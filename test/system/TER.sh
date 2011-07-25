#!/bin/sh 
#

if [ $# -ne 7 ]; then
    echo "TER.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TER.$1.$2.$3.$4.$5.$6.$7

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TER.sh: exact restart test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TER.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TER.sh: exact restart test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TER.sh: this exact restart test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

rundir=${CLM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -r ${rundir}
fi
mkdir -p ${rundir}
if [ $? -ne 0 ]; then
    echo "TER.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

initial_length=${6%+*}
restart_length=${6#*+}
if [ ${initial_length} = $6 ] || [ ${restart_length} = $6 ]; then
    echo "TER.sh: error processing input argument for run lengths" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

full_length=`expr $initial_length + $restart_length`

echo "TER.sh: calling TSM.sh for smoke test of full length ${full_length}" 
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4 $5 $full_length $7
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

echo "TER.sh: calling TSM.sh for smoke test of initial length ${initial_length}" 
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4 $5 $initial_length $7
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

echo "TER.sh: calling TSM.sh for smoke test of restart length ${restart_length}" 
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4 $5 ${restart_length} \
                         continue+$1.$2.$3.$4.$5.${initial_length}.$7
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

mv ${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.${restart_length}.continue/*.clm?.h*.nc .
echo "TER.sh: starting b4b comparisons " 
files_to_compare=`ls *.clm?.h*.nc`
if [ -z " ${files_to_compare}" ]  && [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
    echo "TER.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10
fi

if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
   status="PASS"
else
   status="GEN"
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    if [ ! -f $compare_file ]; then
       echo "TER.sh: trouble finding file $compare_file"
       exit 11
    fi
    if [ ! -f ${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.${full_length}.$7/$compare_file ]; then
       echo "TER.sh: trouble finding file ${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.${full_length}.$7/$compare_file"
       exit 12
    fi
    ${CLM_SCRIPTDIR}/CLM_compare.sh \
	${compare_file} \
	${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.${full_length}.$7/${compare_file}
    rc=$?
    mv cprnc.out cprnc.${compare_file}.out
    if [ $rc -eq 0 ]; then
        echo "TER.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
	echo "TER.sh: error from CLM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	all_comparisons_good="FALSE"
    fi
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TER.sh: exact restart test passed" 
    echo "$status" > TestStatus
    if [ $CLM_RETAIN_FILES != "TRUE" ]; then
        echo "TER.sh: removing some unneeded files to save disc space" 
        rm *.nc
        rm *.r*
    fi
else
    echo "TER.sh: at least one file comparison did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 13
fi

exit 0
