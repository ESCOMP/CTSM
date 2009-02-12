#!/bin/sh 
#

if [ $# -ne 8 ]; then
    echo "TSMpergro.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMpergro.$1.$2.$3.$4.$5.$6.$7.R8

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMpergro.sh: CN spinup test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMpergro.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMpergro.sh: PERGRO test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 1
	else
	    echo "TSMpergro.sh: this PERGRO test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=${CLM_SCRIPTDIR}/../../../../../bld
rundir=${CLM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -r ${rundir}
fi
mkdir -p ${rundir}/timing/checkpoints
if [ $? -ne 0 ]; then
    echo "TSMpergro.sh: error, unable to create work subdirectory" 
    exit 2
fi
cd ${rundir}

echo "TSMpergro.sh: Calling TSM.sh for smoke test"
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $4 $5 $6 $7 $8
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMpergro.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 3
fi

echo "TSMpergro.sh: Calling TSM.sh for smoke test"
${CLM_SCRIPTDIR}/TSM.sh $1 $3 $4 $5 $6 $7 $8
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMpergro.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

mv ${CLM_TESTDIR}/TSM.$1.$2.$4.$5.$6.$7.$8/*.clm?.h*.nc .
echo "TSMpergro.sh: starting comparisons "
files_to_compare=`ls *.clm?.h*.nc`
if [ -z " ${files_to_compare}" ]  && [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
    echo "TSMpergro.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    if [ ! -f $compare_file ]; then
       echo "TSMpergro.sh: trouble finding file $compare_file"
       exit 6
    fi
    if [ ! -f ${CLM_TESTDIR}/TSM.$1.$3.$4.$5.$6.$7.$8/$compare_file ]; then
       echo "TSMpergro.sh: trouble finding file ${CLM_TESTDIR}/TSM.$1.$3.$4.$5.$6.$7.$8/$compare_file"
       exit 7
    fi
    ${CLM_SCRIPTDIR}/CLM_compare.sh \
        ${compare_file} \
        ${CLM_TESTDIR}/TSM.$1.$3.$4.$5.$6.$7.$8/${compare_file}
    rc=$?
    mv cprnc.out cprnc.${compare_file}.out
    echo "TSMpergro.sh: comparison completed; output in ${rundir}/cprnc.${compare_file}.out"
done

exit 0
