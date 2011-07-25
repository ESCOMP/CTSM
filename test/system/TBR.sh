#!/bin/sh 
#

if [ $# -ne 7 ]; then
    echo "TBR.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TBR.$1.$2.$3.$4.$5.$7

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBR.sh: branch test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}"
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBR.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TBR.sh: branch test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TBR.sh: this branch test failed under job ${prev_jobid} - moving those results to "
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
    echo "TBR.sh: error, unable to create work subdirectory" 
    exit 3
fi

cd ${rundir}

initial_length=${6%+*}
branch_length=${6#*+}

if [ ${initial_length} = $6 ] || [ ${branch_length} = $6 ]; then
    echo "TBR.sh: error processing input argument for run lengths" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi
full_length=`expr $initial_length + $branch_length`

echo "TBR.sh: calling TSM.sh for smoke test of full length ${full_length}"
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4 $5 $full_length $7
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBR.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

echo "TBR.sh: calling TSM.sh for smoke test of initial length ${initial_length}"
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4 $5 ${initial_length} $7
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBR.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

echo "TBR.sh: branching clm; output in ${CLM_TESTDIR}/${test_name}/test.log"

start_ymd=`echo $3 | awk -F: '{print $1}'`
dtime=`echo     $3 | awk -F: '{print $2}'`
if [ $branch_length -lt 0 ]; then
  branch_len_days=$((-$branch_length))
else
  branch_len_days=$(($branch_length * $dtime / 86400))
fi
start_ymd=$(( $start_ymd + $branch_len_days))

branch_nlops="$start_ymd:$dtime"

echo "TBR.sh: calling TSM.sh for smoke test of branch length ${branch_length}"
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $branch_nlops $4 $5 ${branch_length} \
            branch+$1.$2.$3.$4.$5.$initial_length.$7
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBR.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi


mv ${CLM_TESTDIR}/TSM.$1.$2.$branch_nlops.$4.$5.${branch_length}.branch/*.clm?.h*.nc .
echo "TBR.sh: starting b4b comparisons " 
files_to_compare=`ls *.clm?.h*.nc`
first_file=`ls *.clm?.h*.nc | head -1`
if [ -z "${files_to_compare}" ] && [ "${debug}" != "YES" ] && [ "$compile_only" != "YES" ]; then
    echo "TBR.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10
fi

if [ "$first_file" = "$files_to_compare" ] && [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
    echo "TBR.sh: only one file to compare -- not enough"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 11
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    if [ ! -f $compare_file ]; then
       echo "TBR.sh: error finding file $compare_file in " `pwd`
       exit 12
    fi
    if [ ! -f ${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.${full_length}.$7/$compare_file ]; then
       echo "TBR.sh: error finding file $compare_file in ${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.${full_length}.$7"
       exit 13
    fi
    if [ "$compare_file" != "$first_file" ]; then
       ${CLM_SCRIPTDIR}/CLM_compare.sh \
	   ${compare_file} \
	   ${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.${full_length}.$7/${compare_file}
       rc=$?
       mv cprnc.out cprnc.${compare_file}.out
       if [ $rc -eq 0 ]; then
           echo "TBR.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
       else
	   echo "TBR.sh: error from CLM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	   all_comparisons_good="FALSE"
       fi
    fi
done

if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
   status="PASS"
else
   status="GEN"
fi

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TBR.sh: branch test passed" 
    echo "$status" > TestStatus
    if [ $CLM_RETAIN_FILES != "TRUE" ]; then
        echo "TBR.sh: removing some unneeded files to save disc space" 
        rm *.nc
        rm *.r*
    fi
else
    echo "TBR.sh: at least one file comparison did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 11
fi

exit 0
