#!/bin/sh 
#

if [ $# -ne 0 ]; then
    echo "TSMruncase.sh: incorrect number of input arguments" 
    exit 1
fi

case_name="clmtestrun"
test_name="TSMruncase"

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMruncase.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMruncase.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMruncase.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMruncase.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi
testdir=${CLM_TESTDIR}/${test_name}
if [ -d ${testdir} ]; then
    rm -r ${testdir}
fi
mkdir -p ${testdir}/timing
if [ $? -ne 0 ]; then
    echo "TSMruncase.sh: error, unable to create work subdirectory"
    exit 3
fi

echo "TSMruncase.sh: calling TCSruncase.sh to create run script" 
${CLM_SCRIPTDIR}/TCSruncase.sh $testdir $case_name
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMruncase.sh: error from TCSruncase.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

# Go to CCSM job scripts directory and run the *.test script you find there
echo "change directory to $testdir/$case_name"
cd $testdir/$case_name
sandboxscript=`/bin/ls -1 *.csh`
chmod +x $sandboxscript

echo "TSMruncase.sh: running CLM run script; output in ${CLM_TESTDIR}/${test_name}/test.log" 

if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
  ./${sandboxscript} > ${CLM_TESTDIR}/${test_name}/test.log 2>&1
  status="PASS"
  rc=$?
else
  status="GEN"
  rc=0
fi
if [ $rc -ne 0 ]; then
    echo "TSMruncase.sh: error from ${sandboxscript} = $rc" 
    echo "FAIL.job${JOBID}" > ${CLM_TESTDIR}/${test_name}/TestStatus
    exit 6
else
    echo "TSMruncase.sh: smoke test passed" 
    echo "$status" > ${CLM_TESTDIR}/${test_name}/TestStatus
fi

exit 0
