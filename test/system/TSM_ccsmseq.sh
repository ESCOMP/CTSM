#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TSM_ccsmseq.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSM_ccsmseq.$1.$2.$3.${CCSM_MACH}

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSM_ccsmseq.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSM_ccsmseq.sh: test was already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSM_ccsmseq.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSM_ccsmseq.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi
testdir=${CLM_TESTDIR}/${test_name}
if [ -d ${testdir} ]; then
    rm -r ${testdir}
fi
mkdir -p ${testdir}
if [ $? -ne 0 ]; then
    echo "TSM.sh: error, unable to create work subdirectory"
    exit 3
fi

echo "TSM_ccsmseq.sh: calling TCT_ccsmseq.sh to prepare sequential CCSM executable" 
${CLM_SCRIPTDIR}/TCT_ccsmseq.sh $1 $2 $3
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSM_ccsmseq.sh: error from TCT_ccsmseq.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

# Go to CCSM job scripts directory and run the *.test script you find there
echo "change directory to ${CLM_TESTDIR}/$1.$2.$3.${CCSM_MACH}.sc.${JOBID}"
cd ${CLM_TESTDIR}/$1.$2.$3.${CCSM_MACH}.sc.${JOBID}
sandboxscript=`/bin/ls -1 *.test`
sandboxbuild=`/bin/ls -1 *.build`

echo "TSM_ccsmseq.sh: building sequential-CCSM; output in ${CLM_TESTDIR}/${test_name}/test.log" 

if [ "$debug" != "YES" ]; then
   ./${sandboxbuild} > ${CLM_TESTDIR}/${test_name}/test.log 2>&1
   rc=$?
else
   rc=0
fi
if [ $rc -ne 0 ]; then
    echo "TSM_ccsmseq.sh: error from ${sandboxbuild} = $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi
 
echo "TSM_ccsmseq.sh: running sequential-CCSM; output in ${CLM_TESTDIR}/${test_name}/test.log" 

if [ "$debug" != "YES" ] || [ "$compile_only" != "YES" ]; then
   ./${sandboxscript} > ${CLM_TESTDIR}/${test_name}/test.log 2>&1
   rc=$?
else
   echo "GEN" > TestStatus
   rc=0
fi
if [ $rc -ne 0 ]; then
    echo "TSM_ccsmseq.sh: error from ${sandboxscript} = $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi
if [ -f TestStatus ]; then
    if grep -c PASS TestStatus > /dev/null; then
       echo "TSM_ccsmseq.sh: smoke test passed" 
       echo "PASS" > ${CLM_TESTDIR}/${test_name}/TestStatus
    elif grep -c GEN TestStatus > /dev/null; then
       echo "TSM_ccsmseq.sh: smoke test passed" 
       echo "GEN" > ${CLM_TESTDIR}/${test_name}/TestStatus
    else
       echo "TSM_ccsmseq.sh: TestStatus reports an error"
       echo "FAIL.job${JOBID}" > ${CLM_TESTDIR}/${test_name}/TestStatus
       exit 7
    fi
else
    echo "TSM_ccsmseq.sh: TestStatus not reported -- must have been an error"
    echo "FAIL.job${JOBID}" > ${CLM_TESTDIR}/${test_name}/TestStatus
    exit 8
fi

exit 0
