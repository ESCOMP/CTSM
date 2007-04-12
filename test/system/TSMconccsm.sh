#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "TSMconccsm.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMconccsm.$1.$2.${CLM_COMPSET}.${CCSM_MACH}

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMconccsm.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMconccsm.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMconccsm.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
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

echo "TSMconccsm.sh: calling TCTconccsm.sh to prepare concurrent CCSM executable" 
${CLM_SCRIPTDIR}/TCTconccsm.sh $1 $2
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMconccsm.sh: error from TCTconccsm.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

# Go to CCSM job scripts directory and run the *.test script you find there
echo "change directory to ${CLM_TESTDIR}/$1.$2.${CLM_COMPSET}.${CCSM_MACH}.${JOBID}"
cd ${CLM_TESTDIR}/$1.$2.${CLM_COMPSET}.${CCSM_MACH}.${JOBID}
sandboxscript=`/bin/ls -1 *.test`
sandboxbuild=`/bin/ls -1 *.build`

echo "TSMconccsm.sh: building concurrent-CCSM; output in ${CLM_TESTDIR}/${test_name}/test.log" 

${sandboxbuild} > ${CLM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMconccsm.sh: error from ${sandboxbuild} = $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi
 
echo "TSMconccsm.sh: running concurrent-CCSM; output in ${CLM_TESTDIR}/${test_name}/test.log" 

${sandboxscript} > ${CLM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMconccsm.sh: error from ${sandboxscript} = $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi
if [ -f TestStatus ]; then
    if grep -c PASS TestStatus > /dev/null; then
       echo "TSMconccsm.sh: smoke test passed" 
       echo "PASS" > ${CLM_TESTDIR}/${test_name}/TestStatus
    else
       echo "TSMconccsm.sh: TestStatus reports an error"
       echo "FAIL.job${JOBID}" > ${CLM_TESTDIR}/${test_name}/TestStatus
       exit 7
    fi
else
    echo "TSMconccsm.sh: TestStatus not reported -- must have been an error"
    echo "FAIL.job${JOBID}" > ${CLM_TESTDIR}/${test_name}/TestStatus
    exit 8
fi

exit 0
