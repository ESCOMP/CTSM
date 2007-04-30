#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TSM.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSM.$1.$2.$3

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSM.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSM.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSM.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=${CLM_SCRIPTDIR}/../../bld
rundir=${CLM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -r ${rundir}
fi
mkdir -p ${rundir} 
if [ $? -ne 0 ]; then
    echo "TSM.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSM.sh: calling TCB.sh to prepare clm executable" 
${CLM_SCRIPTDIR}/TCB.sh $1
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSM.sh: error from TCB.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ ! -f ${CLM_SCRIPTDIR}/nl_files/$2 ]; then
    echo "TSM.sh: namelist options file ${CLM_SCRIPTDIR}/nl_files/$2 not found" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi


echo "TSM.sh: running clm; output in ${CLM_TESTDIR}/${test_name}/test.log" 
echo "TSM.sh: obtaining namelist:" 

run_length=${3}
export run_length

restart_type=0
export restart_type

clm_nrevsn=' '
export $clm_revsn

${CLM_SCRIPTDIR}/nl_files/$2 
rc=$?
if [ $rc -eq 0 ]; then
    echo "TSM.sh: namelist creation was successful" 
    cat ${CLM_SCRIPTDIR}/nl_files/$2
else
    echo "TSM.sh: error building namelist, error= $rc" 
    echo "TSM.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

echo "TSM.sh calling CLM_runcmnd.sh to build run command" 
${CLM_SCRIPTDIR}/CLM_runcmnd.sh $1
rc=$?
if [ $rc -eq 0 ] && [ -f clm_run_command.txt ]; then
    read cmnd < clm_run_command.txt
    echo "TSM.sh: clm run command:" 
    echo "        $cmnd ${CLM_TESTDIR}/TCB.$1/clm"  
    rm clm_run_command.txt
else
    echo "TSM.sh: error building run command; error from CLM_runcmnd.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

${cmnd} ${CLM_TESTDIR}/TCB.$1/clm >> test.log 2>&1
rc=$?
if [ $rc -eq 0 ] && grep -c "TERMINATING CLM MODEL" test.log > /dev/null; then
    echo "TSM.sh: smoke test passed" 
    echo "PASS" > TestStatus
    if [ $CLM_RETAIN_FILES != "TRUE" ]; then
        echo "TSM.sh: removing some unneeded files to save disc space" 
        if [ -f *.clm*.i.* ]; then
            rm *.clm*.i.*
	fi
    fi
else
    echo "TSM.sh: error running clm, error= $rc" 
    echo "TSM.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

exit 0
