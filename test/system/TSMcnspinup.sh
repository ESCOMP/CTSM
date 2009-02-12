#!/bin/sh 
#

if [ $# -ne 8 ]; then
    echo "TSMcnspinup.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMcnspinup.$1.$2.$3.$4.$5.$6.$7

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMcnspinup.sh: CN spinup test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMcnspinup.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMcnspinup.sh: CN spinup test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMcnspinup.sh: this CN spinup test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSMcnspinup.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

spinup_length=${8%+*}
exitspinup_length=${8#*+}
normal_length=-1
if [ ${spinup_length} = $8 ] || [ ${exitspinup_length} = $8 ]; then
    echo "TSMcnspinup.sh: error processing input argument for run lengths" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

echo "TSMcnspinup.sh: Calling TSM.sh for smoke test of spinup length ${spinup_length}" 
${CLM_SCRIPTDIR}/TSM.sh $1 $4 $5 $6 $7 $spinup_length arb_ic
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMcnspinup.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

echo "TSMcnspinup.sh: Calling TSM.sh for restart test of exit_spinup length ${exitspinup_length}" 
${CLM_SCRIPTDIR}/TSM.sh $2 $4 $5 $6 $7 $exitspinup_length continue+$1.$4.$5.$6.$7.${spinup_length}.arb_ic
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMcnspinup.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

echo "TSMcnspinup.sh: calling TSM.sh for smoke test of normal length ${normal_length}" 
${CLM_SCRIPTDIR}/TSM.sh $3 $4 $5 $6 $7 ${normal_length} continue+$2.$4.$5.$6.$7.${exitspinup_length}.continue
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMcnspinup.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

exit 0
