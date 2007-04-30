#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TER.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TER.$1.$2.$3

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TER.sh: exact restart test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
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

cfgdir=${CLM_SCRIPTDIR}/../../bld
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

initial_length=${3%+*}
restart_length=${3#*+}
if [ ${initial_length} = $3 ] || [ ${restart_length} = $3 ]; then
    echo "TER.sh: error processing input argument for run lengths" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

full_length=`expr $initial_length + $restart_length`

echo "TER.sh: calling TSM.sh for smoke test of full length ${full_length}" 
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $full_length
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

echo "TER.sh: calling TSM.sh for smoke test of initial length ${initial_length}" 
${CLM_SCRIPTDIR}/TSM.sh $1 $2 ${initial_length}
rc=$?
if [ $rc -ne 0 ]; then
    echo "TER.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi


cp ${CLM_TESTDIR}/TSM.$1.$2.${initial_length}/*clm* ${rundir}/.
cp ${CLM_TESTDIR}/TSM.$1.$2.${initial_length}/rpointer* ${rundir}/.

echo "TER.sh: restarting offline; output in ${CLM_TESTDIR}/${test_name}/test.log" 
echo "TER.sh: call to build-namelist:"

run_length=${restart_length}
export run_length

restart_type=1
export restart_type

clm_nrevsn=' '
export $clm_revsn

${CLM_SCRIPTDIR}/nl_files/$2 
rc=$?
if [ $rc -eq 0 ]; then
    echo "TER.sh: build-namelist was successful"
    cat namelist 
else
    echo "TER.sh: error building namelist, error= $rc"
    echo "TER.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

echo "TER.sh calling CLM_runcmnd.sh to build run command"
## modify the # of tasks/threads for restart
env CLM_THREADS=${CLM_RESTART_THREADS} CLM_TASKS=${CLM_RESTART_TASKS} \
    ${CLM_SCRIPTDIR}/CLM_runcmnd.sh $1
rc=$?
if [ $rc -eq 0 ] && [ -f clm_run_command.txt ]; then
    read cmnd < clm_run_command.txt
    echo "TER.sh: clm run command:"
    echo "        $cmnd ${CLM_TESTDIR}/TCB.$1/clm"
    rm clm_run_command.txt
else
    echo "TER.sh: error building run command; error from CLM_runcmnd.sh= $rc"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

${cmnd} ${CLM_TESTDIR}/TCB.$1/clm >> test.log 2>&1
rc=$?
if [ $rc -eq 0 ] && grep -c "TERMINATING CLM MODEL" test.log > /dev/null; then
    echo "TER.sh: restart of offline clm completed successfully"
else
    echo "TER.sh: error on restart run of offline clm, error= $rc" 
    echo "TER.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 9
fi

echo "TER.sh: starting b4b comparisons " 
files_to_compare=`ls *.clm*h*.nc`
if [ -z "${files_to_compare}" ]; then
    echo "TER.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    ${CLM_SCRIPTDIR}/CLM_compare.sh \
	${compare_file} \
	${CLM_TESTDIR}/TSM.$1.$2.${full_length}/${compare_file}
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
    echo "PASS" > TestStatus
    if [ $CLM_RETAIN_FILES != "TRUE" ]; then
        echo "TER.sh: removing some unneeded files to save disc space" 
        rm *.nc
        rm *.r*
    fi
else
    echo "TER.sh: at least one file comparison did not pass" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 11
fi

exit 0
