#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TBR.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TBR.$1.$2.$3

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBR.sh: branch test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}"
        exit 0
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

cfgdir=${CLM_SCRIPTDIR}/../../bld
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

initial_length=${3%+*}
branch_length=${3#*+}

if [ ${initial_length} = $3 ] || [ ${branch_length} = $3 ]; then
    echo "TBR.sh: error processing input argument for run lengths" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi
full_length=`expr $initial_length + $branch_length`

echo "TBR.sh: calling TSM.sh for smoke test of full length ${full_length}"
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $full_length
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBR.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi

echo "TBR.sh: calling TSM.sh for smoke test of initial length ${initial_length}"
${CLM_SCRIPTDIR}/TSM.sh $1 $2 ${initial_length}
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBR.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

cp ${CLM_TESTDIR}/TSM.$1.$2.${initial_length}/*clm* ${rundir}/.
#cp ${CLM_TESTDIR}/TSM.$1.$2.${initial_length}/rpointer.* ${rundir}/.

echo "TBR.sh: branching clm; output in ${CLM_TESTDIR}/${test_name}/test.log"
echo "TBR.sh: call to build-namelist:"

run_length=${branch_length}
export run_length

restart_type=3
export restart_type

master_clm_restart=`ls -1rt ${CLM_TESTDIR}/TSM.$1.$2.${initial_length}/*.clm*.r.* \
    | tail -1 | head -1`
nrevsn=${master_clm_restart}
export nrevsn

${CLM_SCRIPTDIR}/nl_files/$2 
rc=$?
if [ $rc -eq 0 ]; then
    echo "TBR.sh: build-namelist was successful"
    cat lnd.stdin
else
    echo "TBR.sh: error building namelist, error from build-namelist= $rc"
    echo "TBR.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

echo "TBR.sh calling CLM_runcmnd.sh to build run command"
## modify the # of tasks/threads for branch
env CLM_THREADS=${CLM_RESTART_THREADS} CLM_TASKS=${CLM_RESTART_TASKS} \
    ${CLM_SCRIPTDIR}/CLM_runcmnd.sh $1
rc=$?
if [ $rc -eq 0 ] && [ -f clm_run_command.txt ]; then
    read cmnd < clm_run_command.txt
    echo "TBR.sh: clm run command:"
    echo "        $cmnd ${CLM_TESTDIR}/TCB.$1/clm"
    rm clm_run_command.txt
else
    echo "TBR.sh: error building run command; error from CLM_runcmnd.sh= $rc"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

${cmnd} ${CLM_TESTDIR}/TCB.$1/clm >> test.log 2>&1
rc=$?
if [ $rc -eq 0 ] && grep -c "TERMINATING CLM MODEL" test.log > /dev/null; then
    echo "TBR.sh: branch of clm completed successfully"
else
    echo "TBR.sh: error on branch run of clm, error= $rc"
    echo "TBR.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 9
fi

echo "TBR.sh: starting b4b comparisons " 
files_to_compare=`ls *.clm*h*.nc`
if [ -z "${files_to_compare}" ]; then
    echo "TBR.sh: error locating files to compare"
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
        echo "TBR.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
	echo "TBR.sh: error from CLM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	all_comparisons_good="FALSE"
    fi
done

if [ ${all_comparisons_good} = "TRUE" ]; then
    echo "TBR.sh: branch test passed" 
    echo "PASS" > TestStatus
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
