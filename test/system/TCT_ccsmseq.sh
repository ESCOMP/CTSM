#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TCT_ccsmseq.sh: incorrect number of input arguments" 
    exit 1
fi

echo "TCT.ccsm.sh: creating seq-ccsm test - output from create_test follows..."

cd ${CLM_SCRIPTDIR}/../../../../../scripts

test_name="TCT_ccsmseq.${1}.${2}.${3}.${CCSM_MACH}"

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCT_ccsmseq.sh: sequential CCSM create test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCT_ccsmseq.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCT_ccsmseq.sh: sequential CCSM create test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCT_ccsmseq.sh: this sequential CCSM create test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

if [ -d ${CLM_TESTDIR}/${test_name} ] && [ $CLM_RETAIN_FILES != "TRUE" ]; then
    echo "TCT.ccsm.sh: removing old ccsm test directories"
    rm -rf ${CLM_TESTDIR}/${test_name}
fi


blddir=${CLM_TESTDIR}/${test_name}

echo "blddir: $blddir"

if [ -d ${blddir} ] && [ $CLM_RETAIN_FILES != "TRUE" ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCT_ccsmseq.sh: error, unable to create work subdirectory" 
    exit 3
fi

echo "cat env_pes"

cat > ${CLM_TESTDIR}/env_pes << EOF
#!/bin/csh -f

set ntasks_atm=$CLM_TASKS
set ntasks_lnd=\$ntasks_atm
set ntasks_ice=\$ntasks_atm
set ntasks_ocn=\$ntasks_atm
set ntasks_cpl=\$ntasks_atm
set ntasks_glc=\$ntasks_atm
set nthrds_atm=1
set rootpe_atm=0
set nthrds_lnd=1
set rootpe_lnd=0
set nthrds_ice=1
set rootpe_ice=0
set nthrds_ocn=1
set rootpe_ocn=0
set nthrds_cpl=1
set rootpe_cpl=0
set nthrds_glc=1
set rootpe_glc=0
EOF

echo "./create_test -testname $1.$2.$3.${CCSM_MACH} -testroot ${CLM_TESTDIR} -testid \"sc.${JOBID}\" -clean off -pes_file ${CLM_TESTDIR}/env_pes"
./create_test -testname $1.$2.$3.${CCSM_MACH} -testroot ${CLM_TESTDIR} \
    -testid "sc.${JOBID}" -clean off -pes_file ${CLM_TESTDIR}/env_pes
rc=$?
echo "rc = $rc"
if [ $rc -ne 0 ]; then
    echo "TCT_ccsmseq.sh: create_test failed, error= $rc"
    echo "FAIL.job${JOBID}" > ${CLM_TESTDIR}/${test_name}/TestStatus
    exit 5
else
    echo "TCT_ccsmseq.sh: sequential ccsm create test completed successfully"
    echo "PASS" > ${CLM_TESTDIR}/${test_name}/TestStatus
fi
exit 0
