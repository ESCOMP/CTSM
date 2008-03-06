#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "TCText_ccsmcon.sh: incorrect number of input arguments" 
    exit 1
fi

echo "TCT.ccsm.sh: creating ccsm test - output from create_test follows..."

cd ${CLM_CCSMROOT}/scripts

test_name=TCText_ccsmcon.${1}.${2}.${CLM_COMPSET}.${CCSM_MACH}

if [ -d ${CLM_TESTDIR}/${test_name} ]; then
    echo "TCT.ccsm.sh: removing old ccsm test directories"
    rm -rf ${CLM_TESTDIR}/${test_name}
fi

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCText_ccsmcon.sh: external concurrent CCSM create test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCText_ccsmcon.sh: external concurrent CCSM create test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCTexternal conccsm.sh: this external concurrent CCSM create test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

blddir=${CLM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCText_ccsmcon.sh: error, unable to create work subdirectory" 
    exit 3
fi

bl_opts=""

hostname=`hostname`
case $hostname in
   ##bluevista
    bv* )
##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to temp file vvvvvvvvvvvvvvvvvvv
## 32 8-way tasks -- this MUST be in sync with the tasks/threads in the test_driver.sh script
## this is also setup to be consistent with the compset specified in test_driver.sh
cat > ${CLM_TESTDIR}/${test_name}/pes.tmp << EOF
set ntasks_atm =  1; set nthrds_atm = 1
set ntasks_lnd = 25; set nthrds_lnd = 1
set ntasks_ice =  1; set nthrds_ice = 1
set ntasks_ocn =  1; set nthrds_ocn = 1
set ntasks_cpl =  4; set nthrds_cpl = 1
EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to temp file ^^^^^^^^^^^^^^^^^^^

    ;;
   ##blueice
    bl* )
##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to temp file vvvvvvvvvvvvvvvvvvv
## 8 32-way tasks -- this MUST be in sync with the tasks/threads in the test_driver.sh script
## this is also setup to be consistent with the compset specified in test_driver.sh
cat > ${CLM_TESTDIR}/${test_name}/pes.tmp << EOF
set ntasks_atm =  1; set nthrds_atm = 1
set ntasks_lnd =  4; set nthrds_lnd = 8
set ntasks_ice =  1; set nthrds_ice = 1
set ntasks_ocn =  1; set nthrds_ocn = 1
set ntasks_cpl =  1; set nthrds_cpl = 1
EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to temp file ^^^^^^^^^^^^^^^^^^^

    ;;
   ## default
    * )
    echo "TCText_ccsmcon.sh: bad hostname (= $hostname) being run on -- use one of the hosts in the TCText_ccsmcon.sh script"
    exit 4

esac

./create_test -test $1 -res $2 -compset ${CLM_COMPSET} -testroot ${CLM_TESTDIR} -mach ${CCSM_MACH} \
    -testid ${JOBID} $bl_opts -pes_file ${CLM_TESTDIR}/${test_name}/pes.tmp -clean off
rc=$?
if [ $rc -ne 0 ]; then
    echo "TCText_ccsmcon.sh: create_test failed, error= $rc"
    exit 5
else
    echo "TCText_ccsmcon.sh: ccsm create test completed successfully"
fi
exit 0
