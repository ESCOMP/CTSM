#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TSMseqccsm.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMseqccsm.$1.$2.$3
cfg_name=TCMseqccsm.$1

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMseqccsm.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMseqccsm.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMseqccsm.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSMseqccsm.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSMseqccsm.sh: calling TCBseqccsm.sh to prepare cam executable" 
${CLM_SCRIPTDIR}/TCBseqccsm.sh $1
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMseqccsm.sh: error from TCBseqccsm.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ ! -f ${CLM_SCRIPTDIR}/nl_files/$2 ]; then
    echo "TSMseqccsm.sh: namelist options file ${CLM_SCRIPTDIR}/nl_files/$2 not found" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi


echo "TSMseqccsm.sh: running sequential-CCSM; output in ${CLM_TESTDIR}/${test_name}/test.log" 
echo "TSMseqccsm.sh: obtaining namelist:" 

echo "TSM.sh: call to build-namelist:"
echo "        ${CLM_SEQCCSMROOT}/models/atm/cam/bld/build-namelist -s -runtype startup \
    -ignore_ic_date \
    -config ${CLM_TESTDIR}/TCBseqccsm.$1/config_cache.xml -infile ${CLM_SCRIPTDIR}/nl_files/$2 \
    -cam_cfg ${CLM_SEQCCSMROOT}/models/atm/cam/bld \
    -namelist \"&timemgr_inparm stop_n=$3 stop_option=\'nsteps\' /\""

${CLM_SEQCCSMROOT}/models/atm/cam/bld/build-namelist -s -runtype startup \
    -ignore_ic_date \
    -config ${CLM_TESTDIR}/TCBseqccsm.$1/config_cache.xml \
    -infile ${CLM_SCRIPTDIR}/nl_files/$2 \
    -cam_cfg ${CLM_SEQCCSMROOT}/models/atm/cam/bld \
    -namelist "&timemgr_inparm stop_n=$3 stop_option='nsteps' /" > test.log
rc=$?
if [ $rc -eq 0 ]; then
    echo "TSMseqccsm.sh: build-namelist was successful"
    cat *_in
else
    echo "TSMseqccsm.sh: error building namelist, error from build-namelist= $rc"
    echo "TSMseqccsm.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

echo "TSMseqccsm.sh calling CLM_runcmnd.sh to build run command" 
${CLM_SCRIPTDIR}/CLM_runcmnd.sh $1
rc=$?
if [ $rc -eq 0 ] && [ -f clm_run_command.txt ]; then
    read cmnd < clm_run_command.txt
    echo "TSMseqccsm.sh: cam run command:" 
    echo "        $cmnd ${CLM_TESTDIR}/TCBseqccsm.$1/cam"  
    rm clm_run_command.txt
else
    echo "TSMseqccsm.sh: error building run command; error from CLM_runcmnd.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

${cmnd} ${CLM_TESTDIR}/TCBseqccsm.$1/cam >> test.log 2>&1
rc=$?
if [ $rc -eq 0 ] && grep -c "END OF MODEL RUN" test.log > /dev/null; then
    echo "TSMseqccsm.sh: smoke test passed" 
    echo "PASS" > TestStatus
    if [ $CLM_RETAIN_FILES != "TRUE" ]; then
        echo "TSMseqccsm.sh: removing some unneeded files to save disc space" 
        if [ -f "*.clm*.i.* *.cam*.i.*" ]; then
            rm *.clm*.i.* *.cam*.i.*
	fi
    fi
else
    echo "TSMseqccsm.sh: error running cam, error= $rc" 
    echo "TSMseqccsm.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

exit 0
