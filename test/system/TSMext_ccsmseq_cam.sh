#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TSMext_ccsmseq_cam.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMext_ccsmseq_cam.$1.$2.$3
cfg_name=TCMext_ccsmseq_cam.$1

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMext_ccsmseq_cam.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMext_ccsmseq_cam.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMext_ccsmseq_cam.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSMext_ccsmseq_cam.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSMext_ccsmseq_cam.sh: calling TCBext_ccsmseq_cam.sh to prepare cam executable" 
${CLM_SCRIPTDIR}/TCBext_ccsmseq_cam.sh $1
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMext_ccsmseq_cam.sh: error from TCBext_ccsmseq_cam.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ ! -f ${CLM_SCRIPTDIR}/nl_files/$2 ]; then
    echo "TSMext_ccsmseq_cam.sh: namelist options file ${CLM_SCRIPTDIR}/nl_files/$2 not found" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi


echo "TSMext_ccsmseq_cam.sh: running external sequential-CCSM with CAM; output in ${CLM_TESTDIR}/${test_name}/test.log" 
echo "TSMext_ccsmseq_cam.sh: obtaining namelist:" 

echo "TSMext_ccsmseq_cam.sh.sh: call to build-namelist:"
echo "        ${CLM_SEQCCSMROOT}/models/atm/cam/bld/build-namelist -s -runtype startup \
    -ignore_ic_date \
    -config ${CLM_TESTDIR}/TCBext_ccsmseq_cam.$1/config_cache.xml -infile ${CLM_SCRIPTDIR}/nl_files/$2 \
    -namelist \"&timemgr_inparm stop_n=$3 stop_option=\'nsteps\' /\""

${CLM_SEQCCSMROOT}/models/atm/cam/bld/build-namelist -s -runtype startup \
    -ignore_ic_date \
    -config ${CLM_TESTDIR}/TCBext_ccsmseq_cam.$1/config_cache.xml \
    -infile ${CLM_SCRIPTDIR}/nl_files/$2 \
    -namelist "&timemgr_inparm stop_n=$3 stop_option='nsteps' /" > test.log
rc=$?
if [ $rc -eq 0 ]; then
    echo "TSMext_ccsmseq_cam.sh: build-namelist was successful"
    cat *_in
else
    echo "TSMext_ccsmseq_cam.sh: error building namelist, error from build-namelist= $rc"
    echo "TSMext_ccsmseq_cam.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

echo "TSMext_ccsmseq_cam.sh calling CAM_runcmnd.sh to build run command" 
export CAM_SCRIPTDIR=${CLM_SCRIPTDIR}
export CAM_THREADS=${CLM_THREADS}
export CAM_TASKS=${CLM_TASKS}

echo "CAM_SCRIPTDIR=${CLM_SCRIPTDIR}"
echo "CAM_THREADS=${CLM_THREADS}"
echo "CAM_TASKS=${CLM_TASKS}"
echo "LSB_MCPU_HOSTS=${LSB_MCPU_HOSTS}"

${CLM_SEQCCSMROOT}/models/atm/cam/test/system/CAM_runcmnd.sh $1
rc=$?
if [ $rc -eq 0 ] && [ -f cam_run_command.txt ]; then
    read cmnd < cam_run_command.txt
    echo "TSMext_ccsmseq_cam.sh: cam run command:" 
    echo "        $cmnd ${CLM_TESTDIR}/TCBext_ccsmseq_cam.$1/cam"  
    rm cam_run_command.txt
else
    echo "TSMext_ccsmseq_cam.sh: error building run command; error from CAM_runcmnd.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

if [ "$debug" != "YES" ]; then
  ${cmnd} ${CLM_TESTDIR}/TCBext_ccsmseq_cam.$1/cam >> test.log 2>&1
  rc=$?
else
  echo "END OF MODEL RUN" > test.log
  rc=0
fi
if [ $rc -eq 0 ] && grep -c "END OF MODEL RUN" test.log > /dev/null; then
    echo "TSMext_ccsmseq_cam.sh: smoke test passed" 
    echo "PASS" > TestStatus
    if [ $CLM_RETAIN_FILES != "TRUE" ]; then
        echo "TSMext_ccsmseq_cam.sh: removing some unneeded files to save disc space" 
        if [ -f "*.clm*.i.*" ]; then
            rm *.clm*.i.*
	fi
    fi
else
    echo "TSMext_ccsmseq_cam.sh: error running cam, error= $rc" 
    echo "TSMext_ccsmseq_cam.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

exit 0
