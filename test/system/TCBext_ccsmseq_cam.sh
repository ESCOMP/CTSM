#!/bin/sh 
#

if [ $# -ne 1 ]; then
    echo "TCBext_ccsmseq_cam.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCBext_ccsmseq_cam.$1

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCBext_ccsmseq_cam.sh: configure and build test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCBext_ccsmseq_cam.sh: configure and build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCBext_ccsmseq_cam.sh: this configure and build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=${CLM_SEQCCSMROOT}/models/atm/cam/bld
blddir=${CLM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCBext_ccsmseq_cam.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

if [ ! -f ${CLM_SCRIPTDIR}/config_files/$1 ]; then
    echo "TCBext_ccsmseq_cam.sh: configure options file ${CLM_SCRIPTDIR}/config_files/$1 not found" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

##construct string of args to configure
config_string=$CFG_STRING
while read config_arg; do
    config_string="${config_string}${config_arg} "
done < ${CLM_SCRIPTDIR}/config_files/$1

# Add user-source to CLM directories
config_string="${config_string} -usr_src ${CLM_ROOT}/src/main,"
config_string="${config_string}${CLM_ROOT}/src/csm_share/shr,"
config_string="${config_string}${CLM_ROOT}/src/csm_share/eshr,"
config_string="${config_string}${CLM_ROOT}/src/csm_share/dshr,"
config_string="${config_string}${CLM_ROOT}/src/biogeochem,"
config_string="${config_string}${CLM_ROOT}/src/biogeophys,"
config_string="${config_string}${CLM_ROOT}/src/riverroute"

echo "TCBext_ccsmseq_cam.sh: building external seq-ccsm executable with CAM; output in ${CLM_TESTDIR}/${test_name}/test.log" 

attempt=1
still_compiling="TRUE"
while [ $still_compiling = "TRUE" ]; do

    echo "TCBext_ccsmseq_cam.sh: call to configure:" 
    echo "        ${cfgdir}/configure ${config_string}" 

    ${cfgdir}/configure ${config_string} > test.log 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then
	echo "TCBext_ccsmseq_cam.sh: configure was successful" 
    else
	echo "TCBext_ccsmseq_cam.sh: external seq-ccsm CAM configure failed, error from configure= $rc" 
	echo "TCBext_ccsmseq_cam.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 5
    fi

    echo "TCBext_ccsmseq_cam.sh: call to make:" 
    echo "        ${MAKE_CMD}" 
    if [ "$debug" != "YES" ]; then
       ${MAKE_CMD} >> test.log 2>&1
       rc=$?
    else
       rc=0
    fi
    if [ $rc -eq 0 ]; then
	echo "TCBext_ccsmseq_cam.sh: make was successful" 
	echo "TCBext_ccsmseq_cam.sh: configure and build test passed"
	echo "PASS" > TestStatus
	if [ $CLM_RETAIN_FILES != "TRUE" ]; then
	    echo "TCBext_ccsmseq_cam.sh: removing some unneeded files to save disc space" 
	    rm *.o
	    rm *.mod
	fi
	still_compiling="FALSE"
    elif [ $attempt -lt 10 ] && \
        grep -c "LICENSE MANAGER PROBLEM" test.log > /dev/null; then
        attempt=`expr $attempt + 1`
        echo "TCBext_ccsmseq_cam.sh: encountered License Manager Problem; launching attempt #$attempt"
    else
	echo "TCBext_ccsmseq_cam.sh: external seq-ccsm build with CAM failed, error from make= $rc" 
	echo "TCBext_ccsmseq_cam.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 6
    fi
done

exit 0
