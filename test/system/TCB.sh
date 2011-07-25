#!/bin/sh 
#

if [ $# -ne 1 ]; then
    echo "TCB.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCB.$1

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCB.sh: configure and build test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCB.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCB.sh: configure and build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCB.sh: this configure and build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=${CLM_SCRIPTDIR}/../../bld
blddir=${CLM_TESTDIR}/${test_name}
if [ -d ${blddir} ] && [ $CLM_RETAIN_FILES != "TRUE" ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCB.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

if [ ! -f ${CLM_SCRIPTDIR}/config_files/$1 ]; then
    echo "TCB.sh: configure options file ${CLM_SCRIPTDIR}/config_files/$1 not found" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

##construct string of args to configure
config_string="$CFG_STRING -mach $CESM_MACH -nc_path $NETCDF_PATH "
while read config_arg; do
    config_string="${config_string}${config_arg} "
done < ${CLM_SCRIPTDIR}/config_files/$1

config_string="${config_string} "

echo "TCB.sh: building clm executable; output in ${CLM_TESTDIR}/${test_name}/test.log" 

attempt=1
still_compiling="TRUE"
while [ $still_compiling = "TRUE" ]; do

    echo "TCB.sh: call to configure:" 
    echo "        ${cfgdir}/configure ${config_string}" 

    if [ -f $blddir/config_cache.xml ]; then
       echo "TCB.sh: configure already ran"
    else
       ${cfgdir}/configure ${config_string} > test.log 2>&1
    fi
    rc=$?
    if [ $rc -eq 0 ]; then
	echo "TCB.sh: configure was successful" 
    else
	echo "TCB.sh: clm configure failed, error from configure= $rc" 
	echo "TCB.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 5
    fi

    echo "TCB.sh: call to make:" 
    echo "        ${MAKE_CMD}" 
    if [ "$debug" != "YES" ]; then
      ${MAKE_CMD} >> test.log 2>&1
      status="PASS"
      rc=$?
    else
      status="GEN"
      rc=0
    fi
    if [ $rc -eq 0 ]; then
	echo "TCB.sh: make was successful" 
	echo "TCB.sh: configure and build test passed"
	echo "$status" > TestStatus
	if [ $CLM_RETAIN_FILES != "TRUE" ]; then
	    echo "TCB.sh: removing some unneeded files to save disc space" 
	    rm *.o
	    rm *.mod
	fi
	still_compiling="FALSE"
    elif [ $attempt -lt 10 ] && \
        grep -c "LICENSE MANAGER PROBLEM" test.log > /dev/null; then
        attempt=`expr $attempt + 1`
        echo "TCB.sh: encountered License Manager Problem; launching attempt #$attempt"
    else
	echo "TCB.sh: clm build failed, error from make= $rc" 
	echo "TCB.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 6
    fi
done

exit 0
