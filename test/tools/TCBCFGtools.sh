#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "TCBCFGtools.sh: incorrect number of input arguments" 
    exit 1
fi

tool=$(basename $1)
test_name=TCBCFGtools.$tool.$2

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCBCFGtools.sh: build test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCBCFGtools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCBCFGtools.sh: build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCBCFGtools.sh: this build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=`ls -1d ${CLM_ROOT}/tools/${1}`
if [ $? -ne 0 ];then
   cfgdir=`ls -1d ${CIME_ROOT}/tools/mapping/${1}*`
   echo "use: $cfgdir"
fi
blddir=${CLM_TESTDIR}/${test_name}/src
if [ -d ${blddir} ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir}
if [ $? -ne 0 ]; then
    echo "TCBCFGtools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

echo "TCBCFGtools.sh: building $tool executable; output in ${blddir}/test.log" 
#
# Copy build files over
#
cp $cfgdir/src/Makefile .
cp $cfgdir/src/Filepath .
#
# Add cfgdir path to beginning of each path in Filepath
#
touch Filepath
while read filepath_arg; do
    echo "${cfgdir}/src/${filepath_arg}" >> Filepath
done < ${cfgdir}/src/Filepath

#
# Figure out configuration
#
if [ ! -f ${CLM_SCRIPTDIR}/config_files/$tool ]; then
    echo "TCB.sh: configure options file ${CLM_SCRIPTDIR}/config_files/$tool not found"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

##construct string of args to configure
config_string=" "
while read config_arg; do
    config_string="${config_string}${config_arg} "
done < ${CLM_SCRIPTDIR}/config_files/$tool

if [ "$TOOLSLIBS" != "" ]; then
   export SLIBS=$TOOLSLIBS
fi
echo "env CIMEROOT=$CLM_ROOT/cime COMPILER=$CESM_COMP $config_string $CLM_ROOT/cime/tools/configure --macros-format Makefile --machine $CESM_MACH $TOOLS_CONF_STRING"
env       CIMEROOT=$CLM_ROOT/cime COMPILER=$CESM_COMP $config_string $CLM_ROOT/cime/tools/configure --macros-format Makefile --machine $CESM_MACH $TOOLS_CONF_STRING >> test.log 2>&1
rc=$?
if [ $rc -ne 0 ]; then
   echo "TCBCFGtools.sh: configure failed, error from configure= $rc" 
   echo "TCBCFGtools.sh: see ${blddir}/test.log for details"
   echo "FAIL.job${JOBID}" > TestStatus
   exit 5
fi

. $INITMODULES
. ./.env_mach_specific.sh

attempt=1
still_compiling="TRUE"
while [ $still_compiling = "TRUE" ]; do

    echo "TCBCFGtools.sh: call to make:" 
    echo "        ${MAKE_CMD} USER_CPPDEFS=-DLINUX"
    if [ "$debug" != "YES" ]; then
       ${MAKE_CMD} USER_CPPDEFS=-DLINUX >> test.log 2>&1
       status="PASS"
       rc=$?
    else
       status="GEN"
       rc=0
    fi
    if [ $rc -eq 0 ]; then
	echo "TCBCFGtools.sh: make was successful" 
	echo "TCBCFGtools.sh: configure and build test passed"
	echo "$status" > TestStatus
	if [ $CLM_RETAIN_FILES != "TRUE" ]; then
	    echo "TCBCFGtools.sh: removing some unneeded files to save disc space" 
	    rm *.o
	    rm *.mod
	fi
	still_compiling="FALSE"
    elif [ $attempt -lt 10 ] && \
        grep -c "LICENSE MANAGER PROBLEM" test.log > /dev/null; then
        attempt=`expr $attempt + 1`
        echo "TCBCFGtools.sh: encountered License Manager Problem; launching attempt #$attempt"
    else
	echo "TCBCFGtools.sh: clm build failed, error from make= $rc" 
	echo "TCBCFGtools.sh: see ${blddir}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 6
    fi
done
if [ "$TOOLSLIBS" != "" ]; then
   export -n SLIBS
fi

exit 0
