#!/bin/sh 
#

if [ $# -ne 1 ]; then
    echo "TCBtools.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCBtools.$1

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCBtools.sh: build test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCBtools.sh: build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCBtools.sh: this build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=${CLM_ROOT}/tools/$1
blddir=${CLM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCBtools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

echo "TCBtools.sh: building $1 executable; output in ${CLM_TESTDIR}/${test_name}/test.log" 
#
# Copy build files over
#
cp $cfgdir/Makefile .
cp $cfgdir/Srcfiles .
cp $cfgdir/*.h      .
#
# Add cfgdir path to begining of each path in Filepath
#
touch Filepath
while read filepath_arg; do
    echo "${cfgdir}/${filepath_arg}" >> Filepath
done < ${cfgdir}/Filepath


attempt=1
still_compiling="TRUE"
while [ $still_compiling = "TRUE" ]; do

    echo "TCBtools.sh: call to make:" 
    echo "        ${MAKE_CMD} ${TOOLS_MAKE_STRING}" 
    ${MAKE_CMD} ${TOOLS_MAKE_STRING} >> test.log 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then
	echo "TCBtools.sh: make was successful" 
	echo "TCBtools.sh: configure and build test passed"
	echo "PASS" > TestStatus
	if [ $CLM_RETAIN_FILES != "TRUE" ]; then
	    echo "TCBtools.sh: removing some unneeded files to save disc space" 
	    rm *.o
	    rm *.mod
	fi
	still_compiling="FALSE"
    elif [ $attempt -lt 10 ] && \
        grep -c "LICENSE MANAGER PROBLEM" test.log > /dev/null; then
        attempt=`expr $attempt + 1`
        echo "TCBtools.sh: encountered License Manager Problem; launching attempt #$attempt"
    else
	echo "TCBtools.sh: clm build failed, error from make= $rc" 
	echo "TCBtools.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 6
    fi
done

exit 0
