#!/bin/sh 
#

if [ $# -ne 5 ]; then
    echo "TSCext_ccsmseq_scam.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSCext_ccsmseq_scam.$1.$2.$3.$4.$5
cfg_name=TCMscam.$1

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSCext_ccsmseq_scam.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSCext_ccsmseq_scam.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSCext_ccsmseq_scam.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSCext_ccsmseq_scam.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSCext_ccsmseq_scam.sh: calling TSMext_ccsmseq_cam.sh to generate iop datafiles" 
${CLM_SCRIPTDIR}/TSMext_ccsmseq_cam.sh $1 $2 $5
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSCext_ccsmseq_scam.sh: error from TSMext_ccsmseq_cam.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

#temporarily stage these files one level up in work directory tree
echo "TSCext_ccsmseq_scam.sh: stage output files in $CLM_TESTDIR"
file=${CLM_TESTDIR}/TSMext_ccsmseq_cam.$1.$2.$5/camrun.cam2.i.0000-09-01-00000.nc
echo "cp -f $file ${CLM_TESTDIR}/."
if [ ! -f $file ]; then
    echo "TSCext_ccsmseq_scam.sh: error -- $file does not exist as expected"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi
cp -f $file ${CLM_TESTDIR}/.

file=${CLM_TESTDIR}/TSMext_ccsmseq_cam.$1.$2.$5/camrun.cam2.h1.0000-09-01-00000.nc
echo "cp -f $file ${CLM_TESTDIR}/."
if [ ! -f $file ]; then
    echo "TSCext_ccsmseq_scam.sh: error -- $file does not exist as expected"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi
cp -f $file ${CLM_TESTDIR}/.

echo "TSCext_ccsmseq_scam.sh: now run scam with the generated datafiles as input" 
if [ "$debug" != "YES" ]; then
  ${CLM_SCRIPTDIR}/TSMext_ccsmseq_cam.sh $3 $4 $5
  rc=$?
else
  rc=0
fi

if [ $rc -ne 0 ]; then
    echo "TSCext_ccsmseq_scam.sh: error from TSMext_ccsmseq_cam.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

#remove temporarily staged files
if [ "$CAM_RETAIN_FILES" != "TRUE" ]; then
    rm ../camrun.*.nc
fi

# Now test the output
echo "TSCext_ccsmseq_scam.sh: Comparing answers to ensure SCAM gives bit-for-bit answers as CAM ... "
if [ "$debug" != "YES" ]; then
   myvar=`ncdump -ff -p 9,17 -v QDIFF,TDIFF ${CLM_TESTDIR}/TSMext_ccsmseq_cam.$3.$4.$5/camrun.*.h1.*.nc | egrep //\.\*DIFF | sed s/^\ \*// | sed s/\[,\;\].\*\$// | uniq`
else
   myvar=0
fi
if [ "$myvar" == "0" ]; then
    echo "TSCext_ccsmseq_scam.sh:  scam b4b test passed"
    echo "PASS" > TestStatus
else
    echo "TSCext_ccsmseq_scam.sh: scam b4b test did not pass"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

exit 0
