#!/bin/sh 
#

if [ $# -ne 5 ]; then
    echo "TSCscam.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSCscam.$1.$2.$3.$4.$5
cfg_name=TCMscam.$1

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSCscam.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSCscam.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSCscam.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSCscam.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSCscam.sh: calling TSMseqccsm.sh to generate iop datafiles" 
${CLM_SCRIPTDIR}/TSMseqccsm.sh $1 $2 $5
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSCscam.sh: error from TSMseqccsm.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

#temporarily stage these files one level up in work directory tree
cp ${CLM_TESTDIR}/TSMseqccsm.$1.$2.$5/camrun.*.i.*.nc ../.
cp ${CLM_TESTDIR}/TSMseqccsm.$1.$2.$5/camrun.*.h1.*.nc ../.

echo "TSCscam.sh: now run scam with the generated datafiles as input" 
${CLM_SCRIPTDIR}/TSMseqccsm.sh $3 $4 $5
rc=$?

#remove temporarily staged files
if [ $CAM_RETAIN_FILES != "TRUE" ]; then
    rm ../camrun.*
fi

if [ $rc -ne 0 ]; then
    echo "TSCscam.sh: error from TSMseqccsm.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

# Now test the output
echo "TSCscam.sh: Comparing answers to ensure SCAM gives bit-for-bit answers as CAM ... "
myvar=`ncdump -ff -p 9,17 -v QDIFF,TDIFF ${CLM_TESTDIR}/TSMseqccsm.$3.$4.$5/camrun.*.h1.*.nc | egrep //\.\*DIFF | sed s/^\ \*// | sed s/\[,\;\].\*\$// | uniq`
if [ "$myvar" == "0" ]; then
    echo "TSCscam.sh:  scam b4b test passed"
    echo "PASS" > TestStatus
else
    echo "TSCscam.sh: scam b4b test did not pass"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

exit 0
