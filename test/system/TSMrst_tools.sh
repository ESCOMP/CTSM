#!/bin/sh 
#

if [ $# -ne 9 ]; then
    echo "TSMrst_tools.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSMrst_tools.$1.$2.$3.$4.$5.$6.$7.$8.$9
cascfg=$1
cfgtol=$2
casrun=$3
casdat=$4
frmres=$5
frmmsk=$6
tores=$7
tomsk=$8
caslen=$9
runtol="tools__o"

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMrst_tools.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSMrst_tools.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSMrst_tools.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSMrst_tools.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
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
    echo "TSMrst_tools.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSMrst_tools.sh: calling TCBtools.sh to prepare $cfgtol executable" 
${CLM_SCRIPTDIR}/TCBtools.sh $cfgtol $runtol
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMrst_tools.sh: error from TCBtools.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

echo "TSMrst_tools.sh: calling TSM.sh to run startup case for high resolution to interpolate from"
echo "${CLM_SCRIPTDIR}/TSM.sh $cascfg $casrun $casdat $frmres $frmmsk $caslen startup"
${CLM_SCRIPTDIR}/TSM.sh $cascfg $casrun $casdat $frmres $frmmsk $caslen startup
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMrst_tools.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi
frmcas="$cascfg.$casrun.$casdat.$frmres.$frmmsk.$caslen.startup"
cd ${CLM_TESTDIR}/TSM.$frmcas
copyfile=`ls -1rt *.clm?.r.*.nc | tail -1 | head -1`
fromfile="$frmres.$frmmsk.clm2.r.$casdat.$caslen.nc"
cp -p $copyfile $rundir/$fromfile
rc=$?
cd $rundir
if [ "$debug" = "YES" ] || [ "$compile_only" = "YES" ]; then
    touch $fromfile
elif [ -z "$fromfile" ] || [$rc -ne 0]; then
    echo "TSMrst_tools.sh: error finding restart file from startup case" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10
fi

echo "TSMrst_tools.sh: calling TSM.sh to run cold case for low resolution to interpolate to"
${CLM_SCRIPTDIR}/TSM.sh $cascfg $casrun $casdat $tores $tomsk $caslen cold
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMrst_tools.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi
tocas="$cascfg.$casrun.$casdat.$tores.$tomsk.$caslen.cold"
cd ${CLM_TESTDIR}/TSM.$tocas
copyfile=`ls -1rt *.clm?.r.*.nc | tail -1 | head -1`
tofile="$tores.$tomsk.clm2.r.$casdat.$caslen.nc"
cp -p $copyfile $rundir/$tofile
rc=$?
cd $rundir
if [ "$debug" = "YES" ] || [ "$compile_only" = "YES" ]; then
    touch $tofile
elif [ -z "$tofile" ] || ["$rc" -ne 0]; then
    echo "TSMrst_tools.sh: error finding restart file from cold case" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10
fi

echo "TSMrst_tools.sh: running $cfgtol; output in ${rundir}/test.log" 

toolrun="env OMP_NUM_THREADS=${CLM_THREADS} ${CLM_TESTDIR}/TCBtools.$cfgtol.$runtol/$cfgtol"

runopt="-i $fromfile -o $tofile"
echo "$toolrun $runopt"
if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
   $toolrun  $runopt >> test.log 2>&1
   status="PASS"
   rc=$?
else
   echo "success" > test.log
   status="GEN"
   rc=0
fi

if [ $rc -eq 0 ] && grep -ci "success" test.log > /dev/null; then
    echo "TSMrst_tools.sh: smoke test passed" 
    echo "$status" > TestStatus
else
    echo "TSMrst_tools.sh: error running $cfgtol, error= $rc" 
    echo "TSMrst_tools.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

tocas="$cascfg.$casrun.$casdat.$tores.$tomsk.$caslen.startup"
mkdir -p ${CLM_TESTDIR}/TSM.$tocas
if [ $? -ne 0 ]; then
    echo "TSMrst_tools.sh: error, unable to create work subdirectory" 
    exit 3
fi

finidat="$rundir/$tofile"
export finidat
echo "TSMrst_tools.sh: calling TSM.sh to run startup case for interpolated resolution using: $tofile"
${CLM_SCRIPTDIR}/TSM.sh $cascfg $casrun $casdat $tores $tomsk $caslen startup
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSMrst_tools.sh: error from TSM.sh=$rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

exit 0
