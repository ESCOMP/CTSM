#!/bin/sh 
#

if [ -z "$BL_ROOT" ] && [ -z "$BL_TESTDIR" ]; then
    echo "TBL.sh: no environment variables set for baseline test - will skip" 
    exit 255
fi

if [ $# -ne 7 ]; then
    echo "TBL.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TBL.$1.$2.$3.$4.$5.$6.$7

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBL.sh: baseline test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}"
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TBL.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TBL.sh: baseline test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TBL.sh: this baseline test failed under job ${prev_jobid} - moving those results to "
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
    echo "TBL.sh: error, unable to create work subdirectory" 
    exit 3
fi

cd ${rundir}

echo "TBL.sh: calling TSM.sh for smoke test"
${CLM_SCRIPTDIR}/TSM.sh $1 $2 $3 $4 $5 $6 $7
rc=$?
if [ $rc -ne 0 ]; then
    echo "TBL.sh: error from TSM.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

if [ -n "${BL_ROOT}" ]; then
    if [ -z "$BL_TESTDIR" ]; then
	BL_TESTDIR=${CLM_TESTDIR}.bl
    fi
    echo "TBL.sh: generating baseline data from root $BL_ROOT - results in $BL_TESTDIR"

    echo "TBL.sh: calling ****baseline**** TSM.sh for smoke test"
    bl_dir=`/bin/ls -1d ${BL_ROOT}/models/lnd/clm/test/system`
    env CLM_TESTDIR=${BL_TESTDIR} \
	CLM_SCRIPTDIR=$bl_dir \
	$bl_dir/TSM.sh $1 $2 $3 $4 $5 $6 $7
    rc=$?
    if [ $rc -ne 0 ]; then
	echo "TBL.sh: error from *baseline* TSM.sh= $rc" 
	echo "FAIL.job${JOBID}" > TestStatus
	exit 5
    fi
fi

echo "TBL.sh: starting b4b comparisons " 
files_to_compare=`cd ${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.$6.$7; ls *.clm?.h*.nc`
if [ -z "${files_to_compare}" ] && [ "$debug" != "YES" ]; then
    echo "TBL.sh: error locating files to compare"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

all_comparisons_good="TRUE"
for compare_file in ${files_to_compare}; do

    ${CLM_SCRIPTDIR}/CLM_compare.sh \
	${BL_TESTDIR}/TSM.$1.$2.$3.$4.$5.$6.$7/${compare_file} \
	${CLM_TESTDIR}/TSM.$1.$2.$3.$4.$5.$6.$7/${compare_file}
    rc=$?
    mv cprnc.out cprnc.${compare_file}.out
    if [ $rc -eq 0 ]; then
        echo "TBL.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
    else
	echo "TBL.sh: error from CLM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
	all_comparisons_good="FALSE"
    fi
done

if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
   status="PASS"
else
   status="GEN"
fi

if [ "${all_comparisons_good}" = "TRUE" ]; then
    echo "TBL.sh: baseline test passed" 
    echo "$status" > TestStatus
    if [ "$CLM_RETAIN_FILES" != "TRUE" ]; then
        echo "TBL.sh: removing some unneeded files to save disc space" 
        rm *.nc
        rm *.r*
    fi
else
    echo "RMS error of last file that do not match"
    grep RMS ${rundir}/cprnc.${compare_file}.out | grep -v 0.0000E+00 | tail -200
    echo "TBL.sh: at least one file comparison did not pass" 
    echo "Last bit of non-zero RMS errors of last file that did not match"
    grep RMS ${rundir}/cprnc.${compare_file}.out | grep -v 0.0000E+00 | tail -200
    echo "FAIL.job${JOBID}" > TestStatus
    exit 7
fi

exit 0
