#!/bin/sh 
#

if [ $# -ne 8 ]; then
    echo "TRP.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TRP.$1.$2.$3.$4.$5.$6.$7.${8%+*}

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TRP.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TRP.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TRP.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TRP.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

rundir=${CLM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -r ${rundir}
fi
mkdir -p ${rundir}/timing/checkpoints
if [ $? -ne 0 ]; then
    echo "TRP.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}
mkdir -p $rundir/lastrun

ntest=0;
while [ $ntest -le $2 ] ; do

  ntest=`expr $ntest + 1`

  echo "TRP.sh: calling TSM.sh for smoke test"
  ${CLM_SCRIPTDIR}/TSM.sh $1 $3 $4 $5 $6 $7 $8
  rc=$?
  if [ $rc -ne 0 ]; then
      echo "TRP.sh: error from TSM.sh= $rc" 
      echo "FAIL.job${JOBID}" > TestStatus
      exit 4
  fi

  if [ $ntest -lt $2 ] ; then
     rm ${CLM_TESTDIR}/TSM.$1.$3.$4.$5.$6.$7.$8/TestStatus
  fi
  files_to_compare=`cd ${CLM_TESTDIR}/TSM.$1.$3.$4.$5.$6.$7.$8; ls *.clm*h*.nc`
  if [ -z "${files_to_compare}" ] && [ "$debug" != "YES" ]; then
      echo "TRP.sh: error locating files to compare"
      echo "FAIL.job${JOBID}" > TestStatus
      exit 5
  fi
  echo "Files to compare: $files_to_compare"
  cd ${rundir}

  # Compare history files to previous run
  if [ $ntest -gt 1 ]; then
     echo "TRP.sh: starting b4b comparisons on test # $ntest" 
     all_comparisons_good="TRUE"
     for compare_file in ${files_to_compare}; do

         ${CLM_SCRIPTDIR}/CLM_compare.sh \
             ${CLM_TESTDIR}/TSM.$1.$3.$4.$5.$6.$7.$8/${compare_file} \
             ${rundir}/lastrun/${compare_file}
         rc=$?
         mv cprnc.out cprnc.${compare_file}.out
         if [ $rc -eq 0 ]; then
             echo "TRP.sh: comparison successful; output in ${rundir}/cprnc.${compare_file}.out"
         else
             echo "TRP.sh: error from CLM_compare.sh= $rc; see ${rundir}/cprnc.${compare_file}.out for details" 
             all_comparisons_good="FALSE"
         fi
     done

     if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
        status="PASS"
     else
        status="GEN"
     fi

     if [ "${all_comparisons_good}" = "TRUE" ]; then
         echo "TRP.sh: baseline test passed" 
         echo "$status" > TestStatus
     else
         echo "TRP.sh: at least one file comparison did not pass" 
         echo "FAIL.job${JOBID}" > TestStatus
         exit 6
     fi
  fi
  echo "Copy files to $rundir/lastrun"
  for compare_file in ${files_to_compare}; do
    echo "cp -pf ${CLM_TESTDIR}/TSM.$1.$3.$4.$5.$6.$7.$8/${compare_file} $rundir/lastrun"
    cp -pf ${CLM_TESTDIR}/TSM.$1.$3.$4.$5.$6.$7.$8/${compare_file} $rundir/lastrun
  done

done

if [ "$CLM_RETAIN_FILES" != "TRUE" ]; then
   echo "TRP.sh: removing some unneeded files to save disc space" 
   rm $rundir/lastrun/*.nc
   rm $rundir/lastrun/*.r*
fi

exit 0
