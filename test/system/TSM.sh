#!/bin/sh 
#

if [ $# -ne 7 ]; then
    echo "TSM.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TSM.$1.$2.$3.$4.$5.$6.${7%+*}

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSM.sh: smoke test has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    elif grep -c GEN ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TSM.sh: test already generated"
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TSM.sh: smoke test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TSM.sh: this smoke test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=${CLM_SCRIPTDIR}/../../../../../bld
rundir=${CLM_TESTDIR}/${test_name}
if [ -d ${rundir} ]; then
    rm -r ${rundir}
fi
mkdir -p ${rundir}/timing/checkpoints
if [ $? -ne 0 ]; then
    echo "TSM.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${rundir}

echo "TSM.sh: calling TCB.sh to prepare clm executable" 
${CLM_SCRIPTDIR}/TCB.sh $1
rc=$?
if [ $rc -ne 0 ]; then
    echo "TSM.sh: error from TCB.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

echo "TSM.sh: running clm; output in ${CLM_TESTDIR}/${test_name}/test.log" 
echo "TSM.sh: obtaining namelist:" 

run_length=${6}
export run_length
atm_rpointer="rpointer.atm"
export atm_rpointer
drv_rpointer="rpointer.drv"
export drv_rpointer
lnd_rpointer="rpointer.lnd"
export lnd_rpointer

if     [ "$7" = "arb_ic" ] || [ "$7" = "cold" ]; then

   nrevsn=' '
   export nrevsn
   drv_restart=' '
   export drv_restart
   datm_restfilm='null'
   export datm_restfilm
   datm_restfils='null'
   export datm_restfils

elif   [ "$7" = "startup" ]; then

   nrevsn=' '
   export nrevsn
   drv_restart=' '
   export drv_restart
   datm_restfilm='null'
   export datm_restfilm
   datm_restfils='null'
   export datm_restfils

elif [ "${7%+*}" = "continue" ]; then

   nrevsn=' '
   export nrevsn
   drv_restart=' '
   export drv_restart
   datm_restfilm='null'
   export datm_restfilm
   datm_restfils='null'
   export datm_restfils

   cp ${CLM_TESTDIR}/TSM.${7#*+}/*clm?.*.*     $rundir/.
   cp ${CLM_TESTDIR}/TSM.${7#*+}/*drv.r*       $rundir/.
   cp ${CLM_TESTDIR}/TSM.${7#*+}/*cpl.r*       $rundir/.
   cp ${CLM_TESTDIR}/TSM.${7#*+}/*datm.r*      $rundir/.
   cp ${CLM_TESTDIR}/TSM.${7#*+}/$atm_rpointer $rundir/.
   cp ${CLM_TESTDIR}/TSM.${7#*+}/$lnd_rpointer $rundir/.
   cp ${CLM_TESTDIR}/TSM.${7#*+}/$drv_rpointer $rundir/.

elif [ "${7%+*}" = "branch" ]; then

   if [ "$debug" = "YES" ] || [ "$compile_only" = "YES" ]; then
       touch ${CLM_TESTDIR}/TSM.${7#*+}/clmrun.clm3.r.1967-01-01-00000.nc
   fi
   cp ${CLM_TESTDIR}/TSM.${7#*+}/*.clm?.r.* $rundir/.
   master_clm_restart=`ls -1rt ${CLM_TESTDIR}/TSM.${7#*+}/*.clm?.r.* \
       | tail -1 | head -1`
   rc=$?
   if [ $rc -ne 0 ]; then
      echo "Can not find clm restart file in ${CLM_TESTDIR}/TSM.${7#*+}"
      exit 6
   fi
   nrevsn=${master_clm_restart}
   export nrevsn

   if [ "$debug" = "YES" ] || [ "$compile_only" = "YES" ]; then
       touch ${CLM_TESTDIR}/TSM.${7#*+}/clmrun.cpl.r.1967-01-01-00000
   fi
   cp ${CLM_TESTDIR}/TSM.${7#*+}/*cpl.r* $rundir/.
   drv_restart=`ls -1rt ${CLM_TESTDIR}/TSM.${7#*+}/*.cpl.r.* \
       | tail -1 | head -1`
   export drv_restart

   if [ "$debug" = "YES" ] || [ "$compile_only" = "YES" ]; then
       touch ${CLM_TESTDIR}/TSM.${7#*+}/clmrun.datm.rs1.1967-01-01-00000.bin
   fi
   cp ${CLM_TESTDIR}/TSM.${7#*+}/*datm.rs.* $rundir/.
   datm_restfils=`ls -1rt ${CLM_TESTDIR}/TSM.${7#*+}/*datm.rs1.* \
       | tail -1 | head -1`
   if [ $rc -ne 0 ]; then
      echo "Can not find datm streams restart file in ${CLM_TESTDIR}/TSM.${7#*+}"
      exit 6
   fi
   export datm_restfils

   if [ "$debug" = "YES" ] || [ "$compile_only" = "YES" ]; then
       touch ${CLM_TESTDIR}/TSM.${7#*+}/clmrun.datm.r.1967-01-01-00000.nc
   fi
   cp ${CLM_TESTDIR}/TSM.${7#*+}/*datm.r.* $rundir/.
   datm_restfilm=`ls -1rt ${CLM_TESTDIR}/TSM.${7#*+}/*datm.r.* \
       | tail -1 | head -1`
   export datm_restfilm

else
    echo "TSM.sh: bad start type = $7, can only handle cold, arb_ic, startup, continue or branch"
    exit 7
fi


config_file="${CLM_TESTDIR}/TCB.$1/config_cache.xml"
if [ "$7" = "arb_ic" ] || [ "$7" = "startup" ] || [ "$7" = "cold" ]; then
   ${CLM_SCRIPTDIR}/mknamelist $2 $3 $4 $5 $config_file ${7%+*}
else
   env CLM_THREADS=$CLM_RESTART_THREADS ${CLM_SCRIPTDIR}/mknamelist $2 $3 $4 $5 $config_file ${7%+*}
fi
rc=$?
if [ $rc -eq 0 ]; then
    echo "TSM.sh: namelist creation was successful" 
    cat *_in
else
    echo "TSM.sh: error building namelist, error= $rc" 
    echo "TSM.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 8
fi

echo "TSM.sh calling CLM_runcmnd.sh to build run command" 
if [ "$7" = "arb_ic" ] || [ "$7" = "startup" ] || [ "$7" = "cold" ]; then
   ${CLM_SCRIPTDIR}/CLM_runcmnd.sh $1
else
   env CLM_THREADS=$CLM_RESTART_THREADS CLM_TASKS=$CLM_RESTART_TASKS ${CLM_SCRIPTDIR}/CLM_runcmnd.sh $1
fi
rc=$?
if [ $rc -eq 0 ] && [ -f clm_run_command.txt ]; then
    read cmnd < clm_run_command.txt
    echo "TSM.sh: clm run command:" 
    echo "        $cmnd ${CLM_TESTDIR}/TCB.$1/clm"  
    rm clm_run_command.txt
else
    echo "TSM.sh: error building run command; error from CLM_runcmnd.sh= $rc" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 9
fi

if [ "$debug" != "YES" ] && [ "$compile_only" != "YES" ]; then
   ${cmnd} ${CLM_TESTDIR}/TCB.$1/clm >> test.log 2>&1
   status="PASS"
else
   echo "${cmnd} ${CLM_TESTDIR}/TCB.$1/clm"
   echo "SUCCESSFUL TERMINATION" >> test.log
   status="GEN"
fi
rc=$?
if [ $rc -eq 0 ] && grep -c "SUCCESSFUL TERMINATION" *test.log > /dev/null; then
    echo "TSM.sh: smoke test passed" 
    echo "$status" > TestStatus
else
    echo "TSM.sh: error running clm, error= $rc" 
    echo "TSM.sh: see ${CLM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 10 
fi

exit 0
