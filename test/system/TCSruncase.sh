#!/bin/sh 
#

if [ $# -ne 2 ]; then
    echo "TCSruncase.sh: incorrect number of input arguments" 
    exit 1
fi

test_dir=$1
case=$2
echo "TCS creating run case test - output from create_newcase follows..."

cfgdir=${CLM_SCRIPTDIR}/../../bld

test_name=TCSruncase

if [ -d ${CLM_TESTDIR}/${test_name} ]; then
    echo "TCS.ccsm.sh: removing old CLM create_newcase directories"
    rm -rf ${CLM_TESTDIR}/${test_name}
fi

if [ -f ${CLM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CLM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCSruncase.sh: CLM runcase create_newcase has already passed; results are in "
	echo "        ${CLM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CLM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCSruncase.sh: CLM create_newcase test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CLM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCSruncase.sh:  this CLM create_newcase test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CLM_TESTDIR}/${test_name} ${CLM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

blddir=${CLM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCSruncase.sh: error, unable to create work subdirectory" 
    exit 3
fi

hostname=`hostname`
case $hostname in
   ##bluevista
    bv* )

    script_template="$cfgdir/run-ibm.csh"
    ;;
   ##bluefire
    be* )

    script_template="$cfgdir/run-ibm.csh"
    ;;
   ##bangkok,calgary
    ba* | b0* | ca* | c0* ) 

    script_template="$cfgdir/run-pc.csh"
    ;;

   ##lightning
    ln* )

    script_template="$cfgdir/run-lightning.csh"
    ;;

   ## default
    * )
    echo "TCSruncase.sh: bad hostname (= $hostname) being run on -- use one of the hosts in the TCSruncase.sh script"
    exit 4

esac

test_run_dir="$test_dir/runs"
mkdir -p "$test_dir/runs"
if [ $? -ne 0 ]; then
    echo "TCSruncase.sh: error, unable to create work subdirectory: $test_run_dir" 
    exit 3
fi

$cfgdir/create_newcase -case $case -dirofcases $test_dir \
                       -script_template $script_template -wrkdir $test_run_dir
rc=$?
if [ $rc -ne 0 ]; then
    echo "TCSruncase.sh: create_newcase failed, error= $rc"
    exit 5
else

    echo "TCSruncase.sh: create_newcase test completed successfully"
fi
exit 0
