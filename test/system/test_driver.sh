#!/bin/sh 
#

# test_driver.sh:  driver script for the testing of CLM in Sequential CCSM
#
# usage on bluesky, bangkok, calgary, tempest, bluevista, lightning: 
# ./test_driver.sh
#
# valid arguments: 
# -i    interactive usage
# -f    force batch submission (avoids user prompt)
# -h    displays this help message
#
# **pass environment variables by preceding above commands 
#   with 'env var1=setting var2=setting '
# **more details in the CLM testing user's guide, accessible 
#   from the CLM developers web page


#will attach timestamp onto end of script name to prevent overwriting
cur_time=`date '+%H:%M:%S'`

hostname=`hostname`
case $hostname in

    ##bluesky
    bs* )
    submit_script="test_driver_bluesky_${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
# bluesky
# @ class       = share
# Number of nodes (CHANGE THIS if needed)
# @ node        = 1
# Switch to use (CHANGE THIS if needed)
# @ network.MPI = csss,shared,us
# @ output      = test_dr.o\$(jobid)
# @ error       = test_dr.o\$(jobid)
# @ node_usage  = shared
# @ job_type    = parallel
# @ tasks_per_node = 8
# @ account_no = 93300370
# Export all Environment variables
# @ environment = COPY_ALL
# @ queue
#

if [ -n "\$LOADL_JOB_NAME" ]; then   #batch job
    export JOBID=\`echo \${LOADL_JOB_NAME} | cut -f2 -d'.'\`
    initdir=\${LOADL_STEP_INITDIR}
    interactive="NO"
else
    interactive="YES"
fi

##omp threads
export CLM_THREADS=4
export CLM_RESTART_THREADS=3

##mpi tasks; ignored by load-leveller!
export CLM_TASKS=8
export CLM_RESTART_TASKS=4

export INC_NETCDF=/usr/local/include
export LIB_NETCDF=/usr/local/lib64/r4i4
export AIXTHREAD_SCOPE=S
export MALLOCMULTIHEAP=true
export OMP_DYNAMIC=false
export XLSMPOPTS="stack=40000000"
export MAKE_CMD="gmake -j 32"
export CFG_STRING=""
export CCSM_MACH="bluesky32"
export MACH_WORKSPACE="/ptmp"
export CPRNC_EXE=/contrib/newcprnc3.0/bin/newcprnc
dataroot="/fs/cgd/csm"
echo_arg=""
input_file="tests_pretag_bluesky"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##bluevista
    bv* )
    submit_script="test_driver_bluevista_${cur_time}.sh"

    account_name=`grep -i "^${LOGNAME}:" /etc/project.ncar | cut -f 1 -d "," | cut -f 2 -d ":" `
    if [ ! -n "${account_name}" ]; then
	echo "ERROR: unable to locate an account number to charge for this job under user: $LOGNAME"
	exit 2
    fi

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

#BSUB -a poe                    # use LSF poe elim
#BSUB -x                          # exclusive use of node (not_shared)
#BSUB -n 32                       # total tasks needed
#BSUB -R "span[ptile=16]"          # max number of tasks (MPI) per node
#BSUB -o test_dr.o%J              # output filename
#BSUB -e test_dr.o%J              # error filename
#BSUB -q premium                  # queue
#BSUB -W 4:28                     
#BSUB -P $account_name      

if [ -n "\$LSB_JOBID" ]; then   #batch job
    export JOBID=\${LSB_JOBID}
    initdir=\${LS_SUBCWD}
    interactive="NO"
else
    interactive="YES"
fi

##omp threads
export CLM_THREADS=4
export CLM_RESTART_THREADS=8

##mpi tasks
export CLM_TASKS=8
export CLM_RESTART_TASKS=4

export INC_NETCDF=/usr/local/include
export LIB_NETCDF=/usr/local/lib64/r4i4
export AIXTHREAD_SCOPE=S
export MALLOCMULTIHEAP=true
export OMP_DYNAMIC=false
export XLSMPOPTS="stack=40000000"
export MAKE_CMD="gmake -j"
export CFG_STRING=""
export CCSM_MACH="bluevista16"
export MACH_WORKSPACE="/ptmp"
export CPRNC_EXE=/contrib/newcprnc3.0/bin/newcprnc
dataroot="/fs/cgd/csm"
echo_arg=""
input_file="tests_pretag_bluevista"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##lightning
    ln* )
    submit_script="test_driver_lightning_${cur_time}.sh"

    account_name=`grep -i "^${LOGNAME}:" /etc/project | cut -f 1 -d "," | cut -f 2 -d ":" `
    if [ ! -n "${account_name}" ]; then
	echo "ERROR: unable to locate an account number to charge for this job under user: $LOGNAME"
	exit 2
    fi

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

#BSUB -a mpich_gm            #lightning requirement
#BSUB -x                          # exclusive use of node (not_shared)
#BSUB -n 4                       # total tasks needed
#BSUB -o test_dr.o%J              # output filename
#BSUB -e test_dr.o%J              # error filename
#BSUB -q regular                  # queue
#BSUB -W 1:58                     
#BSUB -P $account_name      

if [ -n "\$LSB_JOBID" ]; then   #batch job
    export JOBID=\${LSB_JOBID}
    initdir=\${LS_SUBCWD}
    interactive="NO"
else
    interactive="YES"
fi

##omp threads
export CLM_THREADS=1
export CLM_RESTART_THREADS=2

##mpi tasks
export CLM_TASKS=4
export CLM_RESTART_TASKS=2

export INC_NETCDF=/contrib/2.6/netcdf/3.6.0-p1-pathscale-2.4-64/include
export LIB_NETCDF=/contrib/2.6/netcdf/3.6.0-p1-pathscale-2.4-64/lib
mpich=/contrib/2.6/mpich-gm/1.2.6..14a-pathscale-2.4-64
export INC_MPI=\${mpich}/include
export LIB_MPI=\${mpich}/lib
export PS=/contrib/2.6/pathscale/2.4
export PATH=\${mpich}/bin:\${PS}/bin:\${PATH}
export LD_LIBRARY_PATH=\${PS}/lib/2.4:\${LD_LIBRARY_PATH}
export MAKE_CMD="gmake -j 2"
export CFG_STRING="-fc pathf90 -linker mpif90 "
export MACH_WORKSPACE="/ptmp"
export CPRNC_EXE=/contrib/newcprnc3.0/bin/newcprnc
dataroot="/fs/cgd/csm"
echo_arg="-e"
input_file="tests_posttag_lightning"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##bangkok,calgary
    ba* | b0* | ca* | c0* ) 
    submit_script="test_driver_bangkok_${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
#PBS -q long
# Number of nodes (CHANGE THIS if needed)
#PBS -l nodes=2:ppn=2
# output file base name
#PBS -N test_dr
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive="NO"
else
    interactive="YES"
fi

##omp threads
export CLM_THREADS=1
export CLM_RESTART_THREADS=2

##mpi tasks
export CLM_TASKS=4
export CLM_RESTART_TASKS=2

if [ "\$CLM_FC" = "PGI" ]; then
    export PGI=/usr/local/pgi-pgcc-pghf-6.1-3
    export INC_NETCDF=/usr/local/netcdf-3.6.1-beta3-pgi-hpf-cc-6.0-5/include
    export LIB_NETCDF=/usr/local/netcdf-3.6.1-beta3-pgi-hpf-cc-6.0-5/lib
    mpich=/usr/local/mpich-1.2.7p1-pgi-pgcc-pghf-6.1-3
    export INC_MPI=\${mpich}/include
    export LIB_MPI=\${mpich}/lib
    export PATH=\${PGI}/linux86/6.1/bin:\${mpich}/bin:\${PATH}
    export LD_LIBRARY_PATH=\${PGI}/linux86/6.1/lib:\${LD_LIBRARY_PATH}
    export CFG_STRING=""
else
    export LAHEY=/usr/local/lf9562
    export INC_NETCDF=/usr/local/netcdf-3.6.1beta3-gcc-4.0.2-g77-lf9562/include
    export LIB_NETCDF=/usr/local/netcdf-3.6.1beta3-gcc-4.0.2-g77-lf9562/lib
    mpich=/usr/local/mpich-1.2.7p1-gcc-g++-4.0.2-8-lf9562
    export INC_MPI=\${mpich}/include
    export LIB_MPI=\${mpich}/lib
    export PATH=\${LAHEY}/bin:\${mpich}/bin:\${PATH}
    export CFG_STRING="-fc lf95 "
fi
export MAKE_CMD="gmake -j"   ##using hyper-threading on calgary
export MACH_WORKSPACE="/scratch/cluster"
export CPRNC_EXE=/contrib/newcprnc3.0/bin/newcprnc
dataroot="/fs/cgd/csm"
echo_arg="-e"
input_file="tests_pretag_bangkok"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;


    ##robin/phoenix
    ro* ) 
    submit_script="test_driver_robin_${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
# #PBS -q debug
# Number of nodes (CHANGE THIS if needed)
#PBS -l walltime=02:58:00,mppe=8
# output file base name
#PBS -N test_dr
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
#PBS -A CLI017
# End of options

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive="NO"
    input_file="tests_posttag_phoenix"
    echo_arg=""
else
    interactive="YES"
    input_file="tests_posttag_robin"
    echo_arg="-e"
fi

##omp threads
export CLM_THREADS=1
export CLM_RESTART_THREADS=2

##mpi tasks
export CLM_TASKS=8
export CLM_RESTART_TASKS=4

. /opt/modules/modules/init/sh
module purge
module load open
module load PrgEnv.5407
module unload mpt
module load mpt.2.4.0.6
module load pbs
module load netcdf/3.5.1_r4

netcdf=\$NETCDF_SV2
export LIB_NETCDF=\${netcdf}/lib
export INC_NETCDF=\${netcdf}/include
export MOD_NETCDF=\${netcdf}/include
export CFG_STRING="-pcols 258 -target_os unicosmp -cppdefs \"-DSYSUNICOS\" "
export CCSM_MACH="phoenix"

export MAKE_CMD="gmake -j 2"   ##using hyper-threading on calgary
export MACH_WORKSPACE="/scratch/scr101"
export CPRNC_EXE=/spin/proj/ccsm/models/atm/cam/bin/newcprnc/cprnc
dataroot="/spin/proj/ccsm"
EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##tempest
    te* )
    submit_script="test_driver_tempest_${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

#QSUB -q ded_16        # Name of the queue (CHANGE THIS if needed)
#QSUB -l mpp_p=20      # Maximum number of processes (CHANGE THIS if needed)
#QSUB -x               # Export all Environment variables
#QSUB -eo              # Put standard error and standard out in same file
#QSUB -J y             # Put job log in its own file
#QSUB                  # End of options

if [ -n "\$QSUB_REQID" ]; then  #batch job
    export JOBID=\`echo \${QSUB_REQID} | cut -f1 -d'.'\`
    initdir=\${QSUB_WORKDIR}
    interactive="NO"
else
    interactive="YES"
fi

##omp threads
export CLM_THREADS=4
export CLM_RESTART_THREADS=8

##mpi tasks
export CLM_TASKS=4
export CLM_RESTART_TASKS=2

export INC_NETCDF=/usr/local/include
export LIB_NETCDF=/usr/local/lib64/r4i4
mpich=/opt/mpt/mpt/usr
export INC_MPI=\${mpich}/include
export LIB_MPI=\${mpich}/lib64
export OMP_DYNAMIC="FALSE"
export _DSM_PLACEMENT="ROUND_ROBIN"
export _DSM_WAIT="SPIN"
export MPC_GANG="OFF"
export MP_SLAVE_STACKSIZE="40000000"
## MIPSpro module required for f90
. /opt/modules/modules/init/sh
module purge
module load MIPSpro mpt
export MAKE_CMD="gmake -j 16"
export CFG_STRING=""
export MACH_WORKSPACE="/ptmp"
export CPRNC_EXE=/contrib/newcprnc3.0/bin/newcprnc
dataroot="/fs/cgd/csm"
echo_arg=""
input_file="tests_pretag_tempest"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    * ) echo "ERROR: machine $hostname not currently supported"; exit 1 ;;
esac

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat >> ./${submit_script} << EOF

##check if interactive job
if [ "\$interactive" = "YES" ]; then

    echo "test_driver.sh: interactive run - setting JOBID to \$\$"
    export JOBID=\$\$
    if [ \$0 = "test_driver.sh" ]; then
	initdir="."
    else
	initdir=\${0%/*}
    fi
fi

##establish script dir and clm_root
if [ -f \${initdir}/test_driver.sh ]; then
    export CLM_SCRIPTDIR=\`cd \${initdir}; pwd \`
    export CLM_ROOT=\`cd \${CLM_SCRIPTDIR}/../../../../.. ; pwd \`
else
    if [ -n "\${CLM_ROOT}" ] && [ -f \${CLM_ROOT}/test/system/test_driver.sh ]; then
	export CLM_SCRIPTDIR=\`cd \${CLM_ROOT}/test/system; pwd \`
    else
	echo "ERROR: unable to determine script directory "
	echo "       if initiating batch job from directory other than the one containing test_driver.sh, "
	echo "       you must set the environment variable CLM_ROOT to the full path of directory containing "
        echo "       <models>. "
	exit 3
    fi
fi

##output files
clm_log=\${initdir}/td.\${JOBID}.log
if [ -f \$clm_log ]; then
    rm \$clm_log
fi
clm_status=\${initdir}/td.\${JOBID}.status
if [ -f \$clm_status ]; then
    rm \$clm_status
fi

##setup test work directory
if [ -z "\$CLM_TESTDIR" ]; then
    export CLM_TESTDIR=\${MACH_WORKSPACE}/\$LOGNAME/test-driver.\${JOBID}
    if [ -d \$CLM_TESTDIR ]; then
        rm -r \$CLM_TESTDIR
    fi
fi
if [ ! -d \$CLM_TESTDIR ]; then
    mkdir -p \$CLM_TESTDIR
    if [ \$? -ne 0 ]; then
	echo "ERROR: unable to create work directory \$CLM_TESTDIR"
	exit 4
    fi
fi

##set our own environment vars
export CSMDATA=\${dataroot}/inputdata/lnd/clm2
export MPI_TYPE_MAX=100000

##process other env vars possibly coming in
if [ -z "\$CLM_RETAIN_FILES" ]; then
    export CLM_RETAIN_FILES=FALSE
fi
if [ -z "\$CLM_CCSMROOT" ]; then
    export CLM_CCSMROOT=\`ls -1td \${dataroot}/collections/ccsm3_1_beta* | head -1\`
fi
if [ -n "\${CLM_INPUT_TESTS}" ]; then
    input_file=\$CLM_INPUT_TESTS
else
    input_file=\${CLM_SCRIPTDIR}/\${input_file}
fi
if [ ! -f \${input_file} ]; then
    echo "ERROR: unable to locate input file \${input_file}"
    exit 5
fi

if [ \$interactive = "YES" ]; then
    echo "reading tests from \${input_file}"
else
    echo "reading tests from \${input_file}" >> \${clm_log}
fi

num_tests=\`wc -w < \${input_file}\`
echo "STATUS OF CLM TESTING UNDER JOB \${JOBID};  scheduled to run \$num_tests tests from:" >> \${clm_status}
echo "\$input_file" >> \${clm_status}
echo "" >> \${clm_status}
echo "tests of CCSM will use source code from:" >> \${clm_status}
echo "\$CLM_CCSMROOT" >> \${clm_status}
if [ \$interactive = "NO" ]; then
    echo "see \${clm_log} for more detailed output" >> \${clm_status}
fi
echo "" >> \${clm_status}

test_list=""
while read input_line; do
    test_list="\${test_list}\${input_line} "
done < \${input_file}

##initialize flags, counter
skipped_tests="NO"
pending_tests="NO"
count=0

##loop through the tests of input file
for test_id in \${test_list}; do
    count=\`expr \$count + 1\`
    while [ \${#count} -lt 3 ]; do
        count="0\${count}"
    done

    master_line=\`grep \$test_id \${CLM_SCRIPTDIR}/input_tests_master\`
    status_out=""
    for arg in \${master_line}; do
        status_out="\${status_out}\${arg} "
    done

    test_cmd=\${status_out#* }

    status_out="\${count} \${status_out}"

    if [ \$interactive = "YES" ]; then
        echo ""
        echo "************************************************************"
        echo "\${status_out}"
        echo "************************************************************"
    else
        echo "" >> \${clm_log}
        echo "************************************************************"\
            >> \${clm_log}
        echo "\$status_out" >> \${clm_log}
        echo "************************************************************"\
            >> \${clm_log}
    fi

    if [ \${#status_out} -gt 64 ]; then
        status_out=\`echo "\${status_out}" | cut -c1-70\`
    fi
    while [ \${#status_out} -lt 67 ]; do
        status_out="\${status_out}."
    done

    echo \$echo_arg "\$status_out\c" >> \${clm_status}

    if [ \$interactive = "YES" ]; then
        \${CLM_SCRIPTDIR}/\${test_cmd}
        rc=\$?
    else
        \${CLM_SCRIPTDIR}/\${test_cmd} >> \${clm_log} 2>&1
        rc=\$?
    fi
    if [ \$rc -eq 0 ]; then
        echo "PASS" >> \${clm_status}
    elif [ \$rc -eq 255 ]; then
        echo "SKIPPED*" >> \${clm_status}
        skipped_tests="YES"
    elif [ \$rc -eq 254 ]; then
        echo "PENDING**" >> \${clm_status}
        pending_tests="YES"
    else
        echo "FAIL! rc= \$rc" >> \${clm_status}
	if [ \$interactive = "YES" ]; then
	    if [ "\$CLM_SOFF" != "FALSE" ]; then
		echo "stopping on first failure"
		echo "stopping on first failure" >> \${clm_status}
		exit 6
	    fi
	else
	    if [ "\$CLM_SOFF" = "TRUE" ]; then
		echo "stopping on first failure" >> \${clm_status}
		echo "stopping on first failure" >> \${clm_log}
		exit 6
	    fi
	fi
    fi
done

echo "end of input" >> \${clm_status}
if [ \$interactive = "YES" ]; then
    echo "end of input"
else
    echo "end of input" >> \${clm_log}
fi

if [ \$skipped_tests = "YES" ]; then
    echo "*  please verify that any skipped tests are not required of your clm commit" >> \${clm_status}
fi
if [ \$pending_tests = "YES" ]; then
    echo "** tests that are pending must be checked manually for a successful completion" >> \${clm_status}
    if [ \$interactive = "NO" ]; then
	echo "   see the test's output in \${clm_log} " >> \${clm_status}
	echo "   for the location of test results" >> \${clm_status}
    fi
fi
exit 0

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^


chmod a+x $submit_script
arg1=${1##*-}
case $arg1 in
    [iI]* )
    ./${submit_script}
    exit 0
    ;;

    [fF]* )
    ;;

    "" )
    echo ""
    echo "**********************"
    echo "$submit_script has been created and will be submitted to the batch queue..."
    echo "(ret) to continue, (a) to abort"
    read ans
    case $ans in
	[aA]* ) 
	echo "aborting...type ./test_driver.sh -h for help message"
	exit 0
	;;
    esac
    ;;

    * )
    echo ""
    echo "**********************"
    echo "usage on bluesky, bangkok, tempest, bluevista, lightning, robin: "
    echo "./test_driver.sh"
    echo ""
    echo "valid arguments: "
    echo "-i    interactive usage"
    echo "-f    force batch submission (avoids user prompt)"
    echo "-h    displays this help message"
    echo ""
    echo "**pass environment variables by preceding above commands "
    echo "  with 'env var1=setting var2=setting '"
    echo "**more details in the CLM testing user's guide, accessible "
    echo "  from the CLM developers web page"
    echo ""
    echo "**********************"
    exit 0
    ;;
esac

echo "submitting..."
case $hostname in
    ##bluesky
    bs* )  llsubmit ${submit_script};;

    ##bluevista
    bv* )  bsub < ${submit_script};;

    ##lightning
    ln* )  bsub < ${submit_script};;

    ##bangkok,calgary
    ba* | b0* | ca* | c0* )  qsub ${submit_script};;

    ##robin/phoenix
    ro* )  qsub ${submit_script};;

    ##tempest
    te* )  qsub ${submit_script};;

esac
exit 0
