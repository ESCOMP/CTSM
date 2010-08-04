#!/bin/sh 
#

# test_driver.sh:  driver script for the offline testing of CLM
#
# usage on breeze, edinburgh, bluefire, intrepid: 
# ./test_driver.sh
#
# usage on jaguar:
# env CLM_JOBID=1001 ./test_driver.sh -c
# env CLM_JOBID=1001 ./test_driver.sh
#
# valid arguments: 
# -i    interactive usage
# -d    debug usage -- display tests that will run -- but do NOT actually execute them
# -f    force batch submission (avoids user prompt)
# -h    displays this help message
#
#
# **pass environment variables by preceding above commands 
#   with 'env var1=setting var2=setting '
# **more details in the CLM testing user's guide, accessible 
#   from the CLM developers web page


#will attach timestamp onto end of script name to prevent overwriting
cur_time=`date '+%H:%M:%S'`

hostname=`hostname`
case $hostname in

    ##bluefire
    be* )
    submit_script="test_driver_bluefire${cur_time}.sh"

    if [ -z "$CLM_ACCOUNT" ]; then
	export CLM_ACCOUNT=`grep -i "^${LOGNAME}:" /etc/project.ncar | cut -f 1 -d "," | cut -f 2 -d ":" `
	if [ -z "${CLM_ACCOUNT}" ]; then
	    echo "ERROR: unable to locate an account number to charge for this job under user: $LOGNAME"
	    exit 2
	fi
    fi
    if [ -z "$CLM_CCSMBLD" ]; then
	export CLM_CCSMBLD="TRUE"
    fi

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

#BSUB -a poe                      # use LSF poe elim
#BSUB -x                          # exclusive use of node (not_shared)
#BSUB -n 192                      # total tasks needed
#BSUB -R "span[ptile=32]"         # max number of tasks (MPI) per node
#BSUB -o test_dr.o%J              # output filename
#BSUB -e test_dr.o%J              # error filename
#BSUB -q regular                  # queue
#BSUB -W 6:00                     
#BSUB -P $CLM_ACCOUNT
#BSUB -J clmtest

if [ -n "\$LSB_JOBID" ]; then   #batch job
    export JOBID=\${LSB_JOBID}
    initdir=\${LS_SUBCWD}
    interactive="NO"
    input_file="tests_pretag_bluefire"
    c_threads=2
    r_threads=4
else
    interactive="YES"
    export LSB_MCPU_HOSTS="\$hostname 8"
    input_file="tests_pretag_bluefire_nompi"
    c_threads=13
    r_threads=25
fi

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=\$c_threads
fi
export CLM_RESTART_THREADS=\$r_threads

##mpi tasks
export CLM_TASKS=96
export CLM_RESTART_TASKS=46

export CLM_COMPSET="I"

export OBJECT_MODE=64
export XLSMPOPTS="stack=256000000"
export OMP_DYNAMIC=FALSE
export AIXTHREAD_SCOPE=S
export MALLOCMULTIHEAP=TRUE
export MP_LABELIO=yes

# MPI Environment
export MP_RC_USE_LMC=yes
export LAPI_DEBUG_RC_WAIT_ON_QP_SETUP=yes
export MP_INFOLEVEL=2
export MP_EUIDEVICE=sn_all
export MP_SHARED_MEMORY=yes
export LAPI_USE_SHM=yes
export MP_EUILIB=us
# commenting out the following line because we believe it will be better to use 
# the defaults, which change with processor count
#export MP_EAGER_LIMIT=32k
export MP_BULK_MIN_MSG_SIZE=64k
export MP_POLLING_INTERVAL=20000000
export MEMORY_AFFINITY=MCM
export LAPI_DEBUG_ENABLE_AFFINITY=YES
export LAPI_DEBUG_BINDPROC_AFFINITY=YES
export MP_SYNC_QP=YES
export MP_RFIFO_SIZE=16777216
export MP_SHM_ATTACH_THRESH=500000
export MP_EUIDEVELOP=min
export MP_USE_BULK_XFER=yes
export MP_BUFFER_MEM=64M

export MP_RC_MAX_QP=8192
export LAPI_DEBUG_RC_DREG_THRESHOLD=1000000
export LAPI_DEBUG_QP_NOTIFICATION=no
export LAPI_DEBUG_RC_INIT_SETUP=no

export INC_NETCDF=/contrib/netcdf-3.6.2/include
export LIB_NETCDF=/contrib/netcdf-3.6.2/lib
export MAKE_CMD="gmake -j 65"
export CCSM_MACH="bluefire"
export CFG_STRING=""
export TOOLS_MAKE_STRING=""
export MACH_WORKSPACE="/ptmp"
CPRNC_EXE="/contrib/newcprnc3.0/bin/newcprnc"
newcprnc="\$MACH_WORKSPACE/\$LOGIN/newcprnc"
/bin/cp -fp \$CPRNC_EXE \$newcprnc
export CPRNC_EXE="\$newcprnc"
export DATM_QIAN_DATA_DIR="/cgd/tss/atm_forcing.datm7.Qian.T62.c080727"
export PFTDATA="/cgd/tss"
dataroot="/fis/cgd/cseg/csm"


echo_arg=""

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##breeze
    breeze | gale | gust | hail )
    submit_script="test_driver_breeze_${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

interactive="YES"

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=2
fi
export CLM_RESTART_THREADS=2

##mpi tasks
export CLM_TASKS=2
export CLM_RESTART_TASKS=1

export CLM_COMPSET="I"

netcdf=/contrib/netcdf-3.6.2/intel
export INC_NETCDF=\$netcdf/include
export LIB_NETCDF=\$netcdf/lib
export intel=/fs/local
export PATH=\${intel}/bin:\${PATH}
export MAKE_CMD="gmake -j5 "
export CCSM_MACH="generic_linux_intel"
export CFG_STRING="-fc ifort -cc icc  -cppdefs '-DFORTRANUNDERSCORE' "
export TOOLS_MAKE_STRING="USER_FC=ifort USER_LINKER=ifort "
export MACH_WORKSPACE="/ptmp"
export CPRNC_EXE=/fs/home/erik/bin/cprnc
export DATM_QIAN_DATA_DIR="/cgd/tss/atm_forcing.datm7.Qian.T62.c080727"
export PFTDATA="/cgd/tss"
dataroot="/fis/cgd/cseg/csm"
echo_arg="-e"
input_file="tests_posttag_breeze"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ## edinburgh
    edinburgh* | e0*) 
    submit_script="test_driver_edinburgh_${cur_time}.sh"
    export PATH=/cluster/torque/bin:${PATH}

    if [ -z "$CLM_CCSMBLD" ]; then
	export CLM_CCSMBLD="TRUE"
    fi

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
#PBS -q long
# Number of nodes (CHANGE THIS if needed)
#PBS -l nodes=1:ppn=8
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
    input_file="tests_pretag_edinburgh"
else
    export CLM_CCSMBLD="FALSE"
    interactive="YES"
    input_file="tests_pretag_edinburgh_nompi"
fi

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=2
fi
export CLM_RESTART_THREADS=1

##mpi tasks
export CLM_TASKS=8
export CLM_RESTART_TASKS=7

export CLM_COMPSET="I"

export PGI=/usr/local/pgi-pgcc-pghf-7.2-5
export LAHEY=/usr/local/lf6481
export INTEL=/usr/local/intel-cluster-3.2.02
export P4_GLOBMEMSIZE=500000000


if [ "\$CLM_FC" = "PGI" ]; then
    netcdf=/usr/local/netcdf-3.6.3-pgi-hpf-cc-7.2-5
    mpich=/usr/local/mpich-1.2.7p1-pgi-hpf-cc-7.2-5
    export LD_LIBRARY_PATH=\${PGI}/linux86/lib:/cluster/torque/lib:\${LD_LIBRARY_PATH}
    export PATH=\${PGI}/linux86/bin:\${mpich}/bin:\${PATH}
    export CCSM_MACH="edinburgh_pgi"
    export CFG_STRING=""
    export TOOLS_MAKE_STRING=""
elif [ "\$CLM_FC" = "INTEL" ]; then
    netcdf=/usr/local/netcdf-3.6.3-intel-3.2.02
    mpich=/usr/local/mpich-1.2.7p1-intel-3.2.02
    export LD_LIBRARY_PATH=/cluster/torque/lib:\${INTEL}/cc/11.0.074/lib/intel64:\${INTEL}/fc/11.0.074/lib/intel64:\${LD_LIBRARY_PATH}
    export PATH=\${INTEL}/fc/11.0.074/bin/intel64:\${INTEL}/cc/11.0.074/bin/intel64:\${mpich}/bin:\${PATH}
    export CCSM_MACH="edinburgh_intel"
    export CFG_STRING="-fc ifort "
    export TOOLS_MAKE_STRING="USER_FC=ifort "
    /usr/local/intel-cluster-3.2.02/intel-login-script.sh
else
    netcdf=/usr/local/netcdf-3.6.3-gcc-4.1.2-lf95-8.0_x86_64
    mpich=/usr/local/mpich-1.2.7p1-gcc-g++-4.1.2-42-lf9581
    export LD_LIBRARY_PATH=\${LAHEY}/lib64:/cluster/torque/lib:\${LD_LIBRARY_PATH}
    export PATH=\${LAHEY}/bin:\${mpich}/bin:\${PATH}
    export CCSM_MACH="edinburgh_lahey"
    export CFG_STRING="-fc lf95 -cc gcc "
    export TOOLS_MAKE_STRING="USER_FC=lf95 USER_LINKER=lf95 "
fi
export INC_NETCDF=\${netcdf}/include
export LIB_NETCDF=\${netcdf}/lib
export INC_MPI=\${mpich}/include
export LIB_MPI=\${mpich}/lib
export MAKE_CMD="gmake -j 5"   ##using hyper-threading on edinburgh
export MACH_WORKSPACE="/scratch/cluster"
export CPRNC_EXE=/fs/cgd/csm/tools/cprnc_64/cprnc
export DATM_QIAN_DATA_DIR="/project/tss/atm_forcing.datm7.Qian.T62.c080727"
export PFTDATA="/project/tss"
dataroot="/fs/cgd/csm"
echo_arg="-e"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##jaguar
    jaguar* ) 
    submit_script="test_driver_jaguar_${cur_time}.sh"

    if [ -z "$CLM_CCSMBLD" ]; then
	export CLM_CCSMBLD="TRUE"
    fi

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
# #PBS -q batch
# Number of nodes (CHANGE THIS if needed)
#PBS -l walltime=04:00:00,size=528
# output file base name
#PBS -N test_dr
# Put standard error and standard out in same file
#PBS -j oe
# Use sh
#PBS -S /bin/sh
# Export all Environment variables
#PBS -V
#PBS -A CLI017dev
# End of options

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

echo_arg="-e"
if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive="NO"
    input_file="tests_pretag_jaguar"
else
    interactive="YES"
    input_file="tests_pretag_jaguar_nompi"
    if [ "\$compile_only" = "YES" ]; then
       input_file="tests_pretag_jaguar"
    fi
fi


##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=2
fi
export CLM_RESTART_THREADS=4

##mpi tasks
export CLM_TASKS=264
export CLM_RESTART_TASKS=129

export CLM_COMPSET="I"

#export PATH="/opt/public/bin:/opt/cray/bin:/usr/bin/X11"
#export PATH="\${PATH}:/usr/bin:/bin:/opt/bin:/sbin:/usr/sbin:/apps/jaguar/bin"
source /opt/modules/default/init/sh
#module load PrgEnv-pgi Base-opts
#module load xtpe-quadcore
#module load torque moab
#module switch pgi pgi/9.0.2               # 9.0.2 tested for bfb on 2009-sep-25
#module switch xt-mpt    xt-mpt/3.2.0      # 3.2.0  is default on 2009-sep-25
#module switch xt-libsci xt-libsci/10.3.5  # 10.3.5 is default on 2009-sep-25
module load pgi/9.0.2                     # 9.0.2 tested for bfb on 2009-sep-25
module load xt-mpt/3.2.0                  # 3.2.0  is default on 2009-sep-25
module load xt-libsci/10.3.5              # 10.3.5 is default on 2009-sep-25
module load   netcdf/3.6.2                # 3.6.2  is default on 2008-sep-03
module load p-netcdf
module load   ncl
module load subversion

export MPICH_MAX_SHORT_MSG_SIZE=32000 # default is 128000 bytes
export MPICH_PTL_UNEX_EVENTS=960000   # default is  90000 (unexpected recv queue size)
export MPICH_UNEX_BUFFER_SIZE=1000M   # default is    60M (unexpected short msgs buff size)
export MPICH_MSGS_PER_PROC=160000     # default is  32768
export MPICH_PTL_SEND_CREDITS=-1

export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1

# These environment variables were suggested by Helen He to help get around compiler issues
# with pgi9
export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=536870912

# The environment variables below produce corefiles and maybe (?) should be
# moved to DEBUG mode at some point
export MPICH_DBMASK=0x200

export LIB_NETCDF=\${NETCDF_DIR}/lib
export INC_NETCDF=\${NETCDF_DIR}/include
export MOD_NETCDF=\${NETCDF_DIR}/include
export INC_PNETCDF=\${PNETCDF_DIR}/include
export LIB_PNETCDF=\${PNETCDF_DIR}/lib
export CCSM_MACH="jaguar"
export CFG_STRING="-fc ftn "
export TOOLS_MAKE_STRING="USER_FC=ftn USER_CC=cc "
export MAKE_CMD="gmake -j 9 "
export MACH_WORKSPACE="/tmp/work"
export CPRNC_EXE=/tmp/proj/ccsm/tools/ccsm_cprnc/cprnc
export DATM_QIAN_DATA_DIR="/tmp/proj/ccsm/inputdata/atm/datm7/atm_forcing.datm7.Qian.T62.c080727"
dataroot="/tmp/proj/ccsm"
export PFTDATA="\$dataroot/inputdata/lnd/clm2/rawdata"
EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##yong
    yong* )
    submit_script="test_driver_yong_${cur_time}.sh"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

interactive="YES"

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=2
fi
export CLM_RESTART_THREADS=1

##mpi tasks
export CLM_TASKS=2
export CLM_RESTART_TASKS=1

export CLM_COMPSET="I"

export CCSM_MACH="yong_g95"
export NETCDF=/usr/local/netcdf-3-6-3
export INC_NETCDF=\$NETCDF/include
export LIB_NETCDF=\$NETCDF/lib
export MAKE_CMD="make -j 4"
export CFG_STRING="-fc g95 -cc gcc "
export TOOLS_MAKE_STRING=""
export MACH_WORKSPACE="$HOME/runs"
export CPRNC_EXE=$HOME/bin/newcprnc
export DATM_QIAN_DATA_DIR="/cgd/tss/atm_forcing.datm7.Qian.T62.c080727"
export PFTDATA="/cgd/tss";
dataroot="$HOME"
echo_arg=""
input_file="tests_posttag_spot1"

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##intrepid
    login* )
    submit_script="test_driver_intrepid_${cur_time}.sh"

    if [ -z "$CLM_CCSMBLD" ]; then
	export CLM_CCSMBLD="TRUE"
    fi

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ./${submit_script} << EOF
#!/bin/sh
#

if [ -n "\$COBALT_JOBID" ]; then    #batch job
    export JOBID=\`echo \${COBALT_JOBID} | cut -f1 -d'.'\`
    initdir=\`pwd\`
    interactive="NO"
    input_file="tests_posttag_intrepid"
else
    interactive="YES"
    input_file="tests_posttag_intrepid_nompi"
fi

##omp threads
if [ -z "\$CLM_THREADS" ]; then   #threads NOT set on command line
   export CLM_THREADS=4
fi
export CLM_RESTART_THREADS=2

##mpi tasks
export CLM_TASKS=256
export CLM_RESTART_TASKS=120

export CLM_COMPSET="I"

export CCSM_MACH="intrepid"

export OBJECT_MODE=32
export OMP_DYNAMIC=FALSE
export AIXTHREAD_SCOPE=S
export MALLOCMULTIHEAP=TRUE
export MPI_TYPE_MAX=100000

export NETCDF=/soft/apps/netcdf-3.6.2
export INC_NETCDF=\$NETCDF/include
export LIB_NETCDF=\$NETCDF/lib
export MAKE_CMD="make -j 5"
export CFG_STRING=""
export TOOLS_MAKE_STRING=""
export MACH_WORKSPACE="/intrepid-fs0/users/$USER/scratch"
dataroot="/gpfs/home/projects/ccsm"
CPRNC_EXE="\$dataroot/tools/cprnc/cprnc"
newcprnc="\$MACH_WORKSPACE/\$LOGIN/newcprnc"
/bin/cp -fp \$CPRNC_EXE \$newcprnc
export CPRNC_EXE="\$newcprnc"
export DATM_QIAN_DATA_DIR="\$dataroot/inputdata/atm/datm7/atm_forcing.datm7.Qian.T62.c080727"
export PFTDATA="\$dataroot/inputdata/lnd/clm2/rawdata"
echo_arg=""

EOF
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    * ) echo "ERROR: machine $hostname not currently supported"; exit 1 ;;
esac

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat >> ./${submit_script} << EOF

if [ -n "\${CLM_JOBID}" ]; then
    export JOBID=\${CLM_JOBID}
fi
##check if interactive job

if [ "\$interactive" = "YES" ]; then

    if [ -z "\${JOBID}" ]; then
       export JOBID=\$\$
    fi
    echo "test_driver.sh: interactive run - setting JOBID to \$JOBID"
    if [ \$0 = "test_driver.sh" ]; then
	initdir="."
    else
	initdir=\${0%/*}
    fi
fi

##establish script dir and clm_root
if [ -f \${initdir}/test_driver.sh ]; then
    export CLM_SCRIPTDIR=\`cd \${initdir}; pwd \`
    export CLM_ROOT=\`cd \${CLM_SCRIPTDIR}/../../../../..; pwd \`
else
    if [ -n "\${CLM_ROOT}" ] && [ -f \${CLM_ROOT}/models/lnd/clm*/test/system/test_driver.sh ]; then
	export CLM_SCRIPTDIR=\`cd \${CLM_ROOT}/models/lnd/clm*/test/system; pwd \`
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
    if [ -d \$CLM_TESTDIR ] && [ \$CLM_RETAIN_FILES != "TRUE" ]; then
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
export CSMDATA=\${dataroot}/inputdata
export DIN_LOC_ROOT=\${CSMDATA}
export MPI_TYPE_MAX=100000

##process other env vars possibly coming in
if [ -z "\$CLM_RETAIN_FILES" ]; then
    export CLM_RETAIN_FILES=FALSE
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
echo " on machine: $hostname" >> \${clm_status}
if [ -n "${BL_ROOT}" ]; then
   echo "tests of baseline will use source code from:" >> \${clm_status}
   echo "\$BL_ROOT" >> \${clm_status}
fi
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

    if [ -z "\$status_out" ]; then
	echo "No test matches \$test_id in \${CLM_SCRIPTDIR}/input_tests_master"
        exit 3
    fi

    test_cmd=\${status_out#* }

    status_out="\${count} \${status_out}"

    if [ \$interactive = "YES" ]; then
        echo ""
        echo "***********************************************************************************"
        echo "\${status_out}"
        echo "***********************************************************************************"
    else
        echo "" >> \${clm_log}
        echo "***********************************************************************************"\
            >> \${clm_log}
        echo "\$status_out" >> \${clm_log}
        echo "***********************************************************************************"\
            >> \${clm_log}
    fi

    if [ \${#status_out} -gt 94 ]; then
        status_out=\`echo "\${status_out}" | cut -c1-100\`
    fi
    while [ \${#status_out} -lt 97 ]; do
        status_out="\${status_out}."
    done

    echo \$echo_arg "\$status_out\c" >> \${clm_status}

    if [   \$interactive = "YES" ]; then
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
if [ ! -z "$CLM_RETAIN_FILES" ]; then
   export CLM_RETAIN_FILES="FALSE"
fi
arg1=${1##*-}
case $arg1 in
    [iI]* )
    debug="NO"
    interactive="YES"
    compile_only="NO"
    export debug
    export interactive
    export compile_only
    ./${submit_script}
    exit 0
    ;;

    [cC]* )
    debug="NO"
    interactive="YES"
    compile_only="YES"
    export debug
    export CLM_RETAIN_FILES="TRUE"
    export interactive
    export compile_only
    export CLM_RETAIN_FILES="TRUE"
    ./${submit_script}
    exit 0
    ;;

    [dD]* )
    debug="YES"
    interactive="YES"
    compile_only="NO"
    export debug
    export interactive
    export compile_only
    ./${submit_script}
    exit 0
    ;;

    [fF]* )
    debug="NO"
    interactive="NO"
    compile_only="NO"
    export debug
    export interactive
    export compile_only
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
    debug="NO"
    interactive="NO"
    compile_only="NO"
    export debug
    export interactive
    export compile_only
    ;;

    * )
    echo ""
    echo "**********************"
    echo "usage on bluefire, edinburgh, breeze, intrepid: "
    echo "./test_driver.sh"
    echo ""
    echo "usage on jaguar: (compile interactively before submitting)"
    echo "env CLM_JOBID=1001 ./test_driver.sh -c"
    echo "env CLM_JOBID=1001 ./test_driver.sh"
    echo ""
    echo "valid arguments: "
    echo "-i    interactive usage"
    echo "-c    compile-only usage (run configure and compile do not run clm)"
    echo "-d    debug-only  usage (run configure and build-namelist do NOT compile or run clm)"
    echo "-f    force batch submission (avoids user prompt)"
    echo "-h    displays this help message"
    echo ""
    echo "**pass environment variables by preceding above commands "
    echo "  with 'env var1=setting var2=setting '"
    echo ""
    echo "**********************"
    exit 0
    ;;
esac

echo "submitting..."
case $hostname in
    ##bluefire
    be* )  bsub < ${submit_script};;

    ##edinburgh
    edinburgh** | e0* )  qsub ${submit_script};;

    ##jaguar
    jaguar* )  qsub ${submit_script};;

    #intrepid
    login* )  qsub -n 256 -t 60 -q prod-devel --mode script ${submit_script};;

    #default
    * )
    echo "no submission capability on this machine"
    exit 0
    ;;

esac
exit 0
