#!/bin/sh 
#

if [ $# -ne 1 ]; then
    echo "CLM_runcmnd.sh: incorrect number of input arguments"
    exit 1
fi

if [ ! -f ${CLM_SCRIPTDIR}/config_files/$1 ]; then
    echo "CLM_runcmnd.sh: configure options file ${CLM_SCRIPTDIR}/config_files/$1 not found"
    exit 2
fi

hostname=`hostname`
case $hostname in

    ##bluesky
    bs* )
    ##search config options file for parallelization info; default on aix is hybrid
    ##note: load-leveller will ignore any MP_ environment vars
    if grep -ic NOSPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""                                   
	else
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	fi
    else
	echo "CLM_runcmnd.sh: NOTE...if running in batch mode, MP_* env vars will be ignored by "
        echo "                load-leveller and processor assignment will be governed by "
        echo "                directives in job command file"
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="env MP_NODES=2 MP_PROCS=${CLM_TASKS} \
		MP_SHARED_MEMORY=yes MP_EUILIB=us MP_RMPOOL=1 poe "
	else
            ##hybrid
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} MP_NODES=1 \
		MP_PROCS=${CLM_TASKS} MP_SHARED_MEMORY=yes MP_EUILIB=us \
		MP_RMPOOL=1 poe "
	fi
    fi ;;

    ##bluevista
    bv* )
    ##search config options file for parallelization info; default on aix is hybrid
    if grep -ic NOSPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""                                   
	else
            ##open-mp only
#	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	    cmnd="env LSB_PJL_TASK_GEOMETRY="\{\(0\)\}" OMP_NUM_THREADS=${CLM_THREADS} "
	fi
    else
	num_nodes=`echo $LSB_MCPU_HOSTS | wc -w`
	num_nodes=`expr $num_nodes / 2`
	tpn=`expr $CLM_TASKS / $num_nodes `
	proc=0
	geo_string="\{"
	count1=$num_nodes
	while [ "$count1" != "0" ]; do
	    geo_string="${geo_string}\("
	    count2=$tpn
	    while [ "$count2" != "0" ]; do
		if [ "$count2" != "$tpn" ]; then
		    geo_string="${geo_string}\,"
		fi
		geo_string="${geo_string}$proc"
		proc=`expr $proc + 1`
		count2=`expr $count2 - 1`
	    done
	    geo_string="${geo_string}\)"
	    count1=`expr $count1 - 1`
	done
	geo_string="${geo_string}\}"

	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} mpirun.lsf "
	else
            ##hybrid
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} OMP_NUM_THREADS=${CLM_THREADS} mpirun.lsf "
	fi
    fi ;;

    ##lightning
    ln* )
    ##search config options file for parallelization info; default on linux is mpi
    if grep -ic NOSPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""
	elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	else
            ##serial
	    cmnd=""
	fi
    else
	num_nodes=`echo $LSB_MCPU_HOSTS | wc -w`
	num_nodes=`expr $num_nodes / 2`
	tpn=`expr $CLM_TASKS / $num_nodes `
	proc=0
	geo_string="\{"
	count1=$num_nodes
	while [ "$count1" != "0" ]; do
	    geo_string="${geo_string}\("
	    count2=$tpn
	    while [ "$count2" != "0" ]; do
		if [ "$count2" != "$tpn" ]; then
		    geo_string="${geo_string}\,"
		fi
		geo_string="${geo_string}$proc"
		proc=`expr $proc + 1`
		count2=`expr $count2 - 1`
	    done
	    geo_string="${geo_string}\)"
	    count1=`expr $count1 - 1`
	done
	geo_string="${geo_string}\}"

        if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} mpirun.lsf "
        elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} OMP_NUM_THREADS=${CLM_THREADS} mpirun.lsf "
        else
            ##mpi only
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} mpirun.lsf "
        fi
    fi ;;


    ##blueice
    bl* )
    ##search config options file for parallelization info; default on aix is hybrid
    if grep -ic NOSPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""                                   
	else
            ##open-mp only
#	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	    cmnd="env LSB_PJL_TASK_GEOMETRY="\{\(0\)\}" OMP_NUM_THREADS=${CLM_THREADS} "
	fi
    else
	num_nodes=`echo $LSB_MCPU_HOSTS | wc -w`
	num_nodes=`expr $num_nodes / 2`
	tpn=`expr $CLM_TASKS / $num_nodes `
	proc=0
	geo_string="\{"
	count1=$num_nodes
	while [ "$count1" != "0" ]; do
	    geo_string="${geo_string}\("
	    count2=$tpn
	    while [ "$count2" != "0" ]; do
		if [ "$count2" != "$tpn" ]; then
		    geo_string="${geo_string}\,"
		fi
		geo_string="${geo_string}$proc"
		proc=`expr $proc + 1`
		count2=`expr $count2 - 1`
	    done
	    geo_string="${geo_string}\)"
	    count1=`expr $count1 - 1`
	done
	geo_string="${geo_string}\}"

	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} mpirun.lsf "
	else
            ##hybrid
	    cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} OMP_NUM_THREADS=${CLM_THREADS} mpirun.lsf "
	fi
    fi ;;


    ##bangkok,calgary
    ba* | b0* | ca* | c0* )
    ##search config options file for parallelization info; default on linux is mpi
    if grep -ic NOSPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""
	elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	else
            ##serial
	    cmnd=""
	fi
    else
        if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
            cmnd="mpirun -np ${CLM_TASKS} "
        elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
            cmnd="env OMP_NUM_THREADS=${CLM_THREADS} mpirun -np ${CLM_TASKS} "
        else
            ##mpi only
            cmnd="mpirun -np ${CLM_TASKS} "
        fi
    fi ;;

    ##phoenix
    ph* )
    ##search config options file for parallelization info; default on cray is mpi
    if grep -ic NOSPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""
	elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	else
            ##serial
	    cmnd=""
	fi
    else
        if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="aprun -c core=unlimited -A -n ${CLM_TASKS} "
        elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} aprun -c core=unlimited -A -n ${CLM_TASKS} "
        else
            ##mpi only
	    cmnd="aprun -c core=unlimited -A -n ${CLM_TASKS} "
        fi
    fi ;;

    ##jaguar
    jaguar* )
    ##search config options file for parallelization info; default on XT4 is mpi
    if grep -ic NOSPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""
	elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	else
            ##serial
	    cmnd=""
	fi
    else
        if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="yod -VN -np ${CLM_TASKS} "
        elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} yod -VN -np ${CLM_TASKS} "
        else
            ##mpi only
	    cmnd="yod -VN -np ${CLM_TASKS} "
        fi
    fi ;;

    ##tempest

    ##tempest
    te* )
    ##search config options file for parallelization info; default on irix64 is open-mp
    if grep -ic NOSPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""
	else
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	fi
    elif grep -ic SPMD ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="mpirun -np ${CLM_TASKS} "
	else
            ##hybrid
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} mpirun -np ${CLM_TASKS} "
	fi
    else
	if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##serial
	    cmnd=""
	else
            ##open-mp only
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
	fi
    fi ;;

    * ) 
    echo "CLM_runcmnd.sh: unable to construct run command for unsupported machine $hostname "
    exit 3;;
esac

#store command in temporary file for calling script to access
echo ${cmnd} > clm_run_command.txt
exit 0
