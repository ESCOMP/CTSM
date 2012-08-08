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

   ##bluefire
   be* )
   ##search config options file for parallelization info; default on aix is hybrid
   if grep -ic NOUSE_MPISERIAL ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
      num_nodes=`echo $LSB_MCPU_HOSTS | wc -w`
      num_nodes=`expr $num_nodes / 2`
      proc=0
      count1=$num_nodes
      tpn=`expr $CLM_TASKS / $num_nodes `
      geo_string="\{"
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
         cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} TARGET_CPU_LIST="-1" mpirun.lsf /usr/local/bin/launch"
      else
         ##hybrid
         cmnd="env LSB_PJL_TASK_GEOMETRY=${geo_string} TARGET_CPU_RANGE="-1" OMP_NUM_THREADS=${CLM_THREADS} mpirun.lsf /usr/local/bin/hybrid_launch"
      fi
   else
      if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
         ##serial
         cmnd=""                                   
	   else
         ##open-mp only
         #	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
         cmnd="env LSB_PJL_TASK_GEOMETRY="\{\(0\)\}" OMP_NUM_THREADS=${CLM_THREADS} "
	   fi
   fi ;;
   
   
    ## edinburgh
    edinburgh* | e0* )
    ##search config options file for parallelization info; default on linux is mpi
    if grep -ic NOUSE_MPISERIAL ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
        if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
            cmnd="/usr/local/mpiexec/bin/mpiexec -n ${CLM_TASKS} "
        elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
            cmnd="env OMP_NUM_THREADS=${CLM_THREADS} /usr/local/mpiexec/bin/mpiexec -n ${CLM_TASKS} "
        else
            ##mpi only
            cmnd="/usr/local/mpiexec/bin/mpiexec -n ${CLM_TASKS} "
        fi
    else
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
    fi ;;

    ##lynx
    lynx* )
    ##search config options file for parallelization info; default on XT4 is mpi
    if grep -ic NOUSE_MPISERIAL ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
        if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="aprun -n ${CLM_TASKS} "
        elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} aprun -n ${CLM_TASKS} -d ${CLM_THREADS}"
        else
            ##mpi only
	    cmnd="aprun -n ${CLM_TASKS} "
        fi
    else
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
    fi ;;

    ##jaguarpf
    jaguarpf* | aprunjag* )
    ##search config options file for parallelization info; default on XT4 is mpi
    if grep -ic NOUSE_MPISERIAL ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
        if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="aprun -n ${CLM_TASKS} "
        elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} aprun -n ${CLM_TASKS} -d ${CLM_THREADS}"
        else
            ##mpi only
	    cmnd="aprun -n ${CLM_TASKS} "
        fi
    else
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
    fi ;;

    ##yong
    yong* | vpn* )
    if grep -ic NOUSE_MPISERIAL ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
        if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##mpi only
	    cmnd="mpirun -n ${CLM_TASKS} "
        elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
            ##hybrid
	    cmnd="env OMP_NUM_THREADS=${CLM_THREADS} mpirun -n ${CLM_TASKS} "
        else
            ##mpi only
	    cmnd="mpirun -n ${CLM_TASKS} "
        fi
    else
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
    fi
    ;;

    ##mirage
    mirage* )
    if grep -ic NOSMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
       ##serial
       cmnd=""
    elif grep -ic SMP ${CLM_SCRIPTDIR}/config_files/$1 > /dev/null; then
       ##open-mp only
       cmnd="env OMP_NUM_THREADS=${CLM_THREADS} "
    fi
    ;;

    * ) 
    echo "CLM_runcmnd.sh: unable to construct run command for unsupported machine $hostname "
    exit 3;;
esac

#store command in temporary file for calling script to access
echo ${cmnd} > clm_run_command.txt
exit 0
