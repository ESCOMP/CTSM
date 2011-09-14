#!/bin/bash
#----------------------------------------------------------------------
#
# mkmapdata.sh
#
# Create needed mapping files for mksurfdata_map and CLM.
# 
# Example to run for an output resolution of 4x5
#
# mkmapdata.sh -r 4x5
#
# valid arguments: 
# -r <res> Output resolution
# -i       interactive usage
# -l       list mapping files required (so can use check_input_data to get them)
# -d       debug usage -- display mkmapdata that will be run but don't execute them
# -v       verbose usage -- log more information on what is happening
# -h       displays this help message
#
# You can also set the following env variables:
#
# SCRIPGRIDFIL - Full pathname of SCRIP grid file to use (in place of resolution)
# ESMFBIN_PATH - Path to ESMF binaries
# CSMDATA ------ Path to CESM input data
# MPIEXEC ------ Name of mpirun executable
# REGRID_PROC -- Number of MPI processors to use
#
#----------------------------------------------------------------------
echo $0
dir=${0%/*}
if [ "$dir" = "$0" ];then
  dir="."
fi
outfilelist="clm.input_data_list"

#----------------------------------------------------------------------
# SET SOME DEFAULTS -- if not set via env variables outside

if [ -z "$CSMDATA" ]; then
   CSMDATA=/fis/cgd/cseg/csm/inputdata
fi
if [ -z "$REGRID_PROC" ]; then
  REGRID_PROC=8
fi
#----------------------------------------------------------------------
# Usage subroutine
usage() {
  echo ""
  echo "**********************"
  echo "usage on bluefire:"
  echo "./mkmapdata.sh"
  echo ""
  echo "valid arguments: "
  echo "-r <res> resolution"
  echo "-i       interactive usage"
  echo "-l       list mapping files required (so can use check_input_data to get them)"
  echo "         (also writes data to $outfilelist)"
  echo "-d       debug-only  usage (don't actually run mkmapdata just echo what would happen)"
  echo "-h       displays this help message"
  echo "-v       verbose usage -- log more information on what is happening"
  echo ""
  echo "**pass environment variables by preceding above commands "
  echo "  with 'env var1=setting var2=setting '"
  echo ""
  echo "**********************"
}
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# runcmd subroutine
runcmd() {
   cmd=$@
   if [ -z "$cmd" ]; then
       echo "No command given to the runcmd function"
       exit 3
   fi
   if [ "$verbose" = "YES" ]; then
       echo "$cmd"
   fi
   if [ "$debug" != "YES" ]; then
       ${cmd}
       rc=$?
   else
       rc=0
   fi
   if [ $rc != 0 ]; then
       echo "Error status returned from mkmapdata script"
       exit 4
   fi
   return 0
}
#----------------------------------------------------------------------

# Process input arguments
interactive="no"
debug="no"
res="10x15"
verbose="no"
list="no"
declare -i narg=1
for arg in $@; do
   ARGV[narg]=$arg
   narg=narg+1
done

declare -i narg=1
until ((narg>${#ARGV[*]})); do
   arg1=${ARGV[narg]}
   arg1=${arg1##*-}
   case $arg1 in
      [iI]* )
      interactive="YES"
      ;;
      [dD]* )
      debug="YES"
      interactive="YES"
      ;;
      [lL]* )
      debug="YES"
      interactive="YES"
      list="YES"
      ;;
      [rR]* )
      narg=narg+1
      res=${ARGV[narg]} 
      ;;
      [vV]* )
      verbose="YES"
      ;;
      [hH]* )
      usage
      exit 0

      ;;
      * )
      echo "ERROR:: invalid argument sent in: ${ARGV[narg]}"
      usage
      exit 1
      ;;
   esac
   narg=narg+1
done


echo "Script to create mapping files for use by mksurfdata_map"
echo "Output resolution will be $res"

# Find the output grid file for this resolution using the XML database
QUERY="$dir/../../bld/queryDefaultNamelist.pl -silent -namelist clmexp -onlyfiles "
QUERY="$QUERY -justvalue -options sim_year=2000 -csmdata $CSMDATA"
if [ ! -z "$SCRIPGRIDFIL" ]; then
   OUTGRID=$SCRIPGRIDFIL
else
   QRYSCRIPGRID="$QUERY -var scripgriddata -res $res -options lmask=nomask"
   if [ "$verbose" = "YES" ]; then
      echo $QRYSCRIPGRID
   fi
   OUTGRID=`$QRYSCRIPGRID`
fi

if [ -z "$OUTGRID" ]; then
   echo "Output grid file was NOT found for this resolution: $res\n";
   exit 1
fi
if [ "$list" = "YES" ]; then
   echo "outgrid = $OUTGRID"
   echo "outgrid = $OUTGRID" > $outfilelist
elif [ ! -f "$OUTGRID" ]; then
   echo "Output SCRIP grid file does NOT exist: $OUTGRID\n";
   echo "Make sure CSMDATA environment variable is set correctly"
   exit 1
fi
#----------------------------------------------------------------------
# Set the creation date for the output files
CDATE="c"`date +%y%m%d`

#
# Set the list of files to create mapping files for
declare -i nfile=1


# List of gri id files for mksurfdata_map
grids=("0.5x0.5_nomask"  "0.5x0.5_USGS"     "0.5x0.5_AVHRR"    "0.5x0.5_MODIS" \
       "5x5min_nomask"   "5x5min_IGBP-GSDP" "10x10min_nomask"                  \
       "10x10min_IGBPmergeICESatGIS")
for gridmask in ${grids[*]}
do
   grd=${gridmask%_*}
   lmask=${gridmask#*_}
   QUERYFIL="$QUERY -var scripgriddata -res $grd -options lmask=$lmask"
   if [ "$verbose" = "YES" ]; then
      echo $QUERYFIL
   fi
   INGRID[nfile]=`$QUERYFIL`
   if [ "$list" = "YES" ]; then
      echo "ingrid = ${INGRID[nfile]}"
      echo "ingrid = ${INGRID[nfile]}" >> $outfilelist
   fi
   OUTGRID[nfile]=$OUTGRID
   OUTFILE[nfile]=map_${grd}_${lmask}_to_${res}_nomask_aave_da_$CDATE.nc
   nfile=nfile+1
done

# RTM grid for clm
grd=0.5x0.5
lmask=nomask
INGRID[nfile]=$OUTGRID
QUERYFIL="$QUERY -var scripgriddata -res $grd -options lmask=$lmask"
if [ "$verbose" = "YES" ]; then
   echo $QUERYFIL
fi
OUTGRID[nfile]=`$QUERYFIL`
OUTFILE[nfile]=map_${res}_${lmask}_to_${grd}_${lmask}_aave_da_$CDATE.nc
if [ "$list" = "YES" ]; then
   echo "ingrid = ${OUTGRID[nfile]}"
   echo "ingrid = ${OUTGRID[nfile]}" >> $outfilelist
   echo "Succesffully found and listed all required mapping files"
   exit 0
fi
nfile=nfile+1

#----------------------------------------------------------------------
# Machine specific env stuff
#
# Bluefire specific environment settings to run MPI
#
hostname=`hostname`
case $hostname in
  ##bluefire
  be* )
  export TARGET_CPU_LIST="-1"

  export MP_RC_USE_LMC=yes

  export LAPI_DEBUG_ENABLE_AFFINITY=YES 
  export LAPI_DEBUG_BINDPROC_AFFINITY=YES 
  export LAPI_DEBUG_MTU_4K=YES
  export MP_SYNC_QP=YES
  export MP_RFIFO_SIZE=16777216
  export MP_SHM_ATTACH_THRESH=500000 # default is better sometimes
  export MP_EUIDEVELOP=min 
  export MP_LABELIO=yes

  export MP_SINGLE_THREAD=NO #do NOT use if more than one application 
                             #(including OpenMP) thread is making comm. calls

  export OMP_NUM_THREADS=1
  export XLSMPOPTS=stack=256000000
  export MP_PULSE=0
  export MP_USE_BULK_XFER=yes

  export MP_BULK_MIN_MSG_SIZE=64k
  export MP_RC_MAX_QP=8192
  export LAPI_DEBUG_RC_DREG_THRESHOLD=1000000
  export LAPI_DEBUG_QP_NOTIFICATION=no
  export LAPI_DEBUG_RC_INIT_SETUP=no	#try both = yes and = no

  if [ -z "$ESMFBIN_PATH" ]; then
     ESMFBIN_PATH=/contrib/esmf-5.2.0r-64-O/bin
  fi
  if [ -z "$MPIEXEC" ]; then
    MPIEXEC="mpirun.lsf"
  fi

  if [ "$interactive" != "no" ]; then
    # Bluefire specific commands to prepare to run interactively
    export MP_PROCS=$REGRID_PROC
    export MP_EUILIB=ip
                                                                                                                                                                          
    hostname > hostfile
    declare -i p=2
    until ((p>$MP_PROCS)); do
       hostname >> hostfile
       p=p+1
    done
    export MP_HOSTFILE=hostfile
  fi
  ;;

  ##jaguarpf
  jaguarpf* )
  if [ -z "$ESMFBIN_PATH" ]; then
     module load esmf/5.2.0-p1_with-netcdf_g
     ESMFBIN_PATH=$ESMF_BINDIR
  fi
  if [ -z "$MPIEXEC" ]; then
    MPIEXEC="aprun -n $REGRID_PROC"
  fi
  ;;

  *)
  echo "Machine $hostname NOT recognized"
  ;;

esac

if [ "$interactive" = "no" ]; then
  echo "Running in batch mode using MPI"
  if [ -z "$MPIEXEC" ]; then
     echo "Name of MPI exec to use was NOT set"
     echo "Set the environment variable: MPIEXEC"
     exit 1
  fi
  if [ -x "$MPIEXEC" ]; then
     echo "The MPIEXEC pathname given is NOT an executable: $MPIEXEC"
     echo "Set the environment variable: MPIEXEC or run in interactive mode without MPI"
     exit 1
  fi
  mpirun=$MPIEXEC
else
  mpirun=""
  echo "Running interactively"
fi
  

#----------------------------------------------------------------------

if [ ! -d "$ESMFBIN_PATH" ]; then
    echo "Path to ESMF binary directory does NOT exist: $ESMFBIN_PATH"
    echo "Set the environment variable: ESMFBIN_PATH"
    exit 1
fi

ESMF_REGRID="$ESMFBIN_PATH/ESMF_RegridWeightGen"

if [ ! -x "$ESMF_REGRID" ]; then
    echo "ESMF_RegridWeightGen does NOT exist in ESMF binary directory: $ESMFBIN_PATH\n"
    echo "Upgrade to a newer version of ESMF with this utility included"
    echo "Set the environment variable: ESMFBIN_PATH"
    exit 1
fi


# Remove previous log files
rm PET*.Log


#
# Now run the mapping for each file, checking that input files exist
# and then afterwards that the output mapping file exists
#
declare -i nfile=1
until ((nfile>${#INGRID[*]})); do
   echo "Creating ${OUTFILE[nfile]}"
   echo "From input grid: ${INGRID[nfile]}"
   echo "For output grid: ${OUTGRID[nfile]}"
   if [ -z "${INGRID[nfile]}" ] || [ -z "${OUTGRID[nfile]}" ] || [ -z "${OUTFILE[nfile]}" ]; then
      echo "Either input or output grid or output mapping file is NOT set"
      exit 3
   fi
   if [ ! -f "${INGRID[nfile]}" ]; then
      echo "Input grid file does NOT exist: ${INGRID[nfile]}"
      exit 2
   fi
   if [ ! -f "${OUTGRID[nfile]}" ]; then
      echo "Output grid file does NOT exist: ${OUTGRID[nfile]}"
      exit 3
   fi
   cmd="$mpirun $ESMF_REGRID --ignore_unmapped -s ${INGRID[nfile]} "
   cmd="$cmd -d ${OUTGRID[nfile]} -m conserve -w ${OUTFILE[nfile]}"
   runcmd $cmd

   if [ "$debug" != "YES" ] && [ ! -f "${OUTFILE[nfile]}" ]; then
      echo "Output mapping file was NOT created: ${OUTFILE[nfile]}"
      exit 4
   fi
   # add some metadata to the file
   HOST=`hostname`
   history="$ESMF_REGRID"
   runcmd "ncatted -a history,global,a,c,"$history"  ${OUTFILE[nfile]}"
   runcmd "ncatted -a hostname,global,a,c,$HOST   -h ${OUTFILE[nfile]}"
   runcmd "ncatted -a logname,global,a,c,$LOGNAME -h ${OUTFILE[nfile]}"

   # check for duplicate mapping weights
   newfile="rmdups_${OUTFILE[nfile]}"
   runcmd "rm -f $newfile"
   runcmd "env MAPFILE=${OUTFILE[nfile]} NEWMAPFILE=$newfile ncl rmdups.ncl"
   if [ -f "$newfile" ]; then
      runcmd "mv $newfile ${OUTFILE[nfile]}"
   fi

   nfile=nfile+1
done

echo "Successffully created needed mapping files for $res"

exit 0
