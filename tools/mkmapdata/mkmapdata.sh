#!/bin/bash
set -e
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
# -f <scripfilename> Input grid filename 
# -t <type>          Output type, supported values are [regional, global]
# -r <res>           Output resolution
# -b                 use batch mode (not default)
# -i                 High resolution mode (Only used with -f)
# -l                 list mapping files required (so can use check_input_data to get them)
# -d                 debug usage -- display mkmapdata that will be run but don't execute them
# -v                 verbose usage -- log more information on what is happening
# -h                 displays this help message
#
# You can also set the following env variables:
#
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
default_res="10x15"

#----------------------------------------------------------------------
# SET SOME DEFAULTS -- if not set via env variables outside

hostname=`hostname`
case $hostname in

  ##cheyenne
  cheyenne* | r* )
  if [ -z "$CSMDATA" ]; then
     CSMDATA=/glade/p/cesm/cseg/inputdata
  fi
  ;;

  ##casper
  casper* | crthc* )
  if [ -z "$CSMDATA" ]; then
     CSMDATA=/glade/p/cesm/cseg/inputdata
  fi
  ;;

  ##hobart/izumi/thorodin
  hobart* | izumi* | thorodin* )
  if [ -z "$CSMDATA" ]; then
     CSMDATA=/fs/cgd/csm/inputdata
  fi
  ;;

esac

if [[ -z ${CSMDATA} ]]; then
    echo "CSMDATA path not known for host ${hostname}. Set manually before calling mkmapdata.sh. E.g., bash: export CSMDATA=/path/to/csmdata"
    exit 7
elif [[ ! -d "${CSMDATA}" ]]; then
    echo "CSMDATA not found: ${CSMDATA}"
    exit 8
fi

#----------------------------------------------------------------------
# Usage subroutine
usage() {
  echo ""
  echo "**********************"
  echo "usage:"
  echo "./mkmapdata.sh"
  echo ""
  echo "valid arguments: "
  echo "[-f|--gridfile <gridname>] "
  echo "     Full pathname of model SCRIP grid file to use "
  echo "     This variable should be set if this is not a supported grid" 
  echo "     This variable will override the automatic generation of the"
  echo "     filename generated from the -res argument "
  echo "     the filename is generated ASSUMING that this is a supported "
  echo "     grid that has entries in the file namelist_defaults_ctsm.xml"
  echo "     the -r|--res argument MUST be specied if this argument is specified" 
  echo "[-r|--res <res>]"
  echo "     Model output resolution (default is $default_res)"
  echo "[-t|--gridtype <type>]"
  echo "     Model output grid type"
  echo "     supported values are [regional,global], (default is global)"
  echo "[-b|--batch]"
  echo "     Toggles batch mode usage (and run with mpi). If you want to run in batch mode"
  echo "     you need to have a separate batch script for a supported machine"
  echo "     that calls this script interactively - you cannot submit this"
  echo "     script directly to the batch system"
  echo "[-i|--hires]"
  echo "     Output maps are high resolution and large file support should be used"
  echo "[-l|--list]"
  echo "     List mapping files required (use check_input_data to get them)"
  echo "     also writes data to $outfilelist"
  echo "[-d|--debug]"
  echo "     Toggles debug-only (don't actually run mkmapdata just echo what would happen)"
  echo "[-h|--help]  "
  echo "     Displays this help message"
  echo "[-v|--verbose]"
  echo "     Toggle verbose usage -- log more information on what is happening "
  echo "[--fast]"
  echo "     Toggle fast maps only -- only create the maps that can be done quickly "
  echo ""
  echo " You can also set the following env variables:"
  echo "  ESMFBIN_PATH - Path to ESMF binaries "
  echo "                 (default is determined by machine running on)"
  echo "  CSMDATA ------ Path to CESM input data"
  echo "                 (default is $CSMDATA)"
  echo "  MPIEXEC ------ Name of mpirun executable"
  echo "                 (default is determined by machine running on)"
  echo "  REGRID_PROC -- Number of MPI processors to use"
  echo "                 (default is $REGRID_PROC)"
  echo ""
  echo "**defaults can be determined on the machines: cheyenne or casper"
  echo ""
  echo "**pass environment variables by preceding above commands "
  echo "  with 'env var1=setting var2=setting '"
  echo "**********************"
}
#----------------------------------------------------------------------
# runcmd subroutine
#----------------------------------------------------------------------

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
undo
   fi
   return 0
}

#----------------------------------------------------------------------
# Process input arguments
#----------------------------------------------------------------------

interactive="YES"
debug="no"
res="default"
type="global"
phys="clm4_5"
verbose="no"
list="no"
outgrid=""
gridfile="default"
fast="no"
netcdfout="none"

while [ $# -gt 0 ]; do
   case $1 in
       -v|-V)
	   verbose="YES"
	   ;;
       -b|--batch) 
	   interactive="NO"
	   ;;
       -d|--debug)
	   debug="YES"
	   ;;
       --fast)
	   fast="YES"
	   ;;
       -i|--hires)
           netcdfout="64bit_offset"
	   ;;
       -l|--list)
	   debug="YES"
	   list="YES"
	   ;;
       -r|--res)
	   res=$2
	   shift
	   ;;
       -f|--gridfile)
	   gridfile=$2
	   shift
	   ;;
       -t|--gridtype)
	   type=$2
	   shift
	   ;;
       -h|--help )
	   usage
	   exit 0
	   ;;
       * )
	   echo "ERROR:: invalid argument sent in: $2"
	   usage
	   exit 1
	   ;;
   esac
   shift
done

echo "Script to create mapping files required by mksurfdata_map"

#----------------------------------------------------------------------
# Determine output scrip grid file
#----------------------------------------------------------------------

# Set general query command used below
QUERY="$dir/../../bld/queryDefaultNamelist.pl -silent -namelist clmexp "
QUERY="$QUERY -justvalue -options sim_year=2000 -csmdata $CSMDATA"
echo "query command is $QUERY"

echo ""
DST_EXTRA_ARGS=""
if [ "$gridfile" != "default" ]; then
    GRIDFILE=$gridfile
    echo "Using user specified scrip grid file: $GRIDFILE" 
    if [ "$res" = "default" ]; then
       echo "When user specified grid file is given you MUST set the resolution (as the name of your grid)\n";
       exit 1
    fi
    
    # For now, maked the assumption about user-specified grids --
    # that they are SCRIP format. In the future we may want to 
    # provide a command-line options to allow the user to
    # override that default.
    DST_LRGFIL=$netcdfout
    DST_TYPE="SCRIP"
else
    if [ "$res" = "default" ]; then
       res=$default_res
    fi

    QUERYARGS="-res $res -options lmask=nomask"

    # Find the output grid file for this resolution using the XML database
    QUERYFIL="$QUERY -var scripgriddata $QUERYARGS -onlyfiles"
    if [ "$verbose" = "YES" ]; then
	echo $QUERYFIL
    fi
    GRIDFILE=`$QUERYFIL`
    echo "Using default scrip grid file: $GRIDFILE" 

    # Determine extra information about the destination grid file
    DST_LRGFIL=`$QUERY -var scripgriddata_lrgfile_needed $QUERYARGS`
    DST_TYPE=`$QUERY -var scripgriddata_type $QUERYARGS`
    if [ "$DST_TYPE" = "UGRID" ]; then
        # For UGRID, we need extra information: the meshname variable
	dst_meshname=`$QUERY -var scripgriddata_meshname $QUERYARGS`
	DST_EXTRA_ARGS="$DST_EXTRA_ARGS --dst_meshname $dst_meshname"
    fi
fi

if [ "$type" = "global" ] && [ `echo "$res" | grep -c "1x1_"` = 1 ]; then
   echo "This is a regional resolution and yet it is being run as global, set type with '-t' option\n";
   exit 1
fi
if [ "$type" = "global" ] && [ `echo "$res" | grep -c "5x5_"` = 1 ]; then
   echo "This is a regional resolution and yet it is being run as global, set type with '-t' option\n";
   exit 1
fi
echo "Output grid resolution is $res"
if [ -z "$GRIDFILE" ]; then
   echo "Output grid file was NOT found for this resolution: $res\n";
   exit 1
fi

if [ "$list" = "YES" ]; then
   echo "outgrid = $GRIDFILE"
   echo "outgrid = $GRIDFILE" > $outfilelist
elif [ ! -f "$GRIDFILE" ]; then
   echo "Input SCRIP grid file does NOT exist: $GRIDFILE\n";
   echo "Make sure CSMDATA environment variable is set correctly"
   exit 1
fi

#----------------------------------------------------------------------
# Determine all input grid files and output file names 
#----------------------------------------------------------------------

if [ "$phys" = "clm4_5" ]; then
    grids=(                    \
           "0.5x0.5_nomask"     \
           "0.25x0.25_nomask"   \
           "0.125x0.125_nomask"   \
           "3x3min_nomask" \
           "5x5min_nomask"     \
           "10x10min_nomask"   \
           "0.9x1.25_nomask" \
           "1km-merge-10min_HYDRO1K-merge-nomask" \
          )

else
    echo "ERROR: Unknown value for phys: $phys"
    exit 1
fi

# Set timestamp for names below 
# The flag `-d "-0 days"` can serve as a time saver as follows:
# If the script aborted without creating all of the map_ files and
# the user resubmits to create the remaining files on the next day,
# the user could change -0 to -1 to prevent the script from
# duplicating files already generated the day before.
# 
CDATE="c"`date -d "-0 days" +%y%m%d`

# Set name of each output mapping file
# First determine the name of the input scrip grid file  
# for each of the above grids
declare -i nfile=1
for gridmask in ${grids[*]}
do
   grid=${gridmask%_*}
   lmask=${gridmask#*_}

   QUERYARGS="-res $grid -options lmask=$lmask,glc_nec=10 "

   QUERYFIL="$QUERY -var scripgriddata $QUERYARGS -onlyfiles"
   if [ "$verbose" = "YES" ]; then
      echo $QUERYFIL
   fi
   INGRID[nfile]=`$QUERYFIL`
   if [ "$list" = "YES" ]; then
      echo "ingrid = ${INGRID[nfile]}"
      echo "ingrid = ${INGRID[nfile]}" >> $outfilelist
   fi

   OUTFILE[nfile]=map_${grid}_${lmask}_to_${res}_nomask_aave_da_$CDATE.nc

   # Determine extra information about the source grid file
   SRC_EXTRA_ARGS[nfile]=""
   SRC_LRGFIL[nfile]=`$QUERY -var scripgriddata_lrgfile_needed $QUERYARGS`
   SRC_TYPE[nfile]=`$QUERY -var scripgriddata_type $QUERYARGS`
   if [ "${SRC_TYPE[nfile]}" = "UGRID" ]; then
       # For UGRID, we need extra information: the meshname variable
       src_meshname=`$QUERY -var scripgriddata_meshname $QUERYARGS`
       SRC_EXTRA_ARGS[nfile]="${SRC_EXTRA_ARGS[nfile]} --src_meshname $src_meshname"
   fi

   nfile=nfile+1
done

#----------------------------------------------------------------------
# Determine supported machine specific stuff
#----------------------------------------------------------------------

if [ -n "$NERSC_HOST" ]; then
   hostname=$NERSC_HOST
fi
echo "Hostname = $hostname"
case $hostname in
  ##cheyenne
  cheyenne* | r* )
  . /glade/u/apps/ch/opt/lmod/8.1.7/lmod/lmod/init/bash
  if [ -z "$REGRID_PROC" ]; then
     REGRID_PROC=36
  fi
  if [ interactive = "YES" ]; then
     REGRID_PROC=1
  fi
  if [ "$verbose" = "YES" ]; then
     echo "Number of processors to regrid with = $REGRID_PROC"
  fi
  esmfvers=8.2.0b13
  intelvers=19.1.1
  module purge
  module load intel/$intelvers
# module load esmf_libs
# module load esmf_libs/$esmfvers
  module load nco

  if [[ $REGRID_PROC > 1 ]]; then
     mpi=mpt
     module load mpt/2.22
  else
     mpi=mpiuni
  fi
# module load esmf-${esmfvers}-ncdfio-${mpi}-O
  module use /glade/p/cesmdata/cseg/PROGS/modulefiles/esmfpkgs/intel/$intelvers
  module load esmf-${esmfvers}-ncdfio-${mpi}-g
  if [ -z "$ESMFBIN_PATH" ]; then
     ESMFBIN_PATH=`grep ESMF_APPSDIR $ESMFMKFILE | awk -F= '{print $2}'`
  fi
  if [ -z "$MPIEXEC" ]; then
     MPIEXEC="mpiexec_mpt -np $REGRID_PROC"
  fi
  if [ "$verbose" = "YES" ]; then
     echo "list of modules"
     module list
  fi
  ;;

  ## Casper
  casper* | crthc* )
  . /glade/u/apps/dav/opt/lmod/8.1.7/lmod/8.1.7/init/bash
  if [ -z "$REGRID_PROC" ]; then
     REGRID_PROC=8
  fi
  if [ interactive = "YES" ]; then
     REGRID_PROC=1
  fi
  echo "REGRID_PROC=$REGRID_PROC"
  esmfvers=7.1.0r
  intelvers=17.0.1
  module purge
  module load intel/$intelvers
  if [ $? != 0 ]; then
    echo "Error doing module load: intel/$intelvers"
    exit 1
  fi
  module load nco
  module load netcdf
  module load ncarcompilers

  module load esmflibs/$esmfvers
  if [ $? != 0 ]; then
    echo "Error doing module load: esmflibs/$esmfvers"
    exit 1
  fi

  if [[ $REGRID_PROC > 1 ]]; then
     mpi=mpi
     echo "MPI option is NOT currently available"
     exit 1
  else
     mpi=uni
  fi
  module load esmf-${esmfvers}-ncdfio-${mpi}-O
  if [ $? != 0 ]; then
    echo "Error doing module load: esmf-${esmfvers}-ncdfio-${mpi}-O"
    exit 1
  fi
  if [ -z "$ESMFBIN_PATH" ]; then
     ESMFBIN_PATH=`grep ESMF_APPSDIR $ESMFMKFILE | awk -F= '{print $2}'`
  fi
  echo "ESMFMKFILE: $ESMFMKFILE"
  echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

  if [ -z "$MPIEXEC" ]; then
     MPIEXEC="mpiexec -n $REGRID_PROC"
  fi
  ;;

  ##no other machine currently supported    
  *)
  echo "Machine $hostname NOT recognized"
  ;;

esac

# Error checks
if [ ! -d "$ESMFBIN_PATH" ]; then
    echo "Path to ESMF binary directory does NOT exist: $ESMFBIN_PATH"
    echo "Set the environment variable: ESMFBIN_PATH"
    exit 1
fi

#----------------------------------------------------------------------
# Generate the mapping files needed for surface dataset generation
#----------------------------------------------------------------------
 
# Resolve interactive or batch mode command
# NOTE - if you want to run in batch mode - you need to have a separate
# batch file that calls this script interactively - you cannot submit
# this script to the batch system

if [ "$interactive" = "NO" ]; then
   echo "Running in batch mode using MPI"
   if [ -z "$MPIEXEC" ]; then
      echo "Name of MPI exec to use was NOT set"
      echo "Set the environment variable: MPIEXEC"
      exit 1
   fi
   if [ ! -x `which ${MPIEXEC%% *}` ]; then
      echo "The MPIEXEC pathname given is NOT an executable: ${MPIEXEC%% *}"
      echo "Set the environment variable: MPIEXEC or run in interactive mode without MPI"
      exit 1
   fi
   mpirun=$MPIEXEC
   echo "Running in batch mode"
else
   mpirun=""
fi
  
ESMF_REGRID="$ESMFBIN_PATH/ESMF_RegridWeightGen"
if [ ! -x "$ESMF_REGRID" ]; then
    echo "ESMF_RegridWeightGen does NOT exist in ESMF binary directory: $ESMFBIN_PATH\n"
    echo "Upgrade to a newer version of ESMF with this utility included"
    echo "Set the environment variable: ESMFBIN_PATH"
    exit 1
fi

# Remove previous log files, if any
rm PET*.Log ||:

#
# Now run the mapping for each file, checking that input files exist
# and then afterwards that the output mapping file exists
#
declare -i nfile=1
until ((nfile>${#INGRID[*]})); do
   echo "Creating mapping file: ${OUTFILE[nfile]}"
   echo "From input grid: ${INGRID[nfile]}"
   echo "For output grid: $GRIDFILE"
   echo " "
   if [ -z "${INGRID[nfile]}" ] || [ -z "$GRIDFILE" ] || [ -z "${OUTFILE[nfile]}" ]; then
      echo "Either input or output grid or output mapping file is NOT set"
      exit 3
   fi
   if [ ! -f "${INGRID[nfile]}" ]; then
      echo "Input grid file does NOT exist: ${INGRID[nfile]}"
      if [ ! "$list" = "YES" ]; then
         exit 2
      fi
   fi
   if [ ! -f "$GRIDFILE" ]; then
      echo "Output grid file does NOT exist: $GRIDFILE"
      exit 3
   fi

   # Determine what (if any) large file support is needed. Use the
   # most extreme large file support needed by either the source file
   # or the destination file.
   if [ "$DST_LRGFIL" = "netcdf4" ] || [ "${SRC_LRGFIL[nfile]}" = "netcdf4" ]; then
       lrgfil="--netcdf4"
   elif [ "$DST_LRGFIL" = "64bit_offset" ] || [ "${SRC_LRGFIL[nfile]}" = "64bit_offset" ]; then
       lrgfil="--64bit_offset"
   elif [ "$DST_LRGFIL" = "none" ] && [ "${SRC_LRGFIL[nfile]}" = "none" ]; then
       lrgfil=""
   else
       echo "Unknown LRGFIL type:"
       echo "DST_LRGFIL = $DST_LRGFIL"
       echo "SRC_LRGFIL = ${SRC_LRGFIL[nfile]}"
       exit 4
   fi

   # Skip if file already exists
   if [ -f "${OUTFILE[nfile]}" ]; then
      echo "Skipping creation of ${OUTFILE[nfile]} as already exists"
   # Skip if large file and Fast mode is on
   elif [ "$fast" = "YES" ] && [ "${SRC_LRGFIL[nfile]}" = "netcdf4" ]; then
      echo "Skipping creation of ${OUTFILE[nfile]} as fast mode is on so skipping large files in NetCDF4 format"
   else

      cmd="$mpirun $ESMF_REGRID --ignore_unmapped -s ${INGRID[nfile]} "
      cmd="$cmd -d $GRIDFILE -m conserve -w ${OUTFILE[nfile]}"
      if [ $type = "regional" ]; then
        cmd="$cmd --dst_regional"
      fi

      cmd="$cmd --src_type ${SRC_TYPE[nfile]} ${SRC_EXTRA_ARGS[nfile]} --dst_type $DST_TYPE $DST_EXTRA_ARGS"
      cmd="$cmd $lrgfil"

      runcmd $cmd

      if [ "$debug" != "YES" ] && [ ! -f "${OUTFILE[nfile]}" ]; then
         echo "Output mapping file was NOT created: ${OUTFILE[nfile]}"
         exit 6
      fi
      # add some metadata to the file
      HOST=`hostname`
      history="$ESMF_REGRID"
      runcmd "ncatted -a history,global,a,c,"$history"  ${OUTFILE[nfile]}"
      runcmd "ncatted -a hostname,global,a,c,$HOST   -h ${OUTFILE[nfile]}"
      runcmd "ncatted -a logname,global,a,c,$LOGNAME -h ${OUTFILE[nfile]}"
   fi

   nfile=nfile+1
done

echo "Successfully created needed mapping files for $res"

exit 0
