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
# -f <scripfilename> Input grid filename 
# -t <type>          Output type, supported values are [regional, global]
# -r <res>           Output resolution
# -b                 use batch mode (not default)
# -l                 list mapping files required (so can use check_input_data to get them)
# -d                 debug usage -- display mkmapdata that will be run but don't execute them
# -o <ogrid>         Also map to the input ocean-grid resolution"
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

if [ -z "$CSMDATA" ]; then
   CSMDATA=/glade/proj3/cseg/inputdata
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
  echo "[-f|--gridfile <gridname>] "
  echo "     Full pathname of model SCRIP grid file to use "
  echo "     This variable should be set if this is not a supported grid" 
  echo "     This variable will override the automatic generation of the"
  echo "     filename generated from the -res argument "
  echo "     the filename is generated ASSUMING that this is a supported "
  echo "     grid that has entries in the file namelist_defaults_clm.xml"
  echo "     the -r|--res argument MUST be specied if this argument is specified" 
  echo "[-r|--res <res>]"
  echo "     Model output resolution (default is $default_res)"
  echo "[-t|--gridtype <type>]"
  echo "     Model output grid type"
  echo "     supported values are [regional,global], (default is global)"
  echo "[-b|--batch]"
  echo "     Toggles batch mode usage. If you want to run in batch mode"
  echo "     you need to have a separate batch script for a supported machine"
  echo "     that calls this script interactively - you cannot submit this"
  echo "     script directory to the batch system"
  echo "[-l|--list]"
  echo "     List mapping files required (use check_input_data to get them)"
  echo "     also writes data to $outfilelist"
  echo "[-d|--debug]"
  echo "     Toggles debug-only (don't actually run mkmapdata just echo what would happen)"
  echo "[-o|--ocn] <ogrid>"
  echo "     Also map to the input ocean-grid resolution"
  echo "[-h|--help]  "
  echo "     Displays this help message"
  echo "[-v|--verbose]"
  echo "     Toggle verbose usage -- log more information on what is happening "
  echo ""
  echo " You can also set the following env variables:"
  echo "  ESMFBIN_PATH - Path to ESMF binaries "
  echo "                 (default is /contrib/esmf-5.3.0-64-O/bin)"
  echo "  CSMDATA ------ Path to CESM input data"
  echo "                 (default is /glade/proj3/cseg/inputdata)"
  echo "  MPIEXEC ------ Name of mpirun executable"
  echo "                 (default is mpirun.lsf)"
  echo "  REGRID_PROC -- Number of MPI processors to use"
  echo "                 (default is 8)"
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
verbose="no"
ocean="no"
list="no"
outgrid=""
gridfile="default"

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
       -l|--list)
	   debug="YES"
	   list="YES"
	   ;;
       -r|--res)
	   res=$2
	   shift
	   ;;
       -o|--ocn)
	   ocn=$2
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
QUERY="$dir/../../bld/queryDefaultNamelist.pl -silent -namelist clmexp -onlyfiles "
QUERY="$QUERY -justvalue -options sim_year=2000 -csmdata $CSMDATA"
echo "query command is $QUERY"

echo ""
if [ "$gridfile" != "default" ]; then
    GRIDFILE=$gridfile
    echo "Using user specified scrip grid file: $GRIDFILE" 
    if [ "$res" = "default" ]; then
       echo "When user specified grid file is given you MUST set the resolution (as the name of your grid)\n";
       exit 1
    fi
else
    if [ "$res" = "default" ]; then
       res=$default_res
    fi
    # Find the output grid file for this resolution using the XML database
    QUERYFIL="$QUERY -var scripgriddata -res $res -options lmask=nomask"
    if [ "$verbose" = "YES" ]; then
	echo $QUERYFIL
    fi
    GRIDFILE=`$QUERYFIL`
    echo "Using default scrip grid file: $GRIDFILE" 
fi
if [ "$type" = "global" ] && [ `echo "$res" | grep -c "1x1_"` = 1 ]; then
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

grids=(                    \
       "1.9x2.5_nomask"    \ 
       "0.5x0.5_USGS"      \
       "0.5x0.5_AVHRR"     \
       "0.5x0.5_MODIS"     \
       "3x3min_LandScan2004" \
       "3x3min_MODIS"      \
       "3x3min_USGS"       \
       "5x5min_nomask"     \
       "5x5min_IGBP-GSDP"  \
       "5x5min_ISRIC-WISE" \
       "10x10min_nomask"   \
       "3x3min_GLOBE-Gardner" \
       "3x3min_GLOBE-Gardner-mergeGIS" )

# Set timestamp for names below 
CDATE="c"`date +%y%m%d`

# Set name of each output mapping file OTHER than RTM
# First determine the name of the input scrip grid file  
# for each of the above grids
declare -i nfile=1
for gridmask in ${grids[*]}
do
   grid=${gridmask%_*}
   lmask=${gridmask#*_}
   QUERYFIL="$QUERY -var scripgriddata -res $grid -options lmask=$lmask"
   if [ "$verbose" = "YES" ]; then
      echo $QUERYFIL
   fi
   INGRID[nfile]=`$QUERYFIL`
   if [ "$list" = "YES" ]; then
      echo "ingrid = ${INGRID[nfile]}"
      echo "ingrid = ${INGRID[nfile]}" >> $outfilelist
   fi
   if [[ "$gridmask" = "3x3min_MODIS" || "$gridmask" = "3x3min_LandScan2004" || "$gridmask" = "3x3min_GLOBE-Gardner" || "$gridmask" = "3x3min_GLOBE-Gardner-mergeGIS" ]]; then
      LRGFIL[nfile]="yes"
   else
      LRGFIL[nfile]="no"
   fi
   GRIDFILE[nfile]=$GRIDFILE
   OUTFILE[nfile]=map_${grid}_${lmask}_to_${res}_nomask_aave_da_$CDATE.nc
   nfile=nfile+1
done

# Add main mapping file from ocean to atmosphere
if [ "$ocean" != "no" ]; then
   INGRID[nfile]=`$QUERY -var scripgriddata -options hgrid=$ocean,lmask=$ocean`
   OUTGRID[nfile]=$GRIDFILE
   OUTFILE[nfile]=map_${ocean}_to_${res}_aave_da_$CDATE.nc
   nfile=nfile+1
fi

# Set name of RTM output mapping file
if [ "$type" = "global" ]; then
   grid=0.5x0.5
   lmask=nomask
   INGRID[nfile]=$GRIDFILE
   QUERYFIL="$QUERY -var scripgriddata -res $grid -options lmask=$lmask"
   if [ "$verbose" = "YES" ]; then
      echo $QUERYFIL
   fi
   GRIDFILE[nfile]=`$QUERYFIL`
   OUTFILE[nfile]=map_${res}_${lmask}_to_${grid}_${lmask}_aave_da_$CDATE.nc
   if [ "$list" = "YES" ]; then
      echo "ingrid = ${GRIDFILE[nfile]}"
      echo "ingrid = ${GRIDFILE[nfile]}" >> $outfilelist
      echo "Succesffully found and listed all required mapping files"
      exit 0
   fi
   nfile=nfile+1
fi

#----------------------------------------------------------------------
# Determine supported machine specific stuff
#----------------------------------------------------------------------

hostname=`hostname`
case $hostname in
  ##bluefire
  be* )
  source /contrib/Modules/3.2.6/init/bash
  module load netcdf/4.1.3_seq

  if [ -z "$ESMFBIN_PATH" ]; then
     if [ "$type" = "global" ]; then
        ESMFBIN_PATH=/contrib/esmf-5.2.0rp1bs09-64/bin
     else
        ESMFBIN_PATH=/contrib/esmf-5.3.0-64-O-mpiuni-netcdf4.1.3/bin
     fi
  fi
  if [ -z "$MPIEXEC" ]; then
     MPIEXEC="mpirun.lsf"
  fi
  if [ "$interactive" != "NO" ]; then
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
  ## NOTE that for jaguarpf there is no batch script for now
  jaguarpf* )
  if [ -z "$ESMFBIN_PATH" ]; then
     module load esmf/5.2.0-p1_with-netcdf_g
     ESMFBIN_PATH=$ESMF_BINDIR
  fi
  if [ -z "$MPIEXEC" ]; then
    MPIEXEC="aprun -n $REGRID_PROC"
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
# and the RTM run time mapping
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
   if [ ! -x `which $MPIEXEC` ]; then
      echo "The MPIEXEC pathname given is NOT an executable: $MPIEXEC"
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

# Remove previous log files
rm PET*.Log

#
# Now run the mapping for each file, checking that input files exist
# and then afterwards that the output mapping file exists
#
declare -i nfile=1
until ((nfile>${#INGRID[*]})); do
   echo "Creating mapping file: ${OUTFILE[nfile]}"
   echo "From input grid: ${INGRID[nfile]}"
   echo "For output grid: ${GRIDFILE[nfile]}"
   echo " "
   if [ -z "${INGRID[nfile]}" ] || [ -z "${GRIDFILE[nfile]}" ] || [ -z "${OUTFILE[nfile]}" ]; then
      echo "Either input or output grid or output mapping file is NOT set"
      exit 3
   fi
   if [ ! -f "${INGRID[nfile]}" ]; then
      echo "Input grid file does NOT exist: ${INGRID[nfile]}"
      exit 2
   fi
   if [ ! -f "${GRIDFILE[nfile]}" ]; then
      echo "Output grid file does NOT exist: ${GRIDFILE[nfile]}"
      exit 3
   fi
   cmd="$mpirun $ESMF_REGRID --ignore_unmapped -s ${INGRID[nfile]} "
   cmd="$cmd -d ${GRIDFILE[nfile]} -m conserve -w ${OUTFILE[nfile]}"
   if [ $type = "regional" ]; then
     cmd="$cmd --dst_regional"
   fi

   if [ "${LRGFIL[nfile]}" = "yes" ]; then
      cmd="$cmd --64bit_offset "
   fi
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
