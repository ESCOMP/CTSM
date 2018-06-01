#!/bin/bash
#
#
# Batch script to submit to create mapping files for all standard
# resolutions.  If you provide a single resolution via "$RES", only
# that resolution will be used. In that case: If it is a regional or
# single point resolution, you should set '#SBATCH -n' to 1, and be sure
# that '-t regional' is specified in cmdargs.
#
# geyser specific batch commands:
#SBATCH -J regrid        # job name
#SBATCH -n 8
#SBATCH --ntasks-per-node=4
#SBATCH --mem=450G
#SBATCH -t 03:00:00
#SBATCH -A P93300606
#SBATCH -p dav
#SBATCH -e regrid.%J.out   # output filename
#SBATCH -o regrid.%J.err   # error filename

#----------------------------------------------------------------------
# Set parameters
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Begin main script
#----------------------------------------------------------------------

if [ -z "$RES" ]; then
   echo "Run for all valid resolutions"
   resols=`../../bld/queryDefaultNamelist.pl -res list -silent`
else
   resols="$RES"
fi
echo "Create mapping files for this list of resolutions: $resols"

#----------------------------------------------------------------------

for res in $resols; do
   echo "Create mapping files for: $res"
#----------------------------------------------------------------------
   cmdargs="-r $res"

   # For single-point and regional resolutions, tell mkmapdata that
   # output type is regional
   if [[ `echo "$res" | grep -c "1x1_"` -gt 0 || `echo "$res" | grep -c "5x5_"` -gt 0 ]]; then
       res_type="regional"
   else
       res_type="global"
   fi

   cmdargs="$cmdargs -t $res_type"

   echo "$res_type"
   if [ "$res_type" = "regional" ]; then
       echo "regional"
       # For regional and (especially) single-point grids, we can get
       # errors when trying to use multiple processors - so just use 1.
       # We also do NOT set batch mode in this case, because some
       # machines (e.g., yellowstone) do not listen to REGRID_PROC, so to
       # get a single processor, we need to run mkmapdata.sh in
       # interactive mode.
       regrid_num_proc=1
   else
       echo "global"
       regrid_num_proc=8
       if [ ! -z "$LSFUSER" ]; then
           echo "batch"
	   cmdargs="$cmdargs -b"
       fi
       if [ ! -z "$PBS_O_WORKDIR" ]; then
           cd $PBS_O_WORKDIR
	   cmdargs="$cmdargs -b"
       fi
   fi

   echo "args: $cmdargs"
   echo "time env REGRID_PROC=$regrid_num_proc ./mkmapdata.sh $cmdargs\n"
   time env REGRID_PROC=$regrid_num_proc ./mkmapdata.sh $cmdargs
done
