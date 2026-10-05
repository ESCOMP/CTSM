#!/usr/bin/env bash
#PBS -N monthly_clim
#PBS  -r n 
#PBS  -j oe
#PBS  -k eod
#PBS  -S /bin/bash  
#PBS  -l select=1:ncpus=128:mpiprocs=128:ompthreads=1:mem=235GB
#PBS  -l walltime=04:00:00
#PBS -A P93300606
#PBS -q main
#
# ==============================================================================================
#
# Build a 12-month climatology from daily CESM coupler history files:
#   $ARCHDIR/$CASE/cpl/hist/$CASE.cpl.hx.atm.24h.avrg.YYYY-MM-DD-00000.nc
# for the years BEGYEAR to ENDYEAR
# Step 1: monthly mean for each year (ncra over the daily files of that month),
#         for the default variable list, a list given with -v, or all variables with -a
# Step 2: average each calendar month over all years (ncra), then ncrcat to 12 records
# Step 3: reset the time axis to mid-month of YEAR (noleap calendar)
#
# ==============================================================================================


# NOTE: This will fail on errors, or if variables are used but not set, and pipefail fails on first exit in a pipe of commands
set -euo pipefail

default_vars="atmImp_Faxa_bcph2,atmImp_Faxa_bcph3,atmImp_Faxa_ocph1,atmImp_Faxa_ocph2,atmImp_Faxa_ocph3,atmImp_Faxa_dstwet1,atmImp_Faxa_dstdry1,atmImp_Faxa_dstwet2,atmImp_Faxa_dstdry2,atmImp_Faxa_dstwet3,atmImp_Faxa_dstdry3,atmImp_Faxa_dstwet4,atmImp_Faxa_dstdry4"
# Grid fields that are always written to the output, whatever variable list is used
always_vars="atmImp_lon,atmImp_lat"
default_year="1850"
default_outdir="ARCHDIR"
# Number of simultaneous ncra jobs: 16 inside a PBS batch job, 1 otherwise
if [[ -n ${PBS_JOBID:-} ]]; then
  default_njobs=16
else
  default_njobs=1
fi

# Print the help: to stdout and exit 0 for -h, to stderr and exit 1 for bad arguments
usage() {
  local status=${1:-1}
  if (( status == 0 )); then exec 3>&1; else exec 3>&2; fi
  cat >&3 <<EOF
Usage: $0 [-b [-A ACCOUNT]] [-k] [-h] [-j NJOBS] [-o DIR] [-y YEAR] [-v VARLIST | -a] ARCHDIR CASE BEGYEAR ENDYEAR
  BEGYEAR ENDYEAR  first and last year of the daily files to average (inclusive)
  -a          All variables in the files
  -A ACCOUNT  PBS project account to charge with -b (default: the #PBS -A line in this script)
  -b          Batch submission. Check the input, then submit this script to the PBS batch queue with the same
              arguments (uses the #PBS settings at the top of this script)
  -h          Exit with this help
  -j NJOBS    Number of ncra jobs to run at the same time in the background
              (default: 16 in a PBS batch job, 1 otherwise; currently $default_njobs)
  -k          Keep the intermediate files in the output CASE/proc/mon directory
  -o DIR      Output directory if different from ARCHDIR
              Files will go under DIR/CASE/proc
              (default: $default_outdir)
  -v VARLIST  Comma-separated list of variables to average
              (default: $default_vars)
              $always_vars are always included
  -y YEAR     Year to set the climatology to
              (default: $default_year)
EOF
  exit "$status"
}

keep=0
varlist=$default_vars
allvars=0
userlist=0
year=$default_year
outdir=$default_outdir
njobs=$default_njobs
batch=0
qsub_opts=()

# Inside a job submitted with -b: restore the original arguments saved by the submitting script
# (qsub can't pass arguments to a job script, so they go through a file named in the environment)
if [[ -n ${MONTHLY_CLIM_ARGFILE:-} ]]; then
  cd "${MONTHLY_CLIM_WORKDIR:-.}"
  eval "set -- $(cat "$MONTHLY_CLIM_ARGFILE")"
  rm -f "$MONTHLY_CLIM_ARGFILE"
fi

# Arguments to pass on to the batch job: everything except -b and -A
job_args=()
while getopts "bA:kj:v:o:y:ah" opt; do
  case $opt in
    b) batch=1 ;;
    A) qsub_opts+=(-A "$OPTARG") ;;
    k) keep=1;                       job_args+=(-k) ;;
    j) njobs=$OPTARG;                job_args+=(-j "$OPTARG") ;;
    v) varlist=$OPTARG; userlist=1;  job_args+=(-v "$OPTARG") ;;
    o) outdir=$OPTARG;               job_args+=(-o "$OPTARG") ;;
    y) year=$OPTARG;                 job_args+=(-y "$OPTARG") ;;
    a) allvars=1;                    job_args+=(-a) ;;
    h) usage 0 ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))
[[ $# -eq 4 ]] || usage
archdir=$1; case=$2; begyr=$3; endyr=$4
job_args+=("$@")
if (( ${#qsub_opts[@]} > 0 && ! batch )); then
    echo "-A is only used with -b" >&2
    exit 1
fi

# Error check input
for y in "$begyr" "$endyr"; do
  if ! [[ $y =~ ^[0-9]{1,4}$ ]]; then
      echo "BEGYEAR and ENDYEAR must be integers from 0 to 9999: $y" >&2
      exit 1
  fi
done
begyr=$((10#$begyr)); endyr=$((10#$endyr))
if (( begyr > endyr )); then
    echo "BEGYEAR ($begyr) must not be after ENDYEAR ($endyr)" >&2
    exit 1
fi
if ! [[ $njobs =~ ^[1-9][0-9]*$ ]]; then
    echo "NJOBS must be a positive integer: $njobs" >&2
    exit 1
fi
if (( allvars && userlist )); then
    echo "Use either -v or -a, not both" >&2
    exit 1
fi
# Variable subset passed to ncra (coordinates and time_bnds are included automatically)
if (( allvars )); then
    vopt=()
    echo "Averaging all variables"
else
    [[ -n $varlist ]] || { echo "Variable list given with -v is empty" >&2; exit 1; }
    # Add the grid fields that are always output, dropping any duplicates
    IFS=',' read -ra vl <<< "$varlist,$always_vars"
    varlist=$(printf '%s\n' "${vl[@]}" | awk 'NF && !seen[$0]++' | paste -sd, -)
    vopt=(-v "$varlist")
    echo "Averaging variables: $varlist"
fi
if ! [[ $year =~ ^[0-9]{1,4}$ ]]; then
    echo "YEAR must be an integer from 0 to 9999: $year" >&2
    exit 1
fi
year=$(printf '%04d' "$((10#$year))")
if [ ! -d "$archdir" ]; then
    echo "Input $archdir directory does not exist" >&2
    exit 1
fi
indir="$archdir/$case/cpl/hist"
if [ "$outdir" = "$default_outdir" ]; then
    outdir="$archdir/$case/proc"
elif [ ! -d "$outdir" ]; then
    echo "Output $outdir directory does not exist" >&2
    exit 1
else
    outdir="$outdir/$case/proc"
fi
if [ ! -d "$indir" ]; then
    echo "$indir directory does not exist" >&2
    exit 1
fi
mkdir -p "$outdir/mon" || { echo "Can not create $outdir/mon" >&2; exit 1; }

echo "Create a monthly climatology from a list of years from daily averaged CPLHIST files for DATM"

if command -v module >/dev/null 2>&1; then
   # Lmod's module function can reference unset variables, so relax set -u around it
   set +u
   module load nco
   set -u
else
   echo "module command not found, so expecting NCO to already be in the PATH" >&2
fi
for cmd in ncra ncrcat ncap2 ncatted ncks; do
  command -v $cmd >/dev/null || { echo "$cmd not found (module load nco)" >&2; exit 1; }
done
prefix="$case.cpl.hx.atm.24h.avrg"
noleap_days=(31 28 31 30 31 30 31 31 30 31 30 31)

# Year-months to process, BEGYEAR-01 to ENDYEAR-12 (files are matched by name, not the time coordinate)
yms=()
for (( y = begyr; y <= endyr; y++ )); do
  for mm in 01 02 03 04 05 06 07 08 09 10 11 12; do
    yms+=("$(printf '%04d' "$y")-$mm")
  done
done
echo "Processing ${#yms[@]} year-months: ${yms[0]} .. ${yms[-1]}"

# Check that every requested variable is in the first file before doing any work
if (( ! allvars )); then
  echo "Checking that all requested variables are present in the first file for ${yms[0]}"
  first=$(find "$indir" -maxdepth 1 -name "${prefix}.${yms[0]}-??-00000.nc" | sort | sed -n 1p)
  [[ -n $first ]] || { echo "ERROR: no files for ${yms[0]} matching ${prefix}.${yms[0]}-??-00000.nc in $indir" >&2; exit 1; }
  missing=()
  IFS=',' read -ra vars <<< "$varlist"
  for v in "${vars[@]}"; do
    ncks -m -v "$v" "$first" >/dev/null 2>&1 || missing+=("$v")
  done
  if (( ${#missing[@]} > 0 )); then
    echo "ERROR: variables not found in $first: ${missing[*]}" >&2
    exit 1
  fi
fi

# Print the year at January of BEGYEAR, BEGYEAR+10, BEGYEAR+20, ... (every 120 year-months)
progress() {
  local y=$((10#${1%%-*}))
  if [[ ${1#*-} == 01 ]] && (( (y - begyr) % 10 == 0 )); then
    echo "  at year $y"
  fi
}

# Run ncra in the background, with at most $njobs running at once
# Usage: bg_ncra OUTFILE INFILE...  (uses the ncra_opts array for extra ncra options)
# The input list goes to ncra on stdin, so a long list never hits the shell argument limit.
# A failed job leaves the $failed marker file, which wait_ncra checks.
failed="$outdir/mon/.ncra_failed"
rm -f "$failed"
bg_ncra() {
  local out=$1; shift
  while (( $(jobs -rp | wc -l) >= njobs )); do wait -n || true; done
  ( printf '%s\n' "$@" | ncra -O "${ncra_opts[@]}" -o "$out" \
      || { echo "ERROR: ncra failed creating $out" >&2; touch "$failed"; } ) &
}
# Wait for all background ncra jobs, and stop if any of them failed
wait_ncra() {
  wait
  if [[ -e $failed ]]; then
    echo "ERROR: one or more ncra jobs failed, stopping" >&2
    exit 1
  fi
}

# On Ctrl+C (INT), or when PBS kills the job with qdel or at the walltime limit (TERM),
# also stop the background ncra jobs. Background jobs in a script ignore Ctrl+C,
# so without this they would keep running after the script exits.
stop_on_signal() {
  local sig=$1 pid
  trap '' INT TERM
  echo "Caught SIG$sig, stopping background ncra jobs" >&2
  for pid in $(jobs -p); do
    pkill -TERM -P "$pid" 2>/dev/null || true   # the ncra started by this background job
    kill -TERM "$pid" 2>/dev/null || true
  done
  wait 2>/dev/null || true
  if [[ $sig == INT ]]; then exit 130; else exit 143; fi
}
trap 'stop_on_signal INT' INT
trap 'stop_on_signal TERM' TERM

# Check every month in the range has a full set of daily files before doing any averaging
bad=0
echo "Checking that every month has the expected number of daily files"
for ym in "${yms[@]}"; do
  mm=${ym#*-}
  progress "$ym"
  n=$(find "$indir" -maxdepth 1 -name "${prefix}.${ym}-??-00000.nc" | wc -l)
  expect=${noleap_days[$((10#$mm - 1))]}
  if (( n != expect )); then
    echo "ERROR: $ym has $n daily files (expected $expect)" >&2
    bad=1
  fi
done
(( ! bad )) || exit 1

# With -b: the input has been checked, so submit this script to PBS and stop here.
# The arguments are saved (shell-quoted, so commas and spaces survive) to a file that the job reads back.
if (( batch )); then
  command -v qsub >/dev/null || { echo "qsub not found, can not submit with -b" >&2; exit 1; }
  argfile="$outdir/monthly_clim.args.$$"
  printf '%q ' "${job_args[@]}" > "$argfile"
  echo "Submitting to PBS with arguments: ${job_args[*]}"
  if ! jobid=$(qsub "${qsub_opts[@]}" -v "MONTHLY_CLIM_ARGFILE=$argfile,MONTHLY_CLIM_WORKDIR=$PWD" "$0"); then
    rm -f "$argfile"
    echo "ERROR: qsub failed" >&2
    exit 1
  fi
  echo "Successfully submitted job $jobid (log file will be monthly_clim.o${jobid%%.*} in $PWD)"
  exit 0
fi

# ---- Step 1: monthly means for each year ----
echo "Computing monthly means for each year ($njobs ncra jobs at a time)"
ncra_opts=("${vopt[@]}")
for ym in "${yms[@]}"; do
  progress "$ym"
  mapfile -t files < <(find "$indir" -maxdepth 1 -name "${prefix}.${ym}-??-00000.nc" | sort)
  bg_ncra "$outdir/mon/mon_${ym}.nc" "${files[@]}"
done
wait_ncra

# ---- Step 2: average each calendar month over the years ----
echo "Averaging each calendar month over the years"
ncra_opts=()
for mm in 01 02 03 04 05 06 07 08 09 10 11 12; do
  # Only handle this run's years, so leftover files from an earlier -k run are not included
  mons=()
  for (( y = begyr; y <= endyr; y++ )); do
    mons+=("$outdir/mon/mon_$(printf '%04d' "$y")-${mm}.nc")
  done
  echo "Month $mm: averaging ${#mons[@]} years ($begyr-$endyr)"
  bg_ncra "$outdir/mon/clim_${mm}.nc" "${mons[@]}"
done
wait_ncra

echo "Put the climatology into one and file and finalize metadata on the file"
out="$outdir/${case}.cpl.hx.atm.monclim.$(printf '%04d-%04d' "$begyr" "$endyr")_c$(date +%Y%m%d).nc"
ncrcat -O "$outdir"/mon/clim_{01..12}.nc "$out"

# ---- Step 3: reset time to mid-month of $year (noleap), days since $year-01-01 ----
ncap2 -O -s 'time[$time]={15.5,45.0,74.5,105.0,135.5,166.0,196.5,227.5,258.0,288.5,319.0,349.5};' "$out" "$out"
ncatted -O -a units,time,o,c,"days since ${year}-01-01 00:00:00" \
           -a calendar,time,o,c,"noleap" "$out"
# If the file has time bounds, set them to the start and end of each month
if ncks --trd -m -v time_bnds "$out" >/dev/null 2>&1; then
  ncap2 -O -s '*lb[$time]={0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334.};
               *ub[$time]={31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334.,365.};
               time_bnds(:,0)=lb; time_bnds(:,1)=ub;' "$out" "$out"
  ncatted -O -a units,time_bnds,o,c,"days since ${year}-01-01 00:00:00" "$out"
fi
echo "Wrote $out (time = mid-month of year $year, noleap)"

# Remove the temporary monthly files unless -k was given
if (( keep )); then
  echo "Kept intermediate files in $outdir/mon"
else
  rm -r "$outdir/mon"
fi

echo ""
echo "Successfully created monthly climatology: $out"
