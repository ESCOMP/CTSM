"""
Create a 12-month climatology from daily CMEPS coupler history (CPLHIST) files, for use by DATM.

The input files are the daily averaged atmosphere coupler history files:

    ARCHDIR/CASE/cpl/hist/CASE.cpl.hx.atm.24h.avrg.YYYY-MM-DD-00000.nc

for the years BEGYEAR to ENDYEAR. The steps are:

  1. Average the daily files of each month of each year (ncra)
  2. Average each calendar month over all the years (ncra), and put the 12 months in one
     file (ncrcat)
  3. Set the time axis to the middle of each month of YEAR (noleap calendar)

The output file is:

    OUTDIR/CASE/proc/CASE.cpl.hx.atm.monclim.BEGYEAR-ENDYEAR_cYYYYMMDD.nc

This needs the NCO commands ncra, ncrcat and ncatted (on Derecho: module load nco).
With --batch, the input is checked and then this tool is submitted to the PBS batch queue.

A previous version of this for CESM2 is in NCL and here on Derecho by Keith Oleson:
    /glade/u/home/oleson/misc_programs/CLM5Dev/CMIP6/CreateAeroDepFile_Task5B.ncl
It did 5 year averages over the historical period and then did a linear interpolation
between the 5-year averages. It also used CPLHIST output from several ensemble members and
averged them together with ncea. This is something that will have to be done in the future
in this tool as well.

Another example for Nitrogen deposition is here from Simone Tilmes and Mike Mills:
   https://svn.code.sf.net/p/codescripts/code/trunk/ncl/cam_forcing/CreateDepositionFile.ncl

Neither were quite right for what's needed here. And both are in NCL which is well after end of life now.

Also there is code in LDF to create monthly climatologies, but it's not setup to bring in ensemble
averages and not specialized to create datafiles for DATM aerosol forcing.
"""

import argparse
import concurrent.futures  # Handles parallel multi-processing by sending ncra jobs on different threads
import datetime
import getpass
import io
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys

from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm.args_utils import comma_separated_list
from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
    log,
)
from ctsm.git_utils import get_ctsm_git_short_hash
from ctsm.os_utils import run_cmd_output_on_error
from ctsm.path_utils import path_to_ctsm_root
from ctsm.toolchain.gen_mksurfdata_jobscript_single import write_runscript_part1
from ctsm.utils import abort

logger = logging.getLogger(__name__)

# Variables averaged by default
DEFAULT_VARS = [
    "atmImp_Faxa_bcph2",
    "atmImp_Faxa_bcph3",
    "atmImp_Faxa_ocph1",
    "atmImp_Faxa_ocph2",
    "atmImp_Faxa_ocph3",
    "atmImp_Faxa_dstwet1",
    "atmImp_Faxa_dstdry1",
    "atmImp_Faxa_dstwet2",
    "atmImp_Faxa_dstdry2",
    "atmImp_Faxa_dstwet3",
    "atmImp_Faxa_dstdry3",
    "atmImp_Faxa_dstwet4",
    "atmImp_Faxa_dstdry4",
]
# Grid fields that are always written to the output, whatever variable list is used
ALWAYS_VARS = ["atmImp_lon", "atmImp_lat"]
DEFAULT_YEAR = 1850
# Days in each month of the noleap calendar
NOLEAP_DAYS = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
# File name tag of the daily averaged atmosphere coupler history files
HIST_TAG = "cpl.hx.atm.24h.avrg"
# NCO commands that are needed
NCO_COMMANDS = ("ncra", "ncrcat", "ncatted")
# Number of simultaneous ncra jobs inside a PBS batch job
PBS_NJOBS = 16
# PBS batch settings
PBS_JOBNAME = "daily_to_monthly_clim"
DEFAULT_ACCOUNT = "P93300606"
DEFAULT_QUEUE = "main"
DEFAULT_WALLTIME = "04:00:00"
# The executable wrapper for this module, run by the batch job
WRAPPER = os.path.join("tools", "cplhist", "daily_to_monthly_clim")


def get_parser():
    """
    Get the parser object for daily_to_monthly_clim
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "archdir", help="Archive directory with the case directories (ARCHDIR/CASE/cpl/hist)"
    )
    parser.add_argument("case", help="Case name")
    parser.add_argument(
        "begyear", type=int, help="First year of the daily files to average (inclusive)"
    )
    parser.add_argument(
        "endyear", type=int, help="Last year of the daily files to average (inclusive)"
    )

    vars_group = parser.add_mutually_exclusive_group()
    vars_group.add_argument(
        "--vars",
        help=f"""Comma-separated list of variables to average
            ({",".join(ALWAYS_VARS)} are always included)
            (default: {",".join(DEFAULT_VARS)})""",
        action="store",
        dest="vars",
        type=comma_separated_list,
        default=DEFAULT_VARS,
    )
    vars_group.add_argument(
        "-a",
        "--all-vars",
        help="Average all variables in the files",
        action="store_true",
        dest="all_vars",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="""Output directory if different from ARCHDIR. Files go under OUTDIR/CASE/proc
            (default: ARCHDIR)""",
        action="store",
        dest="outdir",
        default=None,
    )
    parser.add_argument(
        "-y",
        "--year",
        help="Year to set the climatology time axis to (default: %(default)s)",
        action="store",
        dest="year",
        type=int,
        default=DEFAULT_YEAR,
    )
    parser.add_argument(
        "-k",
        "--keep",
        help="Keep the intermediate files in the OUTDIR/CASE/proc/mon directory",
        action="store_true",
        dest="keep",
    )
    parser.add_argument(
        "-j",
        "--njobs",
        help=f"""Number of ncra jobs to run at the same time
            (default: {PBS_NJOBS} in a PBS batch job, 1 otherwise)""",
        action="store",
        dest="njobs",
        type=int,
        default=None,
    )
    parser.add_argument(
        "-b",
        "--batch",
        help="Check the input, then submit this tool to the PBS batch queue with the same options",
        action="store_true",
        dest="batch",
    )
    parser.add_argument(
        "-A",
        "--account",
        help=f"""PBS project account to charge with --batch
            (default: $ACCOUNT if set, otherwise {DEFAULT_ACCOUNT})""",
        action="store",
        dest="account",
        default=None,
    )
    parser.add_argument(
        "--queue",
        help="PBS queue to use with --batch (default: %(default)s)",
        action="store",
        dest="queue",
        default=DEFAULT_QUEUE,
    )
    parser.add_argument(
        "--walltime",
        help="PBS wall clock time limit to use with --batch (default: %(default)s)",
        action="store",
        dest="walltime",
        default=DEFAULT_WALLTIME,
    )
    add_logging_args(parser)
    return parser


def default_njobs():
    """
    Number of simultaneous ncra jobs when not given: more inside a PBS batch job
    """
    if os.environ.get("PBS_JOBID"):
        return PBS_NJOBS
    return 1


def input_dir(archdir, case):
    """
    Directory with the daily coupler history files
    """
    return os.path.join(archdir, case, "cpl", "hist")


def output_dir(archdir, outdir, case):
    """
    Directory for the output (OUTDIR/CASE/proc, where OUTDIR defaults to ARCHDIR)
    """
    if outdir is None:
        outdir = archdir
    return os.path.join(outdir, case, "proc")


def check_args(args):
    """
    Check the command line arguments, aborting if there is a problem
    """
    for name, year in (("BEGYEAR", args.begyear), ("ENDYEAR", args.endyear), ("YEAR", args.year)):
        if not 0 <= year <= 9999:
            abort(f"{name} must be an integer from 0 to 9999: {year}")
    if args.begyear > args.endyear:
        abort(f"BEGYEAR ({args.begyear}) must not be after ENDYEAR ({args.endyear})")
    if args.njobs is not None and args.njobs < 1:
        abort(f"NJOBS must be a positive integer: {args.njobs}")
    if args.account is not None and not args.batch:
        abort("--account is only used with --batch")
    if not args.all_vars and not build_var_list(args.vars, []):
        abort("Variable list given with --vars is empty")
    if not os.path.isdir(args.archdir):
        abort(f"Input {args.archdir} directory does not exist")
    if args.outdir is not None and not os.path.isdir(args.outdir):
        abort(f"Output {args.outdir} directory does not exist")
    indir = input_dir(args.archdir, args.case)
    if not os.path.isdir(indir):
        abort(f"{indir} directory does not exist")
    #
    # Abort if the output file already exists
    #
    todaysdate = datetime.date.today()
    outdir = output_dir(args.archdir, args.outdir, args.case)
    outfile = output_filename(outdir, args.case, args.begyear, args.endyear, todaysdate)
    if path.exists(outfile):
        abort(f"Output file {outfile} already exists, remove it or use a different output directory")


def build_var_list(variables, always_vars):
    """
    Combine the variables to average with the ones that are always output

    Keeps the order given, drops duplicates and empty names. Returns None when variables is None,
    meaning all variables.
    """
    if variables is None:
        return None
    var_list = []
    for var in list(variables) + list(always_vars):
        if var and var not in var_list:
            var_list.append(var)
    return var_list


def year_months(begyr, endyr):
    """
    List of (year, month) to process, from BEGYEAR-01 to ENDYEAR-12
    """
    return [(year, month) for year in range(begyr, endyr + 1) for month in range(1, 13)]


def index_daily_files(indir, case):
    """
    Find the daily coupler history files for this case, matching them by name (not by the time
    coordinate in the files)

    Returns a dictionary of sorted file path lists, with (year, month) keys
    """
    pattern = re.compile(
        rf"^{re.escape(case)}\.{re.escape(HIST_TAG)}\.(\d{{4}})-(\d{{2}})-(\d{{2}})-00000\.nc$"
    )
    files_by_month = {}
    for filename in sorted(os.listdir(indir)):
        match = pattern.match(filename)
        if match:
            key = (int(match.group(1)), int(match.group(2)))
            files_by_month.setdefault(key, []).append(os.path.join(indir, filename))
    return files_by_month


def check_day_counts(files_by_month, yms):
    """
    Check that every month has a full set of daily files (noleap calendar)

    Returns a list of error messages, empty if all the months are complete
    """
    errors = []
    for year, month in yms:
        nfiles = len(files_by_month.get((year, month), []))
        expect = NOLEAP_DAYS[month - 1]
        if nfiles != expect:
            errors.append(f"{year:04d}-{month:02d} has {nfiles} daily files (expected {expect})")
    return errors


def missing_vars(filename, var_list):
    """
    Return the variables in var_list that are not in the netCDF file
    """
    with Dataset(filename) as ncfile:
        file_vars = set(ncfile.variables)
    return [var for var in var_list if var not in file_vars]


def check_nco():
    """
    Abort if any of the needed NCO commands can't be found
    """
    missing = [cmd for cmd in NCO_COMMANDS if shutil.which(cmd) is None]
    if missing:
        abort(f"NCO commands not found: {' '.join(missing)} (on Derecho: module load nco)")


def check_input(args, var_list):
    """
    Check the input files before doing any work: every requested variable is in the first file,
    and every month in the range has a full set of daily files

    Returns the dictionary of daily files from index_daily_files
    """
    indir = input_dir(args.archdir, args.case)
    yms = year_months(args.begyear, args.endyear)
    log(
        logger,
        f"Processing {len(yms)} year-months: {args.begyear:04d}-01 .. {args.endyear:04d}-12",
    )
    files_by_month = index_daily_files(indir, args.case)

    if var_list is not None:
        first_files = files_by_month.get(yms[0], [])
        if not first_files:
            abort(f"No files for {yms[0][0]:04d}-{yms[0][1]:02d} in {indir}")
        log(logger, f"Checking the requested variables are in {first_files[0]}")
        missing = missing_vars(first_files[0], var_list)
        if missing:
            abort(f"Variables not found in {first_files[0]}: {' '.join(missing)}")

    log(logger, "Checking that every month has the expected number of daily files")
    errors = check_day_counts(files_by_month, yms)
    if errors:
        abort("\n".join(errors))
    return files_by_month


def ncra_command(infiles, outfile, var_list=None):
    """
    The ncra command (as a list) to average the input files into outfile
    If file already exists, returns an empty list and logs a message instead of running ncra
    """
    if path.exists(outfile):
        log(logger, f"File already exists:{outfile} so skipping ncra command")
        return []
    # Use --hst option so that history attribute isn't appended to as it's too long
    cmd = ["ncra", "--hst"]
    if var_list is not None:
        cmd += ["-v", ",".join(var_list)]
    return cmd + list(infiles) + ["-o", outfile]


def run_jobs(commands, njobs, *, progress_every=None, runner=run_cmd_output_on_error):
    """
    Run the commands (each a list, with the output file last), with at most njobs at once

    Stops at the first failure: the jobs that have not started are cancelled and the runner's
    error (normally an abort) is raised. Logs progress every progress_every completed jobs.
    """
    ndone = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=njobs) as executor:
        futures = [executor.submit(runner, cmd, f"Failed creating {cmd[-1]}") for cmd in commands]
        try:
            for future in concurrent.futures.as_completed(futures):
                future.result()
                ndone += 1
                if progress_every and ndone % progress_every == 0:
                    log(logger, f"  completed {ndone} of {len(commands)}")
        except BaseException:
            # Includes KeyboardInterrupt and the SystemExit from abort
            # Will die if the user hits Ctrl-C or if the batch job is killed
            for future in futures:
                future.cancel()
            raise


def mon_file(mondir, year, month):
    """
    Intermediate file with the monthly mean for one year and month
    """
    return os.path.join(mondir, f"mon_{year:04d}-{month:02d}.nc")


def clim_file(mondir, month):
    """
    Intermediate file with the climatology for one calendar month
    """
    return os.path.join(mondir, f"clim_{month:02d}.nc")


def monthly_means(files_by_month, yms, mondir, var_list, njobs):
    """
    Step 1: average the daily files of each month of each year
    """
    log(logger, f"Computing monthly means for each year ({njobs} ncra jobs at a time)")
    commands = [
        ncra_command(files_by_month[(year, month)], mon_file(mondir, year, month), var_list)
        for year, month in yms
    ]
    # Show progress every 10 years
    run_jobs(commands, njobs, progress_every=120)


def monthly_climatology(mondir, begyr, endyr, njobs):
    """
    Step 2: average each calendar month over the years

    Only this run's years are used, so leftover files from an earlier --keep run are not included
    """
    log(logger, f"Averaging each calendar month over the years {begyr}-{endyr}")
    commands = [
        ncra_command(
            [mon_file(mondir, year, month) for year in range(begyr, endyr + 1)],
            clim_file(mondir, month),
        )
        for month in range(1, 13)
    ]
    run_jobs(commands, njobs)


def output_filename(outdir, case, begyr, endyr, todaysdate):
    """
    Name of the output climatology file, with the year range and creation date
    """
    return os.path.join(
        outdir, f"{case}.cpl.hx.atm.monclim.{begyr:04d}-{endyr:04d}_c{todaysdate:%Y%m%d}.nc"
    )


def month_bounds():
    """
    Start and end of each month, in days since the start of a noleap year
    """
    bounds = []
    start = 0
    for ndays in NOLEAP_DAYS:
        bounds.append((start, start + ndays))
        start += ndays
    return bounds


def mid_month_times():
    """
    Middle of each month, in days since the start of a noleap year
    """
    return [(start + end) / 2.0 for start, end in month_bounds()]


def set_time_axis(filename, year):
    """
    Set the time and time_bnds of the 12-month file to the middle (and bounds) of
    each month of year, in the noleap calendar
    """
    units = f"days since {year:04d}-01-01 00:00:00"
    with Dataset(filename, "a") as ncfile:
        time = ncfile.variables["time"]
        time[:] = mid_month_times()
        time.setncattr("climatology_bounds", "climatology_bounds")
        time.setncattr("units", units)
        time.setncattr("calendar", "noleap")
        #
        # Change time_bounds to climatology_bounds or create it if it doesn't exist
        #
        if "time_bounds" in ncfile.variables:
            ncfile.renameVariable("time_bounds", "climatology_bounds")
        if "climatology_bounds" not in ncfile.variables:
            if "ntb" not in ncfile.dimensions:
                ncfile.createDimension("ntb", 2)
            ncfile.createVariable("climatology_bounds", "f8", ("time", "ntb"))
        climatology_bounds = ncfile.variables["climatology_bounds"]
        climatology_bounds[:] = month_bounds()
        climatology_bounds.setncattr("units", units)
        climatology_bounds.setncattr("calendar", "noleap")

def add_var_attributes(filename):
    """
    Add variable attributes to the file
    """
    with Dataset(filename, "a") as ncfile:
        for var in ncfile.variables:
            if var != "time" and var != "climatology_bounds":
                ncfile.variables[var].setncattr("cell_methods", "time: mean")
                ncfile.variables[var].setncattr("coordinates", "time atmImp_lon atmImp_lat")

def write_provenance(case, filename, begyear, endyear, year, todaysdate):
    """
    Add global attributes saying when, by whom and with what the file was created (like
    update_metadata in ctsm/site_and_regional/base_case.py, but using ncatted)
    """
    created_with = f"./{os.path.basename(WRAPPER)} -- {get_ctsm_git_short_hash()}"
    # Use Overwrite option since the file already exists and we want to add attributes to it
    cmd = ["ncatted", "-O"]
    comment = "Monthly climatology created from daily averaged CPLHIST files for DATM"
    title = "Monthly Climatology of Daily Averaged Atmosphere Coupler History Files for cyclical year {year}"
    source = f"Created from daily averaged CPLHIST files for case {case} for years {begyear:04d}-{endyear:04d}"
    institution = "National Science Foundation (NSF) - National Center for Atmospheric Research (NCAR) Community Earth System Model (CESM) project"
    for name, value in (
        ("Created_on", f"{todaysdate}"),
        ("Created_by", getpass.getuser()),
        ("Created_with", created_with),
        ("comment", comment),
        ("title", title),
        ("source", source),
        ("institution", institution),
        ("Conventions", "CF-1.13"),
        ("case", case),
    ):
        cmd += ["-a", f"{name},global,o,c,{value}"]
    run_cmd_output_on_error(cmd + [filename], f"Failed adding provenance to {filename}")


def write_output(case, mondir, outfile, begyear, endyear, year, todaysdate):
    """
    Step 3: put the 12 monthly climatologies in one file, reset its time axis to the middle of
    each month of year, and add provenance attributes
    """
    log(logger, "Putting the climatology into one file and finalizing its metadata")
    clim_files = [clim_file(mondir, month) for month in range(1, 13)]
    run_cmd_output_on_error(["ncrcat"] + clim_files + [outfile], f"Failed creating {outfile}")
    set_time_axis(outfile, year)
    add_var_attributes(outfile)
    write_provenance(case, outfile, begyear, endyear, year, todaysdate)
    log(logger, f"Wrote {outfile} (time = mid-month of year {year:04d}, noleap)")


def job_args(args):
    """
    The command line arguments for the batch job: the same as given, without the batch options
    """
    argv = [args.archdir, args.case, str(args.begyear), str(args.endyear)]
    if args.all_vars:
        argv.append("--all-vars")
    elif args.vars != DEFAULT_VARS:
        argv += ["--vars", ",".join(args.vars)]
    if args.outdir is not None:
        argv += ["--outdir", args.outdir]
    if args.year != DEFAULT_YEAR:
        argv += ["--year", str(args.year)]
    if args.keep:
        argv.append("--keep")
    if args.njobs is not None:
        argv += ["--njobs", str(args.njobs)]
    for option in ("verbose", "silent", "debug"):
        if getattr(args, option):
            argv.append(f"--{option}")
    return argv


def batch_script(argv, *, account, walltime, workdir):
    """
    Text of the PBS batch job script that runs this tool with the arguments argv
    Written out to a variable so it can be submitted to qsub from stdin, instead of writing a temporary file
    """
    runfile = io.StringIO()
    write_runscript_part1(
        number_of_nodes=1,
        tasks_per_node=argv.njobs,
        machine="derecho",
        account=account,
        walltime=walltime,
        runfile=runfile,
        name=PBS_JOBNAME,
        comment="This is a batch script to create a monthly climatology from daily CPLHIST files",
    )
    runfile.write("module load nco\n")
    runfile.write(f"cd {shlex.quote(workdir)}\n")
    wrapper = os.path.join(path_to_ctsm_root(), WRAPPER)
    runfile.write(shlex.join([sys.executable, wrapper] + argv) + "\n")
    return runfile.getvalue()


def submit_batch(script, queue):
    """
    Submit the batch job script to PBS (qsub reads it from stdin), and return the job id
    """
    if shutil.which("qsub") is None:
        abort("qsub not found, can not submit with --batch")
    result = subprocess.run(
        ["qsub", "-q", queue], input=script, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        abort(f"qsub failed:\n{result.stdout}{result.stderr}")
    return result.stdout.strip()


def main():
    """
    Main program for daily_to_monthly_clim
    """
    setup_logging_pre_config()
    args = get_parser().parse_args()
    process_logging_args(args)
    check_args(args)

    njobs = args.njobs if args.njobs is not None else default_njobs()
    var_list = None if args.all_vars else build_var_list(args.vars, ALWAYS_VARS)
    if var_list is None:
        log(logger, "Averaging all variables")
    else:
        log(logger, f"Averaging variables: {','.join(var_list)}")
    log(logger, "Create a monthly climatology from daily averaged CPLHIST files for DATM")

    check_nco()
    files_by_month = check_input(args, var_list)

    # With --batch: the input has been checked, so submit to PBS and stop here
    if args.batch:
        account = args.account or os.environ.get("ACCOUNT") or DEFAULT_ACCOUNT
        argv = job_args(args)
        script = batch_script(argv, account=account, walltime=args.walltime, workdir=os.getcwd())
        logger.debug("Batch job script:\n%s", script)
        jobid = submit_batch(script, args.queue)
        log(
            logger,
            f"Successfully submitted job {jobid} with arguments: {' '.join(argv)}"
            f" (log file will be {PBS_JOBNAME}.o{jobid.split('.')[0]} in {os.getcwd()})",
        )
        return

    outdir = output_dir(args.archdir, args.outdir, args.case)
    mondir = os.path.join(outdir, "mon")
    os.makedirs(mondir, exist_ok=True)
    todaysdate = datetime.date.today()

    monthly_means(files_by_month, year_months(args.begyear, args.endyear), mondir, var_list, njobs)
    monthly_climatology(mondir, args.begyear, args.endyear, njobs)
    outfile = output_filename(outdir, args.case, args.begyear, args.endyear, todaysdate)
    write_output(args.case, mondir, outfile, args.begyear, args.endyear, args.year, todaysdate)

    if args.keep:
        log(logger, f"Kept intermediate files in {mondir}")
    else:
        shutil.rmtree(mondir)

    log(logger, f"Successfully created monthly climatology: {outfile}")
