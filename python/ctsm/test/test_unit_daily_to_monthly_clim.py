#!/usr/bin/env python3

"""Unit tests for daily_to_monthly_clim"""

import datetime
import io
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm import unit_testing
from ctsm.cplhist import daily_to_monthly_clim as d2mc

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

CASE = "TESTCASE"


# A couple utility functions for the tests
def _daily_name(case, year, month, day):
    """Name of a daily coupler history file"""
    return f"{case}.{d2mc.HIST_TAG}.{year:04d}-{month:02d}-{day:02d}-00000.nc"


def _write_nc(filename, variables, ntime=1):
    """
    Write a small netCDF file with the dimensions time (unlimited, ntime records) and n (2),
    and the given variables: a dict of name -> (dimensions, attributes)
    """
    with Dataset(filename, "w") as ncfile:
        ncfile.createDimension("time", None)
        ncfile.createDimension("n", 2)
        for name, (dims, attrs) in variables.items():
            var = ncfile.createVariable(name, "f8", dims)
            var.setncatts(attrs)
            shape = tuple(ntime if dim == "time" else 2 for dim in dims)
            var[...] = np.ones(shape)


#
# Tests for argument parsing
#
class TestGetParser(unittest.TestCase):
    """Tests of get_parser"""

    def test_defaults(self):
        """The defaults when only the positional arguments are given"""
        args = d2mc.get_parser().parse_args(["arch", CASE, "2001", "2002"])
        self.assertEqual(args.archdir, "arch")
        self.assertEqual(args.case, CASE)
        self.assertEqual((args.begyear, args.endyear), (2001, 2002))
        self.assertEqual(args.vars, d2mc.DEFAULT_VARS)
        self.assertFalse(args.all_vars)
        self.assertIsNone(args.outdir)
        self.assertEqual(args.year, d2mc.DEFAULT_YEAR)
        self.assertFalse(args.keep)
        self.assertIsNone(args.njobs)
        self.assertFalse(args.batch)
        self.assertIsNone(args.account)
        self.assertEqual(args.queue, d2mc.DEFAULT_QUEUE)
        self.assertEqual(args.walltime, d2mc.DEFAULT_WALLTIME)

    def test_varsSplitsCommaList(self):
        """--vars is split on commas, with spaces stripped"""
        args = d2mc.get_parser().parse_args(["--vars", "a, b,c", "arch", CASE, "1", "2"])
        self.assertEqual(args.vars, ["a", "b", "c"])

    def test_allVarsWithVars_exits(self):
        """--all-vars and --vars can't be used together"""
        with mock.patch("sys.stderr", new_callable=io.StringIO):
            with self.assertRaises(SystemExit):
                d2mc.get_parser().parse_args(["-a", "--vars", "a", "arch", CASE, "1", "2"])


#
# Tests for the argument error checking
#
class TestCheckArgs(unittest.TestCase):
    """Tests of check_args"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()
        self._archdir = os.path.join(self._testdir, "archive")
        os.makedirs(d2mc.input_dir(self._archdir, CASE))

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def _args(self, *options, archdir=None, years=("2001", "2002")):
        """Parse the options, with the archive directory, case and years"""
        if archdir is None:
            archdir = self._archdir
        return d2mc.get_parser().parse_args(list(options) + [archdir, CASE] + list(years))

    def test_validArgs(self):
        """Valid arguments pass"""
        d2mc.check_args(self._args())

    def test_yearOutOfRange(self):
        """BEGYEAR, ENDYEAR and YEAR must be from 1 to 9999"""
        with self.assertRaisesRegex(SystemExit, "ENDYEAR must be an integer from 1 to 9999"):
            d2mc.check_args(self._args(years=("2001", "10000")))
        with self.assertRaisesRegex(SystemExit, "BEGYEAR must be an integer from 1 to 9999"):
            d2mc.check_args(self._args(years=("0", "2001")))
        with self.assertRaisesRegex(SystemExit, "YEAR must be an integer from 1 to 9999"):
            d2mc.check_args(self._args("--year", "10000"))
        with self.assertRaisesRegex(SystemExit, "YEAR must be an integer from 1 to 9999"):
            d2mc.check_args(self._args("--year", "0"))

    def test_begyearAfterEndyear(self):
        """BEGYEAR can't be after ENDYEAR"""
        with self.assertRaisesRegex(SystemExit, "must not be after ENDYEAR"):
            d2mc.check_args(self._args(years=("2002", "2001")))

    def test_njobsZero(self):
        """--njobs must be positive"""
        with self.assertRaisesRegex(SystemExit, "NJOBS must be a positive integer more than one"):
            d2mc.check_args(self._args("--njobs", "0"))

    def test_accountWithoutBatch(self):
        """--account is only allowed with --batch"""
        with self.assertRaisesRegex(SystemExit, "--account is only used with --batch"):
            d2mc.check_args(self._args("--account", "P1"))
        d2mc.check_args(self._args("--batch", "--account", "P1"))

    def test_onlyAlwaysVars(self):
        """--vars with only the ALWAYS_VARS leaves nothing to average"""
        with self.assertRaisesRegex(SystemExit, "No variables to average"):
            d2mc.check_args(self._args("--vars", ",".join(d2mc.ALWAYS_VARS)))

    def test_missingDirectories(self):
        """The archive, output and input directories must exist"""
        with self.assertRaisesRegex(SystemExit, "Input .* directory does not exist"):
            d2mc.check_args(
                self._args(archdir=os.path.join(self._testdir, "dir_that_does_not_exist"))
            )
        with self.assertRaisesRegex(SystemExit, "Output .* directory does not exist"):
            d2mc.check_args(
                self._args("--outdir", os.path.join(self._testdir, "dir_that_does_not_exist"))
            )
        shutil.rmtree(d2mc.input_dir(self._archdir, CASE))
        with self.assertRaisesRegex(SystemExit, "case archive directory does not exist"):
            d2mc.check_args(self._args())

    def test_outputExists(self):
        """Abort if today's output file already exists"""
        outdir = d2mc.output_dir(self._archdir, None, CASE)
        os.makedirs(outdir)
        outfile = d2mc.output_filename(outdir, CASE, 2001, 2002, datetime.date.today())
        with open(outfile, "w", encoding="utf-8"):
            pass
        with self.assertRaisesRegex(SystemExit, "already exists"):
            d2mc.check_args(self._args())


#
# Tests for batch and number of jobs
#
class TestDefaultNjobs(unittest.TestCase):
    """Tests of default_njobs"""

    def test_inBatchJob(self):
        """Inside a PBS job (PBS_JOBID set) use PBS_NJOBS"""
        with mock.patch.dict(os.environ, {"PBS_JOBID": "1234.desched1"}):
            self.assertEqual(d2mc.default_njobs(), d2mc.PBS_NJOBS)

    def test_notInBatchJob(self):
        """Outside a PBS job use INTERACTIVE_NJOBS"""
        with mock.patch.dict(os.environ):
            os.environ.pop("PBS_JOBID", None)
            self.assertEqual(d2mc.default_njobs(), d2mc.INTERACTIVE_NJOBS)


#
# Tests for the variable lists
#
class TestVarLists(unittest.TestCase):
    """Tests of build_var_list and averaged_vars"""

    def test_buildVarList(self):
        """Keeps the order, drops duplicates and empty names, adds the extra variables at the end"""
        self.assertEqual(d2mc.build_var_list(["b", "a", "", "b"], ["c", "a"]), ["b", "a", "c"])

    def test_buildVarList_none(self):
        """None (all variables) stays None"""
        self.assertIsNone(d2mc.build_var_list(None, ["c"]))

    def test_averagedVars(self):
        """The ALWAYS_VARS are removed from the list to average"""
        variables = ["a", d2mc.ALWAYS_VARS[0], "", "a", "b"]
        self.assertEqual(d2mc.averaged_vars(variables), ["a", "b"])

    def test_averagedVars_none(self):
        """None (all variables) stays None"""
        self.assertIsNone(d2mc.averaged_vars(None))


#
# Tests for the year/month list and the daily file indexing
#
class TestYearMonths(unittest.TestCase):
    """Tests of year_months"""

    def test_twoYears(self):
        """Two years in order"""
        yms = d2mc.year_months(2001, 2002)
        self.assertEqual(len(yms), 24)
        self.assertEqual(yms[0], (2001, 1))
        self.assertEqual(yms[11], (2001, 12))
        self.assertEqual(yms[12], (2002, 1))
        self.assertEqual(yms[-1], (2002, 12))

    def test_singleYear(self):
        """A single year, including year 1 and the default year"""
        year = 1
        self.assertEqual(d2mc.year_months(year, year), [(1, month) for month in range(1, 13)])
        year = d2mc.DEFAULT_YEAR
        self.assertEqual(d2mc.year_months(year, year), [(year, month) for month in range(1, 13)])


#
# Test the list of daily files and the file count
#
class TestDailyFiles(unittest.TestCase):
    """Tests of index_daily_files and check_day_counts"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()
        # Empty files are enough, the files are found by name
        for month, ndays in enumerate(d2mc.NOLEAP_DAYS, start=1):
            for day in range(1, ndays + 1):
                self._create_empty_file(_daily_name(CASE, 2001, month, day))
        # Files that should be ignored: another case, and other kinds of files
        self._create_empty_file(_daily_name("DIFFERENT_CASE", 2001, 1, 1))
        self._create_empty_file(f"{CASE}.cpl.hi.2001-01-01-00000.nc")
        self._create_empty_file(_daily_name(CASE, 2001, 1, 1) + ".base")

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def _create_empty_file(self, filename):
        """Create an empty file in the test directory"""
        with open(os.path.join(self._testdir, filename), "w", encoding="utf-8"):
            pass

    def test_indexDailyFiles(self):
        """
        Check that index_daily_files properly groups and sorts by year and month
        and ignores any other files
        """
        # Files that should be ignored: another case, and other kinds of files
        self._create_empty_file(_daily_name("DIFFERENT_CASE", 2001, 1, 1))
        self._create_empty_file(f"{CASE}.cpl.hi.2001-01-01-00000.nc")
        self._create_empty_file(_daily_name(CASE, 2001, 1, 1) + ".base")
        #
        files_by_month = d2mc.index_daily_files(self._testdir, CASE)
        self.assertEqual(sorted(files_by_month), d2mc.year_months(2001, 2001))
        for (_, month), files in files_by_month.items():
            self.assertEqual(len(files), d2mc.NOLEAP_DAYS[month - 1])
            self.assertEqual(files, sorted(files))
        self.assertEqual(
            files_by_month[(2001, 1)][0], os.path.join(self._testdir, _daily_name(CASE, 2001, 1, 1))
        )

    def test_checkDayCounts_complete(self):
        """Test that checkDayCounts for a complete year is correct"""
        files_by_month = d2mc.index_daily_files(self._testdir, CASE)
        self.assertEqual(d2mc.check_day_counts(files_by_month, d2mc.year_months(2001, 2001)), [])

    def test_checkDayCounts_missing(self):
        """Test that checkDayCounts reports missing days and years correctly"""
        # Remove one day from the first year
        os.remove(os.path.join(self._testdir, _daily_name(CASE, 2001, 2, 14)))
        # NOTE: Setup only created files for the first year, so the second year has no files
        files_by_month = d2mc.index_daily_files(self._testdir, CASE)
        errors = d2mc.check_day_counts(files_by_month, d2mc.year_months(2001, 2002))
        self.assertEqual(len(errors), 13)
        # There will be an error reported for the missing day in the first year
        # and each month of the second year
        self.assertEqual(errors[0], "2001-02 has 27 daily files (expected 28)")
        self.assertEqual(errors[1], "2002-01 has 0 daily files (expected 31)")


#
# Test the netCDF variable handling functions for the monthly files that were kept
# previously -- to make sure they have the needed variables
#
class TestNetcdfVars(unittest.TestCase):
    """Tests of file_vars, missing_vars and check_reused_mon_file"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()
        self._daily = os.path.join(self._testdir, "daily.nc")
        _write_nc(
            self._daily,
            {
                "a": (("time", "n"), {}),
                "b": (("time", "n"), {}),
                d2mc.ALWAYS_VARS[0]: (("time", "n"), {}),
            },
        )

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def _monfile(self, names):
        """Write a monthly file with the given variables (on time and n)"""
        monfile = os.path.join(self._testdir, "mon.nc")
        _write_nc(monfile, {name: (("time", "n"), {}) for name in names})
        return monfile

    def test_fileVars(self):
        """The set of variable names in the file"""
        self.assertEqual(d2mc.file_vars(self._daily), {"a", "b", d2mc.ALWAYS_VARS[0]})

    def test_missingVars(self):
        """Only the variables not in the file are returned, in order"""
        self.assertEqual(d2mc.missing_vars(self._daily, ["c", "a", "d"]), ["c", "d"])

    def test_checkReusedMonFile_ok(self):
        """A monthly file with the needed variables passes, with a list or all variables"""
        monfile = self._monfile(["a", "b"])
        d2mc.check_reused_mon_file(monfile, [self._daily], ["a"])
        d2mc.check_reused_mon_file(monfile, [self._daily], None)

    def test_checkReusedMonFile_missing(self):
        """A monthly file missing a requested variable aborts"""
        monfile = self._monfile(["a"])
        with self.assertRaisesRegex(SystemExit, "is missing variables: c"):
            d2mc.check_reused_mon_file(monfile, [self._daily], ["a", "c"])
        with self.assertRaisesRegex(SystemExit, "is missing variables: b"):
            d2mc.check_reused_mon_file(monfile, [self._daily], None)

    # This checks against an original version that had the ALWAYS_VARS on the monthly files
    # TODO: Remove this since it's not really needed
    def test_checkReusedMonFile_stale(self):
        """A monthly file with the ALWAYS_VARS was made by an older version, so aborts"""
        monfile = self._monfile(["a", "b", d2mc.ALWAYS_VARS[0]])
        with self.assertRaisesRegex(SystemExit, "older version"):
            d2mc.check_reused_mon_file(monfile, [self._daily], ["a"])


#
# Check some of the NCO (NetCDF Climate Operators) functions
#
class TestCheckNco(unittest.TestCase):
    """Tests of check_nco"""

    def test_missingCommand(self):
        """Abort naming the NCO commands that aren't found"""
        # mock d2mc.shutil with the input arguments below, and the side effects listed
        with mock.patch.object(
            # lambda is used here to have a small inline function to define the side effects
            d2mc.shutil,
            "which",
            side_effect=lambda cmd: None if cmd == "ncrcat" else "/bin/x",
        ):
            with self.assertRaisesRegex(SystemExit, "NCO commands not found: ncrcat "):
                d2mc.check_nco()

    def test_allFound(self):
        """No abort when all the commands are found"""
        # mock d2mc.shutil with the input arguments below, and the side effects as listed
        with mock.patch.object(d2mc.shutil, "which", return_value="/bin/x"):
            d2mc.check_nco()


#
# Tests of the ncra command for reacord averaging
#
class TestNcraCommand(unittest.TestCase):
    """Tests of ncra_command"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()
        self._out = os.path.join(self._testdir, "out.nc")

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_varList(self):
        """A variable list gives -v, and exclude is ignored"""
        cmd = d2mc.ncra_command(["f1", "f2"], self._out, ["a", "b"], exclude=["c"])
        self.assertEqual(cmd, ["ncra", "-O", "--hst", "-v", "a,b", "f1", "f2", "-o", self._out])

    def test_allVarsExclude(self):
        """All variables with exclude gives -x -v"""
        cmd = d2mc.ncra_command(["f1"], self._out, None, exclude=["c", "d"])
        self.assertEqual(cmd, ["ncra", "-O", "--hst", "-x", "-v", "c,d", "f1", "-o", self._out])

    def test_allVars(self):
        """All variables without exclude averages everything"""
        cmd = d2mc.ncra_command(["f1"], self._out)
        self.assertEqual(cmd, ["ncra", "-O", "--hst", "f1", "-o", self._out])

    def test_skipExisting(self):
        """With skip_existing, an existing output file gives no command"""
        expected = ["ncra", "-O", "--hst", "f1", "-o", self._out]
        self.assertEqual(d2mc.ncra_command(["f1"], self._out, skip_existing=True), expected)
        with open(self._out, "w", encoding="utf-8"):
            pass
        self.assertEqual(d2mc.ncra_command(["f1"], self._out, skip_existing=True), [])
        self.assertEqual(d2mc.ncra_command(["f1"], self._out), expected)


#
# Tests of running jobs
#
class TestRunJobs(unittest.TestCase):
    """Tests of run_jobs (with a mock runner instead of running the commands)"""

    def test_allRun(self):
        """Every non-empty command is run, with an error message naming its output file"""
        # Create runner as a MagicMock so that it does nothing, but records each call
        runner = mock.MagicMock()
        # Run a job in one call sending it three different commands to run at the same time
        commands = [
            ["a", "o1"],  # run "a" creating file o1
            [],  # empty command (output already exists) is skipped
            ["b", "o2"],  # run "b" creating file o2
        ]
        d2mc.run_jobs(commands, 2, runner=runner)
        # Check that runner had recorded these two calls, and the appropriate error
        runner.assert_has_calls(
            [
                mock.call(["a", "o1"], "Failed creating o1"),
                mock.call(["b", "o2"], "Failed creating o2"),
            ],
            any_order=True,  # running jobs simultaneously we don't know which one will be first
        )
        self.assertEqual(runner.call_count, 2)

    def test_oneJobInOrder(self):
        """With one job at a time the commands run in order"""
        # Create runner as a MagicMock so that it does nothing, but records each call
        runner = mock.MagicMock()
        commands = [[f"c{i}", f"o{i}"] for i in range(5)]
        d2mc.run_jobs(commands, 1, runner=runner)
        # Check that the calls ran from runner match what's expected
        self.assertEqual([call.args[0] for call in runner.call_args_list], commands)

    def test_noCommands(self):
        """Only empty commands runs nothing"""
        # Create runner as a MagicMock so that it does nothing, but records each call
        runner = mock.MagicMock()
        d2mc.run_jobs([[], []], 4, runner=runner)
        runner.assert_not_called()

    def test_failure(self):
        """A failing command (the runner aborts) stops with SystemExit"""

        # Define a runner that aborts if the input command is "bad"
        def runner(cmd, errmsg):
            if cmd[0] == "bad":
                sys.exit(errmsg)

        with self.assertRaisesRegex(SystemExit, "Failed creating o2"):
            d2mc.run_jobs([["good", "o1"], ["bad", "o2"], ["good", "o3"]], 2, runner=runner)


#
# Tests of the different output files
#
class TestFilenames(unittest.TestCase):
    """Tests of output_filename, mon_file and clim_file"""

    def test_outputFilename(self):
        """Zero-padded years and the creation date"""
        self.assertEqual(
            d2mc.output_filename("/out", CASE, 100, 130, datetime.date(2026, 10, 3)),
            os.path.join("/out", f"{CASE}.cpl.hx.atm.monclim.0100-0130_c20261003.nc"),
        )

    def test_intermediateFiles(self):
        """Names of the monthly and climatology files"""
        self.assertEqual(d2mc.mon_file("/m", 100, 2), os.path.join("/m", "mon_0100-02.nc"))
        self.assertEqual(d2mc.clim_file("/m", 2), os.path.join("/m", "clim_02.nc"))


#
# Test file output for time and time and climaatology bounds
#
class TestTimeValues(unittest.TestCase):
    """Tests of month_bounds, mid_month_times and climatology_bounds_values"""

    def test_monthBounds(self):
        """Start and end of each month of a noleap year"""
        bounds = d2mc.month_bounds()
        self.assertEqual(len(bounds), 12)  # Check that there is just 12 bounds for each month
        self.assertEqual(bounds[0], (0, 31))  # January
        self.assertEqual(bounds[1], (31, 59))  # February
        self.assertEqual(bounds[-1], (334, 365))  # December (last month)

    def test_midMonthTimes(self):
        """Middle of each month"""
        times = d2mc.mid_month_times()
        self.assertEqual(times[0], 15.5)  # January
        self.assertEqual(times[1], 45.0)  # Febreary
        self.assertEqual(times[-1], 349.5)  # December (last month)

    def test_climatologyBounds_singleYear(self):
        """Averaging only YEAR gives the bounds of each month of YEAR"""
        self.assertEqual(d2mc.climatology_bounds_values(1850, 1850, 1850), d2mc.month_bounds())

    def test_climatologyBounds_span(self):
        """The bounds go from the month in BEGYEAR to the month in ENDYEAR, negative before YEAR"""
        bounds = d2mc.climatology_bounds_values(100, 130, 1850)
        self.assertEqual(bounds[0], ((100 - 1850) * 365, (130 - 1850) * 365 + 31))
        # December ends at the start of ENDYEAR+1
        self.assertEqual(bounds[-1][1], (131 - 1850) * 365)
        self.assertTrue(all(start < 0 and end < 0 for start, end in bounds))


#
# Test setting of the time axis
#
class TestSetTimeAxis(unittest.TestCase):
    """Tests of set_time_axis"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()
        self._file = os.path.join(self._testdir, "clim.nc")

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def _write_clim(self, with_time_bnds=False):
        """
        Write a 12-record file like the one from ncrcat, optionally with time_bnds
        This is a utility function called both ways in later tests
        """
        with Dataset(self._file, "w") as ncfile:
            ncfile.createDimension("time", None)
            time = ncfile.createVariable("time", "f8", ("time",))
            time.setncatts({"units": "days since 0001-01-01", "bounds": "time_bnds"})
            time.setncattr("cell_methods", "time: mean")
            time[:] = np.arange(12.0)
            if with_time_bnds:
                ncfile.createDimension("ntb", 2)
                time_bnds = ncfile.createVariable("time_bnds", "f8", ("time", "ntb"))
                time_bnds.setncatts({"axis": "T", "calendar": "noleap", "bounds": "time_bnds"})
                time_bnds.setncattr("cell_methods", "time: mean")
                time_bnds[:] = np.zeros((12, 2))

    def _check(self):
        """Check the time axis and climatology bounds after set_time_axis"""
        d2mc.set_time_axis(self._file, 1850, 2001, 2002)
        units = "days since 1850-01-01 00:00:00"
        with Dataset(self._file) as ncfile:
            self.assertNotIn("time_bnds", ncfile.variables)
            time = ncfile.variables["time"]
            np.testing.assert_array_equal(time[:], d2mc.mid_month_times())
            self.assertEqual(time.units, units)
            self.assertEqual(time.calendar, "noleap")
            self.assertEqual(time.climatology, "climatology_bounds")
            self.assertNotIn("bounds", time.ncattrs())
            self.assertNotIn("cell_methods", time.ncattrs())
            bounds = ncfile.variables["climatology_bounds"]
            self.assertEqual(bounds.dimensions, ("time", "ntb"))
            np.testing.assert_array_equal(
                bounds[:], d2mc.climatology_bounds_values(2001, 2002, 1850)
            )
            self.assertEqual(bounds.units, units)
            for attr in ("axis", "calendar", "bounds", "cell_methods"):
                self.assertNotIn(attr, bounds.ncattrs())

    def test_withTimeBnds(self):
        """time_bnds is renamed to climatology_bounds"""
        self._write_clim(with_time_bnds=True)
        self._check()

    def test_withoutTimeBnds(self):
        """climatology_bounds is created when there is no time_bnds"""
        self._write_clim(with_time_bnds=False)
        self._check()


#
# Test the copy of the grid varaibles (the ALWAYS_VARS)
#
class TestCopyGridVars(unittest.TestCase):
    """Tests of copy_grid_vars"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()
        self._src = os.path.join(self._testdir, "daily.nc")
        self._dst = os.path.join(self._testdir, "clim.nc")
        # Create a CPLHIST history file like from CMEPS:
        # float lon/lat on (time, ny, nx) with a fill value and attributes
        with Dataset(self._src, "w") as ncfile:
            ncfile.createDimension("time", None)
            ncfile.createDimension("atmImp_ny", 1)
            ncfile.createDimension("atmImp_nx", 3)
            for name, base in (("atmImp_lon", 10.0), ("atmImp_lat", -60.0)):
                var = ncfile.createVariable(
                    name, "f4", ("time", "atmImp_ny", "atmImp_nx"), fill_value=1.0e30
                )
                var.setncatts({"long_name": name, "units": "degrees"})
                var.setncattr("cell_methods", "time: mean")
                # Different values for the second time, so the first time slice can be checked
                var[0] = base + np.arange(3.0).reshape(1, 3)
                var[1] = -999.0
        # The destination has atmImp_ny but not atmImp_nx
        with Dataset(self._dst, "w") as ncfile:
            ncfile.createDimension("time", None)
            ncfile.createDimension("atmImp_ny", 1)

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_copyGridVars(self):
        """Copied without time, with the first time slice, attributes and fill value"""
        d2mc.copy_grid_vars(self._src, self._dst, ["atmImp_lon", "atmImp_lat"])
        with Dataset(self._dst) as ncfile:
            self.assertEqual(len(ncfile.dimensions["atmImp_nx"]), 3)
            for name, base in (("atmImp_lon", 10.0), ("atmImp_lat", -60.0)):
                var = ncfile.variables[name]
                self.assertEqual(var.dimensions, ("atmImp_ny", "atmImp_nx"))
                self.assertEqual(var.dtype, np.float32)
                self.assertAlmostEqual(float(var.getncattr("_FillValue")), 1.0e30, delta=1.0e24)
                self.assertEqual(var.long_name, name)
                self.assertEqual(var.units, "degrees")
                self.assertNotIn("cell_methods", var.ncattrs())
                np.testing.assert_array_equal(var[:], base + np.arange(3.0).reshape(1, 3))


#
# Test that adding variable attributes works correctly
#
class TestAddVarAttributes(unittest.TestCase):
    """Tests of add_var_attributes"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()
        self._file = os.path.join(self._testdir, "clim.nc")
        lon, lat = d2mc.ALWAYS_VARS
        _write_nc(
            self._file,
            {
                "time": (("time",), {}),
                "climatology_bounds": (("time", "n"), {}),
                lon: (("n",), {"cell_methods": "time: mean"}),
                lat: (("n",), {}),
                "a": (("time", "n"), {"standard_name": "a", "cell_methods": "time: mean"}),
                "b": (("time", "n"), {"long_name": "Keep this name"}),
                "c": (("n",), {}),
            },
        )

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_addVarAttributes(self):
        """Only the data variables on time get the climatology attributes"""
        d2mc.add_var_attributes(self._file)
        cell_methods = "time: mean within years time: mean over years"
        with Dataset(self._file) as ncfile:
            for name in ("a", "b"):
                var = ncfile.variables[name]
                self.assertEqual(var.cell_methods, cell_methods)
                self.assertEqual(var.coordinates, " ".join(d2mc.ALWAYS_VARS))
                self.assertNotIn("standard_name", var.ncattrs())
            self.assertEqual(ncfile.variables["a"].long_name, "a")
            self.assertEqual(ncfile.variables["b"].long_name, "Keep this name")
            for name in ["time", "climatology_bounds", "c"] + d2mc.ALWAYS_VARS:
                self.assertNotIn("cell_methods", ncfile.variables[name].ncattrs())
                self.assertNotIn("coordinates", ncfile.variables[name].ncattrs())


#
# Tests of writing prvenance information is correct
#
class TestWriteProvenance(unittest.TestCase):
    """Tests of write_provenance (with the ncatted command mocked)"""

    def test_oneNcattedCommand(self):
        """One ncatted --hst command sets the global attributes and deletes file_version"""
        # mock run_cmnd_output_on_error function as run_cmd that does nothing, but returns as set
        # mock sys.argv to work as given next
        with mock.patch.object(d2mc, "run_cmd_output_on_error") as run_cmd, mock.patch.object(
            d2mc, "get_ctsm_git_short_hash", return_value="abc1234"
        ), mock.patch.object(sys, "argv", ["daily_to_monthly_clim", "--vars", "a,b", "x"]):
            d2mc.write_provenance(CASE, "out.nc", 100, 130, 1850, datetime.date(2026, 10, 3))
        run_cmd.assert_called_once()
        cmd = run_cmd.call_args.args[0]
        self.assertEqual(cmd[:3], ["ncatted", "-O", "--hst"])
        self.assertEqual(cmd[-1], "out.nc")
        attrs = cmd[4:-1:2]
        self.assertEqual(cmd[3:-1:2], ["-a"] * len(attrs))
        self.assertIn("Conventions,global,o,c,CF-1.11", attrs)
        self.assertIn("Created_on,global,o,c,2026-10-03", attrs)
        self.assertIn("Created_with,global,o,c,./daily_to_monthly_clim -- abc1234", attrs)
        self.assertIn("command_used,global,o,c,daily_to_monthly_clim --vars a,b x", attrs)
        self.assertIn(f"made_from_cesm_case,global,o,c,{CASE}", attrs)
        self.assertIn("file_version,global,d,,", attrs)
        title = [attr for attr in attrs if attr.startswith("title,")]
        self.assertTrue(title[0].endswith("for cyclical year 1850"))


#
# Test the batch options
#
class TestBatch(unittest.TestCase):
    """Tests of job_args, batch_script and submit_batch"""

    def test_jobArgs_default(self):
        """Only the positional arguments when the defaults are used"""
        args = d2mc.get_parser().parse_args(["arch", CASE, "2001", "2002"])
        self.assertEqual(d2mc.job_args(args), ["arch", CASE, "2001", "2002"])

    def test_jobArgs_removes_batch_options(self):
        """Set all options, but show that job_args returns only the non-batch options"""
        args = d2mc.get_parser().parse_args(
            [
                "--batch",
                "--account",
                "P1",
                "--queue",
                "develop",
                "--walltime",
                "01:00:00",
                "--vars",
                "a,b",
                "-o",
                "out",
                "-y",
                "2000",
                "-k",
                "-j",
                "8",
                "-v",
                "arch",
                CASE,
                "2001",
                "2002",
            ]
        )
        self.assertEqual(
            d2mc.job_args(args),
            ["arch", CASE, "2001", "2002", "--vars", "a,b", "--outdir", "out", "--year", "2000"]
            + ["--keep", "--njobs", "8", "--verbose"],
        )

    def test_jobArgs_allVars(self):
        """Test that --all-vars is passed on when it's called as '-a'"""
        args = d2mc.get_parser().parse_args(["-a", "arch", CASE, "2001", "2002"])
        self.assertEqual(d2mc.job_args(args), ["arch", CASE, "2001", "2002", "--all-vars"])

    def test_batchScript(self):
        """
        Test that batch_script creates: PBS settings, NCO, the directory and the quoted command
        """
        argv = ["arch", CASE, "2001", "2002", "--vars", "a,b", "--outdir", "my dir"]
        script = d2mc.batch_script(
            argv, account="P1", walltime="01:00:00", workdir="/work dir", njobs=8
        )
        lines = script.splitlines()
        self.assertEqual(lines[0], "#!/bin/bash")
        self.assertIn(f"#PBS -N {d2mc.PBS_JOBNAME}", lines)
        self.assertIn("#PBS -A P1", lines)
        self.assertIn("#PBS -l walltime=01:00:00", lines)
        self.assertIn("mpiprocs=8", script)
        self.assertIn("module load nco", lines)
        self.assertIn("cd '/work dir'", lines)
        command = shlex.split(lines[-1])
        self.assertEqual(command[0], sys.executable)
        self.assertTrue(command[1].endswith(d2mc.WRAPPER))
        self.assertEqual(command[2:], argv)

    def test_submitBatch(self):
        """
        Test submitBatch that qsub gets the queue and the script on stdin,
        and the job id is returned
        """
        result = subprocess.CompletedProcess([], 0, stdout="1234.desched1\n", stderr="")
        # mock d2mc.shutil for which to just return qsub path
        # Then mock d2mc.subprocess as run to return the result
        with mock.patch.object(d2mc.shutil, "which", return_value="/bin/qsub"), mock.patch.object(
            d2mc.subprocess, "run", return_value=result
        ) as run:
            self.assertEqual(d2mc.submit_batch("the script", "develop"), "1234.desched1")
        self.assertEqual(run.call_args.args[0], ["qsub", "-q", "develop"])
        self.assertEqual(run.call_args.kwargs["input"], "the script")

    def test_submitBatch_fails(self):
        """test that submitBatch fails if the qsub aborts"""
        # Create a result that mimics a qsub error
        result = subprocess.CompletedProcess([], 1, stdout="", stderr="qsub: bad account")
        # mock d2mc.shutil on which to return qsub path
        # Then mock d2mc.subprocess run to return the result
        with mock.patch.object(d2mc.shutil, "which", return_value="/bin/qsub"), mock.patch.object(
            d2mc.subprocess, "run", return_value=result
        ):
            with self.assertRaisesRegex(SystemExit, "qsub failed:\nqsub: bad account"):
                d2mc.submit_batch("the script", "main")

    def test_submitBatch_noQsub(self):
        """Abort if qsub isn't found"""
        # mock d2mc.shutil for which, to return nothing as it would if the command wasn't found
        with mock.patch.object(d2mc.shutil, "which", return_value=None):
            with self.assertRaisesRegex(SystemExit, "qsub not found"):
                d2mc.submit_batch("the script", "main")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
