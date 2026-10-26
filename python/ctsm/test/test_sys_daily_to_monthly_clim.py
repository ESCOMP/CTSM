#!/usr/bin/env python3

"""System tests for daily_to_monthly_clim (needs the NCO commands, on Derecho: module load nco)"""

import contextlib
import datetime
import io
import os
import shlex
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

# -- add python/ctsm  to path (needed if we want to run test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.cplhist import daily_to_monthly_clim as d2mc

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

CASE = "TC"
BEGYEAR = 2001
ENDYEAR = 2002
# A variable that isn't in the default list
EXTRA_VAR = "atmImp_Sa_extra"
NX = 3
LON = 10.0 + np.arange(NX, dtype="f4")
LAT = -60.0 + np.arange(NX, dtype="f4")


def _value(year, month):
    """Value of the default variables in the daily files for a year and month"""
    return (year - 2000) * 100 + month


def _write_daily_file(filename, year, month, day_of_year):
    """
    Write a daily file like the CMEPS 24h average coupler history files: 2 time samples, time_bnds,
    the default variables and EXTRA_VAR on (time, atmImp_ny, atmImp_nx), and lon/lat with time too
    """
    with Dataset(filename, "w") as ncfile:
        ncfile.createDimension("time", None)
        ncfile.createDimension("ntb", 2)
        ncfile.createDimension("atmImp_ny", 1)
        ncfile.createDimension("atmImp_nx", NX)
        ncfile.file_version = "cmeps0"
        start = (year - 1) * 365 + day_of_year
        time = ncfile.createVariable("time", "f8", ("time",))
        time.setncatts({"units": "days since 0001-01-01 00:00:00", "calendar": "noleap"})
        time[:] = [start + 0.25, start + 0.75]
        time_bnds = ncfile.createVariable("time_bnds", "f8", ("time", "ntb"))
        time_bnds.setncatts({"bounds": "time_bnds", "axis": "T", "calendar": "noleap"})
        time_bnds[:] = [[start, start + 0.5], [start + 0.5, start + 1.0]]
        dims = ("time", "atmImp_ny", "atmImp_nx")
        value = _value(year, month)
        for name in d2mc.DEFAULT_VARS + [EXTRA_VAR]:
            var = ncfile.createVariable(name, "f4", dims, fill_value=1.0e30)
            var.setncatts({"units": "kg m-2 s-1", "standard_name": name})
            var[:] = (2.0 * value if name == EXTRA_VAR else value) * np.ones((2, 1, NX))
        for name, values, units in (
            ("atmImp_lon", LON, "degrees_E"),
            ("atmImp_lat", LAT, "degrees_N"),
        ):
            var = ncfile.createVariable(name, "f4", dims, fill_value=1.0e30)
            var.setncatts({"long_name": name, "units": units})
            var[:] = values * np.ones((2, 1, NX))


@unittest.skipUnless(shutil.which("ncra"), "NCO not found (on Derecho: module load nco)")
class TestSysDailyToMonthlyClim(unittest.TestCase):
    """System tests for daily_to_monthly_clim, running the NCO commands"""

    @classmethod
    def setUpClass(cls):
        """Write BEGYEAR to ENDYEAR of daily files, once for all the tests"""
        cls._indir_root = tempfile.mkdtemp()
        cls._archdir = os.path.join(cls._indir_root, "archive")
        indir = d2mc.input_dir(cls._archdir, CASE)
        os.makedirs(indir)
        cls._first_daily = None
        for year in range(BEGYEAR, ENDYEAR + 1):
            day_of_year = 0
            for month, ndays in enumerate(d2mc.NOLEAP_DAYS, start=1):
                for day in range(1, ndays + 1):
                    filename = os.path.join(
                        indir, f"{CASE}.{d2mc.HIST_TAG}.{year:04d}-{month:02d}-{day:02d}-00000.nc"
                    )
                    _write_daily_file(filename, year, month, day_of_year)
                    day_of_year += 1
                    if cls._first_daily is None:
                        cls._first_daily = filename

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._indir_root, ignore_errors=True)

    def setUp(self):
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
        os.chdir(self._tempdir)

    def tearDown(self):
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def _run(self, *options, outdir=None):
        """
        Run main with the options (and 4 ncra jobs), writing under outdir (a new directory under
        the test directory by default). Returns the output file name and what was printed.
        """
        if outdir is None:
            outdir = tempfile.mkdtemp(dir=self._tempdir)
        argv = ["daily_to_monthly_clim", "-j", "4", "-o", outdir] + list(options)
        argv += [self._archdir, CASE, str(BEGYEAR), str(ENDYEAR)]
        output = io.StringIO()
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(output):
            d2mc.main()
        outfile = d2mc.output_filename(
            d2mc.output_dir(self._archdir, outdir, CASE),
            CASE,
            BEGYEAR,
            ENDYEAR,
            datetime.date.today(),
        )
        self.assertTrue(os.path.exists(outfile), msg=f"{outfile} not created")
        return outfile, output.getvalue(), argv

    def _check_output(self, outfile, data_vars, year=d2mc.DEFAULT_YEAR):
        """Check the variables, values, grid, time axis and metadata of the output file"""
        units = f"days since {year:04d}-01-01 00:00:00"
        with Dataset(outfile) as ncfile:
            expected = set(data_vars) | set(d2mc.ALWAYS_VARS) | {"time", "climatology_bounds"}
            self.assertEqual(set(ncfile.variables), expected)
            # Each month is the average over the years of that month
            for name in data_vars:
                var = ncfile.variables[name]
                self.assertEqual(var.dimensions, ("time", "atmImp_ny", "atmImp_nx"))
                for month in range(1, 13):
                    expect = np.mean([_value(yr, month) for yr in range(BEGYEAR, ENDYEAR + 1)])
                    if name == EXTRA_VAR:
                        expect *= 2.0
                    np.testing.assert_allclose(var[month - 1], expect, rtol=1.0e-6)
                self.assertEqual(var.cell_methods, "time: mean within years time: mean over years")
                self.assertEqual(var.coordinates, "atmImp_lon atmImp_lat")
                self.assertNotIn("standard_name", var.ncattrs())
            # The grid without the time dimension
            for name, values in (("atmImp_lon", LON), ("atmImp_lat", LAT)):
                var = ncfile.variables[name]
                self.assertEqual(var.dimensions, ("atmImp_ny", "atmImp_nx"))
                np.testing.assert_array_equal(var[0], values)
                self.assertNotIn("cell_methods", var.ncattrs())
            # The time axis and the climatology bounds
            time = ncfile.variables["time"]
            np.testing.assert_array_equal(time[:], d2mc.mid_month_times())
            self.assertEqual(time.units, units)
            self.assertEqual(time.calendar, "noleap")
            self.assertEqual(time.climatology, "climatology_bounds")
            np.testing.assert_array_equal(
                ncfile.variables["climatology_bounds"][:],
                d2mc.climatology_bounds_values(BEGYEAR, ENDYEAR, year),
            )
            self.assertEqual(ncfile.Conventions, "CF-1.11")
            self.assertEqual(ncfile.made_from_cesm_case, CASE)
            self.assertNotIn("file_version", ncfile.ncattrs())
            return ncfile.command_used

    def test_defaultVars(self):
        """The default variables, with the intermediate files removed at the end"""
        outfile, _, argv = self._run()
        command_used = self._check_output(outfile, d2mc.DEFAULT_VARS)
        self.assertEqual(command_used, shlex.join(argv))
        self.assertFalse(os.path.exists(os.path.join(os.path.dirname(outfile), "mon")))

    def test_varsAndYear(self):
        """A variable list (including a grid variable) and the time axis set to another year"""
        outfile, _, argv = self._run("--vars", f"{EXTRA_VAR},atmImp_lat", "--year", "2001")
        command_used = self._check_output(outfile, [EXTRA_VAR], year=2001)
        # The comma in the --vars list makes it through ncatted
        self.assertEqual(command_used, shlex.join(argv))

    def test_allVars(self):
        """All variables: time_bnds is renamed to climatology_bounds, lon/lat aren't averaged"""
        outfile, _, _ = self._run("--all-vars")
        self._check_output(outfile, d2mc.DEFAULT_VARS + [EXTRA_VAR])

    def test_keepAndReuse(self):
        """With --keep the monthly files are kept, and reused on a rerun"""
        outdir = tempfile.mkdtemp(dir=self._tempdir)
        outfile, _, _ = self._run("--keep", outdir=outdir)
        mondir = os.path.join(os.path.dirname(outfile), "mon")
        nyears = ENDYEAR - BEGYEAR + 1
        self.assertEqual(len(os.listdir(mondir)), 12 * nyears + 12)
        os.rename(outfile, os.path.join(self._tempdir, "first.nc"))
        outfile, printed, _ = self._run("--keep", outdir=outdir)
        self.assertEqual(printed.count("skipping ncra command"), 12 * nyears)
        self._check_output(outfile, d2mc.DEFAULT_VARS)

    def test_reuseWithDifferentVars_aborts(self):
        """Monthly files kept from a run with fewer variables abort a run with all variables"""
        outdir = tempfile.mkdtemp(dir=self._tempdir)
        outfile, _, _ = self._run("--keep", outdir=outdir)
        os.remove(outfile)
        with self.assertRaisesRegex(SystemExit, f"is missing variables: {EXTRA_VAR}"):
            self._run("--all-vars", outdir=outdir)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
