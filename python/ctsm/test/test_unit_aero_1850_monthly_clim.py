#!/usr/bin/env python3

"""
Unit tests for aero_1850_monthly_clim, and the CplhistMonthlyClim class methods it uses

The module-level functions of cplhist_monthly_clim are tested in test_unit_cplhist_monthly_clim.py
"""

import datetime
import io
import os
import shlex
import shutil
import sys
import tempfile
import unittest
from unittest import mock

from ctsm import unit_testing
from ctsm.cplhist import aero_1850_monthly_clim as aero
from ctsm.cplhist import cplhist_monthly_clim as mclim

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

CASE = "TESTCASE"
# The product class tested here
CLIM = aero.Aero1850MonthlyClim


def _check_args(args):
    """Run check_args for the parsed arguments"""
    CLIM(args).check_args()


#
# Tests for argument parsing
#
class TestGetParser(unittest.TestCase):
    """Tests of get_parser"""

    def test_defaults(self):
        """The defaults when only the positional arguments are given"""
        args = CLIM.get_parser().parse_args(["arch", CASE, "2001", "2002"])
        self.assertEqual(args.archdir, "arch")
        self.assertEqual(args.case, CASE)
        self.assertEqual((args.begyear, args.endyear), (2001, 2002))
        self.assertEqual(args.vars, CLIM.DEFAULT_VARS)
        self.assertFalse(args.all_vars)
        self.assertIsNone(args.outdir)
        self.assertEqual(args.year, CLIM.DEFAULT_YEAR)
        self.assertFalse(args.keep)
        self.assertIsNone(args.njobs)
        self.assertFalse(args.batch)
        self.assertIsNone(args.account)
        self.assertEqual(args.queue, mclim.DEFAULT_QUEUE)
        self.assertEqual(args.walltime, CLIM.DEFAULT_WALLTIME)

    def test_varsSplitsCommaList(self):
        """--vars is split on commas, with spaces stripped"""
        args = CLIM.get_parser().parse_args(["--vars", "a, b,c", "arch", CASE, "1", "2"])
        self.assertEqual(args.vars, ["a", "b", "c"])

    def test_allVarsWithVars_exits(self):
        """--all-vars and --vars can't be used together"""
        with mock.patch("sys.stderr", new_callable=io.StringIO):
            with self.assertRaises(SystemExit):
                CLIM.get_parser().parse_args(["-a", "--vars", "a", "arch", CASE, "1", "2"])


#
# Tests for the argument error checking
#
class TestCheckArgs(unittest.TestCase):
    """Tests of check_args"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()
        self._archdir = os.path.join(self._testdir, "archive")
        os.makedirs(mclim.input_dir(self._archdir, CASE))

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def _args(self, *options, archdir=None, years=("2001", "2002")):
        """Parse the options, with the archive directory, case and years"""
        if archdir is None:
            archdir = self._archdir
        return CLIM.get_parser().parse_args(list(options) + [archdir, CASE] + list(years))

    def test_validArgs(self):
        """Valid arguments pass"""
        _check_args(self._args())

    def test_yearOutOfRange(self):
        """BEGYEAR, ENDYEAR and YEAR must be from 1 to 9999"""
        with self.assertRaisesRegex(SystemExit, "ENDYEAR must be an integer from 1 to 9999"):
            _check_args(self._args(years=("2001", "10000")))
        with self.assertRaisesRegex(SystemExit, "BEGYEAR must be an integer from 1 to 9999"):
            _check_args(self._args(years=("0", "2001")))
        with self.assertRaisesRegex(SystemExit, "YEAR must be an integer from 1 to 9999"):
            _check_args(self._args("--year", "10000"))
        with self.assertRaisesRegex(SystemExit, "YEAR must be an integer from 1 to 9999"):
            _check_args(self._args("--year", "0"))

    def test_begyearAfterEndyear(self):
        """BEGYEAR can't be after ENDYEAR"""
        with self.assertRaisesRegex(SystemExit, "must not be after ENDYEAR"):
            _check_args(self._args(years=("2002", "2001")))

    def test_njobsZero(self):
        """--njobs must be positive"""
        with self.assertRaisesRegex(SystemExit, "NJOBS must be a positive integer"):
            _check_args(self._args("--njobs", "0"))

    def test_accountWithoutBatch(self):
        """--account is only allowed with --batch"""
        with self.assertRaisesRegex(SystemExit, "--account is only used with --batch"):
            _check_args(self._args("--account", "P1"))
        _check_args(self._args("--batch", "--account", "P1"))

    def test_onlyAlwaysVars(self):
        """--vars with only the ALWAYS_VARS leaves nothing to average"""
        with self.assertRaisesRegex(SystemExit, "No variables to average"):
            _check_args(self._args("--vars", ",".join(mclim.ALWAYS_VARS)))

    def test_missingDirectories(self):
        """The archive, output and input directories must exist"""
        with self.assertRaisesRegex(SystemExit, "Input archive directory .* does not exist"):
            _check_args(self._args(archdir=os.path.join(self._testdir, "dir_that_does_not_exist")))
        with self.assertRaisesRegex(SystemExit, "Output archive directory .* does not exist"):
            _check_args(
                self._args("--outdir", os.path.join(self._testdir, "dir_that_does_not_exist"))
            )
        shutil.rmtree(mclim.input_dir(self._archdir, CASE))
        with self.assertRaisesRegex(SystemExit, "Input case archive directory .* does not exist"):
            _check_args(self._args())

    def test_outputExists(self):
        """Abort if today's output file already exists"""
        outdir = mclim.output_dir(self._archdir, None, CASE)
        os.makedirs(outdir)
        outfile = mclim.output_filename(outdir, CASE, 2001, 2002, datetime.date.today())
        with open(outfile, "w", encoding="utf-8"):
            pass
        with self.assertRaisesRegex(SystemExit, "already exists"):
            _check_args(self._args())


#
# Tests of writing provenance information is correct
#
class TestWriteProvenance(unittest.TestCase):
    """Tests of write_provenance (with the ncatted command mocked)"""

    def test_oneNcattedCommand(self):
        """One ncatted --hst command sets the global attributes and deletes file_version"""
        # Replace run_cmd_output_on_error with a mock (run_cmd) that does nothing but record
        # the call, get_ctsm_git_short_hash to return a fixed hash, and sys.argv with the
        # command line given
        clim = CLIM(CLIM.get_parser().parse_args(["arch", CASE, "100", "130"]))
        clim.todaysdate = datetime.date(2026, 10, 3)
        with mock.patch.object(mclim, "run_cmd_output_on_error") as run_cmd, mock.patch.object(
            mclim, "get_ctsm_git_short_hash", return_value="abc1234"
        ), mock.patch.object(sys, "argv", ["aero_1850_monthly_clim", "--vars", "a,b", "x"]):
            clim.write_provenance("out.nc")
        run_cmd.assert_called_once()
        cmd = run_cmd.call_args.args[0]
        self.assertEqual(cmd[:3], ["ncatted", "-O", "--hst"])
        self.assertEqual(cmd[-1], "out.nc")
        attrs = cmd[4:-1:2]
        self.assertEqual(cmd[3:-1:2], ["-a"] * len(attrs))
        self.assertIn("Conventions,global,o,c,CF-1.11", attrs)
        self.assertIn("Created_on,global,o,c,2026-10-03", attrs)
        self.assertIn("Created_with,global,o,c,./aero_1850_monthly_clim -- abc1234", attrs)
        self.assertIn("command_used,global,o,c,aero_1850_monthly_clim --vars a,b x", attrs)
        self.assertIn(f"made_from_cesm_case,global,o,c,{CASE}", attrs)
        self.assertIn("file_version,global,d,,", attrs)
        title = [attr for attr in attrs if attr.startswith("title,")]
        self.assertTrue(title[0].endswith("for cyclical year 1850"))


#
# Test the batch options
#
class TestBatch(unittest.TestCase):
    """Tests of job_args and batch_script"""

    def test_jobArgs_default(self):
        """Only the positional arguments when the defaults are used"""
        args = CLIM.get_parser().parse_args(["arch", CASE, "2001", "2002"])
        self.assertEqual(CLIM(args).job_args(), ["arch", CASE, "2001", "2002"])

    def test_jobArgs_removes_batch_options(self):
        """Set all options, but show that job_args returns only the non-batch options"""
        args = CLIM.get_parser().parse_args(
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
            CLIM(args).job_args(),
            ["arch", CASE, "2001", "2002", "--vars", "a,b", "--outdir", "out", "--year", "2000"]
            + ["--keep", "--njobs", "8", "--verbose"],
        )

    def test_jobArgs_allVars(self):
        """Test that --all-vars is passed on when it's called as '-a'"""
        args = CLIM.get_parser().parse_args(["-a", "arch", CASE, "2001", "2002"])
        self.assertEqual(CLIM(args).job_args(), ["arch", CASE, "2001", "2002", "--all-vars"])

    def test_batchScript(self):
        """
        Test that batch_script creates: PBS settings, NCO, the directory and the quoted command
        """
        argv = ["arch", CASE, "2001", "2002", "--vars", "a,b", "--outdir", "my dir"]
        clim = CLIM(CLIM.get_parser().parse_args(["arch", CASE, "2001", "2002"]))
        script = clim.batch_script(
            argv, account="P1", walltime="01:00:00", workdir="/work dir", njobs=8
        )
        lines = script.splitlines()
        self.assertEqual(lines[0], "#!/bin/bash")
        self.assertIn(f"#PBS -N {CLIM.PBS_JOBNAME}", lines)
        self.assertIn("#PBS -A P1", lines)
        self.assertIn("#PBS -l walltime=01:00:00", lines)
        self.assertIn("mpiprocs=8", script)
        self.assertIn("module load nco", lines)
        self.assertIn("cd '/work dir'", lines)
        command = shlex.split(lines[-1])
        self.assertEqual(command[0], sys.executable)
        self.assertTrue(command[1].endswith(CLIM.WRAPPER))
        self.assertEqual(command[2:], argv)


#
# Test that the class attributes of a product subclass are used (this is what a new product,
# for example ozone, will set)
#
class _TestProductClim(mclim.CplhistMonthlyClim):
    """A product with different settings than the aerosol one"""

    DEFAULT_VARS = ["atmImp_Sa_o3"]
    DEFAULT_YEAR = 2000
    OUTPUT_TAG = "testclim"
    PBS_JOBNAME = "test_product_clim"
    WRAPPER = os.path.join("tools", "cplhist", "test_product_clim")
    DEFAULT_WALLTIME = "01:00:00"
    TITLE = "Test product title"
    COMMENT = "Test product comment"


class TestProductSubclass(unittest.TestCase):
    """Tests that a product subclass's class attributes flow through the base class"""

    def setUp(self):
        self._args = _TestProductClim.get_parser().parse_args(["arch", CASE, "2001", "2002"])
        self._clim = _TestProductClim(self._args)

    def test_parserDefaults(self):
        """The parser defaults come from the subclass"""
        self.assertEqual(self._args.vars, ["atmImp_Sa_o3"])
        self.assertEqual(self._args.year, 2000)
        self.assertEqual(self._args.walltime, "01:00:00")

    def test_jobArgs(self):
        """The subclass defaults aren't passed on to the batch job, the aerosol ones are"""
        self.assertEqual(self._clim.job_args(), ["arch", CASE, "2001", "2002"])
        args = _TestProductClim.get_parser().parse_args(
            ["--vars", ",".join(CLIM.DEFAULT_VARS), "-y", "1850", "arch", CASE, "2001", "2002"]
        )
        self.assertEqual(
            _TestProductClim(args).job_args(),
            ["arch", CASE, "2001", "2002", "--vars", ",".join(CLIM.DEFAULT_VARS)]
            + ["--year", "1850"],
        )

    def test_outputFile(self):
        """The output file name uses OUTPUT_TAG"""
        self._clim.todaysdate = datetime.date(2026, 10, 3)
        self.assertEqual(
            self._clim.output_file(),
            os.path.join(
                "arch", CASE, "proc", f"{CASE}.cpl.hx.atm.testclim.2001-2002_c20261003.nc"
            ),
        )

    def test_batchScript(self):
        """The batch script uses the subclass's job name and wrapper"""
        script = self._clim.batch_script(
            ["arch"], account="P1", walltime="01:00:00", workdir="/work", njobs=8
        )
        self.assertIn("#PBS -N test_product_clim", script.splitlines())
        self.assertTrue(shlex.split(script.splitlines()[-1])[1].endswith(_TestProductClim.WRAPPER))

    def test_provenanceAttrs(self):
        """The title, comment and Created_with come from the subclass"""
        with mock.patch.object(mclim, "get_ctsm_git_short_hash", return_value="abc1234"):
            attrs = dict(self._clim.provenance_attrs())
        self.assertEqual(attrs["title"], "Test product title for cyclical year 2000")
        self.assertEqual(attrs["comment"], "Test product comment")
        self.assertEqual(attrs["Created_with"], "./test_product_clim -- abc1234")

    def test_missingRequiredAttribute(self):
        """A subclass that doesn't set a required class attribute is rejected when it's defined"""
        with self.assertRaisesRegex(TypeError, "must set the class attributes: TITLE, COMMENT"):

            class _IncompleteClim(mclim.CplhistMonthlyClim):  # pylint: disable=unused-variable
                """A product missing TITLE and COMMENT"""

                DEFAULT_VARS = ["a"]
                PBS_JOBNAME = "incomplete"
                WRAPPER = "incomplete"


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
