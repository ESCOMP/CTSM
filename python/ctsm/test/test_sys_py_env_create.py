#!/usr/bin/env python3

"""
System tests for py_env_create
"""
from __future__ import annotations

import os
import unittest
import random
import string
import subprocess
import json

from ctsm import unit_testing
from ctsm.path_utils import path_to_ctsm_root

# Allow test names that pylint doesn't like; otherwise hard to make them readable
# pylint: disable=invalid-name


def get_conda_envs():
    """
    List the user's conda environments
    """
    try:
        out = subprocess.run(
            ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(e.stderr) from e
    return json.loads(out.stdout).get("envs", [])


def does_env_exist(this_env, env_paths):
    """
    Checks whether the user has a certain conda env
    """
    return any(os.path.split(path)[-1] == this_env for path in env_paths)


def get_unique_env_name(length):
    """
    Generates a unique name for a conda environment
    """
    env_name = "pec_test_" + "".join(random.choices(string.ascii_letters, k=length))

    if does_env_exist(env_name, get_conda_envs()):
        env_name = get_unique_env_name(length)
    return env_name


class TestSysPyEnvCreate(unittest.TestCase):
    """System tests for py_env_create"""

    def setUp(self):
        """Run this before each test"""

        # Get path to py_env_create
        self.py_env_create = os.path.join(path_to_ctsm_root(), "py_env_create")
        assert os.path.exists(self.py_env_create)

        # Get path to testing condafile
        self.empty_condafile = os.path.join(path_to_ctsm_root(), "python", "empty.txt")
        assert os.path.exists(self.empty_condafile)

        # Set up other variables
        self.env_names = []

    def tearDown(self):
        """Run this after each test, whether it succeeded or failed"""
        # Remove test envs
        envs = get_conda_envs()
        for env_name in self.env_names:
            if does_env_exist(env_name, envs):
                cmd = ["conda", "remove", "--all", "--name", env_name, "--yes"]
                subprocess.run(cmd, capture_output=True, text=True, check=False)

    def _create_empty_env(self, check=None, extra_args=None, expect_error=False, new_env_name=True):
        """Run py_env_create once, optionally making sure it was created"""
        if expect_error:
            if check:
                raise RuntimeError("check=True incompatible with expect_error=True")
            check = False
        elif check is None:
            check = True

        if new_env_name:
            self.env_names.append(get_unique_env_name(5))

        # Form and run command
        cmd = [self.py_env_create, "-n", self.env_names[-1], "-f", self.empty_condafile, "--yes"]
        if extra_args:
            cmd += extra_args
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)

        print("cmd: " + " ".join(cmd))
        print("STDOUT:\n" + out.stdout)
        print("STDERR:\n" + out.stderr)
        print(f"RETURN CODE: {out.returncode}")

        # Check result
        if expect_error:
            self.assertNotEqual(
                out.returncode, 0, "Unexpected success of py_env_create call:\n" + out.stdout
            )
        else:

            self.assertEqual(
                out.returncode, 0, "Unexpected failure of py_env_create call:\n" + out.stderr
            )

        # Check that env was successfully created
        if check:
            assert does_env_exist(self.env_names[-1], get_conda_envs())

        # If you renamed something, make sure it gets deleted in tearDown
        if extra_args is not None and "-r" in extra_args:
            self.env_names.append(extra_args[extra_args.index("-r") + 1])

        return cmd, out

    def _check_dry_run(self, stdout):
        self.assertTrue("to install the python environment" in stdout)
        self.assertTrue("Using the file" in stdout)
        self.assertTrue(
            stdout.endswith("Exiting before any actions taken, due to -d/--dry-run option.\n")
        )

    def test_py_env_create_dry_run_short(self):
        """
        Ensure py_env_create doesn't do anything if -d given
        """

        # Run py_env_create
        _, out = self._create_empty_env(extra_args=["-d"], check=False)

        # Check stdout
        self._check_dry_run(out.stdout)

    def test_py_env_create_dry_run_long(self):
        """
        Ensure py_env_create doesn't do anything if --dry-run given
        """

        # Run py_env_create
        _, out = self._create_empty_env(extra_args=["--dry-run"], check=False)

        # Check stdout
        self._check_dry_run(out.stdout)

    def test_py_env_create_dry_run_o(self):
        """
        Ensure py_env_create doesn't do anything if -d and -o given
        """

        # Create env
        self._create_empty_env(check=True)

        # Run py_env_create in dry run mode with overwrite
        _, out = self._create_empty_env(extra_args=["-d", "-o"], check=False, new_env_name=False)

        # Check stdout
        self._check_dry_run(out.stdout)
        self.assertTrue("Overwriting existing" in out.stdout)

    def test_py_env_create_dry_run_r(self):
        """
        Ensure py_env_create doesn't do anything if -d and -r given
        """

        # Create env
        self._create_empty_env(check=True)

        # Run py_env_create in dry run mode with rename existing
        _, out = self._create_empty_env(
            extra_args=["-d", "-r", "abc123"], check=False, new_env_name=False
        )

        # Check stdout
        self._check_dry_run(out.stdout)
        self.assertTrue("Renaming existing" in out.stdout)

    def test_py_env_create_help_r(self):
        """
        Ensure that name of new env is in the help text for -r/--rename
        """

        # Run py_env_create in help mode with rename existing
        _, out = self._create_empty_env(extra_args=["-h"], check=False)

        # Check stdout
        expected_str = f"the environment you're installing ({self.env_names[-1]})"
        self.assertTrue(expected_str in out.stdout)

    def test_py_env_create_old(self):
        """
        Ensure py_env_create works with --old
        """

        # Generate env name
        self.env_names.append(get_unique_env_name(5))

        # Run py_env_create
        cmd = [self.py_env_create, "--old", "-n", self.env_names[-1], "--dry-run"]
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)
        print("STDOUT")
        print(out.stdout)
        print("STDERR")
        print(out.stderr)

        # Check stdout
        self._check_dry_run(out.stdout)

    def test_py_env_create_error_both_r_and_o(self):
        """
        Ensure py_env_create errors if both -r and -o are specified
        """

        # Run py_env_create
        _, out = self._create_empty_env(extra_args=["-o", "-r", "abc123"], expect_error=True)

        # Check error
        self.assertTrue("Only specify one of -o/--overwrite or -r/--rename-existing." in out.stderr)

    def test_py_env_create_error_both_r_and_o_longnames(self):
        """
        Ensure py_env_create errors if both --rename-existing and --overwrite are specified
        """

        # Run py_env_create
        _, out = self._create_empty_env(
            extra_args=["--overwrite", "--rename-existing", "abc123"],
            expect_error=True,
        )

        # Check error
        self.assertTrue("Only specify one of -o/--overwrite or -r/--rename-existing." in out.stderr)

    def test_py_env_create_error_exists_without_r_or_o(self):
        """
        Ensure py_env_create errors if environment already exists without -o or -r
        """

        # Run py_env_create once, making sure it was created
        self._create_empty_env()

        # Try doing it again without specifying -o or -r
        _, out = self._create_empty_env(expect_error=True, new_env_name=False)

        # Check error
        self.assertTrue(f"Environment {self.env_names[0]} already exists." in out.stderr)
        self.assertTrue("Try again using one of:" in out.stderr)

    def test_py_env_create_error_exists_and_so_does_r(self):
        """
        Ensure py_env_create errors if environment already exists and -r name does too
        """

        # Run py_env_create once, making sure it was created
        self._create_empty_env()

        # Try doing it again with an existing name in -r
        _, out = self._create_empty_env(
            expect_error=True, new_env_name=False, extra_args=["-r", self.env_names[0]]
        )

        # Check error
        self.assertTrue(f"{self.env_names[0]} also already exists" in out.stderr)

    def test_py_env_create_error_renaming_current(self):
        """
        Ensure py_env_create errors if trying to rename current env
        """

        # Run py_env_create once, making sure it was created
        cmd, _ = self._create_empty_env()

        # Try doing it again in that conda env with its name in -r
        cmd = ["conda", "run", "-n", self.env_names[0]] + cmd + ["-r", self.env_names[0]]
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)

        # Check error
        self.assertTrue("Not going to let you rename the currently active env" in out.stderr)

    def test_py_env_create_error_overwriting_current(self):
        """
        Ensure py_env_create errors if trying to overwrite current env
        """

        # Run py_env_create once, making sure it was created
        cmd, _ = self._create_empty_env()

        # Try doing it again in that conda env with -o
        cmd = ["conda", "run", "-n", self.env_names[0]] + cmd + ["-o"]
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)

        # Check error
        self.assertTrue("Not going to let you overwrite the currently active env" in out.stderr)

    def test_complete_py_env_create(self):
        """
        A few calls of py_env_create to ensure it's working right.
        """

        # Run py_env_create once, making sure it was created
        self._create_empty_env(check=True)

        # Run py_env_create again, overwriting existing
        _, out = self._create_empty_env(new_env_name=False, extra_args=["-o"], check=True)
        # Ensure we only added one
        try:
            self.assertEqual(len(self.env_names), 1)
        except AssertionError as e:
            print(out.stdout)
            raise e

        # Run py_env_create again, renaming existing
        _, out = self._create_empty_env(
            new_env_name=False, extra_args=["-r", self.env_names[-1] + "_orig"]
        )
        # Ensure both exist now
        try:
            self.assertEqual(len(self.env_names), 2)
        except AssertionError as e:
            print(out.stdout)
            raise e
        env_list = get_conda_envs()
        for env_name in self.env_names:
            assert does_env_exist(env_name, env_list)

    def test_complete_py_env_create_mamba(self):
        """
        A few calls of py_env_create to ensure it's working right with mamba instead of conda
        """

        # Run py_env_create once, making sure it was created
        self._create_empty_env(check=True, extra_args=["-m"])

        # Run py_env_create again, overwriting existing
        _, out = self._create_empty_env(new_env_name=False, extra_args=["-o", "-m"], check=True)
        # Ensure we only added one
        try:
            self.assertEqual(len(self.env_names), 1)
        except AssertionError as e:
            print(out.stdout)
            raise e

        # Run py_env_create again, renaming existing
        _, out = self._create_empty_env(
            new_env_name=False, extra_args=["-r", self.env_names[-1] + "_orig", "--mamba"]
        )
        # Ensure both exist now
        try:
            self.assertEqual(len(self.env_names), 2)
        except AssertionError as e:
            print(out.stdout)
            raise e
        env_list = get_conda_envs()
        for env_name in self.env_names:
            assert does_env_exist(env_name, env_list)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
