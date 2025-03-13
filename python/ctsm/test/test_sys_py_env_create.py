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


def get_unique_env_name(length):
    """
    Generates a unique name for a conda environment
    """
    env_name = "pec_test_" + "".join(random.choices(string.ascii_letters, k=length))

    if any(os.path.split(path)[-1] == env_name for path in get_conda_envs()):
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
            if any(os.path.split(path)[-1] == env_name for path in envs):
                cmd = ["conda", "remove", "--all", "--name", env_name, "--yes"]
                subprocess.run(cmd, capture_output=True, text=True, check=False)

    def _create_empty_env(self, check=None, extra_args=None, expect_error=False, new_env_name=True):
        """Run py_env_create once, optionally making sure it was created"""
        if expect_error:
            if check:
                raise RuntimeError("check=True incompatible with expect_error=True")
            check = False

        if extra_args is not None and "-r" in extra_args:
            self.env_names.append(extra_args[extra_args.index("-r") + 1])
        if new_env_name:
            self.env_names.append(get_unique_env_name(5))

        # Form and run command
        cmd = [self.py_env_create, "-n", self.env_names[-1], "-f", self.empty_condafile, "--yes"]
        if extra_args:
            cmd += extra_args
        out = subprocess.run(cmd, capture_output=True, text=True, check=False)

        # Check result
        if expect_error:
            self.assertNotEqual(
                out.returncode, 0, "Unexpected success of py_env_create call:\n" + out.stdout
            )
        else:
            self.assertEqual(
                out.returncode, 0, "Unexpected failure of py_env_create call:\n" + out.stderr
            )

        if check:
            assert any(os.path.split(path)[-1] == self.env_names[-1] for path in get_conda_envs())

        return cmd, out

    def test_py_env_create_error_both_r_and_o(self):
        """
        Ensure py_env_create errors if both -r and -o are specified
        """

        # Run py_env_create
        _, out = self._create_empty_env(extra_args=["-o", "-r", "abc123"], expect_error=True)

        # Check error
        print(f"out.stderr: {out.stderr}")
        for x in out.stderr:
            print(x)
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
        self.assertTrue(f"Conda environment {self.env_names[0]} already exists." in out.stderr)
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
        self.assertTrue("Not going to let you rename the currently active conda env" in out.stderr)

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
        self.assertTrue(
            "Not going to let you overwrite the currently active conda env" in out.stderr
        )

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
        for env_name in self.env_names:
            assert any(os.path.split(path)[-1] == env_name for path in get_conda_envs())


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
