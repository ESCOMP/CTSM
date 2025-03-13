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


def get_unique_env_name(length):
    """
    Generates a unique name for a conda environment
    """
    env_name = "pec_test_" + "".join(random.choices(string.ascii_letters, k=length))
    out = subprocess.run(
        ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
    )
    if any(os.path.split(path)[-1] == env_name for path in json.loads(out.stdout).get("envs", [])):
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
        for env_name in self.env_names:
            cmd = ["conda", "remove", "--all", "--name", env_name, "--yes"]
            subprocess.run(cmd, capture_output=True, text=True, check=True)

    def test_py_env_create(self):
        """
        A few calls of py_env_create to ensure it's working right.
        """

        # Run py_env_create once, making sure it was created
        self.env_names.append(get_unique_env_name(5))
        cmd = [self.py_env_create, "-n", self.env_names[0], "-f", self.empty_condafile, "--yes"]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        out = subprocess.run(
            ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
        )
        assert any(
            os.path.split(path)[-1] == self.env_names[0]
            for path in json.loads(out.stdout).get("envs", [])
        )

        # Run py_env_create again, renaming existing
        self.env_names.append(self.env_names[0])
        self.env_names[0] += "_orig"
        cmd = [
            self.py_env_create,
            "-n",
            self.env_names[1],
            "-f",
            self.empty_condafile,
            "-r",
            self.env_names[0],
            "--yes",
        ]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(e.stderr) from e
        # Ensure both exist now
        out = subprocess.run(
            ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
        )
        envs = json.loads(out.stdout).get("envs", [])
        for env_name in self.env_names:
            assert any(os.path.split(path)[-1] == env_name for path in envs)

        # Run py_env_create again, overwriting existing
        cmd = [
            self.py_env_create,
            "-n",
            self.env_names[1],
            "-f",
            self.empty_condafile,
            "-o",
            "--yes",
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Ensure both still exist
        out = subprocess.run(
            ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
        )
        envs = json.loads(out.stdout).get("envs", [])
        for env_name in self.env_names:
            assert any(os.path.split(path)[-1] == env_name for path in envs)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
