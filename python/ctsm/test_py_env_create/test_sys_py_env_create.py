#!/usr/bin/env python3

"""
System tests for py_env_create
"""
from __future__ import annotations

import os
import random
import string
import subprocess


import json
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from conda.testing.fixtures import CondaCLIFixture

# Allow test names that pylint doesn't like; otherwise hard to make them readable
# pylint: disable=invalid-name

pytest_plugins = "conda.testing.fixtures"


def get_unique_env_path(location, length):
    """
    Generates a unique name and path for a conda environment
    """
    env_name = "pec_test_" + "".join(random.choices(string.ascii_letters, k=length))
    env_path = os.path.join(location, env_name)
    if os.path.exists(env_path):
        env_name, env_path = get_unique_env_path(location, length)
    return env_name, env_path


def test_conda_create(conda_cli: CondaCLIFixture, tmp_path: Path):
    """
    A simple test to ensure that the Conda API is working right.
    """
    _, env_path = get_unique_env_path(tmp_path, 5)
    out, err, code = conda_cli("create", "--prefix", env_path, "--yes")

    assert f"conda activate {env_path}" in out
    assert not err  # no errors
    assert not code  # success!

    # verify everything worked using the `conda env list` command
    out, err, code = conda_cli("env", "list", "--json")
    assert any(Path(env_path).samefile(path) for path in json.loads(out).get("envs", []))
    assert not err  # no errors
    assert not code  # success!

    # cleanup, remove environment
    out, err, code = conda_cli("remove", "--all", "--prefix", env_path, "--yes")

    assert out
    assert not err  # no errors
    assert not code  # success!


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


def test_py_env_create(conda_cli: CondaCLIFixture):
    """
    A simple test to ensure that the Conda API is working right.
    """

    # Get path to py_env_create
    parent_dir = os.path.join(
        os.path.dirname(__file__),
        os.pardir,
        os.pardir,
        os.pardir,
    )
    py_env_create = os.path.join(parent_dir, "py_env_create")
    assert os.path.exists(py_env_create)

    # Get path to testing condafile
    condafile = os.path.join(parent_dir, "python", "empty.txt")
    assert os.path.exists(condafile)

    env_names = []

    # Run py_env_create once, making sure it was created
    env_names.append(get_unique_env_name(5))
    cmd = [py_env_create, "-n", env_names[0], "-f", condafile]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    out = subprocess.run(
        ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
    )
    assert any(
        os.path.split(path)[-1] == env_names[0] for path in json.loads(out.stdout).get("envs", [])
    )

    # Run py_env_create again, renaming existing
    env_names.append(env_names[0])
    env_names[0] += "_orig"
    cmd = [py_env_create, "-n", env_names[1], "-f", condafile, "-r", env_names[0]]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    # Ensure both exist now
    out = subprocess.run(
        ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
    )
    envs = json.loads(out.stdout).get("envs", [])
    for env_name in env_names:
        assert any(
            os.path.split(path)[-1] == env_name for path in envs
        )

    # Run py_env_create again, overwriting existing
    cmd = [py_env_create, "-n", env_names[1], "-f", condafile, "-o"]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    # Ensure both still exist
    out = subprocess.run(
        ["conda", "env", "list", "--json"], capture_output=True, text=True, check=True
    )
    envs = json.loads(out.stdout).get("envs", [])
    for env_name in env_names:
        assert any(
            os.path.split(path)[-1] == env_name for path in envs
        )

    # Remove test envs
    for env_name in env_names:
        conda_cli("remove", "--all", "--name", env_name, "--yes")
