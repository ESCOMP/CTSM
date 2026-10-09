#!/usr/bin/env python3
"""Unit tests for check-derecho-versions.py.

Run as `python docker/ctsm-ci-derecho-gnu/test_check_derecho_versions.py`, or
with `-v`. Standard library only (unittest, not pytest) -- same constraint as
the checker itself, so CI needs no environment beyond a python3.

Every test builds a self-contained fixture under a temporary directory: a fake
Dockerfile, a fake config_machines.xml, a fake intel_derecho.cmake, a fake ini
and a fake two-root module tree of empty .lua files. Nothing here reads /glade
or the repo's real files, so the tests say the same thing on a hosted runner as
on casper.

The fake config_machines.xml is trimmed from the real one but keeps every block
whose *exclusion* the checker depends on -- the intel/cray/nvhpc compiler
blocks, the mpilib="mpich" block and both mpilib="!mpi-serial" blocks -- with
versions that differ from the serial ones, so a broken filter shows up as a
failure rather than as an accidental pass.

Each test mutates one thing. Where a mutation should produce exactly one
failure, that count is asserted: the design calls for each check to report
distinctly rather than collapsing into a neighbour's message or passing.
"""

import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CHECKER_PATH = os.path.join(SCRIPT_DIR, "check-derecho-versions.py")


def _load_checker():
    """Import the checker despite its dashed, non-importable filename."""
    spec = importlib.util.spec_from_file_location(
        "check_derecho_versions", CHECKER_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


checker = _load_checker()

# MPICH_VERSION and PNETCDF_VERSION are present but no longer checked; they stay
# in the real Dockerfile until the image is stripped to serial-only, and their
# presence here guards against a check reappearing for them by accident.
BASE_DOCKERFILE = """\
FROM fedora:40
ARG GCC_VERSION=14.3.0
ARG MPICH_VERSION=3.4.3
ARG HDF5_VERSION=1.14.6
ARG NETCDF_C_VERSION=4.9.3
ARG NETCDF_FORTRAN_VERSION=4.6.2
ARG PNETCDF_VERSION=1.14.1
ARG ESMF_VERSION=8.9.1
ARG PFUNIT_VERSION=4.8.0
ARG MPI_SERIAL_VERSION=2.5.3
ARG PIO_VERSION=2.6.8
"""

# Trimmed from ccs_config/machines/derecho/config_machines.xml, with two
# deliberate tripwires so that a broken filter fails an ARG comparison rather
# than only the module-existence check -- the merge gate runs
# --skip-module-tree, so a filter tested only through the tree is untested
# where it matters most:
#
#   - the non-gnu compiler blocks carry *different* cray-libsci versions, so a
#     broken compiler filter makes cray-libsci resolve to conflicting versions
#     and the deviation guard goes red;
#   - the mpilib="!mpi-serial" block loads a module named plain `esmf` (derecho
#     loads `esmf-mpi` there), so a broken mpilib filter makes ESMF_VERSION
#     resolve to conflicting versions and goes red.
BASE_CONFIG = """\
<machine MACH="derecho">
    <COMPILERS>intel,gnu,nvhpc,cray</COMPILERS>
    <MPILIBS>mpich,openmpi</MPILIBS>
    <module_system type="module" allow_error="true">
      <modules>
        <command name="load">ncarenv/25.10</command>
        <command name="load">cesmdev/1.0</command>
        <command name="purge"/>
        <command name="load">conda/latest</command>
        <command name="load">nco</command>
        <command name="load">craype</command>
        <command name="load">cmake</command>
      </modules>
      <modules compiler="intel">
        <command name="load">intel/2025.3.2</command>
        <command name="load">mkl</command>
      </modules>
      <modules compiler="cray">
        <command name="load">cce/19.0.0</command>
        <command name="load">cray-libsci/24.01.0</command>
      </modules>
      <modules compiler="gnu">
        <command name="load">gcc/14.3.0</command>
        <command name="load">cray-libsci/25.03.0</command>
      </modules>
      <modules compiler="nvhpc">
        <command name="load">nvhpc/25.9</command>
        <command name="load">cray-libsci/23.09.0</command>
      </modules>
      <modules>
        <command name="load">ncarcompilers/1.2.0</command>
      </modules>
      <modules mpilib="mpich">
        <command name="load">cray-mpich/8.1.32</command>
      </modules>
      <modules mpilib="mpich" compiler="gnu" gpu_type="!none">
        <command name="load">cuda/12.9</command>
      </modules>
      <modules mpilib="mpi-serial" DEBUG="FALSE">
        <command name="load">netcdf/4.9.3</command>
        <command name="load">mpi-serial/2.5.3</command>
        <command name="load">parallelio-serial/2.6.8</command>
        <command name="load">esmf/8.9.1</command>
      </modules>
      <modules mpilib="mpi-serial" DEBUG="TRUE">
        <command name="load">netcdf/4.9.3</command>
        <command name="load">mpi-serial/2.5.3</command>
        <command name="load">parallelio-serial/2.6.8-debug</command>
        <command name="load">esmf/8.9.1</command>
      </modules>
      <modules mpilib="!mpi-serial" DEBUG="FALSE">
        <command name="load">netcdf-mpi/4.9.4</command>
        <command name="load">parallel-netcdf/1.14.1</command>
        <command name="load">parallelio/2.6.9</command>
        <command name="load">esmf-mpi/8.9.2</command>
        <command name="load">esmf/9.9.9</command>
      </modules>
      <modules mpilib="!mpi-serial" DEBUG="TRUE">
        <command name="load">netcdf-mpi/4.9.4-debug</command>
        <command name="load">parallel-netcdf/1.14.1-debug</command>
        <command name="load">parallelio/2.6.9-debug</command>
        <command name="load">esmf-mpi/8.9.2-debug</command>
      </modules>
    </module_system>
</machine>
"""

BASE_INI = """\
[snapshot]
hdf5 = 1.14.6
netcdf_fortran = 4.6.2
measured_against_netcdf = 4.9.3

[deviation_guard]
cray_libsci = 25.03.0
"""

BASE_CMAKE = """\
if (MPILIB STREQUAL mpi-serial)
endif()
set(PFUNIT_PATH "$ENV{CESMDATAROOT}/tools/pFUnit/\\
pFUnit4.8.0_derecho_Intel2023.2.1_noMPI_noOpenMP")
"""

# The modulefiles a gnu + mpi-serial run would find on MODULEPATH, split across
# the two roots as derecho has them today -- minus conda/latest, which derecho
# does have but UNCHECKABLE_PINS excludes. Leaving it out means every test here
# would go red if the exclusion stopped working.
BASE_SYSTEM_MODULES = [
    "environment/ncarenv/25.10.lua",
    "environment/cesmdev/1.0.lua",
    "25.10/Core/gcc/14.3.0.lua",
    "25.10/gcc/14.3.0/ncarcompilers/1.2.0.lua",
    "25.10/gcc/14.3.0/cray-libsci/25.03.0.lua",
    "25.10/gcc/14.3.0/netcdf/4.9.3.lua",
    "25.10/gcc/14.3.0/mpi-serial/2.5.3.lua",
    "25.10/gcc/14.3.0/parallelio-serial/2.6.8.lua",
    "25.10/gcc/14.3.0/esmf/8.9.1.lua",
]
# Only here, never under the system root -- the real arrangement, and the reason
# the checker searches two roots.
BASE_CSEG_MODULES = [
    "25.10/gcc/14.3.0/parallelio-serial/2.6.8-debug.lua",
]


def failures(lines):
    return [line for line in lines if "❌" in line]


class CheckerTestCase(unittest.TestCase):
    """Builds the fixture; each test mutates an attribute before running."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmpdir = tmp.name
        self.dockerfile_text = BASE_DOCKERFILE
        self.config_text = BASE_CONFIG
        self.ini_text = BASE_INI
        self.cmake_text = BASE_CMAKE
        self.system_modules = list(BASE_SYSTEM_MODULES)
        self.cseg_modules = list(BASE_CSEG_MODULES)
        # Flipping either of these to False leaves the root nonexistent, which
        # is how the checker sees a runner with no /glade.
        self.make_system_root = True
        self.make_cseg_root = True
        # Names of inputs to leave unwritten, for the setup-problem cases.
        self.omit = set()

    def _write(self, path, text):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)

    def _write_tree(self, root, relpaths, create_root):
        if not create_root:
            return root
        os.makedirs(root, exist_ok=True)
        for rel in relpaths:
            self._write(os.path.join(root, rel), "-- fake modulefile\n")
        return root

    def run_checks(self, **kwargs):
        dockerfile = os.path.join(self.tmpdir, "Dockerfile")
        config_xml = os.path.join(self.tmpdir, "config_machines.xml")
        ini = os.path.join(self.tmpdir, "derecho-versions.ini")
        cmake = os.path.join(self.tmpdir, "intel_derecho.cmake")
        for name, path, text in (
            ("dockerfile", dockerfile, self.dockerfile_text),
            ("config_xml", config_xml, self.config_text),
            ("snapshot", ini, self.ini_text),
            ("pfunit_cmake", cmake, self.cmake_text),
        ):
            if name not in self.omit:
                self._write(path, text)
        system_root = self._write_tree(
            os.path.join(self.tmpdir, "sysmodules"),
            self.system_modules,
            self.make_system_root,
        )
        cseg_root = self._write_tree(
            os.path.join(self.tmpdir, "csegmodules"),
            self.cseg_modules,
            self.make_cseg_root,
        )
        return checker.run_checks(
            dockerfile=dockerfile,
            config_xml=config_xml,
            pfunit_cmake=cmake,
            snapshot=ini,
            module_roots=(system_root, cseg_root),
            **kwargs,
        )


class TestBaseline(CheckerTestCase):
    def test_everything_agreeing_passes(self):
        ok, lines = self.run_checks()
        self.assertEqual(failures(lines), [])
        self.assertTrue(ok)

    def test_debug_pio_resolves_in_the_cseg_root(self):
        """parallelio-serial/<ver>-debug exists only under the CSEG root.

        Regression guard, not a hypothetical: this is derecho's real layout, and
        a checker that walks only the system root fails on its first run.
        """
        ok, lines = self.run_checks()
        self.assertTrue(ok)
        resolved = [
            line for line in lines if "parallelio-serial/2.6.8-debug" in line
        ]
        self.assertEqual(len(resolved), 1)
        self.assertIn("csegmodules", resolved[0])

    def test_mpi_stack_modules_are_ignored(self):
        """The MPI flavor is not the container's subject and must not leak in.

        Asserted in both modes, because they fail differently. With the tree
        checked, a leaked module goes red by resolving nowhere; with the tree
        skipped -- which is how the merge gate runs -- nothing but an ARG
        comparison can catch it, which is what the fixture's `esmf/9.9.9` in the
        !mpi-serial block is for.
        """
        self.config_text = (
            self.config_text.replace("netcdf-mpi/4.9.4", "netcdf-mpi/9.9.9")
            .replace("esmf-mpi/8.9.2", "esmf-mpi/9.9.9")
            .replace("parallel-netcdf/1.14.1", "parallel-netcdf/9.9.9")
            .replace("parallelio/2.6.9", "parallelio/9.9.9")
            .replace("cray-mpich/8.1.32", "cray-mpich/9.9.9")
        )
        for skip in (False, True):
            with self.subTest(skip_module_tree=skip):
                ok, lines = self.run_checks(skip_module_tree=skip)
                self.assertEqual(failures(lines), [])
                self.assertTrue(ok)

    def test_conda_latest_is_skipped_and_said_so(self):
        """"latest" is not a version, so resolving it would prove nothing.

        Reported, not dropped: a list that quietly got shorter is the failure
        mode this check exists against.
        """
        ok, lines = self.run_checks()
        self.assertTrue(ok)
        skipped = [line for line in lines if "conda/latest" in line]
        self.assertEqual(len(skipped), 1)
        self.assertIn("not checked", skipped[0])
        self.assertNotIn("✅", skipped[0])
        self.assertNotIn("❌", skipped[0])

    def test_a_real_conda_version_is_still_existence_checked(self):
        """The exclusion keys on "conda/latest", not on "conda".

        Keyed on the name, a derecho that started pinning conda/<version> would
        silently stop being checked -- which is the opposite of the point.
        """
        self.config_text = self.config_text.replace(
            "conda/latest", "conda/25.1.0"
        )
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("conda/25.1.0", bad[0])

    def test_unpinned_loads_are_not_existence_checked(self):
        """nco, craype and cmake carry no version, so nothing can go stale."""
        ok, lines = self.run_checks()
        self.assertTrue(ok)
        for name in ("nco", "craype", "cmake"):
            self.assertEqual(
                [line for line in lines if line.startswith(f"✅ module {name}/")],
                [],
            )


class TestArgChecks(CheckerTestCase):
    def test_perturbed_arg_fails_alone(self):
        self.dockerfile_text = self.dockerfile_text.replace(
            "NETCDF_C_VERSION=4.9.3", "NETCDF_C_VERSION=4.9.2"
        )
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("NETCDF_C_VERSION", bad[0])

    def test_missing_arg_fails(self):
        self.dockerfile_text = self.dockerfile_text.replace(
            "ARG ESMF_VERSION=8.9.1\n", ""
        )
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("ESMF_VERSION", bad[0])
        self.assertIn("ARG not found", bad[0])

    def test_perturbed_config_module_fails(self):
        """esmf, not netcdf: perturbing netcdf also trips the provenance check.

        The module tree gains the new version too, so the single ❌ here is the
        direct comparison and not the existence check.
        """
        self.config_text = self.config_text.replace("esmf/8.9.1", "esmf/8.9.2")
        self.system_modules.append("25.10/gcc/14.3.0/esmf/8.9.2.lua")
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("ESMF_VERSION", bad[0])
        self.assertIn("esmf/8.9.2", bad[0])

    def test_pfunit_mismatch_fails(self):
        self.cmake_text = self.cmake_text.replace("pFUnit4.8.0", "pFUnit4.9.0")
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("PFUNIT_VERSION", bad[0])


class TestSnapshot(CheckerTestCase):
    def test_perturbed_snapshot_value_fails(self):
        self.ini_text = self.ini_text.replace("hdf5 = 1.14.6", "hdf5 = 1.14.5")
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("HDF5_VERSION", bad[0])

    def test_stale_provenance_refuses_both_snapshot_args(self):
        """A snapshot measured against a superseded netcdf is not comparable.

        Three distinct failures, not one: the provenance itself, and each ARG
        reported as *not compared* rather than silently matching a stale number.
        """
        self.ini_text = self.ini_text.replace(
            "measured_against_netcdf = 4.9.3", "measured_against_netcdf = 4.9.2"
        )
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 3)
        self.assertIn("[snapshot] provenance", bad[0])
        self.assertIn("HDF5_VERSION", bad[1])
        self.assertIn("not compared", bad[1])
        self.assertIn("NETCDF_FORTRAN_VERSION", bad[2])
        self.assertIn("not compared", bad[2])


class TestDeviationGuard(CheckerTestCase):
    def test_cray_libsci_move_fails(self):
        """The guard has no version ARG; it compares the config to the ini."""
        self.config_text = self.config_text.replace(
            "cray-libsci/25.03.0", "cray-libsci/25.04.0"
        )
        self.system_modules.append("25.10/gcc/14.3.0/cray-libsci/25.04.0.lua")
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("cray-libsci", bad[0])
        self.assertIn("25.03.0", bad[0])


class TestModuleTree(CheckerTestCase):
    def test_missing_pinned_module_fails(self):
        self.system_modules.remove("25.10/gcc/14.3.0/esmf/8.9.1.lua")
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("esmf", bad[0])
        self.assertIn("8.9.1", bad[0])

    def test_module_under_cray_mpich_branch_does_not_count(self):
        """That directory only joins MODULEPATH once cray-mpich is loaded.

        A serial run never loads it, so finding the module there would hide
        exactly the drift this check exists for.
        """
        self.system_modules.remove("25.10/gcc/14.3.0/esmf/8.9.1.lua")
        self.system_modules.append(
            "25.10/cray-mpich/8.1.32/gcc/14.3.0/esmf/8.9.1.lua"
        )
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("esmf", bad[0])
        self.assertIn("8.9.1", bad[0])

    def test_unreadable_system_root_fails_loudly(self):
        """No /glade is a failure, not a silent skip and not an exception."""
        self.make_system_root = False
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 1)
        self.assertIn("sysmodules", bad[0])
        self.assertIn("--skip-module-tree", bad[0])
        # Per-module lines would all say the same thing, so they are suppressed.
        self.assertEqual(
            [line for line in lines if line.startswith("✅ module ")], []
        )

    def test_unreadable_cseg_root_still_reports_per_module(self):
        """The system root alone still answers truthfully for most modules."""
        self.make_cseg_root = False
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertIn("csegmodules", bad[0])
        # The -debug variant lives only in the CSEG root, so it is named as
        # missing rather than the whole tier going unreported.
        self.assertTrue(
            any("parallelio-serial/2.6.8-debug" in line for line in bad)
        )
        self.assertTrue(
            any(line.startswith("✅ module netcdf/4.9.3") for line in lines)
        )

    def test_skip_module_tree_says_not_checked(self):
        ok, lines = self.run_checks(skip_module_tree=True)
        self.assertTrue(ok)
        self.assertTrue(any("NOT CHECKED" in line for line in lines))
        self.assertEqual(
            [line for line in lines if line.startswith("✅ module ")], []
        )


class TestSetupProblemsRaise(CheckerTestCase):
    """Unreadable inputs raise; only drift is reported as a ❌ line.

    The distinction is load-bearing. A ❌ asserts that derecho and the
    Dockerfile disagree, and a run that could not read its inputs has
    established no such thing -- it must stop, not report. A refactor that
    caught these and returned (False, [...]) would pass every other test here.
    """

    def test_missing_dockerfile_raises(self):
        self.omit = {"dockerfile"}
        with self.assertRaises(FileNotFoundError):
            self.run_checks()

    def test_missing_config_machines_raises(self):
        self.omit = {"config_xml"}
        with self.assertRaises(FileNotFoundError):
            self.run_checks()

    def test_missing_ini_raises(self):
        self.omit = {"snapshot"}
        with self.assertRaises(FileNotFoundError):
            self.run_checks()

    def test_missing_pfunit_cmake_raises(self):
        self.omit = {"pfunit_cmake"}
        with self.assertRaises(FileNotFoundError) as caught:
            self.run_checks()
        message = str(caught.exception)
        self.assertIn("intel_derecho.cmake", message)
        self.assertIn("bin/git-fleximod update ccs_config", message)

    def test_missing_pfunit_cmake_raises_even_when_no_check_reads_it(self):
        """Pins the hoist: the file's existence is asserted before the loop.

        With PFUNIT_VERSION absent from the Dockerfile the pfunit check never
        opens intel_derecho.cmake, so a version that raised only where the file
        is read would return a ❌ report for a setup problem instead.
        """
        self.omit = {"pfunit_cmake"}
        self.dockerfile_text = self.dockerfile_text.replace(
            "ARG PFUNIT_VERSION=4.8.0\n", ""
        )
        with self.assertRaises(FileNotFoundError):
            self.run_checks()

    def test_config_without_module_system_raises(self):
        self.config_text = '<machine MACH="derecho"><OS>CNL</OS></machine>\n'
        with self.assertRaises(ValueError):
            self.run_checks()

    def test_missing_ini_section_raises_naming_the_file(self):
        """Validated up front, so the error names the file and the keys.

        Left to configparser inside the loop, this surfaced as a bare
        NoSectionError naming neither, and discarded every line gathered so far.
        """
        self.ini_text = self.ini_text.split("[deviation_guard]")[0]
        with self.assertRaises(ValueError) as caught:
            self.run_checks()
        message = str(caught.exception)
        self.assertIn("derecho-versions.ini", message)
        self.assertIn("[deviation_guard] cray_libsci", message)


class TestUnsafeModuleNames(CheckerTestCase):
    """A load command naming a path must not be searched as a module.

    os.path.join drops everything before an absolute component, so an absolute
    "name" would collapse every searched path onto itself and could report a
    file outside derecho's tree as resolved. Needs a garbled config_machines.xml
    to happen, so this is robustness, not exposure.
    """

    def test_path_like_module_names_are_reported_not_searched(self):
        marker = '<command name="load">esmf/8.9.1</command>'
        self.config_text = self.config_text.replace(
            marker,
            marker
            + '\n<command name="load">/elsewhere/escaped/1.0</command>'
            + '\n<command name="load">../escaped/2.0</command>',
            1,
        )
        ok, lines = self.run_checks()
        self.assertFalse(ok)
        bad = failures(lines)
        self.assertEqual(len(bad), 2)
        self.assertTrue(any("is an absolute path" in line for line in bad))
        self.assertTrue(any('contains ".."' in line for line in bad))
        # Reported, and never reported as resolved.
        self.assertEqual([line for line in lines if "escaped" in line], bad)


class TestExitCodes(unittest.TestCase):
    """The exit code is the only thing CI reads, so run the real script.

    Every other test calls run_checks() directly and would survive a tail
    rewritten as `raise SystemExit(main())` -- an inverted gate that reports
    every mismatch as success while the suite stays green. These cases run the
    script in a subprocess and assert the codes, and in passing they pin the
    --skip-module-tree spelling that the workflow hardcodes.

    The fixture is a miniature repo rather than the real one: the checker
    derives every input path from its own location, so copying it beside a fake
    Dockerfile and ini, with a fake ccs_config two levels up, redirects all of
    them without the script needing a test-only way to be told where to look.
    """

    def _fixture_script(self, ini_text=BASE_INI, omit=()):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        script_dir = os.path.join(tmp.name, "docker", "ctsm-ci-derecho-gnu")
        config_dir = os.path.join(tmp.name, "ccs_config", "machines", "derecho")
        os.makedirs(script_dir)
        os.makedirs(config_dir)
        shutil.copy(CHECKER_PATH, script_dir)
        for name, path, text in (
            ("Dockerfile", os.path.join(script_dir, "Dockerfile"),
             BASE_DOCKERFILE),
            ("derecho-versions.ini",
             os.path.join(script_dir, "derecho-versions.ini"), ini_text),
            ("config_machines.xml",
             os.path.join(config_dir, "config_machines.xml"), BASE_CONFIG),
            ("intel_derecho.cmake",
             os.path.join(config_dir, "intel_derecho.cmake"), BASE_CMAKE),
        ):
            if name in omit:
                continue
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
        return os.path.join(script_dir, os.path.basename(CHECKER_PATH))

    def _run(self, script):
        # --skip-module-tree so the fixture needs no /glade and these cases say
        # the same thing on a hosted runner as on casper.
        return subprocess.run(
            [sys.executable, script, "--skip-module-tree"],
            capture_output=True,
            text=True,
            check=False,
        )

    def test_exits_zero_when_everything_agrees(self):
        proc = self._run(self._fixture_script())
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertIn("OK: all version checks passed", proc.stdout)
        self.assertIn("NOT CHECKED", proc.stdout)

    def test_exits_one_on_a_mismatch(self):
        proc = self._run(
            self._fixture_script(
                ini_text=BASE_INI.replace("hdf5 = 1.14.6", "hdf5 = 1.14.5")
            )
        )
        self.assertEqual(proc.returncode, 1, proc.stdout + proc.stderr)
        self.assertIn("FAILED: version mismatch", proc.stdout)
        self.assertIn("HDF5_VERSION", proc.stdout)

    def test_setup_problem_does_not_masquerade_as_a_mismatch(self):
        proc = self._run(self._fixture_script(omit=("Dockerfile",)))
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("FileNotFoundError", proc.stderr)
        self.assertNotIn("FAILED: version mismatch", proc.stdout)


if __name__ == "__main__":
    unittest.main()
