#!/usr/bin/env python3
"""Check that the ctsm-ci-derecho-gnu Dockerfile versions match derecho's gnu stack.

The container replicates derecho's gnu software stack (declared in
ccs_config/machines/derecho/config_machines.xml). This check compares the
Dockerfile's version ARGs to that config, in three modes:

  direct    - the ARG must equal the derecho gnu module version, read live from
              config_machines.xml (gcc, netcdf-mpi, parallel-netcdf, esmf-mpi,
              mpi-serial, parallelio-serial).
  deviation - the ARG is an intentional open-source stand-in that is NOT
              compared to derecho; instead the derecho module version (read
              live) must still equal a recorded value, so a derecho change
              trips the check (MPICH_VERSION <-> cray-mpich).
  snapshot  - the derecho version has NO standalone entry in
              config_machines.xml -- netCDF-Fortran is bundled inside the
              netcdf-mpi module, HDF5 is the hdf5-mpi module netcdf-mpi pulls
              in via depends_on -- so the ARG is compared to a hand-recorded
              value (HDF5, netCDF-Fortran).
  pfunit    - derecho has no pFUnit module at all; the version is embedded in
              the PFUNIT_PATH set by intel_derecho.cmake, read live from there.

derecho's gnu stack comes in two flavors -- the MPI one (MPILIB=mpich, the
machine default, asserted against <MPILIBS>) and the serial one
(MPILIB=mpi-serial) -- and since ccs_config_cesm1.0.88 the same library can
appear in both under different module names and versions (esmf-mpi vs esmf,
parallelio vs parallelio-serial). Every module-reading check therefore names the
flavor it is asking about.

Where one ARG builds both flavors in the container -- ESMF_VERSION and
NETCDF_C_VERSION each compile twice, once MPI-linked and once serial -- the
check also asserts derecho's two modules still agree with each other. If they
ever diverge, no single ARG can match both and the Dockerfile needs a second
one, so that is reported instead of an arbitrary half-truth.

Recorded values (snapshot + deviation guard) live in derecho-versions.ini.

Standard library only. Needs ccs_config populated (bin/git-fleximod update
ccs_config). A version mismatch prints a per-component report and exits 1; a
setup problem (missing/malformed inputs) raises an exception.
"""

import configparser
import os
import re
import xml.etree.ElementTree as ET

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DOCKERFILE = os.path.join(SCRIPT_DIR, "Dockerfile")
SNAPSHOT = os.path.join(SCRIPT_DIR, "derecho-versions.ini")
CONFIG_XML = os.path.join(
    REPO_ROOT, "ccs_config", "machines", "derecho", "config_machines.xml"
)
# derecho defines PFUNIT_PATH only for intel; see get_derecho_pfunit_version.
PFUNIT_CMAKE = os.path.join(
    REPO_ROOT, "ccs_config", "machines", "derecho", "intel_derecho.cmake"
)

COMPILER = "gnu"

# The two mpilib flavors of derecho's gnu stack. MPI_STACK is derecho's default
# MPILIB (the first entry of <MPILIBS>mpich,openmpi</MPILIBS>) and is what the
# container's MPICH stand-in replaces; SERIAL_STACK is the mpi-serial build.
# check_default_mpilib() asserts MPI_STACK is still that first entry: the
# per-stack module names below (netcdf-mpi vs netcdf, esmf-mpi vs esmf) are
# hand-written for mpich, so if derecho changes its default the right answer is
# to fail and have someone re-derive them, not to follow it silently.
MPI_STACK = "mpich"
SERIAL_STACK = "mpi-serial"

# Dockerfile ARG -> how to check it.
#   "direct":    ARG must equal the config gnu module version, in the "stack"
#                named by the entry.
#   "deviation": config gnu module version must equal the recorded value (the
#                ARG is a deliberate stand-in and is NOT compared).
#   "snapshot":  ARG must equal the recorded value (module absent from config).
# "strip_suffix" collapses derecho's "-debug" module variants onto the version
# they are a debug build of; "also_equal" asserts a second module, usually the
# other stack's twin, carries the same version as the primary one.
CHECKS = [
    # The gcc block carries no mpilib attribute, so it applies to both stacks;
    # MPI_STACK is named here only because every module check names one.
    {"arg": "GCC_VERSION", "mode": "direct", "module": "gcc", "stack": MPI_STACK},
    {
        # Like ESMF below: one ARG builds both of the container's netCDF-Cs
        # (the MPICH-linked one and the static serial one under
        # /usr/local/serial), so derecho's two must agree before either can be
        # compared to it.
        "arg": "NETCDF_C_VERSION",
        "mode": "direct",
        "module": "netcdf-mpi",
        "stack": MPI_STACK,
        "strip_suffix": "-debug",
        "also_equal": {"module": "netcdf", "stack": SERIAL_STACK},
    },
    {
        "arg": "PNETCDF_VERSION",
        "mode": "direct",
        "module": "parallel-netcdf",
        "stack": MPI_STACK,
        "strip_suffix": "-debug",
    },
    {
        # One ARG builds both of the container's ESMFs (the MPI one and the
        # mpiuni one), so derecho's two must agree before either can be
        # compared to it.
        "arg": "ESMF_VERSION",
        "mode": "direct",
        "module": "esmf-mpi",
        "stack": MPI_STACK,
        "strip_suffix": "-debug",
        "also_equal": {"module": "esmf", "stack": SERIAL_STACK},
    },
    {
        "arg": "MPICH_VERSION",
        "mode": "deviation",
        "module": "cray-mpich",
        "stack": MPI_STACK,
        "snap": ("deviation_guard", "cray_mpich"),
    },
    {
        "arg": "MPI_SERIAL_VERSION",
        "mode": "direct",
        "module": "mpi-serial",
        "stack": SERIAL_STACK,
    },
    {
        "arg": "PIO_VERSION",
        "mode": "direct",
        "module": "parallelio-serial",
        "stack": SERIAL_STACK,
        "strip_suffix": "-debug",
    },
    {"arg": "HDF5_VERSION", "mode": "snapshot", "snap": ("snapshot", "hdf5")},
    {
        "arg": "NETCDF_FORTRAN_VERSION",
        "mode": "snapshot",
        "snap": ("snapshot", "netcdf_fortran"),
    },
    {"arg": "PFUNIT_VERSION", "mode": "pfunit"},
]


def parse_dockerfile_args(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Dockerfile not found at {path}")
    args = {}
    pat = re.compile(r"^\s*ARG\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(\S+)")
    with open(path, encoding="utf-8") as f:
        for line in f:
            m = pat.match(line)
            if m:
                args[m.group(1)] = m.group(2).strip('"').strip("'")
    return args


def attr_applies(attr_value, target):
    """Does a <modules> attribute value select `target`?

    Mirrors CIME's EnvMachSpecific._match (cime/CIME/XML/env_mach_specific.py):
    an absent attribute matches anything, a leading "!" negates, and otherwise
    the value is an anchored regex -- so alternations ("gnu|nvhpc") work and a
    comma would be a literal, not a separator. No ccs_config machine file uses
    comma-separated values, so none is special-cased here.
    """
    if attr_value is None:
        return True  # unconstrained block, e.g. <modules mpilib="mpich">
    pattern = attr_value.strip()
    if pattern.startswith("!"):
        return re.match(pattern[1:] + "$", target) is None
    return re.match(pattern + "$", target) is not None


def check_default_mpilib(path):
    """Return (None) if derecho's default MPILIB is still MPI_STACK, else why not."""
    root = ET.parse(path).getroot()  # <machine MACH="derecho">
    el = root.find("MPILIBS")
    if el is None or not (el.text or "").strip():
        return f"no <MPILIBS> element in {path}; cannot confirm the default mpilib"
    first = el.text.split(",")[0].strip()
    if first != MPI_STACK:
        return (
            f"derecho's default mpilib is now {first!r}, not {MPI_STACK!r} "
            f"(<MPILIBS>{el.text.strip()}</MPILIBS>). Every MPI-stack check "
            "below reads modules named for mpich (netcdf-mpi, esmf-mpi, "
            "cray-mpich); re-derive them for the new default before trusting "
            "this check."
        )
    return None


def get_gnu_module_versions(path):
    """Return {mpilib: {module_name: {version, ...}}} for gnu load commands.

    Keyed by MPI_STACK / SERIAL_STACK, because since ccs_config_cesm1.0.88
    derecho's library modules live in <modules> blocks keyed by mpilib (and
    DEBUG) with no compiler attribute at all; only the gcc block is still
    compiler="gnu". Attributes other than compiler and mpilib are not filtered:
    DEBUG deliberately, so a module's debug and non-debug variants collapse via
    strip_suffix, and gpu_type harmlessly, since the only extra modules it
    admits (cuda) are ones no check reads.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"config_machines.xml not found at {path}. In CI make sure the "
            "'bin/git-fleximod update ccs_config' step ran; locally, run that "
            "first."
        )
    root = ET.parse(path).getroot()  # <machine MACH="derecho">
    module_system = root.find("module_system")
    if module_system is None:
        raise ValueError(f"no <module_system> element in {path}")
    stacks = {MPI_STACK: {}, SERIAL_STACK: {}}
    for modules in module_system.findall("modules"):
        # compiler="gnu" or a negation that spares gnu ("!intel") or no
        # compiler attribute at all; this still rejects "intel", "cray",
        # "nvhpc" and "!gnu", so the non-gnu twins never participate.
        if not attr_applies(modules.get("compiler"), COMPILER):
            continue
        applies_to = [s for s in stacks if attr_applies(modules.get("mpilib"), s)]
        if not applies_to:
            continue
        for cmd in modules.findall("command"):
            if cmd.get("name") != "load":
                continue
            text = (cmd.text or "").strip()
            if "/" not in text:
                continue  # unversioned load, e.g. "nco", "cmake"
            # Module names here are never path-like, so split name/version on
            # the last "/" (e.g. "esmf/8.6.0-debug" -> "esmf", "8.6.0-debug").
            name, version = text.rsplit("/", 1)
            for stack in applies_to:
                stacks[stack].setdefault(name, set()).add(version)
    return stacks


def resolve_config_version(stacks, module, stack, strip_suffix=None):
    """Return (version, None) or (None, reason) for a gnu module in `stack`.

    An absent module or conflicting versions are reported as a check finding
    (a returned reason -> per-component ❌), not raised: they are exactly the
    kind of derecho drift this check exists to flag.
    """
    where = (
        f'<modules> blocks of config_machines.xml applying to compiler="gnu", '
        f'mpilib="{stack}"'
    )
    versions = stacks[stack]
    if module not in versions:
        return None, f"module '{module}' not loaded by any of the {where}"
    vers = set(versions[module])
    if strip_suffix:
        vers = {
            v[: -len(strip_suffix)] if v.endswith(strip_suffix) else v
            for v in vers
        }
    if len(vers) != 1:
        return None, (
            f"module '{module}' has conflicting versions "
            f"{sorted(versions[module])} across the {where}; cannot pick one"
        )
    return vers.pop(), None


def get_derecho_pfunit_version(path):
    """Return (version, None) or (None, reason) for derecho's pFUnit.

    derecho has no pFUnit module, so there is nothing to read from
    config_machines.xml. CTSM's unit tests locate pFUnit through PFUNIT_PATH,
    whose value embeds the version:

        set(PFUNIT_PATH "$ENV{CESMDATAROOT}/tools/pFUnit/\
            pFUnit4.8.0_derecho_Intel2023.2.1_noMPI_noOpenMP")

    ccs_config sets that only in the intel macros -- there is no gnu
    PFUNIT_PATH for derecho -- so the container's gnu-built pFUnit is matched
    to the intel one by version. The compiler necessarily differs; the version
    is the part worth guarding, and this reads it live rather than recording a
    snapshot that could silently rot.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"intel_derecho.cmake not found at {path}. In CI make sure the "
            "'bin/git-fleximod update ccs_config' step ran; locally, run that "
            "first."
        )
    with open(path, encoding="utf-8") as f:
        text = f.read()
    matches = re.findall(
        r"""set\s*\(\s*PFUNIT_PATH\s+"[^"]*?pFUnit([0-9]+(?:\.[0-9]+)*)_""", text
    )
    if not matches:
        return None, (
            "no set(PFUNIT_PATH ... pFUnit<version>_ ...) found in "
            f"{os.path.basename(path)}; derecho may have moved or renamed its "
            "pFUnit install"
        )
    if len(set(matches)) != 1:
        return None, (
            f"conflicting pFUnit versions {sorted(set(matches))} in "
            f"{os.path.basename(path)}; cannot pick one"
        )
    return matches[0], None


def load_snapshot(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"snapshot file not found at {path}")
    # interpolation=None: version strings never contain "%", and this avoids a
    # latent InterpolationSyntaxError footgun if one ever did.
    parser = configparser.ConfigParser(interpolation=None)
    parser.read(path, encoding="utf-8")  # raises configparser.Error if malformed
    return parser


def snap_get(parser, section, key):
    # Raises configparser.NoSectionError / NoOptionError if absent.
    return parser.get(section, key).strip()


def main():
    args = parse_dockerfile_args(DOCKERFILE)
    config = get_gnu_module_versions(CONFIG_XML)
    mpilib_problem = check_default_mpilib(CONFIG_XML)
    snap = load_snapshot(SNAPSHOT)

    ok = True
    if mpilib_problem:
        print(f"\u274c default mpilib: {mpilib_problem}")
        ok = False
    else:
        print(
            f'\u2705 derecho\'s default mpilib is still "{MPI_STACK}"; the '
            "MPI-stack module names below apply"
        )
    for chk in CHECKS:
        arg = chk["arg"]
        mode = chk["mode"]
        if arg not in args:
            print(f"❌ {arg}: ARG not found in Dockerfile")
            ok = False
            continue
        arg_val = args[arg]

        if mode == "direct":
            stack = chk["stack"]
            cfg_val, reason = resolve_config_version(
                config, chk["module"], stack, chk.get("strip_suffix")
            )
            twin = chk.get("also_equal")
            twin_val, twin_reason = (None, None)
            if cfg_val is not None and twin:
                twin_val, twin_reason = resolve_config_version(
                    config, twin["module"], twin["stack"], chk.get("strip_suffix")
                )
            if cfg_val is None:
                print(f"❌ {arg}: {reason}")
                ok = False
            elif twin and twin_val is None:
                print(f"❌ {arg}: {twin_reason}")
                ok = False
            elif twin and twin_val != cfg_val:
                print(
                    f"❌ {arg}: derecho {twin['module']}/{twin_val} "
                    f"(mpilib=\"{twin['stack']}\") != {chk['module']}/{cfg_val} "
                    f'(mpilib="{stack}"). Derecho\'s two stacks have diverged, '
                    f"so the single {arg} that builds both of the container's "
                    "flavors can no longer match them; the Dockerfile needs a "
                    "second ARG."
                )
                ok = False
            elif arg_val == cfg_val:
                print(
                    f"✅ {arg}={arg_val} matches derecho {chk['module']}/{cfg_val} "
                    f'(gnu, mpilib="{stack}")'
                )
            else:
                print(
                    f"❌ {arg}={arg_val} != derecho {chk['module']}/{cfg_val} "
                    f'(config_machines.xml, gnu, mpilib="{stack}"). Update the '
                    "Dockerfile ARG, or the module changed on derecho."
                )
                ok = False

        elif mode == "deviation":
            live, reason = resolve_config_version(
                config, chk["module"], chk["stack"]
            )
            recorded = snap_get(snap, *chk["snap"])
            if live is None:
                print(f"❌ {arg} guard: {reason}")
                ok = False
            elif live == recorded:
                print(
                    f"✅ derecho {chk['module']} still {recorded} (recorded); "
                    f"{arg}={arg_val} is an intentional open-source stand-in, "
                    "not compared"
                )
            else:
                print(
                    f"❌ derecho {chk['module']} changed {recorded} (recorded) "
                    f"-> {live}. Re-evaluate the {arg}={arg_val} stand-in and "
                    "update [deviation_guard] in derecho-versions.ini."
                )
                ok = False

        elif mode == "pfunit":
            derecho_val, reason = get_derecho_pfunit_version(PFUNIT_CMAKE)
            if derecho_val is None:
                print(f"❌ {arg}: {reason}")
                ok = False
            elif arg_val == derecho_val:
                print(
                    f"✅ {arg}={arg_val} matches derecho pFUnit {derecho_val} "
                    "(from PFUNIT_PATH in intel_derecho.cmake; the container "
                    "builds it with gnu, a known compiler deviation)"
                )
            else:
                print(
                    f"❌ {arg}={arg_val} != derecho pFUnit {derecho_val} "
                    "(PFUNIT_PATH in intel_derecho.cmake). Update the "
                    "Dockerfile ARG, or derecho changed its pFUnit."
                )
                ok = False

        elif mode == "snapshot":
            recorded = snap_get(snap, *chk["snap"])
            if arg_val == recorded:
                print(
                    f"✅ {arg}={arg_val} matches recorded derecho {recorded} "
                    "(recorded from netcdf-mpi; not in config_machines.xml)"
                )
            else:
                print(
                    f"❌ {arg}={arg_val} != recorded derecho {recorded}. "
                    "Neither has a standalone module in config_machines.xml: "
                    "netCDF-Fortran is bundled in netcdf-mpi, HDF5 is the "
                    "hdf5-mpi module it pulls in. Verify on derecho "
                    "(module show netcdf-mpi | grep -i hdf5 / nf-config "
                    "--version) and update the Dockerfile ARG or [snapshot] "
                    "in derecho-versions.ini."
                )
                ok = False

        else:  # pragma: no cover - guards a typo in the CHECKS table above
            raise ValueError(f"unknown check mode {mode!r} for {arg}")

    # cray-libsci -> reference LAPACK/BLAS is a known deviation with no
    # Dockerfile version ARG (dnf installs lapack-devel/blas-devel unversioned),
    # so there is nothing to check for it.
    print(
        "note: cray-libsci -> reference lapack/blas is a known deviation with "
        "no version ARG; not checked."
    )
    print("OK: all version checks passed" if ok else "FAILED: version mismatch")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
