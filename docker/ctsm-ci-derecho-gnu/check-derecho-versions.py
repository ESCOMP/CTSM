#!/usr/bin/env python3
"""Check that the ctsm-ci-derecho-gnu Dockerfile versions match derecho's gnu stack.

The container replicates derecho's gnu software stack (declared in
ccs_config/machines/derecho/config_machines.xml). This check compares the
Dockerfile's version ARGs to that config, in three modes:

  direct    - the ARG must equal the derecho gnu module version, read live from
              config_machines.xml (gcc, netcdf-mpi, parallel-netcdf, esmf).
  deviation - the ARG is an intentional open-source stand-in that is NOT
              compared to derecho; instead the derecho module version (read
              live) must still equal a recorded value, so a derecho change
              trips the check (MPICH_VERSION <-> cray-mpich).
  snapshot  - the derecho version is bundled inside another Spack module and has
              NO standalone entry in config_machines.xml, so the ARG is compared
              to a hand-recorded value (HDF5, netCDF-Fortran).

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

COMPILER = "gnu"

# Dockerfile ARG -> how to check it.
#   "direct":    ARG must equal the config gnu module version.
#   "deviation": config gnu module version must equal the recorded value (the
#                ARG is a deliberate stand-in and is NOT compared).
#   "snapshot":  ARG must equal the recorded value (module absent from config).
CHECKS = [
    {"arg": "GCC_VERSION", "mode": "direct", "module": "gcc"},
    {"arg": "NETCDF_C_VERSION", "mode": "direct", "module": "netcdf-mpi"},
    {"arg": "PNETCDF_VERSION", "mode": "direct", "module": "parallel-netcdf"},
    {
        "arg": "ESMF_VERSION",
        "mode": "direct",
        "module": "esmf",
        "strip_suffix": "-debug",
    },
    {
        "arg": "MPICH_VERSION",
        "mode": "deviation",
        "module": "cray-mpich",
        "snap": ("deviation_guard", "cray_mpich"),
    },
    {"arg": "HDF5_VERSION", "mode": "snapshot", "snap": ("snapshot", "hdf5")},
    {
        "arg": "NETCDF_FORTRAN_VERSION",
        "mode": "snapshot",
        "snap": ("snapshot", "netcdf_fortran"),
    },
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


def get_gnu_module_versions(path):
    """Return {module_name: {version, ...}} for compiler="gnu" load commands."""
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
    versions = {}
    for modules in module_system.findall("modules"):
        # Exact "gnu" match: excludes compiler="!gnu" (and no-compiler) blocks,
        # so the non-gnu twins (cray-mpich/8.1.29, netcdf-mpi/4.9.3,
        # parallel-netcdf/1.14.0) never participate.
        if modules.get("compiler") != COMPILER:
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
            versions.setdefault(name, set()).add(version)
    return versions


def resolve_config_version(versions, module, strip_suffix=None):
    """Return (version, None) or (None, reason) for a gnu module.

    An absent module or conflicting versions are reported as a check finding
    (a returned reason -> per-component ❌), not raised: they are exactly the
    kind of derecho drift this check exists to flag.
    """
    if module not in versions:
        return None, (
            f"module '{module}' not found in any compiler=\"gnu\" <modules> "
            "block of config_machines.xml"
        )
    vers = set(versions[module])
    if strip_suffix:
        vers = {
            v[: -len(strip_suffix)] if v.endswith(strip_suffix) else v
            for v in vers
        }
    if len(vers) != 1:
        return None, (
            f"module '{module}' has conflicting gnu versions "
            f"{sorted(versions[module])}; cannot pick one"
        )
    return vers.pop(), None


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
    snap = load_snapshot(SNAPSHOT)

    ok = True
    for chk in CHECKS:
        arg = chk["arg"]
        mode = chk["mode"]
        if arg not in args:
            print(f"❌ {arg}: ARG not found in Dockerfile")
            ok = False
            continue
        arg_val = args[arg]

        if mode == "direct":
            cfg_val, reason = resolve_config_version(
                config, chk["module"], chk.get("strip_suffix")
            )
            if cfg_val is None:
                print(f"❌ {arg}: {reason}")
                ok = False
            elif arg_val == cfg_val:
                print(
                    f"✅ {arg}={arg_val} matches derecho {chk['module']}/{cfg_val}"
                )
            else:
                print(
                    f"❌ {arg}={arg_val} != derecho {chk['module']}/{cfg_val} "
                    "(config_machines.xml). Update the Dockerfile ARG, or the "
                    "module changed on derecho."
                )
                ok = False

        elif mode == "deviation":
            live, reason = resolve_config_version(config, chk["module"])
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

        elif mode == "snapshot":
            recorded = snap_get(snap, *chk["snap"])
            if arg_val == recorded:
                print(
                    f"✅ {arg}={arg_val} matches recorded derecho {recorded} "
                    "(bundled in netcdf-mpi; not in config_machines.xml)"
                )
            else:
                print(
                    f"❌ {arg}={arg_val} != recorded derecho {recorded}. HDF5 "
                    "and netCDF-Fortran are bundled in derecho's netcdf-mpi "
                    "module (no standalone module in config). Verify on derecho "
                    "(module show netcdf-mpi / nf-config --version) and update "
                    "the Dockerfile ARG or [snapshot] in derecho-versions.ini."
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
