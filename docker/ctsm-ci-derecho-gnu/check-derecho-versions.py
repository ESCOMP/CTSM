#!/usr/bin/env python3
"""Check the ctsm-ci-derecho-gnu Dockerfile against derecho's gnu+mpi-serial stack.

The container replicates the software a **gnu + mpi-serial** standalone run on
derecho links. That is the subject of every check here: not derecho's gnu stack
considered in the abstract, and not its MPI flavor, which the container does not
have. Concretely, the modules read are the ones in the `<modules>` blocks of
ccs_config/machines/derecho/config_machines.xml whose `compiler` attribute
admits "gnu" and whose `mpilib` attribute admits "mpi-serial" -- so
netcdf/mpi-serial/parallelio-serial/esmf, never netcdf-mpi/parallelio/esmf-mpi,
and never cray-mpich or parallel-netcdf.

The Dockerfile's version ARGs are compared in four modes:

  direct    - the ARG must equal the derecho gnu+mpi-serial module version, read
              live from config_machines.xml (gcc, netcdf, esmf, mpi-serial,
              parallelio-serial).
  deviation - the container ships an intentional open-source stand-in that is
              NOT compared to derecho; instead the derecho module version (read
              live) must still equal a recorded value, so a derecho change trips
              the check and the stand-in gets re-examined (cray-libsci ->
              reference lapack/blas). A deviation need not have a version ARG at
              all: the lapack/blas stand-in is dnf-installed and unversioned.
  snapshot  - the derecho version has NO standalone entry in
              config_machines.xml -- netCDF-Fortran is bundled inside the netcdf
              module, HDF5 is the hdf5 module netcdf pulls in via depends_on --
              so the ARG is compared to a hand-recorded value (HDF5,
              netCDF-Fortran). Because that value cannot be read live, the ini
              also records which netcdf it was measured against; that IS read
              live, so a netcdf *version* bump that left the record behind
              fails instead of comparing stale against stale. That is the whole
              of what it covers: derecho can repoint an unchanged netcdf/<ver>
              at a different hdf5, and nothing here would notice. Closing that
              needs a live read of the module tree's depends_on and belongs to
              the scheduled drift check (NEXT_STEPS item 5), not to a check the
              PR gate runs with --skip-module-tree.
  pfunit    - CTSM does not use derecho's pFUnit module (there is one under
              ncarenv/25.10, built gnu). Its unit tests take pFUnit from the
              PFUNIT_PATH set by intel_derecho.cmake, which names a hand-built
              tree under CESMDATAROOT and embeds the version; that path is what
              is read live, so the module's version is deliberately not
              consulted.

Beyond the ARGs, every pinned module config_machines.xml names for this path is
checked for existence in derecho's module tree under /glade -- bar those listed
in UNCHECKABLE_PINS, whose "version" is not one. config_machines.xml
is CTSM's *claim* about derecho, never verified: if it pins a version derecho
has removed, a derecho run dies at `module load`, there is no run to compare
against, and comparing an ARG to the pinned number reports a match that means
nothing. That check resolves modulefiles only -- it never loads a module and
never invokes Lmod. Pass --skip-module-tree where /glade is unreachable.

Recorded values (snapshot + deviation guard) live in derecho-versions.ini.

Standard library only. Needs ccs_config populated (bin/git-fleximod update
ccs_config). A version mismatch prints a per-component report and exits 1; a
setup problem (missing/malformed inputs) raises an exception.
"""

import argparse
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

# The one path this check is about. Blocks carrying any other compiler or
# mpilib are not the container's subject and are filtered out; in particular
# mpilib="!mpi-serial" (the MPI flavor) must not leak in, which is what
# attr_applies()'s negation handling is for.
COMPILER = "gnu"
MPILIB = "mpi-serial"

# derecho's modulefiles live under two roots, and both are needed. The system
# root is NCAR's; the CSEG root is what cesmdev/1.0 adds -- its modulefile
# appends NCAR_MODULEROOT_CSEG to NCAR_VARS_MODULEROOT, and ncarenv then puts
# that root's Core and gcc/<ver> directories on MODULEPATH alongside the system
# root's. This is not a nicety: parallelio-serial/<ver>-debug, which
# config_machines.xml loads for DEBUG="TRUE", exists ONLY under the CSEG root, so
# a check that walks the system root alone fails on its very first run.
SYSTEM_MODULE_ROOT = "/glade/u/apps/derecho/modules"
CSEG_MODULE_ROOT = "/glade/u/apps/cesmdev/modules"

# Pinned loads whose "version" is not a version. conda/latest resolves --
# Core/conda/latest.lua is a real file -- but all the check could confirm is
# that the name exists, never that what it names is stable, so a pass would
# claim more than it knows. Matched on the full "name/version" string, so a
# future conda/<real version> is existence-checked like anything else.
UNCHECKABLE_PINS = {"conda/latest"}

# The [snapshot] values are read out of one derecho module: netcdf bundles
# netCDF-Fortran and pulls in the hdf5 module. Nothing about those versions can
# be read live, but the version of the module they came from can be -- recording
# it is what lets a stale snapshot fail instead of pass. If a snapshot ever comes
# from some other module, this grows a per-check field.
SNAPSHOT_SOURCE_MODULE = "netcdf"
SNAPSHOT_PROVENANCE_KEY = "measured_against_netcdf"

# Dockerfile ARG -> how to check it. See the module docstring for the modes.
# "strip_suffix" collapses derecho's "-debug" module variants onto the version
# they are a debug build of. "arg" is required except in deviation mode, where
# the stand-in may carry no version ARG.
CHECKS = [
    {"arg": "GCC_VERSION", "mode": "direct", "module": "gcc"},
    {"arg": "NETCDF_C_VERSION", "mode": "direct", "module": "netcdf"},
    {"arg": "ESMF_VERSION", "mode": "direct", "module": "esmf"},
    {"arg": "MPI_SERIAL_VERSION", "mode": "direct", "module": "mpi-serial"},
    {
        # DEBUG="TRUE" loads parallelio-serial/<ver>-debug, the debug build of
        # the same <ver>; one ARG covers both.
        "arg": "PIO_VERSION",
        "mode": "direct",
        "module": "parallelio-serial",
        "strip_suffix": "-debug",
    },
    {"arg": "HDF5_VERSION", "mode": "snapshot", "snap": ("snapshot", "hdf5")},
    {
        "arg": "NETCDF_FORTRAN_VERSION",
        "mode": "snapshot",
        "snap": ("snapshot", "netcdf_fortran"),
    },
    {"arg": "PFUNIT_VERSION", "mode": "pfunit"},
    {
        # No "arg": the stand-in is dnf's lapack-devel/blas-devel, installed
        # unversioned, so there is no Dockerfile ARG to compare or to report.
        "mode": "deviation",
        "module": "cray-libsci",
        "stand_in": "reference lapack/blas (dnf, unversioned)",
        "snap": ("deviation_guard", "cray_libsci"),
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


def require_ccs_config_file(path):
    """Raise unless a ccs_config input is present.

    Both callers want the same message, and the existence half has to be
    callable on its own: the pfunit check reads intel_derecho.cmake from inside
    the per-check loop, and raising there would discard the lines already
    gathered.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"{os.path.basename(path)} not found at {path}. In CI make sure "
            "the 'bin/git-fleximod update ccs_config' step ran; locally, run "
            "that first."
        )


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


def check_snapshot_provenance(config, snap):
    """Return None if [snapshot] still describes derecho's current module.

    The recorded values have no entry in config_machines.xml, so they cannot be
    read live and the check can only compare an ARG to a number someone typed.
    What CAN be read live is the version of the module they were measured from.
    If that has moved, the recorded numbers describe a bundle derecho no longer
    has -- and an equally stale ARG then compares equal to them and prints a
    pass, which is exactly how HDF5 and netCDF-Fortran stayed green across a
    netcdf bump.

    It catches a netcdf *version* bump and nothing else. derecho can repoint an
    unchanged netcdf/<ver> at a different hdf5: ncarenv/24.12 and ncarenv/25.10
    both carry netcdf/4.9.3, with depends_on("hdf5/1.12.3") and
    depends_on("hdf5/1.14.6") respectively. A ccs_config ncarenv bump that kept
    netcdf/4.9.3 would leave this guard green on a real HDF5 change, and hdf5
    is not named in config_machines.xml so the module-tree check never looks at
    it either. That gap is NEXT_STEPS item 5's to close.
    """
    recorded = snap_get(snap, "snapshot", SNAPSHOT_PROVENANCE_KEY)
    live, reason = resolve_config_version(config, SNAPSHOT_SOURCE_MODULE)
    if live is None:
        return f"cannot read derecho {SNAPSHOT_SOURCE_MODULE}: {reason}"
    if live != recorded:
        return (
            f"derecho {SNAPSHOT_SOURCE_MODULE} is {live}, but [snapshot] was "
            f"measured against {recorded}, so the recorded HDF5 and "
            "netCDF-Fortran versions describe a bundle derecho no longer has. "
            "Re-measure on derecho (recipe in derecho-versions.ini) and set "
            f"{SNAPSHOT_PROVENANCE_KEY} = {live}."
        )
    return None


def get_gnu_module_versions(path):
    """Return {module_name: {version, ...}} for the gnu + mpi-serial loads.

    One flat mapping, because one stack is in scope. Since ccs_config_cesm1.0.88
    derecho's library modules live in <modules> blocks keyed by mpilib (and
    DEBUG) with no compiler attribute at all; only the gcc block is still
    compiler="gnu". Attributes other than compiler and mpilib are not filtered:
    DEBUG deliberately, so a module's debug and non-debug variants collapse via
    strip_suffix, and gpu_type harmlessly, since the extra modules it admits
    (cuda) sit in mpilib="mpich" blocks this never reaches.

    Unversioned loads (nco, craype, cmake) are dropped. Only a pinned version
    can go stale in the way this check is about -- an unpinned load follows
    whatever derecho makes default, which is explicitly not drift.
    """
    require_ccs_config_file(path)
    root = ET.parse(path).getroot()  # <machine MACH="derecho">
    module_system = root.find("module_system")
    if module_system is None:
        raise ValueError(f"no <module_system> element in {path}")
    versions = {}
    for modules in module_system.findall("modules"):
        # compiler="gnu" or a negation that spares gnu ("!intel") or no compiler
        # attribute at all; this still rejects "intel", "cray", "nvhpc" and
        # "!gnu", so the non-gnu twins never participate. Likewise mpilib:
        # "mpi-serial" and "!mpich" are in, "mpich" and "!mpi-serial" are out.
        if not attr_applies(modules.get("compiler"), COMPILER):
            continue
        if not attr_applies(modules.get("mpilib"), MPILIB):
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
    """Return (version, None) or (None, reason) for a gnu + mpi-serial module.

    An absent module or conflicting versions are reported as a check finding
    (a returned reason -> per-component ❌), not raised: they are exactly the
    kind of derecho drift this check exists to flag.
    """
    where = (
        f'<modules> blocks of config_machines.xml applying to compiler="gnu", '
        f'mpilib="{MPILIB}"'
    )
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

    pFUnit is absent from config_machines.xml, so there is nothing to read
    there. derecho does carry a pfunit module under ncarenv/25.10, but CTSM
    does not use it: its unit tests locate pFUnit through PFUNIT_PATH, whose
    value embeds the version:

        set(PFUNIT_PATH "$ENV{CESMDATAROOT}/tools/pFUnit/\
            pFUnit4.8.0_derecho_Intel2023.2.1_noMPI_noOpenMP")

    ccs_config sets that only in the intel macros -- there is no gnu
    PFUNIT_PATH for derecho -- so the container's gnu-built pFUnit is matched
    to the intel one by version. The compiler necessarily differs; the version
    is the part worth guarding, and this reads it live rather than recording a
    snapshot that could silently rot.
    """
    require_ccs_config_file(path)  # also asserted up front; see run_checks
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
    # Safe after check_snapshot_keys(); raises configparser.NoSectionError /
    # NoOptionError if called without it.
    return parser.get(section, key).strip()


def check_snapshot_keys(parser, path):
    """Raise unless the ini carries every key the CHECKS table will ask for.

    Up front rather than lazily inside the loop. A missing key is a setup
    problem, not drift, so it must raise -- but configparser's own
    NoSectionError names neither the file nor the check that wanted the key,
    and raising mid-loop also throws away every line gathered so far.
    """
    required = [("snapshot", SNAPSHOT_PROVENANCE_KEY)]
    required += [chk["snap"] for chk in CHECKS if "snap" in chk]
    missing = [
        f"[{section}] {key}"
        for section, key in required
        if not parser.has_option(section, key)
    ]
    if missing:
        raise ValueError(
            f"{path} does not carry every hand-recorded value the checks "
            "need. Missing: " + ", ".join(missing) + ". See the HOW TO "
            "UPDATE recipe at the top of that file."
        )


def unsafe_module_name(name, version):
    """Return why a module name must not be joined onto a root, or None.

    Load commands come from XML and are split on the last "/", so a command
    naming an absolute path yields an absolute "name" -- and os.path.join drops
    everything before an absolute component, collapsing every searched path onto
    that one path and reporting a file outside derecho's tree as resolved. ".."
    does the same thing a level at a time. Neither can name a real derecho
    module, so both are reported as findings rather than searched.
    """
    for part in (name, version):
        if os.path.isabs(part):
            return "is an absolute path"
        if ".." in part.split("/"):
            return 'contains ".."'
    return None


def module_search_paths(name, version, env_version, gcc_version, roots):
    """Where a pinned gnu + mpi-serial module may legitimately live.

    These three tiers are exactly what lands on MODULEPATH for such a run, and
    deliberately nothing else: ncarenv's modulefile appends <root>/<env>/Core,
    gcc's appends <root>/<env>/gcc/<gccver>, and the environment tree holds
    ncarenv and cesmdev themselves (they are not under any ncarenv version).
    A module found only under <root>/<env>/cray-mpich/<ver>/gcc/<gccver>/ must
    NOT count as resolved -- that directory reaches MODULEPATH only once
    cray-mpich is loaded, which a serial run never does, so counting it would
    hide exactly the drift this check is for.
    """
    paths = [os.path.join(roots[0], "environment", name, version + ".lua")]
    for root in roots:
        paths.append(os.path.join(root, env_version, "Core", name, version + ".lua"))
    for root in roots:
        paths.append(
            os.path.join(
                root, env_version, "gcc", gcc_version, name, version + ".lua"
            )
        )
    return paths


def check_module_tree(config, roots):
    """Return (ok, lines): does every pinned module still exist on derecho?

    Every pinned module except those in UNCHECKABLE_PINS, which are reported as
    skipped rather than dropped.

    Resolves modulefiles under /glade; never loads a module and never invokes
    Lmod. Loading derecho's stack off derecho is impossible anyway (cray-mpich
    needs Cray PE) and is not needed to answer "does this version still exist".
    """
    lines = []
    ok = True
    if not os.path.isdir(roots[0]):
        lines.append(
            f"❌ module tree: cannot read {roots[0]}. This check resolves "
            "derecho's modulefiles directly, so it must run somewhere /glade "
            "is readable (the Cirrus runner, derecho or casper). Pass "
            "--skip-module-tree on a runner that has no /glade."
        )
        return False, lines  # per-module lines would all be the same noise
    if not os.path.isdir(roots[1]):
        lines.append(
            f"❌ module tree: cannot read {roots[1]}, the module root "
            "cesmdev/1.0 adds. Modules that live only there (e.g. the "
            "parallelio-serial -debug variants) cannot be found; the "
            "per-module results below are against the system root only."
        )
        ok = False
        roots = (roots[0],)

    # The ncarenv and gcc versions name two of the three tiers, so read them
    # from the parsed config rather than hardcoding them: when ccs_config bumps
    # either, the search must follow it in the same commit.
    env_version, env_reason = resolve_config_version(config, "ncarenv")
    gcc_version, gcc_reason = resolve_config_version(config, "gcc")
    if env_version is None or gcc_version is None:
        lines.append(
            "❌ module tree: cannot build the module search paths: "
            f"{env_reason or gcc_reason}. Both ncarenv and gcc must resolve, "
            "since their versions name the directories MODULEPATH gets."
        )
        return False, lines

    for name in sorted(config):
        for version in sorted(config[name]):
            if f"{name}/{version}" in UNCHECKABLE_PINS:
                # Said out loud, because a list that quietly got shorter is
                # the failure mode this whole check is against.
                lines.append(
                    f"⚠️  module {name}/{version}: not checked. "
                    f'"{version}" is not a version, so resolving the file '
                    "would say only that the name exists, not that what it "
                    "names is stable."
                )
                continue
            unsafe = unsafe_module_name(name, version)
            if unsafe:
                lines.append(
                    f"❌ module {name}/{version}: the module name {unsafe} and "
                    "so cannot be resolved relative to a module root. Fix the "
                    "<command name=\"load\"> in config_machines.xml."
                )
                ok = False
                continue
            paths = module_search_paths(
                name, version, env_version, gcc_version, roots
            )
            found = next((p for p in paths if os.path.isfile(p)), None)
            if found:
                # Name the file, not just "found": a module migrating between
                # the system and CSEG roots is then visible in the log rather
                # than invisible.
                lines.append(f"✅ module {name}/{version} resolves: {found}")
            else:
                searched = "\n     ".join(paths)
                lines.append(
                    f"❌ module {name}/{version} is named by "
                    "config_machines.xml but resolves nowhere on derecho's "
                    "MODULEPATH, so there is no derecho run to compare "
                    "against. Searched:\n     " + searched
                )
                ok = False
    return ok, lines


def run_checks(
    dockerfile=DOCKERFILE,
    config_xml=CONFIG_XML,
    pfunit_cmake=PFUNIT_CMAKE,
    snapshot=SNAPSHOT,
    module_roots=(SYSTEM_MODULE_ROOT, CSEG_MODULE_ROOT),
    skip_module_tree=False,
):
    """Return (ok, lines)."""
    args = parse_dockerfile_args(dockerfile)
    config = get_gnu_module_versions(config_xml)
    # The pfunit check reads intel_derecho.cmake from inside the loop below, so
    # assert it is there first: a missing input is a setup problem and must
    # raise before any check line exists, not halfway through the report. Only
    # a partially populated ccs_config reaches this, since config_machines.xml
    # comes from the same external and is read just above.
    require_ccs_config_file(pfunit_cmake)
    snap = load_snapshot(snapshot)
    check_snapshot_keys(snap, snapshot)
    snapshot_problem = check_snapshot_provenance(config, snap)

    ok = True
    lines = []
    if snapshot_problem:
        lines.append(f"❌ [snapshot] provenance: {snapshot_problem}")
        ok = False
    else:
        lines.append(
            "✅ [snapshot] was measured against the "
            f"{SNAPSHOT_SOURCE_MODULE} derecho has today"
        )
    for chk in CHECKS:
        arg = chk.get("arg")
        mode = chk["mode"]
        arg_val = None
        if arg is not None:
            if arg not in args:
                lines.append(f"❌ {arg}: ARG not found in Dockerfile")
                ok = False
                continue
            arg_val = args[arg]

        if mode == "direct":
            cfg_val, reason = resolve_config_version(
                config, chk["module"], chk.get("strip_suffix")
            )
            if cfg_val is None:
                lines.append(f"❌ {arg}: {reason}")
                ok = False
            elif arg_val == cfg_val:
                lines.append(
                    f"✅ {arg}={arg_val} matches derecho "
                    f"{chk['module']}/{cfg_val} (gnu, mpilib=\"{MPILIB}\")"
                )
            else:
                lines.append(
                    f"❌ {arg}={arg_val} != derecho {chk['module']}/{cfg_val} "
                    f'(config_machines.xml, gnu, mpilib="{MPILIB}"). Update '
                    "the Dockerfile ARG, or the module changed on derecho."
                )
                ok = False

        elif mode == "deviation":
            live, reason = resolve_config_version(config, chk["module"])
            recorded = snap_get(snap, *chk["snap"])
            stand_in = chk["stand_in"]
            # The stand-in may have no version ARG (dnf installs it
            # unversioned), so word both outcomes without naming one.
            which = f"{arg}={arg_val}" if arg else stand_in
            if live is None:
                lines.append(f"❌ {chk['module']} guard: {reason}")
                ok = False
            elif live == recorded:
                lines.append(
                    f"✅ derecho {chk['module']} still {recorded} (recorded); "
                    f"the container's {which} is an intentional open-source "
                    "stand-in, not compared"
                )
            else:
                lines.append(
                    f"❌ derecho {chk['module']} changed {recorded} (recorded) "
                    f"-> {live}. Re-evaluate the {which} stand-in and update "
                    "[deviation_guard] in derecho-versions.ini."
                )
                ok = False

        elif mode == "pfunit":
            derecho_val, reason = get_derecho_pfunit_version(pfunit_cmake)
            if derecho_val is None:
                lines.append(f"❌ {arg}: {reason}")
                ok = False
            elif arg_val == derecho_val:
                lines.append(
                    f"✅ {arg}={arg_val} matches derecho pFUnit {derecho_val} "
                    "(from PFUNIT_PATH in intel_derecho.cmake; the container "
                    "builds it with gnu, a known compiler deviation)"
                )
            else:
                lines.append(
                    f"❌ {arg}={arg_val} != derecho pFUnit {derecho_val} "
                    "(PFUNIT_PATH in intel_derecho.cmake). Update the "
                    "Dockerfile ARG, or derecho changed its pFUnit."
                )
                ok = False

        elif mode == "snapshot":
            recorded = snap_get(snap, *chk["snap"])
            if snapshot_problem:
                lines.append(
                    f"❌ {arg}={arg_val}: not compared. The recorded derecho "
                    f"value ({recorded}) came from a superseded "
                    f"{SNAPSHOT_SOURCE_MODULE}; see the [snapshot] provenance "
                    "failure above."
                )
                ok = False
            elif arg_val == recorded:
                lines.append(
                    f"✅ {arg}={arg_val} matches recorded derecho {recorded} "
                    f"(recorded from {SNAPSHOT_SOURCE_MODULE}; not in "
                    "config_machines.xml)"
                )
            else:
                lines.append(
                    f"❌ {arg}={arg_val} != recorded derecho {recorded}. "
                    "Neither has a standalone module in config_machines.xml: "
                    "netCDF-Fortran is bundled in netcdf, HDF5 is the hdf5 "
                    "module it pulls in. Verify on derecho (module show "
                    "netcdf | grep -i hdf5 / nf-config --version) and update "
                    "the Dockerfile ARG or [snapshot] in derecho-versions.ini."
                )
                ok = False

        else:  # pragma: no cover - guards a typo in the CHECKS table above
            raise ValueError(f"unknown check mode {mode!r}")

    if skip_module_tree:
        lines.append(
            "⚠️  derecho module tree existence: NOT CHECKED "
            "(--skip-module-tree). Nothing above confirms that the versions "
            "config_machines.xml pins still exist on derecho."
        )
    else:
        tree_ok, tree_lines = check_module_tree(config, module_roots)
        lines.extend(tree_lines)
        ok = ok and tree_ok

    lines.append("OK: all version checks passed" if ok else "FAILED: version mismatch")
    return ok, lines


def main(argv=None):
    """Parse args, run the checks, print the lines, return ok."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--skip-module-tree",
        action="store_true",
        help="Skip resolving derecho's modulefiles under /glade. This exists "
        "for the hosted PR-gate runner, which has no /glade and can only "
        "compare in-repo files; everywhere else the check runs.",
    )
    opts = parser.parse_args(argv)
    ok, lines = run_checks(skip_module_tree=opts.skip_module_tree)
    for line in lines:
        print(line)
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
