"""Implementation of the CIME LILACSMOKE (LILAC smoke) test.

This is a CTSM-specific test. It tests the building and running of CTSM via LILAC. Compset
is ignored, but grid is important. Also, it's important that this test use the nuopc
driver, both for the sake of the build and for extracting some runtime settings. This test
should also use the lilac testmod (or a testmod that derives from it) in order to
establish the user_nl_ctsm file correctly.

Important directories under CASEROOT are:
- lilac_build: this contains the build and the runtime inputs for the lilac run
- lilac_atm_driver: this contains the build of the test driver as well as the run
  directory in which the test is actually run

Note that namelists for this test are generated in the build phase; they are NOT
regenerated when the test is submitted / run. This means that, if you have made any
changes that will impact namelists, you will need to rebuild this test.

Note that this test is tied to a specific resolution (10x15) and has a hard-coded domain
file for this resolution: see the setting of lnd_domain_file below.

"""

import glob
import os
import shutil

from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.utils import run_cmd, run_cmd_no_fail, symlink_force, new_lid, safe_copy
from CIME.status import append_testlog
from CIME.build import post_build
from CIME.test_status import (
    NAMELIST_PHASE,
    GENERATE_PHASE,
    BASELINE_PHASE,
    TEST_PASS_STATUS,
    TEST_FAIL_STATUS,
)
from CIME.XML.standard_module_setup import *

logger = logging.getLogger(__name__)

_LILAC_RUNTIME_FILES = ["lnd_in", "lnd_modelio.nml", "drv_flds_in", "lilac_in"]


class LILACSMOKE(SystemTestsCommon):
    def __init__(self, case):
        SystemTestsCommon.__init__(self, case)

    def build_phase(self, sharedlib_only=False, model_only=False):
        if not sharedlib_only:
            lndroot = self._case.get_value("COMP_ROOT_DIR_LND")
            exeroot = self._case.get_value("EXEROOT")
            build_dir = self._lilac_build_dir()
            script_path = os.path.abspath(os.path.join(lndroot, "lilac", "build_ctsm"))

            # We only run the initial build command if the build_dir doesn't exist
            # yet. This is to support rebuilding the test case. (The first time through,
            # the build_dir won't exist yet; subsequent times, it will already exist, so
            # we skip to the rebuild command.)
            if not os.path.isdir(build_dir):
                machine = self._case.get_value("MACH")
                compiler = self._case.get_value("COMPILER")
                debug = self._case.get_value("DEBUG")
                # It would be possible to do this testing via the python interface rather
                # than through a separate subprocess. However, we do it through a
                # subprocess in order to test the full build_ctsm script, including
                # command-line parsing.
                cmd = "{script_path} {build_dir} --machine {machine} --compiler {compiler}".format(
                    script_path=script_path, build_dir=build_dir, machine=machine, compiler=compiler
                )
                # It isn't straightforward to determine if pnetcdf is available on a
                # machine. To keep things simple, always run without pnetcdf.
                cmd += " --no-pnetcdf"
                if debug:
                    cmd += " --build-debug"
                self._run_build_cmd(cmd, exeroot, "build_ctsm.bldlog")

            # We call the build script with --rebuild even for an initial build. This is
            # so we make sure to test the code path for --rebuild. (This is also needed if
            # the user rebuilds the test case, in which case this will be the only command
            # run, since the build_dir will already exist.)
            cmd = "{script_path} {build_dir} --rebuild".format(
                script_path=script_path, build_dir=build_dir
            )
            self._run_build_cmd(cmd, exeroot, "rebuild_ctsm.bldlog")

            self._build_atm_driver()

            self._create_link_to_atm_driver()

            self._create_runtime_inputs()

            self._setup_atm_driver_rundir()

            self._cmpgen_namelists()

            # Setting logs=[] implies that we don't bother gzipping any of the build log
            # files; that seems fine for these purposes (and it keeps the above code
            # simpler).
            post_build(self._case, logs=[], build_complete=True)

    def _build_atm_driver(self):
        caseroot = self._case.get_value("CASEROOT")
        lndroot = self._case.get_value("COMP_ROOT_DIR_LND")
        blddir = os.path.join(caseroot, "lilac_atm_driver", "bld")

        if not os.path.exists(blddir):
            os.makedirs(blddir)
        symlink_force(
            os.path.join(lndroot, "lilac", "atm_driver", "Makefile"),
            os.path.join(blddir, "Makefile"),
        )
        symlink_force(
            os.path.join(lndroot, "lilac", "atm_driver", "atm_driver.F90"),
            os.path.join(blddir, "atm_driver.F90"),
        )
        symlink_force(
            os.path.join(self._lilac_build_dir(), "case", "Macros.make"),
            os.path.join(blddir, "Macros.make"),
        )

        makevars = "COMPILER={compiler} DEBUG={debug} CTSM_MKFILE={ctsm_mkfile}".format(
            compiler=self._case.get_value("COMPILER"),
            debug=str(self._case.get_value("DEBUG")).upper(),
            ctsm_mkfile=os.path.join(self._lilac_build_dir(), "ctsm.mk"),
        )
        makecmd = "make {makevars} atm_driver".format(makevars=makevars)

        # Normally the user will source either ctsm_build_environment.sh or
        # ctsm_build_environment.csh before building the atmosphere model. In the context
        # of this test case, case.load_env does the equivalent.
        self._case.load_env()

        self._run_build_cmd(makecmd, blddir, "atm_driver.bldlog")

    def _create_link_to_atm_driver(self):
        caseroot = self._case.get_value("CASEROOT")
        run_exe = (self._case.get_value("run_exe")).strip()

        # Make a symlink to the atm_driver executable so that the case's run command finds
        # it in the expected location
        symlink_force(os.path.join(caseroot, "lilac_atm_driver", "bld", "atm_driver.exe"), run_exe)

    def _create_runtime_inputs(self):
        caseroot = self._case.get_value("CASEROOT")
        runtime_inputs = self._runtime_inputs_dir()

        # NOTE: *** this test is currently tied to this single 10x15 grid resolution ***
        lnd_grid = self._case.get_value("LND_GRID")
        expect(
            lnd_grid == "10x15", "this test is currently tied to this single 10x15 grid resolution"
        )
        lnd_domain_file = os.path.join(
            self._case.get_value("DIN_LOC_ROOT"),
            "share",
            "domains",
            "domain.lnd.fv10x15_gx3v7.180321.nc",
        )

        # Cheat a bit here: Get the fsurdat file from the already-generated lnd_in file in
        # the host test case - i.e., from the standard cime-based preview_namelists. But
        # this isn't really a morally-objectionable cheat, because in the real workflow,
        # we expect the user to identify fsurdat manually; in this testing situation, we
        # need to come up with some way to replace this manual identification, so cheating
        # feels acceptable.
        self._case.create_namelists(component="lnd")
        fsurdat = self._extract_var_from_namelist(
            nl_filename=os.path.join(caseroot, "CaseDocs", "lnd_in"), varname="fsurdat"
        )

        self._fill_in_variables_in_file(
            filepath=os.path.join(runtime_inputs, "ctsm.cfg"),
            replacements={"lnd_domain_file": lnd_domain_file, "fsurdat": fsurdat},
        )

        # The user_nl_ctsm in the case directory is set up based on the standard testmods
        # mechanism. We use that one in place of the standard user_nl_ctsm, since the one
        # in the case directory may contain test-specific modifications.
        shutil.copyfile(
            src=os.path.join(caseroot, "user_nl_ctsm"),
            dst=os.path.join(runtime_inputs, "user_nl_ctsm"),
        )

        script_to_run = os.path.join(runtime_inputs, "make_runtime_inputs")
        self._run_build_cmd(
            "{} --rundir {}".format(script_to_run, runtime_inputs),
            runtime_inputs,
            "make_runtime_inputs.log",
        )

        # In lilac_in, we intentionally use the land mesh file for both atm_mesh_filename
        # and lnd_mesh_filename
        lnd_mesh = self._case.get_value("LND_DOMAIN_MESH")
        casename = self._case.get_value("CASE")
        self._fill_in_variables_in_file(
            filepath=os.path.join(runtime_inputs, "lilac_in"),
            replacements={
                "caseid": casename,
                "atm_mesh_filename": lnd_mesh,
                "lnd_mesh_filename": lnd_mesh,
                "lilac_histfreq_option": "ndays",
            },
            placeholders={"caseid": "ctsm_lilac", "lilac_histfreq_option": "never"},
        )

        # We run download_input_data partly because it may be needed and partly to test
        # this script.
        script_to_run = os.path.join(runtime_inputs, "download_input_data")
        self._run_build_cmd(
            "{} --rundir {}".format(script_to_run, runtime_inputs),
            runtime_inputs,
            "download_input_data.log",
        )

    def _setup_atm_driver_rundir(self):
        """Set up the directory from which we will actually do the run"""
        lndroot = self._case.get_value("COMP_ROOT_DIR_LND")
        rundir = self._atm_driver_rundir()

        if not os.path.exists(rundir):
            os.makedirs(rundir)
            shutil.copyfile(
                src=os.path.join(lndroot, "lilac", "atm_driver", "atm_driver_in"),
                dst=os.path.join(rundir, "atm_driver_in"),
            )

        # As elsewhere: assume the land variables also apply to the atmosphere
        lnd_mesh = self._case.get_value("LND_DOMAIN_MESH")
        lnd_nx = self._case.get_value("LND_NX")
        lnd_ny = self._case.get_value("LND_NY")
        expect(
            self._case.get_value("STOP_OPTION") == "ndays",
            "LILAC testing currently assumes STOP_OPTION of ndays, not {}".format(
                self._case.get_value("STOP_OPTION")
            ),
        )
        stop_n = self._case.get_value("STOP_N")
        casename = self._case.get_value("CASE")
        self._fill_in_variables_in_file(
            filepath=os.path.join(rundir, "atm_driver_in"),
            replacements={
                "caseid": casename,
                "atm_mesh_file": lnd_mesh,
                "atm_global_nx": str(lnd_nx),
                "atm_global_ny": str(lnd_ny),
                "atm_stop_day": str(stop_n + 1),
                "atm_ndays_all_segs": str(stop_n),
            },
        )

        for file_to_link in _LILAC_RUNTIME_FILES:
            symlink_force(
                os.path.join(self._runtime_inputs_dir(), file_to_link),
                os.path.join(rundir, file_to_link),
            )

        init_generated_files_dir = os.path.join(rundir, "init_generated_files")
        if not os.path.exists(init_generated_files_dir):
            os.mkdir(init_generated_files_dir)

    def _cmpgen_namelists(self):
        """Redoes the namelist comparison & generation with appropriate namelists

        The standard namelist comparison & generation is done with the CaseDocs directory
        from the test case. That isn't appropriate here, because those namelists aren't
        actually used in this test. Instead, we want to compare & generate the namelists
        used by the atm_driver-lilac-ctsm execution. Here, we do some file copies and then
        re-call the namelist comparison & generation script in order to accomplish
        this. This will overwrite the namelists generated earlier in the test, and will
        also replace the results of the NLCOMP phase.

        Note that we expect a failure in the NLCOMP phase that is run earlier in the test,
        because that one will have compared the test's standard CaseDocs with the files
        generated from here - and those two sets of namelists can be quite different.
        """
        caseroot = self._case.get_value("CASEROOT")
        casedocs = os.path.join(caseroot, "CaseDocs")
        if os.path.exists(casedocs):
            shutil.rmtree(casedocs)
        os.makedirs(casedocs)

        # case_cmpgen_namelists uses the existence of drv_in to decide whether namelists
        # need to be regenerated. We do NOT want it to regenerate namelists, so we give it
        # the file it wants.
        with open(os.path.join(casedocs, "drv_in"), "a") as drv_in:
            pass

        for onefile in _LILAC_RUNTIME_FILES + ["atm_driver_in"]:
            safe_copy(
                os.path.join(self._atm_driver_rundir(), onefile), os.path.join(casedocs, onefile)
            )

        success = self._case.case_cmpgen_namelists()
        # The setting of the NLCOMP phase status in case_cmpgen_namelists doesn't work
        # here (probably because the test object has a saved version of the test status
        # and so, when it goes to write the status of the build phase, it ends up
        # overwriting whatever was set by case_cmpgen_namelists). So we need to set it
        # here.
        with self._test_status:
            self._test_status.set_status(
                NAMELIST_PHASE,
                TEST_PASS_STATUS if success else TEST_FAIL_STATUS,
                comments="(used lilac namelists)",
            )

    def _extract_var_from_namelist(self, nl_filename, varname):
        """Tries to find a variable named varname in the given file; returns its value

        If not found, aborts
        """
        with open(nl_filename) as nl_file:
            for line in nl_file:
                match = re.search(r'^ *{} *= *[\'"]([^\'"]+)'.format(varname), line)
                if match:
                    return match.group(1)
        expect(False, "{} not found in {}".format(varname, nl_filename))

    def _fill_in_variables_in_file(self, filepath, replacements, placeholders=None):
        """For the given file, make the given replacements

        replacements should be a dictionary mapping variable names to their values

        If placeholders is given, it should be a dictionary mapping some subset of
        variable names to their placeholders. Anything not given here uses a placeholder
        of 'FILL_THIS_IN'.
        """
        if placeholders is None:
            placeholders = {}

        orig_filepath = "{}.orig".format(filepath)
        if not os.path.exists(orig_filepath):
            shutil.copyfile(src=filepath, dst=orig_filepath)
        os.remove(filepath)

        counts = dict.fromkeys(replacements, 0)
        with open(orig_filepath) as orig_file:
            with open(filepath, "w") as new_file:
                for orig_line in orig_file:
                    line = orig_line
                    for varname in replacements:
                        if varname in placeholders:
                            this_placeholder = placeholders[varname]
                        else:
                            this_placeholder = "FILL_THIS_IN"
                        line, replacement_done = self._fill_in_variable(
                            line=line,
                            varname=varname,
                            value=replacements[varname],
                            placeholder=this_placeholder,
                        )
                        if replacement_done:
                            counts[varname] += 1
                            break
                    new_file.write(line)

        for varname in counts:
            expect(
                counts[varname] > 0,
                "Did not find any instances of <{}> to replace in {}".format(varname, filepath),
            )

    def _fill_in_variable(self, line, varname, value, placeholder):
        """Fill in a placeholder variable in a config or namelist file

        Returns a tuple: (newline, replacement_done)
        - newline is the line with the given placeholder replaced with the given value if this
          line is for varname; otherwise returns line unchanged
        - replacement_done is True if the replacement was done, otherwise False
        """
        if re.search(r"^ *{} *=".format(varname), line):
            expect(
                placeholder in line,
                "Placeholder to replace ({}) not found in <{}>".format(placeholder, line.strip()),
            )
            newline = line.replace(placeholder, value)
            replacement_done = True
        else:
            newline = line
            replacement_done = False
        return (newline, replacement_done)

    def _lilac_build_dir(self):
        return os.path.join(self._case.get_value("CASEROOT"), "lilac_build")

    def _runtime_inputs_dir(self):
        return os.path.join(self._lilac_build_dir(), "runtime_inputs")

    def _atm_driver_rundir(self):
        return os.path.join(self._case.get_value("CASEROOT"), "lilac_atm_driver", "run")

    @staticmethod
    def _run_build_cmd(cmd, exeroot, logfile):
        """
        Runs the given build command, with output to the given logfile

        Args:
        cmd: str (command to run)
        exeroot: str (path to exeroot)
        logfile: str (path to logfile)
        """
        append_testlog(cmd)
        run_cmd_no_fail(cmd, arg_stdout=logfile, combine_output=True, from_dir=exeroot)
        with open(os.path.join(exeroot, logfile)) as lf:
            append_testlog(lf.read())

    def run_phase(self):
        # This mimics a bit of what's done in the typical case.run. Note that
        # case.get_mpirun_cmd creates a command that runs the executable given by
        # case.run_exe. So it's important that (elsewhere in this test script) we create a
        # link pointing from that to the atm_driver.exe executable.
        #
        # 2024/5/28 slevis: We added the load_env here to replace the
        # behavior of the PBS -V directive that was removed from
        # /ccs_config/machines/config_batch.xml
        self._case.load_env(reset=True)
        lid = new_lid()
        os.environ["OMP_NUM_THREADS"] = str(self._case.thread_count)
        cmd = self._case.get_mpirun_cmd(allow_unresolved_envvars=False)
        run_cmd_no_fail(cmd, from_dir=self._atm_driver_rundir())

        self._link_to_output_files()

    def _link_to_output_files(self):
        """Make links to the output files so that they appear in the directory expected by the test case

        Note: We do the run from a different directory in order to ensure that the run
        isn't using any of the files that are staged by the test case in the standard run
        directory. But then we need to create these links afterwards for the sake of
        baseline generation / comparison.
        """
        casename = self._case.get_value("CASE")
        rundir = self._case.get_value("RUNDIR")
        pattern = "{}*.nc".format(casename)

        # First remove any old files from the run directory
        old_files = glob.glob(os.path.join(rundir, pattern))
        for one_file in old_files:
            os.remove(one_file)

        # Now link to new files
        output_files = glob.glob(os.path.join(self._atm_driver_rundir(), pattern))
        for one_file in output_files:
            file_basename = os.path.basename(one_file)
            symlink_force(one_file, os.path.join(rundir, file_basename))
