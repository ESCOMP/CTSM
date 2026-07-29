.. include:: ../substitutions.rst

.. _modifying-surface-datasets:

===========================
 Modifying surface datasets
===========================

fsurdat_modifier is a tool that modifies fsurdat files. It reads a surface dataset (fsurdat file) and outputs a modified copy of the same file. Current applications are limited to the simplest CTSM SP mode, so bgc, fire, urban, vic, lake, transient, and crop-related variables in the fsurdat file remain unchanged.

This tool differs from `modify_singlept_site_neon` in that the latter specifically modifies soil properties of single-point surface datasets. It also differs from the `subset_data` tool in that the latter subsets fsurdat files to regional or single-point domains; `subset_data` does include some "modify" functionality when subsetting, but such functionality is more prescriptive in `subset_data`. In particular:

===================================  ======================================
`fsurdat_modifier` options           `subset_data` options                   
===================================  ======================================
std_elev (user sets STD_ELEV value)  uniform-snowpack (sets STD_ELEV to 20)
max_sat_area (user sets FMAX value)  cap-saturation (sets FMAX to zero)    
===================================  ======================================

Files involved
--------------

::

  python/ctsm/modify_input_files/fsurdat_modifier.py
  python/ctsm/modify_input_files/modify_fsurdat.py
  tools/modify_input_files/fsurdat_modifier
  tools/modify_input_files/modify_fsurdat_template.cfg

Instructions
------------

1) Activate conda (however you do this on your system and if not already active), run py_env_create (if necessary), and activate ctsm_pylib:

::

   ./py_env_create  # once per machine, unless needing to update the ctsm_pylib environment
   conda activate ctsm_pylib  # every time you come to this step

(Use "deactivate" to reverse the latter.)

2) Copy, then modify the configure file named `modify_fsurdat_template.cfg`, which contains all the arguments needed by the script.

3) Run the script `./fsurdat_modifier` pointing to the copied/modified `.cfg` file, e.g.

::

  ./fsurdat_modifier modify_users_copy.cfg

See `modify_fsurdat_template.cfg` for required and optional settings.

4) Use the `--verbose option` to see progress output on your screen.

