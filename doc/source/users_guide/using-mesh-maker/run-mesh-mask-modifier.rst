.. include:: ../substitutions.rst

.. _run-mesh-mask-modifier:

========================================
 Modifying the mask of an ESMF mesh file
========================================

The mesh_mask_modifier tool modifies the mask in ESMF mesh files. It reads a mesh file and outputs a modified copy of the same file.

Files involved
--------------

::

   python/ctsm/modify_input_files/mesh_mask_modifier.py
   python/ctsm/modify_input_files/modify_mesh_mask.py
   tools/modify_input_files/mesh_mask_modifier
   tools/modify_input_files/modify_mesh_template.cfg

Instructions
------------

1) Activate conda (however you do this on your system and if not already active), run py_env_create (if necessary), and activate ctsm_pylib:

::

   ./py_env_create  # once per machine, unless needing to update the ctsm_pylib environment
   conda activate ctsm_pylib  # every time you come to this step

(Use "deactivate" to reverse the latter.)

2) Copy, then modify the configure file named ``modify_mesh_template.cfg``, which contains all the arguments needed by the script.

3) Run the script ``mesh_mask_modifier`` pointing to the copied/modified ``.cfg`` file, e.g. ``./mesh_mask_modifier modify_users_copy.cfg``


================================================
Example: F-Case, modify the continental geometry
================================================

User wants to make the Indian Ocean into grassland. In a netcdf file, they specify their own land mask on the CESM 1-degree grid, as well as the area to be specified as grassland. This has been obtained by modifying the default land fraction of CESM. The file contains two arrays:

- ``landmask`` = the new landmask
- ``mod_lnd_props`` = set to 1 where the new land surface has been specified (i.e., where grassland needs to be specified) and zero elsewhere

This use-case requires modification to both the fsurdat file and the mesh file. To modify the former, use the ``modify_fsurdat`` tool (see section :numref:`modifying-surface-datasets`). Here are the steps to modify the mesh file:

In your copy of the CTSM (say, ``~<user>/ctsm``), go to the appropriate tool:

::

   cd tools/modify_input_files
   cp modify_mesh_template.cfg modify_fill_indianocean.cfg

Enter the following (or similar) selections in ``modify_fill_indianocean.cfg``:

::

   mesh_mask_in = /glade/campaign/cesm/cesmdata/cseg/inputdata/share/meshes/fv0.9x1.25_141008_polemod_ESMFmesh.nc
   mesh_mask_out = fv0.9x1.25_141008_polemod_ESMFmesh_modified.nc
   landmask_file = .../path_to_your_copy_of/fill_indianocean.nc

Run the tool

::

   ./mesh_mask_modifier modify_fill_indianocean.cfg

A modified mesh file should appear in the directory where you ran. Point to this file in your case's ``env_run.xml`` in the line that sets ``MASK_MESH``.

.. note:: If for some reason this fails, hardwire the ocean domain mesh file name in ``~<user>/ctsm/ccs_config/component_grids_nuopc.xml`` before starting your CTSM or CESM simulation. In the specific example shown here, hardwire the mesh file name for domain name ``gx1v7``. (This note has not been updated since around 2021.)
