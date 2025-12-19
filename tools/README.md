# CTSM Tools for Preprocessing of Input Datasets or Postprocessing of History Output
#### $CTSMROOT/tools/README

CTSM tools for analysis of CTSM history files -- or for creation or
modification of CTSM input files.
	
I.  General directory structure:

    `$CTSMROOT/tools`
	mksurfdata_esmf -- Create surface datasets.

        crop_calendars --- Regrid and process GGCMI sowing and harvest date files for use in CTSM.

        mkmapgrids ------- Create regular lat/lon SCRIP grid files

        site_and_regional  Scripts for handling input datasets for site and regional cases.
                           These scripts both help with creation of datasets using the
                           standard process as well as subsetting existing datasets and overwriting
                           some aspects for a specific case.

        modify_input_files Scripts to modify CTSM input files. Specifically modifying the surface
                           datasets and mesh files.

        contrib ---------- Miscellaneous tools for pre or post processing of CTSM.
                           Typically these are contributed by anyone who has something
                           they think might be helpful to the community. They may not
                           be as well tested or supported as other tools.

II. Notes on building/running for each of the above tools:

    mksurfdata_esmf has a cime configure and CMake based build using the following files:

        gen_mksurfdata_build ---- Build mksurfdata_esmf
        src/CMakeLists.txt ------ Tells CMake how to build the source code
        Makefile ---------------- GNU makefile to link the program together
        cmake ------------------- CMake macros for finding libraries

    mkmapgrids, and site_and_regional only contain scripts so don't have the above build files.

    Some tools have copies of files from other directories -- see the README.filecopies
    file for more information on this.

    Tools may also have files with the directory name followed by namelist to provide sample namelists.

	<directory>.namelist ------ Namelist to create a global file.

    These files are also used by the test scripts to test the tools (see the
    README.testing) file.

[!NOTE] Be sure to change the path of the datasets references by these namelists to
    point to where you have exported your CESM inputdata datasets.

III. Process sequence to create input datasets needed to run CTSM

    1. Create ESMF MESH grid files (if needed)

       a. For standard resolutions these files will already be created. (done)

       b. Run `tools/site_and_regional/subset_data point` to create single-point datasets

       This creates just the fsurdat file as MESH files are NOT needed for single-point cases.

       c. Run `tools/site_and_regional/subset_data region` to create regional datasets subset from a global dataset

       This creates both the fsurdat file and MESH file needed to run.

       d. General custom grid

        You'll need to convert or create MESH grid files on your own (using scripts
        or other tools) for the general case where you have an unstructured grid, or
        a grid that is not regular in latitude and longitude. And that grid is custom
        and not merely subset from one of the global grids.

    2. Create surface datasets with mksurfdata_esmf on Derecho
        (See mksurfdata_esmf/README.md for more help on doing this)

       - gen_mksurfdata_build to build
       - gen_mksurfdata_namelist to build the namelist
       - gen_mksurfdata_jobscript_single to build a batch script to run on Derecho
       - Submit the batch script just created above

       - This step uses the results of step 1. entered into the XML database.
       - If datasets were NOT entered into the XML database, set the resolution
         by entering the mesh file using the options: --model-mesh --model-mesh-nx --model-mesh-ny

       Example: for 0.9x1.25 resolution for 1850

``` shell
       # On Derecho
       cd mksurfdata_esmf
       ./gen_mksurfdata_build
       ./gen_mksurfdata_namelist --res 0.9x1.25 --start-year 1850 --end-year 1850
       ./gen_mksurfdata_jobscript_single --number-of-nodes 24 --tasks-per-node 12 --namelist-file target.namelist
       qsub mksurfdata_jobscript_single.sh
```

    3. Add new files to XML data or using user_nl_clm (optional)

       See notes on doing this in step (1) above.

IV.  Notes on which input datasets are needed for CTSM

       global or regional grids
         - need fsurdata
         - need mesh files in env_run.xml ATM_DOMAIN_MESH and LND_DOMAIN_MESH

       single-point grids
        - Just need fsurdata
