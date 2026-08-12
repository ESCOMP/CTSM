.. include:: ../substitutions.rst

.. _adding-resolutions:

========================
 Adding New Resolutions
========================

In the last chapter we gave the details on how to create new files for input into CLM. These files could be either global resolutions, regional-grids or even a single grid point. If you want to easily have these files available for continued use in your development you will then want to include them in the build-namelist database so that build-namelist can easily find them for you. You can deal with them, just by putting the settings in the ``user_nl_clm`` namelist file, or by using ``CLM_USRDAT_NAME``. Another way to deal with them is to enter them into the database for build-namelist, so that build-namelist can find them for you. This keeps one central database for all your files, rather than having multiple locations to keep track of files. If you have a LOT of files to keep track of it also might be easier than keeping track by hand, especially if you have to periodically update your files. If you just have a few quick experiments to try, for a short time period you might be best off using the other methods mentioned above.

There are two parts to adding files to the build-namelist database. The first part is adding new resolution names which is done in the ``$CTSMROOT/bld/namelist_files/namelist_definition_ctsm.xml`` file. You can then use the new resolution by using ``CLM_USRDAT_NAME``. If you also want to be able to give the resolution to ``$CTSMROOT/cime/scripts/create_newcase`` -- you'll need to add the grid to the ``$CIMEROOT/config/cesm/config_grid.xml`` file.

The second part is actually adding the new filenames which is done in the ``$CTSMROOT/bld/namelist_files/namelist_definition_ctsm.xml`` file. If you aren't adding any new resolutions, and you are just changing the files for existing resolutions, you don't need to edit the namelist_definition file.

