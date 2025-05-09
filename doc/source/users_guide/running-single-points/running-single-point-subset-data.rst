.. include:: ../substitutions.rst

.. _single_point_subset_data:

****************************************
Running a single point using global data
****************************************

``subset_data`` enables you to run the model using global datasets, but just picking a single point from those datasets and operating on it. It can be a very quick way to do fast simulations and get a quick turnaround.

Subset the data
------------------

For single-point cases, you need to subset a surface dataset and (optionally) DATM data. The Python script to subset this data can be found in the CTSM repository at ``tools/site_and_regional/subset_data``.

Note that you will need to have a python environment set up that includes the packages ``scipy``, ``xarray``, and ``numpy``. If you have conda or miniconda installed, you can create a conda environment for this and other CTSM python tools using the script ``py_env_create`` at the top level of your CTSM checkout. See :ref:`using-ctsm-pylib` for more information.

To subset surface data and climate forcings (DATM) for a single point, use the command:

.. code:: shell

   tools/site_and_regional/subset_data point \
      --lat $my_lat --lon $my_lon --lon-type $my_lon_type --site $my_site_name \
      --create-surface --create-datm \
      --datm-syr $my_start_year --datm-eyr $my_end_year \
      --create-user-mods --outdir $my_output_dir

-  ``$my_lat``: latitude of point, *must be between -90 and 90 degrees*. E.g., Boulder, CO, USA: 40.
-  ``$my_lon``: longitude of point. *Must be between -180 and 360 degrees.* E.g., Boulder, CO, USA: 255 or -105.
-  ``$my_lon_type``: 180 if your longitude is in the [-180, 180] format (i.e., centered at the Prime/0th Meridian); 360 if it's in the [0, 360] format (i.e., centered at the 180th Meridian). Note that ``--lon-type $my_lon_type`` is not necessary if your longitude is unambiguous---i.e., it's only needed if your longitude is in the range [0, 180].
-  ``$my_site_name``: name of site, *used for file naming*
-  ``$my_start_year``: start year for DATM data to subset, *default between 1901 and 2014*
-  ``$my_end_year``: end year for DATM data to subset, *default between 1901 and 2014; the default CRUJRA2024 DATM data ends in 2023, while the old default GSWP3 ends in 2015; see note below about switching the default DATM data*
-  ``$my_output_dir``: output directory to place the subset data and user_mods directory. This should be something specific to *just* your data for ``$my_site_name``.

You can also have the script subset land-use data. See the help (``tools/site_and_regional/subset_data --help``) for all argument options.

.. note::
   This script defaults to subsetting specific surface, domain, and land-use files and the CRUJRA2024 DATM data, and can currently only be run as-is on Derecho. If you're not on Derecho, use ``--inputdata-dir`` to specify where the top level of your CESM input data is. Also, to subset GSWP3 instead of CRUJRA2024 DATM data, you currently need to hardwire ``datm_type = "datm_gswp3"`` (instead of the default ``"datm_crujra"``) in ``python/ctsm/subset_data.py``.

The ``--create-user-mods`` command tells the script to set up a user mods directory in your specified ``$my_output_dir`` and to specify the required ``PTS_LAT`` and ``PTS_LON`` settings. You can then use this user mods directory to set up your CTSM case, as described below.

Create the case
------------------

You can use the user mods directory set up in the previous subset data step to tell CIME/CTSM where your subset files are located.

.. code:: shell

   cime/scripts/create_newcase --case $my_case_name --res CLM_USRDAT \
      --compset $compset --run-unsupported \
      --user-mods-dirs $my_output_dir/user_mods

-  ``$my_case_name``: the path of the case directory you want to create
-  ``$compset``: the compset you would like to use (for example, ``I2000Clm60Bgc``)
-  Note the use of ``$my_output_dir/user_mods`` which is the ``user_mods/`` directory that the subset data script set up within your specified ``$my_output_dir``.

Note that ``./case.setup`` on Derecho will automatically set queue to ``develop`` and walltime to one hour. You might need a longer walltime, but the maximum walltime for ``develop`` is one hour. To change it to two hours on Derecho:

.. code:: shell

   ./xmlchange --subgroup case.run JOB_QUEUE=main,JOB_WALLCLOCK_TIME=2:00:00
