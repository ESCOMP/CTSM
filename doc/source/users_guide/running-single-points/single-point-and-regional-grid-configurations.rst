.. include:: ../substitutions.rst

.. _single-point-configurations:

*****************************************
Single and Regional Grid Configurations
*****************************************

CLM allows you to set up and run cases with a single-point or a local region as well as global resolutions. This is often useful for running quick cases for testing, evaluating specific vegetation types, or land-units, or running with observed data for a specific site.

There are two different ways to do this for normal-supported site

``subset_data``
  runs for a single point using global datasets.

``CLM_USRDAT_NAME``
  runs using your own datasets (single-point or regional).

.. _options-for-single-points:

=========================================
 Choosing the right single point options
=========================================

Running for a *normal supported site* is a great solution, if one of the supported single-point/regional datasets, is your region of interest (see :ref:`running-single-point-datasets`). All the datasets are created for you, and you can easily select one and run, out of the box with it using a supported resolution from the top level of the CESM scripts. The problem is that there is a very limited set of supported datasets. You can also use this method for your own datasets, but you have to create the datasets, and add them to the XML database in scripts, CLM and to the DATM. This is worthwhile if you want to repeat many multiple cases for a given point or region.

In general :ref:`single_point_subset_data` is the quick and dirty method that gets you started, but it has limitations. It's good for an initial attempt at seeing results for a point of interest, but since you can NOT restart with it, its usage is limited. It is the quickest method as you can create a case for it directly from ``cime/scripts/create_newcase``. Although you can't restart, running a single point is very fast, and you can run for long simulation times even without restarts.

Next, ``CLM_USRDAT_NAME`` using ``subset_data`` is the best way to setup cases quickly where you have a simple tool to create your own datasets (see :ref:`single_point_subset_data`). With this method you don't have to change DATM or add files to the XML database. ``subset_data`` will create a usermod directory where you can store your files and the files needed to directly run a case.

Finally, if you also have meteorology data that you want to force your CLM simulations with you'll need to setup cases as described in :ref:`creating-your-own-singlepoint-dataset`. You'll need to create CLM datasets either according to ``CLM_USRDAT_NAME``. You may also need to modify DATM to use your forcing data. And you'll need to change your forcing data to be in a format that DATM can use.

