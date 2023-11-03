.. include:: ../substitutions.rst

.. _single-point-configurations:

*****************************************
Single and Regional Grid Configurations
*****************************************

CLM allows you to set up and run cases with a single-point or a local region as well as global resolutions. This is often useful for running quick cases for testing, evaluating specific vegetation types, or land-units, or running with observed data for a specific site.

There are three different ways to do this for normal-supported site

``PTS_MODE``
  runs for a single point using global datasets.

``CLM_USRDAT_NAME``
  runs using your own datasets (single-point or regional).

``PTCLMmkdata``
  easily setup simulations to run for tower sites..

.. note:: ``PTS_MODE`` and ``PTCLMmkdata`` only works for a single point, while the other two options can also work for regional datasets as well.

.. _options-for-single-points:

=========================================
 Choosing the right single point options
=========================================

Running for a *normal supported site* is a great solution, if one of the supported single-point/regional datasets, is your region of interest (see :ref:`running-single-point-datasets`). All the datasets are created for you, and you can easily select one and run, out of the box with it using a supported resolution from the top level of the CESM scripts. The problem is that there is a very limited set of supported datasets. You can also use this method for your own datasets, but you have to create the datasets, and add them to the XML database in scripts, CLM and to the DATM. This is worthwhile if you want to repeat many multiple cases for a given point or region.

In general :ref:`pts_mode` is the quick and dirty method that gets you started without having to create datasets -- but has limitations. It's good for an initial attempt at seeing results for a point of interest, but since you can NOT restart with it, it's usage is limited. It is the quickest method as you can create a case for it directly from **create_newcase**. Although you can't restart, running a single point is very fast, and you can run for long simulation times even without restarts.

Next, ``CLM_USRDAT_NAME`` is the best way to setup cases quickly where you have to create your own datasets (see :ref:`running-single-point-datasets`). With this method you don't have to change DATM or add files to the XML database -- but you have to follow a strict naming convention for files. However, once the files are named and in the proper location, you can easily setup new cases that use these datasets. This is good for treating all the required datasets as a "group" and for a particular model version. For advanced CLM developers who need to track dataset changes with different model versions you would be best off adding these datasets as supported datasets with the "normal supported datasets" method.

Lastly *PTCLMmkdata* is a great way to easily create datasets, setup simulations and run simulations for tower sites. It takes advantage of both normal supported site functionality and CLM_USRDAT_NAME internally. A big advantage to it, is that it's one-stop shopping, it runs tools to create datasets, and runs **create_newcase** and sets the appropriate env variables for you. So you only have to learn how to run one tool, rather than work with many different ones. PTCLMmkdata is described in the next chapter, :ref:`running-PTCLM`.

Finally, if you also have meteorology data that you want to force your CLM simulations with you'll need to setup cases as described in :ref:`creating-your-own-singlepoint-dataset`. You'll need to create CLM datasets either according to ``CLM_USRDAT_NAME``. You may also need to modify DATM to use your forcing data. And you'll need to change your forcing data to be in a format that DATM can use. :ref:`converting-ameriflux-for-ptclmmkdata` tells you how to use AmeriFlux data for atmospheric forcing.

