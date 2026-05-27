.. include:: ../substitutions.rst

.. _observational-sites-datasets:

*******************************
Observational Sites Datasets
*******************************

There are two ways to customize datasets for a particular observational site. The first is to customize the input to the tools that create the dataset, and the second is to overwrite the default data after you've created a given dataset. Depending on the tool it might be easier to do it one way or the other. Files that you may customize include ``fatmlndfrc``, ``fsurdat``, ``faerdep`` (for DATM), and ``stream_fldfilename_ndep``. To generate custom ``fsurdat`` files, you may modify the inputs needed by the ``mksurfdata_esmf`` tool. Also see :ref:`generic_single_point_runs` for relevant information pertaining to single-point simuations.

Another aspect of customizing your input datasets is customizing the input atmospheric forcing datasets; see :ref:`generic_single_point_runs` for information on this.
