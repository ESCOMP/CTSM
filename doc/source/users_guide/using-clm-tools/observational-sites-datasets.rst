.. include:: ../substitutions.rst

.. _observational-sites-datasets:

*******************************
Observational Sites Datasets
*******************************

A way to customize the model input datasets is to customize the inputs to the tools that create the datasets. One can also overwrite the data of already generated datasets. Depending on the tool, the dataset, and the planned simulation, one way or the other may be simpler, or a combination of both methods may make sense.

Files that you may customize include ``fatmlndfrc``, ``fsurdat``, ``faerdep`` (for DATM), and ``stream_fldfilename_ndep``. To customize ``fsurdat`` files, one may modify the inputs needed by the ``mksurfdata_esmf`` tool. In addition (or instead) we strongly recommend using the ``subset_data`` tool for single-point and regional simulations (see :ref:`generic_single_point_runs`). A combination of methods may make the most sense in some cases.

Another aspect of customizing your input datasets is customizing the input atmospheric forcing datasets; see :ref:`generic_single_point_runs` for information on this.
