.. include:: ../substitutions.rst

.. _observational-sites-datasets:

*******************************
Observational Sites Datasets
*******************************

A way to customize the model input datasets for single-point simulations is to customize the inputs to the tools that create the datasets. Another way is to overwrite the data after you have generated the datasets with default inputs. Depending on the tool and the dataset, one way or the other may be simpler.

Files that you may customize include ``fatmlndfrc``, ``fsurdat``, ``faerdep`` (for DATM), and ``stream_fldfilename_ndep``. Though generally to customize ``fsurdat`` files, you could modify the inputs needed by the ``mksurfdata_esmf`` tool, for single-point simulations we recommend a different method in :ref:`generic_single_point_runs`.

Another aspect of customizing your input datasets is customizing the input atmospheric forcing datasets; see :ref:`generic_single_point_runs` for information on this.
