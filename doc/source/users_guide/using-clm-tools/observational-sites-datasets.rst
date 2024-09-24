.. include:: ../substitutions.rst

.. _observational-sites-datasets:

*******************************
Observational Sites Datasets
*******************************

.. todo::
    Update this.

There are two ways to customize datasets for a particular observational site. The first is to customize the input to the tools that create the dataset, and the second is to overwrite the default data after you've created a given dataset. Depending on the tool it might be easier to do it one way or the other. In Table :numref:`reqd-files-table` we list the files that are most likely to be customized and the way they might be customized. Of those files, the ones you are most likely to customize are: ``fatmlndfrc``, ``fsurdat``, ``faerdep`` (for DATM), and ``stream_fldfilename_ndep``. Note ``mksurfdata_esmf`` as documented previously has options to overwrite the vegetation and soil types. For more information on this also see :ref:`creating-your-own-singlepoint-dataset`.

Another aspect of customizing your input datasets is customizing the input atmospheric forcing datasets; see :ref:`creating-your-own-singlepoint-dataset` for more information on this.
