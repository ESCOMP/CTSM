.. _observational-sites-datasets:

.. include:: ../substitutions.rst

*******************************
Observational Sites Datasets
*******************************

There are two ways to customize datasets for a particular observational site. 
The first is to customize the input to the tools that create the dataset, and the second is to over-write the default data after you've created a given dataset. 
Depending on the tool it might be easier to do it one way or the other. 
In `Table 3-1 <CLM-URL>`_ we list the files that are most likely to be customized and the way they might be customized. 
Of those files, the ones you are most likely to customize are: fatmlndfrc, fsurdat, faerdep (for DATM), and stream_fldfilename_ndep. 
Note **mksurfdata_map** as documented previously has options to overwrite the vegetation and soil types. 
For more information on this also see `the Section called Creating your own single-point/regional surface datasets in Chapter 5 <CLM-URL>`_. 
And PTCLM uses these methods to customize datasets see `Chapter 6 <CLM-URL>`_.

Another aspect of customizing your input datasets is customizing the input atmospheric forcing datasets. 
See `the Section called Running with your own atmosphere forcing in Chapter 5 <CLM-URL>`_ for more information on this. 
Also the chapter on PTCLM in `the Section called Converting AmeriFlux Data for use by PTCLM in Chapter 6 <CLM-URL>`_ has information on using the AmeriFlux tower site data as atmospheric forcing.
