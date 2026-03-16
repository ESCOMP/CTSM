.. include:: ../substitutions.rst

.. _customizing-the-mizuRoute-namelist:

===================================
 Customizing the mizuRoute namelist
===================================

When running "I" compsets with CLM you can use the mizuRoute ROF model to model the river flow. Compsets with "Mz" in the alias name or ``_MIZUROUTE_`` in the long compset name use the mizuRoute ROF model.

1. **mizuRoute control file** (``mizuRoute.control``)
2. **mizuRoute Namelist** (``mizuRoute_in``)

The `mizuRoute User's Guide Documentation <https://mizuroute.readthedocs.io/en/latest/users_guide/index.html>`_ gives the details of all the options for mizuRoute namelist and control file.

To change options in the mizuRoute namelist you edit the user_nl_mizuroute file in your case directory similar to other namelist files.
To change options in the mizuRoute control file you edit the user_nl_mizuroute_control file in your case directory. The difference is the format requires the variable to be in ``<>`` brackets, the value given after whitespace, and a comment about the variable at the end.
So for example to set the history file output and new file frequency to daily you would add the following lines to your user_nl_mizuroute_control file:

.. code-block::

    <newFileFrequency>           daily   !  frequency for new output files (daily, monthly, yearly, single)
    <outputFrequency>            daily   !  time frequency used for temporal aggregation of output variables - numeric or daily, monthyly, or yearly

Normally, those are the only two mizuRoute control variables you would be likely to change. The other settings are set for you based on the mizuRoute grid you are running with. Some of the CTSM case settings also get placed for you in the mizuRoute control file, and those are things you should leave alone.

