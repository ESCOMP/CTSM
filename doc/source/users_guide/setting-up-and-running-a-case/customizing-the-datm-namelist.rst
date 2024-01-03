.. include:: ../substitutions.rst

.. _customizing-the-datm-namelist:

===============================
 Customizing the DATM namelist
===============================

When running "I" compsets with CLM you use the DATM model to give atmospheric forcing data to CLM. There are two ways to customize DATM:

1. **DATM Main Namelist and Stream Namlist gorup** (``datm_in``)
2. **DATM stream files**

The `Data Model Documentation <http://esmci.github.io/cime/data_models/data-atm.html>`_ gives the details of all the options for the data models and for DATM specifically. It goes into detail on all namelist items both for DATM and for DATM streams. So here we won't list ALL of the DATM namelist options, nor go into great details about stream files. But, we will talk about a few of the different options that are relevant for running with CLM. All of the options for changing the namelists or stream files is done by editing the ``user_nl_datm`` file.

Because, they aren't useful for work with CLM we will NOT discuss any of the options for the main DATM namelist. Use the DATM Users Guide at the link above to find details of that. For the streams namelist we will discuss three items:

1. ``mapalgo``
#. ``taxmode``
#. ``tintalgo``

And for the streams file itself we will discuss:

offset
  Again everything else (and including the above items) are discussed in the Data Model User's Guide. Of the above the last three: offset, taxmode and tintalgo are all closely related and have to do with the time interpolation of the DATM data.

mapalgo
  ``mapalgo`` sets the spatial interpolation method to go from the DATM input data to the output DATM model grid. The default is ``bilinear``. For ``CLM1PT`` we set it to ``nn`` to just select the nearest neighbor. This saves time and we also had problems running the interpolation for single-point mode.

taxmode
  ``taxmode`` is the time axis mode. For CLM we usually have it set to ``cycle`` which means that once the end of the data is reached it will start over at the beginning. The extend modes is used have it use the last time-step of the forcing data once it reaches the end of forcing data (or use the first time-step before it reaches where the forcing data starts). See the warning below about the extend mode.

.. warning:: *THE extend OPTION NEEDS TO BE USED WITH CAUTION!* It is only invoked by default for the CLM1PT mode and is only intended for the supported urban datasets to extend the data for a single time-step. If you have the model *run extensively through periods in this mode you will effectively be repeating that last time-step over that entire period*. This means the output of your simulation will be worthless.

offset (in the stream file)
  ``offset`` is the time offset in seconds to give to each stream of data. Normally it is NOT used because the time-stamps for data is set correctly for each stream of data. Note, the ``offset`` may NEED to be adjusted depending on the ``taxmode`` described above, or it may need to be adjusted to account for data that is time-stamped at the END of an interval rather than the middle or beginning of interval. The ``offset`` can is set in the stream file rather than on the stream namelist. For data with a ``taxmode`` method of ``coszen`` the time-stamp needs to be for the beginning of the interval, while for other data it should be the midpoint. The ``offset`` can be used to adjust the time-stamps to get the data to line up correctly.

tintalgo
  ``tintalgo`` is the time interpolation algorithm. For CLM we usually use one of three modes: ``coszen``, ``nearest``, or ``linear``. We use ``coszen`` for solar data, nearest for precipitation data, and linear for everything else. If your data is half-hourly or hourly, nearest will work fine for everything. The ``coszen`` scaling is useful for longer periods (three hours or more) to try to get the solar to match the cosine of the solar zenith angle over that longer period of time. If you use linear for longer intervals, the solar will cut out at night-time anyway, and the straight line will be a poor approximation of the cosine of the solar zenith angle of actual solar data. nearest likewise would be bad for longer periods where it would be much higher than the actual values.

.. note:: For ``coszen`` the time-stamps of the data should correspond to the beginning of the interval the data is measured for. Either make sure the time-stamps on the datafiles is set this way, or use the ``offset`` described above to set it.

.. note:: For nearest and linear the time-stamps of the data should correspond to the middle of the interval the data is measured for. Either make sure the time-stamps on the datafiles is set this way, or use the ``offset`` described above to set it.

In the sections below we go over each of the relevant ``DATM_MODE`` options and what the above DATM settings are for each. This gives you examples of actual usage for the settings. We also describe in what ways you might want to customize them for your own case.

--------------------------------------
CLMGSWP3v1 mode and it's DATM settings
--------------------------------------

In ``CLMGSWP3v1`` mode the GSWP3 NCEP forcing dataset is used and all of it's data is on a 3-hourly interval. Like ``CLM_QIAN`` the dataset is divided into those three data streams: solar, precipitation, and everything else (temperature, pressure, humidity, Long-Wave down and wind). The time-stamps of the data were also adjusted so that they are the beginning of the interval for solar, and the middle for the other two. Because, of this the ``offset`` is set to zero, and the ``tintalgo`` is: ``coszen``, ``nearest``, and ``linear`` for the solar, precipitation and other data respectively. ``taxmode`` is set to ``cycle`` and ``mapalgo`` is set to ``bilinear`` so that the data is spatially interpolated from the input exact half degree grid to the grid the atmosphere model is being run at (to run at this same model resolution use the 360x720cru_360x720cru resolution).

----------------------------------------
CLMCRUNCEPv7 mode and it's DATM settings
----------------------------------------

In ``CLMCRUNCEPv7`` mode the CRUNCEP dataset is used and all of it's data is on a 6-hourly interval. Like ``CLM_QIAN`` the dataset is divided into those three data streams: solar, precipitation, and everything else (temperature, pressure, humidity and wind). The time-stamps of the data were also adjusted so that they are the beginning of the interval for solar, and the middle for the other two. Because, of this the ``offset`` is set to zero, and the ``tintalgo`` is: ``coszen``, ``nearest``, and ``linear`` for the solar, precipitation and other data respectively. ``taxmode`` is set to ``cycle`` and ``mapalgo`` is set to ``bilinear`` so that the data is spatially interpolated from the input exact half degree grid to the grid the atmosphere model is being run at (to run at this same model resolution use the 360x720cru_360x720cru resolution)... note:: The "everything else" data stream (of temperature, pressure, humidity and wind) also includes the data for longwave downward forcing as well. Our simulations showed sensitivity to this field, so we backed off in using it, and let DATM calculate longwave down from the other fields.

For more information on CRUNCEP forcing see `http://dods.extra.cea.fr/data/p529viov/cruncep/ <http://dods.extra.cea.fr/data/p529viov/cruncep/>`_.

.. _clmcruncep-and-its-datm:

--------------------------------------
CLMCRUNCEP mode and it's DATM settings
--------------------------------------

``CLMCRUNCEP`` is similar to the ``CLMCRUNCEPv7`` mode above, except it uses Version 4 of the CRUNCEP data rather than version 7.

.. _clmqian-and-its-datm:

------------------------------------
CLM_QIAN mode and it's DATM settings
------------------------------------

In ``CLM_QIAN`` mode the Qian dataset is used which has 6-hourly solar and precipitation data, and 3-hourly for everything else. The dataset is divided into those three data streams: solar, precipitation, and everything else (temperature, pressure, humidity and wind). The time-stamps of the data were also adjusted so that they are the beginning of the interval for solar, and the middle for the other two. Because, of this the ``offset`` is set to zero, and the ``tintalgo`` is: ``coszen``, ``nearest``, and ``linear`` for the solar, precipitation and other data respectively. ``taxmode`` is set to ``cycle`` and ``mapalgo`` is set to ``bilinear`` so that the data is spatially interpolated from the input T62 grid to the grid the atmosphere model is being run at.

Normally you wouldn't customize the ``CLM_QIAN`` settings, but you might replicate it's use for your own global data that had similar temporal characteristics.

.. _clm1pt-and-its-datm:

----------------------------------
CLM1PT mode and it's DATM settings
----------------------------------

In ``CLM1PT`` mode the model is assumed to have half-hourly or hourly data for a single-point. For the supported datasets that is exactly what it has. But, if you add your own data you may need to make adjustments accordingly. Using the ``CLM_USRDAT_NAME`` resolution you can easily extend this mode for your own datasets that may be regional or even global and could be at different temporal frequencies. If you do so you'll need to make adjustments to your DATM settings. The dataset has all data in a single stream file. The time-stamps of the data were also adjusted so that they are at the middle of the interval. Because, of this the ``offset`` is set to zero, and the ``tintalgo`` is set to ``nearest``. ``taxmode`` is set to ``extend`` and ``mapalgo`` is set to ``nn`` so that simply the nearest point is used.

If you are using your own data for this mode and it's not at least hourly you'll want to adjust the DATM settings for it. If the data is three or six hourly, you'll need to divide it up into separate streams like in ``CLM_QIAN`` mode which will require fairly extensive changes to the DATM namelist and streams files. For an example of doing this see :ref:`eg-sim-data-from-prev-sim`.

.. _cplhistforcing:

------------------------------------------
CPLHISTForcing mode and it's DATM settings
------------------------------------------

In ``CPLHISTForcing`` mode the model is assumed to have 3-hourly for a global grid from a previous CESM simulation. Like ``CLM_QIAN`` mode the data is divided into three streams: one for precipitation, one for solar, and one for everything else. The time-stamps for Coupler history files for CESM is at the end of the interval, so the ``offset`` needs to be set in order to adjust the time-stamps to what it needs to be for the ``tintalgo`` settings. For precipitation ``taxmode`` is set to ``nearest`` so the ``offset`` is set to ``-5400`` seconds so that the ending time-step is adjusted by an hour and half to the middle of the interval. For solar ``taxmode`` is set to ``coszen`` so the offset is set to ``-10800`` seconds so that the ending time-step is adjust by three hours to the beginning of the interval. For everything else ``taxmode`` is set to ``linear`` so the offset is set to ``-5400`` seconds so that the ending time-step is adjusted by an hour and half to the middle of the interval. For an example of such a case see :ref:`running-with-moar-data`.

Normally you wouldn't modify the DATM settings for this mode. However, if you had data at a different frequency than 3-hours you would need to modify the ``offset`` and possibly the ``taxmode``. The other two things that you might modify would be the path to the data or the domain file for the resolution (which is currently hardwired to f09). For data at a different input resolution you would need to change the domain file in the streams file to use a domain file to the resolution that the data comes in on.
