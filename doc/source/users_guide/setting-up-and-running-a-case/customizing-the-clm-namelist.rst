.. include:: ../substitutions.rst

.. _customizing-a-case:

============================
 Customizing CLM's namelist
============================

Once a case has run ``case.setup``, we can then customize the case further, by editing the run-time namelist for CLM. First let's list the definition of each namelist item and their valid values, and then we'll list the default values for them. Next for some of the most used or tricky namelist items we'll give examples of their use, and give you example namelists that highlight these features.

In the following, various examples of namelists are provided that feature the use of different namelist options to customize a case for particular uses. Most the examples revolve around how to customize the output history fields. This should give you a good basis for setting up your own CLM namelist.

.. _def-nl-items-and-defaults:

-----------------------------------------------------
Definition of Namelist items and their default values
-----------------------------------------------------

Here we point to you where you can find the definition of each namelist item and separately the default values for them. The default values may change depending on the resolution, land-mask, simulation-year and other attributes. Both of these files are viewable in your web browser, and then expand each in turn.

1. `Definition of Namelists <https://github.com/ESCOMP/CTSM/blob/master/bld/namelist_files/namelist_definition_ctsm.xml>`_

2. `Default values of each Namelist Item <https://github.com/ESCOMP/CTSM/blob/master/bld/namelist_files/namelist_defaults_ctsm.xml>`_

List of fields that can be added to your output history files by namelist
-------------------------------------------------------------------------

One set of the namelist items allows you to add fields to the output history files: ``hist_fincl1``, ``hist_fincl2``, ``hist_fincl3``, ``hist_fincl4``, ``hist_fincl5``, and ``hist_fincl6``. The :doc:`history_fields_nofates` and :doc:`history_fields_fates` files list all of the history fields available and gives the long-name and units for each.

---------------------------------------------
Examples of using different namelist features
---------------------------------------------

Below we will give examples of user namelists that activate different commonly used namelist features. We will discuss the namelist features in different examples and then show a user namelist that includes an example of the use of these features. First we will show the default namelist that doesn't activate any user options.

The default namelist
--------------------

Here we give the default namelist as it would be created for an "I1850Clm50BgcCropCru" compset at 0.9x1.25 resolution with a gx1v7 land-mask on cheyenne. To edit the namelist you would edit the ``user_nl_clm`` user namelist with just the items you want to change. For simplicity we will just show the CLM namelist and NOT the entire file. In the sections below, for simplicity we will just show the user namelist (``user_nl_clm``) that will add (or modify existing) namelist items to the namelist.

Example 1-2. Default CLM Namelist
---------------------------------
::

     &clm_inparm
      albice = 0.60,0.40
      co2_ppmv = 284.7
      co2_type = 'constant'
      create_crop_landunit = .false.
      dtime = 1800
      fatmlndfrc = '/glade/p/cesm/cseg/inputdata/share/domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc'
      finidat = '$DIN_LOC_ROOT/lnd/clm2/initdata_map/clmi.I1850Clm50BgcCropCru.0241-01-01.0.9x1.25_g1v6_simyr1850_c130531.nc'
      fpftcon = '/glade/p/cesm/cseg/inputdata/lnd/clm2/pftdata/pft-physiology.c130503.nc'
      fsnowaging = '/glade/p/cesm/cseg/inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc'
      fsnowoptics = '/glade/p/cesm/cseg/inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c090915.nc'
      fsurdat = '/glade/p/cesm/cseg/inputdata/lnd/clm2/surfdata_map/surfdata_0.9x1.25_simyr1850_c130415.nc'
      maxpatch_glc = 0
      more_vertlayers = .false.
      nsegspc = 20
      spinup_state = 0
      urban_hac = 'ON'
      urban_traffic = .false.
     /
     &ndepdyn_nml
      ndepmapalgo = 'bilinear'
      stream_fldfilename_ndep = '/glade/p/cesm/cseg/inputdata/lnd/clm2/ndepdata/fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc'
      stream_year_first_ndep = 1850
      stream_year_last_ndep = 1850
     /
     &popd_streams
      popdensmapalgo = 'bilinear'
      stream_fldfilename_popdens = '$DIN_LOC_ROOT/lnd/clm2/firedata/clmforc.Li_2012_hdm_0.5x0.5_AVHRR_simyr1850-2010_c130401.nc'
      stream_year_first_popdens = 1850
      stream_year_last_popdens = 1850
     /
     &light_streams
      lightngmapalgo = 'bilinear'
      stream_fldfilename_lightng = '/glade/p/cesm/cseg/inputdata/atm/datm7/NASA_LIS/clmforc.Li_2012_climo1995-2011.T62.lnfm_c130327.nc'
      stream_year_first_lightng = 0001
      stream_year_last_lightng = 0001
     /
     &clm_hydrology1_inparm
     /
     &clm_soilhydrology_inparm
     /
     &ch4par_in
      fin_use_fsat = .true.
     /

Adding/removing fields on your primary history file
---------------------------------------------------

The primary history files are output monthly, and contain an extensive list of fieldnames, but the list of fieldnames can be added to using ``hist_fincl1`` or removed from by adding fieldnames to ``hist_fexcl1``. A sample user namelist ``user_nl_clm`` adding few new fields (cosine of solar zenith angle, and solar declination) and excluding a few standard fields is (ground temperature, vegetation temperature, soil temperature and soil water).:

Example 1-3. Example user_nl_clm namelist adding and removing fields on primary history file
--------------------------------------------------------------------------------------------
::

   hist_fincl1 = 'COSZEN', 'DECL'
   hist_fexcl1 = 'TG', 'TV', 'TSOI', 'H2OSOI'

Adding auxiliary history files and changing output frequency
------------------------------------------------------------

The ``hist_fincl2`` through ``hist_fincl6`` set of namelist variables add given history fieldnames to auxiliary history file "streams", and ``hist_fexcl2`` through ``hist_fexcl6`` set of namelist variables remove given history fieldnames from history file auxiliary "streams". A history "stream" is a set of history files that are produced at a given frequency. By default there is only one stream of monthly data files. To add more streams you add history fieldnames to ``hist_fincl2`` through ``hist_fincl6``. The output frequency and the way averaging is done can be different for each history file stream. By default the primary history files are monthly and any others are daily. You can have up to six active history streams, but you need to activate them in order. So if you activate stream "6" by setting ``hist_fincl6``, but if any of ``hist_fincl2`` through ``hist_fincl5`` are unset, only the history streams up to the first blank one will be activated.

The frequency of the history file streams is given by the namelist variable ``hist_nhtfrq`` which is an array of rank six for each history stream. The values of the array ``hist_nhtfrq`` must be integers, where the following values have the given meaning:

*Positive value* means the output frequency is the number of model steps between output. *Negative value* means the output frequency is the absolute value in hours given (i.e -1 would mean an hour and -24 would mean a full day). Daily (-24) is the default value for all auxiliary files. *Zero* means the output frequency is monthly. This is the default for the primary history files.

The number of samples on each history file stream is given by the namelist variable ``hist_mfilt`` which is an array of rank six for each history stream. The values of the array ``hist_mfilt`` must be positive integers. By default the primary history file stream has one time sample on it (i.e. output is to separate monthly files), and all other streams have thirty time samples on them.

A sample user namelist ``user_nl_clm`` turning on four extra file streams for output: daily, six-hourly, hourly, and every time-step, leaving the primary history files as monthly, and changing the number of samples on the streams to: yearly (12), thirty, weekly (28), daily (24), and daily (48) is:

Example: user_nl_clm namelist adding auxiliary history files and changing output frequency
------------------------------------------------------------------------------------------
::

   hist_fincl2 = 'TG', 'TV'
   hist_fincl3 = 'TG', 'TV'
   hist_fincl4 = 'TG', 'TV'
   hist_fincl5 = 'TG', 'TV'
   hist_nhtfrq = 0, -24, -6, -1, 1
   hist_mfilt  = 12, 30, 28, 24, 48

Removing all history fields
---------------------------

Sometimes for various reasons you want to remove all the history fields either because you want to do testing without any output, or you only want a very small custom list of output fields rather than the default extensive list of fields. By default only the primary history files are active, so technically using ``hist_fexcl1`` explained in the first example, you could list ALL of the history fields that are output in ``hist_fexcl1`` and then you wouldn't get any output. However, as the list is very extensive this would be a cumbersome thing to do. So to facilitate this ``hist_empty_htapes`` allows you to turn off all default output. You can still use ``hist_fincl1`` to turn your own list of fields on, but you then start from a clean slate. A sample user namelist ``user_nl_clm`` turning off all history fields and then activating just a few selected fields (ground and vegetation temperatures and absorbed solar radiation) is:

Example 1-5. Example user_nl_clm namelist removing all history fields
---------------------------------------------------------------------
::

   hist_empty_htapes = .true.
   hist_fincl1 = 'TG', 'TV', 'FSA'

Various ways to change history output averaging flags
-----------------------------------------------------

There are two ways to change the averaging of output history fields. The first is using ``hist_avgflag_pertape`` which gives a default value for each history stream, the second is when you add fields using ``hist_fincl*``, you add an averaging flag to the end of the field name after a colon (for example ``TSOI:X`` would output the maximum of ``TSOI``). The types of averaging that can be done are:

- ``A`` Average, over the output interval.
- ``I`` Instantaneous, output the value at the output interval.
- ``X`` Maximum, over the output interval.
- ``M`` Minimum, over the output interval.

The default averaging depends on the specific fields, but for most fields is an average. A sample user namelist ``user_nl_clm`` making the monthly output fields all averages (except ``TSOI`` for the first two streams and ``FIRE`` for the 5th stream), and adding auxiliary file streams for instantaneous (6-hourly), maximum (daily), minimum (daily), and average (daily). For some of the fields we diverge from the per-tape value given and customize to some different type of optimization.

Example: user_nl_clm namelist with various ways to average history fields
-------------------------------------------------------------------------------------
::

   hist_empty_htapes = .true.
   hist_fincl1 = 'TSOI:X', 'TG',   'TV',   'FIRE',   'FSR', 'FSH',
		 'EFLX_LH_TOT', 'WT'
   hist_fincl2 = 'TSOI:X', 'TG',   'TV',   'FIRE',   'FSR', 'FSH',
		 'EFLX_LH_TOT', 'WT'
   hist_fincl3 = 'TSOI',   'TG:I', 'TV',   'FIRE',   'FSR', 'FSH',
		 'EFLX_LH_TOT', 'WT'
   hist_fincl4 = 'TSOI',   'TG',   'TV:I', 'FIRE',   'FSR', 'FSH',
		 'EFLX_LH_TOT', 'WT'
   hist_fincl5 = 'TSOI',   'TG',   'TV',   'FIRE:I', 'FSR', 'FSH',
		 'EFLX_LH_TOT', 'WT'
   hist_avgflag_pertape = 'A', 'I', 'X',   'M', 'A'
   hist_nhtfrq = 0, -6, -24, -24, -24

In the example we put the same list of fields on each of the tapes: soil-temperature, ground temperature, vegetation temperature, emitted longwave radiation, reflected solar radiation, sensible heat, total latent-heat, and total water storage. We also modify the soil temperature for the primary and secondary auxiliary tapes by outputting them for a maximum instead of the prescribed per-tape of average and instantaneous respectively. For the tertiary auxiliary tape we output ground temperature instantaneous instead of as a maximum, and for the fourth auxiliary tape we output vegetation temperature instantaneous instead of as a minimum. Finally, for the fifth auxiliary tapes we output ``FIRE`` instantaneously instead of as an average.

.. note:: We also use ``hist_empty_htapes`` as in the previous example, so we can list ONLY the fields that we want on the primary history tapes.

Outputting history files as a vector in order to analyze the plant function types within gridcells
--------------------------------------------------------------------------------------------------

By default the output to history files are the grid-cell average of all land-units, and vegetation types within that grid-cell, and output is on the full 2D latitude/longitude grid with ocean masked out. Sometimes it's important to understand how different land-units or vegetation types are acting within a grid-cell. The way to do this is to output history files as a 1D-vector of all land-units and vegetation types. In order to display this, you'll need to do extensive post-processing to make sense of the output. Often you may only be interested in a few points, so once you figure out the 1D indices for the grid-cells of interest, you can easily view that data. 1D vector output can also be useful for single point datasets, since it's then obvious that all data is for the same grid cell.

To do this you use ``hist_dov2xy`` which is an array of rank six for each history stream. Set it to ``.false.`` if you want one of the history streams to be a 1D vector. You can also use ``hist_type1d_pertape`` if you want to average over all the: Plant-Function-Types, columns, land-units, or grid-cells. A sample user namelist ``user_nl_clm`` leaving the primary monthly files as 2D, and then doing grid-cell (GRID), column (COLS), and no averaging over auxiliary tapes output daily for a single field (ground temperature) is:

Example: user_nl_clm namelist outputting some files in 1D Vector format
-----------------------------------------------------------------------
::

   hist_fincl2 = 'TG'
   hist_fincl3 = 'TG'
   hist_fincl4 = 'TG'
   hist_fincl5 = 'TG'
   hist_fincl6 = 'TG'
   hist_dov2xy = .true., .false., .false., .false.
   hist_type2d_pertape = ' ', 'GRID', 'COLS', ' '
   hist_nhtfrq = 0, -24, -24, -24

.. warning:: ``LAND`` and ``COLS`` are also options to the pertape averaging, but currently there is a bug with them and they fail to work.

.. note:: Technically the default for ``hist_nhtfrq`` is for primary files output monthly and the other auxiliary tapes for daily, so we don't actually have to include ``hist_nhtfrq``, we could use the default for it. Here we specify it for clarity.

Visualizing global 1D vector files will take effort. You'll probably want to do some post-processing and possibly just extract out single points of interest to see what is going on. Since the output is a 1D vector of only land points, traditional plots won't be helpful. The number of points per grid-cell will also vary for anything but grid-cell averaging. You'll need to use the output fields ``pfts1d_ixy``, and ``pfts1d_jxy``, to get the mapping of the fields to the global 2D array. ``pfts1d_itype_veg`` gives you the PFT number for each PFT. Most likely you'll want to do this analysis in a data processing tool (such as NCL, Matlab, Mathmatica, IDL, etc. that is able to read and process NetCDF data files).
