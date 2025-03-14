.. running-with-custom-crop-calendars:

.. include:: ../substitutions.rst

=======================================
 Running with custom crop calendars
=======================================

Since CLM5.1, functionality has been added to enable the customization of crop sowing window and maturity requirements.

Sowing window
-------------
Crops are allowed to be sown only within certain date windows. By default, these are static in time and specific to each crop in each hemisphere. These values can be inspected by looking at (and changed by modifying) the following variables on the parameter file:

- ``min_NH_planting_date`` (start of sowing window, northern hemisphere)
- ``max_NH_planting_date`` (end of sowing window, northern hemisphere)
- ``min_SH_planting_date`` (start of sowing window, southern hemisphere)
- ``max_SH_planting_date`` (end of sowing window, southern hemisphere)

However, geographically- and temporally-varying maps can also be used to prescribe sowing window like so (in ``user_nl_clm``):
::

   ! Input files with maps of the start and end date (1-365) of every crop
   stream_fldFileName_swindow_start = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/cropdata/calendars/processed/swindow_starts_ggcmi_crop_calendar_phase3_v1.01.2000-2000.20231005_145103.nc'
   stream_fldFileName_swindow_end = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/cropdata/calendars/processed/swindow_ends_ggcmi_crop_calendar_phase3_v1.01.2000-2000.20231005_145103.nc'

   ! A mesh file matching the resolution of the sowing window datasets
   stream_meshfile_cropcal = '/glade/p/cesmdata/cseg/inputdata/share/meshes/360x720_120830_ESMFmesh_c20210507_cdf5.nc'

   ! First and last years on the sowing window datasets
   stream_year_first_cropcal_swindows = 2000
   stream_year_last_cropcal_swindows = 2000

Sowing date
-----------
Specific sowing *dates* can be prescribed for any crop in any gridcell by setting the start and end dates of its sowing windows to the same value. The simplest way to do this for all crops in all gridcells is to specify the same file for both ``stream_fldFileName_swindow_start`` and ``stream_fldFileName_swindow_end``.

.. note:: In cells with prescribed sowing dates, the usual weather- and climate-based criteria for determining whether planting is allowed are ignored. The crop will be planted on the prescribed day no matter what.

Maturity requirements
---------------------
The heat unit accumulation required for a crop to reach maturity (and thus be harvested) is typically determined by a formula with crop-specific parameters that are specified on the parameter file. However, geographically- and temporally-varying maps of maturity requirement (in units of degree-days) can also be specified using the ``user_nl_clm`` input variable ``stream_fldFileName_cultivar_gdds``. (Note that ``stream_meshfile_cropcal``, ``stream_year_first_cropcal_cultivar_gdds``, and ``stream_year_last_cropcal_cultivar_gdds``---see above---are all also required.)

Generating maturity requirements
--------------------------------
For phase 3 of the Global Gridded Crop Model Intercomparison (GGCMI), maturity requirements should be the average (over the 1980--2009 growing seasons) growing degree-days accumulated between specified observation-derived sowing and harvest dates. In CLM, this requires the use of a special "GDD-generating" run and some postprocessing.

In a GDD-generating run, crops are planted on the specified sowing dates and are then left in the field---regardless of when they reach maturity and ignoring maximum growing season length---for 364 days. This is set up like so in ``user_nl_clm``:
::

   ! Variables that we introduced above
   stream_fldFileName_swindow_start = '/path/to/sowing_date_file.nc'
   stream_fldFileName_swindow_end = '/path/to/sowing_date_file.nc'
   stream_meshfile_cropcal = '/path/to/mesh_file.nc'
   stream_year_first_cropcal_swindows = YEAR
   stream_year_last_cropcal_swindows = YEAR

   ! Special settings for "GDD-generating" run
   generate_crop_gdds = .true.
   use_mxmat = .false.

   ! (h0) Save default variables monthly instead of daily to save space
   hist_nhtfrq = 0
   hist_mfilt = 12

   ! (h1) Annual outputs for GDD generation
   hist_fincl2 = 'GRAINC_TO_FOOD_PERHARV', 'GRAINC_TO_FOOD_ANN', 'SDATES', 'SDATES_PERHARV', 'SYEARS_PERHARV', 'HDATES', 'GDDHARV_PERHARV', 'GDDACCUM_PERHARV', 'HUI_PERHARV', 'SOWING_REASON_PERHARV', 'HARVEST_REASON_PERHARV'
   hist_nhtfrq(2) = 17520
   hist_mfilt(2) = 999
   hist_type1d_pertape(2) = 'PFTS'
   hist_dov2xy(2) = .false.
   
   ! (h2) Daily outputs for GDD generation
   hist_fincl3 = 'GDDACCUM', 'GDDHARV'
   hist_nhtfrq(3) = -24
   hist_mfilt(3) = 365
   hist_type1d_pertape(3) = 'PFTS'
   hist_dov2xy(3) = .false.

Once the GDD-generating run completes, calling the following Python script will generate a file that can serve as ``stream_fldFileName_cultivar_gdds`` in subsequent runs (along with maps illustrating the results):
::

   python3 python/ctsm/crop_calendars/generate_gdds.py \
   --input-dir /path/to/output/dir/from/gdd-generating-run \
   --first-season 1980 \
   --last-season 2009 \
   --sdates-file '/path/to/sowing_date_file.nc' \
   --hdates-file '/path/to/harvest_date_file.nc' \
   --output-dir '/path/where/you/want/results/saved' \
   --skip-crops miscanthus,irrigated_miscanthus

The entire process can be illustrated with the RXCROPMATURITY system test. E.g.:

::

   run_sys_tests -t RXCROPMATURITY_Lm61.f10_f10_mg37.IHistClm60BgcCrop.cheyenne_intel.clm-cropMonthOutput --skip-generate --skip-compare
