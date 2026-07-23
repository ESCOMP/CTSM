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

Here we give the default namelist as it would be created for an "I1850Clm60BgcCrop" compset at 0.9x1.25 resolution with a t232 land-mask on derecho. To edit the namelist you would edit the ``user_nl_clm`` user namelist with just the items you want to change. For simplicity we remove the namelist groups that are empty or not relevant to this compset. In the sections below, for simplicity we will just show the user namelist (``user_nl_clm``) that will add (or modify existing) namelist items to the namelist.

Example 1-2. Default CLM Namelist
---------------------------------
.. I copied the following file from the baseline for ctsm5.4.047 for the SMS_D_Ld5_PS.f09_f09_mt232.I1850Clm60BgcCrop.derecho_gnu.clm-f09_ObscureStreamOpts test.
.. As said above I removed the empty namelist groups and namelist groups that aren't relevant for this case. And then added leading spaces.

::

     &clm_inparm
      albice = 0.50,0.30
      co2_ppmv = 284.7
      co2_type = 'constant'
      collapse_urban = .false.
      compname = 'clm2'
      convert_ocean_to_land = .true.
      create_crop_landunit = .true.
      crop_residue_removal_frac = 0.5d00
      do_sno_oc = .false.
      downscale_hillslope_meteorology = .true.
      finidat = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/initdata_esmf/ctsm5.4/clmi.f09_interp_from.ctsm5.4.CMIP7_ciso_ctsm5.3.075_f09_124_pSASU.clm2.r.0161_c251118.nc'
      flush_gdd20 = .false.
      for_testing_no_crop_seed_replenishment = .false.
      for_testing_run_ncdiopio_tests = .false.
      for_testing_use_repr_structure_pool = .false.
      for_testing_use_second_grain_pool = .false.
      fsnowaging = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc'
      fsnowoptics = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c013122.nc'
      fsurdat = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/surfdata_esmf/ctsm5.4.0/surfdata_0.9x1.25_hist_1850_78pfts_c251022.nc'
      glc_do_dynglacier = .false.
      glc_snow_persistence_max_days = 0
      h2osno_max = 10000.0
      hillslope_fsat_equals_zero = .false.
      hist_dov2xy = .true.,.false.
      hist_fields_list_file = .false.
      hist_fincl1 = 'RC13_CANAIR','RC14_CANAIR','HDM','LNFM','TBUILD_MAX','P_AC','EXCESS_ICE','LND_MBL'
     hist_fincl2 = 'TG', 'TBOT', 'FIRE', 'FIRA', 'FLDS', 'FSDS', 'FSR', 'FSA', 'FGEV', 'FSH', 'FGR',
             'TSOI', 'ERRSOI', 'SABV', 'SABG', 'FSDSVD', 'FSDSND', 'FSDSVI', 'FSDSNI', 'FSRVD', 'FSRND', 'FSRVI',
             'FSRNI', 'TSA', 'FCTR', 'FCEV', 'QBOT', 'RH2M', 'H2OSOI', 'H2OSNO', 'SOILLIQ', 'SOILICE', 'TSA_U',
             'TSA_R', 'TREFMNAV_U', 'TREFMNAV_R', 'TREFMXAV_U', 'TREFMXAV_R', 'TG_U', 'TG_R', 'RH2M_U', 'RH2M_R', 'QRUNOFF_U', 'QRUNOFF_R',
             'SoilAlpha_U', 'SWup', 'LWup', 'URBAN_AC', 'URBAN_HEAT'
     hist_mfilt = 1,1
     hist_ndens = 1,1,1,1,1,1
     hist_nhtfrq = -24,-8
     hist_wrt_matrixcn_diag = .false.
     irrigate = .false.
     maxpatch_glc = 10
     n_dom_landunits = 0
     n_dom_pfts = 0
     nlevsno = 12
     nsegspc = 35
     o3_veg_stress_method = 'unset'
     paramfile = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/paramdata/ctsm60_params.c260518.nc'
     run_zero_weight_urban = .false.
     snicar_dust_optics = 'sahara'
     snicar_numrad_snw = 5
     snicar_snobc_intmix = .true.
     snicar_snodst_intmix = .false.
     snicar_snw_shape = 'hexagonal_plate'
     snicar_solarspec = 'mid_latitude_winter'
     snicar_use_aerosol = .true.
     snow_cover_fraction_method = 'SwensonLawrence2012'
     snow_thermal_cond_glc_method = 'Sturm1997'
     snow_thermal_cond_lake_method = 'Sturm1997'
     snow_thermal_cond_method = 'Sturm1997'
     soil_layerstruct_predefined = '20SL_8.5m'
     spinup_state = 0
     suplnitro = 'NONE'
     toosmall_crop = 0.d00
     toosmall_glacier = 0.d00
     toosmall_lake = 0.d00
     toosmall_soil = 0.d00
     toosmall_urban = 0.d00
     toosmall_wetland = 0.d00
     use_bedrock = .true.
     use_c13 = .true.
     use_c13_timeseries = .true.
     use_c14 = .true.
     use_c14_bombspike = .true.
     use_cn = .true.
     use_crop = .true.
     use_excess_ice = .true.
     use_fates = .false.
     use_fertilizer = .true.
     use_flexiblecn = .true.
     use_fun = .true.
     use_grainproduct = .true.
     use_hillslope = .false.
     use_hillslope_routing = .false.
     use_hydrstress = .true.
     use_lai_streams = .false.
     use_lch4 = .true.
     use_luna = .true.
     use_matrixcn = .false.
     use_nguardrail = .true.
     use_nitrif_denitrif = .true.
     use_nvmovement = .false.
     use_snicar_frc = .false.
     use_soil_matrixcn = .false.
     use_soil_moisture_streams = .false.
     use_ssre = .true.
     use_subgrid_fluxes = .true.
     use_z0m_snowmelt = .true.
     z0param_method = 'Meier2022'
    /
    &ndepdyn_nml
     ndep_taxmode = 'cycle'
     ndep_tintalgo = 'lower'
     ndep_varlist = 'NDEP_month'
     ndepmapalgo = 'redist'
     stream_fldfilename_ndep = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/ndepdata/fndep_clm_WACCM6_CMIP6piControl001_y21-50avg_1850monthly_0.95x1.25_c180802.nc'
     stream_meshfile_ndep = '/glade/campaign/cesm/cesmdata/inputdata/share/meshes/fv0.9x1.25_141008_polemod_ESMFmesh.nc'
     stream_year_first_ndep = 1850
     stream_year_last_ndep = 1850
    /
    &popd_streams
     popdens_tintalgo = 'linear'
     popdensmapalgo = 'consd'
     stream_fldfilename_popdens = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/firedata/clmforc.Li_2025_CMIP7_hdm_0.5x0.5_simyr1850-2025_c251013.nc'
     stream_meshfile_popdens = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/firedata/clmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0_ESMFmesh_cdf5_100621.nc'
     stream_year_first_popdens = 1850
     stream_year_last_popdens = 1850
    /
    &urbantv_streams
     stream_fldfilename_urbantv = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/urbandata/CTSM52_urbantv_Li_2024_0.9x1.25_simyr1849-2106_c20260217.nc'
     stream_meshfile_urbantv = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/urbandata/CLM50_tbuildmax_Oleson_2016_0.9x1_ESMFmesh_cdf5_100621.nc'
     stream_year_first_urbantv = 1850
     stream_year_last_urbantv = 1850
     urbantv_tintalgo = 'upper'
     urbantvmapalgo = 'redist'
    /
    &light_streams
     lightng_tintalgo = 'nearest'
     lightngmapalgo = 'consf'
     stream_fldfilename_lightng = '/glade/campaign/cesm/cesmdata/inputdata/atm/datm7/NASA_LIS/clmforc.Li_2016_climo1995-2013.360x720.lnfm_Total_c160825.nc'
     stream_meshfile_lightng = '/glade/campaign/cesm/cesmdata/inputdata/atm/datm7/NASA_LIS/clmforc.Li_2016_climo1995-2013.360x720_ESMFmesh_cdf5_150621.nc'
     stream_year_first_lightng = 0001
     stream_year_last_lightng = 0001
    /
    &atm2lnd_inparm
     glcmec_downscale_longwave = .false.
     lapse_rate = 0.006
     repartition_rain_snow = .true.
    /
    &lnd2atm_inparm
     melt_non_icesheet_ice_runoff = .true.
    /
    &clm_canopyhydrology_inparm
     use_clm5_fpi = .true.
    /
    &cnphenology
     generate_crop_gdds = .false.
     initial_seed_at_planting = 3.d00
     min_critical_dayl_method = 'DependsOnLat'
     onset_thresh_depends_on_veg = .true.
     use_mxmat = .true.
    /
    &cropcal_streams
     cropcals_rx = .false.
     cropcals_rx_adapt = .true.
     model_year_align_cropcal_cultivar_gdds = 2000
     model_year_align_cropcal_swindows = 2000
     stream_fldfilename_cultivar_gdds = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/cropdata/calendars/processed/gdds_20230829_161011.tweaked_latlons.no_nan_fill.nc'
     stream_fldfilename_gdd20_baseline = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/cropdata/calendars/processed/20230714_cropcals_pr2_1deg.actually2deg.1980-2009.from_GDDB20.interpd_halfdeg.tweaked_latlons.no_nan_fill.nc'
     stream_fldfilename_swindow_end = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/cropdata/calendars/processed/swindow_ends_ggcmi_crop_calendar_phase3_v1.01.2000-2000.20231005_145103.tweaked_latlons.no_nan_fill.nc'
     stream_fldfilename_swindow_start = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/cropdata/calendars/processed/swindow_starts_ggcmi_crop_calendar_phase3_v1.01.2000-2000.20231005_145103.tweaked_latlons.no_nan_fill.nc'
     stream_gdd20_seasons = .false.
     stream_meshfile_cropcal = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/cropdata/calendars/processed/360x720_120830_ESMFmesh_c20210507_cdf5.tweaked_latlons.no_nan_fill.nc'
     stream_year_first_cropcal_cultivar_gdds = 2000
     stream_year_first_cropcal_swindows = 2000
     stream_year_last_cropcal_cultivar_gdds = 2000
     stream_year_last_cropcal_swindows = 2000
    /
    &cnvegcarbonstate
     initial_vegc = 100.d00
    /
    &friction_velocity
     zetamaxstable = 2.0d00
    /
    &mineral_nitrogen_dynamics
     freelivfix_intercept = 0.0117d00
     freelivfix_slope_wet = 0.0006d00
    /
    &soilwater_movement_inparm
     dtmin = 60.
     expensive = 42
     flux_calculation = 1
     inexpensive = 1
     lower_boundary_condition = 2
     soilwater_movement_method = 1
     upper_boundary_condition = 1
     verysmall = 1.e-8
     xtolerlower = 1.e-2
     xtolerupper = 1.e-1
    /
    &rooting_profile_inparm
     rooting_profile_method_carbon = 1
     rooting_profile_method_water = 1
    /
    &soil_resis_inparm
     soil_resis_method = 1
    /
    &bgc_shared
     constrain_stress_deciduous_onset = .true.
    /
    &canopyfluxes_inparm
     itmax_canopy_fluxes = 40
     use_biomass_heat_storage = .true.
     use_undercanopy_stability = .false.
    /
    &clmu_inparm
     building_temp_method = 1
     urban_explicit_ac = .true.
     urban_hac = 'ON_WASTEHEAT'
     urban_traffic = .false.
    /
    &clm_soilstate_inparm
     organic_frac_squared = .false.
    /
    &clm_nitrogen
     carbon_resp_opt = 0
     cn_evergreen_phenology_opt = 1
     cnratio_floating = .true.
     lnc_opt = .true.
     mm_nuptake_opt = .true.
     reduce_dayl_factor = .false.
     vcmax_opt = 3
    /
    &clm_snowhydrology_inparm
     lotmp_snowdensity_method = 'Slater2017'
     reset_snow = .false.
     reset_snow_glc = .false.
     reset_snow_glc_ela = 1.e9
     snow_dzmax_l_1 = 0.03d00
     snow_dzmax_l_2 = 0.07d00
     snow_dzmax_u_1 = 0.02d00
     snow_dzmax_u_2 = 0.05d00
     snow_dzmin_1 = 0.010d00
     snow_dzmin_2 = 0.015d00
     snow_overburden_compaction_method = 'Vionnet2012'
     wind_dependent_snow_density = .true.
    /
    &cnprecision_inparm
     cnegcrit = -6.d+1
     ncrit = 1.d-9
     nnegcrit = -6.d+0
    /
    &clm_glacier_behavior
     glacier_region_behavior = 'single_at_atm_topo','UNSET','virtual','multiple'
     glacier_region_ice_runoff_behavior = 'melted','UNSET','remains_ice','remains_ice'
     glacier_region_melt_behavior = 'remains_in_place','UNSET','replaced_by_ice','replaced_by_ice'
    /
    &crop_inparm
     baset_latvary_intercept = 12.0d00
     baset_latvary_slope = 0.4d00
     baset_mapping = 'varytropicsbylat'
    /
    &surfacealbedo_inparm
     snowveg_affects_radiation = .true.
    /
    &tillage_inparm
     tillage_mode = 'low'
    /
    &ch4par_in
     finundation_method = 'TWS_inversion'
     use_aereoxid_prog = .true.
    /
    &clm_humanindex_inparm
     calc_human_stress_indices = 'FAST'
    /
    &cnmresp_inparm
     br_root = 0.83d-06
    /
    &cnfun_inparm
     nfix_method = 'Bytnerowicz'
    /
    &photosyns_inparm
     leafresp_method = 2
     light_inhibit = .true.
     modifyphoto_and_lmr_forcrop = .true.
     rootstem_acc = .false.
     stomatalcond_method = 'Medlyn2011'
    /
    &cnfire_inparm
     fire_method = 'li2024crujra'
    /
    &cn_general
     dribble_crophrv_xsmrpool_2atm = .true.
    /
    &lifire_inparm
     boreal_peatfire_c = 0.58d-4
     borpeat_fire_soilmoist_denom =  0.3d00
     bt_max = 0.98d00
     bt_min = 0.85d00
     cli_scale = 0.03d00
     cmb_cmplt_fact_cwd = 0.28d00
     cmb_cmplt_fact_litter = 0.5d00
     cropfire_a1 = 0.34d00
     defo_fire_precip_thresh_bdt = 0.6d00
     defo_fire_precip_thresh_bet = 3.0d00
     lfuel = 75.d00
     max_rh30_affecting_fuel = 95.
     non_boreal_peatfire_c = 0.75d-4
     nonborpeat_fire_precip_denom = 6.5d00
     occur_hi_gdp_tree = 0.33d00
     pot_hmn_ign_counts_alpha = 0.010d00
     rh_hgh = 85.0d00
     rh_low = 30.0d00
     ufuel = 825.d00
    /
    &ch4finundated
     ch4finundatedmapalgo = 'redist'
     stream_fldfilename_ch4finundated = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/paramdata/finundated_inversiondata_0.9x1.25_c170706.nc'
     stream_meshfile_ch4finundated = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/paramdata/finundated_inversiondata_0.9x1_ESMFmesh_cdf5_130621.nc'
    /
    &exice_streams
     stream_fldfilename_exice = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/paramdata/exice_init_0.125x0.125_c20220516.nc'
     stream_mapalgo_exice = 'nn'
     stream_meshfile_exice = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/paramdata/exice_init_0.125x0.125_ESMFmesh_cdf5_c20220802.nc'
     use_excess_ice_streams = .true.
    /
    &clm_temperature_inparm
     excess_ice_coldstart_depth = 0.5
     excess_ice_coldstart_temp = -3.15
    /
    &soilbgc_decomp
     soil_decomp_method = 'CENTURYKoven2013'
    /
    &clm_canopy_inparm
     leaf_mr_vcm = 0.015d00
    /
    &prigentroughness
     prigentroughnessmapalgo = 'consf'
     stream_fldfilename_prigentroughness = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/dustemisdata/Prigent_2005_roughness_0.25x0.25_cdf5_c260218.nc'
     stream_meshfile_prigentroughness = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/dustemisdata/dust_0.25x0.25_ESMFmesh_cdf5_c240222.nc'
     use_prigent_roughness = .true.
    /
    &carbon_isotope_streams
     stream_fldfilename_atm_c13 = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/isotopes/ctsmforc.Graven.atm_delta_C13_CMIP7_global_1700-2023_yearly_v3.0_c251013.nc'
     stream_fldfilename_atm_c14 = '/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/isotopes/ctsmforc.Graven.atm_delta_C14_CMIP7_360x720_1700-2023_yearly_v3.0_tweaked_latlons_c260108.no_nan_fill.nc'
     stream_meshfile_atm_c14 = '/glade/campaign/cesm/cesmdata/inputdata/share/meshes/360x720_120830_ESMFmesh_tweaked_latlons_c20260108.nc'
     stream_year_first_atm_c13 = 1850
     stream_year_first_atm_c14 = 1850
     stream_year_last_atm_c13 = 1850
     stream_year_last_atm_c14 = 1850
    /
    &scf_swenson_lawrence_2012_inparm
     int_snow_max = 2000.
     n_melt_glcmec = 1.0d00
    /

Adding/removing fields on your primary history file
---------------------------------------------------

The primary history files are output monthly and contain an extensive list of fieldnames, but the list of fieldnames can be added to using ``hist_fincl1`` or removed from by adding fieldnames to ``hist_fexcl1``. For maximum output, ``hist_all_fields`` will enable all fieldnames on the primary history file.

A sample user namelist ``user_nl_clm`` adding few new fields (cosine of solar zenith angle, and solar declination) and excluding a few standard fields is (ground temperature, vegetation temperature, soil temperature and soil water).:

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

The default averaging depends on the specific fields, but for most fields is an average. A sample user namelist ``user_nl_clm`` follows making the monthly output fields all averages (except ``TSOI``), and adding auxiliary file streams for instantaneous (6-hourly), maximum (daily), minimum (daily), and average (daily). For some of the fields we diverge from the per-tape value given and customize to some different type of averaging.

.. note:: As of ctsm5.4, history files that used to be labeled with hX (where X is an integer from 0 to 4 in the example) will be labeled with hXi and hXa (as separate files) to indicate instantaneous versus non-instantaneous (average, etc.) files. The change intends to prevent confusion associated with the time corresponding to instantaneous history fields by now putting them on separate files than non-instantaneous fields. The separate instantaneous history files represent the exact time step when they were written and do not include a time_bounds variable. Conversely, non-instantaneous history files represent the period of their time_bounds variable. As a result, time data on non-instantaneous history files are now read correctly during post processing (e.g. by xarray). Special handling may still be needed for instantaneous history files, whose timestamps represent the date and time at the END of the history timestep. So, e.g., an instantaneous variable saved at the end of year 2023 will get the timestamp 2024-01-01 00:00:00.

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
   hist_avgflag_pertape = 'A', 'I', 'X', 'M', 'A'
   hist_nhtfrq = 0, -6, -24, -24, -24

In the example we put the same list of fields on each of the tapes: soil-temperature, ground temperature, vegetation temperature, emitted longwave radiation, reflected solar radiation, sensible heat, total latent-heat, and total water storage. We also modify the soil temperature for the primary and secondary auxiliary tapes by outputting them for a maximum instead of the prescribed per-tape of average and instantaneous respectively. For the tertiary auxiliary tape we output ground temperature as instantaneous instead of as maximum, and for the fourth auxiliary tape we output vegetation temperature as instantaneous instead of as minimum. Finally, for the fifth auxiliary tapes we output ``FIRE`` as instantaneous instead of as average.

.. note:: We also use ``hist_empty_htapes`` as in the previous example, so we can list ONLY the fields that we want on the primary history tape.

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
   hist_type1d_pertape = ' ', 'GRID', 'COLS', ' '
   hist_nhtfrq = 0, -24, -24, -24

.. warning:: ``LAND`` and ``COLS`` are also options to the pertape averaging, but currently there is a bug with them and they fail to work.

.. note:: Technically the default for ``hist_nhtfrq`` is for primary files output monthly and the other auxiliary tapes for daily, so we don't actually have to include ``hist_nhtfrq``, we could use the default for it. Here we specify it for clarity.

Visualizing global 1D vector files will take effort. You'll probably want to do some post-processing and possibly just extract out single points of interest to see what is going on. Since the output is a 1D vector of only land points, traditional plots won't be helpful. The number of points per grid-cell will also vary for anything but grid-cell averaging. You'll need to use the output fields ``pfts1d_ixy``, and ``pfts1d_jxy``, to get the mapping of the fields to the global 2D array. ``pfts1d_itype_veg`` gives you the PFT number for each PFT. Most likely you'll want to do this analysis in a data processing tool (such as NCL, Matlab, Mathmatica, IDL, etc. that is able to read and process NetCDF data files).
