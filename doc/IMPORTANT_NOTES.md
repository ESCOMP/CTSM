# Important Notes on Experimental Features of CTSM

Namelist items that are not regularly tested or used. Some aren't even implemented.

    See

    '../bld/namelist_files/namelist_definition_ctsm.xml' -- for definitions of all namelist variables

## CTSM experimental namelist items

    The following are tested but not on by default (for any physics)

    - all_active
    - allow_invalid_gdd20_season_inputs
    - h2osfcflag (deprecated)
    - use_nvmovement
    - use_soil_moisture_streams

   The following are NOT currently tested nor turned on by default:

    - allowlakeprod
    - allow_invalid_swindow_inputs
    - carbon_resp_opt
    - ch4offline
    - CN_evergreen_phenology_opt
    - CNratio_floating
    - do_sno_oc
    - finundation_method = h2osfc
    - maxpatch_glc /= 10
    - megan_use_gamma_sm
    - no_frozen_nitrif_denitrif
    - overburden_compress_tfactor
    - override_bgc_restart_mismatch_dump
    - perchroot
    - perchroot_alt
    - pertlim (deprecated)
    - reduce_dayl_factor (not implemented, it's commented out in the code)
    - replenishlakec
    - rooting_profile_method_soilcarbon (DELETE)
    - rooting_profile_varindex_carbon /= 2
    - rooting_profile_varindex_water /= 1
    - snicar_dust_optics /= sahara
    - snicar_numrad_snw /= 5
    - snicar_snobc_intmix /= TRUE
    - snicar_snodst_intmix /= TRUE
    - snicar_snw_shape /= hexagonal_plate
    - snicar_solarspec /= mid_latitude_winter
    - snicar_use_aerosol /= FALSE
    - urban_traffic (not implemented)
    - use_cndv (deprecated)
    - use_extralakelayers
    - usefrootc
    - usephfact
    - use_vichydro (deprecated)
    - vcmax_opt = 4

## FATES experimental namelist items

    FATES is a relatively new subcomponent of CTSM
    Almost all FATES options include "fates" in the name

    The following are tested, but not turned on by default:

   - fates_seeddisp_cadence > 0
   - fates_parteh_mode > 1
   - use_fates_planthydro
   - use_fates_managed_fire
   - use_fates_tree_damage
   - use_fates_cohort_age_tracking
   - use_fates_lupft
   - use_fates_potentialveg
   - use_fates_ed_st3

   The following are NOT currently tested nor turned on by default:

   - fates_spitfire_mode == 2
   - fates_spitfire_mode == 5
   - fates_hydro_solver /= 2DPicard
   - fates_photosynth_acclimation == kumarathunge2019
   - fates_cstarvation_model == exponential
   - fates_stomatal_assimilation = = gross
   - fates_stomatal_model ==  medlyn2011
   - fates_leafresp_model ==  atkin2017
   - use_fates_potentialveg
   - use_fates_daylength_factor == FALSE
   - fates_history_dimlevel == 0

