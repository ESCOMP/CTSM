  ! ========================================================================
  ! We want these to have separate instances for each water isotope (and any other
  ! tracers), in addition to the bulk instance.
  !
  ! For the first step, Matt is doing this by moving these variables to a new type
  ! (though in the end we may end up with the name WaterStateType for this new type).
  ! ========================================================================
  
  ! ------------------------------------------------------------------------
  ! These had separate water isotope instances in the old, work-in-progress water isotope
  ! branch, and so we're pretty sure need separate instances.
  ! ------------------------------------------------------------------------
  
  real(r8), pointer :: snowice_col            (:)   ! col average snow ice lens
  real(r8), pointer :: snowliq_col            (:)   ! col average snow liquid water
  real(r8), pointer :: h2osno_col             (:)   ! col snow water (mm H2O)
  real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
  real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
  real(r8), pointer :: h2osoi_liqice_10cm_col (:)   ! col liquid water + ice lens in top 10cm of soil (kg/m2)
  real(r8), pointer :: h2osoi_vol_col         (:,:) ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
  real(r8), pointer :: h2ocan_patch           (:)   ! patch canopy water (mm H2O)

  real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)
  real(r8), pointer :: qg_col                 (:)   ! col ground specific humidity [kg/kg]

  ! ------------------------------------------------------------------------
  ! These variables are new to this version of the code (i.e., weren't in the clm4-based
  ! water isotope code), but we're pretty sure we need separate instances for each
  ! isotope / tracer.
  ! ------------------------------------------------------------------------

  real(r8), pointer :: h2osfc_col             (:)   ! col surface water (mm H2O)
  real(r8), pointer :: snocan_patch           (:)   ! patch canopy snow water (mm H2O)
  real(r8), pointer :: liqcan_patch           (:)   ! patch canopy liquid water (mm H2O)

  ! The old code had separate instances for wt_col, documented as "total water storage
  ! (unsaturated soil water + groundwater) (mm)". tws_grc looks like the closest analog
  ! in the new code.
  real(r8), pointer :: tws_grc                (:)   ! grc total water storage (mm H2O)

  ! ------------------------------------------------------------------------
  ! I lean towards thinking that these should be moved over, but I'm not sure about them.
  ! These are all variables that have been introduced since clm4.  It would be fine to
  ! either move them now or wait and move them once we know we need to.
  ! ------------------------------------------------------------------------

  ! This looks like the kind of variable that should have separate instances
  real(r8), pointer :: snounload_patch        (:)   ! Canopy snow unloading (mm H2O)

  ! These look similar to qg, which had a separate water isotope instance in the old
  ! branch
  real(r8), pointer :: qg_snow_col            (:)   ! col ground specific humidity [kg/kg]
  real(r8), pointer :: qg_soil_col            (:)   ! col ground specific humidity [kg/kg]
  real(r8), pointer :: qg_h2osfc_col          (:)   ! col ground specific humidity [kg/kg]
  real(r8), pointer :: qaf_lun                (:)   ! lun urban canopy air specific humidity (kg/kg)


  ! ========================================================================
  ! Introduce a new WaterBalanceType. This will also need separate instances for each
  ! isotope / tracer (if we want to do balance checks and related adjustments on the
  ! isotopes/tracers), but it seems useful to separate out these variables that are
  ! just needed for the sake of balance checks.
  ! ========================================================================

  real(r8), pointer :: h2osno_old_col         (:)   ! col snow mass for previous time step (kg/m2) (new)
  real(r8), pointer :: liq1_grc               (:)   ! grc initial gridcell total h2o liq content
  real(r8), pointer :: liq2_grc               (:)   ! grc post land cover change total liq content
  real(r8), pointer :: ice1_grc               (:)   ! grc initial gridcell total h2o ice content
  real(r8), pointer :: ice2_grc               (:)   ! grc post land cover change total ice content

  real(r8), pointer :: begwb_col              (:)   ! water mass begining of the time step
  real(r8), pointer :: endwb_col              (:)   ! water mass end of the time step
  real(r8), pointer :: errh2o_patch           (:)   ! water conservation error (mm H2O)
  real(r8), pointer :: errh2o_col             (:)   ! water conservation error (mm H2O)
  real(r8), pointer :: errh2osno_col          (:)   ! snow water conservation error(mm H2O)

  ! ========================================================================
  ! Other variables in WaterStateType
  !
  ! For some of these, we're pretty sure we don't need separate instances. For others,
  ! we're unsure. Some of these are summary variables (e.g., sum over column) that are
  ! probably mostly/entirely for diagnostic purposes.
  !
  ! If someone knows we'll need separate instances for each water isotope / tracer,
  ! we'll move them over now; otherwise, we'll wait to move them until it becomes
  ! clearly necessary.
  !
  ! We may rename this something like WaterAuxiliaryType (or just WaterAuxType).
  ! ========================================================================

  real(r8), pointer :: snow_depth_col         (:)   ! col snow height of snow covered area (m)
  real(r8), pointer :: snow_persistence_col   (:)   ! col length of time that ground has had non-zero snow thickness (sec)
  real(r8), pointer :: snowdp_col             (:)   ! col area-averaged snow height (m)
  real(r8), pointer :: int_snow_col           (:)   ! col integrated snowfall (mm H2O)
  real(r8), pointer :: snow_layer_unity_col   (:,:) ! value 1 for each snow layer, used for history diagnostics
  real(r8), pointer :: bw_col                 (:,:) ! col partial density of water in the snow pack (ice + liquid) [kg/m3] 

  real(r8), pointer :: h2osoi_liq_tot_col     (:)   ! vertically summed col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
  real(r8), pointer :: h2osoi_ice_tot_col     (:)   ! vertically summed col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
  real(r8), pointer :: air_vol_col            (:,:) ! col air filled porosity
  real(r8), pointer :: h2osoi_liqvol_col      (:,:) ! col volumetric liquid water content (v/v)
  real(r8), pointer :: swe_old_col            (:,:) ! col initial snow water

  real(r8), pointer :: total_plant_stored_h2o_col(:) ! col water that is bound in plants, including roots, sapwood, leaves, etc
  ! in most cases, the vegetation scheme does not have a dynamic
  ! water storage in plants, and thus 0.0 is a suitable for the trivial case.
  ! When FATES is coupled in with plant hydraulics turned on, this storage
  ! term is set to non-zero. (kg/m2 H2O)

  real(r8), pointer :: snw_rds_col            (:,:) ! col snow grain radius (col,lyr)    [m^-6, microns]
  real(r8), pointer :: snw_rds_top_col        (:)   ! col snow grain radius (top layer)  [m^-6, microns]
  real(r8), pointer :: h2osno_top_col         (:)   ! col top-layer mass of snow  [kg]
  real(r8), pointer :: sno_liq_top_col        (:)   ! col snow liquid water fraction (mass), top layer  [fraction]

  real(r8), pointer :: rh_ref2m_patch         (:)   ! patch 2 m height surface relative humidity (%)
  real(r8), pointer :: rh_ref2m_r_patch       (:)   ! patch 2 m height surface relative humidity - rural (%)
  real(r8), pointer :: rh_ref2m_u_patch       (:)   ! patch 2 m height surface relative humidity - urban (%)
  real(r8), pointer :: rh_af_patch            (:)   ! patch fractional humidity of canopy air (dimensionless) ! private
  real(r8), pointer :: rh10_af_patch          (:)   ! 10-day mean patch fractional humidity of canopy air (dimensionless)
  real(r8), pointer :: dqgdT_col              (:)   ! col d(qg)/dT

  ! Fractions
  real(r8), pointer :: frac_sno_col           (:)   ! col fraction of ground covered by snow (0 to 1)
  real(r8), pointer :: frac_sno_eff_col       (:)   ! col fraction of ground covered by snow (0 to 1)
  real(r8), pointer :: frac_iceold_col        (:,:) ! col fraction of ice relative to the tot water (new) (-nlevsno+1:nlevgrnd) 
  real(r8), pointer :: frac_h2osfc_col        (:)   ! col fractional area with surface water greater than zero
  real(r8), pointer :: frac_h2osfc_nosnow_col (:)   ! col fractional area with surface water greater than zero (if no snow present)
  real(r8), pointer :: wf_col                 (:)   ! col soil water as frac. of whc for top 0.05 m (0-1) 
  real(r8), pointer :: wf2_col                (:)   ! col soil water as frac. of whc for top 0.17 m (0-1) 
  real(r8), pointer :: fwet_patch             (:)   ! patch canopy fraction that is wet (0 to 1)
  real(r8), pointer :: fcansno_patch          (:)   ! patch canopy fraction that is snow covered (0 to 1)
  real(r8), pointer :: fdry_patch             (:)   ! patch canopy fraction of foliage that is green and dry [-] (new)

  ! ========================================================================
  ! Variables that had separate wtr_ instances in the in-progress water isotope
  ! branch, but are missing from the new code
  ! ========================================================================

  real(r8), pointer :: wtr_h2ocan_loss_col        (:,:)   ! wiso mass balance correction term for dynamic weights
