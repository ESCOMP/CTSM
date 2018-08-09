module WaterFluxBulkType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water fluxes that just apply to bulk water. Note
  ! that this type extends the base waterflux_type, so the full waterfluxbulk_type
  ! contains the union of the fields defined here and the fields defined in
  ! waterflux_type.
  !
  ! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevsoi
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch   
  use CNSharedParamsMod           , only : use_fun
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  use WaterFluxType     , only : waterflux_type
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, extends(waterflux_type), public :: waterfluxbulk_type

     ! water fluxes are in units or mm/s

     real(r8), pointer :: qflx_irrig_patch         (:)   ! patch irrigation flux (mm H2O/s) [+]
     real(r8), pointer :: qflx_irrig_col           (:)   ! col irrigation flux (mm H2O/s) [+]
     real(r8), pointer :: qflx_phs_neg_col         (:)   ! col sum of negative hydraulic redistribution fluxes (mm H2O/s) [+]

     ! In the snow capping parametrization excess mass above h2osno_max is removed.  A breakdown of mass into liquid 
     real(r8), pointer :: qflx_snwcp_discarded_liq_col(:) ! col excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_discarded_ice_col(:) ! col excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
     real(r8), pointer :: qflx_glcice_dyn_water_flux_col(:) ! col water flux needed for balance check due to glc_dyn_runoff_routing (mm H2O/s) (positive means addition of water to the system); valid for all columns


     real(r8), pointer :: qflx_snowindunload_patch (:)   ! patch canopy snow wind unloading (mm H2O /s)
     real(r8), pointer :: qflx_snowindunload_col   (:)   ! col canopy snow wind unloading (mm H2O /s)
     real(r8), pointer :: qflx_snotempunload_patch (:)   ! patch canopy snow temp unloading (mm H2O /s) 
     real(r8), pointer :: qflx_snotempunload_col   (:)   ! col canopy snow temp unloading (mm H2O /s) 

     real(r8), pointer :: qflx_ev_snow_patch       (:)   ! patch evaporation heat flux from snow       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_snow_col         (:)   ! col evaporation heat flux from snow         (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_soil_patch       (:)   ! patch evaporation heat flux from soil       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_soil_col         (:)   ! col evaporation heat flux from soil         (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_h2osfc_patch     (:)   ! patch evaporation heat flux from soil       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_h2osfc_col       (:)   ! col evaporation heat flux from soil         (mm H2O/s) [+ to atm]

     real(r8), pointer :: qflx_adv_col             (:,:) ! col advective flux across different soil layer interfaces [mm H2O/s] [+ downward]
     real(r8), pointer :: qflx_rootsoi_col         (:,:) ! col root and soil water exchange [mm H2O/s] [+ into root]
     real(r8), pointer :: qflx_sat_excess_surf_col (:)   ! col surface runoff due to saturated surface (mm H2O /s)
     real(r8), pointer :: qflx_infl_excess_col     (:)   ! col infiltration excess runoff (mm H2O /s)
     real(r8), pointer :: qflx_infl_excess_surf_col(:)   ! col surface runoff due to infiltration excess (mm H2O /s)
     real(r8), pointer :: qflx_h2osfc_surf_col     (:)   ! col surface water runoff (mm H2O /s)
     real(r8), pointer :: qflx_rain_plus_snomelt_col(:)  ! col rain plus snow melt falling on the soil (mm/s)
     real(r8), pointer :: qflx_in_soil_col         (:)   ! col surface input to soil (mm/s)
     real(r8), pointer :: qflx_in_soil_limited_col (:)   ! col surface input to soil, limited by max infiltration rate (mm/s)
     real(r8), pointer :: qflx_h2osfc_drain_col    (:)   ! col bottom drainage from h2osfc (mm/s)
     real(r8), pointer :: qflx_top_soil_to_h2osfc_col(:) ! col portion of qflx_top_soil going to h2osfc, minus evaporation (mm/s)
     real(r8), pointer :: qflx_in_h2osfc_col(:)          ! col total surface input to h2osfc
     real(r8), pointer :: qflx_h2osfc_to_ice_col   (:)   ! col conversion of h2osfc to ice
     real(r8), pointer :: qflx_snow_h2osfc_col     (:)   ! col snow falling on surface water
     real(r8), pointer :: qflx_drain_perched_col   (:)   ! col sub-surface runoff from perched wt (mm H2O /s)
     real(r8), pointer :: qflx_deficit_col         (:)   ! col water deficit to keep non-negative liquid water content (mm H2O)   
     real(r8), pointer :: qflx_snomelt_lyr_col     (:,:) ! col snow melt in each layer (mm H2O /s)
     real(r8), pointer :: qflx_snow_drain_col      (:)   ! col drainage from snow pack
     real(r8), pointer :: qflx_ice_runoff_snwcp_col(:)   ! col solid runoff from snow capping (mm H2O /s)
     real(r8), pointer :: qflx_ice_runoff_xs_col   (:)   ! col solid runoff from excess ice in soil (mm H2O /s)
     real(r8), pointer :: qflx_drain_vr_col        (:,:) ! col liquid water losted as drainage (m /time step)
     real(r8), pointer :: snow_sources_col         (:)   ! col snow sources (mm H2O/s)
     real(r8), pointer :: snow_sinks_col           (:)   ! col snow sinks (mm H2O/s)


     ! Objects that help convert once-per-year dynamic land cover changes into fluxes
     ! that are dribbled throughout the year
     type(annual_flux_dribbler_type) :: qflx_liq_dynbal_dribbler
     type(annual_flux_dribbler_type) :: qflx_ice_dynbal_dribbler

     ! ET accumulation
     real(r8), pointer :: AnnEt                    (:)   ! Annual average ET flux mmH20/s                                     

   contains
 
     
     
     procedure, public  :: InitBulk
     procedure, public  :: RestartBulk     
     procedure, private :: InitBulkAllocate 
     procedure, private :: InitBulkHistory  
     procedure, private :: InitBulkCold     
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars

  end type waterfluxbulk_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine InitBulk(this, bounds, info)

    class(waterfluxbulk_type) :: this
    type(bounds_type), intent(in)    :: bounds  
    class(water_info_base_type), intent(in), target :: info

    call this%Init(bounds, info)
    call this%InitBulkAllocate(bounds) ! same as "call initAllocate_type(hydro, bounds)"
    call this%InitBulkHistory(bounds)
    call this%InitBulkCold(bounds)

  end subroutine InitBulk

  !------------------------------------------------------------------------
  subroutine InitBulkAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterfluxbulk_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%qflx_irrig_patch         (begp:endp))              ; this%qflx_irrig_patch         (:)   = nan

    allocate(this%qflx_snowindunload_patch (begp:endp))              ; this%qflx_snowindunload_patch (:)   = nan
    allocate(this%qflx_snowindunload_col   (begp:endp))              ; this%qflx_snowindunload_col   (:)   = nan
    allocate(this%qflx_snotempunload_patch (begp:endp))              ; this%qflx_snotempunload_patch (:)   = nan
    allocate(this%qflx_snotempunload_col   (begp:endp))              ; this%qflx_snotempunload_col   (:)   = nan
	

    allocate(this%qflx_irrig_col           (begc:endc))              ; this%qflx_irrig_col           (:)   = nan
    allocate(this%qflx_snwcp_discarded_liq_col(begc:endc))           ; this%qflx_snwcp_discarded_liq_col(:) = nan
    allocate(this%qflx_snwcp_discarded_ice_col(begc:endc))           ; this%qflx_snwcp_discarded_ice_col(:) = nan
    allocate(this%qflx_glcice_dyn_water_flux_col(begc:endc))         ; this%qflx_glcice_dyn_water_flux_col (:) = nan
    allocate(this%qflx_phs_neg_col         (begc:endc))              ; this%qflx_phs_neg_col       (:)   = nan

    allocate( this%qflx_ev_snow_patch      (begp:endp))              ; this%qflx_ev_snow_patch       (:)   = nan
    allocate( this%qflx_ev_snow_col        (begc:endc))              ; this%qflx_ev_snow_col         (:)   = nan
    allocate( this%qflx_ev_soil_patch      (begp:endp))              ; this%qflx_ev_soil_patch       (:)   = nan
    allocate( this%qflx_ev_soil_col        (begc:endc))              ; this%qflx_ev_soil_col         (:)   = nan
    allocate( this%qflx_ev_h2osfc_patch    (begp:endp))              ; this%qflx_ev_h2osfc_patch     (:)   = nan
    allocate( this%qflx_ev_h2osfc_col      (begc:endc))              ; this%qflx_ev_h2osfc_col       (:)   = nan

    allocate(this%qflx_drain_vr_col      (begc:endc,1:nlevsoi))      ; this%qflx_drain_vr_col        (:,:) = nan
    allocate(this%qflx_adv_col             (begc:endc,0:nlevsoi))    ; this%qflx_adv_col             (:,:) = nan
    allocate(this%qflx_rootsoi_col         (begc:endc,1:nlevsoi))    ; this%qflx_rootsoi_col         (:,:) = nan
    allocate(this%qflx_sat_excess_surf_col (begc:endc))              ; this%qflx_sat_excess_surf_col (:)   = nan
    allocate(this%qflx_infl_excess_col     (begc:endc))              ; this%qflx_infl_excess_col     (:)   = nan
    allocate(this%qflx_rain_plus_snomelt_col(begc:endc))             ; this%qflx_rain_plus_snomelt_col(:)  = nan
    allocate(this%qflx_in_soil_col         (begc:endc))              ; this%qflx_in_soil_col         (:)   = nan
    allocate(this%qflx_in_soil_limited_col (begc:endc))              ; this%qflx_in_soil_limited_col (:)   = nan
    allocate(this%qflx_h2osfc_drain_col    (begc:endc))              ; this%qflx_h2osfc_drain_col    (:)   = nan
    allocate(this%qflx_top_soil_to_h2osfc_col(begc:endc))            ; this%qflx_top_soil_to_h2osfc_col(:) = nan
    allocate(this%qflx_in_h2osfc_col       (begc:endc))              ; this%qflx_in_h2osfc_col(:)          = nan
    allocate(this%qflx_h2osfc_to_ice_col   (begc:endc))              ; this%qflx_h2osfc_to_ice_col   (:)   = nan
    allocate(this%qflx_infl_excess_surf_col(begc:endc))              ; this%qflx_infl_excess_surf_col(:)   = nan
    allocate(this%qflx_h2osfc_surf_col     (begc:endc))              ; this%qflx_h2osfc_surf_col     (:)   = nan
    allocate(this%qflx_snow_h2osfc_col     (begc:endc))              ; this%qflx_snow_h2osfc_col     (:)   = nan
    allocate(this%qflx_snomelt_lyr_col     (begc:endc,-nlevsno+1:0)) ; this%qflx_snomelt_lyr_col     (:,:) = nan
    allocate(this%qflx_snow_drain_col      (begc:endc))              ; this%qflx_snow_drain_col      (:)   = nan
    allocate(this%qflx_drain_perched_col   (begc:endc))              ; this%qflx_drain_perched_col   (:)   = nan
    allocate(this%qflx_deficit_col         (begc:endc))              ; this%qflx_deficit_col         (:)   = nan
    allocate(this%qflx_ice_runoff_snwcp_col(begc:endc))              ; this%qflx_ice_runoff_snwcp_col(:)   = nan
    allocate(this%qflx_ice_runoff_xs_col   (begc:endc))              ; this%qflx_ice_runoff_xs_col   (:)   = nan
    allocate(this%snow_sources_col         (begc:endc))              ; this%snow_sources_col         (:)   = nan   
    allocate(this%snow_sinks_col           (begc:endc))              ; this%snow_sinks_col           (:)   = nan   
    allocate(this%AnnET                    (begc:endc))              ; this%AnnET                    (:)   = nan


    this%qflx_liq_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = 'qflx_liq_dynbal', &
         units = 'mm H2O')

    this%qflx_ice_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = 'qflx_ice_dynbal', &
         units = 'mm H2O')

  end subroutine InitBulkAllocate

  !------------------------------------------------------------------------
  subroutine InitBulkHistory(this, bounds)
    !
    ! !USES:
    use clm_varctl  , only : use_cn
    use histFileMod , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(waterfluxbulk_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    this%qflx_snow_drain_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SNOW_DRAIN'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('drainage from snow pack'), &
         ptr_col=this%qflx_snow_drain_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SNOW_DRAIN_ICE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('drainage from snow pack melt (ice landunits only)'), &
         ptr_col=this%qflx_snow_drain_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_snomelt_lyr_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%qflx_snomelt_lyr_col(begc:endc,-nlevsno+1:0)
    call hist_addfld2d ( &
         fname=this%info%fname('SNO_MELT'),  &
         units='mm/s', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('snow melt rate in each snow layer'), &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d ( &
         fname=this%info%fname('SNO_MELT_ICE'),  &
         units='mm/s', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('snow melt rate in each snow layer (ice landunits only)'), &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')


    this%qflx_h2osfc_to_ice_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QH2OSFC_TO_ICE'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('surface water converted to ice'), &
         ptr_col=this%qflx_h2osfc_to_ice_col, default='inactive')

    call hist_addfld2d ( &
         fname=this%info%fname('QROOTSINK'),  &
         units='mm/s', type2d='levsoi', &
         avgflag='A', &
         long_name=this%info%lname('water flux from soil to root in each soil-layer'), &
         ptr_col=this%qflx_rootsoi_col, set_spec=spval, l2g_scale_type='veg', default='inactive')

    this%qflx_ev_snow_patch(begp:endp) = spval

    this%qflx_snowindunload_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNO_WINDUNLOAD'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('canopy snow wind unloading'), &
         ptr_patch=this%qflx_snowindunload_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_snotempunload_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNO_TEMPUNLOAD'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('canopy snow temp unloading'), &
         ptr_patch=this%qflx_snotempunload_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_irrig_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QIRRIG'), &
         units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('water added through irrigation'), &
         ptr_patch=this%qflx_irrig_patch)

    this%qflx_h2osfc_surf_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QH2OSFC'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('surface water runoff'), &
         ptr_col=this%qflx_h2osfc_surf_col)

    this%qflx_drain_perched_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRAI_PERCH'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('perched wt drainage'), &
         ptr_col=this%qflx_drain_perched_col, c2l_scale_type='urbanf')

    this%qflx_phs_neg_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QPHSNEG'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('net negative hydraulic redistribution flux'), &
         ptr_col=this%qflx_phs_neg_col, default='inactive')

    ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at any
    ! given time step but only if there is at least one snow layer (for all landunits 
    ! except lakes).  Also note that monthly average files of snow_sources and snow_sinks
    ! sinks must be weighted by number of days in the month to diagnose, for example, an 
    ! annual value of the change in h2osno. 

    this%snow_sources_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOW_SOURCES'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('snow sources (liquid water)'), &
         ptr_col=this%snow_sources_col, c2l_scale_type='urbanf')

    this%snow_sinks_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOW_SINKS'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('snow sinks (liquid water)'), &
         ptr_col=this%snow_sinks_col, c2l_scale_type='urbanf')

    this%AnnET(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('AnnET'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('Annual ET'), &
         ptr_col=this%AnnET, c2l_scale_type='urbanf', default='inactive')
         
  end subroutine InitBulkHistory
  
  
  


  !-----------------------------------------------------------------------
  subroutine InitBulkCold(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterfluxbulk_type) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l
    !-----------------------------------------------------------------------

    this%qflx_phs_neg_col(bounds%begc:bounds%endc)   = 0.0_r8

    this%qflx_h2osfc_surf_col(bounds%begc:bounds%endc) = 0._r8
    this%qflx_snow_drain_col(bounds%begc:bounds%endc)  = 0._r8

    ! This variable only gets set in the hydrology filter; need to initialize it to 0 for
    ! the sake of columns outside this filter
    this%qflx_ice_runoff_xs_col(bounds%begc:bounds%endc) = 0._r8

    ! Initialize qflx_glcice_dyn_water_flux_col to 0 for all columns because we want this
    ! flux to remain 0 for columns where is is never set, including non-glacier columns.
    !
    ! Other qflx_glcice fluxes intentionally remain unset (spval) outside the do_smb
    ! filter, so that they are flagged as missing value outside that filter.
    this%qflx_glcice_dyn_water_flux_col(bounds%begc:bounds%endc) = 0._r8

    this%AnnEt(bounds%begc:bounds%endc)                 = 0._r8


  end subroutine InitBulkCold

   !-----------------------------------------------------------------------
    subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use clm_varcon  , only : spval
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(waterfluxbulk_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !---------------------------------------------------------------------

    if (use_fun) then
   
       call init_accum_field (name='AnnET', units='MM H2O/S', &
            desc='365-day running mean of total ET', accum_type='runmean', accum_period=-365, &
            subgrid_type='column', numlev=1, init_value=0._r8)

    end if

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
    !
     subroutine InitAccVars (this, bounds)
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(waterfluxbulk_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begc, endc
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    !---------------------------------------------------------------------
    begc = bounds%begc; endc = bounds%endc

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begc:endc), stat=ier)

    ! Determine time step
    nstep = get_nstep()

    if (use_fun) then
       call extract_accum_field ('AnnET', rbufslp, nstep)
       this%qflx_evap_tot_col(begc:endc) = rbufslp(begc:endc)
    end if

    deallocate(rbufslp)

  end subroutine InitAccVars
  
  
  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(waterfluxbulk_type)                 :: this
    type(bounds_type)      , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p                     ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begc, endc
    real(r8), pointer :: rbufslp(:)      ! temporary single level - patch level
    !---------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level patch field

    allocate(rbufslp(begc:endc), stat=ier)
    
    do c = begc,endc
       rbufslp(c) = this%qflx_evap_tot_col(c)
    end do
    if (use_fun) then
       ! Accumulate and extract AnnET (accumulates total ET as 365-day running mean)
       call update_accum_field  ('AnnET', rbufslp, nstep)
       call extract_accum_field ('AnnET', this%AnnET, nstep)
    
    end if

    deallocate(rbufslp)
    
  end subroutine UpdateAccVars

  !------------------------------------------------------------------------
  subroutine RestartBulk(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterfluxbulk_type)            :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call this%restart ( bounds, ncid, flag=flag )
    ! needed for SNICAR

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('qflx_snow_drain')//':'//this%info%fname('qflx_snow_melt'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('drainage from snow column'), &
         units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snow_drain_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snow_drain to zero
       this%qflx_snow_drain_col(bounds%begc:bounds%endc) = 0._r8
    endif
    
    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('AnnET'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('Annual ET'), &
         units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%AnnET)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snow_drain to zero
       this%AnnET(bounds%begc:bounds%endc) = 0._r8
    endif

    call this%qflx_liq_dynbal_dribbler%Restart(bounds, ncid, flag)
    call this%qflx_ice_dynbal_dribbler%Restart(bounds, ncid, flag)

  end subroutine RestartBulk

end module WaterFluxBulkType
