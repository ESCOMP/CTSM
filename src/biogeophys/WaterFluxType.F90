module WaterFluxType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water fluxes that apply to both bulk water and
  ! water tracers.
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
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterflux_type

     class(water_info_base_type), pointer :: info

     ! water fluxes are in units or mm/s

     real(r8), pointer :: qflx_prec_grnd_patch     (:)   ! patch water onto ground including canopy runoff [kg/(m2 s)]
     real(r8), pointer :: qflx_prec_grnd_col       (:)   ! col water onto ground including canopy runoff [kg/(m2 s)]
     real(r8), pointer :: qflx_rain_grnd_patch     (:)   ! patch rain on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_rain_grnd_col       (:)   ! col rain on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_snow_grnd_patch     (:)   ! patch snow on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_snow_grnd_col       (:)   ! col snow on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_sub_snow_patch      (:)   ! patch sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_sub_snow_col        (:)   ! col sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_evap_soi_patch      (:)   ! patch soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_soi_col        (:)   ! col soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_veg_patch      (:)   ! patch vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_veg_col        (:)   ! col vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_can_patch      (:)   ! patch evaporation from leaves and stems (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_can_col        (:)   ! col evaporation from leaves and stems (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_tot_patch      (:)   ! patch pft_qflx_evap_soi + pft_qflx_evap_veg + qflx_tran_veg
     real(r8), pointer :: qflx_evap_tot_col        (:)   ! col col_qflx_evap_soi + col_qflx_evap_veg + qflx_tran_veg
     real(r8), pointer :: qflx_evap_grnd_patch     (:)   ! patch ground surface evaporation rate (mm H2O/s) [+]
     real(r8), pointer :: qflx_evap_grnd_col       (:)   ! col ground surface evaporation rate (mm H2O/s) [+]

     ! In the snow capping parametrization excess mass above h2osno_max is removed.  A breakdown of mass into liquid 
     ! and solid fluxes is done, these are represented by qflx_snwcp_liq_col and qflx_snwcp_ice_col. 
     real(r8), pointer :: qflx_snwcp_liq_col       (:)   ! col excess liquid h2o due to snow capping (outgoing) (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_ice_col       (:)   ! col excess solid h2o due to snow capping (outgoing) (mm H2O /s)
     real(r8), pointer :: qflx_glcice_col(:)              ! col net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC; only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_frz_col (:)         ! col ice growth (positive definite) (mm H2O/s); only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_melt_col(:)         ! col ice melt (positive definite) (mm H2O/s); only valid inside the do_smb_c filter


     real(r8), pointer :: qflx_tran_veg_patch      (:)   ! patch vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_tran_veg_col        (:)   ! col vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_dew_snow_patch      (:)   ! patch surface dew added to snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_snow_col        (:)   ! col surface dew added to snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_patch      (:)   ! patch ground surface dew formation (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_col        (:)   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
     real(r8), pointer :: qflx_prec_intr_patch     (:)   ! patch interception of precipitation [mm/s]
     real(r8), pointer :: qflx_prec_intr_col       (:)   ! col interception of precipitation [mm/s]

     real(r8), pointer :: qflx_infl_col            (:)   ! col infiltration (mm H2O /s)
     real(r8), pointer :: qflx_surf_col            (:)   ! col total surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_drain_col           (:)   ! col sub-surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_top_soil_col        (:)   ! col net water input into soil from top (mm/s)
     real(r8), pointer :: qflx_floodc_col          (:)   ! col flood water flux at column level
     real(r8), pointer :: qflx_sl_top_soil_col     (:)   ! col liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
     real(r8), pointer :: qflx_snomelt_col         (:)   ! col snow melt (mm H2O /s)
     real(r8), pointer :: qflx_qrgwl_col           (:)   ! col qflx_surf at glaciers, wetlands, lakes
     real(r8), pointer :: qflx_runoff_col          (:)   ! col total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
     real(r8), pointer :: qflx_runoff_r_col        (:)   ! col Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
     real(r8), pointer :: qflx_runoff_u_col        (:)   ! col urban total runoff (qflx_drain+qflx_surf) (mm H2O /s) 
     real(r8), pointer :: qflx_rsub_sat_col        (:)   ! col soil saturation excess [mm/s]
     real(r8), pointer :: qflx_snofrz_lyr_col      (:,:) ! col snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
     real(r8), pointer :: qflx_snofrz_col          (:)   ! col column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]

     ! Dynamic land cover change
     real(r8), pointer :: qflx_liq_dynbal_grc      (:)   ! grc liq dynamic land cover change conversion runoff flux
     real(r8), pointer :: qflx_ice_dynbal_grc      (:)   ! grc ice dynamic land cover change conversion runoff flux

   contains
 
     
     
     procedure, public  :: Init
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type waterflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info)

    class(waterflux_type) :: this
    type(bounds_type), intent(in)    :: bounds  
    class(water_info_base_type), intent(in), target :: info

    this%info => info

    call this%InitAllocate(bounds) ! same as "call initAllocate_type(hydro, bounds)"
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
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

    allocate(this%qflx_prec_intr_patch     (begp:endp))              ; this%qflx_prec_intr_patch     (:)   = nan
    allocate(this%qflx_prec_grnd_patch     (begp:endp))              ; this%qflx_prec_grnd_patch     (:)   = nan
    allocate(this%qflx_rain_grnd_patch     (begp:endp))              ; this%qflx_rain_grnd_patch     (:)   = nan
    allocate(this%qflx_snow_grnd_patch     (begp:endp))              ; this%qflx_snow_grnd_patch     (:)   = nan
    allocate(this%qflx_sub_snow_patch      (begp:endp))              ; this%qflx_sub_snow_patch      (:)   = 0.0_r8
    allocate(this%qflx_tran_veg_patch      (begp:endp))              ; this%qflx_tran_veg_patch      (:)   = nan

    allocate(this%qflx_dew_grnd_patch      (begp:endp))              ; this%qflx_dew_grnd_patch      (:)   = nan
    allocate(this%qflx_dew_snow_patch      (begp:endp))              ; this%qflx_dew_snow_patch      (:)   = nan

    allocate(this%qflx_prec_intr_col       (begc:endc))              ; this%qflx_prec_intr_col       (:)   = nan
    allocate(this%qflx_prec_grnd_col       (begc:endc))              ; this%qflx_prec_grnd_col       (:)   = nan
    allocate(this%qflx_rain_grnd_col       (begc:endc))              ; this%qflx_rain_grnd_col       (:)   = nan
    allocate(this%qflx_snow_grnd_col       (begc:endc))              ; this%qflx_snow_grnd_col       (:)   = nan
    allocate(this%qflx_sub_snow_col        (begc:endc))              ; this%qflx_sub_snow_col        (:)   = 0.0_r8
    allocate(this%qflx_snwcp_liq_col       (begc:endc))              ; this%qflx_snwcp_liq_col       (:)   = nan
    allocate(this%qflx_snwcp_ice_col       (begc:endc))              ; this%qflx_snwcp_ice_col       (:)   = nan
    allocate(this%qflx_glcice_col          (begc:endc))              ; this%qflx_glcice_col          (:)   = nan
    allocate(this%qflx_glcice_frz_col      (begc:endc))              ; this%qflx_glcice_frz_col      (:)   = nan
    allocate(this%qflx_glcice_melt_col     (begc:endc))              ; this%qflx_glcice_melt_col     (:)   = nan
    allocate(this%qflx_tran_veg_col        (begc:endc))              ; this%qflx_tran_veg_col        (:)   = nan
    allocate(this%qflx_evap_veg_col        (begc:endc))              ; this%qflx_evap_veg_col        (:)   = nan
    allocate(this%qflx_evap_can_col        (begc:endc))              ; this%qflx_evap_can_col        (:)   = nan
    allocate(this%qflx_evap_soi_col        (begc:endc))              ; this%qflx_evap_soi_col        (:)   = nan
    allocate(this%qflx_evap_tot_col        (begc:endc))              ; this%qflx_evap_tot_col        (:)   = nan
    allocate(this%qflx_evap_grnd_col       (begc:endc))              ; this%qflx_evap_grnd_col       (:)   = nan
    allocate(this%qflx_dew_grnd_col        (begc:endc))              ; this%qflx_dew_grnd_col        (:)   = nan
    allocate(this%qflx_dew_snow_col        (begc:endc))              ; this%qflx_dew_snow_col        (:)   = nan
    allocate(this%qflx_evap_veg_patch      (begp:endp))              ; this%qflx_evap_veg_patch      (:)   = nan
    allocate(this%qflx_evap_can_patch      (begp:endp))              ; this%qflx_evap_can_patch      (:)   = nan
    allocate(this%qflx_evap_soi_patch      (begp:endp))              ; this%qflx_evap_soi_patch      (:)   = nan
    allocate(this%qflx_evap_tot_patch      (begp:endp))              ; this%qflx_evap_tot_patch      (:)   = nan
    allocate(this%qflx_evap_grnd_patch     (begp:endp))              ; this%qflx_evap_grnd_patch     (:)   = nan


    allocate(this%qflx_infl_col            (begc:endc))              ; this%qflx_infl_col            (:)   = nan
    allocate(this%qflx_surf_col            (begc:endc))              ; this%qflx_surf_col            (:)   = nan
    allocate(this%qflx_drain_col           (begc:endc))              ; this%qflx_drain_col           (:)   = nan
    allocate(this%qflx_top_soil_col        (begc:endc))              ; this%qflx_top_soil_col        (:)   = nan
    allocate(this%qflx_snomelt_col         (begc:endc))              ; this%qflx_snomelt_col         (:)   = nan
    allocate(this%qflx_snofrz_col          (begc:endc))              ; this%qflx_snofrz_col          (:)   = nan
    allocate(this%qflx_snofrz_lyr_col      (begc:endc,-nlevsno+1:0)) ; this%qflx_snofrz_lyr_col      (:,:) = nan
    allocate(this%qflx_qrgwl_col           (begc:endc))              ; this%qflx_qrgwl_col           (:)   = nan
    allocate(this%qflx_floodc_col          (begc:endc))              ; this%qflx_floodc_col          (:)   = nan
    allocate(this%qflx_sl_top_soil_col     (begc:endc))              ; this%qflx_sl_top_soil_col     (:)   = nan
    allocate(this%qflx_runoff_col          (begc:endc))              ; this%qflx_runoff_col          (:)   = nan
    allocate(this%qflx_runoff_r_col        (begc:endc))              ; this%qflx_runoff_r_col        (:)   = nan
    allocate(this%qflx_runoff_u_col        (begc:endc))              ; this%qflx_runoff_u_col        (:)   = nan
    allocate(this%qflx_rsub_sat_col        (begc:endc))              ; this%qflx_rsub_sat_col        (:)   = nan

    allocate(this%qflx_liq_dynbal_grc      (begg:endg))              ; this%qflx_liq_dynbal_grc      (:)   = nan
    allocate(this%qflx_ice_dynbal_grc      (begg:endg))              ; this%qflx_ice_dynbal_grc      (:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use clm_varctl  , only : use_cn
    use histFileMod , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
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

    this%qflx_top_soil_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QTOPSOIL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('water input to surface'), &
         ptr_col=this%qflx_top_soil_col, c2l_scale_type='urbanf', default='inactive')

    this%qflx_infl_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QINFL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('infiltration'), &
         ptr_col=this%qflx_infl_col, c2l_scale_type='urbanf')

    this%qflx_surf_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QOVER'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('total surface runoff (includes QH2OSFC)'), &
         ptr_col=this%qflx_surf_col, c2l_scale_type='urbanf')

    this%qflx_qrgwl_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QRGWL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname( &
         'surface runoff at glaciers (liquid only), wetlands, lakes; also includes melted ice runoff from QSNWCPICE'), &
         ptr_col=this%qflx_qrgwl_col, c2l_scale_type='urbanf')

    this%qflx_drain_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRAI'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('sub-surface drainage'), &
         ptr_col=this%qflx_drain_col, c2l_scale_type='urbanf')

    this%qflx_liq_dynbal_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_LIQ_DYNBAL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('liq dynamic land cover change conversion runoff flux'), &
         ptr_lnd=this%qflx_liq_dynbal_grc)     

    this%qflx_ice_dynbal_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_ICE_DYNBAL'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('ice dynamic land cover change conversion runoff flux'), &
         ptr_lnd=this%qflx_ice_dynbal_grc)

    this%qflx_runoff_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QRUNOFF'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('total liquid runoff not including correction for land use change'), &
         ptr_col=this%qflx_runoff_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QRUNOFF_ICE'), &
         units='mm/s', avgflag='A', &
         long_name=this%info%lname('total liquid runoff not incl corret for LULCC (ice landunits only)'), &
         ptr_col=this%qflx_runoff_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_runoff_u_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QRUNOFF_U'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('Urban total runoff'), &
         ptr_col=this%qflx_runoff_u_col, set_nourb=spval, c2l_scale_type='urbanf', default='inactive')

    this%qflx_runoff_r_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QRUNOFF_R'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('Rural total runoff'), &
         ptr_col=this%qflx_runoff_r_col, set_spec=spval, default='inactive')

    this%qflx_snomelt_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNOMELT'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('snow melt rate'), &
         ptr_col=this%qflx_snomelt_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QSNOMELT_ICE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('snow melt (ice landunits only)'), &
         ptr_col=this%qflx_snomelt_col, c2l_scale_type='urbanf', l2g_scale_type='ice')


    this%qflx_snofrz_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNOFRZ'), &
         units='kg/m2/s', &
         avgflag='A', &
         long_name=this%info%lname('column-integrated snow freezing rate'), &
         ptr_col=this%qflx_snofrz_col, set_lake=spval, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QSNOFRZ_ICE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('column-integrated snow freezing rate (ice landunits only)'), &
         ptr_col=this%qflx_snofrz_col, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_snofrz_lyr_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%qflx_snofrz_lyr_col(begc:endc,-nlevsno+1:0)
    call hist_addfld2d ( &
         fname=this%info%fname('SNO_FRZ'),  &
         units='kg/m2/s', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('snow freezing rate in each snow layer'), &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d ( &
         fname=this%info%fname('SNO_FRZ_ICE'),  &
         units='mm/s', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('snow freezing rate in each snow layer (ice landunits only)'), &
         ptr_col=data2dptr, c2l_scale_type='urbanf',no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    this%qflx_prec_intr_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QINTR'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('interception'), &
         ptr_patch=this%qflx_prec_intr_patch, set_lake=0._r8)

    this%qflx_prec_grnd_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRIP'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('throughfall'), &
         ptr_patch=this%qflx_prec_grnd_patch, c2l_scale_type='urbanf')

    this%qflx_evap_soi_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSOIL'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname( 'Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)'), &
         ptr_patch=this%qflx_evap_soi_patch, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QSOIL_ICE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('Ground evaporation (ice landunits only)'), &
         ptr_patch=this%qflx_evap_soi_patch, c2l_scale_type='urbanf', l2g_scale_type='ice')

    this%qflx_evap_can_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QVEGE'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('canopy evaporation'), &
         ptr_patch=this%qflx_evap_can_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_tran_veg_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QVEGT'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('canopy transpiration'), &
         ptr_patch=this%qflx_tran_veg_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QSNOEVAP'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('evaporation from snow'), &
         ptr_patch=this%qflx_tran_veg_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_snwcp_liq_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNOCPLIQ'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('excess liquid h2o due to snow capping not including correction for land use change'), &
         ptr_col=this%qflx_snwcp_liq_col, c2l_scale_type='urbanf')

    this%qflx_snwcp_ice_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QSNWCPICE'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('excess solid h2o due to snow capping not including correction for land use change'), &
         ptr_col=this%qflx_snwcp_ice_col, c2l_scale_type='urbanf')

    this%qflx_glcice_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QICE'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('ice growth/melt'), &
         ptr_col=this%qflx_glcice_col, l2g_scale_type='ice')

    this%qflx_glcice_frz_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QICE_FRZ'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('ice growth'), &
         ptr_col=this%qflx_glcice_frz_col, l2g_scale_type='ice')

    this%qflx_glcice_melt_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QICE_MELT'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('ice melt'), &
         ptr_col=this%qflx_glcice_melt_col, l2g_scale_type='ice')

    this%qflx_rain_grnd_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_RAIN_GRND'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('rain on ground after interception'), &
         ptr_patch=this%qflx_rain_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_snow_grnd_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SNOW_GRND'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('snow on ground after interception'), &
         ptr_patch=this%qflx_snow_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_grnd_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_EVAP_GRND'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('ground surface evaporation'), &
         ptr_patch=this%qflx_evap_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_veg_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_EVAP_VEG'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('vegetation evaporation'), &
         ptr_patch=this%qflx_evap_veg_patch, default='inactive', c2l_scale_type='urbanf')

    this%qflx_evap_tot_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_EVAP_TOT'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('qflx_evap_soi + qflx_evap_can + qflx_tran_veg'), &
         ptr_patch=this%qflx_evap_tot_patch, c2l_scale_type='urbanf')

    this%qflx_dew_grnd_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_DEW_GRND'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('ground surface dew formation'), &
         ptr_patch=this%qflx_dew_grnd_patch, c2l_scale_type='urbanf')

    this%qflx_sub_snow_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SUB_SNOW'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('sublimation rate from snow pack (also includes bare ice sublimation from glacier columns)'), &
         ptr_patch=this%qflx_sub_snow_patch, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_SUB_SNOW_ICE'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('sublimation rate from snow pack (also includes bare ice sublimation from glacier columns) '// &
         '(ice landunits only)'), &
         ptr_patch=this%qflx_sub_snow_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%qflx_dew_snow_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QFLX_DEW_SNOW'), &
         units='mm H2O/s', &
         avgflag='A', &
         long_name=this%info%lname('surface dew added to snow pacK'), &
         ptr_patch=this%qflx_dew_snow_patch, c2l_scale_type='urbanf')

    this%qflx_rsub_sat_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QDRAI_XS'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('saturation excess drainage'), &
         ptr_col=this%qflx_rsub_sat_col, c2l_scale_type='urbanf')

  end subroutine InitHistory
  
  
  

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l
    !-----------------------------------------------------------------------

    this%qflx_evap_grnd_patch(bounds%begp:bounds%endp) = 0.0_r8
    this%qflx_dew_grnd_patch (bounds%begp:bounds%endp) = 0.0_r8
    this%qflx_dew_snow_patch (bounds%begp:bounds%endp) = 0.0_r8

    this%qflx_evap_grnd_col(bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_grnd_col (bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_snow_col (bounds%begc:bounds%endc) = 0.0_r8

    ! needed for CNNLeaching 
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%qflx_drain_col(c) = 0._r8
          this%qflx_surf_col(c)  = 0._r8
       end if
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterflux_type)            :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    ! needed for SNICAR
    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('qflx_snofrz_lyr'), &
         xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name=this%info%lname('snow layer ice freezing rate'), &
         units='kg m-2 s-1', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snofrz_lyr_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snofrz_lyr to zero
       this%qflx_snofrz_lyr_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

  end subroutine Restart

end module WaterFluxType
