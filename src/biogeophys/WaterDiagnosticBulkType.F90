module WaterDiagnosticBulkType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water diagnostic variables that just apply to bulk
  ! water. Diagnostic variables are neither fundamental state variables nor fluxes
  ! between those fundamental states, but are typically derived from those states and/or
  ! fluxes. Note that this type extends the base waterdiagnostic_type, so the full
  ! waterdiagnosticbulk_type contains the union of the fields defined here and the fields
  ! defined in waterdiagnostic_type.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use clm_varctl     , only : use_cn, iulog, use_luna
  use clm_varpar     , only : nlevgrnd, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use filterColMod   , only : filter_col_type, col_filter_from_ltypes
  use WaterDiagnosticType, only : waterdiagnostic_type
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  use WaterStateType, only : waterstate_type
  use WaterStateBulkType, only : waterstatebulk_type
  use WaterFluxType, only : waterflux_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, extends(waterdiagnostic_type), public :: waterdiagnosticbulk_type

     real(r8), pointer :: h2osno_total_col       (:)   ! col total snow water (mm H2O)
     real(r8), pointer :: snow_depth_col         (:)   ! col snow height of snow covered area (m)
     real(r8), pointer :: snowdp_col             (:)   ! col area-averaged snow height (m)
     real(r8), pointer :: snow_layer_unity_col   (:,:) ! value 1 for each snow layer, used for history diagnostics
     real(r8), pointer :: bw_col                 (:,:) ! col partial density of water in the snow pack (ice + liquid) [kg/m3] 

     real(r8), pointer :: h2osoi_liq_tot_col     (:)   ! vertically summed col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_tot_col     (:)   ! vertically summed col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: air_vol_col            (:,:) ! col air filled porosity
     real(r8), pointer :: h2osoi_liqvol_col      (:,:) ! col volumetric liquid water content (v/v)
     real(r8), pointer :: swe_old_col            (:,:) ! col initial snow water

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
     real(r8), pointer :: frac_sno_eff_col       (:)   ! col fraction of ground covered by snow (0 to 1) (note: this can be 1 even if there is no snow, but should be ignored in the no-snow case)
     real(r8), pointer :: frac_iceold_col        (:,:) ! col fraction of ice relative to the tot water (new) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: frac_h2osfc_col        (:)   ! col fractional area with surface water greater than zero
     real(r8), pointer :: frac_h2osfc_nosnow_col (:)   ! col fractional area with surface water greater than zero (if no snow present)
     real(r8), pointer :: wf_col                 (:)   ! col soil water as frac. of whc for top 0.05 m (0-1) 
     real(r8), pointer :: wf2_col                (:)   ! col soil water as frac. of whc for top 0.17 m (0-1) 
     real(r8), pointer :: fwet_patch             (:)   ! patch canopy fraction that is wet (0 to 1)
     real(r8), pointer :: fcansno_patch          (:)   ! patch canopy fraction that is snow covered (0 to 1)
     real(r8), pointer :: fdry_patch             (:)   ! patch canopy fraction of foliage that is green and dry [-] (new)

     ! Summed fluxes
     real(r8), pointer :: qflx_prec_intr_patch   (:)   ! patch interception of precipitation (mm H2O/s)
     real(r8), pointer :: qflx_prec_grnd_col     (:)   ! col water onto ground including canopy runoff (mm H2O/s)

   contains

     procedure, public  :: InitBulk
     procedure, public  :: RestartBulk
     procedure, public  :: Summary
     procedure, public  :: ResetBulkFilter
     procedure, public  :: ResetBulk
     procedure, private :: InitBulkAllocate 
     procedure, private :: InitBulkHistory  
     procedure, private :: InitBulkCold     
     procedure, private :: RestartBackcompatIssue783

  end type waterdiagnosticbulk_type

   ! PUBLIC MEMBER FUNCTIONS
  public :: readParams

  type, private :: params_type
      real(r8) :: zlnd  ! Momentum roughness length for soil, glacier, wetland (m)
  end type params_type
  type(params_type), private ::  params_inst

  ! minimum allowed snow effective radius (also "fresh snow" value) [microns]
  real(r8), public, parameter :: snw_rds_min = 54.526_r8    

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

 subroutine readParams( ncid )
    !
    ! !USES:
    use ncdio_pio, only: file_desc_t
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'readParams_WaterDiagnosticBulk'
    !--------------------------------------------------------------------

    ! Momentum roughness length for soil, glacier, wetland (m)
    call readNcdioScalar(ncid, 'zlnd', subname, params_inst%zlnd)

  end subroutine readParams

  !------------------------------------------------------------------------
  subroutine InitBulk(this, bounds, info, vars, &
       snow_depth_input_col, h2osno_input_col)

    class(waterdiagnosticbulk_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds  
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: vars
    real(r8)          , intent(in) :: snow_depth_input_col(bounds%begc:)
    real(r8)          , intent(in) :: h2osno_input_col(bounds%begc:)  ! Initial total snow water (mm H2O)


    call this%Init(bounds, info, vars)

    call this%InitBulkAllocate(bounds) 

    call this%InitBulkHistory(bounds)

    call this%InitBulkCold(bounds, &
       snow_depth_input_col, h2osno_input_col)

  end subroutine InitBulk

  !------------------------------------------------------------------------
  subroutine InitBulkAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate(this%h2osno_total_col       (begc:endc))                     ; this%h2osno_total_col       (:)   = nan
    allocate(this%snow_depth_col         (begc:endc))                     ; this%snow_depth_col         (:)   = nan
    allocate(this%snowdp_col             (begc:endc))                     ; this%snowdp_col             (:)   = nan
    allocate(this%snow_layer_unity_col   (begc:endc,-nlevsno+1:0))        ; this%snow_layer_unity_col   (:,:) = nan
    allocate(this%bw_col                 (begc:endc,-nlevsno+1:0))        ; this%bw_col                 (:,:) = nan   
    allocate(this%air_vol_col            (begc:endc, 1:nlevgrnd))         ; this%air_vol_col            (:,:) = nan
    allocate(this%h2osoi_liqvol_col      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liqvol_col      (:,:) = nan
    allocate(this%h2osoi_ice_tot_col     (begc:endc))                     ; this%h2osoi_ice_tot_col     (:)   = nan
    allocate(this%h2osoi_liq_tot_col     (begc:endc))                     ; this%h2osoi_liq_tot_col     (:)   = nan
    allocate(this%swe_old_col            (begc:endc,-nlevsno+1:0))        ; this%swe_old_col            (:,:) = nan   

    allocate(this%snw_rds_col            (begc:endc,-nlevsno+1:0))        ; this%snw_rds_col            (:,:) = nan
    allocate(this%snw_rds_top_col        (begc:endc))                     ; this%snw_rds_top_col        (:)   = nan
    allocate(this%h2osno_top_col         (begc:endc))                     ; this%h2osno_top_col         (:)   = nan
    allocate(this%sno_liq_top_col        (begc:endc))                     ; this%sno_liq_top_col        (:)   = nan

    allocate(this%dqgdT_col              (begc:endc))                     ; this%dqgdT_col              (:)   = nan   
    allocate(this%rh_ref2m_patch         (begp:endp))                     ; this%rh_ref2m_patch         (:)   = nan
    allocate(this%rh_ref2m_u_patch       (begp:endp))                     ; this%rh_ref2m_u_patch       (:)   = nan
    allocate(this%rh_ref2m_r_patch       (begp:endp))                     ; this%rh_ref2m_r_patch       (:)   = nan
    allocate(this%rh_af_patch            (begp:endp))                     ; this%rh_af_patch            (:)   = nan
    allocate(this%rh10_af_patch          (begp:endp))                     ; this%rh10_af_patch          (:)   = spval

    allocate(this%frac_sno_col           (begc:endc))                     ; this%frac_sno_col           (:)   = nan
    allocate(this%frac_sno_eff_col       (begc:endc))                     ; this%frac_sno_eff_col       (:)   = nan
    allocate(this%frac_iceold_col        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%frac_iceold_col        (:,:) = nan
    allocate(this%frac_h2osfc_col        (begc:endc))                     ; this%frac_h2osfc_col        (:)   = nan 
    allocate(this%frac_h2osfc_nosnow_col (begc:endc))                     ; this%frac_h2osfc_nosnow_col        (:)   = nan 
    allocate(this%wf_col                 (begc:endc))                     ; this%wf_col                 (:)   = nan
    allocate(this%wf2_col                (begc:endc))                     ; this%wf2_col                (:)   = nan
    allocate(this%fwet_patch             (begp:endp))                     ; this%fwet_patch             (:)   = nan
    allocate(this%fcansno_patch          (begp:endp))                     ; this%fcansno_patch          (:)   = nan
    allocate(this%fdry_patch             (begp:endp))                     ; this%fdry_patch             (:)   = nan
    allocate(this%qflx_prec_intr_patch   (begp:endp))                     ; this%qflx_prec_intr_patch   (:)   = nan
    allocate(this%qflx_prec_grnd_col     (begc:endc))                     ; this%qflx_prec_grnd_col     (:)   = nan

  end subroutine InitBulkAllocate

  !------------------------------------------------------------------------
  subroutine InitBulkHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal, no_snow_zero
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    this%h2osno_total_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('H2OSNO'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('snow depth (liquid water)'), &
         ptr_col=this%h2osno_total_col, c2l_scale_type='urbanf')
    call hist_addfld1d ( &
         fname=this%info%fname('H2OSNO_ICE'), &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('snow depth (liquid water, ice landunits only)'), &
         ptr_col=this%h2osno_total_col, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%h2osoi_liq_tot_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('TOTSOILLIQ'),  &
         units='kg/m2', &
         avgflag='A', &
         long_name=this%info%lname('vertically summed soil liquid water (veg landunits only)'), &
         ptr_col=this%h2osoi_liq_tot_col, l2g_scale_type='veg')

    this%h2osoi_ice_tot_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('TOTSOILICE'),  &
         units='kg/m2', &
         avgflag='A', &
         long_name=this%info%lname('vertically summed soil cie (veg landunits only)'), &
         ptr_col=this%h2osoi_ice_tot_col, l2g_scale_type='veg')

    this%rh_ref2m_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('RH2M'), &
         units='%',  &
         avgflag='A', &
         long_name=this%info%lname('2m relative humidity'), &
         ptr_patch=this%rh_ref2m_patch)

    this%rh_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('RH2M_R'), &
         units='%',  &
         avgflag='A', &
         long_name=this%info%lname('Rural 2m specific humidity'), &
         ptr_patch=this%rh_ref2m_r_patch, set_spec=spval, default='inactive')

    this%rh_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('RH2M_U'), &
         units='%',  &
         avgflag='A', &
         long_name=this%info%lname('Urban 2m relative humidity'), &
         ptr_patch=this%rh_ref2m_u_patch, set_nourb=spval, default='inactive')

    this%rh_af_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('RHAF'), &
         units='fraction', &
         avgflag='A', &
         long_name=this%info%lname('fractional humidity of canopy air'), &
         ptr_patch=this%rh_af_patch, set_spec=spval, default='inactive')

    if(use_luna)then
       call hist_addfld1d ( &
            fname=this%info%fname('RHAF10'), &
            units='fraction', &
            avgflag='A', &
            long_name=this%info%lname('10 day running mean of fractional humidity of canopy air'), &
            ptr_patch=this%rh10_af_patch, set_spec=spval, default='inactive')
    endif

    ! Fractions

    this%frac_h2osfc_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('FH2OSFC'),  &
         units='unitless',  &
         avgflag='A', &
         long_name=this%info%lname('fraction of ground covered by surface water'), &
         ptr_col=this%frac_h2osfc_col)

    this%frac_h2osfc_nosnow_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('FH2OSFC_NOSNOW'),  &
         units='unitless',  &
         avgflag='A', &
         long_name=this%info%lname('fraction of ground covered by surface water (if no snow present)'), &
         ptr_col=this%frac_h2osfc_nosnow_col, default='inactive')

    this%frac_sno_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('FSNO'),  &
         units='unitless',  &
         avgflag='A', &
         long_name=this%info%lname('fraction of ground covered by snow'), &
         ptr_col=this%frac_sno_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('FSNO_ICE'),  &
         units='unitless',  &
         avgflag='A', &
         long_name=this%info%lname('fraction of ground covered by snow (ice landunits only)'), &
         ptr_col=this%frac_sno_col, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%frac_sno_eff_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('FSNO_EFF'),  &
         units='unitless',  &
         avgflag='A', &
         long_name=this%info%lname('effective fraction of ground covered by snow'), &
         ptr_col=this%frac_sno_eff_col, c2l_scale_type='urbanf')!, default='inactive')

    if (use_cn) then
       this%fwet_patch(begp:endp) = spval
       call hist_addfld1d ( &
            fname=this%info%fname('FWET'), &
            units='proportion', &
            avgflag='A', &
            long_name=this%info%lname('fraction of canopy that is wet'), &
            ptr_patch=this%fwet_patch, default='inactive')
    end if

    if (use_cn) then
       this%fcansno_patch(begp:endp) = spval
       call hist_addfld1d ( &
            fname=this%info%fname('FCANSNO'), &
            units='proportion', &
            avgflag='A', &
            long_name=this%info%lname('fraction of canopy that is wet'), &
            ptr_patch=this%fcansno_patch, default='inactive')
    end if

    if (use_cn) then
       this%fdry_patch(begp:endp) = spval
       call hist_addfld1d ( &
            fname=this%info%fname('FDRY'), &
            units='proportion', &
            avgflag='A', &
            long_name=this%info%lname('fraction of foliage that is green and dry'), &
            ptr_patch=this%fdry_patch, default='inactive')
    end if

    if (use_cn)then
       this%frac_iceold_col(begc:endc,:) = spval
       call hist_addfld2d ( &
            fname=this%info%fname('FRAC_ICEOLD'), &
            units='proportion', type2d='levgrnd', &
            avgflag='A', &
            long_name=this%info%lname('fraction of ice relative to the tot water'), &
            ptr_col=this%frac_iceold_col, default='inactive')
    end if

    ! Snow properties - these will be vertically averaged over the snow profile

    this%snow_depth_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOW_DEPTH'),  &
         units='m',  &
         avgflag='A', &
         long_name=this%info%lname('snow height of snow covered area'), &
         ptr_col=this%snow_depth_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('SNOW_DEPTH_ICE'), &
         units='m',  &
         avgflag='A', &
         long_name=this%info%lname('snow height of snow covered area (ice landunits only)'), &
         ptr_col=this%snow_depth_col, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%snowdp_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOWDP'),  &
         units='m',  &
         avgflag='A', &
         long_name=this%info%lname('gridcell mean snow height'), &
         ptr_col=this%snowdp_col, c2l_scale_type='urbanf')

    if (use_cn) then
       this%wf_col(begc:endc) = spval
       call hist_addfld1d ( &
            fname=this%info%fname('WF'), &
            units='proportion', &
            avgflag='A', &
            long_name=this%info%lname('soil water as frac. of whc for top 0.05 m'), &
            ptr_col=this%wf_col, default='inactive')
    end if

    this%h2osno_top_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('H2OSNO_TOP'), &
         units='kg/m2', &
         avgflag='A', &
         long_name=this%info%lname('mass of snow in top snow layer'), &
         ptr_col=this%h2osno_top_col, set_urb=spval)

    this%snw_rds_top_col(begc:endc) = spval 
    call hist_addfld1d ( &
         fname=this%info%fname('SNORDSL'), &
         units='m^-6', &
         avgflag='A', &
         long_name=this%info%lname('top snow layer effective grain radius'), &
         ptr_col=this%snw_rds_top_col, set_urb=spval, default='inactive')

    this%sno_liq_top_col(begc:endc) = spval 
    call hist_addfld1d ( &
         fname=this%info%fname('SNOLIQFL'), &
         units='fraction', &
         avgflag='A', &
         long_name=this%info%lname('top snow layer liquid water fraction (land)'), &
         ptr_col=this%sno_liq_top_col, set_urb=spval, default='inactive')

    ! We determine the fractional time (and fraction of the grid cell) over which each
    ! snow layer existed by running the snow averaging routine on a field whose value is 1
    ! everywhere
    data2dptr => this%snow_layer_unity_col(:,-nlevsno+1:0)
    call hist_addfld2d ( &
         fname=this%info%fname('SNO_EXISTENCE'), &
         units='unitless', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('Fraction of averaging period for which each snow layer existed'), &
         ptr_col=data2dptr, no_snow_behavior=no_snow_zero, default='inactive')

    this%bw_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%bw_col(:,-nlevsno+1:0)
    call hist_addfld2d ( &
         fname=this%info%fname('SNO_BW'), &
         units='kg/m3', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('Partial density of water in the snow pack (ice + liquid)'), &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d ( &
         fname=this%info%fname('SNO_BW_ICE'), &
         units='kg/m3', type2d='levsno', &
         avgflag='A', &
         long_name=this%info%lname('Partial density of water in the snow pack (ice + liquid, ice landunits only)'), &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    this%snw_rds_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%snw_rds_col(:,-nlevsno+1:0)
    call hist_addfld2d ( &
         fname=this%info%fname('SNO_GS'), &
         units='Microns', type2d='levsno',  &
         avgflag='A', &
         long_name=this%info%lname('Mean snow grain size'), &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d ( &
         fname=this%info%fname('SNO_GS_ICE'), &
         units='Microns', type2d='levsno',  &
         avgflag='A', &
         long_name=this%info%lname('Mean snow grain size (ice landunits only)'), &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    ! Summed fluxes

    this%qflx_prec_intr_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('QINTR'), &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('interception'), &
         ptr_patch=this%qflx_prec_intr_patch, set_lake=0._r8)

  end subroutine InitBulkHistory

  !-----------------------------------------------------------------------
  subroutine InitBulkCold(this, bounds, &
       snow_depth_input_col, h2osno_input_col)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type), intent(in) :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: snow_depth_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: h2osno_input_col(bounds%begc:)  ! Initial total snow water (mm H2O)
    !
    ! !LOCAL VARIABLES:
    integer            :: c,l
    real(r8)           :: snowbd      ! temporary calculation of snow bulk density (kg/m3)
    real(r8)           :: fmelt       ! snowbd/100
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(snow_depth_input_col) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(h2osno_input_col) == (/bounds%endc/)), sourcefile, __LINE__)

    do c = bounds%begc,bounds%endc
       this%snow_depth_col(c)         = snow_depth_input_col(c)
       this%snow_layer_unity_col(c,:) = 1._r8
    end do

    do c = bounds%begc,bounds%endc
       this%wf_col(c) = spval
       this%wf2_col(c) = spval
    end do


    associate(snl => col%snl) 

      this%frac_h2osfc_col(bounds%begc:bounds%endc) = 0._r8

      this%fwet_patch(bounds%begp:bounds%endp) = 0._r8
      this%fdry_patch(bounds%begp:bounds%endp) = 0._r8
      this%fcansno_patch(bounds%begp:bounds%endp) = 0._r8

      this%qflx_prec_intr_patch(bounds%begp:bounds%endp) = 0._r8

      !--------------------------------------------
      ! Set snow water
      !--------------------------------------------

      ! Note: Glacier_mec columns are initialized with half the maximum snow cover.
      ! This gives more realistic values of qflx_glcice sooner in the simulation
      ! for columns with net ablation, at the cost of delaying ice formation
      ! in columns with net accumulation.

      do c = bounds%begc, bounds%endc
         l = col%landunit(c)
         if (lun%urbpoi(l)) then
            ! From Bonan 1996 (LSM technical note)
            this%frac_sno_col(c) = min( this%snow_depth_col(c)/0.05_r8, 1._r8)
         else
            this%frac_sno_col(c) = 0._r8
            ! snow cover fraction as in Niu and Yang 2007
            if(this%snow_depth_col(c) > 0.0)  then
               snowbd   = min(400._r8, h2osno_input_col(c)/this%snow_depth_col(c)) !bulk density of snow (kg/m3)
               fmelt    = (snowbd/100.)**1.
               ! 100 is the assumed fresh snow density; 1 is a melting factor that could be
               ! reconsidered, optimal value of 1.5 in Niu et al., 2007
               this%frac_sno_col(c) = tanh( this%snow_depth_col(c) / (2.5 * params_inst%zlnd * fmelt) )
            endif
         end if
      end do

      do c = bounds%begc,bounds%endc
         if (snl(c) < 0) then
            this%snw_rds_col(c,snl(c)+1:0)        = snw_rds_min
            this%snw_rds_col(c,-nlevsno+1:snl(c)) = 0._r8
            this%snw_rds_top_col(c)               = snw_rds_min
         elseif (h2osno_input_col(c) > 0._r8) then
            this%snw_rds_col(c,0)                 = snw_rds_min
            this%snw_rds_col(c,-nlevsno+1:-1)     = 0._r8
            this%snw_rds_top_col(c)               = spval
            this%sno_liq_top_col(c)               = spval
         else
            this%snw_rds_col(c,:)                 = 0._r8
            this%snw_rds_top_col(c)               = spval
            this%sno_liq_top_col(c)               = spval
         endif
      end do


    end associate

  end subroutine InitBulkCold

  !------------------------------------------------------------------------
  subroutine RestartBulk(this, bounds, ncid, flag, writing_finidat_interp_dest_file, waterstatebulk_inst)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod          , only : masterproc
    use clm_varcon       , only : pondmx, watmin, spval, nameg
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varctl       , only : bound_h2osoi
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    logical, intent(in) :: writing_finidat_interp_dest_file ! true if we are writing a finidat_interp_dest file (ignored for flag=='read')
    type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    logical  :: readvar
    !------------------------------------------------------------------------


    call this%Restart(bounds, ncid, flag=flag)

    if(use_luna)then
       call restartvar(ncid=ncid, flag=flag, &
            varname=this%info%fname('rh10'), &
            xtype=ncd_double,  &
            dim1name='pft', &
            long_name=this%info%lname('10-day mean boundary layer relative humidity'), &
            units='unitless', &
            interpinic_flag='interp', readvar=readvar, data=this%rh10_af_patch)
    endif

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('FH2OSFC'), &
         xtype=ncd_double,  &
         dim1name='column',&
         long_name=this%info%lname('fraction of ground covered by h2osfc (0 to 1)'), &
         units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frac_h2osfc_col)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('SNOW_DEPTH'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('snow depth'), &
         units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%snow_depth_col) 

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('frac_sno_eff'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('fraction of ground covered by snow (0 to 1)'),&
         units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%frac_sno_eff_col)
    if (flag == 'read' .and. .not. readvar) then
       this%frac_sno_eff_col(bounds%begc:bounds%endc) = 0.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('frac_sno'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('fraction of ground covered by snow (0 to 1)'),&
         units='unitless',&
         interpinic_flag='interp', readvar=readvar, data=this%frac_sno_col)
    call this%RestartBackcompatIssue783( &
         bounds = bounds, &
         ncid = ncid, &
         flag = flag, &
         writing_finidat_interp_dest_file = writing_finidat_interp_dest_file, &
         waterstatebulk_inst = waterstatebulk_inst)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('FWET'), &
         xtype=ncd_double,  &
         dim1name='pft', &
         long_name=this%info%lname('fraction of canopy that is wet (0 to 1)'), &
         units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fwet_patch)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('FCANSNO'), &
         xtype=ncd_double,  &
         dim1name='pft', &
         long_name=this%info%lname('fraction of canopy that is snow covered (0 to 1)'), &
         units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fcansno_patch)

    ! column type physical state variable - snw_rds
    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('snw_rds'), &
         xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name=this%info%lname('snow layer effective radius'), &
         units='um', &
         interpinic_flag='interp', readvar=readvar, data=this%snw_rds_col)
    if (flag == 'read' .and. .not. readvar) then
       ! NOTE(wjs, 2018-08-03) There was some code here that looked like it was just for
       ! the sake of backwards compatibility, dating back to 2014 or earlier. I was
       ! tempted to just remove it, but on the off-chance that this conditional is still
       ! ever entered, I'm putting an endrun call here to notify users of this removed
       ! code.
       if (masterproc) then
          write(iulog,*) "SNICAR: This is an initial run (not a restart), and grain size/aerosol " // &
               "mass data are not defined in initial condition file. This situation is no longer handled."
       endif
       call endrun(msg = "Absent snw_rds on initial conditions file no longer handled. "// &
            errMsg(sourcefile, __LINE__))
    endif

    if (use_cn) then
       call restartvar(ncid=ncid, flag=flag, &
            varname=this%info%fname('wf'), &
            xtype=ncd_double,  &
            dim1name='column', &
            long_name=this%info%lname(''), &
            units='', &
            interpinic_flag='interp', readvar=readvar, data=this%wf_col) 
    end if



  end subroutine RestartBulk

  !-----------------------------------------------------------------------
  subroutine RestartBackcompatIssue783(this, bounds, ncid, flag, &
       writing_finidat_interp_dest_file, waterstatebulk_inst)
    !
    ! !DESCRIPTION:
    ! Apply backwards compatibility corrections to address issue ESCOMP/ctsm#783
    !
    ! BACKWARDS_COMPATIBILITY(wjs, 2019-10-15) Due to ESCOMP/ctsm#783, old restart files
    ! can have frac_sno == 0 for lake points despite having a snow pack. This can cause
    ! other problems, so fix that here. However, it is apparently possible for frac_sno to
    ! be 0 legitimately when h2osno_total > 0. So if we apply this correction always, then
    ! we sometimes introduce unintentional changes to newer restart files where we don't
    ! actually need to apply this correction. We avoid this by writing metadata to the
    ! restart file indicating that it's new enough to have this correction already in
    ! place, then avoiding doing the correction here if we find we're working with a
    ! new-enough restart file.
    !
    ! This backwards compatibility code can be removed once we can rely on all restart
    ! files being new enough. i.e., we can remove this code once we can rely all restart
    ! files having this new piece of metadata (at which point we can also stop writing
    ! this metadata, as long as we don't need to use newer restart files with older code
    ! versions).
    !
    ! !USES:
    use ncdio_pio        , only : file_desc_t
    use IssueFixedMetadataHandler, only : write_issue_fixed_metadata, read_issue_fixed_metadata
    use landunit_varcon  , only : istdlak
    use clm_time_manager , only : is_restart
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    logical, intent(in) :: writing_finidat_interp_dest_file ! true if this is a finidat_interp_dest file
    type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    integer  :: attribute_value
    logical  :: do_correction
    real(r8) :: h2osno_total(bounds%begc:bounds%endc)  ! total snow water (mm H2O)
    type(filter_col_type) :: filter_lakec  ! filter for lake columns

    integer, parameter :: issue_num = 783

    character(len=*), parameter :: subname = 'RestartBackcompatIssue783'
    !-----------------------------------------------------------------------

    if (flag == 'define') then
       call write_issue_fixed_metadata( &
            ncid = ncid, &
            writing_finidat_interp_dest_file = writing_finidat_interp_dest_file, &
            issue_num = issue_num)

    else if (flag == 'read' .and. .not. is_restart()) then
       call read_issue_fixed_metadata( &
            ncid = ncid, &
            issue_num = issue_num, &
            attribute_value = attribute_value)
       if (attribute_value == 0) then
          do_correction = .true.
       else
          do_correction = .false.
       end if

       if (do_correction) then
          filter_lakec = col_filter_from_ltypes( &
               bounds = bounds, &
               ltypes = [istdlak], &
               include_inactive = .true.)
          call waterstatebulk_inst%CalculateTotalH2osno( &
               bounds = bounds, &
               num_c = filter_lakec%num, &
               filter_c = filter_lakec%indices, &
               caller = 'WaterDiagnosticBulkType_RestartBulk', &
               h2osno_total = h2osno_total(bounds%begc:bounds%endc))
          do fc = 1, filter_lakec%num
             c = filter_lakec%indices(fc)
             if (this%frac_sno_col(c) == 0._r8 .and. h2osno_total(c) > 0._r8) then
                ! Often the value should be between 0 and 1 rather than being 1, but 1 is at
                ! least better than 0 in this case, and it would be tricky or impossible to
                ! figure out the "correct" value.
                this%frac_sno_col(c) = 1._r8
             end if
          end do
       end if
    end if

  end subroutine RestartBackcompatIssue783

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, &
       num_soilp, filter_soilp, &
       num_allc, filter_allc, &
       num_nolakec, filter_nolakec, &
       waterstate_inst, waterflux_inst)
    !
    ! !DESCRIPTION:
    ! Compute end-of-timestep summaries of water diagnostic terms
    !
    ! !USES:
    use clm_varpar   , only : nlevsoi
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilp          ! number of patches in soilp filter
    integer                     , intent(in)    :: filter_soilp(:)    ! filter for soil patches
    integer                     , intent(in)    :: num_allc           ! number of columns in allc filter
    integer                     , intent(in)    :: filter_allc(:)     ! filter for all columns
    integer                     , intent(in)    :: num_nolakec        ! number of columns in no-lake columnc filter
    integer                     , intent(in)    :: filter_nolakec(:)  ! filter for no-lake  columns
    class(waterstate_type)      , intent(in)    :: waterstate_inst
    class(waterflux_type)       , intent(in)    :: waterflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p, j, l, fc, c            ! Indices
    real(r8):: fracl                         ! fraction of soil layer contributing to 10cm total soil water

    character(len=*), parameter :: subname = 'Summary'
    !-----------------------------------------------------------------------
    associate(                                                 &
         dz                 => col%dz                        , & !  Input:  [real(r8) (:,:) ]  layer thickness depth (m)             
         zi                 => col%zi                        , & !  Input:  [real(r8) (:,:) ]  interface depth (m)  

         h2osoi_ice         => waterstate_inst%h2osoi_ice_col, & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq         => waterstate_inst%h2osoi_liq_col, & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)

         h2osoi_ice_tot     => this%h2osoi_ice_tot_col       , & ! Output: [real(r8) (:)   ]  vertically summed ice lens (kg/m2)
         h2osoi_liq_tot     => this%h2osoi_liq_tot_col       , & ! Output: [real(r8) (:)   ]  vertically summed liquid water (kg/m2)   
         h2osoi_liqice_10cm => this%h2osoi_liqice_10cm_col     & ! Output: [real(r8) (:)   ]  liquid water + ice lens in top 10cm of soil (kg/m2)
    )

    call this%waterdiagnostic_type%Summary(bounds, &
         num_soilp, filter_soilp, &
         num_allc, filter_allc, &
         num_nolakec, filter_nolakec, &
         waterstate_inst, waterflux_inst)

    call waterstate_inst%CalculateTotalH2osno(bounds, num_allc, filter_allc, &
         caller = 'WaterDiagnosticBulkType:Summary', &
         h2osno_total = this%h2osno_total_col(bounds%begc:bounds%endc))

    do fp = 1, num_soilp
       p = filter_soilp(fp)
       this%qflx_prec_intr_patch(p) = &
            waterflux_inst%qflx_intercepted_liq_patch(p) + &
            waterflux_inst%qflx_intercepted_snow_patch(p)
    end do

    do fc = 1, num_allc
       c = filter_allc(fc)
       this%qflx_prec_grnd_col(c) = &
            waterflux_inst%qflx_liq_grnd_col(c) + &
            waterflux_inst%qflx_snow_grnd_col(c)
    end do
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       l = col%landunit(c)
       if (.not. lun%urbpoi(l)) then
          h2osoi_liqice_10cm(c) = 0.0_r8
          h2osoi_liq_tot(c) = 0._r8
          h2osoi_ice_tot(c) = 0._r8
       end if
    end do
    do j = 1, nlevsoi
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = col%landunit(c)
          if (.not. lun%urbpoi(l)) then
             if (zi(c,j) <= 0.1_r8) then
                fracl = 1._r8
                h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                     (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                     fracl
             else
                if (zi(c,j) > 0.1_r8 .and. zi(c,j-1) < 0.1_r8) then
                   fracl = (0.1_r8 - zi(c,j-1))/dz(c,j)
                   h2osoi_liqice_10cm(c) = h2osoi_liqice_10cm(c) + &
                        (h2osoi_liq(c,j)+h2osoi_ice(c,j))* &
                        fracl
                end if
             end if
             h2osoi_liq_tot(c) = h2osoi_liq_tot(c) + h2osoi_liq(c,j)
             h2osoi_ice_tot(c) = h2osoi_ice_tot(c) + h2osoi_ice(c,j)
          end if
       end do
    end do

    end associate

  end subroutine Summary

  !-----------------------------------------------------------------------
  subroutine ResetBulkFilter(this, num_c, filter_c)
    !
    ! !DESCRIPTION:
    ! Initialize SNICAR variables for fresh snow columns, for all columns in the given
    ! filter
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type), intent(inout) :: this
    integer, intent(in) :: num_c       ! number of columns in filter_c
    integer, intent(in) :: filter_c(:) ! column filter to operate over
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'ResetBulkFilter'
    !-----------------------------------------------------------------------

    do fc = 1, num_c
       c = filter_c(fc)
       call this%ResetBulk(c)
    end do

  end subroutine ResetBulkFilter

  !-----------------------------------------------------------------------
  subroutine ResetBulk(this, column)
    !
    ! !DESCRIPTION:
    ! Intitialize SNICAR variables for fresh snow column
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type), intent(inout) :: this
    integer , intent(in)   :: column     ! column index
    !-----------------------------------------------------------------------

    this%snw_rds_col(column,0)  = snw_rds_min

  end subroutine ResetBulk

end module WaterDiagnosticBulkType
