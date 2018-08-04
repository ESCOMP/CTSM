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
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use WaterStateBulkType, only : waterstatebulk_type
  use WaterDiagnosticType, only : waterdiagnostic_type
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, extends(waterdiagnostic_type), public :: waterdiagnosticbulk_type

     real(r8), pointer :: snow_depth_col         (:)   ! col snow height of snow covered area (m)
     real(r8), pointer :: snowdp_col             (:)   ! col area-averaged snow height (m)
     real(r8), pointer :: snow_layer_unity_col   (:,:) ! value 1 for each snow layer, used for history diagnostics
     real(r8), pointer :: bw_col                 (:,:) ! col partial density of water in the snow pack (ice + liquid) [kg/m3] 

     real(r8), pointer :: h2osoi_liq_tot_col     (:)   ! vertically summed col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_tot_col     (:)   ! vertically summed col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: air_vol_col            (:,:) ! col air filled porosity
     real(r8), pointer :: h2osoi_liqvol_col      (:,:) ! col volumetric liquid water content (v/v)
     real(r8), pointer :: snounload_patch        (:)   ! Canopy snow unloading (mm H2O)
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


   contains

     procedure          :: InitBulk         
     procedure          :: RestartBulk      
     procedure, public  :: ResetBulk 
     procedure, private :: InitBulkAllocate 
     procedure, private :: InitBulkHistory  
     procedure, private :: InitBulkCold     

  end type waterdiagnosticbulk_type

  ! minimum allowed snow effective radius (also "fresh snow" value) [microns]
  real(r8), public, parameter :: snw_rds_min = 54.526_r8    

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine InitBulk(this, bounds, info, &
       snow_depth_input_col, waterstatebulk_inst)

    class(waterdiagnosticbulk_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    class(water_info_base_type), intent(in), target :: info
    real(r8)          , intent(in) :: snow_depth_input_col(bounds%begc:)
    class(waterstatebulk_type), intent(in) :: waterstatebulk_inst


    call this%Init(bounds, info)

    call this%InitBulkAllocate(bounds) 

    call this%InitBulkHistory(bounds)

    call this%InitBulkCold(bounds, &
       snow_depth_input_col, waterstatebulk_inst)

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
    class(waterdiagnosticbulk_type) :: this
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

    allocate(this%snow_depth_col         (begc:endc))                     ; this%snow_depth_col         (:)   = nan
    allocate(this%snowdp_col             (begc:endc))                     ; this%snowdp_col             (:)   = nan
    allocate(this%snow_layer_unity_col   (begc:endc,-nlevsno+1:0))        ; this%snow_layer_unity_col   (:,:) = nan
    allocate(this%bw_col                 (begc:endc,-nlevsno+1:0))        ; this%bw_col                 (:,:) = nan   
    allocate(this%air_vol_col            (begc:endc, 1:nlevgrnd))         ; this%air_vol_col            (:,:) = nan
    allocate(this%h2osoi_liqvol_col      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liqvol_col      (:,:) = nan
    allocate(this%h2osoi_ice_tot_col     (begc:endc))                     ; this%h2osoi_ice_tot_col     (:)   = nan
    allocate(this%h2osoi_liq_tot_col     (begc:endc))                     ; this%h2osoi_liq_tot_col     (:)   = nan
    allocate(this%snounload_patch        (begp:endp))                     ; this%snounload_patch        (:)   = nan  
    allocate(this%swe_old_col            (begc:endc,-nlevsno+1:0))        ; this%swe_old_col            (:,:) = nan   

    allocate(this%total_plant_stored_h2o_col(begc:endc))                  ; this%total_plant_stored_h2o_col(:) = nan

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
    allocate(this%wf2_col                (begc:endc))                     ; 
    allocate(this%fwet_patch             (begp:endp))                     ; this%fwet_patch             (:)   = nan
    allocate(this%fcansno_patch          (begp:endp))                     ; this%fcansno_patch          (:)   = nan
    allocate(this%fdry_patch             (begp:endp))                     ; this%fdry_patch             (:)   = nan

  end subroutine InitBulkAllocate

  !------------------------------------------------------------------------
  subroutine InitBulkHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_cn, use_lch4
    use clm_varctl     , only : hist_wrtch4diag
    use clm_varpar     , only : nlevsno, nlevsoi
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal, no_snow_zero
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type) :: this
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

    this%h2osoi_liq_tot_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('TOTSOILLIQ'),  &
         units='kg/m2', &
         avgflag='A', &
         long_name=this%info%lname('vertically summed soil liquid water (veg landunits only)'), &
         ptr_col=this%h2osoi_liq_tot_col, set_urb=spval, set_lake=spval, l2g_scale_type='veg')

    this%h2osoi_ice_tot_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('TOTSOILICE'),  &
         units='kg/m2', &
         avgflag='A', &
         long_name=this%info%lname('vertically summed soil cie (veg landunits only)'), &
         ptr_col=this%h2osoi_ice_tot_col, set_urb=spval, set_lake=spval, l2g_scale_type='veg')

    this%snounload_patch(begp:endp) = spval 
    call hist_addfld1d ( &
         fname=this%info%fname('SNOUNLOAD'), &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('Canopy snow unloading'), &
         ptr_patch=this%snounload_patch, set_lake=0._r8)


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

  end subroutine InitBulkHistory

  !-----------------------------------------------------------------------
  subroutine InitBulkCold(this, bounds, &
       snow_depth_input_col, waterstatebulk_inst)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    use shr_const_mod   , only : shr_const_pi
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use shr_spfn_mod    , only : shr_spfn_erf
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use landunit_varcon , only : istwet, istsoil, istdlak, istcrop, istice_mec  
    use column_varcon   , only : icol_shadewall, icol_road_perv
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use clm_varcon      , only : spval, sb, bdsno 
    use clm_varcon      , only : zlnd, tfrz, spval, pc
    use clm_varctl      , only : fsurdat, iulog
    use clm_varctl        , only : use_bedrock
    use spmdMod         , only : masterproc
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type)                :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: snow_depth_input_col(bounds%begc:)
    class(waterstatebulk_type), intent(in)                :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev
    real(r8)           :: maxslope, slopemax, minslope
    real(r8)           :: d, fd, dfdd, slope0,slopebeta
    real(r8) ,pointer  :: std (:)     
    logical            :: readvar 
    type(file_desc_t)  :: ncid        
    character(len=256) :: locfn       
    real(r8)           :: snowbd      ! temporary calculation of snow bulk density (kg/m3)
    real(r8)           :: fmelt       ! snowbd/100
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(snow_depth_input_col) == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))

    do c = bounds%begc,bounds%endc
       this%snow_depth_col(c)         = snow_depth_input_col(c)
       this%snow_layer_unity_col(c,:) = 1._r8
    end do

    do c = bounds%begc,bounds%endc
       this%wf_col(c) = spval
       this%wf2_col(c) = spval
    end do


    ! Water Stored in plants is almost always a static entity, with the exception
    ! of when FATES-hydraulics is used. As such, this is trivially set to 0.0 (rgk 03-2017)
    this%total_plant_stored_h2o_col(bounds%begc:bounds%endc) = 0.0_r8


    associate(snl => col%snl) 

      this%snounload_patch(bounds%begp:bounds%endp) = 0._r8
      this%frac_h2osfc_col(bounds%begc:bounds%endc) = 0._r8

      this%fwet_patch(bounds%begp:bounds%endp) = 0._r8
      this%fdry_patch(bounds%begp:bounds%endp) = 0._r8
      this%fcansno_patch(bounds%begp:bounds%endp) = 0._r8
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
               snowbd   = min(400._r8, waterstatebulk_inst%h2osno_col(c)/this%snow_depth_col(c)) !bulk density of snow (kg/m3)
               fmelt    = (snowbd/100.)**1.
               ! 100 is the assumed fresh snow density; 1 is a melting factor that could be
               ! reconsidered, optimal value of 1.5 in Niu et al., 2007
               this%frac_sno_col(c) = tanh( this%snow_depth_col(c) /(2.5 * zlnd * fmelt) )
            endif
         end if
      end do

      do c = bounds%begc,bounds%endc
         if (snl(c) < 0) then
            this%snw_rds_col(c,snl(c)+1:0)        = snw_rds_min
            this%snw_rds_col(c,-nlevsno+1:snl(c)) = 0._r8
            this%snw_rds_top_col(c)               = snw_rds_min
         elseif (waterstatebulk_inst%h2osno_col(c) > 0._r8) then
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
  subroutine RestartBulk(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod          , only : masterproc
    use clm_varcon       , only : pondmx, watmin, spval, nameg
    use landunit_varcon  , only : istcrop, istdlak, istsoil  
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager , only : is_first_step
    use clm_varctl       , only : bound_h2osoi
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical  :: readvar
    !------------------------------------------------------------------------


    call this%Restart(bounds, ncid, flag=flag)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('SNOUNLOAD'), &
         xtype=ncd_double,  &
         dim1name='pft', &
         long_name=this%info%lname('Canopy snow unloading'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%snounload_patch)


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
    if (flag == 'read' .and. .not. readvar) then
       this%frac_h2osfc_col(bounds%begc:bounds%endc) = 0.0_r8
    end if

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
  subroutine ResetBulk(this, column)
    !
    ! !DESCRIPTION:
    ! Intitialize SNICAR variables for fresh snow column
    !
    ! !ARGUMENTS:
    class(waterdiagnosticbulk_type) :: this
    integer , intent(in)   :: column     ! column index
    !-----------------------------------------------------------------------

    this%snw_rds_col(column,0)  = snw_rds_min

  end subroutine ResetBulk

end module WaterDiagnosticBulkType
