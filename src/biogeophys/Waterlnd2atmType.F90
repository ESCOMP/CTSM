module Waterlnd2atmType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water lnd2atm variables that apply to both bulk water
  ! and water tracers.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varctl     , only : iulog
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use WaterInfoBaseType, only : water_info_base_type
  use WaterIsotopesMod, only : WisoCompareBulkToTracer
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterlnd2atm_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: q_ref2m_grc        (:)   => null() ! 2m surface specific humidity (kg/kg)       
     real(r8), pointer :: h2osno_grc         (:)   => null() ! snow water (mm H2O)
     real(r8), pointer :: qflx_evap_tot_grc  (:)   => null() ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg                        
     real(r8), pointer :: qflx_rofliq_grc         (:)   => null() ! rof liq forcing
     real(r8), pointer :: qflx_rofliq_qsur_grc    (:)   => null() ! rof liq -- surface runoff component
     real(r8), pointer :: qflx_rofliq_qsub_grc    (:)   => null() ! rof liq -- subsurface runoff component
     real(r8), pointer :: qflx_rofliq_qgwl_grc    (:)   => null() ! rof liq -- glacier, wetland and lakes water balance residual component
     real(r8), pointer :: qflx_rofliq_drain_perched_grc    (:)   => null() ! rof liq -- perched water table runoff component        
     real(r8), pointer :: qflx_rofice_grc    (:)   => null() ! rof ice forcing
     real(r8), pointer :: qflx_liq_from_ice_col(:) => null() ! liquid runoff from converted ice runoff                              
     real(r8), pointer :: qirrig_grc         (:)   => null() ! irrigation flux                                                      

   contains

     procedure          :: Init         
     procedure          :: Restart      
     procedure          :: TracerConsistencyCheck
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type waterlnd2atm_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info)

    class(waterlnd2atm_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds  
    class(water_info_base_type), intent(in), target :: info

    this%info => info

    call this%InitAllocate(bounds)

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
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value     
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate(this%q_ref2m_grc        (begg:endg))            ; this%q_ref2m_grc        (:)   =ival
    allocate(this%h2osno_grc         (begg:endg))            ; this%h2osno_grc         (:)   =ival
    allocate(this%qflx_evap_tot_grc  (begg:endg))            ; this%qflx_evap_tot_grc  (:)   =ival
    allocate(this%qflx_rofliq_grc    (begg:endg))            ; this%qflx_rofliq_grc    (:)   =ival
    allocate(this%qflx_rofliq_qsur_grc    (begg:endg))       ; this%qflx_rofliq_qsur_grc    (:)   =ival
    allocate(this%qflx_rofliq_qsub_grc    (begg:endg))       ; this%qflx_rofliq_qsub_grc    (:)   =ival
    allocate(this%qflx_rofliq_qgwl_grc    (begg:endg))       ; this%qflx_rofliq_qgwl_grc    (:)   =ival
    allocate(this%qflx_rofliq_drain_perched_grc    (begg:endg))       ; this%qflx_rofliq_drain_perched_grc    (:)   =ival
    allocate(this%qflx_rofice_grc    (begg:endg))            ; this%qflx_rofice_grc    (:)   =ival
    allocate(this%qflx_liq_from_ice_col(begc:endc))          ; this%qflx_liq_from_ice_col(:) = ival
    allocate(this%qirrig_grc         (begg:endg))            ; this%qirrig_grc         (:)   =ival


  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_lch4
    use clm_varctl     , only : hist_wrtch4diag
    use clm_varpar     , only : nlevsno, nlevsoi
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal, no_snow_zero
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(in) :: this
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


    this%qflx_rofliq_grc(begg:endg) = 0._r8
    call hist_addfld1d (fname=this%info%fname('QRUNOFF_TO_COUPLER'),  units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('total liquid runoff sent to coupler (includes corrections for land use change)'), &
         ptr_lnd=this%qflx_rofliq_grc)

    this%qflx_rofice_grc(begg:endg) = 0._r8
    call hist_addfld1d (fname=this%info%fname('QRUNOFF_ICE_TO_COUPLER'),  units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('total ice runoff sent to coupler (includes corrections for land use change)'), &
         ptr_lnd=this%qflx_rofice_grc)

    this%qflx_liq_from_ice_col(begc:endc) = 0._r8
    call hist_addfld1d (fname=this%info%fname('QRUNOFF_ICE_TO_LIQ'), units='mm/s', &
         avgflag='A', &
         long_name=this%info%lname('liquid runoff from converted ice runoff'), &
         ptr_col=this%qflx_liq_from_ice_col, default='inactive')



  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use column_varcon   , only : icol_shadewall, icol_road_perv
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use clm_varcon      , only : denice, denh2o, spval, sb, bdsno 
    use clm_varcon      , only : zlnd, tfrz, spval, pc, aquifer_water_baseline
    use clm_varctl      , only : fsurdat, iulog
    use clm_varctl        , only : use_bedrock
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(in) :: this
    type(bounds_type)     , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev,nlevs 
    real(r8)           :: maxslope, slopemax, minslope
    real(r8)           :: d, fd, dfdd, slope0,slopebeta
    real(r8) ,pointer  :: std (:)     
    logical            :: readvar 
    type(file_desc_t)  :: ncid        
    character(len=256) :: locfn       
    integer            :: nbedrock
    !-----------------------------------------------------------------------



  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod          , only : masterproc
    use clm_varcon       , only : denice, denh2o, pondmx, watmin, spval, nameg
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager , only : is_first_step
    use clm_varctl       , only : bound_h2osoi
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(in) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,j,nlevs
    logical  :: readvar
    real(r8) :: maxwatsat    ! maximum porosity    
    real(r8) :: excess       ! excess volumetric soil water
    real(r8) :: totwat       ! total soil water (mm)
    !------------------------------------------------------------------------


!    call restartvar(ncid=ncid, flag=flag, &
!         varname=this%info%fname('H2OSFC'), &
!         xtype=ncd_double,  &
!         dim1name='column', &
!         long_name=this%info%lname('surface water'), &
!         units='kg/m2', &
!         interpinic_flag='interp', readvar=readvar, data=this%h2osfc_col)
!    if (flag=='read' .and. .not. readvar) then
!       this%h2osfc_col(bounds%begc:bounds%endc) = 0.0_r8
!    end if



  end subroutine Restart

  !------------------------------------------------------------------------
  subroutine TracerConsistencyCheck(this,bounds,tracer)
    ! !DESCRIPTION:
    ! Check consistency of water tracer with that of bulk water
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    class(waterlnd2atm_type), intent(in) :: tracer
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%q_ref2m_grc(bounds%begg:bounds%endg), &
         tracer%q_ref2m_grc(bounds%begg:bounds%endg), &
         'q_ref2m_grc')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%h2osno_grc(bounds%begg:bounds%endg), &
         tracer%h2osno_grc(bounds%begg:bounds%endg), &
         'h2osno_grc')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%qflx_evap_tot_grc(bounds%begg:bounds%endg), &
         tracer%qflx_evap_tot_grc(bounds%begg:bounds%endg), &
         'qflx_evap_tot_grc')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%qflx_rofliq_grc(bounds%begg:bounds%endg), &
         tracer%qflx_rofliq_grc(bounds%begg:bounds%endg), &
         'qflx_rofliq_grc')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%qflx_rofliq_qsur_grc(bounds%begg:bounds%endg), &
         tracer%qflx_rofliq_qsur_grc(bounds%begg:bounds%endg), &
         'qflx_rofliq_qsur_grc')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%qflx_rofliq_qsub_grc(bounds%begg:bounds%endg), &
         tracer%qflx_rofliq_qsub_grc(bounds%begg:bounds%endg), &
         'qflx_rofliq_qsub_grc')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%qflx_rofliq_qgwl_grc(bounds%begg:bounds%endg), &
         tracer%qflx_rofliq_qgwl_grc(bounds%begg:bounds%endg), &
         'qflx_rofliq_qgwl_grc')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%qflx_rofliq_drain_perched_grc(bounds%begg:bounds%endg), &
         tracer%qflx_rofliq_drain_perched_grc(bounds%begg:bounds%endg), &
         'qflx_rofliq_drain_perched_grc')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%qflx_rofice_grc(bounds%begg:bounds%endg), &
         tracer%qflx_rofice_grc(bounds%begg:bounds%endg), &
         'qflx_rofice_grc')

    call WisoCompareBulkToTracer(bounds%begc, bounds%endc, &
         this%qflx_liq_from_ice_col(bounds%begc:bounds%endc), &
         tracer%qflx_liq_from_ice_col(bounds%begc:bounds%endc), &
         'qflx_liq_from_ice_col')

    call WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
         this%qirrig_grc(bounds%begg:bounds%endg), &
         tracer%qirrig_grc(bounds%begg:bounds%endg), &
         'qirrig_grc')

  end subroutine TracerConsistencyCheck

end module Waterlnd2atmType
