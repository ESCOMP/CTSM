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
  use clm_varcon     , only : spval
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerUtils, only : CompareBulkToTracer
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterlnd2atm_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: q_ref2m_grc        (:)   ! 2m surface specific humidity (kg/kg)
     real(r8), pointer :: h2osno_grc         (:)   ! snow water (mm H2O)
     real(r8), pointer :: qflx_evap_tot_grc  (:)   ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
     real(r8), pointer :: qflx_rofliq_grc         (:)   ! rof liq forcing
     real(r8), pointer :: qflx_rofliq_qsur_grc    (:)   ! rof liq -- surface runoff component
     real(r8), pointer :: qflx_rofliq_qsub_grc    (:)   ! rof liq -- subsurface runoff component
     real(r8), pointer :: qflx_rofliq_qgwl_grc    (:)   ! rof liq -- glacier, wetland and lakes water balance residual component
     real(r8), pointer :: qflx_rofliq_drain_perched_grc    (:)   ! rof liq -- perched water table runoff component
     real(r8), pointer :: qflx_rofice_grc    (:)   ! rof ice forcing
     real(r8), pointer :: qflx_liq_from_ice_col(:) ! liquid runoff from converted ice runoff
     real(r8), pointer :: qirrig_grc         (:)   ! irrigation flux

   contains

     procedure, public  :: Init
     procedure, public  :: TracerConsistencyCheck
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
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
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
    ! Initialize history vars
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begc, endc
    integer           :: begg, endg
    !------------------------------------------------------------------------

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
    ! Initialize cold start conditions
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(inout) :: this
    type(bounds_type)     , intent(in)      :: bounds
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    ! Nothing to do for now

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine TracerConsistencyCheck(this, bounds, bulk, caller_location)
    ! !DESCRIPTION:
    ! Check consistency of water tracer with that of bulk water
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    class(waterlnd2atm_type), intent(in) :: bulk
    character(len=*), intent(in) :: caller_location  ! brief description of where this is called from (for error messages)
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%q_ref2m_grc(bounds%begg:bounds%endg), &
         tracer=this%q_ref2m_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='q_ref2m_grc')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%h2osno_grc(bounds%begg:bounds%endg), &
         tracer=this%h2osno_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='h2osno_grc')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%qflx_evap_tot_grc(bounds%begg:bounds%endg), &
         tracer=this%qflx_evap_tot_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='qflx_evap_tot_grc')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%qflx_rofliq_grc(bounds%begg:bounds%endg), &
         tracer=this%qflx_rofliq_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='qflx_rofliq_grc')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%qflx_rofliq_qsur_grc(bounds%begg:bounds%endg), &
         tracer=this%qflx_rofliq_qsur_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='qflx_rofliq_qsur_grc')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%qflx_rofliq_qsub_grc(bounds%begg:bounds%endg), &
         tracer=this%qflx_rofliq_qsub_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='qflx_rofliq_qsub_grc')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%qflx_rofliq_qgwl_grc(bounds%begg:bounds%endg), &
         tracer=this%qflx_rofliq_qgwl_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='qflx_rofliq_qgwl_grc')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%qflx_rofliq_drain_perched_grc(bounds%begg:bounds%endg), &
         tracer=this%qflx_rofliq_drain_perched_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='qflx_rofliq_drain_perched_grc')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%qflx_rofice_grc(bounds%begg:bounds%endg), &
         tracer=this%qflx_rofice_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='qflx_rofice_grc')

    call CompareBulkToTracer(bounds%begc, bounds%endc, &
         bulk=bulk%qflx_liq_from_ice_col(bounds%begc:bounds%endc), &
         tracer=this%qflx_liq_from_ice_col(bounds%begc:bounds%endc), &
         caller_location=caller_location, &
         name='qflx_liq_from_ice_col')

    call CompareBulkToTracer(bounds%begg, bounds%endg, &
         bulk=bulk%qirrig_grc(bounds%begg:bounds%endg), &
         tracer=this%qirrig_grc(bounds%begg:bounds%endg), &
         caller_location=caller_location, &
         name='qirrig_grc')

  end subroutine TracerConsistencyCheck

end module Waterlnd2atmType
