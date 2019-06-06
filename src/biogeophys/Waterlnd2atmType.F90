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
  use decompMod      , only : BOUNDS_SUBGRID_COLUMN, BOUNDS_SUBGRID_GRIDCELL
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  use WaterTracerUtils, only : AllocateVar1d
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
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type waterlnd2atm_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info, tracer_vars)

    class(waterlnd2atm_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: tracer_vars

    this%info => info

    call this%InitAllocate(bounds, tracer_vars)

    call this%InitHistory(bounds)

    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, tracer_vars)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterlnd2atm_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(water_tracer_container_type), intent(inout) :: tracer_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    !------------------------------------------------------------------------

    call AllocateVar1d(var = this%q_ref2m_grc, name = 'q_ref2m_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%h2osno_grc, name = 'h2osno_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%qflx_evap_tot_grc, name = 'qflx_evap_tot_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%qflx_rofliq_grc, name = 'qflx_rofliq_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%qflx_rofliq_qsur_grc, name = 'qflx_rofliq_qsur_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%qflx_rofliq_qsub_grc, name = 'qflx_rofliq_qsub_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%qflx_rofliq_qgwl_grc, name = 'qflx_rofliq_qgwl_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%qflx_rofliq_drain_perched_grc, name = 'qflx_rofliq_drain_perched_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%qflx_rofice_grc, name = 'qflx_rofice_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%qflx_liq_from_ice_col, name = 'qflx_liq_from_ice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN, &
         ival=ival)
    call AllocateVar1d(var = this%qirrig_grc, name = 'qirrig_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)

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

end module Waterlnd2atmType
