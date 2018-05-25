module WaterType

  ! This is the only type that will appear directly in clm_instMod: clm_inst will have a
  ! single instance of water_type, named water_inst, and will call water_inst%Init and
  ! water_inst%Restart.
  type, public :: water_type
     type(waterstate_type) :: waterstate_bulk_inst
     type(waterdiag_type)  :: waterdiag_bulk_inst
     type(waterbalance_type) :: waterbalance_bulk_inst
     type(waterstate_type), allocatable :: waterstate_tracer_inst(:)
     type(waterdiag_type), allocatable :: waterdiag_tracer_inst(:)
     type(waterbalance_type), allocatable :: waterbalance_tracer_inst(:)
     type(water_bulk_only_state_type) :: water_bulk_only_state_inst
     type(water_bulk_only_diag_type) :: water_bulk_only_diag_inst

     ! Probably move fluxes in here, too - eventually
   contains
     procedure :: Init
     procedure :: Restart
  end type water_type

contains

  subroutine Init(this, bounds, ...)
    call this%waterstate_bulk_inst%Init(bounds, ...)
    call this%waterdiag_bulk_inst%Init(bounds, ...)
    call this%waterbalance_bulk_inst%Init(bounds, ...)

    ! Handle tracer instances here

    call this%water_bulk_only_state_inst%Init(bounds, ...)
    call this%water_bulk_only_diag_inst%Init(bounds, ...)
  end subroutine Init

  ! And similarly for Restart

end module WaterType

! Calls from the driver and associated subroutine headers could then look like one of
! these:

! (a)

call canopy_hydrology(..., water_inst, ...)

subroutine canopy_hydrology(..., water_inst, ...)
  type(water_type), intent(inout) :: water_inst
  associate( &
       foo => water_inst%waterstate_bulk_inst%foo, &  ! Input
       bar => water_inst%waterdiag_bulk_inst%bar   &  ! Output
       )
  end associate
end subroutine canopy_hydrology

! (b)

call canopy_hydrology(..., water_inst%waterstate_bulk_inst, water_inst%waterdiag_bulk_inst, ...)

subroutine canopy_hydrology(..., waterstate_bulk_inst, waterdiag_bulk_inst, ...)
  ! Note that this lets us specify intent separately for state vs. diag
  type(waterstate_type), intent(in) :: waterstate_bulk_inst
  type(waterdiag_type), intent(inout) :: waterdiag_bulk_inst
  associate( &
       foo => waterstate_bulk_inst%foo, &  ! Input
       bar => waterdiag_bulk_inst%bar   &  ! Output
       )
  end associate
end subroutine canopy_hydrology
