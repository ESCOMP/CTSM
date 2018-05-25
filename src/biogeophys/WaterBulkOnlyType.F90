module WaterBulkOnlyType

  ! This type groups together the types that only apply to bulk water. This isn't really
  ! necessary, but it simplifies some argument passing, and keeps things symmetrical with
  ! water_bulk_and_tracer_type.
  type, public :: water_bulk_only_type
     type(water_bulk_only_state_type) :: water_bulk_only_state_inst
     type(water_bulk_only_diag_type) :: water_bulk_only_diag_inst
   contains
     procedure :: Init
     procedure :: Restart
  end type water_bulk_only_type

contains

  subroutine Init(this, bounds, ...)
    call this%water_bulk_only_state_inst%Init(bounds, ...)
    call this%water_bulk_only_diag_inst%Init(bounds, ...)
  end subroutine Init

  subroutine Restart(this, bounds, ...)
    call this%water_bulk_only_state_inst%Restart(bounds, ...)
    call this%water_bulk_only_diag_inst%Restart(bounds, ...)
  end subroutine Restart

end module WaterBulkOnlyType
