module WaterBulkAndTracerType

  ! This type groups together the types that apply to both bulk water and tracers. I like
  ! that this groups together instances relating to the same tracer, and that it
  ! potentially simplifies argument lists.
  type, public :: water_bulk_and_tracer_type
     type(water_state_type) :: water_state_inst
     type(water_diag_type) :: water_diag_inst
     type(water_balance_type) :: water_balance_inst
   contains
     procedure :: Init
     procedure :: Restart
  end type water_bulk_and_tracer_type

contains

  subroutine Init(this, bounds, ...)
    call this%water_state_inst%Init(bounds, ...)
    call this%water_diag_inst%Init(bounds, ...)
    call this%water_balance_inst%Init(bounds, ...)
  end subroutine Init

  subroutine Restart(this, bounds, ...)
    call this%water_state_inst%Restart(bounds, ...)
    call this%water_diag_inst%Restart(bounds, ...)
    call this%water_balance_inst%Restart(bounds, ...)
  end subroutine Restart

end module WaterBulkAndTracerType
