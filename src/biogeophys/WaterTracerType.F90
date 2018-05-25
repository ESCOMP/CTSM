module WaterTracerType

  ! This type groups together the types that apply to water tracers. I like that this
  ! groups together instances relating to the same tracer, and that it potentially
  ! simplifies argument lists.
  !
  ! Mat: You don't need to implement this now: This will come in the next iteration
  type, public :: watertracer_type
     type(waterstate_type) :: waterstate_inst
     type(waterdiag_type) :: waterdiag_inst
     type(waterbalance_type) :: waterbalance_inst
   contains
     procedure :: Init
     procedure :: Restart
  end type watertracer_type

contains

  subroutine Init(this, bounds, ...)
    call this%waterstate_inst%Init(bounds, ...)
    call this%waterdiag_inst%Init(bounds, ...)
    call this%waterbalance_inst%Init(bounds, ...)
  end subroutine Init

  subroutine Restart(this, bounds, ...)
    call this%waterstate_inst%Restart(bounds, ...)
    call this%waterdiag_inst%Restart(bounds, ...)
    call this%waterbalance_inst%Restart(bounds, ...)
  end subroutine Restart

end module WaterTracerType
