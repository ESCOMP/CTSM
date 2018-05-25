module WaterBulkType

  ! This type groups together the types that only apply to bulk water. This isn't really
  ! necessary, but it simplifies some argument passing, and keeps things symmetrical with
  ! watertracer_type.
  type, public :: waterbulk_type
     type(waterstate_bulk_type) :: waterstate_inst
     type(waterdiag_bulk_type) :: waterdiag_inst
     type(waterbalance_type) :: waterbalance_inst
   contains
     procedure :: Init
     procedure :: Restart
  end type waterbulk_type

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

end module WaterBulkType
