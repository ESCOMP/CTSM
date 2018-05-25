module WaterStateBulkType

  ! This type contains everything in waterstate_type, plus water state variables that
  ! only apply to the bulk water state.
  type, extends(waterstate_type), public :: waterstate_bulk_type
     real(r8), pointer :: snow_persistence_col   (:)   ! col length of time that ground has had non-zero snow thickness (sec)
     real(r8), pointer :: int_snow_col           (:)   ! col integrated snowfall (mm H2O)
   contains
     procedure :: Init
     procedure :: Restart
  end type waterstate_bulk_type

contains

  subroutine Init(this, bounds, ...)
    call this%waterstate_type%Init(bounds, ...)

    ! Then do init stuff needed for the variables declared in this module
  end subroutine Init

  subroutine Restart(this, bounds, ...)
    call this%waterstate_type%Restart(bounds, ...)

    ! Then do restart stuff needed for the variables declared in this module
  end subroutine Restart

end module WaterStateBulkType
