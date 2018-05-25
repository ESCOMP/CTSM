module WaterBulkOnlyStateType

  ! This type contains water state variables that only apply to the bulk water state. I
  ! don't love the name.
  type, public :: water_bulk_only_state_type
     real(r8), pointer :: snow_persistence_col   (:)   ! col length of time that ground has had non-zero snow thickness (sec)
     real(r8), pointer :: int_snow_col           (:)   ! col integrated snowfall (mm H2O)
   contains
     procedure :: Init
     procedure :: Restart
  end type water_bulk_only_state_type

contains

  ! Standard infrastructure routines here

end module WaterBulkOnlyStateType
