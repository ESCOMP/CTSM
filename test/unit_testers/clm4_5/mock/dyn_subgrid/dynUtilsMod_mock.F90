module dynUtilsMod

  ! This is a mock replacement for dynUtilsMod. It mocks out dyn_time_weights, returning
  ! a fixed weight rather than going through the clm_time_manager (which would be awkward
  ! from a unit test).

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

contains

  subroutine dyn_time_weights(wt1, offset)
    ! Return a fixed weight, rather than going through the clm_time_manager, which would
    ! be awkward from a unit test

    real(r8) , intent(out)          :: wt1
    integer  , intent(in), optional :: offset
    
    real(r8), parameter :: fixed_weight = 0.25_r8

    wt1 = fixed_weight
  end subroutine dyn_time_weights
    
end module dynUtilsMod
