module lilac_constants

  use shr_kind_mod, only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8

  implicit none
  public

  logical,  parameter :: lilac_constants_statewrite_flag = .false.
  real(R8), parameter :: lilac_constants_spval_init      = 0.0_R8  ! spval for initialization
  real(R8), parameter :: lilac_constants_spval           = 0.0_R8  ! spval
  real(R8), parameter :: lilac_constants_czero           = 0.0_R8  ! spval
  integer,  parameter :: lilac_constants_ispval_mask     = -987987 ! spval for RH mask values
  integer,  parameter :: lilac_constants_SecPerDay       = 86400   ! Seconds per day
  integer             :: lilac_constants_dbug_flag       = 0

end module lilac_constants
