module abortutils

  !-----------------------------------------------------------------------
  ! !MODULE: abortutils
  !
  ! !DESCRIPTION:
  ! Abort the model for abnormal termination
  !-----------------------------------------------------------------------

  private
  save

  public :: endrun

  interface endrun
     module procedure endrun_vanilla
     module procedure endrun_globalindex
  end interface

CONTAINS

  !-----------------------------------------------------------------------
  subroutine endrun_vanilla(msg, additional_msg)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    use shr_sys_mod , only: shr_sys_abort
    use clm_varctl  , only: iulog
    !
    ! !ARGUMENTS:
    implicit none

    ! Generally you want to at least provide msg. The main reason to separate msg from
    ! additional_msg is to supported expected-exception unit testing: you can put
    ! volatile stuff in additional_msg, as in:
    !   call endrun(msg='Informative message', additional_msg=errmsg(__FILE__, __LINE__))
    ! and then just assert against msg.
    character(len=*), intent(in), optional :: msg            ! string to be passed to shr_sys_abort
    character(len=*), intent(in), optional :: additional_msg ! string to be printed, but not passed to shr_sys_abort
    !-----------------------------------------------------------------------

    if (present (additional_msg)) then
       write(iulog,*)'ENDRUN: ', trim(additional_msg)
    else
       write(iulog,*)'ENDRUN:'
    end if

    call shr_sys_abort(msg)

  end subroutine endrun_vanilla

  !-----------------------------------------------------------------------
  subroutine endrun_globalindex(decomp_index, clmlevel, msg, additional_msg)

    !-----------------------------------------------------------------------
    ! Description:
    ! Abort the model for abnormal termination
    !
    use shr_sys_mod       , only: shr_sys_abort
    use clm_varctl        , only: iulog
    use GetGlobalValuesMod, only: GetGlobalWrite
    !
    ! Arguments:
    implicit none
    integer          , intent(in)           :: decomp_index
    character(len=*) , intent(in)           :: clmlevel

    ! Generally you want to at least provide msg. The main reason to separate msg from
    ! additional_msg is to supported expected-exception unit testing: you can put
    ! volatile stuff in additional_msg, as in:
    !   call endrun(msg='Informative message', additional_msg=errmsg(__FILE__, __LINE__))
    ! and then just assert against msg.
    character(len=*), intent(in), optional :: msg            ! string to be passed to shr_sys_abort
    character(len=*), intent(in), optional :: additional_msg ! string to be printed, but not passed to shr_sys_abort
    !
    ! Local Variables:
    integer :: igrc, ilun, icol 
    !-----------------------------------------------------------------------

    write(6,*)'calling getglobalwrite with decomp_index= ',decomp_index,' and clmlevel= ',trim(clmlevel)
    call GetGlobalWrite(decomp_index, clmlevel)

    if (present (additional_msg)) then
       write(iulog,*)'ENDRUN: ', additional_msg
    else
       write(iulog,*)'ENDRUN:'
    end if

    call shr_sys_abort(msg)

  end subroutine endrun_globalindex

end module abortutils
