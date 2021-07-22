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
     module procedure endrun_write_point_context
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
  subroutine endrun_write_point_context(subgrid_index, subgrid_level, msg, additional_msg)

    !-----------------------------------------------------------------------
    ! Description:
    ! Abort the model for abnormal termination
    !
    ! This version also prints additional information about the point causing the error.
    !
    use shr_sys_mod       , only: shr_sys_abort
    use clm_varctl        , only: iulog
    use GetGlobalValuesMod, only: write_point_context
    !
    ! Arguments:
    implicit none
    integer          , intent(in)           :: subgrid_index  ! index of interest (can be at any subgrid level or gridcell level)
    character(len=*) , intent(in)           :: subgrid_level  ! one of nameg, namel, namec or namep

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

    call write_point_context(subgrid_index, subgrid_level)

    if (present (additional_msg)) then
       write(iulog,*)'ENDRUN: ', additional_msg
    else
       write(iulog,*)'ENDRUN:'
    end if

    call shr_sys_abort(msg)

  end subroutine endrun_write_point_context

end module abortutils
