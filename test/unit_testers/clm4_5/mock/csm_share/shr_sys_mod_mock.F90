module shr_sys_mod

  ! This is a mock replacement for shr_sys_mod, which just contains an abort routine
  ! without any dependencies

  implicit none

  public :: shr_sys_abort

contains

  subroutine shr_sys_abort(string,rc)

    implicit none

    character(*) ,optional :: string  ! error message string
    integer      ,optional :: rc      ! error code

    !----- local -----
    integer :: ierr
    logical :: flag

    !----- formats -----
    character(*),parameter :: subName =   '(shr_sys_abort) '
    character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"

    !-------------------------------------------------------------------------------
    ! PURPOSE: consistent stopping mechanism
    !-------------------------------------------------------------------------------

    if (present(string)) then
       if (len_trim(string) > 0) write(*,F00) 'ERROR: '//trim(string)
    end if

    stop

  end subroutine shr_sys_abort

end module shr_sys_mod
