module shr_sys_mod

  use shr_kind_mod  ! defines real & integer kinds

  implicit none
  private

#include <mpif.h>

  public :: shr_sys_abort   ! abort a program

!===============================================================================
CONTAINS
!===============================================================================

  subroutine shr_sys_abort(message, file, line)

    ! Parallel emergency stop

    ! input/output variables
    character(len=*), optional, intent(in) :: message
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: line

    ! Local variables
    integer :: rc
    integer :: ier
    character(len=SHR_KIND_CL):: abort_msg

    if (.not. present(message)) then
       rc = 1001
       call mpi_abort(MPI_COMM_WORLD, rc, ier)
    else
       if (present(file) .and. present(line)) then
          write(abort_msg, '(4a,i0)') trim(message),' at ',trim(file),':',line
       else if (present(file)) then
          write(abort_msg, '(3a)') trim(message),' at ',trim(file)
       else if (present(line)) then
          write(abort_msg, '(2a,i0)') trim(message),' on line ',line
       else
          write(abort_msg, '(a)') trim(message)
       end if

       write(6,*) trim(message)
       rc = 1001
       call mpi_abort(MPI_COMM_WORLD, rc, ier)
    end if

  end subroutine shr_sys_abort

end module shr_sys_mod
