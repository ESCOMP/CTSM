module shr_sys_mod

  use shr_kind_mod  ! defines real & integer kinds

  implicit none
  private

#include <mpif.h>

  public :: shr_sys_getenv  ! get an environment variable
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

  !===============================================================================
  subroutine shr_sys_getenv(name, val, rcode)

    !-------------------------------------------------------------------------------
    ! PURPOSE: an architecture independant system call
    !-------------------------------------------------------------------------------

    ! input/output variables
    character(len=*) ,intent(in)     :: name    ! env var name
    character(len=*) ,intent(inout)  :: val     ! env var value
    integer          ,intent(inout)  :: rcode   ! return code

    !----- local -----
    integer                :: lenname ! length of env var name
    integer                :: lenval  ! length of env var value
    character(SHR_KIND_CL) :: tmpval  ! temporary env var value
    character(*),parameter :: subName =   '(shr_sys_getenv) '
    character(*),parameter :: F00     = "('(shr_sys_getenv) ',4a)"

    lenname=len_trim(name)

#if (defined IRIX64 || defined CRAY || defined UNICOSMP)

    call pxfgetenv(name, lenname, val, lenval, rcode)

#elif (defined AIX || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)

    call getenv(trim(name),tmpval)
    val=trim(tmpval)
    rcode = 0
    if (len_trim(val) ==  0         ) rcode = 1
    if (len_trim(val) >  SHR_KIND_CL) rcode = 2

#else

    write(6,F00) 'ERROR: no implementation of getenv for this architecture'
    call shr_sys_abort(subname//'no implementation of getenv for this machine')

#endif

  end subroutine shr_sys_getenv

end module shr_sys_mod
