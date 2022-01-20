module shr_sys_mod

  use shr_kind_mod  ! defines real & integer kinds

  implicit none
  private

  public :: shr_sys_getenv  ! get an environment variable
  public :: shr_sys_abort   ! abort a program

  integer :: s_loglev = 0
  integer :: s_logunit = 6

!===============================================================================
CONTAINS
!===============================================================================

  subroutine shr_sys_abort(string,rc)

    !-------------------------------------------------------------------------------
    ! PURPOSE: consistent stopping mechanism
    !-------------------------------------------------------------------------------

    character(len=*) , optional :: string  ! error message string
    integer          , optional :: rc      ! error code

    ! local variables
    character(*),parameter :: subName =   '(shr_sys_abort) '
    character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"

    if (len_trim(string) > 0) write(s_logunit,F00) 'ERROR: '//trim(string)
    write(s_logunit,F00) 'WARNING: stopping'
    call abort()
    stop

  end subroutine shr_sys_abort

  !===============================================================================
  subroutine shr_sys_getenv(name, val, rcode)

    !-------------------------------------------------------------------------------
    ! PURPOSE: an architecture independant system call
    !-------------------------------------------------------------------------------

    ! input/output variables
    character(len=*) ,intent(in)  :: name    ! env var name
    character(len=*) ,intent(out) :: val     ! env var value
    integer          ,intent(out) :: rcode   ! return code

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

    write(s_logunit,F00) 'ERROR: no implementation of getenv for this architecture'
    call shr_sys_abort(subname//'no implementation of getenv for this machine')

#endif

  end subroutine shr_sys_getenv

end module shr_sys_mod
