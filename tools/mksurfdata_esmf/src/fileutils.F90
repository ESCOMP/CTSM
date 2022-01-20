module fileutils

  !-----------------------------------------------------------------------
  ! Module containing file I/O utilities
  !-----------------------------------------------------------------------

  use shr_sys_mod, only : shr_sys_abort

  implicit none

  public :: get_filename  !Returns filename given full pathname

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

  character(len=256) function get_filename (fulpath)
    ! Returns filename given full pathname

    ! input/output variables
    character(len=*), intent(in)  :: fulpath !full pathname

    ! local variables:
    integer :: i    !loop index
    integer :: klen !length of fulpath character string
    !------------------------------------------------------------------------

    klen = len_trim(fulpath)
    do i = klen, 1, -1
       if (fulpath(i:i) == '/') go to 10
    end do
    i = 0
10  get_filename = fulpath(i+1:klen)

    return
  end function get_filename

end module fileutils
