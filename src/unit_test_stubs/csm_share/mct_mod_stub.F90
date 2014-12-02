module mct_mod

  ! This is a stub of mct_mod, which only includes the bare minimum needed to build CLM
  ! unit tests

  implicit none

  public :: mct_gsMap
  public :: mct_gsMap_OP

  type mct_gsMap
     ! Empty, dummy type
  end type mct_gsMap

contains

  subroutine mct_gsMap_OP(GSMap, PEno, Points)
    ! Stub routine that simply matches the signature of mct_gsMap_OP
    ! this routine allocates the Points array, to match the documented behavior of the
    ! real routine. This is needed so that a later deallocate will succeed. But note that
    ! it is just allocated to be of size 1, so it cannot be used for any real
    ! calculations.
    type(mct_gsMap), intent(in) :: GSMap
    integer, intent(in) :: PEno
    integer,dimension(:),pointer :: Points

    allocate(Points(1))
  end subroutine mct_gsMap_OP

end module mct_mod
