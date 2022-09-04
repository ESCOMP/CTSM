module rof_comp_esmf

  ! ------------------------------------------------------------------------
  ! This is a stub version of rof_comp_esmf that can be used when we don't have a true
  ! rof component, just to satisfy the necessary interfaces in LILAC.
  ! ------------------------------------------------------------------------

  use ESMF

  implicit none
  private

  public :: rof_register

!===============================================================================
contains
!===============================================================================

  subroutine rof_register(comp, rc)

    ! Stub rof_register routine - shouldn't ever be called!

    ! input/output argumenents
    type(ESMF_GridComp)  :: comp  ! ROF grid component
    integer, intent(out) :: rc    ! return status

    rc = ESMF_RC_NOT_IMPL
  end subroutine rof_register

end module rof_comp_esmf
