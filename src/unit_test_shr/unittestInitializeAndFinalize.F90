module unittestInitializeAndFinalize

  ! Subroutines to do per-executable initialization and finalization (i.e., things that
  ! need to happen once per executable, not once per test)

  implicit none
  private

  public :: unittest_finalize_esmf ! finalization routine for any tests that might have initialized ESMF

contains

  !-----------------------------------------------------------------------
  subroutine unittest_finalize_esmf()
    !
    ! !DESCRIPTION:
    ! Finalization routine for any tests that might have initialized ESMF
    !
    ! A good rule is: for any pFUnit test that links against esmf (i.e., includes esmf in
    ! the LINK_LIBRARIES line in the add_pfunit_ctest call), this finalization routine
    ! should be called. (It is safe to call this even if a test didn't initialize ESMF:
    ! the finalization is done conditionally.)
    !
    ! This can be called by adding the following to the add_pfunit_ctest call:
    !   EXTRA_FINALIZE unittest_finalize_esmf
    !   EXTRA_USE unittestInitializeAndFinalize
    !
    ! !USES:
    use ESMF, only : ESMF_SUCCESS, ESMF_Finalize, ESMF_IsInitialized
    !
    ! !LOCAL VARIABLES:
    logical :: esmf_is_initialized
    integer :: rc

    character(len=*), parameter :: subname = 'unittest_finalize_esmf'
    !-----------------------------------------------------------------------

    esmf_is_initialized = ESMF_IsInitialized(rc=rc)
    if (rc /= ESMF_SUCCESS) then
       ! Calling shr_sys_abort from this finalization routine leads to a pFUnit test
       ! result of pass instead of fail. Doing a 'stop 1' leads to a failure, as desired.
       print *, 'Error in ESMF_IsInitialized'
       stop 1
    end if
    if (esmf_is_initialized) then
       print *, 'Finalizing ESMF'
       call ESMF_Finalize(rc=rc)
       if (rc /= ESMF_SUCCESS) then
          print *, 'Error in ESMF_Finalize'
          stop 1
       end if
    end if
  end subroutine unittest_finalize_esmf

end module unittestInitializeAndFinalize