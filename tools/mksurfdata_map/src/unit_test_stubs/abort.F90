subroutine abort()
    ! Replacement for abort that throws a pfunit exception rather than aborting
    !
    ! This can be used to test expected errors (i.e., failure testing).
    !
    ! If this occurs within a pFUnit-based test:
    !
    ! - If you have code like:
    !
    !   @assertExceptionRaised("ABORTED:")
    !
    ! - If you don't have
    !
    !   @assertExceptionRaised
    !
    !   or
    !
    !   call assertExceptionRaised
    !
    !   then this will result in the given pFUnit test failing.
    use pfunit_mod, only : throw
    implicit none

    call throw("ABORTED:")
end subroutine abort
