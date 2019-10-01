module shr_utils_mod

  use shr_sys_mod, only : shr_sys_abort
  implicit none
  private

  public :: shr_utils_ChkErr

  character(*), parameter :: u_FILE_u = __FILE__

!=========================================================
contains
!=========================================================

  logical function shr_utils_ChkErr(rc, line, file, mpierr)

    use mpi , only : MPI_ERROR_STRING, MPI_MAX_ERROR_STRING, MPI_SUCCESS
    use ESMF, only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO
    use ESMF, only : ESMF_FAILURE, ESMF_LogWrite

    ! input/output arguments
    integer           , intent(in) :: rc
    integer           , intent(in) :: line
    character(len=*)  , intent(in) :: file
    logical, optional , intent(in) :: mpierr

    ! local variables
    character(MPI_MAX_ERROR_STRING) :: lstring
    integer :: dbrc, lrc, len, ierr
    !------------------------------------------

    shr_utils_ChkErr = .false.
    lrc = rc
    if (present(mpierr) .and. mpierr) then
       if (rc == MPI_SUCCESS) return
       call MPI_ERROR_STRING(rc, lstring, len, ierr)
       call ESMF_LogWrite("ERROR: "//trim(lstring), ESMF_LOGMSG_INFO, line=line, file=file, rc=dbrc)
       lrc = ESMF_FAILURE
    endif

    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       shr_utils_ChkErr = .true.
    endif

  end function shr_utils_ChkErr

end module shr_utils_mod
