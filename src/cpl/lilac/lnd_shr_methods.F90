module lnd_shr_methods

  use ESMF         , only : ESMF_State, ESMF_MAXSTR, ESMF_StateGet
  use ESMF         , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_FAILURE
  use ESMF         , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError
  use ESMF         , only : ESMF_SUCCESS, ESMF_Field, ESMF_LOGMSG_INFO
  use shr_kind_mod , only : r8 => shr_kind_r8, cs=>shr_kind_cs
  use shr_sys_mod  , only : shr_sys_abort

  implicit none
  private

  public  :: state_diagnose
  public  :: chkerr

  private :: field_getfldptr

  ! Module data
  character(len=1024)         :: msgString
  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine state_diagnose(State, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    integer                         :: i,j,n
    type(ESMf_Field)                :: lfield
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(r8), pointer               :: dataPtr1d(:)
    real(r8), pointer               :: dataPtr2d(:,:)
    character(len=*),parameter      :: subname='(state_diagnose)'
    ! ----------------------------------------------

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    call ESMF_StateGet(state, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount

       call ESMF_StateGet(state, itemName=lfieldnamelist(n), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call field_getfldptr(lfield, rc=rc, fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data
       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
          endif
       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
          endif
       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    deallocate(lfieldnamelist)

  end subroutine state_diagnose

!===============================================================================

  subroutine field_getfldptr(field, rc, fldptr1, fldptr2, rank, abort)

    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------
    use ESMF          , only : ESMF_FieldStatus_Flag, ESMF_Mesh
    use ESMF          , only : ESMF_FieldGet, ESMF_GEOMTYPE_GRID, ESMF_GEOMTYPE_MESH
    use ESMF          , only : ESMF_MeshGet, ESMF_GeomType_Flag
    use ESMF          , only : ESMF_FIELDSTATUS_COMPLETE
    use ESMF          , only : operator(/=), operator(==)

    ! input/output variables
    type(ESMF_Field)  , intent(in)              :: field
    integer           , intent(out)             :: rc
    real(r8), pointer , intent(inout), optional :: fldptr1(:)
    real(r8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort

    ! local variables
    type(ESMF_GeomType_Flag)    :: geomtype
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Mesh)             :: lmesh
    integer                     :: lrank, nnodes, nelements
    logical                     :: labort
    character(len=*), parameter :: subname='(field_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    labort = .true.
    if (present(abort)) then
       labort = abort
    endif
    lrank = -99

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       lrank = 0
       if (labort) then
          call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       else
          call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       endif
    else

       call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (geomtype == ESMF_GEOMTYPE_GRID) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (geomtype == ESMF_GEOMTYPE_MESH) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (nnodes == 0 .and. nelements == 0) lrank = 0
       else  
          call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", &
               ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       endif ! geomtype

       if (lrank == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
               ESMF_LOGMSG_INFO)
       elseif (lrank == 1) then
          if (.not.present(fldptr1)) then
             call ESMF_LogWrite(trim(subname)//": ERROR missing rank=1 array ", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (lrank == 2) then
          if (.not.present(fldptr2)) then
             call ESMF_LogWrite(trim(subname)//": ERROR missing rank=2 array ", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//": ERROR in rank ", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       endif

    endif  ! status

    if (present(rank)) then
       rank = lrank
    endif

  end subroutine field_getfldptr

!===============================================================================

  logical function chkerr(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

end module lnd_shr_methods
