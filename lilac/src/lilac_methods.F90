module lilac_methods

  !-----------------------------------------------------------------------------
  ! Generic operation methods used by the Mediator Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use mpi             , only : MPI_ERROR_STRING, MPI_MAX_ERROR_STRING, MPI_SUCCESS
  use shr_kind_mod    , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use lilac_constants , only : dbug_flag => lilac_constants_dbug_flag
  use lilac_constants , only : czero => lilac_constants_czero

  implicit none
  private

  interface lilac_methods_FB_accum ; module procedure &
    lilac_methods_FB_accumFB2FB
  end interface

  interface lilac_methods_FB_copy ; module procedure &
    lilac_methods_FB_copyFB2FB
  end interface

  interface lilac_methods_FieldPtr_compare ; module procedure &
    lilac_methods_FieldPtr_compare1, &
    lilac_methods_FieldPtr_compare2
  end interface

  ! used/reused in module

  logical                               :: isPresent
  character(len=1024)                   :: msgString
  type(ESMF_GeomType_Flag)              :: geomtype
  type(ESMF_FieldStatus_Flag)           :: status
  character(*)      , parameter         :: u_FILE_u = &
       __FILE__

  public lilac_methods_FB_copy
  public lilac_methods_FB_accum
  public lilac_methods_FB_average
  public lilac_methods_FB_reset
  public lilac_methods_FB_clean
  public lilac_methods_FB_diagnose
  public lilac_methods_FB_FldChk
  public lilac_methods_FB_GetFldPtr
  public lilac_methods_FB_getNameN
  public lilac_methods_FB_getFieldN
  public lilac_methods_FB_getFieldByName
  public lilac_methods_FB_getNumflds
  public lilac_methods_FB_Field_diagnose
  public lilac_methods_State_diagnose
  public lilac_methods_State_GetFldPtr
  public lilac_methods_State_SetScalar
  public lilac_methods_State_GetScalar
  public lilac_methods_Clock_TimePrint
  public lilac_methods_FieldPtr_compare
  public chkerr

  private lilac_methods_Mesh_Print
  private lilac_methods_Mesh_Write
  private lilac_methods_Field_GetFldPtr
  private lilac_methods_FB_SetFldPtr
  private lilac_methods_FB_copyFB2FB
  private lilac_methods_FB_accumFB2FB
  private lilac_methods_State_getNameN
  private lilac_methods_State_SetFldPtr

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_getNameN(FB, fieldnum, fieldname, rc)

    ! ----------------------------------------------
    ! Get name of field number fieldnum in input field bundle FB
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    character(len=*)      , intent(out)   :: fieldname
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*),parameter      :: subname='(lilac_methods_FB_getNameN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    fieldname = ' '

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (fieldnum > fieldCount) then
      call ESMF_LogWrite(trim(subname)//": ERROR fieldnum > fieldCount ", ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
    endif

    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    fieldname = lfieldnamelist(fieldnum)

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_getNameN

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_getFieldN(FB, fieldnum, field, rc)

    ! ----------------------------------------------
    ! Get field with number fieldnum in input field bundle FB
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    type(ESMF_Field)      , intent(inout) :: field
    integer               , intent(out)   :: rc

    ! local variables
    character(len=ESMF_MAXSTR) :: name
    character(len=*),parameter :: subname='(lilac_methods_FB_getFieldN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call lilac_methods_FB_getNameN(FB, fieldnum, name, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(FB, fieldName=name, field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_getFieldN

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_getFieldByName(FB, fieldname, field, rc)

    ! ----------------------------------------------
    ! Get field associated with fieldname out of FB
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    character(len=*)      , intent(in)    :: fieldname
    type(ESMF_Field)      , intent(inout) :: field
    integer               , intent(out)   :: rc

    ! local variables
    character(len=*),parameter :: subname='(lilac_methods_FB_getFieldByName)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldName=fieldname, field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_getFieldByName

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_State_getNameN(State, fieldnum, fieldname, rc)

    ! ----------------------------------------------
    ! Get field number fieldnum name out of State
    ! ----------------------------------------------

    type(ESMF_State), intent(in)    :: State
    integer         , intent(in)    :: fieldnum
    character(len=*), intent(out)   :: fieldname
    integer         , intent(out)   :: rc

    ! local variables
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*),parameter      :: subname='(lilac_methods_State_getNameN)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    fieldname = ' '

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (fieldnum > fieldCount) then
      call ESMF_LogWrite(trim(subname)//": ERROR fieldnum > fieldCount ", ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
    endif

    allocate(lfieldnamelist(fieldCount))
    call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    fieldname = lfieldnamelist(fieldnum)

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_State_getNameN

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_clean(FB, rc)

    ! ----------------------------------------------
    ! Destroy fields in FB and FB
    ! ----------------------------------------------

    type(ESMF_FieldBundle), intent(inout) :: FB
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    type(ESMF_Field)                :: field
    character(len=*),parameter      :: subname='(lilac_methods_FB_clean)'
    ! ----------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
      call ESMF_FieldBundleGet(FB, fieldName=lfieldnamelist(n), field=field, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldDestroy(field, rc=rc, noGarbage=.true.)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    enddo

    call ESMF_FieldBundleDestroy(FB, rc=rc, noGarbage=.true.)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    deallocate(lfieldnamelist)
    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine lilac_methods_FB_clean

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_reset(FB, value, rc)
    ! ----------------------------------------------
    ! Set all fields to value in FB
    ! If value is not provided, reset to 0.0
    ! ----------------------------------------------

    ! intput/output variables
    type(ESMF_FieldBundle), intent(inout)        :: FB
    real(R8)    , intent(in), optional :: value
    integer               , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(R8)              :: lvalue
    character(len=*),parameter      :: subname='(lilac_methods_FB_reset)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lvalue = czero
    if (present(value)) then
      lvalue = value
    endif

    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
      call lilac_methods_FB_SetFldPtr(FB, lfieldnamelist(n), lvalue, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_reset

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_average(FB, count, rc)

    ! ----------------------------------------------
    ! Set all fields to zero in FB
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout) :: FB
    integer               , intent(in)    :: count
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(R8), pointer               :: dataPtr1(:)
    real(R8), pointer               :: dataPtr2(:,:)
    character(len=*),parameter      :: subname='(lilac_methods_FB_average)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    if (count == 0) then

       if (dbug_flag > 10) then
          call ESMF_LogWrite(trim(subname)//": WARNING count is 0", ESMF_LOGMSG_INFO)
       end if
       !call ESMF_LogWrite(trim(subname)//": WARNING count is 0 set avg to spval", ESMF_LOGMSG_INFO)
       !call lilac_methods_FB_reset(FB, value=spval, rc=rc)
       !if (chkerr(rc,__LINE__,u_FILE_u)) return

    else

      call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldnamelist(fieldCount))
      call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      do n = 1, fieldCount
        call lilac_methods_FB_GetFldPtr(FB, lfieldnamelist(n), dataPtr1, dataPtr2, lrank, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        if (lrank == 0) then
          ! no local data
        elseif (lrank == 1) then
          do i=lbound(dataptr1,1),ubound(dataptr1,1)
            dataptr1(i) = dataptr1(i) / real(count, R8)
          enddo
        elseif (lrank == 2) then
          do j=lbound(dataptr2,2),ubound(dataptr2,2)
          do i=lbound(dataptr2,1),ubound(dataptr2,1)
            dataptr2(i,j) = dataptr2(i,j) / real(count, R8)
          enddo
          enddo
        else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
        endif
      enddo
      deallocate(lfieldnamelist)

    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_average

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_diagnose(FB, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of FB
    ! ----------------------------------------------

    type(ESMF_FieldBundle) , intent(inout)        :: FB
    character(len=*)       , intent(in), optional :: string
    integer                , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR), pointer :: lfieldnamelist(:)
    character(len=CL)               :: lstring
    real(R8), pointer               :: dataPtr1d(:)
    real(R8), pointer               :: dataPtr2d(:,:)
    character(len=*), parameter     :: subname='(lilac_methods_FB_diagnose)'
    ! ----------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string) // ' '
    endif

    ! Determine number of fields in field bundle and allocate memory for lfieldnamelist
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    ! Get the fields in the field bundle
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! For each field in the bundle, get its memory location and print out the field
    do n = 1, fieldCount
       call lilac_methods_FB_GetFldPtr(FB, lfieldnamelist(n), &
            fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data

       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), " no data"
          endif

       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    ! Deallocate memory
    deallocate(lfieldnamelist)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine lilac_methods_FB_diagnose

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_State_diagnose(State, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    type(ESMF_State), intent(in)           :: State
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=CS)               :: lstring
    real(R8), pointer               :: dataPtr1d(:)
    real(R8), pointer               :: dataPtr2d(:,:)
    character(len=*),parameter      :: subname='(lilac_methods_State_diagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    endif

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    call ESMF_StateGet(State, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount

       call lilac_methods_State_GetFldPtr(State, lfieldnamelist(n), &
            fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data

       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif

       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
   endif

  end subroutine lilac_methods_State_diagnose

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_Field_diagnose(FB, fieldname, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout)  :: FB
    character(len=*), intent(in)           :: fieldname
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    integer           :: lrank
    character(len=CS) :: lstring
    real(R8), pointer :: dataPtr1d(:)
    real(R8), pointer :: dataPtr2d(:,:)
    character(len=*),parameter      :: subname='(lilac_methods_FB_FieldDiagnose)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call lilac_methods_FB_GetFldPtr(FB, fieldname, dataPtr1d, dataPtr2d, lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
       ! no local data
    elseif (lrank == 1) then
       if (size(dataPtr1d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    elseif (lrank == 2) then
       if (size(dataPtr2d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    else
       call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    endif
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_Field_diagnose
  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_copyFB2FB(FBout, FBin, rc)

    ! ----------------------------------------------
    ! Copy common field names from FBin to FBout
    ! ----------------------------------------------

    type(ESMF_FieldBundle), intent(inout) :: FBout
    type(ESMF_FieldBundle), intent(in)    :: FBin
    integer               , intent(out)   :: rc
    character(len=*), parameter :: subname='(lilac_methods_FB_copyFB2FB)'
    ! ----------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    call lilac_methods_FB_accum(FBout, FBin, copy=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_copyFB2FB

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_accumFB2FB(FBout, FBin, copy, rc)

    ! ----------------------------------------------
    ! Accumulate common field names from FBin to FBout
    ! If copy is passed in and true, the this is a copy
    ! ----------------------------------------------

    type(ESMF_FieldBundle), intent(inout) :: FBout
    type(ESMF_FieldBundle), intent(in)    :: FBin
    logical, optional     , intent(in)    :: copy
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lranki, lranko
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    logical                         :: exists
    logical                         :: lcopy
    real(R8), pointer               :: dataPtri1(:)  , dataPtro1(:)
    real(R8), pointer               :: dataPtri2(:,:), dataPtro2(:,:)
    character(len=*), parameter     :: subname='(lilac_methods_FB_accumFB2FB)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lcopy = .false.  ! accumulate by default
    if (present(copy)) then
      lcopy = copy
    endif

    call ESMF_FieldBundleGet(FBout, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FBout, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
      call ESMF_FieldBundleGet(FBin, fieldName=lfieldnamelist(n), isPresent=exists, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (exists) then
        call lilac_methods_FB_GetFldPtr(FBin,  lfieldnamelist(n), dataPtri1, dataPtri2, lranki, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        call lilac_methods_FB_GetFldPtr(FBout, lfieldnamelist(n), dataPtro1, dataPtro2, lranko, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        if (lranki == 1 .and. lranko == 1) then

          if (.not.lilac_methods_FieldPtr_Compare(dataPtro1, dataPtri1, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR in dataPtr1 size ", ESMF_LOGMSG_ERROR)
            rc = ESMF_FAILURE
            return
          endif

          if (lcopy) then
            do i=lbound(dataPtri1,1),ubound(dataPtri1,1)
              dataPtro1(i) = dataPtri1(i)
            enddo
          else
            do i=lbound(dataPtri1,1),ubound(dataPtri1,1)
              dataPtro1(i) = dataPtro1(i) + dataPtri1(i)
            enddo
          endif

        elseif (lranki == 2 .and. lranko == 2) then

          if (.not.lilac_methods_FieldPtr_Compare(dataPtro2, dataPtri2, subname, rc)) then
            call ESMF_LogWrite(trim(subname)//": ERROR in dataPtr2 size ", ESMF_LOGMSG_ERROR)
            rc = ESMF_FAILURE
            return
          endif

          if (lcopy) then
            do j=lbound(dataPtri2,2),ubound(dataPtri2,2)
            do i=lbound(dataPtri2,1),ubound(dataPtri2,1)
              dataPtro2(i,j) = dataPtri2(i,j)
            enddo
            enddo
          else
            do j=lbound(dataPtri2,2),ubound(dataPtri2,2)
            do i=lbound(dataPtri2,1),ubound(dataPtri2,1)
              dataPtro2(i,j) = dataPtro2(i,j) + dataPtri2(i,j)
            enddo
            enddo
          endif

        else

          write(msgString,'(a,2i8)') trim(subname)//": ranki, ranko = ",lranki,lranko
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
          call ESMF_LogWrite(trim(subname)//": ERROR ranki ranko not supported "//trim(lfieldnamelist(n)), &
               ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return

        endif

      endif
    enddo

    deallocate(lfieldnamelist)

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_accumFB2FB

  !-----------------------------------------------------------------------------

  logical function lilac_methods_FB_FldChk(FB, fldname, rc)

    ! ----------------------------------------------
    ! Determine if field with fldname is in input field bundle
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*)      , intent(in)  :: fldname
    integer               , intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(lilac_methods_FB_FldChk)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    ! If field bundle is not created then set return to .false.
    if (.not. ESMF_FieldBundleIsCreated(FB)) then
       lilac_methods_FB_FldChk = .false.
       return
    end if

    ! If field bundle is created determine if fldname is present in field bundle
    lilac_methods_FB_FldChk = .false.

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) then
       call ESMF_LogWrite(trim(subname)//" Error checking field: "//trim(fldname), & 
            ESMF_LOGMSG_ERROR)
       return
    endif
    if (isPresent) then
       lilac_methods_FB_FldChk = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end function lilac_methods_FB_FldChk

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_Field_GetFldPtr(field, fldptr1, fldptr2, rank, abort, rc)

    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_Field)  , intent(in)              :: field
    real(R8), pointer , intent(inout), optional :: fldptr1(:)
    real(R8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort
    integer           , intent(out)  , optional :: rc

    ! local variables
    type(ESMF_Mesh) :: lmesh
    integer         :: lrank, nnodes, nelements
    logical         :: labort
    character(len=*), parameter :: subname='(lilac_methods_Field_GetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    if (.not.present(rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR rc not present ", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

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

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_Field_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_GetFldPtr(FB, fldname, fldptr1, fldptr2, rank, field, rc)

    ! ----------------------------------------------
    ! Get pointer to a field bundle field
    ! ----------------------------------------------

    type(ESMF_FieldBundle) , intent(in)              :: FB
    character(len=*)       , intent(in)              :: fldname
    real(R8), pointer      , intent(inout), optional :: fldptr1(:)
    real(R8), pointer      , intent(inout), optional :: fldptr2(:,:)
    integer                , intent(out),   optional :: rank
    integer                , intent(out),   optional :: rc
    type(ESMF_Field)       , intent(out),   optional :: field

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    character(len=*), parameter :: subname='(lilac_methods_FB_GetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    if (.not.present(rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR rc not present "//trim(fldname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    rc = ESMF_SUCCESS

    if (.not. lilac_methods_FB_FldChk(FB, trim(fldname), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR field "//trim(fldname)//" not in FB ", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call lilac_methods_Field_GetFldPtr(lfield, &
         fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (present(rank)) then
      rank = lrank
    endif
    if (present(field)) then
       field = lfield
    endif
    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_SetFldPtr(FB, fldname, val, rc)

    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*)      , intent(in)  :: fldname
    real(R8)    , intent(in)  :: val
    integer               , intent(out) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    real(R8), pointer :: fldptr1(:)
    real(R8), pointer :: fldptr2(:,:)
    character(len=*), parameter :: subname='(lilac_methods_FB_SetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call lilac_methods_FB_GetFldPtr(FB, fldname, fldptr1, fldptr2, lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
      ! no local data
    elseif (lrank == 1) then
      fldptr1 = val
    elseif (lrank == 2) then
      fldptr2 = val
    else
       call ESMF_LogWrite(trim(subname)//": ERROR in rank "//trim(fldname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_FB_SetFldPtr

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_State_GetFldPtr(ST, fldname, fldptr1, fldptr2, rank, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    type(ESMF_State),            intent(in)              :: ST
    character(len=*),            intent(in)              :: fldname
    real(R8), pointer, intent(inout), optional :: fldptr1(:)
    real(R8), pointer, intent(inout), optional :: fldptr2(:,:)
    integer         ,            intent(out),   optional :: rank
    integer         ,            intent(out),   optional :: rc

    ! local variables
    type(ESMF_Field)           :: lfield
    integer                    :: lrank
    character(len=*), parameter :: subname='(lilac_methods_State_GetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    if (.not.present(rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR rc not present "//trim(fldname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    rc = ESMF_SUCCESS

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call lilac_methods_Field_GetFldPtr(lfield, &
         fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (present(rank)) then
      rank = lrank
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_State_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_State_SetFldPtr(ST, fldname, val, rc)

    type(ESMF_State)  , intent(in)  :: ST
    character(len=*)  , intent(in)  :: fldname
    real(R8), intent(in)  :: val
    integer           , intent(out) :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    real(R8), pointer :: fldptr1(:)
    real(R8), pointer :: fldptr2(:,:)
    character(len=*), parameter :: subname='(lilac_methods_State_SetFldPtr)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call lilac_methods_State_GetFldPtr(ST, fldname, fldptr1, fldptr2, lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
      ! no local data
    elseif (lrank == 1) then
      fldptr1 = val
    elseif (lrank == 2) then
      fldptr2 = val
    else
       call ESMF_LogWrite(trim(subname)//": ERROR in rank "//trim(fldname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
      rc = ESMF_FAILURE
      return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_State_SetFldPtr

  !-----------------------------------------------------------------------------

  logical function lilac_methods_FieldPtr_Compare1(fldptr1, fldptr2, cstring, rc)

    real(R8), pointer, intent(in)  :: fldptr1(:)
    real(R8), pointer, intent(in)  :: fldptr2(:)
    character(len=*)           , intent(in)  :: cstring
    integer                    , intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(lilac_methods_FieldPtr_Compare1)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lilac_methods_FieldPtr_Compare1 = .false.
    if (lbound(fldptr2,1) /= lbound(fldptr1,1) .or. &
        ubound(fldptr2,1) /= ubound(fldptr1,1)) then
      call ESMF_LogWrite(trim(subname)//": ERROR in data size "//trim(cstring), ESMF_LOGMSG_ERROR, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      write(msgString,*) trim(subname)//': fldptr1 ',lbound(fldptr1),ubound(fldptr1)
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
      write(msgString,*) trim(subname)//': fldptr2 ',lbound(fldptr2),ubound(fldptr2)
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    else
      lilac_methods_FieldPtr_Compare1 = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end function lilac_methods_FieldPtr_Compare1

  !-----------------------------------------------------------------------------

  logical function lilac_methods_FieldPtr_Compare2(fldptr1, fldptr2, cstring, rc)

    real(R8), pointer, intent(in)  :: fldptr1(:,:)
    real(R8), pointer, intent(in)  :: fldptr2(:,:)
    character(len=*)           , intent(in)  :: cstring
    integer                    , intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(lilac_methods_FieldPtr_Compare2)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    lilac_methods_FieldPtr_Compare2 = .false.
    if (lbound(fldptr2,2) /= lbound(fldptr1,2) .or. &
        lbound(fldptr2,1) /= lbound(fldptr1,1) .or. &
        ubound(fldptr2,2) /= ubound(fldptr1,2) .or. &
        ubound(fldptr2,1) /= ubound(fldptr1,1)) then
      call ESMF_LogWrite(trim(subname)//": ERROR in data size "//trim(cstring), ESMF_LOGMSG_ERROR, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      write(msgString,*) trim(subname)//': fldptr2 ',lbound(fldptr2),ubound(fldptr2)
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
      write(msgString,*) trim(subname)//': fldptr1 ',lbound(fldptr1),ubound(fldptr1)
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    else
      lilac_methods_FieldPtr_Compare2 = .true.
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end function lilac_methods_FieldPtr_Compare2

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_Mesh_Print(mesh, string, rc)

    type(ESMF_Mesh) , intent(in)  :: mesh
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    type(ESMF_Distgrid)         :: distgrid
    type(ESMF_DELayout)         :: delayout
    integer                     :: pdim, sdim, nnodes, nelements
    integer                     :: localDeCount
    integer                     :: DeCount
    integer                     :: dimCount, tileCount
    integer, allocatable        :: minIndexPTile(:,:), maxIndexPTile(:,:)
    type(ESMF_MeshStatus_Flag)  :: meshStatus
    logical                     :: elemDGPresent, nodeDGPresent
    character(len=*),parameter  :: subname='(lilac_methods_Mesh_Print)'
    ! ----------------------------------------------

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS

    call ESMF_MeshGet(mesh, elementDistGridIsPresent=elemDGPresent, &
         nodalDistgridIsPresent=nodeDGPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(mesh, status=meshStatus, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! first get the distgrid, which should be available
    if (elemDGPresent) then
       call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": distGrid=element"
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DistGridGet(distgrid, deLayout=deLayout, dimCount=dimCount, &
            tileCount=tileCount, deCount=deCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    dimCount=", dimCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    tileCount=", tileCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    deCount=", deCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DELayoutGet(deLayout, localDeCount=localDeCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    localDeCount=", localDeCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
       allocate(minIndexPTile(dimCount, tileCount), &
            maxIndexPTile(dimCount, tileCount))

       ! get minIndex and maxIndex arrays
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    minIndexPTile=", minIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    maxIndexPTile=", maxIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       deallocate(minIndexPTile, maxIndexPTile)

    endif

    if (nodeDGPresent) then
       call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": distGrid=nodal"
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DistGridGet(distgrid, deLayout=deLayout, dimCount=dimCount, &
            tileCount=tileCount, deCount=deCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    dimCount=", dimCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    tileCount=", tileCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    deCount=", deCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DELayoutGet(deLayout, localDeCount=localDeCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    localDeCount=", localDeCount
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! allocate minIndexPTile and maxIndexPTile accord. to dimCount and tileCount
       allocate(minIndexPTile(dimCount, tileCount), &
            maxIndexPTile(dimCount, tileCount))

       ! get minIndex and maxIndex arrays
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    minIndexPTile=", minIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//":    maxIndexPTile=", maxIndexPTile
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       deallocate(minIndexPTile, maxIndexPTile)

    endif

    if (.not. elemDGPresent .and. .not. nodeDGPresent) then
       call ESMF_LogWrite(trim(subname)//": cannot print distgrid from mesh", &
            ESMF_LOGMSG_WARNING, rc=rc)
       return
    endif

    ! if mesh is complete, also get additional parameters
    if (meshStatus==ESMF_MESHSTATUS_COMPLETE) then
       ! access localDeCount to show this is a real Grid
       call ESMF_MeshGet(mesh, parametricDim=pdim, spatialDim=sdim, &
            numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       write (msgString,*) trim(subname)//":"//trim(string)//": parametricDim=", pdim
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       write (msgString,*) trim(subname)//":"//trim(string)//": spatialDim=", sdim
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       write (msgString,*) trim(subname)//":"//trim(string)//": numOwnedNodes=", nnodes
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       write (msgString,*) trim(subname)//":"//trim(string)//": numOwnedElements=", nelements
       call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_Mesh_Print


!-----------------------------------------------------------------------------
  subroutine lilac_methods_Clock_TimePrint(clock,string,rc)

    ! input/output variables
    type(ESMF_Clock) , intent(in)          :: clock
    character(len=*) , intent(in),optional :: string
    integer          , intent(out)         :: rc

    ! local variables
    type(ESMF_Time)         :: time
    type(ESMF_TimeInterval) :: timeStep
    character(len=CS)       :: timestr
    character(len=CL)       :: lstring
    character(len=*), parameter :: subname='(lilac_methods_Clock_TimePrint)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    if (present(string)) then
      lstring = trim(subname)//":"//trim(string)
    else
      lstring = trim(subname)
    endif

    call ESMF_ClockGet(clock,currtime=time,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": currtime = "//trim(timestr), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock,starttime=time,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(time,timestring=timestr,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": startime = "//trim(timestr), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock,timestep=timestep,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(timestep,timestring=timestr,rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(lstring)//": timestep = "//trim(timestr), ESMF_LOGMSG_INFO)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_Clock_TimePrint

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_Mesh_Write(mesh, string, rc)

    type(ESMF_Mesh) ,intent(in)  :: mesh
    character(len=*),intent(in)  :: string
    integer         ,intent(out) :: rc

    ! local
    integer             :: n,l,i,lsize,ndims
    character(len=CS)   :: name
    type(ESMF_DISTGRID) :: distgrid
    type(ESMF_Array)    :: array
    real(R8), pointer   :: rawdata(:)
    real(R8), pointer   :: coord(:)
    character(len=*),parameter  :: subname='(lilac_methods_Mesh_Write)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

#if (1 == 0)
    !--- elements ---

    call ESMF_MeshGet(mesh, spatialDim=ndims, numownedElements=lsize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(rawdata(ndims*lsize))
    allocate(coord(lsize))

    call ESMF_MeshGet(mesh, elementDistgrid=distgrid, ownedElemCoords=rawdata, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1,ndims
      name = "unknown"
      if (n == 1) name = "lon_element"
      if (n == 2) name = "lat_element"
    do l = 1,lsize
      i = 2*(l-1) + n
      coord(l) = rawdata(i)
      array = ESMF_ArrayCreate(distgrid, farrayPtr=coord, name=name, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      call lilac_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    enddo
    enddo

    deallocate(rawdata,coord)

    !--- nodes ---

    call ESMF_MeshGet(mesh, spatialDim=ndims, numownedNodes=lsize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(rawdata(ndims*lsize))
    allocate(coord(lsize))

    call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, ownedNodeCoords=rawdata, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1,ndims
      name = "unknown"
      if (n == 1) name = "lon_nodes"
      if (n == 2) name = "lat_nodes"
    do l = 1,lsize
      i = 2*(l-1) + n
      coord(l) = rawdata(i)
      array = ESMF_ArrayCreate(distgrid, farrayPtr=coord, name=name, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      call lilac_methods_Array_diagnose(array, string=trim(string)//"_"//trim(name), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(array, trim(string)//"_"//trim(name)//".nc", overwrite=.true., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
    enddo
    enddo

    deallocate(rawdata,coord)
#else
      call ESMF_LogWrite(trim(subname)//": turned off right now", ESMF_LOGMSG_INFO)
#endif

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine lilac_methods_Mesh_Write

  !-----------------------------------------------------------------------------

!================================================================================

  subroutine lilac_methods_State_GetScalar(state, scalar_id, scalar_value, flds_scalar_name, flds_scalar_num, rc)

    ! ----------------------------------------------
    ! Get scalar data from State for a particular name and broadcast it to all other pets
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)     :: state
    integer,          intent(in)     :: scalar_id
    real(R8),         intent(out)    :: scalar_value
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask, ierr, len, icount
    type(ESMF_VM)     :: vm
    type(ESMF_Field)  :: field
    real(R8), pointer :: farrayptr(:,:)
    real(r8)          :: tmp(1)
    character(len=*), parameter :: subname='(lilac_methods_State_GetScalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! check item exist or not?
    call ESMF_StateGet(State, itemSearch=trim(flds_scalar_name), itemCount=icount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (icount > 0) then
      call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      if (mytask == 0) then
        call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
        if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
        endif
        tmp(:) = farrayptr(scalar_id,:)
      endif
      call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      scalar_value = tmp(1)
    else
      call ESMF_LogWrite(trim(subname)//": no ESMF_Field found named: "//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
    end if 

  end subroutine lilac_methods_State_GetScalar

!================================================================================

  subroutine lilac_methods_State_SetScalar(scalar_value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)

    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------

    ! input/output arguments
    real(R8),         intent(in)     :: scalar_value
    integer,          intent(in)     :: scalar_id
    type(ESMF_State), intent(inout)  :: State
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_Field)  :: field
    type(ESMF_VM)     :: vm
    real(R8), pointer :: farrayptr(:,:)
    character(len=*), parameter :: subname='(lilac_methods_State_SetScalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
        call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO)
        rc = ESMF_FAILURE
        return
      endif
      farrayptr(scalar_id,1) = scalar_value
    endif

  end subroutine lilac_methods_State_SetScalar

  !-----------------------------------------------------------------------------

  subroutine lilac_methods_FB_getNumFlds(FB, string, nflds, rc)

    ! ---------------------------------------------- 
    ! Determine if fieldbundle is created and if so, the number of non-scalar
    ! fields in the field bundle
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: string
    integer                , intent(out)   :: nflds
    integer                , intent(inout) :: rc
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. ESMF_FieldBundleIsCreated(FB)) then
       call ESMF_LogWrite(trim(string)//": has not been created, returning", ESMF_LOGMSG_INFO)
       nflds = 0 
    else
       ! Note - the scalar field has been removed from all mediator
       ! field bundles - so this is why we check if the fieldCount is 0 and not 1 here

       call ESMF_FieldBundleGet(FB, fieldCount=nflds, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (nflds == 0) then
          call ESMF_LogWrite(trim(string)//": only has scalar data, returning", ESMF_LOGMSG_INFO)
       end if
    end if

  end subroutine lilac_methods_FB_getNumFlds

!===============================================================================

  logical function ChkErr(rc, line, file, mpierr)

    integer, intent(in) :: rc
    integer, intent(in) :: line

    character(len=*), intent(in) :: file
    logical, optional, intent(in) :: mpierr

    character(MPI_MAX_ERROR_STRING) :: lstring
    integer :: dbrc, lrc, len, ierr

    ChkErr = .false.
    lrc = rc
    if (present(mpierr)) then
       if (mpierr) then
          if (rc == MPI_SUCCESS) return
          call MPI_ERROR_STRING(rc, lstring, len, ierr)
          call ESMF_LogWrite("ERROR: "//trim(lstring), ESMF_LOGMSG_INFO, line=line, file=file, rc=dbrc)
          lrc = ESMF_FAILURE
       endif
    endif

    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       ChkErr = .true.
    endif

  end function ChkErr

end module lilac_methods

