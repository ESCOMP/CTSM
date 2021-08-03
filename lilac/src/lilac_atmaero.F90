module lilac_atmaero

  !-----------------------------------------------------------------------
  ! Contains methods for reading in atmosphere aerosal data
  ! This will be done on the CTSM grid with the CTSM decomposition
  ! (after the redistribution from atm-> lnd)
  !-----------------------------------------------------------------------

  use ESMF

  ! share code uses
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_sys_mod      , only : shr_sys_abort
  use shr_nl_mod       , only : shr_nl_find_group_name
  use shr_log_mod      , only : shr_log_errMsg
  use shr_mpi_mod      , only : shr_mpi_bcast
  use shr_cal_mod      , only : shr_cal_ymd2date

  ! cdeps uses
  use dshr_methods_mod , only : dshr_fldbun_getfldptr
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_init_from_inline, shr_strdata_advance

  ! lilac uses
  use lilac_methods    , only : chkerr
  use lilac_methods    , only : lilac_methods_FB_getFieldN
  use lilac_constants  , only : field_index_unset, logunit
  use ctsm_LilacCouplingFields, only : a2l_fields, lilac_atm2lnd
  use ctsm_LilacCouplingFieldIndices

  implicit none
  private

  type, private :: field_mapping_type
     character(len=:), allocatable :: field_name
     integer :: field_index = field_index_unset
  end type field_mapping_type

  public :: lilac_atmaero_init   ! initialize stream data type sdat
  public :: lilac_atmaero_interp ! interpolates between two years of ndep file data

  ! module data
  type(shr_strdata_type) :: sdat  ! input data stream

  ! The first num_fields_to_read in the fields_to_read list are the fields that this
  ! module will read from data. This is set up to have the same ordering as the fields in
  ! sdat.
  integer :: num_fields_to_read
  type(field_mapping_type), allocatable :: fields_to_read(:)

  character(*),parameter :: u_file_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine lilac_atmaero_init(atm2cpl_state, lilac_clock, rc)

    ! ----------------------------------------
    ! Initialize data stream information.
    ! ----------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(inout) :: atm2cpl_state
    type(ESMF_Clock) , intent(inout) :: lilac_clock
    integer          , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_Mesh)        :: mesh
    type(ESMF_FieldBundle) :: lfieldbundle
    type(ESMF_Field)       :: lfield
    integer                :: mytask                       ! mpi task number
    integer                :: mpicom                       ! mpi communicator
    integer                :: n, nfld                      ! indices
    integer                :: field_index                  ! field index
    integer                :: nunit                        ! namelist input unit
    integer                :: ierr                         ! namelist i/o error flag
    character(len=cl)      :: stream_fldfilename           ! name of input stream datafile
    character(len=cl)      :: stream_meshfile              ! name of input stream meshfile
    character(len=CL)      :: mapalgo = 'bilinear'         ! type of 2d mapping
    character(len=CS)      :: taxmode = 'extend'           ! time extrapolation
    integer                :: stream_year_first            ! first year in stream to use
    integer                :: stream_year_last             ! last year in stream to use
    integer                :: model_year_align             ! align stream_year_first with model year
    type(field_mapping_type), allocatable :: all_fields(:) ! all fields that can possibly be read from data
    character(len=CS), allocatable :: fldlistFile(:)   
    character(len=CS), allocatable :: fldlistModel(:)   
    !-----------------------------------------------------------------------

    namelist /atmaero_stream/ &
         stream_year_first, stream_year_last, model_year_align,  &
         stream_fldfilename, stream_meshfile

    rc = ESMF_SUCCESS

    all_fields = [ &
         field_mapping_type('BCDEPWET', lilac_a2l_Faxa_bcphiwet), &
         field_mapping_type('BCPHODRY', lilac_a2l_Faxa_bcphodry), &
         field_mapping_type('BCPHIDRY', lilac_a2l_Faxa_bcphidry), &
         field_mapping_type('OCDEPWET', lilac_a2l_Faxa_ocphiwet), &
         field_mapping_type('OCPHIDRY', lilac_a2l_Faxa_ocphidry), &
         field_mapping_type('OCPHODRY', lilac_a2l_Faxa_ocphodry), &
         field_mapping_type('DSTX01WD', lilac_a2l_Faxa_dstwet1), &
         field_mapping_type('DSTX01DD', lilac_a2l_Faxa_dstdry1), &
         field_mapping_type('DSTX02WD', lilac_a2l_Faxa_dstwet2), &
         field_mapping_type('DSTX02DD', lilac_a2l_Faxa_dstdry2), &
         field_mapping_type('DSTX03WD', lilac_a2l_Faxa_dstwet3), &
         field_mapping_type('DSTX03DD', lilac_a2l_Faxa_dstdry3), &
         field_mapping_type('DSTX04WD', lilac_a2l_Faxa_dstwet4), &
         field_mapping_type('DSTX04DD', lilac_a2l_Faxa_dstdry4)]

    nfld = 0
    do n = 1, size(all_fields)
       field_index = all_fields(n)%field_index
       if (a2l_fields%is_needed_from_data(field_index)) then
          nfld = nfld + 1
       end if
    end do

    num_fields_to_read = nfld

    if (num_fields_to_read == 0) then
       return
    end if

    allocate(fields_to_read(num_fields_to_read))
    allocate(fldlistFile(num_fields_to_read))
    allocate(fldlistModel(num_fields_to_read))
    nfld = 0
    do n = 1, size(all_fields)
       field_index = all_fields(n)%field_index
       if (a2l_fields%is_needed_from_data(field_index)) then
          nfld = nfld + 1
          fields_to_read(nfld) = all_fields(n)
          fldListFile(nfld) = trim(fields_to_read(nfld)%field_name)
          fldListModel(nfld) = trim(a2l_fields%get_fieldname(field_index))
       end if
    end do
          
    ! default values for namelist
    stream_year_first  = 1                ! first year in stream to use
    stream_year_last   = 1                ! last  year in stream to use
    model_year_align   = 1                ! align stream_year_first with this model year
    stream_fldFileName = ' '
    stream_meshfile    = ' '

    ! get mytask and mpicom
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, localPet=mytask, mpiCommunicator=mpicom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Read namelist
    if (mytask == 0) then
       open(newunit=nunit, file='lilac_in', status='old', iostat=ierr )
       call shr_nl_find_group_name(nunit, 'atmaero_stream', status=ierr)
       if (ierr == 0) then
          read(nunit, atmaero_stream, iostat=ierr)
          if (ierr /= 0) then
             call shr_sys_abort(' ERROR reading namelist '//shr_log_errMsg(u_file_u, __LINE__))
          end if
       else
          call shr_sys_abort(' ERROR finding namelist '//shr_log_errMsg(u_file_u, __LINE__))
       end if
       close(nunit)
    endif
    call shr_mpi_bcast(stream_year_first , mpicom)
    call shr_mpi_bcast(stream_year_last  , mpicom)
    call shr_mpi_bcast(model_year_align  , mpicom)
    call shr_mpi_bcast(stream_fldfilename, mpicom)
    call shr_mpi_bcast(stream_meshfile   , mpicom)

    if (mytask == 0) then
       print *, ' '
       write(logunit,'(a)') 'atmaero stream settings:'
       write(logunit,'(a,i8)')'  stream_year_first  = ',stream_year_first
       write(logunit,'(a,i8)')'  stream_year_last   = ',stream_year_last
       write(logunit,'(a,i8)')'  model_year_align   = ',model_year_align
       write(logunit,'(a)'   )'  stream_fldFileName = ',trim(stream_fldFileName)
       write(logunit,'(a)'   )'  stream_meshfile    = ',trim(stream_meshfile)
       print *, ' '
    endif

    ! ------------------------------
    ! obtain atm mesh
    ! ------------------------------

    call ESMF_StateGet(atm2cpl_state, 'a2c_fb', lfieldbundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call lilac_methods_FB_getFieldN(lfieldbundle, fieldnum=1, field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, mesh=mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------
    ! create the stream data sdat
    ! ------------------------------
    call shr_strdata_init_from_inline(sdat,                  &
         my_task             = mytask,                       &
         logunit             = logunit,                      &
         compname            = 'LND',                        &
         model_clock         = lilac_clock,                  &
         model_mesh          = mesh,                         &
         stream_meshfile     = trim(stream_meshfile),        &
         stream_lev_dimname  = 'null',                       & 
         stream_mapalgo      = trim(mapalgo),                &
         stream_filenames    = (/trim(stream_fldfilename)/), &
         stream_fldlistFile  = fldlistFile,                  &
         stream_fldListModel = fldlistModel,                 &
         stream_yearFirst    = stream_year_first,            &
         stream_yearLast     = stream_year_last,             &
         stream_yearAlign    = model_year_align,             &
         stream_offset       = 0,                            &
         stream_taxmode      = taxmode,                      &
         stream_dtlimit      = 1.5_r8,                       &
         stream_tintalgo     = 'linear',                     &
         stream_name         = 'ATMAERO data ',              &  
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine lilac_atmaero_init

  !================================================================

  subroutine lilac_atmaero_interp(clock, rc)

    ! input/output variables
    type(ESMF_Clock)       :: clock
    integer, intent(out)   :: rc

    ! local variables
    type(ESMF_Time)   :: currTime
    integer           :: yy, mm, dd, sec, curr_ymd
    integer           :: n
    integer           :: field_index
    character(len=CS) :: stream_varname
    real(r8), pointer :: dataptr1d(:)
    character(len=*), parameter :: subname='lilac_atmaero: [lilac_atmaero_interp]'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (num_fields_to_read == 0) then
       return
    end if

    ! get current time info
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    ! advance the streams
    call shr_strdata_advance(sdat, ymd=curr_ymd, tod=sec, logunit=logunit, istr='atmaero', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! obtain the stream data
    do n = 1, num_fields_to_read
       field_index = fields_to_read(n)%field_index
       stream_varname = a2l_fields%get_fieldname(field_index)
       write(6,*)'DEBUG: stream_varname = ',trim(stream_varname)
       call dshr_fldbun_getFldPtr(sdat%pstrm(1)%fldbun_model, trim(stream_varname), fldptr1=dataptr1d, rc=rc)

       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       call lilac_atm2lnd(field_index, dataptr1d)
    end do

  end subroutine lilac_atmaero_interp

end module lilac_atmaero
