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
  use shr_strdata_mod  , only : shr_strdata_type, shr_strdata_create
  use shr_strdata_mod  , only : shr_strdata_print, shr_strdata_advance
  use shr_string_mod   , only : shr_string_listAppend
  use shr_cal_mod      , only : shr_cal_ymd2date
  use shr_pio_mod      , only : shr_pio_getiotype
  use mct_mod          , only : mct_avect_indexra, mct_gsmap, mct_ggrid
  use mct_mod          , only : mct_gsmap_init, mct_gsmap_orderedpoints
  use mct_mod          , only : mct_ggrid_init, mct_ggrid_importIAttr, mct_ggrid_importRattr

  ! ctsm uses
  use ncdio_pio        , only : pio_subsystem
  use domainMod        , only : ldomain
  use clm_time_manager , only : get_calendar

  ! lilac uses
  use lilac_atmcap     , only : gindex_atm
  use lilac_methods    , only : chkerr
  use lilac_methods    , only : lilac_methods_FB_getFieldN
  use lilac_constants  , only : field_index_unset
  use ctsm_LilacCouplingFields, only : a2l_fields, lilac_atm2lnd
  use ctsm_LilacCouplingFieldIndices

  implicit none
  private

  type, private :: field_mapping_type
     character(len=:), allocatable :: field_name
     integer :: field_index = field_index_unset
  end type field_mapping_type

  public :: lilac_atmaero_init  ! initialize stream data type sdat
  public :: lilac_atmaero_interp   ! interpolates between two years of ndep file data

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

  subroutine lilac_atmaero_init(atm2cpl_state, rc)

    ! ----------------------------------------
    ! Initialize data stream information.
    ! ----------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(inout) :: atm2cpl_state
    integer          , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_Mesh)        :: lmesh
    type(ESMF_FieldBundle) :: lfieldbundle
    type(ESMF_Field)       :: lfield
    type(mct_ggrid)        :: ggrid_atm                  ! domain information
    type(mct_gsmap)        :: gsmap_atm                  ! decompositoin info
    type(field_mapping_type), allocatable :: all_fields(:) ! all fields that can possibly be read from data
    integer                :: mytask                     ! mpi task number
    integer                :: mpicom                     ! mpi communicator
    integer                :: n                          ! index
    integer                :: field_index
    integer                :: lsize                      ! local size
    integer                :: gsize                      ! global size
    integer                :: nunit                      ! namelist input unit
    integer                :: ierr                       ! namelist i/o error flag
    character(len=cl)      :: stream_fldfilename ! name of input stream file
    character(len=CL)      :: mapalgo = 'bilinear'       ! type of 2d mapping
    character(len=CS)      :: taxmode = 'extend'         ! time extrapolation
    character(len=CL)      :: fldlistFile                ! name of fields in input stream file
    character(len=CL)      :: fldlistModel               ! name of fields in model
    integer                :: stream_year_first  ! first year in stream to use
    integer                :: stream_year_last   ! last year in stream to use
    integer                :: model_year_align   ! align stream_year_first with model year
    integer                :: spatialDim
    integer                :: numOwnedElements
    real(r8), pointer      :: ownedElemCoords(:)
    real(r8), pointer      :: mesh_lons(:)
    real(r8), pointer      :: mesh_lats(:)
    real(r8), pointer      :: mesh_areas(:)
    real(r8), pointer      :: rdata(:)
    integer , pointer      :: idata(:)
    !-----------------------------------------------------------------------

    namelist /atmaero_stream/ &
         stream_year_first, stream_year_last, model_year_align,  &
         stream_fldfilename

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

    num_fields_to_read = 0
    allocate(fields_to_read(size(all_fields)))
    fldlistFile = ' '
    fldlistModel = ' '
    do n = 1, size(all_fields)
       field_index = all_fields(n)%field_index
       if (a2l_fields%is_needed_from_data(field_index)) then
          num_fields_to_read = num_fields_to_read + 1
          fields_to_read(num_fields_to_read) = all_fields(n)
          call shr_string_listAppend(fldlistFile, fields_to_read(num_fields_to_read)%field_name)
          call shr_string_listAppend(fldlistModel, a2l_fields%get_fieldname(field_index))
       end if
    end do

    if (num_fields_to_read == 0) then
       return
    end if

    ! default values for namelist
    stream_year_first  = 1                ! first year in stream to use
    stream_year_last   = 1                ! last  year in stream to use
    model_year_align   = 1                ! align stream_year_first with this model year
    stream_fldFileName = ' '

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

    if (mytask == 0) then
       print *, ' '
       print *, 'atmaero stream settings:'
       print *, '  stream_year_first  = ',stream_year_first
       print *, '  stream_year_last   = ',stream_year_last
       print *, '  model_year_align   = ',model_year_align
       print *, '  stream_fldFileName = ',stream_fldFileName
       print *, ' '
    endif

    ! ------------------------------
    ! create the mct gsmap
    ! ------------------------------
    lsize = size(gindex_atm)
    gsize = ldomain%ni * ldomain%nj
    call mct_gsmap_init( gsmap_atm, gindex_atm, mpicom, 1, lsize, gsize )

    ! ------------------------------
    ! obtain mesh lats, lons and areas
    ! ------------------------------

    call ESMF_StateGet(atm2cpl_state, 'a2c_fb', lfieldbundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call lilac_methods_FB_getFieldN(lfieldbundle, fieldnum=1, field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(lmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (numOwnedElements /= lsize) then
       call shr_sys_abort('ERROR: numOwnedElements is not equal to lsize')
    end if
    allocate(ownedElemCoords(spatialDim*numOwnedElements))

    call ESMF_MeshGet(lmesh, ownedElemCoords=ownedElemCoords, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(mesh_lons(numOwnedElements))
    allocate(mesh_lats(numOwnedElements))
    allocate(mesh_areas(numOwnedElements))
    do n = 1,numOwnedElements
       mesh_lons(n) = ownedElemCoords(2*n-1)
       mesh_lats(n) = ownedElemCoords(2*n)
       mesh_areas(n) = 1.e36 ! hard-wire for now for testing
    end do

    ! ------------------------------
    ! create the mct ggrid
    ! ------------------------------
    call mct_ggrid_init( ggrid=ggrid_atm, CoordChars='lat:lon:hgt', OtherChars='area:aream:mask:frac', lsize=lsize)
    call mct_gsmap_orderedpoints(gsmap_atm, mytask, idata)
    call mct_gGrid_importIAttr(ggrid_atm,'GlobGridNum', idata, lsize)
    call mct_gGrid_importRattr(ggrid_atm,"lon" , mesh_lons , lsize)
    call mct_gGrid_importRattr(ggrid_atm,"lat" , mesh_lats , lsize)
    call mct_gGrid_importRattr(ggrid_atm,"area", mesh_areas, lsize)
    allocate(rdata(lsize))
    rdata(:) = 1._R8
    call mct_gGrid_importRattr(ggrid_atm,"mask", rdata, lsize)
    deallocate(mesh_lons, mesh_lats, mesh_areas, rdata)

    ! ------------------------------
    ! create the stream data sdat
    ! ------------------------------
    call shr_strdata_create(sdat,&
         name          = "atmaero",                            &
         pio_subsystem = pio_subsystem,                        &
         pio_iotype    = shr_pio_getiotype(compid= 1),         &
         mpicom        = mpicom,                               &
         compid        = 1,                                    &
         gsmap         = gsmap_atm,                            &
         ggrid         = ggrid_atm,                            &
         nxg           = ldomain%ni,                           &
         nyg           = ldomain%nj,                           &
         yearFirst     = stream_year_first,            &
         yearLast      = stream_year_last,             &
         yearAlign     = model_year_align,             &
         offset        = 0,                                    &
         domFilePath   = '',                                   &
         domfilename   = trim(stream_fldfilename),     &
         domTvarName   = 'time',                               &
         domXvarName   = 'lon' ,                               &
         domYvarName   = 'lat' ,                               &
         domAreaName   = 'area',                               &
         domMaskName   = 'mask',                               &
         filePath      = '',                                   &
         filename      = (/trim(stream_fldfilename)/), &
         fldListFile   = trim(fldlistFile),                    &
         fldListModel  = trim(fldlistModel),                   &
         fillalgo      = 'none',                               &
         mapalgo       = mapalgo,                              &
         calendar      = get_calendar(),                       &
         taxmode       = taxmode                               )

    if (mytask == 0) then
       call shr_strdata_print(sdat,'ATMAERO data')
    endif

  end subroutine lilac_atmaero_init

  !================================================================

  subroutine lilac_atmaero_interp(clock, rc)

    ! input/output variables
    type(ESMF_Clock)       :: clock
    integer, intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    integer                :: mpicom ! mpi communicator
    integer                :: mytask ! mpi task number
    type(ESMF_FieldBundle) :: lfieldbundle
    type(ESMF_Time)        :: currTime
    integer                :: yy, mm, dd, sec, curr_ymd
    integer                :: n
    character(len=*), parameter :: subname='lilac_atmaero: [lilac_atmaero_interp]'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (num_fields_to_read == 0) then
       return
    end if

    ! get mytask and mpicom
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, localPet=mytask, mpiCommunicator=mpicom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! get current time info
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    ! advance the streams
    call shr_strdata_advance(sdat, curr_ymd, sec, mpicom, 'atmaero')

    do n = 1, num_fields_to_read
       call set_field(n)
    end do

  end subroutine lilac_atmaero_interp

  !==============================================================================

  subroutine set_field(fieldnum)

    ! input/output data
    integer, intent(in) :: fieldnum  ! index into fields_to_read and sdat (which are assumed to have the same ordering)

    ! local data
    integer :: field_index ! index in a2l_fields
    !-----------------------------------------------------------------------

    field_index = fields_to_read(fieldnum)%field_index
    call lilac_atm2lnd(field_index, sdat%avs(1)%rAttr(fieldnum,:))

  end subroutine set_field

end module lilac_atmaero
