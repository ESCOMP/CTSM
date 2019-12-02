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
  use lilac_utils      , only : gindex_atm
  use lilac_methods    , only : chkerr
  use lilac_methods    , only : lilac_methods_FB_getFieldN

  implicit none
  private

  public :: lilac_atmaero_init  ! initialize stream data type sdat
  public :: lilac_atmaero_interp   ! interpolates between two years of ndep file data

  ! module data
  type(shr_strdata_type) :: sdat  ! input data stream

  character(*),parameter :: u_file_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine lilac_atmaero_init(atm2lnd_a_state, rc)

    ! ----------------------------------------
    ! Initialize data stream information.
    ! ----------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(inout) :: atm2lnd_a_state
    integer          , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_Mesh)        :: lmesh
    type(ESMF_FieldBundle) :: lfieldbundle
    type(ESMF_Field)       :: lfield
    type(mct_ggrid)        :: ggrid_atm                  ! domain information
    type(mct_gsmap)        :: gsmap_atm                  ! decompositoin info
    integer                :: mytask                     ! mpi task number
    integer                :: mpicom                     ! mpi communicator
    integer                :: n,i,j                      ! index
    integer                :: lsize                      ! local size
    integer                :: gsize                      ! global size
    integer                :: nunit                      ! namelist input unit
    integer                :: ierr                       ! namelist i/o error flag
    character(len=cl)      :: stream_fldfilename_atmaero ! name of input stream file
    character(len=CL)      :: mapalgo = 'bilinear'       ! type of 2d mapping
    character(len=CS)      :: taxmode = 'extend'         ! time extrapolation
    character(len=CL)      :: fldlistFile                ! name of fields in input stream file
    character(len=CL)      :: fldlistModel               ! name of fields in data stream code
    integer                :: stream_year_first_atmaero  ! first year in stream to use
    integer                :: stream_year_last_atmaero   ! last year in stream to use
    integer                :: model_year_align_atmaero   ! align stream_year_first with model year
    integer                :: spatialDim
    integer                :: numOwnedElements
    real(r8), pointer      :: ownedElemCoords(:)
    real(r8), pointer      :: mesh_lons(:)
    real(r8), pointer      :: mesh_lats(:)
    real(r8), pointer      :: mesh_areas(:)
    real(r8), pointer      :: rdata(:)
    integer , pointer      :: idata(:)
    !-----------------------------------------------------------------------

    namelist /atmaero_stream/      &
         stream_year_first_atmaero, &
         stream_year_last_atmaero,  &
         model_year_align_atmaero,  &
         stream_fldfilename_atmaero

    ! default values for namelist
    stream_year_first_atmaero  = 1                ! first year in stream to use
    stream_year_last_atmaero   = 1                ! last  year in stream to use
    model_year_align_atmaero   = 1                ! align stream_year_first_atmaero with this model year
    stream_fldFileName_atmaero = ' '

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

    call shr_mpi_bcast(stream_year_first_atmaero , mpicom)
    call shr_mpi_bcast(stream_year_last_atmaero  , mpicom)
    call shr_mpi_bcast(model_year_align_atmaero  , mpicom)
    call shr_mpi_bcast(stream_fldfilename_atmaero, mpicom)

    if (mytask == 0) then
       print *, ' '
       print *, 'atmaero stream settings:'
       print *, '  stream_year_first_atmaero  = ',stream_year_first_atmaero
       print *, '  stream_year_last_atmaero   = ',stream_year_last_atmaero
       print *, '  model_year_align_atmaero   = ',model_year_align_atmaero
       print *, '  stream_fldFileName_atmaero = ',stream_fldFileName_atmaero
       print *, ' '
    endif

    ! ------------------------------
    ! create the field list for these urbantv fields...use in shr_strdata_create
    ! ------------------------------
    fldlistFile =                      'BCDEPWET:BCPHODRY:BCPHIDRY:'
    fldlistFile = trim(fldlistFile) // 'OCDEPWET:OCPHIDRY:OCPHODRY:DSTX01WD:'
    fldlistFile = trim(fldlistFile) // 'DSTX01DD:DSTX02WD:DSTX02DD:DSTX03WD:'
    fldlistFile = trim(fldlistFile) // 'DSTX03DD:DSTX04WD:DSTX04DD'

    fldlistModel =                       'Faxa_bcphiwet:Faxa_bcphodry:Faxa_bcphidry:'
    fldlistModel = trim(fldlistModel) // 'Faxa_ocphiwet:Faxa_ocphidry:Faxa_ocphodry:'
    fldlistModel = trim(fldlistModel) // 'Faxa_dstwet1:Faxa_dstdry1:Faxa_dstwet2:Faxa_dstdry2:'
    fldlistModel = trim(fldlistModel) // 'Faxa_dstwet3:Faxa_dstdry3:Faxa_dstwet4:Faxa_dstdry4'

    ! ------------------------------
    ! create the mct gsmap
    ! ------------------------------
    lsize = size(gindex_atm)
    gsize = ldomain%ni * ldomain%nj
    call mct_gsmap_init( gsmap_atm, gindex_atm, mpicom, 1, lsize, gsize )

    ! ------------------------------
    ! obtain mesh lats, lons and areas
    ! ------------------------------

    call ESMF_StateGet(atm2lnd_a_state, 'a2c_fb', lfieldbundle, rc=rc)
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
         yearFirst     = stream_year_first_atmaero,            &
         yearLast      = stream_year_last_atmaero,             &
         yearAlign     = model_year_align_atmaero,             &
         offset        = 0,                                    &
         domFilePath   = '',                                   &
         domfilename   = trim(stream_fldfilename_atmaero),     &
         domTvarName   = 'time',                               &
         domXvarName   = 'lon' ,                               &
         domYvarName   = 'lat' ,                               &
         domAreaName   = 'area',                               &
         domMaskName   = 'mask',                               &
         filePath      = '',                                   &
         filename      = (/trim(stream_fldfilename_atmaero)/), &
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

  subroutine lilac_atmaero_interp(atm2lnd_a_state, clock, rc)

    ! input/output variables
    type(ESMF_State)       :: atm2lnd_a_state
    type(ESMF_Clock)       :: clock
    integer, intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    integer                :: mpicom ! mpi communicator
    integer                :: mytask ! mpi task number
    type(ESMF_FieldBundle) :: lfieldbundle
    type(ESMF_Time)        :: currTime
    integer                :: yy, mm, dd, sec, curr_ymd
    character(len=*), parameter :: subname='lilac_atmaero: [lilac_atmaero_interp]'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

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

    ! set field bundle data
    call ESMF_StateGet(atm2lnd_a_state, "a2c_fb", lfieldbundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call set_fieldbundle_data('Faxa_bcphidry' , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_bcphodry' , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_bcphiwet' , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_ocphidry' , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_ocphodry' , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_ocphiwet' , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstwet1'  , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstdry1'  , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstwet2'  , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstdry2'  , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstwet3'  , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstdry3'  , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstwet4'  , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstdry4'  , lfieldbundle, rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lilac_atmaero_interp

  !==============================================================================

  subroutine set_fieldbundle_data(fldname, fieldbundle, rc)

    ! input/output data
    character(len=*)       , intent(in)    :: fldname
    type(ESMF_FieldBundle) , intent(inout) :: fieldbundle
    integer                , intent(out)   :: rc

    ! local data
    type(ESMF_field)  :: lfield
    integer           :: nfld, i
    real(r8), pointer :: fldptr1d(:)
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fieldBundle, fieldName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! error check
    if (size(fldptr1d) /= size(sdat%avs(1)%rAttr, dim=2)) then
       call shr_sys_abort("ERROR: size of fldptr1d and sdat%avs(1)%rattr dim2 are not equal")
    end if

    nfld = mct_avect_indexra(sdat%avs(1),trim(fldname))
    do i = 1, size(fldptr1d)
       fldptr1d(i)= sdat%avs(1)%rAttr(nfld,i)
    end do

  end subroutine set_fieldbundle_data

end module lilac_atmaero
