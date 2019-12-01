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
  use mct_mod          , only : mct_avect_indexra, mct_ggrid

  ! ctsm uses
  use ncdio_pio        , only : pio_subsystem
  use decompMod        , only : bounds_type, get_proc_bounds, gsmap_lnd_gdc2glo
  use domainMod        , only : ldomain
  use spmdMod          , only : mpicom, masterproc, comp_id
  use ndepStreamMod    , only : clm_domain_mct
  use clm_time_manager , only : get_calendar

  ! lilac share
  use lilac_methods    , only : chkerr

  implicit none
  private

  public :: lilac_atmaero_init   ! initialize stream data type sdat
  public :: lilac_atmaero_interp ! interpolates between two years of ndep file data

  ! module data
  type(shr_strdata_type) :: sdat  ! input data stream

  character(*),parameter :: u_file_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine lilac_atmaero_init()

    ! ----------------------------------------
    ! Initialize data stream information.
    ! ----------------------------------------

    ! local variables
    integer           :: nunit
    integer           :: ierr                       ! namelist i/o error flag
    type(mct_ggrid)   :: domain_mct                 ! domain information
    character(len=cl) :: stream_fldfilename_atmaero ! name of input stream file
    character(len=CL) :: mapalgo = 'bilinear'       ! type of 2d mapping
    character(len=CS) :: taxmode = 'extend'         ! time extrapolation
    character(len=CL) :: fldlistFile                ! name of fields in input stream file
    character(len=CL) :: fldlistModel               ! name of fields in data stream code
    integer           :: stream_year_first_atmaero  ! first year in stream to use
    integer           :: stream_year_last_atmaero   ! last year in stream to use
    integer           :: model_year_align_atmaero   ! align stream_year_first with model year
    type(bounds_type) :: bounds          
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

    ! Read namelist
    if (masterproc) then
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

    if (masterproc) then
       print *, ' '
       print *, 'atmaero stream settings:'
       print *, '  stream_year_first_atmaero  = ',stream_year_first_atmaero
       print *, '  stream_year_last_atmaero   = ',stream_year_last_atmaero
       print *, '  model_year_align_atmaero   = ',model_year_align_atmaero
       print *, '  stream_fldFileName_atmaero = ',stream_fldFileName_atmaero
       print *, ' '
    endif

    ! Create the mct domain 
    call get_proc_bounds(bounds)
    call clm_domain_mct (bounds, domain_mct)

    ! Create the field list for these urbantv fields...use in shr_strdata_create
    fldlistFile =                      'BCDEPWET:BCPHODRY:BCPHIDRY:'
    fldlistFile = trim(fldlistFile) // 'OCDEPWET:OCPHIDRY:OCPHODRY:DSTX01WD:'
    fldlistFile = trim(fldlistFile) // 'DSTX01DD:DSTX02WD:DSTX02DD:DSTX03WD:'
    fldlistFile = trim(fldlistFile) // 'DSTX03DD:DSTX04WD:DSTX04DD'

    fldlistModel = 'bcphiwet:bcphodry:bcphidry:'
    fldlistModel = trim(fldlistModel) // 'ocphiwet:ocphidry:ocphodry:'
    fldlistModel = trim(fldlistModel) // 'dstwet1:dstdry1:dstwet2:dstdry2'
    fldlistModel = trim(fldlistModel) // 'dstwet3:dstdry3:dstwet4:dstdry4'

    call shr_strdata_create(sdat,&
         name          = "atmaero",                            &
         pio_subsystem = pio_subsystem,                        &
         pio_iotype    = shr_pio_getiotype(compid= 1),         &
         mpicom        = mpicom,                               &
         compid        = comp_id,                              &
         gsmap         = gsmap_lnd_gdc2glo,                    &
         ggrid         = domain_mct,                           &
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

    if (masterproc) then
       call shr_strdata_print(sdat,'ATMAERO data')
    endif

  end subroutine lilac_atmaero_init

  !================================================================

  subroutine lilac_atmaero_interp(c2l_fb, clock, rc)

    ! input/output variables
    type(ESMF_FieldBundle) :: c2l_fb
    type(ESMF_Clock)       :: clock
    integer, intent(out)   :: rc  

    ! local variables
    type(ESMF_Time) :: currTime
    integer :: yy, mm, dd, sec, curr_ymd
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call shr_strdata_advance(sdat, curr_ymd, sec, mpicom, 'atmaero')

    ! Set field bundle data
    call set_fieldbundle_data('Faxa_bcphidry' , c2l_fb, rc)    ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_bcphodry' , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_bcphiwet' , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_ocphidry' , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_ocphodry' , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_ocphiwet' , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstwet1'  , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstdry1'  , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstwet2'  , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstdry2'  , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstwet3'  , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstdry3'  , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstwet4'  , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call set_fieldbundle_data('Faxa_dstdry4'  , c2l_fb, rc=rc) ; if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lilac_atmaero_interp

  !==============================================================================

  subroutine set_fieldbundle_data(fldname, fieldbundle, rc)

    ! input/output data
    character(len=*)       , intent(in)    :: fldname
    type(ESMF_FieldBundle) , intent(inout) :: fieldbundle
    integer                , intent(out)   :: rc

    ! local data
    type(ESMF_field)  :: lfield
    integer           :: n, nfld, indx
    real(r8), pointer :: fldptr1d(:)
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(fieldBundle, fieldName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    nfld = mct_avect_indexra(sdat%avs(1),trim(fldname))
    do indx = 1, size(fldptr1d)
       fldptr1d(n)= sdat%avs(1)%rAttr(nfld,indx)
    end do

  end subroutine set_fieldbundle_data

end module lilac_atmaero
