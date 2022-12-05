module ndepStreamMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in nitrogen deposition data file
  ! Also includes functions for dynamic ndep file handling and
  ! interpolation.
  !
  ! !USES
  use ESMF             , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_Finalize, ESMF_END_ABORT
  use dshr_strdata_mod , only : shr_strdata_type
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use spmdMod          , only : mpicom, masterproc, iam
  use decompMod        , only : bounds_type
  use clm_varctl       , only : iulog, inst_name
  use abortutils       , only : endrun

  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ndep_init      ! position datasets for dynamic ndep
  public :: ndep_interp    ! interpolates between two years of ndep file data

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: check_units   ! Check the units and make sure they can be used

  ! ! PRIVATE TYPES
  type(shr_strdata_type) :: sdat_ndep                      ! input data stream
  logical                :: divide_by_secs_per_yr = .true. ! divide by the number of seconds per year
  character(len=CS)      :: stream_varnames(1)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine ndep_init(bounds, NLFilename)
    !
    ! Initialize data stream information.
    !
    ! Uses:
    use shr_nl_mod       , only : shr_nl_find_group_name
    use shr_string_mod   , only : shr_string_listGetName, shr_string_listGetNum
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    use shr_mpi_mod      , only : shr_mpi_bcast
    use lnd_comp_shr     , only : mesh, model_clock
    use dshr_strdata_mod , only : shr_strdata_init_from_inline
    !
    ! arguments
    type(bounds_type), intent(in) :: bounds
    character(len=*),  intent(in) :: NLFilename   ! Namelist filename
    !
                                                      ! local variables
    integer                 :: nu_nml                 ! unit for namelist file
    integer                 :: nml_error              ! namelist i/o error flag
    integer                 :: stream_ndep_year_first ! first year in stream to use
    integer                 :: stream_ndep_year_last  ! last year in stream to use
    integer                 :: stream_ndep_year_align ! align stream_year_firstndep with
    real(r8)                :: stream_ndep_dtlimit = 1.0e30_r8
    character(len=CL)       :: stream_ndep_mapalgo = 'bilinear'
    character(len=CL)       :: stream_ndep_tintalgo = 'linear'
    character(len=CS)       :: stream_ndep_taxmode = 'extend'
    character(len=CL)       :: stream_ndep_data_varlist = 'NDEP_year'
    integer                 :: stream_ndep_offset = 0 ! Offset in time for dataset (sec)
    character(len=CL)       :: stream_ndep_data_filename
    character(len=CL)       :: stream_ndep_mesh_filename
    integer                 :: stream_nflds
    integer                 :: rc
    character(*), parameter :: subName = "('ndepdyn_init')"
    !-----------------------------------------------------------------------

    namelist /ndepdyn_nml/          &
         stream_ndep_year_first,    &
         stream_ndep_year_last,     &
         stream_ndep_year_align,    &
         stream_ndep_mapalgo,       &
         stream_ndep_taxmode,       &
         stream_ndep_data_varlist,  &
         stream_ndep_tintalgo,      &
         stream_ndep_data_filename, &
         stream_ndep_mesh_filename

    ! Default values for namelist
    stream_ndep_year_first    = 1 ! first year in stream to use
    stream_ndep_year_last     = 1 ! last  year in stream to use
    stream_ndep_year_align    = 1 ! align stream_ndep_year_first with this model year
    stream_ndep_data_filename = ' '
    stream_ndep_mesh_filename = ' '

    ! Read ndepdyn_nml namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call shr_nl_find_group_name(nu_nml, 'ndepdyn_nml', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=ndepdyn_nml,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg=' ERROR reading ndepdyn_nml namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg=' ERROR finding ndepdyn_nml namelist'//errMsg(sourcefile, __LINE__))
       end if
       close(nu_nml)
    endif

    call shr_mpi_bcast(stream_ndep_year_first    , mpicom)
    call shr_mpi_bcast(stream_ndep_year_last     , mpicom)
    call shr_mpi_bcast(stream_ndep_year_align    , mpicom)
    call shr_mpi_bcast(stream_ndep_data_varlist  , mpicom)
    call shr_mpi_bcast(stream_ndep_taxmode       , mpicom)
    call shr_mpi_bcast(stream_ndep_mapalgo       , mpicom)
    call shr_mpi_bcast(stream_ndep_tintalgo      , mpicom)
    call shr_mpi_bcast(stream_ndep_data_filename , mpicom)
    call shr_mpi_bcast(stream_ndep_mesh_filename , mpicom)

    stream_nflds = shr_string_listGetNum(stream_ndep_data_varlist)      ! Get number of fields in list, fn
    if (stream_nflds /= 1) then
       call endrun(msg=' ERROR stream_nflds is not 1 for '//errMsg(sourcefile, __LINE__))
    end if
    call shr_string_listGetName(stream_ndep_data_varlist, 1, stream_varnames(1))

    if (masterproc) then
       write(iulog,'(a)'   ) ' '
       write(iulog,'(a,i8)') 'ndepdyn stream settings:'
       write(iulog,'(a,i8)') '  stream_ndep_year_first    = ',stream_ndep_year_first
       write(iulog,'(a,i8)') '  stream_ndep_year_last     = ',stream_ndep_year_last
       write(iulog,'(a,i8)') '  stream_ndep_year_align    = ',stream_ndep_year_align
       write(iulog,'(a,a)' ) '  stream_ndep_data_filename = ',trim(stream_ndep_data_filename)
       write(iulog,'(a,a)' ) '  stream_ndep_mesh_filename = ',trim(stream_ndep_mesh_filename)
       write(iulog,'(a,a)' ) '  stream_varnames           = ',trim(stream_varnames(1))
       write(iulog,'(a,a)' ) '  stream_ndep_taxmode       = ',trim(stream_ndep_taxmode)
       write(iulog,'(a,a)' ) '  stream_ndep_tintalgo      = ',trim(stream_ndep_tintalgo)
       write(iulog,'(a)'   ) ' '
    endif

    ! Read in units
    call check_units( stream_ndep_data_filename )

    ! Initialize the cdeps data type sdat_ndep
    call shr_strdata_init_from_inline(sdat_ndep,                    &
         my_task             = iam,                                 &
         logunit             = iulog,                               &
         compname            = 'LND',                               &
         model_clock         = model_clock,                         &
         model_mesh          = mesh,                                &
         stream_meshfile     = trim(stream_ndep_mesh_filename),     &
         stream_lev_dimname  = 'null',                              &
         stream_mapalgo      = trim(stream_ndep_mapalgo),           &
         stream_filenames    = (/trim(stream_ndep_data_filename)/), &
         stream_fldlistFile  = stream_varnames,                     &
         stream_fldListModel = stream_varnames,                     &
         stream_yearFirst    = stream_ndep_year_first,              &
         stream_yearLast     = stream_ndep_year_last,               &
         stream_yearAlign    = stream_ndep_year_align,              &
         stream_offset       = stream_ndep_offset,                  &
         stream_taxmode      = stream_ndep_taxmode,                 &
         stream_dtlimit      = stream_ndep_dtlimit,                 &
         stream_tintalgo     = stream_ndep_tintalgo,                &
         stream_name         = 'Nitrogen deposition data ',         &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine ndep_init

  !================================================================
  subroutine check_units( stream_ndep_data_filename)

    !-------------------------------------------------------------------
    ! Check that units are correct on the file and if need any conversion

    use ncdio_pio     , only : ncd_pio_openfile, ncd_inqvid, ncd_getatt, ncd_pio_closefile, ncd_nowrite
    use ncdio_pio     , only : file_desc_t, var_desc_t
    use shr_log_mod   , only : errMsg => shr_log_errMsg

    ! Arguments
    character(len=*), intent(in)  :: stream_ndep_data_filename  ! ndep filename
    !
    ! Local variables
    type(file_desc_t) :: ncid     ! NetCDF filehandle for ndep file
    type(var_desc_t)  :: vardesc  ! variable descriptor
    integer           :: varid    ! variable index
    logical           :: readvar  ! If variable was read
    character(len=CS) :: ndepunits! ndep units
    !-----------------------------------------------------------------------

    call ncd_pio_openfile( ncid, trim(stream_ndep_data_filename), ncd_nowrite )
    call ncd_inqvid(ncid, stream_varnames(1), varid, vardesc, readvar=readvar)
    if ( readvar ) then
       call ncd_getatt(ncid, varid, "units", ndepunits)
    else
       call endrun(msg=' ERROR finding variable: '//trim(stream_varnames(1))//" in file: "// &
            trim(stream_ndep_data_filename)//errMsg(sourcefile, __LINE__))
    end if
    call ncd_pio_closefile( ncid )

    ! Now check to make sure they are correct
    if (trim(ndepunits) == "g(N)/m2/s"  )then
       divide_by_secs_per_yr = .false.
    else if ( trim(ndepunits) == "g(N)/m2/yr" )then
       divide_by_secs_per_yr = .true.
    else
       call endrun(msg=' ERROR in units for nitrogen deposition equal to: '//trim(ndepunits)//" not units expected"// &
            errMsg(sourcefile, __LINE__))
    end if

  end subroutine check_units

  !================================================================
  subroutine ndep_interp(bounds, atm2lnd_inst)

    !-----------------------------------------------------------------------
    use clm_time_manager , only : get_curr_date, get_curr_days_per_year
    use clm_varcon       , only : secspday
    use atm2lndType      , only : atm2lnd_type
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    !
    ! Arguments
    type(bounds_type) , intent(in)    :: bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_inst
    !
    ! Local variables
    integer :: g, ig
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: dayspyr ! days per year
    integer :: rc
    real(r8), pointer :: dataptr1d(:)
    !-----------------------------------------------------------------------

    ! Advance sdat stream
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    call shr_strdata_advance(sdat_ndep, ymd=mcdate, tod=sec, logunit=iulog, istr='ndepdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Get pointer for stream data that is time and spatially interpolate to model time and grid
    call dshr_fldbun_getFldPtr(sdat_ndep%pstrm(1)%fldbun_model, stream_varnames(1), fldptr1=dataptr1d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Fill in atm2lnd_inst%forc_ndep_grc
    if ( divide_by_secs_per_yr )then
       ig = 0
       dayspyr = get_curr_days_per_year( )
       do g = bounds%begg,bounds%endg
          ig = ig+1
          atm2lnd_inst%forc_ndep_grc(g) = dataptr1d(ig) / (secspday * dayspyr)
       end do
    else
       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          atm2lnd_inst%forc_ndep_grc(g) = dataptr1d(ig)
       end do
    end if

  end subroutine ndep_interp

end module ndepStreamMod
