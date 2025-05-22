module cropcalStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read crop calendars from streams
  !
  ! !USES:
  use ESMF             , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_Finalize
  use ESMF             , only : ESMF_END_ABORT
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_CL, CS => shr_kind_CS
  use dshr_strdata_mod , only : shr_strdata_type
  use decompMod        , only : bounds_type
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog
  use clm_varctl       , only : FL => fname_len
  use clm_varctl       , only : use_crop
  use clm_varctl       , only : use_cropcal_rx_swindows, use_cropcal_rx_cultivar_gdds, use_cropcal_streams
  use clm_varctl       , only : adapt_cropcal_rx_cultivar_gdds
  use clm_varpar       , only : mxpft
  use clm_varpar       , only : mxsowings
  use perf_mod         , only : t_startf, t_stopf
  use spmdMod          , only : masterproc, mpicom, iam
  use pftconMod        , only : npcropmin
  use CNPhenologyMod  , only : generate_crop_gdds
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: cropcal_init    ! position datasets for crop calendars
  public :: cropcal_advance ! Advance the crop calendar streams (outside of a Open-MP threading loop)
  public :: cropcal_interp  ! interpolates between two years of crop calendar data

  ! !PRIVATE MEMBER DATA:
  integer, allocatable        :: g_to_ig(:)         ! Array matching gridcell index to data index
  type(shr_strdata_type)      :: sdat_cropcal_swindow_start   ! sowing window start input data stream
  type(shr_strdata_type)      :: sdat_cropcal_swindow_end     ! sowing window end input data stream
  type(shr_strdata_type)      :: sdat_cropcal_cultivar_gdds   ! maturity requirement input data stream
  type(shr_strdata_type)      :: sdat_cropcal_gdd20_baseline  ! GDD20 baseline input data stream
  type(shr_strdata_type)      :: sdat_cropcal_gdd20_season_start   ! gdd20 season start input data stream
  type(shr_strdata_type)      :: sdat_cropcal_gdd20_season_end     ! gdd20 season end input data stream
  character(len=CS), allocatable :: stream_varnames_sdate(:) ! used for both start and end dates
  character(len=CS), allocatable :: stream_varnames_cultivar_gdds(:)
  character(len=CS), allocatable :: stream_varnames_gdd20_baseline(:)
  character(len=CS), allocatable :: stream_varnames_gdd20_season_enddate(:) ! start uses stream_varnames_sdate
  integer                     :: ncft               ! Number of crop functional types (excl. generic crops)
  logical                     :: allow_invalid_swindow_inputs ! Fall back on paramfile sowing windows in cases of invalid values in stream_fldFileName_swindow_start and _end?
  character(len=FL)       :: stream_fldFileName_swindow_start ! sowing window start stream filename to read
  character(len=FL)       :: stream_fldFileName_swindow_end   ! sowing window end stream filename to read
  character(len=FL)       :: stream_fldFileName_cultivar_gdds ! cultivar growing degree-days stream filename to read
  character(len=FL)       :: stream_fldFileName_gdd20_baseline ! GDD20 baseline stream filename to read
  logical                 :: cropcals_rx ! Used only for setting input files in namelist; does nothing in code, but needs to be here so namelist read doesn't crash
  logical                 :: cropcals_rx_adapt ! Used only for setting input files in namelist; does nothing in code, but needs to be here so namelist read doesn't crash
  logical :: stream_gdd20_seasons  ! Read start and end dates for gdd20 seasons from streams instead of using hemisphere-specific values
  logical                     :: allow_invalid_gdd20_season_inputs ! Fall back on hemisphere "warm periods" in cases of invalid values in stream_fldFileName_gdd20_season_start and _end?
  character(len=FL)       :: stream_fldFileName_gdd20_season_start ! Stream filename to read for start of gdd20 season
  character(len=FL)       :: stream_fldFileName_gdd20_season_end ! Stream filename to read for end of gdd20 season

  character(len=*), parameter :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine cropcal_init(bounds)
    !
    ! Initialize data stream information for crop calendars.
    !
    ! !USES:
    use shr_mpi_mod      , only : shr_mpi_bcast
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use lnd_comp_shr     , only : mesh, model_clock
    use dshr_strdata_mod , only : shr_strdata_init_from_inline
    use controlMod       , only : NLFilename
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds          ! bounds
    !
    ! !LOCAL VARIABLES:
    integer                 :: i,n,ivt                    ! index
    integer                 :: stream_year_first_cropcal_swindows  ! first year in sowing window streams to use
    integer                 :: stream_year_last_cropcal_swindows   ! last year in sowing window streams to use
    integer                 :: model_year_align_cropcal_swindows   ! alignment year for sowing window streams
    integer                 :: stream_year_first_cropcal_cultivar_gdds  ! first year in cultivar gdd stream to use
    integer                 :: stream_year_last_cropcal_cultivar_gdds   ! last year in cultivar gdd stream to use
    integer                 :: model_year_align_cropcal_cultivar_gdds   ! alignment year for cultivar gdd stream
    integer                 :: nu_nml                     ! unit for namelist file
    integer                 :: nml_error                  ! namelist i/o error flag
    character(len=FL)       :: stream_meshfile_cropcal    ! crop calendar stream meshfile
    character(len=CL)       :: cropcal_mapalgo  = 'nn'        ! Mapping alogrithm
    character(len=CL)       :: cropcal_tintalgo = 'nearest'   ! Time interpolation alogrithm
    integer                 :: cropcal_offset = 0             ! Offset in time for dataset (sec)
    integer                 :: rc
    character(*), parameter :: subName = "('cropcaldyn_init')"
    !-----------------------------------------------------------------------
    !
    ! deal with namelist variables here in init
    !
    namelist /cropcal_streams/         &
         stream_year_first_cropcal_swindows, &
         stream_year_last_cropcal_swindows,  &
         model_year_align_cropcal_swindows,  &
         stream_year_first_cropcal_cultivar_gdds, &
         stream_year_last_cropcal_cultivar_gdds,  &
         model_year_align_cropcal_cultivar_gdds,  &
         allow_invalid_swindow_inputs, &
         stream_fldFileName_swindow_start, &
         stream_fldFileName_swindow_end,   &
         stream_fldFileName_cultivar_gdds, &
         stream_fldFileName_gdd20_baseline, &
         stream_meshfile_cropcal, &
         cropcals_rx, &
         cropcals_rx_adapt, &
         stream_gdd20_seasons, &
         allow_invalid_gdd20_season_inputs, &
         stream_fldFileName_gdd20_season_start, &
         stream_fldFileName_gdd20_season_end

    ! Default values for namelist
    stream_year_first_cropcal_swindows  = 1      ! first year in sowing window streams to use
    stream_year_last_cropcal_swindows   = 1      ! last  year in sowing window streams to use
    model_year_align_cropcal_swindows   = 1      ! alignment year for sowing window streams
    stream_year_first_cropcal_cultivar_gdds  = 1 ! first year in cultivar gdd stream to use
    stream_year_last_cropcal_cultivar_gdds   = 1 ! last  year in cultivar gdd stream to use
    model_year_align_cropcal_cultivar_gdds   = 1 ! alignment year for cultivar gdd stream
    allow_invalid_swindow_inputs = .false.
    stream_meshfile_cropcal    = ''
    stream_fldFileName_swindow_start = ''
    stream_fldFileName_swindow_end   = ''
    stream_fldFileName_cultivar_gdds = ''
    stream_fldFileName_gdd20_baseline = ''
    stream_gdd20_seasons = .false.
    allow_invalid_gdd20_season_inputs = .false.
    stream_fldFileName_gdd20_season_start = ''
    stream_fldFileName_gdd20_season_end = ''
    ! Will need modification to work with mxsowings > 1
    ncft = mxpft - npcropmin + 1 ! Ignores generic crops
    allocate(stream_varnames_sdate(ncft))
    allocate(stream_varnames_cultivar_gdds(ncft))
    allocate(stream_varnames_gdd20_baseline(ncft))
    allocate(stream_varnames_gdd20_season_enddate(ncft))
    do n = 1,ncft
       ivt = npcropmin + n - 1
       write(stream_varnames_sdate(n),'(a,i0)') "sdate1_",ivt
       write(stream_varnames_cultivar_gdds(n),'(a,i0)') "gdd1_",ivt
       write(stream_varnames_gdd20_baseline(n),'(a,i0)') "gdd20bl_",ivt
       write(stream_varnames_gdd20_season_enddate(n),'(a,i0)') "hdate1_",ivt
    end do

    ! Read cropcal_streams namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'cropcal_streams', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=cropcal_streams,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading cropcal_streams namelist')
          end if
       else
          call endrun(subname // ':: ERROR finding cropcal_streams namelist')
       end if
       close(nu_nml)
    endif
    call shr_mpi_bcast(stream_year_first_cropcal_swindows  , mpicom)
    call shr_mpi_bcast(stream_year_last_cropcal_swindows   , mpicom)
    call shr_mpi_bcast(model_year_align_cropcal_swindows   , mpicom)
    call shr_mpi_bcast(stream_year_first_cropcal_cultivar_gdds, mpicom)
    call shr_mpi_bcast(stream_year_last_cropcal_cultivar_gdds , mpicom)
    call shr_mpi_bcast(model_year_align_cropcal_cultivar_gdds , mpicom)
    call shr_mpi_bcast(allow_invalid_swindow_inputs, mpicom)
    call shr_mpi_bcast(stream_fldFileName_swindow_start, mpicom)
    call shr_mpi_bcast(stream_fldFileName_swindow_end  , mpicom)
    call shr_mpi_bcast(stream_fldFileName_cultivar_gdds, mpicom)
    call shr_mpi_bcast(stream_fldFileName_gdd20_baseline, mpicom)
    call shr_mpi_bcast(stream_meshfile_cropcal    , mpicom)
    call shr_mpi_bcast(stream_gdd20_seasons, mpicom)
    call shr_mpi_bcast(allow_invalid_gdd20_season_inputs, mpicom)
    call shr_mpi_bcast(stream_fldFileName_gdd20_season_start, mpicom)
    call shr_mpi_bcast(stream_fldFileName_gdd20_season_end, mpicom)

    if (masterproc) then
       write(iulog,*)
       write(iulog,*) 'cropcal_stream settings:'
       write(iulog,'(a,i8)') '  stream_year_first_cropcal_swindows  = ',stream_year_first_cropcal_swindows
       write(iulog,'(a,i8)') '  stream_year_last_cropcal_swindows   = ',stream_year_last_cropcal_swindows
       write(iulog,'(a,i8)') '  model_year_align_cropcal_swindows   = ',model_year_align_cropcal_swindows
       write(iulog,'(a,i8)') '  stream_year_first_cropcal_cultivar_gdds  = ',stream_year_first_cropcal_cultivar_gdds
       write(iulog,'(a,i8)') '  stream_year_last_cropcal_cultivar_gdds   = ',stream_year_last_cropcal_cultivar_gdds
       write(iulog,'(a,i8)') '  model_year_align_cropcal_cultivar_gdds   = ',model_year_align_cropcal_cultivar_gdds
       write(iulog,'(a,l1)') '  allow_invalid_swindow_inputs = ',allow_invalid_swindow_inputs
       write(iulog,'(a,a)' ) '  stream_fldFileName_swindow_start   = ',trim(stream_fldFileName_swindow_start)
       write(iulog,'(a,a)' ) '  stream_fldFileName_swindow_end     = ',trim(stream_fldFileName_swindow_end)
       write(iulog,'(a,a)' ) '  stream_fldFileName_cultivar_gdds   = ',trim(stream_fldFileName_cultivar_gdds)
       write(iulog,'(a,a)' ) '  stream_fldFileName_gdd20_baseline  = ',trim(stream_fldFileName_gdd20_baseline)
       write(iulog,'(a,a)' ) '  stream_meshfile_cropcal    = ',trim(stream_meshfile_cropcal)
       write(iulog,'(a,l1)') '  stream_gdd20_seasons  = ',stream_gdd20_seasons
       write(iulog,'(a,l1)') '  allow_invalid_gdd20_season_inputs = ',allow_invalid_gdd20_season_inputs
       write(iulog,'(a,a)' ) '  stream_fldFileName_gdd20_season_start  = ',stream_fldFileName_gdd20_season_start
       write(iulog,'(a,a)' ) '  stream_fldFileName_gdd20_season_end    = ',stream_fldFileName_gdd20_season_end
       do n = 1,ncft
          write(iulog,'(a,a)' ) '  stream_varnames_sdate  = ',trim(stream_varnames_sdate(n))
          write(iulog,'(a,a)' ) '  stream_varnames_cultivar_gdds  = ',trim(stream_varnames_cultivar_gdds(n))
          write(iulog,'(a,a)' ) '  stream_varnames_gdd20_season_enddate  = ',trim(stream_varnames_gdd20_season_enddate(n))
          write(iulog,'(a,a)' ) '  stream_varnames_gdd20_baseline  = ',trim(stream_varnames_gdd20_baseline(n))
       end do
       write(iulog,*)
    endif

    ! CLMBuildNamelist checks that both start and end files are provided if either is
    use_cropcal_rx_swindows      = stream_fldFileName_swindow_start /= ''
    use_cropcal_rx_cultivar_gdds = stream_fldFileName_cultivar_gdds /= ''
    adapt_cropcal_rx_cultivar_gdds = stream_fldFileName_gdd20_baseline /= ''
    use_cropcal_streams = .false.  ! Will be set to true if any file is read

    if (use_cropcal_rx_swindows) then
       use_cropcal_streams = .true.

       ! Initialize the cdeps data type sdat_cropcal_swindow_start
       ! NOTE: stream_dtlimit 1.5 didn't work for some reason
       call shr_strdata_init_from_inline(sdat_cropcal_swindow_start,  &
            my_task             = iam,                                &
            logunit             = iulog,                              &
            compname            = 'LND',                              &
            model_clock         = model_clock,                        &
            model_mesh          = mesh,                               &
            stream_meshfile     = trim(stream_meshfile_cropcal),      &
            stream_lev_dimname  = 'null',                             &
            stream_mapalgo      = trim(cropcal_mapalgo),              &
            stream_filenames    = (/trim(stream_fldFileName_swindow_start)/), &
            stream_fldlistFile  = stream_varnames_sdate,              &
            stream_fldListModel = stream_varnames_sdate,              &
            stream_yearFirst    = stream_year_first_cropcal_swindows, &
            stream_yearLast     = stream_year_last_cropcal_swindows,  &
            stream_yearAlign    = model_year_align_cropcal_swindows,           &
            stream_offset       = cropcal_offset,                     &
            stream_taxmode      = 'extend',                           &
            stream_dtlimit      = 1.0e30_r8,                          &
            stream_tintalgo     = cropcal_tintalgo,                   &
            stream_name         = 'sowing window start data',         &
            rc                  = rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if

       ! Initialize the cdeps data type sdat_cropcal_swindow_end
       ! NOTE: stream_dtlimit 1.5 didn't work for some reason
       call shr_strdata_init_from_inline(sdat_cropcal_swindow_end,    &
            my_task             = iam,                                &
            logunit             = iulog,                              &
            compname            = 'LND',                              &
            model_clock         = model_clock,                        &
            model_mesh          = mesh,                               &
            stream_meshfile     = trim(stream_meshfile_cropcal),      &
            stream_lev_dimname  = 'null',                             &
            stream_mapalgo      = trim(cropcal_mapalgo),              &
            stream_filenames    = (/trim(stream_fldFileName_swindow_end)/), &
            stream_fldlistFile  = stream_varnames_sdate,              &
            stream_fldListModel = stream_varnames_sdate,              &
            stream_yearFirst    = stream_year_first_cropcal_swindows, &
            stream_yearLast     = stream_year_last_cropcal_swindows,  &
            stream_yearAlign    = model_year_align_cropcal_swindows,  &
            stream_offset       = cropcal_offset,                     &
            stream_taxmode      = 'extend',                           &
            stream_dtlimit      = 1.0e30_r8,                          &
            stream_tintalgo     = cropcal_tintalgo,                   &
            stream_name         = 'sowing window end data',           &
            rc                  = rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! Initialize the cdeps data type sdat_cropcal_cultivar_gdds
    ! NOTE: stream_dtlimit 1.5 didn't work for some reason
    if (use_cropcal_rx_cultivar_gdds) then
       use_cropcal_streams = .true.
       call shr_strdata_init_from_inline(sdat_cropcal_cultivar_gdds,  &
            my_task             = iam,                                &
            logunit             = iulog,                              &
            compname            = 'LND',                              &
            model_clock         = model_clock,                        &
            model_mesh          = mesh,                               &
            stream_meshfile     = trim(stream_meshfile_cropcal),      &
            stream_lev_dimname  = 'null',                             &
            stream_mapalgo      = trim(cropcal_mapalgo),              &
            stream_filenames    = (/trim(stream_fldFileName_cultivar_gdds)/), &
            stream_fldlistFile  = stream_varnames_cultivar_gdds,      &
            stream_fldListModel = stream_varnames_cultivar_gdds,      &
            stream_yearFirst    = stream_year_first_cropcal_cultivar_gdds,&
            stream_yearLast     = stream_year_last_cropcal_cultivar_gdds, &
            stream_yearAlign    = model_year_align_cropcal_cultivar_gdds, &
            stream_offset       = cropcal_offset,                     &
            stream_taxmode      = 'extend',                           &
            stream_dtlimit      = 1.0e30_r8,                          &
            stream_tintalgo     = cropcal_tintalgo,                   &
            stream_name         = 'cultivar gdd data',                &
            rc                  = rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! Initialize the cdeps data type sdat_cropcal_gdd20_baseline
    ! NOTE: Hard-coded to one particular year because it should NOT vary over time. Note that the
    ! particular year chosen doesn't matter. Users can base their file on whatever baseline they
    ! want; they just need to put 2000 on the time axis.
    if (adapt_cropcal_rx_cultivar_gdds) then
       use_cropcal_streams = .true.
       call shr_strdata_init_from_inline(sdat_cropcal_gdd20_baseline,  &
            my_task             = iam,                                &
            logunit             = iulog,                              &
            compname            = 'LND',                              &
            model_clock         = model_clock,                        &
            model_mesh          = mesh,                               &
            stream_meshfile     = trim(stream_meshfile_cropcal),      &
            stream_lev_dimname  = 'null',                             &
            stream_mapalgo      = 'nn',                         &
            stream_filenames    = (/trim(stream_fldFileName_gdd20_baseline)/), &
            stream_fldlistFile  = stream_varnames_gdd20_baseline,     &
            stream_fldListModel = stream_varnames_gdd20_baseline,     &
            stream_yearFirst    = 2000,                               &
            stream_yearLast     = 2000,                               &
            stream_yearAlign    = 2000,                               &
            stream_offset       = cropcal_offset,                     &
            stream_taxmode      = 'extend',                           &
            stream_dtlimit      = 1.0e30_r8,                          &
            stream_tintalgo     = cropcal_tintalgo,                   &
            stream_name         = 'GDD20 baseline data',              &
            rc                  = rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    if (stream_gdd20_seasons) then
       use_cropcal_streams = .true.

       ! Initialize the cdeps data type sdat_cropcal_gdd20_season_start
      ! NOTE: Hard-coded to one particular year because it should NOT vary over time. Note that the
      ! particular year chosen doesn't matter.
      call shr_strdata_init_from_inline(sdat_cropcal_gdd20_season_start,  &
           my_task             = iam,                                &
           logunit             = iulog,                              &
           compname            = 'LND',                              &
           model_clock         = model_clock,                        &
           model_mesh          = mesh,                               &
           stream_meshfile     = trim(stream_meshfile_cropcal),      &
           stream_lev_dimname  = 'null',                             &
           stream_mapalgo      = trim(cropcal_mapalgo),              &
           stream_filenames    = (/trim(stream_fldFileName_gdd20_season_start)/), &
           stream_fldlistFile  = stream_varnames_sdate,              &
           stream_fldListModel = stream_varnames_sdate,              &
           stream_yearFirst    = 2000,                               &
           stream_yearLast     = 2000,                               &
           stream_yearAlign    = 2000,                               &
           stream_offset       = cropcal_offset,                     &
           stream_taxmode      = 'extend',                           &
           stream_dtlimit      = 1.0e30_r8,                          &
           stream_tintalgo     = cropcal_tintalgo,                   &
           stream_name         = 'gdd20 season start data',         &
           rc                  = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Initialize the cdeps data type sdat_cropcal_gdd20_season_end
      ! NOTE: Hard-coded to one particular year because it should NOT vary over time. Note that the
      ! particular year chosen doesn't matter.
      call shr_strdata_init_from_inline(sdat_cropcal_gdd20_season_end,  &
           my_task             = iam,                                &
           logunit             = iulog,                              &
           compname            = 'LND',                              &
           model_clock         = model_clock,                        &
           model_mesh          = mesh,                               &
           stream_meshfile     = trim(stream_meshfile_cropcal),      &
           stream_lev_dimname  = 'null',                             &
           stream_mapalgo      = trim(cropcal_mapalgo),              &
           stream_filenames    = (/trim(stream_fldFileName_gdd20_season_end)/), &
           stream_fldlistFile  = stream_varnames_gdd20_season_enddate, &
           stream_fldListModel = stream_varnames_gdd20_season_enddate, &
           stream_yearFirst    = 2000,                               &
           stream_yearLast     = 2000,                               &
           stream_yearAlign    = 2000,                               &
           stream_offset       = cropcal_offset,                     &
           stream_taxmode      = 'extend',                           &
           stream_dtlimit      = 1.0e30_r8,                          &
           stream_tintalgo     = cropcal_tintalgo,                   &
           stream_name         = 'gdd20 season start data',         &
           rc                  = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

   end if

   if (masterproc) then
      write(iulog,*)
      write(iulog,*) 'cropcal_stream DERIVED settings:'
      write(iulog,'(a,l1)') '  use_cropcal_rx_swindows  = ',use_cropcal_rx_swindows
      write(iulog,'(a,l1)') '  use_cropcal_rx_cultivar_gdds   = ',use_cropcal_rx_cultivar_gdds
      write(iulog,'(a,l1)') '  adapt_cropcal_rx_cultivar_gdds   = ',adapt_cropcal_rx_cultivar_gdds
      write(iulog,'(a,l1)') '  use_cropcal_streams  = ',use_cropcal_streams
      write(iulog,*)
   endif

  end subroutine cropcal_init

  !================================================================
  subroutine cropcal_advance( bounds )
    !
    ! Advance crop calendar streams
    !
    ! !USES:
    use clm_time_manager , only : get_curr_date
    use dshr_strdata_mod , only : shr_strdata_advance
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g, ig   ! Indices
    integer :: begg, endg  ! gridcell bounds
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: rc
    !-----------------------------------------------------------------------

    begg = bounds%begg
    endg = bounds%endg

    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    if (use_cropcal_rx_swindows) then
       call shr_strdata_advance(sdat_cropcal_swindow_start, ymd=mcdate, tod=sec, logunit=iulog, istr='cropcaldyn', rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       call shr_strdata_advance(sdat_cropcal_swindow_end, ymd=mcdate, tod=sec, logunit=iulog, istr='cropcaldyn', rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if
    if (use_cropcal_rx_cultivar_gdds) then
       call shr_strdata_advance(sdat_cropcal_cultivar_gdds, ymd=mcdate, tod=sec, logunit=iulog, istr='cropcaldyn', rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! The following should not have an associated time axis, but still need to be here
    ! - GDD20 baseline values
    ! - GDD20 season start dates
    ! - GDD20 season end dates
    if (adapt_cropcal_rx_cultivar_gdds) then
       call shr_strdata_advance(sdat_cropcal_gdd20_baseline, ymd=mcdate, tod=sec, logunit=iulog, istr='cropcaldyn', rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if
    if (stream_gdd20_seasons) then
       call shr_strdata_advance(sdat_cropcal_gdd20_season_start, ymd=mcdate, tod=sec, logunit=iulog, istr='cropcaldyn', rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       call shr_strdata_advance(sdat_cropcal_gdd20_season_end, ymd=mcdate, tod=sec, logunit=iulog, istr='cropcaldyn', rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    if ( .not. allocated(g_to_ig) )then
       allocate (g_to_ig(begg:endg) )
       ig = 0
       do g = begg,endg
          ig = ig+1
          g_to_ig(g) = ig
       end do
    end if

  end subroutine cropcal_advance

  !================================================================

  subroutine cropcal_interp(bounds, num_pcropp, filter_pcropp, init, crop_inst)
    !
    ! Interpolate data stream information for crop calendars.
    !
    ! !USES:
    use CropType        , only : crop_type
    use PatchType       , only : patch
    use clm_time_manager, only : get_curr_days_per_year
    use pftconMod       , only : pftname
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                , intent(in)    :: init  ! is this being called as initialization?
    type(crop_type)        , intent(inout) :: crop_inst
    !
    ! !LOCAL VARIABLES:
    integer :: ivt, p, ip, ig
    integer :: nc, fp
    integer :: dayspyr
    integer           :: n, g
    integer           :: rc
    integer           :: begg, endg
    integer           :: begp, endp
    real(r8), pointer :: dataptr1d_swindow_start(:)
    real(r8), pointer :: dataptr1d_swindow_end  (:)
    real(r8), pointer :: dataptr1d_cultivar_gdds(:)
    real(r8), pointer :: dataptr1d_gdd20_baseline(:)
    real(r8), pointer :: dataptr1d_gdd20_season_start(:)
    real(r8), pointer :: dataptr1d_gdd20_season_end  (:)
    real(r8), pointer :: dataptr2d_swindow_start(:,:)
    real(r8), pointer :: dataptr2d_swindow_end  (:,:)
    real(r8), pointer :: dataptr2d_cultivar_gdds(:,:)
    real(r8), pointer :: dataptr2d_gdd20_baseline(:,:)
    real(r8), pointer :: dataptr2d_gdd20_season_start(:,:)
    real(r8), pointer :: dataptr2d_gdd20_season_end  (:,:)
    !-----------------------------------------------------------------------

    associate( &
         swindow_starts => crop_inst%rx_swindow_starts_thisyr_patch, &
         swindow_ends   => crop_inst%rx_swindow_ends_thisyr_patch,   &
         gdd20_season_starts => crop_inst%gdd20_season_start_patch, &
         gdd20_season_ends   => crop_inst%gdd20_season_end_patch    &
         )

    begg = bounds%begg
    endg = bounds%endg
    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= endg ), sourcefile, __LINE__)

    ! Get pointer for stream data that is time and spatially interpolate to model time and grid
    ! Place all data from each type into a temporary 2d array

    begp = bounds%begp
    endp = bounds%endp

    dayspyr = get_curr_days_per_year()

    ! Read prescribed sowing window start dates from input files
    allocate(dataptr2d_swindow_start(begg:endg, ncft))
    dataptr2d_swindow_start(begg:endg,:) = -1._r8
    allocate(dataptr2d_swindow_end  (begg:endg, ncft))
    dataptr2d_swindow_end(begg:endg,:) = -1._r8
    if (use_cropcal_rx_swindows) then
       ! Starting with npcropmin will skip generic crops
       do n = 1, ncft
          call dshr_fldbun_getFldPtr(sdat_cropcal_swindow_start%pstrm(1)%fldbun_model, trim(stream_varnames_sdate(n)), &
               fldptr1=dataptr1d_swindow_start,  rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
          call dshr_fldbun_getFldPtr(sdat_cropcal_swindow_end%pstrm(1)%fldbun_model, trim(stream_varnames_sdate(n)), &
               fldptr1=dataptr1d_swindow_end,  rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
          ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
          ! So an explicit loop is required here
          do g = begg, endg

             ! If read-in value is invalid, set to -1. Will be handled later in this subroutine.
             if (dataptr1d_swindow_start(g) <= 0 .or. dataptr1d_swindow_start(g) > dayspyr &
                 .or. dataptr1d_swindow_end(g) <= 0 .or. dataptr1d_swindow_end(g) > dayspyr) then
                dataptr1d_swindow_start(g) = -1
                dataptr1d_swindow_end  (g) = -1
             end if

            dataptr2d_swindow_start(g,n) = dataptr1d_swindow_start(g)
            dataptr2d_swindow_end  (g,n) = dataptr1d_swindow_end  (g)
          end do
       end do

       ! Set sowing window for each gridcell/patch combination
       do fp = 1, num_pcropp
          p = filter_pcropp(fp)
          ivt = patch%itype(p)
          ! Will skip generic crops
          if (ivt >= npcropmin) then
             n = ivt - npcropmin + 1
             ! vegetated pft
             ig = g_to_ig(patch%gridcell(p))
             swindow_starts(p,1) = dataptr2d_swindow_start(ig,n)
             swindow_ends(p,1)   = dataptr2d_swindow_end  (ig,n)
         else
             write(iulog,'(a,i0)') 'cropcal_interp(), prescribed sowing windows: Crop patch has ivt ',ivt
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
       end do

       ! Ensure that, if mxsowings > 1, sowing windows are ordered such that ENDS are monotonically increasing. This is necessary because of how get_swindow() works.
       if (mxsowings > 1) then
           if (any(swindow_ends(begp:endp,2:mxsowings) <= swindow_ends(begp:endp,1:mxsowings-1) .and. &
                   swindow_ends(begp:endp,2:mxsowings) >= 1)) then
               write(iulog, *) 'Sowing window inputs must be ordered such that end dates are monotonically increasing.'
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
           end if
       end if

       ! Handle invalid sowing window values
       if (any(swindow_starts(begp:endp,:) < 1 .or. swindow_ends(begp:endp,:) < 1)) then
           ! Fail if not allowing fallback to paramfile sowing windows
           if ((.not. allow_invalid_swindow_inputs) .and. any(all(swindow_starts(begp:endp,:) < 1, dim=2) .and. patch%wtgcell(begp:endp) > 0._r8 .and. patch%itype(begp:endp) >= npcropmin)) then
               write(iulog, *) 'At least one crop in one gridcell has invalid prescribed sowing window start date(s). To ignore and fall back to paramfile sowing windows, set allow_invalid_swindow_inputs to .true.'
               write(iulog, *) 'Affected crops:'
               do ivt = npcropmin, mxpft
                   do fp = 1, num_pcropp
                       p = filter_pcropp(fp)
                       if (ivt == patch%itype(p) .and. patch%wtgcell(p) > 0._r8 .and. all(swindow_starts(p,:) < 1)) then
                           write(iulog, *) '    ',pftname(ivt),'  (',ivt,')'
                           exit  ! Stop looking for patches of this type
                       end if
                   end do
               end do
               call ESMF_Finalize(endflag=ESMF_END_ABORT)

           ! Fail if a sowing window start date is prescribed without an end date (or vice versa)
           else if (any((swindow_starts(begp:endp,:) >= 1 .and. swindow_ends(begp:endp,:) < 1) .or. (swindow_starts(begp:endp,:) < 1 .and. swindow_ends(begp:endp,:) >= 1))) then
               write(iulog, *) 'Every prescribed sowing window start date must have a corresponding end date.'
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
           end if
       end if

    end if ! use_cropcal_rx_swindows
    deallocate(dataptr2d_swindow_start)
    deallocate(dataptr2d_swindow_end)
   
    allocate(dataptr2d_cultivar_gdds(begg:endg, ncft))
    if (use_cropcal_rx_cultivar_gdds) then
       ! Read prescribed cultivar GDDs from input files
       ! Starting with npcropmin will skip generic crops
       do n = 1, ncft
          call dshr_fldbun_getFldPtr(sdat_cropcal_cultivar_gdds%pstrm(1)%fldbun_model, trim(stream_varnames_cultivar_gdds(n)), &
               fldptr1=dataptr1d_cultivar_gdds,  rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if

          ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
          ! So an explicit loop is required here
          do g = begg, endg
   
             !  If read-in value is invalid, have PlantCrop() set gddmaturity to PFT-default value.
             if (dataptr1d_cultivar_gdds(g) < 0 .or. dataptr1d_cultivar_gdds(g) > 1000000._r8) then
                dataptr1d_cultivar_gdds(g) = -1
             end if
            
             dataptr2d_cultivar_gdds(g,n) = dataptr1d_cultivar_gdds(g)
          end do
       end do
   
       ! Set rx_cultivar_gdd for each gridcell/patch combination
       do fp = 1, num_pcropp
          p = filter_pcropp(fp)

          ivt = patch%itype(p)
          ! Will skip generic crops
          if (ivt >= npcropmin) then
             n = ivt - npcropmin + 1

             if (n > ncft) then
                 write(iulog,'(a,i0,a,i0,a)') 'n (',n,') > ncft (',ncft,')'
                 call ESMF_Finalize(endflag=ESMF_END_ABORT)
             end if

             ! vegetated pft
             ig = g_to_ig(patch%gridcell(p))

             if (ig < begg .or. ig > endg) then
                 write(iulog,'(a,i0,a,i0,a)') 'ig (',ig,')  < begg (',begg,') or > endg (',endg,')'
                 call ESMF_Finalize(endflag=ESMF_END_ABORT)
             end if

             crop_inst%rx_cultivar_gdds_thisyr_patch(p,1) = dataptr2d_cultivar_gdds(ig,n)
   
          else
             write(iulog,'(a,i0)') 'cropcal_interp(), rx_cultivar_gdds: Crop patch has ivt ',ivt
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
       end do
   end if ! use_cropcal_rx_cultivar_gdds

   deallocate(dataptr2d_cultivar_gdds)

   allocate(dataptr2d_gdd20_baseline(begg:endg, ncft))
   if (adapt_cropcal_rx_cultivar_gdds) then
      ! Read GDD20 baselines from input files
      ! Starting with npcropmin will skip generic crops
      do n = 1, ncft
         call dshr_fldbun_getFldPtr(sdat_cropcal_gdd20_baseline%pstrm(1)%fldbun_model, trim(stream_varnames_gdd20_baseline(n)), &
              fldptr1=dataptr1d_gdd20_baseline,  rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
         ! So an explicit loop is required here
         do g = begg, endg
            dataptr2d_gdd20_baseline(g,n) = dataptr1d_gdd20_baseline(g)
         end do
      end do
  
      ! Set gdd20_baseline_patch for each gridcell/patch combination
      do fp = 1, num_pcropp
         p = filter_pcropp(fp)

         ivt = patch%itype(p)
         ! Will skip generic crops
         if (ivt >= npcropmin) then
            n = ivt - npcropmin + 1

            if (n > ncft) then
                write(iulog,'(a,i0,a,i0,a)') 'n (',n,') > ncft (',ncft,')'
                call ESMF_Finalize(endflag=ESMF_END_ABORT)
            end if

            ! vegetated pft
            ig = g_to_ig(patch%gridcell(p))

            if (ig < begg .or. ig > endg) then
                write(iulog,'(a,i0,a,i0,a)') 'ig (',ig,')  < begg (',begg,') or > endg (',endg,')'
                call ESMF_Finalize(endflag=ESMF_END_ABORT)
            end if

            crop_inst%gdd20_baseline_patch(p) = dataptr2d_gdd20_baseline(ig,n)
  
         else
            write(iulog,'(a,i0)') 'cropcal_interp(), rx_gdd20_baseline: Crop patch has ivt ',ivt
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         endif
      end do
  end if ! adapt_cropcal_rx_cultivar_gdds

  deallocate(dataptr2d_gdd20_baseline)


  ! Read prescribed gdd20 season start dates from input files
  allocate(dataptr2d_gdd20_season_start(begg:endg, ncft))
  dataptr2d_gdd20_season_start(begg:endg,:) = -1._r8
  allocate(dataptr2d_gdd20_season_end  (begg:endg, ncft))
  dataptr2d_gdd20_season_end(begg:endg,:) = -1._r8
  if (stream_gdd20_seasons) then
     ! Starting with npcropmin will skip generic crops
     do n = 1, ncft
        call dshr_fldbun_getFldPtr(sdat_cropcal_gdd20_season_start%pstrm(1)%fldbun_model, trim(stream_varnames_sdate(n)), &
             fldptr1=dataptr1d_gdd20_season_start,  rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
        call dshr_fldbun_getFldPtr(sdat_cropcal_gdd20_season_end%pstrm(1)%fldbun_model, trim(stream_varnames_gdd20_season_enddate(n)), &
             fldptr1=dataptr1d_gdd20_season_end,  rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
        ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
        ! So an explicit loop is required here
        do g = begg, endg

           ! If read-in value is invalid, set to -1. Will be handled later in this subroutine.
           if (dataptr1d_gdd20_season_start(g) <= 0 .or. dataptr1d_gdd20_season_start(g) > 366 &
               .or. dataptr1d_gdd20_season_start(g) /= dataptr1d_gdd20_season_start(g)) then
              dataptr1d_gdd20_season_start(g) = -1
           end if
           if (dataptr1d_gdd20_season_end(g) <= 0 .or. dataptr1d_gdd20_season_end(g) > 366 &
               .or. dataptr1d_gdd20_season_end(g) /= dataptr1d_gdd20_season_end(g)) then
              dataptr1d_gdd20_season_end  (g) = -1
           end if

          dataptr2d_gdd20_season_start(g,n) = dataptr1d_gdd20_season_start(g)
          dataptr2d_gdd20_season_end  (g,n) = dataptr1d_gdd20_season_end  (g)
        end do
     end do

     ! Set gdd20 season for each gridcell/patch combination
     do fp = 1, num_pcropp
        p = filter_pcropp(fp)
        ivt = patch%itype(p)
        ! Will skip generic crops
        if (ivt >= npcropmin) then
           n = ivt - npcropmin + 1
           ! vegetated pft
           ig = g_to_ig(patch%gridcell(p))

           gdd20_season_starts(p) = real(dataptr2d_gdd20_season_start(ig,n), r8)
           gdd20_season_ends(p)   = real(dataptr2d_gdd20_season_end  (ig,n), r8)
       else
           write(iulog,'(a,i0)') 'cropcal_interp(), gdd20 seasons: Crop patch has ivt ',ivt
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
        endif
     end do

     ! Handle invalid gdd20 season values
     if (any(gdd20_season_starts(begp:endp) < 1._r8 .or. gdd20_season_ends(begp:endp) < 1._r8)) then
         ! Fail if not allowing fallback to paramfile sowing windows. Only need to check for
         ! values < 1 because values outside [1, 366] are set to -1 above.
         if ((.not. allow_invalid_gdd20_season_inputs) .and. any(gdd20_season_starts(begp:endp) < 1._r8 .and. patch%wtgcell(begp:endp) > 0._r8 .and. patch%itype(begp:endp) >= npcropmin)) then
             write(iulog, *) 'At least one crop in one gridcell has invalid gdd20 season start and/or end date(s). To ignore and fall back to paramfile sowing windows for such crop-gridcells, set allow_invalid_gdd20_season_inputs to .true.'
             write(iulog, *) 'Affected crops:'
             do ivt = npcropmin, mxpft
                 do fp = 1, num_pcropp
                     p = filter_pcropp(fp)
                     if (ivt == patch%itype(p) .and. patch%wtgcell(p) > 0._r8 .and. gdd20_season_starts(p) < 1._r8) then
                         write(iulog, *) '    ',pftname(ivt),'  (',ivt,')'
                         exit  ! Stop looking for patches of this type
                     end if
                 end do
             end do
             call ESMF_Finalize(endflag=ESMF_END_ABORT)

         ! Fail if a gdd20 season start date is given without an end date (or vice versa)
         else if (any((gdd20_season_starts(begp:endp) >= 1._r8 .and. gdd20_season_ends(begp:endp) < 1._r8) .or. (gdd20_season_starts(begp:endp) < 1._r8 .and. gdd20_season_ends(begp:endp) >= 1._r8))) then
             write(iulog, *) 'Every gdd20 season start date must have a corresponding end date.'
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if
     end if

  end if ! stream_gdd20_seasons
  deallocate(dataptr2d_gdd20_season_start)
  deallocate(dataptr2d_gdd20_season_end)


   end associate

  end subroutine cropcal_interp

end module cropcalStreamMod
