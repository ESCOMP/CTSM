module cropcalStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read crop calendars from streams
  !
  ! !USES:
  use ESMF
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_CL, CS => shr_kind_CS
  use dshr_strdata_mod , only : shr_strdata_type
  use decompMod        , only : bounds_type
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog
  use clm_varctl       , only : use_cropcal_rx_sdates, use_cropcal_rx_cultivar_gdds, use_cropcal_streams
  use clm_varpar       , only : mxpft
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

!  ! !PUBLIC MEMBER DATA
!  logical, public, save :: use_cropcal_rx_sdates = .false.
!  logical, public, save :: use_cropcal_rx_cultivar_gdds = .false.
!  logical, public, save :: use_cropcal_streams = .false.

  ! !PRIVATE MEMBER DATA:
  integer, allocatable        :: g_to_ig(:)         ! Array matching gridcell index to data index
  type(shr_strdata_type)      :: sdat_cropcal_sdate           ! sdate input data stream
  type(shr_strdata_type)      :: sdat_cropcal_cultivar_gdds   ! sdate input data stream
  character(len=CS), allocatable :: stream_varnames_sdate(:)
  character(len=CS), allocatable :: stream_varnames_cultivar_gdds(:)
  integer                     :: ncft               ! Number of crop functional types (excl. generic crops)
  character(len=CL)       :: stream_fldFileName_sdate   ! sdate stream filename to read
  character(len=CL)       :: stream_fldFileName_cultivar_gdds ! cultivar growing degree-days stream filename to read

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
    integer                 :: stream_year_first_cropcal  ! first year in crop calendar streams to use
    integer                 :: stream_year_last_cropcal   ! last year in crop calendar streams to use
    integer                 :: model_year_align_cropcal   ! align stream_year_first_cropcal with
    integer                 :: nu_nml                     ! unit for namelist file
    integer                 :: nml_error                  ! namelist i/o error flag
    character(len=CL)       :: stream_meshfile_cropcal    ! crop calendar stream meshfile
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
         stream_year_first_cropcal,    &
         stream_year_last_cropcal,     &
         model_year_align_cropcal,     &
         stream_fldFileName_sdate,     &
         stream_fldFileName_cultivar_gdds, &
         stream_meshfile_cropcal

    ! Default values for namelist
    stream_year_first_cropcal  = 1      ! first year in stream to use
    stream_year_last_cropcal   = 1      ! last  year in stream to use
    model_year_align_cropcal   = 1      ! align stream_year_first_cropcal with this model year
    stream_meshfile_cropcal    = ''
    stream_fldFileName_sdate = ''
    stream_fldFileName_cultivar_gdds = ''
    ! Will need modification to work with mxsowings > 1
    ncft = mxpft - npcropmin + 1 ! Ignores generic crops
    allocate(stream_varnames_sdate(ncft))
    allocate(stream_varnames_cultivar_gdds(ncft))
    do n = 1,ncft
       ivt = npcropmin + n - 1
       write(stream_varnames_sdate(n),'(a,i0)') "sdate1_",ivt
       write(stream_varnames_cultivar_gdds(n),'(a,i0)') "gdd1_",ivt
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
    call shr_mpi_bcast(stream_year_first_cropcal  , mpicom)
    call shr_mpi_bcast(stream_year_last_cropcal   , mpicom)
    call shr_mpi_bcast(model_year_align_cropcal   , mpicom)
    call shr_mpi_bcast(stream_fldFileName_sdate   , mpicom)
    call shr_mpi_bcast(stream_fldFileName_cultivar_gdds, mpicom)
    call shr_mpi_bcast(stream_meshfile_cropcal    , mpicom)

    if (masterproc) then
       write(iulog,*)
       write(iulog,*) 'cropcal_stream settings:'
       write(iulog,'(a,i8)') '  stream_year_first_cropcal  = ',stream_year_first_cropcal
       write(iulog,'(a,i8)') '  stream_year_last_cropcal   = ',stream_year_last_cropcal
       write(iulog,'(a,i8)') '  model_year_align_cropcal   = ',model_year_align_cropcal
       write(iulog,'(a,a)' ) '  stream_fldFileName_sdate   = ',trim(stream_fldFileName_sdate)
       write(iulog,'(a,a)' ) '  stream_fldFileName_cultivar_gdds   = ',trim(stream_fldFileName_cultivar_gdds)
       write(iulog,'(a,a)' ) '  stream_meshfile_cropcal    = ',trim(stream_meshfile_cropcal)
       do n = 1,ncft
          write(iulog,'(a,a)' ) '  stream_varnames_sdate  = ',trim(stream_varnames_sdate(n))
          write(iulog,'(a,a)' ) '  stream_varnames_cultivar_gdds  = ',trim(stream_varnames_cultivar_gdds(n))
       end do
       write(iulog,*)
    endif

    use_cropcal_rx_sdates = stream_fldFileName_sdate /= ''
    use_cropcal_rx_cultivar_gdds = stream_fldFileName_cultivar_gdds /= ''
    if (use_cropcal_rx_cultivar_gdds .and. generate_crop_gdds) then
        write(iulog,*) 'Do not specify both generate_crop_gdds=.true. and stream_fldFileName_cultivar_gdds'
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    use_cropcal_streams = use_cropcal_rx_sdates .or. use_cropcal_rx_cultivar_gdds

    ! Initialize the cdeps data type sdat_cropcal_sdate
    if (use_cropcal_rx_sdates) then
       call shr_strdata_init_from_inline(sdat_cropcal_sdate,          &
            my_task             = iam,                                &
            logunit             = iulog,                              &
            compname            = 'LND',                              &
            model_clock         = model_clock,                        &
            model_mesh          = mesh,                               &
            stream_meshfile     = trim(stream_meshfile_cropcal),      &
            stream_lev_dimname  = 'null',                             &
            stream_mapalgo      = trim(cropcal_mapalgo),              &
            stream_filenames    = (/trim(stream_fldFileName_sdate)/), &
            stream_fldlistFile  = stream_varnames_sdate,              &
            stream_fldListModel = stream_varnames_sdate,              &
            stream_yearFirst    = stream_year_first_cropcal,          &
            stream_yearLast     = stream_year_last_cropcal,           &
            stream_yearAlign    = model_year_align_cropcal,           &
            stream_offset       = cropcal_offset,                     &
            stream_taxmode      = 'cycle',                            &
            stream_dtlimit      = 1.5_r8,                             &
            stream_tintalgo     = cropcal_tintalgo,                   &
            stream_name         = 'sowing date data',                 &
            rc                  = rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

    ! Initialize the cdeps data type sdat_cropcal_cultivar_gdds
    if (use_cropcal_rx_cultivar_gdds) then
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
            stream_yearFirst    = stream_year_first_cropcal,          &
            stream_yearLast     = stream_year_last_cropcal,           &
            stream_yearAlign    = model_year_align_cropcal,           &
            stream_offset       = cropcal_offset,                     &
            stream_taxmode      = 'cycle',                            &
            stream_dtlimit      = 1.5_r8,                             &
            stream_tintalgo     = cropcal_tintalgo,                   &
            stream_name         = 'cultivar gdd data',                &
            rc                  = rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
    end if

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
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: rc
    logical           :: verbose = .false.
    !-----------------------------------------------------------------------

    if (verbose) write(iulog,*) 'cropcal_advance(): Beginning'

    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    if (use_cropcal_rx_sdates) then
       call shr_strdata_advance(sdat_cropcal_sdate, ymd=mcdate, tod=sec, logunit=iulog, istr='cropcaldyn', rc=rc)
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

    if ( .not. allocated(g_to_ig) )then
       allocate (g_to_ig(bounds%begg:bounds%endg) )
       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          g_to_ig(g) = ig
       end do
    end if

    if (verbose) write(iulog,*) 'cropcal_advance(): Ending'

  end subroutine cropcal_advance

  !================================================================

  subroutine cropcal_interp(bounds, num_pcropp, filter_pcropp, crop_inst)
    !
    ! Interpolate data stream information for crop calendars.
    !
    ! !USES:
    use CropType        , only : crop_type
    use PatchType       , only : patch
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    type(crop_type)        , intent(inout) :: crop_inst
    !
    ! !LOCAL VARIABLES:
    integer :: ivt, p, ip, ig
    integer :: nc, fp
    integer           :: n, g
    integer           :: lsize
    integer           :: rc
    real(r8), pointer :: dataptr1d_sdate(:)
    real(r8), pointer :: dataptr2d_sdate(:,:)
    real(r8), pointer :: dataptr1d_cultivar_gdds(:)
    real(r8), pointer :: dataptr2d_cultivar_gdds(:,:)
    logical           :: verbose = .false.
    !-----------------------------------------------------------------------

    if (verbose) write(iulog,*) 'cropcal_interp(): Beginning'

    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)

    ! Get pointer for stream data that is time and spatially interpolate to model time and grid
    ! Place all data from each type into a temporary 2d array
    lsize = bounds%endg - bounds%begg + 1

    ! Read prescribed sowing dates from input files
    allocate(dataptr2d_sdate(lsize, ncft))
    dataptr2d_sdate(:,:) = -5
    if (use_cropcal_rx_sdates) then
       ! Starting with npcropmin will skip generic crops
       if (verbose) write(iulog,*) 'cropcal_interp(): Reading sdate file'
       do n = 1, ncft
          ivt = n + npcropmin - 1
          call dshr_fldbun_getFldPtr(sdat_cropcal_sdate%pstrm(1)%fldbun_model, trim(stream_varnames_sdate(n)), &
               fldptr1=dataptr1d_sdate,  rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
          ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
          ! So an explicit loop is required here
          do g = 1,lsize
   
             ! If read-in value is invalid, allow_unprescribed_planting in CropPhenology()
             if (dataptr1d_sdate(g) <= 0 .or. dataptr1d_sdate(g) > 365) then
                dataptr1d_sdate(g) = -1
             end if
   
            dataptr2d_sdate(g,n) = dataptr1d_sdate(g)
          end do
       end do
   
       ! Set rx_sdate for each gridcell/patch combination
       if (verbose) write(iulog,*) 'cropcal_interp(): Set rx_sdate for each gridcell/patch combination'
       do fp = 1, num_pcropp
          p = filter_pcropp(fp)
          ivt = patch%itype(p)
          ! Will skip generic crops
          if (ivt >= npcropmin) then
             n = ivt - npcropmin + 1
             ! vegetated pft
             ig = g_to_ig(patch%gridcell(p))
             crop_inst%rx_sdates_thisyr(p,1) = dataptr2d_sdate(ig,n)
   
             ! Sanity check: Should only read in valid values
             if (crop_inst%rx_sdates_thisyr(p,1) > 365) then
                 write(iulog,'(a,i0,a,i0)') 'cropcal_interp(): Crop patch (ivt ',ivt,') has dataptr2d prescribed sowing date ',crop_inst%rx_sdates_thisyr(p,1)
                 call ESMF_Finalize(endflag=ESMF_END_ABORT)
             end if
   
             ! Only for first sowing date of the year
             ! The conditional here is to ensure nothing weird happens if it's called incorrectly on day 365
             if (crop_inst%sdates_thisyr(p,1) <= 0) then
                 crop_inst%next_rx_sdate(p) = crop_inst%rx_sdates_thisyr(p,1)
             end if
         else
             write(iulog,'(a,i0)') 'cropcal_interp(), rx_sdates: Crop patch has ivt ',ivt
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
       end do
    end if ! use_cropcal_rx_sdates
    deallocate(dataptr2d_sdate)
   
    allocate(dataptr2d_cultivar_gdds(lsize, ncft))
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
          do g = 1,lsize
   
             !  If read-in value is invalid, have PlantCrop() set gddmaturity to PFT-default value.
             if (dataptr1d_cultivar_gdds(g) < 0 .or. dataptr1d_cultivar_gdds(g) > 1000000._r8) then
                dataptr1d_cultivar_gdds(g) = -1
             end if
            
             dataptr2d_cultivar_gdds(g,n) = dataptr1d_cultivar_gdds(g)
          end do
       end do
   
       ! Set rx_cultivar_gdd for each gridcell/patch combination
       if (verbose) write(iulog,*) 'cropcal_interp(): Set rx_cultivar_gdd for each gridcell/patch combination'
       do fp = 1, num_pcropp
          p = filter_pcropp(fp)

!          if (.not. patch%active(p)) then
!              continue
!          end if

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

             if (ig > lsize) then
                 write(iulog,'(a,i0,a,i0,a)') 'ig (',ig,') > lsize (',lsize,')'
                 call ESMF_Finalize(endflag=ESMF_END_ABORT)
             end if

             crop_inst%rx_cultivar_gdds_thisyr(p,1) = dataptr2d_cultivar_gdds(ig,n)
   
             ! Sanity check: Try to catch uninitialized values
             if (crop_inst%rx_cultivar_gdds_thisyr(p,1) > 1000000._r8) then
                 write(iulog,'(a,i0,a,f20.9)') 'cropcal_interp(): Crop patch (ivt ',ivt,') has rx_cultivar_gdds_thisyr(p,1) HUGE ',crop_inst%rx_cultivar_gdds_thisyr(p,1)
                 call ESMF_Finalize(endflag=ESMF_END_ABORT)
             end if
          else
             write(iulog,'(a,i0)') 'cropcal_interp(), rx_cultivar_gdds: Crop patch has ivt ',ivt
             call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
       end do
      write(iulog,*) 'cropcal_interp(): Reading cultivar_gdds file DONE'
   end if ! use_cropcal_rx_cultivar_gdds

   deallocate(dataptr2d_cultivar_gdds)

   if (verbose) write(iulog,*) 'cropcal_interp(): All done!'


  end subroutine cropcal_interp

end module cropcalStreamMod
