module laiStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read LAI from stream
  !
  ! !USES:
  use ESMF             , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_Finalize, ESMF_END_ABORT
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_CL, CS => shr_kind_CS
  use dshr_strdata_mod , only : shr_strdata_type
  use decompMod        , only : bounds_type
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog, FL => fname_len
  use perf_mod         , only : t_startf, t_stopf
  use spmdMod          , only : masterproc, mpicom, iam
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lai_init    ! position datasets for LAI
  public :: lai_advance ! Advance the LAI streams (outside of a Open-MP threading loop)
  public :: lai_interp  ! interpolates between two years of LAI data (when LAI streams

  ! !PRIVATE MEMBER DATA:
  integer, allocatable        :: g_to_ig(:)         ! Array matching gridcell index to data index
  type(shr_strdata_type)      :: sdat_lai           ! LAI input data stream
  character(*), parameter     :: laiString = "LAI_" ! base string for field string
  integer     , parameter     :: numLaiFields = 16  ! number of fields to build field string
  character(len=CS)           :: stream_varnames(numLaiFields)

  character(len=*), parameter :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine lai_init(bounds)
    !
    ! Initialize data stream information for LAI.
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
    integer                 :: i,n, ig, g                 ! index
    integer                 :: stream_year_first_lai      ! first year in Lai stream to use
    integer                 :: stream_year_last_lai       ! last year in Lai stream to use
    integer                 :: model_year_align_lai       ! align stream_year_first_lai with
    integer                 :: nu_nml                     ! unit for namelist file
    integer                 :: nml_error                  ! namelist i/o error flag
    character(len=FL)       :: stream_fldFileName_lai     ! lai stream filename to read
    character(len=FL)       :: stream_meshfile_lai        ! lai stream meshfile
    real(r8)                :: lai_dtlimit = 1.5_r8       ! dlimit for lai stream to use
    character(len=CL)       :: lai_mapalgo = 'bilinear'   ! Mapping alogrithm
    character(len=CL)       :: lai_tintalgo = 'linear'    ! Time interpolation alogrithm
    integer                 :: lai_offset = 0             ! Offset in time for dataset (sec)
    integer                 :: rc
    character(*), parameter :: subName = "('laidyn_init')"
    !-----------------------------------------------------------------------
    !
    ! deal with namelist variables here in init
    !
    namelist /lai_streams/         &
         stream_year_first_lai,    &
         stream_year_last_lai,     &
         model_year_align_lai,     &
         lai_dtlimit,              &
         lai_mapalgo,              &
         stream_fldFileName_lai,   &
         stream_meshfile_lai,      &
         lai_tintalgo

    ! Default values for namelist
    stream_year_first_lai  = 1      ! first year in stream to use
    stream_year_last_lai   = 1      ! last  year in stream to use
    model_year_align_lai   = 1      ! align stream_year_first_lai with this model year
    stream_fldFileName_lai = ''
    stream_meshfile_lai    = ''
    do n = 1,numLaiFields
       write(stream_varnames(n),'(a,i0)') laiString,n
    end do

    ! Read lai_streams namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'lai_streams', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=lai_streams,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading lai_streams namelist')
          end if
       else
          call endrun(subname // ':: ERROR finding lai_streams namelist')
       end if
       close(nu_nml)
    endif
    call shr_mpi_bcast(stream_year_first_lai  , mpicom)
    call shr_mpi_bcast(stream_year_last_lai   , mpicom)
    call shr_mpi_bcast(model_year_align_lai   , mpicom)
    call shr_mpi_bcast(stream_fldFileName_lai , mpicom)
    call shr_mpi_bcast(stream_meshfile_lai    , mpicom)
    call shr_mpi_bcast(lai_tintalgo           , mpicom)
    call shr_mpi_bcast(lai_dtlimit            , mpicom)
    call shr_mpi_bcast(lai_mapalgo            , mpicom)

    if (masterproc) then
       write(iulog,*)
       write(iulog,'(a)') 'lai_stream settings:'
       write(iulog,'(a,i8)') '  stream_year_first_lai  = ',stream_year_first_lai
       write(iulog,'(a,i8)') '  stream_year_last_lai   = ',stream_year_last_lai
       write(iulog,'(a,i8)') '  model_year_align_lai   = ',model_year_align_lai
       write(iulog,'(a,a)' ) '  stream_fldFileName_lai = ',trim(stream_fldFileName_lai)
       write(iulog,'(a,a)' ) '  stream_meshfile_lai    = ',trim(stream_meshfile_lai)
       write(iulog,'(a,a)' ) '  lai_tintalgo           = ',trim(lai_tintalgo)
       write(iulog,'(a,a)' ) '  lai_mapalgo            = ',trim(lai_mapalgo)
       write(iulog,'(a,a)' ) '  lai_dtlimit            = ', lai_dtlimit

       do n = 1,numLaiFields
          write(iulog,'(a,a)' ) '  stream_varname         = ',trim(stream_varnames(n))
       end do
       write(iulog,*)
    endif

    ! Initialize the cdeps data type sdat_lai
    call shr_strdata_init_from_inline(sdat_lai,                  &
         my_task             = iam,                              &
         logunit             = iulog,                            &
         compname            = 'LND',                            &
         model_clock         = model_clock,                      &
         model_mesh          = mesh,                             &
         stream_meshfile     = trim(stream_meshfile_lai),        &
         stream_lev_dimname  = 'null',                           &
         stream_mapalgo      = trim(lai_mapalgo),                &
         stream_filenames    = (/trim(stream_fldfilename_lai)/), &
         stream_fldlistFile  = stream_varnames,                  &
         stream_fldListModel = stream_varnames,                  &
         stream_yearFirst    = stream_year_first_lai,            &
         stream_yearLast     = stream_year_last_lai,             &
         stream_yearAlign    = model_year_align_lai,             &
         stream_offset       = lai_offset,                       &
         stream_taxmode      = 'cycle',                          &
         stream_dtlimit      = lai_dtlimit,                      &
         stream_tintalgo     = lai_tintalgo,                     &
         stream_name         = 'LAI data',                       &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    
    if ( .not. allocated(g_to_ig) )then
       allocate (g_to_ig(bounds%begg:bounds%endg) )
       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          g_to_ig(g) = ig
       end do
    end if

  end subroutine lai_init

  !================================================================
  subroutine lai_advance( bounds )
    !
    ! Advance LAI streams
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
    !-----------------------------------------------------------------------

    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    call shr_strdata_advance(sdat_lai, ymd=mcdate, tod=sec, logunit=iulog, istr='laidyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    if ( .not. allocated(g_to_ig) )then
       allocate (g_to_ig(bounds%begg:bounds%endg) )
       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          g_to_ig(g) = ig
       end do
    end if

  end subroutine lai_advance

  !================================================================
  subroutine lai_interp(bounds, canopystate_inst)
    !
    ! Interpolate data stream information for Lai.
    !
    ! !USES:
    use pftconMod        , only : noveg
    use PatchType        , only : patch
    use CanopyStateType  , only : canopystate_type
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    type(canopystate_type) , intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer           :: ivt, p, ip, ig, n, g
    integer           :: lsize
    integer           :: rc
    real(r8), pointer :: dataptr1d(:)
    real(r8), pointer :: dataptr2d(:,:)
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)

    ! Get pointer for stream data that is time and spatially interpolate to model time and grid
    ! Place all lai data from each type into a temporary 2d array
    lsize = bounds%endg - bounds%begg + 1
    allocate(dataptr2d(lsize, numLaiFields))
    do n = 1,numLaiFields
       call dshr_fldbun_getFldPtr(sdat_lai%pstrm(1)%fldbun_model, trim(stream_varnames(n)), &
            fldptr1=dataptr1d,  rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       ! Note that the size of dataptr1d includes ocean points so it will be around 3x larger than lsize
       ! So an explicit loop is required here
       do g = 1,lsize
          dataptr2d(g,n) = dataptr1d(g)
       end do
    end do

    do p = bounds%begp, bounds%endp
       ivt = patch%itype(p)
       ! Set lai for each gridcell/patch combination
       if (ivt /= noveg) then
          ! vegetated pft
          ig = g_to_ig(patch%gridcell(p))
          canopystate_inst%tlai_input_patch(p) = dataptr2d(ig,ivt)
       else
          ! non-vegetated pft
          canopystate_inst%tlai_input_patch(p) = 0._r8
       endif
    end do
    deallocate(dataptr2d)

  end subroutine lai_interp

end module LaiStreamMod
