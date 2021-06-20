module laiStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read LAI from stream
  !
  ! !USES:
  use shr_strdata_mod , only : shr_strdata_type, shr_strdata_create
  use shr_strdata_mod , only : shr_strdata_print, shr_strdata_advance
  use shr_kind_mod    , only : r8=>shr_kind_r8, CL=>shr_kind_CL, CS=>shr_kind_CS, CXX=>shr_kind_CXX
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : iulog, inst_name
  use perf_mod        , only : t_startf, t_stopf
  use spmdMod         , only : masterproc, mpicom, comp_id
  use ncdio_pio
  use mct_mod
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lai_init    ! position datasets for LAI
  public :: lai_advance ! Advance the LAI streams (outside of a Open-MP threading loop)
  public :: lai_interp  ! interpolates between two years of LAI data (when LAI streams

  ! !PRIVATE MEMBER DATA:
  integer, allocatable    :: g_to_ig(:)         ! Array matching gridcell index to data index
  type(shr_strdata_type)  :: sdat_lai           ! LAI input data stream

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
    use clm_time_manager , only : get_calendar
    use ncdio_pio        , only : pio_subsystem
    use shr_pio_mod      , only : shr_pio_getiotype
    use shr_stream_mod   , only : shr_stream_file_null
    use shr_string_mod   , only : shr_string_listCreateField
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use ndepStreamMod    , only : clm_domain_mct
    use histFileMod      , only : hist_addfld1d
    use domainMod        , only : ldomain
    use decompMod        , only : gsmap_global
    use controlMod       , only : NLFilename
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds          ! bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: stream_year_first_lai      ! first year in Lai stream to use
    integer            :: stream_year_last_lai       ! last year in Lai stream to use
    integer            :: model_year_align_lai       ! align stream_year_first_lai with
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    type(mct_ggrid)    :: dom_clm                    ! domain information
    character(len=CL)  :: stream_fldFileName_lai     ! lai stream filename to read
    character(len=CL)  :: lai_mapalgo = 'bilinear'   ! Mapping alogrithm
    character(len=CL)  :: lai_tintalgo = 'linear'    ! Time interpolation alogrithm
    character(len=CXX) :: fldList                    ! field string
    character(*), parameter :: laiString = "LAI"     ! base string for field string
    integer     , parameter :: numLaiFields = 16     ! number of fields to build field string
    character(*), parameter :: subName = "('laidyn_init')"
    !-----------------------------------------------------------------------
    !
    ! deal with namelist variables here in init
    !
    namelist /lai_streams/         &
         stream_year_first_lai,    &
         stream_year_last_lai,     &
         model_year_align_lai,     &
         lai_mapalgo,              &
         stream_fldFileName_lai,   &
         lai_tintalgo

    ! Default values for namelist
    stream_year_first_lai     = 1      ! first year in stream to use
    stream_year_last_lai      = 1      ! last  year in stream to use
    model_year_align_lai      = 1      ! align stream_year_first_lai with this model year
    stream_fldFileName_lai    = shr_stream_file_null

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
    call shr_mpi_bcast(lai_tintalgo           , mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'lai_stream settings:'
       write(iulog,*) '  stream_year_first_lai  = ',stream_year_first_lai
       write(iulog,*) '  stream_year_last_lai   = ',stream_year_last_lai
       write(iulog,*) '  model_year_align_lai   = ',model_year_align_lai
       write(iulog,*) '  stream_fldFileName_lai = ',trim(stream_fldFileName_lai)
       write(iulog,*) '  lai_tintalgo           = ',trim(lai_tintalgo)
    endif

    call clm_domain_mct (bounds, dom_clm)

    ! create the field list for these lai fields...use in shr_strdata_create
    fldList = shr_string_listCreateField( numLaiFields, laiString )

    call shr_strdata_create(sdat_lai,name="laidyn",    &
         pio_subsystem=pio_subsystem,                  &
         pio_iotype=shr_pio_getiotype(inst_name),      &
         mpicom=mpicom, compid=comp_id,                &
         gsmap=gsmap_global, ggrid=dom_clm,            &
         nxg=ldomain%ni, nyg=ldomain%nj,               &
         yearFirst=stream_year_first_lai,              &
         yearLast=stream_year_last_lai,                &
         yearAlign=model_year_align_lai,               &
         offset=0,                                     &
         domFilePath='',                               &
         domFileName=trim(stream_fldFileName_lai),     &
         domTvarName='time',                           &
         domXvarName='lon' ,                           &
         domYvarName='lat' ,                           &
         domAreaName='area',                           &
         domMaskName='mask',                           &
         filePath='',                                  &
         filename=(/stream_fldFileName_lai/),          &
         fldListFile=fldList,                          &
         fldListModel=fldList,                         &
         fillalgo='none',                              &
         mapalgo=lai_mapalgo,                          &
         tintalgo=lai_tintalgo,                        &
         calendar=get_calendar(),                      &
         taxmode='cycle'                               )

    if (masterproc) then
       call shr_strdata_print(sdat_lai,'LAI data')
    endif

  end subroutine lai_init

  !==============================================================================
  subroutine lai_advance( bounds )
    !
    ! Advance LAI streams
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g, ig   ! Indices
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    !-----------------------------------------------------------------------

    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(sdat_lai, mcdate, sec, mpicom, 'laidyn')
    if ( .not. allocated(g_to_ig) )then
       allocate (g_to_ig(bounds%begg:bounds%endg) )
       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          g_to_ig(g) = ig
       end do
    end if

  end subroutine lai_advance

  !==============================================================================
  subroutine lai_interp(bounds, canopystate_inst)
    !
    ! Interpolate data stream information for Lai.
    !
    ! !USES:
    use pftconMod       , only : noveg
    use CanopyStateType , only : canopystate_type
    use PatchType       , only : patch
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds
    type(canopystate_type) , intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: ivt, p, ip, ig
    character(len=CL)  :: stream_var_name
    !-----------------------------------------------------------------------
    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (lbound(sdat_lai%avs(1)%rAttr,2) <= g_to_ig(bounds%begg) ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(sdat_lai%avs(1)%rAttr,2) >= g_to_ig(bounds%endg) ), sourcefile, __LINE__)

    do p = bounds%begp, bounds%endp
       ivt = patch%itype(p)
       ! Set lai for each gridcell/patch combination
       if (ivt /= noveg) then
          ! vegetated pft
          write(stream_var_name,"(i6)") ivt
          stream_var_name = 'LAI_'//trim(adjustl(stream_var_name))
          ip = mct_aVect_indexRA(sdat_lai%avs(1),trim(stream_var_name))
          ig = g_to_ig(patch%gridcell(p))
          canopystate_inst%tlai_patch(p) = sdat_lai%avs(1)%rAttr(ip,ig)
       else                       
          ! non-vegetated pft
          canopystate_inst%tlai_patch(p) = 0._r8
       endif
    end do

  end subroutine lai_interp

end module LaiStreamMod
