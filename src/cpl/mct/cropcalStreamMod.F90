module cropcalStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read crop calendars from stream
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
  public :: cropcal_init    ! position datasets for crop calendars
  public :: cropcal_advance ! Advance the crop calendar streams (outside of a Open-MP threading loop)
  public :: cropcal_interp  ! interpolates between two years of crop calendar data

  ! !PRIVATE MEMBER DATA:
  integer, allocatable    :: g_to_ig(:)         ! Array matching gridcell index to data index
  ! SSR TODO: Make these work with max_growingseasons_per_year > 1
  type(shr_strdata_type)  :: sdat_sdate        ! sowing date input data stream

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
    use clm_varpar       , only : cft_lb, cft_ub
    use clm_time_manager , only : get_calendar
    use ncdio_pio        , only : pio_subsystem
    use shr_pio_mod      , only : shr_pio_getiotype
    use shr_stream_mod   , only : shr_stream_file_null
    use shr_string_mod   , only : shr_string_listCreateField_range
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use ndepStreamMod    , only : clm_domain_mct
    use histFileMod      , only : hist_addfld1d
    use domainMod        , only : ldomain
    use controlMod       , only : NLFilename
    use lnd_set_decomp_and_domain , only : gsmap_global
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds          ! bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: stream_year_first_cropcal      ! first year in crop calendar streams to use
    integer            :: stream_year_last_cropcal       ! last year in crop calendar streams to use
    integer            :: model_year_align_cropcal       ! align stream_year_first_cropcal with
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    type(mct_ggrid)    :: dom_clm                    ! domain information
    character(len=CL)  :: stream_fldFileName_sdate     ! sowing date stream filename to read
    character(len=CL)  :: cropcal_mapalgo  = 'nn'      ! Mapping alogrithm
    character(len=CL)  :: cropcal_tintalgo = 'nearest' ! Time interpolation alogrithm
    
    ! SSR TODO: Make this work with max_growingseasons_per_year > 1
    character(len=CXX) :: fldList_sdate1                  ! field string for 1st sowing dates
    
    character(*), parameter :: subName = "('cropcaldyn_init')"
    !-----------------------------------------------------------------------
    !
    ! deal with namelist variables here in init
    !
    namelist /cropcal_streams/         &
         stream_year_first_cropcal,    &
         stream_year_last_cropcal,     &
         model_year_align_cropcal,     &
         stream_fldFileName_sdate

    ! Default values for namelist
    stream_year_first_cropcal     = 1      ! first year in stream to use
    stream_year_last_cropcal      = 1      ! last  year in stream to use
    model_year_align_cropcal      = 1      ! align stream_year_first_cropcal with this model year
    stream_fldFileName_sdate      = shr_stream_file_null

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
    call shr_mpi_bcast(stream_year_first_cropcal, mpicom)
    call shr_mpi_bcast(stream_year_last_cropcal , mpicom)
    call shr_mpi_bcast(model_year_align_cropcal , mpicom)
    call shr_mpi_bcast(stream_fldFileName_sdate , mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'cropcal_stream settings:'
       write(iulog,*) '  stream_year_first_cropcal  = ',stream_year_first_cropcal
       write(iulog,*) '  stream_year_last_cropcal   = ',stream_year_last_cropcal
       write(iulog,*) '  model_year_align_cropcal   = ',model_year_align_cropcal
       write(iulog,*) '  stream_fldFileName_sdate = ',trim(stream_fldFileName_sdate)
    endif

    call clm_domain_mct (bounds, dom_clm)

    ! create the field list for these cropcal fields...use in shr_strdata_create
    ! SSR TODO: Make this work with max_growingseasons_per_year > 1
    write(iulog,'(a,i6)') '        cft_lb = ',cft_lb
    write(iulog,'(a,i6)') '        cft_ub = ',cft_ub
    fldList_sdate1 = shr_string_listCreateField_range( cft_lb, cft_ub, "sdate1" )
    write(iulog,*) 'fldList_sdate1 = ',trim(fldList_sdate1)

    ! SSR TODO:
    ! - Delete "area" and "mask"?
    ! - Is this correct taxmode?
    call shr_strdata_create(sdat_sdate,               &
         name="cropcaldyn",                            &
         pio_subsystem=pio_subsystem,                  &
         pio_iotype=shr_pio_getiotype(inst_name),      &
         mpicom=mpicom, compid=comp_id,                &
         gsmap=gsmap_global, ggrid=dom_clm,            &
         nxg=ldomain%ni, nyg=ldomain%nj,               &
         yearFirst=stream_year_first_cropcal,          &
         yearLast=stream_year_last_cropcal,            &
         yearAlign=model_year_align_cropcal,           &
         offset=0,                                     &
         domFilePath='',                               &
         domFileName=trim(stream_fldFileName_sdate),   &
         domTvarName='time',                           &
         domXvarName='lon' ,                           &
         domYvarName='lat' ,                           &
         domAreaName='area',                           &
         domMaskName='mask',                           &
         filePath='',                                  &
         filename=(/stream_fldFileName_sdate/),        &
         fldListFile=fldList_sdate1,                   &
         fldListModel=fldList_sdate1,                  &
         fillalgo='none',                              &
         mapalgo=cropcal_mapalgo,                      &
         tintalgo=cropcal_tintalgo,                    &
         calendar=get_calendar(),                      &
         taxmode='cycle'                               )

    if (masterproc) then
       call shr_strdata_print(sdat_sdate,'sdate1 data')
    endif

  end subroutine cropcal_init

  !==============================================================================
  subroutine cropcal_advance( bounds )
    !
    ! Advance crop calendar streams
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

    ! SSR TODO: Make this work with max_growingseasons_per_year > 1
    call shr_strdata_advance(sdat_sdate, mcdate, sec, mpicom, 'cropcaldyn')

    if ( .not. allocated(g_to_ig) )then
       allocate (g_to_ig(bounds%begg:bounds%endg) )
       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          g_to_ig(g) = ig
       end do
    end if

  end subroutine cropcal_advance

  !==============================================================================
  subroutine cropcal_interp(bounds, num_pcropp, filter_pcropp, crop_inst)
    !
    ! Interpolate data stream information for crop calendars.
    !
    ! !USES:
    use pftconMod       , only : noveg
    use CropType        , only : crop_type
    use PatchType       , only : patch
    use filterMod       , only : filter
    use decompMod       , only : get_proc_clumps
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
    integer :: patch_gridcell_g
    character(len=CL)  :: stream_var_name
    !-----------------------------------------------------------------------
    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)

    ! SSR TODO: Make this work with max_growingseasons_per_year > 1
    write(iulog,*) 'cropcal_interp() A'
    SHR_ASSERT_FL( (lbound(sdat_sdate%avs(1)%rAttr,2) <= g_to_ig(bounds%begg) ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(sdat_sdate%avs(1)%rAttr,2) >= g_to_ig(bounds%endg) ), sourcefile, __LINE__)

    ! SSR TODO: Make these work with max_growingseasons_per_year > 1
    write(iulog,*) 'cropcal_interp() B'
    do fp = 1, num_pcropp
       write(iulog,*) 'cropcal_interp() C'
       p = filter_pcropp(fp)
       ivt = patch%itype(p)
       ! Set crop calendars for each gridcell/patch combination
       write(stream_var_name,"(i6)") ivt
       
       ! SSR TODO: Add check that variable exists in netCDF
       write(iulog,*) 'cropcal_interp() D1'
       stream_var_name = 'sdate1_'//trim(adjustl(stream_var_name))
       write(iulog,*) 'cropcal_interp() D2'
       ip = mct_aVect_indexRA(sdat_sdate%avs(1),trim(stream_var_name))
       if (ivt /= noveg) then
          write(iulog,*) 'cropcal_interp() D3'
          write(iulog,'(a,i6)')'                    p = ', p
          write(iulog,'(a,i6)')'                  ivt = ', ivt
          write(iulog,'(a,i6)')'    patch%gridcell(p) = ', patch%gridcell(p)
          write(iulog,'(a,i6)')'       g_to_ig(above) = ', g_to_ig(patch%gridcell(p))
          ig = g_to_ig(patch%gridcell(p))
          write(iulog,*) 'cropcal_interp() D4'
          crop_inst%sdates_thisyr(p,1) = sdat_sdate%avs(1)%rAttr(ip,ig)
       endif

       ! SSR TODO: Make this work with max_growingseasons_per_year > 1
       crop_inst%n_growingseasons_thisyear_thispatch(p) = crop_inst%sdates_thisyr(p,1) >= 0

       ! Only for first sowing date of the year
       crop_inst%next_rx_sdate(p) = crop_inst%sdates_thisyr(p,1)

    end do

  end subroutine cropcal_interp

end module cropcalStreamMod
