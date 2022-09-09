module FanStreamMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading in FAN nitrogen deposition (in the form of
  ! manure) data file
  ! Also includes functions for fan stream file handling and 
  ! interpolation.
  !
  ! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl
  use clm_varcon, only : ispval
  use shr_strdata_mod
  use shr_stream_mod
  use shr_string_mod
  use shr_sys_mod
  use shr_mct_mod
  use mct_mod
  use spmdMod     , only: mpicom, masterproc, comp_id, iam
  use clm_varctl  , only: iulog
  use abortutils  , only: endrun
  use decompMod   , only: bounds_type
  use domainMod   , only: ldomain
  use ndepStreamMod, only: clm_domain_mct

  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: fanstream_init     ! position datasets for dynamic ndep2
  public :: fanstream_interp   ! interpolates between two years of ndep2 file data
  public :: set_bcast_fanstream_pars
  

  ! ! PRIVATE TYPES
  type(shr_strdata_type)  :: sdat_grz, sdat_sgrz, sdat_ngrz, sdat_urea, sdat_nitr, sdat_soilph         ! input data streams
  integer :: stream_year_first_fan = ispval      ! first year in stream to use
  integer :: stream_year_last_fan = ispval       ! last year in stream to use
  integer :: model_year_align_fan = ispval       ! align stream_year_firstndep2 with 
  character(len=CL)  :: stream_fldFileName_fan
  character(len=CL)  :: fan_mapalgo = 'bilinear'
  logical :: crop_manure_per_crop
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================

  subroutine set_bcast_fanstream_pars(str_yr_first, str_yr_last, mdl_yr_align, mapalgo, str_filename, crop_man_is_percrop)
    integer, intent(in) :: str_yr_first, str_yr_last, mdl_yr_align
    ! whether manure_sgrz and manure_ngrz are per crop or land area:
    logical, intent(in) :: crop_man_is_percrop 
    character(len=*), intent(in) :: str_filename, mapalgo

    stream_year_first_fan = str_yr_first
    stream_year_last_fan = str_yr_last
    model_year_align_fan = mdl_yr_align
    stream_fldFileName_fan = str_filename
    crop_manure_per_crop = crop_man_is_percrop
    fan_mapalgo = mapalgo

    call shr_mpi_bcast(stream_year_first_fan, mpicom)
    call shr_mpi_bcast(stream_year_last_fan, mpicom)
    call shr_mpi_bcast(model_year_align_fan, mpicom)
    call shr_mpi_bcast(stream_fldFileName_fan, mpicom)
    call shr_mpi_bcast(crop_manure_per_crop, mpicom)
    call shr_mpi_bcast(fan_mapalgo, mpicom)
    
  end subroutine set_bcast_fanstream_pars

  !************************************************************************************
  
  subroutine fanstream_init(bounds, NLFilename)
   !    
   ! Initialize data stream information.  
   !
   ! Uses:
   use clm_varctl                , only : inst_name
   use clm_time_manager          , only : get_calendar
   use ncdio_pio                 , only : pio_subsystem
   use shr_pio_mod               , only : shr_pio_getiotype
   use shr_nl_mod                , only : shr_nl_find_group_name
   use shr_log_mod               , only : errMsg => shr_log_errMsg
   use lnd_set_decomp_and_domain , only : gsmap_global
   !
   ! arguments
   implicit none
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm   ! domain information 
   character(*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(*), parameter :: subName = "('fanstream_init')"
   character(*), parameter :: F00 = "('(fanstream_init) ',4a)"
   character(len=80) :: streamvar
   !-----------------------------------------------------------------------

   if (stream_year_first_fan == ispval) then
      call endrun(msg='ERROR stream_year_first_fan not set at '//errMsg(sourcefile, __LINE__))
   end if
   
   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'ndep2dyn stream settings:'
      write(iulog,*) '  stream_year_first_fan  = ',stream_year_first_fan
      write(iulog,*) '  stream_year_last_fan   = ',stream_year_last_fan   
      write(iulog,*) '  model_year_align_fan   = ',model_year_align_fan   
      write(iulog,*) '  stream_fldFileName_fan = ',stream_fldFileName_fan
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   ! Manure N from year-round grazing livestock
   !
   call shr_strdata_create(sdat_grz,name="clmfangrz",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_global, ggrid=dom_clm,          &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='manure_grz',                   &
        fldListModel='manure_grz',                  &
        fillalgo='none',                            &
        mapalgo=fan_mapalgo,                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_grz,'CLMFAN data')
   endif

   ! Manure N from seasonally grazing livestock
   !
   if (crop_manure_per_crop) then
      streamvar = 'manure_sgrz_crop'
   else
      streamvar = 'manure_sgrz'
   end if

   call shr_strdata_create(sdat_sgrz,name="clmfansgrz",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_global, ggrid=dom_clm,          &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile=streamvar,                    &
        fldListModel=streamvar,                   &
        fillalgo='none',                            &
        mapalgo=fan_mapalgo,                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_sgrz,'CLMFAN data')
   endif

   ! Manure N from non-grazing livestock
   !
   if (crop_manure_per_crop) then
      streamvar = 'manure_ngrz_crop'
   else
      streamvar = 'manure_ngrz'
   end if

   call shr_strdata_create(sdat_ngrz,name="clmfanngrz",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_global, ggrid=dom_clm,          &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile=streamvar,                  &
        fldListModel=streamvar,                 &
        fillalgo='none',                            &
        mapalgo=fan_mapalgo,                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_ngrz,'CLMFAN data')
   endif

   ! Fraction of urea N in synthetic fertilizer N
   !
   call shr_strdata_create(sdat_urea,name="clmfanurea",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_global, ggrid=dom_clm,          &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='fract_urea',                &
        fldListModel='fract_urea',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_urea,'CLMFAN data')
   endif

   ! Fraction of nitrate N in synthetic fertilizer N
   !
   call shr_strdata_create(sdat_nitr,name="clmfannitr",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_global, ggrid=dom_clm,          &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='fract_nitr',                &
        fldListModel='fract_nitr',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   ! Soil pH (to be moved to surface dataset)
   ! 
   call shr_strdata_create(sdat_soilph,name="clmfanph",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_global, ggrid=dom_clm,          &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_fan,          &
        yearLast=stream_year_last_fan,            &
        yearAlign=model_year_align_fan,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_fan), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_fan)/),&
        fldListFile='soilph',                &
        fldListModel='fansoilph',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_soilph,'CLMFAN data')
   endif

   
 end subroutine fanstream_init
  
 !================================================================
 subroutine fanstream_interp(bounds, atm2lnd_inst)

   !-----------------------------------------------------------------------
   use clm_time_manager, only : get_curr_date, get_curr_days_per_year
   use clm_varcon      , only : secspday
   use atm2lndType     , only : atm2lnd_type
   use shr_infnan_mod  , only : isinf => shr_infnan_isinf
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
   !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day
   dayspyr = get_curr_days_per_year( )

   call shr_strdata_advance(sdat_grz, mcdate, sec, mpicom, 'clmfangrz')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_ndep_grz_grc(g) = sdat_grz%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
   end do

   call shr_strdata_advance(sdat_sgrz, mcdate, sec, mpicom, 'clmfansgrz')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      if ( isinf(sdat_ngrz%avs(1)%rAttr(1,ig)) )then
         atm2lnd_inst%forc_ndep_sgrz_grc(g) = 9.0e99_r8
      else
         atm2lnd_inst%forc_ndep_sgrz_grc(g) = sdat_sgrz%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
      end if
   end do

   call shr_strdata_advance(sdat_ngrz, mcdate, sec, mpicom, 'clmfanngrz')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      if ( isinf(sdat_ngrz%avs(1)%rAttr(1,ig)) )then
         atm2lnd_inst%forc_ndep_ngrz_grc(g) = 9.0e99_r8
      else
         atm2lnd_inst%forc_ndep_ngrz_grc(g) = sdat_ngrz%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
      end if
   end do
   
   call shr_strdata_advance(sdat_urea, mcdate, sec, mpicom, 'clmfanurea')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_ndep_urea_grc(g) = sdat_urea%avs(1)%rAttr(1,ig)
   end do

   call shr_strdata_advance(sdat_nitr, mcdate, sec, mpicom, 'clmfannitr')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_ndep_nitr_grc(g) = sdat_nitr%avs(1)%rAttr(1,ig)
   end do

   call shr_strdata_advance(sdat_soilph, mcdate, sec, mpicom, 'clmfanph')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_soilph_grc(g) = sdat_soilph%avs(1)%rAttr(1,ig)
   end do

   
 end subroutine fanstream_interp
    
end module FanStreamMod

