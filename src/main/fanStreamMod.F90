module FanStreamMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading in FAN nitrogen deposition (in the form of
  ! manure) data file
  ! Also includes functions for dynamic ndep2 file handling and 
  ! interpolation.
  !
  ! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_strdata_mod
  use shr_stream_mod
  use shr_string_mod
  use shr_sys_mod
  use shr_mct_mod
  use mct_mod
  use spmdMod     , only: mpicom, masterproc, comp_id, iam
  use clm_varctl  , only: iulog
  use abortutils  , only: endrun
  use fileutils   , only: getavu, relavu
  use decompMod   , only: bounds_type, ldecomp, gsmap_lnd_gdc2glo 
  use domainMod   , only: ldomain
!KO
  use ndepStreamMod, only: clm_domain_mct
!KO

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: fanstream_init     ! position datasets for dynamic ndep2
  public :: fanstream_interp   ! interpolates between two years of ndep2 file data
!KO  public :: clm_domain_mct ! Sets up MCT domain for this resolution

  ! ! PRIVATE TYPES
  type(shr_strdata_type)  :: sdat_past, sdat_mix, sdat_urea, sdat_nitr, sdat_soilph         ! input data streams
  integer :: stream_year_first_ndep2      ! first year in stream to use
  integer :: stream_year_last_ndep2       ! last year in stream to use
  integer :: model_year_align_ndep2       ! align stream_year_firstndep2 with 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================

  subroutine fanstream_init(bounds, NLFilename)
   !    
   ! Initialize data stream information.  
   !
   ! Uses:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
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
   character(len=CL)  :: stream_fldFileName_ndep2
   character(len=CL)  :: ndep2mapalgo = 'bilinear'
   character(*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(*), parameter :: subName = "('ndep2dyn_init')"
   character(*), parameter :: F00 = "('(ndep2dyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /ndep2dyn_nml/        &
        stream_year_first_ndep2,  &
	stream_year_last_ndep2,   &
        model_year_align_ndep2,   &
        ndep2mapalgo,             &
        stream_fldFileName_ndep2

   ! Default values for namelist
    stream_year_first_ndep2  = 1                ! first year in stream to use
    stream_year_last_ndep2   = 1                ! last  year in stream to use
    model_year_align_ndep2   = 1                ! align stream_year_first_ndep2 with this model year
    stream_fldFileName_ndep2 = ' '

   ! Read ndep2dyn_nml namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, 'ndep2dyn_nml', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=ndep2dyn_nml,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading ndep2dyn_nml namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding ndep2dyn_nml namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_ndep2, mpicom)
   call shr_mpi_bcast(stream_year_last_ndep2, mpicom)
   call shr_mpi_bcast(model_year_align_ndep2, mpicom)
   call shr_mpi_bcast(stream_fldFileName_ndep2, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'ndep2dyn stream settings:'
      write(iulog,*) '  stream_year_first_ndep2  = ',stream_year_first_ndep2
      write(iulog,*) '  stream_year_last_ndep2   = ',stream_year_last_ndep2   
      write(iulog,*) '  model_year_align_ndep2   = ',model_year_align_ndep2   
      write(iulog,*) '  stream_fldFileName_ndep2 = ',stream_fldFileName_ndep2
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(sdat_past,name="clmndep2past",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_ndep2,          &
        yearLast=stream_year_last_ndep2,            &
        yearAlign=model_year_align_ndep2,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_ndep2), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_ndep2)/),&
        fldListFile='Nmanure_pastures',             &
        fldListModel='Nmanure_pastures',            &
        fillalgo='none',                            &
        mapalgo=ndep2mapalgo,                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_past,'CLMNDEP2 data')
   endif

   call shr_strdata_create(sdat_mix,name="clmndep2mixed",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_ndep2,          &
        yearLast=stream_year_last_ndep2,            &
        yearAlign=model_year_align_ndep2,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_ndep2), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_ndep2)/),&
        fldListFile='Nmanure_mixed',                &
        fldListModel='Nmanure_mixed',               &
        fillalgo='none',                            &
        mapalgo=ndep2mapalgo,                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_mix,'CLMNDEP2 data')
   endif

   call shr_strdata_create(sdat_urea,name="clmndep2urea",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_ndep2,          &
        yearLast=stream_year_last_ndep2,            &
        yearAlign=model_year_align_ndep2,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_ndep2), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_ndep2)/),&
        fldListFile='fract_urea',                &
        fldListModel='fract_urea',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_urea,'CLMNDEP2 data')
   endif

   call shr_strdata_create(sdat_nitr,name="clmndep2nitr",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_ndep2,          &
        yearLast=stream_year_last_ndep2,            &
        yearAlign=model_year_align_ndep2,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_ndep2), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_ndep2)/),&
        fldListFile='fract_nitr',                &
        fldListModel='fract_nitr',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   call shr_strdata_create(sdat_soilph,name="clmndep2ph",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_ndep2,          &
        yearLast=stream_year_last_ndep2,            &
        yearAlign=model_year_align_ndep2,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_ndep2), &
        domTvarName='time',                         &
        domXvarName='x' ,                           &
        domYvarName='y' ,                           &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_ndep2)/),&
        fldListFile='soilph',                &
        fldListModel='fansoilph',               &
        fillalgo='none',                            &
        mapalgo='nn',                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat_mix,'CLMNDEP2 data')
   endif

   
 end subroutine fanstream_init
  
 !================================================================
 subroutine fanstream_interp(bounds, atm2lnd_inst)

   !-----------------------------------------------------------------------
   use clm_time_manager, only : get_curr_date, get_days_per_year
   use clm_varcon      , only : secspday
   use atm2lndType     , only : atm2lnd_type
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
   dayspyr = get_days_per_year( )

   call shr_strdata_advance(sdat_past, mcdate, sec, mpicom, 'clmndep2pasture')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_ndep3_grc(g) = sdat_past%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
   end do

   call shr_strdata_advance(sdat_mix, mcdate, sec, mpicom, 'clmndep2mixed')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_ndep2_grc(g) = sdat_mix%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
   end do

   call shr_strdata_advance(sdat_urea, mcdate, sec, mpicom, 'clmndep2urea')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_ndep_urea_grc(g) = sdat_urea%avs(1)%rAttr(1,ig)
   end do

   call shr_strdata_advance(sdat_nitr, mcdate, sec, mpicom, 'clmndep2nitr')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_ndep_nitr_grc(g) = sdat_nitr%avs(1)%rAttr(1,ig)
   end do

   call shr_strdata_advance(sdat_soilph, mcdate, sec, mpicom, 'clmndep2ph')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_soilph_grc(g) = sdat_soilph%avs(1)%rAttr(1,ig)
   end do

   
 end subroutine fanstream_interp
    
end module FanStreamMod

