module ndep3StreamMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading in FAN nitrogen deposition (in the form of
  ! fertilizer) data file
  ! Also includes functions for dynamic ndep3 file handling and 
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
  public :: ndep3_init     ! position datasets for dynamic ndep3
  public :: ndep3_interp   ! interpolates between two years of ndep3 file data
!KO  public :: clm_domain_mct ! Sets up MCT domain for this resolution

  ! ! PRIVATE TYPES
  type(shr_strdata_type)  :: sdat         ! input data stream
  integer :: stream_year_first_ndep3      ! first year in stream to use
  integer :: stream_year_last_ndep3       ! last year in stream to use
  integer :: model_year_align_ndep3       ! align stream_year_firstndep3 with 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================

  subroutine ndep3_init(bounds, NLFilename)
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
   character(len=CL)  :: stream_fldFileName_ndep3
   character(len=CL)  :: ndep3mapalgo = 'bilinear'
   character(*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(*), parameter :: subName = "('ndep3dyn_init')"
   character(*), parameter :: F00 = "('(ndep3dyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /ndep3dyn_nml/        &
        stream_year_first_ndep3,  &
	stream_year_last_ndep3,   &
        model_year_align_ndep3,   &
        ndep3mapalgo,             &
        stream_fldFileName_ndep3

   ! Default values for namelist
    stream_year_first_ndep3  = 1                ! first year in stream to use
    stream_year_last_ndep3   = 1                ! last  year in stream to use
    model_year_align_ndep3   = 1                ! align stream_year_first_ndep3 with this model year
    stream_fldFileName_ndep3 = ' '

   ! Read ndep3dyn_nml namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, 'ndep3dyn_nml', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=ndep3dyn_nml,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading ndep3dyn_nml namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding ndep3dyn_nml namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_ndep3, mpicom)
   call shr_mpi_bcast(stream_year_last_ndep3, mpicom)
   call shr_mpi_bcast(model_year_align_ndep3, mpicom)
   call shr_mpi_bcast(stream_fldFileName_ndep3, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'ndep3dyn stream settings:'
      write(iulog,*) '  stream_year_first_ndep3  = ',stream_year_first_ndep3
      write(iulog,*) '  stream_year_last_ndep3   = ',stream_year_last_ndep3   
      write(iulog,*) '  model_year_align_ndep3   = ',model_year_align_ndep3   
      write(iulog,*) '  stream_fldFileName_ndep3 = ',stream_fldFileName_ndep3
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(sdat,name="clmndep3",    &
        pio_subsystem=pio_subsystem,                & 
        pio_iotype=shr_pio_getiotype(inst_name),    &
        mpicom=mpicom, compid=comp_id,              &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,     &
        nxg=ldomain%ni, nyg=ldomain%nj,             &
        yearFirst=stream_year_first_ndep3,          &
        yearLast=stream_year_last_ndep3,            &
        yearAlign=model_year_align_ndep3,           &
        offset=0,                                   &
        domFilePath='',                             &
        domFileName=trim(stream_fldFileName_ndep3), &
        domTvarName='time',                         &
        domXvarName='x' ,                         &
        domYvarName='y' ,                         &  
        domAreaName='area',                         &
        domMaskName='mask',                         &
        filePath='',                                &
        filename=(/trim(stream_fldFileName_ndep3)/),&
        fldListFile='NDEP_year',                   &
        fldListModel='NDEP_year',                  &
        fillalgo='none',                            &
        mapalgo=ndep3mapalgo,                       &
        calendar=get_calendar(),                    &
	taxmode='extend'                            )

   if (masterproc) then
      call shr_strdata_print(sdat,'CLMNDEP3 data')
   endif

 end subroutine ndep3_init
  
 !================================================================
 subroutine ndep3_interp(bounds, atm2lnd_inst)

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

   call shr_strdata_advance(sdat, mcdate, sec, mpicom, 'ndep3dyn')

   ig = 0
   dayspyr = get_days_per_year( )
   do g = bounds%begg,bounds%endg
      ig = ig+1
      atm2lnd_inst%forc_ndep3_grc(g) = sdat%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
   end do
   
 end subroutine ndep3_interp

!!==============================================================================
! subroutine clm_domain_mct(bounds, dom_clm)

!   !-------------------------------------------------------------------
!   ! Set domain data type for internal clm grid
!   use clm_varcon  , only : re
!   use domainMod   , only : ldomain
!   use seq_flds_mod
!   implicit none
!   ! 
!   ! arguments
!   type(bounds_type), intent(in) :: bounds  
!   type(mct_ggrid), intent(out)   :: dom_clm     ! Output domain information for land model
!   !
!   ! local variables
!   integer :: g,i,j              ! index
!   integer :: lsize              ! land model domain data size
!   real(r8), pointer :: data(:)  ! temporary
!   integer , pointer :: idata(:) ! temporary
!   !-------------------------------------------------------------------
!   !
!   ! Initialize mct domain type
!   ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
!   ! Note that in addition land carries around landfrac for the purposes of domain checking
!   ! 
!   lsize = mct_gsMap_lsize(gsmap_lnd_gdc2glo, mpicom)
!   call mct_gGrid_init( GGrid=dom_clm, CoordChars=trim(seq_flds_dom_coord), &
!                        OtherChars=trim(seq_flds_dom_other), lsize=lsize )
!   !
!   ! Allocate memory
!   !
!   allocate(data(lsize))
!   !
!   ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
!   !
!   call mct_gsMap_orderedPoints(gsmap_lnd_gdc2glo, iam, idata)
!   call mct_gGrid_importIAttr(dom_clm,'GlobGridNum',idata,lsize)
!   !
!   ! Determine domain (numbering scheme is: West to East and South to North to South pole)
!   ! Initialize attribute vector with special value
!   !
!   data(:) = -9999.0_R8 
!   call mct_gGrid_importRAttr(dom_clm,"lat"  ,data,lsize) 
!   call mct_gGrid_importRAttr(dom_clm,"lon"  ,data,lsize) 
!   call mct_gGrid_importRAttr(dom_clm,"area" ,data,lsize) 
!   call mct_gGrid_importRAttr(dom_clm,"aream",data,lsize) 
!   data(:) = 0.0_R8     
!   call mct_gGrid_importRAttr(dom_clm,"mask" ,data,lsize) 
!   !
!   ! Determine bounds
!   !
!   ! Fill in correct values for domain components
!   ! Note aream will be filled in in the atm-lnd mapper
!   !
!   do g = bounds%begg,bounds%endg
!      i = 1 + (g - bounds%begg)
!      data(i) = ldomain%lonc(g)
!   end do
!   call mct_gGrid_importRattr(dom_clm,"lon",data,lsize) 

!   do g = bounds%begg,bounds%endg
!      i = 1 + (g - bounds%begg)
!      data(i) = ldomain%latc(g)
!   end do
!   call mct_gGrid_importRattr(dom_clm,"lat",data,lsize) 

!   do g = bounds%begg,bounds%endg
!      i = 1 + (g - bounds%begg)
!      data(i) = ldomain%area(g)/(re*re)
!   end do
!   call mct_gGrid_importRattr(dom_clm,"area",data,lsize) 

!   do g = bounds%begg,bounds%endg
!      i = 1 + (g - bounds%begg)
!      data(i) = real(ldomain%mask(g), r8)
!   end do
!   call mct_gGrid_importRattr(dom_clm,"mask",data,lsize) 

!   do g = bounds%begg,bounds%endg
!      i = 1 + (g - bounds%begg)
!      data(i) = real(ldomain%frac(g), r8)
!   end do
!   call mct_gGrid_importRattr(dom_clm,"frac",data,lsize) 

!   deallocate(data)
!   deallocate(idata)

! end subroutine clm_domain_mct
    
end module ndep3StreamMod

