module ndepStreamMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in nitrogen deposition data file
  ! Also includes functions for dynamic ndep file handling and
  ! interpolation.
  !
  ! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_strdata_mod, only: shr_strdata_type, shr_strdata_create
  use shr_strdata_mod, only: shr_strdata_print, shr_strdata_advance
  use mct_mod     , only: mct_ggrid
  use spmdMod     , only: mpicom, masterproc, comp_id, iam
  use clm_varctl  , only: iulog, inst_name
  use abortutils  , only: endrun
  use decompMod   , only: bounds_type
  use domainMod   , only: ldomain

  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ndep_init      ! position datasets for dynamic ndep
  public :: ndep_interp    ! interpolates between two years of ndep file data
  public :: clm_domain_mct ! Sets up MCT domain for this resolution

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: check_units   ! Check the units and make sure they can be used

  ! ! PRIVATE TYPES
  type(shr_strdata_type)  :: sdat           ! input data stream
  integer :: stream_year_first_ndep         ! first year in stream to use
  integer :: stream_year_last_ndep          ! last year in stream to use
  integer :: model_year_align_ndep          ! align stream_year_firstndep with
  logical :: divide_by_secs_per_yr = .true. ! divide by the number of seconds per year

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
   use shr_kind_mod     , only : CS => shr_kind_cs
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   use lnd_set_decomp_and_domain , only : gsMap_lnd2Dsoi_gdc2glo, gsmap_global
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
   character(len=CL)  :: stream_fldFileName_ndep
   character(len=CL)  :: ndepmapalgo = 'bilinear'
   character(len=CL)  :: ndep_tintalgo = 'linear'
   character(len=CS)  :: ndep_taxmode = 'extend'
   character(len=CL)  :: ndep_varlist = 'NDEP_year'
   character(*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(*), parameter :: subName = "('ndepdyn_init')"
   character(*), parameter :: F00 = "('(ndepdyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /ndepdyn_nml/          &
        stream_year_first_ndep,    &
        stream_year_last_ndep,     &
        model_year_align_ndep,     &
        ndepmapalgo, ndep_taxmode, &
        ndep_varlist,              &
        stream_fldFileName_ndep,   &
        ndep_tintalgo

   ! Default values for namelist
    stream_year_first_ndep  = 1                ! first year in stream to use
    stream_year_last_ndep   = 1                ! last  year in stream to use
    model_year_align_ndep   = 1                ! align stream_year_first_ndep with this model year
    stream_fldFileName_ndep = ' '

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

   call shr_mpi_bcast(stream_year_first_ndep , mpicom)
   call shr_mpi_bcast(stream_year_last_ndep  , mpicom)
   call shr_mpi_bcast(model_year_align_ndep  , mpicom)
   call shr_mpi_bcast(stream_fldFileName_ndep, mpicom)
   call shr_mpi_bcast(ndep_varlist           , mpicom)
   call shr_mpi_bcast(ndep_taxmode           , mpicom)
   call shr_mpi_bcast(ndep_tintalgo          , mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'ndepdyn stream settings:'
      write(iulog,*) '  stream_year_first_ndep  = ',stream_year_first_ndep
      write(iulog,*) '  stream_year_last_ndep   = ',stream_year_last_ndep
      write(iulog,*) '  model_year_align_ndep   = ',model_year_align_ndep
      write(iulog,*) '  stream_fldFileName_ndep = ',stream_fldFileName_ndep
      write(iulog,*) '  ndep_varList            = ',ndep_varList
      write(iulog,*) '  ndep_taxmode            = ',ndep_taxmode
      write(iulog,*) '  ndep_tintalgo           = ',ndep_tintalgo
      write(iulog,*) ' '
   endif
   ! Read in units
   call check_units( stream_fldFileName_ndep, ndep_varList )

   ! Set domain and create streams
   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(sdat,name="clmndep",    &
        pio_subsystem=pio_subsystem,               &
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_global, ggrid=dom_clm,         &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_ndep,          &
        yearLast=stream_year_last_ndep,            &
        yearAlign=model_year_align_ndep,           &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_ndep), &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_ndep)/),&
        fldListFile=ndep_varlist,                  &
        fldListModel=ndep_varlist,                 &
        fillalgo='none',                           &
        mapalgo=ndepmapalgo,                       &
        tintalgo=ndep_tintalgo,                    &
        calendar=get_calendar(),                   &
        taxmode=ndep_taxmode                       )


   if (masterproc) then
      call shr_strdata_print(sdat,'CLMNDEP data')
   endif

 end subroutine ndep_init
 !================================================================

 subroutine check_units( stream_fldFileName_ndep, ndep_varList )
   !-------------------------------------------------------------------
   ! Check that units are correct on the file and if need any conversion
   use ncdio_pio     , only : ncd_pio_openfile, ncd_inqvid, ncd_getatt, ncd_pio_closefile, ncd_nowrite
   use ncdio_pio     , only : file_desc_t, var_desc_t
   use shr_kind_mod  , only : CS => shr_kind_cs
   use shr_log_mod   , only : errMsg => shr_log_errMsg
   use shr_string_mod, only : shr_string_listGetName
   implicit none

   !-----------------------------------------------------------------------
   !
   ! Arguments
   character(len=*), intent(IN)  :: stream_fldFileName_ndep  ! ndep filename
   character(len=*), intent(IN)  :: ndep_varList             ! ndep variable list to examine
   !
   ! Local variables
   type(file_desc_t) :: ncid     ! NetCDF filehandle for ndep file
   type(var_desc_t)  :: vardesc  ! variable descriptor
   integer           :: varid    ! variable index
   logical           :: readvar  ! If variable was read
   character(len=CS) :: ndepunits! ndep units
   character(len=CS) :: fname    ! ndep field name
   !-----------------------------------------------------------------------
   call ncd_pio_openfile( ncid, trim(stream_fldFileName_ndep), ncd_nowrite )
   call shr_string_listGetName( ndep_varList, 1, fname )
   call ncd_inqvid(       ncid, fname, varid, vardesc, readvar=readvar )
   if ( readvar ) then
      call ncd_getatt(    ncid, varid, "units", ndepunits )
   else
      call endrun(msg=' ERROR finding variable: '//trim(fname)//" in file: "// &
                      trim(stream_fldFileName_ndep)//errMsg(sourcefile, __LINE__))
   end if
   call ncd_pio_closefile( ncid )

   ! Now check to make sure they are correct
   if (      trim(ndepunits) == "g(N)/m2/s"  )then
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
   use clm_time_manager, only : get_curr_date, get_curr_days_per_year
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

   call shr_strdata_advance(sdat, mcdate, sec, mpicom, 'ndepdyn')

   if ( divide_by_secs_per_yr )then
      ig = 0
      dayspyr = get_curr_days_per_year( )
      do g = bounds%begg,bounds%endg
         ig = ig+1
         atm2lnd_inst%forc_ndep_grc(g) = sdat%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
      end do
   else
      ig = 0
      do g = bounds%begg,bounds%endg
         ig = ig+1
         atm2lnd_inst%forc_ndep_grc(g) = sdat%avs(1)%rAttr(1,ig)
      end do
   end if

 end subroutine ndep_interp

 !==============================================================================
  subroutine clm_domain_mct(bounds, dom_clm, nlevels)

    !-------------------------------------------------------------------
    ! Set domain data type for internal clm grid
    use clm_varcon                , only : re
    use domainMod                 , only : ldomain
    use mct_mod                   , only : mct_ggrid, mct_gsMap_lsize, mct_gGrid_init
    use mct_mod                   , only : mct_gsMap_orderedPoints, mct_gGrid_importIAttr
    use mct_mod                   , only : mct_gGrid_importRAttr, mct_gsMap
    use lnd_set_decomp_and_domain , only : gsMap_lnd2Dsoi_gdc2glo, gsmap_global
    implicit none
    !
    ! arguments
    type(bounds_type), intent(in) :: bounds
    type(mct_ggrid), intent(out)   :: dom_clm     ! Output domain information for land model
    integer, intent(in), optional :: nlevels      ! Number of levels if this is a 3D field
    !
    ! local variables
    integer :: g,i,j,k            ! index
    integer :: lsize              ! land model domain data size
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    integer :: nlevs              ! Number of vertical levels
    type(mct_gsMap), pointer :: gsmap => null() ! MCT GS map
    !-------------------------------------------------------------------
    ! SEt number of levels, and get the GS map for either the 2D or 3D grid
    nlevs = 1
    if ( present(nlevels) ) nlevs = nlevels
    if ( nlevs == 1 ) then
       gsmap => gsmap_global
    else
       gsmap => gsMap_lnd2Dsoi_gdc2glo
    end if
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    !
    lsize = mct_gsMap_lsize(gsmap, mpicom)
    call mct_gGrid_init( GGrid=dom_clm, &
         CoordChars='lat:lon:hgt', OtherChars='area:aream:mask:frac', lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsmap, iam, idata)
    gsmap => null()
    call mct_gGrid_importIAttr(dom_clm,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(dom_clm,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_clm,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_clm,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_clm,"aream",data,lsize)
    data(:) = 0.0_R8
    call mct_gGrid_importRAttr(dom_clm,"mask" ,data,lsize)
    !
    ! Determine bounds
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%lonc(g)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"lon",data,lsize)

    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%latc(g)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"lat",data,lsize)

    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%area(g)/(re*re)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"area",data,lsize)

    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%mask(g), r8)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"mask",data,lsize)

    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%frac(g), r8)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"frac",data,lsize)

    deallocate(data)
    deallocate(idata)

  end subroutine clm_domain_mct

end module ndepStreamMod
