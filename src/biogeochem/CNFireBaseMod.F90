module CNFireBaseMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for fire dynamics 
  ! created in Nov, 2012  and revised in Apr, 2013 by F. Li and S. Levis
  ! based on Li et al. (2012a,b; 2013)
  ! revised in Apr, 2014 according Li et al.(2014)
  ! Fire-related parameters were calibrated or tuned in Apr, 2013 based on the 
  ! 20th Century transient simulations at f19_g16 with (newfire05_clm45sci15_clm4_0_58) 
  ! a CLM4.5 version, Qian et al. (2006) atmospheric forcing, and
  ! climatological lightning data.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_strdata_mod                    , only : shr_strdata_type, shr_strdata_create, shr_strdata_print
  use shr_strdata_mod                    , only : shr_strdata_advance
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varctl                         , only : iulog
  use spmdMod                            , only : masterproc, mpicom, comp_id
  use fileutils                          , only : getavu, relavu
  use decompMod                          , only : gsmap_lnd_gdc2glo
  use domainMod                          , only : ldomain
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use atm2lndType                        , only : atm2lnd_type
  use CNDVType                           , only : dgvs_type
  use CNVegStateType                     , only : cnveg_state_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType             , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType              , only : cnveg_nitrogenflux_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use EnergyFluxType                     , only : energyflux_type
  use SoilHydrologyType                  , only : soilhydrology_type  
  use WaterstateType                     , only : waterstate_type
  use mct_mod
  use CNFireMethodMod                    , only : cnfire_method_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: cnfire_base_type
  !
  type, extends(cnfire_method_type) :: cnfire_base_type
    private
      ! !PRIVATE MEMBER DATA:
      real(r8), public, pointer :: forc_lnfm(:)        ! Lightning frequency
      real(r8), public, pointer :: forc_hdm(:)         ! Human population density

      type(shr_strdata_type) :: sdat_hdm    ! Human population density input data stream
      type(shr_strdata_type) :: sdat_lnfm   ! Lightning input data stream
    contains
      !
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: CNFireInit        ! Initialization of CNFire
      procedure, public :: CNFireInterp      ! Interpolate fire data
      procedure, public :: CNFireArea        ! Calculate fire area
      procedure, public :: CNFireFluxes      ! Calculate fire fluxes
      !
      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: hdm_init    ! position datasets for dynamic human population density
      procedure, private :: hdm_interp  ! interpolates between two years of human pop. density file data
      procedure, private :: lnfm_init   ! position datasets for Lightning
      procedure, private :: lnfm_interp ! interpolates between two years of Lightning file data
  end type cnfire_base_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNFireInit( this, bounds, NLFilename )
    !
    ! !DESCRIPTION:
    ! Initialize CN Fire module
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(cnfire_base_type) :: this
    type(bounds_type), intent(in) :: bounds  
    character(len=*),  intent(in) :: NLFilename
    !-----------------------------------------------------------------------

    ! Allocate lightning forcing data
    allocate( this%forc_lnfm(bounds%begg:bounds%endg) )
    this%forc_lnfm(bounds%begg:) = nan
    ! Allocate pop dens forcing data
    allocate( this%forc_hdm(bounds%begg:bounds%endg) )
    this%forc_hdm(bounds%begg:) = nan

    call this%hdm_init(bounds, NLFilename)
    call this%hdm_interp(bounds)
    call this%lnfm_init(bounds, NLFilename)
    call this%lnfm_interp(bounds)

  end subroutine CNFireInit

  !-----------------------------------------------------------------------
  subroutine CNFireInterp(this,bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate CN Fire datasets
    !
    ! !ARGUMENTS:
    class(cnfire_base_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !-----------------------------------------------------------------------

    call this%hdm_interp(bounds)
    call this%lnfm_interp(bounds)

  end subroutine CNFireInterp

  !-----------------------------------------------------------------------
  subroutine CNFireArea (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       atm2lnd_inst, energyflux_inst, soilhydrology_inst, waterstate_inst, &
       cnveg_state_inst, cnveg_carbonstate_inst, totlitc_col, decomp_cpools_vr_col, t_soi17cm_col)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cnfire_base_type)                               :: this
    type(bounds_type)                     , intent(in)    :: bounds 
    integer                               , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(atm2lnd_type)                    , intent(in)    :: atm2lnd_inst
    type(energyflux_type)                 , intent(in)    :: energyflux_inst
    type(soilhydrology_type)              , intent(in)    :: soilhydrology_inst
    type(waterstate_type)                 , intent(in)    :: waterstate_inst
    type(cnveg_state_type)                , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)          , intent(inout) :: cnveg_carbonstate_inst
    real(r8)                              , intent(in)    :: totlitc_col(bounds%begc:)
    real(r8)                              , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:)
    real(r8)                              , intent(in)    :: t_soi17cm_col(bounds%begc:)
    !

    call endrun( 'cnfire_base::CNFireArea: this method MUST be implemented!' )

  end subroutine CNFireArea

  !-----------------------------------------------------------------------
  subroutine CNFireFluxes (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      dgvs_inst, cnveg_state_inst,                                                                      &
      cnveg_carbonstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
      leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch, &
      totsomc_col, decomp_cpools_vr_col, decomp_npools_vr_col, somc_fire_col)
   !
   ! !DESCRIPTION:
   ! Fire effects routine for coupled carbon-nitrogen code (CN).
   ! Relies primarily on estimate of fractional area burned, from CNFireArea().
   !
   ! Total fire carbon emissions (g C/m2 land area/yr) 
   !  =avg(COL_FIRE_CLOSS)*seconds_per_year + avg(SOMC_FIRE)*seconds_per_year + 
   !   avg(LF_CONV_CFLUX)*seconds_per_year*min(1.0,avg(LFC2)*seconds_per_year)*0.8
   ! where avg means the temporal average in a year
   ! seconds_per_year is the number of seconds in a year.
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   class(cnfire_base_type)                        :: this
   type(bounds_type)              , intent(in)    :: bounds  
   integer                        , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
   integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
   integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
   type(dgvs_type)                , intent(inout) :: dgvs_inst
   type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
   type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
   type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
   type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
   real(r8)                       , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: totsomc_col(bounds%begc:)                ! (gC/m2) total soil organic matter C
   real(r8)                       , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:) ! (gC/m3)  VR decomp. (litter, cwd, soil)
   real(r8)                       , intent(in)    :: decomp_npools_vr_col(bounds%begc:,1:,1:) ! (gC/m3)  VR decomp. (litter, cwd, soil)
   real(r8)                       , intent(out)   :: somc_fire_col(bounds%begc:)              ! (gC/m2/s) fire C emissions due to peat burning
   !
   !-----------------------------------------------------------------------
   somc_fire_col(bounds%begc:) = 1.0_r8
   call endrun( 'cnfire_base::CNFireFluxes: this method MUST be implemented!' )

  end subroutine CNFireFluxes

  !-----------------------------------------------------------------------
  subroutine hdm_init( this, bounds, NLFilename )
   !
   ! !DESCRIPTION:
   ! Initialize data stream information for population density.
   !
   ! !USES:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use clm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : clm_domain_mct
   use histFileMod      , only : hist_addfld1d
   !
   ! !ARGUMENTS:
   implicit none
   class(cnfire_base_type)       :: this
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! !LOCAL VARIABLES:
   integer            :: stream_year_first_popdens   ! first year in pop. dens. stream to use
   integer            :: stream_year_last_popdens    ! last year in pop. dens. stream to use
   integer            :: model_year_align_popdens    ! align stream_year_first_hdm with 
   integer            :: nu_nml                      ! unit for namelist file
   integer            :: nml_error                   ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm                     ! domain information 
   character(len=CL)  :: stream_fldFileName_popdens  ! population density streams filename
   character(len=CL)  :: popdensmapalgo = 'bilinear' ! mapping alogrithm for population density
   character(*), parameter :: subName = "('hdmdyn_init')"
   character(*), parameter :: F00 = "('(hdmdyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /popd_streams/          &
        stream_year_first_popdens,  &
        stream_year_last_popdens,   &
        model_year_align_popdens,   &
        popdensmapalgo,             &
        stream_fldFileName_popdens

   ! Default values for namelist
   stream_year_first_popdens  = 1       ! first year in stream to use
   stream_year_last_popdens   = 1       ! last  year in stream to use
   model_year_align_popdens   = 1       ! align stream_year_first_popdens with this model year
   stream_fldFileName_popdens = ' '

   ! Read popd_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'popd_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=popd_streams,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading popd_streams namelist'//errMsg(__FILE__, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_popdens, mpicom)
   call shr_mpi_bcast(stream_year_last_popdens, mpicom)
   call shr_mpi_bcast(model_year_align_popdens, mpicom)
   call shr_mpi_bcast(stream_fldFileName_popdens, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'popdens_streams settings:'
      write(iulog,*) '  stream_year_first_popdens  = ',stream_year_first_popdens  
      write(iulog,*) '  stream_year_last_popdens   = ',stream_year_last_popdens   
      write(iulog,*) '  model_year_align_popdens   = ',model_year_align_popdens   
      write(iulog,*) '  stream_fldFileName_popdens = ',stream_fldFileName_popdens
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(this%sdat_hdm,name="clmhdm",     &
        pio_subsystem=pio_subsystem,                   & 
        pio_iotype=shr_pio_getiotype(inst_name),       &
        mpicom=mpicom, compid=comp_id,                 &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,        &
        nxg=ldomain%ni, nyg=ldomain%nj,                &
        yearFirst=stream_year_first_popdens,           &
        yearLast=stream_year_last_popdens,             &
        yearAlign=model_year_align_popdens,            &
        offset=0,                                      &
        domFilePath='',                                &
        domFileName=trim(stream_fldFileName_popdens),  &
        domTvarName='time',                            &
        domXvarName='lon' ,                            &
        domYvarName='lat' ,                            &  
        domAreaName='area',                            &
        domMaskName='mask',                            &
        filePath='',                                   &
        filename=(/trim(stream_fldFileName_popdens)/) , &
        fldListFile='hdm',                             &
        fldListModel='hdm',                            &
        fillalgo='none',                               &
        mapalgo=popdensmapalgo,                        &
        calendar=get_calendar(),                       &
        tintalgo='nearest',                            &
        taxmode='extend'                           )

   if (masterproc) then
      call shr_strdata_print(this%sdat_hdm,'population density data')
   endif

   ! Add history fields
   call hist_addfld1d (fname='HDM', units='counts/km^2',      &
         avgflag='A', long_name='human population density',   &
         ptr_lnd=this%forc_hdm, default='inactive')

  end subroutine hdm_init
  
  !-----------------------------------------------------------------------
  subroutine hdm_interp( this, bounds)
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for population density.
  !
  ! !USES:
  use clm_time_manager, only : get_curr_date
  !
  ! !ARGUMENTS:
  class(cnfire_base_type)       :: this
  type(bounds_type), intent(in) :: bounds  
  !
  ! !LOCAL VARIABLES:
  integer :: g, ig
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(this%sdat_hdm, mcdate, sec, mpicom, 'hdmdyn')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      this%forc_hdm(g) = this%sdat_hdm%avs(1)%rAttr(1,ig)
   end do
   
  end subroutine hdm_interp

  !-----------------------------------------------------------------------
  subroutine lnfm_init( this, bounds, NLFilename )
  !
  ! !DESCRIPTION:
  !
  ! Initialize data stream information for Lightning.
  !
  ! !USES:
  use clm_varctl       , only : inst_name
  use clm_time_manager , only : get_calendar
  use ncdio_pio        , only : pio_subsystem
  use shr_pio_mod      , only : shr_pio_getiotype
  use clm_nlUtilsMod   , only : find_nlgroup_name
  use ndepStreamMod    , only : clm_domain_mct
  use histFileMod      , only : hist_addfld1d
  !
  ! !ARGUMENTS:
  implicit none
  class(cnfire_base_type)       :: this
  type(bounds_type), intent(in) :: bounds  
  character(len=*),  intent(in) :: NLFilename
  !
  ! !LOCAL VARIABLES:
  integer            :: stream_year_first_lightng  ! first year in Lightning stream to use
  integer            :: stream_year_last_lightng   ! last year in Lightning stream to use
  integer            :: model_year_align_lightng   ! align stream_year_first_lnfm with 
  integer            :: nu_nml                     ! unit for namelist file
  integer            :: nml_error                  ! namelist i/o error flag
  type(mct_ggrid)    :: dom_clm                    ! domain information 
  character(len=CL)  :: stream_fldFileName_lightng ! lightning stream filename to read
  character(len=CL)  :: lightngmapalgo = 'bilinear'! Mapping alogrithm
  character(*), parameter :: subName = "('lnfmdyn_init')"
  character(*), parameter :: F00 = "('(lnfmdyn_init) ',4a)"
  !-----------------------------------------------------------------------

   namelist /light_streams/         &
        stream_year_first_lightng,  &
        stream_year_last_lightng,   &
        model_year_align_lightng,   &
        lightngmapalgo,             &
        stream_fldFileName_lightng

   ! Default values for namelist
    stream_year_first_lightng  = 1      ! first year in stream to use
    stream_year_last_lightng   = 1      ! last  year in stream to use
    model_year_align_lightng   = 1      ! align stream_year_first_lnfm with this model year
    stream_fldFileName_lightng = ' '

   ! Read light_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'light_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=light_streams,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading light_streams namelist'//errMsg(__FILE__, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_lightng, mpicom)
   call shr_mpi_bcast(stream_year_last_lightng, mpicom)
   call shr_mpi_bcast(model_year_align_lightng, mpicom)
   call shr_mpi_bcast(stream_fldFileName_lightng, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'light_stream settings:'
      write(iulog,*) '  stream_year_first_lightng  = ',stream_year_first_lightng  
      write(iulog,*) '  stream_year_last_lightng   = ',stream_year_last_lightng   
      write(iulog,*) '  model_year_align_lightng   = ',model_year_align_lightng   
      write(iulog,*) '  stream_fldFileName_lightng = ',stream_fldFileName_lightng
      write(iulog,*) ' '
   endif

   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(this%sdat_lnfm,name="clmlnfm",  &
        pio_subsystem=pio_subsystem,                  & 
        pio_iotype=shr_pio_getiotype(inst_name),      &
        mpicom=mpicom, compid=comp_id,                &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,       &
        nxg=ldomain%ni, nyg=ldomain%nj,               &
        yearFirst=stream_year_first_lightng,          &
        yearLast=stream_year_last_lightng,            &
        yearAlign=model_year_align_lightng,           &
        offset=0,                                     &
        domFilePath='',                               &
        domFileName=trim(stream_fldFileName_lightng), &
        domTvarName='time',                           &
        domXvarName='lon' ,                           &
        domYvarName='lat' ,                           &  
        domAreaName='area',                           &
        domMaskName='mask',                           &
        filePath='',                                  &
        filename=(/trim(stream_fldFileName_lightng)/),&
        fldListFile='lnfm',                           &
        fldListModel='lnfm',                          &
        fillalgo='none',                              &
        mapalgo=lightngmapalgo,                       &
        calendar=get_calendar(),                      &
        taxmode='cycle'                            )

   if (masterproc) then
      call shr_strdata_print(this%sdat_lnfm,'Lightning data')
   endif

   ! Add history fields
   call hist_addfld1d (fname='LNFM', units='counts/km^2/hr',  &
         avgflag='A', long_name='Lightning frequency',        &
         ptr_lnd=this%forc_lnfm, default='inactive')

  end subroutine lnfm_init
  
  !-----------------------------------------------------------------------
  subroutine lnfm_interp(this, bounds )
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for Lightning.
  !
  ! !USES:
  use clm_time_manager, only : get_curr_date
  !
  ! !ARGUMENTS:
  class(cnfire_base_type)       :: this
  type(bounds_type), intent(in) :: bounds  
  !
  ! !LOCAL VARIABLES:
  integer :: g, ig
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(this%sdat_lnfm, mcdate, sec, mpicom, 'lnfmdyn')

   ig = 0
   do g = bounds%begg,bounds%endg
      ig = ig+1
      this%forc_lnfm(g) = this%sdat_lnfm%avs(1)%rAttr(1,ig)
   end do
   
  end subroutine lnfm_interp

end module CNFireBaseMod
