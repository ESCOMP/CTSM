module atm2lndType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd, lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use clm_varpar    , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon    , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl    , only : iulog, use_c13, use_cn, use_lch4, use_cndv, use_fates, use_luna
  use decompMod     , only : bounds_type
  use abortutils    , only : endrun
  use PatchType     , only : patch
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC DATA TYPES:

  type, public :: atm2lnd_params_type
     ! true => repartition rain/snow from atm based on temperature
     logical :: repartition_rain_snow

     ! true => downscale longwave radiation
     logical :: glcmec_downscale_longwave

     ! Surface temperature lapse rate (K m-1)
     real(r8) :: lapse_rate

     ! longwave radiation lapse rate (W m-2 m-1)
     real(r8) :: lapse_rate_longwave

     ! Relative limit for how much longwave downscaling can be done (unitless)
     ! The pre-normalized, downscaled longwave is restricted to be in the range
     ! [lwrad*(1-longwave_downscaling_limit), lwrad*(1+longwave_downscaling_limit)]
     real(r8) :: longwave_downscaling_limit

     ! Rain-snow ramp for glacier landunits
     ! frac_rain = (temp - all_snow_t) * frac_rain_slope
     ! (all_snow_t is in K)
     real(r8) :: precip_repartition_glc_all_snow_t
     real(r8) :: precip_repartition_glc_frac_rain_slope

     ! Rain-snow ramp for non-glacier landunits
     ! frac_rain = (temp - all_snow_t) * frac_rain_slope
     ! (all_snow_t is in K)
     real(r8) :: precip_repartition_nonglc_all_snow_t
     real(r8) :: precip_repartition_nonglc_frac_rain_slope
  end type atm2lnd_params_type

  !----------------------------------------------------
  ! atmosphere -> land variables structure
  !
  ! NOTE:
  ! IF there are forcing variables that are downscaled - then the
  ! non-downscaled versions SHOULD NOT be used in the code. Currently
  ! the non-downscaled versions are only used n a handful of places in
  ! the code (and needs to be used in lnd_import_export and the
  ! downscaling routines), but in general should NOT be used in new
  ! code. Instead use the datatype variables that have a _col suffix
  ! which gives the downscaled versions of these fields.
  !----------------------------------------------------
  type, public :: atm2lnd_type
     type(atm2lnd_params_type) :: params

     ! atm->lnd not downscaled
     real(r8), pointer :: forc_u_grc                    (:)   => null() ! atm wind speed, east direction (m/s)
     real(r8), pointer :: forc_v_grc                    (:)   => null() ! atm wind speed, north direction (m/s)
     real(r8), pointer :: forc_wind_grc                 (:)   => null() ! atmospheric wind speed
     real(r8), pointer :: forc_hgt_grc                  (:)   => null() ! atmospheric reference height (m)
     real(r8), pointer :: forc_topo_grc                 (:)   => null() ! atmospheric surface height (m)
     real(r8), pointer :: forc_hgt_u_grc                (:)   => null() ! obs height of wind [m] (new)
     real(r8), pointer :: forc_hgt_t_grc                (:)   => null() ! obs height of temperature [m] (new)
     real(r8), pointer :: forc_hgt_q_grc                (:)   => null() ! obs height of humidity [m] (new)
     real(r8), pointer :: forc_vp_grc                   (:)   => null() ! atmospheric vapor pressure (Pa)
     real(r8), pointer :: forc_psrf_grc                 (:)   => null() ! surface pressure (Pa)
     real(r8), pointer :: forc_pco2_grc                 (:)   => null() ! CO2 partial pressure (Pa)
     real(r8), pointer :: forc_pco2_240_patch           (:)   => null() ! 10-day mean CO2 partial pressure (Pa)
     real(r8), pointer :: forc_solad_grc                (:,:) => null() ! direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai_grc                (:,:) => null() ! diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
     real(r8), pointer :: forc_solar_grc                (:)   => null() ! incident solar radiation
     real(r8), pointer :: forc_ndep_grc                 (:)   => null() ! nitrogen deposition rate (gN/m2/s)
     real(r8), pointer :: forc_pc13o2_grc               (:)   => null() ! C13O2 partial pressure (Pa)
     real(r8), pointer :: forc_po2_grc                  (:)   => null() ! O2 partial pressure (Pa)
     real(r8), pointer :: forc_po2_240_patch            (:)   => null() ! 10-day mean O2 partial pressure (Pa)
     real(r8), pointer :: forc_aer_grc                  (:,:) => null() ! aerosol deposition array
     real(r8), pointer :: forc_pch4_grc                 (:)   => null() ! CH4 partial pressure (Pa)

     real(r8), pointer :: forc_t_not_downscaled_grc     (:)   => null() ! not downscaled atm temperature (Kelvin)       
     real(r8), pointer :: forc_th_not_downscaled_grc    (:)   => null() ! not downscaled atm potential temperature (Kelvin)    
     real(r8), pointer :: forc_pbot_not_downscaled_grc  (:)   => null() ! not downscaled atm pressure (Pa)   
     real(r8), pointer :: forc_pbot240_downscaled_patch (:)   => null() ! 10-day mean downscaled atm pressure (Pa)           
     real(r8), pointer :: forc_rho_not_downscaled_grc   (:)   => null() ! not downscaled atm density (kg/m**3)                      
     real(r8), pointer :: forc_lwrad_not_downscaled_grc (:)   => null() ! not downscaled atm downwrd IR longwave radiation (W/m**2) 

     ! atm->lnd downscaled
     real(r8), pointer :: forc_t_downscaled_col         (:)   => null() ! downscaled atm temperature (Kelvin)
     real(r8), pointer :: forc_th_downscaled_col        (:)   => null() ! downscaled atm potential temperature (Kelvin)
     real(r8), pointer :: forc_pbot_downscaled_col      (:)   => null() ! downscaled atm pressure (Pa)
     real(r8), pointer :: forc_rho_downscaled_col       (:)   => null() ! downscaled atm density (kg/m**3)
     real(r8), pointer :: forc_lwrad_downscaled_col     (:)   => null() ! downscaled atm downwrd IR longwave radiation (W/m**2)


     ! time averaged quantities
     real(r8) , pointer :: fsd24_patch                  (:)   => null() ! patch 24hr average of direct beam radiation 
     real(r8) , pointer :: fsd240_patch                 (:)   => null() ! patch 240hr average of direct beam radiation 
     real(r8) , pointer :: fsi24_patch                  (:)   => null() ! patch 24hr average of diffuse beam radiation 
     real(r8) , pointer :: fsi240_patch                 (:)   => null() ! patch 240hr average of diffuse beam radiation 
     real(r8) , pointer :: wind24_patch                 (:)   => null() ! patch 24-hour running mean of wind
     real(r8) , pointer :: t_mo_patch                   (:)   => null() ! patch 30-day average temperature (Kelvin)
     real(r8) , pointer :: t_mo_min_patch               (:)   => null() ! patch annual min of t_mo (Kelvin)

   contains

     procedure, public  :: Init
     procedure, public  :: InitForTesting  ! version of Init meant for unit testing
     procedure, private :: ReadNamelist
     procedure, private :: InitAllocate
     procedure, private :: InitHistory  
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars
     procedure, public  :: Restart
     procedure, public  :: Clean

  end type atm2lnd_type

  interface atm2lnd_params_type
     module procedure atm2lnd_params_constructor
  end interface atm2lnd_params_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !----------------------------------------------------

contains

  !-----------------------------------------------------------------------
  function atm2lnd_params_constructor(repartition_rain_snow, glcmec_downscale_longwave, &
       lapse_rate, lapse_rate_longwave, longwave_downscaling_limit, &
       precip_repartition_glc_all_snow_t, precip_repartition_glc_all_rain_t, &
       precip_repartition_nonglc_all_snow_t, precip_repartition_nonglc_all_rain_t) &
       result(params)
    !
    ! !DESCRIPTION:
    ! Creates a new instance of atm2lnd_params_type
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(atm2lnd_params_type) :: params  ! function result
    logical, intent(in) :: repartition_rain_snow
    logical, intent(in) :: glcmec_downscale_longwave

    ! Surface temperature lapse rate (K m-1)
    real(r8), intent(in) :: lapse_rate

    ! Longwave radiation lapse rate (W m-2 m-1)
    ! Must be present if glcmec_downscale_longwave is true; ignored otherwise
    real(r8), intent(in), optional :: lapse_rate_longwave

    ! Relative limit for how much longwave downscaling can be done (unitless)
    ! Must be present if glcmec_downscale_longwave is true; ignored otherwise
    real(r8), intent(in), optional :: longwave_downscaling_limit

    ! End-points of the rain-snow ramp for glacier landunits (degrees C)
    ! Must be present if repartition_rain_snow is true; ignored otherwise
    real(r8), intent(in), optional :: precip_repartition_glc_all_snow_t
    real(r8), intent(in), optional :: precip_repartition_glc_all_rain_t

    ! End-points of the rain-snow ramp for non-glacier landunits (degrees C)
    ! Must be present if repartition_rain_snow is true; ignored otherwise
    real(r8), intent(in), optional :: precip_repartition_nonglc_all_snow_t
    real(r8), intent(in), optional :: precip_repartition_nonglc_all_rain_t
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'atm2lnd_params_constructor'
    !-----------------------------------------------------------------------

    params%repartition_rain_snow = repartition_rain_snow
    params%glcmec_downscale_longwave = glcmec_downscale_longwave

    params%lapse_rate = lapse_rate

    if (glcmec_downscale_longwave) then
       if (.not. present(lapse_rate_longwave)) then
          call endrun(subname // &
               ' ERROR: For glcmec_downscale_longwave true, lapse_rate_longwave must be provided')
       end if
       if (.not. present(longwave_downscaling_limit)) then
          call endrun(subname // &
               ' ERROR: For glcmec_downscale_longwave true, longwave_downscaling_limit must be provided')
       end if

       if (longwave_downscaling_limit < 0._r8 .or. &
            longwave_downscaling_limit > 1._r8) then
          call endrun(subname // &
               ' ERROR: longwave_downscaling_limit must be between 0 and 1')
       end if

       params%lapse_rate_longwave = lapse_rate_longwave
       params%longwave_downscaling_limit = longwave_downscaling_limit
    else
       params%lapse_rate_longwave = nan
       params%longwave_downscaling_limit = nan
    end if

    if (repartition_rain_snow) then

       ! Make sure all of the repartitioning-related parameters are present

       if (.not. present(precip_repartition_glc_all_snow_t)) then
          call endrun(subname // &
               ' ERROR: For repartition_rain_snow true, precip_repartition_glc_all_snow_t must be provided')
       end if
       if (.not. present(precip_repartition_glc_all_rain_t)) then
          call endrun(subname // &
               ' ERROR: For repartition_rain_snow true, precip_repartition_glc_all_rain_t must be provided')
       end if
       if (.not. present(precip_repartition_nonglc_all_snow_t)) then
          call endrun(subname // &
               ' ERROR: For repartition_rain_snow true, precip_repartition_nonglc_all_snow_t must be provided')
       end if
       if (.not. present(precip_repartition_nonglc_all_rain_t)) then
          call endrun(subname // &
               ' ERROR: For repartition_rain_snow true, precip_repartition_nonglc_all_rain_t must be provided')
       end if

       ! Do some other error checking

       if (precip_repartition_glc_all_rain_t <= precip_repartition_glc_all_snow_t) then
          call endrun(subname // &
               ' ERROR: Must have precip_repartition_glc_all_snow_t < precip_repartition_glc_all_rain_t')
       end if

       if (precip_repartition_nonglc_all_rain_t <= precip_repartition_nonglc_all_snow_t) then
          call endrun(subname // &
               ' ERROR: Must have precip_repartition_nonglc_all_snow_t < precip_repartition_nonglc_all_rain_t')
       end if

       ! Convert to the form of the parameters we want for the main code

       call compute_ramp_params( &
            all_snow_t_c = precip_repartition_glc_all_snow_t, &
            all_rain_t_c = precip_repartition_glc_all_rain_t, &
            all_snow_t_k = params%precip_repartition_glc_all_snow_t, &
            frac_rain_slope = params%precip_repartition_glc_frac_rain_slope)

       call compute_ramp_params( &
            all_snow_t_c = precip_repartition_nonglc_all_snow_t, &
            all_rain_t_c = precip_repartition_nonglc_all_rain_t, &
            all_snow_t_k = params%precip_repartition_nonglc_all_snow_t, &
            frac_rain_slope = params%precip_repartition_nonglc_frac_rain_slope)

    else  ! .not. repartition_rain_snow
       params%precip_repartition_glc_all_snow_t = nan
       params%precip_repartition_glc_frac_rain_slope = nan
       params%precip_repartition_nonglc_all_snow_t = nan
       params%precip_repartition_nonglc_frac_rain_slope = nan
    end if

  contains
    subroutine compute_ramp_params(all_snow_t_c, all_rain_t_c, &
         all_snow_t_k, frac_rain_slope)
      real(r8), intent(in)  :: all_snow_t_c  ! Temperature at which precip falls entirely as rain (deg C)
      real(r8), intent(in)  :: all_rain_t_c  ! Temperature at which precip falls entirely as snow (deg C)
      real(r8), intent(out) :: all_snow_t_k  ! Temperature at which precip falls entirely as snow (K)
      real(r8), intent(out) :: frac_rain_slope ! Slope of the frac_rain vs. T relationship

      frac_rain_slope = 1._r8 / (all_rain_t_c - all_snow_t_c)
      all_snow_t_k = all_snow_t_c + tfrz
    end subroutine compute_ramp_params

  end function atm2lnd_params_constructor


  !------------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)

    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in) :: NLFilename ! namelist filename

    call this%InitAllocate(bounds)
    call this%ReadNamelist(NLFilename)
    call this%InitHistory(bounds)
    
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitForTesting(this, bounds, params)
    !
    ! !DESCRIPTION:
    ! Does initialization needed for unit testing. Allows caller to prescribe parameter
    ! values (bypassing the namelist read)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! If params isn't provided, we use default values
    type(atm2lnd_params_type), intent(in), optional :: params
    !
    ! !LOCAL VARIABLES:
    type(atm2lnd_params_type) :: l_params

    character(len=*), parameter :: subname = 'InitForTesting'
    !-----------------------------------------------------------------------

    if (present(params)) then
       l_params = params
    else
       ! Use arbitrary values
       l_params = atm2lnd_params_type( &
            repartition_rain_snow = .false., &
            glcmec_downscale_longwave = .false., &
            lapse_rate = 0.01_r8)
    end if

    call this%InitAllocate(bounds)
    this%params = l_params

  end subroutine InitForTesting


  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the atm2lnd namelist
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    class(atm2lnd_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    ! temporary variables corresponding to the components of atm2lnd_params_type
    logical :: repartition_rain_snow
    logical :: glcmec_downscale_longwave
    real(r8) :: lapse_rate
    real(r8) :: lapse_rate_longwave
    real(r8) :: longwave_downscaling_limit
    real(r8) :: precip_repartition_glc_all_snow_t
    real(r8) :: precip_repartition_glc_all_rain_t
    real(r8) :: precip_repartition_nonglc_all_snow_t
    real(r8) :: precip_repartition_nonglc_all_rain_t

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: nmlname = 'atm2lnd_inparm'

    character(len=*), parameter :: subname = 'ReadNamelist'
    !-----------------------------------------------------------------------

    namelist /atm2lnd_inparm/ repartition_rain_snow, glcmec_downscale_longwave, &
         lapse_rate, lapse_rate_longwave, longwave_downscaling_limit, &
         precip_repartition_glc_all_snow_t, precip_repartition_glc_all_rain_t, &
         precip_repartition_nonglc_all_snow_t, precip_repartition_nonglc_all_rain_t

    ! Initialize namelist variables to defaults
    repartition_rain_snow = .false.
    glcmec_downscale_longwave = .false.
    lapse_rate = nan
    lapse_rate_longwave = nan
    longwave_downscaling_limit = nan
    precip_repartition_glc_all_snow_t = nan
    precip_repartition_glc_all_rain_t = nan
    precip_repartition_nonglc_all_snow_t = nan
    precip_repartition_nonglc_all_rain_t = nan

    if (masterproc) then
       unitn = getavu()
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=atm2lnd_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast(repartition_rain_snow, mpicom)
    call shr_mpi_bcast(glcmec_downscale_longwave, mpicom)
    call shr_mpi_bcast(lapse_rate, mpicom)
    call shr_mpi_bcast(lapse_rate_longwave, mpicom)
    call shr_mpi_bcast(longwave_downscaling_limit, mpicom)
    call shr_mpi_bcast(precip_repartition_glc_all_snow_t, mpicom)
    call shr_mpi_bcast(precip_repartition_glc_all_rain_t, mpicom)
    call shr_mpi_bcast(precip_repartition_nonglc_all_snow_t, mpicom)
    call shr_mpi_bcast(precip_repartition_nonglc_all_rain_t, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       ! Write settings one-by-one rather than with a nml write because some settings may
       ! be NaN if certain options are turned off.
       write(iulog,*) 'repartition_rain_snow = ', repartition_rain_snow
       write(iulog,*) 'glcmec_downscale_longwave = ', glcmec_downscale_longwave
       write(iulog,*) 'lapse_rate = ', lapse_rate
       if (glcmec_downscale_longwave) then
          write(iulog,*) 'lapse_rate_longwave = ', lapse_rate_longwave
          write(iulog,*) 'longwave_downscaling_limit = ', longwave_downscaling_limit
       end if
       if (repartition_rain_snow) then
          write(iulog,*) 'precip_repartition_glc_all_snow_t = ', precip_repartition_glc_all_snow_t
          write(iulog,*) 'precip_repartition_glc_all_rain_t = ', precip_repartition_glc_all_rain_t
          write(iulog,*) 'precip_repartition_nonglc_all_snow_t = ', precip_repartition_nonglc_all_snow_t
          write(iulog,*) 'precip_repartition_nonglc_all_rain_t = ', precip_repartition_nonglc_all_rain_t
       end if
       write(iulog,*) ' '
    end if

    this%params = atm2lnd_params_type( &
         repartition_rain_snow = repartition_rain_snow, &
         glcmec_downscale_longwave = glcmec_downscale_longwave, &
         lapse_rate = lapse_rate, &
         lapse_rate_longwave = lapse_rate_longwave, &
         longwave_downscaling_limit = longwave_downscaling_limit, &
         precip_repartition_glc_all_snow_t = precip_repartition_glc_all_snow_t, &
         precip_repartition_glc_all_rain_t = precip_repartition_glc_all_rain_t, &
         precip_repartition_nonglc_all_snow_t = precip_repartition_nonglc_all_snow_t, &
         precip_repartition_nonglc_all_rain_t = precip_repartition_nonglc_all_rain_t)

  end subroutine ReadNamelist


  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize atm2lnd derived type
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    !------------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    ! atm->lnd
    allocate(this%forc_u_grc                    (begg:endg))        ; this%forc_u_grc                    (:)   = ival
    allocate(this%forc_v_grc                    (begg:endg))        ; this%forc_v_grc                    (:)   = ival
    allocate(this%forc_wind_grc                 (begg:endg))        ; this%forc_wind_grc                 (:)   = ival
    allocate(this%forc_hgt_grc                  (begg:endg))        ; this%forc_hgt_grc                  (:)   = ival
    allocate(this%forc_topo_grc                 (begg:endg))        ; this%forc_topo_grc                 (:)   = ival
    allocate(this%forc_hgt_u_grc                (begg:endg))        ; this%forc_hgt_u_grc                (:)   = ival
    allocate(this%forc_hgt_t_grc                (begg:endg))        ; this%forc_hgt_t_grc                (:)   = ival
    allocate(this%forc_hgt_q_grc                (begg:endg))        ; this%forc_hgt_q_grc                (:)   = ival
    allocate(this%forc_vp_grc                   (begg:endg))        ; this%forc_vp_grc                   (:)   = ival
    allocate(this%forc_psrf_grc                 (begg:endg))        ; this%forc_psrf_grc                 (:)   = ival
    allocate(this%forc_pco2_grc                 (begg:endg))        ; this%forc_pco2_grc                 (:)   = ival
    allocate(this%forc_solad_grc                (begg:endg,numrad)) ; this%forc_solad_grc                (:,:) = ival
    allocate(this%forc_solai_grc                (begg:endg,numrad)) ; this%forc_solai_grc                (:,:) = ival
    allocate(this%forc_solar_grc                (begg:endg))        ; this%forc_solar_grc                (:)   = ival
    allocate(this%forc_ndep_grc                 (begg:endg))        ; this%forc_ndep_grc                 (:)   = ival
    allocate(this%forc_pc13o2_grc               (begg:endg))        ; this%forc_pc13o2_grc               (:)   = ival
    allocate(this%forc_po2_grc                  (begg:endg))        ; this%forc_po2_grc                  (:)   = ival
    allocate(this%forc_aer_grc                  (begg:endg,14))     ; this%forc_aer_grc                  (:,:) = ival
    allocate(this%forc_pch4_grc                 (begg:endg))        ; this%forc_pch4_grc                 (:)   = ival
    if(use_luna)then
     allocate(this%forc_pco2_240_patch          (begp:endp))        ; this%forc_pco2_240_patch           (:)   = ival
     allocate(this%forc_po2_240_patch           (begp:endp))        ; this%forc_po2_240_patch            (:)   = ival
     allocate(this%forc_pbot240_downscaled_patch(begp:endp))        ; this%forc_pbot240_downscaled_patch (:)   = ival
    endif

    ! atm->lnd not downscaled
    allocate(this%forc_t_not_downscaled_grc     (begg:endg))        ; this%forc_t_not_downscaled_grc     (:)   = ival
    allocate(this%forc_pbot_not_downscaled_grc  (begg:endg))        ; this%forc_pbot_not_downscaled_grc  (:)   = ival
    allocate(this%forc_th_not_downscaled_grc    (begg:endg))        ; this%forc_th_not_downscaled_grc    (:)   = ival
    allocate(this%forc_rho_not_downscaled_grc   (begg:endg))        ; this%forc_rho_not_downscaled_grc   (:)   = ival
    allocate(this%forc_lwrad_not_downscaled_grc (begg:endg))        ; this%forc_lwrad_not_downscaled_grc (:)   = ival
    
    ! atm->lnd downscaled
    allocate(this%forc_t_downscaled_col         (begc:endc))        ; this%forc_t_downscaled_col         (:)   = ival
    allocate(this%forc_pbot_downscaled_col      (begc:endc))        ; this%forc_pbot_downscaled_col      (:)   = ival
    allocate(this%forc_th_downscaled_col        (begc:endc))        ; this%forc_th_downscaled_col        (:)   = ival
    allocate(this%forc_rho_downscaled_col       (begc:endc))        ; this%forc_rho_downscaled_col       (:)   = ival
    allocate(this%forc_lwrad_downscaled_col     (begc:endc))        ; this%forc_lwrad_downscaled_col     (:)   = ival

    allocate(this%fsd24_patch                   (begp:endp))        ; this%fsd24_patch                   (:)   = nan
    allocate(this%fsd240_patch                  (begp:endp))        ; this%fsd240_patch                  (:)   = nan
    allocate(this%fsi24_patch                   (begp:endp))        ; this%fsi24_patch                   (:)   = nan
    allocate(this%fsi240_patch                  (begp:endp))        ; this%fsi240_patch                  (:)   = nan
    if (use_fates) then
       allocate(this%wind24_patch               (begp:endp))        ; this%wind24_patch                  (:)   = nan
    end if
    allocate(this%t_mo_patch                    (begp:endp))        ; this%t_mo_patch               (:)   = nan
    allocate(this%t_mo_min_patch                (begp:endp))        ; this%t_mo_min_patch           (:)   = spval ! TODO - initialize this elsewhere

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    !---------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    this%forc_wind_grc(begg:endg) = spval
    call hist_addfld1d (fname='WIND', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_lnd=this%forc_wind_grc)
    ! Rename of WIND for Urban intercomparision project
    call hist_addfld1d (fname='Wind', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_gcell=this%forc_wind_grc, default = 'inactive')

    this%forc_hgt_grc(begg:endg) = spval
    call hist_addfld1d (fname='ZBOT', units='m',  &
         avgflag='A', long_name='atmospheric reference height', &
         ptr_lnd=this%forc_hgt_grc)

    this%forc_topo_grc(begg:endg) = spval
    call hist_addfld1d (fname='ATM_TOPO', units='m', &
         avgflag='A', long_name='atmospheric surface height', &
         ptr_lnd=this%forc_topo_grc)

    this%forc_solar_grc(begg:endg) = spval
    call hist_addfld1d (fname='FSDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_lnd=this%forc_solar_grc)

    this%forc_pco2_grc(begg:endg) = spval
    call hist_addfld1d (fname='PCO2', units='Pa',  &
         avgflag='A', long_name='atmospheric partial pressure of CO2', &
         ptr_lnd=this%forc_pco2_grc)

    this%forc_solar_grc(begg:endg) = spval
    call hist_addfld1d (fname='SWdown', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_gcell=this%forc_solar_grc, default='inactive')

    if (use_lch4) then
       this%forc_pch4_grc(begg:endg) = spval
       call hist_addfld1d (fname='PCH4', units='Pa',  &
            avgflag='A', long_name='atmospheric partial pressure of CH4', &
            ptr_lnd=this%forc_pch4_grc)
    end if

    this%forc_t_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='Tair_from_atm', units='K',  &
         avgflag='A', long_name='atmospheric air temperature received from atmosphere (pre-downscaling)', &
         ptr_gcell=this%forc_t_not_downscaled_grc, default='inactive')

    this%forc_t_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='TBOT', units='K',  &
         avgflag='A', long_name='atmospheric air temperature (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_t_downscaled_col)
    call hist_addfld1d (fname='Tair', units='K', &
         avgflag='A', long_name='atmospheric air temperature (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_t_downscaled_col, default='inactive')

    this%forc_pbot_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='PBOT', units='Pa',  &
         avgflag='A', long_name='atmospheric pressure at surface (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_pbot_downscaled_col)
    call hist_addfld1d (fname='PSurf', units='Pa',  &
         avgflag='A', long_name='atmospheric pressure at surface (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_pbot_downscaled_col, default='inactive')

    this%forc_lwrad_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='FLDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_lwrad_downscaled_col)
    call hist_addfld1d (fname='LWdown', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_lwrad_downscaled_col, default='inactive')

    call hist_addfld1d (fname='FLDS_ICE', units='W/m^2',  &
         avgflag='A', &
         long_name='atmospheric longwave radiation (downscaled to columns in glacier regions) (ice landunits only)', &
         ptr_col=this%forc_lwrad_downscaled_col, l2g_scale_type='ice', &
         default='inactive')

    this%forc_th_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='THBOT', units='K',  &
         avgflag='A', long_name='atmospheric air potential temperature (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_th_downscaled_col)


    ! Time averaged quantities
    this%fsi24_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSI24', units='K',  &
         avgflag='A', long_name='indirect radiation (last 24hrs)', &
         ptr_patch=this%fsi24_patch, default='inactive')

    this%fsi240_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSI240', units='K',  &
         avgflag='A', long_name='indirect radiation (last 240hrs)', &
         ptr_patch=this%fsi240_patch, default='inactive')

    this%fsd24_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSD24', units='K',  &
         avgflag='A', long_name='direct radiation (last 24hrs)', &
         ptr_patch=this%fsd24_patch, default='inactive')

    this%fsd240_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSD240', units='K',  &
         avgflag='A', long_name='direct radiation (last 240hrs)', &
         ptr_patch=this%fsd240_patch, default='inactive')

    if (use_cndv) then
       call hist_addfld1d (fname='TDA', units='K',  &
            avgflag='A', long_name='daily average 2-m temperature', &
            ptr_patch=this%t_mo_patch)
    end if

    if(use_luna)then
       this%forc_pco2_240_patch = spval
       call hist_addfld1d (fname='PCO2_240', units='Pa',  &
            avgflag='A', long_name='10 day running mean of CO2 pressure', &
            ptr_patch=this%forc_pco2_240_patch, default='inactive')
       this%forc_po2_240_patch = spval
      call hist_addfld1d (fname='PO2_240', units='Pa',  &
            avgflag='A', long_name='10 day running mean of O2 pressure', &
            ptr_patch=this%forc_po2_240_patch, default='inactive')
       this%forc_pbot240_downscaled_patch = spval
       call hist_addfld1d (fname='PBOT_240', units='Pa',  &
            avgflag='A', long_name='10 day running mean of air pressure', &
            ptr_patch=this%forc_pbot240_downscaled_patch, default='inactive')
    endif

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use clm_varcon  , only : spval
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !---------------------------------------------------------------------

    this%fsd24_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSD24', units='W/m2',                                             &
         desc='24hr average of direct solar radiation',  accum_type='runmean', accum_period=-1,    &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsd240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSD240', units='W/m2',                                            &
         desc='240hr average of direct solar radiation',  accum_type='runmean', accum_period=-10,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsi24_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSI24', units='W/m2',                                             &
         desc='24hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-1,   &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsi240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSI240', units='W/m2',                                            &
         desc='240hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    if ( use_fates ) then
       call init_accum_field (name='WIND24', units='m', &
            desc='24hr average of wind', accum_type='runmean', accum_period=-1, &
            subgrid_type='pft', numlev=1, init_value=0._r8)
    end if

    if(use_luna) then
      this%forc_po2_240_patch(bounds%begp:bounds%endp) = spval
      call init_accum_field (name='po2_240', units='Pa',                                            &
         desc='10-day running mean of parial O2 pressure',  accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=21223._r8)

      this%forc_pco2_240_patch(bounds%begp:bounds%endp) = spval
      call init_accum_field (name='pco2_240', units='Pa',                                            &
         desc='10-day running mean of parial CO2 pressure',  accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=28._r8)

      this%forc_pbot240_downscaled_patch(bounds%begp:bounds%endp) = spval
      call init_accum_field (name='pbot240', units='Pa',                                            &
         desc='10-day running mean of air pressure',  accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=101325._r8)

    endif

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: begc, endc
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    real(r8), pointer :: rbufslc(:)  ! temporary
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="InitAccVars allocation error for rbufslp"//&
            errMsg(sourcefile, __LINE__))
    endif
    ! Allocate needed dynamic memory for single level col field
    allocate(rbufslc(begc:endc), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="InitAccVars allocation error for rbufslc"//&
            errMsg(sourcefile, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('FSD24', rbufslp, nstep)
    this%fsd24_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSD240', rbufslp, nstep)
    this%fsd240_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSI24', rbufslp, nstep)
    this%fsi24_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSI240', rbufslp, nstep)
    this%fsi240_patch(begp:endp) = rbufslp(begp:endp)

    if (use_cndv) then
       call extract_accum_field ('TDA', rbufslp, nstep) 
       this%t_mo_patch(begp:endp) = rbufslp(begp:endp)
    end if

    if (use_fates) then
       call extract_accum_field ('WIND24', rbufslp, nstep)
       this%wind24_patch(begp:endp) = rbufslp(begp:endp)
    end if

    if(use_luna) then
       call extract_accum_field ('po2_240', rbufslp, nstep)
       this%forc_po2_240_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('pco2_240', rbufslp, nstep)
       this%forc_pco2_240_patch(begp:endp) = rbufslp(begp:endp)   

       call extract_accum_field ('pbot240', rbufslp, nstep)
       this%forc_pbot240_downscaled_patch(begp:endp) = rbufslp(begp:endp)  
 
    endif

    deallocate(rbufslp)
    deallocate(rbufslc)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(atm2lnd_type)                 :: this
    type(bounds_type)      , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p                     ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: rbufslp(:)      ! temporary single level - patch level
    real(r8), pointer :: rbufslc(:)      ! temporary single level - column level
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'UpdateAccVars allocation error for rbufslp'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif
    ! Allocate needed dynamic memory for single level col field
    allocate(rbufslc(begc:endc), stat=ier)
    if (ier/=0) then
       write(iulog,*)'UpdateAccVars allocation error for rbufslc'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    ! Accumulate and extract forc_solad24 & forc_solad240 
    do p = begp,endp
       g = patch%gridcell(p)
       rbufslp(p) = this%forc_solad_grc(g,1)
    end do
    call update_accum_field  ('FSD240', rbufslp               , nstep)
    call extract_accum_field ('FSD240', this%fsd240_patch     , nstep)
    call update_accum_field  ('FSD24' , rbufslp               , nstep)
    call extract_accum_field ('FSD24' , this%fsd24_patch      , nstep)

    ! Accumulate and extract forc_solai24 & forc_solai240 
    do p = begp,endp
       g = patch%gridcell(p)
       rbufslp(p) = this%forc_solai_grc(g,1)
    end do
    call update_accum_field  ('FSI24' , rbufslp               , nstep)
    call extract_accum_field ('FSI24' , this%fsi24_patch      , nstep)
    call update_accum_field  ('FSI240', rbufslp               , nstep)
    call extract_accum_field ('FSI240', this%fsi240_patch     , nstep)


    if (use_cndv) then

       ! Accumulate and extract TDA (accumulates TBOT as 30-day average) and 
       ! also determines t_mo_min
       
       do p = begp,endp
          c = patch%column(p)
          rbufslp(p) = this%forc_t_downscaled_col(c)
       end do
       call update_accum_field  ('TDA', rbufslp, nstep)
       call extract_accum_field ('TDA', rbufslp, nstep)
       do p = begp,endp
          this%t_mo_patch(p) = rbufslp(p)
          this%t_mo_min_patch(p) = min(this%t_mo_min_patch(p), rbufslp(p))
       end do

    end if

    if (use_fates) then
       do p = bounds%begp,bounds%endp
          g = patch%gridcell(p) 
          rbufslp(p) = this%forc_wind_grc(g) 
       end do
       call update_accum_field  ('WIND24', rbufslp, nstep)
       call extract_accum_field ('WIND24', this%wind24_patch, nstep)

    end if

    if(use_luna) then
     do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       rbufslp(p) = this%forc_pco2_grc(g) 
     enddo
     call update_accum_field  ('pco2_240', rbufslp, nstep)
     call extract_accum_field ('pco2_240', this%forc_pco2_240_patch, nstep)

     do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       rbufslp(p) = this%forc_po2_grc(g) 
     enddo
     call update_accum_field  ('po2_240', rbufslp, nstep)
     call extract_accum_field ('po2_240', this%forc_po2_240_patch, nstep)

     do p = bounds%begp,bounds%endp
       c = patch%column(p)
       rbufslp(p) = this%forc_pbot_downscaled_col(c) 
     enddo
     call update_accum_field  ('pbot240', rbufslp, nstep)
     call extract_accum_field ('pbot240', this%forc_pbot240_downscaled_patch, nstep)

    endif

    deallocate(rbufslp)
    deallocate(rbufslc)

  end subroutine UpdateAccVars

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar 
    !------------------------------------------------------------------------

    if (use_cndv) then
       call restartvar(ncid=ncid, flag=flag, varname='T_MO_MIN', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%t_mo_min_patch)
    end if

    if(use_luna)then
       call restartvar(ncid=ncid, flag=flag, varname='pco2_240', xtype=ncd_double,  &
            dim1name='pft', long_name='10-day mean CO2 partial pressure', units='Pa', &
            interpinic_flag='interp', readvar=readvar, data=this%forc_pco2_240_patch )
       call restartvar(ncid=ncid, flag=flag, varname='po2_240', xtype=ncd_double,  &
            dim1name='pft', long_name='10-day mean O2 partial pressure', units='Pa', &
            interpinic_flag='interp', readvar=readvar, data=this%forc_po2_240_patch )
       call restartvar(ncid=ncid, flag=flag, varname='pbot240', xtype=ncd_double,  &
            dim1name='pft', long_name='10 day mean atmospheric pressure(Pa)', units='Pa', &
            interpinic_flag='interp', readvar=readvar, data=this%forc_pbot240_downscaled_patch )
    endif

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Finalize this instance
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(atm2lnd_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    ! atm->lnd
    deallocate(this%forc_u_grc)
    deallocate(this%forc_v_grc)
    deallocate(this%forc_wind_grc)
    deallocate(this%forc_hgt_grc)
    deallocate(this%forc_topo_grc)
    deallocate(this%forc_hgt_u_grc)
    deallocate(this%forc_hgt_t_grc)
    deallocate(this%forc_hgt_q_grc)
    deallocate(this%forc_vp_grc)
    deallocate(this%forc_psrf_grc)
    deallocate(this%forc_pco2_grc)
    deallocate(this%forc_solad_grc)
    deallocate(this%forc_solai_grc)
    deallocate(this%forc_solar_grc)
    deallocate(this%forc_ndep_grc)
    deallocate(this%forc_pc13o2_grc)
    deallocate(this%forc_po2_grc)
    deallocate(this%forc_aer_grc)
    deallocate(this%forc_pch4_grc)

    ! atm->lnd not downscaled
    deallocate(this%forc_t_not_downscaled_grc)
    deallocate(this%forc_pbot_not_downscaled_grc)
    deallocate(this%forc_th_not_downscaled_grc)
    deallocate(this%forc_rho_not_downscaled_grc)
    deallocate(this%forc_lwrad_not_downscaled_grc)
    
    ! atm->lnd downscaled
    deallocate(this%forc_t_downscaled_col)
    deallocate(this%forc_pbot_downscaled_col)
    deallocate(this%forc_th_downscaled_col)
    deallocate(this%forc_rho_downscaled_col)
    deallocate(this%forc_lwrad_downscaled_col)

    deallocate(this%fsd24_patch)
    deallocate(this%fsd240_patch)
    deallocate(this%fsi24_patch)
    deallocate(this%fsi240_patch)
    if (use_fates) then
       deallocate(this%wind24_patch)
    end if
    deallocate(this%t_mo_patch)
    deallocate(this%t_mo_min_patch)

  end subroutine Clean


end module atm2lndType
