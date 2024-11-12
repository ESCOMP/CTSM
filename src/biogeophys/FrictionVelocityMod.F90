module FrictionVelocityMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of the friction velocity, relation for potential
  ! temperature and humidity profiles of surface boundary layer.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use shr_const_mod           , only : SHR_CONST_PI
  use decompMod               , only : bounds_type
  use abortutils              , only : endrun
  use clm_varcon              , only : spval
  use clm_varctl              , only : use_cn, use_luna, z0param_method, use_z0m_snowmelt
  use LandunitType            , only : lun
  use ColumnType              , only : col
  use PatchType               , only : patch
  use landunit_varcon         , only : istsoil, istcrop, istice, istwet
  use ncdio_pio               , only : file_desc_t
  use paramUtilMod            , only : readNcdioScalar
  use atm2lndType             , only : atm2lnd_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use CanopyStateType         , only : canopystate_type
  use WaterFluxBulkType       , only : waterfluxbulk_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: frictionvel_type
     private

     ! Scalar parameters
     real(r8), public :: zetamaxstable = -999._r8  ! Max value zeta ("height" used in Monin-Obukhov theory) can go to under stable conditions
     real(r8) :: zsno = -999._r8  ! Momentum roughness length for snow (m)
     real(r8) :: zlnd = -999._r8  ! Momentum roughness length for soil, glacier, wetland (m)
     real(r8) :: zglc = -999._r8  ! Momentum roughness length for glacier (only used with z0param_method = 'Meier2022') (m)

     ! Roughness length/resistance for friction velocity calculation

     real(r8), pointer, public :: forc_hgt_u_patch (:)   ! patch wind forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: forc_hgt_t_patch (:)   ! patch temperature forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: forc_hgt_q_patch (:)   ! patch specific humidity forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: u10_patch        (:)   ! patch 10-m wind (m/s) (for dust model)
     real(r8), pointer, public :: u10_clm_patch    (:)   ! patch 10-m wind (m/s) (for clm_map2gcell)
     real(r8), pointer, public :: va_patch         (:)   ! patch atmospheric wind speed plus convective velocity (m/s)
     real(r8), pointer, public :: vds_patch        (:)   ! patch deposition velocity term (m/s) (for dry dep SO4, NH4NO3)
     real(r8), pointer, public :: fv_patch         (:)   ! patch friction velocity (m/s) (for dust model)
     real(r8), pointer, public :: rb1_patch        (:)   ! patch aerodynamical resistance (s/m) (for dry deposition of chemical tracers)
     real(r8), pointer, public :: rb10_patch       (:)   ! 10-day mean patch aerodynamical resistance (s/m) (for LUNA model)
     real(r8), pointer, public :: ram1_patch       (:)   ! patch aerodynamical resistance (s/m)
     real(r8), pointer, public :: z0mv_patch       (:)   ! patch roughness length over vegetation, momentum [m]
     real(r8), pointer, public :: z0hv_patch       (:)   ! patch roughness length over vegetation, sensible heat [m]
     real(r8), pointer, public :: z0qv_patch       (:)   ! patch roughness length over vegetation, latent heat [m]
     real(r8), pointer, public :: z0mg_patch       (:)   ! patch roughness length over ground, momentum [m]
     real(r8), pointer, public :: z0hg_patch       (:)   ! patch roughness length over ground, sensible heat [m]
     real(r8), pointer, public :: z0qg_patch       (:)   ! patch roughness length over ground, latent heat [m]
     real(r8), pointer, public :: kbm1_patch       (:)   ! natural logarithm of z0mg_p/z0hg_p [-]
     real(r8), pointer, public :: z0mg_col         (:)   ! col roughness length over ground, momentum  [m]
     real(r8), pointer, public :: z0mg_2D_col      (:)   ! 2-D field of input col roughness length over ground, momentum  [m]
     real(r8), pointer, public :: z0hg_col         (:)   ! col roughness length over ground, sensible heat [m]
     real(r8), pointer, public :: z0qg_col         (:)   ! col roughness length over ground, latent heat [m]
     ! variables to add history output from CanopyFluxesMod
     real(r8), pointer, public :: rah1_patch       (:)   ! patch sensible heat flux resistance [s/m]
     real(r8), pointer, public :: rah2_patch       (:)   ! patch below-canopy sensible heat flux resistance [s/m]
     real(r8), pointer, public :: raw1_patch       (:)   ! patch moisture flux resistance [s/m]
     real(r8), pointer, public :: raw2_patch       (:)   ! patch below-canopy moisture flux resistance [s/m]
     real(r8), pointer, public :: ustar_patch      (:)   ! patch friction velocity [m/s]
     real(r8), pointer, public :: um_patch         (:)   ! patch wind speed including the stablity effect [m/s]
     real(r8), pointer, public :: uaf_patch        (:)   ! patch canopy air speed [m/s]
     real(r8), pointer, public :: taf_patch        (:)   ! patch canopy air temperature [K]
     real(r8), pointer, public :: qaf_patch        (:)   ! patch canopy humidity [kg/kg]
     real(r8), pointer, public :: obu_patch        (:)   ! patch Monin-Obukhov length [m]
     real(r8), pointer, public :: zeta_patch       (:)   ! patch dimensionless stability parameter
     real(r8), pointer, public :: vpd_patch        (:)   ! patch vapor pressure deficit [Pa]
     real(r8), pointer, public :: num_iter_patch   (:)   ! patch number of iterations
     real(r8), pointer, public :: z0m_actual_patch (:)   ! patch roughness length actually used in flux calculations, momentum [m]

   contains

     ! Public procedures
     procedure, public :: Init
     procedure, public :: Restart
     procedure, public :: SetRoughnessLengthsAndForcHeightsNonLake  ! Set roughness lengths and forcing heights for non-lake points
     procedure, public :: SetActualRoughnessLengths ! Set roughness lengths actually used in flux calculations
     procedure, public :: FrictionVelocity       ! Calculate friction velocity
     procedure, public :: MoninObukIni           ! Initialization of the Monin-Obukhov length

     procedure, public  :: InitForTesting        ! version of Init meant for unit testing

     ! Private procedures
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private :: ReadNamelist
     procedure, private :: ReadParams
     procedure, private, nopass :: StabilityFunc1        ! Stability function for rib < 0.
     procedure, private, nopass :: StabilityFunc2        ! Stability function for rib < 0.

  end type frictionvel_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename, params_ncid)

    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in) :: NLFilename  ! file name of namelist file
    type(file_desc_t),intent(inout) :: params_ncid   ! pio netCDF file id

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)
    call this%ReadNamelist(NLFilename)
    call this%ReadParams(params_ncid)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitForTesting(this, bounds)
    ! Initialization for unit testing, hardcodes namelist and parameter file settings
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)
    this%zetamaxstable = 0.5_r8
    this%zsno = 0.00085_r8
    this%zlnd = 0.000775_r8
    this%zglc = 0.00230000005_r8

  end subroutine InitForTesting

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%forc_hgt_u_patch (begp:endp)) ; this%forc_hgt_u_patch (:)   = nan
    allocate(this%forc_hgt_t_patch (begp:endp)) ; this%forc_hgt_t_patch (:)   = nan
    allocate(this%forc_hgt_q_patch (begp:endp)) ; this%forc_hgt_q_patch (:)   = nan
    allocate(this%u10_patch        (begp:endp)) ; this%u10_patch        (:)   = nan
    allocate(this%u10_clm_patch    (begp:endp)) ; this%u10_clm_patch    (:)   = nan
    allocate(this%va_patch         (begp:endp)) ; this%va_patch         (:)   = nan
    allocate(this%vds_patch        (begp:endp)) ; this%vds_patch        (:)   = nan
    allocate(this%fv_patch         (begp:endp)) ; this%fv_patch         (:)   = nan
    allocate(this%rb1_patch        (begp:endp)) ; this%rb1_patch        (:)   = nan
    allocate(this%rb10_patch       (begp:endp)) ; this%rb10_patch       (:)   = spval
    allocate(this%ram1_patch       (begp:endp)) ; this%ram1_patch       (:)   = nan
    allocate(this%z0mv_patch       (begp:endp)) ; this%z0mv_patch       (:)   = nan
    allocate(this%z0hv_patch       (begp:endp)) ; this%z0hv_patch       (:)   = nan
    allocate(this%z0qv_patch       (begp:endp)) ; this%z0qv_patch       (:)   = nan
    allocate(this%z0mg_patch       (begp:endp)) ; this%z0mg_patch       (:)   = nan
    allocate(this%z0hg_patch       (begp:endp)) ; this%z0hg_patch       (:)   = nan
    allocate(this%z0qg_patch       (begp:endp)) ; this%z0qg_patch       (:)   = nan
    allocate(this%kbm1_patch       (begp:endp)) ; this%kbm1_patch       (:)   = nan
    allocate(this%z0mg_col         (begc:endc)) ; this%z0mg_col         (:)   = nan
    allocate(this%z0mg_2D_col      (begc:endc)) ; this%z0mg_2D_col      (:)   = nan
    allocate(this%z0qg_col         (begc:endc)) ; this%z0qg_col         (:)   = nan
    allocate(this%z0hg_col         (begc:endc)) ; this%z0hg_col         (:)   = nan
    allocate(this%rah1_patch       (begp:endp)) ; this%rah1_patch       (:)   = nan
    allocate(this%rah2_patch       (begp:endp)) ; this%rah2_patch       (:)   = nan
    allocate(this%raw1_patch       (begp:endp)) ; this%raw1_patch       (:)   = nan
    allocate(this%raw2_patch       (begp:endp)) ; this%raw2_patch       (:)   = nan
    allocate(this%um_patch         (begp:endp)) ; this%um_patch         (:)   = nan
    allocate(this%uaf_patch        (begp:endp)) ; this%uaf_patch        (:)   = nan
    allocate(this%taf_patch        (begp:endp)) ; this%taf_patch        (:)   = nan
    allocate(this%qaf_patch        (begp:endp)) ; this%qaf_patch        (:)   = nan
    allocate(this%ustar_patch      (begp:endp)) ; this%ustar_patch      (:)   = nan
    allocate(this%obu_patch        (begp:endp)) ; this%obu_patch        (:)   = nan
    allocate(this%zeta_patch       (begp:endp)) ; this%zeta_patch       (:)   = nan
    allocate(this%vpd_patch        (begp:endp)) ; this%vpd_patch        (:)   = nan
    allocate(this%num_iter_patch   (begp:endp)) ; this%num_iter_patch   (:)   = nan
    allocate(this%z0m_actual_patch (begp:endp)) ; this%z0m_actual_patch (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%z0mg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0MG', units='m', &
         avgflag='A', long_name='roughness length over ground, momentum (vegetated landunits only)', &
         ptr_col=this%z0mg_col, default='inactive', l2g_scale_type='veg')

    this%z0hg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0HG', units='m', &
         avgflag='A', long_name='roughness length over ground, sensible heat (vegetated landunits only)', &
         ptr_col=this%z0hg_col, default='inactive', l2g_scale_type='veg')

    this%z0qg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0QG', units='m', &
         avgflag='A', long_name='roughness length over ground, latent heat (vegetated landunits only)', &
         ptr_col=this%z0qg_col, default='inactive', l2g_scale_type='veg')

    this%va_patch(begp:endp) = spval
    call hist_addfld1d (fname='VA', units='m/s', &
         avgflag='A', long_name='atmospheric wind speed plus convective velocity', &
         ptr_patch=this%va_patch, default='inactive')

    this%u10_clm_patch(begp:endp) = spval
    call hist_addfld1d (fname='U10', units='m/s', &
         avgflag='A', long_name='10-m wind', &
         ptr_patch=this%u10_clm_patch)

    call hist_addfld1d (fname='U10_ICE', units='m/s',  &
         avgflag='A', long_name='10-m wind (ice landunits only)', &
         ptr_patch=this%u10_clm_patch, l2g_scale_type='ice', default='inactive')

    this%u10_patch(begp:endp) = spval
    call hist_addfld1d (fname='U10_DUST', units='m/s', &
         avgflag='A', long_name='10-m wind for dust model', &
         ptr_patch=this%u10_patch)

    if (use_cn) then
       this%ram1_patch(begp:endp) = spval
       call hist_addfld1d (fname='RAM1', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%ram1_patch, default='inactive')
    end if

       this%fv_patch(begp:endp) = spval
       call hist_addfld1d (fname='FV', units='m/s', &
            avgflag='A', long_name='friction velocity', &
            ptr_patch=this%fv_patch, default='inactive')

       call hist_addfld1d (fname='RAH1', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%rah1_patch, default='inactive')
       this%rah2_patch(begp:endp) = spval
       call hist_addfld1d (fname='RAH2', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%rah2_patch, default='inactive')
       this%raw1_patch(begp:endp) = spval
       call hist_addfld1d (fname='RAW1', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%raw1_patch, default='inactive')
       this%raw2_patch(begp:endp) = spval
       call hist_addfld1d (fname='RAW2', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%raw2_patch, default='inactive')
       this%ustar_patch(begp:endp) = spval
       call hist_addfld1d (fname='USTAR', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%ustar_patch, default='inactive')
       this%um_patch(begp:endp) = spval
       call hist_addfld1d (fname='UM', units='m/s', &
            avgflag='A', long_name='wind speed plus stability effect', &
            ptr_patch=this%um_patch, default='inactive')
       this%uaf_patch(begp:endp) = spval
       call hist_addfld1d (fname='UAF', units='m/s', &
            avgflag='A', long_name='canopy air speed ', &
            ptr_patch=this%uaf_patch, default='inactive')
       this%taf_patch(begp:endp) = spval
       call hist_addfld1d (fname='TAF', units='K', &
            avgflag='A', long_name='canopy air temperature', &
            ptr_patch=this%taf_patch, default='inactive')
       this%qaf_patch(begp:endp) = spval
       call hist_addfld1d (fname='QAF', units='kg/kg', &
            avgflag='A', long_name='canopy air humidity', &
            ptr_patch=this%qaf_patch, default='inactive')
       this%obu_patch(begp:endp) = spval
       call hist_addfld1d (fname='OBU', units='m', &
            avgflag='A', long_name='Monin-Obukhov length', &
            ptr_patch=this%obu_patch, default='inactive')
       this%zeta_patch(begp:endp) = spval
       call hist_addfld1d (fname='ZETA', units='unitless', &
            avgflag='A', long_name='dimensionless stability parameter', &
            ptr_patch=this%zeta_patch, default='inactive')
       this%vpd_patch(begp:endp) = spval
       call hist_addfld1d (fname='VPD', units='Pa', &
            avgflag='A', long_name='vpd', &
            ptr_patch=this%vpd_patch, default='inactive')
       this%num_iter_patch(begp:endp) = spval
       call hist_addfld1d (fname='num_iter', units='unitless', &
            avgflag='A', long_name='number of iterations', &
            ptr_patch=this%num_iter_patch, default='inactive')
       this%rb1_patch(begp:endp) = spval
       call hist_addfld1d (fname='RB', units='s/m', &
            avgflag='A', long_name='leaf boundary resistance', &
            ptr_patch=this%rb1_patch, default='inactive')

    if (use_cn) then
       this%z0hv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0HV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, sensible heat', &
            ptr_patch=this%z0hv_patch, default='inactive', l2g_scale_type='veg')

       this%z0mv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0MV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, momentum', &
            ptr_patch=this%z0mv_patch, default='inactive', l2g_scale_type='veg')

       this%z0qv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0QV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, latent heat', &
            ptr_patch=this%z0qv_patch, default='inactive', l2g_scale_type='veg')

       this%z0hg_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0HG_P', units='m', &
            avgflag='A', long_name='patch roughness length over ground, sensible heat', &
            ptr_patch=this%z0hg_patch, default='inactive')

       this%z0mg_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0MG_P', units='m', &
            avgflag='A', long_name='patch roughness length over ground, momentum', &
            ptr_patch=this%z0mg_patch, default='inactive')

       this%z0qg_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0QG_P', units='m', &
            avgflag='A', long_name='patch roughness length over ground, latent heat', &
            ptr_patch=this%z0qg_patch, default='inactive')

       this%kbm1_patch(begp:endp) = spval
       call hist_addfld1d (fname='KBM1', units='unitless', &
            avgflag='A', long_name='natural logarithm of Z0MG_P/Z0HG_P', &
            ptr_patch=this%kbm1_patch, default='inactive')
    end if

    if (use_luna) then
       call hist_addfld1d (fname='RB10', units='s/m', &
            avgflag='A', long_name='10 day running mean boundary layer resistance', &
            ptr_patch=this%rb10_patch, default='inactive')
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface albedos to reasonable values
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l                         ! indices
    !-----------------------------------------------------------------------

    ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
    ! since this is not initialized before first call to CNVegStructUpdate,
    ! and it is required to set the upper bound for canopy top height.
    ! Changed 3/21/08, KO: still needed but don't have sufficient information 
    ! to set this properly (e.g., patch-level displacement height and roughness 
    ! length). So leave at 30m.

    if (use_cn) then
       do p = bounds%begp, bounds%endp
          this%forc_hgt_u_patch(p) = 30._r8
       end do
    end if

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%lakpoi(l)) then !lake
          this%z0mg_col(c) = 0.0004_r8
       end if
    end do
   
  end subroutine InitCold

  !------------------------------------------------------------------------------
  subroutine ReadParams( this, params_ncid )
    !
    ! !ARGUMENTS:
    class(frictionvel_type), intent(inout) :: this
    type(file_desc_t),intent(inout) :: params_ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'ReadParams_FrictionVelocity'
    !--------------------------------------------------------------------

    ! Momentum roughness length for snow (m)
    call readNcdioScalar(params_ncid, 'zsno', subname, this%zsno)
    ! Momentum roughness length for soil, glacier, wetland (m)
    call readNcdioScalar(params_ncid, 'zlnd', subname, this%zlnd)

    ! Separated roughness length for glacier if z0param_method == 'Meier2022'
    if (z0param_method == 'Meier2022') then
       call readNcdioScalar(params_ncid, 'zglc', subname, this%zglc)
    end if

  end subroutine ReadParams

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='Z0MG', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground momentum roughness length', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%z0mg_col)

    call restartvar(ncid=ncid, flag=flag, varname='OBU', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Monin-Obukhov length', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%obu_patch)

    if(use_luna)then
       call restartvar(ncid=ncid, flag=flag, varname='rb10', xtype=ncd_double,  &
            dim1name='pft', long_name='10-day mean boundary layer resistance at the pacth', units='s/m', &
            interpinic_flag='interp', readvar=readvar, data=this%rb10_patch)
    endif

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine ReadNamelist( this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for Friction Velocity
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    class(frictionvel_type), intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'FrictionVelocityReadNamelist'
    character(len=*), parameter :: nmlname = 'friction_velocity'
    !-----------------------------------------------------------------------
    real(r8) :: zetamaxstable
    namelist /friction_velocity/ zetamaxstable

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    zetamaxstable = 0.5_r8

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=friction_velocity, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (zetamaxstable, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=friction_velocity)
       write(iulog,*) ' '
    end if

    this%zetamaxstable = zetamaxstable

  end subroutine ReadNamelist

  !-----------------------------------------------------------------------
  subroutine SetRoughnessLengthsAndForcHeightsNonLake(this, bounds, &
       num_nolakec, filter_nolakec, num_nolakep, filter_nolakep, &
       atm2lnd_inst, waterdiagnosticbulk_inst, canopystate_inst, waterfluxbulk_inst)
    !
    ! !DESCRIPTION:
    ! Set roughness lengths and forcing heights for non-lake points
    !
    ! !USES:
    use clm_varcon  , only : rpi, b1_param, b4_param, meier_param1, meier_param2
    ! !ARGUMENTS:
    class(frictionvel_type)        , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds    
    integer                        , intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
    integer                        , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    integer                        , intent(in)    :: num_nolakep       ! number of column non-lake points in patch filter
    integer                        , intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
    type(atm2lnd_type)             , intent(in)    :: atm2lnd_inst
    type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
    type(canopystate_type)         , intent(in)    :: canopystate_inst
    type(waterfluxbulk_type)       , intent(in)    :: waterfluxbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    integer :: fp, p
    integer :: l, g

    character(len=*), parameter :: subname = 'SetRoughnessLengthsAndForcHeightsNonLake'
    !-----------------------------------------------------------------------

    associate( &
         z0mv             =>    this%z0mv_patch                       , & ! Output: [real(r8) (:)   ] roughness length over vegetation, momentum [m]
         z0hv             =>    this%z0hv_patch                       , & ! Output: [real(r8) (:)   ] roughness length over vegetation, sensible heat [m]
         z0qv             =>    this%z0qv_patch                       , & ! Output: [real(r8) (:)   ] roughness length over vegetation, latent heat [m]
         z0mg_p           =>    this%z0mg_patch                       , & ! Output: [real(r8) (:)   ] patch roughness length over ground, momentum [m]
         z0hg_p           =>    this%z0hg_patch                       , & ! Output: [real(r8) (:)   ] patch roughness length over ground, sensible heat [m]
         z0qg_p           =>    this%z0qg_patch                       , & ! Output: [real(r8) (:)   ] patch roughness length over ground, latent heat [m]
         kbm1             =>    this%kbm1_patch                       , & ! Output: [real(r8) (:)   ] natural logarithm of z0mg_p/z0hg_p [-]
         z0hg             =>    this%z0hg_col                         , & ! Output: [real(r8) (:)   ] roughness length over ground, sensible heat [m]
         z0mg             =>    this%z0mg_col                         , & ! Output: [real(r8) (:)   ] roughness length over ground, momentum [m]
         z0qg             =>    this%z0qg_col                         , & ! Output: [real(r8) (:)   ] roughness length over ground, latent heat [m]
         forc_hgt_t_patch =>    this%forc_hgt_t_patch                 , & ! Output: [real(r8) (:)   ] observational height of temperature at patch level [m]
         forc_hgt_q_patch =>    this%forc_hgt_q_patch                 , & ! Output: [real(r8) (:)   ] observational height of specific humidity at patch level [m]
         forc_hgt_u_patch =>    this%forc_hgt_u_patch                 , & ! Output: [real(r8) (:)   ] observational height of wind at patch level [m]
         z0m              =>    canopystate_inst%z0m_patch            , & ! Input: [real(r8) (:)   ] momentum roughness length (m)
         displa           =>    canopystate_inst%displa_patch         , & ! Input: [real(r8) (:)   ] displacement height (m)

         frac_veg_nosno   =>    canopystate_inst%frac_veg_nosno_patch , & ! Input:  [integer  (:)   ] fraction of vegetation not covered by snow (0 OR 1) [-]
         frac_sno         =>    waterdiagnosticbulk_inst%frac_sno_col , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         snomelt_accum    =>    waterdiagnosticbulk_inst%snomelt_accum_col , & ! Input:  [real(r8) (:)   ] accumulated col snow melt for z0m calculation (m H2O)
         urbpoi           =>    lun%urbpoi                            , & ! Input:  [logical  (:)   ] true => landunit is an urban point
         z_0_town         =>    lun%z_0_town                          , & ! Input:  [real(r8) (:)   ] momentum roughness length of urban landunit (m)
         z_d_town         =>    lun%z_d_town                          , & ! Input:  [real(r8) (:)   ] displacement height of urban landunit (m)
         forc_hgt_t       =>    atm2lnd_inst%forc_hgt_t_grc           , & ! Input:  [real(r8) (:)   ] observational height of temperature [m]
         forc_hgt_u       =>    atm2lnd_inst%forc_hgt_u_grc           , & ! Input:  [real(r8) (:)   ] observational height of wind [m]
         forc_hgt_q       =>    atm2lnd_inst%forc_hgt_q_grc           , & ! Input:  [real(r8) (:)   ] observational height of specific humidity [m]
         z0mg_2D          =>    this%z0mg_2D_col                        & ! Input:  [real(r8) (:)   ] 2-D field of input col roughness length over ground, momentum  [m]
         )

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)

       ! Ground roughness lengths over non-lake columns (includes bare ground, ground
       ! underneath canopy, wetlands, etc.)

       select case (z0param_method)
       case ('ZengWang2007')
          if (frac_sno(c) > 0._r8) then
             z0mg(c) = this%zsno
          else
             z0mg(c) = this%zlnd
          end if
       case ('Meier2022')           ! Bare ground and ice have a different value
          l = col%landunit(c)
          if (frac_sno(c) > 0._r8) then ! Do snow first because ice could be snow-covered
             if(use_z0m_snowmelt) then
                if ( snomelt_accum(c) < 1.e-5_r8 )then
                    z0mg(c) = exp(-b1_param * rpi * 0.5_r8 + b4_param) * 1.e-3_r8
                else
                    z0mg(c) = exp(b1_param * (atan((log10(snomelt_accum(c)) + meier_param1) / meier_param2)) + b4_param) * 1.e-3_r8
                end if
             else
                z0mg(c) = this%zsno
             end if
          else if (lun%itype(l) == istice) then
             z0mg(c) = this%zglc
          else
             z0mg(c) = this%zlnd
          end if
       end select

       z0hg(c) = z0mg(c)            ! initial set only
       z0qg(c) = z0mg(c)            ! initial set only


    end do

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)

       ! Roughness lengths over vegetation
       z0mv(p)   = z0m(p)
       z0hv(p)   = z0mv(p)
       z0qv(p)   = z0mv(p)

       ! Set to arbitrary value (will be overwritten by respective modules
       z0mg_p(p)   = spval
       z0hg_p(p)   = spval
       z0qg_p(p)   = spval
       kbm1(p)     = spval

    end do

    ! Make forcing height a patch-level quantity that is the atmospheric forcing 
    ! height plus each patch's z0m+displa
    do fp = 1, num_nolakep
       p = filter_nolakep(fp)
       g = patch%gridcell(p)
       l = patch%landunit(p)
       c = patch%column(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          if (frac_veg_nosno(p) == 0) then
             forc_hgt_u_patch(p) = forc_hgt_u(g) + z0mg(c) + displa(p)
             forc_hgt_t_patch(p) = forc_hgt_t(g) + z0hg(c) + displa(p)
             forc_hgt_q_patch(p) = forc_hgt_q(g) + z0qg(c) + displa(p)
          else
             forc_hgt_u_patch(p) = forc_hgt_u(g) + z0mv(p) + displa(p)
             forc_hgt_t_patch(p) = forc_hgt_t(g) + z0hv(p) + displa(p)
             forc_hgt_q_patch(p) = forc_hgt_q(g) + z0qv(p) + displa(p)
          end if
       else if (lun%itype(l) == istwet .or. lun%itype(l) == istice) then
          forc_hgt_u_patch(p) = forc_hgt_u(g) + z0mg(c) + displa(p)
          forc_hgt_t_patch(p) = forc_hgt_t(g) + z0hg(c) + displa(p)
          forc_hgt_q_patch(p) = forc_hgt_q(g) + z0qg(c) + displa(p)
       else if (urbpoi(l)) then
          forc_hgt_u_patch(p) = forc_hgt_u(g) + z_0_town(l) + z_d_town(l)
          forc_hgt_t_patch(p) = forc_hgt_t(g) + z_0_town(l) + z_d_town(l)
          forc_hgt_q_patch(p) = forc_hgt_q(g) + z_0_town(l) + z_d_town(l)
       end if
    end do

    end associate

  end subroutine SetRoughnessLengthsAndForcHeightsNonLake

  !-----------------------------------------------------------------------
  subroutine SetActualRoughnessLengths(this, bounds, &
       num_exposedvegp, filter_exposedvegp, &
       num_noexposedvegp, filter_noexposedvegp, &
       num_urbanp, filter_urbanp, &
       num_lakep, filter_lakep)
    !
    ! !DESCRIPTION:
    ! Set roughness lengths actually used in flux calculations
    !
    ! !ARGUMENTS:
    class(frictionvel_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_exposedvegp         ! number of points in filter_exposedvegp
    integer                 , intent(in)    :: filter_exposedvegp(:)   ! patch filter for non-snow-covered veg
    integer                 , intent(in)    :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer                 , intent(in)    :: filter_noexposedvegp(:) ! patch filter where frac_veg_nosno is 0 (but does NOT include lake or urban)
    integer                 , intent(in)    :: num_urbanp              ! number of points in filter_urbanp
    integer                 , intent(in)    :: filter_urbanp(:)        ! patch filter for urban
    integer                 , intent(in)    :: num_lakep               ! number of points in filter_lakep
    integer                 , intent(in)    :: filter_lakep(:)         ! patch filter for lake
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p, c, l

    character(len=*), parameter :: subname = 'SetActualRoughnessLengths'
    !-----------------------------------------------------------------------

    associate( &
         z_0_town   => lun%z_0_town          , & ! Input:  [real(r8) (:)] momentum roughness length of urban landunit [m]

         z0mv       => this%z0mv_patch       , & ! Input:  [real(r8) (:)] roughness length over vegetation, momentum [m]
         z0mg       => this%z0mg_col         , & ! Input:  [real(r8) (:)] roughness length over ground, momentum [m]
         z0m_actual => this%z0m_actual_patch   & ! Output: [real(r8) (:)] roughness length actually used in flux calculations, momentum [m]
         )

    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)

       z0m_actual(p) = z0mv(p)
    end do

    do fp = 1, num_noexposedvegp
       p = filter_noexposedvegp(fp)
       c = patch%column(p)

       z0m_actual(p) = z0mg(c)
    end do

    do fp = 1, num_urbanp
       p = filter_urbanp(fp)
       l = patch%landunit(p)

       z0m_actual(p) = z_0_town(l)
    end do

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = patch%column(p)

       z0m_actual(p) = z0mg(c)
    end do

    end associate
  end subroutine SetActualRoughnessLengths

  !------------------------------------------------------------------------------
  subroutine FrictionVelocity(this, lbn, ubn, fn, filtern, &
       displa, z0m, z0h, z0q, &
       obu, iter, ur, um, ustar, &
       temp1, temp2, temp12m, temp22m, fm, landunit_index)
    !
    ! !DESCRIPTION:
    ! Calculation of the friction velocity, relation for potential
    ! temperature and humidity profiles of surface boundary layer.
    ! The scheme is based on the work of Zeng et al. (1998):
    ! Intercomparison of bulk aerodynamic algorithms for the computation
    ! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
    ! Vol. 11, 2628-2644.
    !
    ! !USES:
    use clm_varcon, only : vkc
    use clm_varctl, only : iulog
    !
    ! !ARGUMENTS:
    class(frictionvel_type), intent(inout) :: this
    integer  , intent(in)    :: lbn, ubn                 ! pft/landunit array bounds
    integer  , intent(in)    :: fn                       ! number of filtered pft/landunit elements
    integer  , intent(in)    :: filtern(fn)              ! pft/landunit filter
    real(r8) , intent(in)    :: displa  ( lbn: )         ! displacement height (m) [lbn:ubn]
    real(r8) , intent(in)    :: z0m     ( lbn: )         ! roughness length over vegetation, momentum [m] [lbn:ubn]
    real(r8) , intent(in)    :: z0h     ( lbn: )         ! roughness length over vegetation, sensible heat [m] [lbn:ubn]
    real(r8) , intent(in)    :: z0q     ( lbn: )         ! roughness length over vegetation, latent heat [m] [lbn:ubn]
    real(r8) , intent(in)    :: obu     ( lbn: )         ! monin-obukhov length (m) [lbn:ubn]
    integer  , intent(in)    :: iter                     ! iteration number
    real(r8) , intent(in)    :: ur      ( lbn: )         ! wind speed at reference height [m/s] [lbn:ubn]
    real(r8) , intent(in)    :: um      ( lbn: )         ! wind speed including the stablity effect [m/s] [lbn:ubn]
    real(r8) , intent(out)   :: ustar   ( lbn: )         ! friction velocity [m/s] [lbn:ubn]
    ! temp1, temp12m, temp2, temp22m are "inout" rather than "out" to
    ! prevent returning nan when the code returns from this subroutine
    ! before assigning values to these variables
    real(r8) , intent(inout) :: temp1   ( lbn: )         ! relation for potential temperature profile [lbn:ubn]
    real(r8) , intent(inout) :: temp12m ( lbn: )         ! relation for potential temperature profile applied at 2-m [lbn:ubn]
    real(r8) , intent(inout) :: temp2   ( lbn: )         ! relation for specific humidity profile [lbn:ubn]
    real(r8) , intent(inout) :: temp22m ( lbn: )         ! relation for specific humidity profile applied at 2-m [lbn:ubn]
    real(r8) , intent(inout) :: fm      ( lbn: )         ! diagnose 10m wind (DUST only) [lbn:ubn]
    logical  , intent(in), optional :: landunit_index   ! optional argument that defines landunit or pft level
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: zetam = 1.574_r8 ! transition point of flux-gradient relation (wind profile)
    real(r8), parameter :: zetat = 0.465_r8 ! transition point of flux-gradient relation (temp. profile)
    integer  :: f                           ! pft/landunit filter index
    integer  :: n                           ! pft/landunit index
    integer  :: g                           ! gridcell index
    integer  :: pp                          ! pfti,pftf index
    real(r8) :: zldis(lbn:ubn)              ! reference height "minus" zero displacement heght [m]
    real(r8) :: zeta(lbn:ubn)               ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: tmp1,tmp2,tmp3,tmp4         ! Used to diagnose the 10 meter wind
    real(r8) :: fmnew                       ! Used to diagnose the 10 meter wind
    real(r8) :: fm10                        ! Used to diagnose the 10 meter wind
    real(r8) :: zeta10                      ! Used to diagnose the 10 meter wind
    real(r8) :: vds_tmp                     ! Temporary for dry deposition velocity
    !------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(displa)  == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(z0m)     == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(z0h)     == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(z0q)     == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(obu)     == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(ur)      == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(um)      == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(ustar)   == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(temp1)   == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(temp12m) == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(temp2)   == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(temp22m) == (/ubn/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fm)      == (/ubn/)), sourcefile, __LINE__)

    associate(                                                   & 
         pfti             => lun%patchi                          , & ! Input:  [integer  (:) ] beginning pfti index for landunit         
         pftf             => lun%patchf                          , & ! Input:  [integer  (:) ] final pft index for landunit              
         
         forc_hgt_u_patch => this%forc_hgt_u_patch , & ! Input:  [real(r8) (:) ] observational height of wind at pft level [m]
         forc_hgt_t_patch => this%forc_hgt_t_patch , & ! Input:  [real(r8) (:) ] observational height of temperature at pft level [m]
         forc_hgt_q_patch => this%forc_hgt_q_patch , & ! Input:  [real(r8) (:) ] observational height of specific humidity at pft level [m]
         vds              => this%vds_patch        , & ! Output: [real(r8) (:) ] dry deposition velocity term (m/s) (for SO4 NH4NO3)
         u10              => this%u10_patch        , & ! Output: [real(r8) (:) ] 10-m wind (m/s) (for dust model)        
         u10_clm          => this%u10_clm_patch    , & ! Output: [real(r8) (:) ] 10-m wind (m/s)                         
         va               => this%va_patch         , & ! Output: [real(r8) (:) ] atmospheric wind speed plus convective velocity (m/s)
         fv               => this%fv_patch           & ! Output: [real(r8) (:) ] friction velocity (m/s) (for dust model)
         )

      ! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

      do f = 1, fn
         n = filtern(f)
         if (present(landunit_index)) then
            g = lun%gridcell(n) 
         else
            g = patch%gridcell(n)  
         end if

         ! Wind profile

         if (present(landunit_index)) then
            zldis(n) = forc_hgt_u_patch(pfti(n))-displa(n)
         else
            zldis(n) = forc_hgt_u_patch(n)-displa(n)
         end if
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetam) then
            ustar(n) = vkc*um(n)/(log(-zetam*obu(n)/z0m(n))&
                 - this%StabilityFunc1(-zetam) &
                 + this%StabilityFunc1(z0m(n)/obu(n)) &
                 + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8))
         else if (zeta(n) < 0._r8) then
            ustar(n) = vkc*um(n)/(log(zldis(n)/z0m(n))&
                 - this%StabilityFunc1(zeta(n))&
                 + this%StabilityFunc1(z0m(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            ustar(n) = vkc*um(n)/(log(zldis(n)/z0m(n)) + 5._r8*zeta(n) -5._r8*z0m(n)/obu(n))
         else
            ustar(n) = vkc*um(n)/(log(obu(n)/z0m(n))+5._r8-5._r8*z0m(n)/obu(n) &
                 +(5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if

         if (zeta(n) < 0._r8) then
            vds_tmp = 2.e-3_r8*ustar(n) * ( 1._r8 + (300._r8/(-obu(n)))**0.666_r8)
         else
            vds_tmp = 2.e-3_r8*ustar(n)
         endif

         if (present(landunit_index)) then
            do pp = pfti(n),pftf(n)
               vds(pp) = vds_tmp
            end do
         else
            vds(n) = vds_tmp
         end if

         ! Calculate a 10-m wind (10m + z0m + d)
         ! For now, this will not be the same as the 10-m wind calculated for the dust 
         ! model because the CLM stability functions are used here, not the LSM stability
         ! functions used in the dust model. We will eventually change the dust model to be 
         ! consistent with the following formulation.
         ! Note that the 10-m wind calculated this way could actually be larger than the
         ! atmospheric forcing wind because 1) this includes the convective velocity, 2)
         ! this includes the 1 m/s minimum wind threshold

         ! If forcing height is less than or equal to 10m, then set 10-m wind to um
         if (present(landunit_index)) then
            do pp = pfti(n),pftf(n)
               if (zldis(n)-z0m(n) <= 10._r8) then
                  u10_clm(pp) = um(n)
               else
                  if (zeta(n) < -zetam) then
                     u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(-zetam*obu(n)/(10._r8+z0m(n)))      &
                          - this%StabilityFunc1(-zetam)                              &
                          + this%StabilityFunc1((10._r8+z0m(n))/obu(n))              &
                          + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8)) )
                  else if (zeta(n) < 0._r8) then
                     u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))           &
                          - this%StabilityFunc1(zeta(n))                             &
                          + this%StabilityFunc1((10._r8+z0m(n))/obu(n))) )
                  else if (zeta(n) <=  1._r8) then
                     u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))           &
                          + 5._r8*zeta(n) - 5._r8*(10._r8+z0m(n))/obu(n)) )
                  else
                     u10_clm(pp) = um(n) - ( ustar(n)/vkc*(log(obu(n)/(10._r8+z0m(n)))             &
                          + 5._r8 - 5._r8*(10._r8+z0m(n))/obu(n)                &
                          + (5._r8*log(zeta(n))+zeta(n)-1._r8)) )

                  end if
               end if
               va(pp) = um(n)
            end do
         else
            if (zldis(n)-z0m(n) <= 10._r8) then
               u10_clm(n) = um(n)
            else
               if (zeta(n) < -zetam) then
                  u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(-zetam*obu(n)/(10._r8+z0m(n)))         &
                       - this%StabilityFunc1(-zetam)                                 &
                       + this%StabilityFunc1((10._r8+z0m(n))/obu(n))                 &
                       + 1.14_r8*((-zeta(n))**0.333_r8-(zetam)**0.333_r8)) )
               else if (zeta(n) < 0._r8) then
                  u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))              &
                       - this%StabilityFunc1(zeta(n))                                &
                       + this%StabilityFunc1((10._r8+z0m(n))/obu(n))) )
               else if (zeta(n) <=  1._r8) then
                  u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(zldis(n)/(10._r8+z0m(n)))              &
                       + 5._r8*zeta(n) - 5._r8*(10._r8+z0m(n))/obu(n)) )
               else
                  u10_clm(n) = um(n) - ( ustar(n)/vkc*(log(obu(n)/(10._r8+z0m(n)))                &
                       + 5._r8 - 5._r8*(10._r8+z0m(n))/obu(n)                   &
                       + (5._r8*log(zeta(n))+zeta(n)-1._r8)) )
               end if
            end if
            va(n) = um(n)
         end if

         ! Temperature profile

         if (present(landunit_index)) then
            zldis(n) = forc_hgt_t_patch(pfti(n))-displa(n)
         else
            zldis(n) = forc_hgt_t_patch(n)-displa(n)
         end if
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetat) then
            temp1(n) = vkc/(log(-zetat*obu(n)/z0h(n))&
                 - this%StabilityFunc2(-zetat) &
                 + this%StabilityFunc2(z0h(n)/obu(n)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
         else if (zeta(n) < 0._r8) then
            temp1(n) = vkc/(log(zldis(n)/z0h(n)) &
                 - this%StabilityFunc2(zeta(n)) &
                 + this%StabilityFunc2(z0h(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            temp1(n) = vkc/(log(zldis(n)/z0h(n)) + 5._r8*zeta(n) - 5._r8*z0h(n)/obu(n))
         else
            temp1(n) = vkc/(log(obu(n)/z0h(n)) + 5._r8 - 5._r8*z0h(n)/obu(n) &
                 + (5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if

         ! Humidity profile

         if (present(landunit_index)) then
            if (forc_hgt_q_patch(pfti(n)) == forc_hgt_t_patch(pfti(n)) .and. z0q(n) == z0h(n)) then
               temp2(n) = temp1(n)
            else
               zldis(n) = forc_hgt_q_patch(pfti(n))-displa(n)
               zeta(n) = zldis(n)/obu(n)
               if (zeta(n) < -zetat) then
                  temp2(n) = vkc/(log(-zetat*obu(n)/z0q(n)) &
                       - this%StabilityFunc2(-zetat) &
                       + this%StabilityFunc2(z0q(n)/obu(n)) &
                       + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
               else if (zeta(n) < 0._r8) then
                  temp2(n) = vkc/(log(zldis(n)/z0q(n)) &
                       - this%StabilityFunc2(zeta(n)) &
                       + this%StabilityFunc2(z0q(n)/obu(n)))
               else if (zeta(n) <=  1._r8) then
                  temp2(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
               else
                  temp2(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                       + (5._r8*log(zeta(n))+zeta(n)-1._r8))
               end if
            end if
         else
            if (forc_hgt_q_patch(n) == forc_hgt_t_patch(n) .and. z0q(n) == z0h(n)) then
               temp2(n) = temp1(n)
            else
               zldis(n) = forc_hgt_q_patch(n)-displa(n)
               zeta(n) = zldis(n)/obu(n)
               if (zeta(n) < -zetat) then
                  temp2(n) = vkc/(log(-zetat*obu(n)/z0q(n)) &
                       - this%StabilityFunc2(-zetat) &
                       + this%StabilityFunc2(z0q(n)/obu(n)) &
                       + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
               else if (zeta(n) < 0._r8) then
                  temp2(n) = vkc/(log(zldis(n)/z0q(n)) &
                       - this%StabilityFunc2(zeta(n)) &
                       + this%StabilityFunc2(z0q(n)/obu(n)))
               else if (zeta(n) <=  1._r8) then
                  temp2(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
               else
                  temp2(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                       + (5._r8*log(zeta(n))+zeta(n)-1._r8))
               end if
            endif
         endif

         ! Temperature profile applied at 2-m

         zldis(n) = 2.0_r8 + z0h(n)
         zeta(n) = zldis(n)/obu(n)
         if (zeta(n) < -zetat) then
            temp12m(n) = vkc/(log(-zetat*obu(n)/z0h(n))&
                 - this%StabilityFunc2(-zetat) &
                 + this%StabilityFunc2(z0h(n)/obu(n)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
         else if (zeta(n) < 0._r8) then
            temp12m(n) = vkc/(log(zldis(n)/z0h(n)) &
                 - this%StabilityFunc2(zeta(n))  &
                 + this%StabilityFunc2(z0h(n)/obu(n)))
         else if (zeta(n) <=  1._r8) then
            temp12m(n) = vkc/(log(zldis(n)/z0h(n)) + 5._r8*zeta(n) - 5._r8*z0h(n)/obu(n))
         else
            temp12m(n) = vkc/(log(obu(n)/z0h(n)) + 5._r8 - 5._r8*z0h(n)/obu(n) &
                 + (5._r8*log(zeta(n))+zeta(n)-1._r8))
         end if

         ! Humidity profile applied at 2-m

         if (z0q(n) == z0h(n)) then
            temp22m(n) = temp12m(n)
         else
            zldis(n) = 2.0_r8 + z0q(n)
            zeta(n) = zldis(n)/obu(n)
            if (zeta(n) < -zetat) then
               temp22m(n) = vkc/(log(-zetat*obu(n)/z0q(n)) - &
                    this%StabilityFunc2(-zetat) + this%StabilityFunc2(z0q(n)/obu(n)) &
                    + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(n))**(-0.333_r8)))
            else if (zeta(n) < 0._r8) then
               temp22m(n) = vkc/(log(zldis(n)/z0q(n)) - &
                    this%StabilityFunc2(zeta(n))+this%StabilityFunc2(z0q(n)/obu(n)))
            else if (zeta(n) <=  1._r8) then
               temp22m(n) = vkc/(log(zldis(n)/z0q(n)) + 5._r8*zeta(n)-5._r8*z0q(n)/obu(n))
            else
               temp22m(n) = vkc/(log(obu(n)/z0q(n)) + 5._r8 - 5._r8*z0q(n)/obu(n) &
                    + (5._r8*log(zeta(n))+zeta(n)-1._r8))
            end if
         end if

         ! diagnose 10-m wind for dust model (dstmbl.F)
         ! Notes from C. Zender's dst.F:
         ! According to Bon96 p. 62, the displacement height d (here displa) is
         ! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
         ! Therefore d <= 0.034*z1 and may safely be neglected.
         ! Code from LSM routine SurfaceTemperature was used to obtain u10

         if (present(landunit_index)) then
            zldis(n) = forc_hgt_u_patch(pfti(n))-displa(n)
         else
            zldis(n) = forc_hgt_u_patch(n)-displa(n)
         end if
         zeta(n) = zldis(n)/obu(n)
         if (min(zeta(n), 1._r8) < 0._r8) then
            tmp1 = (1._r8 - 16._r8*min(zeta(n),1._r8))**0.25_r8
            tmp2 = log((1._r8+tmp1*tmp1)/2._r8)
            tmp3 = log((1._r8+tmp1)/2._r8)
            fmnew = 2._r8*tmp3 + tmp2 - 2._r8*atan(tmp1) + 1.5707963_r8
         else
            fmnew = -5._r8*min(zeta(n),1._r8)
         endif
         if (iter == 1) then
            fm(n) = fmnew
         else
            fm(n) = 0.5_r8 * (fm(n)+fmnew)
         end if
         zeta10 = min(10._r8/obu(n), 1._r8)
         if (zeta(n) == 0._r8) zeta10 = 0._r8
         if (zeta10 < 0._r8) then
            tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
            tmp2 = log((1.0_r8 + tmp1*tmp1)/2.0_r8)
            tmp3 = log((1.0_r8 + tmp1)/2.0_r8)
            fm10 = 2.0_r8*tmp3 + tmp2 - 2.0_r8*atan(tmp1) + 1.5707963_r8
         else                ! not stable
            fm10 = -5.0_r8 * zeta10
         end if
         if (present(landunit_index)) then
            tmp4 = log( max( 1.0_r8, forc_hgt_u_patch(pfti(n)) / 10._r8) )
         else 
            tmp4 = log( max( 1.0_r8, forc_hgt_u_patch(n) / 10._r8) )
         end if
         if (present(landunit_index)) then
            do pp = pfti(n),pftf(n)
               u10(pp) = ur(n) - ustar(n)/vkc * (tmp4 - fm(n) + fm10)
               fv(pp)  = ustar(n)
            end do
         else
            u10(n) = ur(n) - ustar(n)/vkc * (tmp4 - fm(n) + fm10)
            fv(n)  = ustar(n)
         end if

      end do

    end associate
  end subroutine FrictionVelocity

  !------------------------------------------------------------------------------
  real(r8) function StabilityFunc1(zeta)
    !
    ! !DESCRIPTION:
    ! Stability function for rib < 0.
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    !
    ! !LOCAL VARIABLES:
    real(r8) :: chik, chik2
    !------------------------------------------------------------------------------

    chik2 = sqrt(1._r8-16._r8*zeta)
    chik = sqrt(chik2)
    StabilityFunc1 = 2._r8*log((1._r8+chik)*0.5_r8) &
         + log((1._r8+chik2)*0.5_r8)-2._r8*atan(chik)+SHR_CONST_PI*0.5_r8

  end function StabilityFunc1

  !------------------------------------------------------------------------------
  real(r8) function StabilityFunc2(zeta)
    !
    ! !DESCRIPTION:
    ! Stability function for rib < 0.
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    !
    ! !LOCAL VARIABLES:
    real(r8) :: chik2
    !------------------------------------------------------------------------------

    chik2 = sqrt(1._r8-16._r8*zeta)
    StabilityFunc2 = 2._r8*log((1._r8+chik2)*0.5_r8)

  end function StabilityFunc2

  !-----------------------------------------------------------------------
  subroutine MoninObukIni (this, ur, thv, dthv, zldis, z0m, um, obu)
    !
    ! !DESCRIPTION:
    ! Initialization of the Monin-Obukhov length.
    ! The scheme is based on the work of Zeng et al. (1998):
    ! Intercomparison of bulk aerodynamic algorithms for the computation
    ! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
    ! Vol. 11, 2628-2644.
    !
    ! !USES:
    use clm_varcon, only : grav
    !
    ! !ARGUMENTS:
    class(frictionvel_type), intent(in) :: this
    real(r8), intent(in)  :: ur    ! wind speed at reference height [m/s]
    real(r8), intent(in)  :: thv   ! virtual potential temperature (kelvin)
    real(r8), intent(in)  :: dthv  ! diff of vir. poten. temp. between ref. height and surface
    real(r8), intent(in)  :: zldis ! reference height "minus" zero displacement heght [m]
    real(r8), intent(in)  :: z0m   ! roughness length, momentum [m]
    real(r8), intent(out) :: um    ! wind speed including the stability effect [m/s]
    real(r8), intent(out) :: obu   ! monin-obukhov length (m)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wc    ! convective velocity [m/s]
    real(r8) :: rib   ! bulk Richardson number
    real(r8) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: ustar ! friction velocity [m/s]
    !-----------------------------------------------------------------------

    ! Initial values of u* and convective velocity

    ustar=0.06_r8
    wc=0.5_r8
    if (dthv >= 0._r8) then
       um=max(ur,0.1_r8)
    else
       um=sqrt(ur*ur+wc*wc)
    endif

    rib=grav*zldis*dthv/(thv*um*um)

    if (rib >= 0._r8) then      ! neutral or stable
       zeta = rib*log(zldis/z0m)/(1._r8-5._r8*min(rib,0.19_r8))
       zeta = min(this%zetamaxstable,max(zeta,0.01_r8 ))
    else                     ! unstable
       zeta=rib*log(zldis/z0m)
       zeta = max(-100._r8,min(zeta,-0.01_r8 ))
    endif

    obu=zldis/zeta

  end subroutine MoninObukIni

end module FrictionVelocityMod
