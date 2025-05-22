module CanopyFluxesMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs calculation of leaf temperature and surface fluxes.
  ! SoilFluxes then determines soil/snow and ground temperatures and updates the surface 
  ! fluxes for the new ground temperature.
  !
  ! !USES:
  use shr_sys_mod           , only : shr_sys_flush
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use abortutils            , only : endrun
  use clm_varctl            , only : iulog, use_cn, use_lch4, use_c13, use_c14, use_cndv, use_fates, &
                                     use_luna, use_hydrstress, use_biomass_heat_storage, z0param_method
  use clm_varpar            , only : nlevgrnd, nlevsno, nlevcan, mxpft
  use pftconMod             , only : pftcon
  use decompMod             , only : bounds_type, subgrid_level_patch
  use ActiveLayerMod        , only : active_layer_type
  use PhotosynthesisMod     , only : Photosynthesis, PhotoSynthesisHydraulicStress, PhotosynthesisTotal, Fractionation
  use EDAccumulateFluxesMod , only : AccumulateFluxes_ED
  use SoilMoistStressMod    , only : calc_effective_soilporosity, calc_volumetric_h2oliq
  use SoilMoistStressMod    , only : calc_root_moist_stress, set_perchroot_opt
  use SimpleMathMod         , only : array_div_vector
  use SurfaceResistanceMod  , only : do_soilevap_beta,do_soil_resistance_sl14
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type
  use EnergyFluxType        , only : energyflux_type
  use FrictionvelocityMod   , only : frictionvel_type
  use OzoneBaseMod          , only : ozone_base_type
  use SoilStateType         , only : soilstate_type
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use TemperatureType       , only : temperature_type
  use WaterFluxBulkType         , only : waterfluxbulk_type
  use WaterStateBulkType        , only : waterstatebulk_type
  use WaterDiagnosticBulkType        , only : waterdiagnosticbulk_type
  use Wateratm2lndBulkType        , only : wateratm2lndbulk_type
  use HumanIndexMod         , only : humanindex_type
  use ch4Mod                , only : ch4_type
  use PhotosynthesisMod     , only : photosyns_type
  use GridcellType          , only : grc                
  use ColumnType            , only : col                
  use PatchType             , only : patch                
  use EDTypesMod            , only : ed_site_type
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  use LunaMod               , only : Update_Photosynthesis_Capacity, Acc24_Climate_LUNA,Acc240_Climate_LUNA,Clear24_Climate_LUNA
  use NumericsMod           , only : truncate_small_values
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyFluxesReadNML     ! Read in namelist settings
  public :: CanopyFluxes            ! Calculate canopy fluxes
  public :: readParams

  type, private :: params_type
     real(r8) :: lai_dl   ! Plant litter area index (m2/m2)
     real(r8) :: z_dl     ! Litter layer thickness (m)
     real(r8) :: a_coef   ! Drag coefficient under less dense canopy (unitless)
     real(r8) :: a_exp    ! Drag exponent under less dense canopy (unitless)
     real(r8) :: csoilc   ! Soil drag coefficient under dense canopy (unitless)
     real(r8) :: cv       ! Turbulent transfer coeff. between canopy surface and canopy air (m/s^(1/2))
     real(r8) :: wind_min ! Minimum wind speed at the atmospheric forcing height (m/s)
  end type params_type
  type(params_type), private ::  params_inst

  !
  ! !PUBLIC DATA MEMBERS:
  ! true => btran is based only on unfrozen soil levels
  logical,  public :: perchroot     = .false.  

  ! true  => btran is based on active layer (defined over two years); 
  ! false => btran is based on currently unfrozen levels
  logical,  public :: perchroot_alt = .false.  
  !
  ! !PRIVATE DATA MEMBERS:
  logical, private :: use_undercanopy_stability = .false.      ! use undercanopy stability term or not
  integer, private :: itmax_canopy_fluxes = -1  ! max # of iterations used in subroutine CanopyFluxes

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine CanopyFluxesReadNML(NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for Canopy Fluxes
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CanopyFluxesReadNML'
    character(len=*), parameter :: nmlname = 'canopyfluxes_inparm'
    !-----------------------------------------------------------------------

    namelist /canopyfluxes_inparm/ use_undercanopy_stability
    namelist /canopyfluxes_inparm/ use_biomass_heat_storage
    namelist /canopyfluxes_inparm/ itmax_canopy_fluxes


    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=canopyfluxes_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if

       if (itmax_canopy_fluxes < 1) then
          call endrun(msg=' ERROR: expecting itmax_canopy_fluxes > 0 ' // &
            errMsg(sourcefile, __LINE__))
       end if

       call relavu( unitn )
    end if

    call shr_mpi_bcast (use_undercanopy_stability, mpicom)
    call shr_mpi_bcast (use_biomass_heat_storage, mpicom)
    call shr_mpi_bcast (itmax_canopy_fluxes, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=canopyfluxes_inparm)
       write(iulog,*) ' '
    end if

  end subroutine CanopyFluxesReadNML

  !------------------------------------------------------------------------------
  subroutine readParams( ncid )
    !
    ! !USES:
    use ncdio_pio, only: file_desc_t
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'readParams_CanopyFluxes'
    !--------------------------------------------------------------------

    !added by K.Sakaguchi for litter resistance: Plant litter area index (m2/m2)
    call readNcdioScalar(ncid, 'lai_dl', subname, params_inst%lai_dl)
    !added by K.Sakaguchi for litter resistance: Litter layer thickness (m)
    call readNcdioScalar(ncid, 'z_dl', subname, params_inst%z_dl)
    ! Drag coefficient under less dense canopy (unitless)
    call readNcdioScalar(ncid, 'a_coef', subname, params_inst%a_coef)
    ! Drag exponent under less dense canopy (unitless)
    call readNcdioScalar(ncid, 'a_exp', subname, params_inst%a_exp)
    ! Drag coefficient for soil under dense canopy (unitless)
    call readNcdioScalar(ncid, 'csoilc', subname, params_inst%csoilc)
    ! Turbulent transfer coeff between canopy surface and canopy air (m/s^(1/2))
    call readNcdioScalar(ncid, 'cv', subname, params_inst%cv)
    ! Minimum wind speed at the atmospheric forcing height (m/s)
    call readNcdioScalar(ncid, 'wind_min', subname, params_inst%wind_min)

  end subroutine readParams

  !------------------------------------------------------------------------------
  subroutine CanopyFluxes(bounds,  num_exposedvegp, filter_exposedvegp,                  &
       clm_fates, nc, active_layer_inst, atm2lnd_inst, canopystate_inst,                 &
       energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst,   &
       temperature_inst, waterfluxbulk_inst, waterstatebulk_inst,                        &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, ch4_inst, ozone_inst,            &
       photosyns_inst, &
       humanindex_inst, soil_water_retention_curve, &
       downreg_patch, leafn_patch, froot_carbon, croot_carbon)
    !
    ! !DESCRIPTION:
    ! 1. Calculates the leaf temperature:
    ! 2. Calculates the leaf fluxes, transpiration, photosynthesis and
    !    updates the dew accumulation due to evaporation.
    !
    ! Method:
    ! Use the Newton-Raphson iteration to solve for the foliage
    ! temperature that balances the surface energy budget:
    !
    ! f(t_veg) = Net radiation - Sensible - Latent = 0
    ! f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
    !
    ! Note:
    ! (1) In solving for t_veg, t_grnd is given from the previous timestep.
    ! (2) The partial derivatives of aerodynamical resistances, which cannot
    !     be determined analytically, are ignored for d(H)/dT and d(LE)/dT
    ! (3) The weighted stomatal resistance of sunlit and shaded foliage is used
    ! (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
    !                                                          => Ec + Eg = Ea
    ! (5) Energy loss is due to: numerical truncation of energy budget equation
    !     (*); and "ecidif" (see the code) which is dropped into the sensible
    !     heat
    ! (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n)
    !     and del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference
    !     of water flux from the leaf between the iteration step (n+1) and (n)
    !     less than 0.1 W/m2; or the iterative steps over 40.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_RGAS, shr_const_pi
    use clm_time_manager   , only : get_step_size_real, get_prev_date, is_near_local_noon
    use clm_varcon         , only : sb, cpair, hvap, vkc, grav, denice, c_to_b
    use clm_varcon         , only : denh2o, tfrz, tlsai_crit, alpha_aero
    use clm_varcon         , only : c14ratio, spval
    use clm_varcon         , only : c_water, c_dry_biomass, c_to_b
    use clm_varcon         , only : nu_param, cd1_param
    use perf_mod           , only : t_startf, t_stopf
    use QSatMod            , only : QSat
    use CLMFatesInterfaceMod, only : hlm_fates_interface_type
    use HumanIndexMod      , only : all_human_stress_indices, fast_human_stress_indices, &
                                    Wet_Bulb, Wet_BulbS, HeatIndex, AppTemp, &
                                    swbgt, hmdex, dis_coi, dis_coiS, THIndex, &
                                    SwampCoolEff, KtoC, VaporPres
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use LunaMod            , only : is_time_to_run_LUNA

    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)            :: bounds 
    integer                                , intent(in)            :: num_exposedvegp        ! number of points in filter_exposedvegp
    integer                                , intent(in)            :: filter_exposedvegp(:)  ! patch filter for non-snow-covered veg
    type(hlm_fates_interface_type)         , intent(inout)         :: clm_fates
    integer                                , intent(in)            :: nc ! clump index
    type(active_layer_type)                , intent(in)            :: active_layer_inst
    type(atm2lnd_type)                     , intent(in)            :: atm2lnd_inst
    type(canopystate_type)                 , intent(inout)         :: canopystate_inst
    type(energyflux_type)                  , intent(inout)         :: energyflux_inst
    type(frictionvel_type)                 , intent(inout)         :: frictionvel_inst
    type(solarabs_type)                    , intent(inout)         :: solarabs_inst
    type(surfalb_type)                     , intent(in)            :: surfalb_inst
    type(soilstate_type)                   , intent(inout)         :: soilstate_inst
    type(temperature_type)                 , intent(inout)         :: temperature_inst
    type(waterstatebulk_type)              , intent(inout)         :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)         , intent(inout)         :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)               , intent(inout)         :: waterfluxbulk_inst
    type(wateratm2lndbulk_type)            , intent(inout)         :: wateratm2lndbulk_inst
    type(ch4_type)                         , intent(inout)         :: ch4_inst
    class(ozone_base_type)                 , intent(inout)         :: ozone_inst
    type(photosyns_type)                   , intent(inout)         :: photosyns_inst
    type(humanindex_type)                  , intent(inout)         :: humanindex_inst
    class(soil_water_retention_curve_type) , intent(in)            :: soil_water_retention_curve
    real(r8), intent(in) :: downreg_patch(bounds%begp:) ! fractional reduction in GPP due to N limitation (dimensionless)
    real(r8), intent(in) :: leafn_patch(bounds%begp:)   ! leaf N (gN/m2)
    real(r8), intent(inout) :: froot_carbon(bounds%begp:)  ! fine root biomass (gC/m2)
    real(r8), intent(inout) :: croot_carbon(bounds%begp:)  ! live coarse root biomass (gC/m2)
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer   :: bsun(:)          ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8), pointer   :: bsha(:)          ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8), parameter :: btran0 = 0.0_r8  ! initial value
    real(r8), parameter :: zii = 1000.0_r8  ! convective boundary layer height [m]
    real(r8), parameter :: beta = 1.0_r8    ! coefficient of conective velocity [-]
    real(r8), parameter :: delmax = 1.0_r8  ! maxchange in  leaf temperature [K]
    real(r8), parameter :: dlemin = 0.1_r8  ! max limit for energy flux convergence [w/m2]
    real(r8), parameter :: dtmin = 0.01_r8  ! max limit for temperature convergence [K]
    integer , parameter :: itmin = 2        ! minimum number of iteration [-]

    !added by K.Sakaguchi for stability formulation
    real(r8), parameter :: ria  = 0.5_r8             ! free parameter for stable formulation (currently = 0.5, "gamma" in Sakaguchi&Zeng,2008)
    real(r8) :: dtime                                ! land model time step (sec)
    real(r8) :: zldis(bounds%begp:bounds%endp)       ! reference height "minus" zero displacement height [m]
    real(r8) :: wc                                   ! convective velocity [m/s]
    real(r8) :: dth(bounds%begp:bounds%endp)         ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv(bounds%begp:bounds%endp)        ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(bounds%begp:bounds%endp)         ! diff of humidity between ref. height and surface
    real(r8) :: ur(bounds%begp:bounds%endp)          ! wind speed at reference height [m/s]
    real(r8) :: temp1(bounds%begp:bounds%endp)       ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begp:bounds%endp)     ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begp:bounds%endp)       ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begp:bounds%endp)     ! relation for specific humidity profile applied at 2-m
    real(r8) :: tstar                                ! temperature scaling parameter
    real(r8) :: qstar                                ! moisture scaling parameter
    real(r8) :: thvstar                              ! virtual potential temperature scaling parameter
    real(r8) :: rpp                                  ! fraction of potential evaporation from leaf [-]
    real(r8) :: rppdry                               ! fraction of potential evaporation through transp [-]
    real(r8) :: cf                                   ! heat transfer coefficient from leaves [-]
    real(r8) :: rb(bounds%begp:bounds%endp)          ! leaf boundary layer resistance [s/m]
    real(r8) :: rah(bounds%begp:bounds%endp,2)       ! thermal resistance [s/m]  (air, ground)
    real(r8) :: raw(bounds%begp:bounds%endp,2)       ! moisture resistance [s/m] (air, ground)
    real(r8) :: wta                                  ! heat conductance for air [m/s]
    real(r8) :: wtg(bounds%begp:bounds%endp)         ! heat conductance for ground [m/s]
    real(r8) :: wtl                                  ! heat conductance for leaf [m/s]
    real(r8) :: wtstem                               ! heat conductance for stem [m/s]
    real(r8) :: wta0(bounds%begp:bounds%endp)        ! normalized heat conductance for air [-]
    real(r8) :: wtl0(bounds%begp:bounds%endp)        ! normalized heat conductance for leaf [-]
    real(r8) :: wtg0                                 ! normalized heat conductance for ground [-]
    real(r8) :: wtstem0(bounds%begp:bounds%endp)     ! normalized heat conductance for stem [-]
    real(r8) :: wtal(bounds%begp:bounds%endp)        ! normalized heat conductance for air and leaf [-]
    real(r8) :: wtga(bounds%begp:bounds%endp)        ! normalized heat cond. for air and ground  [-]
    real(r8) :: wtaq                                 ! latent heat conductance for air [m/s]
    real(r8) :: wtlq                                 ! latent heat conductance for leaf [m/s]
    real(r8) :: wtgq(bounds%begp:bounds%endp)        ! latent heat conductance for ground [m/s]
    real(r8) :: wtaq0(bounds%begp:bounds%endp)       ! normalized latent heat conductance for air [-]
    real(r8) :: wtlq0(bounds%begp:bounds%endp)       ! normalized latent heat conductance for leaf [-]
    real(r8) :: wtgq0                                ! normalized heat conductance for ground [-]
    real(r8) :: wtalq(bounds%begp:bounds%endp)       ! normalized latent heat cond. for air and leaf [-]
    real(r8) :: wtgaq                                ! normalized latent heat cond. for air and ground [-]
    real(r8) :: el(bounds%begp:bounds%endp)          ! vapor pressure on leaf surface [pa]
    real(r8) :: qsatl(bounds%begp:bounds%endp)       ! leaf specific humidity [kg/kg]
    real(r8) :: qsatldT(bounds%begp:bounds%endp)     ! derivative of "qsatl" on "t_veg"
    real(r8) :: e_ref2m                              ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: qsat_ref2m                           ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: gs                                   ! canopy conductance for iwue cal [molH2O/m2ground/s]
    real(r8) :: air(bounds%begp:bounds%endp)         ! atmos. radiation temporay set
    real(r8) :: bir(bounds%begp:bounds%endp)         ! atmos. radiation temporay set
    real(r8) :: cir(bounds%begp:bounds%endp)         ! atmos. radiation temporay set
    real(r8) :: dc1,dc2                              ! derivative of energy flux [W/m2/K]
    real(r8) :: delt                                 ! temporary
    real(r8) :: delq(bounds%begp:bounds%endp)        ! temporary
    real(r8) :: del(bounds%begp:bounds%endp)         ! absolute change in leaf temp in current iteration [K]
    real(r8) :: del2(bounds%begp:bounds%endp)        ! change in leaf temperature in previous iteration [K]
    real(r8) :: dele(bounds%begp:bounds%endp)        ! change in latent heat flux from leaf [K]
    real(r8) :: dels                                 ! change in leaf temperature in current iteration [K]
    real(r8) :: det(bounds%begp:bounds%endp)         ! maximum leaf temp. change in two consecutive iter [K]
    real(r8) :: efeb(bounds%begp:bounds%endp)        ! latent heat flux from leaf (previous iter) [mm/s]
    real(r8) :: efeold                               ! latent heat flux from leaf (previous iter) [mm/s]
    real(r8) :: efpot                                ! potential latent energy flux [kg/m2/s]
    real(r8) :: efe(bounds%begp:bounds%endp)         ! water flux from leaf [mm/s]
    real(r8) :: efsh                                 ! sensible heat from leaf [mm/s]
    real(r8) :: obuold(bounds%begp:bounds%endp)      ! monin-obukhov length from previous iteration
    real(r8) :: tlbef(bounds%begp:bounds%endp)       ! leaf temperature from previous iteration [K]
    real(r8) :: tl_ini(bounds%begp:bounds%endp)      ! leaf temperature from beginning of time step [K]
    real(r8) :: ts_ini(bounds%begp:bounds%endp)      ! stem temperature from beginning of time step [K]
    real(r8) :: ecidif                               ! excess energies [W/m2]
    real(r8) :: err(bounds%begp:bounds%endp)         ! balance error
    real(r8) :: erre                                 ! balance error
    real(r8) :: co2(bounds%begp:bounds%endp)         ! atmospheric co2 partial pressure (pa)
    real(r8) :: c13o2(bounds%begp:bounds%endp)       ! atmospheric c13o2 partial pressure (pa)
    real(r8) :: o2(bounds%begp:bounds%endp)          ! atmospheric o2 partial pressure (pa)
    real(r8) :: svpts(bounds%begp:bounds%endp)       ! saturation vapor pressure at t_veg (pa)
    real(r8) :: eah(bounds%begp:bounds%endp)         ! canopy air vapor pressure (pa)
    real(r8) :: s_node                               ! vol_liq/eff_porosity
    real(r8) :: smp_node                             ! matrix potential
    real(r8) :: smp_node_lf                          ! F. Li and S. Levis
    real(r8) :: vol_liq                              ! partial volume of liquid water in layer
    integer  :: itlef                                ! counter for leaf temperature iteration [-]
    integer  :: nmozsgn(bounds%begp:bounds%endp)     ! number of times stability changes sign
    real(r8) :: w                                    ! exp(-LSAI)
    real(r8) :: csoilcn                              ! interpolated csoilc for less than dense canopies
    real(r8) :: fm(bounds%begp:bounds%endp)          ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: wtshi                                ! sensible heat resistance for air, grnd and leaf [-]
    real(r8) :: wtsqi                                ! latent heat resistance for air, grnd and leaf [-]
    integer  :: j                                    ! soil/snow level index
    integer  :: p                                    ! patch index
    integer  :: c                                    ! column index
    integer  :: l                                    ! landunit index
    integer  :: g                                    ! gridcell index
    integer  :: fn                                   ! number of values in vegetated patch filter
    integer  :: filterp(bounds%endp-bounds%begp+1)   ! vegetated patch filter
    integer  :: fnorig                               ! number of values in patch filter copy
    integer  :: fporig(bounds%endp-bounds%begp+1)    ! temporary filter
    integer  :: fnold                                ! temporary copy of patch count
    integer  :: f                                    ! filter index
    logical  :: found                                ! error flag for canopy above forcing hgt
    integer  :: index                                ! patch index for error
    real(r8) :: egvf                                 ! effective green vegetation fraction
    real(r8) :: lt                                   ! elai+esai
    real(r8) :: U_ustar                              ! wind at canopy height divided by friction velocity (unitless)
    real(r8) :: U_ustar_ini                          ! initial guess of wind at canopy height divided by friction velocity (unitless)
    real(r8) :: U_ustar_prev                         ! wind at canopy height divided by friction velocity from the previous iteration (unitless)
    real(r8) :: ri                                   ! stability parameter for under canopy air (unitless)
    real(r8) :: csoilb                               ! turbulent transfer coefficient over bare soil (unitless)
    real(r8) :: ricsoilc                             ! modified transfer coefficient under dense canopy (unitless)
    real(r8) :: snow_depth_c                         ! critical snow depth to cover plant litter (m)
    real(r8) :: rdl                                  ! dry litter layer resistance for water vapor  (s/m)
    real(r8) :: elai_dl                              ! exposed (dry) plant litter area index
    real(r8) :: fsno_dl                              ! effective snow cover over plant litter
    real(r8) :: dayl_factor(bounds%begp:bounds%endp) ! scalar (0-1) for daylength effect on Vcmax
    ! If no unfrozen layers, put all in the top layer.
    real(r8) :: rootsum(bounds%begp:bounds%endp)
    real(r8) :: delt_snow
    real(r8) :: delt_soil
    real(r8) :: delt_h2osfc
    real(r8) :: lw_grnd
    real(r8) :: delq_snow
    real(r8) :: delq_soil
    real(r8) :: delq_h2osfc
    real(r8) :: dt_veg(bounds%begp:bounds%endp)          ! change in t_veg, last iteration (Kelvin)                              
    integer  :: jtop(bounds%begc:bounds%endc)            ! lbning
    integer  :: filterc_tmp(bounds%endp-bounds%begp+1)   ! temporary variable
    integer  :: ft                                       ! plant functional type index
    real(r8) :: h2ocan                                   ! total canopy water (mm H2O)
    real(r8) :: dt_veg_temp(bounds%begp:bounds%endp)
    integer, parameter :: iv=1                           ! index for first canopy layer (iwue calculation)
    real(r8) :: dbh(bounds%begp:bounds%endp)             ! diameter at breast height of vegetation
    real(r8) :: cp_leaf(bounds%begp:bounds%endp)         ! heat capacity of leaves
    real(r8) :: cp_stem(bounds%begp:bounds%endp)         ! heat capacity of stems
    real(r8) :: rstem(bounds%begp:bounds%endp)           ! stem resistance to heat transfer
    real(r8) :: dt_stem(bounds%begp:bounds%endp)         ! change in stem temperature
    real(r8) :: frac_rad_abs_by_stem(bounds%begp:bounds%endp)     ! fraction of incoming radiation absorbed by stems
    real(r8) :: lw_stem(bounds%begp:bounds%endp)         ! internal longwave stem
    real(r8) :: lw_leaf(bounds%begp:bounds%endp)         ! internal longwave leaf
    real(r8) :: sa_stem(bounds%begp:bounds%endp)         ! surface area stem m2/m2_ground
    real(r8) :: sa_leaf(bounds%begp:bounds%endp)         ! surface area leaf m2/m2_ground
    real(r8) :: sa_internal(bounds%begp:bounds%endp)     ! min(sa_stem,sa_leaf)
    real(r8) :: uuc(bounds%begp:bounds%endp)             ! undercanopy windspeed
    real(r8) :: carea_stem                               ! cross-sectional area of stem
    real(r8) :: dlrad_leaf                               ! Downward longwave radition from leaf
    real(r8) :: snocan_baseline(bounds%begp:bounds%endp)  ! baseline of snocan for use in truncate_small_values

    ! Indices for raw and rah
    integer, parameter :: above_canopy = 1         ! Above canopy
    integer, parameter :: below_canopy = 2         ! Below canopy
    ! Biomass heat storage tuning parameters
    ! These parameters can be used to account for differences
    ! in vegetation shape.
    real(r8), parameter :: k_vert     = 0.1_r8         !vertical distribution of stem
    real(r8), parameter :: k_cyl_vol  = 1.0_r8         !departure from cylindrical volume
    real(r8), parameter :: k_cyl_area = 1.0_r8         !departure from cylindrical area
    real(r8), parameter :: k_internal = 0.0_r8         !self-absorbtion of leaf/stem longwave
    real(r8), parameter :: min_stem_diameter = 0.05_r8 !minimum stem diameter for which to calculate stem interactions

    integer :: dummy_to_make_pgi_happy
    !------------------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(downreg_patch) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(leafn_patch) == (/bounds%endp/)), sourcefile, __LINE__)

    associate(                                                                    & 
         t_stem                 => temperature_inst%t_stem_patch                , & ! Output: [real(r8) (:)   ]  stem temperature (Kelvin)
         dhsdt_canopy           => energyflux_inst%dhsdt_canopy_patch           , & ! Output: [real(r8) (:)   ]  change in heat storage of stem (W/m**2) [+ to atm]
         soilresis              => soilstate_inst%soilresis_col                 , & ! Input:  [real(r8) (:)   ]  soil evaporative resistance
         snl                    => col%snl                                      , & ! Input:  [integer  (:)   ]  number of snow layers
         dayl                   => grc%dayl                                     , & ! Input:  [real(r8) (:)   ]  daylength (s)
         max_dayl               => grc%max_dayl                                 , & ! Input:  [real(r8) (:)   ]  maximum daylength for this grid cell (s)
         is_tree                => pftcon%is_tree                               , & ! Input:  tree patch or not
         is_shrub               => pftcon%is_shrub                              , & ! Input:  shrub patch or not
         dleaf                  => pftcon%dleaf                                 , & ! Input:  characteristic leaf dimension (m)
         dbh_param              => pftcon%dbh                                   , & ! Input:  diameter at brest height (m)
         slatop                 => pftcon%slatop                                , & ! SLA at top of canopy [m^2/gC]
         fbw                    => pftcon%fbw                                   , & ! Input:  fraction of biomass that is water
         nstem                  => pftcon%nstem                                 , & ! Input:  stem number density (#ind/m2)
         woody                  => pftcon%woody                                 , & ! Input:  woody flag
         rstem_per_dbh          => pftcon%rstem_per_dbh                         , & ! Input:  stem resistance per stem diameter (s/m**2)
         wood_density           => pftcon%wood_density                          , & ! Input:  dry wood density (kg/m3)

         z0v_Cr                 => pftcon%z0v_Cr                                , & ! Input:  roughness-element drag coefficient for Raupach92 parameterization (-)
         z0v_Cs                 => pftcon%z0v_Cs                                , & ! Input:  substrate-element drag coefficient for Raupach92 parameterization (-)
         z0v_c                  => pftcon%z0v_c                                 , & ! Input:  c parameter for Raupach92 parameterization (-)
         z0v_cw                 => pftcon%z0v_cw                                , & ! Input:  roughness sublayer depth coefficient for Raupach92 parameterization (-)
         z0v_LAIoff             => pftcon%z0v_LAIoff                            , & ! Input:  leaf area index offset for Raupach92 parameterization (-)
         z0v_LAImax             => pftcon%z0v_LAImax                            , & ! Input:  onset of over-sheltering for Raupach92 parameterization (-)

         forc_lwrad             => atm2lnd_inst%forc_lwrad_downscaled_col       , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)                       
         forc_q                 => wateratm2lndbulk_inst%forc_q_downscaled_col  , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)
         forc_pbot              => atm2lnd_inst%forc_pbot_downscaled_col        , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_th                => atm2lnd_inst%forc_th_downscaled_col          , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)                            
         forc_rho               => atm2lnd_inst%forc_rho_downscaled_col         , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)                                                     
         forc_t                 => atm2lnd_inst%forc_t_downscaled_col           , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)                                      
         forc_u                 => atm2lnd_inst%forc_u_grc                      , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)                        
         forc_v                 => atm2lnd_inst%forc_v_grc                      , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)                       
         forc_pco2              => atm2lnd_inst%forc_pco2_grc                   , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)                                             
         forc_pc13o2            => atm2lnd_inst%forc_pc13o2_grc                 , & ! Input:  [real(r8) (:)   ]  partial pressure c13o2 (Pa)                                           
         forc_po2               => atm2lnd_inst%forc_po2_grc                    , & ! Input:  [real(r8) (:)   ]  partial pressure o2 (Pa)                                              

         tc_ref2m               => humanindex_inst%tc_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (C)
         vap_ref2m              => humanindex_inst%vap_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m height vapor pressure (Pa)
         appar_temp_ref2m       => humanindex_inst%appar_temp_ref2m_patch       , & ! Output: [real(r8) (:)   ]  2 m apparent temperature (C)
         appar_temp_ref2m_r     => humanindex_inst%appar_temp_ref2m_r_patch     , & ! Output: [real(r8) (:)   ]  Rural 2 m apparent temperature (C)
         swbgt_ref2m            => humanindex_inst%swbgt_ref2m_patch            , & ! Output: [real(r8) (:)   ]  2 m Simplified Wetbulb Globe temperature (C)
         swbgt_ref2m_r          => humanindex_inst%swbgt_ref2m_r_patch          , & ! Output: [real(r8) (:)   ]  Rural 2 m Simplified Wetbulb Globe temperature (C)
         humidex_ref2m          => humanindex_inst%humidex_ref2m_patch          , & ! Output: [real(r8) (:)   ]  2 m Humidex (C)
         humidex_ref2m_r        => humanindex_inst%humidex_ref2m_r_patch        , & ! Output: [real(r8) (:)   ]  Rural 2 m Humidex (C)
         wbt_ref2m              => humanindex_inst%wbt_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m Stull Wet Bulb temperature (C)
         wbt_ref2m_r            => humanindex_inst%wbt_ref2m_r_patch            , & ! Output: [real(r8) (:)   ]  Rural 2 m Stull Wet Bulb temperature (C)
         wb_ref2m               => humanindex_inst%wb_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m Wet Bulb temperature (C)
         wb_ref2m_r             => humanindex_inst%wb_ref2m_r_patch             , & ! Output: [real(r8) (:)   ]  Rural 2 m Wet Bulb temperature (C)
         teq_ref2m              => humanindex_inst%teq_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m height Equivalent temperature (K)
         teq_ref2m_r            => humanindex_inst%teq_ref2m_r_patch            , & ! Output: [real(r8) (:)   ]  Rural 2 m Equivalent temperature (K)
         ept_ref2m              => humanindex_inst%ept_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m height Equivalent Potential temperature (K)
         ept_ref2m_r            => humanindex_inst%ept_ref2m_r_patch            , & ! Output: [real(r8) (:)   ]  Rural 2 m height Equivalent Potential temperature (K)
         discomf_index_ref2m    => humanindex_inst%discomf_index_ref2m_patch    , & ! Output: [real(r8) (:)   ]  2 m Discomfort Index temperature (C)
         discomf_index_ref2m_r  => humanindex_inst%discomf_index_ref2m_r_patch  , & ! Output: [real(r8) (:)   ]  Rural 2 m Discomfort Index temperature (C)
         discomf_index_ref2mS   => humanindex_inst%discomf_index_ref2mS_patch   , & ! Output: [real(r8) (:)   ]  2 m height Discomfort Index Stull temperature (C)
         discomf_index_ref2mS_r => humanindex_inst%discomf_index_ref2mS_r_patch , & ! Output: [real(r8) (:)   ]  Rural 2 m Discomfort Index Stull temperature (K)
         nws_hi_ref2m           => humanindex_inst%nws_hi_ref2m_patch           , & ! Output: [real(r8) (:)   ]  2 m NWS Heat Index (C)
         nws_hi_ref2m_r         => humanindex_inst%nws_hi_ref2m_r_patch         , & ! Output: [real(r8) (:)   ]  Rural 2 m NWS Heat Index (C)
         thip_ref2m             => humanindex_inst%thip_ref2m_patch             , & ! Output: [real(r8) (:)   ]  2 m Temperature Humidity Index Physiology (C)
         thip_ref2m_r           => humanindex_inst%thip_ref2m_r_patch           , & ! Output: [real(r8) (:)   ]  Rural 2 m Temperature Humidity Index Physiology (C)
         thic_ref2m             => humanindex_inst%thic_ref2m_patch             , & ! Output: [real(r8) (:)   ]  2 m Temperature Humidity Index Comfort (C)
         thic_ref2m_r           => humanindex_inst%thic_ref2m_r_patch           , & ! Output: [real(r8) (:)   ]  Rural 2 m Temperature Humidity Index Comfort (C)
         swmp65_ref2m           => humanindex_inst%swmp65_ref2m_patch           , & ! Output: [real(r8) (:)   ]  2 m Swamp Cooler temperature 65% effi (C)
         swmp65_ref2m_r         => humanindex_inst%swmp65_ref2m_r_patch         , & ! Output: [real(r8) (:)   ]  Rural 2 m Swamp Cooler temperature 65% effi (C)
         swmp80_ref2m           => humanindex_inst%swmp80_ref2m_patch           , & ! Output: [real(r8) (:)   ]  2 m Swamp Cooler temperature 80% effi (C)
         swmp80_ref2m_r         => humanindex_inst%swmp80_ref2m_r_patch         , & ! Output: [real(r8) (:)   ]  Rural 2 m Swamp Cooler temperature 80% effi (C)

         sabv                   => solarabs_inst%sabv_patch                     , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)                       
         par_z_sun              => solarabs_inst%parsun_z_patch                 , & ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)

         frac_veg_nosno         => canopystate_inst%frac_veg_nosno_patch        , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         elai                   => canopystate_inst%elai_patch                  , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow                        
         esai                   => canopystate_inst%esai_patch                  , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow                        
         laisun                 => canopystate_inst%laisun_patch                , & ! Input:  [real(r8) (:)   ]  sunlit leaf area                                                      
         laisha                 => canopystate_inst%laisha_patch                , & ! Input:  [real(r8) (:)   ]  shaded leaf area                                                      
         displa                 => canopystate_inst%displa_patch                , & ! Input:  [real(r8) (:)   ]  displacement height (m)                                               
         stem_biomass           => canopystate_inst%stem_biomass_patch          , & ! Output: [real(r8) (:)   ]  Aboveground stem biomass  (kg/m**2)
         leaf_biomass           => canopystate_inst%leaf_biomass_patch          , & ! Output: [real(r8) (:)   ]  Aboveground leaf biomass  (kg/m**2)
         htop                   => canopystate_inst%htop_patch                  , & ! Input:  [real(r8) (:)   ]  canopy top(m)                                                         
         dleaf_patch            => canopystate_inst%dleaf_patch                 , & ! Output: [real(r8) (:)   ]  mean leaf diameter for this patch/pft
         
         watsat                 => soilstate_inst%watsat_col                    , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)   (constant)                     
         watdry                 => soilstate_inst%watdry_col                    , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=0                      (constant)                                        
         watopt                 => soilstate_inst%watopt_col                    , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=1                      (constant)                                      
         eff_porosity           => soilstate_inst%eff_porosity_col              , & ! Output: [real(r8) (:,:) ]  effective soil porosity
         soilbeta               => soilstate_inst%soilbeta_col                  , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity                               

         u10_clm                => frictionvel_inst%u10_clm_patch               , & ! Input:  [real(r8) (:)   ]  10 m height winds (m/s)
         forc_hgt_t             => atm2lnd_inst%forc_hgt_t_grc                  , & ! Input:  [real(r8) (:)   ] observational height of temperature [m]
         forc_hgt_u             => atm2lnd_inst%forc_hgt_u_grc                  , & ! Input:  [real(r8) (:)   ] observational height of wind [m]
         forc_hgt_q             => atm2lnd_inst%forc_hgt_q_grc                  , & ! Input:  [real(r8) (:)   ] observational height of specific humidity [m]
         forc_hgt_t_patch       => frictionvel_inst%forc_hgt_t_patch            , & ! Output: [real(r8) (:)   ] observational height of temperature at patch level [m]
         forc_hgt_q_patch       => frictionvel_inst%forc_hgt_q_patch            , & ! Output: [real(r8) (:)   ] observational height of specific humidity at patch level [m]
         forc_hgt_u_patch       => frictionvel_inst%forc_hgt_u_patch            , & ! Output:  [real(r8) (:)   ]  observational height of wind at patch level [m]
         z0mg                   => frictionvel_inst%z0mg_col                    , & ! Input:  [real(r8) (:)   ]  roughness length of ground, momentum [m]                              
         zetamax                => frictionvel_inst%zetamaxstable               , & ! Input:  [real(r8)       ]  max zeta value under stable conditions
         ram1                   => frictionvel_inst%ram1_patch                  , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)                                        
         z0mv                   => frictionvel_inst%z0mv_patch                  , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, momentum [m]                        
         z0hv                   => frictionvel_inst%z0hv_patch                  , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, sensible heat [m]                   
         z0qv                   => frictionvel_inst%z0qv_patch                  , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, latent heat [m]                     
         rb1                    => frictionvel_inst%rb1_patch                   , & ! Output: [real(r8) (:)   ]  boundary layer resistance (s/m)                                       

         t_h2osfc               => temperature_inst%t_h2osfc_col                , & ! Input:  [real(r8) (:)   ]  surface water temperature                                             
         t_soisno               => temperature_inst%t_soisno_col                , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                                           
         t_grnd                 => temperature_inst%t_grnd_col                  , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]                                        
         thv                    => temperature_inst%thv_col                     , & ! Input:  [real(r8) (:)   ]  virtual potential temperature (kelvin)                                
         thm                    => temperature_inst%thm_patch                   , & ! Input:  [real(r8) (:)   ]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)                  
         emv                    => temperature_inst%emv_patch                   , & ! Input:  [real(r8) (:)   ]  vegetation emissivity                                                     
         emg                    => temperature_inst%emg_col                     , & ! Input:  [real(r8) (:)   ]  vegetation emissivity                                                 
         t_veg                  => temperature_inst%t_veg_patch                 , & ! Output: [real(r8) (:)   ]  vegetation temperature (Kelvin)                                       
         t_ref2m                => temperature_inst%t_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)                           
         t_ref2m_r              => temperature_inst%t_ref2m_r_patch             , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface air temperature (Kelvin)                     
         t_skin_patch           => temperature_inst%t_skin_patch                , & ! Output: [real(r8) (:)   ]  patch skin temperature (K)  

         frac_h2osfc            => waterdiagnosticbulk_inst%frac_h2osfc_col              , & ! Input:  [real(r8) (:)   ]  fraction of surface water                                             
         fwet                   => waterdiagnosticbulk_inst%fwet_patch                   , & ! Input:  [real(r8) (:)   ]  fraction of canopy that is wet (0 to 1)                               
         fdry                   => waterdiagnosticbulk_inst%fdry_patch                   , & ! Input:  [real(r8) (:)   ]  fraction of foliage that is green and dry [-]                         
         frac_sno               => waterdiagnosticbulk_inst%frac_sno_eff_col             , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)                           
         snow_depth             => waterdiagnosticbulk_inst%snow_depth_col               , & ! Input:  [real(r8) (:)   ]  snow height (m)                                                       
         qg_snow                => waterdiagnosticbulk_inst%qg_snow_col                  , & ! Input:  [real(r8) (:)   ]  specific humidity at snow surface [kg/kg]                             
         qg_soil                => waterdiagnosticbulk_inst%qg_soil_col                  , & ! Input:  [real(r8) (:)   ]  specific humidity at soil surface [kg/kg]                             
         qg_h2osfc              => waterdiagnosticbulk_inst%qg_h2osfc_col                , & ! Input:  [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]                           
         qg                     => waterdiagnosticbulk_inst%qg_col                       , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]                           
         dqgdT                  => waterdiagnosticbulk_inst%dqgdT_col                    , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"                                        

         h2osoi_ice             => waterstatebulk_inst%h2osoi_ice_col               , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                                    
         h2osoi_vol             => waterstatebulk_inst%h2osoi_vol_col               , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] by F. Li and S. Levis
         h2osoi_liq             => waterstatebulk_inst%h2osoi_liq_col               , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                                                
         h2osoi_liqvol          => waterdiagnosticbulk_inst%h2osoi_liqvol_col            , & ! Output: [real(r8) (:,:) ]  volumetric liquid water (v/v) 
         snocan                 => waterstatebulk_inst%snocan_patch                 , & ! Output: [real(r8) (:)   ]  canopy snow (mm H2O)                                                 
         liqcan                 => waterstatebulk_inst%liqcan_patch                 , & ! Output: [real(r8) (:)   ]  canopy liquid (mm H2O)                                                 

         q_ref2m                => waterdiagnosticbulk_inst%q_ref2m_patch                , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)                          
         rh_ref2m_r             => waterdiagnosticbulk_inst%rh_ref2m_r_patch             , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface relative humidity (%)                        
         rh_ref2m               => waterdiagnosticbulk_inst%rh_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)                              
         rhaf                   => waterdiagnosticbulk_inst%rh_af_patch                  , & ! Output: [real(r8) (:)   ]  fractional humidity of canopy air [dimensionless]                     
         vpd_ref2m              => waterdiagnosticbulk_inst%vpd_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m height surface vapor pressure deficit (Pa)
         iwue_ln                => waterdiagnosticbulk_inst%iwue_ln_patch                , & ! Output: [real(r8) (:)   ]  local noon ecosystem-scale inherent water use efficiency (gC kgH2O-1 hPa)

         qflx_tran_veg          => waterfluxbulk_inst%qflx_tran_veg_patch           , & ! Output: [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)                      
         qflx_evap_veg          => waterfluxbulk_inst%qflx_evap_veg_patch           , & ! Output: [real(r8) (:)   ]  vegetation evaporation (mm H2O/s) (+ = to atm)                        
         qflx_evap_soi          => waterfluxbulk_inst%qflx_evap_soi_patch           , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)                              
         qflx_ev_snow           => waterfluxbulk_inst%qflx_ev_snow_patch            , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (mm H2O/s) [+ to atm]                        
         qflx_ev_soil           => waterfluxbulk_inst%qflx_ev_soil_patch            , & ! Output: [real(r8) (:)   ]  evaporation flux from soil (mm H2O/s) [+ to atm]                        
         qflx_ev_h2osfc         => waterfluxbulk_inst%qflx_ev_h2osfc_patch          , & ! Output: [real(r8) (:)   ]  evaporation flux from h2osfc (mm H2O/s) [+ to atm]                      
         gs_mol_sun             => photosyns_inst%gs_mol_sun_patch              , & ! Input: [real(r8) (:)   ]  patch sunlit leaf stomatal conductance (umol H2O/m**2/s)
         gs_mol_sha             => photosyns_inst%gs_mol_sha_patch              , & ! Input: [real(r8) (:)   ]  patch shaded leaf stomatal conductance (umol H2O/m**2/s)
         rssun                  => photosyns_inst%rssun_patch                   , & ! Output: [real(r8) (:)   ]  leaf sunlit stomatal resistance (s/m) (output from Photosynthesis)
         rssha                  => photosyns_inst%rssha_patch                   , & ! Output: [real(r8) (:)   ]  leaf shaded stomatal resistance (s/m) (output from Photosynthesis)
         fpsn                   => photosyns_inst%fpsn_patch                    , & ! Input:  [real(r8) (:)   ]  photosynthesis (umol CO2 /m**2 /s)

         grnd_ch4_cond          => ch4_inst%grnd_ch4_cond_patch                 , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s] 

         htvp                   => energyflux_inst%htvp_col                     , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) [J/kg] (constant)                      
         btran                  => energyflux_inst%btran_patch                  , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1)                                 
         rresis                 => energyflux_inst%rresis_patch                 , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd)                          
         taux                   => energyflux_inst%taux_patch                   , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)                                  
         tauy                   => energyflux_inst%tauy_patch                   , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)                                  
         canopy_cond            => energyflux_inst%canopy_cond_patch            , & ! Output: [real(r8) (:)   ]  tracer conductance for canopy [m/s] 
         cgrnds                 => energyflux_inst%cgrnds_patch                 , & ! Output: [real(r8) (:)   ]  deriv. of soil sensible heat flux wrt soil temp [w/m2/k]              
         cgrndl                 => energyflux_inst%cgrndl_patch                 , & ! Output: [real(r8) (:)   ]  deriv. of soil latent heat flux wrt soil temp [w/m**2/k]              
         dlrad                  => energyflux_inst%dlrad_patch                  , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy [W/m2]                   
         ulrad                  => energyflux_inst%ulrad_patch                  , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy [W/m2]                     
         cgrnd                  => energyflux_inst%cgrnd_patch                  , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]                  
         eflx_sh_snow           => energyflux_inst%eflx_sh_snow_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]                      
         eflx_sh_h2osfc         => energyflux_inst%eflx_sh_h2osfc_patch         , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]                      
         eflx_sh_soil           => energyflux_inst%eflx_sh_soil_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]                      
         eflx_sh_stem           => energyflux_inst%eflx_sh_stem_patch            , & ! Output: [real(r8) (:)   ]  sensible heat flux from stems (W/m**2) [+ to atm]
         eflx_sh_veg            => energyflux_inst%eflx_sh_veg_patch            , & ! Output: [real(r8) (:)   ]  sensible heat flux from leaves (W/m**2) [+ to atm]                    
         eflx_sh_grnd           => energyflux_inst%eflx_sh_grnd_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]                    
         rah1                   => frictionvel_inst%rah1_patch                  , & ! Output: [real(r8) (:)   ]  aerodynamical  resistance [s/m]
         rah2                   => frictionvel_inst%rah2_patch                  , & ! Output: [real(r8) (:)   ]  aerodynamical  resistance [s/m]
         raw1                   => frictionvel_inst%raw1_patch                  , & ! Output: [real(r8) (:)   ]  aerodynamical  resistance [s/m]
         raw2                   => frictionvel_inst%raw2_patch                  , & ! Output: [real(r8) (:)   ]  aerodynamical resistance [s/m]
         ustar                  => frictionvel_inst%ustar_patch                 , & ! Output: [real(r8) (:)   ]  friction velocity [m/s]
         um                     => frictionvel_inst%um_patch                    , & ! Output: [real(r8) (:)   ]  wind speed including the stablity effect [m/s]
         uaf                    => frictionvel_inst%uaf_patch                   , & ! Output: [real(r8) (:)   ]  canopy air speed [m/s]
         taf                    => frictionvel_inst%taf_patch                   , & ! Output: [real(r8) (:)   ]  canopy air temperature [K]
         qaf                    => frictionvel_inst%qaf_patch                   , & ! Output: [real(r8) (:)   ]  canopy air humidity [kg/kg]
         obu                    => frictionvel_inst%obu_patch                   , & ! Output: [real(r8) (:)   ]  Monin-Obukhov length [m]
         zeta                   => frictionvel_inst%zeta_patch                  , & ! Output: [real(r8) (:)   ]  dimensionless stability parameter 
         vpd                    => frictionvel_inst%vpd_patch                   , & ! Output: [real(r8) (:)   ]  vapor pressure deficit [Pa]
         num_iter               => frictionvel_inst%num_iter_patch              , & ! Output: [real(r8) (:)   ]  number of iterations

         begp                   => bounds%begp                                  , &
         endp                   => bounds%endp                                  , &
         begg                   => bounds%begg                                  , &
         endg                   => bounds%endg                                    &
         )
      if (use_hydrstress) then
        bsun                    => energyflux_inst%bsun_patch                       ! Output: [real(r8) (:)   ]  sunlit canopy transpiration wetness factor (0 to 1)
        bsha                    => energyflux_inst%bsha_patch                       ! Output: [real(r8) (:)   ]  sunlit canopy transpiration wetness factor (0 to 1)
      end if

      
      ! Determine step size

      dtime = get_step_size_real()

      ! Make a local copy of the exposedvegp filter. With the current implementation,
      ! this is needed because the filter is modified in the iteration loop.
      !
      ! TODO(wjs, 2014-09-24) Determine if this is really needed. I suspect that we could
      ! do away with either this temporary fn/filterp, or the temporary fnorig/fporig,
      ! with one of these simply using the passed-in filter (num_exposedvegp /
      ! filter_exposedvegp)

      fn = num_exposedvegp
      filterp(1:fn) = filter_exposedvegp(1:fn)
      
      ! -----------------------------------------------------------------
      ! Time step initialization of photosynthesis variables
      ! -----------------------------------------------------------------

      call photosyns_inst%TimeStepInit(bounds)


      ! -----------------------------------------------------------------
      ! Prep some IO variables and some checks on patch pointers if FATES
      ! is running. 
      ! Filter explanation: The patch filter in this routine identifies all
      ! non-lake, non-urban patches that are not covered by ice. The
      ! filter is set over a few steps:
      !
      ! 1a) for CN: 
      !             clm_drv() -> 
      !             bgc_vegetation_inst%EcosystemDynamicsPostDrainage() ->
      !             CNVegStructUpdate()
      !    if(elai(p)+esai(p)>0) frac_veg_nosno_alb(p) = 1
      !    
      ! 1b) for FATES:
      !              clm_drv() -> 
      !              clm_fates%dynamics_driv() -> 
      !              ed_clm_link() -> 
      !              ed_clm_leaf_area_profile():
      !    if(elai(p)+esai(p)>0) frac_veg_nosno_alb(p) = 1
      !
      ! 2) during clm_drv()->clm_drv_init():
      !    frac_veg_nosno_alb(p) is then combined with the active(p)
      !    flag via union to create frac_veg_nosno_patch(p)
      ! 3) immediately after, during clm_drv()->setExposedvegpFilter()
      !    the list used here "exposedvegp(fe)" is incremented if 
      !    frac_veg_nosno_patch > 0
      ! -----------------------------------------------------------------

      if (use_fates) then
         call clm_fates%prep_canopyfluxes(nc, fn, filterp, photosyns_inst)
      end if

      ! Initialize

      do f = 1, fn
         p = filterp(f)
         del(p)    = 0._r8  ! change in leaf temperature from previous iteration
         efeb(p)   = 0._r8  ! latent head flux from leaf for previous iteration
         wtlq0(p)  = 0._r8
         wtalq(p)  = 0._r8
         wtgq(p)   = 0._r8
         wtaq0(p)  = 0._r8
         obuold(p) = 0._r8
         btran(p)  = btran0
         dhsdt_canopy(p) = 0._r8
         eflx_sh_stem(p) = 0._r8
      end do
      !
      ! Calculate biomass heat capacities
      !
      if(use_biomass_heat_storage) then
bioms:   do f = 1, fn
            p = filterp(f)

            ! fraction of stem receiving incoming radiation
            frac_rad_abs_by_stem(p) = (esai(p))/(elai(p)+esai(p))

            ! when elai = 0, do not multiply by k_vert (i.e. frac_rad_abs_by_stem = 1)
            if(elai(p) > 0._r8) frac_rad_abs_by_stem(p) = k_vert * frac_rad_abs_by_stem(p)

            ! if using Satellite Phenology mode, use values in parameter file
            ! otherwise calculate dbh from stem biomass
            if(use_cn) then
               if(stem_biomass(p) > 0._r8) then 
                  dbh(p) = 2._r8 * sqrt(stem_biomass(p) * (1._r8 - fbw(patch%itype(p))) &
                       / ( shr_const_pi * htop(p) * k_cyl_vol & 
                       * nstem(patch%itype(p)) *  wood_density(patch%itype(p))))
               else
                  dbh(p) = 0._r8
               endif
            else
               dbh(p) = dbh_param(patch%itype(p))
            endif

            ! leaf and stem surface area
            sa_leaf(p) = elai(p)
            ! double in spirit of full surface area for sensible heat
            sa_leaf(p) = 2._r8*sa_leaf(p)

            ! Surface area for stem
            sa_stem(p) = nstem(patch%itype(p))*(htop(p)*shr_const_pi*dbh(p))
            ! adjust for departure of cylindrical stem model
            sa_stem(p) = k_cyl_area * sa_stem(p)

            !
            ! only calculate separate leaf/stem heat capacity for trees
            ! and shrubs if dbh is greater than some minimum value
            ! (set surface area for stem, and fraction absorbed by stem to zero)
            if(.not.(is_tree(patch%itype(p)) .or. is_shrub(patch%itype(p))) &
                 .or. dbh(p) < min_stem_diameter) then
               frac_rad_abs_by_stem(p) = 0.0_r8
               sa_stem(p) = 0.0_r8
               sa_leaf(p) = sa_leaf(p) + esai(p)
            endif

            ! if using Satellite Phenology mode, calculate leaf and stem biomass
            if(.not. use_cn) then
               ! 2gbiomass/gC * (1/SLA) * 1e-3 = kg dry mass/m2(leaf)
               leaf_biomass(p) = (1.e-3_r8*c_to_b/slatop(patch%itype(p))) &
                    * max(0.01_r8, 0.5_r8*sa_leaf(p)) &
                    / (1._r8-fbw(patch%itype(p)))
               ! cross-sectional area of stems
               carea_stem = shr_const_pi * (dbh(p)*0.5_r8)**2
               stem_biomass(p) = carea_stem * htop(p) * k_cyl_vol &
                    * nstem(patch%itype(p)) * wood_density(patch%itype(p)) &
                    /(1._r8-fbw(patch%itype(p)))
            endif

            ! internal longwave fluxes between leaf and stem
            ! (use same area of interaction i.e. ignore leaf <-> leaf)
            sa_internal(p) = min(sa_leaf(p),sa_stem(p))
            sa_internal(p) = k_internal * sa_internal(p)

            ! calculate specify heat capacity of vegetation
            ! as weighted averaged of dry biomass and water
            ! lma_dry has units of kg dry mass/m2 here
            ! (Appendix B of Bonan et al., GMD, 2018) 

            cp_leaf(p)  = leaf_biomass(p) * (c_dry_biomass*(1._r8-fbw(patch%itype(p))) + (fbw(patch%itype(p)))*c_water)

            ! cp-stem will have units J/k/ground_area
            cp_stem(p) = stem_biomass(p) * (c_dry_biomass*(1._r8-fbw(patch%itype(p))) + (fbw(patch%itype(p)))*c_water)
            ! adjust for departure from cylindrical stem model
            cp_stem(p) = k_cyl_vol * cp_stem(p)

            ! resistance between internal stem temperature and canopy air 
            rstem(p) = rstem_per_dbh(patch%itype(p))*dbh(p)

         enddo bioms
      else
        ! Otherwise set biomass heat storage terms to zero
        do f = 1, fn
            p = filterp(f)
            sa_leaf(p)              = (elai(p)+esai(p))
            frac_rad_abs_by_stem(p) = 0._r8
            sa_stem(p)              = 0._r8
            sa_internal(p)          = 0._r8
            cp_leaf(p)              = 0._r8
            cp_stem(p)              = 0._r8
            rstem(p)                = 0._r8
        end do
      end if


      ! calculate daylength control for Vcmax
      do f = 1, fn
         p=filterp(f)
         g=patch%gridcell(p)
         ! calculate dayl_factor as the ratio of (current:max dayl)^2
         ! set a minimum of 0.01 (1%) for the dayl_factor
         dayl_factor(p)=min(1._r8,max(0.01_r8,(dayl(g)*dayl(g))/(max_dayl(g)*max_dayl(g))))
      end do

      rb1(begp:endp) = 0._r8

      !assign the temporary filter
      do f = 1, fn
         p = filterp(f)
         filterc_tmp(f)=patch%column(p)
      enddo
      
      !compute effective soil porosity
      call calc_effective_soilporosity(bounds,                          &
            ubj = nlevgrnd,                                              &
            numf = fn,                                                   &
            filter = filterc_tmp(1:fn),                                  &
            watsat = watsat(bounds%begc:bounds%endc, 1:nlevgrnd),        &
            h2osoi_ice = h2osoi_ice(bounds%begc:bounds%endc,1:nlevgrnd), &
            denice = denice,                                             &
            eff_por=eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd) )
      
      !compute volumetric liquid water content
      jtop(bounds%begc:bounds%endc) = 1
      
      call calc_volumetric_h2oliq(bounds,                                    &
            jtop = jtop(bounds%begc:bounds%endc),                             &
            lbj = 1,                                                          &
            ubj = nlevgrnd,                                                   &
            numf = fn,                                                        &
            filter = filterc_tmp(1:fn),                                       &
            eff_porosity = eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd), &
            h2osoi_liq = h2osoi_liq(bounds%begc:bounds%endc, 1:nlevgrnd),     &
            denh2o = denh2o,                                                  &
            vol_liq = h2osoi_liqvol(bounds%begc:bounds%endc, 1:nlevgrnd) )
      
      !set up perchroot options
      call set_perchroot_opt(perchroot, perchroot_alt)

      ! --------------------------------------------------------------------------
      ! if this is a FATES simulation
      ! ask fates to calculate btran functions and distribution of uptake
      ! this will require boundary conditions from CLM, boundary conditions which
      ! may only be available from a smaller subset of patches that meet the
      ! exposed veg.  
      ! calc_root_moist_stress already calculated root soil water stress 'rresis'
      ! this is the input boundary condition to calculate the transpiration
      ! wetness factor btran and the root weighting factors for FATES.  These
      ! values require knowledge of the belowground root structure.
      ! --------------------------------------------------------------------------
      
      if(use_fates)then
         call clm_fates%wrap_btran(nc, fn, filterc_tmp(1:fn), soilstate_inst, &
               waterdiagnosticbulk_inst, temperature_inst, energyflux_inst, soil_water_retention_curve)
         
      else
         
         !calculate root moisture stress
         call calc_root_moist_stress(bounds,     &
            nlevgrnd = nlevgrnd,               &
            fn = fn,                           &
            filterp = filterp,                 &
            active_layer_inst=active_layer_inst, &
            energyflux_inst=energyflux_inst,   &
            soilstate_inst=soilstate_inst,     &
            temperature_inst=temperature_inst, &
            waterstatebulk_inst=waterstatebulk_inst,   &
            waterdiagnosticbulk_inst=waterdiagnosticbulk_inst,   &
              soil_water_retention_curve=soil_water_retention_curve)
     
      end if

      !! Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)
         g = patch%gridcell(p)

         select case (z0param_method)
         case ('ZengWang2007')
            lt = min(elai(p)+esai(p), tlsai_crit)
            egvf =(1._r8 - alpha_aero * exp(-lt)) / (1._r8 - alpha_aero * exp(-tlsai_crit))
            displa(p) = egvf * displa(p)
            z0mv(p)   = exp(egvf * log(z0mv(p)) + (1._r8 - egvf) * log(z0mg(c)))

         case ('Meier2022')
            lt = max(1.e-5_r8, elai(p) + esai(p))
            displa(p) = htop(p) * (1._r8 - (1._r8 - exp(-(cd1_param * lt)**0.5_r8)) / (cd1_param*lt)**0.5_r8)

            lt = min(lt,z0v_LAImax(patch%itype(p)))
            delt = 2._r8
            ! Reminder that (...)**(-0.5) = 1 / sqrt(...)
            U_ustar_ini = (z0v_Cs(patch%itype(p)) + z0v_Cr(patch%itype(p)) * lt * 0.5_r8)**(-0.5_r8) &
                      *z0v_c(patch%itype(p)) * lt * 0.25_r8
            U_ustar = U_ustar_ini

            do while (delt > 1.e-4_r8)
               U_ustar_prev = U_ustar
               U_ustar = U_ustar_ini * exp(U_ustar_prev)
               delt = abs(U_ustar - U_ustar_prev)
            end do

            U_ustar = 4._r8 * U_ustar / lt / z0v_c(patch%itype(p))

            z0mv(p) = htop(p) * (1._r8 - displa(p) / htop(p)) * exp(-vkc * U_ustar + &
                      log(z0v_cw(patch%itype(p))) - 1._r8 + z0v_cw(patch%itype(p))**(-1._r8))


          case default
            write(iulog,*) 'ERROR: unknown z0para_method: ', z0param_method
            call endrun(msg = 'unknown z0param_method', additional_msg = errMsg(sourcefile, __LINE__))
          end select

          z0hv(p)   = z0mv(p)
          z0qv(p)   = z0mv(p)

          ! Update the forcing heights
          forc_hgt_u_patch(p) = forc_hgt_u(g) + z0mv(p) + displa(p)
          forc_hgt_t_patch(p) = forc_hgt_t(g) + z0hv(p) + displa(p)
          forc_hgt_q_patch(p) = forc_hgt_q(g) + z0qv(p) + displa(p)

      end do

      found = .false.
      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)
         g = patch%gridcell(p)

         ! Net absorbed longwave radiation by canopy and ground
         ! =air+bir*t_veg**4+cir*t_grnd(c)**4

         air(p) =   emv(p) * (1._r8+(1._r8-emv(p))*(1._r8-emg(c))) * forc_lwrad(c)
         bir(p) = - (2._r8-emv(p)*(1._r8-emg(c))) * emv(p) * sb
         cir(p) =   emv(p)*emg(c)*sb

         ! Saturated vapor pressure, specific humidity, and their derivatives
         ! at the leaf surface

         call QSat (t_veg(p), forc_pbot(c), qsatl(p), &
              es = el(p), &
              qsdT = qsatldT(p))

         ! Determine atmospheric co2 and o2

         co2(p) = forc_pco2(g)
         o2(p)  = forc_po2(g)

         if ( use_c13 ) then
            c13o2(p) = forc_pc13o2(g)
         end if

         ! Initialize flux profile

         nmozsgn(p) = 0

         taf(p) = (t_grnd(c) + thm(p))/2._r8
         qaf(p) = (forc_q(c)+qg(c))/2._r8

         ur(p) = max(params_inst%wind_min,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
         dth(p) = thm(p)-taf(p)
         dqh(p) = forc_q(c)-qaf(p)
         delq(p) = qg(c) - qaf(p)
         dthv(p) = dth(p)*(1._r8+0.61_r8*forc_q(c))+0.61_r8*forc_th(c)*dqh(p)
         zldis(p) = forc_hgt_u_patch(p) - displa(p)

         ! Check to see if the forcing height is below the canopy height
         if (zldis(p) < 0._r8) then
            found = .true.
            index = p
         end if

      end do

      if (found) then
         if ( .not. use_fates ) then
            write(iulog,*)'Error: Forcing height is below canopy height for patch index '
            call endrun(subgrid_index=index, subgrid_level=subgrid_level_patch, msg=errmsg(sourcefile, __LINE__))
         end if
      end if

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)

         ! Initialize Monin-Obukhov length and wind speed

         call frictionvel_inst%MoninObukIni(ur(p), thv(c), dthv(p), zldis(p), z0mv(p), um(p), obu(p))
         num_iter(p) = 0

         ! Record initial veg/stem temperatures
         tl_ini(p) = t_veg(p)
         ts_ini(p) = t_stem(p)

      end do

      ! Set counter for leaf temperature iteration (itlef)

      itlef = 0    
      fnorig = fn
      fporig(1:fn) = filterp(1:fn)

      ! Begin stability iteration

      call t_startf('can_iter')
      ITERATION : do while (itlef <= itmax_canopy_fluxes .and. fn > 0)

         ! Determine friction velocity, and potential temperature and humidity
         ! profiles of the surface boundary layer

         call frictionvel_inst%FrictionVelocity (begp, endp, fn, filterp, &
              displa(begp:endp), z0mv(begp:endp), z0hv(begp:endp), z0qv(begp:endp), &
              obu(begp:endp), itlef+1, ur(begp:endp), um(begp:endp), ustar(begp:endp), &
              temp1(begp:endp), temp2(begp:endp), temp12m(begp:endp), temp22m(begp:endp), fm(begp:endp))

         do f = 1, fn
            p = filterp(f)
            c = patch%column(p)
            g = patch%gridcell(p)

            tlbef(p) = t_veg(p)
            del2(p) = del(p)

            ! Determine aerodynamic resistances

            ram1(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
            rah(p,above_canopy) = 1._r8/(temp1(p)*ustar(p))
            raw(p,above_canopy) = 1._r8/(temp2(p)*ustar(p))

            ! Bulk boundary layer resistance of leaves

            uaf(p) = um(p)*sqrt( 1._r8/(ram1(p)*um(p)) )

            ! empirical undercanopy wind speed
            uuc(p) = min(0.4_r8,(0.03_r8*um(p)/ustar(p)))

            ! Use pft parameter for leaf characteristic width
            ! dleaf_patch if this is not an fates patch.
            ! Otherwise, the value has already been loaded
            ! during the FATES dynamics call
            if(.not.patch%is_fates(p)) then  
               dleaf_patch(p) = dleaf(patch%itype(p))
            end if

            cf  = params_inst%cv / (sqrt(uaf(p)) * sqrt(dleaf_patch(p)))
            rb(p)  = 1._r8/(cf*uaf(p))
            rb1(p) = rb(p)

            ! Parameterization for variation of csoilc with canopy density from
            ! X. Zeng, University of Arizona

            w = exp(-(elai(p)+esai(p)))

            ! changed by K.Sakaguchi from here
            ! transfer coefficient over bare soil is changed to a local variable
            ! just for readability of the code (from line 680)
            csoilb = vkc / (params_inst%a_coef * (z0mg(c) * uaf(p) / nu_param)**params_inst%a_exp)

            !compute the stability parameter for ricsoilc  ("S" in Sakaguchi&Zeng,2008)

            ri = ( grav*htop(p) * (taf(p) - t_grnd(c)) ) / (taf(p) * uaf(p) **2.00_r8)

            ! modify csoilc value (0.004) if the under-canopy is in stable condition
            if (use_undercanopy_stability .and. (taf(p) - t_grnd(c) ) > 0._r8) then
               ! decrease the value of csoilc by dividing it with (1+gamma*min(S, 10.0))
               ! ria ("gmanna" in Sakaguchi&Zeng, 2008) is a constant (=0.5)
               ricsoilc = params_inst%csoilc / (1.00_r8 + ria*min( ri, 10.0_r8) )
               csoilcn = csoilb*w + ricsoilc*(1._r8-w)
            else
               csoilcn = csoilb*w + params_inst%csoilc*(1._r8-w)
            end if

            !! Sakaguchi changes for stability formulation ends here

            if (use_biomass_heat_storage) then
               ! use uuc for ground fluxes (keep uaf for canopy terms)
               rah(p,below_canopy) = 1._r8/(csoilcn*uuc(p))
            else
               rah(p,below_canopy) = 1._r8/(csoilcn*uaf(p))
            endif

            raw(p,below_canopy) = rah(p,below_canopy)
            if (use_lch4) then
               grnd_ch4_cond(p) = 1._r8/(raw(p,above_canopy)+raw(p,below_canopy))
            end if

            ! Stomatal resistances for sunlit and shaded fractions of canopy.
            ! Done each iteration to account for differences in eah, tv.

            svpts(p) = el(p)                         ! pa
            eah(p) = forc_pbot(c) * qaf(p) / 0.622_r8   ! pa
            rhaf(p) = eah(p)/svpts(p)
            ! variables for history fields
            rah1(p)  = rah(p,above_canopy)
            raw1(p)  = raw(p,above_canopy)
            rah2(p)  = rah(p,below_canopy)
            raw2(p)  = raw(p,below_canopy)
            vpd(p)  = max((svpts(p) - eah(p)), 50._r8) * 0.001_r8

         end do

         if ( use_fates ) then      
            
            call clm_fates%wrap_photosynthesis(nc, bounds, fn, filterp(1:fn), &
                 svpts(begp:endp), eah(begp:endp), o2(begp:endp), &
                 co2(begp:endp), rb(begp:endp), dayl_factor(begp:endp), &
                 atm2lnd_inst, temperature_inst, canopystate_inst, photosyns_inst)

         else ! not use_fates

            if ( use_hydrstress ) then
               call PhotosynthesisHydraulicStress (bounds, fn, filterp, &
                    svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), bsun(begp:endp), &
                    bsha(begp:endp), btran(begp:endp), dayl_factor(begp:endp), leafn_patch(begp:endp), &
                    qsatl(begp:endp), qaf(begp:endp),     &
                    atm2lnd_inst, temperature_inst, soilstate_inst, waterdiagnosticbulk_inst, surfalb_inst, solarabs_inst, &
                    canopystate_inst, ozone_inst, photosyns_inst, waterfluxbulk_inst, &
                    froot_carbon(begp:endp), croot_carbon(begp:endp))
            else
               call Photosynthesis (bounds, fn, filterp, &
                    svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), btran(begp:endp), &
                    dayl_factor(begp:endp), leafn_patch(begp:endp), &
                    atm2lnd_inst, temperature_inst, surfalb_inst, solarabs_inst, &
                    canopystate_inst, ozone_inst, photosyns_inst, phase='sun')
            endif

            if ( use_cn .and. use_c13 ) then
               call Fractionation (bounds, fn, filterp, downreg_patch(begp:endp), &
                    atm2lnd_inst, canopystate_inst, solarabs_inst, surfalb_inst, photosyns_inst, &
                    phase='sun')
            endif

            if ( .not.(use_hydrstress) ) then
               call Photosynthesis (bounds, fn, filterp, &
                    svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), btran(begp:endp), &
                    dayl_factor(begp:endp), leafn_patch(begp:endp), &
                    atm2lnd_inst, temperature_inst, surfalb_inst, solarabs_inst, &
                    canopystate_inst, ozone_inst, photosyns_inst, phase='sha')
            end if

            if ( use_cn .and. use_c13 ) then
               call Fractionation (bounds, fn, filterp, downreg_patch(begp:endp), &
                    atm2lnd_inst, canopystate_inst, solarabs_inst, surfalb_inst, photosyns_inst, &
                    phase='sha')
            end if

         end if ! end of if use_fates

         do f = 1, fn
            p = filterp(f)
            c = patch%column(p)
            g = patch%gridcell(p)

            ! Sensible heat conductance for air, leaf and ground
            ! Moved the original subroutine in-line...

            wta    = 1._r8/rah(p,above_canopy)     ! air
            wtl    = sa_leaf(p)/rb(p)              ! leaf
            wtg(p) = 1._r8/rah(p,below_canopy)     ! ground
            wtstem = sa_stem(p)/(rstem(p) + rb(p)) ! stem

            wtshi  = 1._r8/(wta+wtl+wtstem+wtg(p)) ! air, leaf, stem and ground

            wtl0(p) = wtl*wtshi         ! leaf
            wtg0    = wtg(p)*wtshi      ! ground
            wta0(p) = wta*wtshi         ! air

            wtstem0(p) = wtstem*wtshi              ! stem
            wtga(p)  = wta0(p)+wtg0+wtstem0(p)     ! ground + air + stem
            wtal(p) = wta0(p)+wtl0(p)+wtstem0(p)   ! air + leaf + stem

            ! internal longwave fluxes between leaf and stem
            lw_stem(p) = sa_internal(p) * emv(p) * sb * t_stem(p)**4
            lw_leaf(p) = sa_internal(p) * emv(p) * sb * t_veg(p)**4

            ! Fraction of potential evaporation from leaf

            if (fdry(p) > 0._r8) then
               rppdry  = fdry(p)*rb(p)*(laisun(p)/(rb(p)+rssun(p)) + laisha(p)/(rb(p)+rssha(p)))/elai(p)
            else
               rppdry = 0._r8
            end if
            
            ! Calculate canopy conductance for methane / oxygen (e.g. stomatal conductance & leaf bdy cond)
            if (use_lch4) then
               canopy_cond(p) = (laisun(p)/(rb(p)+rssun(p)) + laisha(p)/(rb(p)+rssha(p)))/max(elai(p), 0.01_r8)
            end if

            efpot = forc_rho(c)*((elai(p)+esai(p))/rb(p))*(qsatl(p)-qaf(p))
            h2ocan = liqcan(p) + snocan(p)

            ! When the hydraulic stress parameterization is active calculate rpp
            ! but not transpiration
            if ( use_hydrstress ) then
              if (efpot > 0._r8) then
                 if (btran(p) > btran0) then
                   rpp = rppdry + fwet(p)
                 else
                   rpp = fwet(p)
                 end if
                 !Check total evapotranspiration from leaves
                 rpp = min(rpp, (qflx_tran_veg(p)+h2ocan/dtime)/efpot)
              else
                 rpp = 1._r8
              end if
            else
              if (efpot > 0._r8) then
                 if (btran(p) > btran0) then
                    qflx_tran_veg(p) = efpot*rppdry
                    rpp = rppdry + fwet(p)
                 else
                    !No transpiration if btran below 1.e-10
                    rpp = fwet(p)
                    qflx_tran_veg(p) = 0._r8
                 end if
                 !Check total evapotranspiration from leaves
                 rpp = min(rpp, (qflx_tran_veg(p)+h2ocan/dtime)/efpot)
              else
                 !No transpiration if potential evaporation less than zero
                 rpp = 1._r8
                 qflx_tran_veg(p) = 0._r8
              end if
            end if

            ! Update conductances for changes in rpp
            ! Latent heat conductances for ground and leaf.
            ! Air has same conductance for both sensible and latent heat.
            ! Moved the original subroutine in-line...

            wtaq    = frac_veg_nosno(p)/raw(p,above_canopy)             ! air
            wtlq    = frac_veg_nosno(p)*(elai(p)+esai(p))/rb(p) * rpp   ! leaf

            !Litter layer resistance. Added by K.Sakaguchi
            snow_depth_c = params_inst%z_dl ! critical depth for 100% litter burial by snow (=litter thickness)
            fsno_dl = snow_depth(c)/snow_depth_c    ! effective snow cover for (dry)plant litter
            elai_dl = params_inst%lai_dl * (1._r8 - min(fsno_dl,1._r8)) ! exposed (dry)litter area index
            rdl = ( 1._r8 - exp(-elai_dl) ) / ( 0.004_r8*uaf(p)) ! dry litter layer resistance

            ! add litter resistance and Lee and Pielke 1992 beta
            if (delq(p) < 0._r8) then  !dew. Do not apply beta for negative flux (follow old rsoil)
               wtgq(p) = frac_veg_nosno(p)/(raw(p,below_canopy)+rdl)
            else
               if (do_soilevap_beta()) then
                  wtgq(p) = soilbeta(c)*frac_veg_nosno(p)/(raw(p,below_canopy)+rdl)
               endif
               if (do_soil_resistance_sl14()) then
                  wtgq(p) = frac_veg_nosno(p)/(raw(p,below_canopy)+soilresis(c))
               endif
            end if

            wtsqi   = 1._r8/(wtaq+wtlq+wtgq(p))

            wtgq0    = wtgq(p)*wtsqi      ! ground
            wtlq0(p) = wtlq*wtsqi         ! leaf
            wtaq0(p) = wtaq*wtsqi         ! air

            wtgaq    = wtaq0(p)+wtgq0     ! air + ground
            wtalq(p) = wtaq0(p)+wtlq0(p)  ! air + leaf

            dc1 = forc_rho(c)*cpair*wtl
            dc2 = hvap*forc_rho(c)*wtlq

            efsh = dc1*(wtga(p)*t_veg(p)-wtg0*t_grnd(c)-wta0(p)*thm(p)-wtstem0(p)*t_stem(p))
            eflx_sh_stem(p) = forc_rho(c)*cpair*wtstem*((wta0(p)+wtg0+wtl0(p))*t_stem(p)-wtg0*t_grnd(c)-wta0(p)*thm(p)-wtl0(p)*t_veg(p))
            efe(p) = dc2*(wtgaq*qsatl(p)-wtgq0*qg(c)-wtaq0(p)*forc_q(c))

            ! Evaporation flux from foliage

            erre = 0._r8
            if (efe(p)*efeb(p) < 0._r8) then
               efeold = efe(p)
               efe(p)  = 0.1_r8*efeold
               erre = efe(p) - efeold
            end if
            ! fractionate ground emitted longwave
            lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
                 +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
                 +frac_h2osfc(c)*t_h2osfc(c)**4)

            dt_veg(p) = ((1._r8-frac_rad_abs_by_stem(p))*(sabv(p) + air(p) &
                  + bir(p)*t_veg(p)**4 + cir(p)*lw_grnd) &
                  - efsh - efe(p) - lw_leaf(p) + lw_stem(p) &
                  - (cp_leaf(p)/dtime)*(t_veg(p) - tl_ini(p))) &
                  / ((1._r8-frac_rad_abs_by_stem(p))*(- 4._r8*bir(p)*t_veg(p)**3) &
                  + 4._r8*sa_internal(p)*emv(p)*sb*t_veg(p)**3 &
                  +dc1*wtga(p) +dc2*wtgaq*qsatldT(p)+ cp_leaf(p)/dtime)

            t_veg(p) = tlbef(p) + dt_veg(p)

            dels = dt_veg(p)
            del(p)  = abs(dels)
            err(p) = 0._r8
            if (del(p) > delmax) then
               dt_veg(p) = delmax*dels/del(p)
               t_veg(p) = tlbef(p) + dt_veg(p)
               err(p) = (1._r8-frac_rad_abs_by_stem(p))*(sabv(p) + air(p) &
                    + bir(p)*tlbef(p)**3*(tlbef(p) + &
                    4._r8*dt_veg(p)) + cir(p)*lw_grnd) &
                    -sa_internal(p)*emv(p)*sb*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p)) &
                    + lw_stem(p) &
                    - (efsh + dc1*wtga(p)*dt_veg(p)) - (efe(p) + &
                    dc2*wtgaq*qsatldT(p)*dt_veg(p)) &
                    - (cp_leaf(p)/dtime)*(t_veg(p) - tl_ini(p))
            end if

            ! Fluxes from leaves to canopy space
            ! "efe" was limited as its sign changes frequently.  This limit may
            ! result in an imbalance in "hvap*qflx_evap_veg" and
            ! "efe + dc2*wtgaq*qsatdt_veg"

            efpot = forc_rho(c)*((elai(p)+esai(p))/rb(p)) &
                 *(wtgaq*(qsatl(p)+qsatldT(p)*dt_veg(p)) &
                 -wtgq0*qg(c)-wtaq0(p)*forc_q(c))
            qflx_evap_veg(p) = rpp*efpot

            ! Calculation of evaporative potentials (efpot) and
            ! interception losses; flux in kg m**-2 s-1.  ecidif
            ! holds the excess energy if all intercepted water is evaporated
            ! during the timestep.  This energy is later added to the
            ! sensible heat flux.

            ! Note that when the hydraulic stress parameterization is active we don't 
            ! adjust transpiration for the new values of potential evaporation and rppdry
            ! as calculated above because transpiration would then no longer be consistent 
            ! with the vertical transpiration sink terms that are passed to Compute_VertTranSink_PHS, 
            ! thereby causing a water balance error. However, because this adjustment occurs
            ! within the leaf temperature iteration, this ends up being a small inconsistency.
            if ( use_hydrstress ) then
               ecidif = max(0._r8, qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan/dtime)
               qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan/dtime)
            else
               ecidif = 0._r8
               if (efpot > 0._r8 .and. btran(p) > btran0) then
                  qflx_tran_veg(p) = efpot*rppdry
               else
                  qflx_tran_veg(p) = 0._r8
               end if
               ecidif = max(0._r8, qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan/dtime)
               qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan/dtime)
            end if

            ! The energy loss due to above two limits is added to
            ! the sensible heat flux.

            eflx_sh_veg(p) = efsh + dc1*wtga(p)*dt_veg(p) + err(p) + erre + hvap*ecidif

           !  Update SH and lw_leaf for changes in t_veg
            eflx_sh_stem(p) = eflx_sh_stem(p) + forc_rho(c)*cpair*wtstem*(-wtl0(p)*dt_veg(p))
            lw_leaf(p) = sa_internal(p)*emv(p)*sb*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p))

            ! Re-calculate saturated vapor pressure, specific humidity, and their
            ! derivatives at the leaf surface

            call QSat(t_veg(p), forc_pbot(c), qsatl(p), &
                 es = el(p), &
                 qsdT = qsatldT(p))

            ! Update vegetation/ground surface temperature, canopy air
            ! temperature, canopy vapor pressure, aerodynamic temperature, and
            ! Monin-Obukhov stability parameter for next iteration.

            taf(p) = wtg0*t_grnd(c) + wta0(p)*thm(p) + wtl0(p)*t_veg(p) + wtstem0(p)*t_stem(p)
            qaf(p) = wtlq0(p)*qsatl(p) + wtgq0*qg(c) + forc_q(c)*wtaq0(p)

            ! Update Monin-Obukhov length and wind speed including the
            ! stability effect

            dth(p) = thm(p)-taf(p)
            dqh(p) = forc_q(c)-qaf(p)
            delq(p) = wtalq(p)*qg(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(c)

            tstar = temp1(p)*dth(p)
            qstar = temp2(p)*dqh(p)

            thvstar = tstar*(1._r8+0.61_r8*forc_q(c)) + 0.61_r8*forc_th(c)*qstar
            zeta(p) = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

            if (zeta(p) >= 0._r8) then     !stable
               zeta(p) = min(zetamax,max(zeta(p),0.01_r8))
               um(p) = max(ur(p),0.1_r8)
            else                     !unstable
               zeta(p) = max(-100._r8,min(zeta(p),-0.01_r8))
               if ( ustar(p)*thvstar > 0.0d00 )then
                  write(iulog,*) 'ustar*thvstar is positive and has to be negative'
                  write(iulog,*) 'p = ', p
                  write(iulog,*) '-grav*ustar(p)*thvstar*zii/thv(c) = ', -grav*ustar(p)*thvstar*zii/thv(c)
                  write(iulog,*) 'ustar = ', ustar(p)
                  write(iulog,*) 'thvstar = ', thvstar
                  write(iulog,*) 'thv = ', thv(c)
                  write(iulog,*) 'displa= ', displa(p)
                  write(iulog,*) 'z0mg= ', z0mg(c)
                  write(iulog,*) 'zeta= ', zeta(p)
                  write(iulog,*) 'temp1= ', temp1(p)
                  write(iulog,*) 'dth= ', dth(p)
                  write(iulog,*) 'rah(above)= ', rah(p,above_canopy)
                  write(iulog,*) 'rah(below)= ', rah(p,below_canopy)
                  !call endrun(decomp_index=p, clmlevel=namep, msg=errmsg(sourcefile, __LINE__))
                  wc = 0.0_r8
               else
                  wc = beta*(-grav*ustar(p)*thvstar*zii/thv(c))**0.333_r8
               end if
               um(p) = sqrt(ur(p)*ur(p)+wc*wc)
            end if
            obu(p) = zldis(p)/zeta(p)

            if (obuold(p)*obu(p) < 0._r8) nmozsgn(p) = nmozsgn(p)+1
            if (nmozsgn(p) >= 4) obu(p) = zldis(p)/(-0.01_r8)
            obuold(p) = obu(p)

         end do   ! end of filtered patch loop

         ! Test for convergence

         itlef = itlef+1
         if (itlef > itmin) then
            do f = 1, fn
               p = filterp(f)
               dele(p) = abs(efe(p)-efeb(p))
               efeb(p) = efe(p)
               det(p)  = max(del(p),del2(p))
               num_iter(p) = itlef
            end do
            fnold = fn
            fn = 0
            do f = 1, fnold
               p = filterp(f)
               if (.not. (det(p) < dtmin .and. dele(p) < dlemin)) then
                  fn = fn + 1
                  filterp(fn) = p
               end if
            end do
         end if
      end do ITERATION     ! End stability iteration
      call t_stopf('can_iter')

      fn = fnorig
      filterp(1:fn) = fporig(1:fn)

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)
         g = patch%gridcell(p)

         ! Energy balance check in canopy

         lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
              +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
              +frac_h2osfc(c)*t_h2osfc(c)**4)

         err(p) = (1.0_r8-frac_rad_abs_by_stem(p))*(sabv(p) + air(p) + bir(p)*tlbef(p)**3 &
              *(tlbef(p) + 4._r8*dt_veg(p)) + cir(p)*lw_grnd) &
                - lw_leaf(p) + lw_stem(p) - eflx_sh_veg(p) - hvap*qflx_evap_veg(p) &
                - ((t_veg(p)-tl_ini(p))*cp_leaf(p)/dtime)

         !  Update stem temperature; adjust outgoing longwave
         !  does not account for changes in SH or internal LW,
         !  as that would change result for t_veg above
         if (use_biomass_heat_storage) then
            if (stem_biomass(p) > 0._r8) then
               dt_stem(p) = (frac_rad_abs_by_stem(p)*(sabv(p) + air(p) + bir(p)*ts_ini(p)**4 &
                    + cir(p)*lw_grnd) - eflx_sh_stem(p) &
                    + lw_leaf(p)- lw_stem(p))/(cp_stem(p)/dtime &
                    - frac_rad_abs_by_stem(p)*bir(p)*4._r8*ts_ini(p)**3)
            else
               dt_stem(p) = 0._r8
            endif


            dhsdt_canopy(p) = dt_stem(p)*cp_stem(p)/dtime &
                 + (t_veg(p)-tl_ini(p))*cp_leaf(p)/dtime

            t_stem(p) =  t_stem(p) + dt_stem(p)
         else
            dt_stem(p) = 0._r8
         endif

         delt    = wtal(p)*t_grnd(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)-wtstem0(p)*t_stem(p)

         ! Fluxes from ground to canopy space

         taux(p) = -forc_rho(c)*forc_u(g)/ram1(p)
         tauy(p) = -forc_rho(c)*forc_v(g)/ram1(p)
         eflx_sh_grnd(p) = cpair*forc_rho(c)*wtg(p)*delt

         ! compute individual sensible heat fluxes
         delt_snow = wtal(p)*t_soisno(c,snl(c)+1)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)-wtstem0(p)*t_stem(p)
         delt_soil  = wtal(p)*t_soisno(c,1)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)-wtstem0(p)*t_stem(p)
         delt_h2osfc  = wtal(p)*t_h2osfc(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)-wtstem0(p)*t_stem(p)

         eflx_sh_snow(p) = cpair*forc_rho(c)*wtg(p)*delt_snow
         eflx_sh_soil(p) = cpair*forc_rho(c)*wtg(p)*delt_soil
         eflx_sh_h2osfc(p) = cpair*forc_rho(c)*wtg(p)*delt_h2osfc
         qflx_evap_soi(p) = forc_rho(c)*wtgq(p)*delq(p)

         ! compute individual latent heat fluxes
         delq_snow = wtalq(p)*qg_snow(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(c)
         qflx_ev_snow(p) = forc_rho(c)*wtgq(p)*delq_snow

         delq_soil = wtalq(p)*qg_soil(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(c)
         qflx_ev_soil(p) = forc_rho(c)*wtgq(p)*delq_soil

         delq_h2osfc = wtalq(p)*qg_h2osfc(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(c)
         qflx_ev_h2osfc(p) = forc_rho(c)*wtgq(p)*delq_h2osfc

         ! 2 m height air temperature

         t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))
         t_ref2m_r(p) = t_ref2m(p)

         ! 2 m height specific humidity

         q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

         ! 2 m height relative humidity

         call QSat(t_ref2m(p), forc_pbot(c), qsat_ref2m, &
              es = e_ref2m)
         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
         rh_ref2m_r(p) = rh_ref2m(p)

         ! 2m vapor pressure deficit
         vpd_ref2m(p) = e_ref2m*(1._r8-rh_ref2m(p)/100._r8)

         ! Human Heat Stress
         if ( all_human_stress_indices .or. fast_human_stress_indices ) then
            call KtoC(t_ref2m(p), tc_ref2m(p))
            call VaporPres(rh_ref2m(p), e_ref2m, vap_ref2m(p))
            call Wet_BulbS(tc_ref2m(p),rh_ref2m(p), wbt_ref2m(p))
            call HeatIndex(tc_ref2m(p), rh_ref2m(p), nws_hi_ref2m(p))
            call AppTemp(tc_ref2m(p), vap_ref2m(p), u10_clm(p), appar_temp_ref2m(p))
            call swbgt(tc_ref2m(p), vap_ref2m(p), swbgt_ref2m(p))
            call hmdex(tc_ref2m(p), vap_ref2m(p), humidex_ref2m(p))
            call dis_coiS(tc_ref2m(p), rh_ref2m(p), wbt_ref2m(p), discomf_index_ref2mS(p))
            if ( all_human_stress_indices ) then
               call Wet_Bulb(t_ref2m(p), vap_ref2m(p), forc_pbot(c), rh_ref2m(p), q_ref2m(p), &
                             teq_ref2m(p), ept_ref2m(p), wb_ref2m(p))
               call dis_coi(tc_ref2m(p), wb_ref2m(p), discomf_index_ref2m(p))
               call THIndex(tc_ref2m(p), wb_ref2m(p), thic_ref2m(p), thip_ref2m(p))
               call SwampCoolEff(tc_ref2m(p), wb_ref2m(p), swmp80_ref2m(p), swmp65_ref2m(p))
            end if
            wbt_ref2m_r(p)            = wbt_ref2m(p)
            nws_hi_ref2m_r(p)         = nws_hi_ref2m(p)
            appar_temp_ref2m_r(p)     = appar_temp_ref2m(p)
            swbgt_ref2m_r(p)          = swbgt_ref2m(p)
            humidex_ref2m_r(p)        = humidex_ref2m(p)
            discomf_index_ref2mS_r(p) = discomf_index_ref2mS(p)
            if ( all_human_stress_indices ) then
               teq_ref2m_r(p)            = teq_ref2m(p)
               ept_ref2m_r(p)            = ept_ref2m(p)
               wb_ref2m_r(p)             = wb_ref2m(p)
               discomf_index_ref2m_r(p)  = discomf_index_ref2m(p)
               thic_ref2m_r(p)           = thic_ref2m(p)
               thip_ref2m_r(p)           = thip_ref2m(p)
               swmp80_ref2m_r(p)         = swmp80_ref2m(p)
               swmp65_ref2m_r(p)         = swmp65_ref2m(p)
            end if

         end if

         ! Downward longwave radiation below the canopy

         dlrad(p) = (1._r8-emv(p))*emg(c)*forc_lwrad(c) &
              + emv(p)*emg(c)*sb*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p)) &
              *(1.0_r8-frac_rad_abs_by_stem(p)) &
              + emv(p)*emg(c)*sb*ts_ini(p)**3*(ts_ini(p) + 4._r8*dt_stem(p)) &
              *frac_rad_abs_by_stem(p)

         ! Upward longwave radiation above the canopy

         ulrad(p) = ((1._r8-emg(c))*(1._r8-emv(p))*(1._r8-emv(p))*forc_lwrad(c) &
              + emv(p)*(1._r8+(1._r8-emg(c))*(1._r8-emv(p)))*sb &
              *tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p))*(1._r8-frac_rad_abs_by_stem(p)) &
              + emv(p)*(1._r8+(1._r8-emg(c))*(1._r8-emv(p)))*sb &
              *ts_ini(p)**3*(ts_ini(p)+ 4._r8*dt_stem(p))*frac_rad_abs_by_stem(p) &
              + emg(c)*(1._r8-emv(p))*sb*lw_grnd)

         ! Calculate the skin temperature as a weighted sum of all the ground and vegetated fraction
         ! The weight is the so-called vegetation emissivity, but not that emv is actually an attentuation 
         ! function that goes to zero as LAI (ELAI + ESAI) go to zero.

         t_skin_patch(p)  =  emv(p)*t_veg(p)  +  (1._r8 - emv(p))*sqrt(sqrt(lw_grnd))

         ! Derivative of soil energy flux with respect to soil temperature

         cgrnds(p) = cgrnds(p) + cpair*forc_rho(c)*wtg(p)*wtal(p)
         cgrndl(p) = cgrndl(p) + forc_rho(c)*wtgq(p)*wtalq(p)*dqgdT(c)
         cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

         ! save before updating
         snocan_baseline(p) = snocan(p)

         ! Update dew accumulation (kg/m2)
         if (t_veg(p) > tfrz ) then ! above freezing, update accumulation in liqcan
            if ((qflx_evap_veg(p)-qflx_tran_veg(p))*dtime > liqcan(p)) then ! all liq evap
               ! In this case, all liqcan will evap. Take remainder from snocan
               snocan(p) = max(0._r8, &
                  snocan(p) + liqcan(p) + (qflx_tran_veg(p) - qflx_evap_veg(p)) * dtime)
            end if
            liqcan(p) = max(0._r8,liqcan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtime)

         else if (t_veg(p) <= tfrz) then ! below freezing, update accumulation in snocan
            if ((qflx_evap_veg(p)-qflx_tran_veg(p))*dtime > snocan(p)) then ! all sno evap
               ! In this case, all snocan will evap. Take remainder from liqcan
               liqcan(p)=liqcan(p)+snocan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtime
            end if
            snocan(p) = max(0._r8,snocan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtime)
         end if

      end do

      ! Remove snocan that got reduced by more than a factor of rel_epsilon
      ! snocan < rel_epsilon * snocan_baseline will be set to zero
      ! See NumericsMod for rel_epsilon value
      call truncate_small_values(fn, filterp, begp, endp, &
         snocan_baseline(begp:endp), snocan(begp:endp), &
         custom_rel_epsilon=1.e-10_r8)
      
      if ( use_fates ) then
         
         
         call clm_fates%wrap_accumulatefluxes(nc,fn,filterp(1:fn))
         call clm_fates%wrap_hydraulics_drive(bounds,nc,soilstate_inst, &
               waterstatebulk_inst,waterdiagnosticbulk_inst,waterfluxbulk_inst, &
               fn, filterp, solarabs_inst,energyflux_inst)

      else

         ! Determine total photosynthesis
         
         call PhotosynthesisTotal(fn, filterp, &
              atm2lnd_inst, canopystate_inst, photosyns_inst)
         
         ! Calculate water use efficiency
         !   does not support multi-layer canopy
         if (nlevcan == 1) then 
            do f = 1, fn
               p = filterp(f)
               c = patch%column(p)
               g = patch%gridcell(p)
               
               if ( is_near_local_noon( grc%londeg(g), deltasec=3600 ) .and. fpsn(p)>0._r8 )then
                  gs        = 1.e-6_r8*(laisun(p)*gs_mol_sun(p,iv)+laisha(p)*gs_mol_sha(p,iv))    ! 1e-6  converts umolH2O->molH2O
                  if ( gs>0._r8 ) then
                     iwue_ln(p)   = fpsn(p)/gs
                  else
                     iwue_ln(p)   = spval
                  end if
               else
                  iwue_ln(p)   = spval
               end if
            end do
         else
            call endrun(msg=' ERROR: IWUELN calculation not compatible with nlevcan>1 ' // &
                 errMsg(sourcefile, __LINE__))
         end if
         ! Calculate ozone uptake. This needs to be done after rssun and rsshade are
         ! computed by the Photosynthesis routine. The updated ozone uptake computed here
         ! will be used in the next time step to calculate ozone stress for the next time
         ! step's photosynthesis calculations.
         
         ! COMPILER_BUG(wjs, 2014-11-29, pgi 14.7) The following dummy variable assignment is
         ! needed with pgi 14.7 on yellowstone; without it, forc_pbot_downscaled_col gets
         ! resized inappropriately in the following subroutine call, due to a compiler bug.
         dummy_to_make_pgi_happy = ubound(atm2lnd_inst%forc_pbot_downscaled_col, 1)
         call ozone_inst%CalcOzoneUptake( &
              bounds, fn, filterp, &
              forc_pbot = atm2lnd_inst%forc_pbot_downscaled_col(bounds%begc:bounds%endc), &
              forc_th   = atm2lnd_inst%forc_th_downscaled_col(bounds%begc:bounds%endc), &
              rssun     = photosyns_inst%rssun_patch(bounds%begp:bounds%endp), &
              rssha     = photosyns_inst%rssha_patch(bounds%begp:bounds%endp), &
              rb        = frictionvel_inst%rb1_patch(bounds%begp:bounds%endp), &
              ram       = frictionvel_inst%ram1_patch(bounds%begp:bounds%endp), &
              tlai      = canopystate_inst%tlai_patch(bounds%begp:bounds%endp),  &
	      forc_o3   = atm2lnd_inst%forc_o3_grc(bounds%begg:bounds%endg))

         !---------------------------------------------------------
         !update Vc,max and Jmax by LUNA model
         if(use_luna)then
            call Acc24_Climate_LUNA(bounds, fn, filterp, &
                 canopystate_inst, photosyns_inst, &
                 surfalb_inst, solarabs_inst, &
                 temperature_inst)
            
            if(is_time_to_run_LUNA())then
               
               call Acc240_Climate_LUNA(bounds, fn, filterp, &
                    o2(begp:endp), &
                    co2(begp:endp), &
                    rb(begp:endp), &
                    rhaf(begp:endp),&
                    temperature_inst, & 
                    photosyns_inst, &
                    surfalb_inst, &
                    solarabs_inst, &
                    waterdiagnosticbulk_inst,&
                    frictionvel_inst) 
               
               call Update_Photosynthesis_Capacity(bounds, fn, filterp, &
                    dayl_factor(begp:endp), &
                    atm2lnd_inst, &
                    temperature_inst, & 
                    canopystate_inst, &
                    photosyns_inst, &
                    surfalb_inst, &
                    solarabs_inst, &
                    waterdiagnosticbulk_inst,&
                    frictionvel_inst, &
                    ozone_inst)
               
               call Clear24_Climate_LUNA(bounds, fn, filterp, &
                    canopystate_inst, photosyns_inst, &
                    surfalb_inst, solarabs_inst, &
                    temperature_inst)
            endif
            
         endif
      end if

      ! Filter out patches which have small energy balance errors; report others

      fnold = fn
      fn = 0
      do f = 1, fnold
         p = filterp(f)
         if (abs(err(p)) > 0.1_r8) then
            fn = fn + 1
            filterp(fn) = p
         end if

      end do

      do f = 1, fn
         p = filterp(f)
         write(iulog,*) 'energy balance in canopy ',p,', err=',err(p)
      end do

    end associate


  end subroutine CanopyFluxes

end module CanopyFluxesMod

