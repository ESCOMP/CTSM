module MLCanopyFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate  multilayer canopy fluxes
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog
  use decompMod           , only : bounds_type
  use atm2lndType         , only : atm2lnd_type
  use CanopyStateType     , only : canopystate_type
  use ColumnType          , only : col
  use EnergyFluxType      , only : energyflux_type
  use FrictionVelocityMod , only : frictionvel_type
  use GridcellType        , only : grc
  use PatchType           , only : patch
  use pftconMod           , only : pftcon
  use SoilStateType       , only : soilstate_type
  use SolarAbsorbedType   , only : solarabs_type
  use SurfaceAlbedoType   , only : surfalb_type
  use TemperatureType     , only : temperature_type
  use WaterFluxType       , only : waterflux_type
  use WaterStateType      , only : waterstate_type
  use MLCanopyFluxesType  , only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PRIVATE TYPES:
  integer, parameter :: nvar1d = 12     ! Number of single-level fluxes to accumulate over sub-time steps
  integer, parameter :: nvar2d = 4      ! Number of multi-level profile fluxes to accumulate over sub-time steps
  integer, parameter :: nvar3d = 10     ! Number of multi-level leaf fluxes to accumulate over sub-time steps
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MLCanopyFluxes              ! Compute canopy fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SubTimeStepFluxIntegration ! Integrate fluxes over model sub-time steps
  private :: CanopyFluxesSum            ! Sum leaf and soil fluxes
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine MLCanopyFluxes (bounds, num_exposedvegp, filter_exposedvegp, &
  atm2lnd_inst, canopystate_inst, soilstate_inst, temperature_inst, waterstate_inst, &
  waterflux_inst, energyflux_inst, frictionvel_inst, surfalb_inst, solarabs_inst, &
  mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute fluxes for sunlit and shaded leaves at each level
    ! and for soil surface
    !
    ! !USES:
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varcon, only : grav, pi => rpi, spval
    use clm_varorb, only : eccen, obliqr, lambm0, mvelpp
    use clm_varpar, only : ivis, inir
    use shr_orb_mod, only : shr_orb_decl, shr_orb_cosz
    use spmdMod, only : masterproc
    use MLclm_varcon, only : mmh2o, mmdry, cpd, cpw, rgas
    use MLclm_varctl, only : mlcan_to_clm, dtime_substep, ml_vert_init
    use MLclm_varpar, only : isun, isha, nlevmlcan, nleaf
    use MLCanopyNitrogenProfileMod, only : CanopyNitrogenProfile
    use MLCanopyTurbulenceMod, only : CanopyTurbulence
    use MLCanopyWaterMod, only : CanopyInterception, CanopyEvaporation
    use MLinitVerticalMod, only : initVerticalProfiles, initVerticalStructure
    use MLLeafBoundaryLayerMod, only : LeafBoundaryLayer
    use MLLeafHeatCapacityMod, only : LeafHeatCapacity
    use MLLeafPhotosynthesisMod, only : LeafPhotosynthesis
    use MLLongwaveRadiationMod, only : LongwaveRadiation
    use MLPlantHydraulicsMod, only: SoilResistance, PlantResistance
    use MLSolarRadiationMod, only: SolarRadiation
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_exposedvegp           ! Number of non-snow-covered veg points in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:)     ! CLM patch filter for non-snow-covered vegetation

    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(waterstate_type)  , intent(inout) :: waterstate_inst
    type(waterflux_type)   , intent(inout) :: waterflux_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    type(frictionvel_type) , intent(inout) :: frictionvel_inst
    type(surfalb_type)     , intent(inout) :: surfalb_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(mlcanopy_type)    , intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: num_mlcan                               ! Number of vegetated patches for multilayer canopy
    integer  :: filter_mlcan(bounds%endp-bounds%begp+1) ! Patch filter for multilayer canopy
    integer  :: fp                                      ! Filter index
    integer  :: p                                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                       ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                                       ! Gridcell index for CLM g/l/c/p hierarchy
    integer  :: ic                                      ! Aboveground layer index
    integer  :: nstep                                   ! Current model time step number
    integer  :: num_sub_steps                           ! Number of sub-time steps
    integer  :: niter                                   ! Current sub-time step
    real(r8) :: dtime                                   ! Model time step (s)
    real(r8) :: caldaym1                                ! Calendar day for zenith angle (1.000 on 0Z January 1 of current year)
    real(r8) :: declinm1                                ! Solar declination angle for zenith angle (radians)
    real(r8) :: eccf                                    ! Earth orbit eccentricity factor
    real(r8) :: coszen                                  ! Cosine solar zenith angle
    real(r8) :: lat, lon                                ! Latitude and longitude (radians)
    real(r8) :: totpai                                  ! Canopy lai+sai for error check (m2/m2)

    ! These are used to accumulate flux variables over the model sub-time step.
    ! The last dimension is the number of variables

    real(r8) :: flux_accumulator(bounds%begp:bounds%endp,nvar1d)                          ! Single-level fluxes
    real(r8) :: flux_accumulator_profile(bounds%begp:bounds%endp,1:nlevmlcan,nvar2d)      ! Multi-level profile fluxes
    real(r8) :: flux_accumulator_leaf(bounds%begp:bounds%endp,1:nlevmlcan,1:nleaf,nvar3d) ! Multi-level leaf fluxes
    !---------------------------------------------------------------------

    ! Variables used in this subroutine. See README.txt for a complete list of
    ! required CLM variables, and also required CLM files. See MLCanopyFluxesType.F90
    ! for a complete list of multilayer canopy variables.

    associate ( &
                                                                  ! *** CLM variables ***
    forc_u         => atm2lnd_inst%forc_u_grc                , &  ! INPUT: Atmospheric wind speed in east direction (m/s)
    forc_v         => atm2lnd_inst%forc_v_grc                , &  ! INPUT: Atmospheric wind speed in north direction (m/s)
    forc_pco2      => atm2lnd_inst%forc_pco2_grc             , &  ! INPUT: Atmospheric CO2 partial pressure (Pa)
    forc_po2       => atm2lnd_inst%forc_po2_grc              , &  ! INPUT: Atmospheric O2 partial pressure (Pa)
    forc_solad     => atm2lnd_inst%forc_solad_grc            , &  ! INPUT: Atmospheric direct beam radiation (W/m2)
    forc_solai     => atm2lnd_inst%forc_solai_grc            , &  ! INPUT: Atmospheric diffuse radiation (W/m2)
    forc_t         => atm2lnd_inst%forc_t_downscaled_col     , &  ! INPUT: Atmospheric temperature (K)
    forc_q         => atm2lnd_inst%forc_q_downscaled_col     , &  ! INPUT: Atmospheric specific humidity (kg/kg)
    forc_pbot      => atm2lnd_inst%forc_pbot_downscaled_col  , &  ! INPUT: Atmospheric pressure (Pa)
    forc_lwrad     => atm2lnd_inst%forc_lwrad_downscaled_col , &  ! INPUT: Atmospheric longwave radiation (W/m2)
    forc_rain      => atm2lnd_inst%forc_rain_downscaled_col  , &  ! INPUT: Rainfall rate (mm/s)
    forc_snow      => atm2lnd_inst%forc_snow_downscaled_col  , &  ! INPUT: Snowfall rate (mm/s)
    elai           => canopystate_inst%elai_patch            , &  ! INPUT: Leaf area index of canopy (m2/m2)
    esai           => canopystate_inst%esai_patch            , &  ! INPUT: Stem area index of canopy (m2/m2)
    snl            => col%snl                                , &  ! INPUT: Number of snow layers
    z              => col%z                                  , &  ! INPUT: Soil layer depth (m)
    zi             => col%zi                                 , &  ! INPUT: Soil layer depth at layer interface (m)
    btran_patch    => energyflux_inst%btran_patch            , &  ! INPUT: Transpiration wetness factor (0 to 1)
    minlwp_SPA     => pftcon%minlwp_SPA                      , &  ! INPUT: Minimum leaf water potential (MPa)
    soilresis      => soilstate_inst%soilresis_col           , &  ! INPUT: Soil evaporative resistance (s/m)
    thk            => soilstate_inst%thk_col                 , &  ! INPUT: Soil layer thermal conductivity (W/m/K)
    smp_l          => soilstate_inst%smp_l_col               , &  ! INPUT: Soil layer matric potential (mm)
    albgrd         => surfalb_inst%albgrd_col                , &  ! INPUT: Direct beam albedo of ground (soil)
    albgri         => surfalb_inst%albgri_col                , &  ! INPUT: Diffuse albedo of ground (soil)
    t_a10_patch    => temperature_inst%t_a10_patch           , &  ! INPUT: 10-day running mean of the 2-m temperature (K)
    t_soisno       => temperature_inst%t_soisno_col          , &  ! INPUT: Soil temperature (K)
    eflx_lh_tot    => energyflux_inst%eflx_lh_tot_patch      , &  ! OUTPUT: patch total latent heat flux (W/m2)
    eflx_sh_tot    => energyflux_inst%eflx_sh_tot_patch      , &  ! OUTPUT: patch total sensible heat flux (W/m2)
    eflx_lwrad_out => energyflux_inst%eflx_lwrad_out_patch   , &  ! OUTPUT: patch emitted infrared (longwave) radiation (W/m2)
    taux           => energyflux_inst%taux_patch             , &  ! OUTPUT: patch wind (shear) stress: e-w (kg/m/s2)
    tauy           => energyflux_inst%tauy_patch             , &  ! OUTPUT: patch wind (shear) stress: n-s (kg/m/s2)
    fv             => frictionvel_inst%fv_patch              , &  ! OUTPUT: patch friction velocity (m/s)
    u10_clm        => frictionvel_inst%u10_clm_patch         , &  ! OUTPUT: patch 10-m wind (m/s)
    fsa            => solarabs_inst%fsa_patch                , &  ! OUTPUT: patch solar radiation absorbed (total) (W/m2)
    albd           => surfalb_inst%albd_patch                , &  ! OUTPUT: patch surface albedo (direct)
    albi           => surfalb_inst%albi_patch                , &  ! OUTPUT: patch surface albedo (diffuse)
    t_ref2m        => temperature_inst%t_ref2m_patch         , &  ! OUTPUT: patch 2 m height surface air temperature (K)
    qflx_evap_tot  => waterflux_inst%qflx_evap_tot_patch     , &  ! OUTPUT: patch total evapotranspiration flux (kg H2O/m2/s)
    q_ref2m        => waterstate_inst%q_ref2m_patch          , &  ! OUTPUT: patch 2 m height surface specific humidity (kg/kg)

                                                                  ! *** Multilayer canopy variables ***
    zref           => mlcanopy_inst%zref_forcing             , &  ! Atmospheric reference height (m)
    tref           => mlcanopy_inst%tref_forcing             , &  ! Air temperature at reference height (K)
    thref          => mlcanopy_inst%thref_forcing            , &  ! Atmospheric potential temperature at reference height (K)
    thvref         => mlcanopy_inst%thvref_forcing           , &  ! Atmospheric virtual potential temperature at reference height (K)
    qref           => mlcanopy_inst%qref_forcing             , &  ! Specific humidity at reference height (kg/kg)
    eref           => mlcanopy_inst%eref_forcing             , &  ! Vapor pressure at reference height (Pa)
    uref           => mlcanopy_inst%uref_forcing             , &  ! Wind speed at reference height (m/s)
    pref           => mlcanopy_inst%pref_forcing             , &  ! Air pressure at reference height (Pa)
    co2ref         => mlcanopy_inst%co2ref_forcing           , &  ! Atmospheric CO2 at reference height (umol/mol)
    o2ref          => mlcanopy_inst%o2ref_forcing            , &  ! Atmospheric O2 at reference height (mmol/mol)
    rhoair         => mlcanopy_inst%rhoair_forcing           , &  ! Air density at reference height (kg/m3)
    rhomol         => mlcanopy_inst%rhomol_forcing           , &  ! Molar density at reference height (mol/m3)
    mmair          => mlcanopy_inst%mmair_forcing            , &  ! Molecular mass of air at reference height (kg/mol)
    cpair          => mlcanopy_inst%cpair_forcing            , &  ! Specific heat of air (constant pressure) at reference height (J/mol/K)
    solar_zen      => mlcanopy_inst%solar_zen_forcing        , &  ! Solar zenith angle (radians)
    swskyb         => mlcanopy_inst%swskyb_forcing           , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd         => mlcanopy_inst%swskyd_forcing           , &  ! Atmospheric diffuse solar radiation (W/m2)
    lwsky          => mlcanopy_inst%lwsky_forcing            , &  ! Atmospheric longwave radiation (W/m2)
    qflx_rain      => mlcanopy_inst%qflx_rain_forcing        , &  ! Rainfall (mm H2O/s = kg H2O/m2/s)
    qflx_snow      => mlcanopy_inst%qflx_snow_forcing        , &  ! Snowfall (mm H2O/s = kg H2O/m2/s)
    tacclim        => mlcanopy_inst%tacclim_forcing          , &  ! Average air temperature for acclimation (K)
    ncan           => mlcanopy_inst%ncan_canopy              , &  ! Number of layers
    ntop           => mlcanopy_inst%ntop_canopy              , &  ! Index for top leaf layer
    nbot           => mlcanopy_inst%nbot_canopy              , &  ! Index for bottom leaf layer
    lai            => mlcanopy_inst%lai_canopy               , &  ! Leaf area index of canopy (m2/m2)
    sai            => mlcanopy_inst%sai_canopy               , &  ! Stem area index of canopy (m2/m2)
    swveg          => mlcanopy_inst%swveg_canopy             , &  ! Absorbed solar radiation: vegetation (W/m2)
    lwup           => mlcanopy_inst%lwup_canopy              , &  ! Upward longwave radiation above canopy (W/m2)
    shflx          => mlcanopy_inst%shflx_canopy             , &  ! Total sensible heat flux, including soil (W/m2)
    lhflx          => mlcanopy_inst%lhflx_canopy             , &  ! Total latent heat flux, including soil (W/m2)
    etflx          => mlcanopy_inst%etflx_canopy             , &  ! Total water vapor flux, including soil (mol H2O/m2/s)
    ustar          => mlcanopy_inst%ustar_canopy             , &  ! Friction velocity (m/s)
    fracminlwp     => mlcanopy_inst%fracminlwp_canopy        , &  ! Fraction of canopy that is water-stressed
    btran          => mlcanopy_inst%btran_soil               , &  ! Soil wetness factor for stomatal conductance (-)
    albsoib        => mlcanopy_inst%albsoib_soil             , &  ! Direct beam albedo of ground (-)
    albsoid        => mlcanopy_inst%albsoid_soil             , &  ! Diffuse albedo of ground (-)
    swsoi          => mlcanopy_inst%swsoi_soil               , &  ! Absorbed solar radiation: ground (W/m2)
    lwsoi          => mlcanopy_inst%lwsoi_soil               , &  ! Absorbed longwave radiation: ground (W/m2)
    rnsoi          => mlcanopy_inst%rnsoi_soil               , &  ! Net radiation: ground (W/m2)
    tg             => mlcanopy_inst%tg_soil                  , &  ! Soil surface temperature (K)
    tg_bef         => mlcanopy_inst%tg_bef_soil              , &  ! Soil surface temperature for previous timestep (K)
    rhg            => mlcanopy_inst%rhg_soil                 , &  ! Relative humidity of airspace at soil surface (fraction)
    soilres        => mlcanopy_inst%soilres_soil             , &  ! Soil evaporative resistance (s/m)
    soil_t         => mlcanopy_inst%soil_t_soil              , &  ! Temperature of first snow/soil layer (K)
    soil_dz        => mlcanopy_inst%soil_dz_soil             , &  ! Depth to temperature of first snow/soil layer (m)
    soil_tk        => mlcanopy_inst%soil_tk_soil             , &  ! Thermal conductivity of first snow/soil layer (W/m/K)
    dlai           => mlcanopy_inst%dlai_profile             , &  ! Canopy layer leaf area index (m2/m2)
    dsai           => mlcanopy_inst%dsai_profile             , &  ! Canopy layer stem area index (m2/m2)
    dpai           => mlcanopy_inst%dpai_profile             , &  ! Canopy layer plant area index (m2/m2)
    sumpai         => mlcanopy_inst%sumpai_profile           , &  ! Canopy layer cumulative plant area index (m2/m2)
    dlai_frac      => mlcanopy_inst%dlai_frac_profile        , &  ! Canopy layer leaf area index (fraction of canopy total)
    dsai_frac      => mlcanopy_inst%dsai_frac_profile        , &  ! Canopy layer stem area index (fraction of canopy total)
    fracsun        => mlcanopy_inst%fracsun_profile          , &  ! Canopy layer sunlit fraction (-)
    tair           => mlcanopy_inst%tair_profile             , &  ! Canopy layer air temperature (K)
    eair           => mlcanopy_inst%eair_profile             , &  ! Canopy layer vapor pressure (Pa)
    cair           => mlcanopy_inst%cair_profile             , &  ! Canopy layer atmospheric CO2 (umol/mol)
    tair_bef       => mlcanopy_inst%tair_bef_profile         , &  ! Canopy layer air temperature for previous timestep (K)
    eair_bef       => mlcanopy_inst%eair_bef_profile         , &  ! Canopy layer vapor pressure for previous timestep (Pa)
    cair_bef       => mlcanopy_inst%cair_bef_profile         , &  ! Canopy layer atmospheric CO2 for previous timestep (umol/mol)
    lwpveg         => mlcanopy_inst%lwpveg_profile           , &  ! Canopy layer leaf water potential (MPa)
    swleaf         => mlcanopy_inst%swleaf_leaf              , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    lwleaf         => mlcanopy_inst%lwleaf_leaf              , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf         => mlcanopy_inst%rnleaf_leaf              , &  ! Leaf net radiation (W/m2 leaf)
    tleaf          => mlcanopy_inst%tleaf_leaf               , &  ! Leaf temperature (K)
    tleaf_bef      => mlcanopy_inst%tleaf_bef_leaf           , &  ! Leaf temperature for previous timestep (K)
    lwpleaf        => mlcanopy_inst%lwpleaf_leaf               &  ! Leaf water potential (MPa)
    )

    ! Get current step counter (nstep) and step size (dtime)

    nstep = get_nstep()
    dtime = get_step_size()

    ! Set number of sub-time steps for flux calculations

    num_sub_steps = int(dtime / dtime_substep)

    ! Build filter of patches to process with multilayer canopy

    num_mlcan = 0
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       g = patch%gridcell(fp)
!      if (grc%latdeg(g) .gt. -2.9_r8 .and. grc%latdeg(g) .lt. -2.7_r8) then
!      if (grc%londeg(g) .gt. 294.5_r8 .and. grc%londeg(g) .lt. 295.5_r8) then
       num_mlcan = num_mlcan + 1
       filter_mlcan(num_mlcan) = p
!      end if
!      end if
    end do

    ! Initialize canopy vertical structure and profiles. This is only done
    ! once (on first time step) because the forcing height (and therefore
    ! the volume of air in the surface layer) changes between time steps.
    ! The call for initialization is triggered by zref = spval.

    ml_vert_init = 0
    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       if (zref(p) == spval) ml_vert_init = 1
    end do

    if (ml_vert_init == 1) then
       if (masterproc) then
          write (iulog,*) 'Attempting to initialize multilayer canopy vertical structure .....'
       end if

       call initVerticalStructure (bounds, num_mlcan, filter_mlcan, &
       canopystate_inst, frictionvel_inst, mlcanopy_inst)

       call initVerticalProfiles (num_mlcan, filter_mlcan, &
       atm2lnd_inst, mlcanopy_inst)

       if (masterproc) then
          write (iulog,*) 'Successfully initialized multilayer canopy vertical structure'
       end if
    end if

    ! Copy CLM variables to multilayer canopy variables. Note the
    ! distinction between grid cell (g), column (c), and patch (p)
    ! variables. All multilayer canopy variables are for patches.

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       c = patch%column(p)
       g = patch%gridcell(p)

       ! Atmospheric forcing: CLM grid cell (g) variables -> patch (p) variables

       uref(p) = max (1.0_r8, sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       swskyb(p,ivis) = forc_solad(g,ivis) ; swskyb(p,inir) = forc_solad(g,inir)
       swskyd(p,ivis) = forc_solai(g,ivis) ; swskyd(p,inir) = forc_solai(g,inir)

       ! Atmospheric forcing: CLM column (c) variables -> patch (p) variables

       tref(p) = forc_t(c)
       qref(p) = forc_q(c)
       pref(p) = forc_pbot(c)
       lwsky(p) = forc_lwrad(c)
       qflx_rain(p) = forc_rain(c)
       qflx_snow(p) = forc_snow(c)

       ! CO2 and O2: CLM grid cell (g) -> patch (p). Note the unit conversion

       co2ref(p) = forc_pco2(g) / forc_pbot(c) * 1.e06_r8  ! Pa -> umol/mol
       o2ref(p)  = forc_po2(g) / forc_pbot(c) * 1.e03_r8   ! Pa -> mmol/mol

       ! Miscellaneous

       btran(p) = btran_patch(p)
       tacclim(p) = t_a10_patch(p)

       ! Ground (soil) albedos: CLM column (c) -> patch (p)

       albsoib(p,ivis) = albgrd(c,ivis) ; albsoib(p,inir) = albgrd(c,inir)
       albsoid(p,ivis) = albgri(c,ivis) ; albsoid(p,inir) = albgri(c,inir)

       ! Soil evaporative resistance: CLM column (c) -> patch (p)

       soilres(p) = soilresis(c)

       ! Properties of first soil layer: CLM column (c) -> patch (p)

       soil_t(p) = t_soisno(c,snl(c)+1)          ! Temperature of first snow/soil layer (K)
       soil_dz(p) = (z(c,snl(c)+1)-zi(c,snl(c))) ! Depth to temperature of first snow/soil layer (m)
       soil_tk(p) = thk(c,snl(c)+1)              ! Thermal conductivity of first snow/soil layer (W/m/K)

    end do

    ! Solar zenith angle. Need to subtract one time step (-dtime) because
    ! zenith angle is calculated for the beginning of the time step. So use
    ! calendar day at beginning of the time step (nstep-1).

    caldaym1 = get_curr_calday(offset=-int(dtime))
    call shr_orb_decl (caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf)

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       c = patch%column(p)
       g = patch%gridcell(p)
       lat = grc%latdeg(g) * pi / 180._r8
       lon = grc%londeg(g) * pi / 180._r8
       coszen = shr_orb_cosz (caldaym1, lat, lon, declinm1)
       solar_zen(p) = acos(max(0.01_r8,coszen))

       ! Compare coszen to that expected from CLM

       if (abs(coszen-surfalb_inst%coszen_col(c)) .gt. 1.e-03_r8) then
          write (iulog,*) nstep, coszen, surfalb_inst%coszen_col(c)
          call endrun (msg=' ERROR: MLCanopyFluxes: coszen error')
       end if
    end do

    ! Derived atmospheric input

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       eref(p) = qref(p) * pref(p) / (mmh2o / mmdry + (1._r8 - mmh2o / mmdry) * qref(p))
       rhomol(p) = pref(p) / (rgas * tref(p))
       rhoair(p) = rhomol(p) * mmdry * (1._r8 - (1._r8 - mmh2o/mmdry) * eref(p) / pref(p))
       mmair(p) = rhoair(p) / rhomol(p)
       cpair(p) = cpd * (1._r8 + (cpw/cpd - 1._r8) * qref(p)) * mmair(p)
       thref(p) = tref(p) + 0.0098_r8 * zref(p)
       thvref(p) = thref(p) * (1._r8 + 0.61_r8 * qref(p))
    end do

    ! Update leaf and stem area profile for current values

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)

       ! Values for current time step from CLM

       lai(p) = elai(p)
       sai(p) = esai(p)

       ! Vertical profiles

       do ic = 1, ncan(p)
          dlai(p,ic) = dlai_frac(p,ic) * lai(p)
          dsai(p,ic) = dsai_frac(p,ic) * sai(p)
          dpai(p,ic) = dlai(p,ic) + dsai(p,ic)
       end do

       totpai = sum(dpai(p,1:ncan(p)))
       if (abs(totpai - (lai(p)+sai(p))) > 1.e-06_r8) then
          call endrun (msg=' ERROR: MLCanopyFluxes: plant area index not updated correctly')
       end if

       ! Cumulative plant area index: Fill in canopy layers
       ! (at the midpoint) starting from the top

       do ic = ntop(p), 1, -1
          if (ic == ntop(p)) then
             sumpai(p,ic) = 0.5_r8 * dpai(p,ic)
          else
             sumpai(p,ic) = sumpai(p,ic+1) + 0.5_r8 * (dpai(p,ic+1) + dpai(p,ic))
          end if
       end do

       ! Layers above the canopy have no vegetation

       do ic = ntop(p)+1, ncan(p)
          sumpai(p,ic) = 0._r8
       end do

    end do

    ! Solar radiation transfer through the canopy

    call SolarRadiation (bounds, num_mlcan, filter_mlcan, mlcanopy_inst)

    ! Plant hydraulics

    call SoilResistance (num_mlcan, filter_mlcan, &
    soilstate_inst, waterstate_inst, mlcanopy_inst)

    call PlantResistance (num_mlcan, filter_mlcan, mlcanopy_inst)

    ! Canopy profile of photosynthetic capacity

    call CanopyNitrogenProfile (num_mlcan, filter_mlcan, mlcanopy_inst)

    ! Use sub-stepping to calculate fluxes over the full time step

    do niter = 1, num_sub_steps

       ! Save values for previous timestep

       do fp = 1, num_mlcan
          p = filter_mlcan(fp)
          tg_bef(p) = tg(p)
          do ic = 1, ncan(p)
             tleaf_bef(p,ic,isun) = tleaf(p,ic,isun)
             tleaf_bef(p,ic,isha) = tleaf(p,ic,isha)
             tair_bef(p,ic) = tair(p,ic)
             eair_bef(p,ic) = eair(p,ic)
             cair_bef(p,ic) = cair(p,ic)
          end do
       end do

       ! Canopy interception

       call CanopyInterception (num_mlcan, filter_mlcan, mlcanopy_inst)

       ! Longwave radiation transfer through canopy

       call LongwaveRadiation (bounds, num_mlcan, filter_mlcan, mlcanopy_inst)

       ! Net radiation at each layer and at ground

       do fp = 1, num_mlcan
          p = filter_mlcan(fp)
          do ic = 1, ncan(p)
             rnleaf(p,ic,isun) = swleaf(p,ic,isun,ivis) + swleaf(p,ic,isun,inir) + lwleaf(p,ic,isun)
             rnleaf(p,ic,isha) = swleaf(p,ic,isha,ivis) + swleaf(p,ic,isha,inir) + lwleaf(p,ic,isha)
          end do
          rnsoi(p) = swsoi(p,ivis) + swsoi(p,inir) + lwsoi(p)
       end do

       ! Leaf heat capacity

       call LeafHeatCapacity (num_mlcan, filter_mlcan, mlcanopy_inst)

       ! Leaf boundary layer conductance

       call LeafBoundaryLayer (num_mlcan, filter_mlcan, isun, mlcanopy_inst)
       call LeafBoundaryLayer (num_mlcan, filter_mlcan, isha, mlcanopy_inst)

       ! Photosynthesis and stomatal conductance

       call LeafPhotosynthesis (num_mlcan, filter_mlcan, isun, mlcanopy_inst)
       call LeafPhotosynthesis (num_mlcan, filter_mlcan, isha, mlcanopy_inst)

       ! Relative humidity in soil airspace

       do fp = 1, num_mlcan
          p = filter_mlcan(fp)
          c = patch%column(p)
          rhg(p) = exp(grav * mmh2o * smp_l(c,1)*1.e-03_r8 / (rgas * t_soisno(c,1)))
       end do

       ! Canopy turbulence, scalar source/sink fluxes for leaves and soil, and scalar
       ! profiles using above- and within-canopy coupling with the roughness sublayer
       ! (RSL) parameterization

       call CanopyTurbulence (niter, num_mlcan, filter_mlcan, mlcanopy_inst)

       ! Update canopy intercepted water for evaporation and dew

       call CanopyEvaporation (num_mlcan, filter_mlcan, mlcanopy_inst)

       ! Fluxes need to be accumulated over all sub-time steps. Other
       ! variables are instantaneous for the final sub-time step.

       call SubTimeStepFluxIntegration (niter, num_sub_steps, num_mlcan, filter_mlcan, &
       flux_accumulator, flux_accumulator_profile, flux_accumulator_leaf, mlcanopy_inst)

    end do    ! End sub-stepping loop

    ! Sum leaf and soil fluxes over canopy

    call CanopyFluxesSum (num_mlcan, filter_mlcan, mlcanopy_inst)

    ! Need to merge temperature and leaf water potential for sunlit and
    ! shaded leaves because sun/shade fractions change over time

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       do ic = nbot(p), ntop(p)
          tleaf(p,ic,isun) = tleaf(p,ic,isun) * fracsun(p,ic) + tleaf(p,ic,isha) * (1._r8 - fracsun(p,ic))
          tleaf(p,ic,isha) = tleaf(p,ic,isun)
          lwpveg(p,ic) = lwpleaf(p,ic,isun) * fracsun(p,ic) + lwpleaf(p,ic,isha) * (1._r8 - fracsun(p,ic))
       end do
    end do

    ! Diagnose fraction of the canopy that is water stressed

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       fracminlwp(p) = 0._r8

       do ic = nbot(p), ntop(p)
          if (lwpveg(p,ic) <= minlwp_SPA(patch%itype(p))) then
             fracminlwp(p) = fracminlwp(p) + dpai(p,ic)
          end if
       end do

       if ((lai(p) + sai(p)) > 0._r8) then
          fracminlwp(p) = fracminlwp(p) / (lai(p) + sai(p))
       end if
    end do

    ! Copy multilayer canopy variables to CLM variables. These are
    ! passed from CLM to CAM. The variables passed to CAM are found
    ! in the CLM routine: main/lnd2atmType.F90. These are CLM grid
    ! cell variables. The mapping from CLM patch to CLM gridcell is
    ! done in: main/lnd2atmMod.F90

    if (mlcan_to_clm == 1) then
       do fp = 1, num_mlcan
          p = filter_mlcan(fp)
          albd(p,ivis) = 0._r8 ; albd(p,inir) = 0._r8
          albi(p,ivis) = 0._r8 ; albi(p,inir) = 0._r8
          taux(p) = 0._r8
          tauy(p) = 0._r8
          eflx_lh_tot(p) = lhflx(p)
          eflx_sh_tot(p) = shflx(p)
          eflx_lwrad_out(p) = lwup(p)
          qflx_evap_tot(p) = etflx(p) * mmh2o
          fv(p) = ustar(p)
          u10_clm(p) = 0._r8
          t_ref2m(p) = 0._r8
          q_ref2m(p) = 0._r8
          fsa(p) = swveg(p,ivis) + swveg(p,inir) + swsoi(p,ivis) + swsoi(p,inir)
       end do
    end if

    end associate
  end subroutine MLCanopyFluxes

  !-----------------------------------------------------------------------
  subroutine SubTimeStepFluxIntegration (niter, num_sub_steps, num_filter, filter, &
  flux_accumulator, flux_accumulator_profile, flux_accumulator_leaf, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Integrate fluxes over model sub-time steps
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: niter                                ! Current sub-time step
    integer, intent(in) :: num_sub_steps                        ! Number of sub-time steps
    integer, intent(in) :: num_filter                           ! Number of patches in filter
    integer, intent(in) :: filter(:)                            ! Patch filter
    real(r8), intent(inout) :: flux_accumulator(:,:)            ! Single-level flux accumulator variable
    real(r8), intent(inout) :: flux_accumulator_profile(:,:,:)  ! Multi-level profile flux accumulator variable
    real(r8), intent(inout) :: flux_accumulator_leaf(:,:,:,:)   ! Multi-level leaf flux accumulator variable
    type(mlcanopy_type), intent(in) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                              ! Filter index
    integer  :: p                                               ! Patch index for CLM g/l/c/p hierarchy
    integer  :: i,j,k                                           ! Variable index
    !---------------------------------------------------------------------

    associate ( &
    ustar       => mlcanopy_inst%ustar_canopy     , &  ! Friction velocity (m/s)
    lwup        => mlcanopy_inst%lwup_canopy      , &  ! Upward longwave radiation above canopy (W/m2)
    lwsoi       => mlcanopy_inst%lwsoi_soil       , &  ! Absorbed longwave radiation: ground (W/m2)
    rnsoi       => mlcanopy_inst%rnsoi_soil       , &  ! Net radiation: ground (W/m2)
    shsoi       => mlcanopy_inst%shsoi_soil       , &  ! Sensible heat flux: ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi_soil       , &  ! Latent heat flux: ground (W/m2)
    etsoi       => mlcanopy_inst%etsoi_soil       , &  ! Water vapor flux: ground (mol H2O/m2/s)
    gsoi        => mlcanopy_inst%gsoi_soil        , &  ! Soil heat flux (W/m2)
    gac0        => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    qflx_intr   => mlcanopy_inst%qflx_intr_canopy , &  ! Intercepted precipitation (kg H2O/m2/s)
    qflx_tflrain => mlcanopy_inst%qflx_tflrain_canopy , &  ! Total rain throughfall onto ground (kg H2O/m2/s)
    qflx_tflsnow => mlcanopy_inst%qflx_tflsnow_canopy , &  ! Total snow throughfall onto ground (kg H2O/m2/s)
    shair       => mlcanopy_inst%shair_profile    , &  ! Canopy layer air sensible heat flux (W/m2)
    etair       => mlcanopy_inst%etair_profile    , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair       => mlcanopy_inst%stair_profile    , &  ! Canopy layer air storage heat flux (W/m2)
    gac         => mlcanopy_inst%gac_profile      , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    lwleaf      => mlcanopy_inst%lwleaf_leaf      , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf      => mlcanopy_inst%rnleaf_leaf      , &  ! Leaf net radiation (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf_leaf      , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf_leaf      , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf_leaf      , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf_leaf      , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    stleaf      => mlcanopy_inst%stleaf_leaf      , &  ! Leaf storage heat flux (W/m2 leaf)
    an          => mlcanopy_inst%an_leaf          , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    ag          => mlcanopy_inst%ag_leaf          , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    gs          => mlcanopy_inst%gs_leaf            &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Initialize flux variables that are summed over sub-time steps

       if (niter == 1) then
          flux_accumulator(p,:) = 0._r8
          flux_accumulator_profile(p,:,:) = 0._r8
          flux_accumulator_leaf(p,:,:,:) = 0._r8
       end if

       ! Accumulate fluxes over sub-time steps

       i = 0
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + ustar(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lwup(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lwsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + rnsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + shsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lhsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + etsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + gsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + gac0(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_intr(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_tflrain(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_tflsnow(p)

       j = 0
       j = j + 1; flux_accumulator_profile(p,:,j) = flux_accumulator_profile(p,:,j) + shair(p,:)
       j = j + 1; flux_accumulator_profile(p,:,j) = flux_accumulator_profile(p,:,j) + etair(p,:)
       j = j + 1; flux_accumulator_profile(p,:,j) = flux_accumulator_profile(p,:,j) + stair(p,:)
       j = j + 1; flux_accumulator_profile(p,:,j) = flux_accumulator_profile(p,:,j) + gac(p,:)

       k = 0
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + lwleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + rnleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + shleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + lhleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + trleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + evleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + stleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + an(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + ag(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + gs(p,:,:)

       if (i > nvar1d .or. j > nvar2d .or. k > nvar3d) then
          call endrun (msg=' ERROR: SubTimeStepFluxIntegration: nvar error')
       end if

       ! Average fluxes over sub-timesteps

       if (niter == num_sub_steps) then

          ! Time averaging

          flux_accumulator(p,:) = flux_accumulator(p,:) / float(num_sub_steps)
          flux_accumulator_profile(p,:,:) = flux_accumulator_profile(p,:,:) / float(num_sub_steps)
          flux_accumulator_leaf(p,:,:,:) = flux_accumulator_leaf(p,:,:,:) / float(num_sub_steps)

          ! Map fluxes to variables: variables must be in the same order as above

          i = 0
          i = i + 1; ustar(p) = flux_accumulator(p,i)
          i = i + 1; lwup(p) = flux_accumulator(p,i)
          i = i + 1; lwsoi(p) = flux_accumulator(p,i)
          i = i + 1; rnsoi(p) = flux_accumulator(p,i)
          i = i + 1; shsoi(p) = flux_accumulator(p,i)
          i = i + 1; lhsoi(p) = flux_accumulator(p,i)
          i = i + 1; etsoi(p) = flux_accumulator(p,i)
          i = i + 1; gsoi(p) = flux_accumulator(p,i)
          i = i + 1; gac0(p) = flux_accumulator(p,i)
          i = i + 1; qflx_intr(p) = flux_accumulator(p,i)
          i = i + 1; qflx_tflrain(p) = flux_accumulator(p,i)
          i = i + 1; qflx_tflsnow(p) = flux_accumulator(p,i)

          j = 0
          j = j + 1; shair(p,:) = flux_accumulator_profile(p,:,j)
          j = j + 1; etair(p,:) = flux_accumulator_profile(p,:,j)
          j = j + 1; stair(p,:) = flux_accumulator_profile(p,:,j)
          j = j + 1; gac(p,:) = flux_accumulator_profile(p,:,j)

          k = 0
          k = k + 1; lwleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; rnleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; shleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; lhleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; trleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; evleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; stleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; an(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; ag(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; gs(p,:,:) = flux_accumulator_leaf(p,:,:,k)

          if (i > nvar1d .or. j > nvar2d .or. k > nvar3d) then
             call endrun (msg=' ERROR: SubTimeStepFluxIntegration: nvar error')
          end if

       end if

    end do

    end associate
  end subroutine SubTimeStepFluxIntegration

  !-----------------------------------------------------------------------
  subroutine CanopyFluxesSum (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Sum leaf and soil fluxes to get canopy fluxes
    !
    ! !USES:
    use clm_varpar, only : ivis, inir
    use MLclm_varctl, only : turb_type
    use MLclm_varpar, only : isun, isha
    use MLWaterVaporMod, only : LatVap
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                      ! Aboveground layer index
    real(r8) :: err                     ! Energy imbalance (W/m2)
    real(r8) :: radin                   ! Incoming radiation (W/m2)
    real(r8) :: radout                  ! Outgoing radiation (W/m2)
    real(r8) :: avail                   ! Available energy (W/m2)
    real(r8) :: flux                    ! Turbulent fluxes + storage (W/m2)
    real(r8) :: fracgreen               ! Green fraction of plant area index: lai/(lai+sai)
    !---------------------------------------------------------------------

    associate ( &
    tref        => mlcanopy_inst%tref_forcing     , &  ! Air temperature at reference height (K)
    swskyb      => mlcanopy_inst%swskyb_forcing   , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd_forcing   , &  ! Atmospheric diffuse solar radiation (W/m2)
    lwsky       => mlcanopy_inst%lwsky_forcing    , &  ! Atmospheric longwave radiation (W/m2)
    ncan        => mlcanopy_inst%ncan_canopy      , &  ! Number of layers
    ntop        => mlcanopy_inst%ntop_canopy      , &  ! Index for top leaf layer
    swveg       => mlcanopy_inst%swveg_canopy     , &  ! Absorbed solar radiation: vegetation (W/m2)
    albcan      => mlcanopy_inst%albcan_canopy    , &  ! Albedo above canopy (-)
    lwup        => mlcanopy_inst%lwup_canopy      , &  ! Upward longwave radiation above canopy (W/m2)
    shsoi       => mlcanopy_inst%shsoi_soil       , &  ! Sensible heat flux: ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi_soil       , &  ! Latent heat flux: ground (W/m2)
    gsoi        => mlcanopy_inst%gsoi_soil        , &  ! Soil heat flux (W/m2)
    swsoi       => mlcanopy_inst%swsoi_soil       , &  ! Absorbed solar radiation: ground (W/m2)
    lwsoi       => mlcanopy_inst%lwsoi_soil       , &  ! Absorbed longwave radiation: ground (W/m2)
    etsoi       => mlcanopy_inst%etsoi_soil       , &  ! Water vapor flux: ground (mol H2O/m2/s)
    dpai        => mlcanopy_inst%dpai_profile     , &  ! Canopy layer plant area index (m2/m2)
    fwet        => mlcanopy_inst%fwet_profile     , &  ! Canopy layer fraction of plant area index that is wet
    fdry        => mlcanopy_inst%fdry_profile     , &  ! Canopy layer fraction of plant area index that is green and dry
    shair       => mlcanopy_inst%shair_profile    , &  ! Canopy layer air sensible heat flux (W/m2)
    etair       => mlcanopy_inst%etair_profile    , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair       => mlcanopy_inst%stair_profile    , &  ! Canopy layer air storage heat flux (W/m2)
    fracsun     => mlcanopy_inst%fracsun_profile  , &  ! Canopy layer sunlit fraction (-)
    lwleaf      => mlcanopy_inst%lwleaf_leaf      , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf      => mlcanopy_inst%rnleaf_leaf      , &  ! Leaf net radiation (W/m2 leaf)
    stleaf      => mlcanopy_inst%stleaf_leaf      , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf_leaf      , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf_leaf      , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf_leaf      , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf_leaf      , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    swleaf      => mlcanopy_inst%swleaf_leaf      , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    an          => mlcanopy_inst%an_leaf          , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    ag          => mlcanopy_inst%ag_leaf          , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
                                                       ! *** Output ***
    rnet        => mlcanopy_inst%rnet_canopy      , &  ! Total net radiation, including soil (W/m2)
    stflx       => mlcanopy_inst%stflx_canopy     , &  ! Canopy storage heat flux (W/m2)
    shflx       => mlcanopy_inst%shflx_canopy     , &  ! Total sensible heat flux, including soil (W/m2)
    lhflx       => mlcanopy_inst%lhflx_canopy     , &  ! Total latent heat flux, including soil (W/m2)
    etflx       => mlcanopy_inst%etflx_canopy     , &  ! Total water vapor flux, including soil (mol H2O/m2/s)
    lwveg       => mlcanopy_inst%lwveg_canopy     , &  ! Absorbed longwave radiation: vegetation (W/m2)
    lwvegsun    => mlcanopy_inst%lwvegsun_canopy  , &  ! Absorbed longwave radiation: sunlit canopy (W/m2)
    lwvegsha    => mlcanopy_inst%lwvegsha_canopy  , &  ! Absorbed longwave radiation: shaded canopy (W/m2)
    shveg       => mlcanopy_inst%shveg_canopy     , &  ! Sensible heat flux: vegetation (W/m2)
    shvegsun    => mlcanopy_inst%shvegsun_canopy  , &  ! Sensible heat flux: sunlit canopy (W/m2)
    shvegsha    => mlcanopy_inst%shvegsha_canopy  , &  ! Sensible heat flux: shaded canopy (W/m2)
    lhveg       => mlcanopy_inst%lhveg_canopy     , &  ! Latent heat flux: vegetation (W/m2)
    lhvegsun    => mlcanopy_inst%lhvegsun_canopy  , &  ! Latent heat flux: sunlit canopy (W/m2)
    lhvegsha    => mlcanopy_inst%lhvegsha_canopy  , &  ! Latent heat flux: shaded canopy (W/m2)
    etveg       => mlcanopy_inst%etveg_canopy     , &  ! Water vapor flux: vegetation (mol H2O/m2/s)
    etvegsun    => mlcanopy_inst%etvegsun_canopy  , &  ! Water vapor flux: sunlit canopy (mol H2O/m2/s)
    etvegsha    => mlcanopy_inst%etvegsha_canopy  , &  ! Water vapor flux: shaded canopy (mol H2O/m2/s)
    gppveg      => mlcanopy_inst%gppveg_canopy    , &  ! Gross primary production: vegetation (umol CO2/m2/s)
    gppvegsun   => mlcanopy_inst%gppvegsun_canopy , &  ! Gross primary production: sunlit canopy (umol CO2/m2/s)
    gppvegsha   => mlcanopy_inst%gppvegsha_canopy , &  ! Gross primary production: shaded canopy (umol CO2/m2/s)
    swsrc       => mlcanopy_inst%swsrc_profile    , &  ! Canopy layer source/sink flux: absorbed solar radiation (W/m2)
    lwsrc       => mlcanopy_inst%lwsrc_profile    , &  ! Canopy layer source/sink flux: absorbed longwave radiation (W/m2)
    rnsrc       => mlcanopy_inst%rnsrc_profile    , &  ! Canopy layer source/sink flux: net radiation (W/m2)
    stsrc       => mlcanopy_inst%stsrc_profile    , &  ! Canopy layer source/sink flux: storage heat flux (W/m2)
    shsrc       => mlcanopy_inst%shsrc_profile    , &  ! Canopy layer source/sink flux: sensible heat (W/m2)
    lhsrc       => mlcanopy_inst%lhsrc_profile    , &  ! Canopy layer source/sink flux: latent heat (W/m2)
    etsrc       => mlcanopy_inst%etsrc_profile    , &  ! Canopy layer source/sink flux: water vapor (mol H2O/m2/s)
    fco2src     => mlcanopy_inst%fco2src_profile    &  ! Canopy layer source/sink flux: CO2 (umol CO2/m2/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Leaf flux profiles

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             lwsrc(p,ic) = (lwleaf(p,ic,isun)*fracsun(p,ic) + lwleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic)
             swsrc(p,ic,ivis) = (swleaf(p,ic,isun,ivis)*fracsun(p,ic) + swleaf(p,ic,isha,ivis)*(1._r8 - fracsun(p,ic))) * dpai(p,ic)
             swsrc(p,ic,inir) = (swleaf(p,ic,isun,inir)*fracsun(p,ic) + swleaf(p,ic,isha,inir)*(1._r8 - fracsun(p,ic))) * dpai(p,ic)
             rnsrc(p,ic) = (rnleaf(p,ic,isun)*fracsun(p,ic) + rnleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic)
             stsrc(p,ic) = (stleaf(p,ic,isun)*fracsun(p,ic) + stleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic)
             shsrc(p,ic) = (shleaf(p,ic,isun)*fracsun(p,ic) + shleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic)
             lhsrc(p,ic) = (lhleaf(p,ic,isun)*fracsun(p,ic) + lhleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic)
             etsrc(p,ic) = (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic) &
                         + (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             fco2src(p,ic) = (an(p,ic,isun)*fracsun(p,ic) + an(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic) * fracgreen
          else
             lwsrc(p,ic) = 0._r8
             swsrc(p,ic,ivis) = 0._r8
             swsrc(p,ic,inir) = 0._r8
             rnsrc(p,ic) = 0._r8
             stsrc(p,ic) = 0._r8
             shsrc(p,ic) = 0._r8
             lhsrc(p,ic) = 0._r8
             etsrc(p,ic) = 0._r8
             fco2src(p,ic) = 0._r8
          end if
       end do

       ! Sum leaf fluxes

       lwveg(p) = 0._r8
       stflx(p) = 0._r8
       shveg(p) = 0._r8
       lhveg(p) = 0._r8
       etveg(p) = 0._r8
       gppveg(p) = 0._r8

       do ic = 1, ncan(p)
          lwveg(p) = lwveg(p) + lwsrc(p,ic)
          stflx(p) = stflx(p) + stsrc(p,ic)
          shveg(p) = shveg(p) + shsrc(p,ic)
          lhveg(p) = lhveg(p) + lhsrc(p,ic)
          etveg(p) = etveg(p) + etsrc(p,ic)
          if (dpai(p,ic) > 0._r8) then
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             gppveg(p) = gppveg(p) + (ag(p,ic,isun)*fracsun(p,ic) + ag(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic) * fracgreen
          end if
       end do

       ! Check energy balance for conservation

       err = swveg(p,ivis) + swveg(p,inir) + lwveg(p) - shveg(p) - lhveg(p) - stflx(p)
       if (abs(err) >= 1.e-03_r8) then
          call endrun (msg=' ERROR: CanopyFluxesSum: energy conservation error (1)')
       end if

       ! Turbulent fluxes

       select case (turb_type)
       case (0)
          ! Sum of layer fluxes
          shflx(p) = shveg(p) + shsoi(p)
          etflx(p) = etveg(p) + etsoi(p)
          lhflx(p) = lhveg(p) + lhsoi(p)
       case (1)
          ! Turbulent fluxes are at the top of the canopy
          shflx(p) = shair(p,ntop(p))
          etflx(p) = etair(p,ntop(p))
          lhflx(p) = etair(p,ntop(p)) * LatVap(tref(p))
       case default
          call endrun (msg=' ERROR: CanopyFluxesSum: turbulence type not valid')
       end select

       ! Add canopy air heat storage to storage flux

       do ic = 1, ntop(p)
          stflx(p) = stflx(p) + stair(p,ic)
       end do

       ! Overall energy balance check:
       ! radiation in - radiation out - soil heat = available energy = turbulent flux + canopy storage flux

       rnet(p) = swveg(p,ivis) + swveg(p,inir) + swsoi(p,ivis) + swsoi(p,inir) + lwveg(p) + lwsoi(p)
       radin = swskyb(p,ivis) + swskyd(p,ivis) + swskyb(p,inir) + swskyd(p,inir) + lwsky(p)
       radout = albcan(p,ivis)*(swskyb(p,ivis)+swskyd(p,ivis)) + albcan(p,inir)*(swskyb(p,inir)+swskyd(p,inir)) + lwup(p)

       err = rnet(p) - (radin - radout)
       if (abs(err) > 0.01_r8) then
          call endrun (msg=' ERROR: CanopyFluxesSum: energy conservation error (2)')
       end if

       avail = radin - radout - gsoi(p)
       flux = shflx(p) + lhflx(p) + stflx(p)
       err = avail - flux
       if (abs(err) > 0.01_r8) then
          call endrun (msg=' ERROR: CanopyFluxesSum: energy conservation error (3)')
       end if

       ! Sunlit and shaded canopy fluxes

       lwvegsun(p) = 0._r8 ; lwvegsha(p) = 0._r8
       shvegsun(p) = 0._r8 ; shvegsha(p) = 0._r8
       lhvegsun(p) = 0._r8 ; lhvegsha(p) = 0._r8
       etvegsun(p) = 0._r8 ; etvegsha(p) = 0._r8
       gppvegsun(p) = 0._r8 ; gppvegsha(p) = 0._r8

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             lwvegsun(p) = lwvegsun(p) + lwleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             lwvegsha(p) = lwvegsha(p) + lwleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             shvegsun(p) = shvegsun(p) + shleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             shvegsha(p) = shvegsha(p) + shleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             lhvegsun(p) = lhvegsun(p) + lhleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             lhvegsha(p) = lhvegsha(p) + lhleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             etvegsun(p) = etvegsun(p) + (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic)
             etvegsha(p) = etvegsha(p) + (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             gppvegsun(p) = gppvegsun(p) + ag(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) * fracgreen
             gppvegsha(p) = gppvegsha(p) + ag(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) * fracgreen
          end if
       end do

    end do

    end associate
  end subroutine CanopyFluxesSum

end module MLCanopyFluxesMod
