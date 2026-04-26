module NVPLayerDynamicsMod

  ! ---------------------------------------------------------------------------
  ! DESCRIPTION:
  ! [PORTED by Hui Tang: NVP (non-vascular plant, i.e. moss/lichen) layer dynamics]
  !
  ! Updates CLM's vertical column geometry for the NVP layer at index 0,
  ! based on FATES-derived coverage (col%frac_nvp) and thickness (col%dz_nvp).
  !
  ! Also provides subroutines for NVP water physics:
  !   NVPWaterRetentionCurve  — van Genuchten water potential [mm]
  !   NVPHydraulicConductivity — Mualem-van Genuchten hydraulic conductivity [m/s]
  !   NVPEvaporation           — NVP surface evaporation flux [kg m-2 s-1]
  !
  ! Called from clmfates_interfaceMod%wrap_update_hlmfates_dyn each FATES
  ! dynamics timestep, after bc_out%nvp_dz_pa and bc_out%nvp_frac_pa have
  ! been aggregated to the column.
  !
  ! Layer design:
  !   - NVP occupies vertical index 0 (the slot below the bottom-most snow layer).
  !   - col%jbot_sno = -1 when NVP is active (snow loops stop at index -1).
  !   - col%jbot_sno = 0  when NVP is inactive (standard CLM snow behaviour).
  !   - col%dz(c,0) and col%z(c,0) are set from col%dz_nvp(c).
  !   - col%zi(c,0) = 0 (soil surface) is unchanged.
  !   - col%zi(c,-1) = -dz_nvp (top interface of NVP layer).
  !
  ! Energy and mass conservation on layer activation/deactivation:
  !   - Activation  (inactive→active): initialise t_soisno(c,0) = t_soisno(c,1),
  !     h2osoi_liq(c,0) = 0, h2osoi_ice(c,0) = 0.  The NVP layer starts
  !     isothermal with the surface soil and dry; FATES will hydrate it over
  !     subsequent timesteps.
  !   - Deactivation (active→inactive): residual water in layer 0 is merged into
  !     layer 1 with energy-weighted temperature; layer 0 state is then zeroed.
  !   - Resize (thickness change while active): temperature and h2osoi are
  !     unchanged — both are intensive quantities (T in K; h2osoi in kg m-2
  !     of ground area independent of layer thickness).
  ! ---------------------------------------------------------------------------

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use ColumnType             , only : col
  use TemperatureType        , only : temperature_type
  use WaterStateType         , only : waterstate_type
  use WaterFluxBulkType      , only : waterfluxbulk_type
  use WaterDiagnosticBulkType, only : waterdiagnosticbulk_type
  ! [PORTED by Hui Tang: soilstate for bidirectional NVP-soil Darcy flux]
  use SoilStateType          , only : soilstate_type
  use clm_varcon             , only : cpliq, cpice, denh2o, roverg, tfrz, denice
  ! [PORTED by Hui Tang: use_nvp_undersnow flag to deactivate NVP when snow present]
  use clm_varctl             , only : use_nvp_undersnow
  use QSatMod                , only : QSat
  ! [PORTED by Hui Tang: runtime-tunable NVP physics parameters]
  use NVPParamsMod

  implicit none
  private

  public :: UpdateNVPLayer
  public :: NVPWaterRetentionCurve
  public :: NVPHydraulicConductivity
  public :: NVPEvaporation
  public :: NVPWaterBalance_Column
  public :: NVPLayerRestart

  character(len=*), parameter, private :: sourcefile = __FILE__

contains

  ! ===========================================================================

  subroutine UpdateNVPLayer(c, temperature_inst, waterstate_inst)
    ! -------------------------------------------------------------------------
    ! Update NVP layer presence and geometry for column c, maintaining energy
    ! and water conservation across activation / deactivation transitions.
    !
    ! Inputs (from col):
    !   col%dz_nvp(c)   — FATES-derived column-effective NVP thickness [m]
    !   col%frac_nvp(c) — FATES-derived column NVP fractional coverage [0-1]
    !
    ! Outputs (written to col, and optionally temperature_inst, waterstate_inst):
    !   col%nvp_layer_active(c)                   — .true. when NVP layer is active
    !   col%jbot_sno(c)                            — -1 (NVP active) or 0 (inactive)
    !   col%dz(c,0)                                — NVP layer thickness [m]
    !   col%z(c,0)                                 — NVP layer centre depth [m]
    !   col%zi(c,-1)                               — NVP layer top interface [m]
    !
    ! Conservation (applied only when both optional args are present):
    !   temperature_inst%t_soisno_col(c,0)         — initialised on activation
    !   waterstate_inst%h2osoi_liq_col(c,0/1)      — conserved on deactivation
    !   waterstate_inst%h2osoi_ice_col(c,0/1)      — conserved on deactivation
    !
    ! Pass temperature_inst and waterstate_inst during normal timestep calls.
    ! Omit them during restart / cold-start where thermodynamic state is set
    ! independently (restart file or initTemperatureMod).
    ! -------------------------------------------------------------------------
    integer,                             intent(in)    :: c               ! column index
    type(temperature_type),   optional,  intent(inout) :: temperature_inst
    class(waterstate_type),   optional,  intent(inout) :: waterstate_inst

    real(r8) :: dz_nvp
    logical  :: was_active   ! NVP layer state at entry (before this update)
    logical  :: now_active   ! NVP layer state after geometry check
    real(r8) :: cv0          ! heat capacity of NVP layer  [J m-2 K-1]
    real(r8) :: cv1          ! heat capacity of soil layer 1 before merge [J m-2 K-1]
    real(r8) :: cv1_new      ! heat capacity of soil layer 1 after merge  [J m-2 K-1]

    dz_nvp     = col%dz_nvp(c)
    was_active = col%nvp_layer_active(c)

    if (col%frac_nvp(c) > nvp_frac_min .and. dz_nvp > 0._r8) then

       ! --- Active (Appear or Grow/shrink) ---
       now_active              = .true.
       col%nvp_layer_active(c) = .true.
       col%jbot_sno(c)         = -1
       col%dz(c,0)             = dz_nvp
       ! Layer centre is half the thickness above the soil surface (zi(c,0)=0)
       col%z(c,0)              = -0.5_r8 * dz_nvp
       ! Top interface of NVP layer = bottom of the snow layer above
       col%zi(c,-1)            = -dz_nvp

    else

       ! --- Inactive (Disappear or Absent) ---
       now_active              = .false.
       col%nvp_layer_active(c) = .false.
       col%jbot_sno(c)         = 0
       col%dz(c,0)             = 0._r8
       col%z(c,0)              = 0._r8
       ! Restore zi(c,-1) to soil surface when NVP is absent
       col%zi(c,-1)            = 0._r8

    end if

    ! [PORTED by Hui Tang: deactivate NVP under snow when use_nvp_undersnow=.false.]
    ! When snow is present (snl < 0), NVP cannot exchange energy or moisture with the
    ! atmosphere directly. Override the FATES-based activation to treat the column as
    ! standard CLM (jbot_sno=0) so snow layers couple directly to soil.
    if (.not. use_nvp_undersnow .and. col%snl(c) < 0 .and. now_active) then
       now_active              = .false.
       col%nvp_layer_active(c) = .false.
       col%jbot_sno(c)         = 0
       col%dz(c,0)             = 0._r8
       col%z(c,0)              = 0._r8
       col%zi(c,-1)            = 0._r8
    end if

    ! [PORTED by Hui Tang: energy and mass conservation on NVP layer state transitions]
    ! Conservation is only applied when temperature_inst and waterstate_inst are
    ! provided (normal timestep path).  During restart / cold-start the thermo-
    ! dynamic state is restored separately and conservation is skipped.
    if (present(temperature_inst) .and. present(waterstate_inst)) then

       if (.not. was_active .and. now_active) then
          ! --- Appear: initialise layer-0 thermodynamic state from soil layer 1 ---
          ! Temperature inherits from layer 1 to avoid spurious heat flux.
          ! h2osoi starts at zero; FATES provides moisture via its own hydrology.
          temperature_inst%t_soisno_col(c,0)  = temperature_inst%t_soisno_col(c,1)
          waterstate_inst%h2osoi_liq_col(c,0) = 0._r8
          waterstate_inst%h2osoi_ice_col(c,0) = 0._r8

       else if (was_active .and. .not. now_active) then
          ! --- Disappear: merge layer 0 into layer 1 conserving energy and water ---
          ! Compute heat capacities (cv = cpliq*liq + cpice*ice, units J m-2 K-1)
          cv0 = cpliq * waterstate_inst%h2osoi_liq_col(c,0) + &
                cpice * waterstate_inst%h2osoi_ice_col(c,0)
          cv1 = cpliq * waterstate_inst%h2osoi_liq_col(c,1) + &
                cpice * waterstate_inst%h2osoi_ice_col(c,1)
          ! Transfer water mass to layer 1
          waterstate_inst%h2osoi_liq_col(c,1) = waterstate_inst%h2osoi_liq_col(c,1) + &
                                                 waterstate_inst%h2osoi_liq_col(c,0)
          waterstate_inst%h2osoi_ice_col(c,1) = waterstate_inst%h2osoi_ice_col(c,1) + &
                                                 waterstate_inst%h2osoi_ice_col(c,0)
          ! Energy-weighted temperature for layer 1 after merge
          cv1_new = cpliq * waterstate_inst%h2osoi_liq_col(c,1) + &
                    cpice * waterstate_inst%h2osoi_ice_col(c,1)
          if (cv1_new > 0._r8) then
             temperature_inst%t_soisno_col(c,1) = &
                  (cv0 * temperature_inst%t_soisno_col(c,0) + &
                   cv1 * temperature_inst%t_soisno_col(c,1)) / cv1_new
          end if
          ! Zero layer 0; set T to layer-1 value to avoid stale data on reactivation
          waterstate_inst%h2osoi_liq_col(c,0) = 0._r8
          waterstate_inst%h2osoi_ice_col(c,0) = 0._r8
          temperature_inst%t_soisno_col(c,0)  = temperature_inst%t_soisno_col(c,1)

       end if
       ! --- Grow/shrink (active→active): T and h2osoi are per unit ground area,
       !     so no adjustment is needed when dz_nvp changes. ---

    end if

  end subroutine UpdateNVPLayer

  ! ===========================================================================

  subroutine NVPWaterRetentionCurve(vol_liq, eff_porosity, n_van, alpha_van, watsat, watres, smp)
    ! -------------------------------------------------------------------------
    ! Convert NVP volumetric liquid water content to soil water potential [mm]
    ! using the van Genuchten (1980) retention curve formulation.
    !
    ! Ported from Python moss_water_code.py::water_retention_curve.
    !
    ! Arguments:
    !   vol_liq   — volumetric liquid water content [m3 m-3]
    !   eff_porosity — effective porosity [m3 m-3]
    !   n_van     — van Genuchten shape parameter n [-] (> 1)
    !   alpha_van — van Genuchten inverse air-entry pressure [cm-1 * 10]
    !   watsat    — saturated volumetric water content (= theta_nvp_max) [m3 m-3]
    !   watres    — residual volumetric water content [m3 m-3]
    !   smp       — soil/NVP matric potential [mm]  (negative; more negative = drier)
    ! -------------------------------------------------------------------------
    real(r8), intent(in)  :: vol_liq
    real(r8), intent(in)  :: eff_porosity
    real(r8), intent(in)  :: n_van
    real(r8), intent(in)  :: alpha_van
    real(r8), intent(in)  :: watsat
    real(r8), intent(in)  :: watres
    real(r8), intent(out) :: smp

    real(r8) :: m_van         ! van Genuchten m = 1 - 1/n
    real(r8) :: satfrac       ! effective saturation fraction [-]

    m_van        = 1.0_r8 - 1.0_r8 / n_van

    satfrac = (vol_liq - watres) / (eff_porosity - watres)
    satfrac = max(0.0_r8, min(1.0_r8, satfrac))   ! clamp to [0, 1]

    ! van Genuchten retention curve: psi in units of (1/alpha_van)
    smp = -(1.0_r8 / alpha_van) * &
          (satfrac**( 1.0_r8 / (-m_van) ) - 1.0_r8)**(1.0_r8 / n_van)
    smp = smp * 102.0_r8   ! convert to mm

  end subroutine NVPWaterRetentionCurve

  ! ===========================================================================

  subroutine NVPHydraulicConductivity(vol_liq, eff_porosity, n_van, watsat, watres, ksat, khydr)
    ! -------------------------------------------------------------------------
    ! Compute NVP hydraulic conductivity [m s-1] using the Mualem-van Genuchten
    ! (1976/1980) formulation.
    !
    ! Ported from Python moss_water_code.py::hydraulic_conductivity.
    !
    ! Arguments:
    !   vol_liq   — volumetric liquid water content [m3 m-3]
    !   eff_porosity — effective porosity [m3 m-3]
    !   n_van     — van Genuchten shape parameter n [-]
    !   watsat    — saturated volumetric water content [m3 m-3]
    !   watres    — residual volumetric water content [m3 m-3]
    !   ksat      — saturated hydraulic conductivity [m s-1]
    !   khydr     — unsaturated hydraulic conductivity [m s-1]
    ! -------------------------------------------------------------------------
    real(r8), intent(in)  :: vol_liq
    real(r8), intent(in)  :: eff_porosity
    real(r8), intent(in)  :: n_van
    real(r8), intent(in)  :: watsat
    real(r8), intent(in)  :: watres
    real(r8), intent(in)  :: ksat
    real(r8), intent(out) :: khydr

    real(r8) :: m_van         ! van Genuchten m = 1 - 1/n
    real(r8) :: satfrac       ! effective saturation fraction [-]

    m_van        = 1.0_r8 - 1.0_r8 / n_van

    satfrac = (vol_liq - watres) / (eff_porosity - watres)
    satfrac = max(0.0_r8, min(1.0_r8, satfrac))   ! clamp to [0, 1]

    ! Mualem-van Genuchten: khydr = ksat * Se^0.5 * (1 - (1 - Se^(1/m))^m)^2
    khydr = ksat &
          * satfrac**0.5_r8 &
          * (1.0_r8 - (1.0_r8 - satfrac**(1.0_r8 / m_van))**m_van)**2.0_r8

  end subroutine NVPHydraulicConductivity

  ! ===========================================================================

  subroutine NVPEvaporation(theta_nvp, t_nvp, forc_pbot, rho_atm, q_atm, raw, &
       n_van, alpha_van, watsat, watres, &
       evap_nvp, rnvp, psi_nvp, alpha_nvp, q_nvp)
    ! -------------------------------------------------------------------------
    ! Compute NVP surface evaporation flux [kg m-2 s-1].
    !
    ! Ported from the hourly loop in Python moss_water_code.py.
    !
    ! Physics:
    !   1. Surface resistance rnvp increases cubically as NVP dries:
    !        rnvp = rnvp_min + rnvp_amp * (1 - satfrac)^rnvp_exp
    !   2. Water potential psi_nvp [mm] from van Genuchten retention curve.
    !   3. Activity correction alpha_nvp [-] from the Kelvin equation:
    !        alpha_nvp = exp(psi_nvp / (roverg * t_nvp))
    !      where roverg = R_w/g * 1000 [mm K] is the CLM constant from clm_varcon.
    !   4. Specific humidity at NVP surface: q_nvp = alpha_nvp * qs(t_nvp, pbot)
    !   5. Evaporation: E = max(-rho_atm * (q_atm - q_nvp) / (raw + rnvp), 0)
    !      (positive = evaporation from NVP; condensation onto NVP is excluded)
    !
    ! Arguments:
    !   theta_nvp  — volumetric liquid water content of NVP [m3 m-3]
    !   t_nvp      — NVP surface temperature [K]
    !   forc_pbot  — atmospheric pressure [Pa]
    !   rho_atm    — air density [kg m-3]
    !   q_atm      — specific humidity of overlying air [kg kg-1]
    !   raw        — aerodynamic resistance to water vapour transfer [s m-1]
    !   n_van      — van Genuchten n [-]
    !   alpha_van  — van Genuchten alpha [cm-1 * 10]
    !   watsat     — saturated volumetric water content [m3 m-3]
    !   watres     — residual volumetric water content [m3 m-3]
    !   evap_nvp   — NVP evaporation flux [kg m-2 s-1]  (out, >= 0)
    !   rnvp       — NVP surface resistance [s m-1]     (out, diagnostic)
    !   psi_nvp    — NVP matric potential [mm]           (out, diagnostic)
    !   alpha_nvp  — Kelvin activity correction [-]      (out, diagnostic)
    !   q_nvp      — specific humidity at NVP surface [kg kg-1] (out, diagnostic)
    ! -------------------------------------------------------------------------
    real(r8), intent(in)  :: theta_nvp
    real(r8), intent(in)  :: t_nvp
    real(r8), intent(in)  :: forc_pbot
    real(r8), intent(in)  :: rho_atm
    real(r8), intent(in)  :: q_atm
    real(r8), intent(in)  :: raw
    real(r8), intent(in)  :: n_van
    real(r8), intent(in)  :: alpha_van
    real(r8), intent(in)  :: watsat
    real(r8), intent(in)  :: watres
    real(r8), intent(out) :: evap_nvp
    real(r8), intent(out) :: rnvp
    real(r8), intent(out) :: psi_nvp
    real(r8), intent(out) :: alpha_nvp
    real(r8), intent(out) :: q_nvp

    real(r8) :: eff_porosity  ! effective porosity [m3 m-3]
    real(r8) :: satfrac       ! effective saturation fraction [-]
    real(r8) :: qs_nvp        ! saturation specific humidity at t_nvp [kg kg-1]

    ! --- 1. Effective saturation fraction ---
    eff_porosity = max(0.01_r8, watsat)
    satfrac      = (theta_nvp - watres) / (eff_porosity - watres)
    satfrac      = max(0.0_r8, min(1.0_r8, satfrac))

    ! --- 2. Surface resistance (high when dry, low when saturated) ---
    rnvp = rnvp_min + rnvp_amp * (1.0_r8 - satfrac)**rnvp_exp

    ! --- 3. Van Genuchten matric potential [mm] ---
    call NVPWaterRetentionCurve(theta_nvp, eff_porosity, n_van, alpha_van, watsat, watres, psi_nvp)

    ! --- 4. Kelvin activity correction: alpha = exp(psi / (roverg * T)) ---
    ! roverg = R_w/g * 1000 [mm K] so that psi [mm] / (roverg [mm K] * T [K]) is
    ! dimensionless.
    alpha_nvp = exp(psi_nvp / (roverg * t_nvp))

    ! --- 5. Specific humidity at NVP surface ---
    call QSat(t_nvp, forc_pbot, qs_nvp)
    q_nvp = alpha_nvp * qs_nvp

    ! --- 6. Evaporation flux (no condensation: capped at zero) ---
    ! Convention: positive evap_nvp means water leaves NVP surface.
    ! Flux = -rho * (q_atm - q_nvp) / (raw + rnvp)
    !      = rho * (q_nvp - q_atm) / (raw + rnvp)
    evap_nvp = max(-rho_atm * (q_atm - q_nvp) / (raw + rnvp), 0.0_r8)

  end subroutine NVPEvaporation

  ! ===========================================================================

  ! [PORTED by Hui Tang: NVP column water balance — gravity drainage from layer 0 to soil layer 1]
  ! [PORTED by Hui Tang: bidirectional Darcy flux between NVP layer 0 and soil layer 1]
  subroutine NVPWaterBalance_Column(bounds, dtime, waterfluxbulk_inst, waterstate_inst, &
       waterdiagnosticbulk_inst, soilstate_inst, temperature_inst)
    ! -------------------------------------------------------------------------
    ! Update h2osoi_liq(c,0) for the NVP layer and compute the net water
    ! exchange with soil layer 1 via a bidirectional Darcy flux.
    !
    ! Physics:
    !   q01 = K_interface * ( (ψ_nvp - ψ_soil1) / Δz_mm + 1 )
    !
    !   where q01 > 0 means downward (NVP drains to soil) and
    !         q01 < 0 means upward   (NVP absorbs from soil).
    !
    !   ψ_nvp   — van Genuchten matric potential of NVP layer [mm]
    !   ψ_soil1 — CLM Clapp-Hornberger matric potential of soil layer 1 [mm],
    !             from previous timestep (soilstate_inst%smp_l_col)
    !   Δz_mm   — distance between layer centres [mm]
    !           = (0.5*dz(c,0) + 0.5*dz(c,1)) * 1000
    !   K_interface — harmonic mean of K_nvp and K_soil1 [mm/s]
    !
    ! Capping:
    !   Downward (q01 > 0): limited by available NVP liquid water
    !   Upward   (q01 < 0): limited by current soil layer-1 liquid water
    !
    ! qflx_nvp_drain_col is signed (+: NVP→soil, -: soil→NVP) and is added
    ! to qflx_infl in Infiltration so SoilWater sees the net exchange.
    !
    ! Must be called after:
    !   - SnowWater (so qflx_rain_plus_snomelt is finalised)
    !   - clm_drv_patch2col p2c (so qflx_ev_nvp_col is valid)
    ! and before:
    !   - SetQflxInputs / Infiltration (so qflx_nvp_drain_col is available)
    !
    ! Outputs written:
    !   waterfluxbulk_inst%qflx_nvp_infl_col(c)  [mm/s]   water into NVP top
    !   waterfluxbulk_inst%qflx_nvp_drain_col(c) [mm/s]   net NVP-soil exchange (+down)
    !   waterstate_inst%h2osoi_liq_col(c,0)       [kg m-2] updated NVP water store
    !   waterdiagnosticbulk_inst%fwet_nvp_col(c)  [-]      updated NVP wet fraction
    ! -------------------------------------------------------------------------
    type(bounds_type),              intent(in)    :: bounds
    real(r8),                       intent(in)    :: dtime                ! timestep [s]
    type(waterfluxbulk_type),       intent(inout) :: waterfluxbulk_inst
    class(waterstate_type),         intent(inout) :: waterstate_inst
    type(waterdiagnosticbulk_type), intent(inout) :: waterdiagnosticbulk_inst
    type(soilstate_type),           intent(in)    :: soilstate_inst
    ! [PORTED by Hui Tang: temperature_inst needed for ice/liquid partitioning of NVP layer]
    type(temperature_type),         intent(in)    :: temperature_inst

    integer  :: c
    real(r8) :: frac_h2osfc    ! fractional area with surface water [-]
    real(r8) :: frac_nvp_eff   ! effective NVP area fraction (not covered by h2osfc) [-]
    real(r8) :: vol_liq        ! NVP volumetric liquid water content [m3 m-3]
    ! [PORTED by Hui Tang: locals for ice partitioning of NVP layer water]
    real(r8) :: vol_ice        ! NVP volumetric ice content [m3 m-3]
    real(r8) :: eff_porosity   ! effective porosity (watsat - vol_ice) [m3 m-3]
    real(r8) :: khydr_nvp      ! NVP unsaturated hydraulic conductivity [m s-1]
    real(r8) :: K_nvp_mms      ! khydr_nvp converted to mm/s
    real(r8) :: K_soil1        ! soil layer 1 hydraulic conductivity [mm/s]
    real(r8) :: K_interface    ! harmonic-mean interface conductivity [mm/s]
    real(r8) :: psi_nvp        ! NVP van Genuchten matric potential [mm]
    real(r8) :: smp_soil1      ! soil layer 1 matric potential, prev. timestep [mm]
    real(r8) :: dz_iface_mm    ! distance between NVP and soil layer 1 centres [mm]
    real(r8) :: q01            ! Darcy flux NVP→soil (+down, -up) [mm s-1]
    real(r8) :: h2osoi_net     ! h2osoi_liq(c,0) after infl and evap [kg m-2]
    real(r8) :: satfrac        ! NVP effective saturation fraction [-]

    associate( &
         qflx_rain_plus_snomelt => waterfluxbulk_inst%qflx_rain_plus_snomelt_col, &
         qflx_ev_nvp_col        => waterfluxbulk_inst%qflx_ev_nvp_col,            &
         qflx_nvp_infl_col      => waterfluxbulk_inst%qflx_nvp_infl_col,          &
         qflx_nvp_drain_col     => waterfluxbulk_inst%qflx_nvp_drain_col,         &
         h2osoi_liq             => waterstate_inst%h2osoi_liq_col,                &
         h2osoi_ice             => waterstate_inst%h2osoi_ice_col,                &  ! [PORTED by Hui Tang: ice content for porosity reduction]
         h2onvp_col             => waterstate_inst%h2onvp_col,                    &  ! [PORTED by Hui Tang: sync diagnostic copy]
         smp_l                  => soilstate_inst%smp_l_col,                      &
         hk_l                   => soilstate_inst%hk_l_col,                       &
         frac_h2osfc_col        => waterdiagnosticbulk_inst%frac_h2osfc_col,      &
         fwet_nvp_col           => waterdiagnosticbulk_inst%fwet_nvp_col,         &
         vwc_nvp_col            => waterdiagnosticbulk_inst%vwc_nvp_col,          &  ! [PORTED by Hui Tang: volumetric water content]
         t_nvp_col              => temperature_inst%t_nvp_col                     &  ! Input:  [real(r8) (:)   ] NVP (moss/lichen) temperature (Kelvin)
         )

      do c = bounds%begc, bounds%endc

        if (.not. col%nvp_layer_active(c)) then
           qflx_nvp_infl_col(c)  = 0._r8
           qflx_nvp_drain_col(c) = 0._r8
           cycle
        end if

        ! --- Effective NVP area: not covered by surface water ---
        frac_h2osfc  = frac_h2osfc_col(c)
        frac_nvp_eff = min(col%frac_nvp(c), max(0._r8, 1._r8 - frac_h2osfc))

        ! --- Water input to NVP from precipitation / snowmelt ---
        qflx_nvp_infl_col(c) = frac_nvp_eff * qflx_rain_plus_snomelt(c)   ! [mm/s]

        ! --- NVP volumetric water content (clamped to valid range) ---
        if (col%dz(c,0) > 0._r8) then
           if (t_nvp_col(c) >= tfrz) then
             ! For unfrozen soil
              vol_ice = min(watsat_nvp, h2osoi_ice(c,0)/(col%dz(c,0)*denice))
              eff_porosity = watsat_nvp-vol_ice
              vol_liq = min(eff_porosity, h2osoi_liq(c,0)/(col%dz(c,0)*denh2o))
           else
              ! For frozen soil, assume NVP water content is at residual (unavailable for evaporation)
              vol_liq = watres_nvp            
           end if
        else
           vol_liq = 0._r8
        end if

        ! --- NVP van Genuchten matric potential and hydraulic conductivity ---
        call NVPWaterRetentionCurve(vol_liq, eff_porosity, n_van_nvp, alpha_van_nvp, &
             watsat_nvp, watres_nvp, psi_nvp)
        call NVPHydraulicConductivity(vol_liq, eff_porosity, n_van_nvp, watsat_nvp, watres_nvp, &
             ksat_nvp, khydr_nvp)
        K_nvp_mms = khydr_nvp * 1000._r8                                   ! m/s → mm/s

        ! --- Bidirectional Darcy flux with soil layer 1 ---
        ! Uses previous-timestep smp_l and hk_l (explicit time stepping)
        smp_soil1  = smp_l(c,1)                                            ! [mm]
        K_soil1    = hk_l(c,1)                                             ! [mm/s]

        ! Harmonic-mean interface conductivity
        if (K_nvp_mms + K_soil1 > 0._r8) then
           K_interface = 2._r8 * K_nvp_mms * K_soil1 / (K_nvp_mms + K_soil1)
        else
           K_interface = 0._r8
        end if

        ! Distance between layer centres [mm]
        dz_iface_mm = (0.5_r8 * col%dz(c,0) + 0.5_r8 * col%dz(c,1)) * 1000._r8

        ! Darcy flux: q = K * (grad_psi + gravity), positive = downward
        q01 = K_interface * ((psi_nvp - smp_soil1) / dz_iface_mm + 1.0_r8)

        ! --- Update h2osoi_liq(c,0): add infl, subtract evap; cannot go negative ---
        h2osoi_net = h2osoi_liq(c,0) &
             + (qflx_nvp_infl_col(c) - qflx_ev_nvp_col(c)) * dtime        ! [kg m-2]
        h2osoi_net = max(0._r8, h2osoi_net)

        if (q01 >= 0._r8) then
           ! Downward drainage: cap by available NVP liquid water
           qflx_nvp_drain_col(c) = min(q01, h2osoi_net / dtime)
        else
           ! Upward absorption from soil: cap by current soil layer-1 liquid water
           qflx_nvp_drain_col(c) = max(q01, -h2osoi_liq(c,1) / dtime)
        end if

        h2osoi_liq(c,0) = h2osoi_net - qflx_nvp_drain_col(c) * dtime      ! [kg m-2]

        ! --- Step 6: Saturation excess — flush water above watsat to drain ---
        if (col%dz(c,0) > 0._r8) then
           vol_liq = h2osoi_liq(c,0) / (denh2o * col%dz(c,0))
           if (vol_liq > watsat_nvp) then
              ! Excess water above saturation drains immediately to soil layer 1
              satfrac = (vol_liq - watsat_nvp) * denh2o * col%dz(c,0)     ! excess [kg m-2]
              qflx_nvp_drain_col(c) = qflx_nvp_drain_col(c) + satfrac / dtime
              h2osoi_liq(c,0)       = h2osoi_liq(c,0) - satfrac
           end if
        end if

        ! --- Update fwet_nvp (saturation fraction passed to FATES) ---
        if (col%dz(c,0) > 0._r8) then
           vol_liq = h2osoi_liq(c,0) / (denh2o * col%dz(c,0))
           satfrac = (vol_liq - watres_nvp) / max(watsat_nvp - watres_nvp, 1.e-10_r8)
           fwet_nvp_col(c) = max(0._r8, min(1._r8, satfrac))
        else
           fwet_nvp_col(c) = 0._r8
        end if

        ! [PORTED by Hui Tang: sync diagnostic copies for history output]
        ! h2onvp_col mirrors h2osoi_liq(c,0) [kg m-2 = mm H2O]; vwc_nvp_col is volumetric [m3 m-3]
        h2onvp_col(c) = h2osoi_liq(c,0)
        if (col%dz(c,0) > 0._r8) then
           vwc_nvp_col(c) = h2osoi_liq(c,0) / (denh2o * col%dz(c,0))
        else
           vwc_nvp_col(c) = 0._r8
        end if

      end do

    end associate

  end subroutine NVPWaterBalance_Column

  ! ===========================================================================

  subroutine NVPLayerRestart(bounds, ncid, flag)
    ! -------------------------------------------------------------------------
    ! [PORTED by Hui Tang: restart NVP column geometry and layer-state variables]
    !
    ! Saves/restores the four column-level variables that define NVP layer
    ! presence and geometry.  Without these, jbot_sno and col%dz(c,0) would
    ! revert to zero at restart, causing CLM physics to miss the NVP layer
    ! until the first FATES dynamics call.
    !
    ! Note: t_soisno(:,0), h2osoi_liq(:,0), h2osoi_ice(:,0) are already
    !       covered by the standard T_SOISNO / H2OSOI_LIQ / H2OSOI_ICE restart
    !       variables (dim2name='levtot', which spans -nlevsno+1:nlevgrnd).
    !       h2onvp_col and t_nvp_col are restarted in WaterStateType and
    !       TemperatureType respectively.
    ! -------------------------------------------------------------------------
    use ncdio_pio   , only : file_desc_t, ncd_double, ncd_int
    ! [PORTED by Hui Tang: use restUtilMod (not restFileMod) for restartvar — restFileMod
    !  uses clm_instMod, creating clm_instMod ↔ NVPLayerDynamicsMod cycle in build deps]
    use restUtilMod , only : restartvar
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds
    type(file_desc_t) , intent(inout) :: ncid
    character(len=*)  , intent(in)    :: flag   ! 'define', 'write', or 'read'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar
    !-----------------------------------------------------------------------

    ! Column NVP layer thickness [m]
    call restartvar(ncid=ncid, flag=flag, varname='DZ_NVP', xtype=ncd_double, &
         dim1name='column', &
         long_name='NVP (moss/lichen) layer thickness', units='m', &
         interpinic_flag='interp', readvar=readvar, data=col%dz_nvp)
    if (flag == 'read' .and. .not. readvar) then
       col%dz_nvp(bounds%begc:bounds%endc) = 0._r8
    end if

    ! Column NVP fractional coverage [-]
    call restartvar(ncid=ncid, flag=flag, varname='FRAC_NVP', xtype=ncd_double, &
         dim1name='column', &
         long_name='NVP (moss/lichen) fractional coverage', units='1', &
         interpinic_flag='interp', readvar=readvar, data=col%frac_nvp)
    if (flag == 'read' .and. .not. readvar) then
       col%frac_nvp(bounds%begc:bounds%endc) = 0._r8
    end if

    ! [PORTED by Hui Tang: nvp_layer_active is fully redundant with jbot_sno (active iff jbot_sno == -1).
    !  restartvar has no logical-array overload, so we restart only JBOT_SNO and derive the flag below.]

    ! Bottom index of active snow: 0 = no NVP, -1 = NVP present at layer 0
    call restartvar(ncid=ncid, flag=flag, varname='JBOT_SNO', xtype=ncd_int, &
         dim1name='column', &
         long_name='bottom index of active snow layers (0 or -1 with NVP)', units='', &
         interpinic_flag='interp', readvar=readvar, data=col%jbot_sno)
    if (flag == 'read' .and. .not. readvar) then
       col%jbot_sno(bounds%begc:bounds%endc) = 0
    end if

    ! Derive nvp_layer_active from jbot_sno on read (covers both restart and cold-start paths)
    if (flag == 'read') then
       col%nvp_layer_active(bounds%begc:bounds%endc) = &
            (col%jbot_sno(bounds%begc:bounds%endc) == -1)
    end if

  end subroutine NVPLayerRestart

end module NVPLayerDynamicsMod
