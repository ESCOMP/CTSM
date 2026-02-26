module MLclm_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing multilayer canopy model constants
  !
  ! !USES:
  use clm_varcon, only : spval
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  save

  ! Physical constants for multilayer canopy

  real(r8) :: rgas = 8.31446_r8                        ! Universal gas constant (J/K/mol)
  real(r8) :: mmdry = 28.97_r8 * 1.e-03_r8             ! Molecular mass of dry air (kg/mol)
  real(r8) :: mmh2o = 18.02_r8 * 1.e-03_r8             ! Molecular mass of water vapor (kg/mol)
  real(r8) :: cpd = 1005._r8                           ! Specific heat of dry air at constant pressure (J/kg/K)
  real(r8) :: cpw = 1846._r8                           ! Specific heat of water vapor at constant pressure (J/kg/K)
  real(r8) :: visc0 = 13.3e-06_r8                      ! Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
  real(r8) :: dh0 = 18.9e-06_r8                        ! Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
  real(r8) :: dv0 = 21.8e-06_r8                        ! Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
  real(r8) :: dc0 = 13.8e-06_r8                        ! Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
  real(r8) :: lapse_rate = 0.0098_r8                   ! Temperature lapse rate (K/m)

  ! Adjustable parameters for multilayer canopy

  !----------------------------------------------------!
  ! Leaf photosynthesis
  !----------------------------------------------------!
  real(r8) :: kc25 = 404.9_r8                          ! Michaelis-Menten constant for CO2 at 25C (umol/mol)
  real(r8) :: kcha = 79430._r8                         ! Activation energy for kc (J/mol)
  real(r8) :: ko25 = 278.4_r8                          ! Michaelis-Menten constant for O2 at 25C (mmol/mol)
  real(r8) :: koha = 36380._r8                         ! Activation energy for ko (J/mol)
  real(r8) :: cp25 = 42.75_r8                          ! CO2 compensation point at 25C (umol/mol)
  real(r8) :: cpha = 37830._r8                         ! Activation energy for cp (J/mol)

  real(r8) :: vcmaxha_noacclim = 65330._r8             ! Activation energy for vcmax without acclimation (J/mol)
  real(r8) :: vcmaxha_acclim   = 72000._r8             ! Activation energy for vcmax with acclimation (J/mol)
  real(r8) :: vcmaxhd_noacclim = 150000._r8            ! Deactivation energy for vcmax without acclimation (J/mol)
  real(r8) :: vcmaxhd_acclim   = 200000._r8            ! Deactivation energy for vcmax with acclimation (J/mol)
  real(r8) :: vcmaxse_noacclim = 490._r8               ! Entropy term for vcmax without acclimation (J/mol/K)
  real(r8) :: vcmaxse_acclim   = spval                 ! Entropy term for vcmax with acclimation (J/mol/K)

  real(r8) :: jmaxha_noacclim = 43540._r8              ! Activation energy for jmax without acclimation (J/mol)
  real(r8) :: jmaxha_acclim   = 50000._r8              ! Activation energy for jmax with acclimation (J/mol)
  real(r8) :: jmaxhd_noacclim = 150000._r8             ! Deactivation energy for jmax without acclimation (J/mol)
  real(r8) :: jmaxhd_acclim   = 200000._r8             ! Deactivation energy for jmax with acclimation (J/mol)
  real(r8) :: jmaxse_noacclim = 490._r8                ! Entropy term for jmax without acclimation (J/mol/K)
  real(r8) :: jmaxse_acclim   = spval                  ! Entropy term for jmax with acclimation (J/mol/K)

  real(r8) :: rdha = 46390._r8                         ! Activation energy for rd (J/mol)
  real(r8) :: rdhd = 150000._r8                        ! Deactivation energy for rd (J/mol)
  real(r8) :: rdse = 490._r8                           ! Entropy term for rd (J/mol/K)

  real(r8) :: jmax25_to_vcmax25_noacclim = 1.67_r8     ! Ratio of jmax to vcmax at 25C without acclimation (umol/umol)
  real(r8) :: jmax25_to_vcmax25_acclim   = spval       ! Ratio of jmax to vcmax at 25C with acclimation (umol/umol)
  real(r8) :: rd25_to_vcmax25_c3 = 0.015_r8            ! Ratio of rd to vcmax at 25C (C3) (umol/umol)
  real(r8) :: rd25_to_vcmax25_c4 = 0.025_r8            ! Ratio of rd to vcmax at 25C (C4) (umol/umol)
  real(r8) :: kp25_to_vcmax25_c4 = 0.02_r8             ! Ratio of kp to vcmax at 25C (C4) (mol/umol)

  real(r8) :: phi_psII = 0.70_r8                       ! C3: quantum yield of PS II
! real(r8) :: phi_psII = 0.85_r8                       ! C3: quantum yield of PS II
  real(r8) :: theta_j = 0.90_r8                        ! C3: empirical curvature parameter for electron transport rate
  real(r8) :: qe_c4 = 0.05_r8                          ! C4: quantum yield (mol CO2 / mol photons)

  real(r8) :: colim_c3a = 0.98_r8                      ! Empirical curvature parameter for C3 co-limitation (Ac, Aj)
  real(r8) :: colim_c3b = spval                        ! Empirical curvature parameter for C3 co-limitation (Ap)
  real(r8) :: colim_c4a = 0.80_r8                      ! Empirical curvature parameter for C4 co-limitation (Ac, Aj)
  real(r8) :: colim_c4b = 0.95_r8                      ! Empirical curvature parameter for C4 co-limitation (Ap)

  !----------------------------------------------------!
  ! Stomatal conductance
  !----------------------------------------------------!
  real(r8) :: dh2o_to_dco2 = 1.6_r8                    ! Diffusivity H2O / Diffusivity CO2
  real(r8) :: rh_min_BB = 0.2_r8                       ! Minimum relative humidity of air for Ball-Berry stomatal conductance (fraction)
  real(r8) :: vpd_min_MED = 100._r8                    ! Minimum vapor pressure deficit for Medlyn stomatal conductance (Pa)

  !----------------------------------------------------!
  ! Leaf heat capacity
  !----------------------------------------------------!
  real(r8) :: cpbio = 4188._r8 / 3._r8                 ! Specific heat of dry biomass (J/kg/K)
  real(r8) :: fcarbon = 0.5_r8                         ! Fraction of dry biomass that is carbon (kg C / kg DM)
  real(r8) :: fwater = 0.7_r8                          ! Fraction of fresh biomass that is water (kg H2O / kg FM)

  !----------------------------------------------------!
  ! Leaf boundary layer conductance
  !----------------------------------------------------!
  real(r8) :: gb_factor = 1.5_r8                       ! Empirical correction factor for Nu

  !----------------------------------------------------!
  ! Canopy interception
  !----------------------------------------------------!
  real(r8) :: dewmx = 0.1_r8                           ! Maximum allowed interception (kg H2O/m2 leaf)
  real(r8) :: maximum_leaf_wetted_fraction = 0.05_r8   ! Maximum fraction of leaf that can be wet
  real(r8) :: interception_fraction = 1.0_r8           ! Fraction of intercepted precipitation
  real(r8) :: fwet_exponent = 0.67_r8                  ! Exponent for wetted canopy fraction
  real(r8) :: clm45_interception_p1 = 0.25_r8          ! CLM4.5: Interception parameter
  real(r8) :: clm45_interception_p2 = -0.50_r8         ! CLM4.5: Interception parameter

  !----------------------------------------------------!
  ! Solar radiation
  !----------------------------------------------------!
  real(r8) :: chil_min = -0.4_r8                       ! Minimum value for xl leaf angle orientation parameter
  real(r8) :: chil_max = 0.6_r8                        ! Maximum value for xl leaf angle orientation parameter
  real(r8) :: kb_max = 40._r8                          ! Maximum value for direct beam extinction coefficient
  real(r8) :: J_to_umol = 4.6_r8                       ! PAR conversion from W/m2 to umol/m2/s (umol/J)

  !----------------------------------------------------!
  ! Longwave radiation
  !----------------------------------------------------!
  real(r8) :: emg = 0.96_r8                            ! Ground (soil) emissivity

  !----------------------------------------------------!
  ! Roughness sublayer parameterization
  !----------------------------------------------------!
  real(r8) :: cd = 0.25_r8                             ! RSL - drag coefficient for canopy elements (dimensionless)
  real(r8) :: beta_neutral_max = 0.35_r8               ! RSL - maximum value for beta in neutral conditions
  real(r8) :: cr = 0.3_r8                              ! RSL - parameter to calculate beta_neutral
  real(r8) :: c2 = 0.5_r8                              ! RSL - depth scale multiplier
  real(r8) :: Pr0 = 0.5_r8                             ! RSL - neutral value for Pr (Sc)
  real(r8) :: Pr1 = 0.3_r8                             ! RSL - magnitude of variation of Pr (Sc) with stability
  real(r8) :: Pr2 = 2.0_r8                             ! RSL - scale of variation of Pr (Sc) with stability
  real(r8) :: z0mg = 0.01_r8                           ! RSL - roughness length of ground (m)

  ! Limits placed on various variables

  real(r8) :: wind_forc_min = 1.0_r8                   ! Minimum wind speed at forcing height (m/s)
  real(r8) :: eta_max = 20._r8                         ! Maximum value for "eta" parameter (used to constrain lm/beta)
  real(r8) :: zeta_min = -2.0_r8                       ! Minimum value for Monin-Obukhov zeta parameter
  real(r8) :: zeta_max = 1.0_r8                        ! Maximum value for Monin-Obukhov zeta parameter
  real(r8) :: beta_min = 0.2_r8                        ! Minimum value for H&F beta parameter
  real(r8) :: beta_max = 0.5_r8                        ! Maximum value for H&F beta parameter
  real(r8) :: wind_min = 0.1_r8                        ! Minimum wind speed in canopy (m/s)
  real(r8) :: ra_max = 500._r8                         ! Maximum aerodynamic resistance (s/m)

  !----------------------------------------------------!
  ! Constants used in the RSL psihat look-up tables
  ! are initialized in subroutine LookupPsihatINI
  !----------------------------------------------------!

  integer, parameter :: nZ = 276, nL = 41              ! Dimensions of RSL psihat look-up tables
  real(r8) :: zdtgridM(nZ,1)                           ! Grid of zdt on which psihat is given for momentum
  real(r8) :: dtLgridM(1,nL)                           ! Grid of dtL on which psihat is given for momentum
  real(r8) :: psigridM(nZ,nL)                          ! Grid of psihat values for momentum
  real(r8) :: zdtgridH(nZ,1)                           ! Grid of zdt on which psihat is given for heat
  real(r8) :: dtLgridH(1,nL)                           ! Grid of dtL on which psihat is given for heat
  real(r8) :: psigridH(nZ,nL)                          ! Grid of psihat values for heat

end module MLclm_varcon
