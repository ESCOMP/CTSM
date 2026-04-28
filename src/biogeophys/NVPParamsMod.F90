module NVPParamsMod

  ! [PORTED by Hui Tang: centralized CLM-side NVP (moss/lichen) physics parameters]
  !
  ! All tunable parameters for the NVP layer in CLM are declared here with their
  ! default values.  They are readable at runtime via the nvp_inparm namelist group
  ! in user_nl_clm (or the model's standard namelist input file).  Default values
  ! reproduce the original hardcoded constants and leave model behaviour unchanged
  ! when nvp_inparm is absent from the namelist.
  !
  ! Usage:
  !   use NVPParamsMod, only : nvp_frac_min, watsat_nvp, ...
  !
  ! Parameters:
  !   nvp_frac_min    — min fractional coverage to activate the NVP layer
  !   rnvp_min        — minimum surface evaporation resistance (fully saturated) [s m-1]
  !   rnvp_amp        — resistance amplitude as NVP dries [s m-1]
  !   rnvp_exp        — exponent of resistance–dryness curve [-]
  !   ksat_nvp        — saturated hydraulic conductivity [m s-1]
  !   n_van_nvp       — van Genuchten shape parameter n [-]
  !   alpha_van_nvp   — van Genuchten inverse air-entry pressure alpha [cm-1]
  !   watsat_nvp      — saturated volumetric water content (porosity) [m3 m-3]
  !   watres_nvp      — residual volumetric water content [m3 m-3]

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  public

  ! Activation threshold
  real(r8) :: nvp_frac_min    = 1.0e-6_r8  ! min NVP coverage fraction to activate layer [-]

  ! Evaporation resistance (van de Griend & Owe, 1994 / Daamen & Simmonds)
  ! rnvp = rnvp_min + rnvp_amp * (1 - satfrac)^rnvp_exp  [s m-1]
  real(r8) :: rnvp_min        = 10.0_r8    ! minimum surface resistance when saturated [s m-1]
  real(r8) :: rnvp_amp        = 500.0_r8   ! amplitude of resistance increase when dry  [s m-1]
  real(r8) :: rnvp_exp        = 3.0_r8     ! exponent of dryness function               [-]
  ! [PORTED by Hui Tang: surface resistance applied when NVP is frozen — typical literature
  !  range for ice/snow surfaces is ~1000-3000 s/m; default 1500 s/m suppresses but does
  !  not zero sublimation. Set to a very large value (e.g. 1e6) to disable evap when frozen.]
  real(r8) :: rnvp_ice        = 1500.0_r8  ! NVP resistance when frozen [s m-1]

  ! Hydraulic properties (Mualem-van Genuchten)
  real(r8) :: ksat_nvp        = 1.0e-4_r8  ! saturated hydraulic conductivity [m s-1]
  real(r8) :: n_van_nvp       = 1.5_r8     ! van Genuchten shape parameter n  [-]
  real(r8) :: alpha_van_nvp   = 0.01_r8    ! van Genuchten alpha              [cm-1]
  real(r8) :: watsat_nvp      = 0.85_r8    ! saturated volumetric water content [m3 m-3]
  real(r8) :: watres_nvp      = 0.05_r8    ! residual volumetric water content  [m3 m-3]

  ! Thermal properties of the dry NVP matrix (Farouki-style mixing with water/ice)
  ! [PORTED by Hui Tang: NVP dry-matrix thermal parameters for SoilTemperatureMod]
  real(r8) :: thk_dry_nvp     = 0.05_r8    ! dry NVP thermal conductivity       [W m-1 K-1]
  real(r8) :: csol_nvp        = 0.58e6_r8  ! dry NVP volumetric heat capacity   [J m-3 K-1]

end module NVPParamsMod
