module MLclm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing multilayer canopy model run control variables
  !
  ! !USES:
  use clm_varcon, only : ispval
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none

  ! Model run options

  integer  :: light_type = 2            ! Radiative transfer: Norman (1) or two-stream approximation (2)
  integer  :: gs_type = 2               ! Stomatal conductance: Medlyn (0), Ball-Berry (1), or WUE optimization (2)
  integer  :: colim_type = 1            ! Photosynthesis co-limitation: minimum rate (0) or co-limited rate (1)
  integer  :: acclim_type = 1           ! Photosynthesis acclimation: off (0) or on (1)
  real(r8) :: kn_val = -999._r8         ! Canopy nitrogen profile: for a user-specified Kn set kn_val to desired value
  integer  :: turb_type = 1             ! Turbulence parameterization: well-mixed assumption (0) or H&F roughness sublayer (1)
  integer  :: pad_type = 1              ! Plant area density: beta distribution (1) or uniform distribution (2)
  integer  :: mlcan_to_clm = 0          ! Pass multilayer canopy fluxes to CLM for use by CAM: no (0), yes (1)
  integer  :: ml_vert_init = ispval     ! Initialize multilayer canopy vertical structure and profiles: no (0), yes (1)

  ! Model sub-timestep (s) is 5 minutes

  real(r8), save :: dtime_substep = 5._r8 * 60._r8

  ! Archival or permanent disk full pathname for RSL psihat look-up tables

  character(len=256), save :: rslfile = '../rsl_lookup_tables/psihat.nc'

end module MLclm_varctl
