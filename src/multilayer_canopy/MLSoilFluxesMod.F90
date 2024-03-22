module MLSoilFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate soil surface temperature and energy balance
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilFluxes
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilFluxes (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Soil surface temperature and energy balance
    !
    ! !USES:
    use MLWaterVaporMod, only : SatVap, LatVap
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p                  ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(out) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lambda                        ! Latent heat of vaporization (J/mol)
    real(r8) :: gws                           ! Soil conductance for water vapor (mol H2O/m2/s)
    real(r8) :: gw                            ! Total conductance for water vapor (mol H2O/m2/s)
    real(r8) :: esat                          ! Saturation vapor pressure (Pa)
    real(r8) :: desat                         ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(r8) :: qsat                          ! Saturation vapor pressure of air (mol/mol)
    real(r8) :: dqsat                         ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: num1, num2, num3, num4, den   ! Intermediate terms
    real(r8) :: err                           ! Surface energy imbalance (W/m2)
    !---------------------------------------------------------------------

    associate ( &
                                                           ! *** Input ***
    tref        => mlcanopy_inst%tref_forcing         , &  ! Air temperature at reference height (K)
    pref        => mlcanopy_inst%pref_forcing         , &  ! Air pressure at reference height (Pa)
    rhomol      => mlcanopy_inst%rhomol_forcing       , &  ! Molar density at reference height (mol/m3)
    cpair       => mlcanopy_inst%cpair_forcing        , &  ! Specific heat of air (constant pressure) at reference height (J/mol/K)
    rnsoi       => mlcanopy_inst%rnsoi_soil           , &  ! Net radiation: ground (W/m2)
    rhg         => mlcanopy_inst%rhg_soil             , &  ! Relative humidity of airspace at soil surface (fraction)
    soilres     => mlcanopy_inst%soilres_soil         , &  ! Soil evaporative resistance (s/m)
    gac0        => mlcanopy_inst%gac0_soil            , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    soil_t      => mlcanopy_inst%soil_t_soil          , &  ! Temperature of first snow/soil layer (K)
    soil_dz     => mlcanopy_inst%soil_dz_soil         , &  ! Depth to temperature of first snow/soil layer (m)
    soil_tk     => mlcanopy_inst%soil_tk_soil         , &  ! Thermal conductivity of first snow/soil layer (W/m/K)
    tg_bef      => mlcanopy_inst%tg_bef_soil          , &  ! Soil surface temperature for previous timestep (K)
    tair        => mlcanopy_inst%tair_profile         , &  ! Canopy layer air temperature (K)
    eair        => mlcanopy_inst%eair_profile         , &  ! Canopy layer vapor pressure (Pa)
                                                           ! *** Output ***
    shsoi       => mlcanopy_inst%shsoi_soil           , &  ! Sensible heat flux: ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi_soil           , &  ! Latent heat flux: ground (W/m2)
    gsoi        => mlcanopy_inst%gsoi_soil            , &  ! Soil heat flux (W/m2)
    etsoi       => mlcanopy_inst%etsoi_soil           , &  ! Water vapor flux: ground (mol H2O/m2/s)
    tg          => mlcanopy_inst%tg_soil              , &  ! Soil surface temperature (K)
    eg          => mlcanopy_inst%eg_soil                &  ! Soil surface vapor pressure (Pa)
    )

    ! Latent heat of vaporization

    lambda = LatVap(tref(p))

    ! Soil conductance to water vapour diffusion

    gws = 1._r8 / soilres(p)                         ! s/m -> m/s
    gws = gws * rhomol(p)                            ! m/s -> mol H2O/m2/s
    gw = gac0(p) * gws / (gac0(p) + gws)             ! total conductance

    ! Saturation vapor pressure at ground temperature (Pa -> mol/mol)

    call SatVap (tg_bef(p), esat, desat)
    qsat = esat / pref(p) ; dqsat = desat / pref(p)

    ! Calculate soil surface temperature

    num1 = cpair(p) * gac0(p)
    num2 = lambda * gw
    num3 = soil_tk(p) / soil_dz(p)
    num4 = rnsoi(p) - num2 * rhg(p) * (qsat - dqsat * tg_bef(p)) + num3 * soil_t(p)
    den = num1 + num2 * dqsat * rhg(p) + num3
    tg(p) = (num1 * tair(p,1) + num2 * (eair(p,1)/pref(p)) + num4) / den

    ! Sensible heat flux

    shsoi(p) = cpair(p) * (tg(p) - tair(p,1)) * gac0(p)

    ! Latent heat flux

    eg(p) = rhg(p) * (esat + desat * (tg(p) - tg_bef(p)))
    lhsoi(p) = lambda * (eg(p) - eair(p,1)) / pref(p) * gw

    ! Soil heat flux

    gsoi(p) = soil_tk(p) * (tg(p) - soil_t(p)) / soil_dz(p)

    ! Error check

    err = rnsoi(p) - shsoi(p) - lhsoi(p) - gsoi(p)
    if (abs(err) > 0.001_r8) then
       call endrun (msg=' ERROR: SoilFluxes: energy balance error')
    end if

    ! Water vapor flux: W/m2 -> mol H2O/m2/s

    etsoi(p) = lhsoi(p) / lambda

    end associate
  end subroutine SoilFluxes

end module MLSoilFluxesMod
