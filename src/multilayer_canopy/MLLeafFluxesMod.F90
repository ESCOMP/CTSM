module MLLeafFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf temperature and energy fluxes
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
  public :: LeafFluxes
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LeafFluxes (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf temperature and energy fluxes
    !
    ! !USES:
    use MLclm_varctl, only : dtime_substep
    use MLWaterVaporMod, only : SatVap, LatVap
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p            ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic           ! Aboveground layer index
    integer, intent(in) :: il           ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dtime                   ! Model time step (s)
    real(r8) :: esat                    ! Saturation vapor pressure (Pa)
    real(r8) :: desat                   ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(r8) :: qsat                    ! Saturation vapor pressure of air (mol/mol)
    real(r8) :: dqsat                   ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: lambda                  ! Latent heat of vaporization (J/mol)
    real(r8) :: gleaf                   ! Leaf conductance for transpiration (mol H2O/m2 leaf/s)
    real(r8) :: gw                      ! Total conductance including evaporation (mol H2O/m2 leaf/s)
    real(r8) :: num1, num2, num3, den   ! Intermediate calculation
    real(r8) :: err                     ! Energy balance error (W/m2)
    !---------------------------------------------------------------------

    associate ( &
                                                   ! *** Input ***
    tref      => mlcanopy_inst%tref_forcing   , &  ! Air temperature at reference height (K)
    pref      => mlcanopy_inst%pref_forcing   , &  ! Air pressure at reference height (Pa)
    cpair     => mlcanopy_inst%cpair_forcing  , &  ! Specific heat of air (constant pressure) at reference height (J/mol/K)
    dpai      => mlcanopy_inst%dpai_profile   , &  ! Canopy layer plant area index (m2/m2)
    tair      => mlcanopy_inst%tair_profile   , &  ! Canopy layer air temperature (K)
    eair      => mlcanopy_inst%eair_profile   , &  ! Canopy layer vapor pressure (Pa)
    cpleaf    => mlcanopy_inst%cpleaf_profile , &  ! Canopy layer leaf heat capacity (J/m2 leaf/K)
    fwet      => mlcanopy_inst%fwet_profile   , &  ! Canopy layer fraction of plant area index that is wet
    fdry      => mlcanopy_inst%fdry_profile   , &  ! Canopy layer fraction of plant area index that is green and dry
    gbh       => mlcanopy_inst%gbh_leaf       , &  ! Leaf boundary layer conductance: heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv_leaf       , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gs        => mlcanopy_inst%gs_leaf        , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    rnleaf    => mlcanopy_inst%rnleaf_leaf    , &  ! Leaf net radiation (W/m2 leaf)
    tleaf_bef => mlcanopy_inst%tleaf_bef_leaf , &  ! Leaf temperature for previous timestep (K)
                                                   ! *** Output ***
    tleaf     => mlcanopy_inst%tleaf_leaf     , &  ! Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf_leaf    , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf_leaf    , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf    => mlcanopy_inst%lhleaf_leaf    , &  ! Leaf latent heat flux (W/m2 leaf)
    evleaf    => mlcanopy_inst%evleaf_leaf    , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    trleaf    => mlcanopy_inst%trleaf_leaf      &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    )

    ! Time step (s)

    dtime = dtime_substep

    ! Latent heat of vaporization

    lambda = LatVap(tref(p))

    if (dpai(p,ic) > 0._r8) then

       ! Saturation vapor pressure (Pa -> mol/mol)

       call SatVap (tleaf_bef(p,ic,il), esat, desat)
       qsat = esat / pref(p) ; dqsat = desat / pref(p)

       ! Leaf conductance for transpiration

       gleaf = gs(p,ic,il) * gbv(p,ic,il) / (gs(p,ic,il) + gbv(p,ic,il))

       ! Total conductance (including evaporation)

       gw = gleaf * fdry(p,ic) + gbv(p,ic,il) * fwet(p,ic)

       ! Linearized leaf temperature calculation that balances the energy budget

       num1 = 2._r8 * cpair(p) * gbh(p,ic,il)
       num2 = lambda * gw 
       num3 = rnleaf(p,ic,il) - lambda * gw * (qsat - dqsat * tleaf_bef(p,ic,il)) &
            + cpleaf(p,ic) / dtime * tleaf_bef(p,ic,il)
       den = cpleaf(p,ic) / dtime + num1 + num2 * dqsat
       tleaf(p,ic,il) = (num1 * tair(p,ic) + num2 * eair(p,ic) / pref(p) + num3) / den

       ! Storage flux

       stleaf(p,ic,il) = (tleaf(p,ic,il) - tleaf_bef(p,ic,il)) * cpleaf(p,ic) / dtime

       ! Sensible heat flux

       shleaf(p,ic,il) = 2._r8 * cpair(p) * (tleaf(p,ic,il) - tair(p,ic)) * gbh(p,ic,il)

       ! Transpiration and evaporation water fluxes: mol H2O/m2/s

       num1 = qsat + dqsat * (tleaf(p,ic,il) - tleaf_bef(p,ic,il)) - eair(p,ic) / pref(p)
       trleaf(p,ic,il) = gleaf * fdry(p,ic) * num1
       evleaf(p,ic,il) = gbv(p,ic,il) * fwet(p,ic) * num1

       ! Latent heat flux

       lhleaf(p,ic,il) = (trleaf(p,ic,il) + evleaf(p,ic,il)) * lambda

       ! Error check

       err = rnleaf(p,ic,il) - shleaf(p,ic,il) - lhleaf(p,ic,il) - stleaf(p,ic,il)
       if (abs(err) > 1.e-03_r8) then
          call endrun (msg=' ERROR: LeafFluxes: energy balance error')
       end if

    else

       tleaf(p,ic,il) = tair(p,ic)
       stleaf(p,ic,il) = 0._r8
       shleaf(p,ic,il) = 0._r8
       lhleaf(p,ic,il) = 0._r8
       evleaf(p,ic,il) = 0._r8
       trleaf(p,ic,il) = 0._r8

    end if

    end associate
  end subroutine LeafFluxes

end module MLLeafFluxesMod
