module MLLeafBoundaryLayerMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf boundary layer conductance
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
  public :: LeafBoundaryLayer
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LeafBoundaryLayer (num_filter, filter, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf boundary layer conductance
    !
    ! !USES:
    use clm_varcon, only : tfrz, grav
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : visc0, dh0, dv0, dc0, gb_factor
    use MLclm_varctl, only : gb_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter        ! Number of patches in filter
    integer, intent(in) :: filter(:)         ! Patch filter
    integer, intent(in) :: il                ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                           ! Filter index
    integer  :: p                            ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                           ! Aboveground layer index
    real(r8) :: visc                         ! Kinematic viscosity (m2/s)
    real(r8) :: dh                           ! Molecular diffusivity, heat (m2/s)
    real(r8) :: dv                           ! Molecular diffusivity, H2O (m2/s)
    real(r8) :: dc                           ! Molecular diffusivity, CO2 (m2/s)
    real(r8) :: fac                          ! Correction factor for temperature and pressure
    real(r8) :: nu_lam                       ! Forced convection - laminar: Nusselt number (dimensionless)
    real(r8) :: shv_lam                      ! Forced convection - laminar: Sherwood number, H2O (dimensionless)
    real(r8) :: shc_lam                      ! Forced convection - laminar: Sherwood number, CO2 (dimensionless)
    real(r8) :: nu_turb                      ! Forced convection - turbulent: Nusselt number (dimensionless)
    real(r8) :: shv_turb                     ! Forced convection - turbulent: Sherwood number, H2O (dimensionless)
    real(r8) :: shc_turb                     ! Forced convection - turbulent: Sherwood number, CO2 (dimensionless)
    real(r8) :: nu_free                      ! Free convection: Nusselt number (dimensionless)
    real(r8) :: shv_free                     ! Free convection: Sherwood number, H2O (dimensionless)
    real(r8) :: shc_free                     ! Free convection: Sherwood number, CO2 (dimensionless)
    real(r8) :: nu                           ! Nusselt number (dimensionless)
    real(r8) :: shv                          ! Sherwood number, H2O (dimensionless)
    real(r8) :: shc                          ! Sherwood number, CO2 (dimensionless)
    real(r8) :: pr                           ! Prandtl number (dimensionless)
    real(r8) :: scv                          ! Schmidt number, H2O (dimensionless)
    real(r8) :: scc                          ! Schmidt number, CO2 (dimensionless)
    real(r8) :: re                           ! Reynolds number (dimensionless)
    real(r8) :: gr                           ! Grashof number (dimensionless)
    !---------------------------------------------------------------------

    associate ( &
                                                   ! *** Input ***
    dleaf     => pftcon%dleaf                 , &  ! CLM: Leaf dimension (m)
    tref      => mlcanopy_inst%tref_forcing   , &  ! Air temperature at reference height (K)
    pref      => mlcanopy_inst%pref_forcing   , &  ! Air pressure at reference height (Pa)
    rhomol    => mlcanopy_inst%rhomol_forcing , &  ! Molar density at reference height (mol/m3)
    ncan      => mlcanopy_inst%ncan_canopy    , &  ! Number of aboveground layers
    dpai      => mlcanopy_inst%dpai_profile   , &  ! Canopy layer plant area index (m2/m2)
    wind      => mlcanopy_inst%wind_profile   , &  ! Canopy layer wind speed (m/s)
    tair      => mlcanopy_inst%tair_profile   , &  ! Canopy layer air temperature (K)
    tleaf     => mlcanopy_inst%tleaf_leaf     , &  ! Leaf temperature (K)
                                                   ! *** Output ***
    gbh       => mlcanopy_inst%gbh_leaf       , &  ! Leaf boundary layer conductance: heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv_leaf       , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gbc       => mlcanopy_inst%gbc_leaf         &  ! Leaf boundary layer conductance: CO2 (mol CO2/m2 leaf/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Adjust diffusivity for temperature and pressure

       fac = 101325._r8 / pref(p) * (tref(p) / tfrz)**1.81_r8
       visc = visc0 * fac
       dh = dh0 * fac
       dv = dv0 * fac
       dc = dc0 * fac

       do ic = 1, ncan(p)

          if (dpai(p,ic) > 0._r8) then

             select case (gb_type)
             case (0)

                ! Use CLM5 simplification: units are m/s

                gbh(p,ic,il) = 0.005 * sqrt(wind(p,ic) / dleaf(patch%itype(p)))
                gbv(p,ic,il) = gbh(p,ic,il)
                gbc(p,ic,il) = gbv(p,ic,il) / 1.4_r8

             case (1:3)

                ! Use full equation set to calculate boundary layer conductances

                ! a. Reynolds number, Prandtl number, Schmidt numbers, and Grashof number

                re = wind(p,ic) * dleaf(patch%itype(p)) / visc
                pr  = visc / dh
                scv = visc / dv
                scc = visc / dc
                gr = grav * dleaf(patch%itype(p))**3 * max(tleaf(p,ic,il)-tair(p,ic), 0._r8) / (tair(p,ic) * visc * visc)

                ! b. Nusselt and Sherwood numbers depend on convection regime

                ! Forced convection

                ! (i) Laminar flow

                nu_lam  = gb_factor * 0.66_r8 *  pr**0.33_r8 * re**0.5_r8
                shv_lam = gb_factor * 0.66_r8 * scv**0.33_r8 * re**0.5_r8
                shc_lam = gb_factor * 0.66_r8 * scc**0.33_r8 * re**0.5_r8

                ! (ii) Turbulent flow

                nu_turb  = gb_factor * 0.036_r8 *  pr**0.33_r8 * re**0.8_r8
                shv_turb = gb_factor * 0.036_r8 * scv**0.33_r8 * re**0.8_r8
                shc_turb = gb_factor * 0.036_r8 * scc**0.33_r8 * re**0.8_r8

                ! Free convection

                nu_free  = 0.54_r8 *  pr**0.25_r8 * gr**0.25_r8
                shv_free = 0.54_r8 * scv**0.25_r8 * gr**0.25_r8
                shc_free = 0.54_r8 * scc**0.25_r8 * gr**0.25_r8

                ! Choose flow regimes to use

                select case (gb_type)
                case (1)

                   ! Use only laminar flow

                   nu = nu_lam
                   shv = shv_lam
                   shc = shc_lam

                case (2)

                   ! Use laminar and turbulent flow

                   nu = max(nu_lam, nu_turb)
                   shv = max(shv_lam, shv_turb)
                   shc = max(shc_lam, shc_turb)

                case (3)

                   ! Both forced and free convection occur together

                   nu = max(nu_lam, nu_turb) + nu_free
                   shv = max(shv_lam, shv_turb) + shv_free
                   shc = max(shc_lam, shc_turb) + shc_free

                end select

                ! Boundary layer conductances (m/s)

                gbh(p,ic,il) = dh *  nu / dleaf(patch%itype(p))
                gbv(p,ic,il) = dv * shv / dleaf(patch%itype(p))
                gbc(p,ic,il) = dc * shc / dleaf(patch%itype(p))

             case default

                call endrun (msg=' ERROR: LeafBoundaryLayer: gb_type not valid')

             end select

             ! Convert conductances (m/s) to (mol/m2/s)

             gbh(p,ic,il) = gbh(p,ic,il) * rhomol(p)
             gbv(p,ic,il) = gbv(p,ic,il) * rhomol(p)
             gbc(p,ic,il) = gbc(p,ic,il) * rhomol(p)

          else

             gbh(p,ic,il) = 0._r8
             gbv(p,ic,il) = 0._r8
             gbc(p,ic,il) = 0._r8

          end if

       end do
    end do

    end associate
  end subroutine LeafBoundaryLayer

end module MLLeafBoundaryLayerMod
