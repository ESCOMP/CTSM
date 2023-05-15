module MLPlantHydraulicsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate plant hydraulics
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
  public :: PlantResistance        ! Calculate whole-plant resistance
  public :: SoilResistance         ! Calculate soil resistance and water uptake
  !-----------------------------------------------------------------------

contains

  subroutine PlantResistance (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Calculate whole-plant leaf-specific conductance (soil-to-leaf)
    ! 
    ! !USES:
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter       ! Number of patches in filter
    integer, intent(in) :: filter(:)        ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                          ! Filter index
    integer  :: p                           ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                          ! Aboveground layer index
    real(r8) :: rplant                      ! Aboveground plant hydraulic resistance (MPa.s.m2/mmol H2O)
    !---------------------------------------------------------------------

    associate ( &
                                                  ! *** Input ***
    gplant_SPA => pftcon%gplant_SPA          , &  ! CLM (new): Stem (xylem-to-leaf) hydraulic conductance (mmol H2O/m2 leaf area/s/MPa)
    ncan       => mlcanopy_inst%ncan_canopy  , &  ! Number of layers
    rsoil      => mlcanopy_inst%rsoil_soil   , &  ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
    dpai       => mlcanopy_inst%dpai_profile , &  ! Canopy layer plant area index (m2/m2)
    zs         => mlcanopy_inst%zs_profile   , &  ! Canopy layer height for scalar concentration and source (m)
                                                  ! *** Output ***
    lsc        => mlcanopy_inst%lsc_profile    &  ! Canopy layer leaf-specific conductance (mmol H2O/m2 leaf/s/MPa)
    )

    do fp = 1, num_filter
       p = filter(fp)
       do ic = 1, ncan(p)

          if (dpai(p,ic) > 0._r8) then

             ! Aboveground plant stem resistance, xylem-to-leaf (MPa.s.m2/mmol H2O)

!            rplant = zs(p,ic) / gplant_SPA(patch%itype(p))       ! gplant_SPA is conductivity (mmol H2O/m/s/MPa)
             rplant = 1._r8 / gplant_SPA(patch%itype(p))          ! gplant_SPA is conductance (mmol H2O/m2/s/MPa)

             ! Leaf specific conductance, soil-to-leaf (mmol H2O/m2/s/MPa)

             lsc(p,ic) = 1._r8 / (rsoil(p) + rplant)

          else

             lsc(p,ic) = 0._r8

          end if

       end do
    end do

    end associate
  end subroutine PlantResistance

  !-----------------------------------------------------------------------
  subroutine SoilResistance (num_filter, filter, &
  soilstate_inst, waterstatebulk_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Calculate soil hydraulic resistance and water uptake from each soil layer
    ! 
    ! !USES:
    use clm_varcon, only : pi => rpi, denh2o, grav
    use clm_varpar, only : nlevsoi
    use ColumnType, only : col
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : mmh2o
    use SoilStateType, only : soilstate_type
    use WaterStateBulkType, only : waterstatebulk_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter            ! Number of patches in filter
    integer, intent(in) :: filter(:)             ! Patch filter
    type(soilstate_type), intent(in) :: soilstate_inst
    type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                               ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                ! Column index for CLM g/l/c/p hierarchy
    integer  :: j                                ! Soil layer index
    real(r8) :: head                             ! Head of pressure  (MPa/m)
    real(r8) :: root_cross_sec_area              ! Root cross-sectional area (m2 root)
    real(r8) :: root_biomass_density             ! Root biomass density (g biomass / m3 soil) 
    real(r8) :: root_length_density              ! Root length density (m root / m3 soil) 
    real(r8) :: root_dist                        ! Mean distance between roots (m)
    real(r8) :: hk                               ! Hydraulic conductivity (mm H2O/s -> mmol H2O/m/s/MPa)
    real(r8) :: soilr1                           ! Soil-to-root resistance (MPa.s.m2/mmol H2O)
    real(r8) :: soilr2                           ! Root-to-stem resistance (MPa.s.m2/mmol H2O)
    real(r8) :: soilr                            ! Belowground resistance (MPa.s.m2/mmol H2O) 
    real(r8) :: smp_mpa(nlevsoi)                 ! Soil matric potential (MPa)
    real(r8) :: evap(nlevsoi)                    ! Maximum transpiration (mmol H2O/m2/s)
    real(r8) :: totevap                          ! Total maximum transpiration (mmol H2O/m2/s)
    !---------------------------------------------------------------------

    associate ( &
                                                               ! *** Input ***
    minlwp_SPA       => pftcon%minlwp_SPA                 , &  ! CLM (new): Minimum leaf water potential (MPa)
    root_radius_SPA  => pftcon%root_radius_SPA            , &  ! CLM (new): Fine root radius (m)
    root_density_SPA => pftcon%root_density_SPA           , &  ! CLM (new): Fine root density (g biomass / m3 root)
    root_resist_SPA  => pftcon%root_resist_SPA            , &  ! CLM (new): Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O)
    dz               => col%dz                            , &  ! CLM: Soil layer thickness (m)
    smp_l            => soilstate_inst%smp_l_col          , &  ! CLM: Soil layer matric potential (mm)
    hk_l             => soilstate_inst%hk_l_col           , &  ! CLM: Soil layer hydraulic conductivity (mm H2O/s)
    rootfr           => soilstate_inst%rootfr_patch       , &  ! CLM: Fraction of roots in each soil layer
    h2osoi_ice       => waterstatebulk_inst%h2osoi_ice_col    , &  ! CLM: Soil layer ice lens (kg H2O/m2)
    lai              => mlcanopy_inst%lai_canopy          , &  ! Leaf area index of canopy (m2/m2)
    root_biomass     => mlcanopy_inst%root_biomass_canopy , &  ! Fine root biomass (g biomass / m2)
                                                               ! *** Output ***
    psis             => mlcanopy_inst%psis_soil           , &  ! Weighted soil water potential (MPa)
    rsoil            => mlcanopy_inst%rsoil_soil          , &  ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
    soil_et_loss     => mlcanopy_inst%soil_et_loss_soil     &  ! Fraction of total transpiration from each soil layer (-)
    )

    head = denh2o * grav * 1.e-06_r8

    ! Soil and root resistances for each layer

    do fp = 1, num_filter
       p = filter(fp)
       c = patch%column(p)
       root_cross_sec_area = pi * root_radius_SPA(patch%itype(p))**2

       rsoil(p) = 0._r8
       do j = 1, nlevsoi

          ! Hydraulic conductivity and matric potential for each layer

          hk = hk_l(c,j) * (1.e-03_r8 / head)                       ! mm/s -> m/s -> m2/s/MPa
          hk = hk * denh2o / mmh2o * 1000._r8                       ! m2/s/MPa -> mmol/m/s/MPa
          smp_mpa(j) = smp_l(c,j) * 1.e-03_r8 * head                ! mm -> m -> MPa

          ! Root biomass density: g biomass / m3 soil

          root_biomass_density = root_biomass(p) * rootfr(p,j) / dz(c,j)
          root_biomass_density = max(root_biomass_density, 1.e-10_r8)

          ! Root length density: m root per m3 soil

          root_length_density = root_biomass_density / (root_density_SPA(patch%itype(p)) * root_cross_sec_area)

          ! Distance between roots: m

          root_dist = sqrt (1._r8  / (root_length_density * pi))

          ! Soil-to-root resistance (MPa.s.m2/mmol H2O)

          soilr1 = log(root_dist/root_radius_SPA(patch%itype(p))) / (2._r8 * pi * root_length_density * dz(c,j) * hk)

          ! Root-to-stem resistance (MPa.s.m2/mmol H2O)

          soilr2 = root_resist_SPA(patch%itype(p)) / (root_biomass_density * dz(c,j))

          ! Belowground resistance (MPa.s.m2/mmol H2O) 

          soilr = soilr1 + soilr2

          ! Total belowground resistance. Sum the conductances (1/soilr) for
          ! each soil layer and then convert back to a resistance after the
          ! summation.

          rsoil(p) = rsoil(p) + 1._r8 / soilr

          ! Maximum transpiration for each layer (mmol H2O/m2/s). No negative
          ! transpiration and no transpiration from frozen soil.

          evap(j) = (smp_mpa(j) - minlwp_SPA(patch%itype(p))) / soilr
          evap(j) = max (evap(j), 0._r8)
          if (h2osoi_ice(c,j) > 0._r8) evap(j) = 0._r8

       end do

       ! Belowground resistance: resistance = 1 / conductance

       rsoil(p) = lai(p) / rsoil(p)

       ! Weighted soil water potential (MPa) and fractional uptake from soil layers

       totevap = sum(evap)
       psis(p) = 0._r8
       soil_et_loss(p,:) = 0._r8

       do j = 1, nlevsoi
          psis(p) = psis(p) + smp_mpa(j) * evap(j)
          if (totevap > 0._r8) then
             soil_et_loss(p,j) = evap(j) / totevap
          else
             soil_et_loss(p,j) = 1._r8 / nlevsoi
          end if
       end do

       if (totevap > 0._r8) then
          psis(p) = psis(p) / totevap
       else
          psis(p) = minlwp_SPA(patch%itype(p))
       end if

    end do

    end associate
  end subroutine SoilResistance

end module MLPlantHydraulicsMod
