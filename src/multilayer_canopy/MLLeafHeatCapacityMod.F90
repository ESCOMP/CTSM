module MLLeafHeatCapacityMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf heat capacity
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
  public :: LeafHeatCapacity
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LeafHeatCapacity (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf heat capacity
    !
    ! !USES:
    use clm_varcon, only : cpliq
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : cpbio, fcarbon, fwater
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter        ! Number of patches in filter
    integer, intent(in) :: filter(:)         ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                           ! Filter index
    integer  :: p                            ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                           ! Aboveground layer index
    real(r8) :: lma                          ! Leaf carbon mass per area (kg C / m2 leaf)
    real(r8) :: dry_weight                   ! Leaf dry mass per area (kg DM / m2 leaf)
    real(r8) :: fresh_weight                 ! Leaf fresh mass per area (kg FM / m2 leaf)
    real(r8) :: leaf_water                   ! Leaf water (kg H2O / m2 leaf)
    !---------------------------------------------------------------------

    associate ( &
                                                 ! *** Input ***
    slatop    => pftcon%slatop              , &  ! CLM: Specific leaf area at top of canopy (m2/gC)
    ncan      => mlcanopy_inst%ncan_canopy  , &  ! Number of aboveground layers
    dpai      => mlcanopy_inst%dpai_profile , &  ! Canopy layer plant area index (m2/m2)
                                                 ! *** Output ***
    cpleaf    => mlcanopy_inst%cpleaf_profile &  ! Canopy layer leaf heat capacity (J/m2 leaf/K)
    )

    ! Leaf heat capacity - need to convert specific leaf area (m2/gC) to
    ! leaf mass per area (kgC/m2) and then convert to dry weight (assume
    ! carbon is 50% of dry biomass). Then need to convert dry biomass to
    ! fresh biomass (assume 70% of fresh biomass is water). Then remember
    ! that 70% of fresh biomass is water when calculating heat capacity.

    do fp = 1, num_filter
       p = filter(fp)
       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             lma = 1._r8 / slatop(patch%itype(p)) * 0.001_r8                ! m2 / g C -> kg C / m2
             dry_weight = lma / fcarbon                                     ! kg C / m2 -> kg DM / m2
             fresh_weight = dry_weight / (1._r8 - fwater)                   ! kg DM / m2 -> kg FM / m2
             leaf_water = fwater * fresh_weight                             ! Leaf water (kg H2O / m2 leaf)
             cpleaf(p,ic) = cpbio * dry_weight + cpliq * leaf_water         ! Heat capacity (J/K/m2 leaf) 
          else
             cpleaf(p,ic) = 0._r8
          end if
       end do
    end do

    end associate
  end subroutine LeafHeatCapacity

end module MLLeafHeatCapacityMod
