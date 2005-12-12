#include <misc.h>
#include <preproc.h>

module DGVMFireMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FireMod
!
! !DESCRIPTION: FireMod
! Effect of the fire on vegetation structure and litter
! Called once per year
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Fire
!
! !REVISION HISTORY:
! Module created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fire
!
! !INTERFACE:
  subroutine Fire(lbp, ubp, afire_frac, acflux_fire)

!
! !DESCRIPTION:
! Effect of the fire on vegetation structure and litter
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp            ! pft bounds
    real(r8), intent(out) :: afire_frac(lbp:ubp)
    real(r8), intent(out) :: acflux_fire(lbp:ubp)
!
! !CALLED FROM:
! subroutine EcosystemDyn in module EcosystemdynMod
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine fire)
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ivt(:)             ! pft vegetation type
    real(r8), pointer :: lm_ind(:)          ! individual leaf mass
    real(r8), pointer :: sm_ind(:)          ! individual sapwood mass
    real(r8), pointer :: hm_ind(:)          ! individual heartwood mass
    real(r8), pointer :: rm_ind(:)          ! individual root mass
    real(r8), pointer :: fpc_grid(:)        ! foliar projective cover on gridcell (fraction)
    logical , pointer :: present(:)         ! whether this pft present in patch
    real(r8), pointer :: fire_length(:)     ! fire season in days
    logical , pointer :: tree(:)            ! ecophys const - whether this pft is a tree type
    real(r8), pointer :: resist(:)          ! ecophys const - fire resistance index [units?]
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: litter_ag(:)       ! above ground litter
    real(r8), pointer :: nind(:)            ! number of individuals (#/m**2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    real(r8), parameter :: minfuel = 200.0_r8  ! fuel threshold to carry a fire (gC/m2)
    integer  :: p                           ! index
    real(r8) :: fire_index
    real(r8) :: fire_term
    real(r8) :: disturb
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (pft-level)

    ivt         => clm3%g%l%c%p%itype
    lm_ind      => clm3%g%l%c%p%pdgvs%lm_ind
    sm_ind      => clm3%g%l%c%p%pdgvs%sm_ind
    hm_ind      => clm3%g%l%c%p%pdgvs%hm_ind
    rm_ind      => clm3%g%l%c%p%pdgvs%rm_ind
    fpc_grid    => clm3%g%l%c%p%pdgvs%fpcgrid
    present     => clm3%g%l%c%p%pdgvs%present
    fire_length => clm3%g%l%c%p%pdgvs%firelength
    litter_ag   => clm3%g%l%c%p%pdgvs%litterag
    nind        => clm3%g%l%c%p%pdgvs%nind
    tree        => dgv_pftcon%tree
    resist      => dgv_pftcon%resist

!dir$ concurrent
!cdir nodep
    do p = lbp,ubp
       if (ivt(p) > 0) then

          ! slevis: Orig. had a daily loop to calculate fire_length
          !         Now fire_length comes from subroutine FireSeason

          ! Calculate annual fire index

          fire_index = fire_length(p) / 365.0_r8

          ! Calculate the fraction of the grid cell affected by fire

          fire_term = fire_index - 1.0_r8
          afire_frac(p) = max(fire_index * &
               exp(fire_term / (-0.13_r8*fire_term**3 + 0.6_r8*fire_term**2 + 0.8_r8*fire_term + 0.45_r8)), &
               0.001_r8)

          ! Reduce fraction of patch affected by fire when fuel
          ! becomes limiting (reduced carrying capacity)

          if (litter_ag(p) < minfuel * fpc_grid(p)) afire_frac(p) = 0.001_r8

          ! Implement the effect of the fire on vegetation structure and litter
          ! in the disturbed fraction.

          ! Each PFT is assigned a resistance to fire, representing the fraction of
          ! the PFT which survives a fire. Grasses assumed already to have completed
          ! their life cycle and thus are not affected by fire, giving them
          ! a competitive advantage against woody PFTs.

          if (present(p) .and. tree(ivt(p))) then

             ! Calculate the fraction of individuals in grid cell which die
             ! (slevis: 'in grid cell' because nind is grid average)

             disturb = (1.0_r8-resist(ivt(p))) * afire_frac(p)

             ! Calculate carbon flux to atmosphere (gC/m2) due to burnt biomass

             acflux_fire(p) = disturb * nind(p) * (lm_ind(p) + sm_ind(p) + hm_ind(p) + rm_ind(p))

             ! Update the individual density

             nind(p) = nind(p) * (1.0_r8-disturb)

             ! Add combusted litter to carbon flux to atmosphere term

             acflux_fire(p) = acflux_fire(p) + afire_frac(p) * litter_ag(p)

          else

             acflux_fire(p) = afire_frac(p) * litter_ag(p)

          end if

          ! Update the above ground litter term

          litter_ag(p) = (1.0_r8 - afire_frac(p)) * litter_ag(p)

       else

          afire_frac(p) = 0.0_r8
          acflux_fire(p) = 0.0_r8

       end if
    end do

  end subroutine Fire

#endif

end module DGVMFireMod
