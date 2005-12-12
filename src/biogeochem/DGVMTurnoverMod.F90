#include <misc.h>
#include <preproc.h>

module DGVMTurnoverMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: TurnoverMod
!
! !DESCRIPTION:
! Turnover of PFT-specific fraction from each living C pool
! Leaf and root C transferred to litter, sapwood C to heartwood
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
  public :: Turnover
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
! !IROUTINE: Turnover
!
! !INTERFACE:
  subroutine Turnover(lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Turnover of PFT-specific fraction from each living C pool
! Leaf and root C transferred to litter, sapwood C to heartwood
! Called once per year
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                  ! pft bounds
    integer, intent(in) :: num_natvegp               ! number of naturally-vegetated pfts in filter
    integer, intent(in) :: filter_natvegp(ubp-lbp+1) ! pft filter for naturally-vegetated points

!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. turnover)
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ivt(:)             ! pft vegetation type
    real(r8), pointer :: nind(:)            ! number of individuals (#/m**2)
    real(r8), pointer :: l_turn(:)          ! ecophys const - leaf turnover period [years]
    real(r8), pointer :: s_turn(:)          ! ecophys const - sapwood turnover period [years]
    real(r8), pointer :: r_turn(:)          ! ecophys const - root turnover period [years]
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: litter_ag(:)       ! above ground litter
    real(r8), pointer :: litter_bg(:)       ! below ground litter
    real(r8), pointer :: lm_ind(:)          ! individual leaf mass
    real(r8), pointer :: sm_ind(:)          ! individual sapwood mass
    real(r8), pointer :: hm_ind(:)          ! individual heartwood mass
    real(r8), pointer :: rm_ind(:)          ! individual root mass
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: turnover_ind(:)    !
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: p, fp
    real(r8) :: l_torate
    real(r8) :: s_torate
    real(r8) :: r_torate
    real(r8) :: lm_turn
    real(r8) :: sm_turn
    real(r8) :: rm_turn
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (pft-level)

    ivt          => clm3%g%l%c%p%itype
    nind         => clm3%g%l%c%p%pdgvs%nind
    litter_ag    => clm3%g%l%c%p%pdgvs%litterag
    litter_bg    => clm3%g%l%c%p%pdgvs%litterbg
    lm_ind       => clm3%g%l%c%p%pdgvs%lm_ind
    sm_ind       => clm3%g%l%c%p%pdgvs%sm_ind
    hm_ind       => clm3%g%l%c%p%pdgvs%hm_ind
    rm_ind       => clm3%g%l%c%p%pdgvs%rm_ind
    turnover_ind => clm3%g%l%c%p%pdgvs%turnover_ind
    l_turn       => dgv_pftcon%l_turn
    s_turn       => dgv_pftcon%s_turn
    r_turn       => dgv_pftcon%r_turn

    ! Determine turnover of pft-specific fraction from each living
    ! C pool

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)

       ! Turnover rates are reciprocals of tissue longevity

       l_torate = 1.0_r8 / l_turn(ivt(p))
       s_torate = 1.0_r8 / s_turn(ivt(p))
       r_torate = 1.0_r8 / r_turn(ivt(p))

       ! Calculate the biomass turnover in this year

       lm_turn = lm_ind(p) * l_torate
       sm_turn = sm_ind(p) * s_torate
       rm_turn = rm_ind(p) * r_torate

       ! Update the pools

       lm_ind(p) = lm_ind(p) - lm_turn
       sm_ind(p) = sm_ind(p) - sm_turn
       rm_ind(p) = rm_ind(p) - rm_turn

       ! Convert sapwood to heartwood

       hm_ind(p) = hm_ind(p) + sm_turn

       ! Transfer to litter pools

       litter_ag(p) = litter_ag(p) + lm_turn * nind(p)
       litter_bg(p) = litter_bg(p) + rm_turn * nind(p)

       ! Record total turnover

       turnover_ind(p) = lm_turn + sm_turn + rm_turn

    end do

  end subroutine Turnover

#endif

end module DGVMTurnoverMod
