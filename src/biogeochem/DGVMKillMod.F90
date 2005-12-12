#include <misc.h>
#include <preproc.h>

module DGVMKillMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: KillMod
!
! !DESCRIPTION:
! Removal of PFTs with negative annual C increment
! NB: PFTs newly beyond their bioclimatic limits are removed in
! subroutine establishment
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
  public :: Kill
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
! !IROUTINE: Kill
!
! !INTERFACE:
  subroutine Kill(lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Removal of PFTs with negative annual C increment
! NB: PFTs newly beyond their bioclimatic limits are removed in
! subroutine establishment
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
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. kill)
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ivt(:)             ! pft vegetation type
    real(r8), pointer :: nind(:)            ! number of individuals (#/m**2)
    real(r8), pointer :: lm_ind(:)          ! individual leaf mass
    real(r8), pointer :: sm_ind(:)          ! individual sapwood mass
    real(r8), pointer :: hm_ind(:)          ! individual heartwood mass
    real(r8), pointer :: rm_ind(:)          ! individual root mass
    real(r8), pointer :: bm_inc(:)          ! biomass increment
    logical , pointer :: tree(:)            ! ecophys const - whether this pft is a tree
!
! local pointers to implicit inout arguments
!
    logical , pointer :: present(:)         ! whether PFT present in patch
    real(r8), pointer :: litter_ag(:)       ! above ground litter
    real(r8), pointer :: litter_bg(:)       ! below ground litter
!EOP
!
! !LOCAL VARIABLES:
!
    integer :: p,fp
!-----------------------------------------------------------------------

   ! Assign local pointers to derived type members (pft-level)

    ivt       => clm3%g%l%c%p%itype
    present   => clm3%g%l%c%p%pdgvs%present
    nind      => clm3%g%l%c%p%pdgvs%nind
    lm_ind    => clm3%g%l%c%p%pdgvs%lm_ind
    sm_ind    => clm3%g%l%c%p%pdgvs%sm_ind
    hm_ind    => clm3%g%l%c%p%pdgvs%hm_ind
    rm_ind    => clm3%g%l%c%p%pdgvs%rm_ind
    bm_inc    => clm3%g%l%c%p%pdgvs%bm_inc
    litter_ag => clm3%g%l%c%p%pdgvs%litterag
    litter_bg => clm3%g%l%c%p%pdgvs%litterbg
    tree      => dgv_pftcon%tree

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)

       if (bm_inc(p) < 0.0_r8) then !negative C increment this year

          present(p) = .false.   !remove PFT

          ! Transfer killed biomass to litter

          if (tree(ivt(p))) then  ! redundant if block? (slevis)
             litter_ag(p) = litter_ag(p) + (lm_ind(p) + sm_ind(p) + hm_ind(p)) * nind(p)
          else         ! if grass
             litter_ag(p) = litter_ag(p) + lm_ind(p) * nind(p)
          end if

          litter_bg(p) = litter_bg(p) + rm_ind(p) * nind(p)

       end if
    end do

  end subroutine Kill

#endif

end module DGVMKillMod
