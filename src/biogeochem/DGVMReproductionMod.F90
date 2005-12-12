#include <misc.h>
#include <preproc.h>

module DGVMReproductionMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ReproductionMod
!
! !DESCRIPTION:
! Deduction of reproduction costs from annual biomass increment
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
  public :: Reproduction
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
! !IROUTINE: Reproduction
!
! !INTERFACE:
  subroutine Reproduction(lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
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
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. reproduction)
!
! !LOCAL VARIABLES:
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: litter_ag(:)        ! above ground litter
    real(r8), pointer :: bm_inc(:)           ! biomass increment
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8), parameter :: reprod_cost = 0.1_r8 ! proportion of NPP lost to reproduction (Harper 1977)
    integer  :: p, fp                        ! pft index
    real(r8) :: reprod                       ! temporary
!-----------------------------------------------------------------------

   ! Assign local pointers to derived type members (pft-level)

    litter_ag => clm3%g%l%c%p%pdgvs%litterag
    bm_inc    => clm3%g%l%c%p%pdgvs%bm_inc

    ! Compute reproduction costs

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)

       ! Calculate allocation to reproduction
       ! Reproduction costs taken simply as a constant fraction of annual NPP

       reprod = max(bm_inc(p) * reprod_cost, 0.0_r8)

       ! assume the costs go to reproductive structures which will
       ! eventually enter the litter pool

       litter_ag(p) = litter_ag(p) + reprod

       ! Reduce biomass increment by reproductive cost

       bm_inc(p) = bm_inc(p) - reprod

    end do

  end subroutine Reproduction

#endif

end module DGVMReproductionMod
