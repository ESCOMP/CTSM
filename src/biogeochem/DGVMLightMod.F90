#include <misc.h>
#include <preproc.h>

module DGVMLightMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: LightMod
!
! !DESCRIPTION:
! Calculate light competition
! Update fpc (for establishment routine)
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
  public :: Light
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
! !IROUTINE: Light
!
! !INTERFACE:
  subroutine Light(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Calculate light competition
! Update fpc (for establishment routine)
! Called once per year
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg                  ! gridcell bounds
    integer, intent(in) :: lbp, ubp                  ! pft bounds
    integer, intent(in) :: num_natvegp               ! number of naturally-vegetated pfts in filter
    integer, intent(in) :: filter_natvegp(ubp-lbp+1) ! pft filter for naturally-vegetated points
!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine light)
! 3/4/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ivt(:)       ! pft vegetation type
    integer , pointer :: pgridcell(:) ! gridcell index of corresponding pft
    real(r8), pointer :: fpcinc(:)    ! foliar projective cover increment (fraction)
    real(r8), pointer :: sm_ind(:)    ! individual stem mass
    real(r8), pointer :: hm_ind(:)    ! individual heartwood mass
    real(r8), pointer :: crownarea(:) ! area that each individual tree takes up (m^2)
    real(r8), pointer :: sla(:)       ! ecophys const - specific leaf area [m2 leaf g-1 carbon]
    logical , pointer :: tree(:)      ! ecophys const - whether this pft is a tree type
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: fpcgrid(:)   ! foliar projective cover on gridcell (fraction)
    real(r8), pointer :: nind(:)      ! number of individuals
    real(r8), pointer :: litterag(:)  ! above ground litter
    real(r8), pointer :: litterbg(:)  ! below ground litter
    real(r8), pointer :: lm_ind(:)    ! individual leaf mass
    real(r8), pointer :: rm_ind(:)    ! individual root mass
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    real(r8), parameter :: fpc_tree_max = 0.95_r8  !maximum total tree FPC
    integer  :: p,fp, g                         ! indices
    real(r8) :: fpc_tree_total(lbg:ubg)
    real(r8) :: fpc_inc_tree(lbg:ubg)
    real(r8) :: fpc_grass_total(lbg:ubg)
    integer  :: ntree(lbg:ubg)
    real(r8) :: excess
    real(r8) :: nind_kill
    real(r8) :: lm_old
    real(r8) :: lm_kill
    real(r8) :: rm_kill
    real(r8) :: lai_ind
    real(r8) :: fpc_ind

!-----------------------------------------------------------------------

    ! Assign local pointers to derived type scalar members

    ivt       => clm3%g%l%c%p%itype
    pgridcell => clm3%g%l%c%p%gridcell
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    fpcinc    => clm3%g%l%c%p%pdgvs%fpcinc
    nind      => clm3%g%l%c%p%pdgvs%nind
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    fpcinc    => clm3%g%l%c%p%pdgvs%fpcinc
    litterag  => clm3%g%l%c%p%pdgvs%litterag
    litterbg  => clm3%g%l%c%p%pdgvs%litterbg
    lm_ind    => clm3%g%l%c%p%pdgvs%lm_ind
    sm_ind    => clm3%g%l%c%p%pdgvs%sm_ind
    hm_ind    => clm3%g%l%c%p%pdgvs%hm_ind
    rm_ind    => clm3%g%l%c%p%pdgvs%rm_ind
    crownarea => clm3%g%l%c%p%pdgvs%crownarea
    tree      => dgv_pftcon%tree
    sla       => pftcon%sla

    ! Initialize gridcell-level metrics

!dir$ concurrent
!cdir nodep
    do g = lbg, ubg
       fpc_tree_total(g) = 0._r8
       fpc_inc_tree(g) = 0._r8
       fpc_grass_total(g) = 0._r8
       ntree(g) = 0
    end do

    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       g = pgridcell(p)
       if (tree(ivt(p))) then
          ntree(g) = ntree(g) + 1
          fpc_tree_total(g) = fpc_tree_total(g) + fpcgrid(p)
          fpc_inc_tree(g) = fpc_inc_tree(g) + fpcinc(p)
       else    ! if grass
          fpc_grass_total(g) = fpc_grass_total(g) + fpcgrid(p)
       end if
    end do

    ! The gridcell level metrics are now in place, continue...

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       g = pgridcell(p)

       ! light competition

       if (tree(ivt(p))) then

          if (fpc_tree_total(g) > fpc_tree_max) then    ! case (1)

             if (fpc_inc_tree(g) > 0.0_r8) then
                excess = (fpc_tree_total(g) - fpc_tree_max) * &
                     fpcinc(p) / fpc_inc_tree(g)
             else
                excess = (fpc_tree_total(g) - fpc_tree_max) / &
                     real(ntree(g))
             end if

             ! Reduce individual density (and thereby gridcell-level biomass)
             ! so that total tree FPC reduced to 'fpc_tree_max'

             nind_kill = nind(p) * excess / fpcgrid(p)
             nind(p) = nind(p) - nind_kill

             ! Transfer lost biomass to litter

             litterag(p) = litterag(p) + nind_kill * (lm_ind(p) + sm_ind(p) + hm_ind(p))
             litterbg(p) = litterbg(p) + nind_kill * rm_ind(p)

          end if

       else   ! if grass

          if (fpc_grass_total(g) > (1.0_r8-min(fpc_tree_total(g), fpc_tree_max))) then

             ! grass competes with itself if total fpc exceeds 1 (**add comments**)

             excess = (min(fpc_tree_total(g), fpc_tree_max) + &
                  fpc_grass_total(g) - 1.0_r8) * fpcgrid(p) / fpc_grass_total(g)
             lm_old = lm_ind(p)
             lm_ind(p) = -2.0_r8 * log(1.0_r8-(fpcgrid(p) - excess)) / sla(ivt(p))
             lm_kill = lm_old - lm_ind(p)
             rm_kill = rm_ind(p) * lm_kill/lm_old
             rm_ind(p) = rm_ind(p) - rm_kill

             ! Transfer lost biomass to litter

             litterag(p) = litterag(p) + lm_kill
             litterbg(p) = litterbg(p) + rm_kill

          end if

       end if   ! end of if-tree

       ! update fpc (for establishment routine) (slevis: lai_ind is local here)

       if (crownarea(p) > 0.0_r8) then
          lai_ind = lm_ind(p) * sla(ivt(p)) / crownarea(p)
       else
          lai_ind = 0.0_r8
       end if
       fpc_ind = 1.0_r8 - exp(-0.5_r8*lai_ind)
       fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind

    end do

  end subroutine Light

#endif

end module DGVMLightMod
