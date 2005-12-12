#include <misc.h>
#include <preproc.h>

module DGVMMortalityMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: MortalityMod
!
! !DESCRIPTION:
! Tree background and stress mortality
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
  public :: Mortality
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
! !IROUTINE: Mortality
!
! !INTERFACE:
  subroutine Mortality(lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Tree background and stress mortality
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
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. mortality)
!
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ivt(:)          ! pft vegetation type
    real(r8), pointer :: bm_inc(:)       ! biomass increment
    real(r8), pointer :: lm_ind(:)       ! individual leaf mass
    real(r8), pointer :: sm_ind(:)       ! individual sapwood mass
    real(r8), pointer :: hm_ind(:)       ! individual heartwood mass
    real(r8), pointer :: rm_ind(:)       ! individual root mass
    real(r8), pointer :: agddtw(:)       ! accumulated growing degree days above twmax
    real(r8), pointer :: turnover_ind(:) !
    logical , pointer :: tree(:)         ! ecophys const - whether this pft is a tree type
    real(r8), pointer :: sla(:)          ! ecophys const - specific leaf area [m2 leaf g-1 carbon]o
!
! local pointers to implicit inout arguments
!
    logical , pointer :: present(:)      ! whether PFT present in patch
    real(r8), pointer :: nind(:)         ! number of individuals
    real(r8), pointer :: litterag(:)     ! above ground litter
    real(r8), pointer :: litterbg(:)     ! below ground litter
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8), parameter :: k_mort = 0.3_r8 !coefficient of growth efficiency in mortality equation
    real(r8), parameter :: ramp_agddtw = 300.0_r8
    integer  :: p,fp      !index
    real(r8) :: mort_max  ! asymptotic maximum mortality rate (/yr)
    real(r8) :: bm_delta  ! net individual living biomass increment
                          ! (incorporating loss through leaf, root and sapwood turnover) (gC)
    real(r8) :: mort      ! tree mortality rate
    real(r8) :: nind_kill ! reduction in individual density due to mortality (indiv/m2)
    real(r8) :: greffic
    real(r8) :: heatstress
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type scalar members

    ivt          => clm3%g%l%c%p%itype
    bm_inc       => clm3%g%l%c%p%pdgvs%bm_inc
    lm_ind       => clm3%g%l%c%p%pdgvs%lm_ind
    sm_ind       => clm3%g%l%c%p%pdgvs%sm_ind
    hm_ind       => clm3%g%l%c%p%pdgvs%hm_ind
    rm_ind       => clm3%g%l%c%p%pdgvs%rm_ind
    agddtw       => clm3%g%l%c%p%pdgvs%agddtw
    turnover_ind => clm3%g%l%c%p%pdgvs%turnover_ind
    present      => clm3%g%l%c%p%pdgvs%present
    nind         => clm3%g%l%c%p%pdgvs%nind
    litterag     => clm3%g%l%c%p%pdgvs%litterag
    litterbg     => clm3%g%l%c%p%pdgvs%litterbg
    tree         => dgv_pftcon%tree
    sla          => pftcon%sla

    ! Compute stress mortality

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       if (tree(ivt(p))) then

          if (ivt(p)==3 .or. ivt(p) == 8) then
             mort_max = 0.03_r8 !testing diff values
          else
             mort_max = 0.01_r8 !original value for all pfts
          end if

          heatstress = min(1.0_r8, agddtw(p) / ramp_agddtw)

          ! Calculate net individual living biomass increment
          bm_delta = max(0.0_r8, bm_inc(p) / nind(p) - turnover_ind(p))

          ! Calculate growth efficiency (net biomass increment per unit leaf area)

          !BUGFIX
          if (lm_ind(p) == 0) then
             greffic = 0._r8
          else
             greffic = bm_delta / lm_ind(p) / sla(ivt(p))
          end if
          !BUGFIX

          ! Mortality rate inversely related to growth efficiency (Prentice et al 1993)
          mort = mort_max / (1.0_r8 + k_mort * greffic)

          ! Reduce individual density (=> gridcell-level biomass) by mortality rate
          mort = min(1.0_r8, mort + heatstress)
          nind_kill = nind(p) * mort
          nind(p) = nind(p) - nind_kill

          ! Transfer lost biomass to litter
          litterag(p) = litterag(p) + nind_kill * (lm_ind(p) + sm_ind(p) + hm_ind(p))
          litterbg(p) = litterbg(p) + nind_kill * rm_ind(p)
       end if
       if (nind(p) == 0.0_r8) present(p) = .false.
    end do

  end subroutine Mortality

#endif

end module DGVMMortalityMod
