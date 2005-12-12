#include <misc.h>
#include <preproc.h>

module DGVMAllocationMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: DGVMAllocationMod
!
! !DESCRIPTION:
! Performs yearly allocation calculation
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Allocation
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
! !IROUTINE: Allocation
!
! !INTERFACE:
  subroutine Allocation (lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Performs yearly allocation calculation
! Allocation of this year's biomass increment (bm_inc_ind) to the
! three living carbon pools, such that the basic allometric
! relationships (A-C below) are always satisfied.
! -------------------------------------------------------------------
! TREE ALLOCATION
! (A) (leaf area) = latosa * (sapwood xs area)
!       (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
! (B) (leaf mass) = lmtorm * (root mass)
! (C) height = allom2 * (stem diameter)**allom3 (source?)
! (D) (crown area) = min (allom1 * (stem diameter)**reinickerp, crownarea_max)
!
! Mathematical derivation:
!
!   (1) bm_inc_ind = lminc_ind + sminc_ind + rminc_ind
!   (2) leaf_area_new = latosa * sap_xsa_new   [from (A)]
!   (3) leaf_area_new = (lm_ind + lminc_ind) * sla
! from (2) & (3),
!   (4) (lm_ind + lminc_ind) * sla = latosa * sap_xsa_new
! from (4),
!   (5) sap_xsa_new = (lm_ind + lminc_ind) * sla / latosa
!   (6) (lm_ind + lminc_ind) = lmtorm * (rm_ind + rminc_ind) [from (B)]
!   (7) height_new = allom2 * stemdiam_new**allom3  [from (C)]
! from (1),
!   (8) sminc_ind = bm_inc_ind - lminc_ind - rminc_ind
! from (6),
!   (9) rminc_ind=((lm_ind + lminc_ind) / lmtorm) - rm_ind
! from (8) & (9),
!  (10) sminc_ind = bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind)  / lmtorm) + rm_ind
!  (11) wooddens = (sm_ind + sminc_ind + hm_ind) / stemvolume_new
!  (12) stemvolume_new = height_new * pi * stemdiam_new**2 / 4
! from (10), (11) & (12)
!  (13) stemdiam_new = [ ((sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
!         / wooddens) / (height_new * pi / 4) ]**(1/2)
! combining (7) and (13),
!  (14) height_new = allom2 * [ ((sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
!         / wooddens) / (height_new * pi / 4) ]**(1/2 * allom3)
! from (14),
!  (15) height_new**(1 + 2 / allom3) = allom2**(2 / allom3)
!         * ((sm_ind + bm_inc_ind - lminc_ind - ((lm_ind + lminc_ind)
!         / lmtorm) + rm_ind + hm_ind) / wooddens) / (pi / 4)
!  (16) wooddens = (sm_ind + sminc_ind) / sapvolume_new
! from (10) and (16),
!  (17) wooddens = (sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind) / sapvolume_new
!  (18) sapvolume_new = height_new * sap_xsa_new
! from (17) and (18),
!  (19) sap_xsa_new = (sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind)
!         / (height_new * wooddens)
! from (19),
!  (20) height_new = (sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
!         / (sap_xsa_new * wooddens)
! from (5) and (20),
!  (21) height_new**(1 + 2 / allom3) = [ (sm_ind + bm_inc_ind
!         - lminc_ind - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
!         / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
!         **(1 + 2 / allom3)
! -------------------------------------------------------------------
!  (15) and (21) are two alternative expressions for
!       height_new**(1 + 2 / allom3). Combining these,

!  (22) allom2**(2 / allom3) * ((sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind + hm_ind)
!         / wooddens) / (pi / 4) - [ (sm_ind + bm_inc_ind - lminc_ind
!         - ((lm_ind + lminc_ind) / lmtorm) + rm_ind )
!         / ((lm_ind + lminc_ind) * sla * wooddens / latosa) ]
!         **(1 + 2 / allom3)
!         = 0
!
! Equation (22) can be expressed in the form f(lminc_ind)=0.
!
! Numerical methods are used to solve the equation for the
! unknown lminc_ind.
! -------------------------------------------------------------------
!
! Work out minimum leaf production to maintain current sapmass
!
!  (23) sap_xsa = sm_ind / wooddens / height
! from (A) and (23),
!  (24) leaf_mass * sla = latosa * sap_mass / wooddens / height
! from (24),
!  (25) leaf_mass = latosa * sap_mass / (wooddens * height * sla)
! from (25), assuming sminc_ind=0,
!  (26) lm_ind + lminc_ind_min = latosa * sm_ind
!         / (wooddens * height * sla)
! from (26),
!  (27) lminc_ind_min = latosa * sm_ind / (wooddens * height * sla)
!         - lm_ind
! Work out minimum root production to support this leaf mass
! (i.e. lm_ind + lminc_ind_min)
! May be negative following a reduction in soil water limitation
! (increase in lmtorm) relative to last year.

! from (B) and (25),
!  (28) root_mass = latosa * sap_mass / (wooddens * height * sla)
!         / lmtorm
! from (28), assuming sminc_ind=0,
!  (29) rm_ind + rminc_ind_min = latosa * sm_ind
!         / (wooddens * height * sla * lmtorm)
! from (29),
!  (30) rminc_ind_min = latosa * sm_ind
!         / (wooddens * height * sla * lmtorm) - rm_ind
! -------------------------------------------------------------------
!
! Attempt to distribute this year's production among leaves and roots only
!  (31) bm_inc_ind = lminc_ind + rminc_ind
! from (31) and (9),
!  (32) bm_inc_ind = lminc_ind + ((lm_ind + lminc_ind) / lmtorm)
!         - rm_ind
! from (32)
!  (33) lminc_ind = (bm_inc_ind - lm_ind / lmtorm + rm_ind) /
!         (1 + 1 / lmtorm)
!
! -------------------------------------------------------------------
!
! from (25),
!  (34) lm_ind + lminc_ind = latosa * (sm_ind + sminc_ind)
!         / (wooddens * height * sla)
! from (34),
!  (35) sminc_ind = (lm_ind + lminc_ind) * wooddens * height * sla
!         / latosa - sm_ind
! -------------------------------------------------------------------
!
! !USES:
    use clmtype
    use shr_const_mod, ONLY: SHR_CONST_PI
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
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. allocation)
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ivt(:)             ! pft vegetation type
    real(r8), pointer :: sla(:)             ! ecophys const - specific leaf area [m2 leaf g-1 carbon]
    logical , pointer :: tree(:)            ! ecophys const - whether this pft is a tree type
    real(r8), pointer :: allom1(:)          ! ecophys const - parameter in allometric
    real(r8), pointer :: allom2(:)          ! ecophys const - parameter in allometric
    real(r8), pointer :: allom3(:)          ! ecophys const - parameter in allometric
    real(r8), pointer :: latosa(:)          ! ecophys const - ratio of leaf area to sapwood cross-sectional area
    real(r8), pointer :: wooddens(:)        ! ecophys const - wood density (gC/m3)
    real(r8), pointer :: reinickerp(:)      ! ecophys const - parameter in allometric equation
    real(r8), pointer :: crownarea_max(:)   ! ecophys const - tree maximum crown area [m2]
    real(r8), pointer :: init_lmtorm(:)     ! ecophys const - leaf:root ratio under non-water stressed conditions
    real(r8), pointer :: bm_inc(:)          ! biomass increment
    real(r8), pointer :: nind(:)            ! number of individuals (#/m**2)
    real(r8), pointer :: annpsn(:)          ! annual photosynthesis (umol CO2 /m**2)
    real(r8), pointer :: annpsnpot(:)       ! annual potential photosynthesis (same units)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: fpc_grid(:)         ! foliar projective cover on gridcell (fraction)
    real(r8), pointer :: crownarea(:)        ! area that each individual tree takes up (m^2)
    real(r8), pointer :: height(:)           ! canopy top (m)
    real(r8), pointer :: lm_ind(:)           ! individual leaf mass
    real(r8), pointer :: sm_ind(:)           ! individual sapwood mass
    real(r8), pointer :: hm_ind(:)           ! individual heartwood mass
    real(r8), pointer :: rm_ind(:)           ! individual root mass
    real(r8), pointer :: litter_ag(:)        ! above ground litter
    real(r8), pointer :: litter_bg(:)        ! below ground litter
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: lai_ind(:)          ! LAI per individual
    real(r8), pointer :: fpc_inc(:)          ! foliar projective cover increment (fraction)
!
!EOP
!
! !LOCAL VARIABLES:
!
    integer , parameter :: nseg = 20
    real(r8), parameter :: xacc = 0.1_r8         ! threshold x-axis and threshold
    real(r8), parameter :: yacc = 1.0e-10_r8     ! y-axis precision of allocation soln
    real(r8), parameter :: pi = SHR_CONST_PI  ! 3.14159...
    integer  :: p, fp                         ! index
    integer  :: fn,fnold                      ! number of elements in local filter
    integer  :: filterp(num_natvegp)          ! local pft filter of filter_natveg
    real(r8) :: wscal
    real(r8) :: lmtorm(lbp:ubp)
    real(r8) :: bm_inc_ind(lbp:ubp)
    real(r8) :: lminc_ind_min(lbp:ubp)
    real(r8) :: rminc_ind_min(lbp:ubp)
    real(r8) :: lminc_ind(lbp:ubp)
    real(r8) :: rminc_ind(lbp:ubp)
    real(r8) :: sminc_ind(lbp:ubp)
    real(r8) :: fpc_ind(lbp:ubp)
    real(r8) :: fpc_grid_old(lbp:ubp)
    real(r8) :: x1(lbp:ubp), x2(lbp:ubp), dx(lbp:ubp)
    real(r8) :: sap_xsa(lbp:ubp)
    real(r8) :: stemdiam(lbp:ubp)
    real(r8) :: fx1(lbp:ubp)
    real(r8) :: fmid(lbp:ubp)
    real(r8) :: xmid(lbp:ubp)
    real(r8) :: sign(lbp:ubp)
    real(r8) :: rtbis(lbp:ubp)
!-----------------------------------------------------------------------

   ! Assign local pointers to derived type members (pft-level)

    ivt           => clm3%g%l%c%p%itype
    height        => clm3%g%l%c%p%pps%htop
    fpc_grid      => clm3%g%l%c%p%pdgvs%fpcgrid
    fpc_inc       => clm3%g%l%c%p%pdgvs%fpcinc
    bm_inc        => clm3%g%l%c%p%pdgvs%bm_inc
    nind          => clm3%g%l%c%p%pdgvs%nind
    crownarea     => clm3%g%l%c%p%pdgvs%crownarea
    lm_ind        => clm3%g%l%c%p%pdgvs%lm_ind
    sm_ind        => clm3%g%l%c%p%pdgvs%sm_ind
    hm_ind        => clm3%g%l%c%p%pdgvs%hm_ind
    rm_ind        => clm3%g%l%c%p%pdgvs%rm_ind
    lai_ind       => clm3%g%l%c%p%pdgvs%lai_ind
    litter_ag     => clm3%g%l%c%p%pdgvs%litterag
    litter_bg     => clm3%g%l%c%p%pdgvs%litterbg
    annpsnpot     => clm3%g%l%c%p%pdgvs%annpsnpot
    annpsn        => clm3%g%l%c%p%pdgvs%annpsn
    tree          => dgv_pftcon%tree
    allom1        => dgv_pftcon%allom1
    allom2        => dgv_pftcon%allom2
    allom3        => dgv_pftcon%allom3
    latosa        => dgv_pftcon%latosa
    wooddens      => dgv_pftcon%wooddens
    reinickerp    => dgv_pftcon%reinickerp
    crownarea_max => dgv_pftcon%crownarea_max
    init_lmtorm   => dgv_pftcon%lmtorm
    sla           => pftcon%sla

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       bm_inc_ind(p) = bm_inc(p) / nind(p)
    end do

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       if (tree(ivt(p))) then

          ! calculate this year's leaf to fine root mass ratio from mean annual
          ! water scalar and pft specific parameter
          ! slevis: in lpj wscal=awscal(pft)/aleafdays(pft), awscal=SUM(dwscal) (dphen>0)
          !         dwscal=min(supply(pft)/demandpot(pft),1) or =1 when demand(pft)=0 etc
          !         here wscal=annpsn/annpsnpot

          if (annpsnpot(p) > 0.0_r8) then
             wscal = annpsn(p)/annpsnpot(p)
          else
             wscal = 1.0_r8
          end if
          lmtorm(p) = init_lmtorm(ivt(p)) * wscal

          ! tree allocation

          lminc_ind_min(p) = latosa(ivt(p))*sm_ind(p)/(wooddens(ivt(p))*height(p)*sla(ivt(p))) - lm_ind(p)           ! eqn(27)
          rminc_ind_min(p) = latosa(ivt(p))*sm_ind(p)/(wooddens(ivt(p))*height(p)*sla(ivt(p)))*lmtorm(p) - rm_ind(p) ! eqn(30)

          if (rminc_ind_min(p) > 0.0_r8 .and. lminc_ind_min(p) > 0.0_r8 .and. &
               rminc_ind_min(p) + lminc_ind_min(p) <= bm_inc_ind(p)) then

             ! Normal allocation (positive increment to all living C compartments)

             ! Calculation of leaf mass increment (lminc_ind) that satisfies
             ! Eqn (22) using Bisection Method (Press et al 1986, p 346)

             ! Seeking a root for non-negative lminc_ind, rminc_ind and
             ! sminc_ind.  There should be exactly one (no proof presented, but
             ! Steve has managed one) and it should lie between x1=0 and
             ! x2=(bm_inc_ind-(lm_ind/lmtorm-rm_ind))/(1+1/lmtorm).

             x1(p) = 0.0_r8
             x2(p) = (bm_inc_ind(p) - (lm_ind(p)/lmtorm(p) - rm_ind(p))) / (1.0_r8 + 1.0_r8 / lmtorm(p))
             dx(p) = (x2(p)-x1(p)) / real(nseg)

             if (lm_ind(p) == 0.0_r8) x1(p) = x1(p) + dx(p) !to avoid division by zero

             ! evaluate f(x1)=LHS of eqn (22) at x1

             fx1(p) = allom2(ivt(p)) ** (2.0_r8/allom3(ivt(p))) *                             &
                  (sm_ind(p) + bm_inc_ind(p) - x1(p) - (lm_ind(p) + x1(p))/lmtorm(p) + rm_ind(p) + hm_ind(p)) / &
                  wooddens(ivt(p)) / (pi/4.0_r8) -                                                 &
                  ( (sm_ind(p) + bm_inc_ind(p) - x1(p) - (lm_ind(p) + x1(p))/lmtorm(p) + rm_ind(p)) /        &
                  ((lm_ind(p) + x1(p))*sla(ivt(p))*wooddens(ivt(p))/latosa(ivt(p))) )**(1.0_r8+2.0_r8/allom3(ivt(p)))

             ! Find approximate location of leftmost root on the interval
             ! (x1,x2).  Subdivide (x1,x2) into nseg equal segments seeking
             ! change in sign of f(xmid) relative to f(x1).

             fmid(p)= fx1(p)
             xmid(p) = x1(p)
          end if

       end if
    end do

    fn = 0
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       if (tree(ivt(p))) then
          if (rminc_ind_min(p) > 0.0_r8 .and. lminc_ind_min(p) > 0.0_r8 .and. &
              rminc_ind_min(p) + lminc_ind_min(p) <= bm_inc_ind(p)) then
             if (fmid(p)*fx1(p) > 0.0_r8 .and. xmid(p) < x2(p)) then
                fn = fn + 1
                filterp(fn) = p
             end if
          end if
       end if
    end do
    do while (fn > 0)
!dir$ concurrent
!cdir nodep
       do fp = 1,fn
          p = filterp(fp)
          xmid(p) = xmid(p) + dx(p)
          fmid(p) = allom2(ivt(p)) ** (2.0_r8/allom3(ivt(p))) * &
                    (sm_ind(p) + bm_inc_ind(p) - xmid(p) - (lm_ind(p)+xmid(p))/lmtorm(p) + rm_ind(p) + hm_ind(p)) / &
                    wooddens(ivt(p)) / (pi/4.0_r8) - &
                    ( (sm_ind(p) + bm_inc_ind(p) - xmid(p) - (lm_ind(p)+xmid(p))/lmtorm(p) + rm_ind(p)) / &
                    ((lm_ind(p) + xmid(p))*sla(ivt(p))*wooddens(ivt(p))/latosa(ivt(p))) )**(1.0_r8+2.0_r8/allom3(ivt(p)))
       end do
       fnold = fn
       fn = 0
       do fp = 1,fnold
          p = filterp(fp)
          if (fmid(p)*fx1(p) > 0.0_r8 .and. xmid(p) < x2(p)) then
             fn = fn + 1
             filterp(fn) = p
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       if (tree(ivt(p))) then
          if (rminc_ind_min(p) > 0.0_r8 .and. lminc_ind_min(p) > 0.0_r8 .and. &
               rminc_ind_min(p) + lminc_ind_min(p) <= bm_inc_ind(p)) then

             x1(p) = xmid(p) - dx(p)
             x2(p) = xmid(p)

             ! Apply bisection method to find root on the new interval (x1,x2)

             fx1(p) = allom2(ivt(p)) ** (2.0_r8/allom3(ivt(p))) *                                       &
                  (sm_ind(p) + bm_inc_ind(p) - x1(p) - (lm_ind(p) + x1(p))/lmtorm(p) + rm_ind(p) + hm_ind(p)) / &
                  wooddens(ivt(p)) / (pi/4.0_r8) -                                                 &
                  ( (sm_ind(p) + bm_inc_ind(p) - x1(p) - (lm_ind(p) + x1(p))/lmtorm(p) + rm_ind(p)) /        &
                  ((lm_ind(p) + x1(p))*sla(ivt(p))*wooddens(ivt(p))/latosa(ivt(p))) )**(1.0_r8+2.0_r8/allom3(ivt(p)))

             if (fx1(p) >= 0.0_r8) then
                sign(p) = -1.0_r8
             else
                sign(p) = 1.0_r8
             end if

             rtbis(p) = x1(p)
             dx(p) = x2(p)-x1(p)

             ! Bisection loop
             ! Search iterates on value of xmid until xmid lies within
             ! xacc of the root, i.e. until |xmid-x|<xacc where f(x)=0

             fmid(p) = 1.0_r8  !dummy value to guarantee entry to next do-while loop

          end if
       end if
    end do

    fn = 0
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       if (tree(ivt(p))) then
          if (rminc_ind_min(p) > 0.0_r8 .and. lminc_ind_min(p) > 0.0_r8 .and. &
              rminc_ind_min(p) + lminc_ind_min(p) <= bm_inc_ind(p)) then
             if (dx(p) >= xacc .and. abs(fmid(p)) > yacc) then
                fn = fn + 1
                filterp(fn) = p
             end if
          end if
       end if
    end do
    do while (fn > 0)
!dir$ concurrent
!cdir nodep
       do fp = 1,fn
          p = filterp(fp)
          dx(p) = dx(p)*0.5_r8
          xmid(p) = rtbis(p)+dx(p)

          ! calculate fmid=f(xmid) [eqn (22)]

          fmid(p) = allom2(ivt(p)) ** (2.0_r8/allom3(ivt(p))) *  &
                    (sm_ind(p) + bm_inc_ind(p) - xmid(p) - (lm_ind(p)+xmid(p))/lmtorm(p) + rm_ind(p) + hm_ind(p)) / &
                    wooddens(ivt(p)) / (pi/4.0_r8) - &
                    ( (sm_ind(p) + bm_inc_ind(p) - xmid(p) - (lm_ind(p)+xmid(p))/lmtorm(p) + rm_ind(p)) / &
                    ((lm_ind(p) + xmid(p))*sla(ivt(p))*wooddens(ivt(p))/latosa(ivt(p))) )**(1.0_r8+2.0_r8/allom3(ivt(p)))

          if (fmid(p)*sign(p) <= 0.0_r8) rtbis(p) = xmid(p)
       end do
       fnold = fn
       fn = 0
       do fp = 1,fnold
          p = filterp(fp)
          if (dx(p) >= xacc .and. abs(fmid(p)) > yacc) then
             fn = fn + 1
             filterp(fn) = p
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       if (tree(ivt(p))) then
          if (rminc_ind_min(p) > 0.0_r8 .and. lminc_ind_min(p) > 0.0_r8 .and. &
               rminc_ind_min(p) + lminc_ind_min(p) <= bm_inc_ind(p)) then

             ! Now rtbis contains numerical solution for lminc_ind given eqn (22)

             lminc_ind(p) = rtbis(p)

             ! Calculate increments in other compartments using allometry relationships

             rminc_ind(p) = (lm_ind(p) + lminc_ind(p)) / lmtorm(p) - rm_ind(p) !eqn (9)
             sminc_ind(p) = bm_inc_ind(p) - lminc_ind(p) - rminc_ind(p)     !eqn (1)

          else

             ! Abnormal allocation: reduction in some C compartment(s) to satisfy allometry

             lminc_ind(p) = (bm_inc_ind(p) - lm_ind(p)/lmtorm(p) + rm_ind(p)) / &
                  (1.0_r8 + 1.0_r8/lmtorm(p)) !eqn (33)

             if (lminc_ind(p) >= 0.0_r8) then

                ! Positive allocation to leafmass

                rminc_ind(p) = bm_inc_ind(p) - lminc_ind(p)  !eqn (31)

                ! Add killed roots (if any) to below-ground litter

                if (rminc_ind(p) < 0.0_r8) then
                   lminc_ind(p) = bm_inc_ind(p)
                   rminc_ind(p) = (lm_ind(p) + lminc_ind(p)) / lmtorm(p) - rm_ind(p)
                   litter_bg(p) = litter_bg(p) - rminc_ind(p) * nind(p)
                end if

             else

                ! Negative allocation to leaf mass

                rminc_ind(p) = bm_inc_ind(p)
                lminc_ind(p) = (rm_ind(p) + rminc_ind(p)) * lmtorm(p) - lm_ind(p) !from eqn (9)

                ! Add killed leaves to litter

                litter_ag(p) = litter_ag(p) - lminc_ind(p) * nind(p)

             end if

             ! Calculate sminc_ind (must be negative)

             sminc_ind(p) = (lm_ind(p) + lminc_ind(p))*wooddens(ivt(p))*height(p)*sla(ivt(p)) / latosa(ivt(p)) - &
                  sm_ind(p) !eqn (35)

             ! Convert killed sapwood to heartwood

             hm_ind(p) = hm_ind(p) - sminc_ind(p)

          end if

          ! Increment C compartments

          lm_ind(p) = lm_ind(p) + lminc_ind(p)
          rm_ind(p) = rm_ind(p) + rminc_ind(p)
          sm_ind(p) = sm_ind(p) + sminc_ind(p)

          ! Calculate new height, diameter and crown area

          sap_xsa(p) = lm_ind(p) * sla(ivt(p)) / latosa(ivt(p))  !eqn (5)

          !BUGFIX
          if (lm_ind(p) == 0) then
             height(p) = 0._r8
          else
             height(p) = sm_ind(p) / sap_xsa(p) / wooddens(ivt(p))
          end if
          !BUGFIX

          stemdiam(p) = (height(p)/allom2(ivt(p))) ** (1.0_r8/allom3(ivt(p))) !eqn (C)
          crownarea(p) = min(allom1(ivt(p))*stemdiam(p)**reinickerp(ivt(p)), crownarea_max(ivt(p))) !eqn (D)

       else !grasses

          lmtorm(p) = init_lmtorm(ivt(p)) !slevis: eliminate influence from soil H2O

          ! GRASS ALLOCATION
          ! Distribute this year's production among leaves and fine roots
          ! according to leaf to rootmass ratio [eqn (33)]
          ! Relocation of C from one compartment to the other not allowed:
          ! negative increment in either compartment transferred to litter

          lminc_ind(p) = (bm_inc_ind(p) - lm_ind(p)/lmtorm(p) + rm_ind(p)) / (1.0_r8 + 1.0_r8/lmtorm(p))
          rminc_ind(p) = bm_inc_ind(p) - lminc_ind(p)

          if (lminc_ind(p) >= 0.0_r8) then

             ! Add killed roots (if any) to below-ground litter

             ! CHECK: take out if statement because if rminc is negative than
             ! root mass has been translocated to the leaves, therefore mass balance
             ! problem since this carbon stays in the vegetation but is in addition
             ! added to the litter pool. ALLOW translocation from roots to leaves
             ! i.e. assume carbon stores in the roots which can be delivered
             ! to the leaves under times of stress.

             ! if (rminc_ind < 0.0) litter_bg = litter_bg -rminc_ind * nind

          else

             ! Negative allocation to leaf mass

             rminc_ind(p) = bm_inc_ind(p)
             lminc_ind(p) = (rm_ind(p) + rminc_ind(p))*lmtorm(p) - lm_ind(p) !from eqn (9)

             ! Add killed leaves to litter

             litter_ag(p) = litter_ag(p) - lminc_ind(p) * nind(p)

          end if

          ! Increment C compartments

          lm_ind(p) = lm_ind(p) + lminc_ind(p)
          rm_ind(p) = rm_ind(p) + rminc_ind(p)

       end if   ! end of if-tree (or grass)

       ! Update LAI and FPC

       if (crownarea(p) > 0.0_r8) then
          lai_ind(p) = lm_ind(p) * sla(ivt(p)) / crownarea(p)
       else
          lai_ind(p) = 0.0_r8
       end if

       fpc_ind(p) = 1.0_r8 - exp(-0.5_r8*lai_ind(p))
       fpc_grid_old(p) = fpc_grid(p)
       fpc_grid(p) = crownarea(p) * nind(p) * fpc_ind(p)
       fpc_inc(p) = max(0.0_r8, fpc_grid(p) - fpc_grid_old(p))

       ! diagnostic (slevis)
       ! write(15,*)lminc_ind/bm_inc_ind,rminc_ind/bm_inc_ind,sminc_ind/bm_inc_ind
       ! end diagnostic

    end do

  end subroutine Allocation

#endif

end module DGVMAllocationMod
