#include <misc.h>
#include <preproc.h>

module DGVMEstablishmentMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: DGVMEstablishmentMod
!
! !DESCRIPTION:
! Calculates establishment of new pfts
! Called once per year
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils, only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Establishment
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
! !IROUTINE: Establishment
!
! !INTERFACE:
  subroutine Establishment(lbp, ubp, lbg, ubg)
!
! !DESCRIPTION:
! Calculates establishment of new pfts
! Called once per year
!
! !USES:
    use clmtype
    use clm_varpar   , only : numpft
    use clm_varcon   , only : istsoil
    use pftvarcon    , only : noveg
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_PI, SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbp, ubp         ! pft bounds
    integer , intent(in) :: lbg, ubg         ! gridcell bounds
!
! !CALLED FROM:
! subroutine lpj in module DGVMMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subr. establishment)
! 3/4/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: plandunit(:)     ! landunit of corresponding pft
    integer , pointer :: pgridcell(:)     ! gridcell of corresponding pft
    integer , pointer :: ltype(:)         ! landunit type for corresponding pft
    real(r8), pointer :: wtgcell(:)       ! pft weight relative to grid cell
    real(r8), pointer :: tmomin20(:)      ! 20-yr running mean of tmomin
    real(r8), pointer :: agdd20(:)        ! 20-yr running mean of agdd
    real(r8), pointer :: agddtw(:)        ! accumulated growing degree days above twmax
    real(r8), pointer :: prec365(:)       ! 365-day running mean of tot. precipitation
    real(r8), pointer :: sla(:)           ! ecophys const - sp. leaf area [m2 leaf g-1 carbon]
    logical , pointer :: tree(:)          ! ecophys const - true=> tree is present
    real(r8), pointer :: crownarea_max(:) ! ecophys const - tree maximum crown area [m2]
    real(r8), pointer :: lm_sapl(:)       ! ecophys const - leaf mass of sapling
    real(r8), pointer :: sm_sapl(:)       ! ecophys const - stem mass of sapling
    real(r8), pointer :: hm_sapl(:)       ! ecophys const - heartwood mass of sapling
    real(r8), pointer :: rm_sapl(:)       ! ecophys const - root mass of saping
    real(r8), pointer :: reinickerp(:)    ! ecophys const - parameter in allometric equation
    real(r8), pointer :: wooddens(:)      ! ecophys const - wood density (gC/m3)
    real(r8), pointer :: latosa(:)        ! ecophys const - ratio of leaf area to sapwood cross-sectional area
    real(r8), pointer :: allom1(:)        ! ecophys const - parameter in allometric
    real(r8), pointer :: allom2(:)        ! ecophys const - parameter in allometric
    real(r8), pointer :: allom3(:)        ! ecophys const - parameter in allometric
    real(r8), pointer :: tcmin(:)         ! ecophys const - minimum coldest monthly mean temperature
    real(r8), pointer :: tcmax(:)         ! ecophys const - maximum coldest monthly mean temperature
    real(r8), pointer :: gddmin(:)        ! ecophys const - minimum growing degree days (at or above 5 C)
!
! local pointers to implicit in/out arguments
!
    integer , pointer :: ivt(:)           ! vegetation type for this pft
    logical , pointer :: present(:)       ! true=> PFT present in patch
    real(r8), pointer :: nind(:)          ! number of individuals (#/m**2)
    real(r8), pointer :: lm_ind(:)        ! individual leaf mass
    real(r8), pointer :: sm_ind(:)        ! individual sapwood mass
    real(r8), pointer :: hm_ind(:)        ! individual heartwood mass
    real(r8), pointer :: rm_ind(:)        ! individual root mass
    real(r8), pointer :: litterag(:)      ! above ground litter
    real(r8), pointer :: litterbg(:)      ! below ground litter
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: fpcgrid(:)       ! foliar projective cover on gridcell (fraction)
    real(r8), pointer :: htop(:)          ! canopy top (m)
    real(r8), pointer :: lai_ind(:)       ! LAI per individual
    real(r8), pointer :: crownarea(:)     ! area that each individual tree takes up (m^2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: g,l,p,m                        ! indices
    integer  :: fn, filterg(ubg-lbg+1)         ! local gridcell filter for error check
!
! gridcell level variables
!
    logical  :: grid_present(lbg:ubg,0:numpft) !true=>pft is present in gridcell
    logical  :: grid_survive(lbg:ubg,0:numpft) !true=>pft survives in gridcell
    logical  :: grid_estab(lbg:ubg,0:numpft)   !true=>pft is established in grirdcell
    integer  :: ngrass(lbg:ubg)                !counter
    integer  :: npft_estab(lbg:ubg)            !counter
    real(r8) :: fpc_tree_total(lbg:ubg)        !total fractional cover of trees in vegetated portion of gridcell
    real(r8) :: fpc_grass_total(lbg:ubg)       !total fractional cover of grass in vegetated portion of gridcell
    real(r8) :: fpc_total(lbg:ubg)             !old-total fractional vegetated portion of gridcell (without bare ground)
    real(r8) :: fpc_total_new(lbg:ubg)         !new-total fractional vegetated portion of gridcell (without bare ground)
    real(r8) :: grid_tmomin20(lbg:ubg)         !20 year running mean of minimum monthly temperature
    real(r8) :: grid_agdd20(lbg:ubg)           !20 year running mean of growing degree days
    real(r8) :: grid_agddtw(lbg:ubg)           !growing degree base tw
!
! pft level variables
!
    real(r8) :: fpc_ind(lbp:ubp)               !individual foliage projective cover
    real(r8) :: estab_rate(lbp:ubp)            !establishment rate
    real(r8) :: estab_grid(lbp:ubp)            !establishment rate on grid cell
!
! temporaries
!
    real(r8) :: bare_max                       !maximum bare soil
    real(r8) :: bare                           !fractional cover of bare soil
    real(r8) :: nind_old                       !old number of individuals
    real(r8) :: fpcgridtemp                    !temporary
    real(r8) :: sm_ind_temp                    !temporary
    real(r8) :: stemdiam                       !stem diameter
!
! minimum individual density for persistence of PFT (indiv/m2)
!
    real(r8), parameter :: nind_min = 1.0e-10_r8
!
! minimum precip. for establishment (mm/s)
!
    real(r8), parameter :: prec_min_estab = 100._r8/(365._r8*SHR_CONST_CDAY)
!
! maximum sapling establishment rate (indiv/m2)
!
    real(r8), parameter :: estab_max = 0.24_r8
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (landunit-level)

    ltype => clm3%g%l%itype

    ! Assign local pointers to derived type members (pft-level)

    ivt           => clm3%g%l%c%p%itype
    pgridcell     => clm3%g%l%c%p%gridcell
    plandunit     => clm3%g%l%c%p%landunit
    wtgcell       => clm3%g%l%c%p%wtgcell
    htop          => clm3%g%l%c%p%pps%htop
    present       => clm3%g%l%c%p%pdgvs%present
    prec365       => clm3%g%l%c%p%pdgvs%prec365
    present       => clm3%g%l%c%p%pdgvs%present
    nind          => clm3%g%l%c%p%pdgvs%nind
    fpcgrid       => clm3%g%l%c%p%pdgvs%fpcgrid
    litterag      => clm3%g%l%c%p%pdgvs%litterag
    litterbg      => clm3%g%l%c%p%pdgvs%litterbg
    lm_ind        => clm3%g%l%c%p%pdgvs%lm_ind
    sm_ind        => clm3%g%l%c%p%pdgvs%sm_ind
    hm_ind        => clm3%g%l%c%p%pdgvs%hm_ind
    rm_ind        => clm3%g%l%c%p%pdgvs%rm_ind
    crownarea     => clm3%g%l%c%p%pdgvs%crownarea
    lai_ind       => clm3%g%l%c%p%pdgvs%lai_ind
    tmomin20      => clm3%g%l%c%p%pdgvs%tmomin20
    agdd20        => clm3%g%l%c%p%pdgvs%agdd20
    agddtw        => clm3%g%l%c%p%pdgvs%agddtw
    crownarea_max => dgv_pftcon%crownarea_max
    lm_sapl       => dgv_pftcon%lm_sapl
    sm_sapl       => dgv_pftcon%sm_sapl
    hm_sapl       => dgv_pftcon%hm_sapl
    rm_sapl       => dgv_pftcon%rm_sapl
    reinickerp    => dgv_pftcon%reinickerp
    wooddens      => dgv_pftcon%wooddens
    latosa        => dgv_pftcon%latosa
    allom1        => dgv_pftcon%allom1
    allom2        => dgv_pftcon%allom2
    allom3        => dgv_pftcon%allom3
    tree          => dgv_pftcon%tree
    tcmax         => dgv_pftcon%tcmax
    tcmin         => dgv_pftcon%tcmin
    gddmin        => dgv_pftcon%gddmin
    sla           => pftcon%sla

    ! **********************************************************************
    ! Slevis version of LPJ's subr. bioclim
    ! Limits based on 20-year running averages of coldest-month mean
    ! temperature and growing degree days (5 degree base).
    ! For SURVIVAL, coldest month temperature and GDD should be
    ! at least as high as PFT-specific limits.
    ! For REGENERATION, PFT must be able to survive AND coldest month
    ! temperature should be no higher than a PFT-specific limit.
    ! **********************************************************************

    ! Initialize gridcell-level metrics

    do m = 0,numpft
       do g = lbg, ubg
          grid_present(g,m) = .false.
          grid_survive(g,m) = .false.
          grid_estab(g,m)   = .false.
       end do
    end do
    do g = lbg, ubg
       ngrass(g) = 0
       npft_estab(g) = 0
       fpc_tree_total(g) = 0._r8
       fpc_grass_total(g) = 0._r8
       fpc_total(g) = 0._r8
       fpc_total_new(g) = 0._r8
       grid_tmomin20(g) = 0._r8
       grid_agdd20(g) = 0._r8
       grid_agddtw(g) = 0._r8
    end do

    ! Calculate total woody FPC, FPC increment and grass cover (= crown area)

    do p = lbp, ubp
       g = pgridcell(p)

       ! Calculate the grid-average bioclimate variables for
       ! survival and establishment

       grid_tmomin20(g) = grid_tmomin20(g) + tmomin20(p) * wtgcell(p)
       grid_agdd20(g)   = grid_agdd20(g)   + agdd20(p)   * wtgcell(p)
       grid_agddtw(g)   = grid_agddtw(g)   + agddtw(p)   * wtgcell(p)

       ! Set the presence of pft for this gridcell
       ! Note: this modifies the pft-level vegetation type if present is false

       grid_present(g,ivt(p)) = present(p)
       if (.not. present(p)) ivt(p) = 0
    end do

    ! Must go thru all 16 pfts and decide which can/cannot establish or survive
    ! Determine present, survive, estab.  Note - if tmomin20 > tcmax then  crops,
    ! shrubs and 2nd boreal summergreen tree cannot exist yet (see EcosystemDynini)
    ! because they often coexist using up all pft patches and allowing for no bare
    ! ground. Thus they make fpc_grid_total < 1. Solve by allowing one more pft
    ! per grid cell than the number of pfts.

    do m = 1, numpft
!dir$ prefervector
       do g = lbg, ubg
          if (grid_tmomin20(g) >= tcmin(m) + SHR_CONST_TKFRZ ) then
             if (grid_tmomin20(g) <= tcmax(m) + SHR_CONST_TKFRZ  .and. &
                  grid_agdd20(g) >= gddmin(m) .and. nint(grid_agddtw(g)) == 0) then
                grid_estab(g,m) = .true.
             end if
             grid_survive(g,m) = .true.
          end if
       end do
    end do

    do p = lbp, ubp
       g = pgridcell(p)
       l = plandunit(p)

       ! Case 1 -- pft ceases to exist -kill pfts not adapted to current climate

       if (present(p) .and. (.not. grid_survive(g,ivt(p)) .or. nind(p)<nind_min)) then
          present(p) = .false.
          grid_present(g,ivt(p)) = .false.

          ! Add killed biomass to litter
          if (tree(ivt(p))) then
             litterag(p) = litterag(p) + nind(p) * (lm_ind(p) + sm_ind(p) + hm_ind(p))
          else              !if grass
             litterag(p) = litterag(p) + nind(p) * lm_ind(p)
          end if
          litterbg(p) = litterbg(p) + nind(p) * rm_ind(p)
          fpcgrid(p) = 0.0_r8
          ivt(p) = 0
       end if

       ! Case 2 -- pft begins to exist - introduce newly "adapted" pfts

       if (ltype(l) == istsoil) then
          if (.not. present(p) .and. prec365(p) >= prec_min_estab) then
             do m = 1, numpft
                if (ivt(p) /= m .and. .not. present(p)) then
                   if (.not. grid_present(g,m) .and. grid_estab(g,m)) then
                      present(p) = .true.
                      grid_present(g,m) = .true.
                      ivt(p) = m
                      if (tree(ivt(p))) then
                         nind(p) = 0.0_r8
                      else
                         nind(p) = 1.0_r8 !each grass PFT = 1 "individual"
                      end if
                      lm_ind(p) = 0.0_r8
                      sm_ind(p) = 0.0_r8
                      rm_ind(p) = 0.0_r8
                      hm_ind(p) = 0.0_r8
                      fpcgrid(p) = 0.0_r8
                      if (.not. tree(ivt(p))) crownarea(p) = 1.0_r8
                   end if   ! conditions suitable for establishment
                end if   ! no pft present and pft 'm' was not here until now
             end do   ! numpft
          end if   ! more conditions for establishment
       end if   ! if soil

       ! Case 3 -- some pfts continue to exist (no change) and some pfts
       ! continue to not exist (no change). Do nothing for this case.

    end do

    ! Sapling and grass establishment
    ! Calculate total woody FPC and number of woody PFTs present and able to establish

    do p = lbp, ubp
       g = pgridcell(p)
       if (present(p)) then
          if (tree(ivt(p))) then
             fpc_tree_total(g) = fpc_tree_total(g) + fpcgrid(p)
             if (grid_estab(g,ivt(p))) npft_estab(g) = npft_estab(g) + 1
          else if (.not.tree(ivt(p)) .and. ivt(p) > 0) then !grass
             ngrass(g) = ngrass(g) + 1
             fpc_grass_total(g) = fpc_grass_total(g) + fpcgrid(p)
          end if
          fpc_total(g) = fpc_total(g) + fpcgrid(p)
       end if
    end do

    ! Above establishment counters at the grid level required for the next steps.
    ! Note that ngrass, npft_estab, fpc_tree_total, fpc_grass_total, and fpc_total
    ! complete for gridcells.

    do p = lbp, ubp
       g = pgridcell(p)

       ! Prohibit establishment under extreme temperature or water stress.

       if (prec365(p) >= prec_min_estab .and. npft_estab(g) > 0) then

          ! Calculate establishment rate over available space, per tree PFT
          ! Maximum establishment rate reduced by shading as tree FPC approaches 1
          ! Total establishment rate partitioned equally among regenerating woody PFTs

          estab_rate(p) = estab_max * (1.0_r8-exp(5.0_r8*(fpc_tree_total(g)-1.0_r8))) / real(npft_estab(g))

          ! Calculate grid-level establishment rate per woody PFT
          ! Space available for woody PFT establishment is proportion of grid cell
          ! not currently occupied by woody PFTs

          estab_grid(p) = estab_rate(p) * (1.0_r8-fpc_tree_total(g))

       else ! if unsuitable climate for establishment

          estab_grid(p) = 0.0_r8

       end if

       if (present(p) .and. tree(ivt(p)) .and. grid_estab(g,ivt(p))) then

          ! Add new saplings to current population

          nind_old = nind(p)
          nind(p) = nind_old + estab_grid(p)
          lm_ind(p)   = (lm_ind(p) * nind_old + lm_sapl(ivt(p)) * estab_grid(p)) / nind(p)
          sm_ind_temp = (sm_ind(p) * nind_old + sm_sapl(ivt(p)) * estab_grid(p)) / nind(p)
          hm_ind(p)   = (hm_ind(p) * nind_old + hm_sapl(ivt(p)) * estab_grid(p)) / nind(p)
          rm_ind(p)   = (rm_ind(p) * nind_old + rm_sapl(ivt(p)) * estab_grid(p)) / nind(p)

          ! Calculate height, diameter and crown area for new average
          ! individual such that the basic allometric relationships (A-C below)
          ! are satisfied.
          ! (A) (leaf area) = latosa * (sapwood xs area)
          !        (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
          ! (B) (leaf mass) = lmtorm * (root mass)
          ! (C) height = allom2 * (stem diameter)**allom3  (source?)
          ! (D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
          !                                 crownarea_max)
          ! From (A),
          !  (1) sap_xsa = lm_ind * sla / latosa
          !  (2) wooddens = (sm_ind + hm_ind) / stemvolume
          !  (3) stemvolume = stem_xsa * height
          ! From (1), (2) & (3),
          !  (4) stem_xsa = (sm_ind + hm_ind) / wooddens / height
          !  (5) stem_xsa = SHR_CONST_PI * (stemdiam**2) / 4
          ! From (5),
          !  (6) stemdiam = ( 4 * stem_xsa / SHR_CONST_PI )**0.5
          ! From (4) & (6),
          !  (7) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens / height /
          !        SHR_CONST_PI )**0.5
          ! From (C) & (7),
          !  (8) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens /
          !        ( allom2 * stemdiam**allom3 ) / SHR_CONST_PI )**0.5
          ! From (8),
          !  (9) stemdiam = ( 4 * (sm_ind + hm_ind ) / wooddens / SHR_CONST_PI /
          !        allom2 )**( 1 / (2 + allom3) )

          stemdiam = (4.0_r8*(sm_ind_temp + hm_ind(p))/wooddens(ivt(p))/SHR_CONST_PI/ &
             allom2(ivt(p)))**(1.0_r8/(2.0_r8+allom3(ivt(p))))  ! Eqn 9
          htop(p) = allom2(ivt(p)) * stemdiam**allom3(ivt(p))                                    ! Eqn C
          crownarea(p) = min(crownarea_max(ivt(p)), allom1(ivt(p))*stemdiam**reinickerp(ivt(p))) ! Eqn D

          ! Recalculate sapwood mass, transferring excess sapwood to heartwood
          ! compartment, if necessary to satisfy Eqn A

          sm_ind(p) = lm_ind(p) * htop(p) * wooddens(ivt(p)) * sla(ivt(p)) / latosa(ivt(p))
          hm_ind(p) = hm_ind(p) + (sm_ind_temp-sm_ind(p))

          ! Update LAI and FPC

          if (crownarea(p) > 0.0_r8) then
             lai_ind(p) = lm_ind(p) * sla(ivt(p)) / crownarea(p)
          else
             lai_ind(p) = 0.0_r8
          end if

          fpc_ind(p) = 1.0_r8 - exp(-0.5_r8*lai_ind(p))
          fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind(p)
          fpc_total_new(g) = fpc_total_new(g) + fpcgrid(p)

       end if   ! add new saplings block
    end do   ! close loop to update fpc_total_new

    ! Adjustments- don't allow trees to exceed 95% of vegetated landunit

    do p = lbp, ubp
       g = pgridcell(p)
       if (fpc_total_new(g) > 0.95_r8) then
          if (tree(ivt(p)) .and. present(p)) then
             nind_old = nind(p)
             nind(p) = nind(p) / (fpc_total_new(g)/0.95_r8)
             fpcgrid(p) = fpcgrid(p) / (fpc_total_new(g)/0.95_r8)
             litterag(p) = litterag(p) + (nind_old - nind(p)) * (lm_ind(p) + sm_ind(p) + hm_ind(p))
             litterbg(p) = litterbg(p) + (nind_old - nind(p)) * rm_ind(p)
          end if
          fpc_total(g) = 0.95_r8
       end if
    end do

    ! Section for grasses. Grasses can establish in non-vegetated areas

!dir$ concurrent
!cdir nodep
    do p = lbp, ubp
       g = pgridcell(p)
       if (present(p) .and. .not. tree(ivt(p))) then
          if (grid_estab(g,ivt(p))) then
             if (ngrass(g) > 0) then
                bare = (1.0_r8 - fpc_total(g)) / real(ngrass(g))
             else
                bare = 0.0_r8
             end if
             bare_max = (-2.0_r8 * crownarea(p) * log(max(1.0_r8 - bare - fpcgrid(p), 0.000001_r8)) &
                  / sla(ivt(p)) - lm_ind(p)) / lm_sapl(ivt(p))
             bare = max(0.0_r8, min(bare, bare_max))
             lm_ind(p) = lm_ind(p) + bare * lm_sapl(ivt(p))
             rm_ind(p) = rm_ind(p) + bare * rm_sapl(ivt(p))
          end if
          if (lm_ind(p) <= 0.0_r8) then
             present(p) = .false.
             litterbg(p) = litterbg(p) + rm_ind(p) * nind(p)
          end if
       end if

    end do   ! end of pft-loop

    ! Recalculate fraction of vegetated landunit for each pft and do error check

    do g = lbg, ubg
       fpc_total(g) = 0.0_r8
    end do

    do p = lbp, ubp
       g = pgridcell(p)
       if (present(p)) then
          if (crownarea(p) > 0.0_r8) then
             lai_ind(p) = lm_ind(p) * sla(ivt(p)) / crownarea(p)
          else
             lai_ind(p) = 0.0_r8
          end if
          fpc_ind(p) = 1.0_r8 - exp(-0.5_r8*lai_ind(p))
          fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind(p)
          fpc_total(g) = fpc_total(g) + fpcgrid(p)
       else
          fpcgrid(p) = 0.0_r8
       end if
    end do

    ! Check for error in establishment
    fn = 0
    do g = lbg, ubg
       if (fpc_total(g) < 0.0_r8 .or. fpc_total(g) > 1.15_r8) then
          fn = fn + 1
          filterg(fn) = g
       end if
    end do
    ! Just print out the first error
    if (fn > 0) then
       g = filterg(1)
       write(6,*) 'Error in Establishment: fpc_total is ',fpc_total(g),' at gridcell ',g
       call endrun
    end if

    ! Adjustment b/c fpc_total > 1. This can happen, because grasses are allowed
    ! to grow on the bare soil as determined before updating the trees
    ! This may imply that LPJ allows a brief coexistence of some
    ! trees and grasses. LSM cannot allow such coexistence.
    ! ivt >= 12 includes all grasses.

    do p = lbp, ubp
       g = pgridcell(p)
       if (fpc_total(g) > 1.0_r8) then
          if (ivt(p) >= 12 .and. fpcgrid(p) > 0.0_r8) then
             fpcgridtemp = fpcgrid(p)
             fpcgrid(p) = max(0.0_r8, fpcgrid(p) - (fpc_total(g)-1.0_r8))
             if (fpcgrid(p) == 0.0_r8) then
                present(p) = .false.
                litterag(p) = litterag(p) + nind(p) * lm_ind(p)
                litterbg(p) = litterbg(p) + nind(p) * rm_ind(p)
             else  ! recalculate lai_ind for consistency with fpcgrid
                lai_ind(p) = -2.0_r8 * log(1.0_r8 - fpcgrid(p)/(crownarea(p)*nind(p)))
             end if
             fpc_total(g) = fpc_total(g) - fpcgridtemp + fpcgrid(p)
          end if
       end if
    end do

    ! More adjustments. Should C pools be reset, to balance carbon?

    do g = lbg, ubg
       fpc_total_new(g) = 0.0_r8
    end do

    do p = lbp, ubp
       g = pgridcell(p)

       ! Avoid fpcgrid=0 with present=.T.

       if (fpcgrid(p) == 0._r8) present(p) = .false.

       ! Avoid ivt(p)=0 with lai>0.

       if (.not. present(p)) lai_ind(p) = 0.0_r8

       ! Set the fpcgrid for bare ground if there is bare ground in
       ! vegetated landunit and pft is bare ground so that everything
       ! can add up to one.

       if (fpc_total(g) < 1.0_r8 .and. ivt(p) == noveg) then
          fpcgrid(p) = 1.0_r8 - fpc_total(g)
          fpc_total(g) = 1.0_r8
       end if

       ! Update fpc_total_new

       fpc_total_new(g) = fpc_total_new(g) + fpcgrid(p)
    end do

    ! Check for error in establishment
    fn = 0
    do g = lbg, ubg
       if (abs(fpc_total_new(g) - 1.0_r8) > 1.0e-6_r8) then
          fn = fn + 1
          filterg(fn) = g
       end if
    end do
    ! Just print out the first error
    if (fn > 0) then
       g = filterg(1)
       write(6,*) 'Error in Establishment: fpc_total_new =',fpc_total_new(g),&
            ' at gridcell ',g
       call endrun
    end if

  end subroutine Establishment

#endif

end module DGVMEstablishmentMod
