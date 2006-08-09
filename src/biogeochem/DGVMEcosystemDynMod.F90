#include <misc.h>
#include <preproc.h>

module DGVMEcosystemDynMod

#if (defined DGVM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: DGVMEcosystemDynMod
!
! !DESCRIPTION:
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public  :: DGVMEcosystemDynini ! LPJ and other DGVM related initializations
  public  :: DGVMEcosystemDyn    ! Ecosystem dynamics: phenology, vegetation
  public  :: DGVMRespiration     ! Compute plant respiration
!
! !PUBLIC MEMBER FUNCTIONS:
  private :: Phenology     ! Compute summer and drought phenology
  private :: FireSeason    ! Calculate length of fire season in a year
  private :: LitterSOM     ! Calculate litter and soil decomposition
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DGVMEcosystemDynini
!
! !INTERFACE:
  subroutine DGVMEcosystemDynini()
!
! !DESCRIPTION:
! LPJ and other DGVM related initializations
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use nanMod
    use clmtype
    use decompMod    , only : get_proc_bounds, get_proc_global
    use clm_varpar   , only : numpft, npftpar, maxpatch_pft
    use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_TKFRZ
    use pftvarcon    , only : pftpar, tree, summergreen, raingreen, sla, &
                              lm_sapl, sm_sapl, hm_sapl, rm_sapl, latosa, &
                              allom1, allom2, allom3, reinickerp, wooddens, &
                              noveg
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from LPJ initialization subroutines)
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: p,n                  ! indices
    integer  :: begp, endp           ! per-proc beginning and ending pft indices
    integer  :: begc, endc           ! per-proc beginning and ending column indices
    integer  :: begl, endl           ! per-proc beginning and ending landunit indices
    integer  :: begg, endg           ! per-proc gridcell ending gridcell indices
    integer  :: pft
    real(r8) :: x
    real(r8) :: stemdiam
    real(r8) :: height_sapl
    real(r8) :: table(0:numpft,1:npftpar)
    type(pft_type), pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    pptr => clm3%g%l%c%p

    ! PFT PARAMETERS (follows LPJ subroutine pftparameters)

    !  1  fraction of roots in upper soil layer
    !  2  plants with C4 (1) or C3 (0) photosynthetic pathway
    !  3  water scalar value at which leaves shed by drought deciduous PFT
    !  4  canopy conductance component (gmin, mm/s) not associated with
    !     photosynthesis (Haxeltine & Prentice 1996, Table 4)
    !  5  maintenance respiration coefficient
    !  6  flammability threshold
    !  7  maximum foliar N content (mg/g)
    !     (Haxeltine & Prentice 1996a, Fig 4)
    !  8  fire resistance index
    !  9  leaf turnover period (years)
    ! 10  leaf longevity (years)
    ! 11  sapwood turnover period (sapwood converted to heartwood) (years)
    ! 12  root turnover period (years)
    ! 13  leaf C:N mass ratio
    ! 14  sapwood C:N mass ratio
    ! 15  root C:N mass ratio
    ! 16  leaf type: broadleaved (1), needleleaved (2) or grass (3)
    ! 17  phenology type: evergreen (1), summergreen (2), raingreen (3),
    !     any type (4)
    ! 18  leaf to root ratio under non-water stressed conditions
    ! 19  summergreen phenology ramp, GDD5 requirement to grow full leaf canopy
    ! 20  tree maximum crown area (m2)
    ! 21  sapling (or grass on initialisation) LAI
    ! 22  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
    ! 23  boreal pft (1), non-boreal pft (0)
    ! 24  low temperature limit for CO2 uptake
    ! 25  lower range of temperature optimum for photosynthesis
    ! 26  upper range of temperature optimum for photosynthesis
    ! 27  high temperature limit for CO2 unptake

    !BIOCLIMATIC LIMITS

    ! 28 minimum coldest monthly mean temperature
    ! 29 maximum coldest monthly mean temperature
    ! 30 minimum growing degree days (at or above 5 deg C)
    ! 31 upper limit of temperature of the warmest month (twmax)
    ! 32 lower limit of growth efficiency (g/m2)


    ! ---------------------------------------------------------------------
    !      1      2      3      4      5      6      7      8          PFT
    ! ---------------------------------------------------------------------

    data ((table(pft,n),n=1,8),pft=0,numpft) /               &
          inf,   inf,   inf,   inf,   inf,  0.15_r8,   inf,   inf, &      !  0
         0.70_r8,   0.0_r8,  0.00_r8,   0.3_r8,  1.20_r8,  0.15_r8, 100.0_r8,  0.12_r8, &      !  1
         0.90_r8,   0.0_r8,  0.00_r8,   0.3_r8,  0.60_r8,  0.15_r8, 100.0_r8,  0.12_r8, &      !  2 was 1.20
         0.90_r8,   0.0_r8,  0.00_r8,   0.3_r8,  0.60_r8,  0.15_r8, 100.0_r8,  0.12_r8, &      !  3 was 1.20
         0.85_r8,   0.0_r8,  0.00_r8,   0.5_r8,  0.50_r8,  0.15_r8, 100.0_r8,  0.12_r8, &      !  4 was 0.20
         0.70_r8,   0.0_r8,  0.00_r8,   0.5_r8,  1.20_r8,  0.15_r8, 100.0_r8,  0.50_r8, &      !  5
         0.70_r8,   0.0_r8,  0.35_r8,   0.5_r8,  0.50_r8,  0.15_r8, 100.0_r8,  0.50_r8, &      !  6 was 0.20
         0.80_r8,   0.0_r8,  0.00_r8,   0.5_r8,  1.20_r8,  0.15_r8, 120.0_r8,  0.12_r8, &      !  7
         0.90_r8,   0.0_r8,  0.00_r8,   0.3_r8,  0.60_r8,  0.15_r8, 100.0_r8,  0.12_r8, &      !  8 was 1.20
         0.85_r8,   0.0_r8,  0.00_r8,   0.5_r8,  1.20_r8,  0.15_r8, 100.0_r8,  0.12_r8, &      !  9 (=4or5)
         0.80_r8,   0.0_r8,  0.00_r8,   0.5_r8,  1.20_r8,  0.15_r8, 120.0_r8,  0.12_r8, &      ! 10
         0.90_r8,   0.0_r8,  0.00_r8,   0.3_r8,  0.60_r8,  0.15_r8, 100.0_r8,  0.12_r8, &      ! 11 was 1.20
         0.90_r8,   0.0_r8,  0.35_r8,   0.5_r8,  0.60_r8,  0.15_r8, 100.0_r8,  1.00_r8, &      ! 12 was 1.20
         0.90_r8,   0.0_r8,  0.35_r8,   0.5_r8,  0.60_r8,  0.15_r8, 100.0_r8,  1.00_r8, &      ! 13 was 1.20
         0.90_r8,   1.0_r8,  0.35_r8,   0.5_r8,  1.20_r8,  0.15_r8, 100.0_r8,  1.00_r8, &      ! 14
         0.90_r8,   1.0_r8,  0.35_r8,   0.5_r8,  1.20_r8,  0.15_r8, 100.0_r8,  1.00_r8, &      ! 15 (.not.
         0.90_r8,   1.0_r8,  0.35_r8,   0.5_r8,  1.20_r8,  0.15_r8, 100.0_r8,  1.00_r8/        ! 16 present)


    ! ---------------------------------------------------------------------
    !      9     10     11     12     13     14     15     16     17   PFT
    ! ---------------------------------------------------------------------

    data ((table(pft,n),n=9,17),pft=0,numpft) /                     &
         inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf,   inf, & ! 0
         2.0_r8,  2.00_r8,  20.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   2.0_r8,   1.0_r8, & ! 1
         2.0_r8,  2.00_r8,  20.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   2.0_r8,   1.0_r8, & ! 2
         1.0_r8,  0.50_r8,  20.0_r8,   1.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   2.0_r8,   2.0_r8, & ! 3
         2.0_r8,  2.00_r8,  20.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   1.0_r8,   1.0_r8, & ! 4
         1.0_r8,  1.00_r8,  20.0_r8,   1.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   1.0_r8,   1.0_r8, & ! 5
         1.0_r8,  0.50_r8,  20.0_r8,   1.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   1.0_r8,   3.0_r8, & ! 6
         1.0_r8,  0.50_r8,  20.0_r8,   1.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   1.0_r8,   2.0_r8, & ! 7
         1.0_r8,  0.50_r8,  20.0_r8,   1.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   2.0_r8,   2.0_r8, & ! 8
         2.0_r8,  2.00_r8,  20.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   1.0_r8,   1.0_r8, & ! 9
         1.0_r8,  0.50_r8,  20.0_r8,   1.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   1.0_r8,   2.0_r8, & !10
         1.0_r8,  0.50_r8,  20.0_r8,   1.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   2.0_r8,   2.0_r8, & !11
         1.0_r8,  1.00_r8,   1.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   3.0_r8,   4.0_r8, & !12
         1.0_r8,  1.00_r8,   1.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   3.0_r8,   4.0_r8, & !13
         1.0_r8,  1.00_r8,   1.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   3.0_r8,   4.0_r8, & !14
         1.0_r8,  1.00_r8,   1.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   3.0_r8,   4.0_r8, & !15
         1.0_r8,  1.00_r8,   1.0_r8,   2.0_r8,  29.0_r8, 330.0_r8,  29.0_r8,   3.0_r8,   4.0_r8/   !16


    ! ------------------------------------------------------
    !       18      19     20      21    22     23     PFT
    ! ------------------------------------------------------

    data ((table(pft,n),n=18,23),pft=0,numpft) /  &
         inf,    inf,   inf,    inf,  inf,   inf, &  ! 0
         1.0_r8, 1000.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   0.0_r8, &  ! 1
         1.0_r8, 1000.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   1.0_r8, &  ! 2
         1.0_r8,  200.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   1.0_r8, &  ! 3
         1.0_r8, 1000.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   0.0_r8, &  ! 4
         1.0_r8, 1000.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   0.0_r8, &  ! 5
         1.0_r8, 1000.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   0.0_r8, &  ! 6
         1.0_r8,  200.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   0.0_r8, &  ! 7
         1.0_r8,  200.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   1.0_r8, &  ! 8
         1.0_r8, 1000.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   0.0_r8, &  ! 9
         1.0_r8,  200.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   0.0_r8, &  !10
         1.0_r8,  200.0_r8,  15.0_r8,  1.500_r8,  1.2_r8,   1.0_r8, &  !11
        0.75_r8,  100.0_r8,   0.0_r8,  0.001_r8,  1.2_r8,   1.0_r8, &  !12
        0.75_r8,  100.0_r8,   0.0_r8,  0.001_r8,  1.2_r8,   1.0_r8, &  !13
        0.75_r8,  100.0_r8,   0.0_r8,  0.001_r8,  1.2_r8,   0.0_r8, &  !14
        0.75_r8,  100.0_r8,   0.0_r8,  0.001_r8,  1.2_r8,   0.0_r8, &  !15
        0.75_r8,  100.0_r8,   0.0_r8,  0.001_r8,  1.2_r8,   0.0_r8/    !16

    ! -------------------------------------
    !      24     25     26      27    PFT
    ! -------------------------------------
    data ((table(pft,n),n=24,27),pft=0,numpft) / &
          inf,   inf,   inf,    inf, & ! 0
         -4.0_r8,  20.0_r8,  30.0_r8,   42.0_r8, & ! 1
         -4.0_r8,  15.0_r8,  25.0_r8,   38.0_r8, & ! 2
         -4.0_r8,  15.0_r8,  25.0_r8,   38.0_r8, & ! 3
          2.0_r8,  25.0_r8,  30.0_r8,   55.0_r8, & ! 4
         -4.0_r8,  20.0_r8,  30.0_r8,   42.0_r8, & ! 5
          2.0_r8,  25.0_r8,  30.0_r8,   55.0_r8, & ! 6
         -4.0_r8,  20.0_r8,  25.0_r8,   38.0_r8, & ! 7
         -4.0_r8,  15.0_r8,  25.0_r8,   38.0_r8, & ! 8
          2.0_r8,  25.0_r8,  30.0_r8,   55.0_r8, & ! 9
         -4.0_r8,  20.0_r8,  25.0_r8,   38.0_r8, & !10
         -4.0_r8,  15.0_r8,  25.0_r8,   38.0_r8, & !11
         -4.0_r8,  10.0_r8,  30.0_r8,   45.0_r8, & !12
         -4.0_r8,  10.0_r8,  30.0_r8,   45.0_r8, & !13
          6.0_r8,  20.0_r8,  45.0_r8,   55.0_r8, & !14
          6.0_r8,  20.0_r8,  45.0_r8,   55.0_r8, & !15
          6.0_r8,  20.0_r8,  45.0_r8,   55.0_r8/   !16

    ! --------------------------------------------------------
    !      28       29      30       31      32     PFT
    ! --------------------------------------------------------
    data ((table(pft,n),n=28,npftpar),pft=0,numpft) / &
         inf,     inf,    inf,  1000.0_r8,    inf, & !  0
         -2.0_r8,    22.0_r8,  900.0_r8,  1000.0_r8,    0.0_r8, & !  1
        -32.5_r8,    -2.0_r8,  600.0_r8,    23.0_r8,    0.0_r8, & !  2
       9999.9_r8,    -2.0_r8,  350.0_r8,    23.0_r8,    0.0_r8, & !  3 (was -1000.0)
         15.5_r8,  1000.0_r8,    0.0_r8,  1000.0_r8,    0.0_r8, & !  4
          3.0_r8,    18.8_r8, 1200.0_r8,  1000.0_r8,    0.0_r8, & !  5
         15.5_r8,  1000.0_r8,    0.0_r8,  1000.0_r8,    0.0_r8, & !  6
        -17.0_r8,    15.5_r8, 1200.0_r8,  1000.0_r8,    0.0_r8, & !  7
      -1000.0_r8,    -2.0_r8,  350.0_r8,    23.0_r8,    0.0_r8, & !  8
       9999.9_r8,  1000.0_r8,    0.0_r8,  1000.0_r8,    0.0_r8, & !  9 (was    15.5)
       9999.9_r8,    15.5_r8, 1200.0_r8,  1000.0_r8,    0.0_r8, & ! 10 (was   -17.0)
       9999.9_r8,    -2.0_r8,  350.0_r8,    23.0_r8,    0.0_r8, & ! 11 (was -1000.0)
      -1000.0_r8,   -17.0_r8,    0.0_r8,  1000.0_r8,    0.0_r8, & ! 12 (an LSM type)
        -17.0_r8,    15.5_r8,    0.0_r8,  1000.0_r8,    0.0_r8, & ! 13 (was -1000.0)
         15.5_r8,  1000.0_r8,    0.0_r8,  1000.0_r8,    0.0_r8, & ! 14
       9999.9_r8,  1000.0_r8,    0.0_r8,  1000.0_r8,    0.0_r8, & ! 15 (was    15.5)
       9999.9_r8,  1000.0_r8,    0.0_r8,  1000.0_r8,    0.0_r8/   ! 16 (was    15.5)
    !----------------------------------------------------------------------------

    pftpar(noveg,:)    = table(noveg,:)
    sla(noveg)         = inf
    lm_sapl(noveg)     = inf
    rm_sapl(noveg)     = inf
    sm_sapl(noveg)     = inf
    hm_sapl(noveg)     = inf
    tree(noveg)        = .false.
    summergreen(noveg) = .false.
    raingreen(noveg)   = .false.

    do pft = 1,numpft

       ! Transfer parameter values to array pftpar

       do n = 1, npftpar
          pftpar(pft,n) = table(pft,n)
       end do

       ! Assign leaf and phenology logicals

       if (pftpar(pft,16) <= 2.0_r8) then     !woody vegetation: trees, shrubs
          tree(pft) = .true.
       else                                !non woody vegetation: grasses
          tree(pft) = .false.
       end if

       if     (pftpar(pft,17) == 1.0_r8) then !evergreen
          summergreen(pft) = .false.
          raingreen(pft)   = .false.
       elseif (pftpar(pft,17) == 2.0_r8) then !summergreen
          summergreen(pft) = .true.
          raingreen(pft)   = .false.
       elseif (pftpar(pft,17) == 3.0_r8) then !raingreen
          summergreen(pft) = .false.
          raingreen(pft)   = .true.
       else                                !any of the above
          summergreen(pft) = .true.
          raingreen(pft)   = .true.
       end if

       ! Calculate specific leaf area (SLA) for each PFT from leaf longevity
       ! Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
       ! Equation based on Reich et al 1997, Fig 1f:

       ! SLA = 2e-4 * exp(6.15 - 0.46 ln (leaf_longevity * 12))

       ! SLA in m2/gC, leaf_longevity in years

       sla(pft) = 2.0e-4_r8 * exp(6.15_r8 - 0.46_r8*log(pftpar(pft,10)*12.0_r8))

       ! Define initial mass structure

       if (tree(pft)) then !woody PFTs

          ! Calculate leafmass for a sapling individual
          !  (1) lai = leafmass * sla / (crown area)
          !  (2) (leaf area) = latosa * (sapwood xs area)
          !         (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
          !  (3) (crown area) = allom1 * (stem diameter) ** reinickerp
          !         (Reinickes theory)
          ! From (1),
          !  (4) leafmass = lai * (crown area) / sla
          ! From (1) & (3),
          !  (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
          ! From (2),
          !  (6) leafmass = latosa * (sapwood xs area) / sla
          !  (7) (sapwood xs area) = SHR_CONST_PI * (sapwood diameter)**2 / 4
          ! From (6) and (7),
          !  (8) leafmass = latosa * SHR_CONST_PI * (sapwood diameter)**2 / 4 / sla
          ! From (8),
          !  (9) (sapwood diameter) = [ 4 * leafmass * sla / SHR_CONST_PI / latosa ]**0.5
          ! (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
          ! Define x,
          ! (11) x = [ (sapwood diameter)+(heartwood diameter) ] /
          !          (sapwood diameter)
          ! From (10) & (11),
          ! (12) (stem diameter) = x * (sapwood diameter)
          ! From (5), (9) & (12),
          ! (13) leafmass = lai * allom1 * x**reinickerp *
          !               (4*leafmass*sla/SHR_CONST_PI/latosa)**(reinickerp*0.5) / sla
          ! From (13),
          ! (14) leafmass = [ lai * allom1 * x**reinickerp *
          !      (4*sla/SHR_CONST_PI/latosa)**(reinickerp*0.5) / sla ]**(2/(2-reinickerp))

          x = pftpar(pft,22)

          lm_sapl(pft) = (pftpar(pft,21) * allom1 * x**reinickerp *          &
               (4.0_r8 * sla(pft) / SHR_CONST_PI / latosa)**(reinickerp * 0.5_r8) / sla(pft))** &
               (2.0_r8/(2.0_r8-reinickerp)) !eqn 14

          ! Calculate sapling stem diameter
          ! From (9) & (12),
          ! (15) (stem diameter) = x * [ 4 * leafmass * sla / SHR_CONST_PI / latosa ]**0.5

          stemdiam = x * (4.0_r8*lm_sapl(pft)*sla(pft)/SHR_CONST_PI/latosa)**0.5_r8 !Eqn 15

          ! Calculate sapling height
          ! (16) height = allom2 * (stem diameter)**allom3 (source?)

          height_sapl = allom2 * stemdiam**allom3 !Eqn 16

          ! Calculate sapling sapwood mass
          ! (17) (sapwood volume) = height * (sapwood xs area)
          ! (18) (sapwood xs area) = leafmass * sla / latosa
          ! From (17) & (18),
          ! (19) (sapwood volume) = height * leafmass * sla / latosa
          ! (20) (sapwood mass) = (wood density) * (sapwood volume)
          ! From (19) & (20),
          ! (21) (sapwood mass) = (wood density) * height * leafmass * sla / latosa

          sm_sapl(pft)=wooddens*height_sapl*lm_sapl(pft)*sla(pft)/latosa !Eqn 21

          ! Calculate sapling heartwood mass
          ! From (11),
          ! (22) (heartwood mass) = (x-1) * (sapwood mass)

          hm_sapl(pft) = (x-1.0_r8) * sm_sapl(pft) !Eqn 22

       else !grass PFTs

          lm_sapl(pft) = pftpar(pft,21) / sla(pft)

       end if

       ! Calculate sapling or initial grass rootmass
       ! (23) lmtorm = (leafmass) / (rootmass)
       ! where lmtorm=pftpar(pft,18)

       rm_sapl(pft) = lm_sapl(pft) / pftpar(pft,18) !From Eqn 23

    end do ! pft type loop

    ! ---------------------------------------------------------------
    ! Some of the following comes from LPJ subroutine initgrid
    ! ---------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

!dir$ concurrent
!cdir nodep
    do p = begp,endp

       ! used in Phenology
       pptr%pdgvs%t10min(p) = 1.0e+36_r8
       pptr%pdgvs%lai_ind(p) = 0.0_r8

       ! updated in Phenology
       pptr%pdgvs%dphen(p) = 0.0_r8
       pptr%pdgvs%leafon(p) = 0.0_r8
       pptr%pdgvs%leafof(p) = 0.0_r8

       ! accumulated in FireSeason (must reset at end of every year)
       pptr%pdgvs%firelength(p) = 0.0_r8

       ! used in FireSeason; updated in LitterSOM; updated in annual portion of LPJ
       pptr%pdgvs%litterag(p) = 0.0_r8
       pptr%pdgvs%litterbg(p) = 0.0_r8

       ! updated in LitterSOM
       pptr%pdgvs%cpool_fast(p) = 0.0_r8
       pptr%pdgvs%cpool_slow(p) = 0.0_r8
       pptr%pdgvs%k_fast_ave(p) = 0.0_r8
       pptr%pdgvs%k_slow_ave(p) = 0.0_r8
       pptr%pdgvs%litter_decom_ave(p) = 0.0_r8
       pptr%pcf%fmicr(p) = 0.0_r8 !initialize b/c use in Biogeochemistry before LitterSOM

       ! used and updated in annual portion of LPJ
       pptr%pdgvs%present(p)   = .false.
       pptr%pdgvs%nind(p)      = 0._r8
       pptr%pdgvs%lm_ind(p)    = 0._r8
       pptr%pdgvs%sm_ind(p)    = 0._r8
       pptr%pdgvs%hm_ind(p)    = 0._r8
       pptr%pdgvs%rm_ind(p)    = 0._r8
       pptr%pdgvs%tmomin20(p)  = SHR_CONST_TKFRZ - 5._r8 !initialize this way for Phenology code
       pptr%pdgvs%agdd20(p)    = 0._r8
       pptr%pdgvs%t_mo_min(p)  = 1.0e+36_r8
       pptr%pdgvs%crownarea(p) = 0._r8

       ! already a variable in LSM but now updated in annual portion of LPJ
       ! note - fpcgrid is not relevant where ist/=istsoil (consistent with subr. surfrd)
       pptr%pps%htop(p) = 0._r8
       pptr%pps%tsai(p) = 0._r8
       pptr%pdgvs%fpcgrid(p) = 1.0_r8 / maxpatch_pft !
       !
       ! accumulated in Biogeochemistry and used/reset in annual portion of LPJ
       pptr%pdgvs%bm_inc(p) = 0.0_r8
       pptr%pdgvs%afmicr(p) = 0.0_r8

       ! accumulated in Stomata and used/reset in annual portion of LPJ
       pptr%pdgvs%annpsn(p) = 0.0_r8
       pptr%pdgvs%annpsnpot(p) = 0.0_r8

    end do

  end subroutine DGVMEcosystemDynini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DGVMEcosystemDyn
!
! !INTERFACE:
  subroutine DGVMEcosystemDyn(lbp, ubp, num_nolakep, filter_nolakep, &
                              doalb, endofyr)
!
! !DESCRIPTION:
! Ecosystem dynamics: phenology, vegetation
! Calculates leaf areas (tlai, elai),  stem areas (tsai, esai) and
! height (htop)
!
! !USES:
    use clmtype
    use shr_const_mod, only: SHR_CONST_CDAY
    use clm_time_manager, only : get_step_size, get_nstep, get_curr_date
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
    logical, intent(in) :: doalb                       ! true = surface albedo calculation time step
    logical, intent(in) :: endofyr
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 2/1/02, Peter Thornton: Migrated to new data structure.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: ivt(:)      ! vegetation type for this pft
    integer , pointer :: pcolumn(:)  ! column index of corresponding pft
    real(r8), pointer :: lai_ind(:)  ! LAI per individual
    real(r8), pointer :: dphen(:)    ! phenology [0 to 1]
    real(r8), pointer :: fpcgrid(:)  ! foliar projective cover on gridcell (fraction)
    real(r8), pointer :: snowdp(:)   ! snow height (m) (column-level)
!
! local pointers to implicit in/out arguments
!
    real(r8), pointer :: htop(:)     ! canopy top (m)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: tlai(:)     ! one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai(:)     ! one-sided stem area index, no burying by snow
    real(r8), pointer :: hbot(:)     ! canopy bottom (m)
    real(r8), pointer :: elai(:)     ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)     ! one-sided stem area index with burying by snow
    integer , pointer :: frac_veg_nosno_alb(:) ! frac of vegetation not covered by snow [-]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: fp,p,c  ! indices
    real(r8) :: ol      ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb      ! fraction of canopy layer covered by snow
    integer  :: nstep   ! model time step number
    real(r8) :: dtime   ! land model time step (sec)
    integer  :: year    ! current year (0 -> ...)
    integer  :: month   ! current month (1 -> 12)
    integer  :: day     ! current day (1 -> 31)
    integer  :: secs    ! seconds of current date
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type scalar members (column-level)

    snowdp             => clm3%g%l%c%cps%snowdp

    ! Assign local pointers to derived type scalar members (pft-level)

    ivt                => clm3%g%l%c%p%itype
    pcolumn            => clm3%g%l%c%p%column
    tlai               => clm3%g%l%c%p%pps%tlai
    tsai               => clm3%g%l%c%p%pps%tsai
    elai               => clm3%g%l%c%p%pps%elai
    esai               => clm3%g%l%c%p%pps%esai
    htop               => clm3%g%l%c%p%pps%htop
    hbot               => clm3%g%l%c%p%pps%hbot
    htop               => clm3%g%l%c%p%pps%htop
    frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
    lai_ind            => clm3%g%l%c%p%pdgvs%lai_ind
    dphen              => clm3%g%l%c%p%pdgvs%dphen
    fpcgrid            => clm3%g%l%c%p%pdgvs%fpcgrid

    ! today's phenology: returns dphen
    ! with small chg in subr, could be called every tstep if preferred

    nstep = get_nstep()
    dtime = get_step_size()
    if (nstep==0 .or. mod(nstep-1, nint(SHR_CONST_CDAY/dtime))==0._r8 .or. endofyr) then
       call Phenology(lbp, ubp, num_nolakep, filter_nolakep)
    end if

    ! fire season; returns firelength for use at end of yr in subr. Fire
    ! litter and soil decomposition; returns litterag,bg and fmicr

    if (nstep > 1) then
       call FireSeason(lbp, ubp, num_nolakep, filter_nolakep)
       call get_curr_date(year, month, day, secs)
       call LitterSOM(lbp, ubp, num_nolakep, filter_nolakep, year)
    end if

    if (doalb .or. endofyr) then

!dir$ concurrent
!cdir nodep
       do fp = 1, num_nolakep
          p = filter_nolakep(fp)
          c = pcolumn(p)

          ! obtain vegetation structure here (tlai+sai,hbot+top)

          if (ivt(p) > 0 .and. fpcgrid(p) > 0.0_r8) then
             if (ivt(p) <= 11) then          !woody vegetation (following IBIS)
                tlai(p) = lai_ind(p) * dphen(p)
                tsai(p) = 0.250_r8 * lai_ind(p)
                hbot(p) = max(0.0_r8, min(3.00_r8, htop(p)-1.00_r8))
             else                             !grasses (following IBIS)
                tlai(p) = lai_ind(p) * dphen(p) / fpcgrid(p)
                tsai(p) = 0.050_r8 * tlai(p)
                htop(p) = max(0.25_r8, tlai(p) * 0.25_r8)
                hbot(p) = max(0.0_r8, min(0.05_r8, htop(p)-0.20_r8))
             end if
          else
             tlai(p) = 0.0_r8
             tsai(p) = 0.0_r8
             htop(p) = 0.0_r8
             hbot(p) = 0.0_r8
          end if

          ! adjust lai and sai for burying by snow. if exposed lai and sai
          ! are less than 0.05 but non zero, set to 0.05 to prevent numerical
          ! problems associated with very small lai and sai.

          ol = min( max(snowdp(c)-hbot(p),0._r8), htop(p)-hbot(p))
          fb = 1._r8 - ol / max(1.e-06_r8,htop(p)-hbot(p))
          elai(p) = max(tlai(p)*fb, 0.0_r8)
          esai(p) = max(tsai(p)*fb, 0.0_r8)
          if (elai(p) > 0.0_r8 .and. elai(p) < 0.05_r8) elai(p) = 0.05_r8
          if (esai(p) > 0.0_r8 .and. esai(p) < 0.05_r8) esai(p) = 0.05_r8

          ! Fraction of vegetation free of snow
          if ((elai(p) + esai(p)) >= 0.05_r8) then
             frac_veg_nosno_alb(p) = 1
          else
             frac_veg_nosno_alb(p) = 0
          end if

       end do ! end of pft loop

    end if !end of if-doalb block

  end subroutine DGVMEcosystemDyn

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DGVMRespiration
!
! !INTERFACE:
  subroutine DGVMRespiration(lbc, ubc, lbp, ubp, &
                             num_nolakec, filter_nolakec, &
                             num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Calculates surface biogeochemical fluxes
!
! !USES:
    use clmtype
    use shr_const_mod , only: SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_varpar    , only : nlevsoi
    use clm_time_manager  , only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(num_nolakec) ! column filter for non-lake points
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(num_nolakep) ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Gordon Bonan and Sam Levis
! 1/31/02, PET: migrated to new data structures
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)    ! pft's column
    integer , pointer :: ivt(:)        ! vegetation type for current pft
    real(r8), pointer :: t_veg(:)      ! vegetation temperature (Kelvin)
    real(r8), pointer :: fpcgrid(:)    ! foliar projective cover on gridcell (fraction)
    real(r8), pointer :: nind(:)       ! number of individuals
    real(r8), pointer :: dphen(:)      ! phenology [0 to 1]
    real(r8), pointer :: lm_ind(:)     ! individual leaf mass
    real(r8), pointer :: sm_ind(:)     ! individual stem mass
    real(r8), pointer :: rm_ind(:)     ! individual root mass
    real(r8), pointer :: respcoeff(:)  ! respiration coefficient (LPJ)
    real(r8), pointer :: l_cton(:)     ! c/n for leaves (LPJ)
    real(r8), pointer :: s_cton(:)     ! c/n for stems (LPJ)
    real(r8), pointer :: r_cton(:)     ! c/n for roots (LPJ)
    real(r8), pointer :: z(:,:)        ! layer thickness (m)
    real(r8), pointer :: dz(:,:)       ! layer depth (m)
    real(r8), pointer :: t_soisno(:,:) ! soil temperature (Kelvin)
    real(r8), pointer :: fpsn(:)       ! photosynthesis (umol CO2 /m**2 /s)
    real(r8), pointer :: fmicr(:)      ! microbial respiration (umol CO2 /m**2 /s)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: bm_inc(:)     ! biomass increment
    real(r8), pointer :: afmicr(:)     ! annual microbial respiration
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: frmf(:)       ! leaf maintenance respiration  (umol CO2 /m**2 /s)
    real(r8), pointer :: frms(:)       ! stem maintenance respiration  (umol CO2 /m**2 /s)
    real(r8), pointer :: frmr(:)       ! root maintenance respiration  (umol CO2 /m**2 /s)
    real(r8), pointer :: frm(:)        ! total maintenance respiration (umol CO2 /m**2/s)
    real(r8), pointer :: frg(:)        ! growth respiration (umol CO2 /m**2 /s)
    real(r8), pointer :: dmi(:)        ! total dry matter production (ug /m**2 /s)
    real(r8), pointer :: fco2(:)       ! net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
    real(r8), pointer :: tsoi25(:)     ! soil temperature to 0.25 m (Kelvin)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    real(r8), parameter :: k = 0.0548_r8 / SHR_CONST_CDAY ! from [/day] to [/second]
    integer  :: fp,p,fc,c,j   ! indices
    real(r8) :: tf            ! temperature factor
    real(r8) :: dtime         ! land model time step (sec)
    real(r8) :: dmcf          ! co2-to-biomass conversion (ug biomass / umol CO2)
    real(r8) :: tsoi(lbc:ubc) ! temporary
    real(r8) :: dep(lbc:ubc)  ! temporary
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    pcolumn  => clm3%g%l%c%p%column
    z        => clm3%g%l%c%cps%z
    dz       => clm3%g%l%c%cps%dz
    t_soisno => clm3%g%l%c%ces%t_soisno

    ! Assign local pointers to derived type members (pft-level)

    ivt       => clm3%g%l%c%p%itype
    pcolumn   => clm3%g%l%c%p%column
    t_veg     => clm3%g%l%c%p%pes%t_veg
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    nind      => clm3%g%l%c%p%pdgvs%nind
    dphen     => clm3%g%l%c%p%pdgvs%dphen
    lm_ind    => clm3%g%l%c%p%pdgvs%lm_ind
    sm_ind    => clm3%g%l%c%p%pdgvs%sm_ind
    rm_ind    => clm3%g%l%c%p%pdgvs%rm_ind
    bm_inc    => clm3%g%l%c%p%pdgvs%bm_inc
    afmicr    => clm3%g%l%c%p%pdgvs%afmicr
    tsoi25    => clm3%g%l%c%p%pdgvs%tsoi25
    fmicr     => clm3%g%l%c%p%pcf%fmicr
    fpsn      => clm3%g%l%c%p%pcf%fpsn
    frmf      => clm3%g%l%c%p%pcf%frmf
    frms      => clm3%g%l%c%p%pcf%frms
    frmr      => clm3%g%l%c%p%pcf%frmr
    frm       => clm3%g%l%c%p%pcf%frm
    frg       => clm3%g%l%c%p%pcf%frg
    dmi       => clm3%g%l%c%p%pcf%dmi
    fco2      => clm3%g%l%c%p%pcf%fco2

    ! Assign local pointers to ecophysiological constants (pft-level)

    l_cton    => dgv_pftcon%l_cton
    s_cton    => dgv_pftcon%s_cton
    r_cton    => dgv_pftcon%r_cton
    respcoeff => dgv_pftcon%respcoeff

    ! Determine time step and step size

    dtime = get_step_size()

    ! Soil temperature to a depth of 0.25 m.

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       tsoi(c) = 0._r8
       dep(c) = 0._r8
    end do
    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (z(c,j)+0.5_r8*dz(c,j) <= 0.25_r8) then
             tsoi(c) = tsoi(c) + t_soisno(c,j)*dz(c,j)
             dep(c) = dep(c) + dz(c,j)
          end if
       end do
    end do

    ! Loop over pfts

!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)

       if (dep(c) /= 0._r8) then
          tsoi25(p) = tsoi(c)/dep(c)
       else
          tsoi25(p) = t_soisno(c,1)
       end if

       ! Set co2-to-biomass conversion (ug biomass / umol CO2)

       dmcf = 28.5_r8

       ! maintenance respiration: LPJ equations w/ units of [gC m-2 gridcell s-1]
       ! converted to LSM units of [umol CO2 m-2 patch s-1]

       if (ivt(p) > 0 .and. fpcgrid(p) > 0.0_r8) then
          if (t_veg(p) >= SHR_CONST_TKFRZ-40._r8) then
             tf = exp(308.56_r8 * (1.0_r8/56.02_r8 - 1.0_r8/(t_veg(p)-227.13_r8)))
          else
             tf = 0.0_r8
          end if

          frmf(p) = respcoeff(ivt(p)) * k * lm_ind(p) * nind(p) / l_cton(ivt(p)) &
               * tf * dphen(p) * 2.0e6_r8 / dmcf / fpcgrid(p)
          frms(p) = respcoeff(ivt(p)) * k * sm_ind(p) * nind(p) / s_cton(ivt(p)) &
               * tf * 2.0e6_r8 / dmcf / fpcgrid(p)

          if (tsoi25(p) >= SHR_CONST_TKFRZ-40._r8) then
             tf = exp(308.56_r8 * (1.0_r8/56.02_r8 - 1.0_r8/(tsoi25(p)-227.13_r8)))
          else
             tf = 0.0_r8
          end if

          frmr(p) = respcoeff(ivt(p)) * k * rm_ind(p) * nind(p) / r_cton(ivt(p)) &
               * tf * dphen(p) * 2.0e6_r8 / dmcf / fpcgrid(p)
       else
          frmf(p) = 0.0_r8
          frms(p) = 0.0_r8
          frmr(p) = 0.0_r8
       end if

       frm(p)  = frmf(p) + frms(p) + frmr(p)

       ! growth respiration and production

       frg(p) = 0.25_r8 * max(fpsn(p) - frm(p), 0.0_r8)      !changed to match LPJ
       dmi(p) = (fpsn(p) - frm(p) - frg(p)) * dmcf

       ! bm_inc=[gC/m2 patch area] from dmi=[ug dry matter/m2 patch area/s]

       bm_inc(p) = bm_inc(p) + dmi(p) * dtime * 0.5e-6_r8

       ! microbial respiration

       ! DGVM calculates clm%fmicr in LitterSOM; problem with units in relation
       ! to history grid averages; {fmicr}=[gC/m2 gridcell vegetated area] in
       ! LPJ calculation => sum(fmicr) over a gridcell would give the total for
       ! the gridcell; it would be wrong to convert to [/m2 patch area] because
       ! soil carbon continues to exist even when the plants above it die and
       ! fpcgrid goes to 0; how to reconcile with history calculation which
       ! will take the following values and weight them by the corresponding
       ! weights of their patches? Could chg. soil carbon and plant litter to
       ! gridcell level pools; this will affect fco2, as well; could treat this
       ! issue simultaneously with soil water, which we also want converted to
       ! the gridcell level for plant competition. (slevis)

       afmicr(p) = afmicr(p) + fmicr(p) ![gC/m2 gridcell vegetated area]

       ! net CO2 flux

       fco2(p) = -fpsn(p) + frm(p) + frg(p) + fmicr(p)

    end do ! end of pft loop

  end subroutine DGVMRespiration

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Phenology
!
! !INTERFACE:
  subroutine Phenology (lbp, ubp, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Summer and drought phenology. Called once per day.
!
! !USES:
    use clmtype
    use shr_const_mod, only : SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Ecosysdyn in this module
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Jon Foley's IBIS subroutine pheno)
! 2/1/02, Peter Thornton: Migrated to new data structure
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer , pointer :: ivt(:)         ! vegetation type for this pft
    real(r8), pointer :: t10(:)         ! 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: agdd0(:)       ! accumulated growing degree days above 0 deg C
    real(r8), pointer :: agdd5(:)       ! accumulated growing degree days above -5
    real(r8), pointer :: fnpsn10(:)     ! 10-day running mean net photosynthesis
    real(r8), pointer :: tmomin20(:)    ! 20 year running mean of monthly minimum
    real(r8), pointer :: l_long(:)      ! ecophys constant - leaf longevity [years]
    logical , pointer :: tree(:)        ! ecophys constant
    logical , pointer :: raingreen(:)   ! ecophys constant
    logical , pointer :: summergreen(:) ! ecophys constant
!
! local pointers to implicit in/out scalars
!
    real(r8), pointer :: t10min(:)      ! annual minimum of 10-day running mean (K)
    real(r8), pointer :: leafon(:)      ! leafon days
    real(r8), pointer :: leafof(:)      ! leafoff days
!
! local pointers to implicit out scalars
!
    real(r8), pointer :: dphen(:)       ! phenology [0 to 1]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    real(r8), parameter :: ddfacu = 1.0_r8 / 15.0_r8 ! 'drop day factor'
    real(r8), parameter :: ddfacl = 1.0_r8 /  5.0_r8 ! logy to switch in 15 or 5 days
    integer  :: fp,p                           ! indices
    real(r8) :: tthreshold
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (pft-level)

    ivt         => clm3%g%l%c%p%itype
    t10         => clm3%g%l%c%p%pdgvs%t10
    t10min      => clm3%g%l%c%p%pdgvs%t10min
    agdd0       => clm3%g%l%c%p%pdgvs%agdd0
    agdd5       => clm3%g%l%c%p%pdgvs%agdd5
    fnpsn10     => clm3%g%l%c%p%pdgvs%fnpsn10
    leafon      => clm3%g%l%c%p%pdgvs%leafon
    leafof      => clm3%g%l%c%p%pdgvs%leafof
    dphen       => clm3%g%l%c%p%pdgvs%dphen
    tmomin20    => clm3%g%l%c%p%pdgvs%tmomin20
    l_long      => dgv_pftcon%l_long
    tree        => dgv_pftcon%tree
    raingreen   => dgv_pftcon%raingreen
    summergreen => dgv_pftcon%summergreen

!dir$ concurrent
!cdir nodep
    do fp = 1, num_nolakep
       p = filter_nolakep(fp)

       t10min(p) = min (t10min(p), t10(p)) !reset to 1.0e+36 once per year

       if (tree(ivt(p)) .and. .not.summergreen(ivt(p)) .and. .not.raingreen(ivt(p))) then

          dphen(p) = 1.0_r8

       else if (tree(ivt(p)) .and. summergreen(ivt(p))) then

          ! ---------------------------------------------------------------------
          ! * * * upper canopy winter phenology * * *
          ! ---------------------------------------------------------------------
          ! temperature threshold for budburst and senescence
          ! temperature threshold is assumed to be SHR_CONST_TKFRZ
          ! or 5 degrees warmer than the coldest monthly temperature
          ! tmomin20 is initialized to tfrz-5.0

          tthreshold = max (SHR_CONST_TKFRZ, tmomin20(p) + 5.0_r8)

          ! determine if growing degree days are initiated
          ! slevis:*t10 = 10-day running mean air temperature (K)
          ! determine leaf display

          if (t10(p) < tthreshold) then
             dphen(p) = max (0.0_r8, dphen(p) - ddfacu)
          else
             dphen(p) = min (1.0_r8, max (0.0_r8, agdd0(p) - 100.0_r8) / 50.0_r8)
          end if

       else if (ivt(p) > 0 .and. .not.tree(ivt(p))) then !NB: grass has no specific phenology

          ! ---------------------------------------------------------------------
          ! * * * lower canopy phenology * * *
          ! ---------------------------------------------------------------------
          ! temperature threshold for budburst and senescence
          ! temperature threshold is assumed to be SHR_CONST_TKFRZ

          tthreshold = SHR_CONST_TKFRZ

          ! determine leaf display

          if (t10(p) < tthreshold) then                ! cold phenology for grasses
             dphen(p) = max (0.0_r8, dphen(p) - ddfacl)   ! slevis: made ddfacl=1/5
          else if (fnpsn10(p) < 0.0_r8) then                 ! drought phenology for grasses
             dphen(p) = max (0.1_r8, dphen(p) - ddfacl)   ! slevis: made ddfacl=1/5
          else
            !dphen(p) = min (1.0, max (0.0, agdd5(p) - 150.0) / 50.0)
             dphen(p) = min (1.0_r8, dphen(p) + ddfacl)   !slevis: try this line instead
          end if

       end if

       if (tree(ivt(p)) .and. raingreen(ivt(p))) then

          ! ---------------------------------------------------------------------
          ! * * * upper canopy drought phenology * * *
          ! ---------------------------------------------------------------------

          if (fnpsn10(p) <  0.0_r8) dphen(p) = max (0.1_r8, dphen(p) - ddfacu)
          if (fnpsn10(p) >= 0.0_r8) dphen(p) = min (1.0_r8, dphen(p) + ddfacu)

          ! trying out enforced drought phenology

          if (dphen(p) > 0.95_r8) leafon(p) = leafon(p) + 1.0_r8
          if (leafon(p) >= 365.0_r8*l_long(ivt(p))) then
             dphen(p) = 0.1_r8
             leafof(p) = leafof(p) + 1.0_r8
             if (leafof(p) >= 365.0_r8*l_long(ivt(p))) then
                dphen(p)  = 1.0_r8
                leafof(p) = 0.0_r8
                leafon(p) = 1.0_r8
             end if
          end if

       end if

    end do

  end subroutine Phenology

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: FireSeason
!
! !INTERFACE:
  subroutine FireSeason (lbp, ubp, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Calculate length of fire season in a year
! Orig. code was called once per day.
! slevis adapted to call every tstep.
! Orig. code operated on a grid cell basis.
! slevis adapted to operate on a patch basis.
!
! !USES:
    use clmtype
    use clm_time_manager, only : get_step_size
    use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_CDAY, SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Ecosysdyn in this module
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine fire)
!
! !LOCAL VARIABLES:
!
! local pointers toimplicit in arguments
!
    integer , pointer :: pcolumn(:)      ! column index for corresponding pft
    integer , pointer :: ivt(:)          ! vegetation type for this pft
    real(r8), pointer :: t_ref2m(:)      ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: litterag(:)     ! above ground litter
    real(r8), pointer :: wf(:)           ! soil water as frac. of whc for top 0.5 m
    real(r8), pointer :: flam(:)         ! ecophys constant - flammability threshold [units?]
!
! local pointers toimplicit in/out arguments
!
    real(r8), pointer :: firelength(:)   ! fire season in days
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer  :: fp, p, c   ! indices
    real(r8) :: dtime      ! land model time step (sec)
    real(r8) :: fire_prob  ! fire probability (tsteps)
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    wf         => clm3%g%l%c%cps%wf

    ! Assign local pointers to derived type members (pft-level)

    pcolumn    => clm3%g%l%c%p%column
    ivt        => clm3%g%l%c%p%itype
    t_ref2m    => clm3%g%l%c%p%pes%t_ref2m
    litterag   => clm3%g%l%c%p%pdgvs%litterag
    firelength => clm3%g%l%c%p%pdgvs%firelength
    flam       => dgv_pftcon%flam

    ! Get model step size

    dtime = get_step_size()

    ! Calculate the length of the fire season (in days)
    ! Calculate today's fire probability, fire_prob
    ! Assume fire is only possible when temperature is above SHR_CONST_TKFRZ
    ! slevis: *wf is top 0.5 m soil water as a fraction of the whc
    !         *divide fire_prob (days) by tsteps/day to get fire_prob (tsteps)
    !         *else need daily avg t_ref2m and wf to calculate fire_prob

!dir$ concurrent
!cdir nodep
    do fp = 1, num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)

       if (t_ref2m(p) > SHR_CONST_TKFRZ .and. litterag(p) > 0.0_r8) then
          fire_prob = EXP((-SHR_CONST_PI/4.0_r8) * (max(0.0_r8,wf(c))/flam(ivt(p)))**2) &
               * dtime / SHR_CONST_CDAY
       else
          fire_prob = 0.0_r8
       end if

       firelength(p) = firelength(p) + fire_prob  ! reset 1/yr in subroutine lpj
    end do

  end subroutine FireSeason

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: LitterSOM
!
! !INTERFACE:
  subroutine LitterSOM (lbp, ubp, num_nolakep, filter_nolakep, kyr)
!
! !DESCRIPTION:
! Litter and soil decomposition
! Incorporates analytical solution for soil pool sizes
! once litter inputs are (assumed to be) at equilibrium,
! reducing spin-up time for carbon fluxes due to soil respiration.
!
! !USES:
    use clmtype
    use clm_time_manager , only : get_step_size
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
    integer, intent(in) :: kyr                         ! year (0, ...) for nstep+1
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis (adapted from Stephen Sitch's LPJ subroutine littersom)
! 2/1/02, Peter Thornton: Migrate to new data structures
!
! !LOCAL VARIABLES:
!
! local pointers to  implicit in arguments
!
    integer , pointer :: pcolumn(:)          ! column index for corresponding pft
    real(r8), pointer :: wf(:)               ! soil water as frac. of whc for top 0.5 m
    real(r8), pointer :: tsoi25(:)           ! soil temperature to 0.25 m (Kelvin)
!
! local pointers to  implicit in/out arguments
!
    real(r8), pointer :: litterag(:)         ! above ground litter
    real(r8), pointer :: litterbg(:)         ! below ground litter
    real(r8), pointer :: cpool_fast(:)       ! fast carbon pool
    real(r8), pointer :: cpool_slow(:)       ! slow carbon pool
    real(r8), pointer :: k_fast_ave(:)       ! decomposition rate
    real(r8), pointer :: k_slow_ave(:)       ! decomposition rate
    real(r8), pointer :: litter_decom_ave(:) ! decomposition rate
!
! local pointers to  implicit out arguments
!
    real(r8), pointer :: fmicr(:)            ! microbial respiration (umol CO2 /m**2 /s)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer , parameter :: soil_equil_year = 400     ! number of years until pool sizes for soil decomposition solved analytically
    real(r8), parameter :: k_litter10 = 0.5_r8          ! litter decomp. rate at 10 deg C (/year)
    real(r8), parameter :: k_soil_fast10 = 0.03_r8      ! fast pool decomp. rate at 10 deg C (/year)
    real(r8), parameter :: k_soil_slow10 = 0.001_r8     ! slow pool decomp. rate at 10 deg C (/year)
    real(r8), parameter :: atmfrac = 0.7_r8             ! fraction of litter decomp. going directly into the atmosphere
    real(r8), parameter :: soilfrac = 1.0_r8 - atmfrac  ! fraction of litter decomp. going to soil C pools
    real(r8), parameter :: fastfrac = 0.985_r8          ! fraction of litter entering fast soil decomposition pool
    real(r8), parameter :: slowfrac = 1.0_r8 - fastfrac ! fraction of litter entering slow soil decomposition pool
    integer  :: fp,p,c             ! indices
    real(r8) :: dtime              ! land model time step (sec)
    real(r8) :: temp_resp          ! temperature response of decomposition
    real(r8) :: moist_resp         ! moisture response of decomposition
    real(r8) :: k_litter           ! litter decomposition rate (/tstep)
    real(r8) :: k_fast             ! fast pool decomposition rate (/tstep)
    real(r8) :: k_slow             ! slow pool decomposition rate (/tstep)
    real(r8) :: litter_decom       ! litter decomposition
    real(r8) :: litter_decom_ag    ! above-ground litter decomposition
    real(r8) :: litter_decom_bg    ! below-ground litter decomposition
    real(r8) :: cflux_litter_soil  ! litter decomposition flux to soil
    real(r8) :: cflux_litter_atmos ! litter decomposition flux to atmosphere
    real(r8) :: cflux_fast_atmos   ! soil fast pool decomposition flux to atmos.
    real(r8) :: cflux_slow_atmos   ! soil slow pool decomposition flux to atmos.
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    wf               => clm3%g%l%c%cps%wf

    ! Assign local pointers to derived type scalar members

    pcolumn          => clm3%g%l%c%p%column
    fmicr            => clm3%g%l%c%p%pcf%fmicr
    tsoi25           => clm3%g%l%c%p%pdgvs%tsoi25
    litterag         => clm3%g%l%c%p%pdgvs%litterag
    litterbg         => clm3%g%l%c%p%pdgvs%litterbg
    cpool_fast       => clm3%g%l%c%p%pdgvs%cpool_fast
    cpool_slow       => clm3%g%l%c%p%pdgvs%cpool_slow
    k_fast_ave       => clm3%g%l%c%p%pdgvs%k_fast_ave
    k_slow_ave       => clm3%g%l%c%p%pdgvs%k_slow_ave
    litter_decom_ave => clm3%g%l%c%p%pdgvs%litter_decom_ave

    ! Get model step size

    dtime = get_step_size()

    ! Determine litter and soil decomposition

!dir$ concurrent
!cdir nodep
    do fp = 1, num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)

       ! Temperature response function is a modified Q10 relationship
       ! (Lloyd & Taylor 1994)
       ! slevis: Original code used monthly avg soil temp (K); I use tstep value

       if (tsoi25(p) <= SHR_CONST_TKFRZ - 40.0_r8) then !avoid division by zero
          temp_resp=0.0_r8
       else                            !Lloyd & Taylor 1994
          temp_resp=exp(308.56_r8*((1.0_r8/56.02_r8)-(1.0_r8/(tsoi25(p)-227.13_r8))))
       end if

       ! Moisture response based on soil layer 1 moisture content (Foley 1995)
       ! slevis: Orig. code used monthly soil water in upper 0.5 m (fraction of whc)
       !         I use the tstep value

       moist_resp = 0.25_r8 + 0.75_r8 * wf(c)

       ! Original divided by 12 to get monthly decomposition rates (k, /month)
       ! as a function of temperature and moisture
       ! slevis: make rates /tstep by dividing by the number of tsteps per year

       k_litter = k_litter10    * temp_resp * moist_resp * dtime / (SHR_CONST_CDAY * 365._r8)
       k_fast   = k_soil_fast10 * temp_resp * moist_resp * dtime / (SHR_CONST_CDAY * 365._r8)
       k_slow   = k_soil_slow10 * temp_resp * moist_resp * dtime / (SHR_CONST_CDAY * 365._r8)

       ! Calculate monthly litter decomposition using equation
       !   (1) dc/dt = -kc     where c=pool size, t=time, k=decomposition rate
       ! from (1),
       !   (2) c = c0*exp(-kt) where c0=initial pool size
       ! from (2), decomposition in any month given by
       !   (3) delta_c = c0 - c0*exp(-k)
       ! from (3)
       !   (4) delta_c = c0*(1.0-exp(-k))

       litter_decom_ag = litterag(p) * (1.0_r8-exp(-k_litter))  !eqn 4
       litter_decom_bg = litterbg(p) * (1.0_r8-exp(-k_litter))
       litter_decom    = litter_decom_ag + litter_decom_bg

       ! Update the litter pools

       litterag(p) = litterag(p) - litter_decom_ag
       litterbg(p) = litterbg(p) - litter_decom_bg

       ! Calculate carbon flux to atmosphere and soil

       cflux_litter_atmos = atmfrac  * litter_decom
       cflux_litter_soil  = soilfrac * litter_decom

       ! Further subdivide soil fraction between fast and slow soil pools


       cpool_fast(p) = cpool_fast(p) + fastfrac * cflux_litter_soil
       cpool_slow(p) = cpool_slow(p) + slowfrac * cflux_litter_soil

       ! Calculate monthly soil decomposition to the atmosphere

       cflux_fast_atmos = cpool_fast(p) * (1.0_r8-exp(-k_fast))  !eqn 4
       cflux_slow_atmos = cpool_slow(p) * (1.0_r8-exp(-k_slow))  !eqn 4

       ! Update the soil pools

       cpool_fast(p) = cpool_fast(p) - cflux_fast_atmos
       cpool_slow(p) = cpool_slow(p) - cflux_slow_atmos

       ! Calculate heterotrophic respiration (in LSM referred to as microbial)

       fmicr(p) = cflux_litter_atmos + cflux_fast_atmos + cflux_slow_atmos

       ! Empty soil pools below a minimum threshold

       if (cpool_fast(p) < 1.0e-5_r8) cpool_fast(p) = 0.0_r8
       if (cpool_slow(p) < 1.0e-5_r8) cpool_slow(p) = 0.0_r8

       if (kyr <= soil_equil_year) then

          ! Update running average respiration rates and litter input
          ! slevis: had to multiply the denominator to chg units from years to tsteps
          k_fast_ave(p)       = k_fast_ave(p)       + k_fast / &
               (real(soil_equil_year) * 365._r8 * SHR_CONST_CDAY / dtime)
          k_slow_ave(p)       = k_slow_ave(p)       + k_slow / &
               (real(soil_equil_year) * 365._r8 * SHR_CONST_CDAY / dtime)
          litter_decom_ave(p) = litter_decom_ave(p) + litter_decom / &
               (real(soil_equil_year) * 365._r8 * SHR_CONST_CDAY / dtime)

       else if (kyr == soil_equil_year + 1) then

          ! SOIL DECOMPOSITION EQUILIBRIUM CALCULATION
          ! Analytical solution of differential flux equations for fast and slow
          ! soil carbon pools.  Implemented after (soil_equil_year) simulation
          ! years, when annual litter inputs should be close to equilibrium.  Assumes
          ! average climate (temperature and soil moisture) from all years up to
          ! soil_equil_year.
          ! slevis: next could be done once

          ! Analytically calculate pool sizes this year only
          ! Rate of change of soil pool size = litter input - decomposition
          !   (5) dc/dt = litter_decom - kc
          ! At equilibrium,
          !   (6) dc/dt = 0
          ! From (5) & (6),
          !   (7) c = litter_decom / k

          cpool_fast(p) = soilfrac * fastfrac * litter_decom_ave(p) / k_fast_ave(p) !eqn 7
          cpool_slow(p) = soilfrac * slowfrac * litter_decom_ave(p) / k_slow_ave(p) !eqn 7

       end if

    end do

  end subroutine LitterSOM

#endif

end module DGVMEcosystemDynMod
