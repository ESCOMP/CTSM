#include <misc.h>
#include <preproc.h>

module CNVegStructUpdateMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNVegStructUpdateMod
!
! !DESCRIPTION:
! Module for vegetation structure updates (LAI, SAI, htop, hbot)
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varcon  , only: istsoil
    use spmdMod     , only: masterproc
    use clm_varpar  , only: nlevsoi
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: CNVegStructUpdate
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNVegStructUpdate
!
! !INTERFACE:
subroutine CNVegStructUpdate(num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, use C state variables and epc to diagnose
! vegetation structure (LAI, SAI, height)
!
! !USES:
   use clmtype
   use clm_atmlnd   , only: clm_a2l
   use pftvarcon    , only: noveg
   use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilp                 ! number of column soil points in pft filter
   integer, intent(in) :: filter_soilp(:)   ! pft filter for soil points
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 10/28/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)         ! pft vegetation type
   integer , pointer :: pcolumn(:)  ! column index associated with each pft
   integer , pointer :: pgridcell(:)   ! pft's gridcell index
   real(r8), pointer :: snowdp(:)   ! snow height (m)
   real(r8), pointer :: leafc(:)              ! (kgC/m2) leaf C
   real(r8), pointer :: deadstemc(:)          ! (kgC/m2) dead stem C
   real(r8), pointer :: woody(:)                        !binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: slatop(:)    !specific leaf area at top of canopy, projected area basis [m^2/gC]
   real(r8), pointer :: dsladlai(:)  !dSLA/dLAI, projected area basis [m^2/gC]
   real(r8), pointer :: z0mr(:)      !ratio of momentum roughness length to canopy top height (-)
   real(r8), pointer :: displar(:)   !ratio of displacement height to canopy top height (-)
   real(r8), pointer :: forc_hgt_u(:)  ! observational height of wind [m]
!
! local pointers to implicit in/out scalars
!
   integer , pointer :: frac_veg_nosno_alb(:) ! frac of vegetation not covered by snow [-]
   real(r8), pointer :: tlai(:) !one-sided leaf area index, no burying by snow
   real(r8), pointer :: tsai(:) !one-sided stem area index, no burying by snow
   real(r8), pointer :: htop(:) !canopy top (m)
   real(r8), pointer :: hbot(:) !canopy bottom (m)
   real(r8), pointer :: elai(:)     ! one-sided leaf area index with burying by snow
   real(r8), pointer :: esai(:)     ! one-sided stem area index with burying by snow
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p,c,g        !indices
   integer :: fp           !lake filter indices
   real(r8):: taper        ! ratio of height:radius_breast_height (tree allometry)
   real(r8):: stocking     ! #stems / ha (stocking density)
   real(r8):: dwood        ! density of wood (kgC/m^3)
   real(r8):: ol           ! thickness of canopy layer covered by snow (m)
   real(r8):: fb           ! fraction of canopy layer covered by snow

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    pcolumn                        => clm3%g%l%c%p%column
    pgridcell                      => clm3%g%l%c%p%gridcell
    leafc                          => clm3%g%l%c%p%pcs%leafc
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    snowdp                         => clm3%g%l%c%cps%snowdp
    forc_hgt_u                     => clm_a2l%forc_hgt_u
    woody                          => pftcon%woody
    slatop                         => pftcon%slatop
    dsladlai                       => pftcon%dsladlai
    z0mr                           => pftcon%z0mr
    displar                        => pftcon%displar

   ! assign local pointers to derived type arrays (out)
    tlai                           => clm3%g%l%c%p%pps%tlai
    tsai                           => clm3%g%l%c%p%pps%tsai
    htop                           => clm3%g%l%c%p%pps%htop
    hbot                           => clm3%g%l%c%p%pps%hbot
    elai                           => clm3%g%l%c%p%pps%elai
    esai                           => clm3%g%l%c%p%pps%esai
    frac_veg_nosno_alb             => clm3%g%l%c%p%pps%frac_veg_nosno_alb

   ! constant allometric parameters
   taper = 200._r8
   stocking = 1000._r8

   ! convert from stems/ha -> stems/m^2
   stocking = stocking / 10000._r8

   ! a typical value for wood density is 500 kg dry mass/m^3
   dwood = 500._r8

   ! convert from kg -> g
   dwood = dwood * 1000._r8

   ! convert from dry mass -> Carbon
   dwood = dwood * 0.5_r8

   ! pft loop
!dir$ concurrent
!cdir nodep
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)
      g = pgridcell(p)

      if (ivt(p) /= noveg) then
          ! update the leaf area index based on leafC and SLA
          !tlai(p) = (2._r8 * slatop(ivt(p)) * leafc(p))/(2._r8 - dsladlai(ivt(p)) * leafc(p))
          if (dsladlai(ivt(p)) > 0._r8) then
             tlai(p) = (exp(leafc(p)*dsladlai(ivt(p)) + log(slatop(ivt(p)))) - slatop(ivt(p)))/dsladlai(ivt(p))
          else
             tlai(p) = slatop(ivt(p)) * leafc(p)
          end if
          ! added 5/4/04, PET: duringg exp38 debugging
          tlai(p) = max(0._r8, tlai(p))

          ! update the stem area index and height based on LAI, stem mass, and veg type.
          ! With the exception of htop for woody vegetation, this follows the DGVM logic.

          if (woody(ivt(p)) == 1._r8) then
             ! trees and shrubs
             tsai(p) = 0.25_r8 * tlai(p)

             ! trees and shrubs for now have a very simple allometry, with hard-wired
             ! stem taper (height:radius) and hard-wired stocking density (#individuals/area)
             htop(p) = ((3._r8 * deadstemc(p) * taper * taper)/(SHR_CONST_PI * stocking * dwood))**(1._r8/3._r8)

             ! Peter Thornton, 5/3/2004
             ! Adding test to keep htop from getting too close to forcing height for windspeed
             ! Also added for grass, below, although it is not likely to ever be an issue.
             htop(p) = min(htop(p),(forc_hgt_u(g)/(displar(ivt(p))+z0mr(ivt(p))))-3._r8)

             ! Peter Thornton, 8/11/2004
             ! Adding constraint to keep htop from going to 0.0.
             ! This becomes an issue when fire mortality is pushing deadstemc
             ! to 0.0.
             htop(p) = max(htop(p), 0.01_r8)

             hbot(p) = max(0._r8, min(3._r8, htop(p)-1._r8))

          else
             ! grasses
             tsai(p) = 0.05_r8 * tlai(p)

             ! height for grasses depends only on LAI
             htop(p) = max(0.25_r8, tlai(p) * 0.25_r8)

             htop(p) = min(htop(p),(forc_hgt_u(g)/(displar(ivt(p))+z0mr(ivt(p))))-3._r8)

             ! Peter Thornton, 8/11/2004
             ! Adding constraint to keep htop from going to 0.0.
             htop(p) = max(htop(p), 0.01_r8)

             hbot(p) = max(0.0_r8, min(0.05_r8, htop(p)-0.20_r8))
          end if

      else
          tlai(p) = 0._r8
          tsai(p) = 0._r8
          htop(p) = 0._r8
          hbot(p) = 0._r8
      end if
      
      ! adjust lai and sai for burying by snow. 

      ol = min(max(snowdp(c)-hbot(p), 0._r8), htop(p)-hbot(p))
      fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
      elai(p) = max(tlai(p)*fb, 0.0_r8)
      esai(p) = max(tsai(p)*fb, 0.0_r8)

      ! Fraction of vegetation free of snow
      if ((elai(p) + esai(p)) > 0._r8) then
         frac_veg_nosno_alb(p) = 1
      else
         frac_veg_nosno_alb(p) = 0
      end if

   end do

end subroutine CNVegStructUpdate
!-----------------------------------------------------------------------

end module CNVegStructUpdateMod
