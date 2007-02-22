#include <misc.h>
#include <preproc.h>

module VOCEmissionMod

#if (defined VOC)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: VOCEmissionMod
!
! !DESCRIPTION:
! Volatile organic compound emission
!
! !USES:
  use abortutils, only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: VOCEmission
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
! !IROUTINE: VOCEmission
!
! !INTERFACE:
  subroutine VOCEmission (lbp, ubp, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Volatile organic compound emission
! This code simulates volatile organic compound emissions
! following the algorithm presented in Guenther, A., 1999: Modeling
! Biogenic Volatile Organic Compound Emissions to the Atmosphere. In
! Reactive Hydrocarbons in the Atmosphere, Ch. 3
! This model relies on the assumption that 90% of isoprene and monoterpene
! emissions originate from canopy foliage:
!    E = epsilon * gamma * density * delta
! The factor delta (longterm activity factor) applies to isoprene emission
! from deciduous plants only. We neglect this factor at the present time.
! This factor is discussed in Guenther (1997).
! Subroutine written to operate at the patch level.
! IN FINAL IMPLEMENTATION, REMEMBER:
! 1. may wish to call this routine only as freq. as rad. calculations
! 2. may wish to place epsilon values directly in pft-physiology file
! Output: vocflx(nvoc) !VOC flux [ug C m-2 h-1]
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nvoc
    use clm_atmlnd   , only : clm_a2l
    use shr_const_mod, only : SHR_CONST_RGAS
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(num_nolakep) ! pft filter for non-lake points
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
! 2/1/02, Peter Thornton: migration to new data structure
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pgridcell(:)     ! gridcell index of corresponding pft
    integer , pointer :: ivt(:)           ! pft vegetation type for current
    real(r8), pointer :: t_veg(:)         ! pft vegetation temperature (Kelvin)
    real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
    real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
    real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (visible only)
    real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation     (visible only)
    real(r8), pointer :: sla(:)           ! ecophys constant - specific leaf area [m2 leaf g-1 carbon]
!
! local pointers to original implicit out arrays
!
    real(r8), pointer :: vocflx(:,:)      ! VOC flux [ug C m-2 h-1]
    real(r8), pointer :: vocflx_tot(:)    ! VOC flux [ug C m-2 h-1]
    real(r8), pointer :: vocflx_1(:)      ! VOC flux(1) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_2(:)      ! VOC flux(2) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_3(:)      ! VOC flux(3) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_4(:)      ! VOC flux(4) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_5(:)      ! VOC flux(5) [ug C m-2 h-1]
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: fp,p,g,n         ! indices
    real(r8) :: epsilon(lbp:ubp) ! emission factor [ug g-1 h-1]
    real(r8) :: gamma(lbp:ubp)   ! activity factor (instantaneous light and temp. condns)
    real(r8) :: density          ! source density factor [g dry wgt foliar mass/m2 ground]
    real(r8) :: cl               ! temporary
    real(r8) :: ct               ! temporary
    real(r8) :: par              ! temporary
    real(r8) :: reciprod         ! temporary

! Constants

    real(r8), parameter :: R   = SHR_CONST_RGAS*0.001_r8 ! univ. gas constant [J K-1 mol-1]
    real(r8), parameter :: alpha = 0.0027_r8 ! empirical coefficient
    real(r8), parameter :: cl1 = 1.066_r8    ! empirical coefficient
    real(r8), parameter :: ct1 = 95000.0_r8  ! empirical coefficient [J mol-1]
    real(r8), parameter :: ct2 = 230000.0_r8 ! empirical coefficient [J mol-1]
    real(r8), parameter :: ct3 = 0.961_r8    ! empirical coefficient
    real(r8), parameter :: tm  = 314.0_r8    ! empirical coefficient [K]
    real(r8), parameter :: tstd = 303.0_r8   ! std temperature [K]
    real(r8), parameter :: bet = 0.09_r8     ! beta empirical coefficient [K-1]

! These are the values from version of genesis-ibis / 1000.
! With DGVM defined, use LPJ's sla [m2 leaf g-1 carbon]
! Divide by 2 in the equation to get dry weight foliar mass from grams carbon

    real(r8) :: hardwire_sla(0:16)
    real(r8) :: slarea(lbp:ubp)           ! Specific leaf areas [m2 leaf g-1 carbon]
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    forc_solad => clm_a2l%forc_solad
    forc_solai => clm_a2l%forc_solai

    ! Assign local pointers to derived subtypes components (pft-level)

    pgridcell  => clm3%g%l%c%p%gridcell
    ivt        => clm3%g%l%c%p%itype
    t_veg      => clm3%g%l%c%p%pes%t_veg
    fsun       => clm3%g%l%c%p%pps%fsun
    elai       => clm3%g%l%c%p%pps%elai
    vocflx     => clm3%g%l%c%p%pvf%vocflx
    vocflx_tot => clm3%g%l%c%p%pvf%vocflx_tot
    vocflx_1   => clm3%g%l%c%p%pvf%vocflx_1
    vocflx_2   => clm3%g%l%c%p%pvf%vocflx_2
    vocflx_3   => clm3%g%l%c%p%pvf%vocflx_3
    vocflx_4   => clm3%g%l%c%p%pvf%vocflx_4
    vocflx_5   => clm3%g%l%c%p%pvf%vocflx_5
    sla        => pftcon%sla

#if (!defined DGVM)
    hardwire_sla( 0) = 0._r8
    hardwire_sla( 1) = 0.0125_r8 !needleleaf
    hardwire_sla( 2) = 0.0125_r8 !Gordon Bonan suggests NET = 0.0076
    hardwire_sla( 3) = 0.0125_r8 !Gordon Bonan suggests NDT = 0.0200
    hardwire_sla( 4) = 0.0250_r8 !broadleaf
    hardwire_sla( 5) = 0.0250_r8 !Gordon Bonan suggests BET = 0.0178
    hardwire_sla( 6) = 0.0250_r8 !Gordon Bonan suggests BDT = 0.0274
    hardwire_sla( 7) = 0.0250_r8
    hardwire_sla( 8) = 0.0250_r8
    hardwire_sla( 9) = 0.0250_r8
    hardwire_sla(10) = 0.0250_r8
    hardwire_sla(11) = 0.0250_r8
    hardwire_sla(12) = 0.0200_r8 !grass
    hardwire_sla(13) = 0.0200_r8
    hardwire_sla(14) = 0.0200_r8
    hardwire_sla(15) = 0.0200_r8
    hardwire_sla(16) = 0.0200_r8 !numpft = 16
#endif

    ! Determine specific leaf array
!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       g = pgridcell(p)

#if (!defined DGVM)
       slarea(p) = hardwire_sla(ivt(p))
#else
       slarea(p) = sla(ivt(p))
#endif
    end do

    ! Begin loop through voc species

    do n = 1, nvoc
       select case (n)

       case (1)
!dir$ concurrent
!cdir nodep
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             g = pgridcell(p)

             ! epsilon: use values from table 3 in Guenther (1997) which originate in
             ! -------  Guenther et al. (1995). In the comments below, I mention the pft
             !          category as described in table 3. Some values were taken directly
             !          from Guenther et al. (1995). Units: [ug g-1 h-1]
             !          Values were updated on 1/2002 (Guenther, personal communication)

             epsilon(p) = 0._r8

             ! isoprenes:
             if (ivt(p) == 1) then       !needleleaf evergreen temperate
                epsilon(p) = 2._r8
             else if (ivt(p) == 2) then  !needleleaf evergreen boreal
                epsilon(p) = 4._r8
             else if (ivt(p) == 4) then  !broadleaf evergreen tropical
                epsilon(p) = 24._r8
             else if (ivt(p) >= 5 .and. ivt(p) <= 11) then !other woody veg
                epsilon(p) = 24._r8
             end if
          end do

       case (2)
!dir$ concurrent
!cdir nodep
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             g = pgridcell(p)

             epsilon(p) = 0._r8
             ! monoterpenes:
             if (ivt(p) >= 1 .and. ivt(p) <= 2) then        !needleleaf evergreen
                epsilon(p) = 2.0_r8
             else if (ivt(p) == 3) then                     !needleleaf deciduous
                epsilon(p) = 1.6_r8
             else if (ivt(p) == 4) then                     !broadleaf everg trop
                epsilon(p) = 0.4_r8
             else if (ivt(p) >= 5 .and. ivt(p) <= 11) then  !other woody veg
                epsilon(p) = 0.8_r8
             else if (ivt(p) >= 12 .and. ivt(p) <= 17) then !grass & crop
                epsilon(p) = 0.1_r8
             end if
          end do

       case (3)
!dir$ concurrent
!cdir nodep
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             g = pgridcell(p)

             ! other VOCs (OVOCs)
             epsilon(p) = 1.0_r8                 !Guenther (personal communication)
          end do

       case (4)
!dir$ concurrent
!cdir nodep
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             g = pgridcell(p)

             ! other reactive VOCs (ORVOCs)
             epsilon(p) = 1.0_r8                 !Guenther (personal communication)
          end do

       case (5)
!dir$ concurrent
!cdir nodep
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             g = pgridcell(p)

             ! CO
             epsilon(p) = 0.3_r8                 !Guenther (personal communication)
          end do

       case default

          write(6,*)'only nvocs up to index 5 are currently supported'
          call endrun()

       end select

       select case (n)

       case (1)

!dir$ concurrent
!cdir nodep
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             g = pgridcell(p)

             ! gamma: Activity factor. Units [dimensionless]

             ! isoprenes:
             ! scale total incident par by fraction of sunlit leaves (added on 1/2002)
             ! multiply w/m2 by 4.6 to get umol/m2/s for par (added 8/14/02)
             ! got this value from subr. Stomata

             reciprod = 1._r8 / (R * t_veg(p) * tstd)
             ct = exp(ct1 * (t_veg(p) - tstd) * reciprod) / &
                  (ct3 + exp(ct2 * (t_veg(p) - tm) * reciprod))

             par = (forc_solad(g,1) + fsun(p) * forc_solai(g,1)) * 4.6_r8
             cl = alpha * cl1 * par * (1._r8 + alpha * alpha * par * par)**(-0.5_r8)
             gamma(p) = cl * ct !gamma = 1 under std temp & light condns

             par = ((1._r8 - fsun(p)) * forc_solai(g,1)) * 4.6_r8
             cl = alpha * cl1 * par * (1._r8 + alpha * alpha * par * par)**(-0.5_r8)
             gamma(p) = gamma(p) + cl * ct !gamma(sun) + gamma(sha)
          end do

       case (2,3,4,5)

!dir$ concurrent
!cdir nodep
          do fp = 1,num_nolakep
             p = filter_nolakep(fp)
             g = pgridcell(p)

             ! monoterpenes, OVOCs, and ORVOCs (Guenther, 1999 and 1995):
             gamma(p) = exp(bet * (t_veg(p) - tstd))
          end do

       end select

!dir$ concurrent
!cdir nodep
       do fp = 1,num_nolakep
          p = filter_nolakep(fp)
          g = pgridcell(p)

          ! density: Source density factor [g dry weight foliar mass m-2 ground]

          if (ivt(p) > 0) then
             density = elai(p) / (slarea(p) * 0.5_r8)
          else
             density = 0._r8
          end if

          ! calculate the voc flux

          vocflx(p,n) = epsilon(p) * gamma(p) * density
       end do   ! end pft loop

    end do   ! end voc species loop

    ! Calculate total voc flux and individual components for history output

!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       vocflx_tot(p) = 0._r8
    end do
    do n = 1, nvoc
!dir$ concurrent
!cdir nodep
       do fp = 1,num_nolakep
          p = filter_nolakep(fp)
          vocflx_tot(p) = vocflx_tot(p) + vocflx(p,n)
       end do
    end do
!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       vocflx_1(p) = vocflx(p,1)
       vocflx_2(p) = vocflx(p,2)
       vocflx_3(p) = vocflx(p,3)
       vocflx_4(p) = vocflx(p,4)
       vocflx_5(p) = vocflx(p,5)
    end do

  end subroutine VOCEmission

#endif

end module VOCEmissionMod
