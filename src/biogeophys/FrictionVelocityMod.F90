#include <misc.h>
#include <preproc.h>

module FrictionVelocityMod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: FrictionVelocityMod
!
! !DESCRIPTION:
! Calculation of the friction velocity, relation for potential
! temperature and humidity profiles of surface boundary layer.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: FrictionVelocity       ! Calculate friction velocity
  public :: MoninObukIni           ! Initialization of the Monin-Obukhov length
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: StabilityFunc1        ! Stability function for rib < 0.
  private :: StabilityFunc2        ! Stability function for rib < 0.
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FrictionVelocity
!
! !INTERFACE:
  subroutine FrictionVelocity(lbp, ubp, fn, filterp, &
                              displa, z0m, z0h, z0q, &
                              obu, iter, ur, um, ustar, &
                              temp1, temp2, temp12m, temp22m, fm)
!
! !DESCRIPTION:
! Calculation of the friction velocity, relation for potential
! temperature and humidity profiles of surface boundary layer.
! The scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
! Vol. 11, 2628-2644.
!
! !USES:
   use clmtype
   use clm_atmlnd, only : clm_a2l
   use clm_varcon, only : vkc
!
! !ARGUMENTS:
   implicit none
   integer , intent(in)  :: lbp, ubp         ! pft array bounds
   integer , intent(in)  :: fn               ! number of filtered pft elements
   integer , intent(in)  :: filterp(fn)      ! pft filter
   real(r8), intent(in)  :: displa(lbp:ubp)  ! displacement height (m)
   real(r8), intent(in)  :: z0m(lbp:ubp)     ! roughness length over vegetation, momentum [m]
   real(r8), intent(in)  :: z0h(lbp:ubp)     ! roughness length over vegetation, sensible heat [m]
   real(r8), intent(in)  :: z0q(lbp:ubp)     ! roughness length over vegetation, latent heat [m]
   real(r8), intent(in)  :: obu(lbp:ubp)     ! monin-obukhov length (m)
   integer,  intent(in)  :: iter             ! iteration number
   real(r8), intent(in)  :: ur(lbp:ubp)      ! wind speed at reference height [m/s]
   real(r8), intent(in)  :: um(lbp:ubp)      ! wind speed including the stablity effect [m/s]
   real(r8), intent(out) :: ustar(lbp:ubp)   ! friction velocity [m/s]
   real(r8), intent(out) :: temp1(lbp:ubp)   ! relation for potential temperature profile
   real(r8), intent(out) :: temp12m(lbp:ubp) ! relation for potential temperature profile applied at 2-m
   real(r8), intent(out) :: temp2(lbp:ubp)   ! relation for specific humidity profile
   real(r8), intent(out) :: temp22m(lbp:ubp) ! relation for specific humidity profile applied at 2-m
   real(r8), intent(inout) :: fm(lbp:ubp)    ! needed for DGVM only to diagnose 10m wind
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! Added arguments to eliminate passing clm derived type into this function.
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   integer , pointer :: pgridcell(:)   ! pft's gridcell index
   real(r8), pointer :: forc_hgt(:)    ! atmospheric reference height (m)
   real(r8), pointer :: forc_hgt_u(:)  ! observational height of wind [m]
   real(r8), pointer :: forc_hgt_t(:)  ! observational height of temperature [m]
   real(r8), pointer :: forc_hgt_q(:)  ! observational height of humidity [m]
!
! local pointers to implicit out arguments
!
   real(r8), pointer :: u10(:)         ! 10-m wind (m/s) (for dust model)
   real(r8), pointer :: fv(:)          ! friction velocity (m/s) (for dust model)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
   real(r8), parameter :: zetam = 1.574_r8 ! transition point of flux-gradient relation (wind profile)
   real(r8), parameter :: zetat = 0.465_r8 ! transition point of flux-gradient relation (temp. profile)
   integer :: f                         ! pft-filter index
   integer :: p                         ! pft index
   integer :: g                         ! gridcell index
   real(r8):: zldis(lbp:ubp)            ! reference height "minus" zero displacement heght [m]
   real(r8):: zeta(lbp:ubp)             ! dimensionless height used in Monin-Obukhov theory
#if (defined DGVM) || (defined DUST)
   real(r8) :: tmp1,tmp2,tmp3,tmp4      ! Used to diagnose the 10 meter wind
   real(r8) :: fmnew                    ! Used to diagnose the 10 meter wind
   real(r8) :: fm10                     ! Used to diagnose the 10 meter wind
   real(r8) :: zeta10                   ! Used to diagnose the 10 meter wind
#endif
!------------------------------------------------------------------------------

   ! Assign local pointers to derived type members (gridcell-level)

   pgridcell  => clm3%g%l%c%p%gridcell
   forc_hgt   => clm_a2l%forc_hgt
   forc_hgt_u => clm_a2l%forc_hgt_u
   forc_hgt_t => clm_a2l%forc_hgt_t
   forc_hgt_q => clm_a2l%forc_hgt_q

   ! Assign local pointers to derived type members (pft-level)

   u10        => clm3%g%l%c%p%pps%u10
   fv         => clm3%g%l%c%p%pps%fv

   ! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

#if (!defined PERGRO)

!dir$ concurrent
!cdir nodep
   do f = 1, fn
      p = filterp(f)
      g = pgridcell(p)

      ! Wind profile

      zldis(p) = forc_hgt_u(g)-displa(p)
      zeta(p) = zldis(p)/obu(p)
      if (zeta(p) < -zetam) then
         ustar(p) = vkc*um(p)/(log(-zetam*obu(p)/z0m(p))&
              - StabilityFunc1(-zetam) &
              + StabilityFunc1(z0m(p)/obu(p)) &
              + 1.14_r8*((-zeta(p))**0.333_r8-(zetam)**0.333_r8))
      else if (zeta(p) < 0._r8) then
         ustar(p) = vkc*um(p)/(log(zldis(p)/z0m(p))&
              - StabilityFunc1(zeta(p))&
              + StabilityFunc1(z0m(p)/obu(p)))
      else if (zeta(p) <=  1._r8) then
         ustar(p) = vkc*um(p)/(log(zldis(p)/z0m(p)) + 5._r8*zeta(p) -5._r8*z0m(p)/obu(p))
      else
         ustar(p) = vkc*um(p)/(log(obu(p)/z0m(p))+5._r8-5._r8*z0m(p)/obu(p) &
              +(5._r8*log(zeta(p))+zeta(p)-1._r8))
      end if

      ! Temperature profile

      zldis(p) = forc_hgt_t(g)-displa(p)
      zeta(p) = zldis(p)/obu(p)
      if (zeta(p) < -zetat) then
         temp1(p) = vkc/(log(-zetat*obu(p)/z0h(p))&
              - StabilityFunc2(-zetat) &
              + StabilityFunc2(z0h(p)/obu(p)) &
              + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
      else if (zeta(p) < 0._r8) then
         temp1(p) = vkc/(log(zldis(p)/z0h(p)) &
              - StabilityFunc2(zeta(p)) &
              + StabilityFunc2(z0h(p)/obu(p)))
      else if (zeta(p) <=  1._r8) then
         temp1(p) = vkc/(log(zldis(p)/z0h(p)) + 5._r8*zeta(p) - 5._r8*z0h(p)/obu(p))
      else
         temp1(p) = vkc/(log(obu(p)/z0h(p)) + 5._r8 - 5._r8*z0h(p)/obu(p) &
              + (5._r8*log(zeta(p))+zeta(p)-1._r8))
      end if

      ! Humidity profile

      if (forc_hgt_q(g) == forc_hgt_t(g) .and. z0q(p) == z0h(p)) then
         temp2(p) = temp1(p)
      else
         zldis(p) = forc_hgt_q(g)-displa(p)
         zeta(p) = zldis(p)/obu(p)
         if (zeta(p) < -zetat) then
            temp2(p) = vkc/(log(-zetat*obu(p)/z0q(p)) &
                 - StabilityFunc2(-zetat) &
                 + StabilityFunc2(z0q(p)/obu(p)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
         else if (zeta(p) < 0._r8) then
            temp2(p) = vkc/(log(zldis(p)/z0q(p)) &
                 - StabilityFunc2(zeta(p)) &
                 + StabilityFunc2(z0q(p)/obu(p)))
         else if (zeta(p) <=  1._r8) then
            temp2(p) = vkc/(log(zldis(p)/z0q(p)) + 5._r8*zeta(p)-5._r8*z0q(p)/obu(p))
         else
            temp2(p) = vkc/(log(obu(p)/z0q(p)) + 5._r8 - 5._r8*z0q(p)/obu(p) &
                 + (5._r8*log(zeta(p))+zeta(p)-1._r8))
         end if
      endif

      ! Temperature profile applied at 2-m

      zldis(p) = 2.0_r8 + z0h(p)
      zeta(p) = zldis(p)/obu(p)
      if (zeta(p) < -zetat) then
         temp12m(p) = vkc/(log(-zetat*obu(p)/z0h(p))&
              - StabilityFunc2(-zetat) &
              + StabilityFunc2(z0h(p)/obu(p)) &
              + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
      else if (zeta(p) < 0._r8) then
         temp12m(p) = vkc/(log(zldis(p)/z0h(p)) &
              - StabilityFunc2(zeta(p))  &
              + StabilityFunc2(z0h(p)/obu(p)))
      else if (zeta(p) <=  1._r8) then
         temp12m(p) = vkc/(log(zldis(p)/z0h(p)) + 5._r8*zeta(p) - 5._r8*z0h(p)/obu(p))
      else
         temp12m(p) = vkc/(log(obu(p)/z0h(p)) + 5._r8 - 5._r8*z0h(p)/obu(p) &
              + (5._r8*log(zeta(p))+zeta(p)-1._r8))
      end if

      ! Humidity profile applied at 2-m

      if (z0q(p) == z0h(p)) then
         temp22m(p) = temp12m(p)
      else
         zldis(p) = 2.0_r8 + z0q(p)
         zeta(p) = zldis(p)/obu(p)
         if (zeta(p) < -zetat) then
            temp22m(p) = vkc/(log(-zetat*obu(p)/z0q(p)) - &
                 StabilityFunc2(-zetat) + StabilityFunc2(z0q(p)/obu(p)) &
                 + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
         else if (zeta(p) < 0._r8) then
            temp22m(p) = vkc/(log(zldis(p)/z0q(p)) - &
                 StabilityFunc2(zeta(p))+StabilityFunc2(z0q(p)/obu(p)))
         else if (zeta(p) <=  1._r8) then
            temp22m(p) = vkc/(log(zldis(p)/z0q(p)) + 5._r8*zeta(p)-5._r8*z0q(p)/obu(p))
         else
            temp22m(p) = vkc/(log(obu(p)/z0q(p)) + 5._r8 - 5._r8*z0q(p)/obu(p) &
                 + (5._r8*log(zeta(p))+zeta(p)-1._r8))
         end if
      end if

#if (defined DGVM) || (defined DUST)
      ! diagnose 10-m wind for dust model (dstmbl.F)
      ! Notes from C. Zender's dst.F:
      ! According to Bon96 p. 62, the displacement height d (here displa) is
      ! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
      ! Therefore d <= 0.034*z1 and may safely be neglected.
      ! Code from LSM routine SurfaceTemperature was used to obtain u10

      zldis(p) = forc_hgt_u(g)-displa(p)
      zeta(p) = zldis(p)/obu(p)
      if (min(zeta(p), 1._r8) < 0._r8) then
         tmp1 = (1._r8 - 16._r8*min(zeta(p),1._r8))**0.25_r8
         tmp2 = log((1._r8+tmp1*tmp1)/2._r8)
         tmp3 = log((1._r8+tmp1)/2._r8)
         fmnew = 2._r8*tmp3 + tmp2 - 2._r8*atan(tmp1) + 1.5707963_r8
      else
         fmnew = -5._r8*min(zeta(p),1._r8)
      endif
      if (iter == 1) then
         fm(p) = fmnew
      else
         fm(p) = 0.5_r8 * (fm(p)+fmnew)
      end if
      zeta10 = min(10._r8/obu(p), 1._r8)
      if (zeta(p) == 0._r8) zeta10 = 0._r8
      if (zeta10 < 0._r8) then
         tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
         tmp2 = log((1.0_r8 + tmp1*tmp1)/2.0_r8)
         tmp3 = log((1.0_r8 + tmp1)/2.0_r8)
         fm10 = 2.0_r8*tmp3 + tmp2 - 2.0_r8*atan(tmp1) + 1.5707963_r8
      else                ! not stable
         fm10 = -5.0_r8 * zeta10
      end if
      tmp4 = log(forc_hgt(g) / 10._r8)
      u10(p) = ur(p) - ustar(p)/vkc * (tmp4 - fm(p) + fm10)
      fv(p)  = ustar(p)
#endif

   end do
#endif


#if (defined PERGRO)

   !===============================================================================
   ! The following only applies when PERGRO is defined
   !===============================================================================

!dir$ concurrent
!cdir nodep
   do f = 1, fn
      p = filterp(f)
      g = pgridcell(p)

      zldis(p) = forc_hgt_u(g)-displa(p)
      zeta(p) = zldis(p)/obu(p)
      if (zeta(p) < -zetam) then           ! zeta < -1
         ustar(p) = vkc * um(p) / log(-zetam*obu(p)/z0m(p))
      else if (zeta(p) < 0._r8) then         ! -1 <= zeta < 0
         ustar(p) = vkc * um(p) / log(zldis(p)/z0m(p))
      else if (zeta(p) <= 1._r8) then        !  0 <= ztea <= 1
         ustar(p)=vkc * um(p)/log(zldis(p)/z0m(p))
      else                             !  1 < zeta, phi=5+zeta
         ustar(p)=vkc * um(p)/log(obu(p)/z0m(p))
      endif

      zldis(p) = forc_hgt_t(g)-displa(p)
      zeta(p) = zldis(p)/obu(p)
      if (zeta(p) < -zetat) then
         temp1(p)=vkc/log(-zetat*obu(p)/z0h(p))
      else if (zeta(p) < 0._r8) then
         temp1(p)=vkc/log(zldis(p)/z0h(p))
      else if (zeta(p) <= 1._r8) then
         temp1(p)=vkc/log(zldis(p)/z0h(p))
      else
         temp1(p)=vkc/log(obu(p)/z0h(p))
      end if

      zldis(p) = forc_hgt_q(g)-displa(p)
      zeta(p) = zldis(p)/obu(p)
      if (zeta(p) < -zetat) then
         temp2(p)=vkc/log(-zetat*obu(p)/z0q(p))
      else if (zeta(p) < 0._r8) then
         temp2(p)=vkc/log(zldis(p)/z0q(p))
      else if (zeta(p) <= 1._r8) then
         temp2(p)=vkc/log(zldis(p)/z0q(p))
      else
         temp2(p)=vkc/log(obu(p)/z0q(p))
      end if

      zldis(p) = 2.0_r8 + z0h(p)
      zeta(p) = zldis(p)/obu(p)
      if (zeta(p) < -zetat) then
         temp12m(p)=vkc/log(-zetat*obu(p)/z0h(p))
      else if (zeta(p) < 0._r8) then
         temp12m(p)=vkc/log(zldis(p)/z0h(p))
      else if (zeta(p) <= 1._r8) then
         temp12m(p)=vkc/log(zldis(p)/z0h(p))
      else
         temp12m(p)=vkc/log(obu(p)/z0h(p))
      end if

      zldis(p) = 2.0_r8 + z0q(p)
      zeta(p) = zldis(p)/obu(p)
      if (zeta(p) < -zetat) then
         temp22m(p)=vkc/log(-zetat*obu(p)/z0q(p))
      else if (zeta(p) < 0._r8) then
         temp22m(p)=vkc/log(zldis(p)/z0q(p))
      else if (zeta(p) <= 1._r8) then
         temp22m(p)=vkc/log(zldis(p)/z0q(p))
      else
         temp22m(p)=vkc/log(obu(p)/z0q(p))
      end if
#if (defined DGVM) || (defined DUST)
      ! diagnose 10-m wind for dust model (dstmbl.F)
      ! Notes from C. Zender's dst.F:
      ! According to Bon96 p. 62, the displacement height d (here displa) is
      ! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
      ! Therefore d <= 0.034*z1 and may safely be neglected.
      ! Code from LSM routine SurfaceTemperature was used to obtain u10

      zldis(p) = forc_hgt_u(g)-displa(p)
      zeta(p) = zldis(p)/obu(p)
      if (min(zeta(p), 1._r8) < 0._r8) then
         tmp1 = (1._r8 - 16._r8*min(zeta(p),1._r8))**0.25_r8
         tmp2 = log((1._r8+tmp1*tmp1)/2._r8)
         tmp3 = log((1._r8+tmp1)/2._r8)
         fmnew = 2._r8*tmp3 + tmp2 - 2._r8*atan(tmp1) + 1.5707963_r8
      else
         fmnew = -5._r8*min(zeta(p),1._r8)
      endif
      if (iter == 1) then
         fm(p) = fmnew
      else
         fm(p) = 0.5_r8 * (fm(p)+fmnew)
      end if
      zeta10 = min(10._r8/obu(p), 1._r8)
      if (zeta(p) == 0._r8) zeta10 = 0._r8
      if (zeta10 < 0._r8) then
         tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
         tmp2 = log((1.0_r8 + tmp1*tmp1)/2.0_r8)
         tmp3 = log((1.0_r8 + tmp1)/2.0_r8)
         fm10 = 2.0_r8*tmp3 + tmp2 - 2.0_r8*atan(tmp1) + 1.5707963_r8
      else                ! not stable
         fm10 = -5.0_r8 * zeta10
      end if
      tmp4 = log(forc_hgt(g) / 10._r8)
      u10(p) = ur(p) - ustar(p)/vkc * (tmp4 - fm(p) + fm10)
      fv(p)  = ustar(p)
#endif
   end do

#endif

   end subroutine FrictionVelocity

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StabilityFunc
!
! !INTERFACE:
   real(r8) function StabilityFunc1(zeta)
!
! !DESCRIPTION:
! Stability function for rib < 0.
!
! !USES:
      use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
!
! !CALLED FROM:
! subroutine FrictionVelocity in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
      real(r8) :: chik, chik2
!------------------------------------------------------------------------------

      chik2 = sqrt(1._r8-16._r8*zeta)
      chik = sqrt(chik2)
      StabilityFunc1 = 2._r8*log((1._r8+chik)*0.5_r8) &
           + log((1._r8+chik2)*0.5_r8)-2._r8*atan(chik)+SHR_CONST_PI*0.5_r8

    end function StabilityFunc1

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StabilityFunc2
!
! !INTERFACE:
   real(r8) function StabilityFunc2(zeta)
!
! !DESCRIPTION:
! Stability function for rib < 0.
!
! !USES:
     use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
     implicit none
     real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
!
! !CALLED FROM:
! subroutine FrictionVelocity in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
     real(r8) :: chik2
!------------------------------------------------------------------------------

     chik2 = sqrt(1._r8-16._r8*zeta)
     StabilityFunc2 = 2._r8*log((1._r8+chik2)*0.5_r8)

   end function StabilityFunc2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: MoninObukIni
!
! !INTERFACE:
  subroutine MoninObukIni (ur, thv, dthv, zldis, z0m, um, obu)
!
! !DESCRIPTION:
! Initialization of the Monin-Obukhov length.
! The scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
! Vol. 11, 2628-2644.
!
! !USES:
    use clm_varcon, only : grav
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: ur    ! wind speed at reference height [m/s]
    real(r8), intent(in)  :: thv   ! virtual potential temperature (kelvin)
    real(r8), intent(in)  :: dthv  ! diff of vir. poten. temp. between ref. height and surface
    real(r8), intent(in)  :: zldis ! reference height "minus" zero displacement heght [m]
    real(r8), intent(in)  :: z0m   ! roughness length, momentum [m]
    real(r8), intent(out) :: um    ! wind speed including the stability effect [m/s]
    real(r8), intent(out) :: obu   ! monin-obukhov length (m)
!
! !CALLED FROM:
! subroutine BareGroundFluxes in module BareGroundFluxesMod.F90
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod.F90
! subroutine CanopyFluxes in module CanopyFluxesMod.F90
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8) :: wc    ! convective velocity [m/s]
    real(r8) :: rib   ! bulk Richardson number
    real(r8) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: ustar ! friction velocity [m/s]
!-----------------------------------------------------------------------

    ! Initial values of u* and convective velocity

    ustar=0.06_r8
    wc=0.5_r8
    if (dthv >= 0._r8) then
       um=max(ur,0.1_r8)
    else
       um=sqrt(ur*ur+wc*wc)
    endif

    rib=grav*zldis*dthv/(thv*um*um)
#if (defined PERGRO)
    rib = 0._r8
#endif

    if (rib >= 0._r8) then      ! neutral or stable
       zeta = rib*log(zldis/z0m)/(1._r8-5._r8*min(rib,0.19_r8))
       zeta = min(2._r8,max(zeta,0.01_r8 ))
    else                     ! unstable
       zeta=rib*log(zldis/z0m)
       zeta = max(-100._r8,min(zeta,-0.01_r8 ))
    endif

    obu=zldis/zeta

  end subroutine MoninObukIni

end module FrictionVelocityMod
