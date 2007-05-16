#include <misc.h>
#include <preproc.h>

module CASASummaryMod

#if (defined CASA)
#if (defined CLAMP)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CASASummaryMod
!
! !DESCRIPTION:
! Module for carbon summary calculations for CASA'
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: CASASummary
!
! !REVISION HISTORY:
! 22 Sept 2006: Created by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CASASummary
!
! !INTERFACE:
subroutine CASASummary(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Perform pft carbon summary calculations
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use CASAMod, only: gppfact, LEAF, WOOD, FROOT, CWD, SURFMET, SURFSTR, &
                      SURFMIC, SOILMET, SOILSTR, SOILMIC, SLOW, PASSIVE, &
                      CWD_TYPE, LITTER_TYPE, SOIL_TYPE
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbp, ubp       ! pft bounds
   integer, intent(in) :: num_soilp      ! number of soil points in pft filter
   integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 22 Sept 2006: Created by Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   real(r8), pointer :: fpsn(:)              !photosynthesis (umol CO2 /m**2 /s)
   real(r8), pointer :: Tpool_C(:,:)         ! Total C by pool
   real(r8), pointer :: Resp_C(:,:)          ! Respired C by pool
   real(r8), pointer :: Closs(:,:)           ! Lost C by pool
   real(r8), pointer :: Ctrans(:,:)          ! Transferred C by pool type
   real(r8), pointer :: Cflux(:)             ! C flux
   real(r8), pointer :: livefr(:,:)          ! Live fraction
   real(r8), pointer :: fnpp(:)              ! NPP (gC/m2/sec)
   real(r8), pointer :: co2flux(:)           ! net CO2 flux (gC/m2/s) [+= atm]
!
!
! local pointers to implicit in/out scalars
!
!
! local pointers to implicit out scalars
   real(r8), pointer :: casa_agnpp(:)        ! above-ground net primary production [gC/m2/s]
   real(r8), pointer :: casa_ar(:)           ! autotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_bgnpp(:)        ! below-ground net primary production [gC/m2/s]
   real(r8), pointer :: casa_cwdc(:)         ! coarse woody debris C [gC/m2]
   real(r8), pointer :: casa_cwdc_hr(:)      ! cwd heterotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_cwdc_loss(:)    ! cwd C loss [gC/m2/s]
   real(r8), pointer :: casa_frootc(:)       ! fine root C [gC/m2]
   real(r8), pointer :: casa_frootc_alloc(:) ! fine root C allocation [gC/m2/s]
   real(r8), pointer :: casa_frootc_loss(:)  ! fine root C loss [gC/m2/s]
   real(r8), pointer :: casa_gpp(:)          ! gross primary production [gC/m2/s]
   real(r8), pointer :: casa_hr(:)           ! total heterotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_leafc(:)        ! leaf C [gC/m2]
   real(r8), pointer :: casa_leafc_alloc(:)  ! leaf C allocation [gC/m2/s]
   real(r8), pointer :: casa_leafc_loss(:)   ! leaf C loss [gC/m2/s]
   real(r8), pointer :: casa_litterc(:)      ! total litter C (excluding cwd C) [gC/m2]
   real(r8), pointer :: casa_litterc_hr(:)   ! litter heterotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_litterc_loss(:) ! litter C loss [gC/m2/s]
   real(r8), pointer :: casa_nee(:)          ! net ecosystem exchange [gC/m2/s]
   real(r8), pointer :: casa_nep(:)          ! net ecosystem production [gC/m2/s]
   real(r8), pointer :: casa_npp(:)          ! net primary production [gC/m2/s]
   real(r8), pointer :: casa_soilc(:)        ! total soil organic matter C (excluding cwd and litter C) [gC/m2]
   real(r8), pointer :: casa_soilc_hr(:)     ! soil heterotrophic respiration [gC/m2/s]
   real(r8), pointer :: casa_soilc_loss(:)   ! total soil organic matter C loss [gC/m2/s]
   real(r8), pointer :: casa_woodc(:)        ! wood C [gC/m2]
   real(r8), pointer :: casa_woodc_alloc(:)  ! wood C allocation [gC/m2/s]
   real(r8), pointer :: casa_woodc_loss(:)   ! wood C loss [gC/m2/s]
!
!
! !OTHER LOCAL VARIABLES:
   integer :: f,p          ! indices
   real(r8):: dtime        ! land model time step (sec)

!EOP
!-----------------------------------------------------------------------

   ! assign local pointers at the pft level
   fpsn                           => clm3%g%l%c%p%pcf%fpsn
   Tpool_C                        => clm3%g%l%c%p%pps%Tpool_C
   Resp_C                         => clm3%g%l%c%p%pps%Resp_C
   Closs                          => clm3%g%l%c%p%pps%Closs
   Ctrans                         => clm3%g%l%c%p%pps%Ctrans
   Cflux                          => clm3%g%l%c%p%pps%Cflux  
   livefr                         => clm3%g%l%c%p%pps%livefr
   fnpp                           => clm3%g%l%c%p%pps%fnpp
   co2flux                        => clm3%g%l%c%p%pps%co2flux

   casa_agnpp                     => clm3%g%l%c%p%pps%casa_agnpp
   casa_ar                        => clm3%g%l%c%p%pps%casa_ar
   casa_bgnpp                     => clm3%g%l%c%p%pps%casa_bgnpp
   casa_cwdc                      => clm3%g%l%c%p%pps%casa_cwdc
   casa_cwdc_hr                   => clm3%g%l%c%p%pps%casa_cwdc_hr
   casa_cwdc_loss                 => clm3%g%l%c%p%pps%casa_cwdc_loss
   casa_frootc                    => clm3%g%l%c%p%pps%casa_frootc
   casa_frootc_alloc              => clm3%g%l%c%p%pps%casa_frootc_alloc
   casa_frootc_loss               => clm3%g%l%c%p%pps%casa_frootc_loss
   casa_gpp                       => clm3%g%l%c%p%pps%casa_gpp
   casa_hr                        => clm3%g%l%c%p%pps%casa_hr
   casa_leafc                     => clm3%g%l%c%p%pps%casa_leafc
   casa_leafc_alloc               => clm3%g%l%c%p%pps%casa_leafc_alloc
   casa_leafc_loss                => clm3%g%l%c%p%pps%casa_leafc_loss
   casa_litterc                   => clm3%g%l%c%p%pps%casa_litterc
   casa_litterc_hr                => clm3%g%l%c%p%pps%casa_litterc_hr
   casa_litterc_loss              => clm3%g%l%c%p%pps%casa_litterc_loss
   casa_nee                       => clm3%g%l%c%p%pps%casa_nee
   casa_nep                       => clm3%g%l%c%p%pps%casa_nep
   casa_npp                       => clm3%g%l%c%p%pps%casa_npp
   casa_soilc                     => clm3%g%l%c%p%pps%casa_soilc
   casa_soilc_hr                  => clm3%g%l%c%p%pps%casa_soilc_hr
   casa_soilc_loss                => clm3%g%l%c%p%pps%casa_soilc_loss
   casa_woodc                     => clm3%g%l%c%p%pps%casa_woodc
   casa_woodc_alloc               => clm3%g%l%c%p%pps%casa_woodc_alloc
   casa_woodc_loss                => clm3%g%l%c%p%pps%casa_woodc_loss

   ! Get step size
   dtime = get_step_size()

   ! pft loop
!dir$ concurrent
!cdir nodep
   do f = 1, num_soilp
      p = filter_soilp(f)

      ! calculate pft-level summary carbon fluxes and states

      ! autotrophic respiration
      casa_ar(p) = gppfact * fpsn(p) * 12.0_r8 * 1.e-6_r8 ! umolC/m2/s to gC/m2/s
      ! gross primary production
      casa_gpp(p) = fpsn(p) * 12.0_r8 * 1.e-6_r8 ! umolC/m2/s to gC/m2/s

      ! total heterotrophic respiration
      casa_hr(p) = Cflux(p)

      ! net ecosystem exchange
      casa_nee(p) = co2flux(p)

      ! net ecosystem production
      casa_nep(p) =  -1._r8 * co2flux(p)

      ! net primary production
      casa_npp(p) = fnpp(p)

      ! above-ground and below-ground NPP
      casa_agnpp(p) = (livefr(p,LEAF)  + 0.8_r8 * livefr(p,WOOD)) * fnpp(p)
      casa_bgnpp(p) = (livefr(p,FROOT) + 0.2_r8 * livefr(p,WOOD)) * fnpp(p)

      ! leaf C
      casa_leafc(p)       = Tpool_C(p,LEAF)
      casa_leafc_alloc(p) = livefr(p,LEAF) * fnpp(p) * dtime
      casa_leafc_loss(p)  = Closs(p,LEAF)

      ! wood C
      casa_woodc(p)       = Tpool_C(p,WOOD)
      casa_woodc_alloc(p) = livefr(p,WOOD) * fnpp(p) * dtime
      casa_woodc_loss(p)  = Closs(p,WOOD)

      ! fine root C
      casa_frootc(p)       = Tpool_C(p,FROOT)
      casa_frootc_alloc(p) = livefr(p,FROOT) * fnpp(p) * dtime
      casa_frootc_loss(p)  = Closs(p,FROOT)

      ! coarse woody debris C
      casa_cwdc(p)      = Tpool_C(p,CWD)
      casa_cwdc_hr(p)   = Resp_C(p,CWD)
      casa_cwdc_loss(p) = Ctrans(p,CWD_TYPE)

      ! litter C
      casa_litterc(p)      = Tpool_C(p,SURFMET) + Tpool_C(p,SURFSTR) + &
                             Tpool_C(p,SURFMIC) + Tpool_C(p,SOILMET) + &
                             Tpool_C(p,SOILSTR) + Tpool_C(p,SOILMIC)
      casa_litterc_hr(p)   = Resp_C(p,SURFMET) + Resp_C(p,SURFSTR) + &
                             Resp_C(p,SURFMIC) + Resp_C(p,SOILMET) + &
                             Resp_C(p,SOILSTR) + Resp_C(p,SOILMIC)
      casa_litterc_loss(p) = Ctrans(p,LITTER_TYPE)

      ! soil C
      casa_soilc(p)      = Tpool_C(p,SLOW) + Tpool_C(p,PASSIVE)
      casa_soilc_hr(p)   = Resp_C(p,SLOW) + Resp_C(p,PASSIVE)
      casa_soilc_loss(p) = Ctrans(p,SOIL_TYPE)

   end do  ! end of pfts loop

end subroutine CASASummary
!-----------------------------------------------------------------------

#endif
#endif

end module CASASummaryMod

