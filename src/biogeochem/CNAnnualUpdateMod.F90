#include <misc.h>
#include <preproc.h>

module CNAnnualUpdateMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNAnnualUpdateMod
!
! !DESCRIPTION:
! Module for updating annual summation variables
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
    public:: CNAnnualUpdate
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
! !IROUTINE: CNAnnualUpdate
!
! !INTERFACE:
subroutine CNAnnualUpdate(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update annual summation variables
!
! !USES:
   use clmtype
   use clm_varctl, only: irad
   use clm_time_manager, only: get_step_size
   use pft2colMod, only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 10/1/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: pcolumn(:)               ! index into column level
                                                 ! quantities
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: annsum_counter(:)        ! seconds since last annual accumulator turnover
   real(r8), pointer :: tempsum_plant_ndemand(:) ! temporary annual sum of plant_ndemand
   real(r8), pointer :: annsum_plant_ndemand(:)  ! annual sum of plant_ndemand
   real(r8), pointer :: tempsum_retransn(:)      ! temporary annual sum of N retranslocation
   real(r8), pointer :: annsum_retransn(:)       ! annual sum of N retranslocation
   real(r8), pointer :: tempavg_t2m(:)           ! temporary average 2m air temperature (K)
   real(r8), pointer :: annavg_t2m(:)            ! annual average 2m air temperature (K)
   real(r8), pointer :: tempsum_npp(:)           ! temporary sum NPP (gC/m2/yr)
   real(r8), pointer :: annsum_npp(:)            ! annual sum NPP (gC/m2/yr)
   real(r8), pointer :: cannsum_npp(:)           ! column annual sum NPP (gC/m2/yr)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p          ! indices
   integer :: fp,fc        ! lake filter indices
   integer :: dtime        ! time step (s)
   real(r8):: dt           ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays
   annsum_counter        => clm3%g%l%c%cps%annsum_counter
   tempsum_plant_ndemand => clm3%g%l%c%p%pepv%tempsum_plant_ndemand
   annsum_plant_ndemand  => clm3%g%l%c%p%pepv%annsum_plant_ndemand
   tempsum_retransn      => clm3%g%l%c%p%pepv%tempsum_retransn
   annsum_retransn       => clm3%g%l%c%p%pepv%annsum_retransn
   tempavg_t2m           => clm3%g%l%c%p%pepv%tempavg_t2m
   annavg_t2m            => clm3%g%l%c%p%pepv%annavg_t2m
   tempsum_npp           => clm3%g%l%c%p%pepv%tempsum_npp
   annsum_npp            => clm3%g%l%c%p%pepv%annsum_npp
   cannsum_npp           => clm3%g%l%c%cps%cannsum_npp
   pcolumn               => clm3%g%l%c%p%column

   ! set time steps
   dtime = get_step_size()
   dt = float(irad)*dtime

   ! column loop
!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      annsum_counter(c) = annsum_counter(c) + dt
   end do

   ! pft loop
!dir$ concurrent
!cdir nodep
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)
      if (annsum_counter(c) > 365._r8 * 86400._r8) then
         ! update annual plant ndemand accumulator
         annsum_plant_ndemand(p)  = tempsum_plant_ndemand(p)
         tempsum_plant_ndemand(p) = 0._r8

         ! update annual total N retranslocation accumulator
         annsum_retransn(p)  = tempsum_retransn(p)
         tempsum_retransn(p) = 0._r8

         ! update annual average 2m air temperature accumulator
         annavg_t2m(p)  = tempavg_t2m(p)
         tempavg_t2m(p) = 0._r8

         ! update annual NPP accumulator, convert to annual total
         annsum_npp(p) = tempsum_npp(p) * dt
         tempsum_npp(p) = 0._r8
      end if
   end do

   ! use p2c routine to get selected column-average pft-level fluxes and states
   call p2c(num_soilc, filter_soilc, annsum_npp, cannsum_npp)

   ! column loop
!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      if (annsum_counter(c) > 365._r8 * 86400._r8) annsum_counter(c) = 0._r8
   end do

end subroutine CNAnnualUpdate
!-----------------------------------------------------------------------

end module CNAnnualUpdateMod
