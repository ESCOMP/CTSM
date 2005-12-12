#include <misc.h>
#include <preproc.h>

module CNC13FluxMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: C13FluxMod
!
! !DESCRIPTION:
! Module for 13-carbon flux variable update, non-mortality fluxes.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use shr_sys_mod , only: shr_sys_flush
    use clm_varcon  , only: istsoil,c13ratio
    use spmdMod     , only: masterproc
    use clm_varpar  , only: nlevsoi
    implicit none
    save
    private
!
! !PUBLIC MEMBER FUNCTIONS:
    public:: C13Flux1
    private:: C13FluxCalc
!
! !REVISION HISTORY:
! 4/21/2005: Created by Peter Thornton and Neil Suits
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Flux1
!
! !INTERFACE:
subroutine C13Flux1(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, set the 13-carbon flux
! variables (except for gap-phase mortality and fire fluxes)
!
! !USES:
   use clmtype
   use time_manager, only: get_step_size
   use clm_varctl, only: irad
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p

   call C13FluxCalc(p%pc13f%leafc_xfer_to_leafc, p%pcf%leafc_xfer_to_leafc, &
                    p%pc13s%leafc_xfer, p%pcs%leafc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%frootc_xfer_to_frootc, p%pcf%frootc_xfer_to_frootc, &
                    p%pc13s%frootc_xfer, p%pcs%frootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livestemc_xfer_to_livestemc, p%pcf%livestemc_xfer_to_livestemc, &
                    p%pc13s%livestemc_xfer, p%pcs%livestemc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%deadstemc_xfer_to_deadstemc, p%pcf%deadstemc_xfer_to_deadstemc, &
                    p%pc13s%deadstemc_xfer, p%pcs%deadstemc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livecrootc_xfer_to_livecrootc, p%pcf%livecrootc_xfer_to_livecrootc, &
                    p%pc13s%livecrootc_xfer, p%pcs%livecrootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%deadcrootc_xfer_to_deadcrootc, p%pcf%deadcrootc_xfer_to_deadcrootc, &
                    p%pc13s%deadcrootc_xfer, p%pcs%deadcrootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%leafc_to_litter, p%pcf%leafc_to_litter, &
                    p%pc13s%leafc, p%pcs%leafc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%frootc_to_litter, p%pcf%frootc_to_litter, &
                    p%pc13s%frootc, p%pcs%frootc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livestemc_to_deadstemc, p%pcf%livestemc_to_deadstemc, &
                    p%pc13s%livestemc, p%pcs%livestemc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livecrootc_to_deadcrootc, p%pcf%livecrootc_to_deadcrootc, &
                    p%pc13s%livecrootc, p%pcs%livecrootc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%leaf_curmr, p%pcf%leaf_curmr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%froot_curmr, p%pcf%froot_curmr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livestem_curmr, p%pcf%livestem_curmr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livecroot_curmr, p%pcf%livecroot_curmr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%leaf_xsmr, p%pcf%leaf_xsmr, &
                    p%pc13s%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%froot_xsmr, p%pcf%froot_xsmr, &
                    p%pc13s%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livestem_xsmr, p%pcf%livestem_xsmr, &
                    p%pc13s%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livecroot_xsmr, p%pcf%livecroot_xsmr, &
                    p%pc13s%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_xsmrpool, p%pcf%cpool_to_xsmrpool, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_leafc, p%pcf%cpool_to_leafc, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_leafc_storage, p%pcf%cpool_to_leafc_storage, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_frootc, p%pcf%cpool_to_frootc, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_frootc_storage, p%pcf%cpool_to_frootc_storage, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_livestemc, p%pcf%cpool_to_livestemc, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_livestemc_storage, p%pcf%cpool_to_livestemc_storage, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_deadstemc, p%pcf%cpool_to_deadstemc, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_deadstemc_storage, p%pcf%cpool_to_deadstemc_storage, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_livecrootc, p%pcf%cpool_to_livecrootc, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_livecrootc_storage, p%pcf%cpool_to_livecrootc_storage, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_deadcrootc, p%pcf%cpool_to_deadcrootc, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_deadcrootc_storage, p%pcf%cpool_to_deadcrootc_storage, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_leaf_gr, p%pcf%cpool_leaf_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_froot_gr, p%pcf%cpool_froot_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_livestem_gr, p%pcf%cpool_livestem_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_deadstem_gr, p%pcf%cpool_deadstem_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_livecroot_gr, p%pcf%cpool_livecroot_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_deadcroot_gr, p%pcf%cpool_deadcroot_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_leaf_storage_gr, p%pcf%cpool_leaf_storage_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_froot_storage_gr, p%pcf%cpool_froot_storage_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_livestem_storage_gr, p%pcf%cpool_livestem_storage_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_deadstem_storage_gr, p%pcf%cpool_deadstem_storage_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_livecroot_storage_gr, p%pcf%cpool_livecroot_storage_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_deadcroot_storage_gr, p%pcf%cpool_deadcroot_storage_gr, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%cpool_to_gresp_storage, p%pcf%cpool_to_gresp_storage, &
                    p%pc13s%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%transfer_leaf_gr, p%pcf%transfer_leaf_gr, &
                    p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%transfer_froot_gr, p%pcf%transfer_froot_gr, &
                    p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%transfer_livestem_gr, p%pcf%transfer_livestem_gr, &
                    p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%transfer_deadstem_gr, p%pcf%transfer_deadstem_gr, &
                    p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%transfer_livecroot_gr, p%pcf%transfer_livecroot_gr, &
                    p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%transfer_deadcroot_gr, p%pcf%transfer_deadcroot_gr, &
                    p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%leafc_storage_to_xfer, p%pcf%leafc_storage_to_xfer, &
                    p%pc13s%leafc_storage, p%pcs%leafc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%frootc_storage_to_xfer, p%pcf%frootc_storage_to_xfer, &
                    p%pc13s%frootc_storage, p%pcs%frootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livestemc_storage_to_xfer, p%pcf%livestemc_storage_to_xfer, &
                    p%pc13s%livestemc_storage, p%pcs%livestemc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%deadstemc_storage_to_xfer, p%pcf%deadstemc_storage_to_xfer, &
                    p%pc13s%deadstemc_storage, p%pcs%deadstemc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%livecrootc_storage_to_xfer, p%pcf%livecrootc_storage_to_xfer, &
                    p%pc13s%livecrootc_storage, p%pcs%livecrootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%deadcrootc_storage_to_xfer, p%pcf%deadcrootc_storage_to_xfer, &
                    p%pc13s%deadcrootc_storage, p%pcs%deadcrootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(p%pc13f%gresp_storage_to_xfer, p%pcf%gresp_storage_to_xfer, &
                    p%pc13s%gresp_storage, p%pcs%gresp_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
!	call C13FluxCalc(p%pc13f%fx, p%pcf%fx, &
!                    p%pc13s%sx, p%pcs%sx, &
!                    num_soilp, filter_soilp, 1._r8, 0)
                    
end subroutine C13Flux1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13FluxCalc
!
! !INTERFACE:
subroutine C13FluxCalc(c13_flux, ctot_flux, c13_state, ctot_state, &
	                    num, filter, frax, diag)
!
! !DESCRIPTION:
! On the radiation time step, set the 13-carbon flux
! variables (except for gap-phase mortality and fire fluxes)
!
! !USES:
   use clmtype
   use time_manager, only: get_step_size
   use clm_varctl, only: irad
!
! !ARGUMENTS:
   implicit none
   real(r8), pointer   :: c13_flux(:)   !OUTPUT 13C flux
   real(r8), pointer   :: ctot_flux(:)  !INPUT  totC flux
   real(r8), pointer   :: c13_state(:)  !INPUT  13C state, upstream pool
   real(r8), pointer   :: ctot_state(:) !INPUT  totC state, upstream pool
   real(r8), intent(in):: frax          !fractionation factor (1 = no fractionation)
   integer, intent(in) :: num           ! number of filter members
   integer, intent(in) :: filter(:)     ! filter indices
   integer, intent(in) :: diag          !0=no diagnostics, 1=print diagnostics
!
! !CALLED FROM:
! subroutine C13Flux1
!
! !REVISION HISTORY:
!
! !OTHER LOCAL VARIABLES:
   integer :: i,f     ! indices
   real(r8) :: temp
!
   ! loop over the supplied filter
!dir$ concurrent
!cdir nodep
   do f = 1,num
      i = filter(f)
      if (ctot_state(i) /= 0._r8) then
      	c13_flux(i) = ctot_flux(i) * (c13_state(i)/ctot_state(i)) * frax
      else
      	c13_flux(i) = 0._r8
      end if
      
      if (diag == 1) then
      	if (ctot_state(i) /= 0._r8) then
         	temp = ((c13_state(i)/(ctot_state(i)*c13ratio))-1._r8)*1000._r8
         else
         	temp = 0._r8
         end if
      	write(6,*) f,i,c13_flux(i), ctot_flux(i), c13_state(i), ctot_state(i), temp
      end if
   end do
end subroutine C13FluxCalc
!-----------------------------------------------------------------------


end module CNC13FluxMod
 
