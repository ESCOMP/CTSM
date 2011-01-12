
module CNC13FluxMod
#if (defined CN) && (defined C13)

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
    implicit none
    save
    private
!
! !PUBLIC MEMBER FUNCTIONS:
    public:: C13Flux1
    public:: C13Flux2
    public:: C13Flux2h
    public:: C13Flux3
    private:: CNC13LitterToColumn
    private:: CNC13GapPftToColumn
    private:: CNC13HarvestPftToColumn
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
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(column_type), pointer :: c
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p
   c => clm3%g%l%c
	
   ! pft-level non-mortality fluxes
   
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
                    
	! call routine to shift pft-level litterfall fluxes to column, for isotopes
   ! the non-isotope version of this routine is called in CNPhenologyMod.F90
   ! For later clean-up, it would be possible to generalize this function to operate on a single 
   ! pft-to-column flux.
   
   call CNC13LitterToColumn(num_soilc, filter_soilc)
   
   ! column-level non-mortality fluxes
   
   call C13FluxCalc(c%cc13f%cwdc_to_litr2c, c%ccf%cwdc_to_litr2c, &
                     c%cc13s%cwdc, c%ccs%cwdc, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%cwdc_to_litr3c, c%ccf%cwdc_to_litr3c, &
                     c%cc13s%cwdc, c%ccs%cwdc, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%litr1_hr, c%ccf%litr1_hr, &
                     c%cc13s%litr1c, c%ccs%litr1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%litr1c_to_soil1c, c%ccf%litr1c_to_soil1c, &
                     c%cc13s%litr1c, c%ccs%litr1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%litr2_hr, c%ccf%litr2_hr, &
                     c%cc13s%litr2c, c%ccs%litr2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%litr2c_to_soil2c, c%ccf%litr2c_to_soil2c, &
                     c%cc13s%litr2c, c%ccs%litr2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%litr3_hr, c%ccf%litr3_hr, &
                     c%cc13s%litr3c, c%ccs%litr3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%litr3c_to_soil3c, c%ccf%litr3c_to_soil3c, &
                     c%cc13s%litr3c, c%ccs%litr3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%soil1_hr, c%ccf%soil1_hr, &
                     c%cc13s%soil1c, c%ccs%soil1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%soil1c_to_soil2c, c%ccf%soil1c_to_soil2c, &
                     c%cc13s%soil1c, c%ccs%soil1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%soil2_hr, c%ccf%soil2_hr, &
                     c%cc13s%soil2c, c%ccs%soil2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%soil2c_to_soil3c, c%ccf%soil2c_to_soil3c, &
                     c%cc13s%soil2c, c%ccs%soil2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%soil3_hr, c%ccf%soil3_hr, &
                     c%cc13s%soil3c, c%ccs%soil3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%soil3c_to_soil4c, c%ccf%soil3c_to_soil4c, &
                     c%cc13s%soil3c, c%ccs%soil3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(c%cc13f%soil4_hr, c%ccf%soil4_hr, &
                     c%cc13s%soil4c, c%ccs%soil4c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   
!	call C13FluxCalc(p%pc13f%fx, p%pcf%fx, &
!                    p%pc13s%sx, p%pcs%sx, &
!                    num_soilp, filter_soilp, 1._r8, 0)
                    
end subroutine C13Flux1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Flux2
!
! !INTERFACE:
subroutine C13Flux2(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, set the 13-carbon fluxes for gap mortality
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(column_type), pointer :: c
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p
   c => clm3%g%l%c
	
   ! pft-level gap mortality fluxes
   
   call C13FluxCalc(p%pc13f%m_leafc_to_litter, p%pcf%m_leafc_to_litter, &
                     p%pc13s%leafc, p%pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_leafc_storage_to_litter, p%pcf%m_leafc_storage_to_litter, &
                     p%pc13s%leafc_storage, p%pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_leafc_xfer_to_litter, p%pcf%m_leafc_xfer_to_litter, &
                     p%pc13s%leafc_xfer, p%pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_frootc_to_litter, p%pcf%m_frootc_to_litter, &
                     p%pc13s%frootc, p%pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_frootc_storage_to_litter, p%pcf%m_frootc_storage_to_litter, &
                     p%pc13s%frootc_storage, p%pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_frootc_xfer_to_litter, p%pcf%m_frootc_xfer_to_litter, &
                     p%pc13s%frootc_xfer, p%pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livestemc_to_litter, p%pcf%m_livestemc_to_litter, &
                     p%pc13s%livestemc, p%pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livestemc_storage_to_litter, p%pcf%m_livestemc_storage_to_litter, &
                     p%pc13s%livestemc_storage, p%pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livestemc_xfer_to_litter, p%pcf%m_livestemc_xfer_to_litter, &
                     p%pc13s%livestemc_xfer, p%pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadstemc_to_litter, p%pcf%m_deadstemc_to_litter, &
                     p%pc13s%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadstemc_storage_to_litter, p%pcf%m_deadstemc_storage_to_litter, &
                     p%pc13s%deadstemc_storage, p%pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadstemc_xfer_to_litter, p%pcf%m_deadstemc_xfer_to_litter, &
                     p%pc13s%deadstemc_xfer, p%pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livecrootc_to_litter, p%pcf%m_livecrootc_to_litter, &
                     p%pc13s%livecrootc, p%pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livecrootc_storage_to_litter, p%pcf%m_livecrootc_storage_to_litter, &
                     p%pc13s%livecrootc_storage, p%pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livecrootc_xfer_to_litter, p%pcf%m_livecrootc_xfer_to_litter, &
                     p%pc13s%livecrootc_xfer, p%pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadcrootc_to_litter, p%pcf%m_deadcrootc_to_litter, &
                     p%pc13s%deadcrootc, p%pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadcrootc_storage_to_litter, p%pcf%m_deadcrootc_storage_to_litter, &
                     p%pc13s%deadcrootc_storage, p%pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadcrootc_xfer_to_litter, p%pcf%m_deadcrootc_xfer_to_litter, &
                     p%pc13s%deadcrootc_xfer, p%pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_gresp_storage_to_litter, p%pcf%m_gresp_storage_to_litter, &
                     p%pc13s%gresp_storage, p%pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_gresp_xfer_to_litter, p%pcf%m_gresp_xfer_to_litter, &
                     p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
	! call routine to shift pft-level gap mortality fluxes to column, for isotopes
   ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

	call CNC13GapPftToColumn(num_soilc, filter_soilc)
   
end subroutine C13Flux2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Flux2h
!
! !INTERFACE:
subroutine C13Flux2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! set the 13-carbon fluxes for harvest mortality
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(column_type), pointer :: c
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p
   c => clm3%g%l%c
	
   ! pft-level gap mortality fluxes
   
   call C13FluxCalc(p%pc13f%hrv_leafc_to_litter, p%pcf%hrv_leafc_to_litter, &
                     p%pc13s%leafc, p%pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_leafc_storage_to_litter, p%pcf%hrv_leafc_storage_to_litter, &
                     p%pc13s%leafc_storage, p%pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_leafc_xfer_to_litter, p%pcf%hrv_leafc_xfer_to_litter, &
                     p%pc13s%leafc_xfer, p%pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_frootc_to_litter, p%pcf%hrv_frootc_to_litter, &
                     p%pc13s%frootc, p%pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_frootc_storage_to_litter, p%pcf%hrv_frootc_storage_to_litter, &
                     p%pc13s%frootc_storage, p%pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_frootc_xfer_to_litter, p%pcf%hrv_frootc_xfer_to_litter, &
                     p%pc13s%frootc_xfer, p%pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_livestemc_to_litter, p%pcf%hrv_livestemc_to_litter, &
                     p%pc13s%livestemc, p%pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_livestemc_storage_to_litter, p%pcf%hrv_livestemc_storage_to_litter, &
                     p%pc13s%livestemc_storage, p%pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_livestemc_xfer_to_litter, p%pcf%hrv_livestemc_xfer_to_litter, &
                     p%pc13s%livestemc_xfer, p%pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_deadstemc_to_prod10c, p%pcf%hrv_deadstemc_to_prod10c, &
                     p%pc13s%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_deadstemc_to_prod100c, p%pcf%hrv_deadstemc_to_prod100c, &
                     p%pc13s%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_deadstemc_storage_to_litter, p%pcf%hrv_deadstemc_storage_to_litter, &
                     p%pc13s%deadstemc_storage, p%pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_deadstemc_xfer_to_litter, p%pcf%hrv_deadstemc_xfer_to_litter, &
                     p%pc13s%deadstemc_xfer, p%pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_livecrootc_to_litter, p%pcf%hrv_livecrootc_to_litter, &
                     p%pc13s%livecrootc, p%pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_livecrootc_storage_to_litter, p%pcf%hrv_livecrootc_storage_to_litter, &
                     p%pc13s%livecrootc_storage, p%pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_livecrootc_xfer_to_litter, p%pcf%hrv_livecrootc_xfer_to_litter, &
                     p%pc13s%livecrootc_xfer, p%pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_deadcrootc_to_litter, p%pcf%hrv_deadcrootc_to_litter, &
                     p%pc13s%deadcrootc, p%pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_deadcrootc_storage_to_litter, p%pcf%hrv_deadcrootc_storage_to_litter, &
                     p%pc13s%deadcrootc_storage, p%pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_deadcrootc_xfer_to_litter, p%pcf%hrv_deadcrootc_xfer_to_litter, &
                     p%pc13s%deadcrootc_xfer, p%pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_gresp_storage_to_litter, p%pcf%hrv_gresp_storage_to_litter, &
                     p%pc13s%gresp_storage, p%pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%hrv_gresp_xfer_to_litter, p%pcf%hrv_gresp_xfer_to_litter, &
                     p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
	call C13FluxCalc(p%pc13f%hrv_xsmrpool_to_atm, p%pcf%hrv_xsmrpool_to_atm, &
                    p%pc13s%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	! call routine to shift pft-level gap mortality fluxes to column, for isotopes
   ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

	call CNC13HarvestPftToColumn(num_soilc, filter_soilc)
   
end subroutine C13Flux2h
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Flux3
!
! !INTERFACE:
subroutine C13Flux3(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, set the 13-carbon fluxes for fire mortality
!
! !USES:
   use clmtype
   use pft2colMod, only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(column_type), pointer :: c
   integer :: fp,pi
   real(r8), pointer :: ptrp(:)         ! pointer to input pft array
   real(r8), pointer :: ptrc(:)         ! pointer to output column array
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p
   c => clm3%g%l%c
	
   ! pft-level fire mortality fluxes
   
   call C13FluxCalc(p%pc13f%m_leafc_to_fire, p%pcf%m_leafc_to_fire, &
                     p%pc13s%leafc, p%pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_leafc_storage_to_fire, p%pcf%m_leafc_storage_to_fire, &
                     p%pc13s%leafc_storage, p%pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_leafc_xfer_to_fire, p%pcf%m_leafc_xfer_to_fire, &
                     p%pc13s%leafc_xfer, p%pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_frootc_to_fire, p%pcf%m_frootc_to_fire, &
                     p%pc13s%frootc, p%pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_frootc_storage_to_fire, p%pcf%m_frootc_storage_to_fire, &
                     p%pc13s%frootc_storage, p%pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_frootc_xfer_to_fire, p%pcf%m_frootc_xfer_to_fire, &
                     p%pc13s%frootc_xfer, p%pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livestemc_to_fire, p%pcf%m_livestemc_to_fire, &
                     p%pc13s%livestemc, p%pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livestemc_storage_to_fire, p%pcf%m_livestemc_storage_to_fire, &
                     p%pc13s%livestemc_storage, p%pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livestemc_xfer_to_fire, p%pcf%m_livestemc_xfer_to_fire, &
                     p%pc13s%livestemc_xfer, p%pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadstemc_to_fire, p%pcf%m_deadstemc_to_fire, &
                     p%pc13s%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadstemc_to_litter_fire, p%pcf%m_deadstemc_to_litter_fire, &
                     p%pc13s%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadstemc_storage_to_fire, p%pcf%m_deadstemc_storage_to_fire, &
                     p%pc13s%deadstemc_storage, p%pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadstemc_xfer_to_fire, p%pcf%m_deadstemc_xfer_to_fire, &
                     p%pc13s%deadstemc_xfer, p%pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livecrootc_to_fire, p%pcf%m_livecrootc_to_fire, &
                     p%pc13s%livecrootc, p%pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livecrootc_storage_to_fire, p%pcf%m_livecrootc_storage_to_fire, &
                     p%pc13s%livecrootc_storage, p%pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_livecrootc_xfer_to_fire, p%pcf%m_livecrootc_xfer_to_fire, &
                     p%pc13s%livecrootc_xfer, p%pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadcrootc_to_fire, p%pcf%m_deadcrootc_to_fire, &
                     p%pc13s%deadcrootc, p%pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadcrootc_to_litter_fire, p%pcf%m_deadcrootc_to_litter_fire, &
                     p%pc13s%deadcrootc, p%pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadcrootc_storage_to_fire, p%pcf%m_deadcrootc_storage_to_fire, &
                     p%pc13s%deadcrootc_storage, p%pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_deadcrootc_xfer_to_fire, p%pcf%m_deadcrootc_xfer_to_fire, &
                     p%pc13s%deadcrootc_xfer, p%pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_gresp_storage_to_fire, p%pcf%m_gresp_storage_to_fire, &
                     p%pc13s%gresp_storage, p%pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(p%pc13f%m_gresp_xfer_to_fire, p%pcf%m_gresp_xfer_to_fire, &
                     p%pc13s%gresp_xfer, p%pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
	! use routine p2c to calculate the column-level flux of deadstem and deadcrootc to
   ! cwdc as the result of fire mortality.
   call p2c(num_soilc, filter_soilc, p%pc13f%m_deadstemc_to_litter_fire, c%cc13f%m_deadstemc_to_cwdc_fire)
   call p2c(num_soilc, filter_soilc, p%pc13f%m_deadcrootc_to_litter_fire, c%cc13f%m_deadcrootc_to_cwdc_fire)

   call C13FluxCalc(c%cc13f%m_litr1c_to_fire, c%ccf%m_litr1c_to_fire, &
                     c%cc13s%litr1c, c%ccs%litr1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
                    
   call C13FluxCalc(c%cc13f%m_litr2c_to_fire, c%ccf%m_litr2c_to_fire, &
                     c%cc13s%litr2c, c%ccs%litr2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
                    
   call C13FluxCalc(c%cc13f%m_litr3c_to_fire, c%ccf%m_litr3c_to_fire, &
                     c%cc13s%litr3c, c%ccs%litr3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
                    
   call C13FluxCalc(c%cc13f%m_cwdc_to_fire, c%ccf%m_cwdc_to_fire, &
                     c%cc13s%cwdc, c%ccs%cwdc, &
                     num_soilc, filter_soilc, 1._r8, 0)
                    
!	call C13FluxCalc(p%pc13f%fx, p%pcf%fx, &
!                    p%pc13s%sx, p%pcs%sx, &
!                    num_soilc, filter_soilc, 1._r8, 0)
                    
end subroutine C13Flux3
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNC13LitterToColumn
!
! !INTERFACE:
subroutine CNC13LitterToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! called at the end of cn_phenology to gather all pft-level litterfall fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)          ! pft vegetation type
   real(r8), pointer :: wtcol(:)        ! weight (relative to column) for this pft (0-1)
   real(r8), pointer :: pwtgcell(:)     ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: lf_flab(:)      ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)      ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)      ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)      ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)      ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)      ! fine root litter lignin fraction
   integer , pointer :: npfts(:)        ! number of pfts for each column
   integer , pointer :: pfti(:)         ! beginning pft index for each column
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: leafc_to_litr1c(:)
   real(r8), pointer :: leafc_to_litr2c(:)
   real(r8), pointer :: leafc_to_litr3c(:)
   real(r8), pointer :: frootc_to_litr1c(:)
   real(r8), pointer :: frootc_to_litr2c(:)
   real(r8), pointer :: frootc_to_litr3c(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
    integer :: fc,c,pi,p
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    wtcol                          => clm3%g%l%c%p%wtcol
    pwtgcell                       => clm3%g%l%c%p%wtgcell  
    leafc_to_litter                => clm3%g%l%c%p%pc13f%leafc_to_litter
    frootc_to_litter               => clm3%g%l%c%p%pc13f%frootc_to_litter
    npfts                          => clm3%g%l%c%npfts
    pfti                           => clm3%g%l%c%pfti
    lf_flab                        => pftcon%lf_flab
    lf_fcel                        => pftcon%lf_fcel
    lf_flig                        => pftcon%lf_flig
    fr_flab                        => pftcon%fr_flab
    fr_fcel                        => pftcon%fr_fcel
    fr_flig                        => pftcon%fr_flig

   ! assign local pointers to derived type arrays (out)
    leafc_to_litr1c                => clm3%g%l%c%cc13f%leafc_to_litr1c
    leafc_to_litr2c                => clm3%g%l%c%cc13f%leafc_to_litr2c
    leafc_to_litr3c                => clm3%g%l%c%cc13f%leafc_to_litr3c
    frootc_to_litr1c               => clm3%g%l%c%cc13f%frootc_to_litr1c
    frootc_to_litr2c               => clm3%g%l%c%cc13f%frootc_to_litr2c
    frootc_to_litr3c               => clm3%g%l%c%cc13f%frootc_to_litr3c

   do pi = 1,max_pft_per_col
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if ( pi <=  npfts(c) ) then
            p = pfti(c) + pi - 1
            if (pwtgcell(p)>0._r8) then

               ! leaf litter carbon fluxes
               leafc_to_litr1c(c) = leafc_to_litr1c(c) + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               leafc_to_litr2c(c) = leafc_to_litr2c(c) + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               leafc_to_litr3c(c) = leafc_to_litr3c(c) + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter carbon fluxes
               frootc_to_litr1c(c) = frootc_to_litr1c(c) + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               frootc_to_litr2c(c) = frootc_to_litr2c(c) + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               frootc_to_litr3c(c) = frootc_to_litr3c(c) + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNC13LitterToColumn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNC13GapPftToColumn
!
! !INTERFACE:
subroutine CNC13GapPftToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! gather all pft-level gap mortality fluxes
! to the column level and assign them to the three litter pools (+ cwd pool)
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col, maxpatch_pft
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! soil column filter
!
! !CALLED FROM:
! subroutine CNphenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)      ! pft vegetation type
   real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
   real(r8), pointer :: pwtgcell(:) ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: lf_flab(:)  ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)  ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)  ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)  ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)  ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)  ! fine root litter lignin fraction
   integer , pointer :: npfts(:)    ! number of pfts for each column
   integer , pointer :: pfti(:)     ! beginning pft index for each column
   real(r8), pointer :: m_leafc_to_litter(:)
   real(r8), pointer :: m_frootc_to_litter(:)
   real(r8), pointer :: m_livestemc_to_litter(:)
   real(r8), pointer :: m_deadstemc_to_litter(:)
   real(r8), pointer :: m_livecrootc_to_litter(:)
   real(r8), pointer :: m_deadcrootc_to_litter(:)
   real(r8), pointer :: m_leafc_storage_to_litter(:)
   real(r8), pointer :: m_frootc_storage_to_litter(:)
   real(r8), pointer :: m_livestemc_storage_to_litter(:)
   real(r8), pointer :: m_deadstemc_storage_to_litter(:)
   real(r8), pointer :: m_livecrootc_storage_to_litter(:)
   real(r8), pointer :: m_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: m_gresp_storage_to_litter(:)
   real(r8), pointer :: m_leafc_xfer_to_litter(:)
   real(r8), pointer :: m_frootc_xfer_to_litter(:)
   real(r8), pointer :: m_livestemc_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: m_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: m_gresp_xfer_to_litter(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: m_leafc_to_litr1c(:)
   real(r8), pointer :: m_leafc_to_litr2c(:)
   real(r8), pointer :: m_leafc_to_litr3c(:)
   real(r8), pointer :: m_frootc_to_litr1c(:)
   real(r8), pointer :: m_frootc_to_litr2c(:)
   real(r8), pointer :: m_frootc_to_litr3c(:)
   real(r8), pointer :: m_livestemc_to_cwdc(:)
   real(r8), pointer :: m_deadstemc_to_cwdc(:)
   real(r8), pointer :: m_livecrootc_to_cwdc(:)
   real(r8), pointer :: m_deadcrootc_to_cwdc(:)
   real(r8), pointer :: m_leafc_storage_to_litr1c(:)
   real(r8), pointer :: m_frootc_storage_to_litr1c(:)
   real(r8), pointer :: m_livestemc_storage_to_litr1c(:)
   real(r8), pointer :: m_deadstemc_storage_to_litr1c(:)
   real(r8), pointer :: m_livecrootc_storage_to_litr1c(:)
   real(r8), pointer :: m_deadcrootc_storage_to_litr1c(:)
   real(r8), pointer :: m_gresp_storage_to_litr1c(:)
   real(r8), pointer :: m_leafc_xfer_to_litr1c(:)
   real(r8), pointer :: m_frootc_xfer_to_litr1c(:)
   real(r8), pointer :: m_livestemc_xfer_to_litr1c(:)
   real(r8), pointer :: m_deadstemc_xfer_to_litr1c(:)
   real(r8), pointer :: m_livecrootc_xfer_to_litr1c(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_litr1c(:)
   real(r8), pointer :: m_gresp_xfer_to_litr1c(:)
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p               ! indices
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti
   m_leafc_to_litr1c              => clm3%g%l%c%cc13f%m_leafc_to_litr1c
   m_leafc_to_litr2c              => clm3%g%l%c%cc13f%m_leafc_to_litr2c
   m_leafc_to_litr3c              => clm3%g%l%c%cc13f%m_leafc_to_litr3c
   m_frootc_to_litr1c             => clm3%g%l%c%cc13f%m_frootc_to_litr1c
   m_frootc_to_litr2c             => clm3%g%l%c%cc13f%m_frootc_to_litr2c
   m_frootc_to_litr3c             => clm3%g%l%c%cc13f%m_frootc_to_litr3c
   m_livestemc_to_cwdc            => clm3%g%l%c%cc13f%m_livestemc_to_cwdc
   m_deadstemc_to_cwdc            => clm3%g%l%c%cc13f%m_deadstemc_to_cwdc
   m_livecrootc_to_cwdc           => clm3%g%l%c%cc13f%m_livecrootc_to_cwdc
   m_deadcrootc_to_cwdc           => clm3%g%l%c%cc13f%m_deadcrootc_to_cwdc
   m_leafc_storage_to_litr1c      => clm3%g%l%c%cc13f%m_leafc_storage_to_litr1c
   m_frootc_storage_to_litr1c     => clm3%g%l%c%cc13f%m_frootc_storage_to_litr1c
   m_livestemc_storage_to_litr1c  => clm3%g%l%c%cc13f%m_livestemc_storage_to_litr1c
   m_deadstemc_storage_to_litr1c  => clm3%g%l%c%cc13f%m_deadstemc_storage_to_litr1c
   m_livecrootc_storage_to_litr1c => clm3%g%l%c%cc13f%m_livecrootc_storage_to_litr1c
   m_deadcrootc_storage_to_litr1c => clm3%g%l%c%cc13f%m_deadcrootc_storage_to_litr1c
   m_gresp_storage_to_litr1c      => clm3%g%l%c%cc13f%m_gresp_storage_to_litr1c
   m_leafc_xfer_to_litr1c         => clm3%g%l%c%cc13f%m_leafc_xfer_to_litr1c
   m_frootc_xfer_to_litr1c        => clm3%g%l%c%cc13f%m_frootc_xfer_to_litr1c
   m_livestemc_xfer_to_litr1c     => clm3%g%l%c%cc13f%m_livestemc_xfer_to_litr1c
   m_deadstemc_xfer_to_litr1c     => clm3%g%l%c%cc13f%m_deadstemc_xfer_to_litr1c
   m_livecrootc_xfer_to_litr1c    => clm3%g%l%c%cc13f%m_livecrootc_xfer_to_litr1c
   m_deadcrootc_xfer_to_litr1c    => clm3%g%l%c%cc13f%m_deadcrootc_xfer_to_litr1c
   m_gresp_xfer_to_litr1c         => clm3%g%l%c%cc13f%m_gresp_xfer_to_litr1c

   ! assign local pointers to pft-level arrays
   ivt                            => clm3%g%l%c%p%itype
   wtcol                          => clm3%g%l%c%p%wtcol
   pwtgcell                       => clm3%g%l%c%p%wtgcell  
   m_leafc_to_litter              => clm3%g%l%c%p%pc13f%m_leafc_to_litter
   m_frootc_to_litter             => clm3%g%l%c%p%pc13f%m_frootc_to_litter
   m_livestemc_to_litter          => clm3%g%l%c%p%pc13f%m_livestemc_to_litter
   m_deadstemc_to_litter          => clm3%g%l%c%p%pc13f%m_deadstemc_to_litter
   m_livecrootc_to_litter         => clm3%g%l%c%p%pc13f%m_livecrootc_to_litter
   m_deadcrootc_to_litter         => clm3%g%l%c%p%pc13f%m_deadcrootc_to_litter
   m_leafc_storage_to_litter      => clm3%g%l%c%p%pc13f%m_leafc_storage_to_litter
   m_frootc_storage_to_litter     => clm3%g%l%c%p%pc13f%m_frootc_storage_to_litter
   m_livestemc_storage_to_litter  => clm3%g%l%c%p%pc13f%m_livestemc_storage_to_litter
   m_deadstemc_storage_to_litter  => clm3%g%l%c%p%pc13f%m_deadstemc_storage_to_litter
   m_livecrootc_storage_to_litter => clm3%g%l%c%p%pc13f%m_livecrootc_storage_to_litter
   m_deadcrootc_storage_to_litter => clm3%g%l%c%p%pc13f%m_deadcrootc_storage_to_litter
   m_gresp_storage_to_litter      => clm3%g%l%c%p%pc13f%m_gresp_storage_to_litter
   m_leafc_xfer_to_litter         => clm3%g%l%c%p%pc13f%m_leafc_xfer_to_litter
   m_frootc_xfer_to_litter        => clm3%g%l%c%p%pc13f%m_frootc_xfer_to_litter
   m_livestemc_xfer_to_litter     => clm3%g%l%c%p%pc13f%m_livestemc_xfer_to_litter
   m_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pc13f%m_deadstemc_xfer_to_litter
   m_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pc13f%m_livecrootc_xfer_to_litter
   m_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pc13f%m_deadcrootc_xfer_to_litter
   m_gresp_xfer_to_litter         => clm3%g%l%c%p%pc13f%m_gresp_xfer_to_litter

   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1

            if (pwtgcell(p)>0._r8) then

               ! leaf gap mortality carbon fluxes
               m_leafc_to_litr1c(c) = m_leafc_to_litr1c(c) + &
                  m_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               m_leafc_to_litr2c(c) = m_leafc_to_litr2c(c) + &
                  m_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               m_leafc_to_litr3c(c) = m_leafc_to_litr3c(c) + &
                  m_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root gap mortality carbon fluxes
               m_frootc_to_litr1c(c) = m_frootc_to_litr1c(c) + &
                  m_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               m_frootc_to_litr2c(c) = m_frootc_to_litr2c(c) + &
                  m_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               m_frootc_to_litr3c(c) = m_frootc_to_litr3c(c) + &
                  m_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood gap mortality carbon fluxes
               m_livestemc_to_cwdc(c)  = m_livestemc_to_cwdc(c)  + &
                  m_livestemc_to_litter(p)  * wtcol(p)
               m_deadstemc_to_cwdc(c)  = m_deadstemc_to_cwdc(c)  + &
                  m_deadstemc_to_litter(p)  * wtcol(p)
               m_livecrootc_to_cwdc(c) = m_livecrootc_to_cwdc(c) + &
                  m_livecrootc_to_litter(p) * wtcol(p)
               m_deadcrootc_to_cwdc(c) = m_deadcrootc_to_cwdc(c) + &
                  m_deadcrootc_to_litter(p) * wtcol(p)

               ! storage gap mortality carbon fluxes
               m_leafc_storage_to_litr1c(c)      = m_leafc_storage_to_litr1c(c)      + &
                  m_leafc_storage_to_litter(p)      * wtcol(p)
               m_frootc_storage_to_litr1c(c)     = m_frootc_storage_to_litr1c(c)     + &
                  m_frootc_storage_to_litter(p)     * wtcol(p)
               m_livestemc_storage_to_litr1c(c)  = m_livestemc_storage_to_litr1c(c)  + &
                  m_livestemc_storage_to_litter(p)  * wtcol(p)
               m_deadstemc_storage_to_litr1c(c)  = m_deadstemc_storage_to_litr1c(c)  + &
                  m_deadstemc_storage_to_litter(p)  * wtcol(p)
               m_livecrootc_storage_to_litr1c(c) = m_livecrootc_storage_to_litr1c(c) + &
                  m_livecrootc_storage_to_litter(p) * wtcol(p)
               m_deadcrootc_storage_to_litr1c(c) = m_deadcrootc_storage_to_litr1c(c) + &
                  m_deadcrootc_storage_to_litter(p) * wtcol(p)
               m_gresp_storage_to_litr1c(c)      = m_gresp_storage_to_litr1c(c)      + &
                  m_gresp_storage_to_litter(p)      * wtcol(p)

               ! transfer gap mortality carbon fluxes
               m_leafc_xfer_to_litr1c(c)      = m_leafc_xfer_to_litr1c(c)      + &
                  m_leafc_xfer_to_litter(p)      * wtcol(p)
               m_frootc_xfer_to_litr1c(c)     = m_frootc_xfer_to_litr1c(c)     + &
                  m_frootc_xfer_to_litter(p)     * wtcol(p)
               m_livestemc_xfer_to_litr1c(c)  = m_livestemc_xfer_to_litr1c(c)  + &
                  m_livestemc_xfer_to_litter(p)  * wtcol(p)
               m_deadstemc_xfer_to_litr1c(c)  = m_deadstemc_xfer_to_litr1c(c)  + &
                  m_deadstemc_xfer_to_litter(p)  * wtcol(p)
               m_livecrootc_xfer_to_litr1c(c) = m_livecrootc_xfer_to_litr1c(c) + &
                  m_livecrootc_xfer_to_litter(p) * wtcol(p)
               m_deadcrootc_xfer_to_litr1c(c) = m_deadcrootc_xfer_to_litr1c(c) + &
                  m_deadcrootc_xfer_to_litter(p) * wtcol(p)
               m_gresp_xfer_to_litr1c(c)      = m_gresp_xfer_to_litr1c(c)      + &
                  m_gresp_xfer_to_litter(p)      * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNC13GapPftToColumn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNC13HarvestPftToColumn
!
! !INTERFACE:
subroutine CNC13HarvestPftToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! gather all pft-level harvest mortality fluxes
! to the column level and assign them to the litter, cwd, and wood product pools
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col, maxpatch_pft
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! soil column filter
!
! !CALLED FROM:
! subroutine CNphenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)      ! pft vegetation type
   real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
   real(r8), pointer :: pwtgcell(:) ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: lf_flab(:)  ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)  ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)  ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)  ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)  ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)  ! fine root litter lignin fraction
   integer , pointer :: npfts(:)    ! number of pfts for each column
   integer , pointer :: pfti(:)     ! beginning pft index for each column
   real(r8), pointer :: hrv_leafc_to_litter(:)
   real(r8), pointer :: hrv_frootc_to_litter(:)
   real(r8), pointer :: hrv_livestemc_to_litter(:)
   real(r8), pointer :: phrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: phrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_to_litter(:)
   real(r8), pointer :: hrv_leafc_storage_to_litter(:)
   real(r8), pointer :: hrv_frootc_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_gresp_storage_to_litter(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litter(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: hrv_leafc_to_litr1c(:)
   real(r8), pointer :: hrv_leafc_to_litr2c(:)
   real(r8), pointer :: hrv_leafc_to_litr3c(:)
   real(r8), pointer :: hrv_frootc_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_to_litr2c(:)
   real(r8), pointer :: hrv_frootc_to_litr3c(:)
   real(r8), pointer :: hrv_livestemc_to_cwdc(:)
   real(r8), pointer :: chrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: chrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_cwdc(:)
   real(r8), pointer :: hrv_deadcrootc_to_cwdc(:)
   real(r8), pointer :: hrv_leafc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_gresp_storage_to_litr1c(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litr1c(:)
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p               ! indices
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti
   hrv_leafc_to_litr1c              => clm3%g%l%c%cc13f%hrv_leafc_to_litr1c
   hrv_leafc_to_litr2c              => clm3%g%l%c%cc13f%hrv_leafc_to_litr2c
   hrv_leafc_to_litr3c              => clm3%g%l%c%cc13f%hrv_leafc_to_litr3c
   hrv_frootc_to_litr1c             => clm3%g%l%c%cc13f%hrv_frootc_to_litr1c
   hrv_frootc_to_litr2c             => clm3%g%l%c%cc13f%hrv_frootc_to_litr2c
   hrv_frootc_to_litr3c             => clm3%g%l%c%cc13f%hrv_frootc_to_litr3c
   hrv_livestemc_to_cwdc            => clm3%g%l%c%cc13f%hrv_livestemc_to_cwdc
   chrv_deadstemc_to_prod10c        => clm3%g%l%c%cc13f%hrv_deadstemc_to_prod10c
   chrv_deadstemc_to_prod100c       => clm3%g%l%c%cc13f%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_cwdc           => clm3%g%l%c%cc13f%hrv_livecrootc_to_cwdc
   hrv_deadcrootc_to_cwdc           => clm3%g%l%c%cc13f%hrv_deadcrootc_to_cwdc
   hrv_leafc_storage_to_litr1c      => clm3%g%l%c%cc13f%hrv_leafc_storage_to_litr1c
   hrv_frootc_storage_to_litr1c     => clm3%g%l%c%cc13f%hrv_frootc_storage_to_litr1c
   hrv_livestemc_storage_to_litr1c  => clm3%g%l%c%cc13f%hrv_livestemc_storage_to_litr1c
   hrv_deadstemc_storage_to_litr1c  => clm3%g%l%c%cc13f%hrv_deadstemc_storage_to_litr1c
   hrv_livecrootc_storage_to_litr1c => clm3%g%l%c%cc13f%hrv_livecrootc_storage_to_litr1c
   hrv_deadcrootc_storage_to_litr1c => clm3%g%l%c%cc13f%hrv_deadcrootc_storage_to_litr1c
   hrv_gresp_storage_to_litr1c      => clm3%g%l%c%cc13f%hrv_gresp_storage_to_litr1c
   hrv_leafc_xfer_to_litr1c         => clm3%g%l%c%cc13f%hrv_leafc_xfer_to_litr1c
   hrv_frootc_xfer_to_litr1c        => clm3%g%l%c%cc13f%hrv_frootc_xfer_to_litr1c
   hrv_livestemc_xfer_to_litr1c     => clm3%g%l%c%cc13f%hrv_livestemc_xfer_to_litr1c
   hrv_deadstemc_xfer_to_litr1c     => clm3%g%l%c%cc13f%hrv_deadstemc_xfer_to_litr1c
   hrv_livecrootc_xfer_to_litr1c    => clm3%g%l%c%cc13f%hrv_livecrootc_xfer_to_litr1c
   hrv_deadcrootc_xfer_to_litr1c    => clm3%g%l%c%cc13f%hrv_deadcrootc_xfer_to_litr1c
   hrv_gresp_xfer_to_litr1c         => clm3%g%l%c%cc13f%hrv_gresp_xfer_to_litr1c

   ! assign local pointers to pft-level arrays
   ivt                            => clm3%g%l%c%p%itype
   wtcol                          => clm3%g%l%c%p%wtcol
   pwtgcell                       => clm3%g%l%c%p%wtgcell  
   hrv_leafc_to_litter              => clm3%g%l%c%p%pc13f%hrv_leafc_to_litter
   hrv_frootc_to_litter             => clm3%g%l%c%p%pc13f%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => clm3%g%l%c%p%pc13f%hrv_livestemc_to_litter
   phrv_deadstemc_to_prod10c        => clm3%g%l%c%p%pc13f%hrv_deadstemc_to_prod10c
   phrv_deadstemc_to_prod100c       => clm3%g%l%c%p%pc13f%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => clm3%g%l%c%p%pc13f%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => clm3%g%l%c%p%pc13f%hrv_deadcrootc_to_litter
   hrv_leafc_storage_to_litter      => clm3%g%l%c%p%pc13f%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => clm3%g%l%c%p%pc13f%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => clm3%g%l%c%p%pc13f%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => clm3%g%l%c%p%pc13f%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => clm3%g%l%c%p%pc13f%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => clm3%g%l%c%p%pc13f%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => clm3%g%l%c%p%pc13f%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => clm3%g%l%c%p%pc13f%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => clm3%g%l%c%p%pc13f%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => clm3%g%l%c%p%pc13f%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pc13f%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pc13f%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pc13f%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => clm3%g%l%c%p%pc13f%hrv_gresp_xfer_to_litter

   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1

            if (pwtgcell(p)>0._r8) then

               ! leaf harvest mortality carbon fluxes
               hrv_leafc_to_litr1c(c) = hrv_leafc_to_litr1c(c) + &
                  hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr2c(c) = hrv_leafc_to_litr2c(c) + &
                  hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr3c(c) = hrv_leafc_to_litr3c(c) + &
                  hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root harvest mortality carbon fluxes
               hrv_frootc_to_litr1c(c) = hrv_frootc_to_litr1c(c) + &
                  hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr2c(c) = hrv_frootc_to_litr2c(c) + &
                  hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr3c(c) = hrv_frootc_to_litr3c(c) + &
                  hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood harvest mortality carbon fluxes
               hrv_livestemc_to_cwdc(c)  = hrv_livestemc_to_cwdc(c)  + &
                  hrv_livestemc_to_litter(p)  * wtcol(p)
               chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                  phrv_deadstemc_to_prod10c(p)  * wtcol(p)
               chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                  phrv_deadstemc_to_prod100c(p)  * wtcol(p)
               hrv_livecrootc_to_cwdc(c) = hrv_livecrootc_to_cwdc(c) + &
                  hrv_livecrootc_to_litter(p) * wtcol(p)
               hrv_deadcrootc_to_cwdc(c) = hrv_deadcrootc_to_cwdc(c) + &
                  hrv_deadcrootc_to_litter(p) * wtcol(p)

               ! storage harvest mortality carbon fluxes
               hrv_leafc_storage_to_litr1c(c)      = hrv_leafc_storage_to_litr1c(c)      + &
                  hrv_leafc_storage_to_litter(p)      * wtcol(p)
               hrv_frootc_storage_to_litr1c(c)     = hrv_frootc_storage_to_litr1c(c)     + &
                  hrv_frootc_storage_to_litter(p)     * wtcol(p)
               hrv_livestemc_storage_to_litr1c(c)  = hrv_livestemc_storage_to_litr1c(c)  + &
                  hrv_livestemc_storage_to_litter(p)  * wtcol(p)
               hrv_deadstemc_storage_to_litr1c(c)  = hrv_deadstemc_storage_to_litr1c(c)  + &
                  hrv_deadstemc_storage_to_litter(p)  * wtcol(p)
               hrv_livecrootc_storage_to_litr1c(c) = hrv_livecrootc_storage_to_litr1c(c) + &
                  hrv_livecrootc_storage_to_litter(p) * wtcol(p)
               hrv_deadcrootc_storage_to_litr1c(c) = hrv_deadcrootc_storage_to_litr1c(c) + &
                  hrv_deadcrootc_storage_to_litter(p) * wtcol(p)
               hrv_gresp_storage_to_litr1c(c)      = hrv_gresp_storage_to_litr1c(c)      + &
                  hrv_gresp_storage_to_litter(p)      * wtcol(p)

               ! transfer harvest mortality carbon fluxes
               hrv_leafc_xfer_to_litr1c(c)      = hrv_leafc_xfer_to_litr1c(c)      + &
                  hrv_leafc_xfer_to_litter(p)      * wtcol(p)
               hrv_frootc_xfer_to_litr1c(c)     = hrv_frootc_xfer_to_litr1c(c)     + &
                  hrv_frootc_xfer_to_litter(p)     * wtcol(p)
               hrv_livestemc_xfer_to_litr1c(c)  = hrv_livestemc_xfer_to_litr1c(c)  + &
                  hrv_livestemc_xfer_to_litter(p)  * wtcol(p)
               hrv_deadstemc_xfer_to_litr1c(c)  = hrv_deadstemc_xfer_to_litr1c(c)  + &
                  hrv_deadstemc_xfer_to_litter(p)  * wtcol(p)
               hrv_livecrootc_xfer_to_litr1c(c) = hrv_livecrootc_xfer_to_litr1c(c) + &
                  hrv_livecrootc_xfer_to_litter(p) * wtcol(p)
               hrv_deadcrootc_xfer_to_litr1c(c) = hrv_deadcrootc_xfer_to_litr1c(c) + &
                  hrv_deadcrootc_xfer_to_litter(p) * wtcol(p)
               hrv_gresp_xfer_to_litr1c(c)      = hrv_gresp_xfer_to_litr1c(c)      + &
                  hrv_gresp_xfer_to_litter(p)      * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNC13HarvestPftToColumn
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
   do f = 1,num
      i = filter(f)
      if (ctot_state(i) /= 0._r8) then
      	c13_flux(i) = ctot_flux(i) * (c13_state(i)/ctot_state(i)) * frax
      else
      	c13_flux(i) = 0._r8
      end if
      
      if (diag == 1) then
      ! put diagnostic print statements here for 13C flux calculations
      end if
   end do
end subroutine C13FluxCalc
!-----------------------------------------------------------------------

#endif

end module CNC13FluxMod
 
