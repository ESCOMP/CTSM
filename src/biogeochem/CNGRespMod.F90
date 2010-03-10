#include <misc.h>
#include <preproc.h>

module CNGRespMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNGRespMod
!
! !DESCRIPTION:
! Module for growth respiration fluxes,
! for coupled carbon-nitrogen code.
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public :: CNGResp
!
! !REVISION HISTORY:
! 9/12/03: Created by Peter Thornton
! 10/27/03, Peter Thornton: migrated to vector data structures
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNGResp
!
! !INTERFACE:
subroutine CNGResp(num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic carbon state
! variables
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn, in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 8/1/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)         ! pft vegetation type
   real(r8), pointer :: cpool_to_leafc(:)
   real(r8), pointer :: cpool_to_leafc_storage(:)
   real(r8), pointer :: cpool_to_frootc(:)
   real(r8), pointer :: cpool_to_frootc_storage(:)
   real(r8), pointer :: cpool_to_livestemc(:)
   real(r8), pointer :: cpool_to_livestemc_storage(:)
   real(r8), pointer :: cpool_to_deadstemc(:)
   real(r8), pointer :: cpool_to_deadstemc_storage(:)
   real(r8), pointer :: cpool_to_livecrootc(:)
   real(r8), pointer :: cpool_to_livecrootc_storage(:)
   real(r8), pointer :: cpool_to_deadcrootc(:)
   real(r8), pointer :: cpool_to_deadcrootc_storage(:)
   real(r8), pointer :: leafc_xfer_to_leafc(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: woody(:) !binary flag for woody lifeform (1=woody, 0=not woody)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: cpool_leaf_gr(:)
   real(r8), pointer :: cpool_leaf_storage_gr(:)
   real(r8), pointer :: transfer_leaf_gr(:)
   real(r8), pointer :: cpool_froot_gr(:)
   real(r8), pointer :: cpool_froot_storage_gr(:)
   real(r8), pointer :: transfer_froot_gr(:)
   real(r8), pointer :: cpool_livestem_gr(:)
   real(r8), pointer :: cpool_livestem_storage_gr(:)
   real(r8), pointer :: transfer_livestem_gr(:)
   real(r8), pointer :: cpool_deadstem_gr(:)
   real(r8), pointer :: cpool_deadstem_storage_gr(:)
   real(r8), pointer :: transfer_deadstem_gr(:)
   real(r8), pointer :: cpool_livecroot_gr(:)
   real(r8), pointer :: cpool_livecroot_storage_gr(:)
   real(r8), pointer :: transfer_livecroot_gr(:)
   real(r8), pointer :: cpool_deadcroot_gr(:)
   real(r8), pointer :: cpool_deadcroot_storage_gr(:)
   real(r8), pointer :: transfer_deadcroot_gr(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p                ! indices
   integer :: fp               ! lake filter pft index
   real(r8):: grperc, grpnow   ! growth respirarion parameters

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   ivt                           => clm3%g%l%c%p%itype
   cpool_to_leafc                => clm3%g%l%c%p%pcf%cpool_to_leafc
   cpool_to_leafc_storage        => clm3%g%l%c%p%pcf%cpool_to_leafc_storage
   cpool_to_frootc               => clm3%g%l%c%p%pcf%cpool_to_frootc
   cpool_to_frootc_storage       => clm3%g%l%c%p%pcf%cpool_to_frootc_storage
   cpool_to_livestemc            => clm3%g%l%c%p%pcf%cpool_to_livestemc
   cpool_to_livestemc_storage    => clm3%g%l%c%p%pcf%cpool_to_livestemc_storage
   cpool_to_deadstemc            => clm3%g%l%c%p%pcf%cpool_to_deadstemc
   cpool_to_deadstemc_storage    => clm3%g%l%c%p%pcf%cpool_to_deadstemc_storage
   cpool_to_livecrootc           => clm3%g%l%c%p%pcf%cpool_to_livecrootc
   cpool_to_livecrootc_storage   => clm3%g%l%c%p%pcf%cpool_to_livecrootc_storage
   cpool_to_deadcrootc           => clm3%g%l%c%p%pcf%cpool_to_deadcrootc
   cpool_to_deadcrootc_storage   => clm3%g%l%c%p%pcf%cpool_to_deadcrootc_storage
   leafc_xfer_to_leafc           => clm3%g%l%c%p%pcf%leafc_xfer_to_leafc
   frootc_xfer_to_frootc         => clm3%g%l%c%p%pcf%frootc_xfer_to_frootc
   livestemc_xfer_to_livestemc   => clm3%g%l%c%p%pcf%livestemc_xfer_to_livestemc
   deadstemc_xfer_to_deadstemc   => clm3%g%l%c%p%pcf%deadstemc_xfer_to_deadstemc
   livecrootc_xfer_to_livecrootc => clm3%g%l%c%p%pcf%livecrootc_xfer_to_livecrootc
   deadcrootc_xfer_to_deadcrootc => clm3%g%l%c%p%pcf%deadcrootc_xfer_to_deadcrootc
   woody => pftcon%woody

   ! Assign local pointers to derived type arrays (out)
   cpool_leaf_gr                 => clm3%g%l%c%p%pcf%cpool_leaf_gr
   cpool_leaf_storage_gr         => clm3%g%l%c%p%pcf%cpool_leaf_storage_gr
   transfer_leaf_gr              => clm3%g%l%c%p%pcf%transfer_leaf_gr
   cpool_froot_gr                => clm3%g%l%c%p%pcf%cpool_froot_gr
   cpool_froot_storage_gr        => clm3%g%l%c%p%pcf%cpool_froot_storage_gr
   transfer_froot_gr             => clm3%g%l%c%p%pcf%transfer_froot_gr
   cpool_livestem_gr             => clm3%g%l%c%p%pcf%cpool_livestem_gr
   cpool_livestem_storage_gr     => clm3%g%l%c%p%pcf%cpool_livestem_storage_gr
   transfer_livestem_gr          => clm3%g%l%c%p%pcf%transfer_livestem_gr
   cpool_deadstem_gr             => clm3%g%l%c%p%pcf%cpool_deadstem_gr
   cpool_deadstem_storage_gr     => clm3%g%l%c%p%pcf%cpool_deadstem_storage_gr
   transfer_deadstem_gr          => clm3%g%l%c%p%pcf%transfer_deadstem_gr
   cpool_livecroot_gr            => clm3%g%l%c%p%pcf%cpool_livecroot_gr
   cpool_livecroot_storage_gr    => clm3%g%l%c%p%pcf%cpool_livecroot_storage_gr
   transfer_livecroot_gr         => clm3%g%l%c%p%pcf%transfer_livecroot_gr
   cpool_deadcroot_gr            => clm3%g%l%c%p%pcf%cpool_deadcroot_gr
   cpool_deadcroot_storage_gr    => clm3%g%l%c%p%pcf%cpool_deadcroot_storage_gr
   transfer_deadcroot_gr         => clm3%g%l%c%p%pcf%transfer_deadcroot_gr

   ! set some parameters (temporary, these will eventually go into
   ! either pepc, or parameter file
   grperc = 0.3_r8
   grpnow = 1.0_r8

   ! Loop through pfts
   ! start pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)


      ! leaf and fine root growth respiration
      cpool_leaf_gr(p)          = cpool_to_leafc(p) * grperc
      cpool_leaf_storage_gr(p)  = cpool_to_leafc_storage(p) * grperc * grpnow
      transfer_leaf_gr(p)       = leafc_xfer_to_leafc(p) * grperc * (1._r8 - grpnow)
      cpool_froot_gr(p)         = cpool_to_frootc(p) * grperc
      cpool_froot_storage_gr(p) = cpool_to_frootc_storage(p) * grperc * grpnow
      transfer_froot_gr(p)      = frootc_xfer_to_frootc(p) * grperc * (1._r8 - grpnow)

      if (woody(ivt(p)) == 1._r8) then
          cpool_livestem_gr(p)          = cpool_to_livestemc(p) * grperc
          cpool_livestem_storage_gr(p)  = cpool_to_livestemc_storage(p) * grperc * grpnow
          transfer_livestem_gr(p)       = livestemc_xfer_to_livestemc(p) * grperc * (1._r8 - grpnow)
          cpool_deadstem_gr(p)          = cpool_to_deadstemc(p) * grperc
          cpool_deadstem_storage_gr(p)  = cpool_to_deadstemc_storage(p) * grperc * grpnow
          transfer_deadstem_gr(p)       = deadstemc_xfer_to_deadstemc(p) * grperc * (1._r8 - grpnow)
          cpool_livecroot_gr(p)         = cpool_to_livecrootc(p) * grperc
          cpool_livecroot_storage_gr(p) = cpool_to_livecrootc_storage(p) * grperc * grpnow
          transfer_livecroot_gr(p)      = livecrootc_xfer_to_livecrootc(p) * grperc * (1._r8 - grpnow)
          cpool_deadcroot_gr(p)         = cpool_to_deadcrootc(p) * grperc
          cpool_deadcroot_storage_gr(p) = cpool_to_deadcrootc_storage(p) * grperc * grpnow
          transfer_deadcroot_gr(p)      = deadcrootc_xfer_to_deadcrootc(p) * grperc * (1._r8 - grpnow)
      end if

   end do

end subroutine CNGResp

#endif

end module CNGRespMod
