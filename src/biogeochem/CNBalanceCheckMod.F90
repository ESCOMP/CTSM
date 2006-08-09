#include <misc.h>
#include <preproc.h>

module CNBalanceCheckMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNBalanceCheckMod
!
! !DESCRIPTION:
! Module for carbon mass balance checking.
!
! !USES:
    use abortutils  , only: endrun
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varcon  , only: istsoil
    use spmdMod     , only: masterproc
    use clm_varpar  , only: nlevsoi
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: BeginCBalance
    public :: BeginNBalance
    public :: CBalanceCheck
    public :: NBalanceCheck
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
! !IROUTINE: BeginCBalance
!
! !INTERFACE:
subroutine BeginCBalance(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, calculate the beginning carbon balance for mass
! conservation checks.
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
! subroutine driver
!
! !REVISION HISTORY:
! 2/4/05: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real(r8), pointer :: cwdc(:)          ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)        ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)        ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)        ! (gC/m2) litter lignin C
   real(r8), pointer :: soil1c(:)        ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)        ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)        ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)        ! (gC/m2) soil organic matter C (slowest pool)
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: col_ctrunc(:)         ! (gC/m2) column-level sink for C truncation
   real(r8), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_begcb(:)   ! carbon mass, beginning of time step (gC/m**2)
   real(r8), pointer :: pft_begcb(:)   ! carbon mass, beginning of time step (gC/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p     ! indices
   integer :: fp,fc   ! lake filter indices
!
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers at the column level
   col_begcb                      => clm3%g%l%c%ccbal%begcb
   col_ctrunc                     => clm3%g%l%c%ccs%col_ctrunc
   cwdc                           => clm3%g%l%c%ccs%cwdc
   litr1c                         => clm3%g%l%c%ccs%litr1c
   litr2c                         => clm3%g%l%c%ccs%litr2c
   litr3c                         => clm3%g%l%c%ccs%litr3c
   soil1c                         => clm3%g%l%c%ccs%soil1c
   soil2c                         => clm3%g%l%c%ccs%soil2c
   soil3c                         => clm3%g%l%c%ccs%soil3c
   soil4c                         => clm3%g%l%c%ccs%soil4c

   ! assign local pointers at the pft level
   pft_begcb                      => clm3%g%l%c%p%pcbal%begcb
   cpool                          => clm3%g%l%c%p%pcs%cpool
   xsmrpool                       => clm3%g%l%c%p%pcs%xsmrpool
   deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
   deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
   deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
   deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
   deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
   deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
   frootc                         => clm3%g%l%c%p%pcs%frootc
   frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
   frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
   gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
   gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
   leafc                          => clm3%g%l%c%p%pcs%leafc
   leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
   leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
   livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
   livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
   livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
   livestemc                      => clm3%g%l%c%p%pcs%livestemc
   livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
   livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
   pft_ctrunc                     => clm3%g%l%c%p%pcs%pft_ctrunc

   ! column loop
!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c = filter_soilc(fc)
 
      ! calculate beginning column-level carbon balance,
      ! for mass conservation check
 
      col_begcb(c) = cwdc(c) + &
         litr1c(c) + litr2c(c) + litr3c(c) + &
         soil1c(c) + soil2c(c) + soil3c(c) + soil4c(c) + &
         col_ctrunc(c)

   end do ! end of columns loop
 
   ! pft loop
!dir$ concurrent
!cdir nodep
   do fp = 1,num_soilp
      p = filter_soilp(fp)
 
      ! calculate beginning pft-level carbon balance,
      ! for mass conservation check
 
      pft_begcb(p) = &
         leafc(p)      + leafc_storage(p)      + leafc_xfer(p) + &
         frootc(p)     + frootc_storage(p)     + frootc_xfer(p) + &
         livestemc(p)  + livestemc_storage(p)  + livestemc_xfer(p) + &
         deadstemc(p)  + deadstemc_storage(p)  + deadstemc_xfer(p) + &
         livecrootc(p) + livecrootc_storage(p) + livecrootc_xfer(p) + &
         deadcrootc(p) + deadcrootc_storage(p) + deadcrootc_xfer(p) + &
         gresp_storage(p) + gresp_xfer(p)  + cpool(p) +  xsmrpool(p) + &
         pft_ctrunc(p)
 
   end do ! end of pft loop

end subroutine BeginCBalance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BeginNBalance
!
! !INTERFACE:
subroutine BeginNBalance(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, calculate the beginning nitrogen balance for mass
! conservation checks.
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
! subroutine driver
!
! !REVISION HISTORY:
! 2/4/05: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real(r8), pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real(r8), pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real(r8), pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real(r8), pointer :: sminn(:)              ! (gN/m2) soil mineral N
   real(r8), pointer :: soil1n(:)             ! (gN/m2) soil organic matter N (fast pool)
   real(r8), pointer :: soil2n(:)             ! (gN/m2) soil organic matter N (medium pool)
   real(r8), pointer :: soil3n(:)             ! (gN/m2) soil orgainc matter N (slow pool)
   real(r8), pointer :: soil4n(:)             ! (gN/m2) soil orgainc matter N (slowest pool)
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real(r8), pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   real(r8), pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_begnb(:)   ! nitrogen mass, beginning of time step (gN/m**2)
   real(r8), pointer :: pft_begnb(:)   ! nitrogen mass, beginning of time step (gN/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p     ! indices
   integer :: fp,fc   ! lake filter indices
!
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers at the column level
   col_begnb                      => clm3%g%l%c%cnbal%begnb
   col_ntrunc                     => clm3%g%l%c%cns%col_ntrunc
   cwdn                           => clm3%g%l%c%cns%cwdn
   litr1n                         => clm3%g%l%c%cns%litr1n
   litr2n                         => clm3%g%l%c%cns%litr2n
   litr3n                         => clm3%g%l%c%cns%litr3n
   sminn                          => clm3%g%l%c%cns%sminn
   soil1n                         => clm3%g%l%c%cns%soil1n
   soil2n                         => clm3%g%l%c%cns%soil2n
   soil3n                         => clm3%g%l%c%cns%soil3n
   soil4n                         => clm3%g%l%c%cns%soil4n

   ! assign local pointers at the pft level
   pft_begnb                      => clm3%g%l%c%p%pnbal%begnb
   pft_ntrunc                     => clm3%g%l%c%p%pns%pft_ntrunc
   deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
   deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
   deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
   deadstemn                      => clm3%g%l%c%p%pns%deadstemn
   deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
   deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
   frootn                         => clm3%g%l%c%p%pns%frootn
   frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
   frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
   leafn                          => clm3%g%l%c%p%pns%leafn
   leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
   leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
   livecrootn                     => clm3%g%l%c%p%pns%livecrootn
   livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
   livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
   livestemn                      => clm3%g%l%c%p%pns%livestemn
   livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
   livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
   npool                          => clm3%g%l%c%p%pns%npool
   retransn                       => clm3%g%l%c%p%pns%retransn

   ! column loop
!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c = filter_soilc(fc)
 
      ! calculate beginning column-level nitrogen balance,
      ! for mass conservation check
 
      col_begnb(c) = cwdn(c) + &
        litr1n(c) + litr2n(c) + litr3n(c) + &
        soil1n(c) + soil2n(c) + soil3n(c) + soil4n(c) + &
        sminn(c)  + col_ntrunc(c)

   end do ! end of columns loop
 
   ! pft loop
!dir$ concurrent
!cdir nodep
   do fp = 1,num_soilp
      p = filter_soilp(fp)
 
      ! calculate beginning pft-level nitrogen balance,
      ! for mass conservation check
 
      pft_begnb(p) = &
         leafn(p)      + leafn_storage(p)      + leafn_xfer(p) + &
         frootn(p)     + frootn_storage(p)     + frootn_xfer(p) + &
         livestemn(p)  + livestemn_storage(p)  + livestemn_xfer(p) + &
         deadstemn(p)  + deadstemn_storage(p)  + deadstemn_xfer(p) + &
         livecrootn(p) + livecrootn_storage(p) + livecrootn_xfer(p) + &
         deadcrootn(p) + deadcrootn_storage(p) + deadcrootn_xfer(p) + &
         retransn(p)   + npool(p)              + pft_ntrunc(p)
 
   end do ! end of pft loop

end subroutine BeginNBalance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CBalanceCheck
!
! !INTERFACE:
subroutine CBalanceCheck(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, perform carbon mass conservation check for column and pft
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varctl, only: irad
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
! 12/9/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   real(r8), pointer :: cwdc(:)       ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)     ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)     ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)     ! (gC/m2) litter lignin C
   real(r8), pointer :: soil1c(:)     ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)     ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)     ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)     ! (gC/m2) soil organic matter C (slowest pool)
   real(r8), pointer :: col_ctrunc(:) ! (gC/m2) column-level sink for C truncation
   real(r8), pointer :: leafc_to_litr1c(:)
   real(r8), pointer :: leafc_to_litr2c(:)
   real(r8), pointer :: leafc_to_litr3c(:)
   real(r8), pointer :: frootc_to_litr1c(:)
   real(r8), pointer :: frootc_to_litr2c(:)
   real(r8), pointer :: frootc_to_litr3c(:)
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
   real(r8), pointer :: m_deadstemc_to_cwdc_fire(:)
   real(r8), pointer :: m_deadcrootc_to_cwdc_fire(:)
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
   real(r8), pointer :: litr1_hr(:)
   real(r8), pointer :: litr2_hr(:)
   real(r8), pointer :: litr3_hr(:)
   real(r8), pointer :: soil1_hr(:)
   real(r8), pointer :: soil2_hr(:)
   real(r8), pointer :: soil3_hr(:)
   real(r8), pointer :: soil4_hr(:)
   real(r8), pointer :: m_litr1c_to_fire(:)
   real(r8), pointer :: m_litr2c_to_fire(:)
   real(r8), pointer :: m_litr3c_to_fire(:)
   real(r8), pointer :: m_cwdc_to_fire(:)
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   real(r8), pointer :: psnsun_to_cpool(:)
   real(r8), pointer :: psnshade_to_cpool(:)
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: m_leafc_to_litter(:)
   real(r8), pointer :: m_frootc_to_litter(:)
   real(r8), pointer :: m_livestemc_to_litter(:)
   real(r8), pointer :: m_deadstemc_to_litter(:)
   real(r8), pointer :: m_livecrootc_to_litter(:)
   real(r8), pointer :: m_deadcrootc_to_litter(:)
   real(r8), pointer :: m_leafc_to_fire(:)
   real(r8), pointer :: m_frootc_to_fire(:)
   real(r8), pointer :: m_livestemc_to_fire(:)
   real(r8), pointer :: m_deadstemc_to_fire(:)
   real(r8), pointer :: m_livecrootc_to_fire(:)
   real(r8), pointer :: m_deadcrootc_to_fire(:)
   real(r8), pointer :: m_deadstemc_to_litter_fire(:)
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)
   real(r8), pointer :: m_leafc_storage_to_litter(:)
   real(r8), pointer :: m_frootc_storage_to_litter(:)
   real(r8), pointer :: m_livestemc_storage_to_litter(:)
   real(r8), pointer :: m_deadstemc_storage_to_litter(:)
   real(r8), pointer :: m_livecrootc_storage_to_litter(:)
   real(r8), pointer :: m_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: m_gresp_storage_to_litter(:)
   real(r8), pointer :: m_leafc_storage_to_fire(:)
   real(r8), pointer :: m_frootc_storage_to_fire(:)
   real(r8), pointer :: m_livestemc_storage_to_fire(:)
   real(r8), pointer :: m_deadstemc_storage_to_fire(:)
   real(r8), pointer :: m_livecrootc_storage_to_fire(:)
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:)
   real(r8), pointer :: m_gresp_storage_to_fire(:)
   real(r8), pointer :: m_leafc_xfer_to_litter(:)
   real(r8), pointer :: m_frootc_xfer_to_litter(:)
   real(r8), pointer :: m_livestemc_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: m_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: m_gresp_xfer_to_litter(:)
   real(r8), pointer :: m_leafc_xfer_to_fire(:)
   real(r8), pointer :: m_frootc_xfer_to_fire(:)
   real(r8), pointer :: m_livestemc_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:)
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)
   real(r8), pointer :: m_gresp_xfer_to_fire(:)
   real(r8), pointer :: leaf_curmr(:)
   real(r8), pointer :: froot_curmr(:)
   real(r8), pointer :: livestem_curmr(:)
   real(r8), pointer :: livecroot_curmr(:)
   real(r8), pointer :: leaf_xsmr(:)
   real(r8), pointer :: froot_xsmr(:)
   real(r8), pointer :: livestem_xsmr(:)
   real(r8), pointer :: livecroot_xsmr(:)
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
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_cinputs(:)  ! (gC/m2/s) total column-level carbon inputs (for balance check)
   real(r8), pointer :: col_coutputs(:) ! (gC/m2/s) total column-level carbon outputs (for balance check)
   real(r8), pointer :: col_begcb(:)    ! carbon mass, beginning of time step (gC/m**2)
   real(r8), pointer :: col_endcb(:)    ! carbon mass, end of time step (gC/m**2)
   real(r8), pointer :: col_errcb(:)    ! carbon balance error for the timestep (gC/m**2)
   real(r8), pointer :: pft_cinputs(:)  ! (gC/m2/s) pft-level carbon inputs (for balance checking)
   real(r8), pointer :: pft_coutputs(:) ! (gC/m2/s) pft-level carbon outputs (for balance checking)
   real(r8), pointer :: pft_begcb(:)    ! carbon mass, beginning of time step (gC/m**2)
   real(r8), pointer :: pft_endcb(:)    ! carbon mass, end of time step (gC/m**2)
   real(r8), pointer :: pft_errcb(:)    ! carbon balance error for the timestep (gC/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,err_index  ! indices
   integer :: fp,fc          ! lake filter indices
   integer :: dtime          ! time step (s)
   logical :: err_found      ! error flag
   real(r8):: dt             ! radiation time step (seconds)
   real(r8):: t1,t2,t3,t4,t5 ! temporary variables for summation
!EOP
!-----------------------------------------------------------------------

    ! assign local pointers to column-level arrays

    cwdc                           => clm3%g%l%c%ccs%cwdc
    litr1c                         => clm3%g%l%c%ccs%litr1c
    litr2c                         => clm3%g%l%c%ccs%litr2c
    litr3c                         => clm3%g%l%c%ccs%litr3c
    soil1c                         => clm3%g%l%c%ccs%soil1c
    soil2c                         => clm3%g%l%c%ccs%soil2c
    soil3c                         => clm3%g%l%c%ccs%soil3c
    soil4c                         => clm3%g%l%c%ccs%soil4c
    col_ctrunc                     => clm3%g%l%c%ccs%col_ctrunc
    leafc_to_litr1c                => clm3%g%l%c%ccf%leafc_to_litr1c
    leafc_to_litr2c                => clm3%g%l%c%ccf%leafc_to_litr2c
    leafc_to_litr3c                => clm3%g%l%c%ccf%leafc_to_litr3c
    frootc_to_litr1c               => clm3%g%l%c%ccf%frootc_to_litr1c
    frootc_to_litr2c               => clm3%g%l%c%ccf%frootc_to_litr2c
    frootc_to_litr3c               => clm3%g%l%c%ccf%frootc_to_litr3c
    m_leafc_to_litr1c              => clm3%g%l%c%ccf%m_leafc_to_litr1c
    m_leafc_to_litr2c              => clm3%g%l%c%ccf%m_leafc_to_litr2c
    m_leafc_to_litr3c              => clm3%g%l%c%ccf%m_leafc_to_litr3c
    m_frootc_to_litr1c             => clm3%g%l%c%ccf%m_frootc_to_litr1c
    m_frootc_to_litr2c             => clm3%g%l%c%ccf%m_frootc_to_litr2c
    m_frootc_to_litr3c             => clm3%g%l%c%ccf%m_frootc_to_litr3c
    m_livestemc_to_cwdc            => clm3%g%l%c%ccf%m_livestemc_to_cwdc
    m_deadstemc_to_cwdc            => clm3%g%l%c%ccf%m_deadstemc_to_cwdc
    m_livecrootc_to_cwdc           => clm3%g%l%c%ccf%m_livecrootc_to_cwdc
    m_deadcrootc_to_cwdc           => clm3%g%l%c%ccf%m_deadcrootc_to_cwdc
    m_deadstemc_to_cwdc_fire       => clm3%g%l%c%ccf%m_deadstemc_to_cwdc_fire
    m_deadcrootc_to_cwdc_fire      => clm3%g%l%c%ccf%m_deadcrootc_to_cwdc_fire
    m_leafc_storage_to_litr1c      => clm3%g%l%c%ccf%m_leafc_storage_to_litr1c
    m_frootc_storage_to_litr1c     => clm3%g%l%c%ccf%m_frootc_storage_to_litr1c
    m_livestemc_storage_to_litr1c  => clm3%g%l%c%ccf%m_livestemc_storage_to_litr1c
    m_deadstemc_storage_to_litr1c  => clm3%g%l%c%ccf%m_deadstemc_storage_to_litr1c
    m_livecrootc_storage_to_litr1c => clm3%g%l%c%ccf%m_livecrootc_storage_to_litr1c
    m_deadcrootc_storage_to_litr1c => clm3%g%l%c%ccf%m_deadcrootc_storage_to_litr1c
    m_gresp_storage_to_litr1c      => clm3%g%l%c%ccf%m_gresp_storage_to_litr1c
    m_leafc_xfer_to_litr1c         => clm3%g%l%c%ccf%m_leafc_xfer_to_litr1c
    m_frootc_xfer_to_litr1c        => clm3%g%l%c%ccf%m_frootc_xfer_to_litr1c
    m_livestemc_xfer_to_litr1c     => clm3%g%l%c%ccf%m_livestemc_xfer_to_litr1c
    m_deadstemc_xfer_to_litr1c     => clm3%g%l%c%ccf%m_deadstemc_xfer_to_litr1c
    m_livecrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%m_livecrootc_xfer_to_litr1c
    m_deadcrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%m_deadcrootc_xfer_to_litr1c
    m_gresp_xfer_to_litr1c         => clm3%g%l%c%ccf%m_gresp_xfer_to_litr1c
    col_cinputs                    => clm3%g%l%c%ccf%col_cinputs
    col_coutputs                   => clm3%g%l%c%ccf%col_coutputs
    litr1_hr                       => clm3%g%l%c%ccf%litr1_hr
    litr2_hr                       => clm3%g%l%c%ccf%litr2_hr
    litr3_hr                       => clm3%g%l%c%ccf%litr3_hr
    soil1_hr                       => clm3%g%l%c%ccf%soil1_hr
    soil2_hr                       => clm3%g%l%c%ccf%soil2_hr
    soil3_hr                       => clm3%g%l%c%ccf%soil3_hr
    soil4_hr                       => clm3%g%l%c%ccf%soil4_hr
    m_litr1c_to_fire               => clm3%g%l%c%ccf%m_litr1c_to_fire
    m_litr2c_to_fire               => clm3%g%l%c%ccf%m_litr2c_to_fire
    m_litr3c_to_fire               => clm3%g%l%c%ccf%m_litr3c_to_fire
    m_cwdc_to_fire                 => clm3%g%l%c%ccf%m_cwdc_to_fire
    col_begcb                      => clm3%g%l%c%ccbal%begcb
    col_endcb                      => clm3%g%l%c%ccbal%endcb
    col_errcb                      => clm3%g%l%c%ccbal%errcb

    ! assign local pointers to pft-level arrays

    leafc                          => clm3%g%l%c%p%pcs%leafc
    leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    frootc                         => clm3%g%l%c%p%pcs%frootc
    frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    livestemc                      => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
    cpool                          => clm3%g%l%c%p%pcs%cpool
    xsmrpool                          => clm3%g%l%c%p%pcs%xsmrpool
    pft_ctrunc                     => clm3%g%l%c%p%pcs%pft_ctrunc
    psnsun_to_cpool                => clm3%g%l%c%p%pcf%psnsun_to_cpool
    psnshade_to_cpool              => clm3%g%l%c%p%pcf%psnshade_to_cpool
    leafc_to_litter                => clm3%g%l%c%p%pcf%leafc_to_litter
    frootc_to_litter               => clm3%g%l%c%p%pcf%frootc_to_litter
    m_leafc_to_litter              => clm3%g%l%c%p%pcf%m_leafc_to_litter
    m_frootc_to_litter             => clm3%g%l%c%p%pcf%m_frootc_to_litter
    m_livestemc_to_litter          => clm3%g%l%c%p%pcf%m_livestemc_to_litter
    m_deadstemc_to_litter          => clm3%g%l%c%p%pcf%m_deadstemc_to_litter
    m_livecrootc_to_litter         => clm3%g%l%c%p%pcf%m_livecrootc_to_litter
    m_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%m_deadcrootc_to_litter
    m_leafc_to_fire                => clm3%g%l%c%p%pcf%m_leafc_to_fire
    m_frootc_to_fire               => clm3%g%l%c%p%pcf%m_frootc_to_fire
    m_livestemc_to_fire            => clm3%g%l%c%p%pcf%m_livestemc_to_fire
    m_deadstemc_to_fire            => clm3%g%l%c%p%pcf%m_deadstemc_to_fire
    m_livecrootc_to_fire           => clm3%g%l%c%p%pcf%m_livecrootc_to_fire
    m_deadcrootc_to_fire           => clm3%g%l%c%p%pcf%m_deadcrootc_to_fire
    m_deadstemc_to_litter_fire     => clm3%g%l%c%p%pcf%m_deadstemc_to_litter_fire
    m_deadcrootc_to_litter_fire    => clm3%g%l%c%p%pcf%m_deadcrootc_to_litter_fire
    m_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%m_leafc_storage_to_litter
    m_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%m_frootc_storage_to_litter
    m_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter
    m_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter
    m_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter
    m_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter
    m_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%m_gresp_storage_to_litter
    m_leafc_storage_to_fire        => clm3%g%l%c%p%pcf%m_leafc_storage_to_fire
    m_frootc_storage_to_fire       => clm3%g%l%c%p%pcf%m_frootc_storage_to_fire
    m_livestemc_storage_to_fire    => clm3%g%l%c%p%pcf%m_livestemc_storage_to_fire
    m_deadstemc_storage_to_fire    => clm3%g%l%c%p%pcf%m_deadstemc_storage_to_fire
    m_livecrootc_storage_to_fire   => clm3%g%l%c%p%pcf%m_livecrootc_storage_to_fire
    m_deadcrootc_storage_to_fire   => clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_fire
    m_gresp_storage_to_fire        => clm3%g%l%c%p%pcf%m_gresp_storage_to_fire
    m_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter
    m_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter
    m_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter
    m_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter
    m_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter
    m_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter
    m_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter
    m_leafc_xfer_to_fire           => clm3%g%l%c%p%pcf%m_leafc_xfer_to_fire
    m_frootc_xfer_to_fire          => clm3%g%l%c%p%pcf%m_frootc_xfer_to_fire
    m_livestemc_xfer_to_fire       => clm3%g%l%c%p%pcf%m_livestemc_xfer_to_fire
    m_deadstemc_xfer_to_fire       => clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_fire
    m_livecrootc_xfer_to_fire      => clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_fire
    m_deadcrootc_xfer_to_fire      => clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_fire
    m_gresp_xfer_to_fire           => clm3%g%l%c%p%pcf%m_gresp_xfer_to_fire
    leaf_curmr                        => clm3%g%l%c%p%pcf%leaf_curmr
    froot_curmr                       => clm3%g%l%c%p%pcf%froot_curmr
    livestem_curmr                    => clm3%g%l%c%p%pcf%livestem_curmr
    livecroot_curmr                   => clm3%g%l%c%p%pcf%livecroot_curmr
    leaf_xsmr                        => clm3%g%l%c%p%pcf%leaf_xsmr
    froot_xsmr                       => clm3%g%l%c%p%pcf%froot_xsmr
    livestem_xsmr                    => clm3%g%l%c%p%pcf%livestem_xsmr
    livecroot_xsmr                   => clm3%g%l%c%p%pcf%livecroot_xsmr
    cpool_leaf_gr                  => clm3%g%l%c%p%pcf%cpool_leaf_gr
    cpool_leaf_storage_gr          => clm3%g%l%c%p%pcf%cpool_leaf_storage_gr
    transfer_leaf_gr               => clm3%g%l%c%p%pcf%transfer_leaf_gr
    cpool_froot_gr                 => clm3%g%l%c%p%pcf%cpool_froot_gr
    cpool_froot_storage_gr         => clm3%g%l%c%p%pcf%cpool_froot_storage_gr
    transfer_froot_gr              => clm3%g%l%c%p%pcf%transfer_froot_gr
    cpool_livestem_gr              => clm3%g%l%c%p%pcf%cpool_livestem_gr
    cpool_livestem_storage_gr      => clm3%g%l%c%p%pcf%cpool_livestem_storage_gr
    transfer_livestem_gr           => clm3%g%l%c%p%pcf%transfer_livestem_gr
    cpool_deadstem_gr              => clm3%g%l%c%p%pcf%cpool_deadstem_gr
    cpool_deadstem_storage_gr      => clm3%g%l%c%p%pcf%cpool_deadstem_storage_gr
    transfer_deadstem_gr           => clm3%g%l%c%p%pcf%transfer_deadstem_gr
    cpool_livecroot_gr             => clm3%g%l%c%p%pcf%cpool_livecroot_gr
    cpool_livecroot_storage_gr     => clm3%g%l%c%p%pcf%cpool_livecroot_storage_gr
    transfer_livecroot_gr          => clm3%g%l%c%p%pcf%transfer_livecroot_gr
    cpool_deadcroot_gr             => clm3%g%l%c%p%pcf%cpool_deadcroot_gr
    cpool_deadcroot_storage_gr     => clm3%g%l%c%p%pcf%cpool_deadcroot_storage_gr
    transfer_deadcroot_gr          => clm3%g%l%c%p%pcf%transfer_deadcroot_gr
    pft_cinputs                    => clm3%g%l%c%p%pcf%pft_cinputs
    pft_coutputs                   => clm3%g%l%c%p%pcf%pft_coutputs
    pft_begcb                      => clm3%g%l%c%p%pcbal%begcb
    pft_endcb                      => clm3%g%l%c%p%pcbal%endcb
    pft_errcb                      => clm3%g%l%c%p%pcbal%errcb

   ! set time steps
   dtime = get_step_size()
   dt = float(irad)*dtime

   err_found = .false.
   ! column loop
!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! calculate the total column-level carbon storage, for mass conservation check

      col_endcb(c) = cwdc(c) + litr1c(c) + litr2c(c) + litr3c(c) + &
         soil1c(c) + soil2c(c) + soil3c(c) + soil4c(c) + col_ctrunc(c)

      ! calculate total column-level inputs

      t1 = &
         leafc_to_litr1c(c)    + leafc_to_litr2c(c)    + leafc_to_litr3c(c) + &
         frootc_to_litr1c(c)   + frootc_to_litr2c(c)   + frootc_to_litr3c(c) + &
         m_leafc_to_litr1c(c)  + m_leafc_to_litr2c(c)  + m_leafc_to_litr3c(c) + &
         m_frootc_to_litr1c(c) + m_frootc_to_litr2c(c) + m_frootc_to_litr3c(c)
      t2 = &
         m_livestemc_to_cwdc(c)      + m_deadstemc_to_cwdc(c) + &
         m_livecrootc_to_cwdc(c)     + m_deadcrootc_to_cwdc(c) + &
         m_deadstemc_to_cwdc_fire(c) + m_deadcrootc_to_cwdc_fire(c)
      t3 = &
         m_leafc_storage_to_litr1c(c)      + m_frootc_storage_to_litr1c(c) + &
         m_livestemc_storage_to_litr1c(c)  + m_deadstemc_storage_to_litr1c(c) + &
         m_livecrootc_storage_to_litr1c(c) + m_deadcrootc_storage_to_litr1c(c) + &
         m_gresp_storage_to_litr1c(c)
      t4 = &
         m_leafc_xfer_to_litr1c(c)      + m_frootc_xfer_to_litr1c(c) + &
         m_livestemc_xfer_to_litr1c(c)  + m_deadstemc_xfer_to_litr1c(c) + &
         m_livecrootc_xfer_to_litr1c(c) + m_deadcrootc_xfer_to_litr1c(c) + &
         m_gresp_xfer_to_litr1c(c)

      col_cinputs(c) = t1+t2+t3+t4

      ! calculate total column-level outputs

      col_coutputs(c) = &
         litr1_hr(c)         + litr2_hr(c)         + litr3_hr(c) + &
         soil1_hr(c)         + soil2_hr(c)         + soil3_hr(c) + soil4_hr(c) + &
         m_litr1c_to_fire(c) + m_litr2c_to_fire(c) + m_litr3c_to_fire(c) + &
         m_cwdc_to_fire(c)

      ! calculate the total column-level carbon balance error for this time step

      col_errcb(c) = (col_cinputs(c) - col_coutputs(c))*dt - &
         (col_endcb(c) - col_begcb(c))

      ! check for significant errors
      if (abs(col_errcb(c)) > 1e-8_r8) then
         err_found = .true.
         err_index = c
      end if
      
   end do ! end of columns loop

   if (err_found) then
      c = err_index
      write(6,*)'column cbalance error = ', col_errcb(c), c
      write(6,*)'begcb       = ',col_begcb(c)
      write(6,*)'endcb       = ',col_endcb(c)
      write(6,*)'delta store = ',col_endcb(c)-col_begcb(c)
      write(6,*)'input mass  = ',col_cinputs(c)*dt
      write(6,*)'output mass = ',col_coutputs(c)*dt
      write(6,*)'net flux    = ',(col_cinputs(c)-col_coutputs(c))*dt
      call endrun
   end if

   err_found = .false.
   ! pft loop
!dir$ concurrent
!cdir nodep
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! calculate the total pft-level carbon storage, for mass conservation check

      pft_endcb(p) = leafc(p)      + leafc_storage(p)      + leafc_xfer(p) + &
         frootc(p)     + frootc_storage(p)     + frootc_xfer(p) + &
         livestemc(p)  + livestemc_storage(p)  + livestemc_xfer(p) + &
         deadstemc(p)  + deadstemc_storage(p)  + deadstemc_xfer(p) + &
         livecrootc(p) + livecrootc_storage(p) + livecrootc_xfer(p) + &
         deadcrootc(p) + deadcrootc_storage(p) + deadcrootc_xfer(p) + &
         gresp_storage(p) + gresp_xfer(p)  + cpool(p) + xsmrpool(p) + &
         pft_ctrunc(p)

      ! calculate total pft-level inputs

      pft_cinputs(p) = psnsun_to_cpool(p) + psnshade_to_cpool(p)

      ! calculate total pft-level outputs

      t1 = &
         leafc_to_litter(p)        + frootc_to_litter(p) + &
         m_leafc_to_litter(p)      + m_frootc_to_litter(p) + &
         m_livestemc_to_litter(p)  + m_deadstemc_to_litter(p) + &
         m_livecrootc_to_litter(p) + m_deadcrootc_to_litter(p) + &
         m_leafc_to_fire(p)        + m_frootc_to_fire(p) + &
         m_livestemc_to_fire(p)    + m_deadstemc_to_fire(p) + &
         m_livecrootc_to_fire(p)   + m_deadcrootc_to_fire(p) + &
         m_deadstemc_to_litter_fire(p) + m_deadcrootc_to_litter_fire(p)
      t2 = &
         m_leafc_storage_to_litter(p)      + m_frootc_storage_to_litter(p) + &
         m_livestemc_storage_to_litter(p)  + m_deadstemc_storage_to_litter(p) + &
         m_livecrootc_storage_to_litter(p) + m_deadcrootc_storage_to_litter(p) + &
         m_gresp_storage_to_litter(p)      + &
         m_leafc_storage_to_fire(p)        + m_frootc_storage_to_fire(p) + &
         m_livestemc_storage_to_fire(p)    + m_deadstemc_storage_to_fire(p) + &
         m_livecrootc_storage_to_fire(p)   + m_deadcrootc_storage_to_fire(p) + &
         m_gresp_storage_to_fire(p)
     t3 = &
         m_leafc_xfer_to_litter(p)      + m_frootc_xfer_to_litter(p) + &
         m_livestemc_xfer_to_litter(p)  + m_deadstemc_xfer_to_litter(p) + &
         m_livecrootc_xfer_to_litter(p) + m_deadcrootc_xfer_to_litter(p) + &
         m_gresp_xfer_to_litter(p)      + &
         m_leafc_xfer_to_fire(p)        + m_frootc_xfer_to_fire(p) + &
         m_livestemc_xfer_to_fire(p)    + m_deadstemc_xfer_to_fire(p) + &
         m_livecrootc_xfer_to_fire(p)   + m_deadcrootc_xfer_to_fire(p) + &
         m_gresp_xfer_to_fire(p)
      t4 = &
         leaf_curmr(p) + froot_curmr(p) + livestem_curmr(p) + livecroot_curmr(p) + &
         leaf_xsmr(p) + froot_xsmr(p) + livestem_xsmr(p) + livecroot_xsmr(p) + &
         cpool_leaf_gr(p)      + cpool_leaf_storage_gr(p)      + transfer_leaf_gr(p) + &
         cpool_froot_gr(p)     + cpool_froot_storage_gr(p)     + transfer_froot_gr(p)
      t5 = &
         cpool_livestem_gr(p)  + cpool_livestem_storage_gr(p)  + transfer_livestem_gr(p) + &
         cpool_deadstem_gr(p)  + cpool_deadstem_storage_gr(p)  + transfer_deadstem_gr(p) + &
         cpool_livecroot_gr(p) + cpool_livecroot_storage_gr(p) + transfer_livecroot_gr(p) + &
         cpool_deadcroot_gr(p) + cpool_deadcroot_storage_gr(p) + transfer_deadcroot_gr(p)

      pft_coutputs(p) = t1+t2+t3+t4+t5
      ! calculate the total pft-level carbon balance error for this time step

      pft_errcb(p) = (pft_cinputs(p) - pft_coutputs(p))*dt - &
         (pft_endcb(p) - pft_begcb(p))

      ! check for significant errors
      if (abs(pft_errcb(p)) > 1e-8_r8) then
         err_found = .true.
         err_index = p
      end if

   end do  ! end of pfts loop

   if (err_found) then
      p = err_index
      write(6,*)'pft cbalance error = ',pft_errcb(p),p
      write(6,*)'begcb       = ',pft_begcb(p)
      write(6,*)'endcb       = ',pft_endcb(p)
      write(6,*)'delta store = ',pft_endcb(p)-pft_begcb(p)
      write(6,*)'input mass  = ',pft_cinputs(p)*dt
      write(6,*)'output mass = ',pft_coutputs(p)*dt
      write(6,*)'net flux    = ',(pft_cinputs(p)-pft_coutputs(p))*dt
      write(6,*)'veg type    = ',clm3%g%l%c%p%itype(p)
      call endrun
   end if

end subroutine CBalanceCheck
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NBalanceCheck
!
! !INTERFACE:
subroutine NBalanceCheck(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, perform nitrogen mass conservation check
! for column and pft
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varctl, only: irad
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
! 12/9/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real(r8), pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real(r8), pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real(r8), pointer :: soil1n(:)             ! (gN/m2) soil organic matter N (fast pool)
   real(r8), pointer :: soil2n(:)             ! (gN/m2) soil organic matter N (medium pool)
   real(r8), pointer :: soil3n(:)             ! (gN/m2) soil orgainc matter N (slow pool)
   real(r8), pointer :: soil4n(:)             ! (gN/m2) soil orgainc matter N (slowest pool)
   real(r8), pointer :: sminn(:)              ! (gN/m2) soil mineral N
   real(r8), pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   real(r8), pointer :: ndep_to_sminn(:)
   real(r8), pointer :: nfix_to_sminn(:)
   real(r8), pointer :: leafn_to_litr1n(:)
   real(r8), pointer :: leafn_to_litr2n(:)
   real(r8), pointer :: leafn_to_litr3n(:)
   real(r8), pointer :: frootn_to_litr1n(:)
   real(r8), pointer :: frootn_to_litr2n(:)
   real(r8), pointer :: frootn_to_litr3n(:)
   real(r8), pointer :: supplement_to_sminn(:)
   real(r8), pointer :: m_leafn_to_litr1n(:)
   real(r8), pointer :: m_leafn_to_litr2n(:)
   real(r8), pointer :: m_leafn_to_litr3n(:)
   real(r8), pointer :: m_frootn_to_litr1n(:)
   real(r8), pointer :: m_frootn_to_litr2n(:)
   real(r8), pointer :: m_frootn_to_litr3n(:)
   real(r8), pointer :: m_livestemn_to_cwdn(:)
   real(r8), pointer :: m_deadstemn_to_cwdn(:)
   real(r8), pointer :: m_livecrootn_to_cwdn(:)
   real(r8), pointer :: m_deadcrootn_to_cwdn(:)
   real(r8), pointer :: m_deadstemn_to_cwdn_fire(:)
   real(r8), pointer :: m_deadcrootn_to_cwdn_fire(:)
   real(r8), pointer :: m_retransn_to_litr1n(:)
   real(r8), pointer :: m_leafn_storage_to_litr1n(:)
   real(r8), pointer :: m_frootn_storage_to_litr1n(:)
   real(r8), pointer :: m_livestemn_storage_to_litr1n(:)
   real(r8), pointer :: m_deadstemn_storage_to_litr1n(:)
   real(r8), pointer :: m_livecrootn_storage_to_litr1n(:)
   real(r8), pointer :: m_deadcrootn_storage_to_litr1n(:)
   real(r8), pointer :: m_leafn_xfer_to_litr1n(:)
   real(r8), pointer :: m_frootn_xfer_to_litr1n(:)
   real(r8), pointer :: m_livestemn_xfer_to_litr1n(:)
   real(r8), pointer :: m_deadstemn_xfer_to_litr1n(:)
   real(r8), pointer :: m_livecrootn_xfer_to_litr1n(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_litr1n(:)
   real(r8), pointer :: sminn_to_denit_l1s1(:)
   real(r8), pointer :: sminn_to_denit_l2s2(:)
   real(r8), pointer :: sminn_to_denit_l3s3(:)
   real(r8), pointer :: sminn_to_denit_s1s2(:)
   real(r8), pointer :: sminn_to_denit_s2s3(:)
   real(r8), pointer :: sminn_to_denit_s3s4(:)
   real(r8), pointer :: sminn_to_denit_s4(:)
   real(r8), pointer :: sminn_to_denit_excess(:)
   real(r8), pointer :: sminn_to_plant(:)
   real(r8), pointer :: sminn_leached(:) 
   real(r8), pointer :: m_litr1n_to_fire(:)
   real(r8), pointer :: m_litr2n_to_fire(:)
   real(r8), pointer :: m_litr3n_to_fire(:)
   real(r8), pointer :: m_cwdn_to_fire(:)
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real(r8), pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   real(r8), pointer :: sminn_to_npool(:)
   real(r8), pointer :: leafn_to_litter(:)
   real(r8), pointer :: frootn_to_litter(:)
   real(r8), pointer :: m_leafn_to_litter(:)
   real(r8), pointer :: m_frootn_to_litter(:)
   real(r8), pointer :: m_livestemn_to_litter(:)
   real(r8), pointer :: m_deadstemn_to_litter(:)
   real(r8), pointer :: m_livecrootn_to_litter(:)
   real(r8), pointer :: m_deadcrootn_to_litter(:)
   real(r8), pointer :: m_retransn_to_litter(:)
   real(r8), pointer :: m_leafn_to_fire(:)
   real(r8), pointer :: m_frootn_to_fire(:)
   real(r8), pointer :: m_livestemn_to_fire(:)
   real(r8), pointer :: m_deadstemn_to_fire(:)
   real(r8), pointer :: m_livecrootn_to_fire(:)
   real(r8), pointer :: m_deadcrootn_to_fire(:)
   real(r8), pointer :: m_deadstemn_to_litter_fire(:)
   real(r8), pointer :: m_deadcrootn_to_litter_fire(:)
   real(r8), pointer :: m_retransn_to_fire(:)
   real(r8), pointer :: m_leafn_storage_to_litter(:)
   real(r8), pointer :: m_frootn_storage_to_litter(:)
   real(r8), pointer :: m_livestemn_storage_to_litter(:)
   real(r8), pointer :: m_deadstemn_storage_to_litter(:)
   real(r8), pointer :: m_livecrootn_storage_to_litter(:)
   real(r8), pointer :: m_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: m_leafn_storage_to_fire(:)
   real(r8), pointer :: m_frootn_storage_to_fire(:)
   real(r8), pointer :: m_livestemn_storage_to_fire(:)
   real(r8), pointer :: m_deadstemn_storage_to_fire(:)
   real(r8), pointer :: m_livecrootn_storage_to_fire(:)
   real(r8), pointer :: m_deadcrootn_storage_to_fire(:)
   real(r8), pointer :: m_leafn_xfer_to_litter(:)
   real(r8), pointer :: m_frootn_xfer_to_litter(:)
   real(r8), pointer :: m_livestemn_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: m_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_litter(:)
   real(r8), pointer :: m_leafn_xfer_to_fire(:)
   real(r8), pointer :: m_frootn_xfer_to_fire(:)
   real(r8), pointer :: m_livestemn_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemn_xfer_to_fire(:)
   real(r8), pointer :: m_livecrootn_xfer_to_fire(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_fire(:)
!
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_ninputs(:)
   real(r8), pointer :: col_noutputs(:)
   real(r8), pointer :: col_begnb(:)    ! nitrogen mass, beginning of time step (gN/m**2)
   real(r8), pointer :: col_endnb(:)    ! nitrogen mass, end of time step (gN/m**2)
   real(r8), pointer :: col_errnb(:)    ! nitrogen balance error for the timestep (gN/m**2)
   real(r8), pointer :: pft_ninputs(:)
   real(r8), pointer :: pft_noutputs(:)
   real(r8), pointer :: pft_begnb(:)    ! nitrogen mass, beginning of time step (gN/m**2)
   real(r8), pointer :: pft_endnb(:)    ! nitrogen mass, end of time step (gN/m**2)
   real(r8), pointer :: pft_errnb(:)    ! nitrogen balance error for the timestep (gN/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,err_index  ! indices
   integer :: fp,fc          ! lake filter indices
   integer :: dtime          ! time step (s)
   logical :: err_found      ! error flag
   real(r8):: dt             ! radiation time step (seconds)
   real(r8):: t1,t2,t3,t4,t5 ! temporary variables for summation
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers to column-level arrays

    cwdn                           => clm3%g%l%c%cns%cwdn
    litr1n                         => clm3%g%l%c%cns%litr1n
    litr2n                         => clm3%g%l%c%cns%litr2n
    litr3n                         => clm3%g%l%c%cns%litr3n
    soil1n                         => clm3%g%l%c%cns%soil1n
    soil2n                         => clm3%g%l%c%cns%soil2n
    soil3n                         => clm3%g%l%c%cns%soil3n
    soil4n                         => clm3%g%l%c%cns%soil4n
    sminn                          => clm3%g%l%c%cns%sminn
    col_ntrunc                     => clm3%g%l%c%cns%col_ntrunc
    ndep_to_sminn                  => clm3%g%l%c%cnf%ndep_to_sminn
    nfix_to_sminn                  => clm3%g%l%c%cnf%nfix_to_sminn
    leafn_to_litr1n                => clm3%g%l%c%cnf%leafn_to_litr1n
    leafn_to_litr2n                => clm3%g%l%c%cnf%leafn_to_litr2n
    leafn_to_litr3n                => clm3%g%l%c%cnf%leafn_to_litr3n
    frootn_to_litr1n               => clm3%g%l%c%cnf%frootn_to_litr1n
    frootn_to_litr2n               => clm3%g%l%c%cnf%frootn_to_litr2n
    frootn_to_litr3n               => clm3%g%l%c%cnf%frootn_to_litr3n
    supplement_to_sminn            => clm3%g%l%c%cnf%supplement_to_sminn
    m_leafn_to_litr1n              => clm3%g%l%c%cnf%m_leafn_to_litr1n
    m_leafn_to_litr2n              => clm3%g%l%c%cnf%m_leafn_to_litr2n
    m_leafn_to_litr3n              => clm3%g%l%c%cnf%m_leafn_to_litr3n
    m_frootn_to_litr1n             => clm3%g%l%c%cnf%m_frootn_to_litr1n
    m_frootn_to_litr2n             => clm3%g%l%c%cnf%m_frootn_to_litr2n
    m_frootn_to_litr3n             => clm3%g%l%c%cnf%m_frootn_to_litr3n
    m_livestemn_to_cwdn            => clm3%g%l%c%cnf%m_livestemn_to_cwdn
    m_deadstemn_to_cwdn            => clm3%g%l%c%cnf%m_deadstemn_to_cwdn
    m_livecrootn_to_cwdn           => clm3%g%l%c%cnf%m_livecrootn_to_cwdn
    m_deadcrootn_to_cwdn           => clm3%g%l%c%cnf%m_deadcrootn_to_cwdn
    m_deadstemn_to_cwdn_fire       => clm3%g%l%c%cnf%m_deadstemn_to_cwdn_fire
    m_deadcrootn_to_cwdn_fire      => clm3%g%l%c%cnf%m_deadcrootn_to_cwdn_fire
    m_retransn_to_litr1n           => clm3%g%l%c%cnf%m_retransn_to_litr1n
    m_leafn_storage_to_litr1n      => clm3%g%l%c%cnf%m_leafn_storage_to_litr1n
    m_frootn_storage_to_litr1n     => clm3%g%l%c%cnf%m_frootn_storage_to_litr1n
    m_livestemn_storage_to_litr1n  => clm3%g%l%c%cnf%m_livestemn_storage_to_litr1n
    m_deadstemn_storage_to_litr1n  => clm3%g%l%c%cnf%m_deadstemn_storage_to_litr1n
    m_livecrootn_storage_to_litr1n => clm3%g%l%c%cnf%m_livecrootn_storage_to_litr1n
    m_deadcrootn_storage_to_litr1n => clm3%g%l%c%cnf%m_deadcrootn_storage_to_litr1n
    m_leafn_xfer_to_litr1n         => clm3%g%l%c%cnf%m_leafn_xfer_to_litr1n
    m_frootn_xfer_to_litr1n        => clm3%g%l%c%cnf%m_frootn_xfer_to_litr1n
    m_livestemn_xfer_to_litr1n     => clm3%g%l%c%cnf%m_livestemn_xfer_to_litr1n
    m_deadstemn_xfer_to_litr1n     => clm3%g%l%c%cnf%m_deadstemn_xfer_to_litr1n
    m_livecrootn_xfer_to_litr1n    => clm3%g%l%c%cnf%m_livecrootn_xfer_to_litr1n
    m_deadcrootn_xfer_to_litr1n    => clm3%g%l%c%cnf%m_deadcrootn_xfer_to_litr1n
    sminn_to_denit_l1s1            => clm3%g%l%c%cnf%sminn_to_denit_l1s1
    sminn_to_denit_l2s2            => clm3%g%l%c%cnf%sminn_to_denit_l2s2
    sminn_to_denit_l3s3            => clm3%g%l%c%cnf%sminn_to_denit_l3s3
    sminn_to_denit_s1s2            => clm3%g%l%c%cnf%sminn_to_denit_s1s2
    sminn_to_denit_s2s3            => clm3%g%l%c%cnf%sminn_to_denit_s2s3
    sminn_to_denit_s3s4            => clm3%g%l%c%cnf%sminn_to_denit_s3s4
    sminn_to_denit_s4              => clm3%g%l%c%cnf%sminn_to_denit_s4
    sminn_to_denit_excess          => clm3%g%l%c%cnf%sminn_to_denit_excess
    sminn_to_plant                 => clm3%g%l%c%cnf%sminn_to_plant
    sminn_leached                  => clm3%g%l%c%cnf%sminn_leached
    m_litr1n_to_fire               => clm3%g%l%c%cnf%m_litr1n_to_fire
    m_litr2n_to_fire               => clm3%g%l%c%cnf%m_litr2n_to_fire
    m_litr3n_to_fire               => clm3%g%l%c%cnf%m_litr3n_to_fire
    m_cwdn_to_fire                 => clm3%g%l%c%cnf%m_cwdn_to_fire
    col_ninputs                    => clm3%g%l%c%cnf%col_ninputs
    col_noutputs                   => clm3%g%l%c%cnf%col_noutputs
    col_begnb                      => clm3%g%l%c%cnbal%begnb
    col_endnb                      => clm3%g%l%c%cnbal%endnb
    col_errnb                      => clm3%g%l%c%cnbal%errnb

    ! assign local pointers to pft-level arrays

    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    retransn                       => clm3%g%l%c%p%pns%retransn
    npool                          => clm3%g%l%c%p%pns%npool
    pft_ntrunc                     => clm3%g%l%c%p%pns%pft_ntrunc
    sminn_to_npool                 => clm3%g%l%c%p%pnf%sminn_to_npool
    leafn_to_litter                => clm3%g%l%c%p%pnf%leafn_to_litter
    frootn_to_litter               => clm3%g%l%c%p%pnf%frootn_to_litter
    m_leafn_to_litter              => clm3%g%l%c%p%pnf%m_leafn_to_litter
    m_frootn_to_litter             => clm3%g%l%c%p%pnf%m_frootn_to_litter
    m_livestemn_to_litter          => clm3%g%l%c%p%pnf%m_livestemn_to_litter
    m_deadstemn_to_litter          => clm3%g%l%c%p%pnf%m_deadstemn_to_litter
    m_livecrootn_to_litter         => clm3%g%l%c%p%pnf%m_livecrootn_to_litter
    m_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter
    m_retransn_to_litter           => clm3%g%l%c%p%pnf%m_retransn_to_litter
    m_leafn_to_fire                => clm3%g%l%c%p%pnf%m_leafn_to_fire
    m_frootn_to_fire               => clm3%g%l%c%p%pnf%m_frootn_to_fire
    m_livestemn_to_fire            => clm3%g%l%c%p%pnf%m_livestemn_to_fire
    m_deadstemn_to_fire            => clm3%g%l%c%p%pnf%m_deadstemn_to_fire
    m_livecrootn_to_fire           => clm3%g%l%c%p%pnf%m_livecrootn_to_fire
    m_deadcrootn_to_fire           => clm3%g%l%c%p%pnf%m_deadcrootn_to_fire
    m_deadstemn_to_litter_fire     => clm3%g%l%c%p%pnf%m_deadstemn_to_litter_fire
    m_deadcrootn_to_litter_fire    => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter_fire
    m_retransn_to_fire             => clm3%g%l%c%p%pnf%m_retransn_to_fire
    m_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%m_leafn_storage_to_litter
    m_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%m_frootn_storage_to_litter
    m_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter
    m_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter
    m_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter
    m_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter
    m_leafn_storage_to_fire        => clm3%g%l%c%p%pnf%m_leafn_storage_to_fire
    m_frootn_storage_to_fire       => clm3%g%l%c%p%pnf%m_frootn_storage_to_fire
    m_livestemn_storage_to_fire    => clm3%g%l%c%p%pnf%m_livestemn_storage_to_fire
    m_deadstemn_storage_to_fire    => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_fire
    m_livecrootn_storage_to_fire   => clm3%g%l%c%p%pnf%m_livecrootn_storage_to_fire
    m_deadcrootn_storage_to_fire   => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_fire
    m_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter
    m_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter
    m_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter
    m_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter
    m_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter
    m_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter
    m_leafn_xfer_to_fire           => clm3%g%l%c%p%pnf%m_leafn_xfer_to_fire
    m_frootn_xfer_to_fire          => clm3%g%l%c%p%pnf%m_frootn_xfer_to_fire
    m_livestemn_xfer_to_fire       => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_fire
    m_deadstemn_xfer_to_fire       => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_fire
    m_livecrootn_xfer_to_fire      => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_fire
    m_deadcrootn_xfer_to_fire      => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_fire
    pft_ninputs                    => clm3%g%l%c%p%pnf%pft_ninputs
    pft_noutputs                   => clm3%g%l%c%p%pnf%pft_noutputs
    pft_begnb                      => clm3%g%l%c%p%pnbal%begnb
    pft_endnb                      => clm3%g%l%c%p%pnbal%endnb
    pft_errnb                      => clm3%g%l%c%p%pnbal%errnb

   ! set time steps
   dtime = get_step_size()
   dt = float(irad)*dtime

   err_found = .false.
   ! column loop
!dir$ concurrent
!cdir nodep
   do fc = 1,num_soilc
      c=filter_soilc(fc)

      ! calculate the total column-level nitrogen storage, for mass conservation check

      col_endnb(c) = cwdn(c) + &
         litr1n(c) + litr2n(c) + litr3n(c) + &
         soil1n(c) + soil2n(c) + soil3n(c) + soil4n(c) + &
         sminn(c)  + col_ntrunc(c)

      ! calculate total column-level inputs

      t1 = &
         ndep_to_sminn(c) + nfix_to_sminn(c) + &
         leafn_to_litr1n(c)  + leafn_to_litr2n(c)  + leafn_to_litr3n(c)  + &
         frootn_to_litr1n(c) + frootn_to_litr2n(c) + frootn_to_litr3n(c) + &
         supplement_to_sminn(c)
      t2 = &
         m_leafn_to_litr1n(c)  + m_leafn_to_litr2n(c)  + m_leafn_to_litr3n(c) + &
         m_frootn_to_litr1n(c) + m_frootn_to_litr2n(c) + m_frootn_to_litr3n(c)
      t3 = &
         m_livestemn_to_cwdn(c)      + m_deadstemn_to_cwdn(c) + &
         m_livecrootn_to_cwdn(c)     + m_deadcrootn_to_cwdn(c) + &
         m_deadstemn_to_cwdn_fire(c) + m_deadcrootn_to_cwdn_fire(c) + &
         m_retransn_to_litr1n(c)
      t4 = &
         m_leafn_storage_to_litr1n(c)      + m_frootn_storage_to_litr1n(c) + &
         m_livestemn_storage_to_litr1n(c)  + m_deadstemn_storage_to_litr1n(c) + &
         m_livecrootn_storage_to_litr1n(c) + m_deadcrootn_storage_to_litr1n(c)
      t5 = &
         m_leafn_xfer_to_litr1n(c)      + m_frootn_xfer_to_litr1n(c) + &
         m_livestemn_xfer_to_litr1n(c)  + m_deadstemn_xfer_to_litr1n(c) + &
         m_livecrootn_xfer_to_litr1n(c) + m_deadcrootn_xfer_to_litr1n(c)

      col_ninputs(c) = t1+t2+t3+t4+t5

      ! calculate total column-level outputs

      col_noutputs(c) = &
         sminn_to_denit_l1s1(c)   + sminn_to_denit_l2s2(c) + sminn_to_denit_l3s3(c) + &
         sminn_to_denit_s1s2(c)   + sminn_to_denit_s2s3(c) + sminn_to_denit_s3s4(c) + &
         sminn_to_denit_s4(c) + &
         sminn_to_denit_excess(c) + sminn_to_plant(c)      + sminn_leached(c)       + &
         m_litr1n_to_fire(c) + m_litr2n_to_fire(c) + m_litr3n_to_fire(c) + &
         m_cwdn_to_fire(c)

      ! calculate the total column-level nitrogen balance error for this time step

      col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - &
         (col_endnb(c) - col_begnb(c))

      if (abs(col_errnb(c)) > 1e-8_r8) then
         err_found = .true.
         err_index = c
      end if

   end do ! end of columns loop

   if (err_found) then
      c = err_index
      write(6,*)'column nbalance error = ', col_errnb(c), c
      write(6,*)'begnb       = ',col_begnb(c)
      write(6,*)'endnb       = ',col_endnb(c)
      write(6,*)'delta store = ',col_endnb(c)-col_begnb(c)
      write(6,*)'input mass  = ',col_ninputs(c)*dt
      write(6,*)'output mass = ',col_noutputs(c)*dt
      write(6,*)'net flux    = ',(col_ninputs(c)-col_noutputs(c))*dt
      call endrun
   end if

   err_found = .false.
   ! pft loop
!dir$ concurrent
!cdir nodep
   do fp=1,num_soilp
      p=filter_soilp(fp)

      ! calculate the total pft-level nitrogen storage, for mass conservation check

      pft_endnb(p) = &
         leafn(p)      + leafn_storage(p)      + leafn_xfer(p) + &
         frootn(p)     + frootn_storage(p)     + frootn_xfer(p) + &
         livestemn(p)  + livestemn_storage(p)  + livestemn_xfer(p) + &
         deadstemn(p)  + deadstemn_storage(p)  + deadstemn_xfer(p) + &
         livecrootn(p) + livecrootn_storage(p) + livecrootn_xfer(p) + &
         deadcrootn(p) + deadcrootn_storage(p) + deadcrootn_xfer(p) + &
         retransn(p)   + npool(p)              + pft_ntrunc(p)

      ! calculate total pft-level inputs

      pft_ninputs(p) = sminn_to_npool(p)

      ! calculate total pft-level outputs

      t1 = leafn_to_litter(p) + frootn_to_litter(p)
      t2 = &
         m_leafn_to_litter(p)      + m_frootn_to_litter(p) + &
         m_livestemn_to_litter(p)  + m_deadstemn_to_litter(p) + &
         m_livecrootn_to_litter(p) + m_deadcrootn_to_litter(p) + &
         m_retransn_to_litter(p)   + &
         m_leafn_to_fire(p)        + m_frootn_to_fire(p) + &
         m_livestemn_to_fire(p)    + m_deadstemn_to_fire(p) + &
         m_livecrootn_to_fire(p)   + m_deadcrootn_to_fire(p) + &
         m_deadstemn_to_litter_fire(p) + m_deadcrootn_to_litter_fire(p) + &
         m_retransn_to_fire(p)
      t3 = &
         m_leafn_storage_to_litter(p)      + m_frootn_storage_to_litter(p) + &
         m_livestemn_storage_to_litter(p)  + m_deadstemn_storage_to_litter(p) + &
         m_livecrootn_storage_to_litter(p) + m_deadcrootn_storage_to_litter(p) + &
         m_leafn_storage_to_fire(p)        + m_frootn_storage_to_fire(p) + &
         m_livestemn_storage_to_fire(p)    + m_deadstemn_storage_to_fire(p) + &
         m_livecrootn_storage_to_fire(p)   + m_deadcrootn_storage_to_fire(p)
      t4 = &
         m_leafn_xfer_to_litter(p)      + m_frootn_xfer_to_litter(p) + &
         m_livestemn_xfer_to_litter(p)  + m_deadstemn_xfer_to_litter(p) + &
         m_livecrootn_xfer_to_litter(p) + m_deadcrootn_xfer_to_litter(p) + &
         m_leafn_xfer_to_fire(p)        + m_frootn_xfer_to_fire(p) + &
         m_livestemn_xfer_to_fire(p)    + m_deadstemn_xfer_to_fire(p) + &
         m_livecrootn_xfer_to_fire(p)   + m_deadcrootn_xfer_to_fire(p)

      pft_noutputs(p) = t1+t2+t3+t4

      ! calculate the total pft-level nitrogen balance error for this time step

      pft_errnb(p) = (pft_ninputs(p) - pft_noutputs(p))*dt - &
         (pft_endnb(p) - pft_begnb(p))

      if (abs(pft_errnb(p)) > 1e-8_r8) then
         err_found = .true.
         err_index = p
      end if

   end do  ! end of pfts loop

   if (err_found) then
      p = err_index
      write(6,*)'pft nbalance error = ', pft_errnb(p), p
      write(6,*)'begnb       = ',pft_begnb(p)
      write(6,*)'endnb       = ',pft_endnb(p)
      write(6,*)'delta store = ',pft_endnb(p)-pft_begnb(p)
      write(6,*)'input mass  = ',pft_ninputs(p)*dt
      write(6,*)'output mass = ',pft_noutputs(p)*dt
      write(6,*)'net flux    = ',(pft_ninputs(p)-pft_noutputs(p))*dt
      write(6,*)'veg type    = ',clm3%g%l%c%p%itype(p)
      call endrun
   end if

end subroutine NBalanceCheck
!-----------------------------------------------------------------------

end module CNBalanceCheckMod
