module CASAiniTimeVarMod

#if (defined CASA)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CASAiniTimeVarMod
!
! !DESCRIPTION:
! Module for initializing time vary variables used only in CASA'
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: CASAiniTimeVar
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
! !IROUTINE: CASAiniTimeVar
!
! !INTERFACE:
subroutine CASAiniTimeVar()
!
! !DESCRIPTION:
! Initialize time varying variables used only in CASA'
!
! !USES:
   use clmtype
   use shr_kind_mod, only: r8 => shr_kind_r8
   use decompMod   , only: get_proc_bounds
   use clm_varpar  , only: nlevgrnd
!
! !ARGUMENTS:
   implicit none
!
! !CALLED FROM:
! subroutine initialize in file initializeMod.F90
!
! !REVISION HISTORY:
!  22 Sept 2006: Created by Forrest Hoffman
!
!
! local pointers to implicit in arguments
!
!
! local pointers to implicit out arguments
!
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
   real(r8), pointer :: soilpsi(:,:)         ! soil water potential in each
!
! !LOCAL VARIABLES:
!EOP
   integer :: g,l,c,p,j    ! indices
   integer :: begp, endp   ! per-clump/proc beginning and ending pft indices
   integer :: begc, endc   ! per-clump/proc beginning and ending column indices
   integer :: begl, endl   ! per-clump/proc beginning and ending landunit indices
   integer :: begg, endg   ! per-clump/proc gridcell ending gridcell indices
!-----------------------------------------------------------------------

    ! assign local pointers at the gridcell level

    ! assign local pointers at the landunit level

    ! assign local pointers at the column level
    soilpsi                        => clm3%g%l%c%cps%soilpsi

    ! assign local pointers at the pft level
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

   ! Determine subgrid bounds on this processor
   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

   ! initialize column-level variables

   ! initialize pft-level variables
   do p = begp, endp
      casa_agnpp(p) = 0._r8
      casa_ar(p) = 0._r8
      casa_bgnpp(p) = 0._r8
      casa_cwdc(p) = 0._r8
      casa_cwdc_hr(p) = 0._r8
      casa_cwdc_loss(p) = 0._r8
      casa_frootc(p) = 0._r8
      casa_frootc_alloc(p) = 0._r8
      casa_frootc_loss(p) = 0._r8
      casa_gpp(p) = 0._r8
      casa_hr(p) = 0._r8
      casa_leafc(p) = 0._r8
      casa_leafc_alloc(p) = 0._r8
      casa_leafc_loss(p) = 0._r8
      casa_litterc(p) = 0._r8
      casa_litterc_hr(p) = 0._r8
      casa_litterc_loss(p) = 0._r8
      casa_nee(p) = 0._r8
      casa_nep(p) = 0._r8
      casa_npp(p) = 0._r8
      casa_soilc(p) = 0._r8
      casa_soilc_hr(p) = 0._r8
      casa_soilc_loss(p) = 0._r8
      casa_woodc(p) = 0._r8
      casa_woodc_alloc(p) = 0._r8
      casa_woodc_loss(p) = 0._r8
   end do   ! end of loop over pfts  

end subroutine CASAiniTimeVar

#endif

end module CASAiniTimeVarMod
