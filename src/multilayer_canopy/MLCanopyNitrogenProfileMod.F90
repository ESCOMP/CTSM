module MLCanopyNitrogenProfileMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Canopy profile of nitrogen and photosynthetic capacity
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyNitrogenProfile
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine CanopyNitrogenProfile (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy profile of nitrogen and photosynthetic capacity
    !
    ! !USES:
    use clm_varcon, only : tfrz
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : jmax25_to_vcmax25_noacclim, jmax25_to_vcmax25_acclim
    use MLclm_varcon, only : rd25_to_vcmax25_c3, rd25_to_vcmax25_c4, kp25_to_vcmax25_c4
    use MLclm_varctl, only : acclim_type, kn_val, leaf_optics_type
    use MLclm_varpar, only : isun, isha
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                      ! Aboveground layer index
    real(r8) :: jmax25_to_vcmax25       ! Ratio of jmax to vcmax at 25C (umol/umol)
    real(r8) :: vcmax25top              ! Canopy top - Maximum carboxylation rate at 25C (umol/m2/s)
    real(r8) :: jmax25top               ! Canopy top - C3: Maximum electron transport rate at 25C (umol/m2/s)
    real(r8) :: rd25top                 ! Canopy top - Leaf respiration rate at 25C (umol CO2/m2/s)
    real(r8) :: kp25top                 ! Canopy top - C4: Initial slope of CO2 response curve at 25C (mol/m2/s)
    real(r8) :: kn                      ! Leaf nitrogen decay coefficient
    real(r8) :: pai_above               ! Cumulative plant area index above canopy layer
    real(r8) :: fn, fn_sun, fn_sha      ! Nitrogen factor integrated over canopy layer and sun/shade components
    real(r8) :: nscale_sun, nscale_sha  ! Nitrogen scaling coefficient for sun/shade leaves
    real(r8) :: numerical, analytical   ! Numerical and analytical values for vcmax25 canopy integration
    !---------------------------------------------------------------------

    associate ( &
                                                          ! *** Input ***
    c3psn           => pftcon%c3psn                  , &  ! CLM: Photosynthetic pathway (1. = C3 plant, 0. = C4 plant)
    vcmaxpft        => pftcon%vcmaxpft               , &  ! CLMml: Maximum carboxylation rate at 25C (umol/m2/s)
    clump_fac       => pftcon%clump_fac              , &  ! CLMml: Foliage clumping index (-)
    tacclim         => mlcanopy_inst%tacclim_forcing , &  ! Average air temperature for acclimation (K)
    ncan            => mlcanopy_inst%ncan_canopy     , &  ! Number of aboveground layers
    lai             => mlcanopy_inst%lai_canopy      , &  ! Leaf area index of canopy (m2/m2)
    sai             => mlcanopy_inst%sai_canopy      , &  ! Stem area index of canopy (m2/m2)
    dpai            => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    fracsun         => mlcanopy_inst%fracsun_profile , &  ! Canopy layer sunlit fraction (-)
    kb              => mlcanopy_inst%kb_profile      , &  ! Direct beam extinction coefficient (-)
    tbi             => mlcanopy_inst%tbi_profile     , &  ! Cumulative transmittance of direct beam onto canopy layer (-)
                                                          ! *** Output ***
    vcmax25_leaf    => mlcanopy_inst%vcmax25_leaf    , &  ! Leaf maximum carboxylation rate at 25C (umol/m2/s)
    jmax25_leaf     => mlcanopy_inst%jmax25_leaf     , &  ! Leaf C3 maximum electron transport rate at 25C (umol/m2/s)
    rd25_leaf       => mlcanopy_inst%rd25_leaf       , &  ! Leaf respiration rate at 25C (umol CO2/m2/s)
    kp25_leaf       => mlcanopy_inst%kp25_leaf       , &  ! Leaf C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    vcmax25_profile => mlcanopy_inst%vcmax25_profile , &  ! Canopy layer leaf maximum carboxylation rate at 25C (umol/m2/s)
    jmax25_profile  => mlcanopy_inst%jmax25_profile  , &  ! Canopy layer C3 maximum electron transport rate at 25C (umol/m2/s)
    rd25_profile    => mlcanopy_inst%rd25_profile    , &  ! Canopy layer leaf respiration rate at 25C (umol CO2/m2/s)
    kp25_profile    => mlcanopy_inst%kp25_profile      &  ! Canopy layer C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Vcmax and other parameters (at 25C and top of canopy)

       vcmax25top = vcmaxpft(patch%itype(p))

       select case (acclim_type)
       case (0)
          jmax25_to_vcmax25 = jmax25_to_vcmax25_noacclim
       case (1)
          jmax25_to_vcmax25_acclim = 2.59_r8 - 0.035_r8*min(max((tacclim(p)-tfrz),11._r8),35._r8)
          jmax25_to_vcmax25 = jmax25_to_vcmax25_acclim
       case default
          call endrun (msg=' ERROR: CanopyNitrogenProfile: acclim_type not valid')             
       end select

       if (nint(c3psn(patch%itype(p))) == 1) then
          jmax25top = jmax25_to_vcmax25 * vcmax25top
          rd25top = rd25_to_vcmax25_c3 * vcmax25top
          kp25top = 0._r8
       else if (nint(c3psn(patch%itype(p))) == 0) then
          jmax25top = 0._r8
          rd25top = rd25_to_vcmax25_c4 * vcmax25top
          kp25top = kp25_to_vcmax25_c4 * vcmax25top
       end if

       ! Leaf nitrogen decay coefficient based on one of two methods:
       ! (1) calculated from Vcmax (kn_val < 0), or 
       ! (2) a specified value (kn_val > 0)

       if (kn_val < 0._r8) then
          kn = exp(0.00963_r8 * vcmax25top - 2.43_r8)
       else if (kn_val > 0._r8) then
          kn = kn_val
       else
          call endrun (msg='ERROR: CanopyNitrogenProfile: incorrect Kn')
       end if

       ! Layer values - nscale is integrated over the layer to get the average
       ! value for the sunlit and shaded portions of the layer

       pai_above = 0._r8
       do ic = ncan(p), 1, -1

          ! Initialize all layers to zero

          vcmax25_leaf(p,ic,isun) = 0._r8 ; vcmax25_leaf(p,ic,isha) = 0._r8
          jmax25_leaf(p,ic,isun) = 0._r8  ; jmax25_leaf(p,ic,isha) = 0._r8
          rd25_leaf(p,ic,isun) = 0._r8    ; rd25_leaf(p,ic,isha) = 0._r8
          kp25_leaf(p,ic,isun) = 0._r8    ; kp25_leaf(p,ic,isha) = 0._r8

          ! Sunlit and shaded leaves

          if (dpai(p,ic) > 0._r8) then

             ! Nitrogen factor - two different cases:
             !
             ! leaf_optics_type = 0 integrates leaf nitrogen over sunlit and
             ! shaded portions of the canopy. It is equivalent to a one-layer
             ! (big-leaf) canopy, but requires that kb and clump_fac are
             ! invariant with height.
             !
             ! leaf_optics_type = 1 allows kb and clump_fac to vary with height.
             ! It sets sunlit and shaded leaves to have the same leaf nitrogen.
             ! This is true for thin layers (dpai <= 0.1), is approximately true
             ! for dpai <= 1, and is not true for large dpai. Consequently, it
             ! does not give the same answer when the canopy is reduced to one layer.

             fn = exp(-kn * pai_above) * (1._r8 - exp(-kn * dpai(p,ic))) / kn
             select case (leaf_optics_type)
             case (0)
                fn_sun = clump_fac(patch%itype(p)) / (kn + kb(p,ic) * clump_fac(patch%itype(p))) &
                       * exp(-kn * pai_above) * tbi(p,ic) &
                       * (1._r8 - exp(-(kn + kb(p,ic)*clump_fac(patch%itype(p))) * dpai(p,ic)))
                fn_sha = fn - fn_sun
                nscale_sun = fn_sun / (fracsun(p,ic) * dpai(p,ic))
                nscale_sha = fn_sha / ((1._r8 - fracsun(p,ic)) * dpai(p,ic))
             case (1)
                nscale_sun = fn / dpai(p,ic)
                nscale_sha = nscale_sun
             end select

             ! Leaf variables

             vcmax25_leaf(p,ic,isun) = vcmax25top * nscale_sun ; vcmax25_leaf(p,ic,isha) = vcmax25top * nscale_sha
             jmax25_leaf(p,ic,isun)  = jmax25top  * nscale_sun ; jmax25_leaf(p,ic,isha)  = jmax25top  * nscale_sha
             rd25_leaf(p,ic,isun)    = rd25top    * nscale_sun ; rd25_leaf(p,ic,isha)    = rd25top    * nscale_sha
             kp25_leaf(p,ic,isun)    = kp25top    * nscale_sun ; kp25_leaf(p,ic,isha)    = kp25top    * nscale_sha

          end if

          ! Increment cumulative plant area index above the next layer

          pai_above = pai_above + dpai(p,ic)

          ! Layer weighted mean

          vcmax25_profile(p,ic) = vcmax25_leaf(p,ic,isun) * fracsun(p,ic) + vcmax25_leaf(p,ic,isha) * (1._r8-fracsun(p,ic))
          jmax25_profile(p,ic)  = jmax25_leaf(p,ic,isun)  * fracsun(p,ic) + jmax25_leaf(p,ic,isha)  * (1._r8-fracsun(p,ic))
          rd25_profile(p,ic)    = rd25_leaf(p,ic,isun)    * fracsun(p,ic) + rd25_leaf(p,ic,isha)    * (1._r8-fracsun(p,ic))
          kp25_profile(p,ic)    = kp25_leaf(p,ic,isun)    * fracsun(p,ic) + kp25_leaf(p,ic,isha)    * (1._r8-fracsun(p,ic))

       end do

       ! Check that canopy sum of vcmax25 equals the expected value obtained by analytical integration
       ! NB slevis 2024/3/22: The next was a valid error check before introducing minimum dlai and dsai

       numerical = sum(vcmax25_profile(p,1:ncan(p)) * dpai(p,1:ncan(p)))
       analytical = vcmax25top * (1._r8 - exp(-kn*(lai(p) + sai(p)))) / kn
!      if (abs(numerical-analytical) > 1.e-06_r8) then
!         call endrun (msg='ERROR: CanopyNitrogenProfile: canopy integration error')
!      end if

    end do

    end associate
  end subroutine CanopyNitrogenProfile

end module MLCanopyNitrogenProfileMod
