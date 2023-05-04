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
    use MLclm_varctl, only : acclim_type, kn_val
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp              ! Filter index
    integer  :: p               ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic              ! Aboveground layer index
    real(r8) :: jmax25_to_vcmax25 ! Ratio of jmax to vcmax at 25C (umol/umol)
    real(r8) :: vcmax25top      ! Canopy top - Maximum carboxylation rate at 25C (umol/m2/s)
    real(r8) :: jmax25top       ! Canopy top - C3: Maximum electron transport rate at 25C (umol/m2/s)
    real(r8) :: rd25top         ! Canopy top - Leaf respiration rate at 25C (umol CO2/m2/s)
    real(r8) :: kp25top         ! Canopy top - C4: Initial slope of CO2 response curve at 25C (mol/m2/s)
    real(r8) :: kn              ! Leaf nitrogen decay coefficient
    real(r8) :: pai_above       ! Cumulative plant area index above canopy layer
    real(r8) :: nscale          ! Nitrogen scaling coefficient
    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    c3psn     => pftcon%c3psn                  , &  ! CLM: Photosynthetic pathway (1. = C3 plant, 0. = C4 plant)
    vcmaxpft  => pftcon%vcmaxpft               , &  ! CLM (new): Maximum carboxylation rate at 25C (umol/m2/s)
    tacclim   => mlcanopy_inst%tacclim_forcing , &  ! Average air temperature for acclimation (K)
    ncan      => mlcanopy_inst%ncan_canopy     , &  ! Number of layers
    dpai      => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
                                                    ! *** Output ***
    vcmax25   => mlcanopy_inst%vcmax25_profile , &  ! Canopy layer leaf maximum carboxylation rate at 25C (umol/m2/s)
    jmax25    => mlcanopy_inst%jmax25_profile  , &  ! Canopy layer C3 maximum electron transport rate at 25C (umol/m2/s)
    rd25      => mlcanopy_inst%rd25_profile    , &  ! Canopy layer leaf respiration rate at 25C (umol CO2/m2/s)
    kp25      => mlcanopy_inst%kp25_profile      &  ! Canopy layer C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Vcmax and other parameters (at 25C and top of canopy)

       select case (acclim_type)
       case (0)
          jmax25_to_vcmax25 = jmax25_to_vcmax25_noacclim
       case (1)
          jmax25_to_vcmax25_acclim = 2.59_r8 - 0.035_r8*min(max((tacclim(p)-tfrz),11._r8),35._r8)
          jmax25_to_vcmax25 = jmax25_to_vcmax25_acclim
       case default
          call endrun (msg=' ERROR: CanopyNitrogenProfile: acclim_type not valid')             
       end select

       vcmax25top = vcmaxpft(patch%itype(p))

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
       ! (1) calculated from Vcmax or (2) a specified value

       if (kn_val <= 0._r8) then
          kn = exp(0.00963_r8 * vcmax25top - 2.43_r8)
       else if (kn_val > 0._r8) then
          kn = kn_val
       end if

       ! Layer values - nscale is integrated over the layer to get the average
       ! value for the layer

       pai_above = 0._r8
       do ic = ncan(p), 1, -1
          if (dpai(p,ic) > 0._r8) then
             nscale = exp(-kn * pai_above) / (dpai(p,ic) * kn) * (1._r8 - exp(-kn * dpai(p,ic)))
             vcmax25(p,ic) = vcmax25top * nscale
             jmax25(p,ic) = jmax25top * nscale
             rd25(p,ic) = rd25top * nscale
             kp25(p,ic) = kp25top * nscale
          else
             vcmax25(p,ic) = 0._r8
             jmax25(p,ic) = 0._r8
             rd25(p,ic) = 0._r8
             kp25(p,ic) = 0._r8
          end if
          pai_above = pai_above + dpai(p,ic)
       end do

    end do

    end associate
  end subroutine CanopyNitrogenProfile

end module MLCanopyNitrogenProfileMod
