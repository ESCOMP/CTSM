module MLWaterVaporMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate saturation vapor pressure and latent heat of vaporization
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
  public :: SatVap     ! Saturation vapor pressure and derivative
  public :: LatVap     ! Latent heat of vaporization
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SatVap (t, es, desdt)
    !
    ! !DESCRIPTION:
    ! Compute saturation vapor pressure and change in saturation vapor pressure
    ! with respect to temperature. Polynomial approximations are from:
    ! Flatau et al (1992) Polynomial fits to saturation vapor pressure.
    ! Journal of Applied Meteorology 31:1507-1513
    !
    ! !USES:
    use clm_varcon, only : tfrz
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: t        ! Temperature (K)
    real(r8), intent(out) :: es       ! Vapor pressure (Pa)
    real(r8), intent(out) :: desdt    ! d(es)/d(t) (Pa/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: tc                    ! Temperature (C)
    !---------------------------------------------------------------------

    ! For water vapor (temperature range is 0C to 100C)
 
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
 
    ! and for derivative
 
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
 
    ! For ice (temperature range is -75C to 0C)
 
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
 
    ! and for derivative
 
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8

    tc = t - tfrz
    if (tc > 100.0_r8) tc = 100.0_r8
    if (tc < -75.0_r8) tc = -75.0_r8

    if (tc >= 0.0_r8) then
       es    = a0 + tc*(a1 + tc*(a2 + tc*(a3 + tc*(a4 &
             + tc*(a5 + tc*(a6 + tc*(a7 + tc*a8)))))))
       desdt = b0 + tc*(b1 + tc*(b2 + tc*(b3 + tc*(b4 &
             + tc*(b5 + tc*(b6 + tc*(b7 + tc*b8)))))))
    else
       es    = c0 + tc*(c1 + tc*(c2 + tc*(c3 + tc*(c4 &
             + tc*(c5 + tc*(c6 + tc*(c7 + tc*c8)))))))
       desdt = d0 + tc*(d1 + tc*(d2 + tc*(d3 + tc*(d4 &
             + tc*(d5 + tc*(d6 + tc*(d7 + tc*d8)))))))
    end if

    es    = es    * 100._r8            ! Convert from mb to Pa
    desdt = desdt * 100._r8            ! Convert from mb to Pa

  end subroutine SatVap

  !-----------------------------------------------------------------------
  function LatVap (t) result(lambda)
    !
    ! !DESCRIPTION:
    ! Latent heat of vaporization in relation to air temperature, as in CLM
    !
    ! !USES:
    use clm_varcon, only : tfrz, hvap, hsub
    use MLclm_varcon, only : mmh2o
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: t     ! Temperature (K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lambda             ! Molar latent heat of vaporization (J/mol)
    !---------------------------------------------------------------------

    if (t > tfrz) then
       lambda = hvap
    else
       lambda = hsub
    end if
    lambda = lambda * mmh2o

  end function LatVap

end module MLWaterVaporMod
