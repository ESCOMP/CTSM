module SoiWatRetCurveParMod

  !-----------------------------------------------------------------------------
  !DESCRIPTIONS
  !module contains functions to compute soil water retention curve
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:   
  public  :: soil_hk
  public  :: soil_suction
  public :: init_soil_WatRetCurve
  !
  ! !PRIVATE DATA MEMBERS:
  !the two parameters below will be set through namelist
  integer :: soil_suction_method         
  integer :: soil_hk_method
  integer, parameter :: soil_suction_clapphornberg_1978=0
  integer, parameter :: soil_hk_clapphornberg_1978=0
  !-----------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------------
  subroutine init_soil_WatRetCurve()
    !
    !DESCRIPTIONS
    !Initialize the methods for soil water retention curve calculation
    implicit none

    soil_suction_method = soil_suction_clapphornberg_1978         
    soil_hk_method = soil_hk_clapphornberg_1978

  end subroutine init_soil_WatRetCurve

  !-----------------------------------------------------------------------------
  subroutine soil_suction(smpsat, s, bsw, smp, dsmpds) 
    !
    !DESCRIPTION
    ! compute soil suction potential
    !
    ! !USES
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use abortutils    , only : endrun   
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)            :: smpsat   !minimum soil suction, positive [mm]
    real(r8), intent(in)            :: s        !reletive saturation, [0, 1]
    real(r8), intent(in)            :: bsw      !shape parameter
    real(r8), intent(out)           :: smp      !soil suction, negative, [mm]
    real(r8), optional, intent(out) :: dsmpds   !d[smp]/ds, [mm]
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'soil_suction'  ! subroutine name
    !------------------------------------------------------------------------------

    select case (soil_suction_method)

    case (soil_suction_clapphornberg_1978)
       if(present(dsmpds))then
          call soil_suction_clapphornberg1978(smpsat, s, bsw, smp, dsmpds)
       else
          call soil_suction_clapphornberg1978(smpsat, s, bsw, smp)      
       endif
    case default
       call endrun(subname // ':: a soil suction function must be specified!')       
    end select

  end subroutine soil_suction

  !-----------------------------------------------------------------------------
  subroutine soil_hk(hksat, imped, s, bsw, hk, dhkds)   
    !
    !DESCRIPTION
    ! compute soil suction potential
    !   
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use abortutils    , only : endrun      
    implicit none
    real(r8), intent(in) :: hksat    !saturated hydraulic conductivity [mm/s]
    real(r8), intent(in) :: imped    !ice impedance
    real(r8), intent(in) :: s        !reletive saturation, [0, 1]
    real(r8), intent(in) :: bsw      !shape parameter

    real(r8), intent(out):: hk       !hydraulic conductivity [mm/s]
    real(r8), optional, intent(out):: dhkds    !d[hk]/ds   [mm/s]

    character(len=32) :: subname = 'soil_hk'  ! subroutine name

    select case (soil_hk_method)

    case (soil_hk_clapphornberg_1978)
       if(present(dhkds))then
          call soil_hk_clappHornberg1978(hksat, imped, s, bsw, hk, dhkds)
       else
          call soil_hk_clappHornberg1978(hksat, imped, s, bsw, hk)
       endif
    case default
       call endrun(subname // ':: a soil hk function must be specified!')        
    end select

  end subroutine soil_hk

  !-----------------------------------------------------------------------------
  subroutine soil_suction_clappHornberg1978(smpsat, s, bsw, smp, dsmpds)
    !
    !DESCRIPTION
    ! compute the suction pressure using the Clapp-Hornberg parameterization
    ! Note: the effect of grdual air entry is ignored and the parameterization
    ! is isothermal
    !
    ! !USES
    use shr_kind_mod  , only : r8 => shr_kind_r8    
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: smpsat   !minimum soil suction, positive [mm]
    real(r8), intent(in) :: s        !reletive saturation, [0, 1]
    real(r8), intent(in) :: bsw      !shape parameter
    real(r8), intent(out):: smp      !soil suction, negative, [mm]
    real(r8), optional, intent(out):: dsmpds   !d[smp]/ds, [mm]
    !------------------------------------------------------------------------------

    !compute soil suction potential, negative
    smp = -smpsat*s**(-bsw)

    !compute derivative
    if(present(dsmpds))then
       dsmpds=-bsw*smp/s
    endif
  end subroutine soil_suction_clapphornberg1978

  !-----------------------------------------------------------------------------
  subroutine soil_hk_clapphornberg1978(hksat, imped, s, bsw, hk, dhkds)
    !
    ! DESCRIPTIONS
    ! compute hydraulic conductivity using the Clapp-Hornberg parameterization
    ! the ice induced impedance will be parameterized outside the
    ! isothermal computation
    !
    ! !USES
    use shr_kind_mod  , only : r8 => shr_kind_r8   
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)            :: hksat    !saturated hydraulic conductivity [mm/s]
    real(r8), intent(in)            :: imped    !ice impedance
    real(r8), intent(in)            :: s        !reletive saturation, [0, 1]
    real(r8), intent(in)            :: bsw      !shape parameter
    real(r8), intent(out)           :: hk       !hydraulic conductivity [mm/s]
    real(r8), optional, intent(out) :: dhkds    !d[hk]/ds   [mm/s]
    !------------------------------------------------------------------------------

    !compute hydraulic conductivity
    hk=imped*hksat*s**(2._r8*bsw+3._r8)

    !compute the derivative
    if(present(dhkds))then
       dhkds=(2._r8*bsw+3._r8)*hk/s
    endif

  end subroutine soil_hk_clapphornberg1978

end module SoiWatRetCurveParMod
