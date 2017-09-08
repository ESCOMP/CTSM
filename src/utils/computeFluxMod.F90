module computeFluxMod

  ! provides interface for the flux modules

  use shr_kind_mod      , only : r8 => shr_kind_r8
  implicit none

  abstract interface
     subroutine fluxTemplate(x,f,dfdx)
       use shr_kind_mod      , only : r8 => shr_kind_r8
       implicit none
       real(r8), intent(in)            :: x    ! state
       real(r8), intent(out)           :: f    ! f(x)
       real(r8), intent(out), optional :: dfdx ! derivative
     end subroutine fluxTemplate
  end interface

contains

  subroutine getFlux(func, x, f, dfdx)
    procedure(fluxTemplate), pointer :: func
    real(r8), intent(in)           :: x    ! state
    real(r8), intent(out)          :: f    ! f(x)
    real(r8), intent(out),optional :: dfdx ! derivative
    call func(x,f,dfdx)
  end subroutine getFlux

end module computeFluxMod
