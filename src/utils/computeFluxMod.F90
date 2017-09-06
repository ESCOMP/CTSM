module computeFluxMod

 ! provides interface for the flux modules

 USE nrtype
 implicit none

 abstract interface
 subroutine fluxTemplate(x,f,dfdx)
  USE nrtype
  implicit none
  real(dp), intent(in)            :: x    ! state
  real(dp), intent(out)           :: f    ! f(x)
  real(dp), intent(out), optional :: dfdx ! derivative
 end subroutine fluxTemplate
 end interface

 contains

  subroutine getFlux(func, x, f, dfdx)
   procedure(fluxTemplate), pointer :: func
   real(dp), intent(in)           :: x    ! state
   real(dp), intent(out)          :: f    ! f(x)
   real(dp), intent(out),optional :: dfdx ! derivative
   call func(x,f,dfdx)
  end subroutine getFlux

end module computeFluxMod
