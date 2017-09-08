module computeFluxMod

  ! provides interface for the flux modules

  use shr_kind_mod      , only : r8 => shr_kind_r8
  implicit none

  type, abstract :: flux_type
   contains
     procedure(getFlux_interface), deferred :: getFlux
  end type flux_type

  abstract interface
     subroutine getFlux_interface(this,x,f,dfdx)
       use shr_kind_mod      , only : r8 => shr_kind_r8
       import :: flux_type
       implicit none
       class(flux_type), intent(in) :: this
       real(r8), intent(in)            :: x    ! state
       real(r8), intent(out)           :: f    ! f(x)
       real(r8), intent(out), optional :: dfdx ! derivative
     end subroutine getFlux_interface
  end interface

end module computeFluxMod
