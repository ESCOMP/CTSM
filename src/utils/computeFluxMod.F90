module computeFluxMod

  ! provides interface for the flux modules

  use shr_kind_mod      , only : r8 => shr_kind_r8
  implicit none

  type, abstract :: computeFlux_type
     ! flux_name is just used to generate more helpful error messages. This needs to be
     ! set in initialization. (If that's awkward to accomplish, then we could have a
     ! deferred procedure getFluxName that simply returns a character string giving the
     ! name of this flux.)
     character(len=128) :: flux_name
   contains
     procedure(getFlux_interface), deferred :: getFlux
  end type computeFlux_type

  abstract interface
     subroutine getFlux_interface(this,x,f,dfdx)
       use shr_kind_mod      , only : r8 => shr_kind_r8
       import :: computeFlux_type
       implicit none
       class(computeFlux_type), intent(in) :: this
       real(r8), intent(in)            :: x    ! state
       real(r8), intent(out)           :: f    ! f(x)
       real(r8), intent(out), optional :: dfdx ! derivative
     end subroutine getFlux_interface
  end interface

end module computeFluxMod
