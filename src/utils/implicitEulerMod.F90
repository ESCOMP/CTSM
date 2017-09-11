module implicitEulerMod

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use computeFluxMod    , only : computeFlux_type
  implicit none

contains

  subroutine implicitEuler(flux_inst, dt, xInit, xNew, err, message)
    ! dummy variables
    class(computeFlux_type), intent(in) :: flux_inst
    real(r8), intent(in)     :: dt            ! time step
    real(r8), intent(in)     :: xInit         ! initial state
    real(r8), intent(out)    :: xNew          ! updated state
    integer,intent(out) :: err           ! error code
    character(*),intent(out) :: message       ! error message
    ! local variables
    integer             :: iter          ! iteration index
    integer,parameter   :: maxiter=20    ! maximum number of iterations
    real(r8), parameter      :: xTol=1.e-8_r8 ! convergence tolerance
    real(r8)                 :: flux          ! flux
    real(r8)                 :: dfdx          ! derivative
    real(r8)                 :: xRes          ! residual error
    real(r8)                 :: delX          ! state update
    ! initialiuze routine
    err=0; message='implicitEuler/'

    ! initialize
    xNew = xInit

    ! iterate
    do iter=1,maxiter
       ! get the flux and derivative
       call flux_inst%getFlux(xNew, flux, dfdx)
       ! get the residual and the update
       xRes = xNew - (xInit + dt*flux)
       delX = -xRes/(1._r8 - dt*dfdx)
       ! update and check for convergence
       xNew = xNew + delX
       if(abs(delX) < xTol) exit
       ! check for non-convergence
       if(iter==maxiter)then
          message=trim(message)//'the implicit Euler solution did not converge!'
          err=20; return
       endif
    end do
  end subroutine implicitEuler

end module implicitEulerMod
