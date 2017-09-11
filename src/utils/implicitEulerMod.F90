module implicitEulerMod

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use computeFluxMod    , only : computeFlux_type
  use abortutils        , only : endrun
  implicit none

contains

  subroutine implicitEuler(computeFlux_inst, dt, xInit, flux)
    ! dummy variables
    class(computeFlux_type), intent(in) :: computeFlux_inst
    real(r8), intent(in)     :: dt            ! time step
    real(r8), intent(in)     :: xInit         ! initial state
    real(r8), intent(out)    :: flux          ! flux
    ! local variables
    integer             :: iter          ! iteration index
    integer,parameter   :: maxiter=20    ! maximum number of iterations
    real(r8), parameter      :: xTol=1.e-8_r8 ! convergence tolerance
    real(r8)                 :: xNew          ! updated state
    real(r8)                 :: dfdx          ! derivative
    real(r8)                 :: xRes          ! residual error
    real(r8)                 :: delX          ! state update

    ! initialize
    xNew = xInit

    ! iterate
    do iter=1,maxiter
       ! get the flux and derivative
       call computeFlux_inst%getFlux(xNew, flux, dfdx)
       ! get the residual and the update
       xRes = xNew - (xInit + dt*flux)
       delX = -xRes/(1._r8 - dt*dfdx)
       ! update and check for convergence
       xNew = xNew + delX
       if(abs(delX) < xTol) exit
       ! check for non-convergence
       if(iter==maxiter)then
          call endrun(computeFlux_inst%flux_name // ': the implicit Euler solution did not converge!')
       endif
    end do

    ! Get final flux
    call computeFlux_inst%getFlux(xNew, flux)
  end subroutine implicitEuler

end module implicitEulerMod
