module implicitEulerMod

 USE nrtype
 USE computeFluxMod
 implicit none

 contains

 subroutine implicitEuler(funcName, dt, xInit, xNew, err, message)
 ! dummy variables
 procedure(fluxTemplate), pointer :: funcName
 real(dp), intent(in)     :: dt            ! time step
 real(dp), intent(in)     :: xInit         ! initial state
 real(dp), intent(out)    :: xNew          ! updated state
 integer(i4b),intent(out) :: err           ! error code
 character(*),intent(out) :: message       ! error message
 ! local variables
 integer(i4b)             :: iter          ! iteration index
 integer(i4b),parameter   :: maxiter=20    ! maximum number of iterations
 real(dp), parameter      :: xTol=1.e-8_dp ! convergence tolerance
 real(dp)                 :: flux          ! flux
 real(dp)                 :: dfdx          ! derivative
 real(dp)                 :: xRes          ! residual error
 real(dp)                 :: delX          ! state update
 ! initialiuze routine
 err=0; message='implicitEuler/'

 ! initialize
 xNew = xInit

 ! iterate
 do iter=1,maxiter
  ! get the flux and derivative
  call getFlux(funcName, xNew, flux, dfdx)
  ! get the residual and the update
  xRes = xNew - (xInit + dt*flux)
  delX = -xRes/(1._dp - dt*dfdx)
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
