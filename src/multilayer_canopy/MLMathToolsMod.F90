module MLMathToolsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Math tools
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  use MLCanopyFluxesType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: hybrid             ! Solve for the root of a function using secant and Brent's methods
  public :: zbrent             ! Use Brent's method to find the root of a function
  public :: quadratic          ! Solve a quadratic equation for its two roots
  public :: tridiag            ! Solve a tridiagonal system of equations
  public :: beta_function      ! Evaluate the beta function at p and q: B(p,q)
  public :: log_gamma_function ! Evaluate the log natural of the gamma function at x: ln(G(x))

  interface
    subroutine func (p, ic, il, mlcanopy_inst, x, val)
    use shr_kind_mod, only : r8 => shr_kind_r8
    use MLCanopyFluxesType, only : mlcanopy_type
    integer, intent(in) :: p, ic, il
    real(r8), intent(in) :: x
    real(r8), intent(out) :: val
    type(mlcanopy_type) :: mlcanopy_inst
    end subroutine func
  end interface

contains

  !-----------------------------------------------------------------------
  function hybrid (msg, p, ic, il, mlcanopy_inst, func, xa, xb, tol) result(root)
    !
    ! !DESCRIPTION:
    ! Solve for the root of a function, given initial estimates xa and xb.
    ! The root is updated until its accuracy is tol.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*) :: msg           ! String to be printed
    integer, intent(in) :: p          ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic         ! Canopy layer index
    integer, intent(in) :: il         ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: xa, xb    ! Initial estimates of root
    real(r8), intent(in) :: tol       ! Error tolerance
    external :: func                  ! Function to solve
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x0, x1                ! Estimates of root
    real(r8) :: f0, f1                ! Function value for x0 and x1
    real(r8) :: minx                  ! x0 or x1 that gives smallest function value
    real(r8) :: minf                  ! Smallest function value obtained using x0 or x1
    real(r8) :: root                  ! Value for root
    real(r8) :: dx                    ! Change in root
    real(r8) :: x                     ! Updated root
    integer :: iter                   ! Iteration loop index
    integer, parameter :: itmax = 40  ! Maximum number of iterations
    !---------------------------------------------------------------------

    x0 = xa
    call func (p, ic, il, mlcanopy_inst, x0, f0)
    if (f0 == 0._r8) then
       root = x0
       return
    end if

    x1 = xb
    call func (p, ic, il, mlcanopy_inst, x1, f1)
    if (f1 == 0._r8) then
       root = x1
       return
    end if

    if (f1 < f0) then
       minx = x1
       minf = f1
    else
       minx = x0
       minf = f0
    end if

    ! First use the secant method, and then use the brent method as a backup

    iter = 0
    do
       iter = iter + 1
       dx = -f1 * (x1 - x0) / (f1 - f0)
       x = x1 + dx
       if (abs(dx) < tol) then
          x0 = x
          exit
       end if
       x0 = x1
       f0 = f1
       x1 = x
       call func (p, ic, il, mlcanopy_inst, x1, f1)
       if (f1 < minf) then
          minx = x1
          minf = f1
       end if

       ! If a root zone is found, use the brent method for a robust backup strategy

       if (f1 * f0 < 0._r8) then
          x = zbrent (msg, p, ic, il, mlcanopy_inst, func, x0, x1, tol)
          x0 = x
          exit
       end if

       ! In case of failing to converge within itmax iterations stop at the minimum function

       if (iter > itmax) then
          call func (p, ic, il, mlcanopy_inst, minx, f1)
          x0 = minx
          exit
       end if

    end do

    root = x0

  end function hybrid

  !-----------------------------------------------------------------------
  function zbrent (msg, p, ic, il, mlcanopy_inst, func, xa, xb, tol) result(root)
    !
    ! !DESCRIPTION:
    ! Use Brent's method to find the root of a function, which is known to exist
    ! between xa and xb. The root is updated until its accuracy is tol.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*) :: msg           ! String to be printed
    integer, intent(in) :: p          ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic         ! Canopy layer index
    integer, intent(in) :: il         ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: xa, xb    ! Minimum and maximum of the variable domain to search
    real(r8), intent(in) :: tol       ! Error tolerance
    external :: func                  ! Function to solve
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer, parameter :: itmax = 50          ! Maximum number of iterations
    real(r8), parameter :: eps = 1.e-08_r8    ! Relative error tolerance
    integer  :: iter                          ! Iteration loop index
    real(r8) :: a,b,c,d,e,fa,fb,fc,pp,q,r,s,tol1,xm
    real(r8) :: root
    !---------------------------------------------------------------------

    a = xa
    b = xb
    call func (p, ic, il, mlcanopy_inst, a, fa)
    call func (p, ic, il, mlcanopy_inst, b, fb)

    if ((fa > 0._r8 .and. fb > 0._r8) .or. (fa < 0._r8 .and. fb < 0._r8)) then
       write (iulog,*) 'zbrent: Root must be bracketed'
       write (iulog,*) 'called from: ',msg
       write (iulog,*) xa, fa
       write (iulog,*) xb, fb
       call endrun (msg=' ERROR: zbrent error')
    end if
    c = b
    fc = fb
    iter = 0
    do
       if (iter == itmax) exit
       iter = iter + 1
       if ((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8)) then
          c = a
          fc = fa
          d = b - a
          e = d
       end if
       if (abs(fc) < abs(fb)) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
       end if
       tol1 = 2._r8 * eps * abs(b) + 0.5_r8 * tol
       xm = 0.5_r8 * (c - b)
       if (abs(xm) <= tol1 .or. fb == 0._r8) exit
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s = fb / fa
          if (a == c) then
             pp = 2._r8 * xm * s
             q = 1._r8 - s
          else
             q = fa / fc
             r = fb / fc
             pp = s * (2._r8 * xm * q * (q - r) - (b - a) * (r - 1._r8))
             q = (q - 1._r8) * (r - 1._r8) * (s - 1._r8)
          end if
          if (pp > 0._r8) q = -q
          pp = abs(pp)
          if (2._r8*pp < min(3._r8*xm*q-abs(tol1*q),abs(e*q))) then
             e = d
             d = pp / q
          else
             d = xm
             e = d
          end if
       else
          d = xm
          e = d
       end if
       a = b
       fa = fb
       if (abs(d) > tol1) then
          b = b + d
       else
          b = b + sign(tol1,xm)
       end if
       call func (p, ic, il, mlcanopy_inst, b, fb)
       if (fb == 0._r8) exit
    end do
    root = b

    if (iter == itmax) then
       write (iulog,*) 'zbrent: Maximum number of interations exceeded'
       write (iulog,*) 'called from: ',msg
       call endrun (msg=' ERROR: zbrent error')
    end if

  end function zbrent

  !-----------------------------------------------------------------------
  subroutine quadratic (a, b, c, r1, r2)
    !
    ! !DESCRIPTION:
    ! Solve a quadratic equation for its two roots
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: a,b,c       ! Terms for quadratic equation
    real(r8), intent(out) :: r1,r2       ! Roots of quadratic equation
    !
    ! !LOCAL VARIABLES:
    real(r8) :: q                        ! Temporary term for quadratic solution
    !---------------------------------------------------------------------

    if (a == 0._r8) then
       write (iulog,*) 'Quadratic solution error: a = ',a
       call endrun (msg=' ERROR: quadratic error')
    end if

    if (b >= 0._r8) then
       q = -0.5_r8 * (b + sqrt(b*b - 4._r8*a*c))
    else
       q = -0.5_r8 * (b - sqrt(b*b - 4._r8*a*c))
    end if

    r1 = q / a
    if (q /= 0._r8) then
       r2 = c / q
    else
       r2 = 1.e36_r8
    end if

  end subroutine quadratic

  !-----------------------------------------------------------------------
  subroutine tridiag (a, b, c, r, u, n)
    !
    ! !DESCRIPTION:
    ! Solve a tridiagonal system of equations
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    integer,  intent(in)  :: n        ! Number of layers
    real(r8), intent(in)  :: a(n)     ! A vector for tridiagonal solution
    real(r8), intent(in)  :: b(n)     ! B vector for tridiagonal solution
    real(r8), intent(in)  :: c(n)     ! C vector for tridiagonal solution
    real(r8), intent(in)  :: r(n)     ! R vector for tridiagonal solution
    real(r8), intent(out) :: u(n)     ! U vector for tridiagonal solution
    !
    ! !LOCAL VARIABLES:
    real(r8) :: gam(n)                ! Temporary calculation
    real(r8) :: bet                   ! Temporary calculation
    integer :: j                      ! Layer index
    !---------------------------------------------------------------------

    ! Tridiagonal solution:
    !
    ! Solve for U given the set of equations F x U = R, where U is a vector
    ! of length N, R is a vector of length N, and F is an N x N tridiagonal
    ! matrix defined by the vectors A, B, C (each of length N). A(1) and
    ! C(N) are undefined and are not referenced by the subroutine.
    !
    !    | b(1) c(1)   0  ...                      |   | u(1)   |   | r(1)   |
    !    | a(2) b(2) c(2) ...                      |   | u(2)   |   | r(2)   |
    !    |                ...                      | x | ...    | = | ...    |
    !    |                ... a(n-1) b(n-1) c(n-1) |   | u(n-1) |   | r(n-1) |
    !    |                ...   0    a(n)   b(n)   |   | u(n)   |   | r(n)   |
    !

    bet = b(1)
    u(1) = r(1) / bet
    do j = 2, n
       gam(j) = c(j-1) / bet
       bet = b(j) - a(j) * gam(j)
       u(j) = (r(j) - a(j) * u(j-1)) / bet
    end do
    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1) * u(j+1)
    end do

  end subroutine tridiag

  !-----------------------------------------------------------------------
  subroutine tridiag_2eq (a1, b11, b12, c1, d1, a2, b21, b22, c2, d2, t, q, n)
    !
    ! !DESCRIPTION:
    ! Solve a tridiagonal system of equations with two coupled equations
    ! for air temperature and water vapor at each layer
    !
    ! !USES:
    use MLclm_varpar, only : nlevmlcan
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n                ! Number of layers
    real(r8), intent(in)  :: a1(nlevmlcan)    ! Coefficient for air temperature
    real(r8), intent(in)  :: b11(nlevmlcan)   ! Coefficient for air temperature
    real(r8), intent(in)  :: b12(nlevmlcan)   ! Coefficient for air temperature
    real(r8), intent(in)  :: c1(nlevmlcan)    ! Coefficient for air temperature
    real(r8), intent(in)  :: d1(nlevmlcan)    ! Coefficient for air temperature
    real(r8), intent(in)  :: a2(nlevmlcan)    ! Coefficient for water vapor mole fraction
    real(r8), intent(in)  :: b21(nlevmlcan)   ! Coefficient for water vapor mole fraction
    real(r8), intent(in)  :: b22(nlevmlcan)   ! Coefficient for water vapor mole fraction
    real(r8), intent(in)  :: c2(nlevmlcan)    ! Coefficient for water vapor mole fraction
    real(r8), intent(in)  :: d2(nlevmlcan)    ! Coefficient for water vapor mole fraction
    real(r8), intent(out) :: t(nlevmlcan)     ! Air temperature (K)
    real(r8), intent(out) :: q(nlevmlcan)     ! Water vapor (mol/mol)
    !
    ! !LOCAL VARIABLES:
    integer :: i                              ! Layer index
    real(r8) :: ainv, binv                    ! "a" and "b" elements of 2x2 matrix to invert
    real(r8) :: cinv, dinv                    ! "c" and "d" elements of 2x2 matrix to invert
    real(r8) :: det                           ! Determinant of 2x2 matrix

    real(r8) :: e11(0:nlevmlcan)              ! Coefficient for air temperature
    real(r8) :: e12(0:nlevmlcan)              ! Coefficient for air temperature
    real(r8) :: f1(0:nlevmlcan)               ! Coefficient for air temperature
    real(r8) :: e21(0:nlevmlcan)              ! Coefficient for water vapor
    real(r8) :: e22(0:nlevmlcan)              ! Coefficient for water vapor
    real(r8) :: f2(0:nlevmlcan)               ! Coefficient for water vapor
    !---------------------------------------------------------------------

    ! Solve for air temperature and water vapor (mol/mol):
    !
    ! a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
    ! a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
    !
    ! The solution rewrites these equations so that:
    ! T(i) = f1(i) - e11(i)*T(i+1) - e12(i)*q(i+1) 
    ! q(i) = f2(i) - e21(i)*T(i+1) - e22(i)*q(i+1) 

    e11(0) = 0._r8
    e12(0) = 0._r8
    e21(0) = 0._r8
    e22(0) = 0._r8
    f1(0) = 0._r8
    f2(0) = 0._r8

    do i = 1, n

       ! The matrix to invert is:
       !
       ! B(i) - A(i)*E(i-1)
       !
       ! which is a 2x2 matrix. The
       ! elements in the 2x2 matrix are:
       !
       !                      | a b |
       ! B(i) - A(i)*E(i-1) = | c d |
       !
       ! Calculate these elements (denoted by ainv,binv,
       ! cinv,dinv) and the determinant of the matrix

       ainv = b11(i) - a1(i) * e11(i-1)
       binv = b12(i) - a1(i) * e12(i-1)
       cinv = b21(i) - a2(i) * e21(i-1)
       dinv = b22(i) - a2(i) * e22(i-1)
       det = ainv * dinv - binv * cinv

       ! E(i) = [B(i) - A(i)*E(i-1)]^(-1) * C(i)

       e11(i) =  dinv * c1(i) / det
       e12(i) = -binv * c2(i) / det
       e21(i) = -cinv * c1(i) / det
       e22(i) =  ainv * c2(i) / det

       ! F(i) = [B(i) - A(i)*E(i-1)]^(-1) * [D(i) - A(i)*F(i-1)]

       f1(i) = ( dinv*(d1(i) - a1(i)*f1(i-1)) - binv*(d2(i) - a2(i)*f2(i-1))) / det
       f2(i) = (-cinv*(d1(i) - a1(i)*f1(i-1)) + ainv*(d2(i) - a2(i)*f2(i-1))) / det

    end do

    ! Solution for top layer

    i = n
    t(i) = f1(i)
    q(i) = f2(i)

    ! Solution for layers through to bottom

    do i = n-1, 1, -1
       t(i) = f1(i) - e11(i) * t(i+1) - e12(i) * q(i+1) 
       q(i) = f2(i) - e21(i) * t(i+1) - e22(i) * q(i+1) 
    end do

  end subroutine tridiag_2eq

  !-----------------------------------------------------------------------
  function beta_function (p, q) result(beta)
    !
    ! !DESCRIPTION:
    ! Return the value of the beta function evaluated at p and q: B(p,q)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: p      ! Input argument
    real(r8), intent(in) :: q      ! Input argument
    !
    ! !LOCAL VARIABLES:
    real(r8) :: beta               ! Beta function: B(p,q)
    !---------------------------------------------------------------------

    beta = exp(log_gamma_function(p) + log_gamma_function(q) - log_gamma_function(p+q))
  
  end function beta_function

  !-----------------------------------------------------------------------
  function log_gamma_function (x) result(gammaln)
    !
    ! !DESCRIPTION:
    ! Return the value of the log natural of the gamma function evaluated at x: ln(G(x)) 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: x     ! Input argument
    !
    ! !LOCAL VARIABLES:
    real(r8) :: gammaln           ! ln(G(x))
    real(r8) :: y, tmp, ser
    integer :: j

    real(r8), parameter :: coef(6) = (/ 76.18009172947146_r8, -86.50532032941677_r8, &
    24.01409824083091_r8, -1.231739572450155_r8, 0.1208650973866179e-02_r8, -0.5395239384953e-05_r8 /)
    real(r8), parameter :: stp = 2.5066282746310005_r8
    !---------------------------------------------------------------------

    y = x
    tmp = x + 5.5_r8
    tmp = (x + 0.5_r8) * log(tmp) - tmp
    ser = 1.000000000190015_r8
    do j = 1, 6
       y = y + 1._r8
       ser = ser + coef(j) / y
    end do
    gammaln = tmp + log(stp * ser / x)

  end function log_gamma_function

end module MLMathToolsMod
