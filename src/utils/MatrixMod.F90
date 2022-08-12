module MatrixMod
!============================================================
!
! Module for linear alegebra matrix methods
!
!============================================================

#include "shr_assert.h"

  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils  , only : endrun
  use shr_log_mod , only : errMsg => shr_log_errMsg
  implicit none
  private

  !
  ! Public methods:
  !
  public inverse    ! Compute the inverse of a matrix

!============================================================
contains
!============================================================

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
  implicit none
  ! Arguments
  integer,intent(in) :: n           ! Size of matrix
  real(r8),intent(in)  :: a(:,:)    ! Input matrix to fine the inverse of
  real(r8),intent(out) :: c(:,:)    ! Output inverse
  ! Local variables
  real(r8) :: L(n,n)   ! matrix of the elimination coefficient
  real(r8) :: U(n,n)   ! Upper triangular part of input matrix A
  real(r8) :: aa(n,n)  ! Temporary equal to input matrix a
  real(r8) :: b(n)     ! Temporary vector
  real(r8) :: d(n)     ! Temporary vector (solution of L*d)
  real(r8) :: x(n)     ! Temporary vector (U*x = d)
  real(r8) :: coeff    ! coefficient
  integer i, j, k      ! Indices
  character(len=*), parameter :: subname = 'inverse'

  !
  ! Verify input matrix sizes
  !
  SHR_ASSERT((size(a,1) == n), errMsg(subname, __LINE__))
  SHR_ASSERT((size(a,2) == n), errMsg(subname, __LINE__))
  SHR_ASSERT((size(c,1) == n), errMsg(subname, __LINE__))
  SHR_ASSERT((size(c,2) == n), errMsg(subname, __LINE__))
  !
  ! Check that diagonals of input matrix aren't zero
  !
  do k=1,n
     if ( a(k,k) == 0.0_r8 )then
        call endrun( subname//" ERROR: A diagonal element of the input matrix is zero" )
        return
     end if
  end do
  !
  ! step 0: initialization for matrices L and U and b
  ! Fortran 90/95 aloows such operations on matrices
  !
  L=0.0
  U=0.0
  b=0.0

  aa=a
  !
  ! Step 1: forward elimination
  !
  do k=1, n-1
     do i=k+1,n
        ! Already verifieid that divisor isn't zero
        coeff=aa(i,k)/aa(k,k)
        L(i,k) = coeff
        do j=k+1,n
           aa(i,j) = aa(i,j)-coeff*aa(k,j)
        end do
     end do
  end do

  !
  ! Step 2: prepare L and U matrices
  ! L matrix is a matrix of the elimination coefficient
  ! + the diagonal elements are 1.0
  !
  do i=1,n
    L(i,i) = 1.0
  end do
  !
  ! U matrix is the upper triangular part of A
  !
  do j=1,n
    do i=1,j
      U(i,j) = aa(i,j)
    end do
  end do
  !
  ! Step 3: compute columns of the inverse matrix C
  !
  do k=1,n
    b(k)=1.0
    d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
    do i=2,n
      d(i)=b(i)
      do j=1,i-1
        d(i) = d(i) - L(i,j)*d(j)
      end do
    end do
    ! Step 3b: Solve Ux=d using the back substitution
    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
      x(i) = d(i)
      do j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
      end do
      ! Already verifieid that divisor isn't zero
      x(i) = x(i)/u(i,i)
    end do
    ! Step 3c: fill the solutions x(n) into column k of C
    do i=1,n
      c(i,k) = x(i)
    end do
    b(k)=0.0
  end do
end subroutine inverse

!============================================================

end module MatrixMod
