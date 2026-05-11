module SparseMatrixMultiplyMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: SparseMatrixMultiplyMod
  !
  ! !DESCRIPTION:
  ! Sparse matrix multiplication add addition
  !
  ! Author: Xingjie Lu
  !
  !EOP
  !-----------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use decompMod                      , only : bounds_type
  use clm_varctl                     , only : iulog
  use abortutils                     , only : endrun
  implicit none
  private
  !
  ! General Sparse matrix type
  !
  ! !PUBLIC TYPES:
  type, public :: sparse_matrix_type
 
     !sparse matrix is in COO format, Both row index and column index should be in ascending order. 
     !Row index should change faster than Column index to ensure SPMP_AB work properly.
 
     real(r8), pointer :: M(:,:) => null() ! non-zero entries in sparse matrix (unit,sparse matrix index)
     integer , pointer :: RI(:) => null()  ! Row index
     integer , pointer :: CI(:) => null()  ! Column index
     integer NE                            ! Number of nonzero entries
     integer SM                            ! Size of matrix, eg. for nxn matrix, SM=n
     integer num_unit                      ! number of active unit, such as patch, col, or gridcell
     integer begu                          ! begin index of unit in current process
     integer endu                          ! end index of unit in current process
  
  contains
    
    ! !PUBLIC MEMBER FUNCTIONS:
    procedure, public :: InitSM         ! subroutine to initilize sparse matrix type
    procedure, public :: ReleaseSM      ! subroutine to deallocate the sparse matrix type data
    procedure, public :: IsAllocSM      ! return true if the sparse matrix type is allocated (InitSM was called)
    procedure, public :: IsEquivIdxSM   ! return true if the sparse matrix indices are the same for the two sparce matrices
    procedure, public :: SetValueSM     ! subroutine to set values in sparse matrix of any shape
    procedure, public :: SetValueA      ! subroutine to set off-diagonal values in sparse matrix of A
    procedure, public :: SetValueA_diag ! subroutine to set diagonal values in sparse matrix of A
    procedure, public :: SetValueCopySM ! subroutine to copy the input sparse matrix to the output
    procedure, public :: CopyIdxSM      ! subroutine to copy the input indices to the sparse matrix
    procedure, public :: IsValuesSetSM  ! return true if the values are set in the matrix
    procedure, public :: SPMM_AK        ! subroutine to calculate sparse matrix multiplication: A(sparse matrix) = A(sparse matrix) * K(diagonal matrix)
    procedure, public :: SPMP_AB        ! subroutine to calculate sparse matrix addition AB(sparse matrix) = A(sparse matrix) + B(sparse matrix)
    procedure, public :: SPMP_B_ACC     ! subroutine to calculate sparse matrix accumulation: B(sparse matrix) = B(sparse matrix) + A(sparse matrix)
    procedure, public :: SPMP_ABC       ! subroutine to calculate sparse matrix addition ABC(sparse matrix) = A(sparse matrix) + B(sparse matrix) + C(sparse matrix)

  end type sparse_matrix_type
  !
  ! Diagonal matrix type
  !
  type, public :: diag_matrix_type
  
     ! Diaganal matrix only store diagnoal entries

     real(r8), pointer :: DM(:,:) => null() ! entries in diagonal matrix (unit,diagonal matrix index)
     integer SM                             ! Size of matrix, eg. for nxn matrix, SM=n
     integer num_unit                       ! number of active unit, such as patch, col, or gridcell
     integer begu                           ! begin index of unit in current process
     integer endu                           ! end index of unit in current process
  
  contains

    ! !PUBLIC MEMBER FUNCTIONS:
    procedure, public :: InitDM           ! subroutine to initialize diagonal matrix type
    procedure, public :: ReleaseDM        ! subroutine to deallocate the diagonal matrix
    procedure, public :: IsAllocDM        ! return true if the diagonal matrix is allocated (InitDM was called)
    procedure, public :: SetValueDM       ! subroutine to set values in diagonal matrix
    procedure, public :: SetValueCopyDM   ! subroutine to copy the input diagonal matrix to the output
    procedure, public :: IsValuesSetDM    ! return true if the values are set in the matrix

  end type diag_matrix_type
  !
  ! Vector type
  !
  type, public :: vector_type
  
     ! Vector format

     real(r8), pointer :: V(:,:) => null() ! entries in vector (unit,vector index)
     integer SV                            ! Size of vector
     integer num_unit                      ! number of active unit, such as patch, col, or gridcell
     integer begu                          ! begin index of unit in current process
     integer endu                          ! end index of unit in current process

  contains 
  
    ! !PUBLIC MEMBER FUNCTIONS:
    procedure, public :: InitV            ! subroutine to initialize vector type
    procedure, public :: ReleaseV         ! subroutine to deallocate veector type
    procedure, public :: IsAllocV         ! return true if the vector is allocated (InitV was called)
    procedure, public :: SetValueV        ! subroutine to set values in vector
    procedure, public :: SetValueV_scaler ! subroutine to set a constant value to a vector
    procedure, public :: SPMM_AX          ! subroutine to calculate multiplication X(vector)=A(sparse matrix)*X(vector)

  end type vector_type

  ! !PUBLIC DATA:
  integer,  public, parameter :: empty_int  = -9999
  real(r8), public, parameter :: empty_real = -9999._r8

  ! !PRIVATE DATA:
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

! ========================================================================

contains

  ! ========================================================================

  subroutine InitSM(this,SM_in,begu,endu,maxsm)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Initialize the sparse matrix by giving the boundary of landunit/gridcell/column/patch,
    ! the size of matrix. Then allocate the matrix in a sparse matrix format

    ! !ARGUMENTS:
    class(sparse_matrix_type) :: this
    integer,intent(in) :: SM_in
    integer,intent(in) :: begu
    integer,intent(in) :: endu
    integer,optional,intent(in) :: maxsm
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname = 'InitSM'
    !-----------------------------------------------------------------------

    if ( this%IsAllocSM() )then
       call endrun( subname//" ERROR: Sparse Matrix was already allocated" )
       return
    end if
    this%SM = SM_in
    this%begu = begu
    this%endu = endu
    if(present(maxsm))then
       SHR_ASSERT_FL((maxsm >= 1), sourcefile, __LINE__)
       SHR_ASSERT_FL((maxsm <= SM_in*SM_in), sourcefile, __LINE__)
       allocate(this%M(begu:endu,1:maxsm))
    else
       allocate(this%M(begu:endu,1:SM_in*SM_in))
    end if
    allocate(this%RI(1:SM_in*SM_in))
    allocate(this%CI(1:SM_in*SM_in))
    this%M(:,:) = empty_real
    this%RI(:) = empty_int
    this%CI(:) = empty_int
    this%NE    = empty_int

  end subroutine InitSM

  ! ========================================================================

  subroutine ReleaseSM(this)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Release the Sparse Matrix data

    ! !ARGUMENTS:
    class(sparse_matrix_type) :: this
    !-----------------------------------------------------------------------
  
    this%SM   = empty_int
    this%begu = empty_int
    this%endu = empty_int
    if ( associated(this%M) )then
       deallocate(this%M)
    end if
    if ( associated(this%RI) )then
       deallocate(this%RI)
    end if
    if ( associated(this%CI) )then
       deallocate(this%CI)
    end if
    this%M => null()
    this%RI=> null()
    this%CI=> null()

  end subroutine ReleaseSM

  ! ========================================================================

  logical function IsAllocSM(this)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Check if the Sparse Matrix has been allocated (InitSM was called on it)
    ! !ARGUMENTS:
    class(sparse_matrix_type) :: this
    !-----------------------------------------------------------------------

    if ( associated(this%M) .or. associated(this%RI) .or. associated(this%CI) )then
       IsAllocSM = .true.
    else
       IsAllocSM = .false.
    end if

  end function IsAllocSM
  

  ! ========================================================================

  logical function IsEquivIdxSM(this, A)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Check if the Sparse Matrix indices are eqiuivalent
    ! !ARGUMENTS:
    class(sparse_matrix_type) :: this
    type(sparse_matrix_type), intent(in)  :: A   ! Sparse matrix indices to compare to
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname = 'IsEquivIdxSM'
    !-----------------------------------------------------------------------

    ! Start checking easy critera and return if can determine status for sure,
    ! keep checking harder things until everything has been checked for
    if ( this%SM /= A%SM )then
       IsEquivIdxSM = .false.
       return
    end if
    if ( this%NE == A%NE )then
        ! If NE is the same and the row and column indices are identical -- the
        ! indices of the two arrays are identical
        if ( all(this%RI(:this%NE) == A%RI(:this%NE)) .and. all(this%CI(:this%NE) == A%CI(:this%NE)) )then
           IsEquivIdxSM = .true.
           return
        else
           ! This needs more checking! The order could be different
           IsEquivIdxSM = .false.
           return
        end if
    else
       ! This needs more checking! There could be some zerod entries in
       ! non-zero positions
       IsEquivIdxSM = .false.
       return
    end if
    call endrun( subname//" ERROR: it should NOT be possible to reach this point" )
    return

  end function IsEquivIdxSM
  
  ! ========================================================================

  subroutine SetValueSM(this,begu,endu,num_unit,filter_u,M,I,J,NE_in)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Set sparse matrix values by giving all non-zero values and the corresponding row and column indices.
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.
    !
    ! !ARGUMENTS:
    class(sparse_matrix_type) :: this
    integer ,intent(in) :: begu
    integer ,intent(in) :: endu
    integer ,intent(in) :: NE_in
    integer ,intent(in) :: num_unit
    integer ,intent(in) :: filter_u(:)
    real(r8),intent(in) :: M(begu:,1:)
    integer ,intent(in) :: I(:)
    integer ,intent(in) :: J(:)
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname = 'SetValueSM'
    integer k,u,fu
    !-----------------------------------------------------------------------

    if ( .not. this%IsAllocSM() )then
       call endrun( subname//" ERROR: Sparse Matrix was NOT already allocated" )
       return
    end if
    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M, 2) >= NE_in), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(I, 1) >= NE_in), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(J, 1) >= NE_in), sourcefile, __LINE__)
    SHR_ASSERT_FL((lbound(M, 1) == begu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M, 1) == endu), sourcefile, __LINE__)
#ifndef _OPENMP
    ! Without OpenMP array sizes will be identical
    SHR_ASSERT_FL((lbound(M, 1) == this%begu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M, 1) == this%endu), sourcefile, __LINE__)
#else
    ! With OpenMP the allocated array sizes might be larger than the input ones
    SHR_ASSERT_FL((lbound(M, 1) >= this%begu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M, 1) <= this%endu), sourcefile, __LINE__)
#endif
    SHR_ASSERT_FL((maxval(I(:this%NE)) <= this%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((minval(I(:this%NE)) >= 1), sourcefile, __LINE__)
    SHR_ASSERT_FL((maxval(J(:this%NE)) <= this%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((minval(J(:this%NE)) >= 1), sourcefile, __LINE__)
    do k = 1,NE_in
       do fu = 1,num_unit
          u = filter_u(fu)  
          this%M(u,k) = M(u,k)
       end do
    end do

    this%NE = NE_in
    do k = 1,NE_in
       this%RI(k) = I(k)
       this%CI(k) = J(k)
    end do

  end subroutine SetValueSM

  ! ========================================================================

  subroutine SetValueA_diag(this,num_unit,filter_u,scaler)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Set diagonal sparse matrix values by giving a constant scaler.
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.

    ! !ARGUMENTS:
    class(sparse_matrix_type) :: this
    real(r8),intent(in) :: scaler
    integer,intent(in) :: num_unit
    integer,intent(in) :: filter_u(:)
    ! !LOCAL VARIABLES:
    integer i,u,fu
    character(len=*),parameter :: subname = 'SetValueA_diag'
    !-----------------------------------------------------------------------

    if ( .not. this%IsAllocSM() )then
       call endrun( subname//" ERROR: Sparse Matrix was NOT already allocated" )
       return
    end if
    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((lbound(this%M,1) == this%begu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(this%M,1) == this%endu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(this%M,2) >= this%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(this%RI,1) >= this%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(this%CI,1) >= this%SM), sourcefile, __LINE__)
    do i=1,this%SM
       do fu=1,num_unit
          u = filter_u(fu)
          this%M(u,i)  = scaler
       end do
    end do

    do i=1,this%SM
       this%RI(i) = i
       this%CI(i) = i
    end do
    this%NE = this%SM

  end subroutine SetValueA_diag

  ! ========================================================================

  subroutine SetValueA(this,begu,endu,num_unit,filter_u,M,AI,AJ,NE_NON,Init_ready,list,RI_A,CI_A)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Set sparse matrix values by giving values, rows, and columns of non-zero and non-diagonal entries.
    ! Then Set the diagonal entries to -1. The information of active landunit/gridcell/column/patch,
    ! The order and indices of non-diagonal entries in full sparse matrix are memorized to save computational cost,
    ! since these indices are usualy time-independent.

    ! !ARGUMENTS:
    class(sparse_matrix_type) :: this
    integer ,intent(in) :: begu
    integer ,intent(in) :: endu
    integer ,intent(in) :: NE_NON
    integer ,intent(in) :: num_unit
    integer ,intent(in) :: filter_u(:)
    real(r8),intent(in) :: M(begu:,1:)
    integer ,intent(in) :: AI(:)
    integer ,intent(in) :: AJ(:)
    logical ,intent(inout) :: Init_ready !True: diagnoal of A has been set to -1,this%RI, this%CI, this%NE and list has been set up
    integer ,intent(inout),optional :: list(:)
    integer ,intent(inout),optional :: RI_A(:)
    integer ,intent(inout),optional :: CI_A(:)
    ! !LOCAL VARIABLES:
    integer i,j,k,fu,u
    logical list_ready
    type(sparse_matrix_type) :: A_diag, A_nondiag
    character(len=*),parameter :: subname = 'SetValueA'
    !-----------------------------------------------------------------------

    list_ready = .false.

    if ( .not. this%IsAllocSM() )then
       call endrun( subname//" ERROR: Sparse Matrix was NOT already allocated" )
       return
    end if
    if(init_ready .and. .not. (present(list) .and. present(RI_A) .and. present(CI_A)))then
       write(iulog,*) "Error: initialization is ready, but at least one of list, RI_A or CI_A is not presented"
       call endrun( subname//" ERROR: required optional arguments were NOT sent in" )
       return
    end if
    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((lbound(M,1) == begu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M,1) == endu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M,2) >= NE_NON), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(AI,1) >= NE_NON), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(AJ,1) >= NE_NON), sourcefile, __LINE__)
    if ( present(list) )then
       SHR_ASSERT_FL((ubound(list,1) >= NE_NON), sourcefile, __LINE__)
    end if
    if ( present(RI_A) )then
       SHR_ASSERT_FL((ubound(RI_A,1) >= NE_NON+this%SM), sourcefile, __LINE__)
    end if
    if ( present(CI_A) )then
       SHR_ASSERT_FL((ubound(CI_A,1) >= NE_NON+this%SM), sourcefile, __LINE__)
    end if
 
    if(Init_ready)then
       do i = 1,this%SM+NE_NON
          do fu = 1,num_unit
             u = filter_u(fu)
             this%M(u,i) = -1._r8
          end do
       end do
       do i = 1,NE_NON
          do fu = 1,num_unit
             u = filter_u(fu)
             this%M(u,list(i)) = M(u,i)
          end do
       end do
       this%NE = this%SM+NE_NON
       this%RI(1:this%NE) = RI_A(1:this%NE)
       this%CI(1:this%NE) = CI_A(1:this%NE)
    else
       if ( A_diag%IsAllocSM()    ) call A_diag%ReleaseSM()
       if ( A_nondiag%IsAllocSM() ) call A_nondiag%ReleaseSM()
       call A_diag%InitSM(this%SM,begu,endu)
       call A_nondiag%InitSM(this%SM,begu,endu)

       call A_diag%SetValueA_diag(num_unit,filter_u,-1._r8)
       call A_nondiag%SetValueSM(begu,endu,num_unit,filter_u,M,AI,AJ,NE_NON)

       if(present(list))then
          call this%SPMP_AB(num_unit,filter_u,A_nondiag,A_diag,list_ready,list_A=list)
       else
          call this%SPMP_AB(num_unit,filter_u,A_nondiag,A_diag,list_ready)
       end if
       if(present(RI_A))RI_A(1:this%NE) = this%RI(1:this%NE)
       if(present(CI_A))CI_A(1:this%NE) = this%CI(1:this%NE)
   
       Init_ready = .true.
       call A_diag%ReleaseSM()
       call A_nondiag%ReleaseSM()
    end if

  end subroutine SetValueA

  ! ========================================================================

  subroutine SetValueCopySM(this, num_unit, filter_u, matrix)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Set the sparse matrix by copying from another sparse matrix

    ! !ARGUMENTS:
    class(sparse_matrix_type) :: this
    type(sparse_matrix_type), intent(in) :: matrix   ! Sparse Matrix to copy
    integer ,intent(in) :: num_unit
    integer ,intent(in) :: filter_u(:)
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname = 'SetValueCopySM'
    !-----------------------------------------------------------------------

    if ( .not. this%IsAllocSM() )then
       call endrun( subname//" ERROR: Sparse Matrix was NOT already allocated" )
       return
    end if
    if ( .not. matrix%IsValuesSetSM() )then
       call endrun( subname//" ERROR: Sparse Matrix data sent in was NOT already set" )
       return
    end if
    SHR_ASSERT_FL( (this%SM   == matrix%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL( (this%begu == matrix%begu), sourcefile, __LINE__)
    SHR_ASSERT_FL( (this%endu == matrix%endu), sourcefile, __LINE__)
    SHR_ASSERT_FL((maxval(matrix%RI(:this%NE)) <= this%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((minval(matrix%RI(:this%NE)) >= 1), sourcefile, __LINE__)
    SHR_ASSERT_FL((maxval(matrix%CI(:this%NE)) <= this%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((minval(matrix%CI(:this%NE)) >= 1), sourcefile, __LINE__)
    call this%SetValueSM( matrix%begu, matrix%endu, num_unit, filter_u, matrix%M, &
                          matrix%RI, matrix%CI, matrix%NE)

  end subroutine SetValueCopySM

  ! ========================================================================

  subroutine SetValueCopyDM(this, num_unit, filter_u, matrix)
     !-----------------------------------------------------------------------
     !
     ! !DESCRIPTION:
     ! Set the sparse matrix by copying from another sparse matrix

     ! !ARGUMENTS:
     class(diag_matrix_type) :: this
     type(diag_matrix_type), intent(in) :: matrix   ! Diagonal Matrix to copy
     integer ,intent(in) :: num_unit
     integer ,intent(in) :: filter_u(:)
     ! !LOCAL VARIABLES:
     character(len=*),parameter :: subname = 'SetValueCopyDM'
     !-----------------------------------------------------------------------

     if ( .not. this%IsAllocDM() )then
        call endrun( subname//" ERROR: Diagonal Matrix was NOT already allocated" )
        return
     end if
     if ( .not. matrix%IsValuesSetDM() )then
        call endrun( subname//" ERROR: Diagonal Matrix data sent in was NOT already set" )
        return
     end if
     SHR_ASSERT_FL( (this%SM   == matrix%SM), sourcefile, __LINE__)
     SHR_ASSERT_FL( (this%begu == matrix%begu), sourcefile, __LINE__)
     SHR_ASSERT_FL( (this%endu == matrix%endu), sourcefile, __LINE__)
     call this%SetValueDM( matrix%begu, matrix%endu, num_unit, filter_u, matrix%DM)

  end subroutine SetValueCopyDM

  ! ========================================================================

  subroutine CopyIdxSM(this, matrix)
     !-----------------------------------------------------------------------
     !
     ! !DESCRIPTION:
     ! Copy the indices from the input matrix to this sparse matrix
     ! also make sure the sizes are consistent

     ! !ARGUMENTS:
     class(sparse_matrix_type) :: this
     type(sparse_matrix_type), intent(in) :: matrix   ! Sparse Matrix to copy
     ! !LOCAL VARIABLES:
     character(len=*),parameter :: subname = 'CopyIdxSM'
     integer :: i
     !-----------------------------------------------------------------------

     if ( .not. this%IsAllocSM() )then
        call endrun( subname//" ERROR: Sparse Matrix was NOT already allocated" )
        return
     end if
     if ( .not. matrix%IsValuesSetSM() )then
        call endrun( subname//" ERROR: Sparse Matrix data sent in was NOT already set" )
        return
     end if
     SHR_ASSERT_FL( (this%SM   == matrix%SM), sourcefile, __LINE__)
     SHR_ASSERT_FL( (this%begu == matrix%begu), sourcefile, __LINE__)
     SHR_ASSERT_FL( (this%endu == matrix%endu), sourcefile, __LINE__)
     SHR_ASSERT_FL((maxval(matrix%RI(:matrix%NE)) <= this%SM), sourcefile, __LINE__)
     SHR_ASSERT_FL((minval(matrix%RI(:matrix%NE)) >= 1), sourcefile, __LINE__)
     SHR_ASSERT_FL((maxval(matrix%CI(:matrix%NE)) <= this%SM), sourcefile, __LINE__)
     SHR_ASSERT_FL((minval(matrix%CI(:matrix%NE)) >= 1), sourcefile, __LINE__)
     !
     ! Figure out the number of non-empty data values and make sure it's same as input
     !
     this%NE = size(this%M,2)
     do i = 1, this%NE
       if ( all(this%M(:,i) == empty_int) )then
          this%NE = i-1
          exit
       end if
     end do
     if ( this%NE /= matrix%NE )then
        call endrun( subname//" ERROR: Sparse Matrix empty data size is different from input one copying the indices from" )
        return
     end if
     !
     ! Copy indices
     !
     this%RI(:this%NE) = matrix%RI(:matrix%NE)
     this%CI(:this%NE) = matrix%CI(:matrix%NE)

  end subroutine CopyIdxSM

  ! ========================================================================

  logical function IsValuesSetSM(this)
     !-----------------------------------------------------------------------
     !
     ! !DESCRIPTION:
     ! Check if the Sparse Matrix has it's data been set (One of the SetValue* subroutines was called on it)

     ! !ARGUMENTS:
     class(sparse_matrix_type) :: this
     !-----------------------------------------------------------------------

     if ( .not. this%IsAllocSM() )then
        IsValuesSetSM = .false.
     else if ( this%NE == empty_int )then
        IsValuesSetSM = .false.
     else
        IsValuesSetSM = .true.
     end if

  end function IsValuesSetSM

  ! ========================================================================

  logical function IsValuesSetDM(this)
     !-----------------------------------------------------------------------
     !
     ! !DESCRIPTION:
     ! Check if the Sparse Matrix has it's data been set (One of the SetValue* subroutines was called on it)

     ! !ARGUMENTS:
     class(diag_matrix_type) :: this
     !-----------------------------------------------------------------------

     if ( .not. this%IsAllocDM() )then
        IsValuesSetDM = .false.
     else if ( this%SM == empty_int )then
        IsValuesSetDM = .false.
     else
        IsValuesSetDM = .true.
     end if

  end function IsValuesSetDM

  ! ========================================================================

  subroutine InitDM(this,SM_in,begu,endu)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Initialize the diagonal matrix by giving the boundary of landunit/gridcell/column/patch,
    ! the size of matrix. Then allocate the matrix in a diagonal matrix format 

    ! !ARGUMENTS:
    class(diag_matrix_type) :: this
    integer,intent(in) :: SM_in
    integer,intent(in) :: begu
    integer,intent(in) :: endu
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname = 'InitDM'
    !-----------------------------------------------------------------------

    if ( this%IsAllocDM() )then
       call endrun( subname//" ERROR: Diagonal Matrix was already allocated" )
       return
    end if
    this%SM = SM_in
    allocate(this%DM(begu:endu,1:SM_in))
    this%DM(:,:) = empty_real
    this%begu = begu
    this%endu = endu

  end subroutine InitDM

  !-----------------------------------------------------------------------

  subroutine ReleaseDM(this)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Release the Diagonal Matrix data

    ! !ARGUMENTS:
    class(diag_matrix_type) :: this
    !-----------------------------------------------------------------------

    this%SM   = empty_int
    this%begu = empty_int
    this%endu = empty_int
    if ( associated(this%DM) )then
       deallocate(this%DM)
    end if
    this%DM => null()
  end subroutine ReleaseDM

  !-----------------------------------------------------------------------

  logical function IsAllocDM(this)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Check if the Diagonal Matrix is allocated (InitDM was called)

    ! !ARGUMENTS:
    class(diag_matrix_type) :: this
    !-----------------------------------------------------------------------

    if ( associated(this%DM) )then
       IsAllocDM = .true.
    else
       IsAllocDM = .false.
    end if

  end function IsAllocDM

  !-----------------------------------------------------------------------

  subroutine SetValueDM(this,begu,endu,num_unit,filter_u,M)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Set the diagonal matrix values by giving the values of diagonal entries in a right order.
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.

    ! !ARGUMENTS:
    class(diag_matrix_type) :: this
    integer ,intent(in) :: begu
    integer ,intent(in) :: endu
    real(r8),intent(in) :: M(begu:,1:)
    integer ,intent(in) :: num_unit
    integer ,intent(in) :: filter_u(:)
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname = 'SetValueDM'
    integer i,fu,u
    !-----------------------------------------------------------------------

    if ( .not. this%IsAllocDM() )then
       call endrun( subname//" ERROR: Diagonal matrix was NOT already allocated" )
       return
    end if
    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((lbound(M,1)        == begu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M,1)        == endu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M,2)        >= this%SM), sourcefile, __LINE__)
    do i = 1,this%SM
       do fu = 1,num_unit
          u = filter_u(fu)
          this%DM(u,i) = M(u,i)
       end do
    end do

  end subroutine SetValueDM

  !-----------------------------------------------------------------------

  subroutine InitV(this,SV_in,begu,endu)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Initialize the vector by giving the boundary of landunit/gridcell/column/patch,
    ! the size of vector. Then allocate the vector in a vector type

    ! !ARGUMENTS:
    class(vector_type) :: this
    integer,intent(in) :: SV_in
    integer,intent(in) :: begu
    integer,intent(in) :: endu
    ! !LOCAL VARIABLES:
    character(len=*),parameter :: subname = 'InitV'
    !-----------------------------------------------------------------------

    if ( this%IsAllocV() )then
      call endrun( subname//" ERROR: Vector was already allocated" )
      return
    end if
    this%SV = SV_in
    allocate(this%V(begu:endu,1:SV_in))
    this%V(:,:) = empty_real
    this%begu = begu
    this%endu = endu

  end subroutine InitV

  !-----------------------------------------------------------------------

  subroutine ReleaseV(this)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Deallocate vector type

    ! !ARGUMENTS:
    class(vector_type) :: this
    !-----------------------------------------------------------------------

    if ( associated(this%V) )then
       deallocate(this%V)
    end if
    this%V => null()
    this%begu = empty_int
    this%endu = empty_int
    this%SV = empty_int

  end subroutine ReleaseV

  ! ========================================================================

  logical function IsAllocV(this)
     !-----------------------------------------------------------------------
     !
     ! !DESCRIPTION:
     ! Check if the Vector has been allocated (InitV was called on it)
  
     ! !ARGUMENTS:
     class(vector_type) :: this
     !-----------------------------------------------------------------------

     if ( associated(this%V) )then
        IsAllocV = .true.
     else
        IsAllocV = .false.
     end if

  end function IsAllocV

  !-----------------------------------------------------------------------

  subroutine SetValueV_scaler(this,num_unit,filter_u,scaler)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Set the vector values by giving a constant value
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.

    ! !ARGUMENTS:
    class(vector_type) :: this
    real(r8),intent(in) :: scaler
    integer,intent(in) :: num_unit
    integer,intent(in) :: filter_u(:)
    ! !LOCAL VARIABLES:
    integer i,fu,u
    character(len=*),parameter :: subname = 'SetValueV_scaler'
    !-----------------------------------------------------------------------

    if ( .not. this%IsAllocV() )then
       call endrun( subname//" ERROR: Vector was NOT already allocated" )
       return
    end if
    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    do i=1,this%SV
       do fu = 1,num_unit
          u = filter_u(fu)
          this%V(u,i) = scaler
       end do
    end do

  end subroutine SetValueV_scaler

  !-----------------------------------------------------------------------

  subroutine SetValueV(this,begu,endu,num_unit,filter_u,M)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Set the vector values by giving the values in a right order.
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.

    ! !ARGUMENTS:
    integer ,intent(in) :: begu
    integer ,intent(in) :: endu
    class(vector_type)  :: this
    real(r8),intent(in) :: M(begu:,1:)
    integer ,intent(in) :: num_unit
    integer ,intent(in) :: filter_u(:)
    ! !LOCAL VARIABLES:
    integer i,fu,u
    character(len=*),parameter :: subname = 'SetValueV'
    !-----------------------------------------------------------------------

    if ( .not. this%IsAllocV() )then
       call endrun( subname//" ERROR: Vector was NOT already allocated" )
       return
    end if
    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((lbound(M,1)        == begu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M,1)        == endu), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(M,2)        >= this%SV), sourcefile, __LINE__)
    do i=1,this%SV
       do fu = 1,num_unit
          u = filter_u(fu)
          this%V(u,i) = M(u,i)
       end do
    end do

  end subroutine SetValueV

  !-----------------------------------------------------------------------

  subroutine SPMM_AK(this,num_unit,filter_u,K)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Calculate sparse matrix multiplication (SPMM) A(this) = A(this)*K
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.
    ! A is a sparse matrix in Coordinate format (COO). 
    ! K is a diagnoal matrix.

    ! !ARGUMENTS:
    class(sparse_matrix_type)  :: this
    type(diag_matrix_type)  ,intent(in)  :: K
    integer,intent(in) :: num_unit
    integer,intent(in) :: filter_u(:)
    ! !LOCAL VARIABLES:
    integer i,fu,u
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SM            == K%SM), sourcefile, __LINE__)
    do i=1,this%NE
       do fu = 1,num_unit
          u = filter_u(fu)
          this%M(u,i) = this%M(u,i) * K%DM(u,this%CI(i))
       end do
    end do

  end subroutine SPMM_AK

  !-----------------------------------------------------------------------

  subroutine SPMM_AX(this,num_unit,filter_u,A)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Calculate sparse matrix multiplication (SPMM) X(this) = X(this) + A*X(this)
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.
    ! A is a sparse matrix in Coordinate format (COO). 
    ! X is a vector type.

    ! !ARGUMENTS:
    class(vector_type)         :: this
    type(sparse_matrix_type),intent(in)    :: A
    integer,intent(in) :: num_unit
    integer,intent(in) :: filter_u(:)
    ! !LOCAL VARIABLES:
    integer i,fu,u
    real(r8) :: V(this%begu:this%endu,1:this%SV)
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(this%V,2)   == this%SV), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SV            <= A%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(this%V,2)   == this%SV), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(A%M,2)      >= A%NE), sourcefile, __LINE__)
    SHR_ASSERT_FL((maxval(A%RI)       <= this%SV), sourcefile, __LINE__)
    SHR_ASSERT_FL((maxval(A%CI)       <= this%SV), sourcefile, __LINE__)
    do i=1,this%SV
       do fu = 1, num_unit
          u = filter_u(fu)
          V(u,i) = this%V(u,i)
       end do
    end do

    do i=1,A%NE
       do fu = 1, num_unit
          u = filter_u(fu)
          this%V(u,A%RI(i)) = this%V(u,A%RI(i)) + A%M(u,i) * V(u,A%CI(i))
       end do
    end do

  end subroutine SPMM_AX

  !-----------------------------------------------------------------------

  subroutine SPMP_B_ACC(this,num_unit,filter_u,A)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Calculate sparse matrix addition (SPMP) B(this) = B(this) + A
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.
    ! A and B are sparse matrix in Coordinate format (COO). 
    ! Entry locations of A and B should be the same.

    ! !ARGUMENTS:
    class(sparse_matrix_type)     :: this
    type(sparse_matrix_type),intent(in)     :: A
    integer,intent(in) :: num_unit
    integer,intent(in) :: filter_u(:)
    ! !LOCAL VARIABLES:
    integer i,fu,u
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SM            == A%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%NE            == A%NE), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((this%RI        == A%RI), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((this%CI        == A%CI), sourcefile, __LINE__)

    do i=1,A%NE
       do fu = 1, num_unit
          u = filter_u(fu)
          this%M(u,i) = this%M(u,i) + A%M(u,i)
       end do
    end do   

  end subroutine SPMP_B_ACC

  !-----------------------------------------------------------------------

  subroutine SPMP_AB(this,num_unit,filter_u,A,B,list_ready,list_A,list_B,NE_AB,RI_AB,CI_AB)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Calculate sparse matrix addition (SPMP) AB(this) =  A + B
    ! The map of each entry in A and B to AB have been memorized to save the computational cost, 
    ! since they are usually time-independent. 
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.
    ! A is a sparse matrix in Coordinate format (COO)
    ! B is a sparse matrix in Coordinate format (COO)
    ! AB is a sparse matrix in Coordinate format (COO)

    ! !ARGUMENTS:
    class(sparse_matrix_type)     :: this
    type(sparse_matrix_type),intent(in)     :: A
    type(sparse_matrix_type),intent(in)     :: B
    logical,intent(inout)                   :: list_ready 
    integer,intent(in) :: num_unit
    integer,intent(in) :: filter_u(:)

    integer,intent(inout),optional :: list_A(:)
    integer,intent(inout),optional :: list_B(:)
    integer,intent(inout),optional :: NE_AB
    integer,intent(inout),optional :: RI_AB(:)
    integer,intent(inout),optional :: CI_AB(:)

    integer,dimension(:) :: Aindex(A%NE+1),Bindex(B%NE+1)
    integer,dimension(:) :: ABindex(this%SM*this%SM)
    ! !LOCAL VARIABLES:
    integer i_a,i_b,i_ab
    integer i,fu,u
    character(len=*),parameter :: subname = 'SPMP_AB'
    !-----------------------------------------------------------------------

    ! 'list_ready = .true.' means list_A, list_B, NE_AB, RI_AB, and CI_AB have been memorized before.
    ! In this case they all need to be presented. Otherwise, use 'list_ready = .false.' to get those information
    ! for the first time call this subroutine.

    if ( present(list_A) )then
       SHR_ASSERT_FL((ubound(list_A,1) >= A%NE), sourcefile, __LINE__)
    end if
    if ( present(list_B) )then
       SHR_ASSERT_FL((ubound(list_B,1) >= B%NE), sourcefile, __LINE__)
    end if
    if ( present(RI_AB) )then
       SHR_ASSERT_FL((ubound(RI_AB,1) >= A%NE+B%NE), sourcefile, __LINE__)
    end if
    if ( present(CI_AB) )then
       SHR_ASSERT_FL((ubound(CI_AB,1) >= A%NE+B%NE), sourcefile, __LINE__)
    end if
    if(list_ready .and. .not. (present(list_A) .and. present(list_B) .and. present(NE_AB) .and. present(RI_AB) .and. present(CI_AB)))then
       write(iulog,*) "error in SPMP_AB: list_ready is True, but at least one of list_A, list_B, NE_AB, RI_AB and CI_AB are not presented"
       call endrun( subname//" ERROR: missing required optional arguments" )
       return
    end if
    SHR_ASSERT_FL((ubound(filter_u,1) >= num_unit), sourcefile, __LINE__)
    SHR_ASSERT_FL((A%NE               > 0), sourcefile, __LINE__)
    SHR_ASSERT_FL((B%NE               > 0), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SM            > 0), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SM            == A%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SM            == B%SM), sourcefile, __LINE__)

    if(.not. list_ready)then
       i_a=1
       i_b=1
       i_ab=1
       Aindex(1:A%NE)    = (A%CI(1:A%NE)-1)*A%SM + A%RI(1:A%NE)
       Bindex(1:B%NE)    = (B%CI(1:B%NE)-1)*B%SM + B%RI(1:B%NE)
       Aindex(A%NE+1)    = A%SM*A%SM + 1
       Bindex(B%NE+1)    = B%SM*B%SM + 1

       do while (i_a .le. A%NE .or. i_b .le. B%NE)
          if(Aindex(i_a) .lt. Bindex(i_b))then
             do fu = 1, num_unit
                u = filter_u(fu)
                this%M(u,i_ab) = A%M(u,i_a)
             end do
             ABindex(i_ab) = Aindex(i_a)
             if(present(list_A))list_A(i_a) = i_ab
             i_a  = i_a  + 1
             i_ab = i_ab + 1
          else
             if(Aindex(i_a) .gt. Bindex(i_b))then 
                do fu = 1, num_unit
                   u = filter_u(fu)
                   this%M(u,i_ab) = B%M(u,i_b)
                end do
                ABindex(i_ab) = Bindex(i_b)
                if(present(list_B))list_B(i_b) = i_ab
                i_b  = i_b  + 1
                i_ab = i_ab + 1
             else
                do fu = 1, num_unit
                   u = filter_u(fu)
                   this%M(u,i_ab) = A%M(u,i_a) + B%M(u,i_b)
                end do
                ABindex(i_ab) = Aindex(i_a)
                if(present(list_A))list_A(i_a) = i_ab
                if(present(list_B))list_B(i_b) = i_ab
                i_a  = i_a  + 1
                i_b  = i_b  + 1
                i_ab = i_ab + 1
             end if
          end if
       end do

       this%NE = i_ab - 1
       this%CI(1:this%NE) = (ABindex(1:this%NE) - 1) / this%SM + 1
       this%RI(1:this%NE) =  ABindex(1:this%NE) - this%SM * (this%CI(1:this%NE) - 1)
       if(present(NE_AB))NE_AB = this%NE
       if(present(CI_AB))CI_AB(1:this%NE) = this%CI(1:this%NE)
       if(present(RI_AB))RI_AB(1:this%NE) = this%RI(1:this%NE)
       if(present(list_A) .and. present(list_B) .and. present(NE_AB) .and. present(RI_AB) .and. present(CI_AB))list_ready = .true.
    else
       do i = 1, NE_AB
          do fu = 1, num_unit
             u = filter_u(fu)
             this%M(u,i) = 0._r8
          end do
       end do
       do i_a = 1, A%NE
          do fu = 1, num_unit
             u = filter_u(fu)
             this%M(u,list_A(i_a)) = A%M(u,i_a)
          end do
       end do
       do i_b = 1, B%NE
          do fu = 1, num_unit
             u = filter_u(fu)
             this%M(u,list_B(i_b)) = this%M(u,list_B(i_b)) + B%M(u,i_b)
          end do
       end do
       this%NE = NE_AB
       this%CI(1:this%NE) = CI_AB(1:NE_AB)
       this%RI(1:this%NE) = RI_AB(1:NE_AB)
    end if

  end subroutine SPMP_AB

  !-----------------------------------------------------------------------

  subroutine SPMP_ABC(this,num_unit,filter_u,A,B,C,list_ready,list_A,list_B,list_C,NE_ABC,RI_ABC,CI_ABC,&
             use_actunit_list_A,num_actunit_A,filter_actunit_A,use_actunit_list_B,num_actunit_B,filter_actunit_B,&
             use_actunit_list_C,num_actunit_C,filter_actunit_C)
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Calculate sparse matrix addition (SPMP) ABC(this) =  A + B + C
    ! The map of each entry in A, B and C to ABC have been memorized to save the computational cost, 
    ! since they are usually time-independent. 
    ! The information of active landunit/gridcell/column/patch is used to save computational cost.
    ! A is a sparse matrix in Coordinate format (COO)
    ! B is a sparse matrix in Coordinate format (COO)
    ! C is a sparse matrix in Coordinate format (COO)
    ! ABC is a sparse matrix in Coordinate format (COO)

    ! !ARGUMENTS:
    class(sparse_matrix_type)     :: this
    type(sparse_matrix_type),intent(in)     :: A
    type(sparse_matrix_type),intent(in)     :: B
    type(sparse_matrix_type),intent(in)     :: C
    logical,intent(inout)                   :: list_ready 
    integer,intent(in) :: num_unit
    integer,intent(in) :: filter_u(:)
    logical,intent(in),optional :: use_actunit_list_A
    logical,intent(in),optional :: use_actunit_list_B
    logical,intent(in),optional :: use_actunit_list_C
    integer,intent(in),optional :: num_actunit_A
    integer,intent(in),optional :: num_actunit_B
    integer,intent(in),optional :: num_actunit_C
    integer,dimension(:),intent(in),optional :: filter_actunit_A
    integer,dimension(:),intent(in),optional :: filter_actunit_B
    integer,dimension(:),intent(in),optional :: filter_actunit_C

    integer,intent(inout),optional :: list_A(:)
    integer,intent(inout),optional :: list_B(:)
    integer,intent(inout),optional :: list_C(:)
    integer,intent(inout),optional :: NE_ABC
    integer,intent(inout),optional :: RI_ABC(:)
    integer,intent(inout),optional :: CI_ABC(:)

    ! !LOCAL VARIABLES:
    integer,dimension(:) :: Aindex(A%NE+1),Bindex(B%NE+1),Cindex(C%NE+1)
    integer,dimension(:) :: ABCindex(this%SM*this%SM)

    integer i_a,i_b,i_c,i_abc
    integer i,fu,u
    character(len=*),parameter :: subname = 'SPMP_ABC'
    !-----------------------------------------------------------------------

    ! 'list_ready = .true.' means list_A, list_B, list_C, NE_ABC, RI_ABC, and CI_ABC have been memorized before.
    ! In this case they all need to be presented. Otherwise, use 'list_ready = .false.' to get those information
    ! for the first time call this subroutine.
    
    SHR_ASSERT_FL((this%SM            > 0), sourcefile, __LINE__)
    SHR_ASSERT_FL((A%NE               > 0), sourcefile, __LINE__)
    SHR_ASSERT_FL((B%NE               > 0), sourcefile, __LINE__)
    SHR_ASSERT_FL((C%NE               > 0), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SM            == A%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SM            == B%SM), sourcefile, __LINE__)
    SHR_ASSERT_FL((this%SM            == C%SM), sourcefile, __LINE__)
    if( present(list_A) )then
       SHR_ASSERT_FL((size(list_A)    >= A%NE), sourcefile, __LINE__)
    end if
    if( present(list_B) )then
       SHR_ASSERT_FL((size(list_B)    >= B%NE), sourcefile, __LINE__)
    end if
    if( present(list_C) )then
       SHR_ASSERT_FL((size(list_C)    >= C%NE), sourcefile, __LINE__)
    end if
    if( present(RI_ABC) )then
       SHR_ASSERT_FL((size(RI_ABC)    >= A%NE+B%NE+C%NE), sourcefile, __LINE__)
    end if
    if( present(CI_ABC) )then
       SHR_ASSERT_FL((size(CI_ABC)    >= A%NE+B%NE+C%NE), sourcefile, __LINE__)
    end if
    if(list_ready .and. .not. (present(list_A) .and. present(list_B) .and. present(list_C) .and. present(NE_ABC) .and. present(RI_ABC) .and. present(CI_ABC)))then
       write(iulog,*) "error in SPMP_ABC: list_ready is True, but at least one of list_A, list_B, list_C, NE_ABC, RI_ABC and CI_ABC are not presented",&
                  present(list_A),present(list_B),present(list_C),present(NE_ABC),present(RI_ABC),present(CI_ABC)
       call endrun( subname//" ERROR: missing required optional arguments" )
       return
    end if
    if(present(num_actunit_A))then
       if(num_actunit_A < 0)then
          write(iulog,*) "error: num_actunit_A cannot be less than 0"
          call endrun( subname//" ERROR: bad value for num_actunit_A" )
          return
       end if
       if(.not. present(filter_actunit_A))then
          write(iulog,*) "error: num_actunit_A is presented but filter_actunit_A is missing"
          call endrun( subname//" ERROR: missing required optional arguments" )
          return
       end if
       SHR_ASSERT_FL((size(filter_actunit_A) >= num_actunit_A), sourcefile, __LINE__)
    end if
    if(present(num_actunit_B))then
       if(num_actunit_B < 0)then
          write(iulog,*) "error: num_actunit_B cannot be less than 0"
          call endrun( subname//" ERROR: bad value for num_actunit_B" )
          return
       end if
       if(.not. present(filter_actunit_B))then
          write(iulog,*) "error: num_actunit_B is presented but filter_actunit_B is missing"
          call endrun( subname//" ERROR: missing required optional arguments" )
          return
       end if
       SHR_ASSERT_FL((size(filter_actunit_B) >= num_actunit_B), sourcefile, __LINE__)
    end if
    if(present(num_actunit_C))then
       if(num_actunit_C < 0)then
          write(iulog,*) "error: num_actunit_C cannot be less than 0"
          call endrun( subname//" ERROR: bad value for num_actunit_C" )
          return
       end if
       if(.not. present(filter_actunit_C))then
          write(iulog,*) "error: num_actunit_C is presented but filter_actunit_C is missing"
          call endrun( subname//" ERROR: missing required optional arguments" )
          return
       end if
       SHR_ASSERT_FL((size(filter_actunit_C) >= num_actunit_C), sourcefile, __LINE__)
    end if

    if(.not. list_ready)then
       i_a=1
       i_b=1
       i_c=1
       i_abc=1
       Aindex(1:A%NE)    = (A%CI(1:A%NE)-1)*A%SM+A%RI(1:A%NE)
       Bindex(1:B%NE)    = (B%CI(1:B%NE)-1)*B%SM+B%RI(1:B%NE)
       Cindex(1:C%NE)    = (C%CI(1:C%NE)-1)*C%SM+C%RI(1:C%NE)
       Aindex(A%NE+1)    = A%SM*A%SM+1
       Bindex(B%NE+1)    = B%SM*B%SM+1
       Cindex(C%NE+1)    = C%SM*C%SM+1

       do while (i_a .le. A%NE .or. i_b .le. B%NE .or. i_c .le. C%NE)
          if(Aindex(i_a) .lt. Bindex(i_b) .and. Aindex(i_a) .lt. Cindex(i_c))then
             do fu = 1, num_unit
                u = filter_u(fu)
                this%M(u,i_abc) = A%M(u,i_a)
             end do
             ABCindex(i_abc) = Aindex(i_a)
             if(present(list_A))list_A(i_a) = i_abc
             i_a   = i_a  + 1
             i_abc = i_abc + 1
          else
             if(Bindex(i_b) .lt. Aindex(i_a) .and. Bindex(i_b) .lt. Cindex(i_c))then 
                do fu = 1, num_unit
                   u = filter_u(fu)
                   this%M(u,i_abc) = B%M(u,i_b)
                end do
                ABCindex(i_abc) = Bindex(i_b)
                if(present(list_B))list_B(i_b) = i_abc
                i_b   = i_b  + 1
                i_abc = i_abc + 1
             else
                if(Cindex(i_c) .lt. Aindex(i_a) .and. Cindex(i_c) .lt. Bindex(i_b))then
                   do fu = 1, num_unit
                      u = filter_u(fu)
                      this%M(u,i_abc) = C%M(u,i_c)
                   end do
                   ABCindex(i_abc) = Cindex(i_c)
                   if(present(list_C))list_C(i_c) = i_abc
                   i_c   = i_c  + 1
                   i_abc = i_abc + 1
                else
                   if(Aindex(i_a) .eq. Bindex(i_b) .and. Aindex(i_a) .lt. Cindex(i_c))then
                      do fu = 1, num_unit
                         u = filter_u(fu)
                         this%M(u,i_abc) = A%M(u,i_a) + B%M(u,i_b)
                      end do
                      ABCindex(i_abc) = Aindex(i_a)
                      if(present(list_A))list_A(i_a) = i_abc
                      if(present(list_B))list_B(i_b) = i_abc
                      i_a   = i_a  + 1
                      i_b   = i_b  + 1
                      i_abc = i_abc + 1
                   else
                      if(Aindex(i_a) .eq. Cindex(i_c) .and. Aindex(i_a) .lt. Bindex(i_b))then
                         do fu = 1, num_unit
                            u = filter_u(fu)
                            this%M(u,i_abc) = A%M(u,i_a) + C%M(u,i_c)
                         end do
                         ABCindex(i_abc) = Aindex(i_a)
                         if(present(list_A))list_A(i_a) = i_abc
                         if(present(list_C))list_C(i_c) = i_abc
                         i_a   = i_a  + 1
                         i_c   = i_c  + 1
                         i_abc = i_abc + 1
                      else
                         if(Bindex(i_b) .eq. Cindex(i_c) .and. Bindex(i_b) .lt. Aindex(i_a))then
                            do fu = 1, num_unit
                               u = filter_u(fu)
                               this%M(u,i_abc) = B%M(u,i_b) + C%M(u,i_c)
                            end do
                            ABCindex(i_abc) = Bindex(i_b)
                            if(present(list_B))list_B(i_b) = i_abc
                            if(present(list_C))list_C(i_c) = i_abc
                            i_b   = i_b  + 1
                            i_c   = i_c  + 1
                            i_abc = i_abc + 1
                         else
                            if(Aindex(i_a) .eq. Bindex(i_b) .and. Aindex(i_a) .eq. Cindex(i_c))then
                               do fu = 1, num_unit
                                  u = filter_u(fu)
                                  this%M(u,i_abc) = A%M(u,i_a) + B%M(u,i_b) + C%M(u,i_c)
                               end do
                               ABCindex(i_abc) = Bindex(i_b)
                               if(present(list_A))list_A(i_a) = i_abc
                               if(present(list_B))list_B(i_b) = i_abc
                               if(present(list_C))list_C(i_c) = i_abc
                               i_a   = i_a  + 1
                               i_b   = i_b  + 1
                               i_c   = i_c  + 1
                               i_abc = i_abc + 1
                            else
                               write(iulog,*) 'Error in subroutine SPMP_ABC',Aindex(i_a),Bindex(i_b),Cindex(i_c)
                            end if
                         end if
                      end if
                   end if
                end if
             end if
          end if
       end do

       this%NE = i_abc - 1
       this%CI(1:this%NE) = (ABCindex(1:this%NE) - 1) / this%SM + 1
       this%RI(1:this%NE) = ABCindex(1:this%NE) - this%SM * (this%CI(1:this%NE) - 1)
       if(present(NE_ABC))NE_ABC = this%NE
       if(present(CI_ABC))CI_ABC(1:this%NE) = this%CI(1:this%NE)
       if(present(RI_ABC))RI_ABC(1:this%NE) = this%RI(1:this%NE)
       if(present(list_A) .and. present(list_B) .and. present(list_C) .and. present(NE_ABC) .and. present(RI_ABC) .and. present(CI_ABC))list_ready = .true.
    else
       do i = 1, NE_ABC
          do fu = 1, num_unit
             u = filter_u(fu)
             this%M(u,i) = 0._r8
          end do
       end do
       if(present(num_actunit_A))then
          do i_a = 1, A%NE
             do fu = 1, num_actunit_A
                u = filter_actunit_A(fu)
                this%M(u,list_A(i_a)) = A%M(u,i_a)
             end do
          end do
       else
          do i_a = 1, A%NE
             do fu = 1, num_unit
                u = filter_u(fu)
                this%M(u,list_A(i_a)) = A%M(u,i_a)
             end do
          end do
       end if
       if(present(num_actunit_B))then
          do i_b = 1, B%NE
             do fu = 1, num_actunit_B
                u = filter_actunit_B(fu)
                this%M(u,list_B(i_b)) = this%M(u,list_B(i_b)) + B%M(u,i_b)
             end do
          end do
       else
          do i_b = 1, B%NE
             do fu = 1, num_unit
                u = filter_u(fu)
                this%M(u,list_B(i_b)) = this%M(u,list_B(i_b)) + B%M(u,i_b)
             end do
          end do
       end if
       if(present(num_actunit_C))then
          do i_c = 1, C%NE
             do fu = 1, num_actunit_C
                u = filter_actunit_C(fu)
                this%M(u,list_C(i_c)) = this%M(u,list_C(i_c)) + C%M(u,i_c)
             end do
          end do
       else
          do i_c = 1, C%NE
             do fu = 1, num_unit
                u = filter_u(fu)
                this%M(u,list_C(i_c)) = this%M(u,list_C(i_c)) + C%M(u,i_c)
             end do
          end do
       end if
       this%NE = NE_ABC
       this%CI(1:this%NE) = CI_ABC(1:NE_ABC)
       this%RI(1:this%NE) = RI_ABC(1:NE_ABC)
    end if

  end subroutine SPMP_ABC

  !-----------------------------------------------------------------------

end module SparseMatrixMultiplyMod
