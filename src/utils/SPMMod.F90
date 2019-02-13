module SPMMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SPMMod
!
! !DESCRIPTION:
! Sparse matrix multiplication add addition
!
! Author: Xingjie Lu
!
!EOP
!-----------------------------------------------------------------------

  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use decompMod                      , only : bounds_type
  implicit none
  private

  type, public :: sparse_matrix_type
 
  !sparse matrix is in COO format, Both row index and column index should be in ascending order. 
  !Row index should change faster than Column index to ensure SPMP_AB work properly.

     real(r8), pointer :: M(:,:)
     integer , pointer :: RI(:) !Row index
     integer , pointer :: CI(:) !Column index
     integer NE
     integer SM
     integer num_unit
     integer begu
     integer endu
  
  contains
    
    procedure, public :: InitSM
    procedure, public :: SetValueSM
    procedure, public :: SetValueA
    procedure, public :: SetValueA_diag
    procedure, public :: SPMM_AK
    procedure, public :: SPMP_AB   
    procedure, public :: SPMP_B_ACC
    procedure, public :: SPMP_ABC

  end type sparse_matrix_type

  type, public :: diag_matrix_type
  
  !diagnoal matrix only store diagnoal entries

     real(r8), pointer :: DM(:,:)
     integer SM
     integer num_unit
     integer begu
     integer endu
  
  contains

    procedure, public :: InitDM
    procedure, public :: SetValueDM

  end type diag_matrix_type

  type, public :: vector_type
  
  !vector 

     real(r8), pointer :: V(:,:)
     integer SV
     integer num_unit
     integer begu
     integer endu

  contains 
  
    procedure, public :: InitV
    procedure, public :: ReleaseV
    procedure, public :: SetValueV
    procedure, public :: SetValueV_scaler
    procedure, public :: SPMM_AX

  end type vector_type

contains

subroutine InitSM(this,SM_in,begu,endu,maxsm)

class(sparse_matrix_type) :: this
integer,intent(in) :: SM_in
integer,intent(in) :: begu
integer,intent(in) :: endu
integer,optional,intent(in) :: maxsm

this%SM = SM_in
this%begu = begu
this%endu = endu
if(present(maxsm))then
   allocate(this%M(begu:endu,1:maxsm))
else
   allocate(this%M(begu:endu,1:SM_in*SM_in))
end if
allocate(this%RI(1:SM_in*SM_in))
allocate(this%CI(1:SM_in*SM_in))
this%M(:,:) = -9999._r8
this%RI(:) = -9999
this%CI(:) = -9999
this%NE    = -9999

end subroutine InitSM

subroutine SetValueSM(this,begu,endu,num_unit,filter_u,M,I,J,NE_in)!,index_ready)

class(sparse_matrix_type) :: this
integer,intent(in) :: begu
integer,intent(in) :: endu
integer,intent(in) :: NE_in
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)
real(r8),dimension(:,:),intent(in) :: M(begu:endu,1:NE_in)
integer ,dimension(:),intent(in) :: I(1:NE_in)
integer ,dimension(:),intent(in) :: J(1:NE_in)
!logical,optional,intent(inout) :: index_ready

integer k,u,fu

!print*,'NE_in',NE_in
!print*,'M',M(1:NE_in)
do k = 1,NE_in
   do fu = 1,num_unit
      u = filter_u(fu)  
      this%M(u,k) = M(u,k)
   end do
end do

this%NE = NE_in
!if(.not. present(index_ready) .or. .not. index_ready)then
   do k = 1,NE_in
      this%RI(k) = I(k)
      this%CI(k) = J(k)
   end do
!   if(present(index_ready))index_ready = .true.
!end if

end subroutine SetValueSM

subroutine SetValueA_diag(this,num_unit,filter_u,scaler)

class(sparse_matrix_type) :: this
real(r8),intent(in) :: scaler
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)
integer i,u,fu

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

subroutine SetValueA(this,begu,endu,num_unit,filter_u,M,AI,AJ,NE_NON,Init_ready,list,RI_A,CI_A)

class(sparse_matrix_type) :: this
integer,intent(in) :: begu
integer,intent(in) :: endu
integer ,intent(in) :: NE_NON
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)
real(r8),dimension(:),intent(in) :: M(begu:endu,1:NE_NON)
integer ,dimension(:),intent(in) :: AI(1:NE_NON)
integer ,dimension(:),intent(in) :: AJ(1:NE_NON)
logical ,intent(inout) :: Init_ready !True: diagnoal of A has been set to -1,this%RI, this%CI, this%NE and list has been set up
integer ,dimension(:),intent(inout),optional :: list(1:NE_NON)
integer ,dimension(:),intent(inout),optional :: RI_A(1:NE_NON+this%SM)
integer ,dimension(:),intent(inout),optional :: CI_A(1:NE_NON+this%SM)

integer i,j,k,fu,u
logical list_ready
type(sparse_matrix_type) :: A_diag, A_nondiag
!logical :: DiagnolSet
list_ready = .false.
!DiagnolSet = .False.

!    if(NE_NON .eq. 200)then
!       do j=1,NE_NON
!          do fu = 1,num_unit
!             u = filter_u(fu)
!             print*,'begin setvalueSM',u,AI(j),AJ(j),j,M(u,j),num_unit,this%begu,this%endu,this%SM,NE_NON
!          end do
!       end do
!    end if
if(init_ready .and. .not. (present(list) .and. present(RI_A) .and. present(CI_A)))then
   write(*,*),"Error: initialization is ready, but at least one of list, RI_A or CI_A is not presented"
end if
 !i=1
 !k=1
 
 !if(.not. list_ready)then
 !   do j=1,NE_NON
 !      do while(i .lt. AJ(j))
 !         if(.not. DiagnolSet)then
 !            do fu = 1, num_unit
 !               u = filter_u(fu)
 !               this%M(u,k) = -1._r8
 !            end do
 !            this%RI(k) = i
 !            this%CI(k) = i
 !            k = k + 1
 !         end if
 !         i = i + 1
 !         DiagnolSet = .False.
 !      !print*,'first loop',i,AJ(j),k
!!       print*,this%M(1:k-1),this%RI(1:k-1),this%CI(1:k-1)
!       end do
!    print*,'end of loop'
!       if(AI(j) .lt. i)then  !AJ(j) .eq. i)
!       print*,'upper diagnoal'
!          do fu = 1, num_unit
!             u = filter_u(fu)
!             this%M(u,k) = M(j)
!          end do
!          list(j)    = k
!          this%RI(k) = AI(j)
!          this%CI(k) = AJ(j)
!          k = k + 1
!       else
!          if(DiagnolSet)then
!!          print*,'lower diagnoal1'
!             do fu = 1, num_unit
!                u = filter_u(fu)
!                this%M(u,k) = M(j)
!             end do
!             list(j)    = k
!             this%RI(k) = AI(j)
!             this%CI(k) = AJ(j)
!             k = k + 1
!           print*,this%M(1:k-1),this%RI(1:k-1),this%CI(1:k-1)
!          else
!          print*,'lower diagnoal2'
!             do fu = 1, num_unit
!                u = filter_u(fu)
!                this%M(u,k) = -1._r8
!             end do
!             this%RI(k) = i
!             this%CI(k) = i
!             k = k + 1
!          print*,this%M(1:k-1),this%RI(1:k-1),this%CI(1:k-1)
!             DiagnolSet = .True.
!             this%M(k) = M(j)
!             list(j)    = k
!             this%RI(k) = AI(j)
!             this%CI(k) = AJ(j)
!             k = k + 1
!          print*,this%M(1:k-1),this%RI(1:k-1),this%CI(1:k-1)
!          end if
!       end if
!    end do
 
!    do while(i .le. this%SM)
!!    print*,'end loop',i,this%SM,DiagnolSet
!       if(.not. DiagnolSet)then
!          this%M(k) = -1._r8
!          this%RI(k) = i
!          this%CI(k) = i
!          k = k + 1
!          !print*,this%M(1:k-1),this%RI(1:k-1),this%CI(1:k-1)
!       end if
!       i = i + 1
!       DiagnolSet = .False.
!    end do
!    this%NE = k - 1
!    list_ready = .True.
! end if
! print*,'init_ready',Init_ready 
 if(Init_ready)then
!    print*,'this%NE',this%NE
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
!    print*,'setvaluA,ri_2nd',this%NE,this%RI
!    print*,'setvaluA,ci_2nd',this%CI
 else
!    print*,'sm,begu,endu',this%SM,this%begu,this%endu
!    call A%initSM(this%SM,this%begu,this%endu)
    call A_diag%InitSM(this%SM,begu,endu)
    call A_nondiag%InitSM(this%SM,begu,endu)

!    print*,'AI',AI
!    print*,'AJ',AJ
!    print*,'M',NE_NON,M(this%begu:this%endu,1:NE_NON)
    call A_diag%SetValueA_diag(num_unit,filter_u,-1._r8)
!    if(NE_NON .eq. 200)then
!       do j=1,NE_NON
!          do fu = 1,num_unit
!             u = filter_u(fu)
!             print*,'before setvalueSM',u,AI(j),AJ(j),j,M(u,j)
!          end do
!       end do
!    end if
    call A_nondiag%SetValueSM(begu,endu,num_unit,filter_u,M,AI,AJ,NE_NON)
!    if(NE_NON .eq. 200)then
!       do j=1,NE_NON
!          do fu = 1,num_unit
!             u = filter_u(fu)
!             print*,'after setvalueSM',u,A_nondiag%RI(j),A_nondiag%CI(j),j,A_nondiag%M(u,j)
!          end do
!       end do
!    end if
!    print*,'A_nondiag',A_nondiag%M(5,:)
!    do fu = 1,num_unit
!       u = filter_u(fu)
!       if(A_diag%M(u,1) .ge. -100)print*,'A_diag',A_diag%M(u,:)
!    end do
!    print*,'A_diagRI',A_diag%RI
!    print*,'A_diagCI',A_diag%CI

!    print*,'present_list',present(list),list_ready
    if(present(list))then
       call this%SPMP_AB(num_unit,filter_u,A_nondiag,A_diag,list_ready,list_A=list)
    else
       call this%SPMP_AB(num_unit,filter_u,A_nondiag,A_diag,list_ready)
    end if
!    if(NE_NON .eq. 200)then
!       do j=1,this%NE
!          do fu = 1,num_unit
!             u = filter_u(fu)
!             print*,'after SPMP_AB',u,this%RI(j),this%CI(j),j
!             print*,'M',this%M(u,j)
!          end do
!       end do
!    end if
!    if(present(list))print*,'after spmp',present(list),list_ready,list
!    do fu = 1,num_unit
!       u = filter_u(fu)
!       if(A%M(u,1) .ge. -100)print*,'A%M',A%M(u,:)
!    end do
!    print*,'A%RI',A%RI
!    print*,'A%CI',A%CI
!    print*,'NE',A%NE
!    call this%SetValueSM(num_unit,filter_u,A_nondiag%M,A_nondiag%RI,A_nondiag%CI,A_nondiag%NE)
!    print*,'A_nondiagRI',A_nondiag%NE,A_nondiag%RI
!    print*,'A_nondiagCI',A_nondiag%CI
    if(present(RI_A))RI_A(1:this%NE) = this%RI(1:this%NE)
    if(present(CI_A))CI_A(1:this%NE) = this%CI(1:this%NE)
!    print*,'setvaluA,ri_1st',this%RI
!    print*,'setvaluA,ci_1st',this%CI
    
    Init_ready = .true.
 end if
! print*,'end SetValueA'        
end subroutine SetValueA

subroutine InitDM(this,SM_in,begu,endu)

class(diag_matrix_type) :: this
integer,intent(in) :: SM_in
integer,intent(in) :: begu
integer,intent(in) :: endu


this%SM = SM_in
allocate(this%DM(begu:endu,1:SM_in))
this%DM(:,:) = -9999._r8
this%begu = begu
this%endu = endu

end subroutine InitDM

subroutine SetValueDM(this,begu,endu,num_unit,filter_u,M)

class(diag_matrix_type) :: this
integer,intent(in) :: begu
integer,intent(in) :: endu
real(r8),dimension(:),intent(in) :: M(begu:endu,1:this%SM)
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)

integer i,fu,u
!print*,'in SetValuDM',begu,endu
do i = 1,this%SM
   do fu = 1,num_unit
      u = filter_u(fu)
!      if(begu .le. 3998 .and. endu .ge. 3998)print*,'SetValueDM',u,i,M(u,i),this%SM
      this%DM(u,i) = M(u,i)
   end do
end do

end subroutine SetValueDM

subroutine InitV(this,SV_in,begu,endu)

class(vector_type) :: this
integer,intent(in) :: SV_in
integer,intent(in) :: begu
integer,intent(in) :: endu

this%SV = SV_in
allocate(this%V(begu:endu,1:SV_in))
this%V(:,:) = -9999._r8
this%begu = begu
this%endu = endu

end subroutine InitV

subroutine ReleaseV(this)

class(vector_type) :: this

deallocate(this%V)

end subroutine ReleaseV

subroutine SetValueV_scaler(this,num_unit,filter_u,scaler)

class(vector_type) :: this
real(r8),intent(in) :: scaler
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)

integer i,fu,u

do i=1,this%SV
   do fu = 1,num_unit
      u = filter_u(fu)
      this%V(u,i) = scaler
   end do
end do

end subroutine SetValueV_scaler

subroutine SetValueV(this,begu,endu,num_unit,filter_u,M)

integer,intent(in) :: begu
integer,intent(in) :: endu
class(vector_type) :: this
real(r8),dimension(:,:),intent(in) :: M(begu:endu,1:this%SV)
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)

integer i,fu,u

do i=1,this%SV
   do fu = 1,num_unit
      u = filter_u(fu)
      this%V(u,i) = M(u,i)
   end do
end do

end subroutine SetValueV

subroutine SPMM_AK(this,num_unit,filter_u,K)!find_nonzero)

! A is a sparse matrix in Coordinate format (COO). 
! K is a diagnoal matrix in one dimension array.

class(sparse_matrix_type)  :: this
type(diag_matrix_type)  ,intent(in)  :: K
!logical,intent(in) :: index_ready
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)
!logical,optional,intent(in) :: find_nonzero

integer i,fu,u

!if(find_nonzero)then
!   if(index_ready)then
!      AK%M(1:AK%NE) = 0._r8
!      do i=1,A%NE
!         if(A%M(i) .ne. 0 .and. K%DM(A%CI(i)) .ne. 0)then
!            AK%M(i) = A%M(i) * K%DM(A%CI(i))
!         end if
!      end do
!   else
!      AK%M = -9999._r8
!      do i=1,A%NE
!         if(A%M(i) .ne. 0 .and. K%DM(A%CI(i)) .ne. 0)then
!            AK%M(i) = A%M(i) * K%DM(A%CI(i))
!         else
!            AK%M(i) = 0._r8
!         end if
!      end do
!   end if
!else
!AK%RI(1:A%NE) = A%RI(1:A%NE)
!AK%CI(1:A%NE) = A%CI(1:A%NE)
!AK%NE = A%NE

!if(present(find_nonzero) .and. find_nonzero)then
!   do i=1,AK%NE
!      do fu = 1,num_unit
!         u = filter_u(fu)
!         AK%M(u,i) = 0._r8
!      end do
!   end do
      
!   do i=1,A%NE
!      do fu = 1,num_unit
!         u = filter_u(fu)
!         if(A%M(u,i) .ne. 0 .and. K%DM(u,A%CI(i)) .ne. 0)then
!            AK%M(u,i) = A%M(u,i) * K%DM(u,A%CI(i))
!         end if
!      end do
!   end do
!else
   do i=1,this%NE
      do fu = 1,num_unit
         u = filter_u(fu)
         this%M(u,i) = this%M(u,i) * K%DM(u,this%CI(i))
      end do
   end do
!end if

!   if(.not. index_ready)then
!   end if

end subroutine SPMM_AK

subroutine SPMM_AX(this,num_unit,filter_u,A)

! A is a sparse matrix in Coordinate format (COO). 
! K is a diagnoal matrix in one dimension array.

class(vector_type)         :: this
type(sparse_matrix_type),intent(in)    :: A
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)

integer i,fu,u
real(r8),dimension(:,:) :: V(this%begu:this%endu,1:this%SV)

do i=1,this%SV
   do fu = 1, num_unit
      u = filter_u(fu)
      V(u,i) = this%V(u,i)
   end do
end do

do i=1,A%NE
!   if(A%M(i) .ne. 0 .and. X%V(A%CI(i)) .ne. 0)then
   do fu = 1, num_unit
      u = filter_u(fu)
!      print*,'in SPMMAX',i,u,A%RI(i),A%CI(i),A%M(u,i),this%V(u,A%CI(i))
      this%V(u,A%RI(i)) = this%V(u,A%RI(i)) + A%M(u,i) * V(u,A%CI(i))
   end do
!   end if
   !print*,'in SPMMAX1',AX%V(1:3)
end do

end subroutine SPMM_AX

subroutine SPMP_B_ACC(this,num_unit,filter_u,A)

!A is a sparse matrix in Coordinate format (COO)

class(sparse_matrix_type)     :: this
type(sparse_matrix_type),intent(in)     :: A
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)

integer i,fu,u

this%NE=A%NE
this%RI=A%RI
this%CI=A%CI

do i=1,A%NE
   do fu = 1, num_unit
      u = filter_u(fu)
      this%M(u,i) = this%M(u,i) + A%M(u,i)
   end do
end do   

end subroutine SPMP_B_ACC

subroutine SPMP_AB(this,num_unit,filter_u,A,B,list_ready,list_A,list_B,NE_AB,RI_AB,CI_AB)

!A is a sparse matrix in Coordinate format (COO)
!B is a sparse matrix in Coordinate format (COO)
!AB is a sparse matrix in Coordinate format (COO)

class(sparse_matrix_type)     :: this
type(sparse_matrix_type),intent(in)     :: A
type(sparse_matrix_type),intent(in)     :: B
logical,intent(inout)                   :: list_ready 
integer,intent(in) :: num_unit
integer,intent(in) :: filter_u(:)

!'list_ready = .true.' means AB%NE, AB%RI, AB%CI, list_A, list_B are previously assigned.
!Otherwise, use 'list_ready = .false.' to get those information.
integer,dimension(1:A%NE),intent(inout),optional :: list_A
integer,dimension(1:B%NE),intent(inout),optional :: list_B
integer,intent(inout),optional :: NE_AB
integer,dimension(1:A%NE+B%NE),intent(inout),optional :: RI_AB
integer,dimension(1:A%NE+B%NE),intent(inout),optional :: CI_AB

integer,dimension(:) :: Aindex(A%NE+1),Bindex(B%NE+1)
integer,dimension(:) :: ABindex(this%SM*this%SM)

!real(r8),dimension(this%begu:this%endu,this%NE) :: A_M

integer i_a,i_b,i_ab
integer i,fu,u

if(list_ready .and. .not. (present(list_A) .and. present(list_B) .and. present(NE_AB) .and. present(RI_AB) .and. present(CI_AB)))then
   write(*,*),"error in SPMP_AB: list_ready is True, but at least one of list_A, list_B, NE_AB, RI_AB and CI_AB are not presented"
end if
!AB%M = -9999._r8

!do i=1,this%NE
!   do fu = 1,num_unit
!      u = filter_u(fu)
!      A_M(u,i) = this%M(u,i)
!   end do
!end do

if(.not. list_ready)then
!   AB%CI = -9999
!   AB%RI = -9999
   i_a=1
   i_b=1
   i_ab=1
   Aindex(1:A%NE)    = (A%CI-1)*A%SM+A%RI
   Bindex(1:B%NE)    = (B%CI-1)*B%SM+B%RI
   Aindex(A%NE+1)    = A%SM*A%SM+1
   Bindex(B%NE+1)    = B%SM*B%SM+1
!print*,'Aindex',A%NE,A%SM,A%CI,A%RI,this%NE
!print*,'Bindex',B%NE,B%SM,B%CI,B%RI

   do while (i_a .le. A%NE .or. i_b .le. B%NE)
      if(Aindex(i_a) .lt. Bindex(i_b))then
!      print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
         do fu = 1, num_unit
            u = filter_u(fu)
            this%M(u,i_ab) = A%M(u,i_a)
         end do
         ABindex(i_ab) = Aindex(i_a)
         if(present(list_A))list_A(i_a) = i_ab
!      print*,'ABindex',ABindex(i_ab)
         i_a  = i_a  + 1
         i_ab = i_ab + 1
      else
         if(Aindex(i_a) .gt. Bindex(i_b))then 
!         print*,'here2'
!         print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
            do fu = 1, num_unit
               u = filter_u(fu)
               this%M(u,i_ab) = B%M(u,i_b)
            end do
            ABindex(i_ab) = Bindex(i_b)
            if(present(list_B))list_B(i_b) = i_ab
!         print*,'ABindex',ABindex(i_ab)
            i_b  = i_b  + 1
            i_ab = i_ab + 1
         else
!         print*,'here3'
!         print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
            do fu = 1, num_unit
               u = filter_u(fu)
               this%M(u,i_ab) = A%M(u,i_a) + B%M(u,i_b)
            end do
            ABindex(i_ab) = Aindex(i_a)
!         print*,'ABindex',ABindex(i_ab)
            if(present(list_A))list_A(i_a) = i_ab
            if(present(list_B))list_B(i_b) = i_ab
            i_a  = i_a  + 1
            i_b  = i_b  + 1
            i_ab = i_ab + 1
         end if
      end if
   end do

   this%NE = i_ab - 1
!   print*,'AB,this%NE',this%NE
   this%CI(1:this%NE) = (ABindex(1:this%NE) - 1) / this%SM + 1
!   print*,'AB,this%CI',this%CI(1:this%NE)
   this%RI(1:this%NE) = ABindex(1:this%NE) - this%SM * (this%CI(1:this%NE) - 1)
!   print*,'AB,this%RI',this%RI(1:this%NE)
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
!   print*,'A%NE',this%NE,B%NE,list_A(:),list_B(:)
   do i_a = 1, A%NE
      !print*,'i_a',i_a
!      if(A%M(i_a) .ne. 0)then
      do fu = 1, num_unit
         u = filter_u(fu)
!         print*,'in SPMP_AB',u,i_a,list_A(i_a)
         this%M(u,list_A(i_a)) = A%M(u,i_a)
      end do
!      end if
   end do
   do i_b = 1, B%NE
!      if(B%M(i_b) .ne. 0)then
      do fu = 1, num_unit
         u = filter_u(fu)
         this%M(u,list_B(i_b)) = this%M(u,list_B(i_b)) + B%M(u,i_b)
      end do
!      end if
   end do
   this%NE = NE_AB
   this%CI(1:this%NE) = CI_AB(1:NE_AB)
   this%RI(1:this%NE) = RI_AB(1:NE_AB)
end if

end subroutine SPMP_AB

subroutine SPMP_ABC(this,num_unit,filter_u,A,B,C,list_ready,list_A,list_B,list_C,NE_ABC,RI_ABC,CI_ABC,&
           use_actunit_list_A,num_actunit_A,filter_actunit_A,use_actunit_list_B,num_actunit_B,filter_actunit_B,&
           use_actunit_list_C,num_actunit_C,filter_actunit_C)

!A is a sparse matrix in Coordinate format (COO)
!B is a sparse matrix in Coordinate format (COO)
!AB is a sparse matrix in Coordinate format (COO)

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


!'list_ready = .true.' means AB%NE, AB%RI, AB%CI, list_A, list_B are previously assigned.
!Otherwise, use 'list_ready = .false.' to get those information.
integer,dimension(1:A%NE),intent(inout),optional :: list_A
integer,dimension(1:B%NE),intent(inout),optional :: list_B
integer,dimension(1:C%NE),intent(inout),optional :: list_C
integer,intent(inout),optional :: NE_ABC
integer,dimension(1:A%NE+B%NE+C%NE),intent(inout),optional :: RI_ABC
integer,dimension(1:A%NE+B%NE+C%NE),intent(inout),optional :: CI_ABC

integer,dimension(:) :: Aindex(A%NE+1),Bindex(B%NE+1),Cindex(C%NE+1)
integer,dimension(:) :: ABCindex(this%SM*this%SM)

!real(r8),dimension(this%begu:this%endu,this%NE) :: A_M

integer i_a,i_b,i_c,i_abc
integer i,fu,u
!print*,'SPMP_AB',A%NE,B%NE,AB%SM
!print*,'A%RI',A%RI(1:10)
!print*,'A%CI',A%CI(1:10)
!print*,'B%RI',B%RI(1:10)
!print*,'B%CI',B%CI(1:10)

if(list_ready .and. .not. (present(list_A) .and. present(list_B) .and. present(list_C) .and. present(NE_ABC) .and. present(RI_ABC) .and. present(CI_ABC)))then
   write(*,*),"error in SPMP_ABC: list_ready is True, but at least one of list_A, list_B, list_C, NE_ABC, RI_ABC and CI_ABC are not presented",&
              present(list_A),present(list_B),present(list_C),present(NE_ABC),present(RI_ABC),present(CI_ABC)
end if
if(present(num_actunit_A))then
   if(num_actunit_A .eq. 0)then
      write(*,*),"error: num_actunit_A cannot be set to 0"
   end if
   if(.not. present(filter_actunit_A))then
      write(*,*),"error: num_actunit_A is presented but filter_actunit_A is missing"
   end if
end if
if(present(num_actunit_B))then
   if(num_actunit_B .eq. 0)then
      write(*,*),"error: num_actunit_B cannot be set to 0"
   end if
   if(.not. present(filter_actunit_B))then
      write(*,*),"error: num_actunit_B is presented but filter_actunit_B is missing"
   end if
end if
if(present(num_actunit_C))then
   if(num_actunit_C .eq. 0)then
      write(*,*),"error: num_actunit_C cannot be set to 0"
   end if
   if(.not. present(filter_actunit_C))then
      write(*,*),"error: num_actunit_C is presented but filter_actunit_C is missing"
   end if
end if
!print*,'here'
!AB%M = -9999._r8

!do i=1,this%NE
!   do fu = 1,num_unit
!      u = filter_u(fu)
!      A_M(u,i) = this%M(u,i)
!   end do
!end do

if(.not. list_ready)then
!   AB%CI = -9999
!   AB%RI = -9999
   i_a=1
   i_b=1
   i_c=1
   i_abc=1
   Aindex(1:A%NE)    = (A%CI-1)*A%SM+A%RI
   Bindex(1:B%NE)    = (B%CI-1)*B%SM+B%RI
   Cindex(1:C%NE)    = (C%CI-1)*C%SM+C%RI
   Aindex(A%NE+1)    = A%SM*A%SM+1
   Bindex(B%NE+1)    = B%SM*B%SM+1
   Cindex(C%NE+1)    = C%SM*C%SM+1
!print*,'Aindex',A%NE,A%SM,A%CI,A%RI,AB%NE
!print*,'Bindex',B%NE,B%SM,B%CI,B%RI

   do while (i_a .le. A%NE .or. i_b .le. B%NE .or. i_c .le. C%NE)
!      if(Aindex(i_a) .lt. Bindex(i_b))then
      if(Aindex(i_a) .lt. Bindex(i_b) .and. Aindex(i_a) .lt. Cindex(i_c))then
!      print*,'here1'
!      print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
         do fu = 1, num_unit
            u = filter_u(fu)
            this%M(u,i_abc) = A%M(u,i_a)
         end do
         ABCindex(i_abc) = Aindex(i_a)
         if(present(list_A))list_A(i_a) = i_abc
!      print*,'ABindex',ABindex(i_ab)
         i_a   = i_a  + 1
         i_abc = i_abc + 1
      else
!         if(Aindex(i_a) .gt. Bindex(i_b))then 
         if(Bindex(i_b) .lt. Aindex(i_a) .and. Bindex(i_b) .lt. Cindex(i_c))then 
!         print*,'here2'
!         print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
            do fu = 1, num_unit
               u = filter_u(fu)
               this%M(u,i_abc) = B%M(u,i_b)
            end do
            ABCindex(i_abc) = Bindex(i_b)
            if(present(list_B))list_B(i_b) = i_abc
!         print*,'ABindex',ABindex(i_ab)
            i_b   = i_b  + 1
            i_abc = i_abc + 1
         else
            if(Cindex(i_c) .lt. Aindex(i_a) .and. Cindex(i_c) .lt. Bindex(i_b))then
!         print*,'here2'
!         print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
               do fu = 1, num_unit
                  u = filter_u(fu)
                  this%M(u,i_abc) = C%M(u,i_c)
               end do
               ABCindex(i_abc) = Cindex(i_c)
               if(present(list_C))list_C(i_c) = i_abc
!         print*,'ABindex',ABindex(i_ab)
               i_c   = i_c  + 1
               i_abc = i_abc + 1
            else
               if(Aindex(i_a) .eq. Bindex(i_b) .and. Aindex(i_a) .lt. Cindex(i_c))then
!         print*,'here3'
!         print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
                  do fu = 1, num_unit
                     u = filter_u(fu)
                     this%M(u,i_abc) = A%M(u,i_a) + B%M(u,i_b)
                  end do
                  ABCindex(i_abc) = Aindex(i_a)
!         print*,'ABindex',ABindex(i_ab)
                  if(present(list_A))list_A(i_a) = i_abc
                  if(present(list_B))list_B(i_b) = i_abc
                  i_a   = i_a  + 1
                  i_b   = i_b  + 1
                  i_abc = i_abc + 1
               else
                  if(Aindex(i_a) .eq. Cindex(i_c) .and. Aindex(i_a) .lt. Bindex(i_b))then
!         print*,'here3'
!         print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
                     do fu = 1, num_unit
                        u = filter_u(fu)
                        this%M(u,i_abc) = A%M(u,i_a) + C%M(u,i_c)
                     end do
                     ABCindex(i_abc) = Aindex(i_a)
!         print*,'ABindex',ABindex(i_ab)
                     if(present(list_A))list_A(i_a) = i_abc
                     if(present(list_C))list_C(i_c) = i_abc
                     i_a   = i_a  + 1
                     i_c   = i_c  + 1
                     i_abc = i_abc + 1
                  else
                     if(Bindex(i_b) .eq. Cindex(i_c) .and. Bindex(i_b) .lt. Aindex(i_a))then
!         print*,'here3'
!         print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
                        do fu = 1, num_unit
                           u = filter_u(fu)
                           this%M(u,i_abc) = B%M(u,i_b) + C%M(u,i_c)
                        end do
                        ABCindex(i_abc) = Bindex(i_b)
!         print*,'ABindex',ABindex(i_ab)
                        if(present(list_B))list_B(i_b) = i_abc
                        if(present(list_C))list_C(i_c) = i_abc
                        i_b   = i_b  + 1
                        i_c   = i_c  + 1
                        i_abc = i_abc + 1
                     else
                        if(Aindex(i_a) .eq. Bindex(i_b) .and. Aindex(i_a) .eq. Cindex(i_c))then
!         print*,'here3'
!         print*,'i_a,i_b,i_ab',i_a,i_b,i_ab,Aindex(i_a),Bindex(i_b)
                           do fu = 1, num_unit
                              u = filter_u(fu)
                              this%M(u,i_abc) = A%M(u,i_a) + B%M(u,i_b) + C%M(u,i_c)
                           end do
                           ABCindex(i_abc) = Bindex(i_b)
!         print*,'ABindex',ABindex(i_ab)
                           if(present(list_A))list_A(i_a) = i_abc
                           if(present(list_B))list_B(i_b) = i_abc
                           if(present(list_C))list_C(i_c) = i_abc
                           i_a   = i_a  + 1
                           i_b   = i_b  + 1
                           i_c   = i_c  + 1
                           i_abc = i_abc + 1
                        else
                           write(*,*),'Error in subroutine SPMP_ABC',Aindex(i_a),Bindex(i_b),Cindex(i_c)
                        end if
                     end if
                  end if
               end if
            end if
         end if
      end if
   end do

   this%NE = i_abc - 1
!   print*,'1st:this%NE',this%NE
   this%CI(1:this%NE) = (ABCindex(1:this%NE) - 1) / this%SM + 1
!   print*,'1st:this%ci',this%CI(1:this%NE)
   this%RI(1:this%NE) = ABCindex(1:this%NE) - this%SM * (this%CI(1:this%NE) - 1)
!   print*,'1st:this%ri',this%RI(1:this%NE)
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
!   print*,'A%NE',this%NE,B%NE,list_A(:),list_B(:)
   if(present(num_actunit_A))then
      do i_a = 1, A%NE
      !print*,'i_a',i_a
!      if(A%M(i_a) .ne. 0)then
         do fu = 1, num_actunit_A
            u = filter_actunit_A(fu)
            this%M(u,list_A(i_a)) = A%M(u,i_a)
         end do
!      end if
      end do
   else
      do i_a = 1, A%NE
      !print*,'i_a',i_a
!      if(A%M(i_a) .ne. 0)then
         do fu = 1, num_unit
            u = filter_u(fu)
            this%M(u,list_A(i_a)) = A%M(u,i_a)
         end do
!      end if
      end do
   end if
   if(present(num_actunit_B))then
      do i_b = 1, B%NE
!      if(B%M(i_b) .ne. 0)then
         do fu = 1, num_actunit_B
            u = filter_actunit_B(fu)
            this%M(u,list_B(i_b)) = this%M(u,list_B(i_b)) + B%M(u,i_b)
         end do
!      end if
      end do
   else
      do i_b = 1, B%NE
!      if(B%M(i_b) .ne. 0)then
         do fu = 1, num_unit
            u = filter_u(fu)
            this%M(u,list_B(i_b)) = this%M(u,list_B(i_b)) + B%M(u,i_b)
         end do
!      end if
      end do
   end if
   if(present(num_actunit_C))then
      do i_c = 1, C%NE
!      if(B%M(i_b) .ne. 0)then
         do fu = 1, num_actunit_C
            u = filter_actunit_C(fu)
            this%M(u,list_C(i_c)) = this%M(u,list_C(i_c)) + C%M(u,i_c)
         end do
!        end if
      end do
   else
      do i_c = 1, C%NE
!      if(B%M(i_b) .ne. 0)then
         do fu = 1, num_unit
            u = filter_u(fu)
            this%M(u,list_C(i_c)) = this%M(u,list_C(i_c)) + C%M(u,i_c)
         end do
!        end if
      end do
   end if
   this%NE = NE_ABC
!   print*,'2nd:this%ne',NE_ABC
   this%CI(1:this%NE) = CI_ABC(1:NE_ABC)
!   print*,'2nd:this%ci',this%CI(1:this%NE)
   this%RI(1:this%NE) = RI_ABC(1:NE_ABC)
!   print*,'2nd:this%ri',this%RI(1:this%NE)
end if

end subroutine SPMP_ABC

end module SPMMod
