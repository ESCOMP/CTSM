module mkpctPftTypeMod

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: mkpctPftType
  !
  ! !DESCRIPTION:
  ! Derived type and associated methods for operating on pct_pft data
  !
  ! !REVISION HISTORY:
  ! Author: Bill Sacks
  !
  !-----------------------------------------------------------------------

  !!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  ! !PUBLIC TYPES
  public :: pct_pft_type
  
  type :: pct_pft_type
     private
     real(r8), allocatable :: pct_p2l(:) ! pct of each pft on the landunit
     real(r8)              :: pct_l2g    ! pct of landunit on the grid cell
   contains
     ! Public routines:
     ! Query routines:
     procedure :: get_pct_p2l         ! get an array holding % of each pft on the landunit
     procedure :: get_pct_p2g         ! get an array holding % of each pft on the gridcell
     procedure :: get_pct_l2g         ! get % of landunit on the grid cell
     procedure :: get_first_pft_index ! get index of the first pft (lower bound of arrays)
     procedure :: get_one_pct_p2g     ! get % of gridcell for a single pft
     ! Routines that modify the data:
     procedure :: set_pct_l2g         ! set % of landunit on the grid cell
     procedure :: set_one_pct_p2g     ! set % pft for a single pft
     procedure :: merge_pfts          ! merge all area from one PFT into another PFT
     procedure :: remove_small_cover  ! set % cover to 0 for any PFT whose grid cell coverage is less than a threshold

     ! Private routines:
     procedure, private :: convert_from_p2g   ! convert a p2g array into p2l and l2g
     procedure, private :: check_vals         ! perform a sanity check after setting values
  end type pct_pft_type

  ! !PUBLIC MEMBER FUNCTIONS
  public :: update_max_array ! given an array of pct_pft_type variables update the max_p2l values from pct_p2l
  public :: get_pct_p2l_array ! given an array of pct_pft_type variables, return a 2-d array of pct_p2l
  public :: get_pct_l2g_array ! given an array of pct_pft_type variables, return an array of pct_l2g

  interface pct_pft_type
     module procedure constructor       ! initialize a new pct_pft_type object
     module procedure constructor_pong  ! initialize a new pct_pft_type object with all PFT's on the gridcell
     module procedure constructor_empty ! initialize a new pct_pft_type object for an empty landunit
  end interface pct_pft_type

  ! !PRIVATE TYPES:
  real(r8), parameter :: tol = 1.e-12_r8  ! tolerance for checking equality

  !EOP

contains
  
  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor_pong(pct_p2g, first_pft_index, default_pct_p2l) result(this)
    !
    ! !DESCRIPTION:
    ! Given the % of each pft on the grid cell, create a pct_pft_type object.
    !
    ! Note that pct_p2g should just contain the pfts in this landunit.
    !
    ! If all PFTs have 0 weight on the grid cell, we arbitrarily set % of each pft on the
    ! landunit based on default_pct_p2l. Note that:
    ! (1) size of default_pct_p2l must match size of pct_p2g
    ! (2) default_pct_p2l must sum to 100%
    !
    ! !ARGUMENTS:
    type(pct_pft_type) :: this  ! function result
    
    real(r8), intent(in) :: pct_p2g(:)         ! % of each pft on the grid cell
    integer , intent(in) :: first_pft_index    ! index of the first pft (lower bound of arrays)
    real(r8), intent(in) :: default_pct_p2l(:) ! default % of each pft on the landunit, used if total landunit area is 0%
    !
    ! !LOCAL VARIABLES:
    integer :: last_pft_index
    
    character(len=*), parameter :: subname = 'constructor_pong'
    !-----------------------------------------------------------------------

    if (size(default_pct_p2l) /= size(pct_p2g)) then
       write(6,*) subname//' ERROR: size of default_pct_p2l must match size of pct_p2g'
       call abort()
    end if

    last_pft_index = first_pft_index + size(pct_p2g) - 1
    allocate(this%pct_p2l(first_pft_index : last_pft_index))
    call this%convert_from_p2g(pct_p2g, default_pct_p2l)

  end function constructor_pong

  !-----------------------------------------------------------------------
  function constructor(pct_p2l, pct_l2g, first_pft_index) result(this)
    !
    ! !DESCRIPTION:
    ! Given the % of each pft on the land cell and % of land unit on grid cell, 
    ! create a pct_pft_type object.
    !
    ! Note that pct_p2g should just contain the pfts in this landunit.
    !
    ! !ARGUMENTS:
    type(pct_pft_type) :: this  ! function result
    
    real(r8), intent(in) :: pct_p2l(:)         ! % of each pft on the landunit
    real(r8), intent(in) :: pct_l2g            ! % of the landunit on the grid cell
    integer , intent(in) :: first_pft_index    ! index of the first pft (lower bound of arrays)
    !
    ! !LOCAL VARIABLES:
    integer :: last_pft_index
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    last_pft_index = first_pft_index + size(pct_p2l) - 1
    allocate(this%pct_p2l(first_pft_index : last_pft_index))
    this%pct_p2l = pct_p2l
    this%pct_l2g = pct_l2g

  end function constructor

  !-----------------------------------------------------------------------
  function constructor_empty() result(this)
    !
    ! !DESCRIPTION:
    ! Initialize a new pct_pft_type object for an empty landunit - that is, one that has
    ! no PFTs on it, and never can (e.g., the crop landunit when we're running without
    ! prognostic crops, so that the landunit is always empty).
    !
    ! !ARGUMENTS:
    type(pct_pft_type) :: this  ! function result
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'constructor_empty'
    !-----------------------------------------------------------------------
    
    this%pct_l2g = 0._r8
    allocate(this%pct_p2l(0))

  end function constructor_empty

  
  
  ! ========================================================================
  ! Public member functions
  ! ========================================================================

  !-----------------------------------------------------------------------
  function get_pct_p2l(this) result(pct_p2l)
    !
    ! !DESCRIPTION:
    ! Get an array holding % of each pft on the landunit
    !
    ! !ARGUMENTS:
    class(pct_pft_type), intent(in) :: this
    real(r8) :: pct_p2l(size(this%pct_p2l))  ! function result
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_pct_p2l'
    !-----------------------------------------------------------------------
    
    pct_p2l = this%pct_p2l

  end function get_pct_p2l

  !-----------------------------------------------------------------------
  function get_pct_p2g(this) result(pct_p2g)
    !
    ! !DESCRIPTION:
    ! Get an array holding % of each pft on the gridcell
    !
    ! !ARGUMENTS:
    class(pct_pft_type), intent(in) :: this
    real(r8) :: pct_p2g(size(this%pct_p2l))  ! function result
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_pct_p2g'
    !-----------------------------------------------------------------------
    
    pct_p2g(:) = this%pct_p2l(:) * this%pct_l2g / 100._r8

  end function get_pct_p2g

  !-----------------------------------------------------------------------
  function get_pct_l2g(this) result(pct_l2g)
    !
    ! !DESCRIPTION:
    ! Get % of landunit on the grid cell
    !
    ! !ARGUMENTS:
    real(r8) :: pct_l2g  ! function result
    class(pct_pft_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_pct_l2g'
    !-----------------------------------------------------------------------
    
    pct_l2g = this%pct_l2g

  end function get_pct_l2g

  !-----------------------------------------------------------------------
  function get_first_pft_index(this) result(first_pft_index)
    !
    ! !DESCRIPTION:
    ! Get index of the first pft (lower bound of arrays)
    !
    ! !ARGUMENTS:
    integer :: first_pft_index  ! function result
    class(pct_pft_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_first_pft_index'
    !-----------------------------------------------------------------------
    
    first_pft_index = lbound(this%pct_p2l, 1)

  end function get_first_pft_index

  !-----------------------------------------------------------------------
  function get_one_pct_p2g(this, pft_index) result(pct_p2g)
    !
    ! !DESCRIPTION:
    ! Get % of gridcell for a single pft
    !
    ! !ARGUMENTS:
    real(r8) :: pct_p2g  ! function result
    class(pct_pft_type), intent(in) :: this
    integer :: pft_index
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'get_one_pct_p2g'
    !-----------------------------------------------------------------------
    
    pct_p2g = this%pct_p2l(pft_index) * this%pct_l2g / 100._r8

  end function get_one_pct_p2g

  !-----------------------------------------------------------------------
  subroutine set_pct_l2g(this, pct_l2g_new)
    !
    ! !DESCRIPTION:
    ! Set percent of landunit on the grid cell. Keep pct_p2l the same as before.
    !
    ! !ARGUMENTS:
    class(pct_pft_type), intent(inout) :: this
    real(r8), intent(in) :: pct_l2g_new  ! new percent of this landunit with respect to grid cell
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'set_pct_l2g'
    !-----------------------------------------------------------------------
    
    if (pct_l2g_new < 0._r8 .or. pct_l2g_new > (100._r8 + tol)) then
       write(6,*) subname//' ERROR: pct_l2g_new must be between 0 and 100%'
       write(6,*) 'pct_l2g_new = ', pct_l2g_new
       call abort()
    end if

    this%pct_l2g = pct_l2g_new

  end subroutine set_pct_l2g

  !-----------------------------------------------------------------------
  subroutine set_one_pct_p2g(this, pft_index, pct_p2g_new)
    !
    ! !DESCRIPTION:
    ! Set percent pft for a single pft, given its weight on the grid cell.
    !
    ! The landunit percent is adjusted appropriately. In addition, the coverage of other
    ! PFTs are adjusted proportionally so that the total pct_pft adds to 100%.
    !
    ! If the resulting total weight on the grid cell is reduced to 0, then pct_p2l
    ! remains as it was before this subroutine call.
    !
    ! Note about pft_index: Note that the first element of the array has index given by
    ! the first_pft_index value given to the constructor.
    ! 
    ! !ARGUMENTS:
    class(pct_pft_type), intent(inout) :: this
    integer , intent(in) :: pft_index  ! index of the pft to change
    real(r8), intent(in) :: pct_p2g_new    ! new percent of this pft, with respect to grid cell 
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: pct_p2g(:) ! % of each pft on the grid cell

    character(len=*), parameter :: subname = 'set_pct_p2g'
    !-----------------------------------------------------------------------
    
    if (pct_p2g_new < 0._r8 .or. pct_p2g_new > (100._r8 + tol)) then
       write(6,*) subname//' ERROR: pct_p2g_new must be between 0 and 100%'
       write(6,*) 'pct_p2g_new = ', pct_p2g_new
       call abort()
    end if
    
    allocate(pct_p2g(lbound(this%pct_p2l, 1) : ubound(this%pct_p2l, 1)))
    pct_p2g(:) = this%get_pct_p2g()
    pct_p2g(pft_index) = pct_p2g_new

    ! Note that by using this%pct_p2l as the default_pct_pl2 argument, we ensure that, if
    ! the new p2g value brings the total % on the grid cell to 0, then we keep the
    ! previous values for pct_p2l
    call this%convert_from_p2g(pct_p2g, this%pct_p2l)

    deallocate(pct_p2g)

  end subroutine set_one_pct_p2g

  !-----------------------------------------------------------------------
  subroutine merge_pfts(this, source, dest)
    !
    ! !DESCRIPTION:
    ! Merge all area from one PFT into another PFT
    !
    ! !ARGUMENTS:
    class(pct_pft_type), intent(inout) :: this
    integer, intent(in) :: source  ! index of source PFT
    integer, intent(in) :: dest    ! index of dest PFT
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'merge_pfts'
    !-----------------------------------------------------------------------
    
    this%pct_p2l(dest) = this%pct_p2l(dest) + this%pct_p2l(source)
    this%pct_p2l(source) = 0._r8

    call this%check_vals(subname)

  end subroutine merge_pfts

  !-----------------------------------------------------------------------
  subroutine remove_small_cover(this, too_small, nsmall)
    !
    ! !DESCRIPTION:
    ! Remove any small PFTs, defined as those whose grid cell coverage is below some
    ! threshold. Also returns the number of small PFTs found.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(pct_pft_type), intent(inout) :: this
    real(r8), intent(in)  :: too_small ! threshold for considering a PFT too small (% of grid cell)
    integer , intent(out) :: nsmall    ! number of small (but non-zero) PFTs found
    !
    ! !LOCAL VARIABLES:
    integer :: pft_lbound
    integer :: pft_ubound
    integer :: pft_index
    real(r8), allocatable :: pct_p2g(:)  ! % of each pft on the grid cell
    logical , allocatable :: is_small(:) ! whether each PFT is considered too small (but not 0)
    logical , allocatable :: is_zero(:)  ! whether each PFT is exactly 0

    character(len=*), parameter :: subname = 'remove_small_cover'
    !-----------------------------------------------------------------------
    
    pft_lbound = lbound(this%pct_p2l, 1)
    pft_ubound = ubound(this%pct_p2l, 1)
    allocate(pct_p2g (pft_lbound : pft_ubound))
    allocate(is_small(pft_lbound : pft_ubound))
    allocate(is_zero (pft_lbound : pft_ubound))

    pct_p2g(:) = this%get_pct_p2g()
    is_zero(:) = (pct_p2g == 0._r8)
    is_small(:) = (pct_p2g < too_small .and. .not. is_zero(:))
    
    nsmall = count(is_small(:))

    if (nsmall > 0) then

       if (all(is_zero(:) .or. is_small(:))) then
          ! If all PFTs are either 0 or small, then set pct_l2g to 0, but don't touch
          ! pct_p2l(:) (We do NOT set pct_p2l to all 0 in this case, because we need to
          ! maintain sum(pct_p2l) = 100%)
          this%pct_l2g = 0._r8

       else
          ! If there are some big PFTs, then we need to adjust pct_p2l as well as pct_l2g
          ! (setting pct_p2l to 0 for the small elements and renormalizing the others)
          do pft_index = pft_lbound, pft_ubound
             if (is_small(pft_index)) then
                call this%set_one_pct_p2g(pft_index, 0._r8)
             end if
          end do
       end if

       call this%check_vals(subname)
    end if

    deallocate(pct_p2g, is_small, is_zero)
  end subroutine remove_small_cover

  ! ========================================================================
  ! Private member functions
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine convert_from_p2g(this, pct_p2g, default_pct_p2l)
    !
    ! !DESCRIPTION:
    ! Given a p2g array, compute the p2l array and l2g
    !
    ! !ARGUMENTS:
    class(pct_pft_type), intent(inout) :: this
    real(r8), intent(in) :: pct_p2g(:)         ! % of each pft on the grid cell
    real(r8), intent(in) :: default_pct_p2l(:) ! default % of each pft on the landunit, used if total landunit area is 0%
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'convert_from_p2g'
    !-----------------------------------------------------------------------
    
    ! Check pre-conditions

    if (size(pct_p2g) /= size(this%pct_p2l) .or. size(default_pct_p2l) /= size(this%pct_p2l)) then
       write(6,*) subname//' ERROR: array size mismatch: '
       write(6,*) size(pct_p2g), size(default_pct_p2l), size(this%pct_p2l)
       call abort()
    end if

    if (abs(sum(default_pct_p2l) - 100._r8) > tol) then
       write(6,*) subname//' ERROR: default_pct_p2l must sum to 100'
       call abort()
    end if

    if (any(pct_p2g < 0._r8)) then
       write(6,*) subname//' ERROR: negative values found in pct_p2g array'
       write(6,*) pct_p2g
       call abort()
    end if

    if (sum(pct_p2g) < 0._r8 .or. sum(pct_p2g) > (100._r8 + tol)) then
       write(6,*) subname//' ERROR: pct_p2g must be between 0 and 100'
       write(6,*) 'sum(pct_p2g) = ', sum(pct_p2g)
       call abort()
    end if

    ! Done checking pre-conditions

    this%pct_l2g = sum(pct_p2g)
    if (this%pct_l2g > 0._r8) then
       this%pct_p2l = pct_p2g / this%pct_l2g * 100._r8
    else
       this%pct_p2l = default_pct_p2l
    end if

    ! Check post-conditions

    call this%check_vals(subname)

  end subroutine convert_from_p2g


  !-----------------------------------------------------------------------
  subroutine check_vals(this, caller)
    !
    ! !DESCRIPTION:
    ! Perform a sanity check after setting values
    !
    ! !ARGUMENTS:
    class(pct_pft_type), intent(in) :: this
    character(len=*), intent(in) :: caller   ! name of the calling subroutine
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'check_vals'
    !-----------------------------------------------------------------------
    
    if (abs(sum(this%pct_p2l) - 100._r8) > tol) then
       write(6,*) subname//' ERROR from ', caller, ': pct_p2l does not sum to 100'
       write(6,*) 'sum(this%pct_p2l) = ', sum(this%pct_p2l)
       call abort()
    end if

    if (any(this%pct_p2l < 0._r8)) then
       write(6,*) subname//' ERROR from ', caller, ': negative values found in pct_p2l'
       write(6,*) this%pct_p2l
       call abort()
    end if

    if (this%pct_l2g < 0._r8 .or. this%pct_l2g > (100._r8 + tol)) then
       write(6,*) subname//' ERROR from ', caller, ': pct_l2g must be between 0 and 100'
       write(6,*) 'pct_l2g = ', this%pct_l2g
       call abort()
    end if

  end subroutine check_vals
    
  ! ========================================================================
  ! Module-level routines (not member functions)
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine update_max_array(pct_pft_max_arr,pct_pft_arr)
    !
    ! !DESCRIPTION:
    ! Given an array of pct_pft_type variables, update all the max_p2l variables.
    !
    ! Assumes that all elements of pct_pft_max_arr and pct_pft_arr have the same 
    ! size and lower bound for their pct_p2l array.
    !
    ! !ARGUMENTS:
    ! workaround for gfortran bug (58043): declare this 'type' rather than 'class':
    type(pct_pft_type), intent(inout) :: pct_pft_max_arr(:)
    type(pct_pft_type), intent(in) :: pct_pft_arr(:)
    !
    ! !LOCAL VARIABLES:
    integer :: pft_lbound
    integer :: pft_ubound
    integer :: arr_index
    integer :: pft_index
    
    character(len=*), parameter :: subname = 'update_max_array'
    !-----------------------------------------------------------------------
 
       
    pft_lbound = lbound(pct_pft_arr(1)%pct_p2l, 1)
    pft_ubound = ubound(pct_pft_arr(1)%pct_p2l, 1)

    do arr_index = 1, size(pct_pft_arr)
       if (lbound(pct_pft_arr(arr_index)%pct_p2l, 1) /= pft_lbound .or. &
            ubound(pct_pft_arr(arr_index)%pct_p2l, 1) /= pft_ubound) then
          write(6,*) subname//' ERROR: all elements of pct_pft_arr must have'
          write(6,*) 'the same size and lower bound for their pct_p2l array'
          call abort()
       end if
          
       if (pct_pft_arr(arr_index)%pct_l2g > pct_pft_max_arr(arr_index)%pct_l2g) then
          pct_pft_max_arr(arr_index)%pct_l2g = pct_pft_arr(arr_index)%pct_l2g
       end if

       do pft_index = pft_lbound, pft_ubound
          if (pct_pft_arr(arr_index)%pct_p2l(pft_index) > pct_pft_max_arr(arr_index)%pct_p2l(pft_index)) then
             pct_pft_max_arr(arr_index)%pct_p2l(pft_index) = pct_pft_arr(arr_index)%pct_p2l(pft_index)
	  end if
       end do
    end do

  end subroutine update_max_array

  !-----------------------------------------------------------------------
  function get_pct_p2l_array(pct_pft_arr) result(pct_p2l)
    !
    ! !DESCRIPTION:
    ! Given an array of pct_pft_type variables, return a 2-d array of pct_p2l.
    !
    ! Assumes that all elements of pct_pft_arr have the same size and lower bound for
    ! their pct_p2l array.
    !
    ! !ARGUMENTS:
    real(r8), allocatable :: pct_p2l(:,:)  ! function result (n_elements, n_pfts)
    ! workaround for gfortran bug (58043): declare this 'type' rather than 'class':
    type(pct_pft_type), intent(in) :: pct_pft_arr(:)
    !
    ! !LOCAL VARIABLES:
    integer :: pft_lbound
    integer :: pft_ubound
    integer :: arr_index
    integer :: pft_index
    
    character(len=*), parameter :: subname = 'get_pct_p2l_array'
    !-----------------------------------------------------------------------
    
    pft_lbound = lbound(pct_pft_arr(1)%pct_p2l, 1)
    pft_ubound = ubound(pct_pft_arr(1)%pct_p2l, 1)

    allocate(pct_p2l(size(pct_pft_arr), pft_lbound:pft_ubound))

    do arr_index = 1, size(pct_pft_arr)
       if (lbound(pct_pft_arr(arr_index)%pct_p2l, 1) /= pft_lbound .or. &
            ubound(pct_pft_arr(arr_index)%pct_p2l, 1) /= pft_ubound) then
          write(6,*) subname//' ERROR: all elements of pct_pft_arr must have'
          write(6,*) 'the same size and lower bound for their pct_p2l array'
          call abort()
       end if
          
       do pft_index = pft_lbound, pft_ubound
          pct_p2l(arr_index, pft_index) = pct_pft_arr(arr_index)%pct_p2l(pft_index)
       end do
    end do

  end function get_pct_p2l_array

  !-----------------------------------------------------------------------
  function get_pct_l2g_array(pct_pft_arr) result(pct_l2g)
    !
    ! !DESCRIPTION:
    ! Given an array of pct_pft_type variables, return an array of pct_l2g.
    !
    ! !ARGUMENTS:
    real(r8), allocatable :: pct_l2g(:)  ! function result
    class(pct_pft_type), intent(in) :: pct_pft_arr(:)
    !
    ! !LOCAL VARIABLES:
    integer :: arr_index
    
    character(len=*), parameter :: subname = 'get_pct_l2g_array'
    !-----------------------------------------------------------------------
    
    allocate(pct_l2g(size(pct_pft_arr)))
    pct_l2g = pct_pft_arr(:)%pct_l2g

  end function get_pct_l2g_array


end module mkpctPftTypeMod
