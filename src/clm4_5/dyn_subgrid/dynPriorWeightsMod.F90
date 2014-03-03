module dynPriorWeightsMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type and associated methods for working with prior subgrid weights
  ! (i.e., before the weight updates of this time step)
  !
  ! !USES:
  use clmtype
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type, BOUNDS_LEVEL_PROC
  use shr_assert_mod , only : shr_assert
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun

  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: prior_weights_type

  type prior_weights_type
     ! Components are public for ease-of-use and efficiency. However, these components
     ! should be treated as read-only!
     real(r8), allocatable, public :: pwtcol(:)     ! prior pft weight on the column
   contains
     procedure :: set_prior_weights      ! set prior weights to current weights
  end type prior_weights_type

  interface prior_weights_type
     module procedure constructor   ! initialize a prior_weights_type object
  end interface prior_weights_type

contains
  
  ! ======================================================================
  ! Constructors
  ! ======================================================================

  ! ----------------------------------------------------------------------
  type(prior_weights_type) function constructor(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize a prior_weights_type object
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds   ! processor bounds
    !
    ! !LOCAL VARIABLES:
    integer :: ier  ! error code
    
    character(len=*), parameter :: subname = 'prior_weights_type constructor'
    ! ----------------------------------------------------------------------
     
    call shr_assert(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')
    
    allocate(constructor%pwtcol(bounds%begp:bounds%endp), stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for pwtcol'//errMsg(__FILE__, __LINE__))
    end if
  end function constructor


  ! ======================================================================
  ! Public methods
  ! ======================================================================
  
  ! ----------------------------------------------------------------------
  subroutine set_prior_weights(this, bounds)
    !
    ! !DESCRIPTION:
    ! Set prior weights to current weights
    !
    ! !ARGUMENTS:
    class(prior_weights_type) , intent(inout) :: this   ! this object
    type(bounds_type)         , intent(in)    :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p   ! pft index
    ! ----------------------------------------------------------------------
    
    do p = bounds%begp, bounds%endp
       this%pwtcol(p) = pft%wtcol(p)
    end do
  end subroutine set_prior_weights

end module dynPriorWeightsMod
