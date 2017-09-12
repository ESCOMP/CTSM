module computeFluxMod

  ! provides interface for the flux modules

  use shr_kind_mod , only : r8 => shr_kind_r8
  use abortutils   , only : endrun
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use clm_varctl   , only : iulog

  implicit none
  private

  ! Different allowed values for the solution method
  integer, parameter :: SOLUTION_METHOD_UNSET          = 0
  integer, parameter :: SOLUTION_METHOD_EXPLICIT_EULER = 1
  integer, parameter :: SOLUTION_METHOD_IMPLICIT_EULER = 2
  integer, parameter :: SOLUTION_METHOD_ANALYTICAL     = 3

  integer, parameter :: flux_name_maxlen = 128

  public :: computeFlux_type
  type, abstract :: computeFlux_type
     private

     character(len=flux_name_maxlen), public :: flux_name = 'UNSET'  ! name of this flux, for output messages
     integer :: solution_method = SOLUTION_METHOD_UNSET
   contains
     procedure :: setMetadata  ! Set metadata that are common to all type extensions
     procedure :: solveForFlux ! Solve for this flux using the chosen solution method

     ! At least one of the following methods should be overridden by each derived class.
     ! Methods that remain unimplemented will restrict the valid values of
     ! solution_method for the given flux.
     procedure :: computeFluxExplicitEuler
     procedure :: computeFluxImplicitEuler
     procedure :: computeFluxAnalytical

     ! This method should also be overridden if a derived class wants to be able to use a
     ! solution method that requires this method, such as implicitEuler. If it is left
     ! unimplemented, then those solution methods will not be available for the given flux.
     procedure :: getFlux  ! Get a flux estimate for use in an iterative solution procedure
  end type computeFlux_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine setMetadata(this, flux_name, solution_method_str)
    !
    ! !DESCRIPTION:
    ! Set metadata that are common to all type extensions of this base class.
    !
    ! solution_method_str can be one of:
    ! - 'explicit_euler'
    ! - 'implicit_euler'
    ! - 'analytical'
    !
    ! !ARGUMENTS:
    class(computeFlux_type), intent(inout) :: this
    character(len=*), intent(in) :: flux_name            ! name of this flux, for output messages
    character(len=*), intent(in) :: solution_method_str  ! name of this solution method
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'setMetadata'
    !-----------------------------------------------------------------------

    if (len_trim(flux_name) > flux_name_maxlen) then
       write(iulog,*) 'flux_name too long'
       write(iulog,*) trim(flux_name), ' exceeds max length: ', flux_name_maxlen
       call endrun(msg='flux_name too long: ' // errMsg(sourcefile, __LINE__))
    end if
    this%flux_name = flux_name

    select case (solution_method_str)
    case ('explicit_euler')
       this%solution_method = SOLUTION_METHOD_EXPLICIT_EULER
    case ('implicit_euler')
       this%solution_method = SOLUTION_METHOD_IMPLICIT_EULER
    case ('analytical')
       this%solution_method = SOLUTION_METHOD_ANALYTICAL
    case default
       write(iulog,*) 'Unknown solution_method: ', trim(solution_method_str)
       call endrun(msg='Unknown solution method: ' // errMsg(sourcefile, __LINE__))
    end select

  end subroutine setMetadata

  !-----------------------------------------------------------------------
  subroutine solveForFlux(this, x, f)
    !
    ! !DESCRIPTION:
    ! Solve for this flux using the chosen solution method
    !
    ! !ARGUMENTS:
    class(computeFlux_type), intent(in) :: this
    real(r8), intent(in)  :: x ! state
    real(r8), intent(out) :: f ! flux
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'solveForFlux'
    !-----------------------------------------------------------------------

    select case(this%solution_method)
    case(SOLUTION_METHOD_EXPLICIT_EULER)
       call this%computeFluxExplicitEuler(x, f)
    case(SOLUTION_METHOD_IMPLICIT_EULER)
       ! NOTE(wjs, 2017-09-11) At first I thought we could put the call to implicitEuler
       ! right here. (Aside: in order to avoid a circular dependency, that would require
       ! moving the implicitEuler routine into this module, or having this solveForFlux in
       ! a 3rd module that uses both computeFluxMod and implicitEulerMod.) But then I
       ! realized that not all fluxes (e.g., QflxH2osfcSurf) use that iterative
       ! implicitEuler method. So we need to allow each flux to determine how to compute
       ! its own implicit euler-based solution - including whether or not this actually
       ! involves a call to implicitEuler.
       call this%computeFluxImplicitEuler(x, f)
    case(SOLUTION_METHOD_ANALYTICAL)
       call this%computeFluxAnalytical(x, f)
    case default
       write(iulog,*) 'Unknown solution method: ', this%solution_method
       call endrun(msg='Unknown solution method: ' // errMsg(sourcefile, __LINE__))
    end select

  end subroutine solveForFlux

  ! ========================================================================
  ! At least one of the following methods should be overridden by each derived class. We
  ! are providing empty implementations here (which abort the run if called) to allow
  ! derived classes to only provide implementations of some subset of these methods. But
  ! note that methods that remain unimplemented will restrict the valid values of
  ! solution_method for the given flux.
  ! ========================================================================

  subroutine computeFluxExplicitEuler(this, x, f)
    class(computeFlux_type), intent(in) :: this
    real(r8), intent(in)  :: x ! state
    real(r8), intent(out) :: f ! flux, f(x)
    call endrun("computeFluxExplicitEuler not implemented for "//trim(this%flux_name))
  end subroutine computeFluxExplicitEuler

  subroutine computeFluxImplicitEuler(this, x, f)
    class(computeFlux_type), intent(in) :: this
    real(r8), intent(in)  :: x ! state
    real(r8), intent(out) :: f ! flux, f(x)
    call endrun("computeFluxImplicitEuler not implemented for "//trim(this%flux_name))
  end subroutine computeFluxImplicitEuler

  subroutine computeFluxAnalytical(this, x, f)
    class(computeFlux_type), intent(in) :: this
    real(r8), intent(in)  :: x ! state
    real(r8), intent(out) :: f ! flux, f(x)
    call endrun("computeFluxAnalytical not implemented for "//trim(this%flux_name))
  end subroutine computeFluxAnalytical

  ! In contrast to the other unimplemented subroutines in the base class, this one is NOT
  ! called directly by solveForFlux. Rather, it is called by some solution methods, such
  ! as implicitEuler. If it is left unimplemented, then those solution methods will not
  ! be available for the given flux.
  subroutine getFlux(this, x, f, dfdx)
    class(computeFlux_type), intent(in) :: this
    real(r8), intent(in)            :: x    ! state
    real(r8), intent(out)           :: f    ! f(x)
    real(r8), intent(out), optional :: dfdx ! derivative
    call endrun("getFlux not implemented for "//trim(this%flux_name))
  end subroutine getFlux

end module computeFluxMod
