module quadraticMod

  use abortutils  ,   only: endrun
  use shr_kind_mod,   only: r8 => shr_kind_r8, CX => shr_kind_cx
  use clm_varctl  ,   only: iulog

  implicit none

  private

  public :: quadratic_roots ! Solve for the two roots of a quadratic equation

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  subroutine quadratic_roots (a, b, c, root1, root2, file, line)
     !
     ! !DESCRIPTION:
     !==============================================================================!
     !----------------- Solve quadratic equation for its two roots -----------------!
     !==============================================================================!
     ! Implements the numerically stable formulation from Press et al (1986)
     ! Numerical Recipes: The Art of Scientific Computing (Cambridge University Press, Cambridge)
     ! Allows for roots that are technically complex if they are close to rounding to zero.
     !
     ! NOTE: Special handling for these cases...
     !   Will truncate the square root term to zero if it is very small and negative, otherwise will error out
     !   Root2 will be set to spval if it would be undefined by a division by zero
     !
     ! !REVISION HISTORY:
     ! 4/5/10:   Adapted from /home/bonan/ecm/psn/An_gs_iterative.f90 by Keith Oleson
     ! 5/14/18:  Modify endrun handling and handle truncation to zero of small negative square root term EBK
     ! 10/9/26:  Clarify variable names, add comments and pass in file and line number of calling routine EBK
     !
     ! !USES:
     use clm_varcon, only: spval

     implicit none
     !
     ! !ARGUMENTS:
     real(r8), intent(in)  :: a, b, c      ! Coefficients of the quadratic equation x = a*x^2 + b*x + c
     real(r8), intent(out) :: root1, root2 ! The two roots of the quadratic equation
     character(len=*), intent(in), optional :: file  ! File name of calling routine
     integer, intent(in), optional :: line           ! Line number of calling routine
     ! !LOCAL VARIABLES:
     real(r8) :: root1_x_a                 ! Temporary term for the solution first root multiplied by the "a" coefficient
     real(r8) :: discriminant              ! Temporary Term that will have a square root taken
     integer :: pline                      ! Line number to use in endrun calls
     character(len=CX) :: pfile            ! File name to use in endrun calls
     !------------------------------------------------------------------------------

     if ( present(file) )then
        pfile = file
     else
        pfile = sourcefile
     end if
     if ( present(line) )then
        pline = line
     else
        pline = __LINE__
     end if

     ! Initialize the roots to spval in case it exits early due to an error
     root1 = spval
     root2 = spval

     ! If the "a" coefficient is zero, then this is linear rather than quadratic which we assume is a mistake
     if (a == 0._r8) then
        write (iulog,*) 'ERROR: Quadratic solution error: the "a" coefficient is zero = ',a
        call endrun(msg='Quadratic solution error: input equation is linear not quadratic', file=pfile, line=pline)
        return
     end if

     ! Compute the term that will have a square root taken
     discriminant = b*b - 4._r8*a*c
     if ( discriminant < 0.0 )then
        ! Truncate to zero if the term is very small and negative (relative to b), otherwise error out
        if ( -discriminant < 3.0_r8*epsilon(b) )then
           discriminant = 0.0_r8
        else
           write (iulog,*) 'ERROR: Quadratic solution error: root would be complex b^2 - 4*a*c = ', discriminant
           call endrun( msg='Quadratic solution error: root would be complex', file=pfile, line=pline )
           return
        end if
     end if

     ! Compute the term that will be the first root multiplied by the "a" coefficient
     if (b >= 0._r8) then
        root1_x_a = -0.5_r8 * (b + sqrt(discriminant))
     else
        root1_x_a = -0.5_r8 * (b - sqrt(discriminant))
     end if

     ! Solve for the two roots to return
     root1 = root1_x_a / a
     if (root1_x_a /= 0._r8) then
        root2 = c / root1_x_a
     else
        root2 = spval
     end if

  end subroutine quadratic_roots

end module quadraticMod
