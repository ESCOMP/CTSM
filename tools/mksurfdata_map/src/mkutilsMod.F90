module mkutilsMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkutils
!
! !DESCRIPTION:
! General-purpose utilities for mksurfdata_map
!
!
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8
 
   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: slightly_below
   public :: slightly_above
!
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!EOP
!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: slightly_below
!
! !INTERFACE:
logical function slightly_below(a, b, eps)
!
! !DESCRIPTION:
! Returns true if a is slightly below b; false if a is significantly below b or if a is
! greater than or equal to b
!
! !USES:
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: a
   real(r8), intent(in) :: b

   ! if provided, eps gives the relative error allowed for checking the "slightly"
   ! condition; if not provided, the tolerance defaults to the value given by eps_default
   real(r8), intent(in), optional :: eps
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   real(r8) :: l_eps
   real(r8), parameter :: eps_default = 1.e-15_r8  ! default relative error tolerance
!------------------------------------------------------------------------------

   if (present(eps)) then
      l_eps = eps
   else
      l_eps = eps_default
   end if

   if (a < b .and. (b - a)/b < l_eps) then
      slightly_below = .true.
   else
      slightly_below = .false.
   end if

end function slightly_below
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: slightly_above
!
! !INTERFACE:
logical function slightly_above(a, b, eps)
!
! !DESCRIPTION:
! Returns true if a is slightly above b; false if a is significantly above b or if a is
! less than or equal to b
!
! !USES:
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: a
   real(r8), intent(in) :: b

   ! if provided, eps gives the relative error allowed for checking the "slightly"
   ! condition; if not provided, the tolerance defaults to the value given by eps_default
   real(r8), intent(in), optional :: eps
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   real(r8) :: l_eps
   real(r8), parameter :: eps_default = 1.e-15_r8  ! default relative error tolerance
!------------------------------------------------------------------------------

   if (present(eps)) then
      l_eps = eps
   else
      l_eps = eps_default
   end if

   if (a > b .and. (a - b)/b < l_eps) then
      slightly_above = .true.
   else
      slightly_above = .false.
   end if

end function slightly_above
!------------------------------------------------------------------------------



end module mkutilsMod
