!-------------------------------------------------------------------------------
! PURPOSE: provide access to rtc clock rate on Cray X1/X1E systems.
!-------------------------------------------------------------------------------

function rtc_rate()
!----------------------------------------------------------------------- 
! 
! Purpose: Return RTC clock rate. (For use on Cray X1/X1E systems.)
!
! Method: Call IRTC_RATE.
! 
! Author: Pat Worley
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Output arguments
!
   integer(kind=8) :: rtc_rate
!
#ifdef UNICOSMP
!
! Externals
!
   integer(kind=8), external :: irtc_rate
!-----------------------------------------------------------------------

   rtc_rate = irtc_rate()
#else
!-----------------------------------------------------------------------

   rtc_rate = -1
#endif
   return
   end function rtc_rate
