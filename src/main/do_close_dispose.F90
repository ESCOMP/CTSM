#include <misc.h>
#include <preproc.h>

module do_close_dispose

contains

!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_disp
!
! !INTERFACE:
subroutine do_disp (ntapes, hist_ntimes, hist_mfilt, if_stop, if_disphist, rstwr, nlend)
!
! !DESCRIPTION:
! Determine logic for closeing and/or disposing history file
! Sets values for if_disphist, if_stop (arguments)
! Remove history files unless this is end of run or
! history file is not full.
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort
  use clm_time_manager, only : is_last_step
!
! !ARGUMENTS:
  implicit none
  integer, intent(in)  :: ntapes              !actual number of history tapes
  integer, intent(in)  :: hist_ntimes(ntapes) !current numbers of time samples on history tape
  integer, intent(in)  :: hist_mfilt(ntapes)  !maximum number of time samples per tape
  logical, intent(out) :: if_stop             !true => last time step of run
  logical, intent(out) :: if_disphist(ntapes) !true => save and dispose history file
  logical, intent(in)  :: rstwr
  logical, intent(in)  :: nlend	
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: t                   ! history tape index
  logical :: rest_now            ! temporary
  logical :: stop_now            ! temporary
!------------------------------------------------------------------------

  rest_now = .false.
  stop_now = .false.

  if (nlend) stop_now = .true.
  if (rstwr) rest_now = .true.

  if_stop = stop_now

  if (stop_now) then
     ! End of run -  dispose all history files

     if_disphist(1:ntapes) = .true.

  else if (rest_now) then
     ! Restart - dispose all history files

     do t = 1,ntapes
        if_disphist(t) = .true.
     end do
  else
     ! Dispose

     if_disphist(1:ntapes) = .false.
     do t = 1,ntapes
        if (hist_ntimes(t) ==  hist_mfilt(t)) then
           if_disphist(t) = .true.
        endif
     end do
  endif
  
end subroutine do_disp

end module do_close_dispose
