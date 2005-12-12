#include <misc.h>
#include <preproc.h>

!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_close_dispose
!
! !INTERFACE:
subroutine do_close_dispose (ntapes, hist_ntimes, hist_mfilt, &
                             if_stop, if_disphist, if_remvhist)
!
! !DESCRIPTION:
! Determine logic for closeing and/or disposing history file
! Sets values for if_disphist, if_stop, if_remvhist (arguments)
! Remove history files unless this is end of run or
! history file is not full.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
#if (defined COUP_CSM)
  use clm_csmMod  , only : csmstop_next, csmrstrt
#else
  use time_manager, only : is_last_step
#endif
!
! !ARGUMENTS:
  implicit none
  integer, intent(in)  :: ntapes              !actual number of history tapes
  integer, intent(in)  :: hist_ntimes(ntapes) !current numbers of time samples on history tape
  integer, intent(in)  :: hist_mfilt(ntapes)  !maximum number of time samples per tape
  logical, intent(out) :: if_stop             !true => last time step of run
  logical, intent(out) :: if_disphist(ntapes) !true => save and dispose history file
  logical, intent(out) :: if_remvhist(ntapes) !true => remove local history file after dispose
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer t                         !history tape index
!------------------------------------------------------------------------

#if (defined OFFLINE) || (defined COUP_CAM)

  if (is_last_step()) then

     ! End of run -  dispose all history files and write restart

     if_stop = .true.
     if_disphist(1:ntapes) = .true.
     if_remvhist(1:ntapes) = .false.

  else if (hist_ntimes(1) == hist_mfilt(1))  then

     ! Time to dispose master history file and write restart file
     ! dispose all other history files as well

     if_stop = .false.
     if_disphist(1:ntapes) = .true.
     if_remvhist(1) = .true.
     do t = 2,ntapes
        if_remvhist(t) = .false.
        if (hist_ntimes(t)==hist_mfilt(t)) if_remvhist(t) = .true.
     end do

  else

     ! If not end of run or time to dispose master history file
     ! then determine if time to dispose individual auxillary files

     if_stop = .false.
     if_disphist(1:ntapes) = .false.
     if_remvhist(1:ntapes) = .false.
     do t = 2,ntapes
        if (hist_ntimes(t) == hist_mfilt(t)) then
           if_disphist(t) = .true.
           if_remvhist(t) = .true.
        endif
     end do

  endif

#elif (defined COUP_CSM)

  if (csmstop_next) then

     ! If coupler says that next time step is end of run then
     ! dispose all history files and write restart

     if_stop = .true.
     if_disphist(1:ntapes) = .true.
     if_remvhist(1:ntapes) = .false.

  else if (csmrstrt) then

     ! If coupler says to write restart then dispose all history files

     if_stop = .false.
     do t = 1,ntapes
        if_disphist(t) = .true.
        if_remvhist(t) = .false.
        if (hist_ntimes(t)==hist_mfilt(t)) if_remvhist(t) = .true.
     end do

  else

     ! Otherwise check if file is full and dispose if it is

     if_stop = .false.
     if_disphist(1:ntapes) = .false.
     if_remvhist(1:ntapes) = .false.
     do t = 1,ntapes
        if (hist_ntimes(t) ==  hist_mfilt(t)) then
           if_disphist(t) = .true.
           if_remvhist(t) = .true.
        endif
     end do

  endif

#endif

  return

end subroutine do_close_dispose
