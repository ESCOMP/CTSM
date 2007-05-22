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
subroutine do_disp (ntapes, hist_ntimes, hist_mfilt, if_stop, if_disphist, if_remvhist, &
	            rstwr, nlend)
!
! !DESCRIPTION:
! Determine logic for closeing and/or disposing history file
! Sets values for if_disphist, if_stop, if_remvhist (arguments)
! Remove history files unless this is end of run or
! history file is not full.
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort
#if (defined COUP_CSM)
  use clm_csmMod      , only : csmstop_next, csmrstrt
#else
  use clm_time_manager, only : is_last_step
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
  logical, intent(in), optional :: rstwr
  logical, intent(in), optional :: nlend	
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: t                   ! history tape index
  logical :: rest_now            ! temporary
  logical :: stop_now            ! temporary
!------------------------------------------------------------------------

  rest_now = .false.
  stop_now = .false.

#if (defined SEQ_MCT) || (defined SEQ_ESMF)
  if (present(nlend) .and. present(rstwr)) then  
     if (nlend) stop_now = .true.
     if (rstwr) rest_now = .true.
  else
     call shr_sys_abort('do_close_dispose error: must specify nlend and rstwr')
  end if
#elif (defined OFFLINE)
  if (is_last_step()                 ) stop_now = .true.
  if (hist_ntimes(1) == hist_mfilt(1)) rest_now = .true.
#elif (defined COUP_CSM)
  if (csmstop_next) stop_now = .true.
  if (csmrstrt    ) rest_now = .true.
#endif

  if_stop = stop_now

  if (stop_now) then
     ! End of run -  dispose all history files and do not remove any

     if_disphist(1:ntapes) = .true.
     if_remvhist(1:ntapes) = .false.

  else if (rest_now) then
     ! Restart - dispose all history files, remove only the ones that are full

     do t = 1,ntapes
        if_disphist(t) = .true.
        if_remvhist(t) = .false.
        if (hist_ntimes(t)==hist_mfilt(t)) if_remvhist(t) = .true.
     end do
  else
     ! Dispose and remove only full files

     if_disphist(1:ntapes) = .false.
     if_remvhist(1:ntapes) = .false.
     do t = 1,ntapes
        if (hist_ntimes(t) ==  hist_mfilt(t)) then
           if_disphist(t) = .true.
           if_remvhist(t) = .true.
        endif
     end do
  endif
  
end subroutine do_disp

end module do_close_dispose
