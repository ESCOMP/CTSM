#include <misc.h>
#include <preproc.h>

!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_restwrite
!
! !INTERFACE:
logical function do_restwrite()
!
! !DESCRIPTION:
! Determine if restart dataset is to be written at this time step
!
! !USES:
#if (defined COUP_CSM)
  use clm_csmMod  , only : csmstop_next, csmrstrt
#else
  use time_manager, only : is_last_step
  use histFileMod , only : if_writrest
#endif
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

  do_restwrite = .false.

#if (defined OFFLINE) || (defined COUP_CAM)

  ! Write restart if end of run or if time to dispose master history file

  if (is_last_step() .or. if_writrest) do_restwrite = .true.

#elif (defined COUP_CSM)

  ! Write restart only if coupler says to

  if (csmrstrt) do_restwrite = .true.

#endif

end function do_restwrite
