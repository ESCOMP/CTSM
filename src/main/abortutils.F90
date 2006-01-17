#include <misc.h>
#include <preproc.h>

module abortutils

!-----------------------------------------------------------------------
!BOP
! !MODULE: abortutils
!
! !DESCRIPTION:
! Abort the model for abnormal termination
!
! !REVISION HISTORY:
! Author: CCM Core group
!
!EOP
!-----------------------------------------------------------------------

#if (defined OFFLINE) || (defined COUP_CSM)

   private
   save

   public :: endrun

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: endrun
!
! !INTERFACE:
   subroutine endrun(msg)
!
! !DESCRIPTION:
! Abort the model for abnormal termination
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
! $Id: abortutils.F90,v 1.1.12.3 2005/09/09 19:48:29 mvertens Exp $
!-----------------------------------------------------------------------
#if (defined SPMD || defined COUP_CSM)
   use mpiinc
#endif
   use shr_sys_mod, only: shr_sys_flush
!
! !ARGUMENTS:
   implicit none
   character(len=*), intent(in), optional :: msg    ! string to be printed
!
! !REVISION HISTORY:
! Author: CCM Core group
!
!EOP
!-----------------------------------------------------------------------
! $Id: abortutils.F90,v 1.1.12.3 2005/09/09 19:48:29 mvertens Exp $
!-----------------------------------------------------------------------

   if (present (msg)) then
      write(6,*)'ENDRUN:', msg
   else
      write(6,*)'ENDRUN: called without a message string'
   end if

#if defined(NEC_SX)
   call mesput("ENDRUN", len("ENDRUN"), 1)
#elif defined(AIX)
   close(5)    ! needed to prevent batch jobs from hanging in xl__trbk
   call xl__trbk()
#endif

   call shr_sys_flush( 6 )   ! Flush all output to standard output

#if (defined SPMD) || (defined COUP_CSM)
   ! passing an argument of 1 to mpi_abort will lead to a STOPALL output
   ! error code of 257
   call mpi_abort (MPI_COMM_WORLD, 1)
#else
   call abort
#endif

end subroutine endrun

#endif

end module abortutils
