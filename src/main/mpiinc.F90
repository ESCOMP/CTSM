#include <misc.h>
#include <preproc.h>

      module mpiinc
        public

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mpiinc
!
! !DESCRIPTION:
! Data and parameters used for MPI.  Note: The #include of "mpif.h" which
! is typically in f77 fixed format means that this module MUST be in fixed
! format. However, the advantage of this module is that modules which
! use this module (and theremore the mpif.h variables) must not be in f77
! fixed format.
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

#if (defined SPMD || defined COUP_CSM)
#include <mpif.h>
#endif

      end module mpiinc
