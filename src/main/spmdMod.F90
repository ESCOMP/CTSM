#include <misc.h>
#include <preproc.h>

module spmdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spmdMod
!
! !DESCRIPTION:
! SPMD initialization
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
#if (defined COUP_CSM)
  use mpiinc
#endif

#if (!defined SPMD)
! ----------------------- SPMD OFF --------------------------------------

  logical :: masterproc = .true. ! proc 0 logical for printing msgs
  integer :: iam = 0
  integer :: npes = 1
  save
#endif

#if (defined SPMD)
! ----------------------- SPMD ON ---------------------------------------

#if (defined COUP_CAM)
  use mpishorthand
  use spmd_utils, only: npes, masterproc, iam
#elif (defined OFFLINE)
  use mpiinc
  integer, public :: npes        !number of processors
  integer, public :: iam         !proc number
  logical, public :: masterproc  !proc 0 logical for printing msgs
  integer, public :: mpicom
  save
#elif (defined COUP_CSM)
  integer, public :: npes        !number of processors
  integer, public :: iam         !proc number
  logical, public :: masterproc  !proc 0 logical for printing msgs
  integer, public :: mpicom
  save
#endif

#if (defined OFFLINE) || (defined COUP_CSM)

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_init
!
! !INTERFACE:
  subroutine spmd_init
!
! !DESCRIPTION:
! MPI initialization (number of cpus, processes, tids, etc)
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j        ! indices
    integer :: ier        ! return error status
    integer, allocatable :: length(:)
    integer, allocatable :: displ(:)
    character*(MPI_MAX_PROCESSOR_NAME), allocatable :: procname(:)
#if (defined OFFLINE)
    logical mpi_running
#endif
!-----------------------------------------------------------------------

#if (defined OFFLINE)

    ! Initialize mpi and set communication group

    call mpi_initialized (mpi_running, ier)
    if (.not. mpi_running) call mpi_init (ier)

    mpicom  = MPI_COMM_WORLD

#elif (defined COUP_CSM)

    ! Initialize mpi and set communication group
    ! Done in program_csm.F90

#endif

    ! Get my processor id

    call mpi_comm_rank(mpicom, iam, ier)
    if (iam==0) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    ! Get number of processors

    call mpi_comm_size(mpicom, npes, ier)

    ! Get my processor names

    allocate (length(0:npes-1), displ(0:npes-1), procname(0:npes-1))

    call mpi_get_processor_name (procname(iam), length(iam), ier)
    call mpi_allgather (length(iam),1,MPI_INTEGER,length,1,MPI_INTEGER,mpicom,ier)
    do i = 0,npes-1
       displ(i)=i*MPI_MAX_PROCESSOR_NAME
    end do
    call mpi_gatherv (procname(iam),length(iam),MPI_CHARACTER, &
                      procname,length,displ,MPI_CHARACTER,0,mpicom,ier)
    if (masterproc) then
       write(6,100)npes
       write(6,200)
       write(6,220)
       do i=0,npes-1
          write(6,250)i,(procname((i))(j:j),j=1,length(i))
       end do
    endif

    deallocate (length, displ, procname)

100 format(i3," pes participating in computation")
200 format(/,35('-'))
220 format(/,"NODE#",2x,"NAME")
250 format("(",i3,")",2x,100a1)

  end subroutine spmd_init

#endif

#endif

end module spmdMod
