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
  implicit none
  save
  private

  ! Default settings valid even if there is no spmd 

  logical, public :: masterproc      ! proc 0 logical for printing msgs
  integer, public :: iam             ! processor number
  integer, public :: npes            ! number of processors for clm
  integer, public :: mpicom          ! communicator group for clm
  integer, public :: comp_id         ! component id

  !
  ! Public methods
  !
  public :: spmd_init                ! Initialization

  !
  ! Values from mpif.h that can be used
  !
  public :: MPI_INTEGER
  public :: MPI_REAL8
  public :: MPI_LOGICAL
  public :: MPI_SUM
  public :: MPI_MIN
  public :: MPI_MAX
  public :: MPI_STATUS_SIZE
  public :: MPI_ANY_SOURCE
  public :: MPI_CHARACTER
  public :: MPI_COMM_WORLD
  public :: MPI_MAX_PROCESSOR_NAME

#include <mpif.h>  

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_init( clm_mpicom )
!
! !INTERFACE:
  subroutine spmd_init( clm_mpicom )
!
! !DESCRIPTION:
! MPI initialization (number of cpus, processes, tids, etc)
!
! !USES
#if (defined COUP_CSM)
    use cpl_comm_mod, only : cpl_comm_mph_cid
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: clm_mpicom
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j         ! indices
    integer :: ier         ! return error status
    logical :: mpi_running ! temporary
    integer, allocatable :: length(:)
    integer, allocatable :: displ(:)
    character*(MPI_MAX_PROCESSOR_NAME), allocatable :: procname(:)
!-----------------------------------------------------------------------

    ! Initialize mpi communicator group

    mpicom = clm_mpicom

    comp_id = 1
#if (defined COUP_CSM)
    comp_id = cpl_comm_mph_cid
#endif

    ! Initialize mpi

#ifdef OFFLINE
    call mpi_initialized (mpi_running, ier)
    if (.not. mpi_running) call mpi_init (ier)
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

end module spmdMod
