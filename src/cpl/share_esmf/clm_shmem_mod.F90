module clm_shmem_mod
  !-----------------------------------------------------------------------------
  ! Per-node MPI-3 shared-memory helper for large arrays that would otherwise be
  ! replicated identically on every MPI rank.  One physical copy is allocated per
  ! shared-memory node and mapped into every rank on that node, freeing
  ! (ranks_per_node - 1) copies per node.
  !
  ! Ported from CAM's cam_shmem_mod (src/utils/cam_shmem_mod.F90) and specialized
  ! for the CTSM decomposition setup: it provides a default-integer rank-1
  ! allocator (the CAM module only has real r4/r8 2d-5d wrappers) plus a
  ! node-leader sum-reduce that builds a globally-summed array in a node-shared
  ! buffer without every rank holding its own global-sized copy.
  !
  ! Usage (collective over the land communicator mpicom):
  !   call clm_shmem_alloc_i4_1d(ptr, win, n)        ! all ranks
  !   if (clm_shmem_is_leader()) ptr(:) = 0          ! leader owns the storage
  !   call clm_shmem_fence(win)                       ! publish the zeros
  !   <each rank stores into its disjoint indices of ptr>
  !   call clm_shmem_leader_allreduce_sum_i4(ptr,win,n) ! fence; sum across nodes; fence
  !   <all ranks may now read ptr for the rest of its lifetime>
  !   call clm_shmem_free(ptr, win)                   ! collective over the node comm
  !
  ! The MPI-3 shared-memory path is used for real MPI builds.  mpi-serial does not
  ! implement the MPI-2/MPI-3 one-sided / shared-memory interfaces, so for
  ! mpi-serial builds (CPP macro NO_MPI2, set by CIME for MPILIB=mpi-serial) a
  ! single-task fallback is compiled instead: each "node-shared" array is a plain
  ! local allocation (one task is its own node and its own leader, so there is no
  ! cross-rank sharing and the leader sum-reduce is a no-op).  The F90 'mpi' module
  ! (not mpif.h) is used because the TYPE(C_PTR) overloads of MPI_WIN_ALLOCATE_SHARED
  ! / MPI_WIN_SHARED_QUERY are only guaranteed there (MPI-3.0).
  !-----------------------------------------------------------------------------

  use mpi
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
  use spmdMod   , only : mpicom
  use abortutils, only : endrun

  implicit none
  private

  public :: clm_shmem_alloc_i4_1d          ! allocate a node-shared default-integer rank-1 array
  public :: clm_shmem_leader_allreduce_sum_i4 ! sum a node-shared array across nodes, in place
  public :: clm_shmem_free                 ! free a node-shared array (MPI_Win_free)
  public :: clm_shmem_fence                ! synchronize a window (publish writes)
  public :: clm_shmem_is_leader            ! .true. on the leader (rank 0) of this node
  public :: clm_shmem_leader_comm          ! communicator containing only node leaders
  public :: clm_shmem_npes_per_node        ! number of ranks sharing this node

  interface clm_shmem_free
     module procedure clm_shmem_free_i4_1d
  end interface clm_shmem_free

  ! Sentinel window handle used by the mpi-serial fallback (no real MPI window).
  integer, parameter :: SHMEM_WIN_NONE = -1

  logical, save :: initialized = .false.
  integer, save :: node_comm   = MPI_COMM_NULL  ! ranks sharing a node
  integer, save :: leader_comm = MPI_COMM_NULL  ! one rank per node (the leaders)
  integer, save :: node_rank   = 0
  integer, save :: node_size   = 1
  logical, save :: is_leader   = .true.

contains

  !=============================================================================
  subroutine init_comms()
    ! Lazily build the node-local and node-leader communicators.  Collective over
    ! mpicom; safe to call from every shared-memory request.
#ifndef NO_MPI2
    integer :: ierr, color
#endif

    if (initialized) return

#ifndef NO_MPI2
    call mpi_comm_split_type(mpicom, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &
                             node_comm, ierr)
    call mpi_comm_rank(node_comm, node_rank, ierr)
    call mpi_comm_size(node_comm, node_size, ierr)
    is_leader = (node_rank == 0)

    ! Communicator of node leaders only.  masterproc (global rank 0) is a leader.
    if (is_leader) then
       color = 0
    else
       color = MPI_UNDEFINED
    end if
    call mpi_comm_split(mpicom, color, 0, leader_comm, ierr)
#else
    ! mpi-serial: a single task is its own node and its own leader.
    node_rank = 0
    node_size = 1
    is_leader = .true.
#endif

    initialized = .true.
  end subroutine init_comms

  !=============================================================================
  subroutine clm_shmem_alloc_i4_1d(ptr, win, n)
    ! Allocate a node-shared default-integer array of length n.  Only the node
    ! leader requests storage; peers map the leader's contiguous segment.
    integer, pointer, intent(out) :: ptr(:)
    integer,          intent(out) :: win
    integer,          intent(in)  :: n

#ifndef NO_MPI2
    integer(kind=MPI_ADDRESS_KIND) :: winsize, qsize
    integer :: ierr, disp_unit, qdisp
    integer :: itmp
    type(c_ptr) :: baseptr
#else
    integer :: istat
#endif

    call init_comms()

#ifndef NO_MPI2
    disp_unit = storage_size(itmp) / 8   ! bytes per default integer (robust to -i8)
    if (is_leader) then
       winsize = int(n, MPI_ADDRESS_KIND) * int(disp_unit, MPI_ADDRESS_KIND)
    else
       winsize = 0_MPI_ADDRESS_KIND
    end if

    call mpi_win_allocate_shared(winsize, disp_unit, MPI_INFO_NULL, node_comm, &
                                 baseptr, win, ierr)
    if (ierr /= MPI_SUCCESS) call endrun('clm_shmem_mod: MPI_Win_allocate_shared failed')

    ! Non-leaders learn the address of the leader's (rank 0) contiguous segment.
    if (.not. is_leader) then
       call mpi_win_shared_query(win, 0, qsize, qdisp, baseptr, ierr)
       if (ierr /= MPI_SUCCESS) call endrun('clm_shmem_mod: MPI_Win_shared_query failed')
    end if

    call c_f_pointer(baseptr, ptr, [n])
#else
    ! mpi-serial: single task, no shared memory -- a plain local allocation.
    allocate(ptr(n), stat=istat)
    if (istat /= 0) call endrun('clm_shmem_mod: allocate failed (mpi-serial path)')
    win = SHMEM_WIN_NONE
#endif
  end subroutine clm_shmem_alloc_i4_1d

  !=============================================================================
  subroutine clm_shmem_leader_allreduce_sum_i4(ptr, win, n)
    ! Build a globally-summed array in the node-shared buffer ptr(1:n): fence so
    ! every rank's stores are visible, then the node leaders sum their per-node
    ! partials across nodes (over leader_comm) into the shared buffer, then fence
    ! to publish the result to all ranks on the node.  Collective over node_comm;
    ! every rank on the node must call it.
    integer, pointer,    intent(inout) :: ptr(:)
    integer,             intent(in)    :: win
    integer,             intent(in)    :: n

#ifndef NO_MPI2
    integer, allocatable :: tmp(:)
    integer :: ierr
#endif

    call clm_shmem_fence(win)             ! all node stores complete and visible to leader
#ifndef NO_MPI2
    if (is_leader) then
       allocate(tmp(n))
       call mpi_allreduce(ptr, tmp, n, MPI_INTEGER, MPI_SUM, leader_comm, ierr)
       if (ierr /= MPI_SUCCESS) call endrun('clm_shmem_mod: MPI_Allreduce failed')
       ptr(1:n) = tmp(1:n)
       deallocate(tmp)
    end if
#else
    ! mpi-serial: the single task owns the whole domain, so ptr already holds the
    ! global array -- there is nothing to sum across nodes.
#endif
    call clm_shmem_fence(win)             ! publish global result to all node ranks
  end subroutine clm_shmem_leader_allreduce_sum_i4

  !=============================================================================
  subroutine clm_shmem_fence(win)
    ! Collective over the node communicator; synchronizes the window so stores
    ! become visible to all ranks on the node.  A no-op for the mpi-serial path.
    integer, intent(in) :: win
#ifndef NO_MPI2
    integer :: ierr
    call mpi_win_fence(0, win, ierr)
    if (ierr /= MPI_SUCCESS) call endrun('clm_shmem_mod: MPI_Win_fence failed')
#endif
  end subroutine clm_shmem_fence

  !=============================================================================
  subroutine clm_shmem_free_i4_1d(ptr, win)
    ! Free the node-shared window and disassociate the pointer.  Collective over
    ! the node communicator; a no-op when win == MPI_WIN_NULL.
    integer, pointer       :: ptr(:)
    integer, intent(inout) :: win
#ifndef NO_MPI2
    integer :: ierr
    if (win /= MPI_WIN_NULL) then
       call mpi_win_free(win, ierr)
       if (ierr /= MPI_SUCCESS) call endrun('clm_shmem_mod: MPI_Win_free failed')
    end if
    if (associated(ptr)) nullify(ptr)
    win = MPI_WIN_NULL
#else
    ! mpi-serial: ptr was a plain local allocation, so deallocate it.
    if (associated(ptr)) deallocate(ptr)
    win = SHMEM_WIN_NONE
#endif
  end subroutine clm_shmem_free_i4_1d

  !=============================================================================
  logical function clm_shmem_is_leader()
    call init_comms()
    clm_shmem_is_leader = is_leader
  end function clm_shmem_is_leader

  !=============================================================================
  integer function clm_shmem_leader_comm()
    call init_comms()
    clm_shmem_leader_comm = leader_comm
  end function clm_shmem_leader_comm

  !=============================================================================
  integer function clm_shmem_npes_per_node()
    call init_comms()
    clm_shmem_npes_per_node = node_size
  end function clm_shmem_npes_per_node

end module clm_shmem_mod
