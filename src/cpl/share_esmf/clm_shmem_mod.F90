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
  ! Terminology (MPI-3 shared-memory concepts used throughout this module):
  !   window - the shared allocation itself (an "MPI window", type MPI_Win), created once
  !            per node and mapped into the address space of every rank on that node; "win"
  !            is the integer handle used to reference and later free it.
  !   leader - the single rank per node that owns the physical storage (node-local rank 0,
  !            clm_shmem_is_leader()); it is the only rank that requests the allocation,
  !            initializes it, and reduces across nodes.
  !   fence  - an MPI_Win_fence synchronization: a collective call over the node that makes
  !            stores issued by one rank visible to the other ranks sharing the window.
  !   color  - the grouping key passed to MPI_Comm_split; here leaders get color 0 (so they
  !            are gathered into leader_comm) while non-leaders get MPI_UNDEFINED (so they
  !            are excluded from it).
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
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer, c_associated, c_null_ptr
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

  logical, save :: initialized   = .false.
  integer, save :: node_comm     = MPI_COMM_NULL ! ranks sharing a node
  integer, save :: leader_comm   = MPI_COMM_NULL ! one rank per node (the leaders)
  integer, save :: node_rank     = 0             ! this rank's index within node_comm (0 => leader)
  integer, save :: node_size     = 1             ! number of ranks sharing this node
  logical, save :: is_leader     = .true.        ! .true. on node rank 0
  ! Whether the MPI-3 shared-memory path is actually usable on this platform.  Decided
  ! once by a runtime probe in init_comms(): on some MPI stacks (e.g. certain MVAPICH2
  ! builds) MPI_Win_shared_query returns a null base to non-leaders, which would crash at
  ! first use.  When the probe fails, the module runs in a "private" fallback mode: every
  ! rank keeps its own full-size copy and the leader-reduce becomes a plain MPI_Allreduce
  ! over mpicom -- bit-for-bit with the pre-optimization all-rank reduce, just without the
  ! per-node memory saving.  Stays .false. for the mpi-serial (NO_MPI2) build.
  logical, save :: shared_active = .false.

contains

  !=============================================================================
  subroutine init_comms()
    ! Lazily build the node-local and node-leader communicators.  Collective over
    ! mpicom; safe to call from every shared-memory request.
    !
    ! NO_MPI2 is defined by CIME only for MPILIB=mpi-serial, whose stub library lacks the
    ! MPI-2/MPI-3 one-sided and shared-memory routines.  Throughout this module, therefore,
    ! "#ifndef NO_MPI2" selects the real-MPI shared-memory path and the "#else" branch
    ! selects the single-task mpi-serial fallback.
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

    ! Only trust the MPI-3 shared-memory path if a runtime probe confirms it works
    ! end-to-end on this node (see probe_shared_mem and the shared_active comment above).
    shared_active = probe_shared_mem()
#else
    ! mpi-serial: a single task is its own node and its own leader; no shared memory.
    node_rank = 0
    node_size = 1
    is_leader = .true.
    shared_active = .false.
#endif

    initialized = .true.
  end subroutine init_comms

#ifndef NO_MPI2
  !=============================================================================
  logical function probe_shared_mem() result(ok)
    ! Verify that MPI-3 shared memory actually works on this node before relying on it.
    ! The leader allocates a one-element shared window and writes a sentinel; every rank
    ! maps the leader's segment (leader: its own base; peers: MPI_Win_shared_query) and
    ! confirms it can read the sentinel back.  This catches both a null base pointer and a
    ! wrong/garbage mapping.  Collective over mpicom.  The final agreement is AND-reduced
    ! over mpicom (not just node_comm) so the WHOLE job picks one mode: a mix of shared and
    ! private ranks would deadlock the later reduce (shared leaders wait on leader_comm
    ! while private ranks wait on mpicom).  Returns the same value on every rank.
    integer, parameter :: SENTINEL = 1234567
    integer(kind=MPI_ADDRESS_KIND) :: winsize, qsize
    integer          :: ierr, disp_unit, qdisp, itmp, pwin
    integer          :: iok(1), iok_all(1)
    type(c_ptr)      :: baseptr
    integer, pointer :: p(:)
    logical          :: have_ptr, have_win

    disp_unit = storage_size(itmp) / 8
    if (is_leader) then
       winsize = int(disp_unit, MPI_ADDRESS_KIND)
    else
       winsize = 0_MPI_ADDRESS_KIND
    end if

    ! MPI_Win_allocate_shared is collective over node_comm and returns uniformly within a
    ! node; different nodes are independent (each runs its own per-node window ops below).
    call mpi_win_allocate_shared(winsize, disp_unit, MPI_INFO_NULL, node_comm, baseptr, pwin, ierr)
    have_win = (ierr == MPI_SUCCESS)

    iok(1) = 0
    if (have_win) then
       ! Peers replace their (zero-size) base with the leader's segment.
       if (.not. is_leader) then
          call mpi_win_shared_query(pwin, 0, qsize, qdisp, baseptr, ierr)
          if (ierr /= MPI_SUCCESS) baseptr = c_null_ptr
       end if

       have_ptr = c_associated(baseptr)
       if (have_ptr) call c_f_pointer(baseptr, p, [1])

       ! The fences are collective over node_comm and this whole node has a window, so all
       ! its ranks call them; only the sentinel store/load is guarded by have_ptr.
       call mpi_win_fence(0, pwin, ierr)
       if (have_ptr .and. is_leader) p(1) = SENTINEL
       call mpi_win_fence(0, pwin, ierr)

       if (have_ptr) then
          if (p(1) == SENTINEL) iok(1) = 1
       end if

       call mpi_win_free(pwin, ierr)
    end if

    ! Whole-job agreement -- every rank reaches this (no early return): shared memory is
    ! usable only if EVERY rank mapped and read its node's segment correctly.
    call mpi_allreduce(iok, iok_all, 1, MPI_INTEGER, MPI_MIN, mpicom, ierr)
    ok = (iok_all(1) == 1)
  end function probe_shared_mem
#endif

  !=============================================================================
  subroutine clm_shmem_alloc_i4_1d(ptr, win, n)
    ! Allocate a node-shared default-integer array of length n.  Only the node
    ! leader requests storage; peers map the leader's contiguous segment.
    integer, pointer, intent(out) :: ptr(:)  ! Fortran pointer mapped onto the node-shared buffer
    integer,          intent(out) :: win     ! MPI window handle for the allocation (pass to clm_shmem_free)
    integer,          intent(in)  :: n       ! number of elements to allocate (identical on every rank)

    integer :: istat
#ifndef NO_MPI2
    integer(kind=MPI_ADDRESS_KIND) :: winsize, qsize
    integer :: ierr, disp_unit, qdisp
    integer :: itmp
    type(c_ptr) :: baseptr  ! shared-segment base address; the MPI-3 win routines return the mapped
                            ! address as a C pointer, which c_f_pointer then turns into the pointer ptr
#endif

    call init_comms()

#ifndef NO_MPI2
    if (shared_active) then
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

       ! Guard the base address itself: the associated(ptr) post-condition below cannot
       ! detect a null base because c_f_pointer takes the shape from [n], not from baseptr.
       ! init_comms's probe should already have set shared_active=.false. wherever this
       ! would trip, so reaching this endrun is unexpected.
       if (.not. c_associated(baseptr)) &
            call endrun('clm_shmem_mod: clm_shmem_alloc_i4_1d: shared-memory base pointer is null')
       call c_f_pointer(baseptr, ptr, [n])
    else
       ! Real MPI but shared memory is unavailable on this platform: fall back to a
       ! private per-rank allocation (freed with deallocate; summed over mpicom).
       allocate(ptr(n), stat=istat)
       if (istat /= 0) call endrun('clm_shmem_mod: allocate failed (private fallback path)')
       win = SHMEM_WIN_NONE
    end if
#else
    ! mpi-serial: single task, no shared memory -- a plain local allocation.
    allocate(ptr(n), stat=istat)
    if (istat /= 0) call endrun('clm_shmem_mod: allocate failed (mpi-serial path)')
    win = SHMEM_WIN_NONE
#endif

    ! Post-condition: the allocation must have produced an associated pointer of length n.
    ! A disassociated ptr here means the shared-memory allocation returned a null base
    ! address (e.g. the MPI-3 path being unavailable); catching it now gives a clear error
    ! instead of a confusing "reference to disassociated pointer" at the first use of ptr.
    if (.not. associated(ptr)) call endrun('clm_shmem_mod: clm_shmem_alloc_i4_1d: allocation did not associate ptr')
    if (size(ptr) /= n) call endrun('clm_shmem_mod: clm_shmem_alloc_i4_1d: allocated size does not match n')
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

    ! ptr must be the node-shared buffer of length n created by clm_shmem_alloc_i4_1d.
    ! Guard against a caller passing an inconsistent n: the mpi_allreduce and the ptr(1:n)
    ! store below would otherwise read or write past the end of the buffer.
    if (.not. associated(ptr)) call endrun('clm_shmem_mod: clm_shmem_leader_allreduce_sum_i4: ptr is not associated')
    if (size(ptr) /= n) call endrun('clm_shmem_mod: clm_shmem_leader_allreduce_sum_i4: size(ptr) does not match n')

    call clm_shmem_fence(win)             ! all node stores complete and visible to leader
                                          ! (a no-op in private mode; win == SHMEM_WIN_NONE)
#ifndef NO_MPI2
    if (shared_active) then
       ! Shared mode: the node's partial sum lives in the one shared buffer; the leaders
       ! reduce their per-node partials across nodes over leader_comm.
       if (is_leader) then
          allocate(tmp(n))
          call mpi_allreduce(ptr, tmp, n, MPI_INTEGER, MPI_SUM, leader_comm, ierr)
          if (ierr /= MPI_SUCCESS) call endrun('clm_shmem_mod: MPI_Allreduce failed')
          ptr(1:n) = tmp(1:n)
          deallocate(tmp)
       end if
    else
       ! Private mode: every rank filled disjoint indices of its own full-size copy, so
       ! sum across all ranks over mpicom -- bit-for-bit with the original all-rank reduce.
       allocate(tmp(n))
       call mpi_allreduce(ptr, tmp, n, MPI_INTEGER, MPI_SUM, mpicom, ierr)
       if (ierr /= MPI_SUCCESS) call endrun('clm_shmem_mod: MPI_Allreduce failed')
       ptr(1:n) = tmp(1:n)
       deallocate(tmp)
    end if
#else
    ! mpi-serial: the single task owns the whole domain, so ptr already holds the
    ! global array -- there is nothing to sum across nodes.
#endif
    call clm_shmem_fence(win)             ! publish global result to all node ranks
                                          ! (a no-op in private mode)
  end subroutine clm_shmem_leader_allreduce_sum_i4

  !=============================================================================
  subroutine clm_shmem_fence(win)
    ! Collective over the node communicator; synchronizes the window so stores
    ! become visible to all ranks on the node.  A no-op for the mpi-serial path.
    integer, intent(in) :: win
#ifndef NO_MPI2
    integer :: ierr
    if (win == SHMEM_WIN_NONE) return   ! private-mode allocation: no window to synchronize
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
    if (win == SHMEM_WIN_NONE) then
       ! private-mode fallback: ptr was a plain local allocation, so deallocate it.
       if (associated(ptr)) deallocate(ptr)
    else
       ! shared-memory window: free it (never deallocate a c_f_pointer'd buffer).
       if (win /= MPI_WIN_NULL) then
          call mpi_win_free(win, ierr)
          if (ierr /= MPI_SUCCESS) call endrun('clm_shmem_mod: MPI_Win_free failed')
       end if
       if (associated(ptr)) nullify(ptr)
    end if
    win = SHMEM_WIN_NONE
#else
    ! mpi-serial: ptr was a plain local allocation, so deallocate it.
    if (associated(ptr)) deallocate(ptr)
    win = SHMEM_WIN_NONE
#endif
  end subroutine clm_shmem_free_i4_1d

  !=============================================================================
  logical function clm_shmem_is_leader()
    call init_comms()
    ! In private-fallback mode every rank owns its own full-size copy and must initialize
    ! it, so every rank reports as a leader; in shared mode only node rank 0 does.
    clm_shmem_is_leader = (is_leader .or. .not. shared_active)
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
