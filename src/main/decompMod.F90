module decompMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  ! Must use shr_sys_abort rather than endrun here to avoid circular dependency
  use shr_sys_mod , only : shr_sys_abort
  use clm_varctl  , only : iulog
  use clm_varcon  , only : grlnd, nameg, namel, namec, namep, nameCohort
  use mct_mod     , only : mct_gsMap
  !
  ! !PUBLIC TYPES:
  implicit none

  ! mct data type still needed for determining subgrid gindex
  type(mct_gsMap), target, public  :: gsmap_global  ! global seg map

  integer, public :: clump_pproc ! number of clumps per MPI process

  ! Define possible bounds subgrid levels
  integer, parameter, public :: BOUNDS_SUBGRID_GRIDCELL = 1
  integer, parameter, public :: BOUNDS_SUBGRID_LANDUNIT = 2
  integer, parameter, public :: BOUNDS_SUBGRID_COLUMN   = 3
  integer, parameter, public :: BOUNDS_SUBGRID_PATCH    = 4
  integer, parameter, public :: BOUNDS_SUBGRID_COHORT   = 5

  ! Define possible bounds levels
  integer, parameter, public :: BOUNDS_LEVEL_PROC  = 1
  integer, parameter, public :: BOUNDS_LEVEL_CLUMP = 2
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public get_beg            ! get beg bound for a given subgrid level
  public get_end            ! get end bound for a given subgrid level
  public get_proc_clumps    ! number of clumps for this processor
  public get_proc_total     ! total no. of gridcells, landunits, columns and patchs for any processor
  public get_proc_global    ! total gridcells, landunits, columns, patchs across all processors
  public get_clmlevel_gsize ! get global size associated with clmlevel
  public get_clmlevel_gindex! get global size associated with clmlevel
  public get_clump_bounds   ! clump beg and end gridcell,landunit,column,patch
  public get_proc_bounds    ! this processor beg and end gridcell,landunit,column,patch

  ! !PRIVATE MEMBER FUNCTIONS:
  !
  ! !PRIVATE TYPES:
  private  ! (now mostly public for decompinitmod)

  integer,public :: nclumps          ! total number of clumps across all processors
  integer,public :: numg             ! total number of gridcells on all procs
  integer,public :: numl             ! total number of landunits on all procs
  integer,public :: numc             ! total number of columns on all procs
  integer,public :: nump             ! total number of patchs on all procs
  integer,public :: numCohort        ! total number of fates cohorts on all procs

  type bounds_type
     integer :: begg, endg           ! beginning and ending gridcell index
     integer :: begl, endl           ! beginning and ending landunit index
     integer :: begc, endc           ! beginning and ending column index
     integer :: begp, endp           ! beginning and ending patch index
     integer :: begCohort, endCohort ! beginning and ending cohort indices
     integer :: level                ! whether defined on the proc or clump level
     integer :: clump_index          ! if defined on the clump level, this gives the clump index
  end type bounds_type
  public :: bounds_type

  !---global information on each pe
  type processor_type
     integer :: nclumps              ! number of clumps for processor_type iam
     integer,pointer :: cid(:)       ! clump indices
     integer :: ncells               ! number of gridcells in proc
     integer :: nlunits              ! number of landunits in proc
     integer :: ncols                ! number of columns in proc
     integer :: npatches             ! number of patchs in proc
     integer :: nCohorts             ! number of cohorts in proc
     integer :: begg, endg           ! beginning and ending gridcell index
     integer :: begl, endl           ! beginning and ending landunit index
     integer :: begc, endc           ! beginning and ending column index
     integer :: begp, endp           ! beginning and ending patch index
     integer :: begCohort, endCohort ! beginning and ending cohort indices
  end type processor_type
  public processor_type
  type(processor_type),public :: procinfo

  !---global information on each pe
  type clump_type
     integer :: owner            ! process id owning clump
     integer :: ncells           ! number of gridcells in clump
     integer :: nlunits          ! number of landunits in clump
     integer :: ncols            ! number of columns in clump
     integer :: npatches          ! number of patchs in clump
     integer :: nCohorts         ! number of cohorts in proc
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending patch index
     integer :: begCohort, endCohort ! beginning and ending cohort indices
  end type clump_type
  public clump_type
  type(clump_type),public, allocatable :: clumps(:)

  ! NOTE: the following are allocated with a lower bound of 1!
  integer, public, pointer :: gindex_global(:)   => null() ! includes ocean points
  integer, public, pointer :: gindex_grc(:)      => null() ! does not include ocean points
  integer, public, pointer :: gindex_lun(:)      => null()
  integer, public, pointer :: gindex_col(:)      => null()
  integer, public, pointer :: gindex_patch(:)    => null()
  integer, public, pointer :: gindex_cohort(:)   => null()
  integer, public, pointer :: gindex_lnd2Dsoi(:) => null()
  integer, public          :: nglob_x, nglob_y  ! global sizes
  !------------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  pure function get_beg(bounds, subgrid_level) result(beg_index)
    !
    ! !DESCRIPTION:
    ! Get beginning bounds for a given subgrid level
    !
    ! subgrid_level should be one of the constants defined in this module:
    ! BOUNDS_SUBGRID_GRIDCELL, BOUNDS_SUBGRID_LANDUNIT, etc.
    !
    ! Returns -1 for invalid subgrid_level (does not abort in this case, in order to keep
    ! this function pure).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: beg_index  ! function result
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: subgrid_level
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_beg'
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case (BOUNDS_SUBGRID_GRIDCELL)
       beg_index = bounds%begg
    case (BOUNDS_SUBGRID_LANDUNIT)
       beg_index = bounds%begl
    case (BOUNDS_SUBGRID_COLUMN)
       beg_index = bounds%begc
    case (BOUNDS_SUBGRID_PATCH)
       beg_index = bounds%begp
    case (BOUNDS_SUBGRID_COHORT)
       beg_index = bounds%begCohort
    case default
       beg_index = -1
    end select

  end function get_beg

  !-----------------------------------------------------------------------
  pure function get_end(bounds, subgrid_level) result(end_index)
    !
    ! !DESCRIPTION:
    ! Get end bounds for a given subgrid level
    !
    ! subgrid_level should be one of the constants defined in this module:
    ! BOUNDS_SUBGRID_GRIDCELL, BOUNDS_SUBGRID_LANDUNIT, etc.
    !
    ! Returns -1 for invalid subgrid_level (does not abort in this case, in order to keep
    ! this function pure).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer                        :: end_index  ! function result
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: subgrid_level
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'get_end'
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case (BOUNDS_SUBGRID_GRIDCELL)
       end_index = bounds%endg
    case (BOUNDS_SUBGRID_LANDUNIT)
       end_index = bounds%endl
    case (BOUNDS_SUBGRID_COLUMN)
       end_index = bounds%endc
    case (BOUNDS_SUBGRID_PATCH)
       end_index = bounds%endp
    case (BOUNDS_SUBGRID_COHORT)
       end_index = bounds%endCohort
    case default
       end_index = -1
    end select

  end function get_end

  !------------------------------------------------------------------------------
   subroutine get_clump_bounds (n, bounds)
     !
     ! !DESCRIPTION:
     ! Determine clump bounds
     !
     ! !ARGUMENTS:
     integer, intent(in)  :: n                ! processor clump index
     type(bounds_type), intent(out) :: bounds ! clump bounds
     !
     ! !LOCAL VARIABLES:
     character(len=32), parameter :: subname = 'get_clump_bounds'  ! Subroutine name
     integer :: cid                                                ! clump id
#ifdef _OPENMP
     integer, external :: OMP_GET_MAX_THREADS
     integer, external :: OMP_GET_NUM_THREADS
     integer, external :: OMP_GET_THREAD_NUM
#endif
     !------------------------------------------------------------------------------
     !    Make sure this IS being called from a threaded region
#ifdef _OPENMP
     ! FIX(SPM, 090314) - for debugging fates and openMP
     !write(iulog,*) 'SPM omp debug decompMod 1 ', &
          !OMP_GET_NUM_THREADS(),OMP_GET_MAX_THREADS(),OMP_GET_THREAD_NUM()

     if ( OMP_GET_NUM_THREADS() == 1 .and. OMP_GET_MAX_THREADS() > 1 )then
        call shr_sys_abort( trim(subname)//' ERROR: Calling from inside a non-threaded region)')
     end if
#endif

     cid  = procinfo%cid(n)
     bounds%begp = 1
     bounds%endp = clumps(cid)%endp - clumps(cid)%begp + 1
     bounds%begc = 1
     bounds%endc = clumps(cid)%endc - clumps(cid)%begc + 1
     bounds%begl = 1
     bounds%endl = clumps(cid)%endl - clumps(cid)%begl + 1
     bounds%begg = 1
     bounds%endg = clumps(cid)%endg - clumps(cid)%begg + 1
     bounds%begCohort = 1
     bounds%endCohort = clumps(cid)%endCohort - clumps(cid)%begCohort + 1

     bounds%level = BOUNDS_LEVEL_CLUMP
     bounds%clump_index = n

   end subroutine get_clump_bounds

   !------------------------------------------------------------------------------
   subroutine get_proc_bounds (bounds)
     !
     ! !DESCRIPTION:
     ! Retrieve processor bounds
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(out) :: bounds ! processor bounds bounds
     !
     ! !LOCAL VARIABLES:
#ifdef _OPENMP
     integer, external :: OMP_GET_NUM_THREADS
     integer, external :: OMP_GET_MAX_THREADS
     integer, external :: OMP_GET_THREAD_NUM
#endif
     character(len=32), parameter :: subname = 'get_proc_bounds'  ! Subroutine name
     !------------------------------------------------------------------------------
     !    Make sure this is NOT being called from a threaded region
#ifdef _OPENMP
     ! FIX(SPM, 090314) - for debugging fates and openMP
     !write(*,*) 'SPM omp debug decompMod 2 ', &
          !OMP_GET_NUM_THREADS(),OMP_GET_MAX_THREADS(),OMP_GET_THREAD_NUM()
     if ( OMP_GET_NUM_THREADS() > 1 )then
        call shr_sys_abort( trim(subname)//' ERROR: Calling from inside  a threaded region')
     end if
#endif

     bounds%begp = 1
     bounds%endp = procinfo%endp - procinfo%begp + 1
     bounds%begc = 1
     bounds%endc = procinfo%endc - procinfo%begc + 1
     bounds%begl = 1
     bounds%endl = procinfo%endl - procinfo%begl + 1
     bounds%begg = 1
     bounds%endg = procinfo%endg - procinfo%begg + 1
     bounds%begCohort = 1
     bounds%endCohort = procinfo%endCohort - procinfo%begCohort + 1

     bounds%level = BOUNDS_LEVEL_PROC
     bounds%clump_index = -1           ! irrelevant for proc, so assigned a bogus value

   end subroutine get_proc_bounds

   !------------------------------------------------------------------------------
   subroutine get_proc_total(pid, ncells, nlunits, ncols, npatches, nCohorts)
     !
     ! !DESCRIPTION:
     ! Count up gridcells, landunits, columns, and patchs on process.
     !
     ! !ARGUMENTS:
     integer, intent(in)  :: pid     ! proc id
     integer, intent(out) :: ncells  ! total number of gridcells on the processor
     integer, intent(out) :: nlunits ! total number of landunits on the processor
     integer, intent(out) :: ncols   ! total number of columns on the processor
     integer, intent(out) :: npatches ! total number of patchs on the processor
     integer, intent(out) :: nCohorts! total number of cohorts on the processor
     !
     ! !LOCAL VARIABLES:
     integer :: cid       ! clump index
     !------------------------------------------------------------------------------

     npatches  = 0
     nlunits  = 0
     ncols    = 0
     ncells   = 0
     nCohorts = 0
     do cid = 1,nclumps
        if (clumps(cid)%owner == pid) then
           ncells  = ncells    + clumps(cid)%ncells
           nlunits = nlunits   + clumps(cid)%nlunits
           ncols   = ncols     + clumps(cid)%ncols
           npatches = npatches   + clumps(cid)%npatches
           nCohorts = nCohorts + clumps(cid)%nCohorts
        end if
     end do
   end subroutine get_proc_total

   !------------------------------------------------------------------------------
   subroutine get_proc_global(ng, nl, nc, np, nCohorts)
     !
     ! !DESCRIPTION:
     ! Return number of gridcells, landunits, columns, and patchs across all processes.
     !
     ! !ARGUMENTS:
     integer, optional, intent(out) :: ng        ! total number of gridcells across all processors
     integer, optional, intent(out) :: nl        ! total number of landunits across all processors
     integer, optional, intent(out) :: nc        ! total number of columns across all processors
     integer, optional, intent(out) :: np        ! total number of patchs across all processors
     integer, optional, intent(out) :: nCohorts  ! total number fates cohorts
     !------------------------------------------------------------------------------

     if (present(np)) np             = nump
     if (present(nc)) nc             = numc
     if (present(nl)) nl             = numl
     if (present(ng)) ng             = numg
     if (present(nCohorts)) nCohorts = numCohort

   end subroutine get_proc_global

   !------------------------------------------------------------------------------
   integer function get_proc_clumps()
     !
     ! !DESCRIPTION:
     ! Return the number of clumps.
     !------------------------------------------------------------------------------

     get_proc_clumps = procinfo%nclumps

   end function get_proc_clumps

   !-----------------------------------------------------------------------
   integer function get_clmlevel_gsize (clmlevel)
     !
     ! !DESCRIPTION:
     ! Determine 1d size from clmlevel
     !
     ! !USES:
     use domainMod , only : ldomain
     !
     ! !ARGUMENTS:
     character(len=*), intent(in) :: clmlevel    !type of clm 1d array
     !-----------------------------------------------------------------------

     select case (clmlevel)
     case(grlnd)
        get_clmlevel_gsize = ldomain%ns
     case(nameg)
        get_clmlevel_gsize = numg
     case(namel)
        get_clmlevel_gsize = numl
     case(namec)
        get_clmlevel_gsize = numc
     case(namep)
        get_clmlevel_gsize = nump
     case(nameCohort)
        get_clmlevel_gsize = numCohort
     case default
        write(iulog,*) 'get_clmlevel_gsize does not match clmlevel type: ', trim(clmlevel)
        call shr_sys_abort()
     end select

   end function get_clmlevel_gsize

   !-----------------------------------------------------------------------
   subroutine get_clmlevel_gindex (clmlevel, gindex)
     !
     ! !DESCRIPTION:
     ! Compute arguments for gatherv, scatterv for vectors
     !
     ! !ARGUMENTS:
     character(len=*), intent(in) :: clmlevel     ! type of input data
     integer, pointer :: gindex(:)
     !----------------------------------------------------------------------

    select case (clmlevel)
    case(grlnd)
       gindex => gindex_global
    case(nameg)
       gindex => gindex_grc
    case(namel)
       gindex => gindex_lun
    case(namec)
       gindex => gindex_col
    case(namep)
       gindex => gindex_patch
    case(nameCohort)
       gindex => gindex_cohort
    case default
       write(iulog,*) 'get_clmlevel_gindex: Invalid expansion character: ',trim(clmlevel)
       call shr_sys_abort()
    end select

  end subroutine get_clmlevel_gindex

end module decompMod
