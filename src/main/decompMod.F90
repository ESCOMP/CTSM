module decompMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8

  use shr_sys_mod , only : shr_sys_abort ! use shr_sys_abort instead of endrun here to avoid circular dependency
  use clm_varctl  , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none

  ! Define possible bounds subgrid levels
  !
  ! subgrid_level_unspecified can be used in some situations where a subgrid level is
  ! generally needed but it would be hard for the caller to provide it; in this case,
  ! some code that depends on subgrid level will be skipped. (But this should only be
  ! used if you know what you're doing: places where this value is allowed will have a
  ! comment documenting this.)
  integer, parameter, public :: subgrid_level_unspecified = -1
  integer, parameter, public :: subgrid_level_lndgrid     = 0
  integer, parameter, public :: subgrid_level_gridcell    = 1
  integer, parameter, public :: subgrid_level_landunit    = 2
  integer, parameter, public :: subgrid_level_column      = 3
  integer, parameter, public :: subgrid_level_patch       = 4
  integer, parameter, public :: subgrid_level_cohort      = 5

  ! Define possible bounds levels
  integer, parameter, public :: bounds_level_proc  = 1
  integer, parameter, public :: bounds_level_clump = 2
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: get_beg                     ! get beg bound for a given subgrid level
  public :: get_end                     ! get end bound for a given subgrid level
  public :: get_clump_bounds            ! clump beg and end gridcell,landunit,column,patch
  public :: get_proc_bounds             ! this processor beg and end gridcell,landunit,column,patch
  public :: get_proc_total              ! total no. of gridcells, landunits, columns and patchs for any processor
  public :: get_proc_global             ! total gridcells, landunits, columns, patchs across all processors
  public :: get_proc_clumps             ! number of clumps for this processor
  public :: get_global_index            ! Determine global index space value for target point
  public :: get_global_index_array      ! Determine global index space value for target array
  public :: get_subgrid_level_from_name ! Given a name like nameg, return a subgrid level index like subgrid_level_gridcell
  public :: get_subgrid_level_gsize     ! get global size associated with subgrid_level
  public :: get_subgrid_level_gindex    ! get global index array associated with subgrid_level

  ! !PRIVATE MEMBER FUNCTIONS:
  !
  ! !PRIVATE TYPES:
  private  ! (now mostly public for decompinitmod)

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
     integer :: npatches         ! number of patchs in clump
     integer :: nCohorts         ! number of cohorts in proc
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending patch index
     integer :: begCohort, endCohort ! beginning and ending cohort indices
  end type clump_type
  public clump_type
  type(clump_type),public, allocatable :: clumps(:)

  ! ---global sizes
  integer,public :: nclumps          ! total number of clumps across all processors
  integer,public :: numg             ! total number of gridcells on all procs
  integer,public :: numl             ! total number of landunits on all procs
  integer,public :: numc             ! total number of columns on all procs
  integer,public :: nump             ! total number of patchs on all procs
  integer,public :: numCohort        ! total number of fates cohorts on all procs

  ! ---NOTE: the following are allocated with a lower bound of 1!
  integer, public, pointer :: gindex_global(:)   => null() ! includes ocean points
  integer, public, pointer :: gindex_grc(:)      => null() ! does not include ocean points
  integer, public, pointer :: gindex_lun(:)      => null()
  integer, public, pointer :: gindex_col(:)      => null()
  integer, public, pointer :: gindex_patch(:)    => null()
  integer, public, pointer :: gindex_cohort(:)   => null()
  !------------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  pure function get_beg(bounds, subgrid_level) result(beg_index)
    !
    ! !DESCRIPTION:
    ! Get beginning bounds for a given subgrid level
    !
    ! subgrid_level should be one of the constants defined in this module:
    ! subgrid_level_gridcell, subgrid_level_landunit, etc.
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
    case (subgrid_level_gridcell)
       beg_index = bounds%begg
    case (subgrid_level_landunit)
       beg_index = bounds%begl
    case (subgrid_level_column)
       beg_index = bounds%begc
    case (subgrid_level_patch)
       beg_index = bounds%begp
    case (subgrid_level_cohort)
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
    ! subgrid_level_gridcell, subgrid_level_landunit, etc.
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
    case (subgrid_level_gridcell)
       end_index = bounds%endg
    case (subgrid_level_landunit)
       end_index = bounds%endl
    case (subgrid_level_column)
       end_index = bounds%endc
    case (subgrid_level_patch)
       end_index = bounds%endp
    case (subgrid_level_cohort)
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
#endif
    !------------------------------------------------------------------------------
    !    Make sure this IS being called from a threaded region
#ifdef _OPENMP
    if ( OMP_GET_NUM_THREADS() == 1 .and. OMP_GET_MAX_THREADS() > 1 )then
       call shr_sys_abort( trim(subname)//' ERROR: Calling from inside a non-threaded region)')
    end if
#endif

    cid  = procinfo%cid(n)
    bounds%begp      = clumps(cid)%begp - procinfo%begp + 1
    bounds%endp      = clumps(cid)%endp - procinfo%begp + 1
    bounds%begc      = clumps(cid)%begc - procinfo%begc + 1
    bounds%endc      = clumps(cid)%endc - procinfo%begc + 1
    bounds%begl      = clumps(cid)%begl - procinfo%begl + 1
    bounds%endl      = clumps(cid)%endl - procinfo%begl + 1
    bounds%begg      = clumps(cid)%begg - procinfo%begg + 1
    bounds%endg      = clumps(cid)%endg - procinfo%begg + 1
    bounds%begCohort = clumps(cid)%begCohort - procinfo%begCohort + 1
    bounds%endCohort = clumps(cid)%endCohort - procinfo%begCohort + 1

    bounds%level = bounds_level_clump
    bounds%clump_index = n

  end subroutine get_clump_bounds

  !------------------------------------------------------------------------------
  subroutine get_proc_bounds (bounds, allow_call_from_threaded_region)
    !
    ! !DESCRIPTION:
    ! Retrieve processor bounds
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(out) :: bounds ! processor bounds bounds

    ! Normally this routine will abort if it is called from within a threaded region,
    ! because in most cases you should be calling get_clump_bounds in that situation. If
    ! you really want to be using this routine from within a threaded region, then set
    ! allow_call_from_threaded_region to .true.
    logical, intent(in), optional :: allow_call_from_threaded_region
    !
    ! !LOCAL VARIABLES:
    logical :: l_allow_call_from_threaded_region
#ifdef _OPENMP
    integer, external :: OMP_GET_NUM_THREADS
#endif
    character(len=32), parameter :: subname = 'get_proc_bounds'  ! Subroutine name
    !------------------------------------------------------------------------------
    !    Make sure this is NOT being called from a threaded region
    if (present(allow_call_from_threaded_region)) then
       l_allow_call_from_threaded_region = allow_call_from_threaded_region
    else
       l_allow_call_from_threaded_region = .false.
    end if
#ifdef _OPENMP
    if ( OMP_GET_NUM_THREADS() > 1 .and. .not. l_allow_call_from_threaded_region )then
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

    bounds%level = bounds_level_proc
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
  integer function get_global_index(subgrid_index, subgrid_level)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target point at given subgrid level
    !
    ! Uses:
    use shr_log_mod, only: errMsg => shr_log_errMsg
    !
    ! Arguments
    integer , intent(in) :: subgrid_index ! index of interest (can be at any subgrid level or gridcell level)
    integer , intent(in) :: subgrid_level ! one of the subgrid_level_* constants defined above
    !
    ! Local Variables:
    type(bounds_type) :: bounds_proc   ! processor bounds
    integer           :: beg_index     ! beginning proc index for subgrid_level
    integer, pointer  :: gindex(:)
    !----------------------------------------------------------------

    call get_proc_bounds(bounds_proc, allow_call_from_threaded_region=.true.)
    beg_index = get_beg(bounds_proc, subgrid_level)
    if (beg_index == -1) then
       write(iulog,*) 'get_global_index: subgrid_level not supported: ', subgrid_level
       call shr_sys_abort('subgrid_level not supported' // &
            errmsg(sourcefile, __LINE__))
    end if

    call get_subgrid_level_gindex(subgrid_level=subgrid_level, gindex=gindex)
    get_global_index = gindex(subgrid_index - beg_index + 1)

  end function get_global_index

  !-----------------------------------------------------------------------
  function get_global_index_array(subgrid_index, bounds1, bounds2, subgrid_level)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target array at given subgrid level
    !
    ! Example from histFileMod.F90:
    ! ilarr = get_global_index_array(lun%gridcell(bounds%begl:bounds%endl), bounds%begl, bounds%endl, subgrid_level=subgrid_level_gridcell)
    ! Note that the last argument (subgrid_level) is set to nameg, which corresponds
    ! to the "gridcell" not the "lun" of the first argument.
    !
    ! Uses:
#include "shr_assert.h"
    use shr_log_mod, only: errMsg => shr_log_errMsg
    !
    ! Arguments
    integer , intent(in) :: bounds1                 ! lower bound of the input & returned arrays
    integer , intent(in) :: bounds2                 ! upper bound of the input & returned arrays
    integer , intent(in) :: subgrid_index(bounds1:) ! array of indices of interest (can be at any subgrid level or gridcell level)
    integer , intent(in) :: subgrid_level           ! one of the subgrid_level_* constants defined above
    integer              :: get_global_index_array(bounds1:bounds2)
    !
    ! Local Variables:
    type(bounds_type) :: bounds_proc   ! processor bounds
    integer           :: beg_index     ! beginning proc index for subgrid_level
    integer           :: i
    integer , pointer :: gindex(:)
    !----------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(subgrid_index) == (/bounds2/)), sourcefile, __LINE__)
    call get_proc_bounds(bounds_proc, allow_call_from_threaded_region=.true.)
    beg_index = get_beg(bounds_proc, subgrid_level)
    if (beg_index == -1) then
       write(iulog,*) 'get_global_index_array: subgrid_level not supported: ', subgrid_level
       call shr_sys_abort('subgrid_level not supported' // &
            errmsg(sourcefile, __LINE__))
    end if

    call get_subgrid_level_gindex(subgrid_level=subgrid_level, gindex=gindex)
    do i = bounds1, bounds2
       get_global_index_array(i) = gindex(subgrid_index(i) - beg_index + 1)
    end do

  end function get_global_index_array

  !-----------------------------------------------------------------------
  function get_subgrid_level_from_name(subgrid_level_name) result(subgrid_level)
    !
    ! !DESCRIPTION:
    ! Given a name like nameg, return a subgrid level index like subgrid_level_gridcell
    !
    ! !USES:
    use clm_varcon, only : grlnd, nameg, namel, namec, namep, nameCohort
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: subgrid_level_name ! grlnd, nameg, namel, namec, namep, or nameCohort
    integer :: subgrid_level  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_subgrid_level_from_name'
    !-----------------------------------------------------------------------

    select case (subgrid_level_name)
    case(grlnd)
       subgrid_level = subgrid_level_lndgrid
    case(nameg)
       subgrid_level = subgrid_level_gridcell
    case(namel)
       subgrid_level = subgrid_level_landunit
    case(namec)
       subgrid_level = subgrid_level_column
    case(namep)
       subgrid_level = subgrid_level_patch
    case(nameCohort)
       subgrid_level = subgrid_level_cohort
    case default
       write(iulog,*) subname//': unknown subgrid_level_name: ', trim(subgrid_level_name)
       call shr_sys_abort()
    end select

  end function get_subgrid_level_from_name


  !-----------------------------------------------------------------------
  integer function get_subgrid_level_gsize (subgrid_level)
    !
    ! !DESCRIPTION:
    ! Determine 1d size from subgrid_level
    !
    ! !USES:
    use domainMod , only : ldomain
    !
    ! !ARGUMENTS:
    integer, intent(in) :: subgrid_level ! one of the subgrid_level_* constants defined above
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case(subgrid_level_lndgrid)
       get_subgrid_level_gsize = ldomain%ns
    case(subgrid_level_gridcell)
       get_subgrid_level_gsize = numg
    case(subgrid_level_landunit)
       get_subgrid_level_gsize = numl
    case(subgrid_level_column)
       get_subgrid_level_gsize = numc
    case(subgrid_level_patch)
       get_subgrid_level_gsize = nump
    case(subgrid_level_cohort)
       get_subgrid_level_gsize = numCohort
    case default
       write(iulog,*) 'get_subgrid_level_gsize: unknown subgrid_level: ', subgrid_level
       call shr_sys_abort()
    end select

  end function get_subgrid_level_gsize

  !-----------------------------------------------------------------------
  subroutine get_subgrid_level_gindex (subgrid_level, gindex)
    !
    ! !DESCRIPTION:
    ! Get subgrid global index space
    !
    ! !ARGUMENTS:
    integer         , intent(in) :: subgrid_level ! one of the subgrid_level_* constants defined above
    integer         , pointer    :: gindex(:)
    !----------------------------------------------------------------------

    select case (subgrid_level)
    case(subgrid_level_lndgrid)
       gindex => gindex_global
    case(subgrid_level_gridcell)
       gindex => gindex_grc
    case(subgrid_level_landunit)
       gindex => gindex_lun
    case(subgrid_level_column)
       gindex => gindex_col
    case(subgrid_level_patch)
       gindex => gindex_patch
    case(subgrid_level_cohort)
       gindex => gindex_cohort
    case default
       write(iulog,*) 'get_subgrid_level_gindex: unknown subgrid_level: ', subgrid_level
       call shr_sys_abort()
    end select

  end subroutine get_subgrid_level_gindex

end module decompMod
