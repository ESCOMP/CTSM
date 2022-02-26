module CropPoolsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains variables giving information about the different crop-specific
  ! pools (e.g., the number of reproductive grain and structure pools), and subroutines
  ! for working with these variables.
  !
  ! !USES:
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: crop_pools_init      ! initialize module-level data
  public :: get_grain_hist_fname ! get name of a given grain pool for use in history field names
  public :: get_grain_rest_fname ! get name of a given grain pool for use in restart field names
  public :: get_grain_longname   ! get name of a given grain pool for use in history/restart long names

  !
  ! !PUBLIC DATA:
  integer, public :: ngrain ! number of crop reproductive grain pools

  !
  ! !PRIVATE DATA:
  character(len=64), allocatable, private :: grain_hist_fnames(:)  ! names of each grain pool for use in history field names
  character(len=64), allocatable, private :: grain_rest_fnames(:)  ! names of each grain pool for use in restart field names
  character(len=64), allocatable, private :: grain_longnames(:)    ! names of each grain pool for use in history / restart long names

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  !-----------------------------------------------------------------------
  subroutine crop_pools_init()
    !
    ! !DESCRIPTION:
    ! Initialize module-level data
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'crop_pools_init'
    !-----------------------------------------------------------------------

    ! TODO(wjs, 2022-02-14) This will eventually depend on whether we're using AgSys:
    ! AgSys will use ngrain = 2, grain_hist_fnames(1) = 'GRAIN_MEAL', grain_hist_fnames(2) = 'GRAIN_OIL', etc.
    ngrain = 1
    allocate(grain_hist_fnames(ngrain))
    allocate(grain_rest_fnames(ngrain))
    allocate(grain_longnames(ngrain))
    grain_hist_fnames(1) = 'GRAIN'
    grain_rest_fnames(1) = 'reproductive_grain'
    grain_longnames(1) = 'grain'

  end subroutine crop_pools_init

  !-----------------------------------------------------------------------
  function get_grain_hist_fname(grain_pool_num) result(fname)
    !
    ! !DESCRIPTION:
    ! Get trimmed name of a given grain pool for use in history field names
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: fname  ! function result
    integer, intent(in) :: grain_pool_num
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_grain_hist_fname'
    !-----------------------------------------------------------------------

    fname = trim(grain_hist_fnames(grain_pool_num))
  end function get_grain_hist_fname

  !-----------------------------------------------------------------------
  function get_grain_rest_fname(grain_pool_num) result(fname)
    !
    ! !DESCRIPTION:
    ! Get trimmed name of a given grain pool for use in restart field names
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: fname  ! function result
    integer, intent(in) :: grain_pool_num
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_grain_rest_fname'
    !-----------------------------------------------------------------------

    fname = trim(grain_rest_fnames(grain_pool_num))
  end function get_grain_rest_fname

  !-----------------------------------------------------------------------
  function get_grain_longname(grain_pool_num) result(longname)
    !
    ! !DESCRIPTION:
    ! Get trimmed name of a given grain pool for use in history / restart long names
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: longname  ! function result
    integer, intent(in) :: grain_pool_num
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_grain_longname'
    !-----------------------------------------------------------------------

    longname = trim(grain_longnames(grain_pool_num))
  end function get_grain_longname

end module CropPoolsMod
