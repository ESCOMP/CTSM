module CropReprPoolsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains variables giving information about the different crop-specific
  ! reproductive pools, and subroutines for working with these variables.
  !
  ! For most purposes, reproductive pool/flux variables are allocated with an extra
  ! dimension that goes from 1:nrepr, and you will want to loop over all of these indices.
  ! However, the reproductive components can be further subdivided into grain vs.
  ! structural reproductive components, and there are some parts of the code where you
  ! just want to work with one or the other of these sets (i.e., just the grain components
  ! or just the structural reproductive components). For these situations you can loop
  ! over repr_grain_min:repr_grain_max, or repr_structure_min:repr_structure_max. The
  ! naming convention of science (state/flux) variables is: Variables that spell out the
  ! full "reproductive" in their names (e.g., reproductivec_patch) are allocated from
  ! 1:nrepr; variables with names containing "repr_grain" are allocated from
  ! repr_grain_min:repr_grain_max; variables with names containing "repr_structure" are
  ! allocated from repr_structure_min:repr_structure_max.
  !
  ! !USES:
#include "shr_assert.h"
  use clm_varctl, only : for_testing_use_second_grain_pool, for_testing_use_repr_structure_pool
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: crop_repr_pools_init ! initialize module-level data
  public :: get_repr_hist_fname  ! get name of a given reproductive pool for use in history field names
  public :: get_repr_rest_fname  ! get name of a given reproductive pool for use in restart field names
  public :: get_repr_longname    ! get name of a given reproductive pool for use in history/restart long names

  !
  ! !PUBLIC DATA:
  integer, public :: nrepr              ! number of crop reproductive pools (grain & structural reproductive components)
  integer, public :: nrepr_grain        ! number of reproductive grain pools
  integer, public :: nrepr_structure    ! number of reproductive structure pools
  integer, public :: repr_grain_min     ! min index for reproductive grain pools
  integer, public :: repr_grain_max     ! max index for reproductive grain pools
  integer, public :: repr_structure_min ! min index for reproductive structure pools
  integer, public :: repr_structure_max ! max index for reproductive structure pools

  !
  ! !PRIVATE DATA:
  character(len=64), allocatable, private :: repr_hist_fnames(:)  ! names of each reproductive pool for use in history field names
  character(len=64), allocatable, private :: repr_rest_fnames(:)  ! names of each reproductive pool for use in restart field names
  character(len=64), allocatable, private :: repr_longnames(:)    ! names of each reproductive pool for use in history / restart long names

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  !-----------------------------------------------------------------------
  subroutine crop_repr_pools_init()
    !
    ! !DESCRIPTION:
    ! Initialize module-level data
    !
    ! !USES:
    use abortutils, only: endrun
    use shr_log_mod, only: errmsg => shr_log_errMsg
    use CNSharedParamsMod, only: use_matrixcn
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: k

    character(len=*), parameter :: subname = 'crop_repr_pools_init'
    !-----------------------------------------------------------------------

    ! NOTE(wjs, 2022-03-03) The following variables are set up for the current,
    ! AgroIBIS-based crop model, which has a single grain pool and no reproductive
    ! structure pools. The details of how this is done (e.g., the pool called 'GRAIN' is
    ! the nrepr'th pool rather than the first pool) support software testing where we have
    ! two grain pools instead of one, in which situation the second grain pool acts like
    ! the normal single grain pool.
    !
    ! TODO(wjs, 2022-02-14) This will eventually depend on whether we're using AgSys:
    ! AgSys will use nrepr = 4, nrepr_grain = 2, nrepr_struture = 2,
    ! repr_hist_fnames(1) = 'GRAIN_MEAL', grain_hist_fnames(2) = 'GRAIN_OIL', etc.
    if (for_testing_use_second_grain_pool) then
       nrepr_grain = 2
       if (use_matrixcn) then
          call endrun(msg="ERROR: for_testing_use_second_grain_pool should be .false. when use_matrixcn = .true."//errmsg(sourcefile, __LINE__))
       end if
    else
       nrepr_grain = 1
    end if
    if (for_testing_use_repr_structure_pool) then
       nrepr_structure = 2
       if (use_matrixcn) then
          call endrun(msg="ERROR: for_testing_use_repr_structure_pool should be .false. when use_matrixcn = .true."//errMsg(sourcefile, __LINE__))
       end if
    else
       nrepr_structure = 0
    end if

    nrepr = nrepr_grain + nrepr_structure
    ! matrixcn works with nrepr = 1 only
    if (use_matrixcn .and. nrepr /= 1) then
       call endrun(msg="ERROR: nrepr should be 1 when use_matrixcn = .true."//errMsg(sourcefile, __LINE__))
    end if
    allocate(repr_hist_fnames(nrepr))
    allocate(repr_rest_fnames(nrepr))
    allocate(repr_longnames(nrepr))
    ! The following is a little more convoluted than it needs to be with a single pool in
    ! order to support software testing where we have two grain pools and want the second
    ! grain pool to act like the normal single grain pool. (In the typical case with a
    ! single pool, the do loop is never entered and nrepr = 1, so the following block
    ! is equivalent to setting the first element of the repr names to 'GRAIN'/'grain'.)
    do k = 1, nrepr-1
       write(repr_hist_fnames(k), '(a, i0)') 'REPRODUCTIVE', k
       write(repr_rest_fnames(k), '(a, i0)') 'reproductive', k
       write(repr_longnames(k),   '(a, i0)') 'reproductive', k
    end do
    repr_hist_fnames(nrepr) = 'GRAIN'
    repr_rest_fnames(nrepr) = 'grain'
    repr_longnames(nrepr)   = 'grain'

    ! NOTE(wjs, 2022-03-03) In contrast to the above code, the code below this point is
    ! general enough that it should work for AgSys as well as the current crop model
    ! without any modifications.

    repr_grain_min = 1
    repr_grain_max = repr_grain_min + nrepr_grain - 1
    repr_structure_min = repr_grain_max + 1
    repr_structure_max = repr_structure_min + nrepr_structure - 1

    SHR_ASSERT_FL((nrepr_grain + nrepr_structure == nrepr), sourcefile, __LINE__)
    SHR_ASSERT_FL((repr_structure_max == nrepr), sourcefile, __LINE__)

  end subroutine crop_repr_pools_init

  !-----------------------------------------------------------------------
  function get_repr_hist_fname(repr_pool_num) result(fname)
    !
    ! !DESCRIPTION:
    ! Get trimmed name of a given reproductive pool for use in history field names
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: fname  ! function result
    integer, intent(in) :: repr_pool_num
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_repr_hist_fname'
    !-----------------------------------------------------------------------

    fname = trim(repr_hist_fnames(repr_pool_num))
  end function get_repr_hist_fname

  !-----------------------------------------------------------------------
  function get_repr_rest_fname(repr_pool_num) result(fname)
    !
    ! !DESCRIPTION:
    ! Get trimmed name of a given reproductive pool for use in restart field names
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: fname  ! function result
    integer, intent(in) :: repr_pool_num
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_repr_rest_fname'
    !-----------------------------------------------------------------------

    fname = trim(repr_rest_fnames(repr_pool_num))
  end function get_repr_rest_fname

  !-----------------------------------------------------------------------
  function get_repr_longname(repr_pool_num) result(longname)
    !
    ! !DESCRIPTION:
    ! Get trimmed name of a given reproductive pool for use in history / restart long names
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: longname  ! function result
    integer, intent(in) :: repr_pool_num
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_repr_longname'
    !-----------------------------------------------------------------------

    longname = trim(repr_longnames(repr_pool_num))
  end function get_repr_longname

end module CropReprPoolsMod
