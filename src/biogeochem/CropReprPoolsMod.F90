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
  ! over repr_grain_min:repr_grain_max, or repr_structure_min:repr_structure_max.
  !
  ! !USES:
#include "shr_assert.h"
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
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'crop_repr_pools_init'
    !-----------------------------------------------------------------------

    ! TODO(wjs, 2022-02-14) This will eventually depend on whether we're using AgSys:
    ! AgSys will use nrepr = 4, nrepr_grain = 2, nrepr_struture = 2,
    ! repr_hist_fnames(1) = 'GRAIN_MEAL', grain_hist_fnames(2) = 'GRAIN_OIL', etc.
    nrepr = 1
    nrepr_grain = 1
    nrepr_structure = 0
    allocate(repr_hist_fnames(nrepr))
    allocate(repr_rest_fnames(nrepr))
    allocate(repr_longnames(nrepr))
    repr_hist_fnames(1) = 'GRAIN'
    repr_rest_fnames(1) = 'grain'
    repr_longnames(1) = 'grain'

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
