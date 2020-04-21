module NumericsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Utility routines for assisting with model numerics
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8

  implicit none
  save
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: truncate_small_values         ! Truncate relatively small values to 0
  public :: truncate_small_values_one_lev ! Truncate relatively small values to 0, for one level of a multi-level field

  ! !PUBLIC MEMBER DATA:

  ! Relative differences below rel_epsilon are considered to be zero.
  !
  ! Note that double precision machine epsilon is approximately 1e-16, so this value of
  ! 1e-13 allows for 3 orders of magnitude of "slop".
  !
  ! Examples of how to use this:
  !
  ! (1) Rather than checking
  !        if (x == y)
  !     instead check
  !        if (abs(x - y) < rel_epsilon * x)
  !     or
  !        if (abs(x - y) < rel_epsilon * y)
  !
  ! (2) After a state update, you can truncate the state to 0 based on this condition:
  !        if (abs(some_state) < rel_epsilon * abs(some_state_orig)) then
  !           some_state = 0._r8
  !        end if
  !     where some_state_orig is the value of the state variable before the update
  real(r8), public, parameter :: rel_epsilon = 1.e-13_r8  ! Relative differences below this are considered to be zero

  ! !PRIVATE MEMBER DATA:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine truncate_small_values(num_f, filter_f, lb, ub, data_baseline, data)
    !
    ! !DESCRIPTION:
    ! Truncate relatively small values to 0, within the given filter.
    !
    ! "Relatively small" is determined by comparison with some "baseline" version of the
    ! data.
    !
    ! For example, this can be used after doing a state update. In this case,
    ! data_baseline should hold the values before the state update, and data should hold
    ! the values after the state update.
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: num_f              ! number of points in filter_f
    integer  , intent(in)    :: filter_f(:)        ! filter of points in data
    integer  , intent(in)    :: lb                 ! lower bound of data
    integer  , intent(in)    :: ub                 ! upper bound of data
    real(r8) , intent(in)    :: data_baseline(lb:) ! baseline version of data, used to define "relatively close to 0"
    real(r8) , intent(inout) :: data(lb:)          ! data to operate on
    !
    ! !LOCAL VARIABLES:
    integer :: fn  ! index into filter
    integer :: n   ! index into data

    character(len=*), parameter :: subname = 'truncate_small_values'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(data_baseline) == (/ub/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(data) == (/ub/)), sourcefile, __LINE__)

    do fn = 1, num_f
       n = filter_f(fn)
       if (abs(data(n)) < rel_epsilon * abs(data_baseline(n))) then
          data(n) = 0._r8
       end if
    end do

  end subroutine truncate_small_values

  !-----------------------------------------------------------------------
  subroutine truncate_small_values_one_lev(num_f, filter_f, lb, ub, lev_lb, lev, data_baseline, data)
    !
    ! !DESCRIPTION:
    ! Truncate relatively small values to 0, for one level of a multi-level field, within
    ! the given filter.
    !
    ! This version works on one level of a multi-level field. The level of interest can
    ! differ for each point, and is given by the 'lev' input. For point n, we operate
    ! just on data(n,lev(n)).
    !
    ! "Relatively small" is determined by comparison with some "baseline" version of the
    ! data.
    !
    ! For example, this can be used after doing a state update. In this case,
    ! data_baseline should hold the values before the state update, and data should hold
    ! the values after the state update.
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: num_f              ! number of points in filter_f
    integer  , intent(in)    :: filter_f(:)        ! filter of points in data
    integer  , intent(in)    :: lb                 ! lower bound of data
    integer  , intent(in)    :: ub                 ! upper bound of data
    integer  , intent(in)    :: lev_lb             ! lower bound of level dimension of data
    integer  , intent(in)    :: lev(lb:)           ! for each point, which level to work on
    real(r8) , intent(in)    :: data_baseline(lb:) ! baseline version of data, used to define "relatively close to 0" (note that this is only 1-d, giving baselines for the appropriate level for each point)
    real(r8) , intent(inout) :: data(lb:, lev_lb:) ! data to operate on
    !
    ! !LOCAL VARIABLES:
    integer :: fn  ! index into filter
    integer :: n   ! index into data

    character(len=*), parameter :: subname = 'truncate_small_values_one_lev'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(lev, 1) == ub), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(data_baseline, 1) == ub), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(data, 1) == ub), sourcefile, __LINE__)

    do fn = 1, num_f
       n = filter_f(fn)
       if (abs(data(n, lev(n))) < rel_epsilon * abs(data_baseline(n))) then
          data(n, lev(n)) = 0._r8
       end if
    end do

  end subroutine truncate_small_values_one_lev

end module NumericsMod
