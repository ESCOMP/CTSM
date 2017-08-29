module NumericsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Utility routines for assisting with model numerics
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use clm_varcon       , only : rel_epsilon

  implicit none
  save
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: truncate_small_values  ! Truncate relatively small values to 0

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

    SHR_ASSERT_ALL((ubound(data_baseline) == (/ub/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(data) == (/ub/)), errMsg(sourcefile, __LINE__))

    do fn = 1, num_f
       n = filter_f(fn)
       if (abs(data(n)) < rel_epsilon * abs(data_baseline(n))) then
          data(n) = 0._r8
       end if
    end do

  end subroutine truncate_small_values

end module NumericsMod
