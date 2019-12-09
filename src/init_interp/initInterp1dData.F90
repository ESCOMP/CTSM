module initInterp1dData

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains routines for interpolating 1-d data fields. These routines do
  ! not do any i/o.
  ! ------------------------------------------------------------------------

#include "shr_assert.h"
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_infnan_mod , only: shr_infnan_isnan
  use clm_varcon     , only: spval, ispval

  implicit none
  private
  save

  ! Public methods

  public :: interp_1d_data
  interface interp_1d_data
     module procedure interp_1d_data_double
     module procedure interp_1d_data_int
  end interface interp_1d_data

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine interp_1d_data_double(begi, endi, bego, endo, sgridindex, keep_existing, &
       data_in, data_out)
    !
    ! !DESCRIPTION:
    ! Interpolate a 1-d double field.
    !
    ! If keep_existing = .true., then points with sgridindex <= 0 are kept at their
    ! original value; if .false., they are set to spval.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: begi              ! beginning index for input array
    integer  , intent(in)    :: endi              ! ending index for input array
    integer  , intent(in)    :: bego              ! beginning index for output array
    integer  , intent(in)    :: endo              ! ending index for output array
    integer  , intent(in)    :: sgridindex(bego:) ! input index mapping to each outputpoint
    logical  , intent(in)    :: keep_existing     ! whether to keep existing values for points with sgridindex <= 0
    real(r8) , intent(in)    :: data_in(begi:)    ! input data
    real(r8) , intent(inout) :: data_out(bego:)   ! output data
    !
    ! !LOCAL VARIABLES:
    integer :: no,ni          ! indices

    character(len=*), parameter :: subname = 'interp_1d_data_double'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(sgridindex) == (/endo/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(data_in) == (/endi/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(data_out) == (/endo/)), sourcefile, __LINE__)

    if (.not. keep_existing) then
       data_out(bego:endo) = spval
    end if

    do no = bego,endo
       ni = sgridindex(no)
       if (ni > 0) then
          if ( shr_infnan_isnan(data_in(ni)) ) then
             data_out(no) = spval
          else
             data_out(no) = data_in(ni)
          end if
       end if
    end do

  end subroutine interp_1d_data_double

  !-----------------------------------------------------------------------
  subroutine interp_1d_data_int(begi, endi, bego, endo, sgridindex, keep_existing, &
       data_in, data_out)
    !
    ! !DESCRIPTION:
    ! Interpolate a 1-d int field
    !
    ! If keep_existing = .true., then points with sgridindex <= 0 are kept at their
    ! original value; if .false., they are set to ispval.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer , intent(in)    :: begi              ! beginning index for input array
    integer , intent(in)    :: endi              ! ending index for input array
    integer , intent(in)    :: bego              ! beginning index for output array
    integer , intent(in)    :: endo              ! ending index for output array
    integer , intent(in)    :: sgridindex(bego:) ! input index mapping to each outputpoint
    logical , intent(in)    :: keep_existing     ! whether to keep existing values for points with sgridindex <= 0
    integer , intent(in)    :: data_in(begi:)    ! input data
    integer , intent(inout) :: data_out(bego:)   ! output data
    !
    ! !LOCAL VARIABLES:
    integer           :: no,ni          !indices

    character(len=*), parameter :: subname = 'interp_1d_data_int'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(sgridindex) == (/endo/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(data_in) == (/endi/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(data_out) == (/endo/)), sourcefile, __LINE__)

    if (.not. keep_existing) then
       data_out(bego:endo) = ispval
    end if

    do no = bego,endo
       ni = sgridindex(no)
       if (ni > 0) then
          data_out(no) = data_in(ni)  
       end if
    end do

  end subroutine interp_1d_data_int

end module initInterp1dData
