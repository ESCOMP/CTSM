module WaterTracerUtils

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Utility routines for working with water tracers
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use clm_varctl     , only : iulog
  use abortutils     , only : endrun
  use shr_infnan_mod , only : shr_infnan_isnan
  use shr_log_mod    , only : errMsg => shr_log_errMsg

  implicit none
  save
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CompareBulkToTracer

  ! !PRIVATE MEMBER DATA:
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine CompareBulkToTracer(bounds_beg, bounds_end, &
       bulk, tracer, caller_location, name, level)
    !
    ! !DESCRIPTION:
    ! Compare the bulk and tracer quantities; abort if they differ
    !
    ! !ARGUMENTS:
    integer, intent(in) :: bounds_beg
    integer, intent(in) :: bounds_end
    real(r8), intent(in) :: bulk( bounds_beg: )
    real(r8), intent(in) :: tracer( bounds_beg: )
    character(len=*), intent(in) :: caller_location  ! brief description of where this is called from (for error messages)
    character(len=*), intent(in) :: name             ! variable name (for error messages)
    integer, intent(in), optional :: level           ! level being compared (for error messages)
    !
    ! !LOCAL VARIABLES:
    logical :: arrays_equal
    integer :: i
    integer :: diffloc
    real(r8) :: val1, val2

    real(r8), parameter :: tolerance = 1.0e-7_r8

    character(len=*), parameter :: subname = 'CompareBulkToTracer'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(bulk) == [bounds_end]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(tracer) == [bounds_end]), errMsg(sourcefile, __LINE__))

    arrays_equal = .true.
    do i = bounds_beg, bounds_end
       if (.not. shr_infnan_isnan(bulk(i)) .and. .not. shr_infnan_isnan(tracer(i))) then
          ! neither value is nan: check error tolerance                        
          val1 = bulk(i)
          val2 = tracer(i)
          if (val1 == 0.0_r8 .and. val2 == 0.0_r8) then
             ! trap special case were both are zero to avoid division by zero. values equal                                                                    
          else if (abs(val1 - val2) / max(abs(val1), abs(val2)) > tolerance) then
             arrays_equal = .false.
             diffloc = i
             exit
          else
             ! error < tolerance, considered equal.                             
          end if
       else if (shr_infnan_isnan(bulk(i)) .and. shr_infnan_isnan(tracer(i))) then
          ! both values are nan: values are considered equal.                   
       else
          ! only one value is nan, not equal                                    
          arrays_equal = .false.
          diffloc = i
          exit
       end if
    end do

    if (.not. arrays_equal) then
       write(iulog, '(a, a, a)') 'ERROR in ', subname, ': tracer does not agree with bulk water'
       write(iulog, '(a, a)') 'Called from: ', trim(caller_location)
       if (present(level)) then
          write(iulog, '(a, a, a, i0)') 'Variable: ', trim(name), ', level ', level
       else
          write(iulog, '(a, a)') 'Variable: ', trim(name)
       end if
       write(iulog, '(a, i0)') 'First difference at index: ', diffloc
       write(iulog, '(a, e25.17)') 'Bulk  : ', bulk(diffloc)
       write(iulog, '(a, e25.17)') 'Tracer: ', tracer(diffloc)
       call endrun(msg=subname//': tracer does not agree with bulk water')
    end if
  end subroutine CompareBulkToTracer

end module WaterTracerUtils
