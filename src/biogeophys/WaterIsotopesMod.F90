module WaterIsotopesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! !REVISION HISTORY
  !
  !   Created David Noone <dcn@colorado.edu> - Thu Jun 16 18:45:47 MDT 2005
  !   Updated David Noone <dcn@colorado.edu> - Fri Jul 13 15:40:29 MDT 2007
  !   Updated Tony Wong <anthony.e.wong@colorado.edu - Wed Aug 14 10:23:38 MDT 2013
  !       (for CLM4 use)
  !-----------------------------------------------------------------------

  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varctl   , only: iulog  
  use abortutils   , only: endrun
  use shr_log_mod        , only : errMsg => shr_log_errMsg

  implicit none
  private
  save

  !-----------------------------------------------------------------------
  ! Module variables
  !
  !
  ! Public interfaces
  public :: WisoCompareBulkToTracer

  
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !
  ! Module variables
  !

contains

  !=======================================================================
  function WisoCompareBulkToTracer(bounds_beg, bounds_end, &
       bulk, tracer, name) result(wiso_inconsistency)

    use abortutils , only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg

    implicit none

    integer, intent(in) :: bounds_beg
    integer, intent(in) :: bounds_end
    real(r8), intent(in) :: bulk(bounds_beg:bounds_end)
    real(r8), intent(in) :: tracer(bounds_beg:bounds_end)
    character(len=*), intent(in) :: name

    logical :: arrays_equal
    logical :: wiso_inconsistency, inconsistency_is_error
    logical, parameter :: debug = .false.

    inconsistency_is_error = .true.

    wiso_inconsistency = .false.


    if (debug) then
       write(iulog, '(a, a, a, i3, a)') "WISO: comparing ", trim(name), " tracer with bulk water...."
    end if
    arrays_equal = compare_arrays_loop(bounds_beg, bounds_end, &
         bulk(bounds_beg:bounds_end), tracer(bounds_beg:bounds_end))
    if (.not. arrays_equal) then
       write(iulog, '(a, a, a, i3, a)') "WISO: ", trim(name), " tracer does not agree with bulk water."
       wiso_inconsistency = .true.
    else
       if (debug) then
          write(iulog, '(a, a, a, i3, a)') "WISO: ", trim(name), " tracer agrees with bulk water."
       end if
    end if

    if (wiso_inconsistency .and. inconsistency_is_error) then
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
  end function WisoCompareBulkToTracer

  function compare_arrays_loop(bounds_beg, bounds_end, &
       array1, array2) result(arrays_equal)

    use shr_infnan_mod, only : shr_infnan_isnan

    implicit none

    integer, intent(in) :: bounds_beg
    integer, intent(in) :: bounds_end
    real(r8), intent(in) :: array1(bounds_beg:bounds_end), array2(bounds_beg:bounds_end)

    real(r8), parameter :: tolerance = 1.0e-7_r8

    logical :: arrays_equal
    integer :: i
    real(r8) :: val1, val2

    arrays_equal = .true.
    do i = bounds_beg, bounds_end
       if (.not. shr_infnan_isnan(array1(i)) .and. .not. shr_infnan_isnan(array2(i))) then
          ! neither value is nan: check error tolerance                        
          val1 = array1(i)
          val2 = array2(i)
          if (val1 == 0.0_r8 .and. val2 == 0.0_r8) then
             ! trap special case were both are zero to avoid division by zero. values equal                                                                    
          else if (abs(val1 - val2) / max(abs(val1), abs(val2)) > tolerance) then
             arrays_equal = .false.
             exit
          else
             ! error < tolerance, considered equal.                             
          end if
       else if (shr_infnan_isnan(array1(i)) .and. shr_infnan_isnan(array2(i))) then
          ! both values are nan: values are considered equal.                   
       else
          ! only one value is nan, not equal                                    
          arrays_equal = .false.
          exit
       end if
    end do

  end function compare_arrays_loop

  
end module WaterIsotopesMod
