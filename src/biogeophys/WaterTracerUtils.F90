module WaterTracerUtils

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Utility routines for working with water tracers
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type, get_beg, get_end
  use clm_varctl     , only : iulog
  use abortutils     , only : endrun
  use shr_infnan_mod , only : shr_infnan_isnan
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : spval
  use WaterTracerContainerType, only : water_tracer_container_type
  use shr_sys_mod    , only : shr_sys_flush

  implicit none
  save
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: AllocateVar1d
  public :: AllocateVar2d
  public :: CalcTracerFromBulk
  public :: CalcTracerFromBulkFixedRatio
  public :: CompareBulkToTracer

  ! !PRIVATE MEMBER DATA:
  
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine AllocateVar1d(var, name, container, bounds, subgrid_level, ival)
    !
    ! !DESCRIPTION:
    ! Allocate a 1-d water tracer variable and add it to the container
    !
    ! Assumes var is not yet allocated (otherwise there will be a memory leak)
    !
    ! !ARGUMENTS:
    real(r8), pointer                 , intent(out)   :: var(:)
    character(len=*)                  , intent(in)    :: name          ! variable name
    type(water_tracer_container_type) , intent(inout) :: container
    type(bounds_type)                 , intent(in)    :: bounds
    integer                           , intent(in)    :: subgrid_level ! one of the BOUNDS_SUBGRID levels defined in decompMod
    real(r8)                   , intent(in), optional :: ival          ! initial value, if not NaN
    !
    ! !LOCAL VARIABLES:
    integer :: begi, endi

    character(len=*), parameter :: subname = 'AllocateVar1d'
    !-----------------------------------------------------------------------

    begi = get_beg(bounds, subgrid_level)
    endi = get_end(bounds, subgrid_level)
    allocate(var(begi:endi))
    if (present(ival)) then
       var(:) = ival
    else
       var(:) = nan
    end if

    call container%add_var( &
         var = var, &
         begi = begi, &
         description = name, &
         subgrid_level = subgrid_level)

  end subroutine AllocateVar1d

  !-----------------------------------------------------------------------
  subroutine AllocateVar2d(var, name, container, bounds, subgrid_level, dim2beg, dim2end, ival)
    !
    ! !DESCRIPTION:
    ! Allocate a 2-d water tracer variable and add it to the container
    !
    ! A separate variable is added to the container for each level
    !
    ! Assumes var is not yet allocated (otherwise there will be a memory leak)
    !
    ! !ARGUMENTS:
    real(r8), pointer                 , intent(out)   :: var(:,:)
    character(len=*)                  , intent(in)    :: name          ! variable name
    type(water_tracer_container_type) , intent(inout) :: container
    type(bounds_type)                 , intent(in)    :: bounds
    integer                           , intent(in)    :: subgrid_level ! one of the BOUNDS_SUBGRID levels defined in decompMod
    integer                           , intent(in)    :: dim2beg
    integer                           , intent(in)    :: dim2end
    real(r8)                   , intent(in), optional :: ival          ! initial value, if not NaN
    !
    ! !LOCAL VARIABLES:
    integer :: begi, endi
    integer :: lev
    character(len=8) :: lev_str
    character(len=:), allocatable :: description

    character(len=*), parameter :: subname = 'AllocateVar2d'
    !-----------------------------------------------------------------------

    begi = get_beg(bounds, subgrid_level)
    endi = get_end(bounds, subgrid_level)
    allocate(var(begi:endi, dim2beg:dim2end))
    if (present(ival)) then
       var(:,:) = ival
    else
       var(:,:) = nan
    end if

    do lev = dim2beg, dim2end
       write(lev_str, '(i0)') lev
       description = trim(name) // ', level ' // lev_str
       call container%add_var( &
            var = var(:,lev), &
            begi = begi, &
            description = description, &
            subgrid_level = subgrid_level)
    end do

  end subroutine AllocateVar2d

  !-----------------------------------------------------------------------
  subroutine CalcTracerFromBulk(lb, num_pts, filter_pts, &
       bulk_source, bulk_val, tracer_source, tracer_val)
    !
    ! !DESCRIPTION:
    ! Calculate a tracer variable from a corresponding bulk variable when the ratio of
    ! the tracer to the bulk should be based on the ratio in some 'source' variable.
    !
    ! Typically this is used to calculate a tracer flux given a bulk flux. In this case,
    ! the source variable should be the source of the flux - since the tracer ratio of
    ! the flux will be the same as the tracer ratio in the source state variable.
    !
    ! tracer_val will be updated within the given filter, and will remain at its original
    ! values elsewhere
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: lb                 ! lower bound for arrays
    integer  , intent(in)    :: num_pts            ! number of points in the filter
    integer  , intent(in)    :: filter_pts(:)      ! filter in which tracer_val should be updated
    real(r8) , intent(in)    :: bulk_source(lb:)   ! values of the source for this variable, for bulk
    real(r8) , intent(in)    :: bulk_val(lb:)      ! values of the variable of interest, for bulk
    real(r8) , intent(in)    :: tracer_source(lb:) ! values of the source for this variable, for the tracer
    real(r8) , intent(inout) :: tracer_val(lb:)    ! output values of the variable of interest, for the tracer
    !
    ! !LOCAL VARIABLES:
    integer :: num
    integer :: fn, n

    character(len=*), parameter :: subname = 'CalcTracerFromBulk'
    !-----------------------------------------------------------------------

    num = size(bulk_val)
    SHR_ASSERT((size(bulk_source) == num), errMsg(sourcefile, __LINE__))
    SHR_ASSERT((size(tracer_source) == num), errMsg(sourcefile, __LINE__))
    SHR_ASSERT((size(tracer_val) == num), errMsg(sourcefile, __LINE__))

    do fn = 1, num_pts
       n = filter_pts(fn)

       if (bulk_source(n) /= 0._r8) then
          ! Standard case
          tracer_val(n) = bulk_val(n) * (tracer_source(n) / bulk_source(n))
       else if (bulk_val(n) == 0._r8 .and. tracer_source(n) == 0._r8) then
          ! This is acceptable: bulk_source, bulk_val and tracer_source are all 0
          tracer_val(n) = 0._r8
       else if (bulk_val(n) /= 0._r8) then
          write(iulog,*) subname//' ERROR: Non-zero bulk val despite zero bulk source:'
          write(iulog,*) 'bulk_val = ', bulk_val(n)
          write(iulog,*) 'at n = ', n
          write(iulog,*) 'This would lead to an indeterminate tracer val.'
          call endrun(msg=subname//': Non-zero bulk val despite zero bulk source', &
               additional_msg=errMsg(sourcefile, __LINE__))
       else if (tracer_source(n) /= 0._r8) then
          ! NOTE(wjs, 2018-09-28) To avoid this error, we might need code elsewhere that
          ! sets tracer state variables to 0 if the corresponding bulk state variable is 0
          ! and the tracer state is originally near 0 (within roundoff) - in order deal
          ! with roundoff issues arising during state updates. (There's a bit of
          ! discussion of this point in https://github.com/ESCOMP/ctsm/issues/487.)
          write(iulog,*) subname//' ERROR: Non-zero tracer source despite zero bulk source:'
          write(iulog,*) 'tracer_source = ', tracer_source(n)
          write(iulog,*) 'at n = ', n
          write(iulog,*) 'This would lead to an indeterminate tracer val.'
          call endrun(msg=subname//': Non-zero tracer source despite zero bulk source', &
               additional_msg=errMsg(sourcefile, __LINE__))
       else
          write(iulog,*) subname//' ERROR: unhandled condition; we should never get here.'
          write(iulog,*) 'This indicates a programming error in this subroutine.'
          write(iulog,*) 'bulk_val = ', bulk_val(n)
          write(iulog,*) 'bulk_source = ', bulk_source(n)
          write(iulog,*) 'tracer_source = ', tracer_source(n)
          write(iulog,*) 'at n = ', n
          call endrun(msg=subname//': unhandled condition; we should never get here', &
               additional_msg=errMsg(sourcefile, __LINE__))
       end if
    end do

  end subroutine CalcTracerFromBulk


  !-----------------------------------------------------------------------
  subroutine CalcTracerFromBulkFixedRatio(bulk, ratio, tracer)
    !
    ! !DESCRIPTION:
    ! Calculate a tracer variable from a corresponding bulk variable when the tracer
    ! should be a fixed ratio times the bulk
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: bulk(:)
    real(r8), intent(in) :: ratio   ! ratio of tracer to bulk
    real(r8), intent(inout) :: tracer(:)
    !
    ! !LOCAL VARIABLES:
    integer :: num
    integer :: i

    character(len=*), parameter :: subname = 'CalcTracerFromBulkFixedRatio'
    !-----------------------------------------------------------------------

    num = size(bulk)
    SHR_ASSERT((size(tracer) == num), errMsg(sourcefile, __LINE__))
    do i = 1, num
       tracer(i) = bulk(i) * ratio
    end do

  end subroutine CalcTracerFromBulkFixedRatio


  !-----------------------------------------------------------------------
  subroutine CompareBulkToTracer(bounds_beg, bounds_end, &
       bulk, tracer, ratio, caller_location, name)
    !
    ! !DESCRIPTION:
    ! Compare the bulk and tracer quantities; abort if they differ
    !
    ! !ARGUMENTS:
    ! We could get bounds_beg and bounds_end from the lbound and ubound of the bulk or
    ! tracer arrays, but passing them in helps catch any accidental omission of bounds
    ! slicing in the caller (e.g., passing in foo_col rather than
    ! foo_col(bounds%begc:bounds%endc)).
    integer, intent(in) :: bounds_beg
    integer, intent(in) :: bounds_end
    real(r8), intent(in) :: bulk( bounds_beg: )
    real(r8), intent(in) :: tracer( bounds_beg: )
    real(r8), intent(in) :: ratio 
    character(len=*), intent(in) :: caller_location  ! brief description of where this is called from (for error messages)
    character(len=*), intent(in) :: name             ! variable name (for error messages)
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
          if ( val1 /= spval) then
             val1 = val1 * ratio
          end if
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
       write(iulog, '(a, a)') 'Variable: ', trim(name)
       write(iulog, '(a, i0)') 'First difference at index: ', diffloc
       write(iulog, '(a, e25.17)') 'Bulk  : ', bulk(diffloc)
       write(iulog, '(a, e25.17)') 'Tracer: ', tracer(diffloc)
       write(iulog, '(a, e25.17)') 'ratio: ',ratio
       if (.not. shr_infnan_isnan(bulk(diffloc))) then
          write(iulog, '(a, e25.17)') 'Bulk*ratio: ',bulk(diffloc)*ratio
       end if
       call shr_sys_flush(iulog)
       call endrun(msg=subname//': tracer does not agree with bulk water')
    end if
  end subroutine CompareBulkToTracer

end module WaterTracerUtils
