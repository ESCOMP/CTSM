module mkchecksMod

  !-----------------------------------------------------------------------
  ! Generic routines to check validity of output fields
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4

  implicit none
  private

  public :: min_bad        ! check the minimum value of a field
  public :: max_bad        ! check the maximum value of a field

  interface min_bad
     module procedure min_bad_int
     module procedure min_bad_r8
     module procedure min_bad_r4
  end interface min_bad

  interface max_bad
     module procedure max_bad_int
     module procedure max_bad_r8
     module procedure max_bad_r4
  end interface max_bad

!===============================================================
contains
!===============================================================

  logical function min_bad_r8(data, min_allowed, varname)

    ! Confirm that no value of data is less than min_allowed.
    ! Returns true if errors found, false otherwise.
    ! Also prints offending points

    ! !ARGUMENTS:
    real(r8)          , intent(in) :: data(:)      ! array of data to check
    real(r8)          , intent(in) :: min_allowed  ! minimum valid value
    character(len=*)  , intent(in) :: varname      ! name of field

    ! !LOCAL VARIABLES:
    logical :: errors_found        ! true if any errors have been found
    integer :: n                   ! index
    character(len=*), parameter :: subname = 'min_bad_r8'
    !------------------------------------------------------------------------------

    errors_found = .false.

    do n = 1, size(data)
       if (data(n) < min_allowed) then
          write(6,*) subname//' ERROR: ', trim(varname), ' = ', data(n), ' less than ',&
               min_allowed, ' at ', n
          errors_found = .true.
       end if
    end do

    min_bad_r8 = errors_found
  end function min_bad_r8

  !===============================================================
  logical function min_bad_r4(data, min_allowed, varname)

    ! Confirm that no value of data is less than min_allowed.
    ! Returns true if errors found, false otherwise.
    ! Also prints offending points

    ! !ARGUMENTS:
    real(r4)          , intent(in) :: data(:)      ! array of data to check
    real(r4)          , intent(in) :: min_allowed  ! minimum valid value
    character(len=*)  , intent(in) :: varname      ! name of field

    ! !LOCAL VARIABLES:
    logical :: errors_found        ! true if any errors have been found
    integer :: n                   ! index
    character(len=*), parameter :: subname = 'min_bad_r8'
    !------------------------------------------------------------------------------

    errors_found = .false.

    do n = 1, size(data)
       if (data(n) < min_allowed) then
          write(6,*) subname//' ERROR: ', trim(varname), ' = ', data(n), ' less than ',&
               min_allowed, ' at ', n
          errors_found = .true.
       end if
    end do

    min_bad_r4 = errors_found
  end function min_bad_r4

  !===============================================================
  logical function min_bad_int(data, min_allowed, varname)

    ! !DESCRIPTION:
    ! Confirm that no value of data is less than min_allowed.
    ! Returns true if errors found, false otherwise.
    ! Also prints offending points

    ! !ARGUMENTS:
    integer          , intent(in) :: data(:)      ! array of data to check
    integer          , intent(in) :: min_allowed  ! minimum valid value
    character(len=*) , intent(in) :: varname      ! name of field

    ! !LOCAL VARIABLES:
    logical :: errors_found        ! true if any errors have been found
    integer :: n                   ! index
    character(len=*), parameter :: subname = 'min_bad_int'
    !------------------------------------------------------------------------------

    errors_found = .false.

    do n = 1, size(data)
       if (data(n) < min_allowed) then
          write(6,*) subname//' ERROR: ', trim(varname), ' = ', data(n), ' less than ',&
               min_allowed, ' at ', n
          errors_found = .true.
       end if
    end do

    min_bad_int = errors_found
  end function min_bad_int

  !===============================================================
  logical function max_bad_r8(data, max_allowed, varname)

    ! Confirm that no value of data is greate than max_allowed.
    ! Returns true if errors found, false otherwise.
    ! Also prints offending points

    ! !ARGUMENTS:
    real(r8)          , intent(in) :: data(:)      ! array of data to check
    real(r8)          , intent(in) :: max_allowed  ! maximum valid value
    character(len=*)  , intent(in) :: varname      ! name of field

    ! !LOCAL VARIABLES:
    logical :: errors_found        ! true if any errors have been found
    integer :: n                   ! index
    character(len=*), parameter :: subname = 'max_bad_r8'
    !------------------------------------------------------------------------------

    errors_found = .false.

    do n = 1, size(data)
       if (data(n) > max_allowed) then
          write(6,*) subname//' ERROR: ', trim(varname), ' = ', data(n), ' greater than ',&
               max_allowed, ' at ', n
          errors_found = .true.
       end if
    end do

    max_bad_r8 = errors_found
  end function max_bad_r8

  !===============================================================
  logical function max_bad_r4(data, max_allowed, varname)

    ! Confirm that no value of data is greate than max_allowed.
    ! Returns true if errors found, false otherwise.
    ! Also prints offending points

    ! !ARGUMENTS:
    real(r4)          , intent(in) :: data(:)      ! array of data to check
    real(r4)          , intent(in) :: max_allowed  ! maximum valid value
    character(len=*)  , intent(in) :: varname      ! name of field

    ! !LOCAL VARIABLES:
    logical :: errors_found        ! true if any errors have been found
    integer :: n                   ! index
    character(len=*), parameter :: subname = 'max_bad_r8'
    !------------------------------------------------------------------------------

    errors_found = .false.

    do n = 1, size(data)
       if (data(n) > max_allowed) then
          write(6,*) subname//' ERROR: ', trim(varname), ' = ', data(n), ' greater than ',&
               max_allowed, ' at ', n
          errors_found = .true.
       end if
    end do

    max_bad_r4 = errors_found
  end function max_bad_r4

  !===============================================================
  logical function max_bad_int(data, max_allowed, varname)

    ! !DESCRIPTION:
    ! Confirm that no value of data is greate than max_allowed.
    ! Returns true if errors found, false otherwise.
    ! Also prints offending points

    ! !ARGUMENTS:
    integer          , intent(in) :: data(:)      ! array of data to check
    integer          , intent(in) :: max_allowed  ! maximum valid value
    character(len=*) , intent(in) :: varname      ! name of field

    ! !LOCAL VARIABLES:
    logical :: errors_found        ! true if any errors have been found
    integer :: n                   ! index
    character(len=*), parameter :: subname = 'max_bad_int'
    !------------------------------------------------------------------------------

    errors_found = .false.

    do n = 1, size(data)
       if (data(n) > max_allowed) then
          write(6,*) subname//' ERROR: ', trim(varname), ' = ', data(n), ' greater than ',&
               max_allowed, ' at ', n
          errors_found = .true.
       end if
    end do

    max_bad_int = errors_found
  end function max_bad_int

end module mkchecksMod
