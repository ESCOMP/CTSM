module ncdio_utils

  !-----------------------------------------------------------------------
  ! This module provides higher-level netcdf i/o utilities, which build on ncdio_pio.
  !
  ! The main reason for putting these utilities in a separate module (rather than putting
  ! them in ncdio_pio) is to enhance testability: These routines can be unit tested with
  ! a stub version of ncdio_pio.
  use ncdio_pio
  !
  implicit none
  save
  private

  public :: find_var_on_file  ! given a list of possible variables, find the one that exists on the file

contains

  !-----------------------------------------------------------------------
  subroutine find_var_on_file(ncid, varname_list, is_dim, varname_on_file)
    !
    ! !DESCRIPTION:
    ! Given a colon-delimited list of possible variable names, return the first one that
    ! was found on the file.
    !
    ! If none are found, arbitrarily return the first variable in the list. (Doing this
    ! rather than returning a special flag simplifies the logic elsewhere - allowing the
    ! ncd_io call to fail rather than requiring extra error-checking logic.)
    !
    ! !USES:
    use shr_string_mod, only : shr_string_listGetNum, shr_string_listGetName
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(inout) :: ncid            ! netcdf file id
    character(len=*)  , intent(in)    :: varname_list    ! colon-delimited list of possible variable names
    logical           , intent(in)    :: is_dim          ! if .true., then look at dimensions rather than variables
    character(len=*)  , intent(out)   :: varname_on_file ! first variable from the list that was found on file
    !
    ! !LOCAL VARIABLES:
    integer :: num_vars
    integer :: n
    logical :: found
    logical :: readvar
    character(len=len(varname_on_file)) :: cur_varname

    character(len=*), parameter :: subname = 'find_var_on_file'
    !-----------------------------------------------------------------------

    num_vars = shr_string_listGetNum(varname_list)

    found = .false.
    n = 1
    do while ((.not. found) .and. (n <= num_vars))
       call shr_string_listGetName(varname_list, n, cur_varname)
       call check_var_or_dim(ncid, cur_varname, is_dim=is_dim, exists=found)
       n = n + 1
    end do

    if (found) then
       varname_on_file = cur_varname
    else
       ! If none are found, arbitrarily return the first variable in the list
       call shr_string_listGetName(varname_list, 1, varname_on_file)
    end if

  end subroutine find_var_on_file

end module ncdio_utils
