module IssueFixedMetadataHandler

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handles the writing and reading of netCDF metadata saying whether a given issue has
  ! been fixed as of the writing of this file.
  !
  ! !USES:
  use clm_varctl , only : iulog
  use abortutils , only : endrun
  use ncdio_pio  , only : file_desc_t, ncd_global
  use ncdio_pio  , only : ncd_putatt, ncd_getatt, check_att
  use ncdio_pio  , only : ncd_inqnatts, ncd_inqattname
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: write_issue_fixed_metadata ! write metadata to a netCDF file indicating that a given issue has been fixed
  public :: read_issue_fixed_metadata  ! read metadata from a netCDF file indicating whether a given issue has been fixed
  public :: copy_issue_fixed_metadata  ! copy all issue_fixed metadata from one netCDF file to another
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: issue_fixed_attname  ! get the attribute name indicating whether a given issue has been fixed

  !
  ! !PRIVATE DATA:
  character(len=*), parameter, private :: issue_fixed_att_prefix = 'issue_fixed_'  ! all attributes written by this module start with this prefix

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine write_issue_fixed_metadata(ncid, writing_finidat_interp_dest_file, &
       issue_num, attribute_value)
    !
    ! !DESCRIPTION:
    ! Write a global attribute to the given netcdf file, indicating that the given issue has been fixed
    !
    ! The netcdf file must be open and in define mode. If this is a finidat_interp_dest
    ! file, nothing is written: Whether this issue has been fixed will be determined by
    ! the presence of the attribute on the input file to init_interp.
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(inout) :: ncid                             ! netcdf id; must be open in define mode
    logical           , intent(in)    :: writing_finidat_interp_dest_file ! true if this is a finidat_interp_dest file; if so, we don't write anything to it
    integer           , intent(in)    :: issue_num                        ! number of the issue for which we're writing metadata

    ! Value to write; if not specified, write 1
    !
    ! Note that a 0 value is used in read_issue_fixed_metadata to indicate that the
    ! attribute wasn't found on the file at all. So you should generally avoid setting
    ! attribute_value = 0 unless you want the behavior to be the same as if the attribute
    ! wasn't on the file at all.
    integer , intent(in), optional   :: attribute_value
    !
    ! !LOCAL VARIABLES:
    integer :: l_attribute_value  ! local version of attribute_value
    character(len=:), allocatable :: attname

    character(len=*), parameter :: subname = 'write_issue_fixed_metadata'
    !-----------------------------------------------------------------------

    if (present(attribute_value)) then
       l_attribute_value = attribute_value
    else
       l_attribute_value = 1
    end if

    attname = issue_fixed_attname(issue_num)

    ! In principle, the caller could do this check. But we do it here because it can be
    ! very important to avoid writing the metadata to the finidat_interp_dest file, and
    ! we don't want to rely on the caller remembering to do this check.
    if (.not. writing_finidat_interp_dest_file) then
       call ncd_putatt(ncid, ncd_global, attname, l_attribute_value)
    end if

  end subroutine write_issue_fixed_metadata

  !-----------------------------------------------------------------------
  subroutine read_issue_fixed_metadata(ncid, issue_num, attribute_value)
    !
    ! !DESCRIPTION:
    ! Read metadata from a netCDF file indicating whether a given issue has been fixed
    !
    ! If the given attribute isn't found, returns 0; if it is found, returns the attribute's value
    !
    ! Typically, then, the caller can check whether the issue is fixed on the given file
    ! by checking whether the attribute's value is 0: if it is 0, the issue has NOT been
    ! fixed; if it is non-zero, the issue has been fixed. (But a particular attribute may
    ! have more nuanced meanings as well, indicated by different integer values.)
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid ! netcdf id
    integer          , intent(in)    :: issue_num ! number of the issue for which we're reading metadata
    integer          , intent(out)   :: attribute_value ! value of the given issue fixed attribute; 0 if not present on file
    !
    ! !LOCAL VARIABLES:
    character(len=:), allocatable :: attname
    logical :: att_found

    character(len=*), parameter :: subname = 'read_issue_fixed_metadata'
    !-----------------------------------------------------------------------

    attname = issue_fixed_attname(issue_num)
    call check_att(ncid, ncd_global, attname, att_found)
    if (att_found) then
       call ncd_getatt(ncid, ncd_global, attname, attribute_value)
    else
       attribute_value = 0
    end if

  end subroutine read_issue_fixed_metadata

  !-----------------------------------------------------------------------
  subroutine copy_issue_fixed_metadata(ncidi, ncido)
    !
    ! !DESCRIPTION:
    ! Copy all issue_fixed metadata from one netCDF file to another
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncidi ! input file (source)
    type(file_desc_t), intent(inout) :: ncido ! output file (destination); must be open in define mode
    !
    ! !LOCAL VARIABLES:
    integer :: input_natts  ! number of attributes on the input file
    integer :: i
    character(len=256) :: attname
    integer :: attribute_value_input
    integer :: attribute_value_output
    logical :: already_exists_on_output ! whether the given attribute already exists on the output file

    character(len=*), parameter :: subname = 'copy_issue_fixed_metadata'
    !-----------------------------------------------------------------------

    call ncd_inqnatts(ncidi, input_natts)
    do i = 1, input_natts
       call ncd_inqattname(ncidi, ncd_global, i, attname)
       if (attname(1:len_trim(issue_fixed_att_prefix)) == trim(issue_fixed_att_prefix)) then
          call ncd_getatt(ncidi, ncd_global, attname, attribute_value_input)
          call check_att(ncido, ncd_global, attname, already_exists_on_output)
          if (already_exists_on_output) then
             ! Avoid trying to overwrite an existing attribute on the output file. If the
             ! existing value is the same as the value in the input file, there's nothing
             ! to do; if the existing value differed, it's unclear what we should do, so
             ! abort.
             call ncd_getatt(ncido, ncd_global, attname, attribute_value_output)
             if (attribute_value_output /= attribute_value_input) then
                write(iulog,*) subname//' ERROR: Attribute already exists on output file with a different value'
                write(iulog,*) 'Attribute, input value, output value = ', &
                     attname, attribute_value_input, attribute_value_output
                call endrun(subname//' ERROR: Attribute already exists on output file with a different value')
             end if

          else
             ! If the attribute doesn't already exist on the output, create it there
             call ncd_putatt(ncido, ncd_global, attname, attribute_value_input)
          end if
       end if
    end do

  end subroutine copy_issue_fixed_metadata

  !-----------------------------------------------------------------------
  function issue_fixed_attname(issue_num) result(attname)
    !
    ! !DESCRIPTION:
    ! Get the attribute name indicating whether a given issue has been fixed
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: attname  ! function result
    integer, intent(in) :: issue_num
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: issue_num_str

    character(len=*), parameter :: subname = 'issue_fixed_attname'
    !-----------------------------------------------------------------------

    write(issue_num_str, '(i0)') issue_num
    attname = issue_fixed_att_prefix // trim(issue_num_str)

  end function issue_fixed_attname

end module IssueFixedMetadataHandler
