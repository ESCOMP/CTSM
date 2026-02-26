module mkindexmapMod

  !-----------------------------------------------------------------------
  ! Module containing subroutines for making maps of index data.
  !
  ! Includes routines for using an index map as indices into a lookup
  ! table, to essentially paint-by-number some other field.
  !
  ! WJS (2-1-12): There is a lookup_2d subroutine, but not a lookup_1d (or any other
  ! dimensionality). That is simply because I needed lookup_2d, but have not yet needed a
  ! routine of other dimensionalities. In the future, it would probably be helpful to at
  ! least have lookup_1d and lookup_1d_netcdf. If this is done, see my notes under the
  ! lookup_2d_netcdf routine for some thoughts on avoiding duplication.
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod, only : r8 => shr_kind_r8, cs => shr_kind_cs
  use shr_sys_mod , only : shr_sys_abort

  implicit none
  private

  ! public types:
  ! dim_slice_type: stores information about dimensions that we use for slicing a multi-
  ! dimensional variable
  type dim_slice_type
     character(len=CS) :: name  ! name of this dimension
     integer           :: val   ! index to use for the slice
  end type dim_slice_type
  public :: dim_slice_type

  ! public member functions:
  public :: lookup_2d             ! create map based on a 2-d lookup table
  public :: lookup_2d_netcdf      ! wrapper to lookup_2d; first read table from netcdf file
  public :: which_max             ! get index of the maximum value in an array

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  subroutine lookup_2d(index1, index2, lookup_table, fill_val, data, ierr, &
       nodata, valid_entries, invalid_okay)
    !
    ! Creates a data array using a paint-by-number approach according to a lookup table
    !
    ! This routine operates on a 2-d lookup table. There are therefore two index arrays
    ! (index1 and index2); these index arrays are on the same grid as the desired data array
    ! (thus, index1, index2 and data must all have the same length).  Each output point, n, is
    ! then generally determined as:
    !
    ! data(n) = lookup_table(index1(n), index2(n))
    !      
    ! fill_val: value to put in data array where either:
    ! (a) index1 or index2 are equal to nodata (if nodata is given)
    !     Note that this condition does NOT result in ierr being set
    ! (b) valid_entries(index1(n), index2(n)) is false (if valid_entries is given)
    !     Note that this condition also results in ierr being set, unless invalid_okay is
    !     present and .true.
    !     (If valid_entries is not given, it is treated as being .true. everywhere)
    ! (c) index1 or index2 out of range
    !     Note that this condition also results in ierr being set
    ! 
    ! ierr: error return code (if non-0, indicates first error encountered):
    !    0: no error
    !    1: attempt to assign values from the lookup table that are invalid according
    !       to valid_entries (note: this is not considered an error if invalid_okay is
    !       present and .true.)
    !    2: attempt to access an out-of-range index in lookup table
    ! WJS (2-2-12): My main reason for using ierr rather than aborting in case of error
    ! is to facilitate unit testing
    !
    ! input/output variables
    integer , intent(in) :: index1(:)  ! index into dim 1 of lookup_table
    integer , intent(in) :: index2(:)  ! index into dim 2 of lookup_table
    real(r8), intent(in) :: lookup_table(:,:)
    real(r8), intent(in) :: fill_val   ! value to put in data where we don't have a valid value (see above for details)
    real(r8), intent(out):: data(:)    ! output arary
    integer , intent(out):: ierr       ! error return code (0 = no error)

    ! nodata flag in index1 and index2 (see above for details):
    integer, intent(in), optional :: nodata

    ! which entries are considered valid (see above for details):
    logical, intent(in), optional :: valid_entries(:,:)

    ! invalid_okay: if true, then assigning fill_val because valid_entries is false does
    ! NOT raise an error flag (invalid_okay defaults to false, meaning an error is
    ! raised in this case):
    logical, intent(in), optional :: invalid_okay

    ! local variables:
    integer :: n
    integer :: i1, i2
    integer :: data_size          ! size of index1, index2 and data arrays
    integer :: table_n1           ! size of dimension 1 of lookup table
    integer :: table_n2           ! size of dimension 2 of lookup table
    logical :: linvalid_okay      ! local version of invalid_okay
    logical, allocatable :: lvalid_entries(:,:)  ! local version of valid_entries
    character(len=*), parameter :: subname = 'lookup_2d'
    !-----------------------------------------------------------------------

    ierr = 0

    ! Error-check array sizes

    data_size = size(data)
    if (size(index1) /= data_size .or. size(index2) /= data_size) then
       write(6,*) subname//' ERROR: data array sizes do not match'
       write(6,*) 'size(data)   = ', data_size
       write(6,*) 'size(index1) = ', size(index1)
       write(6,*) 'size(index2) = ', size(index2)
       call shr_sys_abort()
    end if

    table_n1 = size(lookup_table,1)
    table_n2 = size(lookup_table,2)
    if (present(valid_entries)) then
       if (size(valid_entries,1) /= table_n1 .or. size(valid_entries,2) /= table_n2) then
          write(6,*) subname//' ERROR: size of valid_entries does not match lookup_table'
          write(6,*) 'size(lookup_table)  = ', table_n1, table_n2
          write(6,*) 'size(valid_entries) = ', size(valid_entries,1), &
               size(valid_entries,2)
          call shr_sys_abort()
       end if
    end if

    ! Set local version of invalid_okay & valid_entries

    if (present(invalid_okay)) then
       linvalid_okay = invalid_okay
    else
       linvalid_okay = .false.
    end if

    allocate(lvalid_entries(table_n1, table_n2))
    if (present(valid_entries)) then
       lvalid_entries(:,:) = valid_entries(:,:)
    else
       lvalid_entries(:,:) = .true.
    end if

    ! Do the lookups

    do n = 1, data_size
       i1 = index1(n)
       i2 = index2(n)

       ! First handle special cases:

       ! index is nodata flag (this is NOT an error)
       if (present(nodata)) then
          if (i1 == nodata .or. i2 == nodata) then
             data(n) = fill_val
             cycle
          end if
       end if

       ! index out of range
       if (i1 <= 0 .or. i1 > table_n1 .or. &
           i2 <= 0 .or. i2 > table_n2) then
          data(n) = fill_val
          if (ierr == 0) ierr = 2
          cycle
       end if

       ! lookup table entry is invalid
       if (.not. lvalid_entries(i1, i2)) then
          data(n) = fill_val
          if (.not. linvalid_okay) then
             if (ierr == 0) ierr = 1
          end if
          cycle
       end if

       ! Finally, the "normal" case, if none of the special cases were triggered:
       data(n) = lookup_table(i1, i2)
    end do

    deallocate(lvalid_entries)

  end subroutine lookup_2d

  !------------------------------------------------------------------------------
  subroutine lookup_2d_netcdf(pioid, tablename, lookup_has_invalid, &
       dimname1, dimname2, n_extra_dims, &
       index1, index2, fill_val, data, ierr, &
       extra_dims, nodata, invalid_okay)
    !
    ! Wrapper to lookup_2d that first reads the lookup table from a netcdf file
    !
    ! If lookup_has_invalid is false, then we treat all lookup table entries as valid data
    ! (i.e., all valid_entries are true in the call to lookup_2d). If lookup_has_invalid is
    ! true, then we read the _FillValue attribute for the lookup table variable, and consider
    ! any table entry with value _FillValue to be an invalid entry, thus putting fill_val in
    ! these data locations (and raising an error flag unless invalid_okay is present and
    ! true).
    !
    ! The dimension given by dimname1 -- with the associated indices given by index1 -- is the
    ! fastest-varying dimension in the lookup table. Dimension dimname2 (associated with
    ! index2) is the second-fastest-varying dimension. Similarly, extra_dims should be ordered
    ! from faster-varying to slowest-varying dimension.  (The first dimension in extra_dims is
    ! the third-fastest-varying dimension in the lookup table.)
    !
    ! n_extra_dims gives the number of extra dimensions (in addition to the first two) in the
    ! lookup table. We take a single 2-d slice of the lookup table, by using a single value of
    ! each of these other dimensions. If n_extra_dims > 0, then extra_dims must be present,
    ! with at least n_extra_dims entries. Each entry in extra_dims gives the name of a
    ! dimension and the dimension index to use for the slice.
    !
    ! If size(extra_dims) > n_extra_dims, then we use the first n_extra_dims entries in
    ! extra_dims. If n_extra_dims = 0, then extra_dims is ignored.
    !
    ! Note that we ignore any coordinate variables associated with the dimensions of the
    ! lookup table; we simply treat the lookup table indices as 1,2,3,...
    !
    ! See the lookup_2d documentation for documentation of some other arguments
    !
    ! WJS (2-1-12): Some thoughts on avoiding duplication if we eventually want similar
    ! routines, lookup_1d_netcdf, lookup_3d_netcdf, etc.:
    !
    ! Much of the code in lookup_2d_netcdf could then be pulled out to a shared subroutine
    ! (e.g., much of the error-checking code).
    !
    ! Or, maybe better: we could try to make a single lookup_netcdf subroutine that handles
    ! 1-d, 2-d and any other dimensionality. To do that, we would (1) make a generic interface
    ! (of which lookup_1d and lookup_2d would be implementations); (2) change the repeated
    ! arguments in lookup_2d_netcdf (*1 and *2) to arrays -- maybe using an array of a derived
    ! type containing these arguments; (3) if possible, initially read the lookup table into a
    ! 1-d array (if the netcdf call allows reading a n-d array into a 1-d array) (if netcdf
    ! doesn't allow this, then I think we could achieve the same thing by reading 1-d slices
    ! of the lookup table in a loop, building the full lookup table as a long 1-d array); (4)
    ! in the call to the generic 'lookup' function, reshape the 1-d lookup table
    ! appropriately. (Note: I think it would be challenging to combine lookup_1d and lookup_2d
    ! (etc.)  into a single routine using a similar method.)
    !
    ! input/output variables
    type(file_desc_t) , intent(inout) :: pioid
    character(len=*)  , intent(in)    :: tablename          ! name of the lookup table variable
    logical           , intent(in)    :: lookup_has_invalid ! should we use _FillValue? (see above)
    character(len=*)  , intent(in)    :: dimname1           ! name of the first (fastest-varying) dimension of the lookup table
    character(len=*)  , intent(in)    :: dimname2           ! name of the second dimension of the lookup table
    integer           , intent(in)    :: n_extra_dims       ! number of extra dimensions in the lookup table

    ! The following arguments are passed directly to lookup_2d:
    integer           , intent(in)    :: index1(:)          ! index into dim 1 of lookup table
    integer           , intent(in)    :: index2(:)          ! index into dim 2 of lookup table
    real(r8)          , intent(in)    :: fill_val           ! value to put in data where we don't have a valid value
    real(r8)          , intent(out)   :: data(:)            ! output array
    integer           , intent(out)   :: ierr               ! error return code from the call to lookup_2d

    ! slice to use if lookup table variable has more than 2 dimensions:
    type(dim_slice_type), intent(in), optional :: extra_dims(:)

    ! nodata flag in index1 and index2, passed directly to lookup_2d:
    integer             , intent(in), optional :: nodata

    ! flag for whether trying to use a lookup table value that is equal to the _FillValue
    ! should raise an error flag
    ! (irrelevant if lookup_has_invalid is .false.)
    ! (passed directly to lookup_2d - see the documentation there for more details)
    logical             , intent(in), optional :: invalid_okay

    ! !LOCAL VARIABLES:
    type(var_desc_t)               :: pio_varid
    integer                        :: ndims              ! total number of dimensions of lookup table
    integer                        :: ndims_expected     ! value we expect for ndims, for error checking
    integer                        :: i                  ! index
    integer                        :: rcode              ! error status
    real(r8)                       :: table_fillval      ! value of the _FillValue attribute for the lookup table
    character(len=CS), allocatable :: dimnames(:)        ! dimension names
    integer , allocatable          :: dimids(:)          ! dimension ids
    integer , allocatable          :: dimlens(:)         ! dimension lengths
    integer , allocatable          :: starts(:)          ! starting indices for reading lookup table
    integer , allocatable          :: counts(:)          ! dimension counts for reading lookup table
    real(r8), allocatable          :: lookup_table(:,:)
    logical , allocatable          :: valid_entries(:,:) ! which entries of the lookup table are considered valid
    character(len=*), parameter :: subname = 'lookup_2d_netcdf'
    !-----------------------------------------------------------------------

    ! Error-check extra_dims
    if (n_extra_dims > 0) then
       if (.not. present(extra_dims)) then
          write(6,*) subname//' ERROR: extra_dims must be present for n_extra_dims > 0'
          call shr_sys_abort()
       end if

       if (size(extra_dims) < n_extra_dims) then
          write(6,*) subname//' ERROR: not enough extra dimensions given'
          write(6,*) 'n_extra_dims =     ', n_extra_dims
          write(6,*) 'size(extra_dims) = ', size(extra_dims)
          call shr_sys_abort()
       end if
    end if

    ! Determine number of expected dimensions in the table, and actual number of
    ! dimensions in the netcdf file

    ndims_expected = 2 + n_extra_dims

    rcode = pio_inq_varid(pioid, trim(tablename), pio_varid)
    rcode = pio_inq_varndims(pioid, pio_varid, ndims)
    if (ndims /= ndims_expected) then
       write(6,*) subname//' ERROR: unexpected number of dimensions in ', &
            trim(tablename)
       write(6,*) 'ndims = ', ndims
       write(6,*) 'expected (based on n_extra_dims): ', ndims_expected
       call shr_sys_abort()
    end if

    ! Get dimension names & sizes, and error-check them
    allocate(dimids(ndims), dimlens(ndims), dimnames(ndims))
    rcode = pio_inq_vardimid (pioid, pio_varid, dimids)
    do i = 1, ndims
       rcode = pio_inq_dimname (pioid, dimids(i), dimnames(i))
       rcode = pio_inq_dimlen (pioid, dimids(i), dimlens(i))
    end do
    call check_dimname(dimnames(1), dimname1, 1)
    call check_dimname(dimnames(2), dimname2, 2)
    do i = 1, n_extra_dims
       call check_dimname(dimnames(2+i), extra_dims(i)%name, 2+i)
       call check_dimsize(dimlens(2+i), extra_dims(i)%val, 2+i)
    end do

    ! Read the lookup table; if the given variable has more than 2 dimensions, we read
    ! a single 2-d slice

    allocate(starts(ndims), counts(ndims))
    allocate(lookup_table(dimlens(1), dimlens(2)))
    starts(1:2) = 1
    counts(1:2) = dimlens(1:2)
    do i = 1, n_extra_dims
       starts(2+i) = extra_dims(i)%val
       counts(2+i) = 1
    end do
    rcode = pio_get_var(pioid, pio_varid, starts, counts, lookup_table)

    !allocate(lookup_table(dimlens(1), dimlens(2)))
    !rcode = pio_get_var(pioid, pio_varid, lookup_table)

    ! Determine which entries are valid
    allocate(valid_entries(size(lookup_table, 1), size(lookup_table, 2)))
    valid_entries(:,:) = .true.
    if (lookup_has_invalid) then
       rcode = pio_get_att(pioid, pio_varid, '_FillValue', table_fillval)
       where (lookup_table == table_fillval)
          valid_entries = .false.
       end where
    end if

    ! Do the lookups
    call lookup_2d(index1, index2, lookup_table, fill_val, data, ierr, nodata=nodata, &
         valid_entries=valid_entries, invalid_okay=invalid_okay)

    deallocate(valid_entries)
    deallocate(lookup_table)
    deallocate(starts, counts)
    deallocate(dimids, dimlens, dimnames)

  contains

    !------------------------------------------------------------------------------
    subroutine check_dimname(actual, expected, i)
      ! Make sure names are equal; if not, stop with an error message

      character(len=*), intent(in) :: actual, expected
      integer         , intent(in) :: i  ! dimension number, for output purposes

      if (actual /= expected) then
         write(6,*) subname//' ERROR: unexpected dimension name in ', trim(tablename)
         write(6,*) 'dimension #', i
         write(6,*) 'actual:   ', trim(actual)
         write(6,*) 'expected: ', trim(expected)
         call shr_sys_abort()
      end if
    end subroutine check_dimname

    !------------------------------------------------------------------------------
    subroutine check_dimsize(length, index, i)
      ! Make sure dimension length is long enough; if not, stop with an error message

      integer, intent(in) :: length, index
      integer, intent(in) :: i  ! dimension number, for output purposes

      if (index > length) then
         write(6,*) subname//' ERROR: desired index exceeds dimension length in ', &
              trim(tablename)
         write(6,*) 'dimension #', i
         write(6,*) 'index:  ', index
         write(6,*) 'length: ', length
         call shr_sys_abort()
      end if
    end subroutine check_dimsize

  end subroutine lookup_2d_netcdf

  !------------------------------------------------------------------------------
  subroutine which_max(arr, maxval, maxindex, lbound)
    !
    ! Returns maximum value in arr along with the index of the maximum value
    ! If multiple values are tied, returns index of the first maximum
    !
    ! input/output variables
    real(r8), intent(in) :: arr(:)
    real(r8), intent(out):: maxval   ! maximum value in arr(:)
    integer , intent(out):: maxindex ! first index of maxval
    integer , intent(in), optional :: lbound ! lower bound of indices of arr; 
                                             ! if not supplied, assumed to be 1:
    ! local variables
    integer :: i
    !-----------------------------------------------------------------------

    maxindex = 1
    maxval = arr(1)
    do i = 2, size(arr)
       if (arr(i) > maxval) then
          maxindex = i
          maxval = arr(i)
       end if
    end do
    if (present(lbound)) then
       maxindex = maxindex + (lbound - 1)
    end if
  end subroutine which_max

end module mkindexmapMod
