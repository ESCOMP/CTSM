module multilevel_interp_factory
  ! Factory module for creating instances of interp_multilevel_type

  use initInterpMultilevelInterp, only : interp_multilevel_interp_type
  use shr_kind_mod             , only : r8 => shr_kind_r8

  implicit none
  private
  save

  public :: create_multilevel_interp_no_levclasses
  public :: create_multilevel_interp_with_levclasses

contains

  ! ========================================================================
  ! Public routines
  ! ========================================================================

  function create_multilevel_interp_no_levclasses(coordinates_source, coordinates_dest, &
       index_dest, npts_dest) &
       result(interpolator)
    ! Arguments:
    type(interp_multilevel_interp_type) :: interpolator  ! function result
    real(r8), intent(in) :: coordinates_source(:)  ! coordinates in source data for index_dest
    real(r8), intent(in) :: coordinates_dest(:)    ! coordinates in dest data for index_dest

    integer, intent(in) :: index_dest    ! dest index of interest in tests
    integer, intent(in) :: npts_dest     ! total number of points wanted in dest

    ! Local variables:
    real(r8), allocatable :: coordinates_source_all(:,:)
    real(r8), allocatable :: coordinates_dest_all(:,:)

    character(len=*), parameter :: coord_varname = 'COORD'
    !-----------------------------------------------------------------------

    call create_coordinate_arrays( &
         coordinates_source = coordinates_source, &
         coordinates_dest = coordinates_dest, &
         coordinates_source_all = coordinates_source_all, &
         coordinates_dest_all = coordinates_dest_all, &
         index_dest = index_dest, &
         npts_dest = npts_dest)

    interpolator = interp_multilevel_interp_type( &
         coordinates_source = coordinates_source_all, &
         coordinates_dest   = coordinates_dest_all, &
         coord_varname      = coord_varname)

  end function create_multilevel_interp_no_levclasses

  function create_multilevel_interp_with_levclasses(coordinates_source, coordinates_dest, &
       level_classes_source, level_classes_dest, &
       index_dest, npts_dest) &
       result(interpolator)
    ! Arguments:
    type(interp_multilevel_interp_type) :: interpolator  ! function result
    real(r8), intent(in) :: coordinates_source(:)   ! coordinates in source data for index_dest
    real(r8), intent(in) :: coordinates_dest(:)     ! coordinates in dest data for index_dest
    integer , intent(in) :: level_classes_source(:) ! class indices in source data for index_dest
    integer , intent(in) :: level_classes_dest(:)   ! class indices in dest data for index_dest

    integer, intent(in) :: index_dest    ! dest index of interest in tests
    integer, intent(in) :: npts_dest     ! total number of points wanted in dest

    ! Local variables:
    real(r8), allocatable :: coordinates_source_all(:,:)
    real(r8), allocatable :: coordinates_dest_all(:,:)
    integer , allocatable :: level_classes_source_all(:,:)
    integer , allocatable :: level_classes_dest_all(:,:)

    character(len=*), parameter :: coord_varname = 'COORD'
    !-----------------------------------------------------------------------

    call create_coordinate_arrays( &
         coordinates_source = coordinates_source, &
         coordinates_dest = coordinates_dest, &
         coordinates_source_all = coordinates_source_all, &
         coordinates_dest_all = coordinates_dest_all, &
         index_dest = index_dest, &
         npts_dest = npts_dest)

    call create_class_arrays( &
         level_classes_source = level_classes_source, &
         level_classes_dest = level_classes_dest, &
         level_classes_source_all = level_classes_source_all, &
         level_classes_dest_all = level_classes_dest_all, &
         index_dest = index_dest, &
         npts_dest = npts_dest)

    interpolator = interp_multilevel_interp_type( &
         coordinates_source = coordinates_source_all, &
         coordinates_dest   = coordinates_dest_all, &
         level_classes_source = level_classes_source_all, &
         level_classes_dest   = level_classes_dest_all, &
         coord_varname      = coord_varname)

  end function create_multilevel_interp_with_levclasses

  ! ========================================================================
  ! Private routines
  ! ========================================================================

  subroutine create_coordinate_arrays(coordinates_source, coordinates_dest, &
       coordinates_source_all, coordinates_dest_all, &
       index_dest, npts_dest)
    ! Arguments:
    real(r8), intent(in) :: coordinates_source(:)  ! coordinates in source data for index_dest
    real(r8), intent(in) :: coordinates_dest(:)    ! coordinates in dest data for index_dest
    real(r8), allocatable, intent(out) :: coordinates_source_all(:,:)
    real(r8), allocatable, intent(out) :: coordinates_dest_all(:,:)

    integer, intent(in) :: index_dest    ! dest index of interest in tests
    integer, intent(in) :: npts_dest     ! total number of points wanted in dest

    ! Local variables:
    integer :: nlevels_source
    integer :: nlevels_dest

    integer :: point, level
    !-----------------------------------------------------------------------

    nlevels_source = size(coordinates_source)
    nlevels_dest = size(coordinates_dest)

    allocate(coordinates_source_all(nlevels_source, npts_dest))
    allocate(coordinates_dest_all(nlevels_dest, npts_dest))

    ! Fill coordinates with garbage
    do point = 1, npts_dest
       do level = 1, nlevels_source
          coordinates_source_all(level, point) = 1000._r8 * level
       end do
    end do

    do point = 1, npts_dest
       do level = 1, nlevels_dest
          coordinates_dest_all(level, point) = 100000._r8 * level
       end do
    end do

    ! But put the passed-in coordinates in index_dest:
    coordinates_source_all(:, index_dest) = coordinates_source(:)
    coordinates_dest_all(:, index_dest) = coordinates_dest(:)

  end subroutine create_coordinate_arrays

  subroutine create_class_arrays(level_classes_source, level_classes_dest, &
       level_classes_source_all, level_classes_dest_all, &
       index_dest, npts_dest)
    ! Arguments:
    integer, intent(in) :: level_classes_source(:)  ! classes in source data for index_dest
    integer, intent(in) :: level_classes_dest(:)    ! classes in dest data for index_dest
    integer, allocatable, intent(out) :: level_classes_source_all(:,:)
    integer, allocatable, intent(out) :: level_classes_dest_all(:,:)

    integer, intent(in) :: index_dest    ! dest index of interest in tests
    integer, intent(in) :: npts_dest     ! total number of points wanted in dest

    ! Local variables:
    integer :: nlevels_source
    integer :: nlevels_dest

    integer, parameter :: default_class = 1
    !-----------------------------------------------------------------------

    nlevels_source = size(level_classes_source)
    nlevels_dest = size(level_classes_dest)

    allocate(level_classes_source_all(nlevels_source, npts_dest), source = default_class)
    allocate(level_classes_dest_all(nlevels_dest, npts_dest), source = default_class)

    level_classes_source_all(:, index_dest) = level_classes_source(:)
    level_classes_dest_all(:, index_dest) = level_classes_dest(:)

  end subroutine create_class_arrays

end module multilevel_interp_factory
