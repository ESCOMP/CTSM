module initInterpMultilevelInterp

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines a class for handling multi-level fields by interpolating based on
  ! a coordinate variable. This interpolation scheme also allows different levels to
  ! belong to different level_classes; the interpolation for a level belonging to
  ! level_class N only considers source levels that also belong to level_class N.
  !
  !
  ! ------------------------------------------------------------------------
  ! NOTE(wjs, 2015-10-22) Note about handling of spval (and NaN) values in the source and
  ! destination data (i.e., the data themselves, as opposed to the coordinate or
  ! level_class arrays): It feels like there are two different possible uses for spval in
  ! multi-level arrays; this class assumes the "static" use of spval (#1 below):
  !
  ! (1) A "static" spval, e.g., based on depth for a given variable: e.g., "given this
  ! vertical discretization, this variable only applies to a depth of 5 m, and so should
  ! always remain spval below 5 m". For this case, if the destination is spval, then it
  ! should remain spval - because spval is more about the grid (vertical discretization,
  ! surface dataset, etc.) than about anything dynamic in time. A spval in the source
  ! data should be ignored (rather than being used in the interpolation), so that the
  ! destination array gets filled for all levels where it was not originally spval, using
  ! extrapolation if necessary (e.g., if the source array had spvals at some levels where
  ! the destination array needs valid data).
  !
  ! (2) A "dynamic" spval, e.g., signaling that this layer should be treated specially
  ! right now, but that may change in time. For this case, an spval in the destination
  ! should be treated no differently from any other value: it should be replaced by
  ! values from the source grid. An spval in the source would need to be handled
  ! specially (e.g.: if both sides of the interpolation are spval then the result is
  ! spval; if one side is spval, then use the other point). (We can't simply ignore /
  ! exclude source levels with spval in this case, because an spval has true meaning that
  ! should be transferred to the nearby levels in the destination array.)
  !
  ! These two uses are fundamentally incompatible: Either we keep spvals in the
  ! destination array as spval (the "static" case), or we replace spvals in the
  ! destination array (the "dynamic" case). The only way we could handle both
  ! possibilities is if variables had metadata saying whether they were using a "dynamic"
  ! or "static" treatment of spvals - but that feels confusing.
  !
  ! Looking through current multi-level variables, it looked like there were a number of
  ! uses of "static" spvals, and one use of "dynamic" spvals (in the ch4 code). I have now
  ! removed the "dynamic" spval.
  !
  ! This class assumes the "static" use of spval - i.e., case (1) above. This is partly to
  ! reflect the current usage in the code, and partly because it's easier to detect
  ! improper use of a "dynamic" spval when we assume "static" spval (through an LII test)
  ! than it is to detect improper use of a "static" spval when we assume "dynamic" spval.
  ! (The latter would tend to just show up as a problem when interpolating from one
  ! vertical grid to another, and it's hard to confirm correctness of that interpolation
  ! in an automated test.)
  !
  ! Thus, multi-level arrays (at least those that are potentially interpolated vertically)
  ! should NOT use the "dynamic" use of spval (case (2)). That is, spval should only be
  ! used in places where that spval will remain spval throughout the run.
  ! ------------------------------------------------------------------------
  !
  ! !USES:
#include "shr_assert.h" 

  use shr_kind_mod             , only : r8 => shr_kind_r8
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use abortutils               , only : endrun
  use clm_varctl               , only : iulog
  use clm_varcon               , only : spval, ispval
  use initInterpMultilevelBase , only : interp_multilevel_type
  use array_utils              , only : pack_wrapper

  implicit none
  private
  save

  ! Public types

  public :: interp_multilevel_interp_type

  type, extends(interp_multilevel_type) :: interp_multilevel_interp_type
     private
     character(len=:), allocatable :: coord_varname     ! name of coordinate variable (just used for identification purposes)

     ! The following 'source' arrays store information about the source grid, but
     ! regridded to the destination grid. e.g., coordinates_source(n,i) gives the
     ! coordinate value for level n on the source grid, for the source grid point that
     ! maps to destination point i.
     real(r8), allocatable :: coordinates_source(:,:)   ! coordinate values on source grid [lev, pt]
     real(r8), allocatable :: coordinates_dest(:,:)     ! coordinate values on dest grid [lev, pt]
     real(r8), allocatable :: dzsoi_source(:,:)         ! soil thicknesses on source grid [lev, pt]
     real(r8), allocatable :: dzsoi_dest(:,:)           ! soil thicknesses on dest grid [lev, pt]
     integer , allocatable :: level_classes_source(:,:) ! class indices on source grid [lev, pt]
     integer , allocatable :: level_classes_dest(:,:)   ! class indices on dest grid [lev, pt]

     integer :: npts         ! number of spatial points on dest grid
     integer :: nlev_source  ! number of levels on source grid
     integer :: nlev_dest    ! number of levels on dest grid
   contains
     ! Public methods from base class
     procedure :: check_npts
     procedure :: interp_multilevel
     procedure :: get_description

     ! Public methods specific to this class
     procedure :: get_nlev_source  ! get number of levels on source grid
     procedure :: get_nlev_dest    ! get number of levels on dest grid

     ! Private methods
     procedure, private :: check_coordinate_array
     procedure, private, nopass :: confirm_monotonically_increasing
     procedure, private, nopass :: interp_onelevel   ! do the interpolation for a single dest level
     procedure, private, nopass :: is_missing        ! true if a value is considered missing
  end type interp_multilevel_interp_type

  interface interp_multilevel_interp_type
     module procedure constructor_with_levclasses
     module procedure constructor_no_levclasses
  end interface interp_multilevel_interp_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor_no_levclasses(coordinates_source, coordinates_dest, &
                                     dzsoi_source, dzsoi_dest, &
                                     coord_varname) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Construct an interp_multilevel_interp_type object, where all levels have the same
    ! level_class.
    !
    ! Coordinates must be monotonically increasing.
    !
    ! coordinates_source gives information about the source grid, but regridded to the
    ! destination grid. So coordinates_source(n,i) gives the coordinate value for level
    ! n on the source grid, for the source grid point that maps to destination point i.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(interp_multilevel_interp_type) :: this  ! function result
    real(r8), intent(in) :: coordinates_source(:,:) ! coordinate values on source grid [lev, pt]
    real(r8), intent(in) :: coordinates_dest(:,:)   ! coordinate values on dest grid [lev, pt]
    real(r8), intent(in) :: dzsoi_source(:,:)       ! dzsoi values on source grid [lev, pt]
    real(r8), intent(in) :: dzsoi_dest(:,:)         ! dzsoi values on dest grid [lev, pt]
    character(len=*), intent(in) :: coord_varname   ! name of coordinate variable (just used for identification purposes)
    !
    ! !LOCAL VARIABLES:
    integer :: npts
    integer :: nlev_source, nlev_dest
    integer, allocatable :: level_classes_source(:,:)
    integer, allocatable :: level_classes_dest(:,:)

    character(len=*), parameter :: subname = 'constructor_no_levclasses'
    !-----------------------------------------------------------------------

    npts = size(coordinates_source, 2)
    SHR_ASSERT_FL((size(coordinates_dest, 2) == npts), sourcefile, __LINE__)
    SHR_ASSERT_FL((size(dzsoi_source, 2) == npts), sourcefile, __LINE__)
    SHR_ASSERT_FL((size(dzsoi_dest, 2) == npts), sourcefile, __LINE__)

    nlev_source = size(coordinates_source, 1)
    nlev_dest = size(coordinates_dest, 1)

    allocate(level_classes_source(nlev_source, npts), source = 1)
    allocate(level_classes_dest(nlev_dest, npts), source = 1)

    this = interp_multilevel_interp_type( &
         coordinates_source = coordinates_source, &
         coordinates_dest = coordinates_dest, &
         dzsoi_source = dzsoi_source, &
         dzsoi_dest = dzsoi_dest, &
         level_classes_source = level_classes_source, &
         level_classes_dest = level_classes_dest, &
         coord_varname = coord_varname)

  end function constructor_no_levclasses


  !-----------------------------------------------------------------------
  function constructor_with_levclasses(coordinates_source, coordinates_dest, &
       dzsoi_source, dzsoi_dest, &
       level_classes_source, level_classes_dest, &
       coord_varname) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Construct an interp_multilevel_interp_type object, where classes are specified for
    ! each level.
    !
    ! Coordinates must be monotonically increasing.
    !
    ! coordinates_source and level_classes_source give information about the source grid,
    ! but regridded to the destination grid. e.g., coordinates_source(n,i) gives the
    ! coordinate value for level n on the source grid, for the source grid point that maps
    ! to destination point i.
    !
    ! For the 'level_classes': The particular values are not important; the important
    ! thing is that, for a given column, levels that are fundamentally different should
    ! have different values. This ensures that data are not copied from one class of
    ! levels to another (e.g., between soil and bedrock).
    !
    ! However, a level whose class is ispval is treated specially: This is treated as a
    ! non-existent level: In the destination data, its value is left unchanged; in the
    ! source data, its value is never used in the interpolation.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(interp_multilevel_interp_type) :: this  ! function result
    real(r8), intent(in) :: coordinates_source(:,:)   ! coordinate values on source grid [lev, pt]
    real(r8), intent(in) :: coordinates_dest(:,:)     ! coordinate values on dest grid [lev, pt]
    real(r8), intent(in) :: dzsoi_source(:,:)         ! dzsoi values on source grid [lev, pt]
    real(r8), intent(in) :: dzsoi_dest(:,:)           ! dzsoi values on dest grid [lev, pt]
    integer , intent(in) :: level_classes_source(:,:) ! class indices on source grid [lev, pt]
    integer , intent(in) :: level_classes_dest(:,:)   ! class indices on dest grid [lev, pt]
    character(len=*), intent(in) :: coord_varname     ! name of coordinate variable (just used for identification purposes)
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'constructor_with_levclasses'
    !-----------------------------------------------------------------------

    this%coord_varname = trim(coord_varname)

    this%npts = size(coordinates_source, 2)
    SHR_ASSERT_FL((size(coordinates_dest, 2) == this%npts), sourcefile, __LINE__)

    this%nlev_source = size(coordinates_source, 1)
    this%nlev_dest = size(coordinates_dest, 1)

    SHR_ASSERT_ALL_FL((shape(level_classes_source) == [this%nlev_source, this%npts]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((shape(level_classes_dest) == [this%nlev_dest, this%npts]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((shape(dzsoi_source) == [this%nlev_source, this%npts]), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((shape(dzsoi_dest) == [this%nlev_dest, this%npts]), sourcefile, __LINE__)

    do i = 1, this%npts
       call this%check_coordinate_array(coordinates_source(:,i), level_classes_source(:,i))
    end do
    do i = 1, this%npts
       call this%check_coordinate_array(coordinates_dest(:,i), level_classes_dest(:,i))
    end do

    this%coordinates_source = coordinates_source
    this%coordinates_dest = coordinates_dest
    this%dzsoi_source = dzsoi_source
    this%dzsoi_dest = dzsoi_dest
    this%level_classes_source = level_classes_source
    this%level_classes_dest = level_classes_dest

  end function constructor_with_levclasses


  ! ========================================================================
  ! Public methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine check_npts(this, npts, varname)
    !
    ! !DESCRIPTION:
    ! Checks the number of destination points, to ensure that this interpolator is
    ! appropriate for this variable. This should be called once for each variable.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_interp_type), intent(in) :: this
    integer, intent(in) :: npts             ! number of dest points (on this processor)
    character(len=*), intent(in) :: varname ! variable name (for diagnostic output)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'check_npts'
    !-----------------------------------------------------------------------

    if (npts /= this%npts) then
       write(iulog,*) subname//' ERROR: mismatch in number of dest points for ', &
            trim(varname)
       write(iulog,*) 'Number of dest points: ', npts
       write(iulog,*) 'Expected number of dest points: ', this%npts
       call endrun(msg=subname//' ERROR: mismatch in number of points for '//&
            trim(varname) // ' ' // errMsg(sourcefile, __LINE__))
    end if

  end subroutine check_npts

  !-----------------------------------------------------------------------
  subroutine interp_multilevel(this, data_dest, data_source, index_dest, scale_by_thickness)
    !
    ! !DESCRIPTION:
    ! Interpolates a multi-level field from source to dest, for a single point.
    !
    ! This version does a true interpolation, using a coordinate variable. The coordinate
    ! variable (along with the group to which each level belongs) can vary for each
    ! spatial point. Thus, index_dest is used in this version, and is must be within the
    ! bounds of the metadata. This index should be 1-based.
    !
    ! If level_classes were provided for this object (i.e., level_classes_source and
    ! level_classes_dest), then: For a given destination level, the interpolation is done
    ! only over source levels whose class matches the destination level's class. In
    ! addition, if the destination level's class is ispval, then the destination data is
    ! left unchanged for that level.
    !
    ! Levels whose data value is spval are treated the same as levels whose class is
    ! ispval: If the destination data begins as spval, then it is left as spval (with the
    ! assumption that this is truly meant to remain as missing data in the destination).
    ! And any source level that has data = spval is ignored. NaN values are treated the
    ! same as spval in this respect.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_interp_type), intent(in) :: this
    real(r8) , intent(inout) :: data_dest(:)
    real(r8) , intent(in)    :: data_source(:)
    integer  , intent(in)    :: index_dest
    logical  , intent(in)    :: scale_by_thickness
    !
    ! !LOCAL VARIABLES:
    integer :: lev_dest
    integer :: level_class_dest
    integer :: lev_source

    ! source information for this index_dest
    real(r8) :: my_level_classes_source(this%nlev_source)
    real(r8) :: my_coordinates_source(this%nlev_source)
    real(r8) :: my_dzsoi_source(this%nlev_source)

    ! whether each source level is in the destination level_class
    logical :: source_levels_in_class(this%nlev_source)

    ! data and coordinates packed to just contain levels in the destination level_class:
    real(r8), allocatable :: data_source_in_class(:)
    real(r8), allocatable :: coordinates_source_in_class(:)
    real(r8), allocatable :: dzsoi_source_in_class(:)

    character(len=*), parameter :: subname = 'interp_multilevel'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((size(data_dest) == this%nlev_dest), sourcefile, __LINE__)
    SHR_ASSERT_FL((size(data_source) == this%nlev_source), sourcefile, __LINE__)
    SHR_ASSERT_FL((index_dest >= 1 .and. index_dest <= this%npts), sourcefile, __LINE__)

    my_level_classes_source(:) = this%level_classes_source(:, index_dest)
    my_coordinates_source(:)   = this%coordinates_source(:, index_dest)
    my_dzsoi_source(:)         = this%dzsoi_source(:, index_dest)

    do lev_dest = 1, this%nlev_dest
       level_class_dest = this%level_classes_dest(lev_dest, index_dest)
       if (level_class_dest /= ispval) then
          if (.not. is_missing(data_dest(lev_dest))) then
             do lev_source = 1, this%nlev_source
                if (my_level_classes_source(lev_source) /= level_class_dest) then
                   source_levels_in_class(lev_source) = .false.
                else if (is_missing(data_source(lev_source))) then
                   source_levels_in_class(lev_source) = .false.
                else
                   source_levels_in_class(lev_source) = .true.
                end if
             end do
             call pack_wrapper(data_source_in_class, data_source, source_levels_in_class)
             call pack_wrapper(coordinates_source_in_class, my_coordinates_source, source_levels_in_class)
             call pack_wrapper(dzsoi_source_in_class, my_dzsoi_source, source_levels_in_class)
             call this%interp_onelevel( &
                  data_dest = data_dest(lev_dest), &
                  dzsoi_dest = this%dzsoi_dest(lev_dest, index_dest), &
                  coordinate_dest = this%coordinates_dest(lev_dest, index_dest), &
                  data_source = data_source_in_class, &
                  dzsoi_source = dzsoi_source_in_class, &
                  coordinates_source = coordinates_source_in_class, &
                  scale_by_thickness = scale_by_thickness)
          end if
       end if
    end do

  end subroutine interp_multilevel

  !-----------------------------------------------------------------------
  pure function get_description(this) result(description)
    !
    ! !DESCRIPTION:
    ! Returns a text description of this interpolator
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: description  ! function result
    class(interp_multilevel_interp_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_description'
    !-----------------------------------------------------------------------
    
    description = 'Interpolate using '//this%coord_varname

  end function get_description

  !-----------------------------------------------------------------------
  pure integer function get_nlev_source(this)
    ! Get number of levels on source grid
    class(interp_multilevel_interp_type), intent(in) :: this

    get_nlev_source = this%nlev_source
  end function get_nlev_source

  !-----------------------------------------------------------------------
  pure integer function get_nlev_dest(this)
    ! Get number of levels on dest grid
    class(interp_multilevel_interp_type), intent(in) :: this

    get_nlev_dest = this%nlev_dest
  end function get_nlev_dest


  ! ========================================================================
  ! Private methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine check_coordinate_array(this, coordinates, level_classes)
    !
    ! !DESCRIPTION:
    ! Check validity of a coordinate array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_interp_type), intent(in) :: this
    real(r8), intent(in) :: coordinates(:)   ! coordinates at all levels, for one point
    integer , intent(in) :: level_classes(:) ! classes corresponding to coordinates
    !
    ! !LOCAL VARIABLES:
    integer :: nlevs

    ! subset of coordinates just in levels that exist (with existence defined based on
    ! level_classes)
    real(r8), allocatable :: coordinates_in_existing_levels(:)

    character(len=*), parameter :: subname = 'check_coordinate_array'
    !-----------------------------------------------------------------------

    nlevs = size(coordinates)
    SHR_ASSERT_FL((size(level_classes) == nlevs), sourcefile, __LINE__)

    call pack_wrapper(coordinates_in_existing_levels, coordinates, level_classes /= ispval)

    if (any(is_missing(coordinates_in_existing_levels))) then
       ! In principle we could handle spvals in coordinate arrays just like we handle
       ! ispval in the level_class or spval in data. However, this likely indicates an
       ! error, so for now we look out for that and abort if it is found.
       call endrun(msg='spvals or NaNs found in coordinate array where level_class /= ispval; ' // &
            'this is currently unhandled ' // errMsg(sourcefile, __LINE__))
    end if

    call this%confirm_monotonically_increasing(coordinates_in_existing_levels)

  end subroutine check_coordinate_array


  !-----------------------------------------------------------------------
  subroutine confirm_monotonically_increasing(data)
    !
    ! !DESCRIPTION:
    ! Confirms that an array is monotonically increasing. Dies if not.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: data(:)
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'confirm_monotonically_increasing'
    !-----------------------------------------------------------------------

    do i = 2, size(data)
       if (data(i-1) >= data(i)) then
          write(iulog,*) subname//' ERROR: array not monotonically increasing: '
          write(iulog,*) data(i-1), data(i)
          call endrun(msg=subname//" ERROR: array not monotonically increasing"// &
               errMsg(sourcefile, __LINE__))
       end if
    end do

  end subroutine confirm_monotonically_increasing

  !-----------------------------------------------------------------------
  subroutine interp_onelevel(data_dest, dzsoi_dest, coordinate_dest, &
                             data_source, dzsoi_source, coordinates_source, &
                             scale_by_thickness)
    !
    ! !DESCRIPTION:
    ! Do the interpolation for a single destination level
    !
    ! !ARGUMENTS:
    real(r8), intent(inout) :: data_dest
    real(r8), intent(in)    :: dzsoi_dest
    real(r8), intent(in)    :: coordinate_dest
    real(r8), intent(in)    :: data_source(:)
    real(r8), intent(in)    :: dzsoi_source(:)
    real(r8), intent(in)    :: coordinates_source(:)
    logical , intent(in)    :: scale_by_thickness
    !
    ! !LOCAL VARIABLES:
    logical :: copylevel
    logical :: found
    integer :: nlev_source
    integer :: lev
    integer :: index_lower
    real(r8) :: wt_lower  ! Weight factor for interpolation of extensive properties
    real(r8) :: wt_lower_plus_1  ! Weight factor for interpolation of extensive properties

    real(r8), parameter :: eps = 1.e-13_r8
    character(len=*), parameter :: subname = 'interp_onelevel'
    !-----------------------------------------------------------------------

    nlev_source = size(coordinates_source)

    if (nlev_source == 0) then
       ! If there is no source information, then leave the destination data at its
       ! original value
       return
    end if

    ! ------------------------------------------------------------------------
    ! Find level(s) to use for interpolation
    ! ------------------------------------------------------------------------

    found = .false.

    if (coordinate_dest < coordinates_source(1)) then
       found = .true.
       copylevel = .true.
       index_lower = 1

    else if (coordinate_dest > coordinates_source(nlev_source)) then
       found = .true.
       copylevel = .true.
       index_lower = nlev_source

    else
       ! See if coordinate_dest matches one of the source coordinates (within roundoff)
       do lev = 1, nlev_source
          if ((abs(coordinate_dest - coordinates_source(lev)) < eps)) then
             found = .true.
             copylevel = .true.
             index_lower = lev
             exit
          end if
       end do

       if (.not. found) then
          ! Find the interval in which coordinate_dest falls
          do lev = 1, nlev_source
             if ( (coordinate_dest  >  coordinates_source(lev)) &
                  .and. (coordinate_dest  <  coordinates_source(lev+1)) ) then
                found = .true.
                copylevel = .false.
                index_lower = lev
                exit
             end if
          end do
       end if

    end if

    if (.not. found) then
       call endrun(subname//' ERROR: Could not find levels to use for interpolation' // &
            errMsg(sourcefile, __LINE__))
    end if

    ! Assemble the weights for scaling certain variables by soil thicknesses.
    ! Such variables represent so-called "extensive" properties and have units
    ! such as kg/m2 or mm or mm/s. The weights by which we scale these variables
    ! are unit converters in reality. The denominators dzsoi_source(index_lower)
    ! and dzsoi_source(index_lower+1) convert data_source from a surface density
    ! to a volume density. The numerator dzsoi_dest converts data_dest back to
    ! a surface density.
    if (scale_by_thickness) then
       ! Set wt_lower
       if (dzsoi_source(index_lower) <= 0._r8) then
          call endrun(msg=subname//' ERROR1: about to divide by zero or negative dzsoi '// &
               errMsg(sourcefile, __LINE__))
       end if
       wt_lower = dzsoi_dest / dzsoi_source(index_lower)

       ! Set wt_lower_plus_1 only if (.not. copylevel)
       if (.not. copylevel) then ! not using wt_lower_plus_1
          if (dzsoi_source(index_lower+1) <= 0._r8) then
             call endrun(msg=subname//' ERROR2: about to divide by zero or negative dzsoi '// &
                  errMsg(sourcefile, __LINE__))
          end if
          wt_lower_plus_1 = dzsoi_dest / dzsoi_source(index_lower+1)
       end if
    else  ! the no scaling option for both copylevel and .not. copylevel
       wt_lower = 1.0_r8
       wt_lower_plus_1 = 1.0_r8
    end if

    ! ------------------------------------------------------------------------
    ! Do the interpolation
    ! ------------------------------------------------------------------------

    if ( copylevel) then
       data_dest = data_source(index_lower) * wt_lower
    else
       data_dest = &
            data_source(index_lower+1) * wt_lower_plus_1 &
            * (coordinate_dest - coordinates_source(index_lower)) &
            / (coordinates_source(index_lower+1) - coordinates_source(index_lower)) + &
            data_source(index_lower) * wt_lower &
            * (coordinates_source(index_lower+1) - coordinate_dest ) &
            / (coordinates_source(index_lower+1) - coordinates_source(index_lower))
    end if

  end subroutine interp_onelevel

  !-----------------------------------------------------------------------
  elemental logical function is_missing(val)
    !
    ! !DESCRIPTION:
    ! Returns true if the given value is considered missing.
    !
    ! spval is treated as missing. NaNs are treated the same as spval. Considering the
    ! destination data, the assumption here is that, if a value was uninitialized before
    ! (so is NaN), then it shouldn't be set in init_interp.
    !
    ! !USES:
    use shr_infnan_mod , only : isnan => shr_infnan_isnan
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: val
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'is_missing'
    !-----------------------------------------------------------------------

    if (isnan(val)) then
       is_missing = .true.
    else if (val == spval) then
       is_missing = .true.
    else
       is_missing = .false.
    end if

  end function is_missing


end module initInterpMultilevelInterp
