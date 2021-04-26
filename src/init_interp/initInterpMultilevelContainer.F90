module initInterpMultilevelContainer

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines a class that contains one instance of each interp_multilevel
  ! type. This class is responsible for:
  !
  ! (1) Constructing each interp_multilevel object
  !
  ! (2) Determining which interp_multilevel object should be used for each multi-level
  !     field (based on the field's level dimension)
  !
  ! !USES:
#include "shr_assert.h" 
  use shr_kind_mod               , only : r8 => shr_kind_r8
  use initInterpBounds           , only : interp_bounds_type
  use initInterpMultilevelBase   , only : interp_multilevel_type
  use initInterpMultilevelCopy   , only : interp_multilevel_copy_type
  use initInterpMultilevelInterp , only : interp_multilevel_interp_type
  use initInterpMultilevelSnow   , only : interp_multilevel_snow_type
  use initInterpMultilevelSplit  , only : interp_multilevel_split_type, create_interp_multilevel_split_type
  use initInterp2dvar            , only : interp_2dvar_type
  use initInterp1dData           , only : interp_1d_data
  use ncdio_pio                  , only : file_desc_t, var_desc_t, check_var, ncd_io, ncd_inqdlen, ncd_inqdid
  use clm_varctl                 , only : iulog
  use abortutils                 , only : endrun
  use shr_log_mod                , only : errMsg => shr_log_errMsg
  use spmdMod                    , only : masterproc
  use array_utils                , only : transpose_wrapper

  implicit none
  private
  save

  ! Public types

  public :: interp_multilevel_container_type

  type :: interp_multilevel_container_type
     private
     ! Components need to be pointers so that we can return pointers to them.
     !
     ! (Components of a derived type cannot have the target attribute, but rather take on
     ! the target attribute from their parent object. So the alternative to making these
     ! pointers would be to require all instances of this derived type to have the target
     ! attribute.)
     type(interp_multilevel_copy_type), pointer   :: interp_multilevel_copy
     type(interp_multilevel_interp_type), pointer :: interp_multilevel_levgrnd_col
     type(interp_multilevel_interp_type), pointer :: interp_multilevel_levmaxurbgrnd_col
     type(interp_multilevel_interp_type), pointer :: interp_multilevel_levgrnd_pft
     type(interp_multilevel_snow_type), pointer   :: interp_multilevel_levsno
     type(interp_multilevel_snow_type), pointer   :: interp_multilevel_levsno1
     type(interp_multilevel_split_type), pointer  :: interp_multilevel_levtot_col
   contains
     procedure :: init
     procedure :: find_interpolator
     final :: destroy_interp_multilevel_container_type
  end type interp_multilevel_container_type

  ! Private routines

  private :: create_levgrnd_col_interpolators
  private :: create_levgrnd_pft_interpolator
  private :: get_levmaxurbgrnd_metadata
  private :: interp_levgrnd_check_source_file
  private :: create_snow_interpolators

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Public methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine init(this, ncid_source, ncid_dest, bounds_source, bounds_dest, &
       pftindex, colindex)
    !
    ! !DESCRIPTION:
    ! Initialize this interp_multilevel_container_type instance
    !
    ! !ARGUMENTS:
    class(interp_multilevel_container_type), intent(inout) :: this
    type(file_desc_t), target, intent(inout) :: ncid_source ! netcdf ID for source file
    type(file_desc_t), target, intent(inout) :: ncid_dest   ! netcdf ID for dest file
    type(interp_bounds_type), intent(in) :: bounds_source
    type(interp_bounds_type), intent(in) :: bounds_dest

    ! The following give mappings from source to dest for pft and col-level variables.
    ! e.g., colindex(i) gives source col corresponding to dest col i.
    integer, intent(in) :: pftindex(:)
    integer, intent(in) :: colindex(:)
    !
    ! !LOCAL VARIABLES:
    integer :: levgrnd_source ! size of the levgrnd dimension in source
    integer :: levgrnd_dest   ! size of the levgrnd dimension in dest
    integer :: dimid          ! dimension id

    character(len=*), parameter :: subname = 'init'
    !-----------------------------------------------------------------------

    allocate(this%interp_multilevel_copy)
    this%interp_multilevel_copy = interp_multilevel_copy_type()

    call ncd_inqdlen(ncid_source, dimid, levgrnd_source, name='levgrnd')
    call ncd_inqdlen(ncid_dest  , dimid, levgrnd_dest  , name='levgrnd')

    ! Note that there are two (often identical) interpolators for column-level levgrnd:
    ! interp_multilevel_levgrnd_col is used for interpolating variables that are
    ! dimensioned by levgrnd; interp_multilevel_levmaxurbgrnd_col is used as part of the
    ! levtot interpolator for interpolating variables that are dimensioned by levtot,
    ! because levtot is snow plus levmaxurbgrnd. If nlevgrnd >= nlevurb (which is often
    ! the case), then these two are identical; however, if nlevgrnd < nlevurb for source
    ! and/or destination, then interp_multilevel_levmaxurbgrnd_col will have additional
    ! levels beyond those in interp_multilevel_levgrnd_col.
    allocate(this%interp_multilevel_levgrnd_col)
    allocate(this%interp_multilevel_levmaxurbgrnd_col)
    call create_levgrnd_col_interpolators( &
         ncid_source = ncid_source, &
         ncid_dest = ncid_dest, &
         bounds_source = bounds_source, &
         bounds_dest = bounds_dest, &
         levgrnd_source = levgrnd_source, &
         levgrnd_dest = levgrnd_dest, &
         colindex = colindex, &
         interp_multilevel_levgrnd_col = this%interp_multilevel_levgrnd_col, &
         interp_multilevel_levmaxurbgrnd_col = this%interp_multilevel_levmaxurbgrnd_col)

    allocate(this%interp_multilevel_levgrnd_pft)
    call create_levgrnd_pft_interpolator( &
         ncid_source = ncid_source, &
         ncid_dest = ncid_dest, &
         bounds_source = bounds_source, &
         bounds_dest = bounds_dest, &
         levgrnd_source = levgrnd_source, &
         levgrnd_dest = levgrnd_dest, &
         pftindex = pftindex, &
         interp_multilevel_levgrnd_pft = this%interp_multilevel_levgrnd_pft)

    allocate(this%interp_multilevel_levsno)
    allocate(this%interp_multilevel_levsno1)
    call create_snow_interpolators( &
         interp_multilevel_levsno = this%interp_multilevel_levsno, &
         interp_multilevel_levsno1 = this%interp_multilevel_levsno1, &
         ncid_source = ncid_source, &
         bounds_source = bounds_source, &
         bounds_dest = bounds_dest, &
         colindex = colindex)

    ! levtot is two sets of levels: first snow, then levmaxurbgrnd (where levmaxurbgrnd =
    ! max(levgrnd, levurb))
    allocate(this%interp_multilevel_levtot_col)
    this%interp_multilevel_levtot_col = create_interp_multilevel_split_type( &
         interpolator_first_levels = this%find_interpolator('levsno', 'column'), &
         interpolator_second_levels = this%interp_multilevel_levmaxurbgrnd_col, &
         num_second_levels_source = this%interp_multilevel_levmaxurbgrnd_col%get_nlev_source(), &
         num_second_levels_dest = this%interp_multilevel_levmaxurbgrnd_col%get_nlev_dest())

  end subroutine init

  !-----------------------------------------------------------------------
  function find_interpolator(this, lev_dimname, vec_dimname) result(interpolator)
    !
    ! !DESCRIPTION:
    ! Given the name of the level dimension and the vector dimension, return a pointer to
    ! an interpolator that is appropriate for this multi-level variable.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_type), pointer :: interpolator  ! function result

    class(interp_multilevel_container_type), intent(in) :: this
    character(len=*), intent(in) :: lev_dimname  ! name of level dimension
    character(len=*), intent(in) :: vec_dimname  ! name of vector dimension
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'find_interpolator'
    !-----------------------------------------------------------------------

    select case (lev_dimname)
    case ('levgrnd')
       select case (vec_dimname)
       case ('column')
          interpolator => this%interp_multilevel_levgrnd_col
       case ('pft')
          interpolator => this%interp_multilevel_levgrnd_pft
       case default
          call error_not_found(subname, lev_dimname, vec_dimname)
       end select
    case ('levmaxurbgrnd')
       select case (vec_dimname)
       case ('column')
          ! NOTE(wjs, 2020-10-23) Currently, no variables use this interpolator, but we
          ! need it as part of the levtot interpolator, so we might as well support it as
          ! a standalone interpolator in case it's needed in the future.
          interpolator => this%interp_multilevel_levmaxurbgrnd_col
       case default
          call error_not_found(subname, lev_dimname, vec_dimname)
       end select
    case ('levtot')
       select case (vec_dimname)
       case ('column')
          interpolator => this%interp_multilevel_levtot_col
       case default
          call error_not_found(subname, lev_dimname, vec_dimname)
       end select
    case ('levsno')
       interpolator => this%interp_multilevel_levsno
    case ('levsno1')
       interpolator => this%interp_multilevel_levsno1
    case default
       interpolator => this%interp_multilevel_copy
    end select

  contains
    subroutine error_not_found(subname, lev_dimname, vec_dimname)
      ! Write an error message and abort
      character(len=*), intent(in) :: subname
      character(len=*), intent(in) :: lev_dimname
      character(len=*), intent(in) :: vec_dimname

      write(iulog,*) subname//' ERROR: no multi-level interpolator found for:'
      write(iulog,*) 'lev_dimname = ', trim(lev_dimname)
      write(iulog,*) 'vec_dimname = ', trim(vec_dimname)
      call endrun(msg='ERROR: no multi-level interpolator found '//errMsg(sourcefile, __LINE__))
    end subroutine error_not_found

  end function find_interpolator

  ! ========================================================================
  ! Private methods and routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine create_levgrnd_col_interpolators(ncid_source, ncid_dest, &
       bounds_source, bounds_dest, colindex, levgrnd_source, levgrnd_dest, &
       interp_multilevel_levgrnd_col, interp_multilevel_levmaxurbgrnd_col)
    !
    ! !DESCRIPTION:
    ! Create the interp_multilevel_levgrnd_col and interp_multilevel_levmaxurbgrnd_col interpolators
    !
    ! !ARGUMENTS:
    type(file_desc_t), target, intent(inout) :: ncid_source
    type(file_desc_t), target, intent(inout) :: ncid_dest
    type(interp_bounds_type), intent(in) :: bounds_source
    type(interp_bounds_type), intent(in) :: bounds_dest
    integer, intent(in) :: colindex(:)  ! mappings from column-level source to dest points
    integer, intent(in) :: levgrnd_source ! size of the levgrnd dimension in source
    integer, intent(in) :: levgrnd_dest   ! size of the levgrnd dimension in dest
    type(interp_multilevel_interp_type), intent(out) :: interp_multilevel_levgrnd_col
    type(interp_multilevel_interp_type), intent(out) :: interp_multilevel_levmaxurbgrnd_col
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: coordinates_source(:,:)   ! [lev, vec]
    real(r8), allocatable :: coordinates_dest(:,:)     ! [lev, vec]
    real(r8), allocatable :: dzsoi_source(:,:)         ! [lev, vec]
    real(r8), allocatable :: dzsoi_dest(:,:)           ! [lev, vec]
    integer , allocatable :: level_classes_source(:,:) ! [lev, vec]
    integer , allocatable :: level_classes_dest(:,:)   ! [lev, vec]

    character(len=*), parameter :: subname = 'create_levgrnd_col_interpolators'
    !-----------------------------------------------------------------------

    call get_levmaxurbgrnd_metadata( &
         ncid_source = ncid_source, &
         ncid_dest = ncid_dest, &
         bounds_source = bounds_source, &
         bounds_dest = bounds_dest, &
         coord_varname = 'COL_Z', &
         dzsoi_varname = 'DZSOI', &
         level_class_varname = 'LEVGRND_CLASS', &
         sgridindex = colindex, &
         coordinates_source = coordinates_source, &
         coordinates_dest = coordinates_dest, &
         dzsoi_source = dzsoi_source, &
         dzsoi_dest = dzsoi_dest, &
         level_classes_source = level_classes_source, &
         level_classes_dest = level_classes_dest)

    interp_multilevel_levmaxurbgrnd_col = interp_multilevel_interp_type( &
         coordinates_source = coordinates_source, &
         coordinates_dest = coordinates_dest, &
         dzsoi_source = dzsoi_source, &
         dzsoi_dest = dzsoi_dest, &
         level_classes_source = level_classes_source, &
         level_classes_dest = level_classes_dest, &
         coord_varname = 'COL_Z')

    interp_multilevel_levgrnd_col = interp_multilevel_interp_type( &
         coordinates_source = coordinates_source(1:levgrnd_source, :), &
         coordinates_dest = coordinates_dest(1:levgrnd_dest, :), &
         dzsoi_source = dzsoi_source(1:levgrnd_source, :), &
         dzsoi_dest = dzsoi_dest(1:levgrnd_dest, :), &
         level_classes_source = level_classes_source(1:levgrnd_source, :), &
         level_classes_dest = level_classes_dest(1:levgrnd_dest, :), &
         coord_varname = 'COL_Z down to levgrnd')

  end subroutine create_levgrnd_col_interpolators

  !-----------------------------------------------------------------------
  subroutine create_levgrnd_pft_interpolator(ncid_source, ncid_dest, &
       bounds_source, bounds_dest, pftindex, levgrnd_source, levgrnd_dest, &
       interp_multilevel_levgrnd_pft)
    !
    ! !DESCRIPTION:
    ! Create the interp_multilevel_levgrnd_pft interpolator
    !
    ! !ARGUMENTS:
    type(file_desc_t), target, intent(inout) :: ncid_source
    type(file_desc_t), target, intent(inout) :: ncid_dest
    type(interp_bounds_type), intent(in) :: bounds_source
    type(interp_bounds_type), intent(in) :: bounds_dest
    integer, intent(in) :: pftindex(:)  ! mappings from patch-level source to dest points
    integer, intent(in) :: levgrnd_source ! size of the levgrnd dimension in source
    integer, intent(in) :: levgrnd_dest   ! size of the levgrnd dimension in dest
    type(interp_multilevel_interp_type), intent(out) :: interp_multilevel_levgrnd_pft
    !
    ! !LOCAL VARIABLES:
    real(r8), allocatable :: coordinates_source(:,:)   ! [lev, vec]
    real(r8), allocatable :: coordinates_dest(:,:)     ! [lev, vec]
    integer , allocatable :: level_classes_source(:,:) ! [lev, vec]
    integer , allocatable :: level_classes_dest(:,:)   ! [lev, vec]
    real(r8), allocatable :: dzsoi_source(:,:)         ! [lev, vec]
    real(r8), allocatable :: dzsoi_dest(:,:)           ! [lev, vec]

    character(len=*), parameter :: subname = 'create_levgrnd_pft_interpolator'
    !-----------------------------------------------------------------------

    call get_levmaxurbgrnd_metadata( &
         ncid_source = ncid_source, &
         ncid_dest = ncid_dest, &
         bounds_source = bounds_source, &
         bounds_dest = bounds_dest, &
         coord_varname = 'COL_Z_p', &
         dzsoi_varname = 'DZSOI_p', &
         level_class_varname = 'LEVGRND_CLASS_p', &
         sgridindex = pftindex, &
         coordinates_source = coordinates_source, &
         coordinates_dest = coordinates_dest, &
         dzsoi_source = dzsoi_source, &
         dzsoi_dest = dzsoi_dest, &
         level_classes_source = level_classes_source, &
         level_classes_dest = level_classes_dest)

    interp_multilevel_levgrnd_pft = interp_multilevel_interp_type( &
         coordinates_source = coordinates_source(1:levgrnd_source, :), &
         coordinates_dest = coordinates_dest(1:levgrnd_dest, :), &
         dzsoi_source = dzsoi_source(1:levgrnd_source, :), &
         dzsoi_dest = dzsoi_dest(1:levgrnd_dest, :), &
         level_classes_source = level_classes_source(1:levgrnd_source, :), &
         level_classes_dest = level_classes_dest(1:levgrnd_dest, :), &
         coord_varname = 'COL_Z_p down to levgrnd')

  end subroutine create_levgrnd_pft_interpolator

  !-----------------------------------------------------------------------
  subroutine get_levmaxurbgrnd_metadata(ncid_source, ncid_dest, &
       bounds_source, bounds_dest, &
       coord_varname, dzsoi_varname, level_class_varname, sgridindex, &
       coordinates_source, coordinates_dest, dzsoi_source, dzsoi_dest, &
       level_classes_source, level_classes_dest)
    !
    ! !DESCRIPTION:
    ! Get coordinate, dzsoi, and level class metadata for the levmaxurbgrnd dimension
    !
    ! !ARGUMENTS:
    type(file_desc_t), target, intent(inout) :: ncid_source
    type(file_desc_t), target, intent(inout) :: ncid_dest
    type(interp_bounds_type), intent(in) :: bounds_source
    type(interp_bounds_type), intent(in) :: bounds_dest
    character(len=*), intent(in) :: coord_varname
    character(len=*), intent(in) :: dzsoi_varname
    character(len=*), intent(in) :: level_class_varname
    integer, intent(in) :: sgridindex(:)  ! mappings from source to dest points for the appropriate subgrid level (e.g., column-level mappings if this interpolator is for column-level data)

    ! The following output arrays are all allocated in this subroutine:
    real(r8), allocatable, intent(out) :: coordinates_source(:,:)   ! [lev, vec]
    real(r8), allocatable, intent(out) :: coordinates_dest(:,:)     ! [lev, vec]
    integer , allocatable, intent(out) :: level_classes_source(:,:) ! [lev, vec]
    integer , allocatable, intent(out) :: level_classes_dest(:,:)   ! [lev, vec]
    real(r8), allocatable, intent(out) :: dzsoi_source(:,:)         ! [lev, vec]
    real(r8), allocatable, intent(out) :: dzsoi_dest(:,:)           ! [lev, vec]
    !
    ! !LOCAL VARIABLES:
    type(interp_2dvar_type) :: coord_source
    type(interp_2dvar_type) :: coord_dest
    type(interp_2dvar_type) :: dzs_source
    type(interp_2dvar_type) :: dzs_dest
    type(interp_2dvar_type) :: level_class_source
    type(interp_2dvar_type) :: level_class_dest
    real(r8), pointer     :: coord_data_source_sgrid_1d(:)          ! [vec] On the source grid
    real(r8), allocatable :: coord_data_source(:,:)                 ! [vec, lev] Interpolated to the dest grid, but source vertical grid
    real(r8), pointer     :: coord_data_dest(:,:)                   ! [vec, lev] Dest horiz & vertical grid
    real(r8), pointer     :: dzsoi_data_source_sgrid_1d(:)          ! [vec] On the source grid
    real(r8), allocatable :: dzsoi_data_source(:,:)                 ! [vec, lev] Interpolated to the dest grid, but source vertical grid
    real(r8), pointer     :: dzsoi_data_dest(:,:)                   ! [vec, lev] Dest horiz & vertical grid
    integer , pointer     :: level_class_data_source_sgrid_1d(:)    ! [vec] On the source grid
    integer , allocatable :: level_class_data_source(:,:)           ! [vec, lev] Interpolated to the dest grid, but source vertical grid
    integer , pointer     :: level_class_data_dest(:,:)             ! [vec, lev] Dest horiz & vertical grid

    integer :: dimid ! netcdf dimension id
    logical :: dimexist ! whether the given dimension exists on file

    integer :: levmaxurbgrnd_source ! length of levmaxurbgrnd dimension on source file
    integer :: levmaxurbgrnd_dest   ! length of levmaxurbgrnd dimension on dest file

    integer :: beg_dest
    integer :: end_dest
    integer :: beg_source
    integer :: end_source

    integer :: level, g  ! loop indices
    integer :: nlev_source
    logical :: on_source

    character(len=*), parameter :: levmaxurbgrnd_name = 'levmaxurbgrnd'

    character(len=*), parameter :: subname = 'get_levmaxurbgrnd_metadata'
    !-----------------------------------------------------------------------

    ! Get levmaxurbgrnd dimension size on source and dest
    ! BACKWARDS_COMPATIBILITY(wjs, 2020-10-22) On older initial conditions files,
    ! levmaxurbgrnd doesn't exist, but we can use levgrnd for the same purpose because
    ! prior to the existence of levmaxurbgrnd, it was always the case that levgrnd >=
    ! levurb.
    call ncd_inqdid(ncid_source, levmaxurbgrnd_name, dimid, dimexist)
    if (dimexist) then
       call ncd_inqdlen(ncid_source, dimid, levmaxurbgrnd_source)
    else
       call ncd_inqdlen(ncid_source, dimid, levmaxurbgrnd_source, name='levgrnd')
    end if
    ! For dest, we can assume this dimension exists
    call ncd_inqdlen(ncid_dest, dimid, levmaxurbgrnd_dest, name=levmaxurbgrnd_name)

    ! Set coord_data_dest
    coord_dest = interp_2dvar_type( &
         varname = coord_varname, &
         ncid = ncid_dest, &
         file_is_dest = .true., &
         bounds = bounds_dest)
    ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
    ! resolving the generic reference here, giving the message: 'No specific
    ! match can be found for the generic subprogram call "READVAR"'. So we
    ! explicitly call the specific routine, rather than calling readvar.
    call coord_dest%readvar_double(coord_data_dest)
    call shr_assert(coord_dest%get_nlev() == levmaxurbgrnd_dest, &
         msg = 'dest '//coord_varname//' dimension length does not match expected levmaxurbgrnd size', &
         file = sourcefile, &
         line = __LINE__)
    beg_dest = coord_dest%get_vec_beg()
    end_dest = coord_dest%get_vec_end()

    ! Set dzsoi_data_dest
    dzs_dest = interp_2dvar_type( &
         varname = dzsoi_varname, &
         ncid = ncid_dest, &
         file_is_dest = .true., &
         bounds = bounds_dest)
    ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
    ! resolving the generic reference here, giving the message: 'No specific
    ! match can be found for the generic subprogram call "READVAR"'. So we
    ! explicitly call the specific routine, rather than calling readvar.
    call dzs_dest%readvar_double(dzsoi_data_dest)
    call shr_assert(dzs_dest%get_nlev() == levmaxurbgrnd_dest, &
         msg = 'dest '//dzsoi_varname//' dimension length does not match expected levmaxurbgrnd size', &
         file = sourcefile, &
         line = __LINE__)
    SHR_ASSERT_FL(dzs_dest%get_vec_beg() == beg_dest, sourcefile, __LINE__)
    SHR_ASSERT_FL(dzs_dest%get_vec_end() == end_dest, sourcefile, __LINE__)

    ! Set level_class_data_dest
    level_class_dest = interp_2dvar_type( &
         varname = level_class_varname, &
         ncid = ncid_dest, &
         file_is_dest = .true., &
         bounds = bounds_dest)
    ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
    ! resolving the generic reference here, giving the message: 'No specific
    ! match can be found for the generic subprogram call "READVAR"'. So we
    ! explicitly call the specific routine, rather than calling readvar.
    call level_class_dest%readvar_int(level_class_data_dest)
    call shr_assert(level_class_dest%get_nlev() == levmaxurbgrnd_dest, &
         msg = 'dest '//level_class_varname//' dimension length does not match expected levmaxurbgrnd size', &
         file = sourcefile, &
         line = __LINE__)
    SHR_ASSERT_FL(level_class_dest%get_vec_beg() == beg_dest, sourcefile, __LINE__)
    SHR_ASSERT_FL(level_class_dest%get_vec_end() == end_dest, sourcefile, __LINE__)

    ! NOTE(wjs, 2015-10-18) The following check is helpful while we still have old initial
    ! conditions files that do not have the necessary metadata. Once these old initial
    ! conditions files have been phased out, we can remove this check. (Without this
    ! check, the run will still abort if it can't find the necessary variables - it just
    ! won't have a very helpful error message.)
    call interp_levgrnd_check_source_file(ncid_source, coord_varname)

    ! Set coord_data_source
    coord_source = interp_2dvar_type( &
         varname = coord_varname, &
         ncid = ncid_source, &
         file_is_dest = .false., &
         bounds = bounds_source)
    call shr_assert(coord_source%get_nlev() == levmaxurbgrnd_source, &
         msg = 'source '//coord_varname//' dimension length does not match expected levmaxurbgrnd size', &
         file = sourcefile, &
         line = __LINE__)
    nlev_source = coord_source%get_nlev()
    beg_source = coord_source%get_vec_beg()
    end_source = coord_source%get_vec_end()
    allocate(coord_data_source(beg_dest:end_dest, nlev_source))
    allocate(coord_data_source_sgrid_1d(beg_source:end_source))
    do level = 1, nlev_source
       ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
       ! resolving the generic reference here, giving the message: 'No specific
       ! match can be found for the generic subprogram call "READLEVEL"'. So we
       ! explicitly call the specific routine, rather than calling readlevel.
       call coord_source%readlevel_double(coord_data_source_sgrid_1d, level)
       call interp_1d_data( &
            begi = beg_source, endi = end_source, &
            bego = beg_dest,   endo = end_dest, &
            sgridindex = sgridindex, &
            keep_existing = .false., &
            data_in = coord_data_source_sgrid_1d, &
            data_out = coord_data_source(:,level))
    end do
    deallocate(coord_data_source_sgrid_1d)

    ! The following call checks whether the finidat source file contains
    ! DZSOI, necessary for the vertical interpolation with scaling by layer
    ! thickness. Scaling by thickness is necessary only if the run uses a
    ! different vertical profile than the one found in the finidata file.
    call check_var(ncid_source, dzsoi_varname, on_source)
    allocate(dzsoi_data_source(beg_dest:end_dest, nlev_source))
    if (on_source) then
       ! Set dzsoi_data_source
       dzs_source = interp_2dvar_type( &
            varname = dzsoi_varname, &
            ncid = ncid_source, &
            file_is_dest = .false., &
            bounds = bounds_source)
       call shr_assert(dzs_source%get_nlev() == levmaxurbgrnd_source, &
            msg = 'source '//dzsoi_varname//' dimension length does not match expected levmaxurbgrnd size', &
            file = sourcefile, &
            line = __LINE__)
       SHR_ASSERT_FL(dzs_source%get_vec_beg() == beg_source, sourcefile, __LINE__)
       SHR_ASSERT_FL(dzs_source%get_vec_end() == end_source, sourcefile, __LINE__)
       allocate(dzsoi_data_source_sgrid_1d(beg_source:end_source))
       do level = 1, nlev_source
          ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
          ! resolving the generic reference here, giving the message: 'No specific
          ! match can be found for the generic subprogram call "READLEVEL"'. So we
          ! explicitly call the specific routine, rather than calling readlevel.
          call dzs_source%readlevel_double(dzsoi_data_source_sgrid_1d, level)
          call interp_1d_data( &
               begi = beg_source, endi = end_source, &
               bego = beg_dest,   endo = end_dest, &
               sgridindex = sgridindex, &
               keep_existing = .false., &
               data_in = dzsoi_data_source_sgrid_1d, &
               data_out = dzsoi_data_source(:,level))
       end do
       deallocate(dzsoi_data_source_sgrid_1d)
    else if (levmaxurbgrnd_source == levmaxurbgrnd_dest) then
       ! BACKWARDS_COMPATIBILITY (slevis, 2021-04-22)
       ! When same number of layers and .not. on_source, we skip the calls to
       ! (a) interp_levgrnd_check_source_file(ncid_source, dzsoi_varname)
       ! (b) dzs_source%readlevel_double(dzsoi_data_source...
       ! to avoid stopping or crashing the model due to missing dzsoi
       ! in the finidat file. We do this for backwards compatibility when the
       ! soil profile in the run is the same as in the finidat file and the
       ! finidat file does not include dzsoi, yet.
       !
       ! Skipping these calls is a problem only in the off chance that the
       ! source and dest soil profiles differ despite having the same number of
       ! layers. The model will handle such cases without scaling
       ! the finidat variables that need scaling by soil layer thickness.
       do level = 1, nlev_source
          do g = beg_dest, end_dest
             dzsoi_data_source(g,level) = dzsoi_data_dest(g,level)
          end do
       end do
    else  ! different number of layers and .not. on_source
       ! BACKWARDS_COMPATIBILITY (slevis, 2021-04-22)
       ! Rather than failing when DZSOI is not present on the finidat file,
       ! construct dzsoi_data_source from coord_data_source
       do g = beg_dest, end_dest
          dzsoi_data_source(g,1) = 2._r8 * coord_data_source(g,1)
          if (nlev_source > 1) then
             do level = 2, nlev_source
                dzsoi_data_source(g,level) = 2._r8 * &
                   (coord_data_source(g,level) - &
                    sum(dzsoi_data_source(g,1:level-1)))
             end do
          end if
       end do
    end if

    ! NOTE(wjs, 2015-10-18) The following check is helpful while we still have old initial
    ! conditions files that do not have the necessary metadata. Once these old initial
    ! conditions files have been phased out, we can remove this check. (Without this
    ! check, the run will still abort if it can't find the necessary variables - it just
    ! won't have a very helpful error message.)
    call interp_levgrnd_check_source_file(ncid_source, level_class_varname)

    ! Set level_class_data_source
    level_class_source = interp_2dvar_type( &
         varname = level_class_varname, &
         ncid = ncid_source, &
         file_is_dest = .false., &
         bounds = bounds_source)
    call shr_assert(level_class_source%get_nlev() == levmaxurbgrnd_source, &
         msg = 'source '//level_class_varname//' dimension length does not match expected levmaxurbgrnd size', &
         file = sourcefile, &    
         line = __LINE__)
    SHR_ASSERT_FL(level_class_source%get_vec_beg() == beg_source, sourcefile, __LINE__)
    SHR_ASSERT_FL(level_class_source%get_vec_end() == end_source, sourcefile, __LINE__)
    allocate(level_class_data_source(beg_dest:end_dest, nlev_source))
    allocate(level_class_data_source_sgrid_1d(beg_source:end_source))
    do level = 1, nlev_source
       ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
       ! resolving the generic reference here, giving the message: 'No specific
       ! match can be found for the generic subprogram call "READLEVEL"'. So we
       ! explicitly call the specific routine, rather than calling readlevel.
       call level_class_source%readlevel_int(level_class_data_source_sgrid_1d, level)
       call interp_1d_data( &
            begi = beg_source, endi = end_source, &
            bego = beg_dest, endo = end_dest, &
            sgridindex = sgridindex, &
            keep_existing = .false., &
            data_in = level_class_data_source_sgrid_1d, &
            data_out = level_class_data_source(:,level))
    end do
    deallocate(level_class_data_source_sgrid_1d)

    ! Set output arrays
    call transpose_wrapper(coordinates_source, coord_data_source)
    call transpose_wrapper(coordinates_dest, coord_data_dest)
    call transpose_wrapper(dzsoi_source, dzsoi_data_source)
    call transpose_wrapper(dzsoi_dest, dzsoi_data_dest)
    call transpose_wrapper(level_classes_source, level_class_data_source)
    call transpose_wrapper(level_classes_dest, level_class_data_dest)

    ! Deallocate pointers (allocatables are automatically deallocated)
    deallocate(coord_data_dest)
    deallocate(dzsoi_data_dest)
    deallocate(level_class_data_dest)

  end subroutine get_levmaxurbgrnd_metadata

  !-----------------------------------------------------------------------
  subroutine interp_levgrnd_check_source_file(ncid_source, varname)
    !
    ! !DESCRIPTION:
    ! Ensure that the necessary variables are present on the source file for the levgrnd
    ! interpolator. In separate calls to this subroutine, we check for
    ! - coord
    ! - dzsoi
    ! - level_class
    !
    ! Aborts the run with a useful error message if variable is missing.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid_source
    character(len=*) , intent(in) :: varname
    !
    ! !LOCAL VARIABLES:
    logical :: on_source
    character(len=:), allocatable :: missing

    character(len=*), parameter :: subname = 'interp_levgrnd_check_source_file'
    !-----------------------------------------------------------------------

    missing = ' '
    call check_var(ncid_source, varname, on_source)
    if (.not. on_source) then
       missing = missing // varname // ' '
    end if
    if (missing /= ' ') then
       if (masterproc) then
          write(iulog,*) subname//&
               ' ERROR: source file for init_interp is missing the necessary variable(s):'
          write(iulog,*) missing
          write(iulog,*) 'To solve this problem, run the model for a short time using this tag,'
          write(iulog,*) 'with a configuration that matches the source file, using the source'
          write(iulog,*) 'file as finidat (with use_init_interp = .false.), in order to'
          write(iulog,*) 'produce a new restart file with the necessary metadata.'
          write(iulog,*) 'Then use that new file as the finidat file for init_interp.'
          write(iulog,*) ' '
          write(iulog,*) 'If that is not possible, then an alternative is to run the model for'
          write(iulog,*) 'a short time using this tag, with cold start initial conditions'
          write(iulog,*) '(finidat = " "). Then use a tool like ncks to copy the misssing fields'
          write(iulog,*) 'onto the original source finidat file. Then use that patched file'
          write(iulog,*) 'as the finidat file for init_interp.'
       end if

       call endrun(subname//' ERROR: source file for init_interp is missing '// &
            missing)
    end if

  end subroutine interp_levgrnd_check_source_file

  !-----------------------------------------------------------------------
  subroutine create_snow_interpolators(interp_multilevel_levsno, interp_multilevel_levsno1, &
       ncid_source, bounds_source, bounds_dest, colindex)
    !
    ! !DESCRIPTION:
    ! Create multi-level interpolators for snow variables
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(interp_multilevel_snow_type), intent(out) :: interp_multilevel_levsno
    type(interp_multilevel_snow_type), intent(out) :: interp_multilevel_levsno1
    type(file_desc_t), intent(inout) :: ncid_source ! netcdf ID for source file
    type(interp_bounds_type), intent(in) :: bounds_source
    type(interp_bounds_type), intent(in) :: bounds_dest
    integer, intent(in) :: colindex(:)  ! mappings from source to dest for column-level arrays
    !
    ! !LOCAL VARIABLES:
    ! snlsno_source needs to be a pointer to satisfy the interface of ncd_io
    integer, pointer     :: snlsno_source_sgrid(:)  ! snlsno in source, on source grid
    integer, allocatable :: snlsno_source(:)        ! snlsno_source interpolated to dest
    integer, allocatable :: snlsno_source_plus_1(:) ! snlsno_source+1 interpolated to dest

    character(len=*), parameter :: subname = 'create_snow_interpolators'
    !-----------------------------------------------------------------------

    ! Read snlsno_source_sgrid
    allocate(snlsno_source_sgrid(bounds_source%get_begc() : bounds_source%get_endc()))
    call ncd_io(ncid=ncid_source, varname='SNLSNO', flag='read', &
         data=snlsno_source_sgrid)
    snlsno_source_sgrid(:) = abs(snlsno_source_sgrid(:))

    ! Interpolate to dest
    allocate(snlsno_source(bounds_dest%get_begc() : bounds_dest%get_endc()))
    call interp_1d_data( &
         begi = bounds_source%get_begc(), endi = bounds_source%get_endc(), &
         bego = bounds_dest%get_begc(), endo = bounds_dest%get_endc(), &
         sgridindex = colindex, &
         keep_existing = .false., &
         data_in = snlsno_source_sgrid, data_out = snlsno_source)
    deallocate(snlsno_source_sgrid)

    ! Set up interp_multilevel_levsno
    interp_multilevel_levsno = interp_multilevel_snow_type( &
         num_snow_layers_source = snlsno_source, &
         num_layers_name = 'SNLSNO')

    ! Set up interp_multilevel_levsno1
    !
    ! For variables dimensioned (levsno+1), we assume they have (snlsno+1) active layers.
    ! Thus, if there are 0 active layers in the source, the bottom layer's value will
    ! still get copied for these (levsno+1) variables.
    allocate(snlsno_source_plus_1(bounds_dest%get_begc() : bounds_dest%get_endc()))
    snlsno_source_plus_1(:) = snlsno_source(:) + 1
    interp_multilevel_levsno1 = interp_multilevel_snow_type( &
         num_snow_layers_source = snlsno_source_plus_1, &
         num_layers_name = 'SNLSNO+1')

    deallocate(snlsno_source)
    deallocate(snlsno_source_plus_1)

  end subroutine create_snow_interpolators

  ! ========================================================================
  ! Finalizers
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine destroy_interp_multilevel_container_type(this)
    !
    ! !DESCRIPTION:
    ! Finalize routine for interp_multilevel_container_type
    !
    ! !ARGUMENTS:
    type(interp_multilevel_container_type) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'destroy_interp_multilevel_container_type'
    !-----------------------------------------------------------------------

    deallocate(this%interp_multilevel_copy)
    deallocate(this%interp_multilevel_levgrnd_col)
    deallocate(this%interp_multilevel_levmaxurbgrnd_col)
    deallocate(this%interp_multilevel_levgrnd_pft)
    deallocate(this%interp_multilevel_levsno)
    deallocate(this%interp_multilevel_levsno1)
    deallocate(this%interp_multilevel_levtot_col)

  end subroutine destroy_interp_multilevel_container_type


end module initInterpMultilevelContainer
