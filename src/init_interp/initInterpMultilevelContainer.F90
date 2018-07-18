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
  use ncdio_pio                  , only : file_desc_t, var_desc_t, check_var, ncd_io
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
     type(interp_multilevel_interp_type), pointer :: interp_multilevel_levgrnd_pft
     type(interp_multilevel_snow_type), pointer   :: interp_multilevel_levsno
     type(interp_multilevel_snow_type), pointer   :: interp_multilevel_levsno1
     type(interp_multilevel_split_type), pointer  :: interp_multilevel_levtot_col
   contains
     procedure :: find_interpolator
  end type interp_multilevel_container_type

  interface interp_multilevel_container_type
     module procedure constructor
  end interface interp_multilevel_container_type

  ! Private routines

  private :: create_interp_multilevel_levgrnd
  private :: interp_levgrnd_check_source_file
  private :: create_snow_interpolators

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor(ncid_source, ncid_dest, bounds_source, bounds_dest, &
       pftindex, colindex) result(this)
    !
    ! !DESCRIPTION:
    ! Create an interp_multilevel_container_type instance.
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t
    ! 
    ! !ARGUMENTS:
    type(interp_multilevel_container_type) :: this  ! function result
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

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    allocate(this%interp_multilevel_copy)
    this%interp_multilevel_copy = interp_multilevel_copy_type()

    allocate(this%interp_multilevel_levgrnd_col)
    this%interp_multilevel_levgrnd_col = create_interp_multilevel_levgrnd( &
         ncid_source = ncid_source, &
         ncid_dest = ncid_dest, &
         bounds_source = bounds_source, &
         bounds_dest = bounds_dest, &
         coord_varname = 'COL_Z', &
         level_class_varname = 'LEVGRND_CLASS', &
         sgridindex = colindex)

    allocate(this%interp_multilevel_levgrnd_pft)
    this%interp_multilevel_levgrnd_pft = create_interp_multilevel_levgrnd( &
         ncid_source = ncid_source, &
         ncid_dest = ncid_dest, &
         bounds_source = bounds_source, &
         bounds_dest = bounds_dest, &
         coord_varname = 'COL_Z_p', &
         level_class_varname = 'LEVGRND_CLASS_p', &
         sgridindex = pftindex)

    allocate(this%interp_multilevel_levsno)
    allocate(this%interp_multilevel_levsno1)
    call create_snow_interpolators( &
         interp_multilevel_levsno = this%interp_multilevel_levsno, &
         interp_multilevel_levsno1 = this%interp_multilevel_levsno1, &
         ncid_source = ncid_source, &
         bounds_source = bounds_source, &
         bounds_dest = bounds_dest, &
         colindex = colindex)

    ! levtot is two sets of levels: first snow, then levgrnd
    allocate(this%interp_multilevel_levtot_col)
    this%interp_multilevel_levtot_col = create_interp_multilevel_split_type( &
         interpolator_first_levels = this%find_interpolator('levsno', 'column'), &
         interpolator_second_levels = this%interp_multilevel_levgrnd_col, &
         num_second_levels_source = this%interp_multilevel_levgrnd_col%get_nlev_source(), &
         num_second_levels_dest = this%interp_multilevel_levgrnd_col%get_nlev_dest())

  end function constructor

  ! ========================================================================
  ! Public methods
  ! ========================================================================

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
  function create_interp_multilevel_levgrnd(ncid_source, ncid_dest, &
       bounds_source, bounds_dest, &
       coord_varname, level_class_varname, &
       sgridindex) &
       result(interpolator)
    !
    ! !DESCRIPTION:
    ! Create the interpolator used to interpolate variables dimensioned by levgrnd
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(interp_multilevel_interp_type) :: interpolator  ! function result
    type(file_desc_t), target, intent(inout) :: ncid_source
    type(file_desc_t), target, intent(inout) :: ncid_dest
    type(interp_bounds_type), intent(in) :: bounds_source
    type(interp_bounds_type), intent(in) :: bounds_dest
    character(len=*), intent(in) :: coord_varname
    character(len=*), intent(in) :: level_class_varname
    integer, intent(in) :: sgridindex(:)  ! mappings from source to dest points for the appropriate subgrid level (e.g., column-level mappings if this interpolator is for column-level data)
    !
    ! !LOCAL VARIABLES:
    type(interp_2dvar_type) :: coord_source
    type(interp_2dvar_type) :: coord_dest
    type(interp_2dvar_type) :: level_class_source
    type(interp_2dvar_type) :: level_class_dest
    real(r8), pointer     :: coord_data_source_sgrid_1d(:)          ! [vec] On the source grid
    real(r8), allocatable :: coord_data_source(:,:)                 ! [vec, lev] Interpolated to the dest grid, but source vertical grid
    real(r8), pointer     :: coord_data_dest(:,:)                   ! [vec, lev] Dest horiz & vertical grid
    integer , pointer     :: level_class_data_source_sgrid_1d(:)    ! [vec] On the source grid
    integer , allocatable :: level_class_data_source(:,:)           ! [vec, lev] Interpolated to the dest grid, but source vertical grid
    integer , pointer     :: level_class_data_dest(:,:)             ! [vec, lev] Dest horiz & vertical grid
    real(r8), allocatable :: coord_data_source_transpose(:,:)       ! [lev, vec]
    real(r8), allocatable :: coord_data_dest_transpose(:,:)         ! [lev, vec]
    integer , allocatable :: level_class_data_source_transpose(:,:) ! [lev, vec]
    integer , allocatable :: level_class_data_dest_transpose(:,:)   ! [lev, vec]

    integer :: beg_dest
    integer :: end_dest
    integer :: beg_source
    integer :: end_source

    integer :: level
    integer :: nlev_source

    character(len=*), parameter :: subname = 'create_interp_multilevel_levgrnd'
    !-----------------------------------------------------------------------

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
    beg_dest = coord_dest%get_vec_beg()
    end_dest = coord_dest%get_vec_end()

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
    SHR_ASSERT(level_class_dest%get_vec_beg() == beg_dest, errMsg(sourcefile, __LINE__))
    SHR_ASSERT(level_class_dest%get_vec_end() == end_dest, errMsg(sourcefile, __LINE__))

    ! NOTE(wjs, 2015-10-18) The following check is helpful while we still have old initial
    ! conditions files that do not have the necessary metadata. Once these old initial
    ! conditions files have been phased out, we can remove this check. (Without this
    ! check, the run will still abort if it can't find the necessary variables - it just
    ! won't have a very helpful error message.)
    call interp_levgrnd_check_source_file(ncid_source, coord_varname, level_class_varname)

    ! Set coord_data_source
    coord_source = interp_2dvar_type( &
         varname = coord_varname, &
         ncid = ncid_source, &
         file_is_dest = .false., &
         bounds = bounds_source)
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

    ! Set level_class_data_source
    level_class_source = interp_2dvar_type( &
         varname = level_class_varname, &
         ncid = ncid_source, &
         file_is_dest = .false., &
         bounds = bounds_source)
    SHR_ASSERT(level_class_source%get_nlev() == nlev_source, errMsg(sourcefile, __LINE__))
    SHR_ASSERT(level_class_source%get_vec_beg() == beg_source, errMsg(sourcefile, __LINE__))
    SHR_ASSERT(level_class_source%get_vec_end() == end_source, errMsg(sourcefile, __LINE__))
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

    ! Create interpolator
    call transpose_wrapper(coord_data_source_transpose, coord_data_source)
    call transpose_wrapper(coord_data_dest_transpose, coord_data_dest)
    call transpose_wrapper(level_class_data_source_transpose, level_class_data_source)
    call transpose_wrapper(level_class_data_dest_transpose, level_class_data_dest)
    interpolator = interp_multilevel_interp_type( &
         coordinates_source = coord_data_source_transpose, &
         coordinates_dest = coord_data_dest_transpose, &
         level_classes_source = level_class_data_source_transpose, &
         level_classes_dest = level_class_data_dest_transpose, &
         coord_varname = coord_varname)

    ! Deallocate pointers (allocatables are automatically deallocated)
    deallocate(coord_data_dest)
    deallocate(level_class_data_dest)

  end function create_interp_multilevel_levgrnd

  !-----------------------------------------------------------------------
  subroutine interp_levgrnd_check_source_file(ncid_source, coord_varname, level_class_varname)
    !
    ! !DESCRIPTION:
    ! Ensure that the necessary variables are present on the source file for the levgrnd
    ! interpolator.
    !
    ! Aborts the run with a useful error message if either variable is missing
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid_source
    character(len=*) , intent(in) :: coord_varname
    character(len=*) , intent(in) :: level_class_varname
    !
    ! !LOCAL VARIABLES:
    logical :: coord_on_source
    logical :: level_class_on_source
    type(var_desc_t) :: coord_source_vardesc  ! unused, but needed for check_var interface
    type(var_desc_t) :: level_class_source_vardesc  ! unused, but needed for check_var interface
    character(len=:), allocatable :: variables_missing

    character(len=*), parameter :: subname = 'interp_levgrnd_check_source_file'
    !-----------------------------------------------------------------------

    variables_missing = ' '
    call check_var(ncid_source, coord_varname, coord_source_vardesc, coord_on_source)
    if (.not. coord_on_source) then
       variables_missing = variables_missing // coord_varname // ' '
    end if
    call check_var(ncid_source, level_class_varname, level_class_source_vardesc, level_class_on_source)
    if (.not. level_class_on_source) then
       variables_missing = variables_missing // level_class_varname // ' '
    end if
    if (variables_missing /= ' ') then
       if (masterproc) then
          write(iulog,*) subname//&
               ' ERROR: source file for init_interp is missing the necessary variable(s):'
          write(iulog,*) variables_missing
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
            variables_missing)
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


end module initInterpMultilevelContainer
