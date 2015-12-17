module glcBehaviorMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Determines a number of aspects of the behavior of glacier_mec classes in each grid
  ! cell.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use landunit_varcon, only : istice_mec
  use clm_instur     , only : wt_lunit, wt_glc_mec
  use decompMod      , only : bounds_type
  use filterColMod   , only : filter_col_type
  use ColumnType     , only : col

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: glc_behavior_type
     private

     ! ------------------------------------------------------------------------
     ! Public data
     ! ------------------------------------------------------------------------

     ! If has_virtual_columns_grc(g) is true, then grid cell g has virtual columns for
     ! all possible glc_mec columns.
     !
     ! In principle, this should only be needed within the icemask, where we need virtual
     ! columns for the sake of coupling with CISM. This is needed in order to (1) provide
     ! SMB in all elevation classes, in case it is being used with 1-way coupling (or to
     ! force a later TG run); (2) even with two-way coupling, provide SMB in the
     ! elevation classes above and below existing elevation classes, for the sake of
     ! vertical interpolation; (3) provide place-holder columns (which are already
     ! spun-up) for dynamic landunits. 
     !
     ! However, by making this part of the user-modifiable "glc behavior", we make it easy
     ! for the user to add virtual columns, if this is desired for diagnostic
     ! purposes. (Also, we cannot use icemask for all purposes, because it isn't known at
     ! initialization.)
     logical, allocatable, public :: has_virtual_columns_grc(:)

     ! ------------------------------------------------------------------------
     ! Private data
     ! ------------------------------------------------------------------------

     ! If collapse_to_atm_topo_grc(g) is true, then grid cell g has at most one glc_mec
     ! column, whose topographic height exactly matches the atmosphere's topographic
     ! height for that grid cell (so that there is no adjustment of atmospheric
     ! forcings).
     !
     ! Note that has_virtual_columns_grc(g) is guaranteed to be false if
     ! collapse_to_atm_topo_grc(g) is true.
     logical, allocatable :: collapse_to_atm_topo_grc(:)

   contains

     ! ------------------------------------------------------------------------
     ! Public routines
     ! ------------------------------------------------------------------------

     procedure, public  :: Init
     procedure, public  :: InitForTesting  ! version of Init meant for unit testing

     ! get number of subgrid units in glc_mec landunit on one grid cell
     procedure, public  :: get_num_glc_mec_subgrid

     ! returns true if memory should be allocated for the given glc_mec column, and its
     ! weight on the landunit
     procedure, public  :: glc_mec_col_exists

     ! returns true if glc_mec columns on the given grid cell have dynamic type (type
     ! potentially changing at runtime)
     procedure, public  :: cols_have_dynamic_type

     ! update topographic height and class of glc_mec columns in regions where these are
     ! collapsed to a single column
     procedure, public  :: update_collapsed_columns

     ! update the column class types of any glc_mec columns that need to be updated
     procedure, public  :: update_glc_classes

     ! ------------------------------------------------------------------------
     ! Private routines
     ! ------------------------------------------------------------------------

     procedure, private :: InitAllocate
     procedure, private :: set_masks

     ! returns a column-level filter of ice_mec columns with the collapse_to_atm_topo
     ! behavior
     procedure, private :: collapse_to_atm_topo_icemec_filterc

     ! update class of glc_mec columns in regions where these are collapsed to a single
     ! column, given a filter
     procedure, private :: update_collapsed_columns_classes

  end type glc_behavior_type

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, begg, endg, glcmask, latdeg)
    !
    ! !DESCRIPTION:
    ! Initialize a glc_behavior_type object
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(inout) :: this
    integer, intent(in) :: begg  ! beginning gridcell index
    integer, intent(in) :: endg  ! ending gridcell index
    integer, intent(in) :: glcmask(begg: )
    real(r8), intent(in) :: latdeg(begg: )
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    call this%InitAllocate(begg, endg)
    call this%set_masks(begg, endg, glcmask(begg:endg), latdeg(begg:endg))

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitForTesting(this, begg, endg, &
       has_virtual_columns, collapse_to_atm_topo)
    !
    ! !DESCRIPTION:
    ! Initialize a glc_behavior_type object for testing
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(inout) :: this
    integer, intent(in) :: begg  ! beginning gridcell index
    integer, intent(in) :: endg  ! ending gridcell index
    logical, intent(in) :: has_virtual_columns(begg:)
    logical, intent(in) :: collapse_to_atm_topo(begg:)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitForTesting'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(has_virtual_columns) == (/endg/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(collapse_to_atm_topo) == (/endg/)), errMsg(__FILE__, __LINE__))

    call this%InitAllocate(begg, endg)
    this%has_virtual_columns_grc(:) = has_virtual_columns(:)
    this%collapse_to_atm_topo_grc(:) = collapse_to_atm_topo(:)

  end subroutine InitForTesting


  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, begg, endg)
    !
    ! !DESCRIPTION:
    ! Allocate variables in this object
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(inout) :: this
    integer, intent(in) :: begg  ! beginning gridcell index
    integer, intent(in) :: endg  ! ending gridcell index
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    allocate(this%has_virtual_columns_grc (begg:endg)); this%has_virtual_columns_grc (:) = .false.
    allocate(this%collapse_to_atm_topo_grc(begg:endg)); this%collapse_to_atm_topo_grc(:) = .false.

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine get_num_glc_mec_subgrid(this, gi, atm_topo, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Get number of subgrid units in glc_mec landunit on one grid cell
    !
    ! !USES:
    use clm_varpar      , only : maxpatch_glcmec
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(in) :: this
    integer , intent(in)  :: gi       ! grid cell index
    real(r8), intent(in)  :: atm_topo ! atmosphere's topographic height for this grid cell (m)
    integer , intent(out) :: npatches ! number of glacier_mec patches in this grid cell
    integer , intent(out) :: ncols    ! number of glacier_mec columns in this grid cell
    integer , intent(out) :: nlunits  ! number of glacier_mec landunits in this grid cell
    !
    ! !LOCAL VARIABLES:
    integer  :: m  ! loop index
    logical  :: col_exists
    real(r8) :: col_wt_lunit

    character(len=*), parameter :: subname = 'get_num_glc_mec_subgrid'
    !-----------------------------------------------------------------------

    ncols = 0

    do m = 1, maxpatch_glcmec
       call this%glc_mec_col_exists(gi = gi, elev_class = m, atm_topo = atm_topo, &
            exists = col_exists, col_wt_lunit = col_wt_lunit)
       if (col_exists) then
          ncols = ncols + 1
       end if
    end do

    if (this%collapse_to_atm_topo_grc(gi) .and. &
         wt_lunit(gi, istice_mec) > 0.0_r8) then
       ! For grid cells with the collapse_to_atm_topo behavior, with a non-zero weight
       ! ice_mec landunit, we expect exactly one column
       SHR_ASSERT(ncols == 1, errMsg(__FILE__, __LINE__))
    end if

    if (ncols > 0) then
       npatches = ncols
       nlunits = 1
    else
       npatches = 0
       nlunits = 0
    end if
  
  end subroutine get_num_glc_mec_subgrid

  !-----------------------------------------------------------------------
  subroutine glc_mec_col_exists(this, gi, elev_class, atm_topo, exists, col_wt_lunit)
    !
    ! !DESCRIPTION:
    ! For the given glc_mec column, with elevation class index elev_class, in grid cell
    ! gi: sets exists to true if memory should be allocated for this column, and sets
    ! col_wt_lunit to the column's weight on the icemec landunit.
    !
    ! If exists is false, then col_wt_lunit is arbitrary and should be ignored.
    !
    ! !USES:
    use glc_elevclass_mod, only : glc_get_elevation_class, GLC_ELEVCLASS_ERR_NONE
    use glc_elevclass_mod, only : GLC_ELEVCLASS_ERR_TOO_LOW, GLC_ELEVCLASS_ERR_TOO_HIGH
    use glc_elevclass_mod, only : glc_errcode_to_string
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(in) :: this
    integer,  intent(in)  :: gi           ! grid cell index
    integer,  intent(in)  :: elev_class   ! elevation class index
    real(r8), intent(in)  :: atm_topo     ! atmosphere's topographic height for this grid cell (m)
    logical,  intent(out) :: exists       ! whether memory should be allocated for this column
    real(r8), intent(out) :: col_wt_lunit ! column's weight on the icemec landunit
    !
    ! !LOCAL VARIABLES:
    integer :: atm_elev_class ! elevation class corresponding to atmosphere topographic height
    integer :: err_code

    character(len=*), parameter :: subname = 'glc_mec_col_exists'
    !-----------------------------------------------------------------------

    ! Set default outputs
    exists = .false.
    col_wt_lunit = wt_glc_mec(gi, elev_class)

    if (this%collapse_to_atm_topo_grc(gi)) then
       if (wt_lunit(gi, istice_mec) > 0.0_r8) then
          call glc_get_elevation_class(atm_topo, atm_elev_class, err_code)
          if ( err_code == GLC_ELEVCLASS_ERR_NONE .or. &
               err_code == GLC_ELEVCLASS_ERR_TOO_LOW .or. &
               err_code == GLC_ELEVCLASS_ERR_TOO_HIGH) then
             ! These are all acceptable "errors" - it is even okay for these purposes if
             ! the elevation is lower than the lower bound of elevation class 1, or
             ! higher than the upper bound of the top elevation class.

             ! Do nothing
          else
             write(iulog,*) subname, ': ERROR getting elevation class for topo = ', atm_topo
             write(iulog,*) glc_errcode_to_string(err_code)
             call endrun(msg=subname//': ERROR getting elevation class')
          end if

          if (elev_class == atm_elev_class) then
             exists = .true.
             col_wt_lunit = 1._r8
          else
             exists = .false.
             col_wt_lunit = 0._r8
          end if
       end if

    else  ! collapse_to_atm_topo_grc .false.
       if (this%has_virtual_columns_grc(gi)) then
          exists = .true.
       else if (wt_lunit(gi, istice_mec) > 0.0_r8 .and. &
            wt_glc_mec(gi, elev_class) > 0.0_r8) then
          ! If the landunit has non-zero weight on the grid cell, and this column has
          ! non-zero weight on the landunit...
          exists = .true.
       end if
    end if

  end subroutine glc_mec_col_exists

  !-----------------------------------------------------------------------
  function cols_have_dynamic_type(this, gi)
    !
    ! !DESCRIPTION:
    ! Returns true if glc_mec columns on the given grid cell have dynamic type (i.e.,
    ! type potentially changing at runtime)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical :: cols_have_dynamic_type  ! function result
    class(glc_behavior_type), intent(in) :: this
    integer, intent(in) :: gi  ! grid cell index
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'cols_have_dynamic_type'
    !-----------------------------------------------------------------------

    if (this%collapse_to_atm_topo_grc(gi)) then
       cols_have_dynamic_type = .true.
    else
       cols_have_dynamic_type = .false.
    end if

  end function cols_have_dynamic_type

  !-----------------------------------------------------------------------
  subroutine update_collapsed_columns(this, bounds, atm_topo)
    !
    ! !DESCRIPTION:
    ! Update topographic height and class of glc_mec columns in regions where these are
    ! collapsed to a single column, whose topographic height matches the atmosphere's
    ! topographic height.
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in) :: atm_topo(bounds%begg:)  ! atmosphere's topographic height (m)
    !
    ! !LOCAL VARIABLES:
    type(filter_col_type) :: collapse_filterc
    integer :: fc         ! filter index
    integer :: c          ! column index
    integer :: g          ! grid cell index

    character(len=*), parameter :: subname = 'update_collapsed_columns'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(atm_topo) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))

    collapse_filterc = this%collapse_to_atm_topo_icemec_filterc(bounds)

    do fc = 1, collapse_filterc%num
       c = collapse_filterc%indices(fc)
       g = col%gridcell(c)

       col%glc_topo(c) = atm_topo(g)
    end do

    call this%update_collapsed_columns_classes(collapse_filterc)

  end subroutine update_collapsed_columns

  !-----------------------------------------------------------------------
  subroutine update_glc_classes(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update the column class types of any glc_mec columns that need to be updated.
    !
    ! Assumes that col%glc_topo has already been set appropriately.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    type(filter_col_type) :: collapse_filterc

    character(len=*), parameter :: subname = 'update_glc_classes'
    !-----------------------------------------------------------------------

    collapse_filterc = this%collapse_to_atm_topo_icemec_filterc(bounds)
    call this%update_collapsed_columns_classes(collapse_filterc)

  end subroutine update_glc_classes

  !-----------------------------------------------------------------------
  subroutine update_collapsed_columns_classes(this, collapse_filterc)
    !
    ! !DESCRIPTION:
    ! Update class of glc_mec columns in regions where these are collapsed to a single
    ! column, given a filter.
    !
    ! Assumes that glc_topo has already been updated appropriately for these columns.
    !
    ! !USES:
    use glc_elevclass_mod, only : glc_get_elevation_class, GLC_ELEVCLASS_ERR_NONE
    use glc_elevclass_mod, only : GLC_ELEVCLASS_ERR_TOO_LOW, GLC_ELEVCLASS_ERR_TOO_HIGH
    use glc_elevclass_mod, only : glc_errcode_to_string
    use column_varcon    , only : icemec_class_to_col_itype
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(in) :: this
    type(filter_col_type), intent(in) :: collapse_filterc
    !
    ! !LOCAL VARIABLES:
    integer :: fc         ! filter index
    integer :: c          ! column index
    integer :: elev_class ! elevation class of the single column on the ice_mec landunit
    integer :: err_code

    character(len=*), parameter :: subname = 'update_collapsed_columns_classes'
    !-----------------------------------------------------------------------

    do fc = 1, collapse_filterc%num
       c = collapse_filterc%indices(fc)

       call glc_get_elevation_class(col%glc_topo(c), elev_class, err_code)
       if ( err_code == GLC_ELEVCLASS_ERR_NONE .or. &
            err_code == GLC_ELEVCLASS_ERR_TOO_LOW .or. &
            err_code == GLC_ELEVCLASS_ERR_TOO_HIGH) then
          ! These are all acceptable "errors" - it is even okay for these purposes if
          ! the elevation is lower than the lower bound of elevation class 1, or
          ! higher than the upper bound of the top elevation class.
          
          ! Do nothing
       else
          write(iulog,*) subname, ': ERROR getting elevation class for topo = ', &
               col%glc_topo(c)
          write(iulog,*) glc_errcode_to_string(err_code)
          call endrun(msg=subname//': ERROR getting elevation class')
       end if

       call col%update_itype(c = c, itype = icemec_class_to_col_itype(elev_class))
    end do
          
  end subroutine update_collapsed_columns_classes

  !-----------------------------------------------------------------------
  function collapse_to_atm_topo_icemec_filterc(this, bounds) result(filter)
    !
    ! !DESCRIPTION:
    ! Returns a column-level filter of ice_mec columns with the collapse_to_atm_topo behavior
    !
    ! !USES:
    use filterColMod, only : filter_col_type, col_filter_from_grcflags_ltypes
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(in) :: this
    type(filter_col_type) :: filter  ! function result
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'collapse_to_atm_topo_icemec_filterc'
    !-----------------------------------------------------------------------

    ! Currently this creates the filter on the fly, recreating it every time this
    ! function is called. Since this is a static filter, we could just compute it once
    ! and save it, returning the already-computed filter when this function is called.
    ! However, the problem with that is the need to have a different filter for each
    ! clump (and potentially another filter for calls from outside a clump loop). This
    ! will become easier to handle if we rework CLM's threading so that there is a
    ! separate instance of each object for each clump: in that case, we'll have multiple
    ! instances of glc_behavior_type, each corresponding to one clump, each with its own
    ! filter.

    filter = col_filter_from_grcflags_ltypes( &
         bounds = bounds, &
         grcflags = this%collapse_to_atm_topo_grc(bounds%begg:bounds%endg), &
         ltypes = [istice_mec])

  end function collapse_to_atm_topo_icemec_filterc


  !-----------------------------------------------------------------------
  subroutine set_masks(this, begg, endg, glcmask, latdeg)
    !
    ! !DESCRIPTION:
    ! Set the mask fields in this object
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(inout) :: this
    integer, intent(in) :: begg  ! beginning gridcell index
    integer, intent(in) :: endg  ! ending gridcell index
    integer, intent(in) :: glcmask(begg: )
    real(r8), intent(in) :: latdeg(begg: )
    !
    ! !LOCAL VARIABLES:
    integer :: g
    logical, parameter :: collapse_mountain_glaciers = .true.

    character(len=*), parameter :: subname = 'set_masks'
    !-----------------------------------------------------------------------

    do g = begg, endg
       if (glcmask(g) > 0) then
          this%has_virtual_columns_grc(g) = .true.
       else
          this%has_virtual_columns_grc(g) = .false.
       end if
    end do

    if (collapse_mountain_glaciers) then
       do g = begg, endg
          if (this%has_virtual_columns_grc(g)) then
             ! Greenland
             this%collapse_to_atm_topo_grc(g) = .false.
          else if (latdeg(g) <= -60._r8) then
             ! Antarctica
             this%collapse_to_atm_topo_grc(g) = .false.
          else
             ! mountain glaciers
             this%collapse_to_atm_topo_grc(g) = .true.
          end if
       end do
    else
       this%collapse_to_atm_topo_grc(begg:endg) = .false.
    end if

  end subroutine set_masks

end module glcBehaviorMod
