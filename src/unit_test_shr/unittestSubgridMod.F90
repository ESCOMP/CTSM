module unittestSubgridMod

  ! Provides routines to aid with the setup of subgrid structure for unit tests that need
  ! it. 
  !
  ! In the setup for a test, the following should be done:
  !
  ! (1) call unittest_subgrid_setup_start
  ! (2) add grid cells, landunits, columns & pfts as desired, using the routines defined in
  !     this module (i.e., using unittest_add_landunit, etc. - NOT directly via add_landunit, etc.)
  ! (3) call unittest_subgrid_setup_end
  !
  !   Example: To add a single grid cell, with two landunits (nat. veg. and icemec), with a
  !   single column on the nat veg landunit, the following can be done:
  !
  !     call unittest_subgrid_setup_start()
  !     call unittest_add_gridcell()
  !     call unittest_add_landunit(my_gi=gi, ltype=istsoil, wtgcell=0.4_r8)
  !     call unittest_add_column(my_li=li, ctype=1, wtlunit=1.0_r8)
  !     c_soil = ci
  !     call unittest_add_landunit(my_gi=gi, ltype=istice_mec, wtgcell=0.6_r8)
  !     call unittest_subgrid_setup_end()
  ! 
  !   A few things to note about this example:
  !   (1) Note the use of gi, li and ci to get the index of the most recently-added grid
  !       cell / landunit / column
  !   (2) Note that not all subgrid information has been filled in: no patches were added
  !       to the soil landunit, and no columns or patches were added to the icemec
  !       landunit. This is because this extra level of detail wasn't needed for this
  !       particular unit test. This omission is perfectly acceptable.
  ! 
  ! In the teardown for a test, the following should be done:
  ! 
  ! (1) call unittest_subgrid_teardown

  use shr_kind_mod , only : r8 => shr_kind_r8
  use decompMod    , only : bounds_type

  implicit none
  private
  save

  ! ------------------------------------------------------------------------
  ! Public entities
  ! ------------------------------------------------------------------------

  ! Public routines
  public :: unittest_subgrid_setup_start ! do the initial setup of subgrid stuff needed for unit testing
  public :: unittest_subgrid_setup_end   ! do the last part of setup
  public :: unittest_subgrid_teardown    ! do any teardown needed for the subgrid stuff
  public :: unittest_add_gridcell        ! add a grid cell
  public :: unittest_add_landunit        ! add a landunit
  public :: unittest_add_column          ! add a column
  public :: unittest_add_patch           ! add a patch

  ! bounds info, which can be passed to routines that need it
  ! Note that the end indices here (endg, endl, endc, endp) will be the final indices in
  ! use, in contrast to the module-level endg, endl, etc., which give the final indices
  ! of the allocated arrays.
  type(bounds_type), public, protected :: bounds
  
  ! Indices of last grid cell / landunit / column / patch added
  integer, public, protected :: gi
  integer, public, protected :: li
  integer, public, protected :: ci
  integer, public, protected :: pi

  ! Maximum array sizes at each level
  integer, parameter, public :: numg = 3
  integer, parameter, public :: numl = 30
  integer, parameter, public :: numc = 50
  integer, parameter, public :: nump = 100

  ! Indices of initial grid cell / landunit / column / patch
  !
  ! Note that we do NOT start at 1, in order to catch any code that assumes indices start
  ! at 1.
  integer, parameter, public :: begg = 11
  integer, parameter, public :: begl = 21
  integer, parameter, public :: begc = 31
  integer, parameter, public :: begp = 41

  ! Indices of final grid cell / landunit / column / patch
  ! Note that these are the final indices of the allocated arrays, which may be greater
  ! than the final index that is actually used for a given test.
  integer, parameter, public :: endg = begg + numg - 1
  integer, parameter, public :: endl = begl + numl - 1
  integer, parameter, public :: endc = begc + numc - 1
  integer, parameter, public :: endp = begp + nump - 1
  
  ! ------------------------------------------------------------------------
  ! Private entities
  ! ------------------------------------------------------------------------

  interface initialize_array
     module procedure initialize_array_int
     module procedure initialize_array_r8
     module procedure initialize_array_logical
  end interface initialize_array

contains
  
  !-----------------------------------------------------------------------
  subroutine unittest_subgrid_setup_start
    !
    ! !DESCRIPTION:
    ! Do the initial setup of subgrid stuff needed for unit testing. This should be
    ! called for each test.
    !
    ! !USES:
    use clm_varpar, only : natpft_lb
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'unittest_subgrid_setup_start'
    !-----------------------------------------------------------------------

    call initialize_arrays

    ! Initialize local module variables

    gi = begg - 1
    li = begl - 1
    ci = begc - 1
    pi = begp - 1
    
    ! Initialize other variables needed for the subgrid setup
    
    natpft_lb = 0
    
    
  end subroutine unittest_subgrid_setup_start

  !-----------------------------------------------------------------------
  subroutine unittest_subgrid_setup_end
    !
    ! !DESCRIPTION:
    ! Do the last part of setup. This should be called after adding all of the landunits,
    ! columns, pfts, etc. for the test.
    !
    ! !USES:
    use initSubgridMod, only : clm_ptrs_compdown
    use subgridWeightsMod, only : compute_higher_order_weights
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'unittest_subgrid_setup_end'
    !-----------------------------------------------------------------------
    
    call set_bounds
    call clm_ptrs_compdown(bounds)
    call compute_higher_order_weights(bounds)

  end subroutine unittest_subgrid_setup_end

  !-----------------------------------------------------------------------
  subroutine set_bounds
    !
    ! !DESCRIPTION:
    ! Create the bounds derived type object
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'set_bounds'
    !-----------------------------------------------------------------------
    
    bounds%begg = begg
    bounds%endg = gi
    bounds%begl = begl
    bounds%endl = li
    bounds%begc = begc
    bounds%endc = ci
    bounds%begp = begp
    bounds%endp = pi

    ! Currently, not setting bounds%level and bounds%clump_index
    
  end subroutine set_bounds



  !-----------------------------------------------------------------------
  subroutine initialize_arrays
    !
    ! !DESCRIPTION:
    ! Allocate arrays in clmtype, and initialize them to default values. Note that we only
    ! do allocation if arrays are not currently allocated; this allows us to avoid
    ! deallocating and reallocating arrays for every test - instead just reinitializing
    ! them.
    !
    ! In principle, we should be able to use the init routines in clmtypeInitMod
    ! (init_pft_type, init_column_type, etc.) instead of doing this all manually. But to
    ! do that, I would want to make a few changes, which I don't want to do until
    ! Mariana's big object-oriented refactor is on the trunk. Once it's acceptable to
    ! make these changes, I would want to:
    !
    ! (1) Change the init_*_type routines to set default values (ispval or nan) in all
    ! arrays
    !
    ! (2) Add destructor routines for each type
    !
    ! (3) Change the code here to just call init_pft_type, init_column_type, etc.
    !
    ! (4) Put some teardown code here that calls the destructor for pft, col, etc. Note
    ! that this will mean that we redo the allocation for each test, but I ran some
    ! timing tests that show this to be a minimal additional cost. 
    !
    ! !USES:
    use landunit_varcon, only : max_lunit
    use clm_varcon, only : ispval
    use clmtype, only : grc, lun, col, pft
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initialize_arrays'
    !-----------------------------------------------------------------------
    
    ! Do the same allocation-if-necessary for landunit_indices as is done for other
    ! arrays... but since it is the only 2-d array, we don't have a generic subroutine
    ! for it.
    if (.not. associated(grc%landunit_indices)) then
       allocate(grc%landunit_indices(1:max_lunit, begg:endg))
    end if
    grc%landunit_indices(1:max_lunit, begg:endg) = ispval

    call initialize_array(lun%gridcell, begl, endl)
    call initialize_array(lun%wtgcell, begl, endl)
    call initialize_array(lun%coli, begl, endl)
    call initialize_array(lun%colf, begl, endl)
    call initialize_array(lun%ncolumns, begl, endl)
    call initialize_array(lun%pfti, begl, endl)
    call initialize_array(lun%pftf, begl, endl)
    call initialize_array(lun%npfts, begl, endl)
    call initialize_array(lun%itype, begl, endl)
    call initialize_array(lun%ifspecial, begl, endl)
    call initialize_array(lun%lakpoi, begl, endl)
    call initialize_array(lun%urbpoi, begl, endl)
    call initialize_array(lun%glcmecpoi, begl, endl)
    call initialize_array(lun%active, begl, endl)

    call initialize_array(col%gridcell, begc, endc)
    call initialize_array(col%wtgcell, begc, endc)
    call initialize_array(col%landunit, begc, endc)
    call initialize_array(col%wtlunit, begc, endc)
    call initialize_array(col%pfti, begc, endc)
    call initialize_array(col%pftf, begc, endc)
    call initialize_array(col%npfts, begc, endc)
    call initialize_array(col%itype, begc, endc)
    call initialize_array(col%active, begc, endc)

    call initialize_array(pft%gridcell, begp, endp)
    call initialize_array(pft%wtgcell, begp, endp)
    call initialize_array(pft%landunit, begp, endp)
    call initialize_array(pft%wtlunit, begp, endp)
    call initialize_array(pft%column, begp, endp)
    call initialize_array(pft%wtcol, begp, endp)
    call initialize_array(pft%itype, begp, endp)
    call initialize_array(pft%mxy, begp, endp)
    call initialize_array(pft%active, begp, endp)


  end subroutine initialize_arrays


  !-----------------------------------------------------------------------
  subroutine initialize_array_int(arr, lb, ub)
    !
    ! !DESCRIPTION:
    ! Initialize an integer array. Allocates arr if necessary, and initializes it to a
    ! default value. Note that, if arr is already allocated, then lb and ub are ignored --
    ! arr is NOT resized.
    !
    ! !USES:
    use clm_varcon, only : ispval
    !
    ! !ARGUMENTS:
    integer , intent(inout), pointer :: arr(:) ! array to initialize
    integer , intent(in)             :: lb     ! lower bound of array, if arr needs allocation
    integer , intent(in)             :: ub     ! upper bound of array, if arr needs allocation
                                               !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initialize_array_int'
    !-----------------------------------------------------------------------
    
    if (.not. associated(arr)) then
       allocate(arr(lb:ub))
    end if

    arr(lb:ub) = ispval

  end subroutine initialize_array_int

  !-----------------------------------------------------------------------
  subroutine initialize_array_r8(arr, lb, ub)
    !
    ! !DESCRIPTION:
    ! Initialize an r8 array. Allocates arr if necessary, and initializes it to a default
    ! value. Note that, if arr is already allocated, then lb and ub are ignored -- arr is
    ! NOT resized.
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    real(r8) , intent(inout), pointer :: arr(:) ! array to initialize
    integer  , intent(in)             :: lb     ! lower bound of array, if arr needs allocation
    integer  , intent(in)             :: ub     ! upper bound of array, if arr needs allocation
                                                !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initialize_array_r8'
    !-----------------------------------------------------------------------

    if (.not. associated(arr)) then
       allocate(arr(lb:ub))
    end if

    arr(lb:ub) = nan

  end subroutine initialize_array_r8

  !-----------------------------------------------------------------------
  subroutine initialize_array_logical(arr, lb, ub)
    !
    ! !DESCRIPTION:
    ! Initialize a logical array. Allocates arr if necessary, and initializes it to a default
    ! value. Note that, if arr is already allocated, then lb and ub are ignored -- arr is
    ! NOT resized.
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    logical  , intent(inout), pointer :: arr(:) ! array to initialize
    integer  , intent(in)             :: lb     ! lower bound of array, if arr needs allocation
    integer  , intent(in)             :: ub     ! upper bound of array, if arr needs allocation
                                                !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initialize_array_logical'
    !-----------------------------------------------------------------------

    if (.not. associated(arr)) then
       allocate(arr(lb:ub))
    end if

    arr(lb:ub) = .false.

  end subroutine initialize_array_logical


  !-----------------------------------------------------------------------
  subroutine unittest_subgrid_teardown
    !
    ! !DESCRIPTION:
    ! Do any teardown needed for the subgrid stuff
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'unittest_subgrid_teardown'
    !-----------------------------------------------------------------------
    
    ! For now, nothing is needed... we currently don't bother with deallocation, and the
    ! initialization of arrays is done in the setup routine

  end subroutine unittest_subgrid_teardown

  !-----------------------------------------------------------------------
  subroutine unittest_add_gridcell()
    !
    ! !DESCRIPTION:
    ! Add a grid cell. The index of the just-added grid cell can be obtained from the
    ! module-level variable, gi.
    !
    ! Unlike add_landunit, add_column and add_patch, this is specific to the unit test
    ! code, because no such routine is needed in the production code
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'unittest_add_gridcell'
    !-----------------------------------------------------------------------
    
    gi = gi + 1

  end subroutine unittest_add_gridcell

  !-----------------------------------------------------------------------
  subroutine unittest_add_landunit(my_gi, ltype, wtgcell)
    !
    ! !DESCRIPTION:
    ! Add a landunit. The index of the just-added landunit can be obtained from the
    ! module-level variable, li.
    !
    ! This is simply a wrapper to the routine in initSubgridMod. We provide this for two
    ! reasons:
    !
    ! (1) To allow the module-level li variable to be protected
    !
    ! (2) To insulate most of the unit test code from any changes in the interface to
    ! add_landunit
    !
    ! !USES:
    use initSubgridMod, only : add_landunit
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: my_gi   ! grid cell index on which this landunit should be placed
    integer  , intent(in)    :: ltype   ! landunit type
    real(r8) , intent(in)    :: wtgcell ! weight of the landunit relative to the grid cell
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'unittest_add_landunit'
    !-----------------------------------------------------------------------

    call add_landunit(li=li, gi=my_gi, ltype=ltype, wtgcell=wtgcell)
    
  end subroutine unittest_add_landunit

  !-----------------------------------------------------------------------
  subroutine unittest_add_column(my_li, ctype, wtlunit)
    !
    ! !DESCRIPTION:
    ! Add a column. The index of the just-added column can be obtained from the
    ! module-level variable, ci.
    !
    ! This is simply a wrapper to the routine in initSubgridMod. We provide this for two
    ! reasons:
    !
    ! (1) To allow the module-level ci variable to be protected
    !
    ! (2) To insulate most of the unit test code from any changes in the interface to
    ! add_column
    !
    ! !USES:
    use initSubgridMod, only : add_column
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: my_li   ! landunit index on which this column should be placed
    integer  , intent(in)    :: ctype   ! column type
    real(r8) , intent(in)    :: wtlunit ! weight of the column relative to the land unit
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'unittest_add_column'
    !-----------------------------------------------------------------------

    call add_column(ci=ci, li=my_li, ctype=ctype, wtlunit=wtlunit)
    
  end subroutine unittest_add_column

  !-----------------------------------------------------------------------
  subroutine unittest_add_patch(my_ci, ptype, wtcol)
    !
    ! !DESCRIPTION:
    ! Add a patch. The index of the just-added patch can be obtained from the
    ! module-level variable, pi.
    !
    ! This is simply a wrapper to the routine in initSubgridMod. We provide this for two
    ! reasons:
    !
    ! (1) To allow the module-level pi variable to be protected
    !
    ! (2) To insulate most of the unit test code from any changes in the interface to
    ! add_patch
    !
    ! !USES:
    use initSubgridMod, only : add_patch
    !
    ! !ARGUMENTS:
    integer  , intent(in)    :: my_ci   ! column index on which this patch should be placed
    integer  , intent(in)    :: ptype   ! patch type
    real(r8) , intent(in)    :: wtcol   ! weight of the patch relative to the column
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'unittest_add_patch'
    !-----------------------------------------------------------------------

    call add_patch(pi=pi, ci=my_ci, ptype=ptype, wtcol=wtcol)
    
  end subroutine unittest_add_patch


end module unittestSubgridMod
