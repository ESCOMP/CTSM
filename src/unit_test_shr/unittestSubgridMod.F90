module unittestSubgridMod

  ! Provides routines to aid with the setup of subgrid structure for unit tests that need
  ! it. 
  !
  ! In the setup for a test, the following should be done:
  !
  ! (1) call unittest_subgrid_setup_start
  !     Note: if explicitly setting nlevsno, that must be done *before* the call to
  !     unittest_subgrid_setup_start
  ! (2) add grid cells, landunits, columns & pfts as desired, using the routines defined in
  !     this module (i.e., using unittest_add_landunit, etc. - NOT directly via add_landunit, etc.)
  ! (3) call unittest_subgrid_setup_end
  !
  !   Example: To add a single grid cell, with two landunits (nat. veg. and ice), with a
  !   single column on the nat veg landunit, the following can be done:
  !
  !     call unittest_subgrid_setup_start()
  !     call unittest_add_gridcell()
  !     call unittest_add_landunit(my_gi=gi, ltype=istsoil, wtgcell=0.4_r8)
  !     call unittest_add_column(my_li=li, ctype=1, wtlunit=1.0_r8)
  !     c_soil = ci
  !     call unittest_add_landunit(my_gi=gi, ltype=istice, wtgcell=0.6_r8)
  !     call unittest_subgrid_setup_end()
  ! 
  !   A few things to note about this example:
  !   (1) Note the use of gi, li and ci to get the index of the most recently-added grid
  !       cell / landunit / column
  !   (2) Note that not all subgrid information has been filled in: no patches were added
  !       to the soil landunit, and no columns or patches were added to the ice
  !       landunit. This is because this extra level of detail wasn't needed for this
  !       particular unit test. This omission is perfectly acceptable.
  ! 
  ! In the teardown for a test, the following should be done:
  ! 
  ! (1) call unittest_subgrid_teardown
  !
  !    Note: This can safely be done even if subgrid stuff was never set up for some
  !    tests (it will only do anything if the subgrid setup was done).

  use shr_kind_mod , only : r8 => shr_kind_r8
  use decompMod    , only : bounds_type, procinfo, get_proc_bounds
  use decompMod    , only : gindex_grc, gindex_lun, gindex_col, gindex_patch
  use GridcellType , only : grc                
  use LandunitType , only : lun                
  use ColumnType   , only : col                
  use PatchType    , only : patch                

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
  public :: get_ltype_special            ! get a landunit type corresponding to a special landunit

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
  integer, parameter, public :: numg = 6
  integer, parameter, public :: numl = 30
  integer, parameter, public :: numc = 50
  integer, parameter, public :: nump = 100

  ! Indices of initial grid cell / landunit / column / patch
  !
  ! Now we do start at 1.
  integer, parameter, public :: begg = 1
  integer, parameter, public :: begl = 1
  integer, parameter, public :: begc = 1
  integer, parameter, public :: begp = 1

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

  integer, private :: nlevsno_orig ! original value of nlevsno, saved so we can restore it later
  logical, private :: nlevsno_set  ! whether we set nlevsno here
  logical, private :: unittest_subgrid_needs_teardown = .false. ! whether subgrid stuff has been initialized

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

    unittest_subgrid_needs_teardown = .true.
    
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

    call set_decomp_info
    call create_bounds_object
    call clm_ptrs_compdown(bounds)
    call compute_higher_order_weights(bounds)

  end subroutine unittest_subgrid_setup_end

  !-----------------------------------------------------------------------
  subroutine set_decomp_info
    !
    ! !DESCRIPTION:
    ! Set up decomp info in decompMod.
    !
    ! We need to do this (in addition to just making sure that the bounds derived type
    ! object is set up correctly) for the sake of callers of get_proc_bounds.
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'set_decomp_info'
    !-----------------------------------------------------------------------

    ! For now, not setting up clump info, because it isn't needed in any unit tests. We
    ! may have to fix this in the future.
    procinfo%nclumps = 1
    allocate(procinfo%cid(procinfo%nclumps))
    procinfo%cid(:) = -1

    procinfo%begg = begg
    procinfo%endg = gi
    procinfo%begl = begl
    procinfo%endl = li
    procinfo%begc = begc
    procinfo%endc = ci
    procinfo%begp = begp
    procinfo%endp = pi

    procinfo%ncells = procinfo%endg - procinfo%begg + 1
    procinfo%nlunits = procinfo%endl - procinfo%begl + 1
    procinfo%ncols = procinfo%endc - procinfo%begc + 1
    procinfo%npatches = procinfo%endp - procinfo%begp + 1

    ! Currently leaving cohort info unset because it isn't needed in any unit tests. We
    ! may have to fix this in the future.


    ! The following are needed in endrun calls, but since they are only used for printing
    ! information, we can fill them with garbage.
    allocate(gindex_grc(procinfo%begg:procinfo%endg))
    allocate(gindex_lun(procinfo%begl:procinfo%endl))
    allocate(gindex_col(procinfo%begc:procinfo%endc))
    allocate(gindex_patch(procinfo%begp:procinfo%endp))
    gindex_grc(:) = 0
    gindex_lun(:) = 0
    gindex_col(:) = 0
    gindex_patch(:) = 0

  end subroutine set_decomp_info

  !-----------------------------------------------------------------------
  subroutine create_bounds_object
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

    ! Some routines want a proc-level bounds. So for now, just making bounds be
    ! proc-level. In the future, we may need both a proc-level and clumps-level bounds
    ! object (if other routines want a clump-level bounds). (For the sake of unit
    ! testing, proc-level and clump-level bounds objects can probably be the same except
    ! for bounds%level and bounds%clump_index.)
    call get_proc_bounds(bounds)

  end subroutine create_bounds_object



  !-----------------------------------------------------------------------
  subroutine initialize_arrays
    !
    ! !DESCRIPTION:
    ! Allocate subgrid arrays, and initialize them to default values.
    !
    ! !USES:
    use landunit_varcon , only : max_lunit
    use clm_varcon      , only : ispval
    use GridcellType    , only : grc
    use LandunitType    , only : lun
    use ColumnType      , only : col
    use PatchType       , only : patch
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initialize_arrays'
    !-----------------------------------------------------------------------

    ! column initialization depends on the nlevsno runtime parameter, so we first need to
    ! set that
    call init_nlevsno()

    call grc%Init(begg, endg)
    call lun%Init(begl, endl)
    call col%Init(begc, endc)
    call patch%init(begp, endp)

  end subroutine initialize_arrays

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

    if (unittest_subgrid_needs_teardown) then
       call grc%clean
       call lun%clean
       call col%clean
       call patch%clean

       call reset_nlevsno()

       unittest_subgrid_needs_teardown = .false.
    end if

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
    ! Add a landunit, and make it active. The index of the just-added landunit can be
    ! obtained from the module-level variable, li.
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
    lun%active(li) = .true.
    
  end subroutine unittest_add_landunit

  !-----------------------------------------------------------------------
  subroutine unittest_add_column(my_li, ctype, wtlunit)
    !
    ! !DESCRIPTION:
    ! Add a column, and make it active. The index of the just-added column can be obtained
    ! from the module-level variable, ci.
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
    col%active(ci) = .true.
    
  end subroutine unittest_add_column

  !-----------------------------------------------------------------------
  subroutine unittest_add_patch(my_ci, ptype, wtcol)
    !
    ! !DESCRIPTION:
    ! Add a patch, and make it active. The index of the just-added patch can be obtained
    ! from the module-level variable, pi.
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
    patch%active(pi) = .true.

  end subroutine unittest_add_patch

  !-----------------------------------------------------------------------
  function get_ltype_special() result(ltype)
    !
    ! !DESCRIPTION:
    ! Returns a landunit type corresponding to a special landunit
    !
    ! !USES:
    use landunit_varcon, only : max_lunit, landunit_is_special
    !
    ! !ARGUMENTS:
    integer :: ltype  ! function result
    !
    ! !LOCAL VARIABLES:
    integer :: ltype_test
    logical :: found

    character(len=*), parameter :: subname = 'get_ltype_special'
    !-----------------------------------------------------------------------

    found = .false.
    ltype_test = 1
    do while (ltype_test <= max_lunit .and. .not. found)
       if (landunit_is_special(ltype_test)) then
          ltype = ltype_test
          found = .true.
       else
          ltype_test = ltype_test + 1
       end if
    end do

    if (.not. found) then
       ! We use a 'stop 1' here instead of shr_sys_abort because, in the context of pFUnit
       ! testing, a shr_sys_abort will allow the code to continue (after raising a pFUnit
       ! exception), which can lead to a cryptic error due to the return value from this
       ! function being invalid, giving an array out of bounds error.
       print *, subname//' ERROR: cannot find a special landunit'
       stop 1
    end if

  end function get_ltype_special


  subroutine init_nlevsno()
    ! Initialize nlevsno to a reasonable value, if it is not already set

    use clm_varpar, only : nlevsno

    if (nlevsno <= 0) then
       nlevsno_orig = nlevsno
       nlevsno = 5
       nlevsno_set = .true.
    else
       nlevsno_set = .false.
    end if
  end subroutine init_nlevsno

  subroutine reset_nlevsno
    ! If we set nlevsno in init_nlevsno, then reset it to its original value

    use clm_varpar, only : nlevsno

    if (nlevsno_set) then
       nlevsno = nlevsno_orig
    end if
  end subroutine reset_nlevsno

end module unittestSubgridMod
