module unittestSimpleSubgridSetupsMod

  ! This module provides wrappers to unittestSubgridMod, which give you a variety of
  ! simple subgrid setups.
#include "shr_assert.h"
  use unittestSubgridMod
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod , only : r8 => shr_kind_r8
  use landunit_varcon, only : istsoil
  use pftconMod, only : noveg

  implicit none
  private
  save

  ! ------------------------------------------------------------------------
  ! Routines that do everything needed with the subgrid setup, including the begin & end
  ! call. Once you call these routines, you cannot add any more gridcells, landunits, etc.
  ! ------------------------------------------------------------------------

  ! Create a grid that has a single gridcell with a single vegetated patch
  public :: setup_single_veg_patch

  ! Create a grid that has a single gridcell with N vegetated patches on one column
  public :: setup_n_veg_patches

  ! Create a grid that has a single gridcell with one landunit of a given type with N
  ! columns, each with a single patch of a given type (default patch type = noveg)
  public :: setup_landunit_ncols

  ! Create a grid that has N grid cells, each with a single vegetated patch
  public :: setup_ncells_single_veg_patch

  ! ------------------------------------------------------------------------
  ! Routines that create a single grid cell with certain properties. You can do other
  ! subgrid setup (creating other grid cells) before and after this.
  ! ------------------------------------------------------------------------

  ! Create a grid cell that is 100% natural veg, with a single patch
  public :: create_gridcell_single_veg_patch

  ! ------------------------------------------------------------------------
  ! Routines that create a single landunit with certain properties. You can do other
  ! subgrid setup (creating other landunits, gridcells, etc.) before and after this.
  ! These assume that unittest_add_gridcell has already been called.
  ! ------------------------------------------------------------------------

  ! Create a landunit of a given type with N columns, each with M patches (default=1) of
  ! a given type (default=noveg)
  public :: create_landunit_ncols

  ! Create a vegetated landunit with N patches
  public :: create_vegetated_landunit_n_patches

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Routines that do everything needed with the subgrid setup, including the begin & end
  ! call. Once you call these routines, you cannot add any more gridcells, landunits, etc.
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine setup_single_veg_patch(pft_type)
    !
    ! !DESCRIPTION:
    ! Create a grid that has a single gridcell with a single vegetated patch, with veg
    ! type given by the pft_type argument
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: pft_type  ! the type of the single vegetated patch
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'setup_single_veg_patch'
    !-----------------------------------------------------------------------
    
    call setup_ncells_single_veg_patch(ncells=1, pft_type=pft_type)

  end subroutine setup_single_veg_patch

  !-----------------------------------------------------------------------
  subroutine setup_n_veg_patches(pwtcol, pft_types)
    !
    ! !DESCRIPTION:
    ! Create a grid that has a single gridcell with N vegetated patches on one column.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: pwtcol(:)

    ! If given, this gives the pft type for each patch. If not given, pft types go 1..N.
    integer, optional, intent(in) :: pft_types(:)

    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'setup_n_veg_patches'
    !-----------------------------------------------------------------------

    call unittest_subgrid_setup_start()
    call unittest_add_gridcell()
    call create_vegetated_landunit_n_patches(lweight = 1._r8, &
         pwtcol = pwtcol, pft_types = pft_types)
    call unittest_subgrid_setup_end()

  end subroutine setup_n_veg_patches

  !-----------------------------------------------------------------------
  subroutine setup_landunit_ncols(ltype, ctypes, cweights, ptype)
    !
    ! !DESCRIPTION:
    ! Create a grid that has a single gridcell with one landunit of a given type with N
    ! columns, each with a single patch of type ptype (or noveg if ptype is not given)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: ltype ! landunit type
    integer, intent(in) :: ctypes(:)  ! array of column types; one column is created for each element in the array
    real(r8), intent(in) :: cweights(:) ! array of column weights on the landunit
    integer, intent(in), optional :: ptype ! patch type (if not given, defaults to noveg)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'setup_landunit_ncols'
    !-----------------------------------------------------------------------

    call unittest_subgrid_setup_start()
    call unittest_add_gridcell()
    call create_landunit_ncols(ltype = ltype, lweight = 1._r8, &
         ctypes = ctypes, cweights = cweights, ptype = ptype)
    call unittest_subgrid_setup_end()

  end subroutine setup_landunit_ncols


  !-----------------------------------------------------------------------
  subroutine setup_ncells_single_veg_patch(ncells, pft_type)
    !
    ! !DESCRIPTION:
    ! Create a grid that has ncells grid cells, each with a single vegetated patch. All
    ! vegetated patches have the same type, given by pft_type.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: ncells    ! number of grid cells
    integer, intent(in) :: pft_type  ! pft type
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'setup_ncells_single_veg_patch'
    !-----------------------------------------------------------------------

    call unittest_subgrid_setup_start()
    do i = 1, ncells
       call create_gridcell_single_veg_patch(pft_type = pft_type)
    end do
    call unittest_subgrid_setup_end()

  end subroutine setup_ncells_single_veg_patch


  ! ========================================================================
  ! Routines that create a single grid cell with certain properties. You can do other
  ! subgrid setup (creating other grid cells) before and after this.
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine create_gridcell_single_veg_patch(pft_type)
    !
    ! !DESCRIPTION:
    ! Create a grid cell that is 100% natural veg, with a single patch
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: pft_type  ! the type of the single vegetated patch
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'create_gridcell_single_veg_patch'
    !-----------------------------------------------------------------------

    call unittest_add_gridcell()
    call unittest_add_landunit(my_gi=gi, ltype=istsoil, wtgcell=1.0_r8)
    call unittest_add_column(my_li=li, ctype=1, wtlunit=1.0_r8)
    call unittest_add_patch(my_ci=ci, ptype=pft_type, wtcol=1.0_r8)

  end subroutine create_gridcell_single_veg_patch

  ! ========================================================================
  ! Routines that create a single landunit with certain properties. You can do other
  ! subgrid setup (creating other landunits, gridcells, etc.) before and after this.
  ! These assume that unittest_add_gridcell has already been called.
  ! ========================================================================


  !-----------------------------------------------------------------------
  subroutine create_landunit_ncols(ltype, lweight, ctypes, cweights, npatches, ptype)
    !
    ! !DESCRIPTION:
    ! Create a landunit of a given type with N columns, each with M patches (default=1) of
    ! type noveg.
    !
    ! Assumes that unittest_add_gridcell has already been called.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: ltype ! landunit type
    real(r8), intent(in) :: lweight ! landunit weight on the grid cell
    integer, intent(in) :: ctypes(:)  ! array of column types; one column is created for each element in the array
    real(r8), intent(in) :: cweights(:) ! array of column weights on the landunit

    ! If npatches is provided, it gives the number of patches on each column. Each patch
    ! has equal weight, and is of the same type (ptype). If not provided, the default is 1
    ! patch per column.
    integer, intent(in), optional :: npatches

    ! If ptype is provided, it gives the pft type of every patch. If not provided, the
    ! default is noveg.
    integer, intent(in), optional :: ptype
    !
    ! !LOCAL VARIABLES:
    integer :: ncols
    integer :: l_npatches ! local version of npatches
    integer :: l_ptype    ! local version of ptype
    integer :: c, p

    character(len=*), parameter :: subname = 'create_landunit_ncols'
    !-----------------------------------------------------------------------

    ncols = size(ctypes)
    SHR_ASSERT((size(cweights) == ncols), errMsg(sourcefile, __LINE__))
    SHR_ASSERT(gi >= begg, 'must call unittest_add_gridcell first: ' // errMsg(sourcefile, __LINE__))

    if (present(npatches)) then
       l_npatches = npatches
    else
       l_npatches = 1
    end if

    if (present(ptype)) then
       l_ptype = ptype
    else
       l_ptype = noveg
    end if

    call unittest_add_landunit(my_gi=gi, ltype=ltype, wtgcell=lweight)
    do c = 1, ncols
       call unittest_add_column(my_li=li, ctype=ctypes(c), wtlunit=cweights(c))
       do p = 1, l_npatches
          call unittest_add_patch(my_ci=ci, ptype=l_ptype, wtcol=1.0_r8/l_npatches)
       end do
    end do

  end subroutine create_landunit_ncols

  !-----------------------------------------------------------------------
  subroutine create_vegetated_landunit_n_patches(lweight, pwtcol, pft_types)
    !
    ! !DESCRIPTION:
    ! Create a vegetated landunit with N patches on one column
    !
    ! Assumes that unittest_add_gridcell has already been called
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: lweight ! landunit weight on the gridcell
    real(r8), intent(in) :: pwtcol(:) ! patch weights on the column

    ! If given, this gives the pft type for each patch. If not given, pft types go 1..N.
    integer, optional, intent(in) :: pft_types(:)
    !
    ! !LOCAL VARIABLES:
    integer :: npatches
    integer :: p
    integer, allocatable :: l_pft_types(:)

    character(len=*), parameter :: subname = 'create_vegetated_landunit_n_patches'
    !-----------------------------------------------------------------------

    npatches = size(pwtcol)
    allocate(l_pft_types(npatches))
    if (present(pft_types)) then
       SHR_ASSERT((size(pft_types) == npatches), errMsg(sourcefile, __LINE__))
       l_pft_types = pft_types
    else
       do p = 1, npatches
          l_pft_types(p) = p
       end do
    end if

    call unittest_add_landunit(my_gi=gi, ltype=istsoil, wtgcell=lweight)
    call unittest_add_column(my_li=li, ctype=1, wtlunit=1._r8)
    do p = 1, npatches
       call unittest_add_patch(my_ci=ci, ptype=l_pft_types(p), wtcol=pwtcol(p))
    end do

  end subroutine create_vegetated_landunit_n_patches



end module unittestSimpleSubgridSetupsMod
