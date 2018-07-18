module unittestGlcMec

  ! This module contains routines that assist unit tests working with glc_mec
  ! (istice_mec) landunits.

  use shr_kind_mod , only : r8 => shr_kind_r8
  use unittestSubgridMod
  use unittestSimpleSubgridSetupsMod
  use landunit_varcon, only : istice_mec
  use column_varcon, only : icemec_class_to_col_itype
  use glc_elevclass_mod, only : glc_elevclass_init, glc_elevclass_clean
  use clm_varpar, only : maxpatch_glcmec

  implicit none
  private
  save

  ! Sets up modules to have the given elevation classes. This should be called early in
  ! the setup of each unit test.
  public :: setup_elevation_classes

  ! Do teardown mirroring what was done in setup_elevation_classes. This should be called
  ! at the end of any unit test that called setup_elevation_classes.
  public :: teardown_elevation_classes

  ! Do all subgrid setup needed for setting up a grid with one grid cell with a single
  ! icemec column.
  public :: setup_single_icemec_column

contains

  subroutine setup_elevation_classes(glc_nec, topomax)
    ! Sets up modules to have the given elevation classes.
    !
    ! This should be called early in the setup of each unit test.
    integer, intent(in) :: glc_nec  ! number of elevation classes
    real(r8), intent(in) :: topomax(:)  ! should be size glc_nec+1

    call glc_elevclass_init(glc_nec, topomax)
    maxpatch_glcmec = glc_nec
  end subroutine setup_elevation_classes

  subroutine teardown_elevation_classes()
    ! Do teardown mirroring what was done in setup_elevation_classes.

    ! This should be called at the end of any unit test that called
    ! setup_elevation_classes.

    call glc_elevclass_clean()
    maxpatch_glcmec = 0
  end subroutine teardown_elevation_classes

  subroutine setup_single_icemec_column(elev_class)
    ! Create a grid cell with a single icemec column.
    !
    ! setup_elevation_classes must already have been called.
    integer, intent(in) :: elev_class

    call unittest_subgrid_setup_start()
    call unittest_add_gridcell()
    call create_landunit_ncols(ltype=istice_mec, lweight=1.0_r8, &
         ctypes=[icemec_class_to_col_itype(elev_class)], cweights=[1.0_r8])
    call unittest_subgrid_setup_end()

  end subroutine setup_single_icemec_column


end module unittestGlcMec
