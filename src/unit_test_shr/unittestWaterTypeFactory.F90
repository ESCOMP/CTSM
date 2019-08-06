module unittestWaterTypeFactory

  ! This module contains a class and associated methods that assist unit tests that need
  ! to create a water_type instance
  !
  ! This is object-oriented in case we want to remember initial values for the sake of
  ! doing a clean teardown, or have any other memory between calls (currently we don't
  ! make use of this, but we could in the future).
  !
  ! Typical usage is:
  !
  !    - Include an instance of this class in the class used for unit testing
  !
  !    - Before calling any other methods: call init()
  !
  !    - Before doing any subgrid setup: call setup_before_subgrid()
  !
  !    - After doing any subgrid setup: call setup_after_subgrid()
  !
  !    - Then get the water_type inst with: call create_water_type()
  !
  !    - In the unit test tearDown method, before unittest_subgrid_teardown: call teardown()

#include "shr_assert.h"
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use clm_varpar, only : nlevsoi, nlevgrnd, nlevsno
  use ColumnType, only : col
  use WaterType, only : water_type, water_params_type
  use unittestArrayMod, only : col_array
  use unittestSubgridMod, only : bounds

  implicit none
  private

  type, public :: unittest_water_type_factory_type
   contains
     procedure, public :: init
     procedure, public :: setup_before_subgrid
     procedure, public :: setup_after_subgrid
     procedure, public :: create_water_type
     procedure, public :: teardown
  end type unittest_water_type_factory_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  subroutine init(this)
    ! Initialize this instance (like a constructor, but Fortran constructors can be
    ! problematic / buggy, so we do it via an init subroutine instead)
    class(unittest_water_type_factory_type), intent(inout) :: this

    ! For now, nothing to do: just a placeholder
  end subroutine init

  subroutine setup_before_subgrid(this, my_nlevsoi, nlevgrnd_additional, my_nlevsno)
    ! Do initializations that need to happen before setting up the subgrid structure.
    !
    ! If you have already initialized any of these variables yourself, then simply pass
    ! in the current values to make this a routine a no-op.
    class(unittest_water_type_factory_type), intent(inout) :: this
    integer, intent(in) :: my_nlevsoi
    integer, intent(in) :: nlevgrnd_additional  ! nlevgrnd = nlevsoi + nlevgrnd_additional
    integer, intent(in) :: my_nlevsno

    nlevsoi = my_nlevsoi
    nlevgrnd = nlevsoi + nlevgrnd_additional
    nlevsno = my_nlevsno
  end subroutine setup_before_subgrid

  subroutine setup_after_subgrid(this, snl, dz)
    ! Do initializations that need to happen after setting up the subgrid structure
    class(unittest_water_type_factory_type), intent(inout) :: this

    ! For now, set all snl and dz values to the single input
    ! If either of these are absent, then we assume that the relevant col values have
    ! already been set, and we do nothing here. (e.g., if dz is absent, we assume that
    ! col%dz has already been set, and don't do anything with col%dz here.)
    integer, intent(in), optional :: snl
    real(r8), intent(in), optional :: dz

    if (present(snl)) then
       col%snl(:) = snl
    end if
    if (present(dz)) then
       col%dz(:,:) = dz
    end if
  end subroutine setup_after_subgrid

  subroutine create_water_type(this, water_inst, &
       t_soisno_col, watsat_col, &
       enable_consistency_checks, enable_isotopes)
    ! Initialize water_inst
    !
    ! Assumes that setup_before_subgrid and setup_after_subgrid have been called and that
    ! subgrid setup has been done.
    !
    ! Arbitrary values are used for some of the inputs to the water_inst init routine
    !
    ! If enable_consistency_checks or enable_isotopes are missing, they are assumed to be
    ! false.
    class(unittest_water_type_factory_type), intent(in) :: this
    type(water_type), intent(inout) :: water_inst

    ! Temperature of each soil and snow layer, for each column and layer. Second
    ! dimension should go from -nlevsno+1..nlevgrnd. If not provided, we use arbitrary
    ! values for this variable.
    real(r8), intent(in), optional :: t_soisno_col( bounds%begc: , -nlevsno+1: )

    ! watsat for each soil layer, for each column. Second dimension should go from
    ! 1..nlevgrnd. If not provided, we use arbitrary values for this variable.
    real(r8), intent(in), optional :: watsat_col( bounds%begc: , 1: )

    logical, intent(in), optional :: enable_consistency_checks
    logical, intent(in), optional :: enable_isotopes

    logical :: l_enable_consistency_checks
    logical :: l_enable_isotopes
    type(water_params_type) :: params
    real(r8) :: l_watsat_col(bounds%begc:bounds%endc, nlevgrnd)
    real(r8) :: l_t_soisno_col(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd)

    if (present(t_soisno_col)) then
       SHR_ASSERT_ALL((ubound(t_soisno_col) == [bounds%endc, nlevgrnd]), errMsg(sourcefile, __LINE__))
    end if
    if (present(watsat_col)) then
       SHR_ASSERT_ALL((ubound(watsat_col) == [bounds%endc, nlevgrnd]), errMsg(sourcefile, __LINE__))
    end if

    if (present(enable_consistency_checks)) then
       l_enable_consistency_checks = enable_consistency_checks
    else
       l_enable_consistency_checks = .false.
    end if

    if (present(enable_isotopes)) then
       l_enable_isotopes = enable_isotopes
    else
       l_enable_isotopes = .false.
    end if

    params = water_params_type( &
         enable_consistency_checks = l_enable_consistency_checks, &
         enable_isotopes = l_enable_isotopes)

    if (present(watsat_col)) then
       l_watsat_col(:,:) = watsat_col(:,:)
    else
       l_watsat_col(:,:) = 0.2_r8
    end if

    if (present(t_soisno_col)) then
       l_t_soisno_col(:,:) = t_soisno_col(:,:)
    else
       l_t_soisno_col(:,:) = 275._r8
    end if

    call water_inst%InitForTesting(bounds, params, &
         h2osno_col = col_array(0._r8), &
         snow_depth_col = col_array(0._r8), &
         watsat_col = l_watsat_col, &
         t_soisno_col = l_t_soisno_col)
  end subroutine create_water_type

  subroutine teardown(this, water_inst)
    ! Should be called from the unittest tearDown method, before unittest_subgrid_teardown
    class(unittest_water_type_factory_type), intent(inout) :: this
    type(water_type), intent(in) :: water_inst

    ! For now, nothing to do: just a placeholder.
    !
    ! Ideally this would call water_inst%Clean. But we don't yet have a clean method in
    ! water_type or most of the types it contains. So for now we have a small memory leak
    ! in each test.
  end subroutine teardown

end module unittestWaterTypeFactory
