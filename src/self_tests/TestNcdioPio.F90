module TestNcdioPio

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains tests of ncdio_pio

#include "shr_assert.h"
  use ncdio_pio
  use shr_kind_mod, only : r8 => shr_kind_r8
  use Assertions, only : assert_equal
  use clm_varcon, only : nameg
  use abortutils, only : endrun
  use array_utils, only : convert_to_logical
  use decompMod, only : bounds_type, get_proc_global
  use GridcellType, only : grc
  use spmdMod, only : masterproc

  implicit none
  private
  save

  ! Public routines

  public :: test_ncdio_pio

  ! Module data used in various tests

  integer, parameter :: nlev1 = 5
  character(len=*), parameter :: lev1_name = 'lev1'

  integer, parameter :: nlev2 = 3
  character(len=*), parameter :: lev2_name = 'lev2'

  character(len=*), parameter :: testfilename = 'test_ncdio_pio.nc'

  real(r8), pointer :: data_double_1d_grc(:)
  integer, pointer :: data_int_1d_grc(:)
  integer, pointer :: data_intlogical_1d_grc(:)  ! only has 0 & 1 values (representing a logical)
  real(r8), pointer :: data_double_2d_grc(:,:)
  real(r8), pointer :: data_double_3d_grc(:,:,:)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine test_ncdio_pio(bounds)
    !
    ! !DESCRIPTION:
    ! Drive tests of ncdio_pio
    !
    ! NOTE(wjs, 2020-10-15) Currently, endrun is called when any test assertion fails. I
    ! thought about changing this so that, instead, a counter is incremented for each
    ! failure, then at the end of the testing (in the higher-level self-test driver),
    ! endrun is called if this counter is greater than 0. The benefit of this is that we'd
    ! see all test failures, not just the first failure. To do that, we'd need to change
    ! the assertions here to increment a counter rather than aborting. However, I'm not
    ! spending the time to make this change for now because (1) I'm not sure how much
    ! value we'd get from it; (2) even if we made that change, it's still very possible
    ! for test code to abort for reasons other than assertions, if something goes wrong
    ! inside ncdio_pio or pio; and (3) some tests here are dependent on earlier tests (for
    ! example, the reads depend on the writes having worked), so a failure in an early
    ! phase could really muck things up for later testing phases. Migrating to a
    ! pFUnit-based unit test would solve this problem, since each pFUnit test is
    ! independent, though would prevent us from being able to have dependent tests the
    ! way we do here (where reads depend on earlier writes), for better or for worse.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'test_ncdio_pio'
    !-----------------------------------------------------------------------

    call write_to_log(subname//': create_vars')
    call create_vars(bounds)

    call write_to_log(subname//': test_write_vars')
    call test_write_vars

    call write_to_log(subname//': test_check_var_or_dim')
    call test_check_var_or_dim

    call write_to_log(subname//': test_read_vars')
    call test_read_vars

    call write_to_log(subname//': test_read_vars_change_type')
    call test_read_vars_change_type

    call clean

  end subroutine test_ncdio_pio

  !-----------------------------------------------------------------------
  subroutine create_vars(bounds)
    !
    ! !DESCRIPTION:
    ! Create module data used in tests here
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g
    integer :: lev1, lev2
    integer :: cur_intlogical

    character(len=*), parameter :: subname = 'create_vars'
    !-----------------------------------------------------------------------

    associate( &
         begg => bounds%begg, &
         endg => bounds%endg  &
         )

    allocate(data_double_1d_grc(begg:endg))
    do g = begg, endg
       data_double_1d_grc(g) = grc%latdeg(g)*1.1_r8
    end do

    allocate(data_int_1d_grc(begg:endg))
    do g = begg, endg
       data_int_1d_grc(g) = int(grc%latdeg(g))*3
    end do

    ! This definition of the logical values will lead to differences based on processor
    ! count. That's okay for now, but we may need to revisit this if it causes problems
    ! in the future.
    allocate(data_intlogical_1d_grc(begg:endg))
    cur_intlogical = 1
    do g = begg, endg
       data_intlogical_1d_grc(g) = cur_intlogical
       cur_intlogical = 1 - cur_intlogical
    end do

    allocate(data_double_2d_grc(begg:endg, 1:nlev1))
    do lev1 = 1, nlev1
       do g = begg, endg
          data_double_2d_grc(g, lev1) = grc%latdeg(g)*4.4_r8 + lev1*5.5_r8
       end do
    end do

    allocate(data_double_3d_grc(begg:endg, 1:nlev1, 1:nlev2))
    do lev2 = 1, nlev2
       do lev1 = 1, nlev1
          do g = begg, endg
             data_double_3d_grc(g, lev1, lev2) = grc%latdeg(g)*4.4_r8 + lev1*5.5_r8 + lev2*6.6_r8
          end do
       end do
    end do

    end associate

  end subroutine create_vars

  !-----------------------------------------------------------------------
  subroutine test_write_vars
    !
    ! !DESCRIPTION:
    ! Write all variables to a test file; confirm that the variables are written and that
    ! they are the expected type
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid
    integer :: numg
    integer :: dimid

    character(len=*), parameter :: subname = 'test_write_vars'
    !-----------------------------------------------------------------------

    call ncd_pio_createfile(ncid, testfilename)

    call get_proc_global(ng=numg)
    call ncd_defdim(ncid, nameg, numg, dimid)
    call ncd_defdim(ncid, lev1_name, nlev1, dimid)
    call ncd_defdim(ncid, lev2_name, nlev2, dimid)

    ! Write latitude since many of the variables' values depend on it (in case we want to
    ! double-check those values)
    call ncd_defvar(ncid=ncid, varname='latdeg', xtype=ncd_double, &
         dim1name=nameg, &
         long_name='latitude', units='degrees')

    call ncd_defvar(ncid=ncid, varname='data_double_1d_grc', xtype=ncd_double, &
         dim1name=nameg, &
         long_name=' ', units=' ')
    call ncd_defvar(ncid=ncid, varname='data_double_1d_grc_float', xtype=ncd_float, &
         dim1name=nameg, &
         long_name=' ', units=' ')
    call ncd_defvar(ncid=ncid, varname='data_int_1d_grc', xtype=ncd_int, &
         dim1name=nameg, &
         long_name=' ', units=' ')
    call ncd_defvar(ncid=ncid, varname='data_intlogical_1d_grc', xtype=ncd_int, &
         dim1name=nameg, &
         long_name=' ', units=' ')
    call ncd_defvar(ncid=ncid, varname='data_double_2d_grc', xtype=ncd_double, &
         dim1name=nameg, dim2name=lev1_name, &
         long_name=' ', units=' ')
    call ncd_defvar(ncid=ncid, varname='data_double_2d_grc_switchdim', xtype=ncd_double, &
         dim1name=lev1_name, dim2name=nameg, &
         long_name=' ', units=' ')
    call ncd_defvar(ncid=ncid, varname='data_double_3d_grc', xtype=ncd_double, &
         dim1name=nameg, dim2name=lev1_name, dim3name=lev2_name)

    call ncd_enddef(ncid)

    call write_to_log('Writing latdeg')
    call ncd_io(varname='latdeg', data=grc%latdeg, &
         dim1name=nameg, ncid=ncid, flag='write')

    call write_to_log('Writing data_double_1d_grc')
    call ncd_io(varname='data_double_1d_grc', data=data_double_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='write')
    call write_to_log('Writing data_double_1d_grc_float')
    call ncd_io(varname='data_double_1d_grc_float', data=data_double_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='write')
    call write_to_log('Writing data_int_1d_grc')
    call ncd_io(varname='data_int_1d_grc', data=data_int_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='write')
    call write_to_log('Writing data_intlogical_1d_grc')
    call ncd_io(varname='data_intlogical_1d_grc', data=data_intlogical_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='write')
    call write_to_log('Writing data_double_2d_grc')
    call ncd_io(varname='data_double_2d_grc', data=data_double_2d_grc, &
         dim1name=nameg, ncid=ncid, flag='write')
    call write_to_log('Writing data_double_2d_grc_switchdim')
    call ncd_io(varname='data_double_2d_grc_switchdim', data=data_double_2d_grc, &
         dim1name=nameg, switchdim=.true., ncid=ncid, flag='write')
    call write_to_log('Writing data_double_3d_grc')
    call ncd_io(varname='data_double_3d_grc', data=data_double_3d_grc, &
         dim1name=nameg, ncid=ncid, flag='write')

    call ncd_pio_closefile(ncid)

    call write_to_log('Confirming that all variables have been written with correct type')
    call ncd_pio_openfile(ncid, testfilename, 0)
    call confirm_var_on_file(ncid, 'data_double_1d_grc', ncd_double)
    call confirm_var_on_file(ncid, 'data_double_1d_grc_float', ncd_float)
    call confirm_var_on_file(ncid, 'data_int_1d_grc', ncd_int)
    call confirm_var_on_file(ncid, 'data_intlogical_1d_grc', ncd_int)
    call confirm_var_on_file(ncid, 'data_double_2d_grc', ncd_double)
    call confirm_var_on_file(ncid, 'data_double_2d_grc_switchdim', ncd_double)
    call confirm_var_on_file(ncid, 'data_double_3d_grc', ncd_double)
    call ncd_pio_closefile(ncid)

  end subroutine test_write_vars

  !-----------------------------------------------------------------------
  subroutine confirm_var_on_file(ncid, varname, expected_type)
    !
    ! !DESCRIPTION:
    ! Confirm that the given variable exists on the file, with the expected type
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    character(len=*), intent(in) :: varname
    integer, intent(in) :: expected_type
    !
    ! !LOCAL VARIABLES:
    logical :: readvar
    integer :: vartype
    type(var_desc_t) :: vardesc

    character(len=*), parameter :: subname = 'check_var_written'
    !-----------------------------------------------------------------------

    call check_var(ncid, varname, readvar, vardesc=vardesc)
    if (.not. readvar) then
       call endrun(trim(varname)//' not found on file')
    end if

    call ncd_inqvtype(ncid, vardesc, vartype)
    if (vartype /= expected_type) then
       call endrun(trim(varname)//' not expected type')
    end if

  end subroutine confirm_var_on_file

  !-----------------------------------------------------------------------
  subroutine test_check_var_or_dim()
    !
    ! !DESCRIPTION:
    ! Test the check_var_or_dim subroutine with variables and dimensions, returning true
    ! and false
    !
    ! This also covers check_var and check_dim
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid
    logical :: exists

    character(len=*), parameter :: var_to_check = 'data_double_1d_grc'
    character(len=*), parameter :: dim_to_check = lev1_name

    character(len=*), parameter :: subname = 'test_check_var_or_dim'
    !-----------------------------------------------------------------------

    call ncd_pio_openfile(ncid, testfilename, 0)

    call check_var_or_dim(ncid, var_to_check, is_dim=.false., exists=exists)
    call shr_assert(exists, 'check_var_or_dim: var exists')

    ! Make sure that check_var_or_dim returns false when the given variable doesn't exist
    ! - even if it is an existing dimension
    call check_var_or_dim(ncid, dim_to_check, is_dim=.false., exists=exists)
    call shr_assert(.not. exists, 'check_var_or_dim: var does not exist')

    call check_var_or_dim(ncid, dim_to_check, is_dim=.true., exists=exists)
    call shr_assert(exists, 'check_var_or_dim: dim exists')

    ! Make sure that check_var_or_dim returns false when the given dimension doesn't
    ! exist - even if it is an existing variable
    call check_var_or_dim(ncid, var_to_check, is_dim=.true., exists=exists)
    call shr_assert(.not. exists, 'check_var_or_dim: dim does not exist')

    call ncd_pio_closefile(ncid)
  end subroutine test_check_var_or_dim

  !-----------------------------------------------------------------------
  subroutine test_read_vars()
    !
    ! !DESCRIPTION:
    ! Test reading the variables from file into variables in memory; ensure these match
    ! the originals.
    !
    ! This just tests reading a variable of a given type on file into a variable of the
    ! same type in memory. Tests involving type conversion happen elsewhere.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid
    real(r8), pointer :: local_double_1d_grc(:)
    integer, pointer :: local_int_1d_grc(:)
    real(r8), pointer :: local_double_2d_grc(:,:)
    real(r8), pointer :: local_double_2d_grc_switchdim(:,:)
    real(r8), pointer :: local_double_3d_grc(:,:,:)

    character(len=*), parameter :: subname = 'test_read_vars'
    !-----------------------------------------------------------------------

    allocate(local_double_1d_grc(size(data_double_1d_grc, 1)))
    allocate(local_int_1d_grc(size(data_int_1d_grc, 1)))
    allocate(local_double_2d_grc(size(data_double_2d_grc, 1), size(data_double_2d_grc, 2)))
    allocate(local_double_2d_grc_switchdim(size(data_double_2d_grc, 1), size(data_double_2d_grc, 2)))
    allocate(local_double_3d_grc(size(data_double_3d_grc, 1), size(data_double_3d_grc, 2), size(data_double_3d_grc, 3)))

    call ncd_pio_openfile(ncid, testfilename, 0)

    call write_to_log(subname//': Reading and comparing data_double_1d_grc')
    call ncd_io(varname='data_double_1d_grc', data=local_double_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    call assert_equal(expected=data_double_1d_grc, actual=local_double_1d_grc, &
         msg='data_double_1d_grc')

    call write_to_log(subname//': Reading and comparing data_int_1d_grc')
    call ncd_io(varname='data_int_1d_grc', data=local_int_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    call assert_equal(expected=data_int_1d_grc, actual=local_int_1d_grc, &
         msg='data_int_1d_grc')

    call write_to_log(subname//': Reading and comparing data_double_2d_grc')
    call ncd_io(varname='data_double_2d_grc', data=local_double_2d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    call assert_equal(expected=data_double_2d_grc, actual=local_double_2d_grc, &
         msg='data_double_2d_grc')

    call write_to_log(subname//': Reading and comparing data_double_2d_grc_switchdim')
    call ncd_io(varname='data_double_2d_grc_switchdim', data=local_double_2d_grc_switchdim, &
         dim1name=nameg, switchdim=.true., ncid=ncid, flag='read')
    call assert_equal(expected=data_double_2d_grc, actual=local_double_2d_grc_switchdim, &
         msg='data_double_2d_grc_switchdim')

    call write_to_log(subname//': Reading and comparing data_double_3d_grc')
    call ncd_io(varname='data_double_3d_grc', data=local_double_3d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    call assert_equal(expected=data_double_3d_grc, actual=local_double_3d_grc, &
         msg='data_double_3d_grc')

    call ncd_pio_closefile(ncid)

    deallocate(local_double_1d_grc)
    deallocate(local_int_1d_grc)
    deallocate(local_double_2d_grc)
    deallocate(local_double_2d_grc_switchdim)
    deallocate(local_double_3d_grc)

  end subroutine test_read_vars

  !-----------------------------------------------------------------------
  subroutine test_read_vars_change_type()
    !
    ! !DESCRIPTION:
    ! Test reading some variables from file into variables in memory with different
    ! types, to ensure this type conversion is done correctly.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid
    real(r8), pointer :: local_double_1d_grc(:)
    integer, pointer :: local_int_1d_grc(:)
    logical, pointer :: local_log_1d_grc(:)
    logical, allocatable :: expected_log_1d(:)

    character(len=*), parameter :: subname = 'test_read_vars_change_type'
    !-----------------------------------------------------------------------

    allocate(local_double_1d_grc(size(data_double_1d_grc, 1)))
    allocate(local_int_1d_grc(size(data_double_1d_grc, 1)))
    allocate(local_log_1d_grc(size(data_double_1d_grc, 1)))

    call ncd_pio_openfile(ncid, testfilename, 0)

    call write_to_log(subname//': Reading float into double')
    local_double_1d_grc(:) = 0._r8
    call ncd_io(varname='data_double_1d_grc_float', data=local_double_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    call assert_equal(expected=data_double_1d_grc, actual=local_double_1d_grc, &
         msg='data_double_1d_grc_float to double', abs_tol=1.e-4_r8)

    call write_to_log(subname//': Reading int into double')
    local_double_1d_grc(:) = 0._r8
    call ncd_io(varname='data_int_1d_grc', data=local_double_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    call assert_equal(expected=real(data_int_1d_grc, r8), actual=local_double_1d_grc, &
         msg='data_int_1d_grc to double')

    call write_to_log(subname//': Reading double into int')
    local_int_1d_grc(:) = 0._r8
    call ncd_io(varname='data_double_1d_grc', data=local_int_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    call assert_equal(expected=int(data_double_1d_grc), actual=local_int_1d_grc, &
         msg='data_double_1d_grc to int')

    call write_to_log(subname//': Reading float into int')
    local_int_1d_grc(:) = 0._r8
    call ncd_io(varname='data_double_1d_grc_float', data=local_int_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    call assert_equal(expected=int(data_double_1d_grc), actual=local_int_1d_grc, &
         msg='data_double_1d_grc_float to int')

    call write_to_log(subname//': Reading int into logical')
    local_log_1d_grc(:) = .false.
    call ncd_io(varname='data_intlogical_1d_grc', data=local_log_1d_grc, &
         dim1name=nameg, ncid=ncid, flag='read')
    allocate(expected_log_1d(size(data_intlogical_1d_grc)))
    call convert_to_logical(data_intlogical_1d_grc, expected_log_1d)
    call assert_equal(expected=expected_log_1d, actual=local_log_1d_grc, &
         msg='data_intlogical_1d_grc to logical')

    call ncd_pio_closefile(ncid)

    deallocate(local_double_1d_grc)
    deallocate(local_int_1d_grc)
    deallocate(local_log_1d_grc)

  end subroutine test_read_vars_change_type


  !-----------------------------------------------------------------------
  subroutine write_to_log(msg)
    !
    ! !DESCRIPTION:
    ! Write a message to the log file, just from the masterproc
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: msg
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'write_to_log'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(*,'(a)') msg
    end if

  end subroutine write_to_log

  !-----------------------------------------------------------------------
  subroutine clean
    !
    ! !DESCRIPTION:
    ! Do end-of-testing cleanup
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'clean'
    !-----------------------------------------------------------------------

    deallocate(data_double_1d_grc)
    deallocate(data_int_1d_grc)
    deallocate(data_intlogical_1d_grc)
    deallocate(data_double_2d_grc)
    deallocate(data_double_3d_grc)

  end subroutine clean


end module TestNcdioPio
