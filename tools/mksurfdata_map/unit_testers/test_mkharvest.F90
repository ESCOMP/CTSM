module test_mkharvest
! Module for testing harvest

  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkharvestMod
  use test_mod

  implicit none
  private

  public :: test_harvest_init
  public :: test_harvest_init_old
  public :: test_harvest_data
  public :: test_harvest_data_all1D

  character(len=*), parameter :: modname = 'test_harvest'

  character(len=128) :: testname
  character(len=128) :: test_prefix
  integer, parameter :: ns_o = 4
  
contains
  
!------------------------------------------------------------------------------
  subroutine test_harvest_init

    use mkncdio
    implicit none
    
    integer            :: ncid
    type(harvestDataType) :: harvdata
    character(len=128) :: varname
    integer            :: varid
    logical            :: varexists
    integer            :: ifld
    character(len=*), parameter :: constfieldname(9) = (/           &
                                       'CONST_HARVEST_VH1      ',   &
                                       'CONST_HARVEST_VH2      ',   &
                                       'CONST_HARVEST_SH1      ',   &
                                       'CONST_HARVEST_SH2      ',   &
                                       'CONST_HARVEST_SH3      ',   &
                                       'CONST_GRAZING          ',   &
                                       'CONST_FERTNITRO_CFT    ',   &
                                       'UNREPRESENTED_PFT_LULCC',   &
                                       'UNREPRESENTED_CFT_LULCC'    &
                                                        /)
    character(len=*), parameter :: units(9) = (/     &
                                       'gC/m2/yr',   &
                                       'gC/m2/yr',   &
                                       'gC/m2/yr',   &
                                       'gC/m2/yr',   &
                                       'gC/m2/yr',   &
                                       'gC/m2/yr',   &
                                       'gN/m2/yr',   &
                                       'unitless',   &
                                       'unitless'    &
                                                            /)
    character(len=*), parameter :: fieldname(9) = (/                                            &
                                 'HARVEST_VH1  ', &
                                 'HARVEST_VH2  ', &
                                 'HARVEST_SH1  ', &
                                 'HARVEST_SH2  ', &
                                 'HARVEST_SH3  ', &
                                 'GRAZING      ', &
                                 'FERTNITRO_CFT', &
                                 'PFT_LULCC    ', &
                                 'CFT_LULCC    '  &
                                                           /)
    character(len=*), parameter :: longname(9) = (/                                            &
                                 'harvest from primary forest                             ', &
                                 'harvest from primary non-forest                         ', &
                                 'harvest from secondary mature-forest                    ', &
                                 'harvest from secondary young-forest                     ', &
                                 'harvest from secondary non-forest                       ', &
                                 'grazing of herbacous pfts                               ', &
                                 'constant background nitrogen fertilizer for each crop   ', &
                                 'constant background unrepresented PFT LULCC transitions ', &
                                 'constant background unrepresented crop LULCC transitions'  &
                                                           /)
    character(len=256) :: string
    character(len=*), parameter :: filename = 'unit_testers/inputs/harvestfields.nc'

    character(len=*), parameter :: subname = 'test_harvest_init'
    integer :: nfields

    testname = 'check harvest_init'
    test_prefix = modname//' -- '//subname//' -- '//trim(testname)//' -- '
    ! Open netcdf file that will be used for most tests
    call check_ret(nf_open(filename, 0, ncid), subname)
    varname = 'GRAZING'
    call check_ret(nf_inq_varid(ncid, varname, varid), subname, varexists=varexists)
    call test_is(varexists, trim(test_prefix)//'existing var')
    call test_is( .not.mkharvest_fieldInBounds( 3 ), trim(test_prefix)//'allfieldsoutofboundsbeforeinit')

    call mkharvest_init( ns_o, 0.0_r8, harvdata, filename )
    call test_is( .not.mkharvest_fieldInBounds( 0 ), trim(test_prefix)//'0 out of bounds')
    nfields = mkharvest_numtypes()
    call test_is( .not.mkharvest_fieldInBounds( nfields+1), trim(test_prefix)//'10 out of bounds')

    ! make sure can now do getter functions

    do ifld = 1, mkharvest_numtypes()
       call test_is(mkharvest_fieldname(ifld,constant=.true.), constfieldname(ifld), trim(test_prefix)//'bad const fieldname')
       call test_is(mkharvest_fieldname(ifld), fieldname(ifld), trim(test_prefix)//trim(testname)//'bad fieldname')
       call test_is(mkharvest_units(ifld), units(ifld), trim(test_prefix)//'bad units')
       call test_is(mkharvest_longname(ifld), longname(ifld), trim(test_prefix)//'bad longname')
    end do
    call harvdata%clean()

  end subroutine test_harvest_init

  subroutine test_harvest_data_all1D()
    implicit none
    type(harvestDataType) :: harvdata
    integer :: dim2nd(9)
    integer :: dsizes(2), nfields, ifld, n, doutsizes(2)
    integer :: dims1D(1), dims2D(2)
    character(len=*), parameter :: subname = 'test_harvest_data'
    character(len=*), parameter :: filename = 'unit_testers/inputs/harvestfields.nc'
    integer, parameter :: indices1D(9) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /)
    integer, parameter :: indices2D(1) = (/ -1 /)
    real(r8), pointer :: data1D(:)
    integer, allocatable :: ind1D(:), ind2D(:)
    integer, parameter :: ns_i = 15, ns_o = 10

    testname = 'check harvest_data_all1D'
    test_prefix = modname//' -- '//subname//' -- '//trim(testname)//' -- '
    dim2nd = 0
    call mkharvest_init( ns_o, 0.0_r8, harvdata, filename )
    call harvdata%clean()
    call harvdata%init( dim2nd, ns_i, ns_o, 0.0_r8 )
    do ifld = 1, mkharvest_numtypes()
         call test_is(harvdata%isField1D(ifld), trim(test_prefix)//'field is 1D' )
         call test_is(.not.harvdata%isField2D(ifld), trim(test_prefix)//'field not 2D' )
    end do
    nfields = mkharvest_numtypes()
    call test_is(harvdata%num1DFields(),nfields,trim(test_prefix)//'num 1D fields')
    call test_is(harvdata%num2DFields(),0,trim(test_prefix)//'num 2D fields')
    call harvdata%getFieldsIdx( ind1D, ind2D )
    call test_is(ind1D,indices1D,trim(test_prefix)//'1D fields indices')
    call test_is(ind2D,indices2D,trim(test_prefix)//'2D fields indices')
    dsizes(1) = ns_i
    doutsizes(1) = ns_o
    do n = 1, harvdata%num1DFields()
       call test_is(harvdata%isField1D(indices1D(n)), trim(test_prefix)//'verify field is 1D' )
       data1D => harvdata%get1DFieldPtr( indices1D(n) )
       dims1D = shape(data1D)
       call test_is(dims1D,dsizes(:),trim(test_prefix)//'1D field dims')
       ! Set data
       data1D(:) = real( n, r8 )
       data1D => null()
       ! Output data
       data1D => harvdata%get1DFieldPtr( indices1D(n), output=.true. )
       dims1D = shape(data1D)
       call test_is(dims1D,doutsizes(:),trim(test_prefix)//'1D Output field dims')
       ! Set data
       data1D(:) = real( n*100, r8 )
       data1D => null()
    end do
    ! Check that data is set from setting above
    do n = 1, harvdata%num1DFields()
       data1D => harvdata%get1DFieldPtr( indices1D(n) )
       call test_is(data1D(1),real( n, r8 ), trim(test_prefix)//'field ')
       data1D => null()
       ! output data
       data1D => harvdata%get1DFieldPtr( indices1D(n), output=.true. )
       call test_is(data1D(1),real( n*100, r8 ), trim(test_prefix)//'field ')
       data1D => null()
    end do
    call harvdata%clean()
  end subroutine test_harvest_data_all1D

!------------------------------------------------------------------------------

  subroutine test_harvest_data()
    implicit none
    type(harvestDataType) :: harvdata
    integer :: dsizes(2), nfields, ifld, n, doutsizes(2)
    integer :: dims1D(1), dims2D(2)
    character(len=*), parameter :: subname = 'test_harvest_data'
    character(len=*), parameter :: filename = 'unit_testers/inputs/harvestfields.nc'
    integer, parameter :: indices1D(6) = (/ 1, 2, 3, 4, 5, 6 /)
    integer, parameter :: indices2D(3) = (/  7,  8,  9 /)
    integer, parameter :: dim2nd(3)    = (/ 64, 15, 64 /)
    character(len=10)  :: dimnames(3)  = (/ "cft", "natpft", "cft" /)
    real(r8), pointer :: data1D(:)
    real(r8), pointer :: data2D(:,:)
    integer, allocatable :: ind1D(:), ind2D(:)
    integer, parameter :: ns_i = 4, ns_o = 20

    testname = 'check harvest_data'
    test_prefix = modname//' -- '//subname//' -- '//trim(testname)//' -- '
    call mkharvest_init( ns_o, 0.0_r8, harvdata, filename )
    call harvdata%getFieldsIdx( ind1D, ind2D )
    call test_is(ind1D,indices1D,trim(test_prefix)//'1D fields indices')
    call test_is(ind2D,indices2D,trim(test_prefix)//'2D fields indices')
    call test_is(harvdata%num1DFields(),size(indices1D),trim(test_prefix)//'num 1D fields')
    call test_is(harvdata%num2DFields(),size(indices2D),trim(test_prefix)//'num 2D fields')
    do n = 1, harvdata%num1DFields()
         ifld = ind1D(n)
         call test_is(harvdata%isField1D(ifld), trim(test_prefix)//'field is 1D' )
         call test_is(.not.harvdata%isField2D(ifld), trim(test_prefix)//'field not 2D' )
    end do
    do n = 1, harvdata%num2DFields()
         ifld = ind2D(n)
         call test_is(.not.harvdata%isField1D(ifld), trim(test_prefix)//'field is not 1D' )
         call test_is(harvdata%isField2D(ifld), trim(test_prefix)//'field is 2D' )
    end do
    dsizes(1) = ns_i
    doutsizes(1) = ns_o
    do n = 1, harvdata%num1DFields()
       call test_is(harvdata%isField1D(indices1D(n)), trim(test_prefix)//'verify field is 1D' )
       data1D => harvdata%get1DFieldPtr( indices1D(n) )
       dims1D = shape(data1D)
       call test_is(dims1D,dsizes(:),trim(test_prefix)//'1D field dims')
       call test_is(harvdata%getFieldsDim(indices1D(n)),"none",trim(test_prefix)//'1D field dimname')
       data1D => null()
    end do
    do n = 1, harvdata%num2DFields()
       dsizes(2) = dim2nd(n)
       call test_is(harvdata%isField2D(indices2D(n)), trim(test_prefix)//'verify field is 2D' )
       data2D => harvdata%get2DFieldPtr( indices2D(n) )
       dims2D = shape(data2D)
       call test_is(dims2D,dsizes(:),trim(test_prefix)//'2D field dims')
       call test_is(harvdata%getFieldsDim(indices2D(n)),dimnames(n),trim(test_prefix)//'1D field dimname')
       data2D => null()
    end do
    call harvdata%clean()
  end subroutine test_harvest_data


!------------------------------------------------------------------------------
  subroutine test_harvest_init_old

    use mkncdio
    implicit none
    
    type(harvestDataType) :: harvdata
    character(len=128) :: testname
    integer            :: ncid
    character(len=128) :: varname
    integer            :: varid
    logical            :: varexists
    integer, parameter :: ns_o = 4
    integer            :: ifld

    character(len=*), parameter :: filename = 'unit_testers/inputs/harvestfieldsold.nc'

    character(len=*), parameter :: subname = 'test_harvest_init'
    character(len=*), parameter :: constfieldname(9) = (/                &
                                       'CONST_HARVEST_VH1      ',   &
                                       'CONST_HARVEST_VH2      ',   &
                                       'CONST_HARVEST_SH1      ',   &
                                       'CONST_HARVEST_SH2      ',   &
                                       'CONST_HARVEST_SH3      ',   &
                                       'CONST_GRAZING          ',   &
                                       'CONST_FERTNITRO_CFT    ',   &
                                       'UNREPRESENTED_PFT_LULCC',   &
                                       'UNREPRESENTED_CFT_LULCC'    &
                                                            /)
    character(len=*), parameter :: units(9) = (/     &
                                       'unitless   ',   &
                                       'unitless   ',   &
                                       'unitless   ',   &
                                       'unitless   ',   &
                                       'unitless   ',   &
                                       'unitless   ',   &
                                       'not_read_in',   &
                                       'not_read_in',   &
                                       'not_read_in'    &
                                                            /)
    character(len=*), parameter :: fieldname(9) = (/                                            &
                                 'HARVEST_VH1  ', &
                                 'HARVEST_VH2  ', &
                                 'HARVEST_SH1  ', &
                                 'HARVEST_SH2  ', &
                                 'HARVEST_SH3  ', &
                                 'GRAZING      ', &
                                 'FERTNITRO_CFT', &
                                 'PFT_LULCC    ', &
                                 'CFT_LULCC    '  &
                                                           /)
    character(len=*), parameter :: longname(9) = (/                             &
                                 'harvest from primary forest         ', &
                                 'harvest from primary non-forest     ', &
                                 'harvest from secondary mature-forest', &
                                 'harvest from secondary young-forest ', &
                                 'harvest from secondary non-forest   ', &
                                 'grazing of herbacous pfts           ', &
                                 'FERTNITRO_CFT (zeroed out)          ', &
                                 'PFT_LULCC (zeroed out)              ', &
                                 'CFT_LULCC (zeroed out)              '  &
                                                           /)
    character(len=256) :: string
    testname = 'check harvest_init_old'
    ! Open netcdf file that will be used for most tests
    call check_ret(nf_open(filename, 0, ncid), subname)
    varname = 'GRAZING'
    call check_ret(nf_inq_varid(ncid, varname, varid), subname, varexists=varexists)
    call test_is(varexists, modname//' -- '//subname//' -- '//trim(testname)//' -- existing var')

    call mkharvest_init( ns_o, 0.0_r8, harvdata, filename )

    ! make sure can now do getter functions

    do ifld = 1, mkharvest_numtypes()
       call test_is(mkharvest_fieldname(ifld,constant=.true.), constfieldname(ifld), modname//' -- '//subname//' -- '//trim(testname)//' -- bad const fieldname')
       call test_is(mkharvest_fieldname(ifld), fieldname(ifld), modname//' -- '//subname//' -- '//trim(testname)//' -- bad fieldname')
       call test_is(mkharvest_units(ifld), units(ifld), modname//' -- '//subname//' -- '//trim(testname)//' -- bad units')
       call test_is(mkharvest_longname(ifld), longname(ifld), modname//' -- '//subname//' -- '//trim(testname)//' -- bad longname')
    end do
    call harvdata%clean()

  end subroutine test_harvest_init_old

end module test_mkharvest
