module CIsoAtmTimeseriesMod

  !-----------------------------------------------------------------------
  ! Module for transient atmospheric boundary to the c13 and c14 codes
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use clm_time_manager , only : get_curr_date,get_days_per_year, get_curr_yearfrac
  use clm_varcon       , only : c14ratio, secspday
  use shr_const_mod    , only : SHR_CONST_PDB                    ! Ratio of C13/C12
  use clm_varctl       , only : iulog
  use abortutils       , only : endrun
  use spmdMod          , only : masterproc
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: C14BombSpike           ! Time series for C14 data
  public:: C14_init_BombSpike     ! Initialize C14 data series and read data in
  public:: C13Timeseries          ! Time series for C13 data
  public:: C13_init_Timeseries    ! Initialize C13 data series and read data in

  ! !PUBLIC TYPES:
  logical            , public :: use_c14_bombspike = .false.  ! do we use time-varying atmospheric C14?
  character(len=256) , public :: atm_c14_filename = ' '       ! file name of C14 input data
  logical            , public :: use_c13_timeseries = .false. ! do we use time-varying atmospheric C13?
  character(len=256) , public :: atm_c13_filename = ' '       ! file name of C13 input data
  integer, parameter , public :: nsectors_c14 = 3             ! Number of latitude sectors the C14 data has

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private:: check_units   ! Check the units of the data on the input file

  ! !PRIVATE TYPES:
  real(r8), allocatable, private :: atm_c14file_time(:)       ! time for C14 data
  real(r8), allocatable, private :: atm_delta_c14(:,:)        ! Delta C14 data
  real(r8), allocatable, private :: atm_c13file_time(:)       ! time for C13 data
  real(r8), allocatable, private :: atm_delta_c13(:)          ! Delta C13 data
  real(r8), parameter :: time_axis_offset = 1850.0_r8         ! Offset in years of time on file

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine C14BombSpike( rc14_atm )
    !
    ! !DESCRIPTION:
    ! for transient simulation, read in an atmospheric timeseries file to impose bomb spike
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(out) :: rc14_atm(nsectors_c14)  ! Ratio of C14 to C12
    !
    ! !LOCAL VARIABLES:
    integer  :: yr, mon, day, tod            ! year, month, day, time-of-day
    real(r8) :: dateyear                     ! Date converted to year
    real(r8) :: delc14o2_atm(nsectors_c14)   ! C14 delta units
    integer  :: fp, p, nt                    ! Indices
    integer  :: ind_below                    ! Time index below current time
    integer  :: ntim_atm_ts                  ! Number of times on file
    real(r8) :: twt_1, twt_2                 ! weighting fractions for interpolating
    integer  :: l                            ! Loop index of sectors
    !-----------------------------------------------------------------------

    ! get current date
    call get_curr_date(yr, mon, day, tod)
    dateyear = real(yr) + get_curr_yearfrac()

    ! find points in atm timeseries to interpolate between
    ntim_atm_ts = size(atm_c14file_time)
    ind_below = 0
    do nt = 1, ntim_atm_ts
       if ((dateyear - time_axis_offset) >= atm_c14file_time(nt) ) then
          ind_below = ind_below+1
       endif
    end do

    ! loop over lat bands to pass all three to photosynthesis
    do l = 1,nsectors_c14
       ! interpolate between nearest two points in atm c14 timeseries
       if (ind_below .eq. 0 ) then 
          delc14o2_atm(l) = atm_delta_c14(l,1)
       elseif (ind_below .eq. ntim_atm_ts ) then
          delc14o2_atm(l) = atm_delta_c14(l,ntim_atm_ts)
       else
          twt_2 = min(1._r8, max(0._r8,((dateyear - time_axis_offset)-atm_c14file_time(ind_below)) &
               / (atm_c14file_time(ind_below+1)-atm_c14file_time(ind_below))))
          twt_1 = 1._r8 - twt_2
          delc14o2_atm(l) = atm_delta_c14(l,ind_below) * twt_1 +  atm_delta_c14(l,ind_below+1) * twt_2
       endif
    
       ! change delta units to ratio
       rc14_atm(l) = (delc14o2_atm(l) * 1.e-3_r8 + 1._r8) * c14ratio
    end do
    
  end subroutine C14BombSpike

  !-----------------------------------------------------------------------
  subroutine C14_init_BombSpike()
    !
    ! !DESCRIPTION:
    ! read netcdf file containing a timeseries of atmospheric delta C14 values; save in module-level array 
    !
    ! !USES:
    use ncdio_pio   , only : ncd_pio_openfile, ncd_pio_closefile, file_desc_t, ncd_inqdlen, ncd_io
    use fileutils   , only : getfil
    !
    ! !LOCAL VARIABLES:
    implicit none
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: nsec                       ! number of input data sectors
    integer :: t                          ! time index
    logical :: readvar                    ! if variable read or not
    character(len=*), parameter :: vname = 'Delta14co2_in_air'  ! Variable name on file
    !-----------------------------------------------------------------------

    call getfil(atm_c14_filename, locfn, 0)

    if ( masterproc ) then
       write(iulog, *) 'C14_init_BombSpike: preparing to open file:'
       write(iulog, *) trim(locfn)
    endif

    call ncd_pio_openfile (ncid, trim(locfn), 0)

    call ncd_inqdlen(ncid,dimid,ntim,'time')
    call ncd_inqdlen(ncid,dimid,nsec,'sector')
    if ( nsec /= nsectors_c14 )then
       call endrun(msg="ERROR: number of sectors on file not what's expected"//errMsg(sourcefile, __LINE__))
    end if

    !! allocate arrays based on size of netcdf timeseries
    allocate(atm_c14file_time(ntim))
    allocate(atm_delta_c14(nsectors_c14,ntim))
    atm_delta_c14(:,:) = 0.0_r8

    call ncd_io(ncid=ncid, varname='time', flag='read', data=atm_c14file_time, &
                readvar=readvar)
    if ( .not. readvar ) then
       call endrun(msg="ERROR: time not on file"//errMsg(sourcefile, __LINE__))
    end if

    call ncd_io(ncid=ncid, varname=vname, flag='read', data=atm_delta_c14, &
                readvar=readvar)
    if ( .not. readvar ) then
       call endrun(msg="ERROR: '//vname//' not on file"//errMsg(sourcefile, __LINE__))
    end if
    ! Check units
    call check_units( ncid, vname, "Modern" )
    call ncd_pio_closefile(ncid)

    ! check to make sure that time dimension is well behaved
    do t = 2, ntim
       if ( atm_c14file_time(t) - atm_c14file_time(t-1) <= 0._r8 ) then
          write(iulog, *) 'C14_init_BombSpike: error.  time axis must be monotonically increasing'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    end do

  end subroutine C14_init_BombSpike


  !-----------------------------------------------------------------------
  subroutine C13TimeSeries( rc13_atm )
    !
    ! !DESCRIPTION:
    ! for transient pulse simulation, impose a time-varying atm boundary condition
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(out) :: rc13_atm    ! Ratio of C13 to C12
    !
    ! !LOCAL VARIABLES:
    integer  :: yr, mon, day, tod          ! year, month, day, time-of-day
    real(r8) :: dateyear                   ! date translated to year
    real(r8) :: delc13o2_atm               ! Delta C13
    integer  :: fp, p, nt                  ! Indices
    integer  :: ind_below                  ! Index of time in file before current time
    integer  :: ntim_atm_ts                ! Number of times on file
    real(r8) :: twt_1, twt_2               ! weighting fractions for interpolating
    !-----------------------------------------------------------------------

    ! get current date
    call get_curr_date(yr, mon, day, tod)
    dateyear = real(yr) + get_curr_yearfrac()

    ! find points in atm timeseries to interpolate between
    ntim_atm_ts = size(atm_c13file_time)
    ind_below = 0
    do nt = 1, ntim_atm_ts
       if ((dateyear - time_axis_offset) >= atm_c13file_time(nt) ) then
          ind_below = ind_below+1
       endif
    end do

    ! interpolate between nearest two points in atm c13 timeseries
    ! cdknotes. for now and for simplicity, just use the northern hemisphere values (sector 1)
    if (ind_below .eq. 0 ) then 
       delc13o2_atm = atm_delta_c13(1)
    elseif (ind_below .eq. ntim_atm_ts ) then
       delc13o2_atm = atm_delta_c13(ntim_atm_ts)
    else
       twt_2 = min(1._r8, max(0._r8,((dateyear - time_axis_offset)-atm_c13file_time(ind_below)) &
            / (atm_c13file_time(ind_below+1)-atm_c13file_time(ind_below))))
       twt_1 = 1._r8 - twt_2
       delc13o2_atm = atm_delta_c13(ind_below) * twt_1 +  atm_delta_c13(ind_below+1) * twt_2
    endif

    ! change delta units to ratio, put on patch loop

    rc13_atm = (delc13o2_atm * 1.e-3_r8 + 1._r8) * SHR_CONST_PDB

  end subroutine C13TimeSeries

  !-----------------------------------------------------------------------
  subroutine C13_init_TimeSeries()
    !
    ! !DESCRIPTION:
    ! read netcdf file containing a timeseries of atmospheric delta C13 values; save in module-level array 
    !
    ! !USES:
    use ncdio_pio   , only : ncd_pio_openfile, ncd_pio_closefile, file_desc_t, ncd_inqdlen, ncd_io
    use fileutils   , only : getfil
    !
    ! !LOCAL VARIABLES:
    implicit none
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: t                          ! Time index
    logical :: readvar                    ! if variable read or not
    character(len=*), parameter :: vname = 'delta13co2_in_air'  ! Variable name on file
    !-----------------------------------------------------------------------

    call getfil(atm_c13_filename, locfn, 0)

    if ( masterproc ) then
       write(iulog, *) 'C13_init_TimeSeries: preparing to open file:'
       write(iulog, *) trim(locfn)
    endif

    call ncd_pio_openfile (ncid, trim(locfn), 0)

    call ncd_inqdlen(ncid,dimid,ntim,'time')

    !! allocate arrays based on size of netcdf timeseries
    allocate(atm_c13file_time(ntim))
    allocate(atm_delta_c13(ntim))

    call ncd_io(ncid=ncid, varname='time', flag='read', data=atm_c13file_time, &
                readvar=readvar)
    if ( .not. readvar ) then
       call endrun(msg="ERROR: time not on file"//errMsg(sourcefile, __LINE__))
    end if

    call ncd_io(ncid=ncid, varname=vname, flag='read', data=atm_delta_c13, &
                readvar=readvar)
    if ( .not. readvar ) then
       call endrun(msg="ERROR: '//vname//' not on file"//errMsg(sourcefile, __LINE__))
    end if

    ! Check units
    call check_units( ncid, vname, "VPDB" )
    call ncd_pio_closefile(ncid)

    ! check to make sure that time dimension is well behaved
    do t = 2, ntim
       if ( atm_c13file_time(t) - atm_c13file_time(t-1) <= 0._r8 ) then
          write(iulog, *) 'C13_init_TimeSeries: error.  time axis must be monotonically increasing'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    end do

  end subroutine C13_init_TimeSeries

  !-----------------------------------------------------------------------
  subroutine check_units( ncid, vname, relativeto )
    !
    ! !DESCRIPTION:
    ! check that time and data units are what's expected or else abort
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_inqvid, ncd_getatt, var_desc_t
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid       ! netcdf id
    character(len=*),  intent(in)    :: vname      ! Variable name
    character(len=*),  intent(in)    :: relativeto ! What are data units relative to
    !
    ! !LOCAL VARIABLES:
    type(var_desc_t)  :: vardesc          ! variable descriptor
    integer           :: varid            ! variable index
    character(len=50) :: units            ! Data units
    character(len=50) :: t_units_expected ! Time units to expect

    call ncd_inqvid( ncid, 'time', varid, vardesc )
    call ncd_getatt( ncid, varid, "units", units )
    write(t_units_expected,'("years since", I5, "-01-01 0:0:0.0")' ) nint(time_axis_offset)
    if ( trim(units) /= t_units_expected )then
       call endrun(msg="ERROR: time units on file are NOT what's expected"//errMsg(sourcefile, __LINE__))
    end if
    call ncd_inqvid( ncid, vname, varid, vardesc )
    call ncd_getatt( ncid, varid, "units", units )
    if ( trim(units) /= "per mil, relative to "//relativeto )then
       call endrun(msg="ERROR: units on file for "//vname//" are NOT what's expected"//errMsg(sourcefile, __LINE__))
    end if
  end subroutine check_units


end module CIsoAtmTimeseriesMod
