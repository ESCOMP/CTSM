module DiurnalOzoneType

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Sets up type for converting multi-day ozone input to sub-daily values using an input
    ! ozone anomaly file
    !
    ! !USES:
  #include "shr_assert.h"
    use shr_kind_mod           , only : r8 => shr_kind_r8
    use decompMod              , only : bounds_type
    use clm_varcon             , only : spval
    use clm_varctl             , only : iulog
    use abortutils             , only : endrun

    implicit none
    save
    private

    ! !PUBLIC TYPES:
    type, public :: diurnal_ozone_anom_type
       private
       ! Private data members
       integer                       :: ntimes               ! size of time dimension
       real(r8), public, allocatable :: o3_anomaly_grc(:,:)  ! o3 anomaly data [grc, ntimes]
       real(r8), public, allocatable :: time_arr(:)          ! time dimension (units = seconds of day)

     contains
       ! Public routines
       procedure, public  :: Init
       procedure, private :: InitAllocate
       procedure, public  :: Interp

    end type diurnal_ozone_anom_type

    character(len=*), parameter, private :: sourcefile = &
    __FILE__


  contains

    ! ========================================================================
    ! Infrastructure routines (initialization, etc.)
    ! ========================================================================

    !-----------------------------------------------------------------------
    subroutine Init(this, bounds, n)
      !
      ! DESCRIPTION:
      ! Initialize ozone anomaly data structures
      !
      !
      ! ARGUMENTS:
      class(diurnal_ozone_anom_type), intent(inout) :: this
      type(bounds_type),              intent(in)    :: bounds
      integer,                        intent(in)    :: n
      !-----------------------------------------------------------------------

      ! set ntimes from input
      this%ntimes = n

      ! allocate arrays
      call this%InitAllocate(bounds)

    end subroutine Init


    !-----------------------------------------------------------------------
    subroutine InitAllocate(this, bounds)
      !
      ! !DESCRIPTION:
      ! Allocate module variables and data structures
      !
      ! !USES:
      use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
      !
      ! !ARGUMENTS:
      class(diurnal_ozone_anom_type), intent(inout) :: this
      type(bounds_type),              intent(in)    :: bounds
      !
      ! LOCAL VARIABLES:
      integer  :: begg, endg
      !---------------------------------------------------------------------

      begg = bounds%begg; endg = bounds%endg

      allocate(this%o3_anomaly_grc(begg:endg,1:this%ntimes)); this%o3_anomaly_grc(:,:) = nan
      allocate(this%time_arr(1:this%ntimes)); this%time_arr(:) = nan

      end subroutine InitAllocate

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine Interp(this, forc_o3, forc_o3_down)
    !
    ! !DESCRIPTION:
    ! Downscale/interpolate multi-day ozone data to subdaily
    !
    ! !USES:
    use clm_time_manager , only : get_curr_time
    use clm_varcon       , only : secspday
    !
    ! !ARGUMENTS:
    class(diurnal_ozone_anom_type), intent(in)  :: this
    real(r8),                       intent(in)  :: forc_o3(bounds%begg:bounds%endg)      ! ozone partial pressure (mol/mol)
    real(r8),                       intent(out) :: forc_o3_down(bounds%begg:bounds%endg) ! ozone partial pressure, downscaled (mol/mol)
    !
    ! LOCAL VARIABLES:
    integer :: t               ! looping index
    integer :: yr              ! year
    integer :: mon             ! month
    integer  :: day            ! day of month
    integer  :: tod            ! time of day (seconds past 0Z)
    integer  :: begg, endg     ! bounds
    integer  :: t_prev         ! previous time index
    real(r8) :: anomaly_val_start
    real(r8) :: anomaly_val_end
    real(r8) :: anomaly_scalar
    real(r8) :: tdiff_end, tdiff_start
    real(r8) :: tdiff
    !-----------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    ! Get current date/time - we really only need seconds
    call get_curr_date(yr, mon, day, tod)

    ! find the time interval we are in
    do t = 1, this%ntimes
      if (real(tod) <= this%time_arr(t)) then
        exit
      end if
    end do

    ! interpolate, checking for edge cases
    if (t == 1) then
      t_prev = this%ntimes
      tdiff_end = secspday - this%time_arr(t_prev) + real(tod)
      tdiff = this%time_arr(t) + secspday - this%time_arr(t_prev)
    else
      t_prev = t - 1
      tdiff_end = real(tod) - this%time_arr(t_prev)
      tdiff = this%time_arr(t) - this%time_arr(t_prev)
    end if

    anomaly_val_start = this%o3_anomaly_grc(begg:endg, t_prev)
    anomaly_val_end = this%o3_anomaly_grc(begg:endg, t)
    tdiff_start = this%time_arr(t) - real(tod)

    ! interpolate
    anomaly_scalar = (anomaly_val_start*tdiff_start + anomaly_val_end*tdiff_end)/tdiff

    forc_o3_down(begg:endg) = forc_o3(begg:endg)*anomaly_scalar

    end subroutine Interp

  !-----------------------------------------------------------------------

  end module DiurnalOzoneType
