module Wateratm2lndType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water atm2lnd variables that apply to both bulk water
  ! and water tracers.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use decompMod      , only : BOUNDS_SUBGRID_COLUMN, BOUNDS_SUBGRID_GRIDCELL
  use clm_varcon     , only : spval
  use ColumnType     , only : col
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  use WaterTracerUtils, only : AllocateVar1d, CalcTracerFromBulk, CalcTracerFromBulkFixedRatio
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: wateratm2lnd_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: forc_q_not_downscaled_grc     (:)   ! not downscaled atm specific humidity (kg/kg)
     real(r8), pointer :: forc_rain_not_downscaled_grc  (:)   ! not downscaled atm rain rate [mm/s]
     real(r8), pointer :: forc_snow_not_downscaled_grc  (:)   ! not downscaled atm snow rate [mm/s]
     real(r8), pointer :: forc_q_downscaled_col         (:)   ! downscaled atm specific humidity (kg/kg)
     real(r8), pointer :: forc_flood_grc                (:)   ! rof flood (mm/s)
     real(r8), pointer :: forc_rain_downscaled_col      (:)   ! downscaled atm rain rate [mm/s]
     real(r8), pointer :: forc_snow_downscaled_col      (:)   ! downscaled atm snow rate [mm/s]

     real(r8), pointer :: rain_to_snow_conversion_col   (:)   ! amount of rain converted to snow via precipitation repartitioning (mm/s)
     real(r8), pointer :: snow_to_rain_conversion_col   (:)   ! amount of snow converted to rain via precipitation repartitioning (mm/s)

   contains

     procedure, public  :: Init
     procedure, public  :: IsCommunicatedWithCoupler
     procedure, public  :: SetNondownscaledTracers
     procedure, public  :: SetDownscaledTracers
     procedure, public  :: Clean
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type wateratm2lnd_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info, tracer_vars)

    class(wateratm2lnd_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: tracer_vars

    this%info => info

    call this%InitAllocate(bounds, tracer_vars)

    call this%InitHistory(bounds)

    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, tracer_vars)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(water_tracer_container_type), intent(inout) :: tracer_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    !------------------------------------------------------------------------

    call AllocateVar1d(var = this%forc_q_not_downscaled_grc, name = 'forc_q_not_downscaled_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%forc_rain_not_downscaled_grc, name = 'forc_rain_not_downscaled_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%forc_snow_not_downscaled_grc, name = 'forc_snow_not_downscaled_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%forc_q_downscaled_col, name = 'forc_q_downscaled_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN, &
         ival=ival)
    call AllocateVar1d(var = this%forc_flood_grc, name = 'forc_flood_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL, &
         ival=ival)
    call AllocateVar1d(var = this%forc_rain_downscaled_col, name = 'forc_rain_downscaled_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN, &
         ival=ival)
    call AllocateVar1d(var = this%forc_snow_downscaled_col, name = 'forc_snow_downscaled_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN, &
         ival=ival)
    call AllocateVar1d(var = this%rain_to_snow_conversion_col, name = 'rain_to_snow_conversion_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%snow_to_rain_conversion_col, name = 'snow_to_rain_conversion_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize history vars
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod    , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begc, endc
    integer           :: begg, endg
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg


    this%forc_rain_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname=this%info%fname('RAIN_FROM_ATM'), units='mm/s',  &
         avgflag='A', long_name=this%info%lname('atmospheric rain received from atmosphere (pre-repartitioning)'), &
         ptr_lnd=this%forc_rain_not_downscaled_grc)

    this%forc_snow_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname=this%info%fname('SNOW_FROM_ATM'), units='mm/s',  &
         avgflag='A', long_name=this%info%lname('atmospheric snow received from atmosphere (pre-repartitioning)'), &
         ptr_lnd=this%forc_snow_not_downscaled_grc)

    this%forc_q_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname=this%info%fname('QBOT'), units='kg/kg',  &
         avgflag='A', long_name=this%info%lname('atmospheric specific humidity (downscaled to columns in glacier regions)'), &
         ptr_col=this%forc_q_downscaled_col)
    ! Rename of QBOT for Urban intercomparison project
    call hist_addfld1d (fname=this%info%fname('Qair'), units='kg/kg',  &
         avgflag='A', long_name=this%info%lname('atmospheric specific humidity (downscaled to columns in glacier regions)'), &
         ptr_col=this%forc_q_downscaled_col, default='inactive')

    this%forc_flood_grc(begg:endg) = spval
    call hist_addfld1d (fname=this%info%fname('QFLOOD'),  units='mm/s',  &
         avgflag='A', long_name=this%info%lname('runoff from river flooding'), &
         ptr_lnd=this%forc_flood_grc)

    this%forc_rain_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname=this%info%fname('RAIN'), units='mm/s',  &
         avgflag='A', long_name=this%info%lname('atmospheric rain, after rain/snow repartitioning based on temperature'), &
         ptr_col=this%forc_rain_downscaled_col)
    call hist_addfld1d (fname=this%info%fname('Rainf'), units='mm/s',  &
         avgflag='A', long_name=this%info%lname('atmospheric rain, after rain/snow repartitioning based on temperature'), &
         ptr_col=this%forc_rain_downscaled_col, default='inactive')

    call hist_addfld1d (fname=this%info%fname('RAIN_ICE'), units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('atmospheric rain, after rain/snow repartitioning based on temperature (ice landunits only)'), &
         ptr_col=this%forc_rain_downscaled_col, l2g_scale_type='ice', &
         default='inactive')

    this%forc_snow_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname=this%info%fname('SNOW'), units='mm/s',  &
         avgflag='A', long_name=this%info%lname('atmospheric snow, after rain/snow repartitioning based on temperature'), &
         ptr_col=this%forc_snow_downscaled_col)

    call hist_addfld1d (fname=this%info%fname('SNOW_ICE'), units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('atmospheric snow, after rain/snow repartitioning based on temperature (ice landunits only)'), &
         ptr_col=this%forc_snow_downscaled_col, l2g_scale_type='ice', &
         default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type), intent(inout) :: this
    type(bounds_type)     , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    ! Nothing to do for now

  end subroutine InitCold

  !------------------------------------------------------------------------
  pure function IsCommunicatedWithCoupler(this) result(coupled)
    !
    ! !DESCRIPTION:
    ! Returns true if this tracer is received from the coupler. Returns false if this
    ! tracer is just used internally in CTSM, and should be set to some fixed ratio times
    ! the bulk water.
    !
    ! !ARGUMENTS:
    logical :: coupled  ! function result
    class(wateratm2lnd_type), intent(in) :: this
    !-----------------------------------------------------------------------

    coupled = this%info%is_communicated_with_coupler()

  end function IsCommunicatedWithCoupler


  !-----------------------------------------------------------------------
  subroutine SetNondownscaledTracers(this, bounds, bulk)
    !
    ! !DESCRIPTION:
    ! Set tracer values for the non-downscaled atm2lnd water quantities from the bulk quantities
    !
    ! This should only be called for tracers that are not communicated with the coupler
    ! (i.e., for which this%IsCommunicatedWithCoupler() is false). Note that the tracer
    ! values are set to a fixed ratio times the bulk (because we don't have any other
    ! information to go on for these fields).
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    class(wateratm2lnd_type), intent(in) :: bulk
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ratio

    character(len=*), parameter :: subname = 'SetNondownscaledTracers'
    !-----------------------------------------------------------------------

    associate( &
         begg => bounds%begg, &
         endg => bounds%endg &
         )

    ratio = this%info%get_ratio()

    call CalcTracerFromBulkFixedRatio( &
         bulk = bulk%forc_q_not_downscaled_grc(begg:endg), &
         ratio = ratio, &
         tracer = this%forc_q_not_downscaled_grc(begg:endg))

    call CalcTracerFromBulkFixedRatio( &
         bulk = bulk%forc_rain_not_downscaled_grc(begg:endg), &
         ratio = ratio, &
         tracer = this%forc_rain_not_downscaled_grc(begg:endg))

    call CalcTracerFromBulkFixedRatio( &
         bulk = bulk%forc_snow_not_downscaled_grc(begg:endg), &
         ratio = ratio, &
         tracer = this%forc_snow_not_downscaled_grc(begg:endg))

    call CalcTracerFromBulkFixedRatio( &
         bulk = bulk%forc_flood_grc(begg:endg), &
         ratio = ratio, &
         tracer = this%forc_flood_grc(begg:endg))

    end associate

  end subroutine SetNondownscaledTracers

  !-----------------------------------------------------------------------
  subroutine SetDownscaledTracers(this, bounds, num_allc, filter_allc, &
       bulk)
    !
    ! !DESCRIPTION:
    ! Set tracer values for the downscaled atm2lnd water quantities from the bulk quantities
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type) , intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_allc       ! number of column points in filter_allc
    integer                  , intent(in)    :: filter_allc(:) ! column filter for all points
    class(wateratm2lnd_type) , intent(in)    :: bulk
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, g

    character(len=*), parameter :: subname = 'SetDownscaledTracers'
    !-----------------------------------------------------------------------

    associate( &
         begg => bounds%begg, &
         endg => bounds%endg, &
         begc => bounds%begc, &
         endc => bounds%endc &
         )

    call SetOneDownscaledTracer( &
         bulk_not_downscaled = bulk%forc_q_not_downscaled_grc(begg:endg), &
         tracer_not_downscaled = this%forc_q_not_downscaled_grc(begg:endg), &
         bulk_downscaled = bulk%forc_q_downscaled_col(begc:endc), &
         tracer_downscaled = this%forc_q_downscaled_col(begc:endc))

    call SetOneDownscaledTracer( &
         bulk_not_downscaled = bulk%forc_rain_not_downscaled_grc(begg:endg), &
         tracer_not_downscaled = this%forc_rain_not_downscaled_grc(begg:endg), &
         bulk_downscaled = bulk%rain_to_snow_conversion_col(begc:endc), &
         tracer_downscaled = this%rain_to_snow_conversion_col(begc:endc))

    call SetOneDownscaledTracer( &
         bulk_not_downscaled = bulk%forc_snow_not_downscaled_grc(begg:endg), &
         tracer_not_downscaled = this%forc_snow_not_downscaled_grc(begg:endg), &
         bulk_downscaled = bulk%snow_to_rain_conversion_col(begc:endc), &
         tracer_downscaled = this%snow_to_rain_conversion_col(begc:endc))

    do fc = 1, num_allc
       c = filter_allc(fc)
       g = col%gridcell(c)
       this%forc_rain_downscaled_col(c) = AdjustPrecipTracer( &
            not_downscaled = this%forc_rain_not_downscaled_grc(g), &
            addition = this%snow_to_rain_conversion_col(c), &
            subtraction = this%rain_to_snow_conversion_col(c))
       this%forc_snow_downscaled_col(c) = AdjustPrecipTracer( &
            not_downscaled = this%forc_snow_not_downscaled_grc(g), &
            addition = this%rain_to_snow_conversion_col(c), &
            subtraction = this%snow_to_rain_conversion_col(c))
    end do

    end associate

  contains
    subroutine SetOneDownscaledTracer(bulk_not_downscaled, tracer_not_downscaled, &
         bulk_downscaled, tracer_downscaled)
      real(r8), intent(in) :: bulk_not_downscaled( bounds%begg: )
      real(r8), intent(in) :: tracer_not_downscaled( bounds%begg: )
      real(r8), intent(in) :: bulk_downscaled( bounds%begc: )
      real(r8), intent(inout) :: tracer_downscaled( bounds%begc: )

      integer  :: fc, c, g
      real(r8) :: bulk_not_downscaled_col(bounds%begc:bounds%endc)
      real(r8) :: tracer_not_downscaled_col(bounds%begc:bounds%endc)

      SHR_ASSERT_ALL((ubound(bulk_not_downscaled) == [bounds%endg]), errMsg(sourcefile, __LINE__))
      SHR_ASSERT_ALL((ubound(tracer_not_downscaled) == [bounds%endg]), errMsg(sourcefile, __LINE__))
      SHR_ASSERT_ALL((ubound(bulk_downscaled) == [bounds%endc]), errMsg(sourcefile, __LINE__))
      SHR_ASSERT_ALL((ubound(tracer_downscaled) == [bounds%endc]), errMsg(sourcefile, __LINE__))

      associate( &
           begc => bounds%begc, &
           endc => bounds%endc &
           )

      do fc = 1, num_allc
         c = filter_allc(fc)
         g = col%gridcell(c)
         ! Note that this copying of bulk_not_downscaled to bulk_not_downscaled_col will
         ! be repeated for every tracer. At some point we might want to optimize this so
         ! that it's just done once and shared for all tracers (probably by doing this
         ! copy outside of the loop over tracers that calls SetDownscaledTracers).
         bulk_not_downscaled_col(c) = bulk_not_downscaled(g)
         tracer_not_downscaled_col(c) = tracer_not_downscaled(g)
      end do

      call CalcTracerFromBulk( &
           lb = begc, &
           num_pts = num_allc, &
           filter_pts = filter_allc, &
           bulk_source = bulk_not_downscaled_col(begc:endc), &
           bulk_val = bulk_downscaled(begc:endc), &
           tracer_source = tracer_not_downscaled_col(begc:endc), &
           tracer_val = tracer_downscaled(begc:endc))

      end associate

    end subroutine SetOneDownscaledTracer

    pure function AdjustPrecipTracer(not_downscaled, addition, subtraction) result(downscaled)
      real(r8) :: downscaled
      real(r8), intent(in) :: not_downscaled
      real(r8), intent(in) :: addition
      real(r8), intent(in) :: subtraction

      real(r8), parameter :: tolerance = 1.e-13_r8

      downscaled = not_downscaled + addition - subtraction
      if (not_downscaled /= 0._r8) then
         if (abs(downscaled / not_downscaled) < tolerance) then
            ! Roundoff correction: If it seems that the downscaled quantity is supposed
            ! to go to exactly 0, then make sure it is indeed exactly 0 rather than
            ! roundoff-level different from 0.
            downscaled = 0._r8
         end if
      end if
    end function AdjustPrecipTracer

  end subroutine SetDownscaledTracers

  !-----------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Deallocate memory associated with this instance
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    deallocate(this%forc_q_not_downscaled_grc)
    deallocate(this%forc_rain_not_downscaled_grc)
    deallocate(this%forc_snow_not_downscaled_grc)
    deallocate(this%forc_q_downscaled_col)
    deallocate(this%forc_flood_grc)
    deallocate(this%forc_rain_downscaled_col)
    deallocate(this%forc_snow_downscaled_col)
    deallocate(this%rain_to_snow_conversion_col)
    deallocate(this%snow_to_rain_conversion_col)

  end subroutine Clean

end module Wateratm2lndType
