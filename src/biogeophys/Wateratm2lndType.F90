module Wateratm2lndType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water atm2lnd variables that apply to both bulk water
  ! and water tracers.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varctl     , only : iulog
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use WaterInfoBaseType, only : water_info_base_type
  use WaterIsotopesMod, only : WisoCompareBulkToTracer
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: wateratm2lnd_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: forc_q_not_downscaled_grc     (:)   => null() ! not downscaled atm specific humidity (kg/kg)              
     real(r8), pointer :: forc_rain_not_downscaled_grc  (:)   => null() ! not downscaled atm rain rate [mm/s]
     real(r8), pointer :: forc_snow_not_downscaled_grc  (:)   => null() ! not downscaled atm snow rate [mm/s]
     real(r8), pointer :: forc_q_downscaled_col         (:)   => null() ! downscaled atm specific humidity (kg/kg)                  
     real(r8), pointer :: forc_flood_grc                (:)   => null() ! rof flood (mm/s)
     real(r8), pointer :: forc_rain_downscaled_col      (:)   => null() ! downscaled atm rain rate [mm/s]                                                                                       
     real(r8), pointer :: forc_snow_downscaled_col      (:)   => null() ! downscaled atm snow rate [mm/s]                                                                                       
                             


   contains

     procedure          :: Init         
     procedure          :: Restart      
     procedure          :: TracerConsistencyCheck
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type wateratm2lnd_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info)

    class(wateratm2lnd_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds  
    class(water_info_base_type), intent(in), target :: info

    this%info => info

    call this%InitAllocate(bounds)

    call this%InitHistory(bounds)

    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
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
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value     
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate(this%forc_q_not_downscaled_grc     (begg:endg))        ; this%forc_q_not_downscaled_grc     (:)   = ival
    allocate(this%forc_rain_not_downscaled_grc  (begg:endg))        ; this%forc_rain_not_downscaled_grc  (:)   = ival
    allocate(this%forc_snow_not_downscaled_grc  (begg:endg))        ; this%forc_snow_not_downscaled_grc  (:)   = ival
    allocate(this%forc_q_downscaled_col         (begc:endc))        ; this%forc_q_downscaled_col         (:)   = ival
    allocate(this%forc_flood_grc                (begg:endg))        ; this%forc_flood_grc                (:)   = ival
    allocate(this%forc_rain_downscaled_col      (begc:endc))        ; this%forc_rain_downscaled_col      (:)   = ival
    allocate(this%forc_snow_downscaled_col      (begc:endc))        ; this%forc_snow_downscaled_col      (:)   = ival


  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_lch4
    use clm_varctl     , only : hist_wrtch4diag
    use clm_varpar     , only : nlevsno, nlevsoi
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal, no_snow_zero
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
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
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use column_varcon   , only : icol_shadewall, icol_road_perv
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use clm_varcon      , only : denice, denh2o, spval, sb, bdsno 
    use clm_varcon      , only : zlnd, tfrz, spval, pc, aquifer_water_baseline
    use clm_varctl      , only : fsurdat, iulog
    use clm_varctl        , only : use_bedrock
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type), intent(in) :: this
    type(bounds_type)     , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev,nlevs 
    real(r8)           :: maxslope, slopemax, minslope
    real(r8)           :: d, fd, dfdd, slope0,slopebeta
    real(r8) ,pointer  :: std (:)     
    logical            :: readvar 
    type(file_desc_t)  :: ncid        
    character(len=256) :: locfn       
    integer            :: nbedrock
    !-----------------------------------------------------------------------



  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod          , only : masterproc
    use clm_varcon       , only : denice, denh2o, pondmx, watmin, spval, nameg
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager , only : is_first_step
    use clm_varctl       , only : bound_h2osoi
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(wateratm2lnd_type), intent(in) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,j,nlevs
    logical  :: readvar
    real(r8) :: maxwatsat    ! maximum porosity    
    real(r8) :: excess       ! excess volumetric soil water
    real(r8) :: totwat       ! total soil water (mm)
    !------------------------------------------------------------------------


    call restartvar(ncid=ncid, flag=flag, varname=this%info%fname('qflx_floodg'), xtype=ncd_double, &
         dim1name='gridcell', &
         long_name=this%info%lname('flood water flux'), units='mm/s', &
         interpinic_flag='skip', readvar=readvar, data=this%forc_flood_grc)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, readvar=readvar, not restart: initialize flood to zero
       this%forc_flood_grc = 0._r8
    endif



  end subroutine Restart

  function TracerConsistencyCheck(this,bounds,tracer) result(wiso_inconsistency)
    !
    ! !DESCRIPTION:
    ! Check consistency of water tracer with that of bulk water
    !
    ! !ARGUMENTS:

    logical :: wiso_inconsistency  ! function result
    class(wateratm2lnd_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    class(wateratm2lnd_type), intent(in) :: tracer
    !
    ! !LOCAL VARIABLES:
    integer l
    !-----------------------------------------------------------------------

    wiso_inconsistency = .false.

    wiso_inconsistency = wiso_inconsistency .or. &
      WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
                              this%forc_q_not_downscaled_grc(bounds%begg:bounds%endg), &
                              tracer%forc_q_not_downscaled_grc(bounds%begg:bounds%endg), &
                              'forc_q_not_downscaled_grc')

    wiso_inconsistency = wiso_inconsistency .or. &
      WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
                              this%forc_rain_not_downscaled_grc(bounds%begg:bounds%endg), &
                              tracer%forc_rain_not_downscaled_grc(bounds%begg:bounds%endg), &
                              'forc_rain_not_downscaled_grc')

    wiso_inconsistency = wiso_inconsistency .or. &
      WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
                              this%forc_snow_not_downscaled_grc(bounds%begg:bounds%endg), &
                              tracer%forc_snow_not_downscaled_grc(bounds%begg:bounds%endg), &
                              'forc_snow_not_downscaled_grc')

    wiso_inconsistency = wiso_inconsistency .or. &
      WisoCompareBulkToTracer(bounds%begc, bounds%endc, &
                              this%forc_q_downscaled_col(bounds%begc:bounds%endc), &
                              tracer%forc_q_downscaled_col(bounds%begc:bounds%endc), &
                              'forc_q_downscaled_col')

    wiso_inconsistency = wiso_inconsistency .or. &
      WisoCompareBulkToTracer(bounds%begg, bounds%endg, &
                              this%forc_flood_grc(bounds%begg:bounds%endg), &
                              tracer%forc_flood_grc(bounds%begg:bounds%endg), &
                              'forc_flood_grc')

    wiso_inconsistency = wiso_inconsistency .or. &
      WisoCompareBulkToTracer(bounds%begc, bounds%endc, &
                              this%forc_rain_downscaled_col(bounds%begc:bounds%endc), &
                              tracer%forc_rain_downscaled_col(bounds%begc:bounds%endc), &
                              'forc_rain_downscaled_col')

    wiso_inconsistency = wiso_inconsistency .or. &
      WisoCompareBulkToTracer(bounds%begc, bounds%endc, &
                              this%forc_snow_downscaled_col(bounds%begc:bounds%endc), &
                              tracer%forc_snow_downscaled_col(bounds%begc:bounds%endc), &
                              'forc_snow_downscaled_col')

  end function TracerConsistencyCheck


end module Wateratm2lndType
