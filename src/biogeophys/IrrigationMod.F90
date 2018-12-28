module IrrigationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates irrigation flux.
  !
  ! Usage:
  !
  !   - Call CalcIrrigationNeeded in order to compute whether and how much irrigation is
  !     needed for the next call to CalcIrrigationFluxes. This should be called once per
  !     timestep.
  ! 
  !   - Call CalcIrrigationFluxes in order to calculate irrigation withdrawal and
  !     application fluxes. It is acceptable for this to be called earlier in the timestep
  !     than CalcIrrigationNeeded.
  !
  ! General design notes:
  !
  !   In principle, CalcIrrigationFluxes and CalcIrrigationNeeded could be combined into a
  !   single routine. Their separation is largely for historical reasons: In the past,
  !   irrigation depended on btran, and qflx_irrig is needed earlier in the driver loop
  !   than when btran becomes available. (And qflx_irrig is also used late in the driver
  !   loop - so it wouldn't work, for example, to calculate qflx_irrig after btran is
  !   computed, and then save it on the restart file for the next iteration of the driver
  !   loop: then the uses of qflx_irrig early and late in the driver loop would be
  !   inconsistent.)
  !
  !   Now that we no longer have a dependency on btran, we could call CalcIrrigationNeeded
  !   before the first time qflx_irrig is needed. Thus, there might be some advantage to
  !   combining CalcIrrigationFluxes and CalcIrrigationNeeded - or at least calling these
  !   two routines from the same place.  In particular: this separation of the irrigation
  !   calculation into two routines that are done at different times in the driver loop
  !   makes it harder and less desirable to nest the irrigation object within some other
  !   object: Doing so might make it harder to do the two separate steps at the right
  !   time, and would lead to less clarity about how these two steps are ordered with
  !   respect to the rest of the driver loop. So if we start trying to create a hierarchy
  !   of objects in CLM, we may want to rework this design. And we may want to rework it
  !   to keep things simpler anyway, even if we aren't nesting objects. Note that this
  !   rework would change behavior slightly, because irrigation would be applied in the
  !   same time step that CalcIrrigationNeeded first determines it's needed, rather than
  !   waiting until the following time step.
  !
  ! Design notes regarding where we calculate supply-based irrigation limitations:
  !
  !   We compute supply-based irrigation limitations differently for rivers (volr) vs.
  !   groundwater.
  !
  !   For volr-based limitation, we apply the limit in CalcIrrigationNeeded rather than in
  !   CalcIrrigationFluxes. This is to avoid problems that arise due to evolving volr -
  !   e.g., if we had calculated irrigation demand such that volr was just barely big
  !   enough, then we apply irrigation for some time steps, then we get an updated volr
  !   which is lower due to these irrigation withdrawals - if we applied the volr
  !   limitation in CalcIrrigationFluxes, this would impose an unnecessary limit on the
  !   irrigation to be applied for the rest of this day. It's hard to see a clean way
  !   around this problem, given that volr is updated in some but not all CLM time
  !   steps. So we impose the limit on irrigation deficit, rather than on the
  !   time step-by-time-step flux. The downside here is that it's somewhat more likely
  !   that irrigation would try to draw volr negative, if volr has been evolving for other
  !   reasons.
  !
  !   For groundwater-based limitation, in contrast, we can rely on having an up-to-date
  !   view of available groundwater, so the argument against having the volr limitation
  !   in CalcIrrigationFluxes doesn't apply, and that seems like a more robust place to
  !   do this limitation.
  !
  !   The key difference is that, for volr, we don't have an easy way to determine if this
  !   is a time step when we have an updated volr value, or if we still have volr from a
  !   few time steps ago (since ROF isn't coupled in every time step). So in that case,
  !   having the possibility that we'd exceed available volr seems like the lesser of
  !   evils.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type, get_proc_global
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use clm_instur       , only : irrig_method
  use pftconMod        , only : pftcon
  use clm_varctl       , only : iulog
  use clm_varcon       , only : isecspday, denh2o, spval, ispval, namep, namec, nameg
  use clm_varpar       , only : nlevsoi, nlevgrnd
  use clm_time_manager , only : get_step_size
  use SoilHydrologyMod , only : CalcIrrigWithdrawals
  use SoilHydrologyType, only : soilhydrology_type
  use SoilStateType    , only : soilstate_type
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  use WaterType        , only : water_type
  use WaterFluxBulkType, only : waterfluxbulk_type
  use WaterFluxType    , only : waterflux_type
  use WaterStateBulkType, only : waterstatebulk_type
  use WaterStateType   , only : waterstate_type
  use WaterTracerUtils , only : CalcTracerFromBulk, CalcTracerFromBulkFixedRatio
  use GridcellType     , only : grc
  use ColumnType       , only : col                
  use PatchType        , only : patch                
  use subgridAveMod    , only : p2c, c2g
  use filterColMod     , only : filter_col_type, col_filter_from_logical_array
  !
  implicit none
  private

  ! !PUBLIC TYPES:
  
  ! This type is public (and its components are public, too) to aid unit testing
  type, public :: irrigation_params_type
     ! Minimum LAI for irrigation
     real(r8) :: irrig_min_lai

     ! Time of day to check whether we need irrigation, seconds (0 = midnight). 
     ! We start applying the irrigation in the time step FOLLOWING this time, 
     ! since we won't begin irrigating until the next call to CalcIrrigationFluxes
     integer  :: irrig_start_time

     ! Desired amount of time to irrigate per day (sec). Actual time may 
     ! differ if this is not a multiple of dtime. Irrigation won't work properly 
     ! if dtime > secsperday
     integer  :: irrig_length

     ! Target soil matric potential for irrigation (mm)
     !
     ! When we irrigate, we aim to bring the total soil moisture in the top (irrig_depth)
     ! m of soil up to this level.
     !
     ! (Note: When we convert this to a relative saturation, we ensure that the relative
     ! saturation target is bounded by [0,1].)
     real(r8) :: irrig_target_smp

     ! Soil depth to which we measure for irrigation (m)
     real(r8) :: irrig_depth

     ! Determines soil moisture threshold at which we irrigate. If
     ! h2osoi_liq_wilting_point is the soil moisture level at wilting point and
     ! h2osoi_liq_target is the soil moisture level at the target irrigation level (given
     ! by irrig_target_smp), then the threshold at which we irrigate is
     !     h2osoi_liq_wilting_point +
     !          irrig_threshold_fraction*(h2osoi_liq_target - h2osoi_liq_wilting_point)
     ! A value of 1 means that we irrigate whenever soil moisture falls below the target
     ! A value of 0 means that we only irrigate when soil moisture falls below the
     ! wilting point
     real(r8) :: irrig_threshold_fraction

     ! Threshold for river water volume below which irrigation is shut off, if
     ! limit_irrigation is .true. (fraction of available river water). A threshold of 0
     ! means allow all river water to be used; a threshold of 0.1 means allow 90% of the
     ! river volume to be used; etc.
     real(r8) :: irrig_river_volume_threshold

     ! Whether irrigation is limited based on river storage. This only applies if ROF is
     ! enabled (i.e., rof_prognostic is .true.) - otherwise we don't limit irrigation,
     ! regardless of the value of this flag.
     logical :: limit_irrigation_if_rof_enabled

     ! use groundwater supply for irrigation (in addition to surface water)
     logical :: use_groundwater_irrigation

     ! Irrigation method to use if not set on the surface dataset
     ! (This won't be needed once we can rely on irrig_method always being on the surface
     ! dataset, but it is useful until then - particularly for testing. See also
     ! https://github.com/ESCOMP/ctsm/issues/565.)
     integer :: irrig_method_default

  end type irrigation_params_type

  type, public :: irrigation_type
     private

     ! Private data members; set in initialization:
     type(irrigation_params_type) :: params
     integer :: dtime                ! land model time step (sec)
     integer :: irrig_nsteps_per_day ! number of time steps per day in which we irrigate
     real(r8), pointer :: relsat_wilting_point_col(:,:) ! relative saturation at which smp = wilting point [col, nlevsoi]
     real(r8), pointer :: relsat_target_col(:,:)        ! relative saturation at which smp is at the irrigation target [col, nlevsoi]
     integer , pointer :: irrig_method_patch          (:) ! patch irrigation application method

     ! Private data members; time-varying:
     real(r8), pointer :: sfc_irrig_rate_patch        (:) ! current irrigation rate from surface water [mm/s]
     real(r8), pointer :: irrig_rate_demand_patch     (:) ! current irrigation rate, neglecting surface water source limitation [mm/s]
     integer , pointer :: n_irrig_steps_left_patch    (:) ! number of time steps for which we still need to irrigate today (if 0, ignore)
     real(r8), pointer :: qflx_irrig_demand_patch     (:) ! irrigation flux neglecting surface water source limitation [mm/s]

   contains
     ! Public routines
     ! COMPILER_BUG(wjs, 2014-10-15, pgi 14.7) Add an "Irrigation" prefix to some  generic routines like "Init"
     ! (without this workaround, pgi compilation fails in restFileMod)
     procedure, public :: Init => IrrigationInit
     procedure, public :: Restart
     procedure, public :: CalcIrrigationFluxes
     procedure, public :: CalcIrrigationNeeded
     procedure, public :: UseGroundwaterIrrigation ! Returns true if groundwater irrigation enabled
     procedure, public :: Clean => IrrigationClean ! deallocate memory

     ! Public simply to support unit testing; should not be used from CLM code
     procedure, public :: InitForTesting ! version of Init meant for unit testing
     procedure, public :: RelsatToH2osoi ! convert from relative saturation to kg/m2 water for a single column and layer
     procedure, public :: WrapCalcIrrigWithdrawals ! Wrap the call to CalcIrrigWithdrawals

     ! Private routines
     procedure, private :: ReadNamelist
     procedure, private :: CheckNamelistValidity   ! Check for validity of input parameters
     procedure, private :: InitAllocate => IrrigationInitAllocate
     procedure, private :: InitHistory => IrrigationInitHistory
     procedure, private :: InitCold => IrrigationInitCold
     procedure, private :: CalcBulkWithdrawals      ! calculate irrigation withdrawals for bulk water
     procedure, private :: CalcOneTracerWithdrawals ! calculate irrigation withdrawals for one water tracer
     procedure, private :: CalcTotalGWUnconIrrig    ! calculate total irrigation withdrawal flux from the unconfined aquifer, for either bulk or one water tracer
     procedure, private :: CalcApplicationFluxes    ! calculate irrigation application fluxes for either bulk water or a single tracer
     procedure, private :: CalcIrrigNstepsPerDay    ! given dtime, calculate irrig_nsteps_per_day
     procedure, private :: SetIrrigMethod           ! set irrig_method_patch based on surface dataset
     procedure, private :: PointNeedsCheckForIrrig  ! whether a given point needs to be checked for irrigation now
     procedure, private :: CalcDeficitVolrLimited   ! calculate deficit limited by river volume for each patch
  end type irrigation_type

  interface irrigation_params_type
     module procedure irrigation_params_constructor
  end interface irrigation_params_type

  ! Soil matric potential at wilting point (mm)
  !
  ! There is no reason to make this a tunable parameter, because the behavior it governs
  ! (the trigger for irrigation) can be tuned via other parameters.
  !
  ! TODO(wjs, 2016-09-08) It looks like there is other code in CLM that also uses an
  ! assumed wilting point (CNRootDynMod, maybe others). We should probably make this a
  ! shared parameter, e.g., in clm_varcon.
  real(r8), parameter, private :: wilting_point_smp = -150000._r8

  ! Conversion factors
  real(r8), parameter :: m3_over_km2_to_mm = 1.e-3_r8
     
  ! Irrigation methods
  integer, parameter, public  :: irrig_method_unset = 0
  ! Drip is defined here as irrigation applied directly to soil surface
  integer, parameter, public :: irrig_method_drip = 1
  ! Sprinkler is applied directly to canopy
  integer, parameter, public :: irrig_method_sprinkler = 2

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function irrigation_params_constructor(irrig_min_lai, &
       irrig_start_time, irrig_length, &
       irrig_target_smp, &
       irrig_depth, irrig_threshold_fraction, irrig_river_volume_threshold, &
       limit_irrigation_if_rof_enabled, use_groundwater_irrigation, &
       irrig_method_default) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Create an irrigation_params instance
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(irrigation_params_type) :: this  ! function result
    real(r8), intent(in) :: irrig_min_lai
    integer , intent(in) :: irrig_start_time
    integer , intent(in) :: irrig_length
    real(r8), intent(in) :: irrig_target_smp
    real(r8), intent(in) :: irrig_depth
    real(r8), intent(in) :: irrig_threshold_fraction
    real(r8), intent(in) :: irrig_river_volume_threshold
    logical , intent(in) :: limit_irrigation_if_rof_enabled
    logical , intent(in) :: use_groundwater_irrigation
    integer , intent(in) :: irrig_method_default
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'irrigation_params_constructor'
    !-----------------------------------------------------------------------
    
    this%irrig_min_lai = irrig_min_lai
    this%irrig_start_time = irrig_start_time
    this%irrig_length = irrig_length
    this%irrig_depth = irrig_depth
    this%irrig_target_smp = irrig_target_smp
    this%irrig_threshold_fraction = irrig_threshold_fraction
    this%irrig_river_volume_threshold = irrig_river_volume_threshold
    this%limit_irrigation_if_rof_enabled = limit_irrigation_if_rof_enabled
    this%use_groundwater_irrigation = use_groundwater_irrigation
    this%irrig_method_default = irrig_method_default

  end function irrigation_params_constructor


  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================
  
  !------------------------------------------------------------------------
  subroutine IrrigationInit(this, bounds, NLFilename, &
       soilstate_inst, soil_water_retention_curve, &
       use_aquifer_layer)
    use SoilStateType , only : soilstate_type

    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    character(len=*)       , intent(in)    :: NLFilename ! Namelist filename
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    logical                , intent(in)    :: use_aquifer_layer ! whether an aquifer layer is used in this run

    call this%ReadNamelist(NLFilename, use_aquifer_layer)
    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds, soilstate_inst, soil_water_retention_curve)
  end subroutine IrrigationInit

  !-----------------------------------------------------------------------
  subroutine InitForTesting(this, bounds, params, dtime, &
       relsat_wilting_point, relsat_target)
    !
    ! !DESCRIPTION:
    ! Does initialization needed for unit testing. Allows caller to prescribe values of
    ! some internal variables.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(irrigation_type)       , intent(inout) :: this
    type(bounds_type)            , intent(in)    :: bounds
    type(irrigation_params_type) , intent(in)    :: params
    integer                      , intent(in)    :: dtime ! model time step (sec)
    real(r8)                     , intent(in)    :: relsat_wilting_point( bounds%begc: , 1: ) ! relative saturation at which smp = irrig_wilting_point_smp [col, nlevsoi]
    real(r8)                     , intent(in)    :: relsat_target( bounds%begc: , 1: ) ! relative saturation at which smp = irrig_target_smp [col, nlevsoi]
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitForTesting'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(relsat_wilting_point) == (/bounds%endc, nlevsoi/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(relsat_target) == (/bounds%endc, nlevsoi/)), errMsg(sourcefile, __LINE__))

    call this%InitAllocate(bounds)
    this%params = params
    this%dtime = dtime
    call this%SetIrrigMethod(bounds)
    this%irrig_nsteps_per_day = this%CalcIrrigNstepsPerDay(dtime)
    this%relsat_wilting_point_col(:,:) = relsat_wilting_point(:,:)
    this%relsat_target_col(:,:) = relsat_target(:,:)

  end subroutine InitForTesting

  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename, use_aquifer_layer)
    !
    ! !DESCRIPTION:
    ! Read the irrigation namelist
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    logical, intent(in) :: use_aquifer_layer    ! whether an aquifer layer is used in this run
    !
    ! !LOCAL VARIABLES:

    ! temporary variables corresponding to the components of irrigation_params_type
    real(r8) :: irrig_min_lai
    integer  :: irrig_start_time
    integer  :: irrig_length
    real(r8) :: irrig_target_smp
    real(r8) :: irrig_depth
    real(r8) :: irrig_threshold_fraction
    real(r8) :: irrig_river_volume_threshold
    logical  :: limit_irrigation_if_rof_enabled
    logical  :: use_groundwater_irrigation
    character(len=64) :: irrig_method_default
    integer  :: irrig_method_default_int

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: nmlname = 'irrigation_inparm'

    character(len=*), parameter :: subname = 'ReadNamelist'
    !-----------------------------------------------------------------------

    namelist /irrigation_inparm/ irrig_min_lai, irrig_start_time, irrig_length, &
         irrig_target_smp, irrig_depth, irrig_threshold_fraction, &
         irrig_river_volume_threshold, limit_irrigation_if_rof_enabled, &
         use_groundwater_irrigation, irrig_method_default
    
    ! Initialize options to garbage defaults, forcing all to be specified explicitly in
    ! order to get reasonable results
    irrig_min_lai = nan
    irrig_start_time = 0
    irrig_length = 0
    irrig_target_smp = nan
    irrig_depth = nan
    irrig_threshold_fraction = nan
    irrig_river_volume_threshold = nan
    limit_irrigation_if_rof_enabled = .false.
    use_groundwater_irrigation = .false.
    irrig_method_default = ' '

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=irrigation_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast(irrig_min_lai, mpicom)
    call shr_mpi_bcast(irrig_start_time, mpicom)
    call shr_mpi_bcast(irrig_length, mpicom)
    call shr_mpi_bcast(irrig_target_smp, mpicom)
    call shr_mpi_bcast(irrig_depth, mpicom)
    call shr_mpi_bcast(irrig_threshold_fraction, mpicom)
    call shr_mpi_bcast(irrig_river_volume_threshold, mpicom)
    call shr_mpi_bcast(limit_irrigation_if_rof_enabled, mpicom)
    call shr_mpi_bcast(use_groundwater_irrigation, mpicom)
    call shr_mpi_bcast(irrig_method_default, mpicom)

    call translate_irrig_method_default

    this%params = irrigation_params_type( &
         irrig_min_lai = irrig_min_lai, &
         irrig_start_time = irrig_start_time, &
         irrig_length = irrig_length, &
         irrig_target_smp = irrig_target_smp, &
         irrig_depth = irrig_depth, &
         irrig_threshold_fraction = irrig_threshold_fraction, &
         irrig_river_volume_threshold = irrig_river_volume_threshold, &
         limit_irrigation_if_rof_enabled = limit_irrigation_if_rof_enabled, &
         use_groundwater_irrigation = use_groundwater_irrigation, &
         irrig_method_default = irrig_method_default_int)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       ! Write settings one-by-one rather than with a nml write because
       ! irrig_river_volume_threshold may be NaN
       write(iulog,*) 'irrig_min_lai = ', irrig_min_lai
       write(iulog,*) 'irrig_start_time = ', irrig_start_time
       write(iulog,*) 'irrig_length = ', irrig_length
       write(iulog,*) 'irrig_target_smp = ', irrig_target_smp
       write(iulog,*) 'irrig_depth = ', irrig_depth
       write(iulog,*) 'irrig_threshold_fraction = ', irrig_threshold_fraction
       write(iulog,*) 'limit_irrigation_if_rof_enabled = ', limit_irrigation_if_rof_enabled
       if (limit_irrigation_if_rof_enabled) then
          write(iulog,*) 'irrig_river_volume_threshold = ', irrig_river_volume_threshold
       end if
       write(iulog,*) 'use_groundwater_irrigation = ', use_groundwater_irrigation
       write(iulog,*) 'irrig_method_default = ', irrig_method_default
       write(iulog,*) ' '

       call this%CheckNamelistValidity(use_aquifer_layer)
    end if

  contains
    subroutine translate_irrig_method_default
      select case (irrig_method_default)
      case ('drip')
         irrig_method_default_int = irrig_method_drip
      case ('sprinkler')
         irrig_method_default_int = irrig_method_sprinkler
      case default
         write(iulog,*) 'ERROR: unknown irrig_method_default: ', trim(irrig_method_default)
         call endrun('Unknown irrig_method_default')
      end select
    end subroutine translate_irrig_method_default

  end subroutine ReadNamelist

  !-----------------------------------------------------------------------
  subroutine CheckNamelistValidity(this, use_aquifer_layer)
    !
    ! !DESCRIPTION:
    ! Check for validity of input parameters.
    !
    ! Assumes that the inputs have already been packed into 'this%params'.
    !
    ! Only needs to be called by the master task, since parameters are the same for all
    ! tasks.
    !
    ! !ARGUMENTS:
    class(irrigation_type), intent(in) :: this
    logical, intent(in) :: use_aquifer_layer    ! whether an aquifer layer is used in this run
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CheckNamelistValidity'
    !-----------------------------------------------------------------------

    associate( &
         irrig_min_lai => this%params%irrig_min_lai, &
         irrig_start_time => this%params%irrig_start_time, &
         irrig_length => this%params%irrig_length, &
         irrig_target_smp => this%params%irrig_target_smp, &
         irrig_depth => this%params%irrig_depth, &
         irrig_threshold_fraction => this%params%irrig_threshold_fraction, &
         irrig_river_volume_threshold => this%params%irrig_river_volume_threshold, &
         use_groundwater_irrigation => this%params%use_groundwater_irrigation, &
         limit_irrigation_if_rof_enabled => this%params%limit_irrigation_if_rof_enabled)

    if (irrig_min_lai < 0._r8) then
       write(iulog,*) ' ERROR: irrig_min_lai must be >= 0'
       write(iulog,*) ' irrig_min_lai = ', irrig_min_lai
       call endrun(msg=' ERROR: irrig_min_lai must be >= 0 ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_start_time < 0 .or. irrig_start_time >= isecspday) then
       write(iulog,*) ' ERROR: irrig_start_time must be >= 0 and < ', isecspday
       write(iulog,*) ' irrig_start_time = ', irrig_start_time
       call endrun(msg=' ERROR: irrig_start_time out of bounds ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_length <= 0 .or. irrig_length > isecspday) then
       write(iulog,*) ' ERROR: irrig_length must be > 0 and <= ', isecspday
       write(iulog,*) ' irrig_length = ', irrig_length
       call endrun(msg=' ERROR: irrig_length out of bounds ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_target_smp >= 0._r8) then
       write(iulog,*) ' ERROR: irrig_target_smp must be negative'
       write(iulog,*) ' irrig_target_smp = ', irrig_target_smp
       call endrun(msg=' ERROR: irrig_target_smp must be negative ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_target_smp < wilting_point_smp) then
       write(iulog,*) ' ERROR: irrig_target_smp must be >= wilting_point_smp'
       write(iulog,*) ' irrig_target_smp (from namelist) = ', irrig_target_smp
       write(iulog,*) ' wilting_point_smp (hard-coded) = ', wilting_point_smp
       call endrun(msg=' ERROR: irrig_target_smp must be >= wilting_point_smp ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_depth < 0._r8) then
       write(iulog,*) ' ERROR: irrig_depth must be > 0'
       write(iulog,*) ' irrig_depth = ', irrig_depth
       call endrun(msg=' ERROR: irrig_depth must be > 0 ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_threshold_fraction < 0._r8 .or. irrig_threshold_fraction > 1._r8) then
       write(iulog,*) ' ERROR: irrig_threshold_fraction must be between 0 and 1'
       write(iulog,*) ' irrig_threshold_fraction = ', irrig_threshold_fraction
       call endrun(msg=' ERROR: irrig_threshold_fraction must be between 0 and 1 ' // &
            errMsg(sourcefile, __LINE__))
    end if

    if (limit_irrigation_if_rof_enabled) then
       if (irrig_river_volume_threshold < 0._r8 .or. irrig_river_volume_threshold > 1._r8) then
          write(iulog,*) ' ERROR: irrig_river_volume_threshold must be between 0 and 1'
          write(iulog,*) ' irrig_river_volume_threshold = ', irrig_river_volume_threshold
          call endrun(msg=' ERROR: irrig_river_volume_threshold must be between 0 and 1 ' // &
               errMsg(sourcefile, __LINE__))
       end if
    end if

    if (use_groundwater_irrigation .and. .not. limit_irrigation_if_rof_enabled) then
       write(iulog,*) ' ERROR: use_groundwater_irrigation only makes sense if limit_irrigation_if_rof_enabled is set.'
       write(iulog,*) '(If limit_irrigation_if_rof_enabled is .false., then groundwater extraction will never be invoked.)'
       call endrun(msg=' ERROR: use_groundwater_irrigation only makes sense if limit_irrigation_if_rof_enabled is set' // &
            errMsg(sourcefile, __LINE__))
    end if

    if (use_aquifer_layer .and. use_groundwater_irrigation) then
          write(iulog,*) ' ERROR: use_groundwater_irrigation and use_aquifer_layer may not be used simultaneously'
          call endrun(msg=' ERROR: use_groundwater_irrigation and use_aquifer_layer cannot both be set to true' // &
               errMsg(sourcefile, __LINE__))
    end if

    end associate

  end subroutine CheckNamelistValidity


  !-----------------------------------------------------------------------
  subroutine IrrigationInitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize irrigation data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%qflx_irrig_demand_patch (begp:endp))              ; this%qflx_irrig_demand_patch  (:)   = nan
    allocate(this%relsat_wilting_point_col    (begc:endc,nlevsoi))  ; this%relsat_wilting_point_col     (:,:) = nan
    allocate(this%relsat_target_col           (begc:endc,nlevsoi))  ; this%relsat_target_col            (:,:) = nan
    allocate(this%sfc_irrig_rate_patch        (begp:endp))          ; this%sfc_irrig_rate_patch         (:)   = nan
    allocate(this%irrig_rate_demand_patch     (begp:endp))          ; this%irrig_rate_demand_patch      (:)   = nan
    allocate(this%irrig_method_patch          (begp:endp))          ; this%irrig_method_patch           (:)   = ispval
    allocate(this%n_irrig_steps_left_patch    (begp:endp))          ; this%n_irrig_steps_left_patch     (:)   = 0

  end subroutine IrrigationInitAllocate

  !-----------------------------------------------------------------------
  subroutine IrrigationInitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize irrigation history fields
    !
    ! !USES:
    use histFileMod  , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%qflx_irrig_demand_patch(begp:endp) = spval
    call hist_addfld1d (fname='QIRRIG_DEMAND', units='mm/s', &
         avgflag='A', long_name='irrigation demand', &
         ptr_patch=this%qflx_irrig_demand_patch, default='inactive')

  end subroutine IrrigationInitHistory

  !-----------------------------------------------------------------------
  subroutine IrrigationInitCold(this, bounds, soilstate_inst, soil_water_retention_curve)
    !
    ! !DESCRIPTION:
    ! Do cold-start initialization for irrigation data structure
    !
    ! !USES:
    use SoilStateType , only : soilstate_type
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    integer :: c ! col index
    integer :: j ! level index

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    do j = 1, nlevsoi
       do c = bounds%begc, bounds%endc
          call soil_water_retention_curve%soil_suction_inverse( &
               c = c, &
               j = j, &
               smp_target = wilting_point_smp, &
               soilstate_inst = soilstate_inst, &
               s_target = this%relsat_wilting_point_col(c,j))

          call soil_water_retention_curve%soil_suction_inverse( &
               c = c, &
               j = j, &
               smp_target = this%params%irrig_target_smp, &
               soilstate_inst = soilstate_inst, &
               s_target = this%relsat_target_col(c,j))

          ! Make sure relative saturation targets are bounded by [0,1]
          !
          ! NOTE(wjs, 2016-11-17) These targets can be > 1 if smp_target is too small of
          ! a negative value; we want to force the target to a relative saturation of 1
          ! in that case. I don't see how these targets could end up negative, though I
          ! had a note that I saw negative values at one point. In practice, for
          ! reasonable parameter values, it seems that these min and max functions have
          ! no effect.
          this%relsat_wilting_point_col(c,j) = min(this%relsat_wilting_point_col(c,j), 1._r8)
          this%relsat_wilting_point_col(c,j) = max(this%relsat_wilting_point_col(c,j), 0._r8)
          this%relsat_target_col(c,j) = min(this%relsat_target_col(c,j), 1._r8)
          this%relsat_target_col(c,j) = max(this%relsat_target_col(c,j), 0._r8)
       end do
    end do

    call this%SetIrrigMethod(bounds)
       
    this%dtime = get_step_size()
    this%irrig_nsteps_per_day = this%CalcIrrigNstepsPerDay(this%dtime)

    this%qflx_irrig_demand_patch(bounds%begp:bounds%endp) = 0._r8
    
  end subroutine IrrigationInitCold

  !-----------------------------------------------------------------------
  pure function CalcIrrigNstepsPerDay(this, dtime) result(irrig_nsteps_per_day)
    !
    ! !DESCRIPTION:
    ! Given dtime (sec), determine number of irrigation steps per day
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: irrig_nsteps_per_day  ! function result
    class(irrigation_type) , intent(in) :: this
    integer                , intent(in) :: dtime ! model time step (sec)
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'CalcIrrigNstepsPerDay'
    !-----------------------------------------------------------------------
    
    irrig_nsteps_per_day = ((this%params%irrig_length + (dtime - 1))/dtime)  ! round up

  end function CalcIrrigNstepsPerDay

  !-----------------------------------------------------------------------
  subroutine SetIrrigMethod(this, bounds)
    !
    ! !DESCRIPTION:
    ! Set this%irrig_method_patch based on surface dataset
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p ! patch index
    integer :: g ! gridcell index
    integer :: m ! patch itype

    character(len=*), parameter :: subname = 'SetIrrigMethod'
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       m = patch%itype(p)
       if (m >= lbound(irrig_method, 2) .and. m <= ubound(irrig_method, 2) &
            .and. pftcon%irrigated(m) == 1._r8) then
          this%irrig_method_patch(p) = irrig_method(g,m)
          ! ensure irrig_method is valid; if not set, use drip irrigation
          if(irrig_method(g,m) == irrig_method_unset) then
             this%irrig_method_patch(p) = this%params%irrig_method_default
          else if (irrig_method(g,m) /= irrig_method_drip .and. irrig_method(g,m) /= irrig_method_sprinkler) then
             write(iulog,*) subname //' invalid irrigation method specified'
             call endrun(decomp_index=g, clmlevel=nameg, msg='bad irrig_method '// &
                  errMsg(sourcefile, __LINE__))
          end if
       else
          this%irrig_method_patch(p) = this%params%irrig_method_default
       end if
    end do

  end subroutine SetIrrigMethod


  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart of irrigation variables
    !
    ! !USES:
    use ncdio_pio        , only : file_desc_t, ncd_inqvdlen, ncd_double, ncd_int
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(irrigation_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    logical :: do_io
    integer :: dimlen       ! dimension length
    integer :: nump_global  ! total number of patchs, globally
    integer :: err_code     ! error code
    logical :: readvar      ! determine if variable is on initial file

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------
    
    ! Get expected total number of points, for later error checks
    call get_proc_global(np=nump_global)

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! BACKWARDS_COMPATIBILITY
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'n_irrig_steps_left', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='n_irrig_steps_left', xtype=ncd_int,  &
            dim1name='pft', &
            long_name='number of irrigation time steps left', units='#', &
            interpinic_flag='interp', readvar=readvar, data=this%n_irrig_steps_left_patch)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%n_irrig_steps_left_patch = 0
    end if

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! BACKWARDS_COMPATIBILITY
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'irrig_rate', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='irrig_rate', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='surface irrigation rate', units='mm/s', &
            interpinic_flag='interp', readvar=readvar, data=this%sfc_irrig_rate_patch)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%sfc_irrig_rate_patch = 0.0_r8
    end if

    ! BACKWARDS_COMPATIBILITY(wjs, 2016-11-23) To support older restart files without an
    ! irrig_rate_demand field, get irrig_rate_demand from irrig_rate. I'm abusing the
    ! capability to specify multiple variable names here (even though irrig_rate isn't
    ! really an old version of irrig_rate_demand), rather than having code like 'if
    ! (.not. readvar) then ...', because the latter doesn't work when using init_interp.
    call restartvar(ncid=ncid, flag=flag, varname='irrig_rate_demand:irrig_rate', &
         xtype=ncd_double,  &
         dim1name='pft', &
         long_name='irrigation rate demand, neglecting surface water source limitation', &
         units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%irrig_rate_demand_patch)
  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine IrrigationClean(this)
    !
    ! !DESCRIPTION:
    ! Deallocate memory
    !
    ! !ARGUMENTS:
    class(irrigation_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------
    
    deallocate(this%qflx_irrig_demand_patch)
    deallocate(this%relsat_wilting_point_col)
    deallocate(this%relsat_target_col)
    deallocate(this%sfc_irrig_rate_patch)
    deallocate(this%irrig_rate_demand_patch)
    deallocate(this%irrig_method_patch)
    deallocate(this%n_irrig_steps_left_patch)

  end subroutine IrrigationClean


  ! ========================================================================
  ! Science routines
  ! ========================================================================
  
  !-----------------------------------------------------------------------
  subroutine CalcIrrigationFluxes(this, bounds, num_soilc, &
       filter_soilc, num_soilp, filter_soilp, &
       soilhydrology_inst, soilstate_inst, &
       water_inst)
    !
    ! !DESCRIPTION:
    ! Apply the irrigation computed by CalcIrrigationNeeded in order to set various fluxes
    !
    ! Sets irrigation withdrawal and application fluxes in the waterflux components of
    ! water_inst, for both bulk and tracers:
    ! - qflx_sfc_irrig_col
    ! - qflx_gw_uncon_irrig_lyr_col
    ! - qflx_gw_uncon_irrig_col
    ! - qflx_gw_con_irrig_col
    ! - qflx_irrig_drip_patch
    ! - qflx_irrig_sprinkler_patch
    !
    ! Should be called once, AND ONLY ONCE, per time step.
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    ! number of points in filter_soilc
    integer, intent(in) :: num_soilc
    ! column filter for soil
    integer, intent(in) :: filter_soilc(:)
    ! number of points in filter_soilp
    integer, intent(in) :: num_soilp
    ! patch filter for soil
    integer, intent(in) :: filter_soilp(:)
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(water_type)         , intent(inout) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,fp  ! patch indices
    integer  :: c,fc  ! column indices
    integer  :: g     ! grid cell index
    integer  :: j     ! level index
    integer  :: i     ! tracer index
    real(r8) :: qflx_sfc_irrig_bulk_patch(bounds%begp:bounds%endp)      ! patch surface irrigation flux for bulk water (mm H2O/s)
    real(r8) :: qflx_gw_demand_bulk_patch(bounds%begp:bounds%endp) ! patch amount of irrigation groundwater demand for bulk water (mm H2O/s)
    real(r8) :: qflx_gw_demand_bulk_col(bounds%begc:bounds%endc) ! col amount of irrigation groundwater demand for bulk water (mm H2O/s)
    logical  :: is_bulk

    character(len=*), parameter :: subname = 'CalcIrrigationFluxes'

    !-----------------------------------------------------------------------

    ! This should be called exactly once per time step, so that the counter decrease
    ! works correctly.

    ! Note that these associated variables refer to the bulk instance
    associate( &
         begp => bounds%begp, &
         endp => bounds%endp, &
         begc => bounds%begc, &
         endc => bounds%endc  &
         )

    call this%CalcBulkWithdrawals(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         soilhydrology_inst, soilstate_inst, water_inst%waterfluxbulk_inst, &
         qflx_sfc_irrig_bulk_patch = qflx_sfc_irrig_bulk_patch(begp:endp), &
         qflx_gw_demand_bulk_patch = qflx_gw_demand_bulk_patch(begp:endp), &
         qflx_gw_demand_bulk_col = qflx_gw_demand_bulk_col(begc:endc))

    do i = water_inst%tracers_beg, water_inst%tracers_end
       call this%CalcOneTracerWithdrawals(bounds, num_soilc, filter_soilc, &
            waterstatebulk_inst = water_inst%waterstatebulk_inst, &
            waterfluxbulk_inst = water_inst%waterfluxbulk_inst, &
            waterstate_tracer_inst = water_inst%bulk_and_tracers(i)%waterstate_inst, &
            waterflux_tracer_inst = water_inst%bulk_and_tracers(i)%waterflux_inst)
    end do

    if (this%params%use_groundwater_irrigation) then
       do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
          call this%CalcTotalGWUnconIrrig(num_soilc, filter_soilc, &
               water_inst%bulk_and_tracers(i)%waterflux_inst)
       end do
    end if

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       is_bulk = (i == water_inst%i_bulk)
       call this%CalcApplicationFluxes(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            waterflux_inst = water_inst%bulk_and_tracers(i)%waterflux_inst, &
            is_bulk = is_bulk, &
            qflx_sfc_irrig_bulk_patch = qflx_sfc_irrig_bulk_patch(begp:endp), &
            qflx_sfc_irrig_bulk_col = water_inst%waterfluxbulk_inst%qflx_sfc_irrig_col(begc:endc), &
            qflx_gw_demand_bulk_patch = qflx_gw_demand_bulk_patch(begp:endp), &
            qflx_gw_demand_bulk_col = qflx_gw_demand_bulk_col(begc:endc))
    end do

    end associate

  end subroutine CalcIrrigationFluxes

  !-----------------------------------------------------------------------
  subroutine CalcBulkWithdrawals(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilhydrology_inst, soilstate_inst, waterfluxbulk_inst, &
       qflx_sfc_irrig_bulk_patch, qflx_gw_demand_bulk_patch, qflx_gw_demand_bulk_col)
    !
    ! !DESCRIPTION:
    ! Calculate irrigation withdrawals for bulk water
    !
    ! Sets / updates the following variables:
    ! - qflx_sfc_irrig_bulk_patch
    ! - qflx_gw_demand_bulk_patch
    ! - qflx_gw_demand_bulk_col
    ! - this%qflx_irrig_demand_patch
    ! - this%n_irrig_steps_left_patch
    ! - waterfluxbulk_inst%qflx_sfc_irrig_col
    ! - waterfluxbulk_inst%qflx_gw_uncon_irrig_lyr_col
    ! - waterfluxbulk_inst%qflx_gw_con_irrig_col
    !
    ! !ARGUMENTS:
    class(irrigation_type)   , intent(inout) :: this
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of points in filter_soilc
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil
    integer                  , intent(in)    :: num_soilp       ! number of points in filter_soilp
    integer                  , intent(in)    :: filter_soilp(:) ! patch filter for soil
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(waterfluxbulk_type) , intent(inout) :: waterfluxbulk_inst
    real(r8)                 , intent(inout) :: qflx_sfc_irrig_bulk_patch( bounds%begp: ) ! patch surface irrigation flux for bulk water (mm H2O/s)
    real(r8)                 , intent(inout) :: qflx_gw_demand_bulk_patch( bounds%begp: ) ! patch amount of irrigation groundwater demand for bulk water (mm H2O/s)
    real(r8)                 , intent(inout) :: qflx_gw_demand_bulk_col( bounds%begc: )   ! col amount of irrigation groundwater demand for bulk water (mm H2O/s)

    !
    ! !LOCAL VARIABLES:
    integer  :: p,fp  ! patch indices

    character(len=*), parameter :: subname = 'CalcBulkWithdrawals'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(qflx_sfc_irrig_bulk_patch) == [bounds%endp]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_gw_demand_bulk_patch) == [bounds%endp]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_gw_demand_bulk_col) == [bounds%endc]), errMsg(sourcefile, __LINE__))

    associate( &
         begp => bounds%begp, &
         endp => bounds%endp, &
         begc => bounds%begc, &
         endc => bounds%endc, &

         qflx_sfc_irrig_col          => waterfluxbulk_inst%qflx_sfc_irrig_col          , & ! Output: [real(r8) (:)] col irrigation flux (mm H2O/s)
         qflx_gw_uncon_irrig_lyr_col => waterfluxbulk_inst%qflx_gw_uncon_irrig_lyr_col , & ! Output: [real(r8) (:,:) ] unconfined groundwater irrigation flux, separated by layer (mm H2O/s)
         qflx_gw_con_irrig_col       => waterfluxbulk_inst%qflx_gw_con_irrig_col         & ! Output: [real(r8) (:)] col confined groundwater irrigation flux (mm H2O/s)
         )

    do fp = 1, num_soilp
       p = filter_soilp(fp)
       
       if (this%n_irrig_steps_left_patch(p) > 0) then
          qflx_sfc_irrig_bulk_patch(p)     = this%sfc_irrig_rate_patch(p)
          this%qflx_irrig_demand_patch(p)  = this%irrig_rate_demand_patch(p)
          qflx_gw_demand_bulk_patch(p)     = &
               this%qflx_irrig_demand_patch(p) - qflx_sfc_irrig_bulk_patch(p)
          
          this%n_irrig_steps_left_patch(p) = this%n_irrig_steps_left_patch(p) - 1
       else
          qflx_sfc_irrig_bulk_patch(p)    = 0._r8
          this%qflx_irrig_demand_patch(p) = 0._r8
          qflx_gw_demand_bulk_patch(p)    = 0._r8
       end if
    end do

    call p2c (bounds, num_soilc, filter_soilc, &
         patcharr = qflx_sfc_irrig_bulk_patch(begp:endp), &
         colarr = qflx_sfc_irrig_col(begc:endc))

    call p2c (bounds, num_soilc, filter_soilc, &
         patcharr = qflx_gw_demand_bulk_patch(begp:endp), &
         colarr = qflx_gw_demand_bulk_col(begc:endc))

    if (this%params%use_groundwater_irrigation) then
       ! Supply as much of the remaining deficit as possible from groundwater irrigation
       call this%WrapCalcIrrigWithdrawals(bounds, num_soilc, filter_soilc, &
            soilhydrology_inst, soilstate_inst, &
            qflx_gw_demand = qflx_gw_demand_bulk_col(begc:endc), &
            qflx_gw_uncon_irrig_lyr = qflx_gw_uncon_irrig_lyr_col(begc:endc, :), &
            qflx_gw_con_irrig = qflx_gw_con_irrig_col(begc:endc))
    end if

    end associate

  end subroutine CalcBulkWithdrawals

  !-----------------------------------------------------------------------
  subroutine CalcOneTracerWithdrawals(this, bounds, num_soilc, filter_soilc, &
       waterstatebulk_inst, waterfluxbulk_inst, &
       waterstate_tracer_inst, waterflux_tracer_inst)
    !
    ! !DESCRIPTION:
    ! Calculate irrigation withdrawals for one water tracer
    !
    ! Sets the following variables:
    ! - waterflux_tracer_inst%qflx_sfc_irrig_col
    ! - waterflux_tracer_inst%qflx_gw_uncon_irrig_lyr_col
    ! - waterflux_tracer_inst%qflx_gw_con_irrig_col
    !
    ! !ARGUMENTS:
    class(irrigation_type)    , intent(in)    :: this
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_soilc       ! number of points in filter_soilc
    integer                   , intent(in)    :: filter_soilc(:) ! column filter for soil
    type(waterstatebulk_type) , intent(in)    :: waterstatebulk_inst
    type(waterfluxbulk_type)  , intent(in)    :: waterfluxbulk_inst
    type(waterstate_type)     , intent(in)    :: waterstate_tracer_inst
    type(waterflux_type)      , intent(inout) :: waterflux_tracer_inst
    !
    ! !LOCAL VARIABLES:
    integer :: j  ! level index

    character(len=*), parameter :: subname = 'CalcOneTracerWithdrawals'
    !-----------------------------------------------------------------------

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc &
         )

    ! BUG(wjs, 2018-12-17, ESCOMP/ctsm#512) Eventually, use actual source ratio for
    ! qflx_sfc_irrig_col rather than assuming a fixed ratio
    call CalcTracerFromBulkFixedRatio( &
         bulk   = waterfluxbulk_inst%qflx_sfc_irrig_col(begc:endc), &
         ratio  = waterflux_tracer_inst%info%get_ratio(), &
         tracer = waterflux_tracer_inst%qflx_sfc_irrig_col(begc:endc))

    if (this%params%use_groundwater_irrigation) then

       call CalcTracerFromBulk( &
            lb            = begc, &
            num_pts       = num_soilc, &
            filter_pts    = filter_soilc, &
            bulk_source   = waterstatebulk_inst%wa_col(begc:endc), &
            bulk_val      = waterfluxbulk_inst%qflx_gw_con_irrig_col(begc:endc), &
            tracer_source = waterstate_tracer_inst%wa_col(begc:endc), &
            tracer_val    = waterflux_tracer_inst%qflx_gw_con_irrig_col(begc:endc))
       do j = 1, nlevsoi
          call CalcTracerFromBulk( &
               lb            = begc, &
               num_pts       = num_soilc, &
               filter_pts    = filter_soilc, &
               bulk_source   = waterstatebulk_inst%h2osoi_liq_col(begc:endc, j), &
               bulk_val      = waterfluxbulk_inst%qflx_gw_uncon_irrig_lyr_col(begc:endc, j), &
               tracer_source = waterstate_tracer_inst%h2osoi_liq_col(begc:endc, j), &
               tracer_val    = waterflux_tracer_inst%qflx_gw_uncon_irrig_lyr_col(begc:endc, j))
       end do

    end if

    end associate

  end subroutine CalcOneTracerWithdrawals

  !-----------------------------------------------------------------------
  subroutine CalcTotalGWUnconIrrig(this, num_soilc, filter_soilc, &
       waterflux_inst)
    !
    ! !DESCRIPTION:
    ! Calculate total irrigation withdrawal flux from the unconfined aquifer, for either bulk
    ! or one water tracer
    !
    ! Sets the following variables:
    ! - waterflux_inst%qflx_gw_uncon_irrig_col
    !
    ! !ARGUMENTS:
    class(irrigation_type)    , intent(in)    :: this
    integer                   , intent(in)    :: num_soilc       ! number of points in filter_soilc
    integer                   , intent(in)    :: filter_soilc(:) ! column filter for soil
    class(waterflux_type)     , intent(inout) :: waterflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c, fc  ! column indices
    integer :: j      ! level index

    character(len=*), parameter :: subname = 'CalcTotalGWUnconIrrig'
    !-----------------------------------------------------------------------

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       waterflux_inst%qflx_gw_uncon_irrig_col(c) = 0._r8
    end do
    do j = 1, nlevsoi
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          waterflux_inst%qflx_gw_uncon_irrig_col(c) = &
               waterflux_inst%qflx_gw_uncon_irrig_col(c) + &
               waterflux_inst%qflx_gw_uncon_irrig_lyr_col(c,j)
       end do
    end do

  end subroutine CalcTotalGWUnconIrrig

  !-----------------------------------------------------------------------
  subroutine CalcApplicationFluxes(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       waterflux_inst, is_bulk, &
       qflx_sfc_irrig_bulk_patch, qflx_sfc_irrig_bulk_col, &
       qflx_gw_demand_bulk_patch, qflx_gw_demand_bulk_col)
    !
    ! !DESCRIPTION:
    ! Calculate irrigation application fluxes for either bulk water or a single tracer
    !
    ! Sets:
    ! - waterflux_inst%qflx_irrig_drip_patch
    ! - waterflux_inst%qflx_irrig_sprinkler_patch
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(in)    :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of points in filter_soilc
    integer                , intent(in)    :: filter_soilc(:) ! column filter for soil
    integer                , intent(in)    :: num_soilp       ! number of points in filter_soilp
    integer                , intent(in)    :: filter_soilp(:) ! patch filter for soil
    class(waterflux_type)  , intent(inout) :: waterflux_inst
    logical                , intent(in)    :: is_bulk         ! whether this call is being made for bulk water (is_bulk true) or a tracer (is_bulk false)
    real(r8)               , intent(in)    :: qflx_sfc_irrig_bulk_patch( bounds%begp: ) ! patch surface irrigation flux for bulk water (mm H2O/s)
    real(r8)               , intent(in)    :: qflx_sfc_irrig_bulk_col( bounds%begc: )   ! col surface irrigation flux for bulk water (mm H2O/s)
    real(r8)               , intent(in)    :: qflx_gw_demand_bulk_patch( bounds%begp: ) ! patch amount of irrigation groundwater demand for bulk water (mm H2O/s)
    real(r8)               , intent(in)    :: qflx_gw_demand_bulk_col( bounds%begc: )   ! col amount of irrigation groundwater demand for bulk water (mm H2O/s)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c  ! column indices
    integer  :: fp, p  ! patch indices
    real(r8) :: qflx_gw_irrig_withdrawn_col(bounds%begc:bounds%endc) ! col total amount of irrigation withdrawn from groundwater (mm H2O/s)
    real(r8) :: qflx_sfc_irrig        ! amount of surface irrigation for a single patch (mm H2O/s)
    real(r8) :: qflx_gw_irrig_applied ! amount of groundwater irrigation for a single patch (mm H2O/s)
    real(r8) :: qflx_irrig_tot        ! total amount of irrigation for a single patch (mm H2O/s)

    character(len=*), parameter :: subname = 'CalcApplicationFluxes'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(qflx_sfc_irrig_bulk_patch) == [bounds%endp]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_sfc_irrig_bulk_col) == [bounds%endc]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_gw_demand_bulk_patch) == [bounds%endp]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_gw_demand_bulk_col) == [bounds%endc]), errMsg(sourcefile, __LINE__))

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       ! Note that qflx_gw_uncon_irrig_col and qflx_gw_con_irrig_col will remain at
       ! their cold start values of 0 if use_groundwater_irrigation is false, so this
       ! sum is safe to do in all cases.
       qflx_gw_irrig_withdrawn_col(c) = &
            waterflux_inst%qflx_gw_uncon_irrig_col(c) + &
            waterflux_inst%qflx_gw_con_irrig_col(c)
    end do

    do fp = 1, num_soilp
       p = filter_soilp(fp)
       c = patch%column(p)

       if (is_bulk) then
          qflx_sfc_irrig = qflx_sfc_irrig_bulk_patch(p)
       else
          if (qflx_sfc_irrig_bulk_col(c) > 0._r8) then
             ! We have col-level sfc irrig for each water tracer, but only have
             ! patch-level for the bulk. To get the patch-level tracer flux, we need
             ! to scale the column-level flux by the ratio of patch-to-column-level
             ! flux for the bulk.
             !
             ! Note that this is similar to the scaling that we do for groundwater
             ! irrigation below.
             qflx_sfc_irrig = waterflux_inst%qflx_sfc_irrig_col(c) * &
                  (qflx_sfc_irrig_bulk_patch(p) / qflx_sfc_irrig_bulk_col(c))
          else
             if (qflx_sfc_irrig_bulk_patch(p) > 0._r8 .or. &
                  waterflux_inst%qflx_sfc_irrig_col(c) > 0._r8) then
                write(iulog,*) 'If qflx_sfc_irrig_bulk_col <= 0, ' // &
                     'expect qflx_sfc_irrig_bulk_patch = waterflux_inst%qflx_sfc_irrig_col = 0'
                write(iulog,*) 'qflx_sfc_irrig_bulk_col = ', qflx_sfc_irrig_bulk_col(c)
                write(iulog,*) 'qflx_sfc_irrig_bulk_patch = ', qflx_sfc_irrig_bulk_patch(p)
                write(iulog,*) 'waterflux_inst%qflx_sfc_irrig_col = ', &
                     waterflux_inst%qflx_sfc_irrig_col(c)
                call endrun(decomp_index=p, clmlevel=namep, &
                     msg = 'If qflx_sfc_irrig_bulk_col <= 0, ' // &
                     'expect qflx_sfc_irrig_bulk_patch = waterflux_inst%qflx_sfc_irrig_col = 0', &
                     additional_msg = errMsg(sourcefile, __LINE__))
             end if

             qflx_sfc_irrig = 0._r8
          end if
       end if

       if (qflx_gw_demand_bulk_col(c) > 0._r8) then
          ! Distribute the withdrawn groundwater irrigation to each patch according to
          ! its demand relative to the total column demand.
          !
          ! Note that this expression could be reformulated to:
          !   qflx_gw_irrig_applied = qflx_gw_demand_bulk_patch(p) * &
          !        (qflx_gw_irrig_withdrawn_col(c) / qflx_gw_demand_bulk_col(c))
          !
          ! It may be more intuitive that the above is correct, at least for bulk water:
          ! this shows that we're down-weighting each patch's demand evenly, based on the
          ! ration of withdrawn water to demand in the given column.
          !
          ! Note that this is similar to the scaling that we do for surface irrigation
          ! above.
          qflx_gw_irrig_applied = qflx_gw_irrig_withdrawn_col(c) * &
               (qflx_gw_demand_bulk_patch(p) / qflx_gw_demand_bulk_col(c))
       else
          if (qflx_gw_demand_bulk_patch(p) > 0._r8 .or. &
               qflx_gw_irrig_withdrawn_col(c) > 0._r8) then
             write(iulog,*) 'If qflx_gw_demand_bulk_col <= 0, expect qflx_gw_demand_bulk_patch = qflx_gw_irrig_withdrawn_col = 0'
             write(iulog,*) 'qflx_gw_demand_bulk_col = ', qflx_gw_demand_bulk_col(c)
             write(iulog,*) 'qflx_gw_demand_bulk_patch = ', qflx_gw_demand_bulk_patch(p)
             write(iulog,*) 'qflx_gw_irrig_withdrawn_col = ', qflx_gw_irrig_withdrawn_col(c)
             call endrun(decomp_index=p, clmlevel=namep, &
                  msg = 'If qflx_gw_demand_bulk_col <= 0, expect qflx_gw_demand_bulk_patch = qflx_gw_irrig_withdrawn_col = 0', &
                  additional_msg = errMsg(sourcefile, __LINE__))
          end if

          qflx_gw_irrig_applied = 0._r8
       end if

       qflx_irrig_tot = qflx_sfc_irrig + qflx_gw_irrig_applied

       ! Set drip/sprinkler irrigation based on irrigation method from input data
       waterflux_inst%qflx_irrig_drip_patch(p)      = 0._r8
       waterflux_inst%qflx_irrig_sprinkler_patch(p) = 0._r8

       if(this%irrig_method_patch(p) == irrig_method_drip) then
          waterflux_inst%qflx_irrig_drip_patch(p)      = qflx_irrig_tot
       else if(this%irrig_method_patch(p) == irrig_method_sprinkler) then
          waterflux_inst%qflx_irrig_sprinkler_patch(p) = qflx_irrig_tot
       else
          call endrun(msg=' ERROR: irrig_method_patch set to invalid value ' // &
               errMsg(sourcefile, __LINE__))
       endif

    end do

  end subroutine CalcApplicationFluxes

  !-----------------------------------------------------------------------
  subroutine WrapCalcIrrigWithdrawals(this, bounds, num_soilc, filter_soilc, &
       soilhydrology_inst, soilstate_inst, &
       qflx_gw_demand, &
       qflx_gw_uncon_irrig_lyr, &
       qflx_gw_con_irrig)
    !
    ! !DESCRIPTION:
    ! Wrap the call to the external subroutine, CalcIrrigWithdrawals
    !
    ! The purpose of this wrapping is to allow unit tests to override the behavior
    !
    ! !ARGUMENTS:
    class(irrigation_type)   , intent(in)    :: this
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of points in filter_soilc
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst

    real(r8) , intent(in)    :: qflx_gw_demand( bounds%begc: )              ! groundwater irrigation demand (mm H2O/s)
    real(r8) , intent(inout) :: qflx_gw_uncon_irrig_lyr( bounds%begc:, 1: ) ! unconfined aquifer groundwater irrigation withdrawal flux, separated by layer (mm H2O/s)
    real(r8) , intent(inout) :: qflx_gw_con_irrig( bounds%begc: )           ! total confined aquifer groundwater irrigation withdrawal flux (mm H2O/s)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'WrapCalcIrrigWithdrawals'
    !-----------------------------------------------------------------------

    call CalcIrrigWithdrawals( &
         bounds = bounds, &
         num_soilc = num_soilc, &
         filter_soilc = filter_soilc, &
         soilhydrology_inst = soilhydrology_inst, &
         soilstate_inst = soilstate_inst, &
         qflx_gw_demand = qflx_gw_demand(bounds%begc:bounds%endc), &
         qflx_gw_uncon_irrig_lyr = qflx_gw_uncon_irrig_lyr(bounds%begc:bounds%endc, :), &
         qflx_gw_con_irrig = qflx_gw_con_irrig(bounds%begc:bounds%endc))

  end subroutine WrapCalcIrrigWithdrawals

  !-----------------------------------------------------------------------
  subroutine CalcIrrigationNeeded(this, bounds, num_exposedvegp, filter_exposedvegp, &
       elai, t_soisno, eff_porosity, h2osoi_liq, volr, rof_prognostic)
    !
    ! !DESCRIPTION:
    ! Calculate whether and how much irrigation is needed for each column. However, this
    ! does NOT actually set the irrigation flux.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds

    ! number of points in filter_exposedvegp
    integer, intent(in) :: num_exposedvegp

    ! patch filter for non-snow-covered veg
    integer, intent(in) :: filter_exposedvegp(:)

    ! one-sided leaf area index with burying by snow [patch]
    real(r8), intent(in) :: elai( bounds%begp: )

    ! col soil temperature (K) [col, nlevgrnd] (note that this does NOT contain the snow levels)
    real(r8), intent(in) :: t_soisno( bounds%begc: , 1: )

    ! effective porosity (0 to 1) [col, nlevgrnd]
    real(r8), intent(in) :: eff_porosity( bounds%begc: , 1: )

    ! column liquid water (kg/m2) [col, nlevgrnd] (note that this does NOT contain the snow levels)
    real(r8), intent(in) :: h2osoi_liq( bounds%begc: , 1: )

    ! river water volume (m3) (ignored if rof_prognostic is .false.)
    real(r8), intent(in) :: volr( bounds%begg: )

    ! whether we're running with a prognostic ROF component; this is needed to determine
    ! whether we can limit irrigation based on river volume.
    logical, intent(in) :: rof_prognostic

    !
    ! !LOCAL VARIABLES:
    integer :: fp   ! patch filter index
    integer :: fc   ! column filter index
    integer :: p    ! patch index
    integer :: c    ! column index
    integer :: g    ! gridcell index
    integer :: j    ! level

    ! Filter for columns where we need to check for irrigation
    type(filter_col_type) :: check_for_irrig_col_filter

    ! target soil moisture for this layer [kg/m2]
    real(r8) :: h2osoi_liq_target

    ! soil moisture at wilting point for this layer [kg/m2]
    real(r8) :: h2osoi_liq_wilting_point

    ! Total of h2osoi down to the depth of irrigation in each column [kg/m2]
    real(r8) :: h2osoi_liq_tot(bounds%begc:bounds%endc)

    ! Total of h2osoi_liq_target down to the depth of irrigation in each column [kg/m2]
    real(r8) :: h2osoi_liq_target_tot(bounds%begc:bounds%endc)

    ! Total of h2osoi_liq at wilting point down to the depth of irrigation in each column
    ! [kg/m2]
    real(r8) :: h2osoi_liq_wilting_point_tot(bounds%begc:bounds%endc)

    ! h2osoi_liq at the threshold for irrigation in this column [kg/m2]
    real(r8) :: h2osoi_liq_at_threshold

    ! difference between desired soil moisture level for each column and current soil
    ! moisture level [kg/m2] [i.e., mm]
    real(r8) :: deficit(bounds%begc:bounds%endc)

    ! deficit limited by river volume [kg/m2] [i.e., mm]
    real(r8) :: deficit_volr_limited(bounds%begc:bounds%endc)

    ! where do we need to check soil moisture to see if we need to irrigate?
    logical  :: check_for_irrig_patch(bounds%begp:bounds%endp)
    logical  :: check_for_irrig_col(bounds%begc:bounds%endc)

    ! set to true once we have reached the max allowable depth for irrigation in a given
    ! column
    logical  :: reached_max_depth(bounds%begc:bounds%endc)

    ! Whether we should limit deficits by available volr
    logical :: limit_irrigation

    character(len=*), parameter :: subname = 'CalcIrrigationNeeded'
    !-----------------------------------------------------------------------
    
    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(elai) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, nlevgrnd/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(eff_porosity) == (/bounds%endc, nlevgrnd/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(h2osoi_liq) == (/bounds%endc, nlevgrnd/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(volr) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))


    ! Determine if irrigation is needed (over irrigated soil columns)
    
    ! First, determine in what grid cells we need to bother 'measuring' soil water, to see
    ! if we need irrigation
    check_for_irrig_col(bounds%begc:bounds%endc) = .false.
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       g = patch%gridcell(p)

       check_for_irrig_patch(p) = this%PointNeedsCheckForIrrig( &
            pft_type=patch%itype(p), elai=elai(p), &
            londeg=grc%londeg(g))
       if (check_for_irrig_patch(p)) then
          c = patch%column(p)
          check_for_irrig_col(c) = .true.
       end if
    end do

    check_for_irrig_col_filter = col_filter_from_logical_array(bounds, &
         check_for_irrig_col(bounds%begc:bounds%endc))

    ! Initialize some variables
    do fc = 1, check_for_irrig_col_filter%num
       c = check_for_irrig_col_filter%indices(fc)

       reached_max_depth(c) = .false.
       h2osoi_liq_tot(c) = 0._r8
       h2osoi_liq_target_tot(c) = 0._r8
       h2osoi_liq_wilting_point_tot(c) = 0._r8
    end do

    ! Now 'measure' soil water for the grid cells identified above and see if the soil is
    ! dry enough to warrant irrigation
    do j = 1,nlevsoi
       do fc = 1, check_for_irrig_col_filter%num
          c = check_for_irrig_col_filter%indices(fc)

          if (.not. reached_max_depth(c)) then
             if (col%z(c,j) > this%params%irrig_depth) then
                reached_max_depth(c) = .true.
             else if (j > col%nbedrock(c)) then
                reached_max_depth(c) = .true.
             else if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                ! if level L was frozen, then we don't look at any levels below L
                reached_max_depth(c) = .true.
             else
                h2osoi_liq_tot(c) = h2osoi_liq_tot(c) + h2osoi_liq(c,j)

                h2osoi_liq_target = this%RelsatToH2osoi( &
                     relsat = this%relsat_target_col(c,j), &
                     eff_porosity = eff_porosity(c,j), &
                     dz = col%dz(c,j))
                h2osoi_liq_target_tot(c) = h2osoi_liq_target_tot(c) + &
                     h2osoi_liq_target

                h2osoi_liq_wilting_point = this%RelsatToH2osoi( &
                     relsat = this%relsat_wilting_point_col(c,j), &
                     eff_porosity = eff_porosity(c,j), &
                     dz = col%dz(c,j))
                h2osoi_liq_wilting_point_tot(c) = h2osoi_liq_wilting_point_tot(c) + &
                     h2osoi_liq_wilting_point
             end if
          end if     ! if (.not. reached_max_depth(c))
       end do        ! do fc
    end do           ! do j

    ! Compute deficits
    ! First initialize deficits to 0 everywhere; this is needed for later averaging up to gridcell
    deficit(bounds%begc:bounds%endc) = 0._r8
    do fc = 1, check_for_irrig_col_filter%num
       c = check_for_irrig_col_filter%indices(fc)

       h2osoi_liq_at_threshold = h2osoi_liq_wilting_point_tot(c) + &
            this%params%irrig_threshold_fraction * &
            (h2osoi_liq_target_tot(c) - h2osoi_liq_wilting_point_tot(c))
       if (h2osoi_liq_tot(c) < h2osoi_liq_at_threshold) then
          deficit(c) = h2osoi_liq_target_tot(c) - h2osoi_liq_tot(c)
          ! deficit shouldn't be less than 0: if it is, that implies that the
          ! irrigation target is less than the irrigation threshold, which is not
          ! supposed to happen
          if (deficit(c) < 0._r8) then
             write(iulog,*) subname//' ERROR: deficit < 0'
             write(iulog,*) 'This implies that irrigation target is less than irrigatio&
                  &n threshold, which should never happen'
             call endrun(decomp_index=c, clmlevel=namec, msg='deficit < 0 '// &
                  errMsg(sourcefile, __LINE__))
          end if
       else
          ! We're above the threshold - so don't irrigate
          deficit(c) = 0._r8
       end if
    end do

    ! Limit deficits by available volr, if desired. Note that we cannot do this limiting
    ! if running without a prognostic river model, since we need river volume for this
    ! limiting.
    !
    ! NOTE(wjs, 2016-11-22) In principle we could base this on rof_present rather than
    ! rof_prognostic, but that would depend on the data runoff (drof) model sending river
    ! volume, which it currently does not.
    limit_irrigation = (this%params%limit_irrigation_if_rof_enabled .and. rof_prognostic)
    if (limit_irrigation) then
       call this%CalcDeficitVolrLimited( &
            bounds = bounds, &
            check_for_irrig_col_filter = check_for_irrig_col_filter, &
            deficit = deficit(bounds%begc:bounds%endc), &
            volr = volr(bounds%begg:bounds%endg), &
            deficit_volr_limited = deficit_volr_limited(bounds%begc:bounds%endc))
    else
       deficit_volr_limited(bounds%begc:bounds%endc) = deficit(bounds%begc:bounds%endc)
    end if

    ! Convert deficits to irrigation rate
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       c = patch%column(p)

       if (check_for_irrig_patch(p)) then

          ! Convert units from mm to mm/sec
          this%sfc_irrig_rate_patch(p) = deficit_volr_limited(c) / &
               (this%dtime*this%irrig_nsteps_per_day)
          this%irrig_rate_demand_patch(p) = deficit(c) / &
               (this%dtime*this%irrig_nsteps_per_day)

          ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
          ! in this case, we'll irrigate by 0 for the given number of time steps
          this%n_irrig_steps_left_patch(p) = this%irrig_nsteps_per_day
       end if
    end do

  end subroutine CalcIrrigationNeeded

  !-----------------------------------------------------------------------
  function PointNeedsCheckForIrrig(this, pft_type, elai, londeg) &
       result(check_for_irrig)
    !
    ! !DESCRIPTION:
    ! Determine whether a given patch needs to be checked for irrigation now.
    !
    ! !USES:
    use clm_time_manager, only : get_local_time
    use pftconMod       , only : pftcon
    !
    ! !ARGUMENTS:
    logical :: check_for_irrig  ! function result
    class(irrigation_type), intent(in) :: this
    integer , intent(in) :: pft_type  ! type of pft in this patch
    real(r8), intent(in) :: elai      ! one-sided leaf area index with burying by snow
    real(r8), intent(in) :: londeg    ! longitude (degrees)
    !
    ! !LOCAL VARIABLES:
    ! number of seconds since the prescribed irrigation start time
    integer  :: seconds_since_irrig_start_time

    character(len=*), parameter :: subname = 'PointNeedsCheckForIrrig'
    !-----------------------------------------------------------------------
    
    if (pftcon%irrigated(pft_type) == 1._r8 .and. &
         elai > this%params%irrig_min_lai) then
       ! see if it's the right time of day to start irrigating:
       seconds_since_irrig_start_time = get_local_time( londeg, starttime=this%params%irrig_start_time, offset=-this%dtime )
       if (seconds_since_irrig_start_time < this%dtime) then
          check_for_irrig         = .true.
       else
          check_for_irrig    = .false.
       end if
    else
       check_for_irrig       = .false.
    end if

  end function PointNeedsCheckForIrrig

  !-----------------------------------------------------------------------
  subroutine CalcDeficitVolrLimited(this, bounds, check_for_irrig_col_filter, &
       deficit, volr, deficit_volr_limited)
    !
    ! !DESCRIPTION:
    ! Calculates deficit limited by river volume for each column.
    !
    ! The output array (deficit_volr_limited) is set to 0 outside the given filter.
    !
    ! If the deficit is lower than the volr-based threshold for a given gridcell, then the
    ! volr-limited deficit is simply the original deficit; if the deficit is higher than
    ! the volr-based threshold for a gridcell, then the volr-limited deficit for each
    ! column in the gridcell is equal to:
    !    (original deficit in the column) * (volr-based threshold)/(gridcell-level deficit)
    !
    ! The logic here relies on the fact that all irrigated columns in a given grid cell
    ! will be irrigated at the same time of day. (As opposed to, say, some being irrigated
    ! at 4 AM and others being irrigated at 6 AM: in that case we couldn't just average
    ! this time step's column-level deficits to the gridcell to get the total daily
    ! irrigation deficit: we'd need to account for the fact that we have some irrigation
    ! commitment from earlier deficits, etc.)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(in) :: this
    type(bounds_type)      , intent(in) :: bounds

    ! filter of columns where we need to check for irrigation
    type(filter_col_type)  , intent(in) :: check_for_irrig_col_filter

    ! original deficit: the amount by which we would irrigate if there were no volr
    ! limitation [kg/m2] [i.e., mm]
    !
    ! Should be set for all columns, even those outside the filter, so that averaging to
    ! grid cell happens properly
    real(r8), intent(in) :: deficit( bounds%begc: )

    ! river water volume [m3]
    real(r8), intent(in) :: volr( bounds%begg: )

    ! deficit limited by river volume [kg/m2] [i.e., mm]
    real(r8), intent(out) :: deficit_volr_limited( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer :: fc ! column filter index
    integer :: c  ! column index
    integer :: g  ! gridcell index

    ! Deficit averaged up to gridcell level [kg/m2] [i.e., mm]
    real(r8) :: deficit_grc(bounds%begg:bounds%endg)

    real(r8) :: available_volr ! volr available for withdrawal [m3]

    real(r8) :: max_deficit_supported_by_volr ! [kg/m2] [i.e., mm]

    ! ratio of deficit_volr_limited to deficit for each grid cell
    real(r8) :: deficit_limited_ratio_grc(bounds%begg:bounds%endg)

    character(len=*), parameter :: subname = 'CalcDeficitVolrLimited'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(deficit) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(volr) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(deficit_volr_limited) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    call c2g(bounds, &
         carr = deficit(bounds%begc:bounds%endc), &
         garr = deficit_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    do g = bounds%begg, bounds%endg
       if (volr(g) > 0._r8) then
          available_volr = volr(g) * (1._r8 - this%params%irrig_river_volume_threshold)
          max_deficit_supported_by_volr = available_volr / grc%area(g) * m3_over_km2_to_mm
       else
          ! Ensure that negative volr is treated the same as 0 volr
          max_deficit_supported_by_volr = 0._r8
       end if

       if (deficit_grc(g) > max_deficit_supported_by_volr) then
          ! inadequate river storage, adjust irrigation demand
          deficit_limited_ratio_grc(g) = max_deficit_supported_by_volr / deficit_grc(g)
       else
          ! adequate river storage, no adjustment to irrigation demand
          deficit_limited_ratio_grc(g) = 1._r8
       end if
    end do

    deficit_volr_limited(bounds%begc:bounds%endc) = 0._r8
    do fc = 1, check_for_irrig_col_filter%num
       c = check_for_irrig_col_filter%indices(fc)

       g = col%gridcell(c)
       deficit_volr_limited(c) = deficit(c) * deficit_limited_ratio_grc(g)
    end do

  end subroutine CalcDeficitVolrLimited


  !-----------------------------------------------------------------------
  pure function RelsatToH2osoi(this, relsat, eff_porosity, dz) result(h2osoi_liq)
    !
    ! !DESCRIPTION:
    ! For a given column and layer, convert relative saturation to kg/m2 water
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: h2osoi_liq   ! function result (kg/m2)
    class(irrigation_type), intent(in) :: this
    real(r8), intent(in) :: relsat       ! relative saturation (0 to 1)
    real(r8), intent(in) :: eff_porosity ! effective porosity (0 to 1)
    real(r8), intent(in) :: dz           ! level thickness (m)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: vol_liq  ! partial volume of liquid water in layer [0 to 1]

    character(len=*), parameter :: subname = 'RelsatToH2osoi'
    !-----------------------------------------------------------------------

    vol_liq    = eff_porosity * relsat
    h2osoi_liq = vol_liq * denh2o * dz

  end function RelsatToH2osoi

  !-----------------------------------------------------------------------
  function UseGroundwaterIrrigation(this)
    !
    ! !DESCRIPTION:
    ! Returns true if we're using groundwater irrigation in this run
    !
    ! !ARGUMENTS:
    implicit none
    class(irrigation_type), intent(in) :: this

    logical :: UseGroundwaterIrrigation  ! function result
    !-----------------------------------------------------------------------

    UseGroundwaterIrrigation = this%params%use_groundwater_irrigation
    
  end function UseGroundwaterIrrigation
  
end module IrrigationMod
