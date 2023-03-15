module IrrigationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates irrigation flux.
  !
  ! Usage:
  !
  !   - Call CalcIrrigationNeeded in order to compute whether and how much irrigation is
  !     needed for the next call to ApplyIrrigation. This should be called once per
  !     timestep.
  ! 
  !   - Call ApplyIrrigation in order to calculate qflx_irrig. This should be called
  !     exactly once per time step, before the first time qflx_irrig is needed by other
  !     parts of the code. It is acceptable for this to be called earlier in the timestep
  !     than CalcIrrigationNeeded.
  !
  !   - Access the timestep's irrigation flux via qflx_irrig_patch or
  !     qflx_irrig_col. These should be treated as read-only.
  !
  ! Design notes:
  !
  !   In principle, ApplyIrrigation and CalcIrrigationNeeded could be combined into a
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
  !   combining ApplyIrrigation and CalcIrrigationNeeded - or at least calling these two
  !   routines from the same place.  In particular: this separation of the irrigation
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
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type, get_proc_global
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog
  use clm_varcon       , only : isecspday, denh2o, spval, namec
  use clm_varpar       , only : nlevsoi, nlevgrnd
  use clm_time_manager , only : get_step_size
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
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
     ! since we won't begin irrigating until the next call to ApplyIrrigation
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

  end type irrigation_params_type


  type, public :: irrigation_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules
     real(r8), pointer, public :: qflx_irrig_patch(:) ! patch irrigation flux (mm H2O/s)
     real(r8), pointer, public :: qflx_irrig_col  (:) ! col irrigation flux (mm H2O/s)

     ! Private data members; set in initialization:
     type(irrigation_params_type) :: params
     integer :: dtime                ! land model time step (sec)
     integer :: irrig_nsteps_per_day ! number of time steps per day in which we irrigate
     real(r8), pointer :: relsat_wilting_point_col(:,:) ! relative saturation at which smp = wilting point [col, nlevsoi]
     real(r8), pointer :: relsat_target_col(:,:)        ! relative saturation at which smp is at the irrigation target [col, nlevsoi]

     ! Private data members; time-varying:
     real(r8), pointer :: irrig_rate_patch            (:) ! current irrigation rate [mm/s]
     real(r8), pointer :: irrig_rate_demand_patch     (:) ! current irrigation rate, neglecting surface water source limitation [mm/s]
     integer , pointer :: n_irrig_steps_left_patch    (:) ! number of time steps for which we still need to irrigate today (if 0, ignore)
     real(r8), pointer :: qflx_irrig_demand_patch     (:) ! irrigation flux neglecting surface water source limitation [mm/s]

   contains
     ! Public routines
     ! COMPILER_BUG(wjs, 2014-10-15, pgi 14.7) Add an "Irrigation" prefix to some  generic routines like "Init"
     ! (without this workaround, pgi compilation fails in restFileMod)
     procedure, public :: Init => IrrigationInit
     procedure, public :: Restart
     procedure, public :: ApplyIrrigation
     procedure, public :: CalcIrrigationNeeded
     procedure, public :: Clean => IrrigationClean ! deallocate memory

     ! Public simply to support unit testing; should not be used from CLM code
     procedure, public :: InitForTesting ! version of Init meant for unit testing
     procedure, public :: RelsatToH2osoi ! convert from relative saturation to kg/m2 water for a single column and layer

     ! Private routines
     procedure, private :: ReadNamelist
     procedure, private :: CheckNamelistValidity   ! Check for validity of input parameters
     procedure, private :: InitAllocate => IrrigationInitAllocate
     procedure, private :: InitHistory => IrrigationInitHistory
     procedure, private :: InitCold => IrrigationInitCold
     procedure, private :: CalcIrrigNstepsPerDay   ! given dtime, calculate irrig_nsteps_per_day
     procedure, private :: PointNeedsCheckForIrrig ! whether a given point needs to be checked for irrigation now
     procedure, private :: CalcDeficitVolrLimited  ! calculate deficit limited by river volume for each patch
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
       limit_irrigation_if_rof_enabled) &
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

  end function irrigation_params_constructor


  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================
  
  !------------------------------------------------------------------------
  subroutine IrrigationInit(this, bounds, NLFilename, &
       soilstate_inst, soil_water_retention_curve)
    use SoilStateType , only : soilstate_type

    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    character(len=*)       , intent(in)    :: NLFilename ! Namelist filename
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

    call this%ReadNamelist(NLFilename)
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
    this%irrig_nsteps_per_day = this%CalcIrrigNstepsPerDay(dtime)
    this%relsat_wilting_point_col(:,:) = relsat_wilting_point(:,:)
    this%relsat_target_col(:,:) = relsat_target(:,:)

  end subroutine InitForTesting

  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename)
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
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    class(irrigation_type) , intent(inout) :: this
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

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: nmlname = 'irrigation_inparm'

    character(len=*), parameter :: subname = 'ReadNamelist'
    !-----------------------------------------------------------------------

    namelist /irrigation_inparm/ irrig_min_lai, irrig_start_time, irrig_length, &
         irrig_target_smp, irrig_depth, irrig_threshold_fraction, &
         irrig_river_volume_threshold, limit_irrigation_if_rof_enabled

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

    this%params = irrigation_params_type( &
         irrig_min_lai = irrig_min_lai, &
         irrig_start_time = irrig_start_time, &
         irrig_length = irrig_length, &
         irrig_target_smp = irrig_target_smp, &
         irrig_depth = irrig_depth, &
         irrig_threshold_fraction = irrig_threshold_fraction, &
         irrig_river_volume_threshold = irrig_river_volume_threshold, &
         limit_irrigation_if_rof_enabled = limit_irrigation_if_rof_enabled)

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
       write(iulog,*) ' '

       call this%CheckNamelistValidity()
    end if

  end subroutine ReadNamelist

  !-----------------------------------------------------------------------
  subroutine CheckNamelistValidity(this)
    !
    ! !DESCRIPTION:
    ! Check for validity of input parameters.
    !
    ! Assumes that the inputs have already been packed into 'this%params'.
    !
    ! Only needs to be called by the master task, since parameters are the same for all
    ! tasks.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(irrigation_type), intent(in) :: this
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

    allocate(this%qflx_irrig_patch         (begp:endp))          ; this%qflx_irrig_patch         (:)   = nan
    allocate(this%qflx_irrig_demand_patch  (begp:endp))          ; this%qflx_irrig_demand_patch  (:)   = nan
    allocate(this%qflx_irrig_col           (begc:endc))          ; this%qflx_irrig_col           (:)   = nan
    allocate(this%relsat_wilting_point_col (begc:endc,nlevsoi)) ; this%relsat_wilting_point_col (:,:) = nan
    allocate(this%relsat_target_col        (begc:endc,nlevsoi)) ; this%relsat_target_col        (:,:) = nan
    allocate(this%irrig_rate_patch         (begp:endp))          ; this%irrig_rate_patch         (:)   = nan
    allocate(this%irrig_rate_demand_patch  (begp:endp))          ; this%irrig_rate_demand_patch  (:)   = nan
    allocate(this%n_irrig_steps_left_patch (begp:endp))          ; this%n_irrig_steps_left_patch (:)   = 0

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

    this%qflx_irrig_patch(begp:endp) = spval
    call hist_addfld1d (fname='QIRRIG', units='mm/s', &
         avgflag='A', long_name='water added through irrigation', &
         ptr_patch=this%qflx_irrig_patch)

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

    this%dtime = get_step_size()
    this%irrig_nsteps_per_day = this%CalcIrrigNstepsPerDay(this%dtime)

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
            long_name='irrigation rate', units='mm/s', &
            interpinic_flag='interp', readvar=readvar, data=this%irrig_rate_patch)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%irrig_rate_patch = 0.0_r8
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
    
    deallocate(this%qflx_irrig_patch)
    deallocate(this%qflx_irrig_demand_patch)
    deallocate(this%qflx_irrig_col)
    deallocate(this%relsat_wilting_point_col)
    deallocate(this%relsat_target_col)
    deallocate(this%irrig_rate_patch)
    deallocate(this%irrig_rate_demand_patch)
    deallocate(this%n_irrig_steps_left_patch)

  end subroutine IrrigationClean


  ! ========================================================================
  ! Science routines
  ! ========================================================================
  
  !-----------------------------------------------------------------------
  subroutine ApplyIrrigation(this, bounds)
    !
    ! !DESCRIPTION:
    ! Apply the irrigation computed by CalcIrrigationNeeded to qflx_irrig.
    !
    ! Should be called once, AND ONLY ONCE, per time step. After this is called, you may
    ! access qflx_irrig_patch or qflx_irrig_col.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p  ! patch index
    integer :: g  ! grid cell index
    
    character(len=*), parameter :: subname = 'ApplyIrrigation'

    !-----------------------------------------------------------------------

    ! This should be called exactly once per time step, so that this counter decrease
    ! works correctly.

    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)

       if (this%n_irrig_steps_left_patch(p) > 0) then
          this%qflx_irrig_patch(p)         = this%irrig_rate_patch(p)
          this%qflx_irrig_demand_patch(p)  = this%irrig_rate_demand_patch(p)
          this%n_irrig_steps_left_patch(p) = this%n_irrig_steps_left_patch(p) - 1
       else
          this%qflx_irrig_patch(p)        = 0._r8
          this%qflx_irrig_demand_patch(p) = 0._r8
       end if

    end do

    call p2c (bounds = bounds, &
         parr = this%qflx_irrig_patch(bounds%begp:bounds%endp), &
         carr = this%qflx_irrig_col(bounds%begc:bounds%endc), &
         p2c_scale_type = 'unity')

  end subroutine ApplyIrrigation


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
          this%irrig_rate_patch(p) = deficit_volr_limited(c) / &
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

end module IrrigationMod
