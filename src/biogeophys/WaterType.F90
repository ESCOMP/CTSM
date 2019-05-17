module WaterType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Container for derived types relating to water, both for bulk water and for isotopes
  ! and other tracers.
  !
  ! Variables pertaining to bulk water can be accessed in two ways:
  !
  !    (1) Using water_inst%water*bulk_inst
  !
  !    (2) As one of the indices in water_inst%bulk_and_tracers(:)%water*_inst
  !
  !    Method (1) is greatly preferable when you are just operating on bulk water. Method
  !    (2) is just meant to be used when you are doing the same operation on bulk water
  !    and all water tracers.
  !
  ! To loop through bulk and all tracers, use code like this:
  !    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
  !       associate( &
  !            waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst, &
  !            [and other associations, as necessary])
  !       [Do calculations involving waterflux_inst, etc.]
  !       end associate
  !    end do
  !
  ! To loop through all tracers (not bulk), use code like this:
  !    do i = water_inst%tracers_beg, water_inst%tracers_end
  !       associate( &
  !            waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst, &
  !            [and other associations, as necessary])
  !       [Do calculations involving waterflux_inst, etc.]
  !       end associate
  !    end do
  !
  ! To loop through all isotopes (not bulk or other water tracers), use code like this:
  !    type(water_info_isotope_type), pointer :: iso_info
  !
  !    do i = water_inst%tracers_beg, water_inst%tracers_end
  !       if (water_inst%IsIsotope(i)) then
  !          call water_inst%GetIsotopeInfo(i, iso_info)
  !          associate( &
  !               waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst, &
  !               [and other associations, as necessary])
  !          [Do calculations involving iso_info, waterflux_inst, etc.]
  !          end associate
  !       end if
  !    end do
  !
  ! The associate statements given above aren't crucial. If the block of code refers to
  ! multiple instances (waterstate, waterflux, etc.), but only refers to each one once or
  ! twice, it can be best to just have:
  !    associate(bulk_or_tracer => water_inst%bulk_and_tracers(i))
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod             , only : r8 => shr_kind_r8
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use abortutils               , only : endrun
  use decompMod                , only : bounds_type
  use clm_varctl               , only : iulog
  use clm_varpar               , only : nlevsno
  use ncdio_pio                , only : file_desc_t
  use WaterFluxBulkType        , only : waterfluxbulk_type
  use WaterFluxType            , only : waterflux_type
  use WaterStateBulkType       , only : waterstatebulk_type
  use WaterStateType           , only : waterstate_type
  use WaterDiagnosticType      , only : waterdiagnostic_type
  use WaterDiagnosticBulkType  , only : waterdiagnosticbulk_type
  use WaterBalanceType         , only : waterbalance_type
  use WaterInfoBaseType        , only : water_info_base_type
  use WaterInfoBulkType        , only : water_info_bulk_type
  use WaterInfoTracerType      , only : water_info_tracer_type
  use WaterInfoIsotopeType     , only : water_info_isotope_type
  use Waterlnd2atmType         , only : waterlnd2atm_type
  use Waterlnd2atmBulkType     , only : waterlnd2atmbulk_type
  use Wateratm2lndType         , only : wateratm2lnd_type
  use Wateratm2lndBulkType     , only : wateratm2lndbulk_type
  use WaterTracerContainerType , only : water_tracer_container_type
  use WaterTracerUtils         , only : CompareBulkToTracer

  implicit none
  private

  !
  ! !PRIVATE TYPES:

  ! This type holds instances needed for bulk water or for a single tracer
  type, private :: bulk_or_tracer_type
     private

     ! ------------------------------------------------------------------------
     ! Public data members
     ! ------------------------------------------------------------------------

     class(waterflux_type)       , pointer, public :: waterflux_inst
     class(waterstate_type)      , pointer, public :: waterstate_inst
     class(waterdiagnostic_type) , pointer, public :: waterdiagnostic_inst
     class(waterbalance_type)    , pointer, public :: waterbalance_inst
     class(waterlnd2atm_type)    , pointer, public :: waterlnd2atm_inst
     class(wateratm2lnd_type)    , pointer, public :: wateratm2lnd_inst

     ! ------------------------------------------------------------------------
     ! Private data members
     ! ------------------------------------------------------------------------

     logical :: is_isotope = .false.
     class(water_info_base_type) , pointer :: info
     type(water_tracer_container_type) :: vars

  end type bulk_or_tracer_type

  !
  ! !PUBLIC TYPES:

  ! water_params_type is public for the sake of unit tests
  type, public :: water_params_type
     private

     ! Whether we add tracers that are used for the tracer consistency checks
     logical :: enable_consistency_checks

     ! Whether we add tracers that are used for isotopes
     logical :: enable_isotopes
  end type water_params_type

  type, public :: water_type
     private

     ! ------------------------------------------------------------------------
     ! Public data members
     ! ------------------------------------------------------------------------

     ! indices into the bulk_and_tracers array
     integer, public :: bulk_and_tracers_beg  ! first index when looping over bulk & tracers
     integer, public :: bulk_and_tracers_end  ! last index when looping over bulk & tracers
     integer, public :: tracers_beg           ! first index when looping over just tracers
     integer, public :: tracers_end           ! last index when looping over just tracers
     integer, public :: i_bulk                ! index of bulk in arrays that contain both bulk and tracers

     type(waterfluxbulk_type), pointer, public       :: waterfluxbulk_inst
     type(waterstatebulk_type), pointer, public      :: waterstatebulk_inst
     type(waterdiagnosticbulk_type), pointer, public :: waterdiagnosticbulk_inst
     type(waterbalance_type), pointer, public        :: waterbalancebulk_inst
     type(waterlnd2atmbulk_type), pointer, public    :: waterlnd2atmbulk_inst
     type(wateratm2lndbulk_type), pointer, public    :: wateratm2lndbulk_inst

     type(bulk_or_tracer_type), allocatable, public  :: bulk_and_tracers(:)

     ! ------------------------------------------------------------------------
     ! Private data members
     ! ------------------------------------------------------------------------

     type(water_params_type) :: params
     integer :: bulk_tracer_index  ! index of the tracer that replicates bulk water (-1 if it doesn't exist)

   contains
     ! Public routines
     procedure, public :: Init
     procedure, public :: InitForTesting  ! Init routine just for unit tests
     procedure, public :: InitAccBuffer
     procedure, public :: InitAccVars
     procedure, public :: UpdateAccVars
     procedure, public :: Restart
     procedure, public :: IsIsotope       ! Return true if a given tracer is an isotope
     procedure, public :: GetIsotopeInfo  ! Get a pointer to the object storing isotope info for a given tracer
     procedure, public :: GetBulkTracerIndex ! Get the index of the tracer that replicates bulk water
     procedure, public :: DoConsistencyCheck ! Whether TracerConsistencyCheck should be called in this run
     procedure, public :: TracerConsistencyCheck
     procedure, public :: Summary            ! Calculate end-of-timestep summaries of water diagnostic terms

     ! Private routines
     procedure, private :: DoInit
     procedure, private :: ReadNamelist
     procedure, private :: SetupTracerInfo
     procedure, private :: AllocateBulk
     procedure, private :: AllocateTracer
  end type water_type

  interface water_params_type
     module procedure water_params_constructor
  end interface water_params_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  function water_params_constructor(enable_consistency_checks, enable_isotopes) &
       result(params)
    !
    ! !DESCRIPTION:
    ! Creates a new instance of water_params_type
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(water_params_type) :: params  ! function result
    logical, intent(in) :: enable_consistency_checks
    logical, intent(in) :: enable_isotopes
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'water_params_constructor'
    !-----------------------------------------------------------------------

    params%enable_consistency_checks = enable_consistency_checks
    params%enable_isotopes = enable_isotopes
  end function water_params_constructor

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename, &
       h2osno_col, snow_depth_col, watsat_col, t_soisno_col, use_aquifer_layer)
    !
    ! !DESCRIPTION:
    ! Initialize all water variables
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    character(len=*) , intent(in)    :: NLFilename ! Namelist filename
    real(r8)          , intent(in) :: h2osno_col(bounds%begc:)
    real(r8)          , intent(in) :: snow_depth_col(bounds%begc:)
    real(r8)          , intent(in) :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)          , intent(in) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    logical           , intent(in) :: use_aquifer_layer ! whether an aquifer layer is used in this run
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    call this%ReadNamelist(NLFilename)
    call this%DoInit(bounds = bounds, &
         h2osno_col = h2osno_col, &
         snow_depth_col = snow_depth_col, &
         watsat_col = watsat_col, &
         t_soisno_col = t_soisno_col, &
         use_aquifer_layer = use_aquifer_layer)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitForTesting(this, bounds, params, &
       h2osno_col, snow_depth_col, watsat_col, t_soisno_col, use_aquifer_layer)
    !
    ! !DESCRIPTION:
    ! Version of Init routine just for unit tests
    !
    ! This version has params passed in directly instead of reading from namelist
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    type(water_params_type), intent(in) :: params
    real(r8)          , intent(in) :: h2osno_col(bounds%begc:)
    real(r8)          , intent(in) :: snow_depth_col(bounds%begc:)
    real(r8)          , intent(in) :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)          , intent(in) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    logical , intent(in), optional :: use_aquifer_layer ! whether an aquifer layer is used in this run (false by default)
    !
    ! !LOCAL VARIABLES:
    logical :: l_use_aquifer_layer

    character(len=*), parameter :: subname = 'InitForTesting'
    !-----------------------------------------------------------------------

    l_use_aquifer_layer = .false.
    if (present(use_aquifer_layer)) then
       l_use_aquifer_layer = use_aquifer_layer
    end if

    this%params = params
    call this%DoInit(bounds = bounds, &
         h2osno_col = h2osno_col, &
         snow_depth_col = snow_depth_col, &
         watsat_col = watsat_col, &
         t_soisno_col = t_soisno_col, &
         use_aquifer_layer = l_use_aquifer_layer)

  end subroutine InitForTesting

  !-----------------------------------------------------------------------
  subroutine DoInit(this, bounds, &
       h2osno_col, snow_depth_col, watsat_col, t_soisno_col, use_aquifer_layer)
    !
    ! !DESCRIPTION:
    ! Actually do the initialization (shared between main Init routine and InitForTesting)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    real(r8)         , intent(in) :: h2osno_col(bounds%begc:)
    real(r8)         , intent(in) :: snow_depth_col(bounds%begc:)
    real(r8)         , intent(in) :: watsat_col(bounds%begc:, 1:)            ! volumetric soil water at saturation (porosity)
    real(r8)         , intent(in) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    logical          , intent(in) :: use_aquifer_layer ! whether an aquifer layer is used in this run
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: i

    character(len=*), parameter :: subname = 'DoInit'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    SHR_ASSERT_ALL((ubound(h2osno_col) == [endc]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(snow_depth_col) == [endc]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(watsat_col, 1) == endc), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno_col, 1) == endc), errMsg(sourcefile, __LINE__))

    call this%SetupTracerInfo()

    call this%AllocateBulk()

    associate( &
         bulk_info => this%bulk_and_tracers(this%i_bulk)%info, &
         bulk_vars => this%bulk_and_tracers(this%i_bulk)%vars &
         )

    call bulk_vars%init()

    call this%waterstatebulk_inst%InitBulk(bounds, &
         bulk_info, &
         bulk_vars, &
         h2osno_input_col = h2osno_col(begc:endc),       &
         watsat_col = watsat_col(begc:endc, 1:),   &
         t_soisno_col = t_soisno_col(begc:endc, -nlevsno+1:), &
         use_aquifer_layer = use_aquifer_layer)

    call this%waterdiagnosticbulk_inst%InitBulk(bounds, &
         bulk_info, &
         bulk_vars, &
         snow_depth_input_col = snow_depth_col(begc:endc),    &
         waterstatebulk_inst = this%waterstatebulk_inst )

    call this%waterbalancebulk_inst%Init(bounds, &
         bulk_info, &
         bulk_vars)

    call this%waterfluxbulk_inst%InitBulk(bounds, &
         bulk_info, &
         bulk_vars)

    call this%waterlnd2atmbulk_inst%InitBulk(bounds, &
         bulk_info, &
         bulk_vars)

    call this%wateratm2lndbulk_inst%InitBulk(bounds, &
         bulk_info, &
         bulk_vars)

    call bulk_vars%complete_setup()

    end associate

    do i = this%tracers_beg, this%tracers_end

       call this%AllocateTracer(i)

       call this%bulk_and_tracers(i)%vars%init()

       call this%bulk_and_tracers(i)%waterstate_inst%Init(bounds, &
            this%bulk_and_tracers(i)%info, &
            this%bulk_and_tracers(i)%vars, &
            h2osno_input_col = h2osno_col(begc:endc),       &
            watsat_col = watsat_col(begc:endc, 1:),   &
            t_soisno_col = t_soisno_col(begc:endc, -nlevsno+1:), &
            use_aquifer_layer = use_aquifer_layer)

       call this%bulk_and_tracers(i)%waterdiagnostic_inst%Init(bounds, &
            this%bulk_and_tracers(i)%info, &
            this%bulk_and_tracers(i)%vars)

       call this%bulk_and_tracers(i)%waterbalance_inst%Init(bounds, &
            this%bulk_and_tracers(i)%info, &
            this%bulk_and_tracers(i)%vars)

       call this%bulk_and_tracers(i)%waterflux_inst%Init(bounds, &
            this%bulk_and_tracers(i)%info, &
            this%bulk_and_tracers(i)%vars)

       call this%bulk_and_tracers(i)%waterlnd2atm_inst%Init(bounds, &
            this%bulk_and_tracers(i)%info, &
            this%bulk_and_tracers(i)%vars)

       call this%bulk_and_tracers(i)%wateratm2lnd_inst%Init(bounds, &
            this%bulk_and_tracers(i)%info, &
            this%bulk_and_tracers(i)%vars)

       call this%bulk_and_tracers(i)%vars%complete_setup()

    end do

  end subroutine DoInit


  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the water_tracers namelist; set this%params
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    class(water_type) , intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    ! temporary local variables corresponding to the namelist items
    logical :: enable_water_tracer_consistency_checks
    logical :: enable_water_isotopes

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: nmlname = 'water_tracers'

    character(len=*), parameter :: subname = 'ReadNamelist'
    !-----------------------------------------------------------------------

    namelist /water_tracers_inparm/ &
         enable_water_tracer_consistency_checks, enable_water_isotopes

    ! Initialize namelist variables to defaults
    enable_water_tracer_consistency_checks = .false.
    enable_water_isotopes = .false.

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=water_tracers_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast(enable_water_tracer_consistency_checks, mpicom)
    call shr_mpi_bcast(enable_water_isotopes, mpicom)

    if (masterproc) then
       write(iulog,*)
       write(iulog,*) nmlname, ' settings'
       write(iulog,nml=water_tracers_inparm)
       write(iulog,*)
    end if

    this%params = water_params_type( &
         enable_consistency_checks = enable_water_tracer_consistency_checks, &
         enable_isotopes = enable_water_isotopes)

  end subroutine ReadNamelist


  !-----------------------------------------------------------------------
  subroutine SetupTracerInfo(this)
    !
    ! !DESCRIPTION:
    ! Setup information on each water tracer
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: num_tracers
    integer :: tracer_num
    logical :: enable_bulk_tracer

    character(len=*), parameter :: subname = 'SetupTracerInfo'
    !-----------------------------------------------------------------------

    this%bulk_tracer_index = -1

    num_tracers = 0
    if (this%params%enable_consistency_checks .or. this%params%enable_isotopes) then
       ! NOTE(wjs, 2018-09-05) From looking at the old water isotope code, it looks like
       ! we may need the bulk tracer even if we're not doing consistency checks, in order
       ! to do some roundoff-level adjustments.
       enable_bulk_tracer = .true.
    else
       enable_bulk_tracer = .false.
    end if

    if (enable_bulk_tracer) then
       num_tracers = num_tracers + 1
    end if
    if (this%params%enable_isotopes) then
       num_tracers = num_tracers + 2
    end if
    if (this%params%enable_consistency_checks) then
       num_tracers = num_tracers + 3
    end if

    this%bulk_and_tracers_beg = 0
    this%tracers_beg = 1
    this%bulk_and_tracers_end = num_tracers
    this%tracers_end = num_tracers
    this%i_bulk = 0

    allocate(this%bulk_and_tracers(this%bulk_and_tracers_beg:this%bulk_and_tracers_end))

    allocate(this%bulk_and_tracers(this%i_bulk)%info, source = water_info_bulk_type())

    tracer_num = 1
    if (enable_bulk_tracer) then
       allocate(this%bulk_and_tracers(tracer_num)%info, source = water_info_isotope_type( &
            tracer_name = 'H2OTR', &
            ratio = 1._r8, &
            included_in_consistency_check = .true., &
            communicated_with_coupler = .false.))
       this%bulk_and_tracers(tracer_num)%is_isotope = .true.
       this%bulk_tracer_index = tracer_num
       tracer_num = tracer_num + 1
    end if
    if (this%params%enable_isotopes) then
       allocate(this%bulk_and_tracers(tracer_num)%info, source = water_info_isotope_type( &
            tracer_name = 'HDO', &
            ratio = 0.9_r8, &
            included_in_consistency_check = .false., &
            communicated_with_coupler = .false.))
       this%bulk_and_tracers(tracer_num)%is_isotope = .true.
       tracer_num = tracer_num + 1

       allocate(this%bulk_and_tracers(tracer_num)%info, source = water_info_isotope_type( &
            tracer_name = 'H218O', &
            ratio = 0.5_r8, &
            included_in_consistency_check = .false., &
            communicated_with_coupler = .false.))
       this%bulk_and_tracers(tracer_num)%is_isotope = .true.
       tracer_num = tracer_num + 1
    end if
    if (this%params%enable_consistency_checks) then
       allocate(this%bulk_and_tracers(tracer_num)%info, source = water_info_isotope_type( &
            tracer_name = 'TESTMED', &
            ratio = 0.1_r8, &
            included_in_consistency_check = .true., &
            communicated_with_coupler = .false.))
       this%bulk_and_tracers(tracer_num)%is_isotope = .true.
       tracer_num = tracer_num + 1

       allocate(this%bulk_and_tracers(tracer_num)%info, source = water_info_isotope_type( &
            tracer_name = 'TESTSMALL', &
            ratio = 1.0e-10_r8, &
            included_in_consistency_check = .true., &
            communicated_with_coupler = .false.))
       this%bulk_and_tracers(tracer_num)%is_isotope = .true.
       tracer_num = tracer_num + 1

       allocate(this%bulk_and_tracers(tracer_num)%info, source = water_info_isotope_type( &
            tracer_name = 'TESTBIG', &
            ratio = 10._r8, &
            included_in_consistency_check = .true., &
            communicated_with_coupler = .false.))
       this%bulk_and_tracers(tracer_num)%is_isotope = .true.
       tracer_num = tracer_num + 1
    end if


    if (tracer_num - 1 /= num_tracers) then
       write(iulog,*) subname//' ERROR: tracer_num discrepancy'
       write(iulog,*) 'num_tracers = ', num_tracers
       write(iulog,*) 'but added ', tracer_num - 1, ' tracers'
       call endrun(msg='tracer_num discrepancy '//errMsg(sourcefile, __LINE__))
    end if

  end subroutine SetupTracerInfo

  !-----------------------------------------------------------------------
  subroutine AllocateBulk(this)
    !
    ! !DESCRIPTION:
    ! Allocate each of the bulk objects
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'AllocateBulk'
    !-----------------------------------------------------------------------

    associate( &
         i_bulk => this%i_bulk &
         )

    allocate(this%waterfluxbulk_inst)
    this%bulk_and_tracers(i_bulk)%waterflux_inst => this%waterfluxbulk_inst

    allocate(this%waterstatebulk_inst)
    this%bulk_and_tracers(i_bulk)%waterstate_inst => this%waterstatebulk_inst

    allocate(this%waterdiagnosticbulk_inst)
    this%bulk_and_tracers(i_bulk)%waterdiagnostic_inst => this%waterdiagnosticbulk_inst

    allocate(this%waterbalancebulk_inst)
    this%bulk_and_tracers(i_bulk)%waterbalance_inst => this%waterbalancebulk_inst

    allocate(this%waterlnd2atmbulk_inst)
    this%bulk_and_tracers(i_bulk)%waterlnd2atm_inst => this%waterlnd2atmbulk_inst

    allocate(this%wateratm2lndbulk_inst)
    this%bulk_and_tracers(i_bulk)%wateratm2lnd_inst => this%wateratm2lndbulk_inst

    end associate

  end subroutine AllocateBulk

  !-----------------------------------------------------------------------
  subroutine AllocateTracer(this, i)
    !
    ! !DESCRIPTION:
    ! Allocate each of the tracer objects for tracer i
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    integer, intent(in) :: i  ! tracer number
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'AllocateTracer'
    !-----------------------------------------------------------------------

    allocate(waterflux_type :: this%bulk_and_tracers(i)%waterflux_inst)
    allocate(waterstate_type :: this%bulk_and_tracers(i)%waterstate_inst)
    allocate(waterdiagnostic_type :: this%bulk_and_tracers(i)%waterdiagnostic_inst)
    allocate(waterbalance_type :: this%bulk_and_tracers(i)%waterbalance_inst)
    allocate(waterlnd2atm_type :: this%bulk_and_tracers(i)%waterlnd2atm_inst)
    allocate(wateratm2lnd_type :: this%bulk_and_tracers(i)%wateratm2lnd_inst)

  end subroutine AllocateTracer

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all water variables
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds

    character(len=*), parameter :: subname = 'InitAccBuffer'
    !-----------------------------------------------------------------------

    call this%waterfluxbulk_inst%InitAccBuffer(bounds)
    call this%wateratm2lndbulk_inst%InitAccBuffer(bounds)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize variables that are associated with accumulated fields
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------

    call this%waterfluxbulk_inst%initAccVars(bounds)
    call this%wateratm2lndbulk_inst%initAccVars(bounds)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Update accumulated variables
    !
    ! Should be called every time step
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'UpdateAccVars'
    !-----------------------------------------------------------------------

    call this%waterfluxbulk_inst%UpdateAccVars(bounds)
    call this%wateratm2lndbulk_inst%UpdateAccVars(bounds)

  end subroutine UpdateAccVars


  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       watsat_col)
    !
    ! !DESCRIPTION:
    ! Read/write information to/from restart file for all water variables
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    real(r8)         , intent(in)    :: watsat_col (bounds%begc:, 1:)  ! volumetric soil water at saturation (porosity)
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(watsat_col, 1) == bounds%endc), errMsg(sourcefile, __LINE__))

    call this%waterfluxbulk_inst%restartBulk (bounds, ncid, flag=flag)

    call this%waterstatebulk_inst%restartBulk (bounds, ncid, flag=flag, &
         watsat_col=watsat_col(bounds%begc:bounds%endc,:))

    call this%waterdiagnosticbulk_inst%restartBulk (bounds, ncid, flag=flag)

    do i = this%tracers_beg, this%tracers_end

       call this%bulk_and_tracers(i)%waterflux_inst%Restart(bounds, ncid, flag=flag)

       call this%bulk_and_tracers(i)%waterstate_inst%Restart(bounds, ncid, flag=flag, &
            watsat_col=watsat_col(bounds%begc:bounds%endc,:))

       call this%bulk_and_tracers(i)%waterdiagnostic_inst%Restart(bounds, ncid, flag=flag)

    end do

  end subroutine Restart

  !-----------------------------------------------------------------------
  function IsIsotope(this, i)
    !
    ! !DESCRIPTION:
    ! Returns true if tracer i is an isotope
    !
    ! i must be >= this%tracers_beg and <= this%tracers_end
    !
    ! !ARGUMENTS:
    logical :: IsIsotope  ! function result
    class(water_type), intent(in) :: this
    integer, intent(in) :: i  ! index of tracer
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'IsIsotope'
    !-----------------------------------------------------------------------

    SHR_ASSERT(i >= this%tracers_beg, errMsg(sourcefile, __LINE__))
    SHR_ASSERT(i <= this%tracers_end, errMsg(sourcefile, __LINE__))

    IsIsotope = this%bulk_and_tracers(i)%is_isotope

  end function IsIsotope

  !-----------------------------------------------------------------------
  subroutine GetIsotopeInfo(this, i, isotope_info)
    !
    ! !DESCRIPTION:
    ! Get a pointer to the object storing isotope info for a given tracer
    !
    ! This provides a mechanism for passing the isotope info to subroutines that need it.
    !
    ! i must be >= this%tracers_beg and <= this%tracers_end, and this%IsIsotope(i) must be
    ! true
    !
    ! Assumes that the 'isotope_info' pointer is not currently allocated. (Otherwise this
    ! will result in a memory leak. It is okay for the isotope_info pointer to be
    ! previously associated with something else, though, as long as it doesn't require
    ! deallocation before being associated with something new.)
    !
    ! !ARGUMENTS:
    class(water_type), intent(in) :: this
    integer, intent(in) :: i  ! index of tracer
    type(water_info_isotope_type), pointer, intent(out) :: isotope_info
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'GetIsotopeInfo'
    !-----------------------------------------------------------------------

    SHR_ASSERT(i >= this%tracers_beg, errMsg(sourcefile, __LINE__))
    SHR_ASSERT(i <= this%tracers_end, errMsg(sourcefile, __LINE__))

    select type(info => this%bulk_and_tracers(i)%info)
    type is(water_info_isotope_type)
       isotope_info => info
    class default
       write(iulog,*) subname, ' ERROR: tracer ', i, ' is not an isotope'
       call endrun(subname//' called on a non-isotope tracer')
    end select

  end subroutine GetIsotopeInfo

  !-----------------------------------------------------------------------
  function GetBulkTracerIndex(this) result(index)
    !
    ! !DESCRIPTION:
    ! Get the index of the tracer that replicates bulk water
    !
    ! Returns -1 if there is no tracer that replicates bulk water in this run
    !
    ! !ARGUMENTS:
    integer :: index  ! function result
    class(water_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'GetBulkTracerIndex'
    !-----------------------------------------------------------------------

    index = this%bulk_tracer_index

  end function GetBulkTracerIndex

  !-----------------------------------------------------------------------
  function DoConsistencyCheck(this) result(do_consistency_check)
    !
    ! !DESCRIPTION:
    ! Returns a logical saying whether TracerConsistencyCheck should be called in this run
    !
    ! !ARGUMENTS:
    logical :: do_consistency_check  ! function result
    class(water_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'DoConsistencyCheck'
    !-----------------------------------------------------------------------

    do_consistency_check = this%params%enable_consistency_checks

  end function DoConsistencyCheck


  !------------------------------------------------------------------------
  subroutine TracerConsistencyCheck(this, bounds, caller_location)
    !
    ! !DESCRIPTION:
    ! Check consistency of water tracers with that of bulk water
    !
    ! This should only be called if this%DoConsistencyCheck() returns .true.
    !
    ! !ARGUMENTS:
    class(water_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: caller_location  ! brief description of where this is called from (for error messages)
    !
    ! !LOCAL VARIABLES:
    integer :: i
    integer :: num_vars
    integer :: var_num
    character(len=:), allocatable :: name
    integer :: begi, endi
    real(r8), pointer :: bulk(:)
    real(r8), pointer :: tracer(:)
    character(len=*), parameter :: subname = 'TracerConsistencyCheck'
    !-----------------------------------------------------------------------

    do i = this%tracers_beg, this%tracers_end

       associate( &
            tracer_vars => this%bulk_and_tracers(i)%vars, &
            tracer_info => this%bulk_and_tracers(i)%info, &
            bulk_vars => this%bulk_and_tracers(this%i_bulk)%vars &
            )


       if (tracer_info%is_included_in_consistency_check()) then
          num_vars = tracer_vars%get_num_vars()
          SHR_ASSERT(num_vars == bulk_vars%get_num_vars(), errMsg(sourcefile, __LINE__))

          do var_num = 1, num_vars
             name = tracer_vars%get_description(var_num)
             SHR_ASSERT(name == bulk_vars%get_description(var_num), errMsg(sourcefile, __LINE__))

             call tracer_vars%get_bounds(var_num, bounds, begi, endi)

             call bulk_vars%get_data(var_num, bulk)
             call tracer_vars%get_data(var_num, tracer)

             call CompareBulkToTracer(begi, endi, &
                  bulk   = bulk(begi:endi), &
                  tracer = tracer(begi:endi), &
                  ratio = tracer_info%get_ratio(), &
                  caller_location = caller_location, &
                  name = name)

          end do
       end if
       end associate

     end do

  end subroutine TracerConsistencyCheck

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, &
       num_soilp, filter_soilp, &
       num_allc, filter_allc)
    !
    ! !DESCRIPTION:
    ! Compute end-of-timestep summaries of water diagnostic terms
    !
    ! !ARGUMENTS:
    class(water_type) , intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds
    integer           , intent(in)    :: num_soilp       ! number of patches in soilp filter
    integer           , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer           , intent(in)    :: num_allc        ! number of columns in allc filter
    integer           , intent(in)    :: filter_allc(:)  ! filter for all columns
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'Summary'
    !-----------------------------------------------------------------------

    do i = this%bulk_and_tracers_beg, this%bulk_and_tracers_end
       associate(bulk_or_tracer => this%bulk_and_tracers(i))
       call bulk_or_tracer%waterdiagnostic_inst%Summary( &
            bounds = bounds, &
            num_soilp = num_soilp, &
            filter_soilp = filter_soilp, &
            num_allc = num_allc, &
            filter_allc = filter_allc, &
            waterstate_inst = bulk_or_tracer%waterstate_inst, &
            waterflux_inst = bulk_or_tracer%waterflux_inst)
       end associate
    end do

  end subroutine Summary


end module WaterType
