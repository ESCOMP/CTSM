module WaterType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Container for derived types relating to water, both for bulk water and for isotopes
  ! and other tracers.
  !
  ! To loop through all tracers, use code like this:
  !    do i = 1, water_inst%num_tracers
  !       call some_subroutine(..., water_inst%waterflux_tracer_inst(i), ...)
  !    end do
  !
  ! To loop through all isotopes, use code like this:
  !    type(water_info_isotope_type), pointer :: iso_info
  !
  !    do i = 1, water_inst%num_tracers
  !       if (water_inst%IsIsotope(i)) then
  !          call water_inst%GetIsotopeInfo(i, iso_info)
  !          call some_subroutine(..., iso_info, water_inst%waterflux_tracer_inst(i), ...)
  !       end if
  !    end do
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use abortutils              , only : endrun
  use decompMod               , only : bounds_type
  use clm_varctl              , only : iulog
  use clm_varpar              , only : nlevsno
  use ncdio_pio               , only : file_desc_t
  use WaterFluxBulkType       , only : waterfluxbulk_type
  use WaterFluxType           , only : waterflux_type
  use WaterStateBulkType      , only : waterstatebulk_type
  use WaterStateType          , only : waterstate_type
  use WaterDiagnosticType     , only : waterdiagnostic_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use WaterBalanceType        , only : waterbalance_type
  use WaterInfoBaseType       , only : water_info_base_type
  use WaterInfoBulkType       , only : water_info_bulk_type
  use WaterInfoIsotopeType    , only : water_info_isotope_type
  use Waterlnd2atmType    , only : waterlnd2atm_type
  use Waterlnd2atmBulkType    , only : waterlnd2atmbulk_type
  use Wateratm2lndType    , only : wateratm2lnd_type
  use Wateratm2lndBulkType    , only : wateratm2lndbulk_type

  implicit none
  private

  !
  ! !PRIVATE TYPES:

  type, private :: water_params_type
     private

     ! Whether we add tracers that are used for the tracer consistency checks
     logical :: enable_consistency_checks

     ! Whether we add tracers that are used for isotopes
     logical :: enable_isotopes
  end type water_params_type

  ! This type is a container for objects of class water_info_base_type, to facilitate
  ! having an array of polymorphic entities.
  type, private :: water_info_container_type
     private
     ! 'info' needs to be a pointer so other pointers can point to it (since a derived
     ! type component cannot have the target attribute)
     class(water_info_base_type), pointer :: info
  end type water_info_container_type

  !
  ! !PUBLIC TYPES:
  type, public :: water_type
     private

     ! ------------------------------------------------------------------------
     ! Public data members
     ! ------------------------------------------------------------------------

     integer, public :: num_tracers

     type(waterfluxbulk_type), public       :: waterfluxbulk_inst
     type(waterstatebulk_type), public      :: waterstatebulk_inst
     type(waterdiagnosticbulk_type), public :: waterdiagnosticbulk_inst
     type(waterbalance_type), public        :: waterbalancebulk_inst
     type(waterlnd2atmbulk_type), public        :: waterlnd2atmbulk_inst
     type(wateratm2lndbulk_type), public        :: wateratm2lndbulk_inst

     type(waterflux_type), allocatable, public       :: waterflux_tracer_inst(:)
     type(waterstate_type), allocatable, public      :: waterstate_tracer_inst(:)
     type(waterdiagnostic_type), allocatable, public :: waterdiagnostic_tracer_inst(:)
     type(waterbalance_type), allocatable, public    :: waterbalance_tracer_inst(:)
     type(waterlnd2atm_type), allocatable, public :: waterlnd2atm_tracer_inst(:)
     type(wateratm2lnd_type), allocatable, public :: wateratm2lnd_tracer_inst(:)

     ! ------------------------------------------------------------------------
     ! Private data members
     ! ------------------------------------------------------------------------

     type(water_params_type) :: params

     ! bulk_info needs to be a pointer so other pointers can point to it (since a derived
     ! type component cannot have the target attribute)
     type(water_info_bulk_type), pointer :: bulk_info
     logical, allocatable :: is_isotope(:)
     type(water_info_container_type), allocatable :: tracer_info(:)
     integer :: bulk_tracer_index  ! index of the tracer that replicates bulk water (-1 if it doesn't exist)

   contains
     ! Public routines
     procedure, public :: Init
     procedure, public :: InitAccBuffer
     procedure, public :: InitAccVars
     procedure, public :: UpdateAccVars
     procedure, public :: Restart
     procedure, public :: IsIsotope       ! Return true if a given tracer is an isotope
     procedure, public :: GetIsotopeInfo  ! Get a pointer to the object storing isotope info for a given tracer
     procedure, public :: GetBulkTracerIndex ! Get the index of the tracer that replicates bulk water
     procedure, public :: DoConsistencyCheck ! Whether TracerConsistencyCheck should be called in this run
     procedure, public :: TracerConsistencyCheck

     ! Private routines
     procedure, private :: ReadNamelist
     procedure, private :: SetupTracerInfo
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
       h2osno_col, snow_depth_col, watsat_col, t_soisno_col)
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
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: i

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    SHR_ASSERT_ALL((ubound(h2osno_col) == [endc]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(snow_depth_col) == [endc]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(watsat_col, 1) == endc), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno_col, 1) == endc), errMsg(sourcefile, __LINE__))

    call this%ReadNamelist(NLFilename)

    allocate(this%bulk_info, source = water_info_bulk_type())

    call this%waterstatebulk_inst%InitBulk(bounds, &
         this%bulk_info, &
         h2osno_input_col = h2osno_col(begc:endc),       &
         watsat_col = watsat_col(begc:endc, 1:),   &
         t_soisno_col = t_soisno_col(begc:endc, -nlevsno+1:) )

    call this%waterdiagnosticbulk_inst%InitBulk(bounds, &
         this%bulk_info, &
         snow_depth_input_col = snow_depth_col(begc:endc),    &
         waterstatebulk_inst = this%waterstatebulk_inst )

    call this%waterbalancebulk_inst%Init(bounds, &
         this%bulk_info)

    call this%waterfluxbulk_inst%InitBulk(bounds, &
         this%bulk_info)

    call this%waterlnd2atmbulk_inst%InitBulk(bounds, &
         this%bulk_info)

    call this%wateratm2lndbulk_inst%InitBulk(bounds, &
         this%bulk_info)

    call this%SetupTracerInfo()

    if (this%num_tracers > 0) then
       allocate(this%waterflux_tracer_inst(this%num_tracers))
       allocate(this%waterstate_tracer_inst(this%num_tracers))
       allocate(this%waterdiagnostic_tracer_inst(this%num_tracers))
       allocate(this%waterbalance_tracer_inst(this%num_tracers))
       allocate(this%waterlnd2atm_tracer_inst(this%num_tracers))
       allocate(this%wateratm2lnd_tracer_inst(this%num_tracers))
    end if

    do i = 1, this%num_tracers

       call this%waterstate_tracer_inst(i)%Init(bounds, &
            this%tracer_info(i)%info, &
            h2osno_input_col = h2osno_col(begc:endc),       &
            watsat_col = watsat_col(begc:endc, 1:),   &
            t_soisno_col = t_soisno_col(begc:endc, -nlevsno+1:) )

       call this%waterdiagnostic_tracer_inst(i)%Init(bounds, &
            this%tracer_info(i)%info)

       call this%waterbalance_tracer_inst(i)%Init(bounds, &
            this%tracer_info(i)%info)

       call this%waterflux_tracer_inst(i)%Init(bounds, &
            this%tracer_info(i)%info)

       call this%waterlnd2atm_tracer_inst(i)%Init(bounds, &
            this%tracer_info(i)%info)

       call this%wateratm2lnd_tracer_inst(i)%Init(bounds, &
            this%tracer_info(i)%info)

    end do

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the water_tracers namelist
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

    this%num_tracers = num_tracers
    if (this%num_tracers > 0) then
       allocate(this%tracer_info(this%num_tracers))
       allocate(this%is_isotope(this%num_tracers))
    end if

    tracer_num = 1
    if (enable_bulk_tracer) then
       allocate(this%tracer_info(tracer_num)%info, source = water_info_isotope_type('H2OTR'))
       this%is_isotope(tracer_num) = .true.
       this%bulk_tracer_index = tracer_num
       tracer_num = tracer_num + 1
    end if
    if (this%params%enable_isotopes) then
       allocate(this%tracer_info(tracer_num)%info, source = water_info_isotope_type('HDO'))
       this%is_isotope(tracer_num) = .true.
       tracer_num = tracer_num + 1

       allocate(this%tracer_info(tracer_num)%info, source = water_info_isotope_type('H218O'))
       this%is_isotope(tracer_num) = .true.
       tracer_num = tracer_num + 1
    end if

    if (tracer_num - 1 /= this%num_tracers) then
       write(iulog,*) subname//' ERROR: tracer_num discrepancy'
       write(iulog,*) 'num_tracers = ', this%num_tracers
       write(iulog,*) 'but added ', tracer_num - 1, ' tracers'
       call endrun(msg='tracer_num discrepancy '//errMsg(sourcefile, __LINE__))
    end if

  end subroutine SetupTracerInfo

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

    call this%wateratm2lndbulk_inst%restartBulk (bounds, ncid, flag=flag)

    do i = 1, this%num_tracers

       call this%waterflux_tracer_inst(i)%Restart(bounds, ncid, flag=flag)

       call this%waterstate_tracer_inst(i)%Restart(bounds, ncid, flag=flag, &
            watsat_col=watsat_col(bounds%begc:bounds%endc,:))

       call this%waterdiagnostic_tracer_inst(i)%Restart(bounds, ncid, flag=flag)

       call this%wateratm2lnd_tracer_inst(i)%Restart(bounds, ncid, flag=flag)

    end do

  end subroutine Restart

  !-----------------------------------------------------------------------
  function IsIsotope(this, i)
    !
    ! !DESCRIPTION:
    ! Returns true if tracer i is an isotope
    !
    ! i must be <= this%num_tracers
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical :: IsIsotope  ! function result
    class(water_type), intent(in) :: this
    integer, intent(in) :: i  ! index of tracer
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'IsIsotope'
    !-----------------------------------------------------------------------

    SHR_ASSERT(i <= this%num_tracers, errMsg(sourcefile, __LINE__))

    IsIsotope = this%is_isotope(i)

  end function IsIsotope

  !-----------------------------------------------------------------------
  subroutine GetIsotopeInfo(this, i, isotope_info)
    !
    ! !DESCRIPTION:
    ! Get a pointer to the object storing isotope info for a given tracer
    !
    ! This provides a mechanism for passing the isotope info to subroutines that need it.
    !
    ! i must be <= this%num_tracers, and this%IsIsotope(i) must be true
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

    SHR_ASSERT(i <= this%num_tracers, errMsg(sourcefile, __LINE__))

    select type(info => this%tracer_info(i)%info)
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
    ! !ARGUMENTS:
    class(water_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*) , intent(in) :: caller_location  ! brief description of where this is called from (for error messages)
    !
    ! !LOCAL VARIABLES:
    integer :: i
    character(len=*), parameter :: subname = 'TracerConsistencyCheck'
    !-----------------------------------------------------------------------

    !for now, simply checking bulk vs the bulk tracer, but may eventually 
    !want to loop over all tracers in some situations
    i = this%GetBulkTracerIndex()
    if (i < 0) then
       write(iulog,*) subname, ' Error: requesting tracer consistency check with bulk tracer not enabled'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    call this%waterbalance_tracer_inst(i)%TracerConsistencyCheck(bounds, &
         this%waterbalancebulk_inst, caller_location)

    call this%waterstate_tracer_inst(i)%TracerConsistencyCheck(bounds, &
         this%waterstatebulk_inst, caller_location)

    call this%waterdiagnostic_tracer_inst(i)%TracerConsistencyCheck(bounds, &
         this%waterdiagnosticbulk_inst, caller_location)

    call this%waterflux_tracer_inst(i)%TracerConsistencyCheck(bounds, &
         this%waterfluxbulk_inst, caller_location)

    call this%waterlnd2atm_tracer_inst(i)%TracerConsistencyCheck(bounds, &
         this%waterlnd2atmbulk_inst, caller_location)

    call this%wateratm2lnd_tracer_inst(i)%TracerConsistencyCheck(bounds, &
         this%wateratm2lndbulk_inst, caller_location)

  end subroutine TracerConsistencyCheck

end module WaterType
