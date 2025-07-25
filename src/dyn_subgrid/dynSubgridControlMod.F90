module dynSubgridControlMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Defines a class for storing and querying control flags related to dynamic subgrid
  ! operation.
  !
  ! Note that this is implemented (essentially) as a singleton, so the only instance of
  ! this class is stored in this module. This is done for convenience, to avoid having to
  ! pass around the single instance just to query these control flags.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use abortutils         , only : endrun
  use clm_varctl         , only : fname_len
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: dynSubgridControl_init
  public :: get_flanduse_timeseries ! return the value of the flanduse_timeseries file name
  public :: get_do_transient_pfts   ! return the value of the do_transient_pfts control flag
  public :: get_do_transient_crops  ! return the value of the do_transient_crops control flag
  public :: get_do_transient_lakes  ! return the value of the do_transient_lakes control flag
  public :: get_do_transient_urban  ! return the value of the do_transient_urban control flag
  public :: get_vars_1dwt_w_time  ! return the value of the vars_1dwt_w_time control flag
  public :: run_has_transient_landcover ! returns true if any aspects of prescribed transient landcover are enabled
  public :: get_do_harvest          ! return the value of the do_harvest control flag
  public :: get_do_grossunrep       ! return the value of the do_grossunrep control flag
  public :: get_reset_dynbal_baselines ! return the value of the reset_dynbal_baselines control flag
  public :: get_for_testing_allow_non_annual_changes ! return true if user has requested to allow area changes at times other than the year boundary, for testing purposes
  public :: get_for_testing_zero_dynbal_fluxes ! return true if user has requested to set the dynbal water and energy fluxes to zero, for testing purposes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: read_namelist              ! read namelist variables
  private :: check_namelist_consistency ! check consistency of namelist settings
  !
  ! !PRIVATE TYPES:
  type dyn_subgrid_control_type
     private
     character(len=fname_len) :: flanduse_timeseries = ' ' ! transient landuse dataset
     logical :: do_transient_pfts  = .false. ! whether to apply transient natural PFTs from dataset
     logical :: do_transient_crops = .false. ! whether to apply transient crops from dataset
     logical :: do_transient_lakes = .false. ! whether to apply transient lakes from dataset 
     logical :: do_transient_urban = .false. ! whether to apply transient urban from dataset
     logical :: do_harvest         = .false. ! whether to apply harvest from dataset
     logical :: do_grossunrep      = .false. ! whether to apply gross unrepresented landcover change from dataset
     logical :: vars_1dwt_w_time = .false. ! whether to add the time dimension to 1dwt variables, e.g. pfts1d_wtcol

     logical :: reset_dynbal_baselines = .false. ! whether to reset baseline values of total column water and energy in the first step of the run

     ! The following is only meant for testing: Whether area changes are allowed at times
     ! other than the year boundary. This should only arise in some test configurations
     ! where we artificially create changes more frequently so that we can run short
     ! tests. This flag is only used for error-checking, not controlling any model
     ! behavior.
     logical :: for_testing_allow_non_annual_changes = .false.

     ! The following is only meant for testing: If .true., set the dynbal water and
     ! energy fluxes to zero. This is needed in some tests where we have daily rather
     ! than annual glacier dynamics: if we allow the true dynbal adjustment fluxes in
     ! those tests, we end up with sensible heat fluxes of thousands of W m-2 or more,
     ! which causes CAM to blow up. However, note that setting it to true will break
     ! water and energy conservation!
     logical :: for_testing_zero_dynbal_fluxes = .false.

     logical :: initialized        = .false. ! whether this object has been initialized
  end type dyn_subgrid_control_type
  
  type(dyn_subgrid_control_type) :: dyn_subgrid_control_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  
  !-----------------------------------------------------------------------
  subroutine dynSubgridControl_init( NLFilename )
    !
    ! !DESCRIPTION:
    ! Initialize the dyn_subgrid_control settings.
    !
    ! !USES:
    use spmdMod           , only : masterproc
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'dynSubgridControl_init'
    !-----------------------------------------------------------------------
    
    call read_namelist( NLFilename )
    if (masterproc) then
       call check_namelist_consistency
    end if

    dyn_subgrid_control_inst%initialized = .true.

  end subroutine dynSubgridControl_init

  !-----------------------------------------------------------------------
  subroutine read_namelist( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read dyn_subgrid_control namelist variables
    !
    ! !USES:
    use fileutils      , only : getavu, relavu
    use clm_nlUtilsMod , only : find_nlgroup_name
    use clm_varctl     , only : iulog
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename   ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    ! temporary variables corresponding to the components of dyn_subgrid_control_type:
    character(len=fname_len) :: flanduse_timeseries
    logical :: do_transient_pfts
    logical :: do_transient_crops
    logical :: do_transient_lakes
    logical :: do_transient_urban
    logical :: do_harvest
    logical :: do_grossunrep
    logical :: vars_1dwt_w_time
    logical :: reset_dynbal_baselines
    logical :: for_testing_allow_non_annual_changes
    logical :: for_testing_zero_dynbal_fluxes
    ! other local variables:
    integer :: nu_nml    ! unit for namelist file
    integer :: nml_error ! namelist i/o error flag
    
    character(len=*), parameter :: subname = 'read_namelist'
    !-----------------------------------------------------------------------
    
    namelist /dynamic_subgrid/ &
         flanduse_timeseries, &
         do_transient_pfts, &
         do_transient_crops, &
         do_transient_lakes, &
         do_transient_urban, &
         do_harvest, &
         do_grossunrep, &
         vars_1dwt_w_time, &
         reset_dynbal_baselines, &
         for_testing_allow_non_annual_changes, &
         for_testing_zero_dynbal_fluxes

    ! Initialize options to default values, in case they are not specified in the namelist
    flanduse_timeseries = ' '
    do_transient_pfts  = .false.
    do_transient_crops = .false.
    do_transient_lakes = .false.
    do_transient_urban = .false.
    do_harvest         = .false.
    do_grossunrep      = .false.
    vars_1dwt_w_time = .false.
    reset_dynbal_baselines = .false.
    for_testing_allow_non_annual_changes = .false.
    for_testing_zero_dynbal_fluxes = .false.

    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'dynamic_subgrid', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=dynamic_subgrid, iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading dynamic_subgrid namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg='ERROR finding dynamic_subgrid namelist'//errMsg(sourcefile, __LINE__))
       end if
       close(nu_nml)
       call relavu( nu_nml )
    endif

    call shr_mpi_bcast (flanduse_timeseries, mpicom)
    call shr_mpi_bcast (do_transient_pfts, mpicom)
    call shr_mpi_bcast (do_transient_crops, mpicom)
    call shr_mpi_bcast (do_transient_lakes, mpicom)
    call shr_mpi_bcast (do_transient_urban, mpicom)
    call shr_mpi_bcast (do_harvest, mpicom)
    call shr_mpi_bcast (do_grossunrep, mpicom)
    call shr_mpi_bcast (vars_1dwt_w_time, mpicom)
    call shr_mpi_bcast (reset_dynbal_baselines, mpicom)
    call shr_mpi_bcast (for_testing_allow_non_annual_changes, mpicom)
    call shr_mpi_bcast (for_testing_zero_dynbal_fluxes, mpicom)

    dyn_subgrid_control_inst = dyn_subgrid_control_type( &
         flanduse_timeseries = flanduse_timeseries, &
         do_transient_pfts = do_transient_pfts, &
         do_transient_crops = do_transient_crops, &
         do_transient_lakes = do_transient_lakes, &
         do_transient_urban = do_transient_urban, &
         do_harvest = do_harvest, &
         do_grossunrep = do_grossunrep, &
         vars_1dwt_w_time = vars_1dwt_w_time, &
         reset_dynbal_baselines = reset_dynbal_baselines, &
         for_testing_allow_non_annual_changes = for_testing_allow_non_annual_changes, &
         for_testing_zero_dynbal_fluxes = for_testing_zero_dynbal_fluxes)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'dynamic_subgrid settings:'
       write(iulog,nml=dynamic_subgrid)
       write(iulog,*) ' '
    end if

  end subroutine read_namelist

  !-----------------------------------------------------------------------
  subroutine check_namelist_consistency
    !
    ! !DESCRIPTION:
    ! Check consistency of namelist settings
    !
    ! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl, only : iulog, use_cndv, use_fates, use_cn, use_crop, &
                           n_dom_pfts, n_dom_landunits, collapse_urban, &
                           toosmall_soil, toosmall_crop, toosmall_glacier, &
                           toosmall_lake, toosmall_wetland, toosmall_urban
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'check_namelist_consistency'
    !-----------------------------------------------------------------------
    
    if (dyn_subgrid_control_inst%flanduse_timeseries == ' ') then
       if (dyn_subgrid_control_inst%do_transient_pfts) then
          write(iulog,*) 'ERROR: do_transient_pfts can only be true if you are running with'
          write(iulog,*) 'a flanduse_timeseries file (currently flanduse_timeseries is blank)'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (dyn_subgrid_control_inst%do_transient_crops) then
          write(iulog,*) 'ERROR: do_transient_crops can only be true if you are running with'
          write(iulog,*) 'a flanduse_timeseries file (currently flanduse_timeseries is blank)'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (dyn_subgrid_control_inst%do_transient_lakes) then
          write(iulog,*) 'ERROR: do_transient_lakes can only be true if you are running with'
          write(iulog,*) 'a flanduse_timeseries file (currently flanduse_timeseries is blank)'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (dyn_subgrid_control_inst%do_transient_urban) then
          write(iulog,*) 'ERROR: do_transient_urban can only be true if you are running with'
          write(iulog,*) 'a flanduse_timeseries file (currently flanduse_timeseries is blank)'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (dyn_subgrid_control_inst%do_harvest) then
          write(iulog,*) 'ERROR: do_harvest can only be true if you are running with'
          write(iulog,*) 'a flanduse_timeseries file (currently flanduse_timeseries is blank)'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (dyn_subgrid_control_inst%do_grossunrep) then
          write(iulog,*) 'ERROR: do_grossunrep can only be true if you are running with'
          write(iulog,*) 'a flanduse_timeseries file (currently flanduse_timeseries is blank)'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    if (dyn_subgrid_control_inst%do_transient_pfts) then
       if (use_cndv) then
          write(iulog,*) 'ERROR: do_transient_pfts is incompatible with use_cndv'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (use_fates) then
          write(iulog,*) 'ERROR: do_transient_pfts is incompatible with use_fates'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    ! NOTE(wjs, 2020-08-23) In the following error checks, I'm treating do_transient_lakes
    ! similar to do_transient_pfts and do_transient_crops. I'm not sure if all of these
    ! checks are truly important for transient lakes (in particular, my guess is that
    ! collapse_urban could probably be done with transient lakes - as well as transient
    ! pfts and transient crops for that matter), but some of the checks probably are
    ! needed, and it seems best to keep transient lakes consistent with other transient
    ! areas in this respect.
    if (dyn_subgrid_control_inst%do_transient_pfts .or. &
         dyn_subgrid_control_inst%do_transient_crops .or. &
         dyn_subgrid_control_inst%do_transient_lakes .or. &
         dyn_subgrid_control_inst%do_transient_urban) then
       if (collapse_urban) then
          write(iulog,*) 'ERROR: do_transient_pfts, do_transient_crops, do_transient_lakes and&
                          do_transient_urban are incompatible with collapse_urban = .true.'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (n_dom_pfts > 0 .or. n_dom_landunits > 0 &
           .or. toosmall_soil > 0._r8 .or. toosmall_crop > 0._r8 &
           .or. toosmall_glacier > 0._r8 .or. toosmall_lake > 0._r8 &
           .or. toosmall_wetland > 0._r8 .or. toosmall_urban > 0._r8) then
          write(iulog,*) 'ERROR: do_transient_pfts, do_transient_crops and do_transient_lakes and &
                          do_transient_urban are incompatible with any of the following set to > 0: &
                          n_dom_pfts > 0, n_dom_landunits > 0, &
                          toosmall_soil > 0._r8, toosmall_crop > 0._r8, &
                          toosmall_glacier > 0._r8, toosmall_lake > 0._r8, &
                          toosmall_wetland > 0._r8, toosmall_urban > 0._r8.'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    if (dyn_subgrid_control_inst%do_transient_crops) then
       if (use_fates) then
          ! NOTE(wjs, 2017-01-13) ED / FATES does not currently have a mechanism for
          ! changing its column areas, with the consequent changes in aboveground biomass
          ! per unit area. See https://github.com/NGEET/ed-clm/issues/173
          write(iulog,*) 'ERROR: do_transient_crops does not currently work with use_fates'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    if (dyn_subgrid_control_inst%do_harvest) then
       if (.not. (use_cn .or. use_fates)) then
          write(iulog,*) 'ERROR: do_harvest can only be true if either use_cn or use_fates are true'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

    if (dyn_subgrid_control_inst%do_grossunrep) then
       ! First check if use_fates. In this case the .not. use_cn error will not
       ! appear. The .not. use_cn error will appea
       ! if .not. use_fates and .not. use_cn.
       if (use_fates) then
          write(iulog,*) 'ERROR: do_grossunrep currently does not work with use_fates'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (.not. use_cn) then
          write(iulog,*) 'ERROR: do_grossunrep can only be true if use_cn is true'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine check_namelist_consistency

  !-----------------------------------------------------------------------
  character(len=fname_len) function get_flanduse_timeseries()
    ! !DESCRIPTION:
    ! Return the value of the flanduse_timeseries file name

    character(len=*), parameter :: subname = 'get_flanduse_timeseries'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_flanduse_timeseries = dyn_subgrid_control_inst%flanduse_timeseries

  end function get_flanduse_timeseries

  !-----------------------------------------------------------------------
  logical function get_do_transient_pfts()
    ! !DESCRIPTION:
    ! Return the value of the do_transient_pfts control flag
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_do_transient_pfts = dyn_subgrid_control_inst%do_transient_pfts

  end function get_do_transient_pfts

  !-----------------------------------------------------------------------
  logical function get_do_transient_crops()
    ! !DESCRIPTION:
    ! Return the value of the do_transient_crops control flag
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_do_transient_crops = dyn_subgrid_control_inst%do_transient_crops

  end function get_do_transient_crops
  
  !-----------------------------------------------------------------------
  logical function get_do_transient_lakes()
    ! !DESCRIPTION:
    ! Return the value of the do_transient_lakes control flag
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_do_transient_lakes = dyn_subgrid_control_inst%do_transient_lakes

  end function get_do_transient_lakes

  !-----------------------------------------------------------------------
  logical function get_do_transient_urban()
    ! !DESCRIPTION:
    ! Return the value of the do_transient_urban control flag
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_do_transient_urban = dyn_subgrid_control_inst%do_transient_urban

  end function get_do_transient_urban
  
  !-----------------------------------------------------------------------
  logical function get_vars_1dwt_w_time()
    ! !DESCRIPTION:
    ! Return the value of the vars_1dwt_w_time control flag
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_vars_1dwt_w_time = dyn_subgrid_control_inst%vars_1dwt_w_time

  end function get_vars_1dwt_w_time

  !-----------------------------------------------------------------------
  logical function run_has_transient_landcover()
    ! !DESCRIPTION:
    ! Returns true if any aspects of prescribed transient landcover are enabled
    !-----------------------------------------------------------------------

    run_has_transient_landcover = &
         (get_do_transient_pfts() .or. &
         get_do_transient_crops() .or. &
         get_do_transient_lakes() .or. &
         get_do_transient_urban())
  end function run_has_transient_landcover

  !-----------------------------------------------------------------------
  logical function get_do_harvest()
    ! !DESCRIPTION:
    ! Return the value of the do_harvest control flag
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_do_harvest = dyn_subgrid_control_inst%do_harvest

  end function get_do_harvest

  !-----------------------------------------------------------------------
  logical function get_do_grossunrep()
    ! !DESCRIPTION:
    ! Return the value of the do_grossunrep control flag
    !-----------------------------------------------------------------------
    
    SHR_ASSERT(dyn_subgrid_control_inst%initialized, errMsg(sourcefile, __LINE__))

    get_do_grossunrep = dyn_subgrid_control_inst%do_grossunrep

  end function get_do_grossunrep

  !-----------------------------------------------------------------------
  logical function get_reset_dynbal_baselines()
    ! !DESCRIPTION:
    ! Return the value of the reset_dynbal_baselines control flag
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_reset_dynbal_baselines = dyn_subgrid_control_inst%reset_dynbal_baselines

  end function get_reset_dynbal_baselines

  !-----------------------------------------------------------------------
  logical function get_for_testing_allow_non_annual_changes()
    !
    ! !DESCRIPTION:
    ! Return true if the user has requested to allow area changes at times other than the
    ! year boundary. (This should typically only be true for testing.) (This only
    ! controls error-checking, not any operation of the code.)
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_for_testing_allow_non_annual_changes = dyn_subgrid_control_inst%for_testing_allow_non_annual_changes

  end function get_for_testing_allow_non_annual_changes

  !-----------------------------------------------------------------------
  logical function get_for_testing_zero_dynbal_fluxes()
    !
    ! !DESCRIPTION:
    ! Return true if the user has requested to set the dynbal water and energy fluxes to
    ! zero. This should typically only be true for testing: This is needed in some tests
    ! where we have daily rather than annual glacier dynamics: if we allow the true dynbal
    ! adjustment fluxes in those tests, we end up with sensible heat fluxes of thousands
    ! of W m-2 or more, which causes CAM to blow up. However, note that setting it to
    ! true will break water and energy conservation!
    ! -----------------------------------------------------------------------

    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_for_testing_zero_dynbal_fluxes = dyn_subgrid_control_inst%for_testing_zero_dynbal_fluxes

  end function get_for_testing_zero_dynbal_fluxes

end module dynSubgridControlMod
