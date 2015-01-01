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
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use abortutils         , only : endrun
  use clm_varctl         , only : fname_len
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: dynSubgridControl_init
  public :: get_flanduse_timeseries     ! return the value of the flanduse_timeseries file name
  public :: get_do_transient_pfts       ! return the value of the do_transient_pfts control flag
  public :: get_do_transient_crops      ! return the value of the do_transient_crops control flag
  public :: get_do_harvest              ! return the value of the do_harvest control flag
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
     logical :: do_harvest         = .false. ! whether to apply harvest from dataset
  end type dyn_subgrid_control_type
  
  type(dyn_subgrid_control_type) :: dyn_subgrid_control_inst

contains
  
  !-----------------------------------------------------------------------
  subroutine dynSubgridControl_init
    !
    ! !DESCRIPTION:
    ! Initialize the dyn_subgrid_control settings.
    !
    ! !USES:
    use spmdMod           , only : masterproc
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'dynSubgridControl_init'
    !-----------------------------------------------------------------------
    
    call read_namelist
    if (masterproc) then
       call check_namelist_consistency
    end if

  end subroutine dynSubgridControl_init

  !-----------------------------------------------------------------------
  subroutine read_namelist
    !
    ! !DESCRIPTION:
    ! Read dyn_subgrid_control namelist variables
    !
    ! !USES:
    use fileutils      , only : getavu, relavu
    use clm_nlUtilsMod , only : find_nlgroup_name
    use controlMod     , only : NLFilename
    use clm_varctl     , only : iulog
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    ! temporary variables corresponding to the components of dyn_subgrid_control_type:
    character(len=fname_len) :: flanduse_timeseries
    logical :: do_transient_pfts
    logical :: do_transient_crops
    logical :: do_harvest
    ! other local variables:
    integer :: nu_nml    ! unit for namelist file
    integer :: nml_error ! namelist i/o error flag
    
    character(len=*), parameter :: subname = 'read_namelist'
    !-----------------------------------------------------------------------
    
    namelist /dynamic_subgrid/ &
         flanduse_timeseries, &
         do_transient_pfts, &
         do_transient_crops, &
         do_harvest

    ! Initialize options to default values, in case they are not specified in the namelist
    flanduse_timeseries = ' '
    do_transient_pfts  = .false.
    do_transient_crops = .false.
    do_harvest         = .false.

    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'dynamic_subgrid', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=dynamic_subgrid, iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading dynamic_subgrid namelist'//errMsg(__FILE__, __LINE__))
          end if
       end if
       close(nu_nml)
       call relavu( nu_nml )
    endif

    call shr_mpi_bcast (flanduse_timeseries, mpicom)
    call shr_mpi_bcast (do_transient_pfts, mpicom)
    call shr_mpi_bcast (do_transient_crops, mpicom)
    call shr_mpi_bcast (do_harvest, mpicom)

    dyn_subgrid_control_inst = dyn_subgrid_control_type( &
         flanduse_timeseries = flanduse_timeseries, &
         do_transient_pfts = do_transient_pfts, &
         do_transient_crops = do_transient_crops, &
         do_harvest = do_harvest)

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
    ! Check consistency of namelist settingsn
    !
    ! !USES:
    use clm_varctl     , only : iulog, use_cndv, use_ed, use_crop, use_cn
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
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       if (dyn_subgrid_control_inst%do_transient_crops) then
          write(iulog,*) 'ERROR: do_transient_crops can only be true if you are running with'
          write(iulog,*) 'a flanduse_timeseries file (currently flanduse_timeseries is blank)'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       if (dyn_subgrid_control_inst%do_harvest) then
          write(iulog,*) 'ERROR: do_harvest can only be true if you are running with'
          write(iulog,*) 'a flanduse_timeseries file (currently flanduse_timeseries is blank)'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end if

    if (dyn_subgrid_control_inst%do_transient_pfts) then
       if (use_cndv) then
          write(iulog,*) 'ERROR: do_transient_pfts is incompatible with use_cndv'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       if (use_ed) then
          write(iulog,*) 'ERROR: do_transient_pfts is incompatible with use_ed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end if

    if (dyn_subgrid_control_inst%do_transient_crops) then
       if (.not. use_crop) then
          write(iulog,*) 'ERROR: do_transient_crops can only be true if use_crop is true'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       if (use_ed) then
          write(iulog,*) 'ERROR: do_transient_crops has not been tested with use_ed;'
          write(iulog,*) 'for now these two options cannot be combined'
       end if
    end if

    if (dyn_subgrid_control_inst%do_harvest) then
       if (.not. use_cn) then
          write(iulog,*) 'ERROR: do_harvest can only be true if use_cn is true'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       if (use_ed) then
          write(iulog,*) 'ERROR: do_harvest currently does not work with use_ed'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end if

  end subroutine check_namelist_consistency

  !-----------------------------------------------------------------------
  character(len=fname_len) function get_flanduse_timeseries()
    ! !DESCRIPTION:
    ! Return the value of the flanduse_timeseries file name
    !-----------------------------------------------------------------------
    
    get_flanduse_timeseries = dyn_subgrid_control_inst%flanduse_timeseries

  end function get_flanduse_timeseries

  !-----------------------------------------------------------------------
  logical function get_do_transient_pfts()
    ! !DESCRIPTION:
    ! Return the value of the do_transient_pfts control flag
    !-----------------------------------------------------------------------
    
    get_do_transient_pfts = dyn_subgrid_control_inst%do_transient_pfts

  end function get_do_transient_pfts

  !-----------------------------------------------------------------------
  logical function get_do_transient_crops()
    ! !DESCRIPTION:
    ! Return the value of the do_transient_crops control flag
    !-----------------------------------------------------------------------
    
    get_do_transient_crops = dyn_subgrid_control_inst%do_transient_crops

  end function get_do_transient_crops

  !-----------------------------------------------------------------------
  logical function get_do_harvest()
    ! !DESCRIPTION:
    ! Return the value of the do_harvest control flag
    !-----------------------------------------------------------------------
    
    get_do_harvest = dyn_subgrid_control_inst%do_harvest

  end function get_do_harvest

end module dynSubgridControlMod
