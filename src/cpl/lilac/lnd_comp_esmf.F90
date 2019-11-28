module lnd_comp_esmf

  !----------------------------------------------------------------------------
  ! This is the ESMF cap for CTSM
  !----------------------------------------------------------------------------

  ! External libraries
  use ESMF
  use mpi               , only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE
  use mct_mod           , only : mct_world_init, mct_world_clean, mct_die
  use shr_pio_mod       , only : shr_pio_init1, shr_pio_init2
  use perf_mod          , only : t_startf, t_stopf, t_barrierf

  ! ctsm and share code
  use shr_kind_mod      , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod       , only : shr_sys_abort
  use shr_file_mod      , only : shr_file_setLogUnit, shr_file_getLogUnit 
  use shr_orb_mod       , only : shr_orb_decl, shr_orb_params
  use shr_cal_mod       , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use spmdMod           , only : masterproc, mpicom, spmd_init
  use decompMod         , only : bounds_type, ldecomp, get_proc_bounds
  use domainMod         , only : ldomain
  use controlMod        , only : control_setNL
  use clm_varorb        , only : eccen, obliqr, lambm0, mvelpp
  use clm_varctl        , only : clm_varctl_set, iulog, finidat 
  use clm_varctl        , only : nsrStartup, nsrContinue, nsrBranch
  use clm_varctl        , only : inst_index, inst_suffix, inst_name
  use clm_time_manager  , only : set_timemgr_init, advance_timestep
  use clm_time_manager  , only : set_nextsw_cday, update_rad_dtime
  use clm_time_manager  , only : get_nstep, get_step_size
  use clm_time_manager  , only : get_curr_date, get_curr_calday, set_nextsw_cday
  use clm_initializeMod , only : initialize1, initialize2
  use clm_driver        , only : clm_drv
  use lnd_import_export , only : import_fields, export_fields
  use lnd_shr_methods   , only : chkerr, state_diagnose
  use spmdMod           , only : masterproc, spmd_init
  use glc_elevclass_mod , only : glc_elevclass_init

  implicit none
  private                         ! By default make data private except

  public :: lnd_register          ! register clm initial, run, final methods
  public :: lnd_init              ! clm initialization
  public :: lnd_run               ! clm run phase
  public :: lnd_final             ! clm finalization/cleanup

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer         , parameter    :: dbug_flag = 6
  type(ESMF_Field), public, save :: field

  integer                        :: glc_nec = 10   ! number of glc elevation classes
  integer,          parameter    :: memdebug_level=1
  integer,          parameter    :: dbug = 1
  character(*)    , parameter    :: modName =  "lnd_comp_esmf"
  character(*),     parameter    :: u_FILE_u = &
       __FILE__!
  type(ESMF_Mesh)         :: Emesh, EMeshTemp, lnd_mesh      ! esmf meshes

!===============================================================================
contains
!===============================================================================

  subroutine lnd_register(comp, rc)

    ! Register the clm initial, run, and final phase methods with ESMF.

    ! input/output argumenents
    type(ESMF_GridComp)  :: comp  ! CLM grid component
    integer, intent(out) :: rc    ! return status

    ! local variables
    character(len=*), parameter :: subname=trim(modname)//': [lnd_register] '
    !-----------------------------------------------------------------------

    print *, "in lnd register routine"
    rc = ESMF_SUCCESS
    call ESMF_LogSet ( flush =.true.)
    call ESMF_LogWrite(subname//"lnd gridcompset entry points setting ....!", ESMF_LOGMSG_INFO)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, lnd_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, lnd_run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, lnd_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_LogWrite(subname//"lnd gridcompset entry points finished!", ESMF_LOGMSG_INFO)
  end subroutine lnd_register

  !===============================================================================

  subroutine lnd_init(comp, import_state, export_state, clock, rc)

    ! Initialize land surface model and obtain relevant atmospheric model arrays
    ! back from (i.e. albedos, surface temperature and snow cover over land).

    ! input/output variables
    type(ESMF_GridComp)  :: comp         ! CLM gridded component
    type(ESMF_State)     :: import_state ! CLM import state
    type(ESMF_State)     :: export_state ! CLM export state
    type(ESMF_Clock)     :: clock        ! ESMF synchronization clock
    integer, intent(out) :: rc           ! Return code

    ! local variable
    integer                        :: mpicom_lnd, mpicom_vm, gsize
    type(ESMF_ArraySpec)           :: arrayspec
    type(ESMF_DistGrid)            :: distgrid
    type(ESMF_Array)               :: dom, l2x, x2l
    type(ESMF_VM)                  :: vm
    integer                        :: lsize                           ! size of attribute vector
    integer                        :: g,i,j                           ! indices
    integer                        :: dtime_sync                      ! coupling time-step from the input synchronization clock
    integer                        :: dtime_clm                       ! clm time-step
    logical                        :: exists                          ! true if file exists
    real(r8)                       :: nextsw_cday                     ! calday from clock of next radiation computation
    character(len=CL)              :: caseid                          ! case identifier name
    character(len=CL)              :: ctitle                          ! case description title
    character(len=CL)              :: starttype                       ! start-type (startup, continue, branch, hybrid)
    character(len=CL)              :: calendar                        ! calendar type name
    character(len=CL)              :: hostname                        ! hostname of machine running on
    character(len=CL)              :: version                         ! Model version
    character(len=CL)              :: username                        ! user running the model
    integer                        :: nsrest                          ! clm restart type
    integer                        :: ref_ymd                         ! reference date (YYYYMMDD)
    integer                        :: ref_tod                         ! reference time of day (sec)
    integer                        :: start_ymd                       ! start date (YYYYMMDD)
    integer                        :: start_tod                       ! start time of day (sec)
    integer                        :: stop_ymd                        ! stop date (YYYYMMDD)
    integer                        :: stop_tod                        ! stop time of day (sec)
    logical                        :: brnch_retain_casename           ! flag if should retain the case name on a branch start type
    logical                        :: atm_aero                        ! Flag if aerosol data sent from atm model
    integer                        :: lbnum                           ! input to memory diagnostic
    integer                        :: shrlogunit                      ! old values for log unit and log level
    integer                        :: logunit                         ! original log unit
    type(bounds_type)              :: bounds                          ! bounds
    integer                        :: nfields
    real(R8), pointer              :: fptr(:, :)
    integer                        :: ierr
    integer                        :: ncomps = 1
    integer, pointer               :: comps(:)                        ! array with component ids
    integer, pointer               :: comms(:)                        ! array with mpicoms
    character(len=32), allocatable :: compLabels(:)
    integer,allocatable            :: comp_id(:)                      ! for pio init2
    character(len=64),allocatable  :: comp_name(:)                    ! for pio init2
    integer,allocatable            :: comp_comm(:)                    ! for pio_init2
    logical,allocatable            :: comp_iamin(:)                   ! for pio init2
    integer,allocatable            :: comp_comm_iam(:)                ! for pio_init2
    integer                        :: ymd                             ! CTSM current date (YYYYMMDD)
    integer                        :: orb_iyear_align                 ! associated with model year
    integer                        :: orb_cyear                       ! orbital year for current orbital computation
    integer                        :: orb_iyear                       ! orbital year for current orbital computation
    integer                        :: orb_eccen                       ! orbital year for current orbital computation
    integer                        :: yy, mm ,dd , curr_tod, curr_ymd ! orbital year for current orbital computation
    type(ESMF_Time)                :: currTime                        ! Current time
    type(ESMF_Time)                :: startTime                       ! Start time
    type(ESMF_Time)                :: stopTime                        ! Stop time
    type(ESMF_Time)                :: refTime                         ! Ref time
    type(ESMF_TimeInterval)        :: timeStep
    type(ESMF_Calendar)            :: esmf_calendar                   ! esmf calendar
    type(ESMF_CalKind_Flag)        :: esmf_caltype                    ! esmf calendar type
    integer, pointer               :: gindex(:)                       ! global index space for land and ocean points
    integer, pointer               :: gindex_lnd(:)                   ! global index space for just land points
    integer, pointer               :: gindex_ocn(:)                   ! global index space for just ocean points
    character(ESMF_MAXSTR)         :: cvalue                          ! config data
    integer                        :: nlnd, nocn                      ! local size ofarrays
    integer                        :: n                               ! indices
    integer                        :: year, month, day
    integer                        :: dtime                           ! time step increment (sec)
    type(ESMF_FieldBundle)         :: c2l_fb
    type(ESMF_FieldBundle)         :: l2c_fb
    type(ESMF_State)               :: importState, exportState
    integer                        :: compid                          ! component id
    character(len=32), parameter   :: sub = 'lnd_init'
    character(len=*),  parameter   :: format = "('("//trim(sub)//") :',A)"
    character(len=*), parameter    :: subname=trim(modName)//': [lnd_init_lilac_cap] '
    !------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' is called!', ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    ! Initialize clm MPI communicator
    !------------------------------------------------------------------------
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call ESMF_LogWrite(subname//"ESMF_VMGetCurrent", ESMF_LOGMSG_INFO)
    call ESMF_VMPrint (vm, rc = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    call ESMF_LogWrite(subname//"ESMF_VMGet", ESMF_LOGMSG_INFO)

    ! duplicate the mpi communicator from the current VM
    call MPI_Comm_dup(mpicom_vm, mpicom_lnd, rc)
    call ESMF_LogWrite(subname//"MPI_Comm_dup...", ESMF_LOGMSG_INFO)

!!!! NS : BOTH MPI_INIT and PIO_INIT1 are in lilac_mod.F90


    !------------------------------------------------------------------------
    ! Initialize mct
    ! (needed for data models and cice prescribed capability)
    ! (needed for data model share code - e.g. nitrogen deposition)
    !------------------------------------------------------------------------
    ! TODO: FIX THIS PLEASE!!!!

    allocate(comms(1), comps(1), compLabels(1), comp_iamin(1), comp_comm_iam(1), comp_name(1),stat=ierr)

    comms(1)      = mpicom_lnd !or call MPI_Comm_dup(mpicom_vm, comms(1), ierr)
    comps(1)      = 1
    compLabels(1) = 'LND'
    comp_iamin(1) = .true.
    comp_name(1)  = 'LND'

    call ESMF_VMGet(vm, mpiCommunicator=comms(1), localPet=comp_comm_iam(1), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call shr_pio_init2(comps, compLabels, comp_iamin, comms, comp_comm_iam)

    call ESMF_LogWrite(subname//"after shr_pio_init2", ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//"Now calling mct_world_init", ESMF_LOGMSG_INFO)
    call mct_world_init(ncomps, mpicom_lnd, comms, comps)
    call ESMF_LogWrite(subname//"mct world initialized! ", ESMF_LOGMSG_INFO)

    !deallocate(comms, comps, compLabels, comp_iamin, comp_comm_iam, comp_name)   ???

    ! Initialize model mpi info
    compid = 1
    call spmd_init( mpicom_lnd, compid)
    call ESMF_LogWrite(subname//"initialized model mpi info using spmd_init", ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    !--- Log File ---
    !------------------------------------------------------------------------

    inst_name  = 'LND'; inst_index  = 1; inst_suffix = ""

    ! Initialize io log unit
    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       write(iulog,format) "CLM land model initialization"
    else
       iulog = shrlogunit
    end if
    call shr_file_setLogUnit (iulog)

    !------------------------------------------------------------------------
    !--- Orbital Values ---
    !------------------------------------------------------------------------

    ! TODO: orbital values should be provided by lilac - but for now lets use defaults
    !! hard wire these these in and we can decide on maybe having a namelist/

    !call shr_cal_date2ymd(ymd,year,month,day)
    !orb_cyear = orb_iyear + (year - orb_iyear_align)

    orb_cyear = 2000
    call shr_orb_params(orb_cyear, eccen, obliqr, mvelpp, &
         obliqr, lambm0, mvelpp, masterproc)

    ! for now hard-coding:
    !nextsw_cday =  1.02083333333333
    !eccen       =  1.670366039276560E-002
    !mvelpp      =  4.93745779048816
    !lambm0      =  -3.247249566152933E-0020
    !obliqr      =  0.409101122579779

    !if ((debug >1) .and. (masterproc)) then
    if (masterproc) then
       write(iulog,*) 'shr_obs_params is setting these:', eccen
       write(iulog,*) 'eccen is : ', eccen
       write(iulog,*) 'mvelpp is : ', mvelpp

       write(iulog,*) 'lambm0 is : ', lambm0
       write(iulog,*) 'obliqr is : ', obliqr
    end if

    !----------------------
    ! Consistency check on namelist filename
    !----------------------
    call control_setNL("lnd_in")

    !----------------------
    ! Get properties from clock
    !----------------------


    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    ! TODO: how do we set case_name and nsrest - should we hardwire for now?
    caseid = 'test_lilac'
    nsrest = nsrStartup
    call ESMF_LogWrite(subname//"time manager Initialized....", ESMF_LOGMSG_INFO)

    !----------------------
    ! Initialize CTSM time manager
    !----------------------

    call set_timemgr_init( &
         calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd,         &
         stop_tod_in=stop_tod)
    call ESMF_LogWrite(subname//"time manager is set now!", ESMF_LOGMSG_INFO)

    !----------------------
    ! Read namelist, grid and surface data
    !----------------------

    ! set default values for run control variables
    call clm_varctl_set(caseid_in=caseid, nsrest_in=nsrest)
    call ESMF_LogWrite(subname//"default values for run control variables are set...", ESMF_LOGMSG_INFO)

    !----------------------
    ! Initialize glc_elevclass module
    !----------------------
    call glc_elevclass_init(glc_nec)

    !----------------------
    ! Initialize1
    !----------------------

    ! note that the memory for gindex_ocn will be allocated in the following call
    call initialize1(gindex_ocn)
    ! call initialize1()

    call ESMF_LogWrite(subname//"initialize1 done...", ESMF_LOGMSG_INFO)

    ! obtain global index array for just land points which includes mask=0 or ocean points
    call get_proc_bounds( bounds )

    nlnd = bounds%endg - bounds%begg + 1
    allocate(gindex_lnd(nlnd))
    !print ,* "nlnd is :", nlnd
    do g = bounds%begg,bounds%endg
       n = 1 + (g - bounds%begg)
       gindex_lnd(n) = ldecomp%gdc2glo(g)
    end do

    call ESMF_LogWrite(subname//"obtained global index", ESMF_LOGMSG_INFO)

    ! create a global index that includes both land and ocean points
    nocn = size(gindex_ocn)
    allocate(gindex(nlnd + nocn))
    do n = 1,nlnd+nocn
       if (n <= nlnd) then
          gindex(n) = gindex_lnd(n)
       else
          gindex(n) = gindex_ocn(n-nlnd)
       end if
    end do

    ! create distGrid from global index array
    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    deallocate(gindex)
    call ESMF_LogWrite(subname//"DistGrid created......", ESMF_LOGMSG_INFO)

    !--------------------------------
    ! generate the mesh on ctsm distribution
    !--------------------------------

    ! TODO: mesh file should come into clm as a namelist for lilac only
    ! for now need to hardwire this in cvalue here
    cvalue = '/glade/p/cesmdata/cseg/inputdata/share/meshes/fv4x5_050615_polemod_ESMFmesh.nc' ! this will need to be filled in to run

    EMeshTemp = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    if (masterproc) then
       write(iulog,*)'mesh file for domain is ',trim(cvalue)
    end if

    ! recreate the mesh using the above distGrid
    EMesh = ESMF_MeshCreate(EMeshTemp, elementDistgrid=Distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_LogWrite(subname//" Create Mesh using distgrid ....", ESMF_LOGMSG_INFO)
    lnd_mesh = EMesh
    !--------------------------------
    ! Finish initializing ctsm
    !--------------------------------
    call ESMF_LogWrite(subname//"before initialize2", ESMF_LOGMSG_INFO)

    call initialize2()

    call ESMF_LogWrite(subname//"initialize2 done...", ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Check that ctsm internal dtime aligns with ctsm coupling interval
    !--------------------------------
    call ESMF_LogWrite(subname//"cheking CTSM dtime and coupling intervals....", ESMF_LOGMSG_INFO)

    call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_TimeIntervalGet( timeStep, s=dtime_sync, rc=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    dtime_clm = get_step_size()

    if (masterproc) then
       write(iulog,*)'dtime_sync= ',dtime_sync,' dtime_ctsm= ',dtime_clm,' mod = ',mod(dtime_sync,dtime_clm)
    end if
    if (mod(dtime_sync,dtime_clm) /= 0) then
       write(iulog,*)'ctsm dtime ',dtime_clm,' and clock dtime ',dtime_sync,' never align'
       rc = ESMF_FAILURE
       return
    end if

    !--------------------------------
    ! Create import state (only assume input from atm - not rof and glc)
    !--------------------------------

    c2l_fb =  ESMF_FieldBundleCreate ( name='c2l_fb', rc=rc)


    call fldbundle_add( 'Sa_z'       , c2l_fb,rc)    !1
    call fldbundle_add( 'Sa_topo'    , c2l_fb,rc)    !2
    call fldbundle_add( 'Sa_u'       , c2l_fb,rc)    !3
    call fldbundle_add( 'Sa_v'       , c2l_fb,rc)    !4
    call fldbundle_add( 'Sa_ptem'    , c2l_fb,rc)    !5
    call fldbundle_add( 'Sa_pbot'    , c2l_fb,rc)    !6
    call fldbundle_add( 'Sa_tbot'    , c2l_fb,rc)    !7
    call fldbundle_add( 'Sa_shum'    , c2l_fb,rc)    !8

    call fldbundle_add( 'Faxa_lwdn'  , c2l_fb,rc)    !9
    call fldbundle_add( 'Faxa_rainc'  , c2l_fb,rc)   !10
    call fldbundle_add( 'Faxa_rainl'  , c2l_fb,rc)   !11
    call fldbundle_add( 'Faxa_snowc'  , c2l_fb,rc)   !12
    call fldbundle_add( 'Faxa_snowl'  , c2l_fb,rc)   !13
    call fldbundle_add( 'Faxa_swndr' , c2l_fb,rc)    !14
    call fldbundle_add( 'Faxa_swvdr' , c2l_fb,rc)    !15
    call fldbundle_add( 'Faxa_swndf' , c2l_fb,rc)    !16
    call fldbundle_add( 'Faxa_swvdf' , c2l_fb,rc)    !17
    call ESMF_StateAdd(import_state, fieldbundleList = (/c2l_fb/), rc=rc)

    ! Create export state

    l2c_fb = ESMF_FieldBundleCreate(name='l2c_fb', rc=rc)
    !call fldbundle_add( 'Sl_lfrint'   ,  l2c_fb,rc)   !1
    call fldbundle_add( 'Sl_lfrin'   ,  l2c_fb,rc)   !1
    call fldbundle_add( 'Sl_t'        ,  l2c_fb,rc)   !2
    call fldbundle_add( 'Sl_tref'     ,  l2c_fb,rc)   !3
    call fldbundle_add( 'Sl_qref'     ,  l2c_fb,rc)   !4
    call fldbundle_add( 'Sl_avsdr'    ,  l2c_fb,rc)   !5
    call fldbundle_add( 'Sl_anidr'    ,  l2c_fb,rc)   !6
    call fldbundle_add( 'Sl_avsdf'    ,  l2c_fb,rc)   !7
    call fldbundle_add( 'Sl_anidf'    ,  l2c_fb,rc)   !8
    call fldbundle_add( 'Sl_snowh'    ,  l2c_fb,rc)   !9
    call fldbundle_add( 'Fall_u10'   ,  l2c_fb,rc)   !10
    call fldbundle_add( 'Fall_fv'    ,  l2c_fb,rc)   !11
    call fldbundle_add( 'Fall_ram1'    ,  l2c_fb,rc)   !12
    !call fldbundle_add( 'Fall_taux'   ,  l2c_fb,rc)   !10
    !call fldbundle_add( 'Fall_lwup'   ,  l2c_fb,rc)   !14
    !call fldbundle_add( 'Fall_evap'   ,  l2c_fb,rc)   !15
    !call fldbundle_add( 'Fall_swniet' ,  l2c_fb,rc)   !16
    call ESMF_StateAdd(export_state, fieldbundleList = (/l2c_fb/), rc=rc)
    !call ESMF_StateAdd(exportState, fieldbundleList = (/l2c_fb/), rc=rc)

    !--------------------------------
    ! Create land export state
    !--------------------------------
    call ESMF_LogWrite(subname//"Creating land export state", ESMF_LOGMSG_INFO)

    ! Fill in export state at end of initialization
    call export_fields(comp, bounds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//"Getting Calendar Day of nextsw calculation...", ESMF_LOGMSG_INFO)

    ! Get calendar day of next sw (shortwave) calculation (nextsw_cday)
    if (nsrest == nsrStartup) then
       call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

       call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    end if

    ! Set nextsw_cday
    call set_nextsw_cday( nextsw_cday )

    if (masterproc) then
       write(iulog,*) 'TimeGet ... nextsw_cday is : ', nextsw_cday
    end if

    ! Set Attributes
    call ESMF_LogWrite(subname//"setting attribute!", ESMF_LOGMSG_INFO)

    call ESMF_AttributeSet(export_state, name="lnd_nx", value=ldomain%ni, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_LogWrite(subname//"setting attribute! lnd_nx", ESMF_LOGMSG_INFO)

    call ESMF_AttributeSet(export_state, name="lnd_ny", value=ldomain%nj, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_LogWrite(subname//"setting attribute-lnd_ny!", ESMF_LOGMSG_INFO)

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call State_diagnose(export_state, subname//':ExportState',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_int:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    call ESMF_LogWrite(subname//' CTSM INITIALIZATION DONE SUCCESSFULLY!!!! ', ESMF_LOGMSG_INFO)

  end subroutine lnd_init

  !---------------------------------------------------------------------------

  !subroutine fldbundle_add(stdname, fldptr, fieldbundle,rc)
  subroutine fldbundle_add(stdname, fieldbundle,rc)
    type(ESMF_Field)                   :: field
    !type(ESMF_Mesh)                    :: lnd_mesh
    character(len=*),            intent(in)    :: stdname
    type (ESMF_FieldBundle)      :: fieldbundle
    integer, intent(out) :: rc

    print *, "in lnd register routine"

    rc = ESMF_SUCCESS

    !field = ESMF_FieldCreate(lnd_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(stdname), rc=rc)
    field = ESMF_FieldCreate(lnd_mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT , name=trim(stdname), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_FieldBundleAdd(fieldbundle, (/field/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
  end subroutine fldbundle_add

  !---------------------------------------------------------------------------
  subroutine lnd_run(gcomp, import_state, export_state, clock, rc)

    !------------------------
    ! Run CTSM
    !------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp           ! CLM gridded component
    type(ESMF_State)     :: import_state    ! CLM import state
    type(ESMF_State)     :: export_state    ! CLM export state
    type(ESMF_Clock)     :: clock           ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code

    ! local variables:
    type(ESMF_Alarm)       :: alarm
    type(ESMF_Time)        :: currTime
    type(ESMF_Time)        :: nextTime
    type(ESMF_State)       :: importState, exportState
    character(ESMF_MAXSTR) :: cvalue
    character(ESMF_MAXSTR) :: case_name      ! case name
    integer                :: ymd            ! CTSM current date (YYYYMMDD)
    integer                :: yr             ! CTSM current year
    integer                :: mon            ! CTSM current month
    integer                :: day            ! CTSM current day
    integer                :: tod            ! CTSM current time of day (sec)
    integer                :: ymd_sync       ! Sync date (YYYYMMDD)
    integer                :: yr_sync        ! Sync current year
    integer                :: mon_sync       ! Sync current month
    integer                :: day_sync       ! Sync current day
    integer                :: tod_sync       ! Sync current time of day (sec)
    integer                :: dtime          ! time step increment (sec)
    integer                :: nstep          ! time step index
    logical                :: rstwr          ! .true. ==> write restart file before returning
    logical                :: nlend          ! .true. ==> last time-step
    logical                :: dosend         ! true => send data back to driver
    logical                :: doalb          ! .true. ==> do albedo calculation on this time step
    real(r8)               :: nextsw_cday    ! calday from clock of next radiation computation
    real(r8)               :: caldayp1       ! ctsm calday plus dtime offset
    integer                :: lbnum          ! input to memory diagnostic
    integer                :: g,i            ! counters
    real(r8)               :: calday         ! calendar day for nstep
    real(r8)               :: declin         ! solar declination angle in radians for nstep
    real(r8)               :: declinp1       ! solar declination angle in radians for nstep+1
    real(r8)               :: eccf           ! earth orbit eccentricity factor
    type(bounds_type)      :: bounds         ! bounds
    character(len=32)      :: rdate          ! date char string for restart file names
    character(len=*),parameter  :: subname=trim(modName)//':[lnd_run] '
    character(*),parameter :: F02 = "('[lnd_comp_esmf] ',a, d26.19)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_ModelAdvance:start::',lbnum)
    endif
#endif

    !----------------------
    ! Obtain orbital values
    !----------------------

    !call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) eccen

    !call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) obliqr

    !call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) lambm0

    !call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) mvelpp

    !--------------------------------
    ! Get processor bounds
    !--------------------------------

    call get_proc_bounds(bounds)

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    call t_startf ('lc_lnd_import')
    call import_fields( gcomp, bounds, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('lc_lnd_import')

    !--------------------------------
    ! Run model
    !--------------------------------

    dtime = get_step_size()
    dosend = .false.
    do while(.not. dosend)

       ! TODO: This is currently hard-wired - is there a better way for nuopc?
       ! Note that the model clock is updated at the end of the time step not at the beginning
       nstep = get_nstep()
       if (nstep > 0) then
          dosend = .true.
       end if

       !--------------------------------
       ! Determine calendar day info
       !--------------------------------

       calday = get_curr_calday()
       caldayp1 = get_curr_calday(offset=dtime)

       !--------------------------------
       ! Get  time of next atmospheric shortwave calculation
       !--------------------------------

       ! TODO(NS): nextsw_cday should come directly from atmosphere!
       ! For now I am setting nextsw_cday to be the same caldayp1

       nextsw_cday = calday
       if (masterproc) then
          write(iulog,*) 'State_GetScalar ... nextsw_cday is : ', nextsw_cday
       end if

       !--------------------------------
       ! Determine doalb based on nextsw_cday sent from atm model
       !--------------------------------

       if (nstep == 0) then
          doalb = .false.
          nextsw_cday = caldayp1
       else if (nstep == 1) then
          !doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8)
          doalb = .false.
       else
          doalb = (nextsw_cday >= -0.5_r8)
       end if

       if (masterproc) then
          write(iulog,*) '------------  LILAC  ----------------'
          write(iulog,*) 'nstep       : ', nstep
          write(iulog,*) 'dtime       : ', dtime
          write(iulog,F02) 'calday      : ', calday
          write(iulog,F02) 'caldayp1    : ', caldayp1
          write(iulog,F02) 'nextsw_cday : ', nextsw_cday
          write(iulog,*) '-------------------------------------'
       end if

       call update_rad_dtime(doalb)

       if (masterproc) then
          write(iulog,*) 'doalb is: ', doalb
       end if

       !--------------------------------
       ! Determine if time to write restart
       !--------------------------------

       !call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
       !if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       !if (ChkErr(rc,__LINE__,u_FILE_u)) return
       !rstwr = .true.
       !call ESMF_AlarmRingerOff( alarm, rc=rc )
       !if (ChkErr(rc,__LINE__,u_FILE_u)) return
       !else
       !rstwr = .false.
       !endif

       ! TODO: for now hardwire rstwr to .false.
       rstwr = .false.

       !--------------------------------
       ! Determine if time to stop
       !--------------------------------

       !call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
       !if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       !if (ChkErr(rc,__LINE__,u_FILE_u)) return
       !nlend = .true.
       !call ESMF_AlarmRingerOff( alarm, rc=rc )
       !if (ChkErr(rc,__LINE__,u_FILE_u)) return
       !else
       !   nlend = .false.
       !endif

       !--------------------------------
       ! Run CTSM
       !--------------------------------

       call t_barrierf('sync_ctsm_run1', mpicom)

       call t_startf ('shr_orb_decl')
       calday = get_curr_calday()

       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )

       if (masterproc) then
          write(iulog,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          write(iulog,*) 'doalb           :  ', doalb
          write(iulog,*) 'call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, decl'
          write(iulog,*) 'call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0,  obliqr, decl'
          write(iulog,F02) 'calday is   :  ', calday
          write(iulog,F02) 'eccen is    :  ', eccen
          write(iulog,F02) 'mvelpp is   :  ', mvelpp
          write(iulog,F02) 'lambm0 is   :  ', lambm0
          write(iulog,F02) 'obliqr is   :  ', obliqr
          write(iulog,F02) 'clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate, rof_prognostic)'
          write(iulog,F02) 'declin is   :  ', declin
          write(iulog,F02) 'declinp1 is :  ', declinp1
          write(iulog,*  ) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       end if

       call t_stopf ('shr_orb_decl')

       ! Restart File - use nexttimestr rather than currtimestr here since that is the time at the end of
       ! the timestep and is preferred for restart file names
       ! TODO: is this correct for lilac?

       call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(nexttime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync, mon_sync, day_sync, tod_sync

       call t_startf ('ctsm_run')
       call clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate, rof_prognostic=.false.)
       call t_stopf ('ctsm_run')

       !--------------------------------
       ! Pack export state
       !--------------------------------

       call t_startf ('lc_lnd_export')
       call export_fields(gcomp, bounds,  rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf ('lc_lnd_export')

       !--------------------------------
       ! Advance ctsm time step
       !--------------------------------

       call t_startf ('lc_ctsm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_ctsm2_adv_timestep')

    end do

    ! Check that internal clock is in sync with master clock
    ! Note that the driver clock has not been updated yet - so at this point
    ! CTSM is actually 1 coupling intervals ahead of the driver clock

    call get_curr_date( yr, mon, day, tod, offset=-2*dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) ) then
       write(iulog,*)'ctsm ymd=',ymd     ,' ctsm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       rc = ESMF_FAILURE
       call ESMF_LogWrite(subname//" CTSM clock not in sync with Master Sync clock",ESMF_LOGMSG_ERROR)
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    !if (dbug > 1) then
    !   call State_diagnose(exportState,subname//':ES',rc=rc)
    !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !   if (masterproc) then
    !      call log_clock_advance(clock, 'CTSM', iulog, rc)
    !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !   end if
    !end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_ModelAdvance:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine lnd_run

  !---------------------------------------------------------------------------

  subroutine lnd_final(comp, import_state, export_state, clock, rc)
    !
    ! !DESCRIPTION:
    ! Finalize land surface model
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp            ! CLM gridded component
    type(ESMF_State)     :: import_state    ! CLM import state
    type(ESMF_State)     :: export_state    ! CLM export state
    type(ESMF_Clock)     :: clock          ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Destroy ESMF objects
    !call esmfshr_util_StateArrayDestroy(export_state,'domain',rc)
    !if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !call esmfshr_util_StateArrayDestroy(export_state,'d2x',rc)
    !if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !call esmfshr_util_StateArrayDestroy(import_state,'x2d',rc)
    !if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine lnd_final

  !===============================================================================

  subroutine log_clock_advance(clock, logunit, rc)

    ! input/output variables
    type(ESMF_Clock)               :: clock
    integer          , intent(in)  :: logunit
    integer          , intent(out) :: rc

    ! local variables
    character(len=CL) :: cvalue, prestring
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    write(prestring, *) "------>Advancing CTSM from: "
    call ESMF_ClockPrint(clock, options="currTime", unit=cvalue, preString=trim(prestring), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

    call ESMF_ClockPrint(clock, options="stopTime", unit=cvalue, &
         preString="--------------------------------> to: ", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

  end subroutine log_clock_advance

!===============================================================================

  subroutine memcheck(string, level, mastertask)

    ! input/output variables
    character(len=*) , intent(in) :: string
    integer          , intent(in) :: level
    logical          , intent(in) :: mastertask

    ! local variables
    integer :: ierr
    integer, external :: GPTLprint_memusage
    !-----------------------------------------------------------------------

    if ((mastertask .and. memdebug_level > level) .or. memdebug_level > level+1) then
       ierr = GPTLprint_memusage(string)
    endif

  end subroutine memcheck

end module lnd_comp_esmf
