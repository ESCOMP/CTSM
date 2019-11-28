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
  use spmdMod           , only : masterproc, spmd_init, mpicom
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
  use glc_elevclass_mod , only : glc_elevclass_init  ! TODO: is this needed?

  implicit none
  private                         ! By default make data private except

  public :: lnd_register          ! register clm initial, run, final methods
  public :: lnd_init              ! clm initialization
  public :: lnd_run               ! clm run phase
  public :: lnd_final             ! clm finalization/cleanup

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer                 :: glc_nec = 10   ! fixed # of glc elevation classes
  integer,      parameter :: memdebug_level=1
  integer,      parameter :: dbug = 1
  character(*), parameter :: modName =  "lnd_comp_esmf"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lnd_register(comp, rc)

    ! Register the clm initial, run, and final phase methods with ESMF.

    ! input/output argumenents
    type(ESMF_GridComp)  :: comp  ! CLM grid component
    integer, intent(out) :: rc    ! return status
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogSet ( flush =.true.)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, lnd_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, lnd_run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, lnd_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite("lnd gridcompset entry points finished!", ESMF_LOGMSG_INFO)

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
    integer                        :: ierr                  ! error code
    integer                        :: n,g,i,j               ! indices
    logical                        :: exists                ! true if file exists
    real(r8)                       :: nextsw_cday           ! calday from clock of next radiation computation
    character(len=CL)              :: caseid                ! case identifier name
    character(len=CL)              :: ctitle                ! case description title
    character(len=CL)              :: starttype             ! start-type (startup, continue, branch, hybrid)
    integer                        :: nsrest                ! clm restart type
    logical                        :: brnch_retain_casename ! flag if should retain the case name on a branch start type
    logical                        :: atm_aero              ! Flag if aerosol data sent from atm model
    integer                        :: lbnum                 ! input to memory diagnostic
    integer                        :: shrlogunit            ! old values for log unit and log level
    type(bounds_type)              :: bounds                ! bounds

    ! generation of field bundles
    type(ESMF_State)               :: importState, exportState
    type(ESMF_FieldBundle)         :: c2l_fb
    type(ESMF_FieldBundle)         :: l2c_fb

    ! mesh generation
    type(ESMF_Mesh)                :: lnd_mesh
    character(ESMF_MAXSTR)         :: lnd_mesh_filename     ! full filepath of land mesh file
    integer                        :: nlnd, nocn            ! local size ofarrays
    integer, pointer               :: gindex(:)             ! global index space for land and ocean points
    integer, pointer               :: gindex_lnd(:)         ! global index space for just land points
    integer, pointer               :: gindex_ocn(:)         ! global index space for just ocean points
    type(ESMF_DistGrid)            :: distgrid

    ! clock info
    character(len=CL)              :: calendar              ! calendar type name
    type(ESMF_CalKind_Flag)        :: caltype               ! calendar type from lilac clock
    integer                        :: curr_tod, curr_ymd    ! current time info
    integer                        :: yy, mm, dd            ! query output from lilac clock
    integer                        :: dtime_lilac           ! coupling time-step from the input lilac clock
    integer                        :: ref_ymd               ! reference date (YYYYMMDD)
    integer                        :: ref_tod               ! reference time of day (sec)
    integer                        :: start_ymd             ! start date (YYYYMMDD)
    integer                        :: start_tod             ! start time of day (sec)
    integer                        :: stop_ymd              ! stop date (YYYYMMDD)
    integer                        :: stop_tod              ! stop time of day (sec)
    type(ESMF_Time)                :: currTime              ! Current time
    type(ESMF_Time)                :: startTime             ! Start time
    type(ESMF_Time)                :: stopTime              ! Stop time
    type(ESMF_Time)                :: refTime               ! Ref time
    type(ESMF_TimeInterval)        :: timeStep              ! time step from lilac clock

    ! orbital info
    integer                        :: orb_iyear_align       ! associated with model year
    integer                        :: orb_cyear             ! orbital year for current orbital computation
    integer                        :: orb_iyear             ! orbital year for current orbital computation
    integer                        :: orb_eccen             ! orbital year for current orbital computation

    ! for pio_init2 and mct
    type(ESMF_VM)               :: vm
    integer                     :: mpicom_vm
    integer                     :: ncomps = 1
    integer, pointer            :: mycomms(:)                 ! for mct
    integer, pointer            :: myids(:)                   ! for mct
    integer                     :: compids(1) = (/1/)         ! for both mct and pio_init2 - array with component ids
    integer                     :: comms(1)                   ! for both mct and pio_init2 - array with mpicoms
    character(len=32)           :: compLabels(1) = (/'LND'/)  ! for pio_init2
    character(len=64)           :: comp_name(1) = (/'LND'/)   ! for pio_init2
    logical                     :: comp_iamin(1) = (/.true./) ! for pio init2
    integer                     :: iam(1)                     ! for pio_init2
    character(len=*), parameter :: subname=trim(modName)//': (lnd_init) '
    !------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' is called!', ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    ! Query VM for local PET and mpi communicator
    !------------------------------------------------------------------------

    ! NOTE : both MPI_INIT and PIO_INIT1 are initialized in lilac_mod.F90

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, localPet=iam(1), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    call ESMF_LogWrite(subname//"ESMF_VMGet", ESMF_LOGMSG_INFO)

    !call ESMF_VMPrint (vm, rc = rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    comms(1) = mpicom_vm

    !------------------------------------------------------------------------
    ! Initialize pio_init2 TODO: is this needed here?
    !------------------------------------------------------------------------

    call shr_pio_init2(compids, compLabels, comp_iamin, comms, iam)
    call ESMF_LogWrite(subname//"initialized shr_pio_init2 ...", ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    ! Initialize mct - needed for data model share code - e.g. nitrogen deposition
    !------------------------------------------------------------------------

    allocate(mycomms(1), myids(1))
    mycomms = (/mpicom_vm/) ; myids = (/1/)

    call mct_world_init(ncomps, mpicom_vm, mycomms, myids)
    call ESMF_LogWrite(subname//"initialized mct ...  ", ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    ! Initialize internal ctsm MPI info
    !------------------------------------------------------------------------

    call spmd_init( clm_mpicom=mpicom_vm, lndid=1)
    call ESMF_LogWrite(subname//"initialized model mpi info using spmd_init", ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    !--- Log File ---
    !------------------------------------------------------------------------

    ! TODO: by default iulog = 6 in clm_varctl - this should be generalized so that we
    ! can control the output log file for ctsm running with a lilac driver

    inst_name  = 'LND'; inst_index  = 1; inst_suffix = ""

    ! Initialize io log unit
    call shr_file_getLogUnit (shrlogunit)
    if (.not. masterproc) then
       iulog = shrlogunit  ! All shr code output will go to iulog for masterproc
    end if
    call shr_file_setLogUnit (iulog)

    if (masterproc) then
       write(iulog,*) "========================================="
       write(iulog,*) " starting (lnd_comp_esmf): lnd_comp_init "
       write(iulog,*) " CLM land model initialization"
    end if

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

    if (masterproc) then
       write(iulog,*) 'shr_obs_params is setting the following:'
       write(iulog,*) 'eccen is  : ', eccen
       write(iulog,*) 'mvelpp is : ', mvelpp
       write(iulog,*) 'lambm0 is : ', lambm0
       write(iulog,*) 'obliqr is : ', obliqr
    end if

    !----------------------
    ! Consistency check on namelist filename
    !----------------------
    call control_setNL("lnd_in")

    ! TODO: how do we set case_name and nsrest - should we hardwire for now?
    caseid = 'test_lilac'
    nsrest = nsrStartup

    !----------------------
    ! Initialize module variables in clm_time_manger.F90
    !----------------------

    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, timeStep=timeStep, rc=rc)
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

    call ESMF_TimeGet( currTime, calkindflag=caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    ! The following sets the module variables in clm_time_mamanger.F90 - BUT DOES NOT intialize the
    ! clock. Routine timemgr_init (called by initialize1) initializes the clock using the module variables
    ! that have been set via calls to set_timemgr_init.

    call set_timemgr_init( &
         calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd, stop_tod_in=stop_tod)

    !----------------------
    ! Read namelist, grid and surface data
    !----------------------

    ! set default values for run control variables
    call clm_varctl_set(caseid_in=caseid, nsrest_in=nsrest)
    call ESMF_LogWrite(subname//"default values for run control variables are set...", ESMF_LOGMSG_INFO)

    !----------------------
    ! Initialize glc_elevclass module
    !----------------------

    call glc_elevclass_init(glc_nec)  ! TODO: is this needed still?

    !----------------------
    ! Call initialize1
    !----------------------

    call ESMF_TimeIntervalGet(timeStep, s=dtime_lilac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (masterproc) then
       write(iulog,*)'dtime_lilac= ',dtime_lilac
    end if

    ! Note that routine controlMod.F90 will initialze the dtime module
    ! variable in clm_time_manager to the dtime_lilac AND NOT the
    ! dtime read in from the clm_inparm namelist in this case.  Note
    ! that the memory for gindex_ocn will be allocated in the following call

    call initialize1(gindex_ocn=gindex_ocn, dtime_driver=dtime_lilac)

    call ESMF_LogWrite(subname//"ctsm time manager initialized....", ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(subname//"ctsm initialize1 done...", ESMF_LOGMSG_INFO)

    !--------------------------------
    ! generate the land mesh on ctsm distribution
    !--------------------------------

    ! TODO: mesh file should come into clm as a namelist for lilac only
    ! for now need to hardwire this in lnd_mesh_filename here

    lnd_mesh_filename = '/glade/p/cesmdata/cseg/inputdata/share/meshes/fv4x5_050615_polemod_ESMFmesh.nc'

    ! obtain global index array for just land points which includes mask=0 or ocean points
    call get_proc_bounds( bounds )

    nlnd = bounds%endg - bounds%begg + 1
    allocate(gindex_lnd(nlnd))
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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    deallocate(gindex)
    call ESMF_LogWrite(subname//"DistGrid created......", ESMF_LOGMSG_INFO)

    lnd_mesh = ESMF_MeshCreate(filename=trim(lnd_mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, elementDistgrid=Distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (masterproc) then
       write(iulog,*)'mesh file for domain is ',trim(lnd_mesh_filename)
    end if
    call ESMF_LogWrite(subname//" Create Mesh using file ...."//trim(lnd_mesh_filename), ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Finish initializing ctsm
    !--------------------------------

    call initialize2()
    call ESMF_LogWrite(subname//"ctsm initialize2 done...", ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Create import state (only assume input from atm - not rof and glc)
    !--------------------------------

    ! First create an empty field bundle
    c2l_fb =  ESMF_FieldBundleCreate ( name='c2l_fb', rc=rc)

    ! Now add fields on lnd_mesh to this field bundle
    call fldbundle_add( 'Sa_z'       , c2l_fb,rc) !1
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sa_topo'    , c2l_fb,rc) !2
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sa_u'       , c2l_fb,rc) !3
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sa_v'       , c2l_fb,rc) !4
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sa_ptem'    , c2l_fb,rc) !5
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sa_pbot'    , c2l_fb,rc) !6
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sa_tbot'    , c2l_fb,rc) !7
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sa_shum'    , c2l_fb,rc) !8
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_lwdn'  , c2l_fb,rc) !9
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_rainc' , c2l_fb,rc) !10
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_rainl' , c2l_fb,rc) !11
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_snowc' , c2l_fb,rc) !12
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_snowl' , c2l_fb,rc) !13
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_swndr' , c2l_fb,rc) !14
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_swvdr' , c2l_fb,rc) !15
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_swndf' , c2l_fb,rc) !16
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Faxa_swvdf' , c2l_fb,rc) !17
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateAdd(import_state, fieldbundleList = (/c2l_fb/), rc=rc)

    !--------------------------------
    ! Create export state
    !--------------------------------

    ! First create an empty field bundle
    l2c_fb = ESMF_FieldBundleCreate(name='l2c_fb', rc=rc)

    ! Now add fields on lnd_mesh to this field bundle
    call fldbundle_add( 'Sl_lfrin'    ,  l2c_fb,rc)   !1
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sl_t'        ,  l2c_fb,rc)   !2
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sl_tref'     ,  l2c_fb,rc)   !3
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sl_qref'     ,  l2c_fb,rc)   !4
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sl_avsdr'    ,  l2c_fb,rc)   !5
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sl_anidr'    ,  l2c_fb,rc)   !6
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sl_avsdf'    ,  l2c_fb,rc)   !7
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sl_anidf'    ,  l2c_fb,rc)   !8
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Sl_snowh'    ,  l2c_fb,rc)   !9
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Fall_u10'    ,  l2c_fb,rc)   !10
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Fall_fv'     ,  l2c_fb,rc)   !11
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add( 'Fall_ram1'   ,  l2c_fb,rc)   !12
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call fldbundle_add( 'Fall_taux'   ,  l2c_fb,rc)   !10
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call fldbundle_add( 'Fall_lwup'   ,  l2c_fb,rc)   !14
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call fldbundle_add( 'Fall_evap'   ,  l2c_fb,rc)   !15
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call fldbundle_add( 'Fall_swnet'  ,  l2c_fb,rc)   !16
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateAdd(export_state, fieldbundleList = (/l2c_fb/), rc=rc)

    !--------------------------------
    ! Create land export state
    !--------------------------------

    call ESMF_LogWrite(subname//"Creating land export state", ESMF_LOGMSG_INFO)

    ! Fill in export state at end of initialization
    call export_fields(comp, bounds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//"Getting Calendar Day of nextsw calculation...", ESMF_LOGMSG_INFO)

    ! Get calendar day of next sw (shortwave) calculation (nextsw_cday)
    if (nsrest == nsrStartup) then
       call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
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

    write(iulog,*) " finished (lnd_comp_esmf): lnd_comp_init "
    write(iulog,*) "========================================="

  !---------------------------
  contains
  !---------------------------

    subroutine fldbundle_add(stdname, fieldbundle,rc)
      !---------------------------
      ! Create an empty input field with name 'stdname' to add to fieldbundle
      !---------------------------

      ! input/output variables
      character(len=*)        , intent(in)    :: stdname
      type (ESMF_FieldBundle) , intent(inout) :: fieldbundle
      integer                 , intent(out)   :: rc
      ! local variables
      type(ESMF_Field) :: field
      !-------------------------------------------------------------------------------
      rc = ESMF_SUCCESS
      field = ESMF_FieldCreate(lnd_mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT , name=trim(stdname), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldBundleAdd(fieldbundle, (/field/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end subroutine fldbundle_add

  end subroutine lnd_init

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
    character(*)    , parameter :: F02 = "('[lnd_comp_esmf] ',a, d26.19)"
    character(len=*), parameter :: subname=trim(modName)//':[lnd_run] '
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
          write(iulog,*) trim(subname) // '... nextsw_cday is : ', nextsw_cday
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
          write(iulog,*) trim(subname) // 'nstep       : ', nstep
          write(iulog,*) trim(subname) // 'dtime       : ', dtime
          write(iulog,*) trim(subname) // 'calday      : ', calday
          write(iulog,*) trim(subname) // 'caldayp1    : ', caldayp1
          write(iulog,*) trim(subname) // 'nextsw_cday : ', nextsw_cday
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
          write(iulog,*)   'doalb       :  ', doalb
          write(iulog,F02) 'calday is   :  ', calday
          write(iulog,F02) 'eccen is    :  ', eccen
          write(iulog,F02) 'mvelpp is   :  ', mvelpp
          write(iulog,F02) 'lambm0 is   :  ', lambm0
          write(iulog,F02) 'obliqr is   :  ', obliqr
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

       call advance_timestep()

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
    !---------------------------------
    ! Finalize land surface model
    !---------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: comp            ! CLM gridded component
    type(ESMF_State)     :: import_state    ! CLM import state
    type(ESMF_State)     :: export_state    ! CLM export state
    type(ESMF_Clock)     :: clock          ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! TODO: Destroy ESMF objects

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
