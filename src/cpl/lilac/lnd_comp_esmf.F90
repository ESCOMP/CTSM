module lnd_comp_esmf

  !----------------------------------------------------------------------------
  ! This is the ESMF cap for CTSM
  ! NOTE : both mpi_init and pio_init1 are initialized in lilac_mod.F90
  !----------------------------------------------------------------------------

  ! external libraries
  use ESMF              , only : ESMF_GridComp, ESMF_SUCCESS, ESMF_LogSet, ESMF_State
  use ESMF              , only : ESMF_Clock, ESMF_FieldBundle, ESMF_MAXSTR, ESMF_Field
  use ESMF              , only : ESMF_FieldCreate, ESMF_AttributeGet
  use ESMF              , only : ESMF_Time, ESMF_LogWrite, ESMF_LogFoundError, ESMF_Finalize
  use ESMF              , only : ESMF_FieldBundleAdd, ESMF_FieldBundleCreate
  use ESMF              , only : ESMF_ClockGet, ESMF_ClockGetAlarm, ESMF_LOGMSG_INFO
  use ESMF              , only : ESMF_TYPEKIND_R8, ESMF_MESHLOC_ELEMENT, ESMF_LOGERR_PASSTHRU
  use ESMF              , only : ESMF_END_ABORT, ESMF_TimeGet, ESMF_LOGMSG_ERROR
  use shr_mpi_mod       , only : shr_mpi_bcast
  use perf_mod          , only : t_startf, t_stopf, t_barrierf

  ! lilac code
  use ctsm_LilacCouplingFields, only : a2l_fields, l2a_fields

  ! cime share code
  use shr_kind_mod      , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod       , only : shr_sys_abort
  use shr_file_mod      , only : shr_file_setLogUnit, shr_file_getLogUnit
  use shr_orb_mod       , only : shr_orb_decl, shr_orb_params
  use shr_cal_mod       , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use shr_nl_mod        , only : shr_nl_find_group_name
  use glc_elevclass_mod , only : glc_elevclass_init

  ! ctsm code
  use spmdMod           , only : masterproc, spmd_init, mpicom
  use decompMod         , only : bounds_type, get_proc_bounds
  use domainMod         , only : ldomain
  use controlMod        , only : control_setNL
  use clm_varorb        , only : eccen, obliqr, lambm0, mvelpp
  use clm_varctl        , only : clm_varctl_set, iulog, finidat
  use clm_varctl        , only : nsrStartup, nsrContinue
  use clm_varctl        , only : inst_index, inst_suffix, inst_name
  use clm_time_manager  , only : set_timemgr_init, advance_timestep
  use clm_time_manager  , only : update_rad_dtime
  use clm_time_manager  , only : get_nstep, get_step_size
  use clm_time_manager  , only : get_curr_date, get_curr_calday
  use clm_initializeMod , only : initialize1, initialize2
  use clm_driver        , only : clm_drv
  use lnd_import_export , only : import_fields, export_fields
  use lnd_shr_methods   , only : chkerr, state_diagnose
  use lnd_comp_shr      , only : mesh, model_meshfile, model_clock
  use lnd_set_decomp_and_domain, only :lnd_set_decomp_and_domain_from_readmesh

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
    use ESMF     , only : ESMF_GridCompSetEntryPoint
    use ESMF     , only : ESMF_METHOD_INITIALIZE, ESMF_METHOD_RUN, ESMF_METHOD_FINALIZE

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
    ! Uses:
    use ESMF              , only : ESMF_VM, ESMF_VMGet, ESMF_VMGetCurrent
    use ESMF              , only : ESMF_DistGrid, ESMF_AttributeSet
    use ESMF              , only : ESMF_CalKind_Flag, ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
    use ESMF              , only : ESMF_TimeInterval, ESMF_TimeIntervalGet
    use ESMF              , only : ESMF_StateAdd
    use ESMF              , only : operator(==)

    use shr_dust_emis_mod , only : shr_dust_emis_readnl

    ! input/output variables
    type(ESMF_GridComp)  :: comp         ! CLM gridded component
    type(ESMF_State)     :: import_state ! CLM import state
    type(ESMF_State)     :: export_state ! CLM export state
    type(ESMF_Clock)     :: clock        ! ESMF synchronization clock
    integer, intent(out) :: rc           ! Return code

    ! local variables
    integer                    :: ierr                       ! error code
    integer                    :: n,g,i,j                    ! indices
    logical                    :: exists                     ! true if file exists
    character(len=CL)          :: caseid                     ! case identifier name
    character(len=CL)          :: starttype                  ! start-type (startup, continue, branch, hybrid)
    integer                    :: nsrest                     ! clm restart type
    integer                    :: lbnum                      ! input to memory diagnostic
    integer                    :: shrlogunit                 ! old values for log unit and log level
    type(bounds_type)          :: bounds                     ! bounds
    character(len=CL)          :: cvalue

    ! communicator info
    type(ESMF_VM)              :: vm
    integer                    :: mpicom_vm

    ! generation of field bundles
    type(ESMF_State)           :: importState, exportState
    type(ESMF_FieldBundle)     :: c2l_fb_atm, c2l_fb_rof     ! field bundles in import state
    type(ESMF_FieldBundle)     :: l2c_fb_atm, l2c_fb_rof     ! field bundles in export state
    type(ESMF_Field)           :: lfield

    ! mesh generation
    character(ESMF_MAXSTR)     :: lnd_mesh_filename          ! full filepath of land mesh file
    integer, pointer           :: gindex(:)                  ! global index space for land and ocean points
    type(ESMF_DistGrid)        :: distgrid
    integer                    :: fileunit
    integer                    :: ni, nj

    ! clock info
    character(len=CL)          :: calendar                   ! calendar type name
    type(ESMF_CalKind_Flag)    :: caltype                    ! calendar type from lilac clock
    integer                    :: curr_tod, curr_ymd         ! current time info
    integer                    :: yy, mm, dd                 ! query output from lilac clock
    integer                    :: dtime_lilac                ! coupling time-step from the input lilac clock
    integer                    :: ref_ymd                    ! reference date (YYYYMMDD)
    integer                    :: ref_tod                    ! reference time of day (sec)
    integer                    :: start_ymd                  ! start date (YYYYMMDD)
    integer                    :: start_tod                  ! start time of day (sec)
    type(ESMF_Time)            :: currTime                   ! Current time
    type(ESMF_Time)            :: startTime                  ! Start time
    type(ESMF_Time)            :: refTime                    ! Ref time
    type(ESMF_TimeInterval)    :: timeStep                   ! time step from lilac clock

    ! orbital info
    integer                    :: orb_iyear_align            ! associated with model year
    integer                    :: orb_cyear                  ! orbital year for current orbital computation
    integer                    :: orb_iyear                  ! orbital year for current orbital computation
    integer                    :: orb_eccen                  ! orbital year for current orbital computation

    character(len=*), parameter :: subname=trim(modName)//': (lnd_init) '
    !------------------------------------------------------------------------

    ! input namelist read for ctsm mesh
    namelist /lilac_lnd_input/ lnd_mesh_filename

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' is called!', ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    ! Query VM for local PET and mpi communicator
    !------------------------------------------------------------------------

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    call ESMF_LogWrite(subname//"ESMF_VMGet", ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    ! Initialize internal ctsm MPI info
    !------------------------------------------------------------------------

    call spmd_init( clm_mpicom=mpicom_vm, lndid=1)
    call ESMF_LogWrite(subname//"initialized model mpi info using spmd_init", ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    ! Initialize output log file
    !------------------------------------------------------------------------

    ! TODO: by default iulog = 6 in clm_varctl - this should be generalized so that we
    ! can control the output log file for ctsm running with a lilac driver
    !
    ! See also https://github.com/ESCOMP/CTSM/issues/861

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
    !
    ! See also https://github.com/ESCOMP/CTSM/issues/865

    !call shr_cal_date2ymd(ymd,year,month,day)
    !orb_cyear = orb_iyear + (year - orb_iyear_align)

    orb_cyear = 2000
    call shr_orb_params(orb_cyear, eccen, obliqr, mvelpp, obliqr, lambm0, mvelpp, masterproc)

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

    !----------------------
    ! Read in lilac_in namelists
    !----------------------

    if (masterproc) then
       open(newunit=fileunit, status="old", file="lilac_in")
       call shr_nl_find_group_name(fileunit, 'lilac_lnd_input', ierr)
       if (ierr == 0) then
          read(fileunit, lilac_lnd_input, iostat=ierr)
          if (ierr > 0) then
             call shr_sys_abort( 'problem on read of lilac_lnd_input')
          end if
       end if
       close(fileunit)
    end if
    call shr_mpi_bcast(lnd_mesh_filename, mpicom)

    ! Fill in the value for model_meshfile in lnd_comp_shr used by the stream routines in share_esmf/
    model_meshfile = trim(lnd_mesh_filename)

    ! Reading in the drv_flds_in namelist is required for dust emissions
    call shr_dust_emis_readnl( mpicom, "drv_flds_in")

    !----------------------
    ! Obtain caseid and start type from attributes in import state
    !----------------------

    call ESMF_AttributeGet(import_state, name="caseid", value=caseid, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_LogWrite(subname//"caseid is "//trim(caseid), ESMF_LOGMSG_INFO)

    call ESMF_AttributeGet(import_state, name="starttype", value=starttype, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_LogWrite(subname//"starttype is "//trim(starttype), ESMF_LOGMSG_INFO)

    if (trim(starttype) == trim('startup')) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim('continue') ) then
       nsrest = nsrContinue
    else
       call shr_sys_abort( subname//' ERROR: unknown starttype'//trim(starttype) )
    end if

    !----------------------
    ! Initialize module variables in clm_time_manger.F90
    !----------------------

    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, refTime=RefTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

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

    call ESMF_TimeIntervalGet(timeStep, s=dtime_lilac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (masterproc) then
       write(iulog,*)'dtime = ',dtime_lilac
    end if

    ! The following sets the module variables in clm_time_mamanger.F90 - BUT DOES NOT intialize the
    ! clock. Routine timemgr_init (called by initialize1) initializes the clock using the module variables
    ! that have been set via calls to set_timemgr_init.

    ! Note that we assume that CTSM's internal dtime matches the coupling time step.
    ! i.e., we currently do NOT allow sub-cycling within a coupling time step.
    call set_timemgr_init( &
         calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, dtime_in=dtime_lilac)
    call ESMF_LogWrite(subname//"ctsm time manager initialized....", ESMF_LOGMSG_INFO)

    ! Set model clock in lnd_comp_shr
    model_clock = clock

    !----------------------
    ! Read namelist, grid and surface data
    !----------------------
    ! set default values for run control variables
    call clm_varctl_set(caseid_in=caseid, nsrest_in=nsrest)
    ! Initialize glc_elevclass module
    call glc_elevclass_init(glc_nec)
    call ESMF_LogWrite(subname//"default values for run control variables are set...", ESMF_LOGMSG_INFO)

    !----------------------
    ! Call initialize1
    !----------------------
    call initialize1(dtime=dtime_lilac)
    call ESMF_LogWrite(subname//"ctsm initialize1 done...", ESMF_LOGMSG_INFO)

    !----------------------
    ! Initialize decomposition and domain (ldomain) types and generate land mesh
    !----------------------
    ! TODO: generalize this so that a mask mesh is read in like for nuopc/cmeps
    ! For now set the meshfile_mask equal to the model_meshfile
    call lnd_set_decomp_and_domain_from_readmesh(driver='lilac', vm=vm, &
         meshfile_lnd=lnd_mesh_filename, meshfile_mask=lnd_mesh_filename, &
         mesh_ctsm=mesh, ni=ni, nj=nj, rc=rc)

    !--------------------------------
    ! Finish initializing ctsm
    !--------------------------------
    call initialize2(ni,nj, currTime)
    call ESMF_LogWrite(subname//"ctsm initialize2 done...", ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Create import state (only assume input from atm - not rof and glc)
    !--------------------------------
    ! create an empty field bundle for import of atm fields
    c2l_fb_atm = ESMF_FieldBundleCreate (name='c2l_fb_atm', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! now add atm import fields on mesh to this field bundle
    do n = 1, a2l_fields%num_fields()
       lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT, &
            name=a2l_fields%get_fieldname(n), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(c2l_fb_atm, (/lfield/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! add the field bundle to the state
    call ESMF_StateAdd(import_state, fieldbundleList = (/c2l_fb_atm/))
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create an empty field bundle for the import of rof fields
    c2l_fb_rof = ESMF_FieldBundleCreate (name='c2l_fb_rof', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldbundle_add('Flrr_flood', c2l_fb_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrr_volr', c2l_fb_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrr_volrmch', c2l_fb_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add the field bundle to the state
    call ESMF_StateAdd(import_state, fieldbundleList = (/c2l_fb_rof/))
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Create ctsm export state
    !--------------------------------

    ! create an empty field bundle for atm export fields
    l2c_fb_atm = ESMF_FieldBundleCreate(name='l2c_fb_atm', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! now add atm export fields on mesh to this field bundle
    do n = 1, l2a_fields%num_fields()
       lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT, &
            name=l2a_fields%get_fieldname(n), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(l2c_fb_atm, (/lfield/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! add the field bundle to the state
    call ESMF_StateAdd(export_state, fieldbundleList = (/l2c_fb_atm/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create an empty field bundle for rof export fields
    l2c_fb_rof = ESMF_FieldBundleCreate(name='l2c_fb_rof', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! now add rof export fields on mesh to this field bundle
    call fldbundle_add('Flrl_rofsur', l2c_fb_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrl_rofgwl', l2c_fb_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrl_rofsub', l2c_fb_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrl_rofi', l2c_fb_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrl_irrig', l2c_fb_rof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateAdd(export_state, fieldbundleList = (/l2c_fb_rof/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Fill in ctsm export state
    !--------------------------------

    call get_proc_bounds( bounds )

    call export_fields(export_state, bounds, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(cvalue,*) ldomain%ni
    call ESMF_AttributeSet(export_state, name="lnd_nx", value=trim(cvalue), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_LogWrite(subname//"set attribute lnd_nx to "//trim(cvalue), ESMF_LOGMSG_INFO)

    write(cvalue,*) ldomain%nj
    call ESMF_AttributeSet(export_state, name="lnd_ny", value=trim(cvalue), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_LogWrite(subname//"set attribute lnd_ny to "//trim(cvalue), ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//"Created land export state", ESMF_LOGMSG_INFO)

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

    if (masterproc) then
       write(iulog,*) " finished (lnd_comp_esmf): lnd_comp_init "
       write(iulog,*) "========================================="
    end if

  !---------------------------
  contains
  !---------------------------

    subroutine fldbundle_add(stdname, fieldbundle, rc)
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
      field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT , name=trim(stdname), rc=rc)
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
    use ESMF       , only : ESMF_Alarm, ESMF_AlarmIsRinging, ESMF_AlarmRingerOff
    use ESMF       , only : ESMF_FAILURE, ESMF_ClockGetNextTime

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
    character(ESMF_MAXSTR) :: cvalue
    integer                :: ymd            ! CTSM current date (YYYYMMDD)
    integer                :: yr             ! CTSM current year
    integer                :: mon            ! CTSM current month
    integer                :: day            ! CTSM current day
    integer                :: tod            ! CTSM current time of day (sec)
    integer                :: ymd_lilac       ! Sync date (YYYYMMDD)
    integer                :: yr_lilac        ! Sync current year
    integer                :: mon_lilac       ! Sync current month
    integer                :: day_lilac       ! Sync current day
    integer                :: tod_lilac       ! Sync current time of day (sec)
    integer                :: dtime          ! time step increment (sec)
    integer                :: nstep          ! time step index
    logical                :: rstwr          ! .true. ==> write restart file before returning
    logical                :: nlend          ! .true. ==> last time-step
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
    logical                :: first_call = .true.  ! true if and only if this is the first time this routine is called in this execution
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

    !call ESMF_AttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) eccen

    !call ESMF_AttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) obliqr

    !call ESMF_AttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) lambm0

    !call ESMF_AttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
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
    call import_fields(import_state, bounds, first_call, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('lc_lnd_import')

    !--------------------------------
    ! Run model
    !--------------------------------

    dtime = get_step_size()

    ! We assume that the land model time step matches the coupling interval.
    nstep = get_nstep()

    !--------------------------------
    ! Determine calendar day info
    !--------------------------------

    calday = get_curr_calday(reuse_day_365_for_day_366=.true.)
    caldayp1 = get_curr_calday(offset=dtime, reuse_day_365_for_day_366=.true.)

    !--------------------------------
    ! Get time of next atmospheric shortwave calculation
    !--------------------------------

    ! TODO(NS): nextsw_cday should come directly from atmosphere!
    ! For now I am setting nextsw_cday to be the same caldayp1
    !
    ! See also https://github.com/ESCOMP/CTSM/issues/860

    nextsw_cday = calday
    if (masterproc) then
       write(iulog,*) trim(subname) // '... nextsw_cday is : ', nextsw_cday
    end if

    !--------------------------------
    ! Obtain orbital values
    !--------------------------------

    call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
    call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )

    if (masterproc) then
       write(iulog,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       write(iulog,F02) 'nextsw_cday is :  ', nextsw_cday
       write(iulog,F02) 'calday is      :  ', calday
       write(iulog,F02) 'eccen is       :  ', eccen
       write(iulog,F02) 'mvelpp is      :  ', mvelpp
       write(iulog,F02) 'lambm0 is      :  ', lambm0
       write(iulog,F02) 'obliqr is      :  ', obliqr
       write(iulog,F02) 'declin is      :  ', declin
       write(iulog,F02) 'declinp1 is    :  ', declinp1
       write(iulog,*  ) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    end if

    !--------------------------------
    ! Determine doalb based on nextsw_cday sent from atm model
    !--------------------------------

    if (nstep == 1) then
       !doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8)
       doalb = .false.
    else
       doalb = (nextsw_cday >= -0.5_r8)
    end if

    if (masterproc) then
       write(iulog,*) '------------  LILAC  ----------------'
       write(iulog,*) 'nstep       : ', nstep
       write(iulog,*) 'calday      : ', calday
       write(iulog,*) 'caldayp1    : ', caldayp1
       write(iulog,*) 'nextsw_cday : ', nextsw_cday
       write(iulog,*) 'doalb       : ', doalb
       write(iulog,*) '-------------------------------------'
    end if

    call update_rad_dtime(doalb)

    !--------------------------------
    ! Determine if time to stop
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='lilac_stop_alarm', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       nlend = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       nlend = .false.
    endif
    if (masterproc) then
       write(iulog,*)' stop alarm is ',nlend
    end if

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='lilac_restart_alarm', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       rstwr = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       rstwr = .false.
    endif
    if (masterproc) then
       write(iulog,*)' restart alarm is ',rstwr
    end if

    !--------------------------------
    ! Run CTSM
    !--------------------------------

    call t_barrierf('sync_ctsm_run1', mpicom)

    ! Restart File - use nexttimestr rather than currtimestr here since that is the time at the end of
    ! the timestep and is preferred for restart file names
    ! TODO: is this correct for lilac?

    call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nexttime, yy=yr_lilac, mm=mon_lilac, dd=day_lilac, s=tod_lilac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_lilac, mon_lilac, day_lilac, tod_lilac

    call t_startf ('ctsm_run')
    call clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate, rof_prognostic=.false.)
    call t_stopf ('ctsm_run')

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call t_startf ('lc_lnd_export')
    call export_fields(export_state, bounds,  rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('lc_lnd_export')

    !--------------------------------
    ! Advance ctsm time step
    !--------------------------------

    call advance_timestep()

    !--------------------------------
    ! Check that internal clock is in sync with lilac driver clock
    !--------------------------------

    ! Get ctsm current time info
    call get_curr_date( yr, mon, day, tod, offset=-2*dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod

    ! Get lilac clock info
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yr_lilac, mm=mon_lilac, dd=day_lilac, s=tod_lilac, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_lilac, mon_lilac, day_lilac, ymd_lilac)

    if (masterproc) then
       write(iulog,*)'lilac ymd=',ymd     ,' lilac tod= ',tod
    end if

    ! Note that the driver clock has not been updated yet - so at this point
    ! CTSM is actually 1 coupling intervals ahead of the driver clock

    if ( (ymd /= ymd_lilac) .or. (tod /= tod_lilac) ) then
       write(iulog,*)'ctsm  ymd=',ymd      ,' ctsm  tod= ',tod
       write(iulog,*)'lilac ymd=',ymd_lilac,' lilac tod= ',tod_lilac
       call ESMF_LogWrite(subname//" CTSM clock not in sync with lilac clock",ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_comp_nuopc_ModelAdvance:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    first_call = .false.

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

    !-----------------------------------------------------------------------
    use ESMF       , only : ESMF_ClockPrint
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
