
program lilac_data_driver

  use seq_flds_mod , only:  &
       seq_flds_x2l_states, seq_flds_x2l_fluxes, seq_flds_x2l_fields, &
       seq_flds_l2x_states, seq_flds_l2x_fluxes, seq_flds_l2x_fields, &
       seq_flds_dom_coord, seq_flds_dom_other, seq_flds_dom_fields
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_putdata, seq_infodata_getdata
  use shr_sys_mod  , only: shr_sys_flush, shr_sys_abort
  use shr_orb_mod  , only: shr_orb_params
  use shr_file_mod , only: shr_file_setlogunit, shr_file_setloglevel
  use shr_pio_mod  , only: shr_pio_init1, shr_pio_init2
  use ESMF

  implicit none

#include <mpif.h>         ! mpi library include file

  !----- Clocks -----
  type(ESMF_Clock)           :: EClock           ! Input synchronization clock
  type(ESMF_Time)            :: CurrTime, StartTime, StopTime
  type(ESMF_TimeInterval)    :: TimeStep
  type(ESMF_Alarm)           :: EAlarm_stop, EAlarm_rest
  type(ESMF_Calendar),target :: Calendar
  integer                    :: yy,mm,dd,sec

  !----- MPI/MCT -----
  integer                       :: mpicom_lilac_drv ! local mpicom
  integer                       :: ID_lilac_drv     ! mct ID
  integer                       :: ncomps        ! number of separate components for MCT
  integer                       :: ntasks,mytask ! mpicom size and rank
  integer                       :: global_comm   ! copy of mpi_comm_world for pio
  integer,allocatable           :: comp_id(:)    ! for pio init2
  logical,allocatable           :: comp_iamin(:) ! for pio init2
  character(len=64),allocatable :: comp_name(:)  ! for pio init2
  integer,allocatable           :: comp_comm(:), comp_comm_iam(:) ! for pio_init2

  !----- Land Coupling Data -----
  type(ESMF_GridComp)  :: gridComp
  type(ESMF_State)     :: a2x_state
  type(ESMF_State)     :: x2a_state

  integer                        :: orb_iyear ! Orbital
  real*8                         :: orb_eccen, orb_obliq, orb_mvelp, orb_obliqr, orb_lambm0, orb_mvelpp
  character(len=128)             :: case_name, case_desc, model_version, hostname, username
  character(len=128)             :: start_type
  logical                        :: brnch_retain_casename, single_column, atm_aero
  real*8                         :: scmlat, scmlon

  !----- Atm Model -----
  integer :: atm_nx, atm_ny
  integer :: gsize, lsize, gstart, gend   ! domain decomp info
  integer, allocatable :: gindex(:)       ! domain decomp info
  type(mct_aVect)    :: x2l_a             ! data for land on atm decomp
  type(mct_aVect)    :: l2x_a             ! data from land on atm decomp
  type(mct_gsMap)    :: gsmap_atm         ! gsmap data for atm
  type(mct_rearr)    :: rearr_atm2lnd     ! rearranger for atm to land
  type(mct_rearr)    :: rearr_lnd2atm     ! rearranger for land to atm

  !----- Other -----
  integer :: n,m                       ! counter
  character(len=128) :: string         ! temporary string
  integer :: ierr, rc                  ! local error status
  integer :: iunit = 250               ! lilac_drv log unit number
  integer :: sunit = 249               ! share log unit number
  character(len=*),parameter :: subname = 'lilac_drv'

  type fld_list_type
     character(len=128) :: stdname
     real*8             :: default_value
     character(len=128) :: units
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToCpl_num = 0
  integer                :: fldsFrCpl_num = 0
  type (fld_list_type)   :: fldsToCpl(fldsMax)
  type (fld_list_type)   :: fldsFrCpl(fldsMax)

  !----------------------------------------------

  !----------------------------------------------
  !--- MPI/MCT ---
  !----------------------------------------------

  call MPI_Init(ierr)
  call MPI_Comm_Dup(MPI_COMM_WORLD, mpicom_lilac_drv, ierr)
  call MPI_COMM_RANK(mpicom_lilac_drv, mytask, ierr)
  call MPI_COMM_SIZE(mpicom_lilac_drv, ntasks, ierr)

  call lilac%start()

  !----------------------------------------------
  !--- Log File and PIO ---
  !----------------------------------------------

  global_comm = MPI_COMM_WORLD
  call shr_pio_init1(ncomps, 'pio_in', global_comm)
  allocate(comp_id(ncomps),comp_name(ncomps),comp_iamin(ncomps),comp_comm(ncomps),comp_comm_iam(ncomps))
  do n = 1,ncomps
     comp_id(n) = ID_lilac_drv
     comp_name(n) = 'LND'
     comp_iamin(n) = .true.
     comp_comm(n) = mpicom_lilac_drv
     comp_comm_iam(n) = mytask
  enddo
  call shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)
  deallocate(comp_id,comp_name,comp_iamin,comp_comm,comp_comm_iam)

  write(string,'(a,i4.4)') 'lilac_drv.log.',mytask
  open(iunit, file=trim(string))
  write(iunit,*) subname,' STARTING'
  call shr_sys_flush(iunit)

  write(iunit,*) subname,' ntasks = ',ntasks
  write(iunit,*) subname,' mytask = ',mytask
  write(iunit,*) subname,' mct ID = ',ID_lilac_drv
  call shr_sys_flush(iunit)
  call shr_file_setLogUnit(sunit)
  call shr_file_setLogLevel(1)

  !----------------------------------------------
  !--- Clocks ---
  !----------------------------------------------
  Calendar = ESMF_CalendarCreate(name='lilac_drv_NOLEAP', calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
  call ESMF_TimeSet(StartTime, yy=2000, mm=1, dd=1, s=0, calendar=Calendar, rc=rc)
  call ESMF_TimeSet(StopTime , yy=2000, mm=1, dd=10, s=0, calendar=Calendar, rc=rc)
  call ESMF_TimeIntervalSet(TimeStep, s=3600, rc=rc)
  EClock = ESMF_ClockCreate(name='lilac_drv_EClock', TimeStep=TimeStep, startTime=StartTime, RefTime=StartTime, stopTime=stopTime, rc=rc)

  EAlarm_stop = ESMF_AlarmCreate(name='seq_timemgr_alarm_stop', clock=EClock, ringTime=StopTime, rc=rc)
  EAlarm_rest = ESMF_AlarmCreate(name='seq_timemgr_alarm_restart', clock=EClock, ringTime=StopTime, rc=rc)

  call ESMF_TimeGet( StartTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
  write(iunit,'(1x,2a,4i6)') subname,' StartTime ymds=',yy,mm,dd,sec
  call ESMF_TimeGet( StopTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
  write(iunit,'(1x,2a,4i6)') subname,' StopTime  ymds=',yy,mm,dd,sec
  call shr_sys_flush(iunit)

  !--- set orbital params
  orb_iyear = 1990
  call shr_orb_params(orb_iyear, orb_eccen, orb_obliq, orb_mvelp, orb_obliqr, orb_lambm0, orb_mvelpp, .true.)

  !--- set case information
  case_name = 'lilac_drv'
  case_desc = 'lilac_drv with clm'
  model_version = 'lilac_drv0.1'
  hostname = 'undefined'
  username = 'undefined'
  start_type = 'startup'
  brnch_retain_casename = .true.
  single_column = .false.
  scmlat = 0.0
  scmlon = 0.0
  atm_aero = .true.
  call seq_infodata_putData(infodata, case_name=case_name,    &
       case_desc=case_desc, single_column=single_column, &
       scmlat=scmlat, scmlon=scmlon,                     &
       brnch_retain_casename=brnch_retain_casename,      &
       start_type=start_type, model_version=model_version, &
       hostname=hostname, username=username,             &
       atm_aero=atm_aero )

  !----------------------------------------------
  !--- lnd_init ---
  !----------------------------------------------

  write(iunit,*) subname,' calling lilac%init'
  call shr_sys_flush(iunit)

  call lilac%init(EClock, x2a_state, a2x_state, rc=rc)

  !----------------------------------------------
  !--- atm and atm/lnd coupling init ---
  !----------------------------------------------

  ! read in the mesh
  ! TODO: set cvalue to filepath of atm mesh
  cvalue = "/path/to/foo"

  if (masterproc) then
     write(iulog,*)'mesh file for domain is ',trim(cvalue)
  end if

  EMesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

  ! import fields
  ! call fldlist_add(fldsFrCpl_num, fldsFrCpl, trim(flds_scalar_name))

  ! land states
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_lfrin'      )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_t'          )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_tref'       )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_qref'       )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_avsdr'      )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_anidr'      )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_avsdf'      )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_anidf'      )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_snowh'      )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_u10'        )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_fv'         )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Sl_ram1'       )

  ! fluxes to atm
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_taux'     )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_tauy'     )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_lat'      )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_sen'      )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_lwup'     )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_evap'     )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_swnet'    )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst1'  )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst2'  )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst3'  )
  call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'Fall_flxdst4'  )

  ! call fldlist_add(fldsToCpl_num, fldsToCpl, trim(flds_scalar_name))

  ! from atm
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_z', default_value=30.0, units='m')
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_topo')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_u', default_value=0.0, units='m/s')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_v', default_value=0.0, units='m/s')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_ptem', default_value=280.0, 'degK')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_pbot', default_value=100100.0, units='Pa')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_tbot', default_value=280.0, units='degK')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_shum', default_value=0.0004, units='kg/kg')
  !call fldlist_add(fldsToCpl_num, fldsToCpl, 'Sa_methane'   )

  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_lwdn', default_value=200.0, units='W/m2')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_rainc', default_value=4.0e-8, units='kg/m2s')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_rainl', default_value=3.0e-8, units='kg/m2s')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_snowc', default_value=1.0e-8, units='kg/m2s')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_snowl', default_value=2.0e-8, units='kg/m2s')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swndr', default_value=100.0, units='W/m2')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swvdr', default_value=90.0, units='W/m2')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swndf', default_value=20.0, units='W/m2')
  call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_swvdf', default_value=40.0, units='W/m2')
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphidry')
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphodry')
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_bcphiwet')
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphidry')
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphodry')
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_ocphiwet')
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry1' )
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry2' )
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry3' )
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstdry4' )
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet1' )
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet2' )
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet3' )
  ! call fldlist_add(fldsToCpl_num, fldsToCpl, 'Faxa_dstwet4' )

  ! more: https://github.com/mvertens/ctsm/blob/ae02ffe25dbc4a85c769c9137b5b3d50f2843e89/src/cpl/nuopc/lnd_import_export.F90#L131

  ! Create States
  x2a_state = ESMF_StateCreate(name="x2a_state", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
  a2x_state = ESMF_StateCreate(name="x2a_state",  stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

  ! Coupler to Atmosphere Fields
  FBout = ESMF_FieldBundleCreate(name="x2a_fields", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
  
  ! Create individual states and add to field bundle
  do n = 1,fldsFrCpl_num
    ! create field
    field = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name=trim(fldsFrCpl(n)%stdname), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    ! add field to field bundle
    call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
  enddo

  ! Add FB to state
  call ESMF_StateAdd(x2a_state, (/FBout/), rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

  ! Atmosphere to Coupler Fields
  FBout = ESMF_FieldBundleCreate(name="a2x_fields", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
  
  ! Create individual states and add to field bundle
  do n = 1,fldsToCpl_num
    ! create field
    field = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name=trim(fldsToCpl(n)%stdname), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    ! initialize with default value
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    ldptr = fldsToCpl(n)%default_value

    ! add field to field bundle
    call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
  enddo

  ! Add FB to state
  call ESMF_StateAdd(a2x_state, (/FBout/), rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out


  !----------------------------------------------
  !--- Time Loop ---
  !----------------------------------------------

  call ESMF_ClockGet(Eclock, currTime=CurrTime, rc=rc)
  do while (CurrTime < StopTime)
    call ESMF_ClockAdvance(EClock, rc=rc)
    call ESMF_ClockGet(EClock, currTime=CurrTime, rc=rc)
    call ESMF_TimeGet( CurrTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
    write(iunit,'(1x,2a,4i6)') subname,' lilac_drv ymds=',yy,mm,dd,sec
    call shr_sys_flush(iunit)

    ! can manually override the alarms as needed
    call ESMF_AlarmRingerOff(EAlarm_rest, rc=rc)
    if (mod(dd,5)==0 .and. sec==0) call ESMF_AlarmRingerOn(EAlarm_rest,rc)

    ! run lilac
    write(iunit,*) subname,' call lilac%run',yy,mm,dd,sec
    call lilac%run(EClock, x2a_state, a2x_state, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
  enddo

  !----------------------------------------------
  !--- lnd_final ---
  !----------------------------------------------

  write(iunit,*) subname,' calling lilac%final()'
  call shr_sys_flush(iunit)
  call lilac%final()

  !----------------------------------------------
  !--- Done ---
  !----------------------------------------------

  write(iunit,*) subname,' DONE'
  call shr_sys_flush(iunit)
  call MPI_Finalize(ierr)

end program lilac_data_driver

subroutine fldlist_add(num, fldlist, stdname, default_value, units)
  integer                    intent(inout) :: num
  type(fld_list_type)        intent(inout) :: fldlist(:)
  character(len=*)           intent(in)    :: stdname
  real, optional             intent(in)    :: default_value
  character(len=*), optional intent(in)    :: units


  ! local variables
  integer :: rc
  character(len=*), parameter :: subname='(fldlist_add)'
  !-------------------------------------------------------------------------------

  ! Set up a list of field information
  num = num + 1
  if (num > fldsMax) then
     call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc) return
  endif
  fldlist(num)%stdname = trim(stdname)
  if(present(default_value)) then
    fldlist(num)%default_value = default_value
  else
    fldlist(num)%default_value = 0.
  end if
  if(present(units)) then
    fldlist(num)%units = trim(units)
  else
    fldlist(num)%units = ""
  end if

end subroutine fldlist_add
