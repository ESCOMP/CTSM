
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
  integer                       :: mpicom_clmdrv ! local mpicom
  integer                       :: ID_clmdrv     ! mct ID
  integer                       :: ncomps        ! number of separate components for MCT
  integer                       :: ntasks,mytask ! mpicom size and rank
  integer                       :: global_comm   ! copy of mpi_comm_world for pio
  integer,allocatable           :: comp_id(:)    ! for pio init2
  logical,allocatable           :: comp_iamin(:) ! for pio init2
  character(len=64),allocatable :: comp_name(:)  ! for pio init2
  integer,allocatable           :: comp_comm(:), comp_comm_iam(:) ! for pio_init2

  !----- Land Coupling Data -----
  !  type(seq_cdata)                :: cdata     ! Input land-model driver data
  !  type(seq_infodata_type),target :: infodata  ! infodata type
  !  type(mct_aVect)                :: x2l, l2x  ! land model import and export states
  !  type(mct_gGrid),target         :: dom_lnd   ! domain data for clm
  !  type(mct_gsMap),target         :: gsmap_lnd ! gsmap data for clm
  integer                        :: orb_iyear ! Orbital
  real*8                         :: orb_eccen, orb_obliq, orb_mvelp, orb_obliqr, orb_lambm0, orb_mvelpp
  character(len=128)             :: case_name, case_desc, model_version, hostname, username
  character(len=128)             :: start_type
  logical                        :: brnch_retain_casename, single_column, atm_aero
  real*8                         :: scmlat, scmlon
  integer :: idx_Sa_z, idx_Sa_u, idx_Sa_v, idx_Sa_tbot, idx_Sa_ptem, &
       idx_Sa_shum, idx_Sa_pbot, idx_Faxa_rainc, idx_Faxa_rainl, &
       idx_Faxa_snowc, idx_Faxa_snowl, idx_Faxa_lwdn, idx_Faxa_swndr, &
       idx_Faxa_swvdr, idx_Faxa_swndf, idx_Faxa_swvdf

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
  integer :: iunit = 250               ! clmdrv log unit number
  integer :: sunit = 249               ! share log unit number
  character(len=*),parameter :: subname = 'clmdrv'

  type fld_list_type
     character(len=128) :: stdname
  end type fld_list_type

  !----------------------------------------------

  !----------------------------------------------
  !--- MPI/MCT ---
  !----------------------------------------------

  call MPI_Init(ierr)
  call MPI_Comm_Dup(MPI_COMM_WORLD, mpicom_clmdrv, ierr)
  call MPI_COMM_RANK(mpicom_clmdrv, mytask, ierr)
  call MPI_COMM_SIZE(mpicom_clmdrv, ntasks, ierr)

  ncomps = 1
  ID_clmdrv = 1
  call mct_world_init(ncomps,MPI_COMM_WORLD,mpicom_clmdrv,ID_clmdrv)

  !----------------------------------------------
  !--- Log File and PIO ---
  !----------------------------------------------

  global_comm = MPI_COMM_WORLD
  call shr_pio_init1(ncomps, 'pio_in', global_comm)
  allocate(comp_id(ncomps),comp_name(ncomps),comp_iamin(ncomps),comp_comm(ncomps),comp_comm_iam(ncomps))
  do n = 1,ncomps
     comp_id(n) = ID_clmdrv
     comp_name(n) = 'LND'
     comp_iamin(n) = .true.
     comp_comm(n) = mpicom_clmdrv
     comp_comm_iam(n) = mytask
  enddo
  call shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)
  deallocate(comp_id,comp_name,comp_iamin,comp_comm,comp_comm_iam)

  write(string,'(a,i4.4)') 'clmdrv.log.',mytask
  open(iunit, file=trim(string))
  write(iunit,*) subname,' STARTING'
  call shr_sys_flush(iunit)

  write(iunit,*) subname,' ntasks = ',ntasks
  write(iunit,*) subname,' mytask = ',mytask
  write(iunit,*) subname,' mct ID = ',ID_clmdrv
  call shr_sys_flush(iunit)
  call shr_file_setLogUnit(sunit)
  call shr_file_setLogLevel(1)

  !----------------------------------------------
  !--- Clocks ---
  !----------------------------------------------

  call ESMF_Initialize(rc=rc)
  Calendar = ESMF_CalendarCreate( name='clmdrv_NOLEAP', &
       calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
  call ESMF_TimeSet(StartTime, yy=2000, mm=1, dd=1, s=0, calendar=Calendar, rc=rc)
  call ESMF_TimeSet(StopTime , yy=2000, mm=1, dd=10, s=0, calendar=Calendar, rc=rc)
  call ESMF_TimeIntervalSet(TimeStep, s=3600, rc=rc)
  EClock = ESMF_ClockCreate(name='clmdrv_EClock', &
       TimeStep=TimeStep, startTime=StartTime, &
       RefTime=StartTime, stopTime=stopTime, rc=rc)

  EAlarm_stop = ESMF_AlarmCreate(name='seq_timemgr_alarm_stop'   , &
       clock=EClock, ringTime=StopTime, rc=rc)
  EAlarm_rest = ESMF_AlarmCreate(name='seq_timemgr_alarm_restart', &
       clock=EClock, ringTime=StopTime, rc=rc)

  call ESMF_TimeGet( StartTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
  write(iunit,'(1x,2a,4i6)') subname,' StartTime ymds=',yy,mm,dd,sec
  call ESMF_TimeGet( StopTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
  write(iunit,'(1x,2a,4i6)') subname,' StopTime  ymds=',yy,mm,dd,sec
  call shr_sys_flush(iunit)

  !----------------------------------------------
  !--- Coupling ---
  !----------------------------------------------

  !--- coupling fields
  seq_flds_dom_coord='lat:lon'
  seq_flds_dom_other='area:aream:mask:frac'
  seq_flds_dom_fields=trim(seq_flds_dom_coord)//':'//trim(seq_flds_dom_other)

  seq_flds_x2l_states= 'Sa_z:Sa_u:Sa_v:Sa_tbot:Sa_ptem:Sa_shum:Sa_pbot:Sg_icemask:Sg_icemask_coupled_fluxes'
  seq_flds_x2l_fluxes= 'Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl:Faxa_lwdn:Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf:Faxa_bcphidry:Faxa_bcphodry:Faxa_bcphiwet:Faxa_ocphidry:Faxa_ocphodry:Faxa_ocphiwet:Faxa_dstwet1:Faxa_dstwet2:Faxa_dstwet3:Faxa_dstwet4:Faxa_dstdry1:Faxa_dstdry2:Faxa_dstdry3:Faxa_dstdry4:Flrr_flood:Flrr_volr'
  seq_flds_x2l_fields= trim(seq_flds_x2l_states)//':'//trim(seq_flds_x2l_fluxes)

  seq_flds_l2x_states= 'Sl_avsdr:Sl_anidr:Sl_avsdf:Sl_anidf:Sl_tref:Sl_qref:Sl_t:Sl_fv:Sl_ram1:Sl_snowh:Sl_u10'
  seq_flds_l2x_fluxes= 'Fall_swnet:Fall_taux:Fall_tauy:Fall_lat:Fall_sen:Fall_lwup:Fall_evap:Fall_flxdst1:Fall_flxdst2:Fall_flxdst3:Fall_flxdst4:Flrl_rofl:Flrl_rofi:Fall_voc001:Fall_voc002:Fall_voc003:Fall_voc004:Fall_voc005:Fall_voc006:Fall_voc007:Fall_voc008'
  seq_flds_l2x_fields= trim(seq_flds_l2x_states)//':'//trim(seq_flds_l2x_fluxes)

  !--- set orbital params
  orb_iyear = 1990
  call shr_orb_params(orb_iyear, orb_eccen, orb_obliq, orb_mvelp, &
       orb_obliqr, orb_lambm0, orb_mvelpp, .true.)
  !  call seq_infodata_putData(infodata, orb_eccen=orb_eccen, orb_mvelpp=orb_mvelpp, &
  !       orb_lambm0=orb_lambm0, orb_obliqr=orb_obliqr )

  !--- set case information
  case_name = 'clmdrv'
  case_desc = 'clmdrv with clm'
  model_version = 'clmdrv0.1'
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

  write(iunit,*) subname,' calling lnd_init_mct'
  call shr_sys_flush(iunit)
  !  call lnd_init_mct(Eclock, cdata, x2l, l2x)

  call diag_avect(l2x,mpicom_clmdrv,'l2x_init')

  idx_Sa_z       = mct_avect_indexra(x2l,'Sa_z')
  idx_Sa_u       = mct_avect_indexra(x2l,'Sa_u')
  idx_Sa_v       = mct_avect_indexra(x2l,'Sa_v')
  idx_Sa_tbot    = mct_avect_indexra(x2l,'Sa_tbot')
  idx_Sa_ptem    = mct_avect_indexra(x2l,'Sa_ptem')
  idx_Sa_shum    = mct_avect_indexra(x2l,'Sa_shum')
  idx_Sa_pbot    = mct_avect_indexra(x2l,'Sa_pbot')
  idx_Faxa_rainc = mct_avect_indexra(x2l,'Faxa_rainc')
  idx_Faxa_rainl = mct_avect_indexra(x2l,'Faxa_rainl')
  idx_Faxa_snowc = mct_avect_indexra(x2l,'Faxa_snowc')
  idx_Faxa_snowl = mct_avect_indexra(x2l,'Faxa_snowl')
  idx_Faxa_lwdn  = mct_avect_indexra(x2l,'Faxa_lwdn')
  idx_Faxa_swndr = mct_avect_indexra(x2l,'Faxa_swndr')
  idx_Faxa_swvdr = mct_avect_indexra(x2l,'Faxa_swvdr')
  idx_Faxa_swndf = mct_avect_indexra(x2l,'Faxa_swndf')
  idx_Faxa_swvdf = mct_avect_indexra(x2l,'Faxa_swvdf')

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


  state = ESMF_StateCreate(name=statename,  &
       stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)

  ! Create Field Bundle
  FBout = ESMF_FieldBundleCreate(name=trim(lname), rc=rc)

  ! Create individual states and add to field bundle
  field = ESMF_FieldCreate(lmesh, ESMF_TYPEKIND_R8, meshloc=meshloc, name="Sa_z", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! Add FB to state
  call ESMF_StateAdd(state, (/FBout/), rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

  ! fill in pointer with data
  call ESMF_StateGet(State, itemName="Sa_z", field=lfield, rc=rc)
  call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)

  ! then
  fldptr = 30.0

  !----------------------------------------------
  !--- Time Loop ---
  !----------------------------------------------

  call ESMF_ClockGet(Eclock, currTime=CurrTime, rc=rc)
  do while (CurrTime < StopTime)
     call ESMF_ClockAdvance(EClock, rc=rc)
     call ESMF_ClockGet(EClock, currTime=CurrTime, rc=rc)
     call ESMF_TimeGet( CurrTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
     write(iunit,'(1x,2a,4i6)') subname,' clmdrv ymds=',yy,mm,dd,sec
     call shr_sys_flush(iunit)

     ! can manually override the alarms as needed
     call ESMF_AlarmRingerOff(EAlarm_rest, rc=rc)
     if (mod(dd,5)==0 .and. sec==0) call ESMF_AlarmRingerOn(EAlarm_rest,rc)

     ! set the coupling data that is sent to the land model, this is on atm decomp
     ! this is just sample test data
     ! these all need to be set in the pointers
     Sa_z       = 30.0   ! m
     Sa_u       = 0.0    ! m/s
     Sa_v       = 0.0    ! m/s
     Sa_tbot    = 280.0  ! degK
     Sa_ptem    = 280.0  ! degK
     Sa_shum    = 0.0004 ! kg/kg
     Sa_pbot    = 100100.0 ! Pa
     Faxa_rainc = 4.0e-8 ! kg/m2s
     Faxa_rainl = 3.0e-8 ! kg/m2s
     Faxa_snowc = 1.0e-8 ! kg/m2s
     Faxa_snowl = 2.0e-8 ! kg/m2s
     Faxa_lwdn  = 200.0  ! W/m2
     Faxa_swndr = 100.0  ! W/m2
     Faxa_swvdr =  90.0  ! W/m2
     Faxa_swndf =  20.0  ! W/m2
     Faxa_swvdf =  40.0  ! W/m2

     ! run clm
     write(iunit,*) subname,' call lilac%run',yy,mm,dd,sec
     call lilac%run(importState, exportState, clock)
     !  call lnd_run_mct(Eclock, cdata, x2l, l2x)

  enddo

  !----------------------------------------------
  !--- lnd_final ---
  !----------------------------------------------

  write(iunit,*) subname,' calling lnd_final_mct'
  call shr_sys_flush(iunit)
  !  call lnd_final_mct(Eclock, cdata, x2l, l2x)

  !----------------------------------------------
  !--- Done ---
  !----------------------------------------------

  write(iunit,*) subname,' DONE'
  call shr_sys_flush(iunit)
  call MPI_Finalize(ierr)

  subroutine fldlist_add(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname


    ! local variables
    integer :: rc
    integer :: dbrc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:fldlist_add)'
    !-------------------------------------------------------------------------------


    ! Set up a list of field information


    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
       return
    endif
    fldlist(num)%stdname = trim(stdname)


  end subroutine fldlist_add

end program lilac_data_driver

