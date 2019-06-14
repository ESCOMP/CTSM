
program lilac_data_driver

  use seq_infodata_mod, only: seq_infodata_putdata
  use shr_sys_mod  , only: shr_sys_flush
  use shr_orb_mod  , only: shr_orb_params
  use shr_file_mod , only: shr_file_setlogunit, shr_file_setloglevel
  use shr_pio_mod  , only: shr_pio_init1, shr_pio_init2
  use ESMF

  use lilac_utils , only create_fldlists, fldsMax


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
  type(LilacGrid)      :: gridComp
  type(LilacFields)    :: a2x_state
  type(LilacFields)    :: x2a_state

  integer                        :: orb_iyear ! Orbitalle
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

  call create_fldlists(fldsFrCpl, fldsToCpl)


  call lilac%init(EClock, x2a_state, a2x_state, rc=rc)

  !----------------------------------------------
  !--- atm and atm/lnd coupling init ---
  !----------------------------------------------

  !----------------------------------------------
  !--- Time Loop ---
  !----------------------------------------------

  call ESMF_ClockGet(Eclock, currTime=CurrTime, rc=rc)
  do while (CurrTime < StopTime)
     call ESMF_ClockAdvance(EClock, rc=rc)
     call ESMF_ClockGet(EClock, currTime=CurrTime, rc=rc)
     call ESMF_TimeGet(CurrTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
     write(iunit,'(1x,2a,4i6)') subname,' lilac_drv ymds=',yy,mm,dd,sec
     call shr_sys_flush(iunit)

     ! can manually override the alarms as needed
     call ESMF_AlarmRingerOff(EAlarm_rest, rc=rc)
     if (mod(dd,5)==0 .and. sec==0) call ESMF_AlarmRingerOn(EAlarm_rest,rc)

     ! run lilac
     write(iunit,*) subname,' call lilac%run',yy,mm,dd,sec
     call lilac%run(EClock, x2a_state, a2x_state, rc=rc)
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
