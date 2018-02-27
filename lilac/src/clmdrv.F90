
PROGRAM clmdrv

    use lnd_comp_mct , only: lnd_init_mct, lnd_run_mct, lnd_final_mct
    use seq_flds_mod , only:  &
        seq_flds_x2l_states, seq_flds_x2l_fluxes, seq_flds_x2l_fields, &
        seq_flds_l2x_states, seq_flds_l2x_fluxes, seq_flds_l2x_fields, &
        seq_flds_dom_coord, seq_flds_dom_other, seq_flds_dom_fields
    use seq_cdata_mod, only: seq_cdata
    use seq_infodata_mod, only: seq_infodata_type, seq_infodata_putdata, seq_infodata_getdata
    use shr_sys_mod  , only: shr_sys_flush, shr_sys_abort
    use shr_orb_mod  , only: shr_orb_params
    use shr_file_mod , only: shr_file_setlogunit, shr_file_setloglevel
    use shr_pio_mod  , only: shr_pio_init1, shr_pio_init2
    use mct_mod
    use ESMF

    implicit none

#include <mpif.h>         ! mpi library include file

    !----- Clocks -----
    type(ESMF_Clock) :: EClock           ! Input synchronization clock
    type(ESMF_Time)  :: CurrTime, StartTime, StopTime
    type(ESMF_TimeInterval) :: TimeStep
    type(ESMF_Alarm) :: EAlarm_stop, EAlarm_rest
    type(ESMF_Calendar),target :: Calendar
    integer :: yy,mm,dd,sec

    !----- MPI/MCT -----
    integer :: mpicom_clmdrv             ! local mpicom
    integer :: ID_clmdrv                 ! mct ID
    integer :: ncomps                    ! number of separate components for MCT
    integer :: ntasks,mytask             ! mpicom size and rank
    integer :: global_comm               ! copy of mpi_comm_world for pio
    integer,allocatable :: comp_id(:)    ! for pio init2
    logical,allocatable :: comp_iamin(:) ! for pio init2
    character(len=64),allocatable :: comp_name(:)  ! for pio init2
    integer,allocatable ::  comp_comm(:), comp_comm_iam(:) ! for pio_init2

    !----- Land Coupling Data -----
    type(seq_cdata)  :: cdata            ! Input land-model driver data
    type(seq_infodata_type),target :: infodata  ! infodata type
    type(mct_aVect)  :: x2l, l2x         ! land model import and export states
    type(mct_gGrid),target  :: dom_lnd   ! domain data for clm
    type(mct_gsMap),target  :: gsmap_lnd ! gsmap data for clm
    integer :: orb_iyear                 ! Orbital
    real*8 :: orb_eccen, orb_obliq, orb_mvelp, orb_obliqr, orb_lambm0, orb_mvelpp
    character(len=128) :: case_name, case_desc, model_version, hostname, username
    character(len=128) :: start_type
    logical :: brnch_retain_casename, single_column, atm_aero
    real*8 :: scmlat, scmlon
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

    !--- set mpicom and cdata memory
    cdata%name     =  'cdata_clmdrv'
    cdata%ID       =  ID_clmdrv
    cdata%mpicom   =  mpicom_clmdrv
    cdata%dom      => dom_lnd
    cdata%gsmap    => gsmap_lnd
    cdata%infodata => infodata

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
    call seq_infodata_putData(infodata, orb_eccen=orb_eccen, orb_mvelpp=orb_mvelpp, &
         orb_lambm0=orb_lambm0, orb_obliqr=orb_obliqr )

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
    call lnd_init_mct(Eclock, cdata, x2l, l2x)

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

    ! set atm grid size to land grid size in this example.  for a real
    ! atmosphere model, the atm and land grids should agree at the outset.
    call seq_infodata_getData(infodata,lnd_nx=atm_nx,lnd_ny=atm_ny)

    ! atm decomp
    gstart = ((mytask * atm_nx * atm_ny) / ntasks) + 1
    gend   = (((mytask+1) * atm_nx * atm_ny) / ntasks)
    lsize = gend - gstart + 1
    gsize = atm_nx * atm_ny
    allocate(gindex(lsize))
    do n = gstart, gend
       m = n-gstart+1
       gindex(m) = n
    end do
    write(iunit,'(1x,2a,5i8)') subname,' atm decomp = ',mytask,gsize,lsize,gstart,gend

    ! initialize land grid on atm decomp
    call mct_gsMap_init(gsmap_atm, gindex, mpicom_clmdrv, ID_clmdrv, lsize, gsize)
    deallocate(gindex)

    ! initialize rearrangers between atm and land decomps
    call mct_rearr_init(gsmap_atm, gsmap_lnd, mpicom_clmdrv, rearr_atm2lnd)
    call mct_rearr_init(gsmap_lnd, gsmap_atm, mpicom_clmdrv, rearr_lnd2atm)

    ! initialize atm avects from land avects with atm lsize
    call mct_avect_init(x2l_a, x2l, lsize)
    call mct_avect_zero(x2l_a)
    call mct_avect_init(l2x_a, l2x, lsize)
    call mct_avect_zero(l2x_a)

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
       x2l_a%rAttr(:,:) = 0.0
       x2l_a%rAttr(idx_Sa_z      ,:) = 30.0   ! m
       x2l_a%rAttr(idx_Sa_u      ,:) = 0.0    ! m/s
       x2l_a%rAttr(idx_Sa_v      ,:) = 0.0    ! m/s
       x2l_a%rAttr(idx_Sa_tbot   ,:) = 280.0  ! degK
       x2l_a%rAttr(idx_Sa_ptem   ,:) = 280.0  ! degK
       x2l_a%rAttr(idx_Sa_shum   ,:) = 0.0004 ! kg/kg
       x2l_a%rAttr(idx_Sa_pbot   ,:) = 100100.0 ! Pa
       x2l_a%rAttr(idx_Faxa_rainc,:) = 4.0e-8 ! kg/m2s
       x2l_a%rAttr(idx_Faxa_rainl,:) = 3.0e-8 ! kg/m2s
       x2l_a%rAttr(idx_Faxa_snowc,:) = 1.0e-8 ! kg/m2s
       x2l_a%rAttr(idx_Faxa_snowl,:) = 2.0e-8 ! kg/m2s
       x2l_a%rAttr(idx_Faxa_lwdn ,:) = 200.0  ! W/m2
       x2l_a%rAttr(idx_Faxa_swndr,:) = 100.0  ! W/m2
       x2l_a%rAttr(idx_Faxa_swvdr,:) =  90.0  ! W/m2
       x2l_a%rAttr(idx_Faxa_swndf,:) =  20.0  ! W/m2
       x2l_a%rAttr(idx_Faxa_swvdf,:) =  40.0  ! W/m2

       ! rearrange data to land decomposition
       call mct_rearr_rearrange(x2l_a, x2l, rearr_atm2lnd)

       ! diagnose
       write(iunit,*) subname,' x2l fields: ',yy,mm,dd,sec
!       call diag_avect(x2l_a,mpicom_clmdrv,'x2l_a')
       call diag_avect(x2l,mpicom_clmdrv,'x2l')

       ! run clm
       write(iunit,*) subname,' call lnd_run_mct',yy,mm,dd,sec
       call lnd_run_mct(Eclock, cdata, x2l, l2x)

       ! rearrange data from land decomposition
       call mct_rearr_rearrange(l2x, l2x_a, rearr_lnd2atm)

       ! diagnose
       write(iunit,*) subname,' l2x fields: ',yy,mm,dd,sec
       call diag_avect(l2x,mpicom_clmdrv,'l2x')
!       call diag_avect(l2x_a,mpicom_clmdrv,'l2x_a')
    enddo

    !----------------------------------------------
    !--- lnd_final ---
    !----------------------------------------------

    write(iunit,*) subname,' calling lnd_final_mct'
    call shr_sys_flush(iunit)
    call lnd_final_mct(Eclock, cdata, x2l, l2x)

    !----------------------------------------------
    !--- Done ---
    !----------------------------------------------

    write(iunit,*) subname,' DONE'
    call shr_sys_flush(iunit)
    call MPI_Finalize(ierr)

contains
!======================================================================

  SUBROUTINE diag_avect(av, mpicom, comment)

     use seq_infodata_mod

     implicit none

!    !INPUT/OUTPUT PARAMETERS:

     type(mct_aVect) , intent(in) :: av
     integer         , intent(in) :: mpicom
     character(len=*), intent(in) :: comment

     !--- local ---
     integer           :: n,k         ! counters
     integer           :: npts,nptsg  ! number of local/global pts in AV
     integer           :: kflds       ! number of fields in AV
     real*8, pointer   :: sumbuf (:)  ! sum buffer
     real*8, pointer   :: sumbufg(:)  ! sum buffer reduced
     integer           :: iam         ! pe number
     type(mct_string)  :: mstring     ! mct char type
     character(len=128):: itemc       ! string converted to char

     !----- formats -----
     character(*),parameter :: subName = '(diag_avect) '

     !----------------------------------------------------------------

     npts  = mct_aVect_lsize(AV)
     kflds = mct_aVect_nRattr(AV)
     allocate(sumbuf(kflds),sumbufg(kflds))

     sumbuf = 0.0

     do k = 1,kflds
     do n = 1,npts
        sumbuf(k) = sumbuf(k) + (AV%rAttr(k,n))
     enddo
     enddo

     call MPI_REDUCE(sumbuf,sumbufg,kflds,MPI_REAL8,MPI_SUM,0,mpicom,ierr)
     call MPI_COMM_RANK(mpicom,iam,ierr)

     if (iam == 0) then
        do k = 1,kflds
           call mct_aVect_getRList(mstring,k,AV)
           itemc = mct_string_toChar(mstring)
           call mct_string_clean(mstring)
           write(iunit,101) trim(comment),k,sumbufg(k),trim(itemc)
        enddo
        call shr_sys_flush(iunit)
     endif

     deallocate(sumbuf,sumbufg)

101  format('comm_diag ',a,1x,i3,es26.19,1x,a)

  end subroutine diag_avect

!======================================================================
end PROGRAM clmdrv

