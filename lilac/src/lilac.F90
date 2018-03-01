module lilac
   !
   ! Public interface to lilac
   !

   implicit none
   private

   type, abstract :: lilac_t
      private
   contains
      ! Public API
      procedure :: Init => lilac_init
      procedure :: Shutdown => lilac_shutdown
      procudure :: AdvanceTime => lilac_advance_time

      ! private initialization routines
      procedure, private :: lilac_init_parallel
      procedure, private :: lilac_init_logging
      procedure, private :: lilac_init_io
      procedure, private :: lilac_init_clocks
      procedure, private :: lilac_init_fields
      procedure, private :: lilac_init_orbit
      procedure, private :: lilac_init_land
      procedure, private :: lilac_init_coupling

      ! private shudown routines
      procedure, private :: lilac_shutdown_land
      procedure, private :: lilac_shutdown_parallel

   end type lilac_t



contains

   !
   ! Public API
   !
   subroutine lilac_init(this)

      use mct_mod, only : mct_world_init

      implicit none

      class(lilac_t), intent(inout) :: this

      call this%lilac_init_parallel()

   end subroutine lilac_init

   subroutine lilac_advance_time(this)

      implicit none

      class(lilac_t), intent(inout) :: this

   end subroutine lilac_advance_time

   subroutine lilac_shutdown(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      ! FIXME(bja, 2018-02) master proc only!
      write(this%iunit, *) 'lilac shutting down...'
      call shr_sys_flush(this%iunit)

      call this%lilac_shutdown_land()
      call this%lilac_shutdown_parallel()

      ! FIXME(bja, 2018-02) master proc only!
      write(this%iunit, *) 'lilac shut down complete.'

   end subroutine lilac_shutdown

   !
   ! Private work functions
   !
   subroutine lilac_init_parallel(this)
      ! Initialize parallel components, e.g. MPI, MCT

      implicit none

      class(lilac_t), intent(inout) :: this

      call MPI_Init(ierr)
      call MPI_Comm_Dup(MPI_COMM_WORLD, mpicom_comp, ierr)
      call MPI_COMM_RANK(mpicom_comp, mytask, ierr)
      call MPI_COMM_SIZE(mpicom_comp, ntasks, ierr)

      num_comps = 1
      ID_comp = 1
      call mct_world_init(num_comps, MPI_COMM_WORLD, mpicom_lilac, ID_comp)

   end subroutine lilac_init_parallel

   subroutine lilac_init_logging(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      write(string,'(a,i4.4)') 'lilac.log.', mytask
      open(iunit, file=trim(string))
      write(iunit,*) subname,' STARTING'
      call shr_sys_flush(iunit)

      write(iunit,*) subname, ' ntasks = ', ntasks
      write(iunit,*) subname, ' mytask = ', mytask
      write(iunit,*) subname, ' mct ID = ', ID_comp
      call shr_sys_flush(iunit)
      call shr_file_setLogUnit(sunit)
      call shr_file_setLogLevel(1)

   end subroutine lilac_init_logging

   subroutine lilac_init_io(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      global_comm = MPI_COMM_WORLD
      call shr_pio_init1(ncomps, 'pio_in', global_comm)
      allocate(comp_id(ncomps), comp_name(ncomps), comp_iamin(ncomps), comp_comm(ncomps), comp_comm_iam(ncomps))
      do n = 1, ncomps
         comp_id(n) = ID_lilac
         comp_name(n) = 'LND'
         comp_iamin(n) = .true.
         comp_comm(n) = mpicom_lilac
         comp_comm_iam(n) = mytask
      enddo
      call shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)
      deallocate(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)

   end subroutine lilac_init_io

   subroutine lilac_init_clocks(this)

      use ESMF

      implicit none

      class(lilac_t), intent(inout) :: this

   type(ESMF_Clock) :: EClock           ! Input synchronization clock
   type(ESMF_Time)  :: CurrTime, StartTime, StopTime
   type(ESMF_TimeInterval) :: TimeStep
   type(ESMF_Alarm) :: EAlarm_stop, EAlarm_rest
   type(ESMF_Calendar), target :: Calendar

      call ESMF_Initialize(rc=rc)
      Calendar = ESMF_CalendarCreate( name='lilac_NOLEAP', &
           calkindflag=ESMF_CALKIND_NOLEAP, rc=rc )
      call ESMF_TimeSet(StartTime, yy=2000, mm=1, dd=1, s=0, calendar=Calendar, rc=rc)
      call ESMF_TimeSet(StopTime , yy=2000, mm=1, dd=10, s=0, calendar=Calendar, rc=rc)
      call ESMF_TimeIntervalSet(TimeStep, s=3600, rc=rc)
      EClock = ESMF_ClockCreate(name='lilac_EClock', &
           TimeStep=TimeStep, startTime=StartTime, &
           RefTime=StartTime, stopTime=stopTime, rc=rc)

      EAlarm_stop = ESMF_AlarmCreate(name='seq_timemgr_alarm_stop'   , &
           clock=EClock, ringTime=StopTime, rc=rc)
      EAlarm_rest = ESMF_AlarmCreate(name='seq_timemgr_alarm_restart', &
           clock=EClock, ringTime=StopTime, rc=rc)

      call ESMF_TimeGet( StartTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
      write(iunit,'(1x,2a,4i6)') subname,' StartTime ymds=', yy, mm, dd, sec
      call ESMF_TimeGet( StopTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
      write(iunit,'(1x,2a,4i6)') subname,' StopTime  ymds=', yy, mm, dd, sec
      call shr_sys_flush(iunit)

   end subroutine lilac_init_clocks

   subroutine lilac_init_fields(this)
      ! Set coupling fields.

      implicit none

      class(lilac_t), intent(inout) :: this

      ! FIXME(bja, 2018-02) this should be dynamically created at runtime
      ! instead of hard coded!

      seq_flds_dom_coord='lat:lon'
      seq_flds_dom_other='area:aream:mask:frac'
      seq_flds_dom_fields=trim(seq_flds_dom_coord)//':'//trim(seq_flds_dom_other)

      seq_flds_x2l_states= 'Sa_z:Sa_u:Sa_v:Sa_tbot:Sa_ptem:Sa_shum:Sa_pbot:Sg_icemask:Sg_icemask_coupled_fluxes'
      seq_flds_x2l_fluxes= 'Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl:Faxa_lwdn:Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf:Faxa_bcphidry:Faxa_bcphodry:Faxa_bcphiwet:Faxa_ocphidry:Faxa_ocphodry:Faxa_ocphiwet:Faxa_dstwet1:Faxa_dstwet2:Faxa_dstwet3:Faxa_dstwet4:Faxa_dstdry1:Faxa_dstdry2:Faxa_dstdry3:Faxa_dstdry4:Flrr_flood:Flrr_volr'
      seq_flds_x2l_fields= trim(seq_flds_x2l_states)//':'//trim(seq_flds_x2l_fluxes)

      seq_flds_l2x_states= 'Sl_avsdr:Sl_anidr:Sl_avsdf:Sl_anidf:Sl_tref:Sl_qref:Sl_t:Sl_fv:Sl_ram1:Sl_snowh:Sl_u10'
      seq_flds_l2x_fluxes= 'Fall_swnet:Fall_taux:Fall_tauy:Fall_lat:Fall_sen:Fall_lwup:Fall_evap:Fall_flxdst1:Fall_flxdst2:Fall_flxdst3:Fall_flxdst4:Flrl_rofl:Flrl_rofi:Fall_voc001:Fall_voc002:Fall_voc003:Fall_voc004:Fall_voc005:Fall_voc006:Fall_voc007:Fall_voc008'
      seq_flds_l2x_fields= trim(seq_flds_l2x_states)//':'//trim(seq_flds_l2x_fluxes)


   end subroutine lilac_init_fields

   subroutine lilac_init_orbit(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      !--- set orbital params
      orb_iyear = 1990
      call shr_orb_params(orb_iyear, orb_eccen, orb_obliq, orb_mvelp, &
           orb_obliqr, orb_lambm0, orb_mvelpp, .true.)
      call seq_infodata_putData(infodata, orb_eccen=orb_eccen, orb_mvelpp=orb_mvelpp, &
           orb_lambm0=orb_lambm0, orb_obliqr=orb_obliqr )


   end subroutine lilac_init_orbit

   subroutine lilac_init_land(this)

      implicit none

      write(iunit,*) subname,' calling lnd_init_mct'
      call shr_sys_flush(iunit)
      call lnd_init_mct(Eclock, cdata, x2l, l2x)

      call diag_avect(l2x, mpicom_lilac,'l2x_init')

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

   end subroutine lilac_init_land

   subroutine lilac_init_coupling(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      ! set atm grid size to land grid size in this example.  for a real
      ! atmosphere model, the atm and land grids should agree at the outset.
      call seq_infodata_getData(infodata, lnd_nx=atm_nx, lnd_ny=atm_ny)

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
      write(iunit,'(1x,2a,5i8)') subname,' atm decomp = ', mytask, gsize, lsize, gstart, gend

      ! initialize land grid on atm decomp
      call mct_gsMap_init(gsmap_atm, gindex, mpicom_lilac, ID_lilac, lsize, gsize)
      deallocate(gindex)

      ! initialize rearrangers between atm and land decomps
      call mct_rearr_init(gsmap_atm, gsmap_lnd, mpicom_lilac, rearr_atm2lnd)
      call mct_rearr_init(gsmap_lnd, gsmap_atm, mpicom_lilac, rearr_lnd2atm)

      ! initialize atm avects from land avects with atm lsize
      call mct_avect_init(x2l_a, x2l, lsize)
      call mct_avect_zero(x2l_a)
      call mct_avect_init(l2x_a, l2x, lsize)
      call mct_avect_zero(l2x_a)

   end subroutine lilac_init_coupling

   subroutine lilac_shutdown_land(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      write(iunit, *) 'lilac shutting down component ', this%comp_name
      call lnd_final_mct(Eclock, cdata, x2l, l2x)

   end subroutine lilac_shutdown_land

   subroutine lilac_shutdown_parallel(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      ! FIXME(bja, 2018-02) need to determine if it is our responsibility to shutdown mpi or the caller!?
      ! call MPI_Finalize(ierr)

   end subroutine lilac_shutdown_parallel

end module lilac
