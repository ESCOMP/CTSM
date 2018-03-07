module lilac
   !
   ! Public interface to lilac
   !

   implicit none

   private

   integer, parameter :: LILAC_MASTER_PROC = 0
   integer, parameter :: LILAC_NUM_COMPONENTS = 1

   type, public :: lilac_t
      private
      character(len=STRING_32) :: component_name
      logical :: debug

      integer :: mpicom_lilac
      integer :: my_mpi_rank_lilac
      integer :: num_mpi_tasks_lilac
      integer :: mct_comp_id

      type(ESMF_Clock) :: lilac_clock
      type(ESMF_Time) :: start_time
      type(ESMF_Time) :: stop_time
      type(ESMF_TimeInterval) :: time_step
      type(ESMF_Alarm) :: alarm_stop
      type(ESMF_Alarm) :: alarm_restart
      
   contains
      ! Public API
      procedure, public :: Init => lilac_init
      procedure, public :: Shutdown => lilac_shutdown
      procudure, public :: AdvanceTime => lilac_advance_time

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
   subroutine lilac_init(this, init_data, clock_data, debug)

      use lilac_api_types, only : lilac_clock_data_t
      use mct_mod, only : mct_world_init

      implicit none

      class(lilac_t), intent(inout) :: this
      type(lilac_init_data_t), intent(in) :: init_data
      type(lilac_clock_data_t), intent(in) :: clock_data
      logical, intent(in) :: debug

      this%debug = debug
      this%component_name = init_data%component_name

      call this%lilac_init_parallel(init_data%mpicom_lilac, &
           init_data%mpicom_component, init_data%mpicom_global_shared)

      call this%lilac_init_logging(init_data%output_unit_lilac, init_data%output_unit_component)
      call this%lilac_init_io()
      call this%lilac_init_clocks(clock_data)
      ! TODO(bja, 2018-03) use init_data%component_name to do some model
      ! specific setup, including getting a list of hard coded input and output
      ! exchange fields.
      call this%lilac_init_fields()
      call this%lilac_init_orbit()
      call this%lilac_init_land()
      call this%lilac_init_coupling()

   end subroutine lilac_init

   subroutine lilac_advance_time(this)

      implicit none

      class(lilac_t), intent(inout) :: this

   end subroutine lilac_advance_time

   subroutine lilac_shutdown(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      if (this%my_mpi_rank_lilac == LILAC_MASTER_PROC) then
         write(this%output_unit, *) 'Shutting down lilac interface for component ', this%component_name, ' ...'
      end if

      call shr_sys_flush(this%output_unit)

      call this%lilac_shutdown_land()
      call this%lilac_shutdown_parallel()

      if (this%my_mpi_rank_lilac == LILAC_MASTER_PROC) then
         write(this%output_unit, *) 'lilac shut down for component ', this%component_name, ' complete.'
      end if

   end subroutine lilac_shutdown

   !
   ! Private work functions
   !
   subroutine lilac_init_parallel(this, mpicom_lilac, mpicom_component, mpicom_global_shared)
      ! Initialize parallel components, e.g. MPI, MCT

      implicit none

      class(lilac_t), intent(inout) :: this

      ! should be safe if previously initialized
      call MPI_Init(ierr)

      this%mpicom_lilac = mpicom_lilac
      this%mpicom_component = mpicom_component
      
      call MPI_COMM_RANK(this%mpicom_lilac, this%my_lilac_mpi_rank, ierr)
      call MPI_COMM_SIZE(this%mpicom_lilac, this%num_lilac_mpi_tasks, ierr)

      ! FIXME(bja, 2018-03) 1 (component | lilac) or two (component & lilac)?
      this%mct_num_comps = 1
      this%mct_comp_id = 1
      ! NOTE(bja, 2018-02) MPI_COMM_WORLD should eventually be initialized on
      ! the union of the lilac and component communicators! If 2, then need arrays?!
      call mct_world_init(this%mct_num_comps, MPI_COMM_WORLD, this%mpicom_lilac, this%mct_comp_id)

   end subroutine lilac_init_parallel

   subroutine lilac_init_logging(this, output_unit_lilac, output_unit_global_shared)

      implicit none

      class(lilac_t), intent(inout) :: this

      character(len=*), parameter :: subname = 'lilac_init_logging'
      ! open logfile for lilac

      this%output_unit = output_unit_lilac

      ! FIXME(bja, 2018-03) do we want a single shared log file, or one per rank?
      write(log_file_name,'(a,i4.4)') 'lilac.log.', this%my_mpi_rank_lilac
      open(this%output_unit, file=trim(log_file_name))
      if (this%my_mpi_rank_lilac == LILAC_MASTER_PROC) then
         write(this%output_unit, *) subname, ': Starting lilac interface for component: ', this%component_name
         write(this%output_unit, *) subname, ':   num lilac tasks = ', this%num_mpi_tasks_lilac
         write(this%output_unit, *) subname, ':   my mpi rank = ', this%my_mpi_rank_lilac
         write(this%output_unit, *) subname, ':   mct component ID = ', this%mct_comp_id_lilac
         call shr_sys_flush(this%output_unit)
      end if

      ! NOTE(bja, 2018-02) these are setting global variables within the shr code!
      call shr_file_setLogUnit(output_unit_global_shared)
      call shr_file_setLogLevel(1)

   end subroutine lilac_init_logging

   subroutine lilac_init_io(this)
      ! NOTE(bja, 2018-02) There is only a *single science component* in each
      ! lilac instance. For now assuming just the science component interacts
      ! with pio, but lilac may have some parallel data I/O needs. If so it
      ! needs to be added to these data structures!

      implicit none

      class(lilac_t), intent(inout) :: this

      !
      call shr_pio_init1(LILAC_NUM_COMPONENTS, 'pio_in', this%mpicom_lilac)
      allocate( &
           comp_id(LILAC_NUM_COMPONENTS), &
           comp_name(LILAC_NUM_COMPONENTS), &
           comp_iamin(LILAC_NUM_COMPONENTS), &
           comp_comm(LILAC_NUM_COMPONENTS), &
           comp_comm_iam(LILAC_NUM_COMPONENTS))

      index = 1
      comp_id(index) = 1
      comp_name(index) = MODEL_NAME_LILAC // '_' // trim(this%component_name)
      comp_iamin(index) = .true.
      comp_comm(index) = this%mpicom_lilac
      comp_comm_iam(index) = this%my_mpi_rank_lilac

      ! TODO(bja, 2018-03) Never have more than one science component, remove loop?
      do n = 1, LILAC_NUM_COMPONENTS
         index = index + n
         comp_id(index) = ID_component
         comp_name(index) = this%component_name
         comp_iamin(index) = .true.
         comp_comm(index) = this%mpicom_component
         comp_comm_iam(index) = mytask ! FIXME(bja, 2018-02) when land and lilac are on different comms??
      enddo

      call shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)

      deallocate(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)

   end subroutine lilac_init_io

   subroutine lilac_init_clocks(this, clock_data)

      use ESMF

      implicit none

      class(lilac_t), intent(inout) :: this
      type(lilac_clock_data_t), intent(in) :: clock_data

      type(ESMF_Calendar), target :: calendar ! FIXME(bja, 2018-02) does not need to be freed?!
      integer :: cal_kind_flag
      integer :: year, month, day, sec

      if (clock_data%calendar_is_leap == .false.) then
         cal_kind_flag = ESMF_CALKIND_NOLEAP
      else
         ! FIXME(bja, 2018-03) not implemented error! ESMF_CALKIND_GREGORIAN?
      end if

      if (ESMF_IsInitialized() /= .true.) then
         ! NOTE(bja, 2018-03) allocates and operates on global data!
         call ESMF_Initialize(rc=rc)
      end if

      calendar = ESMF_CalendarCreate( name='lilac', calkindflag=cal_kind_flag, rc=rc )
      call ESMF_TimeSet(this%start_time, yy=clock_data%start_year, mm=clock_data%start_month, &
           dd=clock_data%start_day, s=clock_data%start_seconds, calendar=calendar, rc=rc)

      call ESMF_TimeSet(this%stop_time , yy=clock_data%stop_year, mm=clock_data%stop_month, &
           dd=clock_data%stop_day, s=clock_data%stop_seconds, calendar=calendar, rc=rc)

      call ESMF_TimeIntervalSet(this%time_step, s=clock_data%time_step_seconds, rc=rc)

      this%lilac_clock = ESMF_ClockCreate(name='lilac_clock', &
           TimeStep=this%time_step, startTime=this%start_time, &
           RefTime=this%start_time, stopTime=this%stop_time, rc=rc)

      this%alarm_stop = ESMF_AlarmCreate(name='alarm_stop'   , &
           clock=this%lilac_clock, ringTime=this%stop_time, rc=rc)
      this%alarm_rest = ESMF_AlarmCreate(name='alarm_restart', &
           clock=this%lilac_clock, ringTime=this%stop_time, rc=rc)

      if (this%debug .and. this%my_mpi_rank_lilac == LILAC_MASTER_PROC) then
         call ESMF_TimeGet( start_time, yy=year, mm=month, dd=day, s=sec, rc=rc )
         write(this%output_unit, '(1x,2a,4i6)') subname,': start time ymds=', year, month, day, sec
         call ESMF_TimeGet( stop_time, yy=year, mm=month, dd=day, s=sec, rc=rc )
         write(this%output_unit, '(1x,2a,4i6)') subname,': stop time  ymds=', year, month, day, sec
         call shr_sys_flush(this%output_unit)
      end if

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

      write(this%output_unit,*) subname,' calling lnd_init_mct'
      call shr_sys_flush(this%output_unit)
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
      write(this%output_unit,'(1x,2a,5i8)') subname,' atm decomp = ', mytask, gsize, lsize, gstart, gend

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

      write(this%output_unit, *) 'lilac shutting down component ', this%comp_name
      call lnd_final_mct(Eclock, cdata, x2l, l2x)

   end subroutine lilac_shutdown_land

   subroutine lilac_shutdown_parallel(this)

      implicit none

      class(lilac_t), intent(inout) :: this

      ! FIXME(bja, 2018-02) need to determine if it is our responsibility to shutdown mpi or the caller!?
      ! call MPI_Finalize(ierr)

   end subroutine lilac_shutdown_parallel

end module lilac
