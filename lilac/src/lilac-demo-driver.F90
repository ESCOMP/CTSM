program lilac_demo_driver

   use lnd_comp_mct , only: lnd_init_mct, lnd_run_mct, lnd_final_mct
   use seq_flds_mod , only:  &
        seq_flds_x2l_states, seq_flds_x2l_fluxes, seq_flds_x2l_fields, &
        seq_flds_l2x_states, seq_flds_l2x_fluxes, seq_flds_l2x_fields, &
        seq_flds_dom_coord, seq_flds_dom_other, seq_flds_dom_fields
   use seq_cdata_mod, only: seq_cdata
   use seq_infodata_mod, only: seq_infodata_type, seq_infodata_putdata, seq_infodata_getdata
   use shr_sys_mod  , only: shr_sys_flush, shr_sys_abort
   use shr_orb_mod  , only: shr_orb_params
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
   type(ESMF_Calendar), target :: Calendar
   integer :: yy, mm, dd, sec

   !----- MPI/MCT -----
   integer :: mpicom_lilac             ! local mpicom
   integer :: ID_lilac                 ! mct ID
   integer :: ncomps                    ! number of separate components for MCT
   integer :: ntasks, mytask             ! mpicom size and rank
   integer :: global_comm               ! copy of mpi_comm_world for pio
   integer, allocatable :: comp_id(:)    ! for pio init2
   logical, allocatable :: comp_iamin(:) ! for pio init2
   character(len=64), allocatable :: comp_name(:)  ! for pio init2
   integer, allocatable ::  comp_comm(:), comp_comm_iam(:) ! for pio_init2

   !----- Land Coupling Data -----
   type(seq_cdata)  :: cdata            ! Input land-model driver data
   type(seq_infodata_type), target :: infodata  ! infodata type
   type(mct_aVect)  :: x2l, l2x         ! land model import and export states
   type(mct_gGrid), target  :: dom_lnd   ! domain data for clm
   type(mct_gsMap), target  :: gsmap_lnd ! gsmap data for clm
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
   integer :: n, m                       ! counter
   character(len=128) :: string         ! temporary string
   integer :: ierr, rc                  ! local error status
   integer :: iunit = 250               ! lilac log unit number
   integer :: sunit = 249               ! share log unit number
   character(len=*), parameter :: subname = 'lilac_demo_driver'

   !----------------------------------------------
   type(lilac_init_data_t) :: lilac_init_data
   class(lilac_t) :: lilac

   !
   ! Initialize lilac
   !

   ! Where should these come from in general? namelist?
   call MPI_Comm_Dup(MPI_COMM_WORLD, lilac_init_data%mpicom_lilac, ierr)
   call MPI_Comm_Dup(MPI_COMM_WORLD, lilac_init_data%mpicom_component, ierr)
   lilac_init_data%output_unit_lilac = 250
   lilac_init_data%output_unit_component = 249


   lilac%Init(lilac_init_data)

   ! FIXME(bja, 2018-02) don't want to use the cdata structure, but we still
   ! need to provide this information to the component?!

   !--- set mpicom and cdata memory
   cdata%name     =  'cdata_lilac'
   cdata%ID       =  ID_lilac
   cdata%mpicom   =  mpicom_lilac
   cdata%dom      => dom_lnd
   cdata%gsmap    => gsmap_lnd
   cdata%infodata => infodata

   !--- set case information
   case_name = 'lilac'
   case_desc = 'lilac with clm'
   model_version = 'lilac-v0.1'
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
      call ESMF_TimeGet( CurrTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
      write(iunit,'(1x,2a,4i6)') subname,' lilac ymds=', yy, mm, dd, sec
      call shr_sys_flush(iunit)

      ! can manually override the alarms as needed
      call ESMF_AlarmRingerOff(EAlarm_rest, rc=rc)
      if (mod(dd, 5)==0 .and. sec==0) call ESMF_AlarmRingerOn(EAlarm_rest, rc)

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
      write(iunit,*) subname,' x2l fields: ', yy, mm, dd, sec
      !       call diag_avect(x2l_a, mpicom_lilac,'x2l_a')
      call diag_avect(x2l, mpicom_lilac,'x2l')

      ! run clm
      write(iunit,*) subname,' call lnd_run_mct', yy, mm, dd, sec
      call lnd_run_mct(Eclock, cdata, x2l, l2x)

      ! rearrange data from land decomposition
      call mct_rearr_rearrange(l2x, l2x_a, rearr_lnd2atm)

      ! diagnose
      write(iunit,*) subname,' l2x fields: ', yy, mm, dd, sec
      call diag_avect(l2x, mpicom_lilac,'l2x')
      !       call diag_avect(l2x_a, mpicom_lilac,'l2x_a')
   enddo

   lilac%Shutdown()

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
      integer           :: n, k         ! counters
      integer           :: npts, nptsg  ! number of local/global pts in AV
      integer           :: kflds       ! number of fields in AV
      real*8, pointer   :: sumbuf (:)  ! sum buffer
      real*8, pointer   :: sumbufg(:)  ! sum buffer reduced
      integer           :: iam         ! pe number
      type(mct_string)  :: mstring     ! mct char type
      character(len=128):: itemc       ! string converted to char

      !----- formats -----
      character(*), parameter :: subName = '(diag_avect) '

      !----------------------------------------------------------------

      npts  = mct_aVect_lsize(AV)
      kflds = mct_aVect_nRattr(AV)
      allocate(sumbuf(kflds), sumbufg(kflds))

      sumbuf = 0.0

      do k = 1, kflds
         do n = 1, npts
            sumbuf(k) = sumbuf(k) + (AV%rAttr(k, n))
         enddo
      enddo

      call MPI_REDUCE(sumbuf, sumbufg, kflds, MPI_REAL8, MPI_SUM, 0, mpicom, ierr)
      call MPI_COMM_RANK(mpicom, iam, ierr)

      if (iam == 0) then
         do k = 1, kflds
            call mct_aVect_getRList(mstring, k, AV)
            itemc = mct_string_toChar(mstring)
            call mct_string_clean(mstring)
            write(iunit, 101) trim(comment), k, sumbufg(k), trim(itemc)
         enddo
         call shr_sys_flush(iunit)
      endif

      deallocate(sumbuf, sumbufg)

101   format('comm_diag ', a, 1x, i3, es26.19, 1x, a)

   end subroutine diag_avect

   !======================================================================
end program lilac_demo_driver

