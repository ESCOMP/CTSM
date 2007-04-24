#include <misc.h>
#include <preproc.h>

module clm_csmMod

#if (defined COUP_CSM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_csmMod
!
! !DESCRIPTION:
! Set of routines that define communication between the
! land model and flux coupler. The order of sends/receives is:
! 1) receive orbital data from coupler
! 2) send control data (grids and masks) to coupler
!    land grid does not have valid data, runoff grid does
! 3) receive valid land grid from flux coupler
! 4) send compressed runoff information to flux coupler
! 5) start normal send/recv communication patterm
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use nanMod
  use clm_varpar
  use clm_atmlnd      , only : clm_mapa2l, clm_mapl2a
  use clm_atmlnd      , only : atm_a2l, clm_a2l, clm_l2a, atm_l2a
  use spmdMod         , only : masterproc, mpicom
  use spmdGathScatMod , only : gather_data_to_master
  use cpl_fields_mod
  use cpl_contract_mod
  use cpl_interface_mod
  use RunoffMod        , only : runoff
  use shr_sys_mod      , only : shr_sys_irtc
  use abortutils       , only : endrun
  use perf_mod
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: csm_setup          ! Setup, mpi_init
  public :: csm_shutdown       ! Shutdown, mpi_finalize
  public :: csm_initialize     ! Initialize contracts, etc
  public :: csm_dosndrcv       ! Logic for determining if send/recv
  public :: csm_recv           ! Receive data from flux coupler
  public :: csm_send           ! Send data to flux coupler
  public :: csm_sendalb        ! Send initial albedos, surface temp and snow data
  public :: csm_flxave         ! Flux averaging rougine
  public :: csm_restart        ! Restart code
  public :: csm_compat         ! Checks compatibility of messages send/received
  public :: compat_check_spval ! Checks that data sent from the coupler is valid

! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!  03.01.15 T. Craig Update for cpl6
!  03.04.27 M. Vertenstein, added qref_2m to communication and
!           generalized global sums to include all fields
!
!EOP
!
  private
!
! PRIVATE MEMBER FUNCTIONS:
!
! PRIVATE TYPES:
!
  integer   :: nsend, nrecv, nroff            ! Buffer sizes
  integer   :: ibuffr(cpl_fields_ibuf_total)  ! Integer buffer from cpl
  integer   :: ibuffs(cpl_fields_ibuf_total)  ! Integer buffer to   cpl
  real(r8)  :: rbuffr(cpl_fields_rbuf_total)  ! Real    buffer from cpl
  real(r8)  :: rbuffs(cpl_fields_rbuf_total)  ! Real    buffer to   cpl
  type(cpl_contract)    :: contractRg         ! Contract for grid recvs from cpl
  type(cpl_contract)    :: contractR          ! Contract for recvs from cpl
  type(cpl_contract)    :: contractS          ! Contract for sends to   cpl
  type(cpl_contract)    :: contractSr         ! Contract for runoff sends to cpl
  real(r8), allocatable :: Gbuf(:,:)          ! Temporary generic buffer
  real(r8), allocatable :: bufS(:,:)          ! Send buffer for land
  real(r8), allocatable :: bufR(:,:)          ! Recv buffer for land
  real(r8), allocatable :: bufSr(:,:)         ! Send buffer for runoff
  real(r8), pointer     :: bufSglob(:,:)      ! Send global sum buffer for land
  real(r8), pointer     :: bufRglob(:,:)      ! Recv global sum buffer for land
  real(r8), pointer     :: bufSloc(:,:)       ! Send local sum buffer for land
  real(r8), pointer     :: bufRloc(:,:)       ! Recv local sum buffer for land
  real(r8), pointer     :: fieldS(:)          ! Global sum send field
  real(r8), pointer     :: fieldR(:)          ! Global sum receive field
  real(r8), pointer     :: area(:)            ! Global area field

  integer :: csm_nptg                         ! Loc sizes, grid coupling buffers
  integer :: csm_nptr                         ! Loc sizes, roff coupling buffers
  integer :: begp, endp                       ! per-proc beginning and ending pft indices
  integer :: begc, endc                       ! per-proc beginning and ending column indices
  integer :: begl, endl                       ! per-proc beginning and ending landunit indices
  integer :: begg, endg                       ! per-proc gridcell ending gridcell indices
  integer :: numg                             ! total number of gridcells across all processors
  integer :: numl                             ! total number of landunits across all processors
  integer :: numc                             ! total number of columns across all processors
  integer :: nump                             ! total number of pfts across all processors
!
! Flux averaging arrays and counters
!
  integer  :: icnt                         ! step counter for flux averager
  integer  :: ncnt                         ! number of steps over which to average output fluxes
  real(r8) :: rncnt                        ! reciprocal of ncnt

  real(r8), allocatable :: taux_ave(:)     ! averaged array
  real(r8), allocatable :: tauy_ave(:)     ! averaged array
  real(r8), allocatable :: lhflx_ave(:)    ! averaged array
  real(r8), allocatable :: shflx_ave(:)    ! averaged array
  real(r8), allocatable :: lwup_ave(:)     ! averaged array
  real(r8), allocatable :: qflx_ave(:)     ! averaged array
  real(r8), allocatable :: swabs_ave(:)    ! averaged array
  real(r8), allocatable :: nee_ave(:)      ! averaged array
!
! When to send/receive messages to coupler and when to make restart and stop
!
  integer, private:: ncpday                ! number of send/recv calls per day
  logical, public :: dorecv                ! receive data from coupler this step
  logical, public :: dosend                ! send data to coupler this step
  logical, public :: csmstop_next          ! received stop at eod signal and will stop on next ts
  logical, public :: csmstop_now           ! received stop now signal from coupler
  logical, public :: csmrstrt              ! restart write signal received from coupler
!
! CCSM timers
!
  logical  :: timer_lnd_sendrecv = .false. ! true => timer is on
  logical  :: timer_lnd_recvsend = .false. ! true => timer is on

! !PRIVATE MEMBER FUNCTIONS:
  private :: global_sum_fld1d    ! global sum of 1d fluxes

! Indices for send/recv fields

  integer :: index_l2c_Sl_t        ! temperature
  integer :: index_l2c_Sl_tref     ! 2m reference temperature
  integer :: index_l2c_Sl_qref     ! 2m reference specific humidity
  integer :: index_l2c_Sl_avsdr    ! albedo: direct , visible
  integer :: index_l2c_Sl_anidr    ! albedo: direct , near-ir
  integer :: index_l2c_Sl_avsdf    ! albedo: diffuse, visible
  integer :: index_l2c_Sl_anidf    ! albedo: diffuse, near-ir
  integer :: index_l2c_Sl_snowh    ! snow height
  integer :: index_l2c_Fall_taux   ! wind stress, zonal
  integer :: index_l2c_Fall_tauy   ! wind stress, meridional
  integer :: index_l2c_Fall_lat    ! latent          heat flux
  integer :: index_l2c_Fall_sen    ! sensible        heat flux
  integer :: index_l2c_Fall_lwup   ! upward longwave heat flux
  integer :: index_l2c_Fall_evap   ! evaporation    water flux
  integer :: index_l2c_Fall_swnet  ! 2m reference temperature
  integer :: index_l2c_Fall_nee    ! co2 flux

  integer :: index_c2l_Sa_co2prog  ! bottom atm prognostic co2
  integer :: index_c2l_Sa_co2diag  ! bottom atm diagnostic co2
  integer :: index_c2l_Sa_z        ! bottom atm level height
  integer :: index_c2l_Sa_u        ! bottom atm level zon wind
  integer :: index_c2l_Sa_v        ! bottom atm level mer wind
  integer :: index_c2l_Sa_tbot     ! bottom atm level temp
  integer :: index_c2l_Sa_ptem     ! bottom atm level pot temp
  integer :: index_c2l_Sa_shum     ! bottom atm level spec hum
  integer :: index_c2l_Sa_dens     ! bottom atm level air dens
  integer :: index_c2l_Sa_pbot     ! bottom atm level pressure
  integer :: index_c2l_Sa_pslv     ! sea level atm pressure     
  integer :: index_c2l_Faxa_lwdn   ! downward longwave heat flux
  integer :: index_c2l_Faxa_rainc  ! precip: liquid, convective
  integer :: index_c2l_Faxa_rainl  ! precip: liquid, large-scale
  integer :: index_c2l_Faxa_snowc  ! precip: frozen, convective
  integer :: index_c2l_Faxa_snowl  ! precip: frozen, large-scale
  integer :: index_c2l_Faxa_swndr  ! shortwave: nir direct  down
  integer :: index_c2l_Faxa_swvdr  ! shortwave: vis direct  down
  integer :: index_c2l_Faxa_swndf  ! shortwave: nir diffuse down
  integer :: index_c2l_Faxa_swvdf  ! shortwave: vis diffuse down

  integer :: index_r2c_Forr_roff   ! runoff to ocean
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_setup
!
! !INTERFACE:
  subroutine csm_setup(mpicom)
!
! !DESCRIPTION:
!  Initialize csm coupling, return the communicator group to the
!  application.
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: mpicom  !MPI group communicator
!
! !REVISION HISTORY:
!  03.01.15 T. Craig: initial version
!
!EOP
!
! !LOCAL VARIABLES:
!------------------------------------------------------------------------

   call cpl_interface_init(cpl_fields_lndname,mpicom)

 end subroutine csm_setup

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_shutdown
!
! !INTERFACE:
  subroutine csm_shutdown
!
! !DESCRIPTION:
!  Finalize csm coupling
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  03.01.15 T. Craig: initial version
!
!EOP
!
! !LOCAL VARIABLES:
!------------------------------------------------------------------------

   call cpl_interface_finalize(cpl_fields_lndname)

 end subroutine csm_shutdown

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_initialize
!
! !INTERFACE:
 subroutine csm_initialize(eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION:
! Initialize send/recv clm and rtm contracts with flux coupler
!
! !USES:
    use domainMod    , only : adomain
    use decompMod    , only : adecomp
    use decompMod    , only : get_proc_bounds, get_proc_global
    use RunoffMod    , only : get_proc_rof_bounds, runoff
    use clm_varctl   , only : csm_doflxave, nsrest, irad
    use clm_varcon   , only : re
    use clm_time_manager , only : get_step_size
    use shr_const_mod, only : SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(out) :: eccen  ! Earth's eccentricity of orbit
    real(r8), intent(out) :: obliqr ! Earth's obliquity in radians
    real(r8), intent(out) :: lambm0 ! Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(out) :: mvelpp ! Earth's moving vernal equinox long of perihelion plus pi (radians)
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!  03.01.15 T.Craig Update for cpl6
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: nsize                     ! attribute vector property	
    integer  :: g,gi,n,ni                 ! indices
    real(r8) :: dtime                     ! step size
    real(r8) :: spval                     ! special value
!------------------------------------------------------------------------

    !---------------------------------------------------------------
    ! Initialize module data for processor bounds
    !---------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)
!    call get_proc_rof_bounds(beg_lnd_rof, end_lnd_rof, beg_ocn_rof, end_ocn_rof)

    !---------------------------------------------------------------
    ! Setup contracts for send/recv communication
    !---------------------------------------------------------------

    ! Determine number of send/recv calls steps per day to flux coupler

    dtime = get_step_size()
    if (csm_doflxave) then
       ncpday = nint(SHR_CONST_CDAY/dtime)/irad
    else
       ncpday = nint(SHR_CONST_CDAY/dtime)
    endif

    ibuffs(:)  = 0                                   ! initialize ibuffs
    ibuffs(cpl_fields_ibuf_gsize  ) = adomain%ni*adomain%nj ! global array size
    ibuffs(cpl_fields_ibuf_gisize ) = adomain%ni     ! global number of lons
    ibuffs(cpl_fields_ibuf_gjsize ) = adomain%nj     ! global number of lats
    ibuffs(cpl_fields_ibuf_lsize  ) = endg-begg+1
    ibuffs(cpl_fields_ibuf_lisize ) = endg-begg+1
    ibuffs(cpl_fields_ibuf_ljsize ) = 1
    ibuffs(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
    ibuffs(cpl_fields_ibuf_ncpl   ) = ncpday         ! number of land send/recv calls per day

    allocate(Gbuf(begg:endg,cpl_fields_grid_total))

    do n = begg, endg	
!        i = adecomp%gdc2i(n)
!        j = adecomp%gdc2j(n)
!        gi = (j-1)*adomain%ni + i
        gi = adecomp%gdc2glo(n)
        Gbuf(n,cpl_fields_grid_lon)   = adomain%lonc(n)
        Gbuf(n,cpl_fields_grid_lat)   = adomain%latc(n)
        Gbuf(n,cpl_fields_grid_area)  = adomain%area(n)/(re*re)
        Gbuf(n,cpl_fields_grid_frac)  = adomain%frac(n)
        Gbuf(n,cpl_fields_grid_mask)  = float(adomain%mask(n))
        Gbuf(n,cpl_fields_grid_index) = gi
    end do

    ! Send ibuff and local grid information to flux coupler.  Initialize cpl6 contracts

    call cpl_interface_contractInit(contractS, cpl_fields_lndname, cpl_fields_cplname, &
         cpl_fields_l2c_fields, ibuffs, Gbuf)
    call cpl_interface_contractInit(contractR, cpl_fields_lndname, cpl_fields_cplname, &
         cpl_fields_c2l_fields, ibuffs, Gbuf)

    deallocate(Gbuf)

    !---------------------------------------------------------------
    ! Setup contracts for  runoff communication
    !---------------------------------------------------------------

    ibuffs(:) = 0
    ibuffs(cpl_fields_ibuf_gsize  ) = rtmlon*rtmlat   ! global array size
    ibuffs(cpl_fields_ibuf_gisize ) = rtmlon          ! global number of lons
    ibuffs(cpl_fields_ibuf_gjsize ) = rtmlat          ! global number of lats
    ibuffs(cpl_fields_ibuf_ncpl   ) = ncpday          ! number of land send/recv calls per day
    ibuffs(cpl_fields_ibuf_lsize  ) = runoff%lnumro   ! local array size
    ibuffs(cpl_fields_ibuf_lisize ) = runoff%lnumro   ! local array size
    ibuffs(cpl_fields_ibuf_ljsize ) = 1               ! local array size

    csm_nptr = ibuffs(cpl_fields_ibuf_lsize)

    allocate(Gbuf(csm_nptr,cpl_fields_grid_total))

    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          gi = runoff%gdc2glo(n)
          if (ni > runoff%lnumro) then
             write(6,*)'clm_csmMod: ERROR runoff count',n,ni,runoff%lnumro,gi
             call endrun()
          endif
          Gbuf(ni,cpl_fields_grid_lon  ) = runoff%lonc(n)
          Gbuf(ni,cpl_fields_grid_lat  ) = runoff%latc(n)
          Gbuf(ni,cpl_fields_grid_area ) = runoff%area(n)*1.0e-6_r8/(re*re)
          Gbuf(ni,cpl_fields_grid_mask ) = 1.0_r8 - float(runoff%mask(n))
          Gbuf(ni,cpl_fields_grid_index) = gi
       endif
    end do
    if (ni /= runoff%lnumro) then
       write(6,*)'clm_csmMod: ERROR runoff total count',ni,runoff%lnumro
       call endrun()
    endif

    call cpl_interface_contractInit(contractSr,cpl_fields_lndname,cpl_fields_cplname,cpl_fields_r2c_fields, ibuffs,Gbuf)

    deallocate(Gbuf)

    !---------------------------------------------------------------
    ! Receive initial ibuf message
    !---------------------------------------------------------------

    call cpl_interface_ibufRecv(cpl_fields_cplname,ibuffr,rbuffr)

    spval  = rbuffr(cpl_fields_rbuf_spval)
    eccen  = rbuffr(cpl_fields_rbuf_eccen)
    obliqr = rbuffr(cpl_fields_rbuf_obliqr)
    lambm0 = rbuffr(cpl_fields_rbuf_lambm0)
    mvelpp = rbuffr(cpl_fields_rbuf_mvelpp)

    ! Check that data is good data and not the special value

    if (masterproc) then
       call compat_check_spval(spval, eccen ,'Eccentricity'     )
       call compat_check_spval(spval, obliqr,'Obliquity'        )
       call compat_check_spval(spval, lambm0,'Long of perhelion')
       call compat_check_spval(spval, mvelpp,'Move long of perh')

       write(6,*)'(CSM_INITIALIZE): eccen:  ', eccen
       write(6,*)'(CSM_INITIALIZE): obliqr: ', obliqr
       write(6,*)'(CSM_INITIALIZE): lambm0: ', lambm0
       write(6,*)'(CSM_INITIALIZE): mvelpp: ', mvelpp
    end if

    write(6,*)'(CSM_INITIALIZE): there will be ',ncpday, &
         ' send/recv calls per day from the land model to the flux coupler'
    write(6,*)'(CSM_INITIALIZE):sent l->d control data '

    ! Allocate memory for grid and runoff communication

    nsend = cpl_interface_contractNumatt(contractS)
    allocate(bufS(begg:endg,nsend))

    nrecv = cpl_interface_contractNumatt(contractR)
    allocate(bufR(begg:endg,nrecv))

    nroff = cpl_interface_contractNumatt(contractSr)
    allocate(bufSr(csm_nptr,nroff))

    !---------------------------------------------------------------
    ! Determine indices
    !---------------------------------------------------------------

    ! Determine send indices

    index_l2c_Sl_t       = cpl_interface_contractIndex(contractS,'Sl_t')
    index_l2c_Sl_snowh   = cpl_interface_contractIndex(contractS,'Sl_snowh')
    index_l2c_Sl_avsdr   = cpl_interface_contractIndex(contractS,'Sl_avsdr')
    index_l2c_Sl_anidr   = cpl_interface_contractIndex(contractS,'Sl_anidr')
    index_l2c_Sl_avsdf   = cpl_interface_contractIndex(contractS,'Sl_avsdf')
    index_l2c_Sl_anidf   = cpl_interface_contractIndex(contractS,'Sl_anidf')
    index_l2c_Sl_tref    = cpl_interface_contractIndex(contractS,'Sl_tref')
    index_l2c_Sl_qref    = cpl_interface_contractIndex(contractS,'Sl_qref')
    index_l2c_Fall_taux  = cpl_interface_contractIndex(contractS,'Fall_taux')
    index_l2c_Fall_tauy  = cpl_interface_contractIndex(contractS,'Fall_tauy')
    index_l2c_Fall_lat   = cpl_interface_contractIndex(contractS,'Fall_lat')
    index_l2c_Fall_sen   = cpl_interface_contractIndex(contractS,'Fall_sen')
    index_l2c_Fall_lwup  = cpl_interface_contractIndex(contractS,'Fall_lwup')
    index_l2c_Fall_evap  = cpl_interface_contractIndex(contractS,'Fall_evap')
    index_l2c_Fall_swnet = cpl_interface_contractIndex(contractS,'Fall_swnet')
    index_l2c_Fall_nee   = cpl_interface_contractIndex(contractS,'Fall_nee', perrwith='quiet')

    ! Determine receive indices

    index_c2l_Sa_z       = cpl_interface_contractIndex(contractR,'Sa_z')
    index_c2l_Sa_u       = cpl_interface_contractIndex(contractR,'Sa_u')
    index_c2l_Sa_v       = cpl_interface_contractIndex(contractR,'Sa_v')
    index_c2l_Sa_ptem    = cpl_interface_contractIndex(contractR,'Sa_ptem')
    index_c2l_Sa_shum    = cpl_interface_contractIndex(contractR,'Sa_shum')
    index_c2l_Sa_pbot    = cpl_interface_contractIndex(contractR,'Sa_pbot')
    index_c2l_Sa_tbot    = cpl_interface_contractIndex(contractR,'Sa_tbot')
    index_c2l_Faxa_lwdn  = cpl_interface_contractIndex(contractR,'Faxa_lwdn')
    index_c2l_Faxa_rainc = cpl_interface_contractIndex(contractR,'Faxa_rainc')
    index_c2l_Faxa_rainl = cpl_interface_contractIndex(contractR,'Faxa_rainl')
    index_c2l_Faxa_snowc = cpl_interface_contractIndex(contractR,'Faxa_snowc')
    index_c2l_Faxa_snowl = cpl_interface_contractIndex(contractR,'Faxa_snowl')
    index_c2l_Faxa_swndr = cpl_interface_contractIndex(contractR,'Faxa_swndr')
    index_c2l_Faxa_swvdr = cpl_interface_contractIndex(contractR,'Faxa_swvdr')
    index_c2l_Faxa_swndf = cpl_interface_contractIndex(contractR,'Faxa_swndf')
    index_c2l_Faxa_swvdf = cpl_interface_contractIndex(contractR,'Faxa_swvdf')
    index_c2l_Sa_co2prog = cpl_interface_contractIndex(contractR,'Sa_co2prog', perrwith='quiet')
    index_c2l_Sa_co2diag = cpl_interface_contractIndex(contractR,'Sa_co2diag', perrwith='quiet')

    ! Determine runoff (send) indices

    index_r2c_Forr_roff  = cpl_interface_contractIndex(contractSR,'Forr_roff')

  end subroutine csm_initialize

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_sendalb
!
! !INTERFACE:
  subroutine csm_sendalb()
!
! !DESCRIPTION:
! Send initial albedos, surface temperature and snow data to the
! flux coupler
!
! !USES:
    use clm_varctl  , only : csm_doflxave, nsrest
    use clm_varcon  , only : sb
    use clm_time_manager, only : get_curr_date, get_prev_date, get_nstep
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,gi,n,ni ! generic indices
    integer :: yr            ! current year
    integer :: mon           ! current month
    integer :: day           ! current day (0, 1, ...)
    integer :: ncsec         ! current seconds of current date (0, ..., 86400)
    integer :: ncdate        ! current date (yymmdd format) (e.g., 021105)
! -----------------------------------------------------------------

    ! Fill send buffer for grid data

    bufS(:,:) = 0.0_r8

    if (nsrest == 0 ) then   !initial run

       call clm_mapl2a(clm_l2a, atm_l2a)

       do g = begg,endg
          bufS(g,index_l2c_Sl_t)     = sqrt(sqrt(atm_l2a%eflx_lwrad_out(g)/sb))
          bufS(g,index_l2c_Sl_snowh) = atm_l2a%h2osno(g)
          bufS(g,index_l2c_Sl_avsdr) = atm_l2a%albd(g,1)
          bufS(g,index_l2c_Sl_anidr) = atm_l2a%albd(g,2)
          bufS(g,index_l2c_Sl_avsdf) = atm_l2a%albi(g,1)
          bufS(g,index_l2c_Sl_anidf) = atm_l2a%albi(g,2)
       end do

    else  ! restart run

       ! On a restart run, no meaningful data is sent to the flux coupler -
       ! this includes ocean runoff (which should only contain zero values)
       ! since the runoff code (riverfluxrtm) has not been called yet

       bufS(:,:) = 1.e30_r8

    endif

    ! Determine time index to send to coupler. Note that for a restart run,
    ! the next time step is nstep+1. But must send current time step to flux couper here.

    if (nsrest == 0) then
       call get_curr_date (yr, mon, day, ncsec)
    else
       call get_prev_date (yr, mon, day, ncsec)
    endif

    ncdate = yr*10000 + mon*100 + day

    ibuffs(:)  = 0
    ibuffs(cpl_fields_ibuf_cdate)  = ncdate      !model date (yyyymmdd)
    ibuffs(cpl_fields_ibuf_sec  )  = ncsec       !elapsed seconds in current date

    ! Send grid data to coupler

    call cpl_interface_contractSend(cpl_fields_cplname,contractS,ibuffs,bufS)

    ! Send runoff data to coupler
    ! Must convert runoff to units of kg/m^2/s from m^3/s

    ni = 0
    bufSr(:,:) = 0._r8
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          bufSr(ni,index_r2c_Forr_roff) = runoff%runoff(n)/(runoff%area(n)*1.0e-6_r8*1000._r8)
          if (ni > runoff%lnumro) then
             write(6,*)'clm_csmMod: ERROR runoff count',n,ni,gi
             call endrun()
          endif
       endif
    end do
    if (ni /= runoff%lnumro) then
       write(6,*)'clm_csmMod: ERROR runoff total count',ni,runoff%lnumro
       call endrun()
    endif

    call cpl_interface_contractSend(cpl_fields_cplname,contractSr,ibuffs,bufSr)

  end subroutine csm_sendalb

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_dosndrcv
!
! !INTERFACE:
  subroutine csm_dosndrcv(doalb)
!
! !DESCRIPTION:
! Determine when to send and receive messages to/from the
! flux coupler on this time-step.
! Determine if send/receive information to/from flux coupler
! Send msgs (land state and fluxes) to the flux coupler only when
! doalb is true (i.e. on time steps before the atm does a solar
! radiation computation). Receive msgs (atm state) from the
! flux coupler only when dorad is true (i.e. on time steps
! when the atm does a solar radiation computation).
! The fluxes are then averaged between the send and receive calls.
!
! !USES:
    use clm_varctl   , only : csm_doflxave
    use clm_time_manager , only : get_step_size, get_nstep
    use shr_const_mod, only : SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
    logical, intent(in) :: doalb  !true=>next timestep a radiation time step
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ntspday           !model steps per day
    real(r8) :: dtime             !step size (seconds)
    integer  :: nstep             !time step
!-----------------------------------------------------------------------

    ! Determine if send/receive information to/from flux coupler

    nstep = get_nstep()
    if (csm_doflxave) then
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = doalb
       else
          dorecv = dosend
          dosend = doalb
       endif
    else
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = .true.
       else
          dorecv = .true.
          dosend = .true.
       endif
    endif

    ! If at end of day: check if should write restart file or stop
    ! at next time step. Note, these statements must appear here since
    ! ibuffr is not received at every time step when flux averaging occurs.

    csmstop_next = .false.
    csmrstrt     = .false.
    dtime        = get_step_size()
    ntspday      = nint(SHR_CONST_CDAY/dtime)
    if (mod(nstep,ntspday) == 0) then
       if (ibuffr(cpl_fields_ibuf_stopeod) /= 0) then  !stop at end of day
          csmstop_next = .true.  !will stop on next time step
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
       if (ibuffr(cpl_fields_ibuf_resteod) /= 0) then !write restart at end of day
          csmrstrt = .true.      !will write restart now
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
    endif

  end subroutine csm_dosndrcv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_recv
!
! !INTERFACE:
  subroutine csm_recv()
!
! !DESCRIPTION:
!  Receive and map data from flux coupler
!
! !USES:
    use clm_varctl, only : co2_type
    use clm_varcon, only : rair, o2_molar_const, co2_ppmv_const, c13ratio
    use clmtype   , only : nameg
    use domainMod , only : adomain
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,gi,n,ni ! generic indices
    real(r8):: forc_rainc    ! rainxy Atm flux mm/s
    real(r8):: forc_rainl    ! rainxy Atm flux mm/s
    real(r8):: forc_snowc    ! snowfxy Atm flux  mm/s
    real(r8):: forc_snowl    ! snowfxl Atm flux  mm/s
    real(r8):: co2_ppmv_diag ! temporary
    real(r8):: co2_ppmv_prog ! temporary
    real(r8):: co2_ppmv      ! temporary
    integer :: ier           ! return error code
!-----------------------------------------------------------------------

    ! Start timers

     if (timer_lnd_sendrecv) then
        call t_stopf ('lnd_sendrecv') ; timer_lnd_sendrecv = .false.
     endif

     call t_startf('lnd_recv')

     ibuffr(:) = 0

     call cpl_interface_contractRecv(cpl_fields_cplname,contractR,ibuffr,bufR)

     ! Do global integrals of fluxes if flagged

     if (ibuffr(cpl_fields_ibuf_infobug) >= 2) then

        if (masterproc) then
           if (.not. associated(bufRglob)) then
              allocate(bufRglob(nrecv,numg), stat=ier)
              if (ier /= 0) then
                 write(6,*)'clm_csmMod: allocation error for bufRglob'; call endrun
              end if
              if (.not. associated(fieldR)) then
                 allocate(fieldR(numg), stat=ier)
                 if (ier /= 0) then
                    write(6,*)'clm_csmMod: allocation error for fieldR'; call endrun
                 end if
              endif
              if (.not. associated(area)) then
                 allocate(area(numg), stat=ier)
                 if (ier /= 0) then
                    write(6,*)'clm_csmMod: allocation error for area'; call endrun
                 end if
              endif
           end if
        end if
        if (.not. associated(bufRloc)) then
           allocate(bufRloc(nrecv,begg:endg), stat=ier)
           if (ier /= 0) then
              write(6,*)'clm_csmMod: allocation error for bufRloc'; call endrun
           end if
        end if

        ! Find attribute vector indices

        do g = begg,endg
           do n = 1,nrecv
              bufRloc(n,g) = bufR(g,n)
           end do
        end do
        call gather_data_to_master(bufRloc, bufRglob, clmlevel=nameg)
        call gather_data_to_master(adomain%area, area, clmlevel=nameg)
        if (masterproc) then
           write(6,*)

           if (index_c2l_Sa_co2prog /= 0) then
              fieldR(:) = bufRglob(index_c2l_Sa_co2prog,:)
              write(6,100) 'lnd','recv', index_c2l_Sa_co2prog, global_sum_fld1d(fieldR,area), ' co2prog'
           end if

           if (index_c2l_Sa_co2diag /= 0) then
              fieldR(:) = bufRglob(index_c2l_Sa_co2diag,:)
              write(6,100) 'lnd','recv', index_c2l_Sa_co2diag, global_sum_fld1d(fieldR,area), ' co2diag'
           end if

           fieldR(:) = bufRglob(index_c2l_Sa_z,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_z, global_sum_fld1d(fieldR,area), ' hgt'

           fieldR(:) = bufRglob(index_c2l_Sa_u,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_u, global_sum_fld1d(fieldR,area), ' u'

           fieldR(:) = bufRglob(index_c2l_Sa_v,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_v, global_sum_fld1d(fieldR,area), ' v'

           fieldR(:) = bufRglob(index_c2l_Sa_ptem,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_ptem, global_sum_fld1d(fieldR,area), ' th'

           fieldR(:) = bufRglob(index_c2l_Sa_shum,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_shum, global_sum_fld1d(fieldR,area), ' q'

           fieldR(:) = bufRglob(index_c2l_Sa_pbot,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_pbot, global_sum_fld1d(fieldR,area), ' pbot'

           fieldR(:) = bufRglob(index_c2l_Sa_tbot,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_tbot, global_sum_fld1d(fieldR,area), ' t'

           fieldR(:) = bufRglob(index_c2l_Faxa_lwdn,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_lwdn, global_sum_fld1d(fieldR,area), ' lwrad'

           fieldR(:) = bufRglob(index_c2l_Faxa_rainc,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_rainc, global_sum_fld1d(fieldR,area), ' rainc'

           fieldR(:) = bufRglob(index_c2l_Faxa_rainl,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_rainl, global_sum_fld1d(fieldR,area), ' rainl'

           fieldR(:) = bufRglob(index_c2l_Faxa_snowc,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_snowc, global_sum_fld1d(fieldR,area), ' snowc'

           fieldR(:) = bufRglob(index_c2l_Faxa_snowl,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_snowl, global_sum_fld1d(fieldR,area), ' snowl'

           fieldR(:) = bufRglob(index_c2l_Faxa_swndr,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_swndr, global_sum_fld1d(fieldR,area),' soll '

           fieldR(:) = bufRglob(index_c2l_Faxa_swvdr,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_swvdr, global_sum_fld1d(fieldR,area),' sols '

           fieldR(:) = bufRglob(index_c2l_Faxa_swndf,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_swndf, global_sum_fld1d(fieldR,area),' solld'

           fieldR(:) = bufRglob(index_c2l_Faxa_swvdf,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_swvdf, global_sum_fld1d(fieldR,area),' solsd'

100        format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
           write(6,*)
        endif
     endif

     ! Stop timer

     call t_stopf('lnd_recv')

     ! Check if end of run now, if so stop (each processor does this)

     csmstop_now = .false.
     if (ibuffr(cpl_fields_ibuf_stopnow) /= 0) then
        csmstop_now = .true.
        if (timer_lnd_recvsend) call t_stopf('lnd_recvsend')
        if (timer_lnd_sendrecv) call t_stopf('lnd_sendrecv')
        write(6,*)'(CSM_RECV) stop now signal from flux coupler'
        write(6,*)'(CSM_RECV) ibuffr(cpl_fields_ibuf_stopnow) = ',ibuffr(cpl_fields_ibuf_stopnow)
        if (masterproc) then
           write(6,9001)
           write(6,9002) ibuffr(cpl_fields_ibuf_cdate)
           write(6,9003)
9001       format(/////' ===========> Terminating CLM Model')
9002       format(     '      Date: ',i8)
9003       format(/////' <=========== CLM Model Terminated')
        endif
        RETURN
     endif

     ! More timer logic

     if (.not. timer_lnd_recvsend) then
        call t_startf('lnd_recvsend') ; timer_lnd_recvsend = .true.
     endif

     ! Split data from coupler into component arrays.
     ! Note that the precipitation fluxes received  from the coupler
     ! are in units of kg/s/m^2. To convert these precipitation rates
     ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
     ! by 1000 mm/m resulting in an overall factor of unity.
     ! Below the units are therefore given in mm/s.

     if (co2_type == 'prognostic' .and. index_c2l_Sa_co2prog == 0) then
        write(6,*)' must have nonzero index_c2l_Sa_co2prog for co2_type equal to prognostic'
        call endrun()
     else if (co2_type == 'diagnostic' .and. index_c2l_Sa_co2diag == 0) then
        write(6,*)' must have nonzero index_c2l_Sa_co2diag for co2_type equal to diagnostic'
        call endrun()
     end if

    
     do g = begg,endg

        ! Determine required receive fields

        atm_a2l%forc_hgt(g)     = bufR(g,index_c2l_Sa_z)         ! zgcmxy  Atm state m
        atm_a2l%forc_u(g)       = bufR(g,index_c2l_Sa_u)         ! forc_uxy  Atm state m/s
        atm_a2l%forc_v(g)       = bufR(g,index_c2l_Sa_v)         ! forc_vxy  Atm state m/s
        atm_a2l%forc_th(g)      = bufR(g,index_c2l_Sa_ptem)      ! forc_thxy Atm state K
        atm_a2l%forc_q(g)       = bufR(g,index_c2l_Sa_shum)      ! forc_qxy  Atm state kg/kg
        atm_a2l%forc_pbot(g)    = bufR(g,index_c2l_Sa_pbot)      ! ptcmxy  Atm state Pa
        atm_a2l%forc_t(g)       = bufR(g,index_c2l_Sa_tbot)      ! forc_txy  Atm state K
        atm_a2l%forc_lwrad(g)   = bufR(g,index_c2l_Faxa_lwdn)    ! flwdsxy Atm flux  W/m^2
        forc_rainc               = bufR(g,index_c2l_Faxa_rainc)   ! mm/s
        forc_rainl               = bufR(g,index_c2l_Faxa_rainl)   ! mm/s
        forc_snowc               = bufR(g,index_c2l_Faxa_snowc)   ! mm/s
        forc_snowl               = bufR(g,index_c2l_Faxa_snowl)   ! mm/s
        atm_a2l%forc_solad(g,2) = bufR(g,index_c2l_Faxa_swndr)   ! forc_sollxy  Atm flux  W/m^2
        atm_a2l%forc_solad(g,1) = bufR(g,index_c2l_Faxa_swvdr)   ! forc_solsxy  Atm flux  W/m^2
        atm_a2l%forc_solai(g,2) = bufR(g,index_c2l_Faxa_swndf)   ! forc_solldxy Atm flux  W/m^2
        atm_a2l%forc_solai(g,1) = bufR(g,index_c2l_Faxa_swvdf)   ! forc_solsdxy Atm flux  W/m^2

        ! Determine optional receive fields

        if (index_c2l_Sa_co2prog /= 0) then
           co2_ppmv_prog = bufR(g,index_c2l_Sa_co2prog)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv_const
        end if
 
        if (index_c2l_Sa_co2diag /= 0) then
           co2_ppmv_diag = bufR(g,index_c2l_Sa_co2diag)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv_const
        end if

        ! Determine derived quantities for required fields

        atm_a2l%forc_hgt_u(g) = atm_a2l%forc_hgt(g)    !observational height of wind [m]
        atm_a2l%forc_hgt_t(g) = atm_a2l%forc_hgt(g)    !observational height of temperature [m]
        atm_a2l%forc_hgt_q(g) = atm_a2l%forc_hgt(g)    !observational height of humidity [m]
        atm_a2l%forc_vp(g)    = atm_a2l%forc_q(g) * atm_a2l%forc_pbot(g) &
                                / (0.622_r8 + 0.378_r8 * atm_a2l%forc_q(g))
        atm_a2l%forc_rho(g)   = (atm_a2l%forc_pbot(g) - 0.378_r8 * atm_a2l%forc_vp(g)) &
                                / (rair * atm_a2l%forc_t(g))
        atm_a2l%forc_po2(g)   = o2_molar_const * atm_a2l%forc_pbot(g)
        atm_a2l%forc_wind(g)  = sqrt(atm_a2l%forc_u(g)**2 + atm_a2l%forc_v(g)**2)
        atm_a2l%forc_solar(g) = atm_a2l%forc_solad(g,1) + atm_a2l%forc_solai(g,1) + &
                                atm_a2l%forc_solad(g,2) + atm_a2l%forc_solai(g,2)
        atm_a2l%forc_rain(g)  = forc_rainc + forc_rainl
        atm_a2l%forc_snow(g)  = forc_snowc + forc_snowl
        
        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        if (co2_type == 'prognostic') then
           co2_ppmv = co2_ppmv_prog
        else if (co2_type == 'diagnostic') then
           co2_ppmv = co2_ppmv_diag 
        else
           co2_ppmv = co2_ppmv_const      
        end if
        atm_a2l%forc_pco2(g) = co2_ppmv * 1.e-6_r8 * atm_a2l%forc_pbot(g) 
        ! 4/14/05: PET
        ! Adding isotope code
        atm_a2l%forc_pc13o2(g) = co2_ppmv * c13ratio * 1.e-6_r8 * atm_a2l%forc_pbot(g)

     end do

     call clm_mapa2l(atm_a2l, clm_a2l)

     ! debug write statements (remove)

!    if (masterproc) write(6,*)'co2_type = ', co2_type, ' co2_ppmv = ', co2_ppmv

  end subroutine csm_recv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_send
!
! !INTERFACE:
  subroutine csm_send()
!
! !DESCRIPTION:
! Send data to the flux coupler
!
! !USES:
    use clm_varctl  , only : csm_doflxave
    use clm_varcon  , only : sb
    use clm_time_manager, only : get_curr_date, get_nstep
    use clmtype     , only : nameg
    use domainMod   , only : adomain
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,n,ni,g,c,p  ! indices
    integer :: yr              ! current year
    integer :: mon             ! current month
    integer :: day             ! current day (0, 1, ...)
    integer :: ncsec           ! current seconds of current date (0, ..., 86400)
    integer :: ncdate          ! current date (yymmdd format) (e.g., 021105)
    integer :: ier             ! error status
!-----------------------------------------------------------------------


    ! Send data to the flux coupler

    if (timer_lnd_recvsend) then
       call t_stopf ('lnd_recvsend') ; timer_lnd_recvsend = .false.
    endif

    ! Start timer

    call t_startf('lnd_send')

    ! Recalculate fluxes to send to atm in flux averaged case
    ! Also recalculate a new t_rad based on flux averged lwup

    if (csm_doflxave) then
       do g = begg,endg
          clm_l2a%taux(g)           = taux_ave(g)  
          clm_l2a%tauy(g)           = tauy_ave(g)  
          clm_l2a%eflx_lh_tot(g)    = lhflx_ave(g) 
          clm_l2a%eflx_sh_tot(g)    = shflx_ave(g) 
          clm_l2a%eflx_lwrad_out(g) = lwup_ave(g)  
          clm_l2a%qflx_evap_tot(g)  = qflx_ave(g)  
          clm_l2a%fsa(g)            = swabs_ave(g) 
          clm_l2a%t_rad(g)          = (abs(clm_l2a%eflx_lwrad_out(g)/sb))**0.25_r8
          if (index_l2c_Fall_nee /= 0) then
             clm_l2a%nee(g)         = nee_ave(g) 
          end if
       end do
    endif

    call clm_mapl2a(clm_l2a, atm_l2a)

    bufS(:,:) = 0.0_r8
    do g = begg,endg
       bufS(g,index_l2c_Sl_t)       =  atm_l2a%t_rad(g)
       bufS(g,index_l2c_Sl_snowh)   =  atm_l2a%h2osno(g)
       bufS(g,index_l2c_Sl_avsdr)   =  atm_l2a%albd(g,1)
       bufS(g,index_l2c_Sl_anidr)   =  atm_l2a%albd(g,2)
       bufS(g,index_l2c_Sl_avsdf)   =  atm_l2a%albi(g,1)
       bufS(g,index_l2c_Sl_anidf)   =  atm_l2a%albi(g,2)
       bufS(g,index_l2c_Sl_tref)    =  atm_l2a%t_ref2m(g)
       bufS(g,index_l2c_Sl_qref)    =  atm_l2a%q_ref2m(g)
       bufS(g,index_l2c_Fall_taux)  = -atm_l2a%taux(g)
       bufS(g,index_l2c_Fall_tauy)  = -atm_l2a%tauy(g)
       bufS(g,index_l2c_Fall_lat)   = -atm_l2a%eflx_lh_tot(g)
       bufS(g,index_l2c_Fall_sen)   = -atm_l2a%eflx_sh_tot(g)
       bufS(g,index_l2c_Fall_lwup)  = -atm_l2a%eflx_lwrad_out(g)
       bufS(g,index_l2c_Fall_evap)  = -atm_l2a%qflx_evap_tot(g)
       bufS(g,index_l2c_Fall_swnet) = -atm_l2a%fsa(g)
       if (index_l2c_Fall_nee /= 0) then
          bufS(g,index_l2c_Fall_nee)  = -atm_l2a%nee(g)
       end if
    end do

    call get_curr_date (yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day

    ibuffs(:)  = 0
    ibuffs(cpl_fields_ibuf_cdate)  = ncdate      ! model date (yyyymmdd)
    ibuffs(cpl_fields_ibuf_sec  )  = ncsec       ! elapsed seconds in current date

    call cpl_interface_contractSend(cpl_fields_cplname,contractS ,ibuffs,bufS)

    ! Must convert runoff to units of kg/m^2/s from m^3/s

    ni = 0
    bufSr(:,:) = 0._r8
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          bufSr(ni,index_r2c_Forr_roff) = runoff%runoff(n)/(runoff%area(n)*1.0e-6_r8*1000._r8)
          if (ni > runoff%lnumro) then
             write(6,*)'clm_csmMod: ERROR runoff count',n,ni,runoff%lnumro
             call endrun()
          endif
       endif
    end do
    if (ni /= runoff%lnumro) then
       write(6,*)'clm_csmMod: ERROR runoff total count',ni,runoff%lnumro
       call endrun()
    endif

    call cpl_interface_contractSend(cpl_fields_cplname, contractSr, ibuffs, bufSr)

    ! Do global integrals if flag is set

    if (ibuffr(cpl_fields_ibuf_infobug) >= 2) then

       if (masterproc) then
          if (.not. associated(bufSglob)) then
             allocate(bufSglob(nsend,numg), stat=ier)
             if (ier /= 0) then
                write(6,*)'clm_csmMod: allocation error for bufSglob'; call endrun
             end if
          end if
          if (.not. associated(fieldS)) then
             allocate(fieldS(numg), stat=ier)
             if (ier /= 0) then
                write(6,*)'clm_csmMod: allocation error for fieldS'; call endrun
             end if
          endif
          if (.not. associated(area)) then
             allocate(area(numg), stat=ier)
             if (ier /= 0) then
                write(6,*)'clm_csmMod: allocation error for area'; call endrun
             end if
          endif
       end if
       if (.not. associated(bufSloc)) then
          allocate(bufSloc(nsend,begg:endg), stat=ier)
          if (ier /= 0) then
             write(6,*)'clm_csmMod: allocation error for bufSloc'; call endrun
          end if
       end if

       do g = begg,endg
          do n = 1,nsend
             bufSloc(n,g) = bufS(g,n)
          end do
       end do
       call gather_data_to_master(bufSloc, bufSglob, clmlevel=nameg)
       call gather_data_to_master(adomain%area, area, clmlevel=nameg)

       if (masterproc) then
          write(6,*)
          
          fieldS(:) = bufSglob(index_l2c_Sl_t,:)
          write(6,100) 'lnd','send', index_l2c_Sl_t, global_sum_fld1d(fieldS,area),' trad'

          fieldS(:) = bufSglob(index_l2c_Sl_avsdr,:)
          write(6,100) 'lnd','send', index_l2c_Sl_avsdr, global_sum_fld1d(fieldS,area),' asdir'

          fieldS(:) = bufSglob(index_l2c_Sl_anidr,:)
          write(6,100) 'lnd','send', index_l2c_Sl_anidr, global_sum_fld1d(fieldS,area),' aldir'

          fieldS(:) = bufSglob(index_l2c_Sl_avsdf,:)
          write(6,100) 'lnd','send', index_l2c_Sl_avsdf, global_sum_fld1d(fieldS,area),' asdif'

          fieldS(:) = bufSglob(index_l2c_Sl_anidf,:)
          write(6,100) 'lnd','send', index_l2c_Sl_anidf, global_sum_fld1d(fieldS,area),' aldif'

          fieldS(:) = bufSglob(index_l2c_Fall_taux,:)
          write(6,100) 'lnd','send', index_l2c_Fall_taux, global_sum_fld1d(fieldS,area),' taux'

          fieldS(:) = bufSglob(index_l2c_Fall_tauy,:)
          write(6,100) 'lnd','send', index_l2c_Fall_tauy, global_sum_fld1d(fieldS,area),' tauy'

          fieldS(:) = bufSglob(index_l2c_Fall_lat,:)
          write(6,100) 'lnd','send', index_l2c_Fall_lat, global_sum_fld1d(fieldS,area),' lhflx'

          fieldS(:) = bufSglob(index_l2c_Fall_sen,:)
          write(6,100) 'lnd','send', index_l2c_Fall_sen, global_sum_fld1d(fieldS,area),' shflx'

          fieldS(:) = bufSglob(index_l2c_Fall_lwup,:)
          write(6,100) 'lnd','send', index_l2c_Fall_lwup, global_sum_fld1d(fieldS,area),' lwup'

          fieldS(:) = bufSglob(index_l2c_Fall_evap,:)
          write(6,100) 'lnd','send', index_l2c_Fall_evap, global_sum_fld1d(fieldS,area),' qflx'

          fieldS(:) = bufSglob(index_l2c_Fall_swnet,:)
          write(6,100) 'lnd','send', index_l2c_Fall_swnet, global_sum_fld1d(fieldS,area),' swabs'

          if (index_l2c_Fall_nee /= 0) then
             fieldS(:) = bufSglob(index_l2c_Fall_nee,:)
             write(6,100) 'lnd','send', index_l2c_Fall_nee, global_sum_fld1d(fieldS,area),' nee'
          end if

          write(6,*)
100       format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
       endif
    endif

    ! Stop timers

    call t_stopf('lnd_send')

    if (.not. timer_lnd_recvsend) then
       call t_startf('lnd_sendrecv') ; timer_lnd_sendrecv = .true.
    endif

  end subroutine csm_send

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_flxave
!
! !INTERFACE:
  subroutine csm_flxave()
!
! !DESCRIPTION:
! Average output fluxes for flux coupler
! Add land surface model output fluxes to accumulators every time step.
! When icnt==ncnt, compute the average flux over the time interval.
!
! !USES:
!
    use clmtype      , only : clm3
    use clm_varctl   , only : irad
    use clm_time_manager , only : get_nstep
    use subgridAveMod, only : p2g, c2g
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g                    ! indices
    integer :: nstep                ! model time step
    real(r8), pointer :: taux_pft(:)
    real(r8), pointer :: taux_gcell(:)
    real(r8), pointer :: tauy_pft(:)
    real(r8), pointer :: tauy_gcell(:)
    real(r8), pointer :: eflx_lh_tot_pft(:)
    real(r8), pointer :: eflx_lh_tot_gcell(:)
    real(r8), pointer :: eflx_sh_tot_pft(:)
    real(r8), pointer :: eflx_sh_tot_gcell(:)
    real(r8), pointer :: eflx_lwrad_out_pft(:)
    real(r8), pointer :: eflx_lwrad_out_gcell(:)
    real(r8), pointer :: qflx_evap_tot_pft(:)
    real(r8), pointer :: qflx_evap_tot_gcell(:)
    real(r8), pointer :: fsa_pft(:)
    real(r8), pointer :: fsa_gcell(:)
    real(r8), pointer :: nee_pft(:)
    real(r8), pointer :: nee_col(:)
    real(r8), pointer :: nee_gcell(:)
!-----------------------------------------------------------------------

    ! Set pointers into derived type (pft level)

    taux_pft             => clm3%g%l%c%p%pmf%taux
    tauy_pft             => clm3%g%l%c%p%pmf%tauy
    eflx_lh_tot_pft      => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_sh_tot_pft      => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_lwrad_out_pft   => clm3%g%l%c%p%pef%eflx_lwrad_out
    qflx_evap_tot_pft    => clm3%g%l%c%p%pwf%qflx_evap_tot
    fsa_pft              => clm3%g%l%c%p%pef%fsa
#if (defined CN)
    nee_col              => clm3%g%l%c%ccf%nee
#elif (defined CASA)
    nee_pft              => clm3%g%l%c%p%pps%co2flux
#else
    nee_pft              => clm3%g%l%c%p%pcf%fco2
#endif

    ! Set pointers into derived type (gridcell level)

    taux_gcell           => clm_l2a%taux
    tauy_gcell           => clm_l2a%tauy
    eflx_lh_tot_gcell    => clm_l2a%eflx_lh_tot
    eflx_sh_tot_gcell    => clm_l2a%eflx_sh_tot
    eflx_lwrad_out_gcell => clm_l2a%eflx_lwrad_out
    qflx_evap_tot_gcell  => clm_l2a%qflx_evap_tot
    fsa_gcell            => clm_l2a%fsa
    nee_gcell            => clm_l2a%nee

    ! Allocate dynamic memory if necessary
    if (.not. allocated(taux_ave)) then
       allocate (taux_ave(begg:endg)) ; taux_ave(:) = nan
    endif
    if (.not. allocated(tauy_ave)) then
       allocate (tauy_ave(begg:endg)) ; tauy_ave(:) = nan
    endif
    if (.not. allocated(lhflx_ave)) then
       allocate (lhflx_ave(begg:endg)); lhflx_ave(:) = nan
    endif
    if (.not. allocated(shflx_ave)) then
       allocate (shflx_ave(begg:endg)); shflx_ave(:) = nan
    endif
    if (.not. allocated(lwup_ave)) then
       allocate (lwup_ave(begg:endg)) ; lwup_ave(:) = nan
    endif
    if (.not. allocated(qflx_ave)) then
       allocate (qflx_ave(begg:endg)) ; qflx_ave(:) = nan
    endif
    if (.not. allocated(swabs_ave)) then
       allocate (swabs_ave(begg:endg)) ; swabs_ave(:) = nan
    endif
    if (index_l2c_Fall_nee /= 0) then
       if (.not. allocated(nee_ave)) then
          allocate (nee_ave(begg:endg)) ; nee_ave(:) = nan
       endif
    end if
       
    ! Determine output flux averaging interval

    nstep = get_nstep()
    if (dorecv) then
       icnt = 1
       if ( nstep==0 ) then
          ncnt = irad + 1
       else
          ncnt = irad
       endif
       rncnt = 1._r8/ncnt
    endif

    ! Do flux averaging

    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, taux_pft, taux_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, tauy_pft, tauy_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, eflx_lh_tot_pft, eflx_lh_tot_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, eflx_sh_tot_pft, eflx_sh_tot_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, eflx_lwrad_out_pft, eflx_lwrad_out_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, qflx_evap_tot_pft, qflx_evap_tot_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, fsa_pft, fsa_gcell, &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    if (index_l2c_Fall_nee /= 0) then
#if (defined CN)
       call c2g(begc, endc, begl, endl, begg, endg, nee_col, nee_gcell, &
            c2l_scale_type= 'unity', l2g_scale_type='unity')
#elif (defined CASA)
       call p2g(begp, endp, begc, endc, begl, endl, begg, endg, nee_pft, nee_gcell, &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
#else
       call p2g(begp, endp, begc, endc, begl, endl, begg, endg, nee_pft, nee_gcell, &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
       ! Note that fco2 in is umolC/m2/sec so units need to be changed to gC/m2/sec
       do g = begg,endg
          nee_gcell(g) = nee_gcell(g)*12.011e-6_r8
       end do
#endif
       ! Convert from gC/m2/s to kgCO2/m2/s
       do g = begg,endg
          nee_gcell(g) = nee_gcell(g)*1.0e-3_r8*(44.0_r8/12.0_r8)
       end do
    end if

    if (icnt == 1) then          ! Initial call of averaging interval, copy data to accumulators
       do g = begg,endg
          taux_ave(g)  = taux_gcell(g)
          tauy_ave(g)  = tauy_gcell(g)
          lhflx_ave(g) = eflx_lh_tot_gcell(g) 
          shflx_ave(g) = eflx_sh_tot_gcell(g)
          lwup_ave(g)  = eflx_lwrad_out_gcell(g)
          qflx_ave(g)  = qflx_evap_tot_gcell(g)
          swabs_ave(g) = fsa_gcell(g)             
          if (index_l2c_Fall_nee /= 0) then
             nee_ave(g)  = nee_gcell(g)             
          end if
       end do
    else if (icnt == ncnt) then   ! Final call of averaging interval, complete averaging
       do g = begg,endg
          taux_ave(g)  = rncnt * (taux_ave(g)  + taux_gcell(g))            
          tauy_ave(g)  = rncnt * (tauy_ave(g)  + tauy_gcell(g))            
          lhflx_ave(g) = rncnt * (lhflx_ave(g) + eflx_lh_tot_gcell(g))
          shflx_ave(g) = rncnt * (shflx_ave(g) + eflx_sh_tot_gcell(g))     
          lwup_ave(g)  = rncnt * (lwup_ave(g)  + eflx_lwrad_out_gcell(g))  
          qflx_ave(g)  = rncnt * (qflx_ave(g)  + qflx_evap_tot_gcell(g))   
          swabs_ave(g) = rncnt * (swabs_ave(g) + fsa_gcell(g))             
          if (index_l2c_Fall_nee /= 0) then
             nee_ave(g)  = rncnt * (nee_ave(g)  + nee_gcell(g))             
          end if
       end do
    else                          ! Intermediate call, add data to accumulators
       do g = begg,endg
          taux_ave(g)  = (taux_ave(g)  + taux_gcell(g))            
          tauy_ave(g)  = (tauy_ave(g)  + tauy_gcell(g))            
          lhflx_ave(g) = (lhflx_ave(g) + eflx_lh_tot_gcell(g))
          shflx_ave(g) = (shflx_ave(g) + eflx_sh_tot_gcell(g))     
          lwup_ave(g)  = (lwup_ave(g)  + eflx_lwrad_out_gcell(g))  
          qflx_ave(g)  = (qflx_ave(g)  + qflx_evap_tot_gcell(g))   
          swabs_ave(g) = (swabs_ave(g) + fsa_gcell(g))             
          if (index_l2c_Fall_nee /= 0) then
             nee_ave(g)  = (nee_ave(g)  + nee_gcell(g))             
          end if
       end do
    end if

    ! Increment counter

    icnt = icnt + 1

  end subroutine csm_flxave
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: compat_check
!
! !INTERFACE:
  subroutine compat_check_spval(spval, data, string)
!
! !DESCRIPTION:
! Check that the given piece of real data sent from the coupler is valid
! data and not the couplers special data flag.  This ensures that the data
! you expect is actually being sent by the coupler.
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: spval
    real(r8), intent(in) :: data
    character(len=*), intent(in) :: string
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    if ( spval == data )then
       write(6,*)'ERROR:(compat_check_spval) msg incompatibility'
       write(6,*)'ERROR: I expect to recieve the data type: ',string
       write(6,*)'from CPL, but all I got was the special data flag'
       write(6,*)'coupler must not be sending this data, you are'
       write(6,*)'running with an incompatable version of the coupler'
       call endrun
    end if

  end subroutine compat_check_spval

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_compat
!
! !INTERFACE:
  subroutine csm_compat(cpl_maj_vers, cpl_min_vers, expect_maj_vers, &
                        expect_min_vers)
!
! !DESCRIPTION:
! Checks that the message recieved from the coupler is compatable
! with the type of message that I expect to recieve.  If the minor
! version numbers differ I print a warning message.  If the major
! numbers differ I abort since that means that the change is
! drastic enough that I can't run with the differences.
! Original Author: Erik Kluzek Dec/97
!
! !PARAMETERS:
    implicit none
    integer, intent(in) :: cpl_maj_vers    ! major version from coupler initial ibuffr array
    integer, intent(in) :: cpl_min_vers    ! minor version from coupler initial ibuffr array
    integer, intent(in) :: expect_maj_vers ! major version of the coupler I'm expecting
    integer, intent(in) :: expect_min_vers ! minor version of the coupler I'm expecting
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    write(6,*)'(cpl_COMPAT): This is revision: $Revision: 1.12.8.19.2.8 $'
    write(6,*)'              Tag: $Name: clm3_expa_48_brnchT_fmesh13 $'
    write(6,*)'              of the message compatability interface:'

    if ( cpl_min_vers /= expect_min_vers )then
       write(6,*) 'WARNING(cpl_compat):: Minor version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_min_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_min_vers
    end if

    if ( cpl_maj_vers /= expect_maj_vers )then
       write(6,*) 'ERROR(cpl_compat):: Major version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_maj_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_maj_vers
       call endrun
    end if

  end subroutine csm_compat

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_restart
!
! !INTERFACE:
  subroutine csm_restart(ncid, flag)
!
! !DESCRIPTION:
!  Read/write restart data needed for running in flux coupled mode
!
! !USES:
    use ncdio
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid           ! netcdf id
    character(len=*), intent(in) :: flag  ! 'read' or 'write'
!
! !REVISION HISTORY:
!  02.09.17  Mariana Vertenstein: moved code to be part of ccsm module
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: dosend_int   ! integer version of dosend 
    logical :: readvar      ! determine if variable is on initial file
!-----------------------------------------------------------------------

    if (flag == 'write') then
       if (dosend) then
          dosend_int = 1
       else
          dosend_int = 0 
       end if
    end if
       
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='ccsm_dosend_int', xtype=nf_int,  &
            long_name='ccsm send flag (integer)', units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_ioglobal(varname='ccsm_dosend_int', data=dosend_int, &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) then
             call endrun()
          end if
       end if
    end if
   
    if (flag == 'read') then
       if (dosend_int == 1) then
          dosend = .true.
       else
          dosend = .false.
       end if
    end if

  end subroutine csm_restart

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: global_sum_fld1d
!
! !INTERFACE:
  real(r8) function global_sum_fld1d(array,area)
!
! !DESCRIPTION:
! Performs a global sum on an input flux array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: array(:) !W/m2, Kg/m2-s or N/m2
    real(r8), intent(in) :: area(:)  !area
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g,i,j,n  ! indices
!------------------------------------------------------------------------

    ! Note: area is in km^2

    global_sum_fld1d = 0._r8
    do g = 1,numg
       global_sum_fld1d = global_sum_fld1d + array(g) * area(g) * 1.e6_r8
    end do

  end function global_sum_fld1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_restart
!
! !INTERFACE:
  logical function is_restart( )
!
! !DESCRIPTION:
! Determine if restart run
!
! !USES:
    use clm_varctl, only : nsrest
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    if (nsrest == 1) then
       is_restart = .true.
    else
       is_restart = .false.
    end if

  end function is_restart

#endif

end module clm_csmMod
