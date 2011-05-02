module clm_atmlnd

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_atmlnd
!
! !DESCRIPTION:
! Handle atm2lnd, lnd2atm mapping/downscaling/upscaling/data
!
! !USES:
  use clm_varpar  , only : numrad, ndst   !ndst = number of dust bins.
  use clm_varcon  , only : rair, grav, cpair, hfus, tfrz
  use clm_varctl  , only : iulog
  use decompMod   , only : get_proc_bounds, get_proc_bounds_atm
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
  use seq_drydep_mod, only : n_drydep, drydep_method, DD_XLND
  use clm_varpar  , only : nvoc
!
! !PUBLIC TYPES:
  implicit none
!----------------------------------------------------
! atmosphere -> land variables structure
!----------------------------------------------------
  type atm2lnd_type
     real(r8), pointer :: forc_t(:)       !atmospheric temperature (Kelvin)
     real(r8), pointer :: forc_u(:)       !atm wind speed, east direction (m/s)
     real(r8), pointer :: forc_v(:)       !atm wind speed, north direction (m/s)
     real(r8), pointer :: forc_wind(:)    !atmospheric wind speed   
     real(r8), pointer :: forc_q(:)       !atmospheric specific humidity (kg/kg)
     real(r8), pointer :: forc_hgt(:)     !atmospheric reference height (m)
     real(r8), pointer :: forc_hgt_u(:)   !obs height of wind [m] (new)
     real(r8), pointer :: forc_hgt_t(:)   !obs height of temperature [m] (new)
     real(r8), pointer :: forc_hgt_q(:)   !obs height of humidity [m] (new)
     real(r8), pointer :: forc_pbot(:)    !atmospheric pressure (Pa)
     real(r8), pointer :: forc_th(:)      !atm potential temperature (Kelvin)
     real(r8), pointer :: forc_vp(:)      !atmospheric vapor pressure (Pa) 
     real(r8), pointer :: forc_rho(:)     !density (kg/m**3)
     real(r8), pointer :: forc_rh(:)      !atmospheric relative humidity (%)
     real(r8), pointer :: forc_psrf(:)    !surface pressure (Pa)
     real(r8), pointer :: forc_pco2(:)    !CO2 partial pressure (Pa)
     real(r8), pointer :: forc_lwrad(:)   !downwrd IR longwave radiation (W/m**2)
     real(r8), pointer :: forc_solad(:,:) !direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai(:,:) !diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
     real(r8), pointer :: forc_solar(:)   !incident solar radiation
     real(r8), pointer :: forc_rain(:)    !rain rate [mm/s]
     real(r8), pointer :: forc_snow(:)    !snow rate [mm/s]
     real(r8), pointer :: forc_ndep(:)    !nitrogen deposition rate (gN/m2/s)
     real(r8), pointer :: rainf(:)        !ALMA rain+snow [mm/s]
#ifdef C13
     real(r8), pointer :: forc_pc13o2(:)  !C13O2 partial pressure (Pa)
#endif
     real(r8), pointer :: forc_po2(:)     !O2 partial pressure (Pa)
     real(r8), pointer :: forc_aer(:,:)   ! aerosol deposition array
  end type atm2lnd_type

!----------------------------------------------------
! land -> atmosphere variables structure
!----------------------------------------------------
  type lnd2atm_type
     real(r8), pointer :: t_rad(:)        !radiative temperature (Kelvin)
     real(r8), pointer :: t_ref2m(:)      !2m surface air temperature (Kelvin)
     real(r8), pointer :: q_ref2m(:)      !2m surface specific humidity (kg/kg)
     real(r8), pointer :: u_ref10m(:)     !10m surface wind speed (m/sec)
     real(r8), pointer :: h2osno(:)       !snow water (mm H2O)
     real(r8), pointer :: albd(:,:)       !(numrad) surface albedo (direct)
     real(r8), pointer :: albi(:,:)       !(numrad) surface albedo (diffuse)
     real(r8), pointer :: taux(:)         !wind stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy(:)         !wind stress: n-s (kg/m/s**2)
     real(r8), pointer :: eflx_lh_tot(:)  !total latent HF (W/m**2)  [+ to atm]
     real(r8), pointer :: eflx_sh_tot(:)  !total sensible HF (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lwrad_out(:) !IR (longwave) radiation (W/m**2)
     real(r8), pointer :: qflx_evap_tot(:)!qflx_evap_soi + qflx_evap_can + qflx_tran_veg
     real(r8), pointer :: fsa(:)          !solar rad absorbed (total) (W/m**2)
     real(r8), pointer :: nee(:)          !net CO2 flux (kg CO2/m**2/s) [+ to atm]
     real(r8), pointer :: ram1(:)         !aerodynamical resistance (s/m)
     real(r8), pointer :: fv(:)           !friction velocity (m/s) (for dust model)
     real(r8), pointer :: flxdst(:,:)       !dust flux (size bins)
     real(r8), pointer :: ddvel(:,:)        !dry deposition velocities
     real(r8), pointer :: flxvoc(:,:)       ! VOC flux (size bins)
  end type lnd2atm_type
  
  type(atm2lnd_type),public,target :: atm_a2l      ! a2l fields on atm grid
  type(lnd2atm_type),public,target :: atm_l2a      ! l2a fields on atm grid

  type(atm2lnd_type),public,target :: clm_a2l      ! a2l fields on clm grid
  type(lnd2atm_type),public,target :: clm_l2a      ! l2a fields on clm grid

  real(r8), pointer, public :: adiag_arain(:)
  real(r8), pointer, public :: adiag_asnow(:)
  real(r8), pointer, public :: adiag_aflux(:)
  real(r8), pointer, public :: adiag_lflux(:)

! !PUBLIC MEMBER FUNCTIONS:
  public :: init_adiag_type
  public :: init_atm2lnd_type
  public :: init_lnd2atm_type
  public :: clm_downscale_a2l
  public :: clm_map2gcell
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein and Tony Craig, 2006-01-10
!
! !PRIVATE MEMBER FUNCTIONS:

!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_adiag_type
!
! !INTERFACE:
  subroutine init_adiag_type
!
! !DESCRIPTION:
! Initialize downscaling diagnostics
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Created by T Craig for downscaling diagnostics, 4/2007
!
!
! !LOCAL VARIABLES:
!EOP
!
  integer :: beg,end
!------------------------------------------------------------------------

  call get_proc_bounds(beg, end)
  allocate(adiag_lflux(beg:end))
  adiag_lflux = 0.0_r8

  call get_proc_bounds_atm(beg, end)
  allocate(adiag_arain(beg:end))
  allocate(adiag_asnow(beg:end))
  allocate(adiag_aflux(beg:end))
  adiag_arain = 0.0_r8
  adiag_asnow = 0.0_r8
  adiag_aflux = 0.0_r8

end subroutine init_adiag_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_atm2lnd_type
!
! !INTERFACE:
  subroutine init_atm2lnd_type(beg, end, a2l)
!
! !DESCRIPTION:
! Initialize atmospheric variables required by the land
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (atm2lnd_type), intent(inout):: a2l
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Modified by T Craig, 11/01/05 for finemesh project
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(a2l%forc_t(beg:end))
  allocate(a2l%forc_u(beg:end))
  allocate(a2l%forc_v(beg:end))
  allocate(a2l%forc_wind(beg:end))
  allocate(a2l%forc_q(beg:end))
  allocate(a2l%forc_rh(beg:end))
  allocate(a2l%forc_hgt(beg:end))
  allocate(a2l%forc_hgt_u(beg:end))
  allocate(a2l%forc_hgt_t(beg:end))
  allocate(a2l%forc_hgt_q(beg:end))
  allocate(a2l%forc_pbot(beg:end))
  allocate(a2l%forc_th(beg:end))
  allocate(a2l%forc_vp(beg:end))
  allocate(a2l%forc_rho(beg:end))
  allocate(a2l%forc_psrf(beg:end))
  allocate(a2l%forc_pco2(beg:end))
  allocate(a2l%forc_lwrad(beg:end))
  allocate(a2l%forc_solad(beg:end,numrad))
  allocate(a2l%forc_solai(beg:end,numrad))
  allocate(a2l%forc_solar(beg:end))
  allocate(a2l%forc_rain(beg:end))
  allocate(a2l%forc_snow(beg:end))
  allocate(a2l%forc_ndep(beg:end))
  allocate(a2l%rainf(beg:end))
#if (defined C13)
  allocate(a2l%forc_pc13o2(beg:end))
#endif
  allocate(a2l%forc_po2(beg:end))
  allocate(a2l%forc_aer(beg:end,14))

  ! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  a2l%forc_t(beg:end) = ival
  a2l%forc_u(beg:end) = ival
  a2l%forc_v(beg:end) = ival
  a2l%forc_wind(beg:end) = ival
  a2l%forc_q(beg:end) = ival
  a2l%forc_rh(beg:end) = ival
  a2l%forc_hgt(beg:end) = ival
  a2l%forc_hgt_u(beg:end) = ival
  a2l%forc_hgt_t(beg:end) = ival
  a2l%forc_hgt_q(beg:end) = ival
  a2l%forc_pbot(beg:end) = ival
  a2l%forc_th(beg:end) = ival
  a2l%forc_vp(beg:end) = ival
  a2l%forc_rho(beg:end) = ival
  a2l%forc_psrf(beg:end) = ival
  a2l%forc_pco2(beg:end) = ival
  a2l%forc_lwrad(beg:end) = ival
  a2l%forc_solad(beg:end,1:numrad) = ival
  a2l%forc_solai(beg:end,1:numrad) = ival
  a2l%forc_solar(beg:end) = ival
  a2l%forc_rain(beg:end) = ival
  a2l%forc_snow(beg:end) = ival
  a2l%forc_ndep(beg:end) = ival
  a2l%rainf(beg:end) = nan
#ifdef C13
  a2l%forc_pc13o2(beg:end) = ival
#endif
  a2l%forc_po2(beg:end) = ival
  a2l%forc_aer(beg:end,:) = ival

end subroutine init_atm2lnd_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2atm_type
!
! !INTERFACE:
  subroutine init_lnd2atm_type(beg, end, l2a)
!
! !DESCRIPTION:
! Initialize land variables required by the atmosphere
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: beg, end
  type (lnd2atm_type), intent(inout):: l2a
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Modified by T Craig, 11/01/05 for finemesh project
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(l2a%t_rad(beg:end))
  allocate(l2a%t_ref2m(beg:end))
  allocate(l2a%q_ref2m(beg:end))
  allocate(l2a%u_ref10m(beg:end))
  allocate(l2a%h2osno(beg:end))
  allocate(l2a%albd(beg:end,1:numrad))
  allocate(l2a%albi(beg:end,1:numrad))
  allocate(l2a%taux(beg:end))
  allocate(l2a%tauy(beg:end))
  allocate(l2a%eflx_lwrad_out(beg:end))
  allocate(l2a%eflx_sh_tot(beg:end))
  allocate(l2a%eflx_lh_tot(beg:end))
  allocate(l2a%qflx_evap_tot(beg:end))
  allocate(l2a%fsa(beg:end))
  allocate(l2a%nee(beg:end))
  allocate(l2a%ram1(beg:end))
  allocate(l2a%fv(beg:end))
  allocate(l2a%flxdst(beg:end,1:ndst))
  allocate(l2a%flxvoc(beg:end,1:nvoc))
  if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
     allocate(l2a%ddvel(beg:end,1:n_drydep))
  end if

  ! ival = nan   ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  l2a%t_rad(beg:end) = ival
  l2a%t_ref2m(beg:end) = ival
  l2a%q_ref2m(beg:end) = ival
  l2a%u_ref10m(beg:end) = ival
  l2a%h2osno(beg:end) = ival
  l2a%albd(beg:end,1:numrad) = ival
  l2a%albi(beg:end,1:numrad) = ival
  l2a%taux(beg:end) = ival
  l2a%tauy(beg:end) = ival
  l2a%eflx_lwrad_out(beg:end) = ival
  l2a%eflx_sh_tot(beg:end) = ival
  l2a%eflx_lh_tot(beg:end) = ival
  l2a%qflx_evap_tot(beg:end) = ival
  l2a%fsa(beg:end) = ival
  l2a%nee(beg:end) = ival
  l2a%ram1(beg:end) = ival
  l2a%fv(beg:end) = ival
  l2a%flxdst(beg:end,1:ndst) = ival
  l2a%flxvoc(beg:end,1:nvoc) = ival
  if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
     l2a%ddvel(beg:end, : ) = ival
  end if

end subroutine init_lnd2atm_type

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_downscale_a2l
!
! !INTERFACE:
subroutine clm_downscale_a2l(a2l_src, a2l_dst)
!
! !DESCRIPTION:
! Maps atm2lnd fields from external grid to clm grid
!
! !USES:
  use downscaleMod, only : map1dl_a2l, map1dl_l2a, map_setptrs
  use decompMod   , only : ldecomp,adecomp
  use domainMod   , only : ldomain,adomain
  use QSatMod     , only : QSat

!
! !ARGUMENTS:
  implicit none
  type(atm2lnd_type), intent(in)  :: a2l_src
  type(atm2lnd_type), intent(out) :: a2l_dst
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: n                     ! loop counter
  integer :: ix                    ! field index
  integer :: nflds                 ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  integer          :: nmap         ! size of map
  integer          :: mo           ! size of map
  integer, pointer :: src(:)       ! map src index
  integer, pointer :: dst(:)       ! map dst index
  real(r8),pointer :: wts(:)       ! map wts values
  integer  :: ns                    !source (atm) indexes
  integer  :: nd                    !destination (lnd) indexes
  ! temporaries for topo downscaling:
  real(r8) :: hsurf_a,hsurf_l,Hbot,Hsrf,lapse
  real(r8) :: zbot_a, tbot_a, pbot_a, thbot_a, qbot_a, qs_a, es_a
  real(r8) :: zbot_l, tbot_l, pbot_l, thbot_l, qbot_l, qs_l, es_l
  real(r8) :: tsrf_l, psrf_l, egcm_l, rhos_l
  real(r8) :: dum1,dum2,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8
  real(r8) :: mu
  real(r8) :: precip_max
  real(r8) :: rand_num(6)             ! random number
  real(r8),allocatable :: rain_l(:),snow_l(:)
  real(r8),allocatable :: pnorm(:)
  real(r8),allocatable :: pchance(:),pamount(:)
  integer ,allocatable :: pchance_min(:)
  real(r8),allocatable :: qsum(:)
  integer  :: ier
  logical  :: first_call = .true.

!------------------------------------------------------------------------------

  if (first_call .and. masterproc) then
    write(iulog,*) 'clm_downscale_a2l subroutine'
  endif

  nradflds = size(a2l_src%forc_solad,dim=2)
  if (nradflds /= numrad) then
    write(iulog,*) 'clm_mapa2l ERROR: nradflds ne numrad ',nradflds,numrad
    call endrun()
  endif

  call get_proc_bounds_atm(begg_s, endg_s)
  call get_proc_bounds    (begg_d, endg_d)

  !-topographic downscaling
  !-only call this if there is more than 1 land cell / atm cell somewhere

  call map_setptrs(map1dl_l2a, dstmo=mo)

  if (mo > 1) then

     if (first_call.and.masterproc) then
        write(iulog,*) 'clm_mapa2l downscaling ON'
     endif

     call map_setptrs(map1dl_a2l,nwts=nmap,src=src,dst=dst,dstmo=mo)
     if (mo /= 1) then
        write(iulog,*)' clm_mapa2l ERROR: map1dl_a2l mo not 1 ',mo
        call endrun()
     endif

     lapse   = 0.0065_r8     ! hardwired in multiple places in cam
     mu = 0.20               ! chance of precip coverage [0.,1.]
     precip_max = 10.0_r8    ! mm/hr
     allocate(pchance(begg_d:endg_d),pamount(begg_d:endg_d))

     do n = 1,nmap
        ns = src(n)
        nd = dst(n)
        
        hsurf_a = adomain%topo(ns)        ! atm elevation
        hsurf_l = ldomain%ntop(nd)        ! lnd elevation
        !tcx DOWNSCALING turns off topo downscaling
        !    hsurf_l = hsurf_a

        if (abs(hsurf_a - hsurf_l) .gt. 0.1_r8) then

           !tcx DOWNSCALING if atm lapse rate is available
           !       lapse  = a2l_src%lapse(ns)         ! atm lapse rate
           tbot_a = a2l_src%forc_t(ns)        ! atm temp at bot
           thbot_a= a2l_src%forc_th(ns)       ! atm pot temp at bot
           pbot_a = a2l_src%forc_pbot(ns)     ! atm press at bot
           qbot_a = a2l_src%forc_q(ns)        ! atm sp humidity at bot
           zbot_a = a2l_src%forc_hgt(ns)      ! atm ref height
           
           zbot_l = zbot_a
           tbot_l = tbot_a-lapse*(hsurf_l-hsurf_a)          ! lnd temp for topo
           
           Hbot   = rair*0.5_r8*(tbot_a+tbot_l)/grav        ! scale ht at avg temp
           pbot_l = pbot_a*exp(-(hsurf_l-hsurf_a)/Hbot)     ! lnd press for topo
           thbot_l= tbot_l*exp((zbot_l/Hbot)*(rair/cpair))  ! pot temp calc
           
           tsrf_l = tbot_l-lapse*(-zbot_l)                  ! lnd temp at surface
           Hsrf   = rair*0.5_r8*(tbot_l+tsrf_l)/grav        ! scale ht at avg temp
           psrf_l = pbot_l*exp(-(zbot_l)/Hsrf)              ! lnd press for topo
           
           call Qsat(tbot_a,pbot_a,es_a,dum1,qs_a,dum2)
           call Qsat(tbot_l,pbot_l,es_l,dum1,qs_l,dum2)
           qbot_l = qbot_a*(qs_l/qs_a)
           
           a2l_dst%forc_hgt(nd)  = zbot_l
           a2l_dst%forc_t(nd)    = tbot_l
           a2l_dst%forc_pbot(nd) = pbot_l
           a2l_dst%forc_th(nd)   = thbot_l
           a2l_dst%forc_q(nd)    = qbot_l
           a2l_dst%forc_vp(nd)   = es_l
           a2l_dst%forc_psrf(nd) = psrf_l
           
        endif
        
        ! DOWNSCALING set random numbers (must be "reproducible") for downscaling
        
        rand_num(1) = mod(abs(a2l_dst%forc_t(nd))            *1.0e3_r8,1.0_r8)
        rand_num(2) = 0.5_r8
        rand_num(3) = mod(abs(sin(ldecomp%gdc2glo(nd)*1._r8))*1.0e3_r8,1.0_r8)
        rand_num(4) = mod(abs(a2l_dst%forc_psrf(nd))         *1.0e1_r8,1.0_r8)
        rand_num(5) = 0.5_r8
        rand_num(6) = mod(abs(cos(ldecomp%gdc2glo(nd)*1._r8))*1.0e3_r8,1.0_r8)
        
        pchance(nd) = mod(rand_num(1)+rand_num(2)+rand_num(3),1.0_r8)
        pamount(nd) = -log(max(1.0_r8-(mod(rand_num(4)+rand_num(5)+rand_num(6),1.0_r8)),1.0e-20_r8))
        
        !tcx DOWNSCALING this will remove random aspects, turn on either or both
        !    pchance(nd) = 0._r8    ! sets all fine gridcells to precip true
        !    pamount(nd) = 1._r8    ! sets all precip gridcells to equal amount
        
     enddo

     !tcx diagnostics
     !  write(iulog,*) 'tcx mm1 ',minval(a2l_dst%forc_t),maxval(a2l_dst%forc_t)
     !  write(iulog,*) 'tcx mm2 ',minval(a2l_dst%forc_psrf),maxval(a2l_dst%forc_psrf)
     !  write(iulog,*) 'tcx mm3 ',minval(a2l_dst%forc_q),maxval(a2l_dst%forc_q)
     !  write(iulog,*) 'tcx mm4 ',minval(a2l_dst%forc_th),maxval(a2l_dst%forc_th)
     !  write(iulog,*) 'tcx mm5 ',minval(pchance),maxval(pchance)
     !  write(iulog,*) 'tcx mm6 ',minval(pamount),maxval(pamount)
     
     call map_setptrs(map1dl_l2a,nwts=nmap,src=src,dst=dst,wts=wts)
     
     ! compute precipitation disaggregation
     !--- find min pchance in finemesh cells of coarse mesh, set that cell's
     !--- pchance to 0. to force it to be below mu and take precip
     allocate(pchance_min(begg_s:endg_s))
     pchance_min = -999
     do n = 1,nmap
        ns = dst(n)
        nd = src(n)
        if (pchance_min(ns) < 1) pchance_min(ns) = nd
        if (pchance(nd) < pchance(pchance_min(ns))) pchance_min(ns) = nd
     enddo
     do n = begg_s,endg_s
        if (pchance_min(n) > 0) pchance(pchance_min(n)) = 0._r8
     enddo
     deallocate(pchance_min)
     
     !--- set rain/snow amounts and compute sums for normalization
     allocate(rain_l(begg_s:endg_s),snow_l(begg_s:endg_s))
     rain_l = 0.0_r8
     snow_l = 0.0_r8
     do n = 1,nmap
        ns = dst(n)
        nd = src(n)
        !tcx DOWNSCALING turn off precip disaggr completely
        if (pchance(nd) < mu) then
           !--- rain/snow refractionation
           !--- set to 100% snow < -5C and 100% rain > 0C
           dum2 = (a2l_dst%forc_t(nd)-tfrz + 5.0_r8)/(5.0_r8)
           dum2 = max(dum2,0.0_r8)
           dum2 = min(dum2,1.0_r8)
           
           !tcx DOWNSCALING turn on refractionation and spatial variability
           a2l_dst%forc_rain(nd) = (         dum2) * pamount(nd)
           a2l_dst%forc_snow(nd) = (1.0_r8 - dum2) * pamount(nd)
           !tcx DOWNSCALING turn on just spatial variability but not refractionation
           ! comment out lines above and comment in the two lines below
           !        a2l_dst%forc_rain(nd) = a2l_dst%forc_rain(nd)* pamount(nd)
           !        a2l_dst%forc_snow(nd) = a2l_dst%forc_snow(nd)* pamount(nd)
           !
        else
           a2l_dst%forc_rain(nd) = 0._r8
           a2l_dst%forc_snow(nd) = 0._r8
        endif
        rain_l(ns) = rain_l(ns) + a2l_dst%forc_rain(nd) * wts(n)
        snow_l(ns) = snow_l(ns) + a2l_dst%forc_snow(nd) * wts(n)
     enddo
     deallocate(pchance,pamount)
     
     !--- compute normalization of disaggregation amounts
     allocate(pnorm(begg_s:endg_s))
     pnorm = 0.0_r8
     do ns = begg_s, endg_s
        if (rain_l(ns) + snow_l(ns) == 0._r8) then
           if (a2l_src%forc_rain(ns) + a2l_src%forc_snow(ns) /= 0._r8) then
              write(iulog,*)' clm_mapa2l ERROR: rain/snow normalization',rain_l(ns),snow_l(ns), &
                   a2l_src%forc_rain(ns),a2l_src%forc_snow(ns)
              call endrun()
           endif
        else
           pnorm(ns) = (a2l_src%forc_rain(ns) + a2l_src%forc_snow(ns)) / (rain_l(ns)+snow_l(ns))
        endif
        !--- set pnorm=1 if close, for bfb consitency when coarse=finemesh cases
        if (abs(pnorm(ns) - 1.0_r8) < 1.0e-12) pnorm(ns) = 1.0_r8
     enddo
     deallocate(rain_l,snow_l)
     
     !--- apply normalization
     do n = 1,nmap
        ns = dst(n)
        nd = src(n)
        a2l_dst%forc_rain(nd) = a2l_dst%forc_rain(nd) * pnorm(ns)
        a2l_dst%forc_snow(nd) = a2l_dst%forc_snow(nd) * pnorm(ns)
     enddo
     deallocate(pnorm)
     
     !--- compute q normalization amounts
     call map_setptrs(map1dl_l2a,nwts=nmap,src=src,dst=dst,wts=wts)
     allocate(qsum(begg_s:endg_s))
     qsum = 0.0_r8
     do n = 1,nmap
        ns = dst(n)
        nd = src(n)
        qsum(ns) = qsum(ns) + wts(n)* a2l_dst%forc_q(nd)
     enddo
     
     call map_setptrs(map1dl_a2l,nwts=nmap,src=src,dst=dst)
     do n = 1,nmap
        ns = src(n)
        nd = dst(n)
        
        qbot_a = a2l_src%forc_q(ns)        ! atm specific humidity
        qbot_l = a2l_dst%forc_q(nd)        ! lnd specific humidity
        pbot_l = a2l_dst%forc_pbot(nd)  
        tbot_l = a2l_dst%forc_t(nd)  
        
        qbot_l = qbot_l - (qsum(ns) - qbot_a)        ! normalize
        egcm_l = qbot_l*pbot_l/(0.622+0.378*qbot_l)
        rhos_l = (pbot_l-0.378*egcm_l) / (rair*tbot_l)
        
        a2l_dst%forc_q(nd)    = qbot_l
        a2l_dst%forc_rho(nd)  = rhos_l
        
     enddo
     
     deallocate(qsum)
     
     ! --- check ---
     call map_setptrs(map1dl_l2a,nwts=nmap,src=src,dst=dst,wts=wts)
     adiag_arain = 0.0_r8
     adiag_asnow = 0.0_r8
     adiag_aflux = 0.0_r8
     adiag_lflux = 0.0_r8
     do ns = begg_s,endg_s
        sum1 = 0.0_r8
        sum2 = 0.0_r8
        sum3 = 0.0_r8
        sum4 = 0.0_r8
        sum5 = 0.0_r8
        sum6 = 0.0_r8
        sum7 = 0.0_r8
        sum8 = 0.0_r8
        do n = 1,nmap
           if (dst(n) == ns) then
              nd = src(n)
              sum1 = sum1 + ldomain%ntop(nd)      * wts(n)
              sum2 = sum2 + a2l_dst%forc_t(nd)    * wts(n)
              sum3 = sum3 + a2l_dst%forc_q(nd)    * wts(n)
              sum4 = sum4 + a2l_dst%forc_hgt(nd)  * wts(n)
              sum5 = sum5 + a2l_dst%forc_pbot(nd) * wts(n)
              sum6 = sum6 + a2l_dst%forc_th(nd)   * wts(n)
              sum7 = sum7 + a2l_dst%forc_rain(nd) * wts(n)
              sum8 = sum8 + a2l_dst%forc_snow(nd) * wts(n)
              adiag_lflux(nd) = (a2l_dst%forc_snow(nd) - a2l_src%forc_snow(ns))*wts(n)*hfus
           endif
        enddo
        adiag_arain(ns) = sum7
        adiag_asnow(ns) = sum8
        adiag_aflux(ns) = (sum8 - a2l_src%forc_snow(ns))*hfus
        
        !--- add up rain and snow, compute relative diff vs source
        dum1 = (sum7 + sum8 - a2l_src%forc_rain(ns) - a2l_src%forc_snow(ns))
        if ((sum7+sum8) == 0.0_r8 .and. &
             (a2l_src%forc_rain(ns)+a2l_src%forc_snow(ns)) == 0.0_r8) then
           dum1 = 0.0_r8
        else
           dum1 = dum1 / (max(sum7+sum8,a2l_src%forc_rain(ns)+a2l_src%forc_snow(ns)))
        endif
        
        if   ((abs(sum1 - adomain%topo(ns))      > 1.0e-8) &
             .or.(abs(sum2 - a2l_src%forc_t(ns))    > 1.0e-3) &
             .or.(abs(sum3 - a2l_src%forc_q(ns))    > 1.0e-8) &
             .or.(abs(sum4 - a2l_src%forc_hgt(ns))  > 1.0e-6) &
             !      .or.(abs(sum5 - a2l_src%forc_pbot(ns)) > 1.0e-6) &
             !      .or.(abs(sum6 - a2l_src%forc_th(ns))   > 1.0e-6) &
             .or.(abs(dum1)                         > 1.0e-8) &
             ) then
           write(iulog,*) 'clm_map2l check ERROR topo ',ns,sum1,adomain%topo(ns)
           write(iulog,*) 'clm_map2l check ERROR t    ',ns,sum2,a2l_src%forc_t(ns)
           write(iulog,*) 'clm_map2l check ERROR q    ',ns,sum3,a2l_src%forc_q(ns)
           write(iulog,*) 'clm_map2l check ERROR hgt  ',ns,sum4,a2l_src%forc_hgt(ns)
           write(iulog,*) 'clm_map2l check ERROR pbot ',ns,sum5,a2l_src%forc_pbot(ns)
           write(iulog,*) 'clm_map2l check ERROR th   ',ns,sum6,a2l_src%forc_th(ns)
           write(iulog,*) 'clm_map2l check ERROR rain ',ns,sum7,a2l_src%forc_rain(ns)
           write(iulog,*) 'clm_map2l check ERROR snow ',ns,sum8,a2l_src%forc_snow(ns)
           write(iulog,*) 'clm_map2l check ERROR adiag',ns,adiag_arain(ns),adiag_asnow(ns)
           call endrun()
        endif
     enddo
     
  else    ! mx_ovr > 1
     
     write(iulog,*)' need mx_ovr > 1 for downscaling'
     call endrun()

  endif   ! mx_ovr > 1

  first_call = .false.

end subroutine clm_downscale_a2l

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_map2gcell
!
! !INTERFACE: subroutine clm_map2gcell(init)
subroutine clm_map2gcell(init)
!
! !DESCRIPTION:
! Compute l2a component of gridcell derived type
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use subgridAveMod
  use clm_varcon  , only : sb
  use clm_varpar  , only : numrad
!
! !ARGUMENTS:
  implicit none
  save
  logical, optional, intent(in) :: init  ! if true=>only set a subset of arguments
!
! !REVISION HISTORY:
! Mariana Vertenstein: created 03/10-25
!
!
! !LOCAL VARIABLES:
  integer :: begp, endp      ! per-proc beginning and ending pft indices
  integer :: begc, endc      ! per-proc beginning and ending column indices
  integer :: begl, endl      ! per-proc beginning and ending landunit indices
  integer :: begg, endg      ! per-proc gridcell ending gridcell indices
!
! !USES:
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
! 03-08-25 : Updated to vector data structure (Mariana Vertenstein)
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: g                           ! indices
  type(gridcell_type), pointer :: gptr   ! pointer to gridcell derived subtype
  type(landunit_type), pointer :: lptr   ! pointer to landunit derived subtype
  type(column_type)  , pointer :: cptr   ! pointer to column derived subtype
  type(pft_type)     , pointer :: pptr   ! pointer to pft derived subtype
  type(pft_dflux_type),pointer :: pdf    ! local pointer to derived subtype
  integer             :: n               ! Loop index over nmap
  real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
  real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
  real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
  ! The following converts g of C to kg of CO2
  real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)

!------------------------------------------------------------------------

  ! Set pointers into derived type

  gptr => clm3%g
  lptr => clm3%g%l
  cptr => clm3%g%l%c
  pptr => clm3%g%l%c%p

  ! Determine processor bounds

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

  ! Compute gridcell averages. 

  if (present(init)) then

     call c2g(begc, endc, begl, endl, begg, endg, &
          cptr%cws%h2osno, clm_l2a%h2osno, &
          c2l_scale_type= 'urbanf', l2g_scale_type='unity')
     do g = begg,endg
        clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
     end do
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albd, clm_l2a%albd,&
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albi, clm_l2a%albi,&
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out,&
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
     do g = begg,endg
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
     end do

  else

     call c2g(begc, endc, begl, endl, begg, endg, cptr%cws%h2osno, clm_l2a%h2osno,&
          c2l_scale_type= 'urbanf', l2g_scale_type='unity')
     do g = begg,endg
        clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
     end do

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albd, clm_l2a%albd, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albi, clm_l2a%albi, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pes%t_ref2m, clm_l2a%t_ref2m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pes%q_ref2m, clm_l2a%q_ref2m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pps%u10_clm, clm_l2a%u_ref10m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pmf%taux, clm_l2a%taux, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pmf%tauy, clm_l2a%tauy, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lh_tot, clm_l2a%eflx_lh_tot, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

!DML note: use new array: clm3%g%gef%eflx_sh_totg
!     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
!          pptr%pef%eflx_sh_tot, clm_l2a%eflx_sh_tot, &
!          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     do g = begg,endg
        clm_l2a%eflx_sh_tot(g) = clm3%g%gef%eflx_sh_totg(g)
     end do

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pwf%qflx_evap_tot, clm_l2a%qflx_evap_tot, & 
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%fsa, clm_l2a%fsa, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
                  
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out, &
          p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')
                  
#if (defined CN)
     call c2g(begc, endc, begl, endl, begg, endg, &
          cptr%ccf%nee, clm_l2a%nee, &
          c2l_scale_type= 'unity', l2g_scale_type='unity')
#elif (defined CASA)
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pps%co2flux, clm_l2a%nee, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
#else
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pcf%fco2, clm_l2a%nee, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
     ! Note that fco2 in is umolC/m2/sec so units need to be changed to gC/m2/sec
     do g = begg,endg
        clm_l2a%nee(g) = clm_l2a%nee(g)*12.011e-6_r8
     end do
#endif

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pps%fv, clm_l2a%fv, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pps%ram1, clm_l2a%ram1, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, ndst, &
           pptr%pdf%flx_mss_vrt_dst, clm_l2a%flxdst, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, nvoc, &
           pptr%pvf%vocflx, clm_l2a%flxvoc, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      if ( n_drydep > 0 .and. drydep_method == DD_XLND ) &
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, n_drydep, &
           pptr%pdd%drydepvel, clm_l2a%ddvel, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     ! Convert from gC/m2/s to kgCO2/m2/s
     do g = begg,endg
        clm_l2a%nee(g) = clm_l2a%nee(g)*convertgC2kgCO2
     end do

     do g = begg,endg
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
     end do

     ! DOWNSCALING add rain snow conversion heat flux to latent heat flux
     do g = begg,endg
        clm_l2a%eflx_lh_tot(g) = clm_l2a%eflx_lh_tot(g) + adiag_lflux(g)
     end do

  end if

end subroutine clm_map2gcell

!------------------------------------------------------------------------
!------------------------------------------------------------------------
end module clm_atmlnd

