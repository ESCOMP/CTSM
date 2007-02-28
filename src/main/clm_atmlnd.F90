#include <misc.h>
#include <preproc.h>

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
                                          !only used # ifdef DUST
  use clm_varcon  , only : rair, grav, cpair
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
  use spmdMod     , only : masterproc
  use abortutils,   only : endrun
!
! !PUBLIC TYPES:
  implicit none
!----------------------------------------------------
! atmosphere -> land variables structure
!----------------------------------------------------
type atm2lnd_type
#if (defined OFFLINE)
  real(r8), pointer :: flfall(:)       !frac of liquid water in falling precip
#endif
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
  real(r8), pointer :: forc_psrf(:)    !surface pressure (Pa)
  real(r8), pointer :: forc_pco2(:)    !CO2 partial pressure (Pa)
  real(r8), pointer :: forc_lwrad(:)   !downwrd IR longwave radiation (W/m**2)
  real(r8), pointer :: forc_solad(:,:) !direct beam radiation (numrad)
                                       !(vis=forc_sols , nir=forc_soll )
  real(r8), pointer :: forc_solai(:,:) !diffuse radiation (numrad)
                                       !(vis=forc_solsd, nir=forc_solld)
  real(r8), pointer :: forc_solar(:)   !incident solar radiation
  real(r8), pointer :: forc_rain(:)    !rain rate [mm/s]
  real(r8), pointer :: forc_snow(:)    !snow rate [mm/s]
  real(r8), pointer :: forc_ndep(:)    !nitrogen deposition rate (gN/m2/s)
  ! 4/14/05: PET
  ! Adding isotope code
  real(r8), pointer :: forc_pc13o2(:)  !C13O2 partial pressure (Pa)
  real(r8), pointer :: forc_po2(:)     !O2 partial pressure (Pa)
end type atm2lnd_type

!----------------------------------------------------
! land -> atmosphere variables structure
!----------------------------------------------------
type lnd2atm_type
  real(r8), pointer :: t_rad(:)        !radiative temperature (Kelvin)
  real(r8), pointer :: t_ref2m(:)      !2m surface air temperature (Kelvin)
  real(r8), pointer :: q_ref2m(:)      !2m surface specific humidity (kg/kg)
  real(r8), pointer :: h2osno(:)       !snow water (mm H2O)
  real(r8), pointer :: albd(:,:)       !(numrad) surface albedo (direct)
  real(r8), pointer :: albi(:,:)       !(numrad) surface albedo (diffuse)
  real(r8), pointer :: taux(:)         !wind stress: e-w (kg/m/s**2)
  real(r8), pointer :: tauy(:)         !wind stress: n-s (kg/m/s**2)
  real(r8), pointer :: eflx_lh_tot(:)  !total latent HF (W/m**2)  [+ to atm]
  real(r8), pointer :: eflx_sh_tot(:)  !total sensible HF (W/m**2) [+ to atm]
  real(r8), pointer :: eflx_lwrad_out(:) !IR (longwave) radiation (W/m**2)
  real(r8), pointer :: qflx_evap_tot(:)!qflx_evap(_soi + _veg) + qflx_tran_veg
  real(r8), pointer :: fsa(:)          !solar rad absorbed (total) (W/m**2)
  real(r8), pointer :: nee(:)          !net CO2 flux (kg C/m**2/s) [+ to atm]
#if (defined DUST || defined  PROGSSLT )
  real(r8), pointer :: ram1(:)         !aerodynamical resistance (s/m)
  real(r8), pointer :: fv(:)           !friction velocity (m/s) (for dust model)
#endif
#if (defined DUST  )
  real(r8), pointer :: flxdst(:,:)       !dust flux (size bins)
#endif
end type lnd2atm_type

  type(atm2lnd_type),public,target :: atm_a2l      ! a2l fields on atm grid
  type(lnd2atm_type),public,target :: atm_l2a      ! l2a fields on atm grid

  type(atm2lnd_type),public,target :: clm_a2l      ! a2l fields on clm grid
  type(lnd2atm_type),public,target :: clm_l2a      ! l2a fields on clm grid

! !PUBLIC MEMBER FUNCTIONS:
  public :: init_atm2lnd_type
  public :: init_lnd2atm_type
  public :: clm_mapa2l
  public :: clm_mapl2a
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
!EOP
!
! !LOCAL VARIABLES:
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

#if (defined OFFLINE)
  allocate(a2l%flfall(beg:end))
#endif
  allocate(a2l%forc_t(beg:end))
  allocate(a2l%forc_u(beg:end))
  allocate(a2l%forc_v(beg:end))
  allocate(a2l%forc_wind(beg:end))
  allocate(a2l%forc_q(beg:end))
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
  ! 4/14/05: PET
  ! Adding isotope code
  allocate(a2l%forc_pc13o2(beg:end))
  allocate(a2l%forc_po2(beg:end))

! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

#if (defined OFFLINE)
  a2l%flfall(beg:end) = ival
#endif
  a2l%forc_t(beg:end) = ival
  a2l%forc_u(beg:end) = ival
  a2l%forc_v(beg:end) = ival
  a2l%forc_wind(beg:end) = ival
  a2l%forc_q(beg:end) = ival
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
  ! 4/14/05: PET
  ! Adding isotope code
  a2l%forc_pc13o2(beg:end) = ival
  a2l%forc_po2(beg:end) = ival

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
!EOP
!
! !LOCAL VARIABLES:
  real(r8) :: ival   ! initial value
!------------------------------------------------------------------------

  allocate(l2a%t_rad(beg:end))
  allocate(l2a%t_ref2m(beg:end))
  allocate(l2a%q_ref2m(beg:end))
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
#if (defined DUST || defined  PROGSSLT )
  allocate(l2a%ram1(beg:end))
  allocate(l2a%fv(beg:end))
#endif
#if (defined DUST )
  allocate(l2a%flxdst(beg:end,1:ndst))
#endif

! ival = nan      ! causes core dump in map_maparray, tcx fix
  ival = 0.0_r8

  l2a%t_rad(beg:end) = ival
  l2a%t_ref2m(beg:end) = ival
  l2a%q_ref2m(beg:end) = ival
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
#if (defined DUST || defined  PROGSSLT )
  l2a%ram1(beg:end) = ival
  l2a%fv(beg:end) = ival
#endif
#if (defined DUST )
  l2a%flxdst(beg:end,1:ndst) = ival
#endif
end subroutine init_lnd2atm_type

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_mapa2l
!
! !INTERFACE:
  subroutine clm_mapa2l(a2l_src, a2l_dst)
!
! !DESCRIPTION:
! Maps atm2lnd fields from external grid to clm grid
!
! !USES:
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
  use areaMod  , only : map_maparrayl, map1dl_a2l, map1dl_l2a, map_setptrs
  use decompMod, only : ldecomp,adecomp
  use domainMod, only : ldomain,adomain
  use QSatMod,   only : QSat
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
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: ix                    ! field index
  integer :: nflds                 ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
  integer          :: nmap         ! size of map
  integer          :: mo           ! size of map
  integer, pointer :: src(:)       ! map src index
  integer, pointer :: dst(:)       ! map dst index
  real(r8),pointer :: wts(:)       ! map wts values
  integer :: ns                    !source (atm) indexes
  integer :: nd                    !destination (lnd) indexes
  ! temporaries for topo downscaling:
  real(r8):: hsurf_a,hsurf_l,Hbot,Hsrf,lapse
  real(r8):: zbot_a, tbot_a, pbot_a, thbot_a, qbot_a, qs_a, es_a
  real(r8):: zbot_l, tbot_l, pbot_l, thbot_l, qbot_l, qs_l, es_l
  real(r8):: tsrf_l, psrf_l, egcm_l, rhos_l
  real(r8):: dum1,dum2,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8
  real(r8),allocatable :: qsum(:)
  logical :: first_call = .true.
!------------------------------------------------------------------------------

  if (first_call .and. masterproc) then
    write(6,*) 'clm_mapa2l subroutine'
  endif

  nradflds = size(a2l_src%forc_solad,dim=2)
  if (nradflds /= numrad) then
    write(6,*) 'clm_mapa2l ERROR: nradflds ne numrad ',nradflds,numrad
    call endrun()
  endif

  !--- allocate temporaries
  call get_proc_bounds_atm(begg_s, endg_s)
  call get_proc_bounds    (begg_d, endg_d)

  nflds = 21+2*numrad

  allocate(asrc(begg_s:endg_s,nflds))
  allocate(adst(begg_d:endg_d,nflds))

  ix = 0
  ix=ix+1; asrc(:,ix) = a2l_src%forc_t(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_u(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_v(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_wind(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_q(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_u(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_t(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_hgt_q(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_pbot(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_th(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_vp(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_rho(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_psrf(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_pco2(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_lwrad(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_solar(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_rain(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_snow(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_pc13o2(:)  
  ix=ix+1; asrc(:,ix) = a2l_src%forc_po2(:)  
  do n = 1,numrad
     ix=ix+1; asrc(:,ix) = a2l_src%forc_solad(:,n)  
     ix=ix+1; asrc(:,ix) = a2l_src%forc_solai(:,n)  
  enddo
!-forc_ndep is not recd from atm,don't know why it's in a2l (TCFIX) ---
!-forc_ndep cannot be updated here, array will be trashed and CN will fail ---
!  asrc(:,xx) = a2l_src%forc_ndep(:)  

#if (defined OFFLINE)
  call map_maparrayl(begg_s, endg_s, begg_d, endg_d, 1, a2l_src%flfall, a2l_dst%flfall, map1dl_a2l)
#endif
  call map_maparrayl(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, map1dl_a2l)

  ix = 0
  ix=ix+1; a2l_dst%forc_t(:)     =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_u(:)     =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_v(:)     =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_wind(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_q(:)     =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_hgt(:)   =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_hgt_u(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_hgt_t(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_hgt_q(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_pbot(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_th(:)    =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_vp(:)    =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_rho(:)   =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_psrf(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_pco2(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_lwrad(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_solar(:) =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_rain(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_snow(:)  =   adst(:,ix)
  ix=ix+1; a2l_dst%forc_pc13o2(:)=   adst(:,ix)
  ix=ix+1; a2l_dst%forc_po2(:)   =   adst(:,ix)
  do n = 1,numrad
     ix=ix+1; a2l_dst%forc_solad(:,n)  = adst(:,ix)
     ix=ix+1; a2l_dst%forc_solai(:,n)  = adst(:,ix)
  enddo

  deallocate(asrc)
  deallocate(adst)

  if (first_call.and.masterproc) then
    write(6,*) 'clm_mapa2l mapping complete'
  endif

!-topographic downscaling
!-only call this if there is more than 1 land cell / atm cell somewhere
  call map_setptrs(map1dl_l2a,dstmo=mo)
  if (mo > 1) then

  if (first_call.and.masterproc) then
    write(6,*) 'clm_mapa2l downscaling ON'
  endif

  call map_setptrs(map1dl_a2l,nwts=nmap,src=src,dst=dst,dstmo=mo)
  if (mo /= 1) then
     write(6,*)' clm_mapa2l ERROR: map1dl_a2l mo not 1 ',mo
     call endrun()
  endif

  lapse   = 0.0065_r8                  ! hardwired in multiple places in cam

  do n = 1,nmap
    ns = src(n)
    nd = dst(n)

    hsurf_a = adomain%topo(ns)        ! atm elevation
    hsurf_l = ldomain%ntop(nd)        ! lnd elevation

    if (abs(hsurf_a - hsurf_l) .gt. 0.1_r8) then

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
  enddo

  allocate(qsum(begg_s:endg_s))
  qsum = 0.0_r8
  call map_setptrs(map1dl_l2a,nwts=nmap,src=src,dst=dst,wts=wts)
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
  do ns = begg_s,endg_s
    sum1 = 0.0_r8
    sum2 = 0.0_r8
    sum3 = 0.0_r8
    sum4 = 0.0_r8
    sum5 = 0.0_r8
    sum6 = 0.0_r8
    do n = 1,nmap
      if (dst(n) == ns) then
        nd = src(n)
        sum1 = sum1 + ldomain%ntop(nd)   * wts(n)
        sum2 = sum2 + a2l_dst%forc_t(nd)    * wts(n)
        sum3 = sum3 + a2l_dst%forc_q(nd)    * wts(n)
        sum4 = sum4 + a2l_dst%forc_hgt(nd)  * wts(n)
        sum5 = sum5 + a2l_dst%forc_pbot(nd) * wts(n)
        sum6 = sum6 + a2l_dst%forc_th(nd)   * wts(n)
      endif
    enddo
    if   ((abs(sum1 - adomain%topo(ns))   > 1.0e-8) &
      .or.(abs(sum2 - a2l_src%forc_t(ns))    > 1.0e-3) &
      .or.(abs(sum3 - a2l_src%forc_q(ns))    > 1.0e-8) &
      .or.(abs(sum4 - a2l_src%forc_hgt(ns))  > 1.0e-6) &
!      .or.(abs(sum5 - a2l_src%forc_pbot(ns)) > 1.0e-6) &
!      .or.(abs(sum6 - a2l_src%forc_th(ns))   > 1.0e-6) &
       ) then
      write(6,*) 'clm_map2l check ERROR topo ',sum1,adomain%topo(ns)
      write(6,*) 'clm_map2l check ERROR t    ',sum2,a2l_src%forc_t(ns)
      write(6,*) 'clm_map2l check ERROR q    ',sum3,a2l_src%forc_q(ns)
      write(6,*) 'clm_map2l check ERROR hgt  ',sum4,a2l_src%forc_hgt(ns)
      write(6,*) 'clm_map2l check ERROR pbot ',sum5,a2l_src%forc_pbot(ns)
      write(6,*) 'clm_map2l check ERROR th   ',sum6,a2l_src%forc_th(ns)
!      call endrun()
    endif
  enddo

  endif   ! mx_ovr > 1

  first_call = .false.

end subroutine clm_mapa2l

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_mapl2a
!
! !INTERFACE:
  subroutine clm_mapl2a(l2a_src, l2a_dst)
!
! !DESCRIPTION:
! Maps lnd2atm fields from clm grid to external grid
!
! !USES:
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
  use areaMod  , only : map_maparrayl, map1dl_l2a
!
! !ARGUMENTS:
  implicit none
  type(lnd2atm_type), intent(in)  :: l2a_src
  type(lnd2atm_type), intent(out) :: l2a_dst
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: ix                    ! field index
  integer :: nflds                 ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
#if (defined DUST )
  integer :: m                     ! loop counter
#endif
!------------------------------------------------------------------------------

  nradflds = size(l2a_src%albd,dim=2)
  if (nradflds /= numrad) then
    write(6,*) 'clm_mapl2a ERROR: nradflds ne numrad ',nradflds,numrad
    call endrun()
  endif

  !--- allocate temporaries
  call get_proc_bounds    (begg_s, endg_s)
  call get_proc_bounds_atm(begg_d, endg_d)

  nflds = 12+2*numrad

#if (defined DUST || defined  PROGSSLT )  
  ! add on fv, ram1 
  nflds = nflds + 2 
#endif

#if (defined DUST  ) 
  ! add on the number of dust bins (for flxdust)
  nflds = nflds + ndst
#endif


  allocate(asrc(begg_s:endg_s,nflds))
  allocate(adst(begg_d:endg_d,nflds))

  ix = 0
  ix=ix+1; asrc(:,ix) = l2a_src%t_rad(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%t_ref2m(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%q_ref2m(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%h2osno(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%taux(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%tauy(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%eflx_lh_tot(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%eflx_sh_tot(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%eflx_lwrad_out(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%qflx_evap_tot(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%fsa(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%nee(:)  
  do n = 1,numrad
     ix=ix+1; asrc(:,ix) = l2a_src%albd(:,n)  
     ix=ix+1; asrc(:,ix) = l2a_src%albi(:,n)  
  enddo
#if (defined DUST || defined  PROGSSLT )
  ix=ix+1; asrc(:,ix) = l2a_src%ram1(:)  
  ix=ix+1; asrc(:,ix) = l2a_src%fv(:)
#endif
#if (defined DUST )
  do m = 1,ndst  ! dust bins
     ix=ix+1; asrc(:,ix) = l2a_src%flxdst(:,m)  
  end do !m
#endif

  call map_maparrayl(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, map1dl_l2a)

  ix = 0
  ix=ix+1; l2a_dst%t_rad(:)          = adst(:,ix)
  ix=ix+1; l2a_dst%t_ref2m(:)        = adst(:,ix)
  ix=ix+1; l2a_dst%q_ref2m(:)        = adst(:,ix)
  ix=ix+1; l2a_dst%h2osno(:)         = adst(:,ix)
  ix=ix+1; l2a_dst%taux(:)           = adst(:,ix)
  ix=ix+1; l2a_dst%tauy(:)           = adst(:,ix)
  ix=ix+1; l2a_dst%eflx_lh_tot(:)    = adst(:,ix)
  ix=ix+1; l2a_dst%eflx_sh_tot(:)    = adst(:,ix)
  ix=ix+1; l2a_dst%eflx_lwrad_out(:) = adst(:,ix)
  ix=ix+1; l2a_dst%qflx_evap_tot(:)  = adst(:,ix)
  ix=ix+1; l2a_dst%fsa(:)            = adst(:,ix)
  ix=ix+1; l2a_dst%nee(:)            = adst(:,ix)
  do n = 1,numrad
     ix=ix+1; l2a_dst%albd(:,n)      = adst(:,ix)
     ix=ix+1; l2a_dst%albi(:,n)      = adst(:,ix)
  enddo

#if (defined DUST || defined  PROGSSLT )
  ix=ix+1; l2a_dst%ram1(:)           = adst(:,ix)
  ix=ix+1; l2a_dst%fv(:)             = adst(:,ix)
#endif
#if (defined DUST  )
  do m = 1,ndst  ! dust bins
     ix=ix+1; l2a_dst%flxdst(:,m)    = adst(:,ix)
  end do !m
#endif

  deallocate(asrc)
  deallocate(adst)

end subroutine clm_mapl2a

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
  use decompMod   , only : get_proc_bounds
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
!EOP
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
!EOP
!
! !LOCAL VARIABLES:
  integer :: g                          ! indices
  type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
  type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
  type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
  type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
#if (defined DUST)
  type(pft_dflux_type),pointer :: pdf   ! local pointer to derived subtype
  integer n
#endif

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
          c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep
     do g = begg,endg
        clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
     end do
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albd, clm_l2a%albd,&
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albi, clm_l2a%albi,&
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
      
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out,&
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep
     do g = begg,endg
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
     end do

  else

     call c2g(begc, endc, begl, endl, begg, endg, cptr%cws%h2osno, clm_l2a%h2osno,&
          c2l_scale_type= 'unity', l2g_scale_type='unity')
!dir$ concurrent
!cdir nodep
     do g = begg,endg
        clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
     end do

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albd, clm_l2a%albd, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, numrad, &
          pptr%pps%albi, clm_l2a%albi, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pes%t_ref2m, clm_l2a%t_ref2m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pes%q_ref2m, clm_l2a%q_ref2m, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pmf%taux, clm_l2a%taux, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pmf%tauy, clm_l2a%tauy, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lh_tot, clm_l2a%eflx_lh_tot, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_sh_tot, clm_l2a%eflx_sh_tot, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pwf%qflx_evap_tot, clm_l2a%qflx_evap_tot, & 
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%fsa, clm_l2a%fsa, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
                  
     call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
          pptr%pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out, &
          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
                  
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

#if (defined DUST || defined  PROGSSLT )
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pps%fv, clm_l2a%fv, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, &
           pptr%pps%ram1, clm_l2a%ram1, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
#endif

#if (defined DUST )
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, ndst, &
           pptr%pdf%flx_mss_vrt_dst, clm_l2a%flxdst, &
           p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
#endif

     ! Convert from gC/m2/s to kgC/m2/s
!dir$ concurrent
!cdir nodep
     do g = begg,endg
        clm_l2a%nee(g) = clm_l2a%nee(g)*1.0e-3_r8
     end do

!dir$ concurrent
!cdir nodep
     do g = begg,endg
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
     end do

  end if

end subroutine clm_map2gcell

!------------------------------------------------------------------------
!------------------------------------------------------------------------
end module clm_atmlnd

