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
  use clm_varpar  , only : numrad
  use areaMod     , only : gridmap_type
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan
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
end type lnd2atm_type

  type(atm2lnd_type),public,target :: atm_a2l      ! a2l fields on atm grid
  type(lnd2atm_type),public,target :: atm_l2a      ! l2a fields on atm grid

  type(atm2lnd_type),public,target :: clm_a2l      ! a2l fields on clm grid
  type(lnd2atm_type),public,target :: clm_l2a      ! l2a fields on clm grid

  type(gridmap_type),public        :: gridmap_a2l  ! mapping from a2l
  type(gridmap_type),public        :: gridmap_l2a  ! mapping from l2a

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

! ival = nan      ! causes core dump in gridmap_maparray, tcx fix
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

! ival = nan      ! causes core dump in gridmap_maparray, tcx fix
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

end subroutine init_lnd2atm_type

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_mapa2l
!
! !INTERFACE:
  subroutine clm_mapa2l(a2l_src,a2l_dst,gridmap)
!
! !DESCRIPTION:
! Maps atm2lnd fields from external grid to clm grid
!
! !USES:
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
  use areaMod  , only : gridmap_maparray
!
! !ARGUMENTS:
  implicit none
  type(atm2lnd_type), intent(in)  :: a2l_src
  type(atm2lnd_type), intent(out) :: a2l_dst
  type(gridmap_type), intent(in)  :: gridmap
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: nflds                 ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  logical :: a2ltrue               ! a2l or l2a map type flag
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
!------------------------------------------------------------------------------

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

  asrc(:,1)  = a2l_src%forc_t(:)  
  asrc(:,2)  = a2l_src%forc_u(:)  
  asrc(:,3)  = a2l_src%forc_v(:)  
  asrc(:,4)  = a2l_src%forc_wind(:)  
  asrc(:,5)  = a2l_src%forc_q(:)  
  asrc(:,6)  = a2l_src%forc_hgt(:)  
  asrc(:,7)  = a2l_src%forc_hgt_u(:)  
  asrc(:,8)  = a2l_src%forc_hgt_t(:)  
  asrc(:,9)  = a2l_src%forc_hgt_q(:)  
  asrc(:,10) = a2l_src%forc_pbot(:)  
  asrc(:,11) = a2l_src%forc_th(:)  
  asrc(:,12) = a2l_src%forc_vp(:)  
  asrc(:,13) = a2l_src%forc_rho(:)  
  asrc(:,14) = a2l_src%forc_psrf(:)  
  asrc(:,15) = a2l_src%forc_pco2(:)  
  asrc(:,16) = a2l_src%forc_lwrad(:)  
  asrc(:,17) = a2l_src%forc_solar(:)  
  asrc(:,18) = a2l_src%forc_rain(:)  
  asrc(:,19) = a2l_src%forc_snow(:)  
  asrc(:,20) = a2l_src%forc_pc13o2(:)  
  asrc(:,21) = a2l_src%forc_po2(:)  
  do n = 1,numrad
     asrc(:,20+2*n) = a2l_src%forc_solad(:,n)  
     asrc(:,21+2*n) = a2l_src%forc_solai(:,n)  
  enddo
!-forc_ndep is not recd from atm,don't know why it's in a2l (TCFIX) ---
!-forc_ndep cannot be updated here, array will be trashed and CN will fail ---
! call gridmap_maparray(begg_s, endg_s, begg_d, endg_d, a2l_src%forc_ndep  ,a2l_dst%forc_ndep  ,gridmap,a2ltrue)
!  asrc(:,xx) = a2l_src%forc_ndep(:)  

  a2ltrue = .true.
#if (defined OFFLINE)
  call gridmap_maparray(begg_s, endg_s, begg_d, endg_d, 1, a2l_src%flfall, a2l_dst%flfall, gridmap, a2ltrue)
#endif
  call gridmap_maparray(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, gridmap, a2ltrue)

  a2l_dst%forc_t(:)     =   adst(:,1)
  a2l_dst%forc_u(:)     =   adst(:,2)
  a2l_dst%forc_v(:)     =   adst(:,3)
  a2l_dst%forc_wind(:)  =   adst(:,4)
  a2l_dst%forc_q(:)     =   adst(:,5)
  a2l_dst%forc_hgt(:)   =   adst(:,6)
  a2l_dst%forc_hgt_u(:) =   adst(:,7)
  a2l_dst%forc_hgt_t(:) =   adst(:,8)
  a2l_dst%forc_hgt_q(:) =   adst(:,9)
  a2l_dst%forc_pbot(:)  =   adst(:,10)
  a2l_dst%forc_th(:)    =   adst(:,11)
  a2l_dst%forc_vp(:)    =   adst(:,12)
  a2l_dst%forc_rho(:)   =   adst(:,13)
  a2l_dst%forc_psrf(:)  =   adst(:,14)
  a2l_dst%forc_pco2(:)  =   adst(:,15)
  a2l_dst%forc_lwrad(:) =   adst(:,16)
  a2l_dst%forc_solar(:) =   adst(:,17)
  a2l_dst%forc_rain(:)  =   adst(:,18)
  a2l_dst%forc_snow(:)  =   adst(:,19)
  a2l_dst%forc_pc13o2(:)=   adst(:,20)
  a2l_dst%forc_po2(:)   =   adst(:,21)
  do n = 1,numrad
     a2l_dst%forc_solad(:,n)  = adst(:,20+2*n)
     a2l_dst%forc_solai(:,n)  = adst(:,21+2*n)
  enddo

  deallocate(asrc)
  deallocate(adst)

end subroutine clm_mapa2l

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_mapl2a
!
! !INTERFACE:
  subroutine clm_mapl2a(l2a_src,l2a_dst,gridmap)
!
! !DESCRIPTION:
! Maps lnd2atm fields from clm grid to external grid
!
! !USES:
  use decompMod, only : get_proc_bounds, get_proc_bounds_atm
  use areaMod  , only : gridmap_maparray
!
! !ARGUMENTS:
  implicit none
  type(lnd2atm_type), intent(in)  :: l2a_src
  type(lnd2atm_type), intent(out) :: l2a_dst
  type(gridmap_type), intent(in)  :: gridmap
!
! !REVISION HISTORY:
! 2005.11.15  T Craig  Creation.
! 2006.3.30   P Worley Restructuring for improved vector performance
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: n                     ! loop counter
  integer :: nflds                 ! number of fields to be mapped
  integer :: nradflds              ! size of 2nd dim in arrays
  integer :: begg_s,endg_s         ! beg,end of input grid
  integer :: begg_d,endg_d         ! beg,end of output grid
  logical :: a2lfalse              ! a2l or l2a map type flag
  real(r8),pointer :: asrc(:,:)    ! temporary source data
  real(r8),pointer :: adst(:,:)    ! temporary dest data
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

  allocate(asrc(begg_s:endg_s,nflds))
  allocate(adst(begg_d:endg_d,nflds))

  asrc(:,1)  = l2a_src%t_rad(:)  
  asrc(:,2)  = l2a_src%t_ref2m(:)  
  asrc(:,3)  = l2a_src%q_ref2m(:)  
  asrc(:,4)  = l2a_src%h2osno(:)  
  asrc(:,5)  = l2a_src%taux(:)  
  asrc(:,6)  = l2a_src%tauy(:)  
  asrc(:,7)  = l2a_src%eflx_lh_tot(:)  
  asrc(:,8)  = l2a_src%eflx_sh_tot(:)  
  asrc(:,9)  = l2a_src%eflx_lwrad_out(:)  
  asrc(:,10)  = l2a_src%qflx_evap_tot(:)  
  asrc(:,11)  = l2a_src%fsa(:)  
  asrc(:,12)  = l2a_src%nee(:)  
  do n = 1,numrad
     asrc(:,11+2*n)  = l2a_src%albd(:,n)  
     asrc(:,12+2*n)  = l2a_src%albi(:,n)  
  enddo

  a2lfalse = .false.
  call gridmap_maparray(begg_s, endg_s, begg_d, endg_d, nflds, asrc, adst, gridmap, a2lfalse)

  l2a_dst%t_rad(:)          = adst(:,1)
  l2a_dst%t_ref2m(:)        = adst(:,2)
  l2a_dst%q_ref2m(:)        = adst(:,3)
  l2a_dst%h2osno(:)         = adst(:,4)
  l2a_dst%taux(:)           = adst(:,5)
  l2a_dst%tauy(:)           = adst(:,6)
  l2a_dst%eflx_lh_tot(:)    = adst(:,7)
  l2a_dst%eflx_sh_tot(:)    = adst(:,8)
  l2a_dst%eflx_lwrad_out(:) = adst(:,9)
  l2a_dst%qflx_evap_tot(:)  = adst(:,10)
  l2a_dst%fsa(:)            = adst(:,11)
  l2a_dst%nee(:)            = adst(:,12)
  do n = 1,numrad
     l2a_dst%albd(:,n)      = adst(:,11+2*n)
     l2a_dst%albi(:,n)      = adst(:,12+2*n)
  enddo

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

