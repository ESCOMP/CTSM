module clm_atmlnd

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd, lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar    , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon    , only : rair, grav, cpair, hfus, tfrz
  use clm_varctl    , only : iulog, use_c13, use_cn, use_lch4
  use seq_drydep_mod, only : n_drydep, drydep_method, DD_XLND
  use shr_megan_mod , only : shr_megan_mechcomps_n
  use decompMod     , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  !----------------------------------------------------
  ! atmosphere -> land variables structure
  !----------------------------------------------------
  type, public :: atm2lnd_type
     real(r8), pointer :: forc_t(:)        => null() !atmospheric temperature (Kelvin)
     real(r8), pointer :: forc_u(:)        => null() !atm wind speed, east direction (m/s)
     real(r8), pointer :: forc_v(:)        => null() !atm wind speed, north direction (m/s)
     real(r8), pointer :: forc_wind(:)     => null() !atmospheric wind speed   
     real(r8), pointer :: forc_q(:)        => null() !atmospheric specific humidity (kg/kg)
     real(r8), pointer :: forc_hgt(:)      => null() !atmospheric reference height (m)
     real(r8), pointer :: forc_hgt_u(:)    => null() !obs height of wind [m] (new)
     real(r8), pointer :: forc_hgt_t(:)    => null() !obs height of temperature [m] (new)
     real(r8), pointer :: forc_hgt_q(:)    => null() !obs height of humidity [m] (new)
     real(r8), pointer :: forc_pbot(:)     => null() !atmospheric pressure (Pa)
     real(r8), pointer :: forc_th(:)       => null() !atm potential temperature (Kelvin)
     real(r8), pointer :: forc_vp(:)       => null() !atmospheric vapor pressure (Pa) 
     real(r8), pointer :: forc_rho(:)      => null() !density (kg/m**3)
     real(r8), pointer :: forc_rh(:)       => null() !atmospheric relative humidity (%)
     real(r8), pointer :: forc_psrf(:)     => null() !surface pressure (Pa)
     real(r8), pointer :: forc_pco2(:)     => null() !CO2 partial pressure (Pa)
     real(r8), pointer :: forc_lwrad(:)    => null() !downwrd IR longwave radiation (W/m**2)
     real(r8), pointer :: forc_solad(:,:)  => null() !direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai(:,:)  => null() !diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
     real(r8), pointer :: forc_solar(:)    => null() !incident solar radiation
     real(r8), pointer :: forc_rain(:)     => null() !rain rate [mm/s]
     real(r8), pointer :: forc_snow(:)     => null() !snow rate [mm/s]
     real(r8), pointer :: forc_ndep(:)     => null() !nitrogen deposition rate (gN/m2/s)
     real(r8), pointer :: rainf(:)         => null() !ALMA rain+snow [mm/s]
     real(r8), pointer :: forc_pc13o2(:)   => null() !C13O2 partial pressure (Pa)
     real(r8), pointer :: forc_po2(:)      => null() !O2 partial pressure (Pa)
     real(r8), pointer :: forc_flood(:)    => null() !rof flood (mm/s)
     real(r8), pointer :: volr(:)          => null() !rof volr (m3)
     real(r8), pointer :: forc_aer(:,:)    => null() !aerosol deposition array
     real(r8), pointer :: forc_pch4(:)     => null() !CH4 partial pressure (Pa)
  end type atm2lnd_type

  !----------------------------------------------------
  ! land -> atmosphere variables structure
  !----------------------------------------------------
  type, public :: lnd2atm_type
     real(r8), pointer :: t_rad(:)         => null() !radiative temperature (Kelvin)
     real(r8), pointer :: t_ref2m(:)       => null() !2m surface air temperature (Kelvin)
     real(r8), pointer :: q_ref2m(:)       => null() !2m surface specific humidity (kg/kg)
     real(r8), pointer :: u_ref10m(:)      => null() !10m surface wind speed (m/sec)
     real(r8), pointer :: h2osno(:)        => null() !snow water (mm H2O)
     real(r8), pointer :: albd(:,:)        => null() !(numrad) surface albedo (direct)
     real(r8), pointer :: albi(:,:)        => null() !(numrad) surface albedo (diffuse)
     real(r8), pointer :: taux(:)          => null() !wind stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy(:)          => null() !wind stress: n-s (kg/m/s**2)
     real(r8), pointer :: eflx_lh_tot(:)   => null() !total latent HF (W/m**2)  [+ to atm]
     real(r8), pointer :: eflx_sh_tot(:)   => null() !total sensible HF (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lwrad_out(:)=> null() !IR (longwave) radiation (W/m**2)
     real(r8), pointer :: qflx_evap_tot(:) => null() !qflx_evap_soi + qflx_evap_can + qflx_tran_veg
     real(r8), pointer :: fsa(:)           => null() !solar rad absorbed (total) (W/m**2)
     real(r8), pointer :: nee(:)           => null() !net CO2 flux (kg CO2/m**2/s) [+ to atm]
     real(r8), pointer :: ram1(:)          => null() !aerodynamical resistance (s/m)
     real(r8), pointer :: fv(:)            => null() !friction velocity (m/s) (for dust model)
     real(r8), pointer :: h2osoi_vol(:,:)  => null() !volumetric soil water (0~watsat, m3/m3, nlevgrnd) (for dust model)
     real(r8), pointer :: rofliq(:)        => null() !rof liq forcing
     real(r8), pointer :: rofice(:)        => null() !rof ice forcing
     real(r8), pointer :: flxdst(:,:)      => null() !dust flux (size bins)
     real(r8), pointer :: ddvel(:,:)       => null() !dry deposition velocities
     real(r8), pointer :: flxvoc(:,:)      => null() !VOC flux (size bins)
     real(r8), pointer :: flux_ch4(:)      => null() !net CH4 flux (kg C/m**2/s) [+ to atm]
  end type lnd2atm_type
  
  type(atm2lnd_type),public,target :: clm_a2l      ! a2l fields on clm grid
  type(lnd2atm_type),public,target :: clm_l2a      ! l2a fields on clm grid

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_atm2lnd_type
  public :: init_lnd2atm_type
  public :: clm_map2gcell
  !----------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_atm2lnd_type(bounds, a2l)
    !
    ! !DESCRIPTION:
    ! Initialize atmospheric variables required by the land
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    type (atm2lnd_type), intent(inout):: a2l
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival   ! initial value
    !------------------------------------------------------------------------

    ! ival = nan      ! causes core dump in map_maparray, tcx fix
    ival = 0.0_r8

    allocate(a2l%forc_t(bounds%begg:bounds%endg))
    a2l%forc_t(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_u(bounds%begg:bounds%endg))
    a2l%forc_u(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_v(bounds%begg:bounds%endg))
    a2l%forc_v(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_wind(bounds%begg:bounds%endg))
    a2l%forc_wind(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_q(bounds%begg:bounds%endg))
    a2l%forc_q(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_rh(bounds%begg:bounds%endg))
    a2l%forc_rh(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_hgt(bounds%begg:bounds%endg))
    a2l%forc_hgt(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_hgt_u(bounds%begg:bounds%endg))
    a2l%forc_hgt_u(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_hgt_t(bounds%begg:bounds%endg))
    a2l%forc_hgt_t(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_hgt_q(bounds%begg:bounds%endg))
    a2l%forc_hgt_q(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_pbot(bounds%begg:bounds%endg))
    a2l%forc_pbot(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_th(bounds%begg:bounds%endg))
    a2l%forc_th(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_vp(bounds%begg:bounds%endg))
    a2l%forc_vp(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_rho(bounds%begg:bounds%endg))
    a2l%forc_rho(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_psrf(bounds%begg:bounds%endg))
    a2l%forc_psrf(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_pco2(bounds%begg:bounds%endg))
    a2l%forc_pco2(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_lwrad(bounds%begg:bounds%endg))
    a2l%forc_lwrad(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_solad(bounds%begg:bounds%endg,numrad))
    a2l%forc_solad(bounds%begg:bounds%endg,numrad)=ival
    allocate(a2l%forc_solai(bounds%begg:bounds%endg,numrad))
    a2l%forc_solai(bounds%begg:bounds%endg,numrad)=ival
    allocate(a2l%forc_solar(bounds%begg:bounds%endg))
    a2l%forc_solar(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_rain(bounds%begg:bounds%endg))
    a2l%forc_rain(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_snow(bounds%begg:bounds%endg))
    a2l%forc_snow(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_ndep(bounds%begg:bounds%endg))
    a2l%forc_ndep(bounds%begg:bounds%endg)=ival
    allocate(a2l%rainf(bounds%begg:bounds%endg))
    a2l%rainf(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_pc13o2(bounds%begg:bounds%endg))
    a2l%forc_pc13o2(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_po2(bounds%begg:bounds%endg))
    a2l%forc_po2(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_flood(bounds%begg:bounds%endg))
    a2l%forc_flood(bounds%begg:bounds%endg)=ival
    allocate(a2l%volr(bounds%begg:bounds%endg))
    a2l%volr(bounds%begg:bounds%endg)=ival
    allocate(a2l%forc_aer(bounds%begg:bounds%endg,14))
    a2l%forc_aer(bounds%begg:bounds%endg,14)=ival
    allocate(a2l%forc_pch4(bounds%begg:bounds%endg))
    a2l%forc_pch4(bounds%begg:bounds%endg)=ival

  end subroutine init_atm2lnd_type

  !------------------------------------------------------------------------
  subroutine init_lnd2atm_type(bounds, l2a)
    !
    ! !DESCRIPTION:
    ! Initialize land variables required by the atmosphere
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    type (lnd2atm_type), intent(inout):: l2a
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival   ! initial value
    !------------------------------------------------------------------------

    ! ival = nan   ! causes core dump in map_maparray, tcx fix
    ival = 0.0_r8

    allocate(l2a%t_rad(bounds%begg:bounds%endg))
    l2a%t_rad(bounds%begg:bounds%endg)=ival
    allocate(l2a%t_ref2m(bounds%begg:bounds%endg))
    l2a%t_ref2m(bounds%begg:bounds%endg)=ival
    allocate(l2a%q_ref2m(bounds%begg:bounds%endg))
    l2a%q_ref2m(bounds%begg:bounds%endg)=ival
    allocate(l2a%u_ref10m(bounds%begg:bounds%endg))
    l2a%u_ref10m(bounds%begg:bounds%endg)=ival
    allocate(l2a%h2osno(bounds%begg:bounds%endg))
    l2a%h2osno(bounds%begg:bounds%endg)=ival
    allocate(l2a%albd(bounds%begg:bounds%endg,1:numrad))
    l2a%albd(bounds%begg:bounds%endg,1:numrad)=ival
    allocate(l2a%albi(bounds%begg:bounds%endg,1:numrad))
    l2a%albi(bounds%begg:bounds%endg,1:numrad)=ival
    allocate(l2a%taux(bounds%begg:bounds%endg))
    l2a%taux(bounds%begg:bounds%endg)=ival
    allocate(l2a%tauy(bounds%begg:bounds%endg))
    l2a%tauy(bounds%begg:bounds%endg)=ival
    allocate(l2a%eflx_lwrad_out(bounds%begg:bounds%endg))
    l2a%eflx_lwrad_out(bounds%begg:bounds%endg)=ival
    allocate(l2a%eflx_sh_tot(bounds%begg:bounds%endg))
    l2a%eflx_sh_tot(bounds%begg:bounds%endg)=ival
    allocate(l2a%eflx_lh_tot(bounds%begg:bounds%endg))
    l2a%eflx_lh_tot(bounds%begg:bounds%endg)=ival
    allocate(l2a%qflx_evap_tot(bounds%begg:bounds%endg))
    l2a%qflx_evap_tot(bounds%begg:bounds%endg)=ival
    allocate(l2a%fsa(bounds%begg:bounds%endg))
    l2a%fsa(bounds%begg:bounds%endg)=ival
    allocate(l2a%nee(bounds%begg:bounds%endg))
    l2a%nee(bounds%begg:bounds%endg)=ival
    allocate(l2a%ram1(bounds%begg:bounds%endg))
    l2a%ram1(bounds%begg:bounds%endg)=ival
    allocate(l2a%fv(bounds%begg:bounds%endg))
    l2a%fv(bounds%begg:bounds%endg)=ival
    allocate(l2a%h2osoi_vol(bounds%begg:bounds%endg,1:nlevgrnd))
    l2a%h2osoi_vol(bounds%begg:bounds%endg,1:nlevgrnd)=ival
    allocate(l2a%rofliq(bounds%begg:bounds%endg))
    l2a%rofliq(bounds%begg:bounds%endg)=ival
    allocate(l2a%rofice(bounds%begg:bounds%endg))
    l2a%rofice(bounds%begg:bounds%endg)=ival
    allocate(l2a%flxdst(bounds%begg:bounds%endg,1:ndst))
    l2a%flxdst(bounds%begg:bounds%endg,1:ndst)=ival
    allocate(l2a%flux_ch4(bounds%begg:bounds%endg))
    l2a%flux_ch4(bounds%begg:bounds%endg)=ival
    if (shr_megan_mechcomps_n>0) then
       allocate(l2a%flxvoc(bounds%begg:bounds%endg,1:shr_megan_mechcomps_n))
       l2a%flxvoc(bounds%begg:bounds%endg,1:shr_megan_mechcomps_n)=ival
    endif
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(l2a%ddvel(bounds%begg:bounds%endg,1:n_drydep))
       l2a%ddvel(bounds%begg:bounds%endg,1:n_drydep)=ival
    end if

  end subroutine init_lnd2atm_type

  !------------------------------------------------------------------------
  
  subroutine clm_map2gcell(bounds, init)
    !
    ! !DESCRIPTION:
    ! Compute l2a component of gridcell derived type
    !
    ! !USES:
    use clmtype
    use subgridAveMod
    use clm_varcon  , only : sb
    use clm_varpar  , only : numrad
    use ch4varcon   , only : ch4offline
    !
    ! !ARGUMENTS:
    implicit none
    save
    type(bounds_type), intent(in) :: bounds  ! bounds
    logical, optional, intent(in) :: init  ! if true=>only set a subset of arguments
    !
    ! !LOCAL VARIABLES:
    integer :: g,n             ! indices
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    if (present(init)) then

       call c2g(bounds, cws%h2osno, clm_l2a%h2osno, &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       do g = bounds%begg,bounds%endg
          clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
       end do

       call c2g(bounds, nlevgrnd, cws%h2osoi_vol, clm_l2a%h2osoi_vol, &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       call p2g(bounds, numrad, pps%albd, clm_l2a%albd,&
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       call p2g(bounds, numrad, pps%albi, clm_l2a%albi,&
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       call p2g(bounds, pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out,&
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       do g = bounds%begg,bounds%endg
          clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
       end do

    else

       call c2g(bounds, cws%h2osno, clm_l2a%h2osno,&
            c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       do g = bounds%begg,bounds%endg
          clm_l2a%h2osno(g) = clm_l2a%h2osno(g)/1000._r8
       end do

       call p2g(bounds, numrad, pps%albd, clm_l2a%albd, &
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       call p2g(bounds, numrad, pps%albi, clm_l2a%albi, &
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       call p2g(bounds, pes%t_ref2m, clm_l2a%t_ref2m, & 
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

       call p2g(bounds, pes%q_ref2m, clm_l2a%q_ref2m, & 
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

       call p2g(bounds, pps%u10_clm, clm_l2a%u_ref10m, & 
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

       call p2g(bounds, pmf%taux, clm_l2a%taux, & 
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

       call p2g(bounds, pmf%tauy, clm_l2a%tauy, & 
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

       call p2g(bounds, pef%eflx_lh_tot, clm_l2a%eflx_lh_tot, &
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       do g = bounds%begg,bounds%endg
          clm_l2a%eflx_sh_tot(g) = gef%eflx_sh_totg(g)
       end do

       call p2g(bounds, pwf%qflx_evap_tot, clm_l2a%qflx_evap_tot, & 
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       call p2g(bounds, pef%fsa, clm_l2a%fsa, &
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       call p2g(bounds, pef%eflx_lwrad_out, clm_l2a%eflx_lwrad_out, &
            p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       if (use_cn) then
          call c2g(bounds, ccf%nee, clm_l2a%nee, &
               c2l_scale_type= 'unity', l2g_scale_type='unity')
       else
          call p2g(bounds, pcf%fco2, clm_l2a%nee, &
               p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
          ! Note that fco2 in is umolC/m2/sec so units need to be changed to gC/m2/sec
          do g = bounds%begg,bounds%endg
             clm_l2a%nee(g) = clm_l2a%nee(g)*12.011e-6_r8
          end do
       end if

       if (use_lch4) then
          if (.not. ch4offline) then
             ! Adjust flux of CO2 by the net conversion of mineralizing C to CH4
             do g = bounds%begg,bounds%endg
                clm_l2a%nee(g) = clm_l2a%nee(g) + gch4%nem(g) ! nem is in g C/m2/sec
                ! nem is calculated in ch4Mod
                ! flux_ch4 is averaged there also.
             end do
          end if
       end if

       call p2g(bounds, pps%fv, clm_l2a%fv, &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

       call p2g(bounds, pps%ram1, clm_l2a%ram1, &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

       do g = bounds%begg,bounds%endg
          clm_l2a%rofliq(g) = gwf%qflx_runoffg(g)
          clm_l2a%rofice(g) = gwf%qflx_snwcp_iceg(g)
       end do

       call p2g(bounds, ndst, &
            pdf%flx_mss_vrt_dst, clm_l2a%flxdst, &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

       if (shr_megan_mechcomps_n>0) then
          call p2g(bounds, shr_megan_mechcomps_n, pvf%vocflx, clm_l2a%flxvoc, &
               p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
       endif

       if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
          call p2g(bounds, n_drydep, pdd%drydepvel, clm_l2a%ddvel, &
               p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
       endif

       ! Convert from gC/m2/s to kgCO2/m2/s
       do g = bounds%begg,bounds%endg
          clm_l2a%nee(g) = clm_l2a%nee(g)*convertgC2kgCO2
       end do

       do g = bounds%begg,bounds%endg
          clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
       end do

    end if

  end subroutine clm_map2gcell

end module clm_atmlnd

