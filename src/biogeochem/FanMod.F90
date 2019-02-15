module FanMod
#ifdef _PYMOD_
  use qsatmod
#else
  use shr_const_mod
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use QSatMod                         , only : QSat
#endif
  implicit none

#ifdef _REALPR_
  integer, parameter :: r8 = 8
#endif
  
#ifdef _PYMOD_
  public

  real(r8), parameter :: SHR_CONST_BOLTZ   = 1.38065e-23
  real(r8), parameter :: SHR_CONST_AVOGAD  = 6.02214e26
  real(r8), parameter :: SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ
  real(r8), parameter :: SHR_CONST_MWDAIR  = 28.966
  real(r8), parameter :: SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR

#else
  private
  public update_org_n
  public eval_fluxes_storage
  public update_npool
  public update_4pool
  public update_urea
#endif
  
  ! Indices in flux arrays, soil:
  integer, parameter, public :: iflx_air = 1, & ! flux to air
       iflx_soild = 2, & ! diffusion to soil
       iflx_no3 = 3, &   ! nitrification
       iflx_soilq = 4, & ! percolation to soil
       iflx_roff = 5, &  ! surface runoff
       iflx_to_tan = 6   ! conversion to tan (from urea) 
  ! Indices in flux arrays, storage:
  integer, parameter, public :: iflx_air_barns = 1, &
       iflx_air_stores = 2, &
       iflx_appl = 3, &
       iflx_to_store = 4
  ! Indices in the organic N pools and fluxes
  integer, parameter, public :: ind_avail = 1, ind_resist = 2, ind_unavail = 3 
  
  ! nominal depth where the soil TAN concentration vanishes:
  real(r8), parameter, public :: soildepth_reservoir = 0.04_r8

  integer, parameter, public :: err_bad_theta = 1, err_negative_tan = 2, err_negative_flux = 3, &
       err_balance_tan = 4, err_balance_nitr = 5, err_nan = 6, err_bad_subst = 7, err_bad_type = 8

  integer, parameter, public :: subst_tan = 1, subst_urea = 2
  
  real(r8), parameter, public :: water_relax_t = 24*3600.0_r8
  
contains

  ! accessor functions are only needed for the python interface... 
  !
  integer function ind_soild() result(ind)
    ind = iflx_soild
  end function ind_soild

  integer function ind_soilq() result(ind)
    ind = iflx_soilq
  end function ind_soilq

  integer function ind_air() result(ind)
    ind = iflx_air
  end function ind_air

  integer function ind_no3() result(ind)
    ind = iflx_no3
  end function ind_no3

  integer function ind_roff() result(ind)
    ind = iflx_roff
  end function ind_roff

  integer function ind_air_barns() result(ind)
    ind = iflx_air_barns
  end function ind_air_barns

  integer function ind_air_stores() result(ind)
    ind = iflx_air_stores
  end function ind_air_stores

  integer function ind_appl() result(ind)
    ind = iflx_appl
  end function ind_appl

  integer function ind_to_store() result(ind)
    ind = iflx_to_store
  end function ind_to_store
  
  function eval_diffusivity_liq_mq(theta, thetasat, tg) result(diff)
    ! Evaluate the aquous phase diffusivity for TAN in soil according to the Millington &
    ! Quirk model.
    implicit none
    real(r8), intent(in) :: theta ! volumetric water content, m3/m3   
    real(r8), intent(in) :: thetasat ! theta at saturation
    real(r8), intent(in) :: tg ! soil temperature, K

    real(r8) :: diff
    
    real(r8) :: kaq_base
    real(r8), parameter :: pw = 7.0_r8 / 3.0_r8
    real(r8), parameter :: gascnst = 8.314, faraday = 96500.0_r8, lp = 73.4_r8, lm_no3 = 71.4_r8, lm_oh=197.6_r8, lm=lm_no3

    ! Base rate by Nernst-Haskell equation, see Poling et al., 2000. The Properties of
    ! Gases and Liquids.
    kaq_base = 1e-4 * (gascnst*tg / (2*faraday**2)) / (1/lp + 1/lm)

    ! Van Der Molen 1990 fit of the base rate.
    !kaq_base = 9.8e-10_r8 * 1.03_r8 ** (Tg-273.0_r8)

    diff = kaq_base * (theta**pw) / (thetasat**2)

  end function eval_diffusivity_liq_mq

  function eval_diffusivity_gas_mq(theta, thetasat, tg) result(diff)
    ! Evaluate the gas phase diffusivity for NH3 in soil according to the Millington &
    ! Quirk model.
    implicit none
    real(r8), intent(in) :: theta ! volumetric water content, m3/m3   
    real(r8), intent(in) :: thetasat ! theta at saturation
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8) :: diff ! diffusivity, m2/s
    
    real(r8) :: soilair, dair
    real(r8), parameter :: pw = 7.0_r8 / 3.0_r8
    real(r8), parameter :: mNH3 = 17., mair = 29, vNH3 = 14.9, vair = 20.1, press = 1.0
    
    soilair = thetasat - theta

    ! Van Der Molen 1990 fit of the base rate.
    !dair = 1.7e-5_r8 * 1.03_r8**(Tg-293.0_r8)

    !dair = 18e-6_r8
    !dair = 1.4e-5
    
    ! Base rate from the Fuller et al. 1966 method.
    dair = (0.001 * tg**1.75 * sqrt(1/mNH3 + 1/mair)) / (press * (vair**(1./3) * vNH3**(1./3))**2) * 1e-4
    diff = dair * (soilair**pw) / (thetasat**2)

  end function eval_diffusivity_gas_mq

  ! The moldrup 2003 are here but not used currently. Check the base rates if use these.
  
  function eval_diffusivity_gas_m03(theta, thetasat, tg, bsw) result(diff)
    ! Evaluate the gas phase diffusivity for NH3 in soil according to the method of
    ! Moldrup (2003).
    implicit none
    real(r8), intent(in) :: theta, thetasat, tg, bsw
    real(r8) :: diff
    
    real(r8) :: soilair, dair, m03_W
    real(r8), parameter :: pw = 10.0_r8 / 3.0_r8
    real(r8), parameter :: m03_T = 2.0
    real(r8), parameter :: mNH3 = 17., mair = 29, vNH3 = 14.9, vair = 20.1, press = 1.0

    m03_W = 3.0 / bsw
    
    soilair = thetasat - theta
    !dair = 1.7e-5 * 1.03**(Tg-293.0_r8)
    dair = (0.001 * tg**1.75 * sqrt(1/mNH3 + 1/mair)) / (press * (vair**(1./3) * vNH3**(1./3))**2) * 1e-4
   
    diff = dair * soilair**m03_T * (soilair/thetasat)**m03_W / soilair
    
  end function eval_diffusivity_gas_m03

  function eval_diffusivity_liq_m03(theta, thetasat, tg, bsw) result(diff)
    ! Evaluate the aquous phase diffusivity for TAN in soil according to the method of Moldrup (2003).
    implicit none
    real(r8), intent(in) :: theta, thetasat, tg, bsw
    real(r8) :: diff
    
    real(r8) :: kaq_base, m03_W
    real(r8), parameter :: pw = 10.0_r8 / 3.0_r8, m03_T = 2.0_r8
    real(r8), parameter :: gascnst = 8.314, faraday = 96500.0_r8, lp = 73.4_r8, lm_no3 = 71.4_r8, lm_oh=197.6_r8, lm=lm_oh
    
    kaq_base = 9.8e-10 * 1.03 ** (Tg-273.0_r8)
    
    !kaq_base = 1e-4 * (gascnst*tg / (2*faraday**2)) / (1/lp + 1/lm)
    
    m03_W = 0.3333_r8*bsw - 1.0
    
    diff = kaq_base * theta**m03_T * (theta/thetasat)**m03_W / theta

  end function eval_diffusivity_liq_m03

  subroutine partition_tan(tg, Hconc, theta, air, kads, KNH3, fract_nh4)
    ! Partition the bulk TAN (NH3 gas/aq + NH4 (aq)) between gas and aqueous. Outputs the
    ! volatility (gas/aq ratio), see below. 
    implicit none
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: Hconc ! H+ concentration, mol / l
    real(r8), intent(in) :: theta ! water volume fraction
    real(r8), intent(in) :: air   ! air volume fraction
    real(r8), intent(in) :: kads  ! adsorption coefficient, kads = NH4(ads) / NH4(aq)

    real(r8), optional, intent(out) :: fract_nh4 ! mass fraction of NH4.
    real(r8), intent(out) :: KNH3 ! volatility, ratio of concentrations [NH3 (gas)] / [NH3+NH4 (aq)] 
    
    real(r8) :: KNH4, HNH3, cnc_aq, fract_aq
    real(r8), parameter :: Tref = 298.15_r8

    KNH4 = 5.67_r8 * 1e-10_r8 * exp(-6286.0_r8 * (1.0_r8/Tg - 1.0_r8/Tref))
    ! HNH3 = [aq] / [gas] -- solubility
    HNH3 = 4.59_r8 * Tg * exp(4092_r8 * (1.0_r8/Tg - 1.0_r8/Tref))
    
    KNH3 = KNH4 / (HNH3*(KNH4 + Hconc))

    fract_aq = theta / (air*KNH3 + theta + kads*(1.0_r8-theta-air))
    if (present(fract_nh4)) fract_nh4 = fract_aq * Hconc / (KNH4 + Hconc)
    
    !if (present(fract_nh3g)) fract_nh3g = air  * KNH3 / (air*KNH3 + theta)
    !if (present(fract_nh3aq)) fract_nh3aq = fract_aq * (1.0_r8 - Hconc / (KNH4 + Hconc))

  end subroutine partition_tan

  real(r8) function eval_no3prod_v2(theta, theta_sat, Tg) result(kNO3)
    ! Evaluate nitrification rate as in the Riddick et al. (2016) paper but for NH4.
    ! Partitioning between TAN forms is not included.
    real(r8), intent(in) :: theta, theta_sat ! volumetric soil water m/m
    real(r8), intent(in) :: Tg    ! soil temperature, K

    real(r8) :: stf, wmr, smrf, mNH4, soil_dens

    !real(r8), parameter :: soil_dens = 1400.0_r8 ! Soil density, kg/m3
    real(r8), parameter :: water_dens = 1000.0_r8
    real(r8), parameter :: rmax = 1.16e-6_r8   ! Maximum rate of nitrification, s-1
    real(r8), parameter :: tmax = 313.0       ! Maximm temperature of microbial activity, K
    real(r8), parameter :: topt = 301.0       ! Optimal temperature of microbial acticity, K
    real(r8), parameter :: asg  = 2.4_8        ! a_sigma, empirical factor
    real(r8), parameter :: wmr_crit = 0.12_r8  ! Critical water content, g/g
    real(r8), parameter :: smrf_b = 2          ! Parameter in soil moisture response function
    real(r8), parameter :: Tref = 298.15_r8
    

    soil_dens = 2650.0_r8 * (1.0_r8-theta_sat)
    mNH4 = 1.0_r8

    ! soil temperature function
    stf = (max(1e-3_r8, tmax-Tg) / (tmax-topt))**asg * exp(asg * (Tg-topt)/(tmax-topt))

    ! gravimetric soil water
    wmr = theta * water_dens / soil_dens

    ! soil moisture response function
    
    smrf = 1.0_r8 - exp(-(wmr/wmr_crit)**smrf_b)
    !if stf < 1e-9 or smrf < 1e-9:
    !    print theta
    !    1/0
    kNO3 = 2.0_r8 * rmax * mNH4 / (1.0_r8/stf + 1.0_r8/smrf)
    
  end function eval_no3prod_v2
  
  subroutine eval_fluxes_slurry(water, mtan, Hconc, tg, ratm, theta, thetasat, perc, runoff, bsw, kads, fluxes)
    ! Evaluate nitrogen fluxes for a partly infiltrated layer of slurry.
    ! The state of infiltration is detemined from the amounts water on surface and in soil.
    ! Positive flux means loss of TAN.
    implicit none
    real(r8), intent(in) :: water(2) ! water in (surface , subsurface), m
    real(r8), intent(in) :: mtan   ! TAN, mass units / m2, surface + subsurface
    real(r8), intent(out) :: fluxes(5) ! TAN fluxes, see top of the module
    real(r8), intent(in) :: Hconc    ! H+ concentration, -log10(pH)
    real(r8), intent(in) :: tg       ! soil temperature, K
    real(r8), intent(in) :: ratm     ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta    ! volumetric soil water in "clean" soil
    real(r8), intent(in) :: thetasat ! volumetric soil water at saturation
    real(r8), intent(in) :: perc     ! percolation water flux thourgh the bottom of volatilization layer, m/s
    real(r8), intent(in) :: runoff   ! surface runoff, m/s
    real(r8), intent(in) :: bsw
    real(r8), intent(in) :: kads     ! dimensionless distribution coefficient, kads = [TAN (s)] / [TAN (aq)]
    !real(r8), intent(in) :: dt       ! timestep

    real(r8) :: water_tot, cnc, air, depth_soilsat, diffusivity_water, diffusivity_satsoil, halfwater, insoil, r1, dz2, inwater
    real(r8) :: r2, volat_rate, kno3, knh3, depth_lower, fract_nh4, r2a, r2b, g3, gdown, rsld, rkl, rkg


    water_tot = water(1) + water(2)

    air = thetasat - theta
    ! depth of the saturated soil layer below the surface pool
    depth_soilsat = water(2) / air 

    cnc = mtan / water_tot

    fluxes(iflx_roff) = cnc * runoff
    fluxes(iflx_soilq) = cnc * perc


    diffusivity_water = 9.8e-10_r8 * 1.03_r8 ** (tg - 273.0_r8)
    diffusivity_satsoil = eval_diffusivity_liq_mq(thetasat, thetasat, tg) * thetasat
    !diffusivity_satsoil = eval_diffusivity_liq_m03(thetasat, thetasat, tg, bsw) * thetasat
    
    halfwater = 0.5_r8 * water_tot

    ! Calculate the internal resistance r1 of the slurry/soil layer by integrating the
    ! diffusivity for distance that covers half of the slurry water.
    if (water(1) < halfwater) then
       ! contribution from both pool and the saturated soil.
       insoil = (halfwater - water(1)) / thetasat
       inwater = water(1)
    else
       ! pool only
       inwater = halfwater
       insoil = 0.0
       !r1 = halfwater / diffusivity_water
    end if

    r1 = inwater / diffusivity_water + insoil /diffusivity_satsoil

    depth_lower = max(soildepth_reservoir, depth_soilsat*1.5)

    call partition_tan(tg, Hconc, 1.0_r8, 0.0_r8, 0.0_r8, knh3, fract_nh4=fract_nh4)
    volat_rate = &
         knh3/(-ratm*kads*theta + ratm*kads + ratm*thetasat - r1*kads*knh3*theta + r1*kads*knh3 + r1*knh3*thetasat)
    
    !volat_rate = 1.0_r8 / (r1 + ratm / knh3) ! conductance from aqueous TAN in slurry to NH3 in atmosphere
    fluxes(iflx_air) = volat_rate*cnc
        
    ! lower soil resistance consists of liquid diffusion slurry, in saturated layer, and
    ! parallel liquid/gas diffusion below the saturated layer. 
    r2a = (water(1)-inwater) / diffusivity_water
    r2b = (depth_soilsat-insoil) / diffusivity_satsoil
    dz2 = depth_lower-depth_soilsat

    Rkl = dz2 / (eval_diffusivity_liq_mq(theta, thetasat, tg)*theta)
    Rkg = dz2 / (eval_diffusivity_gas_mq(theta, thetasat, tg)*(thetasat-theta))
    Rsld = r2a + r2b

    gdown = -(Rkg + Rkl*knh3)/((Rkg*(Rkl + Rsld) + Rkl*Rsld*knh3)*(kads*(theta - 1) - thetasat))

    ! conductance = inverse resistance
    !g3 = (eval_diffusivity_liq_mq(theta, thetasat, tg)*theta &
    !     + knh3*eval_diffusivity_gas_mq(theta, thetasat, tg)*(thetasat-theta)) &
    !     / dz2
    !r2 = r2a + r2b + 1.0/g3

    fluxes(iflx_soild) = cnc * gdown
    
    ! nitrification
    kno3 = eval_no3prod_v2(thetasat, thetasat, tg)
    fluxes(iflx_no3) = kno3 * mtan * fract_nh4

    !fluxes(3:) = 0
    
  end subroutine eval_fluxes_slurry

  subroutine eval_fluxes_soilroff_ads(mtan, water_manure, Hconc, tg, ratm, theta, thetasat, perc, &
       & runoff, bsw, kads_nh4, soildepth, fluxes, substance, status)
    !
    ! Evaluate nitrogen fluxes from a soil layer. Use for all cases except the partly
    ! infiltrated slurry (above). Fluxes can be evaluated either for urea or TAN: for
    ! urea, only the aqueous phase fluxes are evaluated and nitrification is set to zero.
    ! 
    implicit none
    real(r8), intent(in) :: mtan ! TAN (=NH4 (aq) + NH3 (g) + NH3 (aq)), mass units / m2
    real(r8), intent(in) :: water_manure ! water in the soil pool *in addition to* background soil water
    real(r8), intent(out) :: fluxes(5)   ! nitrogen fluxes, mass units / m2 / s, see top of module
    real(r8), intent(in) :: Hconc        ! Hydrogen ion concentration, mol/l
    real(r8), intent(in) :: tg           ! soil temperature, K
    real(r8), intent(in) :: ratm         ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta        ! volumetric soil water in "clean" soil, m/m
    real(r8), intent(in) :: thetasat     ! volumetric soil water at saturation
    real(r8), intent(in) :: perc         ! downwards water percolation rate at the bottom of layer, m/s, > 0
    real(r8), intent(in) :: runoff       ! runoff water flux, m / s
    real(r8), intent(in) :: bsw          ! b in the soilwater retention curve; needed if the Moldrup 2003 diffusivities are used.
    real(r8), intent(in) :: kads_nh4     ! distribution coefficient kads = [TAN (s)] / [TAN (aq)]. Unit m3(water) / m3(soil).
    real(r8), intent(in) :: soildepth    ! thickness of the volatlization layer
    integer, intent(in) :: substance     ! subst_tan or subst_urea.
    integer, intent(out) :: status       ! error flag
    
    real(r8) :: water_tot, cnc, air, henry_eff, dsl, dsg, dstot, dz2, no3_rate, volat_rate, theta_tot, beta
    real(r8) :: cnc_srfg, cnc_srfaq, cnc_soilaq, cnc_soilg, dz, rsl, rsg
    real(r8) :: fract_gas, fract_nh3aq, fract_nh4, fract_aq, volatility

    real(r8) :: kads ! distribution coefficient, unitless ((g NH4 adsorbed / m3 soil solid) / (g NH4 dissolved / m3 soil water))

    water_tot = water_manure + theta*soildepth
    if (water_tot < 1e-9) then
       fluxes = 0.0
       return
    end if
    
    theta_tot = water_tot / soildepth
    if (theta_tot > thetasat) then
       status = err_bad_theta
       return
    end if

    cnc = mtan / soildepth
    air = thetasat - theta_tot

    dz = 0.5*soildepth 
    dsl = eval_diffusivity_liq_mq(theta_tot, thetasat, Tg)
    dsg = eval_diffusivity_gas_mq(theta_tot, thetasat, Tg)

    rsg = dz / (dsg*air)
    rsl = dz / (dsl*theta_tot)

    if (substance == subst_tan) then
       kads = kads_nh4 
       call partition_tan(tg, Hconc, theta_tot, air, kads, volatility, fract_nh4=fract_nh4)
       no3_rate = eval_no3prod_v2(theta_tot, thetasat, tg)
    else if (substance == subst_urea) then
       volatility = 0.0
       no3_rate = 0.0_r8
       kads = 0.0
    else
       status = err_bad_subst
       return
    end if
    
    call get_srf_cnc(volatility, cnc, 0.0_r8, rsg, rsl, Ratm, runoff, theta_tot, air, cnc_srfg, cnc_srfaq)

    fluxes(iflx_air) = cnc_srfg / ratm
    fluxes(iflx_roff) = runoff * cnc_srfaq

    dz2 = soildepth_reservoir - 0.5 * soildepth
    cnc_soilaq = cnc / (theta_tot + air*volatility + (1-theta_tot-air)*kads)
    cnc_soilg = cnc_soilaq*volatility

    fluxes(iflx_soild) = (cnc_soilaq * dsl*theta_tot) / dz2 + (cnc_soilg * dsg*air) / dz2
    fluxes(iflx_soild) = fluxes(iflx_soild) + mtan / (24*3600.0*365)
    fluxes(iflx_no3) = mtan * no3_rate
    fluxes(iflx_soilq) = cnc_soilaq * perc

    status = 0
    
  contains
    
    subroutine get_srf_cnc(knh3, xs, xag, Rsg, Rsl, Rag, qr, theta, air, cnc_gas, cnc_aq)
      ! automatically generated by compost.py
      real(r8), intent(in) :: knh3 ! volatility
      real(r8), intent(in) :: xag  ! cnc_atm_gas
      real(r8), intent(in) :: qr   ! runoff m/s
      real(r8), intent(in) :: Rag  ! Ratm
      real(r8), intent(in) :: xs   ! mass / m3 soil
      real(r8), intent(in) :: Rsg, Rsl, theta, air

      real(r8), intent(out) :: cnc_gas, cnc_aq

      real(r8) :: x0, x1, x2, x3, x4, x5, x6

      x0 = air*knh3
      x1 = air*kads
      x2 = kads*theta
      x3 = Rag + Rsg
      x4 = Rsl*knh3*x3
      x5 = Rag*Rsg*(Rsl*qr + 1)
      x6 = (Rag*xs*(Rsg + Rsl*knh3) + Rsg*Rsl*xag*(kads + theta + x0 - x1 - x2)) &
           /(Rsl*air*knh3**2*x3 + kads*x4 + kads*x5 + theta*x4 + theta*x5 + x0*x5 - x1*x4 - x1*x5 - x2*x4 - x2*x5)
      cnc_gas = knh3*x6
      cnc_aq = x6

    end subroutine get_srf_cnc

  end subroutine eval_fluxes_soilroff_ads

  subroutine partition_to_layer(water, theta, thetasat, soildepth, fraction_in, fraction_down, fraction_runoff)
    ! Evaluate the fraction of water volume that can be accommodated (before saturation)
    ! by soil layer with current water content theta.
    implicit none
    real(r8), intent(in) :: water ! water to be added to the layer, m
    real(r8), intent(in) :: theta, thetasat ! vol. soil water, current and saturation, m/m
    real(r8), intent(in) :: soildepth ! thickness of the layer, m
    real(r8), intent(out) :: fraction_in ! fraction that fits in
    real(r8), intent(out) :: fraction_down ! fraction excess
    real(r8), intent(out) :: fraction_runoff ! = 1 if theta is > 99% saturation otherwise 0

    real(r8) :: watvol, watvol_sat, vol_to_layer, vol_to_deeper, vol_avail
    
    watvol = theta*soildepth
    watvol_sat = thetasat*soildepth

    if (watvol < 0.99*watvol_sat) then
       vol_avail = watvol_sat - watvol
       vol_to_layer = min(water, vol_avail)
       vol_to_deeper = water - vol_to_layer
       fraction_down = vol_to_deeper / water
       fraction_runoff = 0.0
       fraction_in = 1.0_r8 - fraction_down
    else
       fraction_down = 0.0
       fraction_in = 0.0
       fraction_runoff = 1.0_r8
    end if
    
  end subroutine partition_to_layer

  subroutine age_pools_soil(ndep, dt, pools, mtan, garbage)
    implicit none
    real(r8), intent(in) :: ndep ! flux of TAN input, gN/m2/s
    real(r8), intent(inout) :: mtan(:) ! TAN pools for each age range. gN/m2
    real(r8), intent(in) :: dt ! timestep, s
    real(r8), intent(in) :: pools(:) ! age spans covered by each bin, seconds. size(mtan) = size(pools)
    real(r8), intent(out) :: garbage ! TAN removed from the oldest pool. gN/m2

    real(r8) :: flux_out(size(mtan))

    flux_out = mtan / pools

    ! new nitrogen
    mtan(1) = mtan(1) + ndep * dt
    ! transfer nitrogen from fresh to old pools
    mtan = mtan - flux_out * dt
    if (size(mtan) > 1) mtan(2:) = mtan(2:) + flux_out(:size(mtan)-1) * dt
    ! provided that the oldest pool has wide enough age range, the amount transferred out
    ! should be small.
    garbage = flux_out(size(mtan)) * dt
    
  end subroutine age_pools_soil

  subroutine age_pools_slurry(ndep, dt, water_slurry, tan_slurry, tan_soil, pools, garbage)
    implicit none
    real(r8), intent(in) :: ndep ! flux of TAN input, gN/m2/s
    real(r8), intent(in) :: dt ! timestep, s
    ! water in slurry pool, on surface (1) and below surface (2) not including the background water content (theta)
    real(r8), intent(in) :: water_slurry(2)
    real(r8), intent(inout) :: tan_slurry  ! TAN in slurry pool, gN/m2
    real(r8), intent(inout) :: tan_soil(:) ! TAN in soil pools, gN/m2
    real(r8), intent(in) :: pools(:) ! age spans covered by each pool, including the S0 surface slurry. seconds.
    real(r8), intent(out) :: garbage       ! garbage TAN, see above

    real(r8) :: fract_slurry_soil, flux_to_soilpools
    
    fract_slurry_soil = water_slurry(2) / sum(water_slurry)
    flux_to_soilpools = tan_slurry * fract_slurry_soil / pools(1)
    tan_slurry = tan_slurry + (ndep - flux_to_soilpools) * dt

    call age_pools_soil(flux_to_soilpools, dt, pools(2:), tan_soil, garbage)
    
  end subroutine age_pools_slurry

  subroutine update_4pool(tg, ratm, theta, thetasat, precip, evap, qbot, watertend, runoff, tandep, tanprod, bsw, &
       depth_slurry, poolranges, tanpools, Hconc, fluxes, garbage, dt, status)
    !
    ! Evaluate fluxes and integrate states for a 4-stage slurry model with first pool
    ! representing uninfiltrated slurry.
    !
    implicit none
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: ratm ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta ! volumetric soil water in soil column (unaffected by slurry)
    real(r8), intent(in) :: thetasat ! vol. soil water at saturation
    real(r8), intent(in) :: precip   ! precipitation, m/s
    real(r8), intent(in) :: evap   ! ground evaporation, m/s
    real(r8), intent(in) :: qbot   ! specific humidity (kg/kg) at lowest atmospheric model level
    real(r8), intent(in) :: watertend ! time derivative of theta
    real(r8), intent(in) :: runoff    ! surface runoff flux, m/s
    real(r8), intent(in) :: tandep    ! TAN input flux, gN/m2/s
    real(r8), intent(in) :: tanprod   ! TAN produced in the column, added to aged TAN pool
    !real(r8), intent(in) :: infiltr_slurry ! Slurry infiltration rate, m/s
    real(r8), intent(in) :: depth_slurry ! Initial slurry depth, m
    real(r8), intent(in) :: bsw
    real(r8), intent(in) :: poolranges(4)  ! age ranges of TAN pools S0, S1, S2, sec. Slurry infiltration time is inferred from S0. 
    real(r8), intent(inout) :: tanpools(4) ! TAN pools gN/m2
    real(r8), intent(out) :: fluxes(5,4) ! TAN fluxes, gN/m2/s (type of flux, pool)
    real(r8), intent(in) :: Hconc(4) ! H+ concentration
    real(r8), intent(out) :: garbage     ! over-aged TAN occurring during the step, gN/m.
    real(r8), intent(in) :: dt ! timestep, sec, >0
    integer, intent(out) :: status ! return status, 0 = good
    
    real(r8) :: infiltr_slurry, infiltrated, percolated, evap_slurry, water_slurry(2), perc_slurry_mean, waterloss
    real(r8) :: percolation, water_soil, age_prev, water_in_layer, tanpools_old(4)
    integer :: indpl
    
    real(r8), parameter :: dz_layer = 0.02 ! thickness of the volatilization layer, m
    real(r8), parameter :: kads = 0.0_r8   ! distriution coefficient kads = [TAN (s)] / [TAN (aq)], dimensionless

    ! H+ concentration in each pool 
    !real(r8), parameter :: Hconc(4) = (/10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8), 10.0_r8**(-7_r8)/)

    if (theta > thetasat) then
       status = err_bad_theta
       return
    end if

    tanpools_old = tanpools

    ! Pool S0
    !
    evap_slurry = get_evap_pool(tg, ratm, qbot)
    infiltr_slurry = max(depth_slurry / poolranges(1), precip)
    infiltrated = depth_slurry * infiltr_slurry / (infiltr_slurry + evap_slurry)
    ! Slurry water (in addition to soil water, theta) on surface and in soil. Represents
    ! mean over pool S0.
    water_slurry = (/0.5*depth_slurry, 0.5*infiltrated/)
    ! The excess water assumed to have percolated down from the volat. layer.
    percolated = max(infiltrated - dz_layer*(thetasat-theta), 0.0)
    ! Percolation rate out of volat layer, average over the pool S0.
    perc_slurry_mean = percolated / poolranges(1)

    call eval_fluxes_slurry(water_slurry, tanpools(1), Hconc(1), tg, ratm, theta, thetasat, perc_slurry_mean, &
         runoff, bsw, kads, fluxes(:,1))

    if (any(isnan(fluxes))) then
       status = err_nan * 10
    end if

    call update_pools(tanpools(1:1), fluxes(1:5,1:1), dt, 1, 5)

    if (any(tanpools < -1e-15)) then
       status = err_negative_tan
       return
    end if
    
    ! Pool aging & input
    !
    call age_pools_slurry(tandep, dt, water_slurry, tanpools(1), tanpools(2:), poolranges, garbage)
    ! TAN produced (mineralization) goes to directly the old TAN pool. 
    tanpools(4) = tanpools(4) + tanprod*dt

    ! Soil bins S1 and S2
    !
    age_prev = 0 ! for water evaluations, consider beginning of S1 as the starting point
    water_in_layer = infiltrated - percolated ! water in layer just after slurry has infiltrated
    do indpl = 2, 4
       ! water content lost during the aging
       waterloss = water_in_layer * (waterfunction(age_prev) - waterfunction(age_prev+poolranges(indpl)))
       percolation = eval_perc(waterloss, evap, precip, watertend, poolranges(indpl))

       ! water content at the mean age of the pool
       water_soil = water_in_layer * waterfunction(age_prev + 0.5*poolranges(indpl))
       call eval_fluxes_soilroff_ads(tanpools(indpl), water_soil, Hconc(indpl), tg, &
            & ratm, theta, thetasat, percolation, runoff, bsw, kads, &
            & dz_layer, fluxes(:,indpl), subst_tan, status)

       if (status /= 0) return
       age_prev = age_prev + poolranges(indpl)
    end do
    
    call update_pools(tanpools(2:), fluxes(1:5,2:), dt, 3, 5)
    
    if (any(tanpools < -1e-15)) then
       !if (any(tanpools < -1e-3)) then
       status = err_negative_tan * 10
       return
       !end if
    end if

    if (any(isnan(fluxes))) then
       status = err_nan * 100
       return
    end if

    if (abs(sum(tanpools - tanpools_old) - (-sum(fluxes) + tandep + tanprod)*dt + garbage) > max(sum(tanpools_old)*1e-2, 1e-4)) then
       status = err_balance_tan
       return
    end if
    
    status = 0
    
  end subroutine update_4pool
  
  subroutine update_npool(tg, ratm, theta, thetasat, precip, evap, qbot, watertend, runoff, tandep, tanprod, &
       water_init, bsw, poolranges, Hconc, dz_layer, tanpools, fluxes, garbage, dt, status, numpools)
    !
    ! Evaluate fluxes and update soil TAN pools for a model with arbitrary number of pools
    ! divided by age and pH. For slurry use update_4pool.
    ! 
    implicit none
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: ratm ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta ! volumetric soil water in soil column (unaffected by slurry)
    real(r8), intent(in) :: thetasat ! vol. soil water at saturation
    real(r8), intent(in) :: precip   ! precipitation, m/s
    real(r8), intent(in) :: evap   ! ground evaporation, m/s
    real(r8), intent(in) :: qbot   ! specific humidity (kg/kg) at lowest atmospheric model level
    real(r8), intent(in) :: watertend ! time derivative of theta*dz
    real(r8), intent(in) :: runoff    ! surface runoff flux, m/s
    real(r8), intent(in) :: tandep    ! TAN input flux, gN/m2/s
    real(r8), intent(in) :: tanprod(numpools)   ! flux of TAN produced (from urea/organic n) in the column
    real(r8), intent(in) :: water_init ! Initial water volume in the affected patch, m
    real(r8), intent(in) :: bsw
    real(r8), intent(in) :: poolranges(numpools)  ! age ranges of TAN pools (npools)
    real(r8), intent(in) :: Hconc(numpools)       ! H+ concentration, mol/l (npools)
    real(r8), intent(in) :: dz_layer ! thickness of the volatilization layer, m
    real(r8), intent(inout) :: tanpools(numpools) ! TAN pools gN/m2 (npools)
    real(r8), intent(out) :: fluxes(5,numpools) ! TAN fluxes, gN/m2/s (type of flux, pool)
    real(r8), intent(out) :: garbage     ! "over-aged" TAN produced during the step, gN/m.
    real(r8), intent(in) :: dt ! timestep, sec, >0
    integer, intent(out) :: status ! 0 == OK
    integer, intent(in) :: numpools
    
    real(r8) :: fraction_layer, fraction_reservoir, fraction_runoff, waterloss, direct_runoff
    real(r8) :: percolation, water_soil, age_prev, tandep_remaining, direct_percolation, water_into_layer
    real(r8) :: tanpools_old(size(tanpools)), imbalance
    integer :: indpl
    
    real(r8), parameter :: water_relax_t = 24*3600.0_r8
    real(r8), parameter :: kads = 0.0_r8
    logical :: fixed

    tanpools_old = tanpools
    
    if (theta > thetasat) then
       status = err_bad_theta
       return
    end if
    
    ! Initial water excess goes to runoff if the surface is close to saturation, otherwise to the soil.
    !
    call partition_to_layer(water_init, theta, thetasat, dz_layer, fraction_layer, fraction_reservoir, fraction_runoff)
    direct_runoff = fraction_runoff * tandep
    direct_percolation = fraction_reservoir * tandep
    tandep_remaining = tandep - direct_runoff - direct_percolation
    water_into_layer = water_init * (1.0_r8 - fraction_reservoir - fraction_runoff)
    if (tandep_remaining < -1e-15) then
       print *, tandep, direct_runoff, direct_percolation
       status = err_negative_tan + 10
       return
    end if
    if (water_into_layer/dz_layer + theta > thetasat+1e-5) then
       status = err_bad_theta
       return
    end if

    if(any(isnan(tanpools))) then
       status = err_nan 
       return
    end if
    
    ! Pool aging & TAN input
    !
    call age_pools_soil(tandep_remaining, dt, poolranges, tanpools, garbage)
    ! TAN produced (mineralization) goes to directly the old TAN pool. 
    
    if (any(tanpools < 0)) then
       if (any(tanpools < -1e-15)) then
          print *, '<0 2', tanpools_old, 'new', tanpools, 'fx', &
               sum(fluxes(:,1)) * dt, sum(fluxes(:,2)) * dt, sum(fluxes(:,3)) * dt
          status = err_negative_tan + 10000
          return
       else
          where(tanpools<0) tanpools = 0.0
       end if
    end if

    imbalance = abs((sum(tanpools) - sum(tanpools_old)) - ((tandep_remaining)*dt+garbage))
    if (imbalance > max(1e-14, 0.001*sum(tanpools_old))) then
       print *, imbalance, 'old', tanpools_old, 'new', tanpools, 'g', garbage, tandep_remaining, &
            'diff', sum(tanpools_old-tanpools), &
            'depprod', tandep+sum(tanprod), garbage
       status = err_balance_tan*10
       return
    end if
    
    age_prev = 0 ! for water evaluations, consider beginning of S1 as the starting point
    do indpl = 1, size(tanpools)
       ! water content lost during the aging
       waterloss = water_into_layer * (waterfunction(age_prev) - waterfunction(age_prev+poolranges(indpl)))
       percolation = eval_perc(waterloss, evap, precip, watertend, poolranges(indpl))
       ! water content at the middle of the age range
       water_soil = water_into_layer * waterfunction(age_prev + 0.5*poolranges(indpl))
       call eval_fluxes_soilroff_ads(tanpools(indpl), water_soil, Hconc(indpl), tg, &
            & ratm, theta, thetasat, percolation, runoff, bsw, kads, &
            & dz_layer, fluxes(:,indpl), subst_tan, status)
       if (status /= 0) then
          return
       end if
       age_prev = age_prev + poolranges(indpl)
    end do
    
    call update_pools(tanpools, fluxes, dt, numpools, 5, fixed)
    !print *, 'pools just after', tanpools
    tanpools = tanpools + tanprod*dt
    if(any(isnan(tanpools))) then
       status = err_nan+100
       return
    end if
    
    if (any(isnan(fluxes))) then
       status = err_nan + 1000
    end if
    
    if (any(tanpools < -1e-15)) then
       status = err_negative_tan + 1000
       print *, '<0 2', tanpools_old, tanpools, sum(fluxes(:,1)) * dt, sum(fluxes(:,2)) * dt
       
       return
    end if
    
    
    if (abs(sum(tanpools - tanpools_old) + (sum(fluxes)-tandep_remaining-sum(tanprod))*dt + garbage) &
         > max(sum(tanpools_old)*1e-2, 1d-2)) then
       print *, tanpools, tanpools_old, 'fx', fluxes*dt, 'dp', tandep_remaining*dt, tanprod*dt, 'g', garbage, &
            'ib', sum(tanpools-tanpools_old), (sum(fluxes)-tandep_remaining-sum(tanprod))*dt + garbage, 'fix', fixed
       status = err_balance_tan
       return
    end if

    ! Add the "direct" fluxes to the fluxes of the first pool 
    fluxes(iflx_roff, 1) = fluxes(iflx_roff, 1) + direct_runoff
    fluxes(iflx_soilq, 1) = fluxes(iflx_soilq, 1) + direct_percolation
    
    if (any(fluxes < -1e-6)) then
       print *, fluxes
       status = err_negative_flux
       return
    end if
    
    status = 0
    
  end subroutine update_npool

  subroutine update_pools(tanpools, fluxes, dt, np, nf, fixed)
    ! Update tan pools using the fluxes and an ad-hoc scheme against negative TAN masses.
    implicit none
    real(r8), intent(inout) :: tanpools(np), fluxes(nf,np)
    real(r8), intent(in) :: dt
    integer, intent(in) :: np, nf
    logical, intent(out), optional :: fixed

    integer :: ip
    real(r8) :: sumflux, ff
    logical :: fixed_

    fixed_ = .false.
    do ip = 1, np
       sumflux = sum(fluxes(:,ip))*dt
       if (sumflux > tanpools(ip)) then
          if (sumflux > 1e-15) then
             fixed_ = .true.
             ff = tanpools(ip) / sumflux
             fluxes(:,ip) = fluxes(:,ip) * ff
             sumflux = tanpools(ip)
             sumflux = sum(fluxes(:,ip))*dt
          else
             sumflux = 0.0
          end if
       end if
       tanpools(ip) = tanpools(ip) - sumflux
    end do
    if (present(fixed)) fixed = fixed_
    
  end subroutine update_pools
  
  function get_evap_pool(tg, ratm, qbot) result(evap)
    ! Evaluate evaporation rate for surface water given spcific humidity at the reference
    ! height.
    implicit none
    real(r8), intent(in) :: tg, ratm, qbot
    real(r8) :: evap ! m/s

    real(r8) :: es, esdt, qs, qsdt, dens, flux
    real(r8), parameter :: press = 101300.0_r8
    
    call qsat(tg, press, es, esdt, qs, qsdt)
    if (qbot > qs) then
       evap = 0
       return
    end if

    dens = press / (SHR_CONST_RDAIR * tg)
    flux = dens * (qs - qbot) / ratm ! kg/s/m2 == mm/s
    evap = flux*1e-3
    
  end function get_evap_pool

  ! Waterfunction gives the relaxation of the moisture perturbation normalized between 0
  ! and 1. Either exponential or step.
  
  function waterfunction_exp(pool_age) result(water)
    implicit none
    real(r8), intent(in) :: pool_age ! sec
    real(r8) :: water

    water = exp(-pool_age / water_relax_t)
  end function waterfunction_exp

  function waterfunction(pool_age) result(water)
    implicit none
    real(r8), intent(in) :: pool_age ! sec
    real(r8) :: water

    if (pool_age > water_relax_t) then
       water = 0.0_r8
    else
       water = 1.0_r8 - pool_age / water_relax_t
    end if
    
  end function waterfunction
  
  function eval_perc(waterloss, evap, precip, watertend, dt) result(rate)
    ! Evaluate the downwards water flux at the layer bottom given the infiltration and
    ! evaporation fluxes.
    implicit none
    real(r8), intent(in) :: waterloss ! total water loss during dt, m
    real(r8), intent(in) :: evap ! average evaporation rate, m/s
    real(r8), intent(in) :: precip ! average infiltration rate, m/s
    real(r8), intent(in) :: watertend ! background water tendency, m/s
    real(r8), intent(in) :: dt ! timespan, s

    real(r8) :: rate ! percolation rate, m/s
    real(r8) :: perc_base, perc_adj
    
    perc_base = waterloss / dt
    perc_adj = perc_base + precip - evap - watertend
    rate = max(perc_adj, 0.0)

  end function eval_perc

  subroutine eval_fluxes_storage(nitr_input, barntype, tempr_outside, windspeed, fract_direct, &
       volat_coef_barns, volat_coef_stores, &
       tan_fract_excr, fluxes_nitr, fluxes_tan, status)
    !
    ! Evaluate nitrogen fluxes in animal housings and storage. Only volatilization losses
    ! are assumed. The volatilization fluxes are assumed to depend linearly on the TAN
    ! fluxes entering the housings or storage. The base coefficients are given as
    ! arguments and adjusted according the model of Gyldenkaerne et al.
    ! 
    implicit none
    real(r8), intent(in) :: nitr_input ! total nitrogen excreted by animals in housings
    character(len=*), intent(in) :: barntype ! "closed" (pigs, poultry) or "open" (others)
    real(r8), intent(in) :: tempr_outside ! K
    real(r8), intent(in) :: windspeed ! m/s
    real(r8), intent(in) :: fract_direct ! fraction of manure N applied before storage
    real(r8), intent(in) :: volat_coef_barns, volat_coef_stores ! normalization coefficients, unitless
    real(r8), intent(in) :: tan_fract_excr ! fraction of NH4 nitrogen in excreted N
    real(r8), intent(out) :: fluxes_nitr(4), fluxes_tan(4) ! nitrogen and TAN fluxes, gN/s
                                                           ! (/m2). See top of module for
                                                           ! indices.
    integer, intent(out) :: status ! see top of the module.

    ! parameters for the Gyldenkaerne et al. parameterization
    real(r8), parameter :: Tfloor_barns = 4.0_r8, Tfloor_stores = 1.0_r8
    real(r8), parameter :: Tmin_barns = 0.01_r8
    real(r8), parameter :: Tmax_barns = 12.5_r8
    real(r8), parameter :: tempr_D = 3.0_r8
    real(r8), parameter :: Trec = 21.0_r8
    real(r8), parameter :: Vmin_barns = 0.2_r8
    !real(r8), parameter :: Vmax_barns = 0.228_r8
    real(r8), parameter :: pA = 0.89_r8, pB = 0.26_r8  
    real(r8), parameter :: DTlow = 0.5_r8, DThigh = 1.0_r8
    real(r8) :: Vmax_barns ! depends on barntype
    
    real(r8) :: flux_avail, flux_avail_tan, tempr_stores, tempr_barns, vent_barns, flux_direct, flux_direct_tan, &
         & flux_barn, flux_store, tempr_C
    
    fluxes_nitr = 0.0_r8
    fluxes_tan = 0.0_r8
    tempr_C = tempr_outside - 273
    
    select case(barntype)
    case ('open')
       Vmax_barns = 0.228_r8
    tempr_barns = max(tempr_C+tempr_D, Tfloor_barns)
    case ('closed')
       Vmax_barns = 0.40_r8
       if (Trec + DTlow * (tempr_C - Tmin_barns) < Tmin_barns) then
          tempr_barns = Tmin_barns
       else if (tempr_C < Tmin_barns) then
          tempr_barns = Trec + DTlow * (tempr_C - Tmin_barns)
       else if (tempr_C > Tmax_barns) then
          tempr_barns = Trec + DThigh * (tempr_C - Tmax_barns)
       else
          tempr_barns = Trec
       end if
    case default
       status = err_bad_type
       return
    end select

    if (tempr_C < Tmin_barns) then
       vent_barns = Vmin_barns
    else if (tempr_C > Tmax_barns) then
       vent_barns = Vmax_barns
    else
       vent_barns = Vmin_barns + tempr_C * (Vmax_barns-Vmin_barns) / (Tmax_barns - Tmin_barns)
    end if

    flux_avail = nitr_input
    flux_avail_tan = nitr_input * tan_fract_excr

    if (flux_avail < -1e-15 .or. flux_avail_tan < -1e-15) then
       status = err_negative_flux*1000
       return
    end if
    
    flux_barn = flux_avail_tan * volat_coef_barns * tempr_barns**pA * vent_barns**pB
    flux_barn = min(flux_avail_tan, flux_barn) ! hopefully uncommon
    fluxes_tan(iflx_air_barns) = flux_barn
    fluxes_nitr(iflx_air_barns) = flux_barn
    
    flux_avail = flux_avail - flux_barn
    flux_avail_tan = flux_avail_tan - flux_barn

    if (flux_avail < 0 .or. flux_avail_tan < 0) then
       !print *, flux_avail_tan, flux_avail, flux_barn, tempr_barns, vent_barns, tempr_C
       status = err_negative_flux*10000
       return
    end if
    
    flux_direct = fract_direct * flux_avail
    flux_avail = flux_avail - flux_direct
    flux_direct_tan = flux_avail_tan * fract_direct
    flux_avail_tan = flux_avail_tan - flux_direct_tan

    fluxes_tan(iflx_appl) = flux_direct_tan
    fluxes_nitr(iflx_appl) = flux_direct

    tempr_stores = max(Tfloor_stores, tempr_C)
    flux_store = flux_avail_tan &
         & * volat_coef_stores * tempr_stores**pA * windspeed**pB
    flux_store = min(flux_avail_tan, flux_store)

    fluxes_tan(iflx_air_stores) = flux_store
    fluxes_nitr(iflx_air_stores) = flux_store
    
    flux_avail = flux_avail - flux_store
    flux_avail_tan = flux_avail_tan - flux_store
    if (flux_avail < 0) then
       status = err_negative_flux*10
       return
    end if

    fluxes_nitr(iflx_to_store) = flux_avail
    fluxes_tan(iflx_to_store) = flux_avail_tan
    
    if (abs(sum(fluxes_nitr) - nitr_input) > 1e-5*nitr_input) then
       status = err_balance_nitr
       return
    end if

    if (abs(sum(fluxes_tan) - nitr_input*tan_fract_excr) > 1e-5*nitr_input) then
       status = err_balance_tan
       return
    end if

    if (any(fluxes_nitr < 0) .or. any(fluxes_tan < 0)) then
       status = err_negative_flux*100
       return
    end if
    
    status = 0
    
  end subroutine eval_fluxes_storage

  subroutine update_org_n(flux_input, tg, soilpsi, pools, dt, tanprod, soilflux)
    !
    ! Evaluate the decomposition/mineralization N fluxes from the available, resistant and
    ! unavailable N fractions, and update the organic N pools. In addition, evaluate the
    ! flux of organic N into the soil pools according to a fixed time constant set below.
    implicit none
    real(r8), intent(in) :: flux_input(3) ! organic N entering the pools. gN/m2/s. For
                                          ! indices see at top of the module.
    real(r8), intent(in) :: tg          ! ground temperature, K
    real(r8), intent(in) :: soilpsi     ! soil water potential (MPa)
    real(r8), intent(inout) :: pools(3) ! organic N pools
    real(r8), intent(in) :: dt       ! timestep, sec
    real(r8), intent(out) :: tanprod(3)      ! Flux of TAN formed, both pools
    real(r8), intent(out) :: soilflux  ! Flux of organic nitrogen to soil

    real(r8) :: rate_res, rate_avail, TR, rmoist, psi
    real(r8), parameter :: ka1 = 8.94e-7_r8, ka2 = 6.38e-8 ! 1/s
    real(r8), parameter :: tr1 = 0.0106, tr2 = 0.12979 
    real(r8), parameter :: org_to_soil_time = 365*24*3600.0_r8
    real(r8), parameter :: minpsi = -2.5_r8, maxpsi=-0.002_r8
    real(r8) :: soilfluxes(3)
    
    TR = tr1 * exp(tr2 * (tg-273.15_r8))

    ! The moisture scaling taken from CLM5 litter decomposition scheme:
    psi = min(soilpsi, maxpsi)
    ! decomp only if soilpsi is higher than minpsi
    if (psi > minpsi) then
       rmoist = log(minpsi/psi) / log(minpsi/maxpsi)
    else
       rmoist = 0.0
    end if
     
    tanprod(ind_avail) = ka1 * TR * pools(ind_avail)*rmoist
    tanprod(ind_resist) = ka2 * TR * pools(ind_resist)*rmoist
    tanprod(ind_unavail) = 0.0
    soilfluxes = pools * 1.0_r8 / org_to_soil_time

    pools = pools + (flux_input - tanprod - soilfluxes) * dt
    soilflux = sum(soilfluxes)
    
  end subroutine update_org_n

  subroutine update_urea(tg, theta, thetasat, precip, evap, watertend, runoff, &
       ndep, bsw, pools, fluxes, garbage, ranges, dt, status, numpools)
    !
    ! Evaluate fluxes and update the urea pools. The procedure is similar to updating the
    ! soil TAN pools, but NO3 and volatilization fluxes do not occur.
    ! 
    implicit none
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: theta ! volumetric soil water in soil column (background)
    real(r8), intent(in) :: thetasat ! vol. soil water at saturation
    real(r8), intent(in) :: precip   ! precipitation, m/s
    real(r8), intent(in) :: evap   ! ground evaporation, m/s
    real(r8), intent(in) :: watertend ! time derivative of theta*dz
    real(r8), intent(in) :: runoff    ! surface runoff flux, m/s
    real(r8), intent(in) :: ndep ! nitrogen input, mass unit / s
    real(r8), intent(in) :: bsw  ! b in the soil water retention curve
    real(r8), intent(inout) :: pools(numpools) ! nitrogen pools mass / m2
    real(r8), intent(out) :: fluxes(6, numpools) ! one extra for the to_tan flux
    real(r8), intent(in) :: ranges(numpools) ! pool age extents, s
    real(r8), intent(out) :: garbage ! nitrogen in patches aged beyond the oldest pool. mass / m2 
    real(r8), intent(in) :: dt ! time step, s
    integer, intent(out) :: status ! see top of module
    integer, intent(in) :: numpools 

    real(r8), parameter :: rate = 4.83e-6 ! urea decomposition, 1/s
    real(r8), parameter :: missing = 1e36 ! for the parameters not needed for urea fluxes
    real(r8), parameter :: dz_layer = 0.02 ! thickness of the volatilization layer, m
    
    real(r8) :: age_prev, percolation, old_total, balance
    integer :: indpl

    old_total = sum(pools)
    
    call age_pools_soil(ndep, dt, ranges, pools, garbage)

    age_prev = 0 
    do indpl = 1, numpools
       percolation = eval_perc(0.0_r8, evap, precip, watertend, ranges(indpl))
       ! Hconc and Ratm are missing since they do not affect urea.
       call eval_fluxes_soilroff_ads(pools(indpl), 0.0_r8, missing, tg, &
            & missing, theta, thetasat, percolation, runoff, bsw, 0.0_r8, &
            & dz_layer, fluxes(1:5,indpl), subst_urea, status)
       if (status /= 0) then
          return
       end if
       fluxes(iflx_to_tan, indpl) = rate*pools(indpl)
       age_prev = age_prev + ranges(indpl)
    end do

    ! Here goes also flux_tan!
    call update_pools(pools, fluxes, dt, numpools, 6)

    balance = sum(pools) - old_total
    if (abs(balance - (ndep-sum(fluxes))*dt + garbage) > 1e-9) then
       print *, balance, 'f', sum(fluxes)*dt, ndep*dt, (ndep-sum(fluxes))*dt-garbage - balance, &
            'p', pools, 'g', garbage
       status = err_balance_nitr
       return
    end if

    status = 0
    
  end subroutine update_urea
  
  subroutine get_storage_fluxes_tan_ar(manure_excr, tempr_outside, windspeed, fract_direct, &
     & flux_direct, flux_direct_tan, flux_barn, flux_store, flux_resid, flux_resid_tan, &
     & volat_target_barns, volat_target_stores, volat_coef_barns, volat_coef_stores, tan_fract_excr, nn)

    real(8), intent(in), dimension(nn) :: manure_excr, tempr_outside, windspeed, fract_direct
    real(8), intent(out), dimension(nn) :: flux_barn, flux_store, flux_direct, flux_resid, &
         & flux_direct_tan, flux_resid_tan
    real(8), intent(in) :: volat_target_barns, volat_target_stores, volat_coef_barns, volat_coef_stores, tan_fract_excr
    integer, intent(in) :: nn

    integer :: ii, status
    real(r8) :: fluxes_nitr(4), fluxes_tan(4)
    
    do ii = 1, nn
       call eval_fluxes_storage(manure_excr(ii), 'open', tempr_outside(ii), windspeed(ii), fract_direct(ii), &
            & volat_coef_barns, volat_coef_stores, tan_fract_excr, &
            & fluxes_nitr, fluxes_tan, status)

       if (status /= 0) then
          print *, 'Status = ', status
          return
       end if

       flux_direct(ii) = fluxes_nitr(iflx_appl)
       flux_direct_tan(ii) = fluxes_tan(iflx_appl)
       flux_barn(ii) = fluxes_tan(iflx_air_barns)
       flux_store(ii) = fluxes_tan(iflx_air_stores)
       flux_resid(ii) = fluxes_nitr(iflx_to_store)
       flux_resid_tan(ii) = fluxes_tan(iflx_to_store)
       !print *, '1', fluxes_nitr(iflx_appl), flux_direct(ii)
       
    end do
  end subroutine get_storage_fluxes_tan_ar
  
end module FanMod
