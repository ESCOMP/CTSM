module BiogeophysPreFluxCalcsMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialize variables needed by various biogeophysics flux routines (BareGroundFluxes,
  ! CanopyFluxes, etc.)

  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use abortutils              , only : endrun
  use decompMod               , only : bounds_type
  use PatchType               , only : patch
  use ColumnType              , only : col
  use LandunitType            , only : lun
  use clm_varcon              , only : spval
  use clm_varpar              , only : nlevgrnd, nlevsno, nlevurb, nlevmaxurbgrnd
  use clm_varctl              , only : use_fates, z0param_method, iulog
  use pftconMod               , only : pftcon, noveg
  use column_varcon           , only : icol_roof, icol_sunwall, icol_shadewall
  use landunit_varcon         , only : istsoil, istcrop, istice
  use clm_varcon              , only : hvap, hsub
  use CLMFatesInterfaceMod    , only : hlm_fates_interface_type
  use atm2lndType             , only : atm2lnd_type
  use CanopyStateType         , only : canopystate_type
  use FrictionVelocityMod     , only : frictionvel_type
  use EnergyFluxType          , only : energyflux_type
  use SoilStateType           , only : soilstate_type
  use TemperatureType         , only : temperature_type
  use Wateratm2lndBulkType    , only : wateratm2lndbulk_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use WaterStateBulkType      , only : waterstatebulk_type
  use SurfaceResistanceMod    , only : calc_soilevap_resis
  use WaterFluxBulkType     , only : waterfluxbulk_type

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BiogeophysPreFluxCalcs ! Do various calculations that need to happen before the main biogeophysics flux calculations
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SetZ0mDisp
  private :: CalcInitialTemperatureAndEnergyVars

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BiogeophysPreFluxCalcs(bounds, &
       num_nolakec, filter_nolakec, num_nolakep, filter_nolakep, &
       num_urbanc, filter_urbanc, &
       clm_fates, atm2lnd_inst, canopystate_inst, energyflux_inst, frictionvel_inst, &
       soilstate_inst, temperature_inst, &
       wateratm2lndbulk_inst, waterdiagnosticbulk_inst, waterstatebulk_inst, waterfluxbulk_inst)
    !
    ! !DESCRIPTION:
    ! Do various calculations that need to happen before the main biogeophysics flux calculations
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds    
    integer                        , intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
    integer                        , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    integer                        , intent(in)    :: num_nolakep       ! number of column non-lake points in patch filter
    integer                        , intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
    integer                        , intent(in)    :: num_urbanc        ! number of urban columns in clump
    integer                        , intent(in)    :: filter_urbanc(:)  ! urban column filter
    type(hlm_fates_interface_type) , intent(in)    :: clm_fates
    type(atm2lnd_type)             , intent(in)    :: atm2lnd_inst
    type(canopystate_type)         , intent(inout) :: canopystate_inst
    type(energyflux_type)          , intent(inout) :: energyflux_inst
    type(frictionvel_type)         , intent(inout) :: frictionvel_inst
    type(soilstate_type)           , intent(inout) :: soilstate_inst
    type(temperature_type)         , intent(inout) :: temperature_inst
    type(wateratm2lndbulk_type)    , intent(in)    :: wateratm2lndbulk_inst
    type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
    type(waterstatebulk_type)      , intent(in)    :: waterstatebulk_inst
    type(waterfluxbulk_type)       , intent(in)    :: waterfluxbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p

    character(len=*), parameter :: subname = 'BiogeophysPreFluxCalcs'
    !-----------------------------------------------------------------------

    call SetZ0mDisp(bounds, num_nolakep, filter_nolakep, &
         clm_fates, frictionvel_inst, canopystate_inst)

    call frictionvel_inst%SetRoughnessLengthsAndForcHeightsNonLake(bounds, &
         num_nolakec, filter_nolakec,                       &
         num_nolakep, filter_nolakep,                       &
         atm2lnd_inst, waterdiagnosticbulk_inst, canopystate_inst, &
          waterfluxbulk_inst)

    call CalcInitialTemperatureAndEnergyVars(bounds, &
         num_nolakec, filter_nolakec,                       &
         num_nolakep, filter_nolakep,                       &
         num_urbanc, filter_urbanc,                         &
         atm2lnd_inst, canopystate_inst, frictionvel_inst, &
         wateratm2lndbulk_inst, &
         waterdiagnosticbulk_inst, waterstatebulk_inst, &
         temperature_inst, energyflux_inst)

    ! calculate moisture stress/resistance for soil evaporation
    call calc_soilevap_resis(bounds, &
         num_nolakec, filter_nolakec,                       &
         soilstate_inst, &
         waterstatebulk_inst, waterdiagnosticbulk_inst, &
         temperature_inst)

  end subroutine BiogeophysPreFluxCalcs

  !-----------------------------------------------------------------------
  subroutine SetZ0mDisp(bounds, num_nolakep, filter_nolakep, &
       clm_fates, frictionvel_inst, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Set z0m and displa
    !
    ! !USES:
    use clm_time_manager, only : is_first_step, get_nstep, is_beg_curr_year
    use clm_varcon      , only : cd1_param
    use decompMod       , only : subgrid_level_patch
    use BalanceCheckMod , only : GetBalanceCheckSkipSteps
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds    
    integer                        , intent(in)    :: num_nolakep       ! number of column non-lake points in patch filter
    integer                        , intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
    type(hlm_fates_interface_type) , intent(in)    :: clm_fates
    type(frictionvel_type)         , intent(in)    :: frictionvel_inst
    type(canopystate_type)         , intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p, c

    character(len=*), parameter :: subname = 'SetZ0mDisp'
    real(r8) :: U_ustar                                                 ! wind at canopy height divided by friction velocity (unitless)

    !-----------------------------------------------------------------------

    associate( &
         z0mg             =>    frictionvel_inst%z0mg_col             , & ! Input:  [real(r8) (:)   ] roughness length of ground, momentum [m]
         htop             =>    canopystate_inst%htop_patch           , & ! Input:  [real(r8) (:)   ] canopy top (m)                           
         z0m              =>    canopystate_inst%z0m_patch            , & ! Output: [real(r8) (:)   ] momentum roughness length (m)
         displa           =>    canopystate_inst%displa_patch           & ! Output: [real(r8) (:)   ] displacement height (m)
         )

    ! Set roughness and displacement
    ! Note that FATES passes back z0m and displa at the end
    ! of its dynamics call.  If and when crops are
    ! enabled simultaneously with FATES, we will 
    ! have to apply a filter here.
    if(use_fates) then
       call clm_fates%TransferZ0mDisp(bounds, &
            z0m_patch = z0m(bounds%begp:bounds%endp), &
            displa_patch = displa(bounds%begp:bounds%endp))
    end if

    do fp = 1, num_nolakep
       p = filter_nolakep(fp)
       c = patch%column(p)

       if( .not.(patch%is_fates(p))) then
         select case (z0param_method)
         case ('ZengWang2007')

            z0m(p)    = pftcon%z0mr(patch%itype(p)) * htop(p)
            displa(p) = pftcon%displar(patch%itype(p)) * htop(p)

         case ('Meier2022')

            ! Don't set on first few steps of a simulation, since htop isn't set yet, need to wait until after first do_alb time
            if ( is_first_step() .or. get_nstep() <= GetBalanceCheckSkipSteps()-1 ) then
               z0m(p)    = 0._r8
               displa(p) = 0._r8
               cycle
            ! If a crop type and it's the start of the year, htop gets reset to
            ! zero...
            else if ( is_beg_curr_year() .and. pftcon%crop(patch%itype(p)) /= 0.0_r8 )then
               z0m(p)    = 0._r8
               displa(p) = 0._r8
            end if

            if (patch%itype(p) == noveg) then
               z0m(p)    = 0._r8
               displa(p) = 0._r8

            else
               ! Compute as if elai+esai = LAImax in CanopyFluxes
               displa(p) = htop(p) * (1._r8 - (1._r8 - exp(-(cd1_param * (pftcon%z0v_LAImax(patch%itype(p))))**0.5_r8)) &
                           / (cd1_param*(pftcon%z0v_LAImax(patch%itype(p)) ))**0.5_r8)

               U_ustar = 4._r8 * (pftcon%z0v_Cs(patch%itype(p)) + pftcon%z0v_Cr(patch%itype(p)) *  (pftcon%z0v_LAImax(patch%itype(p))) &
                         / 2._r8)**(-0.5_r8) /  (pftcon%z0v_LAImax(patch%itype(p))) / pftcon%z0v_c(patch%itype(p))

               if ( htop(p) <= 1.e-10_r8 )then
                  z0m(p) = z0mg(c)
               else
                  z0m(p) = htop(p) * (1._r8 - displa(p) / htop(p)) * exp(-0.4_r8 * U_ustar + &
                              log(pftcon%z0v_cw(patch%itype(p))) - 1._r8 + pftcon%z0v_cw(patch%itype(p))**(-1._r8))
               end if

            end if


         end select

       end if
    end do

    end associate

  end subroutine SetZ0mDisp


  !-----------------------------------------------------------------------
  subroutine CalcInitialTemperatureAndEnergyVars(bounds, &
       num_nolakec, filter_nolakec, num_nolakep, filter_nolakep, &
       num_urbanc, filter_urbanc, &
       atm2lnd_inst, canopystate_inst, frictionvel_inst, &
       wateratm2lndbulk_inst, waterdiagnosticbulk_inst, waterstatebulk_inst, &
       temperature_inst, energyflux_inst)
    !
    ! !DESCRIPTION:
    ! Near the start of the time step, calculate initial temperature and energy variables
    ! needed to calculate various biogeophys fluxes, and for the sake of starting with
    ! zeroed-out fluxes everywhere.
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds    
    integer                        , intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
    integer                        , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    integer                        , intent(in)    :: num_nolakep       ! number of column non-lake points in patch filter
    integer                        , intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
    integer                        , intent(in)    :: num_urbanc        ! number of urban columns in clump
    integer                        , intent(in)    :: filter_urbanc(:)  ! urban column filter
    type(atm2lnd_type)             , intent(in)    :: atm2lnd_inst
    type(canopystate_type)         , intent(in)    :: canopystate_inst
    type(frictionvel_type)         , intent(in)    :: frictionvel_inst
    type(wateratm2lndbulk_type)    , intent(in)    :: wateratm2lndbulk_inst
    type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
    type(waterstatebulk_type)      , intent(in)    :: waterstatebulk_inst
    type(temperature_type)         , intent(inout) :: temperature_inst
    type(energyflux_type)          , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: l
    integer  :: fc, c
    integer  :: fp, p
    integer  :: j
    real(r8) :: avmuir       ! ir inverse optical depth per unit leaf area

    character(len=*), parameter :: subname = 'CalcInitialTemperatureAndEnergyVars'
    !-----------------------------------------------------------------------

    associate( &
         snl              =>    col%snl                               , & ! Input:  [integer  (:)   ] number of snow layers
         zii              =>    col%zii                               , & ! Output: [real(r8) (:)   ] convective boundary height [m]
         urbpoi           =>    lun%urbpoi                            , & ! Input:  [logical  (:)   ] true => landunit is an urban point
         forc_t           =>    atm2lnd_inst%forc_t_downscaled_col    , & ! Input:  [real(r8) (:)   ] atmospheric temperature (Kelvin)         
         forc_th          =>    atm2lnd_inst%forc_th_downscaled_col   , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)
         elai             =>    canopystate_inst%elai_patch           , & ! Input:  [real(r8) (:)   ] one-sided leaf area index with burying by snow
         esai             =>    canopystate_inst%esai_patch           , & ! Input:  [real(r8) (:)   ] one-sided stem area index with burying by snow
         forc_hgt_t_patch =>    frictionvel_inst%forc_hgt_t_patch     , & ! Input: [real(r8) (:)   ] observational height of temperature at patch level [m]
         frac_sno_eff     =>    waterdiagnosticbulk_inst%frac_sno_eff_col      , & ! Input:  [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)
         frac_sno         =>    waterdiagnosticbulk_inst%frac_sno_col          , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         frac_h2osfc      =>    waterdiagnosticbulk_inst%frac_h2osfc_col       , & ! Input:  [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Input:  [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)
         forc_q           =>    wateratm2lndbulk_inst%forc_q_downscaled_col    , & ! Input:  [real(r8) (:)   ] atmospheric specific humidity (kg/kg)    
         t_soisno         =>    temperature_inst%t_soisno_col         , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
         t_h2osfc         =>    temperature_inst%t_h2osfc_col         , & ! Input:  [real(r8) (:)   ] surface water temperature
         tssbef           =>    temperature_inst%t_ssbef_col          , & ! Output: [real(r8) (:,:) ] soil/snow temperature before update
         t_h2osfc_bef     =>    temperature_inst%t_h2osfc_bef_col     , & ! Output: [real(r8) (:)   ] saved surface water temperature
         t_grnd           =>    temperature_inst%t_grnd_col           , & ! Output: [real(r8) (:)   ] ground temperature (Kelvin)
         emg              =>    temperature_inst%emg_col              , & ! Output: [real(r8) (:)   ] ground emissivity
         emv              =>    temperature_inst%emv_patch            , & ! Output: [real(r8) (:)   ] vegetation emissivity                    
         beta             =>    temperature_inst%beta_col             , & ! Output: [real(r8) (:)   ] coefficient of convective velocity [-]
         thv              =>    temperature_inst%thv_col              , & ! Output: [real(r8) (:)   ] virtual potential temperature (kelvin)
         thm              =>    temperature_inst%thm_patch            , & ! Output: [real(r8) (:)   ] intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
         htvp             =>    energyflux_inst%htvp_col              , & ! Output: [real(r8) (:)   ] latent heat of vapor of water (or sublimation) [j/kg]
         eflx_sh_tot      =>    energyflux_inst%eflx_sh_tot_patch     , & ! Output: [real(r8) (:)   ] total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_tot_r    =>    energyflux_inst%eflx_sh_tot_r_patch   , & ! Output: [real(r8) (:)   ] rural total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_tot_u    =>    energyflux_inst%eflx_sh_tot_u_patch   , & ! Output: [real(r8) (:)   ] urban total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_veg      =>    energyflux_inst%eflx_sh_veg_patch     , & ! Output: [real(r8) (:)   ] sensible heat flux from leaves (W/m**2) [+ to atm]
         eflx_lh_tot      =>    energyflux_inst%eflx_lh_tot_patch     , & ! Output: [real(r8) (:)   ] total latent heat flux (W/m**2)  [+ to atm]
         eflx_lh_tot_r    =>    energyflux_inst%eflx_lh_tot_r_patch   , & ! Output: [real(r8) (:)   ] rural total latent heat flux (W/m**2)  [+ to atm]
         eflx_lh_tot_u    =>    energyflux_inst%eflx_lh_tot_u_patch   , & ! Output: [real(r8) (:)   ] urban total latent heat flux (W/m**2)  [+ to atm]
         cgrnd            =>    energyflux_inst%cgrnd_patch           , & ! Output: [real(r8) (:)   ] deriv. of soil energy flux wrt to soil temp [w/m2/k]
         cgrnds           =>    energyflux_inst%cgrnds_patch          , & ! Output: [real(r8) (:)   ] deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
         cgrndl           =>    energyflux_inst%cgrndl_patch            & ! Output: [real(r8) (:)   ] deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
         )

    do j = -nlevsno+1, nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
               .and. col%itype(c) /= icol_roof) then
             tssbef(c,j) = t_soisno(c,j)
          end if
       end do
    end do

    do j = -nlevsno+1, nlevmaxurbgrnd
       do fc = 1,num_urbanc
          c = filter_urbanc(fc)
          if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
               .or. col%itype(c) == icol_roof) then
             if (j > nlevurb) then
                tssbef(c,j) = spval 
             else
                tssbef(c,j) = t_soisno(c,j)
             end if
          end if
       end do
    end do

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       l = col%landunit(c)

       ! record t_h2osfc prior to updating
       t_h2osfc_bef(c) = t_h2osfc(c)

       ! ground temperature is weighted average of exposed soil, snow, and h2osfc
       if (snl(c) < 0) then
          t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
               + (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) * t_soisno(c,1) &
               + frac_h2osfc(c) * t_h2osfc(c)
       else
          t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
       end if

       ! Ground emissivity - only calculate for non-urban landunits 
       ! Urban emissivities are currently read in from data file
       if (.not. urbpoi(l)) then
          if (lun%itype(l)==istice) then
             emg(c) = 0.97_r8
          else
             emg(c) = (1._r8-frac_sno(c))*0.96_r8 + frac_sno(c)*0.97_r8
          end if
       end if

       ! Latent heat. We arbitrarily assume that the sublimation occurs
       ! only as h2osoi_liq = 0
       if (h2osoi_liq(c,snl(c)+1) <= 0._r8 .and. h2osoi_ice(c,snl(c)+1) > 0._r8) then
          htvp(c) = hsub
       else
          htvp(c) = hvap
       end if

       ! Potential, virtual potential temperature, and wind speed at the
       ! reference height
       beta(c) = 1._r8
       zii(c)  = 1000._r8
       thv(c)  = forc_th(c)*(1._r8+0.61_r8*forc_q(c))
    end do

    do fp = 1, num_nolakep
       p = filter_nolakep(fp)
       c = patch%column(p)
       l = patch%landunit(p)

       ! Initial set (needed for history tape fields)

       eflx_sh_tot(p) = 0._r8
       if (urbpoi(l)) then
          eflx_sh_tot_u(p) = 0._r8
       else if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then 
          eflx_sh_tot_r(p) = 0._r8
       end if
       eflx_lh_tot(p) = 0._r8
       if (urbpoi(l)) then
          eflx_lh_tot_u(p) = 0._r8
       else if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then 
          eflx_lh_tot_r(p) = 0._r8
       end if
       eflx_sh_veg(p) = 0._r8

       ! Initial set for calculation

       cgrnd(p)  = 0._r8
       cgrnds(p) = 0._r8
       cgrndl(p) = 0._r8

       ! Vegetation Emissivity

       avmuir = 1._r8
       emv(p) = 1._r8-exp(-(elai(p)+esai(p))/avmuir)

       ! thm
       thm(p)  = forc_t(c) + 0.0098_r8*forc_hgt_t_patch(p)
    end do

    end associate

  end subroutine CalcInitialTemperatureAndEnergyVars

end module BiogeophysPreFluxCalcsMod
