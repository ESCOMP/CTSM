module lnd2atmMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle lnd2atm mapping
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use shr_megan_mod        , only : shr_megan_mechcomps_n
  use shr_fire_emis_mod    , only : shr_fire_emis_mechcomps_n
  use clm_varpar           , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon           , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl           , only : iulog, use_lch4
  use seq_drydep_mod       , only : n_drydep, drydep_method, DD_XLND
  use decompMod            , only : bounds_type
  use subgridAveMod        , only : p2g, c2g 
  use lnd2atmType          , only : lnd2atm_type
  use atm2lndType          , only : atm2lnd_type
  use ch4Mod               , only : ch4_type
  use DUSTMod              , only : dust_type
  use DryDepVelocity       , only : drydepvel_type
  use VocEmissionMod       , only : vocemis_type
  use CNFireEmissionsMod   , only : fireemis_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityMod  , only : frictionvel_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use TemperatureType      , only : temperature_type
  use WaterFluxBulkType        , only : waterfluxbulk_type
  use WaterStateBulkType       , only : waterstatebulk_type
  use WaterDiagnosticBulkType       , only : waterdiagnosticbulk_type
  use WaterBalanceType       , only : waterbalance_type
  use glcBehaviorMod       , only : glc_behavior_type
  use glc2lndMod           , only : glc2lnd_type
  use ColumnType           , only : col
  use LandunitType         , only : lun
  use GridcellType         , only : grc                
  use landunit_varcon      , only : istice_mec
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm
  public :: lnd2atm_minimal

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: handle_ice_runoff

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine lnd2atm_minimal(bounds, &
      waterstatebulk_inst, surfalb_inst, energyflux_inst, lnd2atm_inst)
    !
    ! !DESCRIPTION:
    ! Compute clm_l2a_inst component of gridcell derived type. This routine computes
    ! the bare minimum of components necessary to get the first step of a
    ! run started.
    !
    ! !USES:
    use clm_varcon, only : sb
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds  
    type(waterstatebulk_type) , intent(in)    :: waterstatebulk_inst
    type(surfalb_type)    , intent(in)    :: surfalb_inst
    type(energyflux_type) , intent(in)    :: energyflux_inst
    type(lnd2atm_type)    , intent(inout) :: lnd2atm_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: g                                    ! index
    real(r8), parameter :: amC   = 12.0_r8          ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8          ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    call c2g(bounds, &
         waterstatebulk_inst%h2osno_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%h2osno_grc    (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_inst%h2osno_grc(g) = lnd2atm_inst%h2osno_grc(g)/1000._r8
    end do

    call c2g(bounds, nlevgrnd, &
         waterstatebulk_inst%h2osoi_vol_col (bounds%begc:bounds%endc, :), &
         lnd2atm_inst%h2osoi_vol_grc    (bounds%begg:bounds%endg, :), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_inst%albd_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_inst%albd_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_inst%albi_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_inst%albi_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_inst%eflx_lwrad_out_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%eflx_lwrad_out_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_inst%t_rad_grc(g) = sqrt(sqrt(lnd2atm_inst%eflx_lwrad_out_grc(g)/sb))
    end do

  end subroutine lnd2atm_minimal

  !------------------------------------------------------------------------
  subroutine lnd2atm(bounds, &
       atm2lnd_inst, surfalb_inst, temperature_inst, frictionvel_inst, &
       waterstatebulk_inst, waterdiagnosticbulk_inst, waterbalancebulk_inst, waterfluxbulk_inst, energyflux_inst, &
       solarabs_inst, drydepvel_inst,  &
       vocemis_inst, fireemis_inst, dust_inst, ch4_inst, glc_behavior, &
       lnd2atm_inst, &
       net_carbon_exchange_grc) 
    !
    ! !DESCRIPTION:
    ! Compute lnd2atm_inst component of gridcell derived type
    !
    ! !USES:
    use ch4varcon  , only : ch4offline
    !
    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds  
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_inst
    type(surfalb_type)          , intent(in)    :: surfalb_inst
    type(temperature_type)      , intent(in)    :: temperature_inst
    type(frictionvel_type)      , intent(in)    :: frictionvel_inst
    type(waterstatebulk_type)       , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)       , intent(inout) :: waterdiagnosticbulk_inst
    type(waterbalance_type)       , intent(inout) :: waterbalancebulk_inst
    type(waterfluxbulk_type)        , intent(inout) :: waterfluxbulk_inst
    type(energyflux_type)       , intent(in)    :: energyflux_inst
    type(solarabs_type)         , intent(in)    :: solarabs_inst
    type(drydepvel_type)        , intent(in)    :: drydepvel_inst
    type(vocemis_type)          , intent(in)    :: vocemis_inst
    type(fireemis_type)         , intent(in)    :: fireemis_inst
    type(dust_type)             , intent(in)    :: dust_inst
    type(ch4_type)              , intent(in)    :: ch4_inst
    type(glc_behavior_type)     , intent(in)    :: glc_behavior
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_inst 
    real(r8)                    , intent(in)    :: net_carbon_exchange_grc( bounds%begg: )  ! net carbon exchange between land and atmosphere, positive for source (gC/m2/s)
    !
    ! !LOCAL VARIABLES:
    integer  :: c, g  ! indices
    real(r8) :: qflx_ice_runoff_col(bounds%begc:bounds%endc) ! total column-level ice runoff
    real(r8) :: eflx_sh_ice_to_liq_grc(bounds%begg:bounds%endg) ! sensible heat flux generated from the ice to liquid conversion, averaged to gridcell
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(net_carbon_exchange_grc) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))

    call handle_ice_runoff(bounds, waterfluxbulk_inst, glc_behavior, &
         melt_non_icesheet_ice_runoff = lnd2atm_inst%params%melt_non_icesheet_ice_runoff, &
         qflx_ice_runoff_col = qflx_ice_runoff_col(bounds%begc:bounds%endc), &
         qflx_liq_from_ice_col = lnd2atm_inst%qflx_liq_from_ice_col(bounds%begc:bounds%endc), &
         eflx_sh_ice_to_liq_col = lnd2atm_inst%eflx_sh_ice_to_liq_col(bounds%begc:bounds%endc))

    !----------------------------------------------------
    ! lnd -> atm
    !----------------------------------------------------
    
    ! First, compute the "minimal" set of fields.
    call lnd2atm_minimal(bounds, &
         waterstatebulk_inst, surfalb_inst, energyflux_inst, lnd2atm_inst)

    call p2g(bounds, &
         temperature_inst%t_ref2m_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%t_ref2m_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         waterdiagnosticbulk_inst%q_ref2m_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%q_ref2m_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_inst%u10_clm_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%u_ref10m_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_inst%taux_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%taux_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_inst%tauy_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%tauy_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         waterfluxbulk_inst%qflx_evap_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%qflx_evap_tot_grc     (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         solarabs_inst%fsa_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%fsa_grc    (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_inst%fv_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%fv_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_inst%ram1_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%ram1_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g( bounds, &
         energyflux_inst%eflx_sh_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%eflx_sh_tot_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity',c2l_scale_type='urbanf',l2g_scale_type='unity')
    call c2g( bounds, &
         energyflux_inst%eflx_sh_precip_conversion_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%eflx_sh_precip_conversion_grc    (bounds%begg:bounds%endg), &
         c2l_scale_type='urbanf', l2g_scale_type='unity')
    call c2g( bounds, &
         lnd2atm_inst%eflx_sh_ice_to_liq_col(bounds%begc:bounds%endc), &
         eflx_sh_ice_to_liq_grc(bounds%begg:bounds%endg), &
         c2l_scale_type='urbanf', l2g_scale_type='unity')
    do g = bounds%begg, bounds%endg
       lnd2atm_inst%eflx_sh_tot_grc(g) =  lnd2atm_inst%eflx_sh_tot_grc(g) + &
            lnd2atm_inst%eflx_sh_precip_conversion_grc(g) + &
            eflx_sh_ice_to_liq_grc(g) - &
            energyflux_inst%eflx_dynbal_grc(g) 
    enddo

    call p2g(bounds, &
         energyflux_inst%eflx_lh_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%eflx_lh_tot_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg, bounds%endg
       lnd2atm_inst%net_carbon_exchange_grc(g) = &
            net_carbon_exchange_grc(g)
    end do
    if (use_lch4) then
       if (.not. ch4offline) then
          ! Adjust flux of CO2 by the net conversion of mineralizing C to CH4
          do g = bounds%begg,bounds%endg
             ! nem is in g C/m2/sec
             lnd2atm_inst%net_carbon_exchange_grc(g) = &
                  lnd2atm_inst%net_carbon_exchange_grc(g) + lnd2atm_inst%nem_grc(g)
          end do
       end if
    end if
    ! Convert from gC/m2/s to kgCO2/m2/s
    do g = bounds%begg,bounds%endg
       lnd2atm_inst%net_carbon_exchange_grc(g) = &
            lnd2atm_inst%net_carbon_exchange_grc(g)*convertgC2kgCO2
    end do

    ! drydepvel
    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
       call p2g(bounds, n_drydep, &
            drydepvel_inst%velocity_patch (bounds%begp:bounds%endp, :), &
            lnd2atm_inst%ddvel_grc        (bounds%begg:bounds%endg, :), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    endif

    ! voc emission flux
    if (shr_megan_mechcomps_n>0) then
       call p2g(bounds, shr_megan_mechcomps_n, &
            vocemis_inst%vocflx_patch(bounds%begp:bounds%endp,:), &
            lnd2atm_inst%flxvoc_grc  (bounds%begg:bounds%endg,:), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    end if

    ! fire emissions fluxes
     if (shr_fire_emis_mechcomps_n>0) then
        call p2g(bounds, shr_fire_emis_mechcomps_n, &
            -fireemis_inst%fireflx_patch(bounds%begp:bounds%endp,:), &
             lnd2atm_inst%fireflx_grc   (bounds%begg:bounds%endg,:), &
             p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
        call p2g(bounds, &
             fireemis_inst%ztop_patch (bounds%begp:bounds%endp), &
             lnd2atm_inst%fireztop_grc(bounds%begg:bounds%endg), &
             p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
     endif

    ! dust emission flux
    call p2g(bounds, ndst, &
         dust_inst%flx_mss_vrt_dst_patch(bounds%begp:bounds%endp, :), &
         lnd2atm_inst%flxdst_grc        (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')


    ! ch4 flux
    if (use_lch4) then
       call c2g( bounds,     &
            ch4_inst%ch4_surf_flux_tot_col (bounds%begc:bounds%endc), &
            lnd2atm_inst%flux_ch4_grc      (bounds%begg:bounds%endg), &
            c2l_scale_type= 'unity', l2g_scale_type='unity' )
    end if

    !----------------------------------------------------
    ! lnd -> rof
    !----------------------------------------------------

    call c2g( bounds, &
         waterfluxbulk_inst%qflx_surf_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%qflx_rofliq_qsur_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         waterfluxbulk_inst%qflx_drain_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%qflx_rofliq_qsub_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    do c = bounds%begc, bounds%endc
       if (col%active(c)) then
          ! It's not entirely appropriate to put qflx_liq_from_ice_col into
          ! qflx_qrgwl_col, since this isn't necessarily just glaciers, wetlands and
          ! lakes. But since we put the liquid portion of snow capping into
          ! qflx_qrgwl_col, it seems reasonable to put qflx_liq_from_ice_col there as
          ! well.
          waterfluxbulk_inst%qflx_qrgwl_col(c) = waterfluxbulk_inst%qflx_qrgwl_col(c) + &
               lnd2atm_inst%qflx_liq_from_ice_col(c)

          ! qflx_runoff is the sum of a number of terms, including qflx_qrgwl. Since we
          ! are adjusting qflx_qrgwl above, we need to adjust qflx_runoff analogously.
          waterfluxbulk_inst%qflx_runoff_col(c) = waterfluxbulk_inst%qflx_runoff_col(c) + &
               lnd2atm_inst%qflx_liq_from_ice_col(c)
       end if
    end do

    call c2g( bounds, &
         waterfluxbulk_inst%qflx_qrgwl_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%qflx_rofliq_qgwl_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         waterfluxbulk_inst%qflx_runoff_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%qflx_rofliq_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    do g = bounds%begg, bounds%endg
       lnd2atm_inst%qflx_rofliq_qgwl_grc(g) = lnd2atm_inst%qflx_rofliq_qgwl_grc(g) - waterfluxbulk_inst%qflx_liq_dynbal_grc(g)
       lnd2atm_inst%qflx_rofliq_grc(g) = lnd2atm_inst%qflx_rofliq_grc(g) - waterfluxbulk_inst%qflx_liq_dynbal_grc(g)
    enddo

    call c2g( bounds, &
         waterfluxbulk_inst%qflx_drain_perched_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%qflx_rofliq_drain_perched_grc(bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         waterfluxbulk_inst%qflx_irrig_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%qirrig_grc(bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         qflx_ice_runoff_col(bounds%begc:bounds%endc),  &
         lnd2atm_inst%qflx_rofice_grc(bounds%begg:bounds%endg),  & 
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       lnd2atm_inst%qflx_rofice_grc(g) = lnd2atm_inst%qflx_rofice_grc(g) - waterfluxbulk_inst%qflx_ice_dynbal_grc(g)          
    enddo

    ! calculate total water storage for history files
    ! first set tws to gridcell total endwb
    ! second add river storage as gridcell average depth (1.e-3 converts [m3/km2] to [mm])
    ! TODO - this was in BalanceCheckMod - not sure where it belongs?

    call c2g( bounds, &
         waterbalancebulk_inst%endwb_col(bounds%begc:bounds%endc), &
         waterdiagnosticbulk_inst%tws_grc  (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       waterdiagnosticbulk_inst%tws_grc(g) = waterdiagnosticbulk_inst%tws_grc(g) + atm2lnd_inst%volr_grc(g) / grc%area(g) * 1.e-3_r8
    enddo

  end subroutine lnd2atm

  !-----------------------------------------------------------------------
  subroutine handle_ice_runoff(bounds, waterfluxbulk_inst, glc_behavior, &
       melt_non_icesheet_ice_runoff, &
       qflx_ice_runoff_col, qflx_liq_from_ice_col, eflx_sh_ice_to_liq_col)
    !
    ! !DESCRIPTION:
    ! Take column-level ice runoff and divide it between (a) ice runoff, and (b) liquid
    ! runoff with a compensating negative sensible heat flux.
    !
    ! The rationale here is: Ice runoff is largely meant to represent a crude
    ! parameterization of iceberg calving. Iceberg calving is mainly appropriate in
    ! regions where an ice sheet terminates at the land-ocean boundary. Elsewhere, in
    ! reality, we expect most ice runoff to flow downstream and melt before it reaches the
    ! ocean. Furthermore, sending ice runoff directly to the ocean can lead to runaway sea
    ! ice growth in some regions (around the Canadian archipelago, and possibly in more
    ! wide-spread regions of the Arctic Ocean); melting this ice before it reaches the
    ! ocean avoids this problem.
    !
    ! If the river model were able to melt ice, then we might not need this routine.
    !
    ! Note that this routine does NOT handle ice runoff generated via the dynamic
    ! landunits adjustment fluxes (i.e., the fluxes that compensate for a difference in
    ! ice content between the pre- and post-dynamic landunit areas). This is partly
    ! because those gridcell-level dynamic landunits adjustment fluxes do not fit well
    ! with this column-based infrastructure, and partly because either method of handling
    ! these fluxes (i.e., sending an ice runoff or sending a liquid runoff with a
    ! negative sensible heat flux) seems equally justifiable.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(waterfluxbulk_type), intent(in) :: waterfluxbulk_inst
    type(glc_behavior_type), intent(in) :: glc_behavior
    logical, intent(in) :: melt_non_icesheet_ice_runoff
    real(r8), intent(out) :: qflx_ice_runoff_col( bounds%begc: ) ! total column-level ice runoff (mm H2O /s)
    real(r8), intent(out) :: qflx_liq_from_ice_col( bounds%begc: ) ! liquid runoff from converted ice runoff (mm H2O /s)
    real(r8), intent(out) :: eflx_sh_ice_to_liq_col( bounds%begc: ) ! sensible heat flux generated from the ice to liquid conversion (W/m2) (+ to atm)

    !
    ! !LOCAL VARIABLES:
    integer :: c, l, g
    logical :: do_conversion

    character(len=*), parameter :: subname = 'handle_ice_runoff'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(qflx_ice_runoff_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_liq_from_ice_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(eflx_sh_ice_to_liq_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    do c = bounds%begc, bounds%endc
       if (col%active(c)) then
          qflx_ice_runoff_col(c) = waterfluxbulk_inst%qflx_ice_runoff_snwcp_col(c) + &
               waterfluxbulk_inst%qflx_ice_runoff_xs_col(c)
          qflx_liq_from_ice_col(c) = 0._r8
          eflx_sh_ice_to_liq_col(c) = 0._r8
       end if
    end do

    if (melt_non_icesheet_ice_runoff) then
       do c = bounds%begc, bounds%endc
          if (col%active(c)) then
             l = col%landunit(c)
             g = col%gridcell(c)
             do_conversion = .false.
             if (lun%itype(l) /= istice_mec) then
                do_conversion = .true.
             else  ! istice_mec
                if (glc_behavior%ice_runoff_melted_grc(g)) then
                   do_conversion = .true.
                else
                   do_conversion = .false.
                end if
             end if
             if (do_conversion) then
                ! ice to liquid absorbs energy, so results in a negative heat flux to atm
                ! Note that qflx_ice_runoff_col is in mm H2O/s, which is the same as kg
                ! m-2 s-1, so we can simply multiply by hfus.
                eflx_sh_ice_to_liq_col(c) = -qflx_ice_runoff_col(c) * hfus
                qflx_liq_from_ice_col(c) = qflx_ice_runoff_col(c)
                qflx_ice_runoff_col(c) = 0._r8
             end if
          end if
       end do
    end if

  end subroutine handle_ice_runoff


end module lnd2atmMod
