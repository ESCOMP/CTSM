module lnd2atmMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle lnd2atm mapping
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_megan_mod        , only : shr_megan_mechcomps_n
  use shr_fire_emis_mod    , only : shr_fire_emis_mechcomps_n
  use clm_varpar           , only : numrad, ndst, nlevgrnd, nlevmaxurbgrnd !ndst = number of dust bins.
  use clm_varcon           , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl           , only : iulog, use_lch4
  use shr_drydep_mod       , only : n_drydep
  use decompMod            , only : bounds_type
  use subgridAveMod        , only : p2g, c2g, l2g
  use filterColMod         , only : filter_col_type, col_filter_from_logical_array
  use lnd2atmType          , only : lnd2atm_type
  use atm2lndType          , only : atm2lnd_type
  use ch4Mod               , only : ch4_type
  use DustEmisBase         , only : dust_emis_base_type
  use DryDepVelocity       , only : drydepvel_type
  use VocEmissionMod       , only : vocemis_type
  use CNFireEmissionsMod   , only : fireemis_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityMod  , only : frictionvel_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use TemperatureType      , only : temperature_type
  use WaterFluxBulkType    , only : waterfluxbulk_type
  use WaterType            , only : water_type
  use glcBehaviorMod       , only : glc_behavior_type
  use glc2lndMod           , only : glc2lnd_type
  use ColumnType           , only : col
  use LandunitType         , only : lun
  use GridcellType         , only : grc
  use landunit_varcon      , only : istice
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
      water_inst, surfalb_inst, energyflux_inst, lnd2atm_inst)
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
    type(water_type)      , intent(inout) :: water_inst
    type(surfalb_type)    , intent(in)    :: surfalb_inst
    type(energyflux_type) , intent(in)    :: energyflux_inst
    type(lnd2atm_type)    , intent(inout) :: lnd2atm_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: i,g                                   ! index
    type(filter_col_type) :: filter_active_c          ! filter for active columns
    real(r8) :: h2osno_total(bounds%begc:bounds%endc) ! total snow water (mm H2O)
    real(r8), parameter :: amC   = 12.0_r8          ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8          ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    ! TODO(wjs, 2019-06-12) We could use the standard allc filter for this. However,
    ! currently lnd2atm_minimal is called from outside a clump loop, both in
    ! initialization and in the driver loop (via the call to lnd2atm), and filters are
    ! only available inside a clump loop. At a glance, it looks like these calls could be
    ! put inside clump loops. But I want to look a little more carefully and/or do some
    ! additional testing (e.g., writing coupler AVGHIST files for one test, and maybe
    ! examining fields that are sent in initialization - and comparing these with
    ! baselines) to help ensure I'm not breaking anything, before making this change. So
    ! for now we build the necessary filters on the fly.
    filter_active_c = col_filter_from_logical_array(bounds, &
         col%active(bounds%begc:bounds%endc))

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(bulk_or_tracer => water_inst%bulk_and_tracers(i))

       call bulk_or_tracer%waterstate_inst%CalculateTotalH2osno( &
            bounds, &
            filter_active_c%num, &
            filter_active_c%indices, &
            caller = 'lnd2atm_minimal: '//water_inst%GetBulkOrTracerName(i), &
            h2osno_total = h2osno_total(bounds%begc:bounds%endc))
       call c2g(bounds, &
            h2osno_total(bounds%begc:bounds%endc), &
            bulk_or_tracer%waterlnd2atm_inst%h2osno_grc(bounds%begg:bounds%endg), &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity')

       do g = bounds%begg,bounds%endg
          bulk_or_tracer%waterlnd2atm_inst%h2osno_grc(g) = &
               bulk_or_tracer%waterlnd2atm_inst%h2osno_grc(g)/1000._r8
       end do
       end associate
    end do

    call c2g(bounds, nlevmaxurbgrnd, &
         water_inst%waterstatebulk_inst%h2osoi_vol_col (bounds%begc:bounds%endc, :), &
         water_inst%waterlnd2atmbulk_inst%h2osoi_vol_grc    (bounds%begg:bounds%endg, :), &
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
       water_inst, &
       energyflux_inst, solarabs_inst, drydepvel_inst,  &
       vocemis_inst, fireemis_inst, dust_emis_inst, ch4_inst, glc_behavior, &
       lnd2atm_inst, &
       net_carbon_exchange_grc)
    !
    ! !DESCRIPTION:
    ! Compute lnd2atm_inst component of gridcell derived type
    !
    ! !USES:
    use ch4varcon  , only : ch4offline
    use clm_varctl , only : use_hillslope_routing
    !
    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_inst
    type(surfalb_type)          , intent(in)    :: surfalb_inst
    type(temperature_type)      , intent(in)    :: temperature_inst
    type(frictionvel_type)      , intent(in)    :: frictionvel_inst
    type(water_type)            , intent(inout) :: water_inst
    type(energyflux_type)       , intent(in)    :: energyflux_inst
    type(solarabs_type)         , intent(in)    :: solarabs_inst
    type(drydepvel_type)        , intent(in)    :: drydepvel_inst
    type(vocemis_type)          , intent(in)    :: vocemis_inst
    type(fireemis_type)         , intent(in)    :: fireemis_inst
    class(dust_emis_base_type)  , intent(in)    :: dust_emis_inst
    type(ch4_type)              , intent(in)    :: ch4_inst
    type(glc_behavior_type)     , intent(in)    :: glc_behavior
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_inst
    real(r8)                    , intent(in)    :: net_carbon_exchange_grc( bounds%begg: )  ! net carbon exchange between land and atmosphere, positive for source (gC/m2/s)
    !
    ! !LOCAL VARIABLES:
    integer  :: c, l, g  ! indices
    real(r8) :: eflx_sh_ice_to_liq_grc(bounds%begg:bounds%endg) ! sensible heat flux generated from the ice to liquid conversion, averaged to gridcell
    real(r8), allocatable :: qflx_surf_col_to_rof(:)          ! surface runoff that is sent directly to rof
    real(r8), allocatable :: qflx_drain_col_to_rof(:)         ! drainagec that is sent directly to rof
    real(r8), allocatable :: qflx_drain_perched_col_to_rof(:) ! perched drainage that is sent directly to rof
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(net_carbon_exchange_grc) == (/bounds%endg/)), sourcefile, __LINE__)

    call handle_ice_runoff(bounds, water_inst%waterfluxbulk_inst, glc_behavior, &
         melt_non_icesheet_ice_runoff = lnd2atm_inst%params%melt_non_icesheet_ice_runoff, &
         qflx_ice_runoff_col = water_inst%waterlnd2atmbulk_inst%qflx_ice_runoff_col(bounds%begc:bounds%endc), &
         qflx_liq_from_ice_col = water_inst%waterlnd2atmbulk_inst%qflx_liq_from_ice_col(bounds%begc:bounds%endc), &
         eflx_sh_ice_to_liq_col = lnd2atm_inst%eflx_sh_ice_to_liq_col(bounds%begc:bounds%endc))

    !----------------------------------------------------
    ! lnd -> atm
    !----------------------------------------------------

    ! First, compute the "minimal" set of fields.
    call lnd2atm_minimal(bounds, &
         water_inst, surfalb_inst, energyflux_inst, lnd2atm_inst)

    call p2g(bounds, &
         temperature_inst%t_ref2m_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%t_ref2m_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         water_inst%waterdiagnosticbulk_inst%q_ref2m_patch (bounds%begp:bounds%endp), &
         water_inst%waterlnd2atmbulk_inst%q_ref2m_grc      (bounds%begg:bounds%endg), &
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
         water_inst%waterfluxbulk_inst%qflx_evap_tot_patch (bounds%begp:bounds%endp), &
         water_inst%waterlnd2atmbulk_inst%qflx_evap_tot_grc     (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         solarabs_inst%fsa_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%fsa_grc    (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_inst%z0m_actual_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%z0m_grc              (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbans', l2g_scale_type='unity')

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
    if ( n_drydep > 0 ) then
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
         dust_emis_inst%flx_mss_vrt_dst_patch(bounds%begp:bounds%endp, :), &
         lnd2atm_inst%flxdst_grc        (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    !----------------------------------------------------
    ! lnd -> rof
    !----------------------------------------------------

    if (use_hillslope_routing) then
       ! streamflow is volume/time, so sum over landunits (do not weight)
       water_inst%waterlnd2atmbulk_inst%qflx_rofliq_stream_grc(bounds%begg:bounds%endg) = 0._r8
       do l = bounds%begl, bounds%endl
          if(lun%active(l)) then
             g = lun%gridcell(l)
             water_inst%waterlnd2atmbulk_inst%qflx_rofliq_stream_grc(g) = &
                  water_inst%waterlnd2atmbulk_inst%qflx_rofliq_stream_grc(g) &
                  +  water_inst%waterfluxbulk_inst%volumetric_streamflow_lun(l) &
                  *1e3_r8/(grc%area(g)*1.e6_r8)
          endif
       enddo

       ! If hillslope routing is used, exclude inputs to stream channel from gridcell averages to avoid double counting
       allocate( &
            qflx_surf_col_to_rof(bounds%begc:bounds%endc), &
            qflx_drain_col_to_rof(bounds%begc:bounds%endc), &
            qflx_drain_perched_col_to_rof(bounds%begc:bounds%endc))

       qflx_surf_col_to_rof(bounds%begc:bounds%endc)  = 0._r8
       qflx_drain_col_to_rof(bounds%begc:bounds%endc) = 0._r8
       qflx_drain_perched_col_to_rof(bounds%begc:bounds%endc) = 0._r8
       
       do c = bounds%begc, bounds%endc
          ! Exclude hillslope columns from gridcell average
          ! hillslope runoff is sent to stream rather than directly
          ! to rof, and is accounted for in qflx_rofliq_stream_grc
          if (col%active(c) .and. .not. col%is_hillslope_column(c)) then
             qflx_surf_col_to_rof(c) = qflx_surf_col_to_rof(c) &
                  + water_inst%waterfluxbulk_inst%qflx_surf_col(c)
             qflx_drain_col_to_rof(c) = qflx_drain_col_to_rof(c) &
                  + water_inst%waterfluxbulk_inst%qflx_drain_col(c)
             qflx_drain_perched_col_to_rof(c) = &
                  qflx_drain_perched_col_to_rof(c) &
                  + water_inst%waterfluxbulk_inst%qflx_drain_perched_col(c)
          endif
       enddo
          
       call c2g( bounds, &
            qflx_surf_col_to_rof  (bounds%begc:bounds%endc), &
            water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc   (bounds%begg:bounds%endg), &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity')
       
       call c2g( bounds, &
            qflx_drain_col_to_rof (bounds%begc:bounds%endc), &
            water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsub_grc   (bounds%begg:bounds%endg), &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity')
       
       call c2g( bounds, &
            qflx_drain_perched_col_to_rof (bounds%begc:bounds%endc), &
            water_inst%waterlnd2atmbulk_inst%qflx_rofliq_drain_perched_grc(bounds%begg:bounds%endg), &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity')
       
       deallocate(qflx_surf_col_to_rof,qflx_drain_col_to_rof, &
            qflx_drain_perched_col_to_rof)

    else
    
       call c2g( bounds, &
            water_inst%waterfluxbulk_inst%qflx_surf_col (bounds%begc:bounds%endc), &
            water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc   (bounds%begg:bounds%endg), &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
       
       call c2g( bounds, &
            water_inst%waterfluxbulk_inst%qflx_drain_col (bounds%begc:bounds%endc), &
            water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qsub_grc   (bounds%begg:bounds%endg), &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
       
       call c2g( bounds, &
            water_inst%waterfluxbulk_inst%qflx_drain_perched_col (bounds%begc:bounds%endc), &
            water_inst%waterlnd2atmbulk_inst%qflx_rofliq_drain_perched_grc(bounds%begg:bounds%endg), &
            c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    endif

    do c = bounds%begc, bounds%endc
       if (col%active(c)) then
          ! It's not entirely appropriate to put qflx_liq_from_ice_col into
          ! qflx_qrgwl_col, since this isn't necessarily just glaciers, wetlands and
          ! lakes. But since we put the liquid portion of snow capping into
          ! qflx_qrgwl_col, it seems reasonable to put qflx_liq_from_ice_col there as
          ! well.
          water_inst%waterfluxbulk_inst%qflx_qrgwl_col(c) = water_inst%waterfluxbulk_inst%qflx_qrgwl_col(c) + &
               water_inst%waterlnd2atmbulk_inst%qflx_liq_from_ice_col(c)

          ! qflx_runoff is the sum of a number of terms, including qflx_qrgwl. Since we
          ! are adjusting qflx_qrgwl above, we need to adjust qflx_runoff analogously.
          water_inst%waterfluxbulk_inst%qflx_runoff_col(c) = &
            water_inst%waterfluxbulk_inst%qflx_runoff_col(c) + &
            water_inst%waterlnd2atmbulk_inst%qflx_liq_from_ice_col(c)
       end if
    end do

    call c2g( bounds, &
         water_inst%waterfluxbulk_inst%qflx_qrgwl_col (bounds%begc:bounds%endc), &
         water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qgwl_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         water_inst%waterfluxbulk_inst%qflx_runoff_col (bounds%begc:bounds%endc), &
         water_inst%waterlnd2atmbulk_inst%qflx_rofliq_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    do g = bounds%begg, bounds%endg
       water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qgwl_grc(g) = &
         water_inst%waterlnd2atmbulk_inst%qflx_rofliq_qgwl_grc(g) - &
         water_inst%waterfluxbulk_inst%qflx_liq_dynbal_grc(g)
       water_inst%waterlnd2atmbulk_inst%qflx_rofliq_grc(g) = &
         water_inst%waterlnd2atmbulk_inst%qflx_rofliq_grc(g) - &
         water_inst%waterfluxbulk_inst%qflx_liq_dynbal_grc(g)
    enddo

    call c2g( bounds, &
         water_inst%waterfluxbulk_inst%qflx_sfc_irrig_col (bounds%begc:bounds%endc), &
         water_inst%waterlnd2atmbulk_inst%qirrig_grc(bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         water_inst%waterlnd2atmbulk_inst%qflx_ice_runoff_col(bounds%begc:bounds%endc),  &
         water_inst%waterlnd2atmbulk_inst%qflx_rofice_grc(bounds%begg:bounds%endg),  &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       water_inst%waterlnd2atmbulk_inst%qflx_rofice_grc(g) = &
         water_inst%waterlnd2atmbulk_inst%qflx_rofice_grc(g) - &
         water_inst%waterfluxbulk_inst%qflx_ice_dynbal_grc(g)
    enddo

    ! calculate total water storage for history files
    ! first set tws to gridcell total endwb
    ! second add river storage as gridcell average depth (1.e-3 converts [m3/km2] to [mm])
    ! TODO - this was in BalanceCheckMod - not sure where it belongs?

    call c2g( bounds, &
         water_inst%waterbalancebulk_inst%endwb_col(bounds%begc:bounds%endc), &
         water_inst%waterdiagnosticbulk_inst%tws_grc(bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       water_inst%waterdiagnosticbulk_inst%tws_grc(g) = &
         water_inst%waterdiagnosticbulk_inst%tws_grc(g) + &
         water_inst%wateratm2lndbulk_inst%volr_grc(g) / grc%area(g) * 1.e-3_r8
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

    SHR_ASSERT_ALL_FL((ubound(qflx_ice_runoff_col) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(qflx_liq_from_ice_col) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(eflx_sh_ice_to_liq_col) == (/bounds%endc/)), sourcefile, __LINE__)

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
             if (lun%itype(l) /= istice) then
                do_conversion = .true.
             else  ! istice
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
