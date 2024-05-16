module lnd_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use abortutils   , only: endrun
  use decompmod    , only: bounds_type, subgrid_level_gridcell
  use lnd2atmType  , only: lnd2atm_type
  use lnd2glcMod   , only: lnd2glc_type
  use atm2lndType  , only: atm2lnd_type
  use glc2lndMod   , only: glc2lnd_type 
  use Waterlnd2atmBulkType , only: waterlnd2atmbulk_type
  use Wateratm2lndBulkType , only: wateratm2lndbulk_type
  use clm_cpl_indices
  use GridcellType      , only : grc
  !
  implicit none
  !===============================================================================

contains

  !===============================================================================
  subroutine lnd_import( bounds, x2l, glc_present, atm2lnd_inst, glc2lnd_inst, wateratm2lndbulk_inst)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model 
    !
    ! !USES:
    use seq_flds_mod    , only: seq_flds_x2l_fields
    use clm_varctl      , only: co2_type, co2_ppmv, iulog, use_c13
    use clm_varctl      , only: ndep_from_cpl 
    use clm_varcon      , only: c13ratio
    use domainMod       , only: ldomain
    use lnd_import_export_utils, only : derive_quantities, check_for_errors, check_for_nans
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    logical            , intent(in)    :: glc_present       ! .true. => running with a non-stub GLC model
    type(atm2lnd_type) , intent(inout) :: atm2lnd_inst      ! clm internal input data type
    type(glc2lnd_type) , intent(inout) :: glc2lnd_inst      ! clm internal input data type
    type(wateratm2lndbulk_type), intent(inout) :: wateratm2lndbulk_inst   ! clm internal input data type
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg           ! bounds
    integer  :: g,i,k,nstep,ier      ! indices, number of steps, and error code
    real(r8) :: qsat_kg_kg           ! saturation specific humidity (kg/kg)
    real(r8) :: forc_pbot            ! atmospheric pressure (Pa)
    real(r8) :: forc_rainc(bounds%begg:bounds%endg)  ! rainxy Atm flux mm/s
    real(r8) :: forc_rainl(bounds%begg:bounds%endg)  ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc(bounds%begg:bounds%endg)  ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl(bounds%begg:bounds%endg)  ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag        ! temporary
    real(r8) :: co2_ppmv_prog        ! temporary
    real(r8) :: co2_ppmv_val         ! temporary
    integer  :: co2_type_idx         ! integer flag for co2_type options
    character(len=32) :: fname       ! name of field that is NaN
    character(len=32), parameter :: sub = 'lnd_import'

    !---------------------------------------------------------------------------

    ! Set bounds
    begg = bounds%begg; endg = bounds%endg

    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    do g = begg,endg
       i = 1 + (g - begg)

       ! Determine flooding input, sign convention is positive downward and
       ! hierarchy is atm/glc/lnd/rof/ice/ocn.  so water sent from rof to land is negative,
       ! change the sign to indicate addition of water to system.

       wateratm2lndbulk_inst%forc_flood_grc(g)   = -x2l(index_x2l_Flrr_flood,i)  

       wateratm2lndbulk_inst%volr_grc(g)   = x2l(index_x2l_Flrr_volr,i) * (ldomain%area(g) * 1.e6_r8)
       wateratm2lndbulk_inst%volrmch_grc(g)= x2l(index_x2l_Flrr_volrmch,i) * (ldomain%area(g) * 1.e6_r8)

       ! Determine required receive fields

       atm2lnd_inst%forc_hgt_grc(g)                  = x2l(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
       atm2lnd_inst%forc_topo_grc(g)                 = x2l(index_x2l_Sa_topo,i)      ! Atm surface height (m)
       atm2lnd_inst%forc_u_grc(g)                    = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
       atm2lnd_inst%forc_v_grc(g)                    = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
       atm2lnd_inst%forc_solad_not_downscaled_grc(g,2) = x2l(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
       atm2lnd_inst%forc_solad_not_downscaled_grc(g,1) = x2l(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
       atm2lnd_inst%forc_solai_grc(g,2)              = x2l(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
       atm2lnd_inst%forc_solai_grc(g,1)              = x2l(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

       atm2lnd_inst%forc_th_not_downscaled_grc(g)    = x2l(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
       wateratm2lndbulk_inst%forc_q_not_downscaled_grc(g)     = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
       atm2lnd_inst%forc_pbot_not_downscaled_grc(g)  = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
       atm2lnd_inst%forc_t_not_downscaled_grc(g)     = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
       atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) = x2l(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2

       forc_rainc(g)                                 = x2l(index_x2l_Faxa_rainc,i)   ! mm/s
       forc_rainl(g)                                 = x2l(index_x2l_Faxa_rainl,i)   ! mm/s
       forc_snowc(g)                                 = x2l(index_x2l_Faxa_snowc,i)   ! mm/s
       forc_snowl(g)                                 = x2l(index_x2l_Faxa_snowl,i)   ! mm/s

       ! atmosphere coupling, for prognostic/prescribed aerosols
       atm2lnd_inst%forc_aer_grc(g,1)                = x2l(index_x2l_Faxa_bcphidry,i)
       atm2lnd_inst%forc_aer_grc(g,2)                = x2l(index_x2l_Faxa_bcphodry,i)
       atm2lnd_inst%forc_aer_grc(g,3)                = x2l(index_x2l_Faxa_bcphiwet,i)
       atm2lnd_inst%forc_aer_grc(g,4)                = x2l(index_x2l_Faxa_ocphidry,i)
       atm2lnd_inst%forc_aer_grc(g,5)                = x2l(index_x2l_Faxa_ocphodry,i)
       atm2lnd_inst%forc_aer_grc(g,6)                = x2l(index_x2l_Faxa_ocphiwet,i)
       atm2lnd_inst%forc_aer_grc(g,7)                = x2l(index_x2l_Faxa_dstwet1,i)
       atm2lnd_inst%forc_aer_grc(g,8)                = x2l(index_x2l_Faxa_dstdry1,i)
       atm2lnd_inst%forc_aer_grc(g,9)                = x2l(index_x2l_Faxa_dstwet2,i)
       atm2lnd_inst%forc_aer_grc(g,10)               = x2l(index_x2l_Faxa_dstdry2,i)
       atm2lnd_inst%forc_aer_grc(g,11)               = x2l(index_x2l_Faxa_dstwet3,i)
       atm2lnd_inst%forc_aer_grc(g,12)               = x2l(index_x2l_Faxa_dstdry3,i)
       atm2lnd_inst%forc_aer_grc(g,13)               = x2l(index_x2l_Faxa_dstwet4,i)
       atm2lnd_inst%forc_aer_grc(g,14)               = x2l(index_x2l_Faxa_dstdry4,i)

       if (index_x2l_Sa_methane /= 0) then
          atm2lnd_inst%forc_pch4_grc(g) = x2l(index_x2l_Sa_methane,i)
       endif

       !--------------------------
       ! Check for nans from coupler
       !--------------------------

       call check_for_nans(x2l(:,i), fname, begg, "x2l")

    end do

    !--------------------------
    ! Derived quantities for required fields
    ! and corresponding error checks
    !--------------------------

    call derive_quantities(bounds, atm2lnd_inst, wateratm2lndbulk_inst, &
       forc_rainc, forc_rainl, forc_snowc, forc_snowl)

    call check_for_errors(bounds, atm2lnd_inst, wateratm2lndbulk_inst)

    ! Determine derived quantities for optional fields
    ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
    ! Note that forc_pbot is in Pa

    do g = begg,endg
       i = 1 + (g - begg)

       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)

       ! Determine optional receive fields
       if (index_x2l_Sa_co2prog /= 0) then
          co2_ppmv_prog = x2l(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
       else
          co2_ppmv_prog = co2_ppmv
       end if
       if (index_x2l_Sa_co2diag /= 0) then
          co2_ppmv_diag = x2l(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
       else
          co2_ppmv_diag = co2_ppmv
       end if

       if (co2_type_idx == 1) then
          co2_ppmv_val = co2_ppmv_prog
       else if (co2_type_idx == 2) then
          co2_ppmv_val = co2_ppmv_diag 
       else
          co2_ppmv_val = co2_ppmv
       end if
       if ( (co2_ppmv_val < 10.0_r8) .or. (co2_ppmv_val > 15000.0_r8) )then
          call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, &
               msg = sub//' ERROR: CO2 is outside of an expected range' )
       end if
       atm2lnd_inst%forc_pco2_grc(g) = co2_ppmv_val * 1.e-6_r8 * forc_pbot 
       if (use_c13) then
          atm2lnd_inst%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
       end if

       if (ndep_from_cpl) then
          ! The coupler is sending ndep in units if kgN/m2/s - and clm uses units of gN/m2/sec - so the
          ! following conversion needs to happen
          atm2lnd_inst%forc_ndep_grc(g) = (x2l(index_x2l_Faxa_nhx, i) + x2l(index_x2l_faxa_noy, i))*1000._r8
       end if

    end do

    call glc2lnd_inst%set_glc2lnd_fields_mct( &
         bounds = bounds, &
         glc_present = glc_present, &
         ! NOTE(wjs, 2017-12-13) the x2l argument doesn't have the typical bounds
         ! subsetting (bounds%begg:bounds%endg). This mirrors the lack of these bounds in
         ! the call to lnd_import from lnd_run_mct. This is okay as long as this code is
         ! outside a clump loop.
         x2l = x2l, &
         index_x2l_Sg_ice_covered = index_x2l_Sg_ice_covered, &
         index_x2l_Sg_topo = index_x2l_Sg_topo, &
         index_x2l_Flgg_hflx = index_x2l_Flgg_hflx, &
         index_x2l_Sg_icemask = index_x2l_Sg_icemask, &
         index_x2l_Sg_icemask_coupled_fluxes = index_x2l_Sg_icemask_coupled_fluxes)

  end subroutine lnd_import

  !===============================================================================

  subroutine lnd_export( bounds, waterlnd2atmbulk_inst, lnd2atm_inst, lnd2glc_inst, l2x)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the clm model to the coupler 
    ! 
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use seq_flds_mod       , only : seq_flds_l2x_fields
    use clm_varctl         , only : iulog
    use shr_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    use shr_fire_emis_mod  , only : shr_fire_emis_mechcomps_n
    use lnd_import_export_utils, only : check_for_nans
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds  ! bounds
    type(lnd2atm_type), intent(inout) :: lnd2atm_inst ! clm land to atmosphere exchange data type
    type(lnd2glc_type), intent(inout) :: lnd2glc_inst ! clm land to atmosphere exchange data type
    type(waterlnd2atmbulk_type), intent(in) :: waterlnd2atmbulk_inst
    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg  ! bounds
    integer  :: g,i,k ! indices
    integer  :: ier   ! error status
    integer  :: nstep ! time step index
    integer  :: dtime ! time step   
    integer  :: num   ! counter
    character(len=32) :: fname       ! name of field that is NaN
    character(len=32), parameter :: sub = 'lnd_export'
    !---------------------------------------------------------------------------

    ! Set bounds
    begg = bounds%begg; endg = bounds%endg

    ! cesm sign convention is that fluxes are positive downward

    l2x(:,:) = 0.0_r8

    do g = begg,endg
       i = 1 + (g-begg)
       l2x(index_l2x_Sl_t,i)        =  lnd2atm_inst%t_rad_grc(g)
       l2x(index_l2x_Sl_snowh,i)    =  waterlnd2atmbulk_inst%h2osno_grc(g)
       l2x(index_l2x_Sl_avsdr,i)    =  lnd2atm_inst%albd_grc(g,1)
       l2x(index_l2x_Sl_anidr,i)    =  lnd2atm_inst%albd_grc(g,2)
       l2x(index_l2x_Sl_avsdf,i)    =  lnd2atm_inst%albi_grc(g,1)
       l2x(index_l2x_Sl_anidf,i)    =  lnd2atm_inst%albi_grc(g,2)
       l2x(index_l2x_Sl_tref,i)     =  lnd2atm_inst%t_ref2m_grc(g)
       l2x(index_l2x_Sl_qref,i)     =  waterlnd2atmbulk_inst%q_ref2m_grc(g)
       l2x(index_l2x_Sl_u10,i)      =  lnd2atm_inst%u_ref10m_grc(g)
       l2x(index_l2x_Fall_taux,i)   = -lnd2atm_inst%taux_grc(g)
       l2x(index_l2x_Fall_tauy,i)   = -lnd2atm_inst%tauy_grc(g)
       l2x(index_l2x_Fall_lat,i)    = -lnd2atm_inst%eflx_lh_tot_grc(g)
       l2x(index_l2x_Fall_sen,i)    = -lnd2atm_inst%eflx_sh_tot_grc(g)
       l2x(index_l2x_Fall_lwup,i)   = -lnd2atm_inst%eflx_lwrad_out_grc(g)
       l2x(index_l2x_Fall_evap,i)   = -waterlnd2atmbulk_inst%qflx_evap_tot_grc(g)
       l2x(index_l2x_Fall_swnet,i)  =  lnd2atm_inst%fsa_grc(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          l2x(index_l2x_Fall_fco2_lnd,i) = -lnd2atm_inst%net_carbon_exchange_grc(g)  
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x(index_l2x_Sl_ram1,i)     =  lnd2atm_inst%ram1_grc(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x(index_l2x_Sl_fv,i)       =  lnd2atm_inst%fv_grc(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x(index_l2x_Sl_soilw,i)    =  waterlnd2atmbulk_inst%h2osoi_vol_grc(g,1)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x(index_l2x_Fall_flxdst1,i)= -lnd2atm_inst%flxdst_grc(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x(index_l2x_Fall_flxdst2,i)= -lnd2atm_inst%flxdst_grc(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x(index_l2x_Fall_flxdst3,i)= -lnd2atm_inst%flxdst_grc(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x(index_l2x_Fall_flxdst4,i)= -lnd2atm_inst%flxdst_grc(g,4)


       ! for dry dep velocities
       if (index_l2x_Sl_ddvel     /= 0 )  then
          l2x(index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1,i) = &
               lnd2atm_inst%ddvel_grc(g,:n_drydep)
       end if

       ! for MEGAN VOC emis fluxes
       if (index_l2x_Fall_flxvoc  /= 0 ) then
          l2x(index_l2x_Fall_flxvoc:index_l2x_Fall_flxvoc+shr_megan_mechcomps_n-1,i) = &
               -lnd2atm_inst%flxvoc_grc(g,:shr_megan_mechcomps_n)
       end if


       ! for fire emis fluxes
       if (index_l2x_Fall_flxfire  /= 0 ) then
          l2x(index_l2x_Fall_flxfire:index_l2x_Fall_flxfire+shr_fire_emis_mechcomps_n-1,i) = &
               -lnd2atm_inst%fireflx_grc(g,:shr_fire_emis_mechcomps_n)
          l2x(index_l2x_Sl_ztopfire,i) = lnd2atm_inst%fireztop_grc(g)
       end if

       if (index_l2x_Fall_methane /= 0) then
          l2x(index_l2x_Fall_methane,i) = -lnd2atm_inst%ch4_surf_flux_tot_grc(g)
       endif

       ! sign convention is positive downward with 
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  
       ! I.e. water sent from land to rof is positive

       l2x(index_l2x_Flrl_rofsur,i) = waterlnd2atmbulk_inst%qflx_rofliq_qsur_grc(g)

       !  subsurface runoff is the sum of qflx_drain and qflx_perched_drain
       l2x(index_l2x_Flrl_rofsub,i) = waterlnd2atmbulk_inst%qflx_rofliq_qsub_grc(g) &
            + waterlnd2atmbulk_inst%qflx_rofliq_drain_perched_grc(g)

       !  qgwl sent individually to coupler
       l2x(index_l2x_Flrl_rofgwl,i) = waterlnd2atmbulk_inst%qflx_rofliq_qgwl_grc(g)

       ! ice  sent individually to coupler
       l2x(index_l2x_Flrl_rofi,i) = waterlnd2atmbulk_inst%qflx_rofice_grc(g)

       ! irrigation flux to be removed from main channel storage (negative)
       l2x(index_l2x_Flrl_irrig,i) = - waterlnd2atmbulk_inst%qirrig_grc(g)

       ! glc coupling
       ! We could avoid setting these fields if glc_present is .false., if that would
       ! help with performance. (The downside would be that we wouldn't have these fields
       ! available for diagnostic purposes or to force a later T compset with dlnd.)
       do num = 0,glc_nec
          l2x(index_l2x_Sl_tsrf(num),i)   = lnd2glc_inst%tsrf_grc(g,num)
          l2x(index_l2x_Sl_topo(num),i)   = lnd2glc_inst%topo_grc(g,num)
          l2x(index_l2x_Flgl_qice(num),i) = lnd2glc_inst%qice_grc(g,num)
       end do

       !--------------------------
       ! Check for nans to coupler
       !--------------------------

       call check_for_nans(l2x(:,i), fname, begg, "l2x")

    end do

  end subroutine lnd_export

end module lnd_import_export
