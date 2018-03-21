module CNNDynamicsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for mineral nitrogen dynamics (deposition, fixation, leaching)
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use decompMod                       , only : bounds_type
  use clm_varcon                      , only : dzsoi_decomp, zisoi
  use clm_varctl                      , only : use_nitrif_denitrif, use_vertsoilc, nfix_timeconst
  use subgridAveMod                   , only : p2c
  use atm2lndType                     , only : atm2lnd_type
  use CNVegStateType                  , only : cnveg_state_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
!KO
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use TemperatureType                 , only : temperature_type
  use FrictionVelocityMod             , only : frictionvel_type
  use clm_varctl                      , only : iulog
  use shr_infnan_mod                  , only : isnan => shr_infnan_isnan
!KO
  use CNVegNitrogenStateType	      , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType	      , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use WaterStateType                  , only : waterstate_type
  use WaterFluxType                   , only : waterflux_type
  !JV
  use SoilStateType                   , only : soilstate_type
  !JV
  
  use CropType                        , only : crop_type
  use ColumnType                      , only : col                
  use PatchType                       , only : patch                
  use perf_mod                        , only : t_startf, t_stopf
  use FanMod
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNNDynamicsReadNML          ! Read in namelist for Mineral Nitrogen Dynamics
  public :: CNNDeposition               ! Update N deposition rate from atm forcing
  public :: CNNFixation                 ! Update N Fixation rate
  public :: CNNFert                     ! Update N fertilizer for crops
  public :: CNSoyfix                    ! N Fixation for soybeans
  public :: CNFreeLivingFixation        ! N free living fixation

  !
  ! !PRIVATE DATA MEMBERS:
  type, private :: params_type
     real(r8) :: freelivfix_intercept   ! intercept of line of free living fixation with annual ET
     real(r8) :: freelivfix_slope_wET   ! slope of line of free living fixation with annual ET
  end type params_type
  type(params_type) :: params_inst
  
  logical, private, parameter :: debug_fan = .false.
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNNDynamicsReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for Mineral Nitrogen Dynamics
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use abortutils     , only : endrun
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CNNDynamicsReadNML'
    character(len=*), parameter :: nmlname = 'mineral_nitrogen_dynamics'
    !-----------------------------------------------------------------------
    real(r8) :: freelivfix_intercept   ! intercept of line of free living fixation with annual ET
    real(r8) :: freelivfix_slope_wET   ! slope of line of free living fixation with annual ET
    namelist /mineral_nitrogen_dynamics/ freelivfix_slope_wET, freelivfix_intercept

    ! Initialize options to default values, in case they are not specified in
    ! the namelist


    freelivfix_intercept = 0.0117_r8
    freelivfix_slope_wET = 0.0006_r8
    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=mineral_nitrogen_dynamics, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (freelivfix_intercept, mpicom)
    call shr_mpi_bcast (freelivfix_slope_wET, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=mineral_nitrogen_dynamics)
       write(iulog,*) ' '
    end if
    params_inst%freelivfix_intercept = freelivfix_intercept
    params_inst%freelivfix_slope_wET = freelivfix_slope_wET

  end subroutine CNNDynamicsReadNML

  subroutine CNNDeposition(bounds, num_soilc, filter_soilc, &
       atm2lnd_inst, soilbiogeochem_nitrogenflux_inst, cnveg_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_carbonflux_inst, &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       waterstate_inst, soilstate_inst, temperature_inst, &
       waterflux_inst, frictionvel_inst)
    use CNSharedParamsMod    , only: use_fun
    !KO
    use clm_varctl           , only: use_fan
!   use subgridAveMod        , only: p2c
    use clm_time_manager     , only: get_step_size, get_curr_date, get_curr_calday, get_nstep

    use clm_varpar           , only: max_patch_per_col
    use LandunitType         , only: lun
    use shr_sys_mod      , only : shr_sys_flush
!KO    use ColumnType           , only: col
    use GridcellType         , only: grc
    use FanMod
    use clm_varctl     , only : iulog
    use abortutils     , only : endrun
    use pftconMod, only : nc4_grass, nc3_nonarctic_grass
    use landunit_varcon,      only:  istsoil, istcrop
    use clm_varcon, only : spval, ispval

    type(bounds_type)        , intent(in)    :: bounds  
!KO
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
!KO
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
!KO
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)         , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)          , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_carbonflux_type)   , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(waterstate_type)                  , intent(inout) :: waterstate_inst
    type(soilstate_type)                   , intent(in)    :: soilstate_inst
    type(temperature_type)                 , intent(inout) :: temperature_inst
    type(waterflux_type)                   , intent(inout) :: waterflux_inst
    type(frictionvel_type)                 , intent(inout) :: frictionvel_inst

    integer, parameter :: num_substeps = 4, balance_check_freq = 1000
    integer :: c, g, patchcounter, p, status, c1, c2, l, fc, ind_substep
    real(r8) :: dt, ndep_org(3), orgpools(3), tanprod(3), watertend, fluxes(6,3), tanpools3(3), ratm, tandep, &
         fluxes2(6,2), fluxes3(6,3), tanpools2(2), fluxes_tmp(6), garbage_total
    real(r8), parameter :: water_init_grz = 0.005_r8, cnc_nh3_air = 0.0_r8, depth_slurry = 0.005_r8
    real(r8), parameter :: fract_resist=0.225_r8, fract_unavail=0.025_r8, fract_avail=0.25_r8, fract_tan=0.5_r8
    real(r8), parameter :: dz_layer_fert = 0.02_r8, dz_layer_grz = 0.02_r8
    !real(r8), parameter :: fract_resist=0._r8, fract_unavail=0._r8, fract_avail=0._r8, fract_tan=1.0_r8
    
    real(r8), parameter :: slurry_infiltr_time = 6*3600.0_r8, water_init_fert = 1e-6
    real(r8), parameter :: poolranges_grz(2) = (/24*3600.0_r8, 360*24*3600.0_r8/), &
         poolranges_fert(3) = (/2*24*3600.0_r8, 24*3600.0_r8, 360*24*3600.0_r8/), &
         poolranges_slr(3) = (/slurry_infiltr_time, 24*3600.0_r8, 360*24*3600.0_r8/), &
         Hconc_grz(2) = (/10**(-8.5_r8), 10**(-8.0_r8)/), &
         Hconc_fert(3) = (/10**-7.0_r8, 10**(-8.0_r8), 10**(-8.0_r8)/)
    !logical, parameter :: do_balance_checks = .false.
    logical :: do_balance_checks
    real(r8) :: tg, garbage, theta, thetasat, infiltr_m_s, evap_m_s, runoff_m_s, org_n_tot, &
         nstored_old, nsoilman_old, nsoilfert_old, fert_to_air, fert_to_soil, fert_total, fert_urea, fert_tan, &
         soilflux_org, urea_resid
    real(r8) :: tanprod_from_urea(3), ureapools(2), fert_no3, fert_generic
    real(r8), parameter :: fract_urea=0.545, fract_no3=0.048
    integer, parameter :: ind_region = 1
    
    dt = real( get_step_size(), r8 )
    do_balance_checks = mod(get_nstep(), balance_check_freq) == 0
    associate(                                                                &
         ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)
         forc_ndep     =>  atm2lnd_inst%forc_ndep_grc ,                       &
         ! Output: [real(r8) (:)]  atmospheric N deposition to soil mineral N (gN/m2/s)
         ndep_to_sminn =>  soilbiogeochem_nitrogenflux_inst%ndep_to_sminn_col & 
         )
      ! Loop through columns
      do c = bounds%begc, bounds%endc
         g = col%gridcell(c)
         ndep_to_sminn(c) = forc_ndep(g)
      end do
    end associate
    
    associate(&
         ngrz => soilbiogeochem_nitrogenflux_inst%man_n_grz_col, &
         man_u_grz => soilbiogeochem_nitrogenstate_inst%man_u_grz_col, &
         man_a_grz => soilbiogeochem_nitrogenstate_inst%man_a_grz_col, &
         man_r_grz => soilbiogeochem_nitrogenstate_inst%man_r_grz_col, &
         man_u_app => soilbiogeochem_nitrogenstate_inst%man_u_app_col, &
         man_a_app => soilbiogeochem_nitrogenstate_inst%man_a_app_col, &
         man_r_app => soilbiogeochem_nitrogenstate_inst%man_r_app_col, &
         ns => soilbiogeochem_nitrogenstate_inst, &
         nf => soilbiogeochem_nitrogenflux_inst, &
         cnv_nf => cnveg_nitrogenflux_inst, &
         ram1 => frictionvel_inst%ram1_patch, &
         rb1 => frictionvel_inst%rb1_patch)

    nf%fert_n_appl_col(bounds%begc:bounds%endc) = 0.0
    nf%man_n_appl_col(bounds%begc:bounds%endc) = 0.0
    nf%man_tan_appl_col(bounds%begc:bounds%endc) = 0.0

    call p2c(bounds, num_soilc, filter_soilc, &
         cnv_nf%fert_patch(bounds%begp:bounds%endp), &
         nf%fert_n_appl_col(bounds%begc:bounds%endc))
    call p2c(bounds, num_soilc, filter_soilc, &
         cnv_nf%manu_patch(bounds%begp:bounds%endp), &
         nf%man_n_appl_col(bounds%begc:bounds%endc))

    if (any(nf%man_n_appl_col > 100)) then
       write(iulog, *) maxval(nf%man_n_appl_col)
       call endrun('bad man_n_appl_col')
    end if
    if (do_balance_checks) then
       nstored_old = get_total_n(ns, nf, 'pools_storage')
       nsoilman_old = get_total_n(ns, nf, 'pools_manure')
       nsoilfert_old = get_total_n(ns, nf, 'pools_fertilizer')
    end if

    ! Assign the "pastoral" manure entire to the natural vegetation column
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       l = col%landunit(c)
       if (.not. (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)) cycle
       if (.not. col%active(c) .or. col%wtgcell(c) < 1e-6) cycle
       g = col%gridcell(c)
       if (lun%itype(l) == istsoil) then
          ngrz(c) = atm2lnd_inst%forc_ndep3_grc(g) / col%wtgcell(c) * 1e3 ! kg to g 
          if (debug_fan) then
             if (ngrz(c) > 1e12 .or. (isnan(ngrz(c)))) then
                write(iulog, *) 'bad ngrz', atm2lnd_inst%forc_ndep3_grc(g), col%wtgcell(c)
                call endrun('bad ngrz 1')
             end if
          end if
          if (nf%man_n_appl_col(c) > 0) then
             write(iulog, *) nf%man_n_appl_col(c)
             call endrun(msg='Found fertilizer in soil column')
          end if
       else
          ngrz(c) = 0.0
       end if

    end do

    if(debug_fan) then
       write(iulog, *) 'nan count of storage 1', count(isnan(ns%man_n_stored_col))
       if (any(isnan(nf%man_n_appl_col))) then
          call endrun('nan nh3 appl b')
       end if
    end if

    call handle_storage(bounds, temperature_inst, frictionvel_inst, dt, &
         atm2lnd_inst%forc_ndep2_grc, &
         ns%man_n_stored_col, ns%man_tan_stored_col, &
         nf%man_n_appl_col, nf%man_tan_appl_col, &
         nf%man_n_grz_col, nf%man_n_mix_col, &
         nf%nh3_stores_col, nf%nh3_barns_col, &
         nf%man_n_transf_col, filter_soilc, num_soilc)

    if (debug_fan) then
       if (any(isnan(nf%nh3_stores_col))) then
          call endrun('nan nh3 stores')
       end if
       if (any(isnan(nf%nh3_barns_col))) then
          call endrun('nan nh3 barns')
       end if
       if (any(isnan(nf%man_n_appl_col))) then
          call endrun('nan nh3 appl')
       end if
       if (any(isnan(nf%man_n_mix_col))) then
          call endrun('nan nh3 appl')
       end if
    end if

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       l = col%landunit(c)
       if (.not. (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)) cycle
       if (.not. col%active(c) .or. col%wtgcell(c) < 1e-6) cycle

       if (nf%man_n_appl_col(c) > 1e12 .or. ngrz(c) > 1e12) then
          write(iulog, *) c, nf%man_n_appl_col(c), ngrz(c), cnv_nf%fert_patch(col%patchi(c):col%patchf(c)), &
               cnv_nf%manu_patch(col%patchi(c):col%patchf(c))
          call endrun('nf%man_n_appl_col(c) is spval')
       end if

       ! Find and average the atmospheric resistances Rb and Ra.
       ! 
       if (lun%itype(col%landunit(c)) == istcrop) then
          ! for column, only one patch
          p = col%patchi(c)
          if (p /= col%patchf(c)) call endrun(msg='Strange patch for crop')
          ratm = ram1(p) + rb1(p)
       else
          ! if natural, find average over grasses
          ratm = 0.0
          patchcounter = 0
          do p = col%patchi(c), col%patchf(c)
             if (patch%itype(p) == nc4_grass .or. patch%itype(p) == nc3_nonarctic_grass) then
                if (.not. patch%active(p) .or. ram1(p) == spval .or. rb1(p) == spval) cycle
                ratm = ratm + ram1(p) + rb1(p)
                patchcounter = patchcounter + 1
             end if
          end do
          if (patchcounter > 0) then
             ratm = ratm / patchcounter
          else
             ! grass not found, take something.
             do p = col%patchi(c), col%patchf(c)
                if (.not. patch%active(p) .or. ram1(p) == spval .or. rb1(p) == spval) cycle
                ratm = ram1(p) + rb1(p)
                exit
             end do
             if (p == col%patchf(c) + 1) then
                ratm = 150.0_r8
             end if
          end if
       end if

       ! Calculation of the water fluxes should include the background soil moisture
       ! tendency. However, it's unclear how to do this in a numerically consistent
       ! way. Following a naive finite differencing approach led to poorer agreement in
       ! stand-alone simulations so the term is currenltly neglected here.
       watertend = 0.0_r8 
       tg = temperature_inst%t_grnd_col(c)
       theta = waterstate_inst%h2osoi_vol_col(c,1)
       thetasat = soilstate_inst%watsat_col(c,1)
       theta = min(theta, 0.98_r8*thetasat)
       infiltr_m_s = max(waterflux_inst%qflx_infl_col(c), 0.0) * 1e-3 
       evap_m_s = waterflux_inst%qflx_evap_grnd_col(c) * 1e-3
       runoff_m_s = max(waterflux_inst%qflx_runoff_col(c), 0.0) * 1e-3

       !
       ! grazing
       !

       ndep_org(ind_avail) = ngrz(c) * fract_avail
       ndep_org(ind_resist) = ngrz(c) * fract_resist
       ndep_org(ind_unavail) = ngrz(c) * fract_unavail
       tandep = ngrz(c) * fract_tan

       orgpools(ind_avail) = man_a_grz(c)
       orgpools(ind_resist) = man_r_grz(c)
       orgpools(ind_unavail) = man_u_grz(c)
       call update_org_n(ndep_org, tg, orgpools, dt, tanprod, soilflux_org)
       man_a_grz(c) = orgpools(ind_avail)
       man_r_grz(c) = orgpools(ind_resist) 
       man_u_grz(c) = orgpools(ind_unavail)

       tanpools2(1) = ns%tan_g1_col(c)
       tanpools2(2) = ns%tan_g2_col(c)
       if (any(isnan(tanpools2))) then
          call endrun('nan1')
       end if

       fluxes_tmp = 0.0
       garbage_total = 0.0
       fluxes2 = 0.0
       garbage = 0
       do ind_substep = 1, num_substeps
          call update_npool(tg, ratm, &
               theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, tandep, (/0.0_r8, sum(tanprod)/), water_init_grz, &
               cnc_nh3_air, poolranges_grz, Hconc_grz, dz_layer_grz, tanpools2, &
               fluxes2(1:5,:), garbage, dt/num_substeps, status, 2)
          if (status /= 0) then
             write(iulog, *) 'status = ', status, tanpools2, ratm, theta, thetasat, tandep, tanprod
             call endrun(msg='update_npool status /= 0')
          end if
          if (debug_fan .and. any(isnan(tanpools2))) then
             call endrun('nan2')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes2, dim=2)
          garbage_total = garbage_total + garbage
       end do
       fluxes_tmp = fluxes_tmp / num_substeps

       ns%tan_g1_col(c) = tanpools2(1)
       ns%tan_g2_col(c) = tanpools2(2)
       if (debug_fan .and. any(isnan(fluxes2))) then
          write(iulog, *) fluxes2
          call endrun('nan3')
       end if

       nf%nh3_grz_col(c) = fluxes_tmp(iflx_air)
       nf%manure_runoff_col(c) = fluxes_tmp(iflx_roff)
       nf%manure_no3_prod_col(c) = fluxes_tmp(iflx_no3)
       nf%manure_nh4_to_soil_col(c) &
            = fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) + garbage_total / dt + soilflux_org

       !
       ! Manure application
       !

       org_n_tot = nf%man_n_appl_col(c) - nf%man_tan_appl_col(c)
       ! Use the the same fractionation of organic N as for grazing, after removing the
       ! "explicitly" calculated TAN.
       if (1-fract_tan > 1e-6) then
          ndep_org(ind_avail) = org_n_tot * fract_avail / (1-fract_tan)
          ndep_org(ind_resist) = org_n_tot * fract_resist / (1-fract_tan)
          ndep_org(ind_unavail) = org_n_tot * fract_unavail / (1-fract_tan)
       else
          ndep_org = 0.0
       end if
       tandep = nf%man_tan_appl_col(c)

       orgpools(ind_avail) = man_a_app(c)
       orgpools(ind_resist) = man_r_app(c)
       orgpools(ind_unavail) = man_u_app(c)
       call update_org_n(ndep_org, tg, orgpools, dt, tanprod, soilflux_org)
       man_a_app(c) = orgpools(ind_avail)
       man_r_app(c) = orgpools(ind_resist)
       man_u_app(c) = orgpools(ind_unavail)
       tanpools3(1) = ns%tan_s0_col(c)
       tanpools3(2) = ns%tan_s1_col(c)
       tanpools3(3) = ns%tan_s2_col(c)

       if (debug_fan .and. any(isnan(tanpools3))) then
          call endrun('nan31')
       end if

       fluxes_tmp = 0.0
       garbage_total = 0.0
       fluxes3 = 0.0
       do ind_substep = 1, num_substeps
          if (debug_fan .and. any(abs(tanpools3) > 1e12)) then
             write(iulog, *) ind_substep, tanpools3, tandep, nf%fert_n_appl_col(c), &
                  nf%man_n_appl_col(c), ns%man_n_stored_col(c), ns%man_tan_stored_col(c)
             call endrun('bad tanpools (manure app)')
          end if

          call update_3pool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, tandep, sum(tanprod), cnc_nh3_air, depth_slurry, &
               poolranges_slr, tanpools3, fluxes3(1:5,:), garbage, dt / num_substeps, status)
          if (status /= 0) then
             write(iulog, *) 'status = ', status, tanpools3, tg, ratm, 'th', theta, &
                  thetasat, tandep, 'tp', tanprod, 'fx', fluxes
             call endrun(msg='update_3pool status /= 0')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes3, dim=2)
          garbage_total = garbage_total + garbage
       end do
       fluxes_tmp = fluxes_tmp / num_substeps

       ns%tan_s0_col(c) = tanpools3(1)
       ns%tan_s1_col(c) = tanpools3(2)
       ns%tan_s2_col(c) = tanpools3(3)

       if (debug_fan .and. any(isnan(fluxes3))) then
          write(iulog, *) fluxes3, tanpools3,ratm, theta, thetasat, infiltr_m_s, tandep, tanprod
          call endrun('nan4')
       end if

       nf%nh3_man_app_col(c) = fluxes_tmp(iflx_air)
       nf%manure_runoff_col(c) = nf%manure_runoff_col(c) + fluxes_tmp(iflx_roff)
       nf%manure_no3_prod_col(c) = nf%manure_no3_prod_col(c) + fluxes_tmp(iflx_no3)
       nf%manure_nh4_to_soil_col(c) &
            = nf%manure_nh4_to_soil_col(c) + fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) &
            + garbage_total / dt + soilflux_org

       !
       ! Fertilizer
       !

       fert_total = nf%fert_n_appl_col(c)

       fert_urea = fert_total * fract_urea
       fert_no3 = fert_total * fract_no3
       fert_generic = fert_total - fert_urea - fert_no3

       ! Urea decomposition 
       ! 
       ureapools(1) = ns%fert_u0_col(c)
       ureapools(2) = ns%fert_u1_col(c)
       fluxes2 = 0.0
       call update_urea(tg, theta, thetasat, infiltr_m_s, evap_m_s, watertend, &
            runoff_m_s, fert_urea, ureapools,  fluxes2, urea_resid, poolranges_fert(1:2), &
            dt, status, numpools=2)
       if (status /= 0) then
          call endrun(msg='Bad status after update_urea for fertilizer')
       end if
       ! Nitrogen fluxes from urea pool. Be sure to not zero below!
       fluxes_tmp = sum(fluxes2, dim=2)

       ns%fert_u0_col(c) = ureapools(1)
       ns%fert_u1_col(c) = ureapools(2)
       ! Collect the formed ammonia for updating the TAN pools
       tanprod_from_urea(1:2) = fluxes2(iflx_to_tan, 1:2)
       tanprod_from_urea(2) = tanprod_from_urea(2)
       ! There is no urea pool corresponding to tan_f2, because most of the urea will
       ! have decomposed. Here whatever remains gets sent to tan_f2. 
       tanprod_from_urea(3) = urea_resid / dt 

       tanpools3(1) = ns%tan_f0_col(c)
       tanpools3(2) = ns%tan_f1_col(c)
       tanpools3(3) = ns%tan_f2_col(c)         
       garbage_total = 0.0
       fluxes3 = 0.0
       do ind_substep = 1, num_substeps
          ! Fertilizer pools f0...f2
          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, 0.0_r8, tanprod_from_urea, water_init_fert, cnc_nh3_air, &
               poolranges_fert, Hconc_fert, dz_layer_fert, tanpools3, fluxes3(1:5,:), &
               garbage, dt/num_substeps, status, numpools=3)
          if (status /= 0) then
             write(iulog, *) 'status:', status, tanpools3, nf%fert_n_appl_col(c)
             call endrun(msg='Bad status after npool for fertilizer')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes3, dim=2) / num_substeps
          garbage_total = garbage_total + garbage

          ! Fertilizer pool f3
          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, fert_generic, (/0.0_r8/), water_init_fert, cnc_nh3_air, &
               (/360*24*3600.0_r8/), (/10**(-6.5_r8)/), dz_layer_fert, ns%tan_f3_col(c:c), fluxes3(1:5,1:1), &
               garbage, dt/num_substeps, status, numpools=1)
          if (status /= 0) then
             write(iulog, *) 'status:', status, tanpools3, nf%fert_n_appl_col(c)
             call endrun(msg='Bad status after npool for generic')
          end if
          fluxes_tmp = fluxes_tmp + fluxes3(:, 1) / num_substeps
          garbage_total = garbage_total + garbage
       end do

       ns%tan_f0_col(c) = tanpools3(1)
       ns%tan_f1_col(c) = tanpools3(2)
       ns%tan_f2_col(c) = tanpools3(3)
       ! !!tan_f3_col already updated above by update_npool!!

       nf%nh3_fert_col(c) = fluxes_tmp(iflx_air)
       nf%fert_runoff_col(c) = fluxes_tmp(iflx_roff)
       nf%fert_no3_prod_col(c) = fluxes_tmp(iflx_no3) + fert_no3
       nf%fert_nh4_to_soil_col(c) &
            = fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) + fert_to_soil + garbage_total/dt 

       ! Total flux
       ! 
       nf%nh3_total_col(c) = nf%nh3_fert_col(c) + nf%nh3_man_app_col(c) &
            + nf%nh3_grz_col(c) + nf%nh3_stores_col(c) +  nf%nh3_barns_col(c)
       if (nh%nh3_total_col(c) < -1e15) then
          call endrun(msg='ERROR: FAN, negative total emission')
       end if
    end do

    if (do_balance_checks) then
       call balance_check('Storage', nstored_old, &
            get_total_n(ns, nf, 'pools_storage'), get_total_n(ns, nf, 'fluxes_storage'))
       call balance_check('Manure', nsoilman_old, &
            get_total_n(ns, nf, 'pools_manure'), get_total_n(ns, nf, 'fluxes_manure'))
       call balance_check('Fertilizer', nsoilfert_old, &
            get_total_n(ns, nf, 'pools_fertilizer'), get_total_n(ns, nf, 'fluxes_fertilizer'))
    end if

    end associate

  contains

    real(r8) function get_total_n(ns, nf, which) result(total)
      type(soilbiogeochem_nitrogenstate_type), intent(in) :: ns
      type(soilbiogeochem_nitrogenflux_type), intent(in) :: nf
      character(len=*), intent(in) :: which
      
      total = 0

      associate(soilc => filter_soilc(1:num_soilc))
        
      select case(which)
      case('pools_storage')
         total = sum(ns%man_n_stored_col(soilc))
         
      case('fluxes_storage')
         total = sum(nf%man_n_mix_col(soilc))
         total = total - sum(nf%nh3_stores_col(soilc))
         total = total - sum(nf%nh3_barns_col(soilc)) - sum(nf%man_n_transf_col(soilc))
         
      case('pools_manure')
         total = total + sum(ns%tan_g1_col(soilc)) + sum(ns%tan_g2_col(soilc)) 
         total = total + sum(ns%man_u_grz_col(soilc)) &
              + sum(ns%man_a_grz_col(soilc)) + sum(ns%man_r_grz_col(soilc))
         total = total + sum(ns%tan_s0_col(soilc)) &
              + sum(ns%tan_s1_col(soilc)) + sum(ns%tan_s2_col(soilc))
         total = total + sum(ns%man_u_app_col(soilc)) &
              + sum(ns%man_a_app_col(soilc)) + sum(ns%man_r_app_col(soilc))
         
      case('fluxes_manure')
         total = sum(nf%man_n_grz_col(soilc)) + sum(nf%man_n_appl_col(soilc)) 
         total = total - sum(nf%nh3_man_app_col(soilc)) &
              - sum(nf%nh3_grz_col(soilc)) - sum(nf%manure_runoff_col(soilc))
         total = total - sum(nf%manure_no3_prod_col(soilc)) - sum(nf%manure_nh4_to_soil_col(soilc))
         
      case('pools_fertilizer')
         total = sum(ns%tan_f0_col((soilc))) + sum(ns%tan_f1_col((soilc))) + sum(ns%tan_f2_col(soilc)) &
              + sum(ns%tan_f3_col(soilc))
         total = total + sum(ns%fert_u0_col(soilc)) + sum(ns%fert_u1_col(soilc))
         
      case('fluxes_fertilizer')
         total = sum(nf%fert_n_appl_col(soilc))
         total = total - sum(nf%nh3_fert_col(soilc)) - sum(nf%fert_runoff_col(soilc))
         total = total - sum(nf%fert_no3_prod_col(soilc)) - sum(nf%fert_nh4_to_soil_col(soilc))
         
      case default
         call endrun(msg='Bad argument to get_total_n')
         
      end select
      
      end associate

    end function get_total_n
    
    subroutine balance_check(label, total_old, total_new, flux)
      ! Check and report that the net flux equals the accumulated mass in pools. The
      ! total pools and fluxes can be evaluated by the function get_total_n.
      character(len=*), intent(in) :: label
      real(r8), intent(in) :: total_old, total_new, flux

      real(r8) :: diff, accflux
      real(r8) :: tol = 1e-6_r8
      
      diff = total_new - total_old
      accflux = flux*dt
      write(iulog, *) 'Balance check:', label, diff, accflux
      
    end subroutine balance_check
    
  end subroutine CNNDeposition

  
  subroutine handle_storage(bounds, temperature_inst, frictionvel_inst, dt,  &
       ndep_mixed_grc, n_stored_col, tan_stored_col, &
       n_manure_spread_col, tan_manure_spread_col, &
       n_manure_graze_col, n_manure_mixed_col, &
       nh3_flux_stores, nh3_flux_barns, man_n_transf, &
       filter_soilc, num_soilc)
    use landunit_varcon, only : max_lunit
    use pftconMod, only : nc4_grass, nc3_nonarctic_grass
    use clm_varcon, only : ispval
    use landunit_varcon,      only:  istsoil, istcrop
    use abortutils     , only : endrun
    use LandunitType   , only: lun
    use GridcellType   , only: grc
    use clm_varctl     , only : iulog
    use ColumnType     , only : col                

    implicit none
    type(bounds_type), intent(in)    :: bounds
    type(temperature_type) , intent(in) :: temperature_inst
    type(frictionvel_type) , intent(in) :: frictionvel_inst
    real(r8), intent(in) :: dt
    
    ! N excreted in manure, mixed/pastoral systems, gN/m2:
    real(r8), intent(in) :: ndep_mixed_grc(bounds%begg:bounds%endg)
    real(r8), intent(inout) :: n_stored_col(bounds%begc:bounds%endc), tan_stored_col(bounds%begc:bounds%endc) ! N, TAN currently stored, gN/m2
    ! N, TAN spread on grasslands, gN/m2/s:
    real(r8), intent(inout) :: n_manure_spread_col(bounds%begc:bounds%endc) ! for crops, input, determined by crop model, otherwise output
    real(r8), intent(out) :: tan_manure_spread_col(bounds%begc:bounds%endc) ! output, calculated from the above and stored manure
    ! N excreted by animals allocated to mixed production systems temporarily grazing on grasslands:
    real(r8), intent(inout) :: n_manure_graze_col(bounds%begc:bounds%endc)
    ! N excreted by animals in mixed systems, total
    real(r8), intent(out) :: n_manure_mixed_col(bounds%begc:bounds%endc)
    ! NH3 emission fluxes from manure storage and housings, gN/m2/s
    real(r8), intent(out) :: nh3_flux_stores(bounds%begc:bounds%endc), nh3_flux_barns(bounds%begc:bounds%endc)
    ! total nitrogen flux transferred out of a crop column
    real(r8), intent(out) :: man_n_transf(bounds%begc:bounds%endc)
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns

    integer :: begg, endg, g, l, c, il, counter, col_grass, status, p
    real(r8) :: flux_avail, flux_grazing
    real(r8) :: tempr_ave, windspeed_ave ! windspeed and temperature averaged over agricultural patches
    real(r8) :: tempr_barns, tempr_stores, vent_barns, flux_grass_crop, tempr_min_10day, &
         flux_grass_graze, flux_grass_spread, flux_grass_spread_tan, flux_grass_crop_tan
    real(r8) :: cumflux, totalinput
    real(r8) :: fluxes_nitr(4), fluxes_tan(4)
    ! The fraction of manure applied continuously on grasslands (if present in the gridcell)
    real(r8), parameter :: fract_continuous = 0.1_r8, kg_to_g = 0.001_r8, max_grazing_fract = 0.3_r8, &
         tan_fract_excr = 0.5_r8, volat_coef_barns = 0.02_r8, volat_coef_stores = 0.01_r8, &
         tempr_min_grazing = 283.0_r8!!!!

    begg = bounds%begg; endg = bounds%endg
    nh3_flux_stores(bounds%begc:bounds%endc) = 0_r8
    nh3_flux_barns(bounds%begc:bounds%endc) = 0_r8

    totalinput = 0.0
    cumflux = 0.0
    
    do g = begg, endg
       !totalinput = totalinput + ndep_mixed_grc(g)
       
       ! First find out if there are grasslands in this cell. If yes, a fraction of
       ! manure can be diverted to them before storage.
       col_grass = ispval
       do il = 1, max_lunit
          l = grc%landunit_indices(il, g)
          if (lun%itype(l) == istsoil) then
             do p = lun%patchi(l), lun%patchf(l)
                if (patch%itype(p) == nc4_grass .or. patch%itype(p) == nc3_nonarctic_grass) then
                   col_grass = patch%column(p)
                   exit
                end if
             end do
          end if
          if (col_grass /= ispval) exit
       end do
       if (col%wtgcell(col_grass) < 1e-6) col_grass = ispval
       ! Transfer of manure from all crop columns to the natural vegetation column:
       flux_grass_graze = 0_r8
       flux_grass_spread = 0_r8
       flux_grass_spread_tan = 0_r8

       do il = 1, max_lunit
          l = grc%landunit_indices(il, g)
          if (l == ispval) cycle
          if (lun%itype(l) == istcrop) then
             ! flux_avail = manure excreted per m2 of crops (ndep_mixed_grc = per m2 / all land units)
             do c = lun%coli(l), lun%colf(l)
                if (.not. col%active(c)) cycle
                if (col%wtgcell(c) < 1e-6) cycle

                if (col%landunit(c) /= l) then
                   write(iulog, *) g, il, c, col%landunit(c)
                   call endrun('something wrong')
                end if
                if (.not. any(c==filter_soilc(1:num_soilc))) then
                   write(iulog, *) c, n_manure_spread_col(c)
                   call endrun('column not in soilfilter')
                end if

                flux_avail = ndep_mixed_grc(g) * kg_to_g / lun%wtgcell(l)
                if (flux_avail > 1e12 .or. isnan(flux_avail)) then
                   write(iulog, *) 'bad flux_avail', ndep_mixed_grc(g), lun%wtgcell(l)
                   call endrun('bad flux_avail')
                end if
                n_manure_mixed_col(c) = flux_avail
                totalinput = totalinput + flux_avail
                !manure_input(c) = flux_avail

                !tempr_ave = 0_r8
                !windspeed_ave = 0_r8
                counter = 0
                if (col_grass == c) call endrun('Something wrong with the indices')
                if (col%patchi(c) /= col%patchf(c)) then
                   call endrun(msg="ERROR crop column has multiple patches")
                end if

                tempr_ave = temperature_inst%t_ref2m_patch(col%patchi(c))
                windspeed_ave = frictionvel_inst%u10_patch(col%patchi(c))

                tempr_min_10day = temperature_inst%t_a10min_patch(col%patchi(c))
                if (tempr_min_10day > tempr_min_grazing) then
                   ! fraction of animals grazing -> allocate some manure to grasslands before barns
                   flux_grazing = max_grazing_fract * flux_avail
                   flux_avail = flux_avail - flux_grazing
                else
                   flux_grazing = 0_r8
                end if
                flux_grass_graze = flux_grass_graze + flux_grazing*col%wtgcell(c)

                call eval_fluxes_storage(flux_avail, tempr_ave, windspeed_ave, fract_continuous, &
                     volat_coef_barns, volat_coef_stores, tan_fract_excr, fluxes_nitr, fluxes_tan, status)
                if (any(fluxes_nitr > 1e12)) then
                   write(iulog, *) 'bad fluxes', fluxes_nitr
                end if
                if (status /=0) then 
                   write(iulog, *) 'status = ', status
                   call endrun(msg='eval_fluxes_storage failed')
                end if
                cumflux = cumflux + sum(fluxes_nitr)
                
                !flux_grass_spread = flux_grass_spread + flux_grass_crop*col%wtgcell(c)
                flux_grass_spread = flux_grass_spread + fluxes_nitr(iflx_appl)*col%wtgcell(c)
                !flux_grass_spread_tan = flux_grass_spread_tan + flux_grass_crop_tan*col%wtgcell(c)
                flux_grass_spread_tan = flux_grass_spread_tan + fluxes_tan(iflx_appl)*col%wtgcell(c)

                
                !man_n_transf(c) = flux_grazing

                if (fluxes_tan(iflx_to_store) < 0) then
                   call endrun(msg="ERROR too much manure lost")
                end if

                if (n_stored_col(c) < 0) then
                   call endrun(msg='n_stored_col is negative')
                end if
                
                if (n_stored_col(c) > 0_r8) then
                   tan_manure_spread_col(c) = n_manure_spread_col(c) * tan_stored_col(c)/n_stored_col(c)
                else if (n_manure_spread_col(c) > 1e-15_r8) then
                   write(iulog, *) 'stored, spread', n_stored_col(c), n_manure_spread_col(c)
                   call endrun(msg='Inconsistent manure application')
                else
                   tan_manure_spread_col(c) = 0_r8
                end if

                if (tan_manure_spread_col(c) > 1) then
                   write(iulog, *) 'bad tan_manure', tan_manure_spread_col(c), tan_stored_col(c), n_stored_col(c), n_manure_spread_col(c)
                end if
                

                if (n_manure_spread_col(c) > 1) then
                   write(iulog, *) 'bad n_manure', tan_manure_spread_col(c), tan_stored_col(c), n_stored_col(c), n_manure_spread_col(c)
                end if

                n_stored_col(c) = n_stored_col(c) + (fluxes_nitr(iflx_to_store) - n_manure_spread_col(c)) * dt
                tan_stored_col(c) = tan_stored_col(c) &
                     + (fluxes_tan(iflx_to_store) - tan_manure_spread_col(c)) * dt
                if (n_stored_col(c) > 1e6) then
                   call endrun(msg='ERROR bad n_stored_col')
                end if
                if (n_stored_col(c) < 0) then
                   if (n_stored_col(c) > -1e-6_r8) then
                      n_stored_col(c) = 0_r8
                   else
                      call endrun(msg="ERROR negative n_stored_col")
                   end if
                end if

                man_n_transf(c) = fluxes_nitr(iflx_appl) + flux_grazing + n_manure_spread_col(c)
                
                nh3_flux_stores(c) = fluxes_nitr(iflx_air_stores)
                nh3_flux_barns(c) = fluxes_nitr(iflx_air_barns)
                
             end do ! column
          end if ! crop land unit
       end do ! landunit

       if (col_grass /= ispval) then
          if (tan_manure_spread_col(col_grass) > 1) then
             write(iulog, *) 'bad tan_manure col_grass before adding', n_manure_spread_col(col_grass), &
                  tan_manure_spread_col(col_grass)
          end if
          n_manure_spread_col(col_grass) = n_manure_spread_col(col_grass) &
               + flux_grass_spread / col%wtgcell(col_grass)
          tan_manure_spread_col(col_grass) = tan_manure_spread_col(col_grass) &
               + flux_grass_spread_tan / col%wtgcell(col_grass)
          n_manure_graze_col(col_grass) = n_manure_graze_col(col_grass) + flux_grass_graze / col%wtgcell(col_grass)
          !write(iulog, *) 'to grass:', n_manure_spread(col_grass), col_grass
          if (tan_manure_spread_col(col_grass) > 1) then
             write(iulog, *) 'bad tan_manure col_grass', flux_grass_spread_tan, col%wtgcell(col_grass)
          end if

       end if

    end do ! grid

  end subroutine handle_storage

  
  !-----------------------------------------------------------------------
!KO  subroutine CNNDeposition( bounds, &
!KO       atm2lnd_inst, soilbiogeochem_nitrogenflux_inst )
!KO
  subroutine CNNDeposition_old( bounds, num_soilc, filter_soilc, &
       atm2lnd_inst, soilbiogeochem_nitrogenflux_inst, cnveg_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_carbonflux_inst, &
       cnveg_nitrogenstate_inst, waterstate_inst, temperature_inst, &
       waterflux_inst, frictionvel_inst )
!KO
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen deposition rate
    ! from atmospheric forcing. For now it is assumed that all the atmospheric
    ! N deposition goes to the soil mineral N pool.
    ! This could be updated later to divide the inputs between mineral N absorbed
    ! directly into the canopy and mineral N entering the soil pool.
    !
    ! !USES:
    use CNSharedParamsMod    , only: use_fun
!KO
    use clm_varctl           , only: use_fan
!   use subgridAveMod        , only: p2c
    use clm_time_manager     , only: get_step_size, get_curr_date, get_curr_calday
    use clm_varpar           , only: max_patch_per_col
    use LandunitType         , only: lun
    use shr_sys_mod      , only : shr_sys_flush
!KO    use ColumnType           , only: col
    use GridcellType         , only: grc
!KO    use PatchType            , only: patch
!KO
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
!KO
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
!KO
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
!KO
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)         , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenstate_type), intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_carbonflux_type)   , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(waterstate_type)                  , intent(inout) :: waterstate_inst
    type(temperature_type)                 , intent(inout) :: temperature_inst
    type(waterflux_type)                   , intent(inout) :: waterflux_inst
    type(frictionvel_type)                 , intent(inout) :: frictionvel_inst
!KO
    !
    ! !LOCAL VARIABLES:
!KO    integer :: g,c                                    ! indices
!KO
    integer :: g,c,p,fc,pi                            ! Indices
    integer :: yr, mon, day, sec                      ! Outputs from get_curr_date
    real(r8) :: fmm                                   ! Fraction of manure that changes to methane (0.055; Lerner et al 1988)
    real(r8) :: na                                    ! Fraction of N assimilated by cows (0)
    real(r8) :: ca                                    ! Fraction of C assimilated by cows (0)
    real(r8) :: pl                                    ! Fraction of cows feed consumed from live shoots (0.9; Holland et al, 1992)
    real(r8) :: pd                                    ! Fraction of cows feed consumed from dead shoots (0.1; Holland et al, 1992)
    real(r8) :: cn                                    ! Carbon:nitrogen ratio (30)
    real(r8) :: fua                                   ! Fraction of urine that volatilizes (0.2; Holland et al, 1992)
    real(r8) :: fu                                    ! Fraction of manure that is urine (0.5; Parton et al, 2001)
    real(r8) :: ffa                                   ! Fraction of feces that volatilizes pre-industrial (0.27; Bouwman et al, 2002)
    real(r8) :: ff                                    ! Fraction of manure that is feces (0.5; Parton et al, 2001)
    real(r8) :: f1850                                 ! Fraction of ndep_fert in 1850 (0.19; Potter et al, 2010)
    real(r8) :: m1850                                 ! Fraction of ndep_manure in 1850 (0.19; Potter et al, 2010)
    real(r8) :: dt                                    ! Radiation time step (seconds)
    real(r8) :: NFactor                               !KO
    real(r8) :: denit                                 !KO
    real(r8) :: jday                                  !KO
    real(r8) :: loss_manure_u                         !KO
    real(r8) :: loss_manure_n                         !KO
    real(r8) :: loss_manure_a                         !KO
    real(r8) :: loss_manure_r                         !KO
    real(r8) :: TF                                    !KO
    real(r8) :: ndep_fertilizer                       !KO
    real(r8) :: loss_fert_u                           !KO
    real(r8) :: R_nox_to_n2o                          !KO
    real(r8) :: kh                                    !KO
    real(r8) :: nh3_gas_conc                          !KO
    real(r8) :: fert_nh3_gas_conc                     !KO
    real(r8) :: knh4                                  !KO
    real(r8) :: nh3_aq_conc                           !KO
    real(r8) :: nh3_aq_sat                            !KO
    real(r8) :: NH4_fert                              !KO
    real(r8) :: NH4_manu                              !KO
    real(r8) :: kf                                    ! Decay rate of urea
    real(r8) :: dnh4,dno3                             ! Base rate diffusion for nh4 and no3
    real(r8) :: porosity                              !KO
    real(r8) :: ratenh4tosom,rateno3tosom             !KO
    real(r8) :: canopy_frac                           ! Fraction of NH3 emissions captured by canopy
    real(r8) :: manu_inc,fert_inc                     ! N mechanically incorporated into soil on a timescale of one year
    real(r8) :: k_relax                               ! Timescale for manure/fert water to relax to soil water (top 5cm)

    ! local pointers to implicit in scalars
    integer , pointer :: gridcell                (:)  ! Index into gridcell level quantities
    integer , pointer :: npfts                   (:)  ! Number of pfts for each column
    integer , pointer :: pfti                    (:)  ! Beginning pft index for each column
    real(r8), pointer :: latdeg                  (:)  ! Latitude in degrees

    real(r8), pointer :: forc_wind               (:)  ! Wind speed (m s-1)
    real(r8), pointer :: forc_rh                 (:)  ! Relative humidity(%)
    real(r8), pointer :: forc_ndep2              (:)  ! Nitrogen deposition rate (gN/ha/yr)
    real(r8), pointer :: forc_ndep3              (:)  ! Nitrogen deposition rate (gN/ha/yr)

    real(r8), pointer :: ndep_manure             (:)  ! Nitrogen manure deposition rate (gN/m2/s)
    real(r8), pointer :: ndep_fert               (:)  ! Nitrogen fertilizer deposition rate (gN/m2/s)
    real(r8), pointer :: nh3_manure              (:)  ! NH3 emission from manure(gN/m2/s)
    real(r8), pointer :: gamma_nh3               (:)  !KO
    real(r8), pointer :: gamma_nh3_fert          (:)  !KO
    real(r8), pointer :: nh3_fert                (:)  ! NH3 emission from fertilizer(gN/m2/s)
    real(r8), pointer :: nhxdep_to_sminn         (:)  ! NHx deposition from fertilizer & manure NH3 emissions(gN/m2/s) 
    real(r8), pointer :: noydep_to_sminn         (:)  ! NOy deposition (gN/m2/s) 
    real(r8), pointer :: nmanure_to_sminn        (:)  ! N deposition to soil mineral from manure(gN/m2/s)
    real(r8), pointer :: nfert_to_sminn          (:)  ! N deposition to soil mineral from fertilizer(gN/m2/s)
    real(r8), pointer :: N_Run_Off               (:)  ! Nitrogen Run Off from manure (gN/m2/s)
    real(r8), pointer :: N_Run_Off_fert          (:)  ! Nitrogen Run Off from fertilizer(gN/m2/s)
    real(r8), pointer :: manure_f_n2o_nit        (:)  ! N2O emission from nitrification of manure (gN/m2/s)
    real(r8), pointer :: manure_f_n2_denit       (:)  ! N2 emission from denitrification of manure (gN/m2/s)
    real(r8), pointer :: manure_f_nox_nit        (:)  ! NOx emission from nitrification of manure (gN/m2/s)
    real(r8), pointer :: fert_f_n2o_nit          (:)  ! N2O emission from nitrification of fertilizer (gN/m2/s)
    real(r8), pointer :: fert_f_n2_denit         (:)  ! N2 emission from denitrification of fertilizer (gN/m2/s)
    real(r8), pointer :: fert_f_nox_nit          (:)  ! NOx emission from nitrification of fertilizer (gN/m2/s)
    real(r8), pointer :: no3_manure_to_soil      (:)  ! NO3 flow from manure to soil @ 25 % NO3 pool per hour (gN/m2/s)
    real(r8), pointer :: TAN_manure_to_soil      (:)  ! NH4 flow from TAN manure to soil @ 1 % TAN pool per day (gN/m2/s)
    real(r8), pointer :: no3_fert_to_soil        (:)  ! NO3 flow from fertilizer to soil @ 25 % NO3 pool per hour (gN/m2/s)
    real(r8), pointer :: TAN_fert_to_soil        (:)  ! NH4 flow from TAN fertilizer to soil @ 1 % TAN pool per day (gN/m2/s)
    real(r8), pointer :: Nd                      (:)  ! Total Amount of N emitted during denitrification (gN/m2/s)
    real(r8), pointer :: lat_fert                (:)  !KO

    real(r8), pointer :: ndep_total              (:)  ! Total nitrogen deposition rate (gN/m2)
    real(r8), pointer :: TAN_manu                (:)  ! Manure Total Ammoniacal Nitrogen (gN/m2)
    real(r8), pointer :: TAN_fert                (:)  ! Fertilizer Total Ammoniacal Nitrogen (gN/m2)
    real(r8), pointer :: man_water_pool          (:)  ! Volume of water in manure/water solution (m3/m2)              
    real(r8), pointer :: fert_water_pool         (:)  ! Volume of water in fert/water solution (m3/m2)
    real(r8), pointer :: no3_manure              (:)  ! NO3 pool in manure (gN/m2)
    real(r8), pointer :: no3_fert                (:)  ! NO3 pool in fertilizer(gN/m2)
    real(r8), pointer :: nh3_manure_total        (:)  ! Total NH3 created from manure to atmosphere (gN/m2)
    real(r8), pointer :: no3_manure_total        (:)  ! Total NO3 created from manure(gN/m2)
    real(r8), pointer :: nh4_manure_total        (:)  ! Total NH4 created from manure (gN/m2)
    real(r8), pointer :: manure_u                (:)  ! Urine N pool in manure (gN/m2)
    real(r8), pointer :: manure_n                (:)  ! Non-minerizable N pool in manure (gN/m2)
    real(r8), pointer :: manure_a                (:)  ! Available N pool in manure (gN/m2)
    real(r8), pointer :: manure_r                (:)  ! Resistant N pool in manure (gN/m2)
    real(r8), pointer :: fert_u                  (:)  ! Fertilizer N pool (gN/m2)
    real(r8), pointer :: ndep_fert_total         (:)  ! Total nitrogen deposition rate of fertilizer(gN/m2)
    real(r8), pointer :: nh3_fert_total          (:)  ! Total NH3 created from fertilizer to atmosphere (gN/m2)
    real(r8), pointer :: no3_fert_total          (:)  ! Total NO3 created from fertilizer(gN/m2)
    real(r8), pointer :: nh4_fert_total          (:)  ! Total NH4 created from fertilizer (gN/m2)
    real(r8), pointer :: total_ndep              (:)  ! Total nitrogen deposition from Nr (gN/m2)
    real(r8), pointer :: total_nh3               (:)  ! Total NH3 created from Nr to atmosphere (gN/m2)
    real(r8), pointer :: total_N_Run_Off         (:)  ! Total N washed from Nr (gN/m2)
    real(r8), pointer :: total_no3               (:)  ! Total NO3 created from Nr(gN/m2)
    real(r8), pointer :: total_nh4               (:)  ! Total NH4 created from Nr (gN/m2)
    real(r8), pointer :: ra_col                  (:)  ! Aerodynamic resistance for grass pfts (s/m)
    real(r8), pointer :: rb_col                  (:)  ! Leaf boundary layer resistance for grass pfts (s/m)
    real(r8), pointer :: fert_app_jday           (:)  ! Julian day of the first fertilizer application (day)
    real(r8), pointer :: gdd8_col                (:)  ! Col growing degree-days base 8C from planting (ddays)
    real(r8), pointer :: t_a10_col               (:)  ! Col 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: t_a10min_col            (:)  ! Col 10-day running mean of min 2-m temperature (K)
    real(r8), pointer :: N_Run_Off_manure_total  (:)  ! Total N washed from manure (gN/m2)
    real(r8), pointer :: N_Run_Off_fert_total    (:)  ! Total N washed from fertilizer (gN/m2)

    real(r8), pointer :: leafn_manure            (:)  ! Leaf N (gN/m2)
    real(r8), pointer :: deadstemn_manure        (:)  ! Dead stem N eaten by cows (gN/m2)

    real(r8), pointer :: methane_manure          (:)  ! Emission of CH4 from cows (gC/m2/s)
    real(r8), pointer :: cmanure_to_sminn        (:)  ! Deposition of C from manure to soil mineral C (gC/m2/s)
    real(r8), pointer :: somhr                   (:)  ! Soil organic matter heterotrophic respiration (gC/m2/s)   

    real(r8), pointer :: total_leafc             (:)  ! Total C at column-level (gC/m2)
    real(r8), pointer :: leafc                   (:)  ! Leaf carbon (gC/m2)
    real(r8), pointer :: leafc_manure            (:)  ! Leaf C eaten by cows (gC/m2)
    real(r8), pointer :: deadstemc_manure        (:)  ! Dead stem C eaten by cows (gC/m2)

    real(r8), pointer :: t_grnd                  (:)  ! Ground temperature (K)
    real(r8), pointer :: gdd8_patch              (:)  ! patch Growing degree-days base 8C from planting (ddays)
    real(r8), pointer :: t_a10_patch             (:)  ! patch 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: t_a10min_patch          (:)  ! patch 10-day running mean of min 2-m temperature

    real(r8), pointer :: h2osoi_liqice_5cm       (:)  ! Liquid water + ice lens in top 5cm of soil (kg/m2)
    real(r8), pointer :: qflx_runoff             (:)  ! Total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)

    real(r8), pointer :: ram1                    (:)  ! Aerodynamical air resistance (s/m)
    real(r8), pointer :: rb1                     (:)  ! Aerodynamical boundary layer resistance (s/m)
!KO
    !-----------------------------------------------------------------------
    
    associate(                                                                & 
         forc_ndep     =>  atm2lnd_inst%forc_ndep_grc ,                       & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)
         ndep_to_sminn =>  soilbiogeochem_nitrogenflux_inst%ndep_to_sminn_col & ! Output: [real(r8) (:)]  atmospheric N deposition to soil mineral N (gN/m2/s)
         )
!KO
    if ( use_fan ) then

    gridcell           => col%gridcell
    npfts              => col%npatches
    pfti               => col%patchi
    latdeg             => grc%latdeg

    forc_wind          => atm2lnd_inst%forc_wind_grc
    forc_rh            => atm2lnd_inst%forc_rh_grc
    forc_ndep2         => atm2lnd_inst%forc_ndep2_grc
    forc_ndep3         => atm2lnd_inst%forc_ndep3_grc

    ndep_manure        => soilbiogeochem_nitrogenflux_inst%ndep_manure_col
    ndep_fert          => soilbiogeochem_nitrogenflux_inst%ndep_fert_col
    nh3_manure         => soilbiogeochem_nitrogenflux_inst%nh3_manure_col
    gamma_nh3          => soilbiogeochem_nitrogenflux_inst%gamma_nh3_col
    gamma_nh3_fert     => soilbiogeochem_nitrogenflux_inst%gamma_nh3_fert_col
    nh3_fert           => soilbiogeochem_nitrogenflux_inst%nh3_fert_col
    nhxdep_to_sminn    => soilbiogeochem_nitrogenflux_inst%nhxdep_to_sminn_col
    noydep_to_sminn    => soilbiogeochem_nitrogenflux_inst%noydep_to_sminn_col
    nmanure_to_sminn   => soilbiogeochem_nitrogenflux_inst%nmanure_to_sminn_col
    nfert_to_sminn     => soilbiogeochem_nitrogenflux_inst%nfert_to_sminn_col
    N_Run_Off          => soilbiogeochem_nitrogenflux_inst%N_Run_Off_col
    N_Run_Off_fert     => soilbiogeochem_nitrogenflux_inst%N_Run_Off_fert_col
    manure_f_n2o_nit   => soilbiogeochem_nitrogenflux_inst%manure_f_n2o_nit_col
    manure_f_n2_denit  => soilbiogeochem_nitrogenflux_inst%manure_f_n2_denit_col
    manure_f_nox_nit   => soilbiogeochem_nitrogenflux_inst%manure_f_nox_nit_col
    fert_f_n2o_nit     => soilbiogeochem_nitrogenflux_inst%fert_f_n2o_nit_col
    fert_f_n2_denit    => soilbiogeochem_nitrogenflux_inst%fert_f_n2_denit_col
    fert_f_nox_nit     => soilbiogeochem_nitrogenflux_inst%fert_f_nox_nit_col
    no3_manure_to_soil => soilbiogeochem_nitrogenflux_inst%no3_manure_to_soil_col
    TAN_manure_to_soil => soilbiogeochem_nitrogenflux_inst%TAN_manure_to_soil_col
    no3_fert_to_soil   => soilbiogeochem_nitrogenflux_inst%no3_fert_to_soil_col
    TAN_fert_to_soil   => soilbiogeochem_nitrogenflux_inst%TAN_fert_to_soil_col
    Nd                 => soilbiogeochem_nitrogenflux_inst%Nd_col
    lat_fert           => soilbiogeochem_nitrogenflux_inst%lat_fert_col

    ndep_total         => soilbiogeochem_nitrogenstate_inst%ndep_total_col
    TAN_manu           => soilbiogeochem_nitrogenstate_inst%TAN_manu_col
    TAN_fert           => soilbiogeochem_nitrogenstate_inst%TAN_fert_col
    man_water_pool     => soilbiogeochem_nitrogenstate_inst%man_water_pool_col
    fert_water_pool    => soilbiogeochem_nitrogenstate_inst%fert_water_pool_col
    no3_manure         => soilbiogeochem_nitrogenstate_inst%no3_manure_col
    no3_fert           => soilbiogeochem_nitrogenstate_inst%no3_fert_col
    nh3_manure_total   => soilbiogeochem_nitrogenstate_inst%nh3_manure_total_col
    no3_manure_total   => soilbiogeochem_nitrogenstate_inst%no3_manure_total_col
    nh4_manure_total   => soilbiogeochem_nitrogenstate_inst%nh4_manure_total_col
    manure_u           => soilbiogeochem_nitrogenstate_inst%manure_u_col
    manure_n           => soilbiogeochem_nitrogenstate_inst%manure_n_col
    manure_a           => soilbiogeochem_nitrogenstate_inst%manure_a_col
    manure_r           => soilbiogeochem_nitrogenstate_inst%manure_r_col
    fert_u             => soilbiogeochem_nitrogenstate_inst%fert_u_col
    ndep_fert_total    => soilbiogeochem_nitrogenstate_inst%ndep_fert_total_col
    nh3_fert_total     => soilbiogeochem_nitrogenstate_inst%nh3_fert_total_col
    no3_fert_total     => soilbiogeochem_nitrogenstate_inst%no3_fert_total_col
    nh4_fert_total     => soilbiogeochem_nitrogenstate_inst%nh4_fert_total_col
    total_ndep         => soilbiogeochem_nitrogenstate_inst%total_ndep_col
    total_nh3          => soilbiogeochem_nitrogenstate_inst%total_nh3_col
    total_N_Run_Off    => soilbiogeochem_nitrogenstate_inst%total_N_Run_Off_col
    total_no3          => soilbiogeochem_nitrogenstate_inst%total_no3_col
    total_nh4          => soilbiogeochem_nitrogenstate_inst%total_nh4_col
    ra_col             => soilbiogeochem_nitrogenstate_inst%ra_col
    rb_col             => soilbiogeochem_nitrogenstate_inst%rb_col
    fert_app_jday      => soilbiogeochem_nitrogenstate_inst%fert_app_jday_col
    gdd8_col           => soilbiogeochem_nitrogenstate_inst%gdd8_col
    t_a10_col          => soilbiogeochem_nitrogenstate_inst%t_a10_col
    t_a10min_col       => soilbiogeochem_nitrogenstate_inst%t_a10min_col
    N_Run_Off_manure_total => soilbiogeochem_nitrogenstate_inst%N_Run_Off_manure_total_col
    N_Run_Off_fert_total   => soilbiogeochem_nitrogenstate_inst%N_Run_Off_fert_total_col

    leafn_manure       => cnveg_nitrogenstate_inst%leafn_manure_patch
    deadstemn_manure   => cnveg_nitrogenstate_inst%deadstemn_manure_patch

    methane_manure     => soilbiogeochem_carbonflux_inst%methane_manure_col
    cmanure_to_sminn   => soilbiogeochem_carbonflux_inst%cmanure_to_sminn_col
    somhr              => soilbiogeochem_carbonflux_inst%somhr_col

    total_leafc        => cnveg_carbonstate_inst%total_leafc_col
    leafc              => cnveg_carbonstate_inst%leafc_patch
    leafc_manure       => cnveg_carbonstate_inst%leafc_manure_patch
    deadstemc_manure   => cnveg_carbonstate_inst%deadstemc_manure_patch

    t_grnd             => temperature_inst%t_grnd_col
    gdd8_patch         => temperature_inst%gdd8_patch
    t_a10_patch        => temperature_inst%t_a10_patch
    t_a10min_patch     => temperature_inst%t_a10min_patch

    h2osoi_liqice_5cm  => waterstate_inst%h2osoi_liqice_5cm_col

    qflx_runoff        => waterflux_inst%qflx_runoff_col

    ram1               => frictionvel_inst%ram1_patch
    rb1                => frictionvel_inst%rb1_patch

    end if
!KO

!KO
      ! Loop through columns
      do c = bounds%begc, bounds%endc
         g = col%gridcell(c)
         ndep_to_sminn(c) = forc_ndep(g)
      end do
      
      if ( use_fan ) then

      fmm = 0._r8
      na  = 0._r8
      ca  = 0._r8
      pl  = 1.0_r8
      pd  = 0._r8
      cn  = 30._r8          !Carbon:nitrogen
      fua = 0.2_r8
      fu  = 0.5_r8
      ffa = 0.27_r8
      ff  = 0.5_r8

      !*pgmh urea decomposition rate (1/s)
      kf          = 4.83e-06_r8      ! [Agehara and Warncke, 2005
      dnh4        = 9.8e-10_r8       ! from genermont and cellier
      dno3        = 1.3e-08_r8       ! fast x13 than dnh4 consistent with data
      porosity    = 0.5_r8
      !*pgmh
      canopy_frac = 0.6_r8

      call get_curr_date(yr, mon, day, sec)
      m1850 = 0.00000000004_r8 * exp(0.012_r8 * yr)
      f1850 = 1.05_r8/(1._r8+exp(0.13_r8*(1975-yr)))

      dt = real( get_step_size(), r8 )
!KO
      
      ! Loop through columns
      do c = bounds%begc, bounds%endc
         g = col%gridcell(c)
!KO         ndep_to_sminn(c) = forc_ndep(g)
!KO
         !calculate the nitrogen excreted (kg ha-1 s-1) by cows from file
         !Nmanure.nc, then change to (g m-2 s-1)

         ndep_manure(c) = forc_ndep2(g) / (10._r8)  * m1850
         ndep_fert(c) = forc_ndep3(g) / (10._r8)  * f1850

         lat_fert(c) = latdeg(g)

         !Remove any NaN values
         if (ndep_manure(c) .ne. ndep_manure(c)) ndep_manure(c) = 0._r8
         if (ndep_fert(c) .ne. ndep_fert(c)) ndep_fert(c) = 0._r8
!KO

      end do
!KO
      ! Convert pft-level variables to column-level for computing fertilizer application date and 
      ! aerodynamic resistance factors.  These conversions will not be necessary once this code is
      ! implemented in the crop model where crop calculations are done on their own column.  
      call p2c(bounds, num_soilc, filter_soilc, &
           ram1(bounds%begp:bounds%endp), &
           ra_col(bounds%begc:bounds%endc))
      call p2c(bounds, num_soilc, filter_soilc, &
           rb1(bounds%begp:bounds%endp), &
           rb_col(bounds%begc:bounds%endc))
      call p2c(bounds, num_soilc, filter_soilc, &
           gdd8_patch(bounds%begp:bounds%endp), &
           gdd8_col(bounds%begc:bounds%endc))
      call p2c(bounds, num_soilc, filter_soilc, &
           t_a10_patch(bounds%begp:bounds%endp), &
           t_a10_col(bounds%begc:bounds%endc))
      call p2c(bounds, num_soilc, filter_soilc, &
           t_a10min_patch(bounds%begp:bounds%endp), &
           t_a10min_col(bounds%begc:bounds%endc))

      do fc = 1,num_soilc
         c = filter_soilc(fc)
                
         do pi = 1,max_patch_per_col
            if ( pi <=  npfts(c) ) then
               p = pfti(c) + pi - 1
!KO
               if (patch%active(p)) then
!KO
!KO                  total_leafc(c) = total_leafc(c) + leafc(p)
!KO
                  total_leafc(c) = total_leafc(c) + leafc(p) * patch%wtcol(p)
!KO
               end if
            end if
         end do
                        
         ! Fertilize crops 1 time a year when planting conditions are first met.  The conditions
         ! for corn planting are used here for all crops.  The earliest date of planting is April
         ! 1st and the latest date is June 14th for the NHemis, add 180 days for Shemis planting
         jday    = get_curr_calday()

         if (jday .eq. 1) fert_app_jday(c) = 0._r8

         if (t_a10_col(c) >= 283.15_r8 .and. t_a10min_col(c) >= 279.15_r8 .and. &
             gdd8_col(c) >= 50._r8 .and. fert_app_jday(c) .eq. 0._r8 .and. &
             jday .gt. 90._r8 .and. lat_fert(c) .gt. 0._r8) fert_app_jday(c) = jday
         if (t_a10_col(c) >= 283.15_r8 .and. t_a10min_col(c) >= 279.15_r8 .and. &
             gdd8_col(c) >= 50._r8 .and. fert_app_jday(c) .eq. 0._r8 .and. &
             jday .gt. 272._r8 .and. lat_fert(c) .lt. 0._r8) fert_app_jday(c) = jday

         if (jday .eq. 165._r8 .and. lat_fert(c) .gt. 0._r8 .and. fert_app_jday(c) .eq. 0._r8) fert_app_jday(c) = jday
         if (jday .eq. 347._r8 .and. lat_fert(c) .lt. 0._r8 .and. fert_app_jday(c) .eq. 0._r8) fert_app_jday(c) = jday

         ! Add fertilizer (ndep_fertilizer)
         if (jday .eq. fert_app_jday(c)) then
            ndep_fertilizer = (ndep_fert(c) * 3600._r8 * 24._r8 * 365._r8 )

         ! If there is any N in existing fertilizer pools, add this to the soil pool (later on in code)
         ! before adding the new fertilizer N to the TAN_fert pool.
            fert_inc = TAN_fert(c)+fert_u(c)
            fert_u(c) = 0._r8
            TAN_fert(c) = 0._r8
         else
            !KO fert_inc needs to be set here, for now assume 0?
            fert_inc = 0._r8
            ndep_fertilizer = 0._r8
         end if

         fert_u(c) = fert_u(c) + ndep_fertilizer
                
         ! Fertilizer decay into the TAN pool following King and Balogh, 2000 (using one day timescale)
         if (t_grnd(c) <= 273._r8) then
            loss_fert_u = 0._r8
         else
!*pgmh changed to a rate consistent with urea
!           loss_fert_u = fert_u(c)*(1._r8 - exp(-1._r8 * dt * (0.04_r8/(86400._r8)) *(((t_grnd(c) - 273._r8) / 38._r8))**1.33_r8))
            fert_u(c) = fert_u(c) / (1._r8+dt*kf)
            loss_fert_u = dt * kf * fert_u(c)
            TAN_fert(c) = TAN_fert(c) + loss_fert_u
         end if

!        TAN_fert(c) = TAN_fert(c) + loss_fert_u
!        fert_u(c) = fert_u(c) - loss_fert_u
!*pgmh   
         ! Manure decay into TAN pool

         ! Here the N from manure is divided into separate pools.
         ! Urine (urea) is 50 % of the N excreted by the animals (Parton et al., 1987)
         ! Non-mineralizable = 10 % feces N, Available = 44 % feces N & Resistant = 46 % feces N (Andrews, 1996)

         manure_u(c) = manure_u(c) + (ndep_manure(c)* 0.5_r8 *dt)
         manure_n(c) = manure_n(c) + (0.025_r8 * ndep_manure(c) * dt)
         manure_a(c) = manure_a(c) + (0.25_r8 * ndep_manure(c) * dt)
         manure_r(c) = manure_r(c) + (0.225_r8 * ndep_manure(c) * dt)
        
         TF = 0.0106_r8 * exp(0.12979_r8 * (t_grnd(c) - 273._r8))
        
         loss_manure_u = manure_u(c)
         loss_manure_n = 0._r8
         loss_manure_a = manure_a(c) * (1._r8 - exp(TF * (8.938e-7_r8 * (-1._r8) * dt)))
         loss_manure_r = manure_r(c) * (1._r8 - exp(TF * (6.38e-8_r8 * (-1._r8) * dt)))

         TAN_manu(c) = TAN_manu(c) + loss_manure_u + loss_manure_n + loss_manure_a + loss_manure_r
        
         ! Represent mechanical incorporation of manure into soil, timescale of one year
         manu_inc = (dt/(365._r8 * 86400._r8)) * (manure_r(c) + manure_a(c))

         manure_u(c) = manure_u(c) - loss_manure_u
         manure_n(c) = manure_n(c) - loss_manure_n
         manure_r(c) = manure_r(c) - loss_manure_r - (dt/(365._r8 * 86400._r8))*manure_r(c)
         manure_a(c) = manure_a(c) - loss_manure_a - (dt/(365._r8 * 86400._r8))*manure_a(c)

         ! Take N content of manure as 1% dry matter and the water content of manure as 85% of total mass
         ! Thus, the water content added is 567 times the N content added at each timstep.
         k_relax = dt/(3._r8*86400._r8)

         ! Add manure water and relax the manure water content to the water content of the top 5 cm of soil 
         ! on a timescale of 3 days (k_relax)
         if (TAN_manu(c) .gt. 0._r8) then
            man_water_pool(c) = man_water_pool(c)+(1.e-6_r8*ndep_manure(c)*dt*566.6667_r8) - &
            k_relax*(man_water_pool(c)-h2osoi_liqice_5cm(c)*1.e-3_r8)
         else
            man_water_pool(c) = 0._r8
            TAN_manu(c) = 0._r8
         end if

         if (man_water_pool(c) .ne. man_water_pool(c)) man_water_pool(c) = 0._r8
         !Remove Manure TAN pool washed off by rain using the method of the global NEWS model (Harrison et al., 2005)
         if (TAN_manu(c) <= 0._r8) then
            N_Run_Off(c) = 0._r8
            TAN_manu(c) = 0._r8
!*pgmh convert qflx_surf from mm/sec to m/sec
!        elseif ((0.94_r8 * TAN_manu(c) * (qflx_surf(c)))*dt > TAN_manu(c)) then
         elseif ((TAN_manu(c)/man_water_pool(c) * (qflx_runoff(c)*1.e-3_r8))*dt > TAN_manu(c)) then
!*pgmh
            N_Run_Off(c) = TAN_manu(c) / dt
         else
!*pgmh convert qflx_surf from mm/sec to m/sec
!           N_Run_Off(c) = 0.94_r8 * TAN_manu(c) * (qflx_surf(c))
            N_Run_Off(c) = TAN_manu(c)/man_water_pool(c) * (qflx_runoff(c)*1.e-3_r8)
!*pgmh
         end if
            
!KO         TAN_manu(c) = TAN_manu(c) - N_Run_Off(c)*dt    
!KO
         ! Protect against negative values of TAN_manu which can cause divide by
         ! zero errors in NFactor below, and adjust N_Run_Off accordingly
         TAN_manu(c)  = TAN_manu(c) - min(N_Run_Off(c)*dt,TAN_manu(c))
         N_Run_Off(c) = min(N_Run_Off(c),TAN_manu(c)/dt)
!KO

         ! Calculate nh3_gas_conc. (PGMH 29th August 2014).  Note that a neutral pH is assumed (1.e-7) and that 
         ! the units of nh3_gas_conc are g/m3.  If the column is completely dry, the nh3_gas_conc is set to zero.
         kh = 56._r8*exp(4092._r8*((1._r8/t_grnd(c))-(1._r8/298.15_r8)))
         knh4 = 5.67e-10_r8*exp(-6286._r8*((1._r8/t_grnd(c))-(1._r8/298.15_r8))) ! [mole/Liter]

         if (man_water_pool(c) .le. 0._r8) then
            man_water_pool(c) = 0._r8
            nh3_gas_conc = 0._r8
         else
            nh3_gas_conc = (TAN_manu(c)/man_water_pool(c))/ &
                           (1._r8+(kh*t_grnd(c)/12.2_r8)+(kh*t_grnd(c)/12.2_r8)*(1.e-7_r8/knh4))
         endif

         ! Compute nh3 saturation with respect to the aqueous concentration of nh3
         nh3_aq_conc = (t_grnd(c)/12.2_r8)*kh*nh3_gas_conc   ! [g/m3]
         nh3_aq_sat = (10._r8**((966.7475_r8/t_grnd(c))-0.57953_r8))*1.e3_r8    ! [g/m3]

         if (nh3_aq_conc .gt. nh3_aq_sat) then
            nh3_aq_conc=nh3_aq_sat
            nh3_gas_conc=nh3_aq_sat*(12.2_r8/t_grnd(c))*(1._r8/kh)
         endif

         ! Calculate NH3 emission
         if (nh3_gas_conc < 3.e-7_r8) then
            nh3_manure(c) = 0._r8
         else
!*pgmh
!            nh3_manure(c) = ((nh3_gas_conc-3.e-7)/(ra_col(c)+rb_col(c))) * 0.35_r8
            nh3_manure(c) = ((nh3_gas_conc-3.e-7_r8)/(ra_col(c)+rb_col(c))) * (1._r8-canopy_frac)
!*pgmh
         end if

         if (nh3_manure(c) .ne. nh3_manure(c)) nh3_manure(c) = 0._r8
        
         !Calculate the amount of water in TAN of Fertilizer.  Water is added when fertilizer
         ! is applied, and then the total water pool is relaxed to the water in the top 5cm of soil
         if (TAN_fert(c) .gt. 0._r8) then
            fert_water_pool(c) = fert_water_pool(c)+(((1.e-6_r8*(ndep_fertilizer*dt))/0.466_r8)/ &
                                 (0.66_r8 *exp(0.0239_r8*(t_grnd(c)-273._r8))))-&
                                 k_relax*(fert_water_pool(c)-h2osoi_liqice_5cm(c)*1.e-3_r8)
         else
            fert_water_pool(c) = 0._r8
            TAN_fert(c) = 0._r8
         end if
        
         if (fert_water_pool(c) .ne. fert_water_pool(c)) fert_water_pool(c) = 0._r8

         !Remove Fertilizer TAN pool washed off by rain at rate of the global NEWS model (Harrison et al., 2005)
         if (TAN_fert(c) <= 0._r8) then
            N_Run_Off_fert(c) = 0._r8
            TAN_fert(c) = 0._r8
!*pgmh convert qflx_surf from mm/sec to m/sec
!        elseif ((0.94_r8 * TAN_fert(c) *(qflx_surf(c)))*dt > TAN_fert(c)) then
         elseif ((TAN_fert(c)/fert_water_pool(c) *(qflx_runoff(c)*1.e-3_r8))*dt > TAN_fert(c)) then
!*pgmh
            N_Run_Off_fert(c) = TAN_fert(c) / dt
         else
!*pgmh convert qflx_surf from mm/sec to m/sec
!           N_Run_Off_fert(c) = 0.94_r8 * TAN_fert(c) *(qflx_surf(c))
            N_Run_Off_fert(c) = TAN_fert(c)/fert_water_pool(c) *(qflx_runoff(c)*1.e-3_r8)
!*pgmh
         end if
!KO         TAN_fert(c) = TAN_fert(c) - N_Run_Off_fert(c)*dt
!KO
         ! Protect against negative values of TAN_fert which can cause divide by
         ! zero errors in NFactor below, and adjust N_Run_Off_fert accordingly
         TAN_fert(c)  = TAN_fert(c) - min(N_Run_Off_fert(c)*dt,TAN_fert(c))
         N_Run_Off_fert(c) = min(N_Run_Off_fert(c),TAN_fert(c)/dt)
!KO

         !Calculate gamma for fertilizer
         if (fert_water_pool(c) .le. 0._r8) then
            fert_water_pool(c) = 0._r8
            fert_nh3_gas_conc = 0._r8
         else
            fert_nh3_gas_conc = (TAN_fert(c)/fert_water_pool(c))/ &
                                (1._r8+(kh*t_grnd(c)/12.2_r8)+(kh*t_grnd(c)/12.2_r8)*(1.e-7_r8/knh4))
         endif

         ! Compute nh3 saturation with respect to the aqueous concentration of nh3
         nh3_aq_conc = (t_grnd(c)/12.2_r8)*kh*fert_nh3_gas_conc   ! [g/m3]

         if (nh3_aq_conc .gt. nh3_aq_sat) then
            nh3_aq_conc=nh3_aq_sat
            fert_nh3_gas_conc=nh3_aq_sat*(12.2_r8/t_grnd(c))*(1._r8/kh)
         endif

         ! Calculate NH3 emission
         if (fert_nh3_gas_conc < 3.e-7_r8) then
            nh3_fert(c) = 0._r8
         else
!*pgmh
!           nh3_fert(c) = ((fert_nh3_gas_conc-3.e-7)/(ra_col(c)+rb_col(c))) * 0.35_r8
            nh3_fert(c) = ((fert_nh3_gas_conc-3.e-7_r8)/(ra_col(c)+rb_col(c))) * (1._r8-canopy_frac)
!*pgmh
         end if
         if (nh3_fert(c) .ne. nh3_fert(c)) nh3_fert(c)=0._r8

         ! Calculate transfert to soil organic matter
         NH4_manu = (kh*t_grnd(c)/12.2_r8)*(1.e-7_r8/knh4)*nh3_gas_conc*man_water_pool(c)
         if (t_grnd(c) > 313._r8) then
            !nmanure_to_sminn(c) = 0._r8
!*dsw set to canopy captured N for very high temperatures
            nmanure_to_sminn(c) = (nh3_manure(c) / (1._r8-canopy_frac) * canopy_frac)
         else
!*pgmh units incorect for conversion, also assume 70% leave capture
!           nmanure_to_sminn(c) = ((2.*1.16e-6*NH4_manu) / ((1./(1.0 - exp(-1.0 * ((h2osoi_liqice_5cm(c)*9.5e-3)/ 0.12)**2.0))) &
!                          +(1./((((40. - (t_grnd(c)-273.))/12.)**2.4) * exp(0.2*((t_grnd(c)-273.)-28.)))))) &
!                                 +(nh3_manure(c) / 0.35_r8 * 0.65_r8)
            nmanure_to_sminn(c) = ((2._r8*1.16e-6_r8*NH4_manu) / &
                                  ((1._r8/(1._r8 - exp(-1._r8 * ((h2osoi_liqice_5cm(c)*0.019_r8)/ 0.12_r8)**2._r8))) + &
                                  (1._r8/((((40._r8 - (t_grnd(c)-273._r8))/12._r8)**2.4_r8) * &
                                  exp(0.2_r8*((t_grnd(c)-273._r8)-28._r8)))))) &
                                  +(nh3_manure(c) / (1._r8-canopy_frac) * canopy_frac)
!*pgmh units incorect for conversion
         end if
        
         ! Calculate soil organic matter
         NH4_fert = (kh*t_grnd(c)/12.2_r8)*(1.e-7_r8/knh4)*fert_nh3_gas_conc*fert_water_pool(c)
         if (t_grnd(c) > 313._r8) then
!           nfert_to_sminn(c) = 0._r8  
!*dsw set to canopy captured N for very high temperatures
            nfert_to_sminn(c) = (nh3_fert(c) / (1._r8-canopy_frac) * canopy_frac)
         else
!*pgmh units incorect for conversion
!           nfert_to_sminn(c) = ((2.*1.16e-6*NH4_fert) / ((1./(1.0 - exp(-1.0 * ((h2osoi_liqice_5cm(c)*9.5e-3)/0.12)**2.0))) &
!                        +(1./((((40. - (t_grnd(c)-273.))/12.)**2.4) * exp(0.2*((t_grnd(c)-273.)-28.)))))) &
!                               +(nh3_fert(c) / 0.35_r8 * 0.65_r8)

            nfert_to_sminn(c) = ((2._r8*1.16e-6_r8*NH4_fert) / &
                                ((1._r8/(1._r8 - exp(-1._r8 * ((h2osoi_liqice_5cm(c)*0.019_r8)/0.12_r8)**2._r8))) + &
                                (1._r8/((((40._r8 - (t_grnd(c)-273._r8))/12._r8)**2.4_r8) * &
                                exp(0.2_r8*((t_grnd(c)-273._r8)-28._r8)))))) &
                                +(nh3_fert(c) / (1._r8-canopy_frac) * canopy_frac)
!*pgmh units incorect for conversion
         end if

         !Calculate N2O & NOx fluxes from nitrification of manure and fertilizer

         manure_f_n2o_nit(c) = nmanure_to_sminn(c) * 0.0_r8   !0.02_r8 
         fert_f_n2o_nit(c) = nfert_to_sminn(c) * 0.0_r8       !0.02_r8
        
         ! NOx FLUXES
         !------------

         ! Ratio of NOx to N2O from nit or denit: eq 5, Parton '01.  Add soi_gd later for now soi_gd = 0.1
         R_nox_to_n2o = 15.2_r8+(35.5_r8*ATAN(0.68_r8*3.14_r8*(10._r8*0.4_r8 -1.86_r8)))/3.14_r8
        
         manure_f_nox_nit(c) = manure_f_n2o_nit(c) * R_nox_to_n2o
         fert_f_nox_nit(c) = fert_f_n2o_nit(c) * R_nox_to_n2o
        
         nh3_manure(c) = abs(nh3_manure(c))
         nmanure_to_sminn(c) = abs(nmanure_to_sminn(c))
         manure_f_n2o_nit(c) = abs(manure_f_n2o_nit(c))
         manure_f_nox_nit(c) = abs(manure_f_nox_nit(c))
         N_Run_Off(c) = abs(N_Run_Off(c))
         nh3_fert(c) = abs(nh3_fert(c))
         nfert_to_sminn(c) = abs(nfert_to_sminn(c))
         fert_f_n2o_nit(c) = abs(fert_f_n2o_nit(c))
         fert_f_nox_nit(c) = abs(fert_f_nox_nit(c))
         N_Run_Off_fert(c) = abs(N_Run_Off_fert(c))
        
         if (TAN_manu(c) <((nh3_manure(c)+ nmanure_to_sminn(c) + manure_f_n2o_nit(c) +  manure_f_nox_nit(c) )*dt)) then
             NFactor = (nh3_manure(c) + nmanure_to_sminn(c)+ manure_f_n2o_nit(c) +  manure_f_nox_nit(c) )
             nh3_manure(c) = nh3_manure(c) * (TAN_manu(c)/dt) / NFactor
             nmanure_to_sminn(c) = nmanure_to_sminn(c) * (TAN_manu(c)/dt) / NFactor
             manure_f_n2o_nit(c) = manure_f_n2o_nit(c) * (TAN_manu(c)/dt) / NFactor
             manure_f_nox_nit(c) = manure_f_nox_nit(c) * (TAN_manu(c)/dt) / NFactor

             TAN_manu(c) = 0._r8
         else
             TAN_manu(c) = TAN_manu(c) - ((nh3_manure(c) + nmanure_to_sminn(c) + manure_f_n2o_nit(c) +  manure_f_nox_nit(c) )*dt)
         end if
        
         if (TAN_fert(c) <((nh3_fert(c)+ nfert_to_sminn(c) + fert_f_n2o_nit(c) +  fert_f_nox_nit(c) )*dt)) then
             NFactor = (nh3_fert(c) + nfert_to_sminn(c) + fert_f_n2o_nit(c) +  fert_f_nox_nit(c) )
             nh3_fert(c) = nh3_fert(c) * (TAN_fert(c)/dt) / NFactor
             nfert_to_sminn(c) = nfert_to_sminn(c) * (TAN_fert(c)/dt) / NFactor
             fert_f_n2o_nit(c) = fert_f_n2o_nit(c) * (TAN_fert(c)/dt) / NFactor
             fert_f_nox_nit(c) = fert_f_nox_nit(c) * (TAN_fert(c)/dt) / NFactor
        
             TAN_fert(c) = 0._r8
         else
             TAN_fert(c) = TAN_fert(c) - ((nh3_fert(c) + nfert_to_sminn(c) + fert_f_n2o_nit(c) +  fert_f_nox_nit(c) )*dt)
         end if
        
         !Add N from mechanical incorporation into soil N pools
         nmanure_to_sminn(c)=nmanure_to_sminn(c)+manu_inc/dt
         nfert_to_sminn(c)=nfert_to_sminn(c)+fert_inc/dt

!*pgmh!!!!! diffusion times for lengths of 1 cm (dz^2=1.e-4 m)
!         ratenh4tosom=1.e-4*dnh4*(1.03**(t_grnd(c)-273.))*((man_water_pool(c)/5.e-2)**(10./3.))/(porosity**2.)
!         rateno3tosom=1.e-4*dno3*(1.03**(t_grnd(c)-273.))*((fert_water_pool(c)/5.e-2)**(10./3.))/(porosity**2.)
         ratenh4tosom = 1.e4_r8*dnh4*(1.03_r8**(t_grnd(c)-273._r8)) * &
                        ((man_water_pool(c)/5.e-2_r8)**(10._r8/3._r8))/(porosity**2._r8)
         rateno3tosom = 1.e4_r8*dno3*(1.03_r8**(t_grnd(c)-273._r8)) * &
                        ((fert_water_pool(c)/5.e-2_r8)**(10._r8/3._r8))/(porosity**2._r8)
         if (t_grnd(c) < 273.15_r8) then ! Turn vertical diffusion off when temperature is below freezing
           ratenh4tosom = 0._r8
           rateno3tosom = 0._r8
         endif

!*pgmh!!!!!
         !Calculate rate that SOM forms NH4 to the soil 
!*pgmh
!        TAN_manure_to_soil(c) =  1.16e-7 * TAN_manu(c) 
         TAN_manure_to_soil(c) =  ratenh4tosom * TAN_manu(c)
!*pgmh
         TAN_manu(c) = TAN_manu(c) - (TAN_manure_to_soil(c) * dt)
        
!*pgmh
!        TAN_fert_to_soil(c) =  1.16e-7 * TAN_fert(c) 
         TAN_fert_to_soil(c) =  ratenh4tosom * TAN_fert(c)
!*pgmh
         TAN_fert(c) = TAN_fert(c) - (TAN_fert_to_soil(c) * dt)
        
         !Calculate rate that NO3- goes to soil
!* pgmh: at present all no3 goes to the soil in a timestep. nmanure_to_sminn goes back to the no3_pool.
!This includes nitrification and mechanical incorporation of manure (not correct!!)
!*pgmh
!        no3_manure_to_soil(c) = (no3_manure(c) / dt)
!        no3_manure(c) = (nmanure_to_sminn(c) * dt)
!        
!        no3_fert_to_soil(c) = (no3_fert(c) / dt)
!        no3_fert(c) = (nfert_to_sminn(c) * dt)
!pgmh: suggested changes: use diffusion coefficient, correct for no3 pool and 
!         no3_manure_to_soil(c)=rateno3tosom * no3_manure(c)
!         no3_manure(c)=no3_manure(c)+nmanure_to_sminn(c)*dt - (no3_manure_to_soil(c) * dt)
!         no3_fert_to_soil(c)=rateno3tosom * no3_fert(c)
!         no3_fert(c)=no3_fert(c) +nfert_to_sminn(c)*dt - (no3_fert_to_soil(c) * dt)
!*pgmh
         no3_manure_to_soil(c) = rateno3tosom * no3_manure(c)
         no3_fert_to_soil(c) = rateno3tosom * no3_fert(c)
         no3_manure(c) = no3_manure(c) + nmanure_to_sminn(c)*dt
         no3_fert(c) = no3_fert(c) + nfert_to_sminn(c)*dt
         if (no3_manure_to_soil(c)*dt .gt. no3_manure(c)) no3_manure_to_soil(c) = no3_manure(c)/dt
         if (no3_fert_to_soil(c)*dt .gt. no3_fert(c)) no3_fert_to_soil(c) = no3_fert(c)/dt
         no3_fert(c) = no3_fert(c) - (no3_fert_to_soil(c) * dt)
         no3_manure(c) = no3_manure(c) - (no3_manure_to_soil(c) * dt)

         ! Calculate NHx & NOy deposition as 45 & 55 % ndep_to_sminn, respectively  
         ! Then NHx and NOy as a fraction of average Ndep from Peter Hess
        
         call get_curr_date(yr, mon, day, sec)
         if (mon .eq. 1) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 0.986_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 1.185_r8
         elseif (mon .eq. 2) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 0.894_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 1.496_r8
         elseif (mon .eq. 3) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 0.825_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 1.152_r8
         elseif (mon .eq. 4) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 0.792_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 1.013_r8
         elseif (mon .eq. 5) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 0.814_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 0.739_r8
         elseif (mon .eq. 6) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 1.201_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 0.833_r8
         elseif (mon .eq. 7) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 1.242_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 1.038_r8
         elseif (mon .eq. 8) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 1.141_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 0.718_r8
         elseif (mon .eq. 9) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 0.986_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 0.856_r8
         elseif (mon .eq. 10) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 0.964_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 0.808_r8
         elseif (mon .eq. 11) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 1.155_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 1.024_r8
         elseif (mon .eq. 12) then
           noydep_to_sminn(c) = 0.55_r8 * ndep_to_sminn(c) * 1.021_r8
           nhxdep_to_sminn(c) = 0.45_r8 * ndep_to_sminn(c) * 1.137_r8
         endif  
                
         !calculate the weighted amount of C and N eaten at the PFT level
         do pi = 1,max_patch_per_col
            if ( pi <=  npfts(c) ) then
!KO
               p = pfti(c) + pi - 1
               if (patch%active(p)) then
!KO
                if (total_leafc(c) > 0._r8) then
                   leafc_manure(p) = leafc(p) * (cn * ((pl * ndep_manure(c))/(1._r8 - na))) * dt / total_leafc(c)
                   leafn_manure(p) = leafc(p) * ((pl * ndep_manure(c))/(1._r8 - na)) * dt / total_leafc(c)
                   deadstemc_manure(p) = leafc(p) * (cn * ((pd * ndep_manure(c))/(1._r8 - na))) * dt / total_leafc(c)
                   deadstemn_manure(p) = leafc(p) * ((pd * ndep_manure(c))/(1._r8 - na)) * dt / total_leafc(c)
                else
                   leafc_manure(p) = 0._r8
                   leafn_manure(p) = 0._r8
                   deadstemc_manure(p) = 0._r8
                   deadstemn_manure(p) = 0._r8
                end if
               end if
            end if
         end do

      end do  ! End column loop

      end if  ! End use_fan
!KO

    end associate

  end subroutine CNNDeposition_Old

  !-----------------------------------------------------------------------
  subroutine CNFreeLivingFixation(num_soilc, filter_soilc, &
       waterflux_inst, soilbiogeochem_nitrogenflux_inst)


    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
 
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter                                                                                                                     
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns                                                                                                                                  
   
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    type(waterflux_type)                   , intent(inout) :: waterflux_inst 
    !
    ! !LOCAL VARIABLES:                                                                                                                                                                                                           
    integer  :: c,fc            !indices     
    real(r8) :: dayspyr         !days per year 
    real(r8) :: secs_per_year   !seconds per year   

       associate(                                                                        &
                  AnnET            => waterflux_inst%AnnET,                              & ! Input:  [real(:)  ] : Annual average ET flux mmH20/s
                  freelivfix_slope => params_inst%freelivfix_slope_wET,                  & ! Input:  [real     ] : slope of fixation with ET
                  freelivfix_inter => params_inst%freelivfix_intercept,                  & ! Input:  [real     ] : intercept of fixation with ET
                  ffix_to_sminn    => soilbiogeochem_nitrogenflux_inst%ffix_to_sminn_col & ! Output: [real(:)  ] : free living N fixation to soil mineral N (gN/m2/s)
                ) 
       
       dayspyr = get_days_per_year()
       secs_per_year = dayspyr*24_r8*3600_r8

       do fc = 1,num_soilc
           c = filter_soilc(fc)
          ffix_to_sminn(c) = (freelivfix_slope*(max(0._r8,AnnET(c))*secs_per_year) + freelivfix_inter )/secs_per_year !(units g N m-2 s-1)  

       end do

  end associate
  end subroutine CNFreeLivingFixation

  !-----------------------------------------------------------------------
  subroutine CNNFixation(num_soilc, filter_soilc, &
       cnveg_carbonflux_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen fixation rate
    ! as a function of annual total NPP. This rate gets updated once per year.
    ! All N fixation goes to the soil mineral N pool.
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
    use CNSharedParamsMod    , only: use_fun
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_carbonflux_type)            , intent(inout) :: cnveg_carbonflux_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc                  ! indices
    real(r8) :: t                     ! temporary
    real(r8) :: dayspyr               ! days per year
    !-----------------------------------------------------------------------

    associate(                                                                & 
         cannsum_npp    => cnveg_carbonflux_inst%annsum_npp_col ,             & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)                
         col_lag_npp    => cnveg_carbonflux_inst%lag_npp_col    ,             & ! Input: [real(r8) (:)]  (gC/m2/s) lagged net primary production           

         nfix_to_sminn  => soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col & ! Output: [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         )

      dayspyr = get_days_per_year()

      if ( nfix_timeconst > 0._r8 .and. nfix_timeconst < 500._r8 ) then
         ! use exponential relaxation with time constant nfix_timeconst for NPP - NFIX relation
         ! Loop through columns
         do fc = 1,num_soilc
            c = filter_soilc(fc)         

            if (col_lag_npp(c) /= spval) then
               ! need to put npp in units of gC/m^2/year here first
               t = (1.8_r8 * (1._r8 - exp(-0.003_r8 * col_lag_npp(c)*(secspday * dayspyr))))/(secspday * dayspyr)  
               nfix_to_sminn(c) = max(0._r8,t)
            else
               nfix_to_sminn(c) = 0._r8
            endif
         end do
      else
         ! use annual-mean values for NPP-NFIX relation
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            t = (1.8_r8 * (1._r8 - exp(-0.003_r8 * cannsum_npp(c))))/(secspday * dayspyr)
            nfix_to_sminn(c) = max(0._r8,t)
         end do
      endif
      if(use_fun)then
        nfix_to_sminn(c) = 0.0_r8
      end if

    end associate

  end subroutine CNNFixation
 
  !-----------------------------------------------------------------------
  subroutine CNNFert(bounds, num_soilc, filter_soilc, &
       cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen fertilizer for crops
    ! All fertilizer goes into the soil mineral N pool.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)    :: bounds  
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_nitrogenflux_type)          , intent(in)    :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: c,fc                 ! indices
    !-----------------------------------------------------------------------

    associate(                                                                  &   
         fert          =>    cnveg_nitrogenflux_inst%fert_patch ,               & ! Input:  [real(r8) (:)]  nitrogen fertilizer rate (gN/m2/s)                
         fert_to_sminn =>    soilbiogeochem_nitrogenflux_inst%fert_to_sminn_col & ! Output: [real(r8) (:)]                                                    
         )
      
      call p2c(bounds, num_soilc, filter_soilc, &
           fert(bounds%begp:bounds%endp), &
           fert_to_sminn(bounds%begc:bounds%endc))

    end associate

  end subroutine CNNFert

  !-----------------------------------------------------------------------
  subroutine CNSoyfix (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       waterstate_inst, crop_inst, cnveg_state_inst, cnveg_nitrogenflux_inst , &
       soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! This routine handles the fixation of nitrogen for soybeans based on
    ! the EPICPHASE model M. Cabelguenne et al., Agricultural systems 60: 175-196, 1999
    ! N-fixation is based on soil moisture, plant growth phase, and availibility of
    ! nitrogen in the soil root zone.
    !
    ! !USES:
    use pftconMod, only : ntmp_soybean, nirrig_tmp_soybean
    use pftconMod, only : ntrp_soybean, nirrig_trp_soybean
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(crop_type)                         , intent(in)    :: crop_inst
    type(cnveg_state_type)                  , intent(in)    :: cnveg_state_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_state_type)         , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: fp,p,c
    real(r8):: fxw,fxn,fxg,fxr             ! soil water factor, nitrogen factor, growth stage factor
    real(r8):: soy_ndemand                 ! difference between nitrogen supply and demand
    real(r8):: GDDfrac
    real(r8):: sminnthreshold1, sminnthreshold2
    real(r8):: GDDfracthreshold1, GDDfracthreshold2
    real(r8):: GDDfracthreshold3, GDDfracthreshold4
    !-----------------------------------------------------------------------

    associate(                                                                      & 
         wf               =>  waterstate_inst%wf_col                      ,         & ! Input:  [real(r8) (:) ]  soil water as frac. of whc for top 0.5 m          

         hui              =>  crop_inst%gddplant_patch                    ,         & ! Input:  [real(r8) (:) ]  gdd since planting (gddplant)                    
         croplive         =>  crop_inst%croplive_patch                    ,         & ! Input:  [logical  (:) ]  true if planted and not harvested                  

         gddmaturity      =>  cnveg_state_inst%gddmaturity_patch          ,         & ! Input:  [real(r8) (:) ]  gdd needed to harvest                             

         plant_ndemand    =>  cnveg_nitrogenflux_inst%plant_ndemand_patch ,         & ! Input:  [real(r8) (:) ]  N flux required to support initial GPP (gN/m2/s)  
         soyfixn          =>  cnveg_nitrogenflux_inst%soyfixn_patch       ,         & ! Output: [real(r8) (:) ]  nitrogen fixed to each soybean crop               

         fpg              =>  soilbiogeochem_state_inst%fpg_col           ,         & ! Input:  [real(r8) (:) ]  fraction of potential gpp (no units)              

         sminn            =>  soilbiogeochem_nitrogenstate_inst%sminn_col ,         & ! Input:  [real(r8) (:) ]  (kgN/m2) soil mineral N                           
         soyfixn_to_sminn =>  soilbiogeochem_nitrogenflux_inst%soyfixn_to_sminn_col & ! Output: [real(r8) (:) ]                                                    
         )

      sminnthreshold1 = 30._r8
      sminnthreshold2 = 10._r8
      GDDfracthreshold1 = 0.15_r8
      GDDfracthreshold2 = 0.30_r8
      GDDfracthreshold3 = 0.55_r8
      GDDfracthreshold4 = 0.75_r8

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)

         ! if soybean currently growing then calculate fixation

         if (croplive(p) .and. &
              (patch%itype(p) == ntmp_soybean .or. &
               patch%itype(p) == nirrig_tmp_soybean .or. &
               patch%itype(p) == ntrp_soybean .or. &
               patch%itype(p) == nirrig_trp_soybean) ) then

            ! difference between supply and demand

            if (fpg(c) < 1._r8) then
               soy_ndemand = 0._r8
               soy_ndemand = plant_ndemand(p) - plant_ndemand(p)*fpg(c)

               ! fixation depends on nitrogen, soil water, and growth stage

               ! soil water factor

               fxw = 0._r8
               fxw = wf(c)/0.85_r8

               ! soil nitrogen factor (Beth says: CHECK UNITS)

               if (sminn(c) > sminnthreshold1) then
                  fxn = 0._r8
               else if (sminn(c) > sminnthreshold2 .and. sminn(c) <= sminnthreshold1) then
                  fxn = 1.5_r8 - .005_r8 * (sminn(c) * 10._r8)
               else if (sminn(c) <= sminnthreshold2) then
                  fxn = 1._r8
               end if

               ! growth stage factor
               ! slevis: to replace GDDfrac, assume...
               ! Beth's crit_offset_gdd_def is similar to my gddmaturity
               ! Beth's ac_gdd (base 5C) similar to my hui=gddplant (base 10
               ! for soy) 
               ! Ranges below are not firm. Are they lit. based or tuning based?

               GDDfrac = hui(p) / gddmaturity(p)

               if (GDDfrac <= GDDfracthreshold1) then
                  fxg = 0._r8
               else if (GDDfrac > GDDfracthreshold1 .and. GDDfrac <= GDDfracthreshold2) then
                  fxg = 6.67_r8 * GDDfrac - 1._r8
               else if (GDDfrac > GDDfracthreshold2 .and. GDDfrac <= GDDfracthreshold3) then
                  fxg = 1._r8
               else if (GDDfrac > GDDfracthreshold3 .and. GDDfrac <= GDDfracthreshold4) then
                  fxg = 3.75_r8 - 5._r8 * GDDfrac
               else  ! GDDfrac > GDDfracthreshold4
                  fxg = 0._r8
               end if

               ! calculate the nitrogen fixed by the soybean

               fxr = min(1._r8, fxw, fxn) * fxg 
               fxr = max(0._r8, fxr)
               soyfixn(p) =  fxr * soy_ndemand
               soyfixn(p) = min(soyfixn(p), soy_ndemand)

            else ! if nitrogen demand met, no fixation

               soyfixn(p) = 0._r8

            end if

         else ! if not live soybean, no fixation

            soyfixn(p) = 0._r8

         end if
      end do

      call p2c(bounds, num_soilc, filter_soilc, &
           soyfixn(bounds%begp:bounds%endp), &
           soyfixn_to_sminn(bounds%begc:bounds%endc))

    end associate

  end subroutine CNSoyfix

end module CNNDynamicsMod
