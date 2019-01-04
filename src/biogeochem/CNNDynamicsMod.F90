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
         fluxes2(6,2), fluxes3(6,3), fluxes4(6,4), tanpools2(2), tanpools4(4), fluxes_tmp(6), garbage_total
    real(r8), parameter :: water_init_grz = 0.006_r8, cnc_nh3_air = 0.0_r8, depth_slurry = 0.005_r8
    !real(r8), parameter :: fract_resist=0.225_r8, fract_unavail=0.025_r8, fract_avail=0.25_r8, fract_tan=0.6_r8

    real(r8), parameter :: fract_tan=0.5_r8 ! of all N
    real(r8), parameter :: fract_resist=0.45_r8, fract_unavail=0.05_r8, fract_avail=0.5_r8 ! of organic N

    
    real(r8), parameter :: dz_layer_fert = 0.02_r8, dz_layer_grz = 0.02_r8
    !real(r8), parameter :: fract_resist=0._r8, fract_unavail=0._r8, fract_avail=0._r8, fract_tan=1.0_r8
    real(r8), parameter :: fert_incorp_reduct = 0.3_r8
    real(r8), parameter :: slurry_infiltr_time = 12*3600.0_r8, water_init_fert = 1e-6
    real(r8), parameter :: &
         poolranges_grz(3) = (/24*3600.0_r8, 10*24*3600.0_r8, 360*24*3600.0_r8/), &
         poolranges_fert(3) = (/2.36*24*3600.0_r8, 24*3600.0_r8, 360*24*3600.0_r8/), &
         poolranges_slr(4) = (/slurry_infiltr_time, 24*3600.0_r8, 10*24*3600.0_r8, 360*24*3600.0_r8/), &
         !Hconc_grz(3) = (/10**(-8.5_r8), 10**(-8.0_r8), 10**(-7.0_r8)/), &
         Hconc_fert(3) = (/10**(-7.0_r8), 10**(-8.5_r8), 10**(-8.0_r8)/)

    real(r8) :: Hconc_grz(3), Hconc_slr(4), pH_soil, pH_crop
    real(r8) :: fert_inc_tan, fert_inc_no3
    !logical, parameter :: do_balance_checks = .false.
    logical :: do_balance_checks
    real(r8) :: tg, garbage, theta, thetasat, infiltr_m_s, evap_m_s, runoff_m_s, org_n_tot, &
         nstored_old, nsoilman_old, nsoilfert_old, fert_to_air, fert_to_soil, fert_total, fert_urea, fert_tan, &
         soilflux_org, urea_resid
    real(r8) :: tanprod_from_urea(3), ureapools(2), fert_no3, fert_generic, bsw
    !real(r8), parameter :: fract_urea=0.545, fract_no3=0.048
    real(r8) :: fract_urea, fract_no3, soilph_min, soilph_max, soilpsi
    integer, parameter :: ind_region = 1
    integer :: def_ph_count
    
    Hconc_grz(1:2) = (/10**(-8.5_r8), 10**(-8.0_r8)/)
    Hconc_slr(1:3) = (/10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8)/)

    soilph_min = 999
    soilph_max = -999
    def_ph_count = 0
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
          ngrz(c) = atm2lnd_inst%forc_ndep_grz_grc(g) / col%wtgcell(c) * 1e3 ! kg to g 
          if (debug_fan) then
             if (ngrz(c) > 1e12 .or. (isnan(ngrz(c)))) then
                write(iulog, *) 'bad ngrz', atm2lnd_inst%forc_ndep_grz_grc(g), col%wtgcell(c)
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

    call handle_storage_v2(bounds, temperature_inst, frictionvel_inst, dt, &
         atm2lnd_inst%forc_ndep_sgrz_grc, atm2lnd_inst%forc_ndep_ngrz_grc, &
         ns%man_n_stored_col, ns%man_tan_stored_col, &
         nf%man_n_appl_col, nf%man_tan_appl_col, &
         nf%man_n_grz_col, nf%man_n_mix_col, &
         nf%nh3_stores_col, nf%nh3_barns_col, &
         nf%man_n_transf_col, ns%fan_grz_fract_col, &
         nf%man_n_barns_col, &
         fract_tan, &
         filter_soilc, num_soilc)

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
       g = col%gridcell(c)
       if (.not. (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)) cycle
       if (.not. col%active(c) .or. col%wtgcell(c) < 1e-15) cycle

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
          ns%fan_grz_fract_col(c) = 1.0_r8 ! for crops handled by handle_storage
       end if

       ! Calculation of the water fluxes should include the background soil moisture
       ! tendency. However, it is unclear how to do this in a numerically consistent
       ! way. Following a naive finite differencing approach led to worse agreement in
       ! stand-alone simulations so the term is currenltly neglected here.
       watertend = 0.0_r8

       ! use the calculated tend
       watertend = waterstate_inst%h2osoi_tend_tsl_col(c) * 1e-3 ! to meters/sec (ie. m3/m2/s)
       
       tg = temperature_inst%t_grnd_col(c)
       theta = waterstate_inst%h2osoi_vol_col(c,1)
       thetasat = soilstate_inst%watsat_col(c,1)
       bsw = soilstate_inst%bsw_col(c,1)
       theta = min(theta, 0.98_r8*thetasat)
       infiltr_m_s = max(waterflux_inst%qflx_infl_col(c), 0.0) * 1e-3 
       evap_m_s = waterflux_inst%qflx_evap_grnd_col(c) * 1e-3
       runoff_m_s = max(waterflux_inst%qflx_runoff_col(c), 0.0) * 1e-3
       soilpsi = soilstate_inst%soilpsi_col(c,1)

       !
       ! grazing
       !

       ndep_org(ind_avail) = ngrz(c) * (1.0_r8-fract_tan) * fract_avail
       ndep_org(ind_resist) = ngrz(c) * (1.0_r8-fract_tan) * fract_resist
       ndep_org(ind_unavail) = ngrz(c) * (1.0_r8-fract_tan) * fract_unavail
       tandep = ngrz(c) * fract_tan

       orgpools(ind_avail) = man_a_grz(c)
       orgpools(ind_resist) = man_r_grz(c)
       orgpools(ind_unavail) = man_u_grz(c)
       call update_org_n(ndep_org, tg, soilpsi, orgpools, dt, tanprod, soilflux_org)
       man_a_grz(c) = orgpools(ind_avail)
       man_r_grz(c) = orgpools(ind_resist) 
       man_u_grz(c) = orgpools(ind_unavail)

       tanpools3(1) = ns%tan_g1_col(c)
       tanpools3(2) = ns%tan_g2_col(c)
       tanpools3(3) = ns%tan_g3_col(c)
       if (any(isnan(tanpools3))) then
          call endrun('nan1')
       end if

       ph_soil = atm2lnd_inst%forc_soilph_grc(g)
       if (ph_soil < 3.0) then
          ph_soil = 6.5_r8
          def_ph_count = def_ph_count + 1
       end if
       Hconc_grz(3) = 10**-(ph_soil)
       soilph_max = max(soilph_max, ph_soil)
       soilph_min = min(soilph_min, ph_soil)

       fluxes_tmp = 0.0
       garbage_total = 0.0
       fluxes3 = 0.0
       garbage = 0
       do ind_substep = 1, num_substeps
          call update_npool(tg, ratm, &
               theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, tandep, (/0.0_r8, 0.0_r8, sum(tanprod)/), water_init_grz, &
               bsw, poolranges_grz, Hconc_grz, dz_layer_grz, tanpools3, &
               fluxes3(1:5,:), garbage, dt/num_substeps, status, 3)
          if (status /= 0) then
             write(iulog, *) 'status = ', status, tanpools2, ratm, theta, thetasat, tandep, tanprod
             call endrun(msg='update_npool status /= 0')
          end if
          if (debug_fan .and. any(isnan(tanpools2))) then
             call endrun('nan2')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes3, dim=2)
          garbage_total = garbage_total + garbage
       end do
       fluxes_tmp = fluxes_tmp / num_substeps

       ns%tan_g1_col(c) = tanpools3(1)
       ns%tan_g2_col(c) = tanpools3(2)
       ns%tan_g3_col(c) = tanpools3(3)
       if (debug_fan .and. any(isnan(fluxes3))) then
          write(iulog, *) fluxes3
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
       ndep_org(ind_avail) = org_n_tot * fract_avail
       ndep_org(ind_resist) = org_n_tot * fract_resist 
       ndep_org(ind_unavail) = org_n_tot * fract_unavail 
       tandep = nf%man_tan_appl_col(c)

       orgpools(ind_avail) = man_a_app(c)
       orgpools(ind_resist) = man_r_app(c)
       orgpools(ind_unavail) = man_u_app(c)
       call update_org_n(ndep_org, tg, soilpsi, orgpools, dt, tanprod, soilflux_org)
       man_a_app(c) = orgpools(ind_avail)
       man_r_app(c) = orgpools(ind_resist)
       man_u_app(c) = orgpools(ind_unavail)
       tanpools4(1) = ns%tan_s0_col(c)
       tanpools4(2) = ns%tan_s1_col(c)
       tanpools4(3) = ns%tan_s2_col(c)
       tanpools4(4) = ns%tan_s3_col(c)

       ph_crop = min(max(ph_soil, 5.5_r8), 7.5_r8)
       Hconc_slr(4) = 10**-(ph_crop)

       if (debug_fan .and. any(isnan(tanpools4))) then
          call endrun('nan31')
       end if

       fluxes_tmp = 0.0
       garbage_total = 0.0
       fluxes4 = 0.0
       do ind_substep = 1, num_substeps
          if (debug_fan .and. any(abs(tanpools4) > 1e12)) then
             write(iulog, *) ind_substep, tanpools4, tandep, nf%fert_n_appl_col(c), &
                  nf%man_n_appl_col(c), ns%man_n_stored_col(c), ns%man_tan_stored_col(c)
             call endrun('bad tanpools (manure app)')
          end if

          call update_4pool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, tandep, sum(tanprod), bsw, depth_slurry, &
               poolranges_slr, tanpools4, Hconc_slr, fluxes4(1:5,:), garbage, dt / num_substeps, status)
          if (status /= 0) then
             write(iulog, *) 'status = ', status, tanpools4, tg, ratm, 'th', theta, &
                  thetasat, tandep, 'tp', tanprod, 'fx', fluxes4
             call endrun(msg='update_3pool status /= 0')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes4, dim=2)
          garbage_total = garbage_total + garbage
       end do
       fluxes_tmp = fluxes_tmp / num_substeps

       ns%tan_s0_col(c) = tanpools4(1)
       ns%tan_s1_col(c) = tanpools4(2)
       ns%tan_s2_col(c) = tanpools4(3)
       ns%tan_s3_col(c) = tanpools4(4)

       if (debug_fan .and. any(isnan(fluxes4))) then
          write(iulog, *) fluxes3, tanpools4,ratm, theta, thetasat, infiltr_m_s, tandep, tanprod
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

       ! Fraction available for volatilization 
       fert_total = nf%fert_n_appl_col(c)
       
       fract_urea = atm2lnd_inst%forc_ndep_urea_grc(g)
       fract_no3 = atm2lnd_inst%forc_ndep_nitr_grc(g)

       ! Fractions made unavailable by mechanical incorporation, will be added to the
       ! to-soil flux (tan) or no3 production (no3) below.
       fert_inc_tan = fert_total * fert_incorp_reduct * (1.0 - fract_no3)
       fert_inc_no3 = fert_total * fert_incorp_reduct * fract_no3
       
       if (fract_urea < 0 .or. fract_no3 < 0 .or. fract_urea + fract_no3 > 1) then
          call endrun('bad fertilizer fractions')
       end if
       
       fert_urea = fert_total * fract_urea * (1.0_r8 - fert_incorp_reduct)
       
       ! Include the incorporated NO3 fertilizer to the no3 flux
       fert_no3 = fert_total * fract_no3
       
       !fert_generic = (fert_total - fert_urea - fert_no3) * (1.0_r8 - fert_incorp_reduct)
       fert_generic = fert_total * (1.0_r8 - fract_urea - fract_no3) * (1.0_r8 - fert_incorp_reduct)
       
       nf%otherfert_n_appl_col(c) = fert_total * (1.0_r8 - fract_urea) !fert_no3 + fert_generic
       
       ! Urea decomposition 
       ! 
       ureapools(1) = ns%fert_u0_col(c)
       ureapools(2) = ns%fert_u1_col(c)
       fluxes2 = 0.0
       call update_urea(tg, theta, thetasat, infiltr_m_s, evap_m_s, watertend, &
            runoff_m_s, fert_urea, bsw, ureapools,  fluxes2, urea_resid, poolranges_fert(1:2), &
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
       nf%nh3_otherfert_col(c) = 0.0
       do ind_substep = 1, num_substeps
          ! Fertilizer pools f0...f2
          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, 0.0_r8, tanprod_from_urea, water_init_fert, bsw, &
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
               runoff_m_s, fert_generic, (/0.0_r8/), water_init_fert, bsw, &
               !(/360*24*3600.0_r8/), (/10**(-6.0_r8)/), dz_layer_fert, ns%tan_f3_col(c:c), fluxes3(1:5,1:1), &
               (/360*24*3600.0_r8/), (/10**(-ph_crop)/), dz_layer_fert, ns%tan_f3_col(c:c), fluxes3(1:5,1:1), &
               garbage, dt/num_substeps, status, numpools=1)
          if (status /= 0) then
             write(iulog, *) 'status:', status, tanpools3, nf%fert_n_appl_col(c)
             call endrun(msg='Bad status after npool for generic')
          end if
          fluxes_tmp = fluxes_tmp + fluxes3(:, 1) / num_substeps
          garbage_total = garbage_total + garbage
          nf%nh3_otherfert_col(c) = nf%nh3_otherfert_col(c) + fluxes3(iflx_air, 1) / num_substeps
       end do

       ns%tan_f0_col(c) = tanpools3(1)
       ns%tan_f1_col(c) = tanpools3(2)
       ns%tan_f2_col(c) = tanpools3(3)
       ! !!tan_f3_col already updated above by update_npool!!

       nf%nh3_fert_col(c) = fluxes_tmp(iflx_air)
       nf%fert_runoff_col(c) = fluxes_tmp(iflx_roff)
       nf%fert_no3_prod_col(c) = fluxes_tmp(iflx_no3) + fert_no3
       nf%fert_nh4_to_soil_col(c) = fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) + garbage_total/dt + fert_inc_tan

       ! Total flux
       ! 
       nf%nh3_total_col(c) = nf%nh3_fert_col(c) + nf%nh3_man_app_col(c) &
            + nf%nh3_grz_col(c) + nf%nh3_stores_col(c) +  nf%nh3_barns_col(c)
       if (nf%nh3_total_col(c) < -1e15) then
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
       write(iulog, *) 'SoilPH check:', soilph_min, soilph_max, def_ph_count
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
         total = total + sum(ns%tan_g1_col(soilc)) + sum(ns%tan_g2_col(soilc)) + sum(ns%tan_g3_col(soilc)) 
         total = total + sum(ns%man_u_grz_col(soilc)) &
              + sum(ns%man_a_grz_col(soilc)) + sum(ns%man_r_grz_col(soilc))
         total = total + sum(ns%tan_s0_col(soilc)) &
              + sum(ns%tan_s1_col(soilc)) + sum(ns%tan_s2_col(soilc)) + sum(ns%tan_s3_col(soilc))
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

  subroutine handle_storage_v2(bounds, temperature_inst, frictionvel_inst, dt,  &
       ndep_sgrz_grc, ndep_ngrz_grc, n_stored_col, tan_stored_col, &
       n_manure_spread_col, tan_manure_spread_col, &
       n_manure_graze_col, n_manure_mixed_col, &
       nh3_flux_stores, nh3_flux_barns, man_n_transf, &
       grz_fract, man_n_barns, tan_fract_excr, &
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
    
    ! N excreted in manure, gN/m2:
    real(r8), intent(in) :: ndep_sgrz_grc(bounds%begg:bounds%endg) ! seasonally grazing animals
    real(r8), intent(in) :: ndep_ngrz_grc(bounds%begg:bounds%endg) ! non-grazing animals
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
    real(r8), intent(out) :: man_n_barns(bounds%begc:bounds%endc)
    ! fraction of manure excreted when grazing
    real(r8), intent(out) :: grz_fract(bounds%begc:bounds%endc)
    ! TAN fraction in excreted N
    real(r8), intent(in) :: tan_fract_excr
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
    real(r8), parameter :: fract_continuous = 0.1_r8, kg_to_g = 1e3_r8, max_grazing_fract = 0.65_r8, &
         volat_coef_barns = 0.03_r8, volat_coef_stores = 0.025_r8, &
         tempr_min_grazing = 283.0_r8!!!!

    begg = bounds%begg; endg = bounds%endg
    nh3_flux_stores(bounds%begc:bounds%endc) = 0_r8
    nh3_flux_barns(bounds%begc:bounds%endc) = 0_r8
    man_n_barns(bounds%begc:bounds%endc) = 0.0_r8
    
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

                n_manure_mixed_col(c) = (ndep_ngrz_grc(g) + ndep_sgrz_grc(g)) * kg_to_g / lun%wtgcell(l)
                
                tempr_min_10day = temperature_inst%t_a10min_patch(col%patchi(c))
                if (tempr_min_10day > tempr_min_grazing) then
                   ! fraction of animals grazing -> allocate some manure to grasslands before barns
                   flux_grazing = max_grazing_fract * ndep_sgrz_grc(g) * kg_to_g / lun%wtgcell(l)
                   flux_avail = (ndep_ngrz_grc(g) + ndep_sgrz_grc(g)*(1.0_r8 - max_grazing_fract)) * kg_to_g / lun%wtgcell(l)
                   grz_fract(c) = max_grazing_fract
                else
                   flux_grazing = 0.0_r8
                   flux_avail = n_manure_mixed_col(c)
                   grz_fract(c) = 0.0_r8
                end if
                flux_grass_graze = flux_grass_graze + flux_grazing*col%wtgcell(c)

                if (flux_avail > 1e12 .or. isnan(flux_avail)) then
                   write(iulog, *) 'bad flux_avail', ndep_ngrz_grc(g), ndep_sgrz_grc(g), lun%wtgcell(l)
                   call endrun('bad flux_avail')
                end if

                totalinput = totalinput + flux_avail

                counter = 0
                if (col_grass == c) call endrun('Something wrong with the indices')
                if (col%patchi(c) /= col%patchf(c)) then
                   call endrun(msg="ERROR crop column has multiple patches")
                end if

                tempr_ave = temperature_inst%t_ref2m_patch(col%patchi(c))
                windspeed_ave = frictionvel_inst%u10_patch(col%patchi(c))

                man_n_barns(c) = flux_avail
                
                call eval_fluxes_storage(flux_avail, tempr_ave, windspeed_ave, 0.0_r8, &
                     volat_coef_barns, volat_coef_stores, tan_fract_excr, fluxes_nitr, fluxes_tan, status)
                if (any(fluxes_nitr > 1e12)) then
                   write(iulog, *) 'bad fluxes', fluxes_nitr
                end if
                if (status /=0) then 
                   write(iulog, *) 'status = ', status
                   call endrun(msg='eval_fluxes_storage failed')
                end if
                cumflux = cumflux + sum(fluxes_nitr)
                
                if (fluxes_tan(iflx_to_store) < 0) then
                   call endrun(msg="ERROR too much manure lost")
                end if

                flux_grass_spread = flux_grass_spread + fluxes_nitr(iflx_to_store)*col%wtgcell(c)
                flux_grass_spread_tan = flux_grass_spread_tan + fluxes_tan(iflx_to_store)*col%wtgcell(c)

                man_n_transf(c) = flux_grazing + fluxes_nitr(iflx_to_store)
                
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
       else if (flux_grass_spread > 0) then
          call endrun('Cannot spread manure')
       end if

    end do ! grid

  end subroutine handle_storage_v2
  

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
