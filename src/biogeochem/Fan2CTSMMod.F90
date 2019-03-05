module Fan2CTSMMod
  use FanMod
  use shr_kind_mod, only : r8 => shr_kind_r8, CL => shr_kind_cl
  use decompMod                       , only : bounds_type
  use atm2lndType                     , only : atm2lnd_type
  use Wateratm2lndBulkType            , only : wateratm2lndbulk_type
  use TemperatureType                 , only : temperature_type
  use FrictionVelocityMod             , only : frictionvel_type
  use shr_infnan_mod                  , only : isnan => shr_infnan_isnan
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use CNVegNitrogenFluxType	      , only : cnveg_nitrogenflux_type
  use WaterStateBulkType              , only : waterstatebulk_type
  use WaterFluxBulkType               , only : waterfluxbulk_type
  use SoilStateType                   , only : soilstate_type
  use ColumnType                      , only : col                
  use PatchType                       , only : patch                
  use clm_varctl                      , only : iulog

  implicit none

  private

  public fan_readnml
  public fan_eval
  public fan_to_sminn
  
  ! Hydrogen ion concentration in TAN pools, mol/l == 10**-pH
  !
  ! Pastures and slurry. The last age class gets soil pH.
  real(r8), parameter :: Hconc_grz_def(2) = (/10**(-8.5_r8), 10**(-8.0_r8)/)
  real(r8), parameter :: Hconc_slr_def(3) &
       = (/10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8)/)
  ! Urea fertilizer. The other fertilizer (F4) pool gets soil pH.
  real(r8), parameter :: Hconc_fert(3) = (/10**(-7.0_r8), 10**(-8.5_r8), 10**(-8.0_r8)/)

  ! Active layer thickness used by FAN. 
  real(r8), parameter :: dz_layer_fert = 0.02_r8 ! m
  real(r8), parameter :: dz_layer_grz = 0.02_r8  ! m

  ! Manure N composition
  real(r8) :: fract_tan = 0.6_r8 ! fraction of total ammoniacal nitrogen
  ! The following are fractions of non-TAN N
  real(r8), parameter :: &
       fract_resist = 0.45_r8, & ! resistant organic N
       fract_unavail = 0.05_r8, & ! unvavailable organic N 
       fract_avail = 0.5_r8 ! available organic N
  
  ! application rate in meters water:
  real(r8), parameter :: water_init_grz = 0.006_r8 ! urine patch depth (m)
  real(r8), parameter :: depth_slurry = 0.005_r8   ! slurry application rate (m)
  real(r8), parameter :: water_init_fert = 1e-9_r8   ! water in fertilizer (assumed none).
  ! Slurry infiltration time
  real(r8), parameter :: slurry_infiltr_time = 6.0_r8*3600_r8 ! seconds
  ! Reduction factor for fertilizer due to mechanical incorporation.
  ! N available for volatilization becomes multiplied by (1-fert_incorp_reduct).
  real(r8) :: fert_incorp_reduct = 0.25_r8 
  
  ! TAN pool age ranges (sec). 
  real(r8), parameter ::          &
       poolranges_grz(3) = (/24*3600.0_r8, 10*24*3600.0_r8, 360*24*3600.0_r8/), &
       poolranges_fert(3) = (/2.36*24*3600.0_r8, 24*3600.0_r8, 360*24*3600.0_r8/), &
       poolranges_slr(4) = (/slurry_infiltr_time, 24*3600.0_r8, 10*24*3600.0_r8, 360*24*3600.0_r8/), &
       poolrange_otherfert(1) = (/360*24*3600.0_r8/)

  ! soil pH for crops restricted between these limits: 
  real(r8), parameter :: pH_crop_min = 5.5_r8
  real(r8), parameter :: pH_crop_max = 7.5_r8

  ! Parameters for grazing in mixed/landless systems:
  real(r8), parameter :: tempr_min_grazing = 283.0_r8 ! Lowest 10-day daily-min temperature for grazing, K
  ! Fraction of ruminants grazing when permitted by temperature
  real(r8), parameter :: max_grazing_fract = 0.65_r8
  ! Normalization constants for barn and storage emissions.
  real(r8), parameter :: volat_coef_barns_open = 0.03_r8, volat_coef_barns_closed = 0.025, volat_coef_stores = 0.025_r8

  ! Fraction of manure N moved from crop to native columns (manure spreading)
  real(r8) :: fract_spread_grass = 1.0_r8
  
  logical :: fan_to_bgc_crop = .false.
  logical :: fan_to_bgc_veg = .false.
  
  logical, parameter :: debug_fan = .true.
  
contains

  subroutine fan_readnml(NLFilename)
    use spmdMod        , only : masterproc, mpicom
    use fileutils      , only : getavu, relavu, opnfil
    use clm_varctl     , only : use_fan
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use shr_nl_mod     , only : shr_nl_find_group_name
    use abortutils     , only : endrun
    use shr_mpi_mod    , only : shr_mpi_bcast
    use FanStreamMod   , only : set_bcast_fanstream_pars
    
    character(len=*), intent(in) :: NLFilename ! Namelist filename

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: subname = 'fan_readnml'
    character(len=*), parameter :: nmlname = 'fan_nml'
    integer :: stream_year_first_fan      ! first year in stream to use
    integer :: stream_year_last_fan       ! last year in stream to use
    integer :: model_year_align_fan       ! align stream_year_firstndep2 with 
    character(len=CL)  :: stream_fldFileName_fan
    character(len=CL)  :: fan_mapalgo

    namelist /fan_nml/ fan_to_bgc_crop, fan_to_bgc_veg, stream_year_first_fan,  &
         stream_year_last_fan, model_year_align_fan, fan_mapalgo, stream_fldFileName_fan, &
         fract_spread_grass

    if (masterproc) then
       unitn = getavu()
       write(iulog, *) 'Read in ' // nmlname // '  namelist'
       call opnfil(NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=fan_nml, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading " // nmlname // "namelist" // errmsg(__FILE__, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find " // nmlname // "namelist" // errmsg(__FILE__, __LINE__))
       end if
       call relavu(unitn)
    end if

    call set_bcast_fanstream_pars(stream_year_first_fan,  stream_year_last_fan, &
         model_year_align_fan, fan_mapalgo, stream_fldFileName_fan)

    call shr_mpi_bcast(fan_to_bgc_crop, mpicom)
    call shr_mpi_bcast(fan_to_bgc_veg, mpicom)
    call shr_mpi_bcast(fract_spread_grass, mpicom)

    if (fract_spread_grass > 1 .or. fract_spread_grass < 0) then
       call endrun(msg="ERROR invalid fract_spread_grass")
    end if

    
    !call mpi_bcast(fan_to_bgc_crop, 1, MPI_LOGICAL, 0, mpicom, ierr)
    !call mpi_bcast(fan_to_bgc_veg, 1, MPI_LOGICAL, 0, mpicom, ierr)
    !call mpi_bcast(use_fan, 1, MPI_LOGICAL, 0, mpicom, ierr)
    
  end subroutine fan_readnml

  !************************************************************************************
  
  subroutine fan_eval(bounds, num_soilc, filter_soilc, &
       atm2lnd_inst, wateratm2lndbulk_inst, &
       cnveg_nitrogenflux_inst, &
       soilbiogeochem_nitrogenflux_inst, &
       soilbiogeochem_nitrogenstate_inst, &
       waterstatebulk_inst, soilstate_inst, temperature_inst, &
       waterfluxbulk_inst, frictionvel_inst)
    use clm_time_manager, only: get_step_size, get_curr_date, get_curr_calday, get_nstep
    use clm_varpar, only: max_patch_per_col
    use LandunitType, only: lun
    use shr_sys_mod, only : shr_sys_flush
    use GridcellType, only: grc
    use abortutils, only : endrun
    use pftconMod, only : nc4_grass, nc3_nonarctic_grass
    use landunit_varcon, only:  istsoil, istcrop
    use clm_varcon, only : spval, ispval
    use decompMod, only : bounds_type
    use subgridAveMod, only: p2c

    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_inst
    type(wateratm2lndbulk_type), intent(in)  :: wateratm2lndbulk_inst
    type(cnveg_nitrogenflux_type)          , intent(in) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(waterstatebulk_type)              , intent(in)    :: waterstatebulk_inst
    type(soilstate_type)                   , intent(in)    :: soilstate_inst
    type(temperature_type)                 , intent(in) :: temperature_inst
    type(waterfluxbulk_type)               , intent(in)    :: waterfluxbulk_inst
    type(frictionvel_type)                 , intent(in) :: frictionvel_inst

    ! Local variables
    integer, parameter :: num_substeps = 4, balance_check_freq = 1
    integer :: c, g, patchcounter, p, status, c1, c2, l, fc, ind_substep
    real(r8) :: dt, ndep_org(3), orgpools(3), tanprod(3), watertend, fluxes(6,3), tanpools3(3), ratm, tandep, &
         fluxes2(6,2), fluxes3(6,3), fluxes4(6,4), tanpools2(2), tanpools4(4), fluxes_tmp(6), garbage_total

    real(r8) :: Hconc_grz(3), Hconc_slr(4), pH_soil, pH_crop
    real(r8) :: fert_inc_tan, fert_inc_no3, old(num_soilc)
    logical :: do_balance_checks
    real(r8) :: tg, garbage, theta, thetasat, infiltr_m_s, evap_m_s, runoff_m_s, org_n_tot, &
         nstored_old, nsoilman_old, nsoilfert_old, fert_to_air, fert_to_soil, fert_total, fert_urea, fert_tan, &
         soilflux_org, urea_resid
    real(r8) :: tanprod_from_urea(3), ureapools(2), fert_no3, fert_generic, bsw
    real(r8) :: fract_urea, fract_no3, soilph_min, soilph_max, soilpsi
    integer :: def_ph_count
    
    Hconc_grz(1:2) = Hconc_grz_def
    Hconc_slr(1:3) = Hconc_slr_def

    soilph_min = 999
    soilph_max = -999
    def_ph_count = 0
    dt = real(get_step_size(), r8)
    do_balance_checks = mod(get_nstep(), balance_check_freq) == 0

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

    !call p2c(bounds, num_soilc, filter_soilc, &
    !     cnv_nf%manu_patch(bounds%begp:bounds%endp), &
    !     nf%man_n_appl_col(bounds%begc:bounds%endc))
    nf%man_n_appl_col(bounds%begc:bounds%endc) = 0.0_r8
    
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
               cnv_nf%manure_patch(col%patchi(c):col%patchf(c))
          call endrun('nf%man_n_appl_col(c) is spval')
       end if

       ! Find and average the atmospheric resistances Rb and Ra.
       ! 
       if (lun%itype(col%landunit(c)) == istcrop) then
          ! Crop column, only one patch
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
                ! Nothing found. We shouldn't be here.
                ratm = 150.0_r8
             end if
          end if
          ns%fan_grz_fract_col(c) = 1.0_r8 ! for crops handled by handle_storage
       end if

       watertend = waterstatebulk_inst%h2osoi_tend_tsl_col(c) * 1e-3 ! to meters/sec (ie. m3/m2/s)
       
       tg = temperature_inst%t_grnd_col(c)
       theta = waterstatebulk_inst%h2osoi_vol_col(c,1)
       thetasat = soilstate_inst%watsat_col(c,1)
       bsw = soilstate_inst%bsw_col(c,1)
       theta = min(theta, 0.98_r8*thetasat)
       infiltr_m_s = max(waterfluxbulk_inst%qflx_infl_col(c), 0.0) * 1e-3 
       evap_m_s = waterfluxbulk_inst%qflx_evap_grnd_col(c) * 1e-3
       runoff_m_s = max(waterfluxbulk_inst%qflx_runoff_col(c), 0.0) * 1e-3
       soilpsi = soilstate_inst%soilpsi_col(c,1)

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
               wateratm2lndbulk_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, tandep, (/0.0_r8, 0.0_r8, sum(tanprod)/), water_init_grz, &
               bsw, poolranges_grz, Hconc_grz, dz_layer_grz, tanpools3, &
               fluxes3(1:5,:), garbage, dt/num_substeps, status, 3)
          if (status /= 0) then
             write(iulog, *) 'status = ', status, tanpools3, ratm, theta, thetasat, tandep, tanprod
             call endrun(msg='update_npool status /= 0')
          end if
          if (debug_fan .and. any(isnan(tanpools3))) then
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

       ph_crop = min(max(ph_soil, ph_crop_min), ph_crop_max)
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
               wateratm2lndbulk_inst%forc_q_downscaled_col(c), watertend, &
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
               wateratm2lndbulk_inst%forc_q_downscaled_col(c), watertend, &
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
               wateratm2lndbulk_inst%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, fert_generic, (/0.0_r8/), water_init_fert, bsw, &
               poolrange_otherfert, (/10**(-ph_crop)/), dz_layer_fert, ns%tan_f3_col(c:c), fluxes3(1:5,1:1), &
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

    old = ns%fan_totn_col(filter_soilc(1:num_soilc))

    call update_summary(ns, nf, filter_soilc, num_soilc)

    call debug_balance(ns, nf, old, filter_soilc(1:num_soilc))
    
    call debug_balance(ns, nf, old, (/1/))
    
    call debug_balance(ns, nf, old, (/2/))
    
    end associate

  contains

    subroutine debug_balance(ns, nf, oldn, columns)
      type(soilbiogeochem_nitrogenstate_type), intent(in) :: ns
      type(soilbiogeochem_nitrogenflux_type), intent(in) :: nf
      real(r8) :: oldn(:)
      integer :: columns(:)

      real(r8) :: newn(size(old))

      newn = ns%fan_totn_col
      
    
      print *, 'FAN SUMMARY', columns
      print *, 'old total:', sum(oldn(columns))
      print *, 'new total:', sum(newn(columns))
      print *, 'delta:', sum(oldn(columns)) - sum(newn(columns))
      print *, 'new flux:', (sum(nf%fan_totnin(columns)) - sum(nf%fan_totnout(columns)))*dt
      !print *, 'total from funct', get_total_n(ns, nf, 'pools_manure') + get_total_n(ns, nf, 'pools_fertilizer')
      
    end subroutine debug_balance
      
    
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
    
  end subroutine fan_eval

  !************************************************************************************
  
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

    implicit none
    type(bounds_type), intent(in)    :: bounds
    type(temperature_type) , intent(in) :: temperature_inst
    type(frictionvel_type) , intent(in) :: frictionvel_inst
    real(r8), intent(in) :: dt
    
    ! N excreted in manure, gN/m2:
    real(r8), intent(in) :: ndep_sgrz_grc(bounds%begg:bounds%endg) ! seasonally grazing animals
    real(r8), intent(in) :: ndep_ngrz_grc(bounds%begg:bounds%endg) ! non-grazing animals
    real(r8), intent(inout) :: n_stored_col(bounds%begc:bounds%endc), &
         & tan_stored_col(bounds%begc:bounds%endc) ! N, TAN currently stored, gN/m2
    ! N, TAN spread on grasslands, gN/m2/s:
    real(r8), intent(inout) :: n_manure_spread_col(bounds%begc:bounds%endc) 
    real(r8), intent(out) :: tan_manure_spread_col(bounds%begc:bounds%endc) ! output, calculated from the above and stored manure
    ! N excreted by animals allocated to mixed production systems temporarily grazing on grasslands:
    real(r8), intent(inout) :: n_manure_graze_col(bounds%begc:bounds%endc)
    ! N excreted by animals in mixed systems, total
    real(r8), intent(out) :: n_manure_mixed_col(bounds%begc:bounds%endc)
    ! NH3 emission fluxes from manure storage and housings, gN/m2/s
    real(r8), intent(out) :: nh3_flux_stores(bounds%begc:bounds%endc), nh3_flux_barns(bounds%begc:bounds%endc)
    ! total nitrogen flux transferred out of a crop column (manure spreading + temporary grazing)
    real(r8), intent(out) :: man_n_transf(bounds%begc:bounds%endc)
    ! Total nitrogen excreted in barns
    real(r8), intent(out) :: man_n_barns(bounds%begc:bounds%endc)
    ! fraction of manure excreted when grazing
    real(r8), intent(out) :: grz_fract(bounds%begc:bounds%endc)
    ! TAN fraction in excreted N
    real(r8), intent(in) :: tan_fract_excr
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns

    integer :: begg, endg, g, l, c, il, counter, col_grass, status, p
    real(r8) :: flux_avail_rum, flux_avail_mg, flux_grazing
    real(r8) :: tempr_ave, windspeed_ave ! windspeed and temperature averaged over agricultural patches
    real(r8) :: tempr_barns, tempr_stores, vent_barns, flux_grass_crop, tempr_min_10day, &
         flux_grass_graze, flux_grass_spread, flux_grass_spread_tan, flux_grass_crop_tan
    real(r8) :: cumflux, totalinput, total_to_store, total_to_store_tan
    real(r8) :: fluxes_nitr(4,2), fluxes_tan(4,2)
    ! The fraction of manure applied continuously on grasslands (if present in the gridcell)
    real(r8), parameter :: kg_to_g = 1e3_r8
    
    begg = bounds%begg; endg = bounds%endg
    nh3_flux_stores(bounds%begc:bounds%endc) = 0_r8
    nh3_flux_barns(bounds%begc:bounds%endc) = 0_r8
    man_n_barns(bounds%begc:bounds%endc) = 0.0_r8
    
    totalinput = 0.0
    cumflux = 0.0
    
    do g = begg, endg
       ! First find out if there are grasslands in this cell. If yes, a fraction of
       ! manure can be diverted to them before storage.
       col_grass = ispval
       do il = 1, max_lunit
          l = grc%landunit_indices(il, g)
          if (l == ispval) cycle
          if (lun%itype(l) == istsoil) then
             do p = lun%patchi(l), lun%patchf(l)
                if ((patch%itype(p) == nc4_grass .or. patch%itype(p) == nc3_nonarctic_grass) &
                    .and. col%wtgcell(patch%column(p)) > 1e-6) then
                   col_grass = patch%column(p)
                   exit
                end if
             end do
          end if
          if (col_grass /= ispval) exit
       end do
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
                   !flux_avail = (ndep_ngrz_grc(g) + ndep_sgrz_grc(g)*(1.0_r8 - max_grazing_fract)) * kg_to_g / lun%wtgcell(l)
                   flux_avail_rum = (ndep_sgrz_grc(g)*(1.0_r8 - max_grazing_fract)) * kg_to_g / lun%wtgcell(l)
                   grz_fract(c) = max_grazing_fract
                else
                   flux_grazing = 0.0_r8
                   flux_avail_rum = ndep_sgrz_grc(g) * kg_to_g / lun%wtgcell(l)
                   grz_fract(c) = 0.0_r8
                end if
                flux_avail_mg = ndep_ngrz_grc(g) * kg_to_g / lun%wtgcell(l)
                flux_grass_graze = flux_grass_graze + flux_grazing*col%wtgcell(c)

                if (flux_avail_rum > 1e12 .or. flux_avail_mg > 1e12 .or. isnan(flux_avail_mg) .or. isnan(flux_avail_rum)) then
                   write(iulog, *) 'bad flux_avail', ndep_ngrz_grc(g), ndep_sgrz_grc(g), lun%wtgcell(l)
                   call endrun('bad flux_avail')
                end if

                totalinput = totalinput + flux_avail_rum + flux_avail_mg

                counter = 0
                if (col_grass == c) call endrun('Something wrong with the indices')
                if (col%patchi(c) /= col%patchf(c)) then
                   call endrun(msg="ERROR crop column has multiple patches")
                end if

                tempr_ave = temperature_inst%t_ref2m_patch(col%patchi(c))
                windspeed_ave = frictionvel_inst%u10_patch(col%patchi(c))

                man_n_barns(c) = flux_avail_rum + flux_avail_mg
                
                ! Evaluate the NH3 losses, separate for ruminants (open barns) and others
                ! (poultry and pigs, closed barns). Note the slicing of fluxes(:,:) and fluxes_tan(:,:).
                man_n_transf(c) = flux_grazing
                nh3_flux_stores(c) = 0.0

                if (flux_avail_rum < 0) then 
                   write(iulog, *) 'flux:', flux_avail_rum
                   call endrun(msg='negat flux_avail for ruminants')
                end if

                ! Ruminants
                call eval_fluxes_storage(flux_avail_rum, 'open', tempr_ave, windspeed_ave, 0.0_r8, &
                     volat_coef_barns_open, volat_coef_stores, tan_fract_excr, fluxes_nitr(:,1), fluxes_tan(:,1), status)
                if (status /=0) then 
                   write(iulog, *) 'status = ', status
                   call endrun(msg='eval_fluxes_storage failed for ruminants')
                end if

                ! Others
                call eval_fluxes_storage(flux_avail_mg, 'closed', tempr_ave, windspeed_ave, 0.0_r8, &
                     volat_coef_barns_closed, volat_coef_stores, tan_fract_excr, fluxes_nitr(:,2), fluxes_tan(:,2), status)
                if (status /=0) then 
                   write(iulog, *) 'status = ', status
                   call endrun(msg='eval_fluxes_storage failed for other livestock')
                end if
                cumflux = cumflux + sum(fluxes_nitr)
                if (any(isnan(fluxes_nitr))) then
                   write(iulog, *) 'fluxes 1', fluxes_nitr(:,1)
                   write(iulog, *) 'fluxes 2', fluxes_nitr(:,2)
                   call endrun('Nan in fluxes_nitr')
                end if
                if (any(isnan(fluxes_tan))) then
                   write(iulog, *) 'fluxes 1', fluxes_tan(:,1)
                   write(iulog, *) 'fluxes 2', fluxes_tan(:,2)
                   call endrun('Nan in fluxes_tan')
                end if
                
                if (fluxes_tan(iflx_to_store,1) < 0) then
                   call endrun(msg="ERROR too much manure lost")
                end if
                ! Simplification as of 2019: no explicit manure storage. Flux to storage
                ! will be spread "immediately".
                total_to_store = sum(fluxes_nitr(iflx_to_store,:))
                total_to_store_tan = sum(fluxes_tan(iflx_to_store,:))

                n_manure_spread_col(c) = (1.0_r8 - fract_spread_grass) * total_to_store
                tan_manure_spread_col(c) = (1.0_r8 - fract_spread_grass) * total_to_store_tan

                flux_grass_spread = flux_grass_spread + fract_spread_grass * total_to_store*col%wtgcell(c)
                flux_grass_spread_tan = flux_grass_spread_tan + fract_spread_grass * total_to_store_tan*col%wtgcell(c)
                
                man_n_transf(c) = man_n_transf(c) + total_to_store

                nh3_flux_stores(c) = sum(fluxes_nitr(iflx_air_stores,:))
                nh3_flux_barns(c) = sum(fluxes_nitr(iflx_air_barns,:))
                
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
          continue
          !call endrun('Cannot spread manure')
       end if

    end do ! grid

  end subroutine handle_storage_v2

  !************************************************************************************
  
  subroutine update_summary(ns, nf, filter_soilc, num_soilc)
    use ColumnType, only : col
    use LandunitType   , only: lun
    use landunit_varcon, only : istcrop

    type(soilbiogeochem_nitrogenstate_type), intent(inout) :: ns
    type(soilbiogeochem_nitrogenflux_type), intent(inout) :: nf
    integer, intent(in)    :: num_soilc       ! number of soil columns in filter
    integer, intent(in)    :: filter_soilc(:) ! filter for soil columns

    integer :: c, fc
    real(r8) :: total, fluxout, fluxin, flux_loss    
    
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       total = ns%tan_g1_col(c) + ns%tan_g2_col(c) + ns%tan_g3_col(c)
       total = total + ns%man_u_grz_col(c) + ns%man_a_grz_col(c) + ns%man_r_grz_col(c)
       total = total + ns%tan_s0_col(c) + ns%tan_s1_col(c) + ns%tan_s2_col(c) + ns%tan_s3_col(c)
       total = total + ns%man_u_app_col(c) + ns%man_a_app_col(c) + ns%man_r_app_col(c)
       total = total + ns%tan_f0_col(c) + ns%tan_f1_col(c) + ns%tan_f2_col(c) + ns%tan_f3_col(c)
       total = total + ns%fert_u0_col(c) + ns%fert_u1_col(c)
       ns%fan_totn_col(c) = total

       if (lun%itype(col%landunit(c)) == istcrop) then
          ! no grazing, man_n_appl is from the same column and not counted as input
          fluxin = nf%man_n_barns_col(c) + nf%fert_n_appl_col(c)
       else
          ! no barns or fertilization. man_n_appl is transferred from crop columns and not
          ! included in the other inputs.
          fluxin = nf%man_n_grz_col(c) + nf%man_n_appl_col(c)
       end if

       flux_loss = nf%nh3_man_app_col(c) + nf%nh3_grz_col(c) + nf%manure_runoff_col(c) &
            + nf%nh3_stores_col(c) + nf%nh3_barns_col(c) &
            + nf%nh3_fert_col(c) + nf%fert_runoff_col(c)
              
       fluxout = nf%fert_no3_prod_col(c) + nf%fert_nh4_to_soil_col(c) &
            + nf%manure_no3_prod_col(c) + nf%manure_nh4_to_soil_col(c) &
            + nf%man_n_transf_col(c) + flux_loss

       nf%fan_totnin(c) = fluxin
       nf%fan_totnout(c) = fluxout
       
    end do
    
  end subroutine update_summary

  !************************************************************************************
  
  subroutine fan_to_sminn(filter_soilc, num_soilc, sbgc_nf)
    use ColumnType, only : col
    use LandunitType   , only: lun
    use landunit_varcon, only : istcrop, istsoil
    
    integer, intent(in) :: filter_soilc(:)
    integer, intent(in) :: num_soilc
    type(soilbiogeochem_nitrogenflux_type), intent(inout) :: sbgc_nf

    integer :: c, fc
    real(r8) :: fan_nitr
    
    do fc = 1, num_soilc
       c = filter_soilc(c)
       fan_nitr &
            = sbgc_nf%fert_no3_prod_col(c) + sbgc_nf%fert_nh4_to_soil_col(c) &
            + sbgc_nf%manure_no3_prod_col(c) + sbgc_nf%manure_nh4_to_soil_col(c)
       if (lun%itype(col%landunit(c)) == istcrop .and. fan_to_bgc_crop) then
          sbgc_nf%fert_to_sminn_col(c) = fan_nitr
       else if (lun%itype(col%landunit(c)) == istsoil .and. fan_to_bgc_veg) then
          sbgc_nf%fert_to_sminn_col(c) = fan_nitr
       end if
    end do
    
  end subroutine fan_to_sminn
  
end module Fan2CTSMMod
