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
  use clm_varctl                      , only : use_nitrif_denitrif, nfix_timeconst
  use subgridAveMod                   , only : p2c
  use atm2lndType                     , only : atm2lnd_type
  use Wateratm2lndBulkType            , only : wateratm2lndbulk_type
  use CNVegStateType                  , only : cnveg_state_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use TemperatureType                 , only : temperature_type
  use FrictionVelocityMod             , only : frictionvel_type
  use clm_varctl                      , only : iulog
  use shr_infnan_mod                  , only : isnan => shr_infnan_isnan
  use CNVegNitrogenStateType	      , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType	      , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use WaterStateBulkType              , only : waterstatebulk_type
  use WaterFluxBulkType               , only : waterfluxbulk_type
  use SoilStateType                   , only : soilstate_type
  use WaterDiagnosticBulkType         , only : waterdiagnosticbulk_type
  use WaterFluxBulkType               , only : waterfluxbulk_type
  use CropType                        , only : crop_type
  use ColumnType                      , only : col                
  use PatchType                       , only : patch                
  use perf_mod                        , only : t_startf, t_stopf
  use FanMod
  use clm_varctl                      , only: use_fan

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
       atm2lnd_inst, wateratm2lndbulk_inst, &
       soilbiogeochem_nitrogenflux_inst, cnveg_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_carbonflux_inst, &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
       waterstatebulk_inst, soilstate_inst, temperature_inst, &
       waterfluxbulk_inst, frictionvel_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen deposition rate from
    ! atmospheric forcing and if enabled, call FAN to evaluate NH3 volatilization from
    ! fertilizers and manure.
    !
    ! For now it is assumed that all the atmospheric N deposition goes to the soil mineral
    ! N pool.  This could be updated later to divide the inputs between mineral N absorbed
    ! directly into the canopy and mineral N entering the soil pool.
    !
    ! !USES:
    use CNSharedParamsMod    , only: use_fun
    use clm_time_manager     , only: get_step_size, get_curr_date, get_curr_calday, get_nstep
    !
    ! !LOCAL VARIABLES:
    use clm_varpar           , only: max_patch_per_col
    use LandunitType         , only: lun
    use shr_sys_mod          , only : shr_sys_flush
    use GridcellType         , only: grc
    use clm_varctl           , only : iulog
    use abortutils           , only : endrun
    use pftconMod            , only : nc4_grass, nc3_nonarctic_grass
    use landunit_varcon      , only:  istsoil, istcrop
    use clm_varcon           , only : spval, ispval
    use Fan2ctsmMod
    
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_inst
    type(wateratm2lndbulk_type), intent(in)  :: wateratm2lndbulk_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)         , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)          , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type), intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_carbonflux_type)   , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(waterstatebulk_type)              , intent(in)    :: waterstatebulk_inst
    type(soilstate_type)                   , intent(in)    :: soilstate_inst
    type(temperature_type)                 , intent(inout) :: temperature_inst
    type(waterfluxbulk_type)               , intent(in)    :: waterfluxbulk_inst
    type(frictionvel_type)                 , intent(inout) :: frictionvel_inst

    integer :: c, g
    
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
    
    if (use_fan) then
       call fan_eval(bounds, num_soilc, filter_soilc, &
            atm2lnd_inst, wateratm2lndbulk_inst, &
            cnveg_nitrogenflux_inst, &
            soilbiogeochem_nitrogenflux_inst, &
            soilbiogeochem_nitrogenstate_inst, &
            waterstatebulk_inst, soilstate_inst, temperature_inst, &
            waterfluxbulk_inst, frictionvel_inst)
    end if
    
  end subroutine CNNDeposition

  !-----------------------------------------------------------------------
  subroutine CNFreeLivingFixation(num_soilc, filter_soilc, &
       waterfluxbulk_inst, soilbiogeochem_nitrogenflux_inst)


    use clm_time_manager , only : get_days_per_year
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
 
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter                                                                                                                     
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns                                                                                                                                  
   
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    type(waterfluxbulk_type)                   , intent(inout) :: waterfluxbulk_inst 
    !
    ! !LOCAL VARIABLES:                                                                                                                                                                                                           
    integer  :: c,fc            !indices     
    real(r8) :: dayspyr         !days per year 
    real(r8) :: secs_per_year   !seconds per year   

       associate(                                                                        &
                  AnnET            => waterfluxbulk_inst%AnnET,                              & ! Input:  [real(:)  ] : Annual average ET flux mmH20/s
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
    use clm_time_manager , only : get_days_per_year
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
       num_pcropp, filter_pcropp, &
       cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst)
    use Fan2CTSMMod, only : fan_to_sminn, fan_to_bgc_crop
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
    integer                                , intent(in)    :: num_pcropp      ! number of prognostic crop pathches
    integer                                , intent(in)    :: filter_pcropp(:)! filter for prognostic crop patches
    type(cnveg_nitrogenflux_type)          , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: c,fc, p, fp                 ! indices
    real(r8) :: manure_col(bounds%begc:bounds%endc)
    
    !-----------------------------------------------------------------------

    associate(                                                                  &   
         synthfert     =>    cnveg_nitrogenflux_inst%synthfert_patch ,          & ! Input:  [real(r8) (:)]  nitrogen fertilizer rate (gN/m2/s)
         manure        =>    cnveg_nitrogenflux_inst%manure_patch ,             & ! Input:  [real(r8) (:)]  manure nitrogen rate (gN/m2/s)
         totalfert     =>    cnveg_nitrogenflux_inst%nfertilization_patch,      & ! Input:  [real(r8) (:)]  manure nitrogen rate (gN/m2/s)                
         fert_to_sminn =>    soilbiogeochem_nitrogenflux_inst%fert_to_sminn_col & ! Output: [real(r8) (:)]                                                    
         )

      if (.not. fan_to_bgc_crop) then
         ! => Crop columns/patches are not handled by FAN. Use synthfert directly and add
         ! the default CLM manure. No N input to non-crop columns in this case.
         call p2c(bounds, num_soilc, filter_soilc, &
              synthfert(bounds%begp:bounds%endp), &
              fert_to_sminn(bounds%begc:bounds%endc))
         call p2c(bounds, num_soilc, filter_soilc, &
              manure(bounds%begp:bounds%endp), &
              manure_col(bounds%begc:bounds%endc))
         ! Add the manure N processed above:
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            fert_to_sminn(c) = fert_to_sminn(c) + manure_col(c)
         end do
         ! Add up synthetic fertilizer and manure to the nfertilization output variable.
         do fp = 1, num_pcropp
            p = filter_pcropp(fp)
            totalfert(p) = synthfert(p) + manure(p)
         end do
      end if

      ! if fan_to_bgc_crop == .true., FAN fills in the fert_to_sminn and totalfert for
      ! crops. It might also fill in the non-crop columns if enabled.
      call fan_to_sminn(bounds, filter_soilc, num_soilc, soilbiogeochem_nitrogenflux_inst, totalfert)
      
    end associate
    
  end subroutine CNNFert

  !-----------------------------------------------------------------------
  subroutine CNSoyfix (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       waterdiagnosticbulk_inst, crop_inst, cnveg_state_inst, cnveg_nitrogenflux_inst , &
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
    type(waterdiagnosticbulk_type)                   , intent(in)    :: waterdiagnosticbulk_inst
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
         wf               =>  waterdiagnosticbulk_inst%wf_col                      ,         & ! Input:  [real(r8) (:) ]  soil water as frac. of whc for top 0.5 m          

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
