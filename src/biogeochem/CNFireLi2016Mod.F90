module CNFireLi2016Mod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for fire dynamics 
  ! created in Nov, 2012  and revised in Apr, 2013 by F. Li and S. Levis
  ! based on Li et al. (2012a,b; 2013)
  ! revised in Apr, 2014 according to Li et al.(2014)
  ! revised in May, 2015, according to Li et al. (2015, in prep.)
  ! Fire-related parameters were calibrated or tuned in May, 2015 based on the 
  ! 20th Century transient simulations at f19_g16 with a CLM4.5 version
  ! (clm50fire), CRUNCEPv5, and climatological lightning data.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_const_mod                      , only : SHR_CONST_PI,SHR_CONST_TKFRZ
  use shr_infnan_mod                     , only : shr_infnan_isnan
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varctl                         , only : iulog
  use clm_varpar                         , only : nlevdecomp, ndecomp_pools, nlevdecomp_full
  use clm_varcon                         , only : dzsoi_decomp
  use pftconMod                          , only : noveg, pftcon
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use subgridAveMod                      , only : p2c
  use atm2lndType                        , only : atm2lnd_type
  use CNDVType                           , only : dgvs_type
  use CNVegStateType                     , only : cnveg_state_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType             , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType              , only : cnveg_nitrogenflux_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use EnergyFluxType                     , only : energyflux_type
  use SaturatedExcessRunoffMod           , only : saturated_excess_runoff_type
  use WaterDiagnosticBulkType                     , only : waterdiagnosticbulk_type
  use Wateratm2lndBulkType                     , only : wateratm2lndbulk_type
  use GridcellType                       , only : grc                
  use ColumnType                         , only : col                
  use PatchType                          , only : patch                
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  use CNFireMethodMod                    , only : cnfire_method_type
  use CNFireBaseMod                      , only : cnfire_base_type, cnfire_const
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: cnfire_li2016_type
  !
  type, extends(cnfire_base_type) :: cnfire_li2016_type
     private
  contains
     !
     ! !PUBLIC MEMBER FUNCTIONS:
     procedure, public :: CNFireArea    ! Calculate fire area
  end type cnfire_li2016_type

  !
  ! !PRIVATE MEMBER DATA:
  !-----------------------------------------------------------------------

  interface cnfire_li2016_type
     ! initialize a new cnfire_base object
     module procedure constructor
  end interface cnfire_li2016_type
  !-----------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !------------------------------------------------------------------------
  type(cnfire_li2016_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type cnfire_base_type.
    ! !ARGUMENTS:

    constructor%need_lightning_and_popdens = .true.
  end function constructor

  !-----------------------------------------------------------------------
  subroutine CNFireArea (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       atm2lnd_inst, energyflux_inst, saturated_excess_runoff_inst, waterdiagnosticbulk_inst, &
       wateratm2lndbulk_inst, cnveg_state_inst, cnveg_carbonstate_inst, totlitc_col, decomp_cpools_vr_col, t_soi17cm_col)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area 
    !
    ! !USES:
    use clm_time_manager     , only: get_step_size, get_days_per_year, get_curr_date, get_nstep
    use clm_varpar           , only: max_patch_per_col
    use clm_varcon           , only: secspday, secsphr
    use clm_varctl           , only: spinup_state
    use pftconMod            , only: nc4_grass, nc3crop, ndllf_evr_tmp_tree
    use pftconMod            , only: nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree, nbrdlf_evr_shrub
    use dynSubgridControlMod , only : run_has_transient_landcover
    !
    ! !ARGUMENTS:
    class(cnfire_li2016_type)                             :: this
    type(bounds_type)                     , intent(in)    :: bounds 
    integer                               , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(atm2lnd_type)                    , intent(in)    :: atm2lnd_inst
    type(energyflux_type)                 , intent(in)    :: energyflux_inst
    type(saturated_excess_runoff_type)    , intent(in)    :: saturated_excess_runoff_inst
    type(waterdiagnosticbulk_type)                 , intent(in)    :: waterdiagnosticbulk_inst
    type(wateratm2lndbulk_type)                 , intent(in)    :: wateratm2lndbulk_inst
    type(cnveg_state_type)                , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)          , intent(inout) :: cnveg_carbonstate_inst
    real(r8)                              , intent(in)    :: totlitc_col(bounds%begc:)
    real(r8)                              , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:)
    real(r8)                              , intent(in)    :: t_soi17cm_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    !
    integer  :: g,l,c,p,pi,j,fc,fp,kyr, kmo, kda, mcsec   ! index variables
    real(r8) :: dt       ! time step variable (s)
    real(r8) :: dayspyr  ! days per year
    real(r8) :: cli      ! effect of climate on deforestation fires (0-1)
    real(r8) :: cri      ! thresholds used for cli, (mm/d), see Eq.(7) in Li et al.(2013)
    real(r8) :: fb       ! availability of fuel for regs A and C
    real(r8) :: fhd      ! impact of hd on agricultural fire
    real(r8) :: fgdp     ! impact of gdp on agricultural fire
    real(r8) :: fire_m   ! combustability of fuel for fire occurrence
    real(r8) :: spread_m ! combustability of fuel for fire spread
    real(r8) :: Lb_lf    ! length-to-breadth ratio added by Lifang 
    integer  :: i_cwd    ! cwd pool
    real(r8) :: lh       ! anthro. ignitions (count/km2/hr)
    real(r8) :: fs       ! hd-dependent fires suppression (0-1)
    real(r8) :: ig       ! total ignitions (count/km2/hr)
    real(r8) :: hdmlf    ! human density
    real(r8) :: arh, arh30 !combustability of fuel related to RH and RH30
    real(r8) :: afuel    !weight for arh and arh30
    real(r8) :: btran_col(bounds%begc:bounds%endc)
    logical  :: transient_landcover  ! whether this run has any prescribed transient landcover
    real(r8), target  :: prec60_col_target(bounds%begc:bounds%endc)
    real(r8), target  :: prec10_col_target(bounds%begc:bounds%endc)
    real(r8), target  :: rh30_col_target(bounds%begc:bounds%endc) 
    real(r8), pointer :: prec60_col(:)
    real(r8), pointer :: prec10_col(:)
    real(r8), pointer :: rh30_col(:)
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(totlitc_col)           == (/bounds%endc/))                              , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_vr_col)  == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soi17cm_col)         == (/bounds%endc/))                              , errMsg(sourcefile, __LINE__))

    associate(                                                                      & 
         totlitc            => totlitc_col                                     , & ! Input:  [real(r8) (:)     ]  (gC/m2) total lit C (column-level mean)           
         decomp_cpools_vr   => decomp_cpools_vr_col                            , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  VR decomp. (litter, cwd, soil)
         tsoi17             => t_soi17cm_col                                   , & ! Input:  [real(r8) (:)     ]  (K) soil T for top 0.17 m                             

         lfuel              => cnfire_const%lfuel                              , & ! Input:  [real(r8)         ]  (gC/m2) Lower threshold of fuel mass
         ufuel              => cnfire_const%ufuel                              , & ! Input:  [real(r8)         ]  (gC/m2) Upper threshold of fuel mass
         rh_hgh             => cnfire_const%rh_hgh                             , & ! Input:  [real(r8)         ]  (%) High relative humidity
         rh_low             => cnfire_const%rh_low                             , & ! Input:  [real(r8)         ]  (%) Low relative humidity
         bt_min             => cnfire_const%bt_min                             , & ! Input:  [real(r8)         ]  (0-1) Minimum btran
         bt_max             => cnfire_const%bt_max                             , & ! Input:  [real(r8)         ]  (0-1) Maximum btran
         cli_scale          => cnfire_const%cli_scale                          , & ! Input:  [real(r8)         ]  (/d) global constant for deforestation fires
         cropfire_a1        => cnfire_const%cropfire_a1                        , & ! Input:  [real(r8)         ]  (/hr) a1 parameter for cropland fire
         non_boreal_peatfire_c => cnfire_const%non_boreal_peatfire_c           , & ! Input:  [real(r8)         ]  (/hr) c parameter for non-boreal peatland fire
         pot_hmn_ign_counts_alpha => cnfire_const%pot_hmn_ign_counts_alpha     , & ! Input:  [real(r8)         ]  (/person/month) Potential human ignition counts
         boreal_peatfire_c  => cnfire_const%boreal_peatfire_c                  , & ! Input:  [real(r8)         ]  (/hr) c parameter for boreal peatland fire
         
         fsr_pft            => pftcon%fsr_pft                                  , & ! Input:
         fd_pft             => pftcon%fd_pft                                   , & ! Input:

         btran2             => energyflux_inst%btran2_patch                    , & ! Input:  [real(r8) (:)     ]  root zone soil wetness                            
         fsat               => saturated_excess_runoff_inst%fsat_col           , & ! Input:  [real(r8) (:)     ]  fractional area with water table at surface       
         wf2                => waterdiagnosticbulk_inst%wf2_col                         , & ! Input:  [real(r8) (:)     ]  soil water as frac. of whc for top 0.17 m         
         
         is_cwd             => decomp_cascade_con%is_cwd                       , & ! Input:  [logical  (:)     ]  TRUE => pool is a cwd pool                         
         spinup_factor      => decomp_cascade_con%spinup_factor                , & ! Input:  [real(r8) (:)     ]  factor for AD spinup associated with each pool           

         forc_rh            => wateratm2lndbulk_inst%forc_rh_grc                        , & ! Input:  [real(r8) (:)     ]  relative humidity                                 
         forc_wind          => atm2lnd_inst%forc_wind_grc                      , & ! Input:  [real(r8) (:)     ]  atmospheric wind speed (m/s)                       
         forc_t             => atm2lnd_inst%forc_t_downscaled_col              , & ! Input:  [real(r8) (:)     ]  downscaled atmospheric temperature (Kelvin)                  
         forc_rain          => wateratm2lndbulk_inst%forc_rain_downscaled_col           , & ! Input:  [real(r8) (:)     ]  downscaled rain                                              
         forc_snow          => wateratm2lndbulk_inst%forc_snow_downscaled_col           , & ! Input:  [real(r8) (:)     ]  downscaled snow                                              
         prec60             => wateratm2lndbulk_inst%prec60_patch                       , & ! Input:  [real(r8) (:)     ]  60-day running mean of tot. precipitation         
         prec10             => wateratm2lndbulk_inst%prec10_patch                       , & ! Input:  [real(r8) (:)     ]  10-day running mean of tot. precipitation         
         rh30               => wateratm2lndbulk_inst%rh30_patch                         , & ! Input:  [real(r8) (:)     ]  10-day running mean of tot. precipitation 
         dwt_smoothed       => cnveg_state_inst%dwt_smoothed_patch             , & ! Input:  [real(r8) (:)     ]  change in patch weight (-1 to 1) on the gridcell, smoothed over the year
         cropf_col          => cnveg_state_inst%cropf_col                      , & ! Input:  [real(r8) (:)     ]  cropland fraction in veg column                   
         gdp_lf             => cnveg_state_inst%gdp_lf_col                     , & ! Input:  [real(r8) (:)     ]  gdp data                                          
         peatf_lf           => cnveg_state_inst%peatf_lf_col                   , & ! Input:  [real(r8) (:)     ]  peatland fraction data                            
         abm_lf             => cnveg_state_inst%abm_lf_col                     , & ! Input:  [integer  (:)     ]  prescribed crop fire time                          
         baf_crop           => cnveg_state_inst%baf_crop_col                   , & ! Output: [real(r8) (:)     ]  burned area fraction for cropland (/sec)  
         baf_peatf          => cnveg_state_inst%baf_peatf_col                  , & ! Output: [real(r8) (:)     ]  burned area fraction for peatland (/sec)  
         burndate           => cnveg_state_inst%burndate_patch                 , & ! Output: [integer  (:)     ]  burn date for crop                                 
         fbac               => cnveg_state_inst%fbac_col                       , & ! Output: [real(r8) (:)     ]  total burned area out of conversion (/sec)
         fbac1              => cnveg_state_inst%fbac1_col                      , & ! Output: [real(r8) (:)     ]  burned area out of conversion region due to land use fire
         farea_burned       => cnveg_state_inst%farea_burned_col               , & ! Output: [real(r8) (:)     ]  total fractional area burned (/sec)
         nfire              => cnveg_state_inst%nfire_col                      , & ! Output: [real(r8) (:)     ]  fire counts (count/km2/sec), valid only in Reg. C
         fsr_col            => cnveg_state_inst%fsr_col                        , & ! Output: [real(r8) (:)     ]  fire spread rate at column level
         fd_col             => cnveg_state_inst%fd_col                         , & ! Output: [real(r8) (:)     ]  fire duration rate at column level
         lgdp_col           => cnveg_state_inst%lgdp_col                       , & ! Output: [real(r8) (:)     ]  gdp limitation factor for nfire                   
         lgdp1_col          => cnveg_state_inst%lgdp1_col                      , & ! Output: [real(r8) (:)     ]  gdp limitation factor for baf per fire            
         lpop_col           => cnveg_state_inst%lpop_col                       , & ! Output: [real(r8) (:)     ]  pop limitation factor for baf per fire            
         lfwt               => cnveg_state_inst%lfwt_col                       , & ! Output: [real(r8) (:)     ]  fractional coverage of non-crop and non-bare-soil Patches
         trotr1_col         => cnveg_state_inst%trotr1_col                     , & ! Output: [real(r8) (:)     ]  patch weight of BET on the column (0-1)           
         trotr2_col         => cnveg_state_inst%trotr2_col                     , & ! Output: [real(r8) (:)     ]  patch weight of BDT on the column (0-1)           
         dtrotr_col         => cnveg_state_inst%dtrotr_col                     , & ! Output: [real(r8) (:)     ]  decreased frac. coverage of BET+BDT on grid for dt
         lfc                => cnveg_state_inst%lfc_col                        , & ! Output: [real(r8) (:)     ]  conversion area frac. of BET+BDT that haven't burned before
         wtlf               => cnveg_state_inst%wtlf_col                       , & ! Output: [real(r8) (:)     ]  fractional coverage of non-crop Patches              
         
         totvegc            => cnveg_carbonstate_inst%totvegc_col              , & ! Input: [real(r8) (:)     ]  totvegc at column level                           
         deadcrootc         => cnveg_carbonstate_inst%deadcrootc_patch         , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
         deadcrootc_storage => cnveg_carbonstate_inst%deadcrootc_storage_patch , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
         deadcrootc_xfer    => cnveg_carbonstate_inst%deadcrootc_xfer_patch    , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
         deadstemc          => cnveg_carbonstate_inst%deadstemc_patch          , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem root C
         frootc             => cnveg_carbonstate_inst%frootc_patch             , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C                               
         frootc_storage     => cnveg_carbonstate_inst%frootc_storage_patch     , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
         frootc_xfer        => cnveg_carbonstate_inst%frootc_xfer_patch        , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
         livecrootc         => cnveg_carbonstate_inst%livecrootc_patch         , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
         livecrootc_storage => cnveg_carbonstate_inst%livecrootc_storage_patch , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
         livecrootc_xfer    => cnveg_carbonstate_inst%livecrootc_xfer_patch    , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
         leafc              => cnveg_carbonstate_inst%leafc_patch              , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
         leafc_storage      => cnveg_carbonstate_inst%leafc_storage_patch      , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
         leafc_xfer         => cnveg_carbonstate_inst%leafc_xfer_patch         , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
         rootc_col          => cnveg_carbonstate_inst%rootc_col                , & ! Output: [real(r8) (:)     ]  root carbon                                       
         leafc_col          => cnveg_carbonstate_inst%leafc_col                , & ! Output: [real(r8) (:)     ]  leaf carbon at column level                       
         deadstemc_col      => cnveg_carbonstate_inst%deadstemc_col            , & ! Output: [real(r8) (:)     ]  deadstem carbon at column level
         fuelc              => cnveg_carbonstate_inst%fuelc_col                , & ! Output: [real(r8) (:)     ]  fuel load coutside cropland                 
         fuelc_crop         => cnveg_carbonstate_inst%fuelc_crop_col             & ! Output: [real(r8) (:)     ]  fuel load for cropland                 
         )
 
      transient_landcover = run_has_transient_landcover()

      !pft to column average 
      prec10_col =>prec10_col_target
      call p2c(bounds, num_soilc, filter_soilc, &
           prec10(bounds%begp:bounds%endp), &
           prec10_col(bounds%begc:bounds%endc))

      prec60_col =>prec60_col_target
      call p2c(bounds, num_soilc, filter_soilc, &
           prec60(bounds%begp:bounds%endp), &
           prec60_col(bounds%begc:bounds%endc))

      rh30_col =>rh30_col_target
      call p2c(bounds, num_soilc, filter_soilc, &
           rh30(bounds%begp:bounds%endp), &
           rh30_col(bounds%begc:bounds%endc))  
      
      call p2c(bounds, num_soilc, filter_soilc, &
           leafc(bounds%begp:bounds%endp), &
           leafc_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           deadstemc(bounds%begp:bounds%endp), &
           deadstemc_col(bounds%begc:bounds%endc))

     call get_curr_date (kyr, kmo, kda, mcsec)
     dayspyr = get_days_per_year()
     ! Get model step size
     dt      = real( get_step_size(), r8 )
     !
     ! On first time-step, just set area burned to zero and exit
     !
     if ( get_nstep() == 0 )then
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           farea_burned(c) = 0._r8
           baf_crop(c)     = 0._r8
           baf_peatf(c)    = 0._r8
           fbac(c)         = 0._r8
           fbac1(c)        = 0._r8
           cropf_col(c)    = 0._r8 
        end do
        return
     end if
     !
     ! Calculate fraction of crop (cropf_col) and non-crop and non-bare-soil 
     ! vegetation (lfwt) in vegetated column
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        cropf_col(c) = 0._r8 
        lfwt(c)      = 0._r8   
     end do
     do pi = 1,max_patch_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if (pi <=  col%npatches(c)) then
              p = col%patchi(c) + pi - 1
              ! For crop veg types
              if( patch%itype(p) > nc4_grass )then
                 cropf_col(c) = cropf_col(c) + patch%wtcol(p)
              end if
              ! For natural vegetation (non-crop and non-bare-soil)
              if( patch%itype(p) >= ndllf_evr_tmp_tree .and. patch%itype(p) <= nc4_grass )then
                 lfwt(c) = lfwt(c) + patch%wtcol(p)
              end if
           end if
        end do
     end do
     ! 
     ! Calculate crop fuel   
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        fuelc_crop(c)=0._r8
     end do
     do pi = 1,max_patch_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if (pi <=  col%npatches(c)) then
              p = col%patchi(c) + pi - 1
              ! For crop PFTs, fuel load includes leaf and litter; only
              ! column-level litter carbon
              ! is available, so we use leaf carbon to estimate the
              ! litter carbon for crop PFTs
              if( patch%itype(p) > nc4_grass .and. patch%wtcol(p) > 0._r8 .and. leafc_col(c) > 0._r8 )then
                 fuelc_crop(c)=fuelc_crop(c) + (leafc(p) + leafc_storage(p) + &
                      leafc_xfer(p))*patch%wtcol(p)/cropf_col(c)     + &
                      totlitc(c)*leafc(p)/leafc_col(c)*patch%wtcol(p)/cropf_col(c)
              end if
           end if
        end do
     end do
     !   
     ! Calculate noncrop column variables
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        fsr_col(c)   = 0._r8
        fd_col(c)    = 0._r8
        rootc_col(c) = 0._r8
        lgdp_col(c)  = 0._r8
        lgdp1_col(c) = 0._r8
        lpop_col(c)  = 0._r8
        btran_col(c) = 0._r8
        wtlf(c)      = 0._r8
        trotr1_col(c)= 0._r8
        trotr2_col(c)= 0._r8
        if (transient_landcover) then
           dtrotr_col(c)=0._r8
        end if
     end do
     do pi = 1,max_patch_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           g = col%gridcell(c)
           if (pi <=  col%npatches(c)) then
              p = col%patchi(c) + pi - 1

              ! For non-crop -- natural vegetation and bare-soil
              if( patch%itype(p)  <  nc3crop .and. cropf_col(c)  <  1.0_r8 )then
                 if( .not. shr_infnan_isnan(btran2(p))) then
                    if (btran2(p)  <=  1._r8 ) then
                       btran_col(c) = btran_col(c)+btran2(p)*patch%wtcol(p)
                       wtlf(c)      = wtlf(c)+patch%wtcol(p)
                    end if
                 end if

                 ! NOTE(wjs, 2016-12-15) These calculations of the fraction of evergreen
                 ! and deciduous tropical trees (used to determine if a column is
                 ! tropical closed forest) use the current fractions. However, I think
                 ! they are used in code that applies to land cover change. Note that
                 ! land cover change is currently generated on the first time step of the
                 ! year (even though the fire code sees the annually-smoothed dwt). Thus,
                 ! I think that, for this to be totally consistent, this code should
                 ! consider the fractional coverage of each PFT prior to the relevant
                 ! land cover change event. (These fractions could be computed in the
                 ! code that handles land cover change, so that the fire code remains
                 ! agnostic to exactly how and when land cover change happens.)
                 !
                 ! For example, if a year started with fractional coverages of
                 ! nbrdlf_evr_trp_tree = 0.35 and nbrdlf_dcd_trp_tree = 0.35, but then
                 ! the start-of-year land cover change reduced both of these to 0.2: The
                 ! current code would consider the column to NOT be tropical closed
                 ! forest (because nbrdlf_evr_trp_tree+nbrdlf_dcd_trp_tree < 0.6),
                 ! whereas in fact the land cover change occurred when the column *was*
                 ! tropical closed forest.
                 if( patch%itype(p) == nbrdlf_evr_trp_tree .and. patch%wtcol(p)  >  0._r8 )then
                    trotr1_col(c)=trotr1_col(c)+patch%wtcol(p)
                 end if
                 if( patch%itype(p) == nbrdlf_dcd_trp_tree .and. patch%wtcol(p)  >  0._r8 )then
                    trotr2_col(c)=trotr2_col(c)+patch%wtcol(p)
                 end if

                 if (transient_landcover) then
                    if( patch%itype(p) == nbrdlf_evr_trp_tree .or. patch%itype(p) == nbrdlf_dcd_trp_tree )then
                       if(dwt_smoothed(p) < 0._r8)then
                          ! Land cover change in CLM happens all at once on the first time
                          ! step of the year. However, the fire code needs deforestation
                          ! rates throughout the year, in order to combine these
                          ! deforestation rates with the current season's climate. So we
                          ! use a smoothed version of dwt.
                          !
                          ! This isn't ideal, because the carbon stocks that the fire code
                          ! is operating on will have decreased by the full annual amount
                          ! before the fire code does anything. But the biggest effect of
                          ! these deforestation fires is as a trigger for other fires, and
                          ! the C fluxes are merely diagnostic so don't need to be
                          ! conservative, so this isn't a big issue.
                          !
                          ! (Actually, it would be even better if the fire code had a
                          ! realistic breakdown of annual deforestation into the
                          ! different seasons. But having deforestation spread evenly
                          ! throughout the year is much better than having it all
                          ! concentrated on January 1.)
                          dtrotr_col(c)=dtrotr_col(c)-dwt_smoothed(p)
                       end if
                    end if
                 end if
                 if (spinup_state == 2) then         
                    rootc_col(c) = rootc_col(c) + (frootc(p) + frootc_storage(p) + &
                         frootc_xfer(p) + deadcrootc(p) * 10._r8 +       &
                         deadcrootc_storage(p) + deadcrootc_xfer(p) +    &
                         livecrootc(p)+livecrootc_storage(p) +           &
                         livecrootc_xfer(p))*patch%wtcol(p)
                 else
                    rootc_col(c) = rootc_col(c) + (frootc(p) + frootc_storage(p) + &
                         frootc_xfer(p) + deadcrootc(p) +                &
                         deadcrootc_storage(p) + deadcrootc_xfer(p) +    &
                         livecrootc(p)+livecrootc_storage(p) +           &
                         livecrootc_xfer(p))*patch%wtcol(p)
                 endif

                 fsr_col(c) = fsr_col(c) + fsr_pft(patch%itype(p))*patch%wtcol(p)/(1.0_r8-cropf_col(c))

                 hdmlf=this%forc_hdm(g)

                    ! all these constants are in Li et al. BG (2012a,b;2013)

                 if( hdmlf  >  0.1_r8 )then            
                    ! For NOT bare-soil
                    if( patch%itype(p)  /=  noveg )then
                       ! For shrub and grass (crop already excluded above)
                       if( patch%itype(p)  >=  nbrdlf_evr_shrub )then      !for shurb and grass
                           lgdp_col(c)  = lgdp_col(c) + (0.1_r8 + 0.9_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (gdp_lf(c)/8._r8)**0.5_r8))*patch%wtcol(p) &
                                  /(1.0_r8 - cropf_col(c))
                           lgdp1_col(c) = lgdp1_col(c) + (0.2_r8 + 0.8_r8*   &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (gdp_lf(c)/7._r8)))*patch%wtcol(p)/(1._r8 - cropf_col(c))
                           lpop_col(c)  = lpop_col(c) + (0.2_r8 + 0.8_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (hdmlf/450._r8)**0.5_r8))*patch%wtcol(p)/(1._r8 - cropf_col(c))
                        else   ! for trees
                           if( gdp_lf(c)  >  20._r8 )then
                              lgdp_col(c) =lgdp_col(c)+cnfire_const%occur_hi_gdp_tree*patch%wtcol(p)/(1._r8 - cropf_col(c))
                              lgdp1_col(c) =lgdp1_col(c)+0.62_r8*patch%wtcol(p)/(1._r8 - cropf_col(c)) 
                           else
                              if( gdp_lf(c) > 8._r8 )then
                                 lgdp_col(c)=lgdp_col(c)+0.79_r8*patch%wtcol(p)/(1._r8 - cropf_col(c))
                                 lgdp1_col(c)=lgdp1_col(c)+0.83_r8*patch%wtcol(p)/(1._r8 - cropf_col(c))
                              else
                                 lgdp_col(c) = lgdp_col(c)+patch%wtcol(p)/(1._r8 - cropf_col(c))
                                 lgdp1_col(c)=lgdp1_col(c)+patch%wtcol(p)/(1._r8 - cropf_col(c))
                              end if
                           end if
                           lpop_col(c) = lpop_col(c) + (0.4_r8 + 0.6_r8*    &
                                  exp(-1._r8*SHR_CONST_PI* &
                                  (hdmlf/125._r8)))*patch%wtcol(p)/(1._r8 -cropf_col(c))
                         end if
                       end if
                  else
                       lgdp_col(c)  = lgdp_col(c)+patch%wtcol(p)/(1.0_r8 - cropf_col(c))
                       lgdp1_col(c) = lgdp1_col(c)+patch%wtcol(p)/(1.0_r8 -cropf_col(c))
                       lpop_col(c)  = lpop_col(c)+patch%wtcol(p)/(1.0_r8 -cropf_col(c))
                  end if

                 fd_col(c) = fd_col(c) + fd_pft(patch%itype(p)) * patch%wtcol(p) * secsphr / (1.0_r8-cropf_col(c))         
              end if
           end if
        end do
     end do

     ! estimate annual decreased fractional coverage of BET+BDT
     ! land cover conversion in CLM4.5 is the same for each timestep except for the beginning  

     if (transient_landcover) then
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           if( dtrotr_col(c)  >  0._r8 )then
              if( kmo == 1 .and. kda == 1 .and. mcsec == 0)then
                 lfc(c) = 0._r8
              end if
              if( kmo == 1 .and. kda == 1 .and. mcsec == dt)then
                 lfc(c) = dtrotr_col(c)*dayspyr*secspday/dt
              end if
           else
              lfc(c)=0._r8
           end if
        end do
     end if
     !
     ! calculate burned area fraction in cropland
     !
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        baf_crop(c)=0._r8
     end do

     do fp = 1,num_soilp
        p = filter_soilp(fp)  
        if( kmo == 1 .and. kda == 1 .and. mcsec == 0 )then
           burndate(p) = 10000 ! init. value; actual range [0 365]
        end if
     end do

     do pi = 1,max_patch_per_col
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           g= col%gridcell(c)
           hdmlf=this%forc_hdm(g)
           if (pi <=  col%npatches(c)) then
              p = col%patchi(c) + pi - 1
              ! For crop
              if( forc_t(c)  >=  SHR_CONST_TKFRZ .and. patch%itype(p)  >  nc4_grass .and.  &
                   kmo == abm_lf(c) .and. &
                   burndate(p) >= 999 .and. patch%wtcol(p)  >  0._r8 )then ! catch  crop burn time

                 ! calculate human density impact on ag. fire
                 fhd = 0.04_r8+0.96_r8*exp(-1._r8*SHR_CONST_PI*(hdmlf/350._r8)**0.5_r8)

                 ! calculate impact of GDP on ag. fire
                 fgdp = 0.01_r8+0.99_r8*exp(-1._r8*SHR_CONST_PI*(gdp_lf(c)/10._r8))

                 ! calculate burned area
                 fb   = max(0.0_r8,min(1.0_r8,(fuelc_crop(c)-lfuel)/(ufuel-lfuel)))

                 ! crop fire only for generic crop types at this time
                 ! managed crops are treated as grasses if crop model is turned on
                 baf_crop(c) = baf_crop(c) + cropfire_a1/secsphr*fhd*fgdp*patch%wtcol(p)
                 if( fb*fhd*fgdp*patch%wtcol(p)  >  0._r8)then
                    burndate(p)=kda
                 end if
              end if
           end if
        end do
     end do
     !
     ! calculate peatland fire
     !
     do fc = 1, num_soilc
        c = filter_soilc(fc)
        g= col%gridcell(c)
        if(grc%latdeg(g) < cnfire_const%borealat )then
           baf_peatf(c) = non_boreal_peatfire_c/secsphr*max(0._r8, &
                min(1._r8,(4.0_r8-prec60_col(c)*secspday)/ &
                4.0_r8))**2*peatf_lf(c)*(1._r8-fsat(c))
        else
           baf_peatf(c) = boreal_peatfire_c/secsphr*exp(-SHR_CONST_PI*(max(wf2(c),0._r8)/0.3_r8))* &
                max(0._r8,min(1._r8,(tsoi17(c)-SHR_CONST_TKFRZ)/10._r8))*peatf_lf(c)* &
                (1._r8-fsat(c))
        end if
     end do
     !
     ! calculate other fires
     !

     ! Set the number of timesteps for e-folding.
     ! When the simulation has run fewer than this number of steps,
     ! re-scale the e-folding time to get a stable early estimate.

     ! find which pool is the cwd pool
     i_cwd = 0
     do l = 1, ndecomp_pools
        if ( is_cwd(l) ) then
           i_cwd = l
        endif
     end do

     !
     ! begin column loop to calculate fractional area affected by fire
     !
     do fc = 1, num_soilc
        c = filter_soilc(fc)
        g = col%gridcell(c)
        hdmlf=this%forc_hdm(g)
        if( cropf_col(c)  <  1._r8 )then
           fuelc(c) = totlitc(c)+totvegc(c)-rootc_col(c)-fuelc_crop(c)*cropf_col(c)
           if (spinup_state == 2) then
              fuelc(c) = fuelc(c) + ((10._r8 - 1._r8)*deadstemc_col(c))
              do j = 1, nlevdecomp  
                 fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) * dzsoi_decomp(j) * spinup_factor(i_cwd) &
                            * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              end do
           else
              do j = 1, nlevdecomp  
                 fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) * dzsoi_decomp(j)
              end do
           end if
           fuelc(c) = fuelc(c)/(1._r8-cropf_col(c))
           fb       = max(0.0_r8,min(1.0_r8,(fuelc(c)-lfuel)/(ufuel-lfuel)))
           if (trotr1_col(c)+trotr2_col(c)<=0.6_r8) then  
              afuel  =min(1._r8,max(0._r8,(fuelc(c)-2500._r8)/(5000._r8-2500._r8)))
              arh=1._r8-max(0._r8, min(1._r8,(forc_rh(g)-rh_low)/(rh_hgh-rh_low)))
              arh30=1._r8-max(0.7_r8, min(1._r8,rh30_col(c)/90._r8))
              if (forc_rh(g) < rh_hgh.and. wtlf(c) > 0._r8 .and. tsoi17(c)> SHR_CONST_TKFRZ)then
                fire_m   = ((afuel*arh30+(1._r8-afuel)*arh)**1.5_r8)*((1._r8 -max(0._r8,&
                    min(1._r8,(btran_col(c)/wtlf(c)-bt_min)/(bt_max-bt_min))))**0.5_r8)
              else
                fire_m   = 0._r8
              end if
              lh       = pot_hmn_ign_counts_alpha*6.8_r8*hdmlf**(0.43_r8)/30._r8/24._r8
              fs       = 1._r8-(0.01_r8+0.98_r8*exp(-0.025_r8*hdmlf))
              ig       = (lh+this%forc_lnfm(g)/(5.16_r8+2.16_r8*cos(3*min(60._r8,abs(grc%latdeg(g)))))*0.22_r8)  &
                         *(1._r8-fs)*(1._r8-cropf_col(c))
              nfire(c) = ig/secsphr*fb*fire_m*lgdp_col(c) !fire counts/km2/sec
              Lb_lf    = 1._r8+10._r8*(1._r8-EXP(-0.06_r8*forc_wind(g)))
              spread_m = fire_m**0.5_r8
              farea_burned(c) = min(1._r8,(cnfire_const%g0*spread_m*fsr_col(c)* &
                   fd_col(c)/1000._r8)**2*lgdp1_col(c)* &
                   lpop_col(c)*nfire(c)*SHR_CONST_PI*Lb_lf+ &
                   baf_crop(c)+baf_peatf(c))  ! fraction (0-1) per sec
           else
             farea_burned(c)=min(1._r8,baf_crop(c)+baf_peatf(c))
           end if
           !
           ! if landuse change data is used, calculate deforestation fires and 
           ! add it in the total of burned area fraction
           !
           if (transient_landcover) then
              if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 )then
                 if(( kmo == 1 .and. kda == 1 .and. mcsec == 0) .or. &
                      dtrotr_col(c) <=0._r8 )then
                    fbac1(c)        = 0._r8
                    farea_burned(c) = baf_crop(c)+baf_peatf(c)
                 else
                    cri = (4.0_r8*trotr1_col(c)+1.8_r8*trotr2_col(c))/(trotr1_col(c)+trotr2_col(c))
                    cli = (max(0._r8,min(1._r8,(cri-prec60_col(c)*secspday)/cri))**0.5)* &
                         (max(0._r8,min(1._r8,(cri-prec10_col(c)*secspday)/cri))**0.5)* &
                         max(0.0005_r8,min(1._r8,19._r8*dtrotr_col(c)*dayspyr*secspday/dt-0.001_r8))* &
                         max(0._r8,min(1._r8,(0.25_r8-(forc_rain(c)+forc_snow(c))*secsphr)/0.25_r8))
                    farea_burned(c) = cli*(cli_scale/secspday)+baf_crop(c)+baf_peatf(c)
                    ! burned area out of conversion region due to land use fire
                    fbac1(c) = max(0._r8,fb*cli*(cli_scale/secspday) - 2.0_r8*lfc(c)/dt)   
                 end if
                 ! total burned area out of conversion 
                 fbac(c) = fbac1(c)+baf_crop(c)+baf_peatf(c) 
              else
                 fbac(c) = farea_burned(c)
              end if
           end if

        else
           farea_burned(c) = min(1._r8,baf_crop(c)+baf_peatf(c))
        end if

     end do  ! end of column loop

   end associate

 end subroutine CNFireArea

end module CNFireLi2016Mod
