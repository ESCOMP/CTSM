!-----------------------------------------------------------------------
subroutine CNiniTimeVar(bounds)

  ! !DESCRIPTION:
  ! Initializes time varying variables used only in
  ! coupled carbon-nitrogen mode (CN):
  !
  ! !USES:
  use clmtype
  use clm_atmlnd  , only: clm_a2l
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varcon  , only: istsoil, zsoi
  use clm_varpar  , only: nlevdecomp, ndecomp_pools, nlevdecomp_full, crop_prog
  use clm_varcon  , only: istcrop, c13ratio, c14ratio
  use clm_varctl  , only: use_c13, use_c14, use_nitrif_denitrif, use_cndv
  use pftvarcon   , only: noveg
  use pftvarcon   , only: npcropmin
  use decompMod   , only: bounds_type
  !
  ! !ARGUMENTS:
  implicit none
  type(bounds_type), intent(in) :: bounds  ! bounds
  !
  ! !LOCAL VARIABLES:
  integer :: g,l,c,p,j,k  ! indices
  !-----------------------------------------------------------------------

   associate(& 
   annsum_counter                      =>    cps%annsum_counter                          , & ! Output: [real(r8) (:)]  seconds since last annual accumulator turnover    
   cannsum_npp                         =>    cps%cannsum_npp                             , & ! Output: [real(r8) (:)]  annual sum of NPP, averaged from pft-level (gC/m2/yr)
   cannavg_t2m                         =>    cps%cannavg_t2m                             , & ! Output: [real(r8) (:)] annual average of 2m air temperature, averaged from pft-level (K)
   wf                                  =>    cps%wf                                      , & ! Output: [real(r8) (:)]  soil moisture in top 0.05 m                       
   wf2                                 =>    cps%wf2                                     , & ! Output: [real(r8) (:)]                                                    
   nfire                               =>    cps%nfire                                   , & ! Output: [real(r8) (:)]  fire counts/km2/timestep                          
   baf_crop                            =>    cps%baf_crop                                , & ! Output: [real(r8) (:)]  burned area fraction in crop                      
   baf_peatf                           =>    cps%baf_peatf                               , & ! Output: [real(r8) (:)]  burned area fraction in peatland                  
   fbac                                =>    cps%fbac                                    , & ! Output: [real(r8) (:)]                                                    
   fbac1                               =>    cps%fbac1                                   , & ! Output: [real(r8) (:)]                                                    
   farea_burned                        =>    cps%farea_burned                            , & ! Output: [real(r8) (:)]  timestep fractional area burned (proportion)      
   qflx_drain                          =>    cwf%qflx_drain                              , & ! Output: [real(r8) (:)]  sub-surface runoff (mm H2O /s)                    
   qflx_surf                           =>    cwf%qflx_surf                               , & ! Output: [real(r8) (:)]  surface runoff (mm H2O /s)                        
   decomp_cpools                       =>    ccs%decomp_cpools                           , & ! Output: [real(r8) (:,:)]  (gC/m2)  decomposing (litter, cwd, soil) c pools
   decomp_cpools_1m                    =>    ccs%decomp_cpools_1m                        , & ! Output: [real(r8) (:,:)]  (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
   decomp_cpools_vr                    =>    ccs%decomp_cpools_vr                        , & ! Output: [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   decomp_npools                       =>    cns%decomp_npools                           , & ! Output: [real(r8) (:,:)]  (gC/m2)  decomposing (litter, cwd, soil) N pools
   decomp_npools_vr                    =>    cns%decomp_npools_vr                        , & ! Output: [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   decomp_npools_1m                    =>    cns%decomp_npools_1m                        , & ! Output: [real(r8) (:,:)]  (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
   nfixation_prof                      =>    cps%nfixation_prof                          , & ! Output: [real(r8) (:,:)]  (1/m) profile for N fixation additions          
   ndep_prof                           =>    cps%ndep_prof                               , & ! Output: [real(r8) (:,:)]  (1/m) profile for N fixation additions          
   seedc                               =>    ccs%seedc                                   , & ! Output: [real(r8) (:)]  (gC/m2) column-level pool for seeding new PFTs    
   prod10c                             =>    ccs%prod10c                                 , & ! Output: [real(r8) (:)]  (gC/m2) wood product C pool, 10-year lifespan     
   prod100c                            =>    ccs%prod100c                                , & ! Output: [real(r8) (:)]  (gC/m2) wood product C pool, 100-year lifespan    
   totprodc                            =>    ccs%totprodc                                , & ! Output: [real(r8) (:)]  (gC/m2) total wood product C                      
   seedn                               =>    cns%seedn                                   , & ! Output: [real(r8) (:)]  (gN/m2) column-level pool for seeding new PFTs    
   prod10n                             =>    cns%prod10n                                 , & ! Output: [real(r8) (:)]  (gN/m2) wood product N pool, 10-year lifespan     
   prod100n                            =>    cns%prod100n                                , & ! Output: [real(r8) (:)]  (gN/m2) wood product N pool, 100-year lifespan    
   totprodn                            =>    cns%totprodn                                , & ! Output: [real(r8) (:)]  (gN/m2) total wood product N                      
   sminn                               =>    cns%sminn                                   , & ! Output: [real(r8) (:)]  (gN/m2) soil mineral N                            
   col_ctrunc                          =>    ccs%col_ctrunc                              , & ! Output: [real(r8) (:)]  (gC/m2) column-level sink for C truncation (diagnostic)
   sminn_vr                            =>    cns%sminn_vr                                , & ! Output: [real(r8) (:,:)]  (gN/m3) soil mineral N                          
   col_ctrunc_vr                       =>    ccs%col_ctrunc_vr                           , & ! Output: [real(r8) (:,:)]  (gC/m3) column-level sink for C truncation (prognostic)
   col_ntrunc_vr                       =>    cns%col_ntrunc_vr                           , & ! Output: [real(r8) (:,:)]  (gN/m3) column-level sink for N truncation      
   fpi_vr                              =>    cps%fpi_vr                                  , & ! Output: [real(r8) (:,:)]                                                  
   alt                                 =>    cps%alt                                     , & ! Output: [real(r8) (:)]                                                    
   altmax                              =>    cps%altmax                                  , & ! Output: [real(r8) (:)]                                                    
   altmax_lastyear                     =>    cps%altmax_lastyear                         , & ! Output: [real(r8) (:)]                                                    
   som_adv_coef                        =>    cps%som_adv_coef                            , & ! Output: [real(r8) (:,:)]                                                  
   som_diffus_coef                     =>    cps%som_diffus_coef                         , & ! Output: [real(r8) (:,:)]                                                  
   alt_indx                            =>    cps%alt_indx                                , & ! Output: [integer (:)]                                                     
   altmax_indx                         =>    cps%altmax_indx                             , & ! Output: [integer (:)]                                                     
   altmax_lastyear_indx                =>    cps%altmax_lastyear_indx                    , & ! Output: [integer (:)]                                                     
   smin_nh4_vr                         =>    cns%smin_nh4_vr                             , & ! Output: [real(r8) (:,:)]  (gN/m3) soil mineral NH4 pool                   
   smin_no3_vr                         =>    cns%smin_no3_vr                             , & ! Output: [real(r8) (:,:)]  (gN/m3) soil mineral NO3 pool                   
   smin_nh4                            =>    cns%smin_nh4                                , & ! Output: [real(r8) (:)]  (gN/m2) soil mineral NH4 pool                     
   smin_no3                            =>    cns%smin_no3                                , & ! Output: [real(r8) (:)]  (gN/m2) soil mineral NO3 pool                     
   totcolc                             =>    ccs%totcolc                                 , & ! Output: [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool   
   cwdc                                =>    ccs%cwdc                                    , & ! Output: [real(r8) (:)]  (gC/m2) coarse woody debris C                     
   totecosysc                          =>    ccs%totecosysc                              , & ! Output: [real(r8) (:)]  (gC/m2) total ecosystem carbon, incl veg but excl cpool
   totlitc                             =>    ccs%totlitc                                 , & ! Output: [real(r8) (:)]  (gC/m2) total litter carbon                       
   totsomc                             =>    ccs%totsomc                                 , & ! Output: [real(r8) (:)]  (gC/m2) total soil organic matter carbon          
   totlitc_1m                          =>    ccs%totlitc_1m                              , & ! Output: [real(r8) (:)]  (gC/m2) total litter carbon to 1 meter            
   totsomc_1m                          =>    ccs%totsomc_1m                              , & ! Output: [real(r8) (:)]  (gC/m2) total soil organic matter carbon to 1 meter
   totcoln                             =>    cns%totcoln                                 , & ! Output: [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg           
   cwdn                                =>    cns%cwdn                                    , & ! Output: [real(r8) (:)]  (gN/m2) coarse woody debris N                     
   totecosysn                          =>    cns%totecosysn                              , & ! Output: [real(r8) (:)]  (gN/m2) total ecosystem nitrogen, incl veg        
   totlitn                             =>    cns%totlitn                                 , & ! Output: [real(r8) (:)]  (gN/m2) total litter nitrogen                     
   totsomn                             =>    cns%totsomn                                 , & ! Output: [real(r8) (:)]  (gN/m2) total soil organic matter nitrogen        
   totlitn_1m                          =>    cns%totlitn_1m                              , & ! Output: [real(r8) (:)]  (gN/m2) total litter nitrogen to 1 meter          
   totsomn_1m                          =>    cns%totsomn_1m                              , & ! Output: [real(r8) (:)]  (gN/m2) total soil organic matter nitrogen to 1 meter
   seedc13                             =>    cc13s%seedc                                 , & ! Output: [real(r8) (:)]  (gC/m2) column-level pool for seeding new PFTs    
   prod10c13                           =>    cc13s%prod10c                               , & ! Output: [real(r8) (:)]  (gC/m2) wood product C13 pool, 10-year lifespan   
   prod100c13                          =>    cc13s%prod100c                              , & ! Output: [real(r8) (:)]  (gC/m2) wood product C13 pool, 100-year lifespan  
   totprodc13                          =>    cc13s%totprodc                              , & ! Output: [real(r8) (:)]  (gC/m2) total wood product C13                    
   cwdc13                              =>    cc13s%cwdc                                  , & ! Output: [real(r8) (:)]  (gC/m2) coarse woody debris C                     
   decomp_c13pools                     =>    cc13s%decomp_cpools                         , & ! Output: [real(r8) (:,:)]  (gC/m2)  decomposing (litter, cwd, soil) c pools
   decomp_c13pools_vr                  =>    cc13s%decomp_cpools_vr                      , & ! Output: [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   c13_col_ctrunc_vr                   =>    cc13s%col_ctrunc_vr                         , & ! Output: [real(r8) (:,:)]  (gC/m3) C truncation term                       
   decomp_c13pools_1m                  =>    cc13s%decomp_cpools_1m                      , & ! Output: [real(r8) (:,:)]  (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
   c13_psnsun                          =>    pc13f%psnsun                                , & ! Output: [real(r8) (:)]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)    
   c13_psnsha                          =>    pc13f%psnsha                                , & ! Output: [real(r8) (:)]  shaded leaf photosynthesis (umol CO2 /m**2/ s)    
   xsmrpool_c13ratio                   =>    pepv%xsmrpool_c13ratio                      , & ! Output: [real(r8) (:)]  C flux assigned to recovery of negative cpool (gC/m2/s)
   alphapsnsun                         =>    pps%alphapsnsun                             , & ! Output: [real(r8) (:)] sunlit 13c fractionation ([])                      
   alphapsnsha                         =>    pps%alphapsnsha                             , & ! Output: [real(r8) (:)] shaded 13c fractionation ([])                      
   leafc13                             =>    pc13s%leafc                                 , & ! Output: [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc13_storage                     =>    pc13s%leafc_storage                         , & ! Output: [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc13_xfer                        =>    pc13s%leafc_xfer                            , & ! Output: [real(r8) (:)]  (gC/m2) leaf C transfer                           
   frootc13                            =>    pc13s%frootc                                , & ! Output: [real(r8) (:)]  (gC/m2) fine root C                               
   frootc13_storage                    =>    pc13s%frootc_storage                        , & ! Output: [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc13_xfer                       =>    pc13s%frootc_xfer                           , & ! Output: [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livestemc13                         =>    pc13s%livestemc                             , & ! Output: [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc13_storage                 =>    pc13s%livestemc_storage                     , & ! Output: [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc13_xfer                    =>    pc13s%livestemc_xfer                        , & ! Output: [real(r8) (:)]  (gC/m2) live stem C transfer                      
   deadstemc13                         =>    pc13s%deadstemc                             , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc13_storage                 =>    pc13s%deadstemc_storage                     , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc13_xfer                    =>    pc13s%deadstemc_xfer                        , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   livecrootc13                        =>    pc13s%livecrootc                            , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc13_storage                =>    pc13s%livecrootc_storage                    , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc13_xfer                   =>    pc13s%livecrootc_xfer                       , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   deadcrootc13                        =>    pc13s%deadcrootc                            , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc13_storage                =>    pc13s%deadcrootc_storage                    , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc13_xfer                   =>    pc13s%deadcrootc_xfer                       , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   c13_gresp_storage                   =>    pc13s%gresp_storage                         , & ! Output: [real(r8) (:)]  (gC/m2) growth respiration storage                
   c13_gresp_xfer                      =>    pc13s%gresp_xfer                            , & ! Output: [real(r8) (:)]  (gC/m2) growth respiration transfer               
   c13pool                             =>    pc13s%cpool                                 , & ! Output: [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   c13xsmrpool                         =>    pc13s%xsmrpool                              , & ! Output: [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   c13_pft_ctrunc                      =>    pc13s%pft_ctrunc                            , & ! Output: [real(r8) (:)]  (gC/m2) C truncation term                         
   totvegc13                           =>    pc13s%totvegc                               , & ! Output: [real(r8) (:)]  (gC/m2) total vegetation carbon, excluding cpool  
   rc13_canair                         =>    pepv%rc13_canair                            , & ! Output: [real(r8) (:)] C13O2/C12O2 in canopy air                          
   rc13_psnsun                         =>    pepv%rc13_psnsun                            , & ! Output: [real(r8) (:)] C13O2/C12O2 in sunlit canopy psn flux              
   rc13_psnsha                         =>    pepv%rc13_psnsha                            , & ! Output: [real(r8) (:)] C13O2/C12O2 in shaded canopy psn flux              
   seedc14                             =>    cc14s%seedc                                 , & ! Output: [real(r8) (:)]  (gC/m2) column-level pool for seeding new PFTs    
   prod10c14                           =>    cc14s%prod10c                               , & ! Output: [real(r8) (:)]  (gC/m2) wood product C14 pool, 10-year lifespan   
   prod100c14                          =>    cc14s%prod100c                              , & ! Output: [real(r8) (:)]  (gC/m2) wood product C14 pool, 100-year lifespan  
   totprodc14                          =>    cc14s%totprodc                              , & ! Output: [real(r8) (:)]  (gC/m2) total wood product C14                    
   cwdc14                              =>    cc14s%cwdc                                  , & ! Output: [real(r8) (:)]  (gC/m2) coarse woody debris C                     
   decomp_c14pools                     =>    cc14s%decomp_cpools                         , & ! Output: [real(r8) (:,:)]  (gC/m2)  decomposing (litter, cwd, soil) c pools
   decomp_c14pools_vr                  =>    cc14s%decomp_cpools_vr                      , & ! Output: [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   c14_col_ctrunc_vr                   =>    cc14s%col_ctrunc_vr                         , & ! Output: [real(r8) (:,:)]  (gC/m3) C truncation term                       
   decomp_c14pools_1m                  =>    cc14s%decomp_cpools_1m                      , & ! Output: [real(r8) (:,:)]  (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
   c14_psnsun                          =>    pc14f%psnsun                                , & ! Output: [real(r8) (:)]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)    
   c14_psnsha                          =>    pc14f%psnsha                                , & ! Output: [real(r8) (:)]  shaded leaf photosynthesis (umol CO2 /m**2/ s)    
   leafc14                             =>    pc14s%leafc                                 , & ! Output: [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc14_storage                     =>    pc14s%leafc_storage                         , & ! Output: [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc14_xfer                        =>    pc14s%leafc_xfer                            , & ! Output: [real(r8) (:)]  (gC/m2) leaf C transfer                           
   frootc14                            =>    pc14s%frootc                                , & ! Output: [real(r8) (:)]  (gC/m2) fine root C                               
   frootc14_storage                    =>    pc14s%frootc_storage                        , & ! Output: [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc14_xfer                       =>    pc14s%frootc_xfer                           , & ! Output: [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livestemc14                         =>    pc14s%livestemc                             , & ! Output: [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc14_storage                 =>    pc14s%livestemc_storage                     , & ! Output: [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc14_xfer                    =>    pc14s%livestemc_xfer                        , & ! Output: [real(r8) (:)]  (gC/m2) live stem C transfer                      
   deadstemc14                         =>    pc14s%deadstemc                             , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc14_storage                 =>    pc14s%deadstemc_storage                     , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc14_xfer                    =>    pc14s%deadstemc_xfer                        , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   livecrootc14                        =>    pc14s%livecrootc                            , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc14_storage                =>    pc14s%livecrootc_storage                    , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc14_xfer                   =>    pc14s%livecrootc_xfer                       , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   deadcrootc14                        =>    pc14s%deadcrootc                            , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc14_storage                =>    pc14s%deadcrootc_storage                    , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc14_xfer                   =>    pc14s%deadcrootc_xfer                       , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   c14_gresp_storage                   =>    pc14s%gresp_storage                         , & ! Output: [real(r8) (:)]  (gC/m2) growth respiration storage                
   c14_gresp_xfer                      =>    pc14s%gresp_xfer                            , & ! Output: [real(r8) (:)]  (gC/m2) growth respiration transfer               
   c14pool                             =>    pc14s%cpool                                 , & ! Output: [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   c14xsmrpool                         =>    pc14s%xsmrpool                              , & ! Output: [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   c14_pft_ctrunc                      =>    pc14s%pft_ctrunc                            , & ! Output: [real(r8) (:)]  (gC/m2) C truncation term                         
   totvegc14                           =>    pc14s%totvegc                               , & ! Output: [real(r8) (:)]  (gC/m2) total vegetation carbon, excluding cpool  
   rc14_atm                            =>    pepv%rc14_atm                               , & ! Output: [real(r8) (:)] C14O2/C12O2 in atmosphere                          
   soyfixn                             =>    pnf%soyfixn                                 , & ! Output: [real(r8) (:)]                                                    
   fert                                =>    pnf%fert                                    , & ! Output: [real(r8) (:)]                                                    
   fert_counter                        =>    pepv%fert_counter                           , & ! Output: [real(r8) (:)]                                                    
   grain_flag                          =>    pepv%grain_flag                             , & ! Output: [real(r8) (:)]                                                    
   leafc                               =>    pcs%leafc                                   , & ! Output: [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc_storage                       =>    pcs%leafc_storage                           , & ! Output: [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc_xfer                          =>    pcs%leafc_xfer                              , & ! Output: [real(r8) (:)]  (gC/m2) leaf C transfer                           
   grainc                              =>    pcs%grainc                                  , & ! Output: [real(r8) (:)]  (gC/m2) grain C                                   
   grainc_storage                      =>    pcs%grainc_storage                          , & ! Output: [real(r8) (:)]  (gC/m2) grain C storage                           
   grainc_xfer                         =>    pcs%grainc_xfer                             , & ! Output: [real(r8) (:)]  (gC/m2) grain C transfer                          
   frootc                              =>    pcs%frootc                                  , & ! Output: [real(r8) (:)]  (gC/m2) fine root C                               
   frootc_storage                      =>    pcs%frootc_storage                          , & ! Output: [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc_xfer                         =>    pcs%frootc_xfer                             , & ! Output: [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livestemc                           =>    pcs%livestemc                               , & ! Output: [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc_storage                   =>    pcs%livestemc_storage                       , & ! Output: [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc_xfer                      =>    pcs%livestemc_xfer                          , & ! Output: [real(r8) (:)]  (gC/m2) live stem C transfer                      
   deadstemc                           =>    pcs%deadstemc                               , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc_storage                   =>    pcs%deadstemc_storage                       , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc_xfer                      =>    pcs%deadstemc_xfer                          , & ! Output: [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   livecrootc                          =>    pcs%livecrootc                              , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc_storage                  =>    pcs%livecrootc_storage                      , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc_xfer                     =>    pcs%livecrootc_xfer                         , & ! Output: [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   deadcrootc                          =>    pcs%deadcrootc                              , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc_storage                  =>    pcs%deadcrootc_storage                      , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc_xfer                     =>    pcs%deadcrootc_xfer                         , & ! Output: [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   gresp_storage                       =>    pcs%gresp_storage                           , & ! Output: [real(r8) (:)]  (gC/m2) growth respiration storage                
   gresp_xfer                          =>    pcs%gresp_xfer                              , & ! Output: [real(r8) (:)]  (gC/m2) growth respiration transfer               
   cpool                               =>    pcs%cpool                                   , & ! Output: [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   xsmrpool                            =>    pcs%xsmrpool                                , & ! Output: [real(r8) (:)]  (gC/m2) abstract C pool to meet excess MR demand  
   forc_hgt_u_pft                      =>    pps%forc_hgt_u_pft                          , & ! Output: [real(r8) (:)] observational height of wind at pft-level [m]      
   woodc                               =>    pcs%woodc                                   , & ! Output: [real(r8) (:)]  (gC/m2) pft-level wood C                          
   leafn                               =>    pns%leafn                                   , & ! Output: [real(r8) (:)]  (gN/m2) leaf N                                    
   leafn_storage                       =>    pns%leafn_storage                           , & ! Output: [real(r8) (:)]  (gN/m2) leaf N storage                            
   leafn_xfer                          =>    pns%leafn_xfer                              , & ! Output: [real(r8) (:)]  (gN/m2) leaf N transfer                           
   grainn                              =>    pns%grainn                                  , & ! Output: [real(r8) (:)]  (gN/m2) grain N                                   
   grainn_storage                      =>    pns%grainn_storage                          , & ! Output: [real(r8) (:)]  (gN/m2) grain N storage                           
   grainn_xfer                         =>    pns%grainn_xfer                             , & ! Output: [real(r8) (:)]  (gN/m2) grain N transfer                          
   frootn                              =>    pns%frootn                                  , & ! Output: [real(r8) (:)]  (gN/m2) fine root N                               
   frootn_storage                      =>    pns%frootn_storage                          , & ! Output: [real(r8) (:)]  (gN/m2) fine root N storage                       
   frootn_xfer                         =>    pns%frootn_xfer                             , & ! Output: [real(r8) (:)]  (gN/m2) fine root N transfer                      
   livestemn                           =>    pns%livestemn                               , & ! Output: [real(r8) (:)]  (gN/m2) live stem N                               
   livestemn_storage                   =>    pns%livestemn_storage                       , & ! Output: [real(r8) (:)]  (gN/m2) live stem N storage                       
   livestemn_xfer                      =>    pns%livestemn_xfer                          , & ! Output: [real(r8) (:)]  (gN/m2) live stem N transfer                      
   deadstemn                           =>    pns%deadstemn                               , & ! Output: [real(r8) (:)]  (gN/m2) dead stem N                               
   deadstemn_storage                   =>    pns%deadstemn_storage                       , & ! Output: [real(r8) (:)]  (gN/m2) dead stem N storage                       
   deadstemn_xfer                      =>    pns%deadstemn_xfer                          , & ! Output: [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   livecrootn                          =>    pns%livecrootn                              , & ! Output: [real(r8) (:)]  (gN/m2) live coarse root N                        
   livecrootn_storage                  =>    pns%livecrootn_storage                      , & ! Output: [real(r8) (:)]  (gN/m2) live coarse root N storage                
   livecrootn_xfer                     =>    pns%livecrootn_xfer                         , & ! Output: [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   deadcrootn                          =>    pns%deadcrootn                              , & ! Output: [real(r8) (:)]  (gN/m2) dead coarse root N                        
   deadcrootn_storage                  =>    pns%deadcrootn_storage                      , & ! Output: [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   deadcrootn_xfer                     =>    pns%deadcrootn_xfer                         , & ! Output: [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   retransn                            =>    pns%retransn                                , & ! Output: [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   npool                               =>    pns%npool                                   , & ! Output: [real(r8) (:)]  (gN/m2) temporary plant N pool                    
   psnsun                              =>    pcf%psnsun                                  , & ! Output: [real(r8) (:)]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)    
   psnsha                              =>    pcf%psnsha                                  , & ! Output: [real(r8) (:)]  shaded leaf photosynthesis (umol CO2 /m**2/ s)    
   laisun                              =>    pps%laisun                                  , & ! Output: [real(r8) (:)]  sunlit projected leaf area index                  
   laisha                              =>    pps%laisha                                  , & ! Output: [real(r8) (:)]  shaded projected leaf area index                  
   dormant_flag                        =>    pepv%dormant_flag                           , & ! Output: [real(r8) (:)]  dormancy flag                                     
   days_active                         =>    pepv%days_active                            , & ! Output: [real(r8) (:)]  number of days since last dormancy                
   onset_flag                          =>    pepv%onset_flag                             , & ! Output: [real(r8) (:)]  onset flag                                        
   onset_counter                       =>    pepv%onset_counter                          , & ! Output: [real(r8) (:)]  onset days counter                                
   onset_gddflag                       =>    pepv%onset_gddflag                          , & ! Output: [real(r8) (:)]  onset flag for growing degree day sum             
   onset_fdd                           =>    pepv%onset_fdd                              , & ! Output: [real(r8) (:)]  onset freezing degree days counter                
   onset_gdd                           =>    pepv%onset_gdd                              , & ! Output: [real(r8) (:)]  onset growing degree days                         
   onset_swi                           =>    pepv%onset_swi                              , & ! Output: [real(r8) (:)]  onset soil water index                            
   offset_flag                         =>    pepv%offset_flag                            , & ! Output: [real(r8) (:)]  offset flag                                       
   offset_counter                      =>    pepv%offset_counter                         , & ! Output: [real(r8) (:)]  offset days counter                               
   offset_fdd                          =>    pepv%offset_fdd                             , & ! Output: [real(r8) (:)]  offset freezing degree days counter               
   offset_swi                          =>    pepv%offset_swi                             , & ! Output: [real(r8) (:)]  offset soil water index                           
   lgsf                                =>    pepv%lgsf                                   , & ! Output: [real(r8) (:)]  long growing season factor [0-1]                  
   bglfr                               =>    pepv%bglfr                                  , & ! Output: [real(r8) (:)]  background litterfall rate (1/s)                  
   bgtr                                =>    pepv%bgtr                                   , & ! Output: [real(r8) (:)]  background transfer rate (1/s)                    
   annavg_t2m                          =>    pepv%annavg_t2m                             , & ! Output: [real(r8) (:)]  annual average 2m air temperature (K)             
   tempavg_t2m                         =>    pepv%tempavg_t2m                            , & ! Output: [real(r8) (:)]  temporary average 2m air temperature (K)          
   gpp                                 =>    pepv%gpp                                    , & ! Output: [real(r8) (:)]  GPP flux before downregulation (gC/m2/s)          
   availc                              =>    pepv%availc                                 , & ! Output: [real(r8) (:)]  C flux available for allocation (gC/m2/s)         
   xsmrpool_recover                    =>    pepv%xsmrpool_recover                       , & ! Output: [real(r8) (:)]  C flux assigned to recovery of negative cpool (gC/m2/s)
   alloc_pnow                          =>    pepv%alloc_pnow                             , & ! Output: [real(r8) (:)]  fraction of current allocation to display as new growth (DIM)
   c_allometry                         =>    pepv%c_allometry                            , & ! Output: [real(r8) (:)]  C allocation index (DIM)                          
   n_allometry                         =>    pepv%n_allometry                            , & ! Output: [real(r8) (:)]  N allocation index (DIM)                          
   plant_ndemand                       =>    pepv%plant_ndemand                          , & ! Output: [real(r8) (:)]  N flux required to support initial GPP (gN/m2/s)  
   tempsum_potential_gpp               =>    pepv%tempsum_potential_gpp                  , & ! Output: [real(r8) (:)]  temporary annual sum of plant_ndemand             
   annsum_potential_gpp                =>    pepv%annsum_potential_gpp                   , & ! Output: [real(r8) (:)]  annual sum of plant_ndemand                       
   tempmax_retransn                    =>    pepv%tempmax_retransn                       , & ! Output: [real(r8) (:)]  temporary max of retranslocated N pool (gN/m2)    
   annmax_retransn                     =>    pepv%annmax_retransn                        , & ! Output: [real(r8) (:)]  annual max of retranslocated N pool (gN/m2)       
   avail_retransn                      =>    pepv%avail_retransn                         , & ! Output: [real(r8) (:)]  N flux available from retranslocation pool (gN/m2/s)
   plant_nalloc                        =>    pepv%plant_nalloc                           , & ! Output: [real(r8) (:)]  total allocated N flux (gN/m2/s)                  
   plant_calloc                        =>    pepv%plant_calloc                           , & ! Output: [real(r8) (:)]  total allocated C flux (gC/m2/s)                  
   excess_cflux                        =>    pepv%excess_cflux                           , & ! Output: [real(r8) (:)]  C flux not allocated due to downregulation (gC/m2/s)
   downreg                             =>    pepv%downreg                                , & ! Output: [real(r8) (:)]  fractional reduction in GPP due to N limitation (DIM)
   tempsum_npp                         =>    pepv%tempsum_npp                            , & ! Output: [real(r8) (:)]  temporary annual sum of NPP                       
   annsum_npp                          =>    pepv%annsum_npp                             , & ! Output: [real(r8) (:)]  annual sum of NPP                                 
   tempsum_litfall                     =>    pepv%tempsum_litfall                        , & ! Output: [real(r8) (:)]  temporary annual sum of litfall                   
   annsum_litfall                      =>    pepv%annsum_litfall                         , & ! Output: [real(r8) (:)]  annual sum of litfall                             
   dispvegc                            =>    pcs%dispvegc                                , & ! Output: [real(r8) (:)]  (gC/m2) displayed veg carbon, excluding storage and cpool
   pft_ctrunc                          =>    pcs%pft_ctrunc                              , & ! Output: [real(r8) (:)]  (gC/m2) pft-level sink for C truncation           
   storvegc                            =>    pcs%storvegc                                , & ! Output: [real(r8) (:)]  (gC/m2) stored vegetation carbon, excluding cpool 
   totpftc                             =>    pcs%totpftc                                 , & ! Output: [real(r8) (:)]  (gC/m2) total pft-level carbon, including cpool   
   totvegc                             =>    pcs%totvegc                                 , & ! Output: [real(r8) (:)]  (gC/m2) total vegetation carbon, excluding cpool  
   prev_frootc_to_litter               =>    pepv%prev_frootc_to_litter                  , & ! Output: [real(r8) (:)] previous timestep froot C litterfall flux (gC/m2/s)
   prev_leafc_to_litter                =>    pepv%prev_leafc_to_litter                   , & ! Output: [real(r8) (:)] previous timestep leaf C litterfall flux (gC/m2/s) 
   dispvegn                            =>    pns%dispvegn                                , & ! Output: [real(r8) (:)]  (gN/m2) displayed veg nitrogen, excluding storage 
   pft_ntrunc                          =>    pns%pft_ntrunc                              , & ! Output: [real(r8) (:)]  (gN/m2) pft-level sink for N truncation           
   storvegn                            =>    pns%storvegn                                , & ! Output: [real(r8) (:)]  (gN/m2) stored vegetation nitrogen                
   totpftn                             =>    pns%totpftn                                 , & ! Output: [real(r8) (:)]  (gN/m2) total pft-level nitrogen                  
   totvegn                             =>    pns%totvegn                                 , & ! Output: [real(r8) (:)]  (gN/m2) total vegetation nitrogen                 
   evergreen                           =>    pftcon%evergreen                            , & ! Input:  [real(r8) (:)]  binary flag for evergreen leaf habit (0 or 1)     
   woody                               =>    pftcon%woody                                , & ! Input:  [real(r8) (:)]  binary flag for woody lifeform (1=woody, 0=not woody)
   leafcn                              =>    pftcon%leafcn                               , & ! Input:  [real(r8) (:)]  leaf C:N (gC/gN)                                  
   deadwdcn                            =>    pftcon%deadwdcn                             , & ! Input:  [real(r8) (:)]  dead wood (xylem and heartwood) C:N (gC/gN)       
   initial_cn_ratio                    =>    decomp_cascade_con%initial_cn_ratio         , & ! Output: [real(r8) (:)]  c:n ratio for initialization of pools             
   initial_stock                       =>    decomp_cascade_con%initial_stock              & ! Output: [real(r8) (:)]  initial concentration for seeding at spinup       
   )

   ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
   ! since this is not initialized before first call to CNVegStructUpdate,
   ! and it is required to set the upper bound for canopy top height.
   ! Changed 3/21/08, KO: still needed but don't have sufficient information 
   ! to set this properly (e.g., pft-level displacement height and roughness 
   ! length). So leave at 30m.
   do p = bounds%begp,bounds%endp
      forc_hgt_u_pft(p) = 30._r8
   end do

   ! initialize column-level variables
   do c = bounds%begc, bounds%endc
      l = col%landunit(c)
      if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

         ! column physical state variables
         annsum_counter(c) = 0._r8
         cannsum_npp(c)    = 0._r8
         cannavg_t2m(c)    = 280._r8
         ! fire related variables changed by F. Li and S. Levis
         wf(c) = 1.0_r8  ! it needs to be non zero so the first time step has no fires
         wf2(c) = 1.0_r8
         nfire(c) = 0._r8
         baf_crop(c) = 0._r8
         baf_peatf(c) = 0._r8
         fbac(c) = 0._r8
         fbac1(c) = 0._r8
         farea_burned(c) = 0._r8

         
         ! needed for CNNLeaching
         qflx_drain(c) = 0._r8
         qflx_surf(c) = 0._r8

         ! column carbon state variable initialization
         do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
               if (zsoi(j) .lt. 0.3 ) then  !! only initialize upper soil column
                  decomp_cpools_vr(c,j,k) = initial_stock(k)
               else
                  decomp_cpools_vr(c,j,k) = 0._r8
               endif
            end do
            col_ctrunc_vr(c,j) = 0._r8
         end do
         if ( nlevdecomp .gt. 1 ) then
            do j = nlevdecomp+1, nlevdecomp_full
               do k = 1, ndecomp_pools
                  decomp_cpools_vr(c,j,k) = 0._r8
               end do
               col_ctrunc_vr(c,j) = 0._r8
            end do
         end if
         do k = 1, ndecomp_pools
            decomp_cpools(c,k) = initial_stock(k)
            decomp_cpools_1m(c,k) = initial_stock(k)
         end do

         do j = 1, nlevdecomp_full
            ! initialize fpi_vr so that levels below nlevsoi are not nans
            fpi_vr(c,j) = 0._r8
            som_adv_coef(c,j) = 0._r8
            som_diffus_coef(c,j) = 0._r8
            ! here initialize the profiles for converting to vertically resolved carbon pools
            nfixation_prof(c,j) = 0._r8
            ndep_prof(c,j) = 0._r8
         end do

         ! and define alt variables to be zero
         alt(c) = 0._r8
         altmax(c) = 0._r8
         altmax_lastyear(c) = 0._r8
         alt_indx(c) = 0
         altmax_indx(c) = 0
         altmax_lastyear_indx = 0

         cwdc(c)       = 0._r8
         col_ctrunc(c) = 0._r8
         totlitc(c)    = 0._r8
         totsomc(c)    = 0._r8
         totlitc_1m(c) = 0._r8
         totsomc_1m(c) = 0._r8
         totecosysc(c) = 0._r8
         totcolc(c)    = 0._r8

         if ( use_c13 ) then
            do j = 1, nlevdecomp
               do k = 1, ndecomp_pools
                  decomp_c13pools_vr(c,j,k) = decomp_cpools_vr(c,j,k) * c13ratio
               end do
               c13_col_ctrunc_vr(c,j) = col_ctrunc_vr(c,j) * c13ratio
            end do
            if ( nlevdecomp .gt. 1 ) then
               do j = nlevdecomp+1, nlevdecomp_full
                  do k = 1, ndecomp_pools
                     decomp_c13pools_vr(c,j,k) = 0._r8
                  end do
                  c13_col_ctrunc_vr(c,j) = 0._r8
               end do
            end if
            cwdc13(c) = cwdc(c) * c13ratio
            do k = 1, ndecomp_pools
               decomp_c13pools(c,k) = decomp_cpools(c,k) * c13ratio
               decomp_c13pools_1m(c,k) = decomp_cpools_1m(c,k) * c13ratio
            end do
         endif

         if ( use_c14 ) then
            do j = 1, nlevdecomp
               do k = 1, ndecomp_pools
                  decomp_c14pools_vr(c,j,k) = decomp_cpools_vr(c,j,k) * c14ratio
               end do
               c14_col_ctrunc_vr(c,j) = col_ctrunc_vr(c,j) * c14ratio
            end do
            if ( nlevdecomp .gt. 1 ) then
               do j = nlevdecomp+1, nlevdecomp_full
                  do k = 1, ndecomp_pools
                     decomp_c14pools_vr(c,j,k) = 0._r8
                  end do
                  c14_col_ctrunc_vr(c,j) = 0._r8
               end do
            end if
            cwdc14(c) = cwdc(c) * c14ratio
            do k = 1, ndecomp_pools
               decomp_c14pools(c,k) = decomp_cpools(c,k) * c14ratio
               decomp_c14pools_1m(c,k) = decomp_cpools_1m(c,k) * c14ratio
            end do
         endif

         ! column nitrogen state variables
         sminn(c) = 0._r8
         do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
               decomp_npools_vr(c,j,k) = decomp_cpools_vr(c,j,k) / initial_cn_ratio(k)
            end do
            sminn_vr(c,j) = 0._r8
            col_ntrunc_vr(c,j) = 0._r8
         end do
         if ( nlevdecomp .gt. 1 ) then
            do j = nlevdecomp+1, nlevdecomp_full
               do k = 1, ndecomp_pools
                  decomp_npools_vr(c,j,k) = 0._r8
               end do
               sminn_vr(c,j) = 0._r8
               col_ntrunc_vr(c,j) = 0._r8
            end do
         end if
         do k = 1, ndecomp_pools
            decomp_npools(c,k) = decomp_cpools(c,k) / initial_cn_ratio(k)
            decomp_npools_1m(c,k) = decomp_cpools_1m(c,k) / initial_cn_ratio(k)
         end do

         if (use_nitrif_denitrif) then
            do j = 1, nlevdecomp_full
               smin_nh4_vr(c,j) = 0._r8
               smin_no3_vr(c,j) = 0._r8
            end do
            smin_nh4(c) = 0._r8
            smin_no3(c) = 0._r8
         end if
         totlitn(c)    = 0._r8
         totsomn(c)    = 0._r8
         totlitn_1m(c) = 0._r8
         totsomn_1m(c) = 0._r8
         totecosysn(c) = 0._r8
         totcoln(c)    = 0._r8
         cwdn(c)       = 0._r8

	 ! dynamic landcover state variables
     seedc(c)  = 0._r8
	 prod10c(c)    = 0._r8
	 prod100c(c)   = 0._r8
	 totprodc(c)   = 0._r8

         if ( use_c13 ) then
            seedc13(c)    = 0._r8
            prod10c13(c)  = 0._r8
            prod100c13(c) = 0._r8
            totprodc13(c) = 0._r8
         endif

         if ( use_c14 ) then
            seedc14(c)    = 0._r8
            prod10c14(c)  = 0._r8
            prod100c14(c) = 0._r8
            totprodc14(c) = 0._r8
         endif

	 seedn(c)      = 0._r8
	 prod10n(c)    = 0._r8
	 prod100n(c)   = 0._r8
	 totprodn(c)   = 0._r8
	 
	 ! also initialize dynamic landcover fluxes so that they have
	 ! real values on first timestep, prior to calling pftdyn_cnbal
	 ccf%dwt_seedc_to_leaf(c) = 0._r8
	 ccf%dwt_seedc_to_deadstem(c) = 0._r8
	 ccf%dwt_conv_cflux(c) = 0._r8
        ccf%lf_conv_cflux(c) = 0._r8
	 ccf%dwt_prod10c_gain(c) = 0._r8
	 ccf%prod10c_loss(c) = 0._r8
	 ccf%dwt_prod100c_gain(c) = 0._r8
	 ccf%prod100c_loss(c) = 0._r8
         do j = 1, nlevdecomp_full
            ccf%dwt_frootc_to_litr_met_c(c,j) = 0._r8
            ccf%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
            ccf%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
            ccf%dwt_livecrootc_to_cwdc(c,j) = 0._r8
            ccf%dwt_deadcrootc_to_cwdc(c,j) = 0._r8
          end do
	 ccf%dwt_closs(c) = 0._r8

         if ( use_c13 ) then
            cc13f%dwt_seedc_to_leaf(c) = 0._r8
            cc13f%dwt_seedc_to_deadstem(c) = 0._r8
            cc13f%dwt_conv_cflux(c) = 0._r8
            cc13f%dwt_prod10c_gain(c) = 0._r8
            cc13f%prod10c_loss(c) = 0._r8
            cc13f%dwt_prod100c_gain(c) = 0._r8
            cc13f%prod100c_loss(c) = 0._r8
            do j = 1, nlevdecomp_full
               cc13f%dwt_frootc_to_litr_met_c(c,j) = 0._r8
               cc13f%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
               cc13f%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
               cc13f%dwt_livecrootc_to_cwdc(c,j) = 0._r8
               cc13f%dwt_deadcrootc_to_cwdc(c,j) = 0._r8
            end do
            cc13f%dwt_closs(c) = 0._r8
         endif         
         
         if ( use_c14 ) then
            cc14f%dwt_seedc_to_leaf(c) = 0._r8
            cc14f%dwt_seedc_to_deadstem(c) = 0._r8
            cc14f%dwt_conv_cflux(c) = 0._r8
            cc14f%dwt_prod10c_gain(c) = 0._r8
            cc14f%prod10c_loss(c) = 0._r8
            cc14f%dwt_prod100c_gain(c) = 0._r8
            cc14f%prod100c_loss(c) = 0._r8
            do j = 1, nlevdecomp_full
               cc14f%dwt_frootc_to_litr_met_c(c,j) = 0._r8
               cc14f%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
               cc14f%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
               cc14f%dwt_livecrootc_to_cwdc(c,j) = 0._r8
               cc14f%dwt_deadcrootc_to_cwdc(c,j) = 0._r8
            end do
            cc14f%dwt_closs(c) = 0._r8
         endif

	 cnf%dwt_seedn_to_leaf(c) = 0._r8
	 cnf%dwt_seedn_to_deadstem(c) = 0._r8
	 cnf%dwt_conv_nflux(c) = 0._r8
	 cnf%dwt_prod10n_gain(c) = 0._r8
	 cnf%prod10n_loss(c) = 0._r8
	 cnf%dwt_prod100n_gain(c) = 0._r8
	 cnf%prod100n_loss(c) = 0._r8
         do j = 1, nlevdecomp_full
            cnf%dwt_frootn_to_litr_met_n(c,j) = 0._r8
            cnf%dwt_frootn_to_litr_cel_n(c,j) = 0._r8
            cnf%dwt_frootn_to_litr_lig_n(c,j) = 0._r8
            cnf%dwt_livecrootn_to_cwdn(c,j) = 0._r8
            cnf%dwt_deadcrootn_to_cwdn(c,j) = 0._r8
         end do
	 cnf%dwt_nloss(c) = 0._r8
      end if
   end do

   ! initialize pft-level variables
   do p = bounds%begp,bounds%endp
      l = pft%landunit(p)
      if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

         ! carbon state variables
         if (pft%itype(p) == noveg) then
            leafc(p) = 0._r8
            leafc_storage(p) = 0._r8
         else
            if (evergreen(pft%itype(p)) == 1._r8) then
               leafc(p) = 1._r8
               leafc_storage(p) = 0._r8
            else if (pft%itype(p) >= npcropmin) then ! prognostic crop types
               leafc(p) = 0._r8
               leafc_storage(p) = 0._r8
            else
               leafc(p) = 0._r8
               leafc_storage(p) = 1._r8
            end if
         end if
         leafc_xfer(p) = 0._r8
         if ( crop_prog )then
            grainc(p) = 0._r8
            grainc_storage(p) = 0._r8
            grainc_xfer(p) = 0._r8
            fert(p) = 0._r8
            soyfixn(p) = 0._r8
         end if
         frootc(p) = 0._r8
         frootc_storage(p) = 0._r8
         frootc_xfer(p) = 0._r8
         livestemc(p) = 0._r8
         livestemc_storage(p) = 0._r8
         livestemc_xfer(p) = 0._r8

         ! tree types need to be initialized with some stem mass so that
         ! roughness length is not zero in canopy flux calculation

         if (woody(pft%itype(p)) == 1._r8) then
            deadstemc(p) = 0.1_r8
         else
            deadstemc(p) = 0._r8
         end if

         deadstemc_storage(p) = 0._r8
         deadstemc_xfer(p) = 0._r8
         livecrootc(p) = 0._r8
         livecrootc_storage(p) = 0._r8
         livecrootc_xfer(p) = 0._r8
         deadcrootc(p) = 0._r8
         deadcrootc_storage(p) = 0._r8
         deadcrootc_xfer(p) = 0._r8
         gresp_storage(p) = 0._r8
         gresp_xfer(p) = 0._r8
         cpool(p) = 0._r8
         xsmrpool(p) = 0._r8
         pft_ctrunc(p) = 0._r8
         dispvegc(p) = 0._r8
         storvegc(p) = 0._r8
         totpftc(p)  = 0._r8
         ! calculate totvegc explicitly so that it is available for the isotope 
         ! code on the first time step.
         totvegc(p)  = leafc(p) + leafc_storage(p) + leafc_xfer(p) + frootc(p) +  &
            frootc_storage(p) + frootc_xfer(p) + livestemc(p) + livestemc_storage(p) +  &
            livestemc_xfer(p) + deadstemc(p) + deadstemc_storage(p) + deadstemc_xfer(p) +  &
            livecrootc(p) + livecrootc_storage(p) + livecrootc_xfer(p) + deadcrootc(p) +  &
            deadcrootc_storage(p) + deadcrootc_xfer(p) + gresp_storage(p) +  &
            gresp_xfer(p) + cpool(p)

         if ( crop_prog )then
            totvegc(p) = totvegc(p) + grainc(p) + grainc_storage(p) + grainc_xfer(p)
         end if

         woodc(p)    = 0._r8


         if ( use_c13 ) then
            leafc13(p)               = leafc(p)               * c13ratio
            leafc13_storage(p)       = leafc_storage(p)       * c13ratio
            leafc13_xfer(p)          = leafc_xfer(p)          * c13ratio
            frootc13(p)              = frootc(p)              * c13ratio
            frootc13_storage(p)      = frootc_storage(p)      * c13ratio
            frootc13_xfer(p)         = frootc_xfer(p)         * c13ratio
            livestemc13(p)           = livestemc(p)           * c13ratio
            livestemc13_storage(p)   = livestemc_storage(p)   * c13ratio
            livestemc13_xfer(p)      = livestemc_xfer(p)      * c13ratio
            deadstemc13(p)           = deadstemc(p)           * c13ratio
            deadstemc13_storage(p)   = deadstemc_storage(p)   * c13ratio
            deadstemc13_xfer(p)      = deadstemc_xfer(p)      * c13ratio
            livecrootc13(p)          = livecrootc(p)          * c13ratio
            livecrootc13_storage(p)  = livecrootc_storage(p)  * c13ratio
            livecrootc13_xfer(p)     = livecrootc_xfer(p)     * c13ratio
            deadcrootc13(p)          = deadcrootc(p)          * c13ratio
            deadcrootc13_storage(p)  = deadcrootc_storage(p)  * c13ratio
            deadcrootc13_xfer(p)     = deadcrootc_xfer(p)     * c13ratio
            c13_gresp_storage(p)     = gresp_storage(p)       * c13ratio
            c13_gresp_xfer(p)        = gresp_xfer(p)          * c13ratio
            c13pool(p)               = cpool(p)               * c13ratio
            c13xsmrpool(p)           = xsmrpool(p)            * c13ratio
            c13_pft_ctrunc(p)        = pft_ctrunc(p)          * c13ratio
            
            ! calculate totvegc explicitly so that it is available for the isotope 
            ! code on the first time step.
            totvegc13(p)  = leafc13(p) + leafc13_storage(p) + leafc13_xfer(p) + frootc13(p) +  &
                 frootc13_storage(p) + frootc13_xfer(p) + livestemc13(p) + livestemc13_storage(p) +  &
                 livestemc13_xfer(p) + deadstemc13(p) + deadstemc13_storage(p) + deadstemc13_xfer(p) +  &
                 livecrootc13(p) + livecrootc13_storage(p) + livecrootc13_xfer(p) + deadcrootc13(p) +  &
                 deadcrootc13_storage(p) + deadcrootc13_xfer(p) + c13_gresp_storage(p) +  &
                 c13_gresp_xfer(p) + c13pool(p)
         endif
            
         if ( use_c14 ) then
            leafc14(p)               = leafc(p)               * c14ratio
            leafc14_storage(p)       = leafc_storage(p)       * c14ratio
            leafc14_xfer(p)          = leafc_xfer(p)          * c14ratio
            frootc14(p)              = frootc(p)              * c14ratio
            frootc14_storage(p)      = frootc_storage(p)      * c14ratio
            frootc14_xfer(p)         = frootc_xfer(p)         * c14ratio
            livestemc14(p)           = livestemc(p)           * c14ratio
            livestemc14_storage(p)   = livestemc_storage(p)   * c14ratio
            livestemc14_xfer(p)      = livestemc_xfer(p)      * c14ratio
            deadstemc14(p)           = deadstemc(p)           * c14ratio
            deadstemc14_storage(p)   = deadstemc_storage(p)   * c14ratio
            deadstemc14_xfer(p)      = deadstemc_xfer(p)      * c14ratio
            livecrootc14(p)          = livecrootc(p)          * c14ratio
            livecrootc14_storage(p)  = livecrootc_storage(p)  * c14ratio
            livecrootc14_xfer(p)     = livecrootc_xfer(p)     * c14ratio
            deadcrootc14(p)          = deadcrootc(p)          * c14ratio
            deadcrootc14_storage(p)  = deadcrootc_storage(p)  * c14ratio
            deadcrootc14_xfer(p)     = deadcrootc_xfer(p)     * c14ratio
            c14_gresp_storage(p)     = gresp_storage(p)       * c14ratio
            c14_gresp_xfer(p)        = gresp_xfer(p)          * c14ratio
            c14pool(p)               = cpool(p)               * c14ratio
            c14xsmrpool(p)           = xsmrpool(p)            * c14ratio
            c14_pft_ctrunc(p)        = pft_ctrunc(p)          * c14ratio
            
            ! calculate totvegc explicitly so that it is available for the isotope 
            ! code on the first time step.
            totvegc14(p)  = leafc14(p) + leafc14_storage(p) + leafc14_xfer(p) + frootc14(p) +  &
                 frootc14_storage(p) + frootc14_xfer(p) + livestemc14(p) + livestemc14_storage(p) +  &
                 livestemc14_xfer(p) + deadstemc14(p) + deadstemc14_storage(p) + deadstemc14_xfer(p) +  &
                 livecrootc14(p) + livecrootc14_storage(p) + livecrootc14_xfer(p) + deadcrootc14(p) +  &
                 deadcrootc14_storage(p) + deadcrootc14_xfer(p) + c14_gresp_storage(p) +  &
                 c14_gresp_xfer(p) + c14pool(p)
            
            rc14_atm(p) = c14ratio
         endif
                                
         ! nitrogen state variables
         if (pft%itype(p) == noveg) then
            leafn(p) = 0._r8
            leafn_storage(p) = 0._r8
         else
            leafn(p) = leafc(p) / leafcn(pft%itype(p))
            leafn_storage(p) = leafc_storage(p) / leafcn(pft%itype(p))
         end if

         leafn_xfer(p) = 0._r8
         if ( crop_prog )then
            grainn(p) = 0._r8
            grainn_storage(p) = 0._r8
            grainn_xfer(p) = 0._r8
         end if
         frootn(p) = 0._r8
         frootn_storage(p) = 0._r8
         frootn_xfer(p) = 0._r8
         livestemn(p) = 0._r8
         livestemn_storage(p) = 0._r8
         livestemn_xfer(p) = 0._r8

         ! tree types need to be initialized with some stem mass so that
         ! roughness length is not zero in canopy flux calculation

         if (woody(pft%itype(p)) == 1._r8) then
            deadstemn(p) = deadstemc(p) / deadwdcn(pft%itype(p))
         else
            deadstemn(p) = 0._r8
         end if

         deadstemn_storage(p) = 0._r8
         deadstemn_xfer(p) = 0._r8
         livecrootn(p) = 0._r8
         livecrootn_storage(p) = 0._r8
         livecrootn_xfer(p) = 0._r8
         deadcrootn(p) = 0._r8
         deadcrootn_storage(p) = 0._r8
         deadcrootn_xfer(p) = 0._r8
         retransn(p) = 0._r8
         npool(p) = 0._r8
         pft_ntrunc(p) = 0._r8
         dispvegn(p) = 0._r8
         storvegn(p) = 0._r8
         totvegn(p)  = 0._r8
         totpftn(p)  = 0._r8

         ! initialization for psnsun and psnsha required for
         ! proper arbitrary initialization of allocation routine
         ! in initial ecosysdyn call

         psnsun(p) = 0._r8
         psnsha(p) = 0._r8

         if ( use_c13 ) then
            c13_psnsun(p) = 0._r8
            c13_psnsha(p) = 0._r8
         endif
         
         if ( use_c14 ) then
            c14_psnsun(p) = 0._r8
            c14_psnsha(p) = 0._r8
         endif

         laisun(p) = 0._r8
         laisha(p) = 0._r8

         ! ecophysiological variables
         ! phenology variables
         dormant_flag(p) = 1._r8
         days_active(p) = 0._r8
         onset_flag(p) = 0._r8
         onset_counter(p) = 0._r8
         onset_gddflag(p) = 0._r8
         onset_fdd(p) = 0._r8
         onset_gdd(p) = 0._r8
         onset_swi(p) = 0.0_r8
         offset_flag(p) = 0._r8
         offset_counter(p) = 0._r8
         offset_fdd(p) = 0._r8
         offset_swi(p) = 0._r8
         lgsf(p) = 0._r8
         bglfr(p) = 0._r8
         bgtr(p) = 0._r8
         annavg_t2m(p) = 280._r8
         tempavg_t2m(p) = 0._r8
         fert_counter(p) = 0._r8
         grain_flag(p) = 0._r8

         ! non-phenology variables
         gpp(p) = 0._r8
         availc(p) = 0._r8
         xsmrpool_recover(p) = 0._r8
         alloc_pnow(p) = 1._r8
         c_allometry(p) = 0._r8
         n_allometry(p) = 0._r8
         plant_ndemand(p) = 0._r8
         tempsum_potential_gpp(p) = 0._r8
         annsum_potential_gpp(p) = 0._r8
         tempmax_retransn(p) = 0._r8
         annmax_retransn(p) = 0._r8
         avail_retransn(p) = 0._r8
         plant_nalloc(p) = 0._r8
         plant_calloc(p) = 0._r8
         excess_cflux(p) = 0._r8
         downreg(p) = 0._r8
         prev_leafc_to_litter(p) = 0._r8
         prev_frootc_to_litter(p) = 0._r8
         tempsum_npp(p) = 0._r8
         annsum_npp(p) = 0._r8
         if (use_cndv) then
            tempsum_litfall(p) = 0._r8
            annsum_litfall(p) = 0._r8
         end if
         if ( use_c13 ) then
            xsmrpool_c13ratio(p) = c13ratio
            rc13_canair(p) = 0._r8
            rc13_psnsun(p) = 0._r8
            rc13_psnsha(p) = 0._r8
            alphapsnsun(p) = 0._r8
            alphapsnsha(p) = 0._r8
         endif
         
      end if   ! end of if-istsoil block
   end do   ! end of loop over pfts  

    end associate 
 end subroutine CNiniTimeVar
