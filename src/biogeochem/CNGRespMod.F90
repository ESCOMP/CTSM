module CNGRespMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for growth respiration fluxes,
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use pftconMod              , only : npcropmin, pftcon
  use CNVegcarbonfluxType    , only : cnveg_carbonflux_type
  use PatchType              , only : patch    
  use CanopyStateType        , only : canopystate_type              
  use CNVegCarbonStateType   , only : cnveg_carbonstate_type       
  use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type     
  use CropReprPoolsMod           , only : nrepr
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNGResp
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  ! subroutine CNGResp(num_soilp, filter_soilp, cnveg_carbonflux_inst)   
  subroutine CNGResp(num_soilp, filter_soilp, cnveg_carbonflux_inst, canopystate_inst, cnveg_carbonstate_inst, &
       cnveg_nitrogenstate_inst)  
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(canopystate_type)         , intent(in)    :: canopystate_inst            
    type(cnveg_carbonstate_type)   , intent(in)    :: cnveg_carbonstate_inst      
    type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst    
    !
    ! !LOCAL VARIABLES:
    integer :: p,k              ! indices
    integer :: fp               ! lake filter patch index
    real(r8):: respfact_leaf   				
    real(r8):: respfact_froot  				    
    real(r8):: respfact_livecroot  			
    real(r8):: respfact_livestem  			
    real(r8):: respfact_leaf_storage  		
    real(r8):: respfact_froot_storage 		    
    real(r8):: respfact_livecroot_storage  	
    real(r8):: respfact_livestem_storage  	
    !-----------------------------------------------------------------------

    associate(                                                                                           & 
         ivt                           =>    patch%itype                                                 , & ! Input:  [integer (:)]  patch vegetation type                                

         woody                         =>    pftcon%woody                                              , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         grperc                        =>    pftcon%grperc                                             , & ! Input:  growth respiration parameter
         grpnow                        =>    pftcon%grpnow                                             , & ! Input:  growth respiration parameter
         leafcn                        =>    pftcon%leafcn                                             , & ! Input:  leaf C:N (gC/gN)                         
         livewdcn                      =>    pftcon%livewdcn                                           , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN)  
         
         laisun                        =>    canopystate_inst%laisun_patch                             , & ! Input:  [real(r8) (:)]  sunlit projected leaf area index      
         laisha                        =>    canopystate_inst%laisha_patch                             , & ! Input:  [real(r8) (:)]  shaded projected leaf area index  
         
         leafc                 		   =>    cnveg_carbonstate_inst%leafc_patch                        , & ! Input:  [real(r8) (:)]                                       
         frootc                        =>    cnveg_carbonstate_inst%frootc_patch                       , & ! Input:  [real(r8) (:)]                             
         livestemc                     =>    cnveg_carbonstate_inst%livestemc_patch                    , & ! Input:  [real(r8) (:)]                             
         livecrootc                    =>    cnveg_carbonstate_inst%livecrootc_patch                   , & ! Input:  [real(r8) (:)]     
         leafc_storage                 =>    cnveg_carbonstate_inst%leafc_storage_patch                , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                           
         frootc_storage                =>    cnveg_carbonstate_inst%frootc_storage_patch               , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                        
         livestemc_storage             =>    cnveg_carbonstate_inst%livestemc_storage_patch            , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                                          
         livecrootc_storage            =>    cnveg_carbonstate_inst%livecrootc_storage_patch           , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage    
         
         leafn             			   =>    cnveg_nitrogenstate_inst%leafn_patch                      , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N  
         frootn                        =>    cnveg_nitrogenstate_inst%frootn_patch                     , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N 
         livestemn                     =>    cnveg_nitrogenstate_inst%livestemn_patch                  , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                                
         livecrootn                    =>    cnveg_nitrogenstate_inst%livecrootn_patch                 , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                  
         leafn_storage                 =>    cnveg_nitrogenstate_inst%leafn_storage_patch              , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                              
         frootn_storage                =>    cnveg_nitrogenstate_inst%frootn_storage_patch             , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                         
         livestemn_storage             =>    cnveg_nitrogenstate_inst%livestemn_storage_patch          , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                                              
         livecrootn_storage            =>    cnveg_nitrogenstate_inst%livecrootn_storage_patch         , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                  
           
         cpool_to_leafc                =>    cnveg_carbonflux_inst%cpool_to_leafc_patch                , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_leafc_storage        =>    cnveg_carbonflux_inst%cpool_to_leafc_storage_patch        , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_frootc               =>    cnveg_carbonflux_inst%cpool_to_frootc_patch               , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_frootc_storage       =>    cnveg_carbonflux_inst%cpool_to_frootc_storage_patch       , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_livestemc            =>    cnveg_carbonflux_inst%cpool_to_livestemc_patch            , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_livestemc_storage    =>    cnveg_carbonflux_inst%cpool_to_livestemc_storage_patch    , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_deadstemc            =>    cnveg_carbonflux_inst%cpool_to_deadstemc_patch            , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_deadstemc_storage    =>    cnveg_carbonflux_inst%cpool_to_deadstemc_storage_patch    , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_livecrootc           =>    cnveg_carbonflux_inst%cpool_to_livecrootc_patch           , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_livecrootc_storage   =>    cnveg_carbonflux_inst%cpool_to_livecrootc_storage_patch   , & ! Input:  [real(r8) (:)]                                                    
         cpool_to_deadcrootc           =>    cnveg_carbonflux_inst%cpool_to_deadcrootc_patch           , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C (gC/m2/s)        
         cpool_to_deadcrootc_storage   =>    cnveg_carbonflux_inst%cpool_to_deadcrootc_storage_patch   , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C storage (gC/m2/s)
         cpool_to_reproductivec               =>    cnveg_carbonflux_inst%cpool_to_reproductivec_patch               , & ! Input:  [real(r8) (:,:)]  allocation to grain C (gC/m2/s)
         cpool_to_reproductivec_storage       =>    cnveg_carbonflux_inst%cpool_to_reproductivec_storage_patch       , & ! Input:  [real(r8) (:,:)]  allocation to grain C storage (gC/m2/s)
         reproductivec_xfer_to_reproductivec         =>    cnveg_carbonflux_inst%reproductivec_xfer_to_reproductivec_patch         , & ! Input:  [real(r8) (:,:)]  grain C growth from storage (gC/m2/s)
         leafc_xfer_to_leafc           =>    cnveg_carbonflux_inst%leafc_xfer_to_leafc_patch           , & ! Input:  [real(r8) (:)]  leaf C growth from storage (gC/m2/s)              
         frootc_xfer_to_frootc         =>    cnveg_carbonflux_inst%frootc_xfer_to_frootc_patch         , & ! Input:  [real(r8) (:)]  fine root C growth from storage (gC/m2/s)         
         livestemc_xfer_to_livestemc   =>    cnveg_carbonflux_inst%livestemc_xfer_to_livestemc_patch   , & ! Input:  [real(r8) (:)]  live stem C growth from storage (gC/m2/s)         
         deadstemc_xfer_to_deadstemc   =>    cnveg_carbonflux_inst%deadstemc_xfer_to_deadstemc_patch   , & ! Input:  [real(r8) (:)]  dead stem C growth from storage (gC/m2/s)         
         livecrootc_xfer_to_livecrootc =>    cnveg_carbonflux_inst%livecrootc_xfer_to_livecrootc_patch , & ! Input:  [real(r8) (:)]  live coarse root C growth from storage (gC/m2/s)  
         deadcrootc_xfer_to_deadcrootc =>    cnveg_carbonflux_inst%deadcrootc_xfer_to_deadcrootc_patch , & ! Input:  [real(r8) (:)]  dead coarse root C growth from storage (gC/m2/s)  
         cpool_reproductive_gr                =>    cnveg_carbonflux_inst%cpool_reproductive_gr_patch                , & ! Output: [real(r8) (:,:)]
         cpool_reproductive_storage_gr        =>    cnveg_carbonflux_inst%cpool_reproductive_storage_gr_patch        , & ! Output: [real(r8) (:,:)]
         transfer_reproductive_gr             =>    cnveg_carbonflux_inst%transfer_reproductive_gr_patch             , & ! Output: [real(r8) (:,:)]
         cpool_leaf_gr                 =>    cnveg_carbonflux_inst%cpool_leaf_gr_patch                 , & ! Output: [real(r8) (:)]                                                    
         cpool_leaf_storage_gr         =>    cnveg_carbonflux_inst%cpool_leaf_storage_gr_patch         , & ! Output: [real(r8) (:)]                                                    
         transfer_leaf_gr              =>    cnveg_carbonflux_inst%transfer_leaf_gr_patch              , & ! Output: [real(r8) (:)]                                                    
         cpool_froot_gr                =>    cnveg_carbonflux_inst%cpool_froot_gr_patch                , & ! Output: [real(r8) (:)]                                                    
         cpool_froot_storage_gr        =>    cnveg_carbonflux_inst%cpool_froot_storage_gr_patch        , & ! Output: [real(r8) (:)]                                                    
         transfer_froot_gr             =>    cnveg_carbonflux_inst%transfer_froot_gr_patch             , & ! Output: [real(r8) (:)]                                                    
         cpool_livestem_gr             =>    cnveg_carbonflux_inst%cpool_livestem_gr_patch             , & ! Output: [real(r8) (:)]                                                    
         cpool_livestem_storage_gr     =>    cnveg_carbonflux_inst%cpool_livestem_storage_gr_patch     , & ! Output: [real(r8) (:)]                                                    
         transfer_livestem_gr          =>    cnveg_carbonflux_inst%transfer_livestem_gr_patch          , & ! Output: [real(r8) (:)]                                                    
         cpool_deadstem_gr             =>    cnveg_carbonflux_inst%cpool_deadstem_gr_patch             , & ! Output: [real(r8) (:)]                                                    
         cpool_deadstem_storage_gr     =>    cnveg_carbonflux_inst%cpool_deadstem_storage_gr_patch     , & ! Output: [real(r8) (:)]                                                    
         transfer_deadstem_gr          =>    cnveg_carbonflux_inst%transfer_deadstem_gr_patch          , & ! Output: [real(r8) (:)]                                                    
         cpool_livecroot_gr            =>    cnveg_carbonflux_inst%cpool_livecroot_gr_patch            , & ! Output: [real(r8) (:)]                                                    
         cpool_livecroot_storage_gr    =>    cnveg_carbonflux_inst%cpool_livecroot_storage_gr_patch    , & ! Output: [real(r8) (:)]                                                    
         transfer_livecroot_gr         =>    cnveg_carbonflux_inst%transfer_livecroot_gr_patch         , & ! Output: [real(r8) (:)]                                                    
         cpool_deadcroot_gr            =>    cnveg_carbonflux_inst%cpool_deadcroot_gr_patch            , & ! Output: [real(r8) (:)]                                                    
         cpool_deadcroot_storage_gr    =>    cnveg_carbonflux_inst%cpool_deadcroot_storage_gr_patch    , & ! Output: [real(r8) (:)]                                                    
         transfer_deadcroot_gr         =>    cnveg_carbonflux_inst%transfer_deadcroot_gr_patch           & ! Output: [real(r8) (:)]                                                    
         )
      
      ! Loop through patches
      ! start patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
          
         respfact_leaf = 1.0_r8   				
         respfact_froot = 1.0_r8 				    
         respfact_livecroot = 1.0_r8  			
         respfact_livestem = 1.0_r8  			
         respfact_livecroot = 1.0_r8 			   
         respfact_livestem = 1.0_r8 			
         respfact_leaf_storage = 1.0_r8 		
         respfact_froot_storage = 1.0_r8 		    
         respfact_livecroot_storage = 1.0_r8  	
         respfact_livestem_storage = 1.0_r8 	
         respfact_livecroot_storage = 1.0_r8 	
         respfact_livestem_storage = 1.0_r8 	
         
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cpool_livestem_gr(p) = cpool_to_livestemc(p) * grperc(ivt(p)) * respfact_livestem     

            cpool_livestem_storage_gr(p) = cpool_to_livestemc_storage(p) * grperc(ivt(p)) * grpnow(ivt(p)) * &
                 respfact_livestem_storage   

            transfer_livestem_gr(p) = livestemc_xfer_to_livestemc(p) * grperc(ivt(p)) * (1._r8 - grpnow(ivt(p))) * &
                 respfact_livestem_storage   

            do k = 1, nrepr
               cpool_reproductive_gr(p,k) = &
                    cpool_to_reproductivec(p,k) * grperc(ivt(p))
               cpool_reproductive_storage_gr(p,k) = &
                    cpool_to_reproductivec_storage(p,k) * grperc(ivt(p)) * grpnow(ivt(p))
               transfer_reproductive_gr(p,k) = &
                    reproductivec_xfer_to_reproductivec(p,k) * grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
            end do

         end if

         ! leaf and fine root growth respiration
         cpool_leaf_gr(p) = cpool_to_leafc(p) * grperc(ivt(p)) * respfact_leaf   

         cpool_leaf_storage_gr(p) = cpool_to_leafc_storage(p) * grperc(ivt(p)) * grpnow(ivt(p)) * respfact_leaf_storage   

         transfer_leaf_gr(p) = leafc_xfer_to_leafc(p) * grperc(ivt(p)) * (1._r8 - grpnow(ivt(p))) * respfact_leaf_storage   

         cpool_froot_gr(p) = cpool_to_frootc(p) * grperc(ivt(p)) * respfact_froot * respfact_froot   

         cpool_froot_storage_gr(p) = cpool_to_frootc_storage(p) * grperc(ivt(p)) * grpnow(ivt(p)) * respfact_froot_storage   

         transfer_froot_gr(p) = frootc_xfer_to_frootc(p) * grperc(ivt(p)) * (1._r8 - grpnow(ivt(p))) * respfact_froot_storage   

         if (woody(ivt(p)) == 1._r8) then
            cpool_livestem_gr(p) = cpool_to_livestemc(p) * grperc(ivt(p)) * respfact_livestem   

            cpool_livestem_storage_gr(p) = cpool_to_livestemc_storage(p) * grperc(ivt(p)) * grpnow(ivt(p)) * &
respfact_livestem_storage   

            transfer_livestem_gr(p) = livestemc_xfer_to_livestemc(p) * grperc(ivt(p)) * (1._r8 - grpnow(ivt(p))) * &
respfact_livestem_storage   

            cpool_deadstem_gr(p) = cpool_to_deadstemc(p) * grperc(ivt(p)) 

            cpool_deadstem_storage_gr(p) = cpool_to_deadstemc_storage(p) * grperc(ivt(p)) * grpnow(ivt(p)) 

            transfer_deadstem_gr(p) = deadstemc_xfer_to_deadstemc(p) * grperc(ivt(p)) * (1._r8 - grpnow(ivt(p))) 

            cpool_livecroot_gr(p) = cpool_to_livecrootc(p) * grperc(ivt(p)) * respfact_livecroot   

            cpool_livecroot_storage_gr(p) = cpool_to_livecrootc_storage(p) * grperc(ivt(p)) * grpnow(ivt(p)) * &
respfact_livecroot_storage   

            transfer_livecroot_gr(p) = livecrootc_xfer_to_livecrootc(p) * grperc(ivt(p)) * (1._r8 - grpnow(ivt(p))) * &
respfact_livecroot_storage   

            cpool_deadcroot_gr(p) = cpool_to_deadcrootc(p) * grperc(ivt(p)) 

            cpool_deadcroot_storage_gr(p) = cpool_to_deadcrootc_storage(p) * grperc(ivt(p)) * grpnow(ivt(p)) 

            transfer_deadcroot_gr(p) = deadcrootc_xfer_to_deadcrootc(p) * grperc(ivt(p)) * (1._r8 - grpnow(ivt(p))) 
         end if

      end do

    end associate

  end subroutine CNGResp

end module CNGRespMod
