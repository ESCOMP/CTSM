module CNVegNitrogenFluxType

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp, i_litr_min, i_litr_max
  use clm_varpar                         , only : nvegnpool,nnphtrans,nngmtrans,nnfitrans,nnphouttrans,&
                                                  nngmouttrans,nnfiouttrans
  use clm_varpar                         , only : ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                                  ilivestem,ilivestem_st,ilivestem_xf,&
                                                  ideadstem,ideadstem_st,ideadstem_xf,&
                                                  ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                                  ideadcroot,ideadcroot_st,ideadcroot_xf,&
                                                  igrain,igrain_st,igrain_xf,iretransn,ioutn
  use clm_varpar                         , only : mxharvests
  use clm_varcon                         , only : spval, ispval, dzsoi_decomp
  use clm_varctl                         , only : use_nitrif_denitrif, use_crop
  use CNSharedParamsMod                  , only : use_fun, use_matrixcn
  use decompMod                          , only : bounds_type
  use abortutils                         , only : endrun
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use dynSubgridControlMod               , only : get_do_grossunrep
  use CropReprPoolsMod                   , only : nrepr, repr_grain_min, repr_grain_max, repr_structure_min, repr_structure_max
  use CropReprPoolsMod                   , only : get_repr_hist_fname, get_repr_rest_fname, get_repr_longname
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col                
  use PatchType                          , only : patch                
  use SparseMatrixMultiplyMod            , only : sparse_matrix_type, diag_matrix_type, vector_type
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: cnveg_nitrogenflux_type

     ! gap mortality fluxes
     real(r8), pointer :: m_leafn_to_litter_patch                   (:)     ! patch leaf N mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_to_litter_patch                  (:)     ! patch fine root N mortality (gN/m2/s)
     real(r8), pointer :: m_leafn_storage_to_litter_patch           (:)     ! patch leaf N storage mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_storage_to_litter_patch          (:)     ! patch fine root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_storage_to_litter_patch       (:)     ! patch live stem N storage mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_storage_to_litter_patch       (:)     ! patch dead stem N storage mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_storage_to_litter_patch      (:)     ! patch live coarse root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_storage_to_litter_patch      (:)     ! patch dead coarse root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_leafn_xfer_to_litter_patch              (:)     ! patch leaf N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_xfer_to_litter_patch             (:)     ! patch fine root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_xfer_to_litter_patch          (:)     ! patch live stem N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_xfer_to_litter_patch          (:)     ! patch dead stem N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_xfer_to_litter_patch         (:)     ! patch live coarse root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_xfer_to_litter_patch         (:)     ! patch dead coarse root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_to_litter_patch               (:)     ! patch live stem N mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_to_litter_patch               (:)     ! patch dead stem N mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_to_litter_patch              (:)     ! patch live coarse root N mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_to_litter_patch              (:)     ! patch dead coarse root N mortality (gN/m2/s)
     real(r8), pointer :: m_retransn_to_litter_patch                (:)     ! patch retranslocated N pool mortality (gN/m2/s)

     ! harvest fluxes
     real(r8), pointer :: hrv_leafn_to_litter_patch                 (:)     ! patch leaf N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_to_litter_patch                (:)     ! patch fine root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_leafn_storage_to_litter_patch         (:)     ! patch leaf N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_storage_to_litter_patch        (:)     ! patch fine root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_storage_to_litter_patch     (:)     ! patch live stem N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_storage_to_litter_patch     (:)     ! patch dead stem N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_storage_to_litter_patch    (:)     ! patch live coarse root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_storage_to_litter_patch    (:)     ! patch dead coarse root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_leafn_xfer_to_litter_patch            (:)     ! patch leaf N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_xfer_to_litter_patch           (:)     ! patch fine root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_xfer_to_litter_patch        (:)     ! patch live stem N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_xfer_to_litter_patch        (:)     ! patch dead stem N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_xfer_to_litter_patch       (:)     ! patch live coarse root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_xfer_to_litter_patch       (:)     ! patch dead coarse root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_to_litter_patch             (:)     ! patch live stem N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_to_litter_patch            (:)     ! patch live coarse root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_to_litter_patch            (:)     ! patch dead coarse root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_retransn_to_litter_patch              (:)     ! patch retranslocated N pool harvest mortality (gN/m2/s)
     real(r8), pointer :: crop_harvestn_to_cropprodn_patch          (:)     ! patch crop harvest N to crop product pool (gN/m2/s)
     real(r8), pointer :: crop_harvestn_to_cropprodn_col            (:)     ! col crop harvest N to crop product pool (gN/m2/s)
     real(r8), pointer :: m_n_to_litr_fire_col                      (:,:,:) ! col N from leaf, froot, xfer and storage N to litter N by fire (gN/m3/s)
     real(r8), pointer :: harvest_n_to_litr_n_col                   (:,:,:) ! col N fluxes associated with harvest to litter pools (gN/m3/s)
     real(r8), pointer :: harvest_n_to_cwdn_col                     (:,:)   ! col N fluxes associated with harvest to CWD pool (gN/m3/s)

     ! fire N fluxes 
     real(r8), pointer :: m_decomp_npools_to_fire_vr_col            (:,:,:) ! col vertically-resolved decomposing N fire loss (gN/m3/s)
     real(r8), pointer :: m_decomp_npools_to_fire_col               (:,:)   ! col vertically-integrated (diagnostic) decomposing N fire loss (gN/m2/s)
     real(r8), pointer :: m_leafn_to_fire_patch                     (:)     ! patch (gN/m2/s) fire N emissions from leafn 
     real(r8), pointer :: m_leafn_storage_to_fire_patch             (:)     ! patch (gN/m2/s) fire N emissions from leafn_storage            
     real(r8), pointer :: m_leafn_xfer_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from leafn_xfer     
     real(r8), pointer :: m_livestemn_to_fire_patch                 (:)     ! patch (gN/m2/s) fire N emissions from livestemn 
     real(r8), pointer :: m_livestemn_storage_to_fire_patch         (:)     ! patch (gN/m2/s) fire N emissions from livestemn_storage      
     real(r8), pointer :: m_livestemn_xfer_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from livestemn_xfer
     real(r8), pointer :: m_deadstemn_to_fire_patch                 (:)     ! patch (gN/m2/s) fire N emissions from deadstemn
     real(r8), pointer :: m_deadstemn_storage_to_fire_patch         (:)     ! patch (gN/m2/s) fire N emissions from deadstemn_storage         
     real(r8), pointer :: m_deadstemn_xfer_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from deadstemn_xfer
     real(r8), pointer :: m_frootn_to_fire_patch                    (:)     ! patch (gN/m2/s) fire N emissions from frootn
     real(r8), pointer :: m_frootn_storage_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from frootn_storage
     real(r8), pointer :: m_frootn_xfer_to_fire_patch               (:)     ! patch (gN/m2/s) fire N emissions from frootn_xfer
     real(r8), pointer :: m_livecrootn_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from m_livecrootn_to_fire
     real(r8), pointer :: m_livecrootn_storage_to_fire_patch        (:)     ! patch (gN/m2/s) fire N emissions from livecrootn_storage     
     real(r8), pointer :: m_livecrootn_xfer_to_fire_patch           (:)     ! patch (gN/m2/s) fire N emissions from livecrootn_xfer
     real(r8), pointer :: m_deadcrootn_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn
     real(r8), pointer :: m_deadcrootn_storage_to_fire_patch        (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn_storage  
     real(r8), pointer :: m_deadcrootn_xfer_to_fire_patch           (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn_xfer
     real(r8), pointer :: m_retransn_to_fire_patch                  (:)     ! patch (gN/m2/s) fire N emissions from retransn
     real(r8), pointer :: m_leafn_to_litter_fire_patch              (:)     ! patch (gN/m2/s) from leafn to litter N  due to fire               
     real(r8), pointer :: m_leafn_storage_to_litter_fire_patch      (:)     ! patch (gN/m2/s) from leafn_storage to litter N  due to fire                              
     real(r8), pointer :: m_leafn_xfer_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from leafn_xfer to litter N  due to fire                              
     real(r8), pointer :: m_livestemn_to_litter_fire_patch          (:)     ! patch (gN/m2/s) from livestemn to litter N  due to fire                              
     real(r8), pointer :: m_livestemn_storage_to_litter_fire_patch  (:)     ! patch (gN/m2/s) from livestemn_storage to litter N  due to fire                                     
     real(r8), pointer :: m_livestemn_xfer_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from livestemn_xfer to litter N  due to fire                                     
     real(r8), pointer :: m_livestemn_to_deadstemn_fire_patch       (:)     ! patch (gN/m2/s) from livestemn to deadstemn N  due to fire                                     
     real(r8), pointer :: m_deadstemn_to_litter_fire_patch          (:)     ! patch (gN/m2/s) from deadstemn to litter N  due to fire                                     
     real(r8), pointer :: m_deadstemn_storage_to_litter_fire_patch  (:)     ! patch (gN/m2/s) from deadstemn_storage to litter N  due to fire                                               
     real(r8), pointer :: m_deadstemn_xfer_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from deadstemn_xfer to litter N  due to fire                                               
     real(r8), pointer :: m_frootn_to_litter_fire_patch             (:)     ! patch (gN/m2/s) from frootn to litter N  due to fire                                               
     real(r8), pointer :: m_frootn_storage_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from frootn_storage to litter N  due to fire                                               
     real(r8), pointer :: m_frootn_xfer_to_litter_fire_patch        (:)     ! patch (gN/m2/s) from frootn_xfer to litter N  due to fire                                               
     real(r8), pointer :: m_livecrootn_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from livecrootn to litter N  due to fire                                               
     real(r8), pointer :: m_livecrootn_storage_to_litter_fire_patch (:)     ! patch (gN/m2/s) from livecrootn_storage to litter N  due to fire                                                     
     real(r8), pointer :: m_livecrootn_xfer_to_litter_fire_patch    (:)     ! patch (gN/m2/s) from livecrootn_xfer to litter N  due to fire                                                     
     real(r8), pointer :: m_livecrootn_to_deadcrootn_fire_patch     (:)     ! patch (gN/m2/s) from livecrootn_xfer to deadcrootn due to fire                                                     
     real(r8), pointer :: m_deadcrootn_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from deadcrootn to deadcrootn due to fire                                                       
     real(r8), pointer :: m_deadcrootn_storage_to_litter_fire_patch (:)     ! patch (gN/m2/s) from deadcrootn_storage to deadcrootn due to fire                                                        
     real(r8), pointer :: m_deadcrootn_xfer_to_litter_fire_patch    (:)     ! patch (gN/m2/s) from deadcrootn_xfer to deadcrootn due to fire                                                         
     real(r8), pointer :: m_retransn_to_litter_fire_patch           (:)     ! patch (gN/m2/s) from retransn to deadcrootn due to fire                                                         
     real(r8), pointer :: fire_nloss_patch                          (:)     ! patch total patch-level fire N loss (gN/m2/s) 
     real(r8), pointer :: fire_nloss_col                            (:)     ! col total column-level fire N loss (gN/m2/s)
     real(r8), pointer :: fire_nloss_p2c_col                        (:)     ! col patch2col column-level fire N loss (gN/m2/s) (p2c)
     real(r8), pointer :: fire_mortality_n_to_cwdn_col              (:,:)   ! col N fluxes associated with fire mortality to CWD pool (gN/m3/s)

     ! phenology fluxes from transfer pool
     real(r8), pointer :: reproductiven_xfer_to_reproductiven_patch(:,:)     ! patch reproductive (e.g., grain) N growth from storage for prognostic crop model (gN/m2/s)
     real(r8), pointer :: leafn_xfer_to_leafn_patch                 (:)     ! patch leaf N growth from storage (gN/m2/s)
     real(r8), pointer :: frootn_xfer_to_frootn_patch               (:)     ! patch fine root N growth from storage (gN/m2/s)
     real(r8), pointer :: livestemn_xfer_to_livestemn_patch         (:)     ! patch live stem N growth from storage (gN/m2/s)
     real(r8), pointer :: deadstemn_xfer_to_deadstemn_patch         (:)     ! patch dead stem N growth from storage (gN/m2/s)
     real(r8), pointer :: livecrootn_xfer_to_livecrootn_patch       (:)     ! patch live coarse root N growth from storage (gN/m2/s)
     real(r8), pointer :: deadcrootn_xfer_to_deadcrootn_patch       (:)     ! patch dead coarse root N growth from storage (gN/m2/s)

     ! litterfall fluxes
     real(r8), pointer :: livestemn_to_litter_patch                 (:)     ! patch livestem N to litter (gN/m2/s)
     real(r8), pointer :: repr_grainn_to_food_patch               (:,:)     ! patch grain N to food for prognostic crop (gN/m2/s) [patch, repr_grain_min:repr_grain_max]
     real(r8), pointer :: repr_grainn_to_food_perharv_patch       (:,:,:)   ! grain N to food for prognostic crop accumulated by harvest (gN/m2) [patch, harvest, repr_grain_min:repr_grain_max]. Not per-second because this variable represents an accumulation over each growing season, to be instantaneously at the end of each calendar year, to provide output that's easier to work with.
     real(r8), pointer :: repr_grainn_to_food_thisyr_patch        (:,:)     ! grain N to food for prognostic crop accumulated this calendar year (gN/m2) [patch, repr_grain_min:repr_grain_max]. Not per-second because this variable represents an accumulation over an entire calendar year, to be saved instantaneously at the end of each calendar year, to provide output that's easier to work with.
     real(r8), pointer :: repr_structuren_to_cropprod_patch       (:,:)     ! patch reproductive structure N to crop product pool for prognostic crop (gN/m2/s) [patch, repr_structure_min:repr_structure_max]
     real(r8), pointer :: repr_structuren_to_litter_patch         (:,:)     ! patch reproductive structure N to litter for prognostic crop (gN/m2/s) [patch, repr_structure_min:repr_structure_max]
     real(r8), pointer :: leafn_to_biofueln_patch                   (:)     ! patch leaf N to biofuel N (gN/m2/s)
     real(r8), pointer :: livestemn_to_biofueln_patch               (:)     ! patch livestem N to biofuel N (gN/m2/s)
     real(r8), pointer :: leafn_to_removedresiduen_patch            (:)     ! patch leaf N to removed residue N (gN/m2/s)
     real(r8), pointer :: livestemn_to_removedresiduen_patch        (:)     ! patch livestem N to removed residue N (gN/m2/s)
     real(r8), pointer :: repr_grainn_to_seed_patch               (:,:)     ! patch grain N to seed for prognostic crop (gN/m2/s) [patch, repr_grain_min:repr_grain_max]
     real(r8), pointer :: repr_grainn_to_seed_perharv_patch       (:,:,:)   ! grain N to seed for prognostic crop accumulated by harvest (gN/m2) [patch, harvest, repr_grain_min:repr_grain_max]. Not per-second because this variable represents an accumulation over each growing season, to be instantaneously at the end of each calendar year, to provide output that's easier to work with.
     real(r8), pointer :: repr_grainn_to_seed_thisyr_patch        (:,:)     ! grain N to seed for prognostic crop accumulated this calendar year (gN/m2) [patch, repr_grain_min:repr_grain_max]. Not per-second because this variable represents an accumulation over an entire calendar year, to be saved instantaneously at the end of each calendar year, to provide output that's easier to work with.
     real(r8), pointer :: leafn_to_litter_patch                     (:)     ! patch leaf N litterfall (gN/m2/s)
     real(r8), pointer :: leafn_to_retransn_patch                   (:)     ! patch leaf N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: frootn_to_retransn_patch                  (:)     ! patch fine root N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: frootn_to_litter_patch                    (:)     ! patch fine root N litterfall (gN/m2/s)

     ! allocation fluxes
     real(r8), pointer :: retransn_to_npool_patch                   (:)     ! patch deployment of retranslocated N (gN/m2/s)  
     real(r8), pointer :: free_retransn_to_npool_patch              (:)     ! patch deployment of free retranslocated N (gN/m2/s)           
     real(r8), pointer :: sminn_to_npool_patch                      (:)     ! patch deployment of soil mineral N uptake (gN/m2/s)
     real(r8), pointer :: npool_to_reproductiven_patch      (:,:)     ! patch allocation to reproductive (e.g., grain) N for prognostic crop (gN/m2/s)
     real(r8), pointer :: npool_to_reproductiven_storage_patch(:,:)   ! patch allocation to reproductive (e.g., grain) N storage for prognostic crop (gN/m2/s)
     real(r8), pointer :: npool_to_leafn_patch                      (:)     ! patch allocation to leaf N (gN/m2/s)
     real(r8), pointer :: npool_to_leafn_storage_patch              (:)     ! patch allocation to leaf N storage (gN/m2/s)
     real(r8), pointer :: npool_to_frootn_patch                     (:)     ! patch allocation to fine root N (gN/m2/s)
     real(r8), pointer :: npool_to_frootn_storage_patch             (:)     ! patch allocation to fine root N storage (gN/m2/s)
     real(r8), pointer :: npool_to_livestemn_patch                  (:)     ! patch allocation to live stem N (gN/m2/s)
     real(r8), pointer :: npool_to_livestemn_storage_patch          (:)     ! patch allocation to live stem N storage (gN/m2/s)
     real(r8), pointer :: npool_to_deadstemn_patch                  (:)     ! patch allocation to dead stem N (gN/m2/s)
     real(r8), pointer :: npool_to_deadstemn_storage_patch          (:)     ! patch allocation to dead stem N storage (gN/m2/s)
     real(r8), pointer :: npool_to_livecrootn_patch                 (:)     ! patch allocation to live coarse root N (gN/m2/s)
     real(r8), pointer :: npool_to_livecrootn_storage_patch         (:)     ! patch allocation to live coarse root N storage (gN/m2/s)
     real(r8), pointer :: npool_to_deadcrootn_patch                 (:)     ! patch allocation to dead coarse root N (gN/m2/s)
     real(r8), pointer :: npool_to_deadcrootn_storage_patch         (:)     ! patch allocation to dead coarse root N storage (gN/m2/s)

     ! annual turnover of storage to transfer pools           
     real(r8), pointer :: reproductiven_storage_to_xfer_patch(:,:)    ! patch reproductive (e.g., grain) N shift storage to transfer for prognostic crop (gN/m2/s)
     real(r8), pointer :: leafn_storage_to_xfer_patch               (:)     ! patch leaf N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: frootn_storage_to_xfer_patch              (:)     ! patch fine root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: livestemn_storage_to_xfer_patch           (:)     ! patch live stem N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: deadstemn_storage_to_xfer_patch           (:)     ! patch dead stem N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: livecrootn_storage_to_xfer_patch          (:)     ! patch live coarse root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: deadcrootn_storage_to_xfer_patch          (:)     ! patch dead coarse root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: fert_patch                                (:)     ! patch applied fertilizer (gN/m2/s)
     real(r8), pointer :: fert_counter_patch                        (:)     ! patch >0 fertilize; <=0 not
     real(r8), pointer :: soyfixn_patch                             (:)     ! patch soybean fixed N (gN/m2/s)

     ! turnover of livewood to deadwood, with retranslocation 
     real(r8), pointer :: livestemn_to_deadstemn_patch              (:)     ! patch live stem N turnover (gN/m2/s)
     real(r8), pointer :: livestemn_to_retransn_patch               (:)     ! patch live stem N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: livecrootn_to_deadcrootn_patch            (:)     ! patch live coarse root N turnover (gN/m2/s)
     real(r8), pointer :: livecrootn_to_retransn_patch              (:)     ! patch live coarse root N to retranslocated N pool (gN/m2/s)

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: ndeploy_patch                             (:)     ! patch total N deployed to growth and storage (gN/m2/s)
     real(r8), pointer :: wood_harvestn_patch                       (:)     ! patch total N losses to wood product pools (gN/m2/s)
     real(r8), pointer :: wood_harvestn_col                         (:)     ! col total N losses to wood product pools (gN/m2/s) (p2c)
     ! phenology: litterfall and crop fluxes
     real(r8), pointer :: phenology_n_to_litr_n_col                 (:,:,:) ! col N fluxes associated with phenology (litterfall and crop) to litter pools (gN/m3/s)

     ! gap mortality fluxes
     real(r8), pointer :: gap_mortality_n_to_litr_n_col             (:,:,:) ! col N fluxes associated with gap mortality to litter pools (gN/m3/s)
     real(r8), pointer :: gap_mortality_n_to_cwdn_col               (:,:)   ! col N fluxes associated with gap mortality to CWD pool (gN/m3/s)

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seedn_to_leaf_patch                   (:)     ! (gN/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedn_to_leaf_grc                     (:)     ! (gN/m2/s) dwt_seedn_to_leaf_patch summed to the gridcell-level
     real(r8), pointer :: dwt_seedn_to_deadstem_patch               (:)     ! (gN/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedn_to_deadstem_grc                 (:)     ! (gN/m2/s) dwt_seedn_to_deadstem_patch summed to the gridcell-level
     real(r8), pointer :: dwt_conv_nflux_patch                      (:)     ! (gN/m2/s) conversion N flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_conv_nflux_grc                        (:)     ! (gN/m2/s) dwt_conv_nflux_patch summed to the gridcell-level
     real(r8), pointer :: dwt_wood_productn_gain_patch              (:)     ! patch (gN/m2/s) addition to wood product pools from landcover change; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_crop_productn_gain_patch              (:)     ! patch (gN/m2/s) addition to crop product pool from landcover change; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_frootn_to_litr_n_col                  (:,:,:) ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootn_to_cwdn_col                (:,:)   ! col (gN/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootn_to_cwdn_col                (:,:)   ! col (gN/m3/s) dead coarse root to CWD due to landcover change

     ! gross unrepresented landcover fluxes
     real(r8), pointer :: gru_leafn_to_litter_patch                 (:)     ! patch leaf N gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_leafn_storage_to_atm_patch            (:)     ! patch leaf N storage gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_leafn_xfer_to_atm_patch               (:)     ! patch leaf N transfer gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_frootn_to_litter_patch                (:)     ! patch fine root N gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_frootn_storage_to_atm_patch           (:)     ! patch fine root N storage gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_frootn_xfer_to_atm_patch              (:)     ! patch fine root N transfer gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_livestemn_to_atm_patch                (:)     ! patch live stem N gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_livestemn_storage_to_atm_patch        (:)     ! patch live stem N storage gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_livestemn_xfer_to_atm_patch           (:)     ! patch live stem N transfer gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_deadstemn_to_atm_patch                (:)     ! patch dead stem N gross unrepresented landcover change mortality to the atmosphere (gC/m2/s)
     real(r8), pointer :: gru_deadstemn_storage_to_atm_patch        (:)     ! patch dead stem N storage gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_deadstemn_xfer_to_atm_patch           (:)     ! patch dead stem N transfer gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_livecrootn_to_litter_patch            (:)     ! patch live coarse root N gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_livecrootn_storage_to_atm_patch       (:)     ! patch live coarse root N storage gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_livecrootn_xfer_to_atm_patch          (:)     ! patch live coarse root N transfer gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_deadcrootn_to_litter_patch            (:)     ! patch dead coarse root N gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_deadcrootn_storage_to_atm_patch       (:)     ! patch dead coarse root N storage gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_deadcrootn_xfer_to_atm_patch          (:)     ! patch dead coarse root N transfer gross unrepresented landcover change mortality (gN/m2/s)
     real(r8), pointer :: gru_retransn_to_litter_patch              (:)     ! patch retranslocated N pool gross unrepresented landcover change mortality (gN/m2/s)

     real(r8), pointer :: gru_conv_nflux_patch                      (:)     ! (gN/m2/s) conversion N flux (immediate loss to atm)
     real(r8), pointer :: gru_conv_nflux_col                        (:)     ! (gN/m2/s) conversion N flux (immediate loss to atm)
     real(r8), pointer :: gru_conv_nflux_grc                        (:)     ! (gN/m2/s) gru_conv_nflux_patch summed to the gridcell-level
     real(r8), pointer :: gru_wood_productn_gain_patch              (:)     ! patch (gN/m2/s) addition to wood product pools from gross unrepresented landcover change
     real(r8), pointer :: gru_wood_productn_gain_col                (:)     ! column (gN/m2/s) addition to wood product pools from gross unrepresented landcover change
     real(r8), pointer :: gru_wood_productn_gain_grc                (:)     ! gridcell (gN/m2/s) addition to wood product pools from gross unrepresented landcover change
     real(r8), pointer :: gru_n_to_litr_n_col                      (:,:,:)  ! col (gN/m3/s) N to litter due to gross unrepresented landcover change
     real(r8), pointer :: gru_n_to_cwdn_col                         (:,:)   ! col (gN/m3/s) N to CWD due to gross unrepresented landcover change

     ! crop fluxes
     real(r8), pointer :: crop_seedn_to_leaf_patch                  (:)     ! patch (gN/m2/s) seed source to leaf, for crops
     
     ! Misc
     real(r8), pointer :: plant_ndemand_patch                       (:)     ! N flux required to support initial GPP (gN/m2/s)
     real(r8), pointer :: avail_retransn_patch                      (:)     ! N flux available from retranslocation pool (gN/m2/s)
     real(r8), pointer :: plant_nalloc_patch                        (:)     ! total allocated N flux (gN/m2/s)
     real(r8), pointer :: plant_ndemand_retrans_patch               (:)     ! The N demand pool generated for FUN2.0; mainly used for deciduous trees (gN/m2/s)
     real(r8), pointer :: plant_ndemand_season_patch                (:)     ! The N demand pool for seasonal deciduous (gN/m2/s)
     real(r8), pointer :: plant_ndemand_stress_patch                (:)     ! The N demand pool for stress deciduous   (gN/m2/s)
     real(r8), pointer :: Nactive_patch                             (:)     ! N acquired by mycorrhizal uptake  (gN/m2/s)
     real(r8), pointer :: Nnonmyc_patch                             (:)     ! N acquired by non-myc uptake      (gN/m2/s)
     real(r8), pointer :: Nam_patch                                 (:)     ! N acquired by AM plant            (gN/m2/s)
     real(r8), pointer :: Necm_patch                                (:)     ! N acquired by ECM plant           (gN/m2/s)
     real(r8), pointer :: Nactive_no3_patch                         (:)     ! N acquired by mycorrhizal uptake  (gN/m2/s)
     real(r8), pointer :: Nactive_nh4_patch                         (:)     ! N acquired by mycorrhizal uptake  (gN/m2/s)
     real(r8), pointer :: Nnonmyc_no3_patch                         (:)     ! N acquired by non-myc             (gN/m2/s)
     real(r8), pointer :: Nnonmyc_nh4_patch                         (:)     ! N acquired by non-myc             (gN/m2/s)
     real(r8), pointer :: Nam_no3_patch                             (:)     ! N acquired by AM plant            (gN/m2/s)
     real(r8), pointer :: Nam_nh4_patch                             (:)     ! N acquired by AM plant            (gN/m2/s)
     real(r8), pointer :: Necm_no3_patch                            (:)     ! N acquired by ECM plant           (gN/m2/s)
     real(r8), pointer :: Necm_nh4_patch                            (:)     ! N acquired by ECM plant           (gN/m2/s)
     real(r8), pointer :: Nfix_patch                                (:)     ! N acquired by Symbiotic BNF       (gN/m2/s)
     real(r8), pointer :: Npassive_patch                            (:)     ! N acquired by passive uptake      (gN/m2/s)
     real(r8), pointer :: Nretrans_patch                            (:)     ! N acquired by retranslocation     (gN/m2/s)
     real(r8), pointer :: Nretrans_org_patch                        (:)     ! N acquired by retranslocation     (gN/m2/s)
     real(r8), pointer :: Nretrans_season_patch                     (:)     ! N acquired by retranslocation     (gN/m2/s)
     real(r8), pointer :: Nretrans_stress_patch                     (:)     ! N acquired by retranslocation     (gN/m2/s)
     real(r8), pointer :: Nuptake_patch                             (:)     ! Total N uptake of FUN             (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_patch                  (:)     ! Total soil N uptake of FUN        (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_vr_patch               (:,:)   ! Total layer soil N uptake of FUN  (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_no3_vr_patch           (:,:)   ! Total layer no3 uptake of FUN     (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_nh4_vr_patch           (:,:)   ! Total layer nh4 uptake of FUN     (gN/m2/s)
     real(r8), pointer :: cost_nfix_patch                           (:)     ! Average cost of fixation          (gN/m2/s)
     real(r8), pointer :: cost_nactive_patch                        (:)     ! Average cost of active uptake     (gN/m2/s)
     real(r8), pointer :: cost_nretrans_patch                       (:)     ! Average cost of retranslocation   (gN/m2/s)
     real(r8), pointer :: nuptake_npp_fraction_patch                (:)     ! frac of npp spent on N acquisition   (gN/m2/s)
	 ! Matrix
     real(r8), pointer :: matrix_nalloc_patch                       (:,:)   ! B-matrix for nitrogen allocation
     real(r8), pointer :: matrix_Ninput_patch                       (:)     ! I-matrix for nitrogen input
	 
     real(r8), pointer :: matrix_nphtransfer_patch                  (:,:)   ! A-matrix_phenologh for nitrogen
     real(r8), pointer :: matrix_nphturnover_patch                  (:,:)   ! K-matrix_phenologh for nitrogen
     integer,  pointer :: matrix_nphtransfer_doner_patch            (:)     ! A-matrix_phenology non-zero indices (column indices) for nitrogen
     integer,  pointer :: matrix_nphtransfer_receiver_patch         (:)     ! A-matrix_phenology non-zero indices (row indices) for nitrogen
 
     real(r8), pointer :: matrix_ngmtransfer_patch                  (:,:)   ! A-matrix_gap mortality for nitrogen
     real(r8), pointer :: matrix_ngmturnover_patch                  (:,:)   ! K-matrix_gap mortality for nitrogen 
     integer,  pointer :: matrix_ngmtransfer_doner_patch            (:)     ! A-matrix_gap mortality non-zero indices (column indices) for nitrogen
     integer,  pointer :: matrix_ngmtransfer_receiver_patch         (:)     ! A-matrix_gap mortality non-zero indices (row indices) for nitrogen
  
     real(r8), pointer :: matrix_nfitransfer_patch                  (:,:)   ! A-matrix_fire for nitrogen
     real(r8), pointer :: matrix_nfiturnover_patch                  (:,:)   ! K-matrix_fire for nitrogen
     integer,  pointer :: matrix_nfitransfer_doner_patch            (:)     ! A-matrix_fire non-zero indices (column indices) for nitrogen
     integer,  pointer :: matrix_nfitransfer_receiver_patch         (:)     ! A-matrix_fire non-zero indices (row indices) for nitrogen

     integer ileafst_to_ileafxf_ph                    ! Index of phenology related N transfer from leaf storage pool to leaf transfer pool
     integer ileafxf_to_ileaf_ph                      ! Index of phenology related N transfer from leaf transfer pool to leaf pool  
     integer ifrootst_to_ifrootxf_ph                  ! Index of phenology related N transfer from fine root storage pool to fine root transfer pool
     integer ifrootxf_to_ifroot_ph                    ! Index of phenology related N transfer from fine root transfer pool to fine root pool  
     integer ilivestemst_to_ilivestemxf_ph            ! Index of phenology related N transfer from live stem storage pool to live stem transfer pool
     integer ilivestemxf_to_ilivestem_ph              ! Index of phenology related N transfer from live stem transfer pool to live stem pool  
     integer ideadstemst_to_ideadstemxf_ph            ! Index of phenology related N transfer from dead stem storage pool to dead stem transfer pool
     integer ideadstemxf_to_ideadstem_ph              ! Index of phenology related N transfer from dead stem transfer pool to dead stem pool  
     integer ilivecrootst_to_ilivecrootxf_ph          ! Index of phenology related N transfer from live coarse root storage pool to live coarse root transfer pool
     integer ilivecrootxf_to_ilivecroot_ph            ! Index of phenology related N transfer from live coarse root transfer pool to live coarse root pool  
     integer ideadcrootst_to_ideadcrootxf_ph          ! Index of phenology related N transfer from dead coarse root storage pool to dead coarse root transfer pool
     integer ideadcrootxf_to_ideadcroot_ph            ! Index of phenology related N transfer from dead coarse root transfer pool to dead coarse root pool  
     integer ilivestem_to_ideadstem_ph                ! Index of phenology related N transfer from live stem pool to dead stem pool  
     integer ilivecroot_to_ideadcroot_ph              ! Index of phenology related N transfer from live coarse root pool to dead coarse root pool  
     integer iretransn_to_ileaf_ph                    ! Index of phenology related N transfer from retranslocation pool to leaf pool
     integer iretransn_to_ileafst_ph                  ! Index of phenology related N transfer from retranslocation pool to leaf storage pool
     integer iretransn_to_ifroot_ph                   ! Index of phenology related N transfer from retranslocation pool to fine root pool
     integer iretransn_to_ifrootst_ph                 ! Index of phenology related N transfer from retranslocation pool to fine root storage pool
     integer iretransn_to_ilivestem_ph                ! Index of phenology related N transfer from retranslocation pool to live stem pool
     integer iretransn_to_ilivestemst_ph              ! Index of phenology related N transfer from retranslocation pool to live stem storage pool
     integer iretransn_to_ideadstem_ph                ! Index of phenology related N transfer from retranslocation pool to dead stem pool
     integer iretransn_to_ideadstemst_ph              ! Index of phenology related N transfer from retranslocation pool to dead stem storage pool
     integer iretransn_to_ilivecroot_ph               ! Index of phenology related N transfer from retranslocation pool to live coarse root pool
     integer iretransn_to_ilivecrootst_ph             ! Index of phenology related N transfer from retranslocation pool to live coarse root storage pool
     integer iretransn_to_ideadcroot_ph               ! Index of phenology related N transfer from retranslocation pool to dead coarse root pool
     integer iretransn_to_ideadcrootst_ph             ! Index of phenology related N transfer from retranslocation pool to dead coarse root storage pool
     integer iretransn_to_igrain_ph                   ! Index of phenology related N transfer from retranslocation pool to grain pool
     integer iretransn_to_igrainst_ph                 ! Index of phenology related N transfer from retranslocation pool to grain storage pool
     integer ileaf_to_iout_ph                         ! Index of phenology related N transfer from leaf pool to outside of vegetation pools  
     integer ifroot_to_iout_ph                        ! Index of phenology related N transfer from fine root pool to outside of vegetation pools  
     integer ilivestem_to_iout_ph                     ! Index of phenology related N transfer from live stem pool to outside of vegetation pools  
     integer ileaf_to_iretransn_ph                    ! Index of phenology related N transfer from leaf pool to retranslocation pools
     integer ifroot_to_iretransn_ph                   ! Index of phenology related N transfer from fine root pool to retranslocation pools
     integer ilivestem_to_iretransn_ph                ! Index of phenology related N transfer from live stem pool to retranslocation pools
     integer ilivecroot_to_iretransn_ph               ! Index of phenology related N transfer from live coarse root pool to retranslocation pools
     integer igrain_to_iout_ph                        ! Index of phenology related N transfer from grain pool to outside of vegetation pools  
     integer iretransn_to_iout_ph                     ! Index of phenology related N transfer from retranslocation pool to outside of vegetation pools
     integer ileaf_to_iout_gm                         ! Index of gap mortality related N transfer from leaf pool to outside of vegetation pools
     integer ileafst_to_iout_gm                       ! Index of gap mortality related N transfer from leaf storage pool to outside of vegetation pools
     integer ileafxf_to_iout_gm                       ! Index of gap mortality related N transfer from leaf transfer pool to outside of vegetation pools
     integer ifroot_to_iout_gm                        ! Index of gap mortality related N transfer from fine root pool to outside of vegetation pools
     integer ifrootst_to_iout_gm                      ! Index of gap mortality related N transfer from fine root storage pool to outside of vegetation pools
     integer ifrootxf_to_iout_gm                      ! Index of gap mortality related N transfer from fine root transfer pool to outside of vegetation pools
     integer ilivestem_to_iout_gm                     ! Index of gap mortality related N transfer from live stem pool to outside of vegetation pools
     integer ilivestemst_to_iout_gm                   ! Index of gap mortality related N transfer from live stem storage pool to outside of vegetation pools
     integer ilivestemxf_to_iout_gm                   ! Index of gap mortality related N transfer from live stem transfer pool to outside of vegetation pools
     integer ideadstem_to_iout_gm                     ! Index of gap mortality related N transfer from dead stem pool to outside of vegetation pools
     integer ideadstemst_to_iout_gm                   ! Index of gap mortality related N transfer from dead stem storage pool to outside of vegetation pools
     integer ideadstemxf_to_iout_gm                   ! Index of gap mortality related N transfer from dead stem transfer pool to outside of vegetation pools
     integer ilivecroot_to_iout_gm                    ! Index of gap mortality related N transfer from live coarse root pool to outside of vegetation pools
     integer ilivecrootst_to_iout_gm                  ! Index of gap mortality related N transfer from live coarse root storage pool to outside of vegetation pools
     integer ilivecrootxf_to_iout_gm                  ! Index of gap mortality related N transfer from live coarse root transfer pool to outside of vegetation pools
     integer ideadcroot_to_iout_gm                    ! Index of gap mortality related N transfer from dead coarse root pool to outside of vegetation pools
     integer ideadcrootst_to_iout_gm                  ! Index of gap mortality related N transfer from dead coarse root storage pool to outside of vegetation pools
     integer ideadcrootxf_to_iout_gm                  ! Index of gap mortality related N transfer from dead coarse root transfer pool to outside of vegetation pools
     integer iretransn_to_iout_gm                     ! Index of gap mortality related N transfer from retranslocation to outside of vegetation pools
     integer ileaf_to_iout_fi                         ! Index of fire related N transfer from leaf pool to outside of vegetation pools
     integer ileafst_to_iout_fi                       ! Index of fire related N transfer from leaf storage pool to outside of vegetation pools
     integer ileafxf_to_iout_fi                       ! Index of fire related N transfer from leaf transfer pool to outside of vegetation pools
     integer ifroot_to_iout_fi                        ! Index of fire related N transfer from fine root pool to outside of vegetation pools
     integer ifrootst_to_iout_fi                      ! Index of fire related N transfer from fine root storage pool to outside of vegetation pools
     integer ifrootxf_to_iout_fi                      ! Index of fire related N transfer from fine root transfer pool to outside of vegetation pools
     integer ilivestem_to_iout_fi                     ! Index of fire related N transfer from live stem pool to outside of vegetation pools
     integer ilivestemst_to_iout_fi                   ! Index of fire related N transfer from live stem storage pool to outside of vegetation pools
     integer ilivestemxf_to_iout_fi                   ! Index of fire related N transfer from live stem transfer pool to outside of vegetation pools
     integer ideadstem_to_iout_fi                     ! Index of fire related N transfer from dead stem pool to outside of vegetation pools
     integer ideadstemst_to_iout_fi                   ! Index of fire related N transfer from dead stem storage pool to outside of vegetation pools
     integer ideadstemxf_to_iout_fi                   ! Index of fire related N transfer from dead stem transfer pool to outside of vegetation pools
     integer ilivecroot_to_iout_fi                    ! Index of fire related N transfer from live coarse root pool to outside of vegetation pools
     integer ilivecrootst_to_iout_fi                  ! Index of fire related N transfer from live coarse root storage pool to outside of vegetation pools
     integer ilivecrootxf_to_iout_fi                  ! Index of fire related N transfer from live coarse root transfer pool to outside of vegetation pools
     integer ideadcroot_to_iout_fi                    ! Index of fire related N transfer from dead coarse root pool to outside of vegetation pools
     integer ideadcrootst_to_iout_fi                  ! Index of fire related N transfer from dead coarse root storage pool to outside of vegetation pools
     integer ideadcrootxf_to_iout_fi                  ! Index of fire related N transfer from dead coarse root transfer pool to outside of vegetation pools
     integer iretransn_to_iout_fi                     ! Index of fire related N transfer from retranslocation transfer pool to outside of vegetation pools
     integer ilivestem_to_ideadstem_fi                ! Index of fire related N transfer from live stem pool to dead stem pools
     integer ilivecroot_to_ideadcroot_fi              ! Index of fire related N transfer from live coarse root pool to dead coarse root pools

     integer,pointer :: list_phn_phgmn     (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKphn to AKphn+AKgmn
     integer,pointer :: list_gmn_phgmn     (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKgmn to AKphn+AKgmn
     integer,pointer :: list_phn_phgmfin   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKphn to AKphn+AKgmn+AKfin
     integer,pointer :: list_gmn_phgmfin   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKgmn to AKphn+AKgmn+AKfin
     integer,pointer :: list_fin_phgmfin   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKfin to AKphn+AKgmn+AKfin
     integer,pointer :: list_aphn          (:)        ! Indices of non-diagnoal entries in full sparse matrix Aph for N cycle
     integer,pointer :: list_agmn          (:)        ! Indices of non-diagnoal entries in full sparse matrix Agm for N cycle
     integer,pointer :: list_afin          (:)        ! Indices of non-diagnoal entries in full sparse matrix Afi for N cycle

     type(sparse_matrix_type)      :: AKphvegn        ! Aph*Kph for N cycle in sparse matrix format
     type(sparse_matrix_type)      :: AKgmvegn        ! Agm*Kgm for N cycle in sparse matrix format
     type(sparse_matrix_type)      :: AKfivegn        ! Afi*Kfi for N cycle in sparse matrix format
     type(sparse_matrix_type)      :: AKallvegn       ! Aph*Kph + Agm*Kgm + Afi*Kfi for N cycle in sparse matrix format
     integer                       :: NE_AKallvegn    ! Number of entries in AKallvegn
     integer,pointer,dimension(:)  :: RI_AKallvegn    ! Row indices in Akallvegn
     integer,pointer,dimension(:)  :: CI_AKallvegn    ! Column indices in AKallvegn
     integer,pointer,dimension(:)  :: RI_phn          ! Row indices of non-diagonal entires in Aph for N cycle
     integer,pointer,dimension(:)  :: CI_phn          ! Column indices of non-diagonal entries in Aph for N cycle
     integer,pointer,dimension(:)  :: RI_gmn          ! Row indices of non-diagonal entires in Agm for N cycle
     integer,pointer,dimension(:)  :: CI_gmn          ! Column indices of non-diagonal entries in Agm for N cycle
     integer,pointer,dimension(:)  :: RI_fin          ! Row indices of non-diagonal entires in Afi for N cycle
     integer,pointer,dimension(:)  :: CI_fin          ! Column indices of non-diagonal entries in Afi for N cycle
     type(diag_matrix_type)        :: Kvegn           ! Temporary variable of Kph, Kgm or Kfi for N cycle in diagonal matrix format
     type(vector_type)             :: Xvegn           ! Vegetation N of each compartment in a vector format

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: ZeroGRU
     procedure , public  :: Summary => Summary_nitrogenflux
     procedure , private :: InitAllocate 
     procedure , private :: InitTransfer 
     procedure , private :: InitHistory
     procedure , private :: InitCold

  end type cnveg_nitrogenflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, alloc_full_veg)

    class(cnveg_nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    logical,intent(in)            :: alloc_full_veg

    call this%InitAllocate (bounds,alloc_full_veg)
    if(alloc_full_veg)then
       if(use_matrixcn)then
          call this%InitTransfer ()
       end if
       call this%InitHistory (bounds)
       call this%InitCold (bounds)
    end if
    
  end subroutine Init

  subroutine InitTransfer (this)
    !
    ! !AGRUMENTS:
    class (cnveg_nitrogenflux_type) :: this
    
    this%ileaf_to_iretransn_ph           = 1
    this%matrix_nphtransfer_doner_patch(this%ileaf_to_iretransn_ph)              = ileaf
    this%matrix_nphtransfer_receiver_patch(this%ileaf_to_iretransn_ph)           = iretransn
    
    this%ileafst_to_ileafxf_ph           = 2
    this%matrix_nphtransfer_doner_patch(this%ileafst_to_ileafxf_ph)              = ileaf_st
    this%matrix_nphtransfer_receiver_patch(this%ileafst_to_ileafxf_ph)           = ileaf_xf

    this%ileafxf_to_ileaf_ph             = 3
    this%matrix_nphtransfer_doner_patch(this%ileafxf_to_ileaf_ph)                = ileaf_xf
    this%matrix_nphtransfer_receiver_patch(this%ileafxf_to_ileaf_ph)             = ileaf

    this%ifroot_to_iretransn_ph          = 4
    this%matrix_nphtransfer_doner_patch(this%ifroot_to_iretransn_ph)             = ifroot
    this%matrix_nphtransfer_receiver_patch(this%ifroot_to_iretransn_ph)          = iretransn
    
    this%ifrootst_to_ifrootxf_ph         = 5
    this%matrix_nphtransfer_doner_patch(this%ifrootst_to_ifrootxf_ph)            = ifroot_st
    this%matrix_nphtransfer_receiver_patch(this%ifrootst_to_ifrootxf_ph)         = ifroot_xf

    this%ifrootxf_to_ifroot_ph           = 6
    this%matrix_nphtransfer_doner_patch(this%ifrootxf_to_ifroot_ph)              = ifroot_xf
    this%matrix_nphtransfer_receiver_patch(this%ifrootxf_to_ifroot_ph)           = ifroot

    this%ilivestem_to_ideadstem_ph       = 7
    this%matrix_nphtransfer_doner_patch(this%ilivestem_to_ideadstem_ph)          = ilivestem
    this%matrix_nphtransfer_receiver_patch(this%ilivestem_to_ideadstem_ph)       = ideadstem
    
    this%ilivestem_to_iretransn_ph       = 8
    this%matrix_nphtransfer_doner_patch(this%ilivestem_to_iretransn_ph)          = ilivestem
    this%matrix_nphtransfer_receiver_patch(this%ilivestem_to_iretransn_ph)       = iretransn
    
    this%ilivestemst_to_ilivestemxf_ph   = 9
    this%matrix_nphtransfer_doner_patch(this%ilivestemst_to_ilivestemxf_ph)      = ilivestem_st
    this%matrix_nphtransfer_receiver_patch(this%ilivestemst_to_ilivestemxf_ph)   = ilivestem_xf

    this%ilivestemxf_to_ilivestem_ph     = 10
    this%matrix_nphtransfer_doner_patch(this%ilivestemxf_to_ilivestem_ph)        = ilivestem_xf
    this%matrix_nphtransfer_receiver_patch(this%ilivestemxf_to_ilivestem_ph)     = ilivestem

    this%ideadstemst_to_ideadstemxf_ph   = 11
    this%matrix_nphtransfer_doner_patch(this%ideadstemst_to_ideadstemxf_ph)      = ideadstem_st
    this%matrix_nphtransfer_receiver_patch(this%ideadstemst_to_ideadstemxf_ph)   = ideadstem_xf

    this%ideadstemxf_to_ideadstem_ph     = 12
    this%matrix_nphtransfer_doner_patch(this%ideadstemxf_to_ideadstem_ph)        = ideadstem_xf
    this%matrix_nphtransfer_receiver_patch(this%ideadstemxf_to_ideadstem_ph)     = ideadstem

    this%ilivecroot_to_ideadcroot_ph     = 13
    this%matrix_nphtransfer_doner_patch(this%ilivecroot_to_ideadcroot_ph)        = ilivecroot
    this%matrix_nphtransfer_receiver_patch(this%ilivecroot_to_ideadcroot_ph)     = ideadcroot
    
    this%ilivecroot_to_iretransn_ph      = 14
    this%matrix_nphtransfer_doner_patch(this%ilivecroot_to_iretransn_ph)         = ilivecroot
    this%matrix_nphtransfer_receiver_patch(this%ilivecroot_to_iretransn_ph)      = iretransn
    
    this%ilivecrootst_to_ilivecrootxf_ph = 15
    this%matrix_nphtransfer_doner_patch(this%ilivecrootst_to_ilivecrootxf_ph)    = ilivecroot_st
    this%matrix_nphtransfer_receiver_patch(this%ilivecrootst_to_ilivecrootxf_ph) = ilivecroot_xf

    this%ilivecrootxf_to_ilivecroot_ph   = 16
    this%matrix_nphtransfer_doner_patch(this%ilivecrootxf_to_ilivecroot_ph)      = ilivecroot_xf
    this%matrix_nphtransfer_receiver_patch(this%ilivecrootxf_to_ilivecroot_ph)   = ilivecroot

    this%ideadcrootst_to_ideadcrootxf_ph = 17
    this%matrix_nphtransfer_doner_patch(this%ideadcrootst_to_ideadcrootxf_ph)    = ideadcroot_st
    this%matrix_nphtransfer_receiver_patch(this%ideadcrootst_to_ideadcrootxf_ph) = ideadcroot_xf

    this%ideadcrootxf_to_ideadcroot_ph   = 18
    this%matrix_nphtransfer_doner_patch(this%ideadcrootxf_to_ideadcroot_ph)      = ideadcroot_xf
    this%matrix_nphtransfer_receiver_patch(this%ideadcrootxf_to_ideadcroot_ph)   = ideadcroot

    this%iretransn_to_ileaf_ph           = 19
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ileaf_ph)              = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ileaf_ph)           = ileaf
    
    this%iretransn_to_ileafst_ph         = 20
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ileafst_ph)            = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ileafst_ph)         = ileaf_st
    
    this%iretransn_to_ifroot_ph          = 21
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ifroot_ph)             = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ifroot_ph)          = ifroot
    
    this%iretransn_to_ifrootst_ph        = 22
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ifrootst_ph)           = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ifrootst_ph)        = ifroot_st
    
    this%iretransn_to_ilivestem_ph       = 23
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ilivestem_ph)          = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ilivestem_ph)       = ilivestem
    
    this%iretransn_to_ilivestemst_ph     = 24
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ilivestemst_ph)        = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ilivestemst_ph)     = ilivestem_st
    
    this%iretransn_to_ideadstem_ph       = 25
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ideadstem_ph)          = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ideadstem_ph)       = ideadstem
    
    this%iretransn_to_ideadstemst_ph     = 26
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ideadstemst_ph)        = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ideadstemst_ph)     = ideadstem_st
    
    this%iretransn_to_ilivecroot_ph      = 27
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ilivecroot_ph)         = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ilivecroot_ph)      = ilivecroot
    
    this%iretransn_to_ilivecrootst_ph    = 28
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ilivecrootst_ph)       = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ilivecrootst_ph)    = ilivecroot_st
    
    this%iretransn_to_ideadcroot_ph      = 29
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ideadcroot_ph)         = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ideadcroot_ph)      = ideadcroot
    
    this%iretransn_to_ideadcrootst_ph    = 30
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ideadcrootst_ph)       = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ideadcrootst_ph)    = ideadcroot_st
    
    if(.not. use_crop)then
       this%ileaf_to_iout_ph             = 31
       this%matrix_nphtransfer_doner_patch(this%ileaf_to_iout_ph)                = ileaf
       this%matrix_nphtransfer_receiver_patch(this%ileaf_to_iout_ph)             = ioutn
    
       this%ifroot_to_iout_ph            = 32
       this%matrix_nphtransfer_doner_patch(this%ifroot_to_iout_ph)               = ifroot
       this%matrix_nphtransfer_receiver_patch(this%ifroot_to_iout_ph)            = ioutn
    
       this%ilivestem_to_iout_ph         = 33
       this%matrix_nphtransfer_doner_patch(this%ilivestem_to_iout_ph)            = ilivestem
       this%matrix_nphtransfer_receiver_patch(this%ilivestem_to_iout_ph)         = ioutn
    
       this%iretransn_to_iout_ph         = 34
       this%matrix_nphtransfer_doner_patch(this%iretransn_to_iout_ph)            = iretransn
       this%matrix_nphtransfer_receiver_patch(this%iretransn_to_iout_ph)         = ioutn
    else
       this%iretransn_to_igrain_ph       = 31
       this%matrix_nphtransfer_doner_patch(this%iretransn_to_igrain_ph)          = iretransn
       this%matrix_nphtransfer_receiver_patch(this%iretransn_to_igrain_ph)       = igrain

       this%iretransn_to_igrainst_ph     = 32
       this%matrix_nphtransfer_doner_patch(this%iretransn_to_igrainst_ph)        = iretransn
       this%matrix_nphtransfer_receiver_patch(this%iretransn_to_igrainst_ph)     = igrain_st

       this%ileaf_to_iout_ph             = 33
       this%matrix_nphtransfer_doner_patch(this%ileaf_to_iout_ph)                = ileaf
       this%matrix_nphtransfer_receiver_patch(this%ileaf_to_iout_ph)             = ioutn
       
       this%ifroot_to_iout_ph            = 34
       this%matrix_nphtransfer_doner_patch(this%ifroot_to_iout_ph)               = ifroot
       this%matrix_nphtransfer_receiver_patch(this%ifroot_to_iout_ph)            = ioutn
    
       this%ilivestem_to_iout_ph         = 35
       this%matrix_nphtransfer_doner_patch(this%ilivestem_to_iout_ph)            = ilivestem
       this%matrix_nphtransfer_receiver_patch(this%ilivestem_to_iout_ph)         = ioutn

       this%igrain_to_iout_ph            = 36
       this%matrix_nphtransfer_doner_patch(this%igrain_to_iout_ph)               = igrain
       this%matrix_nphtransfer_receiver_patch(this%igrain_to_iout_ph)            = ioutn
    
       this%iretransn_to_iout_ph         = 37
       this%matrix_nphtransfer_doner_patch(this%iretransn_to_iout_ph)            = iretransn
       this%matrix_nphtransfer_receiver_patch(this%iretransn_to_iout_ph)         = ioutn
    end if

    this%ileaf_to_iout_gm                = 1
    this%matrix_ngmtransfer_doner_patch(this%ileaf_to_iout_gm)                   = ileaf 
    this%matrix_ngmtransfer_receiver_patch(this%ileaf_to_iout_gm)                = ioutn
    
    this%ileafst_to_iout_gm              = 2
    this%matrix_ngmtransfer_doner_patch(this%ileafst_to_iout_gm)                 = ileaf_st
    this%matrix_ngmtransfer_receiver_patch(this%ileafst_to_iout_gm)              = ioutn
    
    this%ileafxf_to_iout_gm              = 3
    this%matrix_ngmtransfer_doner_patch(this%ileafxf_to_iout_gm)                 = ileaf_xf
    this%matrix_ngmtransfer_receiver_patch(this%ileafxf_to_iout_gm)              = ioutn
    
    this%ifroot_to_iout_gm               = 4
    this%matrix_ngmtransfer_doner_patch(this%ifroot_to_iout_gm)                  = ifroot 
    this%matrix_ngmtransfer_receiver_patch(this%ifroot_to_iout_gm)               = ioutn
    
    this%ifrootst_to_iout_gm             = 5
    this%matrix_ngmtransfer_doner_patch(this%ifrootst_to_iout_gm)                = ifroot_st
    this%matrix_ngmtransfer_receiver_patch(this%ifrootst_to_iout_gm)             = ioutn
    
    this%ifrootxf_to_iout_gm             = 6
    this%matrix_ngmtransfer_doner_patch(this%ifrootxf_to_iout_gm)                = ifroot_xf
    this%matrix_ngmtransfer_receiver_patch(this%ifrootxf_to_iout_gm)             = ioutn
    
    this%ilivestem_to_iout_gm            = 7
    this%matrix_ngmtransfer_doner_patch(this%ilivestem_to_iout_gm)               = ilivestem 
    this%matrix_ngmtransfer_receiver_patch(this%ilivestem_to_iout_gm)            = ioutn
    
    this%ilivestemst_to_iout_gm          = 8
    this%matrix_ngmtransfer_doner_patch(this%ilivestemst_to_iout_gm)             = ilivestem_st
    this%matrix_ngmtransfer_receiver_patch(this%ilivestemst_to_iout_gm)          = ioutn
    
    this%ilivestemxf_to_iout_gm          = 9
    this%matrix_ngmtransfer_doner_patch(this%ilivestemxf_to_iout_gm)             = ilivestem_xf
    this%matrix_ngmtransfer_receiver_patch(this%ilivestemxf_to_iout_gm)          = ioutn
    
    this%ideadstem_to_iout_gm            = 10
    this%matrix_ngmtransfer_doner_patch(this%ideadstem_to_iout_gm)               = ideadstem 
    this%matrix_ngmtransfer_receiver_patch(this%ideadstem_to_iout_gm)            = ioutn
    
    this%ideadstemst_to_iout_gm          = 11
    this%matrix_ngmtransfer_doner_patch(this%ideadstemst_to_iout_gm)             = ideadstem_st
    this%matrix_ngmtransfer_receiver_patch(this%ideadstemst_to_iout_gm)          = ioutn
    
    this%ideadstemxf_to_iout_gm          = 12
    this%matrix_ngmtransfer_doner_patch(this%ideadstemxf_to_iout_gm)             = ideadstem_xf
    this%matrix_ngmtransfer_receiver_patch(this%ideadstemxf_to_iout_gm)          = ioutn
    
    this%ilivecroot_to_iout_gm           = 13
    this%matrix_ngmtransfer_doner_patch(this%ilivecroot_to_iout_gm)              = ilivecroot 
    this%matrix_ngmtransfer_receiver_patch(this%ilivecroot_to_iout_gm)           = ioutn
    
    this%ilivecrootst_to_iout_gm         = 14
    this%matrix_ngmtransfer_doner_patch(this%ilivecrootst_to_iout_gm)            = ilivecroot_st
    this%matrix_ngmtransfer_receiver_patch(this%ilivecrootst_to_iout_gm)         = ioutn
    
    this%ilivecrootxf_to_iout_gm         = 15
    this%matrix_ngmtransfer_doner_patch(this%ilivecrootxf_to_iout_gm)            = ilivecroot_xf
    this%matrix_ngmtransfer_receiver_patch(this%ilivecrootxf_to_iout_gm)         = ioutn
    
    this%ideadcroot_to_iout_gm           = 16
    this%matrix_ngmtransfer_doner_patch(this%ideadcroot_to_iout_gm)              = ideadcroot 
    this%matrix_ngmtransfer_receiver_patch(this%ideadcroot_to_iout_gm)           = ioutn
    
    this%ideadcrootst_to_iout_gm         = 17
    this%matrix_ngmtransfer_doner_patch(this%ideadcrootst_to_iout_gm)            = ideadcroot_st
    this%matrix_ngmtransfer_receiver_patch(this%ideadcrootst_to_iout_gm)         = ioutn
    
    this%ideadcrootxf_to_iout_gm         = 18
    this%matrix_ngmtransfer_doner_patch(this%ideadcrootxf_to_iout_gm)            = ideadcroot_xf
    this%matrix_ngmtransfer_receiver_patch(this%ideadcrootxf_to_iout_gm)         = ioutn
    
    this%iretransn_to_iout_gm            = 19
    this%matrix_ngmtransfer_doner_patch(this%iretransn_to_iout_gm)               = iretransn
    this%matrix_ngmtransfer_receiver_patch(this%iretransn_to_iout_gm)            = ioutn
    
    this%ilivestem_to_ideadstem_fi       = 1
    this%matrix_nfitransfer_doner_patch(this%ilivestem_to_ideadstem_fi)          = ilivestem 
    this%matrix_nfitransfer_receiver_patch(this%ilivestem_to_ideadstem_fi)       = ideadstem
    
    this%ilivecroot_to_ideadcroot_fi     = 2
    this%matrix_nfitransfer_doner_patch(this%ilivecroot_to_ideadcroot_fi)        = ilivecroot 
    this%matrix_nfitransfer_receiver_patch(this%ilivecroot_to_ideadcroot_fi)     = ideadcroot

    this%ileaf_to_iout_fi                = 3
    this%matrix_nfitransfer_doner_patch(this%ileaf_to_iout_fi)                   = ileaf 
    this%matrix_nfitransfer_receiver_patch(this%ileaf_to_iout_fi)                = ioutn
    
    this%ileafst_to_iout_fi              = 4
    this%matrix_nfitransfer_doner_patch(this%ileafst_to_iout_fi)                 = ileaf_st
    this%matrix_nfitransfer_receiver_patch(this%ileafst_to_iout_fi)              = ioutn
    
    this%ileafxf_to_iout_fi              = 5
    this%matrix_nfitransfer_doner_patch(this%ileafxf_to_iout_fi)                 = ileaf_xf
    this%matrix_nfitransfer_receiver_patch(this%ileafxf_to_iout_fi)              = ioutn
    
    this%ifroot_to_iout_fi               = 6
    this%matrix_nfitransfer_doner_patch(this%ifroot_to_iout_fi)                  = ifroot 
    this%matrix_nfitransfer_receiver_patch(this%ifroot_to_iout_fi)               = ioutn
    
    this%ifrootst_to_iout_fi             = 7
    this%matrix_nfitransfer_doner_patch(this%ifrootst_to_iout_fi)                = ifroot_st
    this%matrix_nfitransfer_receiver_patch(this%ifrootst_to_iout_fi)             = ioutn
    
    this%ifrootxf_to_iout_fi             = 8
    this%matrix_nfitransfer_doner_patch(this%ifrootxf_to_iout_fi)                = ifroot_xf
    this%matrix_nfitransfer_receiver_patch(this%ifrootxf_to_iout_fi)             = ioutn
    
    this%ilivestem_to_iout_fi            = 9
    this%matrix_nfitransfer_doner_patch(this%ilivestem_to_iout_fi)               = ilivestem 
    this%matrix_nfitransfer_receiver_patch(this%ilivestem_to_iout_fi)            = ioutn
    
    this%ilivestemst_to_iout_fi          = 10
    this%matrix_nfitransfer_doner_patch(this%ilivestemst_to_iout_fi)             = ilivestem_st
    this%matrix_nfitransfer_receiver_patch(this%ilivestemst_to_iout_fi)          = ioutn
    
    this%ilivestemxf_to_iout_fi          = 11
    this%matrix_nfitransfer_doner_patch(this%ilivestemxf_to_iout_fi)             = ilivestem_xf
    this%matrix_nfitransfer_receiver_patch(this%ilivestemxf_to_iout_fi)          = ioutn
    
    this%ideadstem_to_iout_fi            = 12
    this%matrix_nfitransfer_doner_patch(this%ideadstem_to_iout_fi)               = ideadstem 
    this%matrix_nfitransfer_receiver_patch(this%ideadstem_to_iout_fi)            = ioutn
    
    this%ideadstemst_to_iout_fi          = 13
    this%matrix_nfitransfer_doner_patch(this%ideadstemst_to_iout_fi)             = ideadstem_st
    this%matrix_nfitransfer_receiver_patch(this%ideadstemst_to_iout_fi)          = ioutn
    
    this%ideadstemxf_to_iout_fi          = 14
    this%matrix_nfitransfer_doner_patch(this%ideadstemxf_to_iout_fi)             = ideadstem_xf
    this%matrix_nfitransfer_receiver_patch(this%ideadstemxf_to_iout_fi)          = ioutn
    
    this%ilivecroot_to_iout_fi           = 15
    this%matrix_nfitransfer_doner_patch(this%ilivecroot_to_iout_fi)              = ilivecroot 
    this%matrix_nfitransfer_receiver_patch(this%ilivecroot_to_iout_fi)           = ioutn
    
    this%ilivecrootst_to_iout_fi         = 16
    this%matrix_nfitransfer_doner_patch(this%ilivecrootst_to_iout_fi)            = ilivecroot_st
    this%matrix_nfitransfer_receiver_patch(this%ilivecrootst_to_iout_fi)         = ioutn
    
    this%ilivecrootxf_to_iout_fi         = 17
    this%matrix_nfitransfer_doner_patch(this%ilivecrootxf_to_iout_fi)            = ilivecroot_xf
    this%matrix_nfitransfer_receiver_patch(this%ilivecrootxf_to_iout_fi)         = ioutn
    
    this%ideadcroot_to_iout_fi           = 18
    this%matrix_nfitransfer_doner_patch(this%ideadcroot_to_iout_fi)              = ideadcroot 
    this%matrix_nfitransfer_receiver_patch(this%ideadcroot_to_iout_fi)           = ioutn
    
    this%ideadcrootst_to_iout_fi         = 19
    this%matrix_nfitransfer_doner_patch(this%ideadcrootst_to_iout_fi)            = ideadcroot_st
    this%matrix_nfitransfer_receiver_patch(this%ideadcrootst_to_iout_fi)         = ioutn
    
    this%ideadcrootxf_to_iout_fi         = 20
    this%matrix_nfitransfer_doner_patch(this%ideadcrootxf_to_iout_fi)            = ideadcroot_xf
    this%matrix_nfitransfer_receiver_patch(this%ideadcrootxf_to_iout_fi)         = ioutn
    
    this%iretransn_to_iout_fi            = 21
    this%matrix_nfitransfer_doner_patch(this%iretransn_to_iout_fi)               = iretransn
    this%matrix_nfitransfer_receiver_patch(this%iretransn_to_iout_fi)            = ioutn
    
  end subroutine InitTransfer 

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, alloc_full_veg)
    !
    ! !DESCRIPTION:
    ! Initialize patch nitrogen flux
    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    logical,intent(in)             :: alloc_full_veg
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    if(alloc_full_veg)then
       begp = bounds%begp; endp = bounds%endp
       begc = bounds%begc; endc = bounds%endc
       begg = bounds%begg; endg = bounds%endg
    else
       begp = 0; endp = 0
       begc = 0; endc = 0
       begg = 0; endg = 0
    end if
       
    allocate(this%m_leafn_to_litter_patch                   (begp:endp)) ; this%m_leafn_to_litter_patch                   (:) = nan
    allocate(this%m_frootn_to_litter_patch                  (begp:endp)) ; this%m_frootn_to_litter_patch                  (:) = nan
    allocate(this%m_leafn_storage_to_litter_patch           (begp:endp)) ; this%m_leafn_storage_to_litter_patch           (:) = nan
    allocate(this%m_frootn_storage_to_litter_patch          (begp:endp)) ; this%m_frootn_storage_to_litter_patch          (:) = nan
    allocate(this%m_livestemn_storage_to_litter_patch       (begp:endp)) ; this%m_livestemn_storage_to_litter_patch       (:) = nan
    allocate(this%m_deadstemn_storage_to_litter_patch       (begp:endp)) ; this%m_deadstemn_storage_to_litter_patch       (:) = nan
    allocate(this%m_livecrootn_storage_to_litter_patch      (begp:endp)) ; this%m_livecrootn_storage_to_litter_patch      (:) = nan
    allocate(this%m_deadcrootn_storage_to_litter_patch      (begp:endp)) ; this%m_deadcrootn_storage_to_litter_patch      (:) = nan
    allocate(this%m_leafn_xfer_to_litter_patch              (begp:endp)) ; this%m_leafn_xfer_to_litter_patch              (:) = nan
    allocate(this%m_frootn_xfer_to_litter_patch             (begp:endp)) ; this%m_frootn_xfer_to_litter_patch             (:) = nan
    allocate(this%m_livestemn_xfer_to_litter_patch          (begp:endp)) ; this%m_livestemn_xfer_to_litter_patch          (:) = nan
    allocate(this%m_deadstemn_xfer_to_litter_patch          (begp:endp)) ; this%m_deadstemn_xfer_to_litter_patch          (:) = nan
    allocate(this%m_livecrootn_xfer_to_litter_patch         (begp:endp)) ; this%m_livecrootn_xfer_to_litter_patch         (:) = nan
    allocate(this%m_deadcrootn_xfer_to_litter_patch         (begp:endp)) ; this%m_deadcrootn_xfer_to_litter_patch         (:) = nan
    allocate(this%m_livestemn_to_litter_patch               (begp:endp)) ; this%m_livestemn_to_litter_patch               (:) = nan
    allocate(this%m_deadstemn_to_litter_patch               (begp:endp)) ; this%m_deadstemn_to_litter_patch               (:) = nan
    allocate(this%m_livecrootn_to_litter_patch              (begp:endp)) ; this%m_livecrootn_to_litter_patch              (:) = nan
    allocate(this%m_deadcrootn_to_litter_patch              (begp:endp)) ; this%m_deadcrootn_to_litter_patch              (:) = nan
    allocate(this%m_retransn_to_litter_patch                (begp:endp)) ; this%m_retransn_to_litter_patch                (:) = nan
    allocate(this%hrv_leafn_to_litter_patch                 (begp:endp)) ; this%hrv_leafn_to_litter_patch                 (:) = nan
    allocate(this%hrv_frootn_to_litter_patch                (begp:endp)) ; this%hrv_frootn_to_litter_patch                (:) = nan
    allocate(this%hrv_leafn_storage_to_litter_patch         (begp:endp)) ; this%hrv_leafn_storage_to_litter_patch         (:) = nan
    allocate(this%hrv_frootn_storage_to_litter_patch        (begp:endp)) ; this%hrv_frootn_storage_to_litter_patch        (:) = nan
    allocate(this%hrv_livestemn_storage_to_litter_patch     (begp:endp)) ; this%hrv_livestemn_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_deadstemn_storage_to_litter_patch     (begp:endp)) ; this%hrv_deadstemn_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_livecrootn_storage_to_litter_patch    (begp:endp)) ; this%hrv_livecrootn_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_deadcrootn_storage_to_litter_patch    (begp:endp)) ; this%hrv_deadcrootn_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_leafn_xfer_to_litter_patch            (begp:endp)) ; this%hrv_leafn_xfer_to_litter_patch            (:) = nan
    allocate(this%hrv_frootn_xfer_to_litter_patch           (begp:endp)) ; this%hrv_frootn_xfer_to_litter_patch           (:) = nan
    allocate(this%hrv_livestemn_xfer_to_litter_patch        (begp:endp)) ; this%hrv_livestemn_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_deadstemn_xfer_to_litter_patch        (begp:endp)) ; this%hrv_deadstemn_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_livecrootn_xfer_to_litter_patch       (begp:endp)) ; this%hrv_livecrootn_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_deadcrootn_xfer_to_litter_patch       (begp:endp)) ; this%hrv_deadcrootn_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_livestemn_to_litter_patch             (begp:endp)) ; this%hrv_livestemn_to_litter_patch             (:) = nan
    allocate(this%hrv_livecrootn_to_litter_patch            (begp:endp)) ; this%hrv_livecrootn_to_litter_patch            (:) = nan
    allocate(this%hrv_deadcrootn_to_litter_patch            (begp:endp)) ; this%hrv_deadcrootn_to_litter_patch            (:) = nan
    allocate(this%hrv_retransn_to_litter_patch              (begp:endp)) ; this%hrv_retransn_to_litter_patch              (:) = nan

    allocate(this%m_leafn_to_fire_patch                     (begp:endp)) ; this%m_leafn_to_fire_patch                     (:) = nan
    allocate(this%m_leafn_storage_to_fire_patch             (begp:endp)) ; this%m_leafn_storage_to_fire_patch             (:) = nan
    allocate(this%m_leafn_xfer_to_fire_patch                (begp:endp)) ; this%m_leafn_xfer_to_fire_patch                (:) = nan
    allocate(this%m_livestemn_to_fire_patch                 (begp:endp)) ; this%m_livestemn_to_fire_patch                 (:) = nan
    allocate(this%m_livestemn_storage_to_fire_patch         (begp:endp)) ; this%m_livestemn_storage_to_fire_patch         (:) = nan
    allocate(this%m_livestemn_xfer_to_fire_patch            (begp:endp)) ; this%m_livestemn_xfer_to_fire_patch            (:) = nan
    allocate(this%m_deadstemn_to_fire_patch                 (begp:endp)) ; this%m_deadstemn_to_fire_patch                 (:) = nan
    allocate(this%m_deadstemn_storage_to_fire_patch         (begp:endp)) ; this%m_deadstemn_storage_to_fire_patch         (:) = nan
    allocate(this%m_deadstemn_xfer_to_fire_patch            (begp:endp)) ; this%m_deadstemn_xfer_to_fire_patch            (:) = nan
    allocate(this%m_frootn_to_fire_patch                    (begp:endp)) ; this%m_frootn_to_fire_patch                    (:) = nan
    allocate(this%m_frootn_storage_to_fire_patch            (begp:endp)) ; this%m_frootn_storage_to_fire_patch            (:) = nan
    allocate(this%m_frootn_xfer_to_fire_patch               (begp:endp)) ; this%m_frootn_xfer_to_fire_patch               (:) = nan
    allocate(this%m_livecrootn_to_fire_patch                (begp:endp)) ;     
    allocate(this%m_livecrootn_storage_to_fire_patch        (begp:endp)) ; this%m_livecrootn_storage_to_fire_patch        (:) = nan
    allocate(this%m_livecrootn_xfer_to_fire_patch           (begp:endp)) ; this%m_livecrootn_xfer_to_fire_patch           (:) = nan
    allocate(this%m_deadcrootn_to_fire_patch                (begp:endp)) ; this%m_deadcrootn_to_fire_patch                (:) = nan
    allocate(this%m_deadcrootn_storage_to_fire_patch        (begp:endp)) ; this%m_deadcrootn_storage_to_fire_patch        (:) = nan
    allocate(this%m_deadcrootn_xfer_to_fire_patch           (begp:endp)) ; this%m_deadcrootn_xfer_to_fire_patch           (:) = nan
    allocate(this%m_retransn_to_fire_patch                  (begp:endp)) ; this%m_retransn_to_fire_patch                  (:) = nan

    allocate(this%m_leafn_to_litter_fire_patch              (begp:endp)) ; this%m_leafn_to_litter_fire_patch              (:) = nan
    allocate(this%m_leafn_storage_to_litter_fire_patch      (begp:endp)) ; this%m_leafn_storage_to_litter_fire_patch      (:) = nan
    allocate(this%m_leafn_xfer_to_litter_fire_patch         (begp:endp)) ; this%m_leafn_xfer_to_litter_fire_patch         (:) = nan
    allocate(this%m_livestemn_to_litter_fire_patch          (begp:endp)) ; this%m_livestemn_to_litter_fire_patch          (:) = nan
    allocate(this%m_livestemn_storage_to_litter_fire_patch  (begp:endp)) ; this%m_livestemn_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_livestemn_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_livestemn_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_livestemn_to_deadstemn_fire_patch       (begp:endp)) ; this%m_livestemn_to_deadstemn_fire_patch       (:) = nan
    allocate(this%m_deadstemn_to_litter_fire_patch          (begp:endp)) ; this%m_deadstemn_to_litter_fire_patch          (:) = nan
    allocate(this%m_deadstemn_storage_to_litter_fire_patch  (begp:endp)) ; this%m_deadstemn_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_deadstemn_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_deadstemn_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootn_to_litter_fire_patch             (begp:endp)) ; this%m_frootn_to_litter_fire_patch             (:) = nan
    allocate(this%m_frootn_storage_to_litter_fire_patch     (begp:endp)) ; this%m_frootn_storage_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootn_xfer_to_litter_fire_patch        (begp:endp)) ; this%m_frootn_xfer_to_litter_fire_patch        (:) = nan
    allocate(this%m_livecrootn_to_litter_fire_patch         (begp:endp)) ; this%m_livecrootn_to_litter_fire_patch         (:) = nan
    allocate(this%m_livecrootn_storage_to_litter_fire_patch (begp:endp)) ; this%m_livecrootn_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_livecrootn_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_livecrootn_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_livecrootn_to_deadcrootn_fire_patch     (begp:endp)) ; this%m_livecrootn_to_deadcrootn_fire_patch     (:) = nan
    allocate(this%m_deadcrootn_to_litter_fire_patch         (begp:endp)) ; this%m_deadcrootn_to_litter_fire_patch         (:) = nan
    allocate(this%m_deadcrootn_storage_to_litter_fire_patch (begp:endp)) ; this%m_deadcrootn_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_deadcrootn_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_deadcrootn_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_retransn_to_litter_fire_patch           (begp:endp)) ; this%m_retransn_to_litter_fire_patch           (:) = nan

    allocate(this%leafn_xfer_to_leafn_patch                 (begp:endp)) ; this%leafn_xfer_to_leafn_patch                 (:) = nan
    allocate(this%frootn_xfer_to_frootn_patch               (begp:endp)) ; this%frootn_xfer_to_frootn_patch               (:) = nan
    allocate(this%livestemn_xfer_to_livestemn_patch         (begp:endp)) ; this%livestemn_xfer_to_livestemn_patch         (:) = nan
    allocate(this%deadstemn_xfer_to_deadstemn_patch         (begp:endp)) ; this%deadstemn_xfer_to_deadstemn_patch         (:) = nan
    allocate(this%livecrootn_xfer_to_livecrootn_patch       (begp:endp)) ; this%livecrootn_xfer_to_livecrootn_patch       (:) = nan
    allocate(this%deadcrootn_xfer_to_deadcrootn_patch       (begp:endp)) ; this%deadcrootn_xfer_to_deadcrootn_patch       (:) = nan
    allocate(this%leafn_to_litter_patch                     (begp:endp)) ; this%leafn_to_litter_patch                     (:) = nan
    allocate(this%leafn_to_retransn_patch                   (begp:endp)) ; this%leafn_to_retransn_patch                   (:) = nan
    allocate(this%frootn_to_retransn_patch                  (begp:endp)) ; this%frootn_to_retransn_patch                  (:) = nan
    allocate(this%frootn_to_litter_patch                    (begp:endp)) ; this%frootn_to_litter_patch                    (:) = nan
    allocate(this%retransn_to_npool_patch                   (begp:endp)) ; this%retransn_to_npool_patch                   (:) = nan
    allocate(this%free_retransn_to_npool_patch              (begp:endp)) ; this%free_retransn_to_npool_patch              (:) = nan
    allocate(this%sminn_to_npool_patch                      (begp:endp)) ; this%sminn_to_npool_patch                      (:) = nan

    allocate(this%npool_to_leafn_patch                      (begp:endp)) ; this%npool_to_leafn_patch                      (:) = nan
    allocate(this%npool_to_leafn_storage_patch              (begp:endp)) ; this%npool_to_leafn_storage_patch              (:) = nan
    allocate(this%npool_to_frootn_patch                     (begp:endp)) ; this%npool_to_frootn_patch                     (:) = nan
    allocate(this%npool_to_frootn_storage_patch             (begp:endp)) ; this%npool_to_frootn_storage_patch             (:) = nan
    allocate(this%npool_to_livestemn_patch                  (begp:endp)) ; this%npool_to_livestemn_patch                  (:) = nan
    allocate(this%npool_to_livestemn_storage_patch          (begp:endp)) ; this%npool_to_livestemn_storage_patch          (:) = nan
    allocate(this%npool_to_deadstemn_patch                  (begp:endp)) ; this%npool_to_deadstemn_patch                  (:) = nan
    allocate(this%npool_to_deadstemn_storage_patch          (begp:endp)) ; this%npool_to_deadstemn_storage_patch          (:) = nan
    allocate(this%npool_to_livecrootn_patch                 (begp:endp)) ; this%npool_to_livecrootn_patch                 (:) = nan
    allocate(this%npool_to_livecrootn_storage_patch         (begp:endp)) ; this%npool_to_livecrootn_storage_patch         (:) = nan
    allocate(this%npool_to_deadcrootn_patch                 (begp:endp)) ; this%npool_to_deadcrootn_patch                 (:) = nan
    allocate(this%npool_to_deadcrootn_storage_patch         (begp:endp)) ; this%npool_to_deadcrootn_storage_patch         (:) = nan
    allocate(this%leafn_storage_to_xfer_patch               (begp:endp)) ; this%leafn_storage_to_xfer_patch               (:) = nan
    allocate(this%frootn_storage_to_xfer_patch              (begp:endp)) ; this%frootn_storage_to_xfer_patch              (:) = nan
    allocate(this%livestemn_storage_to_xfer_patch           (begp:endp)) ; this%livestemn_storage_to_xfer_patch           (:) = nan
    allocate(this%deadstemn_storage_to_xfer_patch           (begp:endp)) ; this%deadstemn_storage_to_xfer_patch           (:) = nan
    allocate(this%livecrootn_storage_to_xfer_patch          (begp:endp)) ; this%livecrootn_storage_to_xfer_patch          (:) = nan
    allocate(this%deadcrootn_storage_to_xfer_patch          (begp:endp)) ; this%deadcrootn_storage_to_xfer_patch          (:) = nan
    allocate(this%livestemn_to_deadstemn_patch              (begp:endp)) ; this%livestemn_to_deadstemn_patch              (:) = nan
    allocate(this%livestemn_to_retransn_patch               (begp:endp)) ; this%livestemn_to_retransn_patch               (:) = nan
    allocate(this%livecrootn_to_deadcrootn_patch            (begp:endp)) ; this%livecrootn_to_deadcrootn_patch            (:) = nan
    allocate(this%livecrootn_to_retransn_patch              (begp:endp)) ; this%livecrootn_to_retransn_patch              (:) = nan
    allocate(this%ndeploy_patch                             (begp:endp)) ; this%ndeploy_patch                             (:) = nan
    allocate(this%wood_harvestn_patch                       (begp:endp)) ; this%wood_harvestn_patch                       (:) = nan
    allocate(this%fire_nloss_patch                          (begp:endp)) ; this%fire_nloss_patch                          (:) = nan
    allocate(this%npool_to_reproductiven_patch       (begp:endp, nrepr)) ; this%npool_to_reproductiven_patch            (:,:) = nan
    allocate(this%npool_to_reproductiven_storage_patch(begp:endp, nrepr)); this%npool_to_reproductiven_storage_patch    (:,:) = nan
    allocate(this%livestemn_to_litter_patch                 (begp:endp)) ; this%livestemn_to_litter_patch                 (:) = nan
    allocate(this%repr_grainn_to_food_patch(begp:endp, repr_grain_min:repr_grain_max)) ; this%repr_grainn_to_food_patch (:,:) = nan
    allocate(this%repr_grainn_to_food_perharv_patch(begp:endp, 1:mxharvests, repr_grain_min:repr_grain_max)) ; this%repr_grainn_to_food_perharv_patch (:,:,:) = nan
    allocate(this%repr_grainn_to_food_thisyr_patch(begp:endp, repr_grain_min:repr_grain_max)) ; this%repr_grainn_to_food_thisyr_patch (:,:) = nan
    allocate(this%repr_structuren_to_cropprod_patch(begp:endp, repr_structure_min:repr_structure_max))
    this%repr_structuren_to_cropprod_patch(:,:) = nan
    allocate(this%repr_structuren_to_litter_patch(begp:endp, repr_structure_min:repr_structure_max))
    this%repr_structuren_to_litter_patch(:,:) = nan
    allocate(this%leafn_to_biofueln_patch                   (begp:endp)) ; this%leafn_to_biofueln_patch                   (:) = nan
    allocate(this%livestemn_to_biofueln_patch               (begp:endp)) ; this%livestemn_to_biofueln_patch               (:) = nan
    allocate(this%leafn_to_removedresiduen_patch            (begp:endp)) ; this%leafn_to_removedresiduen_patch            (:) = nan
    allocate(this%livestemn_to_removedresiduen_patch        (begp:endp)) ; this%livestemn_to_removedresiduen_patch        (:) = nan
    allocate(this%repr_grainn_to_seed_patch(begp:endp, repr_grain_min:repr_grain_max)) ; this%repr_grainn_to_seed_patch (:,:) = nan
    allocate(this%repr_grainn_to_seed_perharv_patch(begp:endp, 1:mxharvests, repr_grain_min:repr_grain_max)) ; this%repr_grainn_to_seed_perharv_patch (:,:,:) = nan
    allocate(this%repr_grainn_to_seed_thisyr_patch(begp:endp, repr_grain_min:repr_grain_max)) ; this%repr_grainn_to_seed_thisyr_patch (:,:) = nan
    allocate(this%reproductiven_xfer_to_reproductiven_patch(begp:endp, nrepr))
    this%reproductiven_xfer_to_reproductiven_patch(:,:) = nan
    allocate(this%reproductiven_storage_to_xfer_patch(begp:endp, nrepr)) ; this%reproductiven_storage_to_xfer_patch     (:,:) = nan
    allocate(this%fert_patch                                (begp:endp)) ; this%fert_patch                                (:) = nan
    allocate(this%fert_counter_patch                        (begp:endp)) ; this%fert_counter_patch                        (:) = nan
    allocate(this%soyfixn_patch                             (begp:endp)) ; this%soyfixn_patch                             (:) = nan

    allocate(this%crop_harvestn_to_cropprodn_patch                 (begp:endp)) ; this%crop_harvestn_to_cropprodn_patch                 (:) = nan
    allocate(this%crop_harvestn_to_cropprodn_col                   (begc:endc)) ; this%crop_harvestn_to_cropprodn_col                   (:) = nan

    allocate(this%fire_nloss_col                            (begc:endc)) ; this%fire_nloss_col                            (:) = nan
    allocate(this%fire_nloss_p2c_col                        (begc:endc)) ; this%fire_nloss_p2c_col                        (:) = nan

    allocate(this%m_n_to_litr_fire_col         (begc:endc,1:nlevdecomp_full,1:ndecomp_pools)) ; this%m_n_to_litr_fire_col     (:,:,:) = nan

    allocate(this%dwt_seedn_to_leaf_patch      (begp:endp))                   ; this%dwt_seedn_to_leaf_patch      (:)   = nan
    allocate(this%dwt_seedn_to_leaf_grc        (begg:endg))                   ; this%dwt_seedn_to_leaf_grc        (:)   = nan
    allocate(this%dwt_seedn_to_deadstem_patch  (begp:endp))                   ; this%dwt_seedn_to_deadstem_patch  (:)   = nan
    allocate(this%dwt_seedn_to_deadstem_grc    (begg:endg))                   ; this%dwt_seedn_to_deadstem_grc    (:)   = nan
    allocate(this%dwt_conv_nflux_patch         (begp:endp))                   ; this%dwt_conv_nflux_patch         (:)   = nan
    allocate(this%dwt_conv_nflux_grc           (begg:endg))                   ; this%dwt_conv_nflux_grc           (:)   = nan
    allocate(this%dwt_wood_productn_gain_patch (begp:endp))                   ; this%dwt_wood_productn_gain_patch (:)   = nan
    allocate(this%dwt_crop_productn_gain_patch (begp:endp))                   ; this%dwt_crop_productn_gain_patch (:)   = nan
    allocate(this%wood_harvestn_col            (begc:endc))                   ; this%wood_harvestn_col            (:)   = nan

    allocate(this%dwt_frootn_to_litr_n_col     (begc:endc,1:nlevdecomp_full,1:i_litr_max)) ; this%dwt_frootn_to_litr_n_col (:,:,:) = nan
    allocate(this%dwt_livecrootn_to_cwdn_col   (begc:endc,1:nlevdecomp_full)) ; this%dwt_livecrootn_to_cwdn_col   (:,:) = nan
    allocate(this%dwt_deadcrootn_to_cwdn_col   (begc:endc,1:nlevdecomp_full)) ; this%dwt_deadcrootn_to_cwdn_col   (:,:) = nan

    allocate(this%gru_leafn_to_litter_patch                 (begp:endp)) ; this%gru_leafn_to_litter_patch                 (:) = nan
    allocate(this%gru_leafn_storage_to_atm_patch            (begp:endp)) ; this%gru_leafn_storage_to_atm_patch            (:) = nan
    allocate(this%gru_leafn_xfer_to_atm_patch               (begp:endp)) ; this%gru_leafn_xfer_to_atm_patch               (:) = nan
    allocate(this%gru_frootn_to_litter_patch                (begp:endp)) ; this%gru_frootn_to_litter_patch                (:) = nan
    allocate(this%gru_frootn_storage_to_atm_patch           (begp:endp)) ; this%gru_frootn_storage_to_atm_patch           (:) = nan
    allocate(this%gru_frootn_xfer_to_atm_patch              (begp:endp)) ; this%gru_frootn_xfer_to_atm_patch              (:) = nan
    allocate(this%gru_livestemn_to_atm_patch                (begp:endp)) ; this%gru_livestemn_to_atm_patch                (:) = nan
    allocate(this%gru_livestemn_storage_to_atm_patch        (begp:endp)) ; this%gru_livestemn_storage_to_atm_patch        (:) = nan
    allocate(this%gru_livestemn_xfer_to_atm_patch           (begp:endp)) ; this%gru_livestemn_xfer_to_atm_patch           (:) = nan
    allocate(this%gru_deadstemn_to_atm_patch                (begp:endp)) ; this%gru_deadstemn_to_atm_patch                (:) = nan
    allocate(this%gru_deadstemn_storage_to_atm_patch        (begp:endp)) ; this%gru_deadstemn_storage_to_atm_patch        (:) = nan
    allocate(this%gru_deadstemn_xfer_to_atm_patch           (begp:endp)) ; this%gru_deadstemn_xfer_to_atm_patch           (:) = nan
    allocate(this%gru_livecrootn_storage_to_atm_patch       (begp:endp)) ; this%gru_livecrootn_storage_to_atm_patch       (:) = nan
    allocate(this%gru_livecrootn_xfer_to_atm_patch          (begp:endp)) ; this%gru_livecrootn_xfer_to_atm_patch          (:) = nan
    allocate(this%gru_livecrootn_to_litter_patch            (begp:endp)) ; this%gru_livecrootn_to_litter_patch            (:) = nan
    allocate(this%gru_deadcrootn_storage_to_atm_patch       (begp:endp)) ; this%gru_deadcrootn_storage_to_atm_patch       (:) = nan
    allocate(this%gru_deadcrootn_xfer_to_atm_patch          (begp:endp)) ; this%gru_deadcrootn_xfer_to_atm_patch          (:) = nan
    allocate(this%gru_deadcrootn_to_litter_patch            (begp:endp)) ; this%gru_deadcrootn_to_litter_patch            (:) = nan
    allocate(this%gru_retransn_to_litter_patch              (begp:endp)) ; this%gru_retransn_to_litter_patch              (:) = nan

    allocate(this%gru_conv_nflux_patch                      (begp:endp))                  ; this%gru_conv_nflux_patch         (:)  =nan
    allocate(this%gru_conv_nflux_col                        (begc:endc))                  ; this%gru_conv_nflux_col           (:)  =nan
    allocate(this%gru_conv_nflux_grc                        (begg:endg))                  ; this%gru_conv_nflux_grc           (:)  =nan
    allocate(this%gru_wood_productn_gain_patch              (begp:endp))                  ; this%gru_wood_productn_gain_patch (:)  =nan
    allocate(this%gru_wood_productn_gain_col                (begc:endc))                  ; this%gru_wood_productn_gain_col   (:)  =nan
    allocate(this%gru_wood_productn_gain_grc                (begg:endg))                  ; this%gru_wood_productn_gain_grc   (:)  =nan
    allocate(this%gru_n_to_litr_n_col                       (begc:endc,1:nlevdecomp_full,1:i_litr_max)); this%gru_n_to_litr_n_col(:,:,:)=nan
    allocate(this%gru_n_to_cwdn_col                         (begc:endc,1:nlevdecomp_full)); this%gru_n_to_cwdn_col            (:,:)=nan

    allocate(this%crop_seedn_to_leaf_patch     (begp:endp))                   ; this%crop_seedn_to_leaf_patch     (:)   = nan

    allocate(this%m_decomp_npools_to_fire_vr_col    (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(this%m_decomp_npools_to_fire_col       (begc:endc,1:ndecomp_pools                  ))

    this%m_decomp_npools_to_fire_vr_col   (:,:,:) = nan
    this%m_decomp_npools_to_fire_col      (:,:)   = nan

    allocate(this%phenology_n_to_litr_n_col         (begc:endc, 1:nlevdecomp_full, 1:ndecomp_pools))
    allocate(this%gap_mortality_n_to_litr_n_col     (begc:endc, 1:nlevdecomp_full, 1:ndecomp_pools))
    allocate(this%gap_mortality_n_to_cwdn_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%fire_mortality_n_to_cwdn_col      (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_litr_n_col           (begc:endc, 1:nlevdecomp_full, 1:ndecomp_pools))
    allocate(this%harvest_n_to_cwdn_col             (begc:endc, 1:nlevdecomp_full))

    this%phenology_n_to_litr_n_col       (:,:,:) = nan
    this%gap_mortality_n_to_litr_n_col   (:,:,:) = nan
    this%gap_mortality_n_to_cwdn_col       (:,:) = nan
    this%fire_mortality_n_to_cwdn_col      (:,:) = nan
    this%harvest_n_to_litr_n_col         (:,:,:) = nan
    this%harvest_n_to_cwdn_col             (:,:) = nan

    allocate(this%plant_ndemand_patch         (begp:endp)) ;    this%plant_ndemand_patch         (:) = nan
    allocate(this%avail_retransn_patch        (begp:endp)) ;    this%avail_retransn_patch        (:) = nan
    allocate(this%plant_nalloc_patch          (begp:endp)) ;    this%plant_nalloc_patch          (:) = nan

    allocate(this%plant_ndemand_retrans_patch (begp:endp)) ;    this%plant_ndemand_retrans_patch (:) = nan
    allocate(this%plant_ndemand_season_patch  (begp:endp)) ;    this%plant_ndemand_season_patch  (:) = nan
    allocate(this%plant_ndemand_stress_patch  (begp:endp)) ;    this%plant_ndemand_stress_patch  (:) = nan
    allocate(this%Nactive_patch               (begp:endp)) ;    this%Nactive_patch               (:) = nan 
    allocate(this%Nnonmyc_patch               (begp:endp)) ;    this%Nnonmyc_patch               (:) = nan
    allocate(this%Nam_patch                   (begp:endp)) ;    this%Nam_patch                   (:) = nan
    allocate(this%Necm_patch                  (begp:endp)) ;    this%Necm_patch                  (:) = nan
    allocate(this%Nactive_no3_patch           (begp:endp)) ;    this%Nactive_no3_patch           (:) = nan
    allocate(this%Nactive_nh4_patch           (begp:endp)) ;    this%Nactive_nh4_patch           (:) = nan
    allocate(this%Nnonmyc_no3_patch           (begp:endp)) ;    this%Nnonmyc_no3_patch           (:) = nan
    allocate(this%Nnonmyc_nh4_patch           (begp:endp)) ;    this%Nnonmyc_nh4_patch           (:) = nan
    allocate(this%Nam_no3_patch               (begp:endp)) ;    this%Nam_no3_patch               (:) = nan
    allocate(this%Nam_nh4_patch               (begp:endp)) ;    this%Nam_nh4_patch               (:) = nan
    allocate(this%Necm_no3_patch              (begp:endp)) ;    this%Necm_no3_patch              (:) = nan
    allocate(this%Necm_nh4_patch              (begp:endp)) ;    this%Necm_nh4_patch              (:) = nan
    allocate(this%Npassive_patch              (begp:endp)) ;    this%Npassive_patch              (:) = nan
    allocate(this%Nfix_patch                  (begp:endp)) ;    this%Nfix_patch                  (:) = nan
    allocate(this%Nretrans_patch              (begp:endp)) ;    this%Nretrans_patch              (:) = nan
    allocate(this%Nretrans_org_patch          (begp:endp)) ;    this%Nretrans_org_patch          (:) = nan
    allocate(this%Nretrans_season_patch       (begp:endp)) ;    this%Nretrans_season_patch       (:) = nan
    allocate(this%Nretrans_stress_patch       (begp:endp)) ;    this%Nretrans_stress_patch       (:) = nan 
    allocate(this%Nuptake_patch               (begp:endp)) ;    this%Nuptake_patch               (:) = nan
    allocate(this%sminn_to_plant_fun_patch    (begp:endp)) ;    this%sminn_to_plant_fun_patch    (:) = nan
    allocate(this%sminn_to_plant_fun_vr_patch (begp:endp,1:nlevdecomp_full)) 
    this%sminn_to_plant_fun_vr_patch          (:,:) = nan
    allocate(this%sminn_to_plant_fun_no3_vr_patch (begp:endp,1:nlevdecomp_full))  
    this%sminn_to_plant_fun_no3_vr_patch      (:,:) = nan
    allocate(this%sminn_to_plant_fun_nh4_vr_patch (begp:endp,1:nlevdecomp_full))  
    this%sminn_to_plant_fun_nh4_vr_patch      (:,:) = nan
    allocate(this%cost_nfix_patch              (begp:endp)) ;    this%cost_nfix_patch            (:) = nan
    allocate(this%cost_nactive_patch           (begp:endp)) ;    this%cost_nactive_patch         (:) = nan
    allocate(this%cost_nretrans_patch          (begp:endp)) ;    this%cost_nretrans_patch        (:) = nan
    allocate(this%nuptake_npp_fraction_patch   (begp:endp)) ;    this%nuptake_npp_fraction_patch            (:) = nan
	! Matrix
    if(use_matrixcn)then
       allocate(this%matrix_Ninput_patch               (begp:endp))               ; this%matrix_Ninput_patch              (:)   =  nan
       allocate(this%matrix_nalloc_patch               (begp:endp,1:nvegnpool))   ; this%matrix_nalloc_patch              (:,:) =  nan

       allocate(this%matrix_nphtransfer_patch          (begp:endp,1:nnphtrans))   ; this%matrix_nphtransfer_patch         (:,:) =  nan
       allocate(this%matrix_nphturnover_patch          (begp:endp,1:nvegnpool))   ; this%matrix_nphturnover_patch         (:,:) =  nan
       allocate(this%matrix_nphtransfer_doner_patch    (1:nnphtrans))             ; this%matrix_nphtransfer_doner_patch   (:)   = -9999
       allocate(this%matrix_nphtransfer_receiver_patch (1:nnphtrans))             ; this%matrix_nphtransfer_receiver_patch(:)   = -9999

       allocate(this%matrix_ngmtransfer_patch          (begp:endp,1:nngmtrans))   ; this%matrix_ngmtransfer_patch         (:,:) =  nan
       allocate(this%matrix_ngmturnover_patch          (begp:endp,1:nvegnpool))   ; this%matrix_ngmturnover_patch         (:,:) =  nan
       allocate(this%matrix_ngmtransfer_doner_patch    (1:nngmtrans))             ; this%matrix_ngmtransfer_doner_patch   (:)   = -9999
       allocate(this%matrix_ngmtransfer_receiver_patch (1:nngmtrans))             ; this%matrix_ngmtransfer_receiver_patch(:)   = -9999

       allocate(this%matrix_nfitransfer_patch          (begp:endp,1:nnfitrans))   ; this%matrix_nfitransfer_patch         (:,:) =  nan
       allocate(this%matrix_nfiturnover_patch          (begp:endp,1:nvegnpool))   ; this%matrix_nfiturnover_patch         (:,:) =  nan
       allocate(this%matrix_nfitransfer_doner_patch    (1:nnfitrans))             ; this%matrix_nfitransfer_doner_patch   (:)   = -9999
       allocate(this%matrix_nfitransfer_receiver_patch (1:nnfitrans))             ; this%matrix_nfitransfer_receiver_patch(:)   = -9999

       allocate(this%list_phn_phgmn                    (1:nnphtrans+nvegnpool))   ; this%list_phn_phgmn   = -9999
       allocate(this%list_gmn_phgmn                    (1:nvegnpool))             ; this%list_gmn_phgmn   = -9999
       allocate(this%list_phn_phgmfin                  (1:nnphtrans+nvegnpool))   ; this%list_phn_phgmfin = -9999
       allocate(this%list_gmn_phgmfin                  (1:nvegnpool))             ; this%list_gmn_phgmfin = -9999
       allocate(this%list_fin_phgmfin                  (1:nnfitrans+nvegnpool))   ; this%list_fin_phgmfin = -9999

       allocate(this%list_aphn                         (1:nnphtrans-nnphouttrans)); this%list_aphn        = -9999
       allocate(this%list_agmn                         (1:nngmtrans-nngmouttrans)); this%list_agmn        = -9999
       allocate(this%list_afin                         (1:nnfitrans-nnfiouttrans)); this%list_afin        = -9999

       call this%AKphvegn%InitSM  (nvegnpool,begp,endp,nnphtrans-nnphouttrans+nvegnpool)
       call this%AKgmvegn%InitSM  (nvegnpool,begp,endp,nngmtrans-nngmouttrans+nvegnpool)
       call this%AKfivegn%InitSM  (nvegnpool,begp,endp,nnfitrans-nnfiouttrans+nvegnpool)

       this%NE_AKallvegn = (nnphtrans-nnphouttrans+nvegnpool) + (nngmtrans-nngmouttrans+nvegnpool) + &
                           nnfitrans-nnfiouttrans+nvegnpool

       call this%AKallvegn%InitSM (nvegnpool,begp,endp,this%NE_AKallvegn)

       allocate(this%RI_AKallvegn                      (1:this%NE_AKallvegn))     ; this%RI_AKallvegn(:) = -9999
       allocate(this%CI_AKallvegn                      (1:this%NE_AKallvegn))     ; this%CI_AKallvegn(:) = -9999
       allocate(this%RI_phn             (1:nnphtrans-nnphouttrans+nvegnpool))     ; this%RI_phn(:)       = -9999
       allocate(this%CI_phn             (1:nnphtrans-nnphouttrans+nvegnpool))     ; this%CI_phn(:)       = -9999
       allocate(this%RI_gmn             (1:nngmtrans-nngmouttrans+nvegnpool))     ; this%RI_gmn(:)       = -9999
       allocate(this%CI_gmn             (1:nngmtrans-nngmouttrans+nvegnpool))     ; this%CI_gmn(:)       = -9999
       allocate(this%RI_fin             (1:nnfitrans-nnfiouttrans+nvegnpool))     ; this%RI_fin(:)       = -9999
       allocate(this%CI_fin             (1:nnfitrans-nnfiouttrans+nvegnpool))     ; this%CI_fin(:)       = -9999

       call this%Kvegn%InitDM (nvegnpool,begp,endp)
       call this%Xvegn%InitV  (nvegnpool,begp,endp)
     end if

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer        :: k,l
    integer        :: begp, endp
    integer        :: begc, endc
    integer        :: begg, endg
    character(10)  :: active
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    this%m_leafn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N mortality', &
         ptr_patch=this%m_leafn_to_litter_patch, default='inactive')

    this%m_frootn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N mortality', &
         ptr_patch=this%m_frootn_to_litter_patch, default='inactive')

    this%m_leafn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage mortality', &
         ptr_patch=this%m_leafn_storage_to_litter_patch, default='inactive')

    this%m_frootn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage mortality', &
         ptr_patch=this%m_frootn_storage_to_litter_patch, default='inactive')

    this%m_livestemn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage mortality', &
         ptr_patch=this%m_livestemn_storage_to_litter_patch, default='inactive')

    this%m_deadstemn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage mortality', &
         ptr_patch=this%m_deadstemn_storage_to_litter_patch, default='inactive')

    this%m_livecrootn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage mortality', &
         ptr_patch=this%m_livecrootn_storage_to_litter_patch, default='inactive')

    this%m_deadcrootn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage mortality', &
         ptr_patch=this%m_deadcrootn_storage_to_litter_patch, default='inactive')

    this%m_leafn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer mortality', &
         ptr_patch=this%m_leafn_xfer_to_litter_patch, default='inactive')

    this%m_frootn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer mortality', &
         ptr_patch=this%m_frootn_xfer_to_litter_patch, default='inactive')

    this%m_livestemn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer mortality', &
         ptr_patch=this%m_livestemn_xfer_to_litter_patch, default='inactive')

    this%m_deadstemn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer mortality', &
         ptr_patch=this%m_deadstemn_xfer_to_litter_patch, default='inactive')

    this%m_livecrootn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer mortality', &
         ptr_patch=this%m_livecrootn_xfer_to_litter_patch, default='inactive')

    this%m_deadcrootn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer mortality', &
         ptr_patch=this%m_deadcrootn_xfer_to_litter_patch, default='inactive')

    this%m_livestemn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N mortality', &
         ptr_patch=this%m_livestemn_to_litter_patch, default='inactive')

    this%m_deadstemn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N mortality', &
         ptr_patch=this%m_deadstemn_to_litter_patch, default='inactive')

    this%m_livecrootn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N mortality', &
         ptr_patch=this%m_livecrootn_to_litter_patch, default='inactive')

    this%m_deadcrootn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N mortality', &
         ptr_patch=this%m_deadcrootn_to_litter_patch, default='inactive')

    this%m_retransn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool mortality', &
         ptr_patch=this%m_retransn_to_litter_patch, default='inactive')

    this%m_leafn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N fire loss', &
         ptr_patch=this%m_leafn_to_fire_patch, default='inactive')

    this%m_frootn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N fire loss ', &
         ptr_patch=this%m_frootn_to_fire_patch, default='inactive')

    this%m_leafn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage fire loss', &
         ptr_patch=this%m_leafn_storage_to_fire_patch, default='inactive')

    this%m_frootn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage fire loss', &
         ptr_patch=this%m_frootn_storage_to_fire_patch, default='inactive')

    this%m_livestemn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage fire loss', &
         ptr_patch=this%m_livestemn_storage_to_fire_patch, default='inactive')

    this%m_deadstemn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage fire loss', &
         ptr_patch=this%m_deadstemn_storage_to_fire_patch, default='inactive')

    this%m_livecrootn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage fire loss', &
         ptr_patch=this%m_livecrootn_storage_to_fire_patch, default='inactive')

    this%m_deadcrootn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage fire loss', &
         ptr_patch=this%m_deadcrootn_storage_to_fire_patch, default='inactive')

    this%m_leafn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer fire loss', &
         ptr_patch=this%m_leafn_xfer_to_fire_patch, default='inactive')

    this%m_frootn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer fire loss', &
         ptr_patch=this%m_frootn_xfer_to_fire_patch, default='inactive')

    this%m_livestemn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer fire loss', &
         ptr_patch=this%m_livestemn_xfer_to_fire_patch, default='inactive')

    this%m_deadstemn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer fire loss', &
         ptr_patch=this%m_deadstemn_xfer_to_fire_patch, default='inactive')

    this%m_livecrootn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer fire loss', &
         ptr_patch=this%m_livecrootn_xfer_to_fire_patch, default='inactive')

    this%m_deadcrootn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer fire loss', &
         ptr_patch=this%m_deadcrootn_xfer_to_fire_patch, default='inactive')

    this%m_livestemn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N fire loss', &
         ptr_patch=this%m_livestemn_to_fire_patch, default='inactive')

    this%m_deadstemn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire loss', &
         ptr_patch=this%m_deadstemn_to_fire_patch, default='inactive')

    this%m_deadstemn_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire mortality to litter', &
         ptr_patch=this%m_deadstemn_to_litter_fire_patch, default='inactive')

    this%m_livecrootn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N fire loss', &
         ptr_patch=this%m_livecrootn_to_fire_patch, default='inactive')

    this%m_deadcrootn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire loss', &
         ptr_patch=this%m_deadcrootn_to_fire_patch, default='inactive')

    this%m_deadcrootn_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire mortality to litter', &
         ptr_patch=this%m_deadcrootn_to_litter_fire_patch, default='inactive')

    this%m_retransn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool fire loss', &
         ptr_patch=this%m_retransn_to_fire_patch, default='inactive')

    this%leafn_xfer_to_leafn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_XFER_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N growth from storage', &
         ptr_patch=this%leafn_xfer_to_leafn_patch, default='inactive')

    this%frootn_xfer_to_frootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_XFER_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N growth from storage', &
         ptr_patch=this%frootn_xfer_to_frootn_patch, default='inactive')

    this%livestemn_xfer_to_livestemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_XFER_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N growth from storage', &
         ptr_patch=this%livestemn_xfer_to_livestemn_patch, default='inactive')

    this%deadstemn_xfer_to_deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_XFER_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N growth from storage', &
         ptr_patch=this%deadstemn_xfer_to_deadstemn_patch, default='inactive')

    this%livecrootn_xfer_to_livecrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_XFER_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N growth from storage', &
         ptr_patch=this%livecrootn_xfer_to_livecrootn_patch, default='inactive')

    this%deadcrootn_xfer_to_deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_XFER_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N growth from storage', &
         ptr_patch=this%deadcrootn_xfer_to_deadcrootn_patch, default='inactive')

    this%leafn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N litterfall', &
         ptr_patch=this%leafn_to_litter_patch)

    this%leafn_to_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N to retranslocated N pool', &
         ptr_patch=this%leafn_to_retransn_patch, default='inactive')

    this%frootn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N litterfall', &
         ptr_patch=this%frootn_to_litter_patch, default='inactive')

    this%retransn_to_npool_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of retranslocated N', &
         ptr_patch=this%retransn_to_npool_patch)
         
    this%free_retransn_to_npool_patch(begp:endp) = spval
    call hist_addfld1d (fname='FREE_RETRANSN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of retranslocated N', &
         ptr_patch=this%free_retransn_to_npool_patch)

    this%sminn_to_npool_patch(begp:endp) = spval
    call hist_addfld1d (fname='SMINN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of soil mineral N uptake', &
         ptr_patch=this%sminn_to_npool_patch)

    this%npool_to_leafn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N', &
         ptr_patch=this%npool_to_leafn_patch, default='inactive')

    this%npool_to_leafn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LEAFN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N storage', &
         ptr_patch=this%npool_to_leafn_storage_patch, default='inactive')

    this%npool_to_frootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N', &
         ptr_patch=this%npool_to_frootn_patch, default='inactive')

    this%npool_to_frootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_FROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N storage', &
         ptr_patch=this%npool_to_frootn_storage_patch, default='inactive')

    this%npool_to_livestemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N', &
         ptr_patch=this%npool_to_livestemn_patch, default='inactive')

    this%npool_to_livestemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N storage', &
         ptr_patch=this%npool_to_livestemn_storage_patch, default='inactive')

    this%npool_to_deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N', &
         ptr_patch=this%npool_to_deadstemn_patch, default='inactive')

    this%npool_to_deadstemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N storage', &
         ptr_patch=this%npool_to_deadstemn_storage_patch, default='inactive')

    this%npool_to_livecrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N', &
         ptr_patch=this%npool_to_livecrootn_patch, default='inactive')

    this%npool_to_livecrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N storage', &
         ptr_patch=this%npool_to_livecrootn_storage_patch, default='inactive')

    this%npool_to_deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N', &
         ptr_patch=this%npool_to_deadcrootn_patch, default='inactive')

    this%npool_to_deadcrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N storage', &
         ptr_patch=this%npool_to_deadcrootn_storage_patch, default='inactive')

    this%leafn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N shift storage to transfer', &
         ptr_patch=this%leafn_storage_to_xfer_patch, default='inactive')

    this%frootn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N shift storage to transfer', &
         ptr_patch=this%frootn_storage_to_xfer_patch, default='inactive')

    this%livestemn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N shift storage to transfer', &
         ptr_patch=this%livestemn_storage_to_xfer_patch, default='inactive')

    this%deadstemn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N shift storage to transfer', &
         ptr_patch=this%deadstemn_storage_to_xfer_patch, default='inactive')

    this%livecrootn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N shift storage to transfer', &
         ptr_patch=this%livecrootn_storage_to_xfer_patch, default='inactive')

    this%deadcrootn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N shift storage to transfer', &
         ptr_patch=this%deadcrootn_storage_to_xfer_patch, default='inactive')

    this%livestemn_to_deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N turnover', &
         ptr_patch=this%livestemn_to_deadstemn_patch, default='inactive')

    this%livestemn_to_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N to retranslocated N pool', &
         ptr_patch=this%livestemn_to_retransn_patch, default='inactive')

    this%livecrootn_to_deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N turnover', &
         ptr_patch=this%livecrootn_to_deadcrootn_patch, default='inactive')

    this%livecrootn_to_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N to retranslocated N pool', &
         ptr_patch=this%livecrootn_to_retransn_patch, default='inactive')

    this%ndeploy_patch(begp:endp) = spval
    call hist_addfld1d (fname='NDEPLOY', units='gN/m^2/s', &
         avgflag='A', long_name='total N deployed in new growth', &
         ptr_patch=this%ndeploy_patch)

    this%wood_harvestn_patch(begp:endp) = spval
    call hist_addfld1d (fname='WOOD_HARVESTN', units='gN/m^2/s', &
         avgflag='A', long_name='wood harvest N (to product pools)', &
         ptr_patch=this%wood_harvestn_patch)

    this%fire_nloss_patch(begp:endp) = spval
    call hist_addfld1d (fname='PFT_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total patch-level fire N loss', &
         ptr_patch=this%fire_nloss_patch)

    if (use_crop) then
       this%fert_patch(begp:endp) = spval
       call hist_addfld1d (fname='NFERTILIZATION', units='gN/m^2/s', &
            avgflag='A', long_name='fertilizer added', &
            ptr_patch=this%fert_patch)
       
       this%repr_grainn_to_food_patch(begp:endp,:) = spval
       this%repr_grainn_to_food_perharv_patch(begp:endp,:,:) = spval
       this%repr_grainn_to_food_thisyr_patch(begp:endp,:) = spval
       this%repr_grainn_to_seed_patch(begp:endp,:) = spval
       this%repr_grainn_to_seed_perharv_patch(begp:endp,:,:) = spval
       this%repr_grainn_to_seed_thisyr_patch(begp:endp,:) = spval
       do k = repr_grain_min, repr_grain_max
          data1dptr => this%repr_grainn_to_food_patch(:,k)
          call hist_addfld1d ( &
               ! e.g., GRAINN_TO_FOOD
               fname=get_repr_hist_fname(k)//'N_TO_FOOD', &
               units='gN/m^2/s', &
               avgflag='A', &
               long_name=get_repr_longname(k)//' N to food (not scientifically supported)', &
               ptr_patch=data1dptr, &
               default='inactive')
          data1dptr => this%repr_grainn_to_seed_patch(:,k)
          call hist_addfld1d ( &
               ! e.g., GRAINN_TO_SEED
               fname=get_repr_hist_fname(k)//'N_TO_SEED', &
               units='gN/m^2/s', &
               avgflag='A', &
               long_name=get_repr_longname(k)//' N to seed (not scientifically supported)', &
               ptr_patch=data1dptr, &
               default='inactive')
          data2dptr => this%repr_grainn_to_food_perharv_patch(:,:,k)
          call hist_addfld2d ( &
               ! e.g., GRAINN_TO_FOOD_PERHARV
               fname=get_repr_hist_fname(k)//'N_TO_FOOD_PERHARV', &
               units='gN/m^2', &
               type2d='mxharvests', &
               avgflag='I', &
               long_name=get_repr_longname(k)//' N to food per harvest; should only be output annually (not scientifically supported)', &
               ptr_patch=data2dptr, &
               default='inactive')
          data1dptr => this%repr_grainn_to_food_thisyr_patch(:,k)
          call hist_addfld1d ( &
               ! e.g., GRAINN_TO_FOOD_ANN
               fname=get_repr_hist_fname(k)//'N_TO_FOOD_ANN', &
               units='gN/m^2', &
               avgflag='I', &
               long_name=get_repr_longname(k)//' N to food harvested per calendar year; should only be output annually (not scientifically supported)', &
               ptr_patch=data1dptr, &
               default='inactive')
          data2dptr => this%repr_grainn_to_seed_perharv_patch(:,:,k)
          call hist_addfld2d ( &
               ! e.g., GRAINN_TO_SEED_PERHARV
               fname=get_repr_hist_fname(k)//'N_TO_SEED_PERHARV', &
               units='gN/m^2', &
               type2d='mxharvests', &
               avgflag='I', &
               long_name=get_repr_longname(k)//' N to seed per harvest; should only be output annually (not scientifically supported)', &
               ptr_patch=data2dptr, &
               default='inactive')
          data1dptr => this%repr_grainn_to_seed_thisyr_patch(:,k)
          call hist_addfld1d ( &
               ! e.g., GRAINN_TO_SEED_ANN
               fname=get_repr_hist_fname(k)//'N_TO_SEED_ANN', &
               units='gN/m^2', &
               avgflag='I', &
               long_name=get_repr_longname(k)//' N to seed harvested per calendar year; should only be output annually (not scientifically supported)', &
               ptr_patch=data1dptr, &
               default='inactive')
       end do
    end if

    if (use_crop .and. .not. use_fun) then
       this%soyfixn_patch(begp:endp) = spval
       call hist_addfld1d (fname='SOYFIXN', units='gN/m^2/s', &
            avgflag='A', long_name='soybean fixation', &
            ptr_patch=this%soyfixn_patch)
    end if

    if (use_crop) then
       this%fert_counter_patch(begp:endp) = spval
       call hist_addfld1d (fname='FERT_COUNTER', units='seconds', &
            avgflag='A', long_name='time left to fertilize', &
            ptr_patch=this%fert_counter_patch, default='inactive')
    end if

    !-------------------------------
    ! N flux variables - native to column
    !-------------------------------

    do k = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
          this%m_decomp_npools_to_fire_col(begc:endc,k) = spval
          data1dptr => this%m_decomp_npools_to_fire_col(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'_N_TO_FIRE'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N fire loss'
          call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          if ( nlevdecomp_full > 1 ) then
             this%m_decomp_npools_to_fire_vr_col(begc:endc,:,k) = spval
             data2dptr => this%m_decomp_npools_to_fire_vr_col(:,:,k)
             fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'_N_TO_FIRE'//trim(vr_suffix)
             longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N fire loss'
             call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       endif
    end do

    this%fire_nloss_col(begc:endc) = spval
    call hist_addfld1d (fname='COL_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total column-level fire N loss', &
         ptr_col=this%fire_nloss_col)

    this%dwt_seedn_to_leaf_grc(begg:endg) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_LEAF', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to patch-level leaf', &
         ptr_gcell=this%dwt_seedn_to_leaf_grc)

    this%dwt_seedn_to_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_LEAF_PATCH', units='gN/m^2/s', &
         avgflag='A', &
         long_name='patch-level seed source to patch-level leaf ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_seedn_to_leaf_patch, default='inactive')

    this%dwt_seedn_to_deadstem_grc(begg:endg) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_DEADSTEM', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to patch-level deadstem', &
         ptr_gcell=this%dwt_seedn_to_deadstem_grc)

    this%dwt_seedn_to_deadstem_patch(begp:endp) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_DEADSTEM_PATCH', units='gN/m^2/s', &
         avgflag='A', &
         long_name='patch-level seed source to patch-level deadstem ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_seedn_to_deadstem_patch, default='inactive')

    this%dwt_conv_nflux_grc(begg:endg) = spval
    call hist_addfld1d (fname='DWT_CONV_NFLUX', units='gN/m^2/s', &
         avgflag='A', &
         long_name='conversion N flux (immediate loss to atm) (0 at all times except first timestep of year)', &
         ptr_gcell=this%dwt_conv_nflux_grc)

    this%dwt_conv_nflux_patch(begp:endp) = spval
    call hist_addfld1d (fname='DWT_CONV_NFLUX_PATCH', units='gN/m^2/s', &
         avgflag='A', &
         long_name='patch-level conversion N flux (immediate loss to atm) ' // &
         '(0 at all times except first timestep of year) ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=this%dwt_conv_nflux_patch, default='inactive')
    
    do k = i_litr_min, i_litr_max
       this%dwt_frootn_to_litr_n_col(begc:endc,:,k) = spval
       data2dptr => this%dwt_frootn_to_litr_n_col(begc:endc,:,k)
       fieldname = 'DWT_FROOTN_TO_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'_N'
       longname =  'fine root N to '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' due to landcover change'
       call hist_addfld_decomp (fname=fieldname, units='gN/m^2/s',  type2d='levdcmp', &
            avgflag='A', long_name=longname, &
            ptr_col=data2dptr, default='inactive')
    end do

    this%dwt_livecrootn_to_cwdn_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_LIVECROOTN_TO_CWDN', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=this%dwt_livecrootn_to_cwdn_col, default='inactive')

    this%dwt_deadcrootn_to_cwdn_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_DEADCROOTN_TO_CWDN', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=this%dwt_deadcrootn_to_cwdn_col, default='inactive')

    if ( get_do_grossunrep() )then
       this%gru_conv_nflux_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRU_CONV_NFLUX', units='gN/m^2/s', &
            avgflag='A', long_name='gross unrepresented conversion N flux (immediate loss to atm) (0 at all times except first timestep of year)', &
            ptr_patch=this%gru_conv_nflux_patch)

       this%gru_wood_productn_gain_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRU_WOODPRODN_GAIN', units='gN/m^2/s', &
            avgflag='A', long_name='gross unrepresented landcover change driven addition to wood product nitrogen pools (0 at all times except first timestep of year)', &
            ptr_patch=this%gru_wood_productn_gain_patch)
    end if

    this%crop_seedn_to_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname='CROP_SEEDN_TO_LEAF', units='gN/m^2/s', &
         avgflag='A', long_name='crop seed source to leaf', &
         ptr_patch=this%crop_seedn_to_leaf_patch, default='inactive')

    this%plant_ndemand_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_NDEMAND', units='gN/m^2/s', &
         avgflag='A', long_name='N flux required to support initial GPP', &
         ptr_patch=this%plant_ndemand_patch)

    this%avail_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='AVAIL_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='N flux available from retranslocation pool', &
         ptr_patch=this%avail_retransn_patch, default='inactive')

    this%plant_nalloc_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_NALLOC', units='gN/m^2/s', &
         avgflag='A', long_name='total allocated N flux', &
         ptr_patch=this%plant_nalloc_patch, default='inactive')
    if (use_matrixcn) then		 
       this%matrix_Ninput_patch(begp:endp) = spval
       call hist_addfld1d (fname='MATRIX PLANT_NALLOC', units='gN/m^2/s', &
            avgflag='A', long_name='total allocated N flux for matrix', &
            ptr_patch=this%matrix_Ninput_patch, default='inactive')
    end if 
    
    if ( use_fun ) then
       this%Nactive_patch(begp:endp)  = spval
       call hist_addfld1d (fname='NACTIVE', units='gN/m^2/s',       &
            avgflag='A', long_name='Mycorrhizal N uptake flux',     &
            ptr_patch=this%Nactive_patch)
   
       this%Nnonmyc_patch(begp:endp)  = spval
       call hist_addfld1d (fname='NNONMYC', units='gN/m^2/s',       &
            avgflag='A', long_name='Non-mycorrhizal N uptake flux', &
            ptr_patch=this%Nnonmyc_patch)    

       this%Nam_patch(begp:endp)      = spval
       call hist_addfld1d (fname='NAM', units='gN/m^2/s',           &
            avgflag='A', long_name='AM-associated N uptake flux',   &
            ptr_patch=this%Nam_patch)

       this%Necm_patch(begp:endp)     = spval
       call hist_addfld1d (fname='NECM', units='gN/m^2/s',          &
            avgflag='A', long_name='ECM-associated N uptake flux',  &
            ptr_patch=this%Necm_patch)

       if (use_nitrif_denitrif) then
          this%Nactive_no3_patch(begp:endp)  = spval
          call hist_addfld1d (fname='NACTIVE_NO3', units='gN/m^2/s',   &
               avgflag='A', long_name='Mycorrhizal N uptake flux',     &
               ptr_patch=this%Nactive_no3_patch)
   
          this%Nactive_nh4_patch(begp:endp)  = spval
          call hist_addfld1d (fname='NACTIVE_NH4', units='gN/m^2/s',   &
               avgflag='A', long_name='Mycorrhizal N uptake flux',     &
               ptr_patch=this%Nactive_nh4_patch)

          this%Nnonmyc_no3_patch(begp:endp)  = spval
          call hist_addfld1d (fname='NNONMYC_NO3', units='gN/m^2/s',   &
               avgflag='A', long_name='Non-mycorrhizal N uptake flux', &
               ptr_patch=this%Nnonmyc_no3_patch)
  
          this%Nnonmyc_nh4_patch(begp:endp)  = spval
          call hist_addfld1d (fname='NNONMYC_NH4', units='gN/m^2/s',   &
               avgflag='A', long_name='Non-mycorrhizal N uptake flux', &
               ptr_patch=this%Nnonmyc_nh4_patch)

          this%Nam_no3_patch(begp:endp)      = spval
          call hist_addfld1d (fname='NAM_NO3', units='gN/m^2/s',       &
               avgflag='A', long_name='AM-associated N uptake flux',   &
               ptr_patch=this%Nam_no3_patch)
 
          this%Nam_nh4_patch(begp:endp)      = spval
          call hist_addfld1d (fname='NAM_NH4', units='gN/m^2/s',       &
               avgflag='A', long_name='AM-associated N uptake flux',   &
               ptr_patch=this%Nam_nh4_patch)

          this%Necm_no3_patch(begp:endp)     = spval
          call hist_addfld1d (fname='NECM_NO3', units='gN/m^2/s',      &
               avgflag='A', long_name='ECM-associated N uptake flux',  &
               ptr_patch=this%Necm_no3_patch) 
  
          this%Necm_nh4_patch(begp:endp)     = spval
          call hist_addfld1d (fname='NECM_NH4', units='gN/m^2/s',      &
               avgflag='A', long_name='ECM-associated N uptake flux',  &
               ptr_patch=this%Necm_nh4_patch)   
       end if 

       this%Npassive_patch(begp:endp) = spval
       call hist_addfld1d (fname='NPASSIVE', units='gN/m^2/s',        &
            avgflag='A', long_name='Passive N uptake flux',           &
            ptr_patch=this%Npassive_patch)

       this%Nfix_patch(begp:endp)     = spval
       call hist_addfld1d (fname='NFIX', units='gN/m^2/s',            &
            avgflag='A', long_name='Symbiotic BNF uptake flux',       &
            ptr_patch=this%Nfix_patch)

       this%Nretrans_patch(begp:endp) = spval
       call hist_addfld1d (fname='NRETRANS', units='gN/m^2/s',        &
            avgflag='A', long_name='Retranslocated N uptake flux',    &
            ptr_patch=this%Nretrans_patch)
  
       this%Nretrans_org_patch(begp:endp) = spval
       call hist_addfld1d (fname='NRETRANS_REG', units='gN/m^2/s',    &
            avgflag='A', long_name='Retranslocated N uptake flux',    &
            ptr_patch=this%Nretrans_org_patch)

       this%Nretrans_season_patch(begp:endp) = spval
       call hist_addfld1d (fname='NRETRANS_SEASON', units='gN/m^2/s', &
            avgflag='A', long_name='Retranslocated N uptake flux',    &
            ptr_patch=this%Nretrans_season_patch)

       this%Nretrans_stress_patch(begp:endp) = spval
       call hist_addfld1d (fname='NRETRANS_STRESS', units='gN/m^2/s', &
            avgflag='A', long_name='Retranslocated N uptake flux',    &
            ptr_patch=this%Nretrans_stress_patch)
  
       this%Nuptake_patch(begp:endp) = spval
       call hist_addfld1d (fname='NUPTAKE', units='gN/m^2/s',         &
            avgflag='A', long_name='Total N uptake of FUN',           &
            ptr_patch=this%Nuptake_patch)

       this%sminn_to_plant_fun_patch(begp:endp) = spval
       call hist_addfld1d (fname='SMINN_TO_PLANT_FUN', units='gN/m^2/s',&
            avgflag='A', long_name='Total soil N uptake of FUN',        &
            ptr_patch=this%sminn_to_plant_fun_patch)      
       
       this%cost_nfix_patch(begp:endp)     = spval
       call hist_addfld1d (fname='COST_NFIX', units='gN/gC',            &
            avgflag='A', long_name='Cost of fixation',       &
            ptr_patch=this%cost_nfix_patch)
            
       this%cost_nactive_patch(begp:endp)     = spval
       call hist_addfld1d (fname='COST_NACTIVE', units='gN/gC',            &
            avgflag='A', long_name='Cost of active uptake',       &
            ptr_patch=this%cost_nactive_patch)
            
       this%cost_nretrans_patch(begp:endp)     = spval
       call hist_addfld1d (fname='COST_NRETRANS', units='gN/gC',            &
            avgflag='A', long_name='Cost of retranslocation',       &
            ptr_patch=this%cost_nretrans_patch)
            
      this%nuptake_npp_fraction_patch(begp:endp)     = spval
       call hist_addfld1d (fname='NUPTAKE_NPP_FRACTION', units='-',            &
            avgflag='A', long_name='frac of NPP used in N uptake',       &
            ptr_patch=this%nuptake_npp_fraction_patch)
                  
     
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenflux_type) :: this 
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,j
    integer :: fp, fc                                    ! filter indices
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !---------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    !-----------------------------------------------
    ! initialize nitrogen flux variables
    !-----------------------------------------------

    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)

       if ( use_crop )then
          this%fert_counter_patch(p)  = spval
          this%fert_patch(p)          = 0._r8 
          this%soyfixn_patch(p)       = 0._r8 
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%fert_counter_patch(p)  = 0._r8
       end if
       if ( use_fun ) then      !previously set to spval for special land units
          if (lun%ifspecial(l)) then
             this%plant_ndemand_patch(p)        = 0._r8
             this%avail_retransn_patch(p)       = 0._r8
             this%plant_nalloc_patch(p)         = 0._r8
             this%Npassive_patch(p)             = 0._r8
             this%Nactive_patch(p)              = 0._r8
             this%Nnonmyc_patch(p)              = 0._r8
             this%Nam_patch(p)                  = 0._r8
             this%Necm_patch(p)                 = 0._r8
             if (use_nitrif_denitrif) then
                this%Nactive_no3_patch(p)       = 0._r8
                this%Nactive_nh4_patch(p)       = 0._r8
                this%Nnonmyc_no3_patch(p)       = 0._r8
                this%Nnonmyc_nh4_patch(p)       = 0._r8
                this%Nam_no3_patch(p)           = 0._r8
                this%Nam_nh4_patch(p)           = 0._r8
                this%Necm_no3_patch(p)          = 0._r8
                this%Necm_nh4_patch(p)          = 0._r8
             end if
             this%Nfix_patch(p)                 = 0._r8
             this%Nretrans_patch(p)             = 0._r8
             this%Nretrans_org_patch(p)         = 0._r8
             this%Nretrans_season_patch(p)      = 0._r8
             this%Nretrans_stress_patch(p)      = 0._r8
             this%Nuptake_patch(p)              = 0._r8
             this%sminn_to_plant_fun_patch(p)   = 0._r8
             this%cost_nfix_patch               = 0._r8
             this%cost_nactive_patch            = 0._r8
             this%cost_nretrans_patch           = 0._r8          
             this%nuptake_npp_fraction_patch    = 0._r8
                             
             do j = 1, nlevdecomp
                this%sminn_to_plant_fun_vr_patch(p,j)       = 0._r8
                this%sminn_to_plant_fun_no3_vr_patch(p,j)   = 0._r8
                this%sminn_to_plant_fun_nh4_vr_patch(p,j)   = 0._r8
             end do 
          end if
       end if
    end do

    ! initialize fields for special filters

    call this%SetValues (nvegnpool=nvegnpool, &
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart (this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenflux_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,k,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=256) :: varname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !------------------------------------------------------------------------

    if (use_crop) then
       call restartvar(ncid=ncid, flag=flag, varname='fert_counter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_counter_patch)

       call restartvar(ncid=ncid, flag=flag, varname='fert', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_patch)
    end if

    if (use_crop) then
       do k = 1, nrepr
          data1dptr => this%reproductiven_xfer_to_reproductiven_patch(:,k)
          ! e.g., grainn_xfer_to_grainn
          varname = get_repr_rest_fname(k)//'n_xfer_to_'//get_repr_rest_fname(k)//'n'
          call restartvar(ncid=ncid, flag=flag,  varname=varname, &
               xtype=ncd_double,  &
               dim1name='pft', &
               long_name=get_repr_longname(k)//' N growth from storage', &
               units='gN/m2/s', &
               interpinic_flag='interp', readvar=readvar, data=data1dptr)
       end do

       ! Read or write variable(s) with mxharvests dimension
       ! BACKWARDS_COMPATIBILITY(wjs/ssr, 2022-06-10) See note in CallRestartvarDimOK()
       if (CallRestartvarDimOK(ncid, flag, 'mxharvests')) then
            do k = repr_grain_min, repr_grain_max
               data2dptr => this%repr_grainn_to_food_perharv_patch(:,:,k)
               ! e.g., grainn_to_food_perharv
               varname = get_repr_rest_fname(k)//'n_to_food_perharv'
               call restartvar(ncid=ncid, flag=flag,  varname=varname, &
                    xtype=ncd_double,  &
                    dim1name='pft', &
                    dim2name='mxharvests', &
                    switchdim=.true., &
                    long_name=get_repr_longname(k)//' N to food per harvest; should only be output annually', &
                    units='gN/m2', &
                    readvar=readvar, &
                    scale_by_thickness=.false., &
                    interpinic_flag='interp', data=data2dptr)
               data2dptr => this%repr_grainn_to_seed_perharv_patch(:,:,k)
               ! e.g., grainn_to_seed_perharv
               varname = get_repr_rest_fname(k)//'n_to_seed_perharv'
               call restartvar(ncid=ncid, flag=flag,  varname=varname, &
                    xtype=ncd_double,  &
                    dim1name='pft', &
                    dim2name='mxharvests', &
                    switchdim=.true., &
                    long_name=get_repr_longname(k)//' N to seed per harvest; should only be output annually', &
                    units='gN/m2', &
                    readvar=readvar, &
                    scale_by_thickness=.false., &
                    interpinic_flag='interp', data=data2dptr)
            end do
       end if
  
       do k = repr_grain_min, repr_grain_max
            data1dptr => this%repr_grainn_to_food_thisyr_patch(:,k)
            ! e.g., grainn_to_food_thisyr
            varname = get_repr_rest_fname(k)//'n_to_food_thisyr'
            call restartvar(ncid=ncid, flag=flag,  varname=varname, &
                 xtype=ncd_double,  &
                 dim1name='pft', &
                 long_name=get_repr_longname(k)//' N to food per calendar year; should only be output annually', &
                 units='gN/m2', &
                 interpinic_flag='interp', readvar=readvar, data=data1dptr)
            data1dptr => this%repr_grainn_to_seed_thisyr_patch(:,k)
            ! e.g., grainn_to_seed_thisyr
            varname = get_repr_rest_fname(k)//'n_to_seed_thisyr'
            call restartvar(ncid=ncid, flag=flag,  varname=varname, &
                 xtype=ncd_double,  &
                 dim1name='pft', &
                 long_name=get_repr_longname(k)//' N to seed per calendar year; should only be output annually', &
                 units='gN/m2', &
                 interpinic_flag='interp', readvar=readvar, data=data1dptr)
       end do
    end if

    if (use_crop) then
       call restartvar(ncid=ncid, flag=flag,  varname='livestemn_to_litter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='livestem N to litter', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemn_to_litter_patch)
    end if

    if (use_crop) then
       do k = repr_grain_min, repr_grain_max
          data1dptr => this%repr_grainn_to_food_patch(:,k)
          ! e.g., grainn_to_food
          varname = get_repr_rest_fname(k)//'n_to_food'
          call restartvar(ncid=ncid, flag=flag,  varname=varname, &
               xtype=ncd_double,  &
               dim1name='pft', &
               long_name=get_repr_longname(k)//' N to food', &
               units='gN/m2/s', &
               interpinic_flag='interp', readvar=readvar, data=data1dptr)
       end do
    end if
    
    if (use_crop) then
       do k = 1, nrepr
          data1dptr => this%npool_to_reproductiven_patch(:,k)
          ! e.g., npool_to_grainn
          varname = 'npool_to_'//get_repr_rest_fname(k)//'n'
          call restartvar(ncid=ncid, flag=flag,  varname=varname, &
               xtype=ncd_double,  &
               dim1name='pft', &
               long_name='allocation to '//get_repr_longname(k)//' N', &
               units='gN/m2/s', &
               interpinic_flag='interp', readvar=readvar, data=data1dptr)
       end do
    end if

    if (use_crop) then
       do k = 1, nrepr
          data1dptr => this%npool_to_reproductiven_storage_patch(:,k)
          ! e.g., npool_to_grainn_storage
          varname = 'npool_to_'//get_repr_rest_fname(k)//'n_storage'
          call restartvar(ncid=ncid, flag=flag,  varname=varname, &
               xtype=ncd_double,  &
               dim1name='pft', &
               long_name='allocation to '//get_repr_longname(k)//' N storage', &
               units='gN/m2/s', &
               interpinic_flag='interp', readvar=readvar, data=data1dptr)
       end do
    end if

    if (use_crop) then
       do k = 1, nrepr
          data1dptr => this%reproductiven_storage_to_xfer_patch(:,k)
          ! e.g., grainn_storage_to_xfer
          varname = get_repr_rest_fname(k)//'n_storage_to_xfer'
          call restartvar(ncid=ncid, flag=flag, varname=varname, &
               xtype=ncd_double,  &
               dim1name='pft', &
               long_name=get_repr_longname(k)//' N shift storage to transfer', &
               units='gN/m2/s', &
               interpinic_flag='interp', readvar=readvar, data=data1dptr)
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, varname='plant_ndemand', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_ndemand_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='avail_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%avail_retransn_patch) 

     if ( use_fun ) then
!       set_missing_vals_to_constant for BACKWARDS_COMPATIBILITY(wrw, 2018-06-28) re. issue #426
!       special land units previously set to spval, not 0
!       modifications here should correct this 
        call restartvar(ncid=ncid, flag=flag, varname='Nactive', xtype=ncd_double,       &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nactive_patch) 
        call set_missing_vals_to_constant(this%Nactive_patch, 0._r8)
    
        call restartvar(ncid=ncid, flag=flag, varname='Nnonmyc', xtype=ncd_double,       &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nnonmyc_patch)
        call set_missing_vals_to_constant(this%Nnonmyc_patch, 0._r8)
 
        call restartvar(ncid=ncid, flag=flag, varname='Nam', xtype=ncd_double,           &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nam_patch)
        call set_missing_vals_to_constant(this%Nam_patch, 0._r8)
 
        call restartvar(ncid=ncid, flag=flag, varname='Necm', xtype=ncd_double,          &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Necm_patch)
        call set_missing_vals_to_constant(this%Necm_patch, 0._r8)
 
        if (use_nitrif_denitrif) then
           call restartvar(ncid=ncid, flag=flag, varname='Nactive_no3', xtype=ncd_double,   &
                dim1name='pft', &
                long_name='', units='', &
                interpinic_flag='interp', readvar=readvar, data=this%Nactive_no3_patch)
           call set_missing_vals_to_constant(this%Nactive_no3_patch, 0._r8)
    
           call restartvar(ncid=ncid, flag=flag, varname='Nactive_nh4', xtype=ncd_double,   &
                dim1name='pft', &
                long_name='', units='', &
                interpinic_flag='interp', readvar=readvar, data=this%Nactive_nh4_patch)   
           call set_missing_vals_to_constant(this%Nactive_nh4_patch, 0._r8)
 
           call restartvar(ncid=ncid, flag=flag, varname='Nnonmyc_no3', xtype=ncd_double,    &
                dim1name='pft', &
                long_name='', units='', &
                interpinic_flag='interp', readvar=readvar, data=this%Nnonmyc_no3_patch)
           call set_missing_vals_to_constant(this%Nnonmyc_no3_patch, 0._r8)
 
           call restartvar(ncid=ncid, flag=flag, varname='Nnonmyc_nh4', xtype=ncd_double,    &
                dim1name='pft', &
                long_name='', units='', &
                interpinic_flag='interp', readvar=readvar, data=this%Nnonmyc_nh4_patch)
           call set_missing_vals_to_constant(this%Nnonmyc_nh4_patch, 0._r8)
 
           call restartvar(ncid=ncid, flag=flag, varname='Nam_no3', xtype=ncd_double,         &
                dim1name='pft', &
                long_name='', units='', &
                interpinic_flag='interp', readvar=readvar, data=this%Nam_no3_patch)
           call set_missing_vals_to_constant(this%Nam_no3_patch, 0._r8)
 
           call restartvar(ncid=ncid, flag=flag, varname='Nam_nh4', xtype=ncd_double,         &
                dim1name='pft', &
                long_name='', units='', &
                interpinic_flag='interp', readvar=readvar, data=this%Nam_nh4_patch)
           call set_missing_vals_to_constant(this%Nam_nh4_patch,  0._r8)
 
           call restartvar(ncid=ncid, flag=flag, varname='Necm_no3', xtype=ncd_double,        &
                dim1name='pft', &
                long_name='', units='', &
                interpinic_flag='interp', readvar=readvar, data=this%Necm_no3_patch)
           call set_missing_vals_to_constant(this%Necm_no3_patch, 0._r8)
 
           call restartvar(ncid=ncid, flag=flag, varname='Necm_nh4', xtype=ncd_double,        &
                dim1name='pft', &
                long_name='', units='', &
                interpinic_flag='interp', readvar=readvar, data=this%Necm_nh4_patch)
           call set_missing_vals_to_constant(this%Necm_nh4_patch, 0._r8)
        end if

        call restartvar(ncid=ncid, flag=flag, varname='Npassive', xtype=ncd_double,      &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Npassive_patch)
        call set_missing_vals_to_constant(this%Npassive_patch, 0._r8)
 
        call restartvar(ncid=ncid, flag=flag, varname='Nfix', xtype=ncd_double,          &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nfix_patch)
        call set_missing_vals_to_constant(this%Nfix_patch, 0._r8)
 
        call restartvar(ncid=ncid, flag=flag, varname='Nretrans', xtype=ncd_double,       &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nretrans_patch)
        call set_missing_vals_to_constant(this%Nretrans_patch, 0._r8)
 
        call restartvar(ncid=ncid, flag=flag, varname='Nretrans_org', xtype=ncd_double,   &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nretrans_org_patch)
        call set_missing_vals_to_constant(this%Nretrans_org_patch, 0._r8)
 
        call restartvar(ncid=ncid, flag=flag, varname='Nretrans_season', xtype=ncd_double, &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nretrans_season_patch)
        call set_missing_vals_to_constant(this%Nretrans_season_patch, 0._r8)
 
        call restartvar(ncid=ncid, flag=flag, varname='Nretrans_stress', xtype=ncd_double, &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nretrans_stress_patch)
        call set_missing_vals_to_constant(this%Nretrans_stress_patch, 0._r8)
 
        call restartvar(ncid=ncid, flag=flag, varname='Nuptake', xtype=ncd_double,            &
             dim1name='pft', &
             long_name='', units='', &
             interpinic_flag='interp', readvar=readvar, data=this%Nuptake_patch)
        call set_missing_vals_to_constant(this%Nuptake_patch, 0._r8)
     end if
! End BACKWARDS_COMPATIBILITY(wrw, 2018-06-28) re. issue #426

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this,nvegnpool, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class (cnveg_nitrogenflux_type) :: this
    integer , intent(in) :: num_patch,nvegnpool
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i=filter_patch(fi)

       this%m_leafn_to_litter_patch(i)                   = value_patch
       this%m_frootn_to_litter_patch(i)                  = value_patch
       this%m_leafn_storage_to_litter_patch(i)           = value_patch
       this%m_frootn_storage_to_litter_patch(i)          = value_patch
       this%m_livestemn_storage_to_litter_patch(i)       = value_patch
       this%m_deadstemn_storage_to_litter_patch(i)       = value_patch
       this%m_livecrootn_storage_to_litter_patch(i)      = value_patch
       this%m_deadcrootn_storage_to_litter_patch(i)      = value_patch
       this%m_leafn_xfer_to_litter_patch(i)              = value_patch
       this%m_frootn_xfer_to_litter_patch(i)             = value_patch
       this%m_livestemn_xfer_to_litter_patch(i)          = value_patch
       this%m_deadstemn_xfer_to_litter_patch(i)          = value_patch
       this%m_livecrootn_xfer_to_litter_patch(i)         = value_patch
       this%m_deadcrootn_xfer_to_litter_patch(i)         = value_patch
       this%m_livestemn_to_litter_patch(i)               = value_patch
       this%m_deadstemn_to_litter_patch(i)               = value_patch
       this%m_livecrootn_to_litter_patch(i)              = value_patch
       this%m_deadcrootn_to_litter_patch(i)              = value_patch
       this%m_retransn_to_litter_patch(i)                = value_patch
       this%hrv_leafn_to_litter_patch(i)                 = value_patch             
       this%hrv_frootn_to_litter_patch(i)                = value_patch            
       this%hrv_leafn_storage_to_litter_patch(i)         = value_patch     
       this%hrv_frootn_storage_to_litter_patch(i)        = value_patch    
       this%hrv_livestemn_storage_to_litter_patch(i)     = value_patch 
       this%hrv_deadstemn_storage_to_litter_patch(i)     = value_patch 
       this%hrv_livecrootn_storage_to_litter_patch(i)    = value_patch
       this%hrv_deadcrootn_storage_to_litter_patch(i)    = value_patch
       this%hrv_leafn_xfer_to_litter_patch(i)            = value_patch        
       this%hrv_frootn_xfer_to_litter_patch(i)           = value_patch       
       this%hrv_livestemn_xfer_to_litter_patch(i)        = value_patch    
       this%hrv_deadstemn_xfer_to_litter_patch(i)        = value_patch    
       this%hrv_livecrootn_xfer_to_litter_patch(i)       = value_patch   
       this%hrv_deadcrootn_xfer_to_litter_patch(i)       = value_patch   
       this%hrv_livestemn_to_litter_patch(i)             = value_patch         
       this%hrv_livecrootn_to_litter_patch(i)            = value_patch        
       this%hrv_deadcrootn_to_litter_patch(i)            = value_patch        
       this%hrv_retransn_to_litter_patch(i)              = value_patch    

       this%gru_leafn_to_litter_patch(i)                 = value_patch             
       this%gru_leafn_storage_to_atm_patch(i)            = value_patch     
       this%gru_leafn_xfer_to_atm_patch(i)               = value_patch        
       this%gru_frootn_to_litter_patch(i)                = value_patch            
       this%gru_frootn_storage_to_atm_patch(i)           = value_patch    
       this%gru_frootn_xfer_to_atm_patch(i)              = value_patch       
       this%gru_livestemn_to_atm_patch(i)                = value_patch         
       this%gru_livestemn_storage_to_atm_patch(i)        = value_patch 
       this%gru_livestemn_xfer_to_atm_patch(i)           = value_patch    
       this%gru_deadstemn_to_atm_patch(i)                = value_patch 
       this%gru_deadstemn_storage_to_atm_patch(i)        = value_patch 
       this%gru_deadstemn_xfer_to_atm_patch(i)           = value_patch    
       this%gru_livecrootn_to_litter_patch(i)            = value_patch        
       this%gru_livecrootn_storage_to_atm_patch(i)       = value_patch
       this%gru_livecrootn_xfer_to_atm_patch(i)          = value_patch   
       this%gru_deadcrootn_storage_to_atm_patch(i)       = value_patch
       this%gru_deadcrootn_xfer_to_atm_patch(i)          = value_patch   
       this%gru_deadcrootn_to_litter_patch(i)            = value_patch        
       this%gru_retransn_to_litter_patch(i)              = value_patch    

       this%m_leafn_to_fire_patch(i)                     = value_patch
       this%m_leafn_storage_to_fire_patch(i)             = value_patch
       this%m_leafn_xfer_to_fire_patch(i)                = value_patch
       this%m_livestemn_to_fire_patch(i)                 = value_patch
       this%m_livestemn_storage_to_fire_patch(i)         = value_patch
       this%m_livestemn_xfer_to_fire_patch(i)            = value_patch
       this%m_deadstemn_to_fire_patch(i)                 = value_patch
       this%m_deadstemn_storage_to_fire_patch(i)         = value_patch
       this%m_deadstemn_xfer_to_fire_patch(i)            = value_patch
       this%m_frootn_to_fire_patch(i)                    = value_patch
       this%m_frootn_storage_to_fire_patch(i)            = value_patch
       this%m_frootn_xfer_to_fire_patch(i)               = value_patch
       this%m_livecrootn_to_fire_patch(i)                = value_patch
       this%m_livecrootn_storage_to_fire_patch(i)        = value_patch
       this%m_livecrootn_xfer_to_fire_patch(i)           = value_patch
       this%m_deadcrootn_to_fire_patch(i)                = value_patch
       this%m_deadcrootn_storage_to_fire_patch(i)        = value_patch
       this%m_deadcrootn_xfer_to_fire_patch(i)           = value_patch
       this%m_retransn_to_fire_patch(i)                  = value_patch

       this%gru_conv_nflux_patch(i)                      = value_patch
       this%gru_wood_productn_gain_patch(i)              = value_patch

       this%m_leafn_to_litter_fire_patch(i)              = value_patch
       this%m_leafn_storage_to_litter_fire_patch(i)      = value_patch
       this%m_leafn_xfer_to_litter_fire_patch(i)         = value_patch
       this%m_livestemn_to_litter_fire_patch(i)          = value_patch
       this%m_livestemn_storage_to_litter_fire_patch(i)  = value_patch
       this%m_livestemn_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_livestemn_to_deadstemn_fire_patch(i)       = value_patch
       this%m_deadstemn_to_litter_fire_patch(i)          = value_patch
       this%m_deadstemn_storage_to_litter_fire_patch(i)  = value_patch
       this%m_deadstemn_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_frootn_to_litter_fire_patch(i)             = value_patch
       this%m_frootn_storage_to_litter_fire_patch(i)     = value_patch
       this%m_frootn_xfer_to_litter_fire_patch(i)        = value_patch
       this%m_livecrootn_to_litter_fire_patch(i)         = value_patch
       this%m_livecrootn_storage_to_litter_fire_patch(i) = value_patch
       this%m_livecrootn_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_livecrootn_to_deadcrootn_fire_patch(i)     = value_patch
       this%m_deadcrootn_to_litter_fire_patch(i)         = value_patch
       this%m_deadcrootn_storage_to_litter_fire_patch(i) = value_patch
       this%m_deadcrootn_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_retransn_to_litter_fire_patch(i)           = value_patch

       this%leafn_xfer_to_leafn_patch(i)                 = value_patch
       this%frootn_xfer_to_frootn_patch(i)               = value_patch
       this%livestemn_xfer_to_livestemn_patch(i)         = value_patch
       this%deadstemn_xfer_to_deadstemn_patch(i)         = value_patch
       this%livecrootn_xfer_to_livecrootn_patch(i)       = value_patch
       this%deadcrootn_xfer_to_deadcrootn_patch(i)       = value_patch
       this%leafn_to_litter_patch(i)                     = value_patch
       this%leafn_to_retransn_patch(i)                   = value_patch
       this%frootn_to_litter_patch(i)                    = value_patch
       this%retransn_to_npool_patch(i)                   = value_patch
       this%free_retransn_to_npool_patch(i)              = value_patch
       this%sminn_to_npool_patch(i)                      = value_patch
       this%npool_to_leafn_patch(i)                      = value_patch
       this%npool_to_leafn_storage_patch(i)              = value_patch
       this%npool_to_frootn_patch(i)                     = value_patch
       this%npool_to_frootn_storage_patch(i)             = value_patch
       this%npool_to_livestemn_patch(i)                  = value_patch
       this%npool_to_livestemn_storage_patch(i)          = value_patch
       this%npool_to_deadstemn_patch(i)                  = value_patch
       this%npool_to_deadstemn_storage_patch(i)          = value_patch
       this%npool_to_livecrootn_patch(i)                 = value_patch
       this%npool_to_livecrootn_storage_patch(i)         = value_patch
       this%npool_to_deadcrootn_patch(i)                 = value_patch
       this%npool_to_deadcrootn_storage_patch(i)         = value_patch
       this%leafn_storage_to_xfer_patch(i)               = value_patch
       this%frootn_storage_to_xfer_patch(i)              = value_patch
       this%livestemn_storage_to_xfer_patch(i)           = value_patch
       this%deadstemn_storage_to_xfer_patch(i)           = value_patch
       this%livecrootn_storage_to_xfer_patch(i)          = value_patch
       this%deadcrootn_storage_to_xfer_patch(i)          = value_patch
       this%livestemn_to_deadstemn_patch(i)              = value_patch
       this%livestemn_to_retransn_patch(i)               = value_patch
       this%livecrootn_to_deadcrootn_patch(i)            = value_patch
       this%livecrootn_to_retransn_patch(i)              = value_patch
       this%ndeploy_patch(i)                             = value_patch
       this%wood_harvestn_patch(i)                       = value_patch
       this%fire_nloss_patch(i)                          = value_patch

       this%crop_seedn_to_leaf_patch(i)                  = value_patch
       this%crop_harvestn_to_cropprodn_patch(i)                 = value_patch
    end do

    if ( use_crop )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%livestemn_to_litter_patch(i)              = value_patch
          this%leafn_to_biofueln_patch(i)                = value_patch
          this%livestemn_to_biofueln_patch(i)            = value_patch
          this%leafn_to_removedresiduen_patch(i)         = value_patch
          this%livestemn_to_removedresiduen_patch(i)     = value_patch
          this%soyfixn_patch(i)                          = value_patch
          this%frootn_to_retransn_patch(i)               = value_patch
       end do

       do k = 1, nrepr
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%reproductiven_xfer_to_reproductiven_patch(i,k) = value_patch
             this%npool_to_reproductiven_patch(i,k)                    = value_patch
             this%npool_to_reproductiven_storage_patch(i,k)            = value_patch
             this%reproductiven_storage_to_xfer_patch(i,k)             = value_patch
          end do
       end do

       do k = repr_grain_min, repr_grain_max
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%repr_grainn_to_food_patch(i,k) = value_patch
             this%repr_grainn_to_seed_patch(i,k) = value_patch
          end do
       end do

       do k = repr_structure_min, repr_structure_max
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%repr_structuren_to_cropprod_patch(i,k) = value_patch
             this%repr_structuren_to_litter_patch(i,k)   = value_patch
          end do
       end do
    end if

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)

          do k = i_litr_min, i_litr_max
             this%phenology_n_to_litr_n_col(i,j,k)       = value_column
             this%gap_mortality_n_to_litr_n_col(i,j,k)   = value_column
             this%harvest_n_to_litr_n_col(i,j,k)         = value_column
             this%m_n_to_litr_fire_col(i,j,k)            = value_column
             ! gross unrepresented landcover change
             this%gru_n_to_litr_n_col(i,j,k)             = value_column
          end do

          this%gap_mortality_n_to_cwdn_col(i,j)          = value_column
          this%fire_mortality_n_to_cwdn_col(i,j)         = value_column
          this%harvest_n_to_cwdn_col(i,j)                = value_column  

          ! gross unrepresented landcover change
          this%gru_n_to_cwdn_col(i,j)                    = value_column  
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)

       this%crop_harvestn_to_cropprodn_col(i)       = value_column
       this%fire_nloss_col(i)                = value_column

       ! Zero p2c column fluxes
       this%fire_nloss_col(i) = value_column
       this%wood_harvestn_col(i) = value_column
       this%gru_conv_nflux_col(i) = value_column 
       this%gru_wood_productn_gain_col(i) = value_column
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%m_decomp_npools_to_fire_col(i,k) = value_column
       end do
    end do
! Matrix
    if(use_matrixcn)then
       do j = 1, nvegnpool
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_nalloc_patch(i,j)       = value_patch
             this%matrix_nphturnover_patch (i,j) = value_patch
             this%matrix_ngmturnover_patch (i,j) = value_patch
             this%matrix_nfiturnover_patch (i,j) = value_patch
          end do
       end do

       do j = 1, nnphtrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_nphtransfer_patch (i,j) = value_patch
          end do
       end do

       do j = 1, nngmtrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_ngmtransfer_patch (i,j) = value_patch
          end do
       end do

       do j = 1, nnfitrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_nfitransfer_patch (i,j) = value_patch
          end do
       end do

    end if
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%m_decomp_npools_to_fire_vr_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenflux_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer :: c, g, j, i       ! indices
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg
       this%dwt_seedn_to_leaf_grc(g)     = 0._r8
       this%dwt_seedn_to_deadstem_grc(g) = 0._r8
       this%dwt_conv_nflux_grc(g)        = 0._r8
    end do

    do j = 1, nlevdecomp_full
       do c = bounds%begc,bounds%endc
          do i = i_litr_min, i_litr_max
             this%dwt_frootn_to_litr_n_col(c,j,i) = 0._r8
          end do
          this%dwt_livecrootn_to_cwdn_col(c,j)   = 0._r8
          this%dwt_deadcrootn_to_cwdn_col(c,j)   = 0._r8
       end do
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine ZeroGru( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenflux_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: g  ! indices
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg
       this%gru_conv_nflux_grc(g)        = 0._r8
       this%gru_wood_productn_gain_grc(g) = 0._r8
    end do

  end subroutine ZeroGru

 !-----------------------------------------------------------------------
  subroutine Summary_nitrogenflux(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !USES:
    use clm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c, c2g 
    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenflux_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in) :: num_soilp       ! number of soil patches in filter
    integer           , intent(in) :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l   ! indices
    integer  :: fp,fc       ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! total N deployment (from sminn and retranslocated N pool) (NDEPLOY)
       this%ndeploy_patch(p) = &
            this%sminn_to_npool_patch(p) + &
            this%retransn_to_npool_patch(p) + &
            this%free_retransn_to_npool_patch(p)  

       ! total patch-level fire N losses
       this%fire_nloss_patch(p) = &
            this%m_leafn_to_fire_patch(p)               + &
            this%m_leafn_storage_to_fire_patch(p)       + &
            this%m_leafn_xfer_to_fire_patch(p)          + &
            this%m_frootn_to_fire_patch(p)              + &
            this%m_frootn_storage_to_fire_patch(p)      + &
            this%m_frootn_xfer_to_fire_patch(p)         + &
            this%m_livestemn_to_fire_patch(p)           + &
            this%m_livestemn_storage_to_fire_patch(p)   + &
            this%m_livestemn_xfer_to_fire_patch(p)      + &
            this%m_deadstemn_to_fire_patch(p)           + &
            this%m_deadstemn_storage_to_fire_patch(p)   + &
            this%m_deadstemn_xfer_to_fire_patch(p)      + &
            this%m_livecrootn_to_fire_patch(p)          + &
            this%m_livecrootn_storage_to_fire_patch(p)  + &
            this%m_livecrootn_xfer_to_fire_patch(p)     + &
            this%m_deadcrootn_to_fire_patch(p)          + &
            this%m_deadcrootn_storage_to_fire_patch(p)  + &
            this%m_deadcrootn_xfer_to_fire_patch(p)     + &
            this%m_retransn_to_fire_patch(p)

       ! (Gross Unrepresented Landcover Change Conversion Flux) - Direct Veg N Loss to Atmosphere
       this%gru_conv_nflux_patch(p) = &
            this%gru_livestemn_to_atm_patch(p)          + &
            this%gru_deadstemn_to_atm_patch(p)          + &
            this%gru_leafn_storage_to_atm_patch(p)      + &
            this%gru_frootn_storage_to_atm_patch(p)     + &
            this%gru_livestemn_storage_to_atm_patch(p)  + &
            this%gru_deadstemn_storage_to_atm_patch(p)  + &
            this%gru_livecrootn_storage_to_atm_patch(p) + &
            this%gru_deadcrootn_storage_to_atm_patch(p) + &
            this%gru_leafn_xfer_to_atm_patch(p)         + &
            this%gru_frootn_xfer_to_atm_patch(p)        + &
            this%gru_livestemn_xfer_to_atm_patch(p)     + &
            this%gru_deadstemn_xfer_to_atm_patch(p)     + &
            this%gru_livecrootn_xfer_to_atm_patch(p)    + &
            this%gru_deadcrootn_xfer_to_atm_patch(p)

    end do

    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_nloss_patch(bounds%begp:bounds%endp), &
         this%fire_nloss_p2c_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%gru_conv_nflux_patch(bounds%begp:bounds%endp), &
         this%gru_conv_nflux_col(bounds%begc:bounds%endc))

    call c2g( bounds = bounds, &
         carr = this%gru_conv_nflux_col(bounds%begc:bounds%endc), &
         garr = this%gru_conv_nflux_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    call c2g( bounds = bounds, &
         carr = this%gru_wood_productn_gain_col(bounds%begc:bounds%endc), &
         garr = this%gru_wood_productn_gain_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    ! vertically integrate column-level fire N losses
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%m_decomp_npools_to_fire_col(c,k) = &
                  this%m_decomp_npools_to_fire_col(c,k) + &
                  this%m_decomp_npools_to_fire_vr_col(c,j,k) * dzsoi_decomp(j)
          end do
       end do
    end do

    ! total column-level fire N losses
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%fire_nloss_col(c) = this%fire_nloss_p2c_col(c)
    end do
    do k = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%fire_nloss_col(c) = &
               this%fire_nloss_col(c) + &
               this%m_decomp_npools_to_fire_col(c,k)
       end do
    end do

  end subroutine Summary_nitrogenflux

end module CNVegNitrogenFluxType

