module CNVegMatrixMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Matrix solution for vegetation C and N cycles
  ! The matrix equation 
  ! Xn+1 = Xn + I*dt + (A*ksi*k)*Xn*dt
  
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use clm_time_manager               , only : get_step_size,is_end_curr_year,is_first_step_of_this_run_segment,&
                                              get_days_per_year,is_beg_curr_year
  use decompMod                      , only : bounds_type 
  use clm_varpar                     , only : nlevdecomp, nvegcpool, nvegnpool
  use clm_varpar                     , only : ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                              ilivestem,ilivestem_st,ilivestem_xf,&
                                              ideadstem,ideadstem_st,ideadstem_xf,&
                                              ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                              ideadcroot,ideadcroot_st,ideadcroot_xf,&
                                              igrain,igrain_st,igrain_xf,iretransn,ioutc,ioutn

  use PatchType                      , only : patch
  use GridcellType                   , only : grc

    use clm_varcon                   , only: secspday
  use pftconMod                      , only : pftcon,npcropmin
  use CNVegCarbonStateType           , only : cnveg_carbonstate_type
  use CNVegNitrogenStateType         , only : cnveg_nitrogenstate_type
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type     !include: callocation,ctransfer, cturnover
  use CNVegNitrogenFluxType          , only : cnveg_nitrogenflux_type
  use CNVegStateType                 , only : cnveg_state_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use clm_varctl                     , only : isspinup, is_outmatrix
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNVegMatrix
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
   subroutine CNVegMatrix(bounds,num_soilp,filter_soilp,cnveg_carbonstate_inst,cnveg_nitrogenstate_inst,&
                        cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,cnveg_state_inst,soilbiogeochem_nitrogenflux_inst)
    ! !DESCRIPTION:
    ! !ARGUMENTS:
     type(bounds_type)                      , intent(in)    :: bounds
     integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
     integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
     type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst  
     type(cnveg_nitrogenstate_type)           , intent(inout) :: cnveg_nitrogenstate_inst
     type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
     type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
     type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
 
     type(cnveg_state_type)    , intent(in) :: cnveg_state_inst
!	
    ! !LOCAL VARIABLES:
     integer :: fc,fp,j,i    ! indices
     integer :: p,c         !  
     real(r8),allocatable,dimension(:,:)   :: vegmatrixc_old       (:,:)
     real(r8),allocatable,dimension(:)     :: vegmatrixc_new       (:)
     real(r8),allocatable,dimension(:,:)   :: vegmatrixn_old       (:,:)
     real(r8),allocatable,dimension(:)     :: vegmatrixn_new       (:)
     real(r8),allocatable,dimension(:)     :: vegmatrixc_input     (:)
     real(r8),allocatable,dimension(:,:)   :: vegmatrixc_transfer  (:,:)
     real(r8),allocatable,dimension(:)     :: vegmatrixn_input     (:)
     real(r8),allocatable,dimension(:,:)   :: vegmatrixn_transfer  (:,:)
     real(r8),allocatable,dimension(:)     :: matrix_calloc_acc    (:)
     real(r8),allocatable,dimension(:)     :: matrix_nalloc_acc    (:)
     real(r8),allocatable,dimension(:,:)   :: matrix_ctransfer_acc (:,:)
     real(r8),allocatable,dimension(:,:)   :: matrix_ntransfer_acc (:,:)
!     real(r8),allocatable,dimension(:,:) :: vegmatrixc_rt(:,:),vegmatrixn_rt(:,:)
!     real(r8),allocatable,dimension(:,:,:) :: matrix_nphturnover(:,:,:),matrix_ngmturnover(:,:,:)!,matrix_nphturnover(:,:,:)
     real(r8),allocatable,dimension(:,:,:) :: matrix_nphtransfer_curr(:,:,:)

! for spinupacc
     real(r8),dimension(1:nvegcpool) :: vegmatrixc_rt,rowonec
     real(r8),dimension(1:nvegnpool) :: vegmatrixn_rt,rowonen,tmp
     real(r8),allocatable,dimension(:,:) :: AKinvc(:,:),AKinvn(:,:)
     logical  ::  end_of_year
     real(r8),save :: counter=0

     real(r8):: dt        ! time step (seconds)
     real(r8):: secspyear        ! time step (seconds)
 
!	
!
    associate(                          &                                                                        
         ivt                          => patch%itype                                               , & ! Input:  [integer  (:) ]  patch vegetation type
         woody                        => pftcon%woody                                              , & ! Input:  binary flag for woody lifeform (1=woody, 0=not woody)
         leafcn                       => pftcon%leafcn                                             , & ! Input:  leaf C:N (gC/gN)                        
         frootcn                      => pftcon%frootcn                                            , & ! Input:  fine root C:N (gC/gN)                   
         livewdcn                     => pftcon%livewdcn                                           , & ! Input:  live wood (phloem and ray parenchyma) C:N (gC/gN)
         deadwdcn                     => pftcon%deadwdcn                                           , & ! Input:  dead wood (xylem and heartwood) C:N (gC/gN)
         c_allometry                  => cnveg_state_inst%c_allometry_patch                        , & ! Output: [real(r8) (:)   ]  C allocation index (DIM) 
         n_allometry                  => cnveg_state_inst%n_allometry_patch                        , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)  
         retransn_to_npool            => cnveg_nitrogenflux_inst%retransn_to_npool_patch           , & !
         plant_nalloc                 => cnveg_nitrogenflux_inst%plant_nalloc_patch                , & ! Output: 
         retransn              => cnveg_nitrogenstate_inst%retransn_patch           , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N
  !
 
         leafc                 => cnveg_carbonstate_inst%leafc_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for matrix calculation                                    
         leafc_storage         => cnveg_carbonstate_inst%leafc_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage C for matrix calculation                                   
         leafc_xfer            => cnveg_carbonstate_inst%leafc_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer C for matrix calcuation                                   
         frootc                => cnveg_carbonstate_inst%frootc_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C for matrix calculation
         frootc_storage        => cnveg_carbonstate_inst%frootc_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C for matrix calculation
         frootc_xfer           => cnveg_carbonstate_inst%frootc_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C for matrix calculation
         livestemc             => cnveg_carbonstate_inst%livestemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem C for matrix calculationleaf C
         livestemc_storage     => cnveg_carbonstate_inst%livestemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage C for matrix calculation
         livestemc_xfer        => cnveg_carbonstate_inst%livestemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer C for matrix calculation
         deadstemc             => cnveg_carbonstate_inst%deadstemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem C for matrix calculationleaf C
         deadstemc_storage     => cnveg_carbonstate_inst%deadstemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage C for matrix calculation
         deadstemc_xfer        => cnveg_carbonstate_inst%deadstemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer C for matrix calculation                                    
         livecrootc            => cnveg_carbonstate_inst%livecrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root C for matrix calculationleaf C 
         livecrootc_storage    => cnveg_carbonstate_inst%livecrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage C for matrix calculation
         livecrootc_xfer       => cnveg_carbonstate_inst%livecrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer C for matrix calculation
         deadcrootc            => cnveg_carbonstate_inst%deadcrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root C for matrix calculationleaf C
         deadcrootc_storage    => cnveg_carbonstate_inst%deadcrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage C for matrix calculation
         deadcrootc_xfer       => cnveg_carbonstate_inst%deadcrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer C for matrix calculation
         grainc                => cnveg_carbonstate_inst%grainc_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root C for matrix calculationleaf C
         grainc_storage        => cnveg_carbonstate_inst%grainc_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage C for matrix calculation
         grainc_xfer           => cnveg_carbonstate_inst%grainc_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer C for matrix calculation
         matrix_cap_leafc                 => cnveg_carbonstate_inst%matrix_cap_leafc_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for matrix calculation                                    
         matrix_cap_leafc_storage         => cnveg_carbonstate_inst%matrix_cap_leafc_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage C for matrix calculation                                   
         matrix_cap_leafc_xfer            => cnveg_carbonstate_inst%matrix_cap_leafc_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer C for matrix calcuation                                   
         matrix_cap_frootc                => cnveg_carbonstate_inst%matrix_cap_frootc_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C for matrix calculation
         matrix_cap_frootc_storage        => cnveg_carbonstate_inst%matrix_cap_frootc_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C for matrix calculation
         matrix_cap_frootc_xfer           => cnveg_carbonstate_inst%matrix_cap_frootc_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C for matrix calculation
         matrix_cap_livestemc             => cnveg_carbonstate_inst%matrix_cap_livestemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem C for matrix calculationleaf C
         matrix_cap_livestemc_storage     => cnveg_carbonstate_inst%matrix_cap_livestemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage C for matrix calculation
         matrix_cap_livestemc_xfer        => cnveg_carbonstate_inst%matrix_cap_livestemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer C for matrix calculation
         matrix_cap_deadstemc             => cnveg_carbonstate_inst%matrix_cap_deadstemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem C for matrix calculationleaf C
         matrix_cap_deadstemc_storage     => cnveg_carbonstate_inst%matrix_cap_deadstemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage C for matrix calculation
         matrix_cap_deadstemc_xfer        => cnveg_carbonstate_inst%matrix_cap_deadstemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer C for matrix calculation                                    
         matrix_cap_livecrootc            => cnveg_carbonstate_inst%matrix_cap_livecrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root C for matrix calculationleaf C 
         matrix_cap_livecrootc_storage    => cnveg_carbonstate_inst%matrix_cap_livecrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage C for matrix calculation
         matrix_cap_livecrootc_xfer       => cnveg_carbonstate_inst%matrix_cap_livecrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer C for matrix calculation
         matrix_cap_deadcrootc            => cnveg_carbonstate_inst%matrix_cap_deadcrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root C for matrix calculationleaf C
         matrix_cap_deadcrootc_storage    => cnveg_carbonstate_inst%matrix_cap_deadcrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage C for matrix calculation
         matrix_cap_deadcrootc_xfer       => cnveg_carbonstate_inst%matrix_cap_deadcrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer C for matrix calculation
         matrix_cap_grainc                => cnveg_carbonstate_inst%matrix_cap_grainc_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for matrix calculation                                    
         matrix_cap_grainc_storage        => cnveg_carbonstate_inst%matrix_cap_grainc_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage C for matrix calculation                                   
         matrix_cap_grainc_xfer           => cnveg_carbonstate_inst%matrix_cap_grainc_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer C for matrix calcuation                                   
         matrix_pot_leafc                 => cnveg_carbonstate_inst%matrix_pot_leafc_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for matrix calculation                                    
         matrix_pot_leafc_storage         => cnveg_carbonstate_inst%matrix_pot_leafc_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage C for matrix calculation                                   
         matrix_pot_leafc_xfer            => cnveg_carbonstate_inst%matrix_pot_leafc_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer C for matrix calcuation                                   
         matrix_pot_frootc                => cnveg_carbonstate_inst%matrix_pot_frootc_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C for matrix calculation
         matrix_pot_frootc_storage        => cnveg_carbonstate_inst%matrix_pot_frootc_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C for matrix calculation
         matrix_pot_frootc_xfer           => cnveg_carbonstate_inst%matrix_pot_frootc_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C for matrix calculation
         matrix_pot_livestemc             => cnveg_carbonstate_inst%matrix_pot_livestemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem C for matrix calculationleaf C
         matrix_pot_livestemc_storage     => cnveg_carbonstate_inst%matrix_pot_livestemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage C for matrix calculation
         matrix_pot_livestemc_xfer        => cnveg_carbonstate_inst%matrix_pot_livestemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer C for matrix calculation
         matrix_pot_deadstemc             => cnveg_carbonstate_inst%matrix_pot_deadstemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem C for matrix calculationleaf C
         matrix_pot_deadstemc_storage     => cnveg_carbonstate_inst%matrix_pot_deadstemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage C for matrix calculation
         matrix_pot_deadstemc_xfer        => cnveg_carbonstate_inst%matrix_pot_deadstemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer C for matrix calculation                                    
         matrix_pot_livecrootc            => cnveg_carbonstate_inst%matrix_pot_livecrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root C for matrix calculationleaf C 
         matrix_pot_livecrootc_storage    => cnveg_carbonstate_inst%matrix_pot_livecrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage C for matrix calculation
         matrix_pot_livecrootc_xfer       => cnveg_carbonstate_inst%matrix_pot_livecrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer C for matrix calculation
         matrix_pot_deadcrootc            => cnveg_carbonstate_inst%matrix_pot_deadcrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root C for matrix calculationleaf C
         matrix_pot_deadcrootc_storage    => cnveg_carbonstate_inst%matrix_pot_deadcrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage C for matrix calculation
         matrix_pot_deadcrootc_xfer       => cnveg_carbonstate_inst%matrix_pot_deadcrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer C for matrix calculation
         matrix_pot_grainc                => cnveg_carbonstate_inst%matrix_pot_grainc_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C for matrix calculation
         matrix_pot_grainc_storage        => cnveg_carbonstate_inst%matrix_pot_grainc_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C for matrix calculation
         matrix_pot_grainc_xfer           => cnveg_carbonstate_inst%matrix_pot_grainc_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C for matrix calculation
!
         leafn                 => cnveg_nitrogenstate_inst%leafn_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf N for matrix calculation                                    
         leafn_storage         => cnveg_nitrogenstate_inst%leafn_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage N for matrix calculation                                   
         leafn_xfer            => cnveg_nitrogenstate_inst%leafn_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer N for  matrix calcuation                                   
         frootn                => cnveg_nitrogenstate_inst%frootn_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         frootn_storage        => cnveg_nitrogenstate_inst%frootn_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
         frootn_xfer           => cnveg_nitrogenstate_inst%frootn_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
         livestemn             => cnveg_nitrogenstate_inst%livestemn_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem N for matrix calculationleaf N
         livestemn_storage     => cnveg_nitrogenstate_inst%livestemn_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage N for matrix calculation
         livestemn_xfer        => cnveg_nitrogenstate_inst%livestemn_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer N for matrix calculation
         deadstemn             => cnveg_nitrogenstate_inst%deadstemn_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem N for matrix calculationleaf N
         deadstemn_storage     => cnveg_nitrogenstate_inst%deadstemn_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage N for matrix calculation
         deadstemn_xfer        => cnveg_nitrogenstate_inst%deadstemn_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer N for matrix calculation                                    
         livecrootn            => cnveg_nitrogenstate_inst%livecrootn_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root N for matrix calculationleaf N 
         livecrootn_storage    => cnveg_nitrogenstate_inst%livecrootn_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage N for matrix calculation
         livecrootn_xfer       => cnveg_nitrogenstate_inst%livecrootn_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer N for matrix calculation
         deadcrootn            => cnveg_nitrogenstate_inst%deadcrootn_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root N for matrix calculationleaf N
         deadcrootn_storage    => cnveg_nitrogenstate_inst%deadcrootn_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage N for matrix calculation
         deadcrootn_xfer       => cnveg_nitrogenstate_inst%deadcrootn_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer N for matrix calculation
         grainn                => cnveg_nitrogenstate_inst%grainn_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         grainn_storage        => cnveg_nitrogenstate_inst%grainn_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
         grainn_xfer           => cnveg_nitrogenstate_inst%grainn_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
         matrix_cap_leafn                 => cnveg_nitrogenstate_inst%matrix_cap_leafn_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf N for matrix calculation                                    
         matrix_cap_leafn_storage         => cnveg_nitrogenstate_inst%matrix_cap_leafn_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage N for matrix calculation                                   
         matrix_cap_leafn_xfer            => cnveg_nitrogenstate_inst%matrix_cap_leafn_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer N for  matrix calcuation                                   
         matrix_cap_frootn                => cnveg_nitrogenstate_inst%matrix_cap_frootn_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         matrix_cap_frootn_storage        => cnveg_nitrogenstate_inst%matrix_cap_frootn_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
         matrix_cap_frootn_xfer           => cnveg_nitrogenstate_inst%matrix_cap_frootn_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
         matrix_cap_livestemn             => cnveg_nitrogenstate_inst%matrix_cap_livestemn_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem N for matrix calculationleaf N
         matrix_cap_livestemn_storage     => cnveg_nitrogenstate_inst%matrix_cap_livestemn_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage N for matrix calculation
         matrix_cap_livestemn_xfer        => cnveg_nitrogenstate_inst%matrix_cap_livestemn_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer N for matrix calculation
         matrix_cap_deadstemn             => cnveg_nitrogenstate_inst%matrix_cap_deadstemn_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem N for matrix calculationleaf N
         matrix_cap_deadstemn_storage     => cnveg_nitrogenstate_inst%matrix_cap_deadstemn_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage N for matrix calculation
         matrix_cap_deadstemn_xfer        => cnveg_nitrogenstate_inst%matrix_cap_deadstemn_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer N for matrix calculation                                    
         matrix_cap_livecrootn            => cnveg_nitrogenstate_inst%matrix_cap_livecrootn_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root N for matrix calculationleaf N 
         matrix_cap_livecrootn_storage    => cnveg_nitrogenstate_inst%matrix_cap_livecrootn_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage N for matrix calculation
         matrix_cap_livecrootn_xfer       => cnveg_nitrogenstate_inst%matrix_cap_livecrootn_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer N for matrix calculation
         matrix_cap_deadcrootn            => cnveg_nitrogenstate_inst%matrix_cap_deadcrootn_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root N for matrix calculationleaf N
         matrix_cap_deadcrootn_storage    => cnveg_nitrogenstate_inst%matrix_cap_deadcrootn_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage N for matrix calculation
         matrix_cap_deadcrootn_xfer       => cnveg_nitrogenstate_inst%matrix_cap_deadcrootn_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer N for matrix calculation
         matrix_cap_grainn                => cnveg_nitrogenstate_inst%matrix_cap_grainn_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         matrix_cap_grainn_storage        => cnveg_nitrogenstate_inst%matrix_cap_grainn_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
         matrix_cap_grainn_xfer           => cnveg_nitrogenstate_inst%matrix_cap_grainn_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
         matrix_pot_leafn                 => cnveg_nitrogenstate_inst%matrix_pot_leafn_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf N for matrix calculation                                    
         matrix_pot_leafn_storage         => cnveg_nitrogenstate_inst%matrix_pot_leafn_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage N for matrix calculation                                   
         matrix_pot_leafn_xfer            => cnveg_nitrogenstate_inst%matrix_pot_leafn_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer N for  matrix calcuation                                   
         matrix_pot_frootn                => cnveg_nitrogenstate_inst%matrix_pot_frootn_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         matrix_pot_frootn_storage        => cnveg_nitrogenstate_inst%matrix_pot_frootn_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
         matrix_pot_frootn_xfer           => cnveg_nitrogenstate_inst%matrix_pot_frootn_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
         matrix_pot_livestemn             => cnveg_nitrogenstate_inst%matrix_pot_livestemn_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem N for matrix calculationleaf N
         matrix_pot_livestemn_storage     => cnveg_nitrogenstate_inst%matrix_pot_livestemn_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage N for matrix calculation
         matrix_pot_livestemn_xfer        => cnveg_nitrogenstate_inst%matrix_pot_livestemn_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer N for matrix calculation
         matrix_pot_deadstemn             => cnveg_nitrogenstate_inst%matrix_pot_deadstemn_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem N for matrix calculationleaf N
         matrix_pot_deadstemn_storage     => cnveg_nitrogenstate_inst%matrix_pot_deadstemn_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage N for matrix calculation
         matrix_pot_deadstemn_xfer        => cnveg_nitrogenstate_inst%matrix_pot_deadstemn_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer N for matrix calculation                                    
         matrix_pot_livecrootn            => cnveg_nitrogenstate_inst%matrix_pot_livecrootn_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root N for matrix calculationleaf N 
         matrix_pot_livecrootn_storage    => cnveg_nitrogenstate_inst%matrix_pot_livecrootn_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage N for matrix calculation
         matrix_pot_livecrootn_xfer       => cnveg_nitrogenstate_inst%matrix_pot_livecrootn_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer N for matrix calculation
         matrix_pot_deadcrootn            => cnveg_nitrogenstate_inst%matrix_pot_deadcrootn_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root N for matrix calculationleaf N
         matrix_pot_deadcrootn_storage    => cnveg_nitrogenstate_inst%matrix_pot_deadcrootn_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage N for matrix calculation
         matrix_pot_deadcrootn_xfer       => cnveg_nitrogenstate_inst%matrix_pot_deadcrootn_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer N for matrix calculation
         matrix_pot_grainn                => cnveg_nitrogenstate_inst%matrix_pot_grainn_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         matrix_pot_grainn_storage        => cnveg_nitrogenstate_inst%matrix_pot_grainn_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
         matrix_pot_grainn_xfer           => cnveg_nitrogenstate_inst%matrix_pot_grainn_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
         leafc0                 => cnveg_carbonstate_inst%leafc0_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for matrix calculation                                    
         leafc0_storage         => cnveg_carbonstate_inst%leafc0_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage C for matrix calculation                                   
         leafc0_xfer            => cnveg_carbonstate_inst%leafc0_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer C for matrix calcuation                                   
         frootc0                => cnveg_carbonstate_inst%frootc0_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C for matrix calculation
         frootc0_storage        => cnveg_carbonstate_inst%frootc0_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C for matrix calculation
         frootc0_xfer           => cnveg_carbonstate_inst%frootc0_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C for matrix calculation
         livestemc0             => cnveg_carbonstate_inst%livestemc0_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem C for matrix calculationleaf C
         livestemc0_storage     => cnveg_carbonstate_inst%livestemc0_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage C for matrix calculation
         livestemc0_xfer        => cnveg_carbonstate_inst%livestemc0_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer C for matrix calculation
         deadstemc0             => cnveg_carbonstate_inst%deadstemc0_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem C for matrix calculationleaf C
         deadstemc0_storage     => cnveg_carbonstate_inst%deadstemc0_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage C for matrix calculation
         deadstemc0_xfer        => cnveg_carbonstate_inst%deadstemc0_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer C for matrix calculation                                    
         livecrootc0            => cnveg_carbonstate_inst%livecrootc0_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root C for matrix calculationleaf C 
         livecrootc0_storage    => cnveg_carbonstate_inst%livecrootc0_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage C for matrix calculation
         livecrootc0_xfer       => cnveg_carbonstate_inst%livecrootc0_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer C for matrix calculation
         deadcrootc0            => cnveg_carbonstate_inst%deadcrootc0_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root C for matrix calculationleaf C
         deadcrootc0_storage    => cnveg_carbonstate_inst%deadcrootc0_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage C for matrix calculation
         deadcrootc0_xfer       => cnveg_carbonstate_inst%deadcrootc0_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer C for matrix calculation
         grainc0                => cnveg_carbonstate_inst%grainc0_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C for matrix calculation
         grainc0_storage        => cnveg_carbonstate_inst%grainc0_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C for matrix calculation
         grainc0_xfer           => cnveg_carbonstate_inst%grainc0_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C for matrix calculation
!
         leafn0                 => cnveg_nitrogenstate_inst%leafn0_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf N for matrix calculation                                    
         leafn0_storage         => cnveg_nitrogenstate_inst%leafn0_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage N for matrix calculation                                   
         leafn0_xfer            => cnveg_nitrogenstate_inst%leafn0_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer N for  matrix calcuation                                   
         frootn0                => cnveg_nitrogenstate_inst%frootn0_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         frootn0_storage        => cnveg_nitrogenstate_inst%frootn0_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
         frootn0_xfer           => cnveg_nitrogenstate_inst%frootn0_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
         livestemn0             => cnveg_nitrogenstate_inst%livestemn0_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem N for matrix calculationleaf N
         livestemn0_storage     => cnveg_nitrogenstate_inst%livestemn0_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage N for matrix calculation
         livestemn0_xfer        => cnveg_nitrogenstate_inst%livestemn0_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer N for matrix calculation
         deadstemn0             => cnveg_nitrogenstate_inst%deadstemn0_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem N for matrix calculationleaf N
         deadstemn0_storage     => cnveg_nitrogenstate_inst%deadstemn0_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage N for matrix calculation
         deadstemn0_xfer        => cnveg_nitrogenstate_inst%deadstemn0_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer N for matrix calculation                                    
         livecrootn0            => cnveg_nitrogenstate_inst%livecrootn0_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root N for matrix calculationleaf N 
         livecrootn0_storage    => cnveg_nitrogenstate_inst%livecrootn0_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage N for matrix calculation
         livecrootn0_xfer       => cnveg_nitrogenstate_inst%livecrootn0_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer N for matrix calculation
         deadcrootn0            => cnveg_nitrogenstate_inst%deadcrootn0_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root N for matrix calculationleaf N
         deadcrootn0_storage    => cnveg_nitrogenstate_inst%deadcrootn0_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage N for matrix calculation
         deadcrootn0_xfer       => cnveg_nitrogenstate_inst%deadcrootn0_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer N for matrix calculation
         grainn0                => cnveg_nitrogenstate_inst%grainn0_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         grainn0_storage        => cnveg_nitrogenstate_inst%grainn0_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
         grainn0_xfer           => cnveg_nitrogenstate_inst%grainn0_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
         retransn0              => cnveg_nitrogenstate_inst%retransn0_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
         matrix_calloc_leaf_acc        => cnveg_carbonstate_inst%matrix_calloc_leaf_acc_patch             , & !
         matrix_calloc_leafst_acc      => cnveg_carbonstate_inst%matrix_calloc_leafst_acc_patch             , & !
         matrix_calloc_froot_acc       => cnveg_carbonstate_inst%matrix_calloc_froot_acc_patch             , & !
         matrix_calloc_frootst_acc     => cnveg_carbonstate_inst%matrix_calloc_frootst_acc_patch             , & !
         matrix_calloc_livestem_acc    => cnveg_carbonstate_inst%matrix_calloc_livestem_acc_patch             , & !
         matrix_calloc_livestemst_acc  => cnveg_carbonstate_inst%matrix_calloc_livestemst_acc_patch             , & !
         matrix_calloc_deadstem_acc    => cnveg_carbonstate_inst%matrix_calloc_deadstem_acc_patch             , & !
         matrix_calloc_deadstemst_acc  => cnveg_carbonstate_inst%matrix_calloc_deadstemst_acc_patch             , & !
         matrix_calloc_livecroot_acc   => cnveg_carbonstate_inst%matrix_calloc_livecroot_acc_patch             , & !
         matrix_calloc_livecrootst_acc => cnveg_carbonstate_inst%matrix_calloc_livecrootst_acc_patch             , & !
         matrix_calloc_deadcroot_acc   => cnveg_carbonstate_inst%matrix_calloc_deadcroot_acc_patch             , & !
         matrix_calloc_deadcrootst_acc => cnveg_carbonstate_inst%matrix_calloc_deadcrootst_acc_patch             , & !
         matrix_calloc_grain_acc       => cnveg_carbonstate_inst%matrix_calloc_grain_acc_patch             , & !
         matrix_calloc_grainst_acc     => cnveg_carbonstate_inst%matrix_calloc_grainst_acc_patch             , & !

         matrix_ctransfer_leafst_to_leafxf_acc            => cnveg_carbonstate_inst%matrix_ctransfer_leafst_to_leafxf_acc_patch          , & !
         matrix_ctransfer_leafxf_to_leaf_acc              => cnveg_carbonstate_inst%matrix_ctransfer_leafxf_to_leaf_acc_patch          , & !
         matrix_ctransfer_frootst_to_frootxf_acc          => cnveg_carbonstate_inst%matrix_ctransfer_frootst_to_frootxf_acc_patch          , & !
         matrix_ctransfer_frootxf_to_froot_acc            => cnveg_carbonstate_inst%matrix_ctransfer_frootxf_to_froot_acc_patch          , & !
         matrix_ctransfer_livestemst_to_livestemxf_acc    => cnveg_carbonstate_inst%matrix_ctransfer_livestemst_to_livestemxf_acc_patch          , & !
         matrix_ctransfer_livestemxf_to_livestem_acc      => cnveg_carbonstate_inst%matrix_ctransfer_livestemxf_to_livestem_acc_patch          , & !
         matrix_ctransfer_deadstemst_to_deadstemxf_acc    => cnveg_carbonstate_inst%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch          , & !
         matrix_ctransfer_deadstemxf_to_deadstem_acc      => cnveg_carbonstate_inst%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch          , & !
         matrix_ctransfer_livecrootst_to_livecrootxf_acc  => cnveg_carbonstate_inst%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch          , & !
         matrix_ctransfer_livecrootxf_to_livecroot_acc    => cnveg_carbonstate_inst%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch          , & !
         matrix_ctransfer_deadcrootst_to_deadcrootxf_acc  => cnveg_carbonstate_inst%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch          , & !
         matrix_ctransfer_deadcrootxf_to_deadcroot_acc    => cnveg_carbonstate_inst%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch          , & !
         matrix_ctransfer_grainst_to_grainxf_acc          => cnveg_carbonstate_inst%matrix_ctransfer_grainst_to_grainxf_acc_patch          , & !
         matrix_ctransfer_grainxf_to_grain_acc            => cnveg_carbonstate_inst%matrix_ctransfer_grainxf_to_grain_acc_patch          , & !
         matrix_ctransfer_livestem_to_deadstem_acc        => cnveg_carbonstate_inst%matrix_ctransfer_livestem_to_deadstem_acc_patch          , & !
         matrix_ctransfer_livecroot_to_deadcroot_acc      => cnveg_carbonstate_inst%matrix_ctransfer_livecroot_to_deadcroot_acc_patch          , & !

         matrix_cturnover_leaf_acc                        => cnveg_carbonstate_inst%matrix_cturnover_leaf_acc_patch               , &
         matrix_cturnover_leafst_acc                      => cnveg_carbonstate_inst%matrix_cturnover_leafst_acc_patch               , &
         matrix_cturnover_leafxf_acc                      => cnveg_carbonstate_inst%matrix_cturnover_leafxf_acc_patch               , &
         matrix_cturnover_froot_acc                       => cnveg_carbonstate_inst%matrix_cturnover_froot_acc_patch               , &
         matrix_cturnover_frootst_acc                     => cnveg_carbonstate_inst%matrix_cturnover_frootst_acc_patch               , &
         matrix_cturnover_frootxf_acc                     => cnveg_carbonstate_inst%matrix_cturnover_frootxf_acc_patch               , &
         matrix_cturnover_livestem_acc                    => cnveg_carbonstate_inst%matrix_cturnover_livestem_acc_patch               , &
         matrix_cturnover_livestemst_acc                  => cnveg_carbonstate_inst%matrix_cturnover_livestemst_acc_patch               , &
         matrix_cturnover_livestemxf_acc                  => cnveg_carbonstate_inst%matrix_cturnover_livestemxf_acc_patch               , &
         matrix_cturnover_deadstem_acc                    => cnveg_carbonstate_inst%matrix_cturnover_deadstem_acc_patch               , &
         matrix_cturnover_deadstemst_acc                  => cnveg_carbonstate_inst%matrix_cturnover_deadstemst_acc_patch               , &
         matrix_cturnover_deadstemxf_acc                  => cnveg_carbonstate_inst%matrix_cturnover_deadstemxf_acc_patch               , &
         matrix_cturnover_livecroot_acc                   => cnveg_carbonstate_inst%matrix_cturnover_livecroot_acc_patch               , &
         matrix_cturnover_livecrootst_acc                 => cnveg_carbonstate_inst%matrix_cturnover_livecrootst_acc_patch               , &
         matrix_cturnover_livecrootxf_acc                 => cnveg_carbonstate_inst%matrix_cturnover_livecrootxf_acc_patch               , &
         matrix_cturnover_deadcroot_acc                   => cnveg_carbonstate_inst%matrix_cturnover_deadcroot_acc_patch               , &
         matrix_cturnover_deadcrootst_acc                 => cnveg_carbonstate_inst%matrix_cturnover_deadcrootst_acc_patch               , &
         matrix_cturnover_deadcrootxf_acc                 => cnveg_carbonstate_inst%matrix_cturnover_deadcrootxf_acc_patch               , &
         matrix_cturnover_grain_acc                       => cnveg_carbonstate_inst%matrix_cturnover_grain_acc_patch               , &
         matrix_cturnover_grainst_acc                     => cnveg_carbonstate_inst%matrix_cturnover_grainst_acc_patch               , &
         matrix_cturnover_grainxf_acc                     => cnveg_carbonstate_inst%matrix_cturnover_grainxf_acc_patch               , &

         matrix_nalloc_leaf_acc        => cnveg_nitrogenstate_inst%matrix_nalloc_leaf_acc_patch             , & !
         matrix_nalloc_leafst_acc      => cnveg_nitrogenstate_inst%matrix_nalloc_leafst_acc_patch             , & !
         matrix_nalloc_froot_acc       => cnveg_nitrogenstate_inst%matrix_nalloc_froot_acc_patch             , & !
         matrix_nalloc_frootst_acc     => cnveg_nitrogenstate_inst%matrix_nalloc_frootst_acc_patch             , & !
         matrix_nalloc_livestem_acc    => cnveg_nitrogenstate_inst%matrix_nalloc_livestem_acc_patch             , & !
         matrix_nalloc_livestemst_acc  => cnveg_nitrogenstate_inst%matrix_nalloc_livestemst_acc_patch             , & !
         matrix_nalloc_deadstem_acc    => cnveg_nitrogenstate_inst%matrix_nalloc_deadstem_acc_patch             , & !
         matrix_nalloc_deadstemst_acc  => cnveg_nitrogenstate_inst%matrix_nalloc_deadstemst_acc_patch             , & !
         matrix_nalloc_livecroot_acc   => cnveg_nitrogenstate_inst%matrix_nalloc_livecroot_acc_patch             , & !
         matrix_nalloc_livecrootst_acc => cnveg_nitrogenstate_inst%matrix_nalloc_livecrootst_acc_patch             , & !
         matrix_nalloc_deadcroot_acc   => cnveg_nitrogenstate_inst%matrix_nalloc_deadcroot_acc_patch             , & !
         matrix_nalloc_deadcrootst_acc => cnveg_nitrogenstate_inst%matrix_nalloc_deadcrootst_acc_patch             , & !
         matrix_nalloc_grain_acc       => cnveg_nitrogenstate_inst%matrix_nalloc_grain_acc_patch             , & !
         matrix_nalloc_grainst_acc     => cnveg_nitrogenstate_inst%matrix_nalloc_grainst_acc_patch             , & !

         matrix_ntransfer_leafst_to_leafxf_acc            => cnveg_nitrogenstate_inst%matrix_ntransfer_leafst_to_leafxf_acc_patch          , & !
         matrix_ntransfer_leafxf_to_leaf_acc              => cnveg_nitrogenstate_inst%matrix_ntransfer_leafxf_to_leaf_acc_patch          , & !
         matrix_ntransfer_frootst_to_frootxf_acc          => cnveg_nitrogenstate_inst%matrix_ntransfer_frootst_to_frootxf_acc_patch          , & !
         matrix_ntransfer_frootxf_to_froot_acc            => cnveg_nitrogenstate_inst%matrix_ntransfer_frootxf_to_froot_acc_patch          , & !
         matrix_ntransfer_livestemst_to_livestemxf_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_livestemst_to_livestemxf_acc_patch          , & !
         matrix_ntransfer_livestemxf_to_livestem_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_livestemxf_to_livestem_acc_patch          , & !
         matrix_ntransfer_deadstemst_to_deadstemxf_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch          , & !
         matrix_ntransfer_deadstemxf_to_deadstem_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_deadstemxf_to_deadstem_acc_patch          , & !
         matrix_ntransfer_livecrootst_to_livecrootxf_acc  => cnveg_nitrogenstate_inst%matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch          , & !
         matrix_ntransfer_livecrootxf_to_livecroot_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_livecrootxf_to_livecroot_acc_patch          , & !
         matrix_ntransfer_deadcrootst_to_deadcrootxf_acc  => cnveg_nitrogenstate_inst%matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch          , & !
         matrix_ntransfer_deadcrootxf_to_deadcroot_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch          , & !
         matrix_ntransfer_grainst_to_grainxf_acc          => cnveg_nitrogenstate_inst%matrix_ntransfer_grainst_to_grainxf_acc_patch          , & !
         matrix_ntransfer_grainxf_to_grain_acc            => cnveg_nitrogenstate_inst%matrix_ntransfer_grainxf_to_grain_acc_patch          , & !
         matrix_ntransfer_livestem_to_deadstem_acc        => cnveg_nitrogenstate_inst%matrix_ntransfer_livestem_to_deadstem_acc_patch          , & !
         matrix_ntransfer_livecroot_to_deadcroot_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_livecroot_to_deadcroot_acc_patch          , & !
 
         matrix_ntransfer_retransn_to_leaf_acc        => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_leaf_acc_patch             , & !
         matrix_ntransfer_retransn_to_leafst_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_leafst_acc_patch             , & !
         matrix_ntransfer_retransn_to_froot_acc       => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_froot_acc_patch             , & !
         matrix_ntransfer_retransn_to_frootst_acc     => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_frootst_acc_patch             , & !
         matrix_ntransfer_retransn_to_livestem_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_livestem_acc_patch             , & !
         matrix_ntransfer_retransn_to_livestemst_acc  => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_livestemst_acc_patch             , & !
         matrix_ntransfer_retransn_to_deadstem_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_deadstem_acc_patch             , & !
         matrix_ntransfer_retransn_to_deadstemst_acc  => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_deadstemst_acc_patch             , & !
         matrix_ntransfer_retransn_to_livecroot_acc   => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_livecroot_acc_patch             , & !
         matrix_ntransfer_retransn_to_livecrootst_acc => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_livecrootst_acc_patch             , & !
         matrix_ntransfer_retransn_to_deadcroot_acc   => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_deadcroot_acc_patch             , & !
         matrix_ntransfer_retransn_to_deadcrootst_acc => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_deadcrootst_acc_patch             , & !
         matrix_ntransfer_retransn_to_grain_acc       => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_grain_acc_patch             , & !
         matrix_ntransfer_retransn_to_grainst_acc     => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_grainst_acc_patch             , & !

         matrix_ntransfer_leaf_to_retransn_acc        => cnveg_nitrogenstate_inst%matrix_ntransfer_leaf_to_retransn_acc_patch             , & !
         matrix_ntransfer_froot_to_retransn_acc       => cnveg_nitrogenstate_inst%matrix_ntransfer_froot_to_retransn_acc_patch             , & !
         matrix_ntransfer_livestem_to_retransn_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_livestem_to_retransn_acc_patch             , & !
         matrix_ntransfer_livecroot_to_retransn_acc   => cnveg_nitrogenstate_inst%matrix_ntransfer_livecroot_to_retransn_acc_patch             , & !

         matrix_nturnover_leaf_acc                    => cnveg_nitrogenstate_inst%matrix_nturnover_leaf_acc_patch               , &
         matrix_nturnover_leafst_acc                  => cnveg_nitrogenstate_inst%matrix_nturnover_leafst_acc_patch               , &
         matrix_nturnover_leafxf_acc                  => cnveg_nitrogenstate_inst%matrix_nturnover_leafxf_acc_patch               , &
         matrix_nturnover_froot_acc                   => cnveg_nitrogenstate_inst%matrix_nturnover_froot_acc_patch               , &
         matrix_nturnover_frootst_acc                 => cnveg_nitrogenstate_inst%matrix_nturnover_frootst_acc_patch               , &
         matrix_nturnover_frootxf_acc                 => cnveg_nitrogenstate_inst%matrix_nturnover_frootxf_acc_patch               , &
         matrix_nturnover_livestem_acc                => cnveg_nitrogenstate_inst%matrix_nturnover_livestem_acc_patch               , &
         matrix_nturnover_livestemst_acc              => cnveg_nitrogenstate_inst%matrix_nturnover_livestemst_acc_patch               , &
         matrix_nturnover_livestemxf_acc              => cnveg_nitrogenstate_inst%matrix_nturnover_livestemxf_acc_patch               , &
         matrix_nturnover_deadstem_acc                => cnveg_nitrogenstate_inst%matrix_nturnover_deadstem_acc_patch               , &
         matrix_nturnover_deadstemst_acc              => cnveg_nitrogenstate_inst%matrix_nturnover_deadstemst_acc_patch               , &
         matrix_nturnover_deadstemxf_acc              => cnveg_nitrogenstate_inst%matrix_nturnover_deadstemxf_acc_patch               , &
         matrix_nturnover_livecroot_acc               => cnveg_nitrogenstate_inst%matrix_nturnover_livecroot_acc_patch               , &
         matrix_nturnover_livecrootst_acc             => cnveg_nitrogenstate_inst%matrix_nturnover_livecrootst_acc_patch               , &
         matrix_nturnover_livecrootxf_acc             => cnveg_nitrogenstate_inst%matrix_nturnover_livecrootxf_acc_patch               , &
         matrix_nturnover_deadcroot_acc               => cnveg_nitrogenstate_inst%matrix_nturnover_deadcroot_acc_patch               , &
         matrix_nturnover_deadcrootst_acc             => cnveg_nitrogenstate_inst%matrix_nturnover_deadcrootst_acc_patch               , &
         matrix_nturnover_deadcrootxf_acc             => cnveg_nitrogenstate_inst%matrix_nturnover_deadcrootxf_acc_patch               , &
         matrix_nturnover_grain_acc                   => cnveg_nitrogenstate_inst%matrix_nturnover_grain_acc_patch               , &
         matrix_nturnover_grainst_acc                 => cnveg_nitrogenstate_inst%matrix_nturnover_grainst_acc_patch               , &
         matrix_nturnover_grainxf_acc                 => cnveg_nitrogenstate_inst%matrix_nturnover_grainxf_acc_patch               , &
         matrix_nturnover_retransn_acc                => cnveg_nitrogenstate_inst%matrix_nturnover_retransn_acc_patch               , &

         matrix_alloc                   => cnveg_carbonflux_inst%matrix_alloc_patch                   , & ! Input:      [real(r8) (:,:)]   (gC/gC) input C allocation matrix    		 
         matrix_nalloc                  => cnveg_nitrogenflux_inst%matrix_nalloc_patch                , & ! Input:      [real(r8) (:,:)]   (gC/gC) input N allocation matrix
 
         matrix_phtransfer              => cnveg_carbonflux_inst%matrix_phtransfer_patch              , & ! In/Output:  [real(r8) (:,:,:)] (gC/gC) turnover N transfer matrix
         matrix_gmtransfer              => cnveg_carbonflux_inst%matrix_gmtransfer_patch              , & ! In/Output:
         matrix_fitransfer              => cnveg_carbonflux_inst%matrix_fitransfer_patch              , & ! In/Output:
         matrix_phturnover              => cnveg_carbonflux_inst%matrix_phturnover_patch              , & ! Output:  [real(r8) (:,:,:)] (gC/gC/s) turnover rate diagonal matrix
         matrix_gmturnover              => cnveg_carbonflux_inst%matrix_gmturnover_patch              , & ! Output: 
         matrix_fiturnover              => cnveg_carbonflux_inst%matrix_fiturnover_patch              , & ! Output: 

         matrix_nphtransfer             => cnveg_nitrogenflux_inst%matrix_nphtransfer_patch           , & ! In/Output:  [real(r8) (:,:,:)] (gC/gC) turnover N transfer matrix
         matrix_ngmtransfer             => cnveg_nitrogenflux_inst%matrix_ngmtransfer_patch           , & ! In/Output:
         matrix_nfitransfer             => cnveg_nitrogenflux_inst%matrix_nfitransfer_patch              , & ! In/Output:
         matrix_nphturnover             => cnveg_nitrogenflux_inst%matrix_nphturnover_patch              , & ! Output:  [real(r8) (:,:,:)] (gC/gC/s) turnover rate diagonal matrix
         matrix_ngmturnover             => cnveg_nitrogenflux_inst%matrix_ngmturnover_patch              , & ! Output: 
         matrix_nfiturnover             => cnveg_nitrogenflux_inst%matrix_nfiturnover_patch              , & ! Output: 

         matrix_Cinput                  => cnveg_carbonflux_inst%matrix_Cinput_patch                  , & !

         matrix_Ninput                  => cnveg_nitrogenflux_inst%matrix_Ninput_patch                  & ! Input:      [real(r8) (:) ]    (gN/m2/s) Input N
         )
    !-----------------------------------------------------------------------
     ! set time steps
      dt = real( get_step_size(), r8 )
      secspyear = get_days_per_year() * secspday

!if(bounds%begp .le. 691 .and. bounds%endp .ge. 691)print*,'begin of VEcmatrix',340,soilbiogeochem_nitrogenflux_inst%matrix_input_col(340,1,1)
!if(bounds%begp .le. 1472 .and. bounds%endc .ge. 1472)print*,'begin of VEcmatrix',743,soilbiogeochem_nitrogenflux_inst%matrix_input_col(743,1,1)

      allocate(vegmatrixc_old       (nvegcpool,nvegcpool)) !save it as two dimensional variables in order to track transfers
      allocate(vegmatrixc_new       (nvegcpool))
      allocate(vegmatrixn_old       (nvegnpool,nvegnpool))
      allocate(vegmatrixn_new       (nvegnpool))
      allocate(vegmatrixc_input     (nvegcpool))
      allocate(vegmatrixc_transfer  (nvegcpool,nvegcpool))
      allocate(vegmatrixn_input     (nvegnpool))
      allocate(vegmatrixn_transfer  (nvegnpool,nvegnpool))
      allocate(matrix_calloc_acc    (nvegcpool))
      allocate(matrix_nalloc_acc    (nvegnpool))
      allocate(matrix_ctransfer_acc (nvegcpool,nvegcpool))
      allocate(matrix_ntransfer_acc (nvegnpool,nvegnpool))
!      allocate(vegmatrixc_rt(bounds%begp:bounds%endp,nvegpool))
!      allocate(vegmatrixn_rt(bounds%begp:bounds%endp,nvegpool))
  
!     allocate(matrix_nphturnover(bounds%begp:bounds%endp,nvegpool,nvegpool))
!     allocate(matrix_ngmturnover(bounds%begp:bounds%endp,nvegpool,nvegpool))
!     allocate(matrix_nfitransfer(bounds%begp:bounds%endp,nvegpool+1,nvegpool))

!if(bounds%begp .le. 691 .and. bounds%endp .ge. 691)print*,'here-2 VEcmatrix',340,soilbiogeochem_nitrogenflux_inst%matrix_input_col(340,1,1)
!if(bounds%begp .le. 1472 .and. bounds%endc .ge. 1472)print*,'here-2 VEcmatrix',743,soilbiogeochem_nitrogenflux_inst%matrix_input_col(743,1,1)

!      allocate(vegmatrix_2d(nvegcpool,nvegcpool))
      allocate(AKinvc(nvegcpool,nvegcpool))

!      allocate(vegmatrixn_2d(nvegnpool,nvegnpool))
      allocate(AKinvn(nvegnpool,nvegnpool))

!if(bounds%begp .le. 691 .and. bounds%endp .ge. 691)print*,'here-1 VEcmatrix',340,soilbiogeochem_nitrogenflux_inst%matrix_input_col(340,1,1)
!if(bounds%begp .le. 1472 .and. bounds%endc .ge. 1472)print*,'here-1 VEcmatrix',743,soilbiogeochem_nitrogenflux_inst%matrix_input_col(743,1,1)
      vegmatrixc_old     (:,:)   = 0._r8
      vegmatrixn_old     (:,:)   = 0._r8
      vegmatrixc_new     (:)     = 0._r8
      vegmatrixn_new     (:)     = 0._r8
      vegmatrixc_input   (:)     = 0._r8
      vegmatrixc_transfer(:,:)   = 0._r8
      vegmatrixn_input   (:)     = 0._r8
      vegmatrixn_transfer(:,:)   = 0._r8
      matrix_calloc_acc    (:)   = 0._r8
      matrix_nalloc_acc    (:)   = 0._r8
      matrix_ctransfer_acc (:,:) = 0._r8
      matrix_ntransfer_acc (:,:) = 0._r8
      AKinvc (:,:) = 0._r8
      AKinvn (:,:) = 0._r8
      rowonec(:) = 1._r8
      rowonen(:) = 1._r8
      
!if(bounds%begp .le. 691 .and. bounds%endp .ge. 691)print*,'here-0.5 VEcmatrix',340,soilbiogeochem_nitrogenflux_inst%matrix_input_col(340,1,1)
!if(bounds%begp .le. 1472 .and. bounds%endc .ge. 1472)print*,'here-0.5 VEcmatrix',743,soilbiogeochem_nitrogenflux_inst%matrix_input_col(743,1,1)

      counter = counter + dt
       if (counter >=1.0 * secspyear) then ! link to the recycling span of climate forcing
          end_of_year = .true.
          counter = 0._r8
       else
          end_of_year = .false.
       end if
!if(bounds%begp .le. 691 .and. bounds%endp .ge. 691)print*,'here-0.25 VEcmatrix',340,soilbiogeochem_nitrogenflux_inst%matrix_input_col(340,1,1)
!if(bounds%begp .le. 1472 .and. bounds%endc .ge. 1472)print*,'here-0.25 VEcmatrix',743,soilbiogeochem_nitrogenflux_inst%matrix_input_col(743,1,1)

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)
!         print*,'p,type',p,ivt(p)
!         if(abs(grc%latdeg(patch%gridcell(p))+40.0) .le. 0.01 .and. abs(grc%londeg(patch%gridcell(p))-150) .le. 0.01)then
!             print*,'really begin of vegcmatrix',c,soilbiogeochem_nitrogenflux_inst%matrix_input_col(c,1,1)
!         end if
!         vegmatrix_2d(1:nvegcpool,1:nvegcpool)=0._r8
!         vegmatrixn_2d(1:nvegnpool,1:nvegnpool)=0._r8
!         AKinv(1:nvegcpool,1:nvegcpool)=0._r8
!         AKinvn(1:nvegnpool,1:nvegnpool)=0._r8
!         if(abs(grc%latdeg(patch%gridcell(p))+40.0) .le. 0.01 .and. abs(grc%londeg(patch%gridcell(p))-150) .le. 0.01)then
!             print*,'begin of vegcmatrix',c,soilbiogeochem_nitrogenflux_inst%matrix_input_col(c,1,1)
!         end if
         vegmatrixc_old(ileaf        ,ileaf)         = leafc(p)
         vegmatrixc_old(ileaf_st     ,ileaf_st)      = leafc_storage(p)
         vegmatrixc_old(ileaf_xf     ,ileaf_xf)      = leafc_xfer(p)
         vegmatrixc_old(ifroot       ,ifroot)        = frootc(p)
         vegmatrixc_old(ifroot_st    ,ifroot_st)     = frootc_storage(p)
         vegmatrixc_old(ifroot_xf    ,ifroot_xf)     = frootc_xfer(p)
         vegmatrixc_old(ilivestem    ,ilivestem)     = livestemc(p)
         vegmatrixc_old(ilivestem_st ,ilivestem_st)  = livestemc_storage(p)
         vegmatrixc_old(ilivestem_xf ,ilivestem_xf)  = livestemc_xfer(p)
         vegmatrixc_old(ideadstem    ,ideadstem)     = deadstemc(p)
         vegmatrixc_old(ideadstem_st ,ideadstem_st)  = deadstemc_storage(p)
         vegmatrixc_old(ideadstem_xf ,ideadstem_xf)  = deadstemc_xfer(p)
         vegmatrixc_old(ilivecroot   ,ilivecroot)    = livecrootc(p)
         vegmatrixc_old(ilivecroot_st,ilivecroot_st) = livecrootc_storage(p)
         vegmatrixc_old(ilivecroot_xf,ilivecroot_xf) = livecrootc_xfer(p)
         vegmatrixc_old(ideadcroot   ,ideadcroot)    = deadcrootc(p)
         vegmatrixc_old(ideadcroot_st,ideadcroot_st) = deadcrootc_storage(p)
         vegmatrixc_old(ideadcroot_xf,ideadcroot_xf) = deadcrootc_xfer(p)
!         if(p .eq. 5228)print*,'before vegmatrixc_old',vegmatrixc_old(p,igrain),vegmatrixc_old(p,igrain_st),vegmatrixc_old(p,igrain_xf),grainc(p),grainc_storage(p),grainc_xfer(p)
!         if(abs(grc%latdeg(patch%gridcell(p))+40.0) .le. 0.01 .and. abs(grc%londeg(patch%gridcell(p))-150) .le. 0.01)then
!             print*,'here0 in vegcmatrix',c,soilbiogeochem_nitrogenflux_inst%matrix_input_col(c,1,1)
!         end if
         if(ivt(p) >= npcropmin)then
            vegmatrixc_old(igrain   ,igrain)     = grainc(p)
            vegmatrixc_old(igrain_st,igrain)     = grainc_storage(p)
            vegmatrixc_old(igrain_xf,igrain)     = grainc_xfer(p)
         end if
!         if(p .eq. 5228)print*,'after vegmatrixc_old',vegmatrixc_old(p,igrain),vegmatrixc_old(p,igrain_st),vegmatrixc_old(p,igrain_xf)
 
         vegmatrixn_old(ileaf        ,ileaf)         = leafn(p)
         vegmatrixn_old(ileaf_st     ,ileaf_st)      = leafn_storage(p)
         vegmatrixn_old(ileaf_xf     ,ileaf_xf)      = leafn_xfer(p)
         vegmatrixn_old(ifroot       ,ifroot)        = frootn(p)
         vegmatrixn_old(ifroot_st    ,ifroot_st)     = frootn_storage(p)
         vegmatrixn_old(ifroot_xf    ,ifroot_xf)     = frootn_xfer(p)
         vegmatrixn_old(ilivestem    ,ilivestem)     = livestemn(p)
         vegmatrixn_old(ilivestem_st ,ilivestem_st)  = livestemn_storage(p)
         vegmatrixn_old(ilivestem_xf ,ilivestem_xf)  = livestemn_xfer(p)
         vegmatrixn_old(ideadstem    ,ideadstem)     = deadstemn(p)
         vegmatrixn_old(ideadstem_st ,ideadstem_st)  = deadstemn_storage(p)
         vegmatrixn_old(ideadstem_xf ,ideadstem_xf)  = deadstemn_xfer(p)
         vegmatrixn_old(ilivecroot   ,ilivecroot)    = livecrootn(p)
         vegmatrixn_old(ilivecroot_st,ilivecroot_st) = livecrootn_storage(p)
         vegmatrixn_old(ilivecroot_xf,ilivecroot_xf) = livecrootn_xfer(p)
         vegmatrixn_old(ideadcroot   ,ideadcroot)    = deadcrootn(p)
         vegmatrixn_old(ideadcroot_st,ideadcroot_st) = deadcrootn_storage(p)
         vegmatrixn_old(ideadcroot_xf,ideadcroot_xf) = deadcrootn_xfer(p)
         if(ivt(p) >= npcropmin)then
            vegmatrixn_old(igrain    ,igrain)        = grainn(p)
            vegmatrixn_old(igrain_st ,igrain_st)     = grainn_storage(p)
            vegmatrixn_old(igrain_xf ,igrain_xf)     = grainn_xfer(p)
         end if
         vegmatrixn_old(iretransn    ,iretransn)     = retransn(p)

         if (is_beg_curr_year() .or. is_first_step_of_this_run_segment() )then
            leafc0(p)                = max(leafc(p),1.e-15_r8)
            leafc0_storage(p)        = max(leafc_storage(p), 1.e-15_r8)
            leafc0_xfer(p)           = max(leafc_xfer(p), 1.e-15_r8)
            frootc0(p)               = max(frootc(p), 1.e-15_r8)
            frootc0_storage(p)       = max(frootc_storage(p), 1.e-15_r8)
            frootc0_xfer(p)          = max(frootc_xfer(p), 1.e-15_r8)
            livestemc0(p)            = max(livestemc(p), 1.e-15_r8)
            livestemc0_storage(p)    = max(livestemc_storage(p), 1.e-15_r8)
            livestemc0_xfer(p)       = max(livestemc_xfer(p), 1.e-15_r8)
            deadstemc0(p)            = max(deadstemc(p), 1.e-15_r8)
            deadstemc0_storage(p)    = max(deadstemc_storage(p), 1.e-15_r8)
            deadstemc0_xfer(p)       = max(deadstemc_xfer(p), 1.e-15_r8)
            livecrootc0(p)           = max(livecrootc(p), 1.e-15_r8)
            livecrootc0_storage(p)   = max(livecrootc_storage(p), 1.e-15_r8)
            livecrootc0_xfer(p)      = max(livecrootc_xfer(p), 1.e-15_r8)
            deadcrootc0(p)           = max(deadcrootc(p), 1.e-15_r8)
            deadcrootc0_storage(p)   = max(deadcrootc_storage(p), 1.e-15_r8)
            deadcrootc0_xfer(p)      = max(deadcrootc_xfer(p), 1.e-15_r8)
            if(ivt(p) >= npcropmin)then
            grainc0(p)               = max(grainc(p), 1.e-15_r8)
            grainc0_storage(p)       = max(grainc_storage(p), 1.e-15_r8)
            grainc0_xfer(p)          = max(grainc_xfer(p), 1.e-15_r8)
            end if
            leafn0(p)                = max(leafn(p),1.e-15_r8)
            leafn0_storage(p)        = max(leafn_storage(p), 1.e-15_r8)
            leafn0_xfer(p)           = max(leafn_xfer(p), 1.e-15_r8)
            frootn0(p)               = max(frootn(p), 1.e-15_r8)
            frootn0_storage(p)       = max(frootn_storage(p), 1.e-15_r8)
            frootn0_xfer(p)          = max(frootn_xfer(p), 1.e-15_r8)
            livestemn0(p)            = max(livestemn(p), 1.e-15_r8)
            livestemn0_storage(p)    = max(livestemn_storage(p), 1.e-15_r8)
            livestemn0_xfer(p)       = max(livestemn_xfer(p), 1.e-15_r8)
            deadstemn0(p)            = max(deadstemn(p), 1.e-15_r8)
            deadstemn0_storage(p)    = max(deadstemn_storage(p), 1.e-15_r8)
            deadstemn0_xfer(p)       = max(deadstemn_xfer(p), 1.e-15_r8)
            livecrootn0(p)           = max(livecrootn(p), 1.e-15_r8)
            livecrootn0_storage(p)   = max(livecrootn_storage(p), 1.e-15_r8)
            livecrootn0_xfer(p)      = max(livecrootn_xfer(p), 1.e-15_r8)
            deadcrootn0(p)           = max(deadcrootn(p), 1.e-15_r8)
            deadcrootn0_storage(p)   = max(deadcrootn_storage(p), 1.e-15_r8)
            deadcrootn0_xfer(p)      = max(deadcrootn_xfer(p), 1.e-15_r8)
            if(ivt(p) >= npcropmin)then
            grainn0(p)               = max(grainn(p), 1.e-15_r8)
            grainn0_storage(p)       = max(grainn_storage(p), 1.e-15_r8)
            grainn0_xfer(p)          = max(grainn_xfer(p), 1.e-15_r8)
            end if
            retransn0(p)             = max(retransn(p), 1.e-15_r8)
         end if
!
!         if(abs(grc%latdeg(patch%gridcell(p))+40.0) .le. 0.01 .and. abs(grc%londeg(patch%gridcell(p))-150) .le. 0.01)then
!             print*,'here1 in vegcmatrix',c,soilbiogeochem_nitrogenflux_inst%matrix_input_col(c,1,1)
!         end if
            
!          vegmatrix_2d(ileaf,ileaf)               = leafc(p)         / leafc0(p)
!          vegmatrix_2d(ileaf_st,ileaf_st)         = leafc_storage(p) / leafc0_storage(p)
!          vegmatrix_2d(ileaf_xf,ileaf_xf)         = leafc_xfer(p)    / leafc0_xfer(p)
!          vegmatrix_2d(ifroot,ifroot)             = frootc(p)        / frootc0(p)
!          vegmatrix_2d(ifroot_st,ifroot_st)       = frootc_storage(p)/ frootc0_storage(p)
!          vegmatrix_2d(ifroot_xf,ifroot_xf)       = frootc_xfer(p)   / frootc0_xfer(p)
!          vegmatrix_2d(ilivestem,ilivestem)       = livestemc(p)     / livestemc0(p)
!          vegmatrix_2d(ilivestem_st,ilivestem_st) = livestemc_storage(p) / livestemc0_storage(p)
!          vegmatrix_2d(ilivestem_xf,ilivestem_xf) = livestemc_xfer(p)    / livestemc0_xfer(p)
!          vegmatrix_2d(ideadstem,ideadstem)       = deadstemc(p)     / deadstemc0(p)
!          vegmatrix_2d(ideadstem_st,ideadstem_st) = deadstemc_storage(p)/ deadstemc0_storage(p)
!          vegmatrix_2d(ideadstem_xf,ideadstem_xf) = deadstemc_xfer(p) / deadstemc0_xfer(p)
!          vegmatrix_2d(ilivecroot,ilivecroot)     = livecrootc(p)    / livecrootc0(p)
!          vegmatrix_2d(ilivecroot_st,ilivecroot_st) = livecrootc_storage(p)/livecrootc0_storage(p)
!          vegmatrix_2d(ilivecroot_xf,ilivecroot_xf) = livecrootc_xfer(p)/ livecrootc0_xfer(p)
!          vegmatrix_2d(ideadcroot,ideadcroot)       = deadcrootc(p)  / deadcrootc0(p)
!          vegmatrix_2d(ideadcroot_st,ideadcroot_st) = deadcrootc_storage(p)/deadcrootc0_storage(p)
!          vegmatrix_2d(ideadcroot_xf,ideadcroot_xf) = deadcrootc_xfer(p)/deadcrootc0_xfer(p)
!          if(ivt(p) >= npcropmin)then
!             vegmatrix_2d(igrain,igrain)             = grainc(p)        / grainc0(p)
!             vegmatrix_2d(igrain_st,igrain_st)       = grainc_storage(p)/ grainc0_storage(p)
!             vegmatrix_2d(igrain_xf,igrain_xf)       = grainc_xfer(p)   / grainc0_xfer(p)
!          end if

!          vegmatrixn_2d(ileaf,ileaf)               = leafn(p)         / leafn0(p)
!          vegmatrixn_2d(ileaf_st,ileaf_st)         = leafn_storage(p) / leafn0_storage(p)
!          vegmatrixn_2d(ileaf_xf,ileaf_xf)         = leafn_xfer(p)    / leafn0_xfer(p)
!          vegmatrixn_2d(ifroot,ifroot)             = frootn(p)        / frootn0(p)
!          vegmatrixn_2d(ifroot_st,ifroot_st)       = frootn_storage(p)/ frootn0_storage(p)
!          vegmatrixn_2d(ifroot_xf,ifroot_xf)       = frootn_xfer(p)   / frootn0_xfer(p)
!          vegmatrixn_2d(ilivestem,ilivestem)       = livestemn(p)     / livestemn0(p)
!          vegmatrixn_2d(ilivestem_st,ilivestem_st) = livestemn_storage(p) / livestemn0_storage(p)
!          vegmatrixn_2d(ilivestem_xf,ilivestem_xf) = livestemn_xfer(p)    / livestemn0_xfer(p)
!          vegmatrixn_2d(ideadstem,ideadstem)       = deadstemn(p)     / deadstemn0(p)
!          vegmatrixn_2d(ideadstem_st,ideadstem_st) = deadstemn_storage(p)/ deadstemn0_storage(p)
!          vegmatrixn_2d(ideadstem_xf,ideadstem_xf) = deadstemn_xfer(p) / deadstemn0_xfer(p)
!          vegmatrixn_2d(ilivecroot,ilivecroot)     = livecrootn(p)    / livecrootn0(p)
!          vegmatrixn_2d(ilivecroot_st,ilivecroot_st) = livecrootn_storage(p)/livecrootn0_storage(p)
!          vegmatrixn_2d(ilivecroot_xf,ilivecroot_xf) = livecrootn_xfer(p)/livecrootn0_xfer(p)
!          vegmatrixn_2d(ideadcroot,ideadcroot)       = deadcrootn(p)  / deadcrootn0(p)
!          vegmatrixn_2d(ideadcroot_st,ideadcroot_st) = deadcrootn_storage(p)/deadcrootn0_storage(p)
!          vegmatrixn_2d(ideadcroot_xf,ideadcroot_xf) = deadcrootn_xfer(p)/deadcrootn0_xfer(p)
!          if(ivt(p) >= npcropmin)then
!             vegmatrixn_2d(igrain,igrain)             = grainn(p)        / grainn0(p)
!             vegmatrixn_2d(igrain_st,igrain_st)       = grainn_storage(p)/ grainn0_storage(p)
!             vegmatrixn_2d(igrain_xf,igrain_xf)       = grainn_xfer(p)   / grainn0_xfer(p)
!          end if
!         if(abs(grc%latdeg(patch%gridcell(p))+40.0) .le. 0.01 .and. abs(grc%londeg(patch%gridcell(p))-150) .le. 0.01)then
!             print*,'here2 in vegcmatrix',c,soilbiogeochem_nitrogenflux_inst%matrix_input_col(c,1,1)
!         end if
!
         matrix_phturnover(p,1:nvegcpool,1:nvegcpool) = 0._r8
         matrix_gmturnover(p,1:nvegcpool,1:nvegcpool) = 0._r8
         matrix_fiturnover(p,1:nvegcpool,1:nvegcpool) = 0._r8
         matrix_nphturnover(p,1:nvegcpool,1:nvegcpool) = 0._r8
         matrix_ngmturnover(p,1:nvegcpool,1:nvegcpool) = 0._r8
         matrix_nfiturnover(p,1:nvegcpool,1:nvegcpool) = 0._r8

!         matrix_nphtransfer_curr(p,:,:) = matrix_phtransfer(p,:,:)
!         matrix_nphtransfer_curr(p,ileaf_xf,ileaf_st)     = matrix_nphtransfer(p,ileaf_xf,ileaf_st)
!         matrix_nphtransfer_curr(p,ifroot_xf,ifroot_st)   = matrix_nphtransfer(p,ifroot_xf,ifroot_st)
   
!         matrix_nphtransfer_curr(p,ideadstem,ilivestem)    =  matrix_nphtransfer(p,ideadstem,ilivestem)
!         matrix_nphtransfer_curr(p,ideadcroot,ilivecroot)  =  matrix_nphtransfer(p,ideadcroot,ilivecroot)

!         if(p .eq. 12 .or. p .eq. 13)write(517,"(A,324F10.2)"),'tranfer matrix',p,matrix_phtransfer(p,1:nvegcpool,1:nvegcpool) * dt
!         if(p .eq. 12 .or. p .eq. 13)write(518,"(A,18F10.2)"),'pool sizes',p,vegmatrixc_old(p,1:nvegcpool)
!         if(p .eq. 12 .or. p .eq. 13)write(518,"(A,18F10.2)"),'pool sizes',p,vegmatrixc_old(p,1:nvegcpool)

         do j=1,nvegcpool
            do i=1,nvegcpool+1
               if(i .ne. j)then
                  matrix_phturnover(p,j,j) = matrix_phturnover(p,j,j) + matrix_phtransfer(p,i,j)
                  matrix_gmturnover(p,j,j) = matrix_gmturnover(p,j,j) + matrix_gmtransfer(p,i,j)
                  matrix_fiturnover(p,j,j) = matrix_fiturnover(p,j,j) + matrix_fitransfer(p,i,j)
               end if
            end do
         end do
         
         do j=1,nvegnpool
            do i=1,nvegnpool+1
               if(i .ne. j)then
                 matrix_nphturnover(p,j,j) = matrix_nphturnover(p,j,j) + matrix_nphtransfer(p,i,j)
                 matrix_ngmturnover(p,j,j) = matrix_ngmturnover(p,j,j) + matrix_ngmtransfer(p,i,j)
                 matrix_nfiturnover(p,j,j) = matrix_nfiturnover(p,j,j) + matrix_nfitransfer(p,i,j)
               end if
            end do
         end do
!         if( p .eq. 2)print*,'phtransfer',matrix_phtransfer(p,ideadcroot,:)
!         if(abs(grc%latdeg(patch%gridcell(p))+40.0) .le. 0.01 .and. abs(grc%londeg(patch%gridcell(p))-150) .le. 0.01)then
!             print*,'here3 in vegcmatrix',c,soilbiogeochem_nitrogenflux_inst%matrix_input_col(c,1,1)
!         end if

         do j=1,nvegcpool
            if(matrix_phturnover(p,j,j) .ne. 0)then
               matrix_phtransfer(p,:,j) = matrix_phtransfer(p,:,j) / matrix_phturnover(p,j,j)
            else
               matrix_phtransfer(p,:,j) = 0._r8
            end if
            matrix_phtransfer(p,j,j) = -1._r8

            if(matrix_gmturnover(p,j,j) .ne. 0)then
               matrix_gmtransfer(p,:,j) = matrix_gmtransfer(p,:,j) / matrix_gmturnover(p,j,j)
            else
               matrix_gmtransfer(p,:,j) = 0._r8
            end if
            matrix_gmtransfer(p,j,j) = -1._r8

            if(matrix_fiturnover(p,j,j) .ne. 0)then
               matrix_fitransfer(p,:,j) = matrix_fitransfer(p,:,j) / matrix_fiturnover(p,j,j)
            else
               matrix_fitransfer(p,:,j) = 0._r8
            end if
            matrix_fitransfer(p,j,j) = -1._r8
         end do
!N
!         if(p .eq. 8)print*,'before generate turnover',matrix_ngmturnover(p,ideadstem,ideadstem),matrix_ngmtransfer(p,ioutn,ideadstem)
         do j=1,nvegnpool
            if(matrix_nphturnover(p,j,j) .ne. 0)then
               matrix_nphtransfer(p,:,j) = matrix_nphtransfer(p,:,j) / matrix_nphturnover(p,j,j)
            else
               matrix_nphtransfer(p,:,j) = 0._r8
            end if
            matrix_nphtransfer(p,j,j) = -1._r8

            if(matrix_ngmturnover(p,j,j) .ne. 0)then
               matrix_ngmtransfer(p,:,j) = matrix_ngmtransfer(p,:,j) / matrix_ngmturnover(p,j,j)
            else
               matrix_ngmtransfer(p,:,j) = 0._r8
            end if
            matrix_ngmtransfer(p,j,j) = -1._r8
!
            if(matrix_nfiturnover(p,j,j) .ne. 0)then
               matrix_nfitransfer(p,:,j) = matrix_nfitransfer(p,:,j) / matrix_nfiturnover(p,j,j)
            else
               matrix_nfitransfer(p,:,j) = 0._r8
            end if
            matrix_nfitransfer(p,j,j) = -1._r8
         end do
!         if(p .eq. 5228)print *, 'before1NNN',matrix_alloc(p,:) * matrix_Cinput(p) * dt, grainc(p),grainc_storage(p),grainc_xfer(p),vegmatrixc_old(p,igrain),vegmatrixc_old(p,igrain_st),vegmatrixc_old(p,igrain_xf)
!         if(p .eq. 5228)then
!             tmp=(matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!             print*,'before matrix_new,phenology',cnveg_carbonstate_inst%deadcrootc_patch(p),tmp(ideadcroot)
!             do i=1,nvegcpool
!                print*,'pool,i',i,matrix_phtransfer(p,i,:),matrix_phturnover(p,i,i),vegmatrixc_old(p,i)
!             end do
!             tmp=(matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!             print*,'before matrix_new,gap mortality',cnveg_carbonstate_inst%deadcrootc_patch(p),tmp(ideadcroot)
!             tmp=(matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!             print*,'before matrix_new,fire',cnveg_carbonstate_inst%deadcrootc_patch(p),tmp(ideadcroot)
!         end if
!         if(p .eq. 8428)then
!            print*,'before matrix_new',vegmatrixc_new(p,ideadcroot)
!         end if
!         if(abs(grc%latdeg(patch%gridcell(p))+40.0) .le. 0.01 .and. abs(grc%londeg(patch%gridcell(p))-150) .le. 0.01)then
!            print*,'before matrix,leafc',p,vegmatrixc_old(p,ileaf),vegmatrixc_old(p,ileaf_st),vegmatrixc_old(p,ileaf_xf),matrix_Cinput(p),matrix_alloc(p,ileaf),matrix_alloc(p,ileaf_st),matrix_alloc(p,ileaf_xf),matrix_phtransfer(p,ileaf,ileaf_xf),matrix_phturnover(p,ileaf_xf,ileaf_xf),matrix_phtransfer(p,ileaf_xf,ileaf_st),matrix_phturnover(p,ileaf_st,ileaf_st),matrix_phturnover(p,ileaf_xf,:),c,soilbiogeochem_nitrogenflux_inst%matrix_input_col(c,1,1)
!         end if
!         print*,'vegmatrixc_old',p,matmul(vegmatrixc_old(:,:),rowonec(:))
         vegmatrixc_input(:) = matrix_alloc(p,:) * matrix_Cinput(p) * dt
         vegmatrixc_transfer(:,:) = (matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc_old(:,:)) &
                                   + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc_old(:,:)) &
                                   + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc_old(:,:))) * dt
         vegmatrixc_new(:) = matmul(vegmatrixc_old(:,:),rowonec(:))  +  vegmatrixc_input(:) + matmul(vegmatrixc_transfer(:,:),rowonec(:)) 
!         if(abs(grc%latdeg(patch%gridcell(p))+40.0) .le. 0.01 .and. abs(grc%londeg(patch%gridcell(p))-150) .le. 0.01)then
!            print*,'after matrix,leafc',p,vegmatrixc_new(p,ileaf),vegmatrixc_new(p,ileaf_st),vegmatrixc_new(p,ileaf_xf)
!         end if
         !CiPEHR
!         if(p .eq. 12 .or. p .eq. 13)write(517,"(A,I,324E17.9)"),'tranfer matrix',p,vegmatrixc_transfer(1:nvegcpool,1:nvegcpool)
!         if(p .eq. 12 .or. p .eq. 13)write(518,"(A,I,18E17.9)"),'vegmatrixc_old',p,matmul(vegmatrixc_old(:,:),rowonec(:))
!         if(p .eq. 12 .or. p .eq. 13)write(513,"(A,I,19E17.9)"),'Cinput,alloc',p,matrix_Cinput(p)*dt,matrix_alloc(p,:) 
         !SPRUCE
!         if(p .eq. 2 .or. p .eq. 8)write(517,"(A,I,324E17.9)"),'tranfer matrix',p,vegmatrixc_transfer(1:nvegcpool,1:nvegcpool)
!         if(p .eq. 2 .or. p .eq. 8)write(518,"(A,I,18E17.9)"),'vegmatrixc_old',p,matmul(vegmatrixc_old(:,:),rowonec(:))
!         if(p .eq. 2 .or. p .eq. 8)write(513,"(A,I,19E17.9)"),'Cinput,alloc',p,matrix_Cinput(p)*dt,matrix_alloc(p,:) 
         !SEV
!         if(p .eq. 1 .or. p .eq. 15)write(517,"(A,I,324E17.9)"),'tranfer matrix',p,vegmatrixc_transfer(1:nvegcpool,1:nvegcpool)
!         if(p .eq. 1 .or. p .eq. 15)write(518,"(A,I,18E17.9)"),'vegmatrixc_old',p,matmul(vegmatrixc_old(:,:),rowonec(:))
!         if(p .eq. 1 .or. p .eq. 15)write(513,"(A,I,19E17.9)"),'Cinput,alloc',p,matrix_Cinput(p)*dt,matrix_alloc(p,:) 
                            
         if(p .eq. -9999)then
            !print*,'before matrix_old deadcrootc',p,vegmatrixc_old(p,ideadcroot)
!             print*,'before matrix_old retransn',p,vegmatrixn_old(iretransn,iretransn)
!             print*,'before matrix_old frootn',p,vegmatrixn_old(ifroot,ifroot)
!            print*,'before matrix_old leafn',p,vegmatrixn_old(p,ileaf)
!         print*,'before matrix_old',p,vegmatrixn_old(p,:)
!         end if
!N matrix
!            print*,'Cinput to deadcrootn',matrix_Cinput(p)*matrix_alloc(p,ideadcroot)*dt
!            print*,'Ninput to grainn',matrix_Ninput(p)*matrix_nalloc(p,ileaf)*dt,matrix_Ninput(p),matrix_nalloc(p,ileaf)*dt
!            print*,'Ninput to deadstemn',matrix_Ninput(p)*matrix_nalloc(p,ideadstem)*dt,matrix_Ninput(p),matrix_nalloc(p,ideadstem)*dt
!            print*,'Transfer to deadcrootc phenology from livecroot',matrix_phtransfer(p,ideadcroot,ilivecroot)*matrix_phturnover(p,ilivecroot,ilivecroot) &
!                                                                 * vegmatrixc_old(p,ilivecroot) * dt
!            print*,'Transfer to deadcrootc phenology from transfer pool',matrix_phtransfer(p,ideadcroot,ideadcroot_xf)*matrix_phturnover(p,ideadcroot_xf,ideadcroot_xf) &
!  * vegmatrixc_old(p,ideadcroot_xf) * dt,vegmatrixc_old(p,ideadcroot_xf),matrix_phturnover(p,ideadcroot_xf,ideadcroot_xf),matrix_phtransfer(p,ideadcroot,ideadcroot_xf) 
!            print*,'Transfer to leafc phenology from transfer pool',matrix_phtransfer(p,ileaf,ileaf_xf)*matrix_phturnover(p,ileaf_xf,ileaf_xf) &
!  * vegmatrixc_old(p,ileaf_xf) * dt,vegmatrixc_old(p,ileaf_xf),matrix_phturnover(p,ileaf_xf,ileaf_xf),matrix_phtransfer(p,ileaf,ileaf_xf) 
!            print*,'Transfer from deadcrootn phenology',matrix_phtransfer(p,ideadcroot,ideadcroot) * matrix_phturnover(p,ideadcroot,ideadcroot) &
!                                                                    * vegmatrixc_old(p,ideadcroot) * dt
!            print*,'Transfer from deadstemn_xfer phenology',matrix_nphtransfer(p,ideadstem,ideadstem_xf) * matrix_nphturnover(p,ideadstem_xf,ideadstem_xf) &
!                                                                    * vegmatrixn_old(p,ideadstem_xf) * dt
!            print*,'Transfer from livestemn phenology',matrix_nphtransfer(p,ideadstem,ilivestem) * matrix_nphturnover(p,ilivestem,ilivestem) &
!                                                                    * vegmatrixn_old(p,ilivestem) * dt
!            print*,'Transfer from leafn to retransn phenology',matrix_nphtransfer(p,iretransn,ileaf) * matrix_nphturnover(p,ileaf,ileaf) &
!                                                                    * vegmatrixn_old(p,ileaf) * dt
!            print*,'Transfer from frootn_xfer to frootn phenology',matrix_nphtransfer(p,ifroot,ifroot_xf) * matrix_nphturnover(p,ifroot_xf,ifroot_xf) &
!                                                                    * vegmatrixn_old(p,ifroot_xf) * dt
!            print*,'Transfer from npool to frootn phenology',matrix_Ninput(p)*matrix_nalloc(p,ifroot)*dt
!            print*,'Transfer from retransn to frootn phenology',matrix_nphtransfer(p,ifroot,iretransn) * matrix_nphturnover(p,iretransn,iretransn) &
!                                                                    * vegmatrixn_old(p,iretransn) * dt
!            print*,'Transfer from frootn to retransn phenology',matrix_nphtransfer(p,iretransn,ifroot) * matrix_nphturnover(p,ifroot,ifroot) &
!                                                                    * vegmatrixn_old(p,ifroot) * dt
!            print*,'Transfer from frootn to litter phenology',sum(matrix_nphtransfer(p,1:nvegnpool,ifroot)) * matrix_nphturnover(p,ifroot,ifroot) &
!                                                                    * vegmatrixn_old(p,ifroot) * dt 
!            print*,'Transfer from livestemn to retransn phenology',matrix_nphtransfer(p,iretransn,ilivestem) * matrix_nphturnover(p,ilivestem,ilivestem) &
!                                                                    * vegmatrixn_old(p,ilivestem) * dt
!            print*,'Transfer from livecrootn to retransn phenology',matrix_nphtransfer(p,iretransn,ilivecroot) * matrix_nphturnover(p,ilivecroot,ilivecroot) &
!                                                                    * vegmatrixn_old(p,ilivecroot) * dt
!            print*,'Transfer from retransn to free',matrix_nphtransfer(p,ioutn,iretransn) * matrix_nphturnover(p,iretransn,iretransn) &
!                                                                    * vegmatrixn_old(p,iretransn) * dt
!            print*,'Transfer from retransn to veg',sum(matrix_nphtransfer(p,1:nvegnpool-1,iretransn)) * matrix_nphturnover(p,iretransn,iretransn) &
!                                                                    * vegmatrixn_old(p,iretransn) * dt
!            print*,'Transfer from retransn to out',sum(matrix_nphtransfer(p,1:nvegnpool,iretransn)) * matrix_nphturnover(p,iretransn,iretransn) &
!                                                                    * vegmatrixn_old(p,iretransn) * dt
            
!            print*,'Transfer from deadcrootn gap mort',matrix_gmtransfer(p,ideadcroot,ideadcroot) * matrix_gmturnover(p,ideadcroot,ideadcroot) &
!                                                                    * vegmatrixc_old(p,ideadcroot) * dt
!            print*,'Transfer from deadcrootn fire',matrix_fitransfer(p,ideadcroot,ideadcroot) * matrix_fiturnover(p,ideadcroot,ideadcroot) &
!                                                                    * vegmatrixc_old(p,ideadcroot) * dt
!            print*,'Transfer to deadcrootn fire',matrix_fitransfer(p,ideadcroot,ilivecroot) * matrix_fiturnover(p,ideadcroot,ilivecroot) &
!                                                                    * vegmatrixc_old(p,ilivecroot) * dt
!            print*,'Transfer from leafn gap mort',matrix_ngmtransfer(p,ileaf,ileaf) * matrix_ngmturnover(p,ileaf,ileaf) &
!                                                                    * vegmatrixn_old(p,ileaf) * dt
!            print*,'Transfer from leafn fire',matrix_fitransfer(p,ileaf,ileaf) * matrix_fiturnover(p,ileaf,ileaf) &
!                                                                    * vegmatrixn_old(p,ileaf) * dt
!            print*,'Transfer from deadstemn gap mort',matrix_ngmtransfer(p,ioutn,ideadstem) * matrix_ngmturnover(p,ideadstem,ideadstem) &
!                                                                    * vegmatrixn_old(p,ideadstem) * dt, matrix_ngmtransfer(p,ioutn,ideadstem)
!            print*,'Transfer from deadstemn fire',matrix_nfitransfer(p,ideadstem,ideadstem) * matrix_nfiturnover(p,ideadstem,ideadstem) &
!                                                                    * vegmatrixn_old(p,ideadstem) * dt
!            print*,'Transfer from livestemn fire',matrix_nfitransfer(p,ideadstem,ilivestem) * matrix_nfiturnover(p,ideadstem,ilivestem) &
!                                                                    * vegmatrixn_old(p,ilivestem) * dt
         end if
!         if(p .eq. 16)print*,'input retransn',matrix_nalloc(p,iretransn) * matrix_Ninput(p) * dt
!         tmp=(matmul(matmul(matrix_nphtransfer(p,1:nvegnpool,:),matrix_nphturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!         if(p .eq. 16)print*,'nphtrans retransn',tmp(iretransn)
!         if(p .eq. 16)print*,'nphtrans frootn',tmp(ifroot)
!         tmp=(matmul(matmul(matrix_ngmtransfer(p,1:nvegnpool,:),matrix_ngmturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!         if(p .eq. 16)print*,'ngmtrans retransn',tmp(iretransn)
!         if(p .eq. 16)print*,'ngmtrans frootn',tmp(ifroot)
!         tmp=(matmul(matmul(matrix_nfitransfer(p,1:nvegnpool,:),matrix_nfiturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!         if(p .eq. 16)print*,'nfitrans retransn',tmp(iretransn)
!         if(p .eq. 16)print*,'nfitrans frootn',tmp(ifroot)
         vegmatrixn_input(:) = matrix_nalloc(p,:) * matrix_Ninput(p) * dt
         vegmatrixn_transfer(:,:) = (matmul(matmul(matrix_nphtransfer(p,1:nvegnpool,:),matrix_nphturnover(p,:,:)),vegmatrixn_old(:,:)) &
                                   + matmul(matmul(matrix_ngmtransfer(p,1:nvegnpool,:),matrix_ngmturnover(p,:,:)),vegmatrixn_old(:,:)) &
                                   + matmul(matmul(matrix_nfitransfer(p,1:nvegnpool,:),matrix_nfiturnover(p,:,:)),vegmatrixn_old(:,:))) * dt
         vegmatrixn_new(:) = matmul(vegmatrixn_old(:,:),rowonen(:)) +  vegmatrixn_input(:) + matmul(vegmatrixn_transfer(:,:),rowonen(:))

!         if(p .eq. 16)then
!            print*,'after matrix_new',p,vegmatrixn_new(p,iretransn)
!            print*,'after matrix_new',p,vegmatrixn_new(p,ifroot)
!         end if
 
! 

         if(isspinup .or. is_outmatrix)then
            matrix_calloc_leaf_acc(p)        = matrix_calloc_leaf_acc(p)        + vegmatrixc_input(ileaf)
            matrix_calloc_leafst_acc(p)      = matrix_calloc_leafst_acc(p)      + vegmatrixc_input(ileaf_st)
            matrix_calloc_froot_acc(p)       = matrix_calloc_froot_acc(p)       + vegmatrixc_input(ifroot)
            matrix_calloc_frootst_acc(p)     = matrix_calloc_frootst_acc(p)     + vegmatrixc_input(ifroot_st)
            matrix_calloc_livestem_acc(p)    = matrix_calloc_livestem_acc(p)    + vegmatrixc_input(ilivestem)
            matrix_calloc_livestemst_acc(p)  = matrix_calloc_livestemst_acc(p)  + vegmatrixc_input(ilivestem_st)
            matrix_calloc_deadstem_acc(p)    = matrix_calloc_deadstem_acc(p)    + vegmatrixc_input(ideadstem)
            matrix_calloc_deadstemst_acc(p)  = matrix_calloc_deadstemst_acc(p)  + vegmatrixc_input(ideadstem_st)
            matrix_calloc_livecroot_acc(p)   = matrix_calloc_livecroot_acc(p)   + vegmatrixc_input(ilivecroot)
            matrix_calloc_livecrootst_acc(p) = matrix_calloc_livecrootst_acc(p) + vegmatrixc_input(ilivecroot_st)
            matrix_calloc_deadcroot_acc(p)   = matrix_calloc_deadcroot_acc(p)   + vegmatrixc_input(ideadcroot)
            matrix_calloc_deadcrootst_acc(p) = matrix_calloc_deadcrootst_acc(p) + vegmatrixc_input(ideadcroot_st)
            if(ivt(p) >= npcropmin)then
               matrix_calloc_grain_acc(p)    = matrix_calloc_grain_acc(p)       + vegmatrixc_input(igrain)
               matrix_calloc_grainst_acc(p)  = matrix_calloc_grainst_acc(p)     + vegmatrixc_input(igrain_st)
            end if

!            print*,'before accumulated deadcrootc_acc',matrix_calloc_deadcroot_acc(p),vegmatrixc_input(ideadcroot),matrix_cturnover_deadcroot_acc(p),vegmatrixc_transfer(ideadcroot,ideadcroot)

            matrix_ctransfer_leafst_to_leafxf_acc(p)           = matrix_ctransfer_leafst_to_leafxf_acc(p) &
                                                               + vegmatrixc_transfer(ileaf_xf,ileaf_st)
            matrix_ctransfer_leafxf_to_leaf_acc(p)             = matrix_ctransfer_leafxf_to_leaf_acc(p) &
                                                               + vegmatrixc_transfer(ileaf,ileaf_xf)
            matrix_ctransfer_frootst_to_frootxf_acc(p)         = matrix_ctransfer_frootst_to_frootxf_acc(p) &
                                                               + vegmatrixc_transfer(ifroot_xf,ifroot_st)
            matrix_ctransfer_frootxf_to_froot_acc(p)           = matrix_ctransfer_frootxf_to_froot_acc(p) &
                                                               + vegmatrixc_transfer(ifroot,ifroot_xf)
            matrix_ctransfer_livestemst_to_livestemxf_acc(p)   = matrix_ctransfer_livestemst_to_livestemxf_acc(p) &
                                                               + vegmatrixc_transfer(ilivestem_xf,ilivestem_st)
            matrix_ctransfer_livestemxf_to_livestem_acc(p)     = matrix_ctransfer_livestemxf_to_livestem_acc(p) &
                                                               + vegmatrixc_transfer(ilivestem,ilivestem_xf)
            matrix_ctransfer_deadstemst_to_deadstemxf_acc(p)   = matrix_ctransfer_deadstemst_to_deadstemxf_acc(p) &
                                                               + vegmatrixc_transfer(ideadstem_xf,ideadstem_st)
            matrix_ctransfer_deadstemxf_to_deadstem_acc(p)     = matrix_ctransfer_deadstemxf_to_deadstem_acc(p) &
                                                               + vegmatrixc_transfer(ideadstem,ideadstem_xf)
            matrix_ctransfer_livecrootst_to_livecrootxf_acc(p) = matrix_ctransfer_livecrootst_to_livecrootxf_acc(p) &
                                                               + vegmatrixc_transfer(ilivecroot_xf,ilivecroot_st)
            matrix_ctransfer_livecrootxf_to_livecroot_acc(p)   = matrix_ctransfer_livecrootxf_to_livecroot_acc(p) &
                                                               + vegmatrixc_transfer(ilivecroot,ilivecroot_xf)
            matrix_ctransfer_deadcrootst_to_deadcrootxf_acc(p) = matrix_ctransfer_deadcrootst_to_deadcrootxf_acc(p) &
                                                               + vegmatrixc_transfer(ideadcroot_xf,ideadcroot_st)
            matrix_ctransfer_deadcrootxf_to_deadcroot_acc(p)   = matrix_ctransfer_deadcrootxf_to_deadcroot_acc(p) &
                                                               + vegmatrixc_transfer(ideadcroot,ideadcroot_xf)
            matrix_ctransfer_livestem_to_deadstem_acc(p)       = matrix_ctransfer_livestem_to_deadstem_acc(p) &
                                                               + vegmatrixc_transfer(ideadstem,ilivestem)
            matrix_ctransfer_livecroot_to_deadcroot_acc(p)     = matrix_ctransfer_livecroot_to_deadcroot_acc(p) &
                                                               + vegmatrixc_transfer(ideadcroot,ilivecroot)
            if(ivt(p) >= npcropmin)then
               matrix_ctransfer_grainst_to_grainxf_acc(p)      = matrix_ctransfer_grainst_to_grainxf_acc(p) &
                                                               + vegmatrixc_transfer(igrain_xf,igrain_st)
               matrix_ctransfer_grainxf_to_grain_acc(p)        = matrix_ctransfer_grainxf_to_grain_acc(p) &
                                                               + vegmatrixc_transfer(igrain,igrain_xf)
            end if
            matrix_cturnover_leaf_acc(p)                       = matrix_cturnover_leaf_acc(p) &
                                                               + vegmatrixc_transfer(ileaf,ileaf)
            matrix_cturnover_leafst_acc(p)                     = matrix_cturnover_leafst_acc(p) &
                                                               + vegmatrixc_transfer(ileaf_st,ileaf_st)
            matrix_cturnover_leafxf_acc(p)                     = matrix_cturnover_leafxf_acc(p) &
                                                               + vegmatrixc_transfer(ileaf_xf,ileaf_xf)
            matrix_cturnover_froot_acc(p)                      = matrix_cturnover_froot_acc(p) &
                                                               + vegmatrixc_transfer(ifroot,ifroot)
            matrix_cturnover_frootst_acc(p)                    = matrix_cturnover_frootst_acc(p) &
                                                               + vegmatrixc_transfer(ifroot_st,ifroot_st)
            matrix_cturnover_frootxf_acc(p)                    = matrix_cturnover_frootxf_acc(p) &
                                                               + vegmatrixc_transfer(ifroot_xf,ifroot_xf)
            matrix_cturnover_livestem_acc(p)                   = matrix_cturnover_livestem_acc(p) &
                                                               + vegmatrixc_transfer(ilivestem,ilivestem)
            matrix_cturnover_livestemst_acc(p)                 = matrix_cturnover_livestemst_acc(p) &
                                                               + vegmatrixc_transfer(ilivestem_st,ilivestem_st)
            matrix_cturnover_livestemxf_acc(p)                 = matrix_cturnover_livestemxf_acc(p) &
                                                               + vegmatrixc_transfer(ilivestem_xf,ilivestem_xf)
            matrix_cturnover_deadstem_acc(p)                   = matrix_cturnover_deadstem_acc(p) &
                                                               + vegmatrixc_transfer(ideadstem,ideadstem)
            matrix_cturnover_deadstemst_acc(p)                 = matrix_cturnover_deadstemst_acc(p) &
                                                               + vegmatrixc_transfer(ideadstem_st,ideadstem_st)
            matrix_cturnover_deadstemxf_acc(p)                 = matrix_cturnover_deadstemxf_acc(p) &
                                                               + vegmatrixc_transfer(ideadstem_xf,ideadstem_xf)
            matrix_cturnover_livecroot_acc(p)                  = matrix_cturnover_livecroot_acc(p) &
                                                               + vegmatrixc_transfer(ilivecroot,ilivecroot)
            matrix_cturnover_livecrootst_acc(p)                = matrix_cturnover_livecrootst_acc(p) &
                                                               + vegmatrixc_transfer(ilivecroot_st,ilivecroot_st)
            matrix_cturnover_livecrootxf_acc(p)                = matrix_cturnover_livecrootxf_acc(p) &
                                                               + vegmatrixc_transfer(ilivecroot_xf,ilivecroot_xf)
            matrix_cturnover_deadcroot_acc(p)                  = matrix_cturnover_deadcroot_acc(p) &
                                                               + vegmatrixc_transfer(ideadcroot,ideadcroot)
            matrix_cturnover_deadcrootst_acc(p)                = matrix_cturnover_deadcrootst_acc(p) &
                                                               + vegmatrixc_transfer(ideadcroot_st,ideadcroot_st)
            matrix_cturnover_deadcrootxf_acc(p)                = matrix_cturnover_deadcrootxf_acc(p) &
                                                               + vegmatrixc_transfer(ideadcroot_xf,ideadcroot_xf)
            if(ivt(p) >= npcropmin)then
               matrix_cturnover_grain_acc(p)                   = matrix_cturnover_grain_acc(p) &
                                                               + vegmatrixn_transfer(igrain,igrain)
               matrix_cturnover_grainst_acc(p)                 = matrix_cturnover_grainst_acc(p) &
                                                               + vegmatrixn_transfer(igrain_st,igrain_st)
               matrix_cturnover_grainxf_acc(p)                 = matrix_cturnover_grainxf_acc(p) &
                                                               + vegmatrixn_transfer(igrain_xf,igrain_xf)
            end if

            matrix_nalloc_leaf_acc(p)        = matrix_nalloc_leaf_acc(p)        + vegmatrixn_input(ileaf)
            matrix_nalloc_leafst_acc(p)      = matrix_nalloc_leafst_acc(p)      + vegmatrixn_input(ileaf_st)
            matrix_nalloc_froot_acc(p)       = matrix_nalloc_froot_acc(p)       + vegmatrixn_input(ifroot)
            matrix_nalloc_frootst_acc(p)     = matrix_nalloc_frootst_acc(p)     + vegmatrixn_input(ifroot_st)
            matrix_nalloc_livestem_acc(p)    = matrix_nalloc_livestem_acc(p)    + vegmatrixn_input(ilivestem)
            matrix_nalloc_livestemst_acc(p)  = matrix_nalloc_livestemst_acc(p)  + vegmatrixn_input(ilivestem_st)
            matrix_nalloc_deadstem_acc(p)    = matrix_nalloc_deadstem_acc(p)    + vegmatrixn_input(ideadstem)
            matrix_nalloc_deadstemst_acc(p)  = matrix_nalloc_deadstemst_acc(p)  + vegmatrixn_input(ideadstem_st)
            matrix_nalloc_livecroot_acc(p)   = matrix_nalloc_livecroot_acc(p)   + vegmatrixn_input(ilivecroot)
            matrix_nalloc_livecrootst_acc(p) = matrix_nalloc_livecrootst_acc(p) + vegmatrixn_input(ilivecroot_st)
            matrix_nalloc_deadcroot_acc(p)   = matrix_nalloc_deadcroot_acc(p)   + vegmatrixn_input(ideadcroot)
            matrix_nalloc_deadcrootst_acc(p) = matrix_nalloc_deadcrootst_acc(p) + vegmatrixn_input(ideadcroot_st)
            if(ivt(p) >= npcropmin)then
               matrix_nalloc_grain_acc(p)    = matrix_nalloc_grain_acc(p)       + vegmatrixn_input(igrain)
               matrix_nalloc_grainst_acc(p)  = matrix_nalloc_grainst_acc(p)     + vegmatrixn_input(igrain_st)
            end if

            matrix_ntransfer_leafst_to_leafxf_acc(p)           = matrix_ntransfer_leafst_to_leafxf_acc(p) &
                                                               + vegmatrixn_transfer(ileaf_xf,ileaf_st)
            matrix_ntransfer_leafxf_to_leaf_acc(p)             = matrix_ntransfer_leafxf_to_leaf_acc(p) &
                                                               + vegmatrixn_transfer(ileaf,ileaf_xf)
            matrix_ntransfer_frootst_to_frootxf_acc(p)         = matrix_ntransfer_frootst_to_frootxf_acc(p) &
                                                               + vegmatrixn_transfer(ifroot_xf,ifroot_st)
            matrix_ntransfer_frootxf_to_froot_acc(p)           = matrix_ntransfer_frootxf_to_froot_acc(p) &
                                                               + vegmatrixn_transfer(ifroot,ifroot_xf)
            matrix_ntransfer_livestemst_to_livestemxf_acc(p)   = matrix_ntransfer_livestemst_to_livestemxf_acc(p) &
                                                               + vegmatrixn_transfer(ilivestem_xf,ilivestem_st)
            matrix_ntransfer_livestemxf_to_livestem_acc(p)     = matrix_ntransfer_livestemxf_to_livestem_acc(p) &
                                                               + vegmatrixn_transfer(ilivestem,ilivestem_xf)
            matrix_ntransfer_deadstemst_to_deadstemxf_acc(p)   = matrix_ntransfer_deadstemst_to_deadstemxf_acc(p) &
                                                               + vegmatrixn_transfer(ideadstem_xf,ideadstem_st)
            matrix_ntransfer_deadstemxf_to_deadstem_acc(p)     = matrix_ntransfer_deadstemxf_to_deadstem_acc(p) &
                                                               + vegmatrixn_transfer(ideadstem,ideadstem_xf)
            matrix_ntransfer_livecrootst_to_livecrootxf_acc(p) = matrix_ntransfer_livecrootst_to_livecrootxf_acc(p) &
                                                               + vegmatrixn_transfer(ilivecroot_xf,ilivecroot_st)
            matrix_ntransfer_livecrootxf_to_livecroot_acc(p)   = matrix_ntransfer_livecrootxf_to_livecroot_acc(p) &
                                                               + vegmatrixn_transfer(ilivecroot,ilivecroot_xf)
            matrix_ntransfer_deadcrootst_to_deadcrootxf_acc(p) = matrix_ntransfer_deadcrootst_to_deadcrootxf_acc(p) &
                                                               + vegmatrixn_transfer(ideadcroot_xf,ideadcroot_st)
            matrix_ntransfer_deadcrootxf_to_deadcroot_acc(p)   = matrix_ntransfer_deadcrootxf_to_deadcroot_acc(p) &
                                                               + vegmatrixn_transfer(ideadcroot,ideadcroot_xf)
            matrix_ntransfer_livestem_to_deadstem_acc(p)       = matrix_ntransfer_livestem_to_deadstem_acc(p) &
                                                               + vegmatrixn_transfer(ideadstem,ilivestem)
            matrix_ntransfer_livecroot_to_deadcroot_acc(p)     = matrix_ntransfer_livecroot_to_deadcroot_acc(p) &
                                                               + vegmatrixn_transfer(ideadcroot,ilivecroot)
            if(ivt(p) >= npcropmin)then
               matrix_ntransfer_grainst_to_grainxf_acc(p)      = matrix_ntransfer_grainst_to_grainxf_acc(p) &
                                                               + vegmatrixn_transfer(igrain_xf,igrain_st)
               matrix_ntransfer_grainxf_to_grain_acc(p)        = matrix_ntransfer_grainxf_to_grain_acc(p) &
                                                               + vegmatrixn_transfer(igrain,igrain_xf)
            end if

            matrix_ntransfer_retransn_to_leaf_acc(p)           = matrix_ntransfer_retransn_to_leaf_acc(p) &
                                                               + vegmatrixn_transfer(ileaf,iretransn)
            matrix_ntransfer_retransn_to_leafst_acc(p)         = matrix_ntransfer_retransn_to_leafst_acc(p) &
                                                               + vegmatrixn_transfer(ileaf_st,iretransn)
            matrix_ntransfer_retransn_to_froot_acc(p)          = matrix_ntransfer_retransn_to_froot_acc(p) &
                                                               + vegmatrixn_transfer(ifroot,iretransn)
            matrix_ntransfer_retransn_to_frootst_acc(p)        = matrix_ntransfer_retransn_to_frootst_acc(p) &
                                                               + vegmatrixn_transfer(ifroot_st,iretransn)
            matrix_ntransfer_retransn_to_livestem_acc(p)       = matrix_ntransfer_retransn_to_livestem_acc(p) &
                                                               + vegmatrixn_transfer(ilivestem,iretransn)
            matrix_ntransfer_retransn_to_livestemst_acc(p)     = matrix_ntransfer_retransn_to_livestemst_acc(p) &
                                                               + vegmatrixn_transfer(ilivestem_st,iretransn)
            matrix_ntransfer_retransn_to_deadstem_acc(p)       = matrix_ntransfer_retransn_to_deadstem_acc(p) &
                                                               + vegmatrixn_transfer(ideadstem,iretransn)
            matrix_ntransfer_retransn_to_deadstemst_acc(p)     = matrix_ntransfer_retransn_to_deadstemst_acc(p) &
                                                               + vegmatrixn_transfer(ideadstem_st,iretransn)
            matrix_ntransfer_retransn_to_livecroot_acc(p)      = matrix_ntransfer_retransn_to_livecroot_acc(p) &
                                                               + vegmatrixn_transfer(ilivecroot,iretransn)
            matrix_ntransfer_retransn_to_livecrootst_acc(p)    = matrix_ntransfer_retransn_to_livecrootst_acc(p) &
                                                               + vegmatrixn_transfer(ilivecroot_st,iretransn)
            matrix_ntransfer_retransn_to_deadcroot_acc(p)      = matrix_ntransfer_retransn_to_deadcroot_acc(p) &
                                                               + vegmatrixn_transfer(ideadcroot,iretransn)
            matrix_ntransfer_retransn_to_deadcrootst_acc(p)    = matrix_ntransfer_retransn_to_deadcrootst_acc(p) &
                                                               + vegmatrixn_transfer(ideadcroot_st,iretransn)
            matrix_ntransfer_leaf_to_retransn_acc(p)           = matrix_ntransfer_leaf_to_retransn_acc(p) &
                                                               + vegmatrixn_transfer(iretransn,ileaf)
            matrix_ntransfer_froot_to_retransn_acc(p)          = matrix_ntransfer_froot_to_retransn_acc(p) &
                                                               + vegmatrixn_transfer(iretransn,ifroot)
            matrix_ntransfer_livestem_to_retransn_acc(p)       = matrix_ntransfer_livestem_to_retransn_acc(p) &
                                                               + vegmatrixn_transfer(iretransn,ilivestem)
            matrix_ntransfer_livecroot_to_retransn_acc(p)      = matrix_ntransfer_livecroot_to_retransn_acc(p) &
                                                               + vegmatrixn_transfer(iretransn,ilivecroot)
            if(ivt(p) >= npcropmin)then
               matrix_ntransfer_retransn_to_grain_acc(p)       = matrix_ntransfer_retransn_to_grain_acc(p) &
                                                               + vegmatrixn_transfer(igrain,iretransn)
               matrix_ntransfer_retransn_to_grainst_acc(p)     = matrix_ntransfer_retransn_to_grainst_acc(p) &
                                                               + vegmatrixn_transfer(igrain_st,iretransn)
            end if

            matrix_nturnover_leaf_acc(p)                       = matrix_nturnover_leaf_acc(p) &
                                                               + vegmatrixn_transfer(ileaf,ileaf)
            matrix_nturnover_leafst_acc(p)                     = matrix_nturnover_leafst_acc(p) &
                                                               + vegmatrixn_transfer(ileaf_st,ileaf_st)
            matrix_nturnover_leafxf_acc(p)                     = matrix_nturnover_leafxf_acc(p) &
                                                               + vegmatrixn_transfer(ileaf_xf,ileaf_xf)
            matrix_nturnover_froot_acc(p)                      = matrix_nturnover_froot_acc(p) &
                                                               + vegmatrixn_transfer(ifroot,ifroot)
            matrix_nturnover_frootst_acc(p)                    = matrix_nturnover_frootst_acc(p) &
                                                               + vegmatrixn_transfer(ifroot_st,ifroot_st)
            matrix_nturnover_frootxf_acc(p)                    = matrix_nturnover_frootxf_acc(p) &
                                                               + vegmatrixn_transfer(ifroot_xf,ifroot_xf)
            matrix_nturnover_livestem_acc(p)                   = matrix_nturnover_livestem_acc(p) &
                                                               + vegmatrixn_transfer(ilivestem,ilivestem)
            matrix_nturnover_livestemst_acc(p)                 = matrix_nturnover_livestemst_acc(p) &
                                                               + vegmatrixn_transfer(ilivestem_st,ilivestem_st)
            matrix_nturnover_livestemxf_acc(p)                 = matrix_nturnover_livestemxf_acc(p) &
                                                               + vegmatrixn_transfer(ilivestem_xf,ilivestem_xf)
            matrix_nturnover_deadstem_acc(p)                   = matrix_nturnover_deadstem_acc(p) &
                                                               + vegmatrixn_transfer(ideadstem,ideadstem)
            matrix_nturnover_deadstemst_acc(p)                 = matrix_nturnover_deadstemst_acc(p) &
                                                               + vegmatrixn_transfer(ideadstem_st,ideadstem_st)
            matrix_nturnover_deadstemxf_acc(p)                 = matrix_nturnover_deadstemxf_acc(p) &
                                                               + vegmatrixn_transfer(ideadstem_xf,ideadstem_xf)
            matrix_nturnover_livecroot_acc(p)                  = matrix_nturnover_livecroot_acc(p) &
                                                               + vegmatrixn_transfer(ilivecroot,ilivecroot)
            matrix_nturnover_livecrootst_acc(p)                = matrix_nturnover_livecrootst_acc(p) &
                                                               + vegmatrixn_transfer(ilivecroot_st,ilivecroot_st)
            matrix_nturnover_livecrootxf_acc(p)                = matrix_nturnover_livecrootxf_acc(p) &
                                                               + vegmatrixn_transfer(ilivecroot_xf,ilivecroot_xf)
            matrix_nturnover_deadcroot_acc(p)                  = matrix_nturnover_deadcroot_acc(p) &
                                                               + vegmatrixn_transfer(ideadcroot,ideadcroot)
            matrix_nturnover_deadcrootst_acc(p)                = matrix_nturnover_deadcrootst_acc(p) &
                                                               + vegmatrixn_transfer(ideadcroot_st,ideadcroot_st)
            matrix_nturnover_deadcrootxf_acc(p)                = matrix_nturnover_deadcrootxf_acc(p) &
                                                               + vegmatrixn_transfer(ideadcroot_xf,ideadcroot_xf)
            if(ivt(p) >= npcropmin)then
               matrix_nturnover_grain_acc(p)                   = matrix_nturnover_grain_acc(p) &
                                                               + vegmatrixn_transfer(igrain,igrain)
               matrix_nturnover_grainst_acc(p)                 = matrix_nturnover_grainst_acc(p) &
                                                               + vegmatrixn_transfer(igrain_st,igrain_st)
               matrix_nturnover_grainxf_acc(p)                 = matrix_nturnover_grainxf_acc(p) &
                                                               + vegmatrixn_transfer(igrain_xf,igrain_xf)
            end if
            matrix_nturnover_retransn_acc(p)                   = matrix_nturnover_retransn_acc(p) &
                                                               + vegmatrixn_transfer(iretransn,iretransn)
         end if 

         leafc(p)               = vegmatrixc_new(ileaf)
         leafc_storage(p)       = vegmatrixc_new(ileaf_st)
         leafc_xfer(p)          = vegmatrixc_new(ileaf_xf)
         frootc(p)              = vegmatrixc_new(ifroot)
         frootc_storage(p)      = vegmatrixc_new(ifroot_st)
         frootc_xfer(p)         = vegmatrixc_new(ifroot_xf)
         livestemc(p)           = vegmatrixc_new(ilivestem)
         livestemc_storage(p)   = vegmatrixc_new(ilivestem_st)
         livestemc_xfer(p)      = vegmatrixc_new(ilivestem_xf)
         deadstemc(p)           = vegmatrixc_new(ideadstem)
         deadstemc_storage(p)   = vegmatrixc_new(ideadstem_st)
         deadstemc_xfer(p)      = vegmatrixc_new(ideadstem_xf)
         livecrootc(p)          = vegmatrixc_new(ilivecroot)
         livecrootc_storage(p)  = vegmatrixc_new(ilivecroot_st)
         livecrootc_xfer(p)     = vegmatrixc_new(ilivecroot_xf)
         deadcrootc(p)          = vegmatrixc_new(ideadcroot)
         deadcrootc_storage(p)  = vegmatrixc_new(ideadcroot_st)
         deadcrootc_xfer(p)     = vegmatrixc_new(ideadcroot_xf)
         if(ivt(p) >= npcropmin)then
         grainc(p)              = vegmatrixc_new(igrain)
         grainc_storage(p)      = vegmatrixc_new(igrain_st)
         grainc_xfer(p)         = vegmatrixc_new(igrain_xf)
         end if
 
         leafn(p)               = vegmatrixn_new(ileaf)
         leafn_storage(p)       = vegmatrixn_new(ileaf_st)
         leafn_xfer(p)          = vegmatrixn_new(ileaf_xf)
         frootn(p)              = vegmatrixn_new(ifroot)
         frootn_storage(p)      = vegmatrixn_new(ifroot_st)
         frootn_xfer(p)         = vegmatrixn_new(ifroot_xf)
         livestemn(p)           = vegmatrixn_new(ilivestem)
         livestemn_storage(p)   = vegmatrixn_new(ilivestem_st)
         livestemn_xfer(p)      = vegmatrixn_new(ilivestem_xf)
         deadstemn(p)           = vegmatrixn_new(ideadstem)
         deadstemn_storage(p)   = vegmatrixn_new(ideadstem_st)
         deadstemn_xfer(p)      = vegmatrixn_new(ideadstem_xf)
         livecrootn(p)          = vegmatrixn_new(ilivecroot)
         livecrootn_storage(p)  = vegmatrixn_new(ilivecroot_st)
         livecrootn_xfer(p)     = vegmatrixn_new(ilivecroot_xf)
         deadcrootn(p)          = vegmatrixn_new(ideadcroot)
         deadcrootn_storage(p)  = vegmatrixn_new(ideadcroot_st)
         deadcrootn_xfer(p)     = vegmatrixn_new(ideadcroot_xf)
         if(ivt(p) >= npcropmin)then
            grainn(p)              = vegmatrixn_new(igrain)
            grainn_storage(p)      = vegmatrixn_new(igrain_st)
            grainn_xfer(p)         = vegmatrixn_new(igrain_xf)
         end if
         retransn(p)            = vegmatrixn_new(iretransn)
         if(isspinup .or. is_outmatrix)then
!print *, 'count',counter, end_of_year
            if(end_of_year)then
               matrix_calloc_acc(ileaf)         = matrix_calloc_leaf_acc(p)               
               matrix_calloc_acc(ileaf_st)      = matrix_calloc_leafst_acc(p)               
               matrix_calloc_acc(ifroot)        = matrix_calloc_froot_acc(p)               
               matrix_calloc_acc(ifroot_st)     = matrix_calloc_frootst_acc(p)               
               matrix_calloc_acc(ilivestem)     = matrix_calloc_livestem_acc(p)               
               matrix_calloc_acc(ilivestem_st)  = matrix_calloc_livestemst_acc(p)               
               matrix_calloc_acc(ideadstem)     = matrix_calloc_deadstem_acc(p)               
               matrix_calloc_acc(ideadstem_st)  = matrix_calloc_deadstemst_acc(p)               
               matrix_calloc_acc(ilivecroot)    = matrix_calloc_livecroot_acc(p)               
               matrix_calloc_acc(ilivecroot_st) = matrix_calloc_livecrootst_acc(p)               
               matrix_calloc_acc(ideadcroot)    = matrix_calloc_deadcroot_acc(p)               
               matrix_calloc_acc(ideadcroot_st) = matrix_calloc_deadcrootst_acc(p)               
               if(ivt(p) >= npcropmin)then
                  matrix_calloc_acc(igrain)     = matrix_calloc_grain_acc(p)               
                  matrix_calloc_acc(igrain_st)  = matrix_calloc_grainst_acc(p)               
               end if

               matrix_ctransfer_acc(ileaf_xf,ileaf_st)           = matrix_ctransfer_leafst_to_leafxf_acc(p)
               matrix_ctransfer_acc(ileaf,ileaf_xf)              = matrix_ctransfer_leafxf_to_leaf_acc(p)
               matrix_ctransfer_acc(ifroot_xf,ifroot_st)         = matrix_ctransfer_frootst_to_frootxf_acc(p)
               matrix_ctransfer_acc(ifroot,ifroot_xf)            = matrix_ctransfer_frootxf_to_froot_acc(p)
               matrix_ctransfer_acc(ilivestem_xf,ilivestem_st)   = matrix_ctransfer_livestemst_to_livestemxf_acc(p)
               matrix_ctransfer_acc(ilivestem,ilivestem_xf)      = matrix_ctransfer_livestemxf_to_livestem_acc(p)
               matrix_ctransfer_acc(ideadstem_xf,ideadstem_st)   = matrix_ctransfer_deadstemst_to_deadstemxf_acc(p)
               matrix_ctransfer_acc(ideadstem,ideadstem_xf)      = matrix_ctransfer_deadstemxf_to_deadstem_acc(p)
               matrix_ctransfer_acc(ilivecroot_xf,ilivecroot_st) = matrix_ctransfer_livecrootst_to_livecrootxf_acc(p)
               matrix_ctransfer_acc(ilivecroot,ilivecroot_xf)    = matrix_ctransfer_livecrootxf_to_livecroot_acc(p)
               matrix_ctransfer_acc(ideadcroot_xf,ideadcroot_st) = matrix_ctransfer_deadcrootst_to_deadcrootxf_acc(p)
               matrix_ctransfer_acc(ideadcroot,ideadcroot_xf)    = matrix_ctransfer_deadcrootxf_to_deadcroot_acc(p)
               if(ivt(p) >= npcropmin)then
                  matrix_ctransfer_acc(igrain_xf,igrain_st)      = matrix_ctransfer_grainst_to_grainxf_acc(p)
                  matrix_ctransfer_acc(igrain,igrain_xf)         = matrix_ctransfer_grainxf_to_grain_acc(p)
               end if
               matrix_ctransfer_acc(ideadstem,ilivestem)         = matrix_ctransfer_livestem_to_deadstem_acc(p)
               matrix_ctransfer_acc(ideadcroot,ilivecroot)       = matrix_ctransfer_livecroot_to_deadcroot_acc(p)

               matrix_ctransfer_acc(ileaf,ileaf)                 = matrix_cturnover_leaf_acc(p)               
               matrix_ctransfer_acc(ileaf_st,ileaf_st)           = matrix_cturnover_leafst_acc(p)               
               matrix_ctransfer_acc(ileaf_xf,ileaf_xf)           = matrix_cturnover_leafxf_acc(p)               
               matrix_ctransfer_acc(ifroot,ifroot)               = matrix_cturnover_froot_acc(p)               
               matrix_ctransfer_acc(ifroot_st,ifroot_st)         = matrix_cturnover_frootst_acc(p)               
               matrix_ctransfer_acc(ifroot_xf,ifroot_xf)         = matrix_cturnover_frootxf_acc(p)               
               matrix_ctransfer_acc(ilivestem,ilivestem)         = matrix_cturnover_livestem_acc(p)               
               matrix_ctransfer_acc(ilivestem_st,ilivestem_st)   = matrix_cturnover_livestemst_acc(p)               
               matrix_ctransfer_acc(ilivestem_xf,ilivestem_xf)   = matrix_cturnover_livestemxf_acc(p)               
               matrix_ctransfer_acc(ideadstem,ideadstem)         = matrix_cturnover_deadstem_acc(p)               
               matrix_ctransfer_acc(ideadstem_st,ideadstem_st)   = matrix_cturnover_deadstemst_acc(p)               
               matrix_ctransfer_acc(ideadstem_xf,ideadstem_xf)   = matrix_cturnover_deadstemxf_acc(p)               
               matrix_ctransfer_acc(ilivecroot,ilivecroot)       = matrix_cturnover_livecroot_acc(p)               
               matrix_ctransfer_acc(ilivecroot_st,ilivecroot_st) = matrix_cturnover_livecrootst_acc(p)               
               matrix_ctransfer_acc(ilivecroot_xf,ilivecroot_xf) = matrix_cturnover_livecrootxf_acc(p)               
               matrix_ctransfer_acc(ideadcroot,ideadcroot)       = matrix_cturnover_deadcroot_acc(p)               
               matrix_ctransfer_acc(ideadcroot_st,ideadcroot_st) = matrix_cturnover_deadcrootst_acc(p)               
               matrix_ctransfer_acc(ideadcroot_xf,ideadcroot_xf) = matrix_cturnover_deadcrootxf_acc(p)               
               if(ivt(p) >= npcropmin)then
                  matrix_ctransfer_acc(igrain,igrain)            = matrix_cturnover_grain_acc(p)               
                  matrix_ctransfer_acc(igrain_st,igrain_st)      = matrix_cturnover_grainst_acc(p)               
                  matrix_ctransfer_acc(igrain_xf,igrain_xf)      = matrix_cturnover_grainxf_acc(p)               
               end if

               matrix_nalloc_acc(ileaf)         = matrix_nalloc_leaf_acc(p)               
               matrix_nalloc_acc(ileaf_st)      = matrix_nalloc_leafst_acc(p)               
               matrix_nalloc_acc(ifroot)        = matrix_nalloc_froot_acc(p)               
               matrix_nalloc_acc(ifroot_st)     = matrix_nalloc_frootst_acc(p)               
               matrix_nalloc_acc(ilivestem)     = matrix_nalloc_livestem_acc(p)               
               matrix_nalloc_acc(ilivestem_st)  = matrix_nalloc_livestemst_acc(p)               
               matrix_nalloc_acc(ideadstem)     = matrix_nalloc_deadstem_acc(p)               
               matrix_nalloc_acc(ideadstem_st)  = matrix_nalloc_deadstemst_acc(p)               
               matrix_nalloc_acc(ilivecroot)    = matrix_nalloc_livecroot_acc(p)               
               matrix_nalloc_acc(ilivecroot_st) = matrix_nalloc_livecrootst_acc(p)               
               matrix_nalloc_acc(ideadcroot)    = matrix_nalloc_deadcroot_acc(p)               
               matrix_nalloc_acc(ideadcroot_st) = matrix_nalloc_deadcrootst_acc(p)               
               if(ivt(p) >= npcropmin)then
                  matrix_nalloc_acc(igrain)     = matrix_nalloc_grain_acc(p)               
                  matrix_nalloc_acc(igrain_st)  = matrix_nalloc_grainst_acc(p)               
               end if

               matrix_ntransfer_acc(ileaf_xf,ileaf_st)           = matrix_ntransfer_leafst_to_leafxf_acc(p)
               matrix_ntransfer_acc(ileaf,ileaf_xf)              = matrix_ntransfer_leafxf_to_leaf_acc(p)
               matrix_ntransfer_acc(ifroot_xf,ifroot_st)         = matrix_ntransfer_frootst_to_frootxf_acc(p)
               matrix_ntransfer_acc(ifroot,ifroot_xf)            = matrix_ntransfer_frootxf_to_froot_acc(p)
               matrix_ntransfer_acc(ilivestem_xf,ilivestem_st)   = matrix_ntransfer_livestemst_to_livestemxf_acc(p)
               matrix_ntransfer_acc(ilivestem,ilivestem_xf)      = matrix_ntransfer_livestemxf_to_livestem_acc(p)
               matrix_ntransfer_acc(ideadstem_xf,ideadstem_st)   = matrix_ntransfer_deadstemst_to_deadstemxf_acc(p)
               matrix_ntransfer_acc(ideadstem,ideadstem_xf)      = matrix_ntransfer_deadstemxf_to_deadstem_acc(p)
               matrix_ntransfer_acc(ilivecroot_xf,ilivecroot_st) = matrix_ntransfer_livecrootst_to_livecrootxf_acc(p)
               matrix_ntransfer_acc(ilivecroot,ilivecroot_xf)    = matrix_ntransfer_livecrootxf_to_livecroot_acc(p)
               matrix_ntransfer_acc(ideadcroot_xf,ideadcroot_st) = matrix_ntransfer_deadcrootst_to_deadcrootxf_acc(p)
               matrix_ntransfer_acc(ideadcroot,ideadcroot_xf)    = matrix_ntransfer_deadcrootxf_to_deadcroot_acc(p)
               if(ivt(p) >= npcropmin)then
                  matrix_ntransfer_acc(igrain_xf,igrain_st)      = matrix_ntransfer_grainst_to_grainxf_acc(p)
                  matrix_ntransfer_acc(igrain,igrain_xf)         = matrix_ntransfer_grainxf_to_grain_acc(p)
               end if
               matrix_ntransfer_acc(ideadstem,ilivestem)         = matrix_ntransfer_livestem_to_deadstem_acc(p)
               matrix_ntransfer_acc(ideadcroot,ilivecroot)       = matrix_ntransfer_livecroot_to_deadcroot_acc(p)

               matrix_ntransfer_acc(ileaf,iretransn)         = matrix_ntransfer_retransn_to_leaf_acc(p)               
               matrix_ntransfer_acc(ileaf_st,iretransn)      = matrix_ntransfer_retransn_to_leafst_acc(p)               
               matrix_ntransfer_acc(ifroot,iretransn)        = matrix_ntransfer_retransn_to_froot_acc(p)               
               matrix_ntransfer_acc(ifroot_st,iretransn)     = matrix_ntransfer_retransn_to_frootst_acc(p)               
               matrix_ntransfer_acc(ilivestem,iretransn)     = matrix_ntransfer_retransn_to_livestem_acc(p)               
               matrix_ntransfer_acc(ilivestem_st,iretransn)  = matrix_ntransfer_retransn_to_livestemst_acc(p)               
               matrix_ntransfer_acc(ideadstem,iretransn)     = matrix_ntransfer_retransn_to_deadstem_acc(p)               
               matrix_ntransfer_acc(ideadstem_st,iretransn)  = matrix_ntransfer_retransn_to_deadstemst_acc(p)               
               matrix_ntransfer_acc(ilivecroot,iretransn)    = matrix_ntransfer_retransn_to_livecroot_acc(p)               
               matrix_ntransfer_acc(ilivecroot_st,iretransn) = matrix_ntransfer_retransn_to_livecrootst_acc(p)               
               matrix_ntransfer_acc(ideadcroot,iretransn)    = matrix_ntransfer_retransn_to_deadcroot_acc(p)               
               matrix_ntransfer_acc(ideadcroot_st,iretransn) = matrix_ntransfer_retransn_to_deadcrootst_acc(p)               
               if(ivt(p) >= npcropmin)then
                  matrix_ntransfer_acc(igrain,iretransn)     = matrix_ntransfer_retransn_to_grain_acc(p)               
                  matrix_ntransfer_acc(igrain_st,iretransn)  = matrix_ntransfer_retransn_to_grainst_acc(p)               
               end if
               matrix_ntransfer_acc(iretransn,ileaf)         = matrix_ntransfer_leaf_to_retransn_acc(p)               
               matrix_ntransfer_acc(iretransn,ifroot)        = matrix_ntransfer_froot_to_retransn_acc(p)               
               matrix_ntransfer_acc(iretransn,ilivestem)     = matrix_ntransfer_livestem_to_retransn_acc(p)               
               matrix_ntransfer_acc(iretransn,ilivecroot)    = matrix_ntransfer_livecroot_to_retransn_acc(p)               

               matrix_ntransfer_acc(ileaf,ileaf)                 = matrix_nturnover_leaf_acc(p)               
               matrix_ntransfer_acc(ileaf_st,ileaf_st)           = matrix_nturnover_leafst_acc(p)               
               matrix_ntransfer_acc(ileaf_xf,ileaf_xf)           = matrix_nturnover_leafxf_acc(p)               
               matrix_ntransfer_acc(ifroot,ifroot)               = matrix_nturnover_froot_acc(p)               
               matrix_ntransfer_acc(ifroot_st,ifroot_st)         = matrix_nturnover_frootst_acc(p)               
               matrix_ntransfer_acc(ifroot_xf,ifroot_xf)         = matrix_nturnover_frootxf_acc(p)               
               matrix_ntransfer_acc(ilivestem,ilivestem)         = matrix_nturnover_livestem_acc(p)               
               matrix_ntransfer_acc(ilivestem_st,ilivestem_st)   = matrix_nturnover_livestemst_acc(p)               
               matrix_ntransfer_acc(ilivestem_xf,ilivestem_xf)   = matrix_nturnover_livestemxf_acc(p)               
               matrix_ntransfer_acc(ideadstem,ideadstem)         = matrix_nturnover_deadstem_acc(p)               
               matrix_ntransfer_acc(ideadstem_st,ideadstem_st)   = matrix_nturnover_deadstemst_acc(p)               
               matrix_ntransfer_acc(ideadstem_xf,ideadstem_xf)   = matrix_nturnover_deadstemxf_acc(p)               
               matrix_ntransfer_acc(ilivecroot,ilivecroot)       = matrix_nturnover_livecroot_acc(p)               
               matrix_ntransfer_acc(ilivecroot_st,ilivecroot_st) = matrix_nturnover_livecrootst_acc(p)               
               matrix_ntransfer_acc(ilivecroot_xf,ilivecroot_xf) = matrix_nturnover_livecrootxf_acc(p)               
               matrix_ntransfer_acc(ideadcroot,ideadcroot)       = matrix_nturnover_deadcroot_acc(p)               
               matrix_ntransfer_acc(ideadcroot_st,ideadcroot_st) = matrix_nturnover_deadcrootst_acc(p)               
               matrix_ntransfer_acc(ideadcroot_xf,ideadcroot_xf) = matrix_nturnover_deadcrootxf_acc(p)               
               if(ivt(p) >= npcropmin)then
                  matrix_ntransfer_acc(igrain,igrain)            = matrix_nturnover_grain_acc(p)               
                  matrix_ntransfer_acc(igrain_st,igrain_st)      = matrix_nturnover_grainst_acc(p)               
                  matrix_ntransfer_acc(igrain_xf,igrain_xf)      = matrix_nturnover_grainxf_acc(p)               
               end if
               matrix_ntransfer_acc(iretransn,iretransn)         = matrix_nturnover_retransn_acc(p)               

               do i=1,nvegcpool
                 if(matrix_ctransfer_acc(i,i) .eq. 0)then
                    matrix_ctransfer_acc(i,i) = 1.e+36
                 end if
               end do
               do i=1,nvegnpool
                 if(matrix_ntransfer_acc(i,i) .eq. 0)then
                    matrix_ntransfer_acc(i,i) = 1.e+36
                 end if
               end do
!              end if
               !print*,'before divided by retransn0,matrix_ntransfer_acc,matrix_nalloc_acc',deadcrootc0(p),matrix_ctransfer_acc(1:nvegnpool,ideadcroot),matrix_calloc_acc(ideadcroot)
               matrix_ctransfer_acc(1:nvegcpool,ileaf)         = matrix_ctransfer_acc(1:nvegcpool,ileaf)         / leafc0(p)
               matrix_ctransfer_acc(1:nvegcpool,ileaf_st)      = matrix_ctransfer_acc(1:nvegcpool,ileaf_st)      / leafc0_storage(p)
               matrix_ctransfer_acc(1:nvegcpool,ileaf_xf)      = matrix_ctransfer_acc(1:nvegcpool,ileaf_xf)      / leafc0_xfer(p)
               matrix_ctransfer_acc(1:nvegcpool,ifroot)        = matrix_ctransfer_acc(1:nvegcpool,ifroot)        / frootc0(p)
               matrix_ctransfer_acc(1:nvegcpool,ifroot_st)     = matrix_ctransfer_acc(1:nvegcpool,ifroot_st)     / frootc0_storage(p)
               matrix_ctransfer_acc(1:nvegcpool,ifroot_xf)     = matrix_ctransfer_acc(1:nvegcpool,ifroot_xf)     / frootc0_xfer(p)
               matrix_ctransfer_acc(1:nvegcpool,ilivestem)     = matrix_ctransfer_acc(1:nvegcpool,ilivestem)     / livestemc0(p)
               matrix_ctransfer_acc(1:nvegcpool,ilivestem_st)  = matrix_ctransfer_acc(1:nvegcpool,ilivestem_st)  / livestemc0_storage(p)
               matrix_ctransfer_acc(1:nvegcpool,ilivestem_xf)  = matrix_ctransfer_acc(1:nvegcpool,ilivestem_xf)  / livestemc0_xfer(p)
               matrix_ctransfer_acc(1:nvegcpool,ideadstem)     = matrix_ctransfer_acc(1:nvegcpool,ideadstem)     / deadstemc0(p)
               matrix_ctransfer_acc(1:nvegcpool,ideadstem_st)  = matrix_ctransfer_acc(1:nvegcpool,ideadstem_st)  / deadstemc0_storage(p)
               matrix_ctransfer_acc(1:nvegcpool,ideadstem_xf)  = matrix_ctransfer_acc(1:nvegcpool,ideadstem_xf)  / deadstemc0_xfer(p)
               matrix_ctransfer_acc(1:nvegcpool,ilivecroot)    = matrix_ctransfer_acc(1:nvegcpool,ilivecroot)    / livecrootc0(p)
               matrix_ctransfer_acc(1:nvegcpool,ilivecroot_st) = matrix_ctransfer_acc(1:nvegcpool,ilivecroot_st) / livecrootc0_storage(p)
               matrix_ctransfer_acc(1:nvegcpool,ilivecroot_xf) = matrix_ctransfer_acc(1:nvegcpool,ilivecroot_xf) / livecrootc0_xfer(p)
               matrix_ctransfer_acc(1:nvegcpool,ideadcroot)    = matrix_ctransfer_acc(1:nvegcpool,ideadcroot)    / deadcrootc0(p)
               matrix_ctransfer_acc(1:nvegcpool,ideadcroot_st) = matrix_ctransfer_acc(1:nvegcpool,ideadcroot_st) / deadcrootc0_storage(p)
               matrix_ctransfer_acc(1:nvegcpool,ideadcroot_xf) = matrix_ctransfer_acc(1:nvegcpool,ideadcroot_xf) / deadcrootc0_xfer(p)
               if(ivt(p) >= npcropmin)then
                  matrix_ctransfer_acc(1:nvegcpool,igrain)     = matrix_ctransfer_acc(1:nvegcpool,igrain)    / grainc0(p)
                  matrix_ctransfer_acc(1:nvegcpool,igrain_st)  = matrix_ctransfer_acc(1:nvegcpool,igrain_st) / grainc0_storage(p)
                  matrix_ctransfer_acc(1:nvegcpool,igrain_xf)  = matrix_ctransfer_acc(1:nvegcpool,igrain_xf) / grainc0_xfer(p)
               end if

               matrix_ntransfer_acc(1:nvegnpool,ileaf)         = matrix_ntransfer_acc(1:nvegnpool,ileaf)         / leafn0(p)
               matrix_ntransfer_acc(1:nvegnpool,ileaf_st)      = matrix_ntransfer_acc(1:nvegnpool,ileaf_st)      / leafn0_storage(p)
               matrix_ntransfer_acc(1:nvegnpool,ileaf_xf)      = matrix_ntransfer_acc(1:nvegnpool,ileaf_xf)      / leafn0_xfer(p)
               matrix_ntransfer_acc(1:nvegnpool,ifroot)        = matrix_ntransfer_acc(1:nvegnpool,ifroot)        / frootn0(p)
               matrix_ntransfer_acc(1:nvegnpool,ifroot_st)     = matrix_ntransfer_acc(1:nvegnpool,ifroot_st)     / frootn0_storage(p)
               matrix_ntransfer_acc(1:nvegnpool,ifroot_xf)     = matrix_ntransfer_acc(1:nvegnpool,ifroot_xf)     / frootn0_xfer(p)
               matrix_ntransfer_acc(1:nvegnpool,ilivestem)     = matrix_ntransfer_acc(1:nvegnpool,ilivestem)     / livestemn0(p)
               matrix_ntransfer_acc(1:nvegnpool,ilivestem_st)  = matrix_ntransfer_acc(1:nvegnpool,ilivestem_st)  / livestemn0_storage(p)
               matrix_ntransfer_acc(1:nvegnpool,ilivestem_xf)  = matrix_ntransfer_acc(1:nvegnpool,ilivestem_xf)  / livestemn0_xfer(p)
               matrix_ntransfer_acc(1:nvegnpool,ideadstem)     = matrix_ntransfer_acc(1:nvegnpool,ideadstem)     / deadstemn0(p)
               matrix_ntransfer_acc(1:nvegnpool,ideadstem_st)  = matrix_ntransfer_acc(1:nvegnpool,ideadstem_st)  / deadstemn0_storage(p)
               matrix_ntransfer_acc(1:nvegnpool,ideadstem_xf)  = matrix_ntransfer_acc(1:nvegnpool,ideadstem_xf)  / deadstemn0_xfer(p)
               matrix_ntransfer_acc(1:nvegnpool,ilivecroot)    = matrix_ntransfer_acc(1:nvegnpool,ilivecroot)    / livecrootn0(p)
               matrix_ntransfer_acc(1:nvegnpool,ilivecroot_st) = matrix_ntransfer_acc(1:nvegnpool,ilivecroot_st) / livecrootn0_storage(p)
               matrix_ntransfer_acc(1:nvegnpool,ilivecroot_xf) = matrix_ntransfer_acc(1:nvegnpool,ilivecroot_xf) / livecrootn0_xfer(p)
               matrix_ntransfer_acc(1:nvegnpool,ideadcroot)    = matrix_ntransfer_acc(1:nvegnpool,ideadcroot)    / deadcrootn0(p)
               matrix_ntransfer_acc(1:nvegnpool,ideadcroot_st) = matrix_ntransfer_acc(1:nvegnpool,ideadcroot_st) / deadcrootn0_storage(p)
               matrix_ntransfer_acc(1:nvegnpool,ideadcroot_xf) = matrix_ntransfer_acc(1:nvegnpool,ideadcroot_xf) / deadcrootn0_xfer(p)
               if(ivt(p) >= npcropmin)then
                  matrix_ntransfer_acc(1:nvegnpool,igrain)     = matrix_ntransfer_acc(1:nvegnpool,igrain)    / grainn0(p)
                  matrix_ntransfer_acc(1:nvegnpool,igrain_st)  = matrix_ntransfer_acc(1:nvegnpool,igrain_st) / grainn0_storage(p)
                  matrix_ntransfer_acc(1:nvegnpool,igrain_xf)  = matrix_ntransfer_acc(1:nvegnpool,igrain_xf) / grainn0_xfer(p)
               end if
               matrix_ntransfer_acc(1:nvegnpool,iretransn)  = matrix_ntransfer_acc(1:nvegnpool,iretransn) / retransn0(p)


!               print*,'matrix_ctransfer_acc,matrix_nalloc_acc',matrix_ctransfer_acc,matrix_calloc_acc(:)
               call inverse(matrix_ctransfer_acc(1:nvegcpool,1:nvegcpool),AKinvc(1:nvegcpool,1:nvegcpool),nvegcpool)
!               print*,'AKinvc',AKinvc(ideadcroot,1:nvegnpool),matrix_calloc_acc(1:nvegnpool)
               vegmatrixc_rt(:) = -matmul(AKinvc(1:nvegcpool,1:nvegcpool),matrix_calloc_acc(1:nvegcpool))
!               print*,'after updating vegmatrixn_rt',nvegcpool,vegmatrixc_rt(ideadcroot)
! N  
               call inverse(matrix_ntransfer_acc(1:nvegnpool,1:nvegnpool),AKinvn(1:nvegnpool,1:nvegnpool),nvegnpool)
               vegmatrixn_rt(:) = -matmul(AKinvn(1:nvegnpool,1:nvegnpool),matrix_nalloc_acc(1:nvegnpool))
!            do i=1,nvegpool
!               print*,'AKinv',i, AKinv(i,i), AKinvn(i,i)
!            end do
!          if(p .eq.7) print*,'input', matrix_alloc(p,ideadstem)* matrix_Cinput(p) * dt, livestemc(p)*((matrix_phtransfer(p,ideadstem,ilivestem)*matrix_phturnover(p,ilivestem,ilivestem)) &
!                                                                                             + (matrix_fitransfer(p,ideadstem,ilivestem)*matrix_fiturnover(p,ilivestem,ilivestem)))   
 
!          if(p .eq.7)print*,'output',deadstemc(p)*((matrix_phtransfer(p,iout,ideadstem)*matrix_phturnover(p,ideadstem,ideadstem))+ (matrix_fitransfer(p,iout,ideadstem)*matrix_fiturnover(p,ideadstem,ideadstem)))
               if(isspinup)then
               leafc(p)                  = vegmatrixc_rt(ileaf)
!               leafc_storage(p)          = vegmatrixc_rt(ileaf_st)
!               leafc_xfer(p)             = vegmatrixc_rt(ileaf_xf)
               frootc(p)                 = vegmatrixc_rt(ifroot)
!               frootc_storage(p)         = vegmatrixc_rt(ifroot_st)
!               frootc_xfer(p)            = vegmatrixc_rt(ifroot_xf)
               livestemc(p)              = vegmatrixc_rt(ilivestem)
!               livestemc_storage(p)      = vegmatrixc_rt(ilivestem_st)
!               livestemc_xfer(p)         = vegmatrixc_rt(ilivestem_xf)
               deadstemc(p)              = vegmatrixc_rt(ideadstem)
!               deadstemc_storage(p)      = vegmatrixc_rt(ideadstem_st)
!               deadstemc_xfer(p)         = vegmatrixc_rt(ideadstem_xf)
               livecrootc(p)              = vegmatrixc_rt(ilivecroot)
!               livecrootc_storage(p)      = vegmatrixc_rt(ilivecroot_st)
!               livecrootc_xfer(p)         = vegmatrixc_rt(ilivecroot_xf)   
               deadcrootc(p)              = vegmatrixc_rt(ideadcroot)
!               deadcrootc_storage(p)      = vegmatrixc_rt(ideadcroot_st)
!               deadcrootc_xfer(p)         = vegmatrixc_rt(ideadcroot_xf) 
               if(ivt(p) >= npcropmin)then
               grainc(p)                 = vegmatrixc_rt(igrain)
               end if
!               print *, 'leafacc',p,vegmatrixc_rt(:),vegmatrixn_rt(:)
               leafn(p)                  = vegmatrixn_rt(ileaf)
!               leafn_storage(p)          = vegmatrixn_rt(ileaf_st)
!               leafn_xfer(p)             = vegmatrixn_rt(ileaf_xf)
               frootn(p)                 = vegmatrixn_rt(ifroot)
!               frootn_storage(p)         = vegmatrixn_rt(ifroot_st)
!               frootn_xfer(p)            = vegmatrixn_rt(ifroot_xf)
               livestemn(p)              = vegmatrixn_rt(ilivestem)
!               livestemn_storage(p)      = vegmatrixn_rt(ilivestem_st)
!               livestemn_xfer(p)         = vegmatrixn_rt(ilivestem_xf)
               deadstemn(p)              = vegmatrixn_rt(ideadstem)
!               deadstemn_storage(p)      = vegmatrixn_rt(ideadstem_st)
!               deadstemn_xfer(p)         = vegmatrixn_rt(ideadstem_xf)
               livecrootn(p)              = vegmatrixn_rt(ilivecroot)
!               livecrootn_storage(p)      = vegmatrixn_rt(ilivecroot_st)
!               livecrootn_xfer(p)         = vegmatrixn_rt(ilivecroot_xf)   
               deadcrootn(p)              = vegmatrixn_rt(ideadcroot)
!               deadcrootn_storage(p)      = vegmatrixn_rt(ideadcroot_st)
!               deadcrootn_xfer(p)         = vegmatrixn_rt(ideadcroot_xf)
               if(ivt(p) >= npcropmin)then
               grainn(p)                  = vegmatrixn_rt(igrain)
               end if
               end if
               if(is_outmatrix)then
                  matrix_cap_leafc(p)                  = vegmatrixc_rt(ileaf)
                  matrix_cap_leafc_storage(p)          = vegmatrixc_rt(ileaf_st)
                  matrix_cap_leafc_xfer(p)             = vegmatrixc_rt(ileaf_xf)
                  matrix_cap_frootc(p)                 = vegmatrixc_rt(ifroot)
                  matrix_cap_frootc_storage(p)         = vegmatrixc_rt(ifroot_st)
                  matrix_cap_frootc_xfer(p)            = vegmatrixc_rt(ifroot_xf)
                  matrix_cap_livestemc(p)              = vegmatrixc_rt(ilivestem)
                  matrix_cap_livestemc_storage(p)      = vegmatrixc_rt(ilivestem_st)
                  matrix_cap_livestemc_xfer(p)         = vegmatrixc_rt(ilivestem_xf)
                  matrix_cap_deadstemc(p)              = vegmatrixc_rt(ideadstem)
                  matrix_cap_deadstemc_storage(p)      = vegmatrixc_rt(ideadstem_st)
                  matrix_cap_deadstemc_xfer(p)         = vegmatrixc_rt(ideadstem_xf)
                  matrix_cap_livecrootc(p)             = vegmatrixc_rt(ilivecroot)
                  matrix_cap_livecrootc_storage(p)     = vegmatrixc_rt(ilivecroot_st)
                  matrix_cap_livecrootc_xfer(p)        = vegmatrixc_rt(ilivecroot_xf)   
                  matrix_cap_deadcrootc(p)             = vegmatrixc_rt(ideadcroot)
                  matrix_cap_deadcrootc_storage(p)     = vegmatrixc_rt(ideadcroot_st)
                  matrix_cap_deadcrootc_xfer(p)        = vegmatrixc_rt(ideadcroot_xf) 
                  if(ivt(p) >= npcropmin)then
                     matrix_cap_grainc(p)                 = vegmatrixc_rt(igrain)
                     matrix_cap_grainc_storage(p)         = vegmatrixc_rt(igrain_st)
                     matrix_cap_grainc_xfer(p)            = vegmatrixc_rt(igrain_xf)
                  end if
                  matrix_cap_leafn(p)                  = vegmatrixn_rt(ileaf)
                  matrix_cap_leafn_storage(p)          = vegmatrixn_rt(ileaf_st)
                  matrix_cap_leafn_xfer(p)             = vegmatrixn_rt(ileaf_xf)
                  matrix_cap_frootn(p)                 = vegmatrixn_rt(ifroot)
                  matrix_cap_frootn_storage(p)         = vegmatrixn_rt(ifroot_st)
                  matrix_cap_frootn_xfer(p)            = vegmatrixn_rt(ifroot_xf)
                  matrix_cap_livestemn(p)              = vegmatrixn_rt(ilivestem)
                  matrix_cap_livestemn_storage(p)      = vegmatrixn_rt(ilivestem_st)
                  matrix_cap_livestemn_xfer(p)         = vegmatrixn_rt(ilivestem_xf)
                  matrix_cap_deadstemn(p)              = vegmatrixn_rt(ideadstem)
                  matrix_cap_deadstemn_storage(p)      = vegmatrixn_rt(ideadstem_st)
                  matrix_cap_deadstemn_xfer(p)         = vegmatrixn_rt(ideadstem_xf)
                  matrix_cap_livecrootn(p)             = vegmatrixn_rt(ilivecroot)
                  matrix_cap_livecrootn_storage(p)     = vegmatrixn_rt(ilivecroot_st)
                  matrix_cap_livecrootn_xfer(p)        = vegmatrixn_rt(ilivecroot_xf)   
                  matrix_cap_deadcrootn(p)             = vegmatrixn_rt(ideadcroot)
                  matrix_cap_deadcrootn_storage(p)     = vegmatrixn_rt(ideadcroot_st)
                  if(ivt(p) >= npcropmin)then
                     matrix_cap_grainn(p)                 = vegmatrixn_rt(igrain)
                     matrix_cap_grainn_storage(p)         = vegmatrixn_rt(igrain_st)
                     matrix_cap_grainn_xfer(p)            = vegmatrixn_rt(igrain_xf)
                  end if
                  matrix_pot_leafc(p)                  = matrix_cap_leafc(p) - leafc(p)
                  matrix_pot_leafc_storage(p)          = matrix_cap_leafc_storage(p) - leafc_storage(p)
                  matrix_pot_leafc_xfer(p)             = matrix_cap_leafc_storage(p) - leafc_xfer(p)
                  matrix_pot_frootc(p)                 = matrix_cap_frootc(p) - frootc(p)
                  matrix_pot_frootc_storage(p)         = matrix_cap_frootc_storage(p) - frootc_storage(p)
                  matrix_pot_frootc_xfer(p)            = matrix_cap_frootc_storage(p) - frootc_xfer(p)
                  matrix_pot_livestemc(p)              = matrix_cap_livestemc(p) - livestemc(p)
                  matrix_pot_livestemc_storage(p)      = matrix_cap_livestemc_storage(p) - livestemc_storage(p)
                  matrix_pot_livestemc_xfer(p)         = matrix_cap_livestemc_storage(p) - livestemc_xfer(p)
                  matrix_pot_deadstemc(p)              = matrix_cap_deadstemc(p) - deadstemc(p)
                  matrix_pot_deadstemc_storage(p)      = matrix_cap_deadstemc_storage(p) - deadstemc_storage(p)
                  matrix_pot_deadstemc_xfer(p)         = matrix_cap_deadstemc_storage(p) - deadstemc_xfer(p)
                  matrix_pot_livecrootc(p)             = matrix_cap_livecrootc(p) - livecrootc(p)
                  matrix_pot_livecrootc_storage(p)     = matrix_cap_livecrootc_storage(p) - livecrootc_storage(p)
                  matrix_pot_livecrootc_xfer(p)        = matrix_cap_livecrootc_storage(p) - livecrootc_xfer(p)
                  matrix_pot_deadcrootc(p)             = matrix_cap_deadcrootc(p) - deadcrootc(p)
                  matrix_pot_deadcrootc_storage(p)     = matrix_cap_deadcrootc_storage(p) - deadcrootc_storage(p)
                  matrix_pot_deadcrootc_xfer(p)        = matrix_cap_deadcrootc_storage(p) - deadcrootc_xfer(p)
                  if(ivt(p) >= npcropmin)then
                  matrix_pot_grainc(p)                 = matrix_cap_grainc(p) - grainc(p)
                  matrix_pot_grainc_storage(p)         = matrix_cap_grainc_storage(p) - grainc_storage(p)
                  matrix_pot_grainc_xfer(p)            = matrix_cap_grainc_storage(p) - grainc_xfer(p)
                  end if
                  matrix_pot_leafn(p)                  = matrix_cap_leafn(p) - leafn(p)
                  matrix_pot_leafn_storage(p)          = matrix_cap_leafn_storage(p) - leafn_storage(p)
                  matrix_pot_leafn_xfer(p)             = matrix_cap_leafn_storage(p) - leafn_xfer(p)
                  matrix_pot_frootn(p)                 = matrix_cap_frootn(p) - frootn(p)
                  matrix_pot_frootn_storage(p)         = matrix_cap_frootn_storage(p) - frootn_storage(p)
                  matrix_pot_frootn_xfer(p)            = matrix_cap_frootn_storage(p) - frootn_xfer(p)
                  matrix_pot_livestemn(p)              = matrix_cap_livestemn(p) - livestemn(p)
                  matrix_pot_livestemn_storage(p)      = matrix_cap_livestemn_storage(p) - livestemn_storage(p)
                  matrix_pot_livestemn_xfer(p)         = matrix_cap_livestemn_storage(p) - livestemn_xfer(p)
                  matrix_pot_deadstemn(p)              = matrix_cap_deadstemn(p) - deadstemn(p)
                  matrix_pot_deadstemn_storage(p)      = matrix_cap_deadstemn_storage(p) - deadstemn_storage(p)
                  matrix_pot_deadstemn_xfer(p)         = matrix_cap_deadstemn_storage(p) - deadstemn_xfer(p)
                  matrix_pot_livecrootn(p)             = matrix_cap_livecrootn(p) - livecrootn(p)
                  matrix_pot_livecrootn_storage(p)     = matrix_cap_livecrootn_storage(p) - livecrootn_storage(p)
                  matrix_pot_livecrootn_xfer(p)        = matrix_cap_livecrootn_storage(p) - livecrootn_xfer(p)
                  matrix_pot_deadcrootn(p)             = matrix_cap_deadcrootn(p) - deadcrootn(p)
                  matrix_pot_deadcrootn_storage(p)     = matrix_cap_deadcrootn_storage(p) - deadcrootn_storage(p)
                  matrix_pot_deadcrootn_xfer(p)        = matrix_cap_deadcrootn_storage(p) - deadcrootn_xfer(p)
                  if(ivt(p) >= npcropmin)then
                     matrix_pot_grainn(p)                 = matrix_cap_grainn(p) - grainn(p)
                     matrix_pot_grainn_storage(p)         = matrix_cap_grainn_storage(p) - grainn_storage(p)
                     matrix_pot_grainn_xfer(p)            = matrix_cap_grainn_storage(p) - grainn_xfer(p)
                  end if
               end if
                 
               matrix_calloc_leaf_acc(p)                = 0._r8
               matrix_calloc_leafst_acc(p)               = 0._r8 
               matrix_calloc_froot_acc(p)               = 0._r8 
               matrix_calloc_frootst_acc(p)                = 0._r8
               matrix_calloc_livestem_acc(p)                = 0._r8
               matrix_calloc_livestemst_acc(p)                = 0._r8
               matrix_calloc_deadstem_acc(p)                = 0._r8
               matrix_calloc_deadstemst_acc(p)                = 0._r8
               matrix_calloc_livecroot_acc(p)                = 0._r8
               matrix_calloc_livecrootst_acc(p)                = 0._r8
               matrix_calloc_deadcroot_acc(p)                = 0._r8
               matrix_calloc_deadcrootst_acc(p)                = 0._r8
               if(ivt(p) >= npcropmin)then
                  matrix_calloc_grain_acc(p)                = 0._r8
                  matrix_calloc_grainst_acc(p)                = 0._r8
               end if

               matrix_ctransfer_leafst_to_leafxf_acc(p) = 0._r8
               matrix_ctransfer_leafxf_to_leaf_acc(p) = 0._r8
               matrix_ctransfer_frootst_to_frootxf_acc(p) = 0._r8
               matrix_ctransfer_frootxf_to_froot_acc(p) = 0._r8
               matrix_ctransfer_livestemst_to_livestemxf_acc(p) = 0._r8
               matrix_ctransfer_livestemxf_to_livestem_acc(p) = 0._r8
               matrix_ctransfer_deadstemst_to_deadstemxf_acc(p) = 0._r8
               matrix_ctransfer_deadstemxf_to_deadstem_acc(p) = 0._r8
               matrix_ctransfer_livecrootst_to_livecrootxf_acc(p) = 0._r8
               matrix_ctransfer_livecrootxf_to_livecroot_acc(p) = 0._r8
               matrix_ctransfer_deadcrootst_to_deadcrootxf_acc(p) = 0._r8
               matrix_ctransfer_deadcrootxf_to_deadcroot_acc(p) = 0._r8
               if(ivt(p) >= npcropmin)then
                  matrix_ctransfer_grainst_to_grainxf_acc(p) = 0._r8
                  matrix_ctransfer_grainxf_to_grain_acc(p) = 0._r8
               end if
               matrix_ctransfer_livestem_to_deadstem_acc(p) = 0._r8
               matrix_ctransfer_livecroot_to_deadcroot_acc(p) = 0._r8

               matrix_cturnover_leaf_acc(p)                = 0._r8
               matrix_cturnover_leafst_acc(p)                = 0._r8
               matrix_cturnover_leafxf_acc(p)                = 0._r8
               matrix_cturnover_froot_acc(p)                = 0._r8
               matrix_cturnover_frootst_acc(p)                = 0._r8
               matrix_cturnover_frootxf_acc(p)                = 0._r8
               matrix_cturnover_livestem_acc(p)                = 0._r8
               matrix_cturnover_livestemst_acc(p)                = 0._r8
               matrix_cturnover_livestemxf_acc(p)                = 0._r8
               matrix_cturnover_deadstem_acc(p)                = 0._r8
               matrix_cturnover_deadstemst_acc(p)                = 0._r8
               matrix_cturnover_deadstemxf_acc(p)                = 0._r8
               matrix_cturnover_livecroot_acc(p)                = 0._r8
               matrix_cturnover_livecrootst_acc(p)                = 0._r8
               matrix_cturnover_livecrootxf_acc(p)                = 0._r8
               matrix_cturnover_deadcroot_acc(p)                = 0._r8
               matrix_cturnover_deadcrootst_acc(p)                = 0._r8
               matrix_cturnover_deadcrootxf_acc(p)                = 0._r8
               if(ivt(p) >= npcropmin)then
                  matrix_cturnover_grain_acc(p)               = 0._r8 
                  matrix_cturnover_grainst_acc(p)                = 0._r8
                  matrix_cturnover_grainxf_acc(p)                = 0._r8
               end if

               matrix_nalloc_leaf_acc(p)                = 0._r8
               matrix_nalloc_leafst_acc(p)                = 0._r8
               matrix_nalloc_froot_acc(p)                = 0._r8
               matrix_nalloc_frootst_acc(p)                = 0._r8
               matrix_nalloc_livestem_acc(p)                = 0._r8
               matrix_nalloc_livestemst_acc(p)                = 0._r8
               matrix_nalloc_deadstem_acc(p)                = 0._r8
               matrix_nalloc_deadstemst_acc(p)                = 0._r8
               matrix_nalloc_livecroot_acc(p)                = 0._r8
               matrix_nalloc_livecrootst_acc(p)                = 0._r8
               matrix_nalloc_deadcroot_acc(p)                = 0._r8
               matrix_nalloc_deadcrootst_acc(p)                = 0._r8
               if(ivt(p) >= npcropmin)then
                  matrix_nalloc_grain_acc(p)                = 0._r8
                  matrix_nalloc_grainst_acc(p)                = 0._r8
               end if

               matrix_ntransfer_leafst_to_leafxf_acc(p) = 0._r8
               matrix_ntransfer_leafxf_to_leaf_acc(p) = 0._r8
               matrix_ntransfer_frootst_to_frootxf_acc(p) = 0._r8
               matrix_ntransfer_frootxf_to_froot_acc(p) = 0._r8
               matrix_ntransfer_livestemst_to_livestemxf_acc(p) = 0._r8
               matrix_ntransfer_livestemxf_to_livestem_acc(p) = 0._r8
               matrix_ntransfer_deadstemst_to_deadstemxf_acc(p) = 0._r8
               matrix_ntransfer_deadstemxf_to_deadstem_acc(p) = 0._r8
               matrix_ntransfer_livecrootst_to_livecrootxf_acc(p) = 0._r8
               matrix_ntransfer_livecrootxf_to_livecroot_acc(p) = 0._r8
               matrix_ntransfer_deadcrootst_to_deadcrootxf_acc(p) = 0._r8
               matrix_ntransfer_deadcrootxf_to_deadcroot_acc(p) = 0._r8
               if(ivt(p) >= npcropmin)then
                  matrix_ntransfer_grainst_to_grainxf_acc(p) = 0._r8
                  matrix_ntransfer_grainxf_to_grain_acc(p) = 0._r8
               end if
               matrix_ntransfer_livestem_to_deadstem_acc(p) = 0._r8
               matrix_ntransfer_livecroot_to_deadcroot_acc(p) = 0._r8

               matrix_ntransfer_retransn_to_leaf_acc(p)               = 0._r8 
               matrix_ntransfer_retransn_to_leafst_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_froot_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_frootst_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_livestem_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_livestemst_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_deadstem_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_deadstemst_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_livecroot_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_livecrootst_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_deadcroot_acc(p)                = 0._r8
               matrix_ntransfer_retransn_to_deadcrootst_acc(p)                = 0._r8
               if(ivt(p) >= npcropmin)then
                  matrix_ntransfer_retransn_to_grain_acc(p)                = 0._r8
                  matrix_ntransfer_retransn_to_grainst_acc(p)                = 0._r8
               end if
               matrix_ntransfer_leaf_to_retransn_acc(p)                = 0._r8
               matrix_ntransfer_froot_to_retransn_acc(p)                = 0._r8
               matrix_ntransfer_livestem_to_retransn_acc(p)                = 0._r8
               matrix_ntransfer_livecroot_to_retransn_acc(p)                = 0._r8

               matrix_nturnover_leaf_acc(p)                = 0._r8
               matrix_nturnover_leafst_acc(p)                = 0._r8
               matrix_nturnover_leafxf_acc(p)                = 0._r8
               matrix_nturnover_froot_acc(p)                = 0._r8
               matrix_nturnover_frootst_acc(p)                = 0._r8
               matrix_nturnover_frootxf_acc(p)                = 0._r8
               matrix_nturnover_livestem_acc(p)                = 0._r8
               matrix_nturnover_livestemst_acc(p)                = 0._r8
               matrix_nturnover_livestemxf_acc(p)                = 0._r8
               matrix_nturnover_deadstem_acc(p)                = 0._r8
               matrix_nturnover_deadstemst_acc(p)                = 0._r8
               matrix_nturnover_deadstemxf_acc(p)                = 0._r8
               matrix_nturnover_livecroot_acc(p)                = 0._r8
               matrix_nturnover_livecrootst_acc(p)                = 0._r8
               matrix_nturnover_livecrootxf_acc(p)                = 0._r8
               matrix_nturnover_deadcroot_acc(p)                = 0._r8
               matrix_nturnover_deadcrootst_acc(p)                = 0._r8
               matrix_nturnover_deadcrootxf_acc(p)                = 0._r8
               if(ivt(p) >= npcropmin)then
                  matrix_nturnover_grain_acc(p)                = 0._r8
                  matrix_nturnover_grainst_acc(p)                = 0._r8
                  matrix_nturnover_grainxf_acc(p)                = 0._r8
               end if
               matrix_nturnover_retransn_acc(p)                = 0._r8
!
               matrix_calloc_acc(:) = 0._r8
               matrix_ctransfer_acc(:,:) = 0._r8


               matrix_nalloc_acc(:) = 0._r8
               matrix_ntransfer_acc(:,:) = 0._r8
! C       
           end if
         end if

         matrix_Cinput(p) = 0._r8
         matrix_phtransfer(p,:,:) = 0._r8
         matrix_gmtransfer(p,:,:) = 0._r8
         matrix_fitransfer(p,:,:) = 0._r8
         matrix_nphtransfer(p,:,:) = 0._r8
         matrix_ngmtransfer(p,:,:) = 0._r8
         matrix_Ninput(p) = 0._r8
         
      end do 
    
   end associate 
 end subroutine CNVegMatrix
 subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer,intent(in) :: n
real(r8),intent(in)  :: a(n,n)
real(r8),intent(out) :: c(n,n)
real(r8) :: L(n,n), U(n,n), aa(n,n), b(n), d(n), x(n)
real(r8) :: coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

aa=a

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=aa(i,k)/aa(k,k)
      L(i,k) = coeff
      do j=k+1,n
         aa(i,j) = aa(i,j)-coeff*aa(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = aa(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

end module CNVegMatrixMod
