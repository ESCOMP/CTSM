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
    use clm_varcon                   , only: secspday
  use pftconMod                      , only : pftcon
  use CNVegCarbonStateType           , only : cnveg_carbonstate_type
  use CNVegNitrogenStateType         , only : cnveg_nitrogenstate_type
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type     !include: callocation,ctransfer, cturnover
  use CNVegNitrogenFluxType          , only : cnveg_nitrogenflux_type
  use CNVegStateType                 , only : cnveg_state_type
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
                        cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,cnveg_state_inst)
    ! !DESCRIPTION:
    ! !ARGUMENTS:
     type(bounds_type)                      , intent(in)    :: bounds
     integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
     integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
     type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst  
     type(cnveg_nitrogenstate_type)           , intent(inout) :: cnveg_nitrogenstate_inst
     type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
     type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
 
     type(cnveg_state_type)    , intent(in) :: cnveg_state_inst
!	
    ! !LOCAL VARIABLES:
     integer :: fc,fp,j,i    ! indices
     integer :: p,c         !  
     real(r8),allocatable,dimension(:,:) :: vegmatrixc_old(:,:)
     real(r8),allocatable,dimension(:,:) :: vegmatrixc_new(:,:)
     real(r8),allocatable,dimension(:,:) :: vegmatrixn_old(:,:)
     real(r8),allocatable,dimension(:,:) :: vegmatrixn_new(:,:)
!     real(r8),allocatable,dimension(:,:) :: vegmatrixc_rt(:,:),vegmatrixn_rt(:,:)
     real(r8),allocatable,dimension(:,:,:) :: matrix_nphturnover(:,:,:),matrix_ngmturnover(:,:,:)!,matrix_nphturnover(:,:,:)
     real(r8),allocatable,dimension(:,:,:) :: matrix_nphtransfer_curr(:,:,:)

! for spinupacc
     real(r8),dimension(1:nvegcpool) :: vegmatrixc_rt
     real(r8),dimension(1:nvegnpool) :: vegmatrixn_rt,tmp
     real(r8),allocatable,dimension(:,:) :: vegmatrix_2d(:,:),AKinv,vegmatrixn_2d(:,:),AKinvn
     logical  ::  end_of_year
     integer, save :: counter=0

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
         matrix_alloc_acc       => cnveg_carbonstate_inst%matrix_alloc_acc_patch             , & !
         matrix_transfer_acc    => cnveg_carbonstate_inst%matrix_transfer_acc_patch          , & !
         matrix_nalloc_acc      => cnveg_nitrogenstate_inst%matrix_nalloc_acc_patch             , & !
         matrix_ntransfer_acc   => cnveg_nitrogenstate_inst%matrix_ntransfer_acc_patch          , & !
 
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

      allocate(vegmatrixc_old(bounds%begp:bounds%endp,nvegcpool))
      allocate(vegmatrixc_new(bounds%begp:bounds%endp,nvegcpool))
      allocate(vegmatrixn_old(bounds%begp:bounds%endp,nvegnpool))
      allocate(vegmatrixn_new(bounds%begp:bounds%endp,nvegnpool))
!      allocate(vegmatrixc_rt(bounds%begp:bounds%endp,nvegpool))
!      allocate(vegmatrixn_rt(bounds%begp:bounds%endp,nvegpool))
  
!     allocate(matrix_nphturnover(bounds%begp:bounds%endp,nvegpool,nvegpool))
!     allocate(matrix_ngmturnover(bounds%begp:bounds%endp,nvegpool,nvegpool))
!     allocate(matrix_nfitransfer(bounds%begp:bounds%endp,nvegpool+1,nvegpool))

      allocate(vegmatrix_2d(nvegcpool,nvegcpool))
      allocate(AKinv(nvegcpool,nvegcpool))

      allocate(vegmatrixn_2d(nvegnpool,nvegnpool))
      allocate(AKinvn(nvegnpool,nvegnpool))

      counter = counter + dt
     if (counter >=1.0 * secspyear) then ! link to the recycling span of climate forcing
          end_of_year = .true.
          counter = 0._r8
       else
          end_of_year = .false.
       end if

      do fp = 1,num_soilp
         vegmatrix_2d(:,:)=0._r8
         vegmatrixn_2d(:,:)=0._r8
         AKinv(:,:)=0._r8
         AKinvn(:,:)=0._r8

         p = filter_soilp(fp)
         vegmatrixc_old(p,ileaf)         = leafc(p)
         vegmatrixc_old(p,ileaf_st)      = leafc_storage(p)
         vegmatrixc_old(p,ileaf_xf)      = leafc_xfer(p)
         vegmatrixc_old(p,ifroot)        = frootc(p)
         vegmatrixc_old(p,ifroot_st)     = frootc_storage(p)
         vegmatrixc_old(p,ifroot_xf)     = frootc_xfer(p)
         vegmatrixc_old(p,ilivestem)     = livestemc(p)
         vegmatrixc_old(p,ilivestem_st)  = livestemc_storage(p)
         vegmatrixc_old(p,ilivestem_xf)  = livestemc_xfer(p)
         vegmatrixc_old(p,ideadstem)     = deadstemc(p)
         vegmatrixc_old(p,ideadstem_st)  = deadstemc_storage(p)
         vegmatrixc_old(p,ideadstem_xf)  = deadstemc_xfer(p)
         vegmatrixc_old(p,ilivecroot)    = livecrootc(p)
         vegmatrixc_old(p,ilivecroot_st) = livecrootc_storage(p)
         vegmatrixc_old(p,ilivecroot_xf) = livecrootc_xfer(p)
         vegmatrixc_old(p,ideadcroot)    = deadcrootc(p)
         vegmatrixc_old(p,ideadcroot_st) = deadcrootc_storage(p)
         vegmatrixc_old(p,ideadcroot_xf) = deadcrootc_xfer(p)
         vegmatrixc_old(p,igrain)        = grainc(p)
         vegmatrixc_old(p,igrain_st)     = grainc_storage(p)
         vegmatrixc_old(p,igrain_xf)     = grainc_xfer(p)
 
         vegmatrixn_old(p,ileaf)         = leafn(p)
         vegmatrixn_old(p,ileaf_st)      = leafn_storage(p)
         vegmatrixn_old(p,ileaf_xf)      = leafn_xfer(p)
         vegmatrixn_old(p,ifroot)        = frootn(p)
         vegmatrixn_old(p,ifroot_st)     = frootn_storage(p)
         vegmatrixn_old(p,ifroot_xf)     = frootn_xfer(p)
         vegmatrixn_old(p,ilivestem)     = livestemn(p)
         vegmatrixn_old(p,ilivestem_st)  = livestemn_storage(p)
         vegmatrixn_old(p,ilivestem_xf)  = livestemn_xfer(p)
         vegmatrixn_old(p,ideadstem)     = deadstemn(p)
         vegmatrixn_old(p,ideadstem_st)  = deadstemn_storage(p)
         vegmatrixn_old(p,ideadstem_xf)  = deadstemn_xfer(p)
         vegmatrixn_old(p,ilivecroot)    = livecrootn(p)
         vegmatrixn_old(p,ilivecroot_st) = livecrootn_storage(p)
         vegmatrixn_old(p,ilivecroot_xf) = livecrootn_xfer(p)
         vegmatrixn_old(p,ideadcroot)    = deadcrootn(p)
         vegmatrixn_old(p,ideadcroot_st) = deadcrootn_storage(p)
         vegmatrixn_old(p,ideadcroot_xf) = deadcrootn_xfer(p)
         vegmatrixn_old(p,igrain)        = grainn(p)
         vegmatrixn_old(p,igrain_st)     = grainn_storage(p)
         vegmatrixn_old(p,igrain_xf)     = grainn_xfer(p)
         vegmatrixn_old(p,iretransn)     = retransn(p)
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
            grainc0(p)               = max(grainc(p), 1.e-15_r8)
            grainc0_storage(p)       = max(grainc_storage(p), 1.e-15_r8)
            grainc0_xfer(p)          = max(grainc_xfer(p), 1.e-15_r8)
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
            grainn0(p)               = max(grainn(p), 1.e-15_r8)
            grainn0_storage(p)       = max(grainn_storage(p), 1.e-15_r8)
            grainn0_xfer(p)          = max(grainn_xfer(p), 1.e-15_r8)
         end if
!
          vegmatrix_2d(ileaf,ileaf)               = leafc(p)         / leafc0(p)
          vegmatrix_2d(ileaf_st,ileaf_st)         = leafc_storage(p) / leafc0_storage(p)
          vegmatrix_2d(ileaf_xf,ileaf_xf)         = leafc_xfer(p)    / leafc0_xfer(p)
          vegmatrix_2d(ifroot,ifroot)             = frootc(p)        / frootc0(p)
          vegmatrix_2d(ifroot_st,ifroot_st)       = frootc_storage(p)/ frootc0_storage(p)
          vegmatrix_2d(ifroot_xf,ifroot_xf)       = frootc_xfer(p)   / frootc0_xfer(p)
          vegmatrix_2d(ilivestem,ilivestem)       = livestemc(p)     / livestemc0(p)
          vegmatrix_2d(ilivestem_st,ilivestem_st) = livestemc_storage(p) / livestemc0_storage(p)
          vegmatrix_2d(ilivestem_xf,ilivestem_xf) = livestemc_xfer(p)    / livestemc0_xfer(p)
          vegmatrix_2d(ideadstem,ideadstem)       = deadstemc(p)     / deadstemc0(p)
          vegmatrix_2d(ideadstem_st,ideadstem_st) = deadstemc_storage(p)/ deadstemc0_storage(p)
          vegmatrix_2d(ideadstem_xf,ideadstem_xf) = deadstemc_xfer(p) / deadstemc0_xfer(p)
          vegmatrix_2d(ilivecroot,ilivecroot)     = livecrootc(p)    / livecrootc0(p)
          vegmatrix_2d(ilivecroot_st,ilivecroot_st) = livecrootc_storage(p)/livecrootc0_storage(p)
          vegmatrix_2d(ilivecroot_xf,ilivecroot_xf) = livecrootc_xfer(p)/ livecrootc0_xfer(p)
          vegmatrix_2d(ideadcroot,ideadcroot)       = deadcrootc(p)  / deadcrootc0(p)
          vegmatrix_2d(ideadcroot_st,ideadcroot_st) = deadcrootc_storage(p)/deadcrootc0_storage(p)
          vegmatrix_2d(ideadcroot_xf,ideadcroot_xf) = deadcrootc_xfer(p)/deadcrootc0_xfer(p)
          vegmatrix_2d(igrain,igrain)             = grainc(p)        / grainc0(p)
          vegmatrix_2d(igrain_st,igrain_st)       = grainc_storage(p)/ grainc0_storage(p)
          vegmatrix_2d(igrain_xf,igrain_xf)       = grainc_xfer(p)   / grainc0_xfer(p)

          vegmatrixn_2d(ileaf,ileaf)               = leafn(p)         / leafn0(p)
          vegmatrixn_2d(ileaf_st,ileaf_st)         = leafn_storage(p) / leafn0_storage(p)
          vegmatrixn_2d(ileaf_xf,ileaf_xf)         = leafn_xfer(p)    / leafn0_xfer(p)
          vegmatrixn_2d(ifroot,ifroot)             = frootn(p)        / frootn0(p)
          vegmatrixn_2d(ifroot_st,ifroot_st)       = frootn_storage(p)/ frootn0_storage(p)
          vegmatrixn_2d(ifroot_xf,ifroot_xf)       = frootn_xfer(p)   / frootn0_xfer(p)
          vegmatrixn_2d(ilivestem,ilivestem)       = livestemn(p)     / livestemn0(p)
          vegmatrixn_2d(ilivestem_st,ilivestem_st) = livestemn_storage(p) / livestemn0_storage(p)
          vegmatrixn_2d(ilivestem_xf,ilivestem_xf) = livestemn_xfer(p)    / livestemn0_xfer(p)
          vegmatrixn_2d(ideadstem,ideadstem)       = deadstemn(p)     / deadstemn0(p)
          vegmatrixn_2d(ideadstem_st,ideadstem_st) = deadstemn_storage(p)/ deadstemn0_storage(p)
          vegmatrixn_2d(ideadstem_xf,ideadstem_xf) = deadstemn_xfer(p) / deadstemn0_xfer(p)
          vegmatrixn_2d(ilivecroot,ilivecroot)     = livecrootn(p)    / livecrootn0(p)
          vegmatrixn_2d(ilivecroot_st,ilivecroot_st) = livecrootn_storage(p)/livecrootn0_storage(p)
          vegmatrixn_2d(ilivecroot_xf,ilivecroot_xf) = livecrootn_xfer(p)/livecrootn0_xfer(p)
          vegmatrixn_2d(ideadcroot,ideadcroot)       = deadcrootn(p)  / deadcrootn0(p)
          vegmatrixn_2d(ideadcroot_st,ideadcroot_st) = deadcrootn_storage(p)/deadcrootn0_storage(p)
          vegmatrixn_2d(ideadcroot_xf,ideadcroot_xf) = deadcrootn_xfer(p)/deadcrootn0_xfer(p)
          vegmatrixn_2d(igrain,igrain)             = grainn(p)        / grainn0(p)
          vegmatrixn_2d(igrain_st,igrain_st)       = grainn_storage(p)/ grainn0_storage(p)
          vegmatrixn_2d(igrain_xf,igrain_xf)       = grainn_xfer(p)   / grainn0_xfer(p)
!
         matrix_phturnover(p,:,:) = 0._r8
         matrix_gmturnover(p,:,:) = 0._r8
         matrix_fiturnover(p,:,:) = 0._r8
         matrix_nphturnover(p,:,:) = 0._r8
         matrix_ngmturnover(p,:,:) = 0._r8
         matrix_nfiturnover(p,:,:) = 0._r8

!         matrix_nphtransfer_curr(p,:,:) = matrix_phtransfer(p,:,:)
!         matrix_nphtransfer_curr(p,ileaf_xf,ileaf_st)     = matrix_nphtransfer(p,ileaf_xf,ileaf_st)
!         matrix_nphtransfer_curr(p,ifroot_xf,ifroot_st)   = matrix_nphtransfer(p,ifroot_xf,ifroot_st)
   
!         matrix_nphtransfer_curr(p,ideadstem,ilivestem)    =  matrix_nphtransfer(p,ideadstem,ilivestem)
!         matrix_nphtransfer_curr(p,ideadcroot,ilivecroot)  =  matrix_nphtransfer(p,ideadcroot,ilivecroot)

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
         if(p .eq. 8)print*,'before generate turnover',matrix_ngmturnover(p,ideadstem,ideadstem),matrix_ngmtransfer(p,ioutn,ideadstem)
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
!          if(p.eq.5)print *, 'before1NNN',matrix_alloc(p,:) * matrix_Cinput(p) * dt
!         if(p .eq. 8428)then
!             tmp=(matmul(matmul(matrix_phtransfer(p,1:nvegpool,:),matrix_phturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!             print*,'before matrix_new,phenology',cnveg_carbonstate_inst%deadcrootc_patch(8428),tmp(ideadcroot)
!             tmp=(matmul(matmul(matrix_gmtransfer(p,1:nvegpool,:),matrix_gmturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!             print*,'before matrix_new,gap mortality',cnveg_carbonstate_inst%deadcrootc_patch(8428),tmp(ideadcroot)
!             tmp=(matmul(matmul(matrix_fitransfer(p,1:nvegpool,:),matrix_fiturnover(p,:,:)),vegmatrixc_old(p,:))) * dt
!             print*,'before matrix_new,fire',cnveg_carbonstate_inst%deadcrootc_patch(8428),tmp(ideadcroot)
!         end if
!         if(p .eq. 8428)then
!            print*,'before matrix_new',vegmatrixc_new(p,ideadcroot)
!         end if
         vegmatrixc_new(p,:) = vegmatrixc_old(p,:)  +  matrix_alloc(p,:) * matrix_Cinput(p) * dt &
                            +(matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc_old(p,:))   &
                            + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc_old(p,:))   &
                            + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc_old(p,:))) * dt 
                            
         if(p .eq. 8)then
            !print*,'before matrix_old deadcrootc',p,vegmatrixc_old(p,ideadcroot)
            print*,'before matrix_old deadstemn',p,vegmatrixn_old(p,ideadstem)
!            print*,'before matrix_old leafn',p,vegmatrixn_old(p,ileaf)
!         print*,'before matrix_old',p,vegmatrixn_old(p,:)
!         end if
!N matrix
!            print*,'Cinput to deadcrootn',matrix_Cinput(p)*matrix_alloc(p,ideadcroot)*dt
!            print*,'Ninput to grainn',matrix_Ninput(p)*matrix_nalloc(p,ileaf)*dt,matrix_Ninput(p),matrix_nalloc(p,ileaf)*dt
            print*,'Ninput to deadstemn',matrix_Ninput(p)*matrix_nalloc(p,ideadstem)*dt,matrix_Ninput(p),matrix_nalloc(p,ideadstem)*dt
!            print*,'Transfer to deadcrootc phenology from livecroot',matrix_phtransfer(p,ideadcroot,ilivecroot)*matrix_phturnover(p,ilivecroot,ilivecroot) &
!                                                                 * vegmatrixc_old(p,ilivecroot) * dt
!            print*,'Transfer to deadcrootc phenology from transfer pool',matrix_phtransfer(p,ideadcroot,ideadcroot_xf)*matrix_phturnover(p,ideadcroot_xf,ideadcroot_xf) &
!  * vegmatrixc_old(p,ideadcroot_xf) * dt,vegmatrixc_old(p,ideadcroot_xf),matrix_phturnover(p,ideadcroot_xf,ideadcroot_xf),matrix_phtransfer(p,ideadcroot,ideadcroot_xf) 
!            print*,'Transfer to leafc phenology from transfer pool',matrix_phtransfer(p,ileaf,ileaf_xf)*matrix_phturnover(p,ileaf_xf,ileaf_xf) &
!  * vegmatrixc_old(p,ileaf_xf) * dt,vegmatrixc_old(p,ileaf_xf),matrix_phturnover(p,ileaf_xf,ileaf_xf),matrix_phtransfer(p,ileaf,ileaf_xf) 
!            print*,'Transfer from deadcrootn phenology',matrix_phtransfer(p,ideadcroot,ideadcroot) * matrix_phturnover(p,ideadcroot,ideadcroot) &
!                                                                    * vegmatrixc_old(p,ideadcroot) * dt
            print*,'Transfer from deadstemn_xfer phenology',matrix_nphtransfer(p,ideadstem,ideadstem_xf) * matrix_nphturnover(p,ideadstem_xf,ideadstem_xf) &
                                                                    * vegmatrixn_old(p,ideadstem_xf) * dt
            print*,'Transfer from livestemn phenology',matrix_nphtransfer(p,ideadstem,ilivestem) * matrix_nphturnover(p,ilivestem,ilivestem) &
                                                                    * vegmatrixn_old(p,ilivestem) * dt
            print*,'Transfer from retransn phenology',matrix_nphtransfer(p,ideadstem,iretransn) * matrix_nphturnover(p,iretransn,iretransn) &
                                                                    * vegmatrixn_old(p,iretransn) * dt
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
            print*,'Transfer from deadstemn gap mort',matrix_ngmtransfer(p,ioutn,ideadstem) * matrix_ngmturnover(p,ideadstem,ideadstem) &
                                                                    * vegmatrixn_old(p,ideadstem) * dt, matrix_ngmtransfer(p,ioutn,ideadstem)
            print*,'Transfer from deadstemn fire',matrix_nfitransfer(p,ideadstem,ideadstem) * matrix_nfiturnover(p,ideadstem,ideadstem) &
                                                                    * vegmatrixn_old(p,ideadstem) * dt
            print*,'Transfer from livestemn fire',matrix_nfitransfer(p,ideadstem,ilivestem) * matrix_nfiturnover(p,ideadstem,ilivestem) &
                                                                    * vegmatrixn_old(p,ilivestem) * dt
         end if
         vegmatrixn_new(p,:) = vegmatrixn_old(p,:)  +  matrix_nalloc(p,:) * matrix_Ninput(p) * dt &
                            + (matmul(matmul(matrix_nphtransfer(p,1:nvegnpool,:),matrix_nphturnover(p,:,:)),vegmatrixn_old(p,:)) &
                            +  matmul(matmul(matrix_ngmtransfer(p,1:nvegnpool,:),matrix_ngmturnover(p,:,:)),vegmatrixn_old(p,:)) &
                            +  matmul(matmul(matrix_nfitransfer(p,1:nvegnpool,:),matrix_nfiturnover(p,:,:)),vegmatrixn_old(p,:))) * dt!&
         if(p .eq. 8)then
            print*,'after matrix_new',p,vegmatrixn_new(p,ideadstem)
         end if
 
! 

         if(isspinup .or. is_outmatrix)then
            matrix_alloc_acc(p,1:nvegcpool) = matrix_alloc_acc(p,1:nvegcpool) + matrix_alloc(p,:) * matrix_Cinput(p) * dt!*(dt/secspyear)
            matrix_transfer_acc(p,1:nvegcpool,:) = matrix_transfer_acc(p,1:nvegcpool,:) &
                               + (matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrix_2d(:,:))& !) * dt
                               + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrix_2d(:,:)) &
                               + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrix_2d(:,:)))*dt

            matrix_nalloc_acc(p,1:nvegnpool) = matrix_nalloc_acc(p,1:nvegnpool) + matrix_nalloc(p,:) * matrix_Ninput(p)* dt
            matrix_ntransfer_acc(p,1:nvegnpool,:) = matrix_ntransfer_acc(p,1:nvegnpool,:) &
                               + (matmul(matmul(matrix_nphtransfer(p,1:nvegnpool,:),matrix_nphturnover(p,:,:)),vegmatrixn_2d(:,:))& !) * dt
                               + matmul(matmul(matrix_gmtransfer(p,1:nvegnpool,:),matrix_gmturnover(p,:,:)),vegmatrixn_2d(:,:))) * dt

         end if 

         leafc(p)               = vegmatrixc_new(p,ileaf)
         leafc_storage(p)       = vegmatrixc_new(p,ileaf_st)
         leafc_xfer(p)          = vegmatrixc_new(p,ileaf_xf)
         frootc(p)              = vegmatrixc_new(p,ifroot)
         frootc_storage(p)      = vegmatrixc_new(p,ifroot_st)
         frootc_xfer(p)         = vegmatrixc_new(p,ifroot_xf)
         livestemc(p)           = vegmatrixc_new(p,ilivestem)
         livestemc_storage(p)   = vegmatrixc_new(p,ilivestem_st)
         livestemc_xfer(p)      = vegmatrixc_new(p,ilivestem_xf)
         deadstemc(p)           = vegmatrixc_new(p,ideadstem)
         deadstemc_storage(p)   = vegmatrixc_new(p,ideadstem_st)
         deadstemc_xfer(p)      = vegmatrixc_new(p,ideadstem_xf)
         livecrootc(p)          = vegmatrixc_new(p,ilivecroot)
         livecrootc_storage(p)  = vegmatrixc_new(p,ilivecroot_st)
         livecrootc_xfer(p)     = vegmatrixc_new(p,ilivecroot_xf)
         deadcrootc(p)          = vegmatrixc_new(p,ideadcroot)
         deadcrootc_storage(p)  = vegmatrixc_new(p,ideadcroot_st)
         deadcrootc_xfer(p)     = vegmatrixc_new(p,ideadcroot_xf)
         grainc(p)              = vegmatrixc_new(p,igrain)
         grainc_storage(p)      = vegmatrixc_new(p,igrain_st)
         grainc_xfer(p)         = vegmatrixc_new(p,igrain_xf)
 
         leafn(p)               = vegmatrixn_new(p,ileaf)
         leafn_storage(p)       = vegmatrixn_new(p,ileaf_st)
         leafn_xfer(p)          = vegmatrixn_new(p,ileaf_xf)
         frootn(p)              = vegmatrixn_new(p,ifroot)
         frootn_storage(p)      = vegmatrixn_new(p,ifroot_st)
         frootn_xfer(p)         = vegmatrixn_new(p,ifroot_xf)
         livestemn(p)           = vegmatrixn_new(p,ilivestem)
         livestemn_storage(p)   = vegmatrixn_new(p,ilivestem_st)
         livestemn_xfer(p)      = vegmatrixn_new(p,ilivestem_xf)
         deadstemn(p)           = vegmatrixn_new(p,ideadstem)
         deadstemn_storage(p)   = vegmatrixn_new(p,ideadstem_st)
         deadstemn_xfer(p)      = vegmatrixn_new(p,ideadstem_xf)
         livecrootn(p)          = vegmatrixn_new(p,ilivecroot)
         livecrootn_storage(p)  = vegmatrixn_new(p,ilivecroot_st)
         livecrootn_xfer(p)     = vegmatrixn_new(p,ilivecroot_xf)
         deadcrootn(p)          = vegmatrixn_new(p,ideadcroot)
         deadcrootn_storage(p)  = vegmatrixn_new(p,ideadcroot_st)
         deadcrootn_xfer(p)     = vegmatrixn_new(p,ideadcroot_xf)
         grainn(p)              = vegmatrixn_new(p,igrain)
         grainn_storage(p)      = vegmatrixn_new(p,igrain_st)
         grainn_xfer(p)         = vegmatrixn_new(p,igrain_xf)
         retransn(p)            = vegmatrixn_new(p,iretransn)
         if(isspinup .or. is_outmatrix)then
!print *, 'count',counter, end_of_year
            if(end_of_year)then
               do i=1,nvegcpool
                 if(matrix_transfer_acc(p,i,i) .eq. 0)then
                    matrix_transfer_acc(p,i,i) = 1.e+36
                 end if
               end do
               do i=1,nvegnpool
                 if(matrix_ntransfer_acc(p,i,i) .eq. 0)then
                    matrix_ntransfer_acc(p,i,i) = 1.e+36
                 end if
               end do
!              end if
               call inverse(matrix_transfer_acc(p,1:nvegcpool,:),AKinv(1:nvegcpool,1:nvegcpool),nvegcpool)
               vegmatrixc_rt(:) = -matmul(AKinv(1:nvegcpool,1:nvegcpool),matrix_alloc_acc(p,:))
! N  
               call inverse(matrix_ntransfer_acc(p,1:nvegnpool,:),AKinvn(1:nvegnpool,1:nvegnpool),nvegnpool)
               vegmatrixn_rt(:) = -matmul(AKinvn(1:nvegnpool,1:nvegnpool),matrix_nalloc_acc(p,:))
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
               grainc(p)                 = vegmatrixc_rt(igrain)
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
               grainn(p)                  = vegmatrixn_rt(igrain)
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
                  matrix_cap_grainc(p)                 = vegmatrixc_rt(igrain)
                  matrix_cap_grainc_storage(p)         = vegmatrixc_rt(igrain_st)
                  matrix_cap_grainc_xfer(p)            = vegmatrixc_rt(igrain_xf)
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
                  matrix_cap_grainn(p)                 = vegmatrixn_rt(igrain)
                  matrix_cap_grainn_storage(p)         = vegmatrixn_rt(igrain_st)
                  matrix_cap_grainn_xfer(p)            = vegmatrixn_rt(igrain_xf)
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
                  matrix_pot_grainc(p)                 = matrix_cap_grainc(p) - grainc(p)
                  matrix_pot_grainc_storage(p)         = matrix_cap_grainc_storage(p) - grainc_storage(p)
                  matrix_pot_grainc_xfer(p)            = matrix_cap_grainc_storage(p) - grainc_xfer(p)
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
                  matrix_pot_grainn(p)                 = matrix_cap_grainn(p) - grainn(p)
                  matrix_pot_grainn_storage(p)         = matrix_cap_grainn_storage(p) - grainn_storage(p)
                  matrix_pot_grainn_xfer(p)            = matrix_cap_grainn_storage(p) - grainn_xfer(p)
               end if
                 
!
               matrix_alloc_acc(p,:) = 0._r8
               matrix_transfer_acc(p,:,:) = 0._r8
               matrix_nalloc_acc(p,:) = 0._r8
               matrix_ntransfer_acc(p,:,:) = 0._r8
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
