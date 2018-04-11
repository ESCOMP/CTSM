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
  use pftconMod                      , only : pftcon,npcropmin
  use CNVegCarbonStateType           , only : cnveg_carbonstate_type
  use CNVegNitrogenStateType         , only : cnveg_nitrogenstate_type
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type     !include: callocation,ctransfer, cturnover
  use CNVegNitrogenFluxType          , only : cnveg_nitrogenflux_type
  use CNVegStateType                 , only : cnveg_state_type
  use clm_varctl                     , only : isspinup, is_outmatrix
  use clm_varctl                     , only : use_c13, use_c14 
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
                          cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,cnveg_state_inst, &
                          c13_cnveg_carbonstate_inst,c14_cnveg_carbonstate_inst,c13_cnveg_carbonflux_inst,&
                          c14_cnveg_carbonflux_inst)
    ! !DESCRIPTION:
    ! !ARGUMENTS:
     type(bounds_type)                      , intent(in)    :: bounds
     integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
     integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
     type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst  
     type(cnveg_nitrogenstate_type)         , intent(inout) :: cnveg_nitrogenstate_inst
     type(cnveg_carbonflux_type)            , intent(inout) :: cnveg_carbonflux_inst
     type(cnveg_nitrogenflux_type)          , intent(inout) :: cnveg_nitrogenflux_inst
     type(cnveg_carbonstate_type)           , intent(inout) :: c13_cnveg_carbonstate_inst
     type(cnveg_carbonstate_type)           , intent(inout) :: c14_cnveg_carbonstate_inst
     type(cnveg_carbonflux_type)            , intent(inout) :: c13_cnveg_carbonflux_inst
     type(cnveg_carbonflux_type)            , intent(inout) :: c14_cnveg_carbonflux_inst
     type(cnveg_state_type)                 , intent(in)    :: cnveg_state_inst
!    ! !LOCAL VARIABLES:
     integer :: fc,fp,j,i    ! indices
     integer :: p,c         !  
     real(r8),allocatable,dimension(:,:) :: vegmatrixc_old(:,:),vegmatrixc13_old(:,:),vegmatrixc14_old(:,:)
     real(r8),allocatable,dimension(:,:) :: vegmatrixc_new(:,:),vegmatrixc13_new(:,:),vegmatrixc14_new(:,:)
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
     real(r8):: epsi 
     real(r8):: days_per_year,decay_const,half_life

     real(r8):: dt        ! time step (seconds)
     real(r8):: secspyear        ! time step (seconds)
 
!	
!
    associate(                          &                                                                        
          ivt                   => patch%itype                                               , & !     
          cf13_veg              => c13_cnveg_carbonflux_inst                         , & ! In
          cf14_veg              => c14_cnveg_carbonflux_inst                         , & ! In
          cs13_veg              => c13_cnveg_carbonstate_inst                        , & ! In/Output
          cs14_veg              => c14_cnveg_carbonstate_inst                        , & ! In/Output  
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
         retransn0               => cnveg_nitrogenstate_inst%retransn0_patch                , & !
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
         matrix_C13input                  => cnveg_carbonflux_inst%matrix_C13input_patch                  , & !
         matrix_C14input                  => cnveg_carbonflux_inst%matrix_C14input_patch                  , & !

         matrix_Ninput                  => cnveg_nitrogenflux_inst%matrix_Ninput_patch                  & ! Input:      [real(r8) (:) ]    (gN/m2/s) Input N
         )
    !-----------------------------------------------------------------------
     ! set time steps
      dt = real( get_step_size(), r8 )
      days_per_year = get_days_per_year()
      secspyear = days_per_year * secspday

      allocate(vegmatrixc_old(bounds%begp:bounds%endp,nvegcpool))
      allocate(vegmatrixc_new(bounds%begp:bounds%endp,nvegcpool))
      allocate(vegmatrixn_old(bounds%begp:bounds%endp,nvegnpool))
      allocate(vegmatrixn_new(bounds%begp:bounds%endp,nvegnpool))

      allocate(vegmatrixc13_old(bounds%begp:bounds%endp,nvegcpool))
      allocate(vegmatrixc13_new(bounds%begp:bounds%endp,nvegcpool))
      allocate(vegmatrixc14_old(bounds%begp:bounds%endp,nvegcpool))
      allocate(vegmatrixc14_new(bounds%begp:bounds%endp,nvegcpool))

      allocate(vegmatrix_2d(nvegcpool,nvegcpool))
      allocate(AKinv(nvegcpool,nvegcpool))

      allocate(vegmatrixn_2d(nvegnpool,nvegnpool))
      allocate(AKinvn(nvegnpool,nvegnpool))

      epsi = 1.e-30_r8    ! small number
      counter = counter + dt
     if (counter >= 1.0 * secspyear) then ! link to the recycling span of climate forcing
          end_of_year = .true.
          counter = 0._r8
       else
          end_of_year = .false.
       end if
      half_life = 5730._r8 * secspyear
      decay_const = - log(0.5_r8) / half_life

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
         if(ivt(p) >= npcropmin)then
         vegmatrixc_old(p,igrain)        = grainc(p)
         vegmatrixc_old(p,igrain_st)     = grainc_storage(p)
         vegmatrixc_old(p,igrain_xf)     = grainc_xfer(p)
         end if
        if ( use_c13 )then
         vegmatrixc13_old(p,ileaf)         = cs13_veg%leafc_patch (p)
         vegmatrixc13_old(p,ileaf_st)      = cs13_veg%leafc_storage_patch (p)
         vegmatrixc13_old(p,ileaf_xf)      = cs13_veg%leafc_xfer_patch (p)
         vegmatrixc13_old(p,ifroot)        = cs13_veg%frootc_patch (p)
         vegmatrixc13_old(p,ifroot_st)     = cs13_veg%frootc_storage_patch (p)
         vegmatrixc13_old(p,ifroot_xf)     = cs13_veg%frootc_xfer_patch (p)
         vegmatrixc13_old(p,ilivestem)     = cs13_veg%livestemc_patch(p)
         vegmatrixc13_old(p,ilivestem_st)  = cs13_veg%livestemc_storage_patch (p)
         vegmatrixc13_old(p,ilivestem_xf)  = cs13_veg%livestemc_xfer_patch (p)
         vegmatrixc13_old(p,ideadstem)     = cs13_veg%deadstemc_patch (p)
         vegmatrixc13_old(p,ideadstem_st)  = cs13_veg%deadstemc_storage_patch (p)
         vegmatrixc13_old(p,ideadstem_xf)  = cs13_veg%deadstemc_xfer_patch (p)
         vegmatrixc13_old(p,ilivecroot)    = cs13_veg%livecrootc_patch (p)
         vegmatrixc13_old(p,ilivecroot_st) = cs13_veg%livecrootc_storage_patch (p)
         vegmatrixc13_old(p,ilivecroot_xf) = cs13_veg%livecrootc_xfer_patch (p)
         vegmatrixc13_old(p,ideadcroot)    = cs13_veg%deadcrootc_patch (p)
         vegmatrixc13_old(p,ideadcroot_st) = cs13_veg%deadcrootc_storage_patch (p)
         vegmatrixc13_old(p,ideadcroot_xf) = cs13_veg%deadcrootc_xfer_patch (p)
         if(ivt(p) >= npcropmin)then
           vegmatrixc13_old(p,igrain)        = cs13_veg%grainc_patch(p)
           vegmatrixc13_old(p,igrain_st)     = cs13_veg%grainc_storage_patch(p)
           vegmatrixc13_old(p,igrain_xf)     = cs13_veg%grainc_xfer_patch(p)
         end if
       end if 
       if ( use_c14 )then
         vegmatrixc14_old(p,ileaf)         = cs14_veg%leafc_patch(p)
         vegmatrixc14_old(p,ileaf_st)      = cs14_veg%leafc_storage_patch(p)
         vegmatrixc14_old(p,ileaf_xf)      = cs14_veg%leafc_xfer_patch(p)
         vegmatrixc14_old(p,ifroot)        = cs14_veg%frootc_patch(p)
         vegmatrixc14_old(p,ifroot_st)     = cs14_veg%frootc_storage_patch(p)
         vegmatrixc14_old(p,ifroot_xf)     = cs14_veg%frootc_xfer_patch(p)
         vegmatrixc14_old(p,ilivestem)     = cs14_veg%livestemc_patch(p)
         vegmatrixc14_old(p,ilivestem_st)  = cs14_veg%livestemc_storage_patch(p)
         vegmatrixc14_old(p,ilivestem_xf)  = cs14_veg%livestemc_xfer_patch(p)
         vegmatrixc14_old(p,ideadstem)     = cs14_veg%deadstemc_patch(p)
         vegmatrixc14_old(p,ideadstem_st)  = cs14_veg%deadstemc_storage_patch(p)
         vegmatrixc14_old(p,ideadstem_xf)  = cs14_veg%deadstemc_xfer_patch(p)
         vegmatrixc14_old(p,ilivecroot)    = cs14_veg%livecrootc_patch(p)
         vegmatrixc14_old(p,ilivecroot_st) = cs14_veg%livecrootc_storage_patch(p)
         vegmatrixc14_old(p,ilivecroot_xf) = cs14_veg%livecrootc_xfer_patch(p)
         vegmatrixc14_old(p,ideadcroot)    = cs14_veg%deadcrootc_patch(p)
         vegmatrixc14_old(p,ideadcroot_st) = cs14_veg%deadcrootc_storage_patch(p)
         vegmatrixc14_old(p,ideadcroot_xf) = cs14_veg%deadcrootc_xfer_patch(p)
         if(ivt(p) >= npcropmin)then
           vegmatrixc14_old(p,igrain)        = cs14_veg%grainc_patch(p)
           vegmatrixc14_old(p,igrain_st)     = cs14_veg%grainc_storage_patch(p)
           vegmatrixc14_old(p,igrain_xf)     = cs14_veg%grainc_xfer_patch(p)
         end if
       end if
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
         if(ivt(p) >= npcropmin)then
         vegmatrixn_old(p,igrain)        = grainn(p)
         vegmatrixn_old(p,igrain_st)     = grainn_storage(p)
         vegmatrixn_old(p,igrain_xf)     = grainn_xfer(p)
         end if
         vegmatrixn_old(p,iretransn)     = retransn(p)
         if(isspinup .or. is_outmatrix)then
           if (is_beg_curr_year() .or. is_first_step_of_this_run_segment() )then
!            if (all(vegmatrixc_old(p,:).le.epsi))then
!                 vegmatrixc_rt(:) = 0.0_r8
!            else
            leafc0(p)                = max(leafc(p),epsi)
            leafc0_storage(p)        = max(leafc_storage(p), epsi)
            leafc0_xfer(p)           = max(leafc_xfer(p), epsi)
            frootc0(p)               = max(frootc(p), epsi)
            frootc0_storage(p)       = max(frootc_storage(p), epsi)
            frootc0_xfer(p)          = max(frootc_xfer(p), epsi)
            livestemc0(p)            = max(livestemc(p), epsi)
            livestemc0_storage(p)    = max(livestemc_storage(p), epsi)
            livestemc0_xfer(p)       = max(livestemc_xfer(p), epsi)
            deadstemc0(p)            = max(deadstemc(p), epsi)
            deadstemc0_storage(p)    = max(deadstemc_storage(p), epsi)
            deadstemc0_xfer(p)       = max(deadstemc_xfer(p), epsi)
            livecrootc0(p)           = max(livecrootc(p), epsi)
            livecrootc0_storage(p)   = max(livecrootc_storage(p), epsi)
            livecrootc0_xfer(p)      = max(livecrootc_xfer(p), epsi)
            deadcrootc0(p)           = max(deadcrootc(p), epsi)
            deadcrootc0_storage(p)   = max(deadcrootc_storage(p), epsi)
            deadcrootc0_xfer(p)      = max(deadcrootc_xfer(p), epsi)
            if(ivt(p) >= npcropmin)then
            grainc0(p)               = max(grainc(p), epsi)
            grainc0_storage(p)       = max(grainc_storage(p), epsi)
            grainc0_xfer(p)          = max(grainc_xfer(p), epsi)
            end if
!           end if
!            if (all(vegmatrixn_old(p,:).le.epsi))then
!                 vegmatrixn_rt(:) = 0.0_r8
!            else
            leafn0(p)                = max(leafn(p),epsi)
            leafn0_storage(p)        = max(leafn_storage(p), epsi)
            leafn0_xfer(p)           = max(leafn_xfer(p), epsi)
            frootn0(p)               = max(frootn(p), epsi)
            frootn0_storage(p)       = max(frootn_storage(p), epsi)
            frootn0_xfer(p)          = max(frootn_xfer(p), epsi)
            livestemn0(p)            = max(livestemn(p), epsi)
            livestemn0_storage(p)    = max(livestemn_storage(p), epsi)
            livestemn0_xfer(p)       = max(livestemn_xfer(p), epsi)
            deadstemn0(p)            = max(deadstemn(p), epsi)
            deadstemn0_storage(p)    = max(deadstemn_storage(p), epsi)
            deadstemn0_xfer(p)       = max(deadstemn_xfer(p), epsi)
            livecrootn0(p)           = max(livecrootn(p), epsi)
            livecrootn0_storage(p)   = max(livecrootn_storage(p), epsi)
            livecrootn0_xfer(p)      = max(livecrootn_xfer(p), epsi)
            deadcrootn0(p)           = max(deadcrootn(p), epsi)
            deadcrootn0_storage(p)   = max(deadcrootn_storage(p), epsi)
            deadcrootn0_xfer(p)      = max(deadcrootn_xfer(p), epsi)
            retransn0(p)             = max(retransn(p), epsi)
            if(ivt(p) >= npcropmin)then
            grainn0(p)               = max(grainn(p), epsi)
            grainn0_storage(p)       = max(grainn_storage(p), epsi)
            grainn0_xfer(p)          = max(grainn_xfer(p), epsi)
            end if
!          end if 
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
          if(ivt(p) >= npcropmin)then
             vegmatrix_2d(igrain,igrain)             = grainc(p)        / grainc0(p)
             vegmatrix_2d(igrain_st,igrain_st)       = grainc_storage(p)/ grainc0_storage(p)
             vegmatrix_2d(igrain_xf,igrain_xf)       = grainc_xfer(p)   / grainc0_xfer(p)
          end if

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
          if(ivt(p) >= npcropmin)then
             vegmatrixn_2d(igrain,igrain)             = grainn(p)        / grainn0(p)
             vegmatrixn_2d(igrain_st,igrain_st)       = grainn_storage(p)/ grainn0_storage(p)
             vegmatrixn_2d(igrain_xf,igrain_xf)       = grainn_xfer(p)   / grainn0_xfer(p)
          end if
          vegmatrixn_2d(iretransn,iretransn) = retransn(p)/retransn0(p)
       end if !isspinup
!
         matrix_phturnover(p,:,:) = 0._r8
         matrix_gmturnover(p,:,:) = 0._r8
         matrix_fiturnover(p,:,:) = 0._r8
         matrix_nphturnover(p,:,:) = 0._r8
         matrix_ngmturnover(p,:,:) = 0._r8
         matrix_nfiturnover(p,:,:) = 0._r8

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

         
           vegmatrixc_new(p,:) = vegmatrixc_old(p,:)  +  matrix_alloc(p,:) * matrix_Cinput(p) * dt &
                            +(matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc_old(p,:))   &
                            + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc_old(p,:))   &
                            + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc_old(p,:))) * dt 
        if ( use_c13 ) then

           vegmatrixc13_new(p,:) = vegmatrixc13_old(p,:)  +  matrix_alloc(p,:) * matrix_C13input(p) * dt &
                            +(matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc13_old(p,:))   &
                            + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc13_old(p,:))   &
                            + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc13_old(p,:))) * dt

         cs13_veg%leafc_patch(p)               = vegmatrixc13_new(p,ileaf)
         cs13_veg%leafc_storage_patch(p)       = vegmatrixc13_new(p,ileaf_st)
         cs13_veg%leafc_xfer_patch(p)          = vegmatrixc13_new(p,ileaf_xf)
         cs13_veg%frootc_patch(p)              = vegmatrixc13_new(p,ifroot)
         cs13_veg%frootc_storage_patch(p)      = vegmatrixc13_new(p,ifroot_st)
         cs13_veg%frootc_xfer_patch(p)         = vegmatrixc13_new(p,ifroot_xf)
         cs13_veg%livestemc_patch(p)           = vegmatrixc13_new(p,ilivestem)
         cs13_veg%livestemc_storage_patch(p)   = vegmatrixc13_new(p,ilivestem_st)
         cs13_veg%livestemc_xfer_patch(p)      = vegmatrixc13_new(p,ilivestem_xf)
         cs13_veg%deadstemc_patch(p)           = vegmatrixc13_new(p,ideadstem)
         cs13_veg%deadstemc_storage_patch(p)   = vegmatrixc13_new(p,ideadstem_st)
         cs13_veg%deadstemc_xfer_patch(p)      = vegmatrixc13_new(p,ideadstem_xf)
         cs13_veg%livecrootc_patch(p)          = vegmatrixc13_new(p,ilivecroot)
         cs13_veg%livecrootc_storage_patch(p)  = vegmatrixc13_new(p,ilivecroot_st)
         cs13_veg%livecrootc_xfer_patch(p)     = vegmatrixc13_new(p,ilivecroot_xf)
         cs13_veg%deadcrootc_patch(p)          = vegmatrixc13_new(p,ideadcroot)
         cs13_veg%deadcrootc_storage_patch(p)  = vegmatrixc13_new(p,ideadcroot_st)
         cs13_veg%deadcrootc_xfer_patch(p)     = vegmatrixc13_new(p,ideadcroot_xf)
         if(ivt(p) >= npcropmin)then
           cs13_veg%grainc_patch(p)              = vegmatrixc13_new(p,igrain)
           cs13_veg%grainc_storage_patch(p)      = vegmatrixc13_new(p,igrain_st)
           cs13_veg%grainc_xfer_patch(p)         = vegmatrixc13_new(p,igrain_xf)
         end if
      end if        
     if ( use_c14 ) then
         vegmatrixc14_new(p,:) = vegmatrixc14_old(p,:)  +  matrix_alloc(p,:) * matrix_C14input(p) * dt &
                            +(matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc14_old(p,:))   &
                            + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc14_old(p,:))   &
                            + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc14_old(p,:))) * dt

         cs14_veg%leafc_patch(p)               = vegmatrixc14_new(p,ileaf)* (1._r8 - decay_const)
         cs14_veg%leafc_storage_patch(p)       = vegmatrixc14_new(p,ileaf_st)* (1._r8 - decay_const)
         cs14_veg%leafc_xfer_patch(p)          = vegmatrixc14_new(p,ileaf_xf)* (1._r8 - decay_const)
         cs14_veg%frootc_patch(p)              = vegmatrixc14_new(p,ifroot)* (1._r8 - decay_const)
         cs14_veg%frootc_storage_patch(p)      = vegmatrixc14_new(p,ifroot_st)* (1._r8 - decay_const)
         cs14_veg%frootc_xfer_patch(p)         = vegmatrixc14_new(p,ifroot_xf)* (1._r8 - decay_const)
         cs14_veg%livestemc_patch(p)           = vegmatrixc14_new(p,ilivestem)* (1._r8 - decay_const)
         cs14_veg%livestemc_storage_patch(p)   = vegmatrixc14_new(p,ilivestem_st)* (1._r8 - decay_const)
         cs14_veg%livestemc_xfer_patch(p)      = vegmatrixc14_new(p,ilivestem_xf)* (1._r8 - decay_const)
         cs14_veg%deadstemc_patch(p)           = vegmatrixc14_new(p,ideadstem)* (1._r8 - decay_const)
         cs14_veg%deadstemc_storage_patch(p)   = vegmatrixc14_new(p,ideadstem_st)* (1._r8 - decay_const)
         cs14_veg%deadstemc_xfer_patch(p)      = vegmatrixc14_new(p,ideadstem_xf)* (1._r8 - decay_const)
         cs14_veg%livecrootc_patch(p)          = vegmatrixc14_new(p,ilivecroot)* (1._r8 - decay_const)
         cs14_veg%livecrootc_storage_patch(p)  = vegmatrixc14_new(p,ilivecroot_st)* (1._r8 - decay_const)
         cs14_veg%livecrootc_xfer_patch(p)     = vegmatrixc14_new(p,ilivecroot_xf)* (1._r8 - decay_const)
         cs14_veg%deadcrootc_patch(p)          = vegmatrixc14_new(p,ideadcroot)* (1._r8 - decay_const)
         cs14_veg%deadcrootc_storage_patch(p)  = vegmatrixc14_new(p,ideadcroot_st)* (1._r8 - decay_const)
         cs14_veg%deadcrootc_xfer_patch(p)     = vegmatrixc14_new(p,ideadcroot_xf)* (1._r8 - decay_const)
         if(ivt(p) >= npcropmin)then
           cs14_veg%grainc_patch(p)              = vegmatrixc14_new(p,igrain)* (1._r8 - decay_const)
           cs14_veg%grainc_storage_patch(p)      = vegmatrixc14_new(p,igrain_st)* (1._r8 - decay_const)
           cs14_veg%grainc_xfer_patch(p)         = vegmatrixc14_new(p,igrain_xf)* (1._r8 - decay_const)
         end if
     end if 

           vegmatrixn_new(p,:) = vegmatrixn_old(p,:)  +  matrix_nalloc(p,:) * matrix_Ninput(p) * dt &
                            + (matmul(matmul(matrix_nphtransfer(p,1:nvegnpool,:),matrix_nphturnover(p,:,:)),vegmatrixn_old(p,:)) &
                            +  matmul(matmul(matrix_ngmtransfer(p,1:nvegnpool,:),matrix_ngmturnover(p,:,:)),vegmatrixn_old(p,:)) &
                            +  matmul(matmul(matrix_nfitransfer(p,1:nvegnpool,:),matrix_nfiturnover(p,:,:)),vegmatrixn_old(p,:))) * dt 
         if(isspinup .or. is_outmatrix)then !
            matrix_alloc_acc(p,:) = matrix_alloc_acc(p,:) + matrix_alloc(p,:) * matrix_Cinput(p) * dt
            matrix_transfer_acc(p,1:nvegcpool,:) = matrix_transfer_acc(p,1:nvegcpool,:) &
                               + (matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrix_2d(:,:))& 
                               + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrix_2d(:,:)) &
                               + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrix_2d(:,:)))*dt

            matrix_nalloc_acc(p,:) = matrix_nalloc_acc(p,:) + matrix_nalloc(p,:) * matrix_Ninput(p)* dt
            matrix_ntransfer_acc(p,1:nvegnpool,:) = matrix_ntransfer_acc(p,1:nvegnpool,:) &
                               + (matmul(matmul(matrix_nphtransfer(p,1:nvegnpool,:),matrix_nphturnover(p,:,:)),vegmatrixn_2d(:,:))& 
                               + matmul(matmul(matrix_ngmtransfer(p,1:nvegnpool,:),matrix_ngmturnover(p,:,:)),vegmatrixn_2d(:,:)) &
                               + matmul(matmul(matrix_nfitransfer(p,1:nvegnpool,:),matrix_nfiturnover(p,:,:)),vegmatrixn_2d(:,:)))*dt
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
         if(ivt(p) >= npcropmin)then
         grainc(p)              = vegmatrixc_new(p,igrain)
         grainc_storage(p)      = vegmatrixc_new(p,igrain_st)
         grainc_xfer(p)         = vegmatrixc_new(p,igrain_xf)
         end if
 
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
         if(ivt(p) >= npcropmin)then
            grainn(p)              = vegmatrixn_new(p,igrain)
            grainn_storage(p)      = vegmatrixn_new(p,igrain_st)
            grainn_xfer(p)         = vegmatrixn_new(p,igrain_xf)
         end if
         retransn(p)               = vegmatrixn_new(p,iretransn)
         if(isspinup .or. is_outmatrix)then
            if(end_of_year)then
               do i=1,nvegcpool
                 if(abs(matrix_transfer_acc(p,i,i)) .le. epsi)then
                    matrix_transfer_acc(p,i,i) = 1.e+36
                 end if
               end do
               do i=1,nvegnpool
                 if(abs(matrix_ntransfer_acc(p,i,i)) .le. epsi)then
                    matrix_ntransfer_acc(p,i,i) = 1.e+36
                 end if
               end do
!              end if
               call inverse(matrix_transfer_acc(p,1:nvegcpool,:),AKinv(1:nvegcpool,1:nvegcpool),nvegcpool)
               vegmatrixc_rt(:) = -matmul(AKinv(1:nvegcpool,1:nvegcpool),matrix_alloc_acc(p,:))
! N 
               call inverse(matrix_ntransfer_acc(p,1:nvegnpool,:),AKinvn(1:nvegnpool,1:nvegnpool),nvegnpool)
               vegmatrixn_rt(:) = -matmul(AKinvn(1:nvegnpool,1:nvegnpool),matrix_nalloc_acc(p,:))
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
!                  grainc_storage(p)          = vegmatrixc_rt(igrain_st)
!                  grainc_xfer(p)             = vegmatrixc_rt(igrain_xf)
               end if
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
!                  grainn_storage(p)          = vegmatrixn_rt(igrain_st)
!                  grainn_xfer(p)             = vegmatrixn_rt(igrain_xf)
               end if
!                retransn(p)                = vegmatrixn_rt(iretransn)
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
!
               matrix_alloc_acc(p,:) = 0._r8
               matrix_transfer_acc(p,:,:) = 0._r8
               matrix_nalloc_acc(p,:) = 0._r8
               matrix_ntransfer_acc(p,:,:) = 0._r8
! C       
           end if
         end if

         matrix_Cinput(p) = 0._r8
         matrix_C13input(p) = 0._r8
         matrix_C14input(p) = 0._r8       
         matrix_phtransfer(p,:,:) = 0._r8
         matrix_gmtransfer(p,:,:) = 0._r8
         matrix_fitransfer(p,:,:) = 0._r8
         matrix_nphtransfer(p,:,:) = 0._r8
         matrix_ngmtransfer(p,:,:) = 0._r8
         matrix_nfitransfer(p,:,:) = 0._r8
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
