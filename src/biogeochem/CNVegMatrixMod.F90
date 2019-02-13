module CNVegMatrixMod

  !----------------------------------------------------------------------------------
  ! The matrix model of CLM5.0 was developed by Yiqi Luo EcoLab members, 
  ! Drs. Xingjie Lu, Yuanyuan Huang and Zhengguang Du, at Northern Arizona University
  !----------------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  ! Matrix solution for vegetation C and N cycles
  ! The matrix equation 
  ! Xn+1 = Xn + I*dt + (A*ksi*k)*Xn*dt
  
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use clm_time_manager               , only : get_step_size,is_end_curr_year,is_first_step_of_this_run_segment,&
                                              get_days_per_year,is_beg_curr_year,is_end_curr_year
  use decompMod                      , only : bounds_type 
  use clm_varpar                     , only : nlevdecomp, nvegcpool, nvegnpool
  use clm_varpar                     , only : ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                              ilivestem,ilivestem_st,ilivestem_xf,&
                                              ideadstem,ideadstem_st,ideadstem_xf,&
                                              ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                              ideadcroot,ideadcroot_st,ideadcroot_xf,&
                                              igrain,igrain_st,igrain_xf,iretransn,ioutc,ioutn,&
                                              ncphtrans,nnphtrans,ncgmtrans,nngmtrans,ncfitrans,nnfitrans,&
                                              ncphouttrans,nnphouttrans,ncgmouttrans,nngmouttrans,ncfiouttrans,nnfiouttrans
  use perf_mod                       , only : t_startf, t_stopf
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
  use clm_varctl                     , only : use_c13, use_c14 
  use SPMMod                         , only : sparse_matrix_type,diag_matrix_type,vector_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNVegMatrix
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
   subroutine CNVegMatrix(bounds,num_soilp,filter_soilp,num_actfirep,filter_actfirep,cnveg_carbonstate_inst,cnveg_nitrogenstate_inst,&
                          cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,cnveg_state_inst,soilbiogeochem_nitrogenflux_inst, &
                          c13_cnveg_carbonstate_inst,c14_cnveg_carbonstate_inst,c13_cnveg_carbonflux_inst,&
                          c14_cnveg_carbonflux_inst)
    ! !DESCRIPTION:
    ! !ARGUMENTS:
     type(bounds_type)                      , intent(in)    :: bounds
     integer                                , intent(in)    :: num_soilp       ! number of soil patches in filter
     integer                                , intent(in)    :: filter_soilp(:) ! filter for soil patches
     integer                                , intent(in)    :: num_actfirep       ! number of soil patches in filter
     integer                                , intent(in)    :: filter_actfirep(:) ! filter for soil patches
     type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst  
     type(cnveg_nitrogenstate_type)         , intent(inout) :: cnveg_nitrogenstate_inst
     type(cnveg_carbonflux_type)            , intent(inout) :: cnveg_carbonflux_inst
     type(cnveg_nitrogenflux_type)          , intent(inout) :: cnveg_nitrogenflux_inst
     type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
     type(cnveg_carbonstate_type)           , intent(inout) :: c13_cnveg_carbonstate_inst
     type(cnveg_carbonstate_type)           , intent(inout) :: c14_cnveg_carbonstate_inst
     type(cnveg_carbonflux_type)            , intent(inout) :: c13_cnveg_carbonflux_inst
     type(cnveg_carbonflux_type)            , intent(inout) :: c14_cnveg_carbonflux_inst
     type(cnveg_state_type)                 , intent(in)    :: cnveg_state_inst
!    ! !LOCAL VARIABLES:
!     type(sparse_matrix_type),save         :: Avegc
!     type(sparse_matrix_type),save         :: Avegn
!     type(sparse_matrix_type),save         :: AKphvegc
!     type(sparse_matrix_type),save         :: AKphvegn
!     type(sparse_matrix_type),save         :: AKgmvegc
!     type(sparse_matrix_type),save         :: AKgmvegn
!     type(sparse_matrix_type),save         :: AKfivegc
!     type(sparse_matrix_type),save         :: AKfivegn
!     type(sparse_matrix_type),save         :: AKtmp1vegc
!     type(sparse_matrix_type),save         :: AKtmp2vegc
!     type(sparse_matrix_type),save         :: AKtmp1vegn
!     type(sparse_matrix_type),save         :: AKtmp2vegn
!     logical,save                          :: SMInitialized = .False.
!     type(diag_matrix_type)                :: Kvegc
!     type(diag_matrix_type)                :: Kvegn
!     type(vector_type)                     :: Xoldvegc
!     type(vector_type)                     :: Xnewvegc
!     type(vector_type)                     :: Xoldvegn
!     type(vector_type)                     :: Xnewvegn
!     type(vector_type)                     :: Xoldveg13c
!     type(vector_type)                     :: Xnewveg13c
!     type(vector_type)                     :: Xoldveg14c
!     type(vector_type)                     :: Xnewveg14c
     integer :: fc,fp,j,i,k    ! indices
     integer :: p,c         !  
     real(r8) tmptmp
     real(r8),dimension(:,:)     :: Aphconed(bounds%begp:bounds%endp,ncphtrans-ncphouttrans)
     real(r8),dimension(:,:)     :: Aphnoned(bounds%begp:bounds%endp,nnphtrans-nnphouttrans)
     real(r8),dimension(:,:)     :: Agmconed(bounds%begp:bounds%endp,ncgmtrans-ncgmouttrans)
     real(r8),dimension(:,:)     :: Agmnoned(bounds%begp:bounds%endp,nngmtrans-nngmouttrans)
     real(r8),dimension(:,:)     :: Aficoned(bounds%begp:bounds%endp,ncfitrans-ncfiouttrans)
     real(r8),dimension(:,:)     :: Afinoned(bounds%begp:bounds%endp,nnfitrans-nnfiouttrans)
!     real(r8),allocatable,dimension(:,:)     :: Aphconed(bounds%begp:bounds%endp,ncphtrans-ncphouttrans)
!     real(r8),allocatable,dimension(:,:)     :: Aphnoned(bounds%begp:bounds%endp,nnphtrans-nnphouttrans)
!     real(r8),allocatable,dimension(:,:)     :: Agmconed(bounds%begp:bounds%endp,ncgmtrans-ncgmouttrans)
!     real(r8),allocatable,dimension(:,:)     :: Agmnoned(bounds%begp:bounds%endp,nngmtrans-nngmouttrans)
!     real(r8),allocatable,dimension(:,:)     :: Aficoned(bounds%begp:bounds%endp,ncfitrans-ncfiouttrans)
!     real(r8),allocatable,dimension(:,:)     :: Afinoned(bounds%begp:bounds%endp,nnfitrans-nnfiouttrans)

     integer,dimension(:)      :: AI_phc(ncphtrans-ncphouttrans)
     integer,dimension(:)      :: AI_phn(nnphtrans-nnphouttrans)
     integer,dimension(:)      :: AI_gmc(ncgmtrans-ncgmouttrans)
     integer,dimension(:)      :: AI_gmn(nngmtrans-nngmouttrans)
     integer,dimension(:)      :: AI_fic(ncfitrans-ncfiouttrans)
     integer,dimension(:)      :: AI_fin(nnfitrans-nnfiouttrans)

!     integer,allocatable,dimension(:)      :: AI_phc
!     integer,allocatable,dimension(:)      :: AI_phn
!     integer,allocatable,dimension(:)      :: AI_gmc
!     integer,allocatable,dimension(:)      :: AI_gmn
!     integer,allocatable,dimension(:)      :: AI_fic
!     integer,allocatable,dimension(:)      :: AI_fin

     integer,dimension(:)      :: AJ_phc(ncphtrans-ncphouttrans)
     integer,dimension(:)      :: AJ_phn(nnphtrans-nnphouttrans)
     integer,dimension(:)      :: AJ_gmc(ncgmtrans-ncgmouttrans)
     integer,dimension(:)      :: AJ_gmn(nngmtrans-nngmouttrans)
     integer,dimension(:)      :: AJ_fic(ncfitrans-ncfiouttrans)
     integer,dimension(:)      :: AJ_fin(nnfitrans-nnfiouttrans)

!     integer,allocatable,dimension(:)      :: AJ_phc
!     integer,allocatable,dimension(:)      :: AJ_phn
!     integer,allocatable,dimension(:)      :: AJ_gmc
!     integer,allocatable,dimension(:)      :: AJ_gmn
!     integer,allocatable,dimension(:)      :: AJ_fic
!     integer,allocatable,dimension(:)      :: AJ_fin
     logical,save                          :: list_ready = .false.
!     logical,save                          :: index_ready = .false.

     real(r8),dimension(1:nvegcpool)       :: Kconed
     real(r8),dimension(1:nvegnpool)       :: Knoned

!     type(vector_type)                     :: vegmatrixc_old
!     type(vector_type)                     :: vegmatrixc13_old
!     type(vector_type)                     :: vegmatrixc14_old
!     type(vector_type)                     :: vegmatrixc_new
!     type(vector_type)                     :: vegmatrixc13_new
!     type(vector_type)                     :: vegmatrixc14_new
!     type(vector_type)                     :: vegmatrixn_old       
!     type(vector_type)                     :: vegmatrixn_new       
     type(vector_type)                     :: vegmatrixc_input    
     type(vector_type)                     :: vegmatrixc13_input  
     type(vector_type)                     :: vegmatrixc14_input  
     type(vector_type)                     :: vegmatrixn_input    
!     type(vector_type)                     :: matrix_calloc_acc   
!     type(vector_type)                     :: matrix_nalloc_acc  
     logical, save                         :: init_ready_aphc      = .false.
     logical, save                         :: init_ready_agmc      = .false.
     logical, save                         :: init_ready_afic      = .false.
     logical, save                         :: init_ready_aphn      = .false.
     logical, save                         :: init_ready_agmn      = .false.
     logical, save                         :: init_ready_afin      = .false.
     logical, save                         :: list_ready_phgmfic   = .false.
     logical, save                         :: list_ready_phgmc     = .false.
     logical, save                         :: list_ready_phgmfin   = .false.
     logical, save                         :: list_ready_phgmn     = .false.
!     integer, save, dimension(:)           :: RI_phc(1:ncphtrans-ncphouttrans+nvegcpool) = -9999
!     integer, save, dimension(:)           :: CI_phc(1:ncphtrans-ncphouttrans+nvegcpool) = -9999
!     integer, save, dimension(:)           :: RI_gmc(1:ncgmtrans-ncgmouttrans+nvegcpool) = -9999
!     integer, save, dimension(:)           :: CI_gmc(1:ncgmtrans-ncgmouttrans+nvegcpool) = -9999
!     integer, save, dimension(:)           :: RI_fic(1:ncfitrans-ncfiouttrans+nvegcpool) = -9999
!     integer, save, dimension(:)           :: CI_fic(1:ncfitrans-ncfiouttrans+nvegcpool) = -9999
!     integer, save, dimension(:)           :: RI_phn(1:nnphtrans-nnphouttrans+nvegnpool) = -9999
!     integer, save, dimension(:)           :: CI_phn(1:nnphtrans-nnphouttrans+nvegnpool) = -9999
!     integer, save, dimension(:)           :: RI_gmn(1:nngmtrans-nngmouttrans+nvegnpool) = -9999
!     integer, save, dimension(:)           :: CI_gmn(1:nngmtrans-nngmouttrans+nvegnpool) = -9999
!     integer, save, dimension(:)           :: RI_fin(1:nnfitrans-nnfiouttrans+nvegnpool) = -9999
!     integer, save, dimension(:)           :: CI_fin(1:nnfitrans-nnfiouttrans+nvegnpool) = -9999

!     integer,save,dimension(:)             :: list_aphc(1:ncphtrans-ncphouttrans)
!     integer,save,dimension(:)             :: list_agmc(1:ncgmtrans-ncgmouttrans)
!     integer,save,dimension(:)             :: list_afic(1:ncfitrans-ncfiouttrans)
!     integer,save,dimension(:)             :: list_aphn(1:nnphtrans-nnphouttrans)
!     integer,save,dimension(:)             :: list_agmn(1:nngmtrans-nngmouttrans)
!     integer,save,dimension(:)             :: list_afin(1:nnfitrans-nnfiouttrans)

!     real(r8),allocatable,dimension(:)     :: vegmatrixc_old
!     real(r8),allocatable,dimension(:)     :: vegmatrixc13_old
!     real(r8),allocatable,dimension(:)     :: vegmatrixc14_old
!     real(r8),allocatable,dimension(:)     :: vegmatrixc_new
!     real(r8),allocatable,dimension(:)     :: vegmatrixc13_new
!     real(r8),allocatable,dimension(:)     :: vegmatrixc14_new
!     real(r8),allocatable,dimension(:)     :: vegmatrixn_old       (:)
!     real(r8),allocatable,dimension(:)     :: vegmatrixn_new       (:)
!     real(r8),allocatable,dimension(:)     :: vegmatrixc_input     (:)
!     real(r8),allocatable,dimension(:)     :: vegmatrixc13_input   (:)
!     real(r8),allocatable,dimension(:)     :: vegmatrixc14_input   (:)
     real(r8),dimension(:,:)   :: vegmatrixc_transfer  (1:nvegcpool,1:nvegcpool)
     real(r8),dimension(:,:)   :: vegmatrixc13_transfer  (1:nvegcpool,1:nvegcpool)
     real(r8),dimension(:,:)   :: vegmatrixc14_transfer  (1:nvegcpool,1:nvegcpool)
!     real(r8),allocatable,dimension(:)     :: vegmatrixn_input     (:)
     real(r8),dimension(:,:)   :: vegmatrixn_transfer  (1:nvegnpool,1:nvegnpool)
     real(r8),dimension(:)     :: matrix_calloc_acc    (1:nvegcpool)
     real(r8),dimension(:)     :: matrix_nalloc_acc    (1:nvegnpool)
     real(r8),dimension(:,:)   :: matrix_ctransfer_acc (1:nvegcpool,1:nvegcpool)
     real(r8),dimension(:,:)   :: matrix_ntransfer_acc (1:nvegnpool,1:nvegnpool)


!     real(r8),allocatable,dimension(:,:) :: vegmatrixc_rt(:,:),vegmatrixn_rt(:,:)
!     real(r8),allocatable,dimension(:,:,:) :: matrix_nphturnover(:,:,:),matrix_ngmturnover(:,:,:)!,matrix_nphturnover(:,:,:)

! for spinupacc
     real(r8),dimension(1:nvegcpool) :: vegmatrixc_rt,rowonec
     real(r8),dimension(1:nvegnpool) :: vegmatrixn_rt,rowonen,tmp
     real(r8),dimension(1:nvegcpool,1:nvegcpool) :: A_C,K_C,AK_C
     real(r8),dimension(1:nvegnpool,1:nvegnpool) :: A_N,K_N,AK_N
     real(r8),dimension(:,:) :: AKinvc(1:nvegcpool,1:nvegcpool),AKinvn(1:nvegnpool,1:nvegnpool)
     real(r8):: epsi 
     real(r8):: days_per_year,decay_const,half_life
     
     real(r8):: dt        ! time step (seconds)
     real(r8):: secspyear        ! time step (seconds)
     integer :: yr,mon,day,sec
!KO     real(r8),dimension(1:3*nvegcpool) :: tempdump
 
!	
!
    associate(                          &                                                                        
          ivt                   => patch%itype                                               , & ! Input:  [integer  (:) ]  patch vegetation type
          cf13_veg              => c13_cnveg_carbonflux_inst                         , & ! In
          cf14_veg              => c14_cnveg_carbonflux_inst                         , & ! In
          cs13_veg              => c13_cnveg_carbonstate_inst                        , & ! In/Output
          cs14_veg              => c14_cnveg_carbonstate_inst                        , & ! In/Output  
          retransn              => cnveg_nitrogenstate_inst%retransn_patch           , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N
  !
         fire_closs             => cnveg_carbonflux_inst%fire_closs_patch            , &
 
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
!         matrix_pot_leafc                 => cnveg_carbonstate_inst%matrix_pot_leafc_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for matrix calculation                                    
!         matrix_pot_leafc_storage         => cnveg_carbonstate_inst%matrix_pot_leafc_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage C for matrix calculation                                   
!         matrix_pot_leafc_xfer            => cnveg_carbonstate_inst%matrix_pot_leafc_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer C for matrix calcuation                                   
!         matrix_pot_frootc                => cnveg_carbonstate_inst%matrix_pot_frootc_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C for matrix calculation
!         matrix_pot_frootc_storage        => cnveg_carbonstate_inst%matrix_pot_frootc_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C for matrix calculation
!         matrix_pot_frootc_xfer           => cnveg_carbonstate_inst%matrix_pot_frootc_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C for matrix calculation
!         matrix_pot_livestemc             => cnveg_carbonstate_inst%matrix_pot_livestemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem C for matrix calculationleaf C
!         matrix_pot_livestemc_storage     => cnveg_carbonstate_inst%matrix_pot_livestemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage C for matrix calculation
!         matrix_pot_livestemc_xfer        => cnveg_carbonstate_inst%matrix_pot_livestemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer C for matrix calculation
!         matrix_pot_deadstemc             => cnveg_carbonstate_inst%matrix_pot_deadstemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem C for matrix calculationleaf C
!         matrix_pot_deadstemc_storage     => cnveg_carbonstate_inst%matrix_pot_deadstemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage C for matrix calculation
!         matrix_pot_deadstemc_xfer        => cnveg_carbonstate_inst%matrix_pot_deadstemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer C for matrix calculation                                    
!         matrix_pot_livecrootc            => cnveg_carbonstate_inst%matrix_pot_livecrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root C for matrix calculationleaf C 
!         matrix_pot_livecrootc_storage    => cnveg_carbonstate_inst%matrix_pot_livecrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage C for matrix calculation
!         matrix_pot_livecrootc_xfer       => cnveg_carbonstate_inst%matrix_pot_livecrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer C for matrix calculation
!         matrix_pot_deadcrootc            => cnveg_carbonstate_inst%matrix_pot_deadcrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root C for matrix calculationleaf C
!         matrix_pot_deadcrootc_storage    => cnveg_carbonstate_inst%matrix_pot_deadcrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage C for matrix calculation
!         matrix_pot_deadcrootc_xfer       => cnveg_carbonstate_inst%matrix_pot_deadcrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer C for matrix calculation
!         matrix_pot_grainc                => cnveg_carbonstate_inst%matrix_pot_grainc_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C for matrix calculation
!         matrix_pot_grainc_storage        => cnveg_carbonstate_inst%matrix_pot_grainc_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C for matrix calculation
!         matrix_pot_grainc_xfer           => cnveg_carbonstate_inst%matrix_pot_grainc_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C for matrix calculation
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
!         matrix_pot_leafn                 => cnveg_nitrogenstate_inst%matrix_pot_leafn_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf N for matrix calculation                                    
!         matrix_pot_leafn_storage         => cnveg_nitrogenstate_inst%matrix_pot_leafn_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage N for matrix calculation                                   
!         matrix_pot_leafn_xfer            => cnveg_nitrogenstate_inst%matrix_pot_leafn_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer N for  matrix calcuation                                   
!         matrix_pot_frootn                => cnveg_nitrogenstate_inst%matrix_pot_frootn_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
!         matrix_pot_frootn_storage        => cnveg_nitrogenstate_inst%matrix_pot_frootn_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
!         matrix_pot_frootn_xfer           => cnveg_nitrogenstate_inst%matrix_pot_frootn_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
!         matrix_pot_livestemn             => cnveg_nitrogenstate_inst%matrix_pot_livestemn_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem N for matrix calculationleaf N
!         matrix_pot_livestemn_storage     => cnveg_nitrogenstate_inst%matrix_pot_livestemn_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage N for matrix calculation
!         matrix_pot_livestemn_xfer        => cnveg_nitrogenstate_inst%matrix_pot_livestemn_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer N for matrix calculation
!         matrix_pot_deadstemn             => cnveg_nitrogenstate_inst%matrix_pot_deadstemn_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem N for matrix calculationleaf N
!         matrix_pot_deadstemn_storage     => cnveg_nitrogenstate_inst%matrix_pot_deadstemn_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage N for matrix calculation
!         matrix_pot_deadstemn_xfer        => cnveg_nitrogenstate_inst%matrix_pot_deadstemn_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer N for matrix calculation                                    
!         matrix_pot_livecrootn            => cnveg_nitrogenstate_inst%matrix_pot_livecrootn_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root N for matrix calculationleaf N 
!         matrix_pot_livecrootn_storage    => cnveg_nitrogenstate_inst%matrix_pot_livecrootn_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage N for matrix calculation
!         matrix_pot_livecrootn_xfer       => cnveg_nitrogenstate_inst%matrix_pot_livecrootn_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer N for matrix calculation
!         matrix_pot_deadcrootn            => cnveg_nitrogenstate_inst%matrix_pot_deadcrootn_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root N for matrix calculationleaf N
!         matrix_pot_deadcrootn_storage    => cnveg_nitrogenstate_inst%matrix_pot_deadcrootn_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage N for matrix calculation
!         matrix_pot_deadcrootn_xfer       => cnveg_nitrogenstate_inst%matrix_pot_deadcrootn_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer N for matrix calculation
!         matrix_pot_grainn                => cnveg_nitrogenstate_inst%matrix_pot_grainn_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root N for matrix calculation
!         matrix_pot_grainn_storage        => cnveg_nitrogenstate_inst%matrix_pot_grainn_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage N for matrix calculation
!         matrix_pot_grainn_xfer           => cnveg_nitrogenstate_inst%matrix_pot_grainn_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer N for matrix calculation
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
         matrix_ctransfer_fire_livestem_to_deadstem_acc   => cnveg_carbonstate_inst%matrix_ctransfer_fire_livestem_to_deadstem_acc_patch          , & !
         matrix_ctransfer_fire_livecroot_to_deadcroot_acc => cnveg_carbonstate_inst%matrix_ctransfer_fire_livecroot_to_deadcroot_acc_patch          , & !

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
!         matrix_cturnover_gm_leaf_acc                     => cnveg_carbonstate_inst%matrix_cturnover_gm_leaf_acc_patch               , &
!         matrix_cturnover_gm_leafst_acc                   => cnveg_carbonstate_inst%matrix_cturnover_gm_leafst_acc_patch               , &
!         matrix_cturnover_gm_leafxf_acc                   => cnveg_carbonstate_inst%matrix_cturnover_gm_leafxf_acc_patch               , &
!         matrix_cturnover_gm_froot_acc                    => cnveg_carbonstate_inst%matrix_cturnover_gm_froot_acc_patch               , &
!         matrix_cturnover_gm_frootst_acc                  => cnveg_carbonstate_inst%matrix_cturnover_gm_frootst_acc_patch               , &
!         matrix_cturnover_gm_frootxf_acc                  => cnveg_carbonstate_inst%matrix_cturnover_gm_frootxf_acc_patch               , &
!         matrix_cturnover_gm_livestem_acc                 => cnveg_carbonstate_inst%matrix_cturnover_gm_livestem_acc_patch               , &
!         matrix_cturnover_gm_livestemst_acc               => cnveg_carbonstate_inst%matrix_cturnover_gm_livestemst_acc_patch               , &
!         matrix_cturnover_gm_livestemxf_acc               => cnveg_carbonstate_inst%matrix_cturnover_gm_livestemxf_acc_patch               , &
!         matrix_cturnover_gm_deadstem_acc                 => cnveg_carbonstate_inst%matrix_cturnover_gm_deadstem_acc_patch               , &
!         matrix_cturnover_gm_deadstemst_acc               => cnveg_carbonstate_inst%matrix_cturnover_gm_deadstemst_acc_patch               , &
!         matrix_cturnover_gm_deadstemxf_acc               => cnveg_carbonstate_inst%matrix_cturnover_gm_deadstemxf_acc_patch               , &
!         matrix_cturnover_gm_livecroot_acc                => cnveg_carbonstate_inst%matrix_cturnover_gm_livecroot_acc_patch               , &
!         matrix_cturnover_gm_livecrootst_acc              => cnveg_carbonstate_inst%matrix_cturnover_gm_livecrootst_acc_patch               , &
!         matrix_cturnover_gm_livecrootxf_acc              => cnveg_carbonstate_inst%matrix_cturnover_gm_livecrootxf_acc_patch               , &
!         matrix_cturnover_gm_deadcroot_acc                => cnveg_carbonstate_inst%matrix_cturnover_gm_deadcroot_acc_patch               , &
!         matrix_cturnover_gm_deadcrootst_acc              => cnveg_carbonstate_inst%matrix_cturnover_gm_deadcrootst_acc_patch               , &
!         matrix_cturnover_gm_deadcrootxf_acc              => cnveg_carbonstate_inst%matrix_cturnover_gm_deadcrootxf_acc_patch               , &
!         matrix_cturnover_fire_leaf_acc                   => cnveg_carbonstate_inst%matrix_cturnover_fire_leaf_acc_patch               , &
!         matrix_cturnover_fire_leafst_acc                 => cnveg_carbonstate_inst%matrix_cturnover_fire_leafst_acc_patch               , &
!         matrix_cturnover_fire_leafxf_acc                 => cnveg_carbonstate_inst%matrix_cturnover_fire_leafxf_acc_patch               , &
!         matrix_cturnover_fire_froot_acc                  => cnveg_carbonstate_inst%matrix_cturnover_fire_froot_acc_patch               , &
!         matrix_cturnover_fire_frootst_acc                => cnveg_carbonstate_inst%matrix_cturnover_fire_frootst_acc_patch               , &
!         matrix_cturnover_fire_frootxf_acc                => cnveg_carbonstate_inst%matrix_cturnover_fire_frootxf_acc_patch               , &
!         matrix_cturnover_fire_livestem_acc               => cnveg_carbonstate_inst%matrix_cturnover_fire_livestem_acc_patch               , &
!         matrix_cturnover_fire_livestemst_acc             => cnveg_carbonstate_inst%matrix_cturnover_fire_livestemst_acc_patch               , &
!         matrix_cturnover_fire_livestemxf_acc             => cnveg_carbonstate_inst%matrix_cturnover_fire_livestemxf_acc_patch               , &
!         matrix_cturnover_fire_deadstem_acc               => cnveg_carbonstate_inst%matrix_cturnover_fire_deadstem_acc_patch               , &
!         matrix_cturnover_fire_deadstemst_acc             => cnveg_carbonstate_inst%matrix_cturnover_fire_deadstemst_acc_patch               , &
!         matrix_cturnover_fire_deadstemxf_acc             => cnveg_carbonstate_inst%matrix_cturnover_fire_deadstemxf_acc_patch               , &
!         matrix_cturnover_fire_livecroot_acc              => cnveg_carbonstate_inst%matrix_cturnover_fire_livecroot_acc_patch               , &
!         matrix_cturnover_fire_livecrootst_acc            => cnveg_carbonstate_inst%matrix_cturnover_fire_livecrootst_acc_patch               , &
!         matrix_cturnover_fire_livecrootxf_acc            => cnveg_carbonstate_inst%matrix_cturnover_fire_livecrootxf_acc_patch               , &
!         matrix_cturnover_fire_deadcroot_acc              => cnveg_carbonstate_inst%matrix_cturnover_fire_deadcroot_acc_patch               , &
!         matrix_cturnover_fire_deadcrootst_acc            => cnveg_carbonstate_inst%matrix_cturnover_fire_deadcrootst_acc_patch               , &
!         matrix_cturnover_fire_deadcrootxf_acc            => cnveg_carbonstate_inst%matrix_cturnover_fire_deadcrootxf_acc_patch               , &

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
         matrix_C13input                  => cnveg_carbonflux_inst%matrix_C13input_patch                  , & !
         matrix_C14input                  => cnveg_carbonflux_inst%matrix_C14input_patch                  , & !

         matrix_Ninput                  => cnveg_nitrogenflux_inst%matrix_Ninput_patch        ,          & ! Input:      [real(r8) (:) ]    (gN/m2/s) Input N
         doner_phc                      => cnveg_carbonflux_inst%matrix_phtransfer_doner_patch        , &
         receiver_phc                   => cnveg_carbonflux_inst%matrix_phtransfer_receiver_patch     , &
         doner_gmc                      => cnveg_carbonflux_inst%matrix_gmtransfer_doner_patch        , &
         receiver_gmc                   => cnveg_carbonflux_inst%matrix_gmtransfer_receiver_patch     , &
         doner_fic                      => cnveg_carbonflux_inst%matrix_fitransfer_doner_patch        , &
         receiver_fic                   => cnveg_carbonflux_inst%matrix_fitransfer_receiver_patch     , &
         doner_phn                      => cnveg_nitrogenflux_inst%matrix_nphtransfer_doner_patch        , &
         receiver_phn                   => cnveg_nitrogenflux_inst%matrix_nphtransfer_receiver_patch     , &
         doner_gmn                      => cnveg_nitrogenflux_inst%matrix_ngmtransfer_doner_patch        , &
         receiver_gmn                   => cnveg_nitrogenflux_inst%matrix_ngmtransfer_receiver_patch     , &
         doner_fin                      => cnveg_nitrogenflux_inst%matrix_nfitransfer_doner_patch        , &
         receiver_fin                   => cnveg_nitrogenflux_inst%matrix_nfitransfer_receiver_patch     , &
         ileafst_to_ileafxf_phc  => cnveg_carbonflux_inst%ileafst_to_ileafxf_ph,         &
         ileafxf_to_ileaf_phc    => cnveg_carbonflux_inst%ileafxf_to_ileaf_ph,         &
         ifrootst_to_ifrootxf_phc  => cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph,         &
         ifrootxf_to_ifroot_phc    => cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph,         &
         ilivestemst_to_ilivestemxf_phc  => cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph,         &
         ilivestemxf_to_ilivestem_phc    => cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph,         &
         ideadstemst_to_ideadstemxf_phc  => cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph,         &
         ideadstemxf_to_ideadstem_phc    => cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph,         &
         ilivecrootst_to_ilivecrootxf_phc  => cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph,         &
         ilivecrootxf_to_ilivecroot_phc    => cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph,         &
         ideadcrootst_to_ideadcrootxf_phc  => cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph,         &
         ideadcrootxf_to_ideadcroot_phc    => cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph,         &
         ilivestem_to_ideadstem_phc  => cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph,         &
         ilivecroot_to_ideadcroot_phc  => cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph,         &
         ileaf_to_iout_phc  => cnveg_carbonflux_inst%ileaf_to_iout_ph,         &
         ifroot_to_iout_phc  => cnveg_carbonflux_inst%ifroot_to_iout_ph,         &
         ilivestem_to_iout_phc  => cnveg_carbonflux_inst%ilivestem_to_iout_ph,         &
         igrain_to_iout_phc   =>  cnveg_carbonflux_inst%igrain_to_iout_ph,         &
         ileaf_to_iout_gmc => cnveg_carbonflux_inst%ileaf_to_iout_gm , &
         ileafst_to_iout_gmc => cnveg_carbonflux_inst%ileafst_to_iout_gm , &
         ileafxf_to_iout_gmc => cnveg_carbonflux_inst%ileafxf_to_iout_gm , &
         ifroot_to_iout_gmc => cnveg_carbonflux_inst%ifroot_to_iout_gm , &
         ifrootst_to_iout_gmc => cnveg_carbonflux_inst%ifrootst_to_iout_gm , &
         ifrootxf_to_iout_gmc => cnveg_carbonflux_inst%ifrootxf_to_iout_gm , &
         ilivestem_to_iout_gmc => cnveg_carbonflux_inst%ilivestem_to_iout_gm , &
         ilivestemst_to_iout_gmc => cnveg_carbonflux_inst%ilivestemst_to_iout_gm , &
         ilivestemxf_to_iout_gmc => cnveg_carbonflux_inst%ilivestemxf_to_iout_gm , &
         ideadstem_to_iout_gmc => cnveg_carbonflux_inst%ideadstem_to_iout_gm , &
         ideadstemst_to_iout_gmc => cnveg_carbonflux_inst%ideadstemst_to_iout_gm , &
         ideadstemxf_to_iout_gmc => cnveg_carbonflux_inst%ideadstemxf_to_iout_gm , &
         ilivecroot_to_iout_gmc => cnveg_carbonflux_inst%ilivecroot_to_iout_gm , &
         ilivecrootst_to_iout_gmc => cnveg_carbonflux_inst%ilivecrootst_to_iout_gm , &
         ilivecrootxf_to_iout_gmc => cnveg_carbonflux_inst%ilivecrootxf_to_iout_gm , &
         ideadcroot_to_iout_gmc => cnveg_carbonflux_inst%ideadcroot_to_iout_gm , &
         ideadcrootst_to_iout_gmc => cnveg_carbonflux_inst%ideadcrootst_to_iout_gm , &
         ideadcrootxf_to_iout_gmc => cnveg_carbonflux_inst%ideadcrootxf_to_iout_gm , &
         ileaf_to_iout_fic                   => cnveg_carbonflux_inst%ileaf_to_iout_fi ,   &
         ileafst_to_iout_fic                 => cnveg_carbonflux_inst%ileafst_to_iout_fi,   &
         ileafxf_to_iout_fic                 => cnveg_carbonflux_inst%ileafxf_to_iout_fi,   &
         ifroot_to_iout_fic                  => cnveg_carbonflux_inst%ifroot_to_iout_fi,   &
         ifrootst_to_iout_fic                => cnveg_carbonflux_inst%ifrootst_to_iout_fi,   &
         ifrootxf_to_iout_fic                => cnveg_carbonflux_inst%ifrootxf_to_iout_fi,   &
         ilivestem_to_iout_fic               => cnveg_carbonflux_inst%ilivestem_to_iout_fi,   &
         ilivestemst_to_iout_fic             => cnveg_carbonflux_inst%ilivestemst_to_iout_fi,   &
         ilivestemxf_to_iout_fic             => cnveg_carbonflux_inst%ilivestemxf_to_iout_fi,   &
         ideadstem_to_iout_fic               => cnveg_carbonflux_inst%ideadstem_to_iout_fi,   &
         ideadstemst_to_iout_fic             => cnveg_carbonflux_inst%ideadstemst_to_iout_fi,   &
         ideadstemxf_to_iout_fic             => cnveg_carbonflux_inst%ideadstemxf_to_iout_fi,   &
         ilivecroot_to_iout_fic              => cnveg_carbonflux_inst%ilivecroot_to_iout_fi,   &
         ilivecrootst_to_iout_fic            => cnveg_carbonflux_inst%ilivecrootst_to_iout_fi,   &
         ilivecrootxf_to_iout_fic            => cnveg_carbonflux_inst%ilivecrootxf_to_iout_fi,   &
         ideadcroot_to_iout_fic              => cnveg_carbonflux_inst%ideadcroot_to_iout_fi,   &
         ideadcrootst_to_iout_fic            => cnveg_carbonflux_inst%ideadcrootst_to_iout_fi,   &
         ideadcrootxf_to_iout_fic            => cnveg_carbonflux_inst%ideadcrootxf_to_iout_fi,   &
         ilivestem_to_ideadstem_fic          => cnveg_carbonflux_inst%ilivestem_to_ideadstem_fi,   &
         ilivecroot_to_ideadcroot_fic        => cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_fi,   &

         ileafst_to_ileafxf_phn  => cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph,         &
         ileafxf_to_ileaf_phn    => cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph,         &
         ifrootst_to_ifrootxf_phn  => cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph,         &
         ifrootxf_to_ifroot_phn    => cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph,         &
         ilivestemst_to_ilivestemxf_phn  => cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph,         &
         ilivestemxf_to_ilivestem_phn    => cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph,         &
         ideadstemst_to_ideadstemxf_phn  => cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph,         &
         ideadstemxf_to_ideadstem_phn    => cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph,         &
         ilivecrootst_to_ilivecrootxf_phn  => cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph,         &
         ilivecrootxf_to_ilivecroot_phn    => cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph,         &
         ideadcrootst_to_ideadcrootxf_phn  => cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph,         &
         ideadcrootxf_to_ideadcroot_phn    => cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph,         &
         ilivestem_to_ideadstem_phn  => cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph,         &
         ilivecroot_to_ideadcroot_phn  => cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph,         &
         ileaf_to_iout_phn  => cnveg_nitrogenflux_inst%ileaf_to_iout_ph,         &
         ifroot_to_iout_phn  => cnveg_nitrogenflux_inst%ifroot_to_iout_ph,         &
         ilivestem_to_iout_phn  => cnveg_nitrogenflux_inst%ilivestem_to_iout_ph,         &
         iretransn_to_iout_phn    => cnveg_nitrogenflux_inst%iretransn_to_iout_ph             , &
         igrain_to_iout_phn   =>  cnveg_nitrogenflux_inst%igrain_to_iout_ph        , &
         ileaf_to_iretransn_phn  => cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph,         &
         ifroot_to_iretransn_phn      => cnveg_nitrogenflux_inst%ifroot_to_iretransn_ph                 , &
         ilivestem_to_iretransn_phn  => cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph,         &
         ilivecroot_to_iretransn_phn  => cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph,         &
         iretransn_to_ileaf_phn       => cnveg_nitrogenflux_inst%iretransn_to_ileaf_ph             , &
         iretransn_to_ileafst_phn     => cnveg_nitrogenflux_inst%iretransn_to_ileafst_ph             , &
         iretransn_to_ifroot_phn      => cnveg_nitrogenflux_inst%iretransn_to_ifroot_ph             , &
         iretransn_to_ifrootst_phn    => cnveg_nitrogenflux_inst%iretransn_to_ifrootst_ph             , &
         iretransn_to_ilivestem_phn   => cnveg_nitrogenflux_inst%iretransn_to_ilivestem_ph             , &
         iretransn_to_ilivestemst_phn => cnveg_nitrogenflux_inst%iretransn_to_ilivestemst_ph             , &
         iretransn_to_ideadstem_phn   => cnveg_nitrogenflux_inst%iretransn_to_ideadstem_ph             , &
         iretransn_to_ideadstemst_phn => cnveg_nitrogenflux_inst%iretransn_to_ideadstemst_ph             , &
         iretransn_to_ilivecroot_phn  => cnveg_nitrogenflux_inst%iretransn_to_ilivecroot_ph             , &
         iretransn_to_ilivecrootst_phn=> cnveg_nitrogenflux_inst%iretransn_to_ilivecrootst_ph             , &
         iretransn_to_ideadcroot_phn  => cnveg_nitrogenflux_inst%iretransn_to_ideadcroot_ph             , &
         iretransn_to_ideadcrootst_phn=> cnveg_nitrogenflux_inst%iretransn_to_ideadcrootst_ph             , &
         iretransn_to_igrain_phn      => cnveg_nitrogenflux_inst%iretransn_to_igrain_ph             , &
         iretransn_to_igrainst_phn    => cnveg_nitrogenflux_inst%iretransn_to_igrainst_ph       ,      &
         ileaf_to_iout_gmn => cnveg_nitrogenflux_inst%ileaf_to_iout_gm , &
         ileafst_to_iout_gmn => cnveg_nitrogenflux_inst%ileafst_to_iout_gm , &
         ileafxf_to_iout_gmn => cnveg_nitrogenflux_inst%ileafxf_to_iout_gm , &
         ifroot_to_iout_gmn => cnveg_nitrogenflux_inst%ifroot_to_iout_gm , &
         ifrootst_to_iout_gmn => cnveg_nitrogenflux_inst%ifrootst_to_iout_gm , &
         ifrootxf_to_iout_gmn => cnveg_nitrogenflux_inst%ifrootxf_to_iout_gm , &
         ilivestem_to_iout_gmn => cnveg_nitrogenflux_inst%ilivestem_to_iout_gm , &
         ilivestemst_to_iout_gmn => cnveg_nitrogenflux_inst%ilivestemst_to_iout_gm , &
         ilivestemxf_to_iout_gmn => cnveg_nitrogenflux_inst%ilivestemxf_to_iout_gm , &
         ideadstem_to_iout_gmn => cnveg_nitrogenflux_inst%ideadstem_to_iout_gm , &
         ideadstemst_to_iout_gmn => cnveg_nitrogenflux_inst%ideadstemst_to_iout_gm , &
         ideadstemxf_to_iout_gmn => cnveg_nitrogenflux_inst%ideadstemxf_to_iout_gm , &
         ilivecroot_to_iout_gmn => cnveg_nitrogenflux_inst%ilivecroot_to_iout_gm , &
         ilivecrootst_to_iout_gmn => cnveg_nitrogenflux_inst%ilivecrootst_to_iout_gm , &
         ilivecrootxf_to_iout_gmn => cnveg_nitrogenflux_inst%ilivecrootxf_to_iout_gm , &
         ideadcroot_to_iout_gmn => cnveg_nitrogenflux_inst%ideadcroot_to_iout_gm , &
         ideadcrootst_to_iout_gmn => cnveg_nitrogenflux_inst%ideadcrootst_to_iout_gm , &
         ideadcrootxf_to_iout_gmn => cnveg_nitrogenflux_inst%ideadcrootxf_to_iout_gm , &
         iretransn_to_iout_gmn => cnveg_nitrogenflux_inst%iretransn_to_iout_gm ,&
         ileaf_to_iout_fin                   => cnveg_nitrogenflux_inst%ileaf_to_iout_fi,   &
         ileafst_to_iout_fin                 => cnveg_nitrogenflux_inst%ileafst_to_iout_fi,   &
         ileafxf_to_iout_fin                 => cnveg_nitrogenflux_inst%ileafxf_to_iout_fi,   &
         ifroot_to_iout_fin                  => cnveg_nitrogenflux_inst%ifroot_to_iout_fi,   &
         ifrootst_to_iout_fin                => cnveg_nitrogenflux_inst%ifrootst_to_iout_fi,   &
         ifrootxf_to_iout_fin                => cnveg_nitrogenflux_inst%ifrootxf_to_iout_fi,   &
         ilivestem_to_iout_fin               => cnveg_nitrogenflux_inst%ilivestem_to_iout_fi,   &
         ilivestemst_to_iout_fin             => cnveg_nitrogenflux_inst%ilivestemst_to_iout_fi,   &
         ilivestemxf_to_iout_fin             => cnveg_nitrogenflux_inst%ilivestemxf_to_iout_fi,   &
         ideadstem_to_iout_fin               => cnveg_nitrogenflux_inst%ideadstem_to_iout_fi,   &
         ideadstemst_to_iout_fin             => cnveg_nitrogenflux_inst%ideadstemst_to_iout_fi,   &
         ideadstemxf_to_iout_fin             => cnveg_nitrogenflux_inst%ideadstemxf_to_iout_fi,   &
         ilivecroot_to_iout_fin              => cnveg_nitrogenflux_inst%ilivecroot_to_iout_fi,   &
         ilivecrootst_to_iout_fin            => cnveg_nitrogenflux_inst%ilivecrootst_to_iout_fi,   &
         ilivecrootxf_to_iout_fin            => cnveg_nitrogenflux_inst%ilivecrootxf_to_iout_fi,   &
         ideadcroot_to_iout_fin              => cnveg_nitrogenflux_inst%ideadcroot_to_iout_fi,   &
         ideadcrootst_to_iout_fin            => cnveg_nitrogenflux_inst%ideadcrootst_to_iout_fi,   &
         ideadcrootxf_to_iout_fin            => cnveg_nitrogenflux_inst%ideadcrootxf_to_iout_fi,   &
         ilivestem_to_ideadstem_fin          => cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_fi,   &
         ilivecroot_to_ideadcroot_fin        => cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_fi,   &
         iretransn_to_iout_fin               => cnveg_nitrogenflux_inst%iretransn_to_iout_fi,   &
         AKphvegc                            => cnveg_carbonflux_inst%AKphvegc,  &
         AKgmvegc                            => cnveg_carbonflux_inst%AKgmvegc,  &
         AKfivegc                            => cnveg_carbonflux_inst%AKfivegc,  &
!         AKphvegc                            => cnveg_carbonflux_inst%AKphvegc,  &
!         AKgmvegc                            => cnveg_carbonflux_inst%AKgmvegc,  &
!         AKfivegc                            => cnveg_carbonflux_inst%AKfivegc,  &
!         AKtmp1vegc                          => cnveg_carbonflux_inst%AKtmp1vegc,  &
         AKallvegc                          => cnveg_carbonflux_inst%AKallvegc,  &
         NE_AKphgmc                         => cnveg_carbonflux_inst%NE_AKphgmc,  &
         NE_AKallvegc                         => cnveg_carbonflux_inst%NE_AKallvegc,  &
         RI_AKphgmc                         => cnveg_carbonflux_inst%RI_AKphgmc,  &
         RI_AKallvegc                         => cnveg_carbonflux_inst%RI_AKallvegc,  &
         CI_AKphgmc                         => cnveg_carbonflux_inst%CI_AKphgmc,  &
         CI_AKallvegc                         => cnveg_carbonflux_inst%CI_AKallvegc,  &
         Kvegc                               => cnveg_carbonflux_inst%Kvegc,  &
!         Xoldvegc                            => cnveg_carbonflux_inst%Xoldvegc,  &
         Xvegc                            => cnveg_carbonflux_inst%Xvegc,  &
!         Xfirelossvegc                    => cnveg_carbonflux_inst%Xfirelossvegc,  &
!         Aphvegn                             => cnveg_nitrogenflux_inst%Aphvegn,  &
!         Agmvegn                             => cnveg_nitrogenflux_inst%Agmvegn,  &
!         Afivegn                             => cnveg_nitrogenflux_inst%Afivegn,  &
         AKphvegn                            => cnveg_nitrogenflux_inst%AKphvegn,  &
         AKgmvegn                            => cnveg_nitrogenflux_inst%AKgmvegn,  &
         AKfivegn                            => cnveg_nitrogenflux_inst%AKfivegn,  &
!         AKtmp1vegn                          => cnveg_nitrogenflux_inst%AKtmp1vegn,  &
         AKallvegn                          => cnveg_nitrogenflux_inst%AKallvegn,  &
         NE_AKphgmn                         => cnveg_nitrogenflux_inst%NE_AKphgmn,  &
         NE_AKallvegn                         => cnveg_nitrogenflux_inst%NE_AKallvegn,  &
         RI_AKphgmn                         => cnveg_nitrogenflux_inst%RI_AKphgmn,  &
         CI_AKphgmn                         => cnveg_nitrogenflux_inst%CI_AKphgmn,  &
         RI_AKallvegn                         => cnveg_nitrogenflux_inst%RI_AKallvegn,  &
         CI_AKallvegn                         => cnveg_nitrogenflux_inst%CI_AKallvegn,  &
         RI_phc                               => cnveg_carbonflux_inst%RI_phc, &
         CI_phc                               => cnveg_carbonflux_inst%CI_phc, &
         RI_gmc                               => cnveg_carbonflux_inst%RI_gmc, &
         CI_gmc                               => cnveg_carbonflux_inst%CI_gmc, &
         RI_fic                               => cnveg_carbonflux_inst%RI_fic, &
         CI_fic                               => cnveg_carbonflux_inst%CI_fic, &
         RI_phn                               => cnveg_nitrogenflux_inst%RI_phn, &
         CI_phn                               => cnveg_nitrogenflux_inst%CI_phn, &
         RI_gmn                               => cnveg_nitrogenflux_inst%RI_gmn, &
         CI_gmn                               => cnveg_nitrogenflux_inst%CI_gmn, &
         RI_fin                               => cnveg_nitrogenflux_inst%RI_fin, &
         CI_fin                               => cnveg_nitrogenflux_inst%CI_fin, &
         Kvegn                               => cnveg_nitrogenflux_inst%Kvegn,  &
         Xvegn                            => cnveg_nitrogenflux_inst%Xvegn,  &
!         Xnewvegn                            => cnveg_nitrogenflux_inst%Xnewvegn,  &
         Xveg13c                          => cnveg_carbonflux_inst%Xveg13c,  &
!         Xnewveg13c                          => cnveg_carbonflux_inst%Xnewveg13c,  &
         Xveg14c                          => cnveg_carbonflux_inst%Xveg14c,  &
!         Xnewveg14c                          => cnveg_carbonflux_inst%Xnewveg14c,  &
         list_aphc                           => cnveg_carbonflux_inst%list_aphc,  &
         list_agmc                           => cnveg_carbonflux_inst%list_agmc,  &
         list_afic                           => cnveg_carbonflux_inst%list_afic,  &
         list_aphn                           => cnveg_nitrogenflux_inst%list_aphn,  &
         list_agmn                           => cnveg_nitrogenflux_inst%list_agmn,  &
         list_afin                           => cnveg_nitrogenflux_inst%list_afin,  &
         list_phc_phgm                       => cnveg_carbonflux_inst%list_phc_phgmc, &
         list_gmc_phgm                       => cnveg_carbonflux_inst%list_gmc_phgmc, &
         list_phc_phgmfi                     => cnveg_carbonflux_inst%list_phc_phgmfic, &
         list_gmc_phgmfi                     => cnveg_carbonflux_inst%list_gmc_phgmfic, &
         list_fic_phgmfi                     => cnveg_carbonflux_inst%list_fic_phgmfic, &
         list_phn_phgm                       => cnveg_nitrogenflux_inst%list_phn_phgmn, &
         list_gmn_phgm                       => cnveg_nitrogenflux_inst%list_gmn_phgmn, &
         list_phn_phgmfi                     => cnveg_nitrogenflux_inst%list_phn_phgmfin, &
         list_gmn_phgmfi                     => cnveg_nitrogenflux_inst%list_gmn_phgmfin, &
         list_fin_phgmfi                     => cnveg_nitrogenflux_inst%list_fin_phgmfin &
         )
    !-----------------------------------------------------------------------
     ! set time steps
     !print*,'here0'
      call t_startf('CN veg matrix-init')
      dt = real( get_step_size(), r8 )
      secspyear = get_days_per_year() * secspday
     

!      if(.not. SMInitialized)then
!         call Avegc%InitSM(nvegcpool)
!         call AKphvegc%InitSM(nvegcpool)
!         call AKgmvegc%InitSM(nvegcpool)
!         call AKfivegc%InitSM(nvegcpool)
!         call AKtmp1vegc%InitSM(nvegcpool)
!         call AKtmp2vegc%InitSM(nvegcpool)
!         call Avegn%InitSM(nvegnpool)
!         call AKphvegn%InitSM(nvegnpool)
!         call AKgmvegn%InitSM(nvegnpool)
!         call AKfivegn%InitSM(nvegnpool)
!         call AKtmp1vegn%InitSM(nvegnpool)
!         call AKtmp2vegn%InitSM(nvegnpool)
!         call Kvegc%InitDM(nvegcpool)
!         call Kvegn%InitDM(nvegnpool)
!         call Xoldvegc%InitV(nvegcpool)
!         call Xoldvegn%InitV(nvegnpool)
!         call Xnewvegc%InitV(nvegcpool)
!         call Xnewvegn%InitV(nvegnpool)
!         if(use_c13)then
!            call Xoldveg13c%InitV(nvegcpool)
!            call Xnewveg13c%InitV(nvegcpool)
!         end if
!         if(use_c14)then
!            call Xoldveg14c%InitV(nvegcpool)
!            call Xnewveg14c%InitV(nvegcpool)
!         end if
!         SMInitialized = .True.
!      end if
      call t_stopf('CN veg matrix-init')
      call t_startf('CN veg matrix-alloc')

!      if(ncphtrans-ncphouttrans .gt. 0)then
!         allocate(Aphconed(bounds%begp:bounds%endp,ncphtrans-ncphouttrans))
!         allocate(AI_phc  (ncphtrans-ncphouttrans))
!         allocate(AJ_phc  (ncphtrans-ncphouttrans))
!      end if
!      if(nnphtrans-nnphouttrans .gt. 0)then
!         allocate(Aphnoned(bounds%begp:bounds%endp,nnphtrans-nnphouttrans))
!         allocate(AI_phn  (nnphtrans-nnphouttrans))
!         allocate(AJ_phn  (nnphtrans-nnphouttrans))
!      end if
!      if(ncgmtrans-ncgmouttrans .gt. 0)then
!         allocate(Agmconed(bounds%begp:bounds%endp,ncgmtrans-ncgmouttrans))
!         allocate(AI_gmc  (ncgmtrans-ncgmouttrans))
!         allocate(AJ_gmc  (ncgmtrans-ncgmouttrans))
!      end if
!      if(nngmtrans-nngmouttrans .gt. 0)then
!         allocate(Agmnoned(bounds%begp:bounds%endp,nngmtrans-nngmouttrans))
!         allocate(AI_gmn  (nngmtrans-nngmouttrans))
!         allocate(AJ_gmn  (nngmtrans-nngmouttrans))
!      end if
!      if(ncfitrans-ncfiouttrans .gt. 0)then
!         allocate(Aficoned(bounds%begp:bounds%endp,ncfitrans-ncfiouttrans))
!         allocate(AI_fic  (ncfitrans-ncfiouttrans))
!         allocate(AJ_fic  (ncfitrans-ncfiouttrans))
!      end if
!      if(nnfitrans-nnfiouttrans .gt. 0)then
!         allocate(Afinoned(bounds%begp:bounds%endp,nnfitrans-nnfiouttrans))
!         allocate(AI_fin  (nnfitrans-nnfiouttrans))
!         allocate(AJ_fin  (nnfitrans-nnfiouttrans))
!      end if
!      allocate(vegmatrixc_old       (nvegcpool)) !save it as two dimensional variables in order to track transfers
!      allocate(vegmatrixc_new       (nvegcpool))
!      allocate(vegmatrixc13_old     (nvegcpool)) !save it as two dimensional variables in order to track transfers
!      allocate(vegmatrixc13_new     (nvegcpool))
!      allocate(vegmatrixc14_old     (nvegcpool)) !save it as two dimensional variables in order to track transfers
!      allocate(vegmatrixc14_new     (nvegcpool))
!      allocate(vegmatrixn_old       (nvegnpool))
!      allocate(vegmatrixn_new       (nvegnpool))
!      allocate(vegmatrixc_input     (nvegcpool))
      call vegmatrixc_input%InitV(nvegcpool,bounds%begp,bounds%endp)
      call vegmatrixc13_input%InitV(nvegcpool,bounds%begp,bounds%endp)
      call vegmatrixc14_input%InitV(nvegcpool,bounds%begp,bounds%endp)
      call vegmatrixn_input%InitV(nvegnpool,bounds%begp,bounds%endp)
!      call matrix_calloc_acc%InitV(nvegcpool,bounds%begp,bounds%endp)
!      call matrix_nalloc_acc%InitV(nvegnpool,bounds%begp,bounds%endp)

!      allocate(vegmatrixc13_input   (nvegcpool))
!      allocate(vegmatrixc14_input   (nvegcpool))
!      allocate(vegmatrixn_input     (nvegnpool))
!      allocate(vegmatrixc_transfer  (nvegcpool,nvegcpool))
!      allocate(vegmatrixc13_transfer(nvegcpool,nvegcpool))
!      allocate(vegmatrixc14_transfer(nvegcpool,nvegcpool))
!      allocate(vegmatrixn_transfer  (nvegnpool,nvegnpool))
!      allocate(matrix_calloc_acc    (nvegcpool))
!      allocate(matrix_nalloc_acc    (nvegnpool))
!      allocate(matrix_ctransfer_acc (nvegcpool,nvegcpool))
!      allocate(matrix_ntransfer_acc (nvegnpool,nvegnpool))
      

!      allocate(vegmatrix_2d(nvegcpool,nvegcpool))
!      allocate(AKinvc(nvegcpool,nvegcpool))

!      allocate(vegmatrixn_2d(nvegnpool,nvegnpool))
!      allocate(AKinvn(nvegnpool,nvegnpool))

!      Xoldvegc%V       (:)   = 0._r8
!      if(use_c13)then
!         Xoldveg13c%V    (:)   = 0._r8
!      end if
!      if(use_c14)then
!         Xoldveg14c%V     (:)   = 0._r8
!      end if
!      Xoldvegn%V       (:)   = 0._r8
!      Xnewvegc%V       (:)     = 0._r8
!      vegmatrixn_new       (:)     = 0._r8
!      vegmatrixc_input     (:,:)     = 0._r8
!      vegmatrixc13_input   (:)     = 0._r8
!      vegmatrixc14_input   (:)     = 0._r8
      vegmatrixc_transfer  (:,:)   = 0._r8
      vegmatrixc13_transfer(:,:)   = 0._r8
      vegmatrixc14_transfer(:,:)   = 0._r8
!      vegmatrixn_input   (:)       = 0._r8
      vegmatrixn_transfer(:,:)     = 0._r8
      matrix_calloc_acc    (:)     = 0._r8
      matrix_nalloc_acc    (:)     = 0._r8
      matrix_ctransfer_acc (:,:)   = 0._r8
      matrix_ntransfer_acc (:,:)   = 0._r8
      AKinvc (:,:) = 0._r8
      AKinvn (:,:) = 0._r8
      rowonec(:) = 1._r8
      rowonen(:) = 1._r8
      
      epsi = 1.e-30_r8     ! small number   
      half_life = 5730._r8 * secspyear
      decay_const = - log(0.5_r8) / half_life
      
      call t_stopf('CN veg matrix-alloc')

      call t_startf('CN veg matrix-assigning matrix')

      do i = 1,nvegcpool
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            matrix_phturnover(p,i) = 0._r8
            matrix_gmturnover(p,i) = 0._r8
            matrix_fiturnover(p,i) = 0._r8
            matrix_nphturnover(p,i) = 0._r8
            matrix_ngmturnover(p,i) = 0._r8
            matrix_nfiturnover(p,i) = 0._r8
         end do
      end do

      do k=1,ncphtrans
         do fp = 1,num_soilp
            p = filter_soilp(fp)
!            if(p .eq. 24 .and. ((receiver_phc(k) .eq. ileaf) .or. (doner_phc(k) .eq. ileaf)))print*,'matrix_phtransfer',p,k,doner_phc(k),receiver_phc(k),matrix_phtransfer(p,k) * dt
            matrix_phturnover(p,doner_phc(k)) = matrix_phturnover(p,doner_phc(k)) + matrix_phtransfer(p,k) * dt
!            if(p .eq. 10580 .and. ((doner_phc(k) .eq. ileaf) .or. (doner_phc(k) .eq. ifroot)))print*,'matrix_phturnover',doner_phc(k),receiver_phc(k),matrix_phturnover(p,doner_phc(k)),matrix_phtransfer(p,k) * dt
         end do
      end do

      do k=1,ncgmtrans
         do fp = 1,num_soilp
            p = filter_soilp(fp)
!            if(p .eq. 24 .and. ((receiver_gmc(k) .eq. ileaf) .or. (doner_gmc(k) .eq. ileaf)))print*,'matrix_gmtransfer',p,k,doner_gmc(k),receiver_gmc(k),matrix_gmtransfer(p,k) * dt
            matrix_gmturnover(p,doner_gmc(k)) = matrix_gmturnover(p,doner_gmc(k)) + matrix_gmtransfer(p,k) * dt
         end do
      end do

      do k=1,ncfitrans
         do fp = 1,num_soilp
            p = filter_soilp(fp)
!            if(p .eq. 24 .and. ((receiver_fic(k) .eq. ileaf) .or. (doner_fic(k) .eq. ileaf)))print*,'matrix_fitransfer',p,k,doner_fic(k),receiver_fic(k),matrix_fitransfer(p,k) * dt
            matrix_fiturnover(p,doner_fic(k)) = matrix_fiturnover(p,doner_fic(k)) + matrix_fitransfer(p,k)* dt
         end do
      end do
      
      do k=1,nnphtrans
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            matrix_nphturnover(p,doner_phn(k)) = matrix_nphturnover(p,doner_phn(k)) + matrix_nphtransfer(p,k)* dt
         end do
      end do


      do k=1,nngmtrans
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            matrix_ngmturnover(p,doner_gmn(k)) = matrix_ngmturnover(p,doner_gmn(k)) + matrix_ngmtransfer(p,k)* dt
         end do
      end do


      do k=1,nnfitrans
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            matrix_nfiturnover(p,doner_fin(k)) = matrix_nfiturnover(p,doner_fin(k)) + matrix_nfitransfer(p,k)* dt
         end do
      end do

      
      if(ncphtrans .gt. ncphouttrans)then
         do k=1,ncphtrans-ncphouttrans
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(matrix_phturnover(p,doner_phc(k)) .ne. 0)then
                  Aphconed(p,k) = matrix_phtransfer(p,k) * dt / matrix_phturnover(p,doner_phc(k))
               else
                  Aphconed(p,k) = 0._r8
               end if
            end do
         end do
      end if

      if(ncgmtrans .gt. ncgmouttrans)then
         do k=1,ncgmtrans-ncgmouttrans
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(matrix_gmturnover(p,doner_gmc(k)) .ne. 0)then
                  Agmconed(p,k) = matrix_gmtransfer(p,k) * dt / matrix_gmturnover(p,doner_gmc(k))
               else
                  Agmconed(p,k) = 0._r8
               end if
            end do
         end do
      end if

      if(ncfitrans .gt. ncfiouttrans)then
         do k=1,ncfitrans-ncfiouttrans
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(matrix_fiturnover(p,doner_fic(k)) .ne. 0)then
                  Aficoned(p,k) = matrix_fitransfer(p,k) * dt / matrix_fiturnover(p,doner_fic(k))
               else
                  Aficoned(p,k) = 0._r8
               end if
            end do
         end do
      end if

      if(nnphtrans .gt. nnphouttrans)then
         do k=1,nnphtrans-nnphouttrans
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(matrix_nphturnover(p,doner_phn(k)) .ne. 0)then
                  Aphnoned(p,k) = matrix_nphtransfer(p,k) * dt / matrix_nphturnover(p,doner_phn(k))
               else
                  Aphnoned(p,k) = 0._r8
               end if
            end do
         end do
      end if

      if(nngmtrans .gt. nngmouttrans)then
         do k=1,nngmtrans-nngmouttrans
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(matrix_ngmturnover(p,doner_phn(k)) .ne. 0)then
                  Agmnoned(p,k) = matrix_ngmtransfer(p,k) * dt / matrix_ngmturnover(p,doner_phn(k))
               else
                  Agmnoned(p,k) = 0._r8
               end if
            end do
         end do
      end if

      if(nnfitrans .gt. nnfiouttrans)then
         do k=1,nnfitrans-nnfiouttrans
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(matrix_nfiturnover(p,doner_fin(k)) .ne. 0)then
                  Afinoned(p,k) = matrix_nfitransfer(p,k) * dt / matrix_nfiturnover(p,doner_fin(k))
               else
                  Afinoned(p,k) = 0._r8
               end if
            end do
         end do
      end if

!N
      call t_stopf('CN veg matrix-assigning matrix')

      call t_startf('CN veg matrix-set old value')
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         Xvegc%V(p,ileaf)         = leafc(p)
         Xvegc%V(p,ileaf_st)      = leafc_storage(p)
         Xvegc%V(p,ileaf_xf)      = leafc_xfer(p)
         Xvegc%V(p,ifroot)        = frootc(p)
         Xvegc%V(p,ifroot_st)     = frootc_storage(p)
         Xvegc%V(p,ifroot_xf)     = frootc_xfer(p)
         Xvegc%V(p,ilivestem)     = livestemc(p)
         Xvegc%V(p,ilivestem_st)  = livestemc_storage(p)
         Xvegc%V(p,ilivestem_xf)  = livestemc_xfer(p)
         Xvegc%V(p,ideadstem)     = deadstemc(p)
         Xvegc%V(p,ideadstem_st)  = deadstemc_storage(p)
         Xvegc%V(p,ideadstem_xf)  = deadstemc_xfer(p)
         Xvegc%V(p,ilivecroot)    = livecrootc(p)
         Xvegc%V(p,ilivecroot_st) = livecrootc_storage(p)
         Xvegc%V(p,ilivecroot_xf) = livecrootc_xfer(p)
         Xvegc%V(p,ideadcroot)    = deadcrootc(p)
         Xvegc%V(p,ideadcroot_st) = deadcrootc_storage(p)
         Xvegc%V(p,ideadcroot_xf) = deadcrootc_xfer(p)
!         Xfirelossvegc%V(p,ileaf)         = leafc(p)
!         Xfirelossvegc%V(p,ileaf_st)      = leafc_storage(p)
!         Xfirelossvegc%V(p,ileaf_xf)      = leafc_xfer(p)
!         Xfirelossvegc%V(p,ifroot)        = frootc(p)
!         Xfirelossvegc%V(p,ifroot_st)     = frootc_storage(p)
!         Xfirelossvegc%V(p,ifroot_xf)     = frootc_xfer(p)
!         Xfirelossvegc%V(p,ilivestem)     = livestemc(p)
!         Xfirelossvegc%V(p,ilivestem_st)  = livestemc_storage(p)
!         Xfirelossvegc%V(p,ilivestem_xf)  = livestemc_xfer(p)
!         Xfirelossvegc%V(p,ideadstem)     = deadstemc(p)
!         Xfirelossvegc%V(p,ideadstem_st)  = deadstemc_storage(p)
!         Xfirelossvegc%V(p,ideadstem_xf)  = deadstemc_xfer(p)
!         Xfirelossvegc%V(p,ilivecroot)    = livecrootc(p)
!         Xfirelossvegc%V(p,ilivecroot_st) = livecrootc_storage(p)
!         Xfirelossvegc%V(p,ilivecroot_xf) = livecrootc_xfer(p)
!         Xfirelossvegc%V(p,ideadcroot)    = deadcrootc(p)
!         Xfirelossvegc%V(p,ideadcroot_st) = deadcrootc_storage(p)
!         Xfirelossvegc%V(p,ideadcroot_xf) = deadcrootc_xfer(p)
!         if(p .eq. 6851)print*,'begin of CNVegmatrix',Xvegc%V(p,ifroot_st),frootc_storage(p),ifroot_st
      end do

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         if(ivt(p) >= npcropmin)then
            Xvegc%V(p,igrain)     = grainc(p)
            Xvegc%V(p,igrain_st)  = grainc_storage(p)
            Xvegc%V(p,igrain_xf)  = grainc_xfer(p)
         end if
      end do
            do fp = 1,num_soilp
               p = filter_soilp(fp)
!               if(p .eq. 6851)print*,'after ini-1',p,Xvegc%V(p,1:nvegcpool)
            end do
      if ( use_c13 )then
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            Xveg13c%V(p,ileaf)         = cs13_veg%leafc_patch(p)
            Xveg13c%V(p,ileaf_st)      = cs13_veg%leafc_storage_patch(p)
            Xveg13c%V(p,ileaf_xf)      = cs13_veg%leafc_xfer_patch(p)
            Xveg13c%V(p,ifroot)        = cs13_veg%frootc_patch(p)
            Xveg13c%V(p,ifroot_st)     = cs13_veg%frootc_storage_patch(p)
            Xveg13c%V(p,ifroot_xf)     = cs13_veg%frootc_xfer_patch(p)
            Xveg13c%V(p,ilivestem)     = cs13_veg%livestemc_patch(p)
            Xveg13c%V(p,ilivestem_st)  = cs13_veg%livestemc_storage_patch(p)
            Xveg13c%V(p,ilivestem_xf)  = cs13_veg%livestemc_xfer_patch(p)
            Xveg13c%V(p,ideadstem)     = cs13_veg%deadstemc_patch(p)
            Xveg13c%V(p,ideadstem_st)  = cs13_veg%deadstemc_storage_patch(p)
            Xveg13c%V(p,ideadstem_xf)  = cs13_veg%deadstemc_xfer_patch(p)
            Xveg13c%V(p,ilivecroot)    = cs13_veg%livecrootc_patch(p)
            Xveg13c%V(p,ilivecroot_st) = cs13_veg%livecrootc_storage_patch(p)
            Xveg13c%V(p,ilivecroot_xf) = cs13_veg%livecrootc_xfer_patch(p)
            Xveg13c%V(p,ideadcroot)    = cs13_veg%deadcrootc_patch(p)
            Xveg13c%V(p,ideadcroot_st) = cs13_veg%deadcrootc_storage_patch(p)
            Xveg13c%V(p,ideadcroot_xf) = cs13_veg%deadcrootc_xfer_patch(p)
         end do

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            if(ivt(p) >= npcropmin)then
               Xveg13c%V(p,igrain)     = cs13_veg%grainc_patch(p)
               Xveg13c%V(p,igrain_st)  = cs13_veg%grainc_storage_patch(p)
               Xveg13c%V(p,igrain_xf)  = cs13_veg%grainc_xfer_patch(p)
            end if
         end do
      end if

      if ( use_c14 )then
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            Xveg14c%V(p,ileaf)         = cs14_veg%leafc_patch(p)
            Xveg14c%V(p,ileaf_st)      = cs14_veg%leafc_storage_patch(p)
            Xveg14c%V(p,ileaf_xf)      = cs14_veg%leafc_xfer_patch(p)
            Xveg14c%V(p,ifroot)        = cs14_veg%frootc_patch(p)
            Xveg14c%V(p,ifroot_st)     = cs14_veg%frootc_storage_patch(p)
            Xveg14c%V(p,ifroot_xf)     = cs14_veg%frootc_xfer_patch(p)
            Xveg14c%V(p,ilivestem)     = cs14_veg%livestemc_patch(p)
            Xveg14c%V(p,ilivestem_st)  = cs14_veg%livestemc_storage_patch(p)
            Xveg14c%V(p,ilivestem_xf)  = cs14_veg%livestemc_xfer_patch(p)
            Xveg14c%V(p,ideadstem)     = cs14_veg%deadstemc_patch(p)
            Xveg14c%V(p,ideadstem_st)  = cs14_veg%deadstemc_storage_patch(p)
            Xveg14c%V(p,ideadstem_xf)  = cs14_veg%deadstemc_xfer_patch(p)
            Xveg14c%V(p,ilivecroot)    = cs14_veg%livecrootc_patch(p)
            Xveg14c%V(p,ilivecroot_st) = cs14_veg%livecrootc_storage_patch(p)
            Xveg14c%V(p,ilivecroot_xf) = cs14_veg%livecrootc_xfer_patch(p)
            Xveg14c%V(p,ideadcroot)    = cs14_veg%deadcrootc_patch(p)
            Xveg14c%V(p,ideadcroot_st) = cs14_veg%deadcrootc_storage_patch(p)
            Xveg14c%V(p,ideadcroot_xf) = cs14_veg%deadcrootc_xfer_patch(p)
         end do

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            if(ivt(p) >= npcropmin)then
               Xveg14c%V(p,igrain)     = cs14_veg%grainc_patch(p)
               Xveg14c%V(p,igrain_st)  = cs14_veg%grainc_storage_patch(p)
               Xveg14c%V(p,igrain_xf)  = cs14_veg%grainc_xfer_patch(p)
            end if
         end do
      end if
             
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         Xvegn%V(p,ileaf)         = leafn(p)
         Xvegn%V(p,ileaf_st)      = leafn_storage(p)
         Xvegn%V(p,ileaf_xf)      = leafn_xfer(p)
         Xvegn%V(p,ifroot)        = frootn(p)
         Xvegn%V(p,ifroot_st)     = frootn_storage(p)
         Xvegn%V(p,ifroot_xf)     = frootn_xfer(p)
         Xvegn%V(p,ilivestem)     = livestemn(p)
         Xvegn%V(p,ilivestem_st)  = livestemn_storage(p)
         Xvegn%V(p,ilivestem_xf)  = livestemn_xfer(p)
         Xvegn%V(p,ideadstem)     = deadstemn(p)
         Xvegn%V(p,ideadstem_st)  = deadstemn_storage(p)
         Xvegn%V(p,ideadstem_xf)  = deadstemn_xfer(p)
         Xvegn%V(p,ilivecroot)    = livecrootn(p)
         Xvegn%V(p,ilivecroot_st) = livecrootn_storage(p)
         Xvegn%V(p,ilivecroot_xf) = livecrootn_xfer(p)
         Xvegn%V(p,ideadcroot)    = deadcrootn(p)
         Xvegn%V(p,ideadcroot_st) = deadcrootn_storage(p)
         Xvegn%V(p,ideadcroot_xf) = deadcrootn_xfer(p)
         Xvegn%V(p,iretransn)     = retransn(p)
      end do

            do fp = 1,num_soilp
               p = filter_soilp(fp)
!               if(p .eq. 6851)print*,'after ini0',p,Xvegc%V(p,1:nvegcpool)
            end do

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         if(ivt(p) >= npcropmin)then
            Xvegn%V(p,igrain)        = grainn(p)
            Xvegn%V(p,igrain_st)     = grainn_storage(p)
            Xvegn%V(p,igrain_xf)     = grainn_xfer(p)
         end if
      end do

      if (is_beg_curr_year())then
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            leafc0(p)                = max(leafc(p),  epsi)
            leafc0_storage(p)        = max(leafc_storage(p),  epsi)
            leafc0_xfer(p)           = max(leafc_xfer(p),  epsi)
            frootc0(p)               = max(frootc(p),  epsi)
            frootc0_storage(p)       = max(frootc_storage(p),  epsi)
            frootc0_xfer(p)          = max(frootc_xfer(p),  epsi)
            livestemc0(p)            = max(livestemc(p),  epsi)
            livestemc0_storage(p)    = max(livestemc_storage(p),  epsi)
            livestemc0_xfer(p)       = max(livestemc_xfer(p),  epsi)
            deadstemc0(p)            = max(deadstemc(p),  epsi)
            deadstemc0_storage(p)    = max(deadstemc_storage(p),  epsi)
            deadstemc0_xfer(p)       = max(deadstemc_xfer(p),  epsi)
            livecrootc0(p)           = max(livecrootc(p),  epsi)
            livecrootc0_storage(p)   = max(livecrootc_storage(p),  epsi)
            livecrootc0_xfer(p)      = max(livecrootc_xfer(p),  epsi)
            deadcrootc0(p)           = max(deadcrootc(p),  epsi)
            deadcrootc0_storage(p)   = max(deadcrootc_storage(p),  epsi)
            deadcrootc0_xfer(p)      = max(deadcrootc_xfer(p),  epsi)
         end do

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            if(ivt(p) >= npcropmin)then
               grainc0(p)               = max(grainc(p),  epsi)
               grainc0_storage(p)       = max(grainc_storage(p),  epsi)
               grainc0_xfer(p)          = max(grainc_xfer(p),  epsi)
            end if
         end do

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            leafn0(p)                = max(leafn(p),  epsi)
            leafn0_storage(p)        = max(leafn_storage(p),  epsi)
            leafn0_xfer(p)           = max(leafn_xfer(p),  epsi)
            frootn0(p)               = max(frootn(p),  epsi)
            frootn0_storage(p)       = max(frootn_storage(p),  epsi)
            frootn0_xfer(p)          = max(frootn_xfer(p),  epsi)
            livestemn0(p)            = max(livestemn(p),  epsi)
            livestemn0_storage(p)    = max(livestemn_storage(p),  epsi)
            livestemn0_xfer(p)       = max(livestemn_xfer(p),  epsi)
            deadstemn0(p)            = max(deadstemn(p),  epsi)
            deadstemn0_storage(p)    = max(deadstemn_storage(p),  epsi)
            deadstemn0_xfer(p)       = max(deadstemn_xfer(p),  epsi)
            livecrootn0(p)           = max(livecrootn(p),  epsi)
            livecrootn0_storage(p)   = max(livecrootn_storage(p),  epsi)
            livecrootn0_xfer(p)      = max(livecrootn_xfer(p),  epsi)
            deadcrootn0(p)           = max(deadcrootn(p),  epsi)
            deadcrootn0_storage(p)   = max(deadcrootn_storage(p),  epsi)
            deadcrootn0_xfer(p)      = max(deadcrootn_xfer(p),  epsi)
            retransn0(p)             = max(retransn(p),  epsi)
         end do

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            if(ivt(p) >= npcropmin)then
               grainn0(p)               = max(grainn(p),  epsi)
               grainn0_storage(p)       = max(grainn_storage(p),  epsi)
               grainn0_xfer(p)          = max(grainn_xfer(p),  epsi)
            end if
         end do
      end if
            do fp = 1,num_soilp
               p = filter_soilp(fp)
!               if(p .eq. 6851)print*,'after ini',p,Xvegc%V(p,1:nvegcpool)
            end do
!         leafc0(p)                = max(leafc0(p),epsi)
!         leafc0_storage(p)        = max(leafc0_storage(p), epsi)
!         leafc0_xfer(p)           = max(leafc0_xfer(p), epsi)
!         frootc0(p)               = max(frootc0(p), epsi)
!         frootc0_storage(p)       = max(frootc0_storage(p), epsi)
!         frootc0_xfer(p)          = max(frootc0_xfer(p), epsi)
!         livestemc0(p)            = max(livestemc0(p), epsi)
!         livestemc0_storage(p)    = max(livestemc0_storage(p), epsi)
!         livestemc0_xfer(p)       = max(livestemc0_xfer(p), epsi)
!         deadstemc0(p)            = max(deadstemc0(p), epsi)
!         deadstemc0_storage(p)    = max(deadstemc0_storage(p), epsi)
!         deadstemc0_xfer(p)       = max(deadstemc0_xfer(p), epsi)
!         livecrootc0(p)           = max(livecrootc0(p), epsi)
!         livecrootc0_storage(p)   = max(livecrootc0_storage(p), epsi)
!         livecrootc0_xfer(p)      = max(livecrootc0_xfer(p), epsi)
!         deadcrootc0(p)           = max(deadcrootc0(p), epsi)
!         deadcrootc0_storage(p)   = max(deadcrootc0_storage(p), epsi)
!         deadcrootc0_xfer(p)      = max(deadcrootc0_xfer(p), epsi)
!         if(ivt(p) >= npcropmin)then
!            grainc0(p)               = max(grainc0(p), epsi)
!            grainc0_storage(p)       = max(grainc0_storage(p), epsi)
!            grainc0_xfer(p)          = max(grainc0_xfer(p), epsi)
!         end if
!         leafn0(p)                = max(leafn0(p),epsi)
!         leafn0_storage(p)        = max(leafn0_storage(p), epsi)
!         leafn0_xfer(p)           = max(leafn0_xfer(p), epsi)
!         frootn0(p)               = max(frootn0(p), epsi)
!         frootn0_storage(p)       = max(frootn0_storage(p), epsi)
!         frootn0_xfer(p)          = max(frootn0_xfer(p), epsi)
!         livestemn0(p)            = max(livestemn0(p), epsi)
!         livestemn0_storage(p)    = max(livestemn0_storage(p), epsi)
!         livestemn0_xfer(p)       = max(livestemn0_xfer(p), epsi)
!         deadstemn0(p)            = max(deadstemn0(p), epsi)
!         deadstemn0_storage(p)    = max(deadstemn0_storage(p), epsi)
!         deadstemn0_xfer(p)       = max(deadstemn0_xfer(p), epsi)
!         livecrootn0(p)           = max(livecrootn0(p), epsi)
!         livecrootn0_storage(p)   = max(livecrootn0_storage(p), epsi)
!         livecrootn0_xfer(p)      = max(livecrootn0_xfer(p), epsi)
!         deadcrootn0(p)           = max(deadcrootn0(p), epsi)
!         deadcrootn0_storage(p)   = max(deadcrootn0_storage(p), epsi)
!         deadcrootn0_xfer(p)      = max(deadcrootn0_xfer(p), epsi)
!         retransn0(p)             = max(retransn0(p), epsi)
!         if(ivt(p) >= npcropmin)then
!            grainn0(p)               = max(grainn0(p), epsi)
!            grainn0_storage(p)       = max(grainn0_storage(p), epsi)
!            grainn0_xfer(p)          = max(grainn0_xfer(p), epsi)
!         end if
         call t_stopf('CN veg matrix-set old value')
         !print*,'input'
         call t_startf('CN veg matrix-matrix multi.')
            do fp = 1,num_soilp
               p = filter_soilp(fp)
!               if(p .eq. 6851)print*,'before input',p,Xvegc%V(p,1:nvegcpool)
            end do
         do i=1,nvegcpool
            do fp = 1,num_soilp
               p = filter_soilp(fp)
!               if(p .eq. 24)print*,'Cinput',i,matrix_alloc(p,i),matrix_Cinput(p) * dt
               vegmatrixc_input%V(p,i) = matrix_alloc(p,i) * matrix_Cinput(p) * dt
!               if(p .eq. 10580)print*,'input',i,matrix_alloc(p,i),matrix_Cinput(p) * dt
            end do
         end do

!       do fp = 1,num_soilp
!          p = filter_soilp(fp)
!          if(bounds%begp .eq. 20037)print*,'matrix_Cinput(p)',p,matrix_Cinput(p) * dt
!       end do
!         vegmatrixc_transfer(:,:) = (matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc_old(:,:)) &
!                                   + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc_old(:,:)) &
!                                   + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc_old(:,:))) * dt
!         vegmatrixc_new(:) = matmul(vegmatrixc_old(:,:),rowonec(:))  +  vegmatrixc_input(:) + matmul(vegmatrixc_transfer(:,:),rowonec(:)) 
         if(ncphtrans .gt. ncphouttrans)then
!            do i=1,ncphtrans-ncphouttrans
!               do fp = 1,num_soilp
!                  p = filter_soilp(fp)
!                  Aphconed(p,i) = matrix_phtransfer(p,i)
!               end do
!            end do
!            if(.not. init_ready_aphc)then
            AI_phc = receiver_phc(1:ncphtrans-ncphouttrans)
            AJ_phc = doner_phc   (1:ncphtrans-ncphouttrans)
!            end if
             
!            if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'Aphconed',Aphconed(3998,:)
            call AKphvegc%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Aphconed,&
                 AI_phc,AJ_phc,ncphtrans-ncphouttrans,init_ready_aphc,list_aphc,RI_phc,CI_phc)
!            if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'afteer AKphvegc%SetValueA',AKphvegc%M(3998,:)
         else
            call AKphvegc%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if
!         Kconed = matrix_phturnover(1:nvegcpool,p) * dt
         !print*,'here0.1'
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'before KvegSetValueDM,matrix_phturnover',matrix_phturnover(3998,:)
         call Kvegc%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_phturnover(bounds%begp:bounds%endp,1:nvegcpool))
         !print*,'here0.3'
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'before SPMP_AK,Kvegc',Kvegc%DM(3998,:)
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'before SPMP_AK,AKphvegc',AKphvegc%M(3998,:)
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'before SPMP_AK,AKphvegc%RI',AKphvegc%RI
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'before SPMP_AK,AKphvegc%CI',AKphvegc%CI
!         if(bounds%begp .le. 10580 .and. bounds%endp .ge. 10580)print*,'leafphk and frootphk',Kvegc%DM(10580,1),Kvegc%DM(10580,4),matrix_phturnover(10580,1),matrix_phturnover(10580,4)
         call AKphvegc%SPMM_AK(num_soilp,filter_soilp,Kvegc)
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'after SPMM_AK,AKphvegc',AKphvegc%M(3998,:)
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'after SPMM_AK,AKphvegc%RI',AKphvegc%RI
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'after SPMM_AK,AKphvegc%CI',AKphvegc%CI
!         A_C( = matrix_phtransfer(1:ncphtrans-ncphouttrans,p) 
!         K_C = matrix_phturnover(1:nvegcpool,p)
!         AK_C = 0
!         call dgemm('N','N',nvegcpool,nvegcpool,nvegcpool,dt,A_C,nvegcpool,K_C,nvegcpool,0._r8,AK_C,nvegcpool)
         !print*,'here1'
            do fp = 1,num_soilp
               p = filter_soilp(fp)
            end do
         if(ncgmtrans .gt. ncgmouttrans)then
!            Agmconed = matrix_gmtransfer(1:ncgmtrans-ncgmouttrans,p) 
!            if(.not. init_ready_agmc)then
            AI_gmc = receiver_gmc(1:ncgmtrans-ncgmouttrans)
            AJ_gmc = doner_gmc   (1:ncgmtrans-ncgmouttrans)
!            end if
            call AKgmvegc%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Agmconed,&
                 AI_gmc,AJ_gmc,ncgmtrans-ncgmouttrans,init_ready_agmc,list_agmc,RI_gmc,CI_gmc)
         else
            call AKgmvegc%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if
!         Kconed   = matrix_gmturnover(1:nvegcpool,p) * dt
         call Kvegc%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_gmturnover(bounds%begp:bounds%endp,1:nvegcpool))
         call AKgmvegc%SPMM_AK(num_soilp,filter_soilp,Kvegc)
!         A_C = matrix_gmtransfer(p,1:nvegcpool,:)
!         K_C = matrix_gmturnover(p,:,:)
!         call dgemm('N','N',nvegcpool,nvegcpool,nvegcpool,dt,A_C,nvegcpool,K_C,nvegcpool,1._r8,AK_C,nvegcpool)
         if(ncfitrans .gt. ncfiouttrans)then
!            Aficoned = matrix_fitransfer(1:ncfitrans-ncfiouttrans,p) 
!            if(.not. init_ready_afic)then
            AI_fic = receiver_fic(1:ncfitrans-ncfiouttrans)
            AJ_fic = doner_fic   (1:ncfitrans-ncfiouttrans)
!            end if
            call AKfivegc%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Aficoned,&
                 AI_fic,AJ_fic,ncfitrans-ncfiouttrans,init_ready_afic,list_afic,RI_fic,CI_fic)
         else
            call AKfivegc%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if
!         Kconed   = matrix_fiturnover(1:nvegcpool,p) * dt
         call Kvegc%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_fiturnover(bounds%begp:bounds%endp,1:nvegcpool))
         call AKfivegc%SPMM_AK(num_soilp,filter_soilp,Kvegc)
!         if(num_actfirep .ne. 0)then

!            do fp=1,num_soilp
!               p = filter_soilp(fp)
!               do i=1,AKfivegc%NE
!               end do
!            end do

!            call Xfirelossvegc%SPMM_AX(num_soilp,filter_soilp,AKfivegc)

 !           do fp=1,num_soilp
 !              p = filter_soilp(fp)
 !              do i = 1,nvegcpool
 !                 fire_closs (p) = fire_closs(p) + (Xfirelossvegc%V(p,i) - Xvegc%V(p,i)) / dt
 !              end do
 !           end do
 !        end if   
!         A_C = matrix_fitransfer(p,1:nvegcpool,:)
!         K_C = matrix_fiturnover(p,:,:)
!         call dgemm('N','N',nvegcpool,nvegcpool,nvegcpool,dt,A_C,nvegcpool,K_C,nvegcpool,1._r8,AK_C,nvegcpool)
!         call AKallvegc%SPMP_AB(num_soilp,filter_soilp,Agmvegc,list_ready_tmp1vegc,list_phc_phgmc(1:AKallvegc%NE),list_gmc_phgmc(1:Agmvegc%NE),&
!                                NE_AKphgmc,RI_AKphgmc,CI_AKphgmc)
!         call AKallvegc%SPMP_AB(num_soilp,filter_soilp,Afivegc,list_ready_tmp2vegc,list_phgmc_phgmfic(1:AKallvegc%NE),list_fic_phgmfic(1:Afivegc%NE),&
!                                NE_AKallvegc,RI_AKallvegc,CI_AKallvegc)
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'before SPMP_AB,AKfivegc',num_actfirep,AKfivegc%M(3998,:)
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'before SPMP_AB,AKfivegc%RI',AKfivegc%RI
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'before SPMP_AB,AKfivegc%CI',AKfivegc%CI
!            do fp = 1,num_soilp
!               p = filter_soilp(fp)
!               if(p .eq. 6851)print*,'before SPMP_AB',p,Xvegc%V(p,1:nvegcpool)
!            end do
!         if(bounds%begp .le. 10580 .and. bounds%endp .ge. 10580 .and. AKphvegc%NE .ge. 0)then
!            do i=1,AKphvegc%NE
!               if(AKphvegc%RI(i) .eq. 4)then
!                  print*,'AKphvegc,froot',AKphvegc%RI(i),AKphvegc%CI(i),AKphvegc%M(10580,i)
!               end if
!               if(AKphvegc%RI(i) .eq. 1)then
!                  print*,'AKphvegc,leaf',AKphvegc%RI(i),AKphvegc%CI(i),AKphvegc%M(10580,i)
!               end if
!            end do
!         end if
!         if(bounds%begp .le. 10580 .and. bounds%endp .ge. 10580 .and. AKgmvegc%NE .ge. 0)then
!            do i=1,AKgmvegc%NE
!               if(AKgmvegc%RI(i) .eq. 4)then
!                  print*,'AKgmvegc,froot',AKgmvegc%RI(i),AKgmvegc%CI(i),AKgmvegc%M(10580,i)
!               end if
!               if(AKgmvegc%RI(i) .eq. 1)then
!                  print*,'AKgmvegc,leaf',AKgmvegc%RI(i),AKgmvegc%CI(i),AKgmvegc%M(10580,i)
!               end if
!            end do
!         end if
!         if(bounds%begp .le. 10580 .and. bounds%endp .ge. 10580 .and. AKfivegc%NE .ge. 0)then
!            do i=1,AKfivegc%NE
!               if(AKfivegc%RI(i) .eq. 4)then
!                  print*,'AKfivegc,froot',AKfivegc%RI(i),AKfivegc%CI(i),AKfivegc%M(10580,i)
!               end if
!               if(AKfivegc%RI(i) .eq. 1)then
!                  print*,'AKfivegc,leaf',AKfivegc%RI(i),AKfivegc%CI(i),AKfivegc%M(10580,i)
!               end if
!            end do
!         end if
         if(num_actfirep .eq. 0)then
            call AKallvegc%SPMP_AB(num_soilp,filter_soilp,AKphvegc,AKgmvegc,list_ready_phgmc,list_A=list_phc_phgm,list_B=list_gmc_phgm,&
                 NE_AB=NE_AKallvegc,RI_AB=RI_AKallvegc,CI_AB=CI_AKallvegc)
         else
            call AKallvegc%SPMP_ABC(num_soilp,filter_soilp,AKphvegc,AKgmvegc,AKfivegc,list_ready_phgmfic,list_A=list_phc_phgmfi,&
                 list_B=list_gmc_phgmfi,list_C=list_fic_phgmfi,NE_ABC=NE_AKallvegc,RI_ABC=RI_AKallvegc,CI_ABC=CI_AKallvegc,&
                 use_actunit_list_C=.True.,num_actunit_C=num_actfirep,filter_actunit_C=filter_actfirep)
         end if
!         if(bounds%begp .le. 10580 .and. bounds%endp .ge. 10580)then
!            do i=1,NE_AKallvegc
!               if(AKallvegc%RI(i) .eq. 4)then
!                  print*,'AKallvegc,froot',AKallvegc%RI(i),AKallvegc%CI(i),AKallvegc%M(10580,i)
!!               end if
 !              if(AKallvegc%RI(i) .eq. 1)then
 !                 print*,'AKallvegc,leaf',AKallvegc%RI(i),AKallvegc%CI(i),AKallvegc%M(10580,i)
 !              end if
 !           end do
 !        end if
!         print*,'here5'
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'Xveg before SPMMAX',Xvegc%V(3998,:),nvegcpool,ivt(3998),ivt(3998),npcropmin
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'Akallvegc',AKallvegc%M(3998,:)
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'RI',AKallvegc%RI
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'CI',AKallvegc%CI
!         do i=1,nvegcpool
!            do fp = 1,num_soilp
!               p = filter_soilp(fp)
!               if(p .eq. 6851)print*,'before update Xvegcveg',p,Xvegc%V(p,1:nvegcpool),AKallvegc%M(p,1:10)
!               if(p .eq. 6851)print*,'RI',AKallvegc%RI(1:10)
!               if(p .eq. 6851)print*,'CI',AKallvegc%CI(1:10)
!            end do
!         end do
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'Xveg after SPMMAX',Xvegc%V(3998,:)
!         do i = 1,nvegcpool
!            do fp = 1,num_soilp
!               p = filter_soilp(fp)
!               if(p .eq. 10580)print*,'Xvegc,before AX to froot',p,i,Xvegc%V(p,i)
!            end do
!         end do
         call Xvegc%SPMM_AX(num_soilp,filter_soilp,AKallvegc)
!         do i = 1,nvegcpool
!            do fp = 1,num_soilp
!               p = filter_soilp(fp)
!               if(p .eq. 10580)print*,'Xvegc,after AX to froot',p,i,Xvegc%V(p,i)
!            end do
!         end do
         do i = 1,nvegcpool
            do fp = 1,num_soilp
               p = filter_soilp(fp)
!               if(p .eq. 6851 .and. i .eq. 5)print*,'Cinput to veg',Xvegc%V(p,i),vegmatrixc_input%V(p,i)
               Xvegc%V(p,i) = Xvegc%V(p,i) + vegmatrixc_input%V(p,i) 
!               if(p .eq. 10580)print*,'input to froot',p,i,Xvegc%V(p,i),vegmatrixc_input%V(p,i)
            end do
         end do
         

!         do fp = 1,num_soilp
!            p = filter_soilp(fp)
!            if(bounds%begp .eq. 20037)print*,'after update Xvegcveg',p,Xvegc%V(p,1:nvegcpool),cnveg_carbonstate_inst%cpool_patch(p),cnveg_carbonstate_inst%gresp_storage_patch(p),cnveg_carbonstate_inst%gresp_xfer_patch(p)
!            if(bounds%begp .eq. 20037)print*,'after update Xvegctotal',p,sum(Xvegc%V(p,1:nvegcpool))+cnveg_carbonstate_inst%cpool_patch(p)+cnveg_carbonstate_inst%gresp_storage_patch(p)+cnveg_carbonstate_inst%gresp_xfer_patch(p)
!         end do
!         if(bounds%begp .le. 3998 .and. bounds%endp .ge. 3998)print*,'Xveg after input',Xvegc%V(3998,:)
!         vegmatrixc_new = vegmatrixc_input + vegmatrixc_old
!         call dgemv('N',nvegcpool,nvegcpool,1._r8,AK_C,nvegcpool,vegmatrixc_old,1,1._r8,vegmatrixc_new,1)

         if ( use_c13 ) then

            do i=1,nvegcpool
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  vegmatrixc13_input%V(p,i) = matrix_alloc(p,i) * matrix_C13input(p) * dt
               end do
            end do
!            vegmatrixc13_transfer(:,:) = (matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc13_old(:,:)) &
!                                   + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc13_old(:,:)) &
!                                   + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc13_old(:,:))) * dt
!            vegmatrixc13_new(:) = matmul(vegmatrixc13_old(:,:),rowonec(:))  +  vegmatrixc13_input(:) + matmul(vegmatrixc13_transfer(:,:),rowonec(:))

!            A = matrix_phtransfer(p,1:nvegcpool,:)
!            B = matrix_phturnover(p,:,:)
!            C = 0
!            call dgemm('N','N',nvegcpool,nvegcpool,nvegcpool,1._r8,A,nvegcpool,B,nvegcpool,0._r8,C,nvegcpool)
!            A = matrix_gmtransfer(p,1:nvegcpool,:)
!            B = matrix_gmturnover(p,:,:)
!            call dgemm('N','N',nvegcpool,nvegcpool,nvegcpool,1._r8,A,nvegcpool,B,nvegcpool,1._r8,C,nvegcpool)
!            A = matrix_fitransfer(p,1:nvegcpool,:)
!            B = matrix_fiturnover(p,:,:)
!            call dgemm('N','N',nvegcpool,nvegcpool,nvegcpool,1._r8,A,nvegcpool,B,nvegcpool,1._r8,C,nvegcpool)

            call Xveg13c%SPMM_AX(num_soilp,filter_soilp,AKallvegc)
            do i=1,nvegcpool
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  Xveg13c%V(p,i) = Xveg13c%V(p,i) + vegmatrixc13_input%V(p,i) 
               end do
            end do

!            vegmatrixc13_new = vegmatrixc13_input + vegmatrixc13_old
!            call dgemv('N',nvegcpool,nvegcpool,1._r8,AK_C,nvegcpool,vegmatrixc13_old,1,1._r8,vegmatrixc13_new,1)
         end if
         if ( use_c14 ) then
            do i=1,nvegcpool
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  vegmatrixc14_input%V(p,i) = matrix_alloc(p,i) * matrix_C14input(p) * dt
               end do
            end do
!            vegmatrixc14_transfer(:,:) = (matmul(matmul(matrix_phtransfer(p,1:nvegcpool,:),matrix_phturnover(p,:,:)),vegmatrixc14_old(:,:)) &
!                                   + matmul(matmul(matrix_gmtransfer(p,1:nvegcpool,:),matrix_gmturnover(p,:,:)),vegmatrixc14_old(:,:)) &
!                                   + matmul(matmul(matrix_fitransfer(p,1:nvegcpool,:),matrix_fiturnover(p,:,:)),vegmatrixc14_old(:,:))) * dt
!            vegmatrixc14_new(:) = matmul(vegmatrixc14_old(:,:),rowonec(:))  +  vegmatrixc14_input(:) + matmul(vegmatrixc14_transfer(:,:),rowonec(:))
         !print*,'here6'


            call Xveg14c%SPMM_AX(num_soilp,filter_soilp,AKallvegc)
            do i=1,nvegcpool
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  Xveg14c%V(p,i) = Xveg14c%V(p,i) + vegmatrixc14_input%V(p,i) 
               end do
            end do
!            vegmatrixc14_new = vegmatrixc14_input + vegmatrixc14_old
!            call dgemv('N',nvegcpool,nvegcpool,1._r8,AK_C,nvegcpool,vegmatrixc14_old,1,1._r8,vegmatrixc14_new,1)
         end if
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
                            
!         if(p .eq. -9999)then
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
!         end if
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
         !print*,'here7'
         do i=1,nvegnpool
            do fp = 1,num_soilp
               p = filter_soilp(fp)
!               if(p .eq. 2111 .and. i .eq. 1)print*,'matrix_nalloc',matrix_nalloc(p,i),matrix_Ninput(p) * dt
               vegmatrixn_input%V(p,i) = matrix_nalloc(p,i) * matrix_Ninput(p) * dt
            end do
         end do
         !print*,'here7.1'
!         vegmatrixn_transfer(:,:) = (matmul(matmul(matrix_nphtransfer(p,1:nvegnpool,:),matrix_nphturnover(p,:,:)),vegmatrixn_old(:,:)) &
!                                   + matmul(matmul(matrix_ngmtransfer(p,1:nvegnpool,:),matrix_ngmturnover(p,:,:)),vegmatrixn_old(:,:)) &
!                                   + matmul(matmul(matrix_nfitransfer(p,1:nvegnpool,:),matrix_nfiturnover(p,:,:)),vegmatrixn_old(:,:))) * dt
!         vegmatrixn_new(:) = matmul(vegmatrixn_old(:,:),rowonen(:)) +  vegmatrixn_input(:) + matmul(vegmatrixn_transfer(:,:),rowonen(:))

         if(nnphtrans .gt. nnphouttrans)then
!            Aphnoned = matrix_nphtransfer(1:nnphtrans-nnphouttrans,p)
!            if(.not. init_ready_aphn)then
             AI_phn = receiver_phn(1:nnphtrans-nnphouttrans)
             AJ_phn = doner_phn   (1:nnphtrans-nnphouttrans)
!            end if
            call AKphvegn%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Aphnoned,&
                 AI_phn,AJ_phn,nnphtrans-nnphouttrans,init_ready_aphn,list_aphn,RI_phn,CI_phn)
         else
            call AKphvegn%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if
!         Knoned   = matrix_nphturnover(1:nvegnpool,p) * dt
         !print*,'here7.2'
         call Kvegn%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_nphturnover(bounds%begp:bounds%endp,1:nvegnpool))
         call AKphvegn%SPMM_AK(num_soilp,filter_soilp,Kvegn)
!         if(bounds%begp .le. 2111 .and. bounds%endp .ge. 2111)print*,'AKphvegn',AKphvegn%M(2111,1),AKphvegn%RI(1),AKphvegn%CI(1),AKphvegn%M(2111,2),AKphvegn%RI(2),AKphvegn%CI(2)
         !print*,'here7.3'
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'AKphvegn',AKphvegn%M(6851,1:5)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'RI',AKphvegn%RI(1:5)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'CI',AKphvegn%CI(1:5)

         if(nngmtrans .gt. nngmouttrans)then
!            Agmnoned = matrix_ngmtransfer(1:nngmtrans-nngmouttrans,p)
!            if(init_ready_agmn)then
             AI_gmn = receiver_gmn(1:nngmtrans-nngmouttrans)
             AJ_gmn = doner_gmn   (1:nngmtrans-nngmouttrans)
!            end if
            call AKgmvegn%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Agmnoned,&
                 AI_gmn,AJ_gmn,nngmtrans-nngmouttrans,init_ready_agmn,list_agmn,RI_gmn,CI_gmn)
         else
            call AKgmvegn%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if
!         Knoned   = matrix_ngmturnover(1:nvegnpool,p) * dt
         !print*,'here7.4'
         call Kvegn%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_ngmturnover(bounds%begp:bounds%endp,1:nvegnpool))
         call AKgmvegn%SPMM_AK(num_soilp,filter_soilp,Kvegn)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'AKgmvegn',AKgmvegn%M(6851,1:5)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'RI',AKgmvegn%RI(1:5)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'CI',AKgmvegn%CI(1:5)
!         if(bounds%begp .le. 2111 .and. bounds%endp .ge. 2111)print*,'AKgmvegn',AKgmvegn%M(2111,1),AKphvegn%RI(1),AKphvegn%CI(1)
         !print*,'here7.5'
!         A_C = matrix_gmtransfer(p,1:nvegcpool,:)
!         K_C = matrix_gmturnover(p,:,:)
!         call dgemm('N','N',nvegcpool,nvegcpool,nvegcpool,dt,A_C,nvegcpool,K_C,nvegcpool,1._r8,AK_C,nvegcpool)
         if(nnfitrans .gt. nnfiouttrans)then
!            Afinoned = matrix_nfitransfer(1:nnfitrans-nnfiouttrans,p)
!            if(init_ready_afin)then
             AI_fin = receiver_fin(1:nnfitrans-nnfiouttrans)
             AJ_fin = doner_fin   (1:nnfitrans-nnfiouttrans)
!            end if
            call AKfivegn%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Afinoned,&
                 AI_fin,AJ_fin,nnfitrans-nnfiouttrans,init_ready_afin,list_afin,RI_fin,CI_fin)
         else
            call AKfivegn%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if
         !print*,'here7.6'
!         Knoned   = matrix_nfiturnover(1:nvegnpool,p) * dt
         call Kvegn%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_nfiturnover(bounds%begp:bounds%endp,1:nvegnpool))
         call AKfivegn%SPMM_AK(num_soilp,filter_soilp,Kvegn)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'AKfivegn',AKfivegn%M(6851,1:5)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'RI',AKfivegn%RI(1:5)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'CI',AKfivegn%CI(1:5)
!         if(bounds%begp .le. 2111 .and. bounds%endp .ge. 2111)print*,'AKfivegn',AKfivegn%M(2111,1),AKfivegn%RI(1),AKfivegn%CI(1)
         !print*,'here7.7'
!         A_C = matrix_fitransfer(p,1:nvegcpool,:)
!         K_C = matrix_fiturnover(p,:,:)
!         call dgemm('N','N',nvegcpool,nvegcpool,nvegcpool,dt,A_C,nvegcpool,K_C,nvegcpool,1._r8,AK_C,nvegcpool)
!         call AKallvegn%SPMP_AB(num_soilp,filter_soilp,AKgmvegn,list_ready_tmp1vegn,list_phn_phgmn(1:AKallvegn%NE),list_gmn_phgmn(1:AKgmvegn%NE),&
!                                NE_AKphgmn,RI_AKphgmn,CI_AKphgmn)
!         call AKallvegn%SPMP_AB(num_soilp,filter_soilp,Afivegn,list_ready_tmp2vegn,list_phgmn_phgmfin(1:AKallvegn%NE),list_fin_phgmfin(1:Afivegn%NE),&
!                                NE_AKallvegn,RI_AKallvegn,CI_AKallvegn)
         
         if(num_actfirep .eq. 0)then
            call AKallvegn%SPMP_AB(num_soilp,filter_soilp,AKphvegn,AKgmvegn,list_ready_phgmn,list_A=list_phn_phgm,list_B=list_gmn_phgm,&
                 NE_AB=NE_AKallvegn,RI_AB=RI_AKallvegn,CI_AB=CI_AKallvegn)
         else
            call AKallvegn%SPMP_ABC(num_soilp,filter_soilp,AKphvegn,AKgmvegn,AKfivegn,list_ready_phgmfin,list_A=list_phn_phgmfi,&
                 list_B=list_gmn_phgmfi,list_C=list_fin_phgmfi,NE_ABC=NE_AKallvegn,RI_ABC=RI_AKallvegn,CI_ABC=CI_AKallvegn,&
                 use_actunit_list_C=.True.,num_actunit_C=num_actfirep,filter_actunit_C=filter_actfirep)
         end if
         
!         list_ready = .TRUE.

         !print*,'here7.8'ph
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'Xvegn,before update',Xvegn%V(6851,1:19)
!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'Xvegn,after update',Xvegn%V(6851,1:19)
         call Xvegn%SPMM_AX(num_soilp,filter_soilp,AKallvegn)
          
         do i=1,nvegnpool
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               Xvegn%V(p,i) = Xvegn%V(p,i) + vegmatrixn_input%V(p,i)
            end do
         end do

!         if(bounds%begp .le. 6851 .and. bounds%endp .ge. 6851)print*,'Xvegn,after input',Xvegn%V(6851,1:19),vegmatrixn_input%V(6851,1:19)
         !print*,'here7.9'
!         A_N = matrix_nphtransfer(p,1:nvegnpool,:)
!         K_N = matrix_nphturnover(p,:,:)
!               Xvegn%V(p,i) = Xvegn%V(p,i) + vegmatrixn_input%V(p,i)
!            end do
!         end do
!         A_N = matrix_nphtransfer(p,1:nvegnpool,:)
!         K_N = matrix_nphturnover(p,:,:)
!         AK_N = 0
!         call dgemm('N','N',nvegnpool,nvegnpool,nvegnpool,dt,A_N,nvegnpool,K_N,nvegnpool,0._r8,AK_N,nvegnpool)
!         A_N = matrix_ngmtransfer(p,1:nvegnpool,:)
!         K_N = matrix_ngmturnover(p,:,:)
!         call dgemm('N','N',nvegnpool,nvegnpool,nvegnpool,dt,A_N,nvegnpool,K_N,nvegnpool,1._r8,AK_N,nvegnpool)
!         A_N = matrix_nfitransfer(p,1:nvegcpool,:)
!         K_N = matrix_nfiturnover(p,:,:)
!         call dgemm('N','N',nvegnpool,nvegnpool,nvegnpool,dt,A_N,nvegnpool,K_N,nvegnpool,1._r8,AK_N,nvegnpool)
!
!         vegmatrixn_new = vegmatrixn_input + vegmatrixn_old
!         call dgemv('N',nvegnpool,nvegnpool,1._r8,AK_N,nvegnpool,vegmatrixn_old,1,1._r8,vegmatrixn_new,1)

         call t_stopf('CN veg matrix-matrix multi.')

!         if(p .eq. 16)then
!         end if
 
! 

         !print*,'here8'
         call t_startf('CN veg matrix-accum. trans.')
         if(isspinup .or. is_outmatrix)then
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               matrix_calloc_leaf_acc(p)        = matrix_calloc_leaf_acc(p)        + vegmatrixc_input%V(p,ileaf)
               matrix_calloc_leafst_acc(p)      = matrix_calloc_leafst_acc(p)      + vegmatrixc_input%V(p,ileaf_st)
               matrix_calloc_froot_acc(p)       = matrix_calloc_froot_acc(p)       + vegmatrixc_input%V(p,ifroot)
               matrix_calloc_frootst_acc(p)     = matrix_calloc_frootst_acc(p)     + vegmatrixc_input%V(p,ifroot_st)
               matrix_calloc_livestem_acc(p)    = matrix_calloc_livestem_acc(p)    + vegmatrixc_input%V(p,ilivestem)
               matrix_calloc_livestemst_acc(p)  = matrix_calloc_livestemst_acc(p)  + vegmatrixc_input%V(p,ilivestem_st)
               matrix_calloc_deadstem_acc(p)    = matrix_calloc_deadstem_acc(p)    + vegmatrixc_input%V(p,ideadstem)
               matrix_calloc_deadstemst_acc(p)  = matrix_calloc_deadstemst_acc(p)  + vegmatrixc_input%V(p,ideadstem_st)
               matrix_calloc_livecroot_acc(p)   = matrix_calloc_livecroot_acc(p)   + vegmatrixc_input%V(p,ilivecroot)
               matrix_calloc_livecrootst_acc(p) = matrix_calloc_livecrootst_acc(p) + vegmatrixc_input%V(p,ilivecroot_st)
               matrix_calloc_deadcroot_acc(p)   = matrix_calloc_deadcroot_acc(p)   + vegmatrixc_input%V(p,ideadcroot)
               matrix_calloc_deadcrootst_acc(p) = matrix_calloc_deadcrootst_acc(p) + vegmatrixc_input%V(p,ideadcroot_st)
            end do
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  matrix_calloc_grain_acc(p)    = matrix_calloc_grain_acc(p)       + vegmatrixc_input%V(p,igrain)
                  matrix_calloc_grainst_acc(p)  = matrix_calloc_grainst_acc(p)     + vegmatrixc_input%V(p,igrain_st)
               end if
            end do


         !print*,'here8,5'
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               matrix_ctransfer_leafst_to_leafxf_acc(p)           = matrix_ctransfer_leafst_to_leafxf_acc(p) &
                                                                  + matrix_phtransfer(p,ileafst_to_ileafxf_phc) &
                                                                  * dt * leafc_storage(p) !matrix_phturnover(p,ileaf_st)*leafc_storage(p)
               matrix_ctransfer_leafxf_to_leaf_acc(p)             = matrix_ctransfer_leafxf_to_leaf_acc(p) &
                                                                  + matrix_phtransfer(p,ileafxf_to_ileaf_phc) &
                                                                  * dt * leafc_xfer(p)!matrix_phturnover(p,ileaf_xf)*leafc_xfer(p)
               matrix_ctransfer_frootst_to_frootxf_acc(p)         = matrix_ctransfer_frootst_to_frootxf_acc(p) &
                                                                  + matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) &
                                                                  * dt * frootc_storage(p)!matrix_phturnover(p,ifroot_st)*frootc_storage(p)
               matrix_ctransfer_frootxf_to_froot_acc(p)           = matrix_ctransfer_frootxf_to_froot_acc(p) &
                                                                  + matrix_phtransfer(p,ifrootxf_to_ifroot_phc) &
                                                                  * dt * frootc_xfer(p)!matrix_phturnover(p,ifroot_xf)*frootc_xfer(p)
               matrix_ctransfer_livestemst_to_livestemxf_acc(p)   = matrix_ctransfer_livestemst_to_livestemxf_acc(p) &
                                                                  + matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) &
                                                                  * dt * livestemc_storage(p)!matrix_phturnover(p,ilivestem_st)*livestemc_storage(p)
               matrix_ctransfer_livestemxf_to_livestem_acc(p)     = matrix_ctransfer_livestemxf_to_livestem_acc(p) &
                                                                  + matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc) &
                                                                  * dt * livestemc_xfer(p)!matrix_phturnover(p,ilivestem_xf)*livestemc_xfer(p)
               matrix_ctransfer_deadstemst_to_deadstemxf_acc(p)   = matrix_ctransfer_deadstemst_to_deadstemxf_acc(p) &
                                                                  + matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) &
                                                                  * dt * deadstemc_storage(p)!matrix_phturnover(p,ideadstem_st)*deadstemc_storage(p)
               matrix_ctransfer_deadstemxf_to_deadstem_acc(p)     = matrix_ctransfer_deadstemxf_to_deadstem_acc(p) &
                                                                  + matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc) &
                                                                  * dt * deadstemc_xfer(p)!matrix_phturnover(p,ideadstem_xf)*deadstemc_xfer(p)
               matrix_ctransfer_livecrootst_to_livecrootxf_acc(p) = matrix_ctransfer_livecrootst_to_livecrootxf_acc(p) &
                                                                  + matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) &
                                                                  * dt * livecrootc_storage(p)!matrix_phturnover(p,ilivecroot_st)*livecrootc_storage(p)
               matrix_ctransfer_livecrootxf_to_livecroot_acc(p)   = matrix_ctransfer_livecrootxf_to_livecroot_acc(p) &
                                                                  + matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) &
                                                                  * dt * livecrootc_xfer(p)!matrix_phturnover(p,ilivecroot_xf)*livecrootc_xfer(p)
               matrix_ctransfer_deadcrootst_to_deadcrootxf_acc(p) = matrix_ctransfer_deadcrootst_to_deadcrootxf_acc(p) &
                                                                  + matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) &
                                                                  * dt * deadcrootc_storage(p)!matrix_phturnover(p,ideadcroot_st)*deadcrootc_storage(p)
               matrix_ctransfer_deadcrootxf_to_deadcroot_acc(p)   = matrix_ctransfer_deadcrootxf_to_deadcroot_acc(p) &
                                                                  + matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) &
                                                                  * dt * deadcrootc_xfer(p)!matrix_phturnover(p,ideadcroot_st)*deadcrootc_xfer(p)
               matrix_ctransfer_livestem_to_deadstem_acc(p)       = matrix_ctransfer_livestem_to_deadstem_acc(p) &
                                                                  +(matrix_phtransfer(p,ilivestem_to_ideadstem_phc)&!matrix_phturnover(p,ilivestem) &
                                                                  + matrix_fitransfer(p,ilivestem_to_ideadstem_fic))&!matrix_fiturnover(p,ilivestem))&
                                                                  * dt * livestemc(p) 
               matrix_ctransfer_livecroot_to_deadcroot_acc(p)     = matrix_ctransfer_livecroot_to_deadcroot_acc(p) &
                                                                  +(matrix_phtransfer(p,ilivecroot_to_ideadcroot_phc)&!*matrix_phturnover(p,ilivecroot) &
                                                                  + matrix_fitransfer(p,ilivecroot_to_ideadcroot_fic))&!*matrix_fiturnover(p,ilivecroot))&
                                                                  * dt * livecrootc(p)
!            if(ivt(p) >= npcropmin)then
!               matrix_ctransfer_grainst_to_grainxf_acc(p)      = matrix_ctransfer_grainst_to_grainxf_acc(p) &
!                                                               + matrix_phAK_C(igrain_xf,igrain_st)*Xoldvegc%V(igrain_st)
!                                                               + vegmatrixc_transfer(igrain_xf,igrain_st)
!               matrix_ctransfer_grainxf_to_grain_acc(p)        = matrix_ctransfer_grainxf_to_grain_acc(p) &
!                                                               + AK_C(igrain,igrain_xf)*Xoldvegc%V(igrain_xf)
!                                                               + vegmatrixc_transfer(igrain,igrain_xf)
!            end if
               matrix_cturnover_leaf_acc(p)        = matrix_cturnover_leaf_acc(p) &
                                                   + (matrix_phturnover(p,ileaf)+matrix_gmturnover(p,ileaf)+matrix_fiturnover(p,ileaf)) &
                                                   * leafc(p)
               matrix_cturnover_leafst_acc(p)      = matrix_cturnover_leafst_acc(p) &
                                                   + (matrix_phturnover(p,ileaf_st)+matrix_gmturnover(p,ileaf_st)+matrix_fiturnover(p,ileaf_st)) &
                                                   * leafc_storage(p)
               matrix_cturnover_leafxf_acc(p)      = matrix_cturnover_leafxf_acc(p) &
                                                   + (matrix_phturnover(p,ileaf_xf)+matrix_gmturnover(p,ileaf_xf)+matrix_fiturnover(p,ileaf_xf)) &
                                                   * leafc_xfer(p)
               matrix_cturnover_froot_acc(p)       = matrix_cturnover_froot_acc(p) &
                                                   + (matrix_phturnover(p,ifroot)+matrix_gmturnover(p,ifroot)+matrix_fiturnover(p,ifroot)) &
                                                   * frootc(p)
               matrix_cturnover_frootst_acc(p)     = matrix_cturnover_frootst_acc(p) &
                                                   + (matrix_phturnover(p,ifroot_st)+matrix_gmturnover(p,ifroot_st)+matrix_fiturnover(p,ifroot_st)) &
                                                   * frootc_storage(p)
               matrix_cturnover_frootxf_acc(p)     = matrix_cturnover_frootxf_acc(p) &
                                                   + (matrix_phturnover(p,ifroot_xf)+matrix_gmturnover(p,ifroot_xf)+matrix_fiturnover(p,ifroot_xf)) &
                                                   * frootc_xfer(p)
               matrix_cturnover_livestem_acc(p)    = matrix_cturnover_livestem_acc(p) &
                                                   + (matrix_phturnover(p,ilivestem)+matrix_gmturnover(p,ilivestem)+matrix_fiturnover(p,ilivestem)) &
                                                   * livestemc(p)
               matrix_cturnover_livestemst_acc(p)  = matrix_cturnover_livestemst_acc(p) &
                                                   + (matrix_phturnover(p,ilivestem_st)+matrix_gmturnover(p,ilivestem_st)+matrix_fiturnover(p,ilivestem_st)) &
                                                   * livestemc_storage(p)
               matrix_cturnover_livestemxf_acc(p)  = matrix_cturnover_livestemxf_acc(p) &
                                                   + (matrix_phturnover(p,ilivestem_xf)+matrix_gmturnover(p,ilivestem_xf)+matrix_fiturnover(p,ilivestem_xf)) &
                                                   * livestemc_xfer(p)
               matrix_cturnover_deadstem_acc(p)    = matrix_cturnover_deadstem_acc(p) &
                                                   + (matrix_phturnover(p,ideadstem)+matrix_gmturnover(p,ideadstem)+matrix_fiturnover(p,ideadstem)) &
                                                   * deadstemc(p)
               matrix_cturnover_deadstemst_acc(p)  = matrix_cturnover_deadstemst_acc(p) &
                                                   + (matrix_phturnover(p,ideadstem_st)+matrix_gmturnover(p,ideadstem_st)+matrix_fiturnover(p,ideadstem_st)) &
                                                   * deadstemc_storage(p)
               matrix_cturnover_deadstemxf_acc(p)  = matrix_cturnover_deadstemxf_acc(p) &
                                                   + (matrix_phturnover(p,ideadstem_xf)+matrix_gmturnover(p,ideadstem_xf)+matrix_fiturnover(p,ideadstem_xf)) &
                                                   * deadstemc_xfer(p)
               matrix_cturnover_livecroot_acc(p)   = matrix_cturnover_livecroot_acc(p) &
                                                   + (matrix_phturnover(p,ilivecroot)+matrix_gmturnover(p,ilivecroot)+matrix_fiturnover(p,ilivecroot)) &
                                                   * livecrootc(p)
               matrix_cturnover_livecrootst_acc(p) = matrix_cturnover_livecrootst_acc(p) &
                                                   + (matrix_phturnover(p,ilivecroot_st)+matrix_gmturnover(p,ilivecroot_st)+matrix_fiturnover(p,ilivecroot_st)) &
                                                   * livecrootc_storage(p)
               matrix_cturnover_livecrootxf_acc(p) = matrix_cturnover_livecrootxf_acc(p) &
                                                   + (matrix_phturnover(p,ilivecroot_xf)+matrix_gmturnover(p,ilivecroot_xf)+matrix_fiturnover(p,ilivecroot_xf)) &
                                                   * livecrootc_xfer(p)
               matrix_cturnover_deadcroot_acc(p)   = matrix_cturnover_deadcroot_acc(p) &
                                                   + (matrix_phturnover(p,ideadcroot)+matrix_gmturnover(p,ideadcroot)+matrix_fiturnover(p,ideadcroot)) & 
                                                   * deadcrootc(p)
               matrix_cturnover_deadcrootst_acc(p) = matrix_cturnover_deadcrootst_acc(p) &
                                                   + (matrix_phturnover(p,ideadcroot_st)+matrix_gmturnover(p,ideadcroot_st)+matrix_fiturnover(p,ideadcroot_st)) & 
                                                   * deadcrootc_storage(p)
               matrix_cturnover_deadcrootxf_acc(p) = matrix_cturnover_deadcrootxf_acc(p) &
                                                   + (matrix_phturnover(p,ideadcroot_xf)+matrix_gmturnover(p,ideadcroot_xf)+matrix_fiturnover(p,ideadcroot_xf)) &
                                                   * deadcrootc_xfer(p)
            end do
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  matrix_cturnover_grain_acc(p)    = matrix_cturnover_grain_acc(p) &
                                                   + (matrix_phturnover(p,igrain)+matrix_gmturnover(p,igrain)+matrix_fiturnover(p,igrain)) &
                                                   * grainc(p)
                  matrix_cturnover_grainst_acc(p)  = matrix_cturnover_grainst_acc(p) &
                                                   + (matrix_phturnover(p,igrain_st)+matrix_gmturnover(p,igrain_st)+matrix_fiturnover(p,igrain_st)) &
                                                   * grainc_storage(p)
                  matrix_cturnover_grainxf_acc(p)  = matrix_cturnover_grainxf_acc(p) &
                                                   + (matrix_phturnover(p,igrain_xf)+matrix_gmturnover(p,igrain_xf)+matrix_fiturnover(p,igrain_xf)) &
                                                   * grainc_xfer(p)
               end if
            end do

         !print*,'here9'

            do fp = 1,num_soilp
               p = filter_soilp(fp)
               matrix_nalloc_leaf_acc(p)        = matrix_nalloc_leaf_acc(p)        + vegmatrixn_input%V(p,ileaf)
               matrix_nalloc_leafst_acc(p)      = matrix_nalloc_leafst_acc(p)      + vegmatrixn_input%V(p,ileaf_st)
               matrix_nalloc_froot_acc(p)       = matrix_nalloc_froot_acc(p)       + vegmatrixn_input%V(p,ifroot)
               matrix_nalloc_frootst_acc(p)     = matrix_nalloc_frootst_acc(p)     + vegmatrixn_input%V(p,ifroot_st)
               matrix_nalloc_livestem_acc(p)    = matrix_nalloc_livestem_acc(p)    + vegmatrixn_input%V(p,ilivestem)
               matrix_nalloc_livestemst_acc(p)  = matrix_nalloc_livestemst_acc(p)  + vegmatrixn_input%V(p,ilivestem_st)
               matrix_nalloc_deadstem_acc(p)    = matrix_nalloc_deadstem_acc(p)    + vegmatrixn_input%V(p,ideadstem)
               matrix_nalloc_deadstemst_acc(p)  = matrix_nalloc_deadstemst_acc(p)  + vegmatrixn_input%V(p,ideadstem_st)
               matrix_nalloc_livecroot_acc(p)   = matrix_nalloc_livecroot_acc(p)   + vegmatrixn_input%V(p,ilivecroot)
               matrix_nalloc_livecrootst_acc(p) = matrix_nalloc_livecrootst_acc(p) + vegmatrixn_input%V(p,ilivecroot_st)
               matrix_nalloc_deadcroot_acc(p)   = matrix_nalloc_deadcroot_acc(p)   + vegmatrixn_input%V(p,ideadcroot)
               matrix_nalloc_deadcrootst_acc(p) = matrix_nalloc_deadcrootst_acc(p) + vegmatrixn_input%V(p,ideadcroot_st)
            end do

            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  matrix_nalloc_grain_acc(p)    = matrix_nalloc_grain_acc(p)       + vegmatrixn_input%V(p,igrain)
                  matrix_nalloc_grainst_acc(p)  = matrix_nalloc_grainst_acc(p)     + vegmatrixn_input%V(p,igrain_st)
               end if
            end do

            do fp = 1,num_soilp
               p = filter_soilp(fp)
               matrix_ntransfer_leafst_to_leafxf_acc(p)           = matrix_ntransfer_leafst_to_leafxf_acc(p) &
                                                                  + matrix_nphtransfer(p,ileafst_to_ileafxf_phn) &
                                                                  * dt * leafn_storage(p)!matrix_nphturnover(p,ileaf_st)*leafn_storage(p)
               matrix_ntransfer_leafxf_to_leaf_acc(p)             = matrix_ntransfer_leafxf_to_leaf_acc(p) &
                                                                  + matrix_nphtransfer(p,ileafxf_to_ileaf_phn) & 
                                                                  * dt * leafn_xfer(p)!matrix_nphturnover(p,ileaf_xf)*leafn_xfer(p)
               matrix_ntransfer_frootst_to_frootxf_acc(p)         = matrix_ntransfer_frootst_to_frootxf_acc(p) &
                                                                  + matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) & 
                                                                  * dt * frootn_storage(p)!matrix_nphturnover(p,ifroot_st)*frootn_storage(p)
               matrix_ntransfer_frootxf_to_froot_acc(p)           = matrix_ntransfer_frootxf_to_froot_acc(p) &
                                                                  + matrix_nphtransfer(p,ifrootxf_to_ifroot_phn) & 
                                                                  * dt * frootn_xfer(p)!matrix_nphturnover(p,ifroot_xf)*frootn_xfer(p)
               matrix_ntransfer_livestemst_to_livestemxf_acc(p)   = matrix_ntransfer_livestemst_to_livestemxf_acc(p) &
                                                                  + matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) & 
                                                                  * dt * livestemn_storage(p)!matrix_nphturnover(p,ilivestem_st)*livestemn_storage(p)
               matrix_ntransfer_livestemxf_to_livestem_acc(p)     = matrix_ntransfer_livestemxf_to_livestem_acc(p) &
                                                                  + matrix_nphtransfer(p,ilivestemxf_to_ilivestem_phn) & 
                                                                  * dt * livestemn_xfer(p)!matrix_nphturnover(p,ilivestem_xf)*livestemn_xfer(p)
               matrix_ntransfer_deadstemst_to_deadstemxf_acc(p)   = matrix_ntransfer_deadstemst_to_deadstemxf_acc(p) &
                                                                  + matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) & 
                                                                  * dt * deadstemn_storage(p)!matrix_nphturnover(p,ideadstem_st)*deadstemn_storage(p)
               matrix_ntransfer_deadstemxf_to_deadstem_acc(p)     = matrix_ntransfer_deadstemxf_to_deadstem_acc(p) &
                                                                  + matrix_nphtransfer(p,ideadstemxf_to_ideadstem_phn) & 
                                                                  * dt * deadstemn_xfer(p)!matrix_nphturnover(p,ideadstem_xf)*deadstemn_storage(p)
               matrix_ntransfer_livecrootst_to_livecrootxf_acc(p) = matrix_ntransfer_livecrootst_to_livecrootxf_acc(p) &
                                                                  + matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn) & 
                                                                  * dt * livecrootn_storage(p)!matrix_nphturnover(p,ilivecroot_st)*livecrootn_storage(p)
               matrix_ntransfer_livecrootxf_to_livecroot_acc(p)   = matrix_ntransfer_livecrootxf_to_livecroot_acc(p) &
                                                                  + matrix_nphtransfer(p,ilivecrootxf_to_ilivecroot_phn) & 
                                                                  * dt * livecrootn_xfer(p)!matrix_nphturnover(p,ilivecroot_xf)*livecrootn_xfer(p)
               matrix_ntransfer_deadcrootst_to_deadcrootxf_acc(p) = matrix_ntransfer_deadcrootst_to_deadcrootxf_acc(p) &
                                                                  + matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn) &
                                                                  * dt * deadcrootn_storage(p)!matrix_nphturnover(p,ideadcroot_st)*deadcrootn_storage(p)
               matrix_ntransfer_deadcrootxf_to_deadcroot_acc(p)   = matrix_ntransfer_deadcrootxf_to_deadcroot_acc(p) &
                                                                  + matrix_nphtransfer(p,ideadcrootxf_to_ideadcroot_phn) & 
                                                                  * dt * deadcrootn_xfer(p)!matrix_nphturnover(p,ideadcroot_st)*deadcrootn_xfer(p)
               matrix_ntransfer_livestem_to_deadstem_acc(p)       = matrix_ntransfer_livestem_to_deadstem_acc(p) &
                                                                  +(matrix_nphtransfer(p,ilivestem_to_ideadstem_phn) &!*matrix_nphturnover(p,ilivestem) &
                                                                  + matrix_nfitransfer(p,ilivestem_to_ideadstem_fin)) &!*matrix_nfiturnover(p,ilivestem)) &
                                                                  * dt * livestemn(p)
               matrix_ntransfer_livecroot_to_deadcroot_acc(p)     = matrix_ntransfer_livecroot_to_deadcroot_acc(p) &
                                                                  +(matrix_nphtransfer(p,ilivecroot_to_ideadcroot_phn) &!*matrix_nphturnover(p,ilivecroot) &
                                                                  + matrix_nfitransfer(p,ilivecroot_to_ideadcroot_fin)) &!*matrix_nfiturnover(p,ilivecroot)) &
                                                                  * dt * livecrootn(p)
!            if(ivt(p) >= npcropmin)then
!               matrix_ntransfer_grainst_to_grainxf_acc(p)      = matrix_ntransfer_grainst_to_grainxf_acc(p) &
!                                                               + matrix_phAK_C(igrain_xf,igrain_st)*Xoldvegn%V(igrain_st)
!                                                               + vegmatrixc_transfer(igrain_xf,igrain_st)
!               matrix_ntransfer_grainxf_to_grain_acc(p)        = matrix_ntransfer_grainxf_to_grain_acc(p) &
!                                                               + AK_C(igrain,igrain_xf)*Xoldvegn%V(igrain_xf)
!                                                               + vegmatrixc_transfer(igrain,igrain_xf)
!            end if

               matrix_ntransfer_retransn_to_leaf_acc(p)           = matrix_ntransfer_retransn_to_leaf_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ileaf_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_leafst_acc(p)         = matrix_ntransfer_retransn_to_leafst_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ileafst_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_froot_acc(p)          = matrix_ntransfer_retransn_to_froot_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ifroot_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_frootst_acc(p)        = matrix_ntransfer_retransn_to_frootst_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ifrootst_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_livestem_acc(p)       = matrix_ntransfer_retransn_to_livestem_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ilivestem_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_livestemst_acc(p)     = matrix_ntransfer_retransn_to_livestemst_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ilivestemst_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_deadstem_acc(p)       = matrix_ntransfer_retransn_to_deadstem_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ideadstem_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_deadstemst_acc(p)     = matrix_ntransfer_retransn_to_deadstemst_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ideadstemst_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_livecroot_acc(p)      = matrix_ntransfer_retransn_to_livecroot_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ilivecroot_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_livecrootst_acc(p)    = matrix_ntransfer_retransn_to_livecrootst_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ilivecrootst_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_deadcroot_acc(p)      = matrix_ntransfer_retransn_to_deadcroot_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ideadcroot_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_retransn_to_deadcrootst_acc(p)    = matrix_ntransfer_retransn_to_deadcrootst_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_ideadcrootst_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               matrix_ntransfer_leaf_to_retransn_acc(p)           = matrix_ntransfer_leaf_to_retransn_acc(p) &
                                                                  + matrix_nphtransfer(p,ileaf_to_iretransn_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,ileaf)*leafn(p)
               matrix_ntransfer_froot_to_retransn_acc(p)          = matrix_ntransfer_froot_to_retransn_acc(p) &
                                                                  + matrix_nphtransfer(p,ifroot_to_iretransn_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,ifroot)*frootn(p)
               matrix_ntransfer_livestem_to_retransn_acc(p)       = matrix_ntransfer_livestem_to_retransn_acc(p) &
                                                                  + matrix_nphtransfer(p,ilivestem_to_iretransn_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,ilivestem)*livestemn(p)
               matrix_ntransfer_livecroot_to_retransn_acc(p)      = matrix_ntransfer_livecroot_to_retransn_acc(p) &
                                                                  + matrix_nphtransfer(p,ilivecroot_to_iretransn_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,ilivecroot)*livecrootn(p)
            end do

            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  matrix_ntransfer_retransn_to_grain_acc(p)       = matrix_ntransfer_retransn_to_grain_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_igrain_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
                  matrix_ntransfer_retransn_to_grainst_acc(p)     = matrix_ntransfer_retransn_to_grainst_acc(p) &
                                                                  + matrix_nphtransfer(p,iretransn_to_igrainst_phn) & 
                                                                  * dt * retransn(p)!matrix_nphturnover(p,iretransn)*retransn(p)
               end if
            end do

         !print*,'here10'
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               matrix_nturnover_leaf_acc(p)           = matrix_nturnover_leaf_acc(p) &
                                                      + (matrix_nphturnover(p,ileaf)+matrix_ngmturnover(p,ileaf)+matrix_nfiturnover(p,ileaf)) & 
                                                      * leafn(p)
               matrix_nturnover_leafst_acc(p)         = matrix_nturnover_leafst_acc(p) &
                                                      + (matrix_nphturnover(p,ileaf_st)+matrix_ngmturnover(p,ileaf_st)+matrix_nfiturnover(p,ileaf_st)) & 
                                                      * leafn_storage(p)
               matrix_nturnover_leafxf_acc(p)         = matrix_nturnover_leafxf_acc(p) &
                                                      + (matrix_nphturnover(p,ileaf_xf)+matrix_ngmturnover(p,ileaf_xf)+matrix_nfiturnover(p,ileaf_xf)) &
                                                      * leafn_xfer(p)
               matrix_nturnover_froot_acc(p)          = matrix_nturnover_froot_acc(p) &
                                                      + (matrix_nphturnover(p,ifroot)+matrix_ngmturnover(p,ifroot)+matrix_nfiturnover(p,ifroot)) &
                                                      * frootn(p)
               matrix_nturnover_frootst_acc(p)        = matrix_nturnover_frootst_acc(p) &
                                                      + (matrix_nphturnover(p,ifroot_st)+matrix_ngmturnover(p,ifroot_st)+matrix_nfiturnover(p,ifroot_st)) &
                                                      * frootn_storage(p)
               matrix_nturnover_frootxf_acc(p)        = matrix_nturnover_frootxf_acc(p) &
                                                      + (matrix_nphturnover(p,ifroot_xf)+matrix_ngmturnover(p,ifroot_xf)+matrix_nfiturnover(p,ifroot_xf)) &
                                                      * frootn_xfer(p)
               matrix_nturnover_livestem_acc(p)       = matrix_nturnover_livestem_acc(p) &
                                                      + (matrix_nphturnover(p,ilivestem)+matrix_ngmturnover(p,ilivestem)+matrix_nfiturnover(p,ilivestem)) &
                                                      * livestemn(p)
               matrix_nturnover_livestemst_acc(p)     = matrix_nturnover_livestemst_acc(p) &
                                                      + (matrix_nphturnover(p,ilivestem_st)+matrix_ngmturnover(p,ilivestem_st)+matrix_nfiturnover(p,ilivestem_st)) & 
                                                      * livestemn_storage(p)
               matrix_nturnover_livestemxf_acc(p)     = matrix_nturnover_livestemxf_acc(p) &
                                                      + (matrix_nphturnover(p,ilivestem_xf)+matrix_ngmturnover(p,ilivestem_xf)+matrix_nfiturnover(p,ilivestem_xf)) & 
                                                      * livestemn_xfer(p)
               matrix_nturnover_deadstem_acc(p)       = matrix_nturnover_deadstem_acc(p) &
                                                      + (matrix_nphturnover(p,ideadstem)+matrix_ngmturnover(p,ideadstem)+matrix_nfiturnover(p,ideadstem)) &
                                                      * deadstemn(p)
               matrix_nturnover_deadstemst_acc(p)     = matrix_nturnover_deadstemst_acc(p) &
                                                      + (matrix_nphturnover(p,ideadstem_st)+matrix_ngmturnover(p,ideadstem_st)+matrix_nfiturnover(p,ideadstem_st)) &
                                                      * deadstemn_storage(p)
               matrix_nturnover_deadstemxf_acc(p)     = matrix_nturnover_deadstemxf_acc(p) &
                                                      + (matrix_nphturnover(p,ideadstem_xf)+matrix_ngmturnover(p,ideadstem_xf)+matrix_nfiturnover(p,ideadstem_xf)) & 
                                                      * deadstemn_storage(p)
               matrix_nturnover_livecroot_acc(p)      = matrix_nturnover_livecroot_acc(p) &
                                                      + (matrix_nphturnover(p,ilivecroot)+matrix_ngmturnover(p,ilivecroot)+matrix_nfiturnover(p,ilivecroot)) &
                                                      * livecrootn(p)
               matrix_nturnover_livecrootst_acc(p)    = matrix_nturnover_livecrootst_acc(p) &
                                                      + (matrix_nphturnover(p,ilivecroot_st)+matrix_ngmturnover(p,ilivecroot_st)+matrix_nfiturnover(p,ilivecroot_st)) &
                                                      * livecrootn_storage(p)
               matrix_nturnover_livecrootxf_acc(p)    = matrix_nturnover_livecrootxf_acc(p) &
                                                      + (matrix_nphturnover(p,ilivecroot_xf)+matrix_ngmturnover(p,ilivecroot_xf)+matrix_nfiturnover(p,ilivecroot_xf)) &
                                                      * livecrootn_xfer(p)
               matrix_nturnover_deadcroot_acc(p)      = matrix_nturnover_deadcroot_acc(p) &
                                                      + (matrix_nphturnover(p,ideadcroot)+matrix_ngmturnover(p,ideadcroot)+matrix_nfiturnover(p,ideadcroot)) &
                                                      * deadcrootn(p)
               matrix_nturnover_deadcrootst_acc(p)    = matrix_nturnover_deadcrootst_acc(p) &
                                                      + (matrix_nphturnover(p,ideadcroot_st)+matrix_ngmturnover(p,ideadcroot_st)+matrix_nfiturnover(p,ideadcroot_st)) &
                                                      * deadcrootn_storage(p)
               matrix_nturnover_deadcrootxf_acc(p)    = matrix_nturnover_deadcrootxf_acc(p) &
                                                      + (matrix_nphturnover(p,ideadcroot_xf)+matrix_ngmturnover(p,ideadcroot_xf)+matrix_nfiturnover(p,ideadcroot_xf)) &
                                                      * deadcrootn_xfer(p)
               matrix_nturnover_retransn_acc(p)       = matrix_nturnover_retransn_acc(p) &
                                                      + (matrix_nphturnover(p,iretransn)+matrix_ngmturnover(p,iretransn)+matrix_nfiturnover(p,iretransn)) &
                                                      * retransn(p)
            end do
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  matrix_nturnover_grain_acc(p)       = matrix_nturnover_grain_acc(p) &
                                                      + (matrix_nphturnover(p,igrain)+matrix_ngmturnover(p,igrain)+matrix_nfiturnover(p,igrain)) !&
!                                                      * grainn(p)
                  matrix_nturnover_grainst_acc(p)     = matrix_nturnover_grainst_acc(p) &
                                                      + (matrix_nphturnover(p,igrain_st)+matrix_ngmturnover(p,igrain_st)+matrix_nfiturnover(p,igrain_st)) !&
!                                                      * grainn_storage(p)
                  matrix_nturnover_grainxf_acc(p)     = matrix_nturnover_grainxf_acc(p) &
                                                      + (matrix_nphturnover(p,igrain_xf)+matrix_ngmturnover(p,igrain_xf)+matrix_nfiturnover(p,igrain_xf)) !& 
!                                                      * grainn_xfer(p)
               end if
            end do
         end if
         call t_stopf('CN veg matrix-accum. trans.')

         call t_startf('CN veg matrix-assign new value')
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            leafc(p)               = Xvegc%V(p,ileaf)
            leafc_storage(p)       = Xvegc%V(p,ileaf_st)
            leafc_xfer(p)          = Xvegc%V(p,ileaf_xf)
            frootc(p)              = Xvegc%V(p,ifroot)
            frootc_storage(p)      = Xvegc%V(p,ifroot_st)
            frootc_xfer(p)         = Xvegc%V(p,ifroot_xf)
            livestemc(p)           = Xvegc%V(p,ilivestem)
            livestemc_storage(p)   = Xvegc%V(p,ilivestem_st)
            livestemc_xfer(p)      = Xvegc%V(p,ilivestem_xf)
            deadstemc(p)           = Xvegc%V(p,ideadstem)
            deadstemc_storage(p)   = Xvegc%V(p,ideadstem_st)
            deadstemc_xfer(p)      = Xvegc%V(p,ideadstem_xf)
            livecrootc(p)          = Xvegc%V(p,ilivecroot)
            livecrootc_storage(p)  = Xvegc%V(p,ilivecroot_st)
            livecrootc_xfer(p)     = Xvegc%V(p,ilivecroot_xf)
            deadcrootc(p)          = Xvegc%V(p,ideadcroot)
            deadcrootc_storage(p)  = Xvegc%V(p,ideadcroot_st)
            deadcrootc_xfer(p)     = Xvegc%V(p,ideadcroot_xf)
         end do

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            if(ivt(p) >= npcropmin)then
               grainc(p)              = Xvegc%V(p,igrain)
               grainc_storage(p)      = Xvegc%V(p,igrain_st)
               grainc_xfer(p)         = Xvegc%V(p,igrain_xf)
            end if
         end do
         
!KO         if (is_end_curr_year())then
!KO            tempdump(1)  = leafc(p)
!KO            tempdump(2)  = leafc_storage(p)
!KO            tempdump(3)  = leafc_xfer(p)
!KO            tempdump(4)  = frootc(p)
!KO            tempdump(5)  = frootc_storage(p)
!KO            tempdump(6)  = frootc_xfer(p)
!KO            tempdump(7)  = livestemc(p)
!KO            tempdump(8)  = livestemc_storage(p)
!KO            tempdump(9)  = livestemc_xfer(p)
!KO            tempdump(10) = deadstemc(p)
!KO            tempdump(11) = deadstemc_storage(p)
!KO            tempdump(12) = deadstemc_xfer(p)
!KO            tempdump(13) = livecrootc(p)
!KO            tempdump(14) = livecrootc_storage(p)
!KO            tempdump(15) = livecrootc_xfer(p)
!KO            tempdump(16) = deadcrootc(p)
!KO            tempdump(17) = deadcrootc_storage(p)
!KO            tempdump(18) = deadcrootc_xfer(p)
!KO            where(tempdump .lt. 1.e-8)
!KO               tempdump = 1.e-8
!KO            end where
            
!!!!!->           write(bounds%begp+1105000000,"(I,18E17.9)"),p,(tempdump(i),i=1,18)
!            write(bounds%begp+1105000000,"(I,18E17.9)"),p,leafc(p),leafc_storage(p),leafc_xfer(p),frootc(p),frootc_storage(p),frootc_xfer(p),&
!            livestemc(p),livestemc_storage(p),livestemc_xfer(p),deadstemc(p),deadstemc_storage(p),deadstemc_xfer(p),&
!            livecrootc(p),livecrootc_storage(p),livecrootc_xfer(p),deadcrootc(p),deadcrootc_storage(p),deadcrootc_xfer(p)
!KO         end if

         if ( use_c13 ) then
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               cs13_veg%leafc_patch(p)               = Xveg13c%V(p,ileaf)
               cs13_veg%leafc_storage_patch(p)       = Xveg13c%V(p,ileaf_st)
               cs13_veg%leafc_xfer_patch(p)          = Xveg13c%V(p,ileaf_xf)
               cs13_veg%frootc_patch(p)              = Xveg13c%V(p,ifroot)
               cs13_veg%frootc_storage_patch(p)      = Xveg13c%V(p,ifroot_st)
               cs13_veg%frootc_xfer_patch(p)         = Xveg13c%V(p,ifroot_xf)
               cs13_veg%livestemc_patch(p)           = Xveg13c%V(p,ilivestem)
               cs13_veg%livestemc_storage_patch(p)   = Xveg13c%V(p,ilivestem_st)
               cs13_veg%livestemc_xfer_patch(p)      = Xveg13c%V(p,ilivestem_xf)
               cs13_veg%deadstemc_patch(p)           = Xveg13c%V(p,ideadstem)
               cs13_veg%deadstemc_storage_patch(p)   = Xveg13c%V(p,ideadstem_st)
               cs13_veg%deadstemc_xfer_patch(p)      = Xveg13c%V(p,ideadstem_xf)
               cs13_veg%livecrootc_patch(p)          = Xveg13c%V(p,ilivecroot)
               cs13_veg%livecrootc_storage_patch(p)  = Xveg13c%V(p,ilivecroot_st)
               cs13_veg%livecrootc_xfer_patch(p)     = Xveg13c%V(p,ilivecroot_xf)
               cs13_veg%deadcrootc_patch(p)          = Xveg13c%V(p,ideadcroot)
               cs13_veg%deadcrootc_storage_patch(p)  = Xveg13c%V(p,ideadcroot_st)
               cs13_veg%deadcrootc_xfer_patch(p)     = Xveg13c%V(p,ideadcroot_xf)
           end do
        
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  cs13_veg%grainc_patch(p)              = Xveg13c%V(p,igrain)
                  cs13_veg%grainc_storage_patch(p)      = Xveg13c%V(p,igrain_st)
                  cs13_veg%grainc_xfer_patch(p)         = Xveg13c%V(p,igrain_xf)
               end if
            end do
         end if   
         
         if ( use_c14 ) then
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               cs14_veg%leafc_patch(p)               = Xveg14c%V(p,ileaf)
               cs14_veg%leafc_storage_patch(p)       = Xveg14c%V(p,ileaf_st)
               cs14_veg%leafc_xfer_patch(p)          = Xveg14c%V(p,ileaf_xf)
               cs14_veg%frootc_patch(p)              = Xveg14c%V(p,ifroot)
               cs14_veg%frootc_storage_patch(p)      = Xveg14c%V(p,ifroot_st)
               cs14_veg%frootc_xfer_patch(p)         = Xveg14c%V(p,ifroot_xf)
               cs14_veg%livestemc_patch(p)           = Xveg14c%V(p,ilivestem)
               cs14_veg%livestemc_storage_patch(p)   = Xveg14c%V(p,ilivestem_st)
               cs14_veg%livestemc_xfer_patch(p)      = Xveg14c%V(p,ilivestem_xf)
               cs14_veg%deadstemc_patch(p)           = Xveg14c%V(p,ideadstem)
               cs14_veg%deadstemc_storage_patch(p)   = Xveg14c%V(p,ideadstem_st)
               cs14_veg%deadstemc_xfer_patch(p)      = Xveg14c%V(p,ideadstem_xf)
               cs14_veg%livecrootc_patch(p)          = Xveg14c%V(p,ilivecroot)
               cs14_veg%livecrootc_storage_patch(p)  = Xveg14c%V(p,ilivecroot_st)
               cs14_veg%livecrootc_xfer_patch(p)     = Xveg14c%V(p,ilivecroot_xf)
               cs14_veg%deadcrootc_patch(p)          = Xveg14c%V(p,ideadcroot)
               cs14_veg%deadcrootc_storage_patch(p)  = Xveg14c%V(p,ideadcroot_st)
               cs14_veg%deadcrootc_xfer_patch(p)     = Xveg14c%V(p,ideadcroot_xf)
            end do
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  cs14_veg%grainc_patch(p)              = Xveg14c%V(p,igrain)
                  cs14_veg%grainc_storage_patch(p)      = Xveg14c%V(p,igrain_st)
                  cs14_veg%grainc_xfer_patch(p)         = Xveg14c%V(p,igrain_xf)
               end if
            end do
         end if
 
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            leafn(p)               = Xvegn%V(p,ileaf)
            leafn_storage(p)       = Xvegn%V(p,ileaf_st)
            leafn_xfer(p)          = Xvegn%V(p,ileaf_xf)
            frootn(p)              = Xvegn%V(p,ifroot)
            frootn_storage(p)      = Xvegn%V(p,ifroot_st)
            frootn_xfer(p)         = Xvegn%V(p,ifroot_xf)
            livestemn(p)           = Xvegn%V(p,ilivestem)
            livestemn_storage(p)   = Xvegn%V(p,ilivestem_st)
            livestemn_xfer(p)      = Xvegn%V(p,ilivestem_xf)
            deadstemn(p)           = Xvegn%V(p,ideadstem)
            deadstemn_storage(p)   = Xvegn%V(p,ideadstem_st)
            deadstemn_xfer(p)      = Xvegn%V(p,ideadstem_xf)
            livecrootn(p)          = Xvegn%V(p,ilivecroot)
            livecrootn_storage(p)  = Xvegn%V(p,ilivecroot_st)
            livecrootn_xfer(p)     = Xvegn%V(p,ilivecroot_xf)
            deadcrootn(p)          = Xvegn%V(p,ideadcroot)
            deadcrootn_storage(p)  = Xvegn%V(p,ideadcroot_st)
            deadcrootn_xfer(p)     = Xvegn%V(p,ideadcroot_xf)
            retransn(p)            = Xvegn%V(p,iretransn)
         end do

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            if(ivt(p) >= npcropmin)then
               grainn(p)              = Xvegn%V(p,igrain)
               grainn_storage(p)      = Xvegn%V(p,igrain_st)
               grainn_xfer(p)         = Xvegn%V(p,igrain_xf)
            end if
         end do
         call t_stopf('CN veg matrix-assign new value')
         !print*,'here11'
         if(isspinup .or. is_outmatrix)then
            if(is_end_curr_year())then
               do fp = 1,num_soilp
                  call t_startf('CN veg matrix-prepare AK^-1')
                  p = filter_soilp(fp)
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

!KO               where(matrix_calloc_acc .lt. 1.e-8)
!KO                  tempdump(1:18) = 0
!KO               else where
!KO                  tempdump(1:18) = matrix_calloc_acc(1:18)
!KO               end where

!               write(bounds%begp+1101000000,"(I,18E17.9)"),p,matrix_calloc_acc
!!!!!->               write(bounds%begp+1101000000,"(I,18E17.9)"),p,(tempdump(i),i=1,18)

                 !print*,'before ctransfer',p,matrix_ctransfer_leafxf_to_leaf_acc(p)
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
                  
   !KO               tempdump(1)  = matrix_ctransfer_acc(ileaf_xf,ileaf_st)
   !KO               tempdump(2)  = matrix_ctransfer_acc(ileaf,ileaf_xf)
   !KO               tempdump(3)  = matrix_ctransfer_acc(ifroot_xf,ifroot_st)
!KO               tempdump(4)  = matrix_ctransfer_acc(ifroot,ifroot_xf)
!KO               tempdump(5)  = matrix_ctransfer_acc(ilivestem_xf,ilivestem_st)
!KO               tempdump(6)  = matrix_ctransfer_acc(ilivestem,ilivestem_xf)
!KO               tempdump(7)  = matrix_ctransfer_acc(ideadstem_xf,ideadstem_st)
!KO               tempdump(8)  = matrix_ctransfer_acc(ideadstem,ideadstem_xf)
!KO               tempdump(9)  = matrix_ctransfer_acc(ilivecroot_xf,ilivecroot_st)
!KO               tempdump(10) = matrix_ctransfer_acc(ilivecroot,ilivecroot_xf)
!KO               tempdump(11) = matrix_ctransfer_acc(ideadcroot_xf,ideadcroot_st)
!KO               tempdump(12) = matrix_ctransfer_acc(ideadcroot,ideadcroot_xf)
!KO               if(abs(matrix_cturnover_livestem_acc(p))-abs(matrix_cturnover_gm_livestem_acc(p))-abs(matrix_cturnover_fire_livestem_acc(p)) .lt. 1.e-8)then
!KO                  tempdump(13) = (matrix_ctransfer_acc(ideadstem,ilivestem) - matrix_ctransfer_fire_livestem_to_deadstem_acc(p))/1.e-8
!KO               else
!KO                  tempdump(13) = (matrix_ctransfer_acc(ideadstem,ilivestem) - matrix_ctransfer_fire_livestem_to_deadstem_acc(p)) &
!KO                  / (abs(matrix_cturnover_livestem_acc(p))-abs(matrix_cturnover_gm_livestem_acc(p))-abs(matrix_cturnover_fire_livestem_acc(p)))
!KO               end if
!KO               if(abs(matrix_cturnover_livecroot_acc(p))-abs(matrix_cturnover_gm_livecroot_acc(p))-abs(matrix_cturnover_fire_livecroot_acc(p)) .lt. 1.e-8)then
!KO                  tempdump(14) = (matrix_ctransfer_acc(ideadcroot,ilivecroot) - matrix_ctransfer_fire_livecroot_to_deadcroot_acc(p))/1.e-8
!KO               else
!KO                  tempdump(14) = (matrix_ctransfer_acc(ideadcroot,ilivecroot) - matrix_ctransfer_fire_livecroot_to_deadcroot_acc(p)) &
!KO                  / (abs(matrix_cturnover_livecroot_acc(p))-abs(matrix_cturnover_gm_livecroot_acc(p))-abs(matrix_cturnover_fire_livecroot_acc(p)))
!KO               end if
!KO               if(abs(matrix_cturnover_fire_livestem_acc(p)) .lt. 1.e-8)then
!KO                  tempdump(15) = matrix_ctransfer_fire_livestem_to_deadstem_acc(p)/1.e-8
!KO               else
!KO                  tempdump(15) = matrix_ctransfer_fire_livestem_to_deadstem_acc(p)/abs(matrix_cturnover_fire_livestem_acc(p))
!KO               end if
!KO               if(abs(matrix_cturnover_fire_livecroot_acc(p)) .lt. 1.e-8)then
!KO                  tempdump(16) = matrix_ctransfer_fire_livecroot_to_deadcroot_acc(p)/1.e-8
!KO               else
!KO                  tempdump(16) = matrix_ctransfer_fire_livecroot_to_deadcroot_acc(p)/abs(matrix_cturnover_fire_livecroot_acc(p))
!KO               end if

!KO               where(abs(tempdump) .lt. 1.e-8)
!KO                  tempdump = 0
!KO               end where

!!!!!->               write(bounds%begp+1102000000,"(I,16E17.9)"),p,(tempdump(i),i=1,16)

!               write(bounds%begp+1102000000,"(I,16E17.9)"),p,matrix_ctransfer_acc(ileaf_xf,ileaf_st),matrix_ctransfer_acc(ileaf,ileaf_xf),&
!                               matrix_ctransfer_acc(ifroot_xf,ifroot_st),matrix_ctransfer_acc(ifroot,ifroot_xf),&
!                               matrix_ctransfer_acc(ilivestem_xf,ilivestem_st),matrix_ctransfer_acc(ilivestem,ilivestem_xf),&
!                               matrix_ctransfer_acc(ideadstem_xf,ideadstem_st),matrix_ctransfer_acc(ideadstem,ideadstem_xf),&
!                               matrix_ctransfer_acc(ilivecroot_xf,ilivecroot_st),matrix_ctransfer_acc(ilivecroot,ilivecroot_xf),&
!                               matrix_ctransfer_acc(ideadcroot_xf,ideadcroot_st),matrix_ctransfer_acc(ideadcroot,ideadcroot_xf),&
!                               matrix_ctransfer_acc(ideadstem,ilivestem), matrix_ctransfer_acc(ideadcroot,ilivecroot),&
!                               matrix_ctransfer_fire_livestem_to_deadstem_acc(p),matrix_ctransfer_fire_livecroot_to_deadcroot_acc(p)

                 !print*,'before cturnover',p,matrix_cturnover_leaf_acc(p),matrix_cturnover_leafxf_acc(p)
                  matrix_ctransfer_acc(ileaf,ileaf)                 = -matrix_cturnover_leaf_acc(p)               
                  matrix_ctransfer_acc(ileaf_st,ileaf_st)           = -matrix_cturnover_leafst_acc(p)               
                  matrix_ctransfer_acc(ileaf_xf,ileaf_xf)           = -matrix_cturnover_leafxf_acc(p)               
                  matrix_ctransfer_acc(ifroot,ifroot)               = -matrix_cturnover_froot_acc(p)               
                  matrix_ctransfer_acc(ifroot_st,ifroot_st)         = -matrix_cturnover_frootst_acc(p)               
                  matrix_ctransfer_acc(ifroot_xf,ifroot_xf)         = -matrix_cturnover_frootxf_acc(p)               
                  matrix_ctransfer_acc(ilivestem,ilivestem)         = -matrix_cturnover_livestem_acc(p)               
                  matrix_ctransfer_acc(ilivestem_st,ilivestem_st)   = -matrix_cturnover_livestemst_acc(p)               
                  matrix_ctransfer_acc(ilivestem_xf,ilivestem_xf)   = -matrix_cturnover_livestemxf_acc(p)               
                  matrix_ctransfer_acc(ideadstem,ideadstem)         = -matrix_cturnover_deadstem_acc(p)               
                  matrix_ctransfer_acc(ideadstem_st,ideadstem_st)   = -matrix_cturnover_deadstemst_acc(p)               
                  matrix_ctransfer_acc(ideadstem_xf,ideadstem_xf)   = -matrix_cturnover_deadstemxf_acc(p)               
                  matrix_ctransfer_acc(ilivecroot,ilivecroot)       = -matrix_cturnover_livecroot_acc(p)               
                  matrix_ctransfer_acc(ilivecroot_st,ilivecroot_st) = -matrix_cturnover_livecrootst_acc(p)               
                  matrix_ctransfer_acc(ilivecroot_xf,ilivecroot_xf) = -matrix_cturnover_livecrootxf_acc(p)               
                  matrix_ctransfer_acc(ideadcroot,ideadcroot)       = -matrix_cturnover_deadcroot_acc(p)               
                  matrix_ctransfer_acc(ideadcroot_st,ideadcroot_st) = -matrix_cturnover_deadcrootst_acc(p)               
                  matrix_ctransfer_acc(ideadcroot_xf,ideadcroot_xf) = -matrix_cturnover_deadcrootxf_acc(p)               
                  if(ivt(p) >= npcropmin)then
                     matrix_ctransfer_acc(igrain,igrain)            = -matrix_cturnover_grain_acc(p)               
                     matrix_ctransfer_acc(igrain_st,igrain_st)      = -matrix_cturnover_grainst_acc(p)               
                     matrix_ctransfer_acc(igrain_xf,igrain_xf)      = -matrix_cturnover_grainxf_acc(p)               
                  end if

!KO            if(leafc0(p) < 1.e-8)then
!                 tempdump(19) = matrix_cturnover_gm_leaf_acc(p)/1.e-8
!                 tempdump(37) = matrix_cturnover_fire_leaf_acc(p)/1.e-8
!                 tempdump(1)  = matrix_cturnover_leaf_acc(p)/1.e-8-tempdump(19)-tempdump(37)
!              else
!                 tempdump(19) = matrix_cturnover_gm_leaf_acc(p)/leafc0(p)
!                 tempdump(37) = matrix_cturnover_fire_leaf_acc(p)/leafc0(p)
!                 tempdump(1)  = matrix_cturnover_leaf_acc(p)/leafc0(p)-tempdump(19)-tempdump(37)
!              end if
!              if(leafc0_storage(p) < 1.e-8)then
!                 tempdump(20)  = matrix_cturnover_gm_leafst_acc(p)/1.e-8
!                 tempdump(38)  = matrix_cturnover_fire_leafst_acc(p)/1.e-8
!                 tempdump(2)   = matrix_cturnover_leafst_acc(p)/1.e-8-tempdump(20)-tempdump(38)
!              else
!                 tempdump(20)  = matrix_cturnover_gm_leafst_acc(p)/leafc0_storage(p)
!                 tempdump(38)  = matrix_cturnover_fire_leafst_acc(p)/leafc0_storage(p)
!                 tempdump(2)   = matrix_cturnover_leafst_acc(p)/leafc0_storage(p)-tempdump(20)-tempdump(38)
!              end if
!              if(leafc0_xfer(p) < 1.e-8)then
!                 tempdump(21)  = matrix_cturnover_gm_leafxf_acc(p)/1.e-8
!                 tempdump(39)  = matrix_cturnover_fire_leafxf_acc(p)/1.e-8
!                 tempdump(3)   = matrix_cturnover_leafxf_acc(p)/1.e-8-tempdump(21)-tempdump(39)
!              else
!                 tempdump(21)  = matrix_cturnover_gm_leafxf_acc(p)/leafc0_xfer(p)
!                 tempdump(39)  = matrix_cturnover_fire_leafxf_acc(p)/leafc0_xfer(p)
!                 tempdump(3)   = matrix_cturnover_leafxf_acc(p)/leafc0_xfer(p)-tempdump(21)-tempdump(39)
!              end if
!              if(frootc0(p) < 1.e-8)then
!                 tempdump(22) = matrix_cturnover_gm_froot_acc(p)/1.e-8
!                 tempdump(40) = matrix_cturnover_fire_froot_acc(p)/1.e-8
!                 tempdump(4)  = matrix_cturnover_froot_acc(p)/1.e-8-tempdump(22)-tempdump(40)
!              else
!                 tempdump(22) = matrix_cturnover_gm_froot_acc(p)/frootc0(p)
!                 tempdump(40) = matrix_cturnover_fire_froot_acc(p)/frootc0(p)
!                 tempdump(4)  = matrix_cturnover_froot_acc(p)/frootc0(p)-tempdump(22)-tempdump(40)
!              end if
!              if(frootc0_storage(p) < 1.e-8)then
!                 tempdump(23)  = matrix_cturnover_gm_frootst_acc(p)/1.e-8
!                 tempdump(41)  = matrix_cturnover_fire_frootst_acc(p)/1.e-8
!                 tempdump(5)   = matrix_cturnover_frootst_acc(p)/1.e-8-tempdump(23)-tempdump(41)
!              else
!                 tempdump(23)  = matrix_cturnover_gm_frootst_acc(p)/frootc0_storage(p)
!                 tempdump(41)  = matrix_cturnover_fire_frootst_acc(p)/frootc0_storage(p)
!                 tempdump(5)   = matrix_cturnover_frootst_acc(p)/frootc0_storage(p)-tempdump(23)-tempdump(41)
!              end if
!              if(frootc0_xfer(p) < 1.e-8)then
!                 tempdump(24)  = matrix_cturnover_gm_frootxf_acc(p)/1.e-8
!                 tempdump(42)  = matrix_cturnover_fire_frootxf_acc(p)/1.e-8
!                 tempdump(6)   = matrix_cturnover_frootxf_acc(p)/1.e-8-tempdump(24)-tempdump(42)
!              else
!                 tempdump(24)  = matrix_cturnover_gm_frootxf_acc(p)/frootc0_xfer(p)
!                 tempdump(42)  = matrix_cturnover_fire_frootxf_acc(p)/frootc0_xfer(p)
!                 tempdump(6)   = matrix_cturnover_frootxf_acc(p)/frootc0_xfer(p)-tempdump(24)-tempdump(42)
!              end if
!              if(livestemc0(p) < 1.e-8)then
!                 tempdump(25) = matrix_cturnover_gm_livestem_acc(p)/1.e-8
!                 tempdump(43) = matrix_cturnover_fire_livestem_acc(p)/1.e-8
!                 tempdump(7)  = matrix_cturnover_livestem_acc(p)/1.e-8-tempdump(25)-tempdump(43)
!              else
!                 tempdump(25) = matrix_cturnover_gm_livestem_acc(p)/livestemc0(p)
!                 tempdump(43) = matrix_cturnover_fire_livestem_acc(p)/livestemc0(p)
!                 tempdump(7)  = matrix_cturnover_livestem_acc(p)/livestemc0(p)-tempdump(25)-tempdump(43)
!              end if
!              if(livestemc0_storage(p) < 1.e-8)then
!                 tempdump(26)  = matrix_cturnover_gm_livestemst_acc(p)/1.e-8
!                 tempdump(44)  = matrix_cturnover_fire_livestemst_acc(p)/1.e-8
!                 tempdump(8)   = matrix_cturnover_livestemst_acc(p)/1.e-8-tempdump(26)-tempdump(44)
!              else
!                 tempdump(26)  = matrix_cturnover_gm_livestemst_acc(p)/livestemc0_storage(p)
!                 tempdump(44)  = matrix_cturnover_fire_livestemst_acc(p)/livestemc0_storage(p)
!                 tempdump(8)   = matrix_cturnover_livestemst_acc(p)/livestemc0_storage(p)-tempdump(26)-tempdump(44)
!              end if
!              if(livestemc0_xfer(p) < 1.e-8)then
!                 tempdump(27)  = matrix_cturnover_gm_livestemxf_acc(p)/1.e-8
!                 tempdump(45)  = matrix_cturnover_fire_livestemxf_acc(p)/1.e-8
!                 tempdump(9)   = matrix_cturnover_livestemxf_acc(p)/1.e-8-tempdump(27)-tempdump(45)
!              else
!                 tempdump(27)  = matrix_cturnover_gm_livestemxf_acc(p)/livestemc0_xfer(p)
!                 tempdump(45)  = matrix_cturnover_fire_livestemxf_acc(p)/livestemc0_xfer(p)
!                 tempdump(9)   = matrix_cturnover_livestemxf_acc(p)/livestemc0_xfer(p)-tempdump(27)-tempdump(45)
!              end if
!              if(deadstemc0(p) < 1.e-8)then
!                 tempdump(28) = matrix_cturnover_gm_deadstem_acc(p)/1.e-8
!                 tempdump(46) = matrix_cturnover_fire_deadstem_acc(p)/1.e-8
!                 tempdump(10) = matrix_cturnover_deadstem_acc(p)/1.e-8-tempdump(28)-tempdump(46)
!              else
!                 tempdump(28) = matrix_cturnover_gm_deadstem_acc(p)/deadstemc0(p)
!                 tempdump(46) = matrix_cturnover_fire_deadstem_acc(p)/deadstemc0(p)
!                 tempdump(10) = matrix_cturnover_deadstem_acc(p)/deadstemc0(p)-tempdump(28)-tempdump(46)
!              end if
!              if(deadstemc0_storage(p) < 1.e-8)then
!                 tempdump(29)  = matrix_cturnover_gm_deadstemst_acc(p)/1.e-8
!                 tempdump(47)  = matrix_cturnover_fire_deadstemst_acc(p)/1.e-8
!                 tempdump(11)  = matrix_cturnover_deadstemst_acc(p)/1.e-8-tempdump(29)-tempdump(47)
!              else
!                 tempdump(29)  = matrix_cturnover_gm_deadstemst_acc(p)/deadstemc0_storage(p)
!                 tempdump(47)  = matrix_cturnover_fire_deadstemst_acc(p)/deadstemc0_storage(p)
!                 tempdump(11)  = matrix_cturnover_deadstemst_acc(p)/deadstemc0_storage(p)-tempdump(29)-tempdump(47)
!              end if
!              if(deadstemc0_xfer(p) < 1.e-8)then
!                 tempdump(30)  = matrix_cturnover_gm_deadstemxf_acc(p)/1.e-8
!                 tempdump(48)  = matrix_cturnover_fire_deadstemxf_acc(p)/1.e-8
!                 tempdump(12)  = matrix_cturnover_deadstemxf_acc(p)/1.e-8-tempdump(30)-tempdump(48)
!              else
!                 tempdump(30)  = matrix_cturnover_gm_deadstemxf_acc(p)/deadstemc0_xfer(p)
!                 tempdump(48)  = matrix_cturnover_fire_deadstemxf_acc(p)/deadstemc0_xfer(p)
!                 tempdump(12)  = matrix_cturnover_deadstemxf_acc(p)/deadstemc0_xfer(p)-tempdump(30)-tempdump(48)
!              end if
!              if(livecrootc0(p) < 1.e-8)then
!                 tempdump(31) = matrix_cturnover_gm_livecroot_acc(p)/1.e-8
!                 tempdump(49) = matrix_cturnover_fire_livecroot_acc(p)/1.e-8
!                 tempdump(13) = matrix_cturnover_livecroot_acc(p)/1.e-8-tempdump(31)-tempdump(49)
!              else
!                 tempdump(31) = matrix_cturnover_gm_livecroot_acc(p)/livecrootc0(p)
!                 tempdump(49) = matrix_cturnover_fire_livecroot_acc(p)/livecrootc0(p)
!                 tempdump(13) = matrix_cturnover_livecroot_acc(p)/livecrootc0(p)-tempdump(31)-tempdump(49)
!              end if
!              if(livecrootc0_storage(p) < 1.e-8)then
!                 tempdump(32)  = matrix_cturnover_gm_livecrootst_acc(p)/1.e-8
!                 tempdump(50)  = matrix_cturnover_fire_livecrootst_acc(p)/1.e-8
!                 tempdump(14)  = matrix_cturnover_livecrootst_acc(p)/1.e-8-tempdump(32)-tempdump(50)
!              else
!                 tempdump(32)  = matrix_cturnover_gm_livecrootst_acc(p)/livecrootc0_storage(p)
!                 tempdump(50)  = matrix_cturnover_fire_livecrootst_acc(p)/livecrootc0_storage(p)
!                 tempdump(14)  = matrix_cturnover_livecrootst_acc(p)/livecrootc0_storage(p)-tempdump(32)-tempdump(50)
!              end if
!              if(livecrootc0_xfer(p) < 1.e-8)then
!                 tempdump(33)  = matrix_cturnover_gm_livecrootxf_acc(p)/1.e-8
!                 tempdump(51)  = matrix_cturnover_fire_livecrootxf_acc(p)/1.e-8
!                 tempdump(15)  = matrix_cturnover_livecrootxf_acc(p)/1.e-8-tempdump(33)-tempdump(51)
!              else
!                 tempdump(33)  = matrix_cturnover_gm_livecrootxf_acc(p)/livecrootc0_xfer(p)
!                 tempdump(51)  = matrix_cturnover_fire_livecrootxf_acc(p)/livecrootc0_xfer(p)
!                 tempdump(15)  = matrix_cturnover_livecrootxf_acc(p)/livecrootc0_xfer(p)-tempdump(33)-tempdump(51)
!              end if
!              if(deadcrootc0(p) < 1.e-8)then
!                 tempdump(34) = matrix_cturnover_gm_deadcroot_acc(p)/1.e-8
!                 tempdump(52) = matrix_cturnover_fire_deadcroot_acc(p)/1.e-8
!                 tempdump(16) = matrix_cturnover_deadcroot_acc(p)/1.e-8-tempdump(34)-tempdump(52)
!              else
!                 tempdump(34) = matrix_cturnover_gm_deadcroot_acc(p)/deadcrootc0(p)
!                 tempdump(52) = matrix_cturnover_fire_deadcroot_acc(p)/deadcrootc0(p)
!                 tempdump(16) = matrix_cturnover_deadcroot_acc(p)/deadcrootc0(p)-tempdump(34)-tempdump(52)
!              end if
!              if(deadcrootc0_storage(p) < 1.e-8)then
!                 tempdump(35)  = matrix_cturnover_gm_deadcrootst_acc(p)/1.e-8
!                 tempdump(53)  = matrix_cturnover_fire_deadcrootst_acc(p)/1.e-8
!                 tempdump(17)  = matrix_cturnover_deadcrootst_acc(p)/1.e-8-tempdump(35)-tempdump(53)
!              else
!                 tempdump(35)  = matrix_cturnover_gm_deadcrootst_acc(p)/deadcrootc0_storage(p)
!                 tempdump(53)  = matrix_cturnover_fire_deadcrootst_acc(p)/deadcrootc0_storage(p)
!                 tempdump(17)  = matrix_cturnover_deadcrootst_acc(p)/deadcrootc0_storage(p)-tempdump(35)-tempdump(53)
!              end if
!              if(deadcrootc0_xfer(p) < 1.e-8)then
!                 tempdump(36)  = matrix_cturnover_gm_deadcrootxf_acc(p)/1.e-8
!                 tempdump(54)  = matrix_cturnover_fire_deadcrootxf_acc(p)/1.e-8
!                 tempdump(18)  = matrix_cturnover_deadcrootxf_acc(p)/1.e-8-tempdump(36)-tempdump(54)
!              else
!                 tempdump(36)  = matrix_cturnover_gm_deadcrootxf_acc(p)/deadcrootc0_xfer(p)
!                 tempdump(54)  = matrix_cturnover_fire_deadcrootxf_acc(p)/deadcrootc0_xfer(p)
!                 tempdump(18)  = matrix_cturnover_deadcrootxf_acc(p)/deadcrootc0_xfer(p)-tempdump(36)-tempdump(54)
!              end if
               
!              where(abs(tempdump) .lt. 1.e-8)
!                 tempdump = 0
!              endwhere

!!!!!->               write(bounds%begp+1103000000,"(I,54E17.9)"),p,(tempdump(i),i=1,54)
!               write(bounds%begp+1103000000,"(I,54E17.9)"),p,matrix_cturnover_leaf_acc(p),matrix_cturnover_leafst_acc(p),matrix_cturnover_leafxf_acc(p),&
!                               matrix_cturnover_froot_acc(p),matrix_cturnover_frootst_acc(p),matrix_cturnover_frootxf_acc(p),&
!                               matrix_cturnover_livestem_acc(p),matrix_cturnover_livestemst_acc(p),matrix_cturnover_livestemxf_acc(p),&
!                               matrix_cturnover_deadstem_acc(p),matrix_cturnover_deadstemst_acc(p),matrix_cturnover_deadstemxf_acc(p),&
!                               matrix_cturnover_livecroot_acc(p),matrix_cturnover_livecrootst_acc(p),matrix_cturnover_livecrootxf_acc(p),&
!                               matrix_cturnover_deadcroot_acc(p),matrix_cturnover_deadcrootst_acc(p),matrix_cturnover_deadcrootxf_acc(p),&
!                               matrix_cturnover_gm_leaf_acc(p),matrix_cturnover_gm_leafst_acc(p),matrix_cturnover_gm_leafxf_acc(p),&
!                               matrix_cturnover_gm_froot_acc(p),matrix_cturnover_gm_frootst_acc(p),matrix_cturnover_gm_frootxf_acc(p),&
!                               matrix_cturnover_gm_livestem_acc(p),matrix_cturnover_gm_livestemst_acc(p),matrix_cturnover_gm_livestemxf_acc(p),&
!                               matrix_cturnover_gm_deadstem_acc(p),matrix_cturnover_gm_deadstemst_acc(p),matrix_cturnover_gm_deadstemxf_acc(p),&
!                               matrix_cturnover_gm_livecroot_acc(p),matrix_cturnover_gm_livecrootst_acc(p),matrix_cturnover_gm_livecrootxf_acc(p),&
!                               matrix_cturnover_gm_deadcroot_acc(p),matrix_cturnover_gm_deadcrootst_acc(p),matrix_cturnover_gm_deadcrootxf_acc(p),&
!                               matrix_cturnover_fire_leaf_acc(p),matrix_cturnover_fire_leafst_acc(p),matrix_cturnover_fire_leafxf_acc(p),&
!                               matrix_cturnover_fire_froot_acc(p),matrix_cturnover_fire_frootst_acc(p),matrix_cturnover_fire_frootxf_acc(p),&
!                               matrix_cturnover_fire_livestem_acc(p),matrix_cturnover_fire_livestemst_acc(p),matrix_cturnover_fire_livestemxf_acc(p),&
!                               matrix_cturnover_fire_deadstem_acc(p),matrix_cturnover_fire_deadstemst_acc(p),matrix_cturnover_fire_deadstemxf_acc(p),&
!                               matrix_cturnover_fire_livecroot_acc(p),matrix_cturnover_fire_livecrootst_acc(p),matrix_cturnover_fire_livecrootxf_acc(p),&
!                               matrix_cturnover_fire_deadcroot_acc(p),matrix_cturnover_fire_deadcrootst_acc(p),matrix_cturnover_fire_deadcrootxf_acc(p)            


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
   
                  matrix_ntransfer_acc(ileaf,ileaf)                 = -matrix_nturnover_leaf_acc(p)               
                  matrix_ntransfer_acc(ileaf_st,ileaf_st)           = -matrix_nturnover_leafst_acc(p)               
                  matrix_ntransfer_acc(ileaf_xf,ileaf_xf)           = -matrix_nturnover_leafxf_acc(p)               
                  matrix_ntransfer_acc(ifroot,ifroot)               = -matrix_nturnover_froot_acc(p)               
                  matrix_ntransfer_acc(ifroot_st,ifroot_st)         = -matrix_nturnover_frootst_acc(p)               
                  matrix_ntransfer_acc(ifroot_xf,ifroot_xf)         = -matrix_nturnover_frootxf_acc(p)               
                  matrix_ntransfer_acc(ilivestem,ilivestem)         = -matrix_nturnover_livestem_acc(p)               
                  matrix_ntransfer_acc(ilivestem_st,ilivestem_st)   = -matrix_nturnover_livestemst_acc(p)               
                  matrix_ntransfer_acc(ilivestem_xf,ilivestem_xf)   = -matrix_nturnover_livestemxf_acc(p)               
                  matrix_ntransfer_acc(ideadstem,ideadstem)         = -matrix_nturnover_deadstem_acc(p)               
                  matrix_ntransfer_acc(ideadstem_st,ideadstem_st)   = -matrix_nturnover_deadstemst_acc(p)               
                  matrix_ntransfer_acc(ideadstem_xf,ideadstem_xf)   = -matrix_nturnover_deadstemxf_acc(p)               
                  matrix_ntransfer_acc(ilivecroot,ilivecroot)       = -matrix_nturnover_livecroot_acc(p)               
                  matrix_ntransfer_acc(ilivecroot_st,ilivecroot_st) = -matrix_nturnover_livecrootst_acc(p)               
                  matrix_ntransfer_acc(ilivecroot_xf,ilivecroot_xf) = -matrix_nturnover_livecrootxf_acc(p)               
                  matrix_ntransfer_acc(ideadcroot,ideadcroot)       = -matrix_nturnover_deadcroot_acc(p)               
                  matrix_ntransfer_acc(ideadcroot_st,ideadcroot_st) = -matrix_nturnover_deadcrootst_acc(p)               
                  matrix_ntransfer_acc(ideadcroot_xf,ideadcroot_xf) = -matrix_nturnover_deadcrootxf_acc(p)               
                  if(ivt(p) >= npcropmin)then
                     matrix_ntransfer_acc(igrain,igrain)            = -matrix_nturnover_grain_acc(p)               
                     matrix_ntransfer_acc(igrain_st,igrain_st)      = -matrix_nturnover_grainst_acc(p)               
                     matrix_ntransfer_acc(igrain_xf,igrain_xf)      = -matrix_nturnover_grainxf_acc(p)               
                  end if
                  matrix_ntransfer_acc(iretransn,iretransn)         = -matrix_nturnover_retransn_acc(p)               

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
                
                  call t_stopf('CN veg matrix-prepare AK^-1')
                  call t_startf('CN veg matrix-inv matrix operation')
                  !print*,'before capc',p,leafc0(p),leafc0_storage(p),leafc0_xfer(p),frootc0(p),frootc0_storage(p),frootc0_xfer(p),livestemc0(p),livestemc0_storage(p),livestemc0_xfer(p),deadstemc0(p),deadstemc0_storage(p),deadstemc0_xfer(p),livecrootc0(p),livecrootc0_storage(p),livecrootc0_xfer(p),deadcrootc0(p),deadcrootc0_storage(p),deadcrootc0_xfer(p)
!                  do i=1,nvegcpool
!                     print*,'transfer',p,matrix_ctransfer_acc(i,1:nvegcpool)
!                  end do
!                  print*,'calloc',matrix_calloc_acc(1:nvegcpool)
                  call inverse(matrix_ctransfer_acc(1:nvegcpool,1:nvegcpool),AKinvc(1:nvegcpool,1:nvegcpool),nvegcpool)
                  vegmatrixc_rt(:) = -matmul(AKinvc(1:nvegcpool,1:nvegcpool),matrix_calloc_acc(1:nvegcpool))
!                  print*,'after capc',p,vegmatrixc_rt(1:nvegcpool)
! N  
!                  print*,'before capn',p,leafn0(p),leafn0_storage(p),leafn0_xfer(p),frootn0(p),frootn0_storage(p),frootn0_xfer(p),livestemn0(p),livestemn0_storage(p),livestemn0_xfer(p),deadstemn0(p),deadstemn0_storage(p),deadstemn0_xfer(p),livecrootn0(p),livecrootn0_storage(p),livecrootn0_xfer(p),deadcrootn0(p),deadcrootn0_storage(p),deadcrootn0_xfer(p)
                  call inverse(matrix_ntransfer_acc(1:nvegnpool,1:nvegnpool),AKinvn(1:nvegnpool,1:nvegnpool),nvegnpool)
                  vegmatrixn_rt(:) = -matmul(AKinvn(1:nvegnpool,1:nvegnpool),matrix_nalloc_acc(1:nvegnpool))
!                  print*,'after capn',p,vegmatrixn_rt(1:nvegnpool)
                  call t_stopf('CN veg matrix-inv matrix operation')
!            do i=1,nvegpool
!            end do
!                                                                                             + (matrix_fitransfer(p,ideadstem,ilivestem)*matrix_fiturnover(p,ilivestem,ilivestem)))   
 
                  call t_startf('CN veg matrix-finalize spinup')
                  if(is_first_step_of_this_run_segment())then
!!!!!->                   write(bounds%begp+1104000000,"(2I,4E17.9)"),p,patch%itype(p),patch%wtgcell(p),grc%latdeg(patch%gridcell(p)),grc%londeg(patch%gridcell(p)),grc%area(patch%gridcell(p))
                  end if
                  
                  if(isspinup .and. .not. is_first_step_of_this_run_segment())then
                     leafc(p)                  = vegmatrixc_rt(ileaf)
                     leafc_storage(p)          = vegmatrixc_rt(ileaf_st)
      !               leafc_xfer(p)             = vegmatrixc_rt(ileaf_xf)
                     frootc(p)                 = vegmatrixc_rt(ifroot)
                     frootc_storage(p)         = vegmatrixc_rt(ifroot_st)
      !               frootc_xfer(p)            = vegmatrixc_rt(ifroot_xf)
                     livestemc(p)              = vegmatrixc_rt(ilivestem)
                     livestemc_storage(p)      = vegmatrixc_rt(ilivestem_st)
      !               livestemc_xfer(p)         = vegmatrixc_rt(ilivestem_xf)
                     deadstemc(p)              = vegmatrixc_rt(ideadstem)
                     deadstemc_storage(p)      = vegmatrixc_rt(ideadstem_st)
      !               deadstemc_xfer(p)         = vegmatrixc_rt(ideadstem_xf)
                     livecrootc(p)              = vegmatrixc_rt(ilivecroot)
                     livecrootc_storage(p)      = vegmatrixc_rt(ilivecroot_st)
      !               livecrootc_xfer(p)         = vegmatrixc_rt(ilivecroot_xf)   
                     deadcrootc(p)              = vegmatrixc_rt(ideadcroot)
                     deadcrootc_storage(p)      = vegmatrixc_rt(ideadcroot_st)
      !               deadcrootc_xfer(p)         = vegmatrixc_rt(ideadcroot_xf) 
                     if(ivt(p) >= npcropmin)then
                        grainc(p)                 = vegmatrixc_rt(igrain)
                        grainc_storage(p)         = vegmatrixc_rt(igrain_st)
                     end if
                     leafn(p)                  = vegmatrixn_rt(ileaf)
                     leafn_storage(p)          = vegmatrixn_rt(ileaf_st)
      !               leafn_xfer(p)             = vegmatrixn_rt(ileaf_xf)
                     frootn(p)                 = vegmatrixn_rt(ifroot)
                     frootn_storage(p)         = vegmatrixn_rt(ifroot_st)
      !               frootn_xfer(p)            = vegmatrixn_rt(ifroot_xf)
                     livestemn(p)              = vegmatrixn_rt(ilivestem)
                     livestemn_storage(p)      = vegmatrixn_rt(ilivestem_st)
      !               livestemn_xfer(p)         = vegmatrixn_rt(ilivestem_xf)
                     deadstemn(p)              = vegmatrixn_rt(ideadstem)
                     deadstemn_storage(p)      = vegmatrixn_rt(ideadstem_st)
      !               deadstemn_xfer(p)         = vegmatrixn_rt(ideadstem_xf)
                     livecrootn(p)              = vegmatrixn_rt(ilivecroot)
                     livecrootn_storage(p)      = vegmatrixn_rt(ilivecroot_st)
      !               livecrootn_xfer(p)         = vegmatrixn_rt(ilivecroot_xf)   
                     deadcrootn(p)              = vegmatrixn_rt(ideadcroot)
                     deadcrootn_storage(p)      = vegmatrixn_rt(ideadcroot_st)
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
!                     matrix_pot_leafc(p)                  = matrix_cap_leafc(p) - leafc(p)
!                     matrix_pot_leafc_storage(p)          = matrix_cap_leafc_storage(p) - leafc_storage(p)
!                     matrix_pot_leafc_xfer(p)             = matrix_cap_leafc_storage(p) - leafc_xfer(p)
!                     matrix_pot_frootc(p)                 = matrix_cap_frootc(p) - frootc(p)
!                     matrix_pot_frootc_storage(p)         = matrix_cap_frootc_storage(p) - frootc_storage(p)
!                     matrix_pot_frootc_xfer(p)            = matrix_cap_frootc_storage(p) - frootc_xfer(p)
!                     matrix_pot_livestemc(p)              = matrix_cap_livestemc(p) - livestemc(p)
!                     matrix_pot_livestemc_storage(p)      = matrix_cap_livestemc_storage(p) - livestemc_storage(p)
!                     matrix_pot_livestemc_xfer(p)         = matrix_cap_livestemc_storage(p) - livestemc_xfer(p)
!                     matrix_pot_deadstemc(p)              = matrix_cap_deadstemc(p) - deadstemc(p)
!                     matrix_pot_deadstemc_storage(p)      = matrix_cap_deadstemc_storage(p) - deadstemc_storage(p)
!                     matrix_pot_deadstemc_xfer(p)         = matrix_cap_deadstemc_storage(p) - deadstemc_xfer(p)
!                     matrix_pot_livecrootc(p)             = matrix_cap_livecrootc(p) - livecrootc(p)
!                     matrix_pot_livecrootc_storage(p)     = matrix_cap_livecrootc_storage(p) - livecrootc_storage(p)
!                     matrix_pot_livecrootc_xfer(p)        = matrix_cap_livecrootc_storage(p) - livecrootc_xfer(p)
!                     matrix_pot_deadcrootc(p)             = matrix_cap_deadcrootc(p) - deadcrootc(p)
!                     matrix_pot_deadcrootc_storage(p)     = matrix_cap_deadcrootc_storage(p) - deadcrootc_storage(p)
!                     matrix_pot_deadcrootc_xfer(p)        = matrix_cap_deadcrootc_storage(p) - deadcrootc_xfer(p)
!                     if(ivt(p) >= npcropmin)then
!                        matrix_pot_grainc(p)                 = matrix_cap_grainc(p) - grainc(p)
!                        matrix_pot_grainc_storage(p)         = matrix_cap_grainc_storage(p) - grainc_storage(p)
!                        matrix_pot_grainc_xfer(p)            = matrix_cap_grainc_storage(p) - grainc_xfer(p)
!                     end if
!                     matrix_pot_leafn(p)                  = matrix_cap_leafn(p) - leafn(p)
!                     matrix_pot_leafn_storage(p)          = matrix_cap_leafn_storage(p) - leafn_storage(p)
!                     matrix_pot_leafn_xfer(p)             = matrix_cap_leafn_storage(p) - leafn_xfer(p)
!                     matrix_pot_frootn(p)                 = matrix_cap_frootn(p) - frootn(p)
!                     matrix_pot_frootn_storage(p)         = matrix_cap_frootn_storage(p) - frootn_storage(p)
!                     matrix_pot_frootn_xfer(p)            = matrix_cap_frootn_storage(p) - frootn_xfer(p)
!                     matrix_pot_livestemn(p)              = matrix_cap_livestemn(p) - livestemn(p)
!                     matrix_pot_livestemn_storage(p)      = matrix_cap_livestemn_storage(p) - livestemn_storage(p)
!                     matrix_pot_livestemn_xfer(p)         = matrix_cap_livestemn_storage(p) - livestemn_xfer(p)
!                     matrix_pot_deadstemn(p)              = matrix_cap_deadstemn(p) - deadstemn(p)
!                     matrix_pot_deadstemn_storage(p)      = matrix_cap_deadstemn_storage(p) - deadstemn_storage(p)
!                     matrix_pot_deadstemn_xfer(p)         = matrix_cap_deadstemn_storage(p) - deadstemn_xfer(p)
!                     matrix_pot_livecrootn(p)             = matrix_cap_livecrootn(p) - livecrootn(p)
!                     matrix_pot_livecrootn_storage(p)     = matrix_cap_livecrootn_storage(p) - livecrootn_storage(p)
!                     matrix_pot_livecrootn_xfer(p)        = matrix_cap_livecrootn_storage(p) - livecrootn_xfer(p)
!                     matrix_pot_deadcrootn(p)             = matrix_cap_deadcrootn(p) - deadcrootn(p)
!                     matrix_pot_deadcrootn_storage(p)     = matrix_cap_deadcrootn_storage(p) - deadcrootn_storage(p)
!                     matrix_pot_deadcrootn_xfer(p)        = matrix_cap_deadcrootn_storage(p) - deadcrootn_xfer(p)
!                     if(ivt(p) >= npcropmin)then
!                        matrix_pot_grainn(p)                 = matrix_cap_grainn(p) - grainn(p)
!                        matrix_pot_grainn_storage(p)         = matrix_cap_grainn_storage(p) - grainn_storage(p)
!                        matrix_pot_grainn_xfer(p)            = matrix_cap_grainn_storage(p) - grainn_xfer(p)
!                     end if
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
                  matrix_ctransfer_fire_livestem_to_deadstem_acc(p) = 0._r8
                  matrix_ctransfer_fire_livecroot_to_deadcroot_acc(p) = 0._r8
   
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
   !               matrix_cturnover_gm_leaf_acc(p)                = 0._r8
   !               matrix_cturnover_gm_leafst_acc(p)                = 0._r8
   !               matrix_cturnover_gm_leafxf_acc(p)                = 0._r8
   !               matrix_cturnover_gm_froot_acc(p)                = 0._r8
   !               matrix_cturnover_gm_frootst_acc(p)                = 0._r8
   !               matrix_cturnover_gm_frootxf_acc(p)                = 0._r8
   !               matrix_cturnover_gm_livestem_acc(p)                = 0._r8
   !               matrix_cturnover_gm_livestemst_acc(p)                = 0._r8
   !               matrix_cturnover_gm_livestemxf_acc(p)                = 0._r8
   !               matrix_cturnover_gm_deadstem_acc(p)                = 0._r8
   !               matrix_cturnover_gm_deadstemst_acc(p)                = 0._r8
   !               matrix_cturnover_gm_deadstemxf_acc(p)                = 0._r8
   !               matrix_cturnover_gm_livecroot_acc(p)                = 0._r8
   !               matrix_cturnover_gm_livecrootst_acc(p)                = 0._r8
   !               matrix_cturnover_gm_livecrootxf_acc(p)                = 0._r8
   !               matrix_cturnover_gm_deadcroot_acc(p)                = 0._r8
   !               matrix_cturnover_gm_deadcrootst_acc(p)                = 0._r8
   !               matrix_cturnover_gm_deadcrootxf_acc(p)                = 0._r8
   !               matrix_cturnover_fire_leaf_acc(p)                = 0._r8
   !               matrix_cturnover_fire_leafst_acc(p)                = 0._r8
   !               matrix_cturnover_fire_leafxf_acc(p)                = 0._r8
   !               matrix_cturnover_fire_froot_acc(p)                = 0._r8
   !               matrix_cturnover_fire_frootst_acc(p)                = 0._r8
   !               matrix_cturnover_fire_frootxf_acc(p)                = 0._r8
   !               matrix_cturnover_fire_livestem_acc(p)                = 0._r8
   !               matrix_cturnover_fire_livestemst_acc(p)                = 0._r8
   !               matrix_cturnover_fire_livestemxf_acc(p)                = 0._r8
   !               matrix_cturnover_fire_deadstem_acc(p)                = 0._r8
   !               matrix_cturnover_fire_deadstemst_acc(p)                = 0._r8
   !               matrix_cturnover_fire_deadstemxf_acc(p)                = 0._r8
   !               matrix_cturnover_fire_livecroot_acc(p)                = 0._r8
   !               matrix_cturnover_fire_livecrootst_acc(p)                = 0._r8
   !               matrix_cturnover_fire_livecrootxf_acc(p)                = 0._r8
   !               matrix_cturnover_fire_deadcroot_acc(p)                = 0._r8
   !               matrix_cturnover_fire_deadcrootst_acc(p)                = 0._r8
   !               matrix_cturnover_fire_deadcrootxf_acc(p)                = 0._r8
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
                  call t_stopf('CN veg matrix-finalize spinup')
               end do 
            end if
         end if
   
         !print*,'here12'
!        matrix_Cinput(p) = 0._r8
!        matrix_C13input(p) = 0._r8
!        matrix_C14input(p) = 0._r8       
!         matrix_phtransfer(:,p) = 0._r8
!         matrix_gmtransfer(:,p) = 0._r8
!         matrix_fitransfer(p,:,:) = 0._r8
!         matrix_nphtransfer(p,:,:) = 0._r8
!         matrix_ngmtransfer(p,:,:) = 0._r8
!         matrix_nfitransfer(p,:,:) = 0._r8
!         matrix_Ninput(p) = 0._r8
      call vegmatrixc_input%ReleaseV()
      call vegmatrixc13_input%ReleaseV()
      call vegmatrixc14_input%ReleaseV()
      call vegmatrixn_input%ReleaseV()
         
    
   end associate 
 end subroutine CNVegMatrix
 subroutine inverse(a,c,n)
!============================================================
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
