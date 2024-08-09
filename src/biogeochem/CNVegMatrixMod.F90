module CNVegMatrixMod

  !---------------------------------------------------------------------------------------
  ! The matrix model of CLM5.0 was developed by Yiqi Luo EcoLab members, 
  ! Drs. Xingjie Lu, Yuanyuan Huang and Zhengguang Du, at Northern Arizona University
  !---------------------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  ! Matrix solution for vegetation C and N cycles
  ! The matrix equation 
  ! Xn+1 = Xn + B*I*dt + (Aph*Kph + Agm*Kgm + Afi*Kfi) * Xn*dt
  ! Xn is the state variable of last time step n, and Xn+1 is the state variable of 
  ! the next time step n+1, I is the input to the vegetation, i.e. NPP in this case.
  ! B is allocation fraction vector.
  ! Aph, Agm and Afi represent transfer coefficient matrix A from phenology, gap mortality 
  ! and fire related C and N transfers.
  ! Kph, Kgm and Kfi represent turnover rate matrix K from phenology, gap mortality
  ! and fire related C and N transfers.
  !---------------------------------------------------------------------------------------
  
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use clm_time_manager               , only : get_step_size,is_end_curr_year,is_first_step_of_this_run_segment,&
                                              is_beg_curr_year,update_DA_nstep
  use decompMod                      , only : bounds_type 
  use clm_varcon                     , only : spval
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
  use pftconMod                      , only : pftcon,npcropmin
  use CNVegCarbonStateType           , only : cnveg_carbonstate_type
  use CNVegNitrogenStateType         , only : cnveg_nitrogenstate_type
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type     !include: callocation,ctransfer, cturnover
  use CNVegNitrogenFluxType          , only : cnveg_nitrogenflux_type
  use CNVegStateType                 , only : cnveg_state_type
  use SoilBiogeochemNitrogenFluxType , only : soilbiogeochem_nitrogenflux_type
  use clm_varctl                     , only : spinup_matrixcn, hist_wrt_matrixcn_diag, nyr_forcing, nyr_SASU, iloop_avg 
  use clm_varctl                     , only : use_c13, use_c14 
  use SparseMatrixMultiplyMod        , only : sparse_matrix_type,diag_matrix_type,vector_type
  use MatrixMod                      , only : inverse
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNVegMatrix
  public:: matrix_update_phc,matrix_update_gmc,matrix_update_fic
  public:: matrix_update_phn,matrix_update_gmn,matrix_update_fin
  public:: CNVegMatrixRest

  ! ! PRIVATE MEMBER DATA:
  integer,save, private :: iyr=0   ! Cycling year number into forcing sequence
  integer,save, private :: iloop=0 ! The iloop^th forcing loop
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

!    LOCAL VARIABLES:
     integer :: fc,fp,j,i,k    ! indices
     integer :: p,c         !  

  ! Temporary variables matrix A for different processes
     real(r8),dimension(:,:)   :: Aphconed(bounds%begp:bounds%endp,ncphtrans-ncphouttrans)
     real(r8),dimension(:,:)   :: Aphnoned(bounds%begp:bounds%endp,nnphtrans-nnphouttrans)
     real(r8),dimension(:,:)   :: Agmconed(bounds%begp:bounds%endp,ncgmtrans-ncgmouttrans)
     real(r8),dimension(:,:)   :: Agmnoned(bounds%begp:bounds%endp,nngmtrans-nngmouttrans)
     real(r8),dimension(:,:)   :: Aficoned(bounds%begp:bounds%endp,ncfitrans-ncfiouttrans)
     real(r8),dimension(:,:)   :: Afic14oned(bounds%begp:bounds%endp,ncfitrans-ncfiouttrans)
     real(r8),dimension(:,:)   :: Afinoned(bounds%begp:bounds%endp,nnfitrans-nnfiouttrans)

  ! Temporary variables saving row indices of all transfers in different processes
     integer,dimension(:)      :: AI_phc(ncphtrans-ncphouttrans)
     integer,dimension(:)      :: AI_phn(nnphtrans-nnphouttrans)
     integer,dimension(:)      :: AI_gmc(ncgmtrans-ncgmouttrans)
     integer,dimension(:)      :: AI_gmn(nngmtrans-nngmouttrans)
     integer,dimension(:)      :: AI_fic(ncfitrans-ncfiouttrans)
     integer,dimension(:)      :: AI_fic14(ncfitrans-ncfiouttrans)
     integer,dimension(:)      :: AI_fin(nnfitrans-nnfiouttrans)

  ! Temporary variables saving column indices of all transfers in different processes
     integer,dimension(:)      :: AJ_phc(ncphtrans-ncphouttrans)
     integer,dimension(:)      :: AJ_phn(nnphtrans-nnphouttrans)
     integer,dimension(:)      :: AJ_gmc(ncgmtrans-ncgmouttrans)
     integer,dimension(:)      :: AJ_gmn(nngmtrans-nngmouttrans)
     integer,dimension(:)      :: AJ_fic(ncfitrans-ncfiouttrans)
     integer,dimension(:)      :: AJ_fic14(ncfitrans-ncfiouttrans)
     integer,dimension(:)      :: AJ_fin(nnfitrans-nnfiouttrans)

  ! Temporary variables for matrix operation, which save C and N inputs to different vegetation compartments as a vector type.
     type(vector_type)         :: vegmatrixc_input    
     type(vector_type)         :: vegmatrixc13_input  
     type(vector_type)         :: vegmatrixc14_input  
     type(vector_type)         :: vegmatrixn_input    

  ! "init" indicators indicate whether A matrices have been initialized. 
     logical, save             :: init_ready_aphc      = .false.
     logical, save             :: init_ready_agmc      = .false.
     logical, save             :: init_ready_afic      = .false.
     logical, save             :: init_ready_afic14    = .false.
     logical, save             :: init_ready_aphn      = .false.
     logical, save             :: init_ready_agmn      = .false.
     logical, save             :: init_ready_afin      = .false.

  ! "list" indicators indicate whether operation of sparse matrix plus SPMP_AB or SPMP_ABC has already been saved.
     logical, save             :: list_ready_phgmfic   = .false.
     logical, save             :: list_ready_phgmfic14 = .false.
     logical, save             :: list_ready_phgmc     = .false.
     logical, save             :: list_ready_phgmfin   = .false.
     logical, save             :: list_ready_phgmn     = .false.

  ! Temporary variables are only used at end of the year to calculate C and N storage capacity
     real(r8),dimension(:)     :: matrix_calloc_acc      (1:nvegcpool)
     real(r8),dimension(:)     :: matrix_nalloc_acc      (1:nvegnpool)
     real(r8),dimension(:,:)   :: matrix_ctransfer_acc   (1:nvegcpool,1:nvegcpool)
     real(r8),dimension(:,:)   :: matrix_ntransfer_acc   (1:nvegnpool,1:nvegnpool)
     real(r8),dimension(:)     :: matrix_c13alloc_acc    (1:nvegcpool)
     real(r8),dimension(:,:)   :: matrix_c13transfer_acc (1:nvegcpool,1:nvegcpool)
     real(r8),dimension(:)     :: matrix_c14alloc_acc    (1:nvegcpool)
     real(r8),dimension(:,:)   :: matrix_c14transfer_acc (1:nvegcpool,1:nvegcpool)

  ! Local variables for capacity calculation and spin up
     real(r8),dimension(:)     :: vegmatrixc_rt(1:nvegcpool) ! C storage capacity
     real(r8),dimension(:)     :: vegmatrixc13_rt(1:nvegcpool) ! C13 storage capacity
     real(r8),dimension(:)     :: vegmatrixc14_rt(1:nvegcpool) ! C14 storage capacity
     real(r8),dimension(:)     :: vegmatrixn_rt(1:nvegnpool) ! N storage capacity
     real(r8),dimension(:,:)   :: AKinvc(1:nvegcpool,1:nvegcpool),AKinvn(1:nvegnpool,1:nvegnpool)
     real(r8):: epsi 

     
     real(r8):: dt        ! time step (seconds)
#ifdef _OPENMP
     integer, external :: OMP_GET_MAX_THREADS
     integer :: nthreads  ! Number of threads
#else
     integer, parameter :: nthreads = 0 ! Number of threads
#endif
     integer, parameter :: irepr = 1    ! Reproductive index to use for grain
 
fr:  associate(                          &
     ivt                   => patch%itype                                       , & ! Input:  [integer  (:) ]  patch vegetation type
     cf13_veg              => c13_cnveg_carbonflux_inst                         , & ! In
     cf14_veg              => c14_cnveg_carbonflux_inst                         , & ! In
     cs13_veg              => c13_cnveg_carbonstate_inst                        , & ! In/Output
     cs14_veg              => c14_cnveg_carbonstate_inst                        , & ! In/Output  

    fire_closs            => cnveg_carbonflux_inst%fire_closs_patch            , &

  ! Original vegetation variables are updated by matrix operation in this module
    leafc                 => cnveg_carbonstate_inst%leafc_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C 
    leafc_storage         => cnveg_carbonstate_inst%leafc_storage_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf storage C 
    leafc_xfer            => cnveg_carbonstate_inst%leafc_xfer_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf transfer C 
    frootc                => cnveg_carbonstate_inst%frootc_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root C 
    frootc_storage        => cnveg_carbonstate_inst%frootc_storage_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root storage C 
    frootc_xfer           => cnveg_carbonstate_inst%frootc_xfer_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) fine root transfer C 
    livestemc             => cnveg_carbonstate_inst%livestemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem C 
    livestemc_storage     => cnveg_carbonstate_inst%livestemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem storage C 
    livestemc_xfer        => cnveg_carbonstate_inst%livestemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live stem transfer C
    deadstemc             => cnveg_carbonstate_inst%deadstemc_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem C
    deadstemc_storage     => cnveg_carbonstate_inst%deadstemc_storage_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem storage C
    deadstemc_xfer        => cnveg_carbonstate_inst%deadstemc_xfer_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead stem transfer C
    livecrootc            => cnveg_carbonstate_inst%livecrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root C
    livecrootc_storage    => cnveg_carbonstate_inst%livecrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root storage C
    livecrootc_xfer       => cnveg_carbonstate_inst%livecrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) live coarse root transfer C
    deadcrootc            => cnveg_carbonstate_inst%deadcrootc_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root C
    deadcrootc_storage    => cnveg_carbonstate_inst%deadcrootc_storage_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root storage C
    deadcrootc_xfer       => cnveg_carbonstate_inst%deadcrootc_xfer_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) dead coarse root transfer C
    reproductivec         => cnveg_carbonstate_inst%reproductivec_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) grain C
    reproductivec_storage => cnveg_carbonstate_inst%reproductivec_storage_patch , & ! In/Output:  [real(r8) (:) ]    (gC/m2) grain storage C
    reproductivec_xfer    => cnveg_carbonstate_inst%reproductivec_xfer_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) grain transfer C

    leafn                 => cnveg_nitrogenstate_inst%leafn_patch               , & ! In/Output:  [real(r8) (:) ]    (gN/m2) leaf N
    leafn_storage         => cnveg_nitrogenstate_inst%leafn_storage_patch       , & ! In/Output:  [real(r8) (:) ]    (gN/m2) leaf storage N
    leafn_xfer            => cnveg_nitrogenstate_inst%leafn_xfer_patch          , & ! In/Output:  [real(r8) (:) ]    (gN/m2) leaf transfer N
    frootn                => cnveg_nitrogenstate_inst%frootn_patch              , & ! In/Output:  [real(r8) (:) ]    (gN/m2) fine root N
    frootn_storage        => cnveg_nitrogenstate_inst%frootn_storage_patch      , & ! In/Output:  [real(r8) (:) ]    (gN/m2) fine root storage N
    frootn_xfer           => cnveg_nitrogenstate_inst%frootn_xfer_patch         , & ! In/Output:  [real(r8) (:) ]    (gN/m2) fine root transfer N
    livestemn             => cnveg_nitrogenstate_inst%livestemn_patch           , & ! In/Output:  [real(r8) (:) ]    (gN/m2) live stem N
    livestemn_storage     => cnveg_nitrogenstate_inst%livestemn_storage_patch   , & ! In/Output:  [real(r8) (:) ]    (gN/m2) live stem storage N
    livestemn_xfer        => cnveg_nitrogenstate_inst%livestemn_xfer_patch      , & ! In/Output:  [real(r8) (:) ]    (gN/m2) live stem transfer N
    deadstemn             => cnveg_nitrogenstate_inst%deadstemn_patch           , & ! In/Output:  [real(r8) (:) ]    (gN/m2) dead stem N
    deadstemn_storage     => cnveg_nitrogenstate_inst%deadstemn_storage_patch   , & ! In/Output:  [real(r8) (:) ]    (gN/m2) dead stem storage N
    deadstemn_xfer        => cnveg_nitrogenstate_inst%deadstemn_xfer_patch      , & ! In/Output:  [real(r8) (:) ]    (gN/m2) dead stem transfer N
    livecrootn            => cnveg_nitrogenstate_inst%livecrootn_patch          , & ! In/Output:  [real(r8) (:) ]    (gN/m2) live coarse root N
    livecrootn_storage    => cnveg_nitrogenstate_inst%livecrootn_storage_patch  , & ! In/Output:  [real(r8) (:) ]    (gN/m2) live coarse root storage N
    livecrootn_xfer       => cnveg_nitrogenstate_inst%livecrootn_xfer_patch     , & ! In/Output:  [real(r8) (:) ]    (gN/m2) live coarse root transfer N
    deadcrootn            => cnveg_nitrogenstate_inst%deadcrootn_patch          , & ! In/Output:  [real(r8) (:) ]    (gN/m2) dead coarse root N
    deadcrootn_storage    => cnveg_nitrogenstate_inst%deadcrootn_storage_patch  , & ! In/Output:  [real(r8) (:) ]    (gN/m2) dead coarse root storage N
    deadcrootn_xfer       => cnveg_nitrogenstate_inst%deadcrootn_xfer_patch     , & ! In/Output:  [real(r8) (:) ]    (gN/m2) dead coarse root transfer N
    reproductiven         => cnveg_nitrogenstate_inst%reproductiven_patch              , & ! In/Output:  [real(r8) (:) ]    (gN/m2) grain N
    reproductiven_storage => cnveg_nitrogenstate_inst%reproductiven_storage_patch      , & ! In/Output:  [real(r8) (:) ]    (gN/m2) grain storage N
    reproductiven_xfer    => cnveg_nitrogenstate_inst%reproductiven_xfer_patch         , & ! In/Output:  [real(r8) (:) ]    (gN/m2) grain transfer N
    retransn              => cnveg_nitrogenstate_inst%retransn_patch            , & ! In/Output:  [real(r8) (:) ]    (gN/m2) plant retranslocated N

    leafc_SASUsave                => cnveg_carbonstate_inst%leafc_SASUsave_patch                 , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for SASU
    leafc_storage_SASUsave        => cnveg_carbonstate_inst%leafc_storage_SASUsave_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for SASU
    leafc_xfer_SASUsave           => cnveg_carbonstate_inst%leafc_xfer_SASUsave_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf C for SASU
    frootc_SASUsave               => cnveg_carbonstate_inst%frootc_SASUsave_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) froot C for SASU
    frootc_storage_SASUsave       => cnveg_carbonstate_inst%frootc_storage_SASUsave_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) froot C for SASU
    frootc_xfer_SASUsave          => cnveg_carbonstate_inst%frootc_xfer_SASUsave_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) froot C for SASU
    livestemc_SASUsave            => cnveg_carbonstate_inst%livestemc_SASUsave_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livestem C for SASU
    livestemc_storage_SASUsave    => cnveg_carbonstate_inst%livestemc_storage_SASUsave_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livestem C for SASU
    livestemc_xfer_SASUsave       => cnveg_carbonstate_inst%livestemc_xfer_SASUsave_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livestem C for SASU
    deadstemc_SASUsave            => cnveg_carbonstate_inst%deadstemc_SASUsave_patch             , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadstem C for SASU
    deadstemc_storage_SASUsave    => cnveg_carbonstate_inst%deadstemc_storage_SASUsave_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadstem C for SASU
    deadstemc_xfer_SASUsave       => cnveg_carbonstate_inst%deadstemc_xfer_SASUsave_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadstem C for SASU
    livecrootc_SASUsave           => cnveg_carbonstate_inst%livecrootc_SASUsave_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livecroot C for SASU
    livecrootc_storage_SASUsave   => cnveg_carbonstate_inst%livecrootc_storage_SASUsave_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livecroot C for SASU
    livecrootc_xfer_SASUsave      => cnveg_carbonstate_inst%livecrootc_xfer_SASUsave_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livecroot C for SASU
    deadcrootc_SASUsave           => cnveg_carbonstate_inst%deadcrootc_SASUsave_patch            , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadcroot C for SASU
    deadcrootc_storage_SASUsave   => cnveg_carbonstate_inst%deadcrootc_storage_SASUsave_patch    , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadcroot C for SASU
    deadcrootc_xfer_SASUsave      => cnveg_carbonstate_inst%deadcrootc_xfer_SASUsave_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadcroot C for SASU
    grainc_SASUsave               => cnveg_carbonstate_inst%grainc_SASUsave_patch                , & ! In/Output:  [real(r8) (:) ]    (gC/m2) grain C for SASU
    grainc_storage_SASUsave       => cnveg_carbonstate_inst%grainc_storage_SASUsave_patch        , & ! In/Output:  [real(r8) (:) ]    (gC/m2) grain storage C for SASU

    leafn_SASUsave                => cnveg_nitrogenstate_inst%leafn_SASUsave_patch               , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf N for SASU
    leafn_storage_SASUsave        => cnveg_nitrogenstate_inst%leafn_storage_SASUsave_patch       , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf N for SASU
    leafn_xfer_SASUsave           => cnveg_nitrogenstate_inst%leafn_xfer_SASUsave_patch          , & ! In/Output:  [real(r8) (:) ]    (gC/m2) leaf N for SASU
    frootn_SASUsave               => cnveg_nitrogenstate_inst%frootn_SASUsave_patch              , & ! In/Output:  [real(r8) (:) ]    (gC/m2) froot N for SASU
    frootn_storage_SASUsave       => cnveg_nitrogenstate_inst%frootn_storage_SASUsave_patch      , & ! In/Output:  [real(r8) (:) ]    (gC/m2) froot N for SASU
    frootn_xfer_SASUsave          => cnveg_nitrogenstate_inst%frootn_xfer_SASUsave_patch         , & ! In/Output:  [real(r8) (:) ]    (gC/m2) froot N for SASU
    livestemn_SASUsave            => cnveg_nitrogenstate_inst%livestemn_SASUsave_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livestem N for SASU
    livestemn_storage_SASUsave    => cnveg_nitrogenstate_inst%livestemn_storage_SASUsave_patch   , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livestem N for SASU
    livestemn_xfer_SASUsave       => cnveg_nitrogenstate_inst%livestemn_xfer_SASUsave_patch      , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livestem N for SASU
    deadstemn_SASUsave            => cnveg_nitrogenstate_inst%deadstemn_SASUsave_patch           , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadstem N for SASU
    deadstemn_storage_SASUsave    => cnveg_nitrogenstate_inst%deadstemn_storage_SASUsave_patch   , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadstem N for SASU
    deadstemn_xfer_SASUsave       => cnveg_nitrogenstate_inst%deadstemn_xfer_SASUsave_patch      , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadstem N for SASU
    livecrootn_SASUsave           => cnveg_nitrogenstate_inst%livecrootn_SASUsave_patch          , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livecroot N for SASU
    livecrootn_storage_SASUsave   => cnveg_nitrogenstate_inst%livecrootn_storage_SASUsave_patch  , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livecroot N for SASU
    livecrootn_xfer_SASUsave      => cnveg_nitrogenstate_inst%livecrootn_xfer_SASUsave_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) livecroot N for SASU
    deadcrootn_SASUsave           => cnveg_nitrogenstate_inst%deadcrootn_SASUsave_patch          , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadcroot N for SASU
    deadcrootn_storage_SASUsave   => cnveg_nitrogenstate_inst%deadcrootn_storage_SASUsave_patch  , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadcroot N for SASU
    deadcrootn_xfer_SASUsave      => cnveg_nitrogenstate_inst%deadcrootn_xfer_SASUsave_patch     , & ! In/Output:  [real(r8) (:) ]    (gC/m2) deadcroot N for SASU
    grainn_SASUsave               => cnveg_nitrogenstate_inst%grainn_SASUsave_patch              , & ! In/Output:  [real(r8) (:) ]    (gC/m2) grain N for SASU
    grainn_storage_SASUsave       => cnveg_nitrogenstate_inst%grainn_storage_SASUsave_patch      , & ! In/Output:  [real(r8) (:) ]    (gC/m2) grain storage N for SASU

  ! Vegetation capacity variables "matrix_cap_*", save the capacity of each vegetation compartment.
    matrix_cap_leafc              => cnveg_carbonstate_inst%matrix_cap_leafc_patch               ,&!Output:[real(r8)(:)] (gC/m2) leaf C capacity
    matrix_cap_leafc_storage      => cnveg_carbonstate_inst%matrix_cap_leafc_storage_patch       ,&!Output:[real(r8)(:)] (gC/m2) leaf storage C capacity
    matrix_cap_leafc_xfer         => cnveg_carbonstate_inst%matrix_cap_leafc_xfer_patch          ,&!Output:[real(r8)(:)] (gC/m2) leaf transfer C capacity
    matrix_cap_frootc             => cnveg_carbonstate_inst%matrix_cap_frootc_patch              ,&!Output:[real(r8)(:)] (gC/m2) fine root C capacity
    matrix_cap_frootc_storage     => cnveg_carbonstate_inst%matrix_cap_frootc_storage_patch      ,&!Output:[real(r8)(:)] (gC/m2) fine root storage C capacity
    matrix_cap_frootc_xfer        => cnveg_carbonstate_inst%matrix_cap_frootc_xfer_patch         ,&!Output:[real(r8)(:)] (gC/m2) fine root transfer C capacity
    matrix_cap_livestemc          => cnveg_carbonstate_inst%matrix_cap_livestemc_patch           ,&!Output:[real(r8)(:)] (gC/m2) live stem C capacity
    matrix_cap_livestemc_storage  => cnveg_carbonstate_inst%matrix_cap_livestemc_storage_patch   ,&!Output:[real(r8)(:)] (gC/m2) live stem storage C capacity
    matrix_cap_livestemc_xfer     => cnveg_carbonstate_inst%matrix_cap_livestemc_xfer_patch      ,&!Output:[real(r8)(:)] (gC/m2) live stem transfer C capacity
    matrix_cap_deadstemc          => cnveg_carbonstate_inst%matrix_cap_deadstemc_patch           ,&!Output:[real(r8)(:)] (gC/m2) dead stem C capacity
    matrix_cap_deadstemc_storage  => cnveg_carbonstate_inst%matrix_cap_deadstemc_storage_patch   ,&!Output:[real(r8)(:)] (gC/m2) dead stem storage C capaicty
    matrix_cap_deadstemc_xfer     => cnveg_carbonstate_inst%matrix_cap_deadstemc_xfer_patch      ,&!Output:[real(r8)(:)] (gC/m2) dead stem transfer C capacity
    matrix_cap_livecrootc         => cnveg_carbonstate_inst%matrix_cap_livecrootc_patch          ,&!Output:[real(r8)(:)] (gC/m2) live coarse root C capacity
    matrix_cap_livecrootc_storage => cnveg_carbonstate_inst%matrix_cap_livecrootc_storage_patch  ,&!Output:[real(r8)(:)] (gC/m2) live coarse root storage C capacity
    matrix_cap_livecrootc_xfer    => cnveg_carbonstate_inst%matrix_cap_livecrootc_xfer_patch     ,&!Output:[real(r8)(:)] (gC/m2) live coarse root transfer C capacity
    matrix_cap_deadcrootc         => cnveg_carbonstate_inst%matrix_cap_deadcrootc_patch          ,&!Output:[real(r8)(:)] (gC/m2) dead coarse root C capacity
    matrix_cap_deadcrootc_storage => cnveg_carbonstate_inst%matrix_cap_deadcrootc_storage_patch  ,&!Output:[real(r8)(:)] (gC/m2) dead coarse root storage C capacity
    matrix_cap_deadcrootc_xfer    => cnveg_carbonstate_inst%matrix_cap_deadcrootc_xfer_patch     ,&!Output:[real(r8)(:)] (gC/m2) dead coarse root transfer C capacity
    matrix_cap_reproc             => cnveg_carbonstate_inst%matrix_cap_reproc_patch              ,&!Output:[real(r8)(:)] (gC/m2) grain C capacity
    matrix_cap_reproc_storage     => cnveg_carbonstate_inst%matrix_cap_reproc_storage_patch      ,&!Output:[real(r8)(:)] (gC/m2) grain storage C capacity
    matrix_cap_reproc_xfer        => cnveg_carbonstate_inst%matrix_cap_reproc_xfer_patch         ,&!Output:[real(r8)(:)] (gC/m2) grain transfer C

    matrix_cap_leafn              => cnveg_nitrogenstate_inst%matrix_cap_leafn_patch             ,&!Output:[real(r8)(:)] (gN/m2) leaf N capacity
    matrix_cap_leafn_storage      => cnveg_nitrogenstate_inst%matrix_cap_leafn_storage_patch     ,&!Output:[real(r8)(:)] (gN/m2) leaf storage N capacity
    matrix_cap_leafn_xfer         => cnveg_nitrogenstate_inst%matrix_cap_leafn_xfer_patch        ,&!Output:[real(r8)(:)] (gN/m2) leaf transfer N capacity 
    matrix_cap_frootn             => cnveg_nitrogenstate_inst%matrix_cap_frootn_patch            ,&!Output:[real(r8)(:)] (gN/m2) fine root N capacity
    matrix_cap_frootn_storage     => cnveg_nitrogenstate_inst%matrix_cap_frootn_storage_patch    ,&!Output:[real(r8)(:)] (gN/m2) fine root storage N capacity
    matrix_cap_frootn_xfer        => cnveg_nitrogenstate_inst%matrix_cap_frootn_xfer_patch       ,&!Output:[real(r8)(:)] (gN/m2) fine root transfer N capacity
    matrix_cap_livestemn          => cnveg_nitrogenstate_inst%matrix_cap_livestemn_patch         ,&!Output:[real(r8)(:)] (gN/m2) live stem N capacity
    matrix_cap_livestemn_storage  => cnveg_nitrogenstate_inst%matrix_cap_livestemn_storage_patch ,&!Output:[real(r8)(:)] (gN/m2) live stem storage N capacity
    matrix_cap_livestemn_xfer     => cnveg_nitrogenstate_inst%matrix_cap_livestemn_xfer_patch    ,&!Output:[real(r8)(:)] (gN/m2) live stem transfer N capacity
    matrix_cap_deadstemn          => cnveg_nitrogenstate_inst%matrix_cap_deadstemn_patch         ,&!Output:[real(r8)(:)] (gN/m2) dead stem N capacity
    matrix_cap_deadstemn_storage  => cnveg_nitrogenstate_inst%matrix_cap_deadstemn_storage_patch ,&!Output:[real(r8)(:)] (gN/m2) dead stem storage N capacity
    matrix_cap_deadstemn_xfer     => cnveg_nitrogenstate_inst%matrix_cap_deadstemn_xfer_patch    ,&!Output:[real(r8)(:)] (gN/m2) dead stem transfer N capacity
    matrix_cap_livecrootn         => cnveg_nitrogenstate_inst%matrix_cap_livecrootn_patch        ,&!Output:[real(r8)(:)] (gN/m2) live coarse root N capacity
    matrix_cap_livecrootn_storage => cnveg_nitrogenstate_inst%matrix_cap_livecrootn_storage_patch,&!Output:[real(r8)(:)] (gN/m2) live coarse root storage N capacity
    matrix_cap_livecrootn_xfer    => cnveg_nitrogenstate_inst%matrix_cap_livecrootn_xfer_patch   ,&!Output:[real(r8)(:)] (gN/m2) live coarse root transfer N capacity
    matrix_cap_deadcrootn         => cnveg_nitrogenstate_inst%matrix_cap_deadcrootn_patch        ,&!Output:[real(r8)(:)] (gN/m2) dead coarse root N capacity
    matrix_cap_deadcrootn_storage => cnveg_nitrogenstate_inst%matrix_cap_deadcrootn_storage_patch,&!Output:[real(r8)(:)] (gN/m2) dead coarse root storage N capacity
    matrix_cap_deadcrootn_xfer    => cnveg_nitrogenstate_inst%matrix_cap_deadcrootn_xfer_patch   ,&!Output:[real(r8)(:)] (gN/m2) dead coarse root transfer N capacity
    matrix_cap_repron             => cnveg_nitrogenstate_inst%matrix_cap_repron_patch            ,&!Output:[real(r8)(:)] (gN/m2) grain N capacity
    matrix_cap_repron_storage     => cnveg_nitrogenstate_inst%matrix_cap_repron_storage_patch    ,&!Output:[real(r8)(:)] (gN/m2) grain storage N capacity
    matrix_cap_repron_xfer        => cnveg_nitrogenstate_inst%matrix_cap_repron_xfer_patch       ,&!Output:[real(r8)(:)] (gN/m2) grain transfer N capacity

  ! Variables matrix_calloc_*_acc, matrix_ctransfer_*_acc, and matrix_cturnover_*_acc are used to calculate the C capacity as the C steady state estimates in spin up.
  ! These variables are all state variables, saving accumulated N transfers during the calendar year.
    matrix_calloc_leaf_acc        => cnveg_carbonstate_inst%matrix_calloc_leaf_acc_patch        , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to leaf during this year
    matrix_calloc_leafst_acc      => cnveg_carbonstate_inst%matrix_calloc_leafst_acc_patch      , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to leaf storage during this year
    matrix_calloc_froot_acc       => cnveg_carbonstate_inst%matrix_calloc_froot_acc_patch       , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to fine root during this year
    matrix_calloc_frootst_acc     => cnveg_carbonstate_inst%matrix_calloc_frootst_acc_patch     , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to fine root storage during this year
    matrix_calloc_livestem_acc    => cnveg_carbonstate_inst%matrix_calloc_livestem_acc_patch    , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to live stem during this year
    matrix_calloc_livestemst_acc  => cnveg_carbonstate_inst%matrix_calloc_livestemst_acc_patch  , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to live stem storage during this year
    matrix_calloc_deadstem_acc    => cnveg_carbonstate_inst%matrix_calloc_deadstem_acc_patch    , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to dead stem during this year
    matrix_calloc_deadstemst_acc  => cnveg_carbonstate_inst%matrix_calloc_deadstemst_acc_patch  , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to dead stem storage during this year
    matrix_calloc_livecroot_acc   => cnveg_carbonstate_inst%matrix_calloc_livecroot_acc_patch   , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to live corase root during this year
    matrix_calloc_livecrootst_acc => cnveg_carbonstate_inst%matrix_calloc_livecrootst_acc_patch , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to live corase root storage during this year
    matrix_calloc_deadcroot_acc   => cnveg_carbonstate_inst%matrix_calloc_deadcroot_acc_patch   , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to dead corase root during this year
    matrix_calloc_deadcrootst_acc => cnveg_carbonstate_inst%matrix_calloc_deadcrootst_acc_patch , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to dead corase root storage during this year
    matrix_calloc_grain_acc       => cnveg_carbonstate_inst%matrix_calloc_grain_acc_patch       , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to grain during this year
    matrix_calloc_grainst_acc     => cnveg_carbonstate_inst%matrix_calloc_grainst_acc_patch     , & 
                                     ! In/Output: [real(r8) (:) ] (gC/m2/year) Input C allocated to grain storage during this year
 
    matrix_ctransfer_leafst_to_leafxf_acc            => cnveg_carbonstate_inst%matrix_ctransfer_leafst_to_leafxf_acc_patch            , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from leaf storage to leaf transfer pool during this year
    matrix_ctransfer_leafxf_to_leaf_acc              => cnveg_carbonstate_inst%matrix_ctransfer_leafxf_to_leaf_acc_patch              , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from leaf transfer to leaf pool during this year
    matrix_ctransfer_frootst_to_frootxf_acc          => cnveg_carbonstate_inst%matrix_ctransfer_frootst_to_frootxf_acc_patch          , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from fine root storage to fine root transfer pool during this year
    matrix_ctransfer_frootxf_to_froot_acc            => cnveg_carbonstate_inst%matrix_ctransfer_frootxf_to_froot_acc_patch            , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from fine root transfer to fine root pool during this year
    matrix_ctransfer_livestemst_to_livestemxf_acc    => cnveg_carbonstate_inst%matrix_ctransfer_livestemst_to_livestemxf_acc_patch    , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from live stem storage to live stem transfer pool during this year
    matrix_ctransfer_livestemxf_to_livestem_acc      => cnveg_carbonstate_inst%matrix_ctransfer_livestemxf_to_livestem_acc_patch      , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from live stem transfer to live stem pool during this year
    matrix_ctransfer_deadstemst_to_deadstemxf_acc    => cnveg_carbonstate_inst%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch    , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from dead stem storage to dead stem transfer pool during this year
    matrix_ctransfer_deadstemxf_to_deadstem_acc      => cnveg_carbonstate_inst%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch      , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from dead stem transfer to dead stem pool during this year
    matrix_ctransfer_livecrootst_to_livecrootxf_acc  => cnveg_carbonstate_inst%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch  , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from live coarse root storage to live coarse root transfer pool during this year
    matrix_ctransfer_livecrootxf_to_livecroot_acc    => cnveg_carbonstate_inst%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch    , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from live coarse root transfer to live coarse root pool during this year
    matrix_ctransfer_deadcrootst_to_deadcrootxf_acc  => cnveg_carbonstate_inst%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch  , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from dead coarse root storage to dead coarse root transfer pool during this year
    matrix_ctransfer_deadcrootxf_to_deadcroot_acc    => cnveg_carbonstate_inst%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch    , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from dead coarse root transfer to dead coarse root pool during this year
    matrix_ctransfer_grainst_to_grainxf_acc          => cnveg_carbonstate_inst%matrix_ctransfer_grainst_to_grainxf_acc_patch          , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from grain storage to grain transfer pool during this year
    matrix_ctransfer_grainxf_to_grain_acc            => cnveg_carbonstate_inst%matrix_ctransfer_grainxf_to_grain_acc_patch            , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from grain transfer to grain pool during this year
    matrix_ctransfer_livestem_to_deadstem_acc        => cnveg_carbonstate_inst%matrix_ctransfer_livestem_to_deadstem_acc_patch        , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from live stem to dead stem pool during this year
    matrix_ctransfer_livecroot_to_deadcroot_acc      => cnveg_carbonstate_inst%matrix_ctransfer_livecroot_to_deadcroot_acc_patch      , &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C transfer from live coarse root to dead coarse root pool during this year

    matrix_cturnover_leaf_acc        => cnveg_carbonstate_inst%matrix_cturnover_leaf_acc_patch        , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from leaf
    matrix_cturnover_leafst_acc      => cnveg_carbonstate_inst%matrix_cturnover_leafst_acc_patch      , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from leaf storage 
    matrix_cturnover_leafxf_acc      => cnveg_carbonstate_inst%matrix_cturnover_leafxf_acc_patch      , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from leaf transfer
    matrix_cturnover_froot_acc       => cnveg_carbonstate_inst%matrix_cturnover_froot_acc_patch       , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from fine root 
    matrix_cturnover_frootst_acc     => cnveg_carbonstate_inst%matrix_cturnover_frootst_acc_patch     , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from fine root storage
    matrix_cturnover_frootxf_acc     => cnveg_carbonstate_inst%matrix_cturnover_frootxf_acc_patch     , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from fine root transfer 
    matrix_cturnover_livestem_acc    => cnveg_carbonstate_inst%matrix_cturnover_livestem_acc_patch    , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from live stem
    matrix_cturnover_livestemst_acc  => cnveg_carbonstate_inst%matrix_cturnover_livestemst_acc_patch  , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from live stem storage
    matrix_cturnover_livestemxf_acc  => cnveg_carbonstate_inst%matrix_cturnover_livestemxf_acc_patch  , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from live stem transfer
    matrix_cturnover_deadstem_acc    => cnveg_carbonstate_inst%matrix_cturnover_deadstem_acc_patch    , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from dead stem 
    matrix_cturnover_deadstemst_acc  => cnveg_carbonstate_inst%matrix_cturnover_deadstemst_acc_patch  , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from dead stem storage
    matrix_cturnover_deadstemxf_acc  => cnveg_carbonstate_inst%matrix_cturnover_deadstemxf_acc_patch  , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from dead stem transfer
    matrix_cturnover_livecroot_acc   => cnveg_carbonstate_inst%matrix_cturnover_livecroot_acc_patch   , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from live coarse root
    matrix_cturnover_livecrootst_acc => cnveg_carbonstate_inst%matrix_cturnover_livecrootst_acc_patch , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from live coarse root storage
    matrix_cturnover_livecrootxf_acc => cnveg_carbonstate_inst%matrix_cturnover_livecrootxf_acc_patch , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from live coarse root transfer
    matrix_cturnover_deadcroot_acc   => cnveg_carbonstate_inst%matrix_cturnover_deadcroot_acc_patch   , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from dead coarse root 
    matrix_cturnover_deadcrootst_acc => cnveg_carbonstate_inst%matrix_cturnover_deadcrootst_acc_patch , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from dead coarse root storage
    matrix_cturnover_deadcrootxf_acc => cnveg_carbonstate_inst%matrix_cturnover_deadcrootxf_acc_patch , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from dead coarse root transfer
    matrix_cturnover_grain_acc       => cnveg_carbonstate_inst%matrix_cturnover_grain_acc_patch       , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from grain
    matrix_cturnover_grainst_acc     => cnveg_carbonstate_inst%matrix_cturnover_grainst_acc_patch     , & 
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from grain storage
    matrix_cturnover_grainxf_acc     => cnveg_carbonstate_inst%matrix_cturnover_grainxf_acc_patch       &
                    ! In/Output: [real(r8) (:) ] (gC/m2/year) C turnover from grain transfer
    )
od: associate(                          &

  ! Variables matrix_nalloc_*_acc, matrix_ntransfer_*_acc, and matrix_nturnover_*_acc are used to calculate the N capacity as the N steady state estimates in spin up.
  ! These variables are all state variables, saving accumulated N transfers during the calendar year.
    matrix_nalloc_leaf_acc        => cnveg_nitrogenstate_inst%matrix_nalloc_leaf_acc_patch        , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to leaf during this year
    matrix_nalloc_leafst_acc      => cnveg_nitrogenstate_inst%matrix_nalloc_leafst_acc_patch      , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to leaf storage during this year
    matrix_nalloc_froot_acc       => cnveg_nitrogenstate_inst%matrix_nalloc_froot_acc_patch       , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to fine root during this year
    matrix_nalloc_frootst_acc     => cnveg_nitrogenstate_inst%matrix_nalloc_frootst_acc_patch     , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to fine root storage during this year
    matrix_nalloc_livestem_acc    => cnveg_nitrogenstate_inst%matrix_nalloc_livestem_acc_patch    , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to live stem during this year
    matrix_nalloc_livestemst_acc  => cnveg_nitrogenstate_inst%matrix_nalloc_livestemst_acc_patch  , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to live stem storage during this year
    matrix_nalloc_deadstem_acc    => cnveg_nitrogenstate_inst%matrix_nalloc_deadstem_acc_patch    , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to dead stem during this year
    matrix_nalloc_deadstemst_acc  => cnveg_nitrogenstate_inst%matrix_nalloc_deadstemst_acc_patch  , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to dead stem storage during this year
    matrix_nalloc_livecroot_acc   => cnveg_nitrogenstate_inst%matrix_nalloc_livecroot_acc_patch   , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to live coarse root during this year
    matrix_nalloc_livecrootst_acc => cnveg_nitrogenstate_inst%matrix_nalloc_livecrootst_acc_patch , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to live coarse root storage during this year
    matrix_nalloc_deadcroot_acc   => cnveg_nitrogenstate_inst%matrix_nalloc_deadcroot_acc_patch   , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to dead coarse root during this year
    matrix_nalloc_deadcrootst_acc => cnveg_nitrogenstate_inst%matrix_nalloc_deadcrootst_acc_patch , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to dead coarse root storage during this year
    matrix_nalloc_grain_acc       => cnveg_nitrogenstate_inst%matrix_nalloc_grain_acc_patch       , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to grain during this year
    matrix_nalloc_grainst_acc     => cnveg_nitrogenstate_inst%matrix_nalloc_grainst_acc_patch     , & 
                                     ! In/Output: [real(r8) (:) ] (gN/m2/year) Input N allocated to grain storage during this year

    matrix_ntransfer_leafst_to_leafxf_acc            => cnveg_nitrogenstate_inst%matrix_ntransfer_leafst_to_leafxf_acc_patch           , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from leaf storage to leaf transfer pool during this year
    matrix_ntransfer_leafxf_to_leaf_acc              => cnveg_nitrogenstate_inst%matrix_ntransfer_leafxf_to_leaf_acc_patch             , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from leaf transfer to leaf pool during this year
    matrix_ntransfer_frootst_to_frootxf_acc          => cnveg_nitrogenstate_inst%matrix_ntransfer_frootst_to_frootxf_acc_patch         , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from fine root storage to fine root transfer pool during this year
    matrix_ntransfer_frootxf_to_froot_acc            => cnveg_nitrogenstate_inst%matrix_ntransfer_frootxf_to_froot_acc_patch           , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from fine root transfer to fine root pool during this year
    matrix_ntransfer_livestemst_to_livestemxf_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_livestemst_to_livestemxf_acc_patch   , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from live stem storage to live stem transfer pool during this year
    matrix_ntransfer_livestemxf_to_livestem_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_livestemxf_to_livestem_acc_patch     , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from live stem transfer to live stem pool during this year
    matrix_ntransfer_deadstemst_to_deadstemxf_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch   , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from dead stem storage to dead stem transfer pool during this year
    matrix_ntransfer_deadstemxf_to_deadstem_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_deadstemxf_to_deadstem_acc_patch     , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from dead stem transfer to dead stem pool during this year
    matrix_ntransfer_livecrootst_to_livecrootxf_acc  => cnveg_nitrogenstate_inst%matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from live coarese root storage to live coarese root transfer pool during this year
    matrix_ntransfer_livecrootxf_to_livecroot_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_livecrootxf_to_livecroot_acc_patch   , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from live coarese root transfer to live coarese root pool during this year
    matrix_ntransfer_deadcrootst_to_deadcrootxf_acc  => cnveg_nitrogenstate_inst%matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from dead coarse root storage to dead coarse root transfer pool during this year
    matrix_ntransfer_deadcrootxf_to_deadcroot_acc    => cnveg_nitrogenstate_inst%matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch   , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from dead coarse root transfer to dead coarse root pool during this year
    matrix_ntransfer_grainst_to_grainxf_acc          => cnveg_nitrogenstate_inst%matrix_ntransfer_grainst_to_grainxf_acc_patch         , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from grain storage to grain transfer pool during this year
    matrix_ntransfer_grainxf_to_grain_acc            => cnveg_nitrogenstate_inst%matrix_ntransfer_grainxf_to_grain_acc_patch           , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from grain transfer to grain pool during this year
    matrix_ntransfer_livestem_to_deadstem_acc        => cnveg_nitrogenstate_inst%matrix_ntransfer_livestem_to_deadstem_acc_patch       , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from live stem storage to dead stem transfer pool during this year
    matrix_ntransfer_livecroot_to_deadcroot_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_livecroot_to_deadcroot_acc_patch     , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from live coarse root to dead coarse root pool during this year
 
    matrix_ntransfer_retransn_to_leaf_acc            => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_leaf_acc_patch           , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to leaf pool during this year
    matrix_ntransfer_retransn_to_leafst_acc          => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_leafst_acc_patch         , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to leaf storage pool during this year
    matrix_ntransfer_retransn_to_froot_acc           => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_froot_acc_patch          , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to fine root pool during this year
    matrix_ntransfer_retransn_to_frootst_acc         => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_frootst_acc_patch        , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to fine root storage pool during this year
    matrix_ntransfer_retransn_to_livestem_acc        => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_livestem_acc_patch       , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to live stem pool during this year
    matrix_ntransfer_retransn_to_livestemst_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_livestemst_acc_patch     , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to livestem storage pool during this year
    matrix_ntransfer_retransn_to_deadstem_acc        => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_deadstem_acc_patch       , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to dead stem pool during this year
    matrix_ntransfer_retransn_to_deadstemst_acc      => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_deadstemst_acc_patch     , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to dead stem storage pool during this year
    matrix_ntransfer_retransn_to_livecroot_acc       => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_livecroot_acc_patch      , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to live coarse root pool during this year
    matrix_ntransfer_retransn_to_livecrootst_acc     => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_livecrootst_acc_patch    , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to live coarse root storage pool during this year
    matrix_ntransfer_retransn_to_deadcroot_acc       => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_deadcroot_acc_patch      , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to dead coarse root pool during this year
    matrix_ntransfer_retransn_to_deadcrootst_acc     => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_deadcrootst_acc_patch    , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to dead coarse root storage pool during this year
    matrix_ntransfer_retransn_to_grain_acc           => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_grain_acc_patch          , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to grain pool during this year
    matrix_ntransfer_retransn_to_grainst_acc         => cnveg_nitrogenstate_inst%matrix_ntransfer_retransn_to_grainst_acc_patch        , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from retranslocated N pool to grain storage pool during this year

    matrix_ntransfer_leaf_to_retransn_acc            => cnveg_nitrogenstate_inst%matrix_ntransfer_leaf_to_retransn_acc_patch           , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from leaf pool to retranslocated N pool during this year
    matrix_ntransfer_froot_to_retransn_acc           => cnveg_nitrogenstate_inst%matrix_ntransfer_froot_to_retransn_acc_patch          , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from fine root pool to retranslocated N pool during this year
    matrix_ntransfer_livestem_to_retransn_acc        => cnveg_nitrogenstate_inst%matrix_ntransfer_livestem_to_retransn_acc_patch       , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from live stem pool to retranslocated N pool during this year
    matrix_ntransfer_livecroot_to_retransn_acc       => cnveg_nitrogenstate_inst%matrix_ntransfer_livecroot_to_retransn_acc_patch      , & 
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N transfer from live coarse root pool to retranslocated N pool during this year

    matrix_nturnover_leaf_acc        => cnveg_nitrogenstate_inst%matrix_nturnover_leaf_acc_patch        , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from leaf
    matrix_nturnover_leafst_acc      => cnveg_nitrogenstate_inst%matrix_nturnover_leafst_acc_patch      , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from leaf storage
    matrix_nturnover_leafxf_acc      => cnveg_nitrogenstate_inst%matrix_nturnover_leafxf_acc_patch      , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from leaf transfer
    matrix_nturnover_froot_acc       => cnveg_nitrogenstate_inst%matrix_nturnover_froot_acc_patch       , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from fine root 
    matrix_nturnover_frootst_acc     => cnveg_nitrogenstate_inst%matrix_nturnover_frootst_acc_patch     , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from fine root storage
    matrix_nturnover_frootxf_acc     => cnveg_nitrogenstate_inst%matrix_nturnover_frootxf_acc_patch     , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from fine root transfer
    matrix_nturnover_livestem_acc    => cnveg_nitrogenstate_inst%matrix_nturnover_livestem_acc_patch    , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from live stem
    matrix_nturnover_livestemst_acc  => cnveg_nitrogenstate_inst%matrix_nturnover_livestemst_acc_patch  , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from live stem storage
    matrix_nturnover_livestemxf_acc  => cnveg_nitrogenstate_inst%matrix_nturnover_livestemxf_acc_patch  , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from live stem transfer
    matrix_nturnover_deadstem_acc    => cnveg_nitrogenstate_inst%matrix_nturnover_deadstem_acc_patch    , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from dead stem
    matrix_nturnover_deadstemst_acc  => cnveg_nitrogenstate_inst%matrix_nturnover_deadstemst_acc_patch  , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from dead stem storage
    matrix_nturnover_deadstemxf_acc  => cnveg_nitrogenstate_inst%matrix_nturnover_deadstemxf_acc_patch  , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from dead stem transfer
    matrix_nturnover_livecroot_acc   => cnveg_nitrogenstate_inst%matrix_nturnover_livecroot_acc_patch   , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from live coarse root
    matrix_nturnover_livecrootst_acc => cnveg_nitrogenstate_inst%matrix_nturnover_livecrootst_acc_patch , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from live coarse root storage
    matrix_nturnover_livecrootxf_acc => cnveg_nitrogenstate_inst%matrix_nturnover_livecrootxf_acc_patch , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from live coarse root transfer
    matrix_nturnover_deadcroot_acc   => cnveg_nitrogenstate_inst%matrix_nturnover_deadcroot_acc_patch   , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from dead coarse root
    matrix_nturnover_deadcrootst_acc => cnveg_nitrogenstate_inst%matrix_nturnover_deadcrootst_acc_patch , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from dead coarse root storage
    matrix_nturnover_deadcrootxf_acc => cnveg_nitrogenstate_inst%matrix_nturnover_deadcrootxf_acc_patch , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from dead coarse root transfer
    matrix_nturnover_grain_acc       => cnveg_nitrogenstate_inst%matrix_nturnover_grain_acc_patch       , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from grain
    matrix_nturnover_grainst_acc     => cnveg_nitrogenstate_inst%matrix_nturnover_grainst_acc_patch     , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from grain storage
    matrix_nturnover_grainxf_acc     => cnveg_nitrogenstate_inst%matrix_nturnover_grainxf_acc_patch     , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from grain transfer
    matrix_nturnover_retransn_acc    => cnveg_nitrogenstate_inst%matrix_nturnover_retransn_acc_patch    , &
                    ! In/Output: [real(r8) (:) ] (gN/m2/year) N turnover from retranslocated N pool

  ! *c0* variables save vegetation pool size at beginning of each year as a base for capacity calculation. For examples, 
  ! C turnover rate of pool KC_leaf (yr-1) is calculated by C turnover during the calendar year: matrix_cturnover_leaf_acc (gC/m2/yr) / leafc0 (gC/m2)
    leafc0              => cnveg_carbonstate_inst%leafc0_patch                , & ! In/Output: [real(r8) (:) ] (gC/m2) leaf C at begin of this year     
    leafc0_storage      => cnveg_carbonstate_inst%leafc0_storage_patch        , & ! In/Output: [real(r8) (:) ] (gC/m2) leaf storage C at begin of this year
    leafc0_xfer         => cnveg_carbonstate_inst%leafc0_xfer_patch           , & ! In/Output: [real(r8) (:) ] (gC/m2) leaf transfer C at begin of this year
    frootc0             => cnveg_carbonstate_inst%frootc0_patch               , & ! In/Output: [real(r8) (:) ] (gC/m2) fine root C at begin of this year
    frootc0_storage     => cnveg_carbonstate_inst%frootc0_storage_patch       , & ! In/Output: [real(r8) (:) ] (gC/m2) fine root storage C at begin of this year
    frootc0_xfer        => cnveg_carbonstate_inst%frootc0_xfer_patch          , & ! In/Output: [real(r8) (:) ] (gC/m2) fine root transfer C at begin of this year
    livestemc0          => cnveg_carbonstate_inst%livestemc0_patch            , & ! In/Output: [real(r8) (:) ] (gC/m2) live stem C at begin of this year
    livestemc0_storage  => cnveg_carbonstate_inst%livestemc0_storage_patch    , & ! In/Output: [real(r8) (:) ] (gC/m2) live stem storage C at begin of this year
    livestemc0_xfer     => cnveg_carbonstate_inst%livestemc0_xfer_patch       , & ! In/Output: [real(r8) (:) ] (gC/m2) live stem transfer C at begin of this year
    deadstemc0          => cnveg_carbonstate_inst%deadstemc0_patch            , & ! In/Output: [real(r8) (:) ] (gC/m2) dead stem C at begin of this year
    deadstemc0_storage  => cnveg_carbonstate_inst%deadstemc0_storage_patch    , & ! In/Output: [real(r8) (:) ] (gC/m2) dead stem storage C at begin of this year
    deadstemc0_xfer     => cnveg_carbonstate_inst%deadstemc0_xfer_patch       , & ! In/Output: [real(r8) (:) ] (gC/m2) dead stem transfer C at begin of this year
    livecrootc0         => cnveg_carbonstate_inst%livecrootc0_patch           , & ! In/Output: [real(r8) (:) ] (gC/m2) live coarse root C at begin of this year
    livecrootc0_storage => cnveg_carbonstate_inst%livecrootc0_storage_patch   , & ! In/Output: [real(r8) (:) ] (gC/m2) live coarse root storage C at begin of this year
    livecrootc0_xfer    => cnveg_carbonstate_inst%livecrootc0_xfer_patch      , & ! In/Output: [real(r8) (:) ] (gC/m2) live coarse root transfer C at begin of this year
    deadcrootc0         => cnveg_carbonstate_inst%deadcrootc0_patch           , & ! In/Output: [real(r8) (:) ] (gC/m2) dead coarse root C at begin of this year
    deadcrootc0_storage => cnveg_carbonstate_inst%deadcrootc0_storage_patch   , & ! In/Output: [real(r8) (:) ] (gC/m2) dead coarse root storage C at begin of this year
    deadcrootc0_xfer    => cnveg_carbonstate_inst%deadcrootc0_xfer_patch      , & ! In/Output: [real(r8) (:) ] (gC/m2) dead coarse root transfer C at begin of this year
    reproc0             => cnveg_carbonstate_inst%reproc0_patch               , & ! In/Output: [real(r8) (:) ] (gC/m2) grain C at begin of this year
    reproc0_storage     => cnveg_carbonstate_inst%reproc0_storage_patch       , & ! In/Output: [real(r8) (:) ] (gC/m2) grain storage C at begin of this year
    reproc0_xfer        => cnveg_carbonstate_inst%reproc0_xfer_patch          , & ! In/Output: [real(r8) (:) ] (gC/m2) grain transfer C at begin of this year

  ! *n0* variables save vegetation pool size at beginning of each year as a base for capacity calculation. For examples, 
  ! N turnover rate of pool KN_leaf (yr-1) is calculated by N turnover during the calendar year matrix_nturnover_leaf_acc (gN/m2/yr) / leafn0 (gN/m2)
    leafn0              => cnveg_nitrogenstate_inst%leafn0_patch              , & ! In/Output: [real(r8) (:) ] (gN/m2) leaf N at begin of this year
    leafn0_storage      => cnveg_nitrogenstate_inst%leafn0_storage_patch      , & ! In/Output: [real(r8) (:) ] (gN/m2) leaf storage N at begin of this year
    leafn0_xfer         => cnveg_nitrogenstate_inst%leafn0_xfer_patch         , & ! In/Output: [real(r8) (:) ] (gN/m2) leaf transfer N at begin of this year
    frootn0             => cnveg_nitrogenstate_inst%frootn0_patch             , & ! In/Output: [real(r8) (:) ] (gN/m2) fine root N at begin of this year
    frootn0_storage     => cnveg_nitrogenstate_inst%frootn0_storage_patch     , & ! In/Output: [real(r8) (:) ] (gN/m2) fine root storage N at begin of this year
    frootn0_xfer        => cnveg_nitrogenstate_inst%frootn0_xfer_patch        , & ! In/Output: [real(r8) (:) ] (gN/m2) fine root transfer N at begin of this year
    livestemn0          => cnveg_nitrogenstate_inst%livestemn0_patch          , & ! In/Output: [real(r8) (:) ] (gN/m2) live stem N at begin of this year
    livestemn0_storage  => cnveg_nitrogenstate_inst%livestemn0_storage_patch  , & ! In/Output: [real(r8) (:) ] (gN/m2) live stem storage N at begin of this year
    livestemn0_xfer     => cnveg_nitrogenstate_inst%livestemn0_xfer_patch     , & ! In/Output: [real(r8) (:) ] (gN/m2) live stem transfer N at begin of this year
    deadstemn0          => cnveg_nitrogenstate_inst%deadstemn0_patch          , & ! In/Output: [real(r8) (:) ] (gN/m2) dead stem N at begin of this year
    deadstemn0_storage  => cnveg_nitrogenstate_inst%deadstemn0_storage_patch  , & ! In/Output: [real(r8) (:) ] (gN/m2) dead stem storage N at begin of this year
    deadstemn0_xfer     => cnveg_nitrogenstate_inst%deadstemn0_xfer_patch     , & ! In/Output: [real(r8) (:) ] (gN/m2) dead stem transfer N at begin of this year
    livecrootn0         => cnveg_nitrogenstate_inst%livecrootn0_patch         , & ! In/Output: [real(r8) (:) ] (gN/m2) live coarse root N at begin of this year
    livecrootn0_storage => cnveg_nitrogenstate_inst%livecrootn0_storage_patch , & ! In/Output: [real(r8) (:) ] (gN/m2) live coarse root storage N at begin of this year
    livecrootn0_xfer    => cnveg_nitrogenstate_inst%livecrootn0_xfer_patch    , & ! In/Output: [real(r8) (:) ] (gN/m2) live coarse root transfer N at begin of this year
    deadcrootn0         => cnveg_nitrogenstate_inst%deadcrootn0_patch         , & ! In/Output: [real(r8) (:) ] (gN/m2) dead coarse root N at begin of this year
    deadcrootn0_storage => cnveg_nitrogenstate_inst%deadcrootn0_storage_patch , & ! In/Output: [real(r8) (:) ] (gN/m2) dead coarse root storage N at begin of this year
    deadcrootn0_xfer    => cnveg_nitrogenstate_inst%deadcrootn0_xfer_patch    , & ! In/Output: [real(r8) (:) ] (gN/m2) dead coarse root transfer N at begin of this year
    repron0             => cnveg_nitrogenstate_inst%repron0_patch             , & ! In/Output: [real(r8) (:) ] (gN/m2) grain N at begin of this year
    repron0_storage     => cnveg_nitrogenstate_inst%repron0_storage_patch     , & ! In/Output: [real(r8) (:) ] (gN/m2) grain storage N at begin of this year
    repron0_xfer        => cnveg_nitrogenstate_inst%repron0_xfer_patch        , & ! In/Output: [real(r8) (:) ] (gN/m2) grain transfer N at begin of this year
    retransn0           => cnveg_nitrogenstate_inst%retransn0_patch             & ! In/Output: [real(r8) (:) ] (gN/m2) plant retranslocated N at begin of this year
    )
sd: associate(                          &

  ! Following variables save the C and N transfer rate of different processes at current time step. 
  ! Eg. ph: phenology, gm: gap mortality (including harvest), fi: fire.
    matrix_alloc       => cnveg_carbonflux_inst%matrix_alloc_patch         , & ! Input:  [real(r8) (:,:)] (gC/gC) input C allocation matrix, updated in NutrientCompetitionFlexibleCNMod or NutrientCompetitionCLM45defaultMod
    matrix_nalloc      => cnveg_nitrogenflux_inst%matrix_nalloc_patch      , & ! Input:  [real(r8) (:,:)] (gC/gC) input N allocation matrix, updated in NutrientCompetitionFlexibleCNMod or NutrientCompetitionCLM45defaultMod
    matrix_phtransfer  => cnveg_carbonflux_inst%matrix_phtransfer_patch    , & ! Input:  [real(r8) (:,:)] (gC/m2/s) C transfer rate from phenology processes, updated in CNPhenology
    matrix_gmtransfer  => cnveg_carbonflux_inst%matrix_gmtransfer_patch    , & ! Input:  [real(r8) (:,:)] (gC/m2/s) C transfer rate from gap mortality processes, updated in CNGapMortality
    matrix_fitransfer  => cnveg_carbonflux_inst%matrix_fitransfer_patch    , & ! Input:  [real(r8) (:,:)] (gC/m2/s) C transfer rate from fire processes, updated in CNFireBaseMod or CNFireLi2014Mod
    matrix_phturnover  => cnveg_carbonflux_inst%matrix_phturnover_patch    , & ! Output: [real(r8) (:,:)] (gC/m2/step) C turnover rate from phenology processes, updated in CNVegMatrixMod and dynHarvestMod
    matrix_gmturnover  => cnveg_carbonflux_inst%matrix_gmturnover_patch    , & ! Output: [real(r8) (:,:)] (gC/m2/step) C turnover rate from gap mortality processe, updated in CNVegMatrixMods
    matrix_fiturnover  => cnveg_carbonflux_inst%matrix_fiturnover_patch    , & ! Output: [real(r8) (:,:)] (gC/m2/step) C turnover rate from fire processe, updated in CNVegMatrixMods

    matrix_nphtransfer => cnveg_nitrogenflux_inst%matrix_nphtransfer_patch , & ! Input:  [real(r8) (:,:)] (gN/m2/s) N transfer rate from phenology processes, updated in CNPhenology and (NutrientCompetitionFlexibleCNMod or NutrientCompetitionCLM45defaultMod)
    matrix_ngmtransfer => cnveg_nitrogenflux_inst%matrix_ngmtransfer_patch , & ! Input:  [real(r8) (:,:)] (gN/m2/s) N transfer rate from gap mortality processes, updated in CNGapMortality and dynHarvestMod
    matrix_nfitransfer => cnveg_nitrogenflux_inst%matrix_nfitransfer_patch , & ! Input:  [real(r8) (:,:)] (gN/m2/s) N transfer rate from fire processes, updated in CNFireBaseMod or CNFireLi2014Mod
    matrix_nphturnover => cnveg_nitrogenflux_inst%matrix_nphturnover_patch , & ! Output: [real(r8) (:,:)] (gN/m2/step) N turnover rate from phenology processes, updated in CNVegMatrixMod
    matrix_ngmturnover => cnveg_nitrogenflux_inst%matrix_ngmturnover_patch , & ! Output: [real(r8) (:,:)] (gN/m2/step) N turnover rate from gap mortality processes, updated in CNVegMatrixMod
    matrix_nfiturnover => cnveg_nitrogenflux_inst%matrix_nfiturnover_patch , & ! Output: [real(r8) (:,:)] (gN/m2/step) N turnover rate from fire processes, updated in CNVegMatrixMod

    matrix_Cinput      => cnveg_carbonflux_inst%matrix_Cinput_patch        , & ! Input:  [real(r8) (:)] (gC/m2/s) C input to vegetation, updated in NutrientCompetitionFlexibleCNMod or NutrientCompetitionCLM45defaultMod
    matrix_C13input    => cnveg_carbonflux_inst%matrix_C13input_patch      , & ! Input:  [real(r8) (:)] (gC/m2/s) C13 input to vegetation, updated in NutrientCompetitionFlexibleCNMod or NutrientCompetitionCLM45defaultMod
    matrix_C14input    => cnveg_carbonflux_inst%matrix_C14input_patch      , & ! Input:  [real(r8) (:)] (gC/m2/s) C14 input to vegetation, updated in NutrientCompetitionFlexibleCNMod or NutrientCompetitionCLM45defaultMod
    matrix_Ninput      => cnveg_nitrogenflux_inst%matrix_Ninput_patch      , & ! Input:  [real(r8) (:)] (gN/m2/s) N input to vegetation, updated in NutrientCompetitionFlexibleCNMod or NutrientCompetitionCLM45defaultMod

  ! Doners and receivers of all transfers from different processes have been prescribed in following variables:
    doner_phc          => cnveg_carbonflux_inst%matrix_phtransfer_doner_patch        , & ! Input:  [integer (:)] Doners of phenology related C transfer
    receiver_phc       => cnveg_carbonflux_inst%matrix_phtransfer_receiver_patch     , & ! Input:  [integer (:)] Receiver of phenology related C transfer
    doner_gmc          => cnveg_carbonflux_inst%matrix_gmtransfer_doner_patch        , & ! Input:  [integer (:)] Doners of gap mortality related C transfer
    receiver_gmc       => cnveg_carbonflux_inst%matrix_gmtransfer_receiver_patch     , & ! Input:  [integer (:)] Receiver of gap mortality related C transfer
    doner_fic          => cnveg_carbonflux_inst%matrix_fitransfer_doner_patch        , & ! Input:  [integer (:)] Doners of fire related C transfer
    receiver_fic       => cnveg_carbonflux_inst%matrix_fitransfer_receiver_patch     , & ! Input:  [integer (:)] Receiver of fire related C transfer
    doner_phn          => cnveg_nitrogenflux_inst%matrix_nphtransfer_doner_patch     , & ! Input:  [integer (:)] Doners of phenology related N transfer
    receiver_phn       => cnveg_nitrogenflux_inst%matrix_nphtransfer_receiver_patch  , & ! Input:  [integer (:)] Receiver of phenology related N transfer
    doner_gmn          => cnveg_nitrogenflux_inst%matrix_ngmtransfer_doner_patch     , & ! Input:  [integer (:)] Doners of gap mortality related N transfer
    receiver_gmn       => cnveg_nitrogenflux_inst%matrix_ngmtransfer_receiver_patch  , & ! Input:  [integer (:)] Receiver of gap mortality related N transfer
    doner_fin          => cnveg_nitrogenflux_inst%matrix_nfitransfer_doner_patch     , & ! Input:  [integer (:)] Doners of fire related N transfer
    receiver_fin       => cnveg_nitrogenflux_inst%matrix_nfitransfer_receiver_patch  , & ! Input:  [integer (:)] Receiver of fire related N transfer

  ! Index of each processes related C transfers. See subroutine InitTransfer in CNVegCarbonFluxType.F90 for details.
    ileafst_to_ileafxf_phc            => cnveg_carbonflux_inst%ileafst_to_ileafxf_ph           , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from leaf storage pool to leaf transfer pool
    ileafxf_to_ileaf_phc              => cnveg_carbonflux_inst%ileafxf_to_ileaf_ph             , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from leaf transfer pool to leaf pool
    ifrootst_to_ifrootxf_phc          => cnveg_carbonflux_inst%ifrootst_to_ifrootxf_ph         , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from fine root storage pool to fine root transfer pool
    ifrootxf_to_ifroot_phc            => cnveg_carbonflux_inst%ifrootxf_to_ifroot_ph           , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from fine root transfer pool to fine root pool
    ilivestemst_to_ilivestemxf_phc    => cnveg_carbonflux_inst%ilivestemst_to_ilivestemxf_ph   , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from live stem storage pool to live stem transfer pool
    ilivestemxf_to_ilivestem_phc      => cnveg_carbonflux_inst%ilivestemxf_to_ilivestem_ph     , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from live stem transfer pool to live stem pool
    ideadstemst_to_ideadstemxf_phc    => cnveg_carbonflux_inst%ideadstemst_to_ideadstemxf_ph   , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
    ideadstemxf_to_ideadstem_phc      => cnveg_carbonflux_inst%ideadstemxf_to_ideadstem_ph     , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from dead stem transfer pool to dead stem pool
    ilivecrootst_to_ilivecrootxf_phc  => cnveg_carbonflux_inst%ilivecrootst_to_ilivecrootxf_ph , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
    ilivecrootxf_to_ilivecroot_phc    => cnveg_carbonflux_inst%ilivecrootxf_to_ilivecroot_ph   , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
    ideadcrootst_to_ideadcrootxf_phc  => cnveg_carbonflux_inst%ideadcrootst_to_ideadcrootxf_ph , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
    ideadcrootxf_to_ideadcroot_phc    => cnveg_carbonflux_inst%ideadcrootxf_to_ideadcroot_ph   , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
    ilivestem_to_ideadstem_phc        => cnveg_carbonflux_inst%ilivestem_to_ideadstem_ph       , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to dead stem pool
    ilivecroot_to_ideadcroot_phc      => cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_ph     , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from live coarse root pool to dead coarse root pool
    ileaf_to_iout_phc                 => cnveg_carbonflux_inst%ileaf_to_iout_ph                , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from leaf pool to outside of vegetation pools 
    ifroot_to_iout_phc                => cnveg_carbonflux_inst%ifroot_to_iout_ph               , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from fine root pool to outside of vegetation pools
    ilivestem_to_iout_phc             => cnveg_carbonflux_inst%ilivestem_to_iout_ph            , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from live stem pool to outside of vegetation pools
    igrain_to_iout_phc                => cnveg_carbonflux_inst%igrain_to_iout_ph               , & 
                          ! Input: [integer (:)] Index of phenology related C transfer from grain pool to outside of vegetation pools
    ileaf_to_iout_gmc                 => cnveg_carbonflux_inst%ileaf_to_iout_gm                , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from leaf pool to outside of vegetation pools
    ileafst_to_iout_gmc               => cnveg_carbonflux_inst%ileafst_to_iout_gm              , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from leaf storage pool to outside of vegetation pools
    ileafxf_to_iout_gmc               => cnveg_carbonflux_inst%ileafxf_to_iout_gm              , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from leaf transfer pool to outside of vegetation pools
    ifroot_to_iout_gmc                => cnveg_carbonflux_inst%ifroot_to_iout_gm               , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from fine root pool to outside of vegetation pools
    ifrootst_to_iout_gmc              => cnveg_carbonflux_inst%ifrootst_to_iout_gm             , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from fine root storage pool to outside of vegetation pools
    ifrootxf_to_iout_gmc              => cnveg_carbonflux_inst%ifrootxf_to_iout_gm             , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from fine root transfer pool to outside of vegetation pools
    ilivestem_to_iout_gmc             => cnveg_carbonflux_inst%ilivestem_to_iout_gm            , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from live stem pool to outside of vegetation pools
    ilivestemst_to_iout_gmc           => cnveg_carbonflux_inst%ilivestemst_to_iout_gm          , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from live stem storage pool to outside of vegetation pools
    ilivestemxf_to_iout_gmc           => cnveg_carbonflux_inst%ilivestemxf_to_iout_gm          , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from live stem transfer pool to outside of vegetation pools
    ideadstem_to_iout_gmc             => cnveg_carbonflux_inst%ideadstem_to_iout_gm            , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from dead stem pool to outside of vegetation pools
    ideadstemst_to_iout_gmc           => cnveg_carbonflux_inst%ideadstemst_to_iout_gm          , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from dead stem storage pool to outside of vegetation pools
    ideadstemxf_to_iout_gmc           => cnveg_carbonflux_inst%ideadstemxf_to_iout_gm          , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from dead stem transfer pool to outside of vegetation pools
    ilivecroot_to_iout_gmc            => cnveg_carbonflux_inst%ilivecroot_to_iout_gm           , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from live coarse root pool to outside of vegetation pools
    ilivecrootst_to_iout_gmc          => cnveg_carbonflux_inst%ilivecrootst_to_iout_gm         , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from live coarse root storage pool to outside of vegetation pools
    ilivecrootxf_to_iout_gmc          => cnveg_carbonflux_inst%ilivecrootxf_to_iout_gm         , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from live coarse root transfer pool to outside of vegetation pools
    ideadcroot_to_iout_gmc            => cnveg_carbonflux_inst%ideadcroot_to_iout_gm           , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from dead coarse root pool to outside of vegetation pools
    ideadcrootst_to_iout_gmc          => cnveg_carbonflux_inst%ideadcrootst_to_iout_gm         , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from dead coarse root storage pool to outside of vegetation pools
    ideadcrootxf_to_iout_gmc          => cnveg_carbonflux_inst%ideadcrootxf_to_iout_gm         , & 
                          ! Input: [integer (:)] Index of gap mortality related C transfer from dead coarse root transfer pool to outside of vegetation pools
    ileaf_to_iout_fic                 => cnveg_carbonflux_inst%ileaf_to_iout_fi                , & 
                          ! Input: [integer (:)] Index of fire related C transfer from leaf pool to outside of vegetation pools
    ileafst_to_iout_fic               => cnveg_carbonflux_inst%ileafst_to_iout_fi              , & 
                          ! Input: [integer (:)] Index of fire related C transfer from leaf storage pool to outside of vegetation pools
    ileafxf_to_iout_fic               => cnveg_carbonflux_inst%ileafxf_to_iout_fi              , & 
                          ! Input: [integer (:)] Index of fire related C transfer from leaf transfer pool to outside of vegetation pools
    ifroot_to_iout_fic                => cnveg_carbonflux_inst%ifroot_to_iout_fi               , & 
                          ! Input: [integer (:)] Index of fire related C transfer from fine root pool to outside of vegetation pools
    ifrootst_to_iout_fic              => cnveg_carbonflux_inst%ifrootst_to_iout_fi             , & 
                          ! Input: [integer (:)] Index of fire related C transfer from fine root storage pool to outside of vegetation pools
    ifrootxf_to_iout_fic              => cnveg_carbonflux_inst%ifrootxf_to_iout_fi             , & 
                          ! Input: [integer (:)] Index of fire related C transfer from fine root transfer pool to outside of vegetation pools
    ilivestem_to_iout_fic             => cnveg_carbonflux_inst%ilivestem_to_iout_fi            , & 
                          ! Input: [integer (:)] Index of fire related C transfer from live stem pool to outside of vegetation pools
    ilivestemst_to_iout_fic           => cnveg_carbonflux_inst%ilivestemst_to_iout_fi          , & 
                          ! Input: [integer (:)] Index of fire related C transfer from live stem storage pool to outside of vegetation pools
    ilivestemxf_to_iout_fic           => cnveg_carbonflux_inst%ilivestemxf_to_iout_fi          , & 
                          ! Input: [integer (:)] Index of fire related C transfer from live stem transfer pool to outside of vegetation pools
    ideadstem_to_iout_fic             => cnveg_carbonflux_inst%ideadstem_to_iout_fi            , & 
                          ! Input: [integer (:)] Index of fire related C transfer from dead stem pool to outside of vegetation pools
    ideadstemst_to_iout_fic           => cnveg_carbonflux_inst%ideadstemst_to_iout_fi          , & 
                          ! Input: [integer (:)] Index of fire related C transfer from dead stem storage pool to outside of vegetation pools
    ideadstemxf_to_iout_fic           => cnveg_carbonflux_inst%ideadstemxf_to_iout_fi          , & 
                          ! Input: [integer (:)] Index of fire related C transfer from dead stem transfer pool to outside of vegetation pools
    ilivecroot_to_iout_fic            => cnveg_carbonflux_inst%ilivecroot_to_iout_fi           , & 
                          ! Input: [integer (:)] Index of fire related C transfer from live coarse root pool to outside of vegetation pools
    ilivecrootst_to_iout_fic          => cnveg_carbonflux_inst%ilivecrootst_to_iout_fi         , & 
                          ! Input: [integer (:)] Index of fire related C transfer from live coarse root storage pool to outside of vegetation pools
    ilivecrootxf_to_iout_fic          => cnveg_carbonflux_inst%ilivecrootxf_to_iout_fi         , & 
                          ! Input: [integer (:)] Index of fire related C transfer from live coarse root transfer pool to outside of vegetation pools
    ideadcroot_to_iout_fic            => cnveg_carbonflux_inst%ideadcroot_to_iout_fi           , & 
                          ! Input: [integer (:)] Index of fire related C transfer from dead coarse root pool to outside of vegetation pools
    ideadcrootst_to_iout_fic          => cnveg_carbonflux_inst%ideadcrootst_to_iout_fi         , & 
                          ! Input: [integer (:)] Index of fire related C transfer from dead coarse root storage pool to outside of vegetation pools
    ideadcrootxf_to_iout_fic          => cnveg_carbonflux_inst%ideadcrootxf_to_iout_fi         , & 
                          ! Input: [integer (:)] Index of fire related C transfer from dead coarse root transfer pool to outside of vegetation pools
    ilivestem_to_ideadstem_fic        => cnveg_carbonflux_inst%ilivestem_to_ideadstem_fi       , & 
                          ! Input: [integer (:)] Index of fire related C transfer from live stem pool to dead stem pool
    ilivecroot_to_ideadcroot_fic      => cnveg_carbonflux_inst%ilivecroot_to_ideadcroot_fi     , & 
                          ! Input: [integer (:)] Index of fire related C transfer from live coarse root pool to dead coarse root pool
  ! Index of each processes related N transfers. See subroutine InitTransfer in CNVegNitrogenFluxType.F90 for details.
    ileafst_to_ileafxf_phn            => cnveg_nitrogenflux_inst%ileafst_to_ileafxf_ph          , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from leaf storage pool to leaf transfer pool
    ileafxf_to_ileaf_phn              => cnveg_nitrogenflux_inst%ileafxf_to_ileaf_ph            , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from leaf transfer pool to leaf pool
    ifrootst_to_ifrootxf_phn          => cnveg_nitrogenflux_inst%ifrootst_to_ifrootxf_ph        , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from fine root storage pool to fine root transfer pool
    ifrootxf_to_ifroot_phn            => cnveg_nitrogenflux_inst%ifrootxf_to_ifroot_ph          , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from fine root transfer pool to fine root pool
    ilivestemst_to_ilivestemxf_phn    => cnveg_nitrogenflux_inst%ilivestemst_to_ilivestemxf_ph  , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live stem storage pool to live stem transfer pool
    ilivestemxf_to_ilivestem_phn      => cnveg_nitrogenflux_inst%ilivestemxf_to_ilivestem_ph    , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live stem transfer pool to live stem pool
    ideadstemst_to_ideadstemxf_phn    => cnveg_nitrogenflux_inst%ideadstemst_to_ideadstemxf_ph  , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from dead stem storage pool to dead stem transfer pool
    ideadstemxf_to_ideadstem_phn      => cnveg_nitrogenflux_inst%ideadstemxf_to_ideadstem_ph    , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from dead stem transfer pool to dead stem pool
    ilivecrootst_to_ilivecrootxf_phn  => cnveg_nitrogenflux_inst%ilivecrootst_to_ilivecrootxf_ph, & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live coarse root storage pool to live coarse root transfer pool
    ilivecrootxf_to_ilivecroot_phn    => cnveg_nitrogenflux_inst%ilivecrootxf_to_ilivecroot_ph  , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live coarse root transfer pool to live coarse root pool
    ideadcrootst_to_ideadcrootxf_phn  => cnveg_nitrogenflux_inst%ideadcrootst_to_ideadcrootxf_ph, & 
                          ! Input: [integer (:)] Index of phenology related N transfer from dead coarse root storage pool to dead coarse root transfer pool
    ideadcrootxf_to_ideadcroot_phn    => cnveg_nitrogenflux_inst%ideadcrootxf_to_ideadcroot_ph  , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from dead coarse root transfer pool to dead coarse root pool
    ilivestem_to_ideadstem_phn        => cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_ph      , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live stem pool to dead stem pool
    ilivecroot_to_ideadcroot_phn      => cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_ph    , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live coarse root pool to dead coarse root pool
    ileaf_to_iout_phn                 => cnveg_nitrogenflux_inst%ileaf_to_iout_ph               , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from leaf pool to outside of vegetation pools
    ifroot_to_iout_phn                => cnveg_nitrogenflux_inst%ifroot_to_iout_ph              , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from fine root pool to outside of vegetation pools
    ilivestem_to_iout_phn             => cnveg_nitrogenflux_inst%ilivestem_to_iout_ph           , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live stem pool to outside of vegetation pools
    iretransn_to_iout_phn             => cnveg_nitrogenflux_inst%iretransn_to_iout_ph           , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to outside of vegetation pools
    igrain_to_iout_phn                => cnveg_nitrogenflux_inst%igrain_to_iout_ph              , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from grain pool to outside of vegetation pools
    ileaf_to_iretransn_phn            => cnveg_nitrogenflux_inst%ileaf_to_iretransn_ph          , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from leaf pool to retranslocated N pool
    ifroot_to_iretransn_phn           => cnveg_nitrogenflux_inst%ifroot_to_iretransn_ph         , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from fine root pool to retranslocated N pool
    ilivestem_to_iretransn_phn        => cnveg_nitrogenflux_inst%ilivestem_to_iretransn_ph      , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live stem pool to retranslocated N pool
    ilivecroot_to_iretransn_phn       => cnveg_nitrogenflux_inst%ilivecroot_to_iretransn_ph     , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from live coarse root pool to retranslocated N pool
    iretransn_to_ileaf_phn            => cnveg_nitrogenflux_inst%iretransn_to_ileaf_ph          , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to leaf pool
    iretransn_to_ileafst_phn          => cnveg_nitrogenflux_inst%iretransn_to_ileafst_ph        , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to leaf storage pool
    iretransn_to_ifroot_phn           => cnveg_nitrogenflux_inst%iretransn_to_ifroot_ph         , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to fine root pool
    iretransn_to_ifrootst_phn         => cnveg_nitrogenflux_inst%iretransn_to_ifrootst_ph       , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to fine root storage pool
    iretransn_to_ilivestem_phn        => cnveg_nitrogenflux_inst%iretransn_to_ilivestem_ph      , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to live stem pool
    iretransn_to_ilivestemst_phn      => cnveg_nitrogenflux_inst%iretransn_to_ilivestemst_ph    , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to live stem storage pool
    iretransn_to_ideadstem_phn        => cnveg_nitrogenflux_inst%iretransn_to_ideadstem_ph      , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to dead stem pool
    iretransn_to_ideadstemst_phn      => cnveg_nitrogenflux_inst%iretransn_to_ideadstemst_ph    , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to dead stem storage pool
    iretransn_to_ilivecroot_phn       => cnveg_nitrogenflux_inst%iretransn_to_ilivecroot_ph     , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to live coarse root pool
    iretransn_to_ilivecrootst_phn     => cnveg_nitrogenflux_inst%iretransn_to_ilivecrootst_ph   , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to live coarse root storage pool
    iretransn_to_ideadcroot_phn       => cnveg_nitrogenflux_inst%iretransn_to_ideadcroot_ph     , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to dead coarse root pool
    iretransn_to_ideadcrootst_phn     => cnveg_nitrogenflux_inst%iretransn_to_ideadcrootst_ph   , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to dead coarse root storage pool
    iretransn_to_igrain_phn           => cnveg_nitrogenflux_inst%iretransn_to_igrain_ph         , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to grain pool
    iretransn_to_igrainst_phn         => cnveg_nitrogenflux_inst%iretransn_to_igrainst_ph       , & 
                          ! Input: [integer (:)] Index of phenology related N transfer from retranslocated N pool to grain storage pool
    ileaf_to_iout_gmn                 => cnveg_nitrogenflux_inst%ileaf_to_iout_gm               , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from leaf pool to outside of vegetation pools
    ileafst_to_iout_gmn               => cnveg_nitrogenflux_inst%ileafst_to_iout_gm             , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from leaf storage pool to outside of vegetation pools
    ileafxf_to_iout_gmn               => cnveg_nitrogenflux_inst%ileafxf_to_iout_gm             , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from leaf transfer pool to outside of vegetation pools
    ifroot_to_iout_gmn                => cnveg_nitrogenflux_inst%ifroot_to_iout_gm              , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from fine root pool to outside of vegetation pools
    ifrootst_to_iout_gmn              => cnveg_nitrogenflux_inst%ifrootst_to_iout_gm            , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from fine root storage pool to outside of vegetation pools
    ifrootxf_to_iout_gmn              => cnveg_nitrogenflux_inst%ifrootxf_to_iout_gm            , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from fine root transfer pool to outside of vegetation pools
    ilivestem_to_iout_gmn             => cnveg_nitrogenflux_inst%ilivestem_to_iout_gm           , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from live stem pool to outside of vegetation pools
    ilivestemst_to_iout_gmn           => cnveg_nitrogenflux_inst%ilivestemst_to_iout_gm         , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from live stem storage pool to outside of vegetation pools
    ilivestemxf_to_iout_gmn           => cnveg_nitrogenflux_inst%ilivestemxf_to_iout_gm         , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from live stem transfer pool to outside of vegetation pools
    ideadstem_to_iout_gmn             => cnveg_nitrogenflux_inst%ideadstem_to_iout_gm           , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from dead stem pool to outside of vegetation pools
    ideadstemst_to_iout_gmn           => cnveg_nitrogenflux_inst%ideadstemst_to_iout_gm         , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from dead stem storage pool to outside of vegetation pools
    ideadstemxf_to_iout_gmn           => cnveg_nitrogenflux_inst%ideadstemxf_to_iout_gm         , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from dead stem transfer pool to outside of vegetation pools
    ilivecroot_to_iout_gmn            => cnveg_nitrogenflux_inst%ilivecroot_to_iout_gm          , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from live coarse root pool to outside of vegetation pools
    ilivecrootst_to_iout_gmn          => cnveg_nitrogenflux_inst%ilivecrootst_to_iout_gm        , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from live coarse root storage pool to outside of vegetation pools
    ilivecrootxf_to_iout_gmn          => cnveg_nitrogenflux_inst%ilivecrootxf_to_iout_gm        , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from live coarse root transfer pool to outside of vegetation pools
    ideadcroot_to_iout_gmn            => cnveg_nitrogenflux_inst%ideadcroot_to_iout_gm          , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from dead coarse root pool to outside of vegetation pools
    ideadcrootst_to_iout_gmn          => cnveg_nitrogenflux_inst%ideadcrootst_to_iout_gm        , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from dead coarse root storage pool to outside of vegetation pools
    ideadcrootxf_to_iout_gmn          => cnveg_nitrogenflux_inst%ideadcrootxf_to_iout_gm        , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from dead coarse root transfer pool to outside of vegetation pools
    iretransn_to_iout_gmn             => cnveg_nitrogenflux_inst%iretransn_to_iout_gm           , & 
                          ! Input: [integer (:)] Index of gap mortality related N transfer from retranslocated N pool to outside of vegetation pools
    ileaf_to_iout_fin                 => cnveg_nitrogenflux_inst%ileaf_to_iout_fi               , & 
                          ! Input: [integer (:)] Index of fire related N transfer from leaf pool to outside of vegetation pools
    ileafst_to_iout_fin               => cnveg_nitrogenflux_inst%ileafst_to_iout_fi             , & 
                          ! Input: [integer (:)] Index of fire related N transfer from leaf storage pool to outside of vegetation pools
    ileafxf_to_iout_fin               => cnveg_nitrogenflux_inst%ileafxf_to_iout_fi             , & 
                          ! Input: [integer (:)] Index of fire related N transfer from leaf transfer pool to outside of vegetation pools
    ifroot_to_iout_fin                => cnveg_nitrogenflux_inst%ifroot_to_iout_fi              , & 
                          ! Input: [integer (:)] Index of fire related N transfer from fine root pool to outside of vegetation pools
    ifrootst_to_iout_fin              => cnveg_nitrogenflux_inst%ifrootst_to_iout_fi            , & 
                          ! Input: [integer (:)] Index of fire related N transfer from fine root storage pool to outside of vegetation pools
    ifrootxf_to_iout_fin              => cnveg_nitrogenflux_inst%ifrootxf_to_iout_fi            , & 
                          ! Input: [integer (:)] Index of fire related N transfer from fine transfer pool to outside of vegetation pools
    ilivestem_to_iout_fin             => cnveg_nitrogenflux_inst%ilivestem_to_iout_fi           , & 
                          ! Input: [integer (:)] Index of fire related N transfer from live stem pool to outside of vegetation pools
    ilivestemst_to_iout_fin           => cnveg_nitrogenflux_inst%ilivestemst_to_iout_fi         , & 
                          ! Input: [integer (:)] Index of fire related N transfer from live stem storage pool to outside of vegetation pools
    ilivestemxf_to_iout_fin           => cnveg_nitrogenflux_inst%ilivestemxf_to_iout_fi         , & 
                          ! Input: [integer (:)] Index of fire related N transfer from live stem transfer pool to outside of vegetation pools
    ideadstem_to_iout_fin             => cnveg_nitrogenflux_inst%ideadstem_to_iout_fi           , & 
                          ! Input: [integer (:)] Index of fire related N transfer from dead stem pool to outside of vegetation pools
    ideadstemst_to_iout_fin           => cnveg_nitrogenflux_inst%ideadstemst_to_iout_fi         , & 
                          ! Input: [integer (:)] Index of fire related N transfer from dead stem storage pool to outside of vegetation pools
    ideadstemxf_to_iout_fin           => cnveg_nitrogenflux_inst%ideadstemxf_to_iout_fi         , & 
                          ! Input: [integer (:)] Index of fire related N transfer from dead stem transfer pool to outside of vegetation pools
    ilivecroot_to_iout_fin            => cnveg_nitrogenflux_inst%ilivecroot_to_iout_fi          , & 
                          ! Input: [integer (:)] Index of fire related N transfer from live coarse root pool to outside of vegetation pools
    ilivecrootst_to_iout_fin          => cnveg_nitrogenflux_inst%ilivecrootst_to_iout_fi        , & 
                          ! Input: [integer (:)] Index of fire related N transfer from live coarse root storage pool to outside of vegetation pools
    ilivecrootxf_to_iout_fin          => cnveg_nitrogenflux_inst%ilivecrootxf_to_iout_fi        , & 
                          ! Input: [integer (:)] Index of fire related N transfer from live coarse root transfer pool to outside of vegetation pools
    ideadcroot_to_iout_fin            => cnveg_nitrogenflux_inst%ideadcroot_to_iout_fi          , & 
                          ! Input: [integer (:)] Index of fire related N transfer from dead coarse root pool to outside of vegetation pools
    ideadcrootst_to_iout_fin          => cnveg_nitrogenflux_inst%ideadcrootst_to_iout_fi        , & 
                          ! Input: [integer (:)] Index of fire related N transfer from dead coarse root storage pool to outside of vegetation pools
    ideadcrootxf_to_iout_fin          => cnveg_nitrogenflux_inst%ideadcrootxf_to_iout_fi        , & 
                          ! Input: [integer (:)] Index of fire related N transfer from dead coarse root transfer pool to outside of vegetation pools
    ilivestem_to_ideadstem_fin        => cnveg_nitrogenflux_inst%ilivestem_to_ideadstem_fi      , & 
                          ! Input: [integer (:)] Index of fire related N transfer from live stem to dead stem pool
    ilivecroot_to_ideadcroot_fin      => cnveg_nitrogenflux_inst%ilivecroot_to_ideadcroot_fi    , & 
                          ! Input: [integer (:)] Index of fire related N transfer from live coarse root pool to dead coarse root pool
    iretransn_to_iout_fin             => cnveg_nitrogenflux_inst%iretransn_to_iout_fi              &
                          ! Input: [integer (:)] Index of fire related N transfer from retranslocated N pool to outside of vegetation pools
    )
td: associate(                          &
 
  ! Sparse matrix type of A*K 
    AKphvegc                            => cnveg_carbonflux_inst%AKphvegc           , & ! Aph*Kph for C cycle in sparse matrix format
    AKgmvegc                            => cnveg_carbonflux_inst%AKgmvegc           , & ! Agm*Kgm for C cycle in sparse matrix format
    AKfivegc                            => cnveg_carbonflux_inst%AKfivegc           , & ! Afi*Kfi for C cycle in sparse matrix format
    AKallvegc                           => cnveg_carbonflux_inst%AKallvegc          , & ! Aph*Kph + Agm*Kgm + Afi*Kfi for C cycle in sparse matrix format
    NE_AKallvegc                        => cnveg_carbonflux_inst%NE_AKallvegc       , & ! Number of entries in AKallvegc
    RI_AKallvegc                        => cnveg_carbonflux_inst%RI_AKallvegc       , & ! Row indices in Akallvegc
    CI_AKallvegc                        => cnveg_carbonflux_inst%CI_AKallvegc       , & ! Column indices in AKallvegc
    Kvegc                               => cnveg_carbonflux_inst%Kvegc              , & ! Temporary variable of Kph, Kgm or Kfi for C cycle in diagonal matrix format
    Xvegc                               => cnveg_carbonflux_inst%Xvegc              , & ! Vegetation C of each compartment in a vector format
    AKphvegn                            => cnveg_nitrogenflux_inst%AKphvegn         , & ! Aph*Kph for N cycle in sparse matrix format
    AKgmvegn                            => cnveg_nitrogenflux_inst%AKgmvegn         , & ! Agm*Kgm for N cycle in sparse matrix format
    AKfivegn                            => cnveg_nitrogenflux_inst%AKfivegn         , & ! Afi*Kfi for N cycle in sparse matrix format
    AKallvegn                           => cnveg_nitrogenflux_inst%AKallvegn        , & ! Aph*Kph + Agm*Kgm + Afi*Kfi for N cycle in sparse matrix format
    NE_AKallvegn                        => cnveg_nitrogenflux_inst%NE_AKallvegn     , & ! Number of entries in AKallvegn
    RI_AKallvegn                        => cnveg_nitrogenflux_inst%RI_AKallvegn     , & ! Row indices in Akallvegn
    CI_AKallvegn                        => cnveg_nitrogenflux_inst%CI_AKallvegn     , & ! Column indices in AKallvegn
    Kvegn                               => cnveg_nitrogenflux_inst%Kvegn            , & ! Temporary variable of Kph, Kgm or Kfi for N cycle in diagonal matrix format
    Xvegn                               => cnveg_nitrogenflux_inst%Xvegn            , & ! Vegetation N of each compartment in a vector format
    Xveg13c                             => cnveg_carbonflux_inst%Xveg13c            , & ! Vegetation C13 of each compartment in a vector format
    Xveg14c                             => cnveg_carbonflux_inst%Xveg14c            , & ! Vegetation C14 of each compartment in a vector format

  ! Row and column indices of A matrices
    RI_phc                              => cnveg_carbonflux_inst%RI_phc             , & ! Row indices of non-diagonal entires in Aph for C cycle
    CI_phc                              => cnveg_carbonflux_inst%CI_phc             , & ! Column indices of non-diagonal entries in Aph for C cycle
    RI_gmc                              => cnveg_carbonflux_inst%RI_gmc             , & ! Row indices of non-diagonal entires in Agm for C cycle
    CI_gmc                              => cnveg_carbonflux_inst%CI_gmc             , & ! Column indices of non-diagonal entries in Agm for C cycle
    RI_fic                              => cnveg_carbonflux_inst%RI_fic             , & ! Row indices of non-diagonal entires in Afi for C cycle
    CI_fic                              => cnveg_carbonflux_inst%CI_fic             , & ! Column indices of non-diagonal entries in Afi for C cycle
    RI_phn                              => cnveg_nitrogenflux_inst%RI_phn           , & ! Row indices of non-diagonal entires in Aph for N cycle
    CI_phn                              => cnveg_nitrogenflux_inst%CI_phn           , & ! Column indices of non-diagonal entries in Aph for N cycle
    RI_gmn                              => cnveg_nitrogenflux_inst%RI_gmn           , & ! Row indices of non-diagonal entires in Agm for N cycle
    CI_gmn                              => cnveg_nitrogenflux_inst%CI_gmn           , & ! Column indices of non-diagonal entries in Agm for N cycle
    RI_fin                              => cnveg_nitrogenflux_inst%RI_fin           , & ! Row indices of non-diagonal entires in Afi for N cycle
    CI_fin                              => cnveg_nitrogenflux_inst%CI_fin           , & ! Column indices of non-diagonal entries in Afi for N cycle

  ! Following list contains indices of non-diagonal entries in full sparse matrix
    list_aphc                           => cnveg_carbonflux_inst%list_aphc          , & ! Indices of non-diagnoal entries in full sparse matrix Aph for C cycle
    list_agmc                           => cnveg_carbonflux_inst%list_agmc          , & ! Indices of non-diagnoal entries in full sparse matrix Agm for C cycle
    list_afic                           => cnveg_carbonflux_inst%list_afic          , & ! Indices of non-diagnoal entries in full sparse matrix Afi for C cycle
    list_aphn                           => cnveg_nitrogenflux_inst%list_aphn        , & ! Indices of non-diagnoal entries in full sparse matrix Aph for N cycle
    list_agmn                           => cnveg_nitrogenflux_inst%list_agmn        , & ! Indices of non-diagnoal entries in full sparse matrix Agm for N cycle
    list_afin                           => cnveg_nitrogenflux_inst%list_afin        , & ! Indices of non-diagnoal entries in full sparse matrix Afi for N cycle

  ! For sparse matrix A, B and A + B, following list contains locations of entries in A or B or C mapped into matrix (A+B) or (A+B+C)
    list_phc_phgm                       => cnveg_carbonflux_inst%list_phc_phgmc     , & ! The locations of entries in AKphvegc mapped into (AKphvegc+AKgmvegc)
    list_gmc_phgm                       => cnveg_carbonflux_inst%list_gmc_phgmc     , & ! The locations of entries in AKgmvegc mapped into (AKphvegc+AKgmvegc)
    list_phc_phgmfi                     => cnveg_carbonflux_inst%list_phc_phgmfic   , & ! The locations of entries in AKphvegc mapped into (AKphvegc+AKgmvegc+AKfivegc)
    list_gmc_phgmfi                     => cnveg_carbonflux_inst%list_gmc_phgmfic   , & ! The locations of entries in AKgmvegc mapped into (AKphvegc+AKgmvegc+AKfivegc)
    list_fic_phgmfi                     => cnveg_carbonflux_inst%list_fic_phgmfic   , & ! The locations of entries in AKfivegc mapped into (AKphvegc+AKgmvegc+AKfivegc)
    list_phn_phgm                       => cnveg_nitrogenflux_inst%list_phn_phgmn   , & ! The locations of entries in AKphvegn mapped into (AKphvegn+AKgmvegn)
    list_gmn_phgm                       => cnveg_nitrogenflux_inst%list_gmn_phgmn   , & ! The locations of entries in AKgmvegn mapped into (AKphvegn+AKgmvegn)
    list_phn_phgmfi                     => cnveg_nitrogenflux_inst%list_phn_phgmfin , & ! The locations of entries in AKphvegn mapped into (AKphvegn+AKgmvegn+AKfivegn)
    list_gmn_phgmfi                     => cnveg_nitrogenflux_inst%list_gmn_phgmfin , & ! The locations of entries in AKgmvegn mapped into (AKphvegn+AKgmvegn+AKfivegn)
    list_fin_phgmfi                     => cnveg_nitrogenflux_inst%list_fin_phgmfin   & ! The locations of entries in AKfivegn mapped into (AKphvegn+AKgmvegn+AKfivegn)
    )
#ifdef _OPENMP
     nthreads = OMP_GET_MAX_THREADS()
#endif
   !-----------------------------------------------------------------------
    ! set time steps
      call t_startf('CN veg matrix-init')
      dt = real( get_step_size(), r8 )

    ! Initialize local variables
      call vegmatrixc_input%InitV(nvegcpool,bounds%begp,bounds%endp)
      if(use_c13)then
         call vegmatrixc13_input%InitV(nvegcpool,bounds%begp,bounds%endp)
      end if
      if(use_c14)then
         call vegmatrixc14_input%InitV(nvegcpool,bounds%begp,bounds%endp)
      end if
      call vegmatrixn_input%InitV(nvegnpool,bounds%begp,bounds%endp)
      
      matrix_calloc_acc    (:)     = 0._r8
      matrix_nalloc_acc    (:)     = 0._r8
      matrix_ctransfer_acc (:,:)   = 0._r8
      matrix_ntransfer_acc (:,:)   = 0._r8
      if(use_c13)then
         matrix_c13alloc_acc    (:)     = 0._r8
         matrix_c13transfer_acc (:,:)   = 0._r8
      end if
      if(use_c14)then
         matrix_c14alloc_acc    (:)     = 0._r8
         matrix_c14transfer_acc (:,:)   = 0._r8
      end if

      AKinvc (:,:) = 0._r8
      AKinvn (:,:) = 0._r8
      
      epsi = 1.e-30_r8     ! small number   
      
      call t_stopf('CN veg matrix-init')

      call t_startf('CN veg matrix-assigning matrix')

  ! Calculate A matrices from C transfers and C turnovers      
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
               if(use_c14)then
                  associate( &
                        matrix_c14fitransfer  => c14_cnveg_carbonflux_inst%matrix_fitransfer_patch    , & ! Input:  [real(r8) (:,:)] (gC/m2/s) C transfer rate from fire processes, updated in (CNFireBaseMod or CNFireLi2014Mod) and CNC14decayMod
                        matrix_c14fiturnover  => c14_cnveg_carbonflux_inst%matrix_fiturnover_patch      & ! Output: [real(r8) (:,:)] (gC/m2/step) C turnover rate from fire processe, updated in CNVegMatrixMods
                  )
                  if(matrix_c14fiturnover(p,doner_fic(k)) .ne. 0)then
                     Afic14oned(p,k) = matrix_c14fitransfer(p,k) * dt / matrix_c14fiturnover(p,doner_fic(k))
                  else
                     Afic14oned(p,k) = 0._r8
                  end if
                  end associate
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

      call t_stopf('CN veg matrix-assigning matrix')

  ! Assign old state variables to vector Xveg* 
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
      end do

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         if(ivt(p) >= npcropmin)then
            ! Use one index of the grain reproductive pools to operate on
            Xvegc%V(p,igrain)     = reproductivec(p,irepr)
            Xvegc%V(p,igrain_st)  = reproductivec_storage(p,irepr)
            Xvegc%V(p,igrain_xf)  = reproductivec_xfer(p,irepr)
         end if
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
               ! Use one index of the grain reproductive pools to operate on
               Xveg13c%V(p,igrain)     = cs13_veg%reproductivec_patch(p,irepr)
               Xveg13c%V(p,igrain_st)  = cs13_veg%reproductivec_storage_patch(p,irepr)
               Xveg13c%V(p,igrain_xf)  = cs13_veg%reproductivec_xfer_patch(p,irepr)
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
               ! Use one index of the grain reproductive pools to operate on
               Xveg14c%V(p,igrain)     = cs14_veg%reproductivec_patch(p,irepr)
               Xveg14c%V(p,igrain_st)  = cs14_veg%reproductivec_storage_patch(p,irepr)
               Xveg14c%V(p,igrain_xf)  = cs14_veg%reproductivec_xfer_patch(p,irepr)
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
         if(ivt(p) >= npcropmin)then
            Xvegn%V(p,igrain)        = sum(reproductiven(p,:))
            Xvegn%V(p,igrain_st)     = sum(reproductiven_storage(p,:))
            Xvegn%V(p,igrain_xf)     = sum(reproductiven_xfer(p,:))
         end if
      end do

  ! Save *c0* and *n0* variables at begin of each year.
      if (is_beg_curr_year())then
         iyr = iyr + 1
         if(mod(iyr-1,nyr_forcing) .eq. 0)then
            iloop = iloop + 1
         end if
         if(.not. spinup_matrixcn .or. spinup_matrixcn .and. mod(iyr-1,nyr_SASU) .eq. 0)then
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
                  ! Use one index of the grain reproductive pools to operate on
                  reproc0(p)               = max(reproductivec(p,irepr),  epsi)
                  reproc0_storage(p)       = max(reproductivec_storage(p,irepr),  epsi)
                  reproc0_xfer(p)          = max(reproductivec_xfer(p,irepr),  epsi)
               end if
            end do

            if(use_c13)then
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  cs13_veg%leafc0_patch(p)                = max(cs13_veg%leafc_patch(p),  epsi)
                  cs13_veg%leafc0_storage_patch(p)        = max(cs13_veg%leafc_storage_patch(p),  epsi)
                  cs13_veg%leafc0_xfer_patch(p)           = max(cs13_veg%leafc_xfer_patch(p),  epsi)
                  cs13_veg%frootc0_patch(p)               = max(cs13_veg%frootc_patch(p),  epsi)
                  cs13_veg%frootc0_storage_patch(p)       = max(cs13_veg%frootc_storage_patch(p),  epsi)
                  cs13_veg%frootc0_xfer_patch(p)          = max(cs13_veg%frootc_xfer_patch(p),  epsi)
                  cs13_veg%livestemc0_patch(p)            = max(cs13_veg%livestemc_patch(p),  epsi)
                  cs13_veg%livestemc0_storage_patch(p)    = max(cs13_veg%livestemc_storage_patch(p),  epsi)
                  cs13_veg%livestemc0_xfer_patch(p)       = max(cs13_veg%livestemc_xfer_patch(p),  epsi)
                  cs13_veg%deadstemc0_patch(p)            = max(cs13_veg%deadstemc_patch(p),  epsi)
                  cs13_veg%deadstemc0_storage_patch(p)    = max(cs13_veg%deadstemc_storage_patch(p),  epsi)
                  cs13_veg%deadstemc0_xfer_patch(p)       = max(cs13_veg%deadstemc_xfer_patch(p),  epsi)
                  cs13_veg%livecrootc0_patch(p)           = max(cs13_veg%livecrootc_patch(p),  epsi)
                  cs13_veg%livecrootc0_storage_patch(p)   = max(cs13_veg%livecrootc_storage_patch(p),  epsi)
                  cs13_veg%livecrootc0_xfer_patch(p)      = max(cs13_veg%livecrootc_xfer_patch(p),  epsi)
                  cs13_veg%deadcrootc0_patch(p)           = max(cs13_veg%deadcrootc_patch(p),  epsi)
                  cs13_veg%deadcrootc0_storage_patch(p)   = max(cs13_veg%deadcrootc_storage_patch(p),  epsi)
                  cs13_veg%deadcrootc0_xfer_patch(p)      = max(cs13_veg%deadcrootc_xfer_patch(p),  epsi)
               end do

               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  if(ivt(p) >= npcropmin)then
                     ! Use one index of the grain reproductive pools to operate on
                     cs13_veg%reproc0_patch(p)               = max(cs13_veg%reproductivec_patch(p,irepr),  epsi)
                     cs13_veg%reproc0_storage_patch(p)       = max(cs13_veg%reproductivec_storage_patch(p,irepr), epsi)
                     cs13_veg%reproc0_xfer_patch(p)          = max(cs13_veg%reproductivec_xfer_patch(p,irepr), epsi)
                  end if
               end do
            end if

            if(use_c14)then
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  cs14_veg%leafc0_patch(p)                = max(cs14_veg%leafc_patch(p),  epsi)
                  cs14_veg%leafc0_storage_patch(p)        = max(cs14_veg%leafc_storage_patch(p),  epsi)
                  cs14_veg%leafc0_xfer_patch(p)           = max(cs14_veg%leafc_xfer_patch(p),  epsi)
                  cs14_veg%frootc0_patch(p)               = max(cs14_veg%frootc_patch(p),  epsi)
                  cs14_veg%frootc0_storage_patch(p)       = max(cs14_veg%frootc_storage_patch(p),  epsi)
                  cs14_veg%frootc0_xfer_patch(p)          = max(cs14_veg%frootc_xfer_patch(p),  epsi)
                  cs14_veg%livestemc0_patch(p)            = max(cs14_veg%livestemc_patch(p),  epsi)
                  cs14_veg%livestemc0_storage_patch(p)    = max(cs14_veg%livestemc_storage_patch(p),  epsi)
                  cs14_veg%livestemc0_xfer_patch(p)       = max(cs14_veg%livestemc_xfer_patch(p),  epsi)
                  cs14_veg%deadstemc0_patch(p)            = max(cs14_veg%deadstemc_patch(p),  epsi)
                  cs14_veg%deadstemc0_storage_patch(p)    = max(cs14_veg%deadstemc_storage_patch(p),  epsi)
                  cs14_veg%deadstemc0_xfer_patch(p)       = max(cs14_veg%deadstemc_xfer_patch(p),  epsi)
                  cs14_veg%livecrootc0_patch(p)           = max(cs14_veg%livecrootc_patch(p),  epsi)
                  cs14_veg%livecrootc0_storage_patch(p)   = max(cs14_veg%livecrootc_storage_patch(p),  epsi)
                  cs14_veg%livecrootc0_xfer_patch(p)      = max(cs14_veg%livecrootc_xfer_patch(p),  epsi)
                  cs14_veg%deadcrootc0_patch(p)           = max(cs14_veg%deadcrootc_patch(p),  epsi)
                  cs14_veg%deadcrootc0_storage_patch(p)   = max(cs14_veg%deadcrootc_storage_patch(p),  epsi)
                  cs14_veg%deadcrootc0_xfer_patch(p)      = max(cs14_veg%deadcrootc_xfer_patch(p),  epsi)
               end do

               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  if(ivt(p) >= npcropmin)then
                     ! Use one index of the grain reproductive pools to operate on
                     cs14_veg%reproc0_patch(p)               = max(cs14_veg%reproductivec_patch(p,irepr),  epsi)
                     cs14_veg%reproc0_storage_patch(p)       = max(cs14_veg%reproductivec_storage_patch(p,irepr),  epsi)
                     cs14_veg%reproc0_xfer_patch(p)          = max(cs14_veg%reproductivec_xfer_patch(p,irepr),  epsi)
                  end if
               end do
            end if

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
                  ! Use one index of the grain reproductive pools to operate on
                  repron0(p)               = max(reproductiven(p,irepr),  epsi)
                  repron0_storage(p)       = max(reproductiven_storage(p,irepr),  epsi)
                  repron0_xfer(p)          = max(reproductiven_xfer(p,irepr),  epsi)
               end if
            end do
         end if
      end if    

         call t_stopf('CN veg matrix-set old value')

         call t_startf('CN veg matrix-matrix multi.')

  ! Start matrix operation
  ! Calculate B*I

         do i=1,nvegcpool
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               vegmatrixc_input%V(p,i) = matrix_alloc(p,i) * matrix_Cinput(p) * dt
            end do
         end do

  ! Set up sparse matrix Aph_c from non-diagonal entires Aphconed, diagonal entries are all set to -1.
  ! Note that AKphvegc here only represent A matrix instead of A * K

         if(ncphtrans .gt. ncphouttrans)then
            AI_phc = receiver_phc(1:ncphtrans-ncphouttrans)
            AJ_phc = doner_phc   (1:ncphtrans-ncphouttrans)
            call AKphvegc%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Aphconed,&
                 AI_phc,AJ_phc,ncphtrans-ncphouttrans,init_ready_aphc,list_aphc,RI_phc,CI_phc)
         else
            call AKphvegc%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if

  ! Set up diagonal matrix Kph_c from diagonal entries matrix_phturnover
         call Kvegc%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_phturnover(bounds%begp:bounds%endp,1:nvegcpool))

  ! Calculate Aph_c*Kph_c using SPMM_AK.
         call AKphvegc%SPMM_AK(num_soilp,filter_soilp,Kvegc)



  ! Set up sparse matrix Agm_c from non-diagonal entires Agmconed, diagonal entries are all set to -1.
  ! Note that AKgmvegc here only represent A matrix instead of A * K

         if(ncgmtrans .gt. ncgmouttrans)then
            AI_gmc = receiver_gmc(1:ncgmtrans-ncgmouttrans)
            AJ_gmc = doner_gmc   (1:ncgmtrans-ncgmouttrans)
            call AKgmvegc%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Agmconed,&
                 AI_gmc,AJ_gmc,ncgmtrans-ncgmouttrans,init_ready_agmc,list_agmc,RI_gmc,CI_gmc)
         else
            call AKgmvegc%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if

  ! Set up diagonal matrix Kgm_c from diagonal entries matrix_gmturnover
         call Kvegc%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_gmturnover(bounds%begp:bounds%endp,1:nvegcpool))

  ! Calculate Agm_c*Kgm_c using SPMM_AK.
         call AKgmvegc%SPMM_AK(num_soilp,filter_soilp,Kvegc)



  ! Set up sparse matrix Afi_c from non-diagonal entires Aficoned, diagonal entries are all set to -1.
  ! Note that AKfivegc here only represent A matrix instead of A * K

         if(ncfitrans .gt. ncfiouttrans)then
            AI_fic = receiver_fic(1:ncfitrans-ncfiouttrans)
            AJ_fic = doner_fic   (1:ncfitrans-ncfiouttrans)
            call AKfivegc%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Aficoned,&
                 AI_fic,AJ_fic,ncfitrans-ncfiouttrans,init_ready_afic,list_afic,RI_fic,CI_fic)
            if(use_c14)then
               associate( &
                   AKfivegc14                          => c14_cnveg_carbonflux_inst%AKfivegc     , & ! Afi*Kfi for C14 cycle in sparse matrix format
                   RI_fic14                            => c14_cnveg_carbonflux_inst%RI_fic             , & ! Row indices of non-diagonal entires in Afi for C cycle
                   CI_fic14                            => c14_cnveg_carbonflux_inst%CI_fic             , & ! Column indices of non-diagonal entries in Afi for C cycle
                   list_afic14                         => c14_cnveg_carbonflux_inst%list_afic            & ! Indices of non-diagnoal entries in full sparse matrix Afi for C cycle
               )
               AI_fic14 = receiver_fic(1:ncfitrans-ncfiouttrans)
               AJ_fic14 = doner_fic   (1:ncfitrans-ncfiouttrans)
               call AKfivegc14%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Afic14oned,&
                    AI_fic14,AJ_fic14,ncfitrans-ncfiouttrans,init_ready_afic14,list_afic14,RI_fic14,CI_fic14)
               end associate
            end if
         else
            call AKfivegc%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
            if(use_c14)then
               associate( &
                   AKfivegc14 => c14_cnveg_carbonflux_inst%AKfivegc  & ! Afi*Kfi for C14 cycle in sparse matrix format
               )
               call AKfivegc14%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
               end associate
            end if
         end if

  ! Set up diagonal matrix Kfi_c from diagonal entries matrix_fiturnover
         call Kvegc%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_fiturnover(bounds%begp:bounds%endp,1:nvegcpool))

  ! Calculate Afi_c*Kfi_c using SPMM_AK.
         call AKfivegc%SPMM_AK(num_soilp,filter_soilp,Kvegc)

         if(use_c14)then
            associate( &
                AKfivegc14                          => c14_cnveg_carbonflux_inst%AKfivegc     , & ! Afi*Kfi for C14 cycle in sparse matrix format
                matrix_c14fitransfer  => c14_cnveg_carbonflux_inst%matrix_fitransfer_patch    , & ! Input:  [real(r8) (:,:)] (gC/m2/s) C transfer rate from fire processes, updated in (CNFireBaseMod or CNFireLi2014Mod) and CNC14decayMod
                matrix_c14fiturnover  => c14_cnveg_carbonflux_inst%matrix_fiturnover_patch      & ! Output: [real(r8) (:,:)] (gC/m2/step) C turnover rate from fire processe, updated in CNVegMatrixMods
            )
  ! Set up diagonal matrix Kfi_c from diagonal entries matrix_fiturnover
            call Kvegc%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_c14fiturnover(bounds%begp:bounds%endp,1:nvegcpool))

  ! Calculate Afi_c*Kfi_c using SPMM_AK.
            call AKfivegc14%SPMM_AK(num_soilp,filter_soilp,Kvegc)
            end associate
         end if

  ! Caclulate AKallvegc = Aph_c*Kph_c + Agm_c*Kgm_c + Afi_c*Kfi_c
  ! When no fire, Afi_c*Kfi_c = 0, AKallvegc = Aph_c*Kph_c + Agm_c*Kgm_c
  ! When fire is on, AKallvegc = Aph_c*Kph_c + Agm_c*Kgm_c + Afi_c*Kfi_c

         if(num_actfirep .eq. 0 .and. nthreads < 2)then
            call AKallvegc%SPMP_AB(num_soilp,filter_soilp,AKphvegc,AKgmvegc,list_ready_phgmc,list_A=list_phc_phgm,list_B=list_gmc_phgm,&
                 NE_AB=NE_AKallvegc,RI_AB=RI_AKallvegc,CI_AB=CI_AKallvegc)
         else
            call AKallvegc%SPMP_ABC(num_soilp,filter_soilp,AKphvegc,AKgmvegc,AKfivegc,list_ready_phgmfic,list_A=list_phc_phgmfi,&
                 list_B=list_gmc_phgmfi,list_C=list_fic_phgmfi,NE_ABC=NE_AKallvegc,RI_ABC=RI_AKallvegc,CI_ABC=CI_AKallvegc,&
                 use_actunit_list_C=.True.,num_actunit_C=num_actfirep,filter_actunit_C=filter_actfirep)
         end if
          
         if(use_c14)then
            associate( &
                AKfivegc14                          => c14_cnveg_carbonflux_inst%AKfivegc     , & ! Afi*Kfi for C14 cycle in sparse matrix format
                AKallvegc14                         => c14_cnveg_carbonflux_inst%AKallvegc    , & ! Aph*Kph + Agm*Kgm + Afi*Kfi for C14 cycle in sparse matrix format
                NE_AKallvegc14                      => c14_cnveg_carbonflux_inst%NE_AKallvegc       , & ! Number of entries in AKallvegc
                RI_AKallvegc14                      => c14_cnveg_carbonflux_inst%RI_AKallvegc       , & ! Row indices in Akallvegc
                CI_AKallvegc14                      => c14_cnveg_carbonflux_inst%CI_AKallvegc       , & ! Column indices in AKallvegc
                list_phc14_phgmfi                   => c14_cnveg_carbonflux_inst%list_phc_phgmfic   , & ! The locations of entries in AKphvegc mapped into (AKphvegc+AKgmvegc+AKfivegc)
                list_gmc14_phgmfi                   => c14_cnveg_carbonflux_inst%list_gmc_phgmfic   , & ! The locations of entries in AKgmvegc mapped into (AKphvegc+AKgmvegc+AKfivegc)
                list_fic14_phgmfi                   => c14_cnveg_carbonflux_inst%list_fic_phgmfic     & ! The locations of entries in AKfivegc mapped into (AKphvegc+AKgmvegc+AKfivegc)
            )
            call AKallvegc14%SPMP_ABC(num_soilp,filter_soilp,AKphvegc,AKgmvegc,AKfivegc14,list_ready_phgmfic14,list_A=list_phc14_phgmfi,&
                 list_B=list_gmc14_phgmfi,list_C=list_fic14_phgmfi,NE_ABC=NE_AKallvegc14,RI_ABC=RI_AKallvegc14,CI_ABC=CI_AKallvegc14)
            end associate
         end if


  ! Xvegc_n+1 = (Aph_c*Kph_c + Agm_c*Kgm_c + Afi_c*Kfi_c) * Xvegc_n + Xvegc_n
         call Xvegc%SPMM_AX(num_soilp,filter_soilp,AKallvegc)

  ! Xvegc_n+1 = (Aph_c*Kph_c + Agm_c*Kgm_c + Afi_c*Kfi_c) * Xvegc_n + Xvegc_n + B*I
         do i = 1,nvegcpool
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               Xvegc%V(p,i) = Xvegc%V(p,i) + vegmatrixc_input%V(p,i) 
            end do
         end do
         

         if ( use_c13 ) then
  ! Calculate B*I_C13
            do i=1,nvegcpool
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  vegmatrixc13_input%V(p,i) = matrix_alloc(p,i) * matrix_C13input(p) * dt
               end do
            end do

  ! Xveg13c_n+1 = (Aph_c*Kph_c + Agm_c*Kgm_c + Afi_c*Kfi_c) * Xveg13c_n + Xveg13c_n
            call Xveg13c%SPMM_AX(num_soilp,filter_soilp,AKallvegc)

  ! Xveg13c_n+1 = (Aph_c*Kph_c + Agm_c*Kgm_c + Afi_c*Kfi_c) * Xveg13c_n + Xveg13c_n + B*I_C13
            do i=1,nvegcpool
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  Xveg13c%V(p,i) = Xveg13c%V(p,i) + vegmatrixc13_input%V(p,i) 
               end do
            end do
         end if


         if ( use_c14 ) then
             associate( &
                matrix_C14input => cnveg_carbonflux_inst%matrix_C14input_patch, & ! Input:  [real(r8) (:)] (gC/m2/s) C14 input to vegetation, updated in NutrientCompetitionFlexibleCNMod or NutrientCompetitionCLM45defaultMod
                AKallvegc14     => c14_cnveg_carbonflux_inst%AKallvegc          & ! Aph*Kph + Agm*Kgm + Afi*Kfi for C14 cycle in sparse matrix format
            ) 
  ! Calculate B*I_C14
            do i=1,nvegcpool
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  vegmatrixc14_input%V(p,i) = matrix_alloc(p,i) * matrix_C14input(p) * dt
               end do
            end do

  ! Xveg14c_n+1 = (Aph_c*Kph_c + Agm_c*Kgm_c + Afi_c*Kfi_c) * Xveg14c_n + Xveg14c_n
            call Xveg14c%SPMM_AX(num_soilp,filter_soilp,AKallvegc14)

  ! Xveg14c_n+1 = (Aph_c*Kph_c + Agm_c*Kgm_c + Afi_c*Kfi_c) * Xveg14c_n + Xveg14c_n + B*I_C14
            do i=1,nvegcpool
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  Xveg14c%V(p,i) = Xveg14c%V(p,i) + vegmatrixc14_input%V(p,i) 
               end do
            end do
            end associate
         end if
                            


  ! Calculate B_N*I_N
         do i=1,nvegnpool
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               vegmatrixn_input%V(p,i) = matrix_nalloc(p,i) * matrix_Ninput(p) * dt
            end do
         end do


  ! Set up sparse matrix Aph_n from non-diagonal entires Aficoned, diagonal entries are all set to -1.
  ! Note that AKphvegn here only represent A matrix instead of A * K

         if(nnphtrans .gt. nnphouttrans)then
             AI_phn = receiver_phn(1:nnphtrans-nnphouttrans)
             AJ_phn = doner_phn   (1:nnphtrans-nnphouttrans)
            call AKphvegn%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Aphnoned,&
                 AI_phn,AJ_phn,nnphtrans-nnphouttrans,init_ready_aphn,list_aphn,RI_phn,CI_phn)
         else
            call AKphvegn%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if

  ! Set up diagonal matrix Kph_n from diagonal entries matrix_nphturnover
         call Kvegn%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_nphturnover(bounds%begp:bounds%endp,1:nvegnpool))

  ! Calculate Aph_n*Kph_n using SPMM_AK.
         call AKphvegn%SPMM_AK(num_soilp,filter_soilp,Kvegn)


  ! Set up sparse matrix Agm_n from non-diagonal entires Aficoned, diagonal entries are all set to -1.
  ! Note that AKgmvegn here only represent A matrix instead of A * K

         if(nngmtrans .gt. nngmouttrans)then
             AI_gmn = receiver_gmn(1:nngmtrans-nngmouttrans)
             AJ_gmn = doner_gmn   (1:nngmtrans-nngmouttrans)
            call AKgmvegn%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Agmnoned,&
                 AI_gmn,AJ_gmn,nngmtrans-nngmouttrans,init_ready_agmn,list_agmn,RI_gmn,CI_gmn)
         else
            call AKgmvegn%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if

  ! Set up diagonal matrix Kgm_n from diagonal entries matrix_ngmturnover
         call Kvegn%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_ngmturnover(bounds%begp:bounds%endp,1:nvegnpool))

  ! Calculate Agm_n*Kgm_n using SPMM_AK.
         call AKgmvegn%SPMM_AK(num_soilp,filter_soilp,Kvegn)


  ! Set up sparse matrix Afi_n from non-diagonal entires Aficoned, diagonal entries are all set to -1.
  ! Note that AKfivegn here only represent A matrix instead of A * K

         if(nnfitrans .gt. nnfiouttrans)then
             AI_fin = receiver_fin(1:nnfitrans-nnfiouttrans)
             AJ_fin = doner_fin   (1:nnfitrans-nnfiouttrans)
            call AKfivegn%SetValueA(bounds%begp,bounds%endp,num_soilp,filter_soilp,Afinoned,&
                 AI_fin,AJ_fin,nnfitrans-nnfiouttrans,init_ready_afin,list_afin,RI_fin,CI_fin)
         else
            call AKfivegn%SetValueA_diag(num_soilp,filter_soilp,-1._r8)
         end if

  ! Set up diagonal matrix Kfi_n from diagonal entries matrix_nfiturnover
         call Kvegn%SetValueDM(bounds%begp,bounds%endp,num_soilp,filter_soilp,matrix_nfiturnover(bounds%begp:bounds%endp,1:nvegnpool))

  ! Calculate Afi_n*Kfi_n using SPMM_AK.
         call AKfivegn%SPMM_AK(num_soilp,filter_soilp,Kvegn)

         
  ! Caclulate AKallvegn = Aph_n*Kph_n + Agm_n*Kgm_n + Afi_n*Kfi_n
  ! When no fire, Afi_n*Kfi_n = 0, AKallvegn = Aph_n*Kph_n + Agm_n*Kgm_n
  ! When fire is on, AKallvegn = Aph_n*Kph_n + Agm_n*Kgm_n + Afi_n*Kfi_n

         if(num_actfirep .eq. 0 .and. nthreads < 2)then
            call AKallvegn%SPMP_AB(num_soilp,filter_soilp,AKphvegn,AKgmvegn,list_ready_phgmn,list_A=list_phn_phgm,list_B=list_gmn_phgm,&
                 NE_AB=NE_AKallvegn,RI_AB=RI_AKallvegn,CI_AB=CI_AKallvegn)
         else
            call AKallvegn%SPMP_ABC(num_soilp,filter_soilp,AKphvegn,AKgmvegn,AKfivegn,list_ready_phgmfin,list_A=list_phn_phgmfi,&
                 list_B=list_gmn_phgmfi,list_C=list_fin_phgmfi,NE_ABC=NE_AKallvegn,RI_ABC=RI_AKallvegn,CI_ABC=CI_AKallvegn,&
                 use_actunit_list_C=.True.,num_actunit_C=num_actfirep,filter_actunit_C=filter_actfirep)
         end if
         
  ! Xvegn_n+1 = (Aph_n*Kph_n + Agm_n*Kgm_n + Afi_n*Kfi_n) * Xvegc_n + Xvegc_n
         call Xvegn%SPMM_AX(num_soilp,filter_soilp,AKallvegn)
          
  ! Xvegn_n+1 = (Aph_n*Kph_n + Agm_n*Kgm_n + Afi_n*Kfi_n) * Xvegc_n + Xvegc_n + B_N*I_N
         do i=1,nvegnpool
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               Xvegn%V(p,i) = Xvegn%V(p,i) + vegmatrixn_input%V(p,i)
            end do
         end do

         call t_stopf('CN veg matrix-matrix multi.')


  ! Accumulate transfers during the whole calendar year

         call t_startf('CN veg matrix-accum. trans.')
         if(spinup_matrixcn .or. hist_wrt_matrixcn_diag)then
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
               if(use_c13)then
                  cs13_veg%matrix_calloc_leaf_acc_patch(p)        = cs13_veg%matrix_calloc_leaf_acc_patch(p)        + vegmatrixc13_input%V(p,ileaf)
                  cs13_veg%matrix_calloc_leafst_acc_patch(p)      = cs13_veg%matrix_calloc_leafst_acc_patch(p)      + vegmatrixc13_input%V(p,ileaf_st)
                  cs13_veg%matrix_calloc_froot_acc_patch(p)       = cs13_veg%matrix_calloc_froot_acc_patch(p)       + vegmatrixc13_input%V(p,ifroot)
                  cs13_veg%matrix_calloc_frootst_acc_patch(p)     = cs13_veg%matrix_calloc_frootst_acc_patch(p)     + vegmatrixc13_input%V(p,ifroot_st)
                  cs13_veg%matrix_calloc_livestem_acc_patch(p)    = cs13_veg%matrix_calloc_livestem_acc_patch(p)    + vegmatrixc13_input%V(p,ilivestem)
                  cs13_veg%matrix_calloc_livestemst_acc_patch(p)  = cs13_veg%matrix_calloc_livestemst_acc_patch(p)  + vegmatrixc13_input%V(p,ilivestem_st)
                  cs13_veg%matrix_calloc_deadstem_acc_patch(p)    = cs13_veg%matrix_calloc_deadstem_acc_patch(p)    + vegmatrixc13_input%V(p,ideadstem)
                  cs13_veg%matrix_calloc_deadstemst_acc_patch(p)  = cs13_veg%matrix_calloc_deadstemst_acc_patch(p)  + vegmatrixc13_input%V(p,ideadstem_st)
                  cs13_veg%matrix_calloc_livecroot_acc_patch(p)   = cs13_veg%matrix_calloc_livecroot_acc_patch(p)   + vegmatrixc13_input%V(p,ilivecroot)
                  cs13_veg%matrix_calloc_livecrootst_acc_patch(p) = cs13_veg%matrix_calloc_livecrootst_acc_patch(p) + vegmatrixc13_input%V(p,ilivecroot_st)
                  cs13_veg%matrix_calloc_deadcroot_acc_patch(p)   = cs13_veg%matrix_calloc_deadcroot_acc_patch(p)   + vegmatrixc13_input%V(p,ideadcroot)
                  cs13_veg%matrix_calloc_deadcrootst_acc_patch(p) = cs13_veg%matrix_calloc_deadcrootst_acc_patch(p) + vegmatrixc13_input%V(p,ideadcroot_st)
               end if
               if(use_c14)then
                  cs14_veg%matrix_calloc_leaf_acc_patch(p)        = cs14_veg%matrix_calloc_leaf_acc_patch(p)        + vegmatrixc14_input%V(p,ileaf)
                  cs14_veg%matrix_calloc_leafst_acc_patch(p)      = cs14_veg%matrix_calloc_leafst_acc_patch(p)      + vegmatrixc14_input%V(p,ileaf_st)
                  cs14_veg%matrix_calloc_froot_acc_patch(p)       = cs14_veg%matrix_calloc_froot_acc_patch(p)       + vegmatrixc14_input%V(p,ifroot)
                  cs14_veg%matrix_calloc_frootst_acc_patch(p)     = cs14_veg%matrix_calloc_frootst_acc_patch(p)     + vegmatrixc14_input%V(p,ifroot_st)
                  cs14_veg%matrix_calloc_livestem_acc_patch(p)    = cs14_veg%matrix_calloc_livestem_acc_patch(p)    + vegmatrixc14_input%V(p,ilivestem)
                  cs14_veg%matrix_calloc_livestemst_acc_patch(p)  = cs14_veg%matrix_calloc_livestemst_acc_patch(p)  + vegmatrixc14_input%V(p,ilivestem_st)
                  cs14_veg%matrix_calloc_deadstem_acc_patch(p)    = cs14_veg%matrix_calloc_deadstem_acc_patch(p)    + vegmatrixc14_input%V(p,ideadstem)
                  cs14_veg%matrix_calloc_deadstemst_acc_patch(p)  = cs14_veg%matrix_calloc_deadstemst_acc_patch(p)  + vegmatrixc14_input%V(p,ideadstem_st)
                  cs14_veg%matrix_calloc_livecroot_acc_patch(p)   = cs14_veg%matrix_calloc_livecroot_acc_patch(p)   + vegmatrixc14_input%V(p,ilivecroot)
                  cs14_veg%matrix_calloc_livecrootst_acc_patch(p) = cs14_veg%matrix_calloc_livecrootst_acc_patch(p) + vegmatrixc14_input%V(p,ilivecroot_st)
                  cs14_veg%matrix_calloc_deadcroot_acc_patch(p)   = cs14_veg%matrix_calloc_deadcroot_acc_patch(p)   + vegmatrixc14_input%V(p,ideadcroot)
                  cs14_veg%matrix_calloc_deadcrootst_acc_patch(p) = cs14_veg%matrix_calloc_deadcrootst_acc_patch(p) + vegmatrixc14_input%V(p,ideadcroot_st)
               end if
            end do
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  matrix_calloc_grain_acc(p)    = matrix_calloc_grain_acc(p)       + vegmatrixc_input%V(p,igrain)
                  matrix_calloc_grainst_acc(p)  = matrix_calloc_grainst_acc(p)     + vegmatrixc_input%V(p,igrain_st)
                  if(use_c13)then
                     cs13_veg%matrix_calloc_grain_acc_patch(p)    = cs13_veg%matrix_calloc_grain_acc_patch(p)       + vegmatrixc13_input%V(p,igrain)
                     cs13_veg%matrix_calloc_grainst_acc_patch(p)  = cs13_veg%matrix_calloc_grainst_acc_patch(p)     + vegmatrixc13_input%V(p,igrain_st)
                  end if
                  if(use_c14)then
                     cs14_veg%matrix_calloc_grain_acc_patch(p)    = cs14_veg%matrix_calloc_grain_acc_patch(p)       + vegmatrixc14_input%V(p,igrain)
                     cs14_veg%matrix_calloc_grainst_acc_patch(p)  = cs14_veg%matrix_calloc_grainst_acc_patch(p)     + vegmatrixc14_input%V(p,igrain_st)
                  end if
               end if
            end do

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
               if(use_c13)then
                  cs13_veg%matrix_ctransfer_leafst_to_leafxf_acc_patch(p)           = cs13_veg%matrix_ctransfer_leafst_to_leafxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ileafst_to_ileafxf_phc) &
                                                                                    * dt * cs13_veg%leafc_storage_patch(p) !matrix_phturnover(p,ileaf_st)*leafc_storage(p)
                  cs13_veg%matrix_ctransfer_leafxf_to_leaf_acc_patch(p)             = cs13_veg%matrix_ctransfer_leafxf_to_leaf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ileafxf_to_ileaf_phc) &
                                                                                    * dt * cs13_veg%leafc_xfer_patch(p)!matrix_phturnover(p,ileaf_xf)*leafc_xfer(p)
                  cs13_veg%matrix_ctransfer_frootst_to_frootxf_acc_patch(p)         = cs13_veg%matrix_ctransfer_frootst_to_frootxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) &
                                                                                    * dt * cs13_veg%frootc_storage_patch(p)!matrix_phturnover(p,ifroot_st)*frootc_storage(p)
                  cs13_veg%matrix_ctransfer_frootxf_to_froot_acc_patch(p)           = cs13_veg%matrix_ctransfer_frootxf_to_froot_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ifrootxf_to_ifroot_phc) &
                                                                                    * dt * cs13_veg%frootc_xfer_patch(p)!matrix_phturnover(p,ifroot_xf)*frootc_xfer(p)
                  cs13_veg%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p)   = cs13_veg%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) &
                                                                                    * dt * cs13_veg%livestemc_storage_patch(p)!matrix_phturnover(p,ilivestem_st)*livestemc_storage(p)
                  cs13_veg%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p)     = cs13_veg%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc) &
                                                                                    * dt * cs13_veg%livestemc_xfer_patch(p)!matrix_phturnover(p,ilivestem_xf)*livestemc_xfer(p)
                  cs13_veg%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p)   = cs13_veg%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) &
                                                                                    * dt * cs13_veg%deadstemc_storage_patch(p)!matrix_phturnover(p,ideadstem_st)*deadstemc_storage(p)
                  cs13_veg%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p)     = cs13_veg%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc) &
                                                                                    * dt * cs13_veg%deadstemc_xfer_patch(p)!matrix_phturnover(p,ideadstem_xf)*deadstemc_xfer(p)
                  cs13_veg%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p) = cs13_veg%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) &
                                                                                    * dt * cs13_veg%livecrootc_storage_patch(p)!matrix_phturnover(p,ilivecroot_st)*livecrootc_storage(p)
                  cs13_veg%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p)   = cs13_veg%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) &
                                                                                    * dt * cs13_veg%livecrootc_xfer_patch(p)!matrix_phturnover(p,ilivecroot_xf)*livecrootc_xfer(p)
                  cs13_veg%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p) = cs13_veg%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) &
                                                                                    * dt * cs13_veg%deadcrootc_storage_patch(p)!matrix_phturnover(p,ideadcroot_st)*deadcrootc_storage(p)
                  cs13_veg%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p)   = cs13_veg%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) &
                                                                                    * dt * cs13_veg%deadcrootc_xfer_patch(p)!matrix_phturnover(p,ideadcroot_st)*deadcrootc_xfer(p)
                  cs13_veg%matrix_ctransfer_livestem_to_deadstem_acc_patch(p)       = cs13_veg%matrix_ctransfer_livestem_to_deadstem_acc_patch(p) &
                                                                                    +(matrix_phtransfer(p,ilivestem_to_ideadstem_phc)&!matrix_phturnover(p,ilivestem) &
                                                                                    + matrix_fitransfer(p,ilivestem_to_ideadstem_fic))&!matrix_fiturnover(p,ilivestem))&
                                                                                    * dt * cs13_veg%livestemc_patch(p) 
                  cs13_veg%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p)     = cs13_veg%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p) &
                                                                                    +(matrix_phtransfer(p,ilivecroot_to_ideadcroot_phc)&!*matrix_phturnover(p,ilivecroot) &
                                                                                    + matrix_fitransfer(p,ilivecroot_to_ideadcroot_fic))&!*matrix_fiturnover(p,ilivecroot))&
                                                                                    * dt * cs13_veg%livecrootc_patch(p)
                  cs13_veg%matrix_cturnover_leaf_acc_patch(p)        = cs13_veg%matrix_cturnover_leaf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ileaf)+matrix_gmturnover(p,ileaf)+matrix_fiturnover(p,ileaf)) &
                                                                     * cs13_veg%leafc_patch(p)
                  cs13_veg%matrix_cturnover_leafst_acc_patch(p)      = cs13_veg%matrix_cturnover_leafst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ileaf_st)+matrix_gmturnover(p,ileaf_st)+matrix_fiturnover(p,ileaf_st)) &
                                                                     * cs13_veg%leafc_storage_patch(p)
                  cs13_veg%matrix_cturnover_leafxf_acc_patch(p)      = cs13_veg%matrix_cturnover_leafxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ileaf_xf)+matrix_gmturnover(p,ileaf_xf)+matrix_fiturnover(p,ileaf_xf)) &
                                                                     * cs13_veg%leafc_xfer_patch(p)
                  cs13_veg%matrix_cturnover_froot_acc_patch(p)       = cs13_veg%matrix_cturnover_froot_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ifroot)+matrix_gmturnover(p,ifroot)+matrix_fiturnover(p,ifroot)) &
                                                                     * cs13_veg%frootc_patch(p)
                  cs13_veg%matrix_cturnover_frootst_acc_patch(p)     = cs13_veg%matrix_cturnover_frootst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ifroot_st)+matrix_gmturnover(p,ifroot_st)+matrix_fiturnover(p,ifroot_st)) &
                                                                     * cs13_veg%frootc_storage_patch(p)
                  cs13_veg%matrix_cturnover_frootxf_acc_patch(p)     = cs13_veg%matrix_cturnover_frootxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ifroot_xf)+matrix_gmturnover(p,ifroot_xf)+matrix_fiturnover(p,ifroot_xf)) &
                                                                     * cs13_veg%frootc_xfer_patch(p)
                  cs13_veg%matrix_cturnover_livestem_acc_patch(p)    = cs13_veg%matrix_cturnover_livestem_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivestem)+matrix_gmturnover(p,ilivestem)+matrix_fiturnover(p,ilivestem)) &
                                                                     * cs13_veg%livestemc_patch(p)
                  cs13_veg%matrix_cturnover_livestemst_acc_patch(p)  = cs13_veg%matrix_cturnover_livestemst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivestem_st)+matrix_gmturnover(p,ilivestem_st)+matrix_fiturnover(p,ilivestem_st)) &
                                                                     * cs13_veg%livestemc_storage_patch(p)
                  cs13_veg%matrix_cturnover_livestemxf_acc_patch(p)  = cs13_veg%matrix_cturnover_livestemxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivestem_xf)+matrix_gmturnover(p,ilivestem_xf)+matrix_fiturnover(p,ilivestem_xf)) &
                                                                     * cs13_veg%livestemc_xfer_patch(p)
                  cs13_veg%matrix_cturnover_deadstem_acc_patch(p)    = cs13_veg%matrix_cturnover_deadstem_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadstem)+matrix_gmturnover(p,ideadstem)+matrix_fiturnover(p,ideadstem)) &
                                                                     * cs13_veg%deadstemc_patch(p)
                  cs13_veg%matrix_cturnover_deadstemst_acc_patch(p)  = cs13_veg%matrix_cturnover_deadstemst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadstem_st)+matrix_gmturnover(p,ideadstem_st)+matrix_fiturnover(p,ideadstem_st)) &
                                                                     * cs13_veg%deadstemc_storage_patch(p)
                  cs13_veg%matrix_cturnover_deadstemxf_acc_patch(p)  = cs13_veg%matrix_cturnover_deadstemxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadstem_xf)+matrix_gmturnover(p,ideadstem_xf)+matrix_fiturnover(p,ideadstem_xf)) &
                                                                     * cs13_veg%deadstemc_xfer_patch(p)
                  cs13_veg%matrix_cturnover_livecroot_acc_patch(p)   = cs13_veg%matrix_cturnover_livecroot_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivecroot)+matrix_gmturnover(p,ilivecroot)+matrix_fiturnover(p,ilivecroot)) &
                                                                     * cs13_veg%livecrootc_patch(p)
                  cs13_veg%matrix_cturnover_livecrootst_acc_patch(p) = cs13_veg%matrix_cturnover_livecrootst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivecroot_st)+matrix_gmturnover(p,ilivecroot_st)+matrix_fiturnover(p,ilivecroot_st)) &
                                                                     * cs13_veg%livecrootc_storage_patch(p)
                  cs13_veg%matrix_cturnover_livecrootxf_acc_patch(p) = cs13_veg%matrix_cturnover_livecrootxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivecroot_xf)+matrix_gmturnover(p,ilivecroot_xf)+matrix_fiturnover(p,ilivecroot_xf)) &
                                                                     * cs13_veg%livecrootc_xfer_patch(p)
                  cs13_veg%matrix_cturnover_deadcroot_acc_patch(p)   = cs13_veg%matrix_cturnover_deadcroot_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadcroot)+matrix_gmturnover(p,ideadcroot)+matrix_fiturnover(p,ideadcroot)) & 
                                                                     * cs13_veg%deadcrootc_patch(p)
                  cs13_veg%matrix_cturnover_deadcrootst_acc_patch(p) = cs13_veg%matrix_cturnover_deadcrootst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadcroot_st)+matrix_gmturnover(p,ideadcroot_st)+matrix_fiturnover(p,ideadcroot_st)) & 
                                                                     * cs13_veg%deadcrootc_storage_patch(p)
                  cs13_veg%matrix_cturnover_deadcrootxf_acc_patch(p) = cs13_veg%matrix_cturnover_deadcrootxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadcroot_xf)+matrix_gmturnover(p,ideadcroot_xf)+matrix_fiturnover(p,ideadcroot_xf)) &
                                                                     * cs13_veg%deadcrootc_xfer_patch(p)
               end if
               if(use_c14)then
                  associate( &
                        matrix_c14fitransfer  => c14_cnveg_carbonflux_inst%matrix_fitransfer_patch    , & ! Input:  [real(r8) (:,:)] (gC/m2/s) C transfer rate from fire processes, updated in (CNFireBaseMod or CNFireLi2014Mod) and CNC14decayMod
                        matrix_c14fiturnover  => c14_cnveg_carbonflux_inst%matrix_fiturnover_patch      & ! Output: [real(r8) (:,:)] (gC/m2/step) C turnover rate from fire processe, updated in CNVegMatrixMods
                  )
                  cs14_veg%matrix_ctransfer_leafst_to_leafxf_acc_patch(p)           = cs14_veg%matrix_ctransfer_leafst_to_leafxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ileafst_to_ileafxf_phc) &
                                                                                    * dt * cs14_veg%leafc_storage_patch(p) !matrix_phturnover(p,ileaf_st)*leafc_storage(p)
                  cs14_veg%matrix_ctransfer_leafxf_to_leaf_acc_patch(p)             = cs14_veg%matrix_ctransfer_leafxf_to_leaf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ileafxf_to_ileaf_phc) &
                                                                                    * dt * cs14_veg%leafc_xfer_patch(p)!matrix_phturnover(p,ileaf_xf)*leafc_xfer(p)
                  cs14_veg%matrix_ctransfer_frootst_to_frootxf_acc_patch(p)         = cs14_veg%matrix_ctransfer_frootst_to_frootxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) &
                                                                                    * dt * cs14_veg%frootc_storage_patch(p)!matrix_phturnover(p,ifroot_st)*frootc_storage(p)
                  cs14_veg%matrix_ctransfer_frootxf_to_froot_acc_patch(p)           = cs14_veg%matrix_ctransfer_frootxf_to_froot_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ifrootxf_to_ifroot_phc) &
                                                                                    * dt * cs14_veg%frootc_xfer_patch(p)!matrix_phturnover(p,ifroot_xf)*frootc_xfer(p)
                  cs14_veg%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p)   = cs14_veg%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) &
                                                                                    * dt * cs14_veg%livestemc_storage_patch(p)!matrix_phturnover(p,ilivestem_st)*livestemc_storage(p)
                  cs14_veg%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p)     = cs14_veg%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc) &
                                                                                    * dt * cs14_veg%livestemc_xfer_patch(p)!matrix_phturnover(p,ilivestem_xf)*livestemc_xfer(p)
                  cs14_veg%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p)   = cs14_veg%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) &
                                                                                    * dt * cs14_veg%deadstemc_storage_patch(p)!matrix_phturnover(p,ideadstem_st)*deadstemc_storage(p)
                  cs14_veg%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p)     = cs14_veg%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc) &
                                                                                    * dt * cs14_veg%deadstemc_xfer_patch(p)!matrix_phturnover(p,ideadstem_xf)*deadstemc_xfer(p)
                  cs14_veg%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p) = cs14_veg%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) &
                                                                                    * dt * cs14_veg%livecrootc_storage_patch(p)!matrix_phturnover(p,ilivecroot_st)*livecrootc_storage(p)
                  cs14_veg%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p)   = cs14_veg%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) &
                                                                                    * dt * cs14_veg%livecrootc_xfer_patch(p)!matrix_phturnover(p,ilivecroot_xf)*livecrootc_xfer(p)
                  cs14_veg%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p) = cs14_veg%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) &
                                                                                    * dt * cs14_veg%deadcrootc_storage_patch(p)!matrix_phturnover(p,ideadcroot_st)*deadcrootc_storage(p)
                  cs14_veg%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p)   = cs14_veg%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p) &
                                                                                    + matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) &
                                                                                    * dt * cs14_veg%deadcrootc_xfer_patch(p)!matrix_phturnover(p,ideadcroot_st)*deadcrootc_xfer(p)
                  cs14_veg%matrix_ctransfer_livestem_to_deadstem_acc_patch(p)       = cs14_veg%matrix_ctransfer_livestem_to_deadstem_acc_patch(p) &
                                                                                    +(matrix_phtransfer(p,ilivestem_to_ideadstem_phc)&!matrix_phturnover(p,ilivestem) &
                                                                                    + matrix_c14fitransfer(p,ilivestem_to_ideadstem_fic))&!matrix_fiturnover(p,ilivestem))&
                                                                                    * dt * cs14_veg%livestemc_patch(p) 
                  cs14_veg%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p)     = cs14_veg%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p) &
                                                                                    +(matrix_phtransfer(p,ilivecroot_to_ideadcroot_phc)&!*matrix_phturnover(p,ilivecroot) &
                                                                                    + matrix_c14fitransfer(p,ilivecroot_to_ideadcroot_fic))&!*matrix_fiturnover(p,ilivecroot))&
                                                                                    * dt * cs14_veg%livecrootc_patch(p)
                  cs14_veg%matrix_cturnover_leaf_acc_patch(p)        = cs14_veg%matrix_cturnover_leaf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ileaf)+matrix_gmturnover(p,ileaf)+matrix_c14fiturnover(p,ileaf)) &
                                                                     * cs14_veg%leafc_patch(p)
                  cs14_veg%matrix_cturnover_leafst_acc_patch(p)      = cs14_veg%matrix_cturnover_leafst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ileaf_st)+matrix_gmturnover(p,ileaf_st)+matrix_c14fiturnover(p,ileaf_st)) &
                                                                     * cs14_veg%leafc_storage_patch(p)
                  cs14_veg%matrix_cturnover_leafxf_acc_patch(p)      = cs14_veg%matrix_cturnover_leafxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ileaf_xf)+matrix_gmturnover(p,ileaf_xf)+matrix_c14fiturnover(p,ileaf_xf)) &
                                                                     * cs14_veg%leafc_xfer_patch(p)
                  cs14_veg%matrix_cturnover_froot_acc_patch(p)       = cs14_veg%matrix_cturnover_froot_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ifroot)+matrix_gmturnover(p,ifroot)+matrix_c14fiturnover(p,ifroot)) &
                                                                     * cs14_veg%frootc_patch(p)
                  cs14_veg%matrix_cturnover_frootst_acc_patch(p)     = cs14_veg%matrix_cturnover_frootst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ifroot_st)+matrix_gmturnover(p,ifroot_st)+matrix_c14fiturnover(p,ifroot_st)) &
                                                                     * cs14_veg%frootc_storage_patch(p)
                  cs14_veg%matrix_cturnover_frootxf_acc_patch(p)     = cs14_veg%matrix_cturnover_frootxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ifroot_xf)+matrix_gmturnover(p,ifroot_xf)+matrix_c14fiturnover(p,ifroot_xf)) &
                                                                     * cs14_veg%frootc_xfer_patch(p)
                  cs14_veg%matrix_cturnover_livestem_acc_patch(p)    = cs14_veg%matrix_cturnover_livestem_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivestem)+matrix_gmturnover(p,ilivestem)+matrix_c14fiturnover(p,ilivestem)) &
                                                                     * cs14_veg%livestemc_patch(p)
                  cs14_veg%matrix_cturnover_livestemst_acc_patch(p)  = cs14_veg%matrix_cturnover_livestemst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivestem_st)+matrix_gmturnover(p,ilivestem_st)+matrix_c14fiturnover(p,ilivestem_st)) &
                                                                     * cs14_veg%livestemc_storage_patch(p)
                  cs14_veg%matrix_cturnover_livestemxf_acc_patch(p)  = cs14_veg%matrix_cturnover_livestemxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivestem_xf)+matrix_gmturnover(p,ilivestem_xf)+matrix_c14fiturnover(p,ilivestem_xf)) &
                                                                     * cs14_veg%livestemc_xfer_patch(p)
                  cs14_veg%matrix_cturnover_deadstem_acc_patch(p)    = cs14_veg%matrix_cturnover_deadstem_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadstem)+matrix_gmturnover(p,ideadstem)+matrix_c14fiturnover(p,ideadstem)) &
                                                                     * cs14_veg%deadstemc_patch(p)
                  cs14_veg%matrix_cturnover_deadstemst_acc_patch(p)  = cs14_veg%matrix_cturnover_deadstemst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadstem_st)+matrix_gmturnover(p,ideadstem_st)+matrix_c14fiturnover(p,ideadstem_st)) &
                                                                     * cs14_veg%deadstemc_storage_patch(p)
                  cs14_veg%matrix_cturnover_deadstemxf_acc_patch(p)  = cs14_veg%matrix_cturnover_deadstemxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadstem_xf)+matrix_gmturnover(p,ideadstem_xf)+matrix_c14fiturnover(p,ideadstem_xf)) &
                                                                     * cs14_veg%deadstemc_xfer_patch(p)
                  cs14_veg%matrix_cturnover_livecroot_acc_patch(p)   = cs14_veg%matrix_cturnover_livecroot_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivecroot)+matrix_gmturnover(p,ilivecroot)+matrix_c14fiturnover(p,ilivecroot)) &
                                                                     * cs14_veg%livecrootc_patch(p)
                  cs14_veg%matrix_cturnover_livecrootst_acc_patch(p) = cs14_veg%matrix_cturnover_livecrootst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivecroot_st)+matrix_gmturnover(p,ilivecroot_st)+matrix_c14fiturnover(p,ilivecroot_st)) &
                                                                     * cs14_veg%livecrootc_storage_patch(p)
                  cs14_veg%matrix_cturnover_livecrootxf_acc_patch(p) = cs14_veg%matrix_cturnover_livecrootxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ilivecroot_xf)+matrix_gmturnover(p,ilivecroot_xf)+matrix_c14fiturnover(p,ilivecroot_xf)) &
                                                                     * cs14_veg%livecrootc_xfer_patch(p)
                  cs14_veg%matrix_cturnover_deadcroot_acc_patch(p)   = cs14_veg%matrix_cturnover_deadcroot_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadcroot)+matrix_gmturnover(p,ideadcroot)+matrix_c14fiturnover(p,ideadcroot)) & 
                                                                     * cs14_veg%deadcrootc_patch(p)
                  cs14_veg%matrix_cturnover_deadcrootst_acc_patch(p) = cs14_veg%matrix_cturnover_deadcrootst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadcroot_st)+matrix_gmturnover(p,ideadcroot_st)+matrix_c14fiturnover(p,ideadcroot_st)) & 
                                                                     * cs14_veg%deadcrootc_storage_patch(p)
                  cs14_veg%matrix_cturnover_deadcrootxf_acc_patch(p) = cs14_veg%matrix_cturnover_deadcrootxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,ideadcroot_xf)+matrix_gmturnover(p,ideadcroot_xf)+matrix_c14fiturnover(p,ideadcroot_xf)) &
                                                                     * cs14_veg%deadcrootc_xfer_patch(p)
                  end associate
               end if
            end do
            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if(ivt(p) >= npcropmin)then
                  matrix_cturnover_grain_acc(p)    = matrix_cturnover_grain_acc(p) &
                                                   + (matrix_phturnover(p,igrain)+matrix_gmturnover(p,igrain)+matrix_fiturnover(p,igrain)) &
                                                   * reproductivec(p,irepr)
                  matrix_cturnover_grainst_acc(p)  = matrix_cturnover_grainst_acc(p) &
                                                   + (matrix_phturnover(p,igrain_st)+matrix_gmturnover(p,igrain_st)+matrix_fiturnover(p,igrain_st)) &
                                                   * reproductivec_storage(p,irepr)
                  matrix_cturnover_grainxf_acc(p)  = matrix_cturnover_grainxf_acc(p) &
                                                   + (matrix_phturnover(p,igrain_xf)+matrix_gmturnover(p,igrain_xf)+matrix_fiturnover(p,igrain_xf)) &
                                                   * reproductivec_xfer(p,irepr)
                  if(use_c13)then
                     cs13_veg%matrix_cturnover_grain_acc_patch(p)    = cs13_veg%matrix_cturnover_grain_acc_patch(p) &
                                                                     + (matrix_phturnover(p,igrain)+matrix_gmturnover(p,igrain)+matrix_fiturnover(p,igrain)) &
                                                                     * cs13_veg%reproductivec_patch(p,irepr)
                     cs13_veg%matrix_cturnover_grainst_acc_patch(p)  = cs13_veg%matrix_cturnover_grainst_acc_patch(p) &
                                                                     + (matrix_phturnover(p,igrain_st)+matrix_gmturnover(p,igrain_st)+matrix_fiturnover(p,igrain_st)) &
                                                                     * cs13_veg%reproductivec_storage_patch(p,irepr)
                     cs13_veg%matrix_cturnover_grainxf_acc_patch(p)  = cs13_veg%matrix_cturnover_grainxf_acc_patch(p) &
                                                                     + (matrix_phturnover(p,igrain_xf)+matrix_gmturnover(p,igrain_xf)+matrix_fiturnover(p,igrain_xf)) &
                                                                     * cs13_veg%reproductivec_xfer_patch(p,irepr)
                  end if
                  if(use_c14)then
                     associate( &
                           matrix_c14fitransfer  => c14_cnveg_carbonflux_inst%matrix_fitransfer_patch    , & ! Input:  [real(r8) (:,:)] (gC/m2/s) C transfer rate from fire processes, updated in (CNFireBaseMod or CNFireLi2014Mod) and CNC14decayMod
                           matrix_c14fiturnover  => c14_cnveg_carbonflux_inst%matrix_fiturnover_patch      & ! Output: [real(r8) (:,:)] (gC/m2/step) C turnover rate from fire processe, updated in CNVegMatrixMods
                     )
                     cs14_veg%matrix_cturnover_grain_acc_patch(p)    = cs14_veg%matrix_cturnover_grain_acc_patch(p) &
                                                               + (matrix_phturnover(p,igrain)+matrix_gmturnover(p,igrain)+matrix_c14fiturnover(p,igrain)) &
                                                               * cs14_veg%reproductivec_patch(p,irepr)
                     cs14_veg%matrix_cturnover_grainst_acc_patch(p)  = cs14_veg%matrix_cturnover_grainst_acc_patch(p) &
                                                               + (matrix_phturnover(p,igrain_st)+matrix_gmturnover(p,igrain_st)+matrix_c14fiturnover(p,igrain_st)) &
                                                               * cs14_veg%reproductivec_storage_patch(p,irepr)
                     cs14_veg%matrix_cturnover_grainxf_acc_patch(p)  = cs14_veg%matrix_cturnover_grainxf_acc_patch(p) &
                                                               + (matrix_phturnover(p,igrain_xf)+matrix_gmturnover(p,igrain_xf)+matrix_c14fiturnover(p,igrain_xf)) &
                                                               * cs14_veg%reproductivec_xfer_patch(p,irepr)
                     end associate
                  end if
               end if
            end do

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
                                                      * deadstemn_xfer(p)
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
                                                      + (matrix_nphturnover(p,igrain)+matrix_ngmturnover(p,igrain)+matrix_nfiturnover(p,igrain)) &
                                                      * reproductiven(p,irepr)
                  matrix_nturnover_grainst_acc(p)     = matrix_nturnover_grainst_acc(p) &
                                                      + (matrix_nphturnover(p,igrain_st)+matrix_ngmturnover(p,igrain_st)+matrix_nfiturnover(p,igrain_st)) &
                                                      * reproductiven_storage(p,irepr)
                  matrix_nturnover_grainxf_acc(p)     = matrix_nturnover_grainxf_acc(p) &
                                                      + (matrix_nphturnover(p,igrain_xf)+matrix_ngmturnover(p,igrain_xf)+matrix_nfiturnover(p,igrain_xf)) & 
                                                      * reproductiven_xfer(p,irepr)
               end if
            end do
         end if
         call t_stopf('CN veg matrix-accum. trans.')

  ! Update state variables
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
               ! NOTE: This assumes only a single grain pool! (i.e nrepr is
               ! fixed at 1)!
               reproductivec(p,:)            = Xvegc%V(p,igrain)
               reproductivec_storage(p,:)    = Xvegc%V(p,igrain_st)
               reproductivec_xfer(p,:)       = Xvegc%V(p,igrain_xf)
            end if
         end do
         
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
                  ! NOTE: This assumes only a single grain pool! (i.e nrepr is
                  ! fixed at 1)!
                  cs13_veg%reproductivec_patch(p,:)            = Xveg13c%V(p,igrain)
                  cs13_veg%reproductivec_storage_patch(p,:)    = Xveg13c%V(p,igrain_st)
                  cs13_veg%reproductivec_xfer_patch(p,:)       = Xveg13c%V(p,igrain_xf)
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
                  ! NOTE: This assumes only a single grain pool! (i.e nrepr is
                  ! fixed at 1)!
                  cs14_veg%reproductivec_patch(p,:)            = Xveg14c%V(p,igrain)
                  cs14_veg%reproductivec_storage_patch(p,:)    = Xveg14c%V(p,igrain_st)
                  cs14_veg%reproductivec_xfer_patch(p,:)       = Xveg14c%V(p,igrain_xf)
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
               reproductiven(p,:)            = Xvegn%V(p,igrain)
               reproductiven_storage(p,:)    = Xvegn%V(p,igrain_st)
               reproductiven_xfer(p,:)       = Xvegn%V(p,igrain_xf)
            end if
         end do
         call t_stopf('CN veg matrix-assign new value')

  ! Calculate C storage capacity. 2D matrix instead of sparse matrix is still used when calculating the inverse
         if(spinup_matrixcn .or. hist_wrt_matrixcn_diag)then
            if((.not. spinup_matrixcn .and. is_end_curr_year()) .or. (spinup_matrixcn .and. is_end_curr_year() .and. mod(iyr,nyr_SASU) .eq. 0))then
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

                  if(use_c13)then
                     matrix_c13alloc_acc(ileaf)         = cs13_veg%matrix_calloc_leaf_acc_patch(p)               
                     matrix_c13alloc_acc(ileaf_st)      = cs13_veg%matrix_calloc_leafst_acc_patch(p)               
                     matrix_c13alloc_acc(ifroot)        = cs13_veg%matrix_calloc_froot_acc_patch(p)               
                     matrix_c13alloc_acc(ifroot_st)     = cs13_veg%matrix_calloc_frootst_acc_patch(p)               
                     matrix_c13alloc_acc(ilivestem)     = cs13_veg%matrix_calloc_livestem_acc_patch(p)               
                     matrix_c13alloc_acc(ilivestem_st)  = cs13_veg%matrix_calloc_livestemst_acc_patch(p)               
                     matrix_c13alloc_acc(ideadstem)     = cs13_veg%matrix_calloc_deadstem_acc_patch(p)               
                     matrix_c13alloc_acc(ideadstem_st)  = cs13_veg%matrix_calloc_deadstemst_acc_patch(p)               
                     matrix_c13alloc_acc(ilivecroot)    = cs13_veg%matrix_calloc_livecroot_acc_patch(p)               
                     matrix_c13alloc_acc(ilivecroot_st) = cs13_veg%matrix_calloc_livecrootst_acc_patch(p)               
                     matrix_c13alloc_acc(ideadcroot)    = cs13_veg%matrix_calloc_deadcroot_acc_patch(p)               
                     matrix_c13alloc_acc(ideadcroot_st) = cs13_veg%matrix_calloc_deadcrootst_acc_patch(p)               
                     if(ivt(p) >= npcropmin)then
                        matrix_c13alloc_acc(igrain)     = cs13_veg%matrix_calloc_grain_acc_patch(p)               
                        matrix_c13alloc_acc(igrain_st)  = cs13_veg%matrix_calloc_grainst_acc_patch(p)               
                     end if

                     matrix_c13transfer_acc(ileaf_xf,ileaf_st)           = cs13_veg%matrix_ctransfer_leafst_to_leafxf_acc_patch(p)
                     matrix_c13transfer_acc(ileaf,ileaf_xf)              = cs13_veg%matrix_ctransfer_leafxf_to_leaf_acc_patch(p)
                     matrix_c13transfer_acc(ifroot_xf,ifroot_st)         = cs13_veg%matrix_ctransfer_frootst_to_frootxf_acc_patch(p)
                     matrix_c13transfer_acc(ifroot,ifroot_xf)            = cs13_veg%matrix_ctransfer_frootxf_to_froot_acc_patch(p)
                     matrix_c13transfer_acc(ilivestem_xf,ilivestem_st)   = cs13_veg%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p)
                     matrix_c13transfer_acc(ilivestem,ilivestem_xf)      = cs13_veg%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p)
                     matrix_c13transfer_acc(ideadstem_xf,ideadstem_st)   = cs13_veg%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p)
                     matrix_c13transfer_acc(ideadstem,ideadstem_xf)      = cs13_veg%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p)
                     matrix_c13transfer_acc(ilivecroot_xf,ilivecroot_st) = cs13_veg%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p)
                     matrix_c13transfer_acc(ilivecroot,ilivecroot_xf)    = cs13_veg%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p)
                     matrix_c13transfer_acc(ideadcroot_xf,ideadcroot_st) = cs13_veg%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p)
                     matrix_c13transfer_acc(ideadcroot,ideadcroot_xf)    = cs13_veg%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p)
                     if(ivt(p) >= npcropmin)then
                        matrix_c13transfer_acc(igrain_xf,igrain_st)      = cs13_veg%matrix_ctransfer_grainst_to_grainxf_acc_patch(p)
                        matrix_c13transfer_acc(igrain,igrain_xf)         = cs13_veg%matrix_ctransfer_grainxf_to_grain_acc_patch(p)
                     end if
                     matrix_c13transfer_acc(ideadstem,ilivestem)         = cs13_veg%matrix_ctransfer_livestem_to_deadstem_acc_patch(p)
                     matrix_c13transfer_acc(ideadcroot,ilivecroot)       = cs13_veg%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p)
                  
                     matrix_c13transfer_acc(ileaf,ileaf)                 = -cs13_veg%matrix_cturnover_leaf_acc_patch(p)               
                     matrix_c13transfer_acc(ileaf_st,ileaf_st)           = -cs13_veg%matrix_cturnover_leafst_acc_patch(p)               
                     matrix_c13transfer_acc(ileaf_xf,ileaf_xf)           = -cs13_veg%matrix_cturnover_leafxf_acc_patch(p)               
                     matrix_c13transfer_acc(ifroot,ifroot)               = -cs13_veg%matrix_cturnover_froot_acc_patch(p)               
                     matrix_c13transfer_acc(ifroot_st,ifroot_st)         = -cs13_veg%matrix_cturnover_frootst_acc_patch(p)               
                     matrix_c13transfer_acc(ifroot_xf,ifroot_xf)         = -cs13_veg%matrix_cturnover_frootxf_acc_patch(p)               
                     matrix_c13transfer_acc(ilivestem,ilivestem)         = -cs13_veg%matrix_cturnover_livestem_acc_patch(p)               
                     matrix_c13transfer_acc(ilivestem_st,ilivestem_st)   = -cs13_veg%matrix_cturnover_livestemst_acc_patch(p)               
                     matrix_c13transfer_acc(ilivestem_xf,ilivestem_xf)   = -cs13_veg%matrix_cturnover_livestemxf_acc_patch(p)               
                     matrix_c13transfer_acc(ideadstem,ideadstem)         = -cs13_veg%matrix_cturnover_deadstem_acc_patch(p)               
                     matrix_c13transfer_acc(ideadstem_st,ideadstem_st)   = -cs13_veg%matrix_cturnover_deadstemst_acc_patch(p)               
                     matrix_c13transfer_acc(ideadstem_xf,ideadstem_xf)   = -cs13_veg%matrix_cturnover_deadstemxf_acc_patch(p)               
                     matrix_c13transfer_acc(ilivecroot,ilivecroot)       = -cs13_veg%matrix_cturnover_livecroot_acc_patch(p)               
                     matrix_c13transfer_acc(ilivecroot_st,ilivecroot_st) = -cs13_veg%matrix_cturnover_livecrootst_acc_patch(p)               
                     matrix_c13transfer_acc(ilivecroot_xf,ilivecroot_xf) = -cs13_veg%matrix_cturnover_livecrootxf_acc_patch(p)               
                     matrix_c13transfer_acc(ideadcroot,ideadcroot)       = -cs13_veg%matrix_cturnover_deadcroot_acc_patch(p)               
                     matrix_c13transfer_acc(ideadcroot_st,ideadcroot_st) = -cs13_veg%matrix_cturnover_deadcrootst_acc_patch(p)               
                     matrix_c13transfer_acc(ideadcroot_xf,ideadcroot_xf) = -cs13_veg%matrix_cturnover_deadcrootxf_acc_patch(p)               
                     if(ivt(p) >= npcropmin)then
                        matrix_c13transfer_acc(igrain,igrain)            = -cs13_veg%matrix_cturnover_grain_acc_patch(p)               
                        matrix_c13transfer_acc(igrain_st,igrain_st)      = -cs13_veg%matrix_cturnover_grainst_acc_patch(p)               
                        matrix_c13transfer_acc(igrain_xf,igrain_xf)      = -cs13_veg%matrix_cturnover_grainxf_acc_patch(p)               
                     end if
                  end if

                  if(use_c14)then
                     matrix_c14alloc_acc(ileaf)         = cs14_veg%matrix_calloc_leaf_acc_patch(p)               
                     matrix_c14alloc_acc(ileaf_st)      = cs14_veg%matrix_calloc_leafst_acc_patch(p)               
                     matrix_c14alloc_acc(ifroot)        = cs14_veg%matrix_calloc_froot_acc_patch(p)               
                     matrix_c14alloc_acc(ifroot_st)     = cs14_veg%matrix_calloc_frootst_acc_patch(p)               
                     matrix_c14alloc_acc(ilivestem)     = cs14_veg%matrix_calloc_livestem_acc_patch(p)               
                     matrix_c14alloc_acc(ilivestem_st)  = cs14_veg%matrix_calloc_livestemst_acc_patch(p)               
                     matrix_c14alloc_acc(ideadstem)     = cs14_veg%matrix_calloc_deadstem_acc_patch(p)               
                     matrix_c14alloc_acc(ideadstem_st)  = cs14_veg%matrix_calloc_deadstemst_acc_patch(p)               
                     matrix_c14alloc_acc(ilivecroot)    = cs14_veg%matrix_calloc_livecroot_acc_patch(p)               
                     matrix_c14alloc_acc(ilivecroot_st) = cs14_veg%matrix_calloc_livecrootst_acc_patch(p)               
                     matrix_c14alloc_acc(ideadcroot)    = cs14_veg%matrix_calloc_deadcroot_acc_patch(p)               
                     matrix_c14alloc_acc(ideadcroot_st) = cs14_veg%matrix_calloc_deadcrootst_acc_patch(p)               
                     if(ivt(p) >= npcropmin)then
                        matrix_c14alloc_acc(igrain)     = cs14_veg%matrix_calloc_grain_acc_patch(p)               
                        matrix_c14alloc_acc(igrain_st)  = cs14_veg%matrix_calloc_grainst_acc_patch(p)               
                     end if

                     matrix_c14transfer_acc(ileaf_xf,ileaf_st)           = cs14_veg%matrix_ctransfer_leafst_to_leafxf_acc_patch(p)
                     matrix_c14transfer_acc(ileaf,ileaf_xf)              = cs14_veg%matrix_ctransfer_leafxf_to_leaf_acc_patch(p)
                     matrix_c14transfer_acc(ifroot_xf,ifroot_st)         = cs14_veg%matrix_ctransfer_frootst_to_frootxf_acc_patch(p)
                     matrix_c14transfer_acc(ifroot,ifroot_xf)            = cs14_veg%matrix_ctransfer_frootxf_to_froot_acc_patch(p)
                     matrix_c14transfer_acc(ilivestem_xf,ilivestem_st)   = cs14_veg%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p)
                     matrix_c14transfer_acc(ilivestem,ilivestem_xf)      = cs14_veg%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p)
                     matrix_c14transfer_acc(ideadstem_xf,ideadstem_st)   = cs14_veg%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p)
                     matrix_c14transfer_acc(ideadstem,ideadstem_xf)      = cs14_veg%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p)
                     matrix_c14transfer_acc(ilivecroot_xf,ilivecroot_st) = cs14_veg%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p)
                     matrix_c14transfer_acc(ilivecroot,ilivecroot_xf)    = cs14_veg%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p)
                     matrix_c14transfer_acc(ideadcroot_xf,ideadcroot_st) = cs14_veg%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p)
                     matrix_c14transfer_acc(ideadcroot,ideadcroot_xf)    = cs14_veg%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p)
                     if(ivt(p) >= npcropmin)then
                        matrix_c14transfer_acc(igrain_xf,igrain_st)      = cs14_veg%matrix_ctransfer_grainst_to_grainxf_acc_patch(p)
                        matrix_c14transfer_acc(igrain,igrain_xf)         = cs14_veg%matrix_ctransfer_grainxf_to_grain_acc_patch(p)
                     end if
                     matrix_c14transfer_acc(ideadstem,ilivestem)         = cs14_veg%matrix_ctransfer_livestem_to_deadstem_acc_patch(p)
                     matrix_c14transfer_acc(ideadcroot,ilivecroot)       = cs14_veg%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p)
                  
                     matrix_c14transfer_acc(ileaf,ileaf)                 = -cs14_veg%matrix_cturnover_leaf_acc_patch(p)               
                     matrix_c14transfer_acc(ileaf_st,ileaf_st)           = -cs14_veg%matrix_cturnover_leafst_acc_patch(p)               
                     matrix_c14transfer_acc(ileaf_xf,ileaf_xf)           = -cs14_veg%matrix_cturnover_leafxf_acc_patch(p)               
                     matrix_c14transfer_acc(ifroot,ifroot)               = -cs14_veg%matrix_cturnover_froot_acc_patch(p)               
                     matrix_c14transfer_acc(ifroot_st,ifroot_st)         = -cs14_veg%matrix_cturnover_frootst_acc_patch(p)               
                     matrix_c14transfer_acc(ifroot_xf,ifroot_xf)         = -cs14_veg%matrix_cturnover_frootxf_acc_patch(p)               
                     matrix_c14transfer_acc(ilivestem,ilivestem)         = -cs14_veg%matrix_cturnover_livestem_acc_patch(p)               
                     matrix_c14transfer_acc(ilivestem_st,ilivestem_st)   = -cs14_veg%matrix_cturnover_livestemst_acc_patch(p)               
                     matrix_c14transfer_acc(ilivestem_xf,ilivestem_xf)   = -cs14_veg%matrix_cturnover_livestemxf_acc_patch(p)               
                     matrix_c14transfer_acc(ideadstem,ideadstem)         = -cs14_veg%matrix_cturnover_deadstem_acc_patch(p)               
                     matrix_c14transfer_acc(ideadstem_st,ideadstem_st)   = -cs14_veg%matrix_cturnover_deadstemst_acc_patch(p)               
                     matrix_c14transfer_acc(ideadstem_xf,ideadstem_xf)   = -cs14_veg%matrix_cturnover_deadstemxf_acc_patch(p)               
                     matrix_c14transfer_acc(ilivecroot,ilivecroot)       = -cs14_veg%matrix_cturnover_livecroot_acc_patch(p)               
                     matrix_c14transfer_acc(ilivecroot_st,ilivecroot_st) = -cs14_veg%matrix_cturnover_livecrootst_acc_patch(p)               
                     matrix_c14transfer_acc(ilivecroot_xf,ilivecroot_xf) = -cs14_veg%matrix_cturnover_livecrootxf_acc_patch(p)               
                     matrix_c14transfer_acc(ideadcroot,ideadcroot)       = -cs14_veg%matrix_cturnover_deadcroot_acc_patch(p)               
                     matrix_c14transfer_acc(ideadcroot_st,ideadcroot_st) = -cs14_veg%matrix_cturnover_deadcrootst_acc_patch(p)               
                     matrix_c14transfer_acc(ideadcroot_xf,ideadcroot_xf) = -cs14_veg%matrix_cturnover_deadcrootxf_acc_patch(p)               
                     if(ivt(p) >= npcropmin)then
                        matrix_c14transfer_acc(igrain,igrain)            = -cs14_veg%matrix_cturnover_grain_acc_patch(p)               
                        matrix_c14transfer_acc(igrain_st,igrain_st)      = -cs14_veg%matrix_cturnover_grainst_acc_patch(p)               
                        matrix_c14transfer_acc(igrain_xf,igrain_xf)      = -cs14_veg%matrix_cturnover_grainxf_acc_patch(p)               
                     end if
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
                    if(matrix_ctransfer_acc(i,i) == 0)then
                       matrix_ctransfer_acc(i,i) = spval
                    end if
                  end do
                  if(use_c13)then
                     do i=1,nvegcpool
                       if(matrix_c13transfer_acc(i,i) == 0)then
                          matrix_c13transfer_acc(i,i) = spval
                       end if
                     end do
                  end if
                  if(use_c14)then
                     do i=1,nvegcpool
                       if(matrix_c14transfer_acc(i,i) == 0)then
                          matrix_c14transfer_acc(i,i) = spval
                       end if
                     end do
                  end if
                  do i=1,nvegnpool
                    if(matrix_ntransfer_acc(i,i) == 0)then
                       matrix_ntransfer_acc(i,i) = spval
                    end if
                  end do

  ! Calculate the transfer rate based on the initial value of the calendar year.
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
                     matrix_ctransfer_acc(1:nvegcpool,igrain)     = matrix_ctransfer_acc(1:nvegcpool,igrain)    / reproc0(p)
                     matrix_ctransfer_acc(1:nvegcpool,igrain_st)  = matrix_ctransfer_acc(1:nvegcpool,igrain_st) / reproc0_storage(p)
                     matrix_ctransfer_acc(1:nvegcpool,igrain_xf)  = matrix_ctransfer_acc(1:nvegcpool,igrain_xf) / reproc0_xfer(p)
                  end if
   
                  if(use_c13)then
                     matrix_c13transfer_acc(1:nvegcpool,ileaf)         = matrix_c13transfer_acc(1:nvegcpool,ileaf)         / cs13_veg%leafc0_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ileaf_st)      = matrix_c13transfer_acc(1:nvegcpool,ileaf_st)      / cs13_veg%leafc0_storage_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ileaf_xf)      = matrix_c13transfer_acc(1:nvegcpool,ileaf_xf)      / cs13_veg%leafc0_xfer_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ifroot)        = matrix_c13transfer_acc(1:nvegcpool,ifroot)        / cs13_veg%frootc0_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ifroot_st)     = matrix_c13transfer_acc(1:nvegcpool,ifroot_st)     / cs13_veg%frootc0_storage_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ifroot_xf)     = matrix_c13transfer_acc(1:nvegcpool,ifroot_xf)     / cs13_veg%frootc0_xfer_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ilivestem)     = matrix_c13transfer_acc(1:nvegcpool,ilivestem)     / cs13_veg%livestemc0_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ilivestem_st)  = matrix_c13transfer_acc(1:nvegcpool,ilivestem_st)  / cs13_veg%livestemc0_storage_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ilivestem_xf)  = matrix_c13transfer_acc(1:nvegcpool,ilivestem_xf)  / cs13_veg%livestemc0_xfer_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ideadstem)     = matrix_c13transfer_acc(1:nvegcpool,ideadstem)     / cs13_veg%deadstemc0_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ideadstem_st)  = matrix_c13transfer_acc(1:nvegcpool,ideadstem_st)  / cs13_veg%deadstemc0_storage_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ideadstem_xf)  = matrix_c13transfer_acc(1:nvegcpool,ideadstem_xf)  / cs13_veg%deadstemc0_xfer_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ilivecroot)    = matrix_c13transfer_acc(1:nvegcpool,ilivecroot)    / cs13_veg%livecrootc0_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ilivecroot_st) = matrix_c13transfer_acc(1:nvegcpool,ilivecroot_st) / cs13_veg%livecrootc0_storage_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ilivecroot_xf) = matrix_c13transfer_acc(1:nvegcpool,ilivecroot_xf) / cs13_veg%livecrootc0_xfer_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ideadcroot)    = matrix_c13transfer_acc(1:nvegcpool,ideadcroot)    / cs13_veg%deadcrootc0_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ideadcroot_st) = matrix_c13transfer_acc(1:nvegcpool,ideadcroot_st) / cs13_veg%deadcrootc0_storage_patch(p)
                     matrix_c13transfer_acc(1:nvegcpool,ideadcroot_xf) = matrix_c13transfer_acc(1:nvegcpool,ideadcroot_xf) / cs13_veg%deadcrootc0_xfer_patch(p)
                     if(ivt(p) >= npcropmin)then
                        matrix_c13transfer_acc(1:nvegcpool,igrain)     = matrix_c13transfer_acc(1:nvegcpool,igrain)    / cs13_veg%reproc0_patch(p)
                        matrix_c13transfer_acc(1:nvegcpool,igrain_st)  = matrix_c13transfer_acc(1:nvegcpool,igrain_st) / cs13_veg%reproc0_storage_patch(p)
                        matrix_c13transfer_acc(1:nvegcpool,igrain_xf)  = matrix_c13transfer_acc(1:nvegcpool,igrain_xf) / cs13_veg%reproc0_xfer_patch(p)
                     end if
                  end if
   
                  if(use_c14)then
                     matrix_c14transfer_acc(1:nvegcpool,ileaf)         = matrix_c14transfer_acc(1:nvegcpool,ileaf)         / cs14_veg%leafc0_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ileaf_st)      = matrix_c14transfer_acc(1:nvegcpool,ileaf_st)      / cs14_veg%leafc0_storage_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ileaf_xf)      = matrix_c14transfer_acc(1:nvegcpool,ileaf_xf)      / cs14_veg%leafc0_xfer_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ifroot)        = matrix_c14transfer_acc(1:nvegcpool,ifroot)        / cs14_veg%frootc0_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ifroot_st)     = matrix_c14transfer_acc(1:nvegcpool,ifroot_st)     / cs14_veg%frootc0_storage_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ifroot_xf)     = matrix_c14transfer_acc(1:nvegcpool,ifroot_xf)     / cs14_veg%frootc0_xfer_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ilivestem)     = matrix_c14transfer_acc(1:nvegcpool,ilivestem)     / cs14_veg%livestemc0_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ilivestem_st)  = matrix_c14transfer_acc(1:nvegcpool,ilivestem_st)  / cs14_veg%livestemc0_storage_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ilivestem_xf)  = matrix_c14transfer_acc(1:nvegcpool,ilivestem_xf)  / cs14_veg%livestemc0_xfer_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ideadstem)     = matrix_c14transfer_acc(1:nvegcpool,ideadstem)     / cs14_veg%deadstemc0_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ideadstem_st)  = matrix_c14transfer_acc(1:nvegcpool,ideadstem_st)  / cs14_veg%deadstemc0_storage_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ideadstem_xf)  = matrix_c14transfer_acc(1:nvegcpool,ideadstem_xf)  / cs14_veg%deadstemc0_xfer_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ilivecroot)    = matrix_c14transfer_acc(1:nvegcpool,ilivecroot)    / cs14_veg%livecrootc0_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ilivecroot_st) = matrix_c14transfer_acc(1:nvegcpool,ilivecroot_st) / cs14_veg%livecrootc0_storage_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ilivecroot_xf) = matrix_c14transfer_acc(1:nvegcpool,ilivecroot_xf) / cs14_veg%livecrootc0_xfer_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ideadcroot)    = matrix_c14transfer_acc(1:nvegcpool,ideadcroot)    / cs14_veg%deadcrootc0_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ideadcroot_st) = matrix_c14transfer_acc(1:nvegcpool,ideadcroot_st) / cs14_veg%deadcrootc0_storage_patch(p)
                     matrix_c14transfer_acc(1:nvegcpool,ideadcroot_xf) = matrix_c14transfer_acc(1:nvegcpool,ideadcroot_xf) / cs14_veg%deadcrootc0_xfer_patch(p)
                     if(ivt(p) >= npcropmin)then
                        matrix_c14transfer_acc(1:nvegcpool,igrain)     = matrix_c14transfer_acc(1:nvegcpool,igrain)    / cs14_veg%reproc0_patch(p)
                        matrix_c14transfer_acc(1:nvegcpool,igrain_st)  = matrix_c14transfer_acc(1:nvegcpool,igrain_st) / cs14_veg%reproc0_storage_patch(p)
                        matrix_c14transfer_acc(1:nvegcpool,igrain_xf)  = matrix_c14transfer_acc(1:nvegcpool,igrain_xf) / cs14_veg%reproc0_xfer_patch(p)
                     end if
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
                     matrix_ntransfer_acc(1:nvegnpool,igrain)     = matrix_ntransfer_acc(1:nvegnpool,igrain)    / repron0(p)
                     matrix_ntransfer_acc(1:nvegnpool,igrain_st)  = matrix_ntransfer_acc(1:nvegnpool,igrain_st) / repron0_storage(p)
                     matrix_ntransfer_acc(1:nvegnpool,igrain_xf)  = matrix_ntransfer_acc(1:nvegnpool,igrain_xf) / repron0_xfer(p)
                  end if
                  matrix_ntransfer_acc(1:nvegnpool,iretransn)  = matrix_ntransfer_acc(1:nvegnpool,iretransn) / retransn0(p)
                
                  call t_stopf('CN veg matrix-prepare AK^-1')
                  call t_startf('CN veg matrix-inv matrix operation')

  ! Calculate the residence time and C storage capacity
                  call inverse(matrix_ctransfer_acc(1:nvegcpool,1:nvegcpool),AKinvc(1:nvegcpool,1:nvegcpool),nvegcpool)
                  vegmatrixc_rt(:) = -matmul(AKinvc(1:nvegcpool,1:nvegcpool),matrix_calloc_acc(1:nvegcpool))

  ! Calculate the residence time and C13 storage capacity
                  if(use_c13)then
                     call inverse(matrix_c13transfer_acc(1:nvegcpool,1:nvegcpool),AKinvc(1:nvegcpool,1:nvegcpool),nvegcpool)
                     vegmatrixc13_rt(:) = -matmul(AKinvc(1:nvegcpool,1:nvegcpool),matrix_c13alloc_acc(1:nvegcpool))
                  end if

  ! Calculate the residence time and C14 storage capacity
                  if(use_c14)then
                     call inverse(matrix_c14transfer_acc(1:nvegcpool,1:nvegcpool),AKinvc(1:nvegcpool,1:nvegcpool),nvegcpool)
                     vegmatrixc14_rt(:) = -matmul(AKinvc(1:nvegcpool,1:nvegcpool),matrix_c14alloc_acc(1:nvegcpool))
                  end if

  ! Calculate the residence time and N storage capacity
                  call inverse(matrix_ntransfer_acc(1:nvegnpool,1:nvegnpool),AKinvn(1:nvegnpool,1:nvegnpool),nvegnpool)
                  vegmatrixn_rt(:) = -matmul(AKinvn(1:nvegnpool,1:nvegnpool),matrix_nalloc_acc(1:nvegnpool))

                  do i=1,nvegcpool
                     if(vegmatrixc_rt(i) .lt. 0)vegmatrixc_rt(i) = epsi
                  end do
                  if(use_c13)then
                     do i=1,nvegcpool
                        if(vegmatrixc13_rt(i) .lt. 0)vegmatrixc13_rt(i) = epsi
                     end do
                  end if
                  if(use_c14)then
                     do i=1,nvegcpool
                        if(vegmatrixc14_rt(i) .lt. 0)vegmatrixc14_rt(i) = epsi
                     end do
                  end if
                  do i=1,nvegnpool
                     if(vegmatrixn_rt(i) .lt. 0)vegmatrixn_rt(i) = epsi
                  end do

                  call t_stopf('CN veg matrix-inv matrix operation')
 
                  call t_startf('CN veg matrix-finalize spinup')
                  
                  if(spinup_matrixcn .and. .not. is_first_step_of_this_run_segment())then
                     deadstemc(p)              = vegmatrixc_rt(ideadstem)
                     deadstemc_storage(p)      = vegmatrixc_rt(ideadstem_st)
                     deadcrootc(p)             = vegmatrixc_rt(ideadcroot)
                     deadcrootc_storage(p)     = vegmatrixc_rt(ideadcroot_st)
                     if(use_c13)then
                        cs13_veg%deadstemc_patch(p)              = vegmatrixc13_rt(ideadstem)
                        cs13_veg%deadstemc_storage_patch(p)      = vegmatrixc13_rt(ideadstem_st)
                        cs13_veg%deadcrootc_patch(p)             = vegmatrixc13_rt(ideadcroot)
                        cs13_veg%deadcrootc_storage_patch(p)     = vegmatrixc13_rt(ideadcroot_st)
                     end if
                     if(use_c14)then
                        cs14_veg%deadstemc_patch(p)              = vegmatrixc14_rt(ideadstem)
                        cs14_veg%deadstemc_storage_patch(p)      = vegmatrixc14_rt(ideadstem_st)
                        cs14_veg%deadcrootc_patch(p)             = vegmatrixc14_rt(ideadcroot)
                        cs14_veg%deadcrootc_storage_patch(p)     = vegmatrixc14_rt(ideadcroot_st)
                     end if
                     deadstemn(p)              = vegmatrixn_rt(ideadstem)
                     deadstemn_storage(p)      = vegmatrixn_rt(ideadstem_st)
                     deadcrootn(p)             = vegmatrixn_rt(ideadcroot)
                     deadcrootn_storage(p)     = vegmatrixn_rt(ideadcroot_st)

                     if(iloop .eq. iloop_avg)then
                        leafc_SASUsave(p)              = leafc_SASUsave(p)              + leafc(p)
                        leafc_storage_SASUsave(p)      = leafc_storage_SASUsave(p)      + leafc_storage(p)
                        leafc_xfer_SASUsave(p)         = leafc_xfer_SASUsave(p)         + leafc_xfer(p)
                        frootc_SASUsave(p)             = frootc_SASUsave(p)             + frootc(p)
                        frootc_storage_SASUsave(p)     = frootc_storage_SASUsave(p)     + frootc_storage(p)
                        frootc_xfer_SASUsave(p)        = frootc_xfer_SASUsave(p)        + frootc_xfer(p)
                        livestemc_SASUsave(p)          = livestemc_SASUsave(p)          + livestemc(p)
                        livestemc_storage_SASUsave(p)  = livestemc_storage_SASUsave(p)  + livestemc_storage(p)
                        livestemc_xfer_SASUsave(p)     = livestemc_xfer_SASUsave(p)     + livestemc_xfer(p)
                        deadstemc_SASUsave(p)          = deadstemc_SASUsave(p)          + deadstemc(p)
                        deadstemc_storage_SASUsave(p)  = deadstemc_storage_SASUsave(p)  + deadstemc_storage(p)
                        deadstemc_xfer_SASUsave(p)     = deadstemc_xfer_SASUsave(p)     + deadstemc_xfer(p)
                        livecrootc_SASUsave(p)         = livecrootc_SASUsave(p)         + livecrootc(p)
                        livecrootc_storage_SASUsave(p) = livecrootc_storage_SASUsave(p) + livecrootc_storage(p)
                        livecrootc_xfer_SASUsave(p)    = livecrootc_xfer_SASUsave(p)    + livecrootc_xfer(p)
                        deadcrootc_SASUsave(p)         = deadcrootc_SASUsave(p)         + deadcrootc(p)
                        deadcrootc_storage_SASUsave(p) = deadcrootc_storage_SASUsave(p) + deadcrootc_storage(p)
                        deadcrootc_xfer_SASUsave(p)    = deadcrootc_xfer_SASUsave(p)    + deadcrootc_xfer(p)
                        if(ivt(p)  >= npcropmin)then
                           grainc_SASUsave(p)          = grainc_SASUsave(p)             + sum(reproductivec(p,:))
                           grainc_storage_SASUsave(p)  = grainc_storage_SASUsave(p)     + sum(reproductivec_storage(p,:))
                        end if
                        if(use_c13)then
                           cs13_veg%leafc_SASUsave_patch(p)              = cs13_veg%leafc_SASUsave_patch(p)              + cs13_veg%leafc_patch(p)
                           cs13_veg%leafc_storage_SASUsave_patch(p)      = cs13_veg%leafc_storage_SASUsave_patch(p)      + cs13_veg%leafc_storage_patch(p)
                           cs13_veg%leafc_xfer_SASUsave_patch(p)         = cs13_veg%leafc_xfer_SASUsave_patch(p)         + cs13_veg%leafc_xfer_patch(p)
                           cs13_veg%frootc_SASUsave_patch(p)             = cs13_veg%frootc_SASUsave_patch(p)             + cs13_veg%frootc_patch(p)
                           cs13_veg%frootc_storage_SASUsave_patch(p)     = cs13_veg%frootc_storage_SASUsave_patch(p)     + cs13_veg%frootc_storage_patch(p)
                           cs13_veg%frootc_xfer_SASUsave_patch(p)        = cs13_veg%frootc_xfer_SASUsave_patch(p)        + cs13_veg%frootc_xfer_patch(p)
                           cs13_veg%livestemc_SASUsave_patch(p)          = cs13_veg%livestemc_SASUsave_patch(p)          + cs13_veg%livestemc_patch(p)
                           cs13_veg%livestemc_storage_SASUsave_patch(p)  = cs13_veg%livestemc_storage_SASUsave_patch(p)  + cs13_veg%livestemc_storage_patch(p)
                           cs13_veg%livestemc_xfer_SASUsave_patch(p)     = cs13_veg%livestemc_xfer_SASUsave_patch(p)     + cs13_veg%livestemc_xfer_patch(p)
                           cs13_veg%deadstemc_SASUsave_patch(p)          = cs13_veg%deadstemc_SASUsave_patch(p)          + cs13_veg%deadstemc_patch(p)
                           cs13_veg%deadstemc_storage_SASUsave_patch(p)  = cs13_veg%deadstemc_storage_SASUsave_patch(p)  + cs13_veg%deadstemc_storage_patch(p)
                           cs13_veg%deadstemc_xfer_SASUsave_patch(p)     = cs13_veg%deadstemc_xfer_SASUsave_patch(p)     + cs13_veg%deadstemc_xfer_patch(p)
                           cs13_veg%livecrootc_SASUsave_patch(p)         = cs13_veg%livecrootc_SASUsave_patch(p)         + cs13_veg%livecrootc_patch(p)
                           cs13_veg%livecrootc_storage_SASUsave_patch(p) = cs13_veg%livecrootc_storage_SASUsave_patch(p) + cs13_veg%livecrootc_storage_patch(p)
                           cs13_veg%livecrootc_xfer_SASUsave_patch(p)    = cs13_veg%livecrootc_xfer_SASUsave_patch(p)    + cs13_veg%livecrootc_xfer_patch(p)
                           cs13_veg%deadcrootc_SASUsave_patch(p)         = cs13_veg%deadcrootc_SASUsave_patch(p)         + cs13_veg%deadcrootc_patch(p)
                           cs13_veg%deadcrootc_storage_SASUsave_patch(p) = cs13_veg%deadcrootc_storage_SASUsave_patch(p) + cs13_veg%deadcrootc_storage_patch(p)
                           cs13_veg%deadcrootc_xfer_SASUsave_patch(p)    = cs13_veg%deadcrootc_xfer_SASUsave_patch(p)    + cs13_veg%deadcrootc_xfer_patch(p)
                           if(ivt(p)  >= npcropmin)then
                              cs13_veg%grainc_SASUsave_patch(p)          = cs13_veg%grainc_SASUsave_patch(p)             + cs13_veg%reproductivec_patch(p,irepr)
                              cs13_veg%grainc_storage_SASUsave_patch(p)  = cs13_veg%grainc_storage_SASUsave_patch(p)     + cs13_veg%reproductivec_storage_patch(p,irepr)
                           end if
                        end if
                        if(use_c14)then
                           cs14_veg%leafc_SASUsave_patch(p)              = cs14_veg%leafc_SASUsave_patch(p)              + cs14_veg%leafc_patch(p)
                           cs14_veg%leafc_storage_SASUsave_patch(p)      = cs14_veg%leafc_storage_SASUsave_patch(p)      + cs14_veg%leafc_storage_patch(p)
                           cs14_veg%leafc_xfer_SASUsave_patch(p)         = cs14_veg%leafc_xfer_SASUsave_patch(p)         + cs14_veg%leafc_xfer_patch(p)
                           cs14_veg%frootc_SASUsave_patch(p)             = cs14_veg%frootc_SASUsave_patch(p)             + cs14_veg%frootc_patch(p)
                           cs14_veg%frootc_storage_SASUsave_patch(p)     = cs14_veg%frootc_storage_SASUsave_patch(p)     + cs14_veg%frootc_storage_patch(p)
                           cs14_veg%frootc_xfer_SASUsave_patch(p)        = cs14_veg%frootc_xfer_SASUsave_patch(p)        + cs14_veg%frootc_xfer_patch(p)
                           cs14_veg%livestemc_SASUsave_patch(p)          = cs14_veg%livestemc_SASUsave_patch(p)          + cs14_veg%livestemc_patch(p)
                           cs14_veg%livestemc_storage_SASUsave_patch(p)  = cs14_veg%livestemc_storage_SASUsave_patch(p)  + cs14_veg%livestemc_storage_patch(p)
                           cs14_veg%livestemc_xfer_SASUsave_patch(p)     = cs14_veg%livestemc_xfer_SASUsave_patch(p)     + cs14_veg%livestemc_xfer_patch(p)
                           cs14_veg%deadstemc_SASUsave_patch(p)          = cs14_veg%deadstemc_SASUsave_patch(p)          + cs14_veg%deadstemc_patch(p)
                           cs14_veg%deadstemc_storage_SASUsave_patch(p)  = cs14_veg%deadstemc_storage_SASUsave_patch(p)  + cs14_veg%deadstemc_storage_patch(p)
                           cs14_veg%deadstemc_xfer_SASUsave_patch(p)     = cs14_veg%deadstemc_xfer_SASUsave_patch(p)     + cs14_veg%deadstemc_xfer_patch(p)
                           cs14_veg%livecrootc_SASUsave_patch(p)         = cs14_veg%livecrootc_SASUsave_patch(p)         + cs14_veg%livecrootc_patch(p)
                           cs14_veg%livecrootc_storage_SASUsave_patch(p) = cs14_veg%livecrootc_storage_SASUsave_patch(p) + cs14_veg%livecrootc_storage_patch(p)
                           cs14_veg%livecrootc_xfer_SASUsave_patch(p)    = cs14_veg%livecrootc_xfer_SASUsave_patch(p)    + cs14_veg%livecrootc_xfer_patch(p)
                           cs14_veg%deadcrootc_SASUsave_patch(p)         = cs14_veg%deadcrootc_SASUsave_patch(p)         + cs14_veg%deadcrootc_patch(p)
                           cs14_veg%deadcrootc_storage_SASUsave_patch(p) = cs14_veg%deadcrootc_storage_SASUsave_patch(p) + cs14_veg%deadcrootc_storage_patch(p)
                           cs14_veg%deadcrootc_xfer_SASUsave_patch(p)    = cs14_veg%deadcrootc_xfer_SASUsave_patch(p)    + cs14_veg%deadcrootc_xfer_patch(p)
                           if(ivt(p)  >= npcropmin)then
                              cs14_veg%grainc_SASUsave_patch(p)          = cs14_veg%grainc_SASUsave_patch(p)             + cs14_veg%reproductivec_patch(p,irepr)
                              cs14_veg%grainc_storage_SASUsave_patch(p)  = cs14_veg%grainc_storage_SASUsave_patch(p)     + cs14_veg%reproductivec_storage_patch(p,irepr)
                           end if
                        end if
                        leafn_SASUsave(p)              = leafn_SASUsave(p)              + leafn(p)
                        leafn_storage_SASUsave(p)      = leafn_storage_SASUsave(p)      + leafn_storage(p)
                        leafn_xfer_SASUsave(p)         = leafn_xfer_SASUsave(p)         + leafn_xfer(p)
                        frootn_SASUsave(p)             = frootn_SASUsave(p)             + frootn(p)
                        frootn_storage_SASUsave(p)     = frootn_storage_SASUsave(p)     + frootn_storage(p)
                        frootn_xfer_SASUsave(p)        = frootn_xfer_SASUsave(p)        + frootn_xfer(p)
                        livestemn_SASUsave(p)          = livestemn_SASUsave(p)          + livestemn(p)
                        livestemn_storage_SASUsave(p)  = livestemn_storage_SASUsave(p)  + livestemn_storage(p)
                        livestemn_xfer_SASUsave(p)     = livestemn_xfer_SASUsave(p)     + livestemn_xfer(p)
                        deadstemn_SASUsave(p)          = deadstemn_SASUsave(p)          + deadstemn(p)
                        deadstemn_storage_SASUsave(p)  = deadstemn_storage_SASUsave(p)  + deadstemn_storage(p)
                        deadstemn_xfer_SASUsave(p)     = deadstemn_xfer_SASUsave(p)     + deadstemn_xfer(p)
                        livecrootn_SASUsave(p)         = livecrootn_SASUsave(p)         + livecrootn(p)
                        livecrootn_storage_SASUsave(p) = livecrootn_storage_SASUsave(p) + livecrootn_storage(p)
                        livecrootn_xfer_SASUsave(p)    = livecrootn_xfer_SASUsave(p)    + livecrootn_xfer(p)
                        deadcrootn_SASUsave(p)         = deadcrootn_SASUsave(p)         + deadcrootn(p)
                        deadcrootn_storage_SASUsave(p) = deadcrootn_storage_SASUsave(p) + deadcrootn_storage(p)
                        deadcrootn_xfer_SASUsave(p)    = deadcrootn_xfer_SASUsave(p)    + deadcrootn_xfer(p)
                        if(ivt(p)  >= npcropmin)then
                           grainn_SASUsave(p)          = grainn_SASUsave(p)             + reproductiven(p,irepr)
                        end if
                        if(iyr .eq. nyr_forcing)then
                           leafc(p)              = leafc_SASUsave(p)                    / (nyr_forcing/nyr_SASU)
                           leafc_storage(p)      = leafc_storage_SASUsave(p)            / (nyr_forcing/nyr_SASU)
                           leafc_xfer(p)         = leafc_xfer_SASUsave(p)               / (nyr_forcing/nyr_SASU)
                           frootc(p)             = frootc_SASUsave(p)                   / (nyr_forcing/nyr_SASU)
                           frootc_storage(p)     = frootc_storage_SASUsave(p)           / (nyr_forcing/nyr_SASU)
                           frootc_xfer(p)        = frootc_xfer_SASUsave(p)              / (nyr_forcing/nyr_SASU)
                           livestemc(p)          = livestemc_SASUsave(p)                / (nyr_forcing/nyr_SASU)
                           livestemc_storage(p)  = livestemc_storage_SASUsave(p)        / (nyr_forcing/nyr_SASU)
                           livestemc_xfer(p)     = livestemc_xfer_SASUsave(p)           / (nyr_forcing/nyr_SASU)
                           deadstemc(p)          = deadstemc_SASUsave(p)                / (nyr_forcing/nyr_SASU)
                           deadstemc_storage(p)  = deadstemc_storage_SASUsave(p)        / (nyr_forcing/nyr_SASU)
                           deadstemc_xfer(p)     = deadstemc_xfer_SASUsave(p)           / (nyr_forcing/nyr_SASU)
                           livecrootc(p)         = livecrootc_SASUsave(p)               / (nyr_forcing/nyr_SASU)
                           livecrootc_storage(p) = livecrootc_storage_SASUsave(p)       / (nyr_forcing/nyr_SASU)
                           livecrootc_xfer(p)    = livecrootc_xfer_SASUsave(p)          / (nyr_forcing/nyr_SASU)
                           deadcrootc(p)         = deadcrootc_SASUsave(p)               / (nyr_forcing/nyr_SASU)
                           deadcrootc_storage(p) = deadcrootc_storage_SASUsave(p)       / (nyr_forcing/nyr_SASU)
                           deadcrootc_xfer(p)    = deadcrootc_xfer_SASUsave(p)          / (nyr_forcing/nyr_SASU)
                           if(ivt(p)  >= npcropmin)then
                              reproductivec(p,:)          = grainc_SASUsave(p)                   / (nyr_forcing/nyr_SASU)
                              reproductivec_storage(p,:)  = grainc_storage_SASUsave(p)           / (nyr_forcing/nyr_SASU)
                           end if
                           if(use_c13)then
                              cs13_veg%leafc_patch(p)              = cs13_veg%leafc_SASUsave_patch(p)                    / (nyr_forcing/nyr_SASU)
                              cs13_veg%leafc_storage_patch(p)      = cs13_veg%leafc_storage_SASUsave_patch(p)            / (nyr_forcing/nyr_SASU)
                              cs13_veg%leafc_xfer_patch(p)         = cs13_veg%leafc_xfer_SASUsave_patch(p)               / (nyr_forcing/nyr_SASU)
                              cs13_veg%frootc_patch(p)             = cs13_veg%frootc_SASUsave_patch(p)                   / (nyr_forcing/nyr_SASU)
                              cs13_veg%frootc_storage_patch(p)     = cs13_veg%frootc_storage_SASUsave_patch(p)           / (nyr_forcing/nyr_SASU)
                              cs13_veg%frootc_xfer_patch(p)        = cs13_veg%frootc_xfer_SASUsave_patch(p)              / (nyr_forcing/nyr_SASU)
                              cs13_veg%livestemc_patch(p)          = cs13_veg%livestemc_SASUsave_patch(p)                / (nyr_forcing/nyr_SASU)
                              cs13_veg%livestemc_storage_patch(p)  = cs13_veg%livestemc_storage_SASUsave_patch(p)        / (nyr_forcing/nyr_SASU)
                              cs13_veg%livestemc_xfer_patch(p)     = cs13_veg%livestemc_xfer_SASUsave_patch(p)           / (nyr_forcing/nyr_SASU)
                              cs13_veg%deadstemc_patch(p)          = cs13_veg%deadstemc_SASUsave_patch(p)                / (nyr_forcing/nyr_SASU)
                              cs13_veg%deadstemc_storage_patch(p)  = cs13_veg%deadstemc_storage_SASUsave_patch(p)        / (nyr_forcing/nyr_SASU)
                              cs13_veg%deadstemc_xfer_patch(p)     = cs13_veg%deadstemc_xfer_SASUsave_patch(p)           / (nyr_forcing/nyr_SASU)
                              cs13_veg%livecrootc_patch(p)         = cs13_veg%livecrootc_SASUsave_patch(p)               / (nyr_forcing/nyr_SASU)
                              cs13_veg%livecrootc_storage_patch(p) = cs13_veg%livecrootc_storage_SASUsave_patch(p)       / (nyr_forcing/nyr_SASU)
                              cs13_veg%livecrootc_xfer_patch(p)    = cs13_veg%livecrootc_xfer_SASUsave_patch(p)          / (nyr_forcing/nyr_SASU)
                              cs13_veg%deadcrootc_patch(p)         = cs13_veg%deadcrootc_SASUsave_patch(p)               / (nyr_forcing/nyr_SASU)
                              cs13_veg%deadcrootc_storage_patch(p) = cs13_veg%deadcrootc_storage_SASUsave_patch(p)       / (nyr_forcing/nyr_SASU)
                              cs13_veg%deadcrootc_xfer_patch(p)    = cs13_veg%deadcrootc_xfer_SASUsave_patch(p)          / (nyr_forcing/nyr_SASU)
                              if(ivt(p)  >= npcropmin)then
                                 cs13_veg%reproductivec_patch(p,:)          = cs13_veg%grainc_SASUsave_patch(p)                   / (nyr_forcing/nyr_SASU)
                                 cs13_veg%reproductivec_storage_patch(p,:)  = cs13_veg%grainc_storage_SASUsave_patch(p)           / (nyr_forcing/nyr_SASU)
                              end if
                           end if
                           if(use_c14)then
                              cs14_veg%leafc_patch(p)              = cs14_veg%leafc_SASUsave_patch(p)                    / (nyr_forcing/nyr_SASU)
                              cs14_veg%leafc_storage_patch(p)      = cs14_veg%leafc_storage_SASUsave_patch(p)            / (nyr_forcing/nyr_SASU)
                              cs14_veg%leafc_xfer_patch(p)         = cs14_veg%leafc_xfer_SASUsave_patch(p)               / (nyr_forcing/nyr_SASU)
                              cs14_veg%frootc_patch(p)             = cs14_veg%frootc_SASUsave_patch(p)                   / (nyr_forcing/nyr_SASU)
                              cs14_veg%frootc_storage_patch(p)     = cs14_veg%frootc_storage_SASUsave_patch(p)           / (nyr_forcing/nyr_SASU)
                              cs14_veg%frootc_xfer_patch(p)        = cs14_veg%frootc_xfer_SASUsave_patch(p)              / (nyr_forcing/nyr_SASU)
                              cs14_veg%livestemc_patch(p)          = cs14_veg%livestemc_SASUsave_patch(p)                / (nyr_forcing/nyr_SASU)
                              cs14_veg%livestemc_storage_patch(p)  = cs14_veg%livestemc_storage_SASUsave_patch(p)        / (nyr_forcing/nyr_SASU)
                              cs14_veg%livestemc_xfer_patch(p)     = cs14_veg%livestemc_xfer_SASUsave_patch(p)           / (nyr_forcing/nyr_SASU)
                              cs14_veg%deadstemc_patch(p)          = cs14_veg%deadstemc_SASUsave_patch(p)                / (nyr_forcing/nyr_SASU)
                              cs14_veg%deadstemc_storage_patch(p)  = cs14_veg%deadstemc_storage_SASUsave_patch(p)        / (nyr_forcing/nyr_SASU)
                              cs14_veg%deadstemc_xfer_patch(p)     = cs14_veg%deadstemc_xfer_SASUsave_patch(p)           / (nyr_forcing/nyr_SASU)
                              cs14_veg%livecrootc_patch(p)         = cs14_veg%livecrootc_SASUsave_patch(p)               / (nyr_forcing/nyr_SASU)
                              cs14_veg%livecrootc_storage_patch(p) = cs14_veg%livecrootc_storage_SASUsave_patch(p)       / (nyr_forcing/nyr_SASU)
                              cs14_veg%livecrootc_xfer_patch(p)    = cs14_veg%livecrootc_xfer_SASUsave_patch(p)          / (nyr_forcing/nyr_SASU)
                              cs14_veg%deadcrootc_patch(p)         = cs14_veg%deadcrootc_SASUsave_patch(p)               / (nyr_forcing/nyr_SASU)
                              cs14_veg%deadcrootc_storage_patch(p) = cs14_veg%deadcrootc_storage_SASUsave_patch(p)       / (nyr_forcing/nyr_SASU)
                              cs14_veg%deadcrootc_xfer_patch(p)    = cs14_veg%deadcrootc_xfer_SASUsave_patch(p)          / (nyr_forcing/nyr_SASU)
                              if(ivt(p)  >= npcropmin)then
                                 cs14_veg%reproductivec_patch(p,:)          = cs14_veg%grainc_SASUsave_patch(p)                   / (nyr_forcing/nyr_SASU)
                                 cs14_veg%reproductivec_storage_patch(p,:)  = cs14_veg%grainc_storage_SASUsave_patch(p)           / (nyr_forcing/nyr_SASU)
                              end if
                           end if
                           leafn(p)              = leafn_SASUsave(p)                    / (nyr_forcing/nyr_SASU)
                           leafn_storage(p)      = leafn_storage_SASUsave(p)            / (nyr_forcing/nyr_SASU)
                           leafn_xfer(p)         = leafn_xfer_SASUsave(p)               / (nyr_forcing/nyr_SASU)
                           frootn(p)             = frootn_SASUsave(p)                   / (nyr_forcing/nyr_SASU)
                           frootn_storage(p)     = frootn_storage_SASUsave(p)           / (nyr_forcing/nyr_SASU)
                           frootn_xfer(p)        = frootn_xfer_SASUsave(p)              / (nyr_forcing/nyr_SASU)
                           livestemn(p)          = livestemn_SASUsave(p)                / (nyr_forcing/nyr_SASU)
                           livestemn_storage(p)  = livestemn_storage_SASUsave(p)        / (nyr_forcing/nyr_SASU)
                           livestemn_xfer(p)     = livestemn_xfer_SASUsave(p)           / (nyr_forcing/nyr_SASU)
                           deadstemn(p)          = deadstemn_SASUsave(p)                / (nyr_forcing/nyr_SASU)
                           deadstemn_storage(p)  = deadstemn_storage_SASUsave(p)        / (nyr_forcing/nyr_SASU)
                           deadstemn_xfer(p)     = deadstemn_xfer_SASUsave(p)           / (nyr_forcing/nyr_SASU)
                           livecrootn(p)         = livecrootn_SASUsave(p)               / (nyr_forcing/nyr_SASU)
                           livecrootn_storage(p) = livecrootn_storage_SASUsave(p)       / (nyr_forcing/nyr_SASU)
                           livecrootn_xfer(p)    = livecrootn_xfer_SASUsave(p)          / (nyr_forcing/nyr_SASU)
                           deadcrootn(p)         = deadcrootn_SASUsave(p)               / (nyr_forcing/nyr_SASU)
                           deadcrootn_storage(p) = deadcrootn_storage_SASUsave(p)       / (nyr_forcing/nyr_SASU)
                           deadcrootn_xfer(p)    = deadcrootn_xfer_SASUsave(p)          / (nyr_forcing/nyr_SASU)
                           if(ivt(p)  >= npcropmin)then
                              reproductiven(p,:)         = grainn_SASUsave(p)                   / (nyr_forcing/nyr_SASU)
                           end if
                           leafc_SASUsave(p)              = 0
                           leafc_storage_SASUsave(p)      = 0
                           leafc_xfer_SASUsave(p)         = 0
                           frootc_SASUsave(p)             = 0
                           frootc_storage_SASUsave(p)     = 0
                           frootc_xfer_SASUsave(p)        = 0
                           livestemc_SASUsave(p)          = 0
                           livestemc_storage_SASUsave(p)  = 0
                           livestemc_xfer_SASUsave(p)     = 0
                           deadstemc_SASUsave(p)          = 0
                           deadstemc_storage_SASUsave(p)  = 0
                           deadstemc_xfer_SASUsave(p)     = 0
                           livecrootc_SASUsave(p)         = 0
                           livecrootc_storage_SASUsave(p) = 0
                           livecrootc_xfer_SASUsave(p)    = 0
                           deadcrootc_SASUsave(p)         = 0
                           deadcrootc_storage_SASUsave(p) = 0
                           deadcrootc_xfer_SASUsave(p)    = 0
                           if(ivt(p)  >= npcropmin)then
                              grainc_SASUsave(p)          = 0
                              grainc_storage_SASUsave(p)  = 0
                           end if
                           if(use_c13)then
                              cs13_veg%leafc_SASUsave_patch(p)              = 0
                              cs13_veg%leafc_storage_SASUsave_patch(p)      = 0
                              cs13_veg%leafc_xfer_SASUsave_patch(p)         = 0
                              cs13_veg%frootc_SASUsave_patch(p)             = 0
                              cs13_veg%frootc_storage_SASUsave_patch(p)     = 0
                              cs13_veg%frootc_xfer_SASUsave_patch(p)        = 0
                              cs13_veg%livestemc_SASUsave_patch(p)          = 0
                              cs13_veg%livestemc_storage_SASUsave_patch(p)  = 0
                              cs13_veg%livestemc_xfer_SASUsave_patch(p)     = 0
                              cs13_veg%deadstemc_SASUsave_patch(p)          = 0
                              cs13_veg%deadstemc_storage_SASUsave_patch(p)  = 0
                              cs13_veg%deadstemc_xfer_SASUsave_patch(p)     = 0
                              cs13_veg%livecrootc_SASUsave_patch(p)         = 0
                              cs13_veg%livecrootc_storage_SASUsave_patch(p) = 0
                              cs13_veg%livecrootc_xfer_SASUsave_patch(p)    = 0
                              cs13_veg%deadcrootc_SASUsave_patch(p)         = 0
                              cs13_veg%deadcrootc_storage_SASUsave_patch(p) = 0
                              cs13_veg%deadcrootc_xfer_SASUsave_patch(p)    = 0
                              if(ivt(p)  >= npcropmin)then
                                 cs13_veg%grainc_SASUsave_patch(p)          = 0
                                 cs13_veg%grainc_storage_SASUsave_patch(p)  = 0
                              end if
                           end if
                           if(use_c14)then
                              cs14_veg%leafc_SASUsave_patch(p)              = 0
                              cs14_veg%leafc_storage_SASUsave_patch(p)      = 0
                              cs14_veg%leafc_xfer_SASUsave_patch(p)         = 0
                              cs14_veg%frootc_SASUsave_patch(p)             = 0
                              cs14_veg%frootc_storage_SASUsave_patch(p)     = 0
                              cs14_veg%frootc_xfer_SASUsave_patch(p)        = 0
                              cs14_veg%livestemc_SASUsave_patch(p)          = 0
                              cs14_veg%livestemc_storage_SASUsave_patch(p)  = 0
                              cs14_veg%livestemc_xfer_SASUsave_patch(p)     = 0
                              cs14_veg%deadstemc_SASUsave_patch(p)          = 0
                              cs14_veg%deadstemc_storage_SASUsave_patch(p)  = 0
                              cs14_veg%deadstemc_xfer_SASUsave_patch(p)     = 0
                              cs14_veg%livecrootc_SASUsave_patch(p)         = 0
                              cs14_veg%livecrootc_storage_SASUsave_patch(p) = 0
                              cs14_veg%livecrootc_xfer_SASUsave_patch(p)    = 0
                              cs14_veg%deadcrootc_SASUsave_patch(p)         = 0
                              cs14_veg%deadcrootc_storage_SASUsave_patch(p) = 0
                              cs14_veg%deadcrootc_xfer_SASUsave_patch(p)    = 0
                              if(ivt(p)  >= npcropmin)then
                                 cs14_veg%grainc_SASUsave_patch(p)          = 0
                                 cs14_veg%grainc_storage_SASUsave_patch(p)  = 0
                              end if
                           end if
                           leafn_SASUsave(p)              = 0
                           leafn_storage_SASUsave(p)      = 0
                           leafn_xfer_SASUsave(p)         = 0
                           frootn_SASUsave(p)             = 0
                           frootn_storage_SASUsave(p)     = 0
                           frootn_xfer_SASUsave(p)        = 0
                           livestemn_SASUsave(p)          = 0
                           livestemn_storage_SASUsave(p)  = 0
                           livestemn_xfer_SASUsave(p)     = 0
                           deadstemn_SASUsave(p)          = 0
                           deadstemn_storage_SASUsave(p)  = 0
                           deadstemn_xfer_SASUsave(p)     = 0
                           livecrootn_SASUsave(p)         = 0
                           livecrootn_storage_SASUsave(p) = 0
                           livecrootn_xfer_SASUsave(p)    = 0
                           deadcrootn_SASUsave(p)         = 0
                           deadcrootn_storage_SASUsave(p) = 0
                           deadcrootn_xfer_SASUsave(p)    = 0
                           if(ivt(p)  >= npcropmin)then
                              grainn_SASUsave(p)          = 0
                           end if
                        end if
                     end if
                     call update_DA_nstep()
                  end if

  ! Save C storage capacity from temporary variables to module variables
                  if(hist_wrt_matrixcn_diag)then
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
                        matrix_cap_reproc(p)              = vegmatrixc_rt(igrain)
                        matrix_cap_reproc_storage(p)      = vegmatrixc_rt(igrain_st)
                        matrix_cap_reproc_xfer(p)         = vegmatrixc_rt(igrain_xf)
                     end if
                     if(use_c13)then
                        cs13_veg%matrix_cap_leafc_patch(p)                  = vegmatrixc13_rt(ileaf)
                        cs13_veg%matrix_cap_leafc_storage_patch(p)          = vegmatrixc13_rt(ileaf_st)
                        cs13_veg%matrix_cap_leafc_xfer_patch(p)             = vegmatrixc13_rt(ileaf_xf)
                        cs13_veg%matrix_cap_frootc_patch(p)                 = vegmatrixc13_rt(ifroot)
                        cs13_veg%matrix_cap_frootc_storage_patch(p)         = vegmatrixc13_rt(ifroot_st)
                        cs13_veg%matrix_cap_frootc_xfer_patch(p)            = vegmatrixc13_rt(ifroot_xf)
                        cs13_veg%matrix_cap_livestemc_patch(p)              = vegmatrixc13_rt(ilivestem)
                        cs13_veg%matrix_cap_livestemc_storage_patch(p)      = vegmatrixc13_rt(ilivestem_st)
                        cs13_veg%matrix_cap_livestemc_xfer_patch(p)         = vegmatrixc13_rt(ilivestem_xf)
                        cs13_veg%matrix_cap_deadstemc_patch(p)              = vegmatrixc13_rt(ideadstem)
                        cs13_veg%matrix_cap_deadstemc_storage_patch(p)      = vegmatrixc13_rt(ideadstem_st)
                        cs13_veg%matrix_cap_deadstemc_xfer_patch(p)         = vegmatrixc13_rt(ideadstem_xf)
                        cs13_veg%matrix_cap_livecrootc_patch(p)             = vegmatrixc13_rt(ilivecroot)
                        cs13_veg%matrix_cap_livecrootc_storage_patch(p)     = vegmatrixc13_rt(ilivecroot_st)
                        cs13_veg%matrix_cap_livecrootc_xfer_patch(p)        = vegmatrixc13_rt(ilivecroot_xf)   
                        cs13_veg%matrix_cap_deadcrootc_patch(p)             = vegmatrixc13_rt(ideadcroot)
                        cs13_veg%matrix_cap_deadcrootc_storage_patch(p)     = vegmatrixc13_rt(ideadcroot_st)
                        cs13_veg%matrix_cap_deadcrootc_xfer_patch(p)        = vegmatrixc13_rt(ideadcroot_xf) 
                        if(ivt(p) >= npcropmin)then
                           cs13_veg%matrix_cap_reproc_patch(p)              = vegmatrixc13_rt(igrain)
                           cs13_veg%matrix_cap_reproc_storage_patch(p)      = vegmatrixc13_rt(igrain_st)
                           cs13_veg%matrix_cap_reproc_xfer_patch(p)         = vegmatrixc13_rt(igrain_xf)
                        end if
                     end if
                     if(use_c14)then
                        cs14_veg%matrix_cap_leafc_patch(p)                  = vegmatrixc14_rt(ileaf)
                        cs14_veg%matrix_cap_leafc_storage_patch(p)          = vegmatrixc14_rt(ileaf_st)
                        cs14_veg%matrix_cap_leafc_xfer_patch(p)             = vegmatrixc14_rt(ileaf_xf)
                        cs14_veg%matrix_cap_frootc_patch(p)                 = vegmatrixc14_rt(ifroot)
                        cs14_veg%matrix_cap_frootc_storage_patch(p)         = vegmatrixc14_rt(ifroot_st)
                        cs14_veg%matrix_cap_frootc_xfer_patch(p)            = vegmatrixc14_rt(ifroot_xf)
                        cs14_veg%matrix_cap_livestemc_patch(p)              = vegmatrixc14_rt(ilivestem)
                        cs14_veg%matrix_cap_livestemc_storage_patch(p)      = vegmatrixc14_rt(ilivestem_st)
                        cs14_veg%matrix_cap_livestemc_xfer_patch(p)         = vegmatrixc14_rt(ilivestem_xf)
                        cs14_veg%matrix_cap_deadstemc_patch(p)              = vegmatrixc14_rt(ideadstem)
                        cs14_veg%matrix_cap_deadstemc_storage_patch(p)      = vegmatrixc14_rt(ideadstem_st)
                        cs14_veg%matrix_cap_deadstemc_xfer_patch(p)         = vegmatrixc14_rt(ideadstem_xf)
                        cs14_veg%matrix_cap_livecrootc_patch(p)             = vegmatrixc14_rt(ilivecroot)
                        cs14_veg%matrix_cap_livecrootc_storage_patch(p)     = vegmatrixc14_rt(ilivecroot_st)
                        cs14_veg%matrix_cap_livecrootc_xfer_patch(p)        = vegmatrixc14_rt(ilivecroot_xf)   
                        cs14_veg%matrix_cap_deadcrootc_patch(p)             = vegmatrixc14_rt(ideadcroot)
                        cs14_veg%matrix_cap_deadcrootc_storage_patch(p)     = vegmatrixc14_rt(ideadcroot_st)
                        cs14_veg%matrix_cap_deadcrootc_xfer_patch(p)        = vegmatrixc14_rt(ideadcroot_xf) 
                        if(ivt(p) >= npcropmin)then
                           cs14_veg%matrix_cap_reproc_patch(p)              = vegmatrixc14_rt(igrain)
                           cs14_veg%matrix_cap_reproc_storage_patch(p)      = vegmatrixc14_rt(igrain_st)
                           cs14_veg%matrix_cap_reproc_xfer_patch(p)         = vegmatrixc14_rt(igrain_xf)
                        end if
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
                        matrix_cap_repron(p)              = vegmatrixn_rt(igrain)
                        matrix_cap_repron_storage(p)      = vegmatrixn_rt(igrain_st)
                        matrix_cap_repron_xfer(p)         = vegmatrixn_rt(igrain_xf)
                     end if
                  end if
                    
  ! Reset accumulated variables to 0 at end of each year after calculating capacity 
                  matrix_calloc_leaf_acc(p)                = 0._r8
                  matrix_calloc_leafst_acc(p)              = 0._r8 
                  matrix_calloc_froot_acc(p)               = 0._r8 
                  matrix_calloc_frootst_acc(p)             = 0._r8
                  matrix_calloc_livestem_acc(p)            = 0._r8
                  matrix_calloc_livestemst_acc(p)          = 0._r8
                  matrix_calloc_deadstem_acc(p)            = 0._r8
                  matrix_calloc_deadstemst_acc(p)          = 0._r8
                  matrix_calloc_livecroot_acc(p)           = 0._r8
                  matrix_calloc_livecrootst_acc(p)         = 0._r8
                  matrix_calloc_deadcroot_acc(p)           = 0._r8
                  matrix_calloc_deadcrootst_acc(p)         = 0._r8
                  if(ivt(p) >= npcropmin)then
                     matrix_calloc_grain_acc(p)            = 0._r8
                     matrix_calloc_grainst_acc(p)          = 0._r8
                  end if
   
                  matrix_ctransfer_leafst_to_leafxf_acc(p)           = 0._r8
                  matrix_ctransfer_leafxf_to_leaf_acc(p)             = 0._r8
                  matrix_ctransfer_frootst_to_frootxf_acc(p)         = 0._r8
                  matrix_ctransfer_frootxf_to_froot_acc(p)           = 0._r8
                  matrix_ctransfer_livestemst_to_livestemxf_acc(p)   = 0._r8
                  matrix_ctransfer_livestemxf_to_livestem_acc(p)     = 0._r8
                  matrix_ctransfer_deadstemst_to_deadstemxf_acc(p)   = 0._r8
                  matrix_ctransfer_deadstemxf_to_deadstem_acc(p)     = 0._r8
                  matrix_ctransfer_livecrootst_to_livecrootxf_acc(p) = 0._r8
                  matrix_ctransfer_livecrootxf_to_livecroot_acc(p)   = 0._r8
                  matrix_ctransfer_deadcrootst_to_deadcrootxf_acc(p) = 0._r8
                  matrix_ctransfer_deadcrootxf_to_deadcroot_acc(p)   = 0._r8
                  if(ivt(p) >= npcropmin)then
                     matrix_ctransfer_grainst_to_grainxf_acc(p)      = 0._r8
                     matrix_ctransfer_grainxf_to_grain_acc(p)        = 0._r8
                  end if
                  matrix_ctransfer_livestem_to_deadstem_acc(p)       = 0._r8
                  matrix_ctransfer_livecroot_to_deadcroot_acc(p)     = 0._r8
   
                  matrix_cturnover_leaf_acc(p)                       = 0._r8
                  matrix_cturnover_leafst_acc(p)                     = 0._r8
                  matrix_cturnover_leafxf_acc(p)                     = 0._r8
                  matrix_cturnover_froot_acc(p)                      = 0._r8
                  matrix_cturnover_frootst_acc(p)                    = 0._r8
                  matrix_cturnover_frootxf_acc(p)                    = 0._r8
                  matrix_cturnover_livestem_acc(p)                   = 0._r8
                  matrix_cturnover_livestemst_acc(p)                 = 0._r8
                  matrix_cturnover_livestemxf_acc(p)                 = 0._r8
                  matrix_cturnover_deadstem_acc(p)                   = 0._r8
                  matrix_cturnover_deadstemst_acc(p)                 = 0._r8
                  matrix_cturnover_deadstemxf_acc(p)                 = 0._r8
                  matrix_cturnover_livecroot_acc(p)                  = 0._r8
                  matrix_cturnover_livecrootst_acc(p)                = 0._r8
                  matrix_cturnover_livecrootxf_acc(p)                = 0._r8
                  matrix_cturnover_deadcroot_acc(p)                  = 0._r8
                  matrix_cturnover_deadcrootst_acc(p)                = 0._r8
                  matrix_cturnover_deadcrootxf_acc(p)                = 0._r8
                  if(ivt(p) >= npcropmin)then
                     matrix_cturnover_grain_acc(p)                   = 0._r8 
                     matrix_cturnover_grainst_acc(p)                 = 0._r8
                     matrix_cturnover_grainxf_acc(p)                 = 0._r8
                  end if

                  if(use_c13)then
                     cs13_veg%matrix_calloc_leaf_acc_patch(p)                = 0._r8
                     cs13_veg%matrix_calloc_leafst_acc_patch(p)              = 0._r8 
                     cs13_veg%matrix_calloc_froot_acc_patch(p)               = 0._r8 
                     cs13_veg%matrix_calloc_frootst_acc_patch(p)             = 0._r8
                     cs13_veg%matrix_calloc_livestem_acc_patch(p)            = 0._r8
                     cs13_veg%matrix_calloc_livestemst_acc_patch(p)          = 0._r8
                     cs13_veg%matrix_calloc_deadstem_acc_patch(p)            = 0._r8
                     cs13_veg%matrix_calloc_deadstemst_acc_patch(p)          = 0._r8
                     cs13_veg%matrix_calloc_livecroot_acc_patch(p)           = 0._r8
                     cs13_veg%matrix_calloc_livecrootst_acc_patch(p)         = 0._r8
                     cs13_veg%matrix_calloc_deadcroot_acc_patch(p)           = 0._r8
                     cs13_veg%matrix_calloc_deadcrootst_acc_patch(p)         = 0._r8
                     if(ivt(p) >= npcropmin)then
                        cs13_veg%matrix_calloc_grain_acc_patch(p)            = 0._r8
                        cs13_veg%matrix_calloc_grainst_acc_patch(p)          = 0._r8
                     end if
   
                     cs13_veg%matrix_ctransfer_leafst_to_leafxf_acc_patch(p)           = 0._r8
                     cs13_veg%matrix_ctransfer_leafxf_to_leaf_acc_patch(p)             = 0._r8
                     cs13_veg%matrix_ctransfer_frootst_to_frootxf_acc_patch(p)         = 0._r8
                     cs13_veg%matrix_ctransfer_frootxf_to_froot_acc_patch(p)           = 0._r8
                     cs13_veg%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p)   = 0._r8
                     cs13_veg%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p)     = 0._r8
                     cs13_veg%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p)   = 0._r8
                     cs13_veg%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p)     = 0._r8
                     cs13_veg%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p) = 0._r8
                     cs13_veg%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p)   = 0._r8
                     cs13_veg%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p) = 0._r8
                     cs13_veg%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p)   = 0._r8
                     if(ivt(p) >= npcropmin)then
                        cs13_veg%matrix_ctransfer_grainst_to_grainxf_acc_patch(p)      = 0._r8
                        cs13_veg%matrix_ctransfer_grainxf_to_grain_acc_patch(p)        = 0._r8
                     end if
                     cs13_veg%matrix_ctransfer_livestem_to_deadstem_acc_patch(p)       = 0._r8
                     cs13_veg%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p)     = 0._r8
      
                     cs13_veg%matrix_cturnover_leaf_acc_patch(p)                       = 0._r8
                     cs13_veg%matrix_cturnover_leafst_acc_patch(p)                     = 0._r8
                     cs13_veg%matrix_cturnover_leafxf_acc_patch(p)                     = 0._r8
                     cs13_veg%matrix_cturnover_froot_acc_patch(p)                      = 0._r8
                     cs13_veg%matrix_cturnover_frootst_acc_patch(p)                    = 0._r8
                     cs13_veg%matrix_cturnover_frootxf_acc_patch(p)                    = 0._r8
                     cs13_veg%matrix_cturnover_livestem_acc_patch(p)                   = 0._r8
                     cs13_veg%matrix_cturnover_livestemst_acc_patch(p)                 = 0._r8
                     cs13_veg%matrix_cturnover_livestemxf_acc_patch(p)                 = 0._r8
                     cs13_veg%matrix_cturnover_deadstem_acc_patch(p)                   = 0._r8
                     cs13_veg%matrix_cturnover_deadstemst_acc_patch(p)                 = 0._r8
                     cs13_veg%matrix_cturnover_deadstemxf_acc_patch(p)                 = 0._r8
                     cs13_veg%matrix_cturnover_livecroot_acc_patch(p)                  = 0._r8
                     cs13_veg%matrix_cturnover_livecrootst_acc_patch(p)                = 0._r8
                     cs13_veg%matrix_cturnover_livecrootxf_acc_patch(p)                = 0._r8
                     cs13_veg%matrix_cturnover_deadcroot_acc_patch(p)                  = 0._r8
                     cs13_veg%matrix_cturnover_deadcrootst_acc_patch(p)                = 0._r8
                     cs13_veg%matrix_cturnover_deadcrootxf_acc_patch(p)                = 0._r8
                     if(ivt(p) >= npcropmin)then
                        cs13_veg%matrix_cturnover_grain_acc_patch(p)                   = 0._r8 
                        cs13_veg%matrix_cturnover_grainst_acc_patch(p)                 = 0._r8
                        cs13_veg%matrix_cturnover_grainxf_acc_patch(p)                 = 0._r8
                     end if
                  end if

                  if(use_c14)then
                     cs14_veg%matrix_calloc_leaf_acc_patch(p)                = 0._r8
                     cs14_veg%matrix_calloc_leafst_acc_patch(p)              = 0._r8 
                     cs14_veg%matrix_calloc_froot_acc_patch(p)               = 0._r8 
                     cs14_veg%matrix_calloc_frootst_acc_patch(p)             = 0._r8
                     cs14_veg%matrix_calloc_livestem_acc_patch(p)            = 0._r8
                     cs14_veg%matrix_calloc_livestemst_acc_patch(p)          = 0._r8
                     cs14_veg%matrix_calloc_deadstem_acc_patch(p)            = 0._r8
                     cs14_veg%matrix_calloc_deadstemst_acc_patch(p)          = 0._r8
                     cs14_veg%matrix_calloc_livecroot_acc_patch(p)           = 0._r8
                     cs14_veg%matrix_calloc_livecrootst_acc_patch(p)         = 0._r8
                     cs14_veg%matrix_calloc_deadcroot_acc_patch(p)           = 0._r8
                     cs14_veg%matrix_calloc_deadcrootst_acc_patch(p)         = 0._r8
                     if(ivt(p) >= npcropmin)then
                        cs14_veg%matrix_calloc_grain_acc_patch(p)            = 0._r8
                        cs14_veg%matrix_calloc_grainst_acc_patch(p)          = 0._r8
                     end if
   
                     cs14_veg%matrix_ctransfer_leafst_to_leafxf_acc_patch(p)           = 0._r8
                     cs14_veg%matrix_ctransfer_leafxf_to_leaf_acc_patch(p)             = 0._r8
                     cs14_veg%matrix_ctransfer_frootst_to_frootxf_acc_patch(p)         = 0._r8
                     cs14_veg%matrix_ctransfer_frootxf_to_froot_acc_patch(p)           = 0._r8
                     cs14_veg%matrix_ctransfer_livestemst_to_livestemxf_acc_patch(p)   = 0._r8
                     cs14_veg%matrix_ctransfer_livestemxf_to_livestem_acc_patch(p)     = 0._r8
                     cs14_veg%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch(p)   = 0._r8
                     cs14_veg%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch(p)     = 0._r8
                     cs14_veg%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch(p) = 0._r8
                     cs14_veg%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch(p)   = 0._r8
                     cs14_veg%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch(p) = 0._r8
                     cs14_veg%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch(p)   = 0._r8
                     if(ivt(p) >= npcropmin)then
                        cs14_veg%matrix_ctransfer_grainst_to_grainxf_acc_patch(p)      = 0._r8
                        cs14_veg%matrix_ctransfer_grainxf_to_grain_acc_patch(p)        = 0._r8
                     end if
                     cs14_veg%matrix_ctransfer_livestem_to_deadstem_acc_patch(p)       = 0._r8
                     cs14_veg%matrix_ctransfer_livecroot_to_deadcroot_acc_patch(p)     = 0._r8
      
                     cs14_veg%matrix_cturnover_leaf_acc_patch(p)                       = 0._r8
                     cs14_veg%matrix_cturnover_leafst_acc_patch(p)                     = 0._r8
                     cs14_veg%matrix_cturnover_leafxf_acc_patch(p)                     = 0._r8
                     cs14_veg%matrix_cturnover_froot_acc_patch(p)                      = 0._r8
                     cs14_veg%matrix_cturnover_frootst_acc_patch(p)                    = 0._r8
                     cs14_veg%matrix_cturnover_frootxf_acc_patch(p)                    = 0._r8
                     cs14_veg%matrix_cturnover_livestem_acc_patch(p)                   = 0._r8
                     cs14_veg%matrix_cturnover_livestemst_acc_patch(p)                 = 0._r8
                     cs14_veg%matrix_cturnover_livestemxf_acc_patch(p)                 = 0._r8
                     cs14_veg%matrix_cturnover_deadstem_acc_patch(p)                   = 0._r8
                     cs14_veg%matrix_cturnover_deadstemst_acc_patch(p)                 = 0._r8
                     cs14_veg%matrix_cturnover_deadstemxf_acc_patch(p)                 = 0._r8
                     cs14_veg%matrix_cturnover_livecroot_acc_patch(p)                  = 0._r8
                     cs14_veg%matrix_cturnover_livecrootst_acc_patch(p)                = 0._r8
                     cs14_veg%matrix_cturnover_livecrootxf_acc_patch(p)                = 0._r8
                     cs14_veg%matrix_cturnover_deadcroot_acc_patch(p)                  = 0._r8
                     cs14_veg%matrix_cturnover_deadcrootst_acc_patch(p)                = 0._r8
                     cs14_veg%matrix_cturnover_deadcrootxf_acc_patch(p)                = 0._r8
                     if(ivt(p) >= npcropmin)then
                        cs14_veg%matrix_cturnover_grain_acc_patch(p)                   = 0._r8 
                        cs14_veg%matrix_cturnover_grainst_acc_patch(p)                 = 0._r8
                        cs14_veg%matrix_cturnover_grainxf_acc_patch(p)                 = 0._r8
                     end if
                  end if

                  matrix_nalloc_leaf_acc(p)                          = 0._r8
                  matrix_nalloc_leafst_acc(p)                        = 0._r8
                  matrix_nalloc_froot_acc(p)                         = 0._r8
                  matrix_nalloc_frootst_acc(p)                       = 0._r8
                  matrix_nalloc_livestem_acc(p)                      = 0._r8
                  matrix_nalloc_livestemst_acc(p)                    = 0._r8
                  matrix_nalloc_deadstem_acc(p)                      = 0._r8
                  matrix_nalloc_deadstemst_acc(p)                    = 0._r8
                  matrix_nalloc_livecroot_acc(p)                     = 0._r8
                  matrix_nalloc_livecrootst_acc(p)                   = 0._r8
                  matrix_nalloc_deadcroot_acc(p)                     = 0._r8
                  matrix_nalloc_deadcrootst_acc(p)                   = 0._r8
                  if(ivt(p) >= npcropmin)then
                     matrix_nalloc_grain_acc(p)                      = 0._r8
                     matrix_nalloc_grainst_acc(p)                    = 0._r8
                  end if
   
                  matrix_ntransfer_leafst_to_leafxf_acc(p)           = 0._r8
                  matrix_ntransfer_leafxf_to_leaf_acc(p)             = 0._r8
                  matrix_ntransfer_frootst_to_frootxf_acc(p)         = 0._r8
                  matrix_ntransfer_frootxf_to_froot_acc(p)           = 0._r8
                  matrix_ntransfer_livestemst_to_livestemxf_acc(p)   = 0._r8
                  matrix_ntransfer_livestemxf_to_livestem_acc(p)     = 0._r8
                  matrix_ntransfer_deadstemst_to_deadstemxf_acc(p)   = 0._r8
                  matrix_ntransfer_deadstemxf_to_deadstem_acc(p)     = 0._r8
                  matrix_ntransfer_livecrootst_to_livecrootxf_acc(p) = 0._r8
                  matrix_ntransfer_livecrootxf_to_livecroot_acc(p)   = 0._r8
                  matrix_ntransfer_deadcrootst_to_deadcrootxf_acc(p) = 0._r8
                  matrix_ntransfer_deadcrootxf_to_deadcroot_acc(p)   = 0._r8
                  if(ivt(p) >= npcropmin)then
                     matrix_ntransfer_grainst_to_grainxf_acc(p)      = 0._r8
                     matrix_ntransfer_grainxf_to_grain_acc(p)        = 0._r8
                  end if
                  matrix_ntransfer_livestem_to_deadstem_acc(p)       = 0._r8
                  matrix_ntransfer_livecroot_to_deadcroot_acc(p)     = 0._r8
   
                  matrix_ntransfer_retransn_to_leaf_acc(p)           = 0._r8 
                  matrix_ntransfer_retransn_to_leafst_acc(p)         = 0._r8
                  matrix_ntransfer_retransn_to_froot_acc(p)          = 0._r8
                  matrix_ntransfer_retransn_to_frootst_acc(p)        = 0._r8
                  matrix_ntransfer_retransn_to_livestem_acc(p)       = 0._r8
                  matrix_ntransfer_retransn_to_livestemst_acc(p)     = 0._r8
                  matrix_ntransfer_retransn_to_deadstem_acc(p)       = 0._r8
                  matrix_ntransfer_retransn_to_deadstemst_acc(p)     = 0._r8
                  matrix_ntransfer_retransn_to_livecroot_acc(p)      = 0._r8
                  matrix_ntransfer_retransn_to_livecrootst_acc(p)    = 0._r8
                  matrix_ntransfer_retransn_to_deadcroot_acc(p)      = 0._r8
                  matrix_ntransfer_retransn_to_deadcrootst_acc(p)    = 0._r8
                  if(ivt(p) >= npcropmin)then
                     matrix_ntransfer_retransn_to_grain_acc(p)       = 0._r8
                     matrix_ntransfer_retransn_to_grainst_acc(p)     = 0._r8
                  end if
                  matrix_ntransfer_leaf_to_retransn_acc(p)           = 0._r8
                  matrix_ntransfer_froot_to_retransn_acc(p)          = 0._r8
                  matrix_ntransfer_livestem_to_retransn_acc(p)       = 0._r8
                  matrix_ntransfer_livecroot_to_retransn_acc(p)      = 0._r8
   
                  matrix_nturnover_leaf_acc(p)                       = 0._r8
                  matrix_nturnover_leafst_acc(p)                     = 0._r8
                  matrix_nturnover_leafxf_acc(p)                     = 0._r8
                  matrix_nturnover_froot_acc(p)                      = 0._r8
                  matrix_nturnover_frootst_acc(p)                    = 0._r8
                  matrix_nturnover_frootxf_acc(p)                    = 0._r8
                  matrix_nturnover_livestem_acc(p)                   = 0._r8
                  matrix_nturnover_livestemst_acc(p)                 = 0._r8
                  matrix_nturnover_livestemxf_acc(p)                 = 0._r8
                  matrix_nturnover_deadstem_acc(p)                   = 0._r8
                  matrix_nturnover_deadstemst_acc(p)                 = 0._r8
                  matrix_nturnover_deadstemxf_acc(p)                 = 0._r8
                  matrix_nturnover_livecroot_acc(p)                  = 0._r8
                  matrix_nturnover_livecrootst_acc(p)                = 0._r8
                  matrix_nturnover_livecrootxf_acc(p)                = 0._r8
                  matrix_nturnover_deadcroot_acc(p)                  = 0._r8
                  matrix_nturnover_deadcrootst_acc(p)                = 0._r8
                  matrix_nturnover_deadcrootxf_acc(p)                = 0._r8
                  if(ivt(p) >= npcropmin)then
                     matrix_nturnover_grain_acc(p)                   = 0._r8
                     matrix_nturnover_grainst_acc(p)                 = 0._r8
                     matrix_nturnover_grainxf_acc(p)                 = 0._r8
                  end if
                  matrix_nturnover_retransn_acc(p)                   = 0._r8
                  matrix_calloc_acc(:) = 0._r8
                  matrix_ctransfer_acc(:,:) = 0._r8
                  matrix_nalloc_acc(:) = 0._r8
                  matrix_ntransfer_acc(:,:) = 0._r8

                  call t_stopf('CN veg matrix-finalize spinup')
               end do 
               if(iloop .eq. iloop_avg .and. iyr .eq. nyr_forcing)iloop = 0
               if(iyr .eq. nyr_forcing)iyr=0
            end if
         end if
   
      call vegmatrixc_input%ReleaseV()
      if ( use_c13 )then
         call vegmatrixc13_input%ReleaseV()
      end if
      if ( use_c14 )then
         call vegmatrixc14_input%ReleaseV()
      end if
      call vegmatrixn_input%ReleaseV()
    
   end associate td
   end associate sd
   end associate od
   end associate fr
 end subroutine CNVegMatrix

 function matrix_update_phc(p,itransfer,rate,dt,cnveg_carbonflux_inst,matrixcheck,acc)

   integer ,intent(in) :: p
   integer ,intent(in) :: itransfer
   real(r8),intent(in) :: rate
   real(r8),intent(in) :: dt
   type(cnveg_carbonflux_type)  , intent(inout) :: cnveg_carbonflux_inst 
   logical ,intent(in),optional :: matrixcheck
   logical ,intent(in),optional :: acc
   real(r8)            :: matrix_update_phc

   associate(                          &
     matrix_phtransfer                   => cnveg_carbonflux_inst%matrix_phtransfer_patch , &
     matrix_phturnover                   => cnveg_carbonflux_inst%matrix_phturnover_patch , &
     doner_phc                           => cnveg_carbonflux_inst%matrix_phtransfer_doner_patch&
   )
      if(.not. present(matrixcheck) .or. matrixcheck)then
         if((.not. present(acc) .or. acc) .and. matrix_phturnover(p,doner_phc(itransfer)) + rate * dt .ge. 1)then
            matrix_update_phc = max(0._r8,(1._r8 - matrix_phturnover(p,doner_phc(itransfer))) / dt)
         else
            matrix_update_phc = rate
         end if
      else
         matrix_update_phc = rate
      end if
      if(.not. present(acc) .or. acc)then
         matrix_phturnover(p,doner_phc(itransfer)) = matrix_phturnover(p,doner_phc(itransfer)) + matrix_update_phc * dt
         matrix_phtransfer(p,itransfer) = matrix_phtransfer(p,itransfer) + matrix_update_phc
      else
         matrix_phturnover(p,doner_phc(itransfer)) = matrix_phturnover(p,doner_phc(itransfer)) - matrix_phtransfer(p,itransfer) * dt + matrix_update_phc * dt
         matrix_phtransfer(p,itransfer) = matrix_update_phc
      end if

      return
   end associate

 end function matrix_update_phc

 function matrix_update_gmc(p,itransfer,rate,dt,cnveg_carbonflux_inst,matrixcheck,acc)

   integer,intent(in) :: p
   integer,intent(in) :: itransfer
   real(r8),intent(in) :: rate
   real(r8),intent(in) :: dt
   type(cnveg_carbonflux_type)  , intent(inout) :: cnveg_carbonflux_inst 
   logical ,intent(in),optional :: matrixcheck
   logical ,intent(in),optional :: acc
   real(r8)            :: matrix_update_gmc

   associate(                          &
   matrix_phturnover                   => cnveg_carbonflux_inst%matrix_phturnover_patch , &
   matrix_gmtransfer                   => cnveg_carbonflux_inst%matrix_gmtransfer_patch , &
   matrix_gmturnover                   => cnveg_carbonflux_inst%matrix_gmturnover_patch , &
   doner_gmc                           => cnveg_carbonflux_inst%matrix_gmtransfer_doner_patch  & ! Input:  [integer (:)] Doners of gap mortality related C transfer
   )

      if(.not. present(matrixcheck) .or. matrixcheck)then
         if((.not. present(acc) .or. acc) .and. matrix_phturnover(p,doner_gmc(itransfer)) + matrix_gmturnover(p,doner_gmc(itransfer)) + rate * dt .ge. 1)then
            matrix_update_gmc = max(0._r8,(1._r8 - matrix_phturnover(p,doner_gmc(itransfer)) - matrix_gmturnover(p,doner_gmc(itransfer))) / dt)
         else
            matrix_update_gmc = rate
         end if
      else
         matrix_update_gmc = rate
      end if
      if(.not. present(acc) .or. acc)then
         matrix_gmturnover(p,doner_gmc(itransfer)) = matrix_gmturnover(p,doner_gmc(itransfer)) + matrix_update_gmc * dt
         matrix_gmtransfer(p,itransfer) = matrix_gmtransfer(p,itransfer) + matrix_update_gmc
      else
         matrix_gmturnover(p,doner_gmc(itransfer)) = matrix_gmturnover(p,doner_gmc(itransfer)) - matrix_gmtransfer(p,itransfer) * dt + matrix_update_gmc * dt
         matrix_gmtransfer(p,itransfer) = matrix_update_gmc
      end if
      return
   end associate

 end function matrix_update_gmc


 function matrix_update_fic(p,itransfer,rate,dt,cnveg_carbonflux_inst,matrixcheck,acc)

   integer,intent(in) :: p
   integer,intent(in) :: itransfer
   real(r8),intent(in) :: rate
   real(r8),intent(in) :: dt
   type(cnveg_carbonflux_type)  , intent(inout) :: cnveg_carbonflux_inst 
   logical ,intent(in),optional :: matrixcheck
   logical ,intent(in),optional :: acc
   real(r8)            :: matrix_update_fic

   associate(                          &
   matrix_phturnover                   => cnveg_carbonflux_inst%matrix_phturnover_patch , &
   matrix_gmturnover                   => cnveg_carbonflux_inst%matrix_gmturnover_patch , &
   matrix_fitransfer                   => cnveg_carbonflux_inst%matrix_fitransfer_patch , &
   matrix_fiturnover                   => cnveg_carbonflux_inst%matrix_fiturnover_patch , &
   doner_fic                           => cnveg_carbonflux_inst%matrix_fitransfer_doner_patch &
   )

      if(.not. present(matrixcheck) .or. matrixcheck)then
         if((.not. present(acc) .or. acc) .and. matrix_phturnover(p,doner_fic(itransfer)) + matrix_gmturnover(p,doner_fic(itransfer)) &
           + matrix_fiturnover(p,doner_fic(itransfer)) + rate * dt .ge. 1)then
            matrix_update_fic = max(0._r8,(1._r8 - matrix_phturnover(p,doner_fic(itransfer)) &
                              - matrix_gmturnover(p,doner_fic(itransfer)) - matrix_fiturnover(p,doner_fic(itransfer))) / dt)
         else
            matrix_update_fic = rate
         end if
      else
         matrix_update_fic = rate
      end if
      if(.not. present(acc) .or. acc)then
         matrix_fiturnover(p,doner_fic(itransfer)) = matrix_fiturnover(p,doner_fic(itransfer)) + matrix_update_fic * dt
         matrix_fitransfer(p,itransfer) = matrix_fitransfer(p,itransfer) + matrix_update_fic
      else
         matrix_fiturnover(p,doner_fic(itransfer)) = matrix_fiturnover(p,doner_fic(itransfer)) - matrix_fitransfer(p,itransfer) * dt + matrix_update_fic * dt
         matrix_fitransfer(p,itransfer) = matrix_update_fic
      end if

      return
   end associate

end function matrix_update_fic

 function matrix_update_phn(p,itransfer,rate,dt,cnveg_nitrogenflux_inst,matrixcheck,acc)

   integer,intent(in) :: p
   integer,intent(in) :: itransfer
   real(r8),intent(in) :: rate
   real(r8),intent(in) :: dt
   type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst 
   logical ,intent(in),optional :: matrixcheck
   logical ,intent(in),optional :: acc
   real(r8)            :: matrix_update_phn

   associate(                          &
     matrix_nphtransfer                   => cnveg_nitrogenflux_inst%matrix_nphtransfer_patch , &
     matrix_nphturnover                   => cnveg_nitrogenflux_inst%matrix_nphturnover_patch , &
     doner_phn                            => cnveg_nitrogenflux_inst%matrix_nphtransfer_doner_patch & ! Input:  [integer (:)] Doners of phenology related N transfer
   )

      if(.not. present(matrixcheck) .or. matrixcheck)then
         if((.not. present(acc) .or. acc) .and. matrix_nphturnover(p,doner_phn(itransfer)) + rate * dt .ge. 1)then
            matrix_update_phn = max(0._r8,(1._r8 - matrix_nphturnover(p,doner_phn(itransfer))) / dt)
         else
            matrix_update_phn = rate
         end if
      else
         matrix_update_phn = rate
      end if
      if(.not. present(acc) .or. acc)then
         matrix_nphturnover(p,doner_phn(itransfer)) = matrix_nphturnover(p,doner_phn(itransfer)) + matrix_update_phn * dt
         matrix_nphtransfer(p,itransfer) = matrix_nphtransfer(p,itransfer) + matrix_update_phn
      else
         matrix_nphturnover(p,doner_phn(itransfer)) = matrix_nphturnover(p,doner_phn(itransfer)) - matrix_nphtransfer(p,itransfer) * dt + matrix_update_phn * dt
         matrix_nphtransfer(p,itransfer) = matrix_update_phn
      end if

      return
   end associate

 end function matrix_update_phn

 function matrix_update_gmn(p,itransfer,rate,dt,cnveg_nitrogenflux_inst,matrixcheck,acc)

   integer ,intent(in) :: p
   integer ,intent(in) :: itransfer
   real(r8),intent(in) :: rate
   real(r8),intent(in) :: dt
   type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst 
   logical ,intent(in),optional :: matrixcheck
   logical ,intent(in),optional :: acc
   real(r8)            :: matrix_update_gmn

   associate(                          &
   matrix_nphturnover                   => cnveg_nitrogenflux_inst%matrix_nphturnover_patch , &
   matrix_ngmtransfer                   => cnveg_nitrogenflux_inst%matrix_ngmtransfer_patch , &
   matrix_ngmturnover                   => cnveg_nitrogenflux_inst%matrix_ngmturnover_patch , &
   doner_gmn                            => cnveg_nitrogenflux_inst%matrix_ngmtransfer_doner_patch & ! Input:  [integer (:)] Doners of gap mortality related N transfer
   )

      if(.not. present(matrixcheck) .or. matrixcheck)then
         if((.not. present(acc) .or. acc) .and. matrix_nphturnover(p,doner_gmn(itransfer)) + matrix_ngmturnover(p,doner_gmn(itransfer)) + rate * dt .ge. 1)then
            matrix_update_gmn = max(0._r8,(1._r8 - matrix_nphturnover(p,doner_gmn(itransfer)) - matrix_ngmturnover(p,doner_gmn(itransfer))) / dt)
         else
            matrix_update_gmn = rate
         end if
      else
         matrix_update_gmn = rate
      end if
      if(.not. present(acc) .or. acc)then
         matrix_ngmturnover(p,doner_gmn(itransfer)) = matrix_ngmturnover(p,doner_gmn(itransfer)) + matrix_update_gmn * dt
         matrix_ngmtransfer(p,itransfer) = matrix_ngmtransfer(p,itransfer) + matrix_update_gmn
      else
         matrix_ngmturnover(p,doner_gmn(itransfer)) = matrix_ngmturnover(p,doner_gmn(itransfer)) - matrix_ngmtransfer(p,itransfer) * dt + matrix_update_gmn * dt
         matrix_ngmtransfer(p,itransfer) = matrix_update_gmn
      end if

      return
   end associate

 end function matrix_update_gmn


 function matrix_update_fin(p,itransfer,rate,dt,cnveg_nitrogenflux_inst,matrixcheck,acc)

   integer ,intent(in) :: p
   integer ,intent(in) :: itransfer
   real(r8),intent(in) :: rate
   real(r8),intent(in) :: dt
   type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst 
   logical ,intent(in),optional :: matrixcheck
   logical ,intent(in),optional :: acc
   real(r8)            :: matrix_update_fin

   associate(                          &
   matrix_nphturnover                   => cnveg_nitrogenflux_inst%matrix_nphturnover_patch , &
   matrix_ngmturnover                   => cnveg_nitrogenflux_inst%matrix_ngmturnover_patch , &
   matrix_nfitransfer                   => cnveg_nitrogenflux_inst%matrix_nfitransfer_patch , &
   matrix_nfiturnover                   => cnveg_nitrogenflux_inst%matrix_nfiturnover_patch , &
   doner_fin                           => cnveg_nitrogenflux_inst%matrix_nfitransfer_doner_patch &
   )

      if(.not. present(matrixcheck) .or. matrixcheck)then
         if((.not. present(acc) .or. acc) .and. matrix_nphturnover(p,doner_fin(itransfer)) + matrix_ngmturnover(p,doner_fin(itransfer)) &
           + matrix_nfiturnover(p,doner_fin(itransfer)) + rate * dt .ge. 1)then
            matrix_update_fin = max(0._r8,(1._r8 - matrix_nphturnover(p,doner_fin(itransfer)) &
                              - matrix_ngmturnover(p,doner_fin(itransfer)) - matrix_nfiturnover(p,doner_fin(itransfer))) / dt)
         else
            matrix_update_fin = rate
         end if
      else
         matrix_update_fin = rate
      end if
      if(.not. present(acc) .or. acc)then
         matrix_nfiturnover(p,doner_fin(itransfer)) = matrix_nfiturnover(p,doner_fin(itransfer)) + matrix_update_fin * dt
         matrix_nfitransfer(p,itransfer) = matrix_nfitransfer(p,itransfer) + matrix_update_fin
      else
         matrix_nfiturnover(p,doner_fin(itransfer)) = matrix_nfiturnover(p,doner_fin(itransfer)) - matrix_nfitransfer(p,itransfer) * dt + matrix_update_fin * dt
         matrix_nfitransfer(p,itransfer) = matrix_update_fin
      end if

      return
   end associate

 end function matrix_update_fin

  !-----------------------------------------------------------------------
   subroutine CNVegMatrixRest( ncid, flag )
    ! !DESCRIPTION:
    !
    !    Read/write restart data needed for the CN Matrix model solution
    !
    ! !USES:
    use restUtilMod   , only: restartvar
    use ncdio_pio     , only: file_desc_t, ncd_int
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !------------------------------------------------------------------------
    call restartvar(ncid=ncid, flag=flag, varname='bgc_cycle_year', xtype=ncd_int,  &
            long_name='Year number in spinup cycle sequence', units='years', &
            interpinic_flag='skip', readvar=readvar, data=iyr)

    call restartvar(ncid=ncid, flag=flag, varname='bgc_cycle_loop', xtype=ncd_int,  &
            long_name='Loop number in spinup cycle sequence', units='years', &
            interpinic_flag='skip', readvar=readvar, data=iloop)

    !------------------------------------------------------------------------
   end subroutine CNVegMatrixRest

end module CNVegMatrixMod
