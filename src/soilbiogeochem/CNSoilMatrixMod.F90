module CNSoilMatrixMod

!#include "shr_assert.h"
  !-----------------------------------------------------------------------
  ! The matrix model of CLM5.0 was developed by Yiqi Luo EcoLab members,
  ! Drs. Xingjie Lu, Yuanyuan Huang and Zhengguang Du, at Northern Arizona University
  !----------------------------------------------------------------------------------
  ! 
  ! DESCRIPTION:
  ! Module for CLM5.0BGC matrices
  ! The matrix equation 
  ! Xn+1 = Xn + I*dt + (A*ksi*k - tri/dz)*Xn*dt
  
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use shr_log_mod                    , only : errMsg => shr_log_errMsg
  use decompMod                      , only : bounds_type  
  use abortutils                     , only : endrun
  use clm_time_manager               , only : get_step_size, is_end_curr_month,get_curr_date
  use clm_time_manager               , only : is_first_step_of_this_run_segment,is_beg_curr_year,is_end_curr_year
  use clm_varpar                     , only : ndecomp_pools, nlevdecomp, ndecomp_pools_vr        !number of biogeochemically active soil layers
  use clm_varpar                     , only : ndecomp_cascade_transitions, ndecomp_cascade_outtransitions
  use clm_varpar                     , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varcon                     , only : dzsoi_decomp,zsoi,secspday
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con    
!  use clm_varctl                     , only : use_vertsoilc
  use CNVegCarbonFluxType            , only : cnveg_carbonflux_type
  use CNVegNitrogenFluxType          , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType        , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType  , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType   , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType  , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType   , only : soilbiogeochem_nitrogenflux_type  
  use CNSharedParamsMod                , only : CNParamsShareInst
  use SoilStateType                  , only : soilstate_type  
  use clm_varctl                     , only : isspinup, use_soil_matrixcn, is_outmatrix
  use ColumnType                      , only : col                
  use GridcellType                   , only : grc
  use clm_varctl                     , only : use_c13, use_c14
  use perf_mod                        , only : t_startf, t_stopf
  use SPMMod                         , only : sparse_matrix_type, diag_matrix_type, vector_type
!
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNSoilMatrix
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNSoilMatrix(bounds,num_soilc, filter_soilc, num_actfirec, filter_actfirec,&
       cnveg_carbonflux_inst,soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_carbonflux_inst,soilbiogeochem_state_inst, &
       cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst, &
       soilbiogeochem_nitrogenstate_inst,c13_soilbiogeochem_carbonstate_inst,&
       c13_soilbiogeochem_carbonflux_inst,c14_soilbiogeochem_carbonstate_inst,&
       c14_soilbiogeochem_carbonflux_inst)
    ! !DESCRIPTION:
    ! !ARGUMENTS:
    type(bounds_type)                        , intent(in)    :: bounds
    integer                                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                  , intent(in)    :: num_actfirec       ! number of soil columns in filter
    integer                                  , intent(in)    :: filter_actfirec(:) ! filter for soil columns
    type(cnveg_carbonflux_type)              , intent(inout) :: cnveg_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)    , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)     , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_state_type)          , intent(inout) :: soilbiogeochem_state_inst
    type(cnveg_nitrogenflux_type)            , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenflux_type)   , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type)  , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_carbonstate_type)    , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)     , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)    , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)     , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    ! !LOCAL VARIABLES:
    integer :: fc,j,i, l,k ! indices
    integer :: c     !  
!    integer,parameter:: nspools=7       
    real(r8):: dt                   ! time step (seconds)
    real(r8):: secspyear            ! time step (seconds)
    real(r8):: epsi,fire_delta      !small number
    real(r8):: tmptmp

    integer :: begc,endc                                    ! bounds 
    real(r8),dimension(bounds%begc:bounds%endc,nlevdecomp*(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)) :: a_ma_vr,na_ma_vr
    real(r8),dimension(bounds%begc:bounds%endc,ndecomp_pools_vr) :: kk_ma_vr, kk_fire_vr
!    real(r8),dimension(bounds%begc:bounds%endc,(ndecomp_pools-1)*(3*nlevdecomp-2)) :: tri_ma_vr !,tranvert,ntranvert
    real(r8),dimension(nlevdecomp*(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)) :: Aoned
    real(r8),dimension(ndecomp_pools_vr) :: Koned
    real(r8),dimension((ndecomp_pools-1)*(3*nlevdecomp-2)) :: trioned !,tranvert,ntranvert
    real(r8),dimension(ndecomp_pools_vr) :: AKonedfire
!    real(r8),dimension(ndecomp_pools*(3*nlevdecomp-2)),save :: tri_i,tri_j
!    logical,save :: trij_assigned=.false.
!    real(r8),dimension(ndecomp_pools_vr,ndecomp_pools_vr) :: tranvert_1c,ntranvert_1c,a_ma_vr_1c,kk_ma_vr_1c,na_ma_vr_1c
!    real(r8),dimension(ndecomp_pools,ndecomp_pools) :: a_ma, kk_ma, matrix_soil_cn
!    real(r8),dimension(ndecomp_pools,ndecomp_pools) :: na_ma
!    real(r8),dimension(num_soilc,nlevdecomp) :: depth_scalar,two_scalar
!    real(r8),dimension(num_soilc,nlevdecomp) :: n_scalar_ave,t_scalar_ave,w_scalar_ave, o_scalar_ave
!    real(r8)::  a_tri_ave_c(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
!    real(r8)::  b_tri_ave_c(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
!    real(r8)::  c_tri_ave_c(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
!    real(r8)::  input_ave_c(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
!    real(r8)::  a_tri_ave_n(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
!    real(r8)::  b_tri_ave_n(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
!    real(r8)::  c_tri_ave_n(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
!    real(r8)::  input_ave_n(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)


!    real(r8),dimension(ndecomp_pools_vr,1) :: matrix_Cinter,matrix_Cinter_next,matrix_Cinput_vector,emulator_tmp,emulator_tmp1,emulator_tmpn
!    real(r8),dimension(bounds%begc:bounds%endc,ndecomp_pools_vr) :: matrix_Cinter,matrix_Cinter_next,matrix_Cinput_vector!,emulator_tmp,emulator_tmp1,emulator_tmpn
!    real(r8),dimension(bounds%begc:bounds%endc,ndecomp_pools_vr) :: matrix_Ninter,matrix_Ninter_next,matrix_Ninput_vector

!    real(r8),dimension(bounds%begc:bounds%endc,ndecomp_pools_vr) :: matrix13_Cinter,   matrix14_Cinter
!    real(r8),dimension(bounds%begc:bounds%endc,ndecomp_pools_vr) :: matrix_Cinput13_vector,matrix_Cinput14_vector
!    real(r8),dimension(bounds%begc:bounds%endc,ndecomp_pools_vr) :: matrix_Cinter13_next,matrix_Cinter14_next

!    type(vector_type) :: matrix_Cinter
!    type(vector_type) :: matrix_Cinter_next
!    type(vector_type) :: matrix13_Cinter
!    type(vector_type) :: matrix13_Cinter_next
!    type(vector_type) :: matrix14_Cinter
!    type(vector_type) :: matrix14_Cinter_next
!    type(vector_type) :: matrix_Ninter
!    type(vector_type) :: matrix_Ninter_next
!    type(vector_type) :: matrix_Cinput_vector
!    type(vector_type) :: matrix_Cinput13_vector
!    type(vector_type) :: matrix_Cinput14_vector
!    type(vector_type) :: matrix_Ninput_vector

    real(r8),dimension(bounds%begc:bounds%endc,1:ndecomp_pools_vr,1) ::  soilmatrixc_cap,soilmatrixn_cap
    real(r8), dimension(1:ndecomp_pools_vr,1:ndecomp_pools_vr)   ::  AKinv,AKinvn
!    real(r8),dimension(ndecomp_pools,1) :: emulator_Cinter_1,emulator_Cinter_next_1,emulator_Cinput_1

    real(r8),dimension(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools) :: cn_decomp_pools
!input for 3-step update
!    real(r8), dimension(num_soilc,nlevdecomp) :: tot_c_to_litr_met_c,tot_c_to_litr_cel_c,tot_c_to_litr_lig_c,tot_c_to_cwdc 
!    real(r8), dimension(num_soilc,nlevdecomp) :: tot_n_to_litr_met_n,tot_n_to_litr_cel_n,tot_n_to_litr_lig_n,tot_n_to_cwdn
    
!    real(r8), dimension(num_soilc,nlevdecomp) :: tot_c_to_litr_met_c1,tot_c_to_litr_cel_c1,tot_c_to_litr_lig_c1,tot_c_to_cwdc1 
!    real(r8), dimension(num_soilc,nlevdecomp) :: tot_n_to_litr_met_n1,tot_n_to_litr_cel_n1,tot_n_to_litr_lig_n1,tot_n_to_cwdn1
  
!    real(r8), dimension(num_soilc,nlevdecomp) :: tot_c_to_litr_met_c23,tot_c_to_litr_cel_c23,tot_c_to_litr_lig_c23,tot_c_to_cwdc23 
!    real(r8), dimension(num_soilc,nlevdecomp) :: tot_n_to_litr_met_n23,tot_n_to_litr_cel_n23,tot_n_to_litr_lig_n23,tot_n_to_cwdn23
    integer :: Ntrans
    integer tranlist_a
    integer j_decomp,j_lev
    integer,dimension(:) :: kk(bounds%begc:bounds%endc)
    integer,dimension(:) :: kfire_i(1:ndecomp_pools_vr)
    integer,dimension(:) :: kfire_j(1:ndecomp_pools_vr)
    real(r8),dimension(:,:) :: Cinter_old(bounds%begc:bounds%endc,1:ndecomp_pools_vr)
    real(r8),dimension(:,:) :: Ninter_old(bounds%begc:bounds%endc,1:ndecomp_pools_vr)
    logical,save :: list_ready1 = .False.
    logical,save :: list_ready1_fire   = .False.
    logical,save :: list_ready1_nofire = .False.
    logical,save :: list_ready2_fire   = .False.
    logical,save :: list_ready2_nofire = .False.
    logical,save :: init_readyAsoilc = .False.
    logical,save :: init_readyAsoiln = .False.
    logical,save :: init_readyAVsoil = .False.
    logical,save :: init_readyAKfiresoil = .False.
!    integer,save,dimension(:) :: RI_a (1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp+ndecomp_pools_vr) = -9999
!    integer,save,dimension(:) :: CI_a (1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp+ndecomp_pools_vr) = -9999
!    integer,save,dimension(:) :: RI_na (1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp+ndecomp_pools_vr) = -9999 
!    integer,save,dimension(:) :: CI_na (1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp+ndecomp_pools_vr) = -9999
!   integer,save :: NE_AKVfiresoil=0
!   integer,save :: NE_AKallsoilc=0
!   integer,save :: NE_AKallsoiln=0
!   integer,dimension(:),save :: RI_AKVfiresoil=0


!    integer,dimension(:) :: list_Asoilc(1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp)
!    integer,dimension(:) :: list_Asoiln(1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp)

!    logical,save :: index_ready1 = .False.
!    logical,save :: index_ready2 = .False.
!    logical,save :: index_ready3 = .False.
!    logical,save :: index_ready4 = .False.
    logical isbegofyear

    real(r8) :: hr_adjust

    !-----------------------------------------------------------------------
    begc = bounds%begc; endc = bounds%endc
    
!    SHR_ASSERT_ALL((ubound(cn_decomp_pools)     == (/endc,nlevdecomp,ndecomp_pools/))  , errMsg(sourcefile, __LINE__))
    associate(                                      &
         cf_veg  => cnveg_carbonflux_inst ,          & ! Input
         cs_soil => soilbiogeochem_carbonstate_inst, & ! Output
         cf_soil => soilbiogeochem_carbonflux_inst, & ! Output
         soil_bg => soilbiogeochem_state_inst, &
         nf_veg  => cnveg_nitrogenflux_inst, &
         nf_soil => soilbiogeochem_nitrogenflux_inst, &
         ns_soil => soilbiogeochem_nitrogenstate_inst, &
         cs13_soil => c13_soilbiogeochem_carbonstate_inst, & ! 
         cf13_soil => c13_soilbiogeochem_carbonflux_inst, & ! 
         cs14_soil => c14_soilbiogeochem_carbonstate_inst, & ! 
         cf14_soil => c14_soilbiogeochem_carbonflux_inst, & ! 

!         matrix_decomp_k  => soilbiogeochem_carbonflux_inst%matrix_decomp_k_col, & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         fpi_vr          => soilbiogeochem_state_inst%fpi_vr_col    , & ! Input:  [real(r8) (:,:)   ]  fraction of potential immobilization (no units) 
         cascade_donor_pool => decomp_cascade_con%cascade_donor_pool, &
         cascade_receiver_pool          => decomp_cascade_con%cascade_receiver_pool,& 
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools, &
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio, &
         rf_decomp_cascade => soilbiogeochem_state_inst%rf_decomp_cascade_col,&
         pathfrac_decomp_cascade => soilbiogeochem_state_inst%pathfrac_decomp_cascade_col, &
         is_cwd                         => decomp_cascade_con%is_cwd                                         , & ! Input:  [logical  (:)     ]  TRUE => pool is a litter pool
         is_litter                      => decomp_cascade_con%is_litter                                      , & ! Input:  [logical  (:)     ]  TRUE => pool is a litter pool
         hr               => soilbiogeochem_carbonflux_inst%hr_col, &
         er               => cnveg_carbonflux_inst%er_col, &
         m_decomp_cpools_to_fire => cnveg_carbonflux_inst%m_decomp_cpools_to_fire_col, &
!         matrix_a_tri     => soilbiogeochem_carbonflux_inst%matrix_a_tri_col, &!(begc:endc,1:nlevdecomp),   & ! Output: "A"-for matrix
!         matrix_b_tri     => soilbiogeochem_carbonflux_inst%matrix_b_tri_col, &!(begc:endc,1:nlevdecomp),   & ! Output: "B"-for matrix
!         matrix_c_tri     => soilbiogeochem_carbonflux_inst%matrix_c_tri_col, &!(begc:endc,1:nlevdecomp),   & ! Output: "C"-for matrix 
         tri_ma_vr     => soilbiogeochem_carbonflux_inst%tri_ma_vr, &!(begc:endc,1:nlevdecomp),   & 
!         matrix_Cinput     => soilbiogeochem_carbonflux_inst%matrix_input_col, &!(begc:endc,1:nlevdecomp,1:ndecomp_pools),   & ! Output: "C"-for matrix 
!         matrix_Ninput     => soilbiogeochem_nitrogenflux_inst%matrix_input_col, &!(begc:endc,1:nlevdecomp,1:ndecomp_pools),  & ! Output: "C"-for matrix 
         matrix_decomp_fire_k => soilbiogeochem_carbonflux_inst%matrix_decomp_fire_k_col, &!(begc:endc,1:nlevdecomp,1:ndecomp_pools) &
         all_i                 => decomp_cascade_con%all_i, &
         all_j                 => decomp_cascade_con%all_j, &
         A_i                   => decomp_cascade_con%A_i, &
         A_j                   => decomp_cascade_con%A_j, &
         tri_i                 => decomp_cascade_con%tri_i, &
         tri_j                 => decomp_cascade_con%tri_j, &
         spm_tranlist_a        => decomp_cascade_con%spm_tranlist_a, &
         n_all_entries         => decomp_cascade_con%n_all_entries, &
         Ntri_setup            => decomp_cascade_con%Ntri_setup, &
         Ntrans_setup          => decomp_cascade_con%Ntrans_setup, &
!         Asoilc                => soilbiogeochem_carbonflux_inst%Asoilc, &
!         Asoiln                => soilbiogeochem_nitrogenflux_inst%Asoiln, &
         AKsoilc               => soilbiogeochem_carbonflux_inst%AKsoilc, &
         AKsoiln               => soilbiogeochem_nitrogenflux_inst%AKsoiln, &
         AVsoil                => soilbiogeochem_carbonflux_inst%AVsoil, &
!         AKfiresoil            => soilbiogeochem_carbonflux_inst%AKfiresoil, &
         AKfiresoil            => soilbiogeochem_carbonflux_inst%AKfiresoil, &
!         NE_AKfiresoil        => soilbiogeochem_carbonflux_inst%NE_AKfiresoil, &
!         RI_AKfiresoil        => soilbiogeochem_carbonflux_inst%RI_AKfiresoil, &
!         CI_AKfiresoil        => soilbiogeochem_carbonflux_inst%CI_AKfiresoil, &
         AKallsoilc            => soilbiogeochem_carbonflux_inst%AKallsoilc, &
         NE_AKallsoilc         => soilbiogeochem_carbonflux_inst%NE_AKallsoilc, &
         RI_AKallsoilc         => soilbiogeochem_carbonflux_inst%RI_AKallsoilc, &
         CI_AKallsoilc         => soilbiogeochem_carbonflux_inst%CI_AKallsoilc, &
         AKallsoiln            => soilbiogeochem_nitrogenflux_inst%AKallsoiln, &
         NE_AKallsoiln         => soilbiogeochem_nitrogenflux_inst%NE_AKallsoiln, &
         RI_AKallsoiln         => soilbiogeochem_nitrogenflux_inst%RI_AKallsoiln, &
         CI_AKallsoiln         => soilbiogeochem_nitrogenflux_inst%CI_AKallsoiln, &
         RI_a                  => soilbiogeochem_carbonflux_inst%RI_a, &
         CI_a                  => soilbiogeochem_carbonflux_inst%CI_a, &
         RI_na                 => soilbiogeochem_nitrogenflux_inst%RI_na, &
         CI_na                 => soilbiogeochem_nitrogenflux_inst%CI_na, &
!         AKXcinc               => soilbiogeochem_carbonflux_inst%AKXcinc, &
!         AKXninc               => soilbiogeochem_nitrogenflux_inst%AKXninc, &
         Ksoil                 => soilbiogeochem_carbonflux_inst%Ksoil, &
         Xdiagsoil             => soilbiogeochem_carbonflux_inst%Xdiagsoil, &
         matrix_Cinter      => soilbiogeochem_carbonstate_inst%matrix_Cinter, &
         matrix_Ninter      => soilbiogeochem_nitrogenstate_inst%matrix_Ninter, &
!         matrix_Cinter_next => soilbiogeochem_carbonflux_inst%matrix_Cinter_next, &
!         matrix_Ninter_next => soilbiogeochem_nitrogenflux_inst%matrix_Ninter_next, &
         matrix_Cinter13      => c13_soilbiogeochem_carbonstate_inst%matrix_Cinter, &
         matrix_Cinter14      => c14_soilbiogeochem_carbonstate_inst%matrix_Cinter, &
!         matrix_Cinter13_next => soilbiogeochem_carbonflux_inst%matrix_Cinter13_next, &
!         matrix_Cinter14_next => soilbiogeochem_carbonflux_inst%matrix_Cinter14_next, &
         matrix_Cinput    => soilbiogeochem_carbonflux_inst%matrix_Cinput, &
         matrix_Cinput13  => c13_soilbiogeochem_carbonflux_inst%matrix_Cinput, &
         matrix_Cinput14  => c14_soilbiogeochem_carbonflux_inst%matrix_Cinput, &
         matrix_Ninput    => soilbiogeochem_nitrogenflux_inst%matrix_Ninput, &
         AKXcacc                 => soilbiogeochem_carbonstate_inst%AKXcacc, &
         AKXnacc                 => soilbiogeochem_nitrogenstate_inst%AKXnacc, &
         list_Asoilc             => decomp_cascade_con%list_Asoilc, &
         list_Asoiln             => decomp_cascade_con%list_Asoiln, &
         list_V_AKVfire          => decomp_cascade_con%list_V_AKVfire, &
         list_fire_AKVfire       => decomp_cascade_con%list_fire_AKVfire, &
         list_AK_AKVfire         => decomp_cascade_con%list_AK_AKVfire, &
         list_AK_AKV             => decomp_cascade_con%list_AK_AKV, &
         list_V_AKV              => decomp_cascade_con%list_V_AKV &
         )

     !print*,'here0,NE',AKXcacc%NE
     !print*,'here0'
!            do i=1,ndecomp_pools_vr
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  print*,'here6.22',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!               end do
!            end do
     ! set time steps
      call t_startf('CN Soil matrix-init. matrix')
      dt = real( get_step_size(), r8 )
!         print*,'begin decomp0',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)

      Ntrans = (ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp
      epsi = 1.e-30_r8 

      isbegofyear = is_beg_curr_year()

      ! calculate c:n ratios of applicable pools
      do l = 1, ndecomp_pools
         if ( floating_cn_ratio_decomp_pools(l)) then
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                  if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,'here0.0a',c,j,l,cs_soil%decomp_cpools_vr_col(c,j,l),ns_soil%decomp_npools_vr_col(c,j,l)
                  if ( ns_soil%decomp_npools_vr_col(c,j,l) > 0._r8 ) then
                     cn_decomp_pools(c,j,l) = cs_soil%decomp_cpools_vr_col(c,j,l) / ns_soil%decomp_npools_vr_col(c,j,l)
                  else
                     cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
                  end if
!                  if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,'here0.1a',cn_decomp_pools(c,j,l)
               end do
            end do
         else
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                  if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,'here0.0b',initial_cn_ratio(l)
                  cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
!                  if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,'here0.1b',cn_decomp_pools(c,j,l)
               end do
            end do
         end if 
      end do
!         cn_decomp_pools(:,:,l) = initial_cn_ratio(l)
!         where(spread(spread(floating_cn_ratio_decomp_pools(:),1,nlevdecomp),1,num_soilc) &
!           .and. cs_soil%decomp_cpools_vr_col(:,:,:) > 0._r8 .and. ns_soil%decomp_npools_vr_col(:,:,:) > 0._r8)
!            cn_decomp_pools(:,:,:) = cs_soil%decomp_cpools_vr_col(:,:,:) / ns_soil%decomp_npools_vr_col(:,:,:)
!         elsewhere
!            cn_decomp_pools(:,:,:) = spread(spread(initial_cn_ratio,1,nlevdecomp),1,num_soilc)
!         end where
    
!      a_ma_vr = 0.0_r8
!      na_ma_vr = 0.0_r8
!      kk_ma_vr = 0.0_r8   ! kk_matrix, decay matrix * scalar matrix
!      kk_fire_vr = 0.0_r8   ! kk_matrix, decay matrix * scalar matrix
!      tri_ma_vr = 0.0_r8
!      matrix_Cinput_vector = 0.0_r8
!      matrix_Ninput_vector = 0.0_r8

!      AKinv(:,:)=0._r8    
!      AKinvn(:,:)=0._r8    
!      matrix_Cinter =0.0_r8
!      matrix_Ninter =0.0_r8

      call t_stopf('CN Soil matrix-init. matrix')
!       if (use_vertsoilc) then
!      do j=1,ndecomp_pools_vr   !70
!         do fc = 1,num_soilc
!            c = filter_soilc(fc)
!            a_ma_vr(c,j,j) = -1.0_r8
!            na_ma_vr(c,j,j)= -1.0_r8
!         end do
!      enddo
     !print*,'here1'
            do i=1,ndecomp_pools_vr
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                  print*,'here6.22',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
               end do
            end do
!             cmatrix_in    = 0.0_r8
      call t_startf('CN Soil matrix-assign matrix-a-na')
!      do k = 1, ndecomp_cascade_transitions
!         do j = 1, nlevdecomp
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               if(cascade_receiver_pool(k) .ne. 0)then  !transition to atmosphere
!                  a_ma_vr(c,spm_tranlist_a(j,k)) = (1.0-rf_decomp_cascade(c,j,k))*pathfrac_decomp_cascade(c,j,k)
!                  if( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)))then
!                     na_ma_vr(c,spm_tranlist_a(j,k)) = (1.0-rf_decomp_cascade(c,j,k))* &
!                               (cn_decomp_pools(c,j,cascade_donor_pool(k))/cn_decomp_pools(c,j,cascade_receiver_pool(k)))*pathfrac_decomp_cascade(c,j,k)
!                  else
!                     na_ma_vr(c,spm_tranlist_a(j,k)) = pathfrac_decomp_cascade(c,j,k)
!                  end if
!               end if
!            end do
!         end do
!      end do

      do k = 1, ndecomp_cascade_transitions
         !print*,'here1.1',k
         !print*,cascade_receiver_pool(k)
         if(cascade_receiver_pool(k) .ne. 0)then  !transition to atmosphere
            !print*,'here1.2'
            do j = 1, nlevdecomp
               tranlist_a = spm_tranlist_a(j,k)
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                  if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,'here1.3',c,j,k,tranlist_a,cascade_donor_pool(k)
                  a_ma_vr(c,tranlist_a) = (1.0_r8-rf_decomp_cascade(c,j,k))*pathfrac_decomp_cascade(c,j,k)
                  if( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)))then
!                     print*,'here1.4a'
!                     if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,rf_decomp_cascade(c,j,k)
!                     if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,cn_decomp_pools(c,j,cascade_donor_pool(k))
!                     if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,cn_decomp_pools(c,j,cascade_receiver_pool(k))
!                     if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,pathfrac_decomp_cascade(c,j,k) 
!                     if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,cn_decomp_pools(c,j,:)
                     na_ma_vr(c,tranlist_a) = (1.0_r8-rf_decomp_cascade(c,j,k))* &
                               (cn_decomp_pools(c,j,cascade_donor_pool(k))/cn_decomp_pools(c,j,cascade_receiver_pool(k)))*pathfrac_decomp_cascade(c,j,k)
                  else
!                     print*,'here1.4b'
                     na_ma_vr(c,tranlist_a) = pathfrac_decomp_cascade(c,j,k)
                  end if
               end do
            end do
         end if
      end do

!         do j=1,Ntrans
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'here6.238',c,A_i(j),A_j(j),j,a_ma_vr(c,j),num_actfirec
!            end do
!         end do
      call t_stopf('CN Soil matrix-assign matrix-a-na')
    !print*,'here2'
!            do i=1,ndecomp_pools_vr
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  print*,'here6.22',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!               end do
!            end do
      call t_startf('CN Soil matrix-assign matrix-in')
      do i = 1,ndecomp_pools
         if(is_litter(i))then
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!               matrix_Cinput_vector%V(c,j+(i-1)*nlevdecomp) = max(matrix_Cinput(c,j,i),epsi) 
!               matrix_Ninput_vector%V(c,j+(i-1)*nlevdecomp) = max(matrix_Ninput(c,j,i),epsi)
!               if ( use_c13 ) then
!                  matrix_Cinput13_vector%V(c,j+(i-1)*nlevdecomp) = max(cf13_soil%matrix_input_col(c,j,i),epsi)
!               end if !c13
!               if ( use_c14 ) then
!                  matrix_Cinput14_vector%V(c,j+(i-1)*nlevdecomp) = max(cf14_soil%matrix_input_col(c,j,i),epsi)
!               end if !c14
!               kk_fire_vr(c,(i-1)*nlevdecomp+j) = matrix_decomp_fire_k(c,j,i) 
                  Ksoil%DM(c,(i-1)*nlevdecomp+j)   = Ksoil%DM(c,(i-1)*nlevdecomp+j) * fpi_vr(c,j) 
               end do
            end do
         end if
      end do

!      if(bounds%begc .le. 32397 .and. bounds%endc .ge. 32397)print*,'Ksoil',Ksoil%DM(32397,8)
      call t_stopf('CN Soil matrix-assign matrix-in')
      call t_startf('CN Soil matrix-assign matrix-tri')

!      kk = 1
!      do i = 1,ndecomp_pools
!         if(.not. is_cwd(i))then
!            do j = 1, nlevdecomp
!               if (j > 1) then ! avoid tranfer from for example,soil1_1st layer to litr3_10th layer
!                  do fc = 1,num_soilc
!                     c = filter_soilc(fc)
!                     tri_ma_vr(c,kk(c)) = matrix_c_tri(c,j-1) / dzsoi_decomp(j-1) * (-dt)
!               !      print*,'tri_ma_vr',c,kk(c),j,tri_ma_vr(c,kk(c)),matrix_c_tri(c,j),dzsoi_decomp(j)
!                     kk(c) = kk(c) + 1
!                  end do
!               end if
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  tri_ma_vr(c,kk(c)) = matrix_b_tri(c,j) / dzsoi_decomp(j) * (-dt)
!               !   print*,'tri_ma_vr',c,kk(c),j,tri_ma_vr(c,kk(c)),matrix_b_tri(c,j),dzsoi_decomp(j)
!                  kk(c) = kk(c) + 1
!               end do
!               if (j < nlevdecomp) then ! avoid tranfer from for example, litr3_10th layer to soil1_1st layer
!                  do fc = 1,num_soilc
!                     c = filter_soilc(fc)
!                     tri_ma_vr(c,kk(c)) = matrix_a_tri(c,j+1) / dzsoi_decomp(j+1) * (-dt)
!               !      print*,'tri_ma_vr',c,kk(c),j,tri_ma_vr(c,kk(c)),matrix_a_tri(c,j),dzsoi_decomp(j)
!                     kk(c) = kk(c) + 1
!                  end do
!               end if
!            end do
!         end if
!      end do

      call t_stopf('CN Soil matrix-assign matrix-tri')
      call t_startf('CN Soil matrix-assign matrix-inter')
      do i = 1,ndecomp_pools
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
!               if(cs_soil%decomp_cpools_vr_col(c,j,i) .ne. 0)then
                  matrix_Cinter%V(c,j+(i-1)*nlevdecomp)   = cs_soil%decomp_cpools_vr_col(c,j,i)                
                  Cinter_old(c,j+(i-1)*nlevdecomp)        = cs_soil%decomp_cpools_vr_col(c,j,i)
!               end if
!               if(ns_soil%decomp_npools_vr_col(c,j,i) .ne. 0)then
                  matrix_Ninter%V(c,j+(i-1)*nlevdecomp)   = ns_soil%decomp_npools_vr_col(c,j,i)     
                  Ninter_old(c,j+(i-1)*nlevdecomp)        = ns_soil%decomp_npools_vr_col(c,j,i)
!               end if
            end do
         end do
      end do
!         do j=1,ndecomp_pools_vr
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'matrix_Cinter%V here3.2',c,j,matrix_Cinter%V(c,j)
!            end do
!         end do
      if ( use_c13 )then
         do i = 1,ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                  if(cs13_soil%decomp_cpools_vr_col(c,j,i) .ne. 0)then
                     matrix_Cinter13%V(c,j+(i-1)*nlevdecomp)   = cs13_soil%decomp_cpools_vr_col(c,j,i)
!                  end if
               end do
            end do
         end do
      end if !c13

      if ( use_c14 )then
         do i = 1,ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                  if(cs14_soil%decomp_cpools_vr_col(c,j,i) .ne. 0)then
                     matrix_Cinter14%V(c,j+(i-1)*nlevdecomp)   = cs14_soil%decomp_cpools_vr_col(c,j,i)
!                  end if
               end do
            end do
         end do
      end if !c14
!      if(c .eq. 34208)print*,'matrix_Cinter',matrix_Cinter(c,1:nlevdecomp)
!      if(c .eq. 34208)print*,'matrix_Cinput',matrix_Cinput_vector(c,1:nlevdecomp)
      call t_stopf('CN Soil matrix-assign matrix-inter')

      
      call t_startf('CN Soil matrix-assign matrix-decomp0')

!      print*,'before decomp0',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)!,cs_soil%decomp0_cpools_vr_col(:,9,10)
      if (isbegofyear)then  
         do i = 1,ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs_soil%decomp0_cpools_vr_col(c,j,i)=cs_soil%decomp_cpools_vr_col(c,j,i)
                  ns_soil%decomp0_npools_vr_col(c,j,i)=ns_soil%decomp_npools_vr_col(c,j,i)
               end do
            end do
         end do
!                  cs_soil%decomp0_cpools_vr_col(c,j,i) = max(cs_soil%decomp_cpools_vr_col(c,j,i),epsi)
!                  ns_soil%decomp0_npools_vr_col(c,j,i) = max(ns_soil%decomp_npools_vr_col(c,j,i),epsi)
         where(cs_soil%decomp0_cpools_vr_col .lt. epsi)
            cs_soil%decomp0_cpools_vr_col = epsi
         end where
         where(ns_soil%decomp0_npools_vr_col .lt. epsi)
            ns_soil%decomp0_npools_vr_col = epsi
         end where
!         print*,'after decomp0',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)!,cs_soil%decomp0_cpools_vr_col(:,9,10)
      end if
      call t_stopf('CN Soil matrix-assign matrix-decomp0')
   

     !print*,'here3'
!            do i=1,ndecomp_pools_vr
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  print*,'here6.22',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!               end do
!            end do
!          tranvert = matmul(a_ma_vr,kk_ma_vr)-tri_ma_vr-kk_fire_vr  !intermediate calculatio
!          ntranvert = matmul(na_ma_vr,kk_ma_vr)-tri_ma_vr-kk_fire_vr  !intermediate calculatio

!         do fc = 1,num_soilc
!            c = filter_soilc(fc)

      call t_startf('CN Soil matrix-matrix mult1-lev3-SetValueAK1')

!            do j=1,Ntrans
!               Aoned = a_ma_vr(:,c)
!            end do
!         do j=1,Ntrans
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
               !print*,'here6.2385',a_ma_vr(c,j),c,A_i(j),A_j(j),j
!            end do
!         end do
      !Set C transfer matrix A from a_ma_vr
!      print*,'before AK%setValue',num_soilc,size(filter_soilc),size(a_ma_vr),size(A_i),size(A_j),size(list_Asoilc),size(RI_a),size(CI_a),begc,endc
      call AKsoilc%SetValueA(begc,endc,num_soilc,filter_soilc,a_ma_vr,A_i,A_j,Ntrans,init_readyAsoilc,list_Asoilc,RI_a,CI_a)
!         do j=1,Ntrans
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'here6.239',AKsoilc%M(c,j),c,AKsoilc%RI(j),AKsoilc%CI(j),j
!            end do
!         end do
      !print*,'here3.1'

      call t_stopf('CN Soil matrix-matrix mult1-lev3-SetValueAK1')

      call t_startf('CN Soil matrix-matrix mult1-lev3-SetValueAK2')

      !print*,'here3.2'
!            Aoned = na_ma_vr(c,:)
      !Set N transfer matrix A from a_ma_vr
      call AKsoiln%SetValueA(begc,endc,num_soilc,filter_soilc,na_ma_vr,A_i,A_j,Ntrans,init_readyAsoiln,list_Asoiln,RI_na,CI_na)

      call t_stopf('CN Soil matrix-matrix mult1-lev3-SetValueAK2')

     !print*,'here3.1'
!            Koned = matrix_decomp_k(c,:) * dt
      ! Set decomposition matrix K
!      call Ksoil%SetValueDM(num_soilc,filter_soilc,matrix_decomp_k)

            
      call t_startf('CN Soil matrix-matrix mult1-lev3-SPMM_AK1')
      ! calculate matrix Ac*K for C
      !print*,'here3.3'
      call AKsoilc%SPMM_AK(num_soilc,filter_soilc,Ksoil)
      !print*,'here3.4'

      call t_stopf('CN Soil matrix-matrix mult1-lev3-SPMM_AK1')

      call t_startf('CN Soil matrix-matrix mult1-lev3-SPMM_AK2')
      ! calculate matrix An*K for N
      call AKsoiln%SPMM_AK(num_soilc,filter_soilc,Ksoil)
      !print*,'here3.5'
      call t_stopf('CN Soil matrix-matrix mult1-lev3-SPMM_AK2')
!         index_ready1 = .True.

!         trioned = tri_ma_vr(c,:) * (-dt)
!     print*,'here3.2',Ntri_setup,c
!     if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'tri_i',tri_i
!     if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'tri_j',tri_j
!     if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'tri_ma_vr',Ntri_setup,tri_ma_vr
!     print*,'trioned',trioned
!      init_readyAVsoil = .False.
!      print*,'init_readyAVsoil1',init_readyAVsoil
      call t_startf('CN Soil matrix-matrix mult1-lev3-SetValueAV,AKfire')
      ! Set vertical transfer matrix V from tri_ma_vr
!      print*,'tri_i',tri_i
!      print*,'tri_j',tri_j
!      print*,'Ntri_setup',Ntri_setup
      call AVsoil%SetValueSM(begc,endc,num_soilc,filter_soilc,tri_ma_vr(begc:endc,1:Ntri_setup),tri_i,tri_j,Ntri_setup)!,init_readyAVsoil)
!      print*,'AVsoil%NE',AVsoil%NE
!      print*,'AVsoil%RI',AVsoil%RI
!      print*,'AVsoil%CI',AVsoil%CI
!      print*,'init_readyAVsoil2',init_readyAVsoil
!      print*,'AVsoil',AVsoil%NE,Ntri_setup
!     print*,'here3.25'

      !print*,'here3.6'
!            AKonedfire = matrix_decomp_fire_k(c,:) * (-dt)
!      if(.not. init_readyAKfiresoil)then
      do j=1,ndecomp_pools_vr
         kfire_i(j) = j
         kfire_j(j) = j
      end do
!      end if
!     print*,'here3.3'
      !print*,'here3.7'
      ! Set fire decomposition matrix Kfire from matrix_decomp_fire_k
!      call AKVfiresoil%SetValueSM(num_soilc,filter_soilc,matrix_decomp_fire_k,kfire_i,kfire_j,ndecomp_pools_vr,init_readyAKfiresoil)
      call AKfiresoil%SetValueSM(begc,endc,num_soilc,filter_soilc,matrix_decomp_fire_k(begc:endc,1:ndecomp_pools_vr),kfire_i,kfire_j,ndecomp_pools_vr)

      call t_stopf('CN Soil matrix-matrix mult1-lev3-SetValueAV,AKfire')
!     print*,'here3.4'
!            if(c .eq. 1)print*,'matrix_cinupt',matrix_Cinput_vector(c,1:3)
!            if(c .eq. 1)print*,'matrix_Cinter',matrix_Cinter(c,1:20)
      call t_startf('CN Soil matrix-matrix mult1-lev3-SPMP_AB')

!            if(matrix_decomp_fire_k(c,1) .eq. 0)then
!               call SPMP_AB(AVsoil ,AKfiresoil ,AKVfiresoil,list_ready,list_A=list_V_Vfire,list_B=list_fire_Vfire)
!               print*,'AVsoil',AVsoil%NE,AVsoil%RI(1:AVsoil%NE),AVsoil%CI(1:AVsoil%NE),AVsoil%M(1:AVsoil%NE)
!               print*,'AKVfiresoil',AKVfiresoil%NE,AKVfiresoil%RI(1:AVsoil%NE),AKVfiresoil%CI(1:AVsoil%NE),AKVfiresoil%M(1:AKVfiresoil%NE)-AVsoil%M(1:AVsoil%NE)
!               call SPMP_AB(AKsoilc,AVsoil,AKallsoilc ,list_ready2,list_A=list_AK_AKAV ,list_B=list_AV_AKAV )
!               call SPMP_AB(AKsoiln,AVsoil,AKallsoiln ,list_ready2,list_A=list_AK_AKAV ,list_B=list_AV_AKAV )
!               list_ready2 = .True.
!            else
      ! Calculate Kfire + V and save to variable AKVfiresoil
!      print*,'list_ready1',list_ready1,AKVfiresoil%NE,list_fire_Vfire
!      call AKVfiresoil%SPMP_AB(num_soilc,filter_soilc,AVsoil,list_ready1,list_A=list_fire_Vfire,list_B=list_V_Vfire,&
!                NE_AB=NE_AKVfiresoil,RI_AB=RI_AKVfiresoil,CI_AB=CI_AKVfiresoil)
!      print*,'AKVfiressoil,list_A',AKVfiresoil%NE,list_V_Vfire(:)
!      print*,'AKVfiressoil,list_B',AKVfiresoil%NE,list_fire_Vfire(:)
!      print*,'AKVfiressoil,RI',AKVfiresoil%NE,RI_AKVfiresoil
!      print*,'AKVfiressoil,CI',AKVfiresoil%NE,CI_AKVfiresoil
       
      ! Calculate Kfire + V + Ac*K and save to variable AKVallsoilc
!      print*,'num_actfirec',num_actfirec,AKallsoilc%NE,list_AK_AKV
!      print*,'list_V_AKV',list_V_AKV
!      print*,'num_actfirec',num_actfirec,filter_actfirec
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'before SPMP_ABC,AKsoilc',AKsoilc%M(8,:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'RI',AKsoilc%RI(:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'CI',AKsoilc%CI(:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'before SPMP_ABC,AVsoil',AVsoil%M(8,:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'RI',AVsoil%RI(:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'CI',AVsoil%CI(:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'before SPMP_ABC,AKfiresoil',AKfiresoil%M(8,:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'RI',AKfiresoil%RI(:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'CI',AKfiresoil%CI(:)
!         do j=1,Ntrans
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'here6.240',c,A_i(j),A_j(j),AKsoilc%RI(j),AKsoilc%CI(j),j,AKsoilc%M(c,j),num_actfirec
!            end do
!         end do
!         do j=1,n_all_entries
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'here6.241',c,j,all_i(j),all_j(j),AKallsoilc%M(c,j),num_actfirec
!            end do
!         end do
      !print*,'here3.8'
      if(num_actfirec .eq. 0)then
         call AKallsoilc%SPMP_AB(num_soilc,filter_soilc,AKsoilc,AVsoil,list_ready1_nofire,list_A=list_AK_AKV, list_B=list_V_AKV,&
              NE_AB=NE_AKallsoilc,RI_AB=RI_AKallsoilc,CI_AB=CI_AKallsoilc)
         call AKallsoiln%SPMP_AB(num_soilc,filter_soilc,AKsoiln,AVsoil,list_ready2_nofire,list_A=list_AK_AKV, list_B=list_V_AKV,&
              NE_AB=NE_AKallsoiln,RI_AB=RI_AKallsoiln,CI_AB=CI_AKallsoiln)
      else
         call AKallsoilc%SPMP_ABC(num_soilc,filter_soilc,AKsoilc,AVsoil,AKfiresoil,list_ready1_fire,list_A=list_AK_AKVfire,&
              list_B=list_V_AKVfire,list_C=list_fire_AKVfire,NE_ABC=NE_AKallsoilc,RI_ABC=RI_AKallsoilc,CI_ABC=CI_AKallsoilc,&
              use_actunit_list_C=.True.,num_actunit_C=num_actfirec,filter_actunit_C=filter_actfirec)
!         print*,'call AKallsoiln_SPMP_ABC'
         call AKallsoiln%SPMP_ABC(num_soilc,filter_soilc,AKsoiln,AVsoil,AKfiresoil,list_ready2_fire,list_A=list_AK_AKVfire,&
              list_B=list_V_AKVfire,list_C=list_fire_AKVfire,NE_ABC=NE_AKallsoiln,RI_ABC=RI_AKallsoiln,CI_ABC=CI_AKallsoiln,&
              use_actunit_list_C=.True.,num_actunit_C=num_actfirec,filter_actunit_C=filter_actfirec)
      end if
!         do j=1,n_all_entries
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'here6.242',c,j,all_i(j),all_j(j),AKallsoilc%M(c,j)
!            end do
!         end do
!      print*,'after SPMP_ABC,list_V_AKV',list_V_AKV
!      print*,'after SPMP_ABC,num_actfirec',num_actfirec,AKallsoilc%NE,list_AK_AKV
!      print*,'AKsoilc',AKsoilc%M(1,140:145)
!      print*,'RI',AKsoilc%RI(140:145)
!      print*,'CI',AKsoilc%CI(140:145)
!      print*,'AKVfiresoilc',AKVfiresoil%M(1,:)
!      print*,'RI',AKVfiresoil%RI(:)
!      print*,'CI',AKVfiresoil%CI(:)
!      print*,'AKallsoilc',AKallsoilc%M(1,79),AKallsoilc%M(1,82),AKallsoilc%M(1,258:260)
!      print*,'RI',AKallsoilc%RI(79),AKallsoilc%RI(82),AKallsoilc%RI(258:260)
!      print*,'CI',AKallsoilc%CI(79),AKallsoilc%CI(82),AKallsoilc%CI(258:260)
      ! Calculate Kfire + V + An*K and save to variable AKVallsoiln
!               list_ready1 = .True.
!            end if

            !call SPMP_AB(AVsoil,AKfiresoil,AKVfiresoil)
!            if(c .eq. 1)print*,'AKsoilc',AKsoilc%M(1:10)
!            if(c .eq. 1)print*,'AKsoilc_RI',AKsoilc%RI(1:10)
!            if(c .eq. 1)print*,'AKsoilc_CI',AKsoilc%CI(1:10)
!            if(c .eq. 1)print*,'AVsoil',AVsoil%M(1:10)
!            if(c .eq. 1)print*,'AVsoil_RI',AVsoil%RI(1:10)
!            if(c .eq. 1)print*,'AVsoil_CI',AVsoil%CI(1:10)
!            if(c .eq. 1)print*,'AKfiresoil',AKfiresoil%M(1:10)
!            if(c .eq. 1)print*,'AKfiresoil_RI',AKfiresoil%RI(1:10)
!            if(c .eq. 1)print*,'AKfiresoil_CI',AKfiresoil%CI(1:10)
!            if(c .eq. 1)print*,'AKVfiresoil',AKVfiresoil%M(1:10)
!            if(c .eq. 1)print*,'AKVfiresoil_RI',AKVfiresoil%RI(1:10)
!            if(c .eq. 1)print*,'AKVfiresoil_CI',AKVfiresoil%CI(1:10)

!     print*,'here3.5'
!            call SPMP_AB(AKsoilc,AKVfiresoil,AKallsoilc)
!     print*,'here3.6'


      call t_stopf('CN Soil matrix-matrix mult1-lev3-SPMP_AB')
!            call t_stopf('CN Soil matrix-matrix mult1-lev3')

      call t_startf('CN Soil matrix-matrix mult2-lev2')

!     print*,'here3.7'
!            if(c .eq. 1)print*,'AKallsoilc',AKallsoilc%M(1:10)
!            if(c .eq. 1)print*,'AKallsoilc_RI',AKallsoilc%RI(1:10)
!            if(c .eq. 1)print*,'AKallsoilc_CI',AKallsoilc%CI(1:10)
 
!      call matrix_Cinter%SetValueV(num_soilc,filter_soilc,matrix_Cinter)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'CinterV before SPMM_AX',matrix_Cinter%V(8,:),AKallsoilc%NE
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'Akallsoilc',AKallsoilc%M(8,:)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'RI',AKallsoilc%RI
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'CI',AKallsoilc%CI
!         do j=1,n_all_entries
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'here3.3',c,j,all_i(j),all_j(j),AKallsoilc%M(c,j)
!            end do
!         end do
!         do j=1,ndecomp_pools_vr
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'matrix_Cinter%V here3.4',c,j,matrix_Cinter%V(c,j)
!            end do
!         end do
      !print*,'here3.9'
      call matrix_Cinter%SPMM_AX(num_soilc,filter_soilc,AKallsoilc)
!         do j=1,ndecomp_pools_vr
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'matrix_Cinter%V here3.5',c,j,matrix_Cinter%V(c,j)
!            end do
!         end do
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'CinterV after SPMM_AX',matrix_Cinter%V(8,:)
!      print*,'matrix_Cinter',matrix_Cinter%V(1,21:)
!      print*,'matrix_Cinter_next',matrix_Cinter_next%V(1,21:)
            
!      call matrix_Ninter%SetValueV(num_soilc,filter_soilc,matrix_Ninter)
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'NinterV before SPMM_AX',matrix_Ninter%V(8,:)
!      if(bounds%begc .le. 32397 .and. bounds%endc .ge. 32397)print*,'Ninter before SPMMAX',matrix_Ninter%V(32397,8)
!     if(bounds%begc .le. 1411 .and. bounds%endc .ge. 1411)print*,'before SPMMAX',matrix_Ninter%V(1411,21:30),AKallsoiln%M(1411,79),AKallsoiln%M(1411,82),AKallsoiln%M(1411,258),AKallsoiln%RI(79),AKallsoiln%RI(82),AKallsoiln%RI(258)!,AKallsoiln%RI,Akallsoilc%CI
      call matrix_Ninter%SPMM_AX(num_soilc,filter_soilc,AKallsoiln)
!     if(bounds%begc .le. 1411 .and. bounds%endc .ge. 1411)print*,'after SPMMAX',matrix_Ninter%V(1411,21:30)      
!      if(bounds%begc .le. 8 .and. bounds%endc .ge. 8)print*,'NinterV after SPMM_AX',matrix_Ninter%V(8,:)
!      if(bounds%begc .le. 32397 .and. bounds%endc .ge. 32397)print*,'Ninter after SPMMAX',matrix_Ninter%V(32397,8)

!     print*,'here3.8'
      if ( use_c13)then
!         call matrix_Cinter13%SetValueV(num_soilc,filter_soilc,matrix_Cinter13)
         call matrix_Cinter13%SPMM_AX(num_soilc,filter_soilc,AKallsoilc)
      end if

      if ( use_c14)then
!         call matrix_Cinter14%SetValueV(num_soilc,filter_soilc,matrix_Cinter14)
         call matrix_Cinter14%SPMM_AX(num_soilc,filter_soilc,AKallsoilc)
      end if

!      do j=1,ndecomp_pools_vr
!         do fc = 1,num_soilc
!            c = filter_soilc(fc)
!            print*,'matrix_Cinter%V here3.6',c,j,matrix_Cinter%V(c,j)
!         end do
!      end do
      do j = 1, ndecomp_pools_vr
         do fc = 1,num_soilc
            c = filter_soilc(fc)
!            if(j .le. 40 .and. j .ge. 21 .and. c .eq. 1)print*,'matrix_Cinter',c,j,matrix_Cinter_next%V(c,j),matrix_Cinput_vector%V(c,j),matrix_Cinter%V(c,j)
            matrix_Cinter%V(c,j) = matrix_Cinput%V(c,j) + matrix_Cinter%V(c,j)
 
!            if(j .le. 20 .and. j .ge. 21 .and. c .eq. 1)print*,'matrix_Ninter',c,matrix_Ninter_next%V(c,j),matrix_Ninput_vector%V(c,j),matrix_Ninter%V(c,j)
!            if(c .eq. 1411 .and. j .ge. 21 .and. j .le. 30)print*,'Cinter before input',c,j,matrix_Ninter%V(c,j),matrix_Ninput%V(c,j)
            matrix_Ninter%V(c,j) = matrix_Ninput%V(c,j) + matrix_Ninter%V(c,j)
!            if(c .eq. 1411 .and. j .ge. 21 .and. j .le. 30)print*,'Cinter after input',c,j,matrix_Ninter%V(c,j),matrix_Ninput%V(c,j)
         end do
      end do
!      do j=1,ndecomp_pools_vr
!         do fc = 1,num_soilc
!            c = filter_soilc(fc)
!            print*,'matrix_Cinter%V here3.7',c,j,matrix_Cinter%V(c,j)
!         end do
!      end do

      if ( use_c13)then
         do j = 1, ndecomp_pools_vr
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               matrix_Cinter13%V(c,j) = matrix_Cinput13%V(c,j) + matrix_Cinter13%V(c,j)
            end do
         end do
      end if     
 
      if ( use_c14)then
         do j = 1, ndecomp_pools_vr
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               matrix_Cinter14%V(c,j) = matrix_Cinput14%V(c,j) + matrix_Cinter14%V(c,j)
            end do
         end do
      end if !c14

      call t_stopf('CN Soil matrix-matrix mult2-lev2')
!            hr_adjust = 0
            do i=1,ndecomp_pools
!               if(.not. is_cwd(i))then
               do j=1,nlevdecomp
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     hr(c) = hr(c) - (matrix_Cinter%V(c,j+(i-1)*nlevdecomp)-Cinter_old(c,j+(i-1)*nlevdecomp)&
                                    - matrix_Cinput%V(c,j+(i-1)*nlevdecomp))*dzsoi_decomp(j) / dt
                  end do
               end do
!               end if
            end do

!      if(bounds%begc .le. 12285 .and. bounds%endc .ge. 12285)print*,'hr',12285,hr(12285) * dt

!            print*,'filter_actfirec in CNSoilMatrix',num_actfirec,filter_actfirec(1:num_actfirec)
            if(num_actfirec .ne. 0)then
               do i=1,ndecomp_pools
                  do j=1,nlevdecomp
                     do fc = 1,num_actfirec
                        c = filter_actfirec(fc)
                        fire_delta = AKfiresoil%M(c,j+(i-1)*nlevdecomp)*Cinter_old(c,j+(i-1)*nlevdecomp)*dzsoi_decomp(j) / dt
!                           - matrix_decomp_fire_k(c,j+(i-1)*nlevdecomp)*Cinter_old(c,j+(i-1)*nlevdecomp)*dzsoi_decomp(j) / dt
!                           + Cinter_old(c,j+(i-1)*nlevdecomp)*dzsoi_decomp(j) 
                        m_decomp_cpools_to_fire(c,i) = m_decomp_cpools_to_fire(c,i) - fire_delta
                        hr(c) = hr(c) + fire_delta
                     end do
                  end do
               end do
            end if
!      print*,'hr',hr(c)
!            if(c .eq. 34208)print*,'matrix_Cinter_next1',er(c), hr(c)

!            hr(c) = hr_adjust

!            if(c .eq. 34208)print*,'matrix_Cinter_next',matrix_Cinter_next(c,1:nlevdecomp),er(c), hr(c)


     !print*,'here4'
            do i=1,ndecomp_pools_vr
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
               end do
            end do
      call t_startf('CN Soil matrix-assign back')

      do i=1,ndecomp_pools
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cs_soil%decomp_cpools_vr_col(c,j,i) = matrix_Cinter%V(c,j+(i-1)*nlevdecomp)
!               if(month .eq. 10)print*,'decomp_cpools',c,j,i,cs_soil%decomp_cpools_vr_col(c,j,i)
               ns_soil%decomp_npools_vr_col(c,j,i) = matrix_Ninter%V(c,j+(i-1)*nlevdecomp)
            end do
         end do
      end do
          
      if( use_c13 ) then
         do i=1,ndecomp_pools
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs13_soil%decomp_cpools_vr_col(c,j,i) = matrix_Cinter13%V(c,j+(i-1)*nlevdecomp)
               end do
            end do
         end do
      end if

      if( use_c14 ) then
         do i=1,ndecomp_pools
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs14_soil%decomp_cpools_vr_col(c,j,i) = matrix_Cinter14%V(c,j+(i-1)*nlevdecomp)
               end do
            end do
         end do
      end if

      call t_stopf('CN Soil matrix-assign back')

      if(use_soil_matrixcn .and. (is_outmatrix .or. isspinup))then
         call t_startf('CN Soil matrix-spinup & output1')
            
         do j=1,ndecomp_pools*nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cs_soil%in_acc (c,j) = cs_soil%in_acc (c,j) + matrix_Cinput%V(c,j)
               ns_soil%in_nacc(c,j) = ns_soil%in_nacc(c,j) + matrix_Ninput%V(c,j)
    !           if(j .eq. 121 .and. (c .eq. 2852 .or. c .eq. 2856 .or. c .eq. 2857 .or. c .eq. 2748 .or. c .eq. 2770 .or. c .eq. 2771))print*,'soil3n,input',c,ns_soil%in_nacc(c,j), matrix_Ninput%V(c,j),matrix_Ninter%V(c,j),col%wtgcell(c)
            end do
         end do

         call t_stopf('CN Soil matrix-spinup & output1')
!            do i = 1,ndecomp_pools
!               do j = 1,nlevdecomp
!                  do fc = 1,num_soilc
!                     c = filter_soilc(fc)
!                     if(j < nlevdecomp)then
!                        up_tran_rate = matrix_a_tri(c,j) / dzsoi_decomp(j)
!                        cs_soil%vert_up_tran_acc(c,j,i) = cs_soil%vert_up_tran_acc(c,j,i) &
!                             + tranvert(c,j+(i-1)*nlevdecomp,j+1+(i-1)*nlevdecomp) * matrix_Cinter(c,j+1+(i-1)*nlevdecomp) * dt
!                             + up_tran_rate * matrix_Cinter(c,j+1+(i-1)*nlevdecomp) * dt
!                        ns_soil%vert_up_tran_nacc(c,j,i) = ns_soil%vert_up_tran_nacc(c,j,i) &
!                             + ntranvert(c,j+(i-1)*nlevdecomp,j+1+(i-1)*nlevdecomp) * matrix_Ninter(c,j+1+(i-1)*nlevdecomp) * dt
!                             + up_tran_rate * matrix_Ninter(c,j+1+(i-1)*nlevdecomp) * dt
!                     end if
!                     if(j > 1)then
!                        down_tran_rate = matrix_c_tri(c,j) / dzsoi_decomp(j)
!                        cs_soil%vert_down_tran_acc(c,j,i) = cs_soil%vert_down_tran_acc(c,j,i) &
!                             + down_tran_rate * matrix_Cinter(c,j-1+(i-1)*nlevdecomp,1) * dt
!                             + tranvert(c,j+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp) * matrix_Cinter(c,j-1+(i-1)*nlevdecomp,1) * dt
!                        ns_soil%vert_down_tran_nacc(c,j,i) = ns_soil%vert_down_tran_nacc(c,j,i) &
!                             + down_tran_rate * matrix_Ninter(c,j-1+(i-1)*nlevdecomp,1) * dt
!                             + ntranvert(c,j+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp) * matrix_Ninter(c,j-1+(i-1)*nlevdecomp,1) * dt
!                     end if
                    
!                     cs_soil%exit_acc(c,j,i) = cs_soil%exit_acc(c,j,i) &
!                             + tranvert(c,j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp) * matrix_Cinter(c,j+(i-1)*nlevdecomp,1) * dt
!                     ns_soil%exit_nacc(c,j,i) = ns_soil%exit_nacc(c,j,i) &
!                             + ntranvert(c,j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp) * matrix_Ninter(c,j+(i-1)*nlevdecomp,1) * dt
!                  end do
!               end do
!            end do

!            do k = 1,ndecomp_cascade_transitions
!               if(cascade_receiver_pool(k) .ne. 0)then  !transition to atmosphere
!                  do j = 1,nlevdecomp
!                     do fc = 1,num_soilc
!                        c = filter_soilc(fc)
!                        cs_soil%hori_tran_acc(c,j,k) = cs_soil%hori_tran_acc(c,j,k) &
!                          + tranvert(c,j+(cascade_receiver_pool(k)-1)*nlevdecomp,j+(cascade_donor_pool(k)-1)*nlevdecomp) &
!                          * matrix_Cinter(c,j+(cascade_donor_pool(k)-1)*nlevdecomp,1) * dt
!                        ns_soil%hori_tran_nacc(c,j,k) = ns_soil%hori_tran_nacc(c,j,k) &
!                          + ntranvert(c,j+(cascade_receiver_pool(k)-1)*nlevdecomp,j+(cascade_donor_pool(k)-1)*nlevdecomp) &
!                          * matrix_Ninter(c,j+(cascade_donor_pool(k)-1)*nlevdecomp,1) * dt
!                     end do
!                  end do
!               end if
!            end do

     !print*,'here5'!,AKXcacc
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               Koned = matrix_Cinter(c,:)
     !print*,'here5.1'
     !print*,'here5.3'
         call Xdiagsoil%SetValueDM(begc,endc,num_soilc,filter_soilc,Ninter_old(begc:endc,1:ndecomp_pools_vr))
     !print*,'here5.5'
         call t_startf('CN Soil matrix-spinup & output1.5')
         call AKallsoiln%SPMM_AK(num_soilc,filter_soilc,Xdiagsoil)
         call t_stopf('CN Soil matrix-spinup & output1.5')

!            Koned = matrix_Ninter(c,:)
     !print*,'here5.4'
!         do j=1,ndecomp_pools_vr
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'matrix_Cinter%V',c,j,matrix_Cinter%V(c,j)
!            end do
!         end do
!         do fc = 1,num_soilc
!            c = filter_soilc(fc)
!            print*,'lat,lon',c,grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
!         end do
         call Xdiagsoil%SetValueDM(begc,endc,num_soilc,filter_soilc,Cinter_old(begc:endc,1:ndecomp_pools_vr))
     !print*,'here5.2'
         call t_startf('CN Soil matrix-spinup & output1.6')
!         do j=1,n_all_entries
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'here6.2429',c,j,all_i(j),all_j(j),AKallsoilc%M(c,j)
!            end do
!         end do
!         do j=1,ndecomp_pools_vr
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'Xdiagsoil',c,j,Xdiagsoil%DM(c,j)
!            end do
!         end do
         call AKallsoilc%SPMM_AK(num_soilc,filter_soilc,Xdiagsoil)
!         do j=1,n_all_entries
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               print*,'here6.243',c,j,all_i(j),all_j(j),AKallsoilc%M(c,j)
!            end do
!         end do
        
         call t_stopf('CN Soil matrix-spinup & output1.6')
!     print*,'here5.6'
!               index_ready2 = .True.
     !print*,'here5.6,AKXcacc',c,AKXcacc(c,:)
     !print*,'here5.6,AKXcinc',AKXcinc%M(1:AKXcinc%NE)
         call t_startf('CN Soil matrix-spinup & output2')

!        do j=1,n_all_entries
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               if(all_i(j) .eq. 121 .and. c .ge. 2852 .and. c .le. 2852)print*,'here6.244,2852',c,j,all_i(j),all_j(j),AKXnacc%M(c,j),AKallsoiln%M(c,j)
!               if(all_i(j) .eq. 121 .and. c .ge. 2856 .and. c .le. 2856)print*,'here6.244,2856',c,j,all_i(j),all_j(j),AKXnacc%M(c,j),AKallsoiln%M(c,j)
!               if(all_i(j) .eq. 121 .and. c .ge. 2857 .and. c .le. 2857)print*,'here6.244,2857',c,j,all_i(j),all_j(j),AKXnacc%M(c,j),AKallsoiln%M(c,j)
!               if(all_i(j) .eq. 121 .and. c .ge. 2748 .and. c .le. 2748)print*,'here6.244,2748',c,j,all_i(j),all_j(j),AKXnacc%M(c,j),AKallsoiln%M(c,j)
!               if(all_i(j) .eq. 121 .and. c .ge. 2770 .and. c .le. 2770)print*,'here6.244,2770',c,j,all_i(j),all_j(j),AKXnacc%M(c,j),AKallsoiln%M(c,j)
!               if(all_i(j) .eq. 121 .and. c .ge. 2771 .and. c .le. 2771)print*,'here6.244,2771',c,j,all_i(j),all_j(j),AKXnacc%M(c,j),AKallsoiln%M(c,j)
!            end do
!         end do
         call AKXnacc%SPMP_B_ACC(num_soilc,filter_soilc,AKallsoiln)
!         do j=1,AKallsoiln%NE
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               AKXnacc%M(c,j) = AKXnacc%M(c,j) + AKallsoiln%M(c,j)
!            end do
!         end do
         call t_stopf('CN Soil matrix-spinup & output2')
!     print*,'here5.7',c
!     print*,'here5.7',AKXnacc(c,:)
!     print*,'here5.7',AKXninc%M(1:AKXninc%NE)
         call t_startf('CN Soil matrix-spinup & output3')


         call AKXcacc%SPMP_B_ACC(num_soilc,filter_soilc,AKallsoilc)


!         do j=1,n_all_entries
!            if(bounds%begc .le. 7856 .and. bounds%endc .ge. 7873 .and. all_j(j) .eq. 132)print*,'AKXcacc1',j,AKXcacc%M(7856,j),AKXcacc%M(7861,j),AKXcacc%M(7867,j),AKXcacc%M(7869,j),AKXcacc%M(7871,j),AKXcacc%M(7873,j)
!            if(bounds%begc .le. 7856 .and. bounds%endc .ge. 7873 .and. all_i(j) .eq. 132)print*,'AKXcacc2',j,AKXcacc%M(7856,j),AKXcacc%M(7861,j),AKXcacc%M(7867,j),AKXcacc%M(7869,j),AKXcacc%M(7871,j),AKXcacc%M(7873,j)
!         end do
!         print*,'hereNE',AKXcacc%NE,AKallsoilc%NE
!         print*,'hereRI1',AKXcacc%RI
!         print*,'hereRI2',AKallsoilc%RI
!         print*,'hereCI1',AKXcacc%CI
!         print*,'hereCI2',AKallsoilc%CI
!         do j=1,AKallsoilc%NE
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               AKXcacc%M(c,j) = AKXcacc%M(c,j) + AKallsoilc%M(c,j)
!            end do
!         end do
         call t_stopf('CN Soil matrix-spinup & output3')
     !print*,'here5.8'
                  
!     print*,'here6',is_end_curr_year(),ndecomp_pools_vr!,is_end_curr_year()
!            do i=1,ndecomp_pools_vr
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  print*,'here6.22',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!               end do
!            end do
!         do fc = 1,num_soilc
!            c = filter_soilc(fc)
!            if(bounds%begc .le. 6851 .and. bounds%endc .ge. 6851)print*,'c',c,grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
!            if(c .ge. 7856 .and. c .le. 7873)print*,'in cap update',c, cs_soil%matrix_cap_decomp_cpools_vr_col(c,12,7)
         !end do
!            do i=1,ndecomp_pools
!               do j = 1,nlevdecomp
!                  do fc = 1,num_soilc
!                     c = filter_soilc(fc)
!                     if(j .le. 3 .and. i .eq. 7 .and. (c .eq. 2852 .or. c .eq. 2856 .or. c .eq. 2857 .or. c .eq. 2748 .or. c .eq. 2770 .or. c .eq. 2771))print*,'soil3n,cap,each step',c,j,i,ns_soil%matrix_cap_decomp_npools_vr_col(c,j,i),col%wtgcell(c)
!                  end do
!               enddo
!            end do
         call t_startf('CN Soil matrix-calc. C capacity')
         if(is_end_curr_year())then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cs_soil%tran_acc (c,1:ndecomp_pools_vr,1:ndecomp_pools_vr) = 0._r8
               ns_soil%tran_nacc(c,1:ndecomp_pools_vr,1:ndecomp_pools_vr) = 0._r8
            end do
!            if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)then
!               c = 115570
!               if(abs(ns_soil%tran_nacc(c,113,113)) .lt. epsi)then
!                  tmptmp=1._r8
!               end if
!            end if
            !print*,'here6.1'
!            do i=1,ndecomp_pools_vr
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  if(begc .eq. 2277)print*,'here6.23',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!                  print*,'here6.23',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!               end do
!            end do
!            if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)then
!               c = 115570
!               if(abs(ns_soil%tran_nacc(c,113,113)) .lt. epsi)then
!                  tmptmp=1._r8
!               end if
!            end if
!            do i = 1,ndecomp_pools
!               do j = 1,nlevdecomp
!                  do fc = 1,num_soilc
!                     c = filter_soilc(fc)
!                     cs_soil%in_acc(c,j+(i-1)*nlevdecomp)   = cs_soil%in_acc_2d(c,j,i) 
!                     ns_soil%in_nacc(c,j+(i-1)*nlevdecomp)  = ns_soil%in_nacc_2d(c,j,i) 
!                     if(j < nlevdecomp)then
!                           cs_soil%tran_acc(c,j+(i-1)*nlevdecomp,j+1+(i-1)*nlevdecomp) = cs_soil%vert_up_tran_acc(c,j,i) &
!                                   / cs_soil%decomp0_cpools_vr_col(c,j+1,i)
!                           ns_soil%tran_nacc(c,j+(i-1)*nlevdecomp,j+1+(i-1)*nlevdecomp) = ns_soil%vert_up_tran_nacc(c,j,i) &
!                                   / ns_soil%decomp0_npools_vr_col(c,j+1,i)
!                        end if
!                        if(j > 1)then
!                           cs_soil%tran_acc(c,j+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp) = cs_soil%vert_down_tran_acc(c,j,i) &
!                                   / cs_soil%decomp0_cpools_vr_col(c,j-1,i)
!                           ns_soil%tran_nacc(c,j+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp) = ns_soil%vert_down_tran_nacc(c,j,i) &
!                                   / ns_soil%decomp0_npools_vr_col(c,j-1,i)
!                        end if
!                        cs_soil%tran_acc(c,j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp) = cs_soil%exit_acc(c,j,i) &
!                                      / cs_soil%decomp0_cpools_vr_col(c,j,i)
!                        ns_soil%tran_nacc(c,j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp) = ns_soil%exit_nacc(c,j,i) & 
!                                      / ns_soil%decomp0_npools_vr_col(c,j,i)
!                  end do
!               end do
!            end do
!            do i=1,ndecomp_pools_vr
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
                  !if(begc .eq. 2277)print*,'here6.24',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!                  print*,'here6.24',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!               end do
!            end do
!         print*,'before here1',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)!,cs_soil%decomp0_cpools_vr_col(:,9,10)
            do j=1,n_all_entries
               j_lev    = mod(all_j(j)-1,nlevdecomp)+1
               j_decomp = (all_j(j) - j_lev)/nlevdecomp + 1
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                  print*,'tran_acc set',c,j,all_i(j),all_j(j),j_lev,j_decomp,AKXcacc%M(c,j),cs_soil%decomp0_cpools_vr_col(c,j_lev,j_decomp)
                  cs_soil%tran_acc(c,all_i(j),all_j(j))  = AKXcacc%M(c,j) / cs_soil%decomp0_cpools_vr_col(c,j_lev,j_decomp)
                  !print*,'aftertran_acc',cs_soil%tran_acc(c,all_i(j),all_j(j))
                  !print*,'tran_nacc set',c,j,all_i(j),all_j(j),j_lev,j_decomp,AKXnacc%M(c,j),ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
                  ns_soil%tran_nacc(c,all_i(j),all_j(j)) = AKXnacc%M(c,j) / ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
                  !print*,'aftertran_nacc',ns_soil%tran_nacc(c,all_i(j),all_j(j))
!                  print*,'here6.245',c,j,cs_soil%tran_acc(c,all_i(j),all_j(j)),all_i(j),all_j(j),AKXcacc%M(c,j)
!                  print*,'here6.1',ns_soil%tran_nacc(1,120,120)
!                   if(all_i(j) .eq. 132 .and. c .ge. 7856 .and. c .le. 7873)print*,'in cap2',c, j,AKXcacc%M(c,j),all_i(j),all_j(j),cs_soil%decomp0_cpools_vr_col(c,j_lev,j_decomp)
!                   if(all_i(j) .eq. 121 .and. c .ge. 2852 .and. c .le. 2852)print*,'in cap1,2852',c, j,AKXnacc%M(c,j),all_i(j),all_j(j),ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
!                   if(all_i(j) .eq. 121 .and. c .ge. 2856 .and. c .le. 2856)print*,'in cap1,2856',c, j,AKXnacc%M(c,j),all_i(j),all_j(j),ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
!                   if(all_i(j) .eq. 121 .and. c .ge. 2857 .and. c .le. 2857)print*,'in cap1,2857',c, j,AKXnacc%M(c,j),all_i(j),all_j(j),ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
!                   if(all_i(j) .eq. 121 .and. c .ge. 2748 .and. c .le. 2748)print*,'in cap1,2748',c, j,AKXnacc%M(c,j),all_i(j),all_j(j),ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
!                   if(all_i(j) .eq. 121 .and. c .ge. 2770 .and. c .le. 2770)print*,'in cap1,2770',c, j,AKXnacc%M(c,j),all_i(j),all_j(j),ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
!                   if(all_i(j) .eq. 121 .and. c .ge. 2771 .and. c .le. 2771)print*,'in cap1,2771',c, j,AKXnacc%M(c,j),all_i(j),all_j(j),ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
               end do
            end do

!            if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)then
!               c = 115570
!               if(abs(ns_soil%tran_nacc(c,113,113)) .le. epsi)then
!                  tmptmp=1._r8
!               end if
!            end if
!               print*,'here6.2'!,ns_soil%tran_nacc(2,115,115)
!               do k = 1,ndecomp_cascade_transitions
!                  do j = 1,nlevdecomp
!                     do fc = 1,num_soilc
!                        c = filter_soilc(fc)
!                        if(cascade_receiver_pool(k) .ne. 0)then
!                           cs_soil%tran_acc(c,j+(cascade_receiver_pool(k)-1)*nlevdecomp,j+(cascade_donor_pool(k)-1)*nlevdecomp) &
!                                                 = cs_soil%hori_tran_acc(c,j,k) / cs_soil%decomp0_cpools_vr_col(c,j,cascade_donor_pool(k))
!                           ns_soil%tran_nacc(c,j+(cascade_receiver_pool(k)-1)*nlevdecomp,j+(cascade_donor_pool(k)-1)*nlevdecomp) &
!                                                 = ns_soil%hori_tran_nacc(c,j,k) / ns_soil%decomp0_npools_vr_col(c,j,cascade_donor_pool(k))
!                        end if
!                     end do
!                  end do
!               end do
         
            
            !print*,'here6.2'
!            do i=1,ndecomp_pools_vr
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
                  !if(begc .eq. 2277)print*,'here6.25',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!                  print*,'here6.25',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
!               end do
!            end do
!            if(begc .eq. 2277)then
!               if(abs(cs_soil%tran_acc(2282,132,132) .le. epsi))then
!                  cs_soil%tran_acc(c,i,i) = 1.e+36_r8
!               end if
!               print*,'2277',abs(cs_soil%tran_acc(2282,132,132)
!            end if
            do i=1,ndecomp_pools_vr
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                     print*,'here6.26',fc,begc,endc,c,i,abs(cs_soil%tran_acc(c,i,i))
                  if (abs(cs_soil%tran_acc(c,i,i)) .le. epsi)then !avoid inversion nan
                      cs_soil%tran_acc(c,i,i) = 1.e+36_r8
                  end if 
!                  if (abs(ns_soil%tran_nacc(c,i,i)) .le. epsi)then
!                     ns_soil%tran_nacc(c,i,i) = 1.e+36
!                  end if 
               end do
            end do

            !print*,'here6.3'
!            if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)then
!               c = 115570
!               if(abs(ns_soil%tran_nacc(c,113,113)) .le. epsi)then
!                  tmptmp=1._r8
!               end if
!            end if
            !print*,'here6.4'
            do i=1,ndecomp_pools_vr
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
!                  if(bounds%begc .le. 26034 .and. bounds%endc .ge. 26034)print*,'tran_nacc',c,i,epsi
!                  if(bounds%begc .le. 26034 .and. bounds%endc .ge. 26034)print*,ns_soil%tran_nacc(c,i,i)
!                  print*,ns_soil%tran_nacc(1,120,120)
!                  if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,'c,i,i',c,i,i
!                  if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,ns_soil%tran_nacc(c,i,i)
                  if (abs(ns_soil%tran_nacc(c,i,i)) .le. epsi)then
!                     if(bounds%begc .le. 26034 .and. bounds%endc .ge. 26034)print*,'hi',ns_soil%tran_nacc(c,i,i)
!                     if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,'hi',c,bounds%begc,bounds%endc,i
!                     if(bounds%begc .le. 115568 .and. bounds%endc .ge. 115568)print*,ns_soil%tran_nacc(c,i,i)
                     ns_soil%tran_nacc(c,i,i) = 1.e+36_r8
                  end if 
               end do
            end do
!         print*,'here2 decomp0',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)
!            print*,ns_soil%tran_nacc(1,120,120)
!    print*,'here6.3'!,ns_soil%tran_nacc(2,115,115)
     !print*,'here7'

            do fc = 1,num_soilc
               c = filter_soilc(fc)
!               print*,'inacc',cs_soil%in_acc(c,1:ndecomp_pools_vr)
!               do i=1,140
!                  print*,'tranacc',i,cs_soil%tran_acc(c,i,1:ndecomp_pools_vr)
!               end do
               call inverse(cs_soil%tran_acc(c,1:ndecomp_pools_vr,1:ndecomp_pools_vr),AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
               call inverse(ns_soil%tran_nacc(c,1:ndecomp_pools_vr,1:ndecomp_pools_vr),AKinvn(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
               soilmatrixc_cap(c,:,1) = -matmul(AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),cs_soil%in_acc(c,1:ndecomp_pools_vr))
               soilmatrixn_cap(c,:,1) = -matmul(AKinvn(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ns_soil%in_nacc(c,1:ndecomp_pools_vr))
            end do
         
!            print*,'after cap decomp0',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)
    !print*,'here8'
            do i=1,ndecomp_pools
               do j = 1,nlevdecomp
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                  !   if(soilmatrixc_cap(j+(i-1)*nlevdecomp,1) .le. 1.e-25)then ! for model stability
                  !      soilmatrixc_cap(j+(i-1)*nlevdecomp,1) = cs_soil%decomp_cpools_vr_col(c,j,i)
                  !      soilmatrixn_cap(j+(i-1)*nlevdecomp,1) = ns_soil%decomp_npools_vr_col(c,j,i)
                  !   end if
                     if(isspinup)then
                        cs_soil%decomp_cpools_vr_col(c,j,i) =  soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1)
                        ns_soil%decomp_npools_vr_col(c,j,i) =  soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1)
                     end if
!                    print*,'in cap1 decomp0',i,j,c,cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)
                     cs_soil%matrix_cap_decomp_cpools_vr_col(c,j,i) = soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1)
!                    print*,'in cap2 decomp0',i,j,c,cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)
                     ns_soil%matrix_cap_decomp_npools_vr_col(c,j,i) = soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1)
!                     if(j .eq. 12 .and. i .eq. 7 .and. c .le. 7856 .and. c .ge. 7873)print*,'in capcal',c, cs_soil%matrix_cap_decomp_cpools_vr_col(c,j,i)
!                     cs_soil%matrix_pot_decomp_cpools_vr_col(c,j,i) = soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1) - cs_soil%decomp_cpools_vr_col(c,j,i)
!                     ns_soil%matrix_pot_decomp_npools_vr_col(c,j,i) = soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1) - ns_soil%decomp_npools_vr_col(c,j,i)
!                     if(j .eq. 1 .and. i .eq. 7 .and. (c .eq. 2852 .or. c .eq. 2856 .or. c .eq. 2857 .or. c .eq. 2748 .or. c .eq. 2770 .or. c .eq. 2771))print*,'soil3n,cap',c,j,i,ns_soil%matrix_cap_decomp_npools_vr_col(c,j,i),col%wtgcell(c)
                  end do
               end do
            end do
!            print*,'after cap1 decomp0',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)
            do j=1,n_all_entries
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  AKXcacc%M(c,j) = 0._r8
                  AKXnacc%M(c,j) = 0._r8
               end do
            end do

!            print*,'before in decomp0',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cs_soil%in_acc   (c,:)   = 0._r8
               ns_soil%in_nacc  (c,:)   = 0._r8
            end do
!               cs_soil%vert_up_tran_acc(c,:,:)   = 0.0_r8 
!               cs_soil%vert_down_tran_acc(c,:,:) = 0.0_r8 
!               cs_soil%exit_acc(c,:,:)           = 0.0_r8 
!               cs_soil%hori_tran_acc(c,:,:)      = 0.0_r8 

!               ns_soil%vert_up_tran_nacc(c,:,:)   = 0.0_r8 
!               ns_soil%vert_down_tran_nacc(c,:,:) = 0.0_r8 
!               ns_soil%exit_nacc(c,:,:)           = 0.0_r8 
!               ns_soil%hori_tran_nacc(c,:,:)      = 0.0_r8 
         end if
         !print*,'here9'
         call t_stopf('CN Soil matrix-calc. C capacity')
      end if !is out_matrix
!         print*,'end decomp0',cs_soil%decomp0_cpools_vr_col(:,10,1),cs_soil%decomp_cpools_vr_col(:,10,1)

   end associate 
 end subroutine CNSoilMatrix
 
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

end module CNSoilMatrixMod

