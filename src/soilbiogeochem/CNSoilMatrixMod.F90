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
  ! Xn+1 = Xn + I*dt + (A*K(ksi) - Kfire - tri/dz)*Xn*dt
  ! Or
  ! Xn+1 = Xn + I*dt + (A*K(ksi) - Kfire - V)*Xn*dt
  
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use decompMod                          , only : bounds_type  
  use abortutils                         , only : endrun
  use spmdMod                            , only : masterproc
  use clm_time_manager                   , only : get_step_size, is_end_curr_month,get_curr_date,update_DA_nstep
  use clm_time_manager                   , only : is_first_restart_step,is_beg_curr_year,is_end_curr_year,is_first_step_of_this_run_segment
  use clm_varpar                         , only : ndecomp_pools, nlevdecomp, ndecomp_pools_vr        !number of biogeochemically active soil layers
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_cascade_outtransitions
  use clm_varpar                         , only : i_cwd
  use clm_varcon                         , only : dzsoi_decomp,zsoi,secspday,c3_r2,c14ratio
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, use_soil_matrixcn
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use CNVegNitrogenFluxType              , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType    , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType     , only : soilbiogeochem_nitrogenflux_type  
  use CNSharedParamsMod                  , only : CNParamsShareInst
  use SoilStateType                      , only : soilstate_type  
  use clm_varctl                         , only : spinup_matrixcn, hist_wrt_matrixcn_diag, nyr_forcing, nyr_SASU, iloop_avg
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use clm_varctl                         , only : use_c13, use_c14, iulog
  use perf_mod                           , only : t_startf, t_stopf
  use SparseMatrixMultiplyMod            , only : sparse_matrix_type, diag_matrix_type, vector_type
  use MatrixMod                          , only : inverse
!
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNSoilMatrixInit        ! Initialization for CN Soil Matrix solution
  public:: CNSoilMatrix
  public:: CNSoilMatrixRest        ! Restart for CN Soil Matrix solution

  ! ! PRIVATE MEMBER DATA:
  integer,save, private :: iyr=0   ! Cycling year number into forcing sequence
  integer,save, private :: iloop=0 ! The iloop^th forcing loop
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNSoilMatrixInit( )
    ! !DESCRIPTION: Initialization for CN soil Matrix solution
    ! !ARGUMENTS:
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    if ( use_soil_matrixcn .and. masterproc) then
       write(iulog,*) 'CN Soil matrix solution is on'
       write(iulog,*) '*****************************'
       if ( spinup_matrixcn ) then
          write(iulog,*) '   Matrix spinup is on'
          write(iulog,*) '   *******************'
          write(iulog,*) '   nyr_forcing = ', nyr_forcing
          write(iulog,*) '   nyr_SASU    = ', nyr_SASU
          write(iulog,*) '   iloop_avg   = ', iloop_avg
       end if
       if ( hist_wrt_matrixcn_diag )then
          write(iulog,*) '   Extra matrix solution tracability output is turned on'
       else
          write(iulog,*) '   no extra matrix solution tracability output'
       end if
    end if
  end subroutine CNSoilMatrixInit

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
    real(r8):: dt                   ! time step (seconds)
    real(r8):: epsi,fire_delta      ! small number

    integer :: begc,endc                                    ! bounds 
    real(r8),dimension(bounds%begc:bounds%endc,nlevdecomp*(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)) :: a_ma_vr,na_ma_vr

    real(r8),dimension(bounds%begc:bounds%endc,1:ndecomp_pools_vr,1) ::  soilmatrixc_cap,soilmatrixc13_cap,soilmatrixc14_cap,soilmatrixn_cap
    real(r8), dimension(1:ndecomp_pools_vr,1:ndecomp_pools_vr)   ::  AKinv,AKinvn

    real(r8),dimension(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools) :: cn_decomp_pools
    integer :: Ntrans
    integer tranlist_a
    integer j_decomp,j_lev,ilev,idecomp
    integer,dimension(:) :: kfire_i(1:ndecomp_pools_vr)
    integer,dimension(:) :: kfire_j(1:ndecomp_pools_vr)
    real(r8),dimension(:,:) :: Cinter_old(bounds%begc:bounds%endc,1:ndecomp_pools_vr)
    real(r8),dimension(:,:) :: C13inter_old(bounds%begc:bounds%endc,1:ndecomp_pools_vr)
    real(r8),dimension(:,:) :: C14inter_old(bounds%begc:bounds%endc,1:ndecomp_pools_vr)
    real(r8),dimension(:,:) :: Ninter_old(bounds%begc:bounds%endc,1:ndecomp_pools_vr)
    logical,save :: list_ready1_fire   = .False.
    logical,save :: list_ready1_nofire = .False.
    logical,save :: list_ready2_fire   = .False.
    logical,save :: list_ready2_nofire = .False.
    logical,save :: list_ready3_fire   = .False.
    logical,save :: init_readyAsoilc = .False.
    logical,save :: init_readyAsoiln = .False.
    logical isbegofyear

    !-----------------------------------------------------------------------
    begc = bounds%begc; endc = bounds%endc
    
!    SHR_ASSERT_ALL((ubound(cn_decomp_pools)     == (/endc,nlevdecomp,ndecomp_pools/))  , errMsg(sourcefile, __LINE__))
    associate(                                      &
    cs_soil   => soilbiogeochem_carbonstate_inst    , & ! In/Output
    ns_soil   => soilbiogeochem_nitrogenstate_inst  , & ! In/Output
    cs13_soil => c13_soilbiogeochem_carbonstate_inst, & ! In/Output 
    cs14_soil => c14_soilbiogeochem_carbonstate_inst, & ! In/Output 
    cf13_soil => c13_soilbiogeochem_carbonflux_inst,  & ! In/Output 
    cf14_soil => c14_soilbiogeochem_carbonflux_inst,  & ! In/Output 

    fpi_vr                        => soilbiogeochem_state_inst%fpi_vr_col                 ,&!Input:[real(r8)(:,:)]fraction of potential immobilization (no units)
    cascade_donor_pool            => decomp_cascade_con%cascade_donor_pool                ,&!Input:[integer(:)]which pool is C taken from for a given decomposition step
    cascade_receiver_pool         => decomp_cascade_con%cascade_receiver_pool             ,&!Input:[integer(:)]which pool is C added to for a given decomposition step 
    floating_cn_ratio_decomp_pools=> decomp_cascade_con%floating_cn_ratio_decomp_pools    ,&!Input:[logical(:)]TRUE => pool has fixed C:N ratio
    initial_cn_ratio              => decomp_cascade_con%initial_cn_ratio                  ,&!Input:[real(r8)(:)]c:n ratio for initialization of pools
    rf_decomp_cascade             => soilbiogeochem_carbonflux_inst%rf_decomp_cascade_col      ,&!Input:[real(r8)(:,:,:)]respired fraction in decomposition step (frac)
    pathfrac_decomp_cascade       => soilbiogeochem_carbonflux_inst%pathfrac_decomp_cascade_col,&!Input:[real(r8)(:,:,:)]what fraction of C leaving a given pool passes 
                                                                                            !                       through a given transition (frac)
    is_cwd                        => decomp_cascade_con%is_cwd                            ,&!Input:[logical(:)]TRUE => pool is a cwd pool
    is_litter                     => decomp_cascade_con%is_litter                         ,&!Input:[logical(:)]TRUE => pool is a litter pool

    hr                         => soilbiogeochem_carbonflux_inst%hr_col                   ,&!Output:[real(r8)(:)]heterotrophic respiration
    trcr_ctendency             => soilbiogeochem_carbonflux_inst%decomp_cpools_transport_tendency_col,&
    trcr_ntendency             => soilbiogeochem_nitrogenflux_inst%decomp_npools_transport_tendency_col,&
    m_decomp_cpools_to_fire    => cnveg_carbonflux_inst%m_decomp_cpools_to_fire_col       ,&!Output:[real(r8)(:,:)]vertically-integrated decomposing C fire loss
    tri_ma_vr                  => soilbiogeochem_carbonflux_inst%tri_ma_vr                ,&!Input:[real(r8)(:,:)]vertical C transfer rate in sparse matrix format (gC*m3)/(gC*m3*step))
    matrix_decomp_fire_k       => soilbiogeochem_carbonflux_inst%matrix_decomp_fire_k_col ,&!Input:[real(r8)(:,:)]decomposition rate due to fire (gC*m3)/(gC*m3*step))

    AKsoilc        => soilbiogeochem_carbonflux_inst%AKsoilc   ,&!Output:[SparseMatrix] A*K for C transfers between pools
    RI_a           => soilbiogeochem_carbonflux_inst%RI_a      ,&!In/Output:[Integer(:)] Row numbers of all entries from AKsoilc, Automatically generated by SetValueA
    CI_a           => soilbiogeochem_carbonflux_inst%CI_a      ,&!In/Output:[Integer(:)] Column numbers of all entries from AKsoilc, Automatically generated by SetValueA
    AKsoiln        => soilbiogeochem_nitrogenflux_inst%AKsoiln ,&!Output:[SparseMatrix] A*K for N transfers between pools
    RI_na          => soilbiogeochem_nitrogenflux_inst%RI_na   ,&!In/Output:[Integer(:)] Row numbers of all entries from AKsoiln, Automatically generated by SetValueA
    CI_na          => soilbiogeochem_nitrogenflux_inst%CI_na   ,&!In/Output:[Integer(:)] Column numbers of all entries from AKsoiln, Automatically generated by SetValueA

    A_i            => decomp_cascade_con%A_i                   ,&!Input:[integer(:)] Prescribed row number of all elements in a_ma_vr
    A_j            => decomp_cascade_con%A_j                   ,&!Input:[integer(:)] Prescribed column number of all elements in na_ma_vr
    spm_tranlist_a => decomp_cascade_con%spm_tranlist_a        ,&!Input:[integer(:,:)] Prescribed subscripts to map 2D variables (transitions,soil layer) to 1D sparse matrix format in a_ma_vr and na_ma_vr

    AVsoil         => soilbiogeochem_carbonflux_inst%AVsoil    ,&!Output:[SparseMatrix] V for C and N transfers between soil layers
    tri_i          => decomp_cascade_con%tri_i                 ,&!Input:[integer(:)] Prescribed row index of all entries in AVsoil
    tri_j          => decomp_cascade_con%tri_j                 ,&!Input:[integer(:)] Prescribed column index of all entries in AVsoil
    Ntri_setup     => decomp_cascade_con%Ntri_setup            ,&!Input:[integer] Number of non-zero entries in AVsoil

    AKfiresoil     => soilbiogeochem_carbonflux_inst%AKfiresoil,&!Output:[SparseMatrix] Kfire for CN transfers from soil to atm due to fire

    AKallsoilc     => soilbiogeochem_carbonflux_inst%AKallsoilc     ,&!Output:[SparseMatrix] (A*K+V-Kfire) for soil C cycle
    NE_AKallsoilc  => soilbiogeochem_carbonflux_inst%NE_AKallsoilc  ,&!In/Output:[Integer] Number of entries in AKallsoilc, Automatically generated by functions SPMP_*
    RI_AKallsoilc  => soilbiogeochem_carbonflux_inst%RI_AKallsoilc  ,&!In/Output:[Integer(:)] Row numbers of entries in AKallsoilc, Automatically generated by functions SPMP_*
    CI_AKallsoilc  => soilbiogeochem_carbonflux_inst%CI_AKallsoilc  ,&!In/Output:[Integer(:)] Column numbers of entries in AKallsoilc, Automatically generated by functions SPMP_*
    AKallsoiln     => soilbiogeochem_nitrogenflux_inst%AKallsoiln   ,&!Output:[SparseMatrix] (A*K+V-Kfire) for soil N cycle
    NE_AKallsoiln  => soilbiogeochem_nitrogenflux_inst%NE_AKallsoiln,&!In/Output:[Integer] Number of entries in AKallsoilc, Automatically generated by functions SPMP_*
    RI_AKallsoiln  => soilbiogeochem_nitrogenflux_inst%RI_AKallsoiln,&!In/Output:[Integer(:)] Row numbers of entries in AKallsoilc, Automatically generated by functions SPMP_*
    CI_AKallsoiln  => soilbiogeochem_nitrogenflux_inst%CI_AKallsoiln,&!In/Output:[Integer(:)] Column numbers of entries in AKallsoilc, Automatically generated by functions SPMP_*
    AKXcacc        => soilbiogeochem_carbonstate_inst%AKXcacc       ,&!In/Output:[SparseMatrix] Accumulated transfers for soil C cycle
    AKXnacc        => soilbiogeochem_nitrogenstate_inst%AKXnacc     ,&!In/Output:[SparseMatrix] Accumulated transfers for soil N cycle
    n_all_entries  => decomp_cascade_con%n_all_entries              ,&!Input:[integer] Number of all entries in AKallsoilc, AKallsoiln, AKXcacc, and AKXnacc
    all_i          => decomp_cascade_con%all_i                      ,&!Input:[integer(:)] Prescribed row index of all entries in AKallsoilc, AKallsoiln, AKXcacc, and AKXnacc
    all_j          => decomp_cascade_con%all_j                      ,&!Input:[integer(:)] Prescribed column index of all entries in AKallsoilc, AKallsoiln, AKXcacc, and AKXnacc

    Ksoil          => soilbiogeochem_carbonflux_inst%Ksoil             ,&!Output:[DiagonalMatrix] C turnover rate in different soil pools and layers
    Ksoiln         => soilbiogeochem_nitrogenflux_inst%Ksoiln          ,&!Output:[DiagonalMatrix] N turnover rate in different soil pools and layers
    Xdiagsoil      => soilbiogeochem_carbonflux_inst%Xdiagsoil         ,&!Output:[DiagonalMatrix] Temporary C and N state variable to calculate accumulation transfers
    matrix_Cinter  => soilbiogeochem_carbonstate_inst%matrix_Cinter    ,&!In/Output:[Vector] Soil C state variables (gC/m3) in different soil pools and layers
    matrix_Ninter  => soilbiogeochem_nitrogenstate_inst%matrix_Ninter  ,&!In/Output:[Vector] Soil N state variables (gN/m3) in different soil pools and layers
    matrix_Cinter13=> c13_soilbiogeochem_carbonstate_inst%matrix_Cinter,&!In/Output:[Vector] Soil C13 state variables (gC13/m3) in different soil pools and layers
    matrix_Cinter14=> c14_soilbiogeochem_carbonstate_inst%matrix_Cinter,&!In/Output:[Vector] Soil C14 state variables (gC14/m3) in different soil pools and layers
    matrix_Cinput  => soilbiogeochem_carbonflux_inst%matrix_Cinput     ,&!Input:[Vector] C input to different soil compartments (pools and layers) (gC/m3/step)
    matrix_Cinput13=> c13_soilbiogeochem_carbonflux_inst%matrix_Cinput ,&!Input:[Vector] C13 input to different soil compartments (pools and layers) (gC13/m3/step)
    matrix_Cinput14=> c14_soilbiogeochem_carbonflux_inst%matrix_Cinput ,&!Input:[Vector] C14 input to different soil compartments (pools and layers) (gC14/m3/step)
    matrix_Ninput  => soilbiogeochem_nitrogenflux_inst%matrix_Ninput   ,&!Input:[Vector] N input to different soil compartments (pools and layers) (gN/m3/step)

    list_Asoilc      => decomp_cascade_con%list_Asoilc      ,&!In/Output:[Integer(:)] Saves mapping indices from a_ma_vr to AKsoilc
    list_Asoiln      => decomp_cascade_con%list_Asoiln      ,&!In/Output:[Integer(:)] Saves mapping indices from na_ma_vr to AKsoiln
    list_V_AKVfire   => decomp_cascade_con%list_V_AKVfire   ,&!In/Output:[Integer(:)] Saves mapping indices from V to (A*K+V-Kfire) in the addition subroutine SPMP_ABC
    list_fire_AKVfire=> decomp_cascade_con%list_fire_AKVfire,&!In/Output:[Integer(:)] Saves mapping indices from Kfire to (A*K+V-Kfire) in the addition subroutine SPMP_ABC
    list_AK_AKVfire  => decomp_cascade_con%list_AK_AKVfire  ,&!In/Output:[Integer(:)] Saves mapping indices from A*K to (A*K+V-Kfire) in the addition subroutine SPMP_ABC
    list_AK_AKV      => decomp_cascade_con%list_AK_AKV      ,&!In/Output:[Integer(:)] Saves mapping indices from A*K to (A*K+V) in the addition subroutine SPMP_AB
    list_V_AKV       => decomp_cascade_con%list_V_AKV        &!In/Output:[Integer(:)] Saves mapping indices from V to (A*K+V) in the addition subroutine SPMP_AB 
    )

     ! set time steps
      call t_startf('CN Soil matrix-init. matrix')
      dt = real( get_step_size(), r8 )

      Ntrans = (ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp
      epsi = 1.e-8_r8 

      isbegofyear = is_beg_curr_year()

      ! calculate c:n ratios of applicable pools
      do l = 1, ndecomp_pools
         if ( floating_cn_ratio_decomp_pools(l)) then
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if ( ns_soil%decomp_npools_vr_col(c,j,l) > 0._r8 ) then
                     cn_decomp_pools(c,j,l) = cs_soil%decomp_cpools_vr_col(c,j,l) / ns_soil%decomp_npools_vr_col(c,j,l)
                  else
                     cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
                  end if
               end do
            end do
         else
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
               end do
            end do
         end if 
      end do

      call t_stopf('CN Soil matrix-init. matrix')

      call t_startf('CN Soil matrix-assign matrix-a-na')
      ! Calculate non-diagonal entries (a_ma_vr and na_ma_vr) in transfer coefficient matrix A 
      do k = 1, ndecomp_cascade_transitions
         if(cascade_receiver_pool(k) .ne. 0)then  !transition to atmosphere
            do j = 1, nlevdecomp
               tranlist_a = spm_tranlist_a(j,k)
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  a_ma_vr(c,tranlist_a) = (1.0_r8-rf_decomp_cascade(c,j,k))*pathfrac_decomp_cascade(c,j,k)
                  if( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)))then
                     na_ma_vr(c,tranlist_a) = (1.0_r8-rf_decomp_cascade(c,j,k))* &
                               (cn_decomp_pools(c,j,cascade_donor_pool(k))/cn_decomp_pools(c,j,cascade_receiver_pool(k)))*pathfrac_decomp_cascade(c,j,k)
                  else
                     na_ma_vr(c,tranlist_a) = pathfrac_decomp_cascade(c,j,k)
                  end if
               end do
            end do
         end if
      end do

      call t_stopf('CN Soil matrix-assign matrix-a-na')

      ! Update the turnover rate matrix K with N limitation (fpi_vr)

      ! Assign old value to vector, and be ready for matrix operation
      call t_startf('CN Soil matrix-assign matrix-inter')
      do i = 1,ndecomp_pools
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
                  matrix_Cinter%V(c,j+(i-1)*nlevdecomp)   = cs_soil%decomp_cpools_vr_col(c,j,i)                
                  Cinter_old(c,j+(i-1)*nlevdecomp)        = cs_soil%decomp_cpools_vr_col(c,j,i) 
                  !Cinter_old is saved for accumulation of C transfer calculation and C flux (hr and fire) adjustment
                  matrix_Ninter%V(c,j+(i-1)*nlevdecomp)   = ns_soil%decomp_npools_vr_col(c,j,i)     
                  Ninter_old(c,j+(i-1)*nlevdecomp)        = ns_soil%decomp_npools_vr_col(c,j,i)
                  !Ninter_old is saved for accumulation of N transfer calculation
            end do
         end do
      end do
      if ( use_c13 )then
         do i = 1,ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  matrix_Cinter13%V(c,j+(i-1)*nlevdecomp)   = cs13_soil%decomp_cpools_vr_col(c,j,i)
                  C13inter_old(c,j+(i-1)*nlevdecomp)        = cs13_soil%decomp_cpools_vr_col(c,j,i)
               end do
            end do
         end do
      end if !c13

      if ( use_c14 )then
         do i = 1,ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  matrix_Cinter14%V(c,j+(i-1)*nlevdecomp)   = cs14_soil%decomp_cpools_vr_col(c,j,i)
                  C14inter_old(c,j+(i-1)*nlevdecomp)        = cs14_soil%decomp_cpools_vr_col(c,j,i)
               end do
            end do
         end do
      end if !c14
      call t_stopf('CN Soil matrix-assign matrix-inter')

      call Ksoiln%SetValueCopyDM(num_soilc,filter_soilc,Ksoil )
      do i = 1,ndecomp_pools
         if ( .not. floating_cn_ratio_decomp_pools(i) ) then
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if(ns_soil%decomp_npools_vr_col(c,j,i) > 0)then
                     Ksoiln%DM(c,j+(i-1)*nlevdecomp) = Ksoil%DM(c,j+(i-1)*nlevdecomp) * cs_soil%decomp_cpools_vr_col(c,j,i) / ns_soil%decomp_npools_vr_col(c,j,i) / initial_cn_ratio(i)
                  end if
               end do
            end do
         end if
      end do
      ! Save the C and N pool size at begin of each year, which are used to calculate C and N capacity at end of each year.
      call t_startf('CN Soil matrix-assign matrix-decomp0')
      if (is_beg_curr_year())then  
         iyr = iyr + 1
         if(mod(iyr-1,nyr_forcing) .eq. 0)then
            iloop = iloop + 1
         end if
         if(.not. spinup_matrixcn .or. spinup_matrixcn .and. mod(iyr-1,nyr_SASU) .eq. 0)then
            do i = 1,ndecomp_pools
               do j = 1, nlevdecomp
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     cs_soil%decomp0_cpools_vr_col(c,j,i)=cs_soil%decomp_cpools_vr_col(c,j,i)
                     if(use_c13)then
                        cs13_soil%decomp0_cpools_vr_col(c,j,i)=cs13_soil%decomp_cpools_vr_col(c,j,i)
                     end if
                     if(use_c14)then
                        cs14_soil%decomp0_cpools_vr_col(c,j,i)=cs14_soil%decomp_cpools_vr_col(c,j,i)
                     end if
                     ns_soil%decomp0_npools_vr_col(c,j,i)=ns_soil%decomp_npools_vr_col(c,j,i)
                  end do
               end do
            end do
            where(cs_soil%decomp0_cpools_vr_col .lt. epsi)
               cs_soil%decomp0_cpools_vr_col = epsi
            end where
            if(use_c13)then
               where(cs13_soil%decomp0_cpools_vr_col .lt. epsi*c3_r2)
                  cs13_soil%decomp0_cpools_vr_col = epsi*c3_r2
               end where
            end if
            if(use_c14)then
               where(cs14_soil%decomp0_cpools_vr_col .lt. epsi*c14ratio) 
                  cs14_soil%decomp0_cpools_vr_col = epsi*c14ratio
               end where
            end if
            where(ns_soil%decomp0_npools_vr_col .lt. epsi)
               ns_soil%decomp0_npools_vr_col = epsi
            end where
         end if
      end if
      call t_stopf('CN Soil matrix-assign matrix-decomp0')
   
      ! Set C transfer matrix Ac from a_ma_vr
      call t_startf('CN Soil matrix-matrix mult1-lev3-SetValueAK1')
      call AKsoilc%SetValueA(begc,endc,num_soilc,filter_soilc,a_ma_vr,A_i,A_j,Ntrans,init_readyAsoilc,list_Asoilc,RI_a,CI_a)
      call t_stopf('CN Soil matrix-matrix mult1-lev3-SetValueAK1')

      ! Set N transfer matrix An from na_ma_vr
      call t_startf('CN Soil matrix-matrix mult1-lev3-SetValueAK2')
      call AKsoiln%SetValueA(begc,endc,num_soilc,filter_soilc,na_ma_vr,A_i,A_j,Ntrans,init_readyAsoiln,list_Asoiln,RI_na,CI_na)
      call t_stopf('CN Soil matrix-matrix mult1-lev3-SetValueAK2')

      ! calculate matrix Ac*K for C
      call t_startf('CN Soil matrix-matrix mult1-lev3-SPMM_AK1')
      call AKsoilc%SPMM_AK(num_soilc,filter_soilc,Ksoil)
      call t_stopf('CN Soil matrix-matrix mult1-lev3-SPMM_AK1')

      ! calculate matrix An*K for N
      call t_startf('CN Soil matrix-matrix mult1-lev3-SPMM_AK2')
      call AKsoiln%SPMM_AK(num_soilc,filter_soilc,Ksoiln)
      call t_stopf('CN Soil matrix-matrix mult1-lev3-SPMM_AK2')

      ! Set vertical transfer matrix V from tri_ma_vr
      call t_startf('CN Soil matrix-matrix mult1-lev3-SetValueAV,AKfire')
      call AVsoil%SetValueSM(begc,endc,num_soilc,filter_soilc,tri_ma_vr(begc:endc,1:Ntri_setup),tri_i,tri_j,Ntri_setup)

      ! Set fire decomposition matrix Kfire from matrix_decomp_fire_k
      do j=1,ndecomp_pools_vr
         kfire_i(j) = j
         kfire_j(j) = j
      end do
      call AKfiresoil%SetValueSM(begc,endc,num_soilc,filter_soilc,matrix_decomp_fire_k(begc:endc,1:ndecomp_pools_vr),kfire_i,kfire_j,ndecomp_pools_vr)
      if(use_c14)then
         do i = 1,ndecomp_pools
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cf14_soil%matrix_decomp_fire_k_col(c,j+(i-1)*nlevdecomp) = cf14_soil%matrix_decomp_fire_k_col(c,j+(i-1)*nlevdecomp) &
                                                                           + matrix_decomp_fire_k(c,j+(i-1)*nlevdecomp)
               end do
            end do
         end do
         call cf14_soil%AKfiresoil%SetValueSM(begc,endc,num_soilc,filter_soilc,&
                             cf14_soil%matrix_decomp_fire_k_col(begc:endc,1:ndecomp_pools_vr),kfire_i,kfire_j,ndecomp_pools_vr)
      end if
      call t_stopf('CN Soil matrix-matrix mult1-lev3-SetValueAV,AKfire')

      ! Calculate AKallsoilc = A*K + AVsoil + AKfiresoil. (AKfiresoil = -Kfire)
      ! When no fire, AKallsoilc = A*K + AVsoil
      ! When fire is on, AKallsoilc = A*K + AVsoil + AKfiresoil
      ! Here, AKallsoilc represents all soil C transfer rate (gC/(gC*m3*step))
      ! and AKallsoiln represents all soil N transfer rate (gN/(gN*m3*step))

      call t_startf('CN Soil matrix-matrix mult1-lev3-SPMP_AB')
      if(num_actfirec .eq. 0)then
         call AKallsoilc%SPMP_AB(num_soilc,filter_soilc,AKsoilc,AVsoil,list_ready1_nofire,list_A=list_AK_AKV, list_B=list_V_AKV,&
              NE_AB=NE_AKallsoilc,RI_AB=RI_AKallsoilc,CI_AB=CI_AKallsoilc)
         call AKallsoiln%SPMP_AB(num_soilc,filter_soilc,AKsoiln,AVsoil,list_ready2_nofire,list_A=list_AK_AKV, list_B=list_V_AKV,&
              NE_AB=NE_AKallsoiln,RI_AB=RI_AKallsoiln,CI_AB=CI_AKallsoiln)
      else
         call AKallsoilc%SPMP_ABC(num_soilc,filter_soilc,AKsoilc,AVsoil,AKfiresoil,list_ready1_fire,list_A=list_AK_AKVfire,&
              list_B=list_V_AKVfire,list_C=list_fire_AKVfire,NE_ABC=NE_AKallsoilc,RI_ABC=RI_AKallsoilc,CI_ABC=CI_AKallsoilc,&
              use_actunit_list_C=.True.,num_actunit_C=num_actfirec,filter_actunit_C=filter_actfirec)
         call AKallsoiln%SPMP_ABC(num_soilc,filter_soilc,AKsoiln,AVsoil,AKfiresoil,list_ready2_fire,list_A=list_AK_AKVfire,&
              list_B=list_V_AKVfire,list_C=list_fire_AKVfire,NE_ABC=NE_AKallsoiln,RI_ABC=RI_AKallsoiln,CI_ABC=CI_AKallsoiln,&
              use_actunit_list_C=.True.,num_actunit_C=num_actfirec,filter_actunit_C=filter_actfirec)
      end if
      if(use_c14)then
         call cf14_soil%AKallsoilc%SPMP_ABC(num_soilc,filter_soilc,AKsoilc,AVsoil,cf14_soil%AKfiresoil,list_ready3_fire,&
              list_A=list_AK_AKVfire,list_B=list_V_AKVfire,list_C=list_fire_AKVfire,NE_ABC=cf14_soil%NE_AKallsoilc,&
              RI_ABC=cf14_soil%RI_AKallsoilc,CI_ABC=cf14_soil%CI_AKallsoilc)
      end if

      call t_stopf('CN Soil matrix-matrix mult1-lev3-SPMP_AB')

      call t_startf('CN Soil matrix-matrix mult2-lev2')

      ! Update soil C pool size: X(matrix_Cinter) = X(matrix_Cinter) + (A*K + AVsoil + AKfiresoil) * X(matrix_Cinter)
      ! Update soil N pool size: X(matrix_Ninter) = X(matrix_Ninter) + (A*K + AVsoil + AKfiresoil) * X(matrix_Ninter)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         do i=1,AVsoil%NE
            ilev = mod(AVsoil%RI(i)-1,nlevdecomp)+1
            idecomp = (AVsoil%RI(i) - ilev)/nlevdecomp + 1
            trcr_ctendency(c,ilev,idecomp) = trcr_ctendency(c,ilev,idecomp) + AVsoil%M(c,i)*matrix_Cinter%V(c,AVsoil%CI(i)) / dt
            trcr_ntendency(c,ilev,idecomp) = trcr_ntendency(c,ilev,idecomp) + AVsoil%M(c,i)*matrix_Ninter%V(c,AVsoil%CI(i)) / dt
         end do
      end do
     
      call matrix_Cinter%SPMM_AX(num_soilc,filter_soilc,AKallsoilc)
      call matrix_Ninter%SPMM_AX(num_soilc,filter_soilc,AKallsoiln)

      ! Update soil C13 pool size: X(matrix_Cinter13) = X(matrix_Cinter13) + (A*K + AVsoil + AKfiresoil) * X(matrix_Cinter13)
      if ( use_c13)then
         call matrix_Cinter13%SPMM_AX(num_soilc,filter_soilc,AKallsoilc)
      end if

      ! Update soil C14 pool size: X(matrix_Cinter14) = X(matrix_Cinter14) + (A*K + AVsoil + AKfiresoil) * X(matrix_Cinter14)
      if ( use_c14)then
         call matrix_Cinter14%SPMM_AX(num_soilc,filter_soilc,cf14_soil%AKallsoilc)
      end if

      ! Update soil C pool size: X(matrix_Cinter) = X(matrix_Cinter) + (A*K + AVsoil + AKfiresoil) * X(matrix_Cinter) + I(matrix_Cinput)
      ! Update soil N pool size: X(matrix_Ninter) = X(matrix_Ninter) + (A*K + AVsoil + AKfiresoil) * X(matrix_Ninter) + I(matrix_Ninput)
      do j = 1, ndecomp_pools_vr
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            matrix_Cinter%V(c,j) = matrix_Cinput%V(c,j) + matrix_Cinter%V(c,j)
            matrix_Ninter%V(c,j) = matrix_Ninput%V(c,j) + matrix_Ninter%V(c,j)
         end do
      end do

      ! Update soil C13 pool size: X(matrix_Cinter13) = X(matrix_Cinter13) + (A*K + AVsoil + AKfiresoil) * X(matrix_Cinter13) + I(matrix_Cinput13)
      if ( use_c13)then
         do j = 1, ndecomp_pools_vr
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               matrix_Cinter13%V(c,j) = matrix_Cinput13%V(c,j) + matrix_Cinter13%V(c,j)
            end do
         end do
      end if     
 
      ! Update soil C14 pool size: X(matrix_Cinter14) = X(matrix_Cinter14) + (A*K + AVsoil + AKfiresoil) * X(matrix_Cinter14) + I(matrix_Cinput14)
      if ( use_c14)then
         do j = 1, ndecomp_pools_vr
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               matrix_Cinter14%V(c,j) = matrix_Cinput14%V(c,j) + matrix_Cinter14%V(c,j)
            end do
         end do
      end if !c14

      ! Adjust heterotrophic respiration and fire flux because the pool size updating order is different between default and matrix code, balance error will occur
      ! while sudden big changes happen in C pool, eg. crop harvest and fire.
      call t_stopf('CN Soil matrix-matrix mult2-lev2')

      call t_startf('CN Soil matrix-assign back')

      ! Send vector type soil C and N pool size back to decomp_cpools_vr_col and decomp_npools_vr_col
      do i=1,ndecomp_pools
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cs_soil%decomp_cpools_vr_col(c,j,i) = matrix_Cinter%V(c,j+(i-1)*nlevdecomp)
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

      if(use_soil_matrixcn .and. (hist_wrt_matrixcn_diag .or. spinup_matrixcn))then
            
         ! Accumulate C transfers during a whole calendar year to calculate the C and N capacity
         do j=1,ndecomp_pools*nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cs_soil%in_acc (c,j) = cs_soil%in_acc (c,j) + matrix_Cinput%V(c,j)
               ns_soil%in_nacc(c,j) = ns_soil%in_nacc(c,j) + matrix_Ninput%V(c,j)
            end do
         end do
         if(use_c13)then
            do j=1,ndecomp_pools*nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs13_soil%in_acc (c,j) = cs13_soil%in_acc (c,j) + matrix_Cinput13%V(c,j)
               end do
            end do
         end if
         if(use_c14)then
            do j=1,ndecomp_pools*nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs14_soil%in_acc (c,j) = cs14_soil%in_acc (c,j) + matrix_Cinput14%V(c,j)
               end do
            end do
         end if

         if(use_c13)then
            call cf13_soil%AKallsoilc%SetValueCopySM (num_soilc,filter_soilc,AKallsoilc)
         end if

!         if(use_c14)then
!            call cf14_soil%AKallsoilc%SetValueCopySM (num_soilc,filter_soilc,AKallsoilc)
!         end if
         ! Calculate all the soil C transfers in current time step. 
         ! After this step, AKallsoilc represents all the C transfers (gC/m3/step)
         call Xdiagsoil%SetValueDM(begc,endc,num_soilc,filter_soilc,Cinter_old(begc:endc,1:ndecomp_pools_vr))
         call AKallsoilc%SPMM_AK(num_soilc,filter_soilc,Xdiagsoil)
        
         if(use_c13)then
            call Xdiagsoil%SetValueDM(begc,endc,num_soilc,filter_soilc,C13inter_old(begc:endc,1:ndecomp_pools_vr))
            call cf13_soil%AKallsoilc%SPMM_AK(num_soilc,filter_soilc,Xdiagsoil)
         end if
        
         if(use_c14)then
            call Xdiagsoil%SetValueDM(begc,endc,num_soilc,filter_soilc,C14inter_old(begc:endc,1:ndecomp_pools_vr))
            call cf14_soil%AKallsoilc%SPMM_AK(num_soilc,filter_soilc,Xdiagsoil)
         end if
        
         ! Calculate all the soil N transfers in current time step
         ! After this step, AKallsoiln represents all the N transfers (gN/m3/step)
         call Xdiagsoil%SetValueDM(begc,endc,num_soilc,filter_soilc,Ninter_old(begc:endc,1:ndecomp_pools_vr))
         call AKallsoiln%SPMM_AK(num_soilc,filter_soilc,Xdiagsoil)

         !
         ! Accumulate soil C transfers: AKXcacc = AKXcacc + AKallsoilc
         !
         ! Copy indices from AKallsoilc on restart step
         if ( is_first_restart_step() )then
           call AKXcacc%CopyIdxSM( AKallsoilc )
         end if
         if ( AKXcacc%IsValuesSetSM() )then
            call AKXcacc%SPMP_B_ACC(num_soilc,filter_soilc,AKallsoilc)
         else
            ! This should only happen on the first time-step
            call AKXcacc%SetValueCopySM(num_soilc,filter_soilc,AKallsoilc)
         end if

         !
         ! Accumulate soil C13 transfers: cs13_soil%AKXcacc = cs13_soil%AKXcacc + cs13_soil%AKallsoilc
         !
         ! Copy indices from AKallsoilc on restart step
         if(use_c13)then
            if ( is_first_restart_step() )then
               call cs13_soil%AKXcacc%CopyIdxSM( cf13_soil%AKallsoilc )
            end if
            if ( cs13_soil%AKXcacc%IsValuesSetSM() )then
               call cs13_soil%AKXcacc%SPMP_B_ACC(num_soilc,filter_soilc,cf13_soil%AKallsoilc)
            else
            ! This should only happen on the first time-step
               call cs13_soil%AKXcacc%SetValueCopySM(num_soilc,filter_soilc,cf13_soil%AKallsoilc)
            end if
         end if

         !
         ! Accumulate soil C14 transfers: cs14_soil%AKXcacc = cs14_soil%AKXcacc + cs14_soil%AKallsoilc
         !
         ! Copy indices from AKallsoilc on restart step
         if(use_c14)then
            if ( is_first_restart_step() )then
               call cs14_soil%AKXcacc%CopyIdxSM( cf14_soil%AKallsoilc )
            end if
            if ( cs14_soil%AKXcacc%IsValuesSetSM() )then
               call cs14_soil%AKXcacc%SPMP_B_ACC(num_soilc,filter_soilc,cf14_soil%AKallsoilc)
            else
            ! This should only happen on the first time-step
               call cs14_soil%AKXcacc%SetValueCopySM(num_soilc,filter_soilc,cf14_soil%AKallsoilc)
            end if
         end if

         !
         ! Accumulate soil N transfers: AKXnacc = AKXnacc + AKallsoiln
         !
         ! Copy indices from AKallsoiln on restart step
         if ( is_first_restart_step() )then
           call AKXnacc%CopyIdxSM( AKallsoiln )
         end if
         if ( AKXnacc%IsValuesSetSM() )then
            call AKXnacc%SPMP_B_ACC(num_soilc,filter_soilc,AKallsoiln)
         else
            ! This should only happen on the first time-step
            call AKXnacc%SetValueCopySM(num_soilc,filter_soilc,AKallsoiln)
         end if

         call t_startf('CN Soil matrix-calc. C capacity')
         if((.not. spinup_matrixcn .and. is_end_curr_year()) .or. (spinup_matrixcn .and. is_end_curr_year() .and. mod(iyr,nyr_SASU) .eq. 0))then
            ! Copy C transfers from sparse matrix to 2D temporary variables tran_acc and tran_nacc
            ! Calculate the C and N transfer rate by dividing CN transfer by base value saved at begin of each year.
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cs_soil%tran_acc (c,1:ndecomp_pools_vr,1:ndecomp_pools_vr) = 0._r8
               ns_soil%tran_nacc(c,1:ndecomp_pools_vr,1:ndecomp_pools_vr) = 0._r8
               if(use_c13)then
                  cs13_soil%tran_acc (c,1:ndecomp_pools_vr,1:ndecomp_pools_vr) = 0._r8
               end if
               if(use_c14)then
                  cs14_soil%tran_acc (c,1:ndecomp_pools_vr,1:ndecomp_pools_vr) = 0._r8
               end if
            end do
            do j=1,n_all_entries
               j_lev    = mod(all_j(j)-1,nlevdecomp)+1
               j_decomp = (all_j(j) - j_lev)/nlevdecomp + 1
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cs_soil%tran_acc(c,all_i(j),all_j(j))  = AKXcacc%M(c,j) / cs_soil%decomp0_cpools_vr_col(c,j_lev,j_decomp)
                  ns_soil%tran_nacc(c,all_i(j),all_j(j)) = AKXnacc%M(c,j) / ns_soil%decomp0_npools_vr_col(c,j_lev,j_decomp)
                  if(use_c13)then
                     cs13_soil%tran_acc(c,all_i(j),all_j(j))  = cs13_soil%AKXcacc%M(c,j) / cs13_soil%decomp0_cpools_vr_col(c,j_lev,j_decomp)
                  end if
                  if(use_c14)then
                     cs14_soil%tran_acc(c,all_i(j),all_j(j))  = cs14_soil%AKXcacc%M(c,j) / cs14_soil%decomp0_cpools_vr_col(c,j_lev,j_decomp)
                  end if
               end do
            end do

            do i=1,ndecomp_pools_vr
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (abs(cs_soil%tran_acc(c,i,i)) .le. epsi)then !avoid inversion nan
                      cs_soil%tran_acc(c,i,i) = 1.e+36_r8
                  end if 
               end do
            end do
        
            if(use_c13)then
               do i=1,ndecomp_pools_vr
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     if (abs(cs13_soil%tran_acc(c,i,i)) .le. epsi)then !avoid inversion nan
                         cs13_soil%tran_acc(c,i,i) = 1.e+36_r8
                     end if 
                  end do
               end do
            end if

            if(use_c14)then
               do i=1,ndecomp_pools_vr
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     if (abs(cs14_soil%tran_acc(c,i,i)) .le. epsi)then !avoid inversion nan
                         cs14_soil%tran_acc(c,i,i) = 1.e+36_r8
                     end if 
                  end do
               end do
            end if

            do i=1,ndecomp_pools_vr
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (abs(ns_soil%tran_nacc(c,i,i)) .le. epsi)then
                     ns_soil%tran_nacc(c,i,i) = 1.e+36_r8
                  end if 
               end do
            end do

            ! Calculate capacity 
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               call inverse(cs_soil%tran_acc(c,1:ndecomp_pools_vr,1:ndecomp_pools_vr),AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
               soilmatrixc_cap(c,:,1) = -matmul(AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),cs_soil%in_acc(c,1:ndecomp_pools_vr))
               if(use_c13)then
                  call inverse(cs13_soil%tran_acc(c,1:ndecomp_pools_vr,1:ndecomp_pools_vr),AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
                  soilmatrixc13_cap(c,:,1) = -matmul(AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),cs13_soil%in_acc(c,1:ndecomp_pools_vr))
               end if
               if(use_c14)then
                  call inverse(cs14_soil%tran_acc(c,1:ndecomp_pools_vr,1:ndecomp_pools_vr),AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
                  soilmatrixc14_cap(c,:,1) = -matmul(AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),cs14_soil%in_acc(c,1:ndecomp_pools_vr))
               end if
               call inverse(ns_soil%tran_nacc(c,1:ndecomp_pools_vr,1:ndecomp_pools_vr),AKinvn(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
               soilmatrixn_cap(c,:,1) = -matmul(AKinvn(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ns_soil%in_nacc(c,1:ndecomp_pools_vr))
            end do

            do i = 1,ndecomp_pools
               do j = 1,nlevdecomp
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     if(soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1) .lt. 0)then
                        soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1) = 0._r8
                     endif
                     if(use_c13 .and. soilmatrixc13_cap(c,j+(i-1)*nlevdecomp,1) .lt. 0)then
                        soilmatrixc13_cap(c,j+(i-1)*nlevdecomp,1) = 0._r8
                     endif
                     if(use_c14 .and. soilmatrixc14_cap(c,j+(i-1)*nlevdecomp,1) .lt. 0)then
                        soilmatrixc14_cap(c,j+(i-1)*nlevdecomp,1) = 0._r8
                     endif
                     if(soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1) .lt. 0)then
                        soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1) = 0._r8
                     endif
                  end do
               end do
            end do
                        

            do i = 1,ndecomp_pools
               do j = 1,nlevdecomp
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     if(nyr_SASU .eq. nyr_forcing .and. &
                                     (soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1)/cs_soil%decomp0_cpools_vr_col(c,j,i) .gt. 100 .and. soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1) .gt. 1.e+5_r8  &
                                 .or. soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1)/ns_soil%decomp0_npools_vr_col(c,j,i) .gt. 100 .and. soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1) .gt. 1.e+3_r8) &
                   .or. nyr_SASU .lt. nyr_forcing .and. i .eq. i_cwd .and. &
                                     (soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1)/cs_soil%decomp0_cpools_vr_col(c,j,i) .gt. 100 .and. soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1) .gt. 1.e+5_r8  &
                                 .or. soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1)/ns_soil%decomp0_npools_vr_col(c,j,i) .gt. 100 .and. soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1) .gt. 1.e+3_r8) )then
                        soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1) = matrix_Cinter%V(c,j+(i-1)*nlevdecomp)
                        if(use_c13)then
                           soilmatrixc13_cap(c,j+(i-1)*nlevdecomp,1) = matrix_Cinter13%V(c,j+(i-1)*nlevdecomp)
                        end if
                        if(use_c14)then
                           soilmatrixc14_cap(c,j+(i-1)*nlevdecomp,1) = matrix_Cinter14%V(c,j+(i-1)*nlevdecomp)
                        end if
                        soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1) = matrix_Ninter%V(c,j+(i-1)*nlevdecomp)
                     end if
                  end do
               end do
            end do

            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if(any(soilmatrixc_cap(c,:,1) .gt. 1.e+8_r8) .or. any(soilmatrixn_cap(c,:,1) .gt. 1.e+8_r8))then
                  soilmatrixc_cap(c,:,1) = matrix_Cinter%V(c,:)
                  if(use_c13)then
                     soilmatrixc13_cap(c,:,1) = matrix_Cinter13%V(c,:)
                  end if
                  if(use_c14)then
                     soilmatrixc14_cap(c,:,1) = matrix_Cinter14%V(c,:)
                  end if
                  soilmatrixn_cap(c,:,1) = matrix_Ninter%V(c,:)
               end if
            end do

            ! If spin up is on, the capacity replaces the pool size with capacity.
            ! Copy the capacity into a 3D variable, and be ready to write to history files.
            do i=1,ndecomp_pools 
               do j = 1,nlevdecomp
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     if(spinup_matrixcn .and. .not. is_first_step_of_this_run_segment())then
                        cs_soil%decomp_cpools_vr_col(c,j,i) =  soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1)
                        if(use_c13)then
                           cs13_soil%decomp_cpools_vr_col(c,j,i) =  soilmatrixc13_cap(c,j+(i-1)*nlevdecomp,1)
                        end if
                        if(use_c14)then
                           cs14_soil%decomp_cpools_vr_col(c,j,i) =  soilmatrixc14_cap(c,j+(i-1)*nlevdecomp,1)
                        end if
                        if(floating_cn_ratio_decomp_pools(i))then
                           ns_soil%decomp_npools_vr_col(c,j,i) =  soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1)
                        else
                           ns_soil%decomp_npools_vr_col(c,j,i) = cs_soil%decomp_cpools_vr_col(c,j,i) / cn_decomp_pools(c,j,i)
                        end if
                        ! calculate the average of the storage capacity when iloop equals to iloop_avg
                        if(iloop .eq. iloop_avg)then
                           cs_soil%decomp_cpools_vr_SASUsave_col(c,j,i) = cs_soil%decomp_cpools_vr_SASUsave_col(c,j,i) + cs_soil%decomp_cpools_vr_col(c,j,i)
                           if(use_c13)then
                              cs13_soil%decomp_cpools_vr_SASUsave_col(c,j,i) = cs13_soil%decomp_cpools_vr_SASUsave_col(c,j,i) + cs13_soil%decomp_cpools_vr_col(c,j,i)
                           end if
                           if(use_c14)then
                              cs14_soil%decomp_cpools_vr_SASUsave_col(c,j,i) = cs14_soil%decomp_cpools_vr_SASUsave_col(c,j,i) + cs14_soil%decomp_cpools_vr_col(c,j,i)
                           end if
                           ns_soil%decomp_npools_vr_SASUsave_col(c,j,i) = ns_soil%decomp_npools_vr_SASUsave_col(c,j,i) + ns_soil%decomp_npools_vr_col(c,j,i)
                           if(iyr .eq. nyr_forcing)then
                              cs_soil%decomp_cpools_vr_col(c,j,i) = cs_soil%decomp_cpools_vr_SASUsave_col(c,j,i) / (nyr_forcing/nyr_SASU)
                              if(use_c13)then
                                 cs13_soil%decomp_cpools_vr_col(c,j,i) = cs13_soil%decomp_cpools_vr_SASUsave_col(c,j,i) / (nyr_forcing/nyr_SASU)
                              end if
                              if(use_c14)then
                                 cs14_soil%decomp_cpools_vr_col(c,j,i) = cs14_soil%decomp_cpools_vr_SASUsave_col(c,j,i) / (nyr_forcing/nyr_SASU)
                              end if
                              ns_soil%decomp_npools_vr_col(c,j,i) = ns_soil%decomp_npools_vr_SASUsave_col(c,j,i) / (nyr_forcing/nyr_SASU)
                              cs_soil%decomp_cpools_vr_SASUsave_col(c,j,i) = 0._r8
                              if(use_c13)then
                                 cs13_soil%decomp_cpools_vr_SASUsave_col(c,j,i) = 0._r8
                              end if
                              if(use_c14)then
                                 cs14_soil%decomp_cpools_vr_SASUsave_col(c,j,i) = 0._r8
                              end if
                              ns_soil%decomp_npools_vr_SASUsave_col(c,j,i) = 0._r8
                           end if
                        end if
                     end if
                     cs_soil%matrix_cap_decomp_cpools_vr_col(c,j,i) = soilmatrixc_cap(c,j+(i-1)*nlevdecomp,1)
                     if(use_c13)then
                        cs13_soil%matrix_cap_decomp_cpools_vr_col(c,j,i) = soilmatrixc13_cap(c,j+(i-1)*nlevdecomp,1)
                     end if
                     if(use_c14)then
                        cs14_soil%matrix_cap_decomp_cpools_vr_col(c,j,i) = soilmatrixc14_cap(c,j+(i-1)*nlevdecomp,1)
                     end if
                     ns_soil%matrix_cap_decomp_npools_vr_col(c,j,i) = soilmatrixn_cap(c,j+(i-1)*nlevdecomp,1)
                  end do
               end do
            end do
         
            if(spinup_matrixcn)call update_DA_nstep()
            if(iloop .eq. iloop_avg .and. iyr .eq. nyr_forcing)iloop = 0
            if(iyr .eq. nyr_forcing)iyr = 0

            ! Reset to accumulation variables to 0 at end of each year
            do j=1,n_all_entries
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  AKXcacc%M(c,j) = 0._r8
                  if(use_c13)then
                     cs13_soil%AKXcacc%M(c,j) = 0._r8
                  end if
                  if(use_c14)then
                     cs14_soil%AKXcacc%M(c,j) = 0._r8
                  end if
                  AKXnacc%M(c,j) = 0._r8
               end do
            end do

            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cs_soil%in_acc   (c,:)   = 0._r8
               if(use_c13)then
                  cs13_soil%in_acc   (c,:)   = 0._r8
               end if
               if(use_c14)then
                  cs14_soil%in_acc   (c,:)   = 0._r8
               end if
               ns_soil%in_nacc  (c,:)   = 0._r8
            end do
         end if
         call t_stopf('CN Soil matrix-calc. C capacity')
      end if !is out_matrix

   end associate 
 end subroutine CNSoilMatrix

  !-----------------------------------------------------------------------
   subroutine CNSoilMatrixRest( ncid, flag )
    ! !DESCRIPTION:
    !
    !    Read/write restart data needed for the CN soil Matrix model solution
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
    call restartvar(ncid=ncid, flag=flag, varname='soil_cycle_year', xtype=ncd_int,  &
            long_name='Year number in soil spinup cycle sequence', units='years', &
            interpinic_flag='skip', readvar=readvar, data=iyr)

    call restartvar(ncid=ncid, flag=flag, varname='soil_cycle_loop', xtype=ncd_int,  &
            long_name='Loop number in soil spinup cycle sequence', units='years', &
            interpinic_flag='skip', readvar=readvar, data=iloop)

    !------------------------------------------------------------------------
   end subroutine CNSoilMatrixRest
 
end module CNSoilMatrixMod

