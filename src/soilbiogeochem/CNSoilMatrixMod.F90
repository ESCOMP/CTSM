module CNSoilMatrixMod

!#include "shr_assert.h"
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for CLM45BGC/CLM5.0BGC matrices
  
  ! The matrix equation 
  ! Xn+1 = Xn + I*dt + (A*ksi*k - tri/dz)*Xn*dt
  
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use shr_log_mod                    , only : errMsg => shr_log_errMsg
  use decompMod                      , only : bounds_type  
  use abortutils                     , only : endrun
  use clm_time_manager               , only : get_step_size, is_end_curr_month,get_curr_date,get_days_per_year
  use clm_time_manager               , only : is_first_step_of_this_run_segment,is_beg_curr_year
  use clm_varpar                     , only : ndecomp_pools, nlevdecomp, ndecomp_pools_vr        !number of biogeochemically active soil layers
  use clm_varpar                     , only :ndecomp_cascade_transitions
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
  use clm_varctl                     , only : use_c13, use_c14
!
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNSoilMatrix
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNSoilMatrix(bounds,num_soilc, filter_soilc, &
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
    type(cnveg_carbonflux_type)              , intent(in)    :: cnveg_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)    , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)     , intent(in)    :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_state_type)          , intent(inout) :: soilbiogeochem_state_inst
    type(cnveg_nitrogenflux_type)            , intent(in)    :: cnveg_nitrogenflux_inst
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
    real(r8):: epsi                 !small number
    logical  ::  end_of_year
    real(r8):: days_per_year,decay_const,half_life
    integer, save :: counter=0

    integer :: begc,endc                                    ! bounds 
    real(r8),dimension(ndecomp_pools_vr,ndecomp_pools_vr) :: a_ma_vr, kk_ma_vr, kk_fire_vr,tri_ma_vr, tranvert,ntranvert
    real(r8),dimension(ndecomp_pools_vr,ndecomp_pools_vr) :: na_ma_vr
!    real(r8),dimension(ndecomp_pools,ndecomp_pools) :: a_ma, kk_ma, matrix_soil_cn
!    real(r8),dimension(ndecomp_pools,ndecomp_pools) :: na_ma
    real(r8),dimension(num_soilc,nlevdecomp) :: depth_scalar,two_scalar
    real(r8),dimension(num_soilc,nlevdecomp) :: n_scalar_ave,t_scalar_ave,w_scalar_ave, o_scalar_ave
    real(r8)::  a_tri_ave_c(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8)::  b_tri_ave_c(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8)::  c_tri_ave_c(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8)::  input_ave_c(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8)::  a_tri_ave_n(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8)::  b_tri_ave_n(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8)::  c_tri_ave_n(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8)::  input_ave_n(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)


    real(r8),dimension(ndecomp_pools_vr,1) :: matrix_Cinter,matrix_Cinter_next,matrix_Cinput_vector,emulator_tmp,emulator_tmp1,emulator_tmpn

    real(r8),dimension(ndecomp_pools_vr,1) :: matrix13_Cinter,   matrix14_Cinter
    real(r8),dimension(ndecomp_pools_vr,1) :: matrix_Cinput13_vector,matrix_Cinput14_vector
    real(r8),dimension(ndecomp_pools_vr,1) :: matrix_Cinter13_next,matrix_Cinter14_next

    real(r8),dimension(ndecomp_pools_vr,1) ::  soilmatrixc_cap,soilmatrixn_cap
    real(r8),dimension(ndecomp_pools_vr,ndecomp_pools_vr) ::    matrix_Cinter_2d,matrix_Ninter_2d
    real(r8), dimension(ndecomp_pools_vr,ndecomp_pools_vr)   ::  AKinv,AKinvn
    real(r8),dimension(ndecomp_pools,1) :: emulator_Cinter_1,emulator_Cinter_next_1,emulator_Cinput_1

    real(r8),dimension(ndecomp_pools_vr,1) :: matrix_Ninter,matrix_Ninter_next,matrix_Ninput_vector
    real(r8),dimension(ndecomp_pools_vr,1) :: cmatrix_in
    real(r8),dimension(ndecomp_pools,1) :: matrix_Ninter_1,matrix_Ninter_next_1,matrix_Ninput_1

    real(r8),dimension(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools) :: cn_decomp_pools
!input for 3-step update
    real(r8), dimension(num_soilc,nlevdecomp) :: tot_c_to_litr_met_c,tot_c_to_litr_cel_c,tot_c_to_litr_lig_c,tot_c_to_cwdc 
    real(r8), dimension(num_soilc,nlevdecomp) :: tot_n_to_litr_met_n,tot_n_to_litr_cel_n,tot_n_to_litr_lig_n,tot_n_to_cwdn
    
    real(r8), dimension(num_soilc,nlevdecomp) :: tot_c_to_litr_met_c1,tot_c_to_litr_cel_c1,tot_c_to_litr_lig_c1,tot_c_to_cwdc1 
    real(r8), dimension(num_soilc,nlevdecomp) :: tot_n_to_litr_met_n1,tot_n_to_litr_cel_n1,tot_n_to_litr_lig_n1,tot_n_to_cwdn1
  
    real(r8), dimension(num_soilc,nlevdecomp) :: tot_c_to_litr_met_c23,tot_c_to_litr_cel_c23,tot_c_to_litr_lig_c23,tot_c_to_cwdc23 
    real(r8), dimension(num_soilc,nlevdecomp) :: tot_n_to_litr_met_n23,tot_n_to_litr_cel_n23,tot_n_to_litr_lig_n23,tot_n_to_cwdn23
    integer :: i_soil1
    integer :: i_soil2
    integer :: i_soil3
!    integer :: i_soil4	

    real(r8):: k_l1                         ! decomposition rate constant litter 1 (1/sec)
    real(r8):: k_l2_l3                      ! decomposition rate constant litter 2 and litter 3 (1/sec)
    real(r8):: k_s1                         ! decomposition rate constant SOM 1 (1/sec)
    real(r8):: k_s2                         ! decomposition rate constant SOM 2 (1/sec)
    real(r8):: k_s3                         ! decomposition rate constant SOM 3 (1/sec)
    real(r8):: k_frag                       ! fragmentation rate constant CWD (1/sec)
    real(r8):: tau_l1                       ! turnover time of  litter 1 (yr)
    real(r8):: tau_l2_l3                    ! turnover time of  litter 2 and litter 3 (yr)
    real(r8):: tau_l3                       ! turnover time of  litter 3 (yr)
    real(r8):: tau_s1                       ! turnover time of  SOM 1 (yr)
    real(r8):: tau_s2                       ! turnover time of  SOM 2 (yr)
    real(r8):: tau_s3                       ! turnover time of  SOM 3 (yr)
    real(r8):: tau_cwd                      ! corrected fragmentation rate constant CWD
!    real(r8):: days_per_year 

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

         matrix_decomp_k  => soilbiogeochem_carbonflux_inst%matrix_decomp_k_col, & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         fpi_vr          => soilbiogeochem_state_inst%fpi_vr_col    , & ! Input:  [real(r8) (:,:)   ]  fraction of potential immobilization (no units) 
         cascade_donor_pool => decomp_cascade_con%cascade_donor_pool, &
         cascade_receiver_pool          => decomp_cascade_con%cascade_receiver_pool,& 
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools, &
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio, &
         rf_decomp_cascade => soilbiogeochem_state_inst%rf_decomp_cascade_col,&
         pathfrac_decomp_cascade => soilbiogeochem_state_inst%pathfrac_decomp_cascade_col, &
         is_cwd                         => decomp_cascade_con%is_cwd                                         , & ! Input:  [logical  (:)     ]  TRUE => pool is a litter pool
         is_litter                      => decomp_cascade_con%is_litter                                      , & ! Input:  [logical  (:)     ]  TRUE => pool is a litter pool
         matrix_a_tri     => soilbiogeochem_carbonflux_inst%matrix_a_tri_col, &!(begc:endc,1:nlevdecomp),   & ! Output: "A"-for matrix
         matrix_b_tri     => soilbiogeochem_carbonflux_inst%matrix_b_tri_col, &!(begc:endc,1:nlevdecomp),   & ! Output: "B"-for matrix
         matrix_c_tri     => soilbiogeochem_carbonflux_inst%matrix_c_tri_col, &!(begc:endc,1:nlevdecomp),   & ! Output: "C"-for matrix 
         matrix_Cinput     => soilbiogeochem_carbonflux_inst%matrix_input_col, &!(begc:endc,1:nlevdecomp,1:ndecomp_pools),   & ! Output: "C"-for matrix 
         matrix_Ninput     => soilbiogeochem_nitrogenflux_inst%matrix_input_col, &!(begc:endc,1:nlevdecomp,1:ndecomp_pools),  & ! Output: "C"-for matrix 
         matrix_decomp_fire_k => soilbiogeochem_carbonflux_inst%matrix_decomp_fire_k_col &!(begc:endc,1:nlevdecomp,1:ndecomp_pools) &
         )

     ! set time steps
      dt = real( get_step_size(), r8 )
      days_per_year = get_days_per_year()
      secspyear = days_per_year* secspday
      half_life = 5730._r8 * secspyear
      decay_const = - log(0.5_r8) / half_life

      epsi = 1.e-30_r8 
      counter = counter + dt
     if (counter >= 1*secspyear) then ! link to the recycling span of climate forcing
          end_of_year = .true.
          counter = 0._r8
       else
          end_of_year = .false.
       end if

      !! turnover rate and time. Copied from SoilBiogeochemDecompCascadeBGCMod
      ! the belowground parameters from century
      tau_l1 = 1./18.5
      tau_l2_l3 = 1./4.9
      tau_s1 = 1./7.3
      tau_s2 = 1./0.2
      tau_s3 = 1./.0045

      ! century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1
      tau_cwd  = 1./0.3

      ! translate to per-second time constant
      k_l1 = 1._r8    / (secspday * days_per_year * tau_l1)
      k_l2_l3 = 1._r8 / (secspday * days_per_year * tau_l2_l3)
      k_s1 = 1._r8    / (secspday * days_per_year * tau_s1)
      k_s2 = 1._r8    / (secspday * days_per_year * tau_s2)
      k_s3 = 1._r8    / (secspday * days_per_year * tau_s3)
      k_frag = 1._r8  / (secspday * days_per_year * tau_cwd)
    
  
         i_soil1 = 5
         i_soil2 = 6
         i_soil3 = 7  
      ! calculate c:n ratios of applicable pools
      do l = 1, ndecomp_pools
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if ( floating_cn_ratio_decomp_pools(l) .and. cs_soil%decomp_cpools_vr_col(c,j,l) > 0._r8 ) then
                     if ( ns_soil%decomp_npools_vr_col(c,j,l) > 0._r8 ) then
                         cn_decomp_pools(c,j,l) = cs_soil%decomp_cpools_vr_col(c,j,l) / ns_soil%decomp_npools_vr_col(c,j,l)
                     end if
                  else
                     cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
                  end if 
               end do
            end do
      end do
        a_ma_vr = 0.0_r8
        na_ma_vr = 0.0_r8
        kk_ma_vr = 0.0_r8   ! kk_matrix, decay matrix * scalar matrix
        kk_fire_vr = 0.0_r8   ! kk_matrix, decay matrix * scalar matrix
        tri_ma_vr = 0.0_r8
        matrix_Cinput_vector = 0.0_r8
        matrix_Ninput_vector = 0.0_r8

        matrix_Cinter_2d = 0.0_r8
        matrix_Ninter_2d = 0.0_r8
        AKinv(:,:)=0._r8    
        AKinvn(:,:)=0._r8    

    do fc = 1,num_soilc
       c = filter_soilc(fc)
!       if (use_vertsoilc) then
          do j=1,ndecomp_pools_vr   !70
             a_ma_vr(j,j) = -1.0_r8
             na_ma_vr(j,j)= -1.0_r8
          enddo

!             cmatrix_in    = 0.0_r8

          do j = 1, nlevdecomp
             do i = 1,ndecomp_pools
                matrix_Cinput_vector(j+(i-1)*nlevdecomp,1) = max(matrix_Cinput(c,j,i),epsi) 
                matrix_Ninput_vector(j+(i-1)*nlevdecomp,1) = max(matrix_Ninput(c,j,i),epsi)
                if ( use_c13 ) then
                   matrix_Cinput13_vector(j+(i-1)*nlevdecomp,1) = max(cf13_soil%matrix_input_col(c,j,i),epsi)
                end if !c13
                if ( use_c14 ) then
                   matrix_Cinput14_vector(j+(i-1)*nlevdecomp,1) = max(cf14_soil%matrix_input_col(c,j,i),epsi)
                end if !c14
                kk_fire_vr((i-1)*nlevdecomp+j,(i-1)*nlevdecomp+j) = matrix_decomp_fire_k(c,j,i) 
                if(is_litter(i))then
                   kk_ma_vr((i-1)*nlevdecomp+j,(i-1)*nlevdecomp+j)   = matrix_decomp_k(c,j,i) * fpi_vr(c,j) 
                else
                   kk_ma_vr((i-1)*nlevdecomp+j,(i-1)*nlevdecomp+j)   = matrix_decomp_k(c,j,i)
                end if
                if(.not. is_cwd(i))then
                   tri_ma_vr(j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp) = matrix_b_tri(c,j) / dzsoi_decomp(j)
                   if (j < nlevdecomp) then    ! avoid tranfer from for example, litr3_10th layer to soil1_1st layer
                      tri_ma_vr(j+(i-1)*nlevdecomp,j+1+(i-1)*nlevdecomp) = matrix_c_tri(c,j) / dzsoi_decomp(j) 
                   end if
                   if (j > 1) then  ! avoid tranfer from for example,soil1_1st layer to litr3_10th layer
                      tri_ma_vr(j+(i-1)*nlevdecomp,j-1+(i-1)*nlevdecomp)=  matrix_a_tri(c,j) / dzsoi_decomp(j) 
                   endif
                end if
                matrix_Cinter(j+(i-1)*nlevdecomp,1)   = cs_soil%decomp_cpools_vr_col(c,j,i)                
                matrix_Ninter(j+(i-1)*nlevdecomp,1)   = ns_soil%decomp_npools_vr_col(c,j,i)     
               if ( use_c13 )then
                 matrix13_Cinter(j+(i-1)*nlevdecomp,1)   = cs13_soil%decomp_cpools_vr_col(c,j,i)
               end if !c13

               if ( use_c14 )then
                 matrix14_Cinter(j+(i-1)*nlevdecomp,1)   = cs14_soil%decomp_cpools_vr_col(c,j,i)
               end if !c14
       
              if (is_beg_curr_year() .or.is_first_step_of_this_run_segment() )then  
                 cs_soil%decomp0_cpools_vr_col(c,j,i)=max(cs_soil%decomp_cpools_vr_col(c,j,i),epsi)
                 ns_soil%decomp0_npools_vr_col(c,j,i)=max(ns_soil%decomp_npools_vr_col(c,j,i),epsi)
              end if
                 matrix_Cinter_2d(j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp)=cs_soil%decomp_cpools_vr_col(c,j,i)/cs_soil%decomp0_cpools_vr_col(c,j,i)
                 matrix_Ninter_2d(j+(i-1)*nlevdecomp,j+(i-1)*nlevdecomp)=ns_soil%decomp_npools_vr_col(c,j,i)/ns_soil%decomp0_npools_vr_col(c,j,i)
             end do

             do k = 1, ndecomp_cascade_transitions
                a_ma_vr((cascade_receiver_pool(k)-1)*nlevdecomp+j,(cascade_donor_pool(k)-1)*nlevdecomp+j) = (1.0-rf_decomp_cascade(c,j,k))*pathfrac_decomp_cascade(c,j,k)
                if( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)))then
                   na_ma_vr((cascade_receiver_pool(k)-1)*nlevdecomp+j,(cascade_donor_pool(k)-1)*nlevdecomp+j) = (1.0-rf_decomp_cascade(c,j,k))* &
                            (cn_decomp_pools(c,j,cascade_donor_pool(k))/cn_decomp_pools(c,j,cascade_receiver_pool(k)))*pathfrac_decomp_cascade(c,j,k)
                else
                   na_ma_vr((cascade_receiver_pool(k)-1)*nlevdecomp+j,(cascade_donor_pool(k)-1)*nlevdecomp+j) = pathfrac_decomp_cascade(c,j,k)
                end if
             end do
          end do
          tranvert = matmul(a_ma_vr,kk_ma_vr)-tri_ma_vr-kk_fire_vr  !intermediate calculatio
          ntranvert = matmul(na_ma_vr,kk_ma_vr)-tri_ma_vr-kk_fire_vr  !intermediate calculatio
 
          matrix_Cinter_next(:,:) = matrix_Cinter + matrix_Cinput_vector + &
                                   matmul(tranvert, matrix_Cinter)*dt

          matrix_Ninter_next(:,:) = matrix_Ninter + matrix_Ninput_vector + &
                                   matmul(ntranvert,matrix_Ninter)*dt
         if ( use_c13)then
              matrix_Cinter13_next(:,:) = matrix13_Cinter + matrix_Cinput13_vector + &
                                   matmul(tranvert, matrix13_Cinter)*dt
          do j = 1,nlevdecomp
             do i=1,ndecomp_pools
                cs13_soil%decomp_cpools_vr_col(c,j,i) = matrix_Cinter13_next(j+(i-1)*nlevdecomp,1)
             end do
          end do
         end if !c13
         if ( use_c14)then
              matrix_Cinter14_next(:,:) = matrix14_Cinter + matrix_Cinput14_vector + &
                                   matmul(tranvert, matrix14_Cinter)*dt
          do j = 1,nlevdecomp
             do i=1,ndecomp_pools
                cs14_soil%decomp_cpools_vr_col(c,j,i) = matrix_Cinter14_next(j+(i-1)*nlevdecomp,1)
             end do
          end do
         end if !c14


          do j = 1,nlevdecomp
             do i=1,ndecomp_pools
                cs_soil%decomp_cpools_vr_col(c,j,i) = matrix_Cinter_next(j+(i-1)*nlevdecomp,1)
                ns_soil%decomp_npools_vr_col(c,j,i) = matrix_Ninter_next(j+(i-1)*nlevdecomp,1)
             end do
          end do

         if(use_soil_matrixcn .and. (is_outmatrix .or. isspinup))then
            cs_soil%in_acc(c,:) = cs_soil%in_acc(c,:) + matrix_Cinput_vector(:,1) 
            cs_soil%tran_acc(c,1:ndecomp_pools_vr,:) = cs_soil%tran_acc(c,1:ndecomp_pools_vr,:)&
                                         + matmul(tranvert(1:ndecomp_pools_vr,:),matrix_Cinter_2d(:,:))*dt 
         
            ns_soil%in_nacc(c,:) = ns_soil%in_nacc(c,:) + matrix_Ninput_vector(:,1) 
            ns_soil%tran_nacc(c,1:ndecomp_pools_vr,:) = ns_soil%tran_nacc(c,1:ndecomp_pools_vr,:)&
                                         + matmul(ntranvert(1:ndecomp_pools_vr,:),matrix_Ninter_2d(:,:))*dt 
            if(end_of_year)then
               do i=1,ndecomp_pools_vr
                  if (abs(cs_soil%tran_acc(c,i,i)) .le. epsi)then !avoid inversion nan
                      cs_soil%tran_acc(c,i,i) = 1.e+36
                  end if 
               end do
               do i=1,ndecomp_pools_vr
                  if (abs(ns_soil%tran_nacc(c,i,i)) .le. epsi)then
                     ns_soil%tran_nacc(c,i,i) = 1.e+36
                  end if 
               end do
               call inverse(cs_soil%tran_acc(c,1:ndecomp_pools_vr,:),AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
               call inverse(ns_soil%tran_nacc(c,1:ndecomp_pools_vr,:),AKinvn(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ndecomp_pools_vr)
               soilmatrixc_cap(:,1) = -matmul(AKinv(1:ndecomp_pools_vr,1:ndecomp_pools_vr),cs_soil%in_acc(c,:))
               soilmatrixn_cap(:,1) = -matmul(AKinvn(1:ndecomp_pools_vr,1:ndecomp_pools_vr),ns_soil%in_nacc(c,:))
         
               do j = 1,nlevdecomp
                  do i=1,ndecomp_pools
                  !   if(soilmatrixc_cap(j+(i-1)*nlevdecomp,1) .le. 1.e-25)then ! for model stability
                  !      soilmatrixc_cap(j+(i-1)*nlevdecomp,1) = cs_soil%decomp_cpools_vr_col(c,j,i)
                  !      soilmatrixn_cap(j+(i-1)*nlevdecomp,1) = ns_soil%decomp_npools_vr_col(c,j,i)
                  !   end if
                     if(isspinup)then
                        cs_soil%decomp_cpools_vr_col(c,j,i) =  soilmatrixc_cap(j+(i-1)*nlevdecomp,1)
                        ns_soil%decomp_npools_vr_col(c,j,i) =  soilmatrixn_cap(j+(i-1)*nlevdecomp,1)
                     end if
                     cs_soil%matrix_cap_decomp_cpools_vr_col(c,j,i) = soilmatrixc_cap(j+(i-1)*nlevdecomp,1)
                     ns_soil%matrix_cap_decomp_npools_vr_col(c,j,i) = soilmatrixn_cap(j+(i-1)*nlevdecomp,1)
                     cs_soil%matrix_pot_decomp_cpools_vr_col(c,j,i) = soilmatrixc_cap(j+(i-1)*nlevdecomp,1) - cs_soil%decomp_cpools_vr_col(c,j,i)
                     ns_soil%matrix_pot_decomp_npools_vr_col(c,j,i) = soilmatrixn_cap(j+(i-1)*nlevdecomp,1) - ns_soil%decomp_npools_vr_col(c,j,i)
                  end do
               end do
               cs_soil%in_acc(c,:)     = 0.0_r8  
               cs_soil%tran_acc(c,:,:) = 0.0_r8 
               ns_soil%in_nacc(c,:)     = 0.0_r8  
               ns_soil%tran_nacc(c,:,:) = 0.0_r8 
            end if
         end if !is out_matrix
   enddo !fc 
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

