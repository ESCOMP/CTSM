module SoilBiogeochemDecompCascadeConType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Decomposition Cascade Type
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : ndecomp_cascade_transitions, ndecomp_pools, nlevdecomp, &
                              ndecomp_cascade_outtransitions
  use clm_varctl     , only : use_soil_matrixcn
  use SPMMod         , only : sparse_matrix_type, diag_matrix_type, vector_type
  !
  implicit none

  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_decomp_cascade_constants
  public :: InitSoilTransfer
  !
  type, public :: decomp_cascade_type
     !-- properties of each pathway along decomposition cascade 
     character(len=8)  , pointer :: cascade_step_name(:)               ! name of transition
     integer           , pointer :: cascade_donor_pool(:)              ! which pool is C taken from for a given decomposition step
     integer           , pointer :: cascade_receiver_pool(:)           ! which pool is C added to for a given decomposition step

     !-- properties of each decomposing pool
     logical           , pointer  :: floating_cn_ratio_decomp_pools(:) ! TRUE => pool has fixed C:N ratio
     character(len=8)  , pointer  :: decomp_pool_name_restart(:)       ! name of pool for restart files
     character(len=8)  , pointer  :: decomp_pool_name_history(:)       ! name of pool for history files
     character(len=20) , pointer  :: decomp_pool_name_long(:)          ! name of pool for netcdf long names
     character(len=8)  , pointer  :: decomp_pool_name_short(:)         ! name of pool for netcdf short names
     logical           , pointer  :: is_litter(:)                      ! TRUE => pool is a litter pool
     logical           , pointer  :: is_soil(:)                        ! TRUE => pool is a soil pool
     logical           , pointer  :: is_cwd(:)                         ! TRUE => pool is a cwd pool
     real(r8)          , pointer  :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
     real(r8)          , pointer  :: initial_stock(:)                  ! initial concentration for seeding at spinup
     real(r8)                     :: initial_stock_soildepth           ! soil depth for initial concentration for seeding at spinup
     logical           , pointer  :: is_metabolic(:)                   ! TRUE => pool is metabolic material
     logical           , pointer  :: is_cellulose(:)                   ! TRUE => pool is cellulose
     logical           , pointer  :: is_lignin(:)                      ! TRUE => pool is lignin
     real(r8)          , pointer  :: spinup_factor(:)                  ! factor by which to scale AD and relevant processes by

     integer,pointer :: spm_tranlist_a(:,:)                            ! Prescribed subscripts to map 2D variables (transitions,soil layer) to 1D sparse matrix format in a_ma_vr and na_ma_vr
     integer,pointer :: A_i(:)                                         ! Prescribed row number of all elements in a_ma_vr
     integer,pointer :: A_j(:)                                         ! Prescribed column number of all elements in na_ma_vr
     integer,pointer :: tri_i(:)                                       ! Prescribed row index of all entries in AVsoil
     integer,pointer :: tri_j(:)                                       ! Prescribed column index of all entries in AVsoil
     integer,pointer :: all_i(:)                                       ! Prescribed row index of all entries in AKallsoilc, AKallsoiln, AKXcacc, and AKXnacc
     integer,pointer :: all_j(:)                                       ! Prescribed column index of all entries in AKallsoilc, AKallsoiln, AKXcacc, and AKXnacc

     integer,pointer :: list_V_AKVfire (:)                             ! Saves mapping indices from V to (A*K+V-Kfire) in the addition subroutine SPMP_ABC
     integer,pointer :: list_AK_AKVfire(:)                             ! Saves mapping indices from A*K to (A*K+V-Kfire) in the addition subroutine SPMP_ABC
     integer,pointer :: list_fire_AKVfire(:)                           ! Saves mapping indices from Kfire to (A*K+V-Kfire) in the addition subroutine SPMP_ABC
     integer,pointer :: list_AK_AKV    (:)                             ! Saves mapping indices from A*K to (A*K+V) in the addition subroutine SPMP_AB
     integer,pointer :: list_V_AKV     (:)                             ! Saves mapping indices from V to (A*K+V) in the addition subroutine SPMP_AB
     integer,pointer :: list_Asoilc    (:)                             ! Saves mapping indices from a_ma_vr to AKsoilc
     integer,pointer :: list_Asoiln    (:)                             ! Saves mapping indices from na_ma_vr to AKsoiln

     integer, public :: n_all_entries                                  ! Number of all entries in AKallsoilc, AKallsoiln, AKXcacc, and AKXnacc
     integer, public :: Ntrans_setup                                   ! Number of horizontal transfers between soil and litter pools
     integer, public :: Ntri_setup                                     ! Number of non-zero entries in AVsoil

  end type decomp_cascade_type

  type(decomp_cascade_type), public :: decomp_cascade_con
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_decomp_cascade_constants()
    !
    ! !DESCRIPTION:
    ! Initialize decomposition cascade state
    !------------------------------------------------------------------------

    !-- properties of each pathway along decomposition cascade 
    allocate(decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions))

    ! NOTE(bja, 2015-10) according to Dave Lawrence and Charlie Koven,
    ! the indexing of decomposing pools from 0:ndecomp_pools is a
    ! bug. The lower bound should be 1. The index zero data shouldn't
    ! be used.
    
    !-- properties of each decomposing pool
    allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_litter(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_soil(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cwd(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_stock(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_metabolic(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cellulose(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_lignin(0:ndecomp_pools))
    allocate(decomp_cascade_con%spinup_factor(1:ndecomp_pools))
 
    if(use_soil_matrixcn)then
       allocate(decomp_cascade_con%spm_tranlist_a(1:nlevdecomp,1:ndecomp_cascade_transitions)); decomp_cascade_con%spm_tranlist_a(:,:) = -9999
       allocate(decomp_cascade_con%A_i(1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp));decomp_cascade_con%A_i(:) = -9999
       allocate(decomp_cascade_con%A_j(1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp));decomp_cascade_con%A_j(:) = -9999
       allocate(decomp_cascade_con%tri_i(1:(3*nlevdecomp-2)*(ndecomp_pools-1))); decomp_cascade_con%tri_i(:) = -9999
       allocate(decomp_cascade_con%tri_j(1:(3*nlevdecomp-2)*(ndecomp_pools-1))); decomp_cascade_con%tri_j(:) = -9999
    end if

    !-- properties of each pathway along decomposition cascade 
    decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions) = ''
    decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions) = 0
    decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions) = 0

    !-- first initialization of properties of each decomposing pool
    decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools) = .false.
    decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools)       = ''
    decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools)       = ''
    decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools)          = ''
    decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools)         = ''
    decomp_cascade_con%is_litter(0:ndecomp_pools)                      = .false.
    decomp_cascade_con%is_soil(0:ndecomp_pools)                        = .false.
    decomp_cascade_con%is_cwd(0:ndecomp_pools)                         = .false.
    decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools)               = nan
    decomp_cascade_con%initial_stock(0:ndecomp_pools)                  = nan
    decomp_cascade_con%initial_stock_soildepth                         = 0.3
    decomp_cascade_con%is_metabolic(0:ndecomp_pools)                   = .false.
    decomp_cascade_con%is_cellulose(0:ndecomp_pools)                   = .false.
    decomp_cascade_con%is_lignin(0:ndecomp_pools)                      = .false.
    decomp_cascade_con%spinup_factor(1:ndecomp_pools)                  = nan

    decomp_cascade_con%Ntrans_setup = (ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp
    decomp_cascade_con%Ntri_setup   = (3*nlevdecomp-2)*(ndecomp_pools - 1) !exclude one cwd 
  end subroutine init_decomp_cascade_constants

 subroutine InitSoilTransfer()
! Initialize sparse matrix variables and index. Count possible non-zero entries and record their x and y in the matrix.
! Collect those non-zero entry information, and save them into the list.

  use SPMMod         , only : sparse_matrix_type, diag_matrix_type, vector_type

  integer i,j,k,m,n
  integer,dimension(:) :: ntrans_per_donor(1:ndecomp_pools)
  real(r8),dimension(:) :: SM(1:1,1:decomp_cascade_con%Ntrans_setup),TRI(1:1,1:decomp_cascade_con%Ntri_setup)
  type(sparse_matrix_type) :: AK, AV, AKfire, AKallsoil!, AKtmp1, AKtmp2
  logical list_ready,init_readyAK
  integer num_soilc,filter_c(1:1)
  
  init_readyAK = .false.  
  ntrans_per_donor = 0
        
  do k = 1,ndecomp_cascade_transitions
     ntrans_per_donor(decomp_cascade_con%cascade_donor_pool(k)) &
               = ntrans_per_donor(decomp_cascade_con%cascade_donor_pool(k)) + 1
  end do

     k = 0
     n = 1
     do i = 1,ndecomp_pools
        do j=1,nlevdecomp
           do m=1,ntrans_per_donor(i)
              if(decomp_cascade_con%cascade_receiver_pool(m+k) .ne. 0)then
                 decomp_cascade_con%spm_tranlist_a(j,m+k) = n
                 decomp_cascade_con%A_i(n) = (decomp_cascade_con%cascade_receiver_pool(m+k)-1)*nlevdecomp+j
                 decomp_cascade_con%A_j(n) = (decomp_cascade_con%cascade_donor_pool(m+k)-1)*nlevdecomp+j
                 n = n + 1
              end if
           end do
        end do
        k = k + ntrans_per_donor(i)
     end do
     
     SM = 1._r8
     if(n-1 .ne. decomp_cascade_con%Ntrans_setup)then
        write(*,*),'error in InitSoilTransfer: number of transfers is error in count'
     end if

     n = 1
     do i = 1,ndecomp_pools
        do j = 1, nlevdecomp
           if(.not. decomp_cascade_con%is_cwd(i))then
              if (j > 1) then ! avoid tranfer from for example,soil1_1st layer to litr3_10th layer
                 TRI(1,n) = 1._r8
                 decomp_cascade_con%tri_j(n) = (i-1)*nlevdecomp + j
                 decomp_cascade_con%tri_i(n) = (i-1)*nlevdecomp + j - 1
                 n = n + 1
              end if
              TRI(1,n) = 1._r8
              decomp_cascade_con%tri_j(n) = (i-1)*nlevdecomp + j
              decomp_cascade_con%tri_i(n) = (i-1)*nlevdecomp + j
              n = n + 1
              if (j < nlevdecomp) then ! avoid tranfer from for example, litr3_10th layer to soil1_1st layer
                 TRI(1,n) = 1._r8
                 decomp_cascade_con%tri_j(n) = (i-1)*nlevdecomp + j
                 decomp_cascade_con%tri_i(n) = (i-1)*nlevdecomp + j + 1
                 n = n + 1
              end if
           end if
        end do
     end do

     if(n-1 .ne. decomp_cascade_con%Ntri_setup)then
        write(*,*),'error in InitSoilTransfer: number of vertical-transfers is error in count'
     end if

     num_soilc = 1
     filter_c(1:1) = 1
     call AK%InitSM(ndecomp_pools*nlevdecomp,1,1)
     call AK%SetValueA(1,1,num_soilc,filter_c,SM,decomp_cascade_con%A_i,decomp_cascade_con%A_j,decomp_cascade_con%Ntrans_setup,init_readyAK)
     allocate(decomp_cascade_con%list_AK_AKVfire(1:AK%NE))
     allocate(decomp_cascade_con%list_AK_AKV    (1:AK%NE))

     call AV%InitSM(ndecomp_pools*nlevdecomp,1,1)
     call AV%SetValueSM(1,1,num_soilc,filter_c,TRI,decomp_cascade_con%tri_i,decomp_cascade_con%tri_j,decomp_cascade_con%Ntri_setup)
     allocate(decomp_cascade_con%list_V_AKVfire (1:AV%NE))
     allocate(decomp_cascade_con%list_V_AKV     (1:AV%NE))
     
     call AKfire%InitSM(ndecomp_pools*nlevdecomp,1,1)
     call AKfire%SetValueA_diag(num_soilc,filter_c,1._r8)
     allocate(decomp_cascade_con%list_fire_AKVfire(1:AKfire%NE))
   
     list_ready = .false.
     call AKallsoil%InitSM(ndecomp_pools*nlevdecomp,1,1)
     call AKallsoil%SPMP_ABC(num_soilc,filter_c,AK,AV,AKfire,list_ready)

     decomp_cascade_con%n_all_entries = AKallsoil%NE
     allocate(decomp_cascade_con%all_i(1:decomp_cascade_con%n_all_entries))
     allocate(decomp_cascade_con%all_j(1:decomp_cascade_con%n_all_entries))
     decomp_cascade_con%all_i = AKallsoil%RI
     decomp_cascade_con%all_j = AKallsoil%CI

     allocate(decomp_cascade_con%list_Asoilc    (1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp))
     allocate(decomp_cascade_con%list_Asoiln    (1:(ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp))

   end subroutine InitSoilTransfer

end module SoilBiogeochemDecompCascadeConType
