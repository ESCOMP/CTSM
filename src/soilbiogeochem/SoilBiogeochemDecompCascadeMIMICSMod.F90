module SoilBiogeochemDecompCascadeMIMICSMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Sets the coeffiecients used in the decomposition cascade submodel.  
  ! This uses the MIMICS parameters.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_const_mod                      , only : SHR_CONST_TKFRZ
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : nlevdecomp, ndecomp_pools_max
  use clm_varpar                         , only : i_phys_som, i_chem_som, i_str_lit, i_met_lit, i_cop_mic, i_oli_mic, i_cwd
  use clm_varpar                         , only : i_litr_min, i_litr_max, i_cwdl2
  use clm_varctl                         , only : iulog, spinup_state, anoxia, use_lch4, use_fates
  use clm_varcon                         , only : zsoi
  use decompMod                          , only : bounds_type
  use spmdMod                            , only : masterproc
  use abortutils                         , only : endrun
  use CNSharedParamsMod                  , only : CNParamsShareInst, nlev_soildecomp_standard 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, InitSoilTransfer, use_soil_matrixcn
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilStateType                      , only : soilstate_type
  use TemperatureType                    , only : temperature_type 
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use ch4Mod                             , only : ch4_type
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  use CLMFatesInterfaceMod               , only : hlm_fates_interface_type

  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams                      ! Read in parameters from params file
  public :: init_decompcascade_mimics       ! Initialization
  public :: decomp_rates_mimics             ! Figure out decomposition rates
  !
  ! !PUBLIC DATA MEMBERS 
  !
  ! !PRIVATE DATA MEMBERS 
  ! next four 2d vars are dimensioned (columns,nlevdecomp)
  real(r8), private, allocatable :: desorp(:,:)
  real(r8), private, allocatable :: fphys_m1(:,:)
  real(r8), private, allocatable :: fphys_m2(:,:)
  real(r8), private, allocatable :: p_scalar(:,:)
  integer, private :: i_avl_som  ! index of available (aka active) SOM
  integer, private :: i_l1m1  ! indices of transitions, eg l1m1: litter 1 -> first microbial pool
  integer, private :: i_l1m2
  integer, private :: i_l2m1
  integer, private :: i_l2m2
  integer, private :: i_s1m1
  integer, private :: i_s1m2
  integer, private :: i_m1s1
  integer, private :: i_m1s2
  integer, private :: i_m2s1
  integer, private :: i_m2s2
  integer, private :: i_s2s1
  integer, private :: i_s3s1
  integer, private :: i_m1s3
  integer, private :: i_m2s3
  real(r8), private :: rf_l1m1  ! respiration fractions by transition
  real(r8), private :: rf_l1m2
  real(r8), private :: rf_l2m1
  real(r8), private :: rf_l2m2
  real(r8), private :: rf_s1m1
  real(r8), private :: rf_s1m2
  real(r8), private :: vint_l1_m1  ! regression intercepts by transition
  real(r8), private :: vint_l2_m1
  real(r8), private :: vint_s1_m1
  real(r8), private :: vint_l1_m2
  real(r8), private :: vint_l2_m2
  real(r8), private :: vint_s1_m2
  real(r8), private :: kint_l1_m1  ! regression intercepts by transition
  real(r8), private :: kint_l2_m1
  real(r8), private :: kint_s1_m1
  real(r8), private :: kint_l1_m2
  real(r8), private :: kint_l2_m2
  real(r8), private :: kint_s1_m2
  real(r8), private :: vmod_l1_m1  ! vmod = vmod * av from Wieder et al 2015
  real(r8), private :: vmod_l2_m1
  real(r8), private :: vmod_s1_m1
  real(r8), private :: vmod_l1_m2
  real(r8), private :: vmod_l2_m2
  real(r8), private :: vmod_s1_m2
  real(r8), private :: kmod_l1_m1  ! kmod = ak / kmod from Wieder et al 2015
  real(r8), private :: kmod_l2_m1
  real(r8), private :: kmod_s1_m1
  real(r8), private :: kmod_l1_m2
  real(r8), private :: kmod_l2_m2
  real(r8), private :: kmod_s1_m2
  real(r8), private :: vslope_l1_m1  ! regression coefficients by transition
  real(r8), private :: vslope_l2_m1
  real(r8), private :: vslope_s1_m1
  real(r8), private :: vslope_l1_m2
  real(r8), private :: vslope_l2_m2
  real(r8), private :: vslope_s1_m2
  real(r8), private :: kslope_l1_m1  ! regression coefficients by transition
  real(r8), private :: kslope_l2_m1
  real(r8), private :: kslope_s1_m1
  real(r8), private :: kslope_l1_m2
  real(r8), private :: kslope_l2_m2
  real(r8), private :: kslope_s1_m2

  type, private :: params_type
     real(r8) :: mimics_nue_into_mic  ! microbial N use efficiency for N fluxes
     real(r8) :: mimics_desorpQ10
     real(r8) :: mimics_densdep  ! exponent controling the density dependence of microbial turnover
     real(r8) :: mimics_tau_mod_factor  ! (1 / tauModDenom) from testbed code
     real(r8) :: mimics_tau_mod_min
     real(r8) :: mimics_tau_mod_max
     real(r8) :: mimics_ko_r
     real(r8) :: mimics_ko_k
     real(r8) :: mimics_cn_r  ! C:N of MICr
     real(r8) :: mimics_cn_k  ! C:N of MICk
     real(r8) :: mimics_cn_mod_num  ! adjusts microbial CN based on fmet
     real(r8) :: mimics_t_soi_ref  ! reference soil temperature (degC)
     real(r8) :: mimics_initial_Cstocks_depth  ! Soil depth for initial C stocks for a cold-start (m)
     real(r8), allocatable :: mimics_initial_Cstocks(:)  ! Initial C stocks for a cold-start (gC/m3)
     ! The next few vectors are dimensioned by the number of decomposition
     ! transitions that make use of the corresponding parameters, currently
     ! six. The transitions are represented in this order:
     ! l1m1 l2m1 s1m1 l1m2 l2m2 s1m2
     real(r8), allocatable :: mimics_mge(:)  ! Microbial growth efficiency (mg/mg)
     real(r8), allocatable :: mimics_vmod(:)  ! vmod = vmod * av from Wieder et al 2015
     real(r8), allocatable :: mimics_vint(:)  ! regression intercepts (5.47 ln(mg Cs (mg MIC)-1 h-1) )
     real(r8), allocatable :: mimics_vslope(:)  ! regression coeffs (ln(mg Cs (mg MIC)-1 h-1) ¡C-1)
     real(r8), allocatable :: mimics_kmod(:)  ! kmod = ak / kmod from Wieder et al 2015
     real(r8), allocatable :: mimics_kint(:)  ! regression intercepts
     real(r8), allocatable :: mimics_kslope(:)  ! regression coeffs
     ! The next few vectors are dimensioned by the number of parameters with
     ! the same name (eg 2 for mimics_tau_r_p1, mimics_tau_r_p2) used in the
     ! respective formula. In the formulas, we use scalar copies of each
     ! parameter with suffixes _p1, p2, p3, ... to distinguish among them.
     ! See allocate statements below for the size of each of the following
     ! vectors.
     real(r8), allocatable :: mimics_fmet(:)
     real(r8), allocatable :: mimics_p_scalar(:)
     real(r8), allocatable :: mimics_fphys_r(:)
     real(r8), allocatable :: mimics_fphys_k(:)
     real(r8), allocatable :: mimics_fchem_r(:)
     real(r8), allocatable :: mimics_fchem_k(:)
     real(r8), allocatable :: mimics_desorp(:)
     real(r8), allocatable :: mimics_tau_r(:)
     real(r8), allocatable :: mimics_tau_k(:)
  end type params_type
  !
  type(params_type), private :: params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'readMimicsParams'
    character(len=100) :: errCode = 'Error reading MIMICS params '
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! Read off of netcdf file
    tString='mimics_initial_Cstocks_depth'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_initial_Cstocks_depth=tempr

    allocate(params_inst%mimics_initial_Cstocks(ndecomp_pools_max))
    tString='mimics_initial_Cstocks'
    call ncd_io(trim(tString), params_inst%mimics_initial_Cstocks(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_mge(ndecomp_pools_max))
    tString='mimics_mge'
    call ncd_io(trim(tString), params_inst%mimics_mge(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_vmod(ndecomp_pools_max))
    tString='mimics_vmod'
    call ncd_io(trim(tString), params_inst%mimics_vmod(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_vslope(ndecomp_pools_max))
    tString='mimics_vslope'
    call ncd_io(trim(tString), params_inst%mimics_vslope(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_vint(ndecomp_pools_max))
    tString='mimics_vint'
    call ncd_io(trim(tString), params_inst%mimics_vint(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_kmod(ndecomp_pools_max))
    tString='mimics_kmod'
    call ncd_io(trim(tString), params_inst%mimics_kmod(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_kslope(ndecomp_pools_max))
    tString='mimics_kslope'
    call ncd_io(trim(tString), params_inst%mimics_kslope(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_kint(ndecomp_pools_max))
    tString='mimics_kint'
    call ncd_io(trim(tString), params_inst%mimics_kint(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_p_scalar(2))
    tString='mimics_p_scalar'
    call ncd_io(trim(tString), params_inst%mimics_p_scalar(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_desorp(2))
    tString='mimics_desorp'
    call ncd_io(trim(tString), params_inst%mimics_desorp(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_fphys_r(2))
    tString='mimics_fphys_r'
    call ncd_io(trim(tString), params_inst%mimics_fphys_r(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_fphys_k(2))
    tString='mimics_fphys_k'
    call ncd_io(trim(tString), params_inst%mimics_fphys_k(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_fmet(4))
    tString='mimics_fmet'
    call ncd_io(trim(tString), params_inst%mimics_fmet(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_fchem_r(2))
    tString='mimics_fchem_r'
    call ncd_io(trim(tString), params_inst%mimics_fchem_r(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_fchem_k(2))
    tString='mimics_fchem_k'
    call ncd_io(trim(tString), params_inst%mimics_fchem_k(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_tau_r(2))
    tString='mimics_tau_r'
    call ncd_io(trim(tString), params_inst%mimics_tau_r(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mimics_tau_k(2))
    tString='mimics_tau_k'
    call ncd_io(trim(tString), params_inst%mimics_tau_k(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    tString='mimics_nue_into_mic'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_nue_into_mic = tempr

    tString='mimics_tau_mod_factor'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_tau_mod_factor = tempr

    tString='mimics_tau_mod_min'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_tau_mod_min = tempr

    tString='mimics_tau_mod_max'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_tau_mod_max = tempr

    tString='mimics_ko_r'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_ko_r = tempr

    tString='mimics_ko_k'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_ko_k = tempr

    tString='mimics_densdep'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_densdep = tempr

    tString='mimics_desorpQ10'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_desorpQ10 = tempr

    tString='mimics_t_soi_ref'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_t_soi_ref = tempr

    tString='mimics_cn_mod_num'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_cn_mod_num = tempr

    tString='mimics_cn_r'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_cn_r = tempr

    tString='mimics_cn_k'
    call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%mimics_cn_k = tempr

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_mimics(bounds, soilbiogeochem_state_inst, soilstate_inst )
    !
    ! !DESCRIPTION:
    ! initialize rate constants and decomposition pathways following the
    ! decomposition cascade of the MIMICS model.
    !
    ! !USES:
    use clm_varcon, only: pct_to_frac
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds  
    type(soilbiogeochem_state_type) , intent(inout) :: soilbiogeochem_state_inst
    type(soilstate_type)            , intent(in)    :: soilstate_inst
    !
    ! !LOCAL VARIABLES
    !-- properties of each decomposing pool
    real(r8) :: mimics_nue_into_mic
    real(r8) :: mimics_p_scalar_p1
    real(r8) :: mimics_p_scalar_p2
    real(r8) :: mimics_fphys_r_p1
    real(r8) :: mimics_fphys_r_p2
    real(r8) :: mimics_fphys_k_p1
    real(r8) :: mimics_fphys_k_p2
    real(r8) :: mimics_desorp_p1
    real(r8) :: mimics_desorp_p2

    real(r8):: speedup_fac                  ! acceleration factor, higher when vertsoilc = .true.

    real(r8) :: clay_frac  ! local copy of cellclay converted from % (fraction)
    integer  :: c, j  ! indices
    !-----------------------------------------------------------------------

    associate(                                                                                     &
         nue_decomp_cascade             => soilbiogeochem_state_inst%nue_decomp_cascade_col      , & ! Output: [real(r8)          (:)     ]  N use efficiency for a given transition (gN going into microbe / gN decomposed)

         cellclay                       => soilstate_inst%cellclay_col                           , & ! Input:  [real(r8)          (:,:)   ]  column 3D clay (%)
         
         cascade_donor_pool             => decomp_cascade_con%cascade_donor_pool                 , & ! Output: [integer           (:)     ]  which pool is C taken from for a given decomposition step 
         cascade_receiver_pool          => decomp_cascade_con%cascade_receiver_pool              , & ! Output: [integer           (:)     ]  which pool is C added to for a given decomposition step   
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools     , & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:N ratio                          
         is_microbe                     => decomp_cascade_con%is_microbe                         , & ! Output: [logical           (:)     ]  TRUE => pool is microbial
         is_litter                      => decomp_cascade_con%is_litter                          , & ! Output: [logical           (:)     ]  TRUE => pool is a litter pool                             
         is_soil                        => decomp_cascade_con%is_soil                            , & ! Output: [logical           (:)     ]  TRUE => pool is a soil pool                               
         is_cwd                         => decomp_cascade_con%is_cwd                             , & ! Output: [logical           (:)     ]  TRUE => pool is a cwd pool                                
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio                   , & ! Output: [real(r8)          (:)     ]  c:n ratio for initialization of pools                    
         initial_stock                  => decomp_cascade_con%initial_stock                      , & ! Output: [real(r8)          (:)     ]  initial concentration for seeding at spinup              
         initial_stock_soildepth        => decomp_cascade_con%initial_stock_soildepth            , & ! Output: [real(r8)          (:)     ]  soil depth for initial concentration for seeding at spinup              
         is_metabolic                   => decomp_cascade_con%is_metabolic                       , & ! Output: [logical           (:)     ]  TRUE => pool is metabolic material                        
         is_cellulose                   => decomp_cascade_con%is_cellulose                       , & ! Output: [logical           (:)     ]  TRUE => pool is cellulose                                 
         is_lignin                      => decomp_cascade_con%is_lignin                          , & ! Output: [logical           (:)     ]  TRUE => pool is lignin                                    
         spinup_factor                  => decomp_cascade_con%spinup_factor                        & ! Output: [real(r8)          (:)     ]  factor for AD spinup associated with each pool           

         )

      allocate(desorp(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(fphys_m1(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(fphys_m2(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(p_scalar(bounds%begc:bounds%endc,1:nlevdecomp))

      !------- time-constant coefficients ---------- !
      mimics_nue_into_mic = params_inst%mimics_nue_into_mic
      mimics_p_scalar_p1 = params_inst%mimics_p_scalar(1)
      mimics_p_scalar_p2 = params_inst%mimics_p_scalar(2)
      mimics_fphys_r_p1 = params_inst%mimics_fphys_r(1)
      mimics_fphys_r_p2 = params_inst%mimics_fphys_r(2)
      mimics_fphys_k_p1 = params_inst%mimics_fphys_k(1)
      mimics_fphys_k_p2 = params_inst%mimics_fphys_k(2)
      mimics_desorp_p1 = params_inst%mimics_desorp(1)
      mimics_desorp_p2 = params_inst%mimics_desorp(2)

      ! set respiration fractions for fluxes between compartments
      rf_l1m1 = 1.0_r8 - params_inst%mimics_mge(1)
      rf_l2m1 = 1.0_r8 - params_inst%mimics_mge(2)
      rf_s1m1 = 1.0_r8 - params_inst%mimics_mge(3)
      rf_l1m2 = 1.0_r8 - params_inst%mimics_mge(4)
      rf_l2m2 = 1.0_r8 - params_inst%mimics_mge(5)
      rf_s1m2 = 1.0_r8 - params_inst%mimics_mge(6)

      ! vmod = "old" vmod * av  AND  kmod = ak / "old" kmod
      ! Table B1 Wieder et al. (2015) and MIMICS params file give diff
      ! ak and av values. I used the params file values.
      vmod_l1_m1 = params_inst%mimics_vmod(1)
      vmod_l2_m1 = params_inst%mimics_vmod(2)
      vmod_s1_m1 = params_inst%mimics_vmod(3)
      vmod_l1_m2 = params_inst%mimics_vmod(4)
      vmod_l2_m2 = params_inst%mimics_vmod(5)
      vmod_s1_m2 = params_inst%mimics_vmod(6)
      kmod_l1_m1 = params_inst%mimics_kmod(1)
      kmod_l2_m1 = params_inst%mimics_kmod(2)
      kmod_s1_m1 = params_inst%mimics_kmod(3)
      kmod_l1_m2 = params_inst%mimics_kmod(4)
      kmod_l2_m2 = params_inst%mimics_kmod(5)
      kmod_s1_m2 = params_inst%mimics_kmod(6)
      vslope_l1_m1 = params_inst%mimics_vslope(1)
      vslope_l2_m1 = params_inst%mimics_vslope(2)
      vslope_s1_m1 = params_inst%mimics_vslope(3)
      vslope_l1_m2 = params_inst%mimics_vslope(4)
      vslope_l2_m2 = params_inst%mimics_vslope(5)
      vslope_s1_m2 = params_inst%mimics_vslope(6)
      kslope_l1_m1 = params_inst%mimics_kslope(1)
      kslope_l2_m1 = params_inst%mimics_kslope(2)
      kslope_s1_m1 = params_inst%mimics_kslope(3)
      kslope_l1_m2 = params_inst%mimics_kslope(4)
      kslope_l2_m2 = params_inst%mimics_kslope(5)
      kslope_s1_m2 = params_inst%mimics_kslope(6)
      vint_l1_m1 = params_inst%mimics_vint(1)
      vint_l2_m1 = params_inst%mimics_vint(2)
      vint_s1_m1 = params_inst%mimics_vint(3)
      vint_l1_m2 = params_inst%mimics_vint(4)
      vint_l2_m2 = params_inst%mimics_vint(5)
      vint_s1_m2 = params_inst%mimics_vint(6)
      kint_l1_m1 = params_inst%mimics_kint(1)
      kint_l2_m1 = params_inst%mimics_kint(2)
      kint_s1_m1 = params_inst%mimics_kint(3)
      kint_l1_m2 = params_inst%mimics_kint(4)
      kint_l2_m2 = params_inst%mimics_kint(5)
      kint_s1_m2 = params_inst%mimics_kint(6)

      ! some of these are dependent on the soil texture properties
      ! One-time initializations here.
      ! Time-dep params in subr. decomp_rates_mimics.

      do c = bounds%begc, bounds%endc
         do j = 1, nlevdecomp
            ! The parameter values currently in the params files always lead to
            ! positive values for the expressions below, so we do not
            ! need to use the max function to limit these expressions.
            ! We apply the min function on cellclay because we are looping over
            ! some non-soil columns here that contain cellclay = 1e36.
            clay_frac = pct_to_frac * &
                        dmin1(100.0_r8, cellclay(c,j))  ! conv. % to fraction
            desorp(c,j) = mimics_desorp_p1 * dexp(mimics_desorp_p2 * clay_frac)
            fphys_m1(c,j) = dmin1(1.0_r8, mimics_fphys_r_p1 * &
                                          dexp(mimics_fphys_r_p2 * clay_frac))
            fphys_m2(c,j) = dmin1(1.0_r8, mimics_fphys_k_p1 * &
                                          dexp(mimics_fphys_k_p2 * clay_frac))
            p_scalar(c,j) = 1.0_r8 / (mimics_p_scalar_p1 * &
                                      dexp(mimics_p_scalar_p2 * dsqrt(clay_frac)))
         end do
      end do
      initial_stock_soildepth = params_inst%mimics_initial_Cstocks_depth

      !-------------------  list of pools and their attributes  ------------
      i_litr_min = 1
      i_met_lit = i_litr_min
      floating_cn_ratio_decomp_pools(i_met_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_met_lit) = 'litr1'
      decomp_cascade_con%decomp_pool_name_history(i_met_lit) = 'LIT_MET'
      decomp_cascade_con%decomp_pool_name_long(i_met_lit) = 'metabolic litter'
      decomp_cascade_con%decomp_pool_name_short(i_met_lit) = 'L1'
      is_microbe(i_met_lit) = .false.
      is_litter(i_met_lit) = .true.
      is_soil(i_met_lit) = .false.
      is_cwd(i_met_lit) = .false.
      initial_cn_ratio(i_met_lit) = 10._r8  ! 90 in BGC; not used in MIMICS
      initial_stock(i_met_lit) = params_inst%mimics_initial_Cstocks(i_met_lit)
      is_metabolic(i_met_lit) = .true.
      is_cellulose(i_met_lit) = .false.
      is_lignin(i_met_lit) = .false.

      i_str_lit = i_met_lit + 1
      floating_cn_ratio_decomp_pools(i_str_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_str_lit) = 'litr2'
      decomp_cascade_con%decomp_pool_name_history(i_str_lit) = 'LIT_STR'
      decomp_cascade_con%decomp_pool_name_long(i_str_lit) = 'structural litter'
      decomp_cascade_con%decomp_pool_name_short(i_str_lit) = 'L2'
      is_microbe(i_str_lit) = .false.
      is_litter(i_str_lit) = .true.
      is_soil(i_str_lit) = .false.
      is_cwd(i_str_lit) = .false.
      initial_cn_ratio(i_str_lit) = 10._r8  ! 90 in BGC; not used in MIMICS
      initial_stock(i_str_lit) = params_inst%mimics_initial_Cstocks(i_str_lit)
      is_metabolic(i_str_lit) = .false.
      is_cellulose(i_str_lit) = .true.
      is_lignin(i_str_lit) = .true.

      i_litr_max = i_str_lit
      if (i_litr_min /= 1 .or. i_litr_max < 2 .or. i_litr_max > 3) then
         write(iulog,*) 'Expecting i_litr_min = 1 and i_litr_max = 2 or 3.'
         write(iulog,*) 'See pftconMod, SoilBiogeochemCarbonFluxType, and'
         write(iulog,*) 'clmfates_interfaceMod for ramifications of changing'
         write(iulog,*) 'this assumption.'
         call endrun(msg='ERROR: i_litr_min and/or i_litr_max out of range '// &
              errMsg(sourcefile, __LINE__))
      end if

      i_avl_som = i_str_lit + 1
      floating_cn_ratio_decomp_pools(i_avl_som) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_avl_som) = 'soil1'
      decomp_cascade_con%decomp_pool_name_history(i_avl_som) = 'SOM_AVL'
      decomp_cascade_con%decomp_pool_name_long(i_avl_som) = 'available soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_avl_som) = 'S1'
      is_microbe(i_avl_som) = .false.
      is_litter(i_avl_som) = .false.
      is_soil(i_avl_som) = .true.
      is_cwd(i_avl_som) = .false.
      initial_cn_ratio(i_avl_som) = 10._r8  ! cn_s1 in BGC; not used in MIMICS
      initial_stock(i_avl_som) = params_inst%mimics_initial_Cstocks(i_avl_som)
      is_metabolic(i_avl_som) = .false.
      is_cellulose(i_avl_som) = .false.
      is_lignin(i_avl_som) = .false.

      i_chem_som = i_avl_som + 1
      floating_cn_ratio_decomp_pools(i_chem_som) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_chem_som) = 'soil2'
      decomp_cascade_con%decomp_pool_name_history(i_chem_som) = 'SOM_CHEM'
      decomp_cascade_con%decomp_pool_name_long(i_chem_som) = 'chemically protected soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_chem_som) = 'S2'
      is_microbe(i_chem_som) = .false.
      is_litter(i_chem_som) = .false.
      is_soil(i_chem_som) = .true.
      is_cwd(i_chem_som) = .false.
      initial_cn_ratio(i_chem_som) = 10._r8  ! cn_s2 in BGC; not used in MIMICS
      initial_stock(i_chem_som) = params_inst%mimics_initial_Cstocks(i_chem_som)
      is_metabolic(i_chem_som) = .false.
      is_cellulose(i_chem_som) = .false.
      is_lignin(i_chem_som) = .false.

      i_phys_som = i_chem_som + 1
      floating_cn_ratio_decomp_pools(i_phys_som) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_phys_som) = 'soil3'
      decomp_cascade_con%decomp_pool_name_history(i_phys_som) = 'SOM_PHYS'
      decomp_cascade_con%decomp_pool_name_long(i_phys_som) = 'physically protected soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_phys_som) = 'S3'
      is_microbe(i_phys_som) = .false.
      is_litter(i_phys_som) = .false.
      is_soil(i_phys_som) = .true.
      is_cwd(i_phys_som) = .false.
      initial_cn_ratio(i_phys_som) = 10._r8  ! cn_s3 in BGC; not used in MIMICS
      initial_stock(i_phys_som) = params_inst%mimics_initial_Cstocks(i_phys_som)
      is_metabolic(i_phys_som) = .false.
      is_cellulose(i_phys_som) = .false.
      is_lignin(i_phys_som) = .false.

      i_cop_mic = i_phys_som + 1
      floating_cn_ratio_decomp_pools(i_cop_mic) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_cop_mic) = 'micr1'
      decomp_cascade_con%decomp_pool_name_history(i_cop_mic) = 'MIC_COP'
      decomp_cascade_con%decomp_pool_name_long(i_cop_mic) = 'copiotrophic microbes'
      decomp_cascade_con%decomp_pool_name_short(i_cop_mic) = 'M1'
      is_microbe(i_cop_mic) = .true.
      is_litter(i_cop_mic) = .false.
      is_soil(i_cop_mic) = .false.
      is_cwd(i_cop_mic) = .false.
      initial_cn_ratio(i_cop_mic) = 10._r8  ! MIMICS may use this
      initial_stock(i_cop_mic) = params_inst%mimics_initial_Cstocks(i_cop_mic)
      is_metabolic(i_cop_mic) = .false.
      is_cellulose(i_cop_mic) = .false.
      is_lignin(i_cop_mic) = .false.

      i_oli_mic = i_cop_mic + 1
      floating_cn_ratio_decomp_pools(i_oli_mic) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_oli_mic) = 'micr2'
      decomp_cascade_con%decomp_pool_name_history(i_oli_mic) = 'MIC_OLI'
      decomp_cascade_con%decomp_pool_name_long(i_oli_mic) = 'oligotrophic microbes'
      decomp_cascade_con%decomp_pool_name_short(i_oli_mic) = 'M2'
      is_microbe(i_oli_mic) = .true.
      is_litter(i_oli_mic) = .false.
      is_soil(i_oli_mic) = .false.
      is_cwd(i_oli_mic) = .false.
      initial_cn_ratio(i_oli_mic) = 10._r8  ! MIMICS may use this
      initial_stock(i_oli_mic) = params_inst%mimics_initial_Cstocks(i_oli_mic)
      is_metabolic(i_oli_mic) = .false.
      is_cellulose(i_oli_mic) = .false.
      is_lignin(i_oli_mic) = .false.

      if (.not. use_fates) then
         ! CWD
         i_cwd = i_oli_mic + 1
         floating_cn_ratio_decomp_pools(i_cwd) = .true.
         decomp_cascade_con%decomp_pool_name_restart(i_cwd) = 'cwd'
         decomp_cascade_con%decomp_pool_name_history(i_cwd) = 'CWD'
         decomp_cascade_con%decomp_pool_name_long(i_cwd) = 'coarse woody debris'
         decomp_cascade_con%decomp_pool_name_short(i_cwd) = 'CWD'
         is_microbe(i_cwd) = .false.
         is_litter(i_cwd) = .false.
         is_soil(i_cwd) = .false.
         is_cwd(i_cwd) = .true.
         initial_cn_ratio(i_cwd) = 10._r8  ! 90 in BGC; not used in MIMICS
         initial_stock(i_cwd) = params_inst%mimics_initial_Cstocks(i_cwd)
         is_metabolic(i_cwd) = .false.
         is_cellulose(i_cwd) = .false.
         is_lignin(i_cwd) = .false.
      endif

      speedup_fac = 1._r8

      !lit1,2
      spinup_factor(i_met_lit) = 1._r8
      spinup_factor(i_str_lit) = 1._r8
      !CWD
      if (.not. use_fates) then
         spinup_factor(i_cwd) = max(1._r8, (speedup_fac * CNParamsShareInst%tau_cwd * 0.5_r8 ))
      end if
      !som1,2,3
      spinup_factor(i_avl_som) = 1._r8
      spinup_factor(i_chem_som) = 1._r8  ! BGC used cwd formula above but
      spinup_factor(i_phys_som) = 1._r8  ! ...w the respective tau_s values
      ! micr1,2
      spinup_factor(i_cop_mic) = 1._r8
      spinup_factor(i_oli_mic) = 1._r8

      if ( masterproc ) then
         write(iulog,*) 'Spinup_state ',spinup_state
         write(iulog,*) 'Spinup factors ',spinup_factor
      end if

      !----------------  list of transitions and their time-independent coefficients  ---------------!
      i_l1m1 = 1
      decomp_cascade_con%cascade_step_name(i_l1m1) = 'L1M1'
      cascade_donor_pool(i_l1m1) = i_met_lit
      cascade_receiver_pool(i_l1m1) = i_cop_mic
      nue_decomp_cascade(i_l1m1) = mimics_nue_into_mic

      i_l1m2 = 2
      decomp_cascade_con%cascade_step_name(i_l1m2) = 'L1M2'
      cascade_donor_pool(i_l1m2) = i_met_lit
      cascade_receiver_pool(i_l1m2) = i_oli_mic
      nue_decomp_cascade(i_l1m2) = mimics_nue_into_mic

      i_l2m1 = 3
      decomp_cascade_con%cascade_step_name(i_l2m1) = 'L2M1'
      cascade_donor_pool(i_l2m1) = i_str_lit
      cascade_receiver_pool(i_l2m1) = i_cop_mic
      nue_decomp_cascade(i_l2m1) = mimics_nue_into_mic

      i_l2m2 = 4
      decomp_cascade_con%cascade_step_name(i_l2m2) = 'L2M2'
      cascade_donor_pool(i_l2m2) = i_str_lit
      cascade_receiver_pool(i_l2m2) = i_oli_mic
      nue_decomp_cascade(i_l2m2) = mimics_nue_into_mic

      i_s1m1 = 5
      decomp_cascade_con%cascade_step_name(i_s1m1) = 'S1M1'
      cascade_donor_pool(i_s1m1) = i_avl_som
      cascade_receiver_pool(i_s1m1) = i_cop_mic
      nue_decomp_cascade(i_s1m1) = mimics_nue_into_mic

      i_s1m2 = 6
      decomp_cascade_con%cascade_step_name(i_s1m2) = 'S1M2'
      cascade_donor_pool(i_s1m2) = i_avl_som
      cascade_receiver_pool(i_s1m2) = i_oli_mic
      nue_decomp_cascade(i_s1m2) = mimics_nue_into_mic

      i_s2s1 = 7
      decomp_cascade_con%cascade_step_name(i_s2s1) = 'S2S1'
      cascade_donor_pool(i_s2s1) = i_chem_som
      cascade_receiver_pool(i_s2s1) = i_avl_som
      nue_decomp_cascade(i_s2s1) = 1.0_r8

      i_s3s1 = 8
      decomp_cascade_con%cascade_step_name(i_s3s1) = 'S3S1'
      cascade_donor_pool(i_s3s1) = i_phys_som
      cascade_receiver_pool(i_s3s1) = i_avl_som
      nue_decomp_cascade(i_s3s1) = 1.0_r8

      i_m1s1 = 9
      decomp_cascade_con%cascade_step_name(i_m1s1) = 'M1S1'
      cascade_donor_pool(i_m1s1) = i_cop_mic
      cascade_receiver_pool(i_m1s1) = i_avl_som
      nue_decomp_cascade(i_m1s1) = 1.0_r8

      i_m1s2 = 10
      decomp_cascade_con%cascade_step_name(i_m1s2) = 'M1S2'
      cascade_donor_pool(i_m1s2) = i_cop_mic
      cascade_receiver_pool(i_m1s2) = i_chem_som
      nue_decomp_cascade(i_m1s2) = 1.0_r8

      i_m1s3 = 11
      decomp_cascade_con%cascade_step_name(i_m1s3) = 'M1S3'
      cascade_donor_pool(i_m1s3) = i_cop_mic
      cascade_receiver_pool(i_m1s3) = i_phys_som
      nue_decomp_cascade(i_m1s3) = 1.0_r8

      i_m2s1 = 12
      decomp_cascade_con%cascade_step_name(i_m2s1) = 'M2S1'
      cascade_donor_pool(i_m2s1) = i_oli_mic
      cascade_receiver_pool(i_m2s1) = i_avl_som
      nue_decomp_cascade(i_m2s1) = 1.0_r8

      i_m2s2 = 13
      decomp_cascade_con%cascade_step_name(i_m2s2) = 'M2S2'
      cascade_donor_pool(i_m2s2) = i_oli_mic
      cascade_receiver_pool(i_m2s2) = i_chem_som
      nue_decomp_cascade(i_m2s2) = 1.0_r8

      i_m2s3 = 14
      decomp_cascade_con%cascade_step_name(i_m2s3) = 'M2S3'
      cascade_donor_pool(i_m2s3) = i_oli_mic
      cascade_receiver_pool(i_m2s3) = i_phys_som
      nue_decomp_cascade(i_m2s3) = 1.0_r8

      if (.not. use_fates) then
         i_cwdl2 = 15
         decomp_cascade_con%cascade_step_name(i_cwdl2) = 'CWDL2'
         cascade_donor_pool(i_cwdl2) = i_cwd
         cascade_receiver_pool(i_cwdl2) = i_str_lit
         nue_decomp_cascade(i_cwdl2) = 1.0_r8
      end if

      if (use_soil_matrixcn) call InitSoilTransfer()

      deallocate(params_inst%mimics_mge)
      deallocate(params_inst%mimics_vmod)
      deallocate(params_inst%mimics_vint)
      deallocate(params_inst%mimics_vslope)
      deallocate(params_inst%mimics_kmod)
      deallocate(params_inst%mimics_kint)
      deallocate(params_inst%mimics_kslope)
      deallocate(params_inst%mimics_p_scalar)
      deallocate(params_inst%mimics_desorp)
      deallocate(params_inst%mimics_fphys_r)
      deallocate(params_inst%mimics_fphys_k)
      deallocate(params_inst%mimics_initial_Cstocks)

    end associate

  end subroutine init_decompcascade_mimics

  !-----------------------------------------------------------------------
  subroutine decomp_rates_mimics(bounds, num_bgc_soilc, filter_bgc_soilc, &
       num_soilp, filter_soilp, clm_fates, &
       soilstate_inst, temperature_inst, cnveg_carbonflux_inst, &
       ch4_inst, soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
       idop)
    !
    ! !DESCRIPTION:
    ! Calculate rates and decomposition pathways for the MIMICS
    ! decomposition cascade model
    !
    ! !USES:
    use clm_time_manager , only : get_average_days_per_year, get_step_size
    use clm_varcon       , only : secspday, secsphr, tfrz
    use clm_varcon       , only : g_to_mg, cm3_to_m3
    use subgridAveMod    , only : p2c
    use PatchType        , only : patch
    use pftconMod        , only : pftname
    use TillageMod       , only : get_do_tillage
    use TillageMod       , only : get_apply_tillage_multipliers
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds          
    integer                              , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                              , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                              , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
    integer                              , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
    type(soilstate_type)                 , intent(in)    :: soilstate_inst
    type(temperature_type)               , intent(in)    :: temperature_inst
    type(cnveg_carbonflux_type)          , intent(in)    :: cnveg_carbonflux_inst
    type(ch4_type)                       , intent(in)    :: ch4_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type), intent(in)    :: soilbiogeochem_carbonstate_inst
    type(hlm_fates_interface_type)       , intent(inout) :: clm_fates
    integer, optional                    , intent(in)    :: idop(:) ! patch day of planting
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: eps = 1.e-6_r8
    real(r8):: frw(bounds%begc:bounds%endc) ! rooting fraction weight
    real(r8), allocatable:: fr(:,:)         ! column-level rooting fraction by soil depth
    real(r8):: psi                          ! temporary soilpsi for water scalar
    real(r8):: k_frag                       ! fragmentation rate constant CWD (1/sec)
    real(r8):: fmet
    real(r8):: favl
    real(r8):: fchem_m1
    real(r8):: fchem_m2
    real(r8):: desorption
    real(r8):: vmax_l1_m1  !
    real(r8):: vmax_l2_m1  !
    real(r8):: vmax_s1_m1  !
    real(r8):: vmax_l1_m2  !
    real(r8):: vmax_l2_m2  !
    real(r8):: vmax_s1_m2  !
    real(r8):: km_l1_m1  !
    real(r8):: km_l2_m1  !
    real(r8):: km_s1_m1  !
    real(r8):: km_l1_m2  !
    real(r8):: km_l2_m2  !
    real(r8):: km_s1_m2  !
    real(r8):: tau_m1  !
    real(r8):: tau_m2  !
    real(r8):: tau_mod
    real(r8):: m1_conc  !
    real(r8):: m2_conc  !
    real(r8):: term_1  !
    real(r8):: term_2  !
    real(r8):: t_soi_degC
!   real(r8):: decomp_depth_efolding        ! (meters) e-folding depth for reduction in decomposition [
    integer :: p, fp, c, fc, j, k, l, s  ! indices
    integer :: pf  ! fates patch index
    integer :: nc  ! clump index
    real(r8):: dt                           ! decomposition time step
    real(r8):: days_per_year                ! days per year
    real(r8):: depth_scalar(bounds%begc:bounds%endc,1:nlevdecomp) 
    real(r8):: w_d_o_scalars  ! product of w_scalar * depth_scalar * o_scalar
    real(r8):: mino2lim                     !minimum anaerobic decomposition rate
    real(r8):: mimics_fmet_p1
    real(r8):: mimics_fmet_p2
    real(r8):: mimics_fmet_p3
    real(r8):: mimics_fmet_p4
    real(r8):: mimics_fchem_r_p1
    real(r8):: mimics_fchem_r_p2
    real(r8):: mimics_fchem_k_p1
    real(r8):: mimics_fchem_k_p2
    real(r8):: mimics_tau_mod_min
    real(r8):: mimics_tau_mod_max
    real(r8):: mimics_tau_mod_factor
    real(r8):: mimics_tau_r_p1
    real(r8):: mimics_tau_r_p2
    real(r8):: mimics_tau_k_p1
    real(r8):: mimics_tau_k_p2
    real(r8):: mimics_ko_r
    real(r8):: mimics_ko_k
    real(r8):: mimics_densdep
    real(r8):: mimics_desorpQ10
    real(r8):: mimics_t_soi_ref
    real(r8):: mimics_cn_mod_num
    real(r8):: mimics_cn_r
    real(r8):: mimics_cn_k

    real(r8):: spinup_geogterm_l1(bounds%begc:bounds%endc) ! geographically-varying spinup term for l1
    real(r8):: spinup_geogterm_l2(bounds%begc:bounds%endc)  ! geographically-varying spinup term for l2
    real(r8):: spinup_geogterm_cwd(bounds%begc:bounds%endc) ! geographically-varying spinup term for cwd
    real(r8):: spinup_geogterm_s1(bounds%begc:bounds%endc) ! geographically-varying spinup term for s1
    real(r8):: spinup_geogterm_s2(bounds%begc:bounds%endc) ! geographically-varying spinup term for s2
    real(r8):: spinup_geogterm_s3(bounds%begc:bounds%endc) ! geographically-varying spinup term for s3
    real(r8):: spinup_geogterm_m1(bounds%begc:bounds%endc)  ! geographically-varying spinup term for m1
    real(r8):: spinup_geogterm_m2(bounds%begc:bounds%endc)  ! geographically-varying spinup term for m2
    real(r8):: annsum_npp_col_local(bounds%begc:bounds%endc)  ! local annual sum of NPP at the column level
    real(r8):: annsum_npp(bounds%begp:bounds%endp)  ! local annual sum of NPP at the patch level
    real(r8):: annsum_npp_col_scalar  ! annual sum of NPP, scalar in column-level loop

    !-----------------------------------------------------------------------

    associate(                                                           &

         rf_cwdl2       => CNParamsShareInst%rf_cwdl2                  , & ! Input:  [real(r8)         ]  respiration fraction in CWD to litter2 transition (frac)
         minpsi         => CNParamsShareInst%minpsi                    , & ! Input:  [real(r8)         ]  minimum soil suction (mm)
         maxpsi         => CNParamsShareInst%maxpsi                    , & ! Input:  [real(r8)         ]  maximum soil suction (mm)
         soilpsi        => soilstate_inst%soilpsi_col                  , & ! Input:  [real(r8) (:,:)   ]  soil water potential in each soil layer (MPa)          
         t_soisno       => temperature_inst%t_soisno_col               , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       
         o2stress_sat   => ch4_inst%o2stress_sat_col                   , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         o2stress_unsat => ch4_inst%o2stress_unsat_col                 , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         finundated     => ch4_inst%finundated_col                     , & ! Input:  [real(r8) (:)     ]  fractional inundated area                                
         decomp_cpools_vr => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col , &  ! Input: [real(r8) (:,:,:) ] (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) C pools
         pathfrac_decomp_cascade => soilbiogeochem_carbonflux_inst%pathfrac_decomp_cascade_col                                                                    , &  ! Output: [real(r8) (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)
         rf_decomp_cascade       => soilbiogeochem_carbonflux_inst%rf_decomp_cascade_col                                                                          , & ! Input:  [real(r8)          (:,:,:) ]  respired fraction in decomposition step (frac)
         w_scalar       => soilbiogeochem_carbonflux_inst%w_scalar_col , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
         o_scalar       => soilbiogeochem_carbonflux_inst%o_scalar_col , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
         cn_col         => soilbiogeochem_carbonflux_inst%cn_col       , & ! Output: [real(r8) (:,:)   ]  C:N ratio
         ligninNratioAvg => soilbiogeochem_carbonflux_inst%litr_lig_c_to_n_col, &  ! Input: [real(r8) (:) ] C:N ratio of litter lignin
         decomp_k       => soilbiogeochem_carbonflux_inst%decomp_k_col , & ! Output: [real(r8) (:,:,:) ]  rate for decomposition (1./sec)
         Ksoil          => soilbiogeochem_carbonflux_inst%Ksoil        , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         spinup_factor  => decomp_cascade_con%spinup_factor              & ! Input:  [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
         )

      if (get_do_tillage() .and. .not. present(idop)) then
         call endrun("Do not enable tillage without providing idop to decomp_rate_constants_mimics().")
      end if

      mino2lim = CNParamsShareInst%mino2lim

      days_per_year = get_average_days_per_year()
      dt = real( get_step_size(), r8 )

!     ! Set "decomp_depth_efolding" parameter
!     decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding

      ! translate to per-second time constant
      k_frag = 1._r8 / (secspday * days_per_year * CNParamsShareInst%tau_cwd)

     ! calc ref rate
      if ( spinup_state >= 1 ) then
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            !
            if ( abs(spinup_factor(i_met_lit) - 1._r8) .gt. eps) then
               spinup_geogterm_l1(c) = spinup_factor(i_met_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_str_lit) - 1._r8) .gt. eps) then
               spinup_geogterm_l2(c) = spinup_factor(i_str_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l2(c) = 1._r8
            endif
            !
            if ( .not. use_fates ) then
               if ( abs(spinup_factor(i_cwd) - 1._r8) .gt. eps) then
                  spinup_geogterm_cwd(c) = spinup_factor(i_cwd) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
               else
                  spinup_geogterm_cwd(c) = 1._r8
               endif
            endif
            !
            if ( abs(spinup_factor(i_avl_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s1(c) = spinup_factor(i_avl_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_chem_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s2(c) = spinup_factor(i_chem_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s2(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_phys_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s3(c) = spinup_factor(i_phys_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s3(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_cop_mic) - 1._r8) .gt. eps) then
               spinup_geogterm_m1(c) = spinup_factor(i_cop_mic) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_m1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_oli_mic) - 1._r8) .gt. eps) then
               spinup_geogterm_m2(c) = spinup_factor(i_oli_mic) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_m2(c) = 1._r8
            endif
            !
         end do
      else
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            spinup_geogterm_l1(c) = 1._r8
            spinup_geogterm_l2(c) = 1._r8
            spinup_geogterm_cwd(c) = 1._r8
            spinup_geogterm_s1(c) = 1._r8
            spinup_geogterm_s2(c) = 1._r8
            spinup_geogterm_s3(c) = 1._r8
            spinup_geogterm_m1(c) = 1._r8
            spinup_geogterm_m2(c) = 1._r8
         end do
      endif

      !--- time dependent coefficients-----!
      if ( nlevdecomp .eq. 1 ) then

         ! calculate function to weight the temperature and water potential scalars
         ! for decomposition control.  


         ! the following normalizes values in fr so that they
         ! sum to 1.0 across top nlevdecomp levels on a column
         frw(bounds%begc:bounds%endc) = 0._r8
         allocate(fr(bounds%begc:bounds%endc,nlev_soildecomp_standard))
         do j=1,nlev_soildecomp_standard
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               frw(c) = frw(c) + col%dz(c,j)
            end do
         end do
         do j = 1,nlev_soildecomp_standard
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (frw(c) /= 0._r8) then
                  fr(c,j) = col%dz(c,j) / frw(c)
               else
                  fr(c,j) = 0._r8
               end if
            end do
         end do

         ! calculate the rate constant scalar for soil water content.
         ! Uses the log relationship with water potential given in
         ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
         ! a comparison of models. Ecology, 68(5):1190-1200.
         ! and supported by data in
         ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
         ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

         do j = 1,nlev_soildecomp_standard
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (j==1) w_scalar(c,:) = 0._r8
               psi = min(soilpsi(c,j),maxpsi)
               ! decomp only if soilpsi is higher than minpsi
               if (psi > minpsi) then
                  w_scalar(c,1) = w_scalar(c,1) + (log(minpsi/psi)/log(minpsi/maxpsi))*fr(c,j)
               end if
            end do
         end do

         ! Calculate ANOXIA
         ! anoxia = .true. when (use_lch4)

         if (anoxia) then

            do j = 1,nlev_soildecomp_standard
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)

                  if (j==1) o_scalar(c,:) = 0._r8

                  o_scalar(c,1) = o_scalar(c,1) + fr(c,j) * max(o2stress_unsat(c,j), mino2lim)
               end do
            end do
         else
            o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
         end if

         deallocate(fr)

      else

         ! calculate the rate constant scalar for soil water content.
         ! Uses the log relationship with water potential given in
         ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
         ! a comparison of models. Ecology, 68(5):1190-1200.
         ! and supported by data in
         ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
         ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

         do j = 1,nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               psi = min(soilpsi(c,j),maxpsi)
               ! decomp only if soilpsi is higher than minpsi
               if (psi > minpsi) then
                  w_scalar(c,j) = (log(minpsi/psi)/log(minpsi/maxpsi))
               else
                  w_scalar(c,j) = 0._r8
               end if
            end do
         end do

         ! Calculate ANOXIA
         ! anoxia = .true. when (use_lch4)

         if (anoxia) then
            do j = 1,nlevdecomp
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)

                  o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim)
               end do
            end do
         else
            o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
         end if

      end if

      ! Term that reduces decomposition rate at depth
      ! Placeholder. For now depth_scalar = 1.
      do j = 1, nlevdecomp
         do fc = 1, num_bgc_soilc
            c = filter_bgc_soilc(fc)
            ! Using fixed e-folding depth as in
            ! SoilBiogeochemDecompCascadeBGCMod.F90
!           depth_scalar(c,j) = exp(-zsoi(j) / decomp_depth_efolding)
            depth_scalar(c,j) = 1.0_r8
         end do
      end do

      ! TODO @ekluzek suggested possibly making the Left Hand Sides into arrays
      ! and I wonder in that case whether to skip these assignments altogether
      ! and use the Right Hand Sides directly
      mimics_fmet_p1 = params_inst%mimics_fmet(1)
      mimics_fmet_p2 = params_inst%mimics_fmet(2)
      mimics_fmet_p3 = params_inst%mimics_fmet(3)
      mimics_fmet_p4 = params_inst%mimics_fmet(4)
      mimics_fchem_r_p1 = params_inst%mimics_fchem_r(1)
      mimics_fchem_r_p2 = params_inst%mimics_fchem_r(2)
      mimics_fchem_k_p1 = params_inst%mimics_fchem_k(1)
      mimics_fchem_k_p2 = params_inst%mimics_fchem_k(2)
      mimics_tau_mod_min = params_inst%mimics_tau_mod_min
      mimics_tau_mod_max = params_inst%mimics_tau_mod_max
      mimics_tau_mod_factor = params_inst%mimics_tau_mod_factor
      mimics_tau_r_p1 = params_inst%mimics_tau_r(1)
      mimics_tau_r_p2 = params_inst%mimics_tau_r(2)
      mimics_tau_k_p1 = params_inst%mimics_tau_k(1)
      mimics_tau_k_p2 = params_inst%mimics_tau_k(2)
      mimics_ko_r = params_inst%mimics_ko_r
      mimics_ko_k = params_inst%mimics_ko_k
      mimics_densdep = params_inst%mimics_densdep
      mimics_desorpQ10 = params_inst%mimics_desorpQ10
      mimics_t_soi_ref = params_inst%mimics_t_soi_ref
      mimics_cn_mod_num = params_inst%mimics_cn_mod_num
      mimics_cn_r = params_inst%mimics_cn_r
      mimics_cn_k = params_inst%mimics_cn_k

      ! If FATES-MIMICS, then use FATES copy of annsum_npp.
      ! The FATES copy of annsum_npp is available when use_lch4 = .true., so
      ! we limit FATES-MIMICS to if (use_lch4).
      fates_if: if (use_fates) then
         lch4_if: if (use_lch4) then

            ! Loop over p to get FATES copy of annsum_npp
            nc = bounds%clump_index
            do fp = 1, num_soilp

               p = filter_soilp(fp)
               c = patch%column(p)

               pf = p - col%patchi(c)
               s  = clm_fates%f2hmap(nc)%hsites(c)
               annsum_npp(p) = clm_fates%fates(nc)%bc_out(s)%annsum_npp_pa(pf)

               ! Initialize local column-level annsum_npp before averaging
               annsum_npp_col_local(c) = 0._r8

            end do  ! p loop

            ! Calculate the column-level average
            call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
                 annsum_npp(bounds%begp:bounds%endp), &
                 annsum_npp_col_local(bounds%begc:bounds%endc))
         else
            call endrun(msg='ERROR: soil_decomp_method = MIMICSWieder2015 '// &
              'will work with use_fates = .true. only if use_lch4 = .true. '// &
              errMsg(sourcefile, __LINE__))
         end if lch4_if
      end if fates_if

      ! calculate rates for all litter and som pools
      do fc = 1,num_bgc_soilc
         c = filter_bgc_soilc(fc)

         if (use_fates) then
            annsum_npp_col_scalar = max(0._r8, annsum_npp_col_local(c))
         else
            annsum_npp_col_scalar = max(0._r8, cnveg_carbonflux_inst%annsum_npp_col(c))
         end if

         ! Time-dependent params from Wieder et al. 2015 & testbed code

         ! Set limits on Lignin:N to keep fmet > 0.2
         ! Necessary for litter quality in boreal forests with high cwd flux
         ! TODO Check for high-freq variations in ligninNratioAvg. To avoid,
         !      replace pool_to_litter terms with ann or other long term mean
         !      in CNVegCarbonFluxType.
         fmet = mimics_fmet_p1 * (mimics_fmet_p2 - mimics_fmet_p3 * &
            min(mimics_fmet_p4, ligninNratioAvg(c)))
         tau_mod = min(mimics_tau_mod_max, max(mimics_tau_mod_min, &
            sqrt(mimics_tau_mod_factor * annsum_npp_col_scalar)))

         ! tau_m1 is tauR and tau_m2 is tauK in Wieder et al. 2015
         ! tau ends up in units of per hour but is expected
         ! in units of per second, so convert here; alternatively
         ! place the conversion once in w_d_o_scalars
         tau_m1 = mimics_tau_r_p1 * exp(mimics_tau_r_p2 * fmet) * tau_mod / &
                  secsphr
         tau_m2 = mimics_tau_k_p1 * exp(mimics_tau_k_p2 * fmet) * tau_mod / &
                  secsphr

         ! These two get used in SoilBiogeochemPotentialMod.F90
         ! cn(c,i_cop_mic), cn(c,i_oli_mic) are CN_r, CN_k in the testbed code
         ! cn_r, cn_k are CNr, CNk in the testbed code
         ! https://github.com/wwieder/biogeochem_testbed_1.1/blob/Testbed_CN/SOURCE_CODE/mimics_cycle_CN.f90#L1613
         cn_col(c,i_cop_mic) = mimics_cn_r * sqrt(mimics_cn_mod_num / fmet)
         cn_col(c,i_oli_mic) = mimics_cn_k * sqrt(mimics_cn_mod_num / fmet)

         ! Used in the update of certain pathfrac terms that vary with time
         ! in the next loop
         fchem_m1 = min(1._r8, max(0._r8, mimics_fchem_r_p1 * &
                    exp(mimics_fchem_r_p2 * fmet)))
         fchem_m2 = min(1._r8, max(0._r8, mimics_fchem_k_p1 * &
                    exp(mimics_fchem_k_p2 * fmet)))

         do j = 1,nlevdecomp
            ! vmax ends up in units of per hour but is expected
            ! in units of per second, so convert here; alternatively
            ! place the conversion once in w_d_o_scalars
            ! Table B1 Wieder et al. 2015 & MIMICS params file give diff
            ! kslope. I used the params file value(s).
            t_soi_degC = t_soisno(c,j) - tfrz

            vmax_l1_m1 = exp(vslope_l1_m1 * t_soi_degC + vint_l1_m1) * &
                         vmod_l1_m1 / secsphr
            vmax_l1_m2 = exp(vslope_l1_m2 * t_soi_degC + vint_l1_m2) * &
                         vmod_l1_m2 / secsphr
            vmax_l2_m1 = exp(vslope_l2_m1 * t_soi_degC + vint_l2_m1) * &
                         vmod_l2_m1 / secsphr
            vmax_l2_m2 = exp(vslope_l2_m2 * t_soi_degC + vint_l2_m2) * &
                         vmod_l2_m2 / secsphr
            vmax_s1_m1 = exp(vslope_s1_m1 * t_soi_degC + vint_s1_m1) * &
                         vmod_s1_m1 / secsphr
            vmax_s1_m2 = exp(vslope_s1_m2 * t_soi_degC + vint_s1_m2) * &
                         vmod_s1_m2 / secsphr

            km_l1_m1 = exp(kslope_l1_m1 * t_soi_degC + kint_l1_m1) * &
                       kmod_l1_m1
            km_l1_m2 = exp(kslope_l1_m2 * t_soi_degC + kint_l1_m2) * &
                       kmod_l1_m2
            km_l2_m1 = exp(kslope_l2_m1 * t_soi_degC + kint_l2_m1) * &
                       kmod_l2_m1
            km_l2_m2 = exp(kslope_l2_m2 * t_soi_degC + kint_l2_m2) * &
                       kmod_l2_m2
            km_s1_m1 = exp(kslope_s1_m1 * t_soi_degC + kint_s1_m1) * &
                       kmod_s1_m1 * p_scalar(c,j)
            km_s1_m2 = exp(kslope_s1_m2 * t_soi_degC + kint_s1_m2) * &
                       kmod_s1_m2 * p_scalar(c,j)

            ! Desorption a function of soil temperature and
            ! Q10 = 1.1 w/ reference temperature of 25C.
            ! Expected in units of per second, so convert; alternatively
            ! place the conversion once in w_d_o_scalars
            desorption = (desorp(c,j) / secsphr) * mimics_desorpQ10 * &
                         exp((t_soi_degC - mimics_t_soi_ref) / 10.0_r8)

            ! Microbial concentration with necessary unit conversions
            ! mgC/cm3 = gC/m3 * (1e3 mg/g) / (1e6 cm3/m3)
            m1_conc = (decomp_cpools_vr(c,j,i_cop_mic) / col%dz(c,j)) * &
                      g_to_mg * cm3_to_m3
            m2_conc = (decomp_cpools_vr(c,j,i_oli_mic) / col%dz(c,j)) * &
                      g_to_mg * cm3_to_m3

            ! Product of w_scalar * depth_scalar * o_scalar
            w_d_o_scalars = w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)

            ! decomp_k used in SoilBiogeochemPotentialMod.F90
            ! also updating pathfrac terms that vary with time
            term_1 = vmax_l1_m1 * m1_conc / (km_l1_m1 + m1_conc)
            term_2 = vmax_l1_m2 * m2_conc / (km_l1_m2 + m2_conc)
            decomp_k(c,j,i_met_lit) = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2 /= 0._r8) then
               pathfrac_decomp_cascade(c,j,i_l1m1) = term_1 / (term_1 + term_2)
               pathfrac_decomp_cascade(c,j,i_l1m2) = term_2 / (term_1 + term_2)
            else
               pathfrac_decomp_cascade(c,j,i_l1m1) = 0._r8
               pathfrac_decomp_cascade(c,j,i_l1m2) = 0._r8
            end if

            term_1 = vmax_l2_m1 * m1_conc / (km_l2_m1 + m1_conc)
            term_2 = vmax_l2_m2 * m2_conc / (km_l2_m2 + m2_conc)
            decomp_k(c,j,i_str_lit) = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2 /= 0._r8) then
               pathfrac_decomp_cascade(c,j,i_l2m1) = term_1 / (term_1 + term_2)
               pathfrac_decomp_cascade(c,j,i_l2m2) = term_2 / (term_1 + term_2)
            else
               pathfrac_decomp_cascade(c,j,i_l2m1) = 0._r8
               pathfrac_decomp_cascade(c,j,i_l2m2) = 0._r8
            end if

            term_1 = vmax_s1_m1 * m1_conc / (km_s1_m1 + m1_conc)
            term_2 = vmax_s1_m2 * m2_conc / (km_s1_m2 + m2_conc)
            decomp_k(c,j,i_avl_som) = (term_1 + term_2) * w_d_o_scalars
            if (term_1 + term_2 /= 0._r8) then
               pathfrac_decomp_cascade(c,j,i_s1m1) = term_1 / (term_1 + term_2)
               pathfrac_decomp_cascade(c,j,i_s1m2) = term_2 / (term_1 + term_2)
            else
               pathfrac_decomp_cascade(c,j,i_s1m1) = 0._r8
               pathfrac_decomp_cascade(c,j,i_s1m2) = 0._r8
            end if

            decomp_k(c,j,i_phys_som) = desorption * depth_scalar(c,j)

            term_1 = vmax_l2_m1 * m1_conc / (mimics_ko_r * km_l2_m1 + m1_conc)
            term_2 = vmax_l2_m2 * m2_conc / (mimics_ko_k * km_l2_m2 + m2_conc)
            ! The right hand side is OXIDAT in the testbed (line 1145)
            decomp_k(c,j,i_chem_som) = (term_1 + term_2) * w_d_o_scalars

            ! Currently, mimics_densdep = 1 so as to have no effect
            decomp_k(c,j,i_cop_mic) = tau_m1 * m1_conc**(mimics_densdep) 
 
            favl = min(1.0_r8, max(0.0_r8, 1.0_r8 - fphys_m1(c,j) - fchem_m1))
            pathfrac_decomp_cascade(c,j,i_m1s1) = favl
            pathfrac_decomp_cascade(c,j,i_m1s2) = fchem_m1

            decomp_k(c,j,i_oli_mic) = tau_m2 * m2_conc**(mimics_densdep)

            favl = min(1.0_r8, max(0.0_r8, 1.0_r8 - fphys_m2(c,j) - fchem_m2))
            pathfrac_decomp_cascade(c,j,i_m2s1) = favl
            pathfrac_decomp_cascade(c,j,i_m2s2) = fchem_m2

            ! Same for cwd but only if fates not enabled; fates handles cwd on
            ! its own structure
            ! TODO This shows how BGC applies the spinup coefficients
            if (.not. use_fates) then
               decomp_k(c,j,i_cwd) = k_frag * w_d_o_scalars  ! * spinup_geogterm_cwd(c)
            end if

            ! Tillage
            if (get_do_tillage()) then
               call get_apply_tillage_multipliers(idop, c, j, decomp_k(c,j,:))
            end if

! Above into soil matrix
            if(use_soil_matrixcn)then
               Ksoil%DM(c,j+nlevdecomp*(i_met_lit-1)) = decomp_k(c,j,i_met_lit) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_str_lit-1)) = decomp_k(c,j,i_str_lit) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_avl_som-1)) = decomp_k(c,j,i_avl_som) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_phys_som-1)) = decomp_k(c,j,i_phys_som) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_chem_som-1)) = decomp_k(c,j,i_chem_som) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_cop_mic-1)) = decomp_k(c,j,i_cop_mic) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_oli_mic-1)) = decomp_k(c,j,i_oli_mic) * dt
               ! same for cwd but only if fates is not enabled; fates handles
               ! CWD
               ! on its own structure
               if (.not. use_fates) then
                  Ksoil%DM(c,j+nlevdecomp*(i_cwd-1))   = decomp_k(c,j,i_cwd) * dt
               end if
            end if !use_soil_matrixcn
         end do
      end do

      ! pathfrac terms not calculated in the previous loop
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = 1.0_r8
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = 1.0_r8
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m1s3) = fphys_m1(bounds%begc:bounds%endc,1:nlevdecomp)
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m2s3) = fphys_m2(bounds%begc:bounds%endc,1:nlevdecomp)
      if (.not. use_fates) then
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = 1.0_r8
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
      end if
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m1) = rf_l1m1
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m2) = rf_l1m2
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m1) = rf_l2m1
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m2) = rf_l2m2
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1m1) = rf_s1m1
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1m2) = rf_s1m2
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = 0.0_r8
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = 0.0_r8
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m1s1) = 0.0_r8
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m1s2) = 0.0_r8
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m1s3) = 0.0_r8
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m2s1) = 0.0_r8
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m2s2) = 0.0_r8
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m2s3) = 0.0_r8

    end associate

 end subroutine decomp_rates_mimics

end module SoilBiogeochemDecompCascadeMIMICSMod
