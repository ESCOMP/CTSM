module SoilBiogeochemDecompCascadeBGCMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Sets the coeffiecients used in the decomposition cascade submodel.  
  ! This uses the CENTURY/BGC parameters
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_const_mod                      , only : SHR_CONST_TKFRZ
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : nlevdecomp, ndecomp_pools_max
  use clm_varpar                         , only : i_litr_min, i_litr_max, i_met_lit, i_cwd, i_cwdl2
  use clm_varctl                         , only : iulog, spinup_state, anoxia, use_lch4, use_fates
  use clm_varcon                         , only : zsoi
  use decompMod                          , only : bounds_type
  use spmdMod                            , only : masterproc
  use abortutils                         , only : endrun
  use CNSharedParamsMod                  , only : CNParamsShareInst, nlev_soildecomp_standard 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, InitSoilTransfer, use_soil_matrixcn
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilStateType                      , only : soilstate_type
  use TemperatureType                    , only : temperature_type 
  use ch4Mod                             , only : ch4_type
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term

  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams                      ! Read in parameters from params file
  public :: init_decompcascade_bgc          ! Initialization
  public :: decomp_rate_constants_bgc       ! Figure out decomposition rates
  !
  ! !PUBLIC DATA MEMBERS 
  logical , public :: normalize_q10_to_century_tfunc = .true.! do we normalize the century decomp. rates so that they match the CLM Q10 at a given tep?
  logical , public :: use_century_tfunc = .false.
  real(r8), public :: normalization_tref = 15._r8            ! reference temperature for normalizaion (degrees C)
  !
  ! !PRIVATE DATA MEMBERS 
  integer, private :: i_pas_som  ! index of passive (aka protected) Soil Organic Matter (SOM)
  integer, private :: i_slo_som  ! index of slow (aka recalcitrant) SOM
  integer, private :: i_act_som  ! index of active (aka available) SOM
  integer, private :: i_cel_lit  ! index of cellulose litter pool
  integer, private :: i_lig_lit  ! index of lignin litter pool

  real(r8), private :: cwd_fcel
  real(r8), private :: rf_l1s1
  real(r8), private :: rf_l2s1
  real(r8), private :: rf_l3s2
  real(r8), private :: rf_s2s1
  real(r8), private :: rf_s2s3
  real(r8), private :: rf_s3s1
  real(r8), private :: rf_cwdl3
  real(r8), private, allocatable :: rf_s1s2(:,:)
  real(r8), private, allocatable :: rf_s1s3(:,:)
  real(r8), private, allocatable :: f_s1s2(:,:)
  real(r8), private, allocatable :: f_s1s3(:,:)
  real(r8), private :: f_s2s1
  real(r8), private :: f_s2s3

  integer, private :: i_l1s1
  integer, private :: i_l2s1
  integer, private :: i_l3s2
  integer, private :: i_s1s2
  integer, private :: i_s1s3
  integer, private :: i_s2s1
  integer, private :: i_s2s3
  integer, private :: i_s3s1
  integer, private :: i_cwdl3

  type, private :: params_type
     real(r8):: cn_s1_bgc     !C:N for SOM 1
     real(r8):: cn_s2_bgc     !C:N for SOM 2
     real(r8):: cn_s3_bgc     !C:N for SOM 3

     real(r8):: rf_l1s1_bgc   !respiration fraction litter 1 -> SOM 1
     real(r8):: rf_l2s1_bgc
     real(r8):: rf_l3s2_bgc

     real(r8):: rf_s2s1_bgc    
     real(r8):: rf_s2s3_bgc    
     real(r8):: rf_s3s1_bgc    

     real(r8):: rf_cwdl3_bgc

     real(r8):: tau_l1_bgc    ! 1/turnover time of  litter 1 from Century (l/18.5) (1/yr)
     real(r8):: tau_l2_l3_bgc ! 1/turnover time of  litter 2 and litter 3 from Century (1/4.9) (1/yr)
     real(r8):: tau_s1_bgc    ! 1/turnover time of  SOM 1 from Century (1/7.3) (1/yr)
     real(r8):: tau_s2_bgc    ! 1/turnover time of  SOM 2 from Century (1/0.2) (1/yr)
     real(r8):: tau_s3_bgc    ! 1/turnover time of  SOM 3 from Century (1/0.0045) (1/yr)

     real(r8) :: cwd_fcel_bgc !cellulose fraction for CWD

     real(r8), allocatable :: bgc_initial_Cstocks(:)  ! Initial Carbon stocks for a cold-start
     real(r8) :: bgc_initial_Cstocks_depth  ! Soil depth for initial Carbon stocks for a cold-start
     
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
    character(len=32)  :: subname = 'CNDecompBgcParamsType'
    character(len=100) :: errCode = 'Error reading in CN const file '
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! Read off of netcdf file
    tString='bgc_tau_l1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_l1_bgc=tempr

    tString='bgc_tau_l2_l3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_l2_l3_bgc=tempr

    tString='bgc_tau_s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s1_bgc=tempr

    tString='bgc_tau_s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s2_bgc=tempr

    tString='bgc_tau_s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s3_bgc=tempr

    tString='bgc_cn_s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s1_bgc=tempr

    tString='bgc_cn_s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s2_bgc=tempr

    tString='bgc_cn_s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s3_bgc=tempr

    tString='bgc_rf_l1s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l1s1_bgc=tempr

    tString='bgc_rf_l2s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l2s1_bgc=tempr

    tString='bgc_rf_l3s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l3s2_bgc=tempr   

    tString='bgc_rf_s2s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s2s1_bgc=tempr

    tString='bgc_rf_s2s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s2s3_bgc=tempr

    tString='bgc_rf_s3s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s3s1_bgc=tempr

    tString='bgc_rf_cwdl3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl3_bgc=tempr

    tString='bgc_cwd_fcel'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cwd_fcel_bgc=tempr

    allocate(params_inst%bgc_initial_Cstocks(ndecomp_pools_max))
    tString='bgc_initial_Cstocks'
    call ncd_io(trim(tString), params_inst%bgc_initial_Cstocks(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    tString='bgc_initial_Cstocks_depth'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%bgc_initial_Cstocks_depth=tempr

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_bgc(bounds, soilbiogeochem_state_inst, soilstate_inst )
    !
    ! !DESCRIPTION:
    !  initialize rate constants and decomposition pathways following the decomposition cascade of the BGC model.
    !  written by C. Koven 
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds  
    type(soilbiogeochem_state_type) , intent(inout) :: soilbiogeochem_state_inst
    type(soilstate_type)            , intent(in)    :: soilstate_inst
    !
    ! !LOCAL VARIABLES
    !-- properties of each decomposing pool
    real(r8) :: cn_s1
    real(r8) :: cn_s2
    real(r8) :: cn_s3
    real(r8):: speedup_fac                  ! acceleration factor, higher when vertsoilc = .true.

    integer  :: c, j    ! indices
    real(r8) :: t       ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                                     &
         cellsand                       => soilstate_inst%cellsand_col                           , & ! Input:  [real(r8)          (:,:)   ]  column 3D sand                                         
         
         cascade_donor_pool             => decomp_cascade_con%cascade_donor_pool                 , & ! Output: [integer           (:)     ]  which pool is C taken from for a given decomposition step 
         cascade_receiver_pool          => decomp_cascade_con%cascade_receiver_pool              , & ! Output: [integer           (:)     ]  which pool is C added to for a given decomposition step   
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools     , & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:N ratio                          
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

      allocate(rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp))

      !------- time-constant coefficients ---------- !
      ! set soil organic matter compartment C:N ratios
      cn_s1 = params_inst%cn_s1_bgc
      cn_s2 = params_inst%cn_s2_bgc
      cn_s3 = params_inst%cn_s3_bgc

      ! set respiration fractions for fluxes between compartments
      rf_l1s1 = params_inst%rf_l1s1_bgc
      rf_l2s1 = params_inst%rf_l2s1_bgc
      rf_l3s2 = params_inst%rf_l3s2_bgc
      rf_s2s1 = params_inst%rf_s2s1_bgc
      rf_s2s3 = params_inst%rf_s2s3_bgc
      rf_s3s1 = params_inst%rf_s3s1_bgc

      rf_cwdl3 = params_inst%rf_cwdl3_bgc

      ! set the cellulose and lignin fractions for coarse woody debris
      cwd_fcel = params_inst%cwd_fcel_bgc

      ! set path fractions
      f_s2s1 = 0.42_r8/(0.45_r8)
      f_s2s3 = 0.03_r8/(0.45_r8)

      ! some of these are dependent on the soil texture properties
      do c = bounds%begc, bounds%endc
         do j = 1, nlevdecomp
            t = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - cellsand(c,j))
            f_s1s2(c,j) = 1._r8 - .004_r8 / (1._r8 - t)
            f_s1s3(c,j) = .004_r8 / (1._r8 - t)
            rf_s1s2(c,j) = t
            rf_s1s3(c,j) = t
         end do
      end do
      initial_stock_soildepth = params_inst%bgc_initial_Cstocks_depth

      !-------------------  list of pools and their attributes  ------------
      i_litr_min = 1
      i_met_lit = i_litr_min
      floating_cn_ratio_decomp_pools(i_met_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_met_lit) = 'litr1'
      decomp_cascade_con%decomp_pool_name_history(i_met_lit) = 'LIT_MET'
      decomp_cascade_con%decomp_pool_name_long(i_met_lit) = 'metabolic litter'
      decomp_cascade_con%decomp_pool_name_short(i_met_lit) = 'L1'
      is_litter(i_met_lit) = .true.
      is_soil(i_met_lit) = .false.
      is_cwd(i_met_lit) = .false.
      initial_cn_ratio(i_met_lit) = 90._r8
      initial_stock(i_met_lit) = params_inst%bgc_initial_Cstocks(i_met_lit)
      is_metabolic(i_met_lit) = .true.
      is_cellulose(i_met_lit) = .false.
      is_lignin(i_met_lit) = .false.

      i_cel_lit = i_met_lit + 1
      floating_cn_ratio_decomp_pools(i_cel_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_cel_lit) = 'litr2'
      decomp_cascade_con%decomp_pool_name_history(i_cel_lit) = 'LIT_CEL'
      decomp_cascade_con%decomp_pool_name_long(i_cel_lit) = 'cellulosic litter'
      decomp_cascade_con%decomp_pool_name_short(i_cel_lit) = 'L2'
      is_litter(i_cel_lit) = .true.
      is_soil(i_cel_lit) = .false.
      is_cwd(i_cel_lit) = .false.
      initial_cn_ratio(i_cel_lit) = 90._r8
      initial_stock(i_cel_lit) = params_inst%bgc_initial_Cstocks(i_cel_lit)
      is_metabolic(i_cel_lit) = .false.
      is_cellulose(i_cel_lit) = .true.
      is_lignin(i_cel_lit) = .false.

      i_lig_lit = i_cel_lit + 1
      floating_cn_ratio_decomp_pools(i_lig_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_lig_lit) = 'litr3'
      decomp_cascade_con%decomp_pool_name_history(i_lig_lit) = 'LIT_LIG'
      decomp_cascade_con%decomp_pool_name_long(i_lig_lit) = 'lignin litter'
      decomp_cascade_con%decomp_pool_name_short(i_lig_lit) = 'L3'
      is_litter(i_lig_lit) = .true.
      is_soil(i_lig_lit) = .false.
      is_cwd(i_lig_lit) = .false.
      initial_cn_ratio(i_lig_lit) = 90._r8
      initial_stock(i_lig_lit) = params_inst%bgc_initial_Cstocks(i_lig_lit)
      is_metabolic(i_lig_lit) = .false.
      is_cellulose(i_lig_lit) = .false.
      is_lignin(i_lig_lit) = .true.

      i_litr_max = i_lig_lit
      if (i_litr_min /= 1 .or. i_litr_max < 2 .or. i_litr_max > 3) then
         write(iulog,*) 'Expecting i_litr_min = 1 and i_litr_max = 2 or 3.'
         write(iulog,*) 'See pftconMod, SoilBiogeochemCarbonFluxType, and'
         write(iulog,*) 'clmfates_interfaceMod for ramifications of changing'
         write(iulog,*) 'this assumption.'
         call endrun(msg='ERROR: i_litr_min and/or i_litr_max out of range '// &
              errMsg(sourcefile, __LINE__))
      end if

      i_act_som = i_lig_lit + 1
      floating_cn_ratio_decomp_pools(i_act_som) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_act_som) = 'soil1'
      decomp_cascade_con%decomp_pool_name_history(i_act_som) = 'SOM_ACT'
      decomp_cascade_con%decomp_pool_name_long(i_act_som) = 'active soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_act_som) = 'S1'
      is_litter(i_act_som) = .false.
      is_soil(i_act_som) = .true.
      is_cwd(i_act_som) = .false.
      initial_cn_ratio(i_act_som) = cn_s1
      initial_stock(i_act_som) = params_inst%bgc_initial_Cstocks(i_act_som)
      is_metabolic(i_act_som) = .false.
      is_cellulose(i_act_som) = .false.
      is_lignin(i_act_som) = .false.

      i_slo_som = i_act_som + 1
      floating_cn_ratio_decomp_pools(i_slo_som) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_slo_som) = 'soil2'
      decomp_cascade_con%decomp_pool_name_history(i_slo_som) = 'SOM_SLO'
      decomp_cascade_con%decomp_pool_name_long(i_slo_som) = 'slow soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_slo_som) = 'S2'
      is_litter(i_slo_som) = .false.
      is_soil(i_slo_som) = .true.
      is_cwd(i_slo_som) = .false.
      initial_cn_ratio(i_slo_som) = cn_s2
      initial_stock(i_slo_som) = params_inst%bgc_initial_Cstocks(i_slo_som)
      is_metabolic(i_slo_som) = .false.
      is_cellulose(i_slo_som) = .false.
      is_lignin(i_slo_som) = .false.

      i_pas_som = i_slo_som + 1
      floating_cn_ratio_decomp_pools(i_pas_som) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_pas_som) = 'soil3'
      decomp_cascade_con%decomp_pool_name_history(i_pas_som) = 'SOM_PAS'
      decomp_cascade_con%decomp_pool_name_long(i_pas_som) = 'passive soil organic matter'
      decomp_cascade_con%decomp_pool_name_short(i_pas_som) = 'S3'
      is_litter(i_pas_som) = .false.
      is_soil(i_pas_som) = .true.
      is_cwd(i_pas_som) = .false.
      initial_cn_ratio(i_pas_som) = cn_s3
      initial_stock(i_pas_som) = params_inst%bgc_initial_Cstocks(i_pas_som)
      is_metabolic(i_pas_som) = .false.
      is_cellulose(i_pas_som) = .false.
      is_lignin(i_pas_som) = .false.

      if (.not. use_fates) then
         ! CWD
         i_cwd = i_pas_som + 1
         floating_cn_ratio_decomp_pools(i_cwd) = .true.
         decomp_cascade_con%decomp_pool_name_restart(i_cwd) = 'cwd'
         decomp_cascade_con%decomp_pool_name_history(i_cwd) = 'CWD'
         decomp_cascade_con%decomp_pool_name_long(i_cwd) = 'coarse woody debris'
         decomp_cascade_con%decomp_pool_name_short(i_cwd) = 'CWD'
         is_litter(i_cwd) = .false.
         is_soil(i_cwd) = .false.
         is_cwd(i_cwd) = .true.
         initial_cn_ratio(i_cwd) = 90._r8
         initial_stock(i_cwd) = params_inst%bgc_initial_Cstocks(i_cwd)
         is_metabolic(i_cwd) = .false.
         is_cellulose(i_cwd) = .false.
         is_lignin(i_cwd) = .false.
      endif

      speedup_fac = 1._r8

      !lit1
      spinup_factor(i_met_lit) = 1._r8
      !lit2,3
      spinup_factor(i_cel_lit) = 1._r8
      spinup_factor(i_lig_lit) = 1._r8
      !CWD
      if (.not. use_fates) then
         spinup_factor(i_cwd) = max(1._r8, (speedup_fac * CNParamsShareInst%tau_cwd / 2._r8 ))
      end if
      !som1
      spinup_factor(i_act_som) = 1._r8
      !som2,3
      spinup_factor(i_slo_som) = max(1._r8, (speedup_fac * params_inst%tau_s2_bgc))
      spinup_factor(i_pas_som) = max(1._r8, (speedup_fac * params_inst%tau_s3_bgc))

      if ( masterproc ) then
         write(iulog,*) 'Spinup_state ',spinup_state
         write(iulog,*) 'Spinup factors ',spinup_factor
      end if

      !----------------  list of transitions and their time-independent coefficients  ---------------!
      i_l1s1 = 1
      decomp_cascade_con%cascade_step_name(i_l1s1) = 'L1S1'
      cascade_donor_pool(i_l1s1) = i_met_lit
      cascade_receiver_pool(i_l1s1) = i_act_som

      i_l2s1 = 2
      decomp_cascade_con%cascade_step_name(i_l2s1) = 'L2S1'
      cascade_donor_pool(i_l2s1) = i_cel_lit
      cascade_receiver_pool(i_l2s1) = i_act_som

      i_l3s2 = 3
      decomp_cascade_con%cascade_step_name(i_l3s2) = 'L3S2'
      cascade_donor_pool(i_l3s2) = i_lig_lit
      cascade_receiver_pool(i_l3s2) = i_slo_som

      i_s1s2 = 4
      decomp_cascade_con%cascade_step_name(i_s1s2) = 'S1S2'
      cascade_donor_pool(i_s1s2) = i_act_som
      cascade_receiver_pool(i_s1s2) = i_slo_som

      i_s1s3 = 5
      decomp_cascade_con%cascade_step_name(i_s1s3) = 'S1S3'
      cascade_donor_pool(i_s1s3) = i_act_som
      cascade_receiver_pool(i_s1s3) = i_pas_som

      i_s2s1 = 6
      decomp_cascade_con%cascade_step_name(i_s2s1) = 'S2S1'
      cascade_donor_pool(i_s2s1) = i_slo_som
      cascade_receiver_pool(i_s2s1) = i_act_som

      i_s2s3 = 7 
      decomp_cascade_con%cascade_step_name(i_s2s3) = 'S2S3'
      cascade_donor_pool(i_s2s3) = i_slo_som
      cascade_receiver_pool(i_s2s3) = i_pas_som

      i_s3s1 = 8
      decomp_cascade_con%cascade_step_name(i_s3s1) = 'S3S1'
      cascade_donor_pool(i_s3s1) = i_pas_som
      cascade_receiver_pool(i_s3s1) = i_act_som

      if (.not. use_fates) then
         i_cwdl2 = 9
         decomp_cascade_con%cascade_step_name(i_cwdl2) = 'CWDL2'
         cascade_donor_pool(i_cwdl2) = i_cwd
         cascade_receiver_pool(i_cwdl2) = i_cel_lit
         
         i_cwdl3 = 10
         decomp_cascade_con%cascade_step_name(i_cwdl3) = 'CWDL3'
         cascade_donor_pool(i_cwdl3) = i_cwd
         cascade_receiver_pool(i_cwdl3) = i_lig_lit
      end if
 
      if(use_soil_matrixcn) call InitSoilTransfer()

      deallocate(params_inst%bgc_initial_Cstocks)

    end associate

  end subroutine init_decompcascade_bgc

  !-----------------------------------------------------------------------
  subroutine decomp_rate_constants_bgc(bounds, num_bgc_soilc, filter_bgc_soilc, &
       soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst, &
       idop)
    !
    ! !DESCRIPTION:
    !  calculate rate constants and decomposition pathways for the CENTURY decomposition cascade model
    !  written by C. Koven based on original CLM4 decomposition cascade
    !
    ! !USES:
    use clm_time_manager , only : get_average_days_per_year, get_step_size
    use shr_const_mod    , only : SHR_CONST_PI
    use clm_varcon       , only : secspday
    use TillageMod       , only : get_do_tillage
    use TillageMod       , only : get_apply_tillage_multipliers
    use landunit_varcon  , only : istcrop
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds          
    integer                              , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
    integer                              , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
    type(soilstate_type)                 , intent(in)    :: soilstate_inst
    type(temperature_type)               , intent(in)    :: temperature_inst
    type(ch4_type)                       , intent(in)    :: ch4_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    integer, optional                    , intent(in)    :: idop(:) ! patch day of planting
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: eps = 1.e-6_r8
    real(r8):: frw(bounds%begc:bounds%endc) ! rooting fraction weight
    real(r8), allocatable:: fr(:,:)         ! column-level rooting fraction by soil depth
    real(r8):: psi                          ! temporary soilpsi for water scalar
    real(r8):: rate_scalar                  ! combined rate scalar for decomp
    real(r8):: k_l1                         ! decomposition rate constant litter 1 (1/sec)
    real(r8):: k_l2_l3                      ! decomposition rate constant litter 2 and litter 3 (1/sec)
    real(r8):: k_s1                         ! decomposition rate constant SOM 1 (1/sec)
    real(r8):: k_s2                         ! decomposition rate constant SOM 2 (1/sec)
    real(r8):: k_s3                         ! decomposition rate constant SOM 3 (1/sec)
    real(r8):: k_frag                       ! fragmentation rate constant CWD (1/sec)
    real(r8):: Q10                          ! temperature dependence
    real(r8):: froz_q10                     ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
    real(r8):: decomp_depth_efolding        ! (meters) e-folding depth for reduction in decomposition [
    integer :: c, fc, j, k, l
    real(r8):: dt                           ! decomposition time step
    real(r8):: catanf                       ! hyperbolic temperature function from CENTURY
    real(r8):: catanf_30                    ! reference rate at 30C
    real(r8):: t1                           ! temperature argument
    real(r8):: normalization_factor         ! factor by which to offset the decomposition rates frm century to a q10 formulation
    real(r8):: days_per_year                ! days per year
    real(r8):: depth_scalar(bounds%begc:bounds%endc,1:nlevdecomp) 
    real(r8):: mino2lim                     !minimum anaerobic decomposition rate
    real(r8):: spinup_geogterm_l1(bounds%begc:bounds%endc) ! geographically-varying spinup term for l1
    real(r8):: spinup_geogterm_l23(bounds%begc:bounds%endc) ! geographically-varying spinup term for l2 and l3
    real(r8):: spinup_geogterm_cwd(bounds%begc:bounds%endc) ! geographically-varying spinup term for cwd
    real(r8):: spinup_geogterm_s1(bounds%begc:bounds%endc) ! geographically-varying spinup term for s1
    real(r8):: spinup_geogterm_s2(bounds%begc:bounds%endc) ! geographically-varying spinup term for s2
    real(r8):: spinup_geogterm_s3(bounds%begc:bounds%endc) ! geographically-varying spinup term for s3

    !-----------------------------------------------------------------------

    !----- CENTURY T response function
    catanf(t1) = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))

    associate(                                                           &
         cwd_flig       => CNParamsShareInst%cwd_flig                  , & ! Input:  [real(r8)         ]  lignin fraction of coarse woody debris (frac)
         rf_cwdl2       => CNParamsShareInst%rf_cwdl2                  , & ! Input:  [real(r8)         ]  respiration fraction in CWD to litter2 transition (frac)
         minpsi         => CNParamsShareInst%minpsi                    , & ! Input:  [real(r8)         ]  minimum soil suction (mm)
         maxpsi         => CNParamsShareInst%maxpsi                    , & ! Input:  [real(r8)         ]  maximum soil suction (mm)
         soilpsi        => soilstate_inst%soilpsi_col                  , & ! Input:  [real(r8) (:,:)   ]  soil water potential in each soil layer (MPa)          

         t_soisno       => temperature_inst%t_soisno_col               , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       

         o2stress_sat   => ch4_inst%o2stress_sat_col                   , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         o2stress_unsat => ch4_inst%o2stress_unsat_col                 , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         finundated     => ch4_inst%finundated_col                     , & ! Input:  [real(r8) (:)     ]  fractional inundated area                                
         rf_decomp_cascade       => soilbiogeochem_carbonflux_inst%rf_decomp_cascade_col                                                               , & ! Output: [real(r8) (:,:,:) ]  respired fraction in decomposition step (frac)
         pathfrac_decomp_cascade => soilbiogeochem_carbonflux_inst%pathfrac_decomp_cascade_col                                                         , & ! Output: [real(r8) (:,:,:) ]  what fraction of C passes from donor to receiver pool through a given transition (frac)
         t_scalar       => soilbiogeochem_carbonflux_inst%t_scalar_col , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp                     
         w_scalar       => soilbiogeochem_carbonflux_inst%w_scalar_col , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
         o_scalar       => soilbiogeochem_carbonflux_inst%o_scalar_col , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
         decomp_k       => soilbiogeochem_carbonflux_inst%decomp_k_col , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         Ksoil          => soilbiogeochem_carbonflux_inst%Ksoil        , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         spinup_factor  => decomp_cascade_con%spinup_factor            & ! Input:  [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
         )

      mino2lim = CNParamsShareInst%mino2lim

      if ( use_century_tfunc .and. normalize_q10_to_century_tfunc ) then
         call endrun(msg='ERROR: cannot have both use_century_tfunc and normalize_q10_to_century_tfunc set as true'//&
              errMsg(sourcefile, __LINE__))
      endif

      if (get_do_tillage() .and. .not. present(idop)) then
         call endrun("Do not enable tillage without providing idop to decomp_rate_constants_bgc().")
      end if

      days_per_year = get_average_days_per_year()
      dt = real( get_step_size(), r8 )

      ! set "Q10" parameter
      Q10 = CNParamsShareInst%Q10

      ! set "froz_q10" parameter
      froz_q10  = CNParamsShareInst%froz_q10 

      ! Set "decomp_depth_efolding" parameter
      decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding

      ! translate to per-second time constant
      k_l1 = 1._r8    / (secspday * days_per_year * params_inst%tau_l1_bgc)
      k_l2_l3 = 1._r8 / (secspday * days_per_year * params_inst%tau_l2_l3_bgc)
      k_s1 = 1._r8    / (secspday * days_per_year * params_inst%tau_s1_bgc)
      k_s2 = 1._r8    / (secspday * days_per_year * params_inst%tau_s2_bgc)
      k_s3 = 1._r8    / (secspday * days_per_year * params_inst%tau_s3_bgc)
      k_frag = 1._r8  / (secspday * days_per_year * CNParamsShareInst%tau_cwd)

     ! calc ref rate
      catanf_30 = catanf(30._r8)

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
            if ( abs(spinup_factor(i_cel_lit) - 1._r8) .gt. eps) then
               spinup_geogterm_l23(c) = spinup_factor(i_cel_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l23(c) = 1._r8
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
            if ( abs(spinup_factor(i_act_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s1(c) = spinup_factor(i_act_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_slo_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s2(c) = spinup_factor(i_slo_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s2(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_pas_som) - 1._r8) .gt. eps) then
               spinup_geogterm_s3(c) = spinup_factor(i_pas_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s3(c) = 1._r8
            endif
            !
         end do
      else
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            spinup_geogterm_l1(c) = 1._r8
            spinup_geogterm_l23(c) = 1._r8
            spinup_geogterm_cwd(c) = 1._r8
            spinup_geogterm_s1(c) = 1._r8
            spinup_geogterm_s2(c) = 1._r8
            spinup_geogterm_s3(c) = 1._r8
         end do
      endif

      !--- time dependent coefficients-----!
      if ( nlevdecomp .eq. 1 ) then

         ! calculate function to weight the temperature and water potential scalars
         ! for decomposition control.  


         ! the following normalizes values in fr so that they
         ! sum to 1.0 across top nlevdecomp levels on a column
         frw(bounds%begc:bounds%endc) = 0._r8
         nlev_soildecomp_standard=5
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

         if ( .not. use_century_tfunc ) then
            ! calculate rate constant scalar for soil temperature
            ! assuming that the base rate constants are assigned for non-moisture
            ! limiting conditions at 25 C. 

            do j = 1,nlev_soildecomp_standard
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if (j==1) t_scalar(c,:) = 0._r8
                  if (t_soisno(c,j) >= SHR_CONST_TKFRZ) then
                     t_scalar(c,1)=t_scalar(c,1) + &
                          (Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))*fr(c,j)
                  else
                     t_scalar(c,1)=t_scalar(c,1) + &
                          (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))*fr(c,j)
                  endif
               end do
            end do

         else
            ! original century uses an arctangent function to calculate the temperature dependence of decomposition
            do j = 1,nlev_soildecomp_standard
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if (j==1) t_scalar(c,:) = 0._r8

                  t_scalar(c,1)=t_scalar(c,1) +max(catanf(t_soisno(c,j)-SHR_CONST_TKFRZ)/catanf_30*fr(c,j),0.01_r8)
               end do
            end do

         endif

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

         if (use_lch4) then
            ! Calculate ANOXIA
            if (anoxia) then
               ! Check for anoxia w/o LCH4 now done in controlMod.

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
         else
            o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
         end if

         deallocate(fr)

      else

         if ( .not. use_century_tfunc ) then
            ! calculate rate constant scalar for soil temperature
            ! assuming that the base rate constants are assigned for non-moisture
            ! limiting conditions at 25 C. 
            ! Peter Thornton: 3/13/09
            ! Replaced the Lloyd and Taylor function with a Q10 formula, with Q10 = 1.5
            ! as part of the modifications made to improve the seasonal cycle of 
            ! atmospheric CO2 concentration in global simulations. This does not impact
            ! the base rates at 25 C, which are calibrated from microcosm studies.

            do j = 1, nlevdecomp
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if (t_soisno(c,j) >= SHR_CONST_TKFRZ) then
                     t_scalar(c,j)= (Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
                  else
                     t_scalar(c,j)= (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))
                  endif
               end do
            end do

         else

            do j = 1, nlevdecomp
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  t_scalar(c,j)= max(catanf(t_soisno(c,j)-SHR_CONST_TKFRZ)/catanf_30, 0.01_r8)
               end do
            end do

         endif

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

         if (use_lch4) then
            ! Calculate ANOXIA
            ! Check for anoxia w/o LCH4 now done in controlMod.

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
         else
            o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
         end if

      end if

      if ( normalize_q10_to_century_tfunc ) then
         ! scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
         normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10**((normalization_tref-25._r8)/10._r8))
         do j = 1, nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               t_scalar(c,j) = t_scalar(c,j) * normalization_factor
            end do
         end do
      endif

      ! add a term to reduce decomposition rate at depth
      ! for now used a fixed e-folding depth
      do j = 1, nlevdecomp
         do fc = 1, num_bgc_soilc
            c = filter_bgc_soilc(fc)
            depth_scalar(c,j) = exp(-zsoi(j) / decomp_depth_efolding)
         end do
      end do

      ! calculate rate constants for all litter and som pools
      do j = 1,nlevdecomp
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            decomp_k(c,j,i_met_lit) = k_l1    * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l1(c)
            decomp_k(c,j,i_cel_lit) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l23(c)
            decomp_k(c,j,i_lig_lit) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l23(c)
            decomp_k(c,j,i_act_som) = k_s1    * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s1(c)
            decomp_k(c,j,i_slo_som) = k_s2    * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s2(c)
            decomp_k(c,j,i_pas_som) = k_s3    * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s3(c)
            ! same for cwd but only if fates is not enabled; fates handles CWD
            ! on its own structure
            if (.not. use_fates) then
               decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_cwd(c)
            end if

            ! Tillage
            if (get_do_tillage()) then
               call get_apply_tillage_multipliers(idop, c, j, decomp_k(c,j,:))
            end if

            ! Above into soil matrix
            if(use_soil_matrixcn)then
               Ksoil%DM(c,j+nlevdecomp*(i_met_lit-1)) = decomp_k(c,j,i_met_lit) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_cel_lit-1)) = decomp_k(c,j,i_cel_lit) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_lig_lit-1)) = decomp_k(c,j,i_lig_lit) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_act_som-1)) = decomp_k(c,j,i_act_som) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_slo_som-1)) = decomp_k(c,j,i_slo_som) * dt
               Ksoil%DM(c,j+nlevdecomp*(i_pas_som-1)) = decomp_k(c,j,i_pas_som) * dt
               ! same for cwd but only if fates is not enabled; fates handles CWD
               ! on its own structure
               if (.not. use_fates) then
                  Ksoil%DM(c,j+nlevdecomp*(i_cwd-1))   = decomp_k(c,j,i_cwd) * dt
               end if
            end if !use_soil_matrixcn
         end do
      end do

      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = 1.0_r8
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s1) = 1.0_r8
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s2) = 1.0_r8
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = f_s2s1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = f_s2s3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = 1.0_r8
      if (.not. use_fates) then
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = cwd_fcel
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = cwd_flig
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = rf_cwdl3
      end if
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = rf_l1s1
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s1) = rf_l2s1
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s2) = rf_l3s2
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = rf_s2s1
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = rf_s2s3
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = rf_s3s1

    end associate

 end subroutine decomp_rate_constants_bgc

end module SoilBiogeochemDecompCascadeBGCMod
