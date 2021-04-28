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
  use clm_varpar                         , only : nlevsoi, nlevgrnd
  use clm_varpar                         , only : nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools, ndecomp_pools_max
  use clm_varpar                         , only : i_litr_min, i_litr_max, i_met_lit, i_litr2, i_litr3, i_cwd
  use clm_varctl                         , only : iulog, spinup_state, anoxia, use_lch4, use_vertsoilc, use_fates
  use clm_varcon                         , only : zsoi
  use decompMod                          , only : bounds_type
  use spmdMod                            , only : masterproc
  use abortutils                         , only : endrun
  use CNSharedParamsMod                  , only : CNParamsShareInst, nlev_soildecomp_standard 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
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
  integer, private            :: i_soil1   = -9         ! Soil Organic Matter (SOM) first pool
  integer, private            :: i_soil2   = -9         ! SOM second pool
  integer, private            :: i_soil3   = -9         ! SOM third pool
  integer, private, parameter :: i_litr1   = i_met_lit  ! First litter pool, metabolic

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

     real(r8):: rf_cwdl2_bgc 
     real(r8):: rf_cwdl3_bgc

     real(r8):: tau_l1_bgc    ! 1/turnover time of  litter 1 from Century (l/18.5) (1/yr)
     real(r8):: tau_l2_l3_bgc ! 1/turnover time of  litter 2 and litter 3 from Century (1/4.9) (1/yr)
     real(r8):: tau_s1_bgc    ! 1/turnover time of  SOM 1 from Century (1/7.3) (1/yr)
     real(r8):: tau_s2_bgc    ! 1/turnover time of  SOM 2 from Century (1/0.2) (1/yr)
     real(r8):: tau_s3_bgc    ! 1/turnover time of  SOM 3 from Century (1/0.0045) (1/yr)
     real(r8):: tau_cwd_bgc   ! corrected fragmentation rate constant CWD, century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1 (1/0.3) (1/yr)

     real(r8) :: cwd_fcel_bgc !cellulose fraction for CWD
     real(r8) :: cwd_flig

     real(r8) :: minpsi_bgc   !minimum soil water potential for heterotrophic resp
     real(r8) :: maxpsi_bgc   !maximum soil water potential for heterotrophic resp

     real(r8), allocatable :: initial_Cstocks(:)  ! Initial Carbon stocks for a cold-start
     real(r8) :: initial_Cstocks_depth      ! Soil depth for initial Carbon stocks for a cold-start
     
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
    tString='tau_l1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_l1_bgc=tempr

    tString='tau_l2_l3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_l2_l3_bgc=tempr

    tString='tau_s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s1_bgc=tempr

    tString='tau_s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s2_bgc=tempr

    tString='tau_s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s3_bgc=tempr

    tString='tau_cwd_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_cwd_bgc=tempr

    tString='cn_s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s1_bgc=tempr

    tString='cn_s2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s2_bgc=tempr

    tString='cn_s3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s3_bgc=tempr

    tString='rf_l1s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l1s1_bgc=tempr

    tString='rf_l2s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l2s1_bgc=tempr

    tString='rf_l3s2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l3s2_bgc=tempr   

    tString='rf_s2s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s2s1_bgc=tempr

    tString='rf_s2s3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s2s3_bgc=tempr

    tString='rf_s3s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s3s1_bgc=tempr

    tString='rf_cwdl2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl2_bgc=tempr

    tString='rf_cwdl3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl3_bgc=tempr

    tString='cwd_fcel'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cwd_fcel_bgc=tempr

    tString='minpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%minpsi_bgc=tempr 

    tString='maxpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%maxpsi_bgc=tempr 

    tString='cwd_flig'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cwd_flig=tempr
    
    allocate(params_inst%initial_Cstocks(ndecomp_pools_max))
    tString='initial_Cstocks_bgc'
    call ncd_io(trim(tString), params_inst%initial_Cstocks(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    tString='initial_Cstocks_depth_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%initial_Cstocks_depth=tempr

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_bgc(bounds, soilbiogeochem_state_inst, soilstate_inst )
    !
    ! !DESCRIPTION:
    !  initialize rate constants and decomposition pathways following the decomposition cascade of the BGC model.
    !  written by C. Koven 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds  
    type(soilbiogeochem_state_type) , intent(inout) :: soilbiogeochem_state_inst
    type(soilstate_type)            , intent(in)    :: soilstate_inst
    !
    ! !LOCAL VARIABLES
    !-- properties of each decomposing pool
    real(r8) :: rf_l1s1
    real(r8) :: rf_l2s1
    real(r8) :: rf_l3s2
    !real(r8) :: rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
    !real(r8) :: rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8), allocatable :: rf_s1s2(:,:)
    real(r8), allocatable :: rf_s1s3(:,:)
    real(r8) :: rf_s2s1
    real(r8) :: rf_s2s3
    real(r8) :: rf_s3s1
    real(r8) :: rf_cwdl2
    real(r8) :: rf_cwdl3
    real(r8) :: cwd_fcel
    real(r8) :: cwd_flig
    real(r8) :: cn_s1
    real(r8) :: cn_s2
    real(r8) :: cn_s3
    !real(r8) :: f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
    !real(r8) :: f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8), allocatable :: f_s1s2(:,:)
    real(r8), allocatable :: f_s1s3(:,:)
    real(r8) :: f_s2s1
    real(r8) :: f_s2s3

    integer :: i_l1s1
    integer :: i_l2s1
    integer :: i_l3s2
    integer :: i_s1s2
    integer :: i_s1s3
    integer :: i_s2s1
    integer :: i_s2s3
    integer :: i_s3s1
    integer :: i_cwdl2
    integer :: i_cwdl3
    real(r8):: speedup_fac                  ! acceleration factor, higher when vertsoilc = .true.

    integer  :: c, j    ! indices
    real(r8) :: t       ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                                     &
         rf_decomp_cascade              => soilbiogeochem_state_inst%rf_decomp_cascade_col       , & ! Input:  [real(r8)          (:,:,:) ]  respired fraction in decomposition step (frac)       
         pathfrac_decomp_cascade        => soilbiogeochem_state_inst%pathfrac_decomp_cascade_col , & ! Input:  [real(r8)          (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)

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

      rf_cwdl2 = params_inst%rf_cwdl2_bgc
      rf_cwdl3 = params_inst%rf_cwdl3_bgc

      ! set the cellulose and lignin fractions for coarse woody debris
      cwd_fcel = params_inst%cwd_fcel_bgc
      cwd_flig = params_inst%cwd_flig

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
      initial_stock_soildepth = params_inst%initial_Cstocks_depth

      !-------------------  list of pools and their attributes  ------------
      i_litr_min = 1
      floating_cn_ratio_decomp_pools(i_litr1) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_litr1) = 'litr1'
      decomp_cascade_con%decomp_pool_name_history(i_litr1) = 'LITR1'
      decomp_cascade_con%decomp_pool_name_long(i_litr1) = 'litter 1'
      decomp_cascade_con%decomp_pool_name_short(i_litr1) = 'L1'
      is_litter(i_litr1) = .true.
      is_soil(i_litr1) = .false.
      is_cwd(i_litr1) = .false.
      initial_cn_ratio(i_litr1) = 90._r8
      initial_stock(i_litr1) = params_inst%initial_Cstocks(i_litr1)
      is_metabolic(i_litr1) = .true.
      is_cellulose(i_litr1) = .false.
      is_lignin(i_litr1) = .false.

      i_litr2 = i_litr1 + 1
      floating_cn_ratio_decomp_pools(i_litr2) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_litr2) = 'litr2'
      decomp_cascade_con%decomp_pool_name_history(i_litr2) = 'LITR2'
      decomp_cascade_con%decomp_pool_name_long(i_litr2) = 'litter 2'
      decomp_cascade_con%decomp_pool_name_short(i_litr2) = 'L2'
      is_litter(i_litr2) = .true.
      is_soil(i_litr2) = .false.
      is_cwd(i_litr2) = .false.
      initial_cn_ratio(i_litr2) = 90._r8
      initial_stock(i_litr2) = params_inst%initial_Cstocks(i_litr2)
      is_metabolic(i_litr2) = .false.
      is_cellulose(i_litr2) = .true.
      is_lignin(i_litr2) = .false.

      i_litr3 = i_litr2 + 1
      floating_cn_ratio_decomp_pools(i_litr3) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_litr3) = 'litr3'
      decomp_cascade_con%decomp_pool_name_history(i_litr3) = 'LITR3'
      decomp_cascade_con%decomp_pool_name_long(i_litr3) = 'litter 3'
      decomp_cascade_con%decomp_pool_name_short(i_litr3) = 'L3'
      is_litter(i_litr3) = .true.
      is_soil(i_litr3) = .false.
      is_cwd(i_litr3) = .false.
      initial_cn_ratio(i_litr3) = 90._r8
      initial_stock(i_litr3) = params_inst%initial_Cstocks(i_litr3)
      is_metabolic(i_litr3) = .false.
      is_cellulose(i_litr3) = .false.
      is_lignin(i_litr3) = .true.

      i_litr_max = i_litr3
      if (i_litr_max > 3) then
         call endrun(msg='ERROR: expecting i_litr_max <= 3; see pftconMod '//&
              errMsg(sourcefile, __LINE__))
      end if

      i_soil1 = i_litr3 + 1
      floating_cn_ratio_decomp_pools(i_soil1) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_soil1) = 'soil1'
      decomp_cascade_con%decomp_pool_name_history(i_soil1) = 'SOIL1'
      decomp_cascade_con%decomp_pool_name_long(i_soil1) = 'soil 1'
      decomp_cascade_con%decomp_pool_name_short(i_soil1) = 'S1'
      is_litter(i_soil1) = .false.
      is_soil(i_soil1) = .true.
      is_cwd(i_soil1) = .false.
      initial_cn_ratio(i_soil1) = cn_s1
      initial_stock(i_soil1) = params_inst%initial_Cstocks(i_soil1)
      is_metabolic(i_soil1) = .false.
      is_cellulose(i_soil1) = .false.
      is_lignin(i_soil1) = .false.

      i_soil2 = i_soil1 + 1
      floating_cn_ratio_decomp_pools(i_soil2) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_soil2) = 'soil2'
      decomp_cascade_con%decomp_pool_name_history(i_soil2) = 'SOIL2'
      decomp_cascade_con%decomp_pool_name_long(i_soil2) = 'soil 2'
      decomp_cascade_con%decomp_pool_name_short(i_soil2) = 'S2'
      is_litter(i_soil2) = .false.
      is_soil(i_soil2) = .true.
      is_cwd(i_soil2) = .false.
      initial_cn_ratio(i_soil2) = cn_s2
      initial_stock(i_soil2) = params_inst%initial_Cstocks(i_soil2)
      is_metabolic(i_soil2) = .false.
      is_cellulose(i_soil2) = .false.
      is_lignin(i_soil2) = .false.

      i_soil3 = i_soil2 + 1
      floating_cn_ratio_decomp_pools(i_soil3) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_soil3) = 'soil3'
      decomp_cascade_con%decomp_pool_name_history(i_soil3) = 'SOIL3'
      decomp_cascade_con%decomp_pool_name_long(i_soil3) = 'soil 3'
      decomp_cascade_con%decomp_pool_name_short(i_soil3) = 'S3'
      is_litter(i_soil3) = .false.
      is_soil(i_soil3) = .true.
      is_cwd(i_soil3) = .false.
      initial_cn_ratio(i_soil3) = cn_s3
      initial_stock(i_soil3) = params_inst%initial_Cstocks(i_soil3)
      is_metabolic(i_soil3) = .false.
      is_cellulose(i_soil3) = .false.
      is_lignin(i_soil3) = .false.

      if (.not. use_fates) then
         ! CWD
         i_cwd = i_soil3 + 1
         floating_cn_ratio_decomp_pools(i_cwd) = .true.
         decomp_cascade_con%decomp_pool_name_restart(i_cwd) = 'cwd'
         decomp_cascade_con%decomp_pool_name_history(i_cwd) = 'CWD'
         decomp_cascade_con%decomp_pool_name_long(i_cwd) = 'coarse woody debris'
         decomp_cascade_con%decomp_pool_name_short(i_cwd) = 'CWD'
         is_litter(i_cwd) = .false.
         is_soil(i_cwd) = .false.
         is_cwd(i_cwd) = .true.
         initial_cn_ratio(i_cwd) = 90._r8
         initial_stock(i_cwd) = params_inst%initial_Cstocks(i_cwd)
         is_metabolic(i_cwd) = .false.
         is_cellulose(i_cwd) = .false.
         is_lignin(i_cwd) = .false.
      endif

      speedup_fac = 1._r8

      !lit1
      spinup_factor(i_litr1) = 1._r8
      !lit2,3
      spinup_factor(i_litr2) = 1._r8
      spinup_factor(i_litr3) = 1._r8
      !CWD
      if (.not. use_fates) then
         spinup_factor(i_cwd) = max(1._r8, (speedup_fac * params_inst%tau_cwd_bgc / 2._r8 ))
      end if
      !som1
      spinup_factor(i_soil1) = 1._r8
      !som2,3
      spinup_factor(i_soil2) = max(1._r8, (speedup_fac * params_inst%tau_s2_bgc))
      spinup_factor(i_soil3) = max(1._r8, (speedup_fac * params_inst%tau_s3_bgc))

      if ( masterproc ) then
         write(iulog,*) 'Spinup_state ',spinup_state
         write(iulog,*) 'Spinup factors ',spinup_factor
      end if

      !----------------  list of transitions and their time-independent coefficients  ---------------!
      i_l1s1 = 1
      decomp_cascade_con%cascade_step_name(i_l1s1) = 'L1S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = rf_l1s1
      cascade_donor_pool(i_l1s1) = i_litr1
      cascade_receiver_pool(i_l1s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = 1.0_r8

      i_l2s1 = 2
      decomp_cascade_con%cascade_step_name(i_l2s1) = 'L2S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s1) = rf_l2s1
      cascade_donor_pool(i_l2s1) = i_litr2
      cascade_receiver_pool(i_l2s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s1)= 1.0_r8

      i_l3s2 = 3
      decomp_cascade_con%cascade_step_name(i_l3s2) = 'L3S2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s2) = rf_l3s2
      cascade_donor_pool(i_l3s2) = i_litr3
      cascade_receiver_pool(i_l3s2) = i_soil2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s2) = 1.0_r8

      i_s1s2 = 4
      decomp_cascade_con%cascade_step_name(i_s1s2) = 'S1S2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
      cascade_donor_pool(i_s1s2) = i_soil1
      cascade_receiver_pool(i_s1s2) = i_soil2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)

      i_s1s3 = 5
      decomp_cascade_con%cascade_step_name(i_s1s3) = 'S1S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
      cascade_donor_pool(i_s1s3) = i_soil1
      cascade_receiver_pool(i_s1s3) = i_soil3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)

      i_s2s1 = 6
      decomp_cascade_con%cascade_step_name(i_s2s1) = 'S2S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = rf_s2s1
      cascade_donor_pool(i_s2s1) = i_soil2
      cascade_receiver_pool(i_s2s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = f_s2s1

      i_s2s3 = 7 
      decomp_cascade_con%cascade_step_name(i_s2s3) = 'S2S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = rf_s2s3
      cascade_donor_pool(i_s2s3) = i_soil2
      cascade_receiver_pool(i_s2s3) = i_soil3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = f_s2s3

      i_s3s1 = 8
      decomp_cascade_con%cascade_step_name(i_s3s1) = 'S3S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = rf_s3s1
      cascade_donor_pool(i_s3s1) = i_soil3
      cascade_receiver_pool(i_s3s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = 1.0_r8

      if (.not. use_fates) then
         i_cwdl2 = 9
         decomp_cascade_con%cascade_step_name(i_cwdl2) = 'CWDL2'
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
         cascade_donor_pool(i_cwdl2) = i_cwd
         cascade_receiver_pool(i_cwdl2) = i_litr2
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = cwd_fcel
         
         i_cwdl3 = 10
         decomp_cascade_con%cascade_step_name(i_cwdl3) = 'CWDL3'
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = rf_cwdl3
         cascade_donor_pool(i_cwdl3) = i_cwd
         cascade_receiver_pool(i_cwdl3) = i_litr3
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = cwd_flig
      end if

      deallocate(rf_s1s2)
      deallocate(rf_s1s3)
      deallocate(f_s1s2)
      deallocate(f_s1s3)
      deallocate(params_inst%initial_Cstocks)

    end associate

  end subroutine init_decompcascade_bgc

  !-----------------------------------------------------------------------
  subroutine decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
       soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    !  calculate rate constants and decomposition pathways for the CENTURY decomposition cascade model
    !  written by C. Koven based on original CLM4 decomposition cascade
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year
    use shr_const_mod    , only : SHR_CONST_PI
    use clm_varcon       , only : secspday
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds          
    integer                              , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                              , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)                 , intent(in)    :: soilstate_inst
    type(temperature_type)               , intent(in)    :: temperature_inst
    type(ch4_type)                       , intent(in)    :: ch4_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    !
    ! !LOCAL VARIABLES:
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
    real(r8):: cwdc_loss                    ! fragmentation rate for CWD carbon (gC/m2/s)
    real(r8):: cwdn_loss                    ! fragmentation rate for CWD nitrogen (gN/m2/s)
    real(r8):: Q10                          ! temperature dependence
    real(r8):: froz_q10                     ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
    real(r8):: decomp_depth_efolding        ! (meters) e-folding depth for reduction in decomposition [
    integer :: c, fc, j, k, l
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
         minpsi         => params_inst%minpsi_bgc                      , & ! Input:  [real(r8)         ]  minimum soil suction (mm)
         maxpsi         => params_inst%maxpsi_bgc                      , & ! Input:  [real(r8)         ]  maximum soil suction (mm)
         soilpsi        => soilstate_inst%soilpsi_col                  , & ! Input:  [real(r8) (:,:)   ]  soil water potential in each soil layer (MPa)          

         t_soisno       => temperature_inst%t_soisno_col               , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       

         o2stress_sat   => ch4_inst%o2stress_sat_col                   , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         o2stress_unsat => ch4_inst%o2stress_unsat_col                 , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         finundated     => ch4_inst%finundated_col                     , & ! Input:  [real(r8) (:)     ]  fractional inundated area                                
         
         t_scalar       => soilbiogeochem_carbonflux_inst%t_scalar_col , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp                     
         w_scalar       => soilbiogeochem_carbonflux_inst%w_scalar_col , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
         o_scalar       => soilbiogeochem_carbonflux_inst%o_scalar_col , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
         decomp_k       => soilbiogeochem_carbonflux_inst%decomp_k_col , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         spinup_factor  => decomp_cascade_con%spinup_factor              & ! Input:  [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
         )

      mino2lim = CNParamsShareInst%mino2lim

      if ( use_century_tfunc .and. normalize_q10_to_century_tfunc ) then
         call endrun(msg='ERROR: cannot have both use_century_tfunc and normalize_q10_to_century_tfunc set as true'//&
              errMsg(sourcefile, __LINE__))
      endif

      days_per_year = get_days_per_year()

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
      k_frag = 1._r8  / (secspday * days_per_year * params_inst%tau_cwd_bgc)

     ! calc ref rate
      catanf_30 = catanf(30._r8)

      if ( spinup_state >= 1 ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            !
            if ( abs(spinup_factor(i_litr1) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_l1(c) = spinup_factor(i_litr1) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_litr2) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_l23(c) = spinup_factor(i_litr2) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l23(c) = 1._r8
            endif
            !
            if ( .not. use_fates ) then
               if ( abs(spinup_factor(i_cwd) - 1._r8) .gt. .000001_r8) then
                  spinup_geogterm_cwd(c) = spinup_factor(i_cwd) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
               else
                  spinup_geogterm_cwd(c) = 1._r8
               endif
            endif
            !
            if ( abs(spinup_factor(i_soil1) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_s1(c) = spinup_factor(i_soil1) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_soil2) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_s2(c) = spinup_factor(i_soil2) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s2(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_soil3) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_s3(c) = spinup_factor(i_soil3) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s3(c) = 1._r8
            endif
            !
         end do
      else
         do fc = 1,num_soilc
            c = filter_soilc(fc)
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
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               frw(c) = frw(c) + col%dz(c,j)
            end do
         end do
         do j = 1,nlev_soildecomp_standard
            do fc = 1,num_soilc
               c = filter_soilc(fc)
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
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
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
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
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
            do fc = 1,num_soilc
               c = filter_soilc(fc)
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
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)

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
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (t_soisno(c,j) >= SHR_CONST_TKFRZ) then
                     t_scalar(c,j)= (Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
                  else
                     t_scalar(c,j)= (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))
                  endif
               end do
            end do

         else

            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
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
            do fc = 1,num_soilc
               c = filter_soilc(fc)
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
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)

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
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               t_scalar(c,j) = t_scalar(c,j) * normalization_factor
            end do
         end do
      endif

      ! add a term to reduce decomposition rate at depth
      ! for now used a fixed e-folding depth
      do j = 1, nlevdecomp
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            if (use_vertsoilc) then
               depth_scalar(c,j) = exp(-zsoi(j) / decomp_depth_efolding)
            else
               depth_scalar(c,j) = 1.0_r8
            end if
         end do
      end do

      ! calculate rate constants for all litter and som pools
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_k(c,j,i_litr1) = k_l1    * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l1(c)
            decomp_k(c,j,i_litr2) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l23(c)
            decomp_k(c,j,i_litr3) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l23(c)
            decomp_k(c,j,i_soil1) = k_s1    * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s1(c)
            decomp_k(c,j,i_soil2) = k_s2    * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s2(c)
            decomp_k(c,j,i_soil3) = k_s3    * t_scalar(c,j) * w_scalar(c,j) * &
               depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s3(c)
            ! same for cwd but only if fates is not enabled; fates handles CWD
            ! on its own structure
            if (.not. use_fates) then
               decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_cwd(c)
            end if
         end do
      end do

    end associate

 end subroutine decomp_rate_constants_bgc

end module SoilBiogeochemDecompCascadeBGCMod
