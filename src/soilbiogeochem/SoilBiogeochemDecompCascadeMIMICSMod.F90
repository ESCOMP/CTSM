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
  use clm_varpar                         , only : i_litr_min, i_litr_max, i_met_lit, i_cwd
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
  public :: init_decompcascade_mimics       ! Initialization
  public :: decomp_rates_mimics             ! Figure out decomposition rates
  !
  ! !PUBLIC DATA MEMBERS 
  !
  ! !PRIVATE DATA MEMBERS 
  integer, private :: i_pro_som  ! index of protected Soil Organic Matter (SOM)
  integer, private :: i_rec_som  ! index of recalcitrant SOM
  integer, private :: i_avl_som  ! index of available SOM
  integer, private :: i_str_lit  ! index of structural litter pool
  integer, private :: i_cop_mic  ! index of copiotrophic microbial pool
  integer, private :: i_oli_mic  ! index of oligotrophic microbial pool

  type, private :: params_type
     real(r8):: rf_s2s3_bgc    
     ! Respiration fraction cwd --> structural litter (keep comment sep for dif)
     real(r8):: rf_cwdl2_bgc 
     real(r8):: tau_s1_bgc    ! 1/turnover time of  SOM 1 from Century (1/7.3) (1/yr)
     real(r8):: tau_s2_bgc    ! 1/turnover time of  SOM 2 from Century (1/0.2) (1/yr)
     real(r8):: tau_cwd_bgc   ! corrected fragmentation rate constant CWD, century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1 (1/0.3) (1/yr)

     real(r8) :: cwd_fcel_bgc !cellulose fraction for CWD
     real(r8) :: cwd_flig_bgc !

     real(r8) :: minpsi_bgc   !minimum soil water potential for heterotrophic resp
     real(r8) :: maxpsi_bgc   !maximum soil water potential for heterotrophic resp

     ! TODO Need for spin-ups?
     real(r8) :: initial_Cstocks_depth      ! Soil depth for initial Carbon stocks for a cold-start
     real(r8), allocatable :: initial_Cstocks(:)  ! Initial Carbon stocks for a cold-start
     real(r8), allocatable :: mge(:)  ! Microbial growth efficiency (mg/mg)
     real(r8), allocatable :: vmod(:)  ! 
     real(r8), allocatable :: vslope(:)  ! 
     real(r8), allocatable :: vint(:)  ! 
     real(r8), allocatable :: kmod(:)  ! 
     real(r8), allocatable :: kslope(:)  ! 
     real(r8), allocatable :: kint(:)  ! 
     
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
    character(len=32)  :: subname = 'readParams'
    character(len=100) :: errCode = 'Error reading MIMICS params '
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! Read off of netcdf file
    ! TODO Add new params here and in the params file.
    ! TODO tau_s1, tau_s2, etc may need _bgc added/removed here/params file
    ! TODO Read MIMICS-specific params here & shared ones (eg *_bgc) like this:
    !      decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding
    ! TODO If Melannie's testbed values differ from the BGC values,
    !      make separate _mimics params.
    ! TODO When ready for all final params values, talk to @wwieder
    tString='tau_s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s1_bgc=tempr

    tString='tau_s2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s2_bgc=tempr

    tString='tau_cwd_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_cwd_bgc=tempr

    tString='rf_s2s3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s2s3_bgc=tempr

    ! MIMICS value = 0 (same as BGC value)
    tString='rf_cwdl2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl2_bgc=tempr

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
    params_inst%cwd_flig_bgc=tempr 

    tString='initial_Cstocks_depth_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%initial_Cstocks_depth=tempr

    allocate(params_inst%initial_Cstocks(ndecomp_pools_max))
    tString='initial_Cstocks_bgc'
    call ncd_io(trim(tString), params_inst%initial_Cstocks(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%mge(ndecomp_pools_max))
    tString='mge_mimics'
    call ncd_io(trim(tString), params_inst%mge(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%vmod(ndecomp_pools_max))
    tString='vmod_mimics'
    call ncd_io(trim(tString), params_inst%vmod(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%vslope(ndecomp_pools_max))
    tString='vslope_mimics'
    call ncd_io(trim(tString), params_inst%vslope(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%vint(ndecomp_pools_max))
    tString='vint_mimics'
    call ncd_io(trim(tString), params_inst%vint(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%kmod(ndecomp_pools_max))
    tString='kmod_mimics'
    call ncd_io(trim(tString), params_inst%kmod(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%kslope(ndecomp_pools_max))
    tString='kslope_mimics'
    call ncd_io(trim(tString), params_inst%kslope(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

    allocate(params_inst%kint(ndecomp_pools_max))
    tString='kint_mimics'
    call ncd_io(trim(tString), params_inst%kint(:), 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_mimics(bounds, soilbiogeochem_state_inst, soilstate_inst )
    !
    ! !DESCRIPTION:
    ! initialize rate constants and decomposition pathways following the
    ! decomposition cascade of the MIMICS model.
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
    real(r8) :: rf_l1m1
    real(r8) :: rf_l2m1
    real(r8) :: rf_s3m1
    real(r8) :: rf_l1m2
    real(r8) :: rf_l2m2
    real(r8) :: rf_s3m2
    real(r8), allocatable :: rf_s1s3(:,:)
    real(r8) :: rf_s2s3
    real(r8) :: rf_cwdl2
    real(r8) :: cwd_fcel
    real(r8) :: cwd_flig
    real(r8), allocatable :: f_s1s3(:,:)
    real(r8) :: f_s2s3
    real(r8), allocatable :: p_scalar(:,:)

    integer :: i_l1m1
    integer :: i_l2m1
    integer :: i_s3m1
    integer :: i_l1m2
    integer :: i_l2m2
    integer :: i_s3m2
    integer :: i_s1s3
    integer :: i_s2s3
    integer :: i_cwdl2
    real(r8):: speedup_fac                  ! acceleration factor, higher when vertsoilc = .true.

    integer  :: c, j    ! indices
    real(r8) :: t       ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                                     &
         rf_decomp_cascade              => soilbiogeochem_state_inst%rf_decomp_cascade_col       , & ! Input:  [real(r8)          (:,:,:) ]  respired fraction in decomposition step (frac)       
         pathfrac_decomp_cascade        => soilbiogeochem_state_inst%pathfrac_decomp_cascade_col , & ! Input:  [real(r8)          (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)

         cellsand                       => soilstate_inst%cellsand_col                           , & ! Input:  [real(r8)          (:,:)   ]  column 3D sand                                         
         cellclay                       => soilstate_inst%cellclay_col                           , & ! Input:  [real(r8)          (:,:)   ]  column 3D clay
         
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

      allocate(rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(p_scalar(bounds%begc:bounds%endc,1:nlevdecomp))

      !------- time-constant coefficients ---------- !
      ! set respiration fractions for fluxes between compartments
      rf_l1m1 = 1.0_r8 - params_inst%mge(1)
      rf_l2m1 = 1.0_r8 - params_inst%mge(2)
      rf_s3m1 = 1.0_r8 - params_inst%mge(3)
      rf_l1m2 = 1.0_r8 - params_inst%mge(4)
      rf_l2m2 = 1.0_r8 - params_inst%mge(5)
      rf_s3m2 = 1.0_r8 - params_inst%mge(6)

      rf_s2s3 = params_inst%rf_s2s3_bgc
      rf_cwdl2 = params_inst%rf_cwdl2_bgc

      ! set the structural fraction for coarse woody debris
      ! TODO Renamed LHS to cwd_fstr. Make RHS = sum or avg of the two?
      cwd_fstr = min(1.0_r8, &
                     params_inst%cwd_fcel_bgc + params_inst%cwd_flig_bgc)

      ! vmod = "old" vmod * av  AND  kmod = ak / "old" kmod
      ! Table B1 Wieder et al. (2015) and MIMICS params file give different
      ! ak and av values. I used the params file values.
      vmod_l1m1 = params_inst%vmod(1)
      vmod_l2m1 = params_inst%vmod(2)
      vmod_s3m1 = params_inst%vmod(3)
      vmod_l1m2 = params_inst%vmod(4)
      vmod_l2m2 = params_inst%vmod(5)
      vmod_s3m2 = params_inst%vmod(6)
      kmod_l1m1 = params_inst%kmod(1)
      kmod_l2m1 = params_inst%kmod(2)
      kmod_s3m1 = params_inst%kmod(3)
      kmod_l1m2 = params_inst%kmod(4)
      kmod_l2m2 = params_inst%kmod(5)
      kmod_s3m2 = params_inst%kmod(6)
      vslope_l1m1 = params_inst%vslope(1)
      vslope_l2m1 = params_inst%vslope(2)
      vslope_s3m1 = params_inst%vslope(3)
      vslope_l1m2 = params_inst%vslope(4)
      vslope_l2m2 = params_inst%vslope(5)
      vslope_s3m2 = params_inst%vslope(6)
      kslope_l1m1 = params_inst%kslope(1)
      kslope_l2m1 = params_inst%kslope(2)
      kslope_s3m1 = params_inst%kslope(3)
      kslope_l1m2 = params_inst%kslope(4)
      kslope_l2m2 = params_inst%kslope(5)
      kslope_s3m2 = params_inst%kslope(6)
      vint_l1m1 = params_inst%vint(1)
      vint_l2m1 = params_inst%vint(2)
      vint_s3m1 = params_inst%vint(3)
      vint_l1m2 = params_inst%vint(4)
      vint_l2m2 = params_inst%vint(5)
      vint_s3m2 = params_inst%vint(6)
      kint_l1m1 = params_inst%kint(1)
      kint_l2m1 = params_inst%kint(2)
      kint_s3m1 = params_inst%kint(3)
      kint_l1m2 = params_inst%kint(4)
      kint_l2m2 = params_inst%kint(5)
      kint_s3m2 = params_inst%kint(6)

      ! set path fractions
      f_s2s3 = 0.03_r8/(0.45_r8)

      ! some of these are dependent on the soil texture properties
      ! TODO rf_s1s3 -> rf_decomp_cascade as template for p_scalar?
      !      One-time initializations here.
      !      Time-dep params in subr. decomp_rates_mimics.

      do c = bounds%begc, bounds%endc
         do j = 1, nlevdecomp
            t = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - cellsand(c,j))
            f_s1s3(c,j) = .004_r8 / (1._r8 - t)
            rf_s1s3(c,j) = t
            p_scalar(c,j) = (0.8_r8 * exp(-3.0_r8 * sqrt(cellclay(c,j))))** &
                            (-1.0_r8)
         end do
      end do
      initial_stock_soildepth = params_inst%initial_Cstocks_depth

      !-------------------  list of pools and their attributes  ------------
      i_litr_min = 1
      i_met_lit = i_litr_min
      floating_cn_ratio_decomp_pools(i_met_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_met_lit) = 'litr1'
      decomp_cascade_con%decomp_pool_name_history(i_met_lit) = 'LITR1'
      decomp_cascade_con%decomp_pool_name_long(i_met_lit) = 'litter 1'
      decomp_cascade_con%decomp_pool_name_short(i_met_lit) = 'L1'
      is_microbe(i_met_lit) = .false.
      is_litter(i_met_lit) = .true.
      is_soil(i_met_lit) = .false.
      is_cwd(i_met_lit) = .false.
      initial_cn_ratio(i_met_lit) = 10._r8  ! 90 in BGC; not used in MIMICS
      initial_stock(i_met_lit) = params_inst%initial_Cstocks(i_met_lit)
      is_metabolic(i_met_lit) = .true.
      is_cellulose(i_met_lit) = .false.
      is_lignin(i_met_lit) = .false.

      i_str_lit = i_met_lit + 1
      floating_cn_ratio_decomp_pools(i_str_lit) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_str_lit) = 'litr2'
      decomp_cascade_con%decomp_pool_name_history(i_str_lit) = 'LITR2'
      decomp_cascade_con%decomp_pool_name_long(i_str_lit) = 'litter 2'
      decomp_cascade_con%decomp_pool_name_short(i_str_lit) = 'L2'
      is_microbe(i_str_lit) = .false.
      is_litter(i_str_lit) = .true.
      is_soil(i_str_lit) = .false.
      is_cwd(i_str_lit) = .false.
      initial_cn_ratio(i_str_lit) = 10._r8  ! 90 in BGC; not used in MIMICS
      initial_stock(i_str_lit) = params_inst%initial_Cstocks(i_str_lit)
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

      i_pro_som = i_str_lit + 1
      floating_cn_ratio_decomp_pools(i_pro_som) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_pro_som) = 'soil1'
      decomp_cascade_con%decomp_pool_name_history(i_pro_som) = 'SOIL1'
      decomp_cascade_con%decomp_pool_name_long(i_pro_som) = 'soil 1'
      decomp_cascade_con%decomp_pool_name_short(i_pro_som) = 'S1'
      is_microbe(i_pro_som) = .false.
      is_litter(i_pro_som) = .false.
      is_soil(i_pro_som) = .true.
      is_cwd(i_pro_som) = .false.
      initial_cn_ratio(i_pro_som) = 10._r8  ! cn_s1 in BGC; not used in MIMICS
      initial_stock(i_pro_som) = params_inst%initial_Cstocks(i_pro_som)
      is_metabolic(i_pro_som) = .false.
      is_cellulose(i_pro_som) = .false.
      is_lignin(i_pro_som) = .false.

      i_rec_som = i_pro_som + 1
      floating_cn_ratio_decomp_pools(i_rec_som) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_rec_som) = 'soil2'
      decomp_cascade_con%decomp_pool_name_history(i_rec_som) = 'SOIL2'
      decomp_cascade_con%decomp_pool_name_long(i_rec_som) = 'soil 2'
      decomp_cascade_con%decomp_pool_name_short(i_rec_som) = 'S2'
      is_microbe(i_rec_som) = .false.
      is_litter(i_rec_som) = .false.
      is_soil(i_rec_som) = .true.
      is_cwd(i_rec_som) = .false.
      initial_cn_ratio(i_rec_som) = 10._r8  ! cn_s2 in BGC; not used in MIMICS
      initial_stock(i_rec_som) = params_inst%initial_Cstocks(i_rec_som)
      is_metabolic(i_rec_som) = .false.
      is_cellulose(i_rec_som) = .false.
      is_lignin(i_rec_som) = .false.

      i_avl_som = i_rec_som + 1
      floating_cn_ratio_decomp_pools(i_avl_som) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_avl_som) = 'soil3'
      decomp_cascade_con%decomp_pool_name_history(i_avl_som) = 'SOIL3'
      decomp_cascade_con%decomp_pool_name_long(i_avl_som) = 'soil 3'
      decomp_cascade_con%decomp_pool_name_short(i_avl_som) = 'S3'
      is_microbe(i_avl_som) = .false.
      is_litter(i_avl_som) = .false.
      is_soil(i_avl_som) = .true.
      is_cwd(i_avl_som) = .false.
      initial_cn_ratio(i_avl_som) = 10._r8  ! cn_s3 in BGC; not used in MIMICS
      initial_stock(i_avl_som) = params_inst%initial_Cstocks(i_avl_som)
      is_metabolic(i_avl_som) = .false.
      is_cellulose(i_avl_som) = .false.
      is_lignin(i_avl_som) = .false.

      i_cop_mic = i_avl_som + 1
      floating_cn_ratio_decomp_pools(i_cop_mic) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_cop_mic) = 'micr1'
      decomp_cascade_con%decomp_pool_name_history(i_cop_mic) = 'MICR1'
      decomp_cascade_con%decomp_pool_name_long(i_cop_mic) = 'microbial 1'
      decomp_cascade_con%decomp_pool_name_short(i_cop_mic) = 'M1'
      is_microbe(i_cop_mic) = .true.
      is_litter(i_cop_mic) = .false.
      is_soil(i_cop_mic) = .false.
      is_cwd(i_cop_mic) = .false.
      initial_cn_ratio(i_cop_mic) = 10._r8  ! MIMICS may use this
      initial_stock(i_cop_mic) = params_inst%initial_Cstocks(i_cop_mic)
      is_metabolic(i_cop_mic) = .false.
      is_cellulose(i_cop_mic) = .false.
      is_lignin(i_cop_mic) = .false.

      i_oli_mic = i_cop_mic + 1
      floating_cn_ratio_decomp_pools(i_oli_mic) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_oli_mic) = 'micr2'
      decomp_cascade_con%decomp_pool_name_history(i_oli_mic) = 'MICR2'
      decomp_cascade_con%decomp_pool_name_long(i_oli_mic) = 'microbial 2'
      decomp_cascade_con%decomp_pool_name_short(i_oli_mic) = 'M2'
      is_microbe(i_oli_mic) = .true.
      is_litter(i_oli_mic) = .false.
      is_soil(i_oli_mic) = .false.
      is_cwd(i_oli_mic) = .false.
      initial_cn_ratio(i_oli_mic) = 10._r8  ! MIMICS may use this
      initial_stock(i_oli_mic) = params_inst%initial_Cstocks(i_oli_mic)
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
         initial_stock(i_cwd) = params_inst%initial_Cstocks(i_cwd)
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
         spinup_factor(i_cwd) = max(1._r8, (speedup_fac * params_inst%tau_cwd_bgc / 2._r8 ))
      end if
      !som1,2,3
      spinup_factor(i_pro_som) = 1._r8
      spinup_factor(i_rec_som) = 1._r8  ! BGC used cwd formula above but
      spinup_factor(i_avl_som) = 1._r8  ! ...w the respective tau_s values
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
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m1) = rf_l1m1
      cascade_donor_pool(i_l1m1) = i_met_lit
      cascade_receiver_pool(i_l1m1) = i_cop_mic
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m1) = 1.0_r8

      i_l2m1 = 2
      decomp_cascade_con%cascade_step_name(i_l2m1) = 'L2M1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m1) = rf_l2m1
      cascade_donor_pool(i_l2m1) = i_str_lit
      cascade_receiver_pool(i_l2m1) = i_cop_mic
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m1)= 1.0_r8

      i_s3m1 = 3
      decomp_cascade_con%cascade_step_name(i_s3m1) = 'S3M1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3m1) = rf_s3m1
      cascade_donor_pool(i_s3m1) = i_avl_som
      cascade_receiver_pool(i_s3m1) = i_cop_mic
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3m1) = 1.0_r8

      i_l1m2 = 4
      decomp_cascade_con%cascade_step_name(i_l1m2) = 'L1M2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m2) = rf_l1m2
      cascade_donor_pool(i_l1m2) = i_met_lit
      cascade_receiver_pool(i_l1m2) = i_oli_mic
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m2) = 1.0_r8

      i_l2m2 = 5
      decomp_cascade_con%cascade_step_name(i_l2m2) = 'L2M2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m2) = rf_l2m2
      cascade_donor_pool(i_l2m2) = i_str_lit
      cascade_receiver_pool(i_l2m2) = i_oli_mic
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m2)= 1.0_r8

      i_s3m2 = 6
      decomp_cascade_con%cascade_step_name(i_s3m2) = 'S3M2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3m2) = rf_s3m2
      cascade_donor_pool(i_s3m2) = i_avl_som
      cascade_receiver_pool(i_s3m2) = i_oli_mic
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3m2) = 1.0_r8

      ! TODO If keeping s1_s3 & s2_s3 code, is pathfrac correct or should = 1?
      i_s1s3 = 7
      decomp_cascade_con%cascade_step_name(i_s1s3) = 'S1S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
      cascade_donor_pool(i_s1s3) = i_pro_som
      cascade_receiver_pool(i_s1s3) = i_avl_som
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)

      i_s2s3 = 8
      decomp_cascade_con%cascade_step_name(i_s2s3) = 'S2S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = rf_s2s3
      cascade_donor_pool(i_s2s3) = i_rec_som
      cascade_receiver_pool(i_s2s3) = i_avl_som
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = f_s2s3

      if (.not. use_fates) then
         i_cwdl2 = 9
         decomp_cascade_con%cascade_step_name(i_cwdl2) = 'CWDL2'
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
         cascade_donor_pool(i_cwdl2) = i_cwd
         cascade_receiver_pool(i_cwdl2) = i_str_lit
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = cwd_fstr
      end if

      deallocate(rf_s1s3)
      deallocate(f_s1s3)
      deallocate(p_scalar)
      deallocate(params_inst%mge)
      deallocate(params_inst%vmod)
      deallocate(params_inst%vslope)
      deallocate(params_inst%vint)
      deallocate(params_inst%kmod)
      deallocate(params_inst%kslope)
      deallocate(params_inst%kint)
      deallocate(params_inst%initial_Cstocks)

    end associate

  end subroutine init_decompcascade_mimics

  !-----------------------------------------------------------------------
  subroutine decomp_rates_mimics(bounds, num_soilc, filter_soilc, &
       soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! Calculate rates and decomposition pathways for the MIMICS
    ! decomposition cascade model
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year
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
    real(r8):: k_s1_s3  ! decomposition rate SOM 1 (1/sec)
    real(r8):: k_s2_s3  ! decomposition rate SOM 2 (1/sec)
    real(r8):: k_frag                       ! fragmentation rate constant CWD (1/sec)
    real(r8):: vmax_l1_m1  !
    real(r8):: vmax_l2_m1  !
    real(r8):: vmax_s3_m1  !
    real(r8):: vmax_l1_m2  !
    real(r8):: vmax_l2_m2  !
    real(r8):: vmax_s3_m2  !
    real(r8):: km_l1_m1  !
    real(r8):: km_l2_m1  !
    real(r8):: km_s3_m1  !
    real(r8):: km_l1_m2  !
    real(r8):: km_l2_m2  !
    real(r8):: km_s3_m2  !
    real(r8):: tau_m1_s1  !
    real(r8):: tau_m1_s2  !
    real(r8):: tau_m1_s3  !
    real(r8):: tau_m2_s1  !
    real(r8):: tau_m2_s2  !
    real(r8):: tau_m2_s3  !
    real(r8):: decomp_depth_efolding        ! (meters) e-folding depth for reduction in decomposition [
    integer :: c, fc, j, k, l
    real(r8):: days_per_year                ! days per year
    real(r8):: depth_scalar(bounds%begc:bounds%endc,1:nlevdecomp) 
    real(r8):: w_d_o_scalars  ! product of w_scalar * depth_scalar * o_scalar
    real(r8):: mino2lim                     !minimum anaerobic decomposition rate
    real(r8):: spinup_geogterm_l1(bounds%begc:bounds%endc) ! geographically-varying spinup term for l1
    real(r8):: spinup_geogterm_l2(bounds%begc:bounds%endc)  ! geographically-varying spinup term for l2
    real(r8):: spinup_geogterm_cwd(bounds%begc:bounds%endc) ! geographically-varying spinup term for cwd
    real(r8):: spinup_geogterm_s1(bounds%begc:bounds%endc) ! geographically-varying spinup term for s1
    real(r8):: spinup_geogterm_s2(bounds%begc:bounds%endc) ! geographically-varying spinup term for s2
    real(r8):: spinup_geogterm_s3(bounds%begc:bounds%endc) ! geographically-varying spinup term for s3
    real(r8):: spinup_geogterm_m1(bounds%begc:bounds%endc)  ! geographically-varying spinup term for m1
    real(r8):: spinup_geogterm_m2(bounds%begc:bounds%endc)  ! geographically-varying spinup term for m2

    !-----------------------------------------------------------------------

    associate(                                                           &
         minpsi         => params_inst%minpsi_bgc                      , & ! Input:  [real(r8)         ]  minimum soil suction (mm)
         maxpsi         => params_inst%maxpsi_bgc                      , & ! Input:  [real(r8)         ]  maximum soil suction (mm)
         soilpsi        => soilstate_inst%soilpsi_col                  , & ! Input:  [real(r8) (:,:)   ]  soil water potential in each soil layer (MPa)          

         t_soisno       => temperature_inst%t_soisno_col               , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       

         o2stress_sat   => ch4_inst%o2stress_sat_col                   , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         o2stress_unsat => ch4_inst%o2stress_unsat_col                 , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         finundated     => ch4_inst%finundated_col                     , & ! Input:  [real(r8) (:)     ]  fractional inundated area                                
         
         w_scalar       => soilbiogeochem_carbonflux_inst%w_scalar_col , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
         o_scalar       => soilbiogeochem_carbonflux_inst%o_scalar_col , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
         decomp_k       => soilbiogeochem_carbonflux_inst%decomp_k_col , & ! Output: [real(r8) (:,:,:) ]  rate for decomposition (1./sec)
         spinup_factor  => decomp_cascade_con%spinup_factor              & ! Input:  [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
         )

      mino2lim = CNParamsShareInst%mino2lim

      days_per_year = get_days_per_year()

      ! Set "decomp_depth_efolding" parameter
      decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding

      ! translate to per-second time constant
      ! TODO vmax, km, tau terms replaced (most?) k_ terms; CONFIRM s-1 UNITS
      k_s1_s3 = 1._r8 / (secspday * days_per_year * params_inst%tau_s1_bgc)
      k_s2_s3 = 1._r8 / (secspday * days_per_year * params_inst%tau_s2_bgc)
      k_frag = 1._r8 / (secspday * days_per_year * params_inst%tau_cwd_bgc)

     ! calc ref rate
      if ( spinup_state >= 1 ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            !
            if ( abs(spinup_factor(i_met_lit) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_l1(c) = spinup_factor(i_met_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_str_lit) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_l2(c) = spinup_factor(i_str_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_l2(c) = 1._r8
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
            if ( abs(spinup_factor(i_pro_som) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_s1(c) = spinup_factor(i_pro_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_rec_som) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_s2(c) = spinup_factor(i_rec_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s2(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_avl_som) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_s3(c) = spinup_factor(i_avl_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_s3(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_cop_mic) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_m1(c) = spinup_factor(i_cop_mic) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_m1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_oli_mic) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_m2(c) = spinup_factor(i_oli_mic) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_m2(c) = 1._r8
            endif
            !
         end do
      else
         do fc = 1,num_soilc
            c = filter_soilc(fc)
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

         ! TODO @wwieder concerned about breaking ch4 code
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

      ! Term that reduces decomposition rate at depth
      do j = 1, nlevdecomp
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            if (use_vertsoilc) then
               ! Using fixed e-folding depth as in
               ! SoilBiogeochemDecompCascadeBGCMod.F90
               depth_scalar(c,j) = exp(-zsoi(j) / decomp_depth_efolding)
            else
               depth_scalar(c,j) = 1.0_r8
            end if
         end do
      end do

      ! calculate rates for all litter and som pools
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            ! vmax, tau, and km are time-dependent terms.
            ! Table B1 Wieder et al. (2015) & MIMICS params file give different
            ! kslope. I used the params file value(s).
            ! TODO Introduce t_soi24_degC_col following t_veg24_patch (running
            !        mean) or t_mo_patch (30-d avg) or ...
            vmax_l1_m1 = exp(vslope_l1_m1 * t_soi24_degC_col(c,j) + vint_l1_m1) * vmod_l1_m1
            vmax_l1_m2 = exp(vslope_l1_m2 * t_soi24_degC_col(c,j) + vint_l1_m2) * vmod_l1_m2
            vmax_l2_m1 = exp(vslope_l2_m1 * t_soi24_degC_col(c,j) + vint_l2_m1) * vmod_l2_m1
            vmax_l2_m2 = exp(vslope_l2_m2 * t_soi24_degC_col(c,j) + vint_l2_m2) * vmod_l2_m2
            vmax_s3_m1 = exp(vslope_s3_m1 * t_soi24_degC_col(c,j) + vint_s3_m1) * vmod_l2_m1
            vmax_s3_m2 = exp(vslope_s3_m2 * t_soi24_degC_col(c,j) + vint_s3_m2) * vmod_l2_m2

            km_l1_m1 = exp(kslope_l1_m1 * t_soi24_degC_col(c,j) + kint_l1_m1) * kmod_l1_m1
            km_l1_m2 = exp(kslope_l1_m2 * t_soi24_degC_col(c,j) + kint_l1_m2) * kmod_l1_m2
            km_l2_m1 = exp(kslope_l2_m1 * t_soi24_degC_col(c,j) + kint_l2_m1) * kmod_l2_m1
            km_l2_m2 = exp(kslope_l2_m2 * t_soi24_degC_col(c,j) + kint_l2_m2) * kmod_l2_m2
            km_s3_m1 = exp(kslope_s3_m1 * t_soi24_degC_col(c,j) + kint_s3_m1) * kmod_s3_m1 / p_scalar(c,j)
            km_s3_m2 = exp(kslope_s3_m2 * t_soi24_degC_col(c,j) + kint_s3_m2) * kmod_s3_m2 / p_scalar(c,j)

! TODO mapping tau params to Wieder et al. 2015 & testbed code
! --------------------------------------------------------------
! tau* are time-dependent

! Solving for tauR (tau_m1_s* here) & tauK (tau_m2_s* here)

! fmet_p(1) = 0.85, fmet_p(2) = 0.85, fmet_p(3) = 0.013
! ligninNratioAvg(npt) =  &
!  ( ligninNratio(npt,leaf)  * (cleaf2met(npt) + cleaf2str(npt)) &
!  + ligninNratio(npt,froot) * (croot2met(npt) + croot2str(npt)) &
!  + ligninNratio(npt,wood)  * (cwd2str(npt)) ) &
!  / max(0.001, cleaf2met(npt)+cleaf2str(npt)+croot2met(npt)+croot2str(npt)+cwd2str(npt))
! set limits on Lignin:N to keep fmet > 0.2
! necessary for litter quality in boreal forests with high cwd flux
! ligninNratioAvg(npt) = min(40.0, ligninNratioAvg(npt))
! fmet(npt) = fmet_p(1) * (fmet_p(2) - fmet_p(3) * ligninNratioAvg(npt))

! tauModDenom = 100, tauMod_min = 0.8, tauMod_max = 1.2
! tauMod(npt) = sqrt(avg_ann_npp_gC_m2_yr / tauModDenom)
! tauMod(npt) = min(tauMod_max, max(tauMod_min, tauMod(npt)))

! tau_r(1) = 5.2e-4_r8, tau_r(2) = 0.3_r8 or 0.4?
! tauR(npt) = tau_r(1) * exp(tau_r(2) * fmet(npt)) * tauMod(npt)

! tau_k(1) = 2.4e-4_r8, tau_k(2) = 0.1_r8
! tauK(npt) = tau_k(1) * exp(tau_k(2) * fmet(npt)) * tauMod(npt)
! --------------------------------------------------------------
! TODO STILL MISSING N-related stuff: DIN...

            ! Product of w_scalar * depth_scalar * o_scalar
            w_d_o_scalars = w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)

            ! decomp_k used in SoilBiogeochemPotentialMod.F90
            ! TODO Defined correctly here?
            ! TODO Combine multiple terms of same donor pools...
            ! TODO Where to calculate m1_pool and m2_pool?
            ! TODO When to introduce litter/soil biomass to the calculation?
            ! TODO Excluded spinup terms for now. See cwd below how used.
            decomp_k(c,j,i_met_lit) = (vmax_l1_m1 * m1_pool(c,j) / &
                                       (km_l1_m1 + m1_pool(c,j))) * &
                                      w_d_o_scalars
            decomp_k(c,j,i_met_lit) = (vmax_l1_m2 * m2_pool(c,j) / &
                                       (km_l1_m2 + m2_pool(c,j))) * &
                                      w_d_o_scalars

            decomp_k(c,j,i_str_lit) = (vmax_l2_m1 * m1_pool(c,j) / &
                                       (km_l2_m1 + m1_pool(c,j))) * &
                                      w_d_o_scalars
            decomp_k(c,j,i_str_lit) = (vmax_l2_m2 * m2_pool(c,j) / &
                                       (km_l2_m2 + m2_pool(c,j))) * &
                                      w_d_o_scalars

            decomp_k(c,j,i_avl_som) = (vmax_s3_m1 * m1_pool(c,j) / &
                                       (km_s3_m1 + m1_pool(c,j))) * &
                                      w_d_o_scalars
            decomp_k(c,j,i_avl_som) = (vmax_s3_m2 * m2_pool(c,j) / &
                                       (km_s3_m2 + m2_pool(c,j))) * &
                                      w_d_o_scalars

            decomp_k(c,j,i_pro_som) = k_s1_s3 * w_d_o_scalars
            decomp_k(c,j,i_rec_som) = k_s2_s3 * w_d_o_scalars

            decomp_k(c,j,i_cop_mic) = tau_m1_s1 * m1_pool(c,j)**(densdep) * &
                                      fphys(c) * w_d_o_scalars
            decomp_k(c,j,i_cop_mic) = tau_m1_s2 * m1_pool(c,j)**(densdep) * &
                                      fchem(c) * w_d_o_scalars
            decomp_k(c,j,i_cop_mic) = tau_m1_s3 * m1_pool(c,j)**(densdep) * &
                                      faval(c) *w_d_o_scalars
            decomp_k(c,j,i_oli_mic) = tau_m2_s1 * m2_pool(c,j)**(densdep) * &
                                      fphys(c) * w_d_o_scalars
            decomp_k(c,j,i_oli_mic) = tau_m2_s2 * m2_pool(c,j)**(densdep) * &
                                      fchem(c) * w_d_o_scalars
            decomp_k(c,j,i_oli_mic) = tau_m2_s3 * m2_pool(c,j)**(densdep) * &
                                      faval(c) * w_d_o_scalars
            ! Same for cwd but only if fates not enabled; fates handles cwd on
            ! its own structure
            ! TODO Kept cwd bc we do not have a MIMICS replacement for it.
            !      Also this shows what the BGC decomp_k formulas looked like.
            if (.not. use_fates) then
               decomp_k(c,j,i_cwd) = k_frag * w_d_o_scalars * spinup_geogterm_cwd(c)
            end if
         end do
      end do

    end associate

 end subroutine decomp_rates_mimics

end module SoilBiogeochemDecompCascadeMIMICSMod
