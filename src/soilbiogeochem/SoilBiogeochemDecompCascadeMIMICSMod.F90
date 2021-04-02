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
  use clm_varpar                         , only : nlevdecomp
  use clm_varpar                         , only : i_met_lit
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
  public :: DecompCascadeMIMICSreadNML      ! Read in namelist
  public :: readParams                      ! Read in parameters from params file
  public :: init_decompcascade_mimics       ! Initialization
  public :: decomp_rate_constants_mimics    ! Figure out decomposition rates
  !
  ! !PUBLIC DATA MEMBERS 
  !
  ! !PRIVATE DATA MEMBERS 

  integer, private, parameter :: nsompools = 3          ! Number of SOM pools
  integer, private, parameter :: i_litr1   = i_met_lit  ! First litter pool, metobolic
  integer, private :: i_litr2 = -9  ! Second litter pool, lignin in MIMICS
  integer, private :: i_micr1 = -9  ! First microbial pool, copiotrophic, as labeled in Wieder et al. 2015
  integer, private :: i_micr2 = -9  ! Second microbial pool, oligotrophic, as labeled in Wieder et al. 2015
  integer, private :: i_cwd   = -9  ! Coarse woody debris pool
  integer, private :: i_soil1 = -9  ! Soil Organic Matter (SOM) first pool, protected SOM (som_p in Wieder et al. 2015)
  integer, private :: i_soil2 = -9  ! SOM second pool, recalcitrant SOM (som_c in Wieder et al. 2015)
  integer, private :: i_soil3 = -9  ! SOM third pool, available SOM (som_a in Wieder et al. 2015)

  type, private :: params_type
     ! TODO slevis: soil/litter C:N is not defined in MIMICS. Microbial C:N is
     !              a parameter that varies with litter quality (fmet).

     ! TODO slevis: New rf params in MIMICS. Update params file accordingly.
     !              I haven't identified what these correspond to in testbed.
     real(r8):: rf_l1m1_mimics  ! respiration fraction litter 1 -> microbe 1
     real(r8):: rf_l2m1_mimics
     real(r8):: rf_s3m1_mimics  ! respiration fraction soil 3 -> microbe 1
     real(r8):: rf_l1m2_mimics
     real(r8):: rf_l2m2_mimics
     real(r8):: rf_s3m2_mimics

     real(r8):: rf_cwdl2_mimics

     ! TODO slevis: Introduce Vmax and Km. Also tau params that define
     !              microbial turnover to SOM are functions.
     real(r8):: tau_m1_mimics     ! 1/turnover time of microbial 1 (1/yr)
     real(r8):: tau_m2_mimics     ! 1/turnover time of microbial 2 (1/yr)
     real(r8):: tau_cwd_bgc   ! corrected fragmentation rate constant CWD, century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1 (1/0.3) (1/yr)

     real(r8) :: cwd_flig_bgc !

     real(r8) :: k_frag_bgc   !fragmentation rate for CWD
     real(r8) :: minpsi_bgc   !minimum soil water potential for heterotrophic resp
     real(r8) :: maxpsi_bgc   !maximum soil water potential for heterotrophic resp

     real(r8) :: initial_Cstocks(nsompools) ! Initial Carbon stocks for a cold-start
     real(r8) :: initial_Cstocks_depth      ! Soil depth for initial Carbon stocks for a cold-start
     
  end type params_type
  !
  type(params_type), private :: params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine DecompCascadeMIMCSreadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for soil MIMICS Decomposition Cascade
    ! TODO slevis: Modify for MIMICS if necessary
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use abortutils     , only : endrun
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'DecompCascadeMIMICSreadNML'
    character(len=*), parameter :: nmlname = 'CENTURY_soilBGCDecompCascade'
    !-----------------------------------------------------------------------
    real(r8) :: initial_Cstocks(nsompools), initial_Cstocks_depth
    namelist /CENTURY_soilBGCDecompCascade/ initial_Cstocks, initial_Cstocks_depth

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    initial_Cstocks(:)    = 200._r8
    initial_Cstocks_depth = 0.3

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=CENTURY_soilBGCDecompCascade, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (initial_Cstocks      , mpicom)
    call shr_mpi_bcast (initial_Cstocks_depth, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=CENTURY_soilBGCDecompCascade)
       write(iulog,*) ' '
    end if

    params_inst%initial_Cstocks(:)    = initial_Cstocks(:)
    params_inst%initial_Cstocks_depth = initial_Cstocks_depth

  end subroutine DecompCascadeMIMICSreadNML

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
    ! TODO slevis: Add new params here and in the params file
    tString='tau_cwd'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_cwd_bgc=tempr

    tString='rf_l1m1_mimics'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l1m1_mimics=tempr

    tString='rf_l2m1_mimics'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l2m1_mimics=tempr

    tString='rf_s3m1_mimics'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s3m1_mimics=tempr

    tString='rf_l1m2_mimics'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l1m2_mimics=tempr

    tString='rf_l2m2_mimics'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l2m2_mimics=tempr

    tString='rf_s3m2_mimics'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s3m2_mimics=tempr

    tString='rf_cwdl2_bgc'  ! MIMICS value of this same as BGC value (= 0)
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl2_mimics=tempr

    tString='k_frag'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_frag_bgc=tempr

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
    
  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_mimcs(bounds, soilbiogeochem_state_inst, soilstate_inst )
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
    ! TODO slevis: Keeping temporarily as a template
    !              for new params not coming from the params file.
    !              Also keeping other params that I'm not sure about, yet.
    real(r8), allocatable :: rf_s1s3(:,:)
    real(r8) :: rf_cwdl2
    real(r8) :: cwd_flig
    real(r8), allocatable :: f_s1s3(:,:)
    real(r8) :: f_s2s3

    integer :: i_l1m1
    integer :: i_l2m1
    integer :: i_s3m1
    integer :: i_l1m2
    integer :: i_l2m2
    integer :: i_s3m2
    integer :: i_cwdl2
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

      !------- time-constant coefficients ---------- !
      ! set respiration fractions for fluxes between compartments
      rf_l1m1 = params_inst%rf_l1m1_mimics
      rf_l2m1 = params_inst%rf_l2m1_mimics
      rf_s3m1 = params_inst%rf_s3m1_mimics
      rf_l1m2 = params_inst%rf_l1m2_mimics
      rf_l2m2 = params_inst%rf_l2m2_mimics
      rf_s3m2 = params_inst%rf_s3m2_mimics

      rf_cwdl2 = params_inst%rf_cwdl2_mimics

      ! set the lignin fractions for coarse woody debris
      cwd_flig = params_inst%cwd_flig_bgc

      ! set path fractions
      f_s2s3 = 0.03_r8/(0.45_r8)

      ! some of these are dependent on the soil texture properties
      ! TODO slevis: Template for calculated mimics params?
      !              But soil texture remains const with time.
      !              So one-time initializations here and time-dep params
      !              in subr. decomp_rate_constants_mimics.

      do c = bounds%begc, bounds%endc
         do j = 1, nlevdecomp
            t = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - cellsand(c,j))
            f_s1s3(c,j) = .004_r8 / (1._r8 - t)
            rf_s1s3(c,j) = t
         end do
      end do
      initial_stock_soildepth = params_inst%initial_Cstocks_depth

      !-------------------  list of pools and their attributes  ------------
      floating_cn_ratio_decomp_pools(i_litr1) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_litr1) = 'litr1'
      decomp_cascade_con%decomp_pool_name_history(i_litr1) = 'LITR1'
      decomp_cascade_con%decomp_pool_name_long(i_litr1) = 'litter 1'
      decomp_cascade_con%decomp_pool_name_short(i_litr1) = 'L1'
      is_microbe(i_litr1) = .false.
      is_litter(i_litr1) = .true.
      is_soil(i_litr1) = .false.
      is_cwd(i_litr1) = .false.
      initial_cn_ratio(i_litr1) = 90._r8
      initial_stock(i_litr1) = 0._r8
      is_metabolic(i_litr1) = .true.
      is_cellulose(i_litr1) = .false.
      is_lignin(i_litr1) = .false.

      i_litr2 = i_litr1 + 1
      floating_cn_ratio_decomp_pools(i_litr2) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_litr2) = 'litr2'
      decomp_cascade_con%decomp_pool_name_history(i_litr2) = 'LITR2'
      decomp_cascade_con%decomp_pool_name_long(i_litr2) = 'litter 2'
      decomp_cascade_con%decomp_pool_name_short(i_litr2) = 'L2'
      is_microbe(i_litr2) = .false.
      is_litter(i_litr2) = .true.
      is_soil(i_litr2) = .false.
      is_cwd(i_litr2) = .false.
      initial_cn_ratio(i_litr2) = 90._r8
      initial_stock(i_litr2) = 0._r8
      is_metabolic(i_litr2) = .false.
      is_cellulose(i_litr2) = .true.  ! TODO ??
      is_lignin(i_litr2) = .true.

      ! TODO slevis: clm_varpar hardwires i_cwd = 4 or = 0 if (use_fates)
      !              I recommend replacing the hardwired value with one that we
      !              update here on the fly. Same for i_soil1,2,3 & i_micr1,2.
      i_cwd = i_litr2
      if (.not. use_fates) then
         ! CWD
         i_cwd = i_litr2 + 1
         floating_cn_ratio_decomp_pools(i_cwd) = .true.
         decomp_cascade_con%decomp_pool_name_restart(i_cwd) = 'cwd'
         decomp_cascade_con%decomp_pool_name_history(i_cwd) = 'CWD'
         decomp_cascade_con%decomp_pool_name_long(i_cwd) = 'coarse woody debris'
         decomp_cascade_con%decomp_pool_name_short(i_cwd) = 'CWD'
         is_microbe(i_cwd) = .false.
         is_litter(i_cwd) = .false.
         is_soil(i_cwd) = .false.
         is_cwd(i_cwd) = .true.
         initial_cn_ratio(i_cwd) = 90._r8
         initial_stock(i_cwd) = 0._r8
         is_metabolic(i_cwd) = .false.
         is_cellulose(i_cwd) = .false.
         is_lignin(i_cwd) = .false.
      endif

      i_soil1 = i_cwd + 1
      floating_cn_ratio_decomp_pools(i_soil1) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_soil1) = 'soil1'
      decomp_cascade_con%decomp_pool_name_history(i_soil1) = 'SOIL1'
      decomp_cascade_con%decomp_pool_name_long(i_soil1) = 'soil 1'
      decomp_cascade_con%decomp_pool_name_short(i_soil1) = 'S1'
      is_microbe(i_soil1) = .false.
      is_litter(i_soil1) = .false.
      is_soil(i_soil1) = .true.
      is_cwd(i_soil1) = .false.
      initial_cn_ratio(i_soil1) = cn_s1  ! TODO ??
      initial_stock(i_soil1) = params_inst%initial_Cstocks(1)
      is_metabolic(i_soil1) = .false.
      is_cellulose(i_soil1) = .false.
      is_lignin(i_soil1) = .false.

      i_soil2 = i_soil1 + 1
      floating_cn_ratio_decomp_pools(i_soil2) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_soil2) = 'soil2'
      decomp_cascade_con%decomp_pool_name_history(i_soil2) = 'SOIL2'
      decomp_cascade_con%decomp_pool_name_long(i_soil2) = 'soil 2'
      decomp_cascade_con%decomp_pool_name_short(i_soil2) = 'S2'
      is_microbe(i_soil2) = .false.
      is_litter(i_soil2) = .false.
      is_soil(i_soil2) = .true.
      is_cwd(i_soil2) = .false.
      initial_cn_ratio(i_soil2) = cn_s2  ! TODO ??
      initial_stock(i_soil2) = params_inst%initial_Cstocks(2)
      is_metabolic(i_soil2) = .false.
      is_cellulose(i_soil2) = .false.
      is_lignin(i_soil2) = .false.

      i_soil3 = i_soil2 + 1
      floating_cn_ratio_decomp_pools(i_soil3) = .false.
      decomp_cascade_con%decomp_pool_name_restart(i_soil3) = 'soil3'
      decomp_cascade_con%decomp_pool_name_history(i_soil3) = 'SOIL3'
      decomp_cascade_con%decomp_pool_name_long(i_soil3) = 'soil 3'
      decomp_cascade_con%decomp_pool_name_short(i_soil3) = 'S3'
      is_microbe(i_soil3) = .false.
      is_litter(i_soil3) = .false.
      is_soil(i_soil3) = .true.
      is_cwd(i_soil3) = .false.
      initial_cn_ratio(i_soil3) = cn_s3  ! TODO ??
      initial_stock(i_soil3) = params_inst%initial_Cstocks(3)
      is_metabolic(i_soil3) = .false.
      is_cellulose(i_soil3) = .false.
      is_lignin(i_soil3) = .false.

      i_micr1 = i_soil3 + 1
      floating_cn_ratio_decomp_pools(i_micr1) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_micr1) = 'micr1'
      decomp_cascade_con%decomp_pool_name_history(i_micr1) = 'MICR1'
      decomp_cascade_con%decomp_pool_name_long(i_micr1) = 'microbial 1'
      decomp_cascade_con%decomp_pool_name_short(i_micr1) = 'M1'
      is_microbe(i_micr1) = .true.
      is_litter(i_micr1) = .false.
      is_soil(i_micr1) = .false.
      is_cwd(i_micr1) = .false.
      initial_cn_ratio(i_micr1) = 90._r8
      initial_stock(i_micr1) = 0._r8
      is_metabolic(i_micr1) = .false.
      is_cellulose(i_micr1) = .false.
      is_lignin(i_micr1) = .false.

      i_micr2 = i_micr1 + 1
      floating_cn_ratio_decomp_pools(i_micr2) = .true.
      decomp_cascade_con%decomp_pool_name_restart(i_micr2) = 'micr2'
      decomp_cascade_con%decomp_pool_name_history(i_micr2) = 'MICR2'
      decomp_cascade_con%decomp_pool_name_long(i_micr2) = 'microbial 2'
      decomp_cascade_con%decomp_pool_name_short(i_micr2) = 'M2'
      is_microbe(i_micr2) = .true.
      is_litter(i_micr2) = .false.
      is_soil(i_micr2) = .false.
      is_cwd(i_micr2) = .false.
      initial_cn_ratio(i_micr2) = 90._r8
      initial_stock(i_micr2) = 0._r8
      is_metabolic(i_micr2) = .false.
      is_cellulose(i_micr2) = .false.
      is_lignin(i_micr2) = .false.

      speedup_fac = 1._r8

      !lit1
      spinup_factor(i_litr1) = 1._r8
      !lit2
      spinup_factor(i_litr2) = 1._r8
      !CWD
      if (.not. use_fates) then
         spinup_factor(i_cwd) = max(1._r8, (speedup_fac * params_inst%tau_cwd_bgc / 2._r8 ))
      end if
      !som1
      spinup_factor(i_soil1) = 1._r8
      !som2,3
      ! TODO slevis: change these to 1?
      spinup_factor(i_soil2) = max(1._r8, (speedup_fac * params_inst%tau_s2_bgc))
      spinup_factor(i_soil3) = max(1._r8, (speedup_fac * params_inst%tau_s3_bgc))
      ! micr1,2
      spinup_factor(i_micr1) = 1._r8
      spinup_factor(i_micr2) = 1._r8

      if ( masterproc ) then
         write(iulog,*) 'Spinup_state ',spinup_state
         write(iulog,*) 'Spinup factors ',spinup_factor
      end if

      !----------------  list of transitions and their time-independent coefficients  ---------------!
      i_l1m1 = 1
      decomp_cascade_con%cascade_step_name(i_l1m1) = 'L1M1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m1) = rf_l1m1
      cascade_donor_pool(i_l1m1) = i_litr1
      cascade_receiver_pool(i_l1m1) = i_micr1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m1) = 1.0_r8

      i_l2m1 = 2
      decomp_cascade_con%cascade_step_name(i_l2m1) = 'L2M1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m1) = rf_l2m1
      cascade_donor_pool(i_l2m1) = i_litr2
      cascade_receiver_pool(i_l2m1) = i_micr1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m1)= 1.0_r8

      i_s3m1 = 3
      decomp_cascade_con%cascade_step_name(i_s3m1) = 'S3M1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3m1) = rf_s3m1
      cascade_donor_pool(i_s3m1) = i_soil3
      cascade_receiver_pool(i_s3m1) = i_micr1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3m1) = 1.0_r8

      i_l1m2 = 4
      decomp_cascade_con%cascade_step_name(i_l1m2) = 'L1M2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m2) = rf_l1m2
      cascade_donor_pool(i_l1m2) = i_litr1
      cascade_receiver_pool(i_l1m2) = i_micr2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m2) = 1.0_r8

      i_l2m2 = 5
      decomp_cascade_con%cascade_step_name(i_l2m2) = 'L2M2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m2) = rf_l2m2
      cascade_donor_pool(i_l2m2) = i_litr2
      cascade_receiver_pool(i_l2m2) = i_micr2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m2)= 1.0_r8

      i_s3m2 = 6
      decomp_cascade_con%cascade_step_name(i_s3m2) = 'S3M2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3m2) = rf_s3m2
      cascade_donor_pool(i_s3m2) = i_soil3
      cascade_receiver_pool(i_s3m2) = i_micr2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3m2) = 1.0_r8

      if (.not. use_fates) then
         i_cwdl2 = 7
         decomp_cascade_con%cascade_step_name(i_cwdl2) = 'CWDL2'
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
         cascade_donor_pool(i_cwdl2) = i_cwd
         cascade_receiver_pool(i_cwdl2) = i_litr2
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = cwd_flig
      end if

      deallocate(rf_s1s3)
      deallocate(f_s1s3)

    end associate

  end subroutine init_decompcascade_mimics

  !-----------------------------------------------------------------------
  subroutine decomp_rate_constants_mimics(bounds, num_soilc, filter_soilc, &
       soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! Calculate rate constants and decomposition pathways for the MIMICS
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
    real(r8):: k_l1                         ! decomposition rate constant litter 1 (1/sec)
    real(r8):: k_l2                         ! decomposition rate constant litter 2 (1/sec)
    real(r8):: k_s1                         ! decomposition rate constant SOM 1 (1/sec)
    real(r8):: k_s2                         ! decomposition rate constant SOM 2 (1/sec)
    real(r8):: k_s3                         ! decomposition rate constant SOM 3 (1/sec)
    real(r8):: k_m1                         ! decomposition rate constant microbial 1 (1/sec)
    real(r8):: k_m2                         ! decomposition rate constant microbial 2 (1/sec)
    real(r8):: k_frag                       ! fragmentation rate constant CWD (1/sec)
    real(r8):: cwdc_loss                    ! fragmentation rate for CWD carbon (gC/m2/s)
    real(r8):: cwdn_loss                    ! fragmentation rate for CWD nitrogen (gN/m2/s)
    real(r8):: decomp_depth_efolding        ! (meters) e-folding depth for reduction in decomposition [
    integer :: c, fc, j, k, l
    real(r8):: days_per_year                ! days per year
    real(r8):: depth_scalar(bounds%begc:bounds%endc,1:nlevdecomp) 
    real(r8):: mino2lim                     !minimum anaerobic decomposition rate
    real(r8):: spinup_geogterm_l1(bounds%begc:bounds%endc) ! geographically-varying spinup term for l1
    real(r8):: spinup_geogterm_l2(bounds%begc:bounds%endc) ! geographically-varying spinup term for l2
    real(r8):: spinup_geogterm_cwd(bounds%begc:bounds%endc) ! geographically-varying spinup term for cwd
    real(r8):: spinup_geogterm_s1(bounds%begc:bounds%endc) ! geographically-varying spinup term for s1
    real(r8):: spinup_geogterm_s2(bounds%begc:bounds%endc) ! geographically-varying spinup term for s2
    real(r8):: spinup_geogterm_s3(bounds%begc:bounds%endc) ! geographically-varying spinup term for s3
    real(r8):: spinup_geogterm_m1(bounds%begc:bounds%endc) ! geographically-varying spinup term for m1
    real(r8):: spinup_geogterm_m2(bounds%begc:bounds%endc) ! geographically-varying spinup term for m2

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
         decomp_k       => soilbiogeochem_carbonflux_inst%decomp_k_col , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         spinup_factor  => decomp_cascade_con%spinup_factor              & ! Input:  [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
         )

      mino2lim = CNParamsShareInst%mino2lim

      days_per_year = get_days_per_year()

      ! Set "decomp_depth_efolding" parameter
      decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding

! TODO slevis: tau params are changing; map to Wieder et al. 2015
! From testbed code
! ------------------

! Vmax, tau*, and Km are time-dependent

! Table B1 Wieder et al. 2015 gives
! Vslope = 0.063, Vint = 5.47, av = 8e-6
! Vmod_r = 10,2,10, Vmod_K = 3,3,2 (for LITm,s,SOMa --> MICr,K respectively)

! Vmax(npt,R1) = exp(Vslope(R1) * Tsoi_day_avg_C + Vint(R1)) * av * Vmod(R1)
! Vmax(npt,R2) = exp(Vslope(R2) * Tsoi_day_avg_C + Vint(R2)) * av * Vmod(R2)
! Vmax(npt,R3) = exp(Vslope(R3) * Tsoi_day_avg_C + Vint(R3)) * av * Vmod(R3)
! Vmax(npt,K1) = exp(Vslope(K1) * Tsoi_day_avg_C + Vint(K1)) * av * Vmod(K1)
! Vmax(npt,K2) = exp(Vslope(K2) * Tsoi_day_avg_C + Vint(K2)) * av * Vmod(K2)
! Vmax(npt,K3) = exp(Vslope(K3) * Tsoi_day_avg_C + Vint(K3)) * av * Vmod(K3)

! ---
! Solving for tauR (tau_m1 here?) & tauK (tau_m2 here?)

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

! tauR is tau_m1, right? tau_r(1) = 5.2e-4_r8, tau_r(2) = 0.3_r8 or 0.4?
! tauR(npt) = tau_r(1) * exp(tau_r(2) * fmet(npt)) * tauMod(npt)

! tauK is tau_m2, right? tau_k(1) = 2.4e-4_r8, tau_k(2) = 0.1_r8
! tauK(npt) = tau_k(1) * exp(tau_k(2) * fmet(npt)) * tauMod(npt)
! ---

! Table B1 Wieder et al. 2015 gives
! Kslope = 0.017, 0.027, 0.017 (for LITm,s,SOMa --> MICr and K)
! Kint = 3.19, aK = 10
! Kmod_r = 0.125,0.5,0.25*P_scalar, Kmod_K = 0.5,0.25,0.167*P_scalar
! P_scalar = (2 * exp(-2 * sqrt(fclay)))**(-1)

! Km(npt,R1) = exp(Kslope(R1) * Tsoi_day_avg_C + Kint(R1)) * ak / Kmod(npt,R1)
! Km(npt,R2) = exp(Kslope(R2) * Tsoi_day_avg_C + Kint(R2)) * ak / Kmod(npt,R2)
! Km(npt,R3) = exp(Kslope(R3) * Tsoi_day_avg_C + Kint(R3)) * ak / Kmod(npt,R3)
! Km(npt,K1) = exp(Kslope(K1) * Tsoi_day_avg_C + Kint(K1)) * ak / Kmod(npt,K1)
! Km(npt,K2) = exp(Kslope(K2) * Tsoi_day_avg_C + Kint(K2)) * ak / Kmod(npt,K2)
! Km(npt,K3) = exp(Kslope(K3) * Tsoi_day_avg_C + Kint(K3)) * ak / Kmod(npt,K3)

! STILL MISSING N-related stuff: DIN...

      ! translate to per-second time constant
      k_l1 = 1._r8    / (secspday * days_per_year * params_inst%tau_l1_bgc)
      k_l2 = 1._r8    / (secspday * days_per_year * params_inst%tau_l2_bgc)
      k_s1 = 1._r8    / (secspday * days_per_year * params_inst%tau_s1_bgc)
      k_s2 = 1._r8    / (secspday * days_per_year * params_inst%tau_s2_bgc)
      k_s3 = 1._r8    / (secspday * days_per_year * params_inst%tau_s3_bgc)
      k_m1 = 1._r8    / (secspday * days_per_year * params_inst%tau_m1_mimics)
      k_m2 = 1._r8    / (secspday * days_per_year * params_inst%tau_m2_mimics)
      k_frag = 1._r8  / (secspday * days_per_year * params_inst%tau_cwd_bgc)

     ! calc ref rate
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
            if ( abs(spinup_factor(i_micr1) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_m1(c) = spinup_factor(i_micr1) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_m1(c) = 1._r8
            endif
            !
            if ( abs(spinup_factor(i_micr2) - 1._r8) .gt. .000001_r8) then
               spinup_geogterm_m2(c) = spinup_factor(i_micr2) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
               spinup_geogterm_m2(c) = 1._r8
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
            spinup_geogterm_m1(c) = 1._r8
            spinup_geogterm_m2(c) = 1._r8
         end do
      endif

      !--- time dependent coefficients-----!
      ! TODO slevis: Not immediately obvious how to reduce the if-else
      !      repetition bc w_scalar, o_scalar sum all j
      !      to layer 1 in the if, while they do not in the else.
      !      Other changes that would be necessary:
      !      - set fr = 1 in the else,
      !      - use a new nlevd = nlev_soildecomp_standard in the if and
      !                  nlevd = nlevdecomp in the else,
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
      ! TODO slevis: Reduced code repetition by defining
      !              depth_scalar(c,j) = 1.0_r8 in the else instead of
      !              repeating the decomp_k lines.
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

      ! calculate rate constants for all litter and som pools
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_k(c,j,i_litr1) = k_l1 * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l1(c)
            decomp_k(c,j,i_litr2) = k_l2 * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l2(c)

            decomp_k(c,j,i_micr1) = k_m1 * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_m1(c)
            decomp_k(c,j,i_micr2) = k_m2 * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_m2(c)

            decomp_k(c,j,i_soil1) = k_s1 * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s1(c)
            decomp_k(c,j,i_soil2) = k_s2 * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s2(c)
            decomp_k(c,j,i_soil3) = k_s3 * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s3(c)
            ! Same for cwd but only if fates not enabled; fates handles cwd on
            ! its own structure
            if (.not. use_fates) then
               decomp_k(c,j,i_cwd) = k_frag * w_scalar(c,j) * &
                  depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_cwd(c)
            end if
         end do
      end do

    end associate

 end subroutine decomp_rate_constants_mimics

end module SoilBiogeochemDecompCascadeMIMICSMod
