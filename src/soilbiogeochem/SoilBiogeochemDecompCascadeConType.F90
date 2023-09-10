module SoilBiogeochemDecompCascadeConType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Decomposition Cascade Type
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use abortutils     , only : endrun
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : ndecomp_cascade_transitions, ndecomp_pools, nlevdecomp
  use clm_varcon     , only : ispval
  use SparseMatrixMultiplyMod, only : sparse_matrix_type, diag_matrix_type, vector_type
  use clm_varctl     , only : iulog
  !
  implicit none

  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: decomp_cascade_par_init          ! Initialize the parameters
  public :: init_decomp_cascade_constants    ! Initialize the constants
  public :: InitSoilTransfer                 ! Initialize soil transfer (for soil matrix)
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
     logical           , pointer  :: is_microbe(:)                     ! TRUE => pool is a microbe pool
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

     ! Matrix data

  end type decomp_cascade_type

  integer, public, parameter :: i_atm = 0                              ! for terminal pools (i.e. 100% respiration) (only used for CN not for BGC)
  integer, public, parameter :: no_soil_decomp = 0                     ! No soil decomposition is done
  integer, public, parameter :: century_decomp = 1                     ! CENTURY decomposition method type
  integer, public, parameter :: mimics_decomp = 2                      ! MIMICS decomposition method type
  integer, public            :: decomp_method  = ispval                ! Type of decomposition to use
  logical, public, parameter :: use_soil_matrixcn = .false.            ! true => use cn matrix solution for soil BGC
  type(decomp_cascade_type), public :: decomp_cascade_con
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine decomp_cascade_par_init( NLFilename )
    use clm_varctl         , only : use_cn, use_fates_bgc
    use clm_varpar         , only : ndecomp_pools_max
    use spmdMod            , only : masterproc, mpicom
    use clm_nlUtilsMod     , only : find_nlgroup_name
    use shr_mpi_mod        , only : shr_mpi_bcast
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    ! !LOGAL VARIABLES:
    integer           :: nu_nml              ! unit for namelist file
    integer           :: nml_error           ! namelist i/o error flag
    character(len=20) :: soil_decomp_method  ! String description of soil decomposition method to use
    character(len=*), parameter :: nmlname = 'soilbgc_decomp'
    namelist/soilbgc_decomp/ soil_decomp_method

    ! Default value
    soil_decomp_method = 'UNSET'
    ! Read in soil BGC decomposition namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, nmlname, status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=soilbgc_decomp,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun('trouble reading '//nmlname//' namelist')
          end if
       else
          call endrun('trouble finding '//nmlname//' namelist')
       end if
       close(nu_nml)
       select case( soil_decomp_method )
       case( 'None' ) 
          decomp_method = no_soil_decomp
       case( 'CENTURYKoven2013' ) 
          decomp_method = century_decomp
       case( 'MIMICSWieder2015' )
          decomp_method = mimics_decomp
       case default
          call endrun('Bad soil_decomp_method = '//soil_decomp_method )
       end select
    endif
    ! Broadcast namelist items to all processors
    call shr_mpi_bcast(decomp_method, mpicom)
    ! Don't do anything if neither FATES or BGC is on
    if ( use_cn ) then
       if ( decomp_method == no_soil_decomp )then
          call endrun('When running with BGC an active soil_decomp_method must be used')
       end if
    else if ( use_fates_bgc ) then
       if ( decomp_method == no_soil_decomp )then
          call endrun('When running with FATES and without FATES-SP, an active soil_decomp_method must be used')
       end if
    else
       if ( decomp_method /= no_soil_decomp )then
          call endrun('When running without FATES or BGC soil_decomp_method must be None')
       end if
    end if
     
    if ( decomp_method /= no_soil_decomp )then
       ! We hardwire these parameters here because we use them
       ! in InitAllocate (in SoilBiogeochemStateType) which is called earlier than
       ! init_decompcascade_bgc where they might have otherwise been derived on the
       ! fly. For reference, if they were determined in init_decompcascade_bgc:
       ! ndecomp_pools would get the value of i_pas_som or i_cwd and
       ! ndecomp_cascade_transitions would get the value of i_s3s1 or i_cwdl3
       ! depending on how use_fates is set.
       if ( use_fates_bgc ) then
          if (decomp_method == century_decomp) then
             ndecomp_pools = 6
             ndecomp_cascade_transitions = 8
          else if (decomp_method == mimics_decomp) then
             ndecomp_pools = 7
             ndecomp_cascade_transitions = 14
          end if
       else
          if (decomp_method == century_decomp) then
             ndecomp_pools = 7
             ndecomp_cascade_transitions = 10
          else if (decomp_method == mimics_decomp) then
             ndecomp_pools = 8
             ndecomp_cascade_transitions = 15
          end if
       endif
       ! The next param also appears as a dimension in the params files dated
       ! c210418.nc and later
       ndecomp_pools_max = 8  ! largest ndecomp_pools value above
    else
       ndecomp_pools               = 7
       ndecomp_cascade_transitions = 7
       ndecomp_pools_max           = 8
    end if
    ! Set ndecomp_pools_vr needed for Matrix solution

  end subroutine decomp_cascade_par_init

  !------------------------------------------------------------------------
  subroutine init_decomp_cascade_constants( )
    !
    ! !DESCRIPTION:
    ! Initialize decomposition cascade state
    !------------------------------------------------------------------------
    ! !LOGAL VARIABLES:
    integer           :: ibeg                ! Beginning index for allocated arrays

    !
    ! NOTE: Return early if soil decomposition isn't being done
    !
    if ( decomp_method == ispval ) then
       call endrun('soil_decomp_method was not set, but soil decomposition cascade initialization is being called' )
    else if ( decomp_method /= no_soil_decomp ) then

       ibeg = 1

       !-- properties of each pathway along decomposition cascade 
       allocate(decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions))
       allocate(decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions))
       allocate(decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions))
   
       !-- properties of each decomposing pool
       allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%decomp_pool_name_restart(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%decomp_pool_name_history(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%decomp_pool_name_long(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%decomp_pool_name_short(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%is_microbe(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%is_litter(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%is_soil(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%is_cwd(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%initial_cn_ratio(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%initial_stock(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%is_metabolic(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%is_cellulose(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%is_lignin(ibeg:ndecomp_pools))
       allocate(decomp_cascade_con%spinup_factor(1:ndecomp_pools))

       ! Allocate soil matrix data
       if(use_soil_matrixcn)then
       end if
   
       !-- properties of each pathway along decomposition cascade 
       decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions) = ''
       decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions) = 0
       decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions) = 0
   
       !-- first initialization of properties of each decomposing pool
       decomp_cascade_con%floating_cn_ratio_decomp_pools(ibeg:ndecomp_pools) = .false.
       decomp_cascade_con%decomp_pool_name_history(ibeg:ndecomp_pools)       = ''
       decomp_cascade_con%decomp_pool_name_restart(ibeg:ndecomp_pools)       = ''
       decomp_cascade_con%decomp_pool_name_long(ibeg:ndecomp_pools)          = ''
       decomp_cascade_con%decomp_pool_name_short(ibeg:ndecomp_pools)         = ''
       decomp_cascade_con%is_microbe(ibeg:ndecomp_pools)                     = .false.
       decomp_cascade_con%is_litter(ibeg:ndecomp_pools)                      = .false.
       decomp_cascade_con%is_soil(ibeg:ndecomp_pools)                        = .false.
       decomp_cascade_con%is_cwd(ibeg:ndecomp_pools)                         = .false.
       decomp_cascade_con%initial_cn_ratio(ibeg:ndecomp_pools)               = nan
       decomp_cascade_con%initial_stock(ibeg:ndecomp_pools)                  = nan
       decomp_cascade_con%initial_stock_soildepth                            = 0.3
       decomp_cascade_con%is_metabolic(ibeg:ndecomp_pools)                   = .false.
       decomp_cascade_con%is_cellulose(ibeg:ndecomp_pools)                   = .false.
       decomp_cascade_con%is_lignin(ibeg:ndecomp_pools)                      = .false.
       decomp_cascade_con%spinup_factor(1:ndecomp_pools)                     = nan
    end if

    ! Soil matrix sizes
    if(use_soil_matrixcn)then
    else
       ! Set to missing value if not used
    end if
  end subroutine init_decomp_cascade_constants

  !------------------------------------------------------------------------
  subroutine InitSoilTransfer()
    !
    ! !DESCRIPTION:
    ! Initialize sparse matrix variables and index. Count possible non-zero entries and record their x and y in the matrix.
    ! Collect those non-zero entry information, and save them into the list.
    !------------------------------------------------------------------------
    ! !USES:
    use SparseMatrixMultiplyMod, only : sparse_matrix_type, diag_matrix_type, vector_type
    ! !LOGAL VARIABLES:
    integer i,j,k,m,n
  
  end subroutine InitSoilTransfer

end module SoilBiogeochemDecompCascadeConType
