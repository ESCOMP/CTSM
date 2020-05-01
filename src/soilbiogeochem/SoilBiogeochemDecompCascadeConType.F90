module SoilBiogeochemDecompCascadeConType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Decomposition Cascade Type
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : ndecomp_cascade_transitions, ndecomp_pools
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_decomp_cascade_constants
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
  end type decomp_cascade_type

  integer, public, parameter :: i_atm = 0                              ! for terminal pools (i.e. 100% respiration) (only used for CN not for BGC)
  type(decomp_cascade_type), public :: decomp_cascade_con
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_decomp_cascade_constants( use_century_decomp )
    !
    ! !DESCRIPTION:
    ! Initialize decomposition cascade state
    !------------------------------------------------------------------------
    ! !ARGUMENTS:
    logical, intent(IN) :: use_century_decomp
    ! !LOGAL VARIABLES:
    integer :: ibeg    ! Beginning index for allocated arrays

    if ( use_century_decomp ) then
       ibeg = 1
    else
       ibeg = i_atm
    end if
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
    allocate(decomp_cascade_con%is_litter(ibeg:ndecomp_pools))
    allocate(decomp_cascade_con%is_soil(ibeg:ndecomp_pools))
    allocate(decomp_cascade_con%is_cwd(ibeg:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cn_ratio(ibeg:ndecomp_pools))
    allocate(decomp_cascade_con%initial_stock(ibeg:ndecomp_pools))
    allocate(decomp_cascade_con%is_metabolic(ibeg:ndecomp_pools))
    allocate(decomp_cascade_con%is_cellulose(ibeg:ndecomp_pools))
    allocate(decomp_cascade_con%is_lignin(ibeg:ndecomp_pools))
    allocate(decomp_cascade_con%spinup_factor(1:ndecomp_pools))

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

  end subroutine init_decomp_cascade_constants

end module SoilBiogeochemDecompCascadeConType
