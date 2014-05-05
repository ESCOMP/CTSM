module EDInitTimeConstMod

  ! !DESCRIPTION:
  ! Initialize ED time invariant variables that interact with CLM

  implicit none
  save
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: EDInitTimeConst

contains

!-----------------------------------------------------------------------
subroutine EDInitTimeConst( bounds, EDpftcon, EDPftvarcon_inst, numpft )
  !
  ! !DESCRIPTION:
  ! Initialize time invariant clm variables
  ! 1) removed references to shallow lake - since it is not used
  ! 2) ***Make z, zi and dz allocatable depending on if you
  !    have lake or soil
  ! 3) rootfr only initialized for soil points
  !
  ! !USES:
  use EDPftvarcon          , only : EDPftvarcon_type
  use EDClmtype            , only : EDpft_epc_type
  use decompMod            , only : bounds_type

  type(bounds_type),       intent(in)    :: bounds  ! bounds
  type(EDpft_epc_type),    intent(inout) :: EDpftcon
  type(EDPftvarcon_type),  intent(in)    :: EDPftvarcon_inst
  integer,                 intent(in)    :: numpft

  integer :: m

    do m = 0,numpft
         !ED variables
         EDpftcon%max_dbh(m)               = EDPftvarcon_inst%max_dbh(m)
         EDpftcon%freezetol(m)             = EDPftvarcon_inst%freezetol(m)
         EDpftcon%wood_density(m)          = EDPftvarcon_inst%wood_density(m)
         EDpftcon%alpha_stem(m)            = EDPftvarcon_inst%alpha_stem(m)
         EDpftcon%hgt_min(m)               = EDPftvarcon_inst%hgt_min(m)
         EDpftcon%cushion(m)               = EDPftvarcon_inst%cushion(m)
         EDpftcon%leaf_stor_priority(m)    = EDPftvarcon_inst%leaf_stor_priority(m)
         EDpftcon%leafwatermax(m)          = EDPftvarcon_inst%leafwatermax(m)
         EDpftcon%rootresist(m)            = EDPftvarcon_inst%rootresist(m)
         EDpftcon%soilbeta(m)              = EDPftvarcon_inst%soilbeta(m)
         EDpftcon%crown(m)                 = EDPftvarcon_inst%crown(m)
         EDpftcon%bark_scaler(m)           = EDPftvarcon_inst%bark_scaler(m)
         EDpftcon%crown_kill(m)            = EDPftvarcon_inst%crown_kill(m)
         EDpftcon%initd(m)                 = EDPftvarcon_inst%initd(m)
         EDpftcon%sd_mort(m)               = EDPftvarcon_inst%sd_mort(m)
         EDpftcon%seed_rain(m)             = EDPftvarcon_inst%seed_rain(m)
         EDpftcon%bb_slope(m)              = EDPftvarcon_inst%bb_slope(m)
         EDpftcon%root_long(m)             = EDPftvarcon_inst%root_long(m)
         EDpftcon%seed_alloc(m)            = EDPftvarcon_inst%seed_alloc(m)
         EDpftcon%clone_alloc(m)           = EDPftvarcon_inst%clone_alloc(m)
         EDpftcon%sapwood_ratio(m)         = EDPftvarcon_inst%sapwood_ratio(m)
     end do

  end subroutine EDInitTimeConst

end module EDInitTimeConstMod
