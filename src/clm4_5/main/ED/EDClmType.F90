module EDClmtype

!----------------------------------------------------------------------- 
! !DESCRIPTION: 
! Define ED derived type hierarchy as it relates to CLM type.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
!                              

!*******************************************************************************
!----------------------------------------------------
! pft ecophysiological constants structure
!----------------------------------------------------
type, public :: EDpft_epc_type
                                                !ED variables. 
   real(r8), pointer :: max_dbh(:)              !maximum dbh at which height growth ceases... 
   real(r8), pointer :: freezetol(:)            !minimum temperature tolerance... 
   real(r8), pointer :: wood_density(:)         ! wood density  g cm^-3  ...  
   real(r8), pointer :: alpha_stem(:)           ! live stem turnover rate. y-1 
   real(r8), pointer :: hgt_min(:)              ! sapling height m 
   real(r8), pointer :: cushion(:)              ! labile carbon storage target as multiple of leaf pool. 
   real(r8), pointer :: leaf_stor_priority(:)   ! leaf turnover vs labile carbon use prioritisation. 
                                                ! (1 = lose  leaves, 0 = use store). 
   real(r8), pointer :: leafwatermax(:)         ! amount of water allowed on leaf   surfaces
   real(r8), pointer :: rootresist(:)
   real(r8), pointer :: soilbeta(:)
   real(r8), pointer :: crown(:)                !fraction of the height of the plant that is occupied by crown. For fire model. 
   real(r8), pointer :: bark_scaler(:)          !scaler from dbh to bark thickness. For fire model. 
   real(r8), pointer :: crown_kill(:)           !scaler on fire death. For fire model. 
   real(r8), pointer :: initd(:)                !initial seedling density 
   real(r8), pointer :: sd_mort(:)              !rate of death of seeds produced from reproduction. 
   real(r8), pointer :: seed_rain(:)            !seeds that come from outside the gridbox.  
   real(r8), pointer :: BB_slope(:)             !ball berry slope parameter
   real(r8), pointer :: root_long(:)            !root longevity (yrs)
   real(r8), pointer :: clone_alloc(:)          !fraction of carbon balance allocated to clonal reproduction. 
   real(r8), pointer :: seed_alloc(:)           !fraction of carbon balance allocated to seeds. 
   real(r8), pointer :: sapwood_ratio(:)        !amount of sapwood per unit leaf carbon and m height 
end type EDpft_epc_type

type(EDpft_epc_type), public, target, save :: EDpftcon

                                                !----------------------------------------------------
                                                ! pft carbon flux variables structure
                                                !----------------------------------------------------
type, public :: EDpft_cflux_type
                                                !new variables for ED code. 
   real(r8), pointer :: trimming(:) 
   real(r8), pointer :: area_plant(:) 
   real(r8), pointer :: area_trees(:) 
   real(r8), pointer :: canopy_spread(:) 
   real(r8), pointer :: GCcanopy(:)             ! mmol m-2 s-1
   real(r8), pointer :: PFTbiomass(:,:)         ! total biomass of each pft
   real(r8), pointer :: PFTleafbiomass(:,:)     ! total biomass of each pft   
   real(r8), pointer :: PFTstorebiomass(:,:)    ! total biomass of each pft   
   real(r8), pointer :: PFTnindivs(:,:)         ! total biomass of each pft 
  
   real(r8), pointer :: nesterov_fire_danger(:) ! total biomass of each pft 
   real(r8), pointer :: spitfire_ROS(:)         ! total biomass of each pft 
   real(r8), pointer :: effect_wspeed(:)              ! total biomass of each pft 
   real(r8), pointer :: TFC_ROS(:)              ! total biomass of each pft 
   real(r8), pointer :: fire_intensity(:)       ! total biomass of each pft 
   real(r8), pointer :: fire_area(:)            ! total biomass of each pft 
   real(r8), pointer :: scorch_height(:)        ! total biomass of each pft 
   real(r8), pointer :: fire_fuel_bulkd(:)      ! total biomass of each pft 
   real(r8), pointer :: fire_fuel_eff_moist(:)  ! total biomass of each pft 
   real(r8), pointer :: fire_fuel_sav(:)        ! total biomass of each pft       
   real(r8), pointer :: fire_fuel_mef(:)        ! total biomass of each pft 
   real(r8), pointer :: sum_fuel(:)             ! total biomass of each pft 

   real(r8), pointer :: litter_in(:)            ! total biomass of each pft 
   real(r8), pointer :: litter_out(:)           ! total biomass of each pft    
   real(r8), pointer :: efpot(:)                ! potential transpiration
   real(r8), pointer :: rb(:)                   ! boundary layer conductance
   
   real(r8), pointer :: daily_temp(:)           ! daily temperature for fire and phenology models
   real(r8), pointer :: daily_rh(:)             ! daily RH for fire model
   real(r8), pointer :: daily_prec(:)           ! daily rain for fire and phenology models. 
   real(r8), pointer :: prec24(:)               ! 24hr sum of precipitation 
   real(r8), pointer :: RH24(:)                 ! 24hr average of RH
   real(r8), pointer :: wind24(:)               ! 24hr average of wind. 
   
                                                !seed model. Aggregated to gridcell for now. 
   real(r8), pointer :: seed_bank(:)            ! Mass of seeds.                 kGC/m2
   real(r8), pointer :: seeds_in(:)             ! Production of seed mass.       kGC/m2/year
   real(r8), pointer :: seed_decay(:)           ! Decay of seed mass.            kGC/m2/year 
   real(r8), pointer :: seed_germination(:)     ! Germiantion rate of seed mass. kGC/m2/year
end type EDpft_cflux_type

type(EDpft_cflux_type), target :: EDpcf         ! ED pft carbon fluxes

!----------------------------------------------------
! define the ED pft structure
!----------------------------------------------------

type, public :: EDpft_type
   integer , pointer :: ED_patch(:)     ! Is there a patch of vegetation defined for this p?  
   integer , pointer :: ED_bareground(:)! Is this patch the designated bare ground fraction?  
   real(r8), pointer :: wtED(:)         ! What is the weight of each patch in the ED natural vegetation area? 
end type EDpft_type

type(EDpft_type), target :: EDpft  !plant functional type (pft) data structure 

!
!EOP
!----------------------------------------------------------------------- 
end module EDClmtype  
