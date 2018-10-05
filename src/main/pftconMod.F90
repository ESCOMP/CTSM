module pftconMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation constants and method to
  ! read and initialize vegetation (PFT) constants.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar  , only : mxpft, numrad, ivis, inir, cft_lb, cft_ub
  use clm_varctl  , only : iulog, use_cndv, use_vertsoilc, use_crop
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! Vegetation type constants
  !
  integer, public :: noveg                  ! value for not vegetated 
  integer, public :: ndllf_evr_tmp_tree     ! value for Needleleaf evergreen temperate tree
  integer, public :: ndllf_evr_brl_tree     ! value for Needleleaf evergreen boreal tree
  integer, public :: ndllf_dcd_brl_tree     ! value for Needleleaf deciduous boreal tree
  integer, public :: nbrdlf_evr_trp_tree    ! value for Broadleaf evergreen tropical tree
  integer, public :: nbrdlf_evr_tmp_tree    ! value for Broadleaf evergreen temperate tree
  integer, public :: nbrdlf_dcd_trp_tree    ! value for Broadleaf deciduous tropical tree
  integer, public :: nbrdlf_dcd_tmp_tree    ! value for Broadleaf deciduous temperate tree
  integer, public :: nbrdlf_dcd_brl_tree    ! value for Broadleaf deciduous boreal tree
  integer, public :: ntree                  ! value for last type of tree
  integer, public :: nbrdlf_evr_shrub       ! value for Broadleaf evergreen shrub
  integer, public :: nbrdlf_dcd_tmp_shrub   ! value for Broadleaf deciduous temperate shrub
  integer, public :: nbrdlf_dcd_brl_shrub   ! value for Broadleaf deciduous boreal shrub
  integer, public :: nc3_arctic_grass       ! value for C3 arctic grass
  integer, public :: nc3_nonarctic_grass    ! value for C3 non-arctic grass
  integer, public :: nc4_grass              ! value for C4 grass
  integer, public :: npcropmin              ! value for first crop
  integer, public :: ntmp_corn              ! value for temperate corn, rain fed (rf)
  integer, public :: nirrig_tmp_corn        ! value for temperate corn, irrigated (ir)
  integer, public :: nswheat                ! value for spring temperate cereal (rf)
  integer, public :: nirrig_swheat          ! value for spring temperate cereal (ir)
  integer, public :: nwwheat                ! value for winter temperate cereal (rf)
  integer, public :: nirrig_wwheat          ! value for winter temperate cereal (ir)
  integer, public :: ntmp_soybean           ! value for temperate soybean (rf)
  integer, public :: nirrig_tmp_soybean     ! value for temperate soybean (ir)
  integer, public :: nbarley                ! value for spring barley (rf)
  integer, public :: nirrig_barley          ! value for spring barley (ir)
  integer, public :: nwbarley               ! value for winter barley (rf)
  integer, public :: nirrig_wbarley         ! value for winter barley (ir)
  integer, public :: nrye                   ! value for spring rye (rf)
  integer, public :: nirrig_rye             ! value for spring rye (ir)
  integer, public :: nwrye                  ! value for winter rye (rf)
  integer, public :: nirrig_wrye            ! value for winter rye (ir)
  integer, public :: ncassava               ! ...and so on
  integer, public :: nirrig_cassava
  integer, public :: ncitrus
  integer, public :: nirrig_citrus
  integer, public :: ncocoa
  integer, public :: nirrig_cocoa
  integer, public :: ncoffee
  integer, public :: nirrig_coffee
  integer, public :: ncotton
  integer, public :: nirrig_cotton
  integer, public :: ndatepalm
  integer, public :: nirrig_datepalm
  integer, public :: nfoddergrass
  integer, public :: nirrig_foddergrass
  integer, public :: ngrapes
  integer, public :: nirrig_grapes
  integer, public :: ngroundnuts
  integer, public :: nirrig_groundnuts
  integer, public :: nmillet
  integer, public :: nirrig_millet
  integer, public :: noilpalm
  integer, public :: nirrig_oilpalm
  integer, public :: npotatoes
  integer, public :: nirrig_potatoes
  integer, public :: npulses
  integer, public :: nirrig_pulses
  integer, public :: nrapeseed
  integer, public :: nirrig_rapeseed
  integer, public :: nrice
  integer, public :: nirrig_rice
  integer, public :: nsorghum
  integer, public :: nirrig_sorghum
  integer, public :: nsugarbeet
  integer, public :: nirrig_sugarbeet
  integer, public :: nsugarcane
  integer, public :: nirrig_sugarcane
  integer, public :: nsunflower
  integer, public :: nirrig_sunflower
  integer, public :: nmiscanthus
  integer, public :: nirrig_miscanthus
  integer, public :: nswitchgrass
  integer, public :: nirrig_switchgrass
  integer, public :: ntrp_corn              !value for tropical corn (rf)
  integer, public :: nirrig_trp_corn        !value for tropical corn (ir)
  integer, public :: ntrp_soybean           !value for tropical soybean (rf)
  integer, public :: nirrig_trp_soybean     !value for tropical soybean (ir)
  integer, public :: npcropmax              ! value for last prognostic crop in list
  integer, public :: nc3crop                ! value for generic crop (rf)
  integer, public :: nc3irrig               ! value for irrigated generic crop (ir)

  ! Number of crop functional types actually used in the model. This includes each CFT for
  ! which is_pft_known_to_model is true. Note that this includes irrigated crops even if
  ! irrigation is turned off in this run: it just excludes crop types that aren't handled
  ! at all, as given by the mergetoclmpft list.
  integer, public :: num_cfts_known_to_model

  ! !PUBLIC TYPES:
  type, public :: pftcon_type

     integer , allocatable :: noveg         (:)   ! value for not vegetated
     integer , allocatable :: tree          (:)   ! tree or not?

     real(r8), allocatable :: dleaf         (:)   ! characteristic leaf dimension (m)
     real(r8), allocatable :: c3psn         (:)   ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8), allocatable :: xl            (:)   ! leaf/stem orientation index
     real(r8), allocatable :: rhol          (:,:) ! leaf reflectance: 1=vis, 2=nir
     real(r8), allocatable :: rhos          (:,:) ! stem reflectance: 1=vis, 2=nir
     real(r8), allocatable :: taul          (:,:) ! leaf transmittance: 1=vis, 2=nir
     real(r8), allocatable :: taus          (:,:) ! stem transmittance: 1=vis, 2=nir
     real(r8), allocatable :: z0mr          (:)   ! ratio of momentum roughness length to canopy top height (-)
     real(r8), allocatable :: displar       (:)   ! ratio of displacement height to canopy top height (-)
     real(r8), allocatable :: roota_par     (:)   ! CLM rooting distribution parameter [1/m]
     real(r8), allocatable :: rootb_par     (:)   ! CLM rooting distribution parameter [1/m]
     real(r8), allocatable :: crop          (:)   ! crop pft: 0. = not crop, 1. = crop pft
     real(r8), allocatable :: irrigated     (:)   ! irrigated pft: 0. = not, 1. = irrigated
     real(r8), allocatable :: smpso         (:)   ! soil water potential at full stomatal opening (mm)
     real(r8), allocatable :: smpsc         (:)   ! soil water potential at full stomatal closure (mm)
     real(r8), allocatable :: fnitr         (:)   ! foliage nitrogen limitation factor (-)

     !  CN code
     real(r8), allocatable :: dwood         (:)   ! wood density (gC/m3)
     real(r8), allocatable :: slatop        (:)   ! SLA at top of canopy [m^2/gC]
     real(r8), allocatable :: dsladlai      (:)   ! dSLA/dLAI [m^2/gC]
     real(r8), allocatable :: leafcn        (:)   ! leaf C:N [gC/gN]
     real(r8), allocatable :: flnr          (:)   ! fraction of leaf N in Rubisco [no units]
     real(r8), allocatable :: woody         (:)   ! woody lifeform flag (0 or 1)
     real(r8), allocatable :: lflitcn       (:)   ! leaf litter C:N (gC/gN)
     real(r8), allocatable :: frootcn       (:)   ! fine root C:N (gC/gN)
     real(r8), allocatable :: livewdcn      (:)   ! live wood (phloem and ray parenchyma) C:N (gC/gN)
     real(r8), allocatable :: deadwdcn      (:)   ! dead wood (xylem and heartwood) C:N (gC/gN)
     real(r8), allocatable :: grperc        (:)   ! growth respiration parameter
     real(r8), allocatable :: grpnow        (:)   ! growth respiration parameter
     real(r8), allocatable :: rootprof_beta (:,:) ! CLM rooting distribution parameter for C and N inputs [unitless]
     real(r8), allocatable :: root_radius   (:)   ! root radius (m)
     real(r8), allocatable :: root_density  (:)   ! root density (gC/m3)

     !  crop

     ! These arrays give information about the merge of unused crop types to the types CLM
     ! knows about. mergetoclmpft(m) gives the crop type that CLM uses to simulate input
     ! type m (and mergetoclmpft(m) == m implies that CLM simulates crop type m
     ! directly). is_pft_known_to_model(m) is true if CLM simulates crop type m, and false
     ! otherwise. Note that these do NOT relate to whether irrigation is on or off in a
     ! given simulation - that is handled separately.
     integer , allocatable :: mergetoclmpft         (:)
     logical , allocatable :: is_pft_known_to_model (:)

     real(r8), allocatable :: graincn       (:)   ! grain C:N (gC/gN)
     real(r8), allocatable :: mxtmp         (:)   ! parameter used in accFlds
     real(r8), allocatable :: baset         (:)   ! parameter used in accFlds
     real(r8), allocatable :: declfact      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: bfact         (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: aleaff        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: arootf        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: astemf        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: arooti        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: fleafi        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: allconsl      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: allconss      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: ztopmx        (:)   ! parameter used in CNVegStructUpdate
     real(r8), allocatable :: laimx         (:)   ! parameter used in CNVegStructUpdate
     real(r8), allocatable :: gddmin        (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: hybgdd        (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: lfemerg       (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: grnfill       (:)   ! parameter used in CNPhenology
     integer , allocatable :: mxmat         (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: mbbopt        (:)   ! Ball-Berry equation slope used in Photosynthesis
     real(r8), allocatable :: medlynslope   (:)   ! Medlyn equation slope used in Photosynthesis
     real(r8), allocatable :: medlynintercept(:)  ! Medlyn equation intercept used in Photosynthesis
     integer , allocatable :: mnNHplantdate (:)   ! minimum planting date for NorthHemisphere (YYYYMMDD)
     integer , allocatable :: mxNHplantdate (:)   ! maximum planting date for NorthHemisphere (YYYYMMDD)
     integer , allocatable :: mnSHplantdate (:)   ! minimum planting date for SouthHemisphere (YYYYMMDD)
     integer , allocatable :: mxSHplantdate (:)   ! maximum planting date for SouthHemisphere (YYYYMMDD)
     real(r8), allocatable :: planttemp     (:)   ! planting temperature used in CNPhenology (K)
     real(r8), allocatable :: minplanttemp  (:)   ! mininum planting temperature used in CNPhenology (K)
     real(r8), allocatable :: froot_leaf    (:)   ! allocation parameter: new fine root C per new leaf C (gC/gC) 
     real(r8), allocatable :: stem_leaf     (:)   ! allocation parameter: new stem c per new leaf C (gC/gC)
     real(r8), allocatable :: croot_stem    (:)   ! allocation parameter: new coarse root C per new stem C (gC/gC)
     real(r8), allocatable :: flivewd       (:)   ! allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
     real(r8), allocatable :: fcur          (:)   ! allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
     real(r8), allocatable :: fcurdv        (:)   ! alternate fcur for use with cndv
     real(r8), allocatable :: lf_flab       (:)   ! leaf litter labile fraction
     real(r8), allocatable :: lf_fcel       (:)   ! leaf litter cellulose fraction
     real(r8), allocatable :: lf_flig       (:)   ! leaf litter lignin fraction
     real(r8), allocatable :: fr_flab       (:)   ! fine root litter labile fraction
     real(r8), allocatable :: fr_fcel       (:)   ! fine root litter cellulose fraction
     real(r8), allocatable :: fr_flig       (:)   ! fine root litter lignin fraction
     real(r8), allocatable :: leaf_long     (:)   ! leaf longevity (yrs)
     real(r8), allocatable :: evergreen     (:)   ! binary flag for evergreen leaf habit (0 or 1)
     real(r8), allocatable :: stress_decid  (:)   ! binary flag for stress-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: season_decid  (:)   ! binary flag for seasonal-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: pconv         (:)   ! proportion of deadstem to conversion flux
     real(r8), allocatable :: pprod10       (:)   ! proportion of deadstem to 10-yr product pool
     real(r8), allocatable :: pprod100      (:)   ! proportion of deadstem to 100-yr product pool
     real(r8), allocatable :: pprodharv10   (:)   ! harvest mortality proportion of deadstem to 10-yr pool

     ! pft paraemeters for fire code
     real(r8), allocatable :: cc_leaf       (:)
     real(r8), allocatable :: cc_lstem      (:)
     real(r8), allocatable :: cc_dstem      (:)
     real(r8), allocatable :: cc_other      (:)
     real(r8), allocatable :: fm_leaf       (:)
     real(r8), allocatable :: fm_lstem      (:)
     real(r8), allocatable :: fm_dstem      (:)
     real(r8), allocatable :: fm_other      (:)
     real(r8), allocatable :: fm_root       (:)
     real(r8), allocatable :: fm_lroot      (:)
     real(r8), allocatable :: fm_droot      (:)
     real(r8), allocatable :: fsr_pft       (:)
     real(r8), allocatable :: fd_pft        (:)

     ! pft parameters for crop code
     real(r8), allocatable :: manunitro     (:)   ! manure
     real(r8), allocatable :: fleafcn       (:)   ! C:N during grain fill; leaf
     real(r8), allocatable :: ffrootcn      (:)   ! C:N during grain fill; fine root
     real(r8), allocatable :: fstemcn       (:)   ! C:N during grain fill; stem

     real(r8), allocatable :: i_vcad        (:)   
     real(r8), allocatable :: s_vcad        (:)   
     real(r8), allocatable :: i_flnr        (:)   
     real(r8), allocatable :: s_flnr        (:)     

     ! pft parameters for CNDV code (from LPJ subroutine pftparameters)
     real(r8), allocatable :: pftpar20      (:)   ! tree maximum crown area (m2)
     real(r8), allocatable :: pftpar28      (:)   ! min coldest monthly mean temperature
     real(r8), allocatable :: pftpar29      (:)   ! max coldest monthly mean temperature
     real(r8), allocatable :: pftpar30      (:)   ! min growing degree days (>= 5 deg C)
     real(r8), allocatable :: pftpar31      (:)   ! upper limit of temperature of the warmest month (twmax)
     
     ! pft parameters for FUN
     real(r8), allocatable :: a_fix         (:)   ! A BNF parameter
     real(r8), allocatable :: b_fix         (:)   ! A BNF parameter
     real(r8), allocatable :: c_fix         (:)   ! A BNF parameter
     real(r8), allocatable :: s_fix         (:)   ! A BNF parameter
     real(r8), allocatable :: akc_active    (:)   ! A mycorrhizal uptake parameter
     real(r8), allocatable :: akn_active    (:)   ! A mycorrhizal uptake parameter
     real(r8), allocatable :: ekc_active    (:)   ! A mycorrhizal uptake parameter
     real(r8), allocatable :: ekn_active    (:)   ! A mycorrhizal uptake parameter
     real(r8), allocatable :: kc_nonmyc     (:)   ! A non-mycorrhizal uptake parameter
     real(r8), allocatable :: kn_nonmyc     (:)   ! A non-mycorrhizal uptake parameter
     real(r8), allocatable :: kr_resorb     (:)   ! A retrasnlcation parameter
     real(r8), allocatable :: perecm        (:)   ! The fraction of ECM-associated PFT 
     real(r8), allocatable :: fun_cn_flex_a (:)   ! Parameter a of FUN-flexcn link code (def 5)
     real(r8), allocatable :: fun_cn_flex_b (:)   ! Parameter b of FUN-flexcn link code (def 200)
     real(r8), allocatable :: fun_cn_flex_c (:)   ! Parameter b of FUN-flexcn link code (def 80)         
     real(r8), allocatable :: FUN_fracfixers(:)   ! Fraction of C that can be used for fixation.    


     ! pft parameters for dynamic root code
     real(r8), allocatable :: root_dmx(:)     !maximum root depth

   contains

     procedure, public  :: Init
     procedure, public  :: InitForTesting ! version of Init meant for unit testing
     procedure, public  :: Clean
     procedure, private :: InitAllocate   
     procedure, private :: InitRead
     procedure, private :: set_is_pft_known_to_model   ! Set is_pft_known_to_model based on mergetoclmpft
     procedure, private :: set_num_cfts_known_to_model ! Set the module-level variable, num_cfts_known_to_model

  end type pftcon_type

  type(pftcon_type), public :: pftcon ! pft type constants structure

  integer, public, parameter :: pftname_len = 40         ! max length of pftname       
  character(len=pftname_len), public :: pftname(0:mxpft) ! PFT description

  real(r8), public, parameter :: reinickerp = 1.6_r8     ! parameter in allometric equation
  real(r8), public, parameter :: dwood  = 2.5e5_r8       ! cn wood density (gC/m3); lpj:2.0e5
  real(r8), public, parameter :: allom1 = 100.0_r8       ! parameters in
  real(r8), public, parameter :: allom2 =  40.0_r8       ! ...allometric
  real(r8), public, parameter :: allom3 =   0.5_r8       ! ...equations
  real(r8), public, parameter :: allom1s = 250.0_r8      ! modified for shrubs by
  real(r8), public, parameter :: allom2s =   8.0_r8      ! X.D.Z
! root radius, density from Bonan, GMD, 2014
  real(r8), public, parameter :: root_density = 0.31e06_r8 !(g biomass / m3 root)
  real(r8), public, parameter :: root_radius = 0.29e-03_r8 !(m)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this)

    class(pftcon_type) :: this

    call this%InitAllocate()
    call this%InitRead()

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitForTesting(this)
    ! Version of Init meant for unit testing
    !
    ! Allocate arrays, but don't try to read from file.
    !
    ! Values can then be set by tests as needed

    class(pftcon_type) :: this

    call this%InitAllocate()

  end subroutine InitForTesting

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this)
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) constants
    !
    ! !USES:
    use clm_varpar  ,  only: nvariants
    implicit none
    !
    ! !ARGUMENTS:
    class(pftcon_type) :: this
    !-----------------------------------------------------------------------

    allocate( this%noveg         (0:mxpft)); this%noveg (:)   =huge(1)
    allocate( this%tree          (0:mxpft)); this%tree  (:)   =huge(1)

    allocate( this%dleaf         (0:mxpft) )       
    allocate( this%c3psn         (0:mxpft) )       
    allocate( this%xl            (0:mxpft) )          
    allocate( this%rhol          (0:mxpft,numrad) ) 
    allocate( this%rhos          (0:mxpft,numrad) ) 
    allocate( this%taul          (0:mxpft,numrad) ) 
    allocate( this%taus          (0:mxpft,numrad) ) 
    allocate( this%z0mr          (0:mxpft) )        
    allocate( this%displar       (0:mxpft) )     
    allocate( this%roota_par     (0:mxpft) )   
    allocate( this%rootb_par     (0:mxpft) )   
    allocate( this%crop          (0:mxpft) )        
    allocate( this%mergetoclmpft (0:mxpft) )
    allocate( this%is_pft_known_to_model  (0:mxpft) )
    allocate( this%irrigated     (0:mxpft) )   
    allocate( this%smpso         (0:mxpft) )       
    allocate( this%smpsc         (0:mxpft) )       
    allocate( this%fnitr         (0:mxpft) )       
    allocate( this%slatop        (0:mxpft) )      
    allocate( this%dsladlai      (0:mxpft) )    
    allocate( this%leafcn        (0:mxpft) )      
    allocate( this%flnr          (0:mxpft) )        
    allocate( this%woody         (0:mxpft) )       
    allocate( this%lflitcn       (0:mxpft) )      
    allocate( this%frootcn       (0:mxpft) )      
    allocate( this%livewdcn      (0:mxpft) )     
    allocate( this%deadwdcn      (0:mxpft) )     
    allocate( this%grperc        (0:mxpft) )       
    allocate( this%grpnow        (0:mxpft) )       
    allocate( this%rootprof_beta (0:mxpft,nvariants) )
    allocate( this%graincn       (0:mxpft) )      
    allocate( this%mxtmp         (0:mxpft) )        
    allocate( this%baset         (0:mxpft) )        
    allocate( this%declfact      (0:mxpft) )     
    allocate( this%bfact         (0:mxpft) )        
    allocate( this%aleaff        (0:mxpft) )       
    allocate( this%arootf        (0:mxpft) )       
    allocate( this%astemf        (0:mxpft) )       
    allocate( this%arooti        (0:mxpft) )       
    allocate( this%fleafi        (0:mxpft) )       
    allocate( this%allconsl      (0:mxpft) )     
    allocate( this%allconss      (0:mxpft) )     
    allocate( this%ztopmx        (0:mxpft) )       
    allocate( this%laimx         (0:mxpft) )        
    allocate( this%gddmin        (0:mxpft) )       
    allocate( this%hybgdd        (0:mxpft) )       
    allocate( this%lfemerg       (0:mxpft) )      
    allocate( this%grnfill       (0:mxpft) )      
    allocate( this%mbbopt        (0:mxpft) )      
    allocate( this%medlynslope   (0:mxpft) )      
    allocate( this%medlynintercept(0:mxpft) )      
    allocate( this%mxmat         (0:mxpft) )        
    allocate( this%mnNHplantdate (0:mxpft) )
    allocate( this%mxNHplantdate (0:mxpft) )
    allocate( this%mnSHplantdate (0:mxpft) )
    allocate( this%mxSHplantdate (0:mxpft) )
    allocate( this%planttemp     (0:mxpft) )    
    allocate( this%minplanttemp  (0:mxpft) ) 
    allocate( this%froot_leaf    (0:mxpft) )   
    allocate( this%stem_leaf     (0:mxpft) )    
    allocate( this%croot_stem    (0:mxpft) )   
    allocate( this%flivewd       (0:mxpft) )      
    allocate( this%fcur          (0:mxpft) )         
    allocate( this%fcurdv        (0:mxpft) )       
    allocate( this%lf_flab       (0:mxpft) )      
    allocate( this%lf_fcel       (0:mxpft) )      
    allocate( this%lf_flig       (0:mxpft) )      
    allocate( this%fr_flab       (0:mxpft) )      
    allocate( this%fr_fcel       (0:mxpft) )      
    allocate( this%fr_flig       (0:mxpft) )      
    allocate( this%leaf_long     (0:mxpft) )
    allocate( this%evergreen     (0:mxpft) )    
    allocate( this%stress_decid  (0:mxpft) ) 
    allocate( this%season_decid  (0:mxpft) ) 
    allocate( this%dwood         (0:mxpft) )
    allocate( this%root_density  (0:mxpft) )
    allocate( this%root_radius   (0:mxpft) )
    allocate( this%pconv         (0:mxpft) )        
    allocate( this%pprod10       (0:mxpft) )      
    allocate( this%pprod100      (0:mxpft) )     
    allocate( this%pprodharv10   (0:mxpft) )  
    allocate( this%cc_leaf       (0:mxpft) )
    allocate( this%cc_lstem      (0:mxpft) )
    allocate( this%cc_dstem      (0:mxpft) )
    allocate( this%cc_other      (0:mxpft) )
    allocate( this%fm_leaf       (0:mxpft) )
    allocate( this%fm_lstem      (0:mxpft) )
    allocate( this%fm_dstem      (0:mxpft) )
    allocate( this%fm_other      (0:mxpft) )
    allocate( this%fm_root       (0:mxpft) )
    allocate( this%fm_lroot      (0:mxpft) )
    allocate( this%fm_droot      (0:mxpft) )
    allocate( this%fsr_pft       (0:mxpft) )
    allocate( this%fd_pft        (0:mxpft) )
    allocate( this%manunitro     (0:mxpft) )
    allocate( this%fleafcn       (0:mxpft) )  
    allocate( this%ffrootcn      (0:mxpft) ) 
    allocate( this%fstemcn       (0:mxpft) )  
    allocate( this%i_vcad        (0:mxpft) )
    allocate( this%s_vcad        (0:mxpft) )
    allocate( this%i_flnr        (0:mxpft) )
    allocate( this%s_flnr        (0:mxpft) )
    allocate( this%pftpar20      (0:mxpft) )   
    allocate( this%pftpar28      (0:mxpft) )   
    allocate( this%pftpar29      (0:mxpft) )   
    allocate( this%pftpar30      (0:mxpft) )   
    allocate( this%pftpar31      (0:mxpft) )   
    allocate( this%a_fix         (0:mxpft) )
    allocate( this%b_fix         (0:mxpft) )
    allocate( this%c_fix         (0:mxpft) )
    allocate( this%s_fix         (0:mxpft) )
    allocate( this%akc_active    (0:mxpft) )
    allocate( this%akn_active    (0:mxpft) )
    allocate( this%ekc_active    (0:mxpft) )
    allocate( this%ekn_active    (0:mxpft) )  
    allocate( this%kc_nonmyc     (0:mxpft) )
    allocate( this%kn_nonmyc     (0:mxpft) )
    allocate( this%kr_resorb     (0:mxpft) )
    allocate( this%perecm        (0:mxpft) )
    allocate( this%root_dmx      (0:mxpft) )
    allocate( this%fun_cn_flex_a (0:mxpft) )
    allocate( this%fun_cn_flex_b (0:mxpft) )
    allocate( this%fun_cn_flex_c (0:mxpft) )
    allocate( this%FUN_fracfixers(0:mxpft) )
    
 
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitRead(this)
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) constants
    !
    ! !USES:
    use shr_log_mod , only : errMsg => shr_log_errMsg
    use fileutils   , only : getfil
    use ncdio_pio   , only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t
    use ncdio_pio   , only : ncd_inqdid, ncd_inqdlen
    use clm_varctl  , only : paramfile, use_fates, use_flexibleCN, use_dynroot
    use spmdMod     , only : masterproc
    use CLMFatesParamInterfaceMod, only : FatesReadPFTs
    !
    ! !ARGUMENTS:
    class(pftcon_type) :: this
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn                ! local file name
    integer            :: i,n,m                ! loop indices
    integer            :: ier                  ! error code
    type(file_desc_t)  :: ncid                 ! pio netCDF file id
    integer            :: dimid                ! netCDF dimension id
    integer            :: npft                 ! number of pfts on pft-physiology file
    logical            :: readv                ! read variable in or not
    character(len=32)  :: subname = 'InitRead' ! subroutine name
    character(len=pftname_len) :: expected_pftnames(0:mxpft) 
    character(len=512) :: msg
    !-----------------------------------------------------------------------
    !
    ! Expected PFT names: The names expected on the paramfile file and the order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending with broadleaf_deciduous_boreal_tree
    !       then shrubs, ending with broadleaf_deciduous_boreal_shrub, then grasses starting with c3_arctic_grass
    !       and finally crops, ending with irrigated_tropical_soybean
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF THE CODE WHERE THE ORDER MATTERS!

    expected_pftnames( 0) = 'not_vegetated                      '
    expected_pftnames( 1) = 'needleleaf_evergreen_temperate_tree'
    expected_pftnames( 2) = 'needleleaf_evergreen_boreal_tree   '
    expected_pftnames( 3) = 'needleleaf_deciduous_boreal_tree   '
    expected_pftnames( 4) = 'broadleaf_evergreen_tropical_tree  '
    expected_pftnames( 5) = 'broadleaf_evergreen_temperate_tree '
    expected_pftnames( 6) = 'broadleaf_deciduous_tropical_tree  '
    expected_pftnames( 7) = 'broadleaf_deciduous_temperate_tree '
    expected_pftnames( 8) = 'broadleaf_deciduous_boreal_tree    '
    expected_pftnames( 9) = 'broadleaf_evergreen_shrub          '
    expected_pftnames(10) = 'broadleaf_deciduous_temperate_shrub'
    expected_pftnames(11) = 'broadleaf_deciduous_boreal_shrub   '
    expected_pftnames(12) = 'c3_arctic_grass                    '
    expected_pftnames(13) = 'c3_non-arctic_grass                '
    expected_pftnames(14) = 'c4_grass                           '
    expected_pftnames(15) = 'c3_crop                            '
    expected_pftnames(16) = 'c3_irrigated                       '
    expected_pftnames(17) = 'temperate_corn                     '
    expected_pftnames(18) = 'irrigated_temperate_corn           '
    expected_pftnames(19) = 'spring_wheat                       '
    expected_pftnames(20) = 'irrigated_spring_wheat             '
    expected_pftnames(21) = 'winter_wheat                       '
    expected_pftnames(22) = 'irrigated_winter_wheat             '
    expected_pftnames(23) = 'temperate_soybean                  '
    expected_pftnames(24) = 'irrigated_temperate_soybean        '
    expected_pftnames(25) = 'barley                             '
    expected_pftnames(26) = 'irrigated_barley                   '
    expected_pftnames(27) = 'winter_barley                      '
    expected_pftnames(28) = 'irrigated_winter_barley            '
    expected_pftnames(29) = 'rye                                '
    expected_pftnames(30) = 'irrigated_rye                      '
    expected_pftnames(31) = 'winter_rye                         '
    expected_pftnames(32) = 'irrigated_winter_rye               '
    expected_pftnames(33) = 'cassava                            '
    expected_pftnames(34) = 'irrigated_cassava                  '
    expected_pftnames(35) = 'citrus                             '
    expected_pftnames(36) = 'irrigated_citrus                   '
    expected_pftnames(37) = 'cocoa                              '
    expected_pftnames(38) = 'irrigated_cocoa                    '
    expected_pftnames(39) = 'coffee                             '
    expected_pftnames(40) = 'irrigated_coffee                   '
    expected_pftnames(41) = 'cotton                             '
    expected_pftnames(42) = 'irrigated_cotton                   '
    expected_pftnames(43) = 'datepalm                           '
    expected_pftnames(44) = 'irrigated_datepalm                 '
    expected_pftnames(45) = 'foddergrass                        '
    expected_pftnames(46) = 'irrigated_foddergrass              '
    expected_pftnames(47) = 'grapes                             '
    expected_pftnames(48) = 'irrigated_grapes                   '
    expected_pftnames(49) = 'groundnuts                         '
    expected_pftnames(50) = 'irrigated_groundnuts               '
    expected_pftnames(51) = 'millet                             '
    expected_pftnames(52) = 'irrigated_millet                   '
    expected_pftnames(53) = 'oilpalm                            '
    expected_pftnames(54) = 'irrigated_oilpalm                  '
    expected_pftnames(55) = 'potatoes                           '
    expected_pftnames(56) = 'irrigated_potatoes                 '
    expected_pftnames(57) = 'pulses                             '
    expected_pftnames(58) = 'irrigated_pulses                   '
    expected_pftnames(59) = 'rapeseed                           '
    expected_pftnames(60) = 'irrigated_rapeseed                 '
    expected_pftnames(61) = 'rice                               '
    expected_pftnames(62) = 'irrigated_rice                     '
    expected_pftnames(63) = 'sorghum                            '
    expected_pftnames(64) = 'irrigated_sorghum                  '
    expected_pftnames(65) = 'sugarbeet                          '
    expected_pftnames(66) = 'irrigated_sugarbeet                '
    expected_pftnames(67) = 'sugarcane                          '
    expected_pftnames(68) = 'irrigated_sugarcane                '
    expected_pftnames(69) = 'sunflower                          '
    expected_pftnames(70) = 'irrigated_sunflower                '
    expected_pftnames(71) = 'miscanthus                         '
    expected_pftnames(72) = 'irrigated_miscanthus               '
    expected_pftnames(73) = 'switchgrass                        '
    expected_pftnames(74) = 'irrigated_switchgrass              '
    expected_pftnames(75) = 'tropical_corn                      '
    expected_pftnames(76) = 'irrigated_tropical_corn            '
    expected_pftnames(77) = 'tropical_soybean                   '
    expected_pftnames(78) = 'irrigated_tropical_soybean         '

    ! Set specific vegetation type values

    if (masterproc) then
       write(iulog,*) 'Attempting to read PFT physiological data .....'
    end if
    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid, 'pft', dimid)
    call ncd_inqdlen(ncid, dimid, npft)

    if (npft - 1 /= mxpft) then
       ! NOTE(bja, 201503) need to subtract 1 because of indexing.
       ! NOTE(bja, 201503) fail early because one of the io libs
       ! throws a useless abort error message deep inside the stack
       ! instead of returning readv so we can get a useful line
       ! number.
       write(msg, '(a, i4, a, i4, a)') "ERROR: The number of pfts in the input netcdf file (", &
            npft, ") does not equal the expected number of pfts (", mxpft, "). "
       call endrun(msg=trim(msg)//errMsg(sourcefile, __LINE__))
    end if

    call ncd_io('pftname',pftname, 'read', ncid, readvar=readv, posNOTonfile=.true.) 
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('z0mr', this%z0mr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('displar', this%displar, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('dleaf', this%dleaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('c3psn', this%c3psn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rholvis', this%rhol(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rholnir', this%rhol(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rhosvis', this%rhos(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rhosnir', this% rhos(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('taulvis', this%taul(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('taulnir', this%taul(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('tausvis', this%taus(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('tausnir', this%taus(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('xl', this%xl, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('roota_par', this%roota_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rootb_par', this%rootb_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('slatop', this%slatop, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('dsladlai', this%dsladlai, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('leafcn', this%leafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('flnr', this%flnr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('smpso', this%smpso, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('smpsc', this%smpsc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fnitr', this%fnitr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('woody', this%woody, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lflitcn', this%lflitcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('frootcn', this%frootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('livewdcn', this%livewdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('deadwdcn', this%deadwdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('grperc', this%grperc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('grpnow', this%grpnow, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('froot_leaf', this%froot_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('stem_leaf', this%stem_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('croot_stem', this%croot_stem, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('flivewd', this%flivewd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fcur', this%fcur, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fcurdv', this%fcurdv, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lf_flab', this%lf_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lf_fcel', this%lf_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lf_flig', this%lf_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fr_flab', this%fr_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fr_fcel', this%fr_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fr_flig', this%fr_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('leaf_long', this%leaf_long, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('evergreen', this%evergreen, 'read', ncid, readvar=readv, posNOTonfile=.true.)    
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('stress_decid', this%stress_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('season_decid', this%season_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar20', this%pftpar20, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar28', this%pftpar28, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar29', this%pftpar29, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar30', this%pftpar30, 'read', ncid, readvar=readv, posNOTonfile=.true.)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar31', this%pftpar31, 'read', ncid, readvar=readv, posNOTonfile=.true.)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('a_fix', this%a_fix, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
   
    call ncd_io('b_fix', this%b_fix, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    
    call ncd_io('c_fix', this%c_fix, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    
    call ncd_io('s_fix', this%s_fix, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('akc_active', this%akc_active, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('akn_active', this%akn_active, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('ekc_active', this%ekc_active, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('ekn_active', this%ekn_active, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('kc_nonmyc', this%kc_nonmyc, 'read', ncid, readvar=readv,   posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
   
    call ncd_io('kn_nonmyc', this%kn_nonmyc, 'read', ncid, readvar=readv,   posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('kr_resorb', this%kr_resorb, 'read', ncid, readvar=readv,   posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('perecm', this%perecm, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fun_cn_flex_a', this%fun_cn_flex_a, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fun_cn_flex_b', this%fun_cn_flex_b, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fun_cn_flex_c', this%fun_cn_flex_c, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
  
    call ncd_io('FUN_fracfixers', this%FUN_fracfixers, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('manunitro', this%manunitro, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fleafcn', this%fleafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('ffrootcn', this%ffrootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fstemcn', this%fstemcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rootprof_beta', this%rootprof_beta, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pconv', this%pconv, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pprod10', this%pprod10, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pprodharv10', this%pprodharv10, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pprod100', this%pprod100, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('graincn', this%graincn, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('mxtmp', this%mxtmp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('baset', this%baset, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('declfact', this%declfact, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('bfact', this%bfact, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('aleaff', this%aleaff, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('arootf', this%arootf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('astemf', this%astemf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('arooti', this%arooti, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fleafi', this%fleafi, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('allconsl', this%allconsl, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('allconss', this%allconss, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('crop', this%crop, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('mergetoclmpft', this%mergetoclmpft, 'read', ncid, readvar=readv)  
    if ( .not. readv ) then
       call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    end if

    call ncd_io('irrigated', this%irrigated, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('ztopmx', this%ztopmx, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('laimx', this%laimx, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('gddmin', this%gddmin, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('hybgdd', this%hybgdd, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lfemerg', this%lfemerg, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('grnfill', this%grnfill, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('mbbopt', this%mbbopt, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('medlynslope', this%medlynslope, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('medlynintercept', this%medlynintercept, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('mxmat', this%mxmat, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('cc_leaf', this% cc_leaf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('cc_lstem', this%cc_lstem, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('cc_dstem', this%cc_dstem, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('cc_other', this%cc_other, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_leaf', this% fm_leaf, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_lstem', this%fm_lstem, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_dstem', this%fm_dstem, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_other', this%fm_other, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_root', this% fm_root, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_lroot', this%fm_lroot, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_droot', this%fm_droot, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fsr_pft', this% fsr_pft, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fd_pft', this%  fd_pft, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('planting_temp', this%planttemp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('min_planting_temp', this%minplanttemp, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('min_NH_planting_date', this%mnNHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('min_SH_planting_date', this%mnSHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('max_NH_planting_date', this%mxNHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('max_SH_planting_date', this%mxSHplantdate, 'read', ncid, readvar=readv)  
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    !
    ! Constants
    !
    !MV (10-08-14) TODO is this right - used to be maxveg - is it okay to set it to mxpft?
    do m = 0,mxpft 
       this%dwood(m) = dwood
       this%root_radius(m)  = root_radius
       this%root_density(m) = root_density

       if (m <= ntree) then
          this%tree(m) = 1
       else
          this%tree(m) = 0
       end if
    end do
    !
    ! clm 5 nitrogen variables
    !
    if (use_flexibleCN) then
       call ncd_io('i_vcad', this%i_vcad, 'read', ncid, readvar=readv) 
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__)) 
       
       call ncd_io('s_vcad', this%s_vcad, 'read', ncid, readvar=readv) 
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__)) 
       
       call ncd_io('i_flnr', this%i_flnr, 'read', ncid, readvar=readv) 
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__)) 
       
       call ncd_io('s_flnr', this%s_flnr, 'read', ncid, readvar=readv) 
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__)) 
    end if

    !
    ! Dynamic Root variables for crops
    !
    if ( use_crop .and. use_dynroot )then
       call ncd_io('root_dmx', this%root_dmx, 'read', ncid, readvar=readv)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    end if
   
    call ncd_pio_closefile(ncid)

    call FatesReadPFTs()

    do i = 0, mxpft
       if (.not. use_fates)then
          if ( trim(adjustl(pftname(i))) /= trim(expected_pftnames(i)) )then
             write(iulog,*)'pftconrd: pftname is NOT what is expected, name = ', &
                  trim(pftname(i)), ', expected name = ', trim(expected_pftnames(i))
             call endrun(msg='pftconrd: bad name for pft on paramfile dataset'//errMsg(sourcefile, __LINE__))
          end if
       end if

       if ( trim(pftname(i)) == 'not_vegetated'                       ) noveg                = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_temperate_tree' ) ndllf_evr_tmp_tree   = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_boreal_tree'    ) ndllf_evr_brl_tree   = i
       if ( trim(pftname(i)) == 'needleleaf_deciduous_boreal_tree'    ) ndllf_dcd_brl_tree   = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_tropical_tree'   ) nbrdlf_evr_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_temperate_tree'  ) nbrdlf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_tropical_tree'   ) nbrdlf_dcd_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_tree'  ) nbrdlf_dcd_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_tree'     ) nbrdlf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_shrub'           ) nbrdlf_evr_shrub     = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_shrub' ) nbrdlf_dcd_tmp_shrub = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_shrub'    ) nbrdlf_dcd_brl_shrub = i
       if ( trim(pftname(i)) == 'c3_arctic_grass'                     ) nc3_arctic_grass     = i
       if ( trim(pftname(i)) == 'c3_non-arctic_grass'                 ) nc3_nonarctic_grass  = i
       if ( trim(pftname(i)) == 'c4_grass'                            ) nc4_grass            = i
       if ( trim(pftname(i)) == 'c3_crop'                             ) nc3crop              = i
       if ( trim(pftname(i)) == 'c3_irrigated'                        ) nc3irrig             = i
       if ( trim(pftname(i)) == 'temperate_corn'                      ) ntmp_corn            = i
       if ( trim(pftname(i)) == 'irrigated_temperate_corn'            ) nirrig_tmp_corn      = i
       if ( trim(pftname(i)) == 'spring_wheat'                        ) nswheat              = i
       if ( trim(pftname(i)) == 'irrigated_spring_wheat'              ) nirrig_swheat        = i
       if ( trim(pftname(i)) == 'winter_wheat'                        ) nwwheat              = i
       if ( trim(pftname(i)) == 'irrigated_winter_wheat'              ) nirrig_wwheat        = i
       if ( trim(pftname(i)) == 'temperate_soybean'                   ) ntmp_soybean         = i
       if ( trim(pftname(i)) == 'irrigated_temperate_soybean'         ) nirrig_tmp_soybean   = i
       if ( trim(pftname(i)) == 'barley'                              ) nbarley              = i
       if ( trim(pftname(i)) == 'irrigated_barley'                    ) nirrig_barley        = i
       if ( trim(pftname(i)) == 'winter_barley'                       ) nwbarley             = i
       if ( trim(pftname(i)) == 'irrigated_winter_barley'             ) nirrig_wbarley       = i
       if ( trim(pftname(i)) == 'rye'                                 ) nrye                 = i
       if ( trim(pftname(i)) == 'irrigated_rye'                       ) nirrig_rye           = i
       if ( trim(pftname(i)) == 'winter_rye'                          ) nwrye                = i
       if ( trim(pftname(i)) == 'irrigated_winter_rye'                ) nirrig_wrye          = i
       if ( trim(pftname(i)) == 'cassava'                             ) ncassava             = i
       if ( trim(pftname(i)) == 'irrigated_cassava'                   ) nirrig_cassava       = i
       if ( trim(pftname(i)) == 'citrus'                              ) ncitrus              = i
       if ( trim(pftname(i)) == 'irrigated_citrus'                    ) nirrig_citrus        = i
       if ( trim(pftname(i)) == 'cocoa'                               ) ncocoa               = i
       if ( trim(pftname(i)) == 'irrigated_cocoa'                     ) nirrig_cocoa         = i
       if ( trim(pftname(i)) == 'coffee'                              ) ncoffee              = i
       if ( trim(pftname(i)) == 'irrigated_coffee'                    ) nirrig_coffee        = i
       if ( trim(pftname(i)) == 'cotton'                              ) ncotton              = i
       if ( trim(pftname(i)) == 'irrigated_cotton'                    ) nirrig_cotton        = i
       if ( trim(pftname(i)) == 'datepalm'                            ) ndatepalm            = i
       if ( trim(pftname(i)) == 'irrigated_datepalm'                  ) nirrig_datepalm      = i
       if ( trim(pftname(i)) == 'foddergrass'                         ) nfoddergrass         = i
       if ( trim(pftname(i)) == 'irrigated_foddergrass'               ) nirrig_foddergrass   = i
       if ( trim(pftname(i)) == 'grapes'                              ) ngrapes              = i
       if ( trim(pftname(i)) == 'irrigated_grapes'                    ) nirrig_grapes        = i
       if ( trim(pftname(i)) == 'groundnuts'                          ) ngroundnuts          = i
       if ( trim(pftname(i)) == 'irrigated_groundnuts'                ) nirrig_groundnuts    = i
       if ( trim(pftname(i)) == 'millet'                              ) nmillet              = i
       if ( trim(pftname(i)) == 'irrigated_millet'                    ) nirrig_millet        = i
       if ( trim(pftname(i)) == 'oilpalm'                             ) noilpalm             = i
       if ( trim(pftname(i)) == 'irrigated_oilpalm'                   ) nirrig_oilpalm       = i
       if ( trim(pftname(i)) == 'potatoes'                            ) npotatoes            = i
       if ( trim(pftname(i)) == 'irrigated_potatoes'                  ) nirrig_potatoes      = i
       if ( trim(pftname(i)) == 'pulses'                              ) npulses              = i
       if ( trim(pftname(i)) == 'irrigated_pulses'                    ) nirrig_pulses        = i
       if ( trim(pftname(i)) == 'rapeseed'                            ) nrapeseed            = i
       if ( trim(pftname(i)) == 'irrigated_rapeseed'                  ) nirrig_rapeseed      = i
       if ( trim(pftname(i)) == 'rice'                                ) nrice                = i
       if ( trim(pftname(i)) == 'irrigated_rice'                      ) nirrig_rice          = i
       if ( trim(pftname(i)) == 'sorghum'                             ) nsorghum             = i
       if ( trim(pftname(i)) == 'irrigated_sorghum'                   ) nirrig_sorghum       = i
       if ( trim(pftname(i)) == 'sugarbeet'                           ) nsugarbeet           = i
       if ( trim(pftname(i)) == 'irrigated_sugarbeet'                 ) nirrig_sugarbeet     = i
       if ( trim(pftname(i)) == 'sugarcane'                           ) nsugarcane           = i
       if ( trim(pftname(i)) == 'irrigated_sugarcane'                 ) nirrig_sugarcane     = i
       if ( trim(pftname(i)) == 'sunflower'                           ) nsunflower           = i
       if ( trim(pftname(i)) == 'irrigated_sunflower'                 ) nirrig_sunflower     = i
       if ( trim(pftname(i)) == 'miscanthus'                          ) nmiscanthus          = i
       if ( trim(pftname(i)) == 'irrigated_miscanthus'                ) nirrig_miscanthus    = i
       if ( trim(pftname(i)) == 'switchgrass'                         ) nswitchgrass         = i
       if ( trim(pftname(i)) == 'irrigated_switchgrass'               ) nirrig_switchgrass   = i
       if ( trim(pftname(i)) == 'tropical_corn'                       ) ntrp_corn            = i
       if ( trim(pftname(i)) == 'irrigated_tropical_corn'             ) nirrig_trp_corn      = i
       if ( trim(pftname(i)) == 'tropical_soybean'                    ) ntrp_soybean         = i
       if ( trim(pftname(i)) == 'irrigated_tropical_soybean'          ) nirrig_trp_soybean   = i
    end do

    ntree                = nbrdlf_dcd_brl_tree  ! value for last type of tree
    npcropmin            = ntmp_corn            ! first prognostic crop
    npcropmax            = mxpft                ! last prognostic crop in list

    call this%set_is_pft_known_to_model()
    call this%set_num_cfts_known_to_model()

    if (use_cndv) then
       this%fcur(:) = this%fcurdv(:)
    end if
    ! When crop is not on, merge prognostic crop types into either the rainfed
    ! or irrigated C3 generic crop types
    if ( .not. use_crop )then
       do i = npcropmin, ntrp_soybean, 2
         this%mergetoclmpft(i) = nc3crop
       end do
       do i = nirrig_tmp_corn, npcropmax, 2
         this%mergetoclmpft(i) = nc3irrig
       end do
    end if
    !
    ! Do some error checking, but not if fates is on.
    !
    ! FIX(SPM,032414) double check if some of these should be on...

    if( .not. use_fates ) then
       if ( npcropmax /= mxpft )then
          call endrun(msg=' ERROR: npcropmax is NOT the last value'//errMsg(sourcefile, __LINE__))
       end if
       do i = 0, mxpft
          if ( this%irrigated(i) == 1.0_r8 .and.                              &
               (i == nc3irrig               .or.                              &
                i == nirrig_tmp_corn        .or.                              &
                i == nirrig_swheat          .or. i == nirrig_wwheat      .or. &
                i == nirrig_tmp_soybean     .or.                              &
                i == nirrig_barley          .or. i == nirrig_wbarley     .or. &
                i == nirrig_rye             .or. i == nirrig_wrye        .or. &
                i == nirrig_cassava         .or.                              &
                i == nirrig_citrus          .or.                              &
                i == nirrig_cocoa           .or. i == nirrig_coffee      .or. &
                i == nirrig_cotton          .or.                              &
                i == nirrig_datepalm        .or.                              &
                i == nirrig_foddergrass     .or.                              &
                i == nirrig_grapes          .or. i == nirrig_groundnuts  .or. &
                i == nirrig_millet          .or.                              &
                i == nirrig_oilpalm         .or.                              &
                i == nirrig_potatoes        .or. i == nirrig_pulses      .or. &
                i == nirrig_rapeseed        .or. i == nirrig_rice        .or. &
                i == nirrig_sorghum         .or.                              &
                i == nirrig_sugarbeet       .or. i == nirrig_sugarcane   .or. &
                i == nirrig_sunflower       .or.                              &
                i == nirrig_miscanthus      .or. i == nirrig_switchgrass .or. &
                i == nirrig_trp_corn        .or.                              &
                i == nirrig_trp_soybean) )then
             ! correct
          else if ( this%irrigated(i) == 0.0_r8 )then
             ! correct
          else
             call endrun(msg=' ERROR: irrigated has wrong values'//errMsg(sourcefile, __LINE__))
          end if
          if (      this%crop(i) == 1.0_r8 .and. (i >= nc3crop .and. i <= npcropmax) )then
             ! correct
          else if ( this%crop(i) == 0.0_r8 )then
             ! correct
          else
             call endrun(msg=' ERROR: crop has wrong values'//errMsg(sourcefile, __LINE__))
          end if
          if ( (i /= noveg) .and. (i < npcropmin) .and. &
               abs(this%pconv(i) + this%pprod10(i) + this%pprod100(i) - 1.0_r8) > 1.e-7_r8 )then
             call endrun(msg=' ERROR: pconv+pprod10+pprod100 do NOT sum to one.'//errMsg(sourcefile, __LINE__))
          end if
          if ( this%pprodharv10(i) > 1.0_r8 .or. this%pprodharv10(i) < 0.0_r8 )then
             call endrun(msg=' ERROR: pprodharv10 outside of range.'//errMsg(sourcefile, __LINE__))
          end if
       end do
    end if

    if (masterproc) then
       write(iulog,*) 'Successfully read PFT physiological data'
       write(iulog,*)
    end if

  end subroutine InitRead

  !-----------------------------------------------------------------------
  subroutine set_is_pft_known_to_model(this)
    !
    ! !DESCRIPTION:
    ! Set is_pft_known_to_model based on mergetoclmpft
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(pftcon_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: m, merge_type

    character(len=*), parameter :: subname = 'set_is_pft_known_to_model'
    !-----------------------------------------------------------------------

    this%is_pft_known_to_model(:) = .false.

    ! NOTE(wjs, 2015-10-04) Currently, type 0 has mergetoclmpft = _FillValue in the file,
    ! so we can't handle it in the general loop below. But CLM always uses type 0, so
    ! handle it specially here.
    this%is_pft_known_to_model(0) = .true.

    ! NOTE(wjs, 2015-10-04) Currently, mergetoclmpft is only used for crop types.
    ! However, we handle it more generally here (treating ALL pft types), in case its use
    ! is ever extended to work with non-crop types as well.
    do m = 1, mxpft
       merge_type = this%mergetoclmpft(m)
       this%is_pft_known_to_model(merge_type) = .true.
    end do

  end subroutine set_is_pft_known_to_model

  !-----------------------------------------------------------------------
  subroutine set_num_cfts_known_to_model(this)
    !
    ! !DESCRIPTION:
    ! Set the module-level variable, num_cfts_known_to_model
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(pftcon_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: m

    character(len=*), parameter :: subname = 'set_num_cfts_known_to_model'
    !-----------------------------------------------------------------------

    num_cfts_known_to_model = 0
    do m = cft_lb, cft_ub
       if (this%is_pft_known_to_model(m)) then
          num_cfts_known_to_model = num_cfts_known_to_model + 1
       end if
    end do

  end subroutine set_num_cfts_known_to_model

  !-----------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Deallocate memory
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(pftcon_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    deallocate( this%noveg)
    deallocate( this%tree)

    deallocate( this%dleaf)
    deallocate( this%c3psn)
    deallocate( this%xl)
    deallocate( this%rhol)
    deallocate( this%rhos)
    deallocate( this%taul)
    deallocate( this%taus)
    deallocate( this%z0mr)
    deallocate( this%displar)
    deallocate( this%roota_par)
    deallocate( this%rootb_par)
    deallocate( this%crop)
    deallocate( this%mergetoclmpft)
    deallocate( this%is_pft_known_to_model)
    deallocate( this%irrigated)
    deallocate( this%smpso)
    deallocate( this%smpsc)
    deallocate( this%fnitr)
    deallocate( this%slatop)
    deallocate( this%dsladlai)
    deallocate( this%leafcn)
    deallocate( this%flnr)
    deallocate( this%woody)
    deallocate( this%lflitcn)
    deallocate( this%frootcn)
    deallocate( this%livewdcn)
    deallocate( this%deadwdcn)
    deallocate( this%grperc)
    deallocate( this%grpnow)
    deallocate( this%rootprof_beta)
    deallocate( this%graincn)
    deallocate( this%mxtmp)
    deallocate( this%baset)
    deallocate( this%declfact)
    deallocate( this%bfact)
    deallocate( this%aleaff)
    deallocate( this%arootf)
    deallocate( this%astemf)
    deallocate( this%arooti)
    deallocate( this%fleafi)
    deallocate( this%allconsl)
    deallocate( this%allconss)
    deallocate( this%ztopmx)
    deallocate( this%laimx)
    deallocate( this%gddmin)
    deallocate( this%hybgdd)
    deallocate( this%lfemerg)
    deallocate( this%grnfill)
    deallocate( this%mbbopt)
    deallocate( this%medlynslope)
    deallocate( this%medlynintercept)
    deallocate( this%mxmat)
    deallocate( this%mnNHplantdate)
    deallocate( this%mxNHplantdate)
    deallocate( this%mnSHplantdate)
    deallocate( this%mxSHplantdate)
    deallocate( this%planttemp)
    deallocate( this%minplanttemp)
    deallocate( this%froot_leaf)
    deallocate( this%stem_leaf)
    deallocate( this%croot_stem)
    deallocate( this%flivewd)
    deallocate( this%fcur)
    deallocate( this%fcurdv)
    deallocate( this%lf_flab)
    deallocate( this%lf_fcel)
    deallocate( this%lf_flig)
    deallocate( this%fr_flab)
    deallocate( this%fr_fcel)
    deallocate( this%fr_flig)
    deallocate( this%leaf_long)
    deallocate( this%evergreen)
    deallocate( this%stress_decid)
    deallocate( this%season_decid)
    deallocate( this%dwood)
    deallocate( this%root_density)
    deallocate( this%root_radius)
    deallocate( this%pconv)
    deallocate( this%pprod10)
    deallocate( this%pprod100)
    deallocate( this%pprodharv10)
    deallocate( this%cc_leaf)
    deallocate( this%cc_lstem)
    deallocate( this%cc_dstem)
    deallocate( this%cc_other)
    deallocate( this%fm_leaf)
    deallocate( this%fm_lstem)
    deallocate( this%fm_dstem)
    deallocate( this%fm_other)
    deallocate( this%fm_root)
    deallocate( this%fm_lroot)
    deallocate( this%fm_droot)
    deallocate( this%fsr_pft)
    deallocate( this%fd_pft)
    deallocate( this%manunitro)
    deallocate( this%fleafcn)
    deallocate( this%ffrootcn)
    deallocate( this%fstemcn)
    deallocate( this%i_vcad)
    deallocate( this%s_vcad)
    deallocate( this%i_flnr)
    deallocate( this%s_flnr)
    deallocate( this%pftpar20)
    deallocate( this%pftpar28)
    deallocate( this%pftpar29)
    deallocate( this%pftpar30)
    deallocate( this%pftpar31)
    deallocate( this%a_fix)
    deallocate( this%b_fix)
    deallocate( this%c_fix)
    deallocate( this%s_fix)
    deallocate( this%akc_active)
    deallocate( this%akn_active)
    deallocate( this%ekc_active)
    deallocate( this%ekn_active)
    deallocate( this%kc_nonmyc)
    deallocate( this%kn_nonmyc)
    deallocate( this%kr_resorb)
    deallocate( this%perecm)
    deallocate( this%root_dmx)
    deallocate( this%fun_cn_flex_a)
    deallocate( this%fun_cn_flex_b)
    deallocate( this%fun_cn_flex_c)
    deallocate( this%FUN_fracfixers)
    
  end subroutine Clean

end module pftconMod

