module  PhotosynthesisMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf photosynthesis and stomatal conductance calculation as described by
  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
  ! a multi-layer canopy
  !
  ! !USES:
  use shr_sys_mod         , only : shr_sys_flush
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
  use abortutils          , only : endrun
  use clm_varctl          , only : use_c13, use_c14, use_cn, use_cndv, use_fates, use_luna, use_hydrstress
  use clm_varctl          , only : iulog
  use clm_varpar          , only : nlevcan, nvegwcs, mxpft
  use clm_varcon          , only : namep, c14ratio, spval
  use decompMod           , only : bounds_type
  use QuadraticMod        , only : quadratic
  use pftconMod           , only : pftcon
  use CIsoAtmTimeseriesMod, only : C14BombSpike, use_c14_bombspike, C13TimeSeries, use_c13_timeseries, nsectors_c14
  use atm2lndType         , only : atm2lnd_type
  use CanopyStateType     , only : canopystate_type
  use WaterDiagnosticBulkType      , only : waterdiagnosticbulk_type
  use WaterFluxBulkType       , only : waterfluxbulk_type
  use SoilStateType       , only : soilstate_type
  use TemperatureType     , only : temperature_type
  use SolarAbsorbedType   , only : solarabs_type
  use SurfaceAlbedoType   , only : surfalb_type
  use OzoneBaseMod        , only : ozone_base_type
  use LandunitType        , only : lun
  use PatchType           , only : patch
  use GridcellType        , only : grc
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Photosynthesis        ! Leaf stomatal resistance and leaf photosynthesis
  public :: PhotosynthesisTotal   ! Determine of total photosynthesis
  public :: Fractionation         ! C13 fractionation during photosynthesis
  ! For plant hydraulics approach
  public :: PhotosynthesisHydraulicStress ! Leaf stomatal resistance and leaf photosynthesis
                                          ! Simultaneous solution of sunlit/shaded per Pierre
                                          ! Gentine/Daniel Kennedy plant hydraulic stress method
  public :: plc                           ! Return value of vulnerability curve at x

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: hybrid         ! hybrid solver for ci
  private :: ci_func        ! ci function
  private :: brent          ! brent solver for root of a single variable function
  private :: ft             ! photosynthesis temperature response
  private :: fth            ! photosynthesis temperature inhibition
  private :: fth25          ! scaling factor for photosynthesis temperature inhibition
  ! For plant hydraulics approach
  private :: hybrid_PHS     ! hybrid solver for ci
  private :: ci_func_PHS    ! ci function
  private :: brent_PHS      ! brent solver for root of a single variable function
  private :: calcstress     ! compute the root water stress
  private :: getvegwp       ! calculate vegetation water potential (sun, sha, xylem, root)
  private :: getqflx        ! calculate sunlit and shaded transpiration
  private :: spacF          ! flux divergence across each vegetation segment
  private :: spacA          ! the inverse Jacobian matrix relating delta(vegwp) to f, d(vegwp)=A*f
  private :: d1plc          ! compute 1st deriv of conductance attenuation for each segment

  ! !PRIVATE DATA:
  integer, parameter, private :: leafresp_mtd_ryan1991  = 1  ! Ryan 1991 method for lmr25top
  integer, parameter, private :: leafresp_mtd_atkin2015 = 2  ! Atkin 2015 method for lmr25top
  integer, parameter, private :: sun=1     ! index for sunlit
  integer, parameter, private :: sha=2     ! index for shaded
  integer, parameter, private :: xyl=3     ! index for xylem
  integer, parameter, private :: root=4    ! index for root
  integer, parameter, private :: veg=0     ! index for vegetation
  integer, parameter, private :: soil=1    ! index for soil
  integer, parameter, private :: stomatalcond_mtd_bb1987     = 1   ! Ball-Berry 1987 method for photosynthesis
  integer, parameter, private :: stomatalcond_mtd_medlyn2011 = 2   ! Medlyn 2011 method for photosynthesis
  ! !PUBLIC VARIABLES:

  type :: photo_params_type
     real(r8), allocatable, public  :: krmax              (:)
     real(r8), allocatable, private :: kmax               (:,:)
     real(r8), allocatable, private :: psi50              (:,:)
     real(r8), allocatable, private :: ck                 (:,:)
     real(r8), allocatable, public  :: psi_soil_ref       (:)
     real(r8), allocatable, private :: lmr_intercept_atkin(:)
  contains
     procedure, private :: allocParams
  end type photo_params_type
  !
  type(photo_params_type), public, protected :: params_inst  ! params_inst is populated in readParamsMod 

  type, public :: photosyns_type

     logical , pointer, private :: c3flag_patch      (:)   ! patch true if C3 and false if C4
     ! Plant hydraulic stress specific variables
     real(r8), pointer, private :: ac_phs_patch      (:,:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: aj_phs_patch      (:,:,:) ! patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ap_phs_patch      (:,:,:) ! patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ag_phs_patch      (:,:,:) ! patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_sun_patch      (:,:)   ! patch sunlit net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_sha_patch      (:,:)   ! patch shaded net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: vcmax_z_phs_patch (:,:,:) ! patch maximum rate of carboxylation (umol co2/m**2/s)
     real(r8), pointer, private :: kp_z_phs_patch    (:,:,:) ! patch initial slope of CO2 response curve (C4 plants)
     real(r8), pointer, private :: tpu_z_phs_patch   (:,:,:) ! patch triose phosphate utilization rate (umol CO2/m**2/s)
     real(r8), pointer, private :: gs_mol_sun_patch  (:,:) ! patch sunlit leaf stomatal conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: gs_mol_sha_patch  (:,:) ! patch shaded leaf stomatal conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: gs_mol_sun_ln_patch (:,:) ! patch sunlit leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon (umol H2O/m**2/s)
     real(r8), pointer, private :: gs_mol_sha_ln_patch (:,:) ! patch shaded leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon (umol H2O/m**2/s)
     real(r8), pointer, private :: ac_patch          (:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: aj_patch          (:,:) ! patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ap_patch          (:,:) ! patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ag_patch          (:,:) ! patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_patch          (:,:) ! patch net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: vcmax_z_patch     (:,:) ! patch maximum rate of carboxylation (umol co2/m**2/s)
     real(r8), pointer, private :: cp_patch          (:)   ! patch CO2 compensation point (Pa)
     real(r8), pointer, private :: kc_patch          (:)   ! patch Michaelis-Menten constant for CO2 (Pa)
     real(r8), pointer, private :: ko_patch          (:)   ! patch Michaelis-Menten constant for O2 (Pa)
     real(r8), pointer, private :: qe_patch          (:)   ! patch quantum efficiency, used only for C4 (mol CO2 / mol photons)
     real(r8), pointer, private :: tpu_z_patch       (:,:) ! patch triose phosphate utilization rate (umol CO2/m**2/s)
     real(r8), pointer, private :: kp_z_patch        (:,:) ! patch initial slope of CO2 response curve (C4 plants)
     real(r8), pointer, private :: theta_cj_patch    (:)   ! patch empirical curvature parameter for ac, aj photosynthesis co-limitation
     real(r8), pointer, private :: bbb_patch         (:)   ! patch Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: mbb_patch         (:)   ! patch Ball-Berry slope of conductance-photosynthesis relationship
     real(r8), pointer, private :: gs_mol_patch      (:,:) ! patch leaf stomatal conductance       (umol H2O/m**2/s)
     real(r8), pointer, private :: gb_mol_patch      (:)   ! patch leaf boundary layer conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: rh_leaf_patch     (:)   ! patch fractional humidity at leaf surface (dimensionless)

     real(r8), pointer, private :: alphapsnsun_patch (:)   ! patch sunlit 13c fractionation ([])
     real(r8), pointer, private :: alphapsnsha_patch (:)   ! patch shaded 13c fractionation ([])

     real(r8), pointer, public  :: rc13_canair_patch (:)   ! patch C13O2/C12O2 in canopy air
     real(r8), pointer, public  :: rc13_psnsun_patch (:)   ! patch C13O2/C12O2 in sunlit canopy psn flux
     real(r8), pointer, public  :: rc13_psnsha_patch (:)   ! patch C13O2/C12O2 in shaded canopy psn flux

     real(r8), pointer, public  :: psnsun_patch      (:)   ! patch sunlit leaf photosynthesis     (umol CO2/m**2/s)
     real(r8), pointer, public  :: psnsha_patch      (:)   ! patch shaded leaf photosynthesis     (umol CO2/m**2/s)
     real(r8), pointer, public  :: c13_psnsun_patch  (:)   ! patch c13 sunlit leaf photosynthesis (umol 13CO2/m**2/s)
     real(r8), pointer, public  :: c13_psnsha_patch  (:)   ! patch c13 shaded leaf photosynthesis (umol 13CO2/m**2/s)
     real(r8), pointer, public  :: c14_psnsun_patch  (:)   ! patch c14 sunlit leaf photosynthesis (umol 14CO2/m**2/s)
     real(r8), pointer, public  :: c14_psnsha_patch  (:)   ! patch c14 shaded leaf photosynthesis (umol 14CO2/m**2/s)

     real(r8), pointer, private :: psnsun_z_patch    (:,:) ! patch canopy layer: sunlit leaf photosynthesis   (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_z_patch    (:,:) ! patch canopy layer: shaded leaf photosynthesis   (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wc_patch   (:)   ! patch Rubsico-limited sunlit leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wc_patch   (:)   ! patch Rubsico-limited shaded leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wj_patch   (:)   ! patch RuBP-limited sunlit leaf photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wj_patch   (:)   ! patch RuBP-limited shaded leaf photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wp_patch   (:)   ! patch product-limited sunlit leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wp_patch   (:)   ! patch product-limited shaded leaf photosynthesis (umol CO2/m**2/s)

     real(r8), pointer, public  :: fpsn_patch        (:)   ! patch photosynthesis                 (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wc_patch     (:)   ! patch Rubisco-limited photosynthesis (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wj_patch     (:)   ! patch RuBP-limited photosynthesis    (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wp_patch     (:)   ! patch product-limited photosynthesis (umol CO2/m**2 ground/s)

     real(r8), pointer, public  :: lnca_patch        (:)   ! top leaf layer leaf N concentration (gN leaf/m^2)

     real(r8), pointer, public  :: lmrsun_patch      (:)   ! patch sunlit leaf maintenance respiration rate               (umol CO2/m**2/s)
     real(r8), pointer, public  :: lmrsha_patch      (:)   ! patch shaded leaf maintenance respiration rate               (umol CO2/m**2/s)
     real(r8), pointer, private :: lmrsun_z_patch    (:,:) ! patch canopy layer: sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
     real(r8), pointer, private :: lmrsha_z_patch    (:,:) ! patch canopy layer: shaded leaf maintenance respiration rate (umol CO2/m**2/s)

     real(r8), pointer, public  :: cisun_z_patch     (:,:) ! patch intracellular sunlit leaf CO2 (Pa)
     real(r8), pointer, public  :: cisha_z_patch     (:,:) ! patch intracellular shaded leaf CO2 (Pa)

     real(r8), pointer, private :: rssun_z_patch     (:,:) ! patch canopy layer: sunlit leaf stomatal resistance (s/m)
     real(r8), pointer, private :: rssha_z_patch     (:,:) ! patch canopy layer: shaded leaf stomatal resistance (s/m)
     real(r8), pointer, public  :: rssun_patch       (:)   ! patch sunlit stomatal resistance (s/m)
     real(r8), pointer, public  :: rssha_patch       (:)   ! patch shaded stomatal resistance (s/m)
     real(r8), pointer, public  :: luvcmax25top_patch (:)   ! vcmax25 !     (umol/m2/s)
     real(r8), pointer, public  :: lujmax25top_patch  (:)   ! vcmax25 (umol/m2/s)
     real(r8), pointer, public  :: lutpu25top_patch   (:)   ! vcmax25 (umol/m2/s)
!!


     ! LUNA specific variables
     real(r8), pointer, public  :: vcmx25_z_patch    (:,:) ! patch  leaf Vc,max25 (umol CO2/m**2/s) for canopy layer 
     real(r8), pointer, public  :: jmx25_z_patch     (:,:) ! patch  leaf Jmax25 (umol electron/m**2/s) for canopy layer 
     real(r8), pointer, public  :: pnlc_z_patch      (:,:) ! patch proportion of leaf nitrogen allocated for light capture for canopy layer
     real(r8), pointer, public  :: enzs_z_patch      (:,:) ! enzyme decay status 1.0-fully active; 0-all decayed during stress
     real(r8), pointer, public  :: fpsn24_patch      (:)   ! 24 hour mean patch photosynthesis (umol CO2/m**2 ground/day)

     ! Logical switches for different options
     logical, public  :: rootstem_acc                      ! Respiratory acclimation for roots and stems
     logical, private :: light_inhibit                     ! If light should inhibit respiration
     integer, private :: leafresp_method                   ! leaf maintencence respiration at 25C for canopy top method to use
     integer, private :: stomatalcond_mtd                  ! Stomatal conduction method type
     logical, private :: modifyphoto_and_lmr_forcrop       ! Modify photosynthesis and LMR for crop
   contains

     ! Public procedures
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: ReadNML
     procedure, public  :: ReadParams
     procedure, public  :: TimeStepInit
     procedure, public  :: NewPatchInit

     ! Private procedures
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type photosyns_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%c3flag_patch      (begp:endp))             ; this%c3flag_patch      (:)     =.false.
    allocate(this%ac_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ac_phs_patch      (:,:,:) = nan
    allocate(this%aj_phs_patch      (begp:endp,2,1:nlevcan)) ; this%aj_phs_patch      (:,:,:) = nan
    allocate(this%ap_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ap_phs_patch      (:,:,:) = nan
    allocate(this%ag_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ag_phs_patch      (:,:,:) = nan
    allocate(this%an_sun_patch      (begp:endp,1:nlevcan))   ; this%an_sun_patch      (:,:)   = nan
    allocate(this%an_sha_patch      (begp:endp,1:nlevcan))   ; this%an_sha_patch      (:,:)   = nan
    allocate(this%vcmax_z_phs_patch (begp:endp,2,1:nlevcan)) ; this%vcmax_z_phs_patch (:,:,:) = nan
    allocate(this%tpu_z_phs_patch   (begp:endp,2,1:nlevcan)) ; this%tpu_z_phs_patch   (:,:,:) = nan
    allocate(this%kp_z_phs_patch    (begp:endp,2,1:nlevcan)) ; this%kp_z_phs_patch    (:,:,:) = nan
    allocate(this%gs_mol_sun_patch  (begp:endp,1:nlevcan))   ; this%gs_mol_sun_patch  (:,:)   = nan
    allocate(this%gs_mol_sha_patch  (begp:endp,1:nlevcan))   ; this%gs_mol_sha_patch  (:,:)   = nan
    allocate(this%gs_mol_sun_ln_patch (begp:endp,1:nlevcan)) ; this%gs_mol_sun_ln_patch (:,:)   = nan
    allocate(this%gs_mol_sha_ln_patch (begp:endp,1:nlevcan)) ; this%gs_mol_sha_ln_patch (:,:)   = nan
    allocate(this%ac_patch          (begp:endp,1:nlevcan)) ; this%ac_patch          (:,:) = nan
    allocate(this%aj_patch          (begp:endp,1:nlevcan)) ; this%aj_patch          (:,:) = nan
    allocate(this%ap_patch          (begp:endp,1:nlevcan)) ; this%ap_patch          (:,:) = nan
    allocate(this%ag_patch          (begp:endp,1:nlevcan)) ; this%ag_patch          (:,:) = nan
    allocate(this%an_patch          (begp:endp,1:nlevcan)) ; this%an_patch          (:,:) = nan
    allocate(this%vcmax_z_patch     (begp:endp,1:nlevcan)) ; this%vcmax_z_patch     (:,:) = nan
    allocate(this%tpu_z_patch       (begp:endp,1:nlevcan)) ; this%tpu_z_patch       (:,:) = nan
    allocate(this%kp_z_patch        (begp:endp,1:nlevcan)) ; this%kp_z_patch        (:,:) = nan
    allocate(this%gs_mol_patch      (begp:endp,1:nlevcan)) ; this%gs_mol_patch      (:,:) = nan
    allocate(this%cp_patch          (begp:endp))           ; this%cp_patch          (:)   = nan
    allocate(this%kc_patch          (begp:endp))           ; this%kc_patch          (:)   = nan
    allocate(this%ko_patch          (begp:endp))           ; this%ko_patch          (:)   = nan
    allocate(this%qe_patch          (begp:endp))           ; this%qe_patch          (:)   = nan
    allocate(this%theta_cj_patch    (begp:endp))           ; this%theta_cj_patch    (:)   = nan
    allocate(this%bbb_patch         (begp:endp))           ; this%bbb_patch         (:)   = nan
    allocate(this%mbb_patch         (begp:endp))           ; this%mbb_patch         (:)   = nan
    allocate(this%gb_mol_patch      (begp:endp))           ; this%gb_mol_patch      (:)   = nan
    allocate(this%rh_leaf_patch     (begp:endp))           ; this%rh_leaf_patch     (:)   = nan

    allocate(this%psnsun_patch      (begp:endp))           ; this%psnsun_patch      (:)   = nan
    allocate(this%psnsha_patch      (begp:endp))           ; this%psnsha_patch      (:)   = nan
    allocate(this%c13_psnsun_patch  (begp:endp))           ; this%c13_psnsun_patch  (:)   = nan
    allocate(this%c13_psnsha_patch  (begp:endp))           ; this%c13_psnsha_patch  (:)   = nan
    allocate(this%c14_psnsun_patch  (begp:endp))           ; this%c14_psnsun_patch  (:)   = nan
    allocate(this%c14_psnsha_patch  (begp:endp))           ; this%c14_psnsha_patch  (:)   = nan

    allocate(this%psnsun_z_patch    (begp:endp,1:nlevcan)) ; this%psnsun_z_patch    (:,:) = nan
    allocate(this%psnsha_z_patch    (begp:endp,1:nlevcan)) ; this%psnsha_z_patch    (:,:) = nan
    allocate(this%psnsun_wc_patch   (begp:endp))           ; this%psnsun_wc_patch   (:)   = nan
    allocate(this%psnsha_wc_patch   (begp:endp))           ; this%psnsha_wc_patch   (:)   = nan
    allocate(this%psnsun_wj_patch   (begp:endp))           ; this%psnsun_wj_patch   (:)   = nan
    allocate(this%psnsha_wj_patch   (begp:endp))           ; this%psnsha_wj_patch   (:)   = nan
    allocate(this%psnsun_wp_patch   (begp:endp))           ; this%psnsun_wp_patch   (:)   = nan
    allocate(this%psnsha_wp_patch   (begp:endp))           ; this%psnsha_wp_patch   (:)   = nan
    allocate(this%fpsn_patch        (begp:endp))           ; this%fpsn_patch        (:)   = nan
    allocate(this%fpsn_wc_patch     (begp:endp))           ; this%fpsn_wc_patch     (:)   = nan
    allocate(this%fpsn_wj_patch     (begp:endp))           ; this%fpsn_wj_patch     (:)   = nan
    allocate(this%fpsn_wp_patch     (begp:endp))           ; this%fpsn_wp_patch     (:)   = nan
    
    allocate(this%lnca_patch        (begp:endp))           ; this%lnca_patch        (:)   = nan

    allocate(this%lmrsun_z_patch    (begp:endp,1:nlevcan)) ; this%lmrsun_z_patch    (:,:) = nan
    allocate(this%lmrsha_z_patch    (begp:endp,1:nlevcan)) ; this%lmrsha_z_patch    (:,:) = nan
    allocate(this%lmrsun_patch      (begp:endp))           ; this%lmrsun_patch      (:)   = nan
    allocate(this%lmrsha_patch      (begp:endp))           ; this%lmrsha_patch      (:)   = nan

    allocate(this%alphapsnsun_patch (begp:endp))           ; this%alphapsnsun_patch (:)   = nan
    allocate(this%alphapsnsha_patch (begp:endp))           ; this%alphapsnsha_patch (:)   = nan
    allocate(this%rc13_canair_patch (begp:endp))           ; this%rc13_canair_patch (:)   = nan
    allocate(this%rc13_psnsun_patch (begp:endp))           ; this%rc13_psnsun_patch (:)   = nan
    allocate(this%rc13_psnsha_patch (begp:endp))           ; this%rc13_psnsha_patch (:)   = nan

    allocate(this%cisun_z_patch     (begp:endp,1:nlevcan)) ; this%cisun_z_patch     (:,:) = nan
    allocate(this%cisha_z_patch     (begp:endp,1:nlevcan)) ; this%cisha_z_patch     (:,:) = nan

    allocate(this%rssun_z_patch     (begp:endp,1:nlevcan)) ; this%rssun_z_patch     (:,:) = nan
    allocate(this%rssha_z_patch     (begp:endp,1:nlevcan)) ; this%rssha_z_patch     (:,:) = nan
    allocate(this%rssun_patch       (begp:endp))           ; this%rssun_patch       (:)   = nan
    allocate(this%rssha_patch       (begp:endp))           ; this%rssha_patch       (:)   = nan
    allocate(this%luvcmax25top_patch(begp:endp))           ; this%luvcmax25top_patch(:) = nan
    allocate(this%lujmax25top_patch (begp:endp))           ; this%lujmax25top_patch(:)  = nan
    allocate(this%lutpu25top_patch  (begp:endp))           ; this%lutpu25top_patch(:)   = nan
!!
!    allocate(this%psncanopy_patch   (begp:endp))           ; this%psncanopy_patch   (:)   = nan
!    allocate(this%lmrcanopy_patch   (begp:endp))           ; this%lmrcanopy_patch   (:)   = nan
    if(use_luna)then
      ! NOTE(bja, 2015-09) because these variables are only allocated
      ! when luna is turned on, they can not be placed into associate
      ! statements.
      allocate(this%vcmx25_z_patch  (begp:endp,1:nlevcan)) ; this%vcmx25_z_patch    (:,:) = 30._r8
      allocate(this%jmx25_z_patch   (begp:endp,1:nlevcan)) ; this%jmx25_z_patch     (:,:) = 60._r8 
      allocate(this%pnlc_z_patch    (begp:endp,1:nlevcan)) ; this%pnlc_z_patch      (:,:) = 0.01_r8
      allocate(this%fpsn24_patch    (begp:endp))           ; this%fpsn24_patch      (:)   = nan
      allocate(this%enzs_z_patch    (begp:endp,1:nlevcan)) ; this%enzs_z_patch      (:,:) = 1._r8
    endif

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), pointer  :: ptr_1d(:)  ! pointer to 1d patch array
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%rh_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname='RH_LEAF', units='fraction', &
         avgflag='A', long_name='fractional humidity at leaf surface', &
         ptr_patch=this%rh_leaf_patch, set_spec=spval, default='inactive')
    this%lnca_patch(begp:endp) = spval
    call hist_addfld1d (fname='LNC', units='gN leaf/m^2', &
         avgflag='A', long_name='leaf N concentration', &
         ptr_patch=this%lnca_patch, set_spec=spval)

    ! Don't output photosynthesis variables when FATES is on as they aren't calculated
    if (.not. use_fates) then
       this%fpsn_patch(begp:endp) = spval
       call hist_addfld1d (fname='FPSN', units='umol/m2s',  &
            avgflag='A', long_name='photosynthesis', &
            ptr_patch=this%fpsn_patch, set_lake=0._r8, set_urb=0._r8)

       ! Don't by default output this rate limiting step as only makes sense if you are outputing
       ! the others each time-step
       this%fpsn_wc_patch(begp:endp) = spval
       call hist_addfld1d (fname='FPSN_WC', units='umol/m2s',  &
            avgflag='I', long_name='Rubisco-limited photosynthesis', &
            ptr_patch=this%fpsn_wc_patch, set_lake=0._r8, set_urb=0._r8, &
            default='inactive')

       ! Don't by default output this rate limiting step as only makes sense if you are outputing
       ! the others each time-step
       this%fpsn_wj_patch(begp:endp) = spval
       call hist_addfld1d (fname='FPSN_WJ', units='umol/m2s',  &
            avgflag='I', long_name='RuBP-limited photosynthesis', &
            ptr_patch=this%fpsn_wj_patch, set_lake=0._r8, set_urb=0._r8, &
            default='inactive')

       ! Don't by default output this rate limiting step as only makes sense if you are outputing
       ! the others each time-step
       this%fpsn_wp_patch(begp:endp) = spval
       call hist_addfld1d (fname='FPSN_WP', units='umol/m2s',  &
            avgflag='I', long_name='Product-limited photosynthesis', &
            ptr_patch=this%fpsn_wp_patch, set_lake=0._r8, set_urb=0._r8, &
            default='inactive')
    end if

    if (use_cn) then
       this%psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='sunlit leaf photosynthesis', &
            ptr_patch=this%psnsun_patch)

       this%psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='shaded leaf photosynthesis', &
            ptr_patch=this%psnsha_patch)
    end if

    if ( use_c13 ) then
       this%c13_psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C13 sunlit leaf photosynthesis', &
            ptr_patch=this%c13_psnsun_patch)

       this%c13_psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C13 shaded leaf photosynthesis', &
            ptr_patch=this%c13_psnsha_patch)
    end if

    if ( use_c14 ) then
       this%c14_psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C14 sunlit leaf photosynthesis', &
            ptr_patch=this%c14_psnsun_patch)

       this%c14_psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C14 shaded leaf photosynthesis', &
            ptr_patch=this%c14_psnsha_patch)
    end if

    if ( use_c13 ) then
       this%rc13_canair_patch(begp:endp) = spval
       call hist_addfld1d (fname='RC13_CANAIR', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for canopy air', &
            ptr_patch=this%rc13_canair_patch)

       this%rc13_psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='RC13_PSNSUN', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for sunlit photosynthesis', &
            ptr_patch=this%rc13_psnsun_patch)

       this%rc13_psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='RC13_PSNSHA', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for shaded photosynthesis', &
            ptr_patch=this%rc13_psnsha_patch)
    endif

    ! Canopy physiology

    if ( use_c13 ) then
       this%alphapsnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='ALPHAPSNSUN', units='proportion', &
            avgflag='A', long_name='sunlit c13 fractionation', &
            ptr_patch=this%alphapsnsun_patch, default='inactive')

       this%alphapsnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='ALPHAPSNSHA', units='proportion', &
            avgflag='A', long_name='shaded c13 fractionation', &
            ptr_patch=this%alphapsnsha_patch, default='inactive')
    endif

    this%rssun_patch(begp:endp) = spval
    call hist_addfld1d (fname='RSSUN', units='s/m',  &
         avgflag='M', long_name='sunlit leaf stomatal resistance', &
         ptr_patch=this%rssun_patch, set_lake=spval, set_urb=spval)

    this%rssha_patch(begp:endp) = spval
    call hist_addfld1d (fname='RSSHA', units='s/m',  &
         avgflag='M', long_name='shaded leaf stomatal resistance', &
         ptr_patch=this%rssha_patch, set_lake=spval, set_urb=spval)

    this%gs_mol_sun_patch(begp:endp,:) = spval
    this%gs_mol_sha_patch(begp:endp,:) = spval
    if (nlevcan>1) then 
       call hist_addfld2d (fname='GSSUN', units='umol H20/m2/s', type2d='nlevcan', &
          avgflag='A', long_name='sunlit leaf stomatal conductance', &
          ptr_patch=this%gs_mol_sun_patch, set_lake=spval, set_urb=spval)

       call hist_addfld2d (fname='GSSHA', units='umol H20/m2/s', type2d='nlevcan', &
          avgflag='A', long_name='shaded leaf stomatal conductance', &
          ptr_patch=this%gs_mol_sha_patch, set_lake=spval, set_urb=spval)
    else
       ptr_1d => this%gs_mol_sun_patch(begp:endp,1)
       call hist_addfld1d (fname='GSSUN', units='umol H20/m2/s', &
          avgflag='A', long_name='sunlit leaf stomatal conductance', &
          ptr_patch=ptr_1d)

       ptr_1d => this%gs_mol_sha_patch(begp:endp,1)
       call hist_addfld1d (fname='GSSHA', units='umol H20/m2/s', &
          avgflag='A', long_name='shaded leaf stomatal conductance', &
          ptr_patch=ptr_1d)

    endif
    this%gs_mol_sun_ln_patch(begp:endp,:) = spval
    this%gs_mol_sha_ln_patch(begp:endp,:) = spval
    if (nlevcan>1) then
       call hist_addfld2d (fname='GSSUNLN', units='umol H20/m2/s', type2d='nlevcan', &
          avgflag='A', long_name='sunlit leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon', &
          ptr_patch=this%gs_mol_sun_ln_patch, set_lake=spval, set_urb=spval)

       call hist_addfld2d (fname='GSSHALN', units='umol H20/m2/s', type2d='nlevcan', &
          avgflag='A', long_name='shaded leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon', &
          ptr_patch=this%gs_mol_sha_ln_patch, set_lake=spval, set_urb=spval)
    else
       ptr_1d => this%gs_mol_sun_ln_patch(begp:endp,1)
       call hist_addfld1d (fname='GSSUNLN', units='umol H20/m2/s', &
          avgflag='A', long_name='sunlit leaf stomatal conductance at local noon', &
          ptr_patch=ptr_1d)

       ptr_1d => this%gs_mol_sha_ln_patch(begp:endp,1)
       call hist_addfld1d (fname='GSSHALN', units='umol H20/m2/s', &
          avgflag='A', long_name='shaded leaf stomatal conductance at local noon', &
          ptr_patch=ptr_1d)

    endif
    if(use_luna)then  
       if(nlevcan>1)then
         call hist_addfld2d (fname='Vcmx25Z', units='umol/m2/s', type2d='nlevcan', &
            avgflag='A', long_name='canopy profile of vcmax25 predicted by LUNA model', &
            ptr_patch=this%vcmx25_z_patch)
 
         call hist_addfld2d (fname='Jmx25Z', units='umol/m2/s', type2d='nlevcan', &
            avgflag='A', long_name='canopy profile of  vcmax25 predicted by LUNA model', &
            ptr_patch=this%jmx25_z_patch)

         call hist_addfld2d (fname='PNLCZ', units='unitless', type2d='nlevcan', &
            avgflag='A', long_name='Proportion of nitrogen allocated for light capture', &
            ptr_patch=this%pnlc_z_patch,default='inactive')
       else
         ptr_1d => this%vcmx25_z_patch(:,1)
         call hist_addfld1d (fname='Vcmx25Z', units='umol/m2/s',&
            avgflag='A', long_name='canopy profile of vcmax25 predicted by LUNA model', &
            ptr_patch=ptr_1d)
         ptr_1d => this%jmx25_z_patch(:,1)
         call hist_addfld1d (fname='Jmx25Z', units='umol/m2/s',&
            avgflag='A', long_name='canopy profile of  vcmax25 predicted by LUNA model', &
            ptr_patch=ptr_1d)
         ptr_1d => this%pnlc_z_patch(:,1)
         call hist_addfld1d (fname='PNLCZ', units='unitless', &
            avgflag='A', long_name='Proportion of nitrogen allocated for light capture', &
            ptr_patch=ptr_1d,default='inactive')

         this%luvcmax25top_patch(begp:endp) = spval
         call hist_addfld1d (fname='VCMX25T', units='umol/m2/s',  &
            avgflag='M', long_name='canopy profile of vcmax25', &
            ptr_patch=this%luvcmax25top_patch, set_lake=spval, set_urb=spval)

         this%lujmax25top_patch(begp:endp) = spval
         call hist_addfld1d (fname='JMX25T', units='umol/m2/s',  &
            avgflag='M', long_name='canopy profile of jmax', &
            ptr_patch=this%lujmax25top_patch, set_lake=spval, set_urb=spval)

            this%lutpu25top_patch(begp:endp) = spval
            call hist_addfld1d (fname='TPU25T', units='umol/m2/s',  &
            avgflag='M', long_name='canopy profile of tpu', &
            ptr_patch=this%lutpu25top_patch, set_lake=spval, set_urb=spval)

       endif
       this%fpsn24_patch = spval 
       call hist_addfld1d (fname='FPSN24', units='umol CO2/m**2 ground/day',&
           avgflag='A', long_name='24 hour accumulative patch photosynthesis starting from mid-night', &
           ptr_patch=this%fpsn24_patch, default='inactive')
   
    endif

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,l                        ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)

       this%alphapsnsun_patch(p) = spval
       this%alphapsnsha_patch(p) = spval

       if (lun%ifspecial(l)) then
          this%psnsun_patch(p) = 0._r8
          this%psnsha_patch(p) = 0._r8
          if ( use_c13 ) then
             this%c13_psnsun_patch(p) = 0._r8
             this%c13_psnsha_patch(p) = 0._r8
          endif
          if ( use_c14 ) then
             this%c14_psnsun_patch(p) = 0._r8
             this%c14_psnsha_patch(p) = 0._r8
          endif
       end if
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine allocParams ( this )
    !
    implicit none

    ! !ARGUMENTS:
    class(photo_params_type) :: this
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'allocParams'
    !-----------------------------------------------------------------------

    ! allocate parameters

    allocate( this%krmax       (0:mxpft) )          ; this%krmax(:)        = nan
    allocate( this%kmax        (0:mxpft,nvegwcs) )  ; this%kmax(:,:)       = nan
    allocate( this%psi50       (0:mxpft,nvegwcs) )  ; this%psi50(:,:)      = nan
    allocate( this%ck          (0:mxpft,nvegwcs) )  ; this%ck(:,:)         = nan
    allocate( this%psi_soil_ref(0:mxpft) )          ; this%psi_soil_ref(:) = nan

    if ( use_hydrstress .and. nvegwcs /= 4 )then
       call endrun(msg='Error:: the Plant Hydraulics Stress methodology is for the spacA function is hardcoded for nvegwcs==4' &
                   //errMsg(__FILE__, __LINE__))
    end if

  end subroutine allocParams

  !-----------------------------------------------------------------------
  subroutine readParams ( this, ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io
    implicit none

    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'readParams'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: temp1d(0:mxpft) ! temporary to read in parameter
    real(r8)           :: temp2d(0:mxpft,nvegwcs) ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! read in parameters


    call params_inst%allocParams()

    tString = "krmax"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%krmax=temp1d
    tString = "psi_soil_ref"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%psi_soil_ref=temp1d
    tString = "lmr_intercept_atkin"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%lmr_intercept_atkin=temp1d
    tString = "kmax"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%kmax=temp2d
    tString = "psi50"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%psi50=temp2d
    tString = "ck"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%ck=temp2d

  end subroutine readParams


  !------------------------------------------------------------------------
  subroutine ReadNML(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for Photosynthesis
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'Photosyn::ReadNML'
    character(len=*), parameter :: nmlname = 'photosyns_inparm'
    logical :: rootstem_acc    = .false.                     ! Respiratory acclimation for roots and stems
    logical :: light_inhibit   = .false.                     ! If light should inhibit respiration
    integer :: leafresp_method = leafresp_mtd_ryan1991       ! leaf maintencence respiration at 25C for canopy top method to use
    logical :: modifyphoto_and_lmr_forcrop = .false.            ! Modify photosynthesis and LMR for crop
    character(len=50) :: stomatalcond_method = 'Ball-Berry1987' ! Photosynthesis method string
    !-----------------------------------------------------------------------

    namelist /photosyns_inparm/ leafresp_method, light_inhibit, &
              rootstem_acc, stomatalcond_method, modifyphoto_and_lmr_forcrop

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=photosyns_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
       this%rootstem_acc    = rootstem_acc
       this%leafresp_method = leafresp_method
       this%light_inhibit   = light_inhibit
       this%modifyphoto_and_lmr_forcrop = modifyphoto_and_lmr_forcrop
       if (      trim(stomatalcond_method) == 'Ball-Berry1987' ) then
          this%stomatalcond_mtd = stomatalcond_mtd_bb1987
       else if ( trim(stomatalcond_method) == 'Medlyn2011'     ) then
          this%stomatalcond_mtd = stomatalcond_mtd_medlyn2011
       else
          call endrun(msg="ERROR bad value for stomtalcond_method in "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
    end if

    call shr_mpi_bcast (this%rootstem_acc   , mpicom)
    call shr_mpi_bcast (this%leafresp_method, mpicom)
    call shr_mpi_bcast (this%light_inhibit  , mpicom)
    call shr_mpi_bcast (this%stomatalcond_mtd, mpicom)
    call shr_mpi_bcast (this%modifyphoto_and_lmr_forcrop, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=photosyns_inparm)
       write(iulog,*) ' '
    end if

  end subroutine ReadNML

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    if ( use_c13 ) then
       call restartvar(ncid=ncid, flag=flag, varname='rc13_canair', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%rc13_canair_patch)

       call restartvar(ncid=ncid, flag=flag, varname='rc13_psnsun', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%rc13_psnsun_patch)

       call restartvar(ncid=ncid, flag=flag, varname='rc13_psnsha', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%rc13_psnsha_patch)
    endif

    call restartvar(ncid=ncid, flag=flag, varname='GSSUN', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='sunlit leaf stomatal conductance', units='umol H20/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%gs_mol_sun_patch)
    
    call restartvar(ncid=ncid, flag=flag, varname='GSSHA', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='shaded leaf stomatal conductance', units='umol H20/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%gs_mol_sha_patch)

    call restartvar(ncid=ncid, flag=flag, varname='GSSUNLN', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='sunlit leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon', &
         units='umol H20/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%gs_mol_sun_ln_patch)

    call restartvar(ncid=ncid, flag=flag, varname='GSSHALN', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='shaded leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon', &
         units='umol H20/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%gs_mol_sha_ln_patch)
    
    call restartvar(ncid=ncid, flag=flag, varname='lnca', xtype=ncd_double,  &
       dim1name='pft', long_name='leaf N concentration', units='gN leaf/m^2', &
       interpinic_flag='interp', readvar=readvar, data=this%lnca_patch)

    if(use_luna) then
      call restartvar(ncid=ncid, flag=flag, varname='vcmx25_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='Maximum carboxylation rate at 25 celcius for canopy layers', units='umol CO2/m**2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%vcmx25_z_patch)
      call restartvar(ncid=ncid, flag=flag, varname='jmx25_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='Maximum carboxylation rate at 25 celcius for canopy layers', units='umol CO2/m**2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%jmx25_z_patch)
      call restartvar(ncid=ncid, flag=flag, varname='pnlc_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='proportion of leaf nitrogen allocated for light capture', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%pnlc_z_patch )
      call restartvar(ncid=ncid, flag=flag, varname='enzs_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='enzyme decay status during stress: 1.0-fully active; 0.0-all decayed', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%enzs_z_patch )
      call restartvar(ncid=ncid, flag=flag, varname='gpp24', xtype=ncd_double,  &
            dim1name='pft', long_name='accumulative gross primary production', units='umol CO2/m**2 ground/day', &
            interpinic_flag='interp', readvar=readvar, data=this%fpsn24_patch)    
   endif
   call restartvar(ncid=ncid, flag=flag, varname='vcmx25t', xtype=ncd_double,  &
         dim1name='pft', long_name='canopy profile of vcmax25', &
         units='umol/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%luvcmax25top_patch)    

   call restartvar(ncid=ncid, flag=flag, varname='jmx25t', xtype=ncd_double,  &
         dim1name='pft', long_name='canopy profile of jmax', &
         units='umol/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%lujmax25top_patch)    

   call restartvar(ncid=ncid, flag=flag, varname='tpu25t', xtype=ncd_double,  &
         dim1name='pft', long_name='canopy profile of tpu', &
         units='umol/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%lutpu25top_patch)    

  end subroutine Restart

  !------------------------------------------------------------------------------
  subroutine TimeStepInit (this, bounds)
    !
    ! Time step initialization
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop, istice_mec, istwet
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type) , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,l ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       l = patch%landunit(p)
       if (.not. lun%lakpoi(l)) then
          this%psnsun_patch(p)    = 0._r8
          this%psnsun_wc_patch(p) = 0._r8
          this%psnsun_wj_patch(p) = 0._r8
          this%psnsun_wp_patch(p) = 0._r8

          this%psnsha_patch(p)    = 0._r8
          this%psnsha_wc_patch(p) = 0._r8
          this%psnsha_wj_patch(p) = 0._r8
          this%psnsha_wp_patch(p) = 0._r8

          this%fpsn_patch(p)      = 0._r8
          this%fpsn_wc_patch(p)   = 0._r8
          this%fpsn_wj_patch(p)   = 0._r8
          this%fpsn_wp_patch(p)   = 0._r8

          if ( use_c13 ) then
             this%alphapsnsun_patch(p) = 0._r8
             this%alphapsnsha_patch(p) = 0._r8
             this%c13_psnsun_patch(p)  = 0._r8
             this%c13_psnsha_patch(p)  = 0._r8
          endif
          if ( use_c14 ) then
             this%c14_psnsun_patch(p) = 0._r8
             this%c14_psnsha_patch(p) = 0._r8
          endif
       end if
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop &
            .or. lun%itype(l) == istice_mec &
            .or. lun%itype(l) == istwet) then
          if (use_c13) then
             this%rc13_canair_patch(p) = 0._r8
             this%rc13_psnsun_patch(p) = 0._r8
             this%rc13_psnsha_patch(p) = 0._r8
          end if
       end if
    end do

  end subroutine TimeStepInit

  !------------------------------------------------------------------------------
  subroutine NewPatchInit (this, p)
    !
    ! For new run-time pft, modify state and flux variables to maintain
    ! carbon and nitrogen balance with dynamic pft-weights.
    ! Called from dyn_cnbal_patch
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    integer, intent(in) :: p
    !-----------------------------------------------------------------------

    if ( use_c13 ) then
       this%alphapsnsun_patch(p) = 0._r8
       this%alphapsnsha_patch(p) = 0._r8
       this%rc13_canair_patch(p) = 0._r8
       this%rc13_psnsun_patch(p) = 0._r8
       this%rc13_psnsha_patch(p) = 0._r8
    endif

    this%psnsun_patch(p) = 0._r8
    this%psnsha_patch(p) = 0._r8

    if (use_c13) then
       this%c13_psnsun_patch(p) = 0._r8
       this%c13_psnsha_patch(p) = 0._r8
    end if
    if ( use_c14 ) then
       this%c14_psnsun_patch(p) = 0._r8
       this%c14_psnsha_patch(p) = 0._r8
    end if

  end subroutine NewPatchInit

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine Photosynthesis ( bounds, fn, filterp, &
       esat_tv, eair, oair, cair, rb, btran, &
       dayl_factor, leafn, &
       atm2lnd_inst, temperature_inst, surfalb_inst, solarabs_inst, &
       canopystate_inst, ozone_inst, photosyns_inst, phase)
    !
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance calculation as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy
    !
    ! !USES:
    use clm_varcon        , only : rgas, tfrz, spval, degpsec, isecspday
    use GridcellType      , only : grc
    use clm_time_manager  , only : get_curr_date, get_step_size
    use clm_varctl     , only : cnallocate_carbon_only
    use clm_varctl     , only : lnc_opt, reduce_dayl_factor, vcmax_opt    
    use pftconMod      , only : nbrdlf_dcd_tmp_shrub, npcropmin

    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: fn                             ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                    ! patch filter
    real(r8)               , intent(in)    :: esat_tv( bounds%begp: )        ! saturation vapor pressure at t_veg (Pa) [pft]
    real(r8)               , intent(in)    :: eair( bounds%begp: )           ! vapor pressure of canopy air (Pa) [pft]
    real(r8)               , intent(in)    :: oair( bounds%begp: )           ! Atmospheric O2 partial pressure (Pa) [pft]
    real(r8)               , intent(in)    :: cair( bounds%begp: )           ! Atmospheric CO2 partial pressure (Pa) [pft]
    real(r8)               , intent(in)    :: rb( bounds%begp: )             ! boundary layer resistance (s/m) [pft]
    real(r8)               , intent(in)    :: btran( bounds%begp: )          ! transpiration wetness factor (0 to 1) [pft]
    real(r8)               , intent(in)    :: dayl_factor( bounds%begp: )    ! scalar (0-1) for daylength
    real(r8)               , intent(in)    :: leafn( bounds%begp: )          ! leaf N (gN/m2)
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(solarabs_type)    , intent(in)    :: solarabs_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    class(ozone_base_type) , intent(in)    :: ozone_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    character(len=*)       , intent(in)    :: phase                          ! 'sun' or 'sha'

    !
    ! !LOCAL VARIABLES:
    !
    ! Leaf photosynthesis parameters
    real(r8) :: jmax_z(bounds%begp:bounds%endp,nlevcan)  ! maximum electron transport rate (umol electrons/m**2/s)
    !real(r8) :: lnc(bounds%begp:bounds%endp)   ! leaf N concentration (gN leaf/m^2)
    real(r8) :: bbbopt(bounds%begp:bounds%endp)! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
    real(r8) :: kn(bounds%begp:bounds%endp)    ! leaf nitrogen decay coefficient
    real(r8) :: vcmax25top     ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25top      ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25top       ! canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top       ! canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25top        ! canopy top: initial slope of CO2 response curve (C4 plants) at 25C

    real(r8) :: vcmax25        ! leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25         ! leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25          ! leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25          ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25           ! leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
    real(r8) :: kc25           ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(r8) :: ko25           ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(r8) :: cp25           ! CO2 compensation point at 25C (Pa)

    real(r8) :: vcmaxha        ! activation energy for vcmax (J/mol)
    real(r8) :: jmaxha         ! activation energy for jmax (J/mol)
    real(r8) :: tpuha          ! activation energy for tpu (J/mol)
    real(r8) :: lmrha          ! activation energy for lmr (J/mol)
    real(r8) :: kcha           ! activation energy for kc (J/mol)
    real(r8) :: koha           ! activation energy for ko (J/mol)
    real(r8) :: cpha           ! activation energy for cp (J/mol)

    real(r8) :: vcmaxhd        ! deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd         ! deactivation energy for jmax (J/mol)
    real(r8) :: tpuhd          ! deactivation energy for tpu (J/mol)
    real(r8) :: lmrhd          ! deactivation energy for lmr (J/mol)

    real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K)
    real(r8) :: tpuse          ! entropy term for tpu (J/mol/K)
    real(r8) :: lmrse          ! entropy term for lmr (J/mol/K)

    real(r8) :: vcmaxc         ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc          ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: tpuc           ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: lmrc           ! scaling factor for high temperature inhibition (25 C = 1.0)

    real(r8) :: fnps           ! fraction of light absorbed by non-photosynthetic pigments
    real(r8) :: theta_psii     ! empirical curvature parameter for electron transport rate

    real(r8) :: theta_ip          ! empirical curvature parameter for ap photosynthesis co-limitation

    ! Other
    integer  :: f,p,c,iv          ! indices
    real(r8) :: cf                ! s m**2/umol -> s/m
    real(r8) :: rsmax0            ! maximum stomatal resistance [s/m]
    real(r8) :: gb                ! leaf boundary layer conductance (m/s)
    real(r8) :: cs                ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: gs                ! leaf stomatal conductance (m/s)
    real(r8) :: hs                ! fractional humidity at leaf surface (dimensionless)
    real(r8) :: sco               ! relative specificity of rubisco
    real(r8) :: ft                ! photosynthesis temperature response (statement function)
    real(r8) :: fth               ! photosynthesis temperature inhibition (statement function)
    real(r8) :: fth25             ! ccaling factor for photosynthesis temperature inhibition (statement function)
    real(r8) :: tl                ! leaf temperature in photosynthesis temperature function (K)
    real(r8) :: ha                ! activation energy in photosynthesis temperature function (J/mol)
    real(r8) :: hd                ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8) :: se                ! entropy term in photosynthesis temperature function (J/mol/K)
    real(r8) :: scaleFactor       ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: ciold             ! previous value of Ci for convergence check
    real(r8) :: gs_mol_err        ! gs_mol for error check
    real(r8) :: je                ! electron transport rate (umol electrons/m**2/s)
    real(r8) :: qabs              ! PAR absorbed by PS II (umol photons/m**2/s)
    real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
    real(r8) :: r1,r2             ! roots of quadratic equation
    real(r8) :: ceair             ! vapor pressure of air, constrained (Pa)
    real(r8) :: fnr               ! (gRubisco/gN in Rubisco)
    real(r8) :: act25             ! (umol/mgRubisco/min) Rubisco activity at 25 C
    integer  :: niter             ! iteration loop index
    real(r8) :: nscaler           ! leaf nitrogen scaling coefficient

    real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)

    real(r8) :: psn_wc_z(bounds%begp:bounds%endp,nlevcan) ! Rubisco-limited contribution to psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wj_z(bounds%begp:bounds%endp,nlevcan) ! RuBP-limited contribution to psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wp_z(bounds%begp:bounds%endp,nlevcan) ! product-limited contribution to psn_z (umol CO2/m**2/s)

    real(r8) :: psncan            ! canopy sum of psn_z
    real(r8) :: psncan_wc         ! canopy sum of psn_wc_z
    real(r8) :: psncan_wj         ! canopy sum of psn_wj_z
    real(r8) :: psncan_wp         ! canopy sum of psn_wp_z
    real(r8) :: lmrcan            ! canopy sum of lmr_z
    real(r8) :: gscan             ! canopy sum of leaf conductance
    real(r8) :: laican            ! canopy sum of lai_z
    real(r8) :: rh_can
    real(r8) , pointer :: lai_z       (:,:)
    real(r8) , pointer :: par_z       (:,:)
    real(r8) , pointer :: vcmaxcint   (:)
    real(r8) , pointer :: alphapsn    (:)
    real(r8) , pointer :: psn         (:)
    real(r8) , pointer :: psn_wc      (:)
    real(r8) , pointer :: psn_wj      (:)
    real(r8) , pointer :: psn_wp      (:)
    real(r8) , pointer :: psn_z       (:,:)
    real(r8) , pointer :: lmr         (:)
    real(r8) , pointer :: lmr_z       (:,:)
    real(r8) , pointer :: rs          (:)
    real(r8) , pointer :: rs_z        (:,:)
    real(r8) , pointer :: ci_z        (:,:)
    real(r8) , pointer :: o3coefv     (:)  ! o3 coefficient used in photo calculation
    real(r8) , pointer :: o3coefg     (:)  ! o3 coefficient used in rs calculation
    real(r8) , pointer :: alphapsnsun (:)
    real(r8) , pointer :: alphapsnsha (:)

    real(r8) :: sum_nscaler              
    real(r8) :: total_lai                
    integer  :: nptreemax                

    integer  :: local_secp1                     ! seconds into current date in local time
    real(r8) :: dtime                           ! land model time step (sec)
    integer  :: year,month,day,secs             ! calendar info for current time step
    integer  :: g                               ! index
    integer, parameter :: noonsec = isecspday / 2 ! seconds at local noon
    !------------------------------------------------------------------------------

    ! Temperature and soil water response functions

    ft(tl,ha) = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )
    fth(tl,hd,se,scaleFactor) = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )
    fth25(hd,se) = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    ! Enforce expected array sizes

    SHR_ASSERT_ALL((ubound(esat_tv)     == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(eair)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(oair)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(cair)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(rb)          == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(btran)       == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dayl_factor) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(leafn)       == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    associate(                                                 &
         c3psn      => pftcon%c3psn                          , & ! Input:  photosynthetic pathway: 0. = c4, 1. = c3
	 crop       => pftcon%crop                           , & ! Input:  crop or not (0 =not crop and 1 = crop)
         leafcn     => pftcon%leafcn                         , & ! Input:  leaf C:N (gC/gN)
         flnr       => pftcon%flnr                           , & ! Input:  fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
         fnitr      => pftcon%fnitr                          , & ! Input:  foliage nitrogen limitation factor (-)
         slatop     => pftcon%slatop                         , & ! Input:  specific leaf area at top of canopy, projected area basis [m^2/gC]
         dsladlai   => pftcon%dsladlai                       , & ! Input:  change in sla per unit lai  
         i_vcad     => pftcon%i_vcad                         , & ! Input:  [real(r8) (:)   ]  
         s_vcad     => pftcon%s_vcad                         , & ! Input:  [real(r8) (:)   ]  
         i_flnr     => pftcon%i_flnr                         , & ! Input:  [real(r8) (:)   ]  
         s_flnr     => pftcon%s_flnr                         , & ! Input:  [real(r8) (:)   ]  
         mbbopt     => pftcon%mbbopt                         , & ! Input:  [real(r8) (:)   ]  Ball-Berry slope of conduct/photosyn (umol H2O/umol CO2)
         ivt        => patch%itype                           , & ! Input:  [integer  (:)   ]  patch vegetation type
         forc_pbot  => atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)

         t_veg      => temperature_inst%t_veg_patch          , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)
         t10        => temperature_inst%t_a10_patch          , & ! Input:  [real(r8) (:)   ]  10-day running mean of the 2 m temperature (K)
         tgcm       => temperature_inst%thm_patch            , & ! Input:  [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)

         nrad       => surfalb_inst%nrad_patch               , & ! Input:  [integer  (:)   ]  pft number of canopy layers, above snow for radiative transfer
         tlai_z     => surfalb_inst%tlai_z_patch             , & ! Input:  [real(r8) (:,:) ]  pft total leaf area index for canopy layer
         tlai       => canopystate_inst%tlai_patch           , & ! Input:  [real(r8)(:)    ]  one-sided leaf area index, no burying by snow  
         c3flag     => photosyns_inst%c3flag_patch           , & ! Output: [logical  (:)   ]  true if C3 and false if C4
         ac         => photosyns_inst%ac_patch               , & ! Output: [real(r8) (:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
         aj         => photosyns_inst%aj_patch               , & ! Output: [real(r8) (:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)
         ap         => photosyns_inst%ap_patch               , & ! Output: [real(r8) (:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag         => photosyns_inst%ag_patch               , & ! Output: [real(r8) (:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)
         an         => photosyns_inst%an_patch               , & ! Output: [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)
         gb_mol     => photosyns_inst%gb_mol_patch           , & ! Output: [real(r8) (:)   ]  leaf boundary layer conductance (umol H2O/m**2/s)
         gs_mol     => photosyns_inst%gs_mol_patch           , & ! Output: [real(r8) (:,:) ]  leaf stomatal conductance (umol H2O/m**2/s)
         gs_mol_sun_ln => photosyns_inst%gs_mol_sun_ln_patch , & ! Output: [real(r8) (:,:) ]  sunlit leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon (umol H2O/m**2/s)
         gs_mol_sha_ln => photosyns_inst%gs_mol_sha_ln_patch , & ! Output: [real(r8) (:,:) ]  shaded leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon (umol H2O/m**2/s)
         vcmax_z    => photosyns_inst%vcmax_z_patch          , & ! Output: [real(r8) (:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         cp         => photosyns_inst%cp_patch               , & ! Output: [real(r8) (:)   ]  CO2 compensation point (Pa)
         kc         => photosyns_inst%kc_patch               , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for CO2 (Pa)
         ko         => photosyns_inst%ko_patch               , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for O2 (Pa)
         qe         => photosyns_inst%qe_patch               , & ! Output: [real(r8) (:)   ]  quantum efficiency, used only for C4 (mol CO2 / mol photons)
         tpu_z      => photosyns_inst%tpu_z_patch            , & ! Output: [real(r8) (:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)
         kp_z       => photosyns_inst%kp_z_patch             , & ! Output: [real(r8) (:,:) ]  initial slope of CO2 response curve (C4 plants)
         theta_cj   => photosyns_inst%theta_cj_patch         , & ! Output: [real(r8) (:)   ]  empirical curvature parameter for ac, aj photosynthesis co-limitation
         bbb        => photosyns_inst%bbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
         mbb        => photosyns_inst%mbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship
         rh_leaf    => photosyns_inst%rh_leaf_patch          , & ! Output: [real(r8) (:)   ]  fractional humidity at leaf surface (dimensionless)
         lnc        => photosyns_inst%lnca_patch             , & ! Output: [real(r8) (:)   ]  top leaf layer leaf N concentration (gN leaf/m^2)
         light_inhibit=> photosyns_inst%light_inhibit        , & ! Input:  [logical        ]  flag if light should inhibit respiration
         leafresp_method=> photosyns_inst%leafresp_method    , & ! Input:  [integer        ]  method type to use for leaf-maint.-respiration at 25C canopy top
         stomatalcond_mtd=> photosyns_inst%stomatalcond_mtd  , & ! Input:  [integer        ]  method type to use for stomatal conductance.GC.fnlprmsn15_r22845
         leaf_mr_vcm => canopystate_inst%leaf_mr_vcm           & ! Input:  [real(r8)       ]  scalar constant of leaf respiration with Vcmax
         )

      if (phase == 'sun') then
         par_z     =>    solarabs_inst%parsun_z_patch        ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)
         lai_z     =>    canopystate_inst%laisun_z_patch     ! Input:  [real(r8) (:,:) ]  leaf area index for canopy layer, sunlit or shaded
         vcmaxcint =>    surfalb_inst%vcmaxcintsun_patch     ! Input:  [real(r8) (:)   ]  leaf to canopy scaling coefficient
         alphapsn  =>    photosyns_inst%alphapsnsun_patch    ! Input:  [real(r8) (:)   ]  13C fractionation factor for PSN ()
         o3coefv   =>    ozone_inst%o3coefvsun_patch         ! Input:  [real(r8) (:)   ]  O3 coefficient used in photosynthesis calculation
         o3coefg   =>    ozone_inst%o3coefgsun_patch         ! Input:  [real(r8) (:)   ]  O3 coefficient used in rs calculation
         ci_z      =>    photosyns_inst%cisun_z_patch        ! Output: [real(r8) (:,:) ]  intracellular leaf CO2 (Pa)
         rs        =>    photosyns_inst%rssun_patch          ! Output: [real(r8) (:)   ]  leaf stomatal resistance (s/m)
         rs_z      =>    photosyns_inst%rssun_z_patch        ! Output: [real(r8) (:,:) ]  canopy layer: leaf stomatal resistance (s/m)
         lmr       =>    photosyns_inst%lmrsun_patch         ! Output: [real(r8) (:)   ]  leaf maintenance respiration rate (umol CO2/m**2/s)
         lmr_z     =>    photosyns_inst%lmrsun_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
         psn       =>    photosyns_inst%psnsun_patch         ! Output: [real(r8) (:)   ]  foliage photosynthesis (umol co2 /m**2/ s) [always +]
         psn_z     =>    photosyns_inst%psnsun_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
         psn_wc    =>    photosyns_inst%psnsun_wc_patch      ! Output: [real(r8) (:)   ]  Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
         psn_wj    =>    photosyns_inst%psnsun_wj_patch      ! Output: [real(r8) (:)   ]  RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
         psn_wp    =>    photosyns_inst%psnsun_wp_patch      ! Output: [real(r8) (:)   ]  product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      else if (phase == 'sha') then
         par_z     =>    solarabs_inst%parsha_z_patch        ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)
         lai_z     =>    canopystate_inst%laisha_z_patch     ! Input:  [real(r8) (:,:) ]  leaf area index for canopy layer, sunlit or shaded
         vcmaxcint =>    surfalb_inst%vcmaxcintsha_patch     ! Input:  [real(r8) (:)   ]  leaf to canopy scaling coefficient
         alphapsn  =>    photosyns_inst%alphapsnsha_patch    ! Input:  [real(r8) (:)   ]  13C fractionation factor for PSN ()
         o3coefv   =>    ozone_inst%o3coefvsha_patch         ! Input:  [real(r8) (:)   ]  O3 coefficient used in photosynthesis calculation
         o3coefg   =>    ozone_inst%o3coefgsha_patch         ! Input:  [real(r8) (:)   ]  O3 coefficient used in rs calculation
         ci_z      =>    photosyns_inst%cisha_z_patch        ! Output: [real(r8) (:,:) ]  intracellular leaf CO2 (Pa)
         rs        =>    photosyns_inst%rssha_patch          ! Output: [real(r8) (:)   ]  leaf stomatal resistance (s/m)
         rs_z      =>    photosyns_inst%rssha_z_patch        ! Output: [real(r8) (:,:) ]  canopy layer: leaf stomatal resistance (s/m)
         lmr       =>    photosyns_inst%lmrsha_patch         ! Output: [real(r8) (:)   ]  leaf maintenance respiration rate (umol CO2/m**2/s)
         lmr_z     =>    photosyns_inst%lmrsha_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
         psn       =>    photosyns_inst%psnsha_patch         ! Output: [real(r8) (:)   ]  foliage photosynthesis (umol co2 /m**2/ s) [always +]
         psn_z     =>    photosyns_inst%psnsha_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
         psn_wc    =>    photosyns_inst%psnsha_wc_patch      ! Output: [real(r8) (:)   ]  Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
         psn_wj    =>    photosyns_inst%psnsha_wj_patch      ! Output: [real(r8) (:)   ]  RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
         psn_wp    =>    photosyns_inst%psnsha_wp_patch      ! Output: [real(r8) (:)   ]  product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      end if

      !==============================================================================!
      ! Photosynthesis and stomatal conductance parameters, from:
      ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      !==============================================================================!

      ! Determine seconds of current time step

      dtime = get_step_size()
      call get_curr_date (year, month, day, secs)

      ! vcmax25 parameters, from CN

      fnr = 7.16_r8
      act25 = 3.6_r8   !umol/mgRubisco/min
      ! Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
      act25 = act25 * 1000.0_r8 / 60.0_r8

      ! Activation energy, from:
      ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
      ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282

      kcha    = 79430._r8
      koha    = 36380._r8
      cpha    = 37830._r8
      vcmaxha = 72000._r8
      jmaxha  = 50000._r8
      tpuha   = 72000._r8
      lmrha   = 46390._r8

      ! High temperature deactivation, from:
      ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
      ! The factor "c" scales the deactivation to a value of 1.0 at 25C

      vcmaxhd = 200000._r8
      jmaxhd  = 200000._r8
      tpuhd   = 200000._r8
      lmrhd   = 150650._r8
      lmrse   = 490._r8
      lmrc    = fth25 (lmrhd, lmrse)

      ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593

      fnps = 0.15_r8
      theta_psii = 0.7_r8
      theta_ip = 0.95_r8

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)

         ! C3 or C4 photosynthesis logical variable

         if (nint(c3psn(patch%itype(p))) == 1) then
            c3flag(p) = .true.
         else if (nint(c3psn(patch%itype(p))) == 0) then
            c3flag(p) = .false.
         end if

         ! C3 and C4 dependent parameters

         if (c3flag(p)) then
            qe(p) = 0._r8
            theta_cj(p) = 0.98_r8
            bbbopt(p) = 10000._r8
         else
            qe(p) = 0.05_r8
            theta_cj(p) = 0.80_r8
            bbbopt(p) = 40000._r8
         end if

         ! Soil water stress applied to Ball-Berry parameters

         bbb(p) = max (bbbopt(p)*btran(p), 1._r8)
         mbb(p) = mbbopt(patch%itype(p))

         ! kc, ko, cp, from: Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
         !
         !       kc25 = 404.9 umol/mol
         !       ko25 = 278.4 mmol/mol
         !       cp25 = 42.75 umol/mol
         !
         ! Derive sco from cp and O2 using present-day O2 (0.209 mol/mol) and re-calculate
         ! cp to account for variation in O2 using cp = 0.5 O2 / sco
         !

         kc25 = (404.9_r8 / 1.e06_r8) * forc_pbot(c)
         ko25 = (278.4_r8 / 1.e03_r8) * forc_pbot(c)
         sco  = 0.5_r8 * 0.209_r8 / (42.75_r8 / 1.e06_r8)
         cp25 = 0.5_r8 * oair(p) / sco

         kc(p) = kc25 * ft(t_veg(p), kcha)
         ko(p) = ko25 * ft(t_veg(p), koha)
         cp(p) = cp25 * ft(t_veg(p), cpha)

      end do

      ! Multi-layer parameters scaled by leaf nitrogen profile.
      ! Loop through each canopy layer to calculate nitrogen profile using
      ! cumulative lai at the midpoint of the layer

      do f = 1, fn
         p = filterp(f)

         if (lnc_opt .eqv. .false.) then     
            ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
            
           if ( (slatop(patch%itype(p)) *leafcn(patch%itype(p))) .le. 0.0_r8)then
              call endrun( "ERROR: slatop or leafcn is zero" )
           end if
           lnc(p) = 1._r8 / (slatop(patch%itype(p)) * leafcn(patch%itype(p)))
         end if   

         ! Using the actual nitrogen allocated to the leaf after
         ! uptake rather than fixing leaf nitrogen based on SLA and CN
         ! ratio
         if (lnc_opt .eqv. .true.) then                                                     
            ! nlevcan and nrad(p) look like the same variable ?? check this later
            sum_nscaler = 0.0_r8                                                    
            laican = 0.0_r8                                                         
            total_lai = 0.0_r8                                                      

            do iv = 1, nrad(p)                                                      

               if (iv == 1) then                                                    
                  laican = 0.5_r8 * tlai_z(p,iv)                                    
                  total_lai = tlai_z(p,iv)                                          
               else                                                                 
                  laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))          
                  total_lai = total_lai + tlai_z(p,iv)                              
               end if                                                               

               ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
               ! profile. If sun/shade big leaf code, use canopy integrated factor.
               if (nlevcan == 1) then                                               
                  nscaler = 1.0_r8                                                  
               else if (nlevcan > 1) then                                           
                  nscaler = exp(-kn(p) * laican)                                    
               end if                                                               

               sum_nscaler = sum_nscaler + nscaler                                  

            end do                                                                  

            if (tlai(p) > 0.0_r8 .AND. sum_nscaler > 0.0_r8) then
               ! dividing by LAI to convert total leaf nitrogen
               ! from m2 ground to m2 leaf; dividing by sum_nscaler to
               ! convert total leaf N to leaf N at canopy top
               lnc(p) = leafn(p) / (tlai(p) * sum_nscaler)
            else                                                                    
               lnc(p) = 0.0_r8                                                      
            end if                                                                  

         end if                                                                     


         ! reduce_dayl_factor .eqv. .false.  
         if (reduce_dayl_factor .eqv. .true.) then                                          
            if (dayl_factor(p) > 0.25_r8) then
               ! dayl_factor(p) = 1.0_r8  
            end if                                                                  
         end if                                                                     


         ! Default
         if (vcmax_opt == 0) then                                                   
            ! vcmax25 at canopy top, as in CN but using lnc at top of the canopy
            vcmax25top = lnc(p) * flnr(patch%itype(p)) * fnr * act25 * dayl_factor(p)
            if (.not. use_cn) then
               vcmax25top = vcmax25top * fnitr(patch%itype(p))
            else
               if ( CNAllocate_Carbon_only() ) vcmax25top = vcmax25top * fnitr(patch%itype(p))
            end if
         else if (vcmax_opt == 3) then                                                                   
            vcmax25top = ( i_vcad(patch%itype(p)) + s_vcad(patch%itype(p)) * lnc(p) ) * dayl_factor(p)  
         else if (vcmax_opt == 4) then                                                                   
            nptreemax = 9  ! is this number correct? check later 
            if (patch%itype(p) >= nptreemax) then   ! if not tree 
               ! for shrubs and herbs 
               vcmax25top = lnc(p) * ( i_flnr(patch%itype(p)) + s_flnr(patch%itype(p)) * lnc(p) ) * fnr * act25 * &
                    dayl_factor(p)
            else
               ! if tree 
               vcmax25top = lnc(p) * ( i_flnr(patch%itype(p)) * exp(s_flnr(patch%itype(p)) * lnc(p)) ) * fnr * act25 * &
                    dayl_factor(p)
               ! for trees 
            end if     
         end if        


         ! Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
         ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of Experimental Botany 44:907-920.

         jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top
         tpu25top  = 0.167_r8 * vcmax25top
         kp25top   = 20000._r8 * vcmax25top

         ! Nitrogen scaling factor. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
         ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al (2010) Biogeosciences, 7, 1833-1859
         ! Remove daylength factor from vcmax25 so that kn is based on maximum vcmax25
         ! But not used as defined here if using sun/shade big leaf code. Instead,
         ! will use canopy integrated scaling factors from SurfaceAlbedo.

         if (dayl_factor(p)  < 1.0e-12_r8) then
            kn(p) =  0._r8
         else
            kn(p) = exp(0.00963_r8 * vcmax25top/dayl_factor(p) - 2.43_r8)
         end if

         if (use_cn) then
            if ( leafresp_method == leafresp_mtd_ryan1991 ) then
            ! Leaf maintenance respiration to match the base rate used in CN
            ! but with the new temperature functions for C3 and C4 plants.
            !
            ! Base rate for maintenance respiration is from:
            ! M. Ryan, 1991. Effects of climate change on plant respiration.
            ! Ecological Applications, 1(2), 157-167.
            ! Original expression is br = 0.0106 molC/(molN h)
            ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
            !
            ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
            !
            ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
            ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s
            !
            ! Then scale this value at the top of the canopy for canopy depth

               lmr25top = 2.525e-6_r8 * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
               lmr25top = lmr25top * lnc(p) / 12.e-06_r8
            
            else if ( leafresp_method == leafresp_mtd_atkin2015 ) then
               !using new form for respiration base rate from Atkin
               !communication. 
               if ( lnc(p) > 0.0_r8 ) then
                  lmr25top = params_inst%lmr_intercept_atkin(ivt(p)) + (lnc(p) * 0.2061_r8) - (0.0402_r8 * (t10(p)-tfrz))
               else
                  lmr25top = 0.0_r8
               end if
            end if

         else
            ! Leaf maintenance respiration in proportion to vcmax25top

            if (c3flag(p)) then
               lmr25top = vcmax25top * leaf_mr_vcm
            else
               lmr25top = vcmax25top * 0.025_r8
            end if
         end if

         ! Loop through canopy layers (above snow). Respiration needs to be
         ! calculated every timestep. Others are calculated only if daytime

         laican = 0._r8
         do iv = 1, nrad(p)

            ! Cumulative lai at middle of layer

            if (iv == 1) then
               laican = 0.5_r8 * tlai_z(p,iv)
            else
               laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))
            end if

            ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
            ! profile. If sun/shade big leaf code, use canopy integrated factor.

            if (nlevcan == 1) then
               nscaler = vcmaxcint(p)
            else if (nlevcan > 1) then
               nscaler = exp(-kn(p) * laican)
            end if

            ! Maintenance respiration

            lmr25 = lmr25top * nscaler

            if(use_luna.and.c3flag(p).and.crop(patch%itype(p))== 0) then
                if(.not.use_cn)then ! If CN is on, use leaf N to predict respiration (above). Otherwise, use Vcmax term from LUNA.  RF
                  lmr25 = leaf_mr_vcm * photosyns_inst%vcmx25_z_patch(p,iv)
                endif
            endif
          
            if (c3flag(p)) then
               lmr_z(p,iv) = lmr25 * ft(t_veg(p), lmrha) * fth(t_veg(p), lmrhd, lmrse, lmrc)
            else
               lmr_z(p,iv) = lmr25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               lmr_z(p,iv) = lmr_z(p,iv) / (1._r8 + exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
            end if

            if (par_z(p,iv) <= 0._r8) then           ! night time

               vcmax_z(p,iv) = 0._r8
               jmax_z(p,iv) = 0._r8
               tpu_z(p,iv) = 0._r8
               kp_z(p,iv) = 0._r8

               if ( use_c13 ) then
                  alphapsn(p) = 1._r8
               end if

            else                                     ! day time

               if(use_luna.and.c3flag(p).and.crop(patch%itype(p))== 0)then
                  vcmax25 = photosyns_inst%vcmx25_z_patch(p,iv)
                  jmax25 = photosyns_inst%jmx25_z_patch(p,iv)
                  tpu25 = 0.167_r8 * vcmax25 
                  !Implement scaling of Vcmax25 from sunlit average to shaded canopy average value. RF & GBB. 1 July 2016
                  if(phase == 'sha'.and.surfalb_inst%vcmaxcintsun_patch(p).gt.0._r8.and.nlevcan==1) then
                     vcmax25 = vcmax25 * surfalb_inst%vcmaxcintsha_patch(p)/surfalb_inst%vcmaxcintsun_patch(p)
                     jmax25  = jmax25  * surfalb_inst%vcmaxcintsha_patch(p)/surfalb_inst%vcmaxcintsun_patch(p) 
                     tpu25   = tpu25   * surfalb_inst%vcmaxcintsha_patch(p)/surfalb_inst%vcmaxcintsun_patch(p) 
                  end if
                         
               else
                  vcmax25 = vcmax25top * nscaler
                  jmax25 = jmax25top * nscaler
                  tpu25 = tpu25top * nscaler                  
               endif
               kp25 = kp25top * nscaler

               ! Adjust for temperature

               vcmaxse = 668.39_r8 - 1.07_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               jmaxse  = 659.70_r8 - 0.75_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               tpuse = vcmaxse
               vcmaxc = fth25 (vcmaxhd, vcmaxse)
               jmaxc  = fth25 (jmaxhd, jmaxse)
               tpuc   = fth25 (tpuhd, tpuse)
               vcmax_z(p,iv) = vcmax25 * ft(t_veg(p), vcmaxha) * fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
               jmax_z(p,iv) = jmax25 * ft(t_veg(p), jmaxha) * fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
               tpu_z(p,iv) = tpu25 * ft(t_veg(p), tpuha) * fth(t_veg(p), tpuhd, tpuse, tpuc)

               if (.not. c3flag(p)) then
                  vcmax_z(p,iv) = vcmax25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                  vcmax_z(p,iv) = vcmax_z(p,iv) / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
                  vcmax_z(p,iv) = vcmax_z(p,iv) / (1._r8 + exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
               end if

               kp_z(p,iv) = kp25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)

            end if

            ! Adjust for soil water

            vcmax_z(p,iv) = vcmax_z(p,iv) * btran(p)
            lmr_z(p,iv) = lmr_z(p,iv) * btran(p)
            
           ! Change to add in light inhibition of respiration. 0.67 from Lloyd et al. 2010, & Metcalfe et al. 2012 
           ! Also pers. comm from Peter Reich (Nov 2015). Might potentially be updated pending findings of Atkin et al. (in prep)
           ! review of light inhibition database. 
           if ( light_inhibit .and. par_z(p,1) > 0._r8) then ! are the lights on? 
              lmr_z(p,iv) = lmr_z(p,iv) * 0.67_r8 ! inhibit respiration accordingly. 
           end if

         end do       ! canopy layer loop
      end do          ! patch loop

      !==============================================================================!
      ! Leaf-level photosynthesis and stomatal conductance
      !==============================================================================!

      rsmax0 = 2.e4_r8

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)
         g = patch%gridcell(p)

         ! Leaf boundary layer conductance, umol/m**2/s

         cf = forc_pbot(c)/(rgas*1.e-3_r8*tgcm(p))*1.e06_r8
         gb = 1._r8/rb(p)
         gb_mol(p) = gb * cf

         ! Loop through canopy layers (above snow). Only do calculations if daytime

         do iv = 1, nrad(p)

            if (par_z(p,iv) <= 0._r8) then           ! night time

               ac(p,iv) = 0._r8
               aj(p,iv) = 0._r8
               ap(p,iv) = 0._r8
               ag(p,iv) = 0._r8
               an(p,iv) = ag(p,iv) - lmr_z(p,iv)
               psn_z(p,iv) = 0._r8
               psn_wc_z(p,iv) = 0._r8
               psn_wj_z(p,iv) = 0._r8
               psn_wp_z(p,iv) = 0._r8
               rs_z(p,iv) = min(rsmax0, 1._r8/bbb(p) * cf)
               ci_z(p,iv) = 0._r8
               rh_leaf(p) = 0._r8

            else                                     ! day time

               !now the constraint is no longer needed, Jinyun Tang
               ceair = min( eair(p),  esat_tv(p) )
               if (      stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
                  rh_can = ceair / esat_tv(p)
               else if ( stomatalcond_mtd == stomatalcond_mtd_medlyn2011 )then
                  ! Put some constraints on RH in the canopy when Medlyn stomatal conductance is being used
                  rh_can = max((esat_tv(p) - ceair), 50._r8) * 0.001_r8
               end if

               ! Electron transport rate for C3 plants. Convert par from W/m2 to
               ! umol photons/m**2/s using the factor 4.6

               qabs  = 0.5_r8 * (1._r8 - fnps) * par_z(p,iv) * 4.6_r8
               aquad = theta_psii
               bquad = -(qabs + jmax_z(p,iv))
               cquad = qabs * jmax_z(p,iv)
               call quadratic (aquad, bquad, cquad, r1, r2)
               je = min(r1,r2)

               ! Iterative loop for ci beginning with initial guess

               if (c3flag(p)) then
                  ci_z(p,iv) = 0.7_r8 * cair(p)
               else
                  ci_z(p,iv) = 0.4_r8 * cair(p)
               end if

               niter = 0

               ! Increment iteration counter. Stop if too many iterations

               niter = niter + 1

               ! Save old ci

               ciold = ci_z(p,iv)

               !find ci and stomatal conductance
               call hybrid(ciold, p, iv, c, gb_mol(p), je, cair(p), oair(p), &
                    lmr_z(p,iv), par_z(p,iv), rh_can, gs_mol(p,iv), niter, &
                    atm2lnd_inst, photosyns_inst)

               ! End of ci iteration.  Check for an < 0, in which case gs_mol = bbb

               if (an(p,iv) < 0._r8) gs_mol(p,iv) = bbb(p)

               ! Get local noon sunlit and shaded stomatal conductance
               local_secp1 = secs + nint((grc%londeg(g)/degpsec)/dtime)*dtime
               local_secp1 = mod(local_secp1,isecspday)

               ! Use time period 1 hour before and 1 hour after local noon inclusive (11AM-1PM)
               if (local_secp1 >= (isecspday/2 - 3600) .and. local_secp1 <= (isecspday/2 + 3600)) then
                  if (phase == 'sun') then
                     gs_mol_sun_ln(p,iv) = gs_mol(p,iv)
                  else if (phase == 'sha') then
                     gs_mol_sha_ln(p,iv) = gs_mol(p,iv)
                  end if
               else
                  if (phase == 'sun') then
                     gs_mol_sun_ln(p,iv) = spval
                  else if (phase == 'sha') then
                     gs_mol_sha_ln(p,iv) = spval
                  end if
               end if

               ! Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)

               cs = cair(p) - 1.4_r8/gb_mol(p) * an(p,iv) * forc_pbot(c)
               cs = max(cs,1.e-06_r8)
               ci_z(p,iv) = cair(p) - an(p,iv) * forc_pbot(c) * (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv))

               ! Trap for values of ci_z less than 1.e-06.  This is needed for
               ! Megan (which can crash with negative values)
               ci_z(p,iv) = max( ci_z(p,iv), 1.e-06_r8 )

               ! Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)

               gs = gs_mol(p,iv) / cf
               rs_z(p,iv) = min(1._r8/gs, rsmax0)
               rs_z(p,iv) = rs_z(p,iv) / o3coefg(p)

               ! Photosynthesis. Save rate-limiting photosynthesis

               psn_z(p,iv) = ag(p,iv)
               psn_z(p,iv) = psn_z(p,iv) * o3coefv(p)

               psn_wc_z(p,iv) = 0._r8
               psn_wj_z(p,iv) = 0._r8
               psn_wp_z(p,iv) = 0._r8

               if (ac(p,iv) <= aj(p,iv) .and. ac(p,iv) <= ap(p,iv)) then
                  psn_wc_z(p,iv) =  psn_z(p,iv)
               else if (aj(p,iv) < ac(p,iv) .and. aj(p,iv) <= ap(p,iv)) then
                  psn_wj_z(p,iv) =  psn_z(p,iv)
               else if (ap(p,iv) < ac(p,iv) .and. ap(p,iv) < aj(p,iv)) then
                  psn_wp_z(p,iv) =  psn_z(p,iv)
               end if

               ! Make sure iterative solution is correct

               if (gs_mol(p,iv) < 0._r8) then
                  write (iulog,*)'Negative stomatal conductance:'
                  write (iulog,*)'p,iv,gs_mol= ',p,iv,gs_mol(p,iv)
                  call endrun(decomp_index=p, clmlevel=namep, msg=errmsg(sourcefile, __LINE__))
               end if

               ! Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b

               hs = (gb_mol(p)*ceair + gs_mol(p,iv)*esat_tv(p)) / ((gb_mol(p)+gs_mol(p,iv))*esat_tv(p))
               rh_leaf(p) = hs
               gs_mol_err = mbb(p)*max(an(p,iv), 0._r8)*hs/cs*forc_pbot(c) + bbb(p)

               if (abs(gs_mol(p,iv)-gs_mol_err) > 1.e-01_r8) then
                  write (iulog,*) 'Ball-Berry error check - stomatal conductance error:'
                  write (iulog,*) gs_mol(p,iv), gs_mol_err
               end if

            end if    ! night or day if branch
         end do       ! canopy layer loop
      end do          ! patch loop

      !==============================================================================!
      ! Canopy photosynthesis and stomatal conductance
      !==============================================================================!

      ! Sum canopy layer fluxes and then derive effective leaf-level fluxes (per
      ! unit leaf area), which are used in other parts of the model. Here, laican
      ! sums to either laisun or laisha.

      do f = 1, fn
         p = filterp(f)

         psncan = 0._r8
         psncan_wc = 0._r8
         psncan_wj = 0._r8
         psncan_wp = 0._r8
         lmrcan = 0._r8
         gscan = 0._r8
         laican = 0._r8
         do iv = 1, nrad(p)
            psncan = psncan + psn_z(p,iv) * lai_z(p,iv)
            psncan_wc = psncan_wc + psn_wc_z(p,iv) * lai_z(p,iv)
            psncan_wj = psncan_wj + psn_wj_z(p,iv) * lai_z(p,iv)
            psncan_wp = psncan_wp + psn_wp_z(p,iv) * lai_z(p,iv)
            lmrcan = lmrcan + lmr_z(p,iv) * lai_z(p,iv)
            gscan = gscan + lai_z(p,iv) / (rb(p)+rs_z(p,iv))
            laican = laican + lai_z(p,iv)
         end do
         if (laican > 0._r8) then
            psn(p) = psncan / laican
            psn_wc(p) = psncan_wc / laican
            psn_wj(p) = psncan_wj / laican
            psn_wp(p) = psncan_wp / laican
            lmr(p) = lmrcan / laican
            rs(p) = laican / gscan - rb(p)
         else
            psn(p) =  0._r8
            psn_wc(p) =  0._r8
            psn_wj(p) =  0._r8
            psn_wp(p) =  0._r8
            lmr(p) = 0._r8
            rs(p) = 0._r8
         end if
      end do

    end associate

  end subroutine Photosynthesis

  !------------------------------------------------------------------------------
  subroutine PhotosynthesisTotal (fn, filterp, &
       atm2lnd_inst, canopystate_inst, photosyns_inst)
    !
    ! Determine total photosynthesis
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: fn                             ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                    ! patch filter
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    !
    ! !LOCAL VARIABLES:
    integer :: f,fp,p,l,g               ! indices

    real(r8) :: rc14_atm(nsectors_c14), rc13_atm
    integer :: sector_c14
    !-----------------------------------------------------------------------

    associate(                                             &
         forc_pco2   => atm2lnd_inst%forc_pco2_grc       , & ! Input:  [real(r8) (:) ]  partial pressure co2 (Pa)
         forc_pc13o2 => atm2lnd_inst%forc_pc13o2_grc     , & ! Input:  [real(r8) (:) ]  partial pressure c13o2 (Pa)
         forc_po2    => atm2lnd_inst%forc_po2_grc        , & ! Input:  [real(r8) (:) ]  partial pressure o2 (Pa)

         laisun      => canopystate_inst%laisun_patch    , & ! Input:  [real(r8) (:) ]  sunlit leaf area
         laisha      => canopystate_inst%laisha_patch    , & ! Input:  [real(r8) (:) ]  shaded leaf area

         psnsun      => photosyns_inst%psnsun_patch      , & ! Input:  [real(r8) (:) ]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)
         psnsha      => photosyns_inst%psnsha_patch      , & ! Input:  [real(r8) (:) ]  shaded leaf photosynthesis (umol CO2 /m**2/ s)
         rc13_canair => photosyns_inst%rc13_canair_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in canopy air
         rc13_psnsun => photosyns_inst%rc13_psnsun_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in sunlit canopy psn flux
         rc13_psnsha => photosyns_inst%rc13_psnsha_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in shaded canopy psn flux
         alphapsnsun => photosyns_inst%alphapsnsun_patch , & ! Output: [real(r8) (:) ]  fractionation factor in sunlit canopy psn flux
         alphapsnsha => photosyns_inst%alphapsnsha_patch , & ! Output: [real(r8) (:) ]  fractionation factor in shaded canopy psn flux
         psnsun_wc   => photosyns_inst%psnsun_wc_patch   , & ! Output: [real(r8) (:) ]  Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
         psnsun_wj   => photosyns_inst%psnsun_wj_patch   , & ! Output: [real(r8) (:) ]  RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
         psnsun_wp   => photosyns_inst%psnsun_wp_patch   , & ! Output: [real(r8) (:) ]  product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
         psnsha_wc   => photosyns_inst%psnsha_wc_patch   , & ! Output: [real(r8) (:) ]  Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
         psnsha_wj   => photosyns_inst%psnsha_wj_patch   , & ! Output: [real(r8) (:) ]  RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
         psnsha_wp   => photosyns_inst%psnsha_wp_patch   , & ! Output: [real(r8) (:) ]  product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
         c13_psnsun  => photosyns_inst%c13_psnsun_patch  , & ! Output: [real(r8) (:) ]  sunlit leaf photosynthesis (umol 13CO2 /m**2/ s)
         c13_psnsha  => photosyns_inst%c13_psnsha_patch  , & ! Output: [real(r8) (:) ]  shaded leaf photosynthesis (umol 13CO2 /m**2/ s)
         c14_psnsun  => photosyns_inst%c14_psnsun_patch  , & ! Output: [real(r8) (:) ]  sunlit leaf photosynthesis (umol 14CO2 /m**2/ s)
         c14_psnsha  => photosyns_inst%c14_psnsha_patch  , & ! Output: [real(r8) (:) ]  shaded leaf photosynthesis (umol 14CO2 /m**2/ s)
         fpsn        => photosyns_inst%fpsn_patch        , & ! Output: [real(r8) (:) ]  photosynthesis (umol CO2 /m**2 /s)
         fpsn_wc     => photosyns_inst%fpsn_wc_patch     , & ! Output: [real(r8) (:) ]  Rubisco-limited photosynthesis (umol CO2 /m**2 /s)
         fpsn_wj     => photosyns_inst%fpsn_wj_patch     , & ! Output: [real(r8) (:) ]  RuBP-limited photosynthesis (umol CO2 /m**2 /s)
         fpsn_wp     => photosyns_inst%fpsn_wp_patch       & ! Output: [real(r8) (:) ]  product-limited photosynthesis (umol CO2 /m**2 /s)
         )

      if ( use_c14 ) then
         if (use_c14_bombspike) then
            call C14BombSpike(rc14_atm)
         else
            rc14_atm(:) = c14ratio
         end if
      end if

      if ( use_c13 ) then
         if (use_c13_timeseries) then
            call C13TimeSeries(rc13_atm)
         end if
      end if

      do f = 1, fn
         p = filterp(f)
         g = patch%gridcell(p)

         if (.not. use_fates) then
            fpsn(p)    = psnsun(p)   *laisun(p) + psnsha(p)   *laisha(p)
            fpsn_wc(p) = psnsun_wc(p)*laisun(p) + psnsha_wc(p)*laisha(p)
            fpsn_wj(p) = psnsun_wj(p)*laisun(p) + psnsha_wj(p)*laisha(p)
            fpsn_wp(p) = psnsun_wp(p)*laisun(p) + psnsha_wp(p)*laisha(p)
         end if

         if (use_cn) then
            if ( use_c13 ) then
               if (use_c13_timeseries) then
                  rc13_canair(p) = rc13_atm
               else
                  rc13_canair(p) = forc_pc13o2(g)/(forc_pco2(g) - forc_pc13o2(g))
               endif
               rc13_psnsun(p) = rc13_canair(p)/alphapsnsun(p)
               rc13_psnsha(p) = rc13_canair(p)/alphapsnsha(p)
               c13_psnsun(p)  = psnsun(p) * (rc13_psnsun(p)/(1._r8+rc13_psnsun(p)))
               c13_psnsha(p)  = psnsha(p) * (rc13_psnsha(p)/(1._r8+rc13_psnsha(p)))

               ! use fixed c13 ratio with del13C of -25 to test the overall c13 structure
               ! c13_psnsun(p) = 0.01095627 * psnsun(p)
               ! c13_psnsha(p) = 0.01095627 * psnsha(p)
            endif
            if ( use_c14 ) then

               ! determine latitute sector for radiocarbon bomb spike inputs
               if ( grc%latdeg(g) .ge. 30._r8 ) then
                  sector_c14 = 1
               else if ( grc%latdeg(g) .ge. -30._r8 ) then            
                  sector_c14 = 2
               else
                  sector_c14 = 3
               endif

               c14_psnsun(p) = rc14_atm(sector_c14) * psnsun(p)
               c14_psnsha(p) = rc14_atm(sector_c14) * psnsha(p)
            endif
         end if

      end do

    end associate

  end subroutine PhotosynthesisTotal

  !------------------------------------------------------------------------------
  subroutine Fractionation(bounds, fn, filterp, downreg, &
       atm2lnd_inst, canopystate_inst, solarabs_inst, surfalb_inst, photosyns_inst, &
       phase)
    !
    ! !DESCRIPTION:
    ! C13 fractionation during photosynthesis is calculated here after the nitrogen
    ! limitation is taken into account in the CNAllocation module.
    ! 
    ! As of CLM5, nutrient downregulation occurs prior to photosynthesis via leafcn, so we may
    ! ignore the downregulation term in this and assume that the Ci/Ca used in the photosynthesis
    ! calculation is consistent with that in the isotope calculation
    !
    !!USES:
    use clm_varctl     , only : use_hydrstress
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: fn                   ! size of pft filter
    integer                , intent(in)    :: filterp(fn)          ! patch filter
    real(r8)               , intent(in)    :: downreg( bounds%begp: ) ! fractional reduction in GPP due to N limitation (dimensionless)
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(solarabs_type)    , intent(in)    :: solarabs_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(photosyns_type)   , intent(in)    :: photosyns_inst
    character(len=*)       , intent(in)    :: phase                ! 'sun' or 'sha'
    !
    ! !LOCAL VARIABLES:
    real(r8) , pointer :: par_z (:,:)   ! needed for backwards compatiblity
    real(r8) , pointer :: alphapsn (:)  ! needed for backwards compatiblity
    real(r8) , pointer :: gs_mol(:,:)   ! leaf stomatal conductance (umol H2O/m**2/s)
    real(r8) , pointer :: an(:,:)       ! net leaf photosynthesis (umol CO2/m**2/s)
    integer  :: f,p,c,g,iv              ! indices
    real(r8) :: co2(bounds%begp:bounds%endp)  ! atmospheric co2 partial pressure (pa)
    real(r8) :: ci
    !------------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(downreg) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    associate(                                                  &
         forc_pbot   => atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         forc_pco2   => atm2lnd_inst%forc_pco2_grc            , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)

         c3psn       => pftcon%c3psn                          , & ! Input:  photosynthetic pathway: 0. = c4, 1. = c3

         nrad        => surfalb_inst%nrad_patch               , & ! Input:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer

         gb_mol      => photosyns_inst%gb_mol_patch             & ! Input:  [real(r8) (:)   ]  leaf boundary layer conductance (umol H2O/m**2/s)
         )

      if (phase == 'sun') then
         par_z    =>    solarabs_inst%parsun_z_patch     ! Input :  [real(r8) (:,:)] par absorbed per unit lai for canopy layer (w/m**2)
         alphapsn =>    photosyns_inst%alphapsnsun_patch ! Output:  [real(r8) (:)]
         if (use_hydrstress) then
            gs_mol => photosyns_inst%gs_mol_sun_patch    ! Input:   [real(r8) (:,:) ] sunlit leaf stomatal conductance (umol H2O/m**2/s)
            an     => photosyns_inst%an_sun_patch        ! Input:  [real(r8) (:,:) ]  net sunlit leaf photosynthesis (umol CO2/m**2/s)
         else
            gs_mol => photosyns_inst%gs_mol_patch        ! Input:   [real(r8) (:,:) ] leaf stomatal conductance (umol H2O/m**2/s)
            an     => photosyns_inst%an_patch            ! Input:  [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)
         end if
      else if (phase == 'sha') then
         par_z    =>    solarabs_inst%parsha_z_patch     ! Input :  [real(r8) (:,:)] par absorbed per unit lai for canopy layer (w/m**2)
         alphapsn =>    photosyns_inst%alphapsnsha_patch ! Output:  [real(r8) (:)]
         if (use_hydrstress) then
            gs_mol => photosyns_inst%gs_mol_sha_patch    ! Input:   [real(r8) (:,:) ] shaded leaf stomatal conductance (umol H2O/m**2/s)
            an     => photosyns_inst%an_sha_patch        ! Input:  [real(r8) (:,:) ]  net shaded leaf photosynthesis (umol CO2/m**2/s)
         else
            gs_mol => photosyns_inst%gs_mol_patch        ! Input:   [real(r8) (:,:) ] leaf stomatal conductance (umol H2O/m**2/s)
            an     => photosyns_inst%an_patch            ! Input:  [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)
         end if
      end if

      do f = 1, fn
         p = filterp(f)
         c= patch%column(p)
         g= patch%gridcell(p)

         co2(p) = forc_pco2(g)
         do iv = 1,nrad(p)
            if (par_z(p,iv) <= 0._r8) then           ! night time
               alphapsn(p) = 1._r8
            else                                     ! day time
               ci = co2(p) - (an(p,iv) * &
                    forc_pbot(c) * &
                    (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv)))
               alphapsn(p) = 1._r8 + (((c3psn(patch%itype(p)) * &
                    (4.4_r8 + (22.6_r8*(ci/co2(p))))) + &
                    ((1._r8 - c3psn(patch%itype(p))) * 4.4_r8))/1000._r8)
            end if
         end do
      end do

    end associate

  end subroutine Fractionation

  !-------------------------------------------------------------------------------
  subroutine hybrid(x0, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z,&
       rh_can, gs_mol,iter, &
       atm2lnd_inst, photosyns_inst)
    !
    !! DESCRIPTION:
    ! use a hybrid solver to find the root of equation
    ! f(x) = x- h(x),
    !we want to find x, s.t. f(x) = 0.
    !the hybrid approach combines the strength of the newton secant approach (find the solution domain)
    !and the bisection approach implemented with the Brent's method to guarrantee convergence.

    !
    !! REVISION HISTORY:
    !Dec 14/2012: created by Jinyun Tang
    !
    !!USES:
    !
    !! ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: x0              !initial guess and final value of the solution
    real(r8), intent(in) :: lmr_z              ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in) :: par_z              ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8), intent(in) :: rh_can             ! canopy air relative humidity
    real(r8), intent(in) :: gb_mol             ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: je                 ! electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: cair               ! Atmospheric CO2 partial pressure (Pa)
    real(r8), intent(in) :: oair               ! Atmospheric O2 partial pressure (Pa)
    integer,  intent(in) :: p, iv, c           ! pft, c3/c4, and column index
    real(r8), intent(out) :: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)
    integer,  intent(out) :: iter              !number of iterations used, for record only
    type(atm2lnd_type)  , intent(in)    :: atm2lnd_inst
    type(photosyns_type), intent(inout) :: photosyns_inst
    !
    !! LOCAL VARIABLES
    real(r8) :: a, b
    real(r8) :: fa, fb
    real(r8) :: x1, f0, f1
    real(r8) :: x, dx
    real(r8), parameter :: eps = 1.e-2_r8      !relative accuracy
    real(r8), parameter :: eps1= 1.e-4_r8
    integer,  parameter :: itmax = 40          !maximum number of iterations
    real(r8) :: tol,minx,minf

    call ci_func(x0, f0, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
         atm2lnd_inst, photosyns_inst)

    if(f0 == 0._r8)return

    minx=x0
    minf=f0
    x1 = x0 * 0.99_r8

    call ci_func(x1,f1, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
         atm2lnd_inst, photosyns_inst)

    if(f1==0._r8)then
       x0 = x1
       return
    endif
    if(f1<minf)then
       minx=x1
       minf=f1
    endif

    !first use the secant approach, then use the brent approach as a backup
    iter = 0
    do
       iter = iter + 1
       dx = - f1 * (x1-x0)/(f1-f0)
       x = x1 + dx
       tol = abs(x) * eps
       if(abs(dx)<tol)then
          x0 = x
          exit
       endif
       x0 = x1
       f0 = f1
       x1 = x

       call ci_func(x1,f1, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
            atm2lnd_inst, photosyns_inst)

       if(f1<minf)then
          minx=x1
          minf=f1
       endif
       if(abs(f1)<=eps1)then
          x0 = x1
          exit
       endif

       !if a root zone is found, use the brent method for a robust backup strategy
       if(f1 * f0 < 0._r8)then

          call brent(x, x0,x1,f0,f1, tol, p, iv, c, gb_mol, je, cair, oair, &
               lmr_z, par_z, rh_can, gs_mol, &
               atm2lnd_inst, photosyns_inst)

          x0=x
          exit
       endif
       if(iter>itmax)then
          !in case of failing to converge within itmax iterations
          !stop at the minimum function
          !this happens because of some other issues besides the stomatal conductance calculation
          !and it happens usually in very dry places and more likely with c4 plants.

          call ci_func(minx,f1, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
               atm2lnd_inst, photosyns_inst)

          exit
       endif
    enddo

  end subroutine hybrid

  !------------------------------------------------------------------------------
  subroutine brent(x, x1,x2,f1, f2, tol, ip, iv, ic, gb_mol, je, cair, oair,&
       lmr_z, par_z, rh_can, gs_mol, &
       atm2lnd_inst, photosyns_inst)
    !
    !!DESCRIPTION:
    !Use Brent's method to find the root of a single variable function ci_func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol.

    !!REVISION HISTORY:
    !Dec 14/2012: Jinyun Tang, modified from numerical recipes in F90 by press et al. 1188-1189
    !
    !!ARGUMENTS:
    real(r8), intent(out) :: x                ! indepedent variable of the single value function ci_func(x)
    real(r8), intent(in) :: x1, x2, f1, f2    ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in) :: tol               ! the error tolerance
    real(r8), intent(in) :: lmr_z             ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in) :: par_z             ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8), intent(in) :: gb_mol            ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: je                ! electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: cair              ! Atmospheric CO2 partial pressure (Pa)
    real(r8), intent(in) :: oair              ! Atmospheric O2 partial pressure (Pa)
    real(r8), intent(in) :: rh_can            ! inside canopy relative humidity
    integer,  intent(in) :: ip, iv, ic        ! pft, c3/c4, and column index
    real(r8), intent(out) :: gs_mol           ! leaf stomatal conductance (umol H2O/m**2/s)
    type(atm2lnd_type)  , intent(in)    :: atm2lnd_inst
    type(photosyns_type), intent(inout) :: photosyns_inst
    !
    !!LOCAL VARIABLES:
    integer, parameter :: itmax=20            !maximum number of iterations
    real(r8), parameter :: eps=1.e-2_r8       !relative error tolerance
    integer :: iter
    real(r8)  :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    !------------------------------------------------------------------------------

    a=x1
    b=x2
    fa=f1
    fb=f2
    if((fa > 0._r8 .and. fb > 0._r8).or.(fa < 0._r8 .and. fb < 0._r8))then
       write(iulog,*) 'root must be bracketed for brent'
       call endrun(msg=errmsg(sourcefile, __LINE__))
    endif
    c=b
    fc=fb
    iter = 0
    do
       if(iter==itmax)exit
       iter=iter+1
       if((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8))then
          c=a   !Rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
       endif
       if( abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2._r8*eps*abs(b)+0.5_r8*tol  !Convergence check.
       xm=0.5_r8*(c-b)
       if(abs(xm) <= tol1 .or. fb == 0.)then
          x=b
          return
       endif
       if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa !Attempt inverse quadratic interpolation.
          if(a == c) then
             p=2._r8*xm*s
             q=1._r8-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2._r8*xm*q*(q-r)-(b-a)*(r-1._r8))
             q=(q-1._r8)*(r-1._r8)*(s-1._r8)
          endif
          if(p > 0._r8) q=-q !Check whether in bounds.
          p=abs(p)
          if(2._r8*p < min(3._r8*xm*q-abs(tol1*q),abs(e*q))) then
             e=d !Accept interpolation.
             d=p/q
          else
             d=xm  !Interpolation failed, use bisection.
             e=d
          endif
       else !Bounds decreasing too slowly, use bisection.
          d=xm
          e=d
       endif
       a=b !Move last best guess to a.
       fa=fb
       if(abs(d) > tol1) then !Evaluate new trial root.
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif

       call ci_func(b, fb, ip, iv, ic, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
         atm2lnd_inst, photosyns_inst)

       if(fb==0._r8)exit

    enddo

    if(iter==itmax)write(iulog,*) 'brent exceeding maximum iterations', b, fb
    x=b

    return
  end subroutine brent

  !-------------------------------------------------------------------------------
  function ft(tl, ha) result(ans)
    !
    !!DESCRIPTION:
    ! photosynthesis temperature response
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    !!USES
    use clm_varcon  , only : rgas, tfrz
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: ha  ! activation energy in photosynthesis temperature function (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )

    return
  end function ft

  !-------------------------------------------------------------------------------
  function fth(tl,hd,se,scaleFactor) result(ans)
    !
    !!DESCRIPTION:
    !photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    use clm_varcon  , only : rgas, tfrz
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: hd  ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8), intent(in) :: se  ! entropy term in photosynthesis temperature function (J/mol/K)
    real(r8), intent(in) :: scaleFactor  ! scaling factor for high temperature inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )

    return
  end function fth

  !-------------------------------------------------------------------------------
  function fth25(hd,se)result(ans)
    !
    !!DESCRIPTION:
    ! scaling factor for photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY:
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    !!USES
    use clm_varcon  , only : rgas, tfrz
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hd    ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8), intent(in) :: se    ! entropy term in photosynthesis temperature function (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    return
  end function fth25

  !------------------------------------------------------------------------------
  subroutine ci_func(ci, fval, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z,&
       rh_can, gs_mol, atm2lnd_inst, photosyns_inst)
    !
    !! DESCRIPTION:
    ! evaluate the function
    ! f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an
    !
    ! remark:  I am attempting to maintain the original code structure, also
    ! considering one may be interested to output relevant variables for the
    ! photosynthesis model, I have decided to add these relevant variables to
    ! the relevant data types.
    !
    !!ARGUMENTS:
    real(r8)             , intent(in)    :: ci       ! intracellular leaf CO2 (Pa)
    real(r8)             , intent(in)    :: lmr_z    ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8)             , intent(in)    :: par_z    ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8)             , intent(in)    :: gb_mol   ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)             , intent(in)    :: je       ! electron transport rate (umol electrons/m**2/s)
    real(r8)             , intent(in)    :: cair     ! Atmospheric CO2 partial pressure (Pa)
    real(r8)             , intent(in)    :: oair     ! Atmospheric O2 partial pressure (Pa)
    real(r8)             , intent(in)    :: rh_can   ! canopy air realtive humidity
    integer              , intent(in)    :: p, iv, c ! pft, vegetation type and column indexes
    real(r8)             , intent(out)   :: fval     ! return function of the value f(ci)
    real(r8)             , intent(out)   :: gs_mol   ! leaf stomatal conductance (umol H2O/m**2/s)
    type(atm2lnd_type)   , intent(in)    :: atm2lnd_inst
    type(photosyns_type) , intent(inout) :: photosyns_inst
    !
    !local variables
    real(r8) :: ai                  ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
    real(r8) :: cs                  ! CO2 partial pressure at leaf surface (Pa)

    real(r8) :: aquad, bquad, cquad  ! terms for quadratic equations
    real(r8) :: r1, r2               ! roots of quadratic equation
    real(r8) :: fnps                 ! fraction of light absorbed by non-photosynthetic pigments
    real(r8) :: theta_psii           ! empirical curvature parameter for electron transport rate
    real(r8) :: theta_ip             ! empirical curvature parameter for ap photosynthesis co-limitation
    !------------------------------------------------------------------------------

    associate(&
         forc_pbot  => atm2lnd_inst%forc_pbot_downscaled_col   , & ! Output: [real(r8) (:)   ]  atmospheric pressure (Pa)
         c3flag     => photosyns_inst%c3flag_patch             , & ! Output: [logical  (:)   ]  true if C3 and false if C4
         ac         => photosyns_inst%ac_patch                 , & ! Output: [real(r8) (:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
         aj         => photosyns_inst%aj_patch                 , & ! Output: [real(r8) (:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)
         ap         => photosyns_inst%ap_patch                 , & ! Output: [real(r8) (:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag         => photosyns_inst%ag_patch                 , & ! Output: [real(r8) (:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)
         an         => photosyns_inst%an_patch                 , & ! Output: [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)
         vcmax_z    => photosyns_inst%vcmax_z_patch            , & ! Input:  [real(r8) (:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         cp         => photosyns_inst%cp_patch                 , & ! Output: [real(r8) (:)   ]  CO2 compensation point (Pa)
         kc         => photosyns_inst%kc_patch                 , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for CO2 (Pa)
         ko         => photosyns_inst%ko_patch                 , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for O2 (Pa)
         qe         => photosyns_inst%qe_patch                 , & ! Output: [real(r8) (:)   ]  quantum efficiency, used only for C4 (mol CO2 / mol photons)
         tpu_z      => photosyns_inst%tpu_z_patch              , & ! Output: [real(r8) (:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)
         kp_z       => photosyns_inst%kp_z_patch               , & ! Output: [real(r8) (:,:) ]  initial slope of CO2 response curve (C4 plants)
         theta_cj   => photosyns_inst%theta_cj_patch           , & ! Output: [real(r8) (:)   ]  empirical curvature parameter for ac, aj photosynthesis co-limitation
         bbb        => photosyns_inst%bbb_patch                , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
         mbb        => photosyns_inst%mbb_patch                  & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship
         )

      ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      fnps = 0.15_r8
      theta_psii = 0.7_r8
      theta_ip = 0.95_r8

      if (c3flag(p)) then
         ! C3: Rubisco-limited photosynthesis
         ac(p,iv) = vcmax_z(p,iv) * max(ci-cp(p), 0._r8) / (ci+kc(p)*(1._r8+oair/ko(p)))

         ! C3: RuBP-limited photosynthesis
         aj(p,iv) = je * max(ci-cp(p), 0._r8) / (4._r8*ci+8._r8*cp(p))

         ! C3: Product-limited photosynthesis
         ap(p,iv) = 3._r8 * tpu_z(p,iv)

      else

         ! C4: Rubisco-limited photosynthesis
         ac(p,iv) = vcmax_z(p,iv)

         ! C4: RuBP-limited photosynthesis
         aj(p,iv) = qe(p) * par_z * 4.6_r8

         ! C4: PEP carboxylase-limited (CO2-limited)
         ap(p,iv) = kp_z(p,iv) * max(ci, 0._r8) / forc_pbot(c)

      end if

      ! Gross photosynthesis. First co-limit ac and aj. Then co-limit ap

      aquad = theta_cj(p)
      bquad = -(ac(p,iv) + aj(p,iv))
      cquad = ac(p,iv) * aj(p,iv)
      call quadratic (aquad, bquad, cquad, r1, r2)
      ai = min(r1,r2)

      aquad = theta_ip
      bquad = -(ai + ap(p,iv))
      cquad = ai * ap(p,iv)
      call quadratic (aquad, bquad, cquad, r1, r2)
      ag(p,iv) = max(0._r8,min(r1,r2))

      ! Net photosynthesis. Exit iteration if an < 0

      an(p,iv) = ag(p,iv) - lmr_z
      if (an(p,iv) < 0._r8) then
         fval = 0._r8
         return
      endif
      ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
      ! With an <= 0, then gs_mol = bbb

      cs = cair - 1.4_r8/gb_mol * an(p,iv) * forc_pbot(c)
      cs = max(cs,1.e-06_r8)
      aquad = cs
      bquad = cs*(gb_mol - bbb(p)) - mbb(p)*an(p,iv)*forc_pbot(c)
      cquad = -gb_mol*(cs*bbb(p) + mbb(p)*an(p,iv)*forc_pbot(c)*rh_can)
      call quadratic (aquad, bquad, cquad, r1, r2)
      gs_mol = max(r1,r2)

      ! Derive new estimate for ci

      fval =ci - cair + an(p,iv) * forc_pbot(c) * (1.4_r8*gs_mol+1.6_r8*gb_mol) / (gb_mol*gs_mol)

    end associate

  end subroutine ci_func

  !------------------------------------------------------------------------------
  subroutine PhotosynthesisHydraulicStress ( bounds, fn, filterp, &
       esat_tv, eair, oair, cair, rb, bsun, bsha, btran, dayl_factor, leafn, &
       qsatl, qaf, &
       atm2lnd_inst, temperature_inst, soilstate_inst, waterdiagnosticbulk_inst, &
       surfalb_inst, solarabs_inst, canopystate_inst, ozone_inst, &
       photosyns_inst, waterfluxbulk_inst, froot_carbon, croot_carbon)
    !
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance calculation as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy
    ! Here, sunlit and shaded photosynthesis and stomatal conductance are solved 
    ! simultaneously per Pierre Gentine/Daniel Kennedy plant hydraulic stress
    ! method
    !
    ! !USES:
    use clm_varcon        , only : rgas, tfrz, rpi, spval, degpsec, isecspday
    use GridcellType      , only : grc
    use clm_time_manager  , only : get_curr_date, get_step_size
    use clm_varctl        , only : cnallocate_carbon_only
    use clm_varctl        , only : lnc_opt, reduce_dayl_factor, vcmax_opt    
    use clm_varpar        , only : nlevsoi
    use pftconMod         , only : nbrdlf_dcd_tmp_shrub, npcropmin
    use ColumnType        , only : col
    use shr_infnan_mod    , only : shr_infnan_isnan

    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: fn                             ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                    ! patch filter
    real(r8)               , intent(in)    :: esat_tv( bounds%begp: )        ! saturation vapor pressure at t_veg (Pa) [pft]
    real(r8)               , intent(in)    :: eair( bounds%begp: )           ! vapor pressure of canopy air (Pa) [pft]
    real(r8)               , intent(in)    :: oair( bounds%begp: )           ! Atmospheric O2 partial pressure (Pa) [pft]
    real(r8)               , intent(in)    :: cair( bounds%begp: )           ! Atmospheric CO2 partial pressure (Pa) [pft]
    real(r8)               , intent(in)    :: rb( bounds%begp: )             ! boundary layer resistance (s/m) [pft]
    real(r8)               , intent(in)    :: dayl_factor( bounds%begp: )    ! scalar (0-1) for daylength
    real(r8)               , intent(in)    :: qsatl ( bounds%begp: )         ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)    :: qaf ( bounds%begp: )           ! humidity of canopy air [kg/kg]
    real(r8)               , intent(in)    :: leafn( bounds%begp: )          ! leaf N (gN/m2)
    real(r8)               , intent(out)   :: bsun( bounds%begp: )           ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out)   :: bsha( bounds%begp: )           ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out)   :: btran( bounds%begp: )          ! transpiration wetness factor (0 to 1) [pft]
    real(r8)               , intent(in)    :: froot_carbon( bounds%begp: )    ! fine root carbon (gC/m2) [pft]   
    real(r8)               , intent(in)    :: croot_carbon( bounds%begp: )    ! live coarse root carbon (gC/m2) [pft]   

    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(solarabs_type)    , intent(in)    :: solarabs_inst
    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)   , intent(inout) :: waterfluxbulk_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    class(ozone_base_type) , intent(in)    :: ozone_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    !
    ! !LOCAL VARIABLES:
    !
    ! Leaf photosynthesis parameters
    real(r8) :: jmax_z(bounds%begp:bounds%endp,2,nlevcan) ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8) :: bbbopt(bounds%begp:bounds%endp)           ! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
    real(r8) :: kn(bounds%begp:bounds%endp)               ! leaf nitrogen decay coefficient
    real(r8) :: vcmax25top     ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25top      ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25top       ! canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top       ! canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25top        ! canopy top: initial slope of CO2 response curve (C4 plants) at 25C

    real(r8) :: vcmax25_sun    ! sunlit leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: vcmax25_sha    ! shaded leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25_sun     ! sunlit leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: jmax25_sha     ! shaded leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25_sun      ! sunlit leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: tpu25_sha      ! shaded leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25_sun      ! sunlit leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25_sha      ! shaded leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25_sun       ! sunlit leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
    real(r8) :: kp25_sha       ! shaded leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
    real(r8) :: kc25           ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(r8) :: ko25           ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(r8) :: cp25           ! CO2 compensation point at 25C (Pa)

    real(r8) :: vcmaxha        ! activation energy for vcmax (J/mol)
    real(r8) :: jmaxha         ! activation energy for jmax (J/mol)
    real(r8) :: tpuha          ! activation energy for tpu (J/mol)
    real(r8) :: lmrha          ! activation energy for lmr (J/mol)
    real(r8) :: kcha           ! activation energy for kc (J/mol)
    real(r8) :: koha           ! activation energy for ko (J/mol)
    real(r8) :: cpha           ! activation energy for cp (J/mol)

    real(r8) :: vcmaxhd        ! deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd         ! deactivation energy for jmax (J/mol)
    real(r8) :: tpuhd          ! deactivation energy for tpu (J/mol)
    real(r8) :: lmrhd          ! deactivation energy for lmr (J/mol)

    real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K)
    real(r8) :: tpuse          ! entropy term for tpu (J/mol/K)
    real(r8) :: lmrse          ! entropy term for lmr (J/mol/K)

    real(r8) :: vcmaxc         ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc          ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: tpuc           ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: lmrc           ! scaling factor for high temperature inhibition (25 C = 1.0)

    real(r8) :: fnps           ! fraction of light absorbed by non-photosynthetic pigments
    real(r8) :: theta_psii     ! empirical curvature parameter for electron transport rate

    real(r8) :: theta_ip       ! empirical curvature parameter for ap photosynthesis co-limitation

    ! Other
    integer  :: f,p,c,iv          ! indices
    real(r8) :: cf                ! s m**2/umol -> s/m
    real(r8) :: rsmax0            ! maximum stomatal resistance [s/m]
    real(r8) :: gb                ! leaf boundary layer conductance (m/s)
    real(r8) :: cs_sun            ! CO2 partial pressure at sunlit leaf surface (Pa)
    real(r8) :: cs_sha            ! CO2 partial pressure at shaded leaf surface (Pa)
    real(r8) :: gs                ! leaf stomatal conductance (m/s)
    real(r8) :: hs                ! fractional humidity at leaf surface (dimensionless)
    real(r8) :: sco               ! relative specificity of rubisco
    real(r8) :: ft                ! photosynthesis temperature response (statement function)
    real(r8) :: fth               ! photosynthesis temperature inhibition (statement function)
    real(r8) :: fth25             ! ccaling factor for photosynthesis temperature inhibition (statement function)
    real(r8) :: tl                ! leaf temperature in photosynthesis temperature function (K)
    real(r8) :: ha                ! activation energy in photosynthesis temperature function (J/mol)
    real(r8) :: hd                ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8) :: se                ! entropy term in photosynthesis temperature function (J/mol/K)
    real(r8) :: scaleFactor       ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: ciold             ! previous value of Ci for convergence check
    real(r8) :: gs_mol_err        ! gs_mol for error check
    real(r8) :: je_sun            ! sunlit leaf electron transport rate (umol electrons/m**2/s)
    real(r8) :: je_sha            ! shaded leaf electron transport rate (umol electrons/m**2/s)
    real(r8) :: qabs              ! PAR absorbed by PS II (umol photons/m**2/s)
    real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
    real(r8) :: r1,r2             ! roots of quadratic equation
    real(r8) :: ceair             ! vapor pressure of air, constrained (Pa)
    real(r8) :: fnr               ! (gRubisco/gN in Rubisco)
    real(r8) :: act25             ! (umol/mgRubisco/min) Rubisco activity at 25 C
    integer  :: iter1             ! number of iterations used, for record only
    integer  :: iter2             ! number of iterations used, for record only 
    real(r8) :: nscaler           ! leaf nitrogen scaling coefficient
    real(r8) :: nscaler_sun       ! sunlit leaf nitrogen scaling coefficient
    real(r8) :: nscaler_sha       ! shaded leaf nitrogen scaling coefficient

    real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)

    real(r8) :: psn_wc_z_sun(bounds%begp:bounds%endp,nlevcan) ! Rubisco-limited contribution to sunlit psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wj_z_sun(bounds%begp:bounds%endp,nlevcan) ! RuBP-limited contribution to sunlit psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wp_z_sun(bounds%begp:bounds%endp,nlevcan) ! product-limited contribution to sunlit psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wc_z_sha(bounds%begp:bounds%endp,nlevcan) ! Rubisco-limited contribution to shaded psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wj_z_sha(bounds%begp:bounds%endp,nlevcan) ! RuBP-limited contribution to shaded psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wp_z_sha(bounds%begp:bounds%endp,nlevcan) ! product-limited contribution to shaded psn_z (umol CO2/m**2/s)
    real(r8) :: rh_leaf_sun(bounds%begp:bounds%endp)          ! fractional humidity at sunlit leaf surface (dimensionless)
    real(r8) :: rh_leaf_sha(bounds%begp:bounds%endp)          ! fractional humidity at shaded leaf surface (dimensionless)

    real(r8) :: psncan_sun            ! canopy sum of sunlit psn_z
    real(r8) :: psncan_wc_sun         ! canopy sum of sunlit psn_wc_z
    real(r8) :: psncan_wj_sun         ! canopy sum of sunlit psn_wj_z
    real(r8) :: psncan_wp_sun         ! canopy sum of sunlit psn_wp_z
    real(r8) :: lmrcan_sun            ! canopy sum of sunlit lmr_z
    real(r8) :: gscan_sun             ! canopy sum of sunlit leaf conductance
    real(r8) :: laican_sun            ! canopy sum of sunlit lai_z
    real(r8) :: psncan_sha            ! canopy sum of shaded psn_z
    real(r8) :: psncan_wc_sha         ! canopy sum of shaded psn_wc_z
    real(r8) :: psncan_wj_sha         ! canopy sum of shaded psn_wj_z
    real(r8) :: psncan_wp_sha         ! canopy sum of shaded psn_wp_z
    real(r8) :: lmrcan_sha            ! canopy sum of shaded lmr_z
    real(r8) :: gscan_sha             ! canopy sum of shaded leaf conductance
    real(r8) :: laican_sha            ! canopy sum of shaded lai_z
    real(r8) :: laican                ! canopy sum of lai_z
    real(r8) :: rh_can                ! canopy air relative humidity

    real(r8) , pointer :: lai_z_sun       (:,:) ! leaf area index for canopy layer, sunlit
    real(r8) , pointer :: par_z_sun       (:,:) ! par absorbed per unit lai for canopy layer, sunlit (w/m**2)
    real(r8) , pointer :: vcmaxcint_sun   (:)   ! leaf to canopy scaling coefficient, sunlit
    real(r8) , pointer :: alphapsn_sun    (:)   ! 13C fractionation factor for PSN, sunlit ()
    real(r8) , pointer :: psn_sun         (:)   ! foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_wc_sun      (:)   ! Rubisco-limited foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_wj_sun      (:)   ! RuBP-limited foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +] 
    real(r8) , pointer :: psn_wp_sun      (:)   ! product-limited foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_z_sun       (:,:) ! canopy layer: foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: lmr_sun         (:)   ! leaf maintenance respiration rate, sunlit (umol CO2/m**2/s)
    real(r8) , pointer :: lmr_z_sun       (:,:) ! canopy layer: leaf maintenance respiration rate, sunlit (umol CO2/m**2/s)
    real(r8) , pointer :: rs_sun          (:)   ! leaf stomatal resistance, sunlit (s/m)
    real(r8) , pointer :: rs_z_sun        (:,:) ! canopy layer: leaf stomatal resistance, sunlit (s/m)
    real(r8) , pointer :: ci_z_sun        (:,:) ! intracellular leaf CO2, sunlit (Pa)
    real(r8) , pointer :: o3coefv_sun     (:)   ! o3 coefficient used in photo calculation, sunlit
    real(r8) , pointer :: o3coefg_sun     (:)   ! o3 coefficient used in rs calculation, sunlit
    real(r8) , pointer :: lai_z_sha       (:,:) ! leaf area index for canopy layer, shaded
    real(r8) , pointer :: par_z_sha       (:,:) ! par absorbed per unit lai for canopy layer, shaded (w/m**2)
    real(r8) , pointer :: vcmaxcint_sha   (:)   ! leaf to canopy scaling coefficient, shaded
    real(r8) , pointer :: alphapsn_sha    (:)   ! 13C fractionation factor for PSN, shaded ()
    real(r8) , pointer :: psn_sha         (:)   ! foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_wc_sha      (:)   ! Rubisco-limited foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_wj_sha      (:)   ! RuBP-limited foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +] 
    real(r8) , pointer :: psn_wp_sha      (:)   ! product-limited foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_z_sha       (:,:) ! canopy layer: foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: lmr_sha         (:)   ! leaf maintenance respiration rate, shaded (umol CO2/m**2/s)
    real(r8) , pointer :: lmr_z_sha       (:,:) ! canopy layer: leaf maintenance respiration rate, shaded (umol CO2/m**2/s)
    real(r8) , pointer :: rs_sha          (:)   ! leaf stomatal resistance, shaded (s/m)
    real(r8) , pointer :: rs_z_sha        (:,:) ! canopy layer: leaf stomatal resistance, shaded (s/m)
    real(r8) , pointer :: ci_z_sha        (:,:) ! intracellular leaf CO2, shaded (Pa)
    real(r8) , pointer :: o3coefv_sha     (:)   ! o3 coefficient used in photo calculation, shaded
    real(r8) , pointer :: o3coefg_sha     (:)   ! o3 coefficient used in rs calculation, shaded
    real(r8) :: sum_nscaler
    real(r8) :: total_lai                
    integer  :: nptreemax                
    integer  :: local_secp1                     ! seconds into current date in local time
    real(r8) :: dtime                           ! land model time step (sec)
    integer  :: year,month,day,secs             ! calendar info for current time step
    integer  :: j,g                     ! index
    real(r8) :: rs_resis                ! combined soil-root resistance [s]
    real(r8) :: r_soil                  ! root spacing [m]
    real(r8) :: root_biomass_density    ! root biomass density [g/m3]
    real(r8) :: root_cross_sec_area     ! root cross sectional area [m2]
    real(r8) :: root_length_density     ! root length density [m/m3]
    real(r8) :: froot_average_length    ! average coarse root length [m]
    real(r8) :: croot_average_length    ! average coarse root length [m]
    real(r8) :: soil_conductance        ! soil to root hydraulic conductance [1/s]
    real(r8) :: root_conductance        ! root hydraulic conductance [1/s]
    real(r8) :: rai(nlevsoi)            ! root area index [m2/m2]
    real(r8) :: fs(nlevsoi)             ! root conductance scale factor (reduction in conductance due to decreasing (more negative) root water potential)
    real(r8) :: gsminsun                ! Minimum stomatal conductance sunlit
    real(r8) :: gsminsha                ! Minimum stomatal conductance shaded
    real(r8) :: gs_slope_sun            ! Slope stomatal conductance sunlit
    real(r8) :: gs_slope_sha            ! Slope stomatal conductance shaded
    real(r8), parameter :: croot_lateral_length = 0.25_r8   ! specified lateral coarse root length [m]
    real(r8), parameter :: c_to_b = 2.0_r8           !(g biomass /g C)
!Note that root density is for dry biomass not carbon. CLM provides root biomass as carbon. The conversion is 0.5 g C / g biomass
    integer, parameter :: noonsec   = isecspday / 2 ! seconds at local noon

    !------------------------------------------------------------------------------

    ! Temperature and soil water response functions

    ft(tl,ha) = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )
    fth(tl,hd,se,scaleFactor) = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )
    fth25(hd,se) = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    ! Enforce expected array sizes

    SHR_ASSERT_ALL((ubound(esat_tv)     == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(eair)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(oair)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(cair)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(rb)          == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(bsun)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(bsha)        == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(btran)       == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dayl_factor) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qsatl)       == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qaf)         == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    associate(                                                 &
         k_soil_root  => soilstate_inst%k_soil_root_patch    , & ! Input:  [real(r8) (:,:) ]  soil-root interface conductance (mm/s)
         hk_l         =>    soilstate_inst%hk_l_col          , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity (mm/s) 
         hksat        => soilstate_inst%hksat_col            , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         smp          => soilstate_inst%smp_l_col            , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]

         froot_leaf   => pftcon%froot_leaf                   , & ! fine root to leaf ratio 
         root_conductance_patch => soilstate_inst%root_conductance_patch , & ! Output:   [real(r8) (:,:)] root conductance
         soil_conductance_patch => soilstate_inst%soil_conductance_patch , & ! Output:   [real(r8) (:,:)] soil conductance
         rootfr       => soilstate_inst%rootfr_patch         , & ! Input:   [real(r8) (:,:)]
         dz           => col%dz                              , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         z            => col%z                               , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         root_radius  => pftcon%root_radius                  , & ! Input: 0.29e-03_r8 !(m) 
         root_density => pftcon%root_density                 , & ! Input: 0.31e06_r8 !(g biomass / m3 root) 
         tsai         => canopystate_inst%tsai_patch         , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         c3psn      => pftcon%c3psn                          , & ! Input:  photosynthetic pathway: 0. = c4, 1. = c3
         crop       => pftcon%crop                           , & ! Input:  crop or not (0 =not crop and 1 = crop)
         leafcn     => pftcon%leafcn                         , & ! Input:  leaf C:N (gC/gN)
         flnr       => pftcon%flnr                           , & ! Input:  fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
         fnitr      => pftcon%fnitr                          , & ! Input:  foliage nitrogen limitation factor (-)
         slatop     => pftcon%slatop                         , & ! Input:  specific leaf area at top of canopy, projected area basis [m^2/gC]
         dsladlai   => pftcon%dsladlai                       , & ! Input:  change in sla per unit lai
         i_vcad     => pftcon%i_vcad                         , & ! Input:  [real(r8) (:)   ]  
         s_vcad     => pftcon%s_vcad                         , & ! Input:  [real(r8) (:)   ]  
         i_flnr     => pftcon%i_flnr                         , & ! Input:  [real(r8) (:)   ]  
         s_flnr     => pftcon%s_flnr                         , & ! Input:  [real(r8) (:)   ]  
         mbbopt     => pftcon%mbbopt                         , & 
         medlynintercept=> pftcon%medlynintercept            , & ! Input:  [real(r8) (:)   ]  Intercept for Medlyn stomatal conductance model method
         medlynslope=> pftcon%medlynslope                    , & ! Input:  [real(r8) (:)   ]  Slope for Medlyn stomatal conductance model method
         forc_pbot  => atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         ivt        => patch%itype                           , & ! Input:  [integer  (:)   ]  patch vegetation type

         t_veg      => temperature_inst%t_veg_patch          , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)
         t10        => temperature_inst%t_a10_patch          , & ! Input:  [real(r8) (:)   ]  10-day running mean of the 2 m temperature (K)
         tgcm       => temperature_inst%thm_patch            , & ! Input:  [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)
         nrad       => surfalb_inst%nrad_patch               , & ! Input:  [integer  (:)   ]  pft number of canopy layers, above snow for radiative transfer
         tlai_z     => surfalb_inst%tlai_z_patch             , & ! Input:  [real(r8) (:,:) ]  pft total leaf area index for canopy layer
         tlai       => canopystate_inst%tlai_patch           , & ! Input:  [real(r8)(:)    ]  one-sided leaf area index, no burying by snow  
         c3flag     => photosyns_inst%c3flag_patch           , & ! Output: [logical  (:)   ]  true if C3 and false if C4
         ac         => photosyns_inst%ac_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
         aj         => photosyns_inst%aj_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)
         ap         => photosyns_inst%ap_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag         => photosyns_inst%ag_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)
         vcmax_z    => photosyns_inst%vcmax_z_phs_patch      , & ! Output: [real(r8) (:,:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         luvcmax25top => photosyns_inst%luvcmax25top_patch   , & !  Output: [real(r8) (:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         lujmax25top  => photosyns_inst%lujmax25top_patch    , & ! Output: [real(r8) (:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         lutpu25top   => photosyns_inst%lutpu25top_patch     , & ! Output: [real(r8) (:) ]  maximum rate of carboxylation (umol co2/m**2/s)
!!!
         tpu_z      => photosyns_inst%tpu_z_phs_patch        , & ! Output: [real(r8) (:,:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)
         kp_z       => photosyns_inst%kp_z_phs_patch         , & ! Output: [real(r8) (:,:,:) ]  initial slope of CO2 response curve (C4 plants)
         gb_mol     => photosyns_inst%gb_mol_patch           , & ! Output: [real(r8) (:)   ]  leaf boundary layer conductance (umol H2O/m**2/s)
         cp         => photosyns_inst%cp_patch               , & ! Output: [real(r8) (:)   ]  CO2 compensation point (Pa)
         kc         => photosyns_inst%kc_patch               , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for CO2 (Pa)
         ko         => photosyns_inst%ko_patch               , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for O2 (Pa)
         qe         => photosyns_inst%qe_patch               , & ! Output: [real(r8) (:)   ]  quantum efficiency, used only for C4 (mol CO2 / mol photons)
         theta_cj   => photosyns_inst%theta_cj_patch         , & ! Output: [real(r8) (:)   ]  empirical curvature parameter for ac, aj photosynthesis co-limitation
         bbb        => photosyns_inst%bbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
         mbb        => photosyns_inst%mbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship
         rh_leaf    => photosyns_inst%rh_leaf_patch          , & ! Output: [real(r8) (:)   ]  fractional humidity at leaf surface (dimensionless)
         lnc        => photosyns_inst%lnca_patch             , & ! Output: [real(r8) (:)   ]  top leaf layer leaf N concentration (gN leaf/m^2)
         light_inhibit=> photosyns_inst%light_inhibit        , & ! Input:  [logical        ]  flag if light should inhibit respiration
         leafresp_method=> photosyns_inst%leafresp_method    , & ! Input:  [integer        ]  method type to use for leaf-maint.-respiration at 25C canopy top
         stomatalcond_mtd=> photosyns_inst%stomatalcond_mtd  , & ! Input:  [integer        ]  method type to use for stomatal conductance
         modifyphoto_and_lmr_forcrop=> photosyns_inst%modifyphoto_and_lmr_forcrop, & ! Input:  [logical        ] modifyphoto_and_lmr_forcrop
         leaf_mr_vcm => canopystate_inst%leaf_mr_vcm         , & ! Input:  [real(r8)       ]  scalar constant of leaf respiration with Vcmax
         vegwp      => canopystate_inst%vegwp_patch          , & ! Input/Output: [real(r8) (:,:) ]  vegetation water matric potential (mm)
         an_sun     => photosyns_inst%an_sun_patch           , & ! Output: [real(r8) (:,:) ]  net sunlit leaf photosynthesis (umol CO2/m**2/s)
         an_sha     => photosyns_inst%an_sha_patch           , & ! Output: [real(r8) (:,:) ]  net shaded leaf photosynthesis (umol CO2/m**2/s)
         gs_mol_sun => photosyns_inst%gs_mol_sun_patch       , & ! Output: [real(r8) (:,:) ]  sunlit leaf stomatal conductance (umol H2O/m**2/s)
         gs_mol_sha => photosyns_inst%gs_mol_sha_patch       , & ! Output: [real(r8) (:,:) ]  shaded leaf stomatal conductance (umol H2O/m**2/s)
         gs_mol_sun_ln => photosyns_inst%gs_mol_sun_ln_patch , & ! Output: [real(r8) (:,:) ]  sunlit leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon (umol H2O/m**2/s)
         gs_mol_sha_ln => photosyns_inst%gs_mol_sha_ln_patch   & ! Output: [real(r8) (:,:) ]  shaded leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon (umol H2O/m**2/s)
         )

      par_z_sun     =>    solarabs_inst%parsun_z_patch        ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)
      lai_z_sun     =>    canopystate_inst%laisun_z_patch     ! Input:  [real(r8) (:,:) ]  leaf area index for canopy layer, sunlit or shaded
      vcmaxcint_sun =>    surfalb_inst%vcmaxcintsun_patch     ! Input:  [real(r8) (:)   ]  leaf to canopy scaling coefficient
      alphapsn_sun  =>    photosyns_inst%alphapsnsun_patch    ! Input:  [real(r8) (:)   ]  13C fractionation factor for PSN ()
      o3coefv_sun   =>    ozone_inst%o3coefvsun_patch         ! Input:  [real(r8) (:)   ]  O3 coefficient used in photosynthesis calculation
      o3coefg_sun   =>    ozone_inst%o3coefgsun_patch         ! Input:  [real(r8) (:)   ]  O3 coefficient used in rs calculation
      ci_z_sun      =>    photosyns_inst%cisun_z_patch        ! Output: [real(r8) (:,:) ]  intracellular leaf CO2 (Pa)
      rs_sun        =>    photosyns_inst%rssun_patch          ! Output: [real(r8) (:)   ]  leaf stomatal resistance (s/m)
      rs_z_sun      =>    photosyns_inst%rssun_z_patch        ! Output: [real(r8) (:,:) ]  canopy layer: leaf stomatal resistance (s/m)
      lmr_sun       =>    photosyns_inst%lmrsun_patch         ! Output: [real(r8) (:)   ]  leaf maintenance respiration rate (umol CO2/m**2/s)
      lmr_z_sun     =>    photosyns_inst%lmrsun_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
      psn_sun       =>    photosyns_inst%psnsun_patch         ! Output: [real(r8) (:)   ]  foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_z_sun     =>    photosyns_inst%psnsun_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wc_sun    =>    photosyns_inst%psnsun_wc_patch      ! Output: [real(r8) (:)   ]  Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wj_sun    =>    photosyns_inst%psnsun_wj_patch      ! Output: [real(r8) (:)   ]  RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wp_sun    =>    photosyns_inst%psnsun_wp_patch      ! Output: [real(r8) (:)   ]  product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      par_z_sha     =>    solarabs_inst%parsha_z_patch        ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)
      lai_z_sha     =>    canopystate_inst%laisha_z_patch     ! Input:  [real(r8) (:,:) ]  leaf area index for canopy layer, sunlit or shaded
      vcmaxcint_sha =>    surfalb_inst%vcmaxcintsha_patch     ! Input:  [real(r8) (:)   ]  leaf to canopy scaling coefficient
      alphapsn_sha  =>    photosyns_inst%alphapsnsha_patch    ! Input:  [real(r8) (:)   ]  13C fractionation factor for PSN ()
      o3coefv_sha   =>    ozone_inst%o3coefvsha_patch         ! Input:  [real(r8) (:)   ]  O3 coefficient used in photosynthesis calculation
      o3coefg_sha   =>    ozone_inst%o3coefgsha_patch         ! Input:  [real(r8) (:)   ]  O3 coefficient used in rs calculation
      ci_z_sha      =>    photosyns_inst%cisha_z_patch        ! Output: [real(r8) (:,:) ]  intracellular leaf CO2 (Pa)
      rs_sha        =>    photosyns_inst%rssha_patch          ! Output: [real(r8) (:)   ]  leaf stomatal resistance (s/m)
      rs_z_sha      =>    photosyns_inst%rssha_z_patch        ! Output: [real(r8) (:,:) ]  canopy layer: leaf stomatal resistance (s/m)
      lmr_sha       =>    photosyns_inst%lmrsha_patch         ! Output: [real(r8) (:)   ]  leaf maintenance respiration rate (umol CO2/m**2/s)
      lmr_z_sha     =>    photosyns_inst%lmrsha_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
      psn_sha       =>    photosyns_inst%psnsha_patch         ! Output: [real(r8) (:)   ]  foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_z_sha     =>    photosyns_inst%psnsha_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wc_sha    =>    photosyns_inst%psnsha_wc_patch      ! Output: [real(r8) (:)   ]  Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wj_sha    =>    photosyns_inst%psnsha_wj_patch      ! Output: [real(r8) (:)   ]  RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wp_sha    =>    photosyns_inst%psnsha_wp_patch      ! Output: [real(r8) (:)   ]  product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      !==============================================================================!
      ! Photosynthesis and stomatal conductance parameters, from:
      ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      !==============================================================================!

      ! Determine seconds off current time step

      dtime = get_step_size()
      call get_curr_date (year, month, day, secs)

      ! vcmax25 parameters, from CN

      fnr = 7.16_r8
      act25 = 3.6_r8   !umol/mgRubisco/min
      ! Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
      act25 = act25 * 1000.0_r8 / 60.0_r8

      ! Activation energy, from:
      ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
      ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282

      kcha    = 79430._r8
      koha    = 36380._r8
      cpha    = 37830._r8
      vcmaxha = 72000._r8
      jmaxha  = 50000._r8
      tpuha   = 72000._r8
      lmrha   = 46390._r8

      ! High temperature deactivation, from:
      ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
      ! The factor "c" scales the deactivation to a value of 1.0 at 25C

      vcmaxhd = 200000._r8
      jmaxhd  = 200000._r8
      tpuhd   = 200000._r8
      lmrhd   = 150650._r8
      lmrse   = 490._r8
      lmrc    = fth25 (lmrhd, lmrse)

! calculate root-soil interface conductance 
      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)
         
         do j = 1,nlevsoi

! calculate conversion from conductivity to conductance
            root_biomass_density = c_to_b * froot_carbon(p) * rootfr(p,j) / dz(c,j)
! ensure minimum root biomass (using 1gC/m2)
            root_biomass_density = max(c_to_b*1._r8,root_biomass_density)

          ! Root length density: m root per m3 soil
            root_cross_sec_area = rpi*root_radius(ivt(p))**2
            root_length_density = root_biomass_density / (root_density(ivt(p)) * root_cross_sec_area)

            ! Root-area index (RAI)
            rai(j) = (tsai(p)+tlai(p)) * froot_leaf(ivt(p)) * rootfr(p,j)

! fix coarse root_average_length to specified length
            croot_average_length = croot_lateral_length

! calculate r_soil using Gardner/spa equation (Bonan, GMD, 2014)
            r_soil = sqrt(1./(rpi*root_length_density)) 

            ! length scale approach
            soil_conductance = min(hksat(c,j),hk_l(c,j))/(1.e3*r_soil)
            
! use vegetation plc function to adjust root conductance
               fs(j)=  plc(smp(c,j),p,c,root,veg)
            
! krmax is root conductance per area per length
            root_conductance = (fs(j)*rai(j)*params_inst%krmax(ivt(p)))/(croot_average_length + z(c,j))

            soil_conductance = max(soil_conductance, 1.e-16_r8)
            root_conductance = max(root_conductance, 1.e-16_r8)

            root_conductance_patch(p,j) = root_conductance
            soil_conductance_patch(p,j) = soil_conductance

! sum resistances in soil and root
            rs_resis = 1._r8/soil_conductance + 1._r8/root_conductance

! conductance is inverse resistance
! explicitly set conductance to zero for top soil layer
            if(rai(j)*rootfr(p,j) > 0._r8 .and. j > 1) then
               k_soil_root(p,j) =  1._r8/rs_resis
            else
               k_soil_root(p,j) =  0.
            endif
            
         end do
      enddo

      ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593

      fnps = 0.15_r8
      theta_psii = 0.7_r8
      theta_ip = 0.95_r8

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)

         ! C3 or C4 photosynthesis logical variable

         if (nint(c3psn(patch%itype(p))) == 1) then
            c3flag(p) = .true.
         else if (nint(c3psn(patch%itype(p))) == 0) then
            c3flag(p) = .false.
         end if

         ! C3 and C4 dependent parameters

         if (c3flag(p)) then
            qe(p) = 0._r8
            theta_cj(p) = 0.98_r8
            bbbopt(p) = 10000._r8
         else
            qe(p) = 0.05_r8
            theta_cj(p) = 0.80_r8
            bbbopt(p) = 40000._r8
         end if
 
         if ( stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
            ! Soil water stress applied to Ball-Berry parameters later in ci_func_PHS
            bbb(p) = bbbopt(p) 
            mbb(p) = mbbopt(patch%itype(p))
         end if
         ! kc, ko, cp, from: Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
         !
         !       kc25 = 404.9 umol/mol
         !       ko25 = 278.4 mmol/mol
         !       cp25 = 42.75 umol/mol
         !
         ! Derive sco from cp and O2 using present-day O2 (0.209 mol/mol) and re-calculate
         ! cp to account for variation in O2 using cp = 0.5 O2 / sco
         !

         kc25 = (404.9_r8 / 1.e06_r8) * forc_pbot(c)
         ko25 = (278.4_r8 / 1.e03_r8) * forc_pbot(c)
         sco  = 0.5_r8 * 0.209_r8 / (42.75_r8 / 1.e06_r8)
         cp25 = 0.5_r8 * oair(p) / sco

         kc(p) = kc25 * ft(t_veg(p), kcha)
         ko(p) = ko25 * ft(t_veg(p), koha)
         cp(p) = cp25 * ft(t_veg(p), cpha)

      end do

      ! Multi-layer parameters scaled by leaf nitrogen profile.
      ! Loop through each canopy layer to calculate nitrogen profile using
      ! cumulative lai at the midpoint of the layer

      do f = 1, fn
         p = filterp(f)

         if (lnc_opt .eqv. .false.) then     
            ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
            lnc(p) = 1._r8 / (slatop(patch%itype(p)) * leafcn(patch%itype(p)))
         end if   

         ! Using the actual nitrogen allocated to the leaf after
         ! uptake rather than fixing leaf nitrogen based on SLA and CN
         ! ratio
         if (lnc_opt .eqv. .true.) then                                                     
            ! nlevcan and nrad(p) look like the same variable ?? check this later
            sum_nscaler = 0.0_r8                                                    
            laican = 0.0_r8                                                         
            total_lai = 0.0_r8                                                      

            do iv = 1, nrad(p)                                                      

               if (iv == 1) then                                                    
                  laican = 0.5_r8 * tlai_z(p,iv)                                    
                  total_lai = tlai_z(p,iv)                                          
               else                                                                 
                  laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))          
                  total_lai = total_lai + tlai_z(p,iv)                              
               end if                                                               

               ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
               ! profile. If sun/shade big leaf code, use canopy integrated factor.
               if (nlevcan == 1) then                                               
                  nscaler = 1.0_r8                                                  
               else if (nlevcan > 1) then                                           
                  nscaler = exp(-kn(p) * laican)                                    
               end if                                                               

               sum_nscaler = sum_nscaler + nscaler                                  

            end do                                                                  

            if (tlai(p) > 0.0_r8 .AND. sum_nscaler > 0.0_r8) then
               ! dividing by LAI to convert total leaf nitrogen
               ! from m2 ground to m2 leaf; dividing by sum_nscaler to
               ! convert total leaf N to leaf N at canopy top
               lnc(p) = leafn(p) / (tlai(p) * sum_nscaler)
            else                                                                    
               lnc(p) = 0.0_r8                                                      
            end if                                                                  

         end if                                                                     
         lnc(p) = min(lnc(p),10._r8)

         ! reduce_dayl_factor .eqv. .false.  
         if (reduce_dayl_factor .eqv. .true.) then                                          
            if (dayl_factor(p) > 0.25_r8) then
               ! dayl_factor(p) = 1.0_r8  
            end if                                                                  
         end if                                                                     


         ! Default
         if (vcmax_opt == 0) then                                                   
            ! vcmax25 at canopy top, as in CN but using lnc at top of the canopy
            vcmax25top = lnc(p) * flnr(patch%itype(p)) * fnr * act25 * dayl_factor(p)
            if (.not. use_cn) then
               vcmax25top = vcmax25top * fnitr(patch%itype(p))
            else
               if ( CNAllocate_Carbon_only() ) vcmax25top = vcmax25top * fnitr(patch%itype(p))
            end if
         else if (vcmax_opt == 3) then
            vcmax25top = ( i_vcad(patch%itype(p)) + s_vcad(patch%itype(p)) * lnc(p) ) * dayl_factor(p)
         else if (vcmax_opt == 4) then
            nptreemax = 9  ! is this number correct? check later
            if (patch%itype(p) >= nptreemax) then   ! if not tree
               ! for shrubs and herbs
               vcmax25top = lnc(p) * ( i_flnr(patch%itype(p)) + s_flnr(patch%itype(p)) * lnc(p) ) * fnr * act25 * &
                    dayl_factor(p)
            else
               ! if tree
               vcmax25top = lnc(p) * ( i_flnr(patch%itype(p)) * exp(s_flnr(patch%itype(p)) * lnc(p)) ) * fnr * act25 * &
                    dayl_factor(p)
               ! for trees
            end if
         end if

         ! Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
         ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of Experimental Botany 44:907-920.

         jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top
         tpu25top  = 0.167_r8 * vcmax25top
         kp25top   = 20000._r8 * vcmax25top
         luvcmax25top(p) = vcmax25top
         lujmax25top(p) = jmax25top
         lutpu25top(p)=tpu25top

         ! Nitrogen scaling factor. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
         ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al (2010) Biogeosciences, 7, 1833-1859
         ! Remove daylength factor from vcmax25 so that kn is based on maximum vcmax25
         ! But not used as defined here if using sun/shade big leaf code. Instead,
         ! will use canopy integrated scaling factors from SurfaceAlbedo.

         if (dayl_factor(p) .eq. 0._r8) then
            kn(p) =  0._r8
         else
            kn(p) = exp(0.00963_r8 * vcmax25top/dayl_factor(p) - 2.43_r8)
         end if

         if (use_cn) then
            if ( leafresp_method == leafresp_mtd_ryan1991 ) then
            ! Leaf maintenance respiration to match the base rate used in CN
            ! but with the new temperature functions for C3 and C4 plants.
            !
            ! Base rate for maintenance respiration is from:
            ! M. Ryan, 1991. Effects of climate change on plant respiration.
            ! Ecological Applications, 1(2), 157-167.
            ! Original expression is br = 0.0106 molC/(molN h)
            ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
            !
            ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
            !
            ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
            ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s
            !
            ! Then scale this value at the top of the canopy for canopy depth

               lmr25top = 2.525e-6_r8 * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
               lmr25top = lmr25top * lnc(p) / 12.e-06_r8

            else if ( leafresp_method == leafresp_mtd_atkin2015 ) then
               !using new form for respiration base rate from Atkin
               !communication. 
               if ( lnc(p) > 0.0_r8 ) then
                  lmr25top = params_inst%lmr_intercept_atkin(ivt(p)) + (lnc(p) * 0.2061_r8) - (0.0402_r8 * (t10(p)-tfrz))
               else
                  lmr25top = 0.0_r8
               end if
            end if

         else
            ! Leaf maintenance respiration in proportion to vcmax25top

            if (c3flag(p)) then
               lmr25top = vcmax25top * leaf_mr_vcm
            else
               lmr25top = vcmax25top * 0.025_r8
            end if
         end if

         ! Loop through canopy layers (above snow). Respiration needs to be
         ! calculated every timestep. Others are calculated only if daytime

         laican = 0._r8
         do iv = 1, nrad(p)

            ! Cumulative lai at middle of layer

            if (iv == 1) then
               laican = 0.5_r8 * tlai_z(p,iv)
            else
               laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))
            end if

            ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
            ! profile. If sun/shade big leaf code, use canopy integrated factor.

            if (nlevcan == 1) then
               nscaler_sun = vcmaxcint_sun(p)
               nscaler_sha = vcmaxcint_sha(p)
            else if (nlevcan > 1) then
               nscaler_sun = exp(-kn(p) * laican)
               nscaler_sha = exp(-kn(p) * laican)
            end if

            ! Maintenance respiration

            lmr25_sun = lmr25top * nscaler_sun
            lmr25_sha = lmr25top * nscaler_sha

            if(use_luna.and.c3flag(p).and.crop(patch%itype(p))== 0)then
                if(.not.use_cn)then ! If CN is on, use leaf N to predict respiration (above). Otherwise, use Vcmax term from LUNA.  RF
                  lmr25_sun = leaf_mr_vcm * photosyns_inst%vcmx25_z_patch(p,iv)
                  lmr25_sha = leaf_mr_vcm * photosyns_inst%vcmx25_z_patch(p,iv)
                endif
            endif
          
            if (c3flag(p)) then
               lmr_z_sun(p,iv) = lmr25_sun * ft(t_veg(p), lmrha) * fth(t_veg(p), lmrhd, lmrse, lmrc)
               lmr_z_sha(p,iv) = lmr25_sha * ft(t_veg(p), lmrha) * fth(t_veg(p), lmrhd, lmrse, lmrc)
            else
               lmr_z_sun(p,iv) = lmr25_sun * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               lmr_z_sun(p,iv) = lmr_z_sun(p,iv) / (1._r8 + exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
               lmr_z_sha(p,iv) = lmr25_sha * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               lmr_z_sha(p,iv) = lmr_z_sha(p,iv) / (1._r8 + exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
            end if

            ! Reduce lmr w/ low lai
            lmr_z_sun(p,iv)  = lmr_z_sun(p,iv)*min((0.2_r8*exp(3.218_r8*tlai_z(p,iv))),1._r8)
            lmr_z_sha(p,iv)  = lmr_z_sha(p,iv)*min((0.2_r8*exp(3.218_r8*tlai_z(p,iv))),1._r8)

            if (par_z_sun(p,iv) <= 0._r8) then        ! night time

               vcmax_z(p,sun,iv) = 0._r8
               jmax_z(p,sun,iv) = 0._r8
               tpu_z(p,sun,iv) = 0._r8
               kp_z(p,sun,iv) = 0._r8

               vcmax_z(p,sha,iv) = 0._r8
               jmax_z(p,sha,iv) = 0._r8
               tpu_z(p,sha,iv) = 0._r8
               kp_z(p,sha,iv) = 0._r8

               if ( use_c13 ) then
                  alphapsn_sun(p) = 1._r8
                  alphapsn_sha(p) = 1._r8
               end if

            else                                     ! day time

               if(use_luna.and.c3flag(p).and.crop(patch%itype(p))== 0)then
                  vcmax25_sun = photosyns_inst%vcmx25_z_patch(p,iv)
                  vcmax25_sha = photosyns_inst%vcmx25_z_patch(p,iv)
                  jmax25_sun = photosyns_inst%jmx25_z_patch(p,iv)
                  jmax25_sha = photosyns_inst%jmx25_z_patch(p,iv)
                  tpu25_sun = 0.167_r8 * vcmax25_sun        
                  tpu25_sha = 0.167_r8 * vcmax25_sha        
                  if(surfalb_inst%vcmaxcintsun_patch(p).gt.0._r8.and.nlevcan==1) then
                    vcmax25_sha = vcmax25_sun * surfalb_inst%vcmaxcintsha_patch(p)/surfalb_inst%vcmaxcintsun_patch(p)
                    jmax25_sha  = jmax25_sun  * surfalb_inst%vcmaxcintsha_patch(p)/surfalb_inst%vcmaxcintsun_patch(p) 
                    tpu25_sha   = tpu25_sun   * surfalb_inst%vcmaxcintsha_patch(p)/surfalb_inst%vcmaxcintsun_patch(p) 
                  end if
               else
                  vcmax25_sun = vcmax25top * nscaler_sun
                  jmax25_sun = jmax25top * nscaler_sun
                  tpu25_sun = tpu25top * nscaler_sun        
                  vcmax25_sha = vcmax25top * nscaler_sha
                  jmax25_sha = jmax25top * nscaler_sha
                  tpu25_sha = tpu25top * nscaler_sha        
               endif
               kp25_sun = kp25top * nscaler_sun
               kp25_sha = kp25top * nscaler_sha

               ! Adjust for temperature

               vcmaxse = 668.39_r8 - 1.07_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               jmaxse  = 659.70_r8 - 0.75_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               tpuse = vcmaxse
               vcmaxc = fth25 (vcmaxhd, vcmaxse)
               jmaxc  = fth25 (jmaxhd, jmaxse)
               tpuc   = fth25 (tpuhd, tpuse)
               vcmax_z(p,sun,iv) = vcmax25_sun * ft(t_veg(p), vcmaxha) * fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
               jmax_z(p,sun,iv) = jmax25_sun * ft(t_veg(p), jmaxha) * fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
               tpu_z(p,sun,iv) = tpu25_sun * ft(t_veg(p), tpuha) * fth(t_veg(p), tpuhd, tpuse, tpuc)
               vcmax_z(p,sha,iv) = vcmax25_sha * ft(t_veg(p), vcmaxha) * fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
               jmax_z(p,sha,iv) = jmax25_sha * ft(t_veg(p), jmaxha) * fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
               tpu_z(p,sha,iv) = tpu25_sha * ft(t_veg(p), tpuha) * fth(t_veg(p), tpuhd, tpuse, tpuc)

               if (.not. c3flag(p)) then
                  vcmax_z(p,sun,iv) = vcmax25_sun * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                  vcmax_z(p,sun,iv) = vcmax_z(p,sun,iv) / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
                  vcmax_z(p,sun,iv) = vcmax_z(p,sun,iv) / (1._r8 + exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
                  vcmax_z(p,sha,iv) = vcmax25_sha * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                  vcmax_z(p,sha,iv) = vcmax_z(p,sha,iv) / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
                  vcmax_z(p,sha,iv) = vcmax_z(p,sha,iv) / (1._r8 + exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
               end if

               kp_z(p,sun,iv) = kp25_sun * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               kp_z(p,sha,iv) = kp25_sha * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)

            end if

           ! Change to add in light inhibition of respiration. 0.67 from Lloyd et al. 2010, & Metcalfe et al. 2012
           ! Also pers. comm from Peter Reich (Nov 2015). Might potentially be updated pending findings of Atkin et al. (in prep)
           ! review of light inhibition database.
           if ( light_inhibit .and. par_z_sun(p,1) > 0._r8) then ! are the lights on?
              lmr_z_sun(p,iv) = lmr_z_sun(p,iv) * 0.67_r8 ! inhibit respiration accordingly.
           end if
           if ( light_inhibit .and. par_z_sha(p,1) > 0._r8) then ! are the lights on?
              lmr_z_sha(p,iv) = lmr_z_sha(p,iv) * 0.67_r8 ! inhibit respiration accordingly.
           end if

         end do       ! canopy layer loop
      end do          ! patch loop

      !==============================================================================!
      ! Leaf-level photosynthesis and stomatal conductance
      !==============================================================================!

      rsmax0 = 2.e4_r8

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)
         g = patch%gridcell(p)

         ! Leaf boundary layer conductance, umol/m**2/s

         cf = forc_pbot(c)/(rgas*1.e-3_r8*tgcm(p))*1.e06_r8
         gb = 1._r8/rb(p)
         gb_mol(p) = gb * cf

         ! Loop through canopy layers (above snow). Only do calculations if daytime

         do iv = 1, nrad(p)

            if (par_z_sun(p,iv) <= 0._r8) then        ! night time

               !zqz temporary signal for night time
               vegwp(p,sun)=1._r8

               if (      stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
                  gsminsun = bbb(p)
                  gsminsha = bbb(p)
               else if ( stomatalcond_mtd == stomatalcond_mtd_medlyn2011 )then
                  gsminsun = medlynintercept(patch%itype(p))
                  gsminsha = medlynintercept(patch%itype(p))
               else
                  gsminsun = nan
                  gsminsha = nan
               end if
               call calcstress(p,c,vegwp(p,:),bsun(p),bsha(p),gb_mol(p),gsminsun, gsminsha, &
                    qsatl(p),qaf(p), atm2lnd_inst,canopystate_inst,waterdiagnosticbulk_inst, &
                    soilstate_inst,temperature_inst, waterfluxbulk_inst)

               ac(p,sun,iv) = 0._r8
               aj(p,sun,iv) = 0._r8
               ap(p,sun,iv) = 0._r8
               ag(p,sun,iv) = 0._r8
               if(crop(patch%itype(p))== 0 .or. .not. modifyphoto_and_lmr_forcrop) then
                  an_sun(p,iv) = ag(p,sun,iv) - bsun(p) * lmr_z_sun(p,iv)
               else
                  an_sun(p,iv) = ag(p,sun,iv) - lmr_z_sun(p,iv)
               endif
               psn_z_sun(p,iv) = 0._r8
               psn_wc_z_sun(p,iv) = 0._r8
               psn_wj_z_sun(p,iv) = 0._r8
               psn_wp_z_sun(p,iv) = 0._r8
               rs_z_sun(p,iv) = min(rsmax0, 1._r8/(max( bsun(p)*gsminsun, 1._r8 )) * cf)
               ci_z_sun(p,iv) = 0._r8
               rh_leaf_sun(p) = 0._r8

               ac(p,sha,iv) = 0._r8
               aj(p,sha,iv) = 0._r8
               ap(p,sha,iv) = 0._r8
               ag(p,sha,iv) = 0._r8
               if(crop(patch%itype(p))== 0 .or. .not. modifyphoto_and_lmr_forcrop) then
                  an_sha(p,iv) = ag(p,sha,iv) - bsha(p) * lmr_z_sha(p,iv)
               else
                  an_sha(p,iv) = ag(p,sha,iv) - lmr_z_sha(p,iv)
               endif
               psn_z_sha(p,iv) = 0._r8
               psn_wc_z_sha(p,iv) = 0._r8
               psn_wj_z_sha(p,iv) = 0._r8
               psn_wp_z_sha(p,iv) = 0._r8
               rs_z_sha(p,iv) = min(rsmax0, 1._r8/(max( bsha(p)*gsminsha, 1._r8 )) * cf)
               ci_z_sha(p,iv) = 0._r8
               rh_leaf_sha(p) = 0._r8

            else                                     ! day time

               !now the constraint is no longer needed, Jinyun Tang
               ceair = min( eair(p),  esat_tv(p) )
               if ( stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
                  rh_can = ceair / esat_tv(p)
               else if ( stomatalcond_mtd == stomatalcond_mtd_medlyn2011 )then
                  ! Put some constraints on RH in the canopy when Medlyn stomatal conductance is being used
                  rh_can = max((esat_tv(p) - ceair), 50._r8) * 0.001_r8
               end if

               ! Electron transport rate for C3 plants. Convert par from W/m2 to
               ! umol photons/m**2/s using the factor 4.6

               ! sun
               qabs  = 0.5_r8 * (1._r8 - fnps) * par_z_sun(p,iv) * 4.6_r8
               aquad = theta_psii
               bquad = -(qabs + jmax_z(p,sun,iv))
               cquad = qabs * jmax_z(p,sun,iv)
               call quadratic (aquad, bquad, cquad, r1, r2)
               je_sun = min(r1,r2)

               ! sha
               qabs  = 0.5_r8 * (1._r8 - fnps) * par_z_sha(p,iv) * 4.6_r8
               aquad = theta_psii
               bquad = -(qabs + jmax_z(p,sha,iv))
               cquad = qabs * jmax_z(p,sha,iv)
               call quadratic (aquad, bquad, cquad, r1, r2)
               je_sha = min(r1,r2)

               ! Iterative loop for ci beginning with initial guess

               if (c3flag(p)) then
                  ci_z_sun(p,iv) = 0.7_r8 * cair(p)
                  ci_z_sha(p,iv) = 0.7_r8 * cair(p)
               else
                  ci_z_sun(p,iv) = 0.4_r8 * cair(p)
                  ci_z_sha(p,iv) = 0.4_r8 * cair(p)
               end if

               !find ci and stomatal conductance
               call hybrid_PHS(ci_z_sun(p,iv), ci_z_sha(p,iv), p, iv, c, gb_mol(p), bsun(p),bsha(p), je_sun, &
                               je_sha, cair(p), oair(p), lmr_z_sun(p,iv), lmr_z_sha(p,iv), &
                               par_z_sun(p,iv), par_z_sha(p,iv), rh_can, gs_mol_sun(p,iv), gs_mol_sha(p,iv), &
                               qsatl(p), qaf(p), iter1, iter2, atm2lnd_inst, photosyns_inst, &
                               canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, temperature_inst, waterfluxbulk_inst)
               if ( stomatalcond_mtd == stomatalcond_mtd_medlyn2011 )then
                  gsminsun     = medlynintercept(patch%itype(p))
                  gsminsha     = medlynintercept(patch%itype(p))
                  gs_slope_sun = medlynslope(patch%itype(p))
                  gs_slope_sha = medlynslope(patch%itype(p))
               else if ( stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
                  gsminsun     = bbb(p)
                  gsminsha     = bbb(p)
                  gs_slope_sun = mbb(p)
                  gs_slope_sha = mbb(p)
               end if

               ! End of ci iteration.  Check for an < 0, in which case gs_mol = bbb

               if (an_sun(p,iv) < 0._r8) gs_mol_sun(p,iv) = max( bsun(p)*gsminsun, 1._r8 )
               if (an_sha(p,iv) < 0._r8) gs_mol_sha(p,iv) = max( bsha(p)*gsminsha, 1._r8 )
               ! Get local noon sunlit and shaded stomatal conductance
               local_secp1 = secs + nint((grc%londeg(g)/degpsec)/dtime)*dtime
               local_secp1 = mod(local_secp1,isecspday)
               ! Use time period 1 hour before and 1 hour after local noon inclusive (11AM-1PM)
               if (local_secp1 >= (isecspday/2 - 3600) .and. local_secp1 <= (isecspday/2 + 3600)) then
                  gs_mol_sun_ln(p,iv) = gs_mol_sun(p,iv)
                  gs_mol_sha_ln(p,iv) = gs_mol_sha(p,iv)
               else
                  gs_mol_sun_ln(p,iv) = spval
                  gs_mol_sha_ln(p,iv) = spval
               end if

               ! Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)

               cs_sun = cair(p) - 1.4_r8/gb_mol(p) * an_sun(p,iv) * forc_pbot(c)
               cs_sun = max(cs_sun,1.e-06_r8)
               ci_z_sun(p,iv) = cair(p) - an_sun(p,iv) * forc_pbot(c) * &
                                (1.4_r8*gs_mol_sun(p,iv)+1.6_r8*gb_mol(p)) / &
                                (gb_mol(p)*gs_mol_sun(p,iv))

               ! Trap for values of ci_z_sun less than 1.e-06.  This is needed for
               ! Megan (which can crash with negative values)
               ci_z_sun(p,iv) = max( ci_z_sun(p,iv), 1.e-06_r8 )

               cs_sha = cair(p) - 1.4_r8/gb_mol(p) * an_sha(p,iv) * forc_pbot(c)
               cs_sha = max(cs_sha,1.e-06_r8)
               ci_z_sha(p,iv) = cair(p) - an_sha(p,iv) * forc_pbot(c) * &
                                (1.4_r8*gs_mol_sha(p,iv)+1.6_r8*gb_mol(p)) / &
                                (gb_mol(p)*gs_mol_sha(p,iv))

               ! Trap for values of ci_z_sha less than 1.e-06.  This is needed for
               ! Megan (which can crash with negative values)
               ci_z_sha(p,iv) = max( ci_z_sha(p,iv), 1.e-06_r8 )

               ! Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)

               gs = gs_mol_sun(p,iv) / cf
               rs_z_sun(p,iv) = min(1._r8/gs, rsmax0)
               rs_z_sun(p,iv) = rs_z_sun(p,iv) / o3coefg_sun(p)
               gs = gs_mol_sha(p,iv) / cf
               rs_z_sha(p,iv) = min(1._r8/gs, rsmax0)
               rs_z_sha(p,iv) = rs_z_sha(p,iv) / o3coefg_sha(p)

               ! Photosynthesis. Save rate-limiting photosynthesis

               psn_z_sun(p,iv) = ag(p,sun,iv)
               psn_z_sun(p,iv) = psn_z_sun(p,iv) * o3coefv_sun(p)

               psn_wc_z_sun(p,iv) = 0._r8
               psn_wj_z_sun(p,iv) = 0._r8
               psn_wp_z_sun(p,iv) = 0._r8

               if (ac(p,sun,iv) <= aj(p,sun,iv) .and. ac(p,sun,iv) <= ap(p,sun,iv)) then
                  psn_wc_z_sun(p,iv) =  psn_z_sun(p,iv)
               else if (aj(p,sun,iv) < ac(p,sun,iv) .and. aj(p,sun,iv) <= ap(p,sun,iv)) then
                  psn_wj_z_sun(p,iv) =  psn_z_sun(p,iv)
               else if (ap(p,sun,iv) < ac(p,sun,iv) .and. ap(p,sun,iv) < aj(p,sun,iv)) then
                  psn_wp_z_sun(p,iv) =  psn_z_sun(p,iv)
               end if

               psn_z_sha(p,iv) = ag(p,sha,iv)
               psn_z_sha(p,iv) = psn_z_sha(p,iv) * o3coefv_sha(p)

               psn_wc_z_sha(p,iv) = 0._r8
               psn_wj_z_sha(p,iv) = 0._r8
               psn_wp_z_sha(p,iv) = 0._r8

               if (ac(p,sha,iv) <= aj(p,sha,iv) .and. ac(p,sha,iv) <= ap(p,sha,iv)) then
                  psn_wc_z_sha(p,iv) =  psn_z_sha(p,iv)
               else if (aj(p,sha,iv) < ac(p,sha,iv) .and. aj(p,sha,iv) <= ap(p,sha,iv)) then
                  psn_wj_z_sha(p,iv) =  psn_z_sha(p,iv)
               else if (ap(p,sha,iv) < ac(p,sha,iv) .and. ap(p,sha,iv) < aj(p,sha,iv)) then
                  psn_wp_z_sha(p,iv) =  psn_z_sha(p,iv)
               end if

               ! Make sure iterative solution is correct

               if (gs_mol_sun(p,iv) < 0._r8 .or. gs_mol_sha(p,iv) < 0._r8) then
                  write (iulog,*)'Negative stomatal conductance:'
                  write (iulog,*)'p,iv,gs_mol_sun,gs_mol_sha= ',p,iv,gs_mol_sun(p,iv),gs_mol_sha(p,iv)
                  call endrun(decomp_index=p, clmlevel=namep, msg=errmsg(sourcefile, __LINE__))
               end if

               ! Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b

               hs = (gb_mol(p)*ceair + gs_mol_sun(p,iv)*esat_tv(p)) / ((gb_mol(p)+gs_mol_sun(p,iv))*esat_tv(p))
               rh_leaf_sun(p) = hs
               gs_mol_err = gs_slope_sun*max(an_sun(p,iv), 0._r8)*hs/cs_sun*forc_pbot(c) + max( bsun(p)*gsminsun, 1._r8 )

               if (abs(gs_mol_sun(p,iv)-gs_mol_err) > 1.e-01_r8 .and.  (stomatalcond_mtd == stomatalcond_mtd_bb1987) ) then
                  write (iulog,*) 'Ball-Berry error check - sunlit stomatal conductance error:'
                  write (iulog,*) gs_mol_sun(p,iv), gs_mol_err
               end if

               hs = (gb_mol(p)*ceair + gs_mol_sha(p,iv)*esat_tv(p)) / ((gb_mol(p)+gs_mol_sha(p,iv))*esat_tv(p))
               rh_leaf_sha(p) = hs
               gs_mol_err = gs_slope_sha*max(an_sha(p,iv), 0._r8)*hs/cs_sha*forc_pbot(c) + max( bsha(p)*gsminsha, 1._r8)

               if (abs(gs_mol_sha(p,iv)-gs_mol_err) > 1.e-01_r8 .and.  (stomatalcond_mtd == stomatalcond_mtd_bb1987) ) then
                  write (iulog,*) 'Ball-Berry error check - shaded stomatal conductance error:'
                  write (iulog,*) gs_mol_sha(p,iv), gs_mol_err
               end if

            end if    ! night or day if branch
         end do       ! canopy layer loop
      end do          ! patch loop

      !==============================================================================!
      ! Canopy photosynthesis and stomatal conductance
      !==============================================================================!

      ! Sum canopy layer fluxes and then derive effective leaf-level fluxes (per
      ! unit leaf area), which are used in other parts of the model. Here, laican
      ! sums to either laisun or laisha.

      do f = 1, fn
         p = filterp(f)

         psncan_sun = 0._r8
         psncan_wc_sun = 0._r8
         psncan_wj_sun = 0._r8
         psncan_wp_sun = 0._r8
         lmrcan_sun = 0._r8
         gscan_sun = 0._r8
         laican_sun = 0._r8
         do iv = 1, nrad(p)
            psncan_sun = psncan_sun + psn_z_sun(p,iv) * lai_z_sun(p,iv)
            psncan_wc_sun = psncan_wc_sun + psn_wc_z_sun(p,iv) * lai_z_sun(p,iv)
            psncan_wj_sun = psncan_wj_sun + psn_wj_z_sun(p,iv) * lai_z_sun(p,iv)
            psncan_wp_sun = psncan_wp_sun + psn_wp_z_sun(p,iv) * lai_z_sun(p,iv)
            if(crop(patch%itype(p))== 0 .and. modifyphoto_and_lmr_forcrop) then
               lmrcan_sun = lmrcan_sun + lmr_z_sun(p,iv) * lai_z_sun(p,iv) * bsun(p)
            else
               lmrcan_sun = lmrcan_sun + lmr_z_sun(p,iv) * lai_z_sun(p,iv)
            endif
            gscan_sun = gscan_sun + lai_z_sun(p,iv) / (rb(p)+rs_z_sun(p,iv))
            laican_sun = laican_sun + lai_z_sun(p,iv)
         end do
         if (laican_sun > 0._r8) then
            psn_sun(p) = psncan_sun / laican_sun
            psn_wc_sun(p) = psncan_wc_sun / laican_sun
            psn_wj_sun(p) = psncan_wj_sun / laican_sun
            psn_wp_sun(p) = psncan_wp_sun / laican_sun
            lmr_sun(p) = lmrcan_sun / laican_sun
            rs_sun(p) = laican_sun / gscan_sun - rb(p)
         else
            psn_sun(p) =  0._r8
            psn_wc_sun(p) =  0._r8
            psn_wj_sun(p) =  0._r8
            psn_wp_sun(p) =  0._r8
            lmr_sun(p) = 0._r8
            rs_sun(p) = 0._r8
         end if
         psncan_sha = 0._r8
         psncan_wc_sha = 0._r8
         psncan_wj_sha = 0._r8
         psncan_wp_sha = 0._r8
         lmrcan_sha = 0._r8
         gscan_sha = 0._r8
         laican_sha = 0._r8
         do iv = 1, nrad(p)
            psncan_sha = psncan_sha + psn_z_sha(p,iv) * lai_z_sha(p,iv)
            psncan_wc_sha = psncan_wc_sha + psn_wc_z_sha(p,iv) * lai_z_sha(p,iv)
            psncan_wj_sha = psncan_wj_sha + psn_wj_z_sha(p,iv) * lai_z_sha(p,iv)
            psncan_wp_sha = psncan_wp_sha + psn_wp_z_sha(p,iv) * lai_z_sha(p,iv)
            if(crop(patch%itype(p))== 0 .and. modifyphoto_and_lmr_forcrop) then
               lmrcan_sha = lmrcan_sha + lmr_z_sha(p,iv) * lai_z_sha(p,iv) * bsha(p)
            else
               lmrcan_sha = lmrcan_sha + lmr_z_sha(p,iv) * lai_z_sha(p,iv)
            endif
            gscan_sha = gscan_sha + lai_z_sha(p,iv) / (rb(p)+rs_z_sha(p,iv))
            laican_sha = laican_sha + lai_z_sha(p,iv)
         end do
         if (laican_sha > 0._r8) then
            psn_sha(p) = psncan_sha / laican_sha
            psn_wc_sha(p) = psncan_wc_sha / laican_sha
            psn_wj_sha(p) = psncan_wj_sha / laican_sha
            psn_wp_sha(p) = psncan_wp_sha / laican_sha
            lmr_sha(p) = lmrcan_sha / laican_sha
            rs_sha(p) = laican_sha / gscan_sha - rb(p)
         else
            psn_sha(p) =  0._r8
            psn_wc_sha(p) =  0._r8
            psn_wj_sha(p) =  0._r8
            psn_wp_sha(p) =  0._r8
            lmr_sha(p) = 0._r8
            rs_sha(p) = 0._r8
         end if
         
         if ( laican_sha+laican_sun > 0._r8 ) then
            btran(p) = bsun(p) * (laican_sun / (laican_sun + laican_sha)) + &
                       bsha(p) * (laican_sha / (laican_sun + laican_sha))         
         else
            ! In this case, bsun and bsha should have the same value and btran 
            ! can be set to either bsun or bsha.
            btran(p) = bsun(p)
         end if

      end do

    end associate

  end subroutine PhotosynthesisHydraulicStress
  !------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  subroutine hybrid_PHS(x0sun, x0sha, p, iv, c, gb_mol, bsun, bsha, jesun, jesha, &
       cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
       gs_mol_sun, gs_mol_sha, qsatl, qaf, iter1, iter2, atm2lnd_inst, photosyns_inst, &
       canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, temperature_inst, waterfluxbulk_inst)
    !
    !! DESCRIPTION:
    !use a hybrid solver to find the root of the ci_func equation for sunlit and shaded leaves
    ! f(x) = x- h(x)                                                                                                                                               
    !we want to find x, s.t. f(x) = 0.
    !outside loop iterates for bsun/bsha, which are functions of stomatal conductance
    !the hybrid approach combines the strength of the newton secant approach (find the solution domain)
    !and the bisection approach implemented with the Brent's method to guarantee convergence.
    !
    !! REVISION HISTORY:
    !
    !
    !!USES:
    !
    !! ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: x0sun,x0sha              ! initial guess and final value of the solution for cisun/cisha
    integer , intent(in)    :: p                        ! pft index
    integer , intent(in)    :: iv                       ! radiation canopy layer index
    integer , intent(in)    :: c                        ! column index
    real(r8), intent(in)    :: gb_mol                   ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(out)   :: bsun                     ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8), intent(out)   :: bsha                     ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8), intent(in)    :: jesun                    ! sunlit leaf electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in)    :: jesha                    ! shaded leaf electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in)    :: cair                     ! Atmospheric CO2 partial pressure (Pa)
    real(r8), intent(in)    :: oair                     ! Atmospheric O2 partial pressure (Pa)
    real(r8), intent(in)    :: lmr_z_sun                ! sunlit canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in)    :: lmr_z_sha                ! shaded canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in)    :: par_z_sun                ! par absorbed per unit lai for sunlit canopy layer (w/m**2)
    real(r8), intent(in)    :: par_z_sha                ! par absorbed per unit lai for shaded canopy layer (w/m**2)
    real(r8), intent(in)    :: rh_can                   ! canopy air relative humidity
    real(r8), intent(out)   :: gs_mol_sun               ! sunlit leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(out)   :: gs_mol_sha               ! shaded leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(in)    :: qsatl                    ! leaf specific humidity [kg/kg]
    real(r8), intent(in)    :: qaf                      ! humidity of canopy air [kg/kg]
    integer,  intent(out)   :: iter1                    ! number of iterations used to find appropriate bsun/bsha
    integer,  intent(out)   :: iter2                    ! number of iterations used to find cisun/cisha
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)   , intent(inout) :: waterfluxbulk_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    !
    !! LOCAL VARIABLES
    real(r8) :: x(nvegwcs) ! working copy of vegwp(p,:)
    real(r8) :: gs0sun   ! unstressed sunlit stomatal conductance
    real(r8) :: gs0sha   ! unstressed shaded stomatal conductance
    logical  :: havegs   ! signals direction of calculation gs->qflx or qflx->gs
    real(r8) :: soilflux ! total soil column transpiration [mm/s] 
    real(r8) :: x1sun    ! second guess for cisun
    real(r8) :: f0sun    ! error of cifunc(x0sun)
    real(r8) :: f1sun    ! error of cifunc(x1sun)
    real(r8) :: xsun     ! open variable for brent to return cisun solution
    real(r8) :: dxsun    ! delta cisun from iter_i to iter_i+1
    real(r8) :: x1sha    ! second guess for cisha
    real(r8) :: f0sha    ! error of cifunc(x0sha)
    real(r8) :: f1sha    ! error of cifunc(x1sha)
    real(r8) :: xsha     ! open variable for brent to return cisha solution
    real(r8) :: dxsha    ! delta cisha from iter_i to iter_i+1
    real(r8) :: b0sun    ! bsun from previous iter
    real(r8) :: b0sha    ! bsha from previous iter
    real(r8) :: dbsun    ! delta(bsun) from iter_i to iter_i+1
    real(r8) :: dbsha    ! delta(bsun) from iter_i to iter_i+1
    logical  :: bflag    ! signals to call calcstress to recalc bsun/bsha (or not)
    real(r8) :: tolsun   ! error tolerance for cisun solution [Pa]
    real(r8) :: tolsha   ! error tolerance for cisun solution [Pa]
    real(r8) :: minf     ! storage spot for best cisun/cisha solution
    real(r8) :: minxsun  ! cisun associated with minf
    real(r8) :: minxsha  ! cisha associated with minf
    real(r8), parameter :: toldb = 1.e-2_r8  ! tolerance for satisfactory bsun/bsha solution
    real(r8), parameter :: eps = 1.e-2_r8    ! relative accuracy
    real(r8), parameter :: eps1= 1.e-4_r8    ! absolute accuracy threshold for fsun/fsha
    integer , parameter :: itmax = 3         ! maximum number of iterations zqz (increase later)
    !------------------------------------------------------------------------------
    
    associate(                                                    &
         qflx_tran_veg => waterfluxbulk_inst%qflx_tran_veg_patch    , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         vegwp         => canopystate_inst%vegwp_patch            & ! Input/Output: [real(r8) (:,:) ]  vegetation water matric potential (mm)
    )

    
    x1sun = x0sun
    x1sha = x0sha
    bflag = .false.
    b0sun = -1._r8
    b0sha = -1._r8
    gs0sun = 0._r8   ! Initialize to zero as good form, not used on first itteration below because of bflag
    gs0sha = 0._r8   ! Initialize to zero as good form, not used on first itteration below because of bflag
    bsun  = 1._r8
    bsha  = 1._r8
    iter1 = 0
    
    do                       !outer loop updates bsun/bsha and makes two ci_func calls for interpolation
       x=vegwp(p,:)
       iter1=iter1+1
       iter2=0
       x0sun=max(0.1_r8,x1sun)  !need to make sure x0 .neq. x1
       x1sun=0.99_r8*x1sun
       x0sha=max(0.1_r8,x1sha)
       x1sha=0.99_r8*x1sha
       tolsun = abs(x1sun) * eps
       tolsha = abs(x1sha) * eps
       
       ! this ci_func_PHS call updates bsun/bsha (except on first iter)
       call ci_func_PHS(x,x0sun, x0sha, f0sun, f0sha, p, iv, c, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
            gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
            qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, &
            temperature_inst, waterfluxbulk_inst)
       
       ! update bsun/bsha convergence vars
       dbsun=b0sun-bsun
       dbsha=b0sha-bsha
       b0sun=bsun
       b0sha=bsha
       bflag=.false.
       
       ! this ci_func_PHS call creates second point for ci interpolation
       call ci_func_PHS(x,x1sun, x1sha, f1sun, f1sha, p, iv, c, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
            gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
            qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, &
            temperature_inst, waterfluxbulk_inst)
       
       do                !inner loop finds ci
          if ( (abs(f0sun) < eps1) .and. (abs(f0sha) < eps1) ) then
             x1sun=x0sun
             x1sha=x0sha
             exit
          endif
          if ( (abs(f1sun) < eps1) .and. (abs(f1sha) < eps1) ) then
             exit
          endif
          iter2=iter2+1
          
          if ( (f1sun - f0sun) == 0._r8) then
             !makes next x1sun the midpt between current x1 & x0
             dxsun = 0.5_r8*(x1sun+x0sun)-x1sun
          else
             dxsun=-f1sun*(x1sun-x0sun)/(f1sun-f0sun)
          end if
          if ( (f1sha - f0sha) == 0._r8) then
             dxsha = 0.5_r8*(x1sha+x0sha)-x1sha
          else
             dxsha=-f1sha*(x1sha-x0sha)/(f1sha-f0sha)
          end if
          x0sun=x1sun
          x1sun=x1sun+dxsun
          x0sha=x1sha
          x1sha=x1sha+dxsha
          
          call ci_func_PHS(x,x1sun, x1sha, f1sun, f1sha, p, iv, c, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
               gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
               qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, &
               temperature_inst, waterfluxbulk_inst)

          if ( (abs(dxsun) < tolsun ) .and. (abs(dxsha) <tolsha) ) then
             x0sun=x1sun
             x0sha=x1sha
             exit
          endif
          if (iter2 .eq. 1) then
             !initialize best ci vars
             minf=abs(f1sun+f1sha)
             minxsun=x1sun
             minxsha=x1sha
          else
             if (abs(f1sun+f1sha)<minf) then
                !update best ci vars
                minf=abs(f1sun+f1sha)
                minxsun=x1sun
                minxsha=x1sha
             endif
          endif
          
          if ( (abs(f1sun) < eps1) .and. (abs(f1sha) < eps1) ) then
             exit
          endif
          
          if ( (f1sun*f0sun < 0._r8) .and. (f1sha*f0sha < 0._r8) ) then
             
             call brent_PHS(xsun, x0sun, x1sun, f0sun, f1sun, xsha, x0sha, x1sha, f0sha, f1sha, &
                  tolsun, p, iv, c, gb_mol, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,&
                  rh_can, gs_mol_sun, gs_mol_sha, bsun, bsha, qsatl, qaf, atm2lnd_inst, photosyns_inst, &
                  canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, temperature_inst, waterfluxbulk_inst)
             x0sun=xsun
             x0sha=xsha
             exit
          endif
          
          if (iter2 > itmax) then
             x1sun=minxsun
             x1sha=minxsha
             call ci_func_PHS(x,x1sun, x1sha, f1sun, f1sha, p, iv, c, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
                  gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
                  qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, &
                  temperature_inst, waterfluxbulk_inst)
             exit
          endif
          
       enddo
       
       !update unstressed stomatal conductance
       if (bsun>0.01_r8) then
          gs0sun=gs_mol_sun/bsun
       endif
       if (bsha>0.01_r8) then
          gs0sha=gs_mol_sha/bsha
       endif
       
       bflag=.true.
       
       if ( (abs(dbsun) < toldb) .and. (abs(dbsha) < toldb) ) then
          exit
       endif
       
       if (iter1 > itmax) then
          exit
       endif
    
    enddo
    x0sun=x1sun
    x0sha=x1sha
    
    !set vegwp for the final gs_mol solution
    call getvegwp(p, c, x, gb_mol, gs_mol_sun, gs_mol_sha, qsatl, qaf, soilflux, &
         atm2lnd_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, temperature_inst)
    vegwp(p,:)=x
    if (soilflux<0._r8) soilflux = 0._r8
    qflx_tran_veg(p) = soilflux
    
    end associate
    
  end subroutine hybrid_PHS
  !--------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine brent_PHS(xsun, x1sun, x2sun, f1sun, f2sun, xsha, x1sha, x2sha, f1sha, f2sha, &
       tol, ip, iv, ic, gb_mol, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,&
       rh_can, gs_mol_sun, gs_mol_sha, bsun, bsha, qsatl, qaf, atm2lnd_inst, photosyns_inst, &
       canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, temperature_inst, waterfluxbulk_inst)
    !------------------------------------------------------------------------------
    implicit none
    !
    !!DESCRIPTION:
    !Use Brent's method to find the root of a single variable function ci_func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol. Performed for cisun and cisha.
    !
    !!REVISION HISTORY:
    !
    !!ARGUMENTS:
    real(r8), intent(out)   :: xsun                 ! independent variable of the single value function ci_func(x)
    real(r8), intent(in)    :: x1sun, x2sun         ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in)    :: f1sun, f2sun         ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(out)   :: xsha                 ! independent variable of the single value function ci_func(x)
    real(r8), intent(in)    :: x1sha, x2sha         ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in)    :: f1sha, f2sha         ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in)    :: tol                  ! the error tolerance
    integer , intent(in)    :: ip, iv, ic           ! pft, c3/c4, and column index
    real(r8), intent(in)    :: gb_mol               ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in)    :: jesun,jesha          ! electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in)    :: cair                 ! Atmospheric CO2 partial pressure (Pa)
    real(r8), intent(in)    :: oair                 ! Atmospheric O2 partial pressure (Pa)
    real(r8), intent(in)    :: lmr_z_sun, lmr_z_sha ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in)    :: par_z_sun, par_z_sha ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8), intent(in)    :: rh_can               ! inside canopy relative humidity
    real(r8), intent(out)   :: gs_mol_sun           ! sunlit leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(out)   :: gs_mol_sha           ! shaded leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(inout) :: bsun                 ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8), intent(inout) :: bsha                 ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8), intent(in)    :: qsatl                ! leaf specific humidity [kg/kg]
    real(r8), intent(in)    :: qaf                  ! humidity of canopy air [kg/kg]
    type(atm2lnd_type)     , intent(in)       :: atm2lnd_inst
    type(photosyns_type)   , intent(inout)    :: photosyns_inst
    type(canopystate_type) , intent(inout)    :: canopystate_inst
    type(waterdiagnosticbulk_type)  , intent(inout)    :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)   , intent(inout)    :: waterfluxbulk_inst
    type(soilstate_type)   , intent(inout)    :: soilstate_inst
    type(temperature_type) , intent(in)       :: temperature_inst
    !------------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    real(r8)                :: gs0sun               ! sunlit leaf stomatal conductance (umol H2O/m**2/s)
    real(r8)                :: gs0sha               ! shaded leaf stomatal conductance (umol H2O/m**2/s)
    integer                 :: phase                ! sun==1, sha==2
    integer , parameter     :: nphs = 2             ! number of phases for sun/shade
    integer , parameter     :: itmax = 20           ! maximum number of iterations
    real(r8), parameter     :: eps = 1.e-4_r8       ! relative error tolerance
    integer                 :: iter                 !
    real(r8)                :: a(nphs),b(nphs),c(nphs),d(nphs),e(nphs),fa(nphs),fb(nphs),fc(nphs)
    real(r8)                :: p(nphs),q(nphs),r(nphs),s(nphs),tol1(nphs),xm(nphs)
    real(r8)                :: x(nvegwcs)           !dummy variable passed to cifunc
    logical , parameter     :: bflag = .false.      !indicates the cifunc should not call calcstress
    !------------------------------------------------------------------------------
    
    a(:)=(/x1sun,x1sha/)
    b(:)=(/x2sun,x2sha/)
    fa(:)=(/f1sun,f1sha/)
    fb(:)=(/f2sun,f2sha/)
    
    do phase=1, nphs
       if ( (fa(phase) > 0._r8 .and. fb(phase) > 0._r8) .or. (fa(phase) < 0._r8 .and. fb(phase) < 0._r8) ) then
          write(iulog,*) 'root must be bracketed for brent'
          call endrun(msg=errmsg(sourcefile, __LINE__))
       endif
    enddo
    
    c=b
    fc=fb
    iter = 0
    do
       if( iter == itmax ) exit
       iter=iter+1
       
       do phase=1, nphs
          if( (fb(phase) > 0._r8 .and. fc(phase) > 0._r8) .or. (fb(phase) < 0._r8 .and. fc(phase) < 0._r8)) then
             c(phase)=a(phase)   !Rename a, b, c and adjust bounding interval d.
             fc(phase)=fa(phase)
             d(phase)=b(phase)-a(phase)
             e(phase)=d(phase)
          endif
          if( abs(fc(phase)) < abs(fb(phase)) ) then
             a(phase)=b(phase)
             b(phase)=c(phase)
             c(phase)=a(phase)
             fa(phase)=fb(phase)
             fb(phase)=fc(phase)
             fc(phase)=fa(phase)
          endif
       enddo
       tol1=2._r8*eps*abs(b)+0.5_r8*tol  !Convergence check.
       xm=0.5_r8*(c-b)
       
       if( abs(xm(sun)) <= tol1(sun) .or. fb(sun) == 0._r8 ) then
          if( abs(xm(sha)) <= tol1(sha) .or. fb(sha) == 0._r8 ) then
             xsun=b(sun)
             xsha=b(sha)
             return
          endif
       endif
       
       do phase=1, nphs
          if( abs(e(phase)) >= tol1(phase) .and. abs(fa(phase)) > abs(fb(phase)) ) then
             s(phase)=fb(phase)/fa(phase) !Attempt inverse quadratic interpolation.
             if(a(phase) == c(phase)) then
                p(phase)=2._r8*xm(phase)*s(phase)
                q(phase)=1._r8-s(phase)
             else
                q(phase)=fa(phase)/fc(phase)
                r(phase)=fb(phase)/fc(phase)
                p(phase)=s(phase)*(2._r8*xm(phase)*q(phase)*(q(phase)-r(phase))-(b(phase)-a(phase))*(r(phase)-1._r8))
                q(phase)=(q(phase)-1._r8)*(r(phase)-1._r8)*(s(phase)-1._r8)
             endif
             if( p(phase) > 0._r8 ) q(phase)=-q(phase) !Check whether in bounds.
             p(phase)=abs(p(phase))
             if( 2._r8*p(phase) < min(3._r8*xm(phase)*q(phase)-abs(tol1(phase)*q(phase)),abs(e(phase)*q(phase))) ) then
                e(phase)=d(phase) !Accept interpolation.
                d(phase)=p(phase)/q(phase)
             else
                d(phase)=xm(phase)  !Interpolation failed, use bisection.
                e(phase)=d(phase)
             endif
          else !Bounds decreasing too slowly, use bisection.
             d(phase)=xm(phase)
             e(phase)=d(phase)
          endif
          a(phase)=b(phase) !Move last best guess to a.
          fa(phase)=fb(phase)
          if( abs(d(phase)) > tol1(phase) ) then !Evaluate new trial root.
             b(phase)=b(phase)+d(phase)
          else
             b(phase)=b(phase)+sign(tol1(phase),xm(phase))
          endif
       enddo
       
       gs0sun = gs_mol_sun
       gs0sha = gs_mol_sha
       call ci_func_PHS(x,b(sun), b(sha), fb(sun), fb(sha), ip, iv, ic, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha, &
            gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
            qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, &
            temperature_inst, waterfluxbulk_inst)
       
       if( (fb(sun) == 0._r8) .and. (fb(sha) == 0._r8) ) exit
    enddo
    if( iter == itmax) write(iulog,*) 'brent exceeding maximum iterations', b, fb
    xsun=b(sun)
    xsha=b(sha)
    
    return
    
  end subroutine brent_PHS
  !--------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine ci_func_PHS(x,cisun, cisha, fvalsun, fvalsha, p, iv, c, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
       gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
       qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, &
       temperature_inst, waterfluxbulk_inst)
    !------------------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! evaluate the function
    ! f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an for sunlit and shaded leaves
    !
    ! !REVISION HISTORY:
    !
    !
    ! !USES:
    use clm_varpar        , only : nlevsoi
    implicit none
    !
    ! !ARGUMENTS:
    real(r8)               , intent(inout) :: x(nvegwcs)         ! working copy of vegwp(p,:) 
    real(r8)               , intent(in)    :: cisun,cisha        ! intracellular leaf CO2 (Pa)
    real(r8)               , intent(out)   :: fvalsun,fvalsha    ! return function of the value f(ci)
    integer                , intent(in)    :: p,c,iv             ! pft, column, and radiation indexes
    real(r8)               , intent(inout) :: bsun               ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(inout) :: bsha               ! shaded canopy transpiration wetness factor (0 to 1)
    logical                , intent(in)    :: bflag              ! signals to call calcstress to recalc bsun/bsha (or not)
    real(r8)               , intent(in)    :: gb_mol             ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)    :: gs0sun,gs0sha      ! local gs_mol copies
    real(r8)               , intent(inout) :: gs_mol_sun,gs_mol_sha !leaf stomatal conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)    :: jesun, jesha       ! electron transport rate (umol electrons/m**2/s)
    real(r8)               , intent(in)    :: cair               ! Atmospheric CO2 partial pressure (Pa)
    real(r8)               , intent(in)    :: oair               ! Atmospheric O2 partial pressure (Pa)
    real(r8)               , intent(in)    :: lmr_z_sun, lmr_z_sha ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8)               , intent(in)    :: par_z_sun, par_z_sha ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8)               , intent(in)    :: rh_can             ! canopy air relative humidity
    real(r8)               , intent(in)    :: qsatl              ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)    :: qaf                ! humidity of canopy air [kg/kg]
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)   , intent(in)    :: waterfluxbulk_inst
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(temperature_type) , intent(in)    :: temperature_inst

    ! !LOCAL VARIABLES:
    real(r8) :: ai                   ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
    real(r8) :: cs_sun,cs_sha        ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: aquad, bquad, cquad  ! terms for quadratic equations
    real(r8) :: r1, r2               ! roots of quadratic equation
    real(r8) :: fnps                 ! fraction of light absorbed by non-photosynthetic pigments
    real(r8) :: theta_psii           ! empirical curvature parameter for electron transport rate
    real(r8) :: theta_ip             ! empirical curvature parameter for ap photosynthesis co-limitation
    real(r8) :: term                 ! intermediate in Medlyn stomatal model
    !
    !------------------------------------------------------------------------------
    
    associate(                                                 &
         forc_pbot  => atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ]    atmospheric pressure (Pa)
         c3flag     => photosyns_inst%c3flag_patch           , & ! Input:  [logical  (:)   ]    true if C3 and false if C4
         medlynslope=> pftcon%medlynslope                    , & ! Input:  [real(r8) (:)   ]  Slope for Medlyn stomatal conductance model method
         medlynintercept=> pftcon%medlynintercept            , & ! Input:  [real(r8) (:)   ]  Intercept for Medlyn stomatal conductance model method
         stomatalcond_mtd=> photosyns_inst%stomatalcond_mtd  , & ! Input:  [integer        ]  method type to use for stomatal conductance.GC.fnlprmsn15_r22845
         ac         => photosyns_inst%ac_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
         aj         => photosyns_inst%aj_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)
         ap         => photosyns_inst%ap_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag         => photosyns_inst%ag_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)
         vcmax_z    => photosyns_inst%vcmax_z_phs_patch      , & ! Input:  [real(r8) (:,:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         cp         => photosyns_inst%cp_patch               , & ! Output: [real(r8) (:)   ]    CO2 compensation point (Pa)
         kc         => photosyns_inst%kc_patch               , & ! Output: [real(r8) (:)   ]    Michaelis-Menten constant for CO2 (Pa)
         ko         => photosyns_inst%ko_patch               , & ! Output: [real(r8) (:)   ]    Michaelis-Menten constant for O2 (Pa)
         qe         => photosyns_inst%qe_patch               , & ! Output: [real(r8) (:)   ]    quantum efficiency, used only for C4 (mol CO2 / mol photons)
         tpu_z      => photosyns_inst%tpu_z_phs_patch        , & ! Output: [real(r8) (:,:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)
         kp_z       => photosyns_inst%kp_z_phs_patch         , & ! Output: [real(r8) (:,:,:) ]  initial slope of CO2 response curve (C4 plants)
         theta_cj   => photosyns_inst%theta_cj_patch         , & ! Output: [real(r8) (:)   ]    empirical curvature parameter for ac, aj photosynthesis co-limitation
         bbb        => photosyns_inst%bbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
         mbb        => photosyns_inst%mbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship
         an_sun     => photosyns_inst%an_sun_patch           , & ! Output: [real(r8) (:,:) ]  net sunlit leaf photosynthesis (umol CO2/m**2/s)
         an_sha     => photosyns_inst%an_sha_patch             & ! Output: [real(r8) (:,:) ]  net shaded leaf photosynthesis (umol CO2/m**2/s)
         )
    
    !------------------------------------------------------------------------------
    ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
    fnps = 0.15_r8
    theta_psii = 0.7_r8
    theta_ip = 0.95_r8
    
    if (bflag) then   !zqz what if bsun==0 ... doesn't break... but follow up

       call calcstress(p,c,x,bsun,bsha,gb_mol,gs0sun,gs0sha,qsatl,qaf, &
            atm2lnd_inst,canopystate_inst,waterdiagnosticbulk_inst,soilstate_inst, &
            temperature_inst, waterfluxbulk_inst)
    endif
    
    if (c3flag(p)) then
       ! C3: Rubisco-limited photosynthesis
       ac(p,sun,iv) = bsun * vcmax_z(p,sun,iv) * max(cisun-cp(p), 0._r8) / (cisun+kc(p)*(1._r8+oair/ko(p)))
       ac(p,sha,iv) = bsha * vcmax_z(p,sha,iv) * max(cisha-cp(p), 0._r8) / (cisha+kc(p)*(1._r8+oair/ko(p)))
       
       ! C3: RuBP-limited photosynthesis
       aj(p,sun,iv) = jesun * max(cisun-cp(p), 0._r8) / (4._r8*cisun+8._r8*cp(p))
       aj(p,sha,iv) = jesha * max(cisha-cp(p), 0._r8) / (4._r8*cisha+8._r8*cp(p))
       
       ! C3: Product-limited photosynthesis
       ap(p,sun,iv) = 3._r8 * tpu_z(p,sun,iv)
       ap(p,sha,iv) = 3._r8 * tpu_z(p,sha,iv)
       
    else
       ! C4: Rubisco-limited photosynthesis
       ac(p,sun,iv) = bsun * vcmax_z(p,sun,iv)
       ac(p,sha,iv) = bsha * vcmax_z(p,sha,iv)
       
       ! C4: RuBP-limited photosynthesis
       aj(p,sun,iv) = qe(p) * par_z_sun * 4.6_r8
       aj(p,sha,iv) = qe(p) * par_z_sha * 4.6_r8
       
       ! C4: PEP carboxylase-limited (CO2-limited)
       ap(p,sun,iv) = kp_z(p,sun,iv) * max(cisun, 0._r8) / forc_pbot(c)
       ap(p,sha,iv) = kp_z(p,sha,iv) * max(cisha, 0._r8) / forc_pbot(c)
       
    end if
    
    ! Gross photosynthesis. First co-limit ac and aj. Then co-limit ap
    
    ! Sunlit
    aquad = theta_cj(p)
    bquad = -(ac(p,sun,iv) + aj(p,sun,iv))
    cquad = ac(p,sun,iv) * aj(p,sun,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ai = min(r1,r2)
    
    aquad = theta_ip
    bquad = -(ai + ap(p,sun,iv))
    cquad = ai * ap(p,sun,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ag(p,sun,iv) = max(0._r8,min(r1,r2))
    
    ! Shaded
    aquad = theta_cj(p)
    bquad = -(ac(p,sha,iv) + aj(p,sha,iv))
    cquad = ac(p,sha,iv) * aj(p,sha,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ai = min(r1,r2)
    
    aquad = theta_ip
    bquad = -(ai + ap(p,sha,iv))
    cquad = ai * ap(p,sha,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ag(p,sha,iv) = max(0._r8,min(r1,r2))
    
    ! Net photosynthesis. Exit iteration if an < 0
    an_sun(p,iv) = ag(p,sun,iv) - bsun * lmr_z_sun
    an_sha(p,iv) = ag(p,sha,iv) - bsha * lmr_z_sha
    
    if (an_sun(p,iv) < 0._r8) then
       if ( stomatalcond_mtd == stomatalcond_mtd_medlyn2011 )then
          gs_mol_sun = medlynintercept(patch%itype(p))
       else if ( stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
          gs_mol_sun = bbb(p)
       else
          gs_mol_sun = nan
       end if
       gs_mol_sun = max( bsun*gs_mol_sun, 1._r8)
       fvalsun = 0._r8  ! really tho? zqz
    endif
    if (an_sha(p,iv) < 0._r8) then
       if ( stomatalcond_mtd == stomatalcond_mtd_medlyn2011 )then
          gs_mol_sha = medlynintercept(patch%itype(p))
       else if ( stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
          gs_mol_sha = bbb(p)
       else
          gs_mol_sha = nan
       end if
       gs_mol_sha = max( bsha*gs_mol_sha, 1._r8)
       fvalsha = 0._r8
    endif
    if ((an_sun(p,iv) < 0._r8) .AND. (an_sha(p,iv) < 0._r8)) then
       return
    endif
    
    ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
    ! With an <= 0, then gs_mol = bbb
    
    ! Sunlit
    cs_sun = cair - 1.4_r8/gb_mol * an_sun(p,iv) * forc_pbot(c)
    cs_sun = max(cs_sun,10.e-06_r8)

    if ( stomatalcond_mtd == stomatalcond_mtd_medlyn2011 )then
       term = 1.6_r8 * an_sun(p,iv) / (cs_sun / forc_pbot(c) * 1.e06_r8)
       aquad = 1.0_r8
       bquad = -(2.0 * (medlynintercept(patch%itype(p))*1.e-06_r8 + term) + (medlynslope(patch%itype(p)) * term)**2 / &
               (gb_mol*1.e-06_r8 * rh_can))
       cquad = medlynintercept(patch%itype(p))*medlynintercept(patch%itype(p))*1.e-12_r8 + &
               (2.0*medlynintercept(patch%itype(p))*1.e-06_r8 + term * &
               (1.0 - medlynslope(patch%itype(p))* medlynslope(patch%itype(p)) / rh_can)) * term

       call quadratic (aquad, bquad, cquad, r1, r2)
       gs_mol_sun = max(r1,r2) * 1.e06_r8
   
       ! Shaded
       cs_sha = cair - 1.4_r8/gb_mol * an_sha(p,iv) * forc_pbot(c)
       cs_sha = max(cs_sha,10.e-06_r8)
   
       term = 1.6_r8 * an_sha(p,iv) / (cs_sha / forc_pbot(c) * 1.e06_r8)
       aquad = 1.0_r8
       bquad = -(2.0 * (medlynintercept(patch%itype(p))*1.e-06_r8 + term) + (medlynslope(patch%itype(p)) * term)**2 / &
               (gb_mol*1.e-06_r8 * rh_can))
       cquad = medlynintercept(patch%itype(p))*medlynintercept(patch%itype(p))*1.e-12_r8 + &
               (2.0*medlynintercept(patch%itype(p))*1.e-06_r8 + term * (1.0 - medlynslope(patch%itype(p))* &
               medlynslope(patch%itype(p)) / rh_can)) * term

       call quadratic (aquad, bquad, cquad, r1, r2)
       gs_mol_sha = max(r1,r2)* 1.e06_r8
    else if ( stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
       aquad = cs_sun
       bquad = cs_sun*(gb_mol - max(bsun*bbb(p),1._r8)) - mbb(p)*an_sun(p,iv)*forc_pbot(c)
       cquad = -gb_mol*(cs_sun*max(bsun*bbb(p),1._r8) + mbb(p)*an_sun(p,iv)*forc_pbot(c)*rh_can)
       call quadratic (aquad, bquad, cquad, r1, r2)
       gs_mol_sun = max(r1,r2)
    
       ! Shaded
       cs_sha = cair - 1.4_r8/gb_mol * an_sha(p,iv) * forc_pbot(c)
       cs_sha = max(cs_sha,10.e-06_r8)
    
       aquad = cs_sha
       bquad = cs_sha*(gb_mol - max(bsha*bbb(p),1._r8)) - mbb(p)*an_sha(p,iv)*forc_pbot(c)
       cquad = -gb_mol*(cs_sha*max(bsha*bbb(p),1._r8) + mbb(p)*an_sha(p,iv)*forc_pbot(c)*rh_can)
       call quadratic (aquad, bquad, cquad, r1, r2)
       gs_mol_sha = max(r1,r2)
    end if
    
    ! Derive new estimate for cisun,cisha
    if (an_sun(p,iv) >= 0._r8) then
       if (gs_mol_sun > 0._r8) then
          fvalsun =cisun - cair + an_sun(p,iv) * forc_pbot(c) * (1.4_r8*gs_mol_sun+1.6_r8*gb_mol) / (gb_mol*gs_mol_sun)
       else
          fvalsun =cisun - cair
       endif
    endif
    if (an_sha(p,iv) >= 0._r8) then
       if (gs_mol_sha > 0._r8) then
          fvalsha =cisha - cair + an_sha(p,iv) * forc_pbot(c) * (1.4_r8*gs_mol_sha+1.6_r8*gb_mol) / (gb_mol*gs_mol_sha)
       else
          fvalsha =cisha - cair
       endif
    endif
    end associate
  end subroutine ci_func_PHS
  !--------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine calcstress(p,c,x,bsun,bsha,gb_mol,gs_mol_sun,gs_mol_sha,qsatl,qaf, &
       atm2lnd_inst,canopystate_inst,waterdiagnosticbulk_inst,soilstate_inst, &
       temperature_inst, waterfluxbulk_inst)
    !
    ! DESCRIPTIONS
    ! compute the transpiration stress using a plant hydraulics approach
    ! calls spacF, spacA, and getvegwp
    !
    ! USES
    use clm_varpar        , only : nlevsoi
    use clm_varcon        , only : rgas
    !!
    ! !ARGUMENTS:
    integer                , intent(in)  :: p               ! pft index
    integer                , intent(in)  :: c               ! column index
    real(r8)               , intent(inout)  :: x(nvegwcs)   ! working copy of vegwp(p,:)
    real(r8)               , intent(out) :: bsun            ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out) :: bsha            ! shaded sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(in)  :: gb_mol          ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)  :: gs_mol_sun      ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)  :: gs_mol_sha      ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)  :: qsatl           ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)  :: qaf             ! humidity of canopy air [kg/kg]
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(waterdiagnosticbulk_type)  , intent(in)  :: waterdiagnosticbulk_inst
    type(soilstate_type)   , intent(in)  :: soilstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    type(waterfluxbulk_type)   , intent(in)  :: waterfluxbulk_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: A(nvegwcs,nvegwcs)    ! matrix relating d(vegwp) and f: d(vegwp)=A*f 
    real(r8) :: f(nvegwcs)            ! flux divergence (mm/s)
    real(r8) :: dx(nvegwcs)           ! change in vegwp from one iter to the next [mm]
    real(r8) :: efpot                 ! potential latent energy flux [kg/m2/s]
    real(r8) :: rppdry_sun            ! fraction of potential evaporation through transp - sunlit [-]
    real(r8) :: rppdry_sha            ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: qflx_sun              ! [kg/m2/s]
    real(r8) :: qflx_sha              ! [kg/m2/s]
    real(r8) :: gs0sun,gs0sha         ! local gs_mol copies
    real(r8) :: qsun,qsha             ! attenuated transpiration fluxes
    integer  :: j                     ! index
    real(r8) :: cf                    ! s m**2/umol -> s/m
    integer  :: iter                  ! newton's method iteration number
    logical  :: flag                  ! signal that matrix was not invertible
    logical  :: night                 ! signal to store vegwp within this routine, b/c it is night-time and full suite won't be called
    integer, parameter  :: itmax=50   ! exit newton's method if iters>itmax
    real(r8), parameter :: tolf=1.e-6,toldx=1.e-9 !tolerances for a satisfactory solution
    logical  :: havegs                ! signals direction of calculation gs->qflx or qflx->gs 
    real(r8) :: soilflux              ! total soil column transpiration [mm/s] 
    real(r8), parameter :: tol_lai=.001_r8 ! minimum lai where transpiration is calc'd 
    !------------------------------------------------------------------------------
    
    associate(                                                    &
         laisun        => canopystate_inst%laisun_patch         , & ! Input:  [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input:  [real(r8) (:)   ]  shaded leaf area
         elai          => canopystate_inst%elai_patch           , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
         esai          => canopystate_inst%esai_patch           , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow
         tsai          => canopystate_inst%tsai_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         fdry          => waterdiagnosticbulk_inst%fdry_patch            , & ! Input:  [real(r8) (:)   ]  fraction of foliage that is green and dry [-]
         forc_rho      => atm2lnd_inst%forc_rho_downscaled_col  , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)
         forc_pbot     => atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         tgcm          => temperature_inst%thm_patch            , & ! Input:  [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)
         bsw           => soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         qflx_tran_veg => waterfluxbulk_inst%qflx_tran_veg_patch    , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         sucsat        => soilstate_inst%sucsat_col               & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         )

    !temporary flag for night time vegwp(sun)>0  
    if (x(sun)>0._r8) then
       night=.TRUE.
       x(sun)=x(sha)
    else
       night=.FALSE.
    endif
    
    !copy to avoid rewriting gs_mol_sun
    gs0sun=gs_mol_sun
    gs0sha=gs_mol_sha
    
    !compute transpiration demand
    havegs=.true.
    call getqflx(p,c,gb_mol,gs0sun,gs0sha,qflx_sun,qflx_sha,qsatl,qaf,havegs, &
         atm2lnd_inst, canopystate_inst, waterdiagnosticbulk_inst, temperature_inst)
    
    if ((laisun(p)>tol_lai .or. laisha(p)>tol_lai).and.&
         (qflx_sun>0._r8 .or. qflx_sha>0._r8))then

    !newton's method solves for matching fluxes through the spac
    iter=0
    do
       
       iter=iter+1

       call spacF(p,c,x,f,qflx_sun,qflx_sha, &
            atm2lnd_inst,canopystate_inst,soilstate_inst,temperature_inst,waterfluxbulk_inst)
          
       if ( sqrt(sum(f*f)) < tolf*(qflx_sun+qflx_sha) ) then  !fluxes balanced -> exit
          flag = .false.
          exit
       end if
       if ( iter>itmax ) then                                 !exceeds max iters -> exit
          flag = .false.
          exit
       end if
       
       call spacA(p,c,x,A,qflx_sun,qflx_sha,flag, &
            atm2lnd_inst,canopystate_inst,soilstate_inst,temperature_inst,waterfluxbulk_inst)

       if (flag) then
          ! cannot invert the matrix, solve for x algebraically assuming no flux                            
          exit
       end if

       if (laisun(p)>tol_lai.and.laisha(p)>tol_lai)then
          dx = matmul(A,f)
       else
          !reduces to 3x3 system
          !in this case, dx is not always [sun,sha,xyl,root]
          !sun and sha flip depending on which is lai==0
          dx(sun)=0._r8
          dx(sha:root)=matmul(A(sha:root,sha:root),f(sha:root))
       endif
       
       
       if ( maxval(abs(dx)) > 50000._r8) then
          dx = 50000._r8 * dx / maxval(abs(dx))  !rescale step to max of 50000
       end if


       if (laisun(p)>tol_lai.and.laisha(p)>tol_lai)then
          x=x+dx
       elseif (laisha(p)>tol_lai) then
          x=x+dx
          x(sun)=x(xyl) ! psi_sun = psi_xyl because laisun==0
       else
          x(xyl:root)=x(xyl:root)+dx(xyl:root)
          x(sun)=x(sun)+dx(sha)  ! implementation ugly bit, chose to flip dx(sun) and dx(sha) for laisha==0 case
          x(sha)=x(xyl) ! psi_sha = psi_xyl because laisha==0
         
       endif


       if ( sqrt(sum(dx*dx)) < toldx) then
          !step in vegwp small -> exit
          exit
       end if
       
       ! this is a catch to force spac gradient to atmosphere
       if ( x(xyl) > x(root) ) x(xyl) = x(root)
       if ( x(sun) > x(xyl) )  x(sun) = x(xyl)
       if ( x(sha) > x(xyl) )  x(sha) = x(xyl)
       
    end do

    else
       ! both qflxsun and qflxsha==0
	flag=.true.
    end if

    if (flag) then
       ! solve algebraically
       call getvegwp(p, c, x, gb_mol, gs0sun, gs0sha, qsatl, qaf, soilflux, &
               atm2lnd_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, temperature_inst)
       bsun = plc(x(sun),p,c,sun,veg)
       bsha = plc(x(sha),p,c,sha,veg)
    else     
    ! compute attenuated flux
    qsun=qflx_sun*plc(x(sun),p,c,sun,veg)
    qsha=qflx_sha*plc(x(sha),p,c,sha,veg)
    
    ! retrieve stressed stomatal conductance
    havegs=.FALSE.
    call getqflx(p,c,gb_mol,gs0sun,gs0sha,qsun,qsha,qsatl,qaf,havegs, &
         atm2lnd_inst, canopystate_inst, waterdiagnosticbulk_inst, temperature_inst)
    
    ! compute water stress
    ! .. generally -> B= gs_stressed / gs_unstressed
    ! .. when gs=0 -> B= plc( x )
    if (qflx_sun>0._r8) then
       bsun = gs0sun/gs_mol_sun
    else
       bsun = plc(x(sun),p,c,sun,veg)
    endif
    if (qflx_sha>0._r8) then
       bsha = gs0sha/gs_mol_sha
    else
       bsha = plc(x(sha),p,c,sha,veg)
    endif
    endif
    if ( bsun < 0.01_r8 ) bsun = 0._r8
    if ( bsha < 0.01_r8 ) bsha = 0._r8

    !zqz is this the best place to do this?
    ! was looking like qflx_tran_veg/vegwp was not being set at night time
    ! set vegwp for the final gs_mol solution
    if (night) then
       gs0sun=bsun*gs_mol_sun
       gs0sha=bsha*gs_mol_sha
       call getvegwp(p, c, x, gb_mol, gs0sun, gs0sha, qsatl, qaf, soilflux, &
            atm2lnd_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, temperature_inst)
       if (soilflux<0._r8) soilflux = 0._r8
       qflx_tran_veg(p) = soilflux
    endif
    
    
    end associate
  
  end subroutine calcstress
   
   !------------------------------------------------------------------------------
   
  !------------------------------------------------------------------------------
  subroutine spacA(p,c,x,invA,qflx_sun,qflx_sha,flag, &
       atm2lnd_inst,canopystate_inst,soilstate_inst, &
       temperature_inst, waterfluxbulk_inst)
    
    !
    ! DESCRIPTION
    !  Returns invA, the inverse matrix relating delta(vegwp) to f
    !   d(vegwp)=invA*f
    !   evaluated at vegwp(p)
    !
    ! The methodology is currently hardcoded for linear algebra assuming the
    ! number of vegetation segments is four. Thus the matrix A and it's inverse
    ! invA are both 4x4 matrices. A more general method could be done using for
    ! example a LINPACK linear algebra solver.
    !
    ! USES
    use clm_varpar        , only : nlevsoi
    use clm_varcon        , only : rgas
    !
    ! !ARGUMENTS:
    integer                , intent(in)  :: p               ! pft index
    integer                , intent(in)  :: c               ! column index
    real(r8)               , intent(in)  :: x(nvegwcs)      ! working copy of veg water potential for patch p [mm H2O] 
    real(r8)               , intent(out) :: invA(nvegwcs,nvegwcs)   ! matrix relating d(vegwp) and f: d(vegwp)=invA*f
    real(r8)               , intent(in)  :: qflx_sun        ! Sunlit leaf transpiration [kg/m2/s] 
    real(r8)               , intent(in)  :: qflx_sha        ! Shaded leaf transpiration [kg/m2/s]
    logical                , intent(out) :: flag            ! tells calling function that the matrix is not invertible
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(soilstate_type)   , intent(in)  :: soilstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    type(waterfluxbulk_type)   , intent(in)  :: waterfluxbulk_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: fsto1                 ! sunlit transpiration reduction function [-]
    real(r8) :: fsto2                 ! shaded transpiration reduction function [-] 
    real(r8) :: fx                    ! fraction of maximum conductance, xylem-to-leaf [-] 
    real(r8) :: fr                    ! fraction of maximum conductance, root-to-xylem [-] 
    real(r8) :: dfsto1                ! 1st derivative of fsto1 w.r.t. change in vegwp
    real(r8) :: dfsto2                ! 1st derivative of fsto2 w.r.t. change in vegwp
    real(r8) :: dfx                   ! 1st derivative of fx w.r.t. change in vegwp
    real(r8) :: dfr                   ! 1st derivative of fr w.r.t. change in vegwp
    real(r8) :: A(nvegwcs,nvegwcs)    ! matrix relating vegwp to flux divergence f=A*d(vegwp)
    real(r8) :: leading               ! inverse of determiniant
    real(r8) :: determ                ! determinant of matrix
    real(r8) :: grav1                 ! gravitational potential surface to canopy top (mm H2O)
    real(r8) :: invfactor             ! 
    real(r8), parameter :: tol_lai=.001_r8 ! minimum lai where transpiration is calc'd
    integer  :: j                     ! index
    !------------------------------------------------------------------------------
#ifndef NDEBUG
    ! Only execute this code if DEBUG=TRUE
    if ( nvegwcs /= 4 )then
       call endrun(msg='Error:: this function is hardcoded for 4x4 matrices with nvegwcs==4'//errMsg(__FILE__, __LINE__))
    end if
#endif
    
    associate(                                                    &
         k_soil_root  => soilstate_inst%k_soil_root_patch       , & ! Input:  [real(r8) (:,:) ]  soil-root interface conductance (mm/s)
         laisun        => canopystate_inst%laisun_patch         , & ! Input:  [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input:  [real(r8) (:)   ]  shaded leaf area
         htop          => canopystate_inst%htop_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy top (m)
         tsai          => canopystate_inst%tsai_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         ivt           => patch%itype                             & ! Input:  [integer  (:)   ]  patch vegetation type
         )
    
    ! initialize all elements to zero
    A = 0._r8
    invA = 0._r8

    grav1 = htop(p)*1000._r8
    
    !compute conductance attentuation for each segment
    fsto1=  plc(x(sun),p,c,sun,veg)
    fsto2=  plc(x(sha),p,c,sha,veg)
    fx=     plc(x(xyl),p,c,xyl,veg)
    fr=     plc(x(root),p,c,root,veg)
    
    !compute 1st deriv of conductance attenuation for each segment
    dfsto1=  d1plc(x(sun),p,c,sun,veg)
    dfsto2=  d1plc(x(sha),p,c,sha,veg)
    dfx=     d1plc(x(xyl),p,c,xyl,veg)
    dfr=     d1plc(x(root),p,c,root,veg)
    
    !A - f=A*d(vegwp)
    A(1,1)= - laisun(p) * params_inst%kmax(ivt(p),sun) * fx&
         - qflx_sun * dfsto1
    A(1,3)= laisun(p) * params_inst%kmax(ivt(p),sun) * dfx * (x(xyl)-x(sun))&
         + laisun(p) * params_inst%kmax(ivt(p),sun) * fx
    A(2,2)= - laisha(p) * params_inst%kmax(ivt(p),sha) * fx&
         - qflx_sha * dfsto2
    A(2,3)= laisha(p) * params_inst%kmax(ivt(p),sha) * dfx * (x(xyl)-x(sha))&
         + laisha(p) * params_inst%kmax(ivt(p),sha) * fx
    A(3,1)= laisun(p) * params_inst%kmax(ivt(p),sun) * fx
    A(3,2)= laisha(p) * params_inst%kmax(ivt(p),sha) * fx
    A(3,3)= - laisun(p) * params_inst%kmax(ivt(p),sun) * dfx * (x(xyl)-x(sun))&
         - laisun(p) * params_inst%kmax(ivt(p),sun) * fx&
         - laisha(p) * params_inst%kmax(ivt(p),sha) * dfx * (x(xyl)-x(sha))&
         - laisha(p) * params_inst%kmax(ivt(p),sha) * fx&
         - tsai(p) * params_inst%kmax(ivt(p),xyl) / htop(p) * fr
    A(3,4)= tsai(p) * params_inst%kmax(ivt(p),xyl) / htop(p) * dfr * (x(root)-x(xyl)-grav1)&
         + tsai(p) * params_inst%kmax(ivt(p),xyl) / htop(p) * fr
    A(4,3)= tsai(p) * params_inst%kmax(ivt(p),xyl) / htop(p) * fr
    A(4,4)= - tsai(p) * params_inst%kmax(ivt(p),xyl) / htop(p) * fr&
         - tsai(p) * params_inst%kmax(ivt(p),xyl) / htop(p) * dfr * (x(root)-x(xyl)-grav1)&
         - sum(k_soil_root(p,1:nlevsoi))

    invfactor=1._r8
    A=invfactor*A

    !matrix inversion
    if (laisun(p)>tol_lai .and. laisha(p)>tol_lai) then
       ! general case

       determ=A(4,4)*A(2,2)*A(3,3)*A(1,1) - A(4,4)*A(2,2)*A(3,1)*A(1,3)&
            - A(4,4)*A(3,2)*A(2,3)*A(1,1) - A(4,3)*A(1,1)*A(2,2)*A(3,4)
       if ( abs(determ) <= 1.e-50_r8 ) then
          flag = .true.  !tells calling function that the matrix is not invertible
          return
       else
          flag = .false.
       end if       
    
       leading = 1._r8/determ

       !algebraic inversion of the matrix
       invA(1,1)=leading*A(4,4)*A(2,2)*A(3,3) - leading*A(4,4)*A(3,2)*A(2,3) - leading*A(4,3)*A(2,2)*A(3,4)
       invA(2,1)=leading*A(2,3)*A(4,4)*A(3,1)
       invA(3,1)=-leading*A(4,4)*A(2,2)*A(3,1)
       invA(4,1)=leading*A(4,3)*A(2,2)*A(3,1)
       invA(1,2)=leading*A(1,3)*A(4,4)*A(3,2)
       invA(2,2)=leading*A(4,4)*A(3,3)*A(1,1)-leading*A(4,4)*A(3,1)*A(1,3)-leading*A(4,3)*A(1,1)*A(3,4)
       invA(3,2)=-leading*A(1,1)*A(4,4)*A(3,2)
       invA(4,2)=leading*A(4,3)*A(1,1)*A(3,2)
       invA(1,3)=-leading*A(1,3)*A(2,2)*A(4,4)
       invA(2,3)=-leading*A(2,3)*A(1,1)*A(4,4)
       invA(3,3)=leading*A(2,2)*A(1,1)*A(4,4)
       invA(4,3)=-leading*A(4,3)*A(1,1)*A(2,2)
       invA(1,4)=leading*A(1,3)*A(3,4)*A(2,2)
       invA(2,4)=leading*A(2,3)*A(3,4)*A(1,1)
       invA(3,4)=-leading*A(3,4)*A(1,1)*A(2,2)
       invA(4,4)=leading*A(2,2)*A(3,3)*A(1,1)-leading*A(2,2)*A(3,1)*A(1,3)-leading*A(3,2)*A(2,3)*A(1,1)
       invA=invfactor*invA !undo inversion scaling
    else
       ! if laisun or laisha ==0 invert the corresponding 3x3 matrix
       ! if both are zero, this routine is not called
       if (laisha(p)<=tol_lai) then
          ! shift nonzero matrix values so that both 3x3 cases can be inverted with the same code
          A(2,2)=A(1,1)
          A(3,2)=A(3,1)
          A(2,3)=A(1,3)
       endif
       determ=A(2,2)*A(3,3)*A(4,4)-A(3,4)*A(2,2)*A(4,3)-A(2,3)*A(3,2)*A(4,4)
       if ( abs(determ) <= 1.e-50_r8 ) then
          flag = .true.  !tells calling function that the matrix is not invertible
          return
       else
          flag = .false.
       end if
       
       !algebraic inversion of the 3x3 matrix stored in A(2:4,2:4)
       invA(2,2)=A(3,3)*A(4,4)-A(3,4)*A(4,3)
       invA(2,3)=-A(2,3)*A(4,4)
       invA(2,4)=A(3,4)*A(2,3)
       invA(3,2)=-A(3,2)*A(4,4)
       invA(3,3)=A(2,2)*A(4,4)
       invA(3,4)=-A(3,4)*A(2,2)
       invA(4,2)=A(3,2)*A(4,3)
       invA(4,3)=-A(2,2)*A(4,3)
       invA(4,4)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
       invA=1._r8/determ*invA
       
    endif

    end associate
    
  end subroutine spacA
  
  !--------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine spacF(p,c,x,f,qflx_sun,qflx_sha, &
       atm2lnd_inst,canopystate_inst,soilstate_inst, &
       temperature_inst, waterfluxbulk_inst)
    !
    ! DESCRIPTION
    ! Returns f, the flux divergence across each vegetation segment
    !  calculated for vegwp(p,:) as passed in via x
    !
    ! USES
    use clm_varpar        , only : nlevsoi
    use clm_varcon        , only : rgas
    use ColumnType        , only : col
    !
    ! !ARGUMENTS:
    integer                , intent(in)  :: p               ! pft index
    integer                , intent(in)  :: c               ! column index
    real(r8)               , intent(in)  :: x(nvegwcs)      ! working copy of veg water potential for patch p [mm H2O]
    real(r8)               , intent(out) :: f(nvegwcs)      ! water flux divergence [mm/s]
    real(r8)               , intent(in)  :: qflx_sun        ! Sunlit leaf transpiration [kg/m2/s] 
    real(r8)               , intent(in)  :: qflx_sha        ! Shaded leaf transpiration [kg/m2/s] 
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(soilstate_type)   , intent(in)  :: soilstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    type(waterfluxbulk_type)   , intent(in)  :: waterfluxbulk_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: fsto1                 ! sunlit transpiration reduction function [-]
    real(r8) :: fsto2                 ! shaded transpiration reduction function [-]
    real(r8) :: fx                    ! fraction of maximum conductance, xylem-to-leaf [-] 
    real(r8) :: fr                    ! fraction of maximum conductance, root-to-xylem [-]
    real(r8) :: grav1                 ! gravitational potential surface to canopy top (mm H2O) 
    real(r8) :: grav2(nlevsoi)        ! soil layer gravitational potential relative to surface (mm H2O) 
    real(r8) :: temp                  ! used to copy f(sun) to f(sha) for special case
    real(r8), parameter :: tol_lai=.001_r8  ! needs to be the same as in calcstress and spacA (poor form, refactor)<
    integer  :: j                     ! index
    !------------------------------------------------------------------------------
    
    associate(                                              &
         k_soil_root  => soilstate_inst%k_soil_root_patch       , & ! Input:  [real(r8) (:,:) ]  soil-root interface conductance (mm/s)
         laisun        => canopystate_inst%laisun_patch         , & ! Input:  [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input:  [real(r8) (:)   ]  shaded leaf area
         htop          => canopystate_inst%htop_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy top (m)
         tsai          => canopystate_inst%tsai_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         smp           => soilstate_inst%smp_l_col              , & ! Input: [real(r8) (:,:) ]  soil matrix potential [mm]
         ivt           => patch%itype                           , & ! Input:  [integer  (:)   ]  patch vegetation type
         qflx_tran_veg => waterfluxbulk_inst%qflx_tran_veg_patch    , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         z             => col%z                                   & ! Input:  [real(r8) (:,:) ]  layer node depth (m)
         )
    
    grav1 = htop(p) * 1000._r8
    grav2(1:nlevsoi) = z(c,1:nlevsoi) * 1000._r8
    
    fsto1=  plc(x(sun),p,c,sun,veg)
    fsto2=  plc(x(sha),p,c,sha,veg)
    fx=     plc(x(xyl),p,c,xyl,veg)
    fr=     plc(x(root),p,c,root,veg)
    
    !compute flux divergence across each plant segment
    f(sun)= qflx_sun * fsto1 - laisun(p) * params_inst%kmax(ivt(p),sun) * fx * (x(xyl)-x(sun))
    f(sha)= qflx_sha * fsto2 - laisha(p) * params_inst%kmax(ivt(p),sha) * fx * (x(xyl)-x(sha))
    f(xyl)= laisun(p) * params_inst%kmax(ivt(p),sun) * fx * (x(xyl)-x(sun))&
         + laisha(p) * params_inst%kmax(ivt(p),sha) * fx * (x(xyl)-x(sha)) &
         - tsai(p) * params_inst%kmax(ivt(p),xyl) / htop(p) * fr * (x(root)-x(xyl)-grav1)
    f(root)= tsai(p) * params_inst%kmax(ivt(p),xyl) / htop(p) * fr * (x(root)-x(xyl)-grav1) &
         + sum( k_soil_root(p,1:nlevsoi) * (x(root)+grav2(1:nlevsoi)) ) &
         - sum( k_soil_root(p,1:nlevsoi) * smp(c,1:nlevsoi) )

    if (laisha(p)<tol_lai) then
       ! special case for laisha ~ 0
       ! flip sunlit and shade fluxes to match special case handling in spacA
       temp=f(sun)
       f(sun)=f(sha)
       f(sha)=temp
    endif

    end associate

  end subroutine spacF
  
  !--------------------------------------------------------------------------------
  subroutine getvegwp(p, c, x, gb_mol, gs_mol_sun, gs_mol_sha, qsatl, qaf, soilflux, &
       atm2lnd_inst, canopystate_inst, waterdiagnosticbulk_inst, soilstate_inst, temperature_inst)
    ! !DESCRIPTION:
    !  Calculates transpiration and returns corresponding vegwp in x
    !
    ! !USES:
    ! calls getqflx
    use clm_varpar  , only : nlevsoi
    use ColumnType  , only : col
    implicit none
    !
    ! !ARGUMENTS:
    integer                , intent(in)  :: p                ! pft index
    integer                , intent(in)  :: c                ! column index
    real(r8)               , intent(out) :: x(nvegwcs)       ! working copy of veg water potential for patch p
    real(r8)               , intent(in)  :: gb_mol           ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(inout)  :: gs_mol_sun    ! Ball-Berry leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(inout)  :: gs_mol_sha    ! Ball-Berry leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)  :: qsatl            ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)  :: qaf              ! humidity of canopy air [kg/kg]
    real(r8)               , intent(out) :: soilflux         ! total soil column transpiration [mm/s]
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(waterdiagnosticbulk_type)  , intent(in)  :: waterdiagnosticbulk_inst
    type(soilstate_type)   , intent(in)  :: soilstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: qflx_sun                 ! Sunlit leaf transpiration [kg/m2/s]
    real(r8) :: qflx_sha                 ! Shaded leaf transpiration [kg/m2/s] 
    real(r8) :: fx                       ! fraction of maximum conductance, xylem-to-leaf [-]  
    real(r8) :: fr                       ! fraction of maximum conductance, root-to-xylem [-]  
    real(r8) :: grav1                    ! gravitational potential surface to canopy top (mm H2O)
    real(r8) :: grav2(nlevsoi)           ! soil layer gravitational potential relative to surface (mm H2O) 
    integer  :: j                        ! index
    logical  :: havegs                   ! signals direction of calculation gs->qflx or qflx->gs 
    !----------------------------------------------------------------------
    associate(                                                    &
         k_soil_root  => soilstate_inst%k_soil_root_patch       , & ! Input:  [real(r8) (:,:) ]  soil-root interface conductance (mm/s)
         laisun        => canopystate_inst%laisun_patch         , & ! Input: [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input: [real(r8) (:)   ]  shaded leaf area
         htop          => canopystate_inst%htop_patch           , & ! Input: [real(r8) (:)   ]  patch canopy top (m)
         tsai          => canopystate_inst%tsai_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         smp           => soilstate_inst%smp_l_col              , & ! Input: [real(r8) (:,:) ]  soil matrix potential [mm]
         rootfr        => soilstate_inst%rootfr_patch           , & ! Input: [real(r8) (:,:) ]  fraction of roots in each soil layer
         bsw           => soilstate_inst%bsw_col                , & ! Input: [real(r8) (:,:) ]  Clapp and Hornberger "b"
         ivt           => patch%itype                           , & ! Input: [integer  (:)   ]  patch vegetation type
         hk_l          => soilstate_inst%hk_l_col               , & ! Input: [real(r8) (:,:) ]  hydraulic conductivity (mm/s)
         hksat         => soilstate_inst%hksat_col              , & ! Input: [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         sucsat        => soilstate_inst%sucsat_col             , & ! Input: [real(r8) (:,:) ]  minimum soil suction (mm)
         z             => col%z                                   & ! Input: [real(r8) (:,:) ]  layer node depth (m)
         )
    
    grav1 = 1000._r8 *htop(p)
    grav2(1:nlevsoi) = 1000._r8 * z(c,1:nlevsoi)
    
    !compute transpiration demand
    havegs=.true.
    call getqflx(p,c,gb_mol,gs_mol_sun,gs_mol_sha,qflx_sun,qflx_sha,qsatl,qaf,havegs, &
         atm2lnd_inst, canopystate_inst, waterdiagnosticbulk_inst, temperature_inst)
    
    !calculate root water potential
    if ( abs(sum(k_soil_root(p,1:nlevsoi))) == 0._r8 ) then
       x(root) = sum(smp(c,1:nlevsoi) - grav2)/nlevsoi
    else
       x(root) = (sum(k_soil_root(p,1:nlevsoi)*(smp(c,1:nlevsoi)-grav2))-qflx_sun-qflx_sha) &
                  /sum(k_soil_root(p,1:nlevsoi))
    endif
    
    !calculate xylem water potential
    fr = plc(x(root),p,c,root,veg)
    if ( (tsai(p) > 0._r8) .and. (fr > 0._r8) ) then
       x(xyl) = x(root) - grav1 - (qflx_sun+qflx_sha)/(fr*params_inst%kmax(ivt(p),root)/htop(p)*tsai(p))!removed htop conversion
    else
       x(xyl) = x(root) - grav1
    endif
    
    !calculate sun/sha leaf water potential
    fx = plc(x(xyl),p,c,xyl,veg)
    if ( (laisha(p) > 0._r8) .and. (fx > 0._r8) ) then
       x(sha) = x(xyl) - (qflx_sha/(fx*params_inst%kmax(ivt(p),xyl)*laisha(p)))
    else
       x(sha) = x(xyl)
    endif
    if ( (laisun(p) > 0._r8) .and. (fx > 0._r8) ) then
       x(sun) = x(xyl) - (qflx_sun/(fx*params_inst%kmax(ivt(p),xyl)*laisun(p)))
    else
       x(sun) = x(xyl)
    endif

    !calculate soil flux
    soilflux = 0._r8
    do j = 1,nlevsoi
       soilflux = soilflux + k_soil_root(p,j)*(smp(c,j)-x(root)-grav2(j))
    enddo

    end associate

  end subroutine getvegwp
  
  !--------------------------------------------------------------------------------
  subroutine getqflx(p,c,gb_mol,gs_mol_sun,gs_mol_sha,qflx_sun,qflx_sha,qsatl,qaf,havegs, &
       atm2lnd_inst, canopystate_inst, waterdiagnosticbulk_inst, temperature_inst)
    ! !DESCRIPTION:
    !  calculate sunlit and shaded transpiration using gb_MOL and gs_MOL
    ! !USES:
    !
    use clm_varcon        , only : rgas
    implicit none
    !
    ! !ARGUMENTS:
    integer  , intent(in)     :: p          ! pft index
    integer  , intent(in)     :: c          ! column index
    real(r8) , intent(in)     :: gb_mol     ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8) , intent(inout)  :: gs_mol_sun ! Ball-Berry leaf conductance (umol H2O/m**2/s)
    real(r8) , intent(inout)  :: gs_mol_sha ! Ball-Berry leaf conductance (umol H2O/m**2/s)
    real(r8) , intent(inout)  :: qflx_sun   ! Sunlit leaf transpiration [kg/m2/s]
    real(r8) , intent(inout)  :: qflx_sha   ! Shaded leaf transpiration [kg/m2/s]
    real(r8) , intent(in)     :: qsatl      ! leaf specific humidity [kg/kg]
    real(r8) , intent(in)     :: qaf        ! humidity of canopy air [kg/kg]
    logical  , intent(in)     :: havegs     ! signals direction of calculation gs->qflx or qflx->gs
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(waterdiagnosticbulk_type)  , intent(in)  :: waterdiagnosticbulk_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                      ! heat conductance for leaf [m/s]
    real(r8) :: efpot                    ! potential latent energy flux [kg/m2/s]
    real(r8) :: rppdry_sun               ! fraction of potential evaporation through transp - sunlit [-]
    real(r8) :: rppdry_sha               ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: cf                       ! s m**2/umol -> s/m
    !----------------------------------------------------------------------
    
    associate(                                                    &
         laisun        => canopystate_inst%laisun_patch         , & ! Input: [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input: [real(r8) (:)   ]  shaded leaf area
         elai          => canopystate_inst%elai_patch           , & ! Input: [real(r8) (:)   ]  one-sided leaf area index with burying by snow
         esai          => canopystate_inst%esai_patch           , & ! Input: [real(r8) (:)   ]  one-sided stem area index with burying by snow
         fdry          => waterdiagnosticbulk_inst%fdry_patch            , & ! Input: [real(r8) (:)   ]  fraction of foliage that is green and dry [-]
         forc_rho      => atm2lnd_inst%forc_rho_downscaled_col  , & ! Input: [real(r8) (:)   ]  density (kg/m**3)
         forc_pbot     => atm2lnd_inst%forc_pbot_downscaled_col , & ! Input: [real(r8) (:)   ]  atmospheric pressure (Pa)
         tgcm          => temperature_inst%thm_patch              & ! Input: [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)
         )
    
    
    cf       = forc_pbot(c)/(rgas*1.e-3_r8*tgcm(p))*1.e6_r8  ! gb->gbmol conversion factor
    wtl      = (elai(p)+esai(p))*gb_mol
    efpot    = forc_rho(c)*wtl*(qsatl-qaf)
    if (havegs) then

       if ( (efpot > 0._r8) .and. (elai(p) > 0._r8) ) then
          if ( gs_mol_sun > 0._r8 ) then
             rppdry_sun = fdry(p)/gb_mol*(laisun(p)/(1._r8/gb_mol+1._r8/gs_mol_sun))/elai(p)
             qflx_sun   = efpot*rppdry_sun/cf
          else
             qflx_sun   = 0._r8
          end if
          if ( gs_mol_sha > 0._r8 ) then
             rppdry_sha = fdry(p)/gb_mol*(laisha(p)/(1._r8/gb_mol+1._r8/gs_mol_sha))/elai(p)
             qflx_sha   = efpot*rppdry_sha/cf
          else
             qflx_sha   = 0._r8
          end if
       else
          qflx_sun      = 0._r8
          qflx_sha      = 0._r8
       end if
       
    else
       if (qflx_sun > 0._r8) then
          gs_mol_sun=gb_mol*qflx_sun*cf*elai(p)/(efpot*fdry(p)*laisun(p)-qflx_sun*cf*elai(p))
       else
          gs_mol_sun=0._r8
       endif
       if (qflx_sha > 0._r8) then
          gs_mol_sha=gb_mol*qflx_sha*cf*elai(p)/(efpot*fdry(p)*laisha(p)-qflx_sha*cf*elai(p))
       else
          gs_mol_sha=0._r8
       endif
       
    endif

    end associate

  end subroutine getqflx

  !--------------------------------------------------------------------------------
  function plc(x,p,c,level,plc_method)
    ! !DESCRIPTION
    ! Return value of vulnerability curve at x
    !
    ! !ARGUMENTS
    real(r8) , intent(in)  :: x             ! water potential input
    integer  , intent(in)  :: p             ! index for pft
    integer  , intent(in)  :: c             ! index for column
    integer  , intent(in)  :: level         ! veg segment lvl (1:nvegwcs) 
    integer  , intent(in)  :: plc_method    !
    real(r8)               :: plc           ! attenuated conductance [0:1] 0=no flow
    !
    ! !PARAMETERS
    integer , parameter :: vegetation_weibull=0  ! case number
    !------------------------------------------------------------------------------
    associate(                                                    &
         ivt  => patch%itype                             & ! Input: [integer  (:)   ]  patch vegetation type
             )
    
    select case (plc_method)
       !possible to add other methods later
    case (vegetation_weibull)
       plc=2._r8**(-(x/params_inst%psi50(ivt(p),level))**params_inst%ck(ivt(p),level))
       if ( plc < 0.005_r8) plc = 0._r8
    case default
       print *,'must choose plc method'
    end select

    end associate
    
  end function plc
  !--------------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------------
  function d1plc(x,p,c,level,plc_method)
    ! !DESCRIPTION
    ! Return 1st derivative of vulnerability curve at x
    !
    ! !ARGUMENTS
    real(r8) , intent(in) :: x                ! water potential input
    integer  , intent(in) :: p                ! index for pft
    integer  , intent(in) :: c                ! index for column
    integer  , intent(in) :: level            ! veg segment lvl (1:nvegwcs)
    integer  , intent(in) :: plc_method       ! 0 for vegetation, 1 for soil
    real(r8)              :: d1plc            ! first deriv of plc curve at x
    !
    ! !PARAMETERS
    integer , parameter :: vegetation_weibull=0  ! case number
    !------------------------------------------------------------------------------
    associate(                                                    &
         ivt           => patch%itype                             & ! Input: [integer  (:)   ]  patch vegetation type
             )

    select case (plc_method)
       !possible to add other methods later
    case (vegetation_weibull)
       d1plc= -params_inst%ck(ivt(p),level) * log(2._r8) * (2._r8**(-(x/params_inst%psi50(ivt(p),level)) &
              **params_inst%ck(ivt(p),level))) &
              * ((x/params_inst%psi50(ivt(p),level))**params_inst%ck(ivt(p),level)) / x
    case default
       print *,'must choose plc method'
    end select

    end associate
    
  end function d1plc  
  
end module PhotosynthesisMod
