module MLCanopyFluxesType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Multilayer canopy module data structure
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varcon, only : ispval, spval
  use clm_varpar, only : nlevgrnd, numrad
  use decompMod , only : bounds_type
  use shr_kind_mod, only : r8 => shr_kind_r8
  use MLclm_varpar, only : nlevmlcan, nleaf
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC DATA TYPES:

  type, public :: mlcanopy_type

    ! ------------------------------------------------------------------------------------
    ! Variables have several naming conventions:
    !
    ! var_canopy  = single-level canopy variable
    ! var_soil    = single-level soil variable
    ! var_forcing = single-level atmospheric forcing variable
    ! var_profile = multi-level variable at each canopy layer
    ! var_leaf    = multi-level variable at each canopy layer for sunlit and shaded leaves
    ! ------------------------------------------------------------------------------------

    ! Vegetation input variables: dimension is (patch)

    real(r8), pointer :: ztop_canopy(:)          ! Canopy height (m)
    real(r8), pointer :: lai_canopy(:)           ! Leaf area index of canopy (m2/m2)
    real(r8), pointer :: sai_canopy(:)           ! Stem area index of canopy (m2/m2)
    real(r8), pointer :: root_biomass_canopy(:)  ! Fine root biomass (g biomass / m2)

    ! Atmospheric forcing variables required as input: dimension is (patch)

    real(r8), pointer :: zref_forcing(:)         ! Atmospheric reference height (m)
    real(r8), pointer :: tref_forcing(:)         ! Air temperature at reference height (K)
    real(r8), pointer :: qref_forcing(:)         ! Specific humidity at reference height (kg/kg)
    real(r8), pointer :: uref_forcing(:)         ! Wind speed at reference height (m/s)
    real(r8), pointer :: pref_forcing(:)         ! Air pressure at reference height (Pa)
    real(r8), pointer :: co2ref_forcing(:)       ! Atmospheric CO2 at reference height (umol/mol)
    real(r8), pointer :: o2ref_forcing(:)        ! Atmospheric O2 at reference height (mmol/mol)
    real(r8), pointer :: swskyb_forcing(:,:)     ! Atmospheric direct beam solar radiation (W/m2) [for numrad wavebands]
    real(r8), pointer :: swskyd_forcing(:,:)     ! Atmospheric diffuse solar radiation (W/m2) [for numrad wavebands]
    real(r8), pointer :: lwsky_forcing(:)        ! Atmospheric longwave radiation (W/m2)
    real(r8), pointer :: qflx_rain_forcing(:)    ! Rainfall (mm H2O/s = kg H2O/m2/s)
    real(r8), pointer :: qflx_snow_forcing(:)    ! Snowfall (mm H2O/s = kg H2O/m2/s)
    real(r8), pointer :: tacclim_forcing(:)      ! Average air temperature for acclimation (K)

    ! Derived atmospheric forcing variables: dimension is (patch)

    real(r8), pointer :: eref_forcing(:)         ! Vapor pressure at reference height (Pa)
    real(r8), pointer :: thref_forcing(:)        ! Atmospheric potential temperature at reference height (K)
    real(r8), pointer :: thvref_forcing(:)       ! Atmospheric virtual potential temperature at reference height (K)
    real(r8), pointer :: rhoair_forcing(:)       ! Air density at reference height (kg/m3)
    real(r8), pointer :: rhomol_forcing(:)       ! Molar density at reference height (mol/m3)
    real(r8), pointer :: mmair_forcing(:)        ! Molecular mass of air at reference height (kg/mol)
    real(r8), pointer :: cpair_forcing(:)        ! Specific heat of air (constant pressure) at reference height (J/mol/K)
    real(r8), pointer :: solar_zen_forcing(:)    ! Solar zenith angle (radians)

    ! Canopy flux variables (fluxes are per m2 ground area): dimension is (patch)

    real(r8), pointer :: swveg_canopy(:,:)       ! Absorbed solar radiation: vegetation (W/m2) [for numrad wavebands]
    real(r8), pointer :: swvegsun_canopy(:,:)    ! Absorbed solar radiation: sunlit canopy (W/m2) [for numrad wavebands]
    real(r8), pointer :: swvegsha_canopy(:,:)    ! Absorbed solar radiation: shaded canopy (W/m2) [for numrad wavebands]

    real(r8), pointer :: lwveg_canopy(:)         ! Absorbed longwave radiation: vegetation (W/m2)
    real(r8), pointer :: lwvegsun_canopy(:)      ! Absorbed longwave radiation: sunlit canopy (W/m2)
    real(r8), pointer :: lwvegsha_canopy(:)      ! Absorbed longwave radiation: shaded canopy (W/m2)

    real(r8), pointer :: shveg_canopy(:)         ! Sensible heat flux: vegetation (W/m2)
    real(r8), pointer :: shvegsun_canopy(:)      ! Sensible heat flux: sunlit canopy (W/m2)
    real(r8), pointer :: shvegsha_canopy(:)      ! Sensible heat flux: shaded canopy (W/m2)

    real(r8), pointer :: lhveg_canopy(:)         ! Latent heat flux: vegetation (W/m2)
    real(r8), pointer :: lhvegsun_canopy(:)      ! Latent heat flux: sunlit canopy (W/m2)
    real(r8), pointer :: lhvegsha_canopy(:)      ! Latent heat flux: shaded canopy (W/m2)

    real(r8), pointer :: etveg_canopy(:)         ! Water vapor flux: vegetation (mol H2O/m2/s)
    real(r8), pointer :: etvegsun_canopy(:)      ! Water vapor flux: sunlit canopy (mol H2O/m2/s)
    real(r8), pointer :: etvegsha_canopy(:)      ! Water vapor flux: shaded canopy (mol H2O/m2/s)

    real(r8), pointer :: gppveg_canopy(:)        ! Gross primary production: vegetation (umol CO2/m2/s)
    real(r8), pointer :: gppvegsun_canopy(:)     ! Gross primary production: sunlit canopy (umol CO2/m2/s)
    real(r8), pointer :: gppvegsha_canopy(:)     ! Gross primary production: shaded canopy (umol CO2/m2/s)

    real(r8), pointer :: albcan_canopy(:,:)      ! Albedo above canopy (-) [for numrad wavebands]
    real(r8), pointer :: lwup_canopy(:)          ! Upward longwave radiation above canopy (W/m2)
    real(r8), pointer :: rnet_canopy(:)          ! Total net radiation, including soil (W/m2)
    real(r8), pointer :: shflx_canopy(:)         ! Total sensible heat flux, including soil (W/m2)
    real(r8), pointer :: lhflx_canopy(:)         ! Total latent heat flux, including soil (W/m2)
    real(r8), pointer :: etflx_canopy(:)         ! Total water vapor flux, including soil (mol H2O/m2/s)
    real(r8), pointer :: stflx_canopy(:)         ! Canopy storage heat flux (W/m2)
    real(r8), pointer :: ustar_canopy(:)         ! Friction velocity (m/s)
    real(r8), pointer :: gac_to_hc_canopy(:)     ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)

    real(r8), pointer :: qflx_intr_canopy(:)     ! Intercepted precipitation (kg H2O/m2/s)
    real(r8), pointer :: qflx_tflrain_canopy(:)  ! Total rain throughfall onto ground (kg H2O/m2/s)
    real(r8), pointer :: qflx_tflsnow_canopy(:)  ! Total snow throughfall onto ground (kg H2O/m2/s)

    ! Canopy diagnostic variables: dimension is (patch)

    real(r8), pointer :: uaf_canopy(:)           ! Wind speed at canopy top (m/s)
    real(r8), pointer :: taf_canopy(:)           ! Air temperature at canopy top (K)
    real(r8), pointer :: qaf_canopy(:)           ! Specific humidity at canopy top (kg/kg)
    real(r8), pointer :: fracminlwp_canopy(:)    ! Fraction of canopy that is water-stressed

    ! Canopy aerodynamic variables: dimension is (patch)

    real(r8), pointer :: obu_canopy(:)           ! Obukhov length (m)
    real(r8), pointer :: obuold_canopy(:)        ! Obukhov length from previous iteration
    integer,  pointer :: nmozsgn_canopy(:)       ! Number of times stability changes sign during iteration
    real(r8), pointer :: beta_canopy(:)          ! Value of u* / u at canopy top (-)
    real(r8), pointer :: PrSc_canopy(:)          ! Prandtl (Schmidt) number at canopy top (-)
    real(r8), pointer :: Lc_canopy(:)            ! Canopy density length scale (m)
    real(r8), pointer :: zdisp_canopy(:)         ! Displacement height (m)

    ! Canopy stomatal conductance variables: dimension is (patch)

    real(r8), pointer :: g0_canopy(:)            ! Ball-Berry or Medlyn minimum leaf conductance (mol H2O/m2/s)
    real(r8), pointer :: g1_canopy(:)            ! Ball-Berry or Medlyn slope parameter

    ! Soil energy balance variables: dimension is (patch)

    real(r8), pointer :: albsoib_soil(:,:)       ! Direct beam albedo of ground (-) [for numrad wavebands]
    real(r8), pointer :: albsoid_soil(:,:)       ! Diffuse albedo of ground (-) [for numrad wavebands]
    real(r8), pointer :: swsoi_soil(:,:)         ! Absorbed solar radiation: ground (W/m2) [for numrad wavebands]
    real(r8), pointer :: lwsoi_soil(:)           ! Absorbed longwave radiation: ground (W/m2)
    real(r8), pointer :: rnsoi_soil(:)           ! Net radiation: ground (W/m2)
    real(r8), pointer :: shsoi_soil(:)           ! Sensible heat flux: ground (W/m2)
    real(r8), pointer :: lhsoi_soil(:)           ! Latent heat flux: ground (W/m2)
    real(r8), pointer :: etsoi_soil(:)           ! Water vapor flux: ground (mol H2O/m2/s)
    real(r8), pointer :: gsoi_soil(:)            ! Soil heat flux (W/m2)
    real(r8), pointer :: tg_soil(:)              ! Soil surface temperature (K)
    real(r8), pointer :: tg_bef_soil(:)          ! Soil surface temperature for previous timestep (K)
    real(r8), pointer :: eg_soil(:)              ! Soil surface vapor pressure (Pa)
    real(r8), pointer :: rhg_soil(:)             ! Relative humidity of airspace at soil surface (fraction)
    real(r8), pointer :: gac0_soil(:)            ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    real(r8), pointer :: soil_t_soil(:)          ! Temperature of first snow/soil layer (K)
    real(r8), pointer :: soil_dz_soil(:)         ! Depth to temperature of first snow/soil layer (m)
    real(r8), pointer :: soil_tk_soil(:)         ! Thermal conductivity of first snow/soil layer (W/m/K)
    real(r8), pointer :: soilres_soil(:)         ! Soil evaporative resistance (s/m)

    ! Soil moisture variables: dimension is (patch)

    real(r8), pointer :: btran_soil(:)           ! Soil wetness factor for stomatal conductance (-)
    real(r8), pointer :: psis_soil(:)            ! Weighted soil water potential (MPa)
    real(r8), pointer :: rsoil_soil(:)           ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
    real(r8), pointer :: soil_et_loss_soil(:,:)  ! Fraction of total transpiration from each soil layer (-) [for nlevgrnd layers]

    ! Canopy layer indices: dimension is (patch)

    integer , pointer :: ncan_canopy(:)          ! Number of layers
    integer , pointer :: ntop_canopy(:)          ! Index for top leaf layer
    integer , pointer :: nbot_canopy(:)          ! Index for bottom leaf layer

    ! Canopy layer variables (fluxes are per m2 ground area): dimension is (patch, level)

    real(r8), pointer :: dlai_frac_profile(:,:)  ! Canopy layer leaf area index (fraction of canopy total)
    real(r8), pointer :: dsai_frac_profile(:,:)  ! Canopy layer stem area index (fraction of canopy total)
    real(r8), pointer :: dlai_profile(:,:)       ! Canopy layer leaf area index (m2/m2)
    real(r8), pointer :: dsai_profile(:,:)       ! Canopy layer stem area index (m2/m2)
    real(r8), pointer :: dpai_profile(:,:)       ! Canopy layer plant area index (m2/m2)
    real(r8), pointer :: sumpai_profile(:,:)     ! Canopy layer cumulative plant area index (m2/m2)
    real(r8), pointer :: zs_profile(:,:)         ! Canopy layer height for scalar concentration and source (m)
    real(r8), pointer :: dz_profile(:,:)         ! Canopy layer thickness (m)

    real(r8), pointer :: vcmax25_profile(:,:)    ! Canopy layer leaf maximum carboxylation rate at 25C (umol/m2/s)
    real(r8), pointer :: jmax25_profile(:,:)     ! Canopy layer C3 maximum electron transport rate at 25C (umol/m2/s)
    real(r8), pointer :: kp25_profile(:,:)       ! Canopy layer C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    real(r8), pointer :: rd25_profile(:,:)       ! Canopy layer leaf respiration rate at 25C (umol/m2/s)
    real(r8), pointer :: cpleaf_profile(:,:)     ! Canopy layer leaf heat capacity (J/m2 leaf/K)

    real(r8), pointer :: fracsun_profile(:,:)    ! Canopy layer sunlit fraction (-)
    real(r8), pointer :: tb_profile(:,:)         ! Canopy layer transmittance of direct beam radiation (-)
    real(r8), pointer :: td_profile(:,:)         ! Canopy layer transmittance of diffuse radiation (-)

    real(r8), pointer :: swsrc_profile(:,:,:)    ! Canopy layer source/sink flux: absorbed solar radiation (W/m2) [for numrad wavebands]
    real(r8), pointer :: lwsrc_profile(:,:)      ! Canopy layer source/sink flux: absorbed longwave radiation (W/m2)
    real(r8), pointer :: rnsrc_profile(:,:)      ! Canopy layer source/sink flux: net radiation (W/m2)
    real(r8), pointer :: stsrc_profile(:,:)      ! Canopy layer source/sink flux: storage heat flux (W/m2)
    real(r8), pointer :: shsrc_profile(:,:)      ! Canopy layer source/sink flux: sensible heat (W/m2)
    real(r8), pointer :: lhsrc_profile(:,:)      ! Canopy layer source/sink flux: latent heat (W/m2)
    real(r8), pointer :: etsrc_profile(:,:)      ! Canopy layer source/sink flux: water vapor (mol H2O/m2/s)
    real(r8), pointer :: fco2src_profile(:,:)    ! Canopy layer source/sink flux: CO2 (umol CO2/m2/s)

    real(r8), pointer :: wind_profile(:,:)       ! Canopy layer wind speed (m/s)
    real(r8), pointer :: tair_profile(:,:)       ! Canopy layer air temperature (K)
    real(r8), pointer :: eair_profile(:,:)       ! Canopy layer vapor pressure (Pa)
    real(r8), pointer :: cair_profile(:,:)       ! Canopy layer atmospheric CO2 (umol/mol)
    real(r8), pointer :: tair_bef_profile(:,:)   ! Canopy layer air temperature for previous timestep (K)
    real(r8), pointer :: eair_bef_profile(:,:)   ! Canopy layer vapor pressure for previous timestep (Pa)
    real(r8), pointer :: cair_bef_profile(:,:)   ! Canopy layer atmospheric CO2 for previous timestep (umol/mol)

    real(r8), pointer :: shair_profile(:,:)      ! Canopy layer air sensible heat flux (W/m2)
    real(r8), pointer :: etair_profile(:,:)      ! Canopy layer air water vapor flux (mol H2O/m2/s)
    real(r8), pointer :: stair_profile(:,:)      ! Canopy layer air storage heat flux (W/m2)
    real(r8), pointer :: gac_profile(:,:)        ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)

    real(r8), pointer :: lwpveg_profile(:,:)     ! Canopy layer leaf water potential (MPa)
    real(r8), pointer :: lsc_profile(:,:)        ! Canopy layer leaf-specific conductance (mmol H2O/m2 leaf/s/MPa)
    real(r8), pointer :: h2ocan_profile(:,:)     ! Canopy layer intercepted water (kg H2O/m2)
    real(r8), pointer :: fwet_profile(:,:)       ! Canopy layer fraction of plant area index that is wet
    real(r8), pointer :: fdry_profile(:,:)       ! Canopy layer fraction of plant area index that is green and dry

    ! Sunlit/shaded leaf variables for canopy layers (fluxes are per m2 leaf area): dimension is (patch, level, leaf)

    real(r8), pointer :: tleaf_leaf(:,:,:)       ! Leaf temperature (K)
    real(r8), pointer :: tleaf_bef_leaf(:,:,:)   ! Leaf temperature for previous timestep (K)
    real(r8), pointer :: swleaf_leaf(:,:,:,:)    ! Leaf absorbed solar radiation (W/m2 leaf) [for numrad wavebands]
    real(r8), pointer :: lwleaf_leaf(:,:,:)      ! Leaf absorbed longwave radiation (W/m2 leaf)
    real(r8), pointer :: rnleaf_leaf(:,:,:)      ! Leaf net radiation (W/m2 leaf)
    real(r8), pointer :: stleaf_leaf(:,:,:)      ! Leaf storage heat flux (W/m2 leaf)
    real(r8), pointer :: shleaf_leaf(:,:,:)      ! Leaf sensible heat flux (W/m2 leaf)
    real(r8), pointer :: lhleaf_leaf(:,:,:)      ! Leaf latent heat flux (W/m2 leaf)
    real(r8), pointer :: trleaf_leaf(:,:,:)      ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: evleaf_leaf(:,:,:)      ! Leaf evaporation flux (mol H2O/m2 leaf/s)

    real(r8), pointer :: gbh_leaf(:,:,:)         ! Leaf boundary layer conductance: heat (mol/m2 leaf/s)
    real(r8), pointer :: gbv_leaf(:,:,:)         ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    real(r8), pointer :: gbc_leaf(:,:,:)         ! Leaf boundary layer conductance: CO2 (mol CO2/m2 leaf/s)

    real(r8), pointer :: kc_leaf(:,:,:)          ! Leaf Michaelis-Menten constant for CO2 (umol/mol)
    real(r8), pointer :: ko_leaf(:,:,:)          ! Leaf Michaelis-Menten constant for O2 (mmol/mol)
    real(r8), pointer :: cp_leaf(:,:,:)          ! Leaf CO2 compensation point (umol/mol)
    real(r8), pointer :: vcmax_leaf(:,:,:)       ! Leaf maximum carboxylation rate (umol/m2/s)
    real(r8), pointer :: jmax_leaf(:,:,:)        ! Leaf maximum electron transport rate (umol/m2/s)
    real(r8), pointer :: kp_leaf(:,:,:)          ! Leaf C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    real(r8), pointer :: ceair_leaf(:,:,:)       ! Leaf vapor pressure of air, constrained for stomatal conductance (Pa)
    real(r8), pointer :: leaf_esat_leaf(:,:,:)   ! Leaf saturation vapor pressure (Pa)

    real(r8), pointer :: apar_leaf(:,:,:)        ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    real(r8), pointer :: je_leaf(:,:,:)          ! Leaf electron transport rate (umol/m2/s)
    real(r8), pointer :: ac_leaf(:,:,:)          ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: aj_leaf(:,:,:)          ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: ap_leaf(:,:,:)          ! Leaf product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: ag_leaf(:,:,:)          ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: an_leaf(:,:,:)          ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: rd_leaf(:,:,:)          ! Leaf respiration rate (umol CO2/m2 leaf/s)
    real(r8), pointer :: ci_leaf(:,:,:)          ! Leaf intercellular CO2 (umol/mol)
    real(r8), pointer :: cs_leaf(:,:,:)          ! Leaf surface CO2 (umol/mol)

    real(r8), pointer :: lwpleaf_leaf(:,:,:)     ! Leaf water potential (MPa)
    real(r8), pointer :: hs_leaf(:,:,:)          ! Leaf fractional humidity at leaf surface (-)
    real(r8), pointer :: vpd_leaf(:,:,:)         ! Leaf vapor pressure deficit (Pa)
    real(r8), pointer :: gs_leaf(:,:,:)          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    real(r8), pointer :: alphapsn_leaf(:,:,:)    ! Leaf 13C fractionation factor for photosynthesis (-)

  contains

    procedure, public  :: Init              ! CLM initialization of data type
    procedure, private :: InitAllocate      ! CLM initialization: allocate module data structure
    procedure, private :: InitHistory       ! CLM initialization: setup history file variables
    procedure, private :: InitCold          ! CLM initialization: cold-start initialization
    procedure, public  :: Restart           ! CLM restart file

  end type mlcanopy_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)
    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start
    !
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate module data structure
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp      ! Beginning patch index for CLM g/l/c/p hierarchy
    integer :: endp      ! Ending patch index for CLM g/l/c/p hierarchy
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp

    ! Vegetation input variables

    allocate (this%ztop_canopy         (begp:endp))                              ; this%ztop_canopy         (:)       = spval
    allocate (this%lai_canopy          (begp:endp))                              ; this%lai_canopy          (:)       = spval
    allocate (this%sai_canopy          (begp:endp))                              ; this%sai_canopy          (:)       = spval
    allocate (this%root_biomass_canopy (begp:endp))                              ; this%root_biomass_canopy (:)       = spval

    ! Atmospheric forcing variables

    allocate (this%zref_forcing        (begp:endp))                              ; this%zref_forcing        (:)       = spval
    allocate (this%tref_forcing        (begp:endp))                              ; this%tref_forcing        (:)       = spval
    allocate (this%qref_forcing        (begp:endp))                              ; this%qref_forcing        (:)       = spval
    allocate (this%uref_forcing        (begp:endp))                              ; this%uref_forcing        (:)       = spval
    allocate (this%pref_forcing        (begp:endp))                              ; this%pref_forcing        (:)       = spval
    allocate (this%co2ref_forcing      (begp:endp))                              ; this%co2ref_forcing      (:)       = spval
    allocate (this%o2ref_forcing       (begp:endp))                              ; this%o2ref_forcing       (:)       = spval
    allocate (this%swskyb_forcing      (begp:endp,1:numrad))                     ; this%swskyb_forcing      (:,:)     = spval
    allocate (this%swskyd_forcing      (begp:endp,1:numrad))                     ; this%swskyd_forcing      (:,:)     = spval
    allocate (this%lwsky_forcing       (begp:endp))                              ; this%lwsky_forcing       (:)       = spval
    allocate (this%qflx_rain_forcing   (begp:endp))                              ; this%qflx_rain_forcing   (:)       = spval
    allocate (this%qflx_snow_forcing   (begp:endp))                              ; this%qflx_snow_forcing   (:)       = spval
    allocate (this%tacclim_forcing     (begp:endp))                              ; this%tacclim_forcing     (:)       = spval

    ! Derived atmospheric forcing variables

    allocate (this%eref_forcing        (begp:endp))                              ; this%eref_forcing        (:)       = spval
    allocate (this%thref_forcing       (begp:endp))                              ; this%thref_forcing       (:)       = spval
    allocate (this%thvref_forcing      (begp:endp))                              ; this%thvref_forcing      (:)       = spval
    allocate (this%rhoair_forcing      (begp:endp))                              ; this%rhoair_forcing      (:)       = spval
    allocate (this%rhomol_forcing      (begp:endp))                              ; this%rhomol_forcing      (:)       = spval
    allocate (this%mmair_forcing       (begp:endp))                              ; this%mmair_forcing       (:)       = spval
    allocate (this%cpair_forcing       (begp:endp))                              ; this%cpair_forcing       (:)       = spval
    allocate (this%solar_zen_forcing   (begp:endp))                              ; this%solar_zen_forcing   (:)       = spval

    ! Canopy flux variables

    allocate (this%swveg_canopy        (begp:endp,1:numrad))                     ; this%swveg_canopy        (:,:)     = spval
    allocate (this%swvegsun_canopy     (begp:endp,1:numrad))                     ; this%swvegsun_canopy     (:,:)     = spval
    allocate (this%swvegsha_canopy     (begp:endp,1:numrad))                     ; this%swvegsha_canopy     (:,:)     = spval
    allocate (this%lwveg_canopy        (begp:endp))                              ; this%lwveg_canopy        (:)       = spval
    allocate (this%lwvegsun_canopy     (begp:endp))                              ; this%lwvegsun_canopy     (:)       = spval
    allocate (this%lwvegsha_canopy     (begp:endp))                              ; this%lwvegsha_canopy     (:)       = spval
    allocate (this%shveg_canopy        (begp:endp))                              ; this%shveg_canopy        (:)       = spval
    allocate (this%shvegsun_canopy     (begp:endp))                              ; this%shvegsun_canopy     (:)       = spval
    allocate (this%shvegsha_canopy     (begp:endp))                              ; this%shvegsha_canopy     (:)       = spval
    allocate (this%lhveg_canopy        (begp:endp))                              ; this%lhveg_canopy        (:)       = spval
    allocate (this%lhvegsun_canopy     (begp:endp))                              ; this%lhvegsun_canopy     (:)       = spval
    allocate (this%lhvegsha_canopy     (begp:endp))                              ; this%lhvegsha_canopy     (:)       = spval
    allocate (this%etveg_canopy        (begp:endp))                              ; this%etveg_canopy        (:)       = spval
    allocate (this%etvegsun_canopy     (begp:endp))                              ; this%etvegsun_canopy     (:)       = spval
    allocate (this%etvegsha_canopy     (begp:endp))                              ; this%etvegsha_canopy     (:)       = spval
    allocate (this%gppveg_canopy       (begp:endp))                              ; this%gppveg_canopy       (:)       = spval
    allocate (this%gppvegsun_canopy    (begp:endp))                              ; this%gppvegsun_canopy    (:)       = spval
    allocate (this%gppvegsha_canopy    (begp:endp))                              ; this%gppvegsha_canopy    (:)       = spval
    allocate (this%albcan_canopy       (begp:endp,1:numrad))                     ; this%albcan_canopy       (:,:)     = spval
    allocate (this%lwup_canopy         (begp:endp))                              ; this%lwup_canopy         (:)       = spval
    allocate (this%rnet_canopy         (begp:endp))                              ; this%rnet_canopy         (:)       = spval
    allocate (this%shflx_canopy        (begp:endp))                              ; this%shflx_canopy        (:)       = spval
    allocate (this%lhflx_canopy        (begp:endp))                              ; this%lhflx_canopy        (:)       = spval
    allocate (this%etflx_canopy        (begp:endp))                              ; this%etflx_canopy        (:)       = spval
    allocate (this%stflx_canopy        (begp:endp))                              ; this%stflx_canopy        (:)       = spval
    allocate (this%ustar_canopy        (begp:endp))                              ; this%ustar_canopy        (:)       = spval
    allocate (this%gac_to_hc_canopy    (begp:endp))                              ; this%gac_to_hc_canopy    (:)       = spval
    allocate (this%qflx_intr_canopy    (begp:endp))                              ; this%qflx_intr_canopy    (:)       = spval
    allocate (this%qflx_tflrain_canopy (begp:endp))                              ; this%qflx_tflrain_canopy (:)       = spval
    allocate (this%qflx_tflsnow_canopy (begp:endp))                              ; this%qflx_tflsnow_canopy (:)       = spval

    ! Canopy diagnostic variables

    allocate (this%uaf_canopy          (begp:endp))                              ; this%uaf_canopy          (:)       = spval
    allocate (this%taf_canopy          (begp:endp))                              ; this%taf_canopy          (:)       = spval
    allocate (this%qaf_canopy          (begp:endp))                              ; this%qaf_canopy          (:)       = spval
    allocate (this%fracminlwp_canopy   (begp:endp))                              ; this%fracminlwp_canopy   (:)       = spval

    ! Canopy aerodynamic variables

    allocate (this%obu_canopy          (begp:endp))                              ; this%obu_canopy          (:)       = spval
    allocate (this%obuold_canopy       (begp:endp))                              ; this%obuold_canopy       (:)       = spval
    allocate (this%nmozsgn_canopy      (begp:endp))                              ; this%nmozsgn_canopy      (:)       = ispval
    allocate (this%beta_canopy         (begp:endp))                              ; this%beta_canopy         (:)       = spval
    allocate (this%PrSc_canopy         (begp:endp))                              ; this%PrSc_canopy         (:)       = spval
    allocate (this%Lc_canopy           (begp:endp))                              ; this%Lc_canopy           (:)       = spval
    allocate (this%zdisp_canopy        (begp:endp))                              ; this%zdisp_canopy        (:)       = spval

    ! Canopy stomatal conductance variables

    allocate (this%g0_canopy           (begp:endp))                              ; this%g0_canopy           (:)       = spval
    allocate (this%g1_canopy           (begp:endp))                              ; this%g1_canopy           (:)       = spval

    ! Soil energy balance variables

    allocate (this%albsoib_soil        (begp:endp,1:numrad))                     ; this%albsoib_soil        (:,:)     = spval
    allocate (this%albsoid_soil        (begp:endp,1:numrad))                     ; this%albsoid_soil        (:,:)     = spval
    allocate (this%swsoi_soil          (begp:endp,1:numrad))                     ; this%swsoi_soil          (:,:)     = spval
    allocate (this%lwsoi_soil          (begp:endp))                              ; this%lwsoi_soil          (:)       = spval
    allocate (this%rnsoi_soil          (begp:endp))                              ; this%rnsoi_soil          (:)       = spval
    allocate (this%shsoi_soil          (begp:endp))                              ; this%shsoi_soil          (:)       = spval
    allocate (this%lhsoi_soil          (begp:endp))                              ; this%lhsoi_soil          (:)       = spval
    allocate (this%etsoi_soil          (begp:endp))                              ; this%etsoi_soil          (:)       = spval
    allocate (this%gsoi_soil           (begp:endp))                              ; this%gsoi_soil           (:)       = spval
    allocate (this%tg_soil             (begp:endp))                              ; this%tg_soil             (:)       = spval
    allocate (this%tg_bef_soil         (begp:endp))                              ; this%tg_bef_soil         (:)       = spval
    allocate (this%eg_soil             (begp:endp))                              ; this%eg_soil             (:)       = spval
    allocate (this%rhg_soil            (begp:endp))                              ; this%rhg_soil            (:)       = spval
    allocate (this%gac0_soil           (begp:endp))                              ; this%gac0_soil           (:)       = spval
    allocate (this%soil_t_soil         (begp:endp))                              ; this%soil_t_soil         (:)       = spval
    allocate (this%soil_dz_soil        (begp:endp))                              ; this%soil_dz_soil        (:)       = spval
    allocate (this%soil_tk_soil        (begp:endp))                              ; this%soil_tk_soil        (:)       = spval
    allocate (this%soilres_soil        (begp:endp))                              ; this%soilres_soil        (:)       = spval

    ! Soil moisture variables

    allocate (this%btran_soil          (begp:endp))                              ; this%btran_soil          (:)       = spval
    allocate (this%psis_soil           (begp:endp))                              ; this%psis_soil           (:)       = spval
    allocate (this%rsoil_soil          (begp:endp))                              ; this%rsoil_soil          (:)       = spval
    allocate (this%soil_et_loss_soil   (begp:endp,1:nlevgrnd))                   ; this%soil_et_loss_soil   (:,:)     = spval

    ! Canopy layer indices

    allocate (this%ncan_canopy         (begp:endp))                              ; this%ncan_canopy         (:)       = ispval
    allocate (this%ntop_canopy         (begp:endp))                              ; this%ntop_canopy         (:)       = ispval
    allocate (this%nbot_canopy         (begp:endp))                              ; this%nbot_canopy         (:)       = ispval

    ! Canopy layer variables

    allocate (this%dlai_frac_profile   (begp:endp,1:nlevmlcan))                  ; this%dlai_frac_profile   (:,:)     = spval
    allocate (this%dsai_frac_profile   (begp:endp,1:nlevmlcan))                  ; this%dsai_frac_profile   (:,:)     = spval
    allocate (this%dlai_profile        (begp:endp,1:nlevmlcan))                  ; this%dlai_profile        (:,:)     = spval
    allocate (this%dsai_profile        (begp:endp,1:nlevmlcan))                  ; this%dsai_profile        (:,:)     = spval
    allocate (this%dpai_profile        (begp:endp,1:nlevmlcan))                  ; this%dpai_profile        (:,:)     = spval
    allocate (this%sumpai_profile      (begp:endp,1:nlevmlcan))                  ; this%sumpai_profile      (:,:)     = spval
    allocate (this%zs_profile          (begp:endp,1:nlevmlcan))                  ; this%zs_profile          (:,:)     = spval
    allocate (this%dz_profile          (begp:endp,1:nlevmlcan))                  ; this%dz_profile          (:,:)     = spval

    allocate (this%vcmax25_profile     (begp:endp,1:nlevmlcan))                  ; this%vcmax25_profile     (:,:)     = spval
    allocate (this%jmax25_profile      (begp:endp,1:nlevmlcan))                  ; this%jmax25_profile      (:,:)     = spval
    allocate (this%kp25_profile        (begp:endp,1:nlevmlcan))                  ; this%kp25_profile        (:,:)     = spval
    allocate (this%rd25_profile        (begp:endp,1:nlevmlcan))                  ; this%rd25_profile        (:,:)     = spval
    allocate (this%cpleaf_profile      (begp:endp,1:nlevmlcan))                  ; this%cpleaf_profile      (:,:)     = spval

    allocate (this%fracsun_profile     (begp:endp,1:nlevmlcan))                  ; this%fracsun_profile     (:,:)     = spval
    allocate (this%tb_profile          (begp:endp,1:nlevmlcan))                  ; this%tb_profile          (:,:)     = spval
    allocate (this%td_profile          (begp:endp,1:nlevmlcan))                  ; this%td_profile          (:,:)     = spval

    allocate (this%swsrc_profile       (begp:endp,1:nlevmlcan,1:numrad))         ; this%swsrc_profile       (:,:,:)   = spval
    allocate (this%lwsrc_profile       (begp:endp,1:nlevmlcan))                  ; this%lwsrc_profile       (:,:)     = spval
    allocate (this%rnsrc_profile       (begp:endp,1:nlevmlcan))                  ; this%rnsrc_profile       (:,:)     = spval
    allocate (this%stsrc_profile       (begp:endp,1:nlevmlcan))                  ; this%stsrc_profile       (:,:)     = spval
    allocate (this%shsrc_profile       (begp:endp,1:nlevmlcan))                  ; this%shsrc_profile       (:,:)     = spval
    allocate (this%lhsrc_profile       (begp:endp,1:nlevmlcan))                  ; this%lhsrc_profile       (:,:)     = spval
    allocate (this%etsrc_profile       (begp:endp,1:nlevmlcan))                  ; this%etsrc_profile       (:,:)     = spval
    allocate (this%fco2src_profile     (begp:endp,1:nlevmlcan))                  ; this%fco2src_profile     (:,:)     = spval

    allocate (this%wind_profile        (begp:endp,1:nlevmlcan))                  ; this%wind_profile        (:,:)     = spval
    allocate (this%tair_profile        (begp:endp,1:nlevmlcan))                  ; this%tair_profile        (:,:)     = spval
    allocate (this%eair_profile        (begp:endp,1:nlevmlcan))                  ; this%eair_profile        (:,:)     = spval
    allocate (this%cair_profile        (begp:endp,1:nlevmlcan))                  ; this%cair_profile        (:,:)     = spval
    allocate (this%tair_bef_profile    (begp:endp,1:nlevmlcan))                  ; this%tair_bef_profile    (:,:)     = spval
    allocate (this%eair_bef_profile    (begp:endp,1:nlevmlcan))                  ; this%eair_bef_profile    (:,:)     = spval
    allocate (this%cair_bef_profile    (begp:endp,1:nlevmlcan))                  ; this%cair_bef_profile    (:,:)     = spval

    allocate (this%shair_profile       (begp:endp,1:nlevmlcan))                  ; this%shair_profile       (:,:)     = spval
    allocate (this%etair_profile       (begp:endp,1:nlevmlcan))                  ; this%etair_profile       (:,:)     = spval
    allocate (this%stair_profile       (begp:endp,1:nlevmlcan))                  ; this%stair_profile       (:,:)     = spval
    allocate (this%gac_profile         (begp:endp,1:nlevmlcan))                  ; this%gac_profile         (:,:)     = spval

    allocate (this%lwpveg_profile      (begp:endp,1:nlevmlcan))                  ; this%lwpveg_profile      (:,:)     = spval
    allocate (this%lsc_profile         (begp:endp,1:nlevmlcan))                  ; this%lsc_profile         (:,:)     = spval
    allocate (this%h2ocan_profile      (begp:endp,1:nlevmlcan))                  ; this%h2ocan_profile      (:,:)     = spval
    allocate (this%fwet_profile        (begp:endp,1:nlevmlcan))                  ; this%fwet_profile        (:,:)     = spval
    allocate (this%fdry_profile        (begp:endp,1:nlevmlcan))                  ; this%fdry_profile        (:,:)     = spval

    ! Sunlit/shaded leaf variables for canopy layers

    allocate (this%tleaf_leaf          (begp:endp,1:nlevmlcan,1:nleaf))          ; this%tleaf_leaf          (:,:,:)   = spval
    allocate (this%tleaf_bef_leaf      (begp:endp,1:nlevmlcan,1:nleaf))          ; this%tleaf_bef_leaf      (:,:,:)   = spval
    allocate (this%swleaf_leaf         (begp:endp,1:nlevmlcan,1:nleaf,1:numrad)) ; this%swleaf_leaf         (:,:,:,:) = spval
    allocate (this%lwleaf_leaf         (begp:endp,1:nlevmlcan,1:nleaf))          ; this%lwleaf_leaf         (:,:,:)   = spval
    allocate (this%rnleaf_leaf         (begp:endp,1:nlevmlcan,1:nleaf))          ; this%rnleaf_leaf         (:,:,:)   = spval
    allocate (this%stleaf_leaf         (begp:endp,1:nlevmlcan,1:nleaf))          ; this%stleaf_leaf         (:,:,:)   = spval
    allocate (this%shleaf_leaf         (begp:endp,1:nlevmlcan,1:nleaf))          ; this%shleaf_leaf         (:,:,:)   = spval
    allocate (this%lhleaf_leaf         (begp:endp,1:nlevmlcan,1:nleaf))          ; this%lhleaf_leaf         (:,:,:)   = spval
    allocate (this%trleaf_leaf         (begp:endp,1:nlevmlcan,1:nleaf))          ; this%trleaf_leaf         (:,:,:)   = spval
    allocate (this%evleaf_leaf         (begp:endp,1:nlevmlcan,1:nleaf))          ; this%evleaf_leaf         (:,:,:)   = spval

    allocate (this%gbh_leaf            (begp:endp,1:nlevmlcan,1:nleaf))          ; this%gbh_leaf            (:,:,:)   = spval
    allocate (this%gbv_leaf            (begp:endp,1:nlevmlcan,1:nleaf))          ; this%gbv_leaf            (:,:,:)   = spval
    allocate (this%gbc_leaf            (begp:endp,1:nlevmlcan,1:nleaf))          ; this%gbc_leaf            (:,:,:)   = spval

    allocate (this%kc_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%kc_leaf             (:,:,:)   = spval
    allocate (this%ko_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%ko_leaf             (:,:,:)   = spval
    allocate (this%cp_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%cp_leaf             (:,:,:)   = spval
    allocate (this%vcmax_leaf          (begp:endp,1:nlevmlcan,1:nleaf))          ; this%vcmax_leaf          (:,:,:)   = spval
    allocate (this%jmax_leaf           (begp:endp,1:nlevmlcan,1:nleaf))          ; this%jmax_leaf           (:,:,:)   = spval
    allocate (this%kp_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%kp_leaf             (:,:,:)   = spval
    allocate (this%ceair_leaf          (begp:endp,1:nlevmlcan,1:nleaf))          ; this%ceair_leaf          (:,:,:)   = spval
    allocate (this%leaf_esat_leaf      (begp:endp,1:nlevmlcan,1:nleaf))          ; this%leaf_esat_leaf      (:,:,:)   = spval

    allocate (this%apar_leaf           (begp:endp,1:nlevmlcan,1:nleaf))          ; this%apar_leaf           (:,:,:)   = spval
    allocate (this%je_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%je_leaf             (:,:,:)   = spval
    allocate (this%ac_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%ac_leaf             (:,:,:)   = spval
    allocate (this%aj_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%aj_leaf             (:,:,:)   = spval
    allocate (this%ap_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%ap_leaf             (:,:,:)   = spval
    allocate (this%ag_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%ag_leaf             (:,:,:)   = spval
    allocate (this%an_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%an_leaf             (:,:,:)   = spval
    allocate (this%rd_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%rd_leaf             (:,:,:)   = spval
    allocate (this%ci_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%ci_leaf             (:,:,:)   = spval
    allocate (this%cs_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%cs_leaf             (:,:,:)   = spval

    allocate (this%lwpleaf_leaf        (begp:endp,1:nlevmlcan,1:nleaf))          ; this%lwpleaf_leaf        (:,:,:)   = spval
    allocate (this%hs_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%hs_leaf             (:,:,:)   = spval
    allocate (this%vpd_leaf            (begp:endp,1:nlevmlcan,1:nleaf))          ; this%vpd_leaf            (:,:,:)   = spval
    allocate (this%gs_leaf             (begp:endp,1:nlevmlcan,1:nleaf))          ; this%gs_leaf             (:,:,:)   = spval
    allocate (this%alphapsn_leaf       (begp:endp,1:nlevmlcan,1:nleaf))          ; this%alphapsn_leaf       (:,:,:)   = spval

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory (this, bounds)
    !
    ! !DESCRIPTION:
    ! Setup the fields that can be output on history files
    !
    ! !USES:
    use histFileMod, only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp= bounds%endp

    this%gppveg_canopy(begp:endp) = spval
    call hist_addfld1d (fname='GPP_ML', units='umol/m2s', &
         avgflag='A', long_name='Gross primary production', &
         ptr_patch=this%gppveg_canopy, set_lake=spval, set_urb=spval)

    this%lwpveg_profile(begp:endp,1:nlevmlcan) = spval
    call hist_addfld2d (fname='LWP_ML', units='MPa', type2d='nlevmlcan', &
         avgflag='A', long_name='Leaf water potential of canopy layer', &
         ptr_patch=this%lwpveg_profile, set_lake=spval, set_urb=spval)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold (this, bounds)
    !
    ! !DESCRIPTION:
    ! Cold-start initialization for multilayer canopy
    !
    ! !USES:
    use MLclm_varpar, only : nlevmlcan
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: p                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic               ! Canopy leaf layer index
    !---------------------------------------------------------------------

    !!! ---------------------------------------------- !!!
    !!! WARNING: this initialization is overwritten by !!!
    !!! subroutine initVerticalProfiles                !!!
    !!! ---------------------------------------------- !!!

    ! Initialize leaf water potential and intercepted water

    do p = bounds%begp, bounds%endp
       do ic = 1, nlevmlcan
          this%lwpveg_profile(p,ic) = -0.1_r8
          this%h2ocan_profile(p,ic) = 0._r8
       end do
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart (this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod, only : restartvar
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !---------------------------------------------------------------------

    ! Example for 1-d patch variable

    call restartvar(ncid=ncid, flag=flag, varname='taf_ml', xtype=ncd_double,  &
       dim1name='pft', long_name='air temperature at canopy top', units='K', &
       interpinic_flag='interp', readvar=readvar, data=this%taf_canopy)

    ! Example for 2-d patch variable

    call restartvar(ncid=ncid, flag=flag, varname='lwp_ml', xtype=ncd_double,  &
       dim1name='pft', dim2name='nlevmlcan', switchdim=.true., &
       long_name='leaf water potential of canopy layer', units='MPa', &
       interpinic_flag='interp', readvar=readvar, data=this%lwpveg_profile)

  end subroutine Restart

end module MLCanopyFluxesType
