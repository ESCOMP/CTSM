module FATESFireBase

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Abstract base class for FATES fire data object
  !
  ! !USES:
  use CNFireBaseMod                      , only : cnfire_base_type
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use CNVegStateType                     , only : cnveg_state_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type

  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: fates_fire_base_type
  !
  type, abstract, extends(cnfire_base_type) :: fates_fire_base_type
      private

      ! !PRIVATE MEMBER DATA:

    contains
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure(GetLight24_interface),    public, deferred :: GetLight24     ! Return the 24-hour averaged lightning data
      procedure(GetGDP_interface),        public, deferred :: GetGDP         ! Return the global gdp data
      procedure(InitAccBuffer_interface), public, deferred :: InitAccBuffer  ! Initialize accumulation processes
      procedure(InitAccVars_interface),   public, deferred :: InitAccVars    ! Initialize accumulation variables
      procedure(UpdateAccVars_interface), public, deferred :: UpdateAccVars  ! Update/extract accumulations vars
      ! Interfaces that need to be implemented because they are in the base class
      ! They are NOT used when FATES is on
      procedure, public :: CNFireReadParams                                  ! Read in parameters    (NOT USED FOR FATES)
      procedure, public :: CNFireArea                                        ! Calculate fire area   (NOT USED FOR FATES)
      procedure, public :: CNFireFluxes                                      ! Calculate fire fluxes (NOT USED FOR FATES)
      
  end type fates_fire_base_type

  !-----------------------

  abstract interface
  !-----------------------------------------------------------------------

  !------------------------------------------------------------------------
  function GetLight24_interface( this ) result(lnfm24)
    !
    ! !DESCRIPTION: Get the 24-hour averaged lightning data
    ! !USES
    use shr_kind_mod   , only: r8 => shr_kind_r8
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    real(r8), pointer :: lnfm24(:)
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
  end function GetLight24_interface
  
  !------------------------------------------------------------------------
  function GetGDP_interface( this ) result(gdp)
    !
    ! !DESCRIPTION: Get the global gross domestic product data
    ! !USES
    use shr_kind_mod   , only: r8 => shr_kind_r8
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    real(r8), pointer :: gdp(:)
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
  end function GetGDP_interface

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer_interface (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize the accumulation buffers
    !
    ! !USES
    use decompMod      , only: bounds_type
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

  end subroutine InitAccBuffer_interface

  !-----------------------------------------------------------------------
  subroutine InitAccVars_interface(this, bounds)
    !
    ! !DESCRIPTION:
    !      Initialize the accumulation variables
    !
    ! !USES
    use decompMod      , only: bounds_type
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

  end subroutine InitAccVars_interface

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars_interface (this, bounds)
    !
    ! !DESCRIPTION:
    ! Update accumulation variables
    !
    ! !USES
    use decompMod      , only: bounds_type
    import :: fates_fire_base_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

  end subroutine UpdateAccVars_interface

  end interface

  contains

  !-----------------------------------------------------------------------
  ! Implement empty subroutines that are required in the FireMethodType base
  ! class, but are NOT used in the FATES version
  !-----------------------------------------------------------------------
  subroutine CNFireFluxes (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      num_actfirec, filter_actfirec, num_actfirep, filter_actfirep, &
      dgvs_inst, cnveg_state_inst,                                                                      &
      cnveg_carbonstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
      soilbiogeochem_carbonflux_inst,                                       &
      leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch, &
      totsomc_col, decomp_cpools_vr_col, decomp_npools_vr_col, somc_fire_col)
   ! !DESCRIPTION:
   ! Fire effects routine for coupled carbon-nitrogen code (CN).  (NOT USED FOR FATES)
   !
   ! !USES:
   use CNDVType                           , only : dgvs_type
   use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
   use CNVegNitrogenStateType             , only : cnveg_nitrogenstate_type
   use CNVegNitrogenFluxType              , only : cnveg_nitrogenflux_type
   !
   ! !ARGUMENTS:
   class(fates_fire_base_type)                    :: this
   type(bounds_type)              , intent(in)    :: bounds  
   integer                        , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
   integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
   integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
   integer                        , intent(out)   :: num_actfirep    ! number of active patches on fire in filter
   integer                        , intent(out)   :: filter_actfirep(:) ! filter for soil patches
   integer                        , intent(out)   :: num_actfirec    ! number of active columns on fire in filter
   integer                        , intent(out)   :: filter_actfirec(:) ! filter for soil columns
   type(dgvs_type)                , intent(inout) :: dgvs_inst
   type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
   type(soilbiogeochem_carbonflux_type), intent(inout) :: soilbiogeochem_carbonflux_inst  ! only for matrix_decomp_fire_k: (gC/m3/step) VR deomp. C fire loss in matrix representation
   type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
   type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
   type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
   real(r8)                       , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: totsomc_col(bounds%begc:)
   real(r8)                       , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:)
   real(r8)                       , intent(in)    :: decomp_npools_vr_col(bounds%begc:,1:,1:)
   real(r8)                       , intent(out)   :: somc_fire_col(bounds%begc:)
   !
   call endrun( "This subroutine should NEVER be called when FATES is active" )
  end subroutine CNFireFluxes

  !-----------------------------------------------------------------------
  subroutine CNFireReadParams( this, ncid )
    !
    ! Read in the constant parameters from the input NetCDF parameter file  (NOT USED FOR FATES)
    ! !USES:
    use ncdio_pio        , only : file_desc_t
    !
    ! !ARGUMENTS:
    implicit none
    class(fates_fire_base_type)     :: this
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    call endrun( "This subroutine should NEVER be called when FATES is active" )
  end subroutine CNFireReadParams

  !-----------------------------------------------------------------------
  subroutine CNFireArea (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       num_exposedvegp, filter_exposedvegp, num_noexposedvegp, filter_noexposedvegp, &
       atm2lnd_inst, energyflux_inst, saturated_excess_runoff_inst, &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, &
       waterstatebulk_inst, soilstate_inst, soil_water_retention_curve, &
       crop_inst, cnveg_state_inst, cnveg_carbonstate_inst, totlitc_col, decomp_cpools_vr_col, t_soi17cm_col)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area  (NOT USED FOR FATES)
    !
    ! !USES:
    use EnergyFluxType                     , only : energyflux_type
    use SaturatedExcessRunoffMod           , only : saturated_excess_runoff_type
    use WaterDiagnosticBulkType            , only : waterdiagnosticbulk_type
    use Wateratm2lndBulkType               , only : wateratm2lndbulk_type
    use WaterStateBulkType                 , only : waterstatebulk_type
    use SoilStateType                      , only : soilstate_type
    use SoilWaterRetentionCurveMod         , only : soil_water_retention_curve_type
    use atm2lndType                        , only : atm2lnd_type
    use CropType                           , only: crop_type
    !
    ! !ARGUMENTS:
    class(fates_fire_base_type)                           :: this
    type(bounds_type)                     , intent(in)    :: bounds 
    integer                               , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                               , intent(in)    :: num_exposedvegp        ! number of points in filter_exposedvegp
    integer                               , intent(in)    :: filter_exposedvegp(:)  ! patch filter for non-snow-covered veg
    integer                               , intent(in)    :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer                               , intent(in)    :: filter_noexposedvegp(:) ! patch filter where frac_veg_nosno is 0
    type(atm2lnd_type)                    , intent(in)    :: atm2lnd_inst
    type(energyflux_type)                 , intent(in)    :: energyflux_inst
    type(saturated_excess_runoff_type)    , intent(in)    :: saturated_excess_runoff_inst
    type(waterdiagnosticbulk_type)        , intent(in)    :: waterdiagnosticbulk_inst
    type(wateratm2lndbulk_type)           , intent(in)    :: wateratm2lndbulk_inst
    type(waterstatebulk_type)             , intent(in)    :: waterstatebulk_inst
    type(soilstate_type)                  , intent(in)    :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in)    :: soil_water_retention_curve
    type(cnveg_state_type)                , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)          , intent(inout) :: cnveg_carbonstate_inst
    type(crop_type)                       , intent(in) :: crop_inst

    real(r8)                              , intent(in)    :: totlitc_col(bounds%begc:)
    real(r8)                              , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:)
    real(r8)                              , intent(in)    :: t_soi17cm_col(bounds%begc:)
    !
    call endrun( "This subroutine should NEVER be called when FATES is active" )
 end subroutine CNFireArea
 !----------------------------------------------

end module FATESFireBase
