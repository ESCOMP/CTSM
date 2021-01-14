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
  use atm2lndType                        , only : atm2lnd_type
  use EnergyFluxType                     , only : energyflux_type
  use SaturatedExcessRunoffMod           , only : saturated_excess_runoff_type
  use WaterDiagnosticBulkType            , only : waterdiagnosticbulk_type
  use Wateratm2lndBulkType               , only : wateratm2lndbulk_type
  use WaterStateBulkType                 , only : waterstatebulk_type
  use WaterStateBulkType                 , only : waterstatebulk_type
  use SoilStateType                      , only : soilstate_type
  use SoilWaterRetentionCurveMod         , only : soil_water_retention_curve_type
  use CNVegStateType                     , only : cnveg_state_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type

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
      procedure(InitAccBuffer_interface), public, deferred :: InitAccBuffer  ! Initialize accumulation processes
      procedure(InitAccVars_interface),   public, deferred :: InitAccVars    ! Initialize accumulation variables
      procedure(UpdateAccVars_interface), public, deferred :: UpdateAccVars  ! Update/extract accumulations vars
      
      procedure, public :: CNFireArea

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

    subroutine CNFireArea (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      num_exposedvegp, filter_exposedvegp, num_noexposedvegp, filter_noexposedvegp, &
      atm2lnd_inst, energyflux_inst, saturated_excess_runoff_inst, waterdiagnosticbulk_inst, &
      wateratm2lndbulk_inst, waterstatebulk_inst, soilstate_inst, soil_water_retention_curve, &
      cnveg_state_inst, cnveg_carbonstate_inst, totlitc_col, decomp_cpools_vr_col, t_soi17cm_col)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area 
    !
    ! !USES:
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
      real(r8)                              , intent(in)    :: totlitc_col(bounds%begc:)
      real(r8)                              , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:)
      real(r8)                              , intent(in)    :: t_soi17cm_col(bounds%begc:)

      call endrun( 'cnfire_base::CNFireArea: this method MUST be implemented!' )

    end subroutine CNFireArea

----------------------------------------------

end module FATESFireBase
