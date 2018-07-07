module CNFireMethodMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Abstract base class for functions to implement CN and BGC fire model
  !
  ! Created by Erik Kluzek, following Bill Sack's implementation of polymorphism
  ! !USES:
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: cnfire_method_type

  type, abstract :: cnfire_method_type
     logical, public :: need_lightning_and_popdens = .true.   ! If need lightning and population density 
   contains

     ! Initialize the fire datasets
     procedure(CNFireInit_interface)   , public, deferred :: CNFireInit

     ! Read namelist for the fire datasets
     procedure(CNFireReadNML_interface), public, deferred :: CNFireReadNML

     ! Interpolate the fire datasets
     procedure(CNFireInterp_interface) , public, deferred :: CNFireInterp

     ! Figure out the fire area
     procedure(CNFireArea_interface)   , public, deferred :: CNFireArea

     ! Figure out the fire fluxes
     procedure(CNFireFluxes_interface) , public, deferred :: CNFireFluxes

  end type cnfire_method_type

  abstract interface

     ! Note: The following code is adapted based on what Bill Sacks has done for soil water retention curve
     ! polymorphism. Therefore, I also keep some suggestions he gave there.
     !
     ! - Make the interfaces contain all possible inputs that are needed by any
     !   implementation; each implementation will then ignore the inputs it doesn't need.
     !
     ! - For inputs that are needed only by particular implementations - and particularly
     !   for inputs that are constant in time 
     !   pass these into the constructor, and save pointers to these inputs as components
     !   of the child type that needs them. Then they aren't needed as inputs to the
     !   individual routines, allowing the interfaces for these routines to remain more
     !   consistent between different implementations.
     !
     !---------------------------------------------------------------------------
  subroutine CNFireInit_interface(this, bounds, NLFilename )
    !
    ! !DESCRIPTION:
    ! Initialize CN Fire datasets
    !
    ! USES
    use decompMod              , only : bounds_type
    import :: cnfire_method_type
    ! !ARGUMENTS:
    class(cnfire_method_type)     :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*),  intent(in) :: NLFilename
    !-----------------------------------------------------------------------

  end subroutine CNFireInit_interface

  subroutine CNFireReadNML_interface(this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read general fire namelist
    !
    ! USES
    import :: cnfire_method_type
    ! !ARGUMENTS:
    class(cnfire_method_type)     :: this
    character(len=*),  intent(in) :: NLFilename
    !-----------------------------------------------------------------------

  end subroutine CNFireReadNML_interface

  subroutine CNFireInterp_interface(this, bounds)
    !
    ! !DESCRIPTION:
    ! Interpolate CN Fire datasets
    !
    ! USES
    use decompMod              , only : bounds_type
    import :: cnfire_method_type
    ! !ARGUMENTS:
    class(cnfire_method_type)     :: this
    type(bounds_type), intent(in) :: bounds
    !-----------------------------------------------------------------------

  end subroutine CNFireInterp_interface

  !-----------------------------------------------------------------------
  subroutine CNFireArea_interface (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       atm2lnd_inst, energyflux_inst, saturated_excess_runoff_inst, waterdiagnosticbulk_inst, &
       cnveg_state_inst, cnveg_carbonstate_inst, totlitc_col, decomp_cpools_vr_col, t_soi17cm_col)
    !
    ! !DESCRIPTION:
    ! Computes column-level burned area 
    !
    ! !USES:
    use shr_kind_mod                       , only : r8 => shr_kind_r8
    use decompMod                          , only : bounds_type
    use atm2lndType                        , only : atm2lnd_type
    use EnergyFluxType                     , only : energyflux_type
    use SaturatedExcessRunoffMod           , only : saturated_excess_runoff_type
    use WaterDiagnosticBulkType                     , only : waterdiagnosticbulk_type
    use CNVegStateType                     , only : cnveg_state_type
    use CNVegCarbonStateType               , only : cnveg_carbonstate_type
    import :: cnfire_method_type
    !
    ! !ARGUMENTS:
    class(cnfire_method_type)                             :: this
    type(bounds_type)                     , intent(in)    :: bounds
    integer                               , intent(in)    :: num_soilc       !  number of soil columns in filter
    integer                               , intent(in)    :: filter_soilc(:) !  filter for soil columns
    integer                               , intent(in)    :: num_soilp       !  number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) !  filter for soil patches
    type(atm2lnd_type)                    , intent(in)    :: atm2lnd_inst
    type(energyflux_type)                 , intent(in)    :: energyflux_inst
    type(saturated_excess_runoff_type)    , intent(in)    :: saturated_excess_runoff_inst
    type(waterdiagnosticbulk_type)                 , intent(in)    :: waterdiagnosticbulk_inst
    type(cnveg_state_type)                , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)          , intent(inout) :: cnveg_carbonstate_inst
    real(r8)                              , intent(in)    :: totlitc_col(bounds%begc:)
    real(r8)                              , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:)
    real(r8)                              , intent(in)    :: t_soi17cm_col(bounds%begc:)
    !-----------------------------------------------------------------------
  end subroutine CNFireArea_interface

 !-----------------------------------------------------------------------
 subroutine CNFireFluxes_interface (this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      dgvs_inst, cnveg_state_inst, &
      cnveg_carbonstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
      leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch, &
      totsomc_col, decomp_cpools_vr_col, decomp_npools_vr_col, somc_fire_col)
   !
   ! !DESCRIPTION:
   ! Fire effects routine for coupled carbon-nitrogen code (CN).
   ! Relies primarily on estimate of fractional area burned, from CNFireArea().
   !
   ! Total fire carbon emissions (g C/m2 land area/yr) 
   !  =avg(COL_FIRE_CLOSS)*seconds_per_year + avg(SOMC_FIRE)*seconds_per_year + 
   !   avg(LF_CONV_CFLUX)*seconds_per_year*min(1.0,avg(LFC2)*seconds_per_year)*0.8
   ! where avg means the temporal average in a year
   ! seconds_per_year is the number of seconds in a year.
   !
   ! !USES:
   use shr_kind_mod                       , only : r8 => shr_kind_r8
   use decompMod                          , only : bounds_type
   use CNDVType                           , only : dgvs_type
   use CNVegStateType                     , only : cnveg_state_type
   use CNVegCarbonStateType               , only : cnveg_carbonstate_type
   use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
   use CNVegNitrogenStateType             , only : cnveg_nitrogenstate_type
   use CNVegNitrogenFluxType              , only : cnveg_nitrogenflux_type
   import :: cnfire_method_type
   !
   ! !ARGUMENTS:
   class(cnfire_method_type)                      :: this
   type(bounds_type)              , intent(in)    :: bounds
   integer                        , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
   integer                        , intent(in)    :: num_soilp       ! number of soil patches in filter
   integer                        , intent(in)    :: filter_soilp(:) ! filter for soil patches
   type(dgvs_type)                , intent(inout) :: dgvs_inst
   type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
   type(cnveg_carbonstate_type)   , intent(inout) :: cnveg_carbonstate_inst
   type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
   type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
   real(r8)                       , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
   real(r8)                       , intent(in)    :: totsomc_col(bounds%begc:) ! (gC/m2) total soil organic matter C
   real(r8)                       , intent(in)    :: decomp_cpools_vr_col(bounds%begc:,1:,1:) ! (gC/m3)  VR decomp. (litter, cwd, soil)
   real(r8)                       , intent(in)    :: decomp_npools_vr_col(bounds%begc:,1:,1:) ! (gC/m3)  VR decomp. (litter, cwd, soil)
   real(r8)                       , intent(out)   :: somc_fire_col(bounds%begc:) ! (gC/m2/s) fire C emissions due to peat burning
   !-----------------------------------------------------------------------
  end subroutine CNFireFluxes_interface

 !-----------------------------------------------------------------------

  end interface

end module CNFireMethodMod
