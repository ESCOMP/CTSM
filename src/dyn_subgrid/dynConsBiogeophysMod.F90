module dynConsBiogeophysMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeophysical quantities (water & energy) with dynamic land
  ! cover.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use decompMod               , only : bounds_type
  use UrbanParamsType         , only : urbanparams_type
  use EnergyFluxType          , only : energyflux_type
  use SoilHydrologyType       , only : soilhydrology_type
  use SoilStateType           , only : soilstate_type
  use TemperatureType         , only : temperature_type
  use WaterFluxType           , only : waterflux_type
  use WaterStateBulkType      , only : waterstatebulk_type
  use WaterStateType          , only : waterstate_type
  use WaterDiagnosticType     , only : waterdiagnostic_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use WaterBalanceType        , only : waterbalance_type
  use WaterType               , only : water_type
  use TotalWaterAndHeatMod    , only : AccumulateSoilLiqIceMassNonLake
  use TotalWaterAndHeatMod    , only : AccumulateSoilHeatNonLake
  use TotalWaterAndHeatMod    , only : ComputeLiqIceMassNonLake, ComputeLiqIceMassLake
  use TotalWaterAndHeatMod    , only : ComputeHeatNonLake, ComputeHeatLake
  use TotalWaterAndHeatMod    , only : AdjustDeltaHeatForDeltaLiq
  use TotalWaterAndHeatMod    , only : heat_base_temp
  use subgridAveMod           , only : p2c, c2g, c2l
  use GridcellType            , only : grc
  use LandunitType            , only : lun
  use ColumnType              , only : col
  use PatchType               , only : patch
  use clm_varcon              , only : tfrz, cpliq, hfus, ispval
  use landunit_varcon         , only : istsoil, istice_mec
  use dynSubgridControlMod    , only : get_for_testing_zero_dynbal_fluxes
  use filterColMod            , only : filter_col_type, col_filter_from_ltypes
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dyn_hwcontent_set_baselines ! set start-of-run baseline values for heat and water content in some columns
  public :: dyn_hwcontent_init          ! compute grid-level heat and water content, before land cover change
  public :: dyn_hwcontent_final         ! compute grid-level heat and water content and dynbal fluxes after land cover change
  !
  ! !PRIVATE MEMBER FUNCTIONS
  private :: dyn_water_content_set_baselines ! set start-of-run baseline values for water content, for a single water tracer or bulk water
  private :: dyn_heat_content_set_baselines  ! set start-of-run baseline values for heat content
  private :: set_glacier_baselines           ! compute start-of-run baseline values for each glacier column, for some dynbal term
  private :: dyn_water_content_final      ! compute grid-level water content and dynbal fluxes after landcover change, for a single water tracer or bulk water
  private :: dyn_water_content            ! compute gridcell total liquid and ice water contents
  private :: dyn_heat_content             ! compute gridcell total heat contents

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dyn_hwcontent_set_baselines(bounds, num_icemecc, filter_icemecc, &
       urbanparams_inst, soilstate_inst, water_inst, temperature_inst)
    !
    ! !DESCRIPTION:
    ! Set start-of-run baseline values for heat and water content in some columns.
    !
    ! These baseline values will be subtracted from each column's total heat and water
    ! calculations in each time step. This can correct for things like virtual mass (e.g.,
    ! the virtual ice column in glacier columns), or can add states that are not
    ! explicitly represented in the model (e.g., adding typical values for soil heat and
    ! water content in glacier columns, where the soil column is not explicitly
    ! represented). These corrections will typically reduce the fictitious dynbal
    ! conservation fluxes.
    !
    ! At a minimum, this should be called during cold start initialization, to initialize
    ! the baseline values based on cold start states. In addition, this can be called when
    ! starting a new run from existing initial conditions, if the user desires. This
    ! optional resetting can further reduce the dynbal fluxes; however, it can break
    ! conservation. (So, for example, it can be done when transitioning from an offline
    ! spinup to a coupled run, but it should not be done when transitioning from a
    ! coupled historical run to a future scenario.)
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds

    ! The following filter should include inactive as well as active points. This could
    ! be important if an inactive point later becomes active, so that we have an
    ! appropriate baseline value for that point.
    integer, intent(in) :: num_icemecc  ! number of points in filter_icemecc
    integer, intent(in) :: filter_icemecc(:) ! filter for icemec (i.e., glacier) columns

    type(urbanparams_type), intent(in) :: urbanparams_inst
    type(soilstate_type), intent(in) :: soilstate_inst
    type(water_type), intent(inout) :: water_inst
    type(temperature_type), intent(inout) :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer :: i
    type(filter_col_type) :: natveg_and_glc_filterc

    character(len=*), parameter :: subname = 'dyn_hwcontent_set_baselines'
    !-----------------------------------------------------------------------

    ! Note that we include inactive points in the following. This could be important if
    ! an inactive point later becomes active, so that we have an appropriate baseline
    ! value for that point.
    natveg_and_glc_filterc = col_filter_from_ltypes( &
         bounds = bounds, &
         ltypes = [istsoil, istice_mec], &
         include_inactive = .true.)

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(bulk_or_tracer => water_inst%bulk_and_tracers(i))

       call dyn_water_content_set_baselines(bounds, natveg_and_glc_filterc, &
            num_icemecc, filter_icemecc, &
            bulk_or_tracer%waterstate_inst)
       end associate
    end do

    call dyn_heat_content_set_baselines(bounds, natveg_and_glc_filterc, &
         num_icemecc, filter_icemecc, &
         urbanparams_inst, soilstate_inst, water_inst%waterstatebulk_inst, &
         temperature_inst)

  end subroutine dyn_hwcontent_set_baselines

  !-----------------------------------------------------------------------
  subroutine dyn_water_content_set_baselines(bounds, natveg_and_glc_filterc, &
       num_icemecc, filter_icemecc, &
       waterstate_inst)
    !
    ! !DESCRIPTION:
    ! Set start-of-run baseline values for water content, for a single water tracer or
    ! bulk water.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(filter_col_type), intent(in) :: natveg_and_glc_filterc  ! filter for natural veg and glacier columns

    ! The following filter should include inactive as well as active points
    integer, intent(in) :: num_icemecc  ! number of points in filter_icemecc
    integer, intent(in) :: filter_icemecc(:) ! filter for icemec (i.e., glacier) columns

    class(waterstate_type), intent(inout) :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: soil_liquid_mass_col(bounds%begc:bounds%endc)
    real(r8) :: soil_ice_mass_col(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'dyn_water_content_set_baselines'
    !-----------------------------------------------------------------------

    associate( &
         dynbal_baseline_liq => waterstate_inst%dynbal_baseline_liq_col, & ! Output: [real(r8) (:)   ]  baseline liquid water content subtracted from each column's total liquid water calculation (mm H2O)
         dynbal_baseline_ice => waterstate_inst%dynbal_baseline_ice_col  & ! Output: [real(r8) (:)   ]  baseline ice content subtracted from each column's total ice calculation (mm H2O)
         )

    soil_liquid_mass_col(bounds%begc:bounds%endc) = 0._r8
    soil_ice_mass_col   (bounds%begc:bounds%endc) = 0._r8

    call AccumulateSoilLiqIceMassNonLake(bounds, &
         natveg_and_glc_filterc%num, natveg_and_glc_filterc%indices, &
         waterstate_inst, &
         liquid_mass = soil_liquid_mass_col(bounds%begc:bounds%endc), &
         ice_mass = soil_ice_mass_col(bounds%begc:bounds%endc))

    ! For glacier columns, the liquid and ice in the "soil" (i.e., in the glacial ice) is
    ! a virtual pool. So we'll subtract this amount when determining the dynbal
    ! fluxes. (Note that a positive value in these baseline variables indicates an extra
    ! stock that will be subtracted later.) But glacier columns do not represent the soil
    ! under the glacial ice. Let's assume that the soil state under each glacier column is
    ! the same as the soil state in the natural vegetation landunit on that grid cell. We
    ! subtract this from the dynbal baseline variables to indicate a missing stock.
    call set_glacier_baselines(bounds, num_icemecc, filter_icemecc, &
         vals_col = soil_liquid_mass_col(bounds%begc:bounds%endc), &
         baselines_col = dynbal_baseline_liq(bounds%begc:bounds%endc))
    call set_glacier_baselines(bounds, num_icemecc, filter_icemecc, &
         vals_col = soil_ice_mass_col(bounds%begc:bounds%endc), &
         baselines_col = dynbal_baseline_ice(bounds%begc:bounds%endc))

    end associate

  end subroutine dyn_water_content_set_baselines

  !-----------------------------------------------------------------------
  subroutine dyn_heat_content_set_baselines(bounds, natveg_and_glc_filterc, &
       num_icemecc, filter_icemecc, &
       urbanparams_inst, soilstate_inst, waterstatebulk_inst, &
       temperature_inst)
    !
    ! !DESCRIPTION:
    ! Set start-of-run baseline values for heat content.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(filter_col_type), intent(in) :: natveg_and_glc_filterc  ! filter for natural veg and glacier columns

    ! The following filter should include inactive as well as active points
    integer, intent(in) :: num_icemecc  ! number of points in filter_icemecc
    integer, intent(in) :: filter_icemecc(:) ! filter for icemec (i.e., glacier) columns

    type(urbanparams_type), intent(in) :: urbanparams_inst
    type(soilstate_type), intent(in) :: soilstate_inst
    type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
    type(temperature_type), intent(inout) :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: soil_heat_col(bounds%begc:bounds%endc) ! soil heat content [J/m^2]
    real(r8) :: soil_heat_liquid_col(bounds%begc:bounds%endc)  ! unused; just needed for AccumulateSoilHeatNonLake interface
    real(r8) :: soil_cv_liquid_col(bounds%begc:bounds%endc)  ! unused; just needed for AccumulateSoilHeatNonLake interface

    character(len=*), parameter :: subname = 'dyn_heat_content_set_baselines'
    !-----------------------------------------------------------------------

    associate( &
         dynbal_baseline_heat => temperature_inst%dynbal_baseline_heat_col & ! Output: [real(r8) (:)   ]  baseline heat content subtracted from each column's total heat calculation (J/m2)
         )

    soil_heat_col(bounds%begc:bounds%endc) = 0._r8
    soil_heat_liquid_col(bounds%begc:bounds%endc) = 0._r8
    soil_cv_liquid_col(bounds%begc:bounds%endc) = 0._r8

    call AccumulateSoilHeatNonLake(bounds, &
         natveg_and_glc_filterc%num, natveg_and_glc_filterc%indices, &
         urbanparams_inst, soilstate_inst, temperature_inst, waterstatebulk_inst, &
         heat = soil_heat_col(bounds%begc:bounds%endc), &
         heat_liquid = soil_heat_liquid_col(bounds%begc:bounds%endc), &
         cv_liquid = soil_cv_liquid_col(bounds%begc:bounds%endc))

    ! See comments in dyn_water_content_set_baselines for rationale for these glacier
    ! baselines. Even though the heat in glacier ice can interact with the rest of the
    ! system (e.g., giving off heat to the atmosphere or receiving heat from the
    ! atmosphere), it is still a virtual state in the sense that the glacier ice column
    ! is virtual. And, as for water, we subtract the heat of the soil in the associated
    ! natural vegetated landunit to account for the fact that we don't explicitly model
    ! the soil under glacial ice.
    !
    ! Aside from these considerations of what is virtual vs. real, another rationale
    ! driving the use of these baselines is the desire to minimize the dynbal fluxes. For
    ! the sake of conservation, it seems okay to pick any fixed baseline when summing the
    ! energy (or water) content of a given column (as long as that baseline doesn't
    ! change over time). By using the baselines computed here, we reduce the dynbal
    ! fluxes to more reasonable values.
    call set_glacier_baselines(bounds, num_icemecc, filter_icemecc, &
         vals_col = soil_heat_col(bounds%begc:bounds%endc), &
         baselines_col = dynbal_baseline_heat(bounds%begc:bounds%endc))

    end associate

  end subroutine dyn_heat_content_set_baselines

  !-----------------------------------------------------------------------
  subroutine set_glacier_baselines(bounds, num_icemecc, filter_icemecc, &
       vals_col, baselines_col)
    !
    ! !DESCRIPTION:
    ! Compute start-of-run baseline values for each glacier column, for some dynbal term.
    !
    ! The baselines for a given glacier column are set equal to the input value in that
    ! glacier column minus the average value for the vegetated landunit on that column's
    ! gridcell.
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds

    ! The following filter should include inactive as well as active points
    integer, intent(in) :: num_icemecc  ! number of points in filter_icemecc
    integer, intent(in) :: filter_icemecc(:) ! filter for icemec (i.e., glacier) columns

    real(r8), intent(in) :: vals_col( bounds%begc: ) ! values in each column; must be set for at least natural veg and glacier columns
    real(r8), intent(inout) :: baselines_col( bounds%begc: ) ! baselines in each column; will be set for all points in the icemecc filter
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, l, g  ! indices
    integer  :: l_natveg     ! index of natural veg landunit on this grid cell
    real(r8) :: vals_lun(bounds%begl:bounds%endl)

    character(len=*), parameter :: subname = 'set_glacier_baselines'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(vals_col) == [bounds%endc]), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(baselines_col) == [bounds%endc]), errMsg(sourcefile, __LINE__))

    ! Average values from column to landunit for the natural veg landunit, in case there
    ! are multiple columns on the natural veg landunit.
    !
    ! It's somewhat inefficient that we produce landunit averages for all landunits, when
    ! all we really need are averages for the natural veg landunit. But since this
    ! subroutine is only called in initialization, code reuse/simplicity is more important
    ! than performance.
    !
    ! No thought has been given to c2l_scale_type here. This may need to be fixed if we
    ! ever start caring about urban averages in this routine.
    call c2l(bounds = bounds, &
         carr = vals_col(bounds%begc:bounds%endc), &
         larr = vals_lun(bounds%begl:bounds%endl), &
         c2l_scale_type = 'urbanf', &
         include_inactive = .true.)

    do fc = 1, num_icemecc
       c = filter_icemecc(fc)
       g = col%gridcell(c)

       ! Start by setting the baseline for this glacier column equal to the value in this
       ! glacier column.
       baselines_col(c) = vals_col(c)

       ! Now subtract the value in the natural vegetated landunit on this gridcell. This
       ! subtraction represents the missing stock, since glacier columns do not have
       ! explicit stocks for the soil under the glacial ice. So we are assuming here that
       ! the soil state under each glacier column is the same as the soil state in the
       ! natural vegetation landunit on that grid cell.
       !
       ! Note that we don't do this subtraction if the gridcell doesn't have a natural
       ! vegetation landunit. But we shouldn't have any transient glacier areas in that
       ! case anyway, so that shouldn't matter.
       l_natveg = grc%landunit_indices(istsoil, g)
       if (l_natveg /= ispval) then
          baselines_col(c) = baselines_col(c) - vals_lun(l_natveg)
       end if
    end do

  end subroutine set_glacier_baselines


  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_init(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       urbanparams_inst, soilstate_inst, &
       water_inst, temperature_inst)
    !
    ! !DESCRIPTION:
    ! Compute grid cell-level heat and water content before land cover change
    !
    ! Should be called BEFORE any subgrid weight updates this time step
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_nolakec
    integer                  , intent(in)    :: filter_nolakec(:)
    integer                  , intent(in)    :: num_lakec
    integer                  , intent(in)    :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)    :: urbanparams_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(water_type)         , intent(inout) :: water_inst
    type(temperature_type)   , intent(inout) :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer :: i

    !-------------------------------------------------------------------------------

    associate( &
         begg => bounds%begg, &
         endg => bounds%endg  &
         )

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(bulk_or_tracer => water_inst%bulk_and_tracers(i))
       call dyn_water_content(bounds, &
            num_nolakec, filter_nolakec, &
            num_lakec, filter_lakec, &
            bulk_or_tracer%waterstate_inst, &
            bulk_or_tracer%waterdiagnostic_inst, &
            liquid_mass = bulk_or_tracer%waterbalance_inst%liq1_grc(begg:endg), &
            ice_mass = bulk_or_tracer%waterbalance_inst%ice1_grc(begg:endg))
       end associate
    end do

    call dyn_heat_content( bounds, &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         urbanparams_inst, soilstate_inst, &
         temperature_inst, &
         water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
         heat_grc = temperature_inst%heat1_grc(begg:endg), &
         liquid_water_temp_grc = temperature_inst%liquid_water_temp1_grc(begg:endg))

    end associate

  end subroutine dyn_hwcontent_init

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_final(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       urbanparams_inst, soilstate_inst, &
       water_inst, &
       temperature_inst, energyflux_inst)
    !
    ! !DESCRIPTION:
    ! Compute grid cell-level heat and water content and dynbal fluxes after land cover change
    !
    ! Should be called AFTER all subgrid weight updates this time step
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_nolakec
    integer                  , intent(in)    :: filter_nolakec(:)
    integer                  , intent(in)    :: num_lakec
    integer                  , intent(in)    :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)    :: urbanparams_inst
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(water_type)         , intent(inout) :: water_inst
    type(temperature_type)   , intent(inout) :: temperature_inst
    type(energyflux_type)    , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: i
    integer  :: g     ! grid cell index
    real(r8) :: this_delta_liq(bounds%begg:bounds%endg)  ! change in gridcell h2o liq content for bulk or one tracer
    real(r8) :: delta_liq_bulk(bounds%begg:bounds%endg)  ! change in gridcell h2o liq content for bulk water
    real(r8) :: delta_heat(bounds%begg:bounds%endg) ! change in gridcell heat content
    !---------------------------------------------------------------------------

    associate( &
         begg => bounds%begg, &
         endg => bounds%endg)

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       associate(bulk_or_tracer => water_inst%bulk_and_tracers(i))
       call dyn_water_content_final(bounds, &
            num_nolakec, filter_nolakec, &
            num_lakec, filter_lakec, &
            bulk_or_tracer%waterstate_inst, &
            bulk_or_tracer%waterdiagnostic_inst, &
            bulk_or_tracer%waterbalance_inst, &
            bulk_or_tracer%waterflux_inst, &
            delta_liq = this_delta_liq(begg:endg))
       if (i == water_inst%i_bulk) then
          delta_liq_bulk(begg:endg) = this_delta_liq(begg:endg)
       end if
       end associate
    end do

    call dyn_heat_content( bounds,                                &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         urbanparams_inst, soilstate_inst, &
         temperature_inst, &
         water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
         heat_grc = temperature_inst%heat2_grc(begg:endg), &
         liquid_water_temp_grc = temperature_inst%liquid_water_temp2_grc(begg:endg))

    if (get_for_testing_zero_dynbal_fluxes()) then
       do g = begg, endg
          delta_heat(g) = 0._r8
       end do
    else
       do g = begg, endg
          delta_heat(g) = temperature_inst%heat2_grc(g) - temperature_inst%heat1_grc(g)
       end do
    end if

    call AdjustDeltaHeatForDeltaLiq( &
         bounds, &
         delta_liq = delta_liq_bulk(begg:endg), &
         liquid_water_temp1 = temperature_inst%liquid_water_temp1_grc(begg:endg), &
         liquid_water_temp2 = temperature_inst%liquid_water_temp2_grc(begg:endg), &
         delta_heat = delta_heat(begg:endg))

    call energyflux_inst%eflx_dynbal_dribbler%set_curr_delta(bounds, &
         delta_heat(begg:endg))
    call energyflux_inst%eflx_dynbal_dribbler%get_curr_flux(bounds, &
         energyflux_inst%eflx_dynbal_grc(begg:endg))

    end associate

  end subroutine dyn_hwcontent_final

  !-----------------------------------------------------------------------
  subroutine dyn_water_content_final(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       waterstate_inst, waterdiagnostic_inst, &
       waterbalance_inst, waterflux_inst, &
       delta_liq)
    !
    ! !DESCRIPTION:
    ! Compute grid cell-level water content and dynbal fluxes after landcover change, for
    ! a single water tracer or bulk water
    !
    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_nolakec
    integer                     , intent(in)    :: filter_nolakec(:)
    integer                     , intent(in)    :: num_lakec
    integer                     , intent(in)    :: filter_lakec(:)
    class(waterstate_type)      , intent(in)    :: waterstate_inst
    class(waterdiagnostic_type) , intent(in)    :: waterdiagnostic_inst
    class(waterbalance_type)    , intent(inout) :: waterbalance_inst
    class(waterflux_type)       , intent(inout) :: waterflux_inst
    real(r8)                    , intent(out)   :: delta_liq(bounds%begg:)  ! change in gridcell h2o liq content
    !
    ! !LOCAL VARIABLES:
    integer  :: g
    real(r8) :: delta_ice(bounds%begg:bounds%endg)  ! change in gridcell h2o ice content

    character(len=*), parameter :: subname = 'dyn_water_content_final'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(delta_liq) == [bounds%endg]), errMsg(sourcefile, __LINE__))

    associate( &
         begg => bounds%begg, &
         endg => bounds%endg)

    call dyn_water_content(bounds, &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         waterstate_inst, waterdiagnostic_inst, &
         liquid_mass = waterbalance_inst%liq2_grc(bounds%begg:bounds%endg), &
         ice_mass    = waterbalance_inst%ice2_grc(bounds%begg:bounds%endg))

    if (get_for_testing_zero_dynbal_fluxes()) then
       do g = begg, endg
          delta_liq(g) = 0._r8
          delta_ice(g) = 0._r8
       end do
    else
       do g = begg, endg
          delta_liq(g)  = waterbalance_inst%liq2_grc(g) - waterbalance_inst%liq1_grc(g)
          delta_ice(g)  = waterbalance_inst%ice2_grc(g) - waterbalance_inst%ice1_grc(g)
       end do
    end if

    call waterflux_inst%qflx_liq_dynbal_dribbler%set_curr_delta(bounds, &
         delta_liq(begg:endg))
    call waterflux_inst%qflx_liq_dynbal_dribbler%get_curr_flux(bounds, &
         waterflux_inst%qflx_liq_dynbal_grc(begg:endg))

    call waterflux_inst%qflx_ice_dynbal_dribbler%set_curr_delta(bounds, &
         delta_ice(begg:endg))
    call waterflux_inst%qflx_ice_dynbal_dribbler%get_curr_flux(bounds, &
         waterflux_inst%qflx_ice_dynbal_grc(begg:endg))

    end associate

  end subroutine dyn_water_content_final

  !-----------------------------------------------------------------------
  subroutine dyn_water_content(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       waterstate_inst, waterdiagnostic_inst, &
       liquid_mass, ice_mass)
    !
    ! !DESCRIPTION:
    ! Compute gridcell total liquid and ice water contents
    !
    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)  :: bounds  
    integer                     , intent(in)  :: num_nolakec
    integer                     , intent(in)  :: filter_nolakec(:)
    integer                     , intent(in)  :: num_lakec
    integer                     , intent(in)  :: filter_lakec(:)
    class(waterstate_type)      , intent(in)  :: waterstate_inst
    class(waterdiagnostic_type) , intent(in)  :: waterdiagnostic_inst
    real(r8)                    , intent(out) :: liquid_mass( bounds%begg: ) ! kg m-2
    real(r8)                    , intent(out) :: ice_mass( bounds%begg: )    ! kg m-2
    !
    ! !LOCAL VARIABLES:
    real(r8) :: liquid_mass_col(bounds%begc:bounds%endc) ! kg m-2
    real(r8) :: ice_mass_col(bounds%begc:bounds%endc)    ! kg m-2

    character(len=*), parameter :: subname = 'dyn_water_content'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(liquid_mass) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(ice_mass) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))

    call ComputeLiqIceMassNonLake(bounds, num_nolakec, filter_nolakec, &
         waterstate_inst, waterdiagnostic_inst, &
         subtract_dynbal_baselines = .true., &
         liquid_mass = liquid_mass_col(bounds%begc:bounds%endc), &
         ice_mass = ice_mass_col(bounds%begc:bounds%endc))

    call ComputeLiqIceMassLake(bounds, num_lakec, filter_lakec, &
         waterstate_inst, &
         subtract_dynbal_baselines = .true., &
         liquid_mass = liquid_mass_col(bounds%begc:bounds%endc), &
         ice_mass = ice_mass_col(bounds%begc:bounds%endc))

    call c2g(bounds, &
         carr = liquid_mass_col(bounds%begc:bounds%endc), &
         garr = liquid_mass(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

    call c2g(bounds, &
         carr = ice_mass_col(bounds%begc:bounds%endc), &
         garr = ice_mass(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

  end subroutine dyn_water_content


  !---------------------------------------------------------------------------
  subroutine dyn_heat_content(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       urbanparams_inst, soilstate_inst, &
       temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
       heat_grc, liquid_water_temp_grc)

    ! !DESCRIPTION:
    ! Compute grid-level heat and water content to track conservation with respect to
    ! dynamic land cover.
    !
    ! Heat content is computed relative to a baseline of 0 C. So temperatures above 0 C
    ! lead to a positive heat content, temperatures below 0 C lead to a negative heat
    ! content. For water, the baseline is considered to be ice at 0 C, so for liquid water
    ! we include the latent heat of fusion.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)  :: bounds  
    integer                  , intent(in)  :: num_nolakec
    integer                  , intent(in)  :: filter_nolakec(:)
    integer                  , intent(in)  :: num_lakec
    integer                  , intent(in)  :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)  :: urbanparams_inst
    type(soilstate_type)     , intent(in)  :: soilstate_inst
    type(temperature_type)   , intent(in)  :: temperature_inst
    type(waterstatebulk_type)    , intent(in)  :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)    , intent(in)  :: waterdiagnosticbulk_inst

    real(r8)                 , intent(out) :: heat_grc( bounds%begg: ) ! total heat content for each grid cell [J/m^2]
    real(r8)                 , intent(out) :: liquid_water_temp_grc( bounds%begg: ) ! weighted average liquid water temperature for each grid cell (K)

    !
    ! !LOCAL VARIABLES:
    integer  :: g

    real(r8) :: heat_col(bounds%begc:bounds%endc)  ! sum of heat content for all columns [J/m^2]
    real(r8) :: heat_liquid_col(bounds%begc:bounds%endc) ! sum of heat content for all columns: liquid water, excluding latent heat [J/m^2]
    real(r8) :: cv_liquid_col(bounds%begc:bounds%endc) ! sum of liquid heat capacity for all columns [J/(m^2 K)]

    real(r8) :: heat_liquid_grc(bounds%begg:bounds%endg) ! heat_liquid_col averaged to grid cell [J/m^2]
    real(r8) :: cv_liquid_grc(bounds%begg:bounds%endg) ! cv_liquid_col averaged to grid cell [J/(m^2 K)]
    !-------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(heat_grc) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(liquid_water_temp_grc) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))

    call ComputeHeatNonLake(bounds, num_nolakec, filter_nolakec, &
         urbanparams_inst, soilstate_inst, &
         temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
         heat = heat_col(bounds%begc:bounds%endc), &
         heat_liquid = heat_liquid_col(bounds%begc:bounds%endc), &
         cv_liquid = cv_liquid_col(bounds%begc:bounds%endc))

    call ComputeHeatLake(bounds, num_lakec, filter_lakec, &
         soilstate_inst, temperature_inst, waterstatebulk_inst, &
         heat = heat_col(bounds%begc:bounds%endc), &
         heat_liquid = heat_liquid_col(bounds%begc:bounds%endc), &
         cv_liquid = cv_liquid_col(bounds%begc:bounds%endc))

    call c2g(bounds, &
         carr = heat_col(bounds%begc:bounds%endc), &
         garr = heat_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

    call c2g(bounds, &
         carr = heat_liquid_col(bounds%begc:bounds%endc), &
         garr = heat_liquid_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

    call c2g(bounds, &
         carr = cv_liquid_col(bounds%begc:bounds%endc), &
         garr = cv_liquid_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'urbanf', &
         l2g_scale_type = 'unity')

    do g = bounds%begg, bounds%endg
       if (cv_liquid_grc(g) > 0._r8) then
          liquid_water_temp_grc(g) = &
               (heat_liquid_grc(g) / cv_liquid_grc(g)) + heat_base_temp
       else
          ! 0 or negative water mass in this grid cell: set an arbitrary temperature
          liquid_water_temp_grc(g) = tfrz
       end if
    end do

  end subroutine dyn_heat_content

end module dynConsBiogeophysMod
