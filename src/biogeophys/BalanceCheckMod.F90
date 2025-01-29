module BalanceCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Water and energy balance check.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use decompMod          , only : bounds_type, get_global_index
  use decompMod          , only : subgrid_level_gridcell, subgrid_level_column, subgrid_level_patch
  use abortutils         , only : endrun
  use clm_varctl         , only : iulog
  use clm_varctl         , only : use_fates_planthydro
  use clm_varpar         , only : nlevsoi
  use atm2lndType        , only : atm2lnd_type
  use EnergyFluxType     , only : energyflux_type
  use SolarAbsorbedType  , only : solarabs_type
  use SoilHydrologyType  , only : soilhydrology_type
  use WaterStateType     , only : waterstate_type
  use LakestateType      , only : lakestate_type
  use WaterDiagnosticBulkType, only : waterdiagnosticbulk_type
  use WaterDiagnosticType, only : waterdiagnostic_type
  use Wateratm2lndType   , only : wateratm2lnd_type
  use Waterlnd2atmType   , only : waterlnd2atm_type
  use WaterBalanceType   , only : waterbalance_type
  use WaterFluxType      , only : waterflux_type
  use WaterType          , only : water_type
  use TotalWaterAndHeatMod, only : ComputeWaterMassNonLake, ComputeWaterMassLake
  use GridcellType       , only : grc                
  use LandunitType       , only : lun                
  use ColumnType         , only : col                
  use PatchType          , only : patch                
  use landunit_varcon    , only : istdlak, istsoil,istcrop,istwet,istice
  use column_varcon      , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon      , only : icol_road_perv, icol_road_imperv
  use clm_varctl         , only : use_hillslope_routing
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  public :: BalanceCheckInit         ! Initialization of Water and energy balance check
  public :: WaterGridcellBalance     ! Grid cell-level water balance check
  public :: BeginWaterColumnBalance  ! Initialize column-level water balance check
  public :: BalanceCheck             ! Water & energy balance checks
  public :: GetBalanceCheckSkipSteps ! Get the number of steps to skip for the balance check
  public :: BalanceCheckClean        ! Clean up for BalanceCheck

  ! !PRIVATE MEMBER DATA:
  real(r8), private, parameter :: skip_size = 3600.0_r8   ! Time steps to skip the balance check at startup (sec)
  integer,  private            :: skip_steps = -999       ! Number of time steps to skip the balance check for

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: WaterGridcellBalanceSingle  ! Grid cell-level water balance check for bulk or a single tracer
  private :: BeginWaterColumnBalanceSingle  ! Initialize column-level water balance check for bulk or a single tracer

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BalanceCheckInit( )
  !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Initialize balance check
    !
    ! !USES:
    use spmdMod           , only : masterproc
    use clm_time_manager  , only : get_step_size_real
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dtime                    ! land model time step (sec)
    !-----------------------------------------------------------------------
    dtime = get_step_size_real()
    ! Skip a minimum of two time steps, but otherwise skip the number of time-steps in the skip_size rounded to the nearest integer
    ! Add an additional step as now required to be after the hourly radiation time-step see github issue #1563
    skip_steps = max(2, nint( (skip_size / dtime) ) ) + 1

    if ( masterproc ) write(iulog,*) ' Skip balance checking for the first ', skip_steps, ' time steps'

  end subroutine BalanceCheckInit

  !-----------------------------------------------------------------------
  subroutine BalanceCheckClean( )
  !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Clean up BalanceCheck
    !
    ! !USES:
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------
    skip_steps = -999
  end subroutine BalanceCheckClean

  !-----------------------------------------------------------------------
  function GetBalanceCheckSkipSteps( ) result( get_skip )
  !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Get the number of steps to skip for the balance check
    !
    ! !ARGUMENTS:
    integer :: get_skip    ! function result
    ! !LOCAL VARIABLES:
    if ( skip_steps > 0 )then
       get_skip = skip_steps
    else
       get_skip = -999
       call endrun('ERROR:: GetBalanceCheckSkipSteps called before BalanceCheckInit')
    end if
  end function GetBalanceCheckSkipSteps
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine WaterGridcellBalance(bounds, &
       num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
       water_inst, lakestate_inst, use_aquifer_layer, flag)
    !
    ! !DESCRIPTION:
    ! Grid cell-level water balance for bulk water and each water tracer
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_nolakec        ! number of column non-lake points in column filter
    integer                 , intent(in)    :: filter_nolakec(:)  ! column filter for non-lake points
    integer                 , intent(in)    :: num_lakec          ! number of column lake points in column filter
    integer                 , intent(in)    :: filter_lakec(:)    ! column filter for lake points
    type(water_type)        , intent(inout) :: water_inst
    type(lakestate_type)    , intent(in)    :: lakestate_inst
    logical                 , intent(in)    :: use_aquifer_layer  ! whether an aquifer layer is used in this run
    character(len=5)        , intent(in)    :: flag  ! specifies begwb or endwb
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'WaterGridcellBalance'
    !-----------------------------------------------------------------------

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       ! Obtain begwb_grc or endwb_grc
       call WaterGridcellBalanceSingle(bounds, &
            num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
            lakestate_inst, &
            water_inst%bulk_and_tracers(i)%waterstate_inst, &
            water_inst%bulk_and_tracers(i)%waterdiagnostic_inst, &
            water_inst%bulk_and_tracers(i)%waterbalance_inst, &
            water_inst%bulk_and_tracers(i)%waterflux_inst, &
            use_aquifer_layer = use_aquifer_layer, flag = flag)
    end do

  end subroutine WaterGridcellBalance

  !-----------------------------------------------------------------------
  subroutine BeginWaterColumnBalance(bounds, &
       num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
       water_inst, soilhydrology_inst, lakestate_inst, &
       use_aquifer_layer)
    !
    ! !DESCRIPTION:
    ! Initialize column-level water balance at beginning of time step, for bulk water and
    ! each water tracer
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                   , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                   , intent(in)    :: num_lakec            ! number of column lake points in column filter
    integer                   , intent(in)    :: filter_lakec(:)      ! column filter for lake points
    type(water_type)          , intent(inout) :: water_inst
    type(lakestate_type)      , intent(in)    :: lakestate_inst
    type(soilhydrology_type)  , intent(in)    :: soilhydrology_inst
    logical                   , intent(in)    :: use_aquifer_layer    ! whether an aquifer layer is used in this run
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'BeginWaterColumnBalance'
    !-----------------------------------------------------------------------

    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
       call BeginWaterColumnBalanceSingle(bounds, &
            num_nolakec, filter_nolakec, &
            num_lakec, filter_lakec, &
            soilhydrology_inst, &
            lakestate_inst, &
            water_inst%bulk_and_tracers(i)%waterstate_inst, &
            water_inst%bulk_and_tracers(i)%waterdiagnostic_inst, &
            water_inst%bulk_and_tracers(i)%waterbalance_inst, &
            use_aquifer_layer = use_aquifer_layer)
    end do

  end subroutine BeginWaterColumnBalance

  !-----------------------------------------------------------------------
  subroutine WaterGridcellBalanceSingle(bounds, &
       num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
       lakestate_inst, waterstate_inst, waterdiagnostic_inst, &
       waterbalance_inst, waterflux_inst, use_aquifer_layer, flag)
    !
    ! !DESCRIPTION:
    ! Grid cell-level water balance for bulk or a single tracer
    ! at beginning or end of time step as specified by "flag"
    !
    ! !USES:
    use subgridAveMod, only: c2g
    use LandunitType , only : lun
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    integer                    , intent(in)    :: num_nolakec        ! number of column non-lake points in column filter
    integer                    , intent(in)    :: filter_nolakec(:)  ! column filter for non-lake points
    integer                    , intent(in)    :: num_lakec          ! number of column lake points in column filter
    integer                    , intent(in)    :: filter_lakec(:)    ! column filter for lake points
    type(lakestate_type)       , intent(in)    :: lakestate_inst
    class(waterstate_type)     , intent(inout) :: waterstate_inst
    class(waterdiagnostic_type), intent(in)    :: waterdiagnostic_inst
    class(waterbalance_type)   , intent(inout) :: waterbalance_inst
    class(waterflux_type)      , intent(inout) :: waterflux_inst
    logical                    , intent(in)    :: use_aquifer_layer  ! whether an aquifer layer is used in this run
    character(len=5)           , intent(in)    :: flag  ! specifies begwb or endwb
    !
    ! !LOCAL VARIABLES:
    integer :: g, l  ! indices
    integer :: begc, endc, begl, endl, begg, endg  ! bounds
    real(r8) :: wb_col(bounds%begc:bounds%endc)  ! temporary column-level water mass
    real(r8) :: wb_grc(bounds%begg:bounds%endg)  ! temporary grid cell-level water mass
    real(r8) :: qflx_liq_dynbal_left_to_dribble(bounds%begg:bounds%endg)  ! grc liq dynamic land cover change conversion runoff flux
    real(r8) :: qflx_ice_dynbal_left_to_dribble(bounds%begg:bounds%endg)  ! grc ice dynamic land cover change conversion runoff flux
    real(r8) :: wa_reset_nonconservation_gain_grc(bounds%begg:bounds%endg)  ! grc mass gained from resetting water in the unconfined aquifer, wa_col (negative indicates mass lost) (mm)

    character(len=*), parameter :: subname = 'WaterGridcellBalanceSingle'
    !-----------------------------------------------------------------------

    associate(                                     &
         begwb_grc => waterbalance_inst%begwb_grc, &  ! Output: [real(r8) (:)]  grid cell-level water mass begining of the time step
         endwb_grc => waterbalance_inst%endwb_grc, &  ! Output: [real(r8) (:)]  grid cell-level water mass end of the time step
         wa_reset_nonconservation_gain_col => waterbalance_inst%wa_reset_nonconservation_gain_col                                  &  ! Input:  [real(r8) (:)]  col mass gained from resetting water in the unconfined aquifer, wa_col (negative indicates mass lost) (mm)
         )

    begc = bounds%begc
    endc = bounds%endc
    begl = bounds%begl
    endl = bounds%endl
    begg = bounds%begg
    endg = bounds%endg

    call ComputeWaterMassNonLake(bounds, num_nolakec, filter_nolakec, &
         waterstate_inst, waterdiagnostic_inst, &
         subtract_dynbal_baselines = .true., &
         water_mass = wb_col(begc:endc))

    call ComputeWaterMassLake(bounds, num_lakec, filter_lakec, &
         waterstate_inst, lakestate_inst, &
         add_lake_water_and_subtract_dynbal_baselines = .true., &
         water_mass = wb_col(begc:endc))

    call c2g(bounds, wb_col(begc:endc), wb_grc(begg:endg), &
             c2l_scale_type='urbanf', l2g_scale_type='unity')

    ! add landunit level state variable, convert from (m3) to (kg m-2)
    if (use_hillslope_routing) then
       do l = begl, endl
          g = lun%gridcell(l)
          wb_grc(g) = wb_grc(g) +  waterstate_inst%stream_water_volume_lun(l) &
               *1e3_r8/(grc%area(g)*1.e6_r8)
       enddo
    endif
    
    ! Call the beginning or ending version of the subroutine according
    ! to flag value
    if (flag == 'begwb') then
       call waterflux_inst%qflx_liq_dynbal_dribbler%get_amount_left_to_dribble_beg( &
         bounds, &
         qflx_liq_dynbal_left_to_dribble(begg:endg))
       call waterflux_inst%qflx_ice_dynbal_dribbler%get_amount_left_to_dribble_beg( &
         bounds, &
         qflx_ice_dynbal_left_to_dribble(begg:endg))
    else if (flag == 'endwb') then
       call waterflux_inst%qflx_liq_dynbal_dribbler%get_amount_left_to_dribble_end( &
         bounds, &
         qflx_liq_dynbal_left_to_dribble(begg:endg))
       call waterflux_inst%qflx_ice_dynbal_dribbler%get_amount_left_to_dribble_end( &
         bounds, &
         qflx_ice_dynbal_left_to_dribble(begg:endg))
    else
       write(iulog,*) 'Unknown flag passed into this subroutine.'
       write(iulog,*) 'Expecting either begwb or endwb.'
       call endrun(msg=errmsg(sourcefile, __LINE__))
    end if

    ! These dynbal dribblers store the delta state, (end - beg). Thus, the
    ! amount dribbled out is the negative of the amount stored in the
    ! dribblers. Therefore, conservation requires us to subtract the amount
    ! remaining to dribble.
    ! This sign convention is opposite to the convention chosen for the
    ! respective dribble terms used in the carbon balance. At some point
    ! it may be worth making the two conventions consistent.
    do g = begg, endg
       wb_grc(g) = wb_grc(g) - qflx_liq_dynbal_left_to_dribble(g)  &
                             - qflx_ice_dynbal_left_to_dribble(g)
    end do

    ! Map wb_grc to beginning/ending water balance according to flag
    if (flag == 'begwb') then
       do g = begg, endg
          begwb_grc(g) = wb_grc(g)
       end do
    else if (flag == 'endwb') then
       ! endwb_grc requires one more step first
       if (use_aquifer_layer) then
          ! wa_reset_nonconservation_gain may be non-zero only when
          ! use_aquifer_layer is true. We do this c2g call only when needed
          ! to avoid unnecessary calculations; by adding this term only when
          ! use_aquifer_layer is true, we effectively let the balance checks
          ! ensure that this term is zero when use_aquifer_layer is false,
          ! as it should be.
          ! The _col term was determined in BeginWaterColumnBalanceSingle
          ! after any dynamic landuse adjustments.
          call c2g( bounds, &
               wa_reset_nonconservation_gain_col(begc:endc), &
               wa_reset_nonconservation_gain_grc(begg:endg), &
               c2l_scale_type='urbanf', l2g_scale_type='unity' )
       else
          wa_reset_nonconservation_gain_grc(begg:endg) = 0._r8
       end if
       do g = begg, endg
          endwb_grc(g) = wb_grc(g) - wa_reset_nonconservation_gain_grc(g)
       end do
    end if

    end associate

  end subroutine WaterGridcellBalanceSingle

   !-----------------------------------------------------------------------
  subroutine BeginWaterColumnBalanceSingle(bounds, &
       num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
       soilhydrology_inst, lakestate_inst, waterstate_inst, & 
       waterdiagnostic_inst, waterbalance_inst, &
       use_aquifer_layer)
    !
    ! !DESCRIPTION:
    ! Initialize column-level water balance at beginning of time step, for bulk or a
    ! single tracer
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds     
    integer                   , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                   , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                   , intent(in)    :: num_lakec            ! number of column lake points in column filter
    integer                   , intent(in)    :: filter_lakec(:)      ! column filter for lake points
    type(soilhydrology_type)  , intent(in)    :: soilhydrology_inst
    type(lakestate_type)      , intent(in)    :: lakestate_inst
    class(waterstate_type)    , intent(inout) :: waterstate_inst
    class(waterdiagnostic_type), intent(in)   :: waterdiagnostic_inst
    class(waterbalance_type)  , intent(inout) :: waterbalance_inst
    logical                   , intent(in)    :: use_aquifer_layer    ! whether an aquifer layer is used in this run
    !
    ! !LOCAL VARIABLES:
    integer :: c, fc                  ! indices
    !-----------------------------------------------------------------------

    associate(                                               &
         zi         => col%zi                           , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         zwt        => soilhydrology_inst%zwt_col       , & ! Input:  [real(r8) (:)   ]  water table depth (m)
         aquifer_water_baseline => waterstate_inst%aquifer_water_baseline, &  ! Input: [real(r8)] baseline value for water in the unconfined aquifer (wa_col) for this bulk / tracer (mm)
         wa         => waterstate_inst%wa_col           , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)
         wa_reset_nonconservation_gain => waterbalance_inst%wa_reset_nonconservation_gain_col                                           , & ! Output: [real(r8) (:)   ]  mass gained from resetting water in the unconfined aquifer, wa_col (negative indicates mass lost) (mm)
         begwb      => waterbalance_inst%begwb_col      , & ! Output: [real(r8) (:)   ]  water mass begining of the time step
         h2osno_old => waterbalance_inst%h2osno_old_col   & ! Output: [real(r8) (:)   ]  snow water (mm H2O) at previous time step
         )

    ! wa(c) gets added to liquid_mass in ComputeLiqIceMassNonLake called here.
    ! wa_reset_nonconservation_gain is calculated for the grid cell-level
    ! water balance check and may be non-zero only when
    ! use_aquifer_layer is true. The grid cell-level balance check ensures
    ! that this term is zero when use_aquifer_layer is false, as it should be.
    ! In particular, we adjust wa back to the baseline under certain
    ! conditions. The right way to do this might be to use explicit fluxes from
    ! some other state, but in this case we don't have a source to pull from,
    ! so we adjust wa without explicit fluxes. Because we do this before
    ! initializing the column-level balance check, the column-level check is
    ! unaware of the adjustment. However, since this adjustment happens after
    ! initializing the gridcell-level balance check, we have to account for
    ! it in the gridcell-level balance check. The normal way to account for an
    ! adjustment like this would be to include the flux in the balance check.
    ! Here we don't have an explicit flux, so instead we track the
    ! non-conservation state. In principle, we could calculate an explicit flux
    ! and use that, but we don't gain anything from using an explicit flux in
    ! this case.
    if(use_aquifer_layer) then
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          if (col%hydrologically_active(c)) then
             if(zwt(c) <= zi(c,nlevsoi)) then
                wa_reset_nonconservation_gain(c) = aquifer_water_baseline - &
                                                   wa(c)
                wa(c) = aquifer_water_baseline
             else
                wa_reset_nonconservation_gain(c) = 0._r8
             end if
          end if
       end do
    endif

    call ComputeWaterMassNonLake(bounds, num_nolakec, filter_nolakec, &
         waterstate_inst, waterdiagnostic_inst, &
         subtract_dynbal_baselines = .false., &
         water_mass = begwb(bounds%begc:bounds%endc))

    call ComputeWaterMassLake(bounds, num_lakec, filter_lakec, &
         waterstate_inst, lakestate_inst, &
         add_lake_water_and_subtract_dynbal_baselines = .false., &
         water_mass = begwb(bounds%begc:bounds%endc))

    call waterstate_inst%CalculateTotalH2osno(bounds, num_nolakec, filter_nolakec, &
         caller = 'BeginWaterBalanceSingle-nolake', &
         h2osno_total = h2osno_old(bounds%begc:bounds%endc))
    call waterstate_inst%CalculateTotalH2osno(bounds, num_lakec, filter_lakec, &
         caller = 'BeginWaterBalanceSingle-lake', &
         h2osno_total = h2osno_old(bounds%begc:bounds%endc))

    end associate 

  end subroutine BeginWaterColumnBalanceSingle

   !-----------------------------------------------------------------------
   subroutine BalanceCheck( bounds, &
        num_allc, filter_allc, &
        atm2lnd_inst, solarabs_inst, waterflux_inst, waterstate_inst, &
        waterdiagnosticbulk_inst, waterbalance_inst, wateratm2lnd_inst, &
        waterlnd2atm_inst, surfalb_inst, energyflux_inst, canopystate_inst)
     !
     ! !DESCRIPTION:
     ! This subroutine accumulates the numerical truncation errors of the water
     ! and energy balance calculation. It is helpful to see the performance of
     ! the process of integration.
     !
     ! The error for energy balance:
     !
     ! error = abs(Net radiation - change of internal energy - Sensible heat
     !             - Latent heat)
     !
     ! The error for water balance:
     !
     ! error = abs(precipitation - change of water storage - evaporation - runoff)
     !
     ! !USES:
     use clm_varcon        , only : spval
     use clm_varctl        , only : use_soil_moisture_streams
     use clm_time_manager  , only : get_step_size_real, get_nstep
     use clm_time_manager  , only : get_nstep_since_startup_or_lastDA_restart_or_pause
     use CanopyStateType   , only : canopystate_type
     use subgridAveMod     , only : c2g
     use dynSubgridControlMod, only : get_for_testing_zero_dynbal_fluxes
     use SurfaceAlbedoType , only : surfalb_type
     !
     ! !ARGUMENTS:
     type(bounds_type)     , intent(in)    :: bounds  
     integer               , intent(in)    :: num_allc        ! number of columns in allc filter
     integer               , intent(in)    :: filter_allc(:)  ! filter for all columns
     type(atm2lnd_type)    , intent(in)    :: atm2lnd_inst
     type(solarabs_type)   , intent(in)    :: solarabs_inst
     class(waterflux_type) , intent(in)    :: waterflux_inst
     class(waterstate_type), intent(in)    :: waterstate_inst
     type(waterdiagnosticbulk_type), intent(in) :: waterdiagnosticbulk_inst
     class(waterbalance_type), intent(inout) :: waterbalance_inst
     class(waterlnd2atm_type), intent(in) :: waterlnd2atm_inst
     class(wateratm2lnd_type) , intent(in) :: wateratm2lnd_inst
     type(surfalb_type)    , intent(in)    :: surfalb_inst
     type(energyflux_type) , intent(inout) :: energyflux_inst
     type(canopystate_type), intent(inout) :: canopystate_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: p,c,l,g,fc                             ! indices
     real(r8) :: dtime                                  ! land model time step (sec)
     integer  :: nstep                                  ! time step number
     integer  :: DAnstep                                ! time step number since last Data Assimilation (DA)
     integer  :: indexp,indexc,indexl,indexg            ! index of first found in search loop
     integer  :: global_index                           ! index in global index space
     real(r8) :: errh2o_grc(bounds%begg:bounds%endg)    ! grid cell level water conservation error [mm H2O]
     real(r8) :: forc_rain_col(bounds%begc:bounds%endc) ! column level rain rate [mm/s]
     real(r8) :: forc_snow_col(bounds%begc:bounds%endc) ! column level snow rate [mm/s]
     real(r8) :: h2osno_total(bounds%begc:bounds%endc)  ! total snow water [mm H2O]
     real(r8) :: qflx_glcice_dyn_water_flux_grc(bounds%begg:bounds%endg)  ! grid cell-level water flux needed for balance check due to glc_dyn_runoff_routing [mm H2O/s] (positive means addition of water to the system)
     real(r8) :: qflx_snwcp_discarded_liq_grc(bounds%begg:bounds%endg)  ! grid cell-level excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack [mm H2O /s]
     real(r8) :: qflx_snwcp_discarded_ice_grc(bounds%begg:bounds%endg)  ! grid cell-level excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack [mm H2O /s]

     real(r8) :: errh2o_max_val                         ! Maximum value of error in water conservation error  over all columns [mm H2O]
     real(r8) :: errh2osno_max_val                      ! Maximum value of error in h2osno conservation error over all columns [kg m-2]
     real(r8) :: errsol_max_val                         ! Maximum value of error in solar radiation conservation error over all columns [W m-2]
     real(r8) :: errlon_max_val                         ! Maximum value of error in longwave radiation conservation error over all columns [W m-2]
     real(r8) :: errseb_max_val                         ! Maximum value of error in surface energy conservation error over all columns [W m-2]
     real(r8) :: errsoi_col_max_val                     ! Maximum value of column-level soil/lake energy conservation error over all columns [W m-2]

     real(r8), parameter :: h2o_warning_thresh       = 1.e-9_r8                       ! Warning threshhold for error in errh2o and errh2osnow 
     real(r8), parameter :: energy_warning_thresh    = 1.e-7_r8                       ! Warning threshhold for error in errsol, errsol, errseb, errlonv
     real(r8), parameter :: error_thresh             = 1.e-5_r8                       ! Error threshhold for conservation error

     !-----------------------------------------------------------------------

     associate(                                                                   & 
          forc_solad_col    =>    atm2lnd_inst%forc_solad_downscaled_col        , & ! Input:  [real(r8) (:,:) ]  direct beam radiation (vis=forc_sols , nir=forc_soll )
          forc_solad        =>    atm2lnd_inst%forc_solad_not_downscaled_grc    , & ! Input:  [real(r8) (:,:) ]  direct beam radiation (vis=forc_sols , nir=forc_soll )
          forc_solai        =>    atm2lnd_inst%forc_solai_grc                   , & ! Input:  [real(r8) (:,:) ]  diffuse radiation     (vis=forc_solsd, nir=forc_solld)
          forc_rain         =>    wateratm2lnd_inst%forc_rain_downscaled_col    , & ! Input:  [real(r8) (:)   ]  column level rain rate [mm/s]
          forc_rain_grc     =>    wateratm2lnd_inst%forc_rain_not_downscaled_grc, & ! Input:  [real(r8) (:)   ]  grid cell-level rain rate [mm/s]
          forc_snow         =>    wateratm2lnd_inst%forc_snow_downscaled_col    , & ! Input:  [real(r8) (:)   ]  column level snow rate [mm/s]
          forc_snow_grc     =>    wateratm2lnd_inst%forc_snow_not_downscaled_grc, & ! Input:  [real(r8) (:)   ]  grid cell-level snow rate [mm/s]
          forc_lwrad              =>    atm2lnd_inst%forc_lwrad_downscaled_col  , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)

          h2osno_old              =>    waterbalance_inst%h2osno_old_col          , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O) at previous time step
          frac_sno_eff            =>    waterdiagnosticbulk_inst%frac_sno_eff_col        , & ! Input:  [real(r8) (:)   ]  effective snow fraction                 
          frac_sno                =>    waterdiagnosticbulk_inst%frac_sno_col            , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          snow_depth              =>    waterdiagnosticbulk_inst%snow_depth_col          , & ! Input:  [real(r8) (:)   ]  snow height (m)                         
          begwb_grc               =>    waterbalance_inst%begwb_grc             , & ! Input:  [real(r8) (:)   ]  grid cell-level water mass begining of the time step
          endwb_grc               =>    waterbalance_inst%endwb_grc             , & ! Output: [real(r8) (:)   ]  grid cell-level water mass end of the time step
          begwb_col               =>    waterbalance_inst%begwb_col             , & ! Input:  [real(r8) (:)   ]  column-level water mass begining of the time step
          endwb_col               =>    waterbalance_inst%endwb_col             , & ! Output: [real(r8) (:)   ]  column-level water mass end of the time step
          errh2o_col              =>    waterbalance_inst%errh2o_col            , & ! Output: [real(r8) (:)   ]  column-level water conservation error (mm H2O)
          errh2osno               =>    waterbalance_inst%errh2osno_col           , & ! Output: [real(r8) (:)   ]  error in h2osno (kg m-2)                
          snow_sources            =>    waterbalance_inst%snow_sources_col         , & ! Output: [real(r8) (:)   ]  snow sources (mm H2O /s)
          snow_sinks              =>    waterbalance_inst%snow_sinks_col           , & ! Output: [real(r8) (:)   ]  snow sinks (mm H2O /s)
          qflx_liq_grnd_col       =>    waterflux_inst%qflx_liq_grnd_col        , & ! Input:  [real(r8) (:)   ]  liquid on ground after interception (mm H2O/s) [+]
          qflx_snow_grnd_col      =>    waterflux_inst%qflx_snow_grnd_col       , & ! Input:  [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
          qflx_snwcp_liq          =>    waterflux_inst%qflx_snwcp_liq_col       , & ! Input:  [real(r8) (:)   ]  excess liquid h2o due to snow capping (outgoing) (mm H2O /s) [+]`
          qflx_snwcp_ice          =>    waterflux_inst%qflx_snwcp_ice_col       , & ! Input:  [real(r8) (:)   ]  excess solid h2o due to snow capping (outgoing) (mm H2O /s) [+]`
          qflx_snwcp_discarded_liq_col => waterflux_inst%qflx_snwcp_discarded_liq_col, & ! Input: [real(r8) (:)] column level excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s) [+]
          qflx_snwcp_discarded_ice_col => waterflux_inst%qflx_snwcp_discarded_ice_col, & ! Input: [real(r8) (:)] column level excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s) [+]
          qflx_evap_tot_col       =>    waterflux_inst%qflx_evap_tot_col        , & ! Input:  [real(r8) (:)   ]  column level qflx_evap_soi + qflx_evap_can + qflx_tran_veg
          qflx_evap_tot_grc       =>    waterlnd2atm_inst%qflx_evap_tot_grc     , & ! Input:  [real(r8) (:)   ]  grid cell-level qflx_evap_soi + qflx_evap_can + qflx_tran_veg
          qflx_soliddew_to_top_layer    => waterflux_inst%qflx_soliddew_to_top_layer_col   , & ! Input:  [real(r8) (:)   ]  rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s) [+]
          qflx_solidevap_from_top_layer => waterflux_inst%qflx_solidevap_from_top_layer_col, & ! Input:  [real(r8) (:)   ]  rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s) [+]
          qflx_liqevap_from_top_layer   => waterflux_inst%qflx_liqevap_from_top_layer_col  , & ! Input:  [real(r8) (:)   ]  rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]
          qflx_liqdew_to_top_layer      => waterflux_inst%qflx_liqdew_to_top_layer_col     , & ! Input:  [real(r8) (:)   ]  rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s) [+]
          qflx_prec_grnd          =>    waterdiagnosticbulk_inst%qflx_prec_grnd_col, & ! Input:  [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
          qflx_snow_h2osfc        =>    waterflux_inst%qflx_snow_h2osfc_col     , & ! Input:  [real(r8) (:)   ]  snow falling on surface water (mm/s)
          qflx_h2osfc_to_ice      =>    waterflux_inst%qflx_h2osfc_to_ice_col   , & ! Input:  [real(r8) (:)   ]  conversion of h2osfc to ice             
          qflx_drain_perched_col  =>    waterflux_inst%qflx_drain_perched_col   , & ! Input:  [real(r8) (:)   ]  column level sub-surface runoff (mm H2O /s)
          qflx_drain_perched_grc  =>    waterlnd2atm_inst%qflx_rofliq_drain_perched_grc, & ! Input: [real(r8) (:)] grid cell-level sub-surface runoff (mm H2O /s)
          qflx_flood_col          =>    waterflux_inst%qflx_floodc_col          , & ! Input:  [real(r8) (:)   ]  column level total runoff due to flooding
          forc_flood_grc          =>    wateratm2lnd_inst%forc_flood_grc        , & ! Input:  [real(r8) (:)   ]  grid cell-level total grid cell-level runoff from river model
          qflx_snow_drain         =>    waterflux_inst%qflx_snow_drain_col      , & ! Input:  [real(r8) (:)   ]  drainage from snow pack                         
          qflx_surf_col           =>    waterflux_inst%qflx_surf_col            , & ! Input:  [real(r8) (:)   ]  column level surface runoff (mm H2O /s)
          qflx_surf_grc           =>    waterlnd2atm_inst%qflx_rofliq_qsur_grc  , & ! Input:  [real(r8) (:)   ]  grid cell-level surface runoff (mm H20 /s)
          qflx_qrgwl_col          =>    waterflux_inst%qflx_qrgwl_col           , & ! Input:  [real(r8) (:)   ]  column level qflx_surf at glaciers, wetlands, lakes
          qflx_qrgwl_grc          =>    waterlnd2atm_inst%qflx_rofliq_qgwl_grc  , & ! Input:  [real(r8) (:)   ]  grid cell-level qflx_surf at glaciers, wetlands, lakes
          qflx_drain_col          =>    waterflux_inst%qflx_drain_col           , & ! Input:  [real(r8) (:)   ]  column level sub-surface runoff (mm H2O /s)
          qflx_drain_grc          =>    waterlnd2atm_inst%qflx_rofliq_qsub_grc  , & ! Input:  [real(r8) (:)   ]  grid cell-level drainage (mm H20 /s)
          qflx_streamflow_grc     =>    waterlnd2atm_inst%qflx_rofliq_stream_grc, & ! Input: [real(r8) (:)   ] streamflow [mm H2O/s]
          qflx_ice_runoff_col     =>    waterlnd2atm_inst%qflx_ice_runoff_col   , & ! Input:  [real(r8) (:)   ] column level solid runoff from snow capping and from excess ice in soil (mm H2O /s)
          qflx_ice_runoff_grc     =>    waterlnd2atm_inst%qflx_rofice_grc       , & ! Input:  [real(r8) (:)   ] grid cell-level solid runoff from snow capping and from excess ice in soil (mm H2O /s)
          qflx_sl_top_soil        =>    waterflux_inst%qflx_sl_top_soil_col     , & ! Input:  [real(r8) (:)   ]  liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)

          qflx_sfc_irrig_col      =>    waterflux_inst%qflx_sfc_irrig_col       , & ! Input:  [real(r8) (:)   ]  column level irrigation flux (mm H2O /s)
          qflx_sfc_irrig_grc      =>    waterlnd2atm_inst%qirrig_grc            , & ! Input:  [real(r8) (:)   ]  grid cell-level irrigation flux (mm H20 /s)
          qflx_glcice_dyn_water_flux_col => waterflux_inst%qflx_glcice_dyn_water_flux_col, & ! Input: [real(r8) (:)]  column level water flux needed for balance check due to glc_dyn_runoff_routing (mm H2O/s) (positive means addition of water to the system)

          dhsdt_canopy            =>    energyflux_inst%dhsdt_canopy_patch      , & ! Input:  [real(r8) (:)   ]  change in heat content of canopy (W/m**2) [+ to atm]
          eflx_lwrad_out          =>    energyflux_inst%eflx_lwrad_out_patch    , & ! Input:  [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)
          eflx_lwrad_net          =>    energyflux_inst%eflx_lwrad_net_patch    , & ! Input:  [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
          eflx_sh_tot             =>    energyflux_inst%eflx_sh_tot_patch       , & ! Input:  [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
          eflx_lh_tot             =>    energyflux_inst%eflx_lh_tot_patch       , & ! Input:  [real(r8) (:)   ]  total latent heat flux (W/m**2)  [+ to atm]
          eflx_soil_grnd          =>    energyflux_inst%eflx_soil_grnd_patch    , & ! Input:  [real(r8) (:)   ]  soil heat flux (W/m**2) [+ = into soil] 
          eflx_wasteheat_patch    =>    energyflux_inst%eflx_wasteheat_patch    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
          eflx_ventilation_patch  =>    energyflux_inst%eflx_ventilation_patch  , & ! Input:  [real(r8) (:)   ]  sensible heat flux from building ventilation (W/m**2)
          eflx_heat_from_ac_patch =>    energyflux_inst%eflx_heat_from_ac_patch , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
          eflx_traffic_patch      =>    energyflux_inst%eflx_traffic_patch      , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
          eflx_dynbal             =>    energyflux_inst%eflx_dynbal_grc         , & ! Input:  [real(r8) (:)   ]  energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]
          errsoi_col              =>    energyflux_inst%errsoi_col              , & ! Output: [real(r8) (:)   ]  column-level soil/lake energy conservation error (W/m**2)
          errsol                  =>    energyflux_inst%errsol_patch            , & ! Output: [real(r8) (:)   ]  solar radiation conservation error (W/m**2)
          errseb                  =>    energyflux_inst%errseb_patch            , & ! Output: [real(r8) (:)   ]  surface energy conservation error (W/m**2)
          errlon                  =>    energyflux_inst%errlon_patch            , & ! Output: [real(r8) (:)   ]  longwave radiation conservation error (W/m**2)

          sabg_soil               =>    solarabs_inst%sabg_soil_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
          sabg_snow               =>    solarabs_inst%sabg_snow_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
          sabg_chk                =>    solarabs_inst%sabg_chk_patch            , & ! Input:  [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
          fsa                     =>    solarabs_inst%fsa_patch                 , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed (total) (W/m**2)
          fsr                     =>    solarabs_inst%fsr_patch                 , & ! Input:  [real(r8) (:)   ]  solar radiation reflected (W/m**2)      
          sabv                    =>    solarabs_inst%sabv_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)
          sabg                    =>    solarabs_inst%sabg_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
          
          elai                    =>    canopystate_inst%elai_patch             , & ! Input:  [real(r8) (:,:)]  
          esai                    =>    canopystate_inst%esai_patch             , & ! Input:  [real(r8) (:,:)]  

          fabd                    =>    surfalb_inst%fabd_patch                 , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit direct flux
          fabi                    =>    surfalb_inst%fabi_patch                 , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit indirect flux
          albd                    =>    surfalb_inst%albd_patch                 , & ! Input:  [real(r8) (:,:)]  surface albedo (direct)
          albi                    =>    surfalb_inst%albi_patch                 , & ! Input:  [real(r8) (:,:)]  surface albedo (diffuse)
          ftdd                    =>    surfalb_inst%ftdd_patch                 , & ! Input:  [real(r8) (:,:)]  down direct flux below canopy per unit direct flux
          ftid                    =>    surfalb_inst%ftid_patch                 , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit direct flux
          ftii                    =>    surfalb_inst%ftii_patch                 , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit diffuse flux

          netrad                  =>    energyflux_inst%netrad_patch              & ! Output: [real(r8) (:)   ]  net radiation (positive downward) (W/m**2)
          )

       ! Get step size and time step

       nstep = get_nstep()
       DAnstep = get_nstep_since_startup_or_lastDA_restart_or_pause()
       dtime = get_step_size_real()

       ! Determine column level incoming snow and rain
       ! Assume no incident precipitation on urban wall columns (as in CanopyHydrologyMod.F90).

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          l = col%landunit(c)       

          if (col%itype(c) == icol_sunwall .or.  col%itype(c) == icol_shadewall) then
             forc_rain_col(c) = 0._r8
             forc_snow_col(c) = 0._r8
          else
             forc_rain_col(c) = forc_rain(c)
             forc_snow_col(c) = forc_snow(c)
          end if
       end do

       ! Water balance check at the column level

       do c = bounds%begc, bounds%endc

          ! add qflx_drain_perched and qflx_flood
          if (col%active(c)) then

             errh2o_col(c) = endwb_col(c) - begwb_col(c) &
                  - (forc_rain_col(c)        &
                  + forc_snow_col(c)         &
                  + qflx_flood_col(c)        &
                  + qflx_sfc_irrig_col(c)    &
                  + qflx_glcice_dyn_water_flux_col(c) &
                  - qflx_evap_tot_col(c)     &
                  - qflx_surf_col(c)         &
                  - qflx_qrgwl_col(c)        &
                  - qflx_drain_col(c)        &
                  - qflx_drain_perched_col(c) &
                  - qflx_ice_runoff_col(c)   &
                  - qflx_snwcp_discarded_liq_col(c) &
                  - qflx_snwcp_discarded_ice_col(c)) * dtime

          else

             errh2o_col(c) = 0.0_r8

          end if

       end do
       
       errh2o_max_val = maxval(abs(errh2o_col(bounds%begc:bounds%endc)))

       if (errh2o_max_val > h2o_warning_thresh) then

           indexc = maxloc( abs(errh2o_col(bounds%begc:bounds%endc)), 1 ) + bounds%begc - 1
           global_index = get_global_index(subgrid_index=indexc, subgrid_level=subgrid_level_column)
           write(iulog,*)'WARNING:  column-level water balance error ',&
             ' nstep= ',nstep, &
             ' local indexc= ',indexc,&
             ' global indexc= ',global_index, &
             ' errh2o= ',errh2o_col(indexc)
         if ((errh2o_max_val > error_thresh) .and. (DAnstep > skip_steps)) then
              
              write(iulog,*)'CTSM is stopping because errh2o > ', error_thresh, ' mm'
              write(iulog,*)'nstep                     = ',nstep
              write(iulog,*)'errh2o_col                = ',errh2o_col(indexc)
              write(iulog,*)'forc_rain                 = ',forc_rain_col(indexc)*dtime
              write(iulog,*)'forc_snow                 = ',forc_snow_col(indexc)*dtime
              write(iulog,*)'endwb_col                 = ',endwb_col(indexc)
              write(iulog,*)'begwb_col                 = ',begwb_col(indexc)

              write(iulog,*)'qflx_evap_tot             = ',qflx_evap_tot_col(indexc)*dtime
              write(iulog,*)'qflx_sfc_irrig            = ',qflx_sfc_irrig_col(indexc)*dtime
              write(iulog,*)'qflx_surf                 = ',qflx_surf_col(indexc)*dtime
              write(iulog,*)'qflx_qrgwl                = ',qflx_qrgwl_col(indexc)*dtime
              write(iulog,*)'qflx_drain                = ',qflx_drain_col(indexc)*dtime

              write(iulog,*)'qflx_ice_runoff           = ',qflx_ice_runoff_col(indexc)*dtime

              write(iulog,*)'qflx_snwcp_discarded_ice  = ',qflx_snwcp_discarded_ice_col(indexc)*dtime
              write(iulog,*)'qflx_snwcp_discarded_liq  = ',qflx_snwcp_discarded_liq_col(indexc)*dtime
              write(iulog,*)'deltawb                   = ',endwb_col(indexc)-begwb_col(indexc)
              write(iulog,*)'deltawb/dtime             = ',(endwb_col(indexc)-begwb_col(indexc))/dtime

              if (.not.(col%itype(indexc) == icol_roof .or. &
                   col%itype(indexc) == icol_road_imperv .or. &
                   col%itype(indexc) == icol_road_perv)) then
                   write(iulog,*)'qflx_drain_perched         = ',qflx_drain_perched_col(indexc)*dtime
                   write(iulog,*)'qflx_flood                 = ',qflx_flood_col(indexc)*dtime
                   write(iulog,*)'qflx_glcice_dyn_water_flux = ', qflx_glcice_dyn_water_flux_col(indexc)*dtime
              end if
              
              write(iulog,*)'CTSM is stopping'
              call endrun(subgrid_index=indexc, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
         end if
       
       end if
       

       ! Water balance check at the grid cell level

       call c2g( bounds,  &
         qflx_glcice_dyn_water_flux_col(bounds%begc:bounds%endc),  &
         qflx_glcice_dyn_water_flux_grc(bounds%begg:bounds%endg),  &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
       call c2g( bounds,  &
         qflx_snwcp_discarded_liq_col(bounds%begc:bounds%endc),  &
         qflx_snwcp_discarded_liq_grc(bounds%begg:bounds%endg),  &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
       call c2g( bounds,  &
         qflx_snwcp_discarded_ice_col(bounds%begc:bounds%endc),  &
         qflx_snwcp_discarded_ice_grc(bounds%begg:bounds%endg),  &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

       do g = bounds%begg, bounds%endg
          errh2o_grc(g) = endwb_grc(g) - begwb_grc(g)  &
               - (forc_rain_grc(g)  &
               + forc_snow_grc(g)  &
               + forc_flood_grc(g)  &
               + qflx_sfc_irrig_grc(g)  &
               + qflx_glcice_dyn_water_flux_grc(g)  &
               - qflx_evap_tot_grc(g)  &
               - qflx_surf_grc(g)  &
               - qflx_qrgwl_grc(g)  &
               - qflx_drain_grc(g)  &
               - qflx_drain_perched_grc(g)  &
               - qflx_ice_runoff_grc(g)  &
               - qflx_snwcp_discarded_liq_grc(g)  &
               - qflx_snwcp_discarded_ice_grc(g)) * dtime
       end do

       ! add landunit level flux variable, convert from (m3/s) to (kg m-2 s-1)
       if (use_hillslope_routing) then
          ! output water flux from streamflow (+)
          do g = bounds%begg, bounds%endg
             errh2o_grc(g) = errh2o_grc(g) &
                  +  qflx_streamflow_grc(g) * dtime
          enddo
       endif

       errh2o_max_val = maxval(abs(errh2o_grc(bounds%begg:bounds%endg)))

       ! BUG(rgk, 2021-04-13, ESCOMP/CTSM#1314) Temporarily bypassing gridcell-level check with use_fates_planthydro until issue 1314 is resolved
       
       if (errh2o_max_val > h2o_warning_thresh .and. .not.use_fates_planthydro) then

          indexg = maxloc( abs(errh2o_grc(bounds%begg:bounds%endg)), 1 ) + bounds%begg - 1
          write(iulog,*)'WARNING:  grid cell-level water balance error ',&
             ' nstep= ',nstep, &
             ' local indexg= ',indexg,&
             ' errh2o_grc= ',errh2o_grc(indexg)
          if (errh2o_max_val > error_thresh .and. DAnstep > skip_steps .and. &
              .not. use_soil_moisture_streams .and. &
              .not. get_for_testing_zero_dynbal_fluxes()) then

             write(iulog,*)'CTSM is stopping because errh2o > ', error_thresh, ' mm'
             write(iulog,*)'nstep                     = ',nstep
             write(iulog,*)'errh2o_grc                = ',errh2o_grc(indexg)
             write(iulog,*)'forc_rain                 = ',forc_rain_grc(indexg)*dtime
             write(iulog,*)'forc_snow                 = ',forc_snow_grc(indexg)*dtime
             write(iulog,*)'endwb_grc                 = ',endwb_grc(indexg)
             write(iulog,*)'begwb_grc                 = ',begwb_grc(indexg)

             write(iulog,*)'qflx_evap_tot             = ',qflx_evap_tot_grc(indexg)*dtime
             write(iulog,*)'qflx_sfc_irrig            = ',qflx_sfc_irrig_grc(indexg)*dtime
             write(iulog,*)'qflx_surf                 = ',qflx_surf_grc(indexg)*dtime
             write(iulog,*)'qflx_qrgwl                = ',qflx_qrgwl_grc(indexg)*dtime
             write(iulog,*)'qflx_drain                = ',qflx_drain_grc(indexg)*dtime
             write(iulog,*)'qflx_ice_runoff           = ',qflx_ice_runoff_grc(indexg)*dtime
             write(iulog,*)'qflx_snwcp_discarded_ice  = ',qflx_snwcp_discarded_ice_grc(indexg)*dtime
             write(iulog,*)'qflx_snwcp_discarded_liq  = ',qflx_snwcp_discarded_liq_grc(indexg)*dtime
             write(iulog,*)'deltawb                   = ',endwb_grc(indexg)-begwb_grc(indexg)
             write(iulog,*)'deltawb/dtime             = ',(endwb_grc(indexg)-begwb_grc(indexg))/dtime
             write(iulog,*)'qflx_drain_perched        = ',qflx_drain_perched_grc(indexg)*dtime
             write(iulog,*)'forc_flood                = ',forc_flood_grc(indexg)*dtime
             write(iulog,*)'qflx_glcice_dyn_water_flux = ',qflx_glcice_dyn_water_flux_grc(indexg)*dtime

             write(iulog,*)'CTSM is stopping'
             call endrun(subgrid_index=indexg, subgrid_level=subgrid_level_gridcell, msg=errmsg(sourcefile, __LINE__))
          end if

       end if

       ! Snow balance check at the column level.

       call waterstate_inst%CalculateTotalH2osno(bounds, num_allc, filter_allc, &
            caller = 'BalanceCheck', &
            h2osno_total = h2osno_total(bounds%begc:bounds%endc))

       do c = bounds%begc,bounds%endc
          if (col%active(c)) then
             g = col%gridcell(c)
             l = col%landunit(c)

             ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
             ! any given time step but only if there is at least one snow layer.  h2osno 
             ! also includes snow that is part of the soil column (an initial snow layer is 
             ! only created if h2osno > 10mm).

             if (col%snl(c) < 0) then
                snow_sources(c) = qflx_prec_grnd(c) + qflx_soliddew_to_top_layer(c) &
                     + qflx_liqdew_to_top_layer(c)
                snow_sinks(c)  = qflx_solidevap_from_top_layer(c) + qflx_liqevap_from_top_layer(c) &
                     + qflx_snow_drain(c) + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) &
                     + qflx_snwcp_discarded_ice_col(c) + qflx_snwcp_discarded_liq_col(c) &
                     + qflx_sl_top_soil(c)

                if (lun%itype(l) == istdlak) then 
                   snow_sources(c) = qflx_snow_grnd_col(c) &
                        + frac_sno_eff(c) * (qflx_liq_grnd_col(c) &
                        +  qflx_soliddew_to_top_layer(c) + qflx_liqdew_to_top_layer(c) ) 
                   snow_sinks(c)   = frac_sno_eff(c) * (qflx_solidevap_from_top_layer(c) &
                        + qflx_liqevap_from_top_layer(c) ) + qflx_snwcp_ice(c) + qflx_snwcp_liq(c)  &
                        + qflx_snwcp_discarded_ice_col(c) + qflx_snwcp_discarded_liq_col(c)  &
                        + qflx_snow_drain(c)  + qflx_sl_top_soil(c)
                endif

                 if (col%itype(c) == icol_road_perv .or. lun%itype(l) == istsoil .or. &
                      lun%itype(l) == istcrop .or. lun%itype(l) == istwet .or. &
                      lun%itype(l) == istice) then
                   snow_sources(c) = (qflx_snow_grnd_col(c) - qflx_snow_h2osfc(c) ) &
                          + frac_sno_eff(c) * (qflx_liq_grnd_col(c) &
                          + qflx_soliddew_to_top_layer(c) + qflx_liqdew_to_top_layer(c) ) &
                          + qflx_h2osfc_to_ice(c)
                   snow_sinks(c) = frac_sno_eff(c) * (qflx_solidevap_from_top_layer(c) &
                          + qflx_liqevap_from_top_layer(c)) + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) &
                          + qflx_snwcp_discarded_ice_col(c) + qflx_snwcp_discarded_liq_col(c) &
                          + qflx_snow_drain(c) + qflx_sl_top_soil(c)
                endif

                errh2osno(c) = (h2osno_total(c) - h2osno_old(c)) - (snow_sources(c) - snow_sinks(c)) * dtime
             else
                snow_sources(c) = 0._r8
                snow_sinks(c) = 0._r8
                errh2osno(c) = 0._r8
             end if
          else
             errh2osno(c) = 0._r8
          end if
       end do


       errh2osno_max_val = maxval( abs(errh2osno(bounds%begc:bounds%endc)))
       
       if (errh2osno_max_val > h2o_warning_thresh) then
            indexc = maxloc( abs(errh2osno(bounds%begc:bounds%endc)), 1) + bounds%begc -1
            global_index = get_global_index(subgrid_index=indexc, subgrid_level=subgrid_level_column)
            write(iulog,*)'WARNING:  snow balance error '
            write(iulog,*)'nstep= ',nstep, &
                 ' local indexc= ',indexc, &
                 ' global indexc= ',global_index, &
                 ' col%itype= ',col%itype(indexc), &
                 ' lun%itype= ',lun%itype(col%landunit(indexc)), &
                 ' errh2osno= ',errh2osno(indexc)

            if ((errh2osno_max_val > error_thresh) .and. (DAnstep > skip_steps) ) then
                 write(iulog,*)'CTSM is stopping because errh2osno > ', error_thresh, ' mm'
                 write(iulog,*)'nstep              = ',nstep
                 write(iulog,*)'errh2osno          = ',errh2osno(indexc)
                 write(iulog,*)'snl                = ',col%snl(indexc)
                 write(iulog,*)'snow_depth         = ',snow_depth(indexc)
                 write(iulog,*)'frac_sno_eff       = ',frac_sno_eff(indexc)
                 write(iulog,*)'h2osno             = ',h2osno_total(indexc)
                 write(iulog,*)'h2osno_old         = ',h2osno_old(indexc)
                 write(iulog,*)'snow_sources       = ',snow_sources(indexc)*dtime
                 write(iulog,*)'snow_sinks         = ',snow_sinks(indexc)*dtime
                 write(iulog,*)'qflx_prec_grnd     = ',qflx_prec_grnd(indexc)*dtime
                 write(iulog,*)'qflx_snow_grnd_col = ',qflx_snow_grnd_col(indexc)*dtime
                 write(iulog,*)'qflx_liq_grnd_col  = ',qflx_liq_grnd_col(indexc)*dtime
                 write(iulog,*)'qflx_solidevap_from_top_layer = ',qflx_solidevap_from_top_layer(indexc)*dtime
                 write(iulog,*)'qflx_snow_drain    = ',qflx_snow_drain(indexc)*dtime
                 write(iulog,*)'qflx_liqevap_from_top_layer = ',qflx_liqevap_from_top_layer(indexc)*dtime
                 write(iulog,*)'qflx_soliddew_to_top_layer  = ',qflx_soliddew_to_top_layer(indexc)*dtime
                 write(iulog,*)'qflx_liqdew_to_top_layer    = ',qflx_liqdew_to_top_layer(indexc)*dtime
                 write(iulog,*)'qflx_snwcp_ice     = ',qflx_snwcp_ice(indexc)*dtime
                 write(iulog,*)'qflx_snwcp_liq     = ',qflx_snwcp_liq(indexc)*dtime
                 write(iulog,*)'qflx_snwcp_discarded_ice = ',qflx_snwcp_discarded_ice_col(indexc)*dtime
                 write(iulog,*)'qflx_snwcp_discarded_liq = ',qflx_snwcp_discarded_liq_col(indexc)*dtime
                 write(iulog,*)'qflx_sl_top_soil   = ',qflx_sl_top_soil(indexc)*dtime
                 write(iulog,*)'CTSM is stopping'
                 call endrun(subgrid_index=indexc, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
            end if

       end if

       ! Energy balance checks

       do p = bounds%begp, bounds%endp
          if (patch%active(p)) then
             c = patch%column(p)
             l = patch%landunit(p)
             g = patch%gridcell(p)

             ! Solar radiation energy balance
             ! Do not do this check for an urban patch since it will not balance on a per-column
             ! level because of interactions between columns and since a separate check is done
             ! in the urban radiation module
             if (.not. lun%urbpoi(l)) then
                   errsol(p) = fsa(p) + fsr(p) &
                        - (forc_solad_col(c,1) + forc_solad_col(c,2) + forc_solai(g,1) + forc_solai(g,2))
             else
                errsol(p) = spval
             end if

             ! Longwave radiation energy balance
             ! Do not do this check for an urban patch since it will not balance on a per-column
             ! level because of interactions between columns and since a separate check is done
             ! in the urban radiation module
             if (.not. lun%urbpoi(l)) then
                errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(c)
             else
                errlon(p) = spval
             end if

             ! Surface energy balance
             ! Changed to using (eflx_lwrad_net) here instead of (forc_lwrad - eflx_lwrad_out) because
             ! there are longwave interactions between urban columns (and therefore patches). 
             ! For surfaces other than urban, (eflx_lwrad_net) equals (forc_lwrad - eflx_lwrad_out),
             ! and a separate check is done above for these terms.

             if (.not. lun%urbpoi(l)) then
                errseb(p) = sabv(p) + sabg_chk(p) + forc_lwrad(c) - eflx_lwrad_out(p) &
                     - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p) - dhsdt_canopy(p)
             else
                errseb(p) = sabv(p) + sabg(p) &
                     - eflx_lwrad_net(p) &
                     - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p) &
                     + eflx_wasteheat_patch(p) + eflx_heat_from_ac_patch(p) + eflx_traffic_patch(p) &
                     + eflx_ventilation_patch(p)
             end if
             !TODO MV - move this calculation to a better place - does not belong in BalanceCheck 
             netrad(p) = fsa(p) - eflx_lwrad_net(p)
          else
             errsol(p) = 0._r8
             errlon(p) = 0._r8
             errseb(p) = 0._r8
          end if
       end do

       ! Solar radiation energy balance check

       errsol_max_val = maxval( abs(errsol(bounds%begp:bounds%endp)), mask = (errsol(bounds%begp:bounds%endp) /= spval) ) 

       if  ((errsol_max_val > energy_warning_thresh) .and. (DAnstep > skip_steps)) then

           indexp = maxloc( abs(errsol(bounds%begp:bounds%endp)), 1 , mask = (errsol(bounds%begp:bounds%endp) /= spval) ) + bounds%begp -1
           indexg = patch%gridcell(indexp)
           indexc = patch%column(indexp)
           write(iulog,*)'WARNING:: BalanceCheck, solar radiation balance error (W/m2)'
           write(iulog,*)'nstep         = ',nstep
           write(iulog,*)'errsol        = ',errsol(indexp)
           
           if (errsol_max_val > error_thresh) then
               write(iulog,*)'CTSM is stopping because errsol > ', error_thresh, ' W/m2'
               write(iulog,*)'fsa           = ',fsa(indexp)
               write(iulog,*)'fsr           = ',fsr(indexp)
               write(iulog,*)'forc_solad(1) = ',forc_solad(indexg,1)
               write(iulog,*)'forc_solad(2) = ',forc_solad(indexg,2)
               write(iulog,*)'forc_solai(1) = ',forc_solai(indexg,1)
               write(iulog,*)'forc_solai(2) = ',forc_solai(indexg,2)
               write(iulog,*)'forc_tot      = ',forc_solad(indexg,1)+forc_solad(indexg,2) &
                    +forc_solai(indexg,1)+forc_solai(indexg,2)

               write(iulog,*)'coszen_col:',surfalb_inst%coszen_col(indexc)
               write(iulog,*)'fabd:',fabd(indexp,:)
               write(iulog,*)'fabi:',fabi(indexp,:)
               write(iulog,*)'albd:',albd(indexp,:)
               write(iulog,*)'albi:',albi(indexp,:)
               write(iulog,*)'ftdd:',ftdd(indexp,:)
               write(iulog,*)'ftid:',ftid(indexp,:)
               write(iulog,*)'ftii:',ftii(indexp,:)

               
               write(iulog,*)'CTSM is stopping'
               call endrun(subgrid_index=indexp, subgrid_level=subgrid_level_patch, msg=errmsg(sourcefile, __LINE__))
           end if

       end if
       
       ! Longwave radiation energy balance check

       errlon_max_val = maxval( abs(errlon(bounds%begp:bounds%endp)), mask = (errlon(bounds%begp:bounds%endp) /= spval) )

       if ((errlon_max_val > energy_warning_thresh) .and. (DAnstep > skip_steps)) then
            indexp = maxloc( abs(errlon(bounds%begp:bounds%endp)), 1 , mask = (errlon(bounds%begp:bounds%endp) /= spval) ) + bounds%begp -1
            write(iulog,*)'indexp         = ',indexp
            write(iulog,*)'WARNING: BalanceCheck: longwave energy balance error (W/m2)'
            write(iulog,*)'nstep        = ',nstep
            write(iulog,*)'errlon       = ',errlon(indexp)
            if (errlon_max_val > error_thresh ) then
                 write(iulog,*)'CTSM is stopping because errlon > ', error_thresh, ' W/m2'
                 call endrun(subgrid_index=indexp, subgrid_level=subgrid_level_patch, msg=errmsg(sourcefile, __LINE__))
            end if
       end if

       ! Surface energy balance check

       errseb_max_val = maxval( abs(errseb(bounds%begp:bounds%endp)))

       if ((errseb_max_val > energy_warning_thresh) .and. (DAnstep > skip_steps)) then

           indexp = maxloc( abs(errseb(bounds%begp:bounds%endp)), 1 ) + bounds%begp -1
           indexc = patch%column(indexp)
           indexg = patch%gridcell(indexp)

           write(iulog,*)'WARNING: BalanceCheck: surface flux energy balance error (W/m2)'
           write(iulog,*)'nstep          = ' ,nstep
           write(iulog,*)'errseb         = ' ,errseb(indexp)

           if ( errseb_max_val > error_thresh ) then
              write(iulog,*)'CTSM is stopping because errseb > ', error_thresh, ' W/m2'
              write(iulog,*)'sabv           = ' ,sabv(indexp)
              write(iulog,*)'sabg           = ' ,sabg(indexp), ((1._r8- frac_sno(indexc))*sabg_soil(indexp) + &
                   frac_sno(indexc)*sabg_snow(indexp)),sabg_chk(indexp)
              write(iulog,*)'forc_tot      = '  ,forc_solad(indexg,1) + forc_solad(indexg,2) + &
                   forc_solai(indexg,1) + forc_solai(indexg,2)

              write(iulog,*)'eflx_lwrad_net = ' ,eflx_lwrad_net(indexp)
              write(iulog,*)'eflx_sh_tot    = ' ,eflx_sh_tot(indexp)
              write(iulog,*)'eflx_lh_tot    = ' ,eflx_lh_tot(indexp)
              write(iulog,*)'eflx_soil_grnd = ' ,eflx_soil_grnd(indexp)
              write(iulog,*)'dhsdt_canopy   = ' ,dhsdt_canopy(indexp)
              write(iulog,*)'fsa fsr = '        ,fsa(indexp),    fsr(indexp)
              write(iulog,*)'fabd fabi = '      ,fabd(indexp,:), fabi(indexp,:)
              write(iulog,*)'albd albi = '      ,albd(indexp,:), albi(indexp,:)
              write(iulog,*)'ftii ftdd ftid = ' ,ftii(indexp,:), ftdd(indexp,:),ftid(indexp,:)
              write(iulog,*)'elai esai = '      ,elai(indexp),   esai(indexp)
              write(iulog,*)'CTSM is stopping'
              call endrun(subgrid_index=indexp, subgrid_level=subgrid_level_patch, msg=errmsg(sourcefile, __LINE__))
           end if

       end if

       ! Soil energy balance check

       errsoi_col_max_val  =  maxval( abs(errsoi_col(bounds%begc:bounds%endc)) , mask = col%active(bounds%begc:bounds%endc))

       if (errsoi_col_max_val > 1.0e-5_r8 ) then
           indexc = maxloc( abs(errsoi_col(bounds%begc:bounds%endc)), 1 , mask = col%active(bounds%begc:bounds%endc) ) + bounds%begc -1
           write(iulog,*)'WARNING: BalanceCheck: soil balance error (W/m2)'
           write(iulog,*)'nstep         = ',nstep
           write(iulog,*)'errsoi_col    = ',errsoi_col(indexc)

           if ((errsoi_col_max_val > 1.e-4_r8) .and. (DAnstep > skip_steps)) then
              write(iulog,*)'CTSM is stopping'
              call endrun(subgrid_index=indexc, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
           end if
       end if 
       
     end associate

   end subroutine BalanceCheck

end module BalanceCheckMod
