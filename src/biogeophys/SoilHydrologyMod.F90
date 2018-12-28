module SoilHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate soil hydrology
  !
#include "shr_assert.h"
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use decompMod         , only : bounds_type
  use clm_varctl        , only : iulog, use_vichydro
  use clm_varcon        , only : e_ice, denh2o, denice, rpi
  use clm_varcon        , only : pondmx_urban
  use clm_varpar        , only : nlevsoi, nlevgrnd, nlayer, nlayert
  use column_varcon     , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon     , only : icol_road_imperv
  use landunit_varcon   , only : istsoil, istcrop
  use clm_time_manager  , only : get_step_size
  use NumericsMod       , only : truncate_small_values
  use EnergyFluxType    , only : energyflux_type
  use InfiltrationExcessRunoffMod, only : infiltration_excess_runoff_type
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use WaterFluxType     , only : waterflux_type
  use WaterFluxBulkType , only : waterfluxbulk_type
  use WaterStateType    , only : waterstate_type
  use WaterStateBulkType, only : waterstatebulk_type
  use WaterDiagnosticBulkType, only : waterdiagnosticbulk_type
  use TemperatureType   , only : temperature_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : patch                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilHydReadNML       ! Read in the Soil hydrology namelist
  public :: SetSoilWaterFractions ! Set diagnostic variables related to the fraction of water and ice in each layer
  public :: SetQflxInputs        ! Set the flux of water into the soil from the top
  public :: UpdateH2osfc         ! Calculate fluxes out of h2osfc and update the h2osfc state
  public :: Infiltration         ! Calculate total infiltration
  public :: TotalSurfaceRunoff   ! Calculate total surface runoff
  public :: UpdateUrbanPonding   ! Update the state variable representing ponding on urban surfaces
  public :: WaterTable           ! Calculate water table before imposing drainage
  public :: Drainage             ! Calculate subsurface drainage
  public :: CLMVICMap
  public :: PerchedWaterTable    ! Calculate perched water table
  public :: PerchedLateralFlow   ! Calculate lateral flow from perched saturated zone
  public :: ThetaBasedWaterTable ! Calculate water table from soil moisture state
  public :: LateralFlowPowerLaw  ! Calculate lateral flow based on power law drainage function
  public :: RenewCondensation    ! Misc. corrections
  public :: CalcIrrigWithdrawals ! Calculate irrigation withdrawals from groundwater by layer
  public :: WithdrawGroundwaterIrrigation   ! Remove groundwater irrigation from unconfined and confined aquifers
  
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: QflxH2osfcSurf      ! Compute qflx_h2osfc_surf
  private :: QflxH2osfcDrain     ! Compute qflx_h2osfc_drain

  !-----------------------------------------------------------------------
  real(r8), private :: baseflow_scalar = 1.e-2_r8

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine soilHydReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for soil hydrology
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'soilHydReadNML'
    character(len=*), parameter :: nmlname = 'soilhydrology_inparm'
    !-----------------------------------------------------------------------
    namelist /soilhydrology_inparm/ baseflow_scalar

    ! Initialize options to default values, in case they are not specified in
    ! the namelist


    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=soilhydrology_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (baseflow_scalar, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=soilhydrology_inparm)
       write(iulog,*) ' '
    end if

  end subroutine soilhydReadNML

  !-----------------------------------------------------------------------
  subroutine SetSoilWaterFractions(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_inst, soilstate_inst, waterstatebulk_inst)
    !
    ! !DESCRIPTION:
    ! Set diagnostic variables related to the fraction of water and ice in each layer
    !
    ! !USES:
    use clm_varcon, only : denice
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds               
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(inout) :: soilstate_inst
    type(waterstatebulk_type), intent(in)    :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: j, fc, c
    real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevsoi) !partial volume of ice lens in layer
    real(r8) :: icefrac_orig ! original formulation for icefrac

    character(len=*), parameter :: subname = 'SetSoilWaterFractions'
    !-----------------------------------------------------------------------

    associate( &
         dz               =>    col%dz                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 

         watsat           =>    soilstate_inst%watsat_col           , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         eff_porosity     =>    soilstate_inst%eff_porosity_col     , & ! Output: [real(r8) (:,:) ]  effective porosity = porosity - vol_ice

         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col      , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col      , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)

         origflag         =>    soilhydrology_inst%origflag         , & ! Input:  logical
         icefrac          =>    soilhydrology_inst%icefrac_col      , & ! Output: [real(r8) (:,:) ]                                                  
         fracice          =>    soilhydrology_inst%fracice_col        & ! Output: [real(r8) (:,:) ]  fractional impermeability (-)                    
         )

    do j = 1,nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Porosity of soil, partial volume of ice and liquid, fraction of ice in each layer,
          ! fractional impermeability
          vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
          eff_porosity(c,j) = max(0.01_r8,watsat(c,j)-vol_ice(c,j))
          icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))

          ! fracice is only used in code with origflag == 1. For this calculation, we use
          ! the version of icefrac that was used in this original hydrology code.
          if (h2osoi_ice(c,j) == 0._r8) then
             ! Avoid possible divide by zero (in case h2osoi_liq(c,j) is also 0)
             icefrac_orig = 0._r8
          else
             icefrac_orig = min(1._r8,h2osoi_ice(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j)))
          end if
          fracice(c,j) = max(0._r8,exp(-3._r8*(1._r8-icefrac_orig))- exp(-3._r8))/(1.0_r8-exp(-3._r8))
       end do
    end do

    end associate

  end subroutine SetSoilWaterFractions

   !-----------------------------------------------------------------------
   subroutine SetQflxInputs(bounds, num_hydrologyc, filter_hydrologyc, &
        waterfluxbulk_inst, waterdiagnosticbulk_inst)
     !
     ! !DESCRIPTION:
     ! Set various input fluxes of water
     !
     ! !ARGUMENTS:
     type(bounds_type)          , intent(in)    :: bounds
     integer                    , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                    , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(waterfluxbulk_type)       , intent(inout) :: waterfluxbulk_inst
     type(waterdiagnosticbulk_type)     , intent(in) :: waterdiagnosticbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer :: fc, c
     real(r8) :: qflx_evap ! evaporation for this column
     real(r8) :: fsno      ! copy of frac_sno

     character(len=*), parameter :: subname = 'SetQflxInputs'
     !-----------------------------------------------------------------------

     associate( &
          snl                     =>    col%snl                                   , & ! Input:  [integer  (:)   ]  minus number of snow layers

          qflx_top_soil           =>    waterfluxbulk_inst%qflx_top_soil_col          , & ! Output: [real(r8) (:)]  net water input into soil from top (mm/s)
          qflx_in_soil            =>    waterfluxbulk_inst%qflx_in_soil_col           , & ! Output: [real(r8) (:)]  surface input to soil (mm/s)
          qflx_top_soil_to_h2osfc =>   waterfluxbulk_inst%qflx_top_soil_to_h2osfc_col , & ! Output: [real(r8) (:)]  portion of qflx_top_soil going to h2osfc, minus evaporation (mm/s)
          qflx_rain_plus_snomelt  =>    waterfluxbulk_inst%qflx_rain_plus_snomelt_col , & ! Input:  [real(r8) (:)]  rain plus snow melt falling on the soil (mm/s)
          qflx_snow_h2osfc        =>    waterfluxbulk_inst%qflx_snow_h2osfc_col       , & ! Input:  [real(r8) (:)]  snow falling on surface water (mm/s)
          qflx_floodc             =>    waterfluxbulk_inst%qflx_floodc_col            , & ! Input:  [real(r8) (:)]  column flux of flood water from RTM
          qflx_ev_soil            =>    waterfluxbulk_inst%qflx_ev_soil_col           , & ! Input:  [real(r8) (:)]  evaporation flux from soil (W/m**2) [+ to atm]
          qflx_evap_grnd          =>    waterfluxbulk_inst%qflx_evap_grnd_col         , & ! Input:  [real(r8) (:)]  ground surface evaporation rate (mm H2O/s) [+]
          qflx_ev_h2osfc          =>    waterfluxbulk_inst%qflx_ev_h2osfc_col         , & ! Input:  [real(r8) (:)]  evaporation flux from h2osfc (W/m**2) [+ to atm]
          qflx_sat_excess_surf    =>    waterfluxbulk_inst%qflx_sat_excess_surf_col   , & ! Input:  [real(r8) (:)]  surface runoff due to saturated surface (mm H2O /s)

          frac_sno                =>    waterdiagnosticbulk_inst%frac_sno_eff_col          , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          frac_h2osfc             =>    waterdiagnosticbulk_inst%frac_h2osfc_col          & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
          )

     do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)

        qflx_top_soil(c) = qflx_rain_plus_snomelt(c) + qflx_snow_h2osfc(c) + qflx_floodc(c)

        ! ------------------------------------------------------------------------
        ! Partition surface inputs between soil and h2osfc
        ! ------------------------------------------------------------------------

        if (snl(c) >= 0) then
           fsno=0._r8
           ! if no snow layers, sublimation is removed from h2osoi_ice in drainage
           qflx_evap=qflx_evap_grnd(c)
        else
           fsno=frac_sno(c)
           qflx_evap=qflx_ev_soil(c)
        endif

        qflx_in_soil(c) = (1._r8 - frac_h2osfc(c)) * (qflx_top_soil(c)  - qflx_sat_excess_surf(c))
        qflx_top_soil_to_h2osfc(c) = frac_h2osfc(c) * (qflx_top_soil(c)  - qflx_sat_excess_surf(c))

        ! remove evaporation (snow treated in SnowHydrology)
        qflx_in_soil(c) = qflx_in_soil(c) - (1.0_r8 - fsno - frac_h2osfc(c))*qflx_evap
        qflx_top_soil_to_h2osfc(c) =  qflx_top_soil_to_h2osfc(c)  - frac_h2osfc(c) * qflx_ev_h2osfc(c)

     end do

     end associate

   end subroutine SetQflxInputs

   !-----------------------------------------------------------------------
   subroutine RouteInfiltrationExcess(bounds, num_hydrologyc, filter_hydrologyc, &
        waterfluxbulk_inst, soilhydrology_inst)
     !
     ! !DESCRIPTION:
     ! Route the infiltration excess runoff flux
     !
     ! This adjusts infiltration, input to h2osfc, and the associated surface runoff flux.
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
     !
     ! !LOCAL VARIABLES:
     integer :: fc, c

     character(len=*), parameter :: subname = 'RouteInfiltrationExcess'
     !-----------------------------------------------------------------------

     associate( &
          qflx_in_soil_limited    => waterfluxbulk_inst%qflx_in_soil_limited_col              , & ! Output: [real(r8) (:)   ] surface input to soil, limited by max infiltration rate (mm H2O /s)
          qflx_in_h2osfc          => waterfluxbulk_inst%qflx_in_h2osfc_col                    , & ! Output: [real(r8) (:)   ] total surface input to h2osfc
          qflx_infl_excess_surf   => waterfluxbulk_inst%qflx_infl_excess_surf_col             , & ! Output: [real(r8) (:)   ] surface runoff due to infiltration excess (mm H2O /s)
          qflx_in_soil            => waterfluxbulk_inst%qflx_in_soil_col                      , & ! Input:  [real(r8) (:)   ] surface input to soil (mm/s)
          qflx_top_soil_to_h2osfc => waterfluxbulk_inst%qflx_top_soil_to_h2osfc_col           , & ! Input:  [real(r8) (:)   ] portion of qflx_top_soil going to h2osfc, minus evaporation (mm/s)
          qflx_infl_excess        => waterfluxbulk_inst%qflx_infl_excess_col                  , & ! Input:  [real(r8) (:)   ] infiltration excess runoff (mm H2O /s)

          h2osfcflag              => soilhydrology_inst%h2osfcflag                          & ! Input:  integer
          )

     do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        if (lun%itype(col%landunit(c)) == istsoil .or. lun%itype(col%landunit(c))==istcrop) then
           qflx_in_soil_limited(c) = qflx_in_soil(c) - qflx_infl_excess(c)
           if (h2osfcflag /= 0) then
              qflx_in_h2osfc(c) =  qflx_top_soil_to_h2osfc(c) + qflx_infl_excess(c)
              qflx_infl_excess_surf(c) = 0._r8
           else
              ! No h2osfc pool, so qflx_infl_excess goes directly to surface runoff
              qflx_in_h2osfc(c) = qflx_top_soil_to_h2osfc(c)
              qflx_infl_excess_surf(c) = qflx_infl_excess(c)
           end if
        else
           ! non-vegetated landunits (i.e. urban) use original CLM4 code
           qflx_in_soil_limited(c) = qflx_in_soil(c)
           qflx_in_h2osfc(c) = 0._r8
           qflx_infl_excess_surf(c) = 0._r8
        end if
     end do

     end associate

   end subroutine RouteInfiltrationExcess

   !-----------------------------------------------------------------------
   subroutine UpdateH2osfc(bounds, num_hydrologyc, filter_hydrologyc, &
        infiltration_excess_runoff_inst, &
        energyflux_inst, soilhydrology_inst, &
        waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculate fluxes out of h2osfc and update the h2osfc state
     !
     ! !USES:
     use shr_const_mod    , only : shr_const_pi
     use clm_varcon       , only : denh2o, denice, roverg, wimp, tfrz
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_sunwall, icol_shadewall, icol_road_perv
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(infiltration_excess_runoff_type), intent(in) :: infiltration_excess_runoff_inst
     type(energyflux_type)    , intent(in)    :: energyflux_inst
     type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(waterdiagnosticbulk_type)    , intent(in) :: waterdiagnosticbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: c,l,fc                                     ! indices
     real(r8) :: dtime                                      ! land model time step (sec)
     real(r8) :: h2osfc_partial(bounds%begc:bounds%endc)    ! partially-updated h2osfc
     !-----------------------------------------------------------------------

     associate(                                                        & 
          qinmax           =>    infiltration_excess_runoff_inst%qinmax_col , & ! Input:  [real(r8) (:)] maximum infiltration rate (mm H2O /s)

          frac_h2osfc      =>    waterdiagnosticbulk_inst%frac_h2osfc_col     , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
          frac_h2osfc_nosnow  => waterdiagnosticbulk_inst%frac_h2osfc_nosnow_col,    & ! Input: [real(r8) (:)   ] col fractional area with surface water greater than zero (if no snow present)
          h2osfc           =>    waterstatebulk_inst%h2osfc_col          , & ! Output: [real(r8) (:)   ]  surface water (mm)                                

          qflx_in_h2osfc   => waterfluxbulk_inst%qflx_in_h2osfc_col      , & ! Input:  [real(r8) (:)   ] total surface input to h2osfc
          qflx_h2osfc_surf =>    waterfluxbulk_inst%qflx_h2osfc_surf_col , & ! Output: [real(r8) (:)   ]  surface water runoff (mm H2O /s)                       
          qflx_h2osfc_drain => waterfluxbulk_inst%qflx_h2osfc_drain_col  , & ! Output: [real(r8) (:)   ]  bottom drainage from h2osfc (mm H2O /s)

          h2osfc_thresh    =>    soilhydrology_inst%h2osfc_thresh_col, & ! Input:  [real(r8) (:)   ]  level at which h2osfc "percolates"                
          h2osfcflag       =>    soilhydrology_inst%h2osfcflag         & ! Input:  integer
          )

     dtime = get_step_size()

     call QflxH2osfcSurf(bounds, num_hydrologyc, filter_hydrologyc, &
          h2osfcflag = h2osfcflag, &
          h2osfc = h2osfc(bounds%begc:bounds%endc), &
          h2osfc_thresh = h2osfc_thresh(bounds%begc:bounds%endc), &
          frac_h2osfc_nosnow = frac_h2osfc_nosnow(bounds%begc:bounds%endc), &
          topo_slope = col%topo_slope(bounds%begc:bounds%endc), &
          qflx_h2osfc_surf = qflx_h2osfc_surf(bounds%begc:bounds%endc))
     
     ! Update h2osfc prior to calculating bottom drainage from h2osfc.
     !
     ! This could be removed if we wanted to do a straight forward Euler, and/or set
     ! things up for a more sophisticated solution method. In the latter, rather than
     ! having h2osfc_partial, we'd have some more sophisticated state estimate here
     do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        h2osfc_partial(c) = h2osfc(c) + (qflx_in_h2osfc(c) - qflx_h2osfc_surf(c)) * dtime
     end do

     call truncate_small_values(num_f = num_hydrologyc, filter_f = filter_hydrologyc, &
          lb = bounds%begc, ub = bounds%endc, &
          data_baseline = h2osfc(bounds%begc:bounds%endc), &
          data = h2osfc_partial(bounds%begc:bounds%endc))

     call QflxH2osfcDrain(bounds, num_hydrologyc, filter_hydrologyc, &
          h2osfcflag = h2osfcflag, &
          h2osfc = h2osfc_partial(bounds%begc:bounds%endc), &
          frac_h2osfc = frac_h2osfc(bounds%begc:bounds%endc), &
          qinmax = qinmax(bounds%begc:bounds%endc), &
          qflx_h2osfc_drain = qflx_h2osfc_drain(bounds%begc:bounds%endc))

     ! Update h2osfc based on fluxes
     do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        h2osfc(c) = h2osfc_partial(c) - qflx_h2osfc_drain(c) * dtime
     end do

     ! Due to rounding errors, fluxes that should have brought h2osfc to exactly 0 may
     ! have instead left it slightly less than or slightly greater than 0. Correct for
     ! that here.
     call truncate_small_values(num_f = num_hydrologyc, filter_f = filter_hydrologyc, &
          lb = bounds%begc, ub = bounds%endc, &
          data_baseline = h2osfc_partial(bounds%begc:bounds%endc), &
          data = h2osfc(bounds%begc:bounds%endc))

    end associate

   end subroutine UpdateH2osfc

   !-----------------------------------------------------------------------
   subroutine Infiltration(bounds, num_hydrologyc, filter_hydrologyc, &
        waterfluxbulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculate total infiltration
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer :: fc, c

     character(len=*), parameter :: subname = 'Infiltration'
     !-----------------------------------------------------------------------

     associate( &
          qflx_infl            => waterfluxbulk_inst%qflx_infl_col            , & ! Output: [real(r8) (:)   ] infiltration (mm H2O /s)
          qflx_in_soil_limited => waterfluxbulk_inst%qflx_in_soil_limited_col , & ! Input:  [real(r8) (:)   ] surface input to soil, limited by max infiltration rate (mm H2O /s)
          qflx_h2osfc_drain    => waterfluxbulk_inst%qflx_h2osfc_drain_col      & ! Input:  [real(r8) (:)   ]  bottom drainage from h2osfc (mm H2O /s)
          )

     do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        qflx_infl(c) = qflx_in_soil_limited(c) + qflx_h2osfc_drain(c)
     end do

     end associate

   end subroutine Infiltration

   !-----------------------------------------------------------------------
   subroutine QflxH2osfcSurf(bounds, num_hydrologyc, filter_hydrologyc, &
        h2osfcflag, h2osfc, h2osfc_thresh, frac_h2osfc_nosnow, topo_slope, &
        qflx_h2osfc_surf)
     !
     ! !DESCRIPTION:
     ! Compute qflx_h2osfc_surf
     !
     ! !USES:
     use clm_varcon, only : pc, mu
     !
     ! !ARGUMENTS:
     type(bounds_type) , intent(in)    :: bounds
     integer           , intent(in)    :: num_hydrologyc                     ! number of column soil points in column filter
     integer           , intent(in)    :: filter_hydrologyc(:)               ! column filter for soil points
     integer           , intent(in)    :: h2osfcflag                         ! true => surface water is active
     real(r8)          , intent(in)    :: h2osfc( bounds%begc: )             ! surface water (mm)
     real(r8)          , intent(in)    :: h2osfc_thresh( bounds%begc: )      ! level at which h2osfc "percolates"
     real(r8)          , intent(in)    :: frac_h2osfc_nosnow( bounds%begc: ) ! fractional area with surface water greater than zero (if no snow present)
     real(r8)          , intent(in)    :: topo_slope( bounds%begc: )         ! topographic slope
     real(r8)          , intent(inout) :: qflx_h2osfc_surf( bounds%begc: )   ! surface water runoff (mm H2O /s)
     !
     ! !LOCAL VARIABLES:
     integer  :: fc, c
     real(r8) :: dtime         ! land model time step (sec)
     real(r8) :: frac_infclust ! fraction of submerged area that is connected
     real(r8) :: k_wet         ! linear reservoir coefficient for h2osfc

     character(len=*), parameter :: subname = 'QflxH2osfcSurf'
     !-----------------------------------------------------------------------

     SHR_ASSERT_ALL((ubound(h2osfc) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(h2osfc_thresh) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(frac_h2osfc_nosnow) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(topo_slope) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(qflx_h2osfc_surf) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

     dtime = get_step_size()

     do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)

        if (h2osfcflag==1) then
           if (frac_h2osfc_nosnow(c) <= pc) then
              frac_infclust=0.0_r8
           else
              frac_infclust=(frac_h2osfc_nosnow(c)-pc)**mu
           endif
        endif

        ! limit runoff to value of storage above S(pc)
        if(h2osfc(c) > h2osfc_thresh(c) .and. h2osfcflag/=0) then
           ! spatially variable k_wet
           k_wet=1.0e-4_r8 * sin((rpi/180.) * topo_slope(c))
           qflx_h2osfc_surf(c) = k_wet * frac_infclust * (h2osfc(c) - h2osfc_thresh(c))

           qflx_h2osfc_surf(c)=min(qflx_h2osfc_surf(c),(h2osfc(c) - h2osfc_thresh(c))/dtime)
        else
           qflx_h2osfc_surf(c)= 0._r8
        endif

        ! cutoff lower limit
        if ( qflx_h2osfc_surf(c) < 1.0e-8) then
           qflx_h2osfc_surf(c) = 0._r8
        end if

     end do

   end subroutine QflxH2osfcSurf

   !-----------------------------------------------------------------------
   subroutine QflxH2osfcDrain(bounds, num_hydrologyc, filter_hydrologyc, &
        h2osfcflag, h2osfc, frac_h2osfc, qinmax, &
        qflx_h2osfc_drain)
     !
     ! !DESCRIPTION:
     ! Compute qflx_h2osfc_drain
     !
     ! Note that, if h2osfc is negative, then qflx_h2osfc_drain will be negative - acting
     ! to exactly restore h2osfc to 0.
     !
     ! !ARGUMENTS:
     type(bounds_type) , intent(in)    :: bounds
     integer           , intent(in)    :: num_hydrologyc                     ! number of column soil points in column filter
     integer           , intent(in)    :: filter_hydrologyc(:)               ! column filter for soil points
     integer           , intent(in)    :: h2osfcflag                         ! true => surface water is active
     real(r8)          , intent(in)    :: h2osfc( bounds%begc: )             ! surface water (mm)
     real(r8)          , intent(in)    :: frac_h2osfc( bounds%begc: )        ! fraction of ground covered by surface water (0 to 1)
     real(r8)          , intent(in)    :: qinmax( bounds%begc: )             ! maximum infiltration rate (mm H2O /s)
     real(r8)          , intent(inout) :: qflx_h2osfc_drain( bounds%begc: )  ! bottom drainage from h2osfc (mm H2O /s)
     !
     ! !LOCAL VARIABLES:
     integer :: fc, c
     real(r8) :: dtime         ! land model time step (sec)

     character(len=*), parameter :: subname = 'QflxH2osfcDrain'
     !-----------------------------------------------------------------------

     SHR_ASSERT_ALL((ubound(h2osfc) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(frac_h2osfc) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(qinmax) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(qflx_h2osfc_drain) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
     
     dtime = get_step_size()

     do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)

        if (h2osfc(c) < 0.0) then
           qflx_h2osfc_drain(c) = h2osfc(c)/dtime
        else
           qflx_h2osfc_drain(c)=min(frac_h2osfc(c)*qinmax(c),h2osfc(c)/dtime)
           if(h2osfcflag==0) then
              ! ensure no h2osfc
              qflx_h2osfc_drain(c)= max(0._r8,h2osfc(c)/dtime)
           end if
        end if
     end do

   end subroutine QflxH2osfcDrain


   !-----------------------------------------------------------------------
   subroutine TotalSurfaceRunoff(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc, &
        waterfluxbulk_inst, soilhydrology_inst, &
        waterstatebulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculate total surface runoff
     !
     ! For hydrologically-active columns, this is the sum of some already-computed terms
     !
     ! In addition, computes total surface runoff for non-hydrologically-active urban
     ! columns
     !
     ! !ARGUMENTS:
     type(bounds_type)    , intent(in)    :: bounds
     integer              , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer              , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     integer              , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer              , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     type(waterfluxbulk_type) , intent(inout) :: waterfluxbulk_inst
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(waterstatebulk_type), intent(in)    :: waterstatebulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: fc, c
     real(r8) :: dtime ! land model time step (sec)

     character(len=*), parameter :: subname = 'TotalSurfaceRunoff'
     !-----------------------------------------------------------------------

     associate( &
          snl              =>    col%snl                             , & ! Input:  [integer  (:)   ]  minus number of snow layers                        

          qflx_surf        =>    waterfluxbulk_inst%qflx_surf_col        , & ! Output: [real(r8) (:)   ]  total surface runoff (mm H2O /s)
          qflx_infl_excess_surf => waterfluxbulk_inst%qflx_infl_excess_surf_col, & ! Input:  [real(r8) (:)   ]  surface runoff due to infiltration excess (mm H2O /s)
          qflx_h2osfc_surf =>    waterfluxbulk_inst%qflx_h2osfc_surf_col,  & ! Input:  [real(r8) (:)   ]  surface water runoff (mm H2O /s)
          qflx_rain_plus_snomelt => waterfluxbulk_inst%qflx_rain_plus_snomelt_col , & ! Input: [real(r8) (:)   ] rain plus snow melt falling on the soil (mm/s)
          qflx_evap_grnd   =>    waterfluxbulk_inst%qflx_evap_grnd_col   , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]    
          qflx_floodc      =>    waterfluxbulk_inst%qflx_floodc_col      , & ! Input:  [real(r8) (:)   ]  column flux of flood water from RTM               
          qflx_sat_excess_surf => waterfluxbulk_inst%qflx_sat_excess_surf_col , & ! Input:  [real(r8) (:)   ]  surface runoff due to saturated surface (mm H2O /s)

          xs_urban         =>    soilhydrology_inst%xs_urban_col     , & ! Output: [real(r8) (:)   ]  excess soil water above urban ponding limit

          h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col        & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                            
          )

     dtime = get_step_size()

     ! ------------------------------------------------------------------------
     ! Set qflx_surf for hydrologically-active columns
     ! ------------------------------------------------------------------------

     do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        ! Depending on whether h2osfcflag is 0 or 1, one of qflx_infl_excess or
        ! qflx_h2osfc_surf will always be 0. But it's safe to just add them both.
        qflx_surf(c) = qflx_sat_excess_surf(c) + qflx_infl_excess_surf(c) + qflx_h2osfc_surf(c)
     end do

     ! ------------------------------------------------------------------------
     ! Set qflx_surf for non-hydrologically-active urban columns
     !
     ! Determine water in excess of ponding limit for urban roof and impervious road.
     ! Excess goes to surface runoff. No surface runoff for sunwall and shadewall.
     ! ------------------------------------------------------------------------

     do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
           ! If there are snow layers then all qflx_rain_plus_snomelt goes to surface runoff
           if (snl(c) < 0) then
              qflx_surf(c) = max(0._r8,qflx_rain_plus_snomelt(c))
           else
              ! NOTE(wjs, 2017-07-11) It looks to me like this use of xs_urban and
              ! h2osoi_liq(c,1) are roughly analogous to h2osfc in hydrologically-active
              ! columns. Why not use h2osfc for urban columns, too?
              xs_urban(c) = max(0._r8, &
                   h2osoi_liq(c,1)/dtime + qflx_rain_plus_snomelt(c) - qflx_evap_grnd(c) - &
                   pondmx_urban/dtime)
              qflx_surf(c) = xs_urban(c)
           end if
           ! send flood water flux to runoff for all urban columns
           qflx_surf(c) = qflx_surf(c)  + qflx_floodc(c)
        else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
           qflx_surf(c) = 0._r8
        end if
     end do

     end associate

   end subroutine TotalSurfaceRunoff

   !-----------------------------------------------------------------------
   subroutine UpdateUrbanPonding(bounds, num_urbanc, filter_urbanc, &
        waterstatebulk_inst, soilhydrology_inst, waterfluxbulk_inst)
     !
     ! !DESCRIPTION:
     ! Update the state variable representing ponding on urban surfaces
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
     type(waterfluxbulk_type)     , intent(in)    :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: fc, c
     real(r8) :: dtime ! land model time step (sec)

     character(len=*), parameter :: subname = 'UpdateUrbanPonding'
     !-----------------------------------------------------------------------

     associate( &
         snl              =>    col%snl                             , & ! Input:  [integer  (:)   ]  minus number of snow layers                        

         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col      , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)

         xs_urban         =>    soilhydrology_inst%xs_urban_col     , & ! Input:  [real(r8) (:)   ]  excess soil water above urban ponding limit

         qflx_rain_plus_snomelt => waterfluxbulk_inst%qflx_rain_plus_snomelt_col , & ! Input: [real(r8) (:)   ] rain plus snow melt falling on the soil (mm/s)
         qflx_evap_grnd   =>    waterfluxbulk_inst%qflx_evap_grnd_col     & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]    
         )

     dtime = get_step_size()

     do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        
        if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
           if (snl(c) >= 0) then
              if (xs_urban(c) > 0.) then
                 h2osoi_liq(c,1) = pondmx_urban
              else
                 h2osoi_liq(c,1) = max(0._r8,h2osoi_liq(c,1)+ &
                      (qflx_rain_plus_snomelt(c)-qflx_evap_grnd(c))*dtime)
              end if
           end if
        end if
     end do

     end associate

   end subroutine UpdateUrbanPonding

   !-----------------------------------------------------------------------
   subroutine WaterTable(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
        soilhydrology_inst, soilstate_inst, temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst) 
     !
     ! !DESCRIPTION:
     ! Calculate watertable, considering aquifer recharge but no drainage.
     !
     ! !USES:
     use clm_varcon       , only : pondmx, tfrz, watmin,denice,denh2o
     use column_varcon    , only : icol_roof, icol_road_imperv
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds  
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(temperature_type)   , intent(in)    :: temperature_inst
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(waterdiagnosticbulk_type)    , intent(inout) :: waterdiagnosticbulk_inst
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: ka                                      ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: qcharge_tot
     real(r8) :: qcharge_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: dflag=0._r8
     !-----------------------------------------------------------------------

     associate(                                                            & 
          snl                =>    col%snl                               , & ! Input:  [integer  (:)   ]  number of snow layers                              
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           

          t_soisno           =>    temperature_inst%t_soisno_col         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       

          h2osfc             =>    waterstatebulk_inst%h2osfc_col            , & ! Input:  [real(r8) (:)   ]  surface water (mm)                                
          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
          h2osoi_vol         =>    waterstatebulk_inst%h2osoi_vol_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
          frac_h2osfc        =>    waterdiagnosticbulk_inst%frac_h2osfc_col       , & ! Input:  [real(r8) (:)   ]                                                    

          qflx_dew_grnd      =>    waterfluxbulk_inst%qflx_dew_grnd_col      , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]      
          qflx_dew_snow      =>    waterfluxbulk_inst%qflx_dew_snow_col      , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]    

          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
          eff_porosity       =>    soilstate_inst%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         

          zwt                =>    soilhydrology_inst%zwt_col            , & ! Output: [real(r8) (:)   ]  water table depth (m)                             
          zwt_perched        =>    soilhydrology_inst%zwt_perched_col    , & ! Output: [real(r8) (:)   ]  perched water table depth (m)                     
          frost_table        =>    soilhydrology_inst%frost_table_col    , & ! Output: [real(r8) (:)   ]  frost table depth (m)                             
          wa                 =>    waterstatebulk_inst%wa_col             , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
          qcharge            =>    soilhydrology_inst%qcharge_col        , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)                      
          origflag           =>    soilhydrology_inst%origflag           , & ! Input:  logical  
          
          qflx_sub_snow      =>    waterfluxbulk_inst%qflx_sub_snow_col      , & ! Output: [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]   
          qflx_drain         =>    waterfluxbulk_inst%qflx_drain_col         , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
          qflx_drain_perched =>    waterfluxbulk_inst%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ]  perched wt sub-surface runoff (mm H2O /s)         
          qflx_rsub_sat      =>    waterfluxbulk_inst%qflx_rsub_sat_col        & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]                 
          )

       ! Get time step

       dtime = get_step_size()

       ! Convert layer thicknesses from m to mm

       do j = 1,nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             dzmm(c,j) = dz(c,j)*1.e3_r8
          end do
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          qflx_drain(c)    = 0._r8
          qflx_rsub_sat(c) = 0._r8
          qflx_drain_perched(c)  = 0._r8       
       end do

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       end do

       !============================== QCHARGE =========================================
       ! Water table changes due to qcharge
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! use analytical expression for aquifer specific yield
          rous = watsat(c,nlevsoi) &
               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
          rous=max(rous,0.02_r8)

          !--  water table is below the soil column  --------------------------------------
          if(jwt(c) == nlevsoi) then             
             wa(c)  = wa(c) + qcharge(c)  * dtime
             zwt(c) = zwt(c) - (qcharge(c)  * dtime)/1000._r8/rous
          else                                
             !-- water table within soil layers 1-9  -------------------------------------
             ! try to raise water table to account for qcharge
             qcharge_tot = qcharge(c) * dtime
             if(qcharge_tot > 0.) then !rising water table
                do j = jwt(c)+1, 1,-1
                   ! use analytical expression for specific yield
                   s_y = watsat(c,j) &
                        * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                   s_y=max(s_y,0.02_r8)

                   qcharge_layer=min(qcharge_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                   qcharge_layer=max(qcharge_layer,0._r8)

                   if(s_y > 0._r8) zwt(c) = zwt(c) - qcharge_layer/s_y/1000._r8

                   qcharge_tot = qcharge_tot - qcharge_layer
                   if (qcharge_tot <= 0.) exit
                enddo
             else ! deepening water table (negative qcharge)
                do j = jwt(c)+1, nlevsoi
                   ! use analytical expression for specific yield
                   s_y = watsat(c,j) &
                        * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                   s_y=max(s_y,0.02_r8)

                   qcharge_layer=max(qcharge_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                   qcharge_layer=min(qcharge_layer,0._r8)
                   qcharge_tot = qcharge_tot - qcharge_layer
                   if (qcharge_tot >= 0.) then 
                      zwt(c) = zwt(c) - qcharge_layer/s_y/1000._r8
                      exit
                   else
                      zwt(c) = zi(c,j)
                   endif

                enddo
                if (qcharge_tot > 0.) zwt(c) = zwt(c) - qcharge_tot/1000._r8/rous
             endif

             !-- recompute jwt for following calculations  ---------------------------------
             ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevsoi
             do j = 1,nlevsoi
                if(zwt(c) <= zi(c,j)) then
                   jwt(c) = j-1
                   exit
                end if
             enddo
          endif
       enddo


       !==  BASEFLOW ==================================================
       ! perched water table code
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then 
             k_frz=nlevsoi
          else
             k_frz=1
          endif

          do k=2, nlevsoi
             if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
                k_frz=k
                exit
             endif
          enddo

          frost_table(c)=z(c,k_frz)

          ! initialize perched water table to frost table, and qflx_drain_perched(c) to zero
          zwt_perched(c)=frost_table(c)

          !===================  water table above frost table  =============================
          ! if water table is above frost table, do not use topmodel baseflow formulation
          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then
          else 
             !===================  water table below frost table  =============================
             !--  compute possible perched water table *and* groundwater table afterwards
             ! locate perched water table from bottom up starting at frost table
             ! sat_lev is an arbitrary saturation level used to determine perched water table
             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) then 
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             if (k_frz > k_perch) then
                ! interpolate between k_perch and k_perch+1 to find perched water table height
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

             endif !k_frz > k_perch 
          endif
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Renew the ice and liquid mass due to condensation

          if (snl(c)+1 >= 1) then

             ! make consistent with how evap_grnd removed in infiltration
             h2osoi_liq(c,1) = h2osoi_liq(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_grnd(c) * dtime
             h2osoi_ice(c,1) = h2osoi_ice(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_snow(c) * dtime
             if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
                qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
                h2osoi_ice(c,1) = 0._r8
             else
                h2osoi_ice(c,1) = h2osoi_ice(c,1) - (1._r8 - frac_h2osfc(c)) * qflx_sub_snow(c) * dtime
             end if
          end if
       end do


       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          ! Renew the ice and liquid mass due to condensation for urban roof and impervious road

          if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
             if (snl(c)+1 >= 1) then
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + qflx_dew_grnd(c) * dtime
                h2osoi_ice(c,1) = h2osoi_ice(c,1) + (qflx_dew_snow(c) * dtime)
                if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
                   qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
                   h2osoi_ice(c,1) = 0._r8
                else
                   h2osoi_ice(c,1) = h2osoi_ice(c,1) - (qflx_sub_snow(c) * dtime)
                end if
             end if
          end if

       end do

     end associate

   end subroutine WaterTable

   !-----------------------------------------------------------------------
   subroutine Drainage(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc,  &
        temperature_inst, soilhydrology_inst, soilstate_inst, waterstatebulk_inst, waterfluxbulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage
     !
     ! !USES:
     use clm_varcon       , only : pondmx, tfrz, watmin,rpi, secspday, nlvic
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(temperature_type)   , intent(in)    :: temperature_inst
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'Drainage'           ! subroutine name
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)            ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: xsia(bounds%begc:bounds%endc)           ! available pore space at layer i (mm)
     real(r8) :: xs1(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer 1 (mm)
     real(r8) :: smpfz(1:nlevsoi)                        ! matric potential of layer right above water table (mm)
     real(r8) :: wtsub                                   ! summation of hk*dzmm for layers below water table (mm**2/s)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: wh_zwt                                  ! water head at the water table depth (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: ka                                      ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: dza                                     ! fff*(zwt-z(jwt)) (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: rsub_top_max
     real(r8) :: h2osoi_vol
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: qcharge_tot
     real(r8) :: qcharge_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: vol_ice
     real(r8) :: dsmax_tmp(bounds%begc:bounds%endc)       ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rsub_tmp                 ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: frac                     ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rel_moist                ! relative moisture, temporary variable
     real(r8) :: wtsub_vic                ! summation of hk*dzmm for layers in the third VIC layer
     !-----------------------------------------------------------------------

     associate(                                                            & 
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          snl                =>    col%snl                               , & ! Input:  [integer  (:)   ] number of snow layers                              

          t_soisno           =>    temperature_inst%t_soisno_col         , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)                       

          h2osfc             =>    waterstatebulk_inst%h2osfc_col            , & ! Input:  [real(r8) (:)   ] surface water (mm)                                

          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  
          eff_porosity       =>    soilstate_inst%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ] effective porosity = porosity - vol_ice         
          hk_l               =>    soilstate_inst%hk_l_col               , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity (mm/s)                    

          depth              =>    soilhydrology_inst%depth_col          , & ! Input:  [real(r8) (:,:) ] VIC soil depth                                   
          c_param            =>    soilhydrology_inst%c_param_col        , & ! Input:  [real(r8) (:)   ] baseflow exponent (Qb)                             
          Dsmax              =>    soilhydrology_inst%dsmax_col          , & ! Input:  [real(r8) (:)   ] max. velocity of baseflow (mm/day)
          max_moist          =>    soilhydrology_inst%max_moist_col      , & ! Input:  [real(r8) (:,:) ] maximum soil moisture (ice + liq)
          moist              =>    soilhydrology_inst%moist_col          , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)                         
          Ds                 =>    soilhydrology_inst%ds_col             , & ! Input:  [real(r8) (:)   ] fracton of Dsmax where non-linear baseflow begins
          Wsvic              =>    soilhydrology_inst%Wsvic_col          , & ! Input:  [real(r8) (:)   ] fraction of maximum soil moisutre where non-liear base flow occurs
          icefrac            =>    soilhydrology_inst%icefrac_col        , & ! Output: [real(r8) (:,:) ] fraction of ice in layer                         
          hkdepth            =>    soilhydrology_inst%hkdepth_col        , & ! Input:  [real(r8) (:)   ] decay factor (m)                                  
          frost_table        =>    soilhydrology_inst%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)                             
          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)                             
          zwt_perched        =>    soilhydrology_inst%zwt_perched_col    , & ! Input:  [real(r8) (:)   ] perched water table depth (m)                     
          aquifer_water_baseline => waterstatebulk_inst%aquifer_water_baseline, & ! Input: [real(r8)] baseline value for water in the unconfined aquifer (wa_col) (mm)
          wa                 =>    waterstatebulk_inst%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)              
          ice                =>    soilhydrology_inst%ice_col            , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)                         
          qcharge            =>    soilhydrology_inst%qcharge_col        , & ! Input:  [real(r8) (:)   ] aquifer recharge rate (mm/s)                      
          origflag           =>    soilhydrology_inst%origflag           , & ! Input:  logical
          h2osfcflag         =>    soilhydrology_inst%h2osfcflag         , & ! Input:  integer
          
          qflx_snwcp_liq     =>    waterfluxbulk_inst%qflx_snwcp_liq_col     , & ! Output: [real(r8) (:)   ] excess liquid h2o due to snow capping (outgoing) (mm H2O /s) [+]
          qflx_ice_runoff_xs =>    waterfluxbulk_inst%qflx_ice_runoff_xs_col , & ! Output: [real(r8) (:)   ] solid runoff from excess ice in soil (mm H2O /s) [+]
          qflx_dew_grnd      =>    waterfluxbulk_inst%qflx_dew_grnd_col      , & ! Output: [real(r8) (:)   ] ground surface dew formation (mm H2O /s) [+]      
          qflx_dew_snow      =>    waterfluxbulk_inst%qflx_dew_snow_col      , & ! Output: [real(r8) (:)   ] surface dew added to snow pack (mm H2O /s) [+]    
          qflx_sub_snow      =>    waterfluxbulk_inst%qflx_sub_snow_col      , & ! Output: [real(r8) (:)   ] sublimation rate from snow pack (mm H2O /s) [+]   
          qflx_drain         =>    waterfluxbulk_inst%qflx_drain_col         , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)                    
          qflx_qrgwl         =>    waterfluxbulk_inst%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ] qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
          qflx_rsub_sat      =>    waterfluxbulk_inst%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ] soil saturation excess [mm h2o/s]                 
          qflx_drain_perched =>    waterfluxbulk_inst%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ] perched wt sub-surface runoff (mm H2O /s)         

          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)                                
          )

       ! Get time step

       dtime = get_step_size()

       ! Convert layer thicknesses from m to mm

       do j = 1,nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             dzmm(c,j) = dz(c,j)*1.e3_r8

             vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             icefrac(c,j) = min(1._r8,vol_ice/watsat(c,j))          
          end do
       end do

       ! Initial set

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          qflx_drain(c)    = 0._r8 
          rsub_bot(c)      = 0._r8
          qflx_rsub_sat(c) = 0._r8
          rsub_top(c)      = 0._r8
          fracice_rsub(c)  = 0._r8

       end do

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       end do

       rous = 0.2_r8

       !==  BASEFLOW ==================================================
       ! perched water table code
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          !  specify maximum drainage rate
          q_perch_max = 1.e-5_r8 * sin(col%topo_slope(c) * (rpi/180._r8))

          ! if layer containing water table is frozen, compute the following:
          !     frost table, perched water table, and drainage from perched saturated layer

          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then 
             k_frz=nlevsoi
          else
             k_frz=1
          endif

          do k=2, nlevsoi
             if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
                k_frz=k
                exit
             endif
          enddo

          frost_table(c)=z(c,k_frz)

          ! initialize perched water table to frost table, and qflx_drain_perched(c) to zero
          zwt_perched(c)=frost_table(c)
          qflx_drain_perched(c) = 0._r8

          !===================  water table above frost table  =============================
          ! if water table is above frost table, do not use topmodel baseflow formulation

          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then
             ! compute drainage from perched saturated region
             wtsub = 0._r8
             q_perch = 0._r8
             do k = jwt(c)+1, k_frz
                imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                wtsub = wtsub + dzmm(c,k)
             end do
             if (wtsub > 0._r8) q_perch = q_perch/wtsub

             qflx_drain_perched(c) = q_perch_max * q_perch &
                  *(frost_table(c) - zwt(c))

             ! remove drainage from perched saturated layers
             rsub_top_tot = -  qflx_drain_perched(c) * dtime
             do k = jwt(c)+1, k_frz
                rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                rsub_top_layer=min(rsub_top_layer,0._r8)
                rsub_top_tot = rsub_top_tot - rsub_top_layer

                h2osoi_liq(c,k) = h2osoi_liq(c,k) + rsub_top_layer

                if (rsub_top_tot >= 0.) then 
                   zwt(c) = zwt(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                   exit
                else
                   zwt(c) = zi(c,k)
                endif
             enddo

             ! if rsub_top_tot is greater than available water (above frost table), 
             !     then decrease qflx_drain_perched by residual amount for water balance
             qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             !-- recompute jwt  ---------------------------------------------------------
             ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevsoi
             do j = 1,nlevsoi
                if(zwt(c) <= zi(c,j)) then
                   jwt(c) = j-1
                   exit
                end if
             enddo
          else 
             !===================  water table below frost table  =============================
             !--  compute possible perched water table *and* groundwater table afterwards
             ! locate perched water table from bottom up starting at frost table
             ! sat_lev is an arbitrary saturation level used to determine perched water table
             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol/watsat(c,k) <= sat_lev) then 
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             if (k_frz > k_perch) then
                ! interpolate between k_perch and k_perch+1 to find perched water table height
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

                ! compute drainage from perched saturated region
                wtsub = 0._r8
                q_perch = 0._r8
                do k = k_perch, k_frz
                   imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                   q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                   wtsub = wtsub + dzmm(c,k)
                end do
                if (wtsub > 0._r8) q_perch = q_perch/wtsub

                qflx_drain_perched(c) = q_perch_max * q_perch &
                     *(frost_table(c) - zwt_perched(c))

                ! no perched water table drainage if using original formulation
                if(origflag == 1) qflx_drain_perched(c) = 0._r8

                ! remove drainage from perched saturated layers
                rsub_top_tot = -  qflx_drain_perched(c) * dtime
                do k = k_perch+1, k_frz
                   rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                   rsub_top_layer=min(rsub_top_layer,0._r8)
                   rsub_top_tot = rsub_top_tot - rsub_top_layer

                   h2osoi_liq(c,k) = h2osoi_liq(c,k) + rsub_top_layer

                   if (rsub_top_tot >= 0.) then 
                      zwt_perched(c) = zwt_perched(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                      exit
                   else
                      zwt_perched(c) = zi(c,k)
                   endif

                enddo

                ! if rsub_top_tot is greater than available water (above frost table), 
                !     then decrease qflx_drain_perched by residual amount for water balance
                qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             else
                qflx_drain_perched(c) = 0._r8
             endif !k_frz > k_perch 

             !-- Topographic runoff  ----------------------------------------------------------------------
             fff(c)         = 1._r8/ hkdepth(c)
             dzsum = 0._r8
             icefracsum = 0._r8
             do j = max(jwt(c),1), nlevsoi
                dzsum  = dzsum + dzmm(c,j)
                icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
             end do
             ! add ice impedance factor to baseflow
             if(origflag == 1) then 
                if (use_vichydro) then
                   call endrun(msg="VICHYDRO is not available for origflag=1"//errmsg(sourcefile, __LINE__))
                else
                   fracice_rsub(c) = max(0._r8,exp(-3._r8*(1._r8-(icefracsum/dzsum))) &
                        - exp(-3._r8))/(1.0_r8-exp(-3._r8))
                   imped=(1._r8 - fracice_rsub(c))
                   rsub_top_max = 5.5e-3_r8
                end if
             else
                if (use_vichydro) then
                   imped=10._r8**(-e_ice*min(1.0_r8,ice(c,nlayer)/max_moist(c,nlayer)))
                   dsmax_tmp(c) = Dsmax(c) * dtime/ secspday !mm/day->mm/dtime
                   rsub_top_max = dsmax_tmp(c)
                else
                   imped=10._r8**(-e_ice*(icefracsum/dzsum))
                   rsub_top_max = 10._r8 * sin((rpi/180.) * col%topo_slope(c))
                end if
             endif
             if (use_vichydro) then
                ! ARNO model for the bottom soil layer (based on bottom soil layer 
                ! moisture from previous time step
                ! use watmin instead for resid_moist to be consistent with default hydrology
                rel_moist = (moist(c,nlayer) - watmin)/(max_moist(c,nlayer)-watmin) 
                frac = (Ds(c) * rsub_top_max )/Wsvic(c)
                rsub_tmp = (frac * rel_moist)/dtime
                if(rel_moist > Wsvic(c))then
                   frac = (rel_moist - Wsvic(c))/(1.0_r8 - Wsvic(c))
                   rsub_tmp = rsub_tmp + (rsub_top_max * (1.0_r8 - Ds(c)/Wsvic(c)) *frac**c_param(c))/dtime
                end if
                rsub_top(c) = imped * rsub_tmp
                ! make sure baseflow isn't negative
                rsub_top(c) = max(0._r8, rsub_top(c))
             else
                rsub_top(c)    = imped * rsub_top_max* exp(-fff(c)*zwt(c))
             end if

             ! use analytical expression for aquifer specific yield
             rous = watsat(c,nlevsoi) &
                  * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
             rous=max(rous,0.02_r8)

             !--  water table is below the soil column  --------------------------------------
             if(jwt(c) == nlevsoi) then             
                wa(c)  = wa(c) - rsub_top(c) * dtime
                zwt(c)     = zwt(c) + (rsub_top(c) * dtime)/1000._r8/rous
                h2osoi_liq(c,nlevsoi) = h2osoi_liq(c,nlevsoi) + max(0._r8,(wa(c)-aquifer_water_baseline))
                wa(c)  = min(wa(c), aquifer_water_baseline)
             else                                
                !-- water table within soil layers 1-9  -------------------------------------
                !============================== RSUB_TOP =========================================
                !--  Now remove water via rsub_top
                rsub_top_tot = - rsub_top(c) * dtime
                !should never be positive... but include for completeness
                if(rsub_top_tot > 0.) then !rising water table

                   call endrun(msg="RSUB_TOP IS POSITIVE in Drainage!"//errmsg(sourcefile, __LINE__))

                else ! deepening water table
                   if (use_vichydro) then
                      wtsub_vic = 0._r8
                      do j = (nlvic(1)+nlvic(2)+1), nlevsoi
                         wtsub_vic = wtsub_vic + hk_l(c,j)*dzmm(c,j)
                      end do

                      do j = (nlvic(1)+nlvic(2)+1), nlevsoi
                         rsub_top_layer=max(rsub_top_tot, rsub_top_tot*hk_l(c,j)*dzmm(c,j)/wtsub_vic)
                         rsub_top_layer=min(rsub_top_layer,0._r8)
                         h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer
                         rsub_top_tot = rsub_top_tot - rsub_top_layer
                      end do
                   else
                      do j = jwt(c)+1, nlevsoi
                         ! use analytical expression for specific yield
                         s_y = watsat(c,j) &
                              * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                         s_y=max(s_y,0.02_r8)

                         rsub_top_layer=max(rsub_top_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                         rsub_top_layer=min(rsub_top_layer,0._r8)
                         h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer

                         rsub_top_tot = rsub_top_tot - rsub_top_layer

                         if (rsub_top_tot >= 0.) then 
                            zwt(c) = zwt(c) - rsub_top_layer/s_y/1000._r8

                            exit
                         else
                            zwt(c) = zi(c,j)
                         endif
                      enddo
                   end if

                   !--  remove residual rsub_top  ---------------------------------------------
                   zwt(c) = zwt(c) - rsub_top_tot/1000._r8/rous
                   wa(c) = wa(c) + rsub_top_tot
                endif

                !-- recompute jwt  ---------------------------------------------------------
                ! allow jwt to equal zero when zwt is in top layer
                jwt(c) = nlevsoi
                do j = 1,nlevsoi
                   if(zwt(c) <= zi(c,j)) then
                      jwt(c) = j-1
                      exit
                   end if
                enddo
             end if! end of jwt if construct

             zwt(c) = max(0.0_r8,zwt(c))
             zwt(c) = min(80._r8,zwt(c))

          endif

       end do

       !  excessive water above saturation added to the above unsaturated layer like a bucket
       !  if column fully saturated, excess water goes to runoff

       do j = nlevsoi,2,-1
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             xsi(c)            = max(h2osoi_liq(c,j)-eff_porosity(c,j)*dzmm(c,j),0._r8)
             h2osoi_liq(c,j)   = min(eff_porosity(c,j)*dzmm(c,j), h2osoi_liq(c,j))
             h2osoi_liq(c,j-1) = h2osoi_liq(c,j-1) + xsi(c)
          end do
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! watmin addition to fix water balance errors
          xs1(c)          = max(max(h2osoi_liq(c,1)-watmin,0._r8)- &
               max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice(c,1)-watmin)),0._r8)
          h2osoi_liq(c,1) = h2osoi_liq(c,1) - xs1(c)

          if (lun%urbpoi(col%landunit(c))) then
             qflx_rsub_sat(c)     = xs1(c) / dtime
          else
             if(h2osfcflag == 1) then
                ! send this water up to h2osfc rather than sending to drainage
                h2osfc(c) = h2osfc(c) + xs1(c)
                qflx_rsub_sat(c)     = 0._r8
             else
                ! use original code to send water to drainage (non-h2osfc case)
                qflx_rsub_sat(c)     = xs1(c) / dtime
             endif
          endif
          ! add in ice check
          xs1(c)          = max(max(h2osoi_ice(c,1),0._r8)-max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1))),0._r8)
          h2osoi_ice(c,1) = min(max(0._r8,pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1)), h2osoi_ice(c,1))
          qflx_ice_runoff_xs(c) = xs1(c) / dtime
       end do

       ! Limit h2osoi_liq to be greater than or equal to watmin.
       ! Get water needed to bring h2osoi_liq equal watmin from lower layer.
       ! If insufficient water in soil layers, get from aquifer water

       do j = 1, nlevsoi-1
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             if (h2osoi_liq(c,j) < watmin) then
                xs(c) = watmin - h2osoi_liq(c,j)
                ! deepen water table if water is passed from below zwt layer
                if(j == jwt(c)) then 
                   zwt(c) = zwt(c) + xs(c)/eff_porosity(c,j)/1000._r8
                endif
             else
                xs(c) = 0._r8
             end if
             h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
             h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
          end do
       end do

       ! Get water for bottom layer from layers above if possible
       j = nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (h2osoi_liq(c,j) < watmin) then
             xs(c) = watmin-h2osoi_liq(c,j)
             searchforwater: do i = nlevsoi-1, 1, -1
                available_h2osoi_liq = max(h2osoi_liq(c,i)-watmin-xs(c),0._r8)
                if (available_h2osoi_liq >= xs(c)) then
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
                   h2osoi_liq(c,i) = h2osoi_liq(c,i) - xs(c)
                   xs(c) = 0._r8
                   exit searchforwater
                else
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + available_h2osoi_liq
                   h2osoi_liq(c,i) = h2osoi_liq(c,i) - available_h2osoi_liq
                   xs(c) = xs(c) - available_h2osoi_liq
                end if
             end do searchforwater
          else
             xs(c) = 0._r8
          end if
          ! Needed in case there is no water to be found
          h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
          ! Instead of removing water from aquifer where it eventually
          ! shows up as excess drainage to the ocean, take it back out of 
          ! drainage
          rsub_top(c) = rsub_top(c) - xs(c)/dtime

       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Sub-surface runoff and drainage

          qflx_drain(c) = qflx_rsub_sat(c) + rsub_top(c)

          ! Set imbalance for snow capping

          qflx_qrgwl(c) = qflx_snwcp_liq(c)

       end do

       ! No drainage for urban columns (except for pervious road as computed above)

       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col%itype(c) /= icol_road_perv) then
             qflx_drain(c) = 0._r8
             ! This must be done for roofs and impervious road (walls will be zero)
             qflx_qrgwl(c) = qflx_snwcp_liq(c)
          end if
       end do

     end associate

   end subroutine Drainage

  !-----------------------------------------------------------------------
  subroutine CLMVICMap(bounds, numf, filter, &
       soilhydrology_inst, waterstatebulk_inst)
     !
     ! !DESCRIPTION:
     ! Performs  the mapping from CLM layers to VIC layers
     ! CLM hydrologically active soil layers are mapped to three VIC layers
     ! by assigning the first nlvic(1) layers to VIC layer 1
     !              the next nlvic(2) layers  to VIC alyer 2
     !              and the remaining to VIC layer 3
     ! mapping from VIC to CLM layers, M.Huang
     !
     ! !USES:
     use clm_varcon  , only : denh2o, denice, watmin
     use decompMod   , only : bounds_type
     !
     ! !REVISION HISTORY:
     ! Created by Maoyi Huang
     ! 11/13/2012, Maoyi Huang: rewrite the mapping modules in CLM4VIC 
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds    
     integer                  , intent(in)    :: numf      ! number of column soil points in column filter
     integer                  , intent(in)    :: filter(:) ! column filter for soil points
     type(waterstatebulk_type)    , intent(in)    :: waterstatebulk_inst 
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     !
     ! !LOCAL VARIABLES
     real(r8) :: ice0(1:nlayer)            ! last step ice lens (mm)  (new)
     real(r8) :: moist0(1:nlayer)          ! last step soil water (mm)  (new)
     integer  :: i, j, c, fc
     ! note: in CLM4 h2osoil_liq unit is kg/m2, in VIC moist is mm
     ! h2osoi_ice is actually water equivalent ice content.
     !-----------------------------------------------------------------------

     associate(                                                   & 
          dz            => col%dz                               , & ! Input:  [real(r8) (:,:)   ] layer depth (m)                        
          zi            => col%zi                               , & ! Input:  [real(r8) (:,:)   ] interface level below a "z" level (m)  
          z             => col%z                                , & ! Input:  [real(r8) (:,:)   ] layer thickness (m)                    

          h2osoi_liq    => waterstatebulk_inst%h2osoi_liq_col       , & ! Input:  [real(r8) (:,:)   ] liquid water (kg/m2)                   
          h2osoi_ice    => waterstatebulk_inst%h2osoi_ice_col       , & ! Input:  [real(r8) (:,:)   ] ice lens (kg/m2)                       
          h2osoi_vol    => waterstatebulk_inst%h2osoi_vol_col       , & ! Input:  [real(r8) (:,:)   ] volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)

          depth         => soilhydrology_inst%depth_col         , & ! Input:  [real(r8) (:,:)   ] layer depth of upper layer (m)         
          porosity      => soilhydrology_inst%porosity_col      , & ! Input:  [real(r8) (:,:)   ] soil porisity (1-bulk_density/soil_density)
          max_moist     => soilhydrology_inst%max_moist_col     , & ! Input:  [real(r8) (:,:)   ] max layer moist + ice (mm)             
          vic_clm_fract => soilhydrology_inst%vic_clm_fract_col , & ! Input:  [real(r8) (:,:,:) ] fraction of VIC layers in each CLM layer
          moist         => soilhydrology_inst%moist_col         , & ! Output: [real(r8) (:,:)   ] liquid water (mm)                      
          ice           => soilhydrology_inst%ice_col           , & ! Output: [real(r8) (:,:)   ] ice lens (mm)                          
          moist_vol     => soilhydrology_inst%moist_vol_col     , & ! Output: [real(r8) (:,:)   ] volumetric soil moisture for VIC soil layers
          top_moist     =>    soilhydrology_inst%top_moist_col    , & ! Output:  [real(r8) (:)   ]  soil moisture in top VIC layers
          top_max_moist =>    soilhydrology_inst%top_max_moist_col, & ! Output:  [real(r8) (:)   ]  maximum soil moisture in top VIC layers
          top_ice       =>    soilhydrology_inst%top_ice_col      , & ! Output:  [real(r8) (:)   ]  ice len in top VIC layers
          top_moist_limited => soilhydrology_inst%top_moist_limited_col  & ! Output:  [real(r8) (:) ]  soil moisture in top layers, limited to no greater than top_max_moist
          )

       ! map CLM to VIC
       do fc = 1, numf
          c = filter(fc)
          do i = 1, nlayer
             ice0(i) = ice(c,i)
             moist0(i) = moist(c,i)
             ice(c,i) = 0._r8
             moist(c,i) = 0._r8
             do j = 1, nlevsoi
                ice(c,i) = ice(c,i) + h2osoi_ice(c,j) * vic_clm_fract(c,i,j)
                moist(c,i) = moist(c,i) + h2osoi_liq(c,j) * vic_clm_fract(c,i,j)
             end do
             ice(c,i)       = min((moist0(i) + ice0(i)), ice(c,i))
             ice(c,i)       = max(0._r8, ice(c,i))
             moist(c,i)     = max(watmin, moist(c,i))
             moist(c,i)     = min(max_moist(c,i)-ice(c,i), moist(c,i))
             moist_vol(c,i) = moist(c,i)/(depth(c,i)*denice) + ice(c,i)/(depth(c,i)*denh2o)
             moist_vol(c,i) = min(porosity(c,i), moist_vol(c,i))
             moist_vol(c,i) = max(0.01_r8, moist_vol(c,i))
          end do

          ! hydrologic inactive layers
          ice(c, nlayer+1:nlayert)       = h2osoi_ice(c, nlevsoi+1:nlevgrnd)
          moist(c, nlayer+1:nlayert)     = h2osoi_liq(c, nlevsoi+1:nlevgrnd)
          moist_vol(c, nlayer+1:nlayert) = h2osoi_vol(c, nlevsoi+1:nlevgrnd)
       end do

       ! Set values related to top VIC layers

       do fc = 1, numf
          c = filter(fc)
          top_moist(c) = 0._r8
          top_ice(c) = 0._r8
          top_max_moist(c) = 0._r8
       end do

       do j = 1, nlayer - 1
          do fc = 1, numf
             c = filter(fc)
             top_ice(c) = top_ice(c) + ice(c,j)
             top_moist(c) =  top_moist(c) + moist(c,j) + ice(c,j)
             top_max_moist(c) = top_max_moist(c) + max_moist(c,j)
          end do
       end do

       do fc = 1, numf
          c = filter(fc)
          top_moist_limited(c) = min(top_moist(c), top_max_moist(c))
       end do

     end associate

   end subroutine CLMVICMap

!#3
   !-----------------------------------------------------------------------
   subroutine PerchedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
        temperature_inst, waterstatebulk_inst, waterfluxbulk_inst) 
     !
     ! !DESCRIPTION:
     ! Calculate watertable, considering aquifer recharge but no drainage.
     !
     ! !USES:
     use clm_varcon       , only : pondmx, tfrz, watmin,denice,denh2o
     use column_varcon    , only : icol_roof, icol_road_imperv
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds  
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(temperature_type)   , intent(in)    :: temperature_inst
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch,k_zwt
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     integer  :: sat_flag
     !-----------------------------------------------------------------------

     associate(                                                            & 
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          t_soisno           =>    temperature_inst%t_soisno_col         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       

          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
          h2osoi_vol         =>    waterstatebulk_inst%h2osoi_vol_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  
          zwt                =>    soilhydrology_inst%zwt_col            , & ! Output: [real(r8) (:)   ]  water table depth (m)                             
          zwt_perched        =>    soilhydrology_inst%zwt_perched_col    , & ! Output: [real(r8) (:)   ]  perched water table depth (m)                     
          frost_table        =>    soilhydrology_inst%frost_table_col    , & ! Output: [real(r8) (:)   ]  frost table depth (m)                             
          origflag           =>    soilhydrology_inst%origflag            & ! Input:  logical  
          )

       ! calculate perched water table location 
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then 
             k_frz=nlevsoi
          else
             k_frz=1
          endif

          do k=2, nlevsoi
             if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
                k_frz=k
                exit
             endif
          enddo

          frost_table(c)=z(c,k_frz)

          ! initialize perched water table to frost table, and qflx_drain_perched(c) to zero
          zwt_perched(c)=frost_table(c)

          !=======  water table above frost table  ===================
          ! if water table is above frost table, do nothing 
          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then
          else 
             !==========  water table below frost table  ============
             ! locate perched water table from bottom up starting at 
             ! frost table sat_lev is an arbitrary saturation level 
             ! used to determine perched water table

             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) then 
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             ! interpolate between k_perch and k_perch+1 to find 
             ! perched water table height
             if (k_frz > k_perch) then
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

             endif !k_frz > k_perch 
          endif
       end do

     end associate

   end subroutine PerchedWaterTable

!#4
   !-----------------------------------------------------------------------
   subroutine PerchedLateralFlow(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
        waterstatebulk_inst, waterfluxbulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage from perched saturated zone
     !
     ! !USES:
     use clm_varcon       , only : pondmx, tfrz, watmin,rpi, secspday, nlvic
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv

     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'PerchedLateralFlow' ! subroutine name
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi) ! layer thickness (mm)
     real(r8) :: wtsub                                   ! summation of hk*dzmm for layers below water table (mm**2/s)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: h2osoi_vol
     real(r8) :: imped
     real(r8) :: drainage_tot
     real(r8) :: drainage_layer
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1, s2, m, b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: vol_ice
     !-----------------------------------------------------------------------

     associate(                                                            & 
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  

          icefrac            =>    soilhydrology_inst%icefrac_col        , & ! Output: [real(r8) (:,:) ] fraction of ice in layer                         
          frost_table        =>    soilhydrology_inst%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)                             
          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)                             
          zwt_perched        =>    soilhydrology_inst%zwt_perched_col    , & ! Input:  [real(r8) (:)   ] perched water table depth (m)                     
          origflag           =>    soilhydrology_inst%origflag           , & ! Input:  logical
          
          qflx_drain_perched =>    waterfluxbulk_inst%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ] perched wt sub-surface runoff (mm H2O /s)         

          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)                                
          )

       ! Get time step

       dtime = get_step_size()

       ! Compute ice fraction in each layer

       do j = 1,nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             dzmm(c,j) = dz(c,j)*1.e3_r8

             vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             icefrac(c,j) = min(1._r8,vol_ice/watsat(c,j))          
          end do
       end do

       ! compute drainage from perched saturated region
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          qflx_drain_perched(c) = 0._r8

          if ((frost_table(c) > zwt_perched(c)) .and. origflag == 0) then

             !  specify maximum drainage rate
             q_perch_max = 1.e-5_r8 * sin(col%topo_slope(c) * (rpi/180._r8))
             
             ! calculate frost table and perched water table locations
             do k=1, nlevsoi
                if (frost_table(c) >= zi(c,k-1) .and. frost_table(c) <= zi(c,k)) then
                   k_frz=k
                   exit
                endif
             enddo
             
             do k=1, nlevsoi
                if (zwt_perched(c) >= zi(c,k-1) .and. zwt_perched(c) <= zi(c,k)) then
                   k_perch=k
                   exit
                endif
             enddo

             wtsub = 0._r8
             q_perch = 0._r8
             do k = k_perch, k_frz
                imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                wtsub = wtsub + dzmm(c,k)
             end do
             if (wtsub > 0._r8) q_perch = q_perch/wtsub
             
             qflx_drain_perched(c) = q_perch_max * q_perch &
                  *(frost_table(c) - zwt_perched(c))
             
             ! no perched water table drainage if using original formulation
             if(origflag == 1) qflx_drain_perched(c) = 0._r8
             
             ! if perched water table exists
             if (k_frz > k_perch) then
                ! remove drainage from perched saturated layers
                drainage_tot = -  qflx_drain_perched(c) * dtime
                do k = k_perch+1, k_frz
                   drainage_layer=max(drainage_tot,-(h2osoi_liq(c,k)-watmin))
                   drainage_layer=min(drainage_layer,0._r8)
                   drainage_tot = drainage_tot - drainage_layer
                   
                   h2osoi_liq(c,k) = h2osoi_liq(c,k) + drainage_layer
                   
                   s_y = watsat(c,k) &
                        * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,k))**(-1./bsw(c,k)))
                   s_y=max(s_y,0.02_r8)
                   if (drainage_tot >= 0.) then 
                      zwt_perched(c) = zwt_perched(c) - drainage_layer/s_y/1000._r8
                      exit
                   else
                      zwt_perched(c) = zi(c,k)
                   endif
                enddo
          
                ! if drainage_tot is greater than available water 
                ! (above frost table), then decrease qflx_drain_perched 
                ! by residual amount for water balance
                qflx_drain_perched(c) = qflx_drain_perched(c) + drainage_tot/dtime          
             else
                qflx_drain_perched(c) = 0._r8
             endif !k_frz > k_perch 
          endif
       enddo

     end associate

   end subroutine PerchedLateralFlow

!#5
   !-----------------------------------------------------------------------
   subroutine ThetaBasedWaterTable(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
        waterstatebulk_inst, waterfluxbulk_inst) 
     !
     ! !DESCRIPTION:
     ! Calculate watertable, considering aquifer recharge but no drainage.
     !
     ! !USES:
     use clm_varcon       , only : denice,denh2o
     use column_varcon    , only : icol_roof, icol_road_imperv
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds  
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(waterstatebulk_type)     , intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type)      , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i                                ! indices
     integer  :: k,k_zwt
     real(r8) :: sat_lev
     real(r8) :: s1,s2,m,b   ! temporary variables used to interpolate theta
     integer  :: sat_flag
     
     !-----------------------------------------------------------------------

     associate(                                                            & 
          nbedrock           =>    col%nbedrock                          , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
          h2osoi_vol         =>    waterstatebulk_inst%h2osoi_vol_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  
          zwt                =>    soilhydrology_inst%zwt_col              & ! Output: [real(r8) (:)   ]  water table depth (m)                             
          )

       ! calculate water table based on soil moisture state
       ! this is a simple search for 1st layer with soil moisture 
       ! less than specified threshold (sat_lev)

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! initialize to depth of bottom of lowest layer
          zwt(c)=zi(c,nlevsoi)

          ! locate water table from bottom up starting at bottom of soil column
          ! sat_lev is an arbitrary saturation level used to determine water table
          sat_lev=0.9
          
          k_zwt=nbedrock(c)
          sat_flag=1 !will remain unchanged if all layers at saturation
           do k=nbedrock(c),1,-1
             h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                  + h2osoi_ice(c,k)/(dz(c,k)*denice)
             
             if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) then 
                k_zwt=k
                sat_flag=0
                exit
             endif
          enddo
          if (sat_flag == 1) k_zwt=1

          ! if soil column above sat_lev, set water table to lower 
          ! interface of first layer
          if (k_zwt == 1) then
             zwt(c)=zi(c,1)
          else if (k_zwt < nbedrock(c)) then
             ! interpolate between k_zwt and k_zwt+1 to find water table height
             s1 = (h2osoi_liq(c,k_zwt)/(dz(c,k_zwt)*denh2o) &
                  + h2osoi_ice(c,k_zwt)/(dz(c,k_zwt)*denice))/watsat(c,k_zwt)
             s2 = (h2osoi_liq(c,k_zwt+1)/(dz(c,k_zwt+1)*denh2o) &
                  + h2osoi_ice(c,k_zwt+1)/(dz(c,k_zwt+1)*denice))/watsat(c,k_zwt+1)
             
             m=(z(c,k_zwt+1)-z(c,k_zwt))/(s2-s1)
             b=z(c,k_zwt+1)-m*s2
             zwt(c)=max(0._r8,m*sat_lev+b)
          else
             zwt(c)=zi(c,nbedrock(c))
          endif
       end do

     end associate

   end subroutine ThetaBasedWaterTable

!#6
   !-----------------------------------------------------------------------
   subroutine LateralFlowPowerLaw(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc,soilhydrology_inst, soilstate_inst, &
        waterstatebulk_inst, waterfluxbulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage
     !
     ! !USES:
     use clm_varcon       , only : pondmx, watmin,rpi, secspday, nlvic
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv
     use GridcellType     , only : grc                

     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'Drainage'           ! subroutine name
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)            ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: xsia(bounds%begc:bounds%endc)           ! available pore space at layer i (mm)
     real(r8) :: xs1(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer 1 (mm)
     real(r8) :: smpfz(1:nlevsoi)                        ! matric potential of layer right above water table (mm)
     real(r8) :: wtsub                                   ! summation of hk*dzmm for layers below water table (mm**2/s)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: h2osoi_vol
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: vol_ice
     real(r8) :: dsmax_tmp(bounds%begc:bounds%endc)       ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rsub_tmp                 ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: frac                     ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rel_moist                ! relative moisture, temporary variable
     real(r8) :: wtsub_vic                ! summation of hk*dzmm for layers in the third VIC layer
     integer :: g
     real(r8), parameter :: n_baseflow = 1 !drainage power law exponent
     !-----------------------------------------------------------------------

     associate(                                                            & 
          nbedrock           =>    col%nbedrock                          , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)           
          h2osno             =>    waterstatebulk_inst%h2osno_col            , & ! Input:  [real(r8) (:)   ] surface water (mm)                                
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          snl                =>    col%snl                               , & ! Input:  [integer  (:)   ] number of snow layers                              
          h2osfc             =>    waterstatebulk_inst%h2osfc_col            , & ! Input:  [real(r8) (:)   ] surface water (mm)                                
          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  
          eff_porosity       =>    soilstate_inst%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ] effective porosity = porosity - vol_ice         
          hk_l               =>    soilstate_inst%hk_l_col               , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity (mm/s)                    

          depth              =>    soilhydrology_inst%depth_col          , & ! Input:  [real(r8) (:,:) ] VIC soil depth                                   
          c_param            =>    soilhydrology_inst%c_param_col        , & ! Input:  [real(r8) (:)   ] baseflow exponent (Qb)                             
          Dsmax              =>    soilhydrology_inst%dsmax_col          , & ! Input:  [real(r8) (:)   ] max. velocity of baseflow (mm/day)
          max_moist          =>    soilhydrology_inst%max_moist_col      , & ! Input:  [real(r8) (:,:) ] maximum soil moisture (ice + liq)
          moist              =>    soilhydrology_inst%moist_col          , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)                         
          Ds                 =>    soilhydrology_inst%ds_col             , & ! Input:  [real(r8) (:)   ] fracton of Dsmax where non-linear baseflow begins
          Wsvic              =>    soilhydrology_inst%Wsvic_col          , & ! Input:  [real(r8) (:)   ] fraction of maximum soil moisutre where non-liear base flow occurs
          icefrac            =>    soilhydrology_inst%icefrac_col        , & ! Output: [real(r8) (:,:) ] fraction of ice in layer                         
          hkdepth            =>    soilhydrology_inst%hkdepth_col        , & ! Input:  [real(r8) (:)   ] decay factor (m)                                  
          frost_table        =>    soilhydrology_inst%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)                             
          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)                             
          wa                 =>    waterstatebulk_inst%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)              
          ice                =>    soilhydrology_inst%ice_col            , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)                         
          qcharge            =>    soilhydrology_inst%qcharge_col        , & ! Input:  [real(r8) (:)   ] aquifer recharge rate (mm/s)                      
          origflag           =>    soilhydrology_inst%origflag           , & ! Input:  logical
          h2osfcflag         =>    soilhydrology_inst%h2osfcflag         , & ! Input:  integer
          
          qflx_snwcp_liq     =>    waterfluxbulk_inst%qflx_snwcp_liq_col     , & ! Output: [real(r8) (:)   ] excess rainfall due to snow capping (mm H2O /s) [+]
          qflx_ice_runoff_xs =>    waterfluxbulk_inst%qflx_ice_runoff_xs_col , & ! Output: [real(r8) (:)   ] solid runoff from excess ice in soil (mm H2O /s) [+]
          qflx_dew_grnd      =>    waterfluxbulk_inst%qflx_dew_grnd_col      , & ! Output: [real(r8) (:)   ] ground surface dew formation (mm H2O /s) [+]      
          qflx_dew_snow      =>    waterfluxbulk_inst%qflx_dew_snow_col      , & ! Output: [real(r8) (:)   ] surface dew added to snow pack (mm H2O /s) [+]    
          qflx_sub_snow      =>    waterfluxbulk_inst%qflx_sub_snow_col      , & ! Output: [real(r8) (:)   ] sublimation rate from snow pack (mm H2O /s) [+]   
          qflx_drain         =>    waterfluxbulk_inst%qflx_drain_col         , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)                    
          qflx_qrgwl         =>    waterfluxbulk_inst%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ] qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
          qflx_rsub_sat      =>    waterfluxbulk_inst%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ] soil saturation excess [mm h2o/s]                 
          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)                                
          )

       ! Get time step

       dtime = get_step_size()

       ! Convert layer thicknesses from m to mm

       do j = 1,nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             dzmm(c,j) = dz(c,j)*1.e3_r8

             vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             icefrac(c,j) = min(1._r8,vol_ice/watsat(c,j))          
          end do
       end do

       ! Initial set

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          qflx_drain(c)    = 0._r8 
          rsub_bot(c)      = 0._r8
          qflx_rsub_sat(c) = 0._r8
          rsub_top(c)      = 0._r8
          fracice_rsub(c)  = 0._r8
       end do

       ! The layer index of the first unsaturated layer, 
       ! i.e., the layer right above the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          jwt(c) = nlevsoi
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevsoi
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       end do

       !-- Topographic runoff  -------------------------
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          fff(c)         = 1._r8/ hkdepth(c)

          dzsum = 0._r8
          icefracsum = 0._r8
          do j = max(jwt(c),1), nlevsoi
             dzsum  = dzsum + dzmm(c,j)
             icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
          end do
          imped=10._r8**(-e_ice*(icefracsum/dzsum))
          !@@
          ! baseflow is power law expression relative to bedrock layer
          if(zwt(c) <= zi(c,nbedrock(c))) then 
             rsub_top(c)    = imped * baseflow_scalar * tan(rpi/180._r8*col%topo_slope(c))* &
                              (zi(c,nbedrock(c)) - zwt(c))**(n_baseflow)
          else
             rsub_top(c) = 0._r8
          endif

          !--  Now remove water via rsub_top
          rsub_top_tot = - rsub_top(c)* dtime

          !should never be positive... but include for completeness
          if(rsub_top_tot > 0.) then !rising water table
             
             call endrun(msg="RSUB_TOP IS POSITIVE in Drainage!"//errmsg(sourcefile, __LINE__))
             
          else ! deepening water table
             do j = jwt(c)+1, nbedrock(c)
                ! use analytical expression for specific yield
                s_y = watsat(c,j) &
                     * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                s_y=max(s_y,0.02_r8)
                rsub_top_layer=max(rsub_top_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                rsub_top_layer=min(rsub_top_layer,0._r8)
                h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer
                   
                rsub_top_tot = rsub_top_tot - rsub_top_layer

                if (rsub_top_tot >= 0.) then 
                   zwt(c) = zwt(c) - rsub_top_layer/s_y/1000._r8
                   
                   exit
                else
                   zwt(c) = zi(c,j)
                endif
             enddo
             
             !--  remove residual rsub_top  --------------------------------
             ! make sure no extra water removed from soil column
             rsub_top(c) = rsub_top(c) - rsub_top_tot/dtime
          endif
          
          zwt(c) = max(0.0_r8,zwt(c))
          zwt(c) = min(80._r8,zwt(c))
       end do

       !  excessive water above saturation added to the above unsaturated layer like a bucket
       !  if column fully saturated, excess water goes to runoff

       do j = nlevsoi,2,-1
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             xsi(c)            = max(h2osoi_liq(c,j)-eff_porosity(c,j)*dzmm(c,j),0._r8)
             h2osoi_liq(c,j)   = min(eff_porosity(c,j)*dzmm(c,j), h2osoi_liq(c,j))
             h2osoi_liq(c,j-1) = h2osoi_liq(c,j-1) + xsi(c)
          end do
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! watmin addition to fix water balance errors
          xs1(c) = max(max(h2osoi_liq(c,1)-watmin,0._r8)- &
               max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice(c,1)-watmin)),0._r8)
          h2osoi_liq(c,1) = h2osoi_liq(c,1) - xs1(c)

          if (lun%urbpoi(col%landunit(c))) then
             qflx_rsub_sat(c)     = xs1(c) / dtime
          else
             ! send this water up to h2osfc rather than sending to drainage
             h2osfc(c) = h2osfc(c) + xs1(c)
             qflx_rsub_sat(c)     = 0._r8
          endif
          ! add in ice check
          xs1(c)          = max(max(h2osoi_ice(c,1),0._r8)-max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1))),0._r8)
          h2osoi_ice(c,1) = min(max(0._r8,pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1)), h2osoi_ice(c,1))
          qflx_ice_runoff_xs(c) = xs1(c) / dtime
       end do

       ! Limit h2osoi_liq to be greater than or equal to watmin.
       ! Get water needed to bring h2osoi_liq equal watmin from lower layer.
       ! If insufficient water in soil layers, get from aquifer water

       do j = 1, nlevsoi-1
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             if (h2osoi_liq(c,j) < watmin) then
                xs(c) = watmin - h2osoi_liq(c,j)
                ! deepen water table if water is passed from below zwt layer
                if(j == jwt(c)) then 
                   zwt(c) = zwt(c) + xs(c)/eff_porosity(c,j)/1000._r8
                endif
             else
                xs(c) = 0._r8
             end if
             h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
             h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
          end do
       end do

       ! Get water for bottom layer from layers above if possible
       j = nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (h2osoi_liq(c,j) < watmin) then
             xs(c) = watmin-h2osoi_liq(c,j)
             searchforwater: do i = nlevsoi-1, 1, -1
                available_h2osoi_liq = max(h2osoi_liq(c,i)-watmin-xs(c),0._r8)
                if (available_h2osoi_liq >= xs(c)) then
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
                   h2osoi_liq(c,i) = h2osoi_liq(c,i) - xs(c)
                   xs(c) = 0._r8
                   exit searchforwater
                else
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + available_h2osoi_liq
                   h2osoi_liq(c,i) = h2osoi_liq(c,i) - available_h2osoi_liq
                   xs(c) = xs(c) - available_h2osoi_liq
                end if
             end do searchforwater
          else
             xs(c) = 0._r8
          end if
          ! Needed in case there is no water to be found
          h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
          ! Instead of removing water from aquifer where it eventually
          ! shows up as excess drainage to the ocean, take it back out of 
          ! drainage
          qflx_rsub_sat(c) = qflx_rsub_sat(c) - xs(c)/dtime

       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Sub-surface runoff and drainage

          qflx_drain(c) = qflx_rsub_sat(c) + rsub_top(c)

          ! Set imbalance for snow capping

          qflx_qrgwl(c) = qflx_snwcp_liq(c)

       end do

       ! No drainage for urban columns (except for pervious road as computed above)

       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col%itype(c) /= icol_road_perv) then
             qflx_drain(c) = 0._r8
             ! This must be done for roofs and impervious road (walls will be zero)
             qflx_qrgwl(c) = qflx_snwcp_liq(c)
          end if
       end do

     end associate

   end subroutine LateralFlowPowerLaw

!#7
   !-----------------------------------------------------------------------
   subroutine RenewCondensation(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst, &
        waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst) 
     !
     ! !DESCRIPTION:
     ! Calculate watertable, considering aquifer recharge but no drainage.
     !
     ! !USES:
     use column_varcon    , only : icol_roof, icol_road_imperv
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds  
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(waterdiagnosticbulk_type)    , intent(inout) :: waterdiagnosticbulk_inst
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     !-----------------------------------------------------------------------

     associate(                                                            & 
          snl                =>    col%snl                               , & ! Input:  [integer  (:)   ]  number of snow layers                              
          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
          frac_h2osfc        =>    waterdiagnosticbulk_inst%frac_h2osfc_col       , & ! Input:  [real(r8) (:)   ]                                                    
          qflx_dew_grnd      =>    waterfluxbulk_inst%qflx_dew_grnd_col      , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]      
          qflx_dew_snow      =>    waterfluxbulk_inst%qflx_dew_snow_col      , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]    
          qflx_sub_snow      =>    waterfluxbulk_inst%qflx_sub_snow_col       & ! Output: [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]   
          )

       ! Get time step

       dtime = get_step_size()

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Renew the ice and liquid mass due to condensation

          if (snl(c)+1 >= 1) then

             ! make consistent with how evap_grnd removed in infiltration
             h2osoi_liq(c,1) = h2osoi_liq(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_grnd(c) * dtime
             h2osoi_ice(c,1) = h2osoi_ice(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_snow(c) * dtime
             if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
                qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
                h2osoi_ice(c,1) = 0._r8
             else
                h2osoi_ice(c,1) = h2osoi_ice(c,1) - (1._r8 - frac_h2osfc(c)) * qflx_sub_snow(c) * dtime
             end if
          end if

       end do


       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          ! Renew the ice and liquid mass due to condensation for urban roof and impervious road

          if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
             if (snl(c)+1 >= 1) then
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + qflx_dew_grnd(c) * dtime
                h2osoi_ice(c,1) = h2osoi_ice(c,1) + (qflx_dew_snow(c) * dtime)
                if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
                   qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
                   h2osoi_ice(c,1) = 0._r8
                else
                   h2osoi_ice(c,1) = h2osoi_ice(c,1) - (qflx_sub_snow(c) * dtime)
                end if
             end if
          end if

       end do

     end associate

   end subroutine RenewCondensation
!#8
   !-----------------------------------------------------------------------
   subroutine CalcIrrigWithdrawals(bounds, &
        num_soilc, filter_soilc, &
        soilhydrology_inst, soilstate_inst, &
        qflx_gw_demand, &
        qflx_gw_uncon_irrig_lyr, &
        qflx_gw_con_irrig)
     !
     ! !DESCRIPTION:
     ! Calculate irrigation withdrawals from groundwater, given groundwater irrigation demand
     !
     ! This routine is called when use_groundwater_irrigation = .true.
     !
     ! It is not compatible with use_aquifer_layer
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_soilc       ! number of column soil points in column filter
     integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil points
     type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
     type(soilstate_type)     , intent(in)    :: soilstate_inst

     real(r8) , intent(in)    :: qflx_gw_demand( bounds%begc: )              ! groundwater irrigation demand (mm H2O/s)
     real(r8) , intent(inout) :: qflx_gw_uncon_irrig_lyr( bounds%begc:, 1: ) ! unconfined aquifer groundwater irrigation withdrawal flux, separated by layer (mm H2O/s)
     real(r8) , intent(inout) :: qflx_gw_con_irrig( bounds%begc: )           ! total confined aquifer groundwater irrigation withdrawal flux (mm H2O/s)
     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'CalcIrrigWithdrawals'  ! subroutine name
     integer  :: c,j,fc                                  ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: s_y
     real(r8) :: irrig_demand_remaining  ! mm H2O
     real(r8) :: available_water_layer   ! mm H2O
     real(r8) :: irrig_layer             ! mm H2O
     !-----------------------------------------------------------------------

     SHR_ASSERT_ALL((ubound(qflx_gw_demand) == [bounds%endc]), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(qflx_gw_uncon_irrig_lyr) == [bounds%endc, nlevsoi]), errMsg(sourcefile, __LINE__))
     SHR_ASSERT_ALL((ubound(qflx_gw_con_irrig) == [bounds%endc]), errMsg(sourcefile, __LINE__))

     associate(                                                            & 
          nbedrock           =>    col%nbedrock                          , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)           
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  
          zwt                =>    soilhydrology_inst%zwt_col              & ! Input:  [real(r8) (:)   ] water table depth (m)                             
          )

       dtime = get_step_size()

       do j = 1, nlevsoi
          do fc = 1, num_soilc
             c = filter_soilc(fc)
             qflx_gw_uncon_irrig_lyr(c,j) = 0._r8
          end do
       end do

       !-- Remove groundwater from unconfined aquifer  -----------
       do fc = 1, num_soilc
          c = filter_soilc(fc)

          irrig_demand_remaining = qflx_gw_demand(c)*dtime
          
          ! should never be negative... but include for completeness
          if(irrig_demand_remaining < 0.) then
             
             call endrun(msg="negative groundwater irrigation demand! "//errmsg(sourcefile, __LINE__))
             
          else 
             jwt(c) = nlevsoi
             ! allow jwt to equal zero when zwt is in top layer
             do j = 1,nlevsoi
                if(zwt(c) <= zi(c,j)) then
                   jwt(c) = j-1
                   exit
                end if
             enddo
             do j = jwt(c)+1, nbedrock(c)
                ! use analytical expression for specific yield
                s_y = watsat(c,j) &
                     * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                s_y=max(s_y,0.02_r8)

                if (j==jwt(c)+1) then
                   available_water_layer=max(0._r8,(s_y*(zi(c,j) - zwt(c))*1.e3))
                else
                   available_water_layer=max(0._r8,(s_y*(zi(c,j) - zi(c,j-1))*1.e3))
                endif

                irrig_layer = min(irrig_demand_remaining, available_water_layer)
                qflx_gw_uncon_irrig_lyr(c,j) = irrig_layer / dtime

                irrig_demand_remaining = irrig_demand_remaining - irrig_layer

                if (irrig_demand_remaining <= 0.) then 
                   exit
                endif
             end do

             if (irrig_demand_remaining > 0._r8) then
                ! NOTE(wjs, 2018-12-15) Eventually we could limit irrigation from the
                ! confined aquifer based on water availability there.
                qflx_gw_con_irrig(c) = irrig_demand_remaining / dtime
             else
                qflx_gw_con_irrig(c) = 0._r8
             end if
          end if
       end do

     end associate

   end subroutine CalcIrrigWithdrawals
!#9
   !-----------------------------------------------------------------------
   subroutine WithdrawGroundwaterIrrigation(bounds, &
        num_soilc, filter_soilc, &
        waterflux_inst, waterstate_inst)
     !
     ! !DESCRIPTION:
     ! Remove groundwater irrigation from unconfined and confined aquifers
     ! This routine is called when use_groundwater_irrigation = .true.
     ! It is not compatible with use_aquifer_layer
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds               
     integer                , intent(in)    :: num_soilc       ! number of column soil points in column filter
     integer                , intent(in)    :: filter_soilc(:) ! column filter for soil points
     class(waterflux_type)  , intent(in)    :: waterflux_inst
     class(waterstate_type) , intent(inout) :: waterstate_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: fc, c, j
     real(r8) :: dtime        ! land model time step (sec)

     character(len=*), parameter :: subname = 'WithdrawGroundwaterIrrigation'
     !-----------------------------------------------------------------------

     associate( &
          qflx_gw_uncon_irrig_lyr => waterflux_inst%qflx_gw_uncon_irrig_lyr_col, & ! Input: [real(r8) (:,:) ] unconfined groundwater irrigation flux, separated by layer (mm H2O/s)
          qflx_gw_con_irrig  => waterflux_inst%qflx_gw_con_irrig_col     , & ! Input: [real(r8) (:,:) ] confined groundwater irrigation flux (mm H2O /s)

          wa                 =>    waterstate_inst%wa_col                , & ! Output: [real(r8) (:)   ] water in the unconfined aquifer (mm)
          h2osoi_liq         =>    waterstate_inst%h2osoi_liq_col          & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          )

     dtime = get_step_size()

     do j = 1, nlevsoi
        do fc = 1, num_soilc
           c = filter_soilc(fc)
           h2osoi_liq(c,j) = h2osoi_liq(c,j) - qflx_gw_uncon_irrig_lyr(c,j) * dtime
        end do
     end do

     ! zwt is not being updated, as it will be updated after HydrologyNoDrainage
       
     do fc = 1, num_soilc
        c = filter_soilc(fc)
        wa(c) = wa(c) - qflx_gw_con_irrig(c) * dtime
     end do

     end associate

   end subroutine WithdrawGroundwaterIrrigation

!#0
end module SoilHydrologyMod
