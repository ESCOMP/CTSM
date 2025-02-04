module SoilHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate soil hydrology
  !
#include "shr_assert.h"
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use decompMod         , only : bounds_type, subgrid_level_column
  use clm_varctl        , only : iulog, use_vichydro
  use clm_varcon        , only : ispval
  use clm_varcon        , only : denh2o, denice, rpi
  use clm_varcon        , only : pondmx_urban
  use clm_varpar        , only : nlevsoi, nlevgrnd, nlayer, nlayert
  use column_varcon     , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon     , only : icol_road_imperv
  use landunit_varcon   , only : istsoil, istcrop
  use clm_time_manager  , only : get_step_size_real
  use NumericsMod       , only : truncate_small_values
  use EnergyFluxType    , only : energyflux_type
  use InfiltrationExcessRunoffMod, only : infiltration_excess_runoff_type
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use Wateratm2lndBulkType, only : wateratm2lndbulk_type
  use WaterFluxType     , only : waterflux_type
  use WaterFluxBulkType , only : waterfluxbulk_type
  use WaterStateType    , only : waterstate_type
  use WaterStateBulkType, only : waterstatebulk_type
  use WaterDiagnosticBulkType, only : waterdiagnosticbulk_type
  use TemperatureType   , only : temperature_type
  use LandunitType      , only : lun                
  use ColumnType        , only : column_type, col
  use PatchType         , only : patch

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilHydReadNML       ! Read in the Soil hydrology namelist
  public :: SetSoilWaterFractions ! Set diagnostic variables related to the fraction of water and ice in each layer
  public :: SetFloodc            ! Apply gridcell flood water flux to non-lake columns
  public :: SetQflxInputs        ! Set the flux of water into the soil from the top
  public :: Infiltration         ! Calculate total infiltration
  public :: TotalSurfaceRunoff   ! Calculate total surface runoff
  public :: UpdateUrbanPonding   ! Update the state variable representing ponding on urban surfaces
  public :: WaterTable           ! Calculate water table before imposing drainage
  public :: Drainage             ! Calculate subsurface drainage
  public :: CLMVICMap
  public :: PerchedWaterTable    ! Calculate perched water table
  public :: PerchedLateralFlow   ! Calculate lateral flow from perched saturated zone
  public :: ThetaBasedWaterTable ! Calculate water table from soil moisture state
  public :: SubsurfaceLateralFlow ! Calculate subsurface lateral flow from saturated zone
  public :: RenewCondensation    ! Misc. corrections
  public :: CalcIrrigWithdrawals ! Calculate irrigation withdrawals from groundwater by layer
  public :: WithdrawGroundwaterIrrigation   ! Remove groundwater irrigation from unconfined and confined aquifers
  public :: readParams

  type, private :: params_type
     real(r8) :: aq_sp_yield_min         ! Minimum aquifer specific yield (unitless)
     real(r8) :: n_baseflow              ! Drainage power law exponent (unitless)
     real(r8) :: perched_baseflow_scalar ! Scalar multiplier for perched base flow rate (kg/m2/s)
     real(r8) :: e_ice                   ! Soil ice impedance factor (unitless)
  end type params_type
  type(params_type), public ::  params_inst
  
  !-----------------------------------------------------------------------
  real(r8), private   :: baseflow_scalar = 1.e-2_r8
  real(r8), parameter :: tolerance = 1.e-12_r8                   ! tolerance for checking whether sublimation is greater than ice in top soil layer

  integer, private :: head_gradient_method    ! Method for calculating hillslope saturated head gradient
  integer, private :: transmissivity_method   ! Method for calculating transmissivity of hillslope columns

  ! Head gradient methods
  integer, parameter, private :: kinematic = 0
  integer, parameter, private :: darcy     = 1
  ! Transmissivity methods
  integer, parameter, private :: uniform_transmissivity = 0
  integer, parameter, private :: layersum = 1

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine hillslope_hydrology_ReadNML(NLFilename)
    !
    ! DESCRIPTION
    ! read in hillslope hydrology namelist variables related to
    ! subsurface lateral flow
    !
    ! !USES:
    use abortutils      , only : endrun
    use fileutils       , only : getavu, relavu
    use spmdMod         , only : mpicom, masterproc
    use shr_mpi_mod     , only : shr_mpi_bcast
    use clm_varctl      , only : iulog
    use clm_nlUtilsMod  , only : find_nlgroup_name

    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !--------------------------------------------------------------------
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    character(len=*), parameter :: nmlname = 'hillslope_hydrology_inparm'
    character(*), parameter    :: subName = "('hillslope_hydrology_ReadNML')"
    character(len=50) :: hillslope_head_gradient_method  = 'Darcy'    ! head gradient method string
    character(len=50) :: hillslope_transmissivity_method = 'LayerSum' ! transmissivity method string
    !-----------------------------------------------------------------------

! MUST agree with name in namelist and read statement
    namelist /hillslope_hydrology_inparm/ &
         hillslope_head_gradient_method,  &
         hillslope_transmissivity_method
    
    ! Default values for namelist
    head_gradient_method    = darcy
    transmissivity_method   = layersum
    
    ! Read hillslope hydrology namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'hillslope_hydrology_inparm', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=hillslope_hydrology_inparm,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading hillslope hydrology namelist')
          end if
       else
          call endrun(subname // ':: ERROR reading hillslope hydrology namelist')
       end if
       close(nu_nml)
       call relavu( nu_nml )

       ! Convert namelist strings to numerical values
       if (      trim(hillslope_head_gradient_method) == 'Kinematic' ) then
          head_gradient_method = kinematic
       else if ( trim(hillslope_head_gradient_method) == 'Darcy'     ) then
          head_gradient_method = darcy
       else
          call endrun(msg="ERROR bad value for hillslope_head_gradient_method in "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if

       if (      trim(hillslope_transmissivity_method) == 'Uniform' ) then
          transmissivity_method = uniform_transmissivity
       else if ( trim(hillslope_transmissivity_method) == 'LayerSum') then
          transmissivity_method = layersum
       else
          call endrun(msg="ERROR bad value for hillslope_transmissivity_method in "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if

    endif

    call shr_mpi_bcast(head_gradient_method, mpicom)
    call shr_mpi_bcast(transmissivity_method, mpicom)

    if (masterproc) then

       write(iulog,*) ' '
       write(iulog,*) 'hillslope_hydrology lateral flow settings:'
       write(iulog,*) '  hillslope_head_gradient_method  = ',hillslope_head_gradient_method
       write(iulog,*) '  hillslope_transmissivity_method = ',hillslope_transmissivity_method

    endif

  end subroutine hillslope_hydrology_ReadNML

  !-----------------------------------------------------------------------
  subroutine readParams( ncid )
    !
    ! !USES:
    use ncdio_pio, only: file_desc_t
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'readParams_SoilHydrology'
    !--------------------------------------------------------------------

    ! Minimum aquifer specific yield (unitless)
    call readNcdioScalar(ncid, 'aq_sp_yield_min', subname, params_inst%aq_sp_yield_min)
    ! Drainage power law exponent (unitless)
    call readNcdioScalar(ncid, 'n_baseflow', subname, params_inst%n_baseflow)
    ! Scalar multiplier for perched base flow rate (kg/m2/s)
    call readNcdioScalar(ncid, 'perched_baseflow_scalar', subname, params_inst%perched_baseflow_scalar)
    ! Soil ice impedance factor (unitless)
    call readNcdioScalar(ncid, 'e_ice', subname, params_inst%e_ice)

  end subroutine readParams

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
    real(r8) :: dz_ext(bounds%begc:bounds%endc,1:nlevsoi)

    character(len=*), parameter :: subname = 'SetSoilWaterFractions'
    !-----------------------------------------------------------------------

    associate( &
         dz               =>    col%dz                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 

         watsat           =>    soilstate_inst%watsat_col           , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         eff_porosity     =>    soilstate_inst%eff_porosity_col     , & ! Output: [real(r8) (:,:) ]  effective porosity = porosity - vol_ice

         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col  , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col  , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)
         excess_ice       =>    waterstatebulk_inst%excess_ice_col  , & ! Input:  [real(r8) (:,:) ]  excess ice (kg/m2)
         icefrac          =>    soilhydrology_inst%icefrac_col        & ! Output: [real(r8) (:,:) ]                                                  
         )

    do j = 1,nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Porosity of soil, partial volume of ice and liquid, fraction of ice in each layer,
          ! fractional impermeability
          dz_ext(c,j) = dz(c,j) + excess_ice(c,j)/denice ! extended layer thickness, should be good for all the columns
          vol_ice(c,j) = min(watsat(c,j), (h2osoi_ice(c,j) + excess_ice(c,j))/(dz_ext(c,j)*denice))
          eff_porosity(c,j) = max(0.01_r8,watsat(c,j)-vol_ice(c,j))
          icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))

       end do
    end do

    end associate

  end subroutine SetSoilWaterFractions

  !-----------------------------------------------------------------------
  subroutine SetFloodc(num_nolakec, filter_nolakec, col, &
       wateratm2lndbulk_inst, waterfluxbulk_inst)
    !
    ! !DESCRIPTION:
    ! Apply gridcell flood water flux to non-lake columns
    !
    ! !ARGUMENTS:
     integer                     , intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
     integer                     , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
     type(column_type)           , intent(in)    :: col
     type(wateratm2lndbulk_type) , intent(in)    :: wateratm2lndbulk_inst
     type(waterfluxbulk_type)    , intent(inout) :: waterfluxbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, g

    character(len=*), parameter :: subname = 'SetFloodc'
    !-----------------------------------------------------------------------

    associate( &
          qflx_floodg => wateratm2lndbulk_inst%forc_flood_grc , & ! Input:  [real(r8) (:)   ]  gridcell flux of flood water from RTM
          qflx_floodc => waterfluxbulk_inst%qflx_floodc_col     & ! Output: [real(r8) (:)   ]  column flux of flood water from RTM
          )

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       g = col%gridcell(c)
       if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall) then
          qflx_floodc(c) = qflx_floodg(g)
       else
          qflx_floodc(c) = 0._r8
       end if
    end do

    end associate

  end subroutine SetFloodc

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
          qflx_liqevap_from_top_layer => waterfluxbulk_inst%qflx_liqevap_from_top_layer_col, & ! Input:  [real(r8) (:)]  rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]
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
           qflx_evap=qflx_liqevap_from_top_layer(c)
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
          qflx_liqevap_from_top_layer => waterfluxbulk_inst%qflx_liqevap_from_top_layer_col, & ! Input:  [real(r8) (:)   ]  rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]    
          qflx_floodc      =>    waterfluxbulk_inst%qflx_floodc_col      , & ! Input:  [real(r8) (:)   ]  column flux of flood water from RTM               
          qflx_sat_excess_surf => waterfluxbulk_inst%qflx_sat_excess_surf_col , & ! Input:  [real(r8) (:)   ]  surface runoff due to saturated surface (mm H2O /s)

          xs_urban         =>    soilhydrology_inst%xs_urban_col     , & ! Output: [real(r8) (:)   ]  excess soil water above urban ponding limit

          h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col        & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                            
          )

     dtime = get_step_size_real()

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
                   h2osoi_liq(c,1)/dtime + qflx_rain_plus_snomelt(c) - qflx_liqevap_from_top_layer(c) - &
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
         qflx_liqevap_from_top_layer => waterfluxbulk_inst%qflx_liqevap_from_top_layer_col & ! Input:  [real(r8) (:)   ]  rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]    
         )

     dtime = get_step_size_real()

     do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        
        if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
           if (snl(c) >= 0) then
              if (xs_urban(c) > 0.) then
                 h2osoi_liq(c,1) = pondmx_urban
              else
                 h2osoi_liq(c,1) = max(0._r8,h2osoi_liq(c,1)+ &
                      (qflx_rain_plus_snomelt(c)-qflx_liqevap_from_top_layer(c))*dtime)
              end if
           end if
        end if
     end do

     end associate

   end subroutine UpdateUrbanPonding

   !-----------------------------------------------------------------------
   subroutine WaterTable(bounds, num_hydrologyc, filter_hydrologyc, &
        soilhydrology_inst, soilstate_inst, temperature_inst, waterstatebulk_inst, waterfluxbulk_inst)
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
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(temperature_type)   , intent(in)    :: temperature_inst
     type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type)     , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
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
     real(r8) :: qflx_solidevap_from_top_layer_save      ! temporary
     !-----------------------------------------------------------------------

     associate(                                                            & 
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           

          t_soisno           =>    temperature_inst%t_soisno_col         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       

          h2osfc             =>    waterstatebulk_inst%h2osfc_col            , & ! Input:  [real(r8) (:)   ]  surface water (mm)                                
          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
          h2osoi_vol         =>    waterstatebulk_inst%h2osoi_vol_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]

          qflx_ev_snow       =>    waterfluxbulk_inst%qflx_ev_snow_col   , & ! In/Out: [real(r8) (:)   ]  evaporation flux from snow (mm H2O/s) [+ to atm]
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
          
          qflx_drain         =>    waterfluxbulk_inst%qflx_drain_col         , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)
          qflx_drain_perched =>    waterfluxbulk_inst%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ]  perched wt sub-surface runoff (mm H2O /s)         
          qflx_rsub_sat      =>    waterfluxbulk_inst%qflx_rsub_sat_col        & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]                 
          )

       ! Get time step

       dtime = get_step_size_real()

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
          rous=max(rous, params_inst%aq_sp_yield_min)

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
                   s_y=max(s_y, params_inst%aq_sp_yield_min)

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
                   s_y=max(s_y, params_inst%aq_sp_yield_min)

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
          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz) then
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
          h2osfcflag         =>    soilhydrology_inst%h2osfcflag         , & ! Input:  integer
          
          qflx_snwcp_liq     =>    waterfluxbulk_inst%qflx_snwcp_liq_col     , & ! Output: [real(r8) (:)   ] excess liquid h2o due to snow capping (outgoing) (mm H2O /s) [+]
          qflx_ice_runoff_xs =>    waterfluxbulk_inst%qflx_ice_runoff_xs_col , & ! Output: [real(r8) (:)   ] solid runoff from excess ice in soil (mm H2O /s) [+]
          qflx_liqdew_to_top_layer      => waterfluxbulk_inst%qflx_liqdew_to_top_layer_col     , & ! Output: [real(r8) (:)   ] rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s) [+]    
          qflx_soliddew_to_top_layer    => waterfluxbulk_inst%qflx_soliddew_to_top_layer_col   , & ! Output: [real(r8) (:)   ] rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s) [+]      
          qflx_solidevap_from_top_layer => waterfluxbulk_inst%qflx_solidevap_from_top_layer_col, & ! Output: [real(r8) (:)   ] rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s) [+]   
          qflx_drain         =>    waterfluxbulk_inst%qflx_drain_col         , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)                    
          qflx_qrgwl         =>    waterfluxbulk_inst%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ] qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
          qflx_rsub_sat      =>    waterfluxbulk_inst%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ] soil saturation excess [mm h2o/s]                 
          qflx_drain_perched =>    waterfluxbulk_inst%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ] perched wt sub-surface runoff (mm H2O /s)         

          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)                                
          )

       ! Get time step

       dtime = get_step_size_real()

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
          qflx_rsub_sat(c) = 0._r8
          rsub_top(c)      = 0._r8
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
          q_perch_max = params_inst%perched_baseflow_scalar * sin(col%topo_slope(c) * (rpi/180._r8))

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

          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz) then
             ! compute drainage from perched saturated region
             wtsub = 0._r8
             q_perch = 0._r8
             do k = jwt(c)+1, k_frz
                imped=10._r8**(-params_inst%e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
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
                   imped=10._r8**(-params_inst%e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                   q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                   wtsub = wtsub + dzmm(c,k)
                end do
                if (wtsub > 0._r8) q_perch = q_perch/wtsub

                qflx_drain_perched(c) = q_perch_max * q_perch &
                     *(frost_table(c) - zwt_perched(c))

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
             if (use_vichydro) then
                imped=10._r8**(-params_inst%e_ice*min(1.0_r8,ice(c,nlayer)/max_moist(c,nlayer)))
                dsmax_tmp(c) = Dsmax(c) * dtime/ secspday !mm/day->mm/dtime
                rsub_top_max = dsmax_tmp(c)
             else
                imped=10._r8**(-params_inst%e_ice*(icefracsum/dzsum))
                rsub_top_max = 10._r8 * sin((rpi/180.) * col%topo_slope(c))
             end if

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
             rous=max(rous, params_inst%aq_sp_yield_min)

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

                   call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                        msg="RSUB_TOP IS POSITIVE in Drainage!"//errmsg(sourcefile, __LINE__))

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
                         s_y=max(s_y, params_inst%aq_sp_yield_min)

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
     use clm_varcon       , only : tfrz, denice, denh2o
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
     type(waterstatebulk_type), intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type) , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i                       ! indices
     integer  :: k,k_frz,k_perch,k_zwt          ! indices
     real(r8) :: s1, s2                         ! temporary moisture values
     real(r8) :: m, b                           ! slope and intercept
     real(r8), parameter :: sat_lev = 0.9       ! saturation value used to identify saturated layers
     !-----------------------------------------------------------------------

     associate(                                                            & 
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)
          t_soisno           =>    temperature_inst%t_soisno_col         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       

          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
          h2osoi_vol         =>    waterstatebulk_inst%h2osoi_vol_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  
          zwt                =>    soilhydrology_inst%zwt_col            , & ! Output: [real(r8) (:)   ]  water table depth (m)                             
          zwt_perched        =>    soilhydrology_inst%zwt_perched_col    , & ! Output: [real(r8) (:)   ]  perched water table depth (m)                     
          frost_table        =>    soilhydrology_inst%frost_table_col      & ! Output: [real(r8) (:)   ]  frost table depth (m)                             
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

          ! frost table is top of frozen layer
          frost_table(c) = zi(c,k_frz-1)

          ! initialize perched water table to frost table
          zwt_perched(c) = frost_table(c)

          !=======  water table above frost table  ===================
          ! if water table is above frost table, do nothing 
          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz) then
          else if (k_frz > 1) then
             !==========  water table below frost table  ============
             ! locate perched water table from bottom up starting at 
             ! frost table sat_lev is an arbitrary saturation level 
             ! used to determine perched water table

             k_perch = 1
             do k=k_frz,1,-1
                h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) then 
                   k_perch = k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, check temperature of layer, 
             ! and only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists above frost table, 
             ! interpolate between k_perch and k_perch+1 to find 
             ! perched water table height
             if (k_frz > k_perch) then
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)

                if (s1 > s2) then 
                   zwt_perched(c) = zi(c,k_perch-1)
                else
                   m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                   b=z(c,k_perch+1)-m*s2
                   zwt_perched(c)=max(0._r8,m*sat_lev+b)
                endif
             endif
          endif
       end do

     end associate

   end subroutine PerchedWaterTable

!#4   
   !-----------------------------------------------------------------------
   subroutine PerchedLateralFlow(bounds, num_hydrologyc, &
        filter_hydrologyc, soilhydrology_inst, soilstate_inst, &
        waterstatebulk_inst, waterfluxbulk_inst, wateratm2lndbulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage from perched saturated zone
     !
     ! !USES:
     use clm_varcon       , only : pondmx, tfrz, watmin,rpi, secspday, nlvic
     use LandunitType     , only : lun                
     use landunit_varcon  , only : istsoil
     use clm_varctl       , only : use_hillslope_routing

     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)      :: bounds               
     integer                  , intent(in)      :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)      :: filter_hydrologyc(:) ! column filter for soil points
     type(soilstate_type)       , intent(in)    :: soilstate_inst
     type(soilhydrology_type)   , intent(inout) :: soilhydrology_inst
     type(waterstatebulk_type)  , intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type)   , intent(inout) :: waterfluxbulk_inst
     type(wateratm2lndbulk_type), intent(in)    :: wateratm2lndbulk_inst
     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'PerchedLateralFlowHillslope' ! subroutine name
     integer  :: c,fc,k,l,g                       ! indices
     real(r8) :: dtime                            ! land model time step (sec)
     real(r8) :: drainage_tot                     ! total amount of drainage to be removed from the column (mm/s)
     real(r8) :: drainage_layer                   ! amount of drainage to be removed from current layer (mm/s)
     real(r8) :: s_y                              ! specific yield (unitless)
     integer  :: k_frost(bounds%begc:bounds%endc) ! indices identifying frost table layer
     integer  :: k_perch(bounds%begc:bounds%endc) ! indices identifying perched water table layer
     real(r8) :: wtsub                            ! temporary variable
     real(r8) :: q_perch                          ! transmissivity (mm2/s)
     real(r8) :: q_perch_max                      ! baseflow coefficient
     real(r8) :: stream_water_depth               ! depth of water in stream channel (m)
     real(r8) :: stream_channel_depth             ! depth of stream channel (m)

     real(r8) :: transmis                         ! transmissivity (m2/s)
     real(r8) :: head_gradient                    ! head gradient (m/m)
     real(r8), parameter :: k_anisotropic = 1._r8 ! anisotropy factor
     integer  :: c0, c_src, c_dst                 ! indices
     real(r8) :: qflx_drain_perched_vol(bounds%begc:bounds%endc)   ! volumetric lateral subsurface flow through active layer [m3/s]
     real(r8) :: qflx_drain_perched_out(bounds%begc:bounds%endc)   ! lateral subsurface flow through active layer [mm/s]

     associate(                                                            & 
          nbedrock           =>    col%nbedrock                          , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          bsw                =>    soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"                        
          hksat              =>    soilstate_inst%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_inst%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)                       
          watsat             =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  

          frost_table        =>    soilhydrology_inst%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)                             
          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)                             
          zwt_perched        =>    soilhydrology_inst%zwt_perched_col    , & ! Input:  [real(r8) (:)   ] perched water table depth (m)                     
          tdepth             =>    wateratm2lndbulk_inst%tdepth_grc      , & ! Input:  [real(r8) (:)   ]  depth of water in tributary channels (m)
          tdepth_bankfull    =>    wateratm2lndbulk_inst%tdepthmax_grc   , & ! Input:  [real(r8) (:)   ]  bankfull depth of tributary channels (m)
          stream_water_volume =>    waterstatebulk_inst%stream_water_volume_lun , & ! Input:  [real(r8) (:)   ] stream water volume (m3)


          qflx_drain_perched =>    waterfluxbulk_inst%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ] perched wt sub-surface runoff (mm H2O /s)         

          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)                                
          )

       ! Get time step

       dtime = get_step_size_real()

       ! locate frost table and perched water table
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          k_frost(c) = nbedrock(c)
          k_perch(c) = nbedrock(c)
          do k = 1,nbedrock(c)
             if (frost_table(c) >= zi(c,k-1) .and. frost_table(c) < zi(c,k)) then
                k_frost(c) = k
                exit
             endif
          enddo

          do k = 1,nbedrock(c)
             if (zwt_perched(c) >= zi(c,k-1) .and. zwt_perched(c) < zi(c,k)) then
                k_perch(c) = k
                exit
             endif
          enddo
       enddo

       ! compute drainage from perched saturated region
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          l = col%landunit(c)
          g = col%gridcell(c)
          qflx_drain_perched(c)     = 0._r8
          qflx_drain_perched_out(c) = 0._r8
          qflx_drain_perched_vol(c) = 0._r8

          if (frost_table(c) > zwt_perched(c)) then
             ! Hillslope columns
             if (col%is_hillslope_column(c) .and. col%active(c)) then

                ! calculate head gradient

                if (head_gradient_method == kinematic) then
                   ! kinematic wave approximation
                   head_gradient = col%hill_slope(c)
                else if (head_gradient_method == darcy) then
                   ! darcy's law 
                   if (col%cold(c) /= ispval) then
                      head_gradient = (col%hill_elev(c)-zwt_perched(c)) &
                           - (col%hill_elev(col%cold(c))-zwt_perched(col%cold(c)))
                      head_gradient = head_gradient / (col%hill_distance(c) - col%hill_distance(col%cold(c)))
                   else
                      if (use_hillslope_routing) then
                         stream_water_depth = stream_water_volume(l) &
                              /lun%stream_channel_length(l)/lun%stream_channel_width(l)
                         stream_channel_depth = lun%stream_channel_depth(l)
                      else
                         stream_water_depth = tdepth(g)
                         stream_channel_depth = tdepth_bankfull(g)
                      endif

                      ! flow between channel and lowest column
                      ! bankfull height is defined to be zero
                      head_gradient = (col%hill_elev(c)-zwt_perched(c)) &
                           ! ignore overbankfull storage
                           - max(min((stream_water_depth - stream_channel_depth),0._r8), &
                           (col%hill_elev(c)-frost_table(c)))

                      head_gradient = head_gradient / (col%hill_distance(c))

                      ! head_gradient cannot be negative when channel is empty
                      if (stream_water_depth <= 0._r8) then
                         head_gradient = max(head_gradient, 0._r8)
                      endif
                   endif
                else                 
                   call endrun(msg="head_gradient_method must be kinematic or darcy"//errmsg(sourcefile, __LINE__))
                endif

                ! Determine source and destination columns
                if (head_gradient >= 0._r8) then
                   c_src = c
                   c_dst = col%cold(c)
                else
                   c_src = col%cold(c)
                   c_dst = c
                endif

                ! Calculate transmissivity of source column
                transmis = 0._r8

                if (transmissivity_method == layersum) then
                   if (head_gradient_method == kinematic) then
                      if(k_perch(c_src) < k_frost(c_src)) then
                         do k = k_perch(c_src), k_frost(c_src)-1
                            if(k == k_perch(c_src)) then
                               transmis = transmis + 1.e-3_r8*hksat(c_src,k)*(zi(c_src,k) - zwt_perched(c_src))
                            else
                               transmis = transmis + 1.e-3_r8*hksat(c_src,k)*dz(c_src,k)
                            endif
                         enddo
                      endif
                   else if (head_gradient_method == darcy) then
                      if(c_src == ispval) then
                         ! lowland, losing stream (c_src == ispval)
                         ! use hksat of c_dst for transmissivity
                         transmis = (1.e-3_r8*hksat(c,k_perch(c_dst)))*stream_water_depth
                      else
                         ! if k_perch equals k_frost, no perched saturated zone exists
                         if(k_perch(c_src) < k_frost(c_src)) then
                            do k = k_perch(c_src), k_frost(c_src)-1
                               if(k == k_perch(c_src)) then
                                  transmis = transmis + 1.e-3_r8*hksat(c_src,k)*(zi(c_src,k) - zwt_perched(c_src))
                               else
                                  if(c_dst == ispval) then 
                                     ! lowland, gaining stream
                                     ! only include layers above stream channel bottom
                                     if ((col%hill_elev(c_src)-z(c_src,k)) > (-stream_channel_depth)) then
                                        
                                        transmis = transmis + 1.e-3_r8*hksat(c_src,k)*dz(c_src,k)
                                     endif
                                  else
                                     ! uplands
                                     ! only include layers above dst water table elevation
                                     if ((col%hill_elev(c_src)-z(c_src,k)) > (col%hill_elev(c_dst) - zwt_perched(c_dst))) then
                                        
                                        transmis = transmis + 1.e-3_r8*hksat(c_src,k)*dz(c_src,k)
                                     endif
                                  endif
                               endif
                            enddo
                         endif
                      endif
                   endif
                else if (transmissivity_method == uniform_transmissivity) then
                   ! constant conductivity based on shallowest saturated layer hydraulic conductivity
                   transmis = (1.e-3_r8*hksat(c_src,k_perch(c_src))) &
                        *(zi(c_src,k_frost(c_src)) - zwt_perched(c_src) )
                endif

                ! adjust by 'anisotropy factor'
                transmis = k_anisotropic*transmis

                qflx_drain_perched_vol(c) = transmis*col%hill_width(c)*head_gradient
                qflx_drain_perched_out(c) = 1.e3_r8*(qflx_drain_perched_vol(c)/col%hill_area(c))
 
             else
                ! Non-hillslope columns
                ! specify maximum drainage rate
                q_perch_max = params_inst%perched_baseflow_scalar &
                     * sin(col%topo_slope(c) * (rpi/180._r8))

                wtsub = 0._r8
                q_perch = 0._r8
                ! this should be consistent with hillslope and k_perch=k_frost means no
                ! saturated zone; should probably change q_perch to tranmis and change
                ! units and q_perch_max
                do k = k_perch(c), k_frost(c)-1
                   q_perch = q_perch + hksat(c,k)*dz(c,k)
                   wtsub = wtsub + dz(c,k)
                end do
                if (wtsub > 0._r8) q_perch = q_perch/wtsub
                
                qflx_drain_perched_out(c) = q_perch_max * q_perch &
                     *(frost_table(c) - zwt_perched(c))
             endif
          endif

       enddo
             
       ! compute net drainage from perched saturated region
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          ! drainage-out
          qflx_drain_perched(c) = qflx_drain_perched(c) + qflx_drain_perched_out(c)
          if (col%is_hillslope_column(c) .and. col%active(c)) then
             ! drainage-in
             if (col%cold(c) /= ispval) then
                qflx_drain_perched(col%cold(c)) = &
                     qflx_drain_perched(col%cold(c)) - &
                     1.e3_r8*(qflx_drain_perched_vol(c))/col%hill_area(col%cold(c))
             endif
          endif
       enddo

       ! remove drainage from soil moisture storage
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          
          ! remove drainage from perched saturated layers
          drainage_tot = qflx_drain_perched(c) * dtime
          ! ignore frozen layer (k_frost)
          do k = k_perch(c), k_frost(c)-1

             s_y = watsat(c,k) &
                  * ( 1. - (1.+1.e3*zwt_perched(c)/sucsat(c,k))**(-1./bsw(c,k)))
             s_y=max(s_y,params_inst%aq_sp_yield_min)
             if (k==k_perch(c)) then
                drainage_layer=min(drainage_tot,(s_y*(zi(c,k) - zwt_perched(c))*1.e3))
             else
                drainage_layer=min(drainage_tot,(s_y*(dz(c,k))*1.e3))
             endif
             
             drainage_layer=max(drainage_layer,0._r8)
             drainage_tot = drainage_tot - drainage_layer
             h2osoi_liq(c,k) = h2osoi_liq(c,k) - drainage_layer
             
          enddo

          ! if drainage_tot is greater than available water
          ! (above frost table), then decrease qflx_drain_perched
          ! by residual amount for water balance
          qflx_drain_perched(c) = qflx_drain_perched(c) - drainage_tot/dtime
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
   subroutine SubsurfaceLateralFlow(bounds,  & 
        num_hydrologyc, filter_hydrologyc,  &
        num_urbanc, filter_urbanc,soilhydrology_inst, soilstate_inst, &
        waterstatebulk_inst, waterfluxbulk_inst, wateratm2lndbulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage
     !
     ! !USES:
     use clm_time_manager , only : get_step_size
     use clm_varpar       , only : nlevsoi, nlevgrnd, nlayer, nlayert
     use clm_varctl       , only : nhillslope
     use clm_varcon       , only : pondmx, watmin,rpi, secspday
     use column_varcon    , only : icol_road_perv
     use abortutils       , only : endrun
     use GridcellType     , only : grc  
     use landunit_varcon  , only : istsoil, istcrop
     use clm_varctl       , only : use_hillslope_routing

     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilstate_type)     , intent(in)    :: soilstate_inst
     type(wateratm2lndbulk_type)       , intent(in)    :: wateratm2lndbulk_inst
     type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
     type(waterstatebulk_type), intent(inout) :: waterstatebulk_inst
     type(waterfluxbulk_type) , intent(inout) :: waterfluxbulk_inst

     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'SubsurfaceLateralFlow' ! subroutine name
     integer  :: c,j,fc,i,l,g                            ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: drainage(bounds%begc:bounds%endc)       ! subsurface drainage (mm/s)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: xs1(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer 1 (mm)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: ice_imped_col(bounds%begc:bounds%endc)  ! column average hydraulic conductivity reduction due to presence of soil ice (-)
     real(r8) :: ice_imped(bounds%begc:bounds%endc,1:nlevsoi) ! hydraulic conductivity reduction due to presence of soil ice (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: h2osoi_vol                              ! volumetric water content (mm3/mm3)
     real(r8) :: drainage_tot                            ! total drainage to be removed from column (mm)
     real(r8) :: drainage_layer                          ! drainage to be removed from current layer (mm)
     real(r8) :: s_y                                     ! specific yield (unitless)
     real(r8) :: vol_ice                          ! volumetric ice content (mm3/mm3)
     logical, parameter :: no_lateral_flow = .false.    ! flag for testing
     real(r8) :: transmis                         ! transmissivity (m2/s)
     real(r8) :: head_gradient                    ! hydraulic head gradient (m/m)
     real(r8) :: stream_water_depth               ! depth of water in stream channel (m)
     real(r8) :: stream_channel_depth             ! depth of stream channel (m)
     real(r8) :: available_stream_water           ! stream water (m3)
     real(r8), parameter :: n_baseflow = 1        ! drainage power law exponent
     real(r8), parameter :: k_anisotropic = 1._r8 ! anisotropy scalar
     real(r8) :: qflx_latflow_out_vol(bounds%begc:bounds%endc) ! volumetric lateral flow (m3/s)
     real(r8) :: qflx_net_latflow(bounds%begc:bounds%endc)     ! net lateral flow in column (mm/s)
     real(r8) :: qflx_latflow_avg(bounds%begc:bounds%endc)     ! average lateral flow (mm/s)
     real(r8) :: larea                            ! area of hillslope in landunit
     integer  :: c0, c_src, c_dst                 ! indices
     
     !-----------------------------------------------------------------------

     associate(                                                            & 
          nbedrock           =>    col%nbedrock                          , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)           
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
          qflx_latflow_out   =>    waterfluxbulk_inst%qflx_latflow_out_col, & ! Output: [real(r8) (:)   ] lateral saturated outflow (mm/s)
          qflx_latflow_in    =>    waterfluxbulk_inst%qflx_latflow_in_col, & ! Output: [real(r8) (:)   ]  lateral saturated inflow (mm/s)
          volumetric_discharge =>  waterfluxbulk_inst%volumetric_discharge_col , & ! Output: [real(r8) (:)   ]  discharge from column (m3/s)

          tdepth             =>    wateratm2lndbulk_inst%tdepth_grc      , & ! Input:  [real(r8) (:)   ]  depth of water in tributary channels (m)
          tdepth_bankfull    =>    wateratm2lndbulk_inst%tdepthmax_grc   , & ! Input:  [real(r8) (:)   ]  bankfull depth of tributary channels (m)

          depth              =>    soilhydrology_inst%depth_col          , & ! Input:  [real(r8) (:,:) ] VIC soil depth                                   
          icefrac            =>    soilhydrology_inst%icefrac_col        , & ! Output: [real(r8) (:,:) ] fraction of ice in layer                         
          frost_table        =>    soilhydrology_inst%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)                             
          zwt                =>    soilhydrology_inst%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)                             
          stream_water_volume =>    waterstatebulk_inst%stream_water_volume_lun, & ! Input:  [real(r8) (:)   ] stream water volume (m3)
          
          qflx_snwcp_liq     =>    waterfluxbulk_inst%qflx_snwcp_liq_col     , & ! Output: [real(r8) (:)   ] excess rainfall due to snow capping (mm H2O /s) [+]
          qflx_ice_runoff_xs =>    waterfluxbulk_inst%qflx_ice_runoff_xs_col , & ! Output: [real(r8) (:)   ] solid runoff from excess ice in soil (mm H2O /s) [+]
          qflx_drain         =>    waterfluxbulk_inst%qflx_drain_col         , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)                    
          qflx_qrgwl         =>    waterfluxbulk_inst%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ] qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
          qflx_rsub_sat      =>    waterfluxbulk_inst%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ] soil saturation excess [mm h2o/s]                 
          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)                                
          )

       ! Get time step

       dtime = get_step_size_real()

       ! Convert layer thicknesses from m to mm

       do j = 1,nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             dzmm(c,j) = dz(c,j)*1.e3_r8

             vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             icefrac(c,j) = min(1._r8,vol_ice/watsat(c,j))
             ice_imped(c,j)=10._r8**(-params_inst%e_ice*icefrac(c,j))
          end do
       end do

       ! Initial set

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          qflx_drain(c)    = 0._r8 
          qflx_rsub_sat(c) = 0._r8
          drainage(c)      = 0._r8
          qflx_latflow_in(c) = 0._r8
          qflx_latflow_out(c) = 0._r8
          qflx_net_latflow(c) = 0._r8
          volumetric_discharge(c)       = 0._r8
          qflx_latflow_out_vol(c) = 0._r8
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

      ! Calculate ice impedance factor (after jwt calculated)
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         dzsum = 0._r8
         icefracsum = 0._r8
         do j = max(jwt(c),1), nlevsoi
            dzsum  = dzsum + dzmm(c,j)
            icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
         end do
         ice_imped_col(c)=10._r8**(-params_inst%e_ice*(icefracsum/dzsum))
      enddo

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         l = col%landunit(c)
         g = col%gridcell(c)
         ! Hillslope columns
         if (col%is_hillslope_column(c) .and. col%active(c)) then

            ! method for calculating head gradient
            if (head_gradient_method == kinematic) then
               head_gradient = col%hill_slope(c)
            else if (head_gradient_method == darcy) then
               if (col%cold(c) /= ispval) then
                  head_gradient = (col%hill_elev(c)-zwt(c)) &
                       - (col%hill_elev(col%cold(c))-zwt(col%cold(c)))
                  head_gradient = head_gradient / (col%hill_distance(c) - col%hill_distance(col%cold(c)))
               else
                  if (use_hillslope_routing) then
                     stream_water_depth = stream_water_volume(l) &
                          /lun%stream_channel_length(l)/lun%stream_channel_width(l)
                     stream_channel_depth = lun%stream_channel_depth(l)
                  else
                     stream_water_depth = tdepth(g)
                     stream_channel_depth = tdepth_bankfull(g)
                  endif

                  ! flow between channel and lowest column
                  ! bankfull height is defined to be zero
                  head_gradient = (col%hill_elev(c)-zwt(c)) &
                       ! ignore overbankfull storage
                       - min((stream_water_depth - stream_channel_depth),0._r8)

                  head_gradient = head_gradient / (col%hill_distance(c))
                  ! head_gradient cannot be negative when channel is empty
                  if (stream_water_depth <= 0._r8) then
                     head_gradient = max(head_gradient, 0._r8)
                  endif
                  ! add vertical drainage for losing streams
                  ! (this could be a separate term from lateral flow...)
                  if (head_gradient < 0._r8) then
                     ! head_gradient = head_gradient - 1._r8
                     ! adjust lateral gradient w/ k_anisotropic
                     head_gradient = head_gradient - 1._r8/k_anisotropic
                  endif
               endif
            else
               call endrun(msg="head_gradient_method must be kinematic or darcy"//errmsg(sourcefile, __LINE__))  
            end if

            !scs: in cases of bad data, where hand differences in 
            ! adjacent bins are very large, cap maximum head_gradient
            ! should a warning be used instead?
            head_gradient = min(max(head_gradient,-2._r8),2._r8)
            
            ! Determine source and destination columns
            if (head_gradient >= 0._r8) then
               c_src = c
               c_dst = col%cold(c)
            else
               c_src = col%cold(c)
               c_dst = c
            endif
           
            ! Calculate transmissivity of source column
            transmis = 0._r8
            if(c_src /= ispval) then 
               ! transmissivity non-zero only when saturated conditions exist
               if(zwt(c_src) <= zi(c_src,nbedrock(c_src))) then 
                  ! sum of layer transmissivities
                  if (transmissivity_method == layersum) then
                     do j = jwt(c_src)+1, nbedrock(c_src)
                        if(j == jwt(c_src)+1) then
                           transmis = transmis + 1.e-3_r8*ice_imped(c_src,j)*hksat(c_src,j)*(zi(c_src,j) - zwt(c_src))
                        else
                           if(c_dst == ispval) then 
                              ! lowland, gaining stream
                              ! only include layers above stream channel bottom
                              if ((col%hill_elev(c_src)-z(c_src,j)) > (-stream_channel_depth)) then
                                 
                                 transmis = transmis + 1.e-3_r8*ice_imped(c_src,j)*hksat(c_src,j)*dz(c_src,j)
                              endif
                           else
                              ! uplands
                              if ((col%hill_elev(c_src)-z(c_src,j)) > (col%hill_elev(c_dst) - zwt(c_dst))) then
                                 transmis = transmis + 1.e-3_r8*ice_imped(c_src,j)*hksat(c_src,j)*dz(c_src,j)
                              endif
                           endif
                        endif
                     end do
                  ! constant conductivity based on shallowest saturated layer hk
                  else if (transmissivity_method == uniform_transmissivity) then
                     transmis = (1.e-3_r8*ice_imped(c_src,jwt(c_src)+1)*hksat(c_src,jwt(c_src)+1)) &
                          *(zi(c_src,nbedrock(c_src)) - zwt(c_src) )
                  else
                     call endrun(msg="transmissivity_method must be LayerSum or Uniform"//errmsg(sourcefile, __LINE__))
                  endif
               endif
            else
               ! transmissivity of losing stream (c_src == ispval)
               transmis = (1.e-3_r8*ice_imped(c,jwt(c)+1)*hksat(c,jwt(c)+1))*stream_water_depth
            endif
            ! adjust transmissivity by 'anisotropy factor'
            transmis = k_anisotropic*transmis

            ! the qflx_latflow_out_vol calculations use the
            ! transmissivity to determine whether saturated flow
            ! conditions exist, b/c gradients will be nonzero
            ! even when no saturated layers are present
            !          qflx_latflow_out_vol(c) = ice_imped(c)*transmis*col%hill_width(c)*head_gradient
            ! include ice impedance in transmissivity
            qflx_latflow_out_vol(c) = transmis*col%hill_width(c)*head_gradient

            ! When head gradient is negative (losing stream channel), 
            ! limit outflow by available stream channel water
            if (use_hillslope_routing .and. (qflx_latflow_out_vol(c) < 0._r8)) then
               available_stream_water = stream_water_volume(l)/lun%stream_channel_number(l)/nhillslope
               if(abs(qflx_latflow_out_vol(c))*dtime > available_stream_water) then
                  qflx_latflow_out_vol(c) = -available_stream_water/dtime
               endif
            endif

            ! volumetric_discharge from lowest column is qflx_latflow_out_vol
            ! scaled by total area of column in gridcell divided by column area
            if (col%cold(c) == ispval) then
               volumetric_discharge(c) = qflx_latflow_out_vol(c) &
                    *(grc%area(g)*1.e6_r8*col%wtgcell(c)/col%hill_area(c))
            endif

            ! convert volumetric flow to equivalent flux
            qflx_latflow_out(c) = 1.e3_r8*qflx_latflow_out_vol(c)/col%hill_area(c)

            ! hilltop column has no inflow
            if (col%colu(c) == ispval) then
               qflx_latflow_in(c) = 0._r8
            endif

            ! current outflow is inflow to downhill column normalized by downhill area
            if (col%cold(c) /= ispval) then
               qflx_latflow_in(col%cold(c)) = qflx_latflow_in(col%cold(c)) + &
                    1.e3_r8*qflx_latflow_out_vol(c)/col%hill_area(col%cold(c))
            endif

         else
            ! Non-hillslope columns
            ! baseflow is power law expression relative to bedrock layer
            if(zwt(c) <= zi(c,nbedrock(c))) then
               qflx_latflow_out(c) = ice_imped_col(c) * baseflow_scalar &
                    * tan(rpi/180._r8*col%topo_slope(c))* &
                    (zi(c,nbedrock(c)) - zwt(c))**(params_inst%n_baseflow)
            endif
            ! convert flux to volumetric flow
            qflx_latflow_out_vol(c) = 1.e-3_r8*qflx_latflow_out(c)*(grc%area(g)*1.e6_r8*col%wtgcell(c))
            volumetric_discharge(c) = qflx_latflow_out_vol(c)
         endif
      enddo

      ! recalculate average flux for no-lateral flow case
      if(no_lateral_flow) then
         if (head_gradient_method /= kinematic) then
            call endrun(msg="head_gradient_method must be kinematic for no_lateral_flow = .true.! "//errmsg(sourcefile, __LINE__))
         endif
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (col%is_hillslope_column(c) .and. col%active(c)) then
               l = col%landunit(c)
               !need to sum all columns w/ same hillslope id for each column
               qflx_latflow_avg(c) = 0._r8
               larea = 0._r8
               do c0 = lun%coli(l), lun%colf(l)
                  if(col%hillslope_ndx(c0) == col%hillslope_ndx(c)) then
                     qflx_latflow_avg(c) = qflx_latflow_avg(c) + qflx_latflow_out_vol(c0)
                     larea = larea + col%hill_area(c0)
                  endif
               enddo
               qflx_latflow_avg(c) = 1.e3_r8*qflx_latflow_avg(c)/larea
            else
               qflx_latflow_avg(c) = qflx_latflow_out(c)
            endif
         enddo
      endif
         
      !-- Topographic runoff  -------------------------
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         
         ! net lateral flow (positive out)
         qflx_net_latflow(c) = qflx_latflow_out(c) - qflx_latflow_in(c)
         if(no_lateral_flow) then
            qflx_net_latflow(c) = qflx_latflow_avg(c)
         endif
         
         !@@
         ! baseflow 
         if(zwt(c) <= zi(c,nbedrock(c))) then 
            ! apply net lateral flow here
            drainage(c) = qflx_net_latflow(c)
         else
            drainage(c) = 0._r8
         endif
         
         !--  Now remove water via drainage
         drainage_tot = - drainage(c) * dtime
         
         if(drainage_tot > 0.) then !rising water table
            do j = jwt(c)+1,1,-1

               ! ensure water is not added to frozen layers
               if (zi(c,j) < frost_table(c)) then 
                  ! analytical expression for specific yield
                  s_y = watsat(c,j) &
                       * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                  s_y=max(s_y,params_inst%aq_sp_yield_min)

                  drainage_layer=min(drainage_tot,(s_y*dz(c,j)*1.e3))

                  drainage_layer=max(drainage_layer,0._r8)
                  h2osoi_liq(c,j) = h2osoi_liq(c,j) + drainage_layer

                  drainage_tot = drainage_tot - drainage_layer

                  if (drainage_tot <= 0.) then 
                     zwt(c) = zwt(c) - drainage_layer/s_y/1000._r8
                     exit
                  else
                     zwt(c) = zi(c,j-1)
                  endif
               endif
                  
            enddo
            
            !--  remove residual drainage  --------------------------------
            h2osfc(c) = h2osfc(c) + drainage_tot
                    
          else ! deepening water table
             do j = jwt(c)+1, nbedrock(c)
                ! analytical expression for specific yield
                s_y = watsat(c,j) &
                     * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                s_y=max(s_y,params_inst%aq_sp_yield_min)
                
                drainage_layer=max(drainage_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                drainage_layer=min(drainage_layer,0._r8)
                h2osoi_liq(c,j) = h2osoi_liq(c,j) + drainage_layer

                drainage_tot = drainage_tot - drainage_layer
                   
                if (drainage_tot >= 0.) then 
                   zwt(c) = zwt(c) - drainage_layer/s_y/1000._r8
                   exit
                else
                   zwt(c) = zi(c,j)
                endif
             enddo
             
             !--  remove residual drainage  -----------------------
             ! make sure no extra water removed from soil column
             drainage(c) = drainage(c) + drainage_tot/dtime
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
          qflx_drain(c) = qflx_rsub_sat(c) + drainage(c)

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

   end subroutine SubsurfaceLateralFlow

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
     integer  :: c ,j,fc,i                                       ! indices
     real(r8) :: dtime                                           ! land model time step (sec)
     real(r8) :: qflx_solidevap_from_top_layer_save              ! temporary
     integer  :: num_modifiedc                                   ! number of columns in filter_modifiedc
     integer  :: filter_modifiedc(bounds%endc-bounds%begc+1)     ! column filter of points modified in this subroutine
     real(r8) :: h2osoi_ice_before_evap(bounds%begc:bounds%endc) ! h2osoi_ice in layer 1 before applying solidevap
     !-----------------------------------------------------------------------

     associate(                                                            & 
          snl                =>    col%snl                               , & ! Input:  [integer  (:)   ]  number of snow layers                              
          h2osoi_liq         =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
          frac_h2osfc        =>    waterdiagnosticbulk_inst%frac_h2osfc_col       , & ! Input:  [real(r8) (:)   ]                                                    
          qflx_liqdew_to_top_layer      => waterfluxbulk_inst%qflx_liqdew_to_top_layer_col  , & ! Input:  [real(r8) (:)   ]  rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s) [+]    
          qflx_soliddew_to_top_layer    => waterfluxbulk_inst%qflx_soliddew_to_top_layer_col, & ! Input:  [real(r8) (:)   ]  rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s) [+]      
          qflx_ev_snow                  => waterfluxbulk_inst%qflx_ev_snow_col    , & ! In/Out: [real(r8) (:)   ]  evaporation flux from snow (mm H2O/s) [+ to atm]
          qflx_solidevap_from_top_layer => waterfluxbulk_inst%qflx_solidevap_from_top_layer_col & ! Output: [real(r8) (:)   ]  rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s) [+]   
          )

       ! Get time step

       dtime = get_step_size_real()
       num_modifiedc = 0

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Renew the ice and liquid mass due to condensation

          if (snl(c)+1 >= 1) then
             num_modifiedc = num_modifiedc + 1
             filter_modifiedc(num_modifiedc) = c

             ! make consistent with how evap_grnd removed in infiltration
             h2osoi_liq(c,1) = h2osoi_liq(c,1) + (1._r8 - frac_h2osfc(c))*qflx_liqdew_to_top_layer(c) * dtime
             h2osoi_ice(c,1) = h2osoi_ice(c,1) + (1._r8 - frac_h2osfc(c))*qflx_soliddew_to_top_layer(c) * dtime
             h2osoi_ice_before_evap(c) = h2osoi_ice(c,1)
             h2osoi_ice(c,1) = h2osoi_ice(c,1) - (1._r8 - frac_h2osfc(c)) * qflx_solidevap_from_top_layer(c) * dtime
          end if

       end do


       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          ! Renew the ice and liquid mass due to condensation for urban roof and impervious road

          if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
             if (snl(c)+1 >= 1) then
                num_modifiedc = num_modifiedc + 1
                filter_modifiedc(num_modifiedc) = c

                h2osoi_liq(c,1) = h2osoi_liq(c,1) + qflx_liqdew_to_top_layer(c) * dtime
                h2osoi_ice(c,1) = h2osoi_ice(c,1) + (qflx_soliddew_to_top_layer(c) * dtime)
                h2osoi_ice_before_evap(c) = h2osoi_ice(c,1)
                h2osoi_ice(c,1) = h2osoi_ice(c,1) - (qflx_solidevap_from_top_layer(c) * dtime)
             end if
          end if

       end do

       call truncate_small_values( &
            num_f              = num_modifiedc, &
            filter_f           = filter_modifiedc, &
            lb                 = bounds%begc, &
            ub                 = bounds%endc, &
            data_baseline      = h2osoi_ice_before_evap(bounds%begc:bounds%endc), &
            data               = h2osoi_ice(bounds%begc:bounds%endc, 1), &
            custom_rel_epsilon = tolerance)

       do fc = 1, num_modifiedc
          c = filter_modifiedc(fc)

          if (h2osoi_ice(c,1) < 0._r8) then
             write(iulog,*) "ERROR: In RenewCondensation, h2osoi_ice has gone significantly negative"
             write(iulog,*) "c = ", c
             write(iulog,*) "h2osoi_ice_before_evap = ", h2osoi_ice_before_evap(c)
             write(iulog,*) "h2osoi_ice(c,1)        = ", h2osoi_ice(c,1)
             write(iulog,*) "qflx_solidevap_from_top_layer*dtime = ", qflx_solidevap_from_top_layer(c)*dtime
             call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                  msg="In RenewCondensation, h2osoi_ice has gone significantly negative")
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

     SHR_ASSERT_ALL_FL((ubound(qflx_gw_demand) == [bounds%endc]), sourcefile, __LINE__)
     SHR_ASSERT_ALL_FL((ubound(qflx_gw_uncon_irrig_lyr) == [bounds%endc, nlevsoi]), sourcefile, __LINE__)
     SHR_ASSERT_ALL_FL((ubound(qflx_gw_con_irrig) == [bounds%endc]), sourcefile, __LINE__)

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

       dtime = get_step_size_real()

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
             
             call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                  msg="negative groundwater irrigation demand! "//errmsg(sourcefile, __LINE__))
             
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
                s_y=max(s_y, params_inst%aq_sp_yield_min)

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

     dtime = get_step_size_real()

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
