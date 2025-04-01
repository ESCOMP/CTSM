module SoilTemperatureMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates snow and soil temperatures including phase change
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use shr_infnan_mod          , only : nan => shr_infnan_nan, assignment(=)
  use decompMod               , only : bounds_type
  use abortutils              , only : endrun
  use perf_mod                , only : t_startf, t_stopf
  use clm_varctl              , only : iulog
  use UrbanParamsType         , only : urbanparams_type
  use UrbanTimeVarType        , only : urbantv_type
  use atm2lndType             , only : atm2lnd_type
  use CanopyStateType         , only : canopystate_type
  use WaterFluxBulkType       , only : waterfluxbulk_type
  use WaterStateBulkType      , only : waterstatebulk_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use SolarAbsorbedType       , only : solarabs_type
  use SoilStateType           , only : soilstate_type
  use EnergyFluxType          , only : energyflux_type
  use TemperatureType         , only : temperature_type
  use LandunitType            , only : lun
  use ColumnType              , only : col
  use PatchType               , only : patch
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilTemperature 
  !
  !    -> SoilTemperature:      soil/snow and ground temperatures
  !          -> SoilTermProp    thermal conductivities and heat capacities
  !          -> Tridiagonal     tridiagonal matrix solution
  !          -> PhaseChange     phase change of liquid/ice contents
  !
  ! (1) Snow and soil temperatures
  !     o The volumetric heat capacity is calculated as a linear combination
  !       in terms of the volumetric fraction of the constituent phases.
  !     o The thermal conductivity of soil is computed from
  !       the algorithm of Johansen (as reported by Farouki 1981), and the
  !       conductivity of snow is from the formulation used in
  !       Sturm (1997) or Jordan (1991) p. 18 depending on namelist option.
  !     o Boundary conditions:
  !       F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
  !     o Soil / snow temperature is predicted from heat conduction
  !       in 10 soil layers and up to nlevsno snow layers.
  !       The thermal conductivities at the interfaces between two
  !       neighboring layers (j, j+1) are derived from an assumption that
  !       the flux across the interface is equal to that from the node j
  !       to the interface and the flux from the interface to the node j+1.
  !       The equation is solved using the Crank-Nicholson method and
  !       results in a tridiagonal system equation.
  ! (2) Phase change
  !
  ! The following is only public for the sake of unit testing; it should not be called
  ! directly by CLM code outside this module
  public :: ComputeGroundHeatFluxAndDeriv       ! Computes G and dG/dT on surface of standing water, snow and soil
  public :: ComputeHeatDiffFluxAndFactor        ! Heat diffusion at layer interface and factor used in setting up of banded matrix
  public :: SetRHSVec                           ! Sets up the RHS vector for the numerical solution of temperature for snow/standing-water/soil
  public :: SetRHSVec_Snow                      ! Sets up the RHS vector corresponding to snow layers for all columns
  public :: SetRHSVec_Soil                      ! Sets up the RHS vector corresponding to soil layers for all columns
  public :: SetRHSVec_StandingSurfaceWater      ! Sets up the RHS vector corresponding to standing water layers for all columns
  public :: SetMatrix                           ! Sets up the matrix for the numerical solution of temperature for snow/standing-water/soil
  public :: AssembleMatrixFromSubmatrices       ! Assemble the full matrix from submatrices.
  public :: SetMatrix_Snow                      ! Set up the matrix entries corresponding to snow layers for all columns
  public :: SetMatrix_Soil                      ! Set up the matrix entries corresponding to soil layers for all columns
  public :: SetMatrix_StandingSurfaceWater      ! Set up the matrix entries corresponding to standing surface water for all columns
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilThermProp       ! Set therm conduct. and heat cap of snow/soil layers
  private :: PhaseChangeH2osfc   ! When surface water freezes move ice to bottom snow layer
  private :: PhaseChange_beta    ! Calculation of the phase change within snow and soil layers
  private :: BuildingHAC         ! Building Heating and Cooling for simpler method (introduced in CLM4.5)

  real(r8), private, parameter :: thin_sfclayer = 1.0e-6_r8   ! Threshold for thin surface layer
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilTemperature(bounds, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, &
       num_nolakep, filter_nolakep, num_nolakec, filter_nolakec, &
       atm2lnd_inst, urbanparams_inst, canopystate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst,&
       solarabs_inst, soilstate_inst, energyflux_inst,  temperature_inst, urbantv_inst)
    !
    ! !DESCRIPTION:
    ! Snow and soil temperatures including phase change
    ! o The volumetric heat capacity is calculated as a linear combination
    !   in terms of the volumetric fraction of the constituent phases.
    ! o The thermal conductivity of soil is computed from
    !   the algorithm of Johansen (as reported by Farouki 1981), and the
    !   conductivity of snow is from the formulation used in
    !   Sturm (1997) or Jordan (1991) p. 18 depending on namelist option.
    ! o Boundary conditions:
    !   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
    ! o Soil / snow temperature is predicted from heat conduction
    !   in 10 soil layers and up to nlevsno snow layers.
    !   The thermal conductivities at the interfaces between two
    !   neighboring layers (j, j+1) are derived from an assumption that
    !   the flux across the interface is equal to that from the node j
    !   to the interface and the flux from the interface to the node j+1.
    !   The equation is solved using the Crank-Nicholson method and
    !   results in a tridiagonal system equation.
    !
    ! !USES:
    use clm_time_manager         , only : get_step_size_real
    use clm_varpar               , only : nlevsno, nlevgrnd, nlevurb, nlevmaxurbgrnd
    use clm_varctl               , only : iulog, use_excess_ice
    use clm_varcon               , only : cnfac, cpice, cpliq, denh2o, denice
    use landunit_varcon          , only : istsoil, istcrop
    use column_varcon            , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv
    use BandDiagonalMod          , only : BandDiagonal
    use UrbanParamsType          , only : IsSimpleBuildTemp, IsProgBuildTemp
    use UrbBuildTempOleson2015Mod, only : BuildingTemperature
    !
    ! !ARGUMENTS:
    type(bounds_type)              ,  intent(in)    :: bounds
    integer                        ,  intent(in)    :: num_nolakep       ! number of non-lake points in patch filter
    integer                        ,  intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
    integer                        ,  intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
    integer                        ,  intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    integer                        ,  intent(in)    :: num_urbanl        ! number of urban landunits in clump
    integer                        ,  intent(in)    :: filter_urbanl(:)  ! urban landunit filter
    integer                        ,  intent(in)    :: num_urbanc        ! number of urban columns in clump
    integer                        ,  intent(in)    :: filter_urbanc(:)  ! urban column filter
    type(atm2lnd_type)             ,  intent(in)    :: atm2lnd_inst
    type(urbanparams_type)         ,  intent(in)    :: urbanparams_inst
    type(urbantv_type)             ,  intent(in)    :: urbantv_inst
    type(canopystate_type)         ,  intent(in)    :: canopystate_inst
    type(waterstatebulk_type)      ,  intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type) ,  intent(inout) :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)       ,  intent(inout) :: waterfluxbulk_inst
    type(soilstate_type)           ,  intent(inout) :: soilstate_inst
    type(solarabs_type)            ,  intent(inout) :: solarabs_inst
    type(energyflux_type)          ,  intent(inout) :: energyflux_inst
    type(temperature_type)         ,  intent(inout) :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l,g                                                  ! indices
    integer  :: fc, fp                                                   ! lake filtered column & patch indices
    integer  :: fl                                                       ! urban filtered landunit indices
    integer  :: jtop(bounds%begc:bounds%endc)                            ! top level at each column
    real(r8) :: dtime                                                    ! land model time step (sec)
    real(r8) :: cv (bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)   ! heat capacity [J/(m2 K)]
    real(r8) :: tk (bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)   ! thermal conductivity [W/(m K)]
    real(r8) :: fn (bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)   ! heat diffusion through the layer interface [W/m2]
    real(r8) :: fn1(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)   ! heat diffusion through the layer interface [W/m2]
    real(r8) :: dzm                                                      ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                      ! used in computing tridiagonal matrix
    real(r8) :: sabg_lyr_col(bounds%begc:bounds%endc,-nlevsno+1:1)       ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) :: eflx_gnet_top                                            ! net energy flux into surface layer, patch-level [W/m2]
    real(r8) :: hs_top(bounds%begc:bounds%endc)                          ! net energy flux into surface layer (col) [W/m2]
    logical  :: cool_on(bounds%begl:bounds%endl)                         ! is urban air conditioning on?
    logical  :: heat_on(bounds%begl:bounds%endl)                         ! is urban heating on?
    real(r8) :: fn_h2osfc(bounds%begc:bounds%endc)                       ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8) :: dz_h2osfc(bounds%begc:bounds%endc)                       ! height of standing surface water [m]
    integer, parameter :: nband=5
    real(r8) :: bmatrix(bounds%begc:bounds%endc,nband,-nlevsno:nlevmaxurbgrnd) ! banded matrix for numerical solution of temperature
    real(r8) :: tvector(bounds%begc:bounds%endc,-nlevsno:nlevmaxurbgrnd)       ! initial temperature solution [K]
    real(r8) :: rvector(bounds%begc:bounds%endc,-nlevsno:nlevmaxurbgrnd)       ! RHS vector for numerical solution of temperature
    real(r8) :: tk_h2osfc(bounds%begc:bounds%endc)                       ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8) :: dhsdT(bounds%begc:bounds%endc)                           ! temperature derivative of "hs" [col]
    real(r8) :: hs_soil(bounds%begc:bounds%endc)                         ! heat flux on soil [W/m2]
    real(r8) :: hs_top_snow(bounds%begc:bounds%endc)                     ! heat flux on top snow layer [W/m2]
    real(r8) :: hs_h2osfc(bounds%begc:bounds%endc)                       ! heat flux on standing water [W/m2]
    integer  :: jbot(bounds%begc:bounds%endc)                            ! bottom level at each column
    real(r8) :: dz_0(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)                ! original layer thickness [m] 
    real(r8) :: z_0(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)                 ! original layer depth [m]
    real(r8) :: zi_0(bounds%begc:bounds%endc,-nlevsno+0:nlevmaxurbgrnd)                ! original layer interface level bellow layer "z" [m]

    !-----------------------------------------------------------------------

    associate(                                                                &
         snl                     => col%snl                                 , & ! Input:  [integer  (:)   ]  number of snow layers                    
         zi                      => col%zi                                  , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m) 
         dz                      => col%dz                                  , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                       
         z                       => col%z                                   , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                   
         ctype                   => col%itype                               , & ! Input: [integer (:)    ]  column type

         
         t_building_max          => urbantv_inst%t_building_max             , & ! Input:  [real(r8) (:)   ]  maximum internal building air temperature [K]
         t_building_min          => urbanparams_inst%t_building_min         , & ! Input:  [real(r8) (:)   ]  minimum internal building air temperature [K]

         
         forc_lwrad              => atm2lnd_inst%forc_lwrad_downscaled_col  , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)

         
         frac_veg_nosno          => canopystate_inst%frac_veg_nosno_patch   , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]

         
         frac_sno_eff            => waterdiagnosticbulk_inst%frac_sno_eff_col        , & ! Input:  [real(r8) (:)   ]  eff. fraction of ground covered by snow (0 to 1)
         snow_depth              => waterdiagnosticbulk_inst%snow_depth_col          , & ! Input:  [real(r8) (:)   ]  snow height (m)                         
         h2osfc                  => waterstatebulk_inst%h2osfc_col                   , & ! Input:  [real(r8) (:)   ]  surface water (mm)                      
         excess_ice              => waterstatebulk_inst%excess_ice_col               , & ! Input:  [real(r8) (:,:) ]  excess ice (kg/m2) (new) (1:nlevgrnd)
         frac_h2osfc             => waterdiagnosticbulk_inst%frac_h2osfc_col         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)

         
         sabg_soil               => solarabs_inst%sabg_soil_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
         sabg_snow               => solarabs_inst%sabg_snow_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
         sabg_chk                => solarabs_inst%sabg_chk_patch            , & ! Output: [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
         sabg_lyr                => solarabs_inst%sabg_lyr_patch            , & ! Input:  [real(r8) (:,:) ]  absorbed solar radiation (pft,lyr) [W/m2]
         sabg                    => solarabs_inst%sabg_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)

         
         htvp                    => energyflux_inst%htvp_col                , & ! Input:  [real(r8) (:)   ]  latent heat of vapor of water (or sublimation) [j/kg]
         cgrnd                   => energyflux_inst%cgrnd_patch             , & ! Input:  [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]
         dlrad                   => energyflux_inst%dlrad_patch             , & ! Input:  [real(r8) (:)   ]  downward longwave radiation blow the canopy [W/m2]
         eflx_sh_grnd            => energyflux_inst%eflx_sh_grnd_patch      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_lwrad_net          => energyflux_inst%eflx_lwrad_net_patch    , & ! Input:  [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
         eflx_sh_snow            => energyflux_inst%eflx_sh_snow_patch      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]
         eflx_sh_soil            => energyflux_inst%eflx_sh_soil_patch      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_h2osfc          => energyflux_inst%eflx_sh_h2osfc_patch    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from surface water (W/m**2) [+ to atm]
         eflx_bot                => energyflux_inst%eflx_bot_col            , & ! Input:  [real(r8) (:)   ]  heat flux from beneath column (W/m**2) [+ = upward]
         eflx_fgr12              => energyflux_inst%eflx_fgr12_col          , & ! Output: [real(r8) (:)   ]  heat flux between soil layer 1 and 2 (W/m2)
         eflx_fgr                => energyflux_inst%eflx_fgr_col            , & ! Output: [real(r8) (:,:) ]  (rural) soil downward heat flux (W/m2) (1:nlevgrnd)
         eflx_traffic            => energyflux_inst%eflx_traffic_lun        , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
         eflx_traffic_patch      => energyflux_inst%eflx_traffic_patch      , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
         eflx_wasteheat          => energyflux_inst%eflx_wasteheat_lun      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_wasteheat_patch    => energyflux_inst%eflx_wasteheat_patch    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_heat_from_ac       => energyflux_inst%eflx_heat_from_ac_lun   , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         eflx_heat_from_ac_patch => energyflux_inst%eflx_heat_from_ac_patch , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         eflx_anthro             => energyflux_inst%eflx_anthro_patch       , & ! Input:  [real(r8) (:)   ]  total anthropogenic heat flux (W/m**2)  
         dgnetdT                 => energyflux_inst%dgnetdT_patch           , & ! Output: [real(r8) (:)   ]  temperature derivative of ground net heat flux  
         eflx_gnet               => energyflux_inst%eflx_gnet_patch         , & ! Output: [real(r8) (:)   ]  net ground heat flux into the surface (W/m**2)
         eflx_building_heat_errsoi => energyflux_inst%eflx_building_heat_errsoi_col, & ! Output: [real(r8) (:)]  heat flux from urban building interior to walls, roof (W/m**2)
         eflx_urban_ac_col       => energyflux_inst%eflx_urban_ac_col       , & ! Output: [real(r8) (:)   ]  urban air conditioning flux (W/m**2)    
         eflx_urban_heat_col     => energyflux_inst%eflx_urban_heat_col     , & ! Output: [real(r8) (:)   ]  urban heating flux (W/m**2)             

         emg                     => temperature_inst%emg_col                , & ! Input:  [real(r8) (:)   ]  ground emissivity                       
         tssbef                  => temperature_inst%t_ssbef_col            , & ! Input:  [real(r8) (:,:) ]  temperature at previous time step [K] 
         t_h2osfc                => temperature_inst%t_h2osfc_col           , & ! Output: [real(r8) (:)   ]  surface water temperature               
         t_soisno                => temperature_inst%t_soisno_col           , & ! Output: [real(r8) (:,:) ]  soil temperature [K]             
         t_grnd                  => temperature_inst%t_grnd_col             , & ! Output: [real(r8) (:)   ]  ground surface temperature [K]          
         t_building              => temperature_inst%t_building_lun         , & ! Output: [real(r8) (:)   ]  internal building air temperature [K]       
         t_roof_inner            => temperature_inst%t_roof_inner_lun       , & ! Input:  [real(r8) (:)   ]  roof inside surface temperature [K]
         t_sunw_inner            => temperature_inst%t_sunw_inner_lun       , & ! Input:  [real(r8) (:)   ]  sunwall inside surface temperature [K]
         t_shdw_inner            => temperature_inst%t_shdw_inner_lun       , & ! Input:  [real(r8) (:)   ]  shadewall inside surface temperature [K]
         xmf                     => temperature_inst%xmf_col                , & ! Output: [real(r8) (:)   ] melting or freezing within a time step [kg/m2]
         xmf_h2osfc              => temperature_inst%xmf_h2osfc_col         , & ! Output: [real(r8) (:)   ] latent heat of phase change of surface water [col]
         fact                    => temperature_inst%fact_col               , & ! Output: [real(r8) (:)   ] used in computing tridiagonal matrix [col, lev]
         c_h2osfc                => temperature_inst%c_h2osfc_col           , & ! Output: [real(r8) (:)   ] heat capacity of surface water [col] 
         
         begc                    =>    bounds%begc                          , & ! Input:  [integer        ] beginning column index
         endc                    =>    bounds%endc                            & ! Input:  [integer        ] ending column index
         )

      ! Get step size

      dtime = get_step_size_real()

      if ( IsSimpleBuildTemp() ) call BuildingHAC( bounds, num_urbanl, &
                                           filter_urbanl, temperature_inst, &
                                           urbanparams_inst, urbantv_inst, &
                                           cool_on, heat_on )


      ! set up compact matrix for band diagonal solver, requires additional
      !     sub/super diagonals (1 each), and one additional row for t_h2osfc
      jtop = -9999
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         jtop(c) = snl(c)
         ! compute jbot
         if ((col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
              .or. col%itype(c) == icol_roof) ) then
            jbot(c) = nlevurb
         else
            jbot(c) = nlevgrnd
         endif
      end do


      !--------------------------------------------------------------
      ! Vertical coordinates adjustment for excess ice calculations
      !--------------------------------------------------------------
      if ( use_excess_ice ) then
         ! Save original soil depth to get put them back in et the end 
         dz_0(begc:endc,-nlevsno+1:nlevmaxurbgrnd) = dz(begc:endc,-nlevsno+1:nlevmaxurbgrnd)
         zi_0(begc:endc,-nlevsno+0:nlevmaxurbgrnd) = zi(begc:endc,-nlevsno+0:nlevmaxurbgrnd)
         z_0(begc:endc,-nlevsno+1:nlevmaxurbgrnd) = z(begc:endc,-nlevsno+1:nlevmaxurbgrnd)
         ! Adjust column depth for excess ice thickness 
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               dz(c,1:nlevmaxurbgrnd) = dz(c,1:nlevmaxurbgrnd) + excess_ice(c,1:nlevmaxurbgrnd) / denice  ! add extra layer thickness
               do j = 1, nlevmaxurbgrnd ! if excess ice amount dropped to zero there will be no adjustment
                  zi(c,j) = zi(c,j) + sum(excess_ice(c,1:j)) / denice
                  z(c,j) = (zi(c,j-1) + zi(c,j)) * 0.5_r8
               end do
            end if
         end do
      end if

      !------------------------------------------------------
      ! Compute ground surface and soil temperatures
      !------------------------------------------------------

      ! Thermal conductivity and Heat capacity

      tk_h2osfc(begc:endc) = nan
      call SoilThermProp(bounds, num_urbanc, filter_urbanc, num_nolakec, filter_nolakec, &
           tk(begc:endc, :), &
           cv(begc:endc, :), &
           tk_h2osfc(begc:endc), &
           urbanparams_inst, temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, soilstate_inst)

      ! Net ground heat flux into the surface and its temperature derivative
      ! Added a patches loop here to get the average of hs and dhsdT over
      ! all Patches on the column. Precalculate the terms that do not depend on PFT.

      call ComputeGroundHeatFluxAndDeriv(bounds, &
           num_nolakep, filter_nolakep, num_nolakec, filter_nolakec,          &
           hs_h2osfc( begc:endc ),                                            &
           hs_top_snow( begc:endc ),                                          &
           hs_soil( begc:endc ),                                              &
           hs_top( begc:endc ),                                               &
           dhsdT( begc:endc ),                                                &
           sabg_lyr_col( begc:endc, -nlevsno+1: ),                            &
           atm2lnd_inst, urbanparams_inst, canopystate_inst, waterdiagnosticbulk_inst, &
           waterfluxbulk_inst, solarabs_inst, energyflux_inst, temperature_inst)

      ! Determine heat diffusion through the layer interface and factor used in computing
      ! banded diagonal matrix and set up vector r and vectors a, b, c that define banded
      ! diagonal matrix and solve system

      call ComputeHeatDiffFluxAndFactor(bounds, num_nolakec, filter_nolakec, &
           dtime,                                                            &
           tk( begc:endc, -nlevsno+1: ),                                     &
           cv( begc:endc, -nlevsno+1: ),                                     &
           fn( begc:endc, -nlevsno+1: ),                                     &
           fact( begc:endc, -nlevsno+1: ),                                   &
           energyflux_inst, temperature_inst)

      ! compute thermal properties of h2osfc

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         if ( (h2osfc(c) > thin_sfclayer) .and. (frac_h2osfc(c) > thin_sfclayer) ) then
            c_h2osfc(c)  = max(thin_sfclayer, cpliq*h2osfc(c)/frac_h2osfc(c)  )
            dz_h2osfc(c) = max(thin_sfclayer, 1.0e-3*h2osfc(c)/frac_h2osfc(c) )
         else
            c_h2osfc(c)  = thin_sfclayer
            dz_h2osfc(c) = thin_sfclayer
         endif
      enddo


      ! Set up right-hand side vecor (vector r).

      call SetRHSVec(bounds, num_nolakec, filter_nolakec, &
           dtime,                                         &
           hs_h2osfc( begc:endc ),                        &
           hs_top_snow( begc:endc ),                      &
           hs_soil( begc:endc ),                          &
           hs_top( begc:endc ),                           &
           dhsdT( begc:endc ),                            &
           sabg_lyr_col (begc:endc, -nlevsno+1: ),        &
           tk( begc:endc, -nlevsno+1: ),                  &
           tk_h2osfc( begc:endc ),                        &
           fact( begc:endc, -nlevsno+1: ),                &
           fn( begc:endc, -nlevsno+1: ),                  &
           c_h2osfc( begc:endc ),                         &
           dz_h2osfc( begc:endc ),                        &
           temperature_inst,                              &
           waterdiagnosticbulk_inst,                      &
           rvector( begc:endc, -nlevsno: ))
      
      ! Set up the banded diagonal matrix

      call SetMatrix(bounds, num_nolakec, filter_nolakec, &
           dtime,                                         &
           nband,                                         &
           dhsdT( begc:endc ),                            &
           tk( begc:endc, -nlevsno+1: ),                  &
           tk_h2osfc( begc:endc ),                        &
           fact( begc:endc, -nlevsno+1: ),                &
           c_h2osfc( begc:endc ),                         &
           dz_h2osfc( begc:endc ),                        &
           waterdiagnosticbulk_inst,                      &
           bmatrix( begc:endc, 1:, -nlevsno: ))

      ! initialize initial temperature vector

      tvector(begc:endc, :) = nan
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         do j = snl(c)+1, 0
            tvector(c,j-1) = t_soisno(c,j)
         end do

         ! surface water layer has two coefficients
         tvector(c,0) = t_h2osfc(c)

         ! soil layers; top layer will have one offset and one extra coefficient
         tvector(c,1:nlevmaxurbgrnd) = t_soisno(c,1:nlevmaxurbgrnd)

      enddo

      call t_startf( 'SoilTempBandDiag')

      ! Solve the system

      call BandDiagonal(bounds, -nlevsno, nlevmaxurbgrnd, jtop(begc:endc), jbot(begc:endc), &
           num_nolakec, filter_nolakec, nband, bmatrix(begc:endc, :, :), &
           rvector(begc:endc, :), tvector(begc:endc, :))
      call t_stopf( 'SoilTempBandDiag')

      ! return temperatures to original array

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         do j = snl(c)+1, 0
            t_soisno(c,j) = tvector(c,j-1) !snow layers
         end do
         t_soisno(c,1:nlevmaxurbgrnd)   = tvector(c,1:nlevmaxurbgrnd)  !soil layers

         if (frac_h2osfc(c) == 0._r8) then
            t_h2osfc(c)=t_soisno(c,1)
         else
            t_h2osfc(c)              = tvector(c,0)           !surface water
         endif
      enddo

      ! Melting or Freezing

      do j = -nlevsno+1,nlevmaxurbgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            if ((col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                 .or. col%itype(c) == icol_roof) .and. j <= nlevurb) then
               if (j >= snl(c)+1) then
                  if (j <= nlevurb-1) then
                     fn1(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j == nlevurb) then
                     ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                     if ( IsSimpleBuildTemp() )then
                       ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                       ! building temperature. (See Oleson urban notes of 6/18/03).
                       ! Note new formulation for fn, this will be used below in net energey flux computations
                       fn1(c,j) = tk(c,j) * (t_building(l) - t_soisno(c,j))/(zi(c,j) - z(c,j))
                       fn(c,j)  = tk(c,j) * (t_building(l) - tssbef(c,j))/(zi(c,j) - z(c,j))

                     else
                        ! the bottom "soil" layer and the equations are derived assuming a prognostic inner
                        ! surface temperature.
                        if (ctype(c) == icol_sunwall) then
                          fn1(c,j) = tk(c,j) * (t_sunw_inner(l) - t_soisno(c,j))/(zi(c,j) - z(c,j))
                          fn(c,j)  = tk(c,j) * (t_sunw_inner(l) - tssbef(c,j))/(zi(c,j) - z(c,j))
                        else if (ctype(c) == icol_shadewall) then
                          fn1(c,j) = tk(c,j) * (t_shdw_inner(l) - t_soisno(c,j))/(zi(c,j) - z(c,j))
                          fn(c,j)  = tk(c,j) * (t_shdw_inner(l) - tssbef(c,j))/(zi(c,j) - z(c,j))
                        else if (ctype(c) == icol_roof) then
                          fn1(c,j) = tk(c,j) * (t_roof_inner(l) - t_soisno(c,j))/(zi(c,j) - z(c,j))
                          fn(c,j)  = tk(c,j) * (t_roof_inner(l) - tssbef(c,j))/(zi(c,j) - z(c,j))
                        end if
                     end if
                  end if
               end if
            else if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                 .and. col%itype(c) /= icol_roof) then
               if (j >= snl(c)+1) then
                  if (j <= nlevgrnd-1) then
                     fn1(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j == nlevgrnd) then
                     fn1(c,j) = 0._r8
                  end if
               end if
            end if
         end do
      end do

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)
         if (lun%urbpoi(l)) then
            if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall .or. col%itype(c) == icol_roof) then
               eflx_building_heat_errsoi(c) = cnfac*fn(c,nlevurb) + (1._r8-cnfac)*fn1(c,nlevurb)
            else
               eflx_building_heat_errsoi(c) = 0._r8
            end if
            if ( IsSimpleBuildTemp() )then
               if (cool_on(l)) then
                 eflx_urban_ac_col(c) = abs(eflx_building_heat_errsoi(c))
                 eflx_urban_heat_col(c) = 0._r8
               else if (heat_on(l)) then
                 eflx_urban_ac_col(c) = 0._r8
                 eflx_urban_heat_col(c) = abs(eflx_building_heat_errsoi(c))
               else
                 eflx_urban_ac_col(c) = 0._r8
                 eflx_urban_heat_col(c) = 0._r8
               end if
            end if
         end if
      end do

      ! compute phase change of h2osfc

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         xmf_h2osfc(c) = 0.
      end do

      call PhaseChangeH2osfc (bounds, num_nolakec, filter_nolakec, &
           dhsdT(bounds%begc:bounds%endc), &
           waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst, temperature_inst,energyflux_inst)

      call Phasechange_beta (bounds, num_nolakec, filter_nolakec, &
           dhsdT(bounds%begc:bounds%endc), &
           soilstate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst, energyflux_inst, temperature_inst)

      !--------------------------------------------------------------
      ! Vertical coordinates adjustment for excess ice calculations
      !--------------------------------------------------------------
      ! bringing back the soil depth to the original state
      if (use_excess_ice) then
         ! Adjust column depth for excess ice thickness 
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               dz(c,1:nlevmaxurbgrnd)=dz_0(c,1:nlevmaxurbgrnd)
               zi(c,1:nlevmaxurbgrnd)=zi_0(c,1:nlevmaxurbgrnd)
               z(c,1:nlevmaxurbgrnd)=z_0(c,1:nlevmaxurbgrnd)
            end if
         end do
      end if


      if ( IsProgBuildTemp() )then
         call BuildingTemperature(bounds, num_urbanl, filter_urbanl, num_nolakec, filter_nolakec, &
                                  tk(bounds%begc:bounds%endc, :), urbanparams_inst,               &
                                  temperature_inst, energyflux_inst, urbantv_inst, atm2lnd_inst)
      end if

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         ! this expression will (should) work whether there is snow or not
         if (snl(c) < 0) then
            if(frac_h2osfc(c) /= 0._r8) then
               t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                    + (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) * t_soisno(c,1) &
                    + frac_h2osfc(c) * t_h2osfc(c)
            else
               t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                    + (1.0_r8 - frac_sno_eff(c)) * t_soisno(c,1)
            end if

         else
            if(frac_h2osfc(c) /= 0._r8) then
               t_grnd(c) = (1._r8 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
            else
               t_grnd(c) = t_soisno(c,1)
            end if
         endif
      end do

      ! Initialize soil heat content

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)
         eflx_fgr12(c)= 0._r8
      end do

      ! Calculate soil heat content and soil plus snow heat content

      do j = -nlevsno+1,nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)

            if (j == 1) then ! this only needs to be done once
               eflx_fgr12(c) = -cnfac*fn(c,1) - (1._r8-cnfac)*fn1(c,1)
            end if
            if (j > 0 .and. j < nlevgrnd .and. (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)) then
               eflx_fgr(c,j) = -cnfac*fn(c,j) - (1._r8-cnfac)*fn1(c,j)
            else if (j == nlevgrnd .and. (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)) then
               eflx_fgr(c,j) = 0._r8
            end if

         end do
      end do

    end associate

  end subroutine SoilTemperature

  !-----------------------------------------------------------------------
  subroutine SoilThermProp (bounds, num_urbanc, filter_urbanc, num_nolakec, filter_nolakec, &
       tk, cv, tk_h2osfc, &
       urbanparams_inst, temperature_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, soilstate_inst)

    !
    ! !DESCRIPTION:
    ! Calculation of thermal conductivities and heat capacities of
    ! snow/soil layers
    ! (1) The volumetric heat capacity is calculated as a linear combination
    !     in terms of the volumetric fraction of the constituent phases.
    !
    ! (2) The thermal conductivity of soil is computed from the algorithm of
    !     Johansen (as reported by Farouki 1981), and of snow is from the
    !     formulation used in Sturm (1997) or Jordan (1991) p. 18 depending on
    !     namelist option.
    ! The thermal conductivities at the interfaces between two neighboring
    ! layers (j, j+1) are derived from an assumption that the flux across
    ! the interface is equal to that from the node j to the interface and the
    ! flux from the interface to the node j+1.
    !
    ! !USES:
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use clm_varpar      , only : nlevsno, nlevgrnd, nlevurb, nlevsoi, nlevmaxurbgrnd
    use clm_varcon      , only : denh2o, denice, tfrz, tkwat, tkice, tkair, cpice,  cpliq, thk_bedrock, csol_bedrock
    use landunit_varcon , only : istice, istwet
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv
    use clm_varctl      , only : iulog, snow_thermal_cond_method
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds 
    integer                , intent(in)    :: num_urbanc        ! number of urban columns in clump
    integer                , intent(in)    :: filter_urbanc(:)  ! urban column filter
    integer                , intent(in)    :: num_nolakec                      ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_nolakec(:)                ! column filter for non-lake points
    real(r8)               , intent(out)   :: cv( bounds%begc: , -nlevsno+1: ) ! heat capacity [J/(m2 K)                              ] [col, lev]
    real(r8)               , intent(out)   :: tk( bounds%begc: , -nlevsno+1: ) ! thermal conductivity at the layer interface [W/(m K) ] [col, lev]
    real(r8)               , intent(out)   :: tk_h2osfc( bounds%begc: )        ! thermal conductivity of h2osfc [W/(m K)              ] [col]
    type(urbanparams_type) , intent(in)    :: urbanparams_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(waterstatebulk_type)  , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: l,c,j                     ! indices
    integer  :: fc                        ! lake filtered column indices
    real(r8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
    real(r8) :: dke                       ! kersten number
    real(r8) :: fl                        ! volume fraction of liquid or unfrozen water to total water
    real(r8) :: satw                      ! relative total water content of soil.
    real(r8) :: zh2osfc

    character(len=*),parameter :: subname = 'SoilThermProp'
    !-----------------------------------------------------------------------

    call t_startf( 'SoilThermProp' )

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(cv)        == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk)        == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk_h2osfc) == (/bounds%endc/)),           sourcefile, __LINE__)

    associate(                                                 & 
         nbedrock     =>    col%nbedrock                     , & ! Input:  [real(r8) (:,:) ]  depth to bedrock (m)                                 
         snl          =>    col%snl			                 , & ! Input:  [integer  (:)   ]  number of snow layers                    
         dz           =>    col%dz			                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                       
         zi           =>    col%zi			                 , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m) 
         z            =>    col%z			                 , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                   
         
         nlev_improad =>    urbanparams_inst%nlev_improad    , & ! Input:  [integer  (:)   ]  number of impervious road layers         
         tk_wall      =>    urbanparams_inst%tk_wall	     , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban wall    
         tk_roof      =>    urbanparams_inst%tk_roof	     , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban roof    
         tk_improad   =>    urbanparams_inst%tk_improad	     , & ! Input:  [real(r8) (:,:) ]  thermal conductivity of urban impervious road
         cv_wall      =>    urbanparams_inst%cv_wall	     , & ! Input:  [real(r8) (:,:) ]  heat capacity of urban wall    
         cv_roof      =>    urbanparams_inst%cv_roof	     , & ! Input:  [real(r8) (:,:) ]  heat capacity of urban roof    
         cv_improad   =>    urbanparams_inst%cv_improad	     , & ! Input:  [real(r8) (:,:) ]  heat capacity of urban impervious road
         
         t_soisno     =>    temperature_inst%t_soisno_col    , & ! Input:  [real(r8) (:,:) ]  soil temperature [K]             
         
         frac_sno     =>    waterdiagnosticbulk_inst%frac_sno_eff_col , & ! Input:  [real(r8) (:)   ]  fractional snow covered area            
         h2osfc       =>    waterstatebulk_inst%h2osfc_col	     , & ! Input:  [real(r8) (:)   ]  surface (mm H2O)                        
         h2osno_no_layers => waterstatebulk_inst%h2osno_no_layers_col , & ! Input:  [real(r8) (:)   ]  snow not resolved into layers (mm H2O)
         h2osoi_liq   =>    waterstatebulk_inst%h2osoi_liq_col   , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                  
         h2osoi_ice   =>    waterstatebulk_inst%h2osoi_ice_col   , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                      
         excess_ice   => waterstatebulk_inst%excess_ice_col      , & ! Input:  [real(r8) (:,:) ]  excess ice lenses (kg/m2) (new) (1:nlevgrnd)
         bw           =>    waterdiagnosticbulk_inst%bw_col	     , & ! Output: [real(r8) (:,:) ]  partial density of water in the snow pack (ice + liquid) [kg/m3] 
         
         tkmg         =>    soilstate_inst%tkmg_col	         , & ! Input:  [real(r8) (:,:) ]  thermal conductivity, soil minerals  [W/m-K]
         tkdry        =>    soilstate_inst%tkdry_col	     , & ! Input:  [real(r8) (:,:) ]  thermal conductivity, dry soil (W/m/K)
         csol         =>    soilstate_inst%csol_col	         , & ! Input:  [real(r8) (:,:) ]  heat capacity, soil solids (J/m**3/K)
         watsat       =>    soilstate_inst%watsat_col	     , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         tksatu       =>    soilstate_inst%tksatu_col	     , & ! Input:  [real(r8) (:,:) ]  thermal conductivity, saturated soil [W/m-K]
         thk          =>    soilstate_inst%thk_col             & ! Output: [real(r8) (:,:) ]  thermal conductivity of each layer  [W/m-K] 
         )

      ! Thermal conductivity of soil from Farouki (1981)

      do j = -nlevsno+1,nlevgrnd
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)

            ! Only examine levels from 1->nlevgrnd
            if (j >= 1) then    
               l = col%landunit(c)

               ! This will include pervious road for all nlevgrnd layers and impervious road for j > nlev_improad
               if ((lun%itype(l) /= istwet .and. lun%itype(l) /= istice &
                  .and. col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall .and. &
                  col%itype(c) /= icol_roof .and. col%itype(c) /= icol_road_imperv) .or. &
                  (col%itype(c) == icol_road_imperv .and. j > nlev_improad(l))) then

                  satw = (h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice +excess_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
                  satw = min(1._r8, satw)
                  if (satw > .1e-6_r8) then
                     if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                        dke = max(0._r8, log10(satw) + 1.0_r8)
                     else                               ! Frozen soil
                        dke = satw
                     end if
                     fl = (h2osoi_liq(c,j)/(denh2o*dz(c,j))) / (h2osoi_liq(c,j)/(denh2o*dz(c,j)) + &
                          h2osoi_ice(c,j)/(denice*dz(c,j))+excess_ice(c,j)/(denice*dz(c,j)))
                     dksat = tkmg(c,j)*tkwat**(fl*watsat(c,j))*tkice**((1._r8-fl)*watsat(c,j))
                     thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
                  else
                     thk(c,j) = tkdry(c,j)
                  endif
                  if (j > nbedrock(c)) thk(c,j) = thk_bedrock
               else if (lun%itype(l) == istice) then
                  thk(c,j) = tkwat
                  if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
               else if (lun%itype(l) == istwet) then                         
                  if (j > nlevsoi) then 
                     thk(c,j) = thk_bedrock
                  else
                     thk(c,j) = tkwat
                     if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
                  endif
               endif
            endif

            ! Thermal conductivity of snow
            ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
            if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then  
               bw(c,j) = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/(frac_sno(c)*dz(c,j))
               select case (snow_thermal_cond_method)
               case ('Jordan1991')
                  thk(c,j) = tkair + (7.75e-5_r8 *bw(c,j) + 1.105e-6_r8*bw(c,j)*bw(c,j))*(tkice-tkair)
               case ('Sturm1997')
                  ! Implemented by Vicky Dutch (VRD), Nick Rutter, and
                  ! Leanne Wake (LMW)
                  ! https://tc.copernicus.org/articles/16/4201/2022/
                  ! Code provided by Adrien Dams to Will Wieder
                  if (bw(c,j) <= 156) then !LMW or 0.156 ?
                     thk(c,j) = 0.023 + 0.234*(bw(c,j)/1000) !LMW - units changed by VRD
                  else !LMW
                     thk(c,j) = 0.138 - 1.01*(bw(c,j)/1000) +(3.233*((bw(c,j)/1000)*(bw(c,j)/1000))) ! LMW Sturm I think
                  end if
               case default
                  write(iulog,*) subname//' ERROR: unknown snow_thermal_cond_method value: ', snow_thermal_cond_method
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               end select
            end if

         end do
      end do

      do j = 1,nlevurb
         do fc = 1, num_urbanc
            c = filter_urbanc(fc)
            l = col%landunit(c)

            if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
               thk(c,j) = tk_wall(l,j)
            else if (col%itype(c) == icol_roof) then
               thk(c,j) = tk_roof(l,j)
            else if (col%itype(c) == icol_road_imperv .and. j <= nlev_improad(l)) then
               thk(c,j) = tk_improad(l,j)
            end if

         end do
      end do

      ! Thermal conductivity at the layer interface

      do j = -nlevsno+1,nlevmaxurbgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if ((col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                 .or. col%itype(c) == icol_roof) .and. j <= nlevurb) then
               if (j >= snl(c)+1 .AND. j <= nlevurb-1) then
                  tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                       /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
               else if (j == nlevurb) then

                  ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                  ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                  ! building temperature. (See Oleson urban notes of 6/18/03).
                  tk(c,j) = thk(c,j)
               end if
            else if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                 .and. col%itype(c) /= icol_roof) then
               if (j >= snl(c)+1 .AND. j <= nlevgrnd-1) then
                  tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                       /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
               else if (j == nlevgrnd) then
                  tk(c,j) = 0._r8
               end if
            end if
         end do
      end do

      ! calculate thermal conductivity of h2osfc
      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         zh2osfc=1.0e-3*(0.5*h2osfc(c)) !convert to [m] from [mm]
         tk_h2osfc(c)= tkwat*thk(c,1)*(z(c,1)+zh2osfc) &
              /(tkwat*z(c,1)+thk(c,1)*zh2osfc)
      enddo

      ! Soil heat capacity, from de Vires (1963)

      do j = 1, nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            if ((lun%itype(l) /= istwet .and. lun%itype(l) /= istice &
               .and. col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall .and. &
               col%itype(c) /= icol_roof .and. col%itype(c) /= icol_road_imperv) .or. &
               (col%itype(c) == icol_road_imperv .and. j > nlev_improad(l))) then
               cv(c,j) = csol(c,j)*(1._r8-watsat(c,j))*dz(c,j) + (h2osoi_ice(c,j)*cpice + &
                         h2osoi_liq(c,j)*cpliq) + excess_ice(c,j)*cpice
               if (j > nbedrock(c)) cv(c,j) = csol_bedrock*dz(c,j)
            else if (lun%itype(l) == istwet) then 
               cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
               if (j > nbedrock(c)) cv(c,j) = csol_bedrock*dz(c,j)
            else if (lun%itype(l) == istice) then
               cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
            endif
         enddo
      end do

      do j = 1, nlevurb
         do fc = 1,num_urbanc
            c = filter_urbanc(fc)
            l = col%landunit(c)
            if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
               cv(c,j) = cv_wall(l,j) * dz(c,j)
            else if (col%itype(c) == icol_roof) then
               cv(c,j) = cv_roof(l,j) * dz(c,j)
            else if (col%itype(c) == icol_road_imperv .and. j <= nlev_improad(l)) then
               cv(c,j) = cv_improad(l,j) * dz(c,j)
            endif
          end do
      end do

      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         if (h2osno_no_layers(c) > 0._r8) then
            cv(c,1) = cv(c,1) + cpice*h2osno_no_layers(c)
         end if
      end do

      ! Snow heat capacity

      do j = -nlevsno+1,0
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
               if (frac_sno(c) > 0._r8) then
                  cv(c,j) = max(thin_sfclayer,(cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j))/frac_sno(c))
               else
                  cv(c,j) = thin_sfclayer
               endif
            end if
         end do
      end do
      call t_stopf( 'SoilThermProp' )

    end associate

  end subroutine SoilThermProp

  !-----------------------------------------------------------------------
  subroutine PhaseChangeH2osfc (bounds, num_nolakec, filter_nolakec, &
       dhsdT, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst, temperature_inst,energyflux_inst)
    !
    ! !DESCRIPTION:
    ! Only freezing is considered.  When water freezes, move ice to bottom snow layer.
    !
    ! !USES:
    use clm_time_manager , only : get_step_size_real
    use clm_varcon       , only : tfrz, hfus, grav, denice, cnfac, cpice, cpliq
    use clm_varpar       , only : nlevsno, nlevgrnd
    use clm_varctl       , only : iulog
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                         
    integer                , intent(in)    :: num_nolakec                          ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_nolakec(:)                    ! column filter for non-lake points
    real(r8)               , intent(in)    :: dhsdT ( bounds%begc: )               ! temperature derivative of "hs" [col               ]
    type(waterstatebulk_type)  , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)   , intent(inout) :: waterfluxbulk_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(energyflux_type) , intent(inout) :: energyflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,g                       !do loop index
    integer  :: fc                          !lake filtered column indices
    real(r8) :: dtime                       !land model time step (sec)
    real(r8) :: temp1                       !temporary variables [kg/m2                    ]
    real(r8) :: h2osno_total(bounds%begc:bounds%endc)  ! total snow water (mm H2O)
    real(r8) :: hm(bounds%begc:bounds%endc) !energy residual [W/m2                         ]
    real(r8) :: xm(bounds%begc:bounds%endc) !melting or freezing within a time step [kg/m2 ]
    real(r8) :: tinc                        !t(n+1)-t(n) [K]
    real(r8) :: smp                         !frozen water potential (mm)
    real(r8) :: rho_avg                     !average density
    real(r8) :: z_avg                       !average of snow depth 
    real(r8) :: c1                          !weight to use for lowest snow layer
    real(r8) :: c2                          !weight to use for surface water layer
    !-----------------------------------------------------------------------

    call t_startf( 'PhaseChangeH2osfc' )

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(dhsdT) == (/bounds%endc/)), sourcefile, __LINE__)

    associate(                                                                   & 
          eflx_h2osfc_to_snow_col  => energyflux_inst%eflx_h2osfc_to_snow_col  , & ! Output: [real(r8) (:)   ] col snow melt to h2osfc heat flux (W/m**2)
         snl                       =>    col%snl                               , & ! Input:  [integer  (:)   ] number of snow layers                    
         dz                        =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer thickness (m)                    
         
         frac_sno                  =>    waterdiagnosticbulk_inst%frac_sno_eff_col      , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
         frac_h2osfc               =>    waterdiagnosticbulk_inst%frac_h2osfc_col       , & ! Input:  [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         h2osno_no_layers          =>    waterstatebulk_inst%h2osno_no_layers_col  , & ! Output: [real(r8) (:)   ] snow that is not resolved into layers (mm H2O)
         h2osoi_ice                =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Input:  [real(r8) (:,:) ] ice lens (kg/m2) (new)                 
         h2osfc                    =>    waterstatebulk_inst%h2osfc_col            , & ! Output: [real(r8) (:)   ] surface water (mm)                      
         int_snow                  =>    waterstatebulk_inst%int_snow_col          , & ! Output: [real(r8) (:)   ] integrated snowfall [mm]               
         snow_depth                =>    waterdiagnosticbulk_inst%snow_depth_col        , & ! Output: [real(r8) (:)   ] snow height (m)                          
         
         qflx_h2osfc_to_ice        =>    waterfluxbulk_inst%qflx_h2osfc_to_ice_col , & ! Output: [real(r8) (:)   ] conversion of h2osfc to ice             
         
         fact                      =>    temperature_inst%fact_col      , &
         c_h2osfc                  =>    temperature_inst%c_h2osfc_col  , &
         xmf_h2osfc                =>    temperature_inst%xmf_h2osfc_col, &
         t_soisno                  =>    temperature_inst%t_soisno_col         , & ! Output: [real(r8) (:,:) ] soil temperature [K]              
         t_h2osfc                  =>    temperature_inst%t_h2osfc_col           & ! Output: [real(r8) (:)   ] surface water temperature               
         )

      ! Get step size

      dtime = get_step_size_real()

      ! Initialization

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         xmf_h2osfc(c) = 0._r8
         hm(c) = 0._r8
         xm(c) = 0._r8
         qflx_h2osfc_to_ice(c) = 0._r8
         eflx_h2osfc_to_snow_col(c) = 0._r8
      end do

      call waterstatebulk_inst%CalculateTotalH2osno(bounds, num_nolakec, filter_nolakec, &
           caller = 'PhaseChangeH2osfc', &
           h2osno_total = h2osno_total(bounds%begc:bounds%endc))

      ! Freezing identification
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)

         ! If liquid exists below melt point, freeze some to ice.
         if ( frac_h2osfc(c) > 0._r8 .AND. t_h2osfc(c) <= tfrz) then
            tinc = tfrz - t_h2osfc(c)
            t_h2osfc(c) = tfrz

            ! energy absorbed beyond freezing temperature
            hm(c) = frac_h2osfc(c)*(dhsdT(c)*tinc - tinc*c_h2osfc(c)/dtime)

            ! mass of water converted from liquid to ice
            xm(c) = hm(c)*dtime/hfus  
            temp1 = h2osfc(c) + xm(c)

            z_avg=frac_sno(c)*snow_depth(c)
            if (z_avg > 0._r8) then 
               rho_avg=min(800._r8,h2osno_total(c)/z_avg)
            else
               rho_avg=200._r8
            endif

            !=====================  xm < h2osfc  ====================================
            if(temp1 >= 0._r8) then ! add some frozen water to snow column

               ! add ice to snow column
               int_snow(c) = int_snow(c) - xm(c)
               if (snl(c) == 0) then
                  h2osno_no_layers(c) = h2osno_no_layers(c) - xm(c)
               else
                  h2osoi_ice(c,0) = h2osoi_ice(c,0) - xm(c)
               end if
               h2osno_total(c) = h2osno_total(c) - xm(c)

               ! remove ice from h2osfc
               h2osfc(c) = h2osfc(c) + xm(c)

               xmf_h2osfc(c) = hm(c)
               qflx_h2osfc_to_ice(c) = -xm(c)/dtime


               ! update snow depth
               if (frac_sno(c) > 0 .and. snl(c) < 0) then 
                  snow_depth(c)=h2osno_total(c)/(rho_avg*frac_sno(c))
               else
                  snow_depth(c)=h2osno_total(c)/denice
               endif

               ! adjust temperature of lowest snow layer to account for addition of ice
               if (snl(c) == 0) then
                  !initialize for next time step
                  t_soisno(c,0) = t_h2osfc(c)
                  eflx_h2osfc_to_snow_col(c) = 0.
               else
                  if (snl(c) == -1)then
                     c1=frac_sno(c)*(dtime/fact(c,0) - dhsdT(c)*dtime)
                  else
                     c1=frac_sno(c)/fact(c,0)*dtime
                  end if
                  if ( frac_h2osfc(c) /= 0.0_r8 )then
                     c2=(-cpliq*xm(c) - frac_h2osfc(c)*dhsdT(c)*dtime)
                  else
                     c2=0.0_r8
                  end if
                  t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc(c)) &
                       /(c1 + c2)             
                  eflx_h2osfc_to_snow_col(c) =(t_h2osfc(c)-t_soisno(c,0))*c2/dtime
                  
               endif

               !=========================  xm > h2osfc  =============================
            else !all h2osfc converted to ice

               rho_avg=(h2osno_total(c)*rho_avg + h2osfc(c)*denice)/(h2osno_total(c) + h2osfc(c))
               int_snow(c) = int_snow(c) + h2osfc(c)
               if (snl(c) == 0) then
                  h2osno_no_layers(c) = h2osno_no_layers(c) + h2osfc(c)
               else
                  h2osoi_ice(c,0) = h2osoi_ice(c,0) + h2osfc(c)
               end if
               h2osno_total(c) = h2osno_total(c) + h2osfc(c)

               qflx_h2osfc_to_ice(c) = h2osfc(c)/dtime

               ! excess energy is used to cool ice layer
               !
               ! NOTE: should compute and then use the heat capacity of frozen h2osfc layer
               !       rather than using heat capacity of the liquid layer. But this causes 
               !       balance check errors as it doesn't know about it.

               ! cool frozen h2osfc layer with extra heat
               t_h2osfc(c) = t_h2osfc(c) - temp1*hfus/(dtime*dhsdT(c) - c_h2osfc(c))

               xmf_h2osfc(c) = (hm(c) - frac_h2osfc(c)*temp1*hfus/dtime)

               ! next, determine equilibrium temperature of combined ice/snow layer
               if (snl(c) == 0) then
                  !initialize for next time step
                  t_soisno(c,0) = t_h2osfc(c)
               else if (snl(c) == -1) then
                  c1=frac_sno(c)*(dtime/fact(c,0) - dhsdT(c)*dtime)
                  if ( frac_h2osfc(c) /= 0.0_r8 )then
                     c2=frac_h2osfc(c)*(c_h2osfc(c) - dtime*dhsdT(c))

                  else
                     c2=0.0_r8
                  end if
                  t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc(c)) &
                       /(c1 + c2)             
                  t_h2osfc(c) = t_soisno(c,0)

               else
                  c1=frac_sno(c)/fact(c,0)*dtime
                  if ( frac_h2osfc(c) /= 0.0_r8 )then
                     c2=frac_h2osfc(c)*(c_h2osfc(c) - dtime*dhsdT(c))
                  else
                     c2=0.0_r8
                  end if
                  t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc(c)) &
                       /(c1 + c2)             
                  t_h2osfc(c) = t_soisno(c,0)
               endif

               ! set h2osfc to zero (all liquid converted to ice)
               h2osfc(c) = 0._r8

               ! update snow depth
               if (frac_sno(c) > 0 .and. snl(c) < 0) then 
                  snow_depth(c)=h2osno_total(c)/(rho_avg*frac_sno(c))
               else
                  snow_depth(c)=h2osno_total(c)/denice
               endif

            endif
         endif
      enddo
      call t_stopf( 'PhaseChangeH2osfc' )

    end associate

  end subroutine PhaseChangeH2osfc

  !-----------------------------------------------------------------------
  subroutine Phasechange_beta (bounds, num_nolakec, filter_nolakec, dhsdT, &
       soilstate_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, waterfluxbulk_inst, energyflux_inst, temperature_inst)
    !
    ! !DESCRIPTION:
    ! Calculation of the phase change within snow and soil layers:
    ! (1) Check the conditions for which the phase change may take place,
    !     i.e., the layer temperature is great than the freezing point
    !     and the ice mass is not equal to zero (i.e. melting),
    !     or the layer temperature is less than the freezing point
    !     and the liquid water mass is greater than the allowable supercooled 
    !     liquid water calculated from freezing point depression (i.e. freezing).
    ! (2) Assess the rate of phase change from the energy excess (or deficit)
    !     after setting the layer temperature to freezing point.
    ! (3) Re-adjust the ice and liquid mass, and the layer temperature
    !
    ! !USES:
    use clm_time_manager , only : get_step_size_real
    use clm_varpar       , only : nlevsno, nlevgrnd, nlevurb, nlevmaxurbgrnd
    use clm_varctl       , only : iulog
    use clm_varcon       , only : tfrz, hfus, grav, denice
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv
    use landunit_varcon  , only : istsoil, istcrop, istice
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                      
    integer                , intent(in)    :: num_nolakec                          ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_nolakec(:)                    ! column filter for non-lake points
    real(r8)               , intent(in)    :: dhsdT ( bounds%begc: )               ! temperature derivative of "hs" [col]
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(waterstatebulk_type)  , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)   , intent(inout) :: waterfluxbulk_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,g,l                            !do loop index
    integer  :: fc                                 !lake filtered column indices
    real(r8) :: dtime                              !land model time step (sec)
    real(r8) :: heatr                              !energy residual or loss after melting or freezing
    real(r8) :: temp1                              !temporary variables [kg/m2]
    real(r8) :: hm(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)    !energy residual [W/m2]
    real(r8) :: xm(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)    !melting or freezing within a time step [kg/m2]
    real(r8) :: xm2(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)   !additional melting or freezing within a time step [kg/m2] (needed for excess ice melt)
    real(r8) :: wmass0(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)!initial mass of ice and liquid (kg/m2)
    real(r8) :: wice0 (bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)!initial mass of ice (kg/m2)
    real(r8) :: wliq0 (bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)!initial mass of liquid (kg/m2)
    real(r8) :: supercool(bounds%begc:bounds%endc,nlevmaxurbgrnd)        !supercooled water in soil (kg/m2) 
    real(r8) :: propor                                                   !proportionality constant (-)
    real(r8) :: tinc(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)  !t(n+1)-t(n) [K]
    real(r8) :: smp                                                      !frozen water potential (mm)
    real(r8) :: wexice0(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd) !initial mass of excess_ice at the timestep (kg/m2)


    !-----------------------------------------------------------------------

    call t_startf( 'PhaseChangebeta' )

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(dhsdT) == (/bounds%endc/)), sourcefile, __LINE__)

    associate(                                                        & 
         snl              =>    col%snl                             , & ! Input:  [integer  (:)   ] number of snow layers                    
         dz               =>    col%dz                              , & ! Input:  [real(r8) (:,:) ] layer thickness (m)                    
         
         bsw              =>    soilstate_inst%bsw_col              , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"               
         sucsat           =>    soilstate_inst%sucsat_col           , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)              
         watsat           =>    soilstate_inst%watsat_col           , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)
         
         frac_sno_eff     =>    waterdiagnosticbulk_inst%frac_sno_eff_col    , & ! Input:  [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)
         frac_h2osfc      =>    waterdiagnosticbulk_inst%frac_h2osfc_col     , & ! Input:  [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         snow_depth       =>    waterdiagnosticbulk_inst%snow_depth_col      , & ! Input:  [real(r8) (:)   ] snow height (m)                         
         exice_subs_col   =>    waterdiagnosticbulk_inst%exice_subs_col      , & ! Output: [real(r8) (:,:) ]  per layer subsidence due to excess ice melt (mm/s)
         h2osno_no_layers =>    waterstatebulk_inst%h2osno_no_layers_col     , & ! Output: [real(r8) (:)   ] snow not resolved into layers (mm H2O)
         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col           , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)             
         h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col           , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2) (new)                 
         excess_ice       =>    waterstatebulk_inst%excess_ice_col           , & ! Input:  [real(r8) (:,:) ]  excess ice (kg/m2) (new) (1:nlevgrnd)
         
         qflx_snow_drain  =>    waterfluxbulk_inst%qflx_snow_drain_col  , & ! Output: [real(r8) (:)   ] drainage from snow pack                           
         qflx_snofrz_lyr  =>    waterfluxbulk_inst%qflx_snofrz_lyr_col  , & ! Output: [real(r8) (:,:) ] snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
         qflx_snofrz      =>    waterfluxbulk_inst%qflx_snofrz_col      , & ! Output: [real(r8) (:)   ] column-integrated snow freezing rate (positive definite) [kg m-2 s-1]
         qflx_snomelt     =>    waterfluxbulk_inst%qflx_snomelt_col     , & ! Output: [real(r8) (:)   ] snow melt (mm H2O /s)
         snomelt_accum    =>    waterdiagnosticbulk_inst%snomelt_accum_col , & ! Output: [real(r8) (:)   ] accumulated snow melt (m)
         qflx_snomelt_lyr =>    waterfluxbulk_inst%qflx_snomelt_lyr_col , & ! Output: [real(r8) (:)   ] snow melt in each layer (mm H2O /s)
         
         eflx_snomelt     =>    energyflux_inst%eflx_snomelt_col    , & ! Output: [real(r8) (:)   ] snow melt heat flux (W/m**2)
         eflx_snomelt_r   =>    energyflux_inst%eflx_snomelt_r_col  , & ! Output: [real(r8) (:)   ] rural snow melt heat flux (W/m**2)
         eflx_snomelt_u   =>    energyflux_inst%eflx_snomelt_u_col  , & ! Output: [real(r8) (:)   ] urban snow melt heat flux (W/m**2)
         
         xmf              =>    temperature_inst%xmf_col            , &
         fact             =>    temperature_inst%fact_col           , &
         
         imelt            =>    temperature_inst%imelt_col          , & ! Output: [integer  (:,:) ] flag for melting (=1), freezing (=2), Not=0 (new)
         t_soisno         =>    temperature_inst%t_soisno_col         & ! Output: [real(r8) (:,:) ] soil temperature [K]              
         )

      ! Get step size

      dtime = get_step_size_real()

      ! Initialization

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)

         xmf(c) = 0._r8
         qflx_snomelt(c)      = 0._r8
         qflx_snofrz(c)       = 0._r8
         qflx_snow_drain(c)   = 0._r8
      end do

      do j = -nlevsno+1,nlevmaxurbgrnd       ! all layers
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if (j >= snl(c)+1) then

               ! Initialization
               imelt(c,j) = 0
               hm(c,j) = 0._r8
               xm(c,j) = 0._r8
               xm2(c,j) = 0._r8
               wice0(c,j) = h2osoi_ice(c,j)
               wliq0(c,j) = h2osoi_liq(c,j)
               wexice0(c,j) = excess_ice(c,j)
               wmass0(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j) + wexice0(c,j)
               if (j >= 1) then
                  exice_subs_col(c,j) = 0._r8
               endif
            endif   ! end of snow layer if-block

            if (j <= 0) then
               ! Do for all possible snow layers in case snl changes over timestep.
               qflx_snomelt_lyr(c,j) = 0._r8
               qflx_snofrz_lyr(c,j)  = 0._r8
            end if
         end do   ! end of column-loop
      enddo   ! end of level-loop

      !--  snow layers  --------------------------------------------------- 
      do j = -nlevsno+1,0             
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if (j >= snl(c)+1) then

               ! Melting identification
               ! If ice exists above melt point, melt some to liquid.
               if (h2osoi_ice(c,j) > 0._r8 .and. t_soisno(c,j) > tfrz) then
                  imelt(c,j) = 1
                  !                tinc(c,j) = t_soisno(c,j) - tfrz 
                  tinc(c,j) = tfrz - t_soisno(c,j) 
                  t_soisno(c,j) = tfrz
               endif

               ! Freezing identification
               ! If liquid exists below melt point, freeze some to ice.
               if (h2osoi_liq(c,j) > 0._r8 .AND. t_soisno(c,j) < tfrz) then
                  imelt(c,j) = 2
                  !                tinc(c,j) = t_soisno(c,j) - tfrz 
                  tinc(c,j) = tfrz - t_soisno(c,j) 
                  t_soisno(c,j) = tfrz
               endif
            endif   ! end of snow layer if-block
         end do   ! end of column-loop
      enddo   ! end of level-loop

      !-- soil layers   --------------------------------------------------- 
      do j = 1,nlevmaxurbgrnd             
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            supercool(c,j) = 0.0_r8
            ! add in urban condition if-block
            if ((col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                 .and. col%itype(c) /= icol_roof .and. j <= nlevgrnd) .or. &
                 ((col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                 .or. col%itype(c) == icol_roof) .and. j <= nlevurb)) then

               if (h2osoi_ice(c,j) > 0. .AND. t_soisno(c,j) > tfrz) then
                  imelt(c,j) = 1
                  !             tinc(c,j) = t_soisno(c,j) - tfrz 
                  tinc(c,j) = tfrz - t_soisno(c,j)
                  t_soisno(c,j) = tfrz
               endif

               ! melt excess ice after normal ice
               if (excess_ice(c,j) > 0._r8 .AND. t_soisno(c,j) > tfrz) then
                  imelt(c,j) = 1
                  tinc(c,j) = tfrz - t_soisno(c,j)
                  t_soisno(c,j) = tfrz
               endif

               ! from Zhao (1997) and Koren (1999)
               supercool(c,j) = 0.0_r8
               if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop .or. col%itype(c) == icol_road_perv) then
                  if(t_soisno(c,j) < tfrz) then
                     smp = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
                     supercool(c,j) = watsat(c,j)*(smp/sucsat(c,j))**(-1._r8/bsw(c,j))
                     supercool(c,j) = supercool(c,j)*dz(c,j)*1000._r8       ! (mm)
                  endif
               endif

               if (h2osoi_liq(c,j) > supercool(c,j) .AND. t_soisno(c,j) < tfrz) then
                  imelt(c,j) = 2
                  !             tinc(c,j) = t_soisno(c,j) - tfrz
                  tinc(c,j) = tfrz - t_soisno(c,j)
                  t_soisno(c,j) = tfrz
               endif

               ! If snow exists, but its thickness is less than the critical value (0.01 m)
               if (h2osno_no_layers(c) > 0._r8 .AND. j == 1) then
                  if (t_soisno(c,j) > tfrz) then
                     imelt(c,j) = 1
                     !                tincc,j) = t_soisno(c,j) - tfrz
                     tinc(c,j) = tfrz - t_soisno(c,j)
                     t_soisno(c,j) = tfrz
                  endif
               endif

            endif

         end do
      enddo


      do j = -nlevsno+1,nlevmaxurbgrnd       ! all layers
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)

            if ((col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                 .and. col%itype(c) /= icol_roof .and. j <= nlevgrnd) .or. &
                 ((col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                 .or. col%itype(c) == icol_roof) .and. j <= nlevurb)) then

               if (j >= snl(c)+1) then

                  ! Calculate the energy surplus and loss for melting and freezing
                  if (imelt(c,j) > 0) then

                     ! added unique cases for this calculation,
                     ! to account for absorbed solar radiation in each layer

                     !==================================================================
                     if (j == snl(c)+1) then ! top layer                   
                        if(j > 0) then
                           hm(c,j) = dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j)
                        else
                           hm(c,j) = frac_sno_eff(c)*(dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j))
                        endif

                        if ( j==1 .and. frac_h2osfc(c) /= 0.0_r8 ) then
                           hm(c,j) = hm(c,j) - frac_h2osfc(c)*(dhsdT(c)*tinc(c,j))
                        end if
                     else if (j == 1) then
                        hm(c,j) = (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) &
                             *dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j)
                     else ! non-interfacial snow/soil layers                   
                        if(j < 1) then
                           hm(c,j) = - frac_sno_eff(c)*(tinc(c,j)/fact(c,j))
                        else
                           hm(c,j) = - tinc(c,j)/fact(c,j)
                        endif
                     endif
                  endif

                  ! These two errors were checked carefully (Y. Dai).  They result from the
                  ! computed error of "Tridiagonal-Matrix" in subroutine "thermal".
                  if (imelt(c,j) == 1 .AND. hm(c,j) < 0._r8) then
                     hm(c,j) = 0._r8
                     imelt(c,j) = 0
                  endif
                  if (imelt(c,j) == 2 .AND. hm(c,j) > 0._r8) then
                     hm(c,j) = 0._r8
                     imelt(c,j) = 0
                  endif

                  ! The rate of melting and freezing

                  if (imelt(c,j) > 0 .and. abs(hm(c,j)) > 0._r8) then
                     xm(c,j) = hm(c,j)*dtime/hfus                           ! kg/m2

                     ! If snow exists, but its thickness is less than the critical value
                     ! (1 cm). Note: more work is needed to determine how to tune the
                     ! snow depth for this case
                     if (j == 1) then
                        if (h2osno_no_layers(c) > 0._r8 .AND. xm(c,j) > 0._r8) then
                           temp1 = h2osno_no_layers(c)                     ! kg/m2
                           h2osno_no_layers(c) = max(0._r8,temp1-xm(c,j))
                           propor = h2osno_no_layers(c)/temp1
                           snow_depth(c) = propor * snow_depth(c)
                           heatr = hm(c,j) - hfus*(temp1-h2osno_no_layers(c))/dtime   ! W/m2
                           if (heatr > 0._r8) then
                              xm(c,j) = heatr*dtime/hfus                    ! kg/m2
                              hm(c,j) = heatr                               ! W/m2
                           else
                              xm(c,j) = 0._r8
                              hm(c,j) = 0._r8
                           endif
                           qflx_snomelt(c) = max(0._r8,(temp1-h2osno_no_layers(c)))/dtime   ! kg/(m2 s)
                           ! no snow layers, so qflx_snomelt_lyr is not set
                           xmf(c) = hfus*qflx_snomelt(c)
                           qflx_snow_drain(c) = qflx_snomelt(c)
                        endif
                     endif

                     heatr = 0._r8
                     if (xm(c,j) > 0._r8) then !if there is excess heat to melt the ice
                        h2osoi_ice(c,j) = max(0._r8, wice0(c,j)-xm(c,j))
                        heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                        xm2(c,j) = xm(c,j) - h2osoi_ice(c,j) !excess ice melting
                        if (h2osoi_ice(c,j) == 0._r8) then ! this might be redundant 
                           if (excess_ice(c,j) >= 0._r8 .and. xm2(c,j)>0._r8 .and. j>=2) then ! if there is excess ice to melt
                              excess_ice(c,j) = max(0._r8,wexice0(c,j) - xm2(c,j))
                              heatr = hm(c,j) - hfus * (wexice0(c,j)-excess_ice(c,j)+wice0(c,j)-h2osoi_ice(c,j)) / dtime
                           endif
                        endif !end of excess ice block
                     else if (xm(c,j) < 0._r8) then
                        if (j <= 0) then
                           h2osoi_ice(c,j) = min(wmass0(c,j), wice0(c,j)-xm(c,j))  ! snow
                        else
                           if (wmass0(c,j) - wexice0(c,j) < supercool(c,j)) then  ! even if excess ice is present, it cannot refreeze
                              h2osoi_ice(c,j) = 0._r8
                           else
                              h2osoi_ice(c,j) = min(wmass0(c,j) - wexice0(c,j) - supercool(c,j),wice0(c,j)-xm(c,j))
                           endif
                        endif
                        heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                     endif

                     h2osoi_liq(c,j) = max(0._r8,wmass0(c,j)-h2osoi_ice(c,j)-excess_ice(c,j)) !melted excess ice is added to the respective soil layers
                     

                     if (abs(heatr) > 0._r8) then
                        if (j == snl(c)+1) then

                           if(j==1) then
                              t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                                   /(1._r8-(1.0_r8 - frac_h2osfc(c))*fact(c,j)*dhsdT(c))
                           else
                              t_soisno(c,j) = t_soisno(c,j) + (fact(c,j)/frac_sno_eff(c))*heatr &
                                   /(1._r8-fact(c,j)*dhsdT(c))

                           endif

                        else if (j == 1) then

                           t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                                /(1._r8-(1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c))*fact(c,j)*dhsdT(c))
                        else
                           if(j > 0) then
                              t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr
                           else
                              if(frac_sno_eff(c) > 0._r8) t_soisno(c,j) = t_soisno(c,j) + (fact(c,j)/frac_sno_eff(c))*heatr
                           endif
                        endif

                        if (j <= 0) then    ! snow
                           if (h2osoi_liq(c,j)*h2osoi_ice(c,j)>0._r8) t_soisno(c,j) = tfrz
                        end if
                     endif  ! end of heatr > 0 if-block

                     if (j >= 1) then
                        xmf(c) = xmf(c) + hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime + &
                                 hfus*(wexice0(c,j)-excess_ice(c,j))/dtime
                        ! subsidence calculation
                        exice_subs_col(c,j) = max(0._r8, (wexice0(c,j)-excess_ice(c,j))/denice) 
                     else
                        xmf(c) = xmf(c) + hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                     endif
                     
                     if (imelt(c,j) == 1 .AND. j < 1) then
                        qflx_snomelt_lyr(c,j) = max(0._r8,(wice0(c,j)-h2osoi_ice(c,j)))/dtime
                        qflx_snomelt(c)       = qflx_snomelt(c) + qflx_snomelt_lyr(c,j)
                        snomelt_accum(c)      = snomelt_accum(c) + qflx_snomelt_lyr(c,j) * dtime * 1.e-3_r8
                     endif

                     ! layer freezing mass flux (positive):
                     if (imelt(c,j) == 2 .AND. j < 1) then
                        qflx_snofrz_lyr(c,j) = max(0._r8,(h2osoi_ice(c,j)-wice0(c,j)))/dtime
                        qflx_snofrz(c)       = qflx_snofrz(c) + qflx_snofrz_lyr(c,j)
                     endif

                  endif

               endif   ! end of snow layer if-block

            endif

         end do   ! end of column-loop
      enddo   ! end of level-loop


      ! Needed for history file output

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         eflx_snomelt(c) = qflx_snomelt(c) * hfus
         l = col%landunit(c)
         if (lun%urbpoi(l)) then
            eflx_snomelt_u(c) = eflx_snomelt(c)
         else if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
            eflx_snomelt_r(c) = eflx_snomelt(c)
         end if
      end do

      call t_stopf( 'PhaseChangebeta' )
    end associate

  end subroutine Phasechange_beta

  !-----------------------------------------------------------------------
  subroutine ComputeGroundHeatFluxAndDeriv(bounds, &
       num_nolakep, filter_nolakep, num_nolakec, filter_nolakec, &
       hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, &
       atm2lnd_inst, urbanparams_inst, canopystate_inst, waterdiagnosticbulk_inst, &
       waterfluxbulk_inst, solarabs_inst, energyflux_inst, temperature_inst)
    !
    ! !DESCRIPTION:
    ! Computes ground heat flux on:
    ! (1) The surface of standing water,
    ! (2) The surface of snow,
    ! (3) The surface of soil, and
    ! (4) Net energy flux into soil surface.
    ! Additionally, derivative of ground heat flux w.r.t to temeprature
    !
    ! !USES:
    use clm_varcon     , only : sb, hvap
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno
    use UrbanParamsType, only : IsSimpleBuildTemp, IsProgBuildTemp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds                                    ! bounds
    integer                , intent(in)    :: num_nolakep                               ! number of non-lake points in patch filter
    integer                , intent(in)    :: filter_nolakep( : )                       ! patch filter for non-lake points
    integer                , intent(in)    :: num_nolakec                               ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_nolakec( : )                       ! column filter for non-lake points
    real(r8)               , intent(out)   :: hs_h2osfc( bounds%begc: )                 ! heat flux on standing water [W/m2]
    real(r8)               , intent(out)   :: hs_top_snow( bounds%begc: )               ! heat flux on top snow layer [W/m2]
    real(r8)               , intent(out)   :: hs_soil( bounds%begc: )                   ! heat flux on soil [W/m2]
    real(r8)               , intent(out)   :: hs_top (bounds%begc: )                    ! net energy flux into surface layer (col) [W/m2]
    real(r8)               , intent(out)   :: dhsdT( bounds%begc: )                     ! temperature derivative of "hs" [col]
    real(r8)               , intent(out)   :: sabg_lyr_col( bounds%begc:, -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(urbanparams_type) , intent(in)    :: urbanparams_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(waterdiagnosticbulk_type)  , intent(in)    :: waterdiagnosticbulk_inst
    type(waterfluxbulk_type)   , intent(in)    :: waterfluxbulk_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,p,l,g                                              ! indices
    integer  :: fc, fp                                                 ! lake filtered column and patch indices
    real(r8) :: hs(bounds%begc:bounds%endc)                            ! net energy flux into the surface (w/m2)
    real(r8) :: lwrad_emit(bounds%begc:bounds%endc)                    ! emitted longwave radiation
    real(r8) :: dlwrad_emit(bounds%begc:bounds%endc)                   ! time derivative of emitted longwave radiation
    integer  :: lyr_top                                                ! index of top layer of snowpack (-4 to 0) [idx]
    real(r8) :: eflx_gnet_top                                          ! net energy flux into surface layer, patch-level [W/m2]
    real(r8) :: lwrad_emit_snow(bounds%begc:bounds%endc)               !
    real(r8) :: lwrad_emit_soil(bounds%begc:bounds%endc)               !
    real(r8) :: lwrad_emit_h2osfc(bounds%begc:bounds%endc)             !
    real(r8) :: eflx_gnet_snow                                         !
    real(r8) :: eflx_gnet_soil                                         !
    real(r8) :: eflx_gnet_h2osfc                                       !
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(hs_h2osfc)     == (/bounds%endc/)),   sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hs_top_snow)   == (/bounds%endc/)),   sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hs_soil)       == (/bounds%endc/)),   sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hs_top)        == (/bounds%endc/)),   sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dhsdT)         == (/bounds%endc/)),   sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sabg_lyr_col)  == (/bounds%endc,1/)), sourcefile, __LINE__)

    associate(                                                                &
         snl                     => col%snl                                 , & ! Input:  [integer (:)    ]  number of snow layers
         z                       => col%z                                   , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         
         forc_lwrad              => atm2lnd_inst%forc_lwrad_downscaled_col  , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)
         
         frac_veg_nosno          => canopystate_inst%frac_veg_nosno_patch   , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         
         frac_sno_eff            => waterdiagnosticbulk_inst%frac_sno_eff_col        , & ! Input:  [real(r8) (:)   ]  eff. fraction of ground covered by snow (0 to 1)
         
         qflx_ev_snow            => waterfluxbulk_inst%qflx_ev_snow_patch       , & ! Input:  [real(r8) (:)   ]  evaporation flux from snow (mm H2O/s) [+ to atm]
         qflx_ev_soil            => waterfluxbulk_inst%qflx_ev_soil_patch       , & ! Input:  [real(r8) (:)   ]  evaporation flux from soil (mm H2O/s) [+ to atm]
         qflx_ev_h2osfc          => waterfluxbulk_inst%qflx_ev_h2osfc_patch     , & ! Input:  [real(r8) (:)   ]  evaporation flux from h2osfc (mm H2O/s) [+ to atm]
         qflx_evap_soi           => waterfluxbulk_inst%qflx_evap_soi_patch      , & ! Input:  [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_tran_veg           => waterfluxbulk_inst%qflx_tran_veg_patch      , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         
         emg                     => temperature_inst%emg_col                , & ! Input:  [real(r8) (:)   ]  ground emissivity                       
         t_h2osfc                => temperature_inst%t_h2osfc_col           , & ! Input:  [real(r8) (:)   ]  surface water temperature               
         t_grnd                  => temperature_inst%t_grnd_col             , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]          
         t_soisno                => temperature_inst%t_soisno_col           , & ! Input:  [real(r8) (:,:) ]  soil temperature [K]             
         
         htvp                    => energyflux_inst%htvp_col                , & ! Input:  [real(r8) (:)   ]  latent heat of vapor of water (or sublimation) [j/kg]
         cgrnd                   => energyflux_inst%cgrnd_patch             , & ! Input:  [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]
         dlrad                   => energyflux_inst%dlrad_patch             , & ! Input:  [real(r8) (:)   ]  downward longwave radiation blow the canopy [W/m2]
         eflx_traffic            => energyflux_inst%eflx_traffic_lun        , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
         eflx_wasteheat          => energyflux_inst%eflx_wasteheat_lun      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_ventilation        => energyflux_inst%eflx_ventilation_lun    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from building ventilation (W/m**2)
         eflx_heat_from_ac       => energyflux_inst%eflx_heat_from_ac_lun   , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         eflx_sh_snow            => energyflux_inst%eflx_sh_snow_patch      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]
         eflx_sh_soil            => energyflux_inst%eflx_sh_soil_patch      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_h2osfc          => energyflux_inst%eflx_sh_h2osfc_patch    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from surface water (W/m**2) [+ to atm]
         eflx_sh_grnd            => energyflux_inst%eflx_sh_grnd_patch      , & ! Input:  [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_lwrad_net          => energyflux_inst%eflx_lwrad_net_patch    , & ! Input:  [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
         eflx_wasteheat_patch    => energyflux_inst%eflx_wasteheat_patch    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_ventilation_patch  => energyflux_inst%eflx_ventilation_patch  , & ! Input:  [real(r8) (:)   ]  sensible heat flux from building ventilation (W/m**2)
         eflx_heat_from_ac_patch => energyflux_inst%eflx_heat_from_ac_patch , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         eflx_traffic_patch      => energyflux_inst%eflx_traffic_patch      , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
         eflx_anthro             => energyflux_inst%eflx_anthro_patch       , & ! Input:  [real(r8) (:)   ]  total anthropogenic heat flux (W/m**2)  
         eflx_gnet               => energyflux_inst%eflx_gnet_patch         , & ! Output: [real(r8) (:)   ]  net ground heat flux into the surface (W/m**2)
         dgnetdT                 => energyflux_inst%dgnetdT_patch           , & ! Output: [real(r8) (:)   ]  temperature derivative of ground net heat flux  
         
         sabg                    => solarabs_inst%sabg_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
         sabg_soil               => solarabs_inst%sabg_soil_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
         sabg_snow               => solarabs_inst%sabg_snow_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
         sabg_chk                => solarabs_inst%sabg_chk_patch            , & ! Output: [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
         sabg_lyr                => solarabs_inst%sabg_lyr_patch            , & ! Output: [real(r8) (:,:) ]  absorbed solar radiation (pft,lyr) [W/m2]
         
         begc                    => bounds%begc                             , & ! Input:  [integer        ] beginning column index
         endc                    => bounds%endc                               & ! Input:  [integer        ] ending column index
         )

      ! Net ground heat flux into the surface and its temperature derivative
      ! Added a pfts loop here to get the average of hs and dhsdT over
      ! all PFTs on the column. Precalculate the terms that do not depend on PFT.

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         lwrad_emit(c)  =    emg(c) * sb * t_grnd(c)**4
         dlwrad_emit(c) = 4._r8*emg(c) * sb * t_grnd(c)**3

         ! fractionate lwrad_emit; balanced in CanopyFluxes & Biogeophysics2
         lwrad_emit_snow(c)    =    emg(c) * sb * t_soisno(c,snl(c)+1)**4
         lwrad_emit_soil(c)    =    emg(c) * sb * t_soisno(c,1)**4
         lwrad_emit_h2osfc(c)  =    emg(c) * sb * t_h2osfc(c)**4
      end do

      hs_soil(begc:endc)   = 0._r8
      hs_h2osfc(begc:endc) = 0._r8
      hs(begc:endc)        = 0._r8
      dhsdT(begc:endc)     = 0._r8
      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)

         if (.not. lun%urbpoi(l)) then
            eflx_gnet(p) = sabg(p) + dlrad(p) &
                 + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit(c) &
                 - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
            ! save sabg for balancecheck, in case frac_sno is set to zero later
            sabg_chk(p) = frac_sno_eff(c) * sabg_snow(p) + (1._r8 - frac_sno_eff(c) ) * sabg_soil(p)

            eflx_gnet_snow = sabg_snow(p) + dlrad(p) &
                 + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_snow(c) &
                 - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

            eflx_gnet_soil = sabg_soil(p) + dlrad(p) &
                 + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_soil(c) &
                 - (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

            eflx_gnet_h2osfc = sabg_soil(p) + dlrad(p) &
                 + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_h2osfc(c) &
                 - (eflx_sh_h2osfc(p)+qflx_ev_h2osfc(p)*htvp(c))
         else
            ! For urban columns we use the net longwave radiation (eflx_lwrad_net) because of
            ! interactions between urban columns.

            ! All wasteheat and traffic flux goes into canyon floor
            if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then
               ! Note that we divide the following landunit variables by 1-wtlunit_roof which 
               ! essentially converts the flux from W/m2 of urban area to W/m2 of canyon floor area
               eflx_wasteheat_patch(p) = eflx_wasteheat(l)/(1._r8-lun%wtlunit_roof(l))
               if ( IsSimpleBuildTemp() ) then
                  eflx_ventilation_patch(p) = 0._r8
               else if ( IsProgBuildTemp() ) then
                  eflx_ventilation_patch(p) = eflx_ventilation(l)/(1._r8-lun%wtlunit_roof(l))
               end if
               eflx_heat_from_ac_patch(p) = eflx_heat_from_ac(l)/(1._r8-lun%wtlunit_roof(l))
               eflx_traffic_patch(p) = eflx_traffic(l)/(1._r8-lun%wtlunit_roof(l))
            else
               eflx_wasteheat_patch(p) = 0._r8
               eflx_ventilation_patch(p) = 0._r8
               eflx_heat_from_ac_patch(p) = 0._r8
               eflx_traffic_patch(p) = 0._r8
            end if
            ! Include transpiration term because needed for previous road
            ! and include wasteheat and traffic flux
            eflx_gnet(p) = sabg(p) + dlrad(p)  &
                 - eflx_lwrad_net(p) &
                 - (eflx_sh_grnd(p) + qflx_evap_soi(p)*htvp(c) + qflx_tran_veg(p)*hvap) &
                 + eflx_wasteheat_patch(p) + eflx_heat_from_ac_patch(p) + eflx_traffic_patch(p) &
                 + eflx_ventilation_patch(p)
            if ( IsSimpleBuildTemp() ) then
               eflx_anthro(p)   = eflx_wasteheat_patch(p) + eflx_traffic_patch(p)
            end if
            eflx_gnet_snow   = eflx_gnet(p)
            eflx_gnet_soil   = eflx_gnet(p)
            eflx_gnet_h2osfc = eflx_gnet(p)
         end if
         dgnetdT(p) = - cgrnd(p) - dlwrad_emit(c)
         hs(c) = hs(c) + eflx_gnet(p) * patch%wtcol(p)
         dhsdT(c) = dhsdT(c) + dgnetdT(p) * patch%wtcol(p)
         ! separate surface fluxes for soil/snow
         hs_soil(c) = hs_soil(c) + eflx_gnet_soil * patch%wtcol(p)
         hs_h2osfc(c) = hs_h2osfc(c) + eflx_gnet_h2osfc * patch%wtcol(p)
      end do

      ! Additional calculations with SNICAR:
      ! Set up tridiagonal matrix in a new manner. There is now
      ! absorbed solar radiation in each snow layer, instead of
      ! only the surface. Following the current implementation,
      ! absorbed solar flux should be: S + ((delS/delT)*dT),
      ! where S is absorbed radiation, and T is temperature. Now,
      ! assume delS/delT is zero, then it is OK to just add S
      ! to each layer

      ! Initialize:
      sabg_lyr_col(begc:endc,-nlevsno+1:1) = 0._r8
      hs_top(begc:endc)                    = 0._r8
      hs_top_snow(begc:endc)               = 0._r8

      do fp = 1,num_nolakep
         p = filter_nolakep(fp)
         c = patch%column(p)
         l = patch%landunit(p)

         lyr_top = snl(c) + 1

         if (.not. lun%urbpoi(l)) then

            eflx_gnet_top = sabg_lyr(p,lyr_top) + dlrad(p) + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
                 - lwrad_emit(c) - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))

            hs_top(c) = hs_top(c) + eflx_gnet_top*patch%wtcol(p)

            eflx_gnet_snow = sabg_lyr(p,lyr_top) + dlrad(p) + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
                 - lwrad_emit_snow(c) - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

            eflx_gnet_soil = sabg_lyr(p,lyr_top) + dlrad(p) + (1._r8-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
                 - lwrad_emit_soil(c) - (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

            hs_top_snow(c) = hs_top_snow(c) + eflx_gnet_snow*patch%wtcol(p)

            do j = lyr_top,1,1
               sabg_lyr_col(c,j) = sabg_lyr_col(c,j) + sabg_lyr(p,j) * patch%wtcol(p)
            enddo
         else

            hs_top(c)      = hs_top(c) + eflx_gnet(p)*patch%wtcol(p)
            hs_top_snow(c) = hs_top_snow(c) + eflx_gnet(p)*patch%wtcol(p)
            sabg_lyr_col(c,lyr_top) = sabg_lyr_col(c,lyr_top) + sabg(p) * patch%wtcol(p)

         endif
      enddo

    end associate

  end subroutine ComputeGroundHeatFluxAndDeriv

  !-----------------------------------------------------------------------
  subroutine ComputeHeatDiffFluxAndFactor(bounds, num_nolakec, filter_nolakec, dtime, &
       tk, cv, fn, fact, &
       energyflux_inst, temperature_inst)
    !
    ! !DESCRIPTION:
    ! Computes:
    ! (1) Heat diffusion at the interface of layers.
    ! (2) Factor used in computing tridiagonal matrix
    !
    ! !USES:
    use clm_varcon     , only : capr, cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb, nlevmaxurbgrnd
    use UrbanParamsType, only : IsSimpleBuildTemp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)  :: bounds                             ! bounds
    integer                , intent(in)  :: num_nolakec                        ! number of column non-lake points in column filter
    integer                , intent(in)  :: filter_nolakec(:)                  ! column filter for non-lake points
    real(r8)               , intent(in)  :: dtime                              ! land model time step (sec)
    real(r8)               , intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )     ! thermal conductivity [W/(m K)]
    real(r8)               , intent(in)  :: cv (bounds%begc: ,-nlevsno+1: )    ! heat capacity [J/(m2 K)]
    real(r8)               , intent(out) :: fn (bounds%begc: ,-nlevsno+1: )    ! heat diffusion through the layer interface [W/m2]
    real(r8)               , intent(out) :: fact( bounds%begc: , -nlevsno+1: ) ! used in computing tridiagonal matrix [col, lev]
    type(energyflux_type)  , intent(in)  :: energyflux_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                           ! indices
    integer  :: fc                                              ! lake filtered column indices
    real(r8) :: dzm                                             ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(tk)   == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(cv)   == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fact) == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fn)   == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)

    associate(&
         zi         => col%zi                          , & ! Input: [real(r8) (:,:) ] interface level below a "z" level (m)
         dz         => col%dz                          , & ! Input: [real(r8) (:,:) ] layer depth (m)
         z          => col%z                           , & ! Input: [real(r8) (:,:) ] layer thickness (m)
         ctype      => col%itype                       , & ! Input: [integer (:)    ]  column type
         t_building => temperature_inst%t_building_lun , & ! Input: [real(r8) (:)   ] internal building temperature [K]       
         t_roof_inner => temperature_inst%t_roof_inner_lun , & ! Input: [real(r8) (:)   ] roof inside surface temperature [K]
         t_sunw_inner => temperature_inst%t_sunw_inner_lun , & ! Input: [real(r8) (:)   ] sunwall inside surface temperature [K]
         t_shdw_inner => temperature_inst%t_shdw_inner_lun , & ! Input: [real(r8) (:)   ] shadewall inside surface temperature [K]
         t_soisno   => temperature_inst%t_soisno_col   , & ! Input: [real(r8) (:,:) ] soil temperature [K]             
         eflx_bot   => energyflux_inst%eflx_bot_col      & ! Input: [real(r8) (:)   ] heat flux from beneath column (W/m**2) [+ = upward]
         )

      ! Determine heat diffusion through the layer interface and factor used in computing
      ! tridiagonal matrix and set up vector r and vectors a, b, c that define tridiagonal
      ! matrix and solve system

      do j = -nlevsno+1,nlevmaxurbgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            if ((col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                 .or. col%itype(c) == icol_roof) .and. j <= nlevurb) then
               if (j >= col%snl(c)+1) then
                  if (j == col%snl(c)+1) then
                     fact(c,j) = dtime/cv(c,j)
                     fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j <= nlevurb-1) then
                     fact(c,j) = dtime/cv(c,j)
                     fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                     dzm     = (z(c,j)-z(c,j-1))
                  else if (j == nlevurb) then
                     fact(c,j) = dtime/cv(c,j)
                     if ( IsSimpleBuildTemp() )then
                       ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                       ! building temperature. (See Oleson urban notes of 6/18/03).
                       fn(c,j) = tk(c,j) * (t_building(l) - cnfac*t_soisno(c,j))/(zi(c,j) - z(c,j))
                     else
                        ! the bottom "soil" layer and the equations are derived assuming a prognostic inner
                        ! surface temperature.
                        if (ctype(c) == icol_sunwall) then
                           fn(c,j) = tk(c,j) * (t_sunw_inner(l) - cnfac*t_soisno(c,j))/(zi(c,j) - z(c,j))
                        else if (ctype(c) == icol_shadewall) then
                           fn(c,j) = tk(c,j) * (t_shdw_inner(l) - cnfac*t_soisno(c,j))/(zi(c,j) - z(c,j))
                        else if (ctype(c) == icol_roof) then
                           fn(c,j) = tk(c,j) * (t_roof_inner(l) - cnfac*t_soisno(c,j))/(zi(c,j) - z(c,j))
                        end if
                     end if
                  end if
               end if
            else if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                 .and. col%itype(c) /= icol_roof .and. j <= nlevgrnd) then
               if (j >= col%snl(c)+1) then
                  if (j == col%snl(c)+1) then
                     fact(c,j) = dtime/cv(c,j) * dz(c,j) / (0.5_r8*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))
                     fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                  else if (j <= nlevgrnd-1) then
                     fact(c,j) = dtime/cv(c,j)
                     fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                     dzm     = (z(c,j)-z(c,j-1))
                  else if (j == nlevgrnd) then
                     fact(c,j) = dtime/cv(c,j)
                     fn(c,j) = eflx_bot(c)
                  end if
               end if
            end if
         end do
      end do

    end associate

  end subroutine ComputeHeatDiffFluxAndFactor

  !-----------------------------------------------------------------------
  subroutine SetRHSVec(bounds, num_nolakec, filter_nolakec, dtime, &
       hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, tk, &
       tk_h2osfc, fact, fn, c_h2osfc, dz_h2osfc, &
       temperature_inst, waterdiagnosticbulk_inst, rvector)

    !
    ! !DESCRIPTION:
    ! Setup the RHS-Vector for the numerical solution of temperature for snow,
    ! standing surface water and soil layers.
    !
    !           |===========|
    !           |   Snow    |
    !           !===========|
    ! rvector = |    SSW    |
    !           !===========|
    !           !   Soil    |
    !           !===========|
    !
    ! !USES:
    use clm_varcon      , only : cnfac, cpliq
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar      , only : nlevsno, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)  :: bounds                            ! bounds
    integer  , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer  , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8) , intent(in)  :: dtime                                      ! land model time step (sec)
    real(r8) , intent(in)  :: hs_h2osfc( bounds%begc: )                  ! heat flux on standing water [W/m2]
    real(r8) , intent(in)  :: hs_top_snow( bounds%begc: )                ! heat flux on top snow layer [W/m2]
    real(r8) , intent(in)  :: hs_soil( bounds%begc: )                    ! heat flux on soil [W/m2]
    real(r8) , intent(in)  :: hs_top( bounds%begc: )                     ! net energy flux into surface layer (col) [W/m2]
    real(r8) , intent(in)  :: dhsdT( bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8) , intent(in)  :: sabg_lyr_col( bounds%begc: , -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) , intent(in)  :: tk( bounds%begc: , -nlevsno+1: )           ! thermal conductivity [W/(m K)]
    real(r8) , intent(in)  :: tk_h2osfc( bounds%begc: )                  ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8) , intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8) , intent(in)  :: fn( bounds%begc: , -nlevsno+1: )           ! heat diffusion through the layer interface [W/m2]
    real(r8) , intent(in)  :: c_h2osfc( bounds%begc: )                   ! heat capacity of surface water [col]
    real(r8) , intent(in)  :: dz_h2osfc( bounds%begc: )                  ! Thickness of standing water [m]
    real(r8) , intent(out) :: rvector( bounds%begc: , -nlevsno: )        ! RHS vector used in numerical solution of temperature
    type(temperature_type) , intent(in) :: temperature_inst
    type(waterdiagnosticbulk_type) , intent(in) :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                                     ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: rt (bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)  ! "r" vector for tridiagonal solution
    real(r8) :: fn_h2osfc(bounds%begc:bounds%endc)                      ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8) :: rt_snow(bounds%begc:bounds%endc,-nlevsno:-1)            ! RHS vector corresponding to snow layers
    real(r8) :: rt_ssw(bounds%begc:bounds%endc,1)                       ! RHS vector corresponding to standing surface water
    real(r8) :: rt_soil(bounds%begc:bounds%endc,1:nlevmaxurbgrnd)       ! RHS vector corresponding to soil layer
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(hs_h2osfc)    == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hs_top_snow)  == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hs_soil)      == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hs_top)       == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dhsdT)        == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk)           == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk_h2osfc)    == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fact)         == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fn)           == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(c_h2osfc)     == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dz_h2osfc)    == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rvector)      == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)

    associate(                                                       &
         t_soisno     => temperature_inst%t_soisno_col             , & ! Input: [real(r8) (:,:) ] soil temperature [K]
         t_h2osfc     => temperature_inst%t_h2osfc_col             , & ! Input: [real(r8) (:)   ] surface water temperature
         frac_h2osfc  => waterdiagnosticbulk_inst%frac_h2osfc_col  , & ! Input: [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         frac_sno_eff => waterdiagnosticbulk_inst%frac_sno_eff_col , & ! Input: [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)
         begc         => bounds%begc                               , & ! Input: [integer        ] beginning column index
         endc         => bounds%endc                                 & ! Input: [integer        ] ending column index
         )

      ! Initialize
      rvector(begc:endc, :) = nan

      !SetRHSVec_ subroutines must be called in the correct order:
      ! 1) SetRHSVec_Snow 
      ! 2) SetRHSVec_StandingSurfaceWater
      ! 3) SetRHSVec_Soil
      !

      ! Set entries in RHS vector for snow layers
      call SetRHSVec_Snow(bounds, num_nolakec, filter_nolakec, &
           hs_top_snow( begc:endc ),                           &
           hs_top( begc:endc ),                                &
           dhsdT( begc:endc ),                                 &
           sabg_lyr_col (begc:endc, -nlevsno+1: ),             &
           fact( begc:endc, -nlevsno+1: ),                     &
           fn( begc:endc, -nlevsno+1: ),                       &
           t_soisno ( begc:endc, -nlevsno+1: ),                &
           t_h2osfc ( begc:endc ),                             &
           rt_snow( begc:endc, -nlevsno:))

      ! Set entries in RHS vector for surface water layer
      call SetRHSVec_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, &
           dtime,                                                              &
           hs_h2osfc( begc:endc ),                                             &
           dhsdT( begc:endc ),                                                 &
           tk_h2osfc( begc:endc ),                                             &
           c_h2osfc( begc:endc ),                                              &
           dz_h2osfc( begc:endc ),                                             &
           fn_h2osfc( begc:endc ),                                             &
           t_soisno ( begc:endc, -nlevsno+1: ),                                &
           t_h2osfc ( begc:endc),                                              &
           rt_ssw( begc:endc, 1:1))

      ! Set entries in RHS vector for soil layers
      call SetRHSVec_Soil(bounds, num_nolakec, filter_nolakec, &
           hs_top_snow( begc:endc ),                           &
           hs_soil( begc:endc ),                               &
           hs_top( begc:endc ),                                &
           dhsdT( begc:endc ),                                 &
           sabg_lyr_col (begc:endc, -nlevsno+1: ),             &
           fact( begc:endc, -nlevsno+1: ),                     &
           fn( begc:endc, -nlevsno+1: ),                       &
           fn_h2osfc( begc:endc ),                             &
           c_h2osfc( begc:endc ),                              &
           frac_h2osfc ( begc:endc),                           &
           frac_sno_eff( begc:endc),                           &
           t_soisno ( begc:endc, -nlevsno+1: ),                &
           rt_soil( begc:endc, 1: ))

      ! Combine the RHS vector
      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         rvector(c, -nlevsno:-1) = rt_snow(c, -nlevsno:-1)
         rvector(c, 0          ) = rt_ssw(c, 1           )
         rvector(c, 1:nlevmaxurbgrnd ) = rt_soil(c, 1:nlevmaxurbgrnd )
      end do

    end associate

  end subroutine SetRHSVec

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_Snow(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_top, dhsdT, sabg_lyr_col, &
       fact, fn, t_soisno, t_h2osfc, rt)
    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to snow layers for all columns.
    !
    ! !USES:
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                             ! bounds
    integer , intent(in)  :: num_nolakec                                ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                          ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow( bounds%begc: )                ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_top( bounds%begc: )                     ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT( bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col( bounds%begc: , -nlevsno+1: ) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )         ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: , -nlevsno+1: )           ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: t_soisno(bounds%begc:, -nlevsno+1:)        ! soil temperature [K] 
    real(r8), intent(in)  :: t_h2osfc(bounds%begc:)                     ! surface water temperature [K] 
    real(r8), intent(out) :: rt(bounds%begc: , -nlevsno: )              ! rhs vector entries
    !-----------------------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                   ! indices
    integer  :: fc                                                      ! lake filtered column indices
    real(r8) :: dzp, dzm                                                ! used in computing tridiagonal matrix
    real(r8) :: hs_top_lev(bounds%endc)

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(hs_top_snow)  == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(hs_top)       == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dhsdT)        == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fact)         == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fn)           == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_soisno)     == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_h2osfc)     == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rt)           == (/bounds%endc, -1/)),       sourcefile, __LINE__)

    associate(                    &
         begc =>    bounds%begc , & ! Input:  [integer        ] beginning column index
         endc =>    bounds%endc , & ! Input:  [integer        ] ending column index
         z    =>    col%z         & ! Input:  [real(r8) (:,:) ] layer thickness [m]
         )

      ! Initialize
      rt(begc:endc, : ) = nan

      do j = -nlevsno+1,0
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            ! urban road and non-urban columns
            hs_top_lev(c) = hs_top_snow(c) 

            ! urban non-road columns
            if (col%itype(c) == icol_sunwall   .or. &
                col%itype(c) == icol_shadewall .or. &
                col%itype(c) == icol_roof) then

                hs_top_lev(c) = hs_top(c)

            end if
         end do
      end do


      do j = -nlevsno+1,0
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if (j == col%snl(c)+1) then
               dzp     = z(c,j+1)-z(c,j)
               rt(c,j-1) = t_soisno(c,j) +  fact(c,j)*( hs_top_lev(c) &
                      - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )

            else if (j > col%snl(c)+1) then
               dzm     = (z(c,j)-z(c,j-1))
               dzp     = (z(c,j+1)-z(c,j))

               rt(c,j-1) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
               rt(c,j-1) = rt(c,j-1) + fact(c,j)*sabg_lyr_col(c,j)
            end if
         end do
      end do

    end associate

  end subroutine SetRHSVec_Snow

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, dtime, &
       hs_h2osfc, dhsdT, tk_h2osfc, c_h2osfc, dz_h2osfc, fn_h2osfc, &
       t_soisno, t_h2osfc, rt)
    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to all standing surface water
    !
    ! !USES:
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                      ! bounds
    integer , intent(in)  :: num_nolakec                         ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                   ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                               ! land model time step (sec)
    real(r8), intent(in)  :: hs_h2osfc(bounds%begc: )            ! heat flux on standing water [W/m2]
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )            ! thermal conductivity of h2osfc [W/(m K)] [col]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )            ! heat capacity of surface water [col]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )            ! Thickness of standing water [m]
    real(r8), intent(out) :: fn_h2osfc (bounds%begc: )           ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(in)  :: t_soisno(bounds%begc:, -nlevsno+1:) ! soil temperature [K]
    real(r8), intent(in)  :: t_h2osfc(bounds%begc:)              ! surface water temperature temperature [K]
    real(r8), intent(out) :: rt(bounds%begc: , 1: )              ! rhs vector entries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                              ! indices
    integer  :: fc                                               ! lake filtered column indices
    real(r8) :: dzm                                              ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(hs_h2osfc)    == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dhsdT)        == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk_h2osfc)    == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(c_h2osfc)     == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dz_h2osfc)    == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fn_h2osfc)    == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_soisno)     == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_h2osfc)     == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rt)           == (/bounds%endc,1/)),         sourcefile, __LINE__)

    ! Initialize
    rt(bounds%begc:bounds%endc, : ) = nan

    !
    ! surface water ------------------------------------------------------------------
    !
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! surface water layer has two coefficients
       dzm=(0.5*dz_h2osfc(c)+col%z(c,1))

       fn_h2osfc(c)=tk_h2osfc(c)*(t_soisno(c,1)-t_h2osfc(c))/dzm
       rt(c,1)= t_h2osfc(c) +  (dtime/c_h2osfc(c)) &
            *( hs_h2osfc(c) - dhsdT(c)*t_h2osfc(c) + cnfac*fn_h2osfc(c) )!rhs for h2osfc

    enddo

  end subroutine SetRHSVec_StandingSurfaceWater

  !-----------------------------------------------------------------------
  subroutine SetRHSVec_Soil(bounds, num_nolakec, filter_nolakec, &
       hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col, fact, fn, fn_h2osfc, c_h2osfc, &
       frac_h2osfc, frac_sno_eff, t_soisno, rt)
    !
    ! !DESCRIPTION:
    ! Sets up RHS vector corresponding to all soil layers
    !
    ! !USES:
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevurb, nlevgrnd, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                     ! bounds
    integer , intent(in)  :: num_nolakec                                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                  ! column filter for non-lake points
    real(r8), intent(in)  :: hs_top_snow(bounds%begc: )                         ! heat flux on top snow layer [W/m2]
    real(r8), intent(in)  :: hs_soil(bounds%begc: )                             ! heat flux on soil [W/m2]
    real(r8), intent(in)  :: hs_top(bounds%begc: )                              ! net energy flux into surface layer (col) [W/m2]
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                               ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: sabg_lyr_col(bounds%begc:, -nlevsno+1: )           ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                 ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: fn (bounds%begc: ,-nlevsno+1: )                    ! heat diffusion through the layer interface [W/m2]
    real(r8), intent(in)  :: fn_h2osfc (bounds%begc: )                          ! heat diffusion through standing-water/soil interface [W/m2]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                           ! heat capacity of surface water [col]
    real(r8), intent(in)  :: frac_h2osfc(bounds%begc: )                         ! fractional area with surface water greater than zero
    real(r8), intent(in)  :: frac_sno_eff(bounds%begc: )                        ! fraction of ground covered by snow (0 to 1)
    real(r8), intent(in)  :: t_soisno(bounds%begc:, -nlevsno+1:)                ! soil temperature [K] 
    real(r8), intent(out) :: rt(bounds%begc: ,1: )                              ! rhs vector entries
    !-----------------------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                           ! indices
    integer  :: fc                                              ! lake filtered column indices
    !-----------------------------------------------------------------------
    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(hs_soil)      == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dhsdT)        == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(sabg_lyr_col) == (/bounds%endc, 1/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fact)         == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fn)           == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fn_h2osfc)    == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(c_h2osfc)     == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(frac_h2osfc)  == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(frac_sno_eff) == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_soisno)     == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(rt)           == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)

    associate(&
         begc     => bounds%begc  , & ! Input:  [integer ] beginning column index
         endc     => bounds%endc    & ! Input:  [integer ] ending column index
         )

      ! Initialize
      rt(begc:endc, : ) = nan

      ! urban non-road columns ------------------------------------------------------------------
      !
      do j = 1,nlevurb
        do fc = 1,num_nolakec
           c = filter_nolakec(fc)
           if (col%itype(c) == icol_sunwall  .or. &
               col%itype(c) == icol_shadewall .or. &
               col%itype(c) == icol_roof) then

               if (j == col%snl(c)+1) then
                  ! changed hs to hs_top
                  rt(c,j) = t_soisno(c,j) +  fact(c,j)*( hs_top(c) - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
               else if (j <= nlevurb-1) then
                  ! if this is a snow layer or the top soil layer,
                  ! add absorbed solar flux to factor 'rt'
                  if (j == 1) then
                     rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                     rt(c,j) = rt(c,j) + (fact(c,j)*sabg_lyr_col(c,j))
                  else
                     rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                  endif

               else if (j == nlevurb) then
                  ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                  ! the bottom "soil" layer and the equations are derived assuming a prescribed internal
                  ! building temperature. (See Oleson urban notes of 6/18/03).
                  rt(c,j) = t_soisno(c,j) + fact(c,j)*( fn(c,j) - cnfac*fn(c,j-1) )
               end if
           end if
        end do
      end do


      ! non-urban and urban road column
      do j = 1,nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            if ((.not. lun%urbpoi(l)) .or. & 
               (col%itype(c) == icol_road_imperv .or. &
                col%itype(c) == icol_road_perv)) then

               if (j == col%snl(c)+1) then
                  rt(c,j) = t_soisno(c,j) +  fact(c,j)*( hs_top_snow(c) &
                       - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
               else if (j == 1) then
                  ! this is the snow/soil interface layer
                  rt(c,j) = t_soisno(c,j) + fact(c,j) &
                       *((1._r8-frac_sno_eff(c))*(hs_soil(c) - dhsdT(c)*t_soisno(c,j)) &
                       + cnfac*(fn(c,j) - frac_sno_eff(c) * fn(c,j-1)))

                  rt(c,j) = rt(c,j) +  frac_sno_eff(c)*fact(c,j)*sabg_lyr_col(c,j)

               else if (j <= nlevgrnd-1) then
                  rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )

               else if (j == nlevgrnd) then
                  rt(c,j) = t_soisno(c,j) - cnfac*fact(c,j)*fn(c,j-1) + fact(c,j)*fn(c,j)
               end if

            end if
         end do
      end do

     !
     ! surface water  -----------------------------------------------------------------
     !
     do fc = 1,num_nolakec
        c = filter_nolakec(fc)
        if ( frac_h2osfc(c) /= 0.0_r8 )then
            rt(c,1)=rt(c,1) &
                 -frac_h2osfc(c)*fact(c,1)*((hs_soil(c) - dhsdT(c)*t_soisno(c,1)) &
                 +cnfac*fn_h2osfc(c))
         end if
      end do

    end associate

  end subroutine SetRHSVec_Soil

  !-----------------------------------------------------------------------
  subroutine SetMatrix(bounds, num_nolakec, filter_nolakec, dtime, nband, &
       dhsdT, tk, tk_h2osfc, fact, c_h2osfc, dz_h2osfc, waterdiagnosticbulk_inst, bmatrix)
    !
    ! !DESCRIPTION:
    ! Setup the matrix for the numerical solution of temperature for snow,
    ! standing surface water and soil layers.
    !
    !
    !           |===========|===========|===========|
    !           |   Snow    |           | Snow-Soil |
    !           !===========|===========|===========|
    ! bmatrix = |           |    SSW    |  SSW-Soil |
    !           !===========|===========|===========|
    !           ! Soil-Snow | Soil-SSW  |    Soil   |
    !           !===========|===========|===========|
    !
    !
    !
    ! !USES:
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                                    ! bounds
    integer , intent(in)  :: num_nolakec                                       ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                                 ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                                             ! land model time step (sec)
    integer , intent(in)  :: nband                                             ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                              ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )                    ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                          ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )                ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                          ! heat capacity of surface water [col]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                          ! Thickness of standing water [m]
    real(r8), intent(out) :: bmatrix(bounds%begc: , 1:,-nlevsno: )             ! matrix for numerical solution of temperature
    type(waterdiagnosticbulk_type), intent(in) :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                                            ! indices
    integer  :: fc                                                             ! lake filtered column indices
    real(r8) :: dzm                                                            ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                            ! used in computing tridiagonal matrix
    real(r8) :: bmatrix_snow(bounds%begc:bounds%endc,nband,-nlevsno:-1      )  ! block-diagonal matrix for snow layers
    real(r8) :: bmatrix_ssw(bounds%begc:bounds%endc,nband,       0:0       )   ! block-diagonal matrix for standing surface water
    real(r8) :: bmatrix_soil(bounds%begc:bounds%endc,nband,       1:nlevmaxurbgrnd)  ! block-diagonal matrix for soil layers
    real(r8) :: bmatrix_snow_soil(bounds%begc:bounds%endc,nband,-1:-1)         ! off-diagonal matrix for snow-soil interaction
    real(r8) :: bmatrix_ssw_soil(bounds%begc:bounds%endc,nband, 0:0 )          ! off-diagonal matrix for standing surface water-soil interaction
    real(r8) :: bmatrix_soil_snow(bounds%begc:bounds%endc,nband, 1:1 )         ! off-diagonal matrix for soil-snow interaction
    real(r8) :: bmatrix_soil_ssw(bounds%begc:bounds%endc,nband, 1:1 )          ! off-diagonal matrix for soil-standing surface water interaction
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(dhsdT)     == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk)        == (/bounds%endc, nlevmaxurbgrnd/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk_h2osfc) == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fact)      == (/bounds%endc, nlevmaxurbgrnd/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(c_h2osfc)  == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dz_h2osfc) == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix)   == (/bounds%endc, nband, nlevmaxurbgrnd/)), sourcefile, __LINE__)

    associate(                                                       &
         z            => col%z                                     , & ! Input: [real(r8) (:,:) ]  layer thickness [m]
         frac_h2osfc  => waterdiagnosticbulk_inst%frac_h2osfc_col  , & ! Input: [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
         frac_sno_eff => waterdiagnosticbulk_inst%frac_sno_eff_col , & ! Input: [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         begc         => bounds%begc                               , & ! Input: [integer        ]  beginning column index
         endc         => bounds%endc                                 & ! Input: [integer        ]  ending column index
         )

      ! Assemble smaller matrices

      call SetMatrix_Snow(bounds, num_nolakec, filter_nolakec, nband, &
           dhsdT( begc:endc ),                                        &
           tk( begc:endc, -nlevsno+1: ),                              &
           fact( begc:endc, -nlevsno+1: ),                            &
           frac_sno_eff(begc:endc),                                   &
           bmatrix_snow( begc:endc, 1:, -nlevsno: ),                  &
           bmatrix_snow_soil( begc:endc, 1:, -1: ))

      call SetMatrix_Soil(bounds, num_nolakec, filter_nolakec, nband, &
           dhsdT( begc:endc ),                                        &
           tk( begc:endc, -nlevsno+1: ),                              &
           tk_h2osfc( begc:endc ),                                    &
           dz_h2osfc( begc:endc ),                                    &
           fact( begc:endc, -nlevsno+1: ),                            &
           frac_h2osfc(begc:endc),                                    &
           frac_sno_eff(begc:endc),                                   &
           bmatrix_soil( begc:endc, 1:, 1: ),                         &
           bmatrix_soil_snow( begc:endc, 1:, 1: ))

      call SetMatrix_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, dtime, nband, &
           dhsdT( begc:endc ),                                                               &
           tk( begc:endc, -nlevsno+1: ),                                                     &
           tk_h2osfc( begc:endc ),                                                           &
           fact( begc:endc, -nlevsno+1: ),                                                   &
           c_h2osfc( begc:endc ),                                                            &
           dz_h2osfc( begc:endc ),                                                           &
           frac_h2osfc(begc:endc),                                                           &
           bmatrix_ssw( begc:endc, 1:, 0: ),                                                 &
           bmatrix_ssw_soil( begc:endc, 1:, 0: ),                                            &
           bmatrix_soil_ssw( begc:endc, 1:, 1: ))

      call AssembleMatrixFromSubmatrices(bounds, num_nolakec, filter_nolakec, nband, &
           bmatrix_snow( begc:endc, 1:, -nlevsno: ),                                 &
           bmatrix_ssw( begc:endc, 1:, 0: ),                                         &
           bmatrix_soil( begc:endc, 1:, 1: ),                                        &
           bmatrix_snow_soil( begc:endc, 1:, -1: ),                                  &
           bmatrix_ssw_soil( begc:endc, 1:, 0: ),                                    &
           bmatrix_soil_snow( begc:endc, 1:, 1: ),                                   &
           bmatrix_soil_ssw( begc:endc, 1:, 1: ),                                    &
           bmatrix( begc:endc, 1:, -nlevsno: ))

    end associate

  end subroutine SetMatrix

  !-----------------------------------------------------------------------
  subroutine AssembleMatrixFromSubmatrices(bounds, num_nolakec, filter_nolakec, nband, &
       bmatrix_snow, bmatrix_ssw, bmatrix_soil, bmatrix_snow_soil, &
       bmatrix_ssw_soil, bmatrix_soil_snow, bmatrix_soil_ssw, bmatrix)

    !
    ! !DESCRIPTION:
    ! Assemble the full matrix from submatrices.
    !
    ! Non-zero pattern of bmatrix (assuming 5 snow layers):
    !
    !        SNOW-LAYERS
    !            |
    !            |  STANDING-SURFACE-WATER
    !            |         |
    !            |         |              SOIL-LAYERS
    !            |         |                  |
    !            v         v                  v
    !
    !      -5 -4 -3 -2 -1| 0| 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    !      ==============================================================
    !  -5 | x  x         |  |                                            |
    !  -4 | x  x  x      |  |                                            |
    !  -3 |    x  x  x   |  |                                            |
    !  -2 |       x  x  x|  |                                            |
    !  -1 |          x  x|  | x                                          |
    !      ==============================================================
    !   0 |              | x| x                                          |
    !      ==============================================================
    !   1 |             x| x| x  x                                       |
    !   2 |              |  | x  x  x                                    |
    !   3 |              |  |    x  x  x                                 |
    !   4 |              |  |       x  x  x                              |
    !   5 |              |  |          x  x  x                           |
    !   6 |              |  |             x  x  x                        |
    !   7 |              |  |                x  x  x                     |
    !   8 |              |  |                   x  x  x                  |
    !   9 |              |  |                      x  x  x               |
    !  10 |              |  |                         x  x  x            |
    !  11 |              |  |                            x  x  x         |
    !  12 |              |  |                               x  x  x      |
    !  13 |              |  |                                  x  x  x   |
    !  14 |              |  |                                     x  x  x|
    !  15 |              |  |                                        x  x|
    !      ==============================================================
    !
    !
    ! !USES:
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: bmatrix_snow(bounds%begc: , 1: , -nlevsno: ) ! block-diagonal matrix for snow layers [col, nband, nlevsno]
    real(r8), intent(in)  :: bmatrix_ssw(bounds%begc: , 1: , 0: )         ! block-diagonal matrix for standing surface water [col, nband, 0:0]
    real(r8), intent(in)  :: bmatrix_soil(bounds%begc: , 1: , 1: )        ! block-diagonal matrix for soil layers [col, nband, nlevgrnd]
    real(r8), intent(in)  :: bmatrix_snow_soil(bounds%begc: , 1: , -1: )  ! off-diagonal matrix for snow-soil interaction [col, nband, -1:-1]
    real(r8), intent(in)  :: bmatrix_ssw_soil(bounds%begc: , 1: , 0: )    ! off-diagonal matrix for standing surface water-soil interaction [col, nband, 0:0]
    real(r8), intent(in)  :: bmatrix_soil_snow(bounds%begc: , 1: , 1: )   ! off-diagonal matrix for soil-snow interaction [col, nband, 1:1]
    real(r8), intent(in)  :: bmatrix_soil_ssw(bounds%begc: , 1: , 1: )    ! off-diagonal matrix for soil-standing surface water interaction [col, nband, 1:1]
    real(r8), intent(out) :: bmatrix(bounds%begc: , 1: , -nlevsno: )      ! full matrix used in numerical solution of temperature [col, nband, -nlevsno:nlevgrnd]
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c                                                                         ! indices
    integer  :: fc                                                                          ! lake filtered column indices
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(bmatrix_snow)        == (/bounds%endc, nband, -1/)),       sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_ssw)         == (/bounds%endc, nband, 0/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_soil)        == (/bounds%endc, nband, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_snow_soil)   == (/bounds%endc, nband, -1/)),       sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_ssw_soil)    == (/bounds%endc, nband, 0/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_soil_snow)   == (/bounds%endc, nband, 1/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_soil_ssw)    == (/bounds%endc, nband, 1/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix)             == (/bounds%endc, nband, nlevmaxurbgrnd/)), sourcefile, __LINE__)

    ! Assemble the full matrix

    bmatrix(bounds%begc:bounds%endc, :, :) = 0.0_r8
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! Snow
       bmatrix(c,2:3,-nlevsno   )   = bmatrix_snow(c,2:3,-nlevsno   )
       bmatrix(c,2:4,-nlevsno+1:-2) = bmatrix_snow(c,2:4,-nlevsno+1:-2)
       bmatrix(c,3:4,-1   )         = bmatrix_snow(c,3:4,-1   )

       ! Snow-Soil
       bmatrix(c,1,-1) = bmatrix_snow_soil(c,1,-1)

       ! StandingSurfaceWater
       bmatrix(c,3,0) = bmatrix_ssw(c,3,0)

       ! StandingSurfaceWater-Soil
       bmatrix(c,2,0) = bmatrix_ssw_soil(c,2,0)

       ! Soil
       bmatrix(c,2:3,1           )  = bmatrix_soil(c,2:3,1           )
       bmatrix(c,2:4,2:nlevmaxurbgrnd-1)  = bmatrix_soil(c,2:4,2:nlevmaxurbgrnd-1)
       bmatrix(c,3:4,nlevmaxurbgrnd    )  = bmatrix_soil(c,3:4,nlevmaxurbgrnd    )

       ! Soil-Snow
       bmatrix(c,5,1)  = bmatrix_soil_snow(c,5,1)

       ! Soil-StandingSurfaceWater
       bmatrix(c,4,1)  = bmatrix_soil_ssw(c,4,1)

    end do

  end subroutine AssembleMatrixFromSubmatrices

  !-----------------------------------------------------------------------
  subroutine SetMatrix_Snow(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, fact, frac_sno_eff, bmatrix_snow, bmatrix_snow_soil)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries corresponding to internal snow layers
    !
    ! !USES:
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                               ! bounds
    integer , intent(in)  :: num_nolakec                                  ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                            ! column filter for non-lake points
    integer , intent(in)  :: nband                                        ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                         ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )               ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )           ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: frac_sno_eff(bounds%begc: )                  ! fraction of ground covered by snow (0 to 1)
    real(r8), intent(out) :: bmatrix_snow(bounds%begc: , 1:, -nlevsno: )  ! matrix enteries
    real(r8), intent(out) :: bmatrix_snow_soil(bounds%begc: , 1:,-1: )    ! matrix enteries
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                       ! indices
    integer  :: fc                                                          ! lake filtered column indices
    integer  :: nlev_thresh(1:num_nolakec)
    real(r8) :: dzm                                                         ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                         ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(dhsdT)             == (/bounds%endc/)),            sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk)                == (/bounds%endc, nlevmaxurbgrnd/)),  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fact)              == (/bounds%endc, nlevmaxurbgrnd/)),  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(frac_sno_eff)      == (/bounds%endc/)),            sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_snow)      == (/bounds%endc, nband, -1/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_snow_soil) == (/bounds%endc, nband, -1/)), sourcefile, __LINE__)

    associate(                                        &
         begc =>    bounds%begc                     , & ! Input:  [integer       ] beginning column index
         endc =>    bounds%endc                     , & ! Input:  [integer       ] ending column index
         z    =>    col%z                             & ! Input:  [real(r8) (:,:)] layer thickness [m]
         )

      ! Initialize
      bmatrix_snow      (begc:endc, :, :) = 0.0_r8
      bmatrix_snow_soil (begc:endc, :, :) = 0.0_r8

      do j = -nlevsno+1,0
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            ! urban road and non-urban columns
            nlev_thresh(fc) = nlevgrnd

            ! urban non-road columns
            if (col%itype(c) == icol_sunwall  .or. & 
                col%itype(c) == icol_shadewall.or. &
                col%itype(c) == icol_roof) then

               nlev_thresh(fc) = nlevurb

            end if
         enddo
      end do


      do j = -nlevsno+1,0
        do fc = 1,num_nolakec
           c = filter_nolakec(fc)
           if (j >= col%snl(c)+1) then
              dzp     = z(c,j+1)-z(c,j)    
              if (j == col%snl(c)+1) then
                 bmatrix_snow (c,4,j-1) = 0._r8                                                        
                 bmatrix_snow (c,3,j-1) = 1._r8+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
              else
                 dzm     = (z(c,j)-z(c,j-1))
                 bmatrix_snow(c,4,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                 bmatrix_snow(c,3,j-1) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
              end if
              if (j /= 0) then
                 bmatrix_snow(c,2,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp !ct
              else !if ( j == 0)
                 bmatrix_snow_soil(c,1,j-1) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
              end if
           end if
        enddo
      end do

end associate

end subroutine SetMatrix_Snow

!-----------------------------------------------------------------------
  subroutine SetMatrix_Soil(bounds, num_nolakec, filter_nolakec, nband, &
       dhsdT, tk, tk_h2osfc, dz_h2osfc, fact, frac_h2osfc, frac_sno_eff,  bmatrix_soil, bmatrix_soil_snow)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries corresponding to internal soil layers
    ! and soil-snow interaction layer. 
    !
    ! !USES:
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon  , only : icol_road_perv, icol_road_imperv
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevurb, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                            ! bounds
    integer , intent(in)  :: num_nolakec                               ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                         ! column filter for non-lake points
    integer , intent(in)  :: nband                                     ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                      ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )            ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                  ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                  ! Thickness of standing water [m]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )        ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: frac_h2osfc(bounds%begc: )                ! fractional area with surface water greater than zero
    real(r8), intent(in)  :: frac_sno_eff(bounds%begc: )               ! fraction of ground covered by snow (0 to 1)
    real(r8), intent(out) :: bmatrix_soil(bounds%begc: , 1:, 1: )      ! matrix enteries corresponding to internal soil layers
    real(r8), intent(out) :: bmatrix_soil_snow(bounds%begc: , 1: ,1: ) ! matrix enteries corresponding to soil-snow interaction
    !
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,l                                                  ! indices
    integer  :: fc                                                     ! lake filtered column indices
    real(r8) :: dzm                                                    ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                    ! used in computing tridiagonal matrix
    ! -----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(dhsdT)             == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk)                == (/bounds%endc, nlevmaxurbgrnd/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk_h2osfc)         == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dz_h2osfc)         == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fact)              == (/bounds%endc, nlevmaxurbgrnd/)),        sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(frac_h2osfc)       == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(frac_sno_eff)      == (/bounds%endc/)),                  sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_soil)      == (/bounds%endc, nband, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_soil_snow) == (/bounds%endc, nband, 1/)),        sourcefile, __LINE__)

    associate(                      &
         begc    => bounds%begc   , & ! Input:  [integer        ] beginning column index
         endc    => bounds%endc   , & ! Input:  [integer        ] ending column index
         zi      =>    col%zi     , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level [m]
         z       =>    col%z        & ! Input:  [real(r8) (:,:) ] layer thickness [m]
         )
      ! Initialize
      bmatrix_soil      (begc:endc, :, :) = 0.0_r8
      bmatrix_soil_snow (begc:endc, :, :) = 0.0_r8


      !
      ! urban non-road columns ---------------------------------------------------------
      !
      do j = 1,nlevurb
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            if (col%itype(c) == icol_sunwall   .or. &
                col%itype(c) == icol_shadewall .or. &
                col%itype(c) == icol_roof) then

               if (j == col%snl(c)+1) then
                  dzp     = z(c,j+1)-z(c,j)
                  if (j /= 1) then
                     bmatrix_soil(c,4,j) = 0._r8
                   else
                     bmatrix_soil_snow(c,5,j) = 0._r8
                  end if
                  bmatrix_soil(c,3,j) = 1._r8+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                  bmatrix_soil(c,2,j) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
               else if (j <= nlevurb-1) then
                  dzm     = (z(c,j)-z(c,j-1))
                  dzp     = (z(c,j+1)-z(c,j))
                  if (j /= 1) then
                     bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                   else
                     bmatrix_soil_snow(c,5,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                  end if
                  bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                  bmatrix_soil(c,2,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
               else if (j == nlevurb) then
                  ! For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
                  ! the bottom "soil" layer and the equations are derived assuming a prognostic inner
                  ! surface temperature.
                  dzm     = ( z(c,j)-z(c,j-1))
                  dzp     = (zi(c,j)-z(c,j))
                  bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)*(tk(c,j-1)/dzm)
                  bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j-1)/dzm + tk(c,j)/dzp)
                  bmatrix_soil(c,2,j) = 0._r8
               end if
            end if
         enddo
      enddo

      !
      ! non-urban columns and urban road columns ---------------------------------------
      !
      do j = 1,nlevgrnd
         do fc = 1,num_nolakec
            c = filter_nolakec(fc)
            l = col%landunit(c)
            if ((col%itype(c) == icol_road_imperv) .or. &
                (col%itype(c) == icol_road_perv)   .or. &
                (.not. lun%urbpoi(l))) then

                if (j == col%snl(c)+1) then
                   dzp     = z(c,j+1)-z(c,j)
                   if (j /= 1) then
                      bmatrix_soil(c,4,j) = 0._r8
                    else 
                      bmatrix_soil_snow(c,5,j) = 0._r8
                   end if
                   bmatrix_soil(c,3,j) = 1._r8+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                   bmatrix_soil(c,2,j) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                else if (j == 1) then
                   ! this is the snow/soil interface layer
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))
                   bmatrix_soil(c,2,j) = - (1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                   bmatrix_soil(c,3,j) = 1._r8 + (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp &
                        + frac_sno_eff(c) * tk(c,j-1)/dzm) &
                        - (1._r8 - frac_sno_eff(c))*fact(c,j)*dhsdT(c)
                   bmatrix_soil_snow(c,5,j) =   - frac_sno_eff(c) * (1._r8-cnfac) * fact(c,j) &
                        * tk(c,j-1)/dzm
                else if (j <= nlevgrnd-1) then
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))
                   bmatrix_soil(c,2,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                   bmatrix_soil(c,3,j) =   1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                   bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                else if (j == nlevgrnd) then
                   dzm     = (z(c,j)-z(c,j-1))
                   bmatrix_soil(c,2,j) = 0._r8
                   bmatrix_soil(c,3,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                   bmatrix_soil(c,4,j) =   - (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                end if
            end if
         end do
      end do

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)

         ! surface water layer has two coefficients
         dzm=(0.5*dz_h2osfc(c)+col%z(c,1))

         ! diagonal element correction for presence of h2osfc
         if ( frac_h2osfc(c) /= 0.0_r8 ) then
            bmatrix_soil(c,3,1)=bmatrix_soil(c,3,1)+ frac_h2osfc(c) &
                 *((1._r8-cnfac)*fact(c,1)*tk_h2osfc(c)/dzm + fact(c,1)*dhsdT(c))
         end if

      enddo

    end associate

  end subroutine SetMatrix_Soil

  !-----------------------------------------------------------------------
  subroutine SetMatrix_StandingSurfaceWater(bounds, num_nolakec, filter_nolakec, dtime, nband, &
       dhsdT, tk, tk_h2osfc, fact, c_h2osfc, dz_h2osfc, frac_h2osfc, bmatrix_ssw , bmatrix_ssw_soil, bmatrix_soil_ssw)
    !
    ! !DESCRIPTION:
    ! Setup the matrix entries corresponding to internal standing water layers,
    ! soil_standing surface water and standing surface water-soil interaction layers
    !
    ! !USES:
    use clm_varcon     , only : cnfac
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar     , only : nlevsno, nlevmaxurbgrnd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                           ! bounds
    integer , intent(in)  :: num_nolakec                              ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                        ! column filter for non-lake points
    real(r8), intent(in)  :: dtime                                    ! land model time step [sec]
    integer , intent(in)  :: nband                                    ! number of bands of the tridigonal matrix
    real(r8), intent(in)  :: dhsdT(bounds%begc: )                     ! temperature derivative of "hs" [col]
    real(r8), intent(in)  :: tk(bounds%begc: ,-nlevsno+1: )           ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: tk_h2osfc(bounds%begc: )                 ! thermal conductivity [W/(m K)]
    real(r8), intent(in)  :: fact( bounds%begc: , -nlevsno+1: )       ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(in)  :: c_h2osfc( bounds%begc: )                 ! heat capacity of surface water [col]
    real(r8), intent(in)  :: dz_h2osfc(bounds%begc: )                 ! Thickness of standing water [m]
    real(r8), intent(in)  :: frac_h2osfc(bounds%begc: )               ! fractional area with surface water greater than zero
    real(r8), intent(out) :: bmatrix_ssw(bounds%begc: , 1:, 0: )      ! matrix enteries for internal standing water layer
    real(r8), intent(out) :: bmatrix_ssw_soil(bounds%begc: , 1: ,0: ) ! matrix enteries for standing surface water-soil layer interaction
    real(r8), intent(out) :: bmatrix_soil_ssw(bounds%begc: , 1:, 1: ) ! matrix enteries for soil layer-standing surface water interaction
    !
    ! !LOCAL VARIABLES:
    integer  :: c                                                     ! indices
    integer  :: fc                                                    ! lake filtered column indices
    real(r8) :: dzm                                                   ! used in computing tridiagonal matrix
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(dhsdT)            == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk)               == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(tk_h2osfc)        == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(fact)             == (/bounds%endc, nlevmaxurbgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(c_h2osfc)         == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dz_h2osfc)        == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(frac_h2osfc)      == (/bounds%endc/)),           sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_ssw)      == (/bounds%endc, nband, 0/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_ssw_soil) == (/bounds%endc, nband, 0/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(bmatrix_soil_ssw) == (/bounds%endc, nband, 1/)), sourcefile, __LINE__)
    !-----------------------------------------------------------------------

    ! Initialize
    bmatrix_ssw      (bounds%begc:bounds%endc, :, :) = 0.0_r8
    bmatrix_ssw_soil (bounds%begc:bounds%endc, :, :) = 0.0_r8
    bmatrix_soil_ssw (bounds%begc:bounds%endc, :, :) = 0.0_r8

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)

       ! surface water layer has two coefficients
       dzm=(0.5*dz_h2osfc(c)+col%z(c,1))

       bmatrix_ssw(c,3,0)= 1._r8+(1._r8-cnfac)*(dtime/c_h2osfc(c)) &
            *tk_h2osfc(c)/dzm -(dtime/c_h2osfc(c))*dhsdT(c) !interaction from atm

       bmatrix_ssw_soil(c,2,0)= -(1._r8-cnfac)*(dtime/c_h2osfc(c))*tk_h2osfc(c)/dzm !flux to top soil layer

       ! top soil layer has sub coef shifted to 2nd super diagonal
       if ( frac_h2osfc(c) /= 0.0_r8 )then
          bmatrix_soil_ssw(c,4,1)=  - frac_h2osfc(c) * (1._r8-cnfac) * fact(c,1) &
               * tk_h2osfc(c)/dzm !flux from h2osfc
       end if
    enddo

  end subroutine SetMatrix_StandingSurfaceWater

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: BuildingHAC
  !
  ! !INTERFACE:
  subroutine BuildingHAC( bounds, num_urbanl, filter_urbanl, temperature_inst, urbanparams_inst, urbantv_inst, &
                          cool_on, heat_on )
    ! !DESCRIPTION:
    !    Simpler method to manage building temperature (first introduced in CLM4.5). Restricts building
    !    temperature to within bounds, and determine's if heating or cooling is on.
    ! !USES:
    implicit none
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds           ! bounds
    integer          , intent(in) :: num_urbanl       ! number of urban landunits in clump
    integer          , intent(in) :: filter_urbanl(:) ! urban landunit filter
    type(temperature_type)  , intent(inout) :: temperature_inst ! Temperature variables
    type(urbanparams_type)  , intent(in)    :: urbanparams_inst ! urban parameters
    type(urbantv_type)      , intent(in)    :: urbantv_inst     ! urban parameters
    logical, intent(out)  :: cool_on(bounds%begl:)            ! is urban air conditioning on?
    logical, intent(out)  :: heat_on(bounds%begl:)            ! is urban heating on?
    !-----------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer  :: fl,l                       ! indices
    !EOP
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(cool_on)  == (/bounds%endl/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(heat_on)  == (/bounds%endl/)), sourcefile, __LINE__)

    associate(& 
    urbpoi         => lun%urbpoi                         , & ! Input:  [logical (:)]  true => landunit is an urban point       

    t_building     => temperature_inst%t_building_lun    , & ! Input:  [real(r8) (:)]  internal building air temperature [K]       

    t_building_max => urbantv_inst%t_building_max        , & ! Input:  [real(r8) (:)]  maximum internal building air temperature [K]
    t_building_min => urbanparams_inst%t_building_min      & ! Input:  [real(r8) (:)]  minimum internal building air temperature [K]
    )
    ! Restrict internal building temperature to between min and max
    ! and determine if heating or air conditioning is on
    do fl = 1,num_urbanl
       l = filter_urbanl(fl)
       if (urbpoi(l)) then
          cool_on(l) = .false. 
          heat_on(l) = .false. 
          if (t_building(l) > t_building_max(l)) then
            t_building(l) = t_building_max(l)
            cool_on(l) = .true.
            heat_on(l) = .false.
          else if (t_building(l) < t_building_min(l)) then
            t_building(l) = t_building_min(l)
            cool_on(l) = .false.
            heat_on(l) = .true.
          end if
      end if
    end do

    end associate 

  end subroutine BuildingHAC

  !-----------------------------------------------------------------------

end module SoilTemperatureMod
