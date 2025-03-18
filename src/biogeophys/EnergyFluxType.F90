module EnergyFluxType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! Energy flux data structure
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : spval
  use clm_varctl     , only : use_biomass_heat_storage, iulog
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch                
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  !
  implicit none
  save
  private
  !
  type, public :: energyflux_type

     ! Fluxes
     real(r8), pointer :: eflx_sh_stem_patch      (:)   ! patch sensible heat flux from stem (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_h2osfc_to_snow_col (:)   ! col snow melt to h2osfc heat flux (W/m**2)
     real(r8), pointer :: eflx_sh_grnd_patch      (:)   ! patch sensible heat flux from ground (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_veg_patch       (:)   ! patch sensible heat flux from leaves (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_snow_patch      (:)   ! patch sensible heat flux from snow (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_soil_patch      (:)   ! patch sensible heat flux from soil  (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_h2osfc_patch    (:)   ! patch sensible heat flux from surface water (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_patch       (:)   ! patch total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_u_patch     (:)   ! patch urban total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_r_patch     (:)   ! patch rural total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_precip_conversion_col(:) ! col sensible heat flux from precipitation conversion (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_patch       (:)   ! patch total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_u_patch     (:)   ! patch urban total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_r_patch     (:)   ! patch rural total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_vegt_patch      (:)   ! patch transpiration heat flux from veg (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_vege_patch      (:)   ! patch evaporation heat flux from veg (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_grnd_patch      (:)   ! patch evaporation heat flux from ground (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_soil_grnd_patch    (:)   ! patch soil heat flux (W/m**2) [+ = into soil] 
     real(r8), pointer :: eflx_soil_grnd_u_patch  (:)   ! patch urban soil heat flux (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_soil_grnd_r_patch  (:)   ! patch rural soil heat flux (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_lwrad_net_patch    (:)   ! patch net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_net_r_patch  (:)   ! patch rural net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_net_u_patch  (:)   ! patch urban net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_out_patch    (:)   ! patch emitted infrared (longwave) radiation (W/m**2)
     real(r8), pointer :: eflx_lwrad_out_r_patch  (:)   ! patch rural emitted infrared (longwave) rad (W/m**2)
     real(r8), pointer :: eflx_lwrad_out_u_patch  (:)   ! patch urban emitted infrared (longwave) rad (W/m**2)
     real(r8), pointer :: eflx_snomelt_col        (:)   ! col snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_snomelt_r_col      (:)   ! col rural snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_snomelt_u_col      (:)   ! col urban snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_gnet_patch         (:)   ! patch net heat flux into ground  (W/m**2)
     real(r8), pointer :: eflx_grnd_lake_patch    (:)   ! patch net heat flux into lake / snow surface, excluding light transmission (W/m**2)
     real(r8), pointer :: eflx_dynbal_grc         (:)   ! grc dynamic land cover change conversion energy flux (W/m**2)
     real(r8), pointer :: eflx_bot_col            (:)   ! col heat flux from beneath the soil or ice column (W/m**2)
     real(r8), pointer :: eflx_fgr12_col          (:)   ! col ground heat flux between soil layers 1 and 2 (W/m**2)
     real(r8), pointer :: eflx_fgr_col            (:,:) ! col (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)
     real(r8), pointer :: eflx_building_heat_errsoi_col(:) ! col heat flux to interior surface of walls and roof for errsoi check (W m-2)
     real(r8), pointer :: eflx_urban_ac_col       (:)   ! col urban air conditioning flux (W/m**2)
     real(r8), pointer :: eflx_urban_heat_col     (:)   ! col urban heating flux (W/m**2)
     real(r8), pointer :: eflx_anthro_patch       (:)   ! patch total anthropogenic heat flux (W/m**2)
     real(r8), pointer :: eflx_traffic_patch      (:)   ! patch traffic sensible heat flux (W/m**2)
     real(r8), pointer :: eflx_wasteheat_patch    (:)   ! patch sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
     real(r8), pointer :: eflx_ventilation_patch  (:)   ! patch sensible heat flux from building ventilation (W/m**2)
     real(r8), pointer :: eflx_heat_from_ac_patch (:)   ! patch sensible heat flux put back into canyon due to removal by AC (W/m**2)
     real(r8), pointer :: eflx_traffic_lun        (:)   ! lun traffic sensible heat flux (W/m**2)
     real(r8), pointer :: eflx_wasteheat_lun      (:)   ! lun sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
     real(r8), pointer :: eflx_ventilation_lun    (:)   ! lun sensible heat flux from building ventilation (W/m**2)
     real(r8), pointer :: eflx_heat_from_ac_lun   (:)   ! lun sensible heat flux to be put back into canyon due to removal by AC (W/m**2)
     real(r8), pointer :: eflx_building_lun       (:)   ! lun building heat flux from change in interior building air temperature (W/m**2)
     real(r8), pointer :: eflx_urban_ac_lun       (:)   ! lun urban air conditioning flux (W/m**2)
     real(r8), pointer :: eflx_urban_heat_lun     (:)   ! lun urban heating flux (W/m**2)

     ! Derivatives of energy fluxes
     real(r8), pointer :: dgnetdT_patch           (:)   ! patch derivative of net ground heat flux wrt soil temp  (W/m**2 K)
     real(r8), pointer :: netrad_patch            (:)   ! col net radiation (W/m**2) [+ = to sfc]
     real(r8), pointer :: cgrnd_patch             (:)   ! col deriv. of soil energy flux wrt to soil temp [W/m2/k]
     real(r8), pointer :: cgrndl_patch            (:)   ! col deriv. of soil latent heat flux wrt soil temp  [W/m**2/k]
     real(r8), pointer :: cgrnds_patch            (:)   ! col deriv. of soil sensible heat flux wrt soil temp [W/m2/k]

     ! Canopy radiation
     real(r8), pointer :: dlrad_patch             (:)   ! col downward longwave radiation below the canopy [W/m2]
     real(r8), pointer :: ulrad_patch             (:)   ! col upward longwave radiation above the canopy [W/m2]

     ! Wind Stress
     real(r8), pointer :: taux_patch              (:)   ! patch wind (shear) stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy_patch              (:)   ! patch wind (shear) stress: n-s (kg/m/s**2)

     ! Conductance
     real(r8), pointer :: canopy_cond_patch       (:)   ! patch tracer conductance for canopy [m/s] 

     ! Transpiration
     real(r8), pointer :: btran_patch             (:)   ! patch transpiration wetness factor (0 to 1)
     real(r8), pointer :: btran_min_patch         (:)   ! patch daily minimum transpiration wetness factor (0 to 1)
     real(r8), pointer :: btran_min_inst_patch    (:)   ! patch instantaneous daily minimum transpiration wetness factor (0 to 1)
     real(r8), pointer :: bsun_patch              (:)   ! patch sunlit canopy transpiration wetness factor (0 to 1)
     real(r8), pointer :: bsha_patch              (:)   ! patch shaded canopy transpiration wetness factor (0 to 1)

     ! Roots
     real(r8), pointer :: rresis_patch            (:,:) ! patch root resistance by layer (0-1)  (nlevgrnd)

     ! Latent heat
     real(r8), pointer :: htvp_col                (:)   ! latent heat of vapor of water (or sublimation) [j/kg]

     ! Canopy heat
     real(r8), pointer :: dhsdt_canopy_patch      (:)   ! patch change in heat content of canopy (leaf+stem) (W/m**2) [+ to atm]

     ! Balance Checks
     real(r8), pointer :: errsoi_patch            (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errsoi_col              (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errseb_patch            (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errseb_col              (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errsol_patch            (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errsol_col              (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errlon_patch            (:)   ! longwave radiation conservation error (W/m**2)
     real(r8), pointer :: errlon_col              (:)   ! longwave radiation conservation error (W/m**2)

     ! Objects that help convert once-per-year dynamic land cover changes into fluxes
     ! that are dribbled throughout the year
     type(annual_flux_dribbler_type) :: eflx_dynbal_dribbler

   contains

     procedure, public  :: Init            ! Public initialization method
     procedure, private :: InitAllocate    ! initialize/allocate
     procedure, private :: InitHistory     ! setup history fields
     procedure, private :: InitCold        ! initialize for cold start
     procedure, public  :: Restart         ! setup restart fields
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars

  end type energyflux_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, t_grnd_col, is_simple_buildtemp, is_prog_buildtemp )
    !
    ! !DESCRIPTION:
    !    Allocate and initialize the data type and setup history, and initialize for cold-start.
    ! !USES:
    implicit none
    ! !ARGUMENTS:
    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )
    logical           , intent(in) :: is_simple_buildtemp        ! If using simple building temp method
    logical           , intent(in) :: is_prog_buildtemp          ! If using prognostic building temp method

    SHR_ASSERT_ALL_FL((ubound(t_grnd_col) == (/bounds%endc/)), sourcefile, __LINE__)

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds, is_simple_buildtemp, is_prog_buildtemp )
    call this%InitCold ( bounds, t_grnd_col, is_simple_buildtemp, is_prog_buildtemp ) 

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevgrnd
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate( this%eflx_h2osfc_to_snow_col (begc:endc))             ; this%eflx_h2osfc_to_snow_col (:)   = nan
    allocate( this%eflx_sh_snow_patch      (begp:endp))             ; this%eflx_sh_snow_patch      (:)   = nan
    allocate( this%eflx_sh_soil_patch      (begp:endp))             ; this%eflx_sh_soil_patch      (:)   = nan
    allocate( this%eflx_sh_h2osfc_patch    (begp:endp))             ; this%eflx_sh_h2osfc_patch    (:)   = nan
    allocate( this%eflx_sh_tot_patch       (begp:endp))             ; this%eflx_sh_tot_patch       (:)   = nan
    allocate( this%eflx_sh_tot_u_patch     (begp:endp))             ; this%eflx_sh_tot_u_patch     (:)   = nan
    allocate( this%eflx_sh_tot_r_patch     (begp:endp))             ; this%eflx_sh_tot_r_patch     (:)   = nan
    allocate( this%eflx_sh_grnd_patch      (begp:endp))             ; this%eflx_sh_grnd_patch      (:)   = nan
    allocate( this%eflx_sh_stem_patch      (begp:endp))             ; this%eflx_sh_stem_patch      (:)   = nan
    allocate( this%eflx_sh_veg_patch       (begp:endp))             ; this%eflx_sh_veg_patch       (:)   = nan
    allocate( this%eflx_sh_precip_conversion_col(begc:endc))        ; this%eflx_sh_precip_conversion_col(:) = nan
    allocate( this%eflx_lh_tot_u_patch     (begp:endp))             ; this%eflx_lh_tot_u_patch     (:)   = nan
    allocate( this%eflx_lh_tot_patch       (begp:endp))             ; this%eflx_lh_tot_patch       (:)   = nan
    allocate( this%eflx_lh_tot_r_patch     (begp:endp))             ; this%eflx_lh_tot_r_patch     (:)   = nan
    allocate( this%eflx_lh_grnd_patch      (begp:endp))             ; this%eflx_lh_grnd_patch      (:)   = nan
    allocate( this%eflx_lh_vege_patch      (begp:endp))             ; this%eflx_lh_vege_patch      (:)   = nan
    allocate( this%eflx_lh_vegt_patch      (begp:endp))             ; this%eflx_lh_vegt_patch      (:)   = nan
    allocate( this%eflx_soil_grnd_patch    (begp:endp))             ; this%eflx_soil_grnd_patch    (:)   = nan
    allocate( this%eflx_soil_grnd_u_patch  (begp:endp))             ; this%eflx_soil_grnd_u_patch  (:)   = nan
    allocate( this%eflx_soil_grnd_r_patch  (begp:endp))             ; this%eflx_soil_grnd_r_patch  (:)   = nan
    allocate( this%eflx_lwrad_net_patch    (begp:endp))             ; this%eflx_lwrad_net_patch    (:)   = nan
    allocate( this%eflx_lwrad_net_u_patch  (begp:endp))             ; this%eflx_lwrad_net_u_patch  (:)   = nan
    allocate( this%eflx_lwrad_net_r_patch  (begp:endp))             ; this%eflx_lwrad_net_r_patch  (:)   = nan
    allocate( this%eflx_lwrad_out_patch    (begp:endp))             ; this%eflx_lwrad_out_patch    (:)   = nan
    allocate( this%eflx_lwrad_out_u_patch  (begp:endp))             ; this%eflx_lwrad_out_u_patch  (:)   = nan
    allocate( this%eflx_lwrad_out_r_patch  (begp:endp))             ; this%eflx_lwrad_out_r_patch  (:)   = nan
    allocate( this%eflx_gnet_patch         (begp:endp))             ; this%eflx_gnet_patch         (:)   = nan
    allocate( this%eflx_grnd_lake_patch    (begp:endp))             ; this%eflx_grnd_lake_patch    (:)   = nan
    allocate( this%eflx_dynbal_grc         (begg:endg))             ; this%eflx_dynbal_grc         (:)   = nan
    allocate( this%eflx_bot_col            (begc:endc))             ; this%eflx_bot_col            (:)   = nan
    allocate( this%eflx_snomelt_col        (begc:endc))             ; this%eflx_snomelt_col        (:)   = nan
    allocate( this%eflx_snomelt_r_col      (begc:endc))             ; this%eflx_snomelt_r_col      (:)   = nan
    allocate( this%eflx_snomelt_u_col      (begc:endc))             ; this%eflx_snomelt_u_col      (:)   = nan
    allocate( this%eflx_fgr12_col          (begc:endc))             ; this%eflx_fgr12_col          (:)   = nan
    allocate( this%eflx_fgr_col            (begc:endc, 1:nlevgrnd)) ; this%eflx_fgr_col            (:,:) = nan
    allocate( this%eflx_building_heat_errsoi_col  (begc:endc))      ; this%eflx_building_heat_errsoi_col(:)= nan
    allocate( this%eflx_urban_ac_col       (begc:endc))             ; this%eflx_urban_ac_col       (:)   = nan
    allocate( this%eflx_urban_heat_col     (begc:endc))             ; this%eflx_urban_heat_col     (:)   = nan
    allocate( this%eflx_wasteheat_patch    (begp:endp))             ; this%eflx_wasteheat_patch    (:)   = nan
    allocate( this%eflx_ventilation_patch  (begp:endp))             ; this%eflx_ventilation_patch  (:)   = nan
    allocate( this%eflx_traffic_patch      (begp:endp))             ; this%eflx_traffic_patch      (:)   = nan
    allocate( this%eflx_heat_from_ac_patch (begp:endp))             ; this%eflx_heat_from_ac_patch (:)   = nan
    allocate( this%eflx_heat_from_ac_lun   (begl:endl))             ; this%eflx_heat_from_ac_lun   (:)   = nan
    allocate( this%eflx_building_lun       (begl:endl))             ; this%eflx_building_lun       (:)   = nan
    allocate( this%eflx_urban_ac_lun       (begl:endl))             ; this%eflx_urban_ac_lun       (:)   = nan
    allocate( this%eflx_urban_heat_lun     (begl:endl))             ; this%eflx_urban_heat_lun     (:)   = nan
    allocate( this%eflx_traffic_lun        (begl:endl))             ; this%eflx_traffic_lun        (:)   = nan
    allocate( this%eflx_wasteheat_lun      (begl:endl))             ; this%eflx_wasteheat_lun      (:)   = nan
    allocate( this%eflx_ventilation_lun    (begl:endl))             ; this%eflx_ventilation_lun    (:)   = nan
    allocate( this%eflx_anthro_patch       (begp:endp))             ; this%eflx_anthro_patch       (:)   = nan

    allocate( this%dgnetdT_patch           (begp:endp))             ; this%dgnetdT_patch           (:)   = nan
    allocate( this%cgrnd_patch             (begp:endp))             ; this%cgrnd_patch             (:)   = nan
    allocate( this%cgrndl_patch            (begp:endp))             ; this%cgrndl_patch            (:)   = nan
    allocate( this%cgrnds_patch            (begp:endp))             ; this%cgrnds_patch            (:)   = nan
    allocate( this%dlrad_patch             (begp:endp))             ; this%dlrad_patch             (:)   = nan
    allocate( this%ulrad_patch             (begp:endp))             ; this%ulrad_patch             (:)   = nan
    allocate( this%netrad_patch            (begp:endp))             ; this%netrad_patch            (:)   = nan  

    allocate( this%taux_patch              (begp:endp))             ; this%taux_patch              (:)   = nan
    allocate( this%tauy_patch              (begp:endp))             ; this%tauy_patch              (:)   = nan

    allocate( this%canopy_cond_patch       (begp:endp))             ; this%canopy_cond_patch       (:)   = nan

    allocate( this%htvp_col                (begc:endc))             ; this%htvp_col                (:)   = nan

    allocate( this%dhsdt_canopy_patch      (begp:endp))             ; this%dhsdt_canopy_patch      (:)   = nan
 
    allocate(this%rresis_patch             (begp:endp,1:nlevgrnd))  ; this%rresis_patch            (:,:) = nan
    allocate(this%btran_patch              (begp:endp))             ; this%btran_patch             (:)   = nan
    allocate(this%btran_min_patch          (begp:endp))             ; this%btran_min_patch         (:)   = nan
    allocate(this%btran_min_inst_patch     (begp:endp))             ; this%btran_min_inst_patch    (:)   = nan
    allocate( this%bsun_patch              (begp:endp))             ; this%bsun_patch              (:)   = nan
    allocate( this%bsha_patch              (begp:endp))             ; this%bsha_patch              (:)   = nan
    allocate( this%errsoi_patch            (begp:endp))             ; this%errsoi_patch            (:)   = nan
    allocate( this%errsoi_col              (begc:endc))             ; this%errsoi_col              (:)   = nan
    allocate( this%errseb_patch            (begp:endp))             ; this%errseb_patch            (:)   = nan
    allocate( this%errseb_col              (begc:endc))             ; this%errseb_col              (:)   = nan
    allocate( this%errsol_patch            (begp:endp))             ; this%errsol_patch            (:)   = nan
    allocate( this%errsol_col              (begc:endc))             ; this%errsol_col              (:)   = nan
    allocate( this%errlon_patch            (begp:endp))             ; this%errlon_patch            (:)   = nan
    allocate( this%errlon_col              (begc:endc))             ; this%errlon_col              (:)   = nan

    this%eflx_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = 'eflx_dynbal', &
         units = 'J/m**2')

  end subroutine InitAllocate
    
  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, is_simple_buildtemp, is_prog_buildtemp)
    !
    ! !DESCRIPTION:
    ! Setup fields that can be output to history files
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevgrnd
    use clm_varctl     , only : use_cn, use_hydrstress
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    use ncdio_pio      , only : ncd_inqvdlen
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    logical          , intent(in) :: is_simple_buildtemp ! If using simple building temp method
    logical          , intent(in) :: is_prog_buildtemp   ! If using prognostic building temp method
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begl, endl
    integer           :: begg, endg
    integer           :: dimlen
    integer           :: err_code
    logical           :: do_io
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg


    this%eflx_dynbal_grc(begg:endg) = spval 
    call hist_addfld1d (fname='EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_lnd=this%eflx_dynbal_grc)

    this%eflx_snomelt_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM',  units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=this%eflx_snomelt_col, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSM_ICE', units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux (ice landunits only)', &
         ptr_col=this%eflx_snomelt_col, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_snomelt_r_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM_R',  units='W/m^2',  &
         avgflag='A', long_name='Rural snow melt heat flux', &
         ptr_col=this%eflx_snomelt_r_col, set_spec=spval, default='inactive')

    this%eflx_snomelt_u_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM_U',  units='W/m^2',  &
         avgflag='A', long_name='Urban snow melt heat flux', &
         ptr_col=this%eflx_snomelt_u_col, c2l_scale_type='urbanf', set_nourb=spval, default='inactive')

    this%eflx_lwrad_net_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRA', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_patch, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FIRA_ICE', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation (ice landunits only)', &
         ptr_patch=this%eflx_lwrad_net_patch, c2l_scale_type='urbanf', l2g_scale_type='ice',&
         default='inactive')

    this%eflx_lwrad_net_r_patch(begp:endp) = spval 
    call hist_addfld1d (fname='FIRA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_r_patch, set_spec=spval)

    this%eflx_lwrad_out_patch(begp:endp) = spval 
    call hist_addfld1d (fname='FIRE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_patch, c2l_scale_type='urbanf')
    ! Rename of FIRE for Urban intercomparision project
    call hist_addfld1d (fname='LWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling longwave radiation', &
         ptr_patch=this%eflx_lwrad_out_patch, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FIRE_ICE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation (ice landunits only)', &
         ptr_patch=this%eflx_lwrad_out_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_lwrad_out_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_R', units='W/m^2',  &
         avgflag='A', long_name='Rural emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_r_patch, set_spec=spval)

    this%eflx_lh_vegt_patch(begp:endp) = spval
    call hist_addfld1d (fname='FCTR', units='W/m^2',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_patch=this%eflx_lh_vegt_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_lh_vege_patch(begp:endp) = spval
    call hist_addfld1d (fname='FCEV', units='W/m^2',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_patch=this%eflx_lh_vege_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_lh_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGEV', units='W/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_patch=this%eflx_lh_grnd_patch, c2l_scale_type='urbanf') 

    this%eflx_sh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH', units='W/m^2',  &
         avgflag='A', long_name='sensible heat not including correction for land use change and rain/snow conversion', &
         ptr_patch=this%eflx_sh_tot_patch, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FSH_ICE', units='W/m^2',  &
         avgflag='A', &
         long_name='sensible heat not including correction for land use change and rain/snow conversion (ice landunits only)', &
         ptr_patch=this%eflx_sh_tot_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_sh_tot_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_R', units='W/m^2',  &
         avgflag='A', long_name='Rural sensible heat', &
         ptr_patch=this%eflx_sh_tot_r_patch, set_spec=spval)

    this%eflx_sh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qh', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_patch=this%eflx_sh_tot_patch, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qle', units='W/m^2',  &
         avgflag='A', long_name='total evaporation', &
         ptr_patch=this%eflx_lh_tot_patch, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT', units='W/m^2', &
         avgflag='A', long_name='total latent heat flux [+ to atm]', &
         ptr_patch=this%eflx_lh_tot_patch, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='EFLX_LH_TOT_ICE', units='W/m^2',  &
         avgflag='A', long_name='total latent heat flux [+ to atm] (ice landunits only)', &
         ptr_patch=this%eflx_lh_tot_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_lh_tot_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT_R', units='W/m^2',  &
         avgflag='A', long_name='Rural total evaporation', &
         ptr_patch=this%eflx_lh_tot_r_patch, set_spec=spval)

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qstor', units='W/m^2',  &
         avgflag='A', long_name='storage heat flux (includes snowmelt)', &
         ptr_patch=this%eflx_soil_grnd_patch, c2l_scale_type='urbanf', &
         default = 'inactive')
    this%eflx_sh_veg_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_V', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from veg', &
         ptr_patch=this%eflx_sh_veg_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    if (use_biomass_heat_storage) then
       this%eflx_sh_stem_patch(begp:endp) = spval
       call hist_addfld1d (fname='FSH_STEM', units='W/m^2',  &
            avgflag='A', long_name='sensible heat from stem', &
            ptr_patch=this%eflx_sh_stem_patch, c2l_scale_type='urbanf',default = 'inactive')

       this%dhsdt_canopy_patch(begp:endp) = spval
       call hist_addfld1d (fname='DHSDT_CANOPY', units='W/m^2',  &
            avgflag='A', long_name='change in canopy heat storage', &
            ptr_patch=this%dhsdt_canopy_patch, set_lake=0._r8, c2l_scale_type='urbanf',default='active')
    endif

    this%eflx_sh_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_G', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_patch=this%eflx_sh_grnd_patch, c2l_scale_type='urbanf')

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGR', units='W/m^2',  &
         avgflag='A', long_name='heat flux into soil/snow including snow melt and lake / snow light transmission', &
         ptr_patch=this%eflx_soil_grnd_patch, c2l_scale_type='urbanf')

    call hist_addfld1d (fname='FGR_ICE', units='W/m^2',  &
         avgflag='A', &
         long_name='heat flux into soil/snow including snow melt and lake / snow light transmission (ice landunits only)', &
         ptr_patch=this%eflx_soil_grnd_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_soil_grnd_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGR_R', units='W/m^2',  &
         avgflag='A', long_name='Rural heat flux into soil/snow including snow melt and snow light transmission', &
         ptr_patch=this%eflx_soil_grnd_r_patch, set_spec=spval, default='inactive')

    this%eflx_lwrad_net_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_u_patch, c2l_scale_type='urbanf', set_nourb=spval, default='inactive')

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_SOIL_GRND', units='W/m^2', &
         avgflag='A', long_name='soil heat flux [+ into soil]', &
         ptr_patch=this%eflx_soil_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%eflx_lwrad_out_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_U', units='W/m^2',  &
         avgflag='A', long_name='Urban emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_u_patch, c2l_scale_type='urbanf', set_nourb=spval, default='inactive')

    this%eflx_sh_tot_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_U', units='W/m^2',  &
         avgflag='A', long_name='Urban sensible heat', &
         ptr_patch=this%eflx_sh_tot_u_patch, c2l_scale_type='urbanf', set_nourb=spval, default='inactive')

    this%eflx_sh_precip_conversion_col(begc:endc) = spval
    call hist_addfld1d (fname = 'FSH_PRECIP_CONVERSION', units='W/m^2', &
         avgflag='A', long_name='Sensible heat flux from conversion of rain/snow atm forcing', &
         ptr_col=this%eflx_sh_precip_conversion_col, c2l_scale_type='urbanf')

    this%eflx_lh_tot_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT_U', units='W/m^2',  &
         avgflag='A', long_name='Urban total evaporation', &
         ptr_patch=this%eflx_lh_tot_u_patch, c2l_scale_type='urbanf', set_nourb=spval, default='inactive')

    this%eflx_soil_grnd_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGR_U', units='W/m^2',  &
         avgflag='A', long_name='Urban heat flux into soil/snow including snow melt', &
         ptr_patch=this%eflx_soil_grnd_u_patch, c2l_scale_type='urbanf', set_nourb=spval, default='inactive')

    this%netrad_patch(begp:endp) = spval
    call hist_addfld1d (fname='Rnet', units='W/m^2',  &
         avgflag='A', long_name='net radiation', &
         ptr_patch=this%netrad_patch, c2l_scale_type='urbanf', &
         default='inactive')

    if (use_cn) then
       this%dlrad_patch(begp:endp) = spval
       call hist_addfld1d (fname='DLRAD', units='W/m^2', &
            avgflag='A', long_name='downward longwave radiation below the canopy', &
            ptr_patch=this%dlrad_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%ulrad_patch(begp:endp) = spval
       call hist_addfld1d (fname='ULRAD', units='W/m^2', &
            avgflag='A', long_name='upward longwave radiation above the canopy', &
            ptr_patch=this%ulrad_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%cgrnd_patch(begp:endp) = spval
       call hist_addfld1d (fname='CGRND', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil energy flux wrt to soil temp', &
            ptr_patch=this%cgrnd_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%cgrndl_patch(begp:endp) = spval
       call hist_addfld1d (fname='CGRNDL', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil latent heat flux wrt soil temp', &
            ptr_patch=this%cgrndl_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%cgrnds_patch(begp:endp) = spval
       call hist_addfld1d (fname='CGRNDS', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil sensible heat flux wrt soil temp', &
            ptr_patch=this%cgrnds_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    this%eflx_gnet_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_GNET', units='W/m^2', &
         avgflag='A', long_name='net heat flux into ground', &
         ptr_patch=this%eflx_gnet_patch, default='inactive', c2l_scale_type='urbanf')

    this%eflx_grnd_lake_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_GRND_LAKE', units='W/m^2', &
         avgflag='A', long_name='net heat flux into lake/snow surface, excluding light transmission', &
         ptr_patch=this%eflx_grnd_lake_patch, set_nolake=spval)

    if (      is_simple_buildtemp )then
       this%eflx_building_heat_errsoi_col(begc:endc) = spval
       call hist_addfld1d (fname='BUILDHEAT', units='W/m^2',  &
            avgflag='A', long_name='heat flux from urban building interior to walls and roof', &
            ptr_col=this%eflx_building_heat_errsoi_col, set_nourb=0._r8, c2l_scale_type='urbanf')

       this%eflx_urban_ac_col(begc:endc) = spval
       call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
            avgflag='A', long_name='urban air conditioning flux', &
            ptr_col=this%eflx_urban_ac_col, set_nourb=0._r8, c2l_scale_type='urbanf')
   
       this%eflx_urban_heat_col(begc:endc) = spval
       call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
            avgflag='A', long_name='urban heating flux', &
            ptr_col=this%eflx_urban_heat_col, set_nourb=0._r8, c2l_scale_type='urbanf')
    else
       this%eflx_building_lun(begl:endl) = spval
       call hist_addfld1d (fname='EFLXBUILD', units='W/m^2',  &
            avgflag='A', long_name='building heat flux from change in interior building air temperature', &
            ptr_lunit=this%eflx_building_lun, set_nourb=0._r8, l2g_scale_type='unity')

       this%eflx_urban_ac_lun(begl:endl) = spval
       call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
            avgflag='A', long_name='urban air conditioning flux', &
            ptr_lunit=this%eflx_urban_ac_lun, set_nourb=0._r8, l2g_scale_type='unity')

       this%eflx_urban_heat_lun(begl:endl) = spval
       call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
            avgflag='A', long_name='urban heating flux', &
            ptr_lunit=this%eflx_urban_heat_lun, set_nourb=0._r8, l2g_scale_type='unity')
    end if


    this%dgnetdT_patch(begp:endp) = spval
    call hist_addfld1d (fname='DGNETDT', units='W/m^2/K', &
         avgflag='A', long_name='derivative of net ground heat flux wrt soil temp', &
         ptr_patch=this%dgnetdT_patch, default='inactive', c2l_scale_type='urbanf')

    this%eflx_fgr12_col(begc:endc) = spval
    call hist_addfld1d (fname='FGR12',  units='W/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=this%eflx_fgr12_col, set_lake=spval)

    this%eflx_fgr_col(begc:endc,:) = spval
    call hist_addfld2d (fname='FGR_SOIL_R', units='watt/m^2', type2d='levgrnd', &
         avgflag='A', long_name='Rural downward heat flux at interface below each soil layer', &
         ptr_col=this%eflx_fgr_col, set_spec=spval, default='inactive')

    this%eflx_traffic_patch(begp:endp) = spval
    call hist_addfld1d (fname='TRAFFICFLUX', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from urban traffic', &
         ptr_patch=this%eflx_traffic_patch, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    this%eflx_wasteheat_patch(begp:endp) = spval
    call hist_addfld1d (fname='WASTEHEAT', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from heating/cooling sources of urban waste heat', &
         ptr_patch=this%eflx_wasteheat_patch, set_nourb=0._r8, c2l_scale_type='urbanf')

    if ( is_prog_buildtemp )then
       this%eflx_ventilation_patch(begp:endp) = spval
       call hist_addfld1d (fname='VENTILATION', units='W/m^2',  &
            avgflag='A', long_name='sensible heat flux from building ventilation', &
            ptr_patch=this%eflx_ventilation_patch, set_nourb=0._r8, c2l_scale_type='urbanf')
    end if

    this%eflx_heat_from_ac_patch(begp:endp) = spval
    call hist_addfld1d (fname='HEAT_FROM_AC', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux put into canyon due to heat removed from air conditioning', &
         ptr_patch=this%eflx_heat_from_ac_patch, set_nourb=0._r8, c2l_scale_type='urbanf')

    if ( is_simple_buildtemp )then
       this%eflx_anthro_patch(begp:endp) = spval
       call hist_addfld1d (fname='Qanth', units='W/m^2',  &
            avgflag='A', long_name='anthropogenic heat flux', &
            ptr_patch=this%eflx_anthro_patch, set_nourb=0._r8, c2l_scale_type='urbanf', &
            default='inactive')
    end if

    this%taux_patch(begp:endp) = spval
    call hist_addfld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_patch=this%taux_patch)
    ! Rename of TAUX for Urban intercomparision project (when U=V)    
    call hist_addfld1d (fname='Qtau', units='kg/m/s^2',  &  
         avgflag='A', long_name='momentum flux', &
         ptr_patch=this%taux_patch, default='inactive')

    this%tauy_patch(begp:endp) = spval
    call hist_addfld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_patch=this%tauy_patch)

    this%btran_patch(begp:endp) = spval
    if (.not. use_hydrstress) then
       call hist_addfld1d (fname='BTRAN', units='unitless',  &
            avgflag='A', long_name='transpiration beta factor', &
            ptr_patch=this%btran_patch, l2g_scale_type='veg')
    end if

    this%btran_min_patch(begp:endp) = spval
    call hist_addfld1d (fname='BTRANMN', units='unitless',  &
         avgflag='A', long_name='daily minimum of transpiration beta factor', &
         ptr_patch=this%btran_min_patch, l2g_scale_type='veg')

    if (use_cn) then
       this%rresis_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='RRESIS', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='root resistance in each soil layer', &
            ptr_patch=this%rresis_patch, l2g_scale_type='veg', default='inactive')
    end if

    this%errsoi_col(begc:endc) = spval
    call hist_addfld1d (fname='ERRSOI',  units='W/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=this%errsoi_col)

    this%errseb_patch(begp:endp) = spval
    call hist_addfld1d (fname='ERRSEB',  units='W/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_patch=this%errseb_patch)

    this%errsol_patch(begp:endp) = spval
    call hist_addfld1d (fname='ERRSOL',  units='W/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_patch=this%errsol_patch, set_urb=spval)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, t_grnd_col, is_simple_buildtemp, is_prog_buildtemp)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_varpar      , only : nlevgrnd
    use clm_varcon      , only : sb
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )
    logical           , intent(in) :: is_simple_buildtemp ! If using simple building temp method
    logical           , intent(in) :: is_prog_buildtemp   ! If using prognostic building temp method
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p,levs,lev
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(t_grnd_col) == (/bounds%endc/)), sourcefile, __LINE__)

    ! Columns
    if ( is_simple_buildtemp )then
       do c = bounds%begc, bounds%endc
          l = col%landunit(c)

          if (lun%urbpoi(l)) then
             this%eflx_building_heat_errsoi_col(c) = 0._r8
             this%eflx_urban_ac_col(c)             = 0._r8
             this%eflx_urban_heat_col(c)           = 0._r8
          else
             this%eflx_building_heat_errsoi_col(c) = 0._r8
             this%eflx_urban_ac_col(c)             = 0._r8
             this%eflx_urban_heat_col(c)           = 0._r8
          end if

       end do
    end if

    ! Patches
    do p = bounds%begp, bounds%endp 
       c = patch%column(p)
       l = patch%landunit(p)

       if (.not. lun%urbpoi(l)) then ! non-urban
          this%eflx_lwrad_net_u_patch(p) = spval
          this%eflx_lwrad_out_u_patch(p) = spval
          this%eflx_lh_tot_u_patch(p)    = spval
          this%eflx_sh_tot_u_patch(p)    = spval
          this%eflx_soil_grnd_u_patch(p) = spval
       end if

       this%eflx_lwrad_out_patch(p) = sb * (t_grnd_col(c))**4
    end do

    ! patches
    do p = bounds%begp, bounds%endp 
       l = patch%landunit(p)

       if (.not. lun%urbpoi(l)) then
          this%eflx_traffic_lun(l)        = spval
          this%eflx_wasteheat_lun(l)      = spval
          this%eflx_ventilation_lun(l)    = spval
          if ( is_prog_buildtemp )then
             this%eflx_building_lun(l)   = 0._r8
             this%eflx_urban_ac_lun(l)   = 0._r8
             this%eflx_urban_heat_lun(l) = 0._r8
          end if

          this%eflx_wasteheat_patch(p)    = 0._r8
          this%eflx_ventilation_patch(p)  = 0._r8
          this%eflx_heat_from_ac_patch(p) = 0._r8
          this%eflx_traffic_patch(p)      = 0._r8
          if ( is_simple_buildtemp) &
             this%eflx_anthro_patch(p)    = 0._r8
       else
          if ( is_prog_buildtemp )then
             this%eflx_building_lun(l)   = 0._r8
             this%eflx_urban_ac_lun(l)   = 0._r8
             this%eflx_urban_heat_lun(l) = 0._r8
             this%eflx_ventilation_lun(l)= 0._r8
          end if
       end if
    end do

    ! initialize rresis, for use in ecosystemdyn
    do p = bounds%begp,bounds%endp
       do lev = 1,nlevgrnd
          this%rresis_patch(p,lev) = 0._r8
       end do
    end do 

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, is_simple_buildtemp, is_prog_buildtemp)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, &
                            ncd_inqvdlen
    use restUtilMod
    use decompMod      , only : get_proc_global
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    logical          , intent(in)    :: is_simple_buildtemp ! If using simple building temp method
    logical          , intent(in)    :: is_prog_buildtemp   ! If using prognostic building temp method
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    integer :: dimlen
    integer :: err_code
    integer :: numl_global
    logical :: readvar      ! determine if variable is on initial file
    logical :: do_io
    !-----------------------------------------------------------------------

    call get_proc_global(nl=numl_global)
    call restartvar(ncid=ncid, flag=flag, varname='EFLX_LWRAD_OUT', xtype=ncd_double,  & 
         dim1name='pft', &
         long_name='emitted infrared (longwave) radiation', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_lwrad_out_patch)

    ! Restart for building air temperature method
    if ( is_prog_buildtemp )then
       ! landunit urban energy state variable - eflx_urban_ac
       do_io = .true.
       ! On a read, confirm that this variable has the expected size (landunit-level); if not, 
       ! don't read it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size (column-level).
       if (flag == 'read') then
          call ncd_inqvdlen(ncid, 'URBAN_AC_L', 1, dimlen, err_code)
          if (dimlen /= numl_global) then
             do_io = .false.
             readvar = .false.
          end if
       end if
       if (do_io) then
          call restartvar(ncid=ncid, flag=flag, varname='URBAN_AC_L', xtype=ncd_double,  &
               dim1name='landunit',&
               long_name='urban air conditioning flux', units='watt/m^2', &
               interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_ac_lun)
       else
          this%eflx_urban_ac_lun = 0.0_r8
       end if
       ! landunit urban energy state variable - eflx_urban_heat
       do_io = .true.
       ! On a read, confirm that this variable has the expected size (landunit-level); if not, 
       ! don't read it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size (column-level).
       if (flag == 'read') then
          call ncd_inqvdlen(ncid, 'URBAN_HEAT_L', 1, dimlen, err_code)
          if (dimlen /= numl_global) then
             do_io = .false.
             readvar = .false.
          end if
       end if
       if (do_io) then
          call restartvar(ncid=ncid, flag=flag, varname='URBAN_HEAT_L', xtype=ncd_double,  &
               dim1name='landunit',&
               long_name='urban heating flux', units='watt/m^2', &
               interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_heat_lun)
       else
          this%eflx_urban_heat_lun = 0.0_r8
       end if
       call restartvar(ncid=ncid, flag=flag, varname='EFLX_VENTILATION', xtype=ncd_double, &
           dim1name='landunit', &
           long_name='sensible heat flux from building ventilation', units='watt/m^2', &
           interpinic_flag='interp', readvar=readvar, data=this%eflx_ventilation_lun)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "can't find EFLX_VENTILATION in initial file..."
          if (masterproc) write(iulog,*) "Initialize EFLX_VENTILATION to zero"
          this%eflx_ventilation_lun(bounds%begl:bounds%endl) = 0._r8
       end if

    else if ( is_simple_buildtemp )then
        call restartvar(ncid=ncid, flag=flag, varname='URBAN_AC', xtype=ncd_double, &
            dim1name='column', &
            long_name='urban air conditioning flux', units='watt/m^2', &
            interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_ac_col)
        call restartvar(ncid=ncid, flag=flag, varname='URBAN_HEAT', xtype=ncd_double,  &
            dim1name='column', &
            long_name='urban heating flux', units='watt/m^2', &
            interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_heat_col)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='BTRAN_MIN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily minimum of transpiration wetness factor', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran_min_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='BTRAN_MIN_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily minimum of transpiration wetness factor', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran_min_inst_patch) 

    call this%eflx_dynbal_dribbler%Restart(bounds, ncid, flag)

  end subroutine Restart
  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    ! Each interval and accumulation type is unique to each field processed.
    ! Routine [initAccBuffer] defines the fields to be processed
    ! and the type of accumulation. 
    ! Routine [updateAccVars] does the actual accumulation for a given field.
    ! Fields are accumulated by calls to subroutine [update_accum_field]. 
    ! To accumulate a field, it must first be defined in subroutine [initAccVars] 
    ! and then accumulated by calls to [updateAccVars].
    ! Four types of accumulations are possible:
    !   o average over time interval
    !   o running mean over time interval
    !   o running accumulation over time interval
    ! Time average fields are only valid at the end of the averaging interval.
    ! Running means are valid once the length of the simulation exceeds the
    ! averaging interval. Accumulated fields are continuously accumulated.
    ! The trigger value "-99999." resets the accumulation to zero.
    !
    ! !USES 
    use accumulMod       , only : init_accum_field
    use clm_time_manager , only : get_step_size_real
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES: 
    real(r8) :: dtime
    integer, parameter :: not_used = huge(1)
    !---------------------------------------------------------------------

    dtime = get_step_size_real()

    call init_accum_field(name='BTRANAV', units='-', &
         desc='average over an hour of btran', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer
  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use clm_time_manager , only : get_nstep
    use clm_varctl       , only : nsrest, nsrStartup
    use abortutils       , only : endrun
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    ! Initialize variables that are to be time accumulated
    ! Initialize btran min values
    if (nsrest == nsrStartup) then
       this%btran_min_patch(begp:endp)        =  spval

       this%btran_min_inst_patch(begp:endp)   =  spval
    end if

  end subroutine InitAccVars
  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod       , only : update_accum_field, extract_accum_field
    use abortutils       , only : endrun
    !
    ! !ARGUMENTS:
    class(energyflux_type)                :: this
    type(bounds_type)      , intent(in)    :: bounds

    !
    ! !LOCAL VARIABLES:
    integer :: m,g,l,c,p                 ! indices
    integer :: ier                       ! error status
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: year                      ! year (0, ...) for nstep
    integer :: month                     ! month (1, ..., 12) for nstep
    integer :: day                       ! day of month (1, ..., 31) for nstep
    integer :: secs                      ! seconds into current date for nstep
    logical :: end_cd                    ! temporary for is_end_curr_day() value
    integer :: begp, endp
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    dtime = get_step_size()
    nstep = get_nstep()
    call get_curr_date (year, month, day, secs)

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate and extract BTRANAV - hourly average btran
    ! Used to compute minimum of hourly averaged btran
    ! over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('BTRANAV', this%btran_patch, nstep)
    call extract_accum_field ('BTRANAV', rbufslp, nstep)
    end_cd = is_end_curr_day()
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          this%btran_min_inst_patch(p) = min(rbufslp(p), this%btran_min_inst_patch(p))
       endif
       if (end_cd) then
          this%btran_min_patch(p) = this%btran_min_inst_patch(p)
          this%btran_min_inst_patch(p) =  spval
       else if (secs == dtime) then
          this%btran_min_patch(p) = spval
       endif
    end do

    deallocate(rbufslp)

  end subroutine UpdateAccVars

end module EnergyFluxType
