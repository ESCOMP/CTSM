module atm2lndMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd forcing
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varpar     , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon     , only : rair, grav, cpair, hfus, tfrz, denh2o, spval
  use clm_varcon     , only : wv_to_dair_weight_ratio
  use clm_varctl     , only : iulog, use_c13, use_cn, use_lch4, iulog
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use atm2lndType    , only : atm2lnd_type
  use TopoMod        , only : topo_type
  use filterColMod   , only : filter_col_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col
  use landunit_varcon, only : istice_mec
  use WaterType      , only : water_type
  use Wateratm2lndBulkType, only : wateratm2lndbulk_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: set_atm2lnd_water_tracers          ! Set tracer values for the atm2lnd water quantities
  public :: downscale_forcings                 ! Downscale atm forcing fields from gridcell to column

  ! The following routine is public for the sake of unit testing; it should not be
  ! called by production code outside this module
  public :: partition_precip             ! Partition precipitation into rain/snow
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: rhos                             ! calculate atmospheric density
  private :: repartition_rain_snow_one_col    ! Re-partition precipitation for a single column
  private :: sens_heat_from_precip_conversion ! Compute sensible heat flux needed to compensate for rain-snow conversion
  private :: downscale_longwave               ! Downscale longwave radiation from gridcell to column
  private :: build_normalization              ! Compute normalization factors so that downscaled fields are conservative
  private :: check_downscale_consistency      ! Check consistency of downscaling

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine set_atm2lnd_water_tracers(bounds, num_allc, filter_allc, water_inst)
    !
    ! !DESCRIPTION:
    ! Set tracer values for the atm2lnd water quantities
    !
    ! Should be called after downscale_forcings
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_allc       ! number of column points in filter_allc
    integer           , intent(in) :: filter_allc(:) ! column filter for all points
    type(water_type)  , intent(in) :: water_inst
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'set_atm2lnd_water_tracers'
    !-----------------------------------------------------------------------

    do i = water_inst%tracers_beg, water_inst%tracers_end
       associate( &
            wateratm2lnd_inst => water_inst%bulk_and_tracers(i)%wateratm2lnd_inst)

       if (.not. wateratm2lnd_inst%IsCommunicatedWithCoupler()) then
          call wateratm2lnd_inst%SetNondownscaledTracers( &
               bounds, water_inst%wateratm2lndbulk_inst)
       end if

       call wateratm2lnd_inst%SetDownscaledTracers( &
            bounds, num_allc, filter_allc, water_inst%wateratm2lndbulk_inst)

       end associate
    end do

  end subroutine set_atm2lnd_water_tracers

  !-----------------------------------------------------------------------
  subroutine downscale_forcings(bounds, &
       topo_inst, atm2lnd_inst, wateratm2lndbulk_inst, eflx_sh_precip_conversion)
    !
    ! !DESCRIPTION:
    ! Downscale atmospheric forcing fields from gridcell to column.
    !
    ! Downscaling is done based on the difference between each CLM column's elevation and
    ! the atmosphere's surface elevation (which is the elevation at which the atmospheric
    ! forcings are valid).
    !
    ! Note that the downscaling procedure can result in changes in grid cell mean values
    ! compared to what was provided by the atmosphere. We conserve fluxes of mass and
    ! energy, but allow states such as temperature to differ.
    !
    ! For most variables, downscaling is done over columns defined by
    ! topo_inst%DownscaleFilterc. But we also do direct copies of gridcell-level forcings
    ! into column-level forcings over all other active columns. In addition, precipitation
    ! (rain vs. snow partitioning) is adjusted everywhere.
    !
    ! !USES:
    use clm_varcon      , only : rair, cpair, grav
    use QsatMod         , only : Qsat
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds  
    class(topo_type)   , intent(in)    :: topo_inst
    type(atm2lnd_type) , intent(inout) :: atm2lnd_inst
    type(wateratm2lndbulk_type) , intent(inout) :: wateratm2lndbulk_inst
    real(r8)           , intent(out)   :: eflx_sh_precip_conversion(bounds%begc:) ! sensible heat flux from precipitation conversion (W/m**2) [+ to atm]
    !
    ! !LOCAL VARIABLES:
    integer :: g, l, c, fc         ! indices
    integer :: clo, cc
    type(filter_col_type) :: downscale_filter_c

    ! temporaries for topo downscaling
    real(r8) :: hsurf_g,hsurf_c
    real(r8) :: Hbot, zbot
    real(r8) :: tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g, rhos_g
    real(r8) :: tbot_c, pbot_c, thbot_c, qbot_c, qs_c, es_c, rhos_c
    real(r8) :: rhos_c_estimate, rhos_g_estimate
    real(r8) :: dum1,   dum2

    character(len=*), parameter :: subname = 'downscale_forcings'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(eflx_sh_precip_conversion) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate(&
         ! Parameters:
         lapse_rate => atm2lnd_inst%params%lapse_rate              , & ! Input:  [real(r8)] Surface temperature lapse rate (K m-1)

         ! Gridcell-level metadata:
         forc_topo_g  => atm2lnd_inst%forc_topo_grc                , & ! Input:  [real(r8) (:)]  atmospheric surface height (m)

         ! Column-level metadata:
         topo_c       => topo_inst%topo_col                        , & ! Input:  [real(r8) (:)] column surface height (m)

         ! Gridcell-level non-downscaled fields:
         forc_t_g     => atm2lnd_inst%forc_t_not_downscaled_grc    , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_th_g    => atm2lnd_inst%forc_th_not_downscaled_grc   , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
         forc_q_g     => wateratm2lndbulk_inst%forc_q_not_downscaled_grc    , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
         forc_pbot_g  => atm2lnd_inst%forc_pbot_not_downscaled_grc , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
         forc_rho_g   => atm2lnd_inst%forc_rho_not_downscaled_grc  , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
         
         ! Column-level downscaled fields:
         forc_t_c     => atm2lnd_inst%forc_t_downscaled_col        , & ! Output: [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_th_c    => atm2lnd_inst%forc_th_downscaled_col       , & ! Output: [real(r8) (:)]  atmospheric potential temperature (Kelvin)
         forc_q_c     => wateratm2lndbulk_inst%forc_q_downscaled_col        , & ! Output: [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
         forc_pbot_c  => atm2lnd_inst%forc_pbot_downscaled_col     , & ! Output: [real(r8) (:)]  atmospheric pressure (Pa)               
         forc_rho_c   => atm2lnd_inst%forc_rho_downscaled_col        & ! Output: [real(r8) (:)]  atmospheric density (kg/m**3)           
         )
      
      ! Initialize column forcing (needs to be done for ALL active columns)
      do c = bounds%begc,bounds%endc
         if (col%active(c)) then
            g = col%gridcell(c)

            forc_t_c(c)     = forc_t_g(g)
            forc_th_c(c)    = forc_th_g(g)
            forc_q_c(c)     = forc_q_g(g)
            forc_pbot_c(c)  = forc_pbot_g(g)
            forc_rho_c(c)   = forc_rho_g(g)
         end if
      end do

      downscale_filter_c = topo_inst%DownscaleFilterc(bounds)

      ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
      ! For glacier_mec columns the downscaling is based on surface elevation.
      ! For other columns the downscaling is a simple copy (above).
      do fc = 1, downscale_filter_c%num
         c = downscale_filter_c%indices(fc)
         l = col%landunit(c)
         g = col%gridcell(c)

         ! This is a simple downscaling procedure 
         ! Note that forc_hgt, forc_u, and forc_v are not downscaled.

         hsurf_g = forc_topo_g(g)                        ! gridcell sfc elevation
         hsurf_c = topo_c(c)                             ! column sfc elevation
         tbot_g  = forc_t_g(g)                           ! atm sfc temp
         thbot_g = forc_th_g(g)                          ! atm sfc pot temp
         qbot_g  = forc_q_g(g)                           ! atm sfc spec humid
         pbot_g  = forc_pbot_g(g)                        ! atm sfc pressure
         rhos_g  = forc_rho_g(g)                         ! atm density
         zbot    = atm2lnd_inst%forc_hgt_grc(g)          ! atm ref height
         tbot_c  = tbot_g-lapse_rate*(hsurf_c-hsurf_g)   ! sfc temp for column
         Hbot    = rair*0.5_r8*(tbot_g+tbot_c)/grav      ! scale ht at avg temp
         pbot_c  = pbot_g*exp(-(hsurf_c-hsurf_g)/Hbot)   ! column sfc press

         ! Derivation of potential temperature calculation:
         ! 
         ! The textbook definition would be:
         ! thbot_c = tbot_c * (p0/pbot_c)^(rair/cpair)
         ! 
         ! Note that pressure is related to scale height as:
         ! pbot_c = p0 * exp(-zbot/H)
         !
         ! Using Hbot in place of H, we get:
         ! pbot_c = p0 * exp(-zbot/Hbot)
         !
         ! Plugging this in to the textbook definition, then manipulating, we get:
         ! thbot_c = tbot_c * (p0/(p0*exp(-zbot/Hbot)))^(rair/cpair)
         !         = tbot_c * (1/exp(-zbot/Hbot))^(rair/cpair)
         !         = tbot_c * (exp(zbot/Hbot))^(rair/cpair)
         !         = tbot_c * exp((zbot/Hbot) * (rair/cpair))
         !
         ! But we want everything expressed in delta form, resulting in:
         ! thbot_c = thbot_g + (tbot_c - tbot_g)*exp((zbot/Hbot)*(rair/cpair))

         thbot_c= thbot_g + (tbot_c - tbot_g)*exp((zbot/Hbot)*(rair/cpair))  ! pot temp calc

         call Qsat(tbot_g,pbot_g,es_g,dum1,qs_g,dum2)
         call Qsat(tbot_c,pbot_c,es_c,dum1,qs_c,dum2)

         qbot_c = qbot_g*(qs_c/qs_g)

         ! For forc_rho_c: We could simply set:
         !
         !    rhos_c = rhos(pbot_c, egcm_c, tbot_c)
         !
         ! However, we want forc_rho_c to be identical to forc_rho_g when topo_c equals
         ! forc_topo_g. So we compute our own version of forc_rho_g using the rhos
         ! function, and then multiply forc_rho_g by the ratio of (computed column-level
         ! rho) to (computed gridcell-level rho).
         rhos_c_estimate = rhos(qbot=qbot_c, pbot=pbot_c, tbot=tbot_c)
         rhos_g_estimate = rhos(qbot=qbot_g, pbot=pbot_g, tbot=tbot_g)
         rhos_c = rhos_g * (rhos_c_estimate / rhos_g_estimate)

         forc_t_c(c)    = tbot_c
         forc_th_c(c)   = thbot_c
         forc_q_c(c)    = qbot_c
         forc_pbot_c(c) = pbot_c
         forc_rho_c(c)  = rhos_c

      end do

      call partition_precip(bounds, atm2lnd_inst, wateratm2lndbulk_inst, &
           eflx_sh_precip_conversion(bounds%begc:bounds%endc))

      call downscale_longwave(bounds, downscale_filter_c, topo_inst, atm2lnd_inst)

      call check_downscale_consistency(bounds, atm2lnd_inst, wateratm2lndbulk_inst)

    end associate

  end subroutine downscale_forcings

  !-----------------------------------------------------------------------
  pure function rhos(qbot, pbot, tbot)
    !
    ! !DESCRIPTION:
    ! Compute atmospheric density (kg/m**3)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: rhos  ! function result: atmospheric density (kg/m**3)
    real(r8), intent(in) :: qbot  ! atmospheric specific humidity (kg/kg)
    real(r8), intent(in) :: pbot  ! atmospheric pressure (Pa)
    real(r8), intent(in) :: tbot  ! atmospheric temperature (K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: egcm

    character(len=*), parameter :: subname = 'rhos'
    !-----------------------------------------------------------------------

    egcm = qbot*pbot / &
         (wv_to_dair_weight_ratio + (1._r8 - wv_to_dair_weight_ratio)*qbot)
    rhos = (pbot - (1._r8 - wv_to_dair_weight_ratio)*egcm) / (rair*tbot)
    
  end function rhos

  !-----------------------------------------------------------------------
  subroutine partition_precip(bounds, atm2lnd_inst, wateratm2lndbulk_inst, eflx_sh_precip_conversion)
    !
    ! !DESCRIPTION:
    ! Partition precipitation into rain/snow based on temperature.
    !
    ! Note that, unlike the other downscalings done here, this is currently applied over
    ! all points - not just those within the downscale filter.
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds  
    type(atm2lnd_type) , intent(inout) :: atm2lnd_inst
    type(wateratm2lndbulk_type) , intent(inout) :: wateratm2lndbulk_inst
    real(r8), intent(inout) :: eflx_sh_precip_conversion(bounds%begc:) ! sensible heat flux from precipitation conversion (W/m**2) [+ to atm]
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,g      ! indices
    real(r8) :: rain_orig  ! rain before conversion
    real(r8) :: snow_orig  ! snow before conversion
    real(r8) :: all_snow_t ! temperature at which all precip falls as snow (K)
    real(r8) :: frac_rain_slope ! slope of the frac_rain vs. temperature relationship

    character(len=*), parameter :: subname = 'partition_precip'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(eflx_sh_precip_conversion) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate(&
         ! Gridcell-level non-downscaled fields:
         forc_rain_g  => wateratm2lndbulk_inst%forc_rain_not_downscaled_grc , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
         forc_snow_g  => wateratm2lndbulk_inst%forc_snow_not_downscaled_grc , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
         
         ! Column-level downscaled fields:
         forc_t_c                  => atm2lnd_inst%forc_t_downscaled_col                , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_rain_c               => wateratm2lndbulk_inst%forc_rain_downscaled_col    , & ! Output: [real(r8) (:)]  rain rate [mm/s]
         forc_snow_c               => wateratm2lndbulk_inst%forc_snow_downscaled_col    , & ! Output: [real(r8) (:)]  snow rate [mm/s]
         rain_to_snow_conversion_c => wateratm2lndbulk_inst%rain_to_snow_conversion_col , & ! Output: [real(r8) (:)]  amount of rain converted to snow via precipitation repartition [mm/s]
         snow_to_rain_conversion_c => wateratm2lndbulk_inst%snow_to_rain_conversion_col   & ! Output: [real(r8) (:)]  amount of snow converted to rain via precipitation repartition [mm/s]
         )

    ! Initialize column forcing
    do c = bounds%begc,bounds%endc
       if (col%active(c)) then
          g = col%gridcell(c)
          forc_rain_c(c)  = forc_rain_g(g)
          forc_snow_c(c)  = forc_snow_g(g)
          rain_to_snow_conversion_c(c) = 0._r8
          snow_to_rain_conversion_c(c) = 0._r8
          eflx_sh_precip_conversion(c) = 0._r8
       end if
    end do

    ! Optionally, convert rain to snow or vice versa based on forc_t_c
    if (atm2lnd_inst%params%repartition_rain_snow) then
       do c = bounds%begc, bounds%endc
          if (col%active(c)) then
             l = col%landunit(c)
             rain_orig = forc_rain_c(c)
             snow_orig = forc_snow_c(c)
             if (lun%itype(l) == istice_mec) then
                all_snow_t = atm2lnd_inst%params%precip_repartition_glc_all_snow_t
                frac_rain_slope = atm2lnd_inst%params%precip_repartition_glc_frac_rain_slope
             else
                all_snow_t = atm2lnd_inst%params%precip_repartition_nonglc_all_snow_t
                frac_rain_slope = atm2lnd_inst%params%precip_repartition_nonglc_frac_rain_slope
             end if
             call repartition_rain_snow_one_col(&
                  temperature = forc_t_c(c), &
                  all_snow_t = all_snow_t, &
                  frac_rain_slope = frac_rain_slope, &
                  rain = forc_rain_c(c), &
                  snow = forc_snow_c(c))
             if (forc_rain_c(c) > rain_orig) then
                snow_to_rain_conversion_c(c) = forc_rain_c(c) - rain_orig
             end if
             if (forc_snow_c(c) > snow_orig) then
                rain_to_snow_conversion_c(c) = forc_snow_c(c) - snow_orig
             end if
             call sens_heat_from_precip_conversion(&
                  rain_to_snow = rain_to_snow_conversion_c(c), &
                  snow_to_rain = snow_to_rain_conversion_c(c), &
                  sens_heat_flux = eflx_sh_precip_conversion(c))
          end if
       end do
    end if

    end associate

  end subroutine partition_precip

  !-----------------------------------------------------------------------
  subroutine repartition_rain_snow_one_col(temperature, all_snow_t, frac_rain_slope, &
       rain, snow)
    !
    ! !DESCRIPTION:
    ! Re-partition precipitation into rain/snow for a single column.
    !
    ! Rain and snow variables should be set initially, and are updated here
    !
    ! !ARGUMENTS:
    real(r8) , intent(in)    :: temperature ! near-surface temperature (K)
    real(r8) , intent(in)    :: all_snow_t  ! temperature at which precip falls entirely as snow (K)
    real(r8) , intent(in)    :: frac_rain_slope ! slope of the frac_rain vs. T relationship
    real(r8) , intent(inout) :: rain        ! atm rain rate [mm/s]
    real(r8) , intent(inout) :: snow        ! atm snow rate [(mm water equivalent)/s]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: frac_rain     ! fraction of precipitation that should become rain
    real(r8) :: total_precip

    character(len=*), parameter :: subname = 'repartition_rain_snow_one_col'
    !-----------------------------------------------------------------------

    frac_rain = (temperature - all_snow_t) * frac_rain_slope

    ! bound in [0,1]
    frac_rain = min(1.0_r8,max(0.0_r8,frac_rain))

    total_precip = rain + snow
    rain = total_precip * frac_rain
    snow = total_precip - rain

  end subroutine repartition_rain_snow_one_col

  !-----------------------------------------------------------------------
  subroutine sens_heat_from_precip_conversion(rain_to_snow, snow_to_rain, sens_heat_flux)
    !
    ! !DESCRIPTION:
    ! Given conversion fluxes from rain to snow and snow to rain, compute the sensible
    ! heat flux needed to compensate for the rain-snow conversion.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: rain_to_snow    ! amount of rain converted to snow [mm/s]
    real(r8), intent(in)  :: snow_to_rain    ! amount of snow converted to rain [mm/s]
    real(r8), intent(out) :: sens_heat_flux  ! [W/m^2]
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: mm_to_m = 1.e-3_r8  ! multiply by this to convert from mm to m
    real(r8), parameter :: tol = 1.e-13_r8     ! relative tolerance for error checks

    character(len=*), parameter :: subname = 'sens_heat_from_precip_conversion'
    !-----------------------------------------------------------------------

    ! rain to snow releases energy, so results in a positive heat flux to atm
    sens_heat_flux = (rain_to_snow - snow_to_rain) * mm_to_m * denh2o * hfus

  end subroutine sens_heat_from_precip_conversion


  !-----------------------------------------------------------------------
  subroutine downscale_longwave(bounds, downscale_filter_c, &
       topo_inst, atm2lnd_inst)
    !
    ! !DESCRIPTION:
    ! Downscale longwave radiation from gridcell to column
    ! Must be done AFTER temperature downscaling
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    type(filter_col_type) , intent(in)    :: downscale_filter_c
    class(topo_type)      , intent(in)    :: topo_inst
    type(atm2lnd_type)    , intent(inout) :: atm2lnd_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,g,fc     ! indices
    real(r8) :: hsurf_c      ! column-level elevation (m)
    real(r8) :: hsurf_g      ! gridcell-level elevation (m)

    real(r8), dimension(bounds%begg : bounds%endg) :: sum_lwrad_g    ! weighted sum of column-level lwrad
    real(r8), dimension(bounds%begg : bounds%endg) :: sum_wts_g      ! sum of weights that contribute to sum_lwrad_g
    real(r8), dimension(bounds%begg : bounds%endg) :: lwrad_norm_g   ! normalization factors
    real(r8), dimension(bounds%begg : bounds%endg) :: newsum_lwrad_g ! weighted sum of column-level lwrad after normalization

    character(len=*), parameter :: subname = 'downscale_longwave'
    !-----------------------------------------------------------------------

    associate(&
         ! Parameters:
         lapse_rate_longwave => atm2lnd_inst%params%lapse_rate_longwave  , & ! Input:  [real(r8)] longwave radiation lapse rate (W m-2 m-1)
         longwave_downscaling_limit => atm2lnd_inst%params%longwave_downscaling_limit, & ! Input:  [real(r8)] Relative limit for how much longwave downscaling can be done (unitless)

         ! Gridcell-level metadata:
         forc_topo_g  => atm2lnd_inst%forc_topo_grc                , & ! Input:  [real(r8) (:)]  atmospheric surface height (m)

         ! Column-level metadata:
         topo_c       => topo_inst%topo_col                        , & ! Input:  [real(r8) (:)] column surface height (m)

         ! Gridcell-level fields:
         forc_lwrad_g => atm2lnd_inst%forc_lwrad_not_downscaled_grc, & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)
         
         ! Column-level (downscaled) fields:
         forc_lwrad_c => atm2lnd_inst%forc_lwrad_downscaled_col      & ! Output: [real(r8) (:)]  downward longwave (W/m**2)
         )
    
      ! Initialize column forcing (needs to be done for ALL active columns)
      do c = bounds%begc, bounds%endc
         if (col%active(c)) then
            g = col%gridcell(c)
            forc_lwrad_c(c) = forc_lwrad_g(g)
         end if
      end do

      ! Optionally, downscale the longwave radiation, conserving energy
      if (atm2lnd_inst%params%glcmec_downscale_longwave) then

         ! Initialize variables related to normalization
         do g = bounds%begg, bounds%endg
            sum_lwrad_g(g) = 0._r8
            sum_wts_g(g) = 0._r8
            newsum_lwrad_g(g) = 0._r8
         end do

         ! Do the downscaling
         do fc = 1, downscale_filter_c%num
            c = downscale_filter_c%indices(fc)
            l = col%landunit(c)
            g = col%gridcell(c)

            hsurf_g = forc_topo_g(g)
            hsurf_c = topo_c(c)

            ! Assume a linear decrease in downwelling longwave radiation with increasing
            ! elevation, based on Van Tricht et al. (2016, TC) Figure 6,
            ! doi:10.5194/tc-10-2379-2016
            forc_lwrad_c(c) = forc_lwrad_g(g) - lapse_rate_longwave * (hsurf_c-hsurf_g)
            ! But ensure that we don't depart too far from the atmospheric forcing value:
            ! negative values of lwrad are certainly bad, but small positive values might
            ! also be bad. We can especially run into trouble due to the normalization: a
            ! small lwrad value in one column can lead to a big normalization factor,
            ! leading to huge lwrad values in other columns.
            forc_lwrad_c(c) = min(forc_lwrad_c(c), &
                 forc_lwrad_g(g) * (1._r8 + longwave_downscaling_limit))
            forc_lwrad_c(c) = max(forc_lwrad_c(c), &
                 forc_lwrad_g(g) * (1._r8 - longwave_downscaling_limit))

            ! Keep track of the gridcell-level weighted sum for later normalization.
            !
            ! This gridcell-level weighted sum just includes points for which we do the
            ! downscaling (e.g., glc_mec points). Thus the contributing weights
            ! generally do not add to 1. So to do the normalization properly, we also
            ! need to keep track of the weights that have contributed to this sum.
            sum_lwrad_g(g) = sum_lwrad_g(g) + col%wtgcell(c)*forc_lwrad_c(c)
            sum_wts_g(g) = sum_wts_g(g) + col%wtgcell(c)
         end do


         ! Normalize forc_lwrad_c(c) to conserve energy

         call build_normalization(orig_field=forc_lwrad_g(bounds%begg:bounds%endg), &
              sum_field=sum_lwrad_g(bounds%begg:bounds%endg), &
              sum_wts=sum_wts_g(bounds%begg:bounds%endg), &
              norms=lwrad_norm_g(bounds%begg:bounds%endg))

         do fc = 1, downscale_filter_c%num
            c = downscale_filter_c%indices(fc)
            l = col%landunit(c)
            g = col%gridcell(c)

            forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
            newsum_lwrad_g(g) = newsum_lwrad_g(g) + col%wtgcell(c)*forc_lwrad_c(c)
         end do


         ! Make sure that, after normalization, the grid cell mean is conserved

         do g = bounds%begg, bounds%endg
            if (sum_wts_g(g) > 0._r8) then
               if (abs((newsum_lwrad_g(g) / sum_wts_g(g)) - forc_lwrad_g(g)) > 1.e-8_r8) then
                  write(iulog,*) 'g, newsum_lwrad_g, sum_wts_g, forc_lwrad_g: ', &
                       g, newsum_lwrad_g(g), sum_wts_g(g), forc_lwrad_g(g)
                  call endrun(msg=' ERROR: Energy conservation error downscaling longwave'//&
                       errMsg(sourcefile, __LINE__))
               end if
            end if
         end do

      end if    ! glcmec_downscale_longwave

    end associate

  end subroutine downscale_longwave

  !-----------------------------------------------------------------------
  subroutine build_normalization(orig_field, sum_field, sum_wts, norms)
    !
    ! !DESCRIPTION:
    ! Build an array of normalization factors that can be applied to a downscaled forcing
    ! field, in order to force the mean of the new field to be the same as the mean of
    ! the old field (for conservation).
    !
    ! This allows for the possibility that only a subset of columns are downscaled. Only
    ! the columns that are adjusted should be included in the weighted sum, sum_field;
    ! sum_wts gives the sum of contributing weights on the grid cell level. 

    ! For example, if a grid cell has an original forcing value of 1.0, and contains 4
    ! columns with the following weights on the gridcell, and the following values after
    ! normalization:
    !
    !       col #:    1     2     3     4
    !      weight:  0.1   0.2   0.3   0.4
    ! downscaled?:  yes   yes    no    no
    !       value:  0.9   1.1   1.0   1.0
    !
    ! Then we would have:
    ! orig_field(g) = 1.0
    ! sum_field(g) = 0.1*0.9 + 0.2*1.1 = 0.31
    ! sum_wts(g) = 0.1 + 0.2 = 0.3
    ! norms(g) = 1.0 / (0.31 / 0.3) = 0.9677
    !
    ! The field can then be normalized as:
    !              forc_lwrad_c(c) = forc_lwrad_c(c) * lwrad_norm_g(g)
    !   where lwrad_norm_g is the array of norms computed by this routine

    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: orig_field(:)  ! the original field, at the grid cell level
    real(r8), intent(in)  :: sum_field(:)   ! the new weighted sum across columns (dimensioned by grid cell)
    real(r8), intent(in)  :: sum_wts(:)     ! sum of the weights used to create sum_field (dimensioned by grid cell)
    real(r8), intent(out) :: norms(:)       ! computed normalization factors
    !-----------------------------------------------------------------------

    SHR_ASSERT((size(orig_field) == size(norms)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT((size(sum_field) == size(norms)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT((size(sum_wts) == size(norms)), errMsg(sourcefile, __LINE__))

    where (sum_wts == 0._r8)
       ! Avoid divide by zero; if sum_wts is 0, then the normalization doesn't matter,
       ! because the adjusted values won't affect the grid cell mean.
       norms = 1.0_r8

    elsewhere (sum_field == 0._r8)
       ! Avoid divide by zero. If this is because both sum_field and orig_field are 0,
       ! then the normalization doesn't matter. If sum_field == 0 while orig_field /= 0,
       ! then we have a problem: no normalization will allow us to recover the original
       ! gridcell mean. We should probably catch this and abort, but for now we're
       ! relying on error checking in the caller (checking for conservation) to catch
       ! this potential problem.
       norms = 1.0_r8

    elsewhere
       ! The standard case
       norms = orig_field / (sum_field / sum_wts)

    end where

  end subroutine build_normalization


  !-----------------------------------------------------------------------
  subroutine check_downscale_consistency(bounds, atm2lnd_inst, wateratm2lndbulk_inst)
    !
    ! !DESCRIPTION:
    ! Check consistency of downscaling
    !
    ! Note that this operates over more than just the filter used for the downscaling,
    ! because it checks some things outside that filter.
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in) :: bounds  
    type(atm2lnd_type), intent(in) :: atm2lnd_inst
    type(wateratm2lndbulk_type), intent(in) :: wateratm2lndbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: g, l, c    ! indices
    character(len=*), parameter :: subname = 'check_downscale_consistency'
    !-----------------------------------------------------------------------

    associate(&
         ! Gridcell-level fields:
         forc_t_g     => atm2lnd_inst%forc_t_not_downscaled_grc     , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_th_g    => atm2lnd_inst%forc_th_not_downscaled_grc    , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
         forc_q_g     => wateratm2lndbulk_inst%forc_q_not_downscaled_grc     , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
         forc_pbot_g  => atm2lnd_inst%forc_pbot_not_downscaled_grc  , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
         forc_rho_g   => atm2lnd_inst%forc_rho_not_downscaled_grc   , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
         forc_rain_g  => wateratm2lndbulk_inst%forc_rain_not_downscaled_grc  , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
         forc_snow_g  => wateratm2lndbulk_inst%forc_snow_not_downscaled_grc  , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
         forc_lwrad_g => atm2lnd_inst%forc_lwrad_not_downscaled_grc , & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)
         
         ! Column-level (downscaled) fields:
         forc_t_c     => atm2lnd_inst%forc_t_downscaled_col         , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
         forc_th_c    => atm2lnd_inst%forc_th_downscaled_col        , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)
         forc_q_c     => wateratm2lndbulk_inst%forc_q_downscaled_col         , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)   
         forc_pbot_c  => atm2lnd_inst%forc_pbot_downscaled_col      , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)               
         forc_rho_c   => atm2lnd_inst%forc_rho_downscaled_col       , & ! Input:  [real(r8) (:)]  atmospheric density (kg/m**3)           
         forc_rain_c  => wateratm2lndbulk_inst%forc_rain_downscaled_col      , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
         forc_snow_c  => wateratm2lndbulk_inst%forc_snow_downscaled_col      , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
         forc_lwrad_c => atm2lnd_inst%forc_lwrad_downscaled_col       & ! Input:  [real(r8) (:)]  downward longwave (W/m**2)
         )

    ! BUG(wjs, 2016-11-15, bugz 2377)
    !
    ! Make sure that, for urban points, the column-level forcing fields are identical to
    ! the gridcell-level forcing fields. This is needed because the urban-specific code
    ! sometimes uses the gridcell-level forcing fields (and it would take a large
    ! refactor to change this to use column-level fields).
    !
    ! However, do NOT check rain & snow: these ARE downscaled for urban points (as for
    ! all other points), and the urban code does not refer to the gridcell-level versions
    ! of these fields.
    
    do c = bounds%begc, bounds%endc
       if (col%active(c)) then
          l = col%landunit(c)
          g = col%gridcell(c)

          if (lun%urbpoi(l)) then
             if (forc_t_c(c)     /= forc_t_g(g)    .or. &
                  forc_th_c(c)    /= forc_th_g(g)   .or. &
                  forc_q_c(c)     /= forc_q_g(g)    .or. &
                  forc_pbot_c(c)  /= forc_pbot_g(g) .or. &
                  forc_rho_c(c)   /= forc_rho_g(g)  .or. &
                  forc_lwrad_c(c) /= forc_lwrad_g(g)) then
                write(iulog,*) subname//' ERROR: column-level forcing differs from gridcell-level forcing for urban point'
                write(iulog,*) 'c, g = ', c, g
                write(iulog,*) 'forc_t_c, forc_t_g = ', forc_t_c(c), forc_t_g(g)
                write(iulog,*) 'forc_th_c, forc_th_g = ', forc_th_c(c), forc_th_g(g)
                write(iulog,*) 'forc_q_c, forc_q_g = ', forc_q_c(c), forc_q_g(g)
                write(iulog,*) 'forc_pbot_c, forc_pbot_g = ', forc_pbot_c(c), forc_pbot_g(g)
                write(iulog,*) 'forc_rho_c, forc_rho_g = ', forc_rho_c(c), forc_rho_g(g)
                write(iulog,*) 'forc_lwrad_c, forc_lwrad_g = ', forc_lwrad_c(c), forc_lwrad_g(g)
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if  ! inequal
          end if  ! urbpoi
       end if  ! active
    end do

    end associate

  end subroutine check_downscale_consistency

end module atm2lndMod
