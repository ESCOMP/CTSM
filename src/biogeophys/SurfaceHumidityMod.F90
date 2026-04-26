module SurfaceHumidityMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate surface humidities, as well as a few intermediate variables that are needed
  ! in the humidity calculations

  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use decompMod               , only : bounds_type
  use abortutils              , only : endrun
  use clm_varcon              , only : denh2o, denice, roverg, tfrz, spval
  use column_varcon           , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon           , only : icol_road_imperv, icol_road_perv
  use landunit_varcon         , only : istice, istwet, istsoil, istcrop
  use clm_varpar              , only : nlevgrnd
  ! [PORTED by Hui Tang: use_nvp flag for NVP ground evaporation blending]
  use clm_varctl              , only : use_nvp
  use NVPLayerDynamicsMod     , only : NVPWaterRetentionCurve
    ! [PORTED by Hui Tang: runtime-tunable NVP physics parameters]
  use NVPParamsMod            , only : n_van_nvp, alpha_van_nvp, watsat_nvp, watres_nvp
  use atm2lndType             , only : atm2lnd_type
  use SoilStateType           , only : soilstate_type
  use TemperatureType         , only : temperature_type
  use Wateratm2lndBulkType    , only : wateratm2lndbulk_type
  use WaterStateBulkType      , only : waterstatebulk_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use LandunitType            , only : lun                
  use ColumnType              , only : col                
  use QSatMod                 , only : QSat
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CalculateSurfaceHumidity  

  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine CalculateSurfaceHumidity(bounds, &
       num_nolakec, filter_nolakec, &
       atm2lnd_inst, temperature_inst, &
       waterstatebulk_inst, wateratm2lndbulk_inst, &
       soilstate_inst, waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Calculate surface humidities, as well as a few intermediate variables that are
    ! needed in the humidity calculations
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds    
    integer                        , intent(in)    :: num_nolakec       ! number of column non-lake points in column filter
    integer                        , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    type(atm2lnd_type)             , intent(in)    :: atm2lnd_inst
    type(temperature_type)         , intent(in)    :: temperature_inst
    type(waterstatebulk_type)      , intent(in)    :: waterstatebulk_inst
    type(wateratm2lndbulk_type)    , intent(in)    :: wateratm2lndbulk_inst
    type(soilstate_type)           , intent(inout) :: soilstate_inst
    type(waterdiagnosticbulk_type) , intent(inout) :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: g,l,c,p      ! indices
    integer  :: j            ! soil/snow level index
    integer  :: fp           ! lake filter patch index
    integer  :: fc           ! lake filter column index
    real(r8) :: qred         ! soil surface relative humidity
    real(r8) :: qsatg        ! saturated humidity [kg/kg]
    real(r8) :: qsatgdT      ! d(qsatg)/dT
    real(r8) :: qsatgdT_snow ! d(qsatg)/dT, for snow
    real(r8) :: qsatgdT_soil ! d(qsatg)/dT, for soil
    real(r8) :: qsatgdT_h2osfc ! d(qsatg)/dT, for h2osfc
    ! [PORTED by Hui Tang: NVP ground humidity variables]
    real(r8) :: qsatgdT_nvp  ! d(qsatg)/dT, for NVP surface
    real(r8) :: frac_nvp_eff ! effective NVP fraction for ground humidity blend
    real(r8) :: hr_nvp       ! alpha NVP
    real(r8) :: psit_nvp     ! negative potential of NVP
    real(r8) :: fac          ! soil wetness of surface layer
    real(r8) :: psit         ! negative potential of soil
    real(r8) :: hr           ! alpha soil
    real(r8) :: hr_road_perv ! alpha soil for urban pervious road
    real(r8) :: wx           ! partial volume of ice and water of surface layer
    real(r8) :: fac_fc       ! soil wetness of surface layer relative to field capacity
    real(r8) :: eff_porosity ! effective porosity in layer
    real(r8) :: vol_ice      ! partial volume of ice lens in layer
    real(r8) :: vol_liq      ! partial volume of liquid water in layer
    !------------------------------------------------------------------------------

    associate( & 
         snl              =>    col%snl                                     , & ! Input:  [integer  (:)   ] number of snow layers
         dz               =>    col%dz                                      , & ! Input:  [real(r8) (:,:) ] layer depth (m)

         forc_pbot        =>    atm2lnd_inst%forc_pbot_downscaled_col       , & ! Input:  [real(r8) (:)   ] atmospheric pressure (Pa)
         forc_q           =>    wateratm2lndbulk_inst%forc_q_downscaled_col , & ! Input:  [real(r8) (:)   ] atmospheric specific humidity (kg/kg)


         frac_h2osfc      =>    waterdiagnosticbulk_inst%frac_h2osfc_col    , & ! Input:  [real(r8) (:)   ] fraction of ground covered by surface water (0 to 1)
         frac_sno_eff     =>    waterdiagnosticbulk_inst%frac_sno_eff_col   , & ! Input:  [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)
         h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col          , & ! Input:  [real(r8) (:,:) ] ice lens (kg/m2)
         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col          , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)
         qg_snow          =>    waterdiagnosticbulk_inst%qg_snow_col        , & ! Output: [real(r8) (:)   ] specific humidity at snow surface [kg/kg]
         qg_soil          =>    waterdiagnosticbulk_inst%qg_soil_col        , & ! Output: [real(r8) (:)   ] specific humidity at soil surface [kg/kg]
         qg               =>    waterdiagnosticbulk_inst%qg_col             , & ! Output: [real(r8) (:)   ] ground specific humidity [kg/kg]
         qg_h2osfc        =>    waterdiagnosticbulk_inst%qg_h2osfc_col      , & ! Output: [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]
         dqgdT            =>    waterdiagnosticbulk_inst%dqgdT_col          , & ! Output: [real(r8) (:)   ] d(qg)/dT
         ! [PORTED by Hui Tang: NVP ground humidity fields for ground evap blending]
         qg_nvp           =>    waterdiagnosticbulk_inst%qg_nvp_col         , & ! Output: [real(r8) (:)   ] NVP surface specific humidity [kg/kg]
         smpmin           =>    soilstate_inst%smpmin_col                   , & ! Input:  [real(r8) (:)   ] restriction for min of soil potential (mm)
         sucsat           =>    soilstate_inst%sucsat_col                   , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)
         watsat           =>    soilstate_inst%watsat_col                   , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)
         watdry           =>    soilstate_inst%watdry_col                   , & ! Input:  [real(r8) (:,:) ] volumetric soil moisture corresponding to no restriction on ET from urban pervious surface
         watopt           =>    soilstate_inst%watopt_col                   , & ! Input:  [real(r8) (:,:) ] volumetric soil moisture corresponding to no restriction on ET from urban pervious surface
         bsw              =>    soilstate_inst%bsw_col                      , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
         rootfr_road_perv =>    soilstate_inst%rootfr_road_perv_col         , & ! Input:  [real(r8) (:,:) ] fraction of roots in each soil layer for urban pervious road
         rootr_road_perv  =>    soilstate_inst%rootr_road_perv_col          , & ! Output: [real(r8) (:,:) ] effective fraction of roots in each soil layer for urban pervious road
         soilalpha        =>    soilstate_inst%soilalpha_col                , & ! Output: [real(r8) (:)   ] factor that reduces ground saturated specific humidity (-)
         soilalpha_u      =>    soilstate_inst%soilalpha_u_col              , & ! Output: [real(r8) (:)   ] Urban factor that reduces ground saturated specific humidity (-)

         t_h2osfc         =>    temperature_inst%t_h2osfc_col               , & ! Input:  [real(r8) (:)   ] surface water temperature
         t_soisno         =>    temperature_inst%t_soisno_col               , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)
         t_grnd           =>    temperature_inst%t_grnd_col                 , & ! Input:  [real(r8) (:)   ] ground temperature (Kelvin)
         ! [PORTED by Hui Tang: NVP layer temperature for NVP surface humidity]
         t_nvp_col        =>    temperature_inst%t_nvp_col                    & ! Input:  [real(r8) (:)   ] NVP (moss/lichen) temperature (Kelvin)
         )

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)

         if (col%itype(c) == icol_road_perv) then
            hr_road_perv = 0._r8
         end if

         ! Saturated vapor pressure, specific humidity and their derivatives
         ! at ground surface
         qred = 1._r8
         if (lun%itype(l)/=istwet .AND. lun%itype(l)/=istice) then

            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)
               fac  = min(1._r8, wx/watsat(c,1))
               fac  = max( fac, 0.01_r8 )
               psit = -sucsat(c,1) * fac ** (-bsw(c,1))
               psit = max(smpmin(c), psit)
               ! modify qred to account for h2osfc
               hr   = exp(psit/roverg/t_soisno(c,1))

               ! [PORTED by Hui Tang: NVP effective fraction for ground humidity blend]
               ! NVP occupies area not covered by snow or surface water
               if (use_nvp) then
                  ! Compute NVP surface humidity as a function of NVP water retention curve
                  ! --- NVP volumetric water content (clamped to valid range) ---
                  if (dz(c,0) > 0._r8) then
                     if (t_soisno(c,0) >= tfrz) then
                        ! For unfrozen soil
                        vol_ice = min(watsat_nvp, h2osoi_ice(c,0)/(dz(c,0)*denice))
                        eff_porosity = watsat_nvp-vol_ice
                        vol_liq = min(eff_porosity, h2osoi_liq(c,0)/(dz(c,0)*denh2o))
                     else
                        ! For frozen soil, assume NVP water content is at residual (unavailable for evaporation)
                        vol_liq = watres_nvp            
                     end if
                     call NVPWaterRetentionCurve(vol_liq, eff_porosity, n_van_nvp, alpha_van_nvp, &
                              watsat_nvp, watres_nvp, psit_nvp)
                     hr_nvp = exp(psit_nvp/roverg/t_nvp_col(c))
                  else
                     ! If dz(c,0) is not positive, set hr_nvp to 0
                     hr_nvp = 0._r8
                  end if
               else
                  hr_nvp = 0._r8
               end if       

               ! [PORTED by Hui Tang: NVP effective fraction for ground humidity blend]
               ! NVP occupies area not covered by snow or surface water
               if (use_nvp) then
                  frac_nvp_eff = min(col%frac_nvp(c), max(0._r8, 1._r8 - frac_sno_eff(c) - frac_h2osfc(c)))
                  qred = (1._r8 - frac_sno_eff(c) - frac_h2osfc(c) - frac_nvp_eff)*hr &
                       + frac_sno_eff(c) + frac_h2osfc(c) + frac_nvp_eff*hr_nvp
               else
                  frac_nvp_eff = 0._r8
                  qred = (1._r8 - frac_sno_eff(c) - frac_h2osfc(c))*hr &
                       + frac_sno_eff(c) + frac_h2osfc(c)
               end if
               soilalpha(c) = qred

            else if (col%itype(c) == icol_road_perv) then
               ! Pervious road depends on water in total soil column
               do j = 1, nlevgrnd
                  if (t_soisno(c,j) >= tfrz) then
                     vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
                     eff_porosity = watsat(c,j)-vol_ice
                     vol_liq = min(eff_porosity, h2osoi_liq(c,j)/(dz(c,j)*denh2o))
                     fac = min( max(vol_liq-watdry(c,j),0._r8) / (watopt(c,j)-watdry(c,j)), 1._r8 )
                  else
                     fac = 0._r8
                  end if
                  rootr_road_perv(c,j) = rootfr_road_perv(c,j)*fac
                  hr_road_perv = hr_road_perv + rootr_road_perv(c,j)
               end do
               ! Allows for sublimation of snow or dew on snow
               qred = (1.-frac_sno_eff(c))*hr_road_perv + frac_sno_eff(c)

               ! Normalize root resistances to get layer contribution to total ET
               if (hr_road_perv > 0._r8) then
                  do j = 1, nlevgrnd
                     rootr_road_perv(c,j) = rootr_road_perv(c,j)/hr_road_perv
                  end do
               end if
               soilalpha_u(c) = qred

            else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
               qred = 0._r8
               soilalpha_u(c) = spval

            else if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
               qred = 1._r8
               soilalpha_u(c) = spval
            end if

         else
            soilalpha(c) = spval

         end if

         ! compute humidities individually for snow, soil, h2osfc for vegetated landunits
         if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

            call QSat(t_soisno(c,1), forc_pbot(c), qsatg, &
                 qsdT = qsatgdT_soil)
            if (qsatg > forc_q(c) .and. forc_q(c) > hr*qsatg) then
               qsatg = forc_q(c)
               qsatgdT_soil = 0._r8
            end if
            qg_soil(c) = hr*qsatg

            if (snl(c) < 0) then
               call QSat(t_soisno(c,snl(c)+1), forc_pbot(c), qsatg, &
                    qsdT = qsatgdT_snow)
               qg_snow(c) = qsatg
               dqgdT(c) = frac_sno_eff(c)*qsatgdT_snow + &
                    (1._r8 - frac_sno_eff(c) - frac_h2osfc(c))*hr*qsatgdT_soil
            else
               ! To be consistent with hs_top values in SoilTemp, set qg_snow to qg_soil
               ! for snl = 0 case. This ensures hs_top_snow will equal hs_top_soil.
               qg_snow(c) = qg_soil(c)
               dqgdT(c) = (1._r8 - frac_h2osfc(c))*hr*qsatgdT_soil
            endif

            if (frac_h2osfc(c) > 0._r8) then
               call QSat(t_h2osfc(c), forc_pbot(c), qsatg, &
                    qsdT = qsatgdT_h2osfc)
               qg_h2osfc(c) = qsatg
               dqgdT(c) = dqgdT(c) + frac_h2osfc(c) * qsatgdT_h2osfc
            else
               qg_h2osfc(c) = qg_soil(c)
            end if

            qg(c) = frac_sno_eff(c)*qg_snow(c) + (1._r8 - frac_sno_eff(c) - frac_h2osfc(c))*qg_soil(c) &
                 + frac_h2osfc(c) * qg_h2osfc(c)

            ! [PORTED by Hui Tang: NVP ground evaporation blending]
            ! When NVP is active, include NVP surface humidity in qg blend.
            ! NVP occupies area not covered by snow or surface water.
            ! qg_nvp = hr_nvp * qsat(t_nvp): hr_nvp acts as surface RH.
            if (use_nvp) then
               qg_nvp(c) = qg_soil(c)  ! default when no NVP coverage; recomputed below if frac_nvp_eff > 0
               if (frac_nvp_eff > 0._r8) then
                  call QSat(t_nvp_col(c), forc_pbot(c), qsatg, &
                       qsdT = qsatgdT_nvp)
                  qg_nvp(c) = hr_nvp * qsatg
                  ! Adjust qg and dqgdT: reduce bare-soil contribution by frac_nvp_eff, add NVP term
                  qg(c) = qg(c) - frac_nvp_eff * qg_soil(c) + frac_nvp_eff * qg_nvp(c)
                  dqgdT(c) = dqgdT(c) - frac_nvp_eff * hr * qsatgdT_soil &
                             + frac_nvp_eff * hr_nvp * qsatgdT_nvp
               end if
            else
               qg_nvp(c) = qg_soil(c)
            end if

         else
            call QSat(t_grnd(c), forc_pbot(c), qsatg, &
                 qsdT = qsatgdT)
            qg(c) = qred*qsatg
            dqgdT(c) = qred*qsatgdT

            if (qsatg > forc_q(c) .and. forc_q(c) > qred*qsatg) then
               qg(c) = forc_q(c)
               dqgdT(c) = 0._r8
            end if

            qg_snow(c) = qg(c)
            qg_soil(c) = qg(c)
            qg_h2osfc(c) = qg(c)
         endif

      end do ! (end of columns loop)


    end associate

  end subroutine CalculateSurfaceHumidity

end module SurfaceHumidityMod
