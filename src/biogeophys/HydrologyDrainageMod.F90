module HydrologyDrainageMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates soil/snow hydrology with drainage (subsurface runoff)
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_varctl        , only : iulog, use_vichydro
  use clm_varcon        , only : e_ice, denh2o, denice, rpi, spval
  use atm2lndType       , only : atm2lnd_type
  use glc2lndMod        , only : glc2lnd_type
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use TemperatureType   , only : temperature_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use IrrigationMod     , only : irrigation_type
  use GlacierSurfaceMassBalanceMod, only : glacier_smb_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: HydrologyDrainage ! Calculates soil/snow hydrolog with drainage
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine HydrologyDrainage(bounds,               &
       num_nolakec, filter_nolakec,                  &
       num_hydrologyc, filter_hydrologyc,            &
       num_urbanc, filter_urbanc,                    &
       num_do_smb_c, filter_do_smb_c,                &
       atm2lnd_inst, glc2lnd_inst, temperature_inst, &
       soilhydrology_inst, soilstate_inst, waterstate_inst, waterflux_inst, &
       irrigation_inst, glacier_smb_inst)
    !
    ! !DESCRIPTION:
    ! Calculates soil/snow hydrology with drainage (subsurface runoff)
    !
    ! !USES:
    use landunit_varcon  , only : istice, istwet, istsoil, istice_mec, istcrop
    use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, icol_shadewall
    use clm_varcon       , only : denh2o, denice
    use clm_varctl       , only : use_vichydro
    use clm_varpar       , only : nlevgrnd, nlevurb
    use clm_time_manager , only : get_step_size, get_nstep
    use SoilHydrologyMod , only : CLMVICMap, Drainage, PerchedLateralFlow, LateralFlowPowerLaw
    use SoilWaterMovementMod , only : use_aquifer_layer
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds               
    integer                  , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                  , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    integer                  , intent(in)    :: num_do_smb_c         ! number of columns in which SMB is calculated, in column filter    
    integer                  , intent(in)    :: filter_do_smb_c(:)   ! column filter for bare landwhere SMB is calculated
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_inst
    type(glc2lnd_type)       , intent(in)    :: glc2lnd_inst
    type(temperature_type)   , intent(in)    :: temperature_inst
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    type(soilstate_type)     , intent(inout) :: soilstate_inst
    type(waterstate_type)    , intent(inout) :: waterstate_inst
    type(waterflux_type)     , intent(inout) :: waterflux_inst
    type(irrigation_type)    , intent(in)    :: irrigation_inst
    type(glacier_smb_type)   , intent(in)    :: glacier_smb_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: g,l,c,j,fc                 ! indices
    real(r8) :: dtime                      ! land model time step (sec)
    !-----------------------------------------------------------------------
    
    associate(                                                         &    
         dz                 => col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer thickness depth (m)                       
         ctype              => col%itype                             , & ! Input:  [integer  (:)   ]  column type                                        

         qflx_floodg        => atm2lnd_inst%forc_flood_grc           , & ! Input:  [real(r8) (:)   ]  gridcell flux of flood water from RTM             
         forc_rain          => atm2lnd_inst%forc_rain_downscaled_col , & ! Input:  [real(r8) (:)   ]  rain rate [mm/s]                                  
         forc_snow          => atm2lnd_inst%forc_snow_downscaled_col , & ! Input:  [real(r8) (:)   ]  snow rate [mm/s]                                  

         wa                 => soilhydrology_inst%wa_col             , & ! Input:  [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
         
         h2ocan             => waterstate_inst%h2ocan_col            , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O)                             
         h2osfc             => waterstate_inst%h2osfc_col            , & ! Input:  [real(r8) (:)   ]  surface water (mm)                                
         h2osno             => waterstate_inst%h2osno_col            , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                               
         begwb              => waterstate_inst%begwb_col             , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step              
         endwb              => waterstate_inst%endwb_col             , & ! Output: [real(r8) (:)   ]  water mass end of the time step                   
         h2osoi_ice         => waterstate_inst%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
         h2osoi_liq         => waterstate_inst%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
         h2osoi_vol         => waterstate_inst%h2osoi_vol_col        , & ! Output: [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]

         qflx_evap_tot      => waterflux_inst%qflx_evap_tot_col      , & ! Input:  [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg     
         qflx_snwcp_ice     => waterflux_inst%qflx_snwcp_ice_col     , & ! Input: [real(r8) (:)   ]  excess solid h2o due to snow capping (outgoing) (mm H2O /s) [+]`
         qflx_h2osfc_surf   => waterflux_inst%qflx_h2osfc_surf_col   , & ! Output: [real(r8) (:)   ]  surface water runoff (mm/s)                        
         qflx_drain_perched => waterflux_inst%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ]  sub-surface runoff from perched zwt (mm H2O /s)   
         qflx_rsub_sat      => waterflux_inst%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]                 
         qflx_drain         => waterflux_inst%qflx_drain_col         , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf          => waterflux_inst%qflx_surf_col          , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)                        
         qflx_infl          => waterflux_inst%qflx_infl_col          , & ! Output: [real(r8) (:)   ]  infiltration (mm H2O /s)                          
         qflx_qrgwl         => waterflux_inst%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes            
         qflx_runoff        => waterflux_inst%qflx_runoff_col        , & ! Output: [real(r8) (:)   ]  total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
         qflx_runoff_u      => waterflux_inst%qflx_runoff_u_col      , & ! Output: [real(r8) (:)   ]  Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
         qflx_runoff_r      => waterflux_inst%qflx_runoff_r_col      , & ! Output: [real(r8) (:)   ]  Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
         qflx_ice_runoff_snwcp => waterflux_inst%qflx_ice_runoff_snwcp_col, & ! Output: [real(r8) (:)] solid runoff from snow capping (mm H2O /s)
         qflx_irrig         => irrigation_inst%qflx_irrig_col          & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)                       
         )

      ! Determine time step and step size

      dtime = get_step_size()

      if (use_vichydro) then
         call CLMVICMap(bounds, num_hydrologyc, filter_hydrologyc, &
              soilhydrology_inst, waterstate_inst)
      endif

      if (use_aquifer_layer()) then 
         call Drainage(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc,&
              temperature_inst, soilhydrology_inst, soilstate_inst, &
              waterstate_inst, waterflux_inst)
      else
         
         call PerchedLateralFlow(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc,&
              soilhydrology_inst, soilstate_inst, &
              waterstate_inst, waterflux_inst)

         
         call LateralFlowPowerLaw(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc,&
              soilhydrology_inst, soilstate_inst, &
              waterstate_inst, waterflux_inst)

      endif

      do j = 1, nlevgrnd
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof) .and. j > nlevurb) then
            else
               h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
            end if
         end do
      end do

      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)

         if (ctype(c) == icol_roof .or. ctype(c) == icol_sunwall &
              .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_road_imperv) then
            endwb(c) = h2ocan(c) + h2osno(c)
         else
            ! add h2osfc to water balance
            endwb(c) = h2ocan(c) + h2osno(c) + h2osfc(c) + wa(c)

         end if
      end do

      do j = 1, nlevgrnd
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof) .and. j > nlevurb) then

            else
               endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
            end if
         end do
      end do

      ! Determine wetland and land ice hydrology (must be placed here
      ! since need snow updated from CombineSnowLayers)

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)
         g = col%gridcell(c)

         if (lun%itype(l)==istwet .or. lun%itype(l)==istice      &
                                  .or. lun%itype(l)==istice_mec) then

            qflx_drain(c)         = 0._r8
            qflx_drain_perched(c) = 0._r8
            qflx_h2osfc_surf(c)   = 0._r8
            qflx_surf(c)          = 0._r8
            qflx_infl(c)          = 0._r8
            qflx_qrgwl(c) = forc_rain(c) + forc_snow(c) + qflx_floodg(g) - qflx_evap_tot(c) - qflx_snwcp_ice(c) - &
                 (endwb(c)-begwb(c))/dtime

         else if (lun%urbpoi(l) .and. ctype(c) /= icol_road_perv) then

            qflx_drain_perched(c) = 0._r8
            qflx_h2osfc_surf(c)   = 0._r8
            qflx_rsub_sat(c)      = spval

         end if

         qflx_ice_runoff_snwcp(c) = qflx_snwcp_ice(c)
      end do

      ! This call needs to be here so that it comes after the initial calculation of
      ! qflx_qrgwl and qflx_ice_runoff_snwcp, but before the ues of qflx_qrgwl in
      ! qflx_runoff.
      call glacier_smb_inst%AdjustRunoffTerms(bounds, num_do_smb_c, filter_do_smb_c, &
           glc2lnd_inst, &
           qflx_qrgwl = qflx_qrgwl(bounds%begc:bounds%endc), &
           qflx_ice_runoff_snwcp = qflx_ice_runoff_snwcp(bounds%begc:bounds%endc))

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)

         qflx_runoff(c) = qflx_drain(c) + qflx_surf(c)  + qflx_h2osfc_surf(c) + qflx_qrgwl(c) + qflx_drain_perched(c)

         if ((lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) .and. col%active(c)) then
            qflx_runoff(c) = qflx_runoff(c) - qflx_irrig(c)
         end if
         if (lun%urbpoi(l)) then
            qflx_runoff_u(c) = qflx_runoff(c)
         else if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
            qflx_runoff_r(c) = qflx_runoff(c)
         end if

      end do

    end associate

  end subroutine HydrologyDrainage

end module HydrologyDrainageMod
