module CanopyHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of
  ! (1) water storage of intercepted precipitation
  ! (2) direct throughfall and canopy drainage of precipitation
  ! (3) the fraction of foliage covered by water and the fraction
  !     of foliage that is dry and transpiring.
  ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_sys_mod     , only : shr_sys_flush
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_time_manager, only : get_step_size
  use clm_varctl      , only : iulog
  use column_varcon   , only : icol_sunwall, icol_shadewall
  use subgridAveMod   , only : p2c
  use LandunitType    , only : lun                
  use atm2lndType     , only : atm2lnd_type
  use AerosolMod      , only : aerosol_type
  use CanopyStateType , only : canopystate_type
  use TemperatureType , only : temperature_type
  use WaterType       , only : water_type
  use WaterFluxBulkType       , only : waterfluxbulk_type
  use Wateratm2lndBulkType    , only : wateratm2lndbulk_type
  use WaterStateBulkType      , only : waterstatebulk_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use WaterTracerUtils        , only : CalcTracerFromBulk
  use ColumnType      , only : col, column_type
  use PatchType       , only : patch, patch_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyHydrology_readnl ! Read namelist
  public :: readParams
  public :: CanopyInterceptionAndThroughfall

  type, private :: params_type
     real(r8) :: liq_canopy_storage_scalar  ! Canopy-storage-of-liquid-water parameter (kg/m2)
     real(r8) :: snow_canopy_storage_scalar  ! Canopy-storage-of-snow parameter (kg/m2)
  end type params_type
  type(params_type), private ::  params_inst
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SumFlux_TopOfCanopyInputs                   ! Compute patch-level precipitation inputs for bulk water or one tracer
  private :: BulkFlux_CanopyInterceptionAndThroughfall   ! Compute canopy interception and throughfall for bulk water
  private :: TracerFlux_CanopyInterceptionAndThroughfall ! Calculate canopy interception and throughfall for one tracer
  private :: UpdateState_AddInterceptionToCanopy         ! Update snocan and liqcan based on interception, for bulk or one tracer
  private :: BulkFlux_CanopyExcess                       ! Compute runoff from canopy due to exceeding maximum storage, for bulk
  private :: TracerFlux_CanopyExcess                     ! Calculate runoff from canopy due to exceeding maximum storage, for one tracer
  private :: UpdateState_RemoveCanfallFromCanopy         ! Update snocan and liqcan based on canfall, for bulk or one tracer
  private :: BulkFlux_SnowUnloading                      ! Compute snow unloading for bulk
  private :: TracerFlux_SnowUnloading                    ! Compute snow unloading for one tracer
  private :: UpdateState_RemoveSnowUnloading             ! Update snocan based on snow unloading, for bulk or one tracer
  private :: SumFlux_FluxesOntoGround                    ! Compute summed fluxes onto ground, for bulk or one tracer
  private :: BulkDiag_FracWet                            ! Determine fraction of vegetated surface that is wet
  !
  ! !PRIVATE DATA MEMBERS:
  real(r8) :: interception_fraction ! Fraction of intercepted precipitation
  real(r8) :: maximum_leaf_wetted_fraction ! Maximum fraction of leaf that may be wet
  logical, private :: use_clm5_fpi    = .false. ! use clm5 fpi equation

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CanopyHydrology_readnl( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CanopyHydrology
    !
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'CanopyHydrology_readnl'  ! subroutine name
    !-----------------------------------------------------------------------
    namelist /clm_canopyhydrology_inparm/ &
         interception_fraction, &
         maximum_leaf_wetted_fraction, &
         use_clm5_fpi

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input. 
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm_CanopyHydrology_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_canopyhydrology_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_canopyhydrology_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clm_canopyhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR finding clm_canopyhydrology_inparm namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(interception_fraction, mpicom)
    call shr_mpi_bcast(maximum_leaf_wetted_fraction, mpicom)
    call shr_mpi_bcast(use_clm5_fpi, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'canopyhydrology settings:'
       write(iulog,*) '  interception_fraction        = ',interception_fraction
       write(iulog,*) '  maximum_leaf_wetted_fraction = ',maximum_leaf_wetted_fraction
       write(iulog,*) '  use_clm5_fpi                 = ',use_clm5_fpi
    endif

   end subroutine CanopyHydrology_readnl

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
    character(len=*), parameter :: subname = 'readParams_CanopyHydrology'
    !--------------------------------------------------------------------

    ! Canopy-storage-of-liquid-water parameter (kg/m2)
    call readNcdioScalar(ncid, 'liq_canopy_storage_scalar', subname, params_inst%liq_canopy_storage_scalar)
    ! Canopy-storage-of-snow parameter (kg/m2)
    call readNcdioScalar(ncid, 'snow_canopy_storage_scalar', subname, params_inst%snow_canopy_storage_scalar)

   end subroutine readParams

   !-----------------------------------------------------------------------
   subroutine CanopyInterceptionAndThroughfall(bounds, &
        num_soilp, filter_soilp, &
        num_nolakep, filter_nolakep, &
        num_nolakec, filter_nolakec, &
        patch, col, &
        canopystate_inst, atm2lnd_inst, water_inst)
     !
     ! !DESCRIPTION:
     ! Coordinate work related to the calculation of canopy interception and throughfall,
     ! as well as removal of excess water from the canopy and snow unloading. This
     ! includes:
     ! (1) water storage of intercepted precipitation
     ! (2) direct throughfall and canopy drainage of precipitation
     ! (3) the fraction of foliage covered by water and the fraction
     !     of foliage that is dry and transpiring.
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds
     integer                , intent(in)    :: num_soilp         ! number of patches in filter_soilp
     integer                , intent(in)    :: filter_soilp(:)   ! patch filter for soil points
     integer                , intent(in)    :: num_nolakep       ! number of patches in filter_nolakep
     integer                , intent(in)    :: filter_nolakep(:) ! patch filter for non-lake points
     integer                , intent(in)    :: num_nolakec       ! number of columns in filter_nolakec
     integer                , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
     type(patch_type)       , intent(in)    :: patch
     type(column_type)      , intent(in)    :: col
     type(canopystate_type) , intent(in)    :: canopystate_inst
     type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
     type(water_type)       , intent(inout) :: water_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: i     ! index of water tracer or bulk
     real(r8) :: dtime ! land model time step (sec)

     real(r8) :: qflx_liq_above_canopy_patch(bounds%begp:bounds%endp)        ! liquid water input above canopy (rain plus irrigation) [mm/s]
     real(r8) :: tracer_qflx_liq_above_canopy_patch(bounds%begp:bounds%endp) ! For one tracer: liquid water input above canopy (rain plus irrigation) [mm/s]
     real(r8) :: forc_snow_patch(bounds%begp:bounds%endp)                    ! atm snow, patch-level [mm/s]
     real(r8) :: tracer_forc_snow_patch(bounds%begp:bounds%endp)             ! For one tracer: atm snow, patch-level [mm/s]

     logical  :: check_point_for_interception_and_excess(bounds%begp:bounds%endp)

     character(len=*), parameter :: subname = 'CanopyInterceptionAndThroughfall'
     !-----------------------------------------------------------------------

     associate( &
          begp => bounds%begp, &
          endp => bounds%endp, &
          begc => bounds%begc, &
          endc => bounds%endc, &
          begg => bounds%begg, &
          endg => bounds%endg, &
          b_wateratm2lnd_inst    => water_inst%wateratm2lndbulk_inst, &
          b_waterflux_inst       => water_inst%waterfluxbulk_inst, &
          b_waterstate_inst      => water_inst%waterstatebulk_inst, &
          b_waterdiagnostic_inst => water_inst%waterdiagnosticbulk_inst &
          )

     dtime = get_step_size()

     ! Note about filters in this routine: Most of the work here just needs to be done
     ! for patches that have some canopy, so we use the soil filter. However, the final
     ! fluxes of water onto the ground (qflx_snow_grnd_col and qflx_liq_grnd_col) are
     ! needed for all non-lake points. So a few routines use the nolake filter to ensure
     ! that these fluxes are set correctly for all patches.


     ! Compute canopy interception and throughfall for bulk water
     !
     ! Note use of nolake filter: We need qflx_through_snow and qflx_through_liq for all
     ! nolake points because these are inputs into qflx_snow_grnd_col and
     ! qflx_liq_grnd_col, which are needed for all nolake points.
     call SumFlux_TopOfCanopyInputs(bounds, num_nolakep, filter_nolakep, &
          ! Inputs
          patch                 = patch, &
          forc_rain             = b_wateratm2lnd_inst%forc_rain_downscaled_col(begc:endc), &
          qflx_irrig_sprinkler  = b_waterflux_inst%qflx_irrig_sprinkler_patch(begp:endp), &
          forc_snow_col         = b_wateratm2lnd_inst%forc_snow_downscaled_col(begc:endc), &
          ! Outputs
          qflx_liq_above_canopy = qflx_liq_above_canopy_patch(begp:endp), &
          forc_snow_patch       = forc_snow_patch(begp:endp))

     call BulkFlux_CanopyInterceptionAndThroughfall(bounds, num_nolakep, filter_nolakep, &
          ! Inputs
          patch                 = patch, &
          col                   = col, &
          frac_veg_nosno        = canopystate_inst%frac_veg_nosno_patch(begp:endp), &
          elai                  = canopystate_inst%elai_patch(begp:endp), &
          esai                  = canopystate_inst%esai_patch(begp:endp), &
          forc_snow             = forc_snow_patch(begp:endp), &
          qflx_liq_above_canopy = qflx_liq_above_canopy_patch(begp:endp), &
          ! Outputs
          qflx_through_snow     = b_waterflux_inst%qflx_through_snow_patch(begp:endp), &
          qflx_through_liq      = b_waterflux_inst%qflx_through_liq_patch(begp:endp), &
          qflx_intercepted_snow = b_waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
          qflx_intercepted_liq  = b_waterflux_inst%qflx_intercepted_liq_patch(begp:endp), &
          check_point_for_interception_and_excess = check_point_for_interception_and_excess(begp:endp))

     ! Calculate canopy interception and throughfall for each tracer
     !
     ! Note use of nolake filter: We need qflx_through_snow and qflx_through_liq for all
     ! nolake points because these are inputs into qflx_snow_grnd_col and
     ! qflx_liq_grnd_col, which are needed for all nolake points.
     do i = water_inst%tracers_beg, water_inst%tracers_end
        associate(w => water_inst%bulk_and_tracers(i))
        call SumFlux_TopOfCanopyInputs(bounds, num_nolakep, filter_nolakep, &
             ! Inputs
             patch                 = patch, &
             forc_rain             = w%wateratm2lnd_inst%forc_rain_downscaled_col(begc:endc), &
             qflx_irrig_sprinkler  = w%waterflux_inst%qflx_irrig_sprinkler_patch(begp:endp), &
             forc_snow_col         = w%wateratm2lnd_inst%forc_snow_downscaled_col(begc:endc), &
             ! Outputs
             qflx_liq_above_canopy = tracer_qflx_liq_above_canopy_patch(begp:endp), &
             forc_snow_patch       = tracer_forc_snow_patch(begp:endp))

        call TracerFlux_CanopyInterceptionAndThroughfall(bounds, num_nolakep, filter_nolakep, &
             ! Inputs
             bulk_forc_snow             = forc_snow_patch(begp:endp), &
             bulk_qflx_liq_above_canopy = qflx_liq_above_canopy_patch(begp:endp), &
             bulk_qflx_through_snow     = b_waterflux_inst%qflx_through_snow_patch(begp:endp), &
             bulk_qflx_intercepted_snow = b_waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
             bulk_qflx_through_liq      = b_waterflux_inst%qflx_through_liq_patch(begp:endp), &
             bulk_qflx_intercepted_liq  = b_waterflux_inst%qflx_intercepted_liq_patch(begp:endp), &
             trac_forc_snow             = tracer_forc_snow_patch(begp:endp), &
             trac_qflx_liq_above_canopy = tracer_qflx_liq_above_canopy_patch(begp:endp), &
             ! Outputs
             trac_qflx_through_snow     = w%waterflux_inst%qflx_through_snow_patch(begp:endp), &
             trac_qflx_intercepted_snow = w%waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
             trac_qflx_through_liq      = w%waterflux_inst%qflx_through_liq_patch(begp:endp), &
             trac_qflx_intercepted_liq  = w%waterflux_inst%qflx_intercepted_liq_patch(begp:endp))
        end associate
     end do

     ! Update snocan and liqcan based on interception, for bulk water and each tracer
     do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
        associate(w => water_inst%bulk_and_tracers(i))
        call UpdateState_AddInterceptionToCanopy(bounds, num_soilp, filter_soilp, &
             ! Inputs
             dtime                 = dtime, &
             qflx_intercepted_snow = w%waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
             qflx_intercepted_liq  = w%waterflux_inst%qflx_intercepted_liq_patch(begp:endp), &
             ! Outputs
             snocan                = w%waterstate_inst%snocan_patch(begp:endp), &
             liqcan                = w%waterstate_inst%liqcan_patch(begp:endp))
        end associate
     end do

     ! Compute runoff from canopy due to exceeding maximum storage, for bulk
     call BulkFlux_CanopyExcess(bounds, num_soilp, filter_soilp, &
          ! Inputs
          dtime = dtime, &
          elai = canopystate_inst%elai_patch(begp:endp), &
          esai = canopystate_inst%esai_patch(begp:endp), &
          snocan = b_waterstate_inst%snocan_patch(begp:endp), &
          liqcan = b_waterstate_inst%liqcan_patch(begp:endp), &
          check_point_for_interception_and_excess = check_point_for_interception_and_excess(begp:endp), &
          ! Outputs
          qflx_snocanfall = b_waterflux_inst%qflx_snocanfall_patch(begp:endp), &
          qflx_liqcanfall = b_waterflux_inst%qflx_liqcanfall_patch(begp:endp))

     ! Calculate runoff from canopy due to exceeding maximum storage, for each tracer
     do i = water_inst%tracers_beg, water_inst%tracers_end
        associate(w => water_inst%bulk_and_tracers(i))
        call TracerFlux_CanopyExcess(bounds, num_soilp, filter_soilp, &
             ! Inputs
             bulk_liqcan          = b_waterstate_inst%liqcan_patch(begp:endp), &
             bulk_snocan          = b_waterstate_inst%snocan_patch(begp:endp), &
             bulk_qflx_liqcanfall = b_waterflux_inst%qflx_liqcanfall_patch(begp:endp), &
             bulk_qflx_snocanfall = b_waterflux_inst%qflx_snocanfall_patch(begp:endp), &
             trac_liqcan          = w%waterstate_inst%liqcan_patch(begp:endp), &
             trac_snocan          = w%waterstate_inst%snocan_patch(begp:endp), &
             ! Outputs
             trac_qflx_liqcanfall = w%waterflux_inst%qflx_liqcanfall_patch(begp:endp), &
             trac_qflx_snocanfall = w%waterflux_inst%qflx_snocanfall_patch(begp:endp))
        end associate
     end do

     ! Update snocan and liqcan based on canfall, for bulk water and each tracer
     do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
        associate(w => water_inst%bulk_and_tracers(i))
        call UpdateState_RemoveCanfallFromCanopy(bounds, num_soilp, filter_soilp, &
             ! Inputs
             dtime           = dtime, &
             qflx_liqcanfall = w%waterflux_inst%qflx_liqcanfall_patch(begp:endp), &
             qflx_snocanfall = w%waterflux_inst%qflx_snocanfall_patch(begp:endp), &
             ! Outputs
             liqcan          = w%waterstate_inst%liqcan_patch(begp:endp), &
             snocan          = w%waterstate_inst%snocan_patch(begp:endp))
        end associate
     end do

     ! Compute snow unloading for bulk
     call BulkFlux_SnowUnloading(bounds, num_soilp, filter_soilp, &
          ! Inputs
          dtime              = dtime, &
          patch              = patch, &
          frac_veg_nosno     = canopystate_inst%frac_veg_nosno_patch(begp:endp), &
          forc_t             = atm2lnd_inst%forc_t_downscaled_col(begc:endc), &
          forc_wind          = atm2lnd_inst%forc_wind_grc(begg:endg), &
          snocan             = b_waterstate_inst%snocan_patch(begp:endp), &
          ! Outputs
          qflx_snotempunload = b_waterflux_inst%qflx_snotempunload_patch(begp:endp), &
          qflx_snowindunload = b_waterflux_inst%qflx_snowindunload_patch(begp:endp), &
          qflx_snow_unload   = b_waterflux_inst%qflx_snow_unload_patch(begp:endp))

     ! Compute snow unloading for each tracer
     do i = water_inst%tracers_beg, water_inst%tracers_end
        associate(w => water_inst%bulk_and_tracers(i))
        call TracerFlux_SnowUnloading(bounds, num_soilp, filter_soilp, &
             ! Inputs
             bulk_snocan           = b_waterstate_inst%snocan_patch(begp:endp), &
             bulk_qflx_snow_unload = b_waterflux_inst%qflx_snow_unload_patch(begp:endp), &
             trac_snocan           = w%waterstate_inst%snocan_patch(begp:endp), &
             ! Outputs
             trac_qflx_snow_unload = w%waterflux_inst%qflx_snow_unload_patch(begp:endp))
        end associate
     end do

     ! Update snocan based on snow unloading, for bulk water and each tracer
     do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
        associate(w => water_inst%bulk_and_tracers(i))
        call UpdateState_RemoveSnowUnloading(bounds, num_soilp, filter_soilp, &
             ! Inputs
             dtime            = dtime, &
             qflx_snow_unload = w%waterflux_inst%qflx_snow_unload_patch(begp:endp), &
             ! Outputs
             snocan           = w%waterstate_inst%snocan_patch(begp:endp))
        end associate
     end do

     ! Compute summed fluxes onto ground, for bulk water and each tracer
     !
     ! Note use of nolake filter: We need qflx_snow_grnd_col and qflx_liq_grnd_col for
     ! all non-lake points.
     do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
        associate(w => water_inst%bulk_and_tracers(i))
        call SumFlux_FluxesOntoGround(bounds, &
             num_nolakep, filter_nolakep, &
             num_nolakec, filter_nolakec, &
             ! Inputs
             qflx_through_snow  = w%waterflux_inst%qflx_through_snow_patch(begp:endp), &
             qflx_snocanfall    = w%waterflux_inst%qflx_snocanfall_patch(begp:endp), &
             qflx_snow_unload   = w%waterflux_inst%qflx_snow_unload_patch(begp:endp), &
             qflx_through_liq   = w%waterflux_inst%qflx_through_liq_patch(begp:endp), &
             qflx_liqcanfall    = w%waterflux_inst%qflx_liqcanfall_patch(begp:endp), &
             qflx_irrig_drip    = w%waterflux_inst%qflx_irrig_drip_patch(begp:endp), &
             ! Outputs
             qflx_snow_grnd_col = w%waterflux_inst%qflx_snow_grnd_col(begc:endc), &
             qflx_liq_grnd_col  = w%waterflux_inst%qflx_liq_grnd_col(begc:endc))
        end associate
     end do

     ! Determine the fraction of foliage covered by water and the fraction of foliage that
     ! is dry and transpiring.
     call BulkDiag_FracWet(bounds, num_soilp, filter_soilp, &
          ! Inputs
          frac_veg_nosno = canopystate_inst%frac_veg_nosno_patch(begp:endp), &
          elai           = canopystate_inst%elai_patch(begp:endp), &
          esai           = canopystate_inst%esai_patch(begp:endp), &
          snocan         = b_waterstate_inst%snocan_patch(begp:endp), &
          liqcan         = b_waterstate_inst%liqcan_patch(begp:endp), &
          ! Outputs
          fwet           = b_waterdiagnostic_inst%fwet_patch(begp:endp), &
          fdry           = b_waterdiagnostic_inst%fdry_patch(begp:endp), &
          fcansno        = b_waterdiagnostic_inst%fcansno_patch(begp:endp))

     end associate

   end subroutine CanopyInterceptionAndThroughfall

   !-----------------------------------------------------------------------
   subroutine SumFlux_TopOfCanopyInputs(bounds, num_nolakep, filter_nolakep, &
        patch, forc_rain, qflx_irrig_sprinkler, forc_snow_col, &
        qflx_liq_above_canopy, forc_snow_patch)
     !
     ! !DESCRIPTION:
     ! Compute patch-level precipitation inputs for bulk water or one tracer
     !
     ! !ARGUMENTS:
     type(bounds_type) , intent(in)    :: bounds
     integer           , intent(in)    :: num_nolakep
     integer           , intent(in)    :: filter_nolakep(:)
     type(patch_type)  , intent(in)    :: patch
     real(r8)          , intent(in)    :: forc_rain( bounds%begc: )             ! atm rain rate [mm/s]
     real(r8)          , intent(in)    :: qflx_irrig_sprinkler( bounds%begp: )  ! sprinkler irrigation [mm/s]
     real(r8)          , intent(in)    :: forc_snow_col( bounds%begc: )         ! atm snow rate [mm/s]
     real(r8)          , intent(inout) :: qflx_liq_above_canopy( bounds%begp: ) ! total incoming liquid water above canopy [mm/s]
     real(r8)          , intent(inout) :: forc_snow_patch( bounds%begp: )       ! atm snow rate, patch-level [mm/s]
     !
     ! !LOCAL VARIABLES:
     integer :: fp, p, c

     character(len=*), parameter :: subname = 'SumFlux_TopOfCanopyInputs'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(forc_rain, 1) == bounds%endc), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_irrig_sprinkler, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(forc_snow_col, 1) == bounds%endc), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_liq_above_canopy, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(forc_snow_patch, 1) == bounds%endp), sourcefile, __LINE__)

     do fp = 1, num_nolakep
        p = filter_nolakep(fp)
        c = patch%column(p)

        qflx_liq_above_canopy(p) = forc_rain(c) + qflx_irrig_sprinkler(p)
        forc_snow_patch(p) = forc_snow_col(c)
     end do

   end subroutine SumFlux_TopOfCanopyInputs

   !-----------------------------------------------------------------------
   subroutine BulkFlux_CanopyInterceptionAndThroughfall(bounds, num_nolakep, filter_nolakep, &
        patch, col, &
        frac_veg_nosno, elai, esai, forc_snow, qflx_liq_above_canopy, &
        qflx_through_snow, qflx_through_liq, &
        qflx_intercepted_snow, qflx_intercepted_liq, &
        check_point_for_interception_and_excess)
     !
     ! !DESCRIPTION:
     ! Compute canopy interception and throughfall for bulk water
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_nolakep
     integer, intent(in) :: filter_nolakep(:)
     type(patch_type), intent(in) :: patch
     type(column_type), intent(in) :: col

     integer  , intent(in)    :: frac_veg_nosno( bounds%begp: )                          ! fraction of vegetation not covered by snow (0 OR 1)
     real(r8) , intent(in)    :: elai( bounds%begp: )                                    ! canopy one-sided leaf area index with burying by snow
     real(r8) , intent(in)    :: esai( bounds%begp: )                                    ! canopy one-sided stem area index with burying by snow
     real(r8) , intent(in)    :: forc_snow( bounds%begp: )                               ! atm snow (mm H2O/s)
     real(r8) , intent(in)    :: qflx_liq_above_canopy( bounds%begp: )                   ! liquid water input above canopy (rain plus irrigation) (mm H2O/s)

     real(r8) , intent(inout) :: qflx_through_snow( bounds%begp: )                       ! canopy throughfall of snow (mm H2O/s)
     real(r8) , intent(inout) :: qflx_through_liq( bounds%begp: )                        ! canopy throughfall of liquid (mm H2O/s)
     real(r8) , intent(inout) :: qflx_intercepted_snow( bounds%begp: )                   ! canopy interception of snow (mm H2O/s)
     real(r8) , intent(inout) :: qflx_intercepted_liq( bounds%begp: )                    ! canopy interception of liquid (mm H2O/s)
     logical  , intent(inout) :: check_point_for_interception_and_excess( bounds%begp: ) ! whether each patch in the filter needs to have the interception calculations (here) and snow/liquid excess calculations (elsewhere) computed
     !
     ! !LOCAL VARIABLES:
     integer :: fp, p, c
     real(r8) :: fpiliq  ! coefficient of interception for liquid
     real(r8) :: fpisnow ! coefficient of interception for snow

     character(len=*), parameter :: subname = 'BulkFlux_CanopyInterceptionAndThroughfall'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(frac_veg_nosno, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(elai, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(esai, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(forc_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_liq_above_canopy, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_through_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_through_liq, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_intercepted_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_intercepted_liq, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(check_point_for_interception_and_excess, 1) == bounds%endp), sourcefile, __LINE__)

     do fp = 1, num_nolakep
        p = filter_nolakep(fp)
        c = patch%column(p)

        check_point_for_interception_and_excess(p) = &
             (frac_veg_nosno(p) == 1 .and. (forc_snow(p) + qflx_liq_above_canopy(p)) > 0._r8)
        if (check_point_for_interception_and_excess(p)) then
           ! Coefficient of interception
           if (use_clm5_fpi) then
              fpiliq = interception_fraction * tanh(elai(p) + esai(p))
           else
              fpiliq = 0.25_r8*(1._r8 - exp(-0.5_r8*(elai(p) + esai(p))))
           end if

           fpisnow = (1._r8 - exp(-0.5_r8*(elai(p) + esai(p))))  ! max interception of 1

           ! Direct throughfall
           qflx_through_snow(p) = forc_snow(p) * (1._r8-fpisnow)
           qflx_through_liq(p)  = qflx_liq_above_canopy(p) * (1._r8-fpiliq)

           ! Canopy interception
           qflx_intercepted_snow(p) = forc_snow(p) * fpisnow
           qflx_intercepted_liq(p) = qflx_liq_above_canopy(p) * fpiliq

        else
           ! Note that special landunits will be handled here, in addition to soil points
           ! with frac_veg_nosno == 0.
           qflx_intercepted_snow(p) = 0._r8
           qflx_intercepted_liq(p) = 0._r8
           if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
              qflx_through_snow(p) = 0._r8
              qflx_through_liq(p)  = 0._r8
           else
              qflx_through_snow(p) = forc_snow(p)
              qflx_through_liq(p)  = qflx_liq_above_canopy(p)
           end if
        end if
     end do

   end subroutine BulkFlux_CanopyInterceptionAndThroughfall

   !-----------------------------------------------------------------------
   subroutine TracerFlux_CanopyInterceptionAndThroughfall(bounds, num_nolakep, filter_nolakep, &
        bulk_forc_snow, bulk_qflx_liq_above_canopy, &
        bulk_qflx_through_snow, bulk_qflx_intercepted_snow, &
        bulk_qflx_through_liq, bulk_qflx_intercepted_liq, &
        trac_forc_snow, trac_qflx_liq_above_canopy, &
        trac_qflx_through_snow, trac_qflx_intercepted_snow, &
        trac_qflx_through_liq, trac_qflx_intercepted_liq)
     !
     ! !DESCRIPTION:
     ! Calculate canopy interception and throughfall for one tracer
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_nolakep
     integer, intent(in) :: filter_nolakep(:)

     ! For description of arguments, see comments in
     ! BulkFlux_CanopyInterceptionAndThroughfall. Here, bulk_* variables refer to bulk
     ! water and trac_* variables to the given water tracer.
     real(r8), intent(in) :: bulk_forc_snow( bounds%begp: )
     real(r8), intent(in) :: bulk_qflx_liq_above_canopy( bounds%begp: )
     real(r8), intent(in) :: bulk_qflx_through_snow( bounds%begp: )
     real(r8), intent(in) :: bulk_qflx_intercepted_snow( bounds%begp: )
     real(r8), intent(in) :: bulk_qflx_through_liq( bounds%begp: )
     real(r8), intent(in) :: bulk_qflx_intercepted_liq( bounds%begp: )
     real(r8), intent(in) :: trac_forc_snow( bounds%begp: )
     real(r8), intent(in) :: trac_qflx_liq_above_canopy( bounds%begp: )

     real(r8), intent(inout) :: trac_qflx_through_snow( bounds%begp: )
     real(r8), intent(inout) :: trac_qflx_intercepted_snow( bounds%begp: )
     real(r8), intent(inout) :: trac_qflx_through_liq( bounds%begp: )
     real(r8), intent(inout) :: trac_qflx_intercepted_liq( bounds%begp: )
     !
     ! !LOCAL VARIABLES:

     character(len=*), parameter :: subname = 'TracerFlux_CanopyInterceptionAndThroughfall'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(bulk_forc_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_qflx_liq_above_canopy, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_qflx_through_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_qflx_intercepted_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_qflx_through_liq, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_qflx_intercepted_liq, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_forc_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_qflx_liq_above_canopy, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_qflx_through_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_qflx_intercepted_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_qflx_through_liq, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_qflx_intercepted_liq, 1) == bounds%endp), sourcefile, __LINE__)

     associate( &
          begp => bounds%begp, &
          endp => bounds%endp  &
          )

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_nolakep, &
          filter_pts    = filter_nolakep, &
          bulk_source   = bulk_forc_snow(begp:endp), &
          bulk_val      = bulk_qflx_through_snow(begp:endp), &
          tracer_source = trac_forc_snow(begp:endp), &
          tracer_val    = trac_qflx_through_snow(begp:endp))

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_nolakep, &
          filter_pts    = filter_nolakep, &
          bulk_source   = bulk_forc_snow(begp:endp), &
          bulk_val      = bulk_qflx_intercepted_snow(begp:endp), &
          tracer_source = trac_forc_snow(begp:endp), &
          tracer_val    = trac_qflx_intercepted_snow(begp:endp))

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_nolakep, &
          filter_pts    = filter_nolakep, &
          bulk_source   = bulk_qflx_liq_above_canopy(begp:endp), &
          bulk_val      = bulk_qflx_through_liq(begp:endp), &
          tracer_source = trac_qflx_liq_above_canopy(begp:endp), &
          tracer_val    = trac_qflx_through_liq(begp:endp))

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_nolakep, &
          filter_pts    = filter_nolakep, &
          bulk_source   = bulk_qflx_liq_above_canopy(begp:endp), &
          bulk_val      = bulk_qflx_intercepted_liq(begp:endp), &
          tracer_source = trac_qflx_liq_above_canopy(begp:endp), &
          tracer_val    = trac_qflx_intercepted_liq(begp:endp))

     end associate

   end subroutine TracerFlux_CanopyInterceptionAndThroughfall

   !-----------------------------------------------------------------------
   subroutine UpdateState_AddInterceptionToCanopy(bounds, num_soilp, filter_soilp, dtime, &
        qflx_intercepted_snow, qflx_intercepted_liq, snocan, liqcan)
     !
     ! !DESCRIPTION:
     ! Update snocan and liqcan based on interception, for bulk or one tracer
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_soilp
     integer, intent(in) :: filter_soilp(:)
     real(r8), intent(in) :: dtime  ! land model time step (sec)

     real(r8) , intent(in)    :: qflx_intercepted_snow( bounds%begp: ) ! canopy interception of snow (mm H2O/s)
     real(r8) , intent(in)    :: qflx_intercepted_liq( bounds%begp: )  ! canopy interception of liquid (mm H2O/s)

     real(r8) , intent(inout) :: snocan( bounds%begp: )                ! canopy snow water (mm H2O)
     real(r8) , intent(inout) :: liqcan( bounds%begp: )                ! canopy liquid water (mm H2O)
     !
     ! !LOCAL VARIABLES:
     integer :: fp, p

     character(len=*), parameter :: subname = 'UpdateState_AddInterceptionToCanopy'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(qflx_intercepted_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_intercepted_liq, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(snocan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(liqcan, 1) == bounds%endp), sourcefile, __LINE__)

     do fp = 1, num_soilp
        p = filter_soilp(fp)

        snocan(p) = max(0._r8, snocan(p) + dtime * qflx_intercepted_snow(p))
        liqcan(p) = max(0._r8, liqcan(p) + dtime * qflx_intercepted_liq(p))
     end do

   end subroutine UpdateState_AddInterceptionToCanopy

   !-----------------------------------------------------------------------
   subroutine BulkFlux_CanopyExcess(bounds, num_soilp, filter_soilp, &
        dtime, elai, esai, snocan, liqcan, &
        check_point_for_interception_and_excess, &
        qflx_snocanfall, qflx_liqcanfall)
     !
     ! !DESCRIPTION:
     ! Compute runoff from canopy due to exceeding maximum storage, for bulk
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_soilp
     integer, intent(in) :: filter_soilp(:)

     real(r8) , intent(in)    :: dtime                                                   ! land model time step (sec)
     real(r8) , intent(in)    :: elai( bounds%begp: )                                    ! canopy one-sided leaf area index with burying by snow
     real(r8) , intent(in)    :: esai( bounds%begp: )                                    ! canopy one-sided stem area index with burying by snow
     real(r8) , intent(in)    :: snocan( bounds%begp: )                                  ! canopy snow water (mm H2O)
     real(r8) , intent(in)    :: liqcan( bounds%begp: )                                  ! canopy liquid water (mm H2O)

     logical  , intent(in)    :: check_point_for_interception_and_excess( bounds%begp: ) ! whether each patch in the filter needs to have the interception calculations (elsewhere) and snow/liquid excess calculations (here) computed

     real(r8) , intent(inout) :: qflx_snocanfall( bounds%begp: )                         ! rate of excess canopy snow falling off canopy (mm H2O/s)
     real(r8) , intent(inout) :: qflx_liqcanfall( bounds%begp: )                         ! rate of excess canopy liquid falling off canopy (mm H2O/s)
     !
     ! !LOCAL VARIABLES:
     integer :: fp, p
     real(r8) :: snocanmx ! maximum allowed snow on canopy (mm H2O)
     real(r8) :: liqcanmx ! maximum allowed liquid water on canopy (mm H2O)

     character(len=*), parameter :: subname = 'BulkFlux_CanopyExcess'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(elai, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(esai, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(snocan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(liqcan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(check_point_for_interception_and_excess, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_snocanfall, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_liqcanfall, 1) == bounds%endp), sourcefile, __LINE__)

     do fp = 1, num_soilp
        p = filter_soilp(fp)
        qflx_liqcanfall(p) = 0._r8
        qflx_snocanfall(p) = 0._r8

        if (check_point_for_interception_and_excess(p)) then
           liqcanmx = params_inst%liq_canopy_storage_scalar * (elai(p) + esai(p))
           qflx_liqcanfall(p) = max((liqcan(p) - liqcanmx)/dtime, 0._r8)
           snocanmx = params_inst%snow_canopy_storage_scalar * (elai(p) + esai(p))
           qflx_snocanfall(p) = max((snocan(p) - snocanmx)/dtime, 0._r8)
        end if
     end do

   end subroutine BulkFlux_CanopyExcess

   !-----------------------------------------------------------------------
   subroutine TracerFlux_CanopyExcess(bounds, num_soilp, filter_soilp, &
        bulk_liqcan, bulk_snocan, &
        bulk_qflx_liqcanfall, bulk_qflx_snocanfall, &
        trac_liqcan, trac_snocan, &
        trac_qflx_liqcanfall, trac_qflx_snocanfall)
     !
     ! !DESCRIPTION:
     ! Calculate runoff from canopy due to exceeding maximum storage, for one tracer
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_soilp
     integer, intent(in) :: filter_soilp(:)

     ! For description of arguments, see comments in BulkFlux_CanopyExcess. Here, bulk_*
     ! variables refer to bulk water and trac_* variables to the given water tracer.
     real(r8), intent(in) :: bulk_liqcan( bounds%begp: )
     real(r8), intent(in) :: bulk_snocan( bounds%begp: )
     real(r8), intent(in) :: bulk_qflx_liqcanfall( bounds%begp: )
     real(r8), intent(in) :: bulk_qflx_snocanfall( bounds%begp: )
     real(r8), intent(in) :: trac_liqcan( bounds%begp: )
     real(r8), intent(in) :: trac_snocan( bounds%begp: )

     real(r8), intent(inout) :: trac_qflx_liqcanfall( bounds%begp: )
     real(r8), intent(inout) :: trac_qflx_snocanfall( bounds%begp: )
     !
     ! !LOCAL VARIABLES:

     character(len=*), parameter :: subname = 'TracerFlux_CanopyExcess'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(bulk_liqcan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_snocan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_qflx_liqcanfall, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_qflx_snocanfall, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_liqcan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_snocan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_qflx_liqcanfall, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_qflx_snocanfall, 1) == bounds%endp), sourcefile, __LINE__)

     associate( &
          begp => bounds%begp, &
          endp => bounds%endp  &
          )

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk_liqcan(begp:endp), &
          bulk_val      = bulk_qflx_liqcanfall(begp:endp), &
          tracer_source = trac_liqcan(begp:endp), &
          tracer_val    = trac_qflx_liqcanfall(begp:endp))

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk_snocan(begp:endp), &
          bulk_val      = bulk_qflx_snocanfall(begp:endp), &
          tracer_source = trac_snocan(begp:endp), &
          tracer_val    = trac_qflx_snocanfall(begp:endp))

     end associate

   end subroutine TracerFlux_CanopyExcess

   !-----------------------------------------------------------------------
   subroutine UpdateState_RemoveCanfallFromCanopy(bounds, num_soilp, filter_soilp, dtime, &
        qflx_liqcanfall, qflx_snocanfall, liqcan, snocan)
     !
     ! !DESCRIPTION:
     ! Update snocan and liqcan based on canfall, for bulk or one tracer
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_soilp
     integer, intent(in) :: filter_soilp(:)
     real(r8), intent(in) :: dtime  ! land model time step (sec)

     real(r8) , intent(in)    :: qflx_liqcanfall( bounds%begp: ) ! rate of excess canopy liquid falling off canopy (mm H2O/s)
     real(r8) , intent(in)    :: qflx_snocanfall( bounds%begp: ) ! rate of excess canopy snow falling off canopy (mm H2O/s)
     real(r8) , intent(inout) :: liqcan( bounds%begp: )          ! canopy liquid water (mm H2O)
     real(r8) , intent(inout) :: snocan( bounds%begp: )          ! canopy snow water (mm H2O)
     !
     ! !LOCAL VARIABLES:
     integer :: fp, p

     character(len=*), parameter :: subname = 'UpdateState_RemoveCanfallFromCanopy'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(qflx_liqcanfall, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_snocanfall, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(liqcan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(snocan, 1) == bounds%endp), sourcefile, __LINE__)

     do fp = 1, num_soilp
        p = filter_soilp(fp)

        liqcan(p) = liqcan(p) - dtime * qflx_liqcanfall(p)
        snocan(p) = snocan(p) - dtime * qflx_snocanfall(p)
     end do

   end subroutine UpdateState_RemoveCanfallFromCanopy

   !-----------------------------------------------------------------------
   subroutine BulkFlux_SnowUnloading(bounds, num_soilp, filter_soilp, dtime, patch, &
        frac_veg_nosno, forc_t, forc_wind, snocan, &
        qflx_snotempunload, qflx_snowindunload, qflx_snow_unload)
     !
     ! !DESCRIPTION:
     ! Compute snow unloading for bulk
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_soilp
     integer, intent(in) :: filter_soilp(:)
     real(r8), intent(in) :: dtime  ! land model time step (sec)
     type(patch_type), intent(in) :: patch

     integer  , intent(in)    :: frac_veg_nosno( bounds%begp: )     ! fraction of vegetation not covered by snow (0 OR 1)
     real(r8) , intent(in)    :: forc_t( bounds%begc: )             ! atmospheric temperature (Kelvin)
     real(r8) , intent(in)    :: forc_wind( bounds%begg: )          ! atmospheric wind speed (m/s)
     real(r8) , intent(in)    :: snocan( bounds%begp: )             ! canopy snow water (mm H2O)

     real(r8) , intent(inout) :: qflx_snotempunload( bounds%begp: ) ! canopy snow unloading from temperature (mm H2O/s)
     real(r8) , intent(inout) :: qflx_snowindunload( bounds%begp: ) ! canopy snow unloading from wind (mm H2O/s)
     real(r8) , intent(inout) :: qflx_snow_unload( bounds%begp: )   ! total canopy snow unloading (mm H2O/s)
     !
     ! !LOCAL VARIABLES:
     integer :: fp, p, c, g

     character(len=*), parameter :: subname = 'BulkFlux_SnowUnloading'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(frac_veg_nosno, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(forc_t, 1) == bounds%endc), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(forc_wind, 1) == bounds%endg), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(snocan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_snotempunload, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_snowindunload, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_snow_unload, 1) == bounds%endp), sourcefile, __LINE__)

     do fp = 1, num_soilp
        p = filter_soilp(fp)

        if (frac_veg_nosno(p) == 1 .and. snocan(p) > 0._r8) then
           c = patch%column(p)
           g = patch%gridcell(p)

           qflx_snotempunload(p) = max(0._r8,snocan(p)*(forc_t(c)-270.15_r8)/1.87e5_r8)
           qflx_snowindunload(p) = 0.5_r8*snocan(p)*forc_wind(g)/1.56e5_r8
           qflx_snow_unload(p) = min(qflx_snotempunload(p) + qflx_snowindunload(p), snocan(p)/dtime)
        else
           qflx_snotempunload(p) = 0._r8
           qflx_snowindunload(p) = 0._r8
           qflx_snow_unload(p) = 0._r8
        end if
     end do

   end subroutine BulkFlux_SnowUnloading

   !-----------------------------------------------------------------------
   subroutine TracerFlux_SnowUnloading(bounds, num_soilp, filter_soilp, &
        bulk_snocan, bulk_qflx_snow_unload, trac_snocan, trac_qflx_snow_unload)
     !
     ! !DESCRIPTION:
     ! Compute snow unloading for one tracer
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_soilp
     integer, intent(in) :: filter_soilp(:)

     ! For description of arguments, see comments in BulkFlux_SnowUnloading. Here, bulk_*
     ! variables refer to bulk water and trac_* variables to the given water tracer.
     real(r8), intent(in) :: bulk_snocan( bounds%begp: )
     real(r8), intent(in) :: bulk_qflx_snow_unload( bounds%begp: )
     real(r8), intent(in) :: trac_snocan( bounds%begp: )
     real(r8), intent(inout) :: trac_qflx_snow_unload( bounds%begp: )
     !
     ! !LOCAL VARIABLES:

     character(len=*), parameter :: subname = 'TracerFlux_SnowUnloading'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(bulk_snocan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(bulk_qflx_snow_unload, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_snocan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(trac_qflx_snow_unload, 1) == bounds%endp), sourcefile, __LINE__)

     associate( &
          begp => bounds%begp, &
          endp => bounds%endp  &
          )

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk_snocan(begp:endp), &
          bulk_val      = bulk_qflx_snow_unload(begp:endp), &
          tracer_source = trac_snocan(begp:endp), &
          tracer_val    = trac_qflx_snow_unload(begp:endp))

     end associate

   end subroutine TracerFlux_SnowUnloading

   !-----------------------------------------------------------------------
   subroutine UpdateState_RemoveSnowUnloading(bounds, num_soilp, filter_soilp, dtime, &
        qflx_snow_unload, snocan)
     !
     ! !DESCRIPTION:
     ! Update snocan based on snow unloading, for bulk or one tracer
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_soilp
     integer, intent(in) :: filter_soilp(:)
     real(r8), intent(in) :: dtime  ! land model time step (sec)

     real(r8) , intent(in)    :: qflx_snow_unload( bounds%begp: ) ! total canopy snow unloading (mm H2O/s)
     real(r8) , intent(inout) :: snocan( bounds%begp: )           ! canopy snow water (mm H2O)
     !
     ! !LOCAL VARIABLES:
     integer :: fp, p

     character(len=*), parameter :: subname = 'UpdateState_RemoveSnowUnloading'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(qflx_snow_unload, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(snocan, 1) == bounds%endp), sourcefile, __LINE__)

     do fp = 1, num_soilp
        p = filter_soilp(fp)

        ! NOTE(wjs, 2019-05-09) The following check sets snocan to exactly 0 if it would
        ! be set to close to 0 (but possibly not exactly 0, because (snocan/dtime)*dtime
        ! is not always exactly equal to snocan). The use of exact equality of these
        ! floating point values in the conditional is probably not desirable long-term. It
        ! currently works for bulk water (because qflx_snow_unload is set to exactly
        ! snocan/dtime if it would exceed this value), but does not necessarily work for
        ! tracers. However, I'm using this exact equality check for now in order to get
        ! bit-for-bit answers during the big refactor of CanopyHydrology.
        if (qflx_snow_unload(p) == snocan(p)/dtime) then
           snocan(p) = 0._r8
        else
           snocan(p) = snocan(p) - qflx_snow_unload(p) * dtime
        end if
     end do

   end subroutine UpdateState_RemoveSnowUnloading

   !-----------------------------------------------------------------------
   subroutine SumFlux_FluxesOntoGround(bounds, &
        num_nolakep, filter_nolakep, &
        num_nolakec, filter_nolakec, &
        qflx_through_snow, qflx_snocanfall, qflx_snow_unload, &
        qflx_through_liq, qflx_liqcanfall, qflx_irrig_drip, &
        qflx_snow_grnd_col, qflx_liq_grnd_col)
     !
     ! !DESCRIPTION:
     ! Compute summed fluxes onto ground, for bulk or one tracer
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_nolakep
     integer, intent(in) :: filter_nolakep(:)
     integer, intent(in) :: num_nolakec
     integer, intent(in) :: filter_nolakec(:)

     real(r8) , intent(in)    :: qflx_through_snow( bounds%begp: ) ! canopy throughfall of snow (mm H2O/s)
     real(r8) , intent(in)    :: qflx_snocanfall( bounds%begp: )   ! rate of excess canopy snow falling off canopy (mm H2O/s)
     real(r8) , intent(in)    :: qflx_snow_unload( bounds%begp: )  ! total canopy snow unloading (mm H2O/s)
     real(r8) , intent(in)    :: qflx_through_liq( bounds%begp: )  ! canopy throughfall of liquid (mm H2O/s)
     real(r8) , intent(in)    :: qflx_liqcanfall( bounds%begp: )   ! rate of excess canopy liquid falling off canopy (mm H2O/s)
     real(r8) , intent(in)    :: qflx_irrig_drip( bounds%begp: )   ! drip irrigation amount (mm H2O/s)

     real(r8) , intent(inout) :: qflx_snow_grnd_col( bounds%begc: )    ! snow on ground after interception (mm H2O/s)
     real(r8) , intent(inout) :: qflx_liq_grnd_col( bounds%begc: )     ! liquid on ground after interception (mm H2O/s)
     !
     ! !LOCAL VARIABLES:
     integer :: fp, p
     real(r8) :: qflx_snow_grnd_patch(bounds%begp:bounds%endp)  ! snow on ground after interception, patch-level (mm H2O/s)
     real(r8) :: qflx_liq_grnd_patch(bounds%begp:bounds%endp)   ! liquid on ground after interception, patch-level (mm H2O/s)

     character(len=*), parameter :: subname = 'SumFlux_FluxesOntoGround'
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(qflx_through_snow, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_snocanfall, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_snow_unload, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_through_liq, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_liqcanfall, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_irrig_drip, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_snow_grnd_col, 1) == bounds%endc), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(qflx_liq_grnd_col, 1) == bounds%endc), sourcefile, __LINE__)

     associate( &
          begp => bounds%begp, &
          endp => bounds%endp, &
          begc => bounds%begc, &
          endc => bounds%endc &
          )

     do fp = 1, num_nolakep
        p = filter_nolakep(fp)

        qflx_snow_grnd_patch(p) = &
             qflx_through_snow(p) + &
             qflx_snocanfall(p) + &
             qflx_snow_unload(p)

        qflx_liq_grnd_patch(p) = &
             qflx_through_liq(p) + &
             qflx_liqcanfall(p) + &
             qflx_irrig_drip(p)
     end do

     call p2c(bounds, num_nolakec, filter_nolakec, &
          qflx_snow_grnd_patch(begp:endp), &
          qflx_snow_grnd_col(begc:endc))

     call p2c(bounds, num_nolakec, filter_nolakec, &
          qflx_liq_grnd_patch(begp:endp), &
          qflx_liq_grnd_col(begc:endc))

     end associate

   end subroutine SumFlux_FluxesOntoGround


   !-----------------------------------------------------------------------
   subroutine BulkDiag_FracWet(bounds, num_soilp, filter_soilp, &
        frac_veg_nosno, elai, esai, snocan, liqcan, &
        fwet, fdry, fcansno)
     !
     ! !DESCRIPTION:
     ! Determine fraction of vegetated surfaces which are wet and
     ! fraction of elai which is dry. The variable ``fwet'' is the
     ! fraction of all vegetation surfaces which are wet including
     ! stem area which contribute to evaporation. The variable ``fdry''
     ! is the fraction of elai which is dry because only leaves
     ! can transpire.  Adjusted for stem area which does not transpire.
     !
     ! ! USES:
     use clm_varcon         , only : tfrz
     ! !ARGUMENTS:
     type(bounds_type), intent(in) :: bounds
     integer, intent(in) :: num_soilp
     integer, intent(in) :: filter_soilp(:)

     integer  , intent(in)    :: frac_veg_nosno( bounds%begp: ) ! fraction of vegetation not covered by snow (0 OR 1)
     real(r8) , intent(in)    :: elai( bounds%begp: )           ! canopy one-sided leaf area index with burying by snow
     real(r8) , intent(in)    :: esai( bounds%begp: )           ! canopy one-sided stem area index with burying by snow
     real(r8) , intent(in)    :: snocan( bounds%begp: )         ! canopy snow water (mm H2O)
     real(r8) , intent(in)    :: liqcan( bounds%begp: )         ! canopy liquid water (mm H2O)
     real(r8) , intent(inout) :: fwet( bounds%begp: )           ! fraction of canopy that is wet (0 to 1)
     real(r8) , intent(inout) :: fdry( bounds%begp: )           ! fraction of foliage that is green and dry [-]
     real(r8) , intent(inout) :: fcansno( bounds%begp: )        ! fraction of canopy that is snow covered (0 to 1)
     !
     ! !LOCAL VARIABLES:
     integer  :: fp,p             ! indices
     real(r8) :: h2ocan           ! total canopy water (mm H2O)
     real(r8) :: vegt             ! lsai
     !-----------------------------------------------------------------------

     SHR_ASSERT_FL((ubound(frac_veg_nosno, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(elai, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(esai, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(snocan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(liqcan, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(fwet, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(fdry, 1) == bounds%endp), sourcefile, __LINE__)
     SHR_ASSERT_FL((ubound(fcansno, 1) == bounds%endp), sourcefile, __LINE__)

     do fp = 1, num_soilp
        p = filter_soilp(fp)
        if (frac_veg_nosno(p) == 1) then
           h2ocan = snocan(p) + liqcan(p)

           if (h2ocan > 0._r8) then
              vegt    = frac_veg_nosno(p)*(elai(p) + esai(p))
              fwet(p) = (h2ocan / (vegt * params_inst%liq_canopy_storage_scalar))**0.666666666666_r8
              fwet(p) = min (fwet(p),maximum_leaf_wetted_fraction)   ! Check for maximum limit of fwet
              if (snocan(p) > 0._r8) then
                 fcansno(p) = (snocan(p) / (vegt * params_inst%snow_canopy_storage_scalar))**0.15_r8 ! must match snocanmx 
                 fcansno(p) = min (fcansno(p),1.0_r8)
              else
                 fcansno(p) = 0._r8
              end if
           else
              fwet(p) = 0._r8
              fcansno(p) = 0._r8
           end if
           fdry(p) = (1._r8-fwet(p))*elai(p)/(elai(p)+esai(p))
        else
           fwet(p) = 0._r8
           fdry(p) = 0._r8
        end if
     end do

   end subroutine BulkDiag_FracWet

end module CanopyHydrologyMod
