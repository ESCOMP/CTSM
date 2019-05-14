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
  use ColumnType      , only : col                
  use PatchType       , only : patch, patch_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyHydrology_readnl ! Read namelist
  public :: CanopyInterceptionAndThroughfall
  public :: CanopyHydrology        ! Run
  public :: readParams

  type, private :: params_type
     real(r8) :: zlnd  ! Roughness length for soil (m)
     real(r8) :: dewmx  ! Canopy maximum storage of liquid water (kg/m2)
     real(r8) :: sno_stor_max  ! Canopy maximum storage of snow (kg/m2)
     real(r8) :: accum_factor  ! Accumulation constant for fractional snow covered area (unitless)
  end type params_type
  type(params_type), private ::  params_inst
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SumFlux_TopOfCanopyInputs
  private :: BulkFlux_CanopyInterceptionAndThroughfall
  private :: BulkDiag_FracWet    ! Determine fraction of vegetated surface that is wet
  private :: FracH2oSfc ! Determine fraction of land surfaces which are submerged  
  !
  ! !PRIVATE DATA MEMBERS:
  integer :: oldfflag=0  ! use old fsno parameterization (N&Y07) 
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
         oldfflag, &
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
    call shr_mpi_bcast(oldfflag, mpicom)
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

    ! Roughness length for soil (m)
    call readNcdioScalar(ncid, 'zlnd', subname, params_inst%zlnd)
    ! Canopy maximum storage of liquid water (kg/m2)
    call readNcdioScalar(ncid, 'dewmx', subname, params_inst%dewmx)
    ! Canopy maximum storage of snow (kg/m2)
    call readNcdioScalar(ncid, 'sno_stor_max', subname, params_inst%sno_stor_max)
    ! Accumulation constant for fractional snow covered area (unitless)
    call readNcdioScalar(ncid, 'accum_factor', subname, params_inst%accum_factor)

   end subroutine readParams

   !-----------------------------------------------------------------------
   subroutine CanopyInterceptionAndThroughfall(bounds, &
        num_soilp, filter_soilp, &
        num_nolakep, filter_nolakep, &
        patch, canopystate_inst, atm2lnd_inst, water_inst)
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
     type(patch_type)       , intent(in)    :: patch
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
        call SumFlux_FluxesOntoGround(bounds, num_nolakep, filter_nolakep, &
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
     integer :: fp, p
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
           qflx_through_snow(p) = forc_snow(p)
           qflx_through_liq(p)  = qflx_liq_above_canopy(p)
           qflx_intercepted_snow(p) = 0._r8
           qflx_intercepted_liq(p) = 0._r8
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
           liqcanmx = params_inst%dewmx * (elai(p) + esai(p))
           qflx_liqcanfall(p) = max((liqcan(p) - liqcanmx)/dtime, 0._r8)
           snocanmx = params_inst%sno_stor_max * (elai(p) + esai(p))
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

        ! FIXME(wjs, 2019-05-09) Put in place adjustments to get bit-for-bit
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
   subroutine SumFlux_FluxesOntoGround(bounds, num_nolakep, filter_nolakep, &
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

        ! FIXME(wjs, 2019-05-10) To get bit-for-bit, probably need to replace
        ! qflx_snow_unload with (qflx_snow_unload*dtime)/dtime. Probably no need to check
        ! that that is only a roundoff-level change (i.e., no need to confirm that
        ! (qflx_snow_unload*dtime)/dtime equals qflx_snow_unload within roundoff.
        qflx_snow_grnd_patch(p) = &
             qflx_through_snow(p) + &
             qflx_snocanfall(p) + &
             qflx_snow_unload(p)

        ! FIXME(wjs, 2019-05-14) Parenthesize the addition of the first two terms to get
        ! bit-for-bit. (Then remove this set of parentheses when I can tolerate
        ! roundoff-level changes.)
        qflx_liq_grnd_patch(p) = &
             qflx_through_liq(p) + &
             qflx_liqcanfall(p) + &
             qflx_irrig_drip(p)
     end do

     call p2c(bounds, num_nolakep, filter_nolakep, &
          qflx_snow_grnd_patch(begp:endp), &
          qflx_snow_grnd_col(begc:endc))

     call p2c(bounds, num_nolakep, filter_nolakep, &
          qflx_liq_grnd_patch(begp:endp), &
          qflx_liq_grnd_col(begc:endc))

     end associate

   end subroutine SumFlux_FluxesOntoGround


   !-----------------------------------------------------------------------
   subroutine CanopyHydrology(bounds, &
        num_nolakec, filter_nolakec, &
        atm2lnd_inst, temperature_inst, &
        aerosol_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, &
        waterfluxbulk_inst)
     !
     ! !DESCRIPTION:
     ! Calculation of snow layer initialization if the snow accumulation exceeds 10 mm.
     ! Note:  The evaporation loss is taken off after the calculation of leaf
     ! temperature in the subroutine clm\_leaftem.f90, not in this subroutine.
     !
     ! !USES:
     use clm_varcon         , only : hfus, denice, rpi, spval, tfrz, int_snow_max
     use column_varcon      , only : icol_roof, icol_sunwall, icol_shadewall
     use landunit_varcon    , only : istcrop, istwet, istsoil, istice_mec 
     use clm_varctl         , only : subgridflag
     use clm_varpar         , only : nlevsoi,nlevsno
     use SnowHydrologyMod   , only : NewSnowBulkDensity
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds     
     integer                , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
     integer                , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
     type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
     type(temperature_type) , intent(inout) :: temperature_inst
     type(aerosol_type)     , intent(inout) :: aerosol_inst
     type(waterstatebulk_type)      , intent(inout) :: waterstatebulk_inst
     type(waterdiagnosticbulk_type) , intent(inout) :: waterdiagnosticbulk_inst
     type(waterfluxbulk_type)       , intent(inout) :: waterfluxbulk_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: f                                            ! filter index
     integer  :: c                                            ! column index
     integer  :: l                                            ! landunit index
     integer  :: g                                            ! gridcell index
     integer  :: newnode                                      ! flag when new snow node is set, (1=yes, 0=no)
     real(r8) :: dtime                                        ! land model time step (sec)
     real(r8) :: dz_snowf                                     ! layer thickness rate change due to precipitation [mm/s]
     real(r8) :: bifall(bounds%begc:bounds%endc)              ! bulk density of newly fallen dry snow [kg/m3]
     real(r8) :: z_avg                                        ! grid cell average snow depth
     real(r8) :: rho_avg                                      ! avg density of snow column
     real(r8) :: temp_snow_depth,temp_intsnow                 ! temporary variables
     real(r8) :: fmelt
     real(r8) :: smr
     real(r8) :: delf_melt
     real(r8) :: fsno_new
     real(r8) :: int_snow_limited ! integrated snowfall, limited to be no greater than int_snow_max [mm]
     real(r8) :: newsnow(bounds%begc:bounds%endc)
     real(r8) :: snowmelt(bounds%begc:bounds%endc)
     integer  :: j
     !-----------------------------------------------------------------------

     associate(                                                             & 
          snl                  => col%snl                                 , & ! Input:  [integer  (:)   ]  number of snow layers                    
          n_melt               => col%n_melt                              , & ! Input:  [real(r8) (:)   ]  SCA shape parameter                     
          zi                   => col%zi                                  , & ! Output: [real(r8) (:,:) ]  interface level below a "z" level (m) 
          dz                   => col%dz                                  , & ! Output: [real(r8) (:,:) ]  layer depth (m)                       
          z                    => col%z                                   , & ! Output: [real(r8) (:,:) ]  layer thickness (m)                   

          forc_t               => atm2lnd_inst%forc_t_downscaled_col      , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)        

          t_grnd               => temperature_inst%t_grnd_col             , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)             
          t_soisno             => temperature_inst%t_soisno_col           , & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)  

          h2osno               => waterstatebulk_inst%h2osno_col              , & ! Output: [real(r8) (:)   ]  snow water (mm H2O)                     
          snow_depth           => waterdiagnosticbulk_inst%snow_depth_col          , & ! Output: [real(r8) (:)   ]  snow height (m)                         
          int_snow             => waterstatebulk_inst%int_snow_col            , & ! Output: [real(r8) (:)   ]  integrated snowfall [mm]                
          frac_sno_eff         => waterdiagnosticbulk_inst%frac_sno_eff_col        , & ! Output: [real(r8) (:)   ]  eff. fraction of ground covered by snow (0 to 1)
          frac_sno             => waterdiagnosticbulk_inst%frac_sno_col            , & ! Output: [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          frac_iceold          => waterdiagnosticbulk_inst%frac_iceold_col         , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water
          h2osoi_ice           => waterstatebulk_inst%h2osoi_ice_col          , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                      
          h2osoi_liq           => waterstatebulk_inst%h2osoi_liq_col          , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                  
          swe_old              => waterdiagnosticbulk_inst%swe_old_col             , & ! Output: [real(r8) (:,:) ]  snow water before update              

          qflx_snow_drain       => waterfluxbulk_inst%qflx_snow_drain_col     , & ! Input: [real(r8) (:)   ]  drainage from snow pack from previous time step       
          qflx_snow_h2osfc     => waterfluxbulk_inst%qflx_snow_h2osfc_col     , & ! Output: [real(r8) (:)   ]  snow falling on surface water (mm/s)     
          qflx_snow_grnd_col   => waterfluxbulk_inst%qflx_snow_grnd_col         & ! Input:  [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
          )

       ! Compute time step
       dtime = get_step_size()

       ! Determine snow height and snow water

       call NewSnowBulkDensity(bounds, num_nolakec, filter_nolakec, &
            atm2lnd_inst, bifall(bounds%begc:bounds%endc))

       do f = 1, num_nolakec
          c = filter_nolakec(f)
          l = col%landunit(c)
          g = col%gridcell(c)

          ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
          ! U.S.Department of Agriculture Forest Service, Project F,
          ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

          qflx_snow_h2osfc(c) = 0._r8
          ! set temporary variables prior to updating
          temp_snow_depth=snow_depth(c)
          ! save initial snow content
          do j= -nlevsno+1,snl(c)
             swe_old(c,j) = 0.0_r8
          end do
          do j= snl(c)+1,0
             swe_old(c,j)=h2osoi_liq(c,j)+h2osoi_ice(c,j)
          enddo

          ! all snow falls on ground, no snow on h2osfc
          newsnow(c) = qflx_snow_grnd_col(c) * dtime

          ! update int_snow
          int_snow(c) = max(int_snow(c),h2osno(c)) !h2osno could be larger due to frost

          ! snowmelt from previous time step * dtime
          snowmelt(c) = qflx_snow_drain(c) * dtime

          if (h2osno(c) > 0.0) then

             !======================  FSCA PARAMETERIZATIONS  ======================
             ! fsca parameterization based on *changes* in swe
             ! first compute change from melt during previous time step
             if(snowmelt(c) > 0._r8) then

                int_snow_limited = min(int_snow(c), int_snow_max)
                smr=min(1._r8,h2osno(c)/int_snow_limited)

                frac_sno(c) = 1. - (acos(min(1._r8,(2.*smr - 1._r8)))/rpi)**(n_melt(c))

             endif

             ! update fsca by new snow event, add to previous fsca
             if (newsnow(c) > 0._r8) then
                fsno_new = 1._r8 - (1._r8 - tanh(params_inst%accum_factor * newsnow(c))) * (1._r8 - frac_sno(c))
                frac_sno(c) = fsno_new

                ! reset int_snow after accumulation events
                temp_intsnow= (h2osno(c) + newsnow(c)) &
                     / (0.5*(cos(rpi*(1._r8-max(frac_sno(c),1e-6_r8))**(1./n_melt(c)))+1._r8))
                int_snow(c) = min(1.e8_r8,temp_intsnow)
             endif

             !====================================================================

             ! for subgrid fluxes
             if (subgridflag ==1 .and. .not. lun%urbpoi(l)) then
                if (frac_sno(c) > 0._r8)then
                   snow_depth(c)=snow_depth(c) + newsnow(c)/(bifall(c) * frac_sno(c))
                else
                   snow_depth(c)=0._r8
                end if
             else
                ! for uniform snow cover
                snow_depth(c)=snow_depth(c)+newsnow(c)/bifall(c)
             endif

             ! use original fsca formulation (n&y 07)
             if (oldfflag == 1) then 
                ! snow cover fraction in Niu et al. 2007
                if(snow_depth(c) > 0.0_r8)  then
                   frac_sno(c) = tanh(snow_depth(c) / (2.5_r8 * params_inst%zlnd * &
                        (min(800._r8,(h2osno(c)+ newsnow(c))/snow_depth(c))/100._r8)**1._r8) )
                endif
                if(h2osno(c) < 1.0_r8)  then
                   frac_sno(c)=min(frac_sno(c),h2osno(c))
                endif
             endif

          else !h2osno == 0
             ! initialize frac_sno and snow_depth when no snow present initially
             if (newsnow(c) > 0._r8) then 
                z_avg = newsnow(c)/bifall(c)
                fmelt=newsnow(c)
                frac_sno(c) = tanh(params_inst%accum_factor * newsnow(c))

                ! make int_snow consistent w/ new fsno, h2osno
                int_snow(c) = 0. !reset prior to adding newsnow below
                temp_intsnow= (h2osno(c) + newsnow(c)) &
                     / (0.5*(cos(rpi*(1._r8-max(frac_sno(c),1e-6_r8))**(1./n_melt(c)))+1._r8))
                int_snow(c) = min(1.e8_r8,temp_intsnow)

                ! update snow_depth and h2osno to be consistent with frac_sno, z_avg
                if (subgridflag ==1 .and. .not. lun%urbpoi(l)) then
                   snow_depth(c)=z_avg/frac_sno(c)
                else
                   snow_depth(c)=newsnow(c)/bifall(c)
                endif
                ! use n&y07 formulation
                if (oldfflag == 1) then 
                   ! snow cover fraction in Niu et al. 2007
                   if(snow_depth(c) > 0.0_r8)  then
                      frac_sno(c) = tanh(snow_depth(c) / (2.5_r8 * params_inst%zlnd * &
                           (min(800._r8,newsnow(c)/snow_depth(c))/100._r8)**1._r8) )
                   endif
                endif
             else
                z_avg = 0._r8
                snow_depth(c) = 0._r8
                frac_sno(c) = 0._r8
             endif
          endif ! end of h2osno > 0

          ! no snow on surface water
          qflx_snow_h2osfc(c) = 0._r8

          ! update h2osno for new snow
          h2osno(c) = h2osno(c) + newsnow(c) 
          int_snow(c) = int_snow(c) + newsnow(c)

          ! update change in snow depth
          dz_snowf = (snow_depth(c) - temp_snow_depth) / dtime

          ! set frac_sno_eff variable
          if (.not. lun%urbpoi(l)) then
             if (subgridflag ==1) then 
                frac_sno_eff(c) = frac_sno(c)
             else
                frac_sno_eff(c) = 1._r8
             endif
          else
             frac_sno_eff(c) = 1._r8
          endif

          if (lun%itype(l)==istwet .and. t_grnd(c)>tfrz) then
             h2osno(c)=0._r8
             snow_depth(c)=0._r8
          end if

          ! When the snow accumulation exceeds 10 mm, initialize snow layer
          ! Currently, the water temperature for the precipitation is simply set
          ! as the surface air temperature

          newnode = 0    ! flag for when snow node will be initialized
          if (snl(c) == 0 .and. frac_sno(c)*snow_depth(c) >= 0.01_r8) then
             newnode = 1
             snl(c) = -1
             dz(c,0) = snow_depth(c)                       ! meter
             z(c,0) = -0.5_r8*dz(c,0)
             zi(c,-1) = -dz(c,0)
             t_soisno(c,0) = min(tfrz, forc_t(c))      ! K
             h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
             h2osoi_liq(c,0) = 0._r8                   ! kg/m2
             frac_iceold(c,0) = 1._r8

             ! intitialize SNICAR variables for fresh snow:
             call aerosol_inst%Reset(column=c)
             call waterdiagnosticbulk_inst%ResetBulk(column=c)
          end if

          ! The change of ice partial density of surface node due to precipitation.
          ! Only ice part of snowfall is added here, the liquid part will be added
          ! later.

          if (snl(c) < 0 .and. newnode == 0) then
             h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+newsnow(c)
             dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
          end if

       end do

       ! update surface water fraction (this may modify frac_sno)
       call FracH2oSfc(bounds, num_nolakec, filter_nolakec, &
            waterstatebulk_inst, waterdiagnosticbulk_inst)

     end associate 

   end subroutine CanopyHydrology

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
     real(r8) :: dewmxi           ! inverse of maximum allowed dew [1/mm]
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
              dewmxi  = 1.0_r8/params_inst%dewmx  ! wasteful division
              fwet(p) = ((dewmxi/vegt)*h2ocan)**0.666666666666_r8
              fwet(p) = min (fwet(p),maximum_leaf_wetted_fraction)   ! Check for maximum limit of fwet
              if (snocan(p) > 0._r8) then
                 dewmxi  = 1.0_r8/params_inst%dewmx  ! wasteful division
                 fcansno(p) = ((dewmxi / (vegt * params_inst%sno_stor_max * 10.0_r8)) * snocan(p))**0.15_r8 ! must match snocanmx 
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

   !-----------------------------------------------------------------------
   subroutine FracH2OSfc(bounds, num_h2osfc, filter_h2osfc, &
        waterstatebulk_inst, waterdiagnosticbulk_inst, no_update)
     !
     ! !DESCRIPTION:
     ! Determine fraction of land surfaces which are submerged  
     ! based on surface microtopography and surface water storage.
     !
     ! !USES:
     use shr_const_mod   , only : shr_const_pi
     use shr_spfn_mod    , only : erf => shr_spfn_erf
     use landunit_varcon , only : istsoil, istcrop
     !
     ! !ARGUMENTS:
     type(bounds_type)     , intent(in)           :: bounds           
     integer               , intent(in)           :: num_h2osfc       ! number of column points in column filter
     integer               , intent(in)           :: filter_h2osfc(:) ! column filter 
     type(waterstatebulk_type) , intent(inout)        :: waterstatebulk_inst
     type(waterdiagnosticbulk_type) , intent(inout)        :: waterdiagnosticbulk_inst
     integer               , intent(in), optional :: no_update        ! flag to make calculation w/o updating variables
     !
     ! !LOCAL VARIABLES:
     integer :: c,f,l          ! indices
     real(r8):: d,fd,dfdd      ! temporary variable for frac_h2o iteration
     real(r8):: sigma          ! microtopography pdf sigma in mm
     real(r8):: min_h2osfc
     !-----------------------------------------------------------------------

     associate(                                              & 
          micro_sigma  => col%micro_sigma                  , & ! Input:  [real(r8) (:)   ] microtopography pdf sigma (m)                     

          h2osno       => waterstatebulk_inst%h2osno_col       , & ! Input:  [real(r8) (:)   ] snow water (mm H2O)                               
          
          h2osoi_liq   => waterstatebulk_inst%h2osoi_liq_col   , & ! Output: [real(r8) (:,:) ] liquid water (col,lyr) [kg/m2]                  
          h2osfc       => waterstatebulk_inst%h2osfc_col       , & ! Output: [real(r8) (:)   ] surface water (mm)                                
          frac_sno     => waterdiagnosticbulk_inst%frac_sno_col     , & ! Output: [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)       
          frac_sno_eff => waterdiagnosticbulk_inst%frac_sno_eff_col , & ! Output: [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)  
          frac_h2osfc  => waterdiagnosticbulk_inst%frac_h2osfc_col,   & ! Output: [real(r8) (:)   ] col fractional area with surface water greater than zero 
          frac_h2osfc_nosnow  => waterdiagnosticbulk_inst%frac_h2osfc_nosnow_col    & ! Output: [real(r8) (:)   ] col fractional area with surface water greater than zero (if no snow present)
          )

       ! arbitrary lower limit on h2osfc for safer numerics...
       min_h2osfc=1.e-8_r8

       do f = 1, num_h2osfc
          c = filter_h2osfc(f)
          l = col%landunit(c)

          ! h2osfc only calculated for soil vegetated land units
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

             !  Use newton-raphson method to iteratively determine frac_h2osfc
             !  based on amount of surface water storage (h2osfc) and 
             !  microtopography variability (micro_sigma)

             if (h2osfc(c) > min_h2osfc) then
                ! a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
                d=0.0

                sigma=1.0e3 * micro_sigma(c) ! convert to mm
                do l=1,10
                   fd = 0.5*d*(1.0_r8+erf(d/(sigma*sqrt(2.0)))) &
                        +sigma/sqrt(2.0*shr_const_pi)*exp(-d**2/(2.0*sigma**2)) &
                        -h2osfc(c)
                   dfdd = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

                   d = d - fd/dfdd
                enddo
                !--  update the submerged areal fraction using the new d value
                frac_h2osfc(c) = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

             else
                frac_h2osfc(c) = 0._r8
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + h2osfc(c)
                h2osfc(c)=0._r8
             endif

             frac_h2osfc_nosnow(c) = frac_h2osfc(c)


             if (.not. present(no_update)) then

                ! adjust fh2o, fsno when sum is greater than zero
                if (frac_sno(c) > (1._r8 - frac_h2osfc(c)) .and. h2osno(c) > 0) then

                   if (frac_h2osfc(c) > 0.01_r8) then             
                      frac_h2osfc(c) = max(1.0_r8 - frac_sno(c),0.01_r8)
                      frac_sno(c) = 1.0_r8 - frac_h2osfc(c)
                   else
                      frac_sno(c) = 1.0_r8 - frac_h2osfc(c)
                   endif
                   frac_sno_eff(c)=frac_sno(c)

                endif

             endif ! end of no_update construct

          else !if landunit not istsoil/istcrop, set frac_h2osfc to zero

             frac_h2osfc(c) = 0._r8

          endif

       end do

     end associate

   end subroutine FracH2OSfc

end module CanopyHydrologyMod
