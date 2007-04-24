#include <misc.h>
#include <preproc.h>

module driver

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: driver
!
! !DESCRIPTION:
! This module provides the main CLM driver calling sequence.  Most
! computations occurs over ``clumps'' of gridcells (and associated subgrid
! scale entities) assigned to each MPI process. Computation is further
! parallelized by looping over clumps on each process using shared memory
! OpenMP or Cray Streaming Directives.
!
! The main CLM driver calling sequence is as follows:
! \begin{verbatim}
!
! * Recv data from flux coupler [COUP_CSM] 
!
! + interpMonthlyVeg      interpolate monthly vegetation data [!DGVM]
!   + readMonthlyVegetation read vegetation data for two months [!DGVM]
!
! ==== Begin Loop over clumps ====
!  -> DriverInit          save of variables from previous time step
!  -> Hydrology1          canopy interception and precip on ground
!     -> FracWet          fraction of wet vegetated surface and dry elai
!  -> SurfaceRadiation    surface solar radiation
!  -> Biogeophysics1      leaf temperature and surface fluxes
!  -> BareGroundFluxes    surface fluxes for bare soil or snow-covered
!                         vegetation patches
!     -> MoninObukIni     first-guess Monin-Obukhov length and wind speed
!     -> FrictionVelocity friction velocity and potential temperature and
!                         humidity profiles
!  -> CanopyFluxes        leaf temperature and surface fluxes for vegetated
!                         patches
!     -> QSat             saturated vapor pressure, specific humidity, &
!                         derivatives at leaf surface
!     -> MoninObukIni     first-guess Monin-Obukhov length and wind speed
!     -> FrictionVelocity friction velocity and potential temperature and
!                         humidity profiles
!     -> Stomata          stomatal resistance and photosynthesis for
!                         sunlit leaves
!     -> Stomata          stomatal resistance and photosynthesis for
!                         shaded leaves
!     -> QSat             recalculation of saturated vapor pressure,
!                         specific humidity, & derivatives at leaf surface
!  -> Biogeophysics_Lake  lake temperature and surface fluxes
!   + VOCEmission         compute VOC emission [VOC]
!   + DGVMRespiration     CO2 respriation and plant production [DGVM]
!   + DGVMEcosystemDyn    DGVM ecosystem dynamics: vegetation phenology [!DGVM]
!  -> EcosystemDyn        "static" ecosystem dynamics: vegetation phenology
!                         and soil carbon [!DGVM]
!  -> Biogeophysics2      soil/snow & ground temp and update surface fluxes
!  -> pft2col             Average from PFT level to column level
!  -> Hydrology2          surface and soil hydrology
!  -> Hydrology_Lake      lake hydrology
!  -> SnowAge             update snow age for surface albedo calcualtion
!  -> BalanceCheck        check for errors in energy and water balances
!  -> SurfaceAlbedo       albedos for next time step
!  ====  End Loop over clumps  ====
!
! * Average fluxes over time interval and send to flux coupler [COUP_CSM]
!  -> write_diagnostic    output diagnostic if appropriate
!   + Rtmriverflux        calls RTM river routing model [RTM]
!  -> updateAccFlds       update accumulated fields
!  -> update_hbuf         accumulate history fields for time interval
!
!  Begin DGVM calculations at end of model year [DGVM]
!    ==== Begin Loop over clumps ====
!     + lpj                 run LPJ ecosystem dynamics: reproduction, turnover,
!                           kill, allocation, light, mortality, fire
!     + lpjreset            reset variables & initialize for next year
!     + resetWeightsDGVM    reset variables and patch weights
!     + resetTimeConstDGVM  reset time constant variables for new pfts
!    ====  End Loop over clumps  ====
!  End DGVM calculations at end of model year [DGVM]
!
!  -> htapes_wrapup       write history tapes if appropriate
!  -> DGVMhist            write DGVM history file
!  -> restFile            write restart file if appropriate
!  -> inicfile            write initial file if appropriate
! \end{verbatim}
! Optional subroutines are denoted by an plus (+) with the associated
! CPP variable in brackets at the end of the line.  Coupler communication
! when coupled with CCSM components is denoted by an asterisk (*).
!
! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clmtype
  use clm_varctl          , only : wrtdia, fpftdyn, fndepdyn
  use spmdMod             , only : masterproc
  use decompMod           , only : get_proc_clumps, get_clump_bounds
  use filterMod           , only : filter, setFilters
  use pftdynMod           , only : pftdyn_interp, pftdyn_wbal_init, pftdyn_wbal 
  use clm_varcon          , only : zlnd
  use clm_time_manager        , only : get_step_size, get_curr_calday, &
                                   get_curr_date, get_ref_date, get_nstep, is_perpetual
  use histFileMod         , only : update_hbuf, htapes_wrapup
  use restFileMod         , only : restFile_write, restFile_write_binary, restFile_filename
  use inicFileMod         , only : inicfile_perp  
  use accFldsMod          , only : updateAccFlds
  use DriverInitMod       , only : DriverInit
  use BalanceCheckMod     , only : BeginWaterBalance, BalanceCheck
  use SurfaceRadiationMod , only : SurfaceRadiation
  use Hydrology1Mod       , only : Hydrology1
  use Hydrology2Mod       , only : Hydrology2
  use HydrologyLakeMod    , only : HydrologyLake
  use Biogeophysics1Mod   , only : Biogeophysics1
  use BareGroundFluxesMod , only : BareGroundFluxes
  use CanopyFluxesMod     , only : CanopyFluxes
  use Biogeophysics2Mod   , only : Biogeophysics2
  use BiogeophysicsLakeMod, only : BiogeophysicsLake
  use SurfaceAlbedoMod    , only : SurfaceAlbedo, Snowage
  use pft2colMod          , only : pft2col
#if (defined DGVM)
  use DGVMEcosystemDynMod , only : DGVMEcosystemDyn, DGVMRespiration
  use DGVMMod             , only : lpj, lpjreset, histDGVM, &
	                           resetweightsdgvm, resettimeconstdgvm 
#elif (defined CN)
  use CNEcosystemDynMod   , only : CNEcosystemDyn
  use CNBalanceCheckMod   , only : BeginCBalance, BeginNBalance, &
                                   CBalanceCheck, NBalanceCheck
  use ndepFileMod         , only : ndepdyn_interp
#else
  use STATICEcosysDynMod  , only : EcosystemDyn, interpMonthlyVeg
#endif
#if (defined DUST)
  use DUSTMod             , only : DustDryDep, DustEmission
#endif
#if (defined VOC)
  use VOCEmissionMod      , only : VOCEmission
#endif
#if (defined CASA)
  use CASAPhenologyMod    , only : CASAPhenology
  use CASAMod             , only : Casa
#endif
#if (defined RTM)
  use RtmMod              , only : Rtmriverflux
#endif
  use abortutils          , only : endrun
  use perf_mod

!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: driver1
  public :: driver2
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: write_diagnostic
  private :: do_restwrite
  private :: do_inicwrite
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!
! !ROUTINE: driver1
!
! !INTERFACE:
subroutine driver1 (doalb, caldayp1, declinp1)
!
! !ARGUMENTS:
  implicit none
  logical , intent(in) :: doalb    ! true if time for surface albedo calc
  real(r8), intent(in) :: caldayp1 ! calendar day for nstep+1
  real(r8), intent(in) :: declinp1 ! declination angle for next time step
!
! !CALLED FROM:
! program program_off (if COUP_OFFLINE cpp variable is defined)
! program program_csm (if COUP_CSM cpp variable is defined)
! subroutine atm_lnddrv in module atm_lndMod (if SEQ_MCT or SEQ_ESMF cpp variable
!   is defined)
!
! !REVISION HISTORY:
! 2002.10.01  Mariana Vertenstein latest update to new data structures
! 11/26/03, Peter Thornton: Added new call for SurfaceRadiationSunShade when
!  cpp directive SUNSHA is set, for sunlit/shaded canopy radiation.
! 4/25/05, Peter Thornton: Made the sun/shade routine the default, no longer
!  need to have SUNSHA defined.  
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: nc, c         ! indices
  integer  :: nclumps       ! number of clumps on this processor
  integer  :: nstep         ! time step number
  integer  :: begp, endp    ! clump beginning and ending pft indices
  integer  :: begc, endc    ! clump beginning and ending column indices
  integer  :: begl, endl    ! clump beginning and ending landunit indices
  integer  :: begg, endg    ! clump beginning and ending gridcell indices
  type(column_type)  , pointer :: cptr    ! pointer to column derived subtype
!-----------------------------------------------------------------------
  ! Set pointers into derived type

  cptr => clm3%g%l%c

#if (defined DGVM)
   ! no call
#elif (defined CN)
   ! no call
#else
  ! ============================================================================
  ! Determine weights for time interpolation of monthly vegetation data.
  ! This also determines whether it is time to read new monthly vegetation and
  ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
  ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
  ! weights obtained here are used in subroutine ecosystemdyn to obtain time
  ! interpolated values.
  ! ============================================================================

  if (doalb) call interpMonthlyVeg()
#endif

  ! ============================================================================
  ! Loop over clumps
  ! ============================================================================

  nclumps = get_proc_clumps()

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
#if !defined (USE_OMP)
!!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
#endif
  do nc = 1,nclumps

     ! ============================================================================
     ! Determine clump boundaries
     ! ============================================================================

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! Initialize the mass balance checks: water, carbon, and nitrogen
     ! ============================================================================

     call t_startf('begwbal')
     call BeginWaterBalance(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec)
     call t_stopf('begwbal')

#if (defined CN)
     if (doalb) then
        call t_startf('begcnbal')
        
        call BeginCBalance(filter(nc)%num_soilc,filter(nc)%soilc, &
             filter(nc)%num_soilp, filter(nc)%soilp)
             
        call BeginNBalance(filter(nc)%num_soilc,filter(nc)%soilc, &
             filter(nc)%num_soilp, filter(nc)%soilp)
             
        call t_stopf('begcnbal')
     end if 
#endif

  end do
!$OMP END PARALLEL DO
#if !defined (USE_OMP)
!!CSD$ END PARALLEL DO
#endif

  ! ============================================================================
  ! Initialize h2ocan_loss to zero
  ! ============================================================================
  call t_startf('pftdynwts')
  call pftdyn_wbal_init()

#if (!defined DGVM)
  ! ============================================================================
  ! Update weights and reset filters if dynamic land use
  ! This needs to be done outside the clumps loop, but after BeginWaterBalance()
  ! ============================================================================
  if (doalb .and. fpftdyn /= ' ') then
     call pftdyn_interp()
     call pftdyn_wbal()
     call setFilters()
  end if
#endif

#if (defined CN)
  ! ============================================================================
  ! Update dynamic N deposition field, on albedo timestep
  ! currently being done outside clumps loop, but no reason why it couldn't be
  ! re-written to go inside.
  ! ============================================================================
  if (doalb .and. fndepdyn /= ' ') then
     call ndepdyn_interp()
  end if
#endif       
  call t_stopf('pftdynwts')

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
#if !defined (USE_OMP)
!!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
#endif
  do nc = 1,nclumps

     ! ============================================================================
     ! Determine clump boundaries
     ! ============================================================================

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! Initialize variables from previous time step and
     ! Determine canopy interception and precipitation onto ground surface.
     ! Determine the fraction of foliage covered by water and the fraction
     ! of foliage that is dry and transpiring. Initialize snow layer if the
     ! snow accumulation exceeds 10 mm.
     ! ============================================================================

     call t_startf('drvinit')
     call DriverInit(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec)
     call t_stopf('drvinit')

     ! ============================================================================
     ! Hydrology1
     ! ============================================================================

     call t_startf('hydro1')
     call Hydrology1(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('hydro1')

     ! ============================================================================
     ! Surface Radiation
     ! ============================================================================

     call t_startf('surfrad')
     call SurfaceRadiation(begp, endp)
     call t_stopf('surfrad')

     ! ============================================================================
     ! Determine leaf temperature and surface fluxes based on ground
     ! temperature from previous time step.
     ! ============================================================================

     call t_startf('bgp1')
     call Biogeophysics1(begg, endg, begc, endc, begp, endp, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp1')

     ! ============================================================================
     ! Determine bare soil or snow-covered vegetation surface temperature and fluxes
     ! Calculate Ground fluxes (frac_veg_nosno is either 1 or 0)
     ! ============================================================================

     call t_startf('bgflux')
     call BareGroundFluxes(begp, endp, &
                           filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgflux')

     ! ============================================================================
     ! Determine non snow-covered vegetation surface temperature and fluxes
     ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
     ! and leaf water change by evapotranspiration
     ! ============================================================================

     call t_startf('canflux')
     call CanopyFluxes(begg, endg, begc, endc, begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('canflux')

     ! ============================================================================
     ! Determine lake temperature and surface fluxes
     ! ============================================================================

     call t_startf('bgplake')
     call BiogeophysicsLake(begc, endc, begp, endp, &
                            filter(nc)%num_lakec, filter(nc)%lakec, &
                            filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('bgplake')

     ! ============================================================================
     ! Determine VOC, DUST and DGVM Respiration if appropriate
     ! ============================================================================

     call t_startf('bgc')

#if (defined DUST)
     ! Dust mobilization (C. Zender's modified codes)
     call DustEmission(begp, endp, begc, endc, begl, endl, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)

     ! Dust dry deposition (C. Zender's modified codes)
     call DustDryDep(begp, endp)
#endif

#if (defined VOC)
     ! VOC emission (A. Guenther's model)
     call VOCEmission(begp, endp, &
                      filter(nc)%num_nolakep, filter(nc)%nolakep)
#endif

     call t_stopf('bgc')

     ! ============================================================================
     ! Determine soil/snow temperatures including ground temperature and
     ! update surface fluxes for new ground temperature.
     ! ============================================================================

     call t_startf('bgp2')
     call Biogeophysics2(begc, endc, begp, endp, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp2')

     ! ============================================================================
     ! Perform averaging from PFT level to column level
     ! ============================================================================

     call t_startf('pft2col')
     call pft2col(begc, endc, filter(nc)%num_nolakec, filter(nc)%nolakec)
     call t_stopf('pft2col')

     ! ============================================================================
     ! Vertical (column) soil and surface hydrology
     ! ============================================================================

     call t_startf('hydro2')
     call Hydrology2(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     filter(nc)%num_snowc, filter(nc)%snowc, &
                     filter(nc)%num_nosnowc, filter(nc)%nosnowc)
     call t_stopf('hydro2')

     ! ============================================================================
     ! Lake hydrology
     ! ============================================================================

     call t_startf('hylake')
     call HydrologyLake(begp, endp, &
                        filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('hylake')

     ! ============================================================================
     ! Update Snow Age (needed for surface albedo calculation
     ! ============================================================================

     call t_startf('snowage')
     call SnowAge(begc, endc)

     ! ============================================================================
     ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
     ! ============================================================================

!dir$ concurrent
!cdir nodep
     do c = begc,endc
        cptr%cps%frac_sno(c) = cptr%cps%snowdp(c) / (10._r8*zlnd + cptr%cps%snowdp(c))
     end do
     call t_stopf('snowage')

     ! ============================================================================
     ! Ecosystem dynamics: Uses CN, DGVM, or static parameterizations
     ! ============================================================================

#if (defined CASA)
     call t_startf('casa')
     call CASAPhenology(begp, endp, filter(nc)%num_soilp, filter(nc)%soilp)
     call Casa(begp, endp, filter(nc)%num_soilp, filter(nc)%soilp)
     call t_stopf('casa')
#endif

     call t_startf('ecosysdyn')

#if (defined DGVM)
     ! Prognostic biogeography,
     ! surface biogeochemical fluxes: co2 respiration and plant production
     call DGVMRespiration(begc, endc, begp, endp, &
                          filter(nc)%num_nolakec, filter(nc)%nolakec, &
                          filter(nc)%num_nolakep, filter(nc)%nolakep)

     call DGVMEcosystemDyn(begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep, &
                       doalb, endofyr=.false.)
#elif (defined CN)
     ! Prescribed biogeography,
     ! fully prognostic canopy structure and C-N biogeochemistry
     call CNEcosystemDyn(begc,endc,begp,endp,filter(nc)%num_soilc,&
                  filter(nc)%soilc, filter(nc)%num_soilp, &
                  filter(nc)%soilp, doalb)          
#else
     ! Prescribed biogeography,
     ! prescribed canopy structure, some prognostic carbon fluxes
     call EcosystemDyn(begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep, &
                       doalb)
#endif
     call t_stopf('ecosysdyn')

     ! ============================================================================
     ! Check the energy and water balance, also carbon and nitrogen balance
     ! ============================================================================

     call t_startf('balchk')
     call BalanceCheck(begp, endp, begc, endc)
     call t_stopf('balchk')
     
#if (defined CN)
     if (doalb) then
        call t_startf('cnbalchk')
        
        call CBalanceCheck(filter(nc)%num_soilc,filter(nc)%soilc, &
             filter(nc)%num_soilp, filter(nc)%soilp)
          
        call NBalanceCheck(filter(nc)%num_soilc,filter(nc)%soilc, &
             filter(nc)%num_soilp, filter(nc)%soilp)
          
        call t_stopf('cnbalchk')
     end if
#endif
        

     ! ============================================================================
     ! Determine albedos for next time step
     ! ============================================================================

     if (doalb) then
        call t_startf('surfalb')
        call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, &
                           caldayp1, declinp1)
        call t_stopf('surfalb')
     end if

  end do
!$OMP END PARALLEL DO
#if !defined (USE_OMP)
!!CSD$ END PARALLEL DO
#endif

end subroutine driver1

!-----------------------------------------------------------------------
!
! !ROUTINE: driver2
!
! !INTERFACE:
subroutine driver2(caldayp1, declinp1, rstwr)
!
! !ARGUMENTS:
  implicit none
  real(r8),          intent(in) :: caldayp1 ! calendar day for nstep+1
  real(r8),          intent(in) :: declinp1 ! declination angle for next time step
  logical, optional, intent(in) :: rstwr    ! true => write restart file this step
!
! !CALLED FROM:
! program program_off (if COUP_OFFLINE cpp variable is defined)
! program program_csm (if COUP_CSM cpp variable is defined)
! module clm_camMod (if SEQ_MCT or SEQ_ESMF cpp variable is defined)
!
! !REVISION HISTORY:
! 2005.05.22  Mariana Vertenstein creation
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: nstep         ! time step number
  real(r8) :: dtime         ! land model time step (sec)
#if (defined DGVM)
  integer  :: nc, c         ! indices
  integer  :: nclumps       ! number of clumps on this processor
  integer  :: yrp1          ! year (0, ...) for nstep+1
  integer  :: monp1         ! month (1, ..., 12) for nstep+1
  integer  :: dayp1         ! day of month (1, ..., 31) for nstep+1
  integer  :: secp1         ! seconds into current date for nstep+1
  integer  :: yr            ! year (0, ...)
  integer  :: mon           ! month (1, ..., 12)
  integer  :: day           ! day of month (1, ..., 31)
  integer  :: sec           ! seconds of the day
  integer  :: ncdate        ! current date
  integer  :: nbdate        ! base date (reference date)
  integer  :: kyr           ! thousand years, equals 2 at end of first year
  integer  :: begp, endp    ! clump beginning and ending pft indices
  integer  :: begc, endc    ! clump beginning and ending column indices
  integer  :: begl, endl    ! clump beginning and ending landunit indices
  integer  :: begg, endg    ! clump beginning and ending gridcell indices
#endif
  character(len=256) :: filer       ! restart file name
  logical :: write_restart
!-----------------------------------------------------------------------

  ! ============================================================================
  ! Write global average diagnostics to standard output
  ! ============================================================================

  call t_startf('wrtdiag')
  nstep = get_nstep()
  call write_diagnostic(wrtdia, nstep)
  call t_stopf('wrtdiag')

#if (defined RTM)
  ! ============================================================================
  ! Route surface and subsurface runoff into rivers
  ! ============================================================================

  call t_startf('clmrtm')
  call Rtmriverflux()
  call t_stopf('clmrtm')
#endif

#if (defined SEQ_MCT) || (defined SEQ_ESMF)
  ! ============================================================================
  ! Read initial snow and soil moisture data at each time step
  ! ============================================================================

  call t_startf('inicperp')
  if (is_perpetual()) call inicfile_perp()
  call t_stopf('inicperp')
#endif

  ! ============================================================================
  ! Update accumulators
  ! ============================================================================

  call t_startf('accum')
  call updateAccFlds()
  call t_stopf('accum')

  ! ============================================================================
  ! Update history buffer
  ! ============================================================================

  call t_startf('hbuf')
  call update_hbuf()
  call t_stopf('hbuf')

  ! ============================================================================
  ! Call DGVM (Dynamic Global Vegetation Model) if appropriate
  ! LPJ is called at last time step of year. Then reset vegetation distribution
  ! and some counters for the next year.
  ! NOTE: monp1, dayp1, and secp1 correspond to nstep+1
  ! NOTE: lpjreset must be called after update_accum and update_hbuf
  ! in order to have the correct values of the accumulated variables
  ! ============================================================================

#if (defined DGVM)
  call t_startf('d2dgvm')
  dtime = get_step_size()
  call get_curr_date(yrp1, monp1, dayp1, secp1, offset=int(dtime))
  if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then

     ! Get date info.  kyr is used in lpj().  At end of first year, kyr = 2.
     call get_curr_date(yr, mon, day, sec)
     ncdate = yr*10000 + mon*100 + day
     call get_ref_date(yr, mon, day, sec)
     nbdate = yr*10000 + mon*100 + day
     kyr = ncdate/10000 - nbdate/10000 + 1

     if (masterproc) write(6,*) 'End of year. DGVM called now: ncdate=', &
                     ncdate,' nbdate=',nbdate,' kyr=',kyr,' nstep=', nstep

     nclumps = get_proc_clumps()

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
#if !defined (USE_OMP)
!!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
#endif
     do nc = 1,nclumps
        call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
        call lpj(begg, endg, begp, endp, filter(nc)%num_natvegp, filter(nc)%natvegp, kyr)
        call lpjreset(begg, endg, begc, endc, begp, endp, filter(nc)%num_nolakep, filter(nc)%nolakep, &
                      caldayp1, declinp1)
        call resetWeightsDGVM(begg, endg, begc, endc, begp, endp)
        call resetTimeConstDGVM(begp, endp)
     end do
#if !defined (USE_OMP)
!!CSD$ END PARALLEL DO
#endif
!$OMP END PARALLEL DO
  end if
  call t_stopf('d2dgvm')
#endif

  ! ============================================================================
  ! Create history and write history tapes if appropriate
  ! ============================================================================
  call t_startf('clm_driver_io')
  call htapes_wrapup()

  ! ============================================================================
  ! Write to DGVM history buffer if appropriate
  ! ============================================================================

#if (defined DGVM)
  if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
     call histDGVM()
     if (masterproc) write(6,*) 'Annual DGVM calculations are complete'
  end if
#endif

  ! ============================================================================
  ! Write restart/initial files if appropriate
  ! ============================================================================

  if (present(rstwr)) then
     write_restart = rstwr
  else
     write_restart = do_restwrite()
  end if
     
  if (write_restart) then
     filer = restFile_filename(type='netcdf')
     call restFile_write( filer )
     filer = restFile_filename(type='binary')
     call restFile_write_binary( filer )
  else if (do_inicwrite()) then
     dtime = get_step_size()
     filer = restFile_filename(type='netcdf', offset=int(dtime))
     call restFile_write( filer )
  end if
  call t_stopf('clm_driver_io')


end subroutine driver2

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diagnostic
!
! !INTERFACE:
subroutine write_diagnostic (wrtdia, nstep)
!
! !DESCRIPTION:
! Write diagnostic surface temperature output each timestep.  Written to
! be fast but not bit-for-bit because order of summations can change each
! timestep.
!
! !USES:
  use clm_atmlnd , only : clm_l2a
  use decompMod  , only : get_proc_bounds, get_proc_global
  use spmdMod    , only : masterproc, npes, MPI_REAL8, MPI_ANY_SOURCE, &
                          MPI_STATUS_SIZE, mpicom
  use shr_sys_mod, only : shr_sys_flush
  use abortutils , only : endrun
!
! !ARGUMENTS:
  implicit none
  logical, intent(in) :: wrtdia     !true => write diagnostic
  integer, intent(in) :: nstep      !model time step
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: p                       ! loop index
  integer :: begp, endp              ! per-proc beginning and ending pft indices
  integer :: begc, endc              ! per-proc beginning and ending column indices
  integer :: begl, endl              ! per-proc beginning and ending landunit indices
  integer :: begg, endg              ! per-proc gridcell ending gridcell indices
  integer :: numg                    ! total number of gridcells across all processors
  integer :: numl                    ! total number of landunits across all processors
  integer :: numc                    ! total number of columns across all processors
  integer :: nump                    ! total number of pfts across all processors
  integer :: ier                     ! error status
  real(r8):: psum                    ! partial sum of ts
  real(r8):: tsum                    ! sum of ts
  real(r8):: tsxyav                  ! average ts for diagnostic output
  integer :: status(MPI_STATUS_SIZE) ! mpi status
!------------------------------------------------------------------------

#if (!defined SEQ_MCT) && (!defined SEQ_ESMF)

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)

  if (wrtdia) then

     call t_barrierf('sync_write_diag', mpicom)
     psum = sum(clm_l2a%t_rad(begg:endg))
     if (masterproc) then
        tsum = psum
        do p = 1, npes-1
           call mpi_recv(psum, 1, MPI_REAL8, p, 999, mpicom, status, ier)
           if (ier/=0) then
              write(6,*) 'write_diagnostic: Error in mpi_recv()'
              call endrun
           end if
           tsum = tsum + psum
        end do
     else
        call mpi_send(psum, 1, MPI_REAL8, 0, 999, mpicom, ier)
        if (ier/=0) then
           write(6,*) 'write_diagnostic: Error in mpi_send()'
           call endrun
        end if
     end if
     if (masterproc) then
        tsxyav = tsum / numg
        write (6,1000) nstep, tsxyav
#ifndef UNICOSMP
        call shr_sys_flush(6)
#endif
     end if

  else

     if (masterproc) then
        write(6,*)'clm2: completed timestep ',nstep
#ifndef UNICOSMP
        call shr_sys_flush(6)
#endif
     end if

  endif

1000 format (1x,'nstep = ',i10,'   TS = ',e21.15)

#endif

end subroutine write_diagnostic

!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_restwrite
!
! !INTERFACE:
logical function do_restwrite()
!
! !DESCRIPTION:
! Determine if restart dataset is to be written at this time step
!
! !USES:
  use restFileMod , only : rest_flag
#if (defined COUP_CSM)
  use clm_csmMod  , only : csmstop_next, csmrstrt
#else
  use clm_time_manager, only : is_last_step
  use histFileMod , only : if_writrest
#endif
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

  do_restwrite = .false.

#if (defined OFFLINE)

  ! Write restart if end of run or if time to dispose master history file

  if (is_last_step() .or. if_writrest) do_restwrite = .true.
  if (.not.rest_flag) do_restwrite = .false.

#elif (defined COUP_CSM)

  ! Write restart only if coupler says to

  if (csmrstrt) do_restwrite = .true.

#endif

end function do_restwrite

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_inicwrite
!
! !INTERFACE:
  logical function do_inicwrite()
!
! !DESCRIPTION:
! Determine if initial dataset is to be written at this time step
! True implies that the initial file will be written one time step
! before the date contained in the filename.
!
! !USES:
    use clm_time_manager, only : get_curr_date, get_prev_date, get_step_size
    use clm_varctl  , only : hist_crtinic
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein

!EOP
!
! !LOCAL VARIABLES:
    integer :: yr         !nstep year (0 -> ...)
    integer :: yrm1       !nstep-1 year (0 -> ...)
    integer :: daym1      !nstep-1 day (1 -> 31)
    integer :: day        !nstep day (1 -> 31)
    integer :: mon        !nstep month (1 -> 12)
    integer :: monm1      !nstep-1 month (1 -> 12)
    integer :: mcsec      !nstep time of day [seconds]
    integer :: mcsecm1    !nstep-1 time of day [seconds]
    integer :: mcsecp1    !nstep+1 time of day [seconds]
    integer :: dayp1      !nstep+1 day (1 -> 31)
    integer :: monp1      !nstep+1 month (1 -> 12)
    integer :: yrp1       !nstep+1 year (0 -> ...)
    integer :: dtime      !timestep size [seconds]
!-----------------------------------------------------------------------

    ! Set calendar for current, previous, and next time steps

    dtime = get_step_size()
    call get_curr_date (yr  , mon  , day  , mcsec  )
    call get_prev_date (yrm1, monm1, daym1, mcsecm1)
    call get_curr_date (yrp1, monp1, dayp1, mcsecp1, offset=dtime)

    ! Determine if time to write out initial dataset

    do_inicwrite = .false.
    if (hist_crtinic /= 'NONE') then
       if      (hist_crtinic == '6-HOURLY') then
          if (mod(mcsecp1,21600) == 0) do_inicwrite = .true.
       elseif  (hist_crtinic == 'DAILY') then
          if (day /= dayp1)  do_inicwrite = .true.
       else if (hist_crtinic == 'MONTHLY') then
          if (mon /= monp1)  do_inicwrite = .true.
       else if (hist_crtinic == 'YEARLY') then
          if (mon == 12 .and. monp1 == 1)  do_inicwrite = .true.
       endif
    endif

  end function do_inicwrite

end module driver
