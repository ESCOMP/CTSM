#include <misc.h>
#include <preproc.h>

#if (defined OFFLINE)

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: program_off
!
! !INTERFACE:
PROGRAM program_off
!
! !DESCRIPTION:
! "off-line" code to mimic coupling to an atmospheric model.
! This program is an "off-line" driver for clm2.
! This code can be used to run the clm2 uncoupled from any atmospheric model.
! The appropriate atmospheric forcing is provided in module [atmdrvMod.F90]
! o If running as an offline driver, the land surface model may use
!   a different grid than the input atmospheric data. The atmospheric
!   data is then interpolated to the land model grid inside the
!   atmospheric driver module [atmdrvMod.F90].
! o If running as part of cam, the land surface model must use the
!   same grid as the cam.
! o If running through the flux coupler, the land surface model grid
!   is interpolated to the atmospheric grid inside the flux coupler
! o To map from the atmospheric grid to the land grid, the atmospheric
!   model must provide latitudes and longitudes (degrees) for each grid
!   point and the North, East, South, and West edges of atmospheric grid.
!   Comparable data for the land grid are provided by the land model.
!   When mapping from land to atm grid, an atm grid cell that is part
!   land and part ocean (as defined by the land surface grid) will have
!   fluxes only based on the land portion.
! o The zenith angle calculation is for the NEXT time step rather
!   than the current time step. Make sure the calendar day is for
!   the NEXT time step. Make sure the calendar day is for Greenwich
!   time (see next comment).
! o The land surface model calculates its own net solar radiation and
!   net longwave radiation at the surface. The net longwave radiation
!   at the surface will differ somewhat from that calculated in the
!   atmospheric model because the atm model will use the upward
!   longwave flux (or radiative temperature) from the previous time
!   step whereas the land surface model uses the flux for the current
!   time step. The net solar radiation should equal that calculated
!   in the atmospheric model. If not, there is a problem in how
!   the models are coupled.
!
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_orb_mod          
  use clm_varpar
  use clm_varctl    , only : irad, nsrest, finidat
  use initializeMod , only : initialize
  use atmdrvMod     , only : atmdrv
  use time_manager  , only : is_last_step, advance_timestep, get_nstep, get_step_size, &
                             get_curr_calday
  use atmdrvMod     , only : atm_getgrid
  use initSurfAlbMod, only : initSurfAlb, do_initsurfalb 
  use driver        , only : driver1, driver2
  use lnd2atmMod    , only : lnd2atm
  use abortutils    , only : endrun
#if (defined SPMD)
  use spmdMod       , only : masterproc, iam, mpicom, spmd_init
#else
  use spmdMod       , only : masterproc, iam
#endif
!
! !ARGUMENTS:
    implicit none
#include <gptl.inc>
#if (defined HAVE_PAPI)
#include <f77papi.h>
#endif
!
! !REVISION HISTORY:
! Author: Gordon Bonan and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  logical  :: doalb     ! true if surface albedo calculation time step
  integer  :: nstep     ! time step index
  real(r8) :: dtime     ! time step increment (sec)
  real(r8) :: caldayp1  ! calendar day for nstep+1
  real(r8) :: calday    ! calendar day for nstep
  real(r8) :: caldaym1  ! calendar day for nstep-1
  real(r8) :: declinp1  ! solar declination angle in radians for nstep+1
  real(r8) :: declin    ! solar declination angle in radians for nstep
  real(r8) :: declinm1  ! solar declination angle in radians for nstep-1
  integer  :: ier       ! error code

! Earth's orbital characteristics

  integer  :: iyear_AD  ! Year (AD) to simulate above earth's orbital parameters for
  real(r8) :: eccen     ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)
  real(r8) :: obliq     ! Earth's obliquity angle (degree's) (-90 to +90) (typically 22-26)
  real(r8) :: mvelp     ! Earth's moving vernal equinox at perhelion (degree's) (0 to 360.0)

! Orbital information after call to routine shr_orbit_params

  real(r8) :: obliqr    ! Earth's obliquity in radians
  real(r8) :: lambm0    ! Mean longitude (radians) of perihelion at the vernal equinox
  real(r8) :: mvelpp    ! Earth's moving vernal equinox longitude
  logical  :: log_print ! true=> print diagnostics
  real(r8) :: eccf      ! earth orbit eccentricity factor
!-----------------------------------------------------------------------

  ! -----------------------------------------------------------------
  ! Initialize timing library
  ! -----------------------------------------------------------------

  ! Set options and initialize timing library.  Use the new "GPTL"
  ! initialization function calls instead of the old "t_" subroutines to
  ! enable CAM to gracefully abort in case of a failure return.  Failures can
  ! happen particularly when enabling PAPI timers and everything was not set
  ! up exactly correctly when building the PAPI lib, or CAM.
  ! 
  ! For logical settings, 2nd arg 0 to gptlsetoption means disable, 
  ! non-zero means enable
  !
  ! Turn off CPU timing (expensive)

  if (gptlsetoption (gptlcpu, 0) < 0) call endrun ('CLM: gptlsetoption')

  ! Compile-time setting of max timer depth

#ifdef GPTLDEPTHLIMIT
  if (gptlsetoption (gptldepthlimit, GPTLDEPTHLIMIT) < 0) call endrun ('CLM: gptlsetoption')
#endif

  ! Next 2 calls only work if PAPI is enabled.  These examples enable counting
  ! of total cycles and floating point ops, respectively
  !
  !   if (gptlsetoption (PAPI_TOT_CYC, 1) < 0) call endrun ('CLM: gptlsetoption')
  !   if (gptlsetoption (PAPI_FP_OPS, 1)  < 0) call endrun ('CLM: gptlsetoption')
  !
  ! Initialize the timing lib.  This call must occur after all gptlsetoption
  ! calls and before all other timing lib calls.

  if (gptlinitialize () < 0) call endrun ('CLM: gptlinitialize')

  ! -----------------------------------------------------------------
  ! Initialize intra-model MPI communication
  ! -----------------------------------------------------------------

#if (defined SPMD)
  call spmd_init()
#endif

  ! -----------------------------------------------------------------
  ! Initialize Orbital parameters
  ! -----------------------------------------------------------------

  ! obliq, eccen and mvelp are determined based on value of iyear_AD

  if (masterproc) then
     log_print = .true.
  else
     log_print = .false.
  end if
  iyear_AD = 1950
  obliq    = SHR_ORB_UNDEF_REAL
  eccen    = SHR_ORB_UNDEF_REAL
  mvelp    = SHR_ORB_UNDEF_REAL
  call shr_orb_params (iyear_AD, eccen, obliq, mvelp, obliqr, &
                       lambm0, mvelpp, log_print)

  ! -----------------------------------------------------------------
  ! Initialize land model
  ! -----------------------------------------------------------------

  call initialize()

  ! Initialize albedos (correct pft filters are needed)
  ! Only done for arbitrary initialization
  
  if (nsrest == 0) then
     if (finidat == ' ' .or. do_initsurfalb) then
        calday = get_curr_calday()
        call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
        
        dtime = get_step_size()
        caldaym1 = get_curr_calday(offset=-int(dtime))
        call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
        
        call initSurfAlb( calday, declin, declinm1 )
     end if
  end if
  
  ! -----------------------------------------------------------------
  ! Initialize "external" atmospheric forcing
  ! -----------------------------------------------------------------

  ! Read atmospheric forcing dataset one time to obtain the longitudes
  ! and latitudes of the atmospheric dataset, as well as the edges. When
  ! coupled to atm model, these are input variables. If no
  ! atmospheric data files are provided, model uses dummy atmospheric
  ! forcing and sets atmospheric grid to land grid.
  
  if (masterproc) write (6,*) 'Attempting to set up atmospheric grid '
  call atm_getgrid()
  if (masterproc) write (6,*) 'Successfully set up atmospheric grid '
  
  ! -----------------------------------------------------------------
  ! Time stepping loop
  ! -----------------------------------------------------------------

  call t_startf('total')

  do
     ! Current atmospheric state and fluxes for all [atmlon] x [atmlat] points.

     nstep = get_nstep()
     call atmdrv(nstep)

     ! doalb is true when the next time step is a radiation time step

     doalb = (irad==1 .or. (mod(nstep,irad)==0 .and. nstep/=0))

     ! Determin declination angle for next time step

     dtime = get_step_size()
     caldayp1 = get_curr_calday( offset=int(dtime) )
     call shr_orb_decl( caldayp1, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )

     ! Call land surface model driver1

     call t_startf('clm_driver1')
     call driver1(doalb, caldayp1, declinp1)
     call t_stopf('clm_driver1')

     ! Determine fields that would be sent to atm for diagnostic purposes
     ! When not in offline mode, this is called from clm_csmMod and lp_coupling

     call lnd2atm()

     ! Call land surface model driver2

     call t_startf('clm_driver2')
     call driver2(caldayp1, declinp1)
     call t_stopf('clm_driver2')

     ! Determine if time to stop

     if (is_last_step()) exit

     ! Increment time step

     call advance_timestep()

  end do
  call t_stopf('total')

  ! -----------------------------------------------------------------
  ! Exit gracefully
  ! -----------------------------------------------------------------

  if (masterproc) then
     write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
  endif
  call t_prf(iam)
#if (defined SPMD)
  call mpi_barrier (mpicom, ier)
  call mpi_finalize(ier)
#endif

  stop
end program program_off

#else

!The following is only here since empty file won't compile
subroutine program_off_stub
  write(6,*) 'PROGRAM_OFF: this routine should not be called'
  return
end subroutine program_off_stub

#endif
