#include <misc.h>
#include <preproc.h>

#if (defined COUP_CSM)

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: program_csm
!
! !INTERFACE:
PROGRAM program_csm
!
! !DESCRIPTION:
! Driver for CLM as the land component of CCSM.
! This program is the driver for CLM to work as the land component of
! CCSM.  The flux coupler will provide all the appropriate atmospheric
! forcing for the land model to run.
! o the land surface model returns to the CCSM flux coupler surface
!   fluxes, temperatures, and albedos for only the land points on the
!   [lsmlon x lsmlat] grid
! o the land surface model uses its own surface type data set. because
!   it processes only land points, this data set must correctly define
!   what points are land and what are not land
! o the land surface model uses its own grid dimensions (lsmlon and
!   lsmlat). currently these must equal lsmlon and lsmlat so that there
!   is a direct correspondence between the atmosphere and land grids
! o the zenith angle calculation is calculated for the
!   NEXT time step rather than the current time step. make sure
!   the calendar day is for the NEXT time step. make sure the
!   solar declination calculation is the same as in the
!   atmospheric model but for the NEXT time step. make sure the
!   calendar day is for greenwich time (see next comment).
! o subroutine calendr: this generates a julian day (with fraction)
!   based on the time step, which is used to calculate the solar
!   zenith angle. this time must be at the greenwich meridian to
!   get the correct zenith angle. also, output from this subroutine
!   is used to calculate the month (1, ..., 12), day (1, ..., 31),
!   and year (00, ...) of the simulation.
! o the land surface model calculates its own net solar radiation and
!   net longwave radiation at the surface. the net longwave radiation
!   at the surface will differ somewhat from that calculated from the
!   CCSM flux coupler because the cpl model will use the upward
!   longwave flux (or radiative temperature) from the previous time
!   step whereas the land surface model uses the flux for the current
!   time step. the net solar radiation should equal that calculated
!   from the flux coupler. if not, there is a problem.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_orb_mod        
  use shr_msg_mod        
  use clm_varpar    , only : lsmlon, lsmlat     
  use clm_varctl    , only : nsrest, irad, csm_doflxave, finidat
  use clm_varorb    , only : eccen, mvelpp, lambm0, obliqr
  use time_manager  , only : advance_timestep, get_nstep, get_curr_calday, get_step_size
  use clm_csmMod    , only : csmstop_now, csm_setup, csm_shutdown, & 
                             csm_dosndrcv, csm_recv, csm_send, csm_flxave, &
                             csm_initialize, csm_sendalb, dorecv, dosend  
  use clm_comp      , only : clm_init1, clm_init2, clm_run1, clm_run2
#if (defined SPMD)
  use spmdMod       , only : masterproc, iam, spmd_init, mpicom
#else
  use spmdMod       , only : masterproc, iam
#endif
  use abortutils    , only : endrun
!
! !ARGUMENTS:
    implicit none
#include <gptl.inc>
#if (defined HAVE_PAPI)
#include <f77papi.h>
#endif
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: i,j          ! loop indices
  integer  :: nstep        ! time step index
  real(r8) :: dtime        ! time step increment (sec)
  logical  :: doalb        ! true if surface albedo calculation time step
  real(r8) :: caldayp1     ! calendar day for nstep+1
  real(r8) :: calday       ! calendar day for nstep
  real(r8) :: caldaym1     ! calendar day for nstep-1
  real(r8) :: declinp1     ! solar declination angle in radians for nstep+1
  real(r8) :: declin       ! solar declination angle in radians for nstep
  real(r8) :: declinm1     ! solar declination angle in radians for nstep-1
  integer  :: ier          ! error code
  logical  :: log_print    ! true=> print diagnostics
  integer  :: mpicom_dummy ! temporary
!
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
  ! Initailize input/output units
  ! -----------------------------------------------------------------

  call shr_msg_stdio ('lnd')

  ! -----------------------------------------------------------------
  ! Initialize inter-model MPI communication 
  ! -----------------------------------------------------------------

#if (defined SPMD)
  call csm_setup(mpicom)
#else
  call csm_setup(mpicom_dummy)
#endif

  ! -----------------------------------------------------------------
  ! Initialize intra-model MPI communication 
  ! -----------------------------------------------------------------

#if (defined SPMD)
  call spmd_init()
#endif

  ! -----------------------------------------------------------------
  ! Initialize land model
  ! -----------------------------------------------------------------

  ! Initialize land model - initialize communication with flux coupler

  call clm_init1()

  ! Initialize flux coupler communication - need to call get orbital 
  ! parameters before the call to initialize albedos
  
  call csm_initialize(eccen, obliqr, lambm0, mvelpp)

  ! Initialize albedos (correct pft filters are needed)

  call clm_init2() 

  ! Send first land model data to flux coupler.

  call csm_sendalb()

  ! -----------------------------------------------------------------
  ! Time stepping loop
  ! -----------------------------------------------------------------

  call t_startf('lnd_timeloop')
  do

     ! doalb is true when the next time step is a radiation time step

     nstep = get_nstep()
     doalb = ((irad==1 .and. nstep/=0) .or. (mod(nstep,irad)==0 .and. nstep/=0))

     ! Determine if information should be sent/received to/from flux coupler

     call csm_dosndrcv(doalb)

     ! Get atmospheric state and fluxes from flux coupler and 

     if (dorecv) then
        ! Fill in a2ls and a2lf
        call csm_recv()

        ! Determine if time to stop
        if (csmstop_now) then
           call t_stopf('clm_driver1')
           exit
        end if
     end if

     call clm_run1()

     ! Average fluxes over interval if appropriate
     ! Surface states sent to the flux coupler states are not time averaged

     if (csm_doflxave) call csm_flxave()

     ! Send fields to flux coupler
     ! Send states[n] (except for snow[n-1]), time averaged fluxes for [n,n-1,n-2],
     ! albedos[n+1], and ocnrof_vec[n]

     if (dosend) call csm_send()

     call clm_run2()
     
     ! Increment time step

     call advance_timestep()

  end do
  call t_stopf ('lnd_timeloop')

  ! -----------------------------------------------------------------
  ! Exit gracefully
  ! -----------------------------------------------------------------

  if (masterproc) then
     write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
  endif
  call t_prf(iam)
  call csm_shutdown()

  stop
end program program_csm

#else

!The following is only here since empty file won't compile
subroutine program_csm_stub
  write(6,*) 'PROGRAM_CSM: this routine should not be called'
  stop 99
end subroutine program_csm_stub

#endif
