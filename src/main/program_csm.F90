#include <misc.h>
#include <preproc.h>

#if (defined COUP_CSM)

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: CLM
!
! !INTERFACE:
#ifdef SINGLE_EXEC
 subroutine ccsm_lnd()
#else
 program ccsm_lnd
#endif
!
! !DESCRIPTION:
! Driver for CLM as the land component of CCSM.
! This program is the driver for CLM to work as the land component of
! CCSM.  The flux coupler will provide all the appropriate atmospheric
! forcing for the land model to run.
! o the land surface model returns to the CCSM flux coupler surface
!   fluxes, temperatures, and albedos for only the land points on the

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
#ifdef SINGLE_EXEC
  use MPH_module      , only : MPH_get_argument
#endif
  use shr_kind_mod    , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_orb_mod        
  use shr_file_mod        
  use controlMod      , only : control_setNL
  use clm_varctl      , only : nsrest, irad, csm_doflxave, finidat
  use clm_varorb      , only : eccen, mvelpp, lambm0, obliqr
  use clm_time_manager, only : advance_timestep, get_nstep, get_curr_calday, get_step_size
  use clm_csmMod      , only : csmstop_now, csm_setup, csm_shutdown, & 
                               csm_dosndrcv, csm_recv, csm_send, csm_flxave, &
                               csm_initialize, csm_sendalb, dorecv, dosend  
  use clm_comp        , only : clm_init0, clm_init1, clm_init2, clm_run1, clm_run2
  use abortutils      , only : endrun
  use spmdMod       
  use ESMF_Mod
  use perf_mod
!
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
#ifdef SINGLE_EXEC
  integer  :: nThreads
#endif
  integer  :: i,j          ! loop indices
  integer  :: nstep        ! time step index
  real(r8) :: dtime        ! time step increment (sec)
  logical  :: doalb        ! true if surface albedo calculation time step
  integer  :: ier          ! error code
  logical  :: log_print    ! true=> print diagnostics
  integer  :: mpicom_top   ! MPI communicator
  character(len=SHR_KIND_CL) :: nlfilename ! Namelist filename
!
!-----------------------------------------------------------------------

#ifdef SINGLE_EXEC
  call MPH_get_argument("THREADS", nThreads, "lnd")
#ifdef _OPENMP
   call OMP_SET_NUM_THREADS(nThreads)
#endif
#endif

  ! -----------------------------------------------------------------
  ! Initialize MPI
  ! -----------------------------------------------------------------

  call csm_setup(mpicom_top)
  call spmd_init(mpicom_top)

  ! -----------------------------------------------------------------
  ! Initialize ESMF
  ! -----------------------------------------------------------------

  call ESMF_Initialize()

  ! -----------------------------------------------------------------
  ! Initialize input/output units
  ! -----------------------------------------------------------------

   call shr_file_chdir('lnd')                    ! all PE's chdir
   if (masterproc) then
      call shr_file_chStdin('lnd', nlfilename=nlfilename)  ! redir unit 5
      call control_setNL( nlfilename )                     ! Set namelist
      call shr_file_chStdout('lnd')                        ! redir unit 6
   end if

  ! -----------------------------------------------------------------
  ! Initialize timing library
  ! -----------------------------------------------------------------

  call t_initf(nlfilename, LogPrint=masterproc, Mpicom=mpicom_top, &
               MasterTask=masterproc)

  ! -----------------------------------------------------------------
  ! Initialize land model
  ! -----------------------------------------------------------------

  ! Initialize land model - initialize communication with flux coupler

  call clm_init0()
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
  call t_prf('timing_all',mpicom_top)
  call t_finalizef()
  call csm_shutdown()

  stop

#ifdef SINGLE_EXEC
 end subroutine ccsm_lnd
#else
 end program ccsm_lnd
#endif

#else

!The following is only here since empty file won't compile
subroutine program_csm_stub
  write(6,*) 'PROGRAM_CSM: this routine should not be called'
  stop 99
end subroutine program_csm_stub

#endif
