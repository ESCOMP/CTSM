module cpl_mod


    use ESMF
    implicit none

    contains



    subroutine cpl_atm2lnd_register(cplcomp, rc)
     type(ESMF_CplComp)   :: cplcomp
     integer, intent(inout) :: rc

      rc = ESMF_FAILURE

      ! Register the callback routines.

      call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_atm2lnd_init, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN, userRoutine=coupler_run, rc=rc)
      !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE, userRoutine=coupler_final, rc=rc)
      !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)


    end subroutine cpl_atm2lnd_register


    subroutine cpl_atm2lnd_init(cplcomp, importState, exportState, clock, rc)

      type(ESMF_CplComp)   :: cplcomp
      type(ESMF_State)     :: importState, exportState
      type(ESMF_Clock)     :: clock
      integer, intent(inout) :: rc

    end subroutine cpl_atm2lnd_init



end module cpl_mod


