module cpl_mod

    use ESMF
    implicit none

    public cpl_atm2lnd_register
    public cpl_lnd2atm_register

contains

    subroutine cpl_atm2lnd_register(cplcomp, rc)
     type(ESMF_CplComp)   :: cplcomp
     integer, intent(out) :: rc

      rc = ESMF_FAILURE

      ! Register the callback routines.

      !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_atm2lnd_init, rc=rc)
      !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, cpl_atm2lnd_init ,rc=rc)
      call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, my_init, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN, userRoutine=coupler_run, rc=rc)
      !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
      !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE, userRoutine=coupler_final, rc=rc)
      !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    end subroutine cpl_atm2lnd_register

    subroutine cpl_lnd2atm_register(cplcomp, rc)
     type(ESMF_CplComp)   :: cplcomp
     integer, intent(out) :: rc
     rc = ESMF_FAILURE
      call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, cpl_lnd2atm_init, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)


    end subroutine cpl_lnd2atm_register


    subroutine my_init(cplcomp, importState, exportState, clock, rc)
      type(ESMF_CplComp) :: cplcomp
      type(ESMF_State) :: importState
      type(ESMF_State) :: exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      print *, "CPLR initialize routine called"
      rc = ESMF_SUCCESS
    end subroutine my_init


    subroutine cpl_lnd2atm_init(cplcomp, importState, exportState, clock, rc)
      type(ESMF_CplComp) :: cplcomp
      type(ESMF_State) :: importState
      type(ESMF_State) :: exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      print *, "Coupler for land to atmosphere initialize routine called"
      rc = ESMF_SUCCESS
    end subroutine cpl_lnd2atm_init



    subroutine cpl_atm2lnd_init(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp)   :: cplcomp
        type(ESMF_State)     :: importState, exportState
        type(ESMF_Clock)     :: clock
        integer, intent(inout) :: rc

        print *, "Coupler Init starting"
        rc = ESMF_SUCCESS

    end subroutine cpl_atm2lnd_init



end module cpl_mod


