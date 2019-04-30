module cpl_mod

    use ESMF
    implicit none

    public cpl_atm2lnd_register
    public cpl_lnd2atm_register


    character(*), parameter :: modname =  "  cpl_mod"

    contains

    subroutine cpl_atm2lnd_register(cplcomp, rc)
        type(ESMF_CplComp)   :: cplcomp
        integer, intent(out) :: rc
        character(len=*), parameter :: subname=trim(modname)//':[cpl_atm2lnd_register] '

        rc = ESMF_SUCCESS
        print *, "in cpl_atm2lnd_register routine"

        ! Register the callback routines.
        ! Set the entry points for coupler ESMF Component methods
        !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_atm2lnd_init, rc=rc)
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, cpl_atm2lnd_init, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN, userRoutine=coupler_run, rc=rc)
        !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE, userRoutine=coupler_final, rc=rc)
        !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    end subroutine cpl_atm2lnd_register

    subroutine cpl_lnd2atm_register(cplcomp, rc)
        type(ESMF_CplComp)   :: cplcomp
        integer, intent(out) :: rc
        character(len=*), parameter :: subname=trim(modname)//':[cpl_lnd2atm_register] '


        rc = ESMF_SUCCESS
        print *, "in cpl_lnd2atm_register routine"

        ! Register the callback routines.
        ! Set the entry points for coupler ESMF Component methods
        !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_lnd2atm_init, rc=rc)
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, cpl_lnd2atm_init, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    end subroutine cpl_lnd2atm_register

    subroutine cpl_lnd2atm_init(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp) :: cplcomp
        type(ESMF_State) :: importState
        type(ESMF_State) :: exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc

        character(len=*), parameter :: subname=trim(modname)//':[cpl_lnd2atm_init] '

        print *, "CPLR initialize routine called"
        print *, "Coupler for land to atmosphere initialize routine called"
        rc = ESMF_SUCCESS
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)
    end subroutine cpl_lnd2atm_init

    subroutine cpl_atm2lnd_init(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp) :: cplcomp
        type(ESMF_State) :: importState
        type(ESMF_State) :: exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc

        character(len=*), parameter :: subname=trim(modname)//':[cpl_lnd2atm_init] '

        print *, "CPLR initialize routine called"
        print *, "Coupler for atmosphere to land initialize routine called"
        rc = ESMF_SUCCESS
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)
    end subroutine cpl_atm2lnd_init



end module cpl_mod


