module cpl_mod

    use ESMF
    implicit none

    private

    type(ESMF_RouteHandle), save :: rh_atm2lnd, rh_lnd2atm


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
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine= cpl_atm2lnd_init, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN       , userRoutine=cpl_atm2lnd_run  , rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE  , userRoutine=cpl_atm2lnd_final, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    end subroutine cpl_atm2lnd_register

    subroutine cpl_lnd2atm_register(cplcomp, rc)
        type(ESMF_CplComp)   :: cplcomp
        integer, intent(out) :: rc
        character(len=*), parameter :: subname=trim(modname)//' : [cpl_lnd2atm_register] '


        rc = ESMF_SUCCESS
        print *, "in cpl_lnd2atm_register routine"

        ! Register the callback routines.
        ! Set the entry points for coupler ESMF Component methods
        !call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_lnd2atm_init, rc=rc)
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, cpl_lnd2atm_init, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN       , userRoutine=cpl_lnd2atm_run  , rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
        call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE  , userRoutine=cpl_lnd2atm_final, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    end subroutine cpl_lnd2atm_register

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine cpl_lnd2atm_init(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp) :: cplcomp
        type(ESMF_State) :: importState
        type(ESMF_State) :: exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc
        type (ESMF_FieldBundle)      ::  import_fieldbundle, export_fieldbundle

        character(len=*), parameter :: subname=trim(modname)//': [cpl_lnd2atm_init] '

        rc = ESMF_SUCCESS
        print *, "Coupler for land to atmosphere initialize routine called"
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(importState, "l2c_fb", import_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(exportState, "c2a_fb", export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        ! For Redisting
        !call ESMF_FieldBundleRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        ! For ReGridding
        call ESMF_FieldBundleRegridStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        !call ESMF_StateGet(importState, itemname="a2c_fb", item=import_fieldbundle, rc=rc)
        !if (chkerr(rc,__LINE__,u_FILE_u)) return
        !call ESMF_StateGet(exportState, itemname="c2a_fb", item=export_fieldbundle, rc=rc)
        !if (chkerr(rc,__LINE__,u_FILE_u)) return
        !call ESMF_FieldRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        !if (chkerr(rc,__LINE__,u_FILE_u)) return
    end subroutine cpl_lnd2atm_init

    subroutine cpl_atm2lnd_init(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp) :: cplcomp
        type(ESMF_State) :: importState
        type(ESMF_State) :: exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc
        type (ESMF_FieldBundle)      ::  import_fieldbundle, export_fieldbundle

        character(len=*), parameter :: subname=trim(modname)//': [cpl_atm2lnd_init] '

        rc = ESMF_SUCCESS
        print *, "Coupler for atmosphere to land initialize routine called"
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(importState, "a2c_fb", import_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(exportState, "c2l_fb", export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        ! For Redisting
        !call ESMF_FieldBundleRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        ! For ReGridding
        call ESMF_FieldBundleRegridStore(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)
    end subroutine cpl_atm2lnd_init

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine cpl_lnd2atm_run(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp) :: cplcomp
        type(ESMF_State) :: importState
        type(ESMF_State) :: exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc
        type (ESMF_FieldBundle)      ::  import_fieldbundle, export_fieldbundle

        character(len=*), parameter :: subname=trim(modname)//': [cpl_lnd2atm_run] '

        rc = ESMF_SUCCESS

        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)
        print *, "Running cpl_lnd2atm_run"

        call ESMF_StateGet(importState, "l2c_fb", import_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_StateGet(exportState, "c2a_fb", export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        !call ESMF_StateGet(importState, itemname=importStateName, item=srcFieldBundle, rc=rc)
        !if (chkerr(rc,__LINE__,u_FILE_u)) return
        !call ESMF_StateGet(exportState, itemname=exportStateName, item=dstFieldBundle, rc=rc)
        !if (chkerr(rc,__LINE__,u_FILE_u)) return
        !call ESMF_FieldBundleRegrid(srcFieldBundle, dstFieldBundle, rh_lnd2atm, rc=rc)
        call ESMF_FieldBundleRegrid(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//" regridding fieldbundles  from land to atmos!", ESMF_LOGMSG_INFO)
        !routehandle, zeroregion, termorderflag, checkflag, rc)
    end subroutine cpl_lnd2atm_run

    subroutine cpl_atm2lnd_run(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp) :: cplcomp
        type(ESMF_State) :: importState
        type(ESMF_State) :: exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc
        type (ESMF_FieldBundle)      ::  import_fieldbundle, export_fieldbundle

        character(len=*), parameter :: subname=trim(modname)//': [cpl_atm2lnd_run] '

        rc = ESMF_SUCCESS

        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)
        print *, "Running cpl_atm2lnd_run"


        call ESMF_StateGet(importState, trim("a2c_fb"), import_fieldbundle, rc=rc)
        !call ESMF_StateGet(importState, itemName=trim("a2c_fb"), item=import_fieldbundle, rc=rc)     ! this syntax was not working???
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//" got a2c fieldbundle!", ESMF_LOGMSG_INFO)

        call ESMF_StateGet(exportState, trim("c2l_fb"), export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//" got c2l fieldbundle!", ESMF_LOGMSG_INFO)
        call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

        !call ESMF_StateGet(importState, itemname=importStateName, item=srcFieldBundle, rc=rc)
        !if (chkerr(rc,__LINE__,u_FILE_u)) return
        !call ESMF_StateGet(exportState, itemname=exportStateName, item=dstFieldBundle, rc=rc)
        !if (chkerr(rc,__LINE__,u_FILE_u)) return


        call ESMF_FieldBundleRegrid(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//" regridding fieldbundles from atmos to land!", ESMF_LOGMSG_INFO)


        !routehandle, zeroregion, termorderflag, checkflag, rc)
     end subroutine cpl_atm2lnd_run

    subroutine cpl_lnd2atm_final(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp) :: cplcomp
        type(ESMF_State) :: importState
        type(ESMF_State) :: exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc
        type (ESMF_FieldBundle)      ::  import_fieldbundle, export_fieldbundle

        character(len=*), parameter :: subname=trim(modname)//': [cpl_lnd2atm_final] '

        rc = ESMF_SUCCESS

        ! Only thing to do here is release redist (or regrid) and route handles
        call ESMF_FieldBundleRegridRelease (routehandle=rh_lnd2atm , rc=rc)

        call ESMF_LogWrite(subname//"---------------------------------!", ESMF_LOGMSG_INFO)
        call ESMF_LogWrite(subname//" rh_lnd2atm route handle released!", ESMF_LOGMSG_INFO)
    end subroutine cpl_lnd2atm_final


    subroutine cpl_atm2lnd_final(cplcomp, importState, exportState, clock, rc)
        type(ESMF_CplComp) :: cplcomp
        type(ESMF_State) :: importState
        type(ESMF_State) :: exportState
        type(ESMF_Clock) :: clock
        integer, intent(out) :: rc
        type (ESMF_FieldBundle)      ::  import_fieldbundle, export_fieldbundle

        character(len=*), parameter :: subname=trim(modname)//': [cpl_atm2lnd_final] '

        rc = ESMF_SUCCESS

        ! Only thing to do here is release redist (or regrid) and route handles
        call ESMF_FieldBundleRegridRelease (routehandle=rh_atm2lnd, rc=rc)

        call ESMF_LogWrite(subname//"---------------------------------!", ESMF_LOGMSG_INFO)
        call ESMF_LogWrite(subname//" rh_atm2lnd route handle released!", ESMF_LOGMSG_INFO)
    end subroutine cpl_atm2lnd_final






end module cpl_mod

