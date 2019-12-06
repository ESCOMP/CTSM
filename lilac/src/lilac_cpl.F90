module lilac_cpl

  !-----------------------------------------------------------------------
  ! Module containing all routines for couplers
  !-----------------------------------------------------------------------

  use ESMF
  use shr_sys_mod  , only : shr_sys_abort
  use lilac_methods, only : chkerr

  implicit none
  private

  public :: cpl_atm2lnd_register
  public :: cpl_lnd2atm_register
  public :: cpl_lnd2rof_register
  public :: cpl_rof2lnd_register

  type(ESMF_RouteHandle) :: rh_atm2lnd
  type(ESMF_RouteHandle) :: rh_lnd2atm
  type(ESMF_RouteHandle) :: rh_lnd2rof
  type(ESMF_RouteHandle) :: rh_rof2lnd

  integer                :: mytask
  integer,  parameter    :: ispval_mask = -987987 ! spval for RH mask values

  character(*), parameter :: modname =  "lilac_cpl"

  character(*), parameter :: u_FILE_u = &
       __FILE__

!======================================================================
contains
!======================================================================

  subroutine cpl_atm2lnd_register(cplcomp, rc)

    ! input/output variables
    type(ESMF_CplComp   ) :: cplcomp
    integer, intent(out ) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    character(len=*) , parameter :: subname=trim(modname ) //' : [cpl_atm2lnd_register] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
       print *, "in cpl_atm2lnd_register routine"
    end if

    ! Register the callback routines.
    ! Set the entry points for coupler ESMF Component methods
    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_atm2lnd_init, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN       , userRoutine=cpl_atm2lnd_run  , rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE  , userRoutine=cpl_atm2lnd_final, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

  end subroutine cpl_atm2lnd_register

!======================================================================

  subroutine cpl_lnd2atm_register(cplcomp, rc)

    type(ESMF_CplComp)    :: cplcomp
    integer, intent(out ) :: rc

    ! local variables
    character(len=*     ) , parameter :: subname=trim(modname ) //' : [cpl_lnd2atm_register] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "in cpl_lnd2atm_register routine"
    end if

    ! Register the callback routines.
    ! Set the entry points for coupler ESMF Component methods
    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_lnd2atm_init, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN, userRoutine=cpl_lnd2atm_run  , rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE, userRoutine=cpl_lnd2atm_final, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

  end subroutine cpl_lnd2atm_register

!======================================================================

  subroutine cpl_lnd2rof_register(cplcomp, rc)

    ! input/output variables
    type(ESMF_CplComp   ) :: cplcomp
    integer, intent(out ) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    character(len=*) , parameter :: subname=trim(modname ) //' : [cpl_atm2lnd_register] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
       print *, "in cpl_atm2lnd_register routine"
    end if

    ! Register the callback routines.
    ! Set the entry points for coupler ESMF Component methods
    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_lnd2rof_init, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN       , userRoutine=cpl_lnd2rof_run  , rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE  , userRoutine=cpl_lnd2rof_final, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

  end subroutine cpl_lnd2rof_register

!======================================================================

  subroutine cpl_rof2lnd_register(cplcomp, rc)

    type(ESMF_CplComp)    :: cplcomp
    integer, intent(out ) :: rc

    ! local variables
    character(len=*     ) , parameter :: subname=trim(modname ) //' : [cpl_rof2lnd_register] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "in cpl_rof2lnd_register routine"
    end if

    ! Register the callback routines.
    ! Set the entry points for coupler ESMF Component methods
    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine=cpl_rof2lnd_init, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN       , userRoutine=cpl_rof2lnd_run  , rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE  , userRoutine=cpl_rof2lnd_final, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

  end subroutine cpl_rof2lnd_register

!======================================================================

  subroutine cpl_atm2lnd_init(cplcomp, importState, exportState, clock, rc)

    ! input/output variables
    type (ESMF_CplComp    ) :: cplcomp
    type (ESMF_State      ) :: importState
    type (ESMF_State      ) :: exportState
    type (ESMF_Clock      ) :: clock
    integer, intent(out   ) :: rc

    ! local variables
    type (ESMF_FieldBundle)         :: import_fieldbundle
    type (ESMF_FieldBundle)         :: export_fieldbundle
    character(len=*), parameter     :: subname=trim(modname) //': [cpl_atm2lnd_init] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Coupler for atmosphere to land initialize routine called"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call cpl_get_fieldbundle(importState, 'a2c_fb', import_fieldbundle, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call cpl_get_fieldbundle(exportState, 'c2l_fb_atm', export_fieldbundle, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in initializing cpl_atm2lnd')

    call ESMF_LogWrite(subname//"cpl_atm2lnd_init finished!", ESMF_LOGMSG_INFO)

  end subroutine cpl_atm2lnd_init

!======================================================================

  subroutine cpl_lnd2atm_init(cplcomp, importState, exportState, clock, rc)

    type (ESMF_CplComp     ) :: cplcomp
    type (ESMF_State       ) :: importState
    type (ESMF_State       ) :: exportState
    type (ESMF_Clock       ) :: clock
    integer, intent(out    ) :: rc

    ! local variables
    type (ESMF_FieldBundle)         :: import_fieldbundle
    type (ESMF_FieldBundle)         :: export_fieldbundle
    character(len=*) , parameter    :: subname=trim(modname ) //': [cpl_lnd2atm_init] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Coupler for land to atmosphere initialize routine called"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call cpl_get_fieldbundle(importState, 'l2c_fb_atm', import_fieldbundle, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call cpl_get_fieldbundle(exportState, 'c2a_fb', export_fieldbundle, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in initializing cpl_lnd2atm')

    call ESMF_LogWrite(subname//"cpl init finished!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2atm_init

!======================================================================

  subroutine cpl_lnd2rof_init(cplcomp, importState, exportState, clock, rc)

    ! input/output variables
    type (ESMF_CplComp    ) :: cplcomp
    type (ESMF_State      ) :: importState
    type (ESMF_State      ) :: exportState
    type (ESMF_Clock      ) :: clock
    integer, intent(out   ) :: rc

    ! local variables
    type (ESMF_FieldBundle)     :: import_fieldbundle
    type (ESMF_FieldBundle)     :: export_fieldbundle
    integer                     :: srcTermProcessing_Value = 0 ! should this be a module variable?
    character(len=*), parameter :: subname=trim(modname) //': [cpl_lnd2rof_init] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Coupler for atmosphere to land initialize routine called"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call cpl_get_fieldbundle(importState, 'l2c_fb_rof', import_fieldbundle, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call cpl_get_fieldbundle(exportState, 'c2r_fb', export_fieldbundle, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleRegridStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2rof, &
         srcMaskValues=(/ispval_mask/), dstMaskValues=(/ispval_mask/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
         normType=ESMF_NORMTYPE_FRACAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., &
         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in initializing cpl_lnd2rof')

    call ESMF_LogWrite(subname//"cpl init finished!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2rof_init

!======================================================================

  subroutine cpl_rof2lnd_init(cplcomp, importState, exportState, clock, rc)

    type (ESMF_CplComp     ) :: cplcomp
    type (ESMF_State       ) :: importState
    type (ESMF_State       ) :: exportState
    type (ESMF_Clock       ) :: clock
    integer, intent(out    ) :: rc

    ! local variables
    type (ESMF_FieldBundle)      :: import_fieldbundle
    type (ESMF_FieldBundle)      :: export_fieldbundle
    integer                      :: srcTermProcessing_Value = 0 ! should this be a module variable?
    character(len=*) , parameter :: subname=trim(modname ) //': [cpl_rof2lnd_init] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    if (mytask == 0) then
       print *, "Coupler for land to atmosphere initialize routine called"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call cpl_get_fieldbundle(importState, 'r2c_fb', import_fieldbundle, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call cpl_get_fieldbundle(exportState, 'c2l_fb_rof', export_fieldbundle, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleRegridStore(import_fieldbundle, export_fieldbundle, routehandle=rh_rof2lnd, &
         srcMaskValues=(/ispval_mask/), dstMaskValues=(/ispval_mask/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
         normType=ESMF_NORMTYPE_FRACAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., &
         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in initializing cpl_rof2lnd')

    call ESMF_LogWrite(subname//"cpl init finished!", ESMF_LOGMSG_INFO)

  end subroutine cpl_rof2lnd_init

!======================================================================

  subroutine cpl_atm2lnd_run(cplcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_CplComp)   :: cplcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle ) :: import_fieldbundle, export_fieldbundle
    character(len=*) , parameter :: subname=trim(modname       ) //': [cpl_atm2lnd_run] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Running cpl_atm2lnd_run"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(importState, "a2c_fb", import_fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(exportState, "c2l_fb_atm", export_fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleRedist(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//" regridding fieldbundles from atmos to land!", ESMF_LOGMSG_INFO)

  end subroutine cpl_atm2lnd_run

!======================================================================

  subroutine cpl_lnd2atm_run(cplcomp, importState, exportState, clock, rc)

    type(ESMF_CplComp)   :: cplcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle) :: import_fieldbundle, export_fieldbundle
    character(len=*) , parameter :: subname=trim(modname       ) //': [cpl_lnd2atm_run] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Running cpl_lnd2atm_run"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(importState, "l2c_fb_atm", import_fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(exportState, "c2a_fb", export_fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleRedist(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//" regridding fieldbundles  from land to atmos!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2atm_run

!======================================================================

  subroutine cpl_lnd2rof_run(cplcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_CplComp)   :: cplcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle ) :: import_fieldbundle, export_fieldbundle
    character(len=*        ) , parameter :: subname=trim(modname) //': [cpl_lnd2rof_run] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Running cpl_lnd2rof_run"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(importState, "l2c_fb_rof", import_fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(exportState, "c2r_fb", export_fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleRegrid(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2rof, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//" regridding fieldbundles from land to river!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2rof_run

!======================================================================

  subroutine cpl_rof2lnd_run(cplcomp, importState, exportState, clock, rc)

    type(ESMF_CplComp)   :: cplcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle ) :: import_fieldbundle, export_fieldbundle
    character(len=*) , parameter :: subname=trim(modname) //': [cpl_rof2lnd_run] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Running cpl_rof2lnd_run"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(importState, "r2c_fb", import_fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(exportState, "c2l_fb_rof", export_fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleRegrid(import_fieldbundle, export_fieldbundle, routehandle=rh_rof2lnd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//" regridding fieldbundles from river to land!", ESMF_LOGMSG_INFO)

  end subroutine cpl_rof2lnd_run

!======================================================================

  subroutine cpl_atm2lnd_final(cplcomp, importState, exportState, clock, rc)

    ! input/output variables
    type (ESMF_CplComp)  :: cplcomp
    type (ESMF_State)    :: importState
    type (ESMF_State)    :: exportState
    type (ESMF_Clock)    :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle ) :: import_fieldbundle, export_fieldbundle
    character(len=*) , parameter :: subname=trim(modname       ) //': [cpl_atm2lnd_final] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"---------------------------------!", ESMF_LOGMSG_INFO)

    ! Only thing to do here is release redist (or regrid) and route handles
    call ESMF_FieldBundleRegridRelease (routehandle=rh_atm2lnd, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//" rh_atm2lnd route handle released!", ESMF_LOGMSG_INFO)

  end subroutine cpl_atm2lnd_final

!======================================================================

  subroutine cpl_lnd2atm_final(cplcomp, importState, exportState, clock, rc)

    type (ESMF_CplComp)  :: cplcomp
    type (ESMF_State)    :: importState
    type (ESMF_State)    :: exportState
    type (ESMF_Clock)    :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*) , parameter :: subname=trim(modname) //': [cpl_lnd2atm_final] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"---------------------------------!", ESMF_LOGMSG_INFO)

    ! Only thing to do here is release redist (or regrid) and route handles
    call ESMF_FieldBundleRegridRelease (routehandle=rh_lnd2atm , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//" rh_lnd2atm route handle released!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2atm_final

!======================================================================

  subroutine cpl_lnd2rof_final(cplcomp, importState, exportState, clock, rc)

    type (ESMF_CplComp)  :: cplcomp
    type (ESMF_State)    :: importState
    type (ESMF_State)    :: exportState
    type (ESMF_Clock)    :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*) , parameter :: subname=trim(modname) //': [cpl_lnd2rof_final] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"---------------------------------!", ESMF_LOGMSG_INFO)

    ! Only thing to do here is release redist (or regrid) and route handles
    call ESMF_FieldBundleRegridRelease (routehandle=rh_lnd2rof , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//" rh_lnd2rof route handle released!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2rof_final

!======================================================================

  subroutine cpl_rof2lnd_final(cplcomp, importState, exportState, clock, rc)

    type (ESMF_CplComp)  :: cplcomp
    type (ESMF_State)    :: importState
    type (ESMF_State)    :: exportState
    type (ESMF_Clock)    :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*) , parameter :: subname=trim(modname) //': [cpl_rof2lnd_final] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"---------------------------------!", ESMF_LOGMSG_INFO)

    ! Only thing to do here is release redist (or regrid) and route handles
    call ESMF_FieldBundleRegridRelease (routehandle=rh_rof2lnd , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//" rh_rof2lnd route handle released!", ESMF_LOGMSG_INFO)

  end subroutine cpl_rof2lnd_final

!======================================================================

  subroutine cpl_get_fieldbundle(state, fbname, fieldbundle, rc)

    ! input/output variables
    type(ESMF_State)       :: state
    character(len=*)       :: fbname
    type(ESMF_FieldBundle) :: fieldbundle
    integer, intent(out)   :: rc

    ! local variables
    integer                         :: n
    integer                         :: fieldcount
    character(len=128), allocatable :: fieldlist(:)
    character(len=128)              :: cvalue
    character(len=*), parameter     :: subname=trim(modname) //': [cpl_get_fieldbundle] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, trim(fbname), fieldbundle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(fieldbundle, fieldCount=fieldCount, rc=rc) 
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) fieldcount
    call ESMF_LogWrite(subname//" trim(fbname)//' field count = "//trim(cvalue), ESMF_LOGMSG_INFO)
    allocate(fieldlist(fieldcount))

    call ESMF_FieldBundleGet(fieldbundle, fieldNameList=fieldlist, rc=rc) 
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1,fieldCount
       write(cvalue,*) n 
       call ESMF_LogWrite(subname//trim(fbname)//" field "//trim(cvalue)//' = '//trim(fieldlist(n)), &
            ESMF_LOGMSG_INFO)
    end do
    deallocate(fieldlist)
    if (mytask == 0) then
       print *, trim(fbname)//' field count = ',fieldcount
    end if

  end subroutine cpl_get_fieldbundle

end module lilac_cpl
