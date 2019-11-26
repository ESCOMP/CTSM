module lilac_cpl

  !-----------------------------------------------------------------------
  ! Module containing all routines for both couplers
  ! 1- coupler 1 : atm ---> lnd (cpl_atm2lnd)
  ! 2- coupler 2 : lnd ---> atm (cpl_lnd2atm)
  !-----------------------------------------------------------------------

  use ESMF
  implicit none

  include 'mpif.h' !TODO: remove this and use ESMF

  private

  public :: cpl_atm2lnd_register
  public :: cpl_lnd2atm_register

  type(ESMF_RouteHandle) :: rh_atm2lnd
  type(ESMF_RouteHandle) :: rh_lnd2atm
  integer                :: mytask

  character(*), parameter :: modname =  "lilac_cpl"

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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    print *,'mytask= ',mytask
    if (mytask == 0) then
       print *, "in cpl_atm2lnd_register routine"
    end if

    ! Register the callback routines.
    ! Set the entry points for coupler ESMF Component methods
    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, userRoutine= cpl_atm2lnd_init, rc=rc)
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
    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_INITIALIZE, cpl_lnd2atm_init, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_RUN       , userRoutine=cpl_lnd2atm_run  , rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_CplCompSetEntryPoint(cplcomp, ESMF_METHOD_FINALIZE  , userRoutine=cpl_lnd2atm_final, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
  end subroutine cpl_lnd2atm_register

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
    integer                         :: n
    integer                         :: fieldcount
    character(len=128), allocatable :: fieldlist(:)
    character(len=128)              :: cvalue
    character(len=*), parameter     :: subname=trim(modname) //': [cpl_atm2lnd_init] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Coupler for atmosphere to land initialize routine called"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(importState, "a2c_fb", import_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldBundleGet(import_fieldbundle, fieldCount=fieldCount, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write(cvalue,*) fieldcount
    call ESMF_LogWrite(subname//" a2c_fb field count = "//trim(cvalue), ESMF_LOGMSG_INFO)
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(import_fieldbundle, fieldNameList=fieldlist, rc=rc) 
    do n = 1,fieldCount
       write(cvalue,*) n 
       call ESMF_LogWrite(subname//" a2c_fb field "//trim(cvalue)//' = '//trim(fieldlist(n)), ESMF_LOGMSG_INFO)
    end do
    deallocate(fieldlist)
    if (mytask == 0) then
       print *, ' a2c_fb field count = ',fieldcount
    end if

    call ESMF_StateGet(exportState, "c2l_fb", export_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldBundleGet(export_fieldbundle, fieldCount=fieldCount, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write(cvalue,*) fieldcount
    call ESMF_LogWrite(subname//" c2l_fb field count = "//trim(cvalue), ESMF_LOGMSG_INFO)
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(export_fieldbundle, fieldNameList=fieldlist, rc=rc) 
    do n = 1,fieldCount
       write(cvalue,*) n 
       call ESMF_LogWrite(subname//" c2l_fb field "//trim(cvalue)//' = '//trim(fieldlist(n)), ESMF_LOGMSG_INFO)
    end do
    deallocate(fieldlist)
    if (mytask == 0) then
       print *, ' c2l_fb field count = ',fieldcount
    end if

    if (mytask == 0) then
       print *, "PRINTING FIELDBUNDLES from atm->lnd"
       call ESMF_FieldBundlePrint (import_fieldbundle, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_FieldBundlePrint (export_fieldbundle, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    end if

    call ESMF_FieldBundleRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(subname//"cpl init finished!", ESMF_LOGMSG_INFO)

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
    integer                         :: n
    integer                         :: fieldcount
    character(len=128), allocatable :: fieldlist(:)
    character(len=128)              :: cvalue
    character(len=*) , parameter    :: subname=trim(modname ) //': [cpl_lnd2atm_init] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS

    if (mytask == 0) then
       print *, "Coupler for land to atmosphere initialize routine called"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(importState, "l2c_fb", import_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldBundleGet(import_fieldbundle, fieldCount=fieldCount, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write(cvalue,*) fieldcount
    call ESMF_LogWrite(subname//" l2c_fb field count = "//trim(cvalue), ESMF_LOGMSG_INFO)
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(import_fieldbundle, fieldNameList=fieldlist, rc=rc) 
    do n = 1,fieldCount
       write(cvalue,*) n 
       call ESMF_LogWrite(subname//" l2c_fb field "//trim(cvalue)//' = '//trim(fieldlist(n)), ESMF_LOGMSG_INFO)
    end do
    deallocate(fieldlist)
    if (mytask == 0) then
       print *, ' l2c_fb field count = ',fieldcount
    end if

    call ESMF_StateGet(exportState, "c2a_fb", export_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldBundleGet(export_fieldbundle, fieldCount=fieldCount, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write(cvalue,*) fieldcount
    call ESMF_LogWrite(subname//" c2a_fb field count = "//trim(cvalue), ESMF_LOGMSG_INFO)
    allocate(fieldlist(fieldcount))
    call ESMF_FieldBundleGet(export_fieldbundle, fieldNameList=fieldlist, rc=rc) 
    do n = 1,fieldCount
       write(cvalue,*) n 
       call ESMF_LogWrite(subname//" c2a_fb field "//trim(cvalue)//' = '//trim(fieldlist(n)), ESMF_LOGMSG_INFO)
    end do
    deallocate(fieldlist)
    if (mytask == 0) then
       print *, ' c2a_fb field count = ',fieldcount
    end if

    call ESMF_FieldBundleRedistStore(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(subname//"cpl init finished!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2atm_init

!======================================================================

  subroutine cpl_atm2lnd_run(cplcomp, importState, exportState, clock, rc)

    type(ESMF_CplComp)   :: cplcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle ) :: import_fieldbundle, export_fieldbundle
    character(len=*        ) , parameter :: subname=trim(modname       ) //': [cpl_atm2lnd_run] '
    !---------------------------------------------------

    rc = ESMF_SUCCESS
    if (mytask == 0) then
       print *, "Running cpl_atm2lnd_run"
    end if
    call ESMF_LogWrite(subname//"-----------------!", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(importState, trim("a2c_fb"), import_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//" got a2c fieldbundle!", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(exportState, trim("c2l_fb"), export_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite(subname//" got c2l fieldbundle!", ESMF_LOGMSG_INFO)

    call ESMF_FieldBundleRedist(import_fieldbundle, export_fieldbundle, routehandle=rh_atm2lnd, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
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

    call ESMF_StateGet(importState, "l2c_fb", import_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(exportState, "c2a_fb", export_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldBundleRedist(import_fieldbundle, export_fieldbundle, routehandle=rh_lnd2atm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(subname//" regridding fieldbundles  from land to atmos!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2atm_run

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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(subname//" rh_lnd2atm route handle released!", ESMF_LOGMSG_INFO)

  end subroutine cpl_lnd2atm_final

end module lilac_cpl
