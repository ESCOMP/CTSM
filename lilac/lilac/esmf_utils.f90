module esmf_utils

  ! Wrappers and derived types exposing ESMF components to LILAC


#include "ESMF.h"
#include <lilac.h>
  use ESMF

  implicit none
  private

  character(*), parameter :: modname =  "(esmf_utils)"

  ! Consider renaming ESMFInfoType (add lilac to name)
  type, public :: ESMFInfoType
     private
     character(len=ESMF_MAXSTR) :: name

     type(ESMF_VM)                :: vm
     type(ESMF_State)             :: land_import
     type(ESMF_State)             :: land_export
     type(ESMF_State)             :: atmos_import
     type(ESMF_State)             :: atmos_export
     type(ESMF_GridComp)          :: atmos_comp
     type(ESMF_GridComp)          :: land_comp
     type(ESMF_CplComp)           :: cpl_comp

   contains
     procedure, public   :: init      => init
     procedure, public   :: run       => run
     procedure, public   :: final     => final

  end type ESMFInfoType

contains

  subroutine init(self, name)
    implicit none
    class(ESMFInfoType), intent(inout)  :: self
    character(len=ESMF_MAXSTR), intent(in) :: name

    ! TODO define subroutines: https://stackoverflow.com/questions/32809769/how-to-pass-subroutine-names-as-arguments-in-fortran

    ! Local variables
    character(len=ESMF_MAXSTR)            :: cname1, cname2, cplname
    integer                               :: localPet, petCount, localrc, rc=ESMF_SUCCESS, userrc=ESMF_SUCCESS

    character(len=*) :: subname=trim(modname)//':(init) '

    call ESMF_LogWrite(subname//"esmf_info%init()", ESMF_LOGMSG_INFO)

    self%name = trim(name)

    ! Create section
    !-------------------------------------------------------------------------

    ! Initialize framework and get back default global VM

    ! only run if not esmf_isintialized()
    call ESMF_Initialize(vm=self%vm, defaultlogfilename="lilac.log", logkindflag=ESMF_LOGKIND_MULTI, rc=localrc)
    if (return_error(localrc, rc)) return

    ! Get number of PETs we are running with
    call ESMF_VMGet(self%vm, petCount=petCount, localPet=localPet, rc=localrc)
    if (return_error(localrc, rc)) return

    ! Create the 2 model components and a coupler
    cname1 = "land"
    ! use petList to define land on all PET
    self%land_comp = ESMF_GridCompCreate(name=cname1, rc=localrc)
    call ESMF_LogWrite(subname//"Created "//trim(cname1)//" component", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return

    cname2 = "atmosphere"
    ! use petList to define atmosphere on all PET
    self%atmos_comp = ESMF_GridCompCreate(name=cname2, rc=localrc)
    call ESMF_LogWrite(subname//"Created "//trim(cname2)//" component", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return

    cplname = "lilac coupler"
    ! no petList means that coupler component runs on all PETs
    self%cpl_comp = ESMF_CplCompCreate(name=cplname, rc=localrc)
    call ESMF_LogWrite(subname//"Created "//trim(cplname)//" component", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return

    call ESMF_LogWrite(subname//"Comp Creates finished", ESMF_LOGMSG_INFO)

    ! Register section
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(self%atmos_comp, userRoutine=atmos_register, userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"atmos SetServices finished", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

    call ESMF_GridCompSetServices(self%land_comp, userRoutine=land_register, userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"land SetServices finished", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

    call ESMF_CplCompSetServices(self%cpl_comp, userRoutine=cpl_register, userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"Cpl SetServices finished", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

    ! Init section
    !-------------------------------------------------------------------------
    ! land import/export states
    self%land_import = ESMF_StateCreate(name="land import", stateintent=ESMF_STATEINTENT_IMPORT, rc=localrc)
    if (return_error(localrc, rc)) return
    self%land_export = ESMF_StateCreate(name="land export", stateintent=ESMF_STATEINTENT_EXPORT, rc=localrc)
    if (return_error(localrc, rc)) return
    call ESMF_GridCompInitialize(self%land_comp, importState=self%land_import, exportState=self%land_export, userRc=userrc, rc=localrc)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return
    call ESMF_LogWrite(subname//"Land Initialize finished", ESMF_LOGMSG_INFO)

    ! atmosphere import/export state
    self%atmos_import = ESMF_StateCreate(name="atmos import",  &
         stateintent=ESMF_STATEINTENT_IMPORT, rc=localrc)
    if (return_error(localrc, rc)) return

    self%atmos_export = ESMF_StateCreate(name="atmos export",  &
         stateintent=ESMF_STATEINTENT_EXPORT, rc=localrc)
    if (return_error(localrc, rc)) return
    call ESMF_GridCompInitialize(self%atmos_comp, exportState=self%atmos_export, userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"Atmosphere Initialize finished", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

    ! call ESMF_CPLCompInitialize twice (once for each grid comp)

  end subroutine init

  subroutine run(self)
    implicit none
    class(ESMFInfoType), intent(inout)  :: self
    integer :: localrc, rc=ESMF_SUCCESS, userrc=ESMF_SUCCESS
    character(len=*), parameter :: subname=trim(modname)//':(init) '

    call ESMF_LogWrite(subname//"esmf_info%run()", ESMF_LOGMSG_INFO)

    ! TODO: need some help on order of imports/exports/runs and whether the land/atm both need import/export states

    ! atmosphere run
    ! copy the atmos state and put it into atmos export
    call ESMF_GridCompRun(self%atmos_comp, exportState=self%atmos_export, phase=1, userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"Atmosphere Run returned", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

    ! coupler run
    call ESMF_CplCompRun(self%cpl_comp, importState=self%atmos_export, exportState=self%land_import, &
         userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"Coupler Run returned", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

    ! land run
    call ESMF_GridCompRun(self%land_comp, importState=self%land_import, exportState=self%land_export, userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"Land Run returned", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

    ! coupler run
    call ESMF_CplCompRun(self%cpl_comp, importState=self%land_export, exportState=self%atmos_import, &
         userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"Coupler Run returned", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

    call ESMF_GridCompRun(self%atmos_comp, importState=self%atmos_import, phase=2, userRc=userrc, rc=localrc)
    call ESMF_LogWrite(subname//"Atmosphere Run returned", ESMF_LOGMSG_INFO)
    if (return_error(localrc, rc)) return
    if (return_error(userrc, rc)) return

  end subroutine run

  subroutine final(self)
    implicit none
    class(ESMFInfoType), intent(inout)  :: self
    integer :: localrc, rc=ESMF_SUCCESS
    character(len=*), parameter :: subname=trim(modname)//':(final) '

    call ESMF_LogWrite(subname//"esmf_info%final()", ESMF_LOGMSG_INFO)

    ! Destroy section
    call ESMF_GridCompDestroy(self%atmos_comp, rc=localrc)
    if (return_error(localrc, rc)) return
    call ESMF_GridCompDestroy(self%land_comp, rc=localrc)
    if (return_error(localrc, rc)) return
    call ESMF_CplCompDestroy(self%cpl_comp, rc=localrc)
    if (return_error(localrc, rc)) return

    call ESMF_StateDestroy(self%land_export, rc=localrc)
    call ESMF_StateDestroy(self%land_import, rc=localrc)
    if (return_error(localrc, rc)) return
    call ESMF_StateDestroy(self%atmos_export, rc=localrc)
    call ESMF_StateDestroy(self%atmos_import, rc=localrc)
    ! do this everywhere
    if (return_error(localrc, rc)) return

    call ESMF_LogWrite(subname//"All Destroy routines done", ESMF_LOGMSG_INFO)

  end subroutine final

  subroutine atmos_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional
    character(len=*), parameter :: subname=trim(modname)//':(atmos_register) '

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
         userRoutine=atoms_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         userRoutine=atoms_copy_atm_to_lilac, phase=1, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         userRoutine=atoms_copy_lilac_to_atm, phase=2, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
         userRoutine=atoms_final, rc=rc)
    ! TODO: check rcs

    rc = ESMF_SUCCESS

  end subroutine atmos_register

  subroutine land_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional
    character(len=*), parameter :: subname=trim(modname)//':(lnd_register) '

    ! land_* comes from ctsm esmf cap

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=land_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=land_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=land_final, rc=rc)
    ! TODO: check rcs

    rc = ESMF_SUCCESS

  end subroutine land_register

  subroutine cpl_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional
    character(len=*), parameter :: subname=trim(modname)//':(cpl_register) '

    rc = ESMF_FAILURE

    ! Register the callback routines.

    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=coupler_init, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=coupler_run, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)
    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=coupler_final, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rc)

    call ESMF_LogWrite(subname//"CouplerMod: Registered Initialize, Run, and Finalize routines", ESMF_LOGMSG_INFO)

    rc = ESMF_SUCCESS

  end subroutine cpl_register

  function return_error(rc, returnrc) result(error)
    ! fight with this later
    integer, intent(in) :: rc, returnrc
    logical             :: error
    if (ESMF_LogFoundError(rc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=returnrc)) then
       error = .true.
    else
       error = .false.
    endif

  end function return_error

end module esmf_utils
