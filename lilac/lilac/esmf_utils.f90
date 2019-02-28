module esmf_utils

  ! Wrappers and derived types exposing ESMF components to LILAC

#include "ESMF.h"
  use ESMF

  implicit none
  private

  character(*), parameter :: modname =  "(esmf_utils)"

  interface
     subroutine userRoutine(gridcomp, rc)
       use ESMF_CompMod
       implicit none
       type(ESMF_GridComp)        :: gridcomp ! must not be optional
       integer, intent(out)       :: rc       ! must not be optional
     end subroutine userRoutine
  end interface

  interface
     subroutine userCplRoutine(cplcomp, rc)
       use ESMF_CompMod
       implicit none
       type(ESMF_CplComp)         :: cplcomp  ! must not be optional
       integer, intent(out)       :: rc       ! must not be optional
     end subroutine userCplRoutine
  end interface

  ! Consider renaming ESMFInfoType (add lilac to name)
  type, public :: ESMFInfoType
     private

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

  subroutine start(self, rc)
    implicit none
    integer, intent(in) :: rc=ESMF_SUCCESS

    integer             :: localPet, petCount

    character(len=*), parameter :: subname=trim(modname)//':(start) '

    call ESMF_LogWrite(subname//"esmf_info%start()", ESMF_LOGMSG_INFO)

    ! Initialize framework and get back default global VM

    ! only run if not esmf_isintialized()
    call ESMF_Initialize(vm=self%vm, defaultlogfilename="lilac.log", logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    ! Get number of PETs we are running with
    call ESMF_VMGet(self%vm, petCount=petCount, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

  end subroutine start

  subroutine init(self, atmos_register, land_register, cpl_register, rc)
    implicit none
    class(ESMFInfoType), intent(inout)         :: self
    procedure(userRoutine)                     :: atmos_register
    procedure(userRoutine)                     :: land_register
    procedure(userCplRoutine)                  :: cpl_register
    integer, intent(in)                        :: rc=ESMF_SUCCESS

    ! Local variables
    character(len=ESMF_MAXSTR)            :: cname1, cname2, cplname
    integer                               :: rc=ESMF_SUCCESS

    character(len=*), parameter :: subname=trim(modname)//':(init) '

    call ESMF_LogWrite(subname//"esmf_info%init()", ESMF_LOGMSG_INFO)

    ! Create section
    !-------------------------------------------------------------------------
    ! Create the 2 model components and a coupler
    cname1 = "land"
    ! use petList to define land on all PET
    self%land_comp = ESMF_GridCompCreate(name=cname1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname1)//" component", ESMF_LOGMSG_INFO)

    cname2 = "atmosphere"
    ! use petList to define atmosphere on all PET
    self%atmos_comp = ESMF_GridCompCreate(name=cname2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cname2)//" component", ESMF_LOGMSG_INFO)

    cplname = "lilac coupler"
    ! no petList means that coupler component runs on all PETs
    self%cpl_comp = ESMF_CplCompCreate(name=cplname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Created "//trim(cplname)//" component", ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//"Comp Creates finished", ESMF_LOGMSG_INFO)

    ! Register section
    !-------------------------------------------------------------------------
    call ESMF_GridCompSetServices(self%atmos_comp, userRoutine=atmos_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"atmos SetServices finished", ESMF_LOGMSG_INFO)

    call ESMF_GridCompSetServices(self%land_comp, userRoutine=land_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"land SetServices finished", ESMF_LOGMSG_INFO)

    call ESMF_CplCompSetServices(self%cpl_comp, userRoutine=cpl_register, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Cpl SetServices finished", ESMF_LOGMSG_INFO)

    ! Init section
    !-------------------------------------------------------------------------
    ! land import/export states
    self%land_import = ESMF_StateCreate(name="land import", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    self%land_export = ESMF_StateCreate(name="land export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompInitialize(self%land_comp, importState=self%land_import, exportState=self%land_export, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Land Initialize finished", ESMF_LOGMSG_INFO)

    ! atmosphere import/export state
    self%atmos_import = ESMF_StateCreate(name="atmos import",  stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    self%atmos_export = ESMF_StateCreate(name="atmos export", stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompInitialize(self%atmos_comp, exportState=self%atmos_export, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Atmosphere Initialize finished", ESMF_LOGMSG_INFO)

  end subroutine init

  subroutine run(self, rc)
    implicit none
    class(ESMFInfoType), intent(inout)  :: self
    integer :: rc=ESMF_SUCCESS
    character(len=*), parameter :: subname=trim(modname)//':(run) '

    call ESMF_LogWrite(subname//"esmf_info%run()", ESMF_LOGMSG_INFO)

    ! atmosphere run
    ! copy the atmos state and put it into atmos export
    call ESMF_GridCompRun(self%atmos_comp, exportState=self%atmos_export, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Atmosphere Run returned", ESMF_LOGMSG_INFO)

    ! coupler run
    call ESMF_CplCompRun(self%cpl_comp, importState=self%atmos_export, exportState=self%land_import, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Coupler Run returned", ESMF_LOGMSG_INFO)

    ! land run
    call ESMF_GridCompRun(self%land_comp, importState=self%land_import, exportState=self%land_export, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Land Run returned", ESMF_LOGMSG_INFO)

    ! coupler run
    call ESMF_CplCompRun(self%cpl_comp, importState=self%land_export, exportState=self%atmos_import, rc=rc)
    call ESMF_LogWrite(subname//"Coupler Run returned", ESMF_LOGMSG_INFO)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompRun(self%atmos_comp, importState=self%atmos_import, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Atmosphere Run returned", ESMF_LOGMSG_INFO)

  end subroutine run

  subroutine final(self, rc)
    implicit none
    class(ESMFInfoType), intent(inout)  :: self
    integer :: rc=ESMF_SUCCESS
    character(len=*), parameter :: subname=trim(modname)//':(final) '

    call ESMF_LogWrite(subname//"esmf_info%final()", ESMF_LOGMSG_INFO)

    ! Destroy section
    call ESMF_GridCompDestroy(self%atmos_comp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompDestroy(self%land_comp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_CplCompDestroy(self%cpl_comp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_StateDestroy(self%land_export, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_StateDestroy(self%land_import, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_StateDestroy(self%atmos_export, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_StateDestroy(self%atmos_import, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_LogWrite(subname//"All Destroy routines done", ESMF_LOGMSG_INFO)

  end subroutine final

end module esmf_utils
