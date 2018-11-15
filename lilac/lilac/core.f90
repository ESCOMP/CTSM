module lilac

  use ESMF
  use esmf_utils

  implicit none

  character(*), parameter :: modname =  "(core)"
  integer, parameter :: LILAC_SUCCESS = ESMF_SUCCESS

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  public :: init
  public :: run
  public :: final

  private  :: atmos_register 
  private  :: land_register
  private  :: cpl_register

  type, public :: LilacType
     private

     type(ESMFInfoType)             :: esmf_info
     character(len=ESMF_MAXSTR)     :: name

   contains
     procedure, public  :: init  => init
     procedure, public  :: run => run
     procedure, public  :: final => final

     ! register methods
     procedure, nopass, private  :: atmos_register  => atmos_register
     procedure, nopass, private  :: land_register => land_register
     procedure, nopass, private  :: cpl_register => cpl_register

     ! Init methods
     procedure, nopass, private  :: atmos_init  => atmos_init
     procedure, nopass, private  :: land_init  => land_init
     procedure, nopass, private  :: coupler_init  => coupler_init

     ! Run methods
     procedure, nopass, private  :: atmos_copy_atm_to_lilac => atmos_copy_atm_to_lilac
     procedure, nopass, private  :: atmos_copy_lilac_to_atm => atmos_copy_lilac_to_atm
     procedure, nopass, private  :: land_run => land_run
     procedure, nopass, private  :: coupler_run => coupler_run

     ! Final methods
     procedure, nopass, private  :: atmos_final => atmos_final
     procedure, nopass, private  :: land_final => land_final
     procedure, nopass, private  :: coupler_final => coupler_final

  end type LilacType

contains

  subroutine init(self, name)
    implicit none
    class(LilacType), intent(inout)  :: self
    character(len=ESMF_MAXSTR), intent(in)     :: name

    character(len=*), parameter :: subname=trim(modname)//':(init) '

    call ESMF_LogWrite(subname//"Initializing lilac", ESMF_LOGMSG_INFO)

    self%name = trim(name)

    ! Initialize ESMF structures
    call self%esmf_info%init(name, atmos_register, land_register, cpl_register)

  end subroutine init

  subroutine run(self)
    implicit none
    class(LilacType), intent(inout)  :: self

    character(len=*), parameter :: subname=trim(modname)//':(run) '

    call ESMF_LogWrite(subname//"Running lilac", ESMF_LOGMSG_INFO)

    call self%esmf_info%run()

  end subroutine run

  subroutine final(self)
    implicit none
    class(LilacType), intent(inout)  :: self

    character(len=*), parameter :: subname=trim(modname)//':(final) '

    call ESMF_LogWrite(subname//"Finalizing lilac", ESMF_LOGMSG_INFO)

    call self%esmf_info%final()

  end subroutine final

  subroutine atmos_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional
    character(len=*), parameter :: subname=trim(modname)//':(atmos_register) '

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
         userRoutine=atmos_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         userRoutine=atmos_copy_atm_to_lilac, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         userRoutine=atmos_copy_lilac_to_atm, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
         userRoutine=atmos_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    rc = ESMF_SUCCESS

  end subroutine atmos_register

  subroutine land_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional
    character(len=*), parameter :: subname=trim(modname)//':(lnd_register) '

    ! land_* comes from ctsm esmf cap

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=land_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=land_run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=land_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    rc = ESMF_SUCCESS

  end subroutine land_register

  subroutine cpl_register(comp, rc)
    type(ESMF_CplComp)   :: comp   ! must not be optional
    integer, intent(out)  :: rc     ! must not be optional
    character(len=*), parameter :: subname=trim(modname)//':(cpl_register) '

    rc = ESMF_FAILURE

    ! Register the callback routines.
    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=coupler_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=coupler_run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=coupler_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_LogWrite(subname//"CouplerMod: Registered Initialize, Run, and Finalize routines", ESMF_LOGMSG_INFO)

    rc = ESMF_SUCCESS

  end subroutine cpl_register

  subroutine atmos_init(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(atmos_init) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"atmos_init has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine atmos_init

  subroutine land_init(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(land_init) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"land_init has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine land_init

  subroutine coupler_init(comp, importState, exportState, clock, rc)
    type(ESMF_CplComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(coupler_init) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"coupler_init has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine coupler_init

  subroutine atmos_copy_atm_to_lilac(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(atmos_copy_atm_to_lilac) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"atmos_copy_atm_to_lilac has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine atmos_copy_atm_to_lilac

  subroutine atmos_copy_lilac_to_atm(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(atmos_copy_lilac_to_atm) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"atmos_copy_lilac_to_atm has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine atmos_copy_lilac_to_atm

  subroutine land_run(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(land_run) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"land_run has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine land_run

  subroutine coupler_run(comp, importState, exportState, clock, rc)
    type(ESMF_CplComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(coupler_run) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"coupler_run has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine coupler_run

  subroutine atmos_final(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(atmos_final) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"atmos_final has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine atmos_final

  subroutine land_final(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(land_final) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"land_final has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine land_final

  subroutine coupler_final(comp, importState, exportState, clock, rc)
    type(ESMF_CplComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(coupler_final) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"coupler_final has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine coupler_final

end module lilac
