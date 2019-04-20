module lnd_cap
  use ESMF
  use lilac_utils

  implicit none

  character(*), parameter :: modname =  "(land)"

  !!integer, parameter      :: fldsMax = 100

  type(ESMF_Field), public, save :: field
  type(ESMF_Field), public, save :: field_sie, field_u

  type(fld_list_type), allocatable :: x2a_fields(:)  
  type(fld_list_type), allocatable :: a2x_fields(:)  

  !private

  public lnd_register
  !public  :: add_fields
  !public  :: import_fields
  !public  :: export_fields

  contains

  subroutine lnd_register(comp, rc)

    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out) :: rc
    character(len=*), parameter :: subname=trim(modname)//':(atmos_register) '

    print *, "in user register routine"

    rc = ESMF_SUCCESS
    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=lnd_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=lnd_run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=lnd_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

  end subroutine lnd_register

  subroutine lnd_init(comp, importState, exportState, clock, rc)
    type (ESMF_GridComp)     :: comp
    type (ESMF_State)        :: importState, exportState
    type (ESMF_Clock)        :: clock
    integer, intent(out)     :: rc

    print *, "  Empty land is created !!!!"
    rc = ESMF_SUCCESS
    !-------------------------------------------------------------------------
    !    Generate -- Read in  the mesh
    !-------------------------------------------------------------------------
  end subroutine lnd_init


  subroutine lnd_run(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(atmos_copy_lilac_to_atm) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"atmos_copy_lilac_to_atm has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine lnd_run

  subroutine lnd_final(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(lnd_final) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"lnd_final is called but has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine lnd_final
  !===============================================================================





end module lnd_cap
