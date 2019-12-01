module lilac_atmcap

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This is a dummy atmosphere cap for setting up lilac structure.
  !-----------------------------------------------------------------------

  ! !USES
  use ESMF
  use lilac_utils , only : atm2lnd, lnd2atm, gindex_atm, atm_mesh_filename
  implicit none

  public :: lilac_atmos_register

  integer            :: mytask 
  integer, parameter :: debug = 0 ! internal debug level

!========================================================================
contains
!========================================================================

  subroutine lilac_atmos_register (comp, rc)

    type(ESMF_GridComp)          :: comp   ! must not be optional
    integer, intent(out)         :: rc

    ! local variables
    type(ESMF_VM)                :: vm
    character(len=*), parameter  :: subname='(lilac_atmos_register): '
    !-------------------------------------------------------------------------

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (mytask == 0) then
       print *, "in user register routine"
    end if

    ! Initialize return code
    rc = ESMF_SUCCESS

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=lilac_atmos_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=lilac_atmos_run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=lilac_atmos_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine lilac_atmos_register

!========================================================================

  subroutine lilac_atmos_init (comp, lnd2atm_a_state, atm2lnd_a_state, clock, rc)

    ! input/output variables
    type (ESMF_GridComp) :: comp
    type (ESMF_State)    :: lnd2atm_a_state, atm2lnd_a_state
    type (ESMF_Clock)    :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)             :: atm_mesh
    type(ESMF_DistGrid)         :: atm_distgrid
    type(ESMF_Field)            :: field
    type(ESMF_FieldBundle)      :: c2a_fb , a2c_fb
    integer                     :: n, i
    character(len=*), parameter :: subname='(lilac_atmos_init): '
    !-------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"------------------------!", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Read in the atm mesh
    !-------------------------------------------------------------------------

    ! Note that in the call to lilac_atm the host atmospere sent both the gindex_atm and 
    ! the atm_mesh_filename that were then set as module variables in lilac_utils

    atm_distgrid = ESMF_DistGridCreate (arbSeqIndexList=gindex_atm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    atm_mesh = ESMF_MeshCreate(filename=trim(atm_mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistGrid=atm_distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(trim(subname)//"Mesh for atmosphere is created for "//trim(atm_mesh_filename), ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // "Mesh for atmosphere is created for "//trim(atm_mesh_filename)
    end if

    !-------------------------------------------------------------------------
    ! Atmosphere to Coupler (land) Fields --  atmos --> land
    ! - Create empty field bundle -- a2c_fb
    ! - Create  Fields and add them to field bundle
    ! - Add a2c_fb to state (atm2lnd_a_state)
    !-------------------------------------------------------------------------

    ! Create individual fields and add to field bundle -- a2c
    a2c_fb = ESMF_FieldBundleCreate(name="a2c_fb", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! create fields and add to field bundle
    do n = 1, size(atm2lnd)
       field = ESMF_FieldCreate(atm_mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
            name=trim(atm2lnd(n)%fldname), farrayPtr=atm2lnd(n)%dataptr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_FieldBundleAdd(a2c_fb, (/field/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    end do

    call ESMF_LogWrite(subname//"fieldbundleadd is finished .... !", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "!Fields to  Coupler (atmos to  land ) (a2c_fb) Field Bundle Created!"
    end if

    ! Add field bundle to state
    call ESMF_StateAdd(atm2lnd_a_state, (/a2c_fb/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite(subname//"atm2lnd_a_state is filled with dummy_var field bundle!", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "!atm2lnd_a_state is filld with dummy_var field bundle!"
    end if

    !-------------------------------------------------------------------------
    ! Coupler (land) to Atmosphere Fields --  c2a
    ! - Create Field Bundle -- c2a_fb for because we are in atmos
    ! - Create  Fields and add them to field bundle
    ! - Add c2a_fb to state (lnd2atm_a_state)
    !-------------------------------------------------------------------------

    c2a_fb = ESMF_FieldBundleCreate (name="c2a_fb", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! create fields and add to field bundle
    do n = 1, size(lnd2atm)
       field = ESMF_FieldCreate(atm_mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
            name=trim(lnd2atm(n)%fldname), farrayPtr=lnd2atm(n)%dataptr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_FieldBundleAdd(c2a_fb, (/field/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       if (debug > 0) then
          call ESMF_FieldPrint(field,  rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       end if
    end do
    call ESMF_LogWrite(subname//"c2a fieldbundleadd is finished .... !", ESMF_LOGMSG_INFO)

    ! Add field bundle to state
    call ESMF_StateAdd(lnd2atm_a_state, (/c2a_fb/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Set Attributes needed by land
    call ESMF_AttributeSet(lnd2atm_a_state, name="nextsw_cday", value=11, rc=rc)  ! TODO: mv what in the world is this???

  end subroutine lilac_atmos_init

!========================================================================

  subroutine lilac_atmos_run(comp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(lilac_atmos_run):'

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"Should atmos_run ", ESMF_LOGMSG_INFO)

  end subroutine lilac_atmos_run

!========================================================================

  subroutine lilac_atmos_final(comp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle) ::  import_fieldbundle, export_fieldbundle
    character(len=*), parameter :: subname='( lilac_atmos_final): '
    !-------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_StateGet(importState, "c2a_fb", import_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(exportState, "a2c_fb", export_fieldbundle, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_FieldBundleDestroy(import_fieldbundle, rc=rc)
    call ESMF_FieldBundleDestroy(export_fieldbundle, rc=rc)

    call ESMF_LogWrite(subname//"?? Are there any other thing for destroying in atmos_final??", ESMF_LOGMSG_INFO)

  end subroutine lilac_atmos_final

end module lilac_atmcap
