module DummyAtmos
  use ESMF
  use lilac_utils, only: fldlist_add, fld_list_type, fldsMax


  implicit none


  type(ESMF_Field), public, save :: field
  type(ESMF_Field), public, save :: field_sie, field_u

  integer                 :: flds_x2a_num = 0
  integer                 :: flds_a2x_num = 0


  type(fld_list_type), allocatable :: x2a_fields(:)  
  type(fld_list_type), allocatable :: a2x_fields(:)  

  !private
  character(*), parameter :: modname =  "(core)"

  public atmos_register
  !public  :: add_fields
  !public  :: import_fields
  !public  :: export_fields




  !!! Adding import export states stuff here....


  contains

  subroutine atmos_register(comp, rc)
    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out) :: rc
    character(len=*), parameter :: subname=trim(modname)//':(atmos_register) '
    
    print *, "in user register routine"
    
    rc = ESMF_SUCCESS
    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=atmos_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=atmos_copy_atm_to_lilac, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=atmos_copy_lilac_to_atm, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=atmos_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

  end subroutine atmos_register

  subroutine atmos_init(comp, importState, exportState, clock, rc) 
      !, dum_var1, dum_var2)

    type (ESMF_GridComp)     :: comp
    type (ESMF_State)        :: importState, exportState
    type (ESMF_Clock)        :: clock
    integer, intent(out)    :: rc
    
    !!! TODO: Maybe it is better to call these fldsToAtm and fldsFrAtm

    type (fld_list_type)    :: fldsToCpl(fldsMax)
    type (fld_list_type)    :: fldsFrCpl(fldsMax)
    integer                 :: fldsToCpl_num
    integer                :: fldsFrCpl_num

    character(len=*), parameter :: subname=trim(modname)//':(atmos_init) '

    type (ESMF_State)       :: x2a_state ! the coupled flow State
    type (ESMF_State)       :: a2x_state ! the coupled flow State
    type (ESMF_FieldBundle) :: FBout
    integer                 :: n
    
    type(ESMF_Mesh)      :: Emesh
    character(len=ESMF_MAXSTR)            :: atmos_mesh_filepath

    real, pointer :: dum_var1_ptr (:)
    real, pointer :: dum_var2_ptr (:)
    real, dimension(10) :: dum_var1
    real, dimension(10) :: dum_var2



    ! Initialize return code

    rc = ESMF_SUCCESS

    !-------------------------------------------------------------------------
    !    Generate -- Read in  the mesh
    !-------------------------------------------------------------------------

    ! For now this is our dummy mesh: 
    atmos_mesh_filepath='/gpfs/fs1/p/cesmdata/cseg/inputdata/share/meshes/T31_040122_ESMFmesh.nc'
    
    EMesh = ESMF_MeshCreate(filename=trim(atmos_mesh_filepath), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Mesh for atmosphere is created!", ESMF_LOGMSG_INFO)
    print *, "!Mesh for atmosphere is created!"

    !-------------------------------------------------------------------------
    ! Create States -- x2a_state (import) -- a2x_state (export)
    !-------------------------------------------------------------------------
    
    EMesh = ESMF_MeshCreate(filename=trim(atmos_mesh_filepath), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    call ESMF_LogWrite(subname//"Mesh for atmosphere is created!", ESMF_LOGMSG_INFO)
    print *, "!Mesh for atmosphere is created!"

    !-------------------------------------------------------------------------
    ! Create States -- x2a_state (import) -- a2x_state (export)
    !-------------------------------------------------------------------------
    x2a_state = ESMF_StateCreate(name="x2a_state", stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    a2x_state = ESMF_StateCreate(name="a2x_state",  stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    print *, "!empty x2a_state (import) is created!"
    print *, "!empty a2x_state (export) is created!"

    !-------------------------------------------------------------------------
    ! Create Field lists -- Basically create a list of fields and add a default
    ! value to them.
    !-------------------------------------------------------------------------

    !!! FOR NOW LET'S JUST ADD TWO THINGS....
    !!! WE WILL PUT THIS UNDER CREATE_FLDLIST LATER
    fldsFrCpl_num = 1
    fldsToCpl_num = 2

    call fldlist_add(fldsToCpl_num, fldsToCpl, 'dum_var1', default_value=30.0, units='m')
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'dum_var2', default_value=30.0, units='m')


    !-------------------------------------------------------------------------
    ! Coupler (land) to Atmosphere Fields --  x2a
    ! I- Create Field Bundle -- FBout for now-- TODO: negin want to rename to x2a_fieldbundle
    ! II- Create  Fields and add them to field bundle 
    ! III - Add FBout to state (x2a_state) 
    !-------------------------------------------------------------------------
    FBout = ESMF_FieldBundleCreate(name="x2a_fields", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out




    !!! THIS IS INIT! We don't Have anything from coupler

    ! Create individual states and add to field bundle
    fldsFrCpl_num = 0
    !call fldlist_add(fldsFrCpl_num, fldsFrCpl, 'dummy_var_1', default_value=0.0, units='m')
    do n = 1,fldsFrCpl_num
       ! create field
       !!! Here we want to pass pointers
       !!! Create With pointer ? or fieldcreate and then fieldget
       field = ESMF_FieldCreate(Emesh,farrayPtr=dum_var1_ptr,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldsFrCpl(n)%stdname), rc=rc)
       !field = ESMF_FieldCreate(lmesh,farrayPtr=x2a_fields%fields(:, n),  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldsFrCpl(n)%stdname), rc=rc)
       print *, trim(fldsFrCpl(n)%stdname)
       !field = ESMF_FieldCreate(EMesh, farrayPtr ,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldsFrCpl(n)%stdname), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
       ! add field to field bundle
       call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    enddo
    print *, "!Fields For Coupler (fldsFrCpl) Field Bundle Created!"

    ! Add FB to state
    call ESMF_StateAdd(x2a_state, (/FBout/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    ! Atmosphere to Coupler Fields
    FBout = ESMF_FieldBundleCreate(name="a2x_fields", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    ! Create individual states and add to field bundle
    fldsToCpl_num = 1
    call fldlist_add(fldsToCpl_num, fldsToCpl, 'dum_var2'      )
    do n = 1,fldsToCpl_num
       ! create field
       !field = ESMF_FieldCreate(lmesh, farrayPtr=a2x_field%fields(:,n) , meshloc=ESMF_MESHLOC_ELEMENT, name=trim(fldsToCpl(n)%stdname), rc=rc)
       field = ESMF_FieldCreate(EMesh, ESMF_TYPEKIND_R8 ,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldsToCpl(n)%stdname), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
       ! initialize with default value
       !call ESMF_FieldGet(field, farrayPtr=fldptr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    enddo
    call ESMF_StateAdd(a2x_state, (/FBout/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    print *, "!a2x_state is filld with dummy_var field bundle!"

  end subroutine atmos_init

  subroutine atmos_copy_atm_to_lilac(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(atmos_copy_atm_to_lilac) '

    ! Initialize return code
    rc = ESMF_SUCCESS
! get a list of fields of variables we need from atmos....
!
    !call ESMF_LogWrite(subname//"atmos_copy_atm_to_lilac has not been implemented yet", ESMF_LOGMSG_INFO)

    ! loop over fields, copying pointer from import to export state

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
  !===============================================================================

  ! Let's put this in a lilac_utils
  !subroutine fldlist_add(num, fldlist, stdname)
  !  integer,                    intent(inout) :: num
  !  type(fld_list_type),        intent(inout) :: fldlist(:)
  !  character(len=*),           intent(in)    :: stdname

    ! local variables
  !  integer :: rc
  !  integer :: dbrc
  !  character(len=*), parameter :: subname='(lnd_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

  !  num = num + 1
   ! if (num > fldsMax) then
   !    call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
   !         ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
   !    return
   ! endif
   ! fldlist(num)%stdname = trim(stdname)

  !end subroutine fldlist_add




end module DummyAtmos
