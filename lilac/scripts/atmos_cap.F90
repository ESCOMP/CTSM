module atmos_cap

  use ESMF
  use lilac_utils


  implicit none

  character(*), parameter :: modname =  "(atmos_cap)"

  !!integer, parameter      :: fldsMax = 100

  type(ESMF_Field), public, save :: field
  type(ESMF_Field), public, save :: field_sie, field_u

  type(fld_list_type), allocatable :: x2a_fields(:)  
  type(fld_list_type), allocatable :: a2x_fields(:)  

  !private

  public atmos_register
  !public  :: add_fields
  !public  :: import_fields
  !public  :: export_fields

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
    type (ESMF_GridComp)     :: comp
    type (ESMF_State)        :: importState, exportState
    type (ESMF_Clock)        :: clock
    integer, intent(out)     :: rc

    ! local variables
    !!! TODO: Maybe it is better to call these fldsToAtm and fldsFrAtm
    type (fld_list_type)       :: fldsToCpl(fldsMax)
    type (fld_list_type)       :: fldsFrCpl(fldsMax)
    integer                    :: fldsToCpl_num
    integer                    :: fldsFrCpl_num
    type (ESMF_FieldBundle)    :: FBout
    integer                    :: n
    type(ESMF_Mesh)            :: atmos_mesh
    character(len=ESMF_MAXSTR) :: atmos_mesh_filepath
    integer                    :: petCount, localrc, urc
    integer                    :: mid, by2, quart, by4
    type(ESMF_Grid)            :: atmos_grid
    type(ESMF_DistGrid)        :: distgridIN, distgridFS
    logical                    :: mesh_switch
    character(len=*), parameter :: subname=trim(modname)//':(atmos_init) '
    !----------------------

    !integer                    :: regDecomp(:,:)
    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(comp, petcount=petcount, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-------------------------------------------------------------------------
    !    Generate -- Read in  the mesh
    !-------------------------------------------------------------------------

    mesh_switch = .false.
    if(mesh_switch) then
        ! For now this is our dummy mesh: 
        atmos_mesh_filepath='/gpfs/fs1/p/cesmdata/cseg/inputdata/share/meshes/T31_040122_ESMFmesh.nc'

        atmos_mesh = ESMF_MeshCreate(filename=trim(atmos_mesh_filepath), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Mesh for atmosphere is created!", ESMF_LOGMSG_INFO)
        print *, "!Mesh for atmosphere is created!"

    else
        !Grid1= ESMF_GridCreateNoPeriDimUfrmR( maxIndex=(/180,360 /), &
        !      minCornerCoord=(/0._ESMF_KIND_R8, 0._ESMF_KIND_R8/), &
        !      maxCornerCoord=(/180._ESMF_KIND_R8, 360._ESMF_KIND_R8/), &
        !      regDecomp=(/petcount,1/), rc=rc)

        atmos_grid = ESMF_GridCreateNoPeriDimUfrm( minIndex= (/1,1/), maxIndex=(/180,360 /), &
              maxCornerCoord=(/180._ESMF_KIND_R8, 360._ESMF_KIND_R8/), &
              minCornerCoord=(/0._ESMF_KIND_R8, 0._ESMF_KIND_R8/), &
              coordSys=ESMF_COORDSYS_CART,&
              regDecomp=(/1,petcount/),&
                                rc=rc)

        call ESMF_LogWrite(subname//"Grid for atmosphere is created!", ESMF_LOGMSG_INFO)
        print *, "Grid for atmosphere is created!"
    endif

    !-------------------------------------------------------------------------
    ! Coupler (land) to Atmosphere Fields --  x2a
    ! I- Create Field Bundle -- FBout for now-- TODO: negin want to rename to x2a_fieldbundle
    ! II- Create  Fields and add them to field bundle 
    ! III - Add FBout to state (x2a_state)
    !-------------------------------------------------------------------------
    FBout = ESMF_FieldBundleCreate(name="x2a_fields", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    ! Create individual fields and add to field bundle
    do n = 1,fldsFrCpl_num
       ! create field
       !!! Here we want to pass pointers
       !field = ESMF_FieldCreate(Emesh,farrayPtr=dum_var1_ptr,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldsFrCpl(n)%stdname), rc=rc)
       !field = ESMF_FieldCreate(lmesh,farrayPtr=x2a_fields%fields(:, n),  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldsFrCpl(n)%stdname), rc=rc)
       print *, trim(fldsFrCpl(n)%stdname)
       if (mesh_switch) then
          field = ESMF_FieldCreate(atmos_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(fldsFrCpl(n)%stdname), farrayPtr=fldsFrCpl(n)%farrayptr1d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
      else
          field = ESMF_FieldCreate(atmos_grid, name=trim(fldsFrCpl(n)%stdname), farrayPtr=fldsFrCpl(n)%farrayptr2d, rc=rc)
          !field = ESMF_FieldCreate(atmos_mesh, name=trim(fldsFrCpl(n)%stdname), farrayPtr=fldsFrCpl(n)%farrayptr1d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
       end if
       ! add field to field bundle
       call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
       enddo
    print *, "!Fields For Coupler (fldsFrCpl) Field Bundle Created!"

    ! Add FB to state
    call ESMF_StateAdd(exportState, (/FBout/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    ! Atmosphere to Coupler Fields
    FBout = ESMF_FieldBundleCreate(name="a2x_fields", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    ! Create individual states and add to field bundle
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------

    !call fldlist_add(fldsToCpl_num, fldsToCpl, 'dum_var2'      )
    do n = 1,fldsToCpl_num
       ! create field
       !field = ESMF_FieldCreate(lmesh, farrayPtr=a2x_field%fields(:,n) , meshloc=ESMF_MESHLOC_ELEMENT, name=trim(fldsToCpl(n)%stdname), rc=rc)
       field = ESMF_FieldCreate(atmos_mesh, ESMF_TYPEKIND_R8 ,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(fldsToCpl(n)%stdname), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
       ! initialize with default value
       !call ESMF_FieldGet(field, farrayPtr=fldptr, rc=rc)
       !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
       !fldptr = fldsToCpl(n)%default_value

       ! add field to field bundle
       call ESMF_FieldBundleAdd(FBout, (/field/), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
    enddo


    print *, "!Fields to  Coupler (fldstoCpl) Field Bundle Created!"

    ! Add FB to state
    call ESMF_StateAdd(importState, (/FBout/), rc=rc)
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





end module atmos_cap
