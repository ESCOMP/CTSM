module atmos_cap

    use ESMF
    use lilac_utils
    !use lilac_mod,   only :    a2l_fields


    implicit none

    character(*), parameter          ::  modname =  "atmos_cap"

    !!integer, parameter              :: fldsMax = 100

    type(ESMF_Field), public, save   ::  field
    type(ESMF_Field), public, save   ::  field_sie, field_u

    !type(fld_list_type), public, allocatable ::  l2a_fields(:)
    !type(fld_list_type), public, allocatable ::  a2l_fields(:)
    !type(fld_list_type), public, save ::  l2a_fields(:)
    !type(fld_list_type), public, save ::  a2l_fields(:)
    type(fld_list_type), allocatable ::  l2a_fields(:)
    type(fld_list_type), allocatable ::  a2l_fields(:)

    !private

    public  :: atmos_register
    !public  :: add_fields
    !public  :: import_fields
    !public  :: export_fields


    !------------------------------------------------------------------------

    contains

    subroutine atmos_register (comp, rc)

        type(ESMF_GridComp)          :: comp   ! must not be optional
        integer, intent(out)         :: rc
        character(len=*), parameter  :: subname=trim(modname)//':(atmos_register) '

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



    subroutine atmos_init (comp, importState, exportState, clock, rc)

        type (ESMF_GridComp)         ::  comp
        type (ESMF_State)            ::  importState, exportState
        type (ESMF_Clock)            ::  clock
        integer, intent(out)         ::  rc

        ! local variables
        !!! TODO: Maybe it is better to call these fldsToAtm and fldsFrAtm
        type (fld_list_type)         ::  flds_a2l(fldsMax)
        type (fld_list_type)         ::  flds_l2a(fldsMax)
        integer                      ::  flds_a2l_num
        integer                      ::  flds_l2a_num
        type (ESMF_FieldBundle)      ::  l2a_fb , a2l_fb
        integer                      ::  n
        type(ESMF_Mesh)              ::  atmos_mesh
        character(len=ESMF_MAXSTR)   ::  atmos_mesh_filepath
        integer                      ::  petCount, localrc, urc
        integer                      ::  mid, by2, quart, by4
        type(ESMF_Grid)              ::  atmos_grid
        type(ESMF_DistGrid)          ::  distgridIN, distgridFS
        logical                      ::  mesh_switch
        character(len=*), parameter  :: subname=trim(modname)//':[atmos_init] '
        !----------------------

        !integer                    :: regDecomp(:,:)
        ! Initialize return code
        rc = ESMF_SUCCESS

        call ESMF_GridCompGet (comp, petcount=petcount, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

        !-------------------------------------------------------------------------
        !    Generate -- Read in  the mesh
        !-------------------------------------------------------------------------
        mesh_switch = .True.

        if(mesh_switch) then
            ! For now this is our dummy mesh: 
            atmos_mesh_filepath  =   '/gpfs/fs1/p/cesmdata/cseg/inputdata/share/meshes/T31_040122_ESMFmesh.nc'

            atmos_mesh           = ESMF_MeshCreate(filename=trim(atmos_mesh_filepath), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
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
        ! Coupler (land) to Atmosphere Fields --  l2a
        ! I- Create Field Bundle -- l2a_fb for now 
        ! II- Create  Fields and add them to field bundle 
        ! III - Add l2a_fb to state (l2a_state)
        !-------------------------------------------------------------------------

        l2a_fb = ESMF_FieldBundleCreate (name="l2a_fields", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        ! Create individual fields and add to field bundle
        flds_l2a_num = 2
        do n = 1,flds_l2a_num
           ! create field
           !!! Here we want to pass pointers
           print *, trim(flds_l2a(n)%stdname)
           if (mesh_switch) then
              field = ESMF_FieldCreate(atmos_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(flds_l2a(n)%stdname), farrayPtr=flds_l2a(n)%farrayptr1d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
           else
              field = ESMF_FieldCreate(atmos_grid, name=trim(flds_l2a(n)%stdname), farrayPtr=flds_l2a(n)%farrayptr2d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
           end if
           ! add field to field bundle
           call ESMF_FieldBundleAdd(l2a_fb, (/field/), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        enddo
        print *, "!Fields For Coupler (flds_l2a) Field Bundle Created!"

        ! Add FB to state
        call ESMF_StateAdd(exportState, (/l2a_fb/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out


        ! Atmosphere to Coupler Fields
        a2l_fb = ESMF_FieldBundleCreate(name="a2l_fields", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        !-------------------------------------------------------------------------
        ! Create individual states and add to field bundle -- a2l
        !-------------------------------------------------------------------------

        !call fldlist_add(flds_a2l_num, flds_a2l, 'dum_var2'      )

        do n = 1,flds_a2l_num
           ! create field
           field = ESMF_FieldCreate(atmos_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(flds_a2l(n)%stdname), farrayPtr=flds_a2l(n)%farrayptr1d, rc=rc)
           !field = ESMF_FieldCreate(atmos_mesh, ESMF_TYPEKIND_R8 ,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(flds_a2l(n)%stdname), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
           ! initialize with default value
           !call ESMF_FieldGet(field, farrayPtr=fldptr, rc=rc)
           !fldptr = flds_a2l(n)%default_value

           ! add field to field bundle
           call ESMF_FieldBundleAdd(a2l_fb, (/field/), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        enddo


        print *, "!Fields to  Coupler (flds_a2l) Field Bundle Created!"

        ! Add FB to state
        call ESMF_StateAdd(importState, (/a2l_fb/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        print *, "!a2l_state is filld with dummy_var field bundle!"
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

end module atmos_cap
