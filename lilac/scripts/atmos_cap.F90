module atmos_cap

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:

    ! !USES
    use ESMF
    use lilac_utils

    implicit none

    character(*), parameter          ::  modname =  "atmos_cap"

    !!integer, parameter              :: fldsMax = 100

    type(ESMF_Field), public, save   ::  field

    type(fld_list_type), public, allocatable ::  c2a_fldlist(:)
    type(fld_list_type), public, allocatable ::  a2c_fldlist(:)

    !type (fld_list_type)            ::  a2c_fldlist(fldsMax)
    !type (fld_list_type)            ::  c2a_fldlist(fldsMax)

    integer                          ::  a2c_fldlist_num
    integer                          ::  c2a_fldlist_num

    !private

    public  :: atmos_register
    !real(kind=ESMF_KIND_R8), dimension(:), public, pointer, save :: fldptr

    !========================================================================
    contains
    !========================================================================

    subroutine atmos_register (comp, rc)

        type(ESMF_GridComp)          :: comp   ! must not be optional
        integer, intent(out)         :: rc
        character(len=*), parameter  :: subname=trim(modname)//':(atmos_register) '

        print *, "in user register routine"

        rc = ESMF_SUCCESS
        ! Set the entry points for standard ESMF Component methods
        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=atmos_init, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=atmos_run, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=atmos_final, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    end subroutine atmos_register



    subroutine atmos_init (comp, lnd2atm_a_state, atm2lnd_a_state, clock, rc)
        type (ESMF_GridComp)         ::  comp
        type (ESMF_State)            ::  lnd2atm_a_state, atm2lnd_a_state
        type (ESMF_Clock)            ::  clock
        integer, intent(out)         ::  rc

        ! local variables
        type (ESMF_FieldBundle)      ::  c2a_fb , a2c_fb
        integer                      ::  n
        type(ESMF_Mesh)              ::  atmos_mesh
        character(len=ESMF_MAXSTR)   ::  atmos_mesh_filepath
        integer                      ::  petCount, localrc, urc
        integer                      ::  mid, by2, quart, by4
        type(ESMF_Grid)              ::  atmos_grid
        type(ESMF_DistGrid)          ::  distgridIN, distgridFS
        logical                      ::  mesh_switch
        character(len=*), parameter  :: subname=trim(modname)//': [atmos_init] '
        !integer                     :: regDecomp(:,:)
        !-------------------------------------------------------------------------
        ! Initialize return code
        rc = ESMF_SUCCESS
        call ESMF_LogWrite(subname//"------------------------!", ESMF_LOGMSG_INFO)

        call ESMF_GridCompGet (comp, petcount=petcount, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

        !-------------------------------------------------------------------------
        !    Read in  the mesh ----or----- Generate the grid
        !-------------------------------------------------------------------------
        mesh_switch = .True.

        if(mesh_switch) then
            ! For now this is our dummy mesh: 
            !atmos_mesh_filepath  =   '/gpfs/fs1/p/cesmdata/cseg/inputdata/share/meshes/T31_040122_ESMFmesh.nc'  !! Negin: This did not work.... 
            atmos_mesh_filepath  =   '/gpfs/fs1/p/cesmdata/cseg/inputdata/share/meshes/fv1.9x2.5_141008_ESMFmesh.nc'

            atmos_mesh           = ESMF_MeshCreate(filename=trim(atmos_mesh_filepath), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
            call ESMF_LogWrite(subname//"Mesh for atmosphere is created!", ESMF_LOGMSG_INFO)
            print *, "!Mesh for atmosphere is created!"

        else
            !atmos_grid= ESMF_GridCreateNoPeriDimUfrmR( maxIndex=(/180,360 /), &
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
        ! Atmosphere to Coupler (land) Fields --  a2l
        ! I- Create empty field bundle -- a2c_fb
        ! II- Create  Fields and add them to field bundle
        ! III - Add a2c_fb to state (atm2lnd_a_state)
        !-------------------------------------------------------------------------

        a2c_fb = ESMF_FieldBundleCreate(name="a2c_fb", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"field bundle", ESMF_LOGMSG_INFO)

        ! Create individual fields and add to field bundle -- a2l

        !call fldlist_add(a2c_fldlist_num, a2c_fldlist, 'dum_var2'      )
        a2c_fldlist_num = 3

        do n = 1,a2c_fldlist_num

           print *, "**********************************************************"
           print *, "creating field for a2l:"
           print *, trim(a2c_fldlist(n)%stdname)

           ! create field
           !!! Here we want to pass pointers
           !field = ESMF_FieldCreate(atmos_mesh, ESMF_TYPEKIND_R8 ,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(a2c_fldlist(n)%stdname), rc=rc)
           field = ESMF_FieldCreate(atmos_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(a2c_fldlist(n)%stdname), farrayPtr=a2c_fldlist(n)%farrayptr1d, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
           call ESMF_FieldFill(field, dataFillScheme = "sincos" , rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

           ! add field to field bundle
           call ESMF_FieldBundleAdd(a2c_fb, (/field/), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out


           print *, a2c_fldlist(n)%farrayptr1d
           print *, "this field is created"

        enddo

        print *, "!Fields to  Coupler (atmos to  land ) (a2c_fb) Field Bundle Created!"

        ! Add field bundle to state
        call ESMF_StateAdd(atm2lnd_a_state, (/a2c_fb/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"atm2lnd_a_state is filled with dummy_var field bundle!", ESMF_LOGMSG_INFO)
        print *, "!atm2lnd_a_state is filld with dummy_var field bundle!"

        !-------------------------------------------------------------------------
        ! Coupler (land) to Atmosphere Fields --  l2a
        ! I- Create Field Bundle -- c2a_fb for now 
        ! II- Create  Fields and add them to field bundle 
        ! III - Add c2a_fb to state (lnd2atm_a_state)
        !-------------------------------------------------------------------------

        c2a_fb = ESMF_FieldBundleCreate (name="c2a_fb", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        ! Create individual fields and add to field bundle -- l2a
        c2a_fldlist_num = 3

        do n = 1,c2a_fldlist_num

           ! create field
           !!! Here we want to pass pointers
           if (mesh_switch) then
              !field = ESMF_FieldCreate(atmos_mesh, ESMF_TYPEKIND_R8 ,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(c2a_fldlist(n)%stdname), rc=rc)
              !field = ESMF_FieldCreate(atmos_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(c2a_fldlist(n)%stdname), farrayPtr=c2a_fldlist(n)%farrayptr1d, rc=rc)
              field = ESMF_FieldCreate(atmos_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(a2c_fldlist(n)%stdname), farrayPtr=a2c_fldlist(n)%farrayptr1d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
          else
              field = ESMF_FieldCreate(atmos_grid, name=trim(c2a_fldlist(n)%stdname), farrayPtr=c2a_fldlist(n)%farrayptr2d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
           end if

           ! add field to field bundle
           call ESMF_FieldBundleAdd(c2a_fb, (/field/), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

           print *, "**********************************************************"
           print *, "creating field for c2a:"
           print *, trim(c2a_fldlist(n)%stdname)
           !print *, c2a_fldlist(n)%farrayptr1d

        enddo

        print *, "!Fields For Coupler (c2a_fldlist) Field Bundle Created!"

        ! Add field bundle to state
        call ESMF_StateAdd(lnd2atm_a_state, (/c2a_fb/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        print *, "!lnd2atm_a_state is filld with dummy_var field bundle!"


    end subroutine atmos_init

    subroutine atmos_run(comp, importState, exportState, clock, rc)
        type(ESMF_GridComp)  :: comp
        type(ESMF_State)     :: importState, exportState
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        character(len=*), parameter :: subname=trim(modname)//': [atmos_run] '

        ! Initialize return code
        rc = ESMF_SUCCESS

        call ESMF_LogWrite(subname//"atmos run has not been implemented yet", ESMF_LOGMSG_INFO)
    end subroutine atmos_run

    subroutine atmos_final(comp, importState, exportState, clock, rc)
        type(ESMF_GridComp)  :: comp
        type(ESMF_State)     :: importState, exportState
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        character(len=*), parameter :: subname=trim(modname)//': [atmos_final] '
        type (ESMF_FieldBundle)     ::  import_fieldbundle, export_fieldbundle

        ! Initialize return code
        rc = ESMF_SUCCESS

        call ESMF_StateGet(importState, "c2a_fb", import_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        call ESMF_StateGet(exportState, "a2c_fb", export_fieldbundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail ou


        call ESMF_FieldBundleDestroy(import_fieldbundle, rc=rc)
        call ESMF_FieldBundleDestroy(export_fieldbundle, rc=rc)

        call ESMF_LogWrite(subname//"atmos_final has not been implemented yet", ESMF_LOGMSG_INFO)

    end subroutine atmos_final

end module atmos_cap
