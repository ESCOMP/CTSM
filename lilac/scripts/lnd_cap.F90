module lnd_cap
  use ESMF
  use lilac_utils

  implicit none

  character(*), parameter :: modname =  "  lnd_cap"

  !!integer, parameter      :: fldsMax = 100

  type(ESMF_Field), public, save :: field
  type(ESMF_Field), public, save :: field_sie, field_u

  type(fld_list_type), public, allocatable ::  c2l_fldlist(:)
  type(fld_list_type), public, allocatable ::  l2c_fldlist(:)

  !private

  public lnd_register
  !public  :: add_fields
  !public  :: import_fields
  !public  :: export_fields

  contains

!-------------------------------------------------------------------------
!    land register
!-------------------------------------------------------------------------
  subroutine lnd_register(comp, rc)

    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out) :: rc
    character(len=*), parameter :: subname=trim(modname)//': [lnd_register] '

    print *, "in lnd register routine"

    rc = ESMF_SUCCESS
    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=lnd_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=lnd_run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=lnd_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

  end subroutine lnd_register

!-------------------------------------------------------------------------
!    land init
!-------------------------------------------------------------------------

  subroutine lnd_init(comp, atm2lnd_l_state, lnd2atm_l_state, clock, rc)

    type (ESMF_GridComp)             :: comp
    type (ESMF_State)                :: atm2lnd_l_state, lnd2atm_l_state
    type (ESMF_Clock)                :: clock
    integer, intent(out)             :: rc

    type (ESMF_FieldBundle)      ::  l2c_fb , c2l_fb
    integer                      ::  n


    logical mesh_switch
    integer                          :: petCount, localrc, urc
    type(ESMF_Mesh)                  :: lnd_mesh
    character(len=ESMF_MAXSTR)       :: lnd_mesh_filepath

    character(len=*), parameter :: subname=trim(modname)//': [lnd_init] '

    type(ESMF_Grid)                  :: lnd_grid

    integer                          ::  c2l_fldlist_num
    integer                          ::  l2c_fldlist_num


    !integer                    :: regDecomp(:,:)
    ! Initialize return code
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//"------------------------!", ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(comp, petcount=petcount, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)


    print *, "  Empty land is created !!!!"
    print *, "in land routine routine"
    !-------------------------------------------------------------------------
    !    Read in  the mesh ----or----- Generate the grid
    !-------------------------------------------------------------------------
    mesh_switch = .true.
    if(mesh_switch) then
        print *, "creating mesh for land"
        ! For now this is our dummy mesh:
        !lnd_mesh_filepath    =      '/gpfs/fs1/p/cesmdata/cseg/inputdata/share/meshes/T62_040121_ESMFmesh.nc'
        lnd_mesh_filepath    =      '/gpfs/fs1/p/cesmdata/cseg/inputdata/share/meshes/T31_040122_ESMFmesh.nc'

        lnd_mesh             =        ESMF_MeshCreate(filename=trim(lnd_mesh_filepath), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Mesh for land is created!", ESMF_LOGMSG_INFO)
        print *, "!Mesh for land is created!"
    else
        lnd_grid = ESMF_GridCreateNoPeriDimUfrm( minIndex= (/1,1/), maxIndex=(/180,360 /), &
              maxCornerCoord=(/180._ESMF_KIND_R8, 360._ESMF_KIND_R8/), &
              minCornerCoord=(/0._ESMF_KIND_R8, 0._ESMF_KIND_R8/), &
              coordSys=ESMF_COORDSYS_CART,&
              regDecomp=(/petcount,1/),&
                                rc=rc)
        call ESMF_GridCompGet(comp, grid= lnd_grid , petcount=petcount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Grid for land is created!", ESMF_LOGMSG_INFO)
        print *, "Grid for land is created!"
    endif



        !-------------------------------------------------------------------------
        ! Coupler (land) to Atmosphere Fields --  l2a
        ! I- Create Field Bundle -- l2c_fb for now
        ! II- Create  Fields and add them to field bundle
        ! III - Add l2c_fb to state (lnd2atm_l_state)
        !-------------------------------------------------------------------------

        l2c_fb = ESMF_FieldBundleCreate (name="l2c_fb", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        print *, 'l2c_fb is created'
        ! Create individual fields and add to field bundle -- l2a
        l2c_fldlist_num = 3

        do n = 1,l2c_fldlist_num

           ! create field
           !!! Here we want to pass pointers
           if (mesh_switch) then
           field = ESMF_FieldCreate(lnd_mesh, ESMF_TYPEKIND_R8 ,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(l2c_fldlist(n)%stdname), rc=rc)
              !field = ESMF_FieldCreate(lnd_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(l2c_fldlist(n)%stdname), farrayPtr=l2c_fldlist(n)%farrayptr1d, rc=
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
           else
              field = ESMF_FieldCreate(lnd_grid, name=trim(l2c_fldlist(n)%stdname), farrayPtr=l2c_fldlist(n)%farrayptr2d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
           end if
           ! add field to field bundle
           call ESMF_FieldBundleAdd(l2c_fb, (/field/), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

           print *, "**********************************************************"
           print *, "creating field for l2a:"
           print *, trim(l2c_fldlist(n)%stdname)
           print *, l2c_fldlist(n)%farrayptr1d

        enddo

        print *, "!Fields For Coupler (l2c_fldlist) Field Bundle Created!"

        ! Add field bundle to state
        call ESMF_StateAdd(lnd2atm_l_state, (/l2c_fb/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        print *, "!lnd2atm_l_state is filld with dummy_var field bundle!"


        !-------------------------------------------------------------------------
        ! Atmosphere to Coupler (land) Fields --  a2l
        ! I- Create empty field bundle -- c2l_fb
        ! II- Create  Fields and add them to field bundle
        ! III - Add c2l_fb to state (atm2lnd_l_state)
        !-------------------------------------------------------------------------

        c2l_fb = ESMF_FieldBundleCreate(name="c2l_fb", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out

        ! Create individual fields and add to field bundle -- a2l

        !call fldlist_add(c2l_fldlist_num, c2l_fldlist, 'dum_var2'      )
        c2l_fldlist_num = 3

        do n = 1,c2l_fldlist_num

           ! create field
           !!! Here we want to pass pointers
           field = ESMF_FieldCreate(lnd_mesh, ESMF_TYPEKIND_R8 ,  meshloc=ESMF_MESHLOC_ELEMENT , name=trim(c2l_fldlist(n)%stdname), rc=rc)
           !field = ESMF_FieldCreate(lnd_mesh, meshloc=ESMF_MESHLOC_ELEMENT, name=trim(c2l_fldlist(n)%stdname), farrayPtr=c2l_fldlist(n)%farrayptr1d, rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
           !call ESMF_FieldGet(field, farrayPtr=fldptr, rc=rc)
           !fldptr = c2l_fldlist(n)%default_value

           ! add field to field bundle
           call ESMF_FieldBundleAdd(c2l_fb, (/field/), rc=rc)
           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out


           print *, "**********************************************************"
           print *, "creating field for a2l:"
           print *, trim(c2l_fldlist(n)%stdname)
           print *, c2l_fldlist(n)%farrayptr1d

        enddo

        print *, "!Fields to  Coupler (atmos to  land ) (c2l_fb) Field Bundle Created!"

        ! Add field bundle to state
        call ESMF_StateAdd(atm2lnd_l_state, (/c2l_fb/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        print *, "!atm2lnd_l_state is filld with dummy_var field bundle!"



  end subroutine lnd_init

!-------------------------------------------------------------------------
!    land run
!-------------------------------------------------------------------------
  subroutine lnd_run(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//': [lnd_run] '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"lnd_run has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine lnd_run

!-------------------------------------------------------------------------
!    land final
!-------------------------------------------------------------------------
  subroutine lnd_final(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//': [lnd_final] '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"lnd_final is called but has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine lnd_final
  !===============================================================================





end module lnd_cap
