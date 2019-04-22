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

!-------------------------------------------------------------------------
!    land register
!-------------------------------------------------------------------------
  subroutine lnd_register(comp, rc)

    type(ESMF_GridComp)   :: comp   ! must not be optional
    integer, intent(out) :: rc
    character(len=*), parameter :: subname=trim(modname)//':(lnd_register) '

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

  subroutine lnd_init(comp, importState, exportState, clock, rc)
    type (ESMF_GridComp)             :: comp
    type (ESMF_State)                :: importState, exportState
    type (ESMF_Clock)                :: clock
    integer, intent(out)             :: rc

    logical mesh_switch
    integer                          :: petCount, localrc, urc
    type(ESMF_Mesh)                  :: lnd_mesh
    character(len=ESMF_MAXSTR)       :: lnd_mesh_filepath

    character(len=*), parameter :: subname=trim(modname)//':(lnd_register) '

    type(ESMF_Grid)                  :: lnd_grid

    print *, "  Empty land is created !!!!"
    rc = ESMF_SUCCESS


    print *, "in land routine routine"
    !-------------------------------------------------------------------------
    !    Read in  the mesh ----or----- Generate the grid
    !-------------------------------------------------------------------------
    mesh_switch = .false.
    if(mesh_switch) then
        ! For now this is our dummy mesh:
        lnd_mesh_filepath='/gpfs/fs1/p/cesmdata/cseg/inputdata/share/meshes/T31_040122_ESMFmesh.nc'

        lnd_mesh = ESMF_MeshCreate(filename=trim(lnd_mesh_filepath), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Mesh for land is created!", ESMF_LOGMSG_INFO)
        print *, "!Mesh for land is created!"
    else
    call ESMF_GridCompGet(comp, petcount=petcount, rc=rc)
        lnd_grid = ESMF_GridCreateNoPeriDimUfrm( minIndex= (/1,1/), maxIndex=(/180,360 /), &
              maxCornerCoord=(/180._ESMF_KIND_R8, 360._ESMF_KIND_R8/), &
              minCornerCoord=(/0._ESMF_KIND_R8, 0._ESMF_KIND_R8/), &
              coordSys=ESMF_COORDSYS_CART,&
              regDecomp=(/petcount,1/),&
                                rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return  ! bail out
        call ESMF_LogWrite(subname//"Grid for land is created!", ESMF_LOGMSG_INFO)
        print *, "Grid for land is created!"
    endif
  end subroutine lnd_init

!-------------------------------------------------------------------------
!    land run
!-------------------------------------------------------------------------
  subroutine lnd_run(comp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=*), parameter :: subname=trim(modname)//':(lnd_run) '

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

    character(len=*), parameter :: subname=trim(modname)//':(lnd_final) '

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"lnd_final is called but has not been implemented yet", ESMF_LOGMSG_INFO)

  end subroutine lnd_final
  !===============================================================================





end module lnd_cap
