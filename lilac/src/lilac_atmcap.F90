module lilac_atmcap

  !-----------------------------------------------------------------------
  ! This is an ESMF lilac cap for the host atmosphere
  !-----------------------------------------------------------------------

  use ESMF
  use shr_kind_mod  , only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_sys_mod   , only : shr_sys_abort
  use lilac_methods , only : chkerr
  use lilac_constants, only : logunit
  use ctsm_LilacCouplingFields, only : a2l_fields, l2a_fields

  implicit none

  public :: lilac_atmcap_init
  public :: lilac_atmcap_register

  ! Time invariant input from host atmosphere
  integer , public, allocatable :: gindex_atm(:) ! global index space
  real(r8), private, allocatable :: atm_lons(:)   ! local longitudes
  real(r8), private, allocatable :: atm_lats(:)   ! local latitudes
  integer , public              :: atm_global_nx
  integer , public              :: atm_global_ny

  ! Time variant input from host atmosphere
  real(r8) :: nextsw_cday = 1.e36_r8  ! calendar day of the next sw calculation

  integer                 :: mytask
  integer     , parameter :: debug = 0 ! internal debug level
  character(*), parameter :: u_FILE_u = &
       __FILE__

!========================================================================
contains
!========================================================================

  subroutine lilac_atmcap_init_vars(atm_gindex_in, atm_lons_in, atm_lats_in, atm_global_nx_in, atm_global_ny_in)

    ! input/output variables
    integer , intent(in) :: atm_gindex_in(:)
    real(r8), intent(in) :: atm_lons_in(:)
    real(r8), intent(in) :: atm_lats_in(:)
    integer , intent(in) :: atm_global_nx_in
    integer , intent(in) :: atm_global_ny_in

    ! glocal variables
    integer :: n, lsize, fileunit
    ! --------------------------------------------

    lsize = size(atm_gindex_in)
    allocate(gindex_atm(lsize))
    allocate(atm_lons(lsize))
    allocate(atm_lats(lsize))

    ! set module variables
    gindex_atm(:) = atm_gindex_in(:)
    atm_lons(:)   = atm_lons_in(:)
    atm_lats(:)   = atm_lats_in(:)
    atm_global_nx = atm_global_nx_in
    atm_global_ny = atm_global_ny_in

  end subroutine lilac_atmcap_init_vars

  !========================================================================
  subroutine lilac_atmcap_register (comp, rc)

    ! input/output variables
    type(ESMF_GridComp)          :: comp   ! must not be optional
    integer, intent(out)         :: rc

    ! local variables
    type(ESMF_VM)                :: vm
    character(len=*), parameter  :: subname='(lilac_atmcap_register): '
    !-------------------------------------------------------------------------

    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
       print *, "in user register routine"
    end if

    ! Initialize return code
    rc = ESMF_SUCCESS

    ! Set the entry points for standard ESMF Component methods
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=lilac_atmcap_init, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=lilac_atmcap_run, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=lilac_atmcap_final, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lilac_atmcap_register

!========================================================================

  subroutine lilac_atmcap_init (comp, lnd2atm_state, atm2lnd_state, clock, rc)

    ! input/output variables
    type (ESMF_GridComp) :: comp
    type (ESMF_State)    :: lnd2atm_state
    type (ESMF_State)    :: atm2lnd_state
    type (ESMF_Clock)    :: clock
    integer, intent(out) :: rc

    ! local variables
    integer                :: fileunit
    type(ESMF_Mesh)        :: atm_mesh
    type(ESMF_DistGrid)    :: atm_distgrid
    type(ESMF_Field)       :: field
    type(ESMF_FieldBundle) :: c2a_fb , a2c_fb
    integer                :: n, i, ierr
    integer                :: lsize
    character(len=cl)      :: atm_mesh_filename
    character(len=cl)      :: cvalue
    integer                :: spatialDim
    integer                :: numOwnedElements
    real(r8), pointer      :: ownedElemCoords(:)
    real(r8)               :: mesh_lon, mesh_lat
    real(r8), parameter    :: tolerance = 1.e-4_r8
    character(len=*), parameter :: subname='(lilac_atmcap_init): '
    !-------------------------------------------------------------------------

    namelist /lilac_atmcap_input/ atm_mesh_filename

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//"------------------------!", ESMF_LOGMSG_INFO)

    !-------------------------------------------------------------------------
    ! Read in the atm mesh
    !-------------------------------------------------------------------------

    ! read in mesh file name from namelist
    open(newunit=fileunit, status="old", file="lilac_in")
    read(fileunit, lilac_atmcap_input, iostat=ierr)
    if (ierr > 0) then
       call shr_sys_abort(trim(subname) // 'error reading in lilac_atmcap_input')
    end if
    close(fileunit)

    ! Note that in the call to lilac_atm the host atmospere sent both the gindex_atm and
    ! the atm_mesh_filename that were then set as module variables here

    atm_distgrid = ESMF_DistGridCreate (arbSeqIndexList=gindex_atm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: the addUserArea failed for the 4x5 grid - need to have a
    ! more robust approach - unless the area will simply be ignored for now?
    ! atm_mesh = ESMF_MeshCreate(filename=trim(atm_mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
    !      elementDistGrid=atm_distgrid, addUserArea=.true., rc=rc)
    atm_mesh = ESMF_MeshCreate(filename=trim(atm_mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistGrid=atm_distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) then
       call shr_sys_abort(trim(subname) // 'Error: failure in creating lilac atmcap from meshfile '//&
            trim(atm_mesh_filename))
    end if

    call ESMF_LogWrite(trim(subname)//"Mesh for atmosphere is created for "//trim(atm_mesh_filename), ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, trim(subname) // "Mesh for atmosphere is created for "//trim(atm_mesh_filename)
    end if

    !-------------------------------------------------------------------------
    ! Check that lons and lats from the host atmospere match those read
    ! in from the atm mesh file
    !-------------------------------------------------------------------------

    call ESMF_MeshGet(atm_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(gindex_atm)
    if (numOwnedElements /= lsize) then
       print *, 'numOwnedElements in atm_mesh = ',numOwnedElements
       print *, 'local size from gindex_atm from host atm = ',lsize
       call shr_sys_abort('ERROR: numOwnedElements is not equal to lsize')
    end if

    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    call ESMF_MeshGet(atm_mesh, ownedElemCoords=ownedElemCoords, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1, lsize
       mesh_lon = ownedElemCoords(2*n-1)
       mesh_lat = ownedElemCoords(2*n)
       if ( abs(mesh_lon - atm_lons(n)) > tolerance) then
          write(logunit,101),n, atm_lons(n), mesh_lon
101       format('ERROR: lilac_atmcap: n, lon, mesh_lon = ',i6,2(f20.10,2x))
          call shr_sys_abort()
       end if
       if ( abs(mesh_lat - atm_lats(n)) > tolerance) then
          write(logunit,102),n, atm_lats(n), mesh_lat
102       format('ERROR: lilac_atmcap: n, lat, mesh_lat = ',i6,2(f20.10,2x))
          call shr_sys_abort()
       end if
    end do
    deallocate(ownedElemCoords)

    !-------------------------------------------------------------------------
    ! Create a2c_fb field bundle and add to atm2lnd_state
    !-------------------------------------------------------------------------

    ! create empty field bundle "a2c_fb"
    a2c_fb = ESMF_FieldBundleCreate(name="a2c_fb", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create fields and add to field bundle
    do n = 1, a2l_fields%num_fields()
       field = ESMF_FieldCreate(atm_mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
            name=a2l_fields%get_fieldname(n), farrayPtr=a2l_fields%get_dataptr(n), &
            rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldBundleAdd(a2c_fb, (/field/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (debug > 0) then
          call ESMF_FieldPrint(field,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    ! add field bundle to atm2lnd_state
    call ESMF_StateAdd(atm2lnd_state, (/a2c_fb/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"lilac a2c_fb fieldbundle created and added to atm2lnd_state", ESMF_LOGMSG_INFO)
    if (mytask == 0) then
       print *, "lilac a2c_fb fieldbundle created and added to atm2lnd_state"
    end if

    ! add nextsw_cday attributes
    write(cvalue,*) nextsw_cday
    call ESMF_AttributeSet(atm2lnd_state, name="nextsw_cday", value=cvalue, rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------------------------------------
    ! Create c2a_fb field bundle and add to lnd2atm_state
    !-------------------------------------------------------------------------

    ! create empty field bundle "c2a_fb"
    c2a_fb = ESMF_FieldBundleCreate (name="c2a_fb", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create fields and add to field bundle
    do n = 1, l2a_fields%num_fields()
       field = ESMF_FieldCreate(atm_mesh, meshloc=ESMF_MESHLOC_ELEMENT, &
            name=l2a_fields%get_fieldname(n), farrayPtr=l2a_fields%get_dataptr(n), &
            rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldBundleAdd(c2a_fb, (/field/), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (debug > 0) then
          call ESMF_FieldPrint(field,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    ! add field bundle to lnd2atm_state
    call ESMF_StateAdd(lnd2atm_state, (/c2a_fb/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//"lilac c2a_fb fieldbundle is done and added to lnd2atm_state", ESMF_LOGMSG_INFO)

  end subroutine lilac_atmcap_init

!========================================================================
  subroutine lilac_atmcap_run(comp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Initialize return code
    ! This routine does nothing - its all in the atm->lnd coupler
    rc = ESMF_SUCCESS

  end subroutine lilac_atmcap_run

!========================================================================
  subroutine lilac_atmcap_final(comp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: comp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type (ESMF_FieldBundle) ::  import_fieldbundle, export_fieldbundle
    character(len=*), parameter :: subname='( lilac_atmcap_final): '
    !-------------------------------------------------------------------------

    ! Initialize return code
    rc = ESMF_SUCCESS

    call ESMF_StateGet(importState, "c2a_fb", import_fieldbundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(exportState, "a2c_fb", export_fieldbundle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleDestroy(import_fieldbundle, rc=rc)
    call ESMF_FieldBundleDestroy(export_fieldbundle, rc=rc)

    call ESMF_LogWrite(subname//"Finished lilac_atmcap_final", ESMF_LOGMSG_INFO)

  end subroutine lilac_atmcap_final

end module lilac_atmcap
