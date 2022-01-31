module mkesmfMod

  use ESMF
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkUtilsMod   , only : chkerr

  implicit none
  private

  public :: regrid_rawdata
  public :: get_meshareas
  public :: create_routehandle_frac
  public :: create_routehandle_r4
  public :: create_routehandle_r8

  interface create_routehandle_frac
     module procedure create_routehandle_frac_r4
     module procedure create_routehandle_frac_r8
  end interface create_routehandle_frac

  interface get_meshareas
     module procedure get_meshareas_r4
     module procedure get_meshareas_r8
  end interface get_meshareas

  interface regrid_rawdata
     module procedure regrid_rawdata1d_r4
     module procedure regrid_rawdata1d_r8
     module procedure regrid_rawdata2d_r4
     module procedure regrid_rawdata2d_r8
  end interface regrid_rawdata

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine create_routehandle_frac_r8(mesh_i, mesh_o, field_type, routehandle, frac_o, rc)

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    type(ESMF_Mesh)        , intent(in)    :: mesh_o
    character(len=*)       , intent(in)    :: field_type
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    real(r4)               , intent(inout) :: frac_o(:)
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: srcMaskValue = 0          ! ignore source points where the mesh mask is 0
    integer           :: dstMaskValue = -987987    ! don't ingore any destination points
    integer           :: srcTermProcessing_Value = 0
    type(ESMF_Field)  :: field_i
    type(ESMF_Field)  :: field_o
    type(ESMF_Field)  :: dstFracField
    real(r8), pointer :: dataptr(:)
    character(len=*), parameter :: subname = 'create_routehandle'
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dstFracField = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create route handle to map field_model to field_data
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_FRACAREA, &
         srcMaskValues=(/srcMaskValue/), &
         dstMaskValues=(/dstMaskValue/), &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         dstFracField= dstFracField, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    call ESMF_FieldGet(dstFracField, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    frac_o(:) = dataptr(:)

    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(dstFracField, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  end subroutine create_routehandle_frac_r8

  !===============================================================
  subroutine create_routehandle_frac_r4(mesh_i, mesh_o, routehandle, frac_o, rc)

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    type(ESMF_Mesh)        , intent(in)    :: mesh_o
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    real(r4)               , intent(inout) :: frac_o(:)
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: srcMaskValue = 0          ! ignore source points where the mesh mask is 0
    integer           :: dstMaskValue = -987987    ! don't ingore any destination points
    integer           :: srcTermProcessing_Value = 0
    type(ESMF_Field)  :: field_i
    type(ESMF_Field)  :: field_o
    type(ESMF_Field)  :: dstFracField
    real(r4), pointer :: dataptr(:)
    character(len=*), parameter :: subname = 'create_routehandle'
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dstFracField = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create route handle to map field_model to field_data
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_FRACAREA, &
         srcMaskValues=(/srcMaskValue/), &
         dstMaskValues=(/dstMaskValue/), &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         dstFracField= dstFracField, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    call ESMF_FieldGet(dstFracField, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    frac_o(:) = dataptr(:)

    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(dstFracField, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  end subroutine create_routehandle_frac_r4

  !===============================================================
  subroutine create_routehandle_r4(mesh_i, mesh_o, routehandle, rc)

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    type(ESMF_Mesh)        , intent(in)    :: mesh_o
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    integer                , intent(out)   :: rc

    ! local variables
    integer          :: srcMaskValue = 0          ! ignore source points where the mesh mask is 0
    integer          :: dstMaskValue = -987987    ! don't ingore any destination points
    integer          :: srcTermProcessing_Value = 0
    type(ESMF_Field) :: field_i
    type(ESMF_Field) :: field_o
    character(len=*), parameter :: subname = 'create_routehandle_r4'
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create route handle to map field_model to field_data
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_FRACAREA, &
         srcMaskValues=(/srcMaskValue/), &
         dstMaskValues=(/dstMaskValue/), &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  end subroutine create_routehandle_r4

  !===============================================================
  subroutine create_routehandle_r8(mesh_i, mesh_o, routehandle, rc)

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    type(ESMF_Mesh)        , intent(in)    :: mesh_o
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    integer                , intent(out)   :: rc

    ! local variables
    integer          :: srcMaskValue = 0          ! ignore source points where the mesh mask is 0
    integer          :: dstMaskValue = -987987    ! don't ingore any destination points
    integer          :: srcTermProcessing_Value = 0
    type(ESMF_Field) :: field_i
    type(ESMF_Field) :: field_o
    character(len=*), parameter :: subname = 'create_routehandle_r8'
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create route handle to map field_model to field_data
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_FRACAREA, &
         srcMaskValues=(/srcMaskValue/), &
         dstMaskValues=(/dstMaskValue/), &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  end subroutine create_routehandle_r8

  !===============================================================
  subroutine regrid_rawdata1d_r4(mesh_i, mesh_o, routehandle, data_i, data_o, rc)

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    type(ESMF_Mesh)        , intent(in)    :: mesh_o
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    real(r4)               , intent(in)    :: data_i(:)
    real(r4)               , intent(inout) :: data_o(:)
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)  :: field_i
    type(ESMF_Field)  :: field_o
    real(r4), pointer :: dataptr(:)
    logical           :: checkflag = .false.
    character(len=*), parameter :: subname = 'regrid_rawdata1d_r4'
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Interpolate data_i to data_o
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = data_i(:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r8
    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    data_o(:) = dataptr(:)

    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  end subroutine regrid_rawdata1d_r4

  !===============================================================
  subroutine regrid_rawdata1d_r8(mesh_i, mesh_o,  routehandle, data_i, data_o, rc)

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    type(ESMF_Mesh)        , intent(in)    :: mesh_o
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    real(r8)               , intent(in)    :: data_i(:)
    real(r8)               , intent(inout) :: data_o(:)
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)  :: field_i
    type(ESMF_Field)  :: field_o
    real(r8), pointer :: dataptr(:)
    logical           :: checkflag = .false.
    character(len=*), parameter :: subname = 'regrid_rawdata1d_r8'
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Interpolate data_i to data_o
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = data_i(:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r8
    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    data_o(:) = dataptr(:)

    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  end subroutine regrid_rawdata1d_r8

  !===============================================================
  subroutine regrid_rawdata2d_r4(mesh_i, mesh_o,  routehandle, data_i, data_o, lbound, ubound, rc)

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    type(ESMF_Mesh)        , intent(in)    :: mesh_o
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    integer                , intent(in)    :: lbound
    integer                , intent(in)    :: ubound
    real(r4)               , intent(in)    :: data_i(:,:)
    real(r4)               , intent(inout) :: data_o(:,:)
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)  :: field_i
    type(ESMF_Field)  :: field_o
    logical           :: checkflag = .false.
    real(r4), pointer :: dataptr(:,:)
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/lbound/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/lbound/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Interpolate data_i to data_o
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:,:) = data_i(:,:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:,:) = 0._r8
    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    data_o(:,:) = dataptr(:,:)

    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  end subroutine regrid_rawdata2d_r4

  !===============================================================
  subroutine regrid_rawdata2d_r8(mesh_i, mesh_o,  routehandle, data_i, data_o, lbound, ubound, rc)

    ! input/output variables
    type(ESMF_Mesh)        , intent(in)    :: mesh_i
    type(ESMF_Mesh)        , intent(in)    :: mesh_o
    type(ESMF_RouteHandle) , intent(inout) :: routehandle
    integer                , intent(in)    :: lbound
    integer                , intent(in)    :: ubound
    real(r8)               , intent(in)    :: data_i(:,:)
    real(r8)               , intent(inout) :: data_o(:,:)
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)  :: field_i
    type(ESMF_Field)  :: field_o
    logical           :: checkflag = .false.
    real(r8), pointer :: dataptr(:,:)
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/lbound/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
         ungriddedLbound=(/lbound/), ungriddedUbound=(/ubound/), gridToFieldMap=(/2/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Interpolate data_i to data_o
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:,:) = data_i(:,:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:,:) = 0._r8

    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    data_o(:,:) = dataptr(:,:)

    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

  end subroutine regrid_rawdata2d_r8

  !===============================================================
  subroutine get_meshareas_r8(mesh, areas, rc)

    ! input/output variables
    type(ESMF_Mesh) , intent(in)    :: mesh
    real(r8)        , intent(inout) :: areas(:)
    integer         , intent(out)   :: rc

    ! local variables
    real(r8), pointer :: dataptr(:)
    type(ESMF_Field)  :: lfield
    integer           :: ns
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_MeshGet(mesh, numOwnedElements=ns, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    areas(:) = dataptr(:)

    call ESMF_FieldDestroy(lfield, nogarbage = .true., rc=rc)
  end subroutine get_meshareas_r8

  !===============================================================
  subroutine get_meshareas_r4(mesh, areas, rc)

    ! input/output variables
    type(ESMF_Mesh) , intent(in)    :: mesh
    real(r4)        , intent(inout) :: areas(:)
    integer         , intent(out)   :: rc

    ! local variables
    real(r8), pointer :: dataptr(:)
    type(ESMF_Field)  :: lfield
    integer           :: ns 
    ! --------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_MeshGet(mesh, numOwnedElements=ns, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    areas(:) = real(dataptr(:), kind=r4)

    call ESMF_FieldDestroy(lfield, nogarbage = .true., rc=rc)
  end subroutine get_meshareas_r4

end module mkesmfMod
