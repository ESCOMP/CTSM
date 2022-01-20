module mklanwatMod

  !-----------------------------------------------------------------------
  ! make %lake and %wetland from input lake / wetland data
  ! also make lake parameters
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkvarpar     , only : re	
  use mkesmfMod    , only : regrid_data
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag

  implicit none
  private

  public mklakwat    ! make % lake
  public mkwetlnd    ! make % wetland
  public mklakparams ! make lake parameters

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mklakwat(file_mesh_i, file_data_i, mesh_o, zero_out, lake_o, rc)

    ! -------------------
    ! make %lake
    ! -------------------

    ! input/output variables
    type(ESMF_Mesh)   , intent(in)  :: mesh_o
    character(len=*)  , intent(in)  :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i ! input data file name
    logical           , intent(in)  :: zero_out    ! if should zero glacier out
    real(r8)          , intent(out) :: lake_o(:)   ! output grid: %lake
    integer           , intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)       :: mesh_i
    type(ESMF_Field)      :: field_i
    type(ESMF_Field)      :: field_o
    real(r8), allocatable :: lake_i(:)              ! input grid: percent lake
    real(r8)              :: sum_fldi               ! global sum of dummy input fld
    real(r8)              :: sum_fldo               ! global sum of dummy output fld
    real(r8)              :: glake_i                ! input  grid: global lake
    real(r8)              :: garea_i                ! input  grid: global area
    real(r8)              :: glake_o                ! output grid: global lake
    real(r8)              :: garea_o                ! output grid: global area
    integer               :: ni,no,ns_i,ns_o,k      ! indices
    integer               :: rcode                  ! error status
    character(len=32)     :: subname = 'mklakwat'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,*) 'Attempting to make %lake .....'
    end if

    ! Initialize lake_o to 0
    ns_o = size(lake_o)
    do no = 1,ns_o
       lake_o(no) = 0.
    enddo

    if ( .not. zero_out ) then

       ! create field on input mesh (first read in input mesh)
       mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(lake_i(ns_i), stat=rcode)
       if (rcode/=0) call shr_sys_abort()

       ! create field on model mesh
       field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! regrid lake_i to lake_o - this also returns lake_i to be used in the global sums below
       call regrid_data(field_i, field_o, 'PCT_LAKE', trim(file_data_i), lake_i, lake_o, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do no = 1,size(lake_o)
          if (lake_o(no) < 1.) lake_o(no) = 0.
       enddo

       ! -----------------------------------------------------------------
       ! Error check prep
       ! Global sum of output field -- must multiply by fraction of
       ! output grid that is land as determined by input grid
       ! -----------------------------------------------------------------
       ! TODO: not sure if we still want this with the ESMF regridding

       ! -----------------------------------------------------------------
       ! Error check2
       ! Compare global areas on input and output grids
       ! -----------------------------------------------------------------
       ! TODO: implement this
       ! ! Input grid

       deallocate (lake_i)
       ! TODO: destroy route handle and created field
    end if

    if (root_task) then
       write (ndiag,*) 'Successfully made %lake'
       write (ndiag,*)
    end if

  end subroutine mklakwat

  !===============================================================
  subroutine mkwetlnd(file_mesh_i, file_data_i, mesh_o, zero_out, swmp_o, rc)

    ! -------------------
    ! make %wetland
    ! -------------------

    ! input/output variables
    type(ESMF_Mesh)   , intent(in)  :: mesh_o
    character(len=*)  , intent(in)  :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i ! input data file name
    logical           , intent(in)  :: zero_out    ! if should zero glacier out
    real(r8)          , intent(out) :: swmp_o(:)   ! output grid: %lake
    integer           , intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)       :: mesh_i
    type(ESMF_Field)      :: field_i
    type(ESMF_Field)      :: field_o
    real(r8), allocatable :: swmp_i(:)              ! input grid: percent lake
    real(r8)              :: sum_fldi               ! global sum of dummy input fld
    real(r8)              :: sum_fldo               ! global sum of dummy output fld
    real(r8)              :: glake_i                ! input  grid: global lake
    real(r8)              :: garea_i                ! input  grid: global area
    real(r8)              :: glake_o                ! output grid: global lake
    real(r8)              :: garea_o                ! output grid: global area
    integer               :: ni,no,ns_i,ns_o,k      ! indices
    integer               :: rcode                  ! error status
    character(len=32)     :: subname = 'mkwtlnd'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,*) 'Attempting to make %wetland .....'
    end if

    ! Initialize swmp_o to 0
    ns_o = size(swmp_o)
    do no = 1,ns_o
       swmp_o(no) = 0.
    enddo

    if ( .not. zero_out ) then

       ! create field on input mesh (first read in input mesh)
       mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(swmp_i(ns_i), stat=rcode)
       if (rcode/=0) call shr_sys_abort()

       ! create field on model mesh
       field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! regrid swmp_i to swmp_o - this also returns swmp_i to be used in the global sums below
       call regrid_data(field_i, field_o, 'PCT_WETLAND', trim(file_data_i), swmp_i, swmp_o, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do no = 1,ns_o
          if (swmp_o(no) < 1.) swmp_o(no) = 0.
       enddo

       ! -----------------------------------------------------------------
       ! Error check prep
       ! Global sum of output field -- must multiply by fraction of
       ! output grid that is land as determined by input grid
       ! -----------------------------------------------------------------
       ! TODO: not sure if we still want this with the ESMF regridding

       ! -----------------------------------------------------------------
       ! Error check2
       ! Compare global areas on input and output grids
       ! -----------------------------------------------------------------
       ! TODO: implememt this

       deallocate (swmp_i)
       ! TODO: destroy route handle and created field

    end if

    if (root_task) then
       write (ndiag,*) 'Successfully made %wetland'
       write (ndiag,*)
    end if

  end subroutine mkwetlnd

  !===============================================================

  subroutine mklakparams(file_mesh_i, file_data_i, mesh_o, lakedepth_o, rc)

    ! -------------------
    ! make lake parameters (currently just lake depth)
    ! -------------------

    !use mkdiagnosticsMod , only : output_diagnostics_continuous
    use mkchecksMod      , only : min_bad

    ! input/output variables
    character(len=*)  , intent(in)  :: file_mesh_i    ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i    ! input data file name
    type(ESMF_Mesh)   , intent(in)  :: mesh_o
    real(r8)          , intent(out) :: lakedepth_o(:) ! output grid: lake depth (m)
    integer           , intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)       :: mesh_i
    type(ESMF_Field)      :: field_i
    type(ESMF_Field)      :: field_o
    real(r8), allocatable :: lakedepth_i(:)    ! input raw data
    real(r8)              :: sum_fldi          ! global sum of dummy input fld
    real(r8)              :: sum_fldo          ! global sum of dummy output fld
    real(r8)              :: glake_i           ! input  grid: global lake
    real(r8)              :: garea_i           ! input  grid: global area
    real(r8)              :: glake_o           ! output grid: global lake
    real(r8)              :: garea_o           ! output grid: global area
    integer               :: ni,no,ns_i,ns_o,k ! indices
    integer               :: rcode             ! error status
    real(r8), parameter   :: min_valid_lakedepth = 0._r8
    character(len=32)     :: subname = 'mklakparams'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,*) 'Attempting to make lake parameters.....'
    end if

    ! create field on input mesh (first read in input mesh)
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lakedepth_i(ns_i), stat=rcode)
    if (rcode/=0) call shr_sys_abort()

    ! create field on model mesh
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! regrid lakedepth_i to lakedepth_o - this also returns lakedepth_i to be used in the global sums below
    call regrid_data(field_i, field_o, 'LAKEDEPTH', trim(file_data_i), lakedepth_i, lakedepth_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Check validity of output data
    if (min_bad(lakedepth_o, min_valid_lakedepth, 'lakedepth')) then
       call shr_sys_abort()
    end if

    ! TODO: implement the following
    !call output_diagnostics_continuous(data_i, lakedepth_o, tgridmap, "Lake Depth", "m", ndiag, tdomain%mask, frac_dst)

    deallocate (lakedepth_i)

    if (root_task) then
       write (ndiag,*) 'Successfully made lake parameters'
       write (ndiag,*)
    end if

  end subroutine mklakparams

end module mklanwatMod
