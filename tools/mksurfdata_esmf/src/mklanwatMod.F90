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
  use mkpioMod     , only : mkpio_get_rawdata 
  use mkesmfMod    , only : regrid_rawdata
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag

  implicit none
  private

  public :: mklakwat    ! make % lake,  % wetland and lake parameters

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mklakwat(file_mesh_i, file_data_i, mesh_o, &
       zero_out_lake, zero_out_wetland, lake_o, swmp_o, lakedepth_o, rc)

#ifdef TODO
    use mkdiagnosticsMod , only : output_diagnostics_continuous
#endif
    use mkchecksMod      , only : min_bad

    ! -------------------
    ! make %lake, %wetland and lake parameters
    ! -------------------

    ! input/output variables
    type(ESMF_Mesh)   , intent(in)  :: mesh_o
    character(len=*)  , intent(in)  :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i      ! input data file name
    logical           , intent(in)  :: zero_out_lake    ! if should zero glacier out
    logical           , intent(in)  :: zero_out_wetland ! if should zero glacier out
    real(r8)          , intent(out) :: lake_o(:)        ! output grid: %lake
    real(r8)          , intent(out) :: swmp_o(:)        ! output grid: %lake
    real(r8)          , intent(out) :: lakedepth_o(:)   ! output grid: lake depth (m)
    integer           , intent(out) :: rc

    ! local variables
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    real(r8), allocatable  :: lake_i(:)              ! input grid: percent lake
    real(r8), allocatable  :: swmp_i(:)              ! input grid: percent lake
    real(r8), allocatable  :: lakedepth_i(:) 
    integer                :: ni,no,ns_i,ns_o,k      ! indices
    integer                :: rcode                  ! error status
    integer                :: srcMaskValue = 0
    integer                :: dstMaskValue = -987987 ! spval for RH mask values
    integer                :: srcTermProcessing_Value = 0
    real(r8), parameter    :: min_valid_lakedepth = 0._r8
    character(len=32)      :: subname = 'mklakwat'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Initialize lake_o and swmp_o to 0
    ns_o = size(lake_o)
    do no = 1,ns_o
       lake_o(no) = 0.
       swmp_o(no) = 0.
    enddo

    ! create field on model mesh
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------------------------
    ! Create route handle for rawdata mesh to model mesh
    ! ----------------------------------------

    ! create field on input mesh (first read in input mesh)
    call ESMF_VMLogMemInfo("Before create mesh_i in lanwat")
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in lanwat")

    ! create route handle to map field_model to field_data
    call ESMF_VMLogMemInfo("Before regridstore in regrid_data")
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         srcMaskValues=(/srcMaskValue/), dstMaskValues=(/dstMaskValue/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_DSTAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in regrid_data")

    ! ----------------------------------------
    ! Create %lake
    ! ----------------------------------------

    if (.not. zero_out_lake) then

       if (root_task) then
          write (ndiag,*) 'Attempting to make %lake .....'
       end if

       ! read in rawdata
       allocate(lake_i(ns_i), stat=rcode)
       if (rcode/=0) call shr_sys_abort()
       call mkpio_get_rawdata(trim(file_data_i), 'PCT_LAKE', mesh_i, lake_i, rc)  
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! regrid lake_i to lake_o - this also returns lake_i to be used in the global sums below
       call ESMF_VMLogMemInfo("Before regrid_data in lanwat")
       call regrid_rawdata(field_i, field_o, routehandle, lake_i, lake_o, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After regrid_data in lanwat")
       do no = 1,size(lake_o)
          write(6,*)'DEBUG: n,lake_o = ',no,lake_o(no)
          if (lake_o(no) < 1.) lake_o(no) = 0.
       enddo

       deallocate (lake_i)

       if (root_task) then
          write (ndiag,*) 'Successfully made %lake'
          write (ndiag,*)
       end if
    end if

    ! ----------------------------------------
    ! Create %wetland
    ! ----------------------------------------

    if (.not. zero_out_wetland) then
       if (root_task) then
          write (ndiag,*) 'Attempting to make %wetland .....'
       end if

       ! read in rawdata
       allocate(swmp_i(ns_i), stat=rcode)
       if (rcode/=0) call shr_sys_abort()
       call mkpio_get_rawdata(trim(file_data_i), 'PCT_WETLAND', mesh_i, swmp_i, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! regrid swmp_i to swmp_o - this also returns swmp_i to be used in the global sums below
       call ESMF_VMLogMemInfo("Before regrid_data for wetland")
       call regrid_rawdata(field_i, field_o, routehandle, swmp_i, swmp_o, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After regrid_data for wetland")
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do no = 1,ns_o
          if (swmp_o(no) < 1.) swmp_o(no) = 0.
       enddo
       deallocate (swmp_i)

       if (root_task) then
          write (ndiag,*) 'Successfully made %wetland'
          write (ndiag,*)
       end if
    end if

    ! ----------------------------------------
    ! Create lake parameter (lakdepth)
    ! ----------------------------------------

    if (root_task) then
       write (ndiag,*) 'Attempting to make lake parameters.....'
    end if

    ! read in rawdata
    allocate(lakedepth_i(ns_i), stat=rcode)
    if (rcode/=0) call shr_sys_abort()
    call mkpio_get_rawdata(trim(file_data_i), 'LAKEDEPTH', mesh_i, lakedepth_i, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! regrid lakedepth_i to lakedepth_o - this also returns lakedepth_i to be used in the global sums below
    call ESMF_VMLogMemInfo("Before regrid_data for lakedepth")
    call regrid_rawdata(field_i, field_o, routehandle, lakedepth_i, lakedepth_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid_data for lakedepth")
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Check validity of output data
    if (min_bad(lakedepth_o, min_valid_lakedepth, 'lakedepth')) then
       call shr_sys_abort()
    end if

#ifdef TODO
    call output_diagnostics_continuous(data_i, lakedepth_o, tgridmap, "Lake Depth", "m", ndiag, tdomain%mask, frac_dst)
#endif

    deallocate (lakedepth_i)

    if (root_task) then
       write (ndiag,*) 'Successfully made lake parameters'
       write (ndiag,*)
    end if

    ! Release memory for route handle

    call ESMF_VMLogMemInfo("Before destroy operation for lanwat ")
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operation for lanwat ")

  end subroutine mklakwat

end module mklanwatMod
