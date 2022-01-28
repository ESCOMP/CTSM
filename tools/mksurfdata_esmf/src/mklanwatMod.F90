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
  use mkpioMod     , only : mkpio_get_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata, create_routehandle_r8
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
    type(file_desc_t)      :: pioid
    real(r8), allocatable  :: lake_i(:)              ! input grid: percent lake
    real(r8), allocatable  :: swmp_i(:)              ! input grid: percent wetland
    real(r8), allocatable  :: lakedepth_i(:)         ! iput grid: lake depth
    integer                :: ni,no,k                ! indices
    integer                :: ns_i,ns_o              ! local sizes
    integer                :: rcode                  ! error status
    real(r8), parameter    :: min_valid_lakedepth = 0._r8
    character(len=*), parameter :: subname = ' mklakwat '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Initialize lake_o and swmp_o to 0
    ns_o = size(lake_o)
    do no = 1,ns_o
       lake_o(no) = 0.
       swmp_o(no) = 0.
    enddo

    ! Open raw data file
    ! ASSUME for now that have only 1 input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Create a route handle between the input and output mesh
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, rc)
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Determine number of elements in input mesh
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! ----------------------------------------
    ! Create %lake
    ! ----------------------------------------

    if (.not. zero_out_lake) then

       if (root_task) then
          write (ndiag,*) 'Attempting to make %lake .....'
       end if

       ! Read in lake_i
       allocate(lake_i(ns_i), stat=rcode)
       if (rcode/=0) call shr_sys_abort()
       call mkpio_get_rawdata(pioid, 'PCT_LAKE', mesh_i, lake_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Regrid lake_i to lake_o
       call regrid_rawdata(mesh_i, mesh_o, routehandle, lake_i, lake_o, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After regrid_rawdata for PCT_LAKE in "//trim(subname))
       do no = 1,size(lake_o)
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

       ! read in swmp_i
       allocate(swmp_i(ns_i), stat=rcode)
       if (rcode/=0) call shr_sys_abort()
       call mkpio_get_rawdata(pioid, 'PCT_WETLAND', mesh_i, lake_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! regrid swmp_i to swmp_o - this also returns swmp_i to be used in the global sums below
       call regrid_rawdata(mesh_i, mesh_o, routehandle, swmp_i, swmp_o, rc)
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
    ! Create lake parameter (lakedepth)
    ! ----------------------------------------

    if (root_task) then
       write (ndiag,*) 'Attempting to make lake parameters.....'
    end if

    ! lakedepth
    allocate(lakedepth_i(ns_i), stat=rcode)
    if (rcode/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'LAKEDEPTH', mesh_i, lakedepth_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! regrid lakedepth_i to lakedepth_o - this also returns lakedepth_i to be used in the global sums below
    call regrid_rawdata(mesh_i, mesh_o, routehandle, lakedepth_i, lakedepth_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid_rawdata for lakedepth in "//trim(subname))

    ! Check validity of output data
    if (min_bad(lakedepth_o, min_valid_lakedepth, 'lakedepth')) then
       call shr_sys_abort()
    end if

#ifdef TODO
    call output_diagnostics_continuous(data_i, lakedepth_o, tgridmap, "Lake Depth", "m", ndiag, tdomain%mask, frac_dst)
#endif
    deallocate (lakedepth_i)

    ! ----------------------------------------
    ! Close the single input data file and release memory
    ! ----------------------------------------

    call pio_closefile(pioid)

    call ESMF_VMLogMemInfo("Before destroy operation for lanwat ")
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operation for lanwat ")

    if (root_task) then
       write (ndiag,*) 'Successfully made lake parameters'
       write (ndiag,*)
    end if

  end subroutine mklakwat

end module mklanwatMod
