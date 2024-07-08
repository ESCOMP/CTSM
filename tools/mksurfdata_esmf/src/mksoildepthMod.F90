module mksoildepthMod

  !-----------------------------------------------------------------------
  ! make fraction soildepth from input soildepth data
  !-----------------------------------------------------------------------
  !
  use ESMF
  use pio            , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite, pio_syncfile
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4, cs => shr_kind_cs
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8
  use mkdiagnosticsMod , only : output_diagnostics_area
  use mkutilsMod     , only : chkerr
  use mkchecksMod    , only : min_bad, max_bad
  use mkfileMod      , only : mkfile_output
  use mkvarctl

  implicit none
  private

  public mksoildepth           ! regrid soildepth data

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mksoildepth(file_mesh_i, file_data_i, mesh_o, pioid_o, rc)
    !
    ! make soildepth
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle)      :: routehandle
    type(ESMF_Mesh)             :: mesh_i
    type(file_desc_t)           :: pioid_i
    integer                     :: ns_i, ns_o
    integer                     :: ni, no
    integer , allocatable       :: mask_i(:)
    real(r8), allocatable       :: frac_i(:)
    real(r8), allocatable       :: frac_o(:)
    real(r8), allocatable       :: soildepth_i(:)
    real(r8), allocatable       :: soildepth_o(:) ! output grid: fraction soildepth
    character(len=CS)           :: varname
    integer                     :: varnum
    real(r8), parameter         :: min_valid = 0._r8         ! minimum valid value
    real(r8), parameter         :: max_valid = 100.000001_r8 ! maximum valid value
    integer                     :: ier, rcode                ! error status
    character(len=*), parameter :: subname = 'mksoildepth'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make soildepth .....'
    end if

    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o and allocate output data
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate ( soildepth_o(ns_o)); soildepth_o(:) = spval

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(frac_o(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
         routehandle=routehandle, frac_o=frac_o, rc=rc)
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r4) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh and get frac_o
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
         routehandle=routehandle, frac_o=frac_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    do no = 1, ns_o
       if ((frac_o(no) < 0.0) .or. (frac_o(no) > 1.0001)) then
          write(6,*) "ERROR:: frac_o out of range: ", frac_o(no),no
          call shr_sys_abort ()
       end if
    end do

    ! Determine variable name to read in
    varnum = 1
    select case (varnum)
    case(1)
       varname = 'Avg_Depth_Median'
    case(2)
       varname = 'Avg_Depth_Mean'
    case(3)
       varname = 'Upland_Valley_Depth_Median'
    case(4)
       varname = 'Upland_Valley_Depth_Mean'
    case(5)
       varname = 'Upland_Hillslope_Depth_Median'
    case(6)
       varname = 'Upland_Hillslope_Depth_Mean'
    case(7)
       varname = 'Lowland_Depth_Mean'
    case(8)
       varname = 'Lowland_Depth_Mean'
    end select

    ! Read in input soil depth data
    allocate(soildepth_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, trim(varname), mesh_i, soildepth_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid soildepth_i to soildepth_o
    call regrid_rawdata(mesh_i, mesh_o, routehandle, soildepth_i, soildepth_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//'after regrid rawdata in '//trim(subname))

    ! Check validity of output data
    if ( min_bad(soildepth_o, min_valid, 'soildepth') .or. &
         max_bad(soildepth_o, max_valid, 'soildepth')) then
       call shr_sys_abort()
    end if

    ! Write output data
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil depth"
    call mkfile_output(pioid_o,  mesh_o,  'zbedrock', soildepth_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call pio_syncfile(pioid_o)

    ! Output diagnostic info
    call output_diagnostics_area(mesh_i, mesh_o, mask_i, frac_o, &
       soildepth_i, soildepth_o, "Soildepth", percent=.false., ndiag=ndiag, rc=rc)

    ! Close the input file
    call pio_closefile(pioid_i)
    call ESMF_VMLogMemInfo("After pio_closefile in "//trim(subname))

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made soildepth'
       write (ndiag,*)
    end if

  end subroutine mksoildepth

end module mksoildepthMod
