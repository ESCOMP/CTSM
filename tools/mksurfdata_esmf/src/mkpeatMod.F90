module mkpeatMod

  !-----------------------------------------------------------------------
  ! make fraction peat from input peat data
  !-----------------------------------------------------------------------

  use ESMF
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use pio              , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite, pio_syncfile
  use mkpioMod         , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkesmfMod        , only : regrid_rawdata, create_routehandle_r8
  use mkvarctl         , only : ndiag, root_task, mpicom, spval
  use mkchecksMod      , only : min_bad, max_bad
  use mkdiagnosticsMod , only : output_diagnostics_area
  use mkutilsMod       , only : chkerr
  use mkfileMod        , only : mkfile_output

  implicit none
  private

  public :: mkpeat

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkpeat(file_mesh_i, file_data_i, mesh_o, pioid_o, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! input model mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid_i
    integer                :: ni,no,k
    integer                :: ns_i, ns_o
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: peat_i(:)              ! input grid: percent peat
    real(r8), allocatable  :: peat_o(:)   ! output grid: fraction peat
    integer                :: ier, rcode             ! error status
    real(r8), parameter    :: min_valid = 0._r8         ! minimum valid value
    real(r8), parameter    :: max_valid = 100.000001_r8 ! maximum valid value
    character(len=*), parameter :: subname = 'mkpeat'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make peat .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if
    call ESMF_VMLogMemInfo("At start of "//trim(subname))

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o and allocate output data
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate (peat_o(ns_o)) ; peat_o(:) = spval

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, frac_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (frac_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Read in peat_i
    allocate(peat_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'peatf', mesh_i, peat_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Create a route handle between the input and output mesh and get frac_o
    allocate(frac_o(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
         routehandle=routehandle, frac_o=frac_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Regrid peat
    call regrid_rawdata(mesh_i, mesh_o, routehandle, peat_i, peat_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (min_bad(peat_o, min_valid, 'peat') .or. max_bad(peat_o, max_valid, 'peat')) then
       call shr_sys_abort(subname//" peat_o does not fall in range of min_valid/max_valid")
    end if

    ! Write out data to output file
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out peatland fraction"
    call mkfile_output(pioid_o, mesh_o, 'peatf', peat_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call pio_syncfile(pioid_o)

    ! Output diagnostic info
    call output_diagnostics_area(mesh_i, mesh_o, mask_i, frac_o, &
         peat_i, peat_o, "Peat", percent=.false., ndiag=ndiag, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    ! Close the file
    call pio_closefile(pioid_i)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made peat'
    end if
    call ESMF_VMLogMemInfo("At end of "//trim(subname))

  end subroutine mkpeat

end module mkpeatMod
