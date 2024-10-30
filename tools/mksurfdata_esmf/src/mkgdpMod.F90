module mkgdpMod

  !-----------------------------------------------------------------------
  ! make GDP from input GDP data
  !-----------------------------------------------------------------------

  use ESMF
  use pio              , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite, pio_syncfile
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use mkpioMod         , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkesmfMod        , only : regrid_rawdata, create_routehandle_r8
  use mkvarctl         , only : ndiag, root_task, mpicom, spval
  use mkdiagnosticsMod , only : output_diagnostics_continuous
  use mkchecksMod      , only : min_bad
  use mkutilsMod       , only : chkerr
  use mkfileMod        , only : mkfile_output 

  implicit none
  private

  public :: mkgdp            ! regrid gdp data

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkgdp(file_mesh_i, file_data_i, mesh_o, pioid_o, rc)
    !
    ! make GDP from input GDP data
    !
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
    real(r8), allocatable  :: gdp_i(:)          ! input grid: percent gdp
    real(r8), allocatable  :: gdp_o(:)          ! output grid: GDP (x1000 1995 US$ per capita)
    real(r8), parameter    :: min_valid = 0._r8 ! minimum valid value
    integer                :: ier, rcode        ! error status
    character(len=*), parameter :: subname = 'mkgdp'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make GDP.....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if
    call ESMF_VMLogMemInfo("At start of"//trim(subname))

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
    allocate (gdp_o(ns_o)); gdp_o(:) = spval

    ! Get the landmask from the file and reset the mesh mask based on that
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

    ! Read in gdp_i
    allocate(gdp_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'gdp', mesh_i, gdp_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Create a route handle between the input and output mesh and get frac_o
    allocate(frac_o(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
         routehandle=routehandle, frac_o=frac_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Regrid gdp
    call regrid_rawdata(mesh_i, mesh_o, routehandle, gdp_i, gdp_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (min_bad(gdp_o, min_valid, 'gdp')) then
       call shr_sys_abort(subname//' error in reading gdp_i')
    end if

    ! Close the file
    call pio_closefile(pioid_i)

    ! Write output data
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out gdp"
    call mkfile_output(pioid_o, mesh_o, 'gdp', gdp_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call pio_syncfile(pioid_o)

    ! Output diagnostic info
    call output_diagnostics_continuous(mesh_i, mesh_o, gdp_i, gdp_o, &
       "GDP", "x1000 US$ per capita", ndiag=ndiag, rc=rc, mask_i=mask_i, frac_o=frac_o)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    ! Clean up memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made GDP'
    end if
    call ESMF_VMLogMemInfo("At end of "//trim(subname))

  end subroutine mkgdp

end module mkgdpMod
