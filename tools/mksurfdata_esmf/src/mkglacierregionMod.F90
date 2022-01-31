module mkglacierregionMod

  !-----------------------------------------------------------------------
  ! make glacier region ID
  ! Regridding is done by finding the nearest neighbor source cell for each destination cell.
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_nn, get_meshareas
  use mkutilsMod     , only : chkerr
#ifdef TODO
  ! use mkdiagnosticsMod, only : output_diagnostics_index
#endif
  use mkchecksMod, only : min_bad
  use mkvarctl

  implicit none
  private

  public :: mkglacierregion  ! make glacier region ID

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkglacierregion(file_mesh_i, file_data_i, mesh_o, glacier_region_o, rc)
    !
    ! Make glacier region ID
    !
    ! input/output variables
    character(len=*)  , intent(in)  :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)  :: mesh_o      ! output mesh
    integer           , intent(out) :: glacier_region_o(:) ! glacier region
    integer           , intent(out) :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle ! nearest neighbor routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer , allocatable  :: glacier_region_i(:) ! glacier region on input grid
    real(r4), allocatable  :: frac_i(:)           ! input mask
    real(r4), allocatable  :: frac_o(:)           ! output fractions
    real(r4), allocatable  :: data_i(:) 
    real(r4), allocatable  :: data_o(:) 
    integer                :: ier, rcode          ! error status
    integer                :: max_region          ! max region ID
    character(len=*), parameter :: subname = 'mkglacierregion'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make glacier region .....'
    end if

    ! Open input data file
    if (root_task) then
       write (ndiag,'(a)') 'Opening glacier region raw data file: ', trim(file_data_i)
    end if
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Create a nearest neighbor route handle between the input and output mesh
    call create_routehandle_nn(mesh_i, mesh_o, routehandle, rc=rc)
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Determine ns_i and allocate glacier_region_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(glacier_region_i(ns_i), stat=ier)
    if (ier/=0) call abort()

    ! Determine ns_o (glacier_region_o has already been allocated)
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read in input data
    call mkpio_get_rawdata(pioid, 'GLACIER_REGION', mesh_i, glacier_region_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Confirm that no value of glacier_region is less than min_allowed.
    if (min_bad(glacier_region_i, 0, 'GLACIER_REGION')) then
       call shr_sys_abort()
    end if

    ! Convert to real4
    allocate(data_i(ns_i))
    do ni = 1,ns_i
       data_i(ni) = real(glacier_region_i(ni), kind=r4)
    end do
    allocate(data_o(ns_o))

    ! Regrid raw data - yse a nearest neighbor map here
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now convert back to integer
    do no = 1,ns_o
       glacier_region_o(no) = nint(data_o(no))
    end do

    ! call output_diagnostics_index(glacier_region_i, glacier_region_o, tgridmap, &
    !      'Glacier Region ID', 0, max_region, ndiag, mask_src=tdomain%mask, frac_dst=frac_dst)

    ! Close the input file
    call pio_closefile(pioid)
    call ESMF_VMLogMemInfo("After pio_closefile in "//trim(subname))

    ! Release memory
    deallocate (frac_i)
    deallocate (frac_o)
    deallocate(glacier_region_i)
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made glacier region'
       write (ndiag,*)
    end if

  end subroutine mkglacierregion

end module mkglacierregionMod
