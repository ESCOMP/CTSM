module mkgdpMod

  !-----------------------------------------------------------------------
  ! make GDP from input GDP data
  !-----------------------------------------------------------------------
  !
  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl       , only : ndiag, root_task, mpicom

  implicit none
  private

#include <mpif.h>

  public :: mkgdp            ! regrid gdp data

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkgdp(file_mesh_i, file_data_i, mesh_o, gdp_o, rc)
    !
    ! make GDP from input GDP data
    !
    use mkdiagnosticsMod, only : output_diagnostics_continuous
    use mkchecksMod, only : min_bad
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! input model mesh
    real(r8)          , intent(inout) :: gdp_o(:)    ! output grid: GDP (x1000 1995 US$ per capita)
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer                :: ni,no,k
    integer                :: ns_i, ns_o
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    real(r8), allocatable  :: gdp_i(:)             ! input grid: percent gdp
    real(r8)               :: sum_fldi             ! global sum of dummy input fld
    real(r8)               :: sum_fldo             ! global sum of dummy output fld
    real(r8)               :: gglac_i              ! input  grid: global glac
    real(r8)               :: garea_i              ! input  grid: global area
    real(r8)               :: gglac_o              ! output grid: global glac
    real(r8)               :: garea_o              ! output grid: global area
    integer                :: ier, rcode           ! error status
    real(r8), parameter    :: min_valid = 0._r8    ! minimum valid value
    character(len=*), parameter :: subname = 'mkgdp'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make GDP.....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, frac_i, rc=rc)
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

    ! Read in peat_i
    allocate(gdp_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'gdp', mesh_i, gdp_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Create a route handle between the input and output mesh and get frac_o
    allocate(frac_o(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Regrid gdp
    call regrid_rawdata(mesh_i, mesh_o, routehandle, gdp_i, gdp_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (min_bad(gdp_o, min_valid, 'gdp')) then
       call shr_sys_abort(subname//' error in reading gdp_i')
    end if

    ! Close the file
    call pio_closefile(pioid)
    call ESMF_VMLogMemInfo("After pio_closefile in "//trim(subname))

    ! Output diagnostic info
    allocate(area_i(ns_i))
    allocate(area_o(ns_o))
    call get_meshareas(mesh_i, area_i, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call get_meshareas(mesh_o, area_o, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! call  output_diagnostics_continuous(area_i, area_o, mask_i, frac_o, &
    !      gdp_i, gdp_o, "GDP", "x1000 US$ per capita", ndiag)

    ! Clean up memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made GDP'
       write (ndiag,'(a)')
    end if

  end subroutine mkgdp

end module mkgdpMod
