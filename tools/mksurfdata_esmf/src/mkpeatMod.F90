module mkpeatMod

  !-----------------------------------------------------------------------
  ! make fraction peat from input peat data
  !-----------------------------------------------------------------------
  !
  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkpioMod       , only : mkpio_iodesc_rawdata, mkpio_get_rawdata_level
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl       , only : ndiag, root_task, mpicom

  implicit none
  private

#include <mpif.h>

  public :: mkpeat           ! regrid peat data

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkpeat(file_mesh_i, file_data_i, mesh_o, peat_o, rc)

    use mkdiagnosticsMod, only : output_diagnostics_area
    use mkchecksMod     , only : min_bad, max_bad

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i   ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i   ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    real(r8)          , intent(inout) :: peat_o(:)     ! output grid: fraction peat
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
    real(r8), allocatable  :: peat_i(:)            ! input grid: percent glac
    real(r8)               :: sum_fldi             ! global sum of dummy input fld
    real(r8)               :: sum_fldo             ! global sum of dummy output fld
    real(r8)               :: gglac_i              ! input  grid: global glac
    real(r8)               :: garea_i              ! input  grid: global area
    real(r8)               :: gglac_o              ! output grid: global glac
    real(r8)               :: garea_o              ! output grid: global area
    integer                :: ier, rcode           ! error status
    real(r8)               :: relerr = 0.00001     ! max error: sum overlap wts ne 1
    real(r8), allocatable :: data_i(:)          ! data on input grid
    real(r8), parameter :: min_valid = 0._r8          ! minimum valid value
    real(r8), parameter :: max_valid = 100.000001_r8  ! maximum valid value
    character(len=*), parameter :: subname = 'mkpeat'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make peat .....'
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
    allocate(peat_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'peatf', mesh_i, peat_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Create a route handle between the input and output mesh and get frac_o
    allocate(frac_o(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Regrid peat
    call regrid_rawdata(mesh_i, mesh_o, routehandle, peat_i, peat_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (min_bad(peat_o, min_valid, 'peat') .or. max_bad(peat_o, max_valid, 'peat')) then
       call shr_sys_abort(subname//" peat_o does not fall in range of min_valid/max_valid")
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
    call output_diagnostics_area(area_i, area_o, mask_i, peat_i, peat_o, 'peat', .true., ndiag)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made organic matter '
    end if

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made peat'
    end if

  end subroutine mkpeat

end module mkpeatMod
