module mkorganicMod

  ! make organic matter dataset

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod     , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata, create_routehandle_r8
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag
  use mkvarpar     , only : nlevsoi

  implicit none
  private

  public mkorganic      ! Set organic soil

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkorganic(file_mesh_i, file_data_i, mesh_o, organic_o, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i      ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o
    real(r8)          , intent(out)   :: organic_o(:,:)   ! output grid: %lake
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: nlay
    integer                :: n, l  ! indices
    integer                :: rcode, ier             ! error status
    integer                :: ndims
    integer , allocatable  :: dimlengths(:)
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    character(len=*), parameter :: subname = 'mkorganic'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write (ndiag,'(a)') ' Attempting to make organic mater dataset .....'
    end if

    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Get dimensions of raw data. 
    !  - raw data is dimensions (lon,lat,lev) 
    !  - input read from pio has dimensions(n,lev)
    !  - esmf field dataptr has dimensions (lev,n)
    allocate(dimlengths(3))
    call mkpio_get_dimlengths(pioid, 'ORGANIC', ndims, dimlengths)
    nlay = dimlengths(3)
    if (nlay /= nlevsoi) then
       write(6,*)'nlay, nlevsoi= ',nlay,nlevsoi,' do not match'
       call shr_sys_abort()
    end if
    deallocate(dimlengths)

    ! Read in input mesh
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
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
       if (frac_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh and get frac_o
    allocate(frac_o(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Read in input data
    ! - levels are the innermost dimension for esmf fields
    ! - levels are the outermost dimension in pio reads
    ! Input data is read into (ns_i,nlay) array and then transferred to data_i(nlay,ns_i)
    allocate(data_i(nlay,ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'ORGANIC', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid data_i to data_o
    allocate(data_o(nlay,ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 1, nlay, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regrid organic_i in  "//trim(subname))

    ! Set organic_o
    do l = 1,nlevsoi
       do no = 1,ns_o
          organic_o(no,l) = data_o(l,no)
       end do
    end do
    do l = 1,nlevsoi
       do no = 1,ns_o
          if ((organic_o(no,l)) > 130.000001_r8) then
             write (6,*) trim(subname)//' error: organic = ',organic_o(no,l), &
                  ' greater than 130.000001 for n,lev ',no,l
             call shr_sys_abort()
          end if
       enddo
    end do

    ! Close the file 
    call pio_closefile(pioid)
    call ESMF_VMLogMemInfo("After pio_closefile in "//trim(subname))

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a,i8)') 'Successfully made organic matter '
    end if

  end subroutine mkorganic

end module mkorganicMod
