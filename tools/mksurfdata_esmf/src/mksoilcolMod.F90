module mksoilcolMod

  use ESMF
  use pio
  use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_sys_mod  , only : shr_sys_abort
  use mkpioMod     , only : mkpio_get_rawdata
  use mkpioMod     , only : mkpio_iodesc_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod    , only : regrid_rawdata, create_routehandle_r4, get_meshareas
  use mkutilsMod   , only : chkerr
  use mkvarctl     , only : root_task, ndiag, mpicom, MPI_INTEGER, MPI_MAX
  use mkvarctl     , only : soil_color_override, unsetcol

  implicit none
  private

#include <mpif.h>

  public  :: mksoilcol      ! Set soil colors
  private :: mkrank

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoilcol(file_data_i, file_mesh_i, mesh_o, soil_color_o, nsoilcol, rc)

    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i     ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i     ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o          ! model mesho
    integer           , intent(inout) :: soil_color_o(:) ! soil color classes
    integer           , intent(out)   :: nsoilcol        ! number of soil colors
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: n,l,k
    integer , allocatable  :: mask_i(:)
    real(r4), allocatable  :: frac_i(:)
    real(r4), allocatable  :: frac_o(:)
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    real(r4), allocatable  :: data_i(:,:)
    real(r4), allocatable  :: data_o(:,:)
    real(r4), allocatable  :: soil_color_i(:)
    logical                :: has_color     ! whether this grid cell has non-zero color
    integer                :: nsoilcol_local
    integer                :: maxindex(1)
    real(r4), allocatable  :: loc_gast_i(:) ! local global area, by surface type
    real(r4), allocatable  :: loc_gast_o(:) ! local global area, by surface type
    real(r4), allocatable  :: gast_i(:)     ! global area, by surface type
    real(r4), allocatable  :: gast_o(:)     ! global area, by surface type
    integer                :: rcode, ier
    character(len=*), parameter :: subname = 'mksoilcol'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Error check soil_color if it is set
    if ( soil_color_override /= unsetcol )then
       if ( soil_color_override < 0 .or. soil_color_override > 20 )then
          write(6,*)'soil_color is out of range = ', soil_color_override
          call shr_sys_abort()
       end if
       write(6,*) 'Replace soil color for all points with: ', soil_color_override
       do no = 1,size(soil_color_o)
          soil_color_o(no) = soil_color_override
       end do
       RETURN
    end if

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make soil color classes .....'
    end if

    ! Open soil color data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)

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
       if (frac_i(ni) > 0._r4) then
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
    call create_routehandle_r4(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Read in input soil color data
    allocate(soil_color_i(ns_i),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'SOIL_COLOR', mesh_i, soil_color_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Determine maximum number of soil colors across all processors
    ! This will be used for the ungridded dimension of data_o below
    nsoilcol_local = maxval(soil_color_i)
    call mpi_allreduce(nsoilcol_local, nsoilcol, 1, MPI_INTEGER, MPI_MAX, mpicom, rcode)
    ! TODO: MPI_SUCCESS could not be accessed from mkvarctl
    !if (rcode /= MPI_SUCCESS) call shr_sys_abort('error for mpi_allredice in '//trim(subname))

    ! Now determine data_i as a real 2d array - for every possible soil color create a global
    ! field with gridcells equal to 1 for that soil color and zero elsewhere
    allocate(data_i(0:nsoilcol,ns_i))
    data_i(:,:) = 0._r4
    do l = 0,nsoilcol
       do n = 1,ns_i
          if (int(soil_color_i(n)) == l) then
             data_i(l,n) = 1._r4 * mask_i(n)
          end if
       end do
    end do

    ! Regrid data_i to data_o
    allocate(data_o(0:nsoilcol, ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, nsoilcol, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//'after regrid rawdata in '//trim(subname))

    ! Determine output soil color
    soil_color_o(:) = 0
    do no = 1,ns_o
       ! If the output cell has any non-zero-colored inputs, then set the weight of
       ! zero-colored inputs to 0, to ensure that the zero-color is NOT dominant.
       if (any(data_o(1:nsoilcol,no) > 0.)) then
          has_color = .true.
          data_o(0,no) = 0.0
       else
          has_color = .false.
       end if

       ! Rank non-zero weights by color type. wsti(1) is the most extensive color type.
       if (has_color) then
          call mkrank (nsoilcol, data_o(0:nsoilcol,no), maxindex)
          soil_color_o(no) = maxindex(1)
       end if

       ! If land but no color, set color to 15 (in older dataset generic soil color 4)
       if (nsoilcol == 8) then
          if (soil_color_o(no)==0) then
             soil_color_o(no) = 4
          end if
       else if (nsoilcol == 20) then
          if (soil_color_o(no)==0) then
             soil_color_o(no) = 15
          end if
       else
          write(6,*) 'MKSOILCOL error: unhandled nsoilcol: ', nsoilcol
          call shr_sys_abort()
       end if

       ! Error checks
       if (soil_color_o(no) < 0 .or. soil_color_o(no) > nsoilcol) then
          write (6,*) 'MKSOILCOL error: land model soil color = ', &
               soil_color_o(no),' is not valid for lon,lat = ',no
          call shr_sys_abort()
       end if
    end do

    ! Compare global area of each soil color on input and output grids

    allocate(area_i(ns_i))
    allocate(area_o(ns_o))
    call get_meshareas(mesh_i, area_i, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call get_meshareas(mesh_o, area_o, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

!     allocate(loc_gast_i(0:nsoilcol))
!     allocate(loc_gast_o(0:nsoilcol))
!     allocate(gast_i(0:nsoilcol))
!     allocate(gast_o(0:nsoilcol))
!     loc_gast_i(:) = 0.
!     do ni = 1,ns_i
!        k = soil_color_i(ni)
!        loc_gast_i(k) = loc_gast_i(k) + area_i(ni) * mask_i(ni)
!     end do
!     call mpi_reduce(loc_gast_i, gast_i, nsoilcol+1, nsoilcol+1, MPI_REAL4, MPI_SUM, 0, mpicom, ier)
!     loc_gast_o(:) = 0.
!     do no = 1,ns_o
!        k = soil_color_o(no)
!        loc_gast_o(k) = loc_gast_o(k) + area_o(no) * frac_o(no)
!     end do
!     call mpi_reduce(loc_gast_o, gast_o, nsoilcol+1, nsoilcol+1, MPI_REAL4, MPI_SUM, 0, mpicom, ier)

!     write (ndiag,*)
!     write (ndiag,'(1x,70a1)') ('.',k=1,70)
!     write (ndiag,101)
! 101 format (1x,'soil color type',5x,' input grid area output grid area',/ &
!             1x,20x,'     10**6 km**2','      10**6 km**2')
!     write (ndiag,'(1x,70a1)') ('.',k=1,70)
!     write (ndiag,*)
!     do k = 0, nsoilcol
!        write (ndiag,'(1x,a,i3,d16.3,d17.3)') 'class ',k, gast_i(k)*1.e-6, gast_o(k)*1.e-6
!     end do

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made soil color classes'
       write (ndiag,'(a)')
    end if

    ! Clean up memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    call ESMF_LogWrite(subname//' finished routine mksoilcol')

  end subroutine mksoilcol

  !===============================================================
  subroutine mkrank (n, a, iv)
    !
    ! Return indices of largest [num] values in array [a].
    !
    ! input/output variables
    integer , intent(in) :: n      !array length
    real(r4), intent(in) :: a(0:n) !array to be ranked
    integer , intent(out):: iv(1)  !index to [num] largest values in array [a]

    ! local variables:
    real(r4) :: a_max  !maximum value in array
    real(r4) :: delmax !tolerance for finding if larger value
    integer  :: i      !array index
    integer  :: m      !do loop index
    integer  :: k      !do loop index
    integer  :: miss   !missing data value
    !-----------------------------------------------------------------------

    ! Find index of largest non-zero number
    delmax = 1.e-06
    miss = 9999
    iv(1) = miss

    a_max = -9999.
    do i = 0, n
       if (a(i)>0. .and. (a(i)-a_max)>delmax) then
          a_max = a(i)
          iv(1)  = i
       end if
    end do
    ! iv(1) = miss indicates no values > 0. this is an error
    if (iv(1) == miss) then
       write (6,*) 'MKRANK error: iv(1) = missing'
       call shr_sys_abort()
    end if
  end subroutine mkrank

end module mksoilcolMod
