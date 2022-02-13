module mksoilfmaxMod

  !-----------------------------------------------------------------------
  ! Make soil fmax
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarpar
  use mkvarctl

  implicit none
  private           ! By default make data private

  public :: mksoilfmax         ! Make percent fmax

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoilfmax(file_mesh_i, file_data_i, mesh_o, fmax_o, rc)
    !
    ! make percent fmax
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    real(r8)          , intent(inout) :: fmax_o(:)   ! output grid: %fmax
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: n,l,k
    integer , allocatable  :: mask_i(:) 
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: fmax_i(:)             ! input grid: percent fmax
    real(r8), allocatable  :: areas_i(:)            ! input mesh areas
    real(r8), allocatable  :: areas_o(:)            ! output mesh areas
    real(r8)               :: sum_fldi              ! global sum of dummy input fld
    real(r8)               :: sum_fldo              ! global sum of dummy output fld
    real(r8)               :: gfmax_i               ! input  grid: global fmax
    real(r8)               :: garea_i               ! input  grid: global area
    real(r8)               :: gfmax_o               ! output grid: global fmax
    real(r8)               :: garea_o               ! output grid: global area
    integer                :: ier, rcode            ! error status
    character(len=32)      :: subname = 'mksoilfmax'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make %fmax .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Error check soil_fmax if it is set
    if ( soil_fmax_override /= unsetsoil )then
       if ( soil_fmax_override < 0.0 .or. soil_fmax_override > 1.0 )then
          write(6,*)'soil_fmax is out of range = ', soil_fmax_override
          call shr_sys_abort()
       end if
       write(6,*) 'Replace soil fmax for all points with: ', soil_fmax_override
    end if

    ! Open input data file
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
       if (frac_i(ni) > 0._r8) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    allocate(frac_o(ns_o))
    call create_routehandle_r8(mesh_i, mesh_o, routehandle, frac_o=frac_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))
    do n = 1, ns_o
       if ((frac_o(n) < 0.0) .or. (frac_o(n) > 1.0001)) then
          write(6,*) "ERROR:: frac_o out of range: ", frac_o(n),n
       end if
    end do

    ! Read in input data
    allocate(fmax_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'FMAX', mesh_i, fmax_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Regrid fmax_i to fmax_o, in points with no data, use globalAvg
    fmax_i(:) = fmax_i(:) * frac_i(:)
    fmax_o(:) = 0._r8
    call regrid_rawdata(mesh_i, mesh_o, routehandle, fmax_i, fmax_o, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//'after regrid rawdata in '//trim(subname))
    do n = 1,ns_o
       if (frac_o(n) == 0._r8) then
          fmax_o(n) = .365783_r8
       end if
    end do

    ! Check for conservation
    do no = 1, ns_o
       if ((fmax_o(no)) > 1.000001_r8) then
          write (6,*) 'MKFMAX error: fmax = ',fmax_o(no),' greater than 1.000001 for no = ',no
          call shr_sys_abort()
       end if
    enddo

    ! -----------------------------------------------------------------
    ! Error check2
    ! Compare global areas on input and output grids
    ! -----------------------------------------------------------------

!     allocate(areas_i(ns_i))
!     call get_meshareas(mesh_i, areas_i, rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return

!     allocate(areas_o(ns_o))
!     call get_meshareas(mesh_o, areas_o, rc)
!     if (ChkErr(rc,__LINE__,u_FILE_u)) return

!     gfmax_i = 0.
!     garea_i = 0.
!     do ni = 1,ns_i
!        garea_i = garea_i + areas_i(ni)*re**2
!        gfmax_i = gfmax_i + fmax_i(ni)*(areas_i(ni)/100.) * frac_i(ni)*re**2
!     end do

!     gfmax_o = 0.
!     garea_o = 0.
!     do no = 1,ns_o
!        garea_o = garea_o + areas_o(no)*re**2
!        gfmax_o = gfmax_o + fmax_o(no)*(areas_o(no)/100.) * frac_o(no)*re**2
!     end do
!     deallocate(areas_i)
!     deallocate(areas_o)

!     ! Diagnostic output
!     if (root_task) then
!        write (ndiag,*)
!        write (ndiag,'(1x,70a1)') ('=',k=1,70)
!        write (ndiag,*) 'Maximum Fractional Saturated Area Output'
!        write (ndiag,'(1x,70a1)') ('=',k=1,70)

!        write (ndiag,*)
!        write (ndiag,'(1x,70a1)') ('.',k=1,70)
!        write (ndiag,2001)
!        write (ndiag,'(1x,70a1)') ('.',k=1,70)
!        write (ndiag,*)
!        write (ndiag,2002) gfmax_i*1.e-06,gfmax_o*1.e-06
!        write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
! 2001   format (1x,'surface type   input grid area  output grid area'/&
!             1x,'                 10**6 km**2      10**6 km**2   ')
! 2002   format (1x,'fmax        ',f14.3,f17.3)
! 2004   format (1x,'all surface ',f14.3,f17.3)
!     end if

    ! Close the input file
    call pio_closefile(pioid)
    call ESMF_VMLogMemInfo("After pio_closefile in "//trim(subname))

    ! Release memory
    deallocate (fmax_i)
    deallocate (frac_i)
    deallocate (frac_o)
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %fmax'
       write (ndiag,*)
    end if

  end subroutine mksoilfmax

end module mksoilfmaxMod
