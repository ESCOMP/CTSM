module mktopostatsMod

  !-----------------------------------------------------------------------
  ! make various topography statistics
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl       , only : ndiag, root_task, mpicom, std_elev

  implicit none
  private

#include <mpif.h>

  public :: mktopostats            ! make topo stddev & mean slope

  type(ESMF_DynamicMask) :: dynamicMask

  logical :: calculate_stddev = .true.

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mktopostats(file_mesh_i, file_data_i, mesh_o, topo_stddev_o, slope_o, rc)

    ! make various topography statistics
    !
    use mkdiagnosticsMod, only : output_diagnostics_continuous, output_diagnostics_continuous_outonly
    use mkchecksMod     , only : min_bad, max_bad
    !
    ! input/output variables
    character(len=*)  , intent(in)  :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i      ! input data file name
    type(ESMF_Mesh)   , intent(in)  :: mesh_o           ! input model mesh
    real(r4)          , intent(out) :: topo_stddev_o(:) ! output grid: standard deviation of elevation (m)
    real(r4)          , intent(out) :: slope_o(:)       ! output grid: slope (degrees)
    integer           , intent(out) :: rc
    !
    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    type(file_desc_t)      :: pioid
    integer                :: ni,no,k
    integer                :: ns_i, ns_o
    real(r4), allocatable  :: data_i(:)
    real(r4), pointer      :: dataptr(:)
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    integer                :: ier, rcode        ! error status
    integer                :: srcTermProcessing_Value = 0
    real(r4), parameter    :: min_valid = 0._r4 ! minimum valid value
    real(r4), parameter    :: min_valid_topo_stddev = 0._r4
    real(r4), parameter    :: min_valid_slope = 0._r4
    real(r4), parameter    :: max_valid_slope = 90._r4
    character(len=*), parameter :: subname = 'mktopostats'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make Topography statistics.....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    if ( std_elev >= 0.0_r8 )then
       if (root_task) then
          write (ndiag,'(a)')' Bypass the reading and just use global values'
          write (ndiag,'(a)')' Setting std deviation of topography to ', std_elev
          write (ndiag,'(a)')' Setting slope of topography to zero'
       end if
       topo_stddev_o(:) = std_elev
       slope_o = 0.0_r8
       RETURN
    end if

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_UGRID, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read in data_i
    allocate(data_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//' error in allocating data_i')
    call mkpio_get_rawdata(pioid, 'ELEVATION', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Create ESMF fields that will be used below
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle
    if (root_task) then
       write(ndiag,'(a)') subname//' creating a routehandle '
    end if
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (root_task) then
       write(ndiag,'(a)') subname//' finished creating a routehandle '
    end if
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    ! -----------------------------
    ! Obtain the standard deviation
    ! -----------------------------

    ! Create a dynamic mask object
    ! The dynamic mask object further holds a pointer to the routine that will be called in order to
    ! handle dynamically masked elements - in this case its DynMaskProc (see below)
    call ESMF_DynamicMaskSetR4R8R4(dynamicMask, dynamicMaskRoutine=StdDevProc, handleAllElements=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Interpolate data_i to data_o
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = data_i(:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r4

    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, dynamicMask=dynamicMask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, rc=rc)
    ! if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    topo_stddev_o(:) = dataptr(:)

    ! call output_diagnostics_continuous_outonly(topo_stddev_o, tgridmap, "Topo Std Dev", "m", ndiag)
    ! Check validity of output data
    !if (min_bad(topo_stddev_o, min_valid_topo_stddev, 'topo_stddev')) then
       !call shr_sys_abort()
    !end if

    ! TODO: get the output diagnostics working
    ! call  output_diagnostics_continuous_outonly(area_i, area_o, mask_i, frac_o, &
    !      data_i, topo_stddev_o, "Topo Std Dev", "m", ndiag)

    ! -----------------------------
    ! Obtain the slope
    ! -----------------------------

    call mkpio_get_rawdata(pioid, 'SLOPE', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = data_i(:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r4

    calculate_stddev = .false.  ! module variable used by dynamic mask
    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, dynamicMask=dynamicMask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    slope_o(:) = dataptr(:)

    !call output_diagnostics_continuous(data_i, slope_o, tgridmap,
    !"Slope", "degrees", ndiag, tdomain%mask, tgridmap%frac_dst)

    ! Check validity of output data
    !if (min_bad(slope_o, min_valid_slope, 'slope') .or. &
       ! max_bad(slope_o, max_valid_slope, 'slope')) then
       !call shr_sys_abort()
    !end if

    ! Close files and deallocate dynamic memory

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made Topography statistics'
       write (ndiag,'(a)')
    end if

  end subroutine mktopostats

  !================================================================================================
  subroutine StdDevProc(dynamicMaskList, dynamicSrcMaskValue, dynamicDstMaskValue, rc)

    use ESMF, only : ESMF_RC_ARG_BAD

    ! input/output arguments
    type(ESMF_DynamicMaskElementR4R8R4) , pointer              :: dynamicMaskList(:)
    real(ESMF_KIND_R4)                  , intent(in), optional :: dynamicSrcMaskValue
    real(ESMF_KIND_R4)                  , intent(in), optional :: dynamicDstMaskValue
    integer                             , intent(out)          :: rc

    ! local variables
    integer  :: i, j
    real(ESMF_KIND_R4)  :: renorm
    real(ESMF_KIND_R4)  :: mean
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Below - ONLY if you do NOT have the source masked out then do
    ! the regridding (which is done explicitly here)
    ! Below i are the destination points and j are the source points

    if (associated(dynamicMaskList)) then
       do i=1, size(dynamicMaskList)
          dynamicMaskList(i)%dstElement = 0.d0 ! set to zero
          renorm = 0.d0 ! reset

          ! Determine the mean
          do j = 1, size(dynamicMaskList(i)%factor)
             dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement + &
                  (dynamicMaskList(i)%factor(j) * dynamicMaskList(i)%srcElement(j))
             renorm = renorm + dynamicMaskList(i)%factor(j)
          enddo
          if (renorm > 0.d0) then
             dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement / renorm
          else
             rc = ESMF_RC_ARG_BAD  ! error detected
             return
          endif

          ! Now compute the standard deviation
          if (calculate_stddev) then
             mean = dynamicMaskList(i)%dstElement
             dynamicMaskList(i)%dstElement = 0.d0 ! reset to zero
             do j = 1, size(dynamicMaskList(i)%factor)
                dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement + &
                     (dynamicMaskList(i)%factor(j) * (dynamicMaskList(i)%srcElement(j) - mean)**2)
             enddo
             dynamicMaskList(i)%dstElement = sqrt(dynamicMaskList(i)%dstElement/renorm)
          end if
       enddo
    endif

  end subroutine StdDevProc

end module mktopostatsMod
