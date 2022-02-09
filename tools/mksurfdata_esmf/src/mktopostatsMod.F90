module mktopostatsMod

  !-----------------------------------------------------------------------
  ! make various topography statistics
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : regrid_rawdata, create_routehandle_r8, get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl       , only : ndiag, root_task, mpicom

#include <mpif.h>

  implicit none
  private

  public :: mktopostats            ! make topo stddev & mean slope

  type(ESMF_DynamicMask) :: dynamicMask

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mktopostats(file_mesh_i, file_data_i, mesh_o, std_elev, topo_stddev_o, slope_o, rc)

    ! make various topography statistics
    !
    use mkdiagnosticsMod, only : output_diagnostics_continuous, output_diagnostics_continuous_outonly
    use mkchecksMod     , only : min_bad, max_bad
    !
    ! input/output variables
    character(len=*)  , intent(in)  :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)  :: file_data_i      ! input data file name
    type(ESMF_Mesh)   , intent(in)  :: mesh_o           ! input model mesh
    real(r8)          , intent(out) :: topo_stddev_o(:) ! output grid: standard deviation of elevation (m)
    real(r8)          , intent(out) :: slope_o(:)       ! output grid: slope (degrees)
    real(r8)          , intent(in)  :: std_elev         ! standard deviation of elevation (m) to use when not using input file
    integer           , intent(out) :: rc

    character(len=*)  , intent(in) :: mapfname          ! input mapping file name
    character(len=*)  , intent(in) :: datfname          ! input data file name
    integer           , intent(in) :: ndiag             ! unit number for diag out
    !
    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid
    integer                :: ni,no,k
    integer                :: ns_i, ns_o
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: rmask_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    real(r8), allocatable  :: gdp_i(:)          ! input grid: percent gdp
    logical                :: bypass_reading    ! bypass reading dataset and just use a global value
    integer                :: ier, rcode        ! error status
    type(ESMF_Field)  :: field_i
    type(ESMF_Field)  :: field_o
    real(r8), pointer :: dataptr(:)
    real(r8), parameter    :: min_valid = 0._r8 ! minimum valid value
    real(r8), parameter    :: min_valid_topo_stddev = 0._r8
    real(r8), parameter    :: min_valid_slope = 0._r8
    real(r8), parameter    :: max_valid_slope = 90._r8
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
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ns_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Force the mask to be one
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    mask_i(:) = 1
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    deallocate(mask_i)

    ! Compute the standard deviation
    ! Use the dynamic masking could be used to do this.
    ! dst_array(no) = dst_array(no) + wt * (src_array(ni) - weighted_means(no))**2 part,
    ! and then just do a plain sparse matrix multiply on a src Field of all 1.0 to do the
    ! weight sum part, and then do the divide and sqrt locally on each PET.

    ! One issue is how to get the weight_means() in for each
    ! dest. location. I think that you could pass them in via the
    ! dst_array and then pull them out before doing the calculation, but
    ! to do so that you to stop it from zeroing out the
    ! dst array which you can do by setting
    ! zeroregion=ESMF_REGION_EMPTY.This would also depend on the
    ! calculation happening only once for each destination location,
    ! which I would guess is true, but Gerhard can confirm.

    ! Read in data_i
    allocate(data_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//' error in allocating data_i')
    call mkpio_get_rawdata(pioid, 'ELEVATION', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Create a route handle between the input and output mesh and get frac_o
    allocate(frac_o(ns_o),stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Create ESMF fields that will be used below
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, srcTermProcessing=0, &
         srcMaskValues=(/-987987/), dstMaskValues=(/-987987/), &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    ! -----------------------------
    ! Obtain the standard deviation
    ! -----------------------------

    ! Create a dynamic mask object
    ! The dynamic mask object further holds a pointer to the routine that will be called in order to
    ! handle dynamically masked elements - in this case its DynMaskProc (see below)
    call ESMF_DynamicMaskSetR8R8R8(dynamicMask, dynamicMaskRoutine=StdDevProc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Interpolate data_i to data_o
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = data_i(:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r8

    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, &
         zeroregion=ESMF_REGION_EMPTY, dynamicMaskRoutine=StdDevProc, rc=rc)

    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    topo_stddev_o(:) = dataptr(:)

    ! call output_diagnostics_continuous_outonly(topo_stddev_o, tgridmap, "Topo Std Dev", "m", ndiag)
    ! Check validity of output data
    if (min_bad(topo_stddev_o, min_valid_topo_stddev, 'topo_stddev')) then
       call shr_sys_abort()
    end if

    ! -----------------------------
    ! Obtain the slope
    ! -----------------------------

    call mkpio_get_rawdata(pioid, 'SLOPE', mesh_i, data_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = data_i(:)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r8

    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, &
         zeroregion=ESMF_REGION_EMPTY, rc=rc)

    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    slope_o(:) = dataptr(:)

    !call output_diagnostics_continuous(data_i, slope_o, tgridmap,
    !"Slope", "degrees", ndiag, tdomain%mask, tgridmap%frac_dst)

    ! Check validity of output data
    if (min_bad(slope_o, min_valid_slope, 'slope') .or. &
        max_bad(slope_o, max_valid_slope, 'slope')) then
       call shr_sys_abort()
    end if

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
    type(ESMF_DynamicMaskElementR8R8R8) , pointer              :: dynamicMaskList(:)
    real(ESMF_KIND_R8)                  , intent(in), optional :: dynamicSrcMaskValue
    real(ESMF_KIND_R8)                  , intent(in), optional :: dynamicDstMaskValue
    integer                             , intent(out)          :: rc

    ! local variables
    integer  :: i, j
    real(ESMF_KIND_R8)  :: renorm
    real(ESMF_KIND_R8)  :: mean
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Below - ONLY if you do NOT have the source masked out then do
    ! the regridding (which is done explicitly here)
    ! Below i are the destination points and j are the source points

    if (associated(dynamicMaskList)) then
       do i=1, size(dynamicMaskList)
          dynamicMaskList(i)%dstElement = czero ! set to zero
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
          mean = dynamicMaskList(i)%dstElement

          ! Now compute the standard deviation
          dynamicMaskList(i)%dstElement = czero ! reset to zero
          do j = 1, size(dynamicMaskList(i)%factor)
             dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement + &
                  (dynamicMaskList(i)%factor(j) * (dynamicMaskList(i)%srcElement(j)) - mean)**2
          enddo
       enddo
    endif

  end subroutine StdDevProc

end module mktopostatsMod
