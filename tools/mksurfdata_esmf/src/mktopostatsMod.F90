module mktopostatsMod

  !-----------------------------------------------------------------------
  ! make various topography statistics
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, pio_iotype, pio_ioformat, pio_iosystem
  use mkutilsMod     , only : chkerr
  use mkvarctl       , only : ndiag, root_task, mpicom, std_elev, spval
  use mkfileMod      , only : mkfile_output  

  implicit none
  private

  public :: mktopostats            ! make topo stddev & mean slope

  type(ESMF_DynamicMask) :: dynamicMask

  logical :: calculate_stddev = .true.

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mktopostats(file_mesh_i, file_data_i, file_data_i_override, mesh_o, pioid_o, rc)

    ! make various topography statistics
    !
    use mkdiagnosticsMod, only : output_diagnostics_continuous
    use mkchecksMod     , only : min_bad, max_bad
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i          ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i          ! input data file name
    character(len=*)  , intent(in)    :: file_data_i_override ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o               ! input model mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc
    !
    ! local variables:
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    type(file_desc_t)      :: pioid_i
    integer                :: ni,no,k
    integer                :: ns_i, ns_o
    real(r4), allocatable  :: data_i(:)
    real(r4), pointer      :: dataptr(:)
    real(r4), allocatable  :: topo_stddev_o(:) ! output grid: standard deviation of elevation (m)
    real(r4), allocatable  :: slope_o(:)       ! output grid: slope (degrees)
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
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
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

    ! Determine ns_o and allocate output data
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate (topo_stddev_o(ns_o)) ; topo_stddev_o(:) = spval
    allocate (slope_o(ns_o))       ; slope_o(:)       = spval

    ! Read file_data_i_override for data that is assumed to already be on the output grid
    if (file_data_i_override /= ' ' ) then
       if (root_task)  write(ndiag, '(a)') trim(subname)//" reading STD_ELEV and SLOPE from "//trim(file_data_i_override)
       ! TODO: get dimensions and make sure that they match the dimensions of mesh_o
       rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i_override), pio_nowrite)
       call mkpio_get_rawdata(pioid_i, 'STD_ELEV', mesh_o, topo_stddev_o, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call mkpio_get_rawdata(pioid_i, 'SLOPE', mesh_o, slope_o, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (root_task)  write(ndiag, '(a)') trim(subname)//" writing topo_stddev "
       call mkfile_output(pioid_o,  mesh_o, 'STD_ELEV', topo_stddev_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for STD_ELEV')
       call mkfile_output(pioid_o,  mesh_o, 'SLOPE', slope_o, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for SLOPE')
       call pio_syncfile(pioid_o)
       RETURN
    end if

    ! Open input data file
    ! Read in data with PIO_IOTYPE_NETCDF rather than PIO_IOTYPE_PNETCDF since there are problems
    ! with the pnetcdf read of this high resolution data
    rcode = pio_openfile(pio_iosystem, pioid_i, PIO_IOTYPE_NETCDF, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

    ! Read in input mesh
    mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_UGRID, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))

    ! Determine ns_i
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read in input data data_i
    allocate(data_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//' error in allocating data_i')
    call mkpio_get_rawdata(pioid_i, 'ELEVATION', mesh_i, data_i, rc=rc)
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

    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    topo_stddev_o(:) = dataptr(:)

    ! Check validity of output data
    if (min_bad(topo_stddev_o, min_valid_topo_stddev, 'topo_stddev')) then
       call shr_sys_abort()
    end if
    call output_diagnostics_continuous(mesh_i, mesh_o, real(data_i,8), real(topo_stddev_o,8), &
         "Topo Std Dev", "m", ndiag=ndiag, rc=rc, nomask=.true.)

    ! -----------------------------
    ! Obtain the slope
    ! -----------------------------

    call mkpio_get_rawdata(pioid_i, 'SLOPE', mesh_i, data_i, rc=rc)
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

    ! Check validity of output data
    if (min_bad(slope_o, min_valid_slope, 'slope') .or. &
        max_bad(slope_o, max_valid_slope, 'slope')) then
       call shr_sys_abort()
    end if
    call output_diagnostics_continuous(mesh_i, mesh_o, real(data_i,8), real(slope_o,8), &
         "Slope", "degrees", ndiag=ndiag, rc=rc, nomask=.true.)

    ! Write out output data
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing topo_stddev "
    call mkfile_output(pioid_o,  mesh_o, 'STD_ELEV', topo_stddev_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for STD_ELEV')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing slope"
    call mkfile_output(pioid_o,  mesh_o, 'SLOPE', slope_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for SLOPE')
    call pio_syncfile(pioid_o)

    ! Close files and deallocate dynamic memory
    call pio_closefile(pioid_i)
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

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
