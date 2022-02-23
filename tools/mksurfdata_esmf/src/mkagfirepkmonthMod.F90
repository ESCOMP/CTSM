module mkagfirepkmonthMod

  !-----------------------------------------------------------------------
  ! Make agricultural fire peak month data
  !-----------------------------------------------------------------------

  use ESMF
  use pio              , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite, pio_syncfile
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use mkpioMod         , only : mkpio_get_rawdata,  pio_iotype, pio_iosystem
  use mkvarctl         , only : ndiag, root_task
  use mkchecksMod      , only : min_bad, max_bad
  use mkdiagnosticsMod , only : output_diagnostics_index
  use mkutilsMod       , only : chkerr
  use mkfileMod        , only : mkfile_output

  implicit none
  private           ! By default make data private

  public  :: mkagfirepkmon       ! Set agricultural fire peak month

  integer , parameter :: min_valid = 1
  integer , parameter :: max_valid = 12
  integer , parameter :: unsetmon = 13

  type(ESMF_DynamicMask) :: dynamicMask

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================
contains
!===============================================================

  subroutine mkagfirepkmon(file_mesh_i, file_data_i, mesh_o, pioid_o, rc)
    !
    ! Make agricultural fire peak month data from higher resolution data
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i      ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i      ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o           ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc
    !
    ! local variables:
    type(ESMF_RouteHandle)         :: routehandle
    type(ESMF_Mesh)                :: mesh_i
    type(ESMF_Field)               :: field_i
    type(ESMF_Field)               :: field_o
    type(ESMF_Field)               :: field_dstfrac
    type(file_desc_t)              :: pioid_i
    integer                        :: k
    integer                        :: ni,no
    integer                        :: ns_i, ns_o
    integer , allocatable          :: mask_i(:)
    real(r4), allocatable          :: rmask_i(:)
    real(r8), allocatable          :: frac_o(:)
    integer , allocatable          :: idata_i(:)    ! input grid: agricultural fire peak month
    integer , allocatable          :: agfirepkmon_o(:) ! agricultural fire peak month
    real(r4), pointer              :: dataptr(:)
    real(r8), pointer              :: dataptr_r8(:)
    integer                        :: rcode, ier    ! error status
    integer                        :: srcTermProcessing_Value = 0
    character(len=*), parameter    :: subname = 'mkagfirepkmon'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make agricultural fire peak month data .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)
    call ESMF_VMLogMemInfo("After pio_openfile "//trim(file_data_i))

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
    allocate (agfirepkmon_o(ns_o)); agfirepkmon_o(:) = -999

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(rmask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" ERROR in allocating rmask_i")
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" ERROR in allocating mask_i")
    call mkpio_get_rawdata(pioid_i, 'LANDMASK', mesh_i, rmask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (rmask_i(ni) > 0.) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Read in agfirepkmon_i
    allocate(idata_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort(subname//" error in allocating idata_i")
    call mkpio_get_rawdata(pioid_i, 'abm', mesh_i, idata_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

     ! Create ESMF fields that will be used below
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_dstfrac = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         dstFracField= field_dstfrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))

    ! Determine frac_o
    call ESMF_FieldGet(field_dstfrac, farrayptr=dataptr_r8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(frac_o(ns_o))
    frac_o(:) = dataptr_r8(:)

    ! Create a dynamic mask object
    ! The dynamic mask object further holds a pointer to the routine that will be called in order to
    ! handle dynamically masked elements - in this case its DynMaskProc (see below)
    call ESMF_DynamicMaskSetR4R8R4(dynamicMask, dynamicMaskRoutine=get_dominant_indices,  &
         handleAllElements=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine dominant fire month
    ! **NOTE** the use of the dynamicMask argument to the ESMF_FieldRegrid call
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = real(idata_i(:), kind=r4)
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r4
    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, dynamicMask=dynamicMask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       agfirepkmon_o(no) = int(dataptr(no))
    end do

    ! Check validity of output data
    if (min_bad(agfirepkmon_o, min_valid, 'agfirepkmon') .or. &
        max_bad(agfirepkmon_o, unsetmon , 'agfirepkmon')) then
        call shr_sys_abort()
     end if

    ! Close the file 
    call pio_closefile(pioid_i)

    ! Write out data
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing abm (agricultural fire peak month)"
    call mkfile_output(pioid_o,  mesh_o, 'abm', agfirepkmon_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    call pio_syncfile(pioid_o)

    ! Output diagnostics comparing global area of each peak month on input and output grids
    call output_diagnostics_index(mesh_i, mesh_o, mask_i, frac_o, &
         1, 13, idata_i, agfirepkmon_o, 'peak fire month', ndiag, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    
    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_o, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_FieldDestroy(field_dstfrac, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made Agricultural fire peak month'
       write (ndiag,*)
    end if

  end subroutine mkagfirepkmon

  !================================================================================================
  subroutine get_dominant_indices(dynamicMaskList, dynamicSrcMaskValue, dynamicDstMaskValue, rc)

    ! input/output arguments
    type(ESMF_DynamicMaskElementR4R8R4) , pointer              :: dynamicMaskList(:)
    real(ESMF_KIND_R4)                  , intent(in), optional :: dynamicSrcMaskValue
    real(ESMF_KIND_R4)                  , intent(in), optional :: dynamicDstMaskValue
    integer                             , intent(out)          :: rc

    ! local variables
    integer            :: ni, no, n
    real(ESMF_KIND_R4) :: wts_o(min_valid:max_valid)
    integer            :: maxindex(1)
    logical            :: hasdata 
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (associated(dynamicMaskList)) then
       do no = 1, size(dynamicMaskList)
          hasdata = .false.
          wts_o(:) = 0.d0
          do ni = 1, size(dynamicMaskList(no)%factor)
             if (dynamicMaskList(no)%srcElement(ni) > 0.d0) then
                do n = min_valid,max_valid
                   if ( dynamicMaskList(no)%srcElement(ni) == n) then
                      wts_o(n) = wts_o(n) + dynamicMaskList(no)%factor(ni)
                      hasdata = .true.
                   end if
                enddo
             end if
          end do
          
          ! Determine the most dominant index of wts_o
          if (hasdata) then
             maxindex = maxloc(wts_o(:)) 
             dynamicMaskList(no)%dstElement = real(maxindex(1), kind=r4)
          else
             dynamicMaskList(no)%dstElement = real(unsetmon, kind=r4)
          end if
       end do
    end if

  end subroutine get_dominant_indices

end module mkagfirepkmonthMod
