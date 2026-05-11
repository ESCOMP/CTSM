module mkglacierregionMod

  !-----------------------------------------------------------------------
  ! make glacier region ID
  ! Regridding is done by finding the nearest neighbor source cell for each destination cell.
  !-----------------------------------------------------------------------

  use ESMF
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use pio              , only : file_desc_t, pio_openfile, pio_closefile, pio_nowrite
  use mkpioMod         , only : mkpio_get_rawdata, pio_iotype, pio_iosystem
  use mkesmfMod        , only : regrid_rawdata, create_routehandle_r8
  use mkchecksMod      , only : min_bad
  use mkvarctl         , only : ndiag, root_task
  use mkdiagnosticsMod , only : output_diagnostics_index
  use mkutilsMod       , only : chkerr
  use mkfileMod        , only : mkfile_output

  implicit none
  private

  public :: mkglacierregion  ! make glacier region ID

  integer, private :: nglacier_regions = 3

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mkglacierregion(file_mesh_i, file_data_i, mesh_o, pioid_o, rc)
    !
    ! Make glacier region ID
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    integer           , intent(out)   :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: routehandle ! nearest neighbor routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(file_desc_t)      :: pioid_i
    integer                :: ni,no,l,k
    integer                :: ns_i, ns_o
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: frac_i(:)
    real(r8), allocatable  :: frac_o(:)
    real(r8), allocatable  :: data_i(:,:)
    real(r8), allocatable  :: data_o(:,:)
    integer , allocatable  :: glacier_region_i(:) ! glacier region on input grid
    integer , allocatable  :: glacier_region_o(:) ! glacier region on output grid
    integer                :: ier, rcode          ! error status
    character(len=*), parameter :: subname = 'mkglacierregion'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make glacier region .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
       write(ndiag,'(a)') ' Input mesh file is '//trim(file_mesh_i)
    end if
    call ESMF_VMLogMemInfo("At start of "//trim(subname))

    ! Open input data file
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

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
    allocate (glacier_region_o(ns_o)) ; glacier_region_o(:) = -999

    ! Read in input data (Confirm that no value of glacier_region is less than min_allowed)
    allocate(glacier_region_i(ns_i), stat=ier)
    if (ier/=0) call abort()
    call mkpio_get_rawdata(pioid_i, 'GLACIER_REGION', mesh_i, glacier_region_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (min_bad(glacier_region_i, 0, 'GLACIER_REGION')) then
       call shr_sys_abort()
    end if
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Reset mesh mask to zero where glacier_region_i is zero
    allocate(frac_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    do ni = 1,ns_i
       if (glacier_region_i(ni) == 0) then
          mask_i(ni) = 0
       else
          mask_i(ni) = 1
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle between the input and output mesh
    allocate(frac_o(ns_o))
    if (ier/=0) call abort()
    call create_routehandle_r8(mesh_i=mesh_i, mesh_o=mesh_o, norm_by_fracs=.true., &
         routehandle=routehandle, frac_o=frac_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create routehandle in "//trim(subname))

    ! Now determine data_i as a real 2d array - for every possible glacier region create a global
    ! field with gridcells equal to 1 for that region and zero elsewhere
    allocate(data_i(0:nglacier_regions,ns_i))
    data_i(:,:) = 0._r4
    do l = 0,nglacier_regions
       do ni = 1,ns_i
          if (glacier_region_i(ni) == l) then
             data_i(l,ni) = 1._r4
          end if
       end do
    end do

    ! Regrid data_i to data_o
    allocate(data_o(0:nglacier_regions, ns_o), stat=ier)
    if (ier/=0) call shr_sys_abort('error allocating data_i(max_regions, ns_o)')
    call regrid_rawdata(mesh_i, mesh_o, routehandle, data_i, data_o, 0, nglacier_regions, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine glacier_region_o
    !
    ! We take the maximum glacier region index that has > 0 coverage in the output grid
    ! cell. The frequency of occurrence of the different glacier regions is irrelevant.
    ! For example, if the value 1 appears in 99% of source cells overlapping a given
    ! destination cell and the value 2 appears in just 1%, we'll put 2 in the destination
    ! cell because it is the maximum value. (This treatment is important so that we give
    ! priority to non-zero indices, and more generally allow setting up the glacier region
    ! indices so that higher values take precedence - i.e., if a CTSM grid cell has any
    ! overlap with a higher-valued region, we use that region.)
    glacier_region_o(:) = 0
    do no = 1,ns_o
       ! Note that we loop starting at the highest index so we give priority to higher
       ! indices. Also note that we stop looking at index 1 (i.e., we don't explicitly
       ! look at index 0): if we haven't found a region with greater than 0 coverage from
       ! looking at the indices greater than 0, then we'll default to assigning this to
       ! glacier region 0, regardless of the coverage of glacier region 0 (though, in
       ! practice, glacier region 0 should have 100% coverage in this case).
       do l = nglacier_regions, 1, -1
          if (data_o(l,no) > 0._r8) then
             glacier_region_o(no) = l
             exit
          end if
       end do
    end do
    
    ! Write output data
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out glacier_region"
    call mkfile_output(pioid_o,  mesh_o, 'GLACIER_REGION', glacier_region_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

    ! Determine global diagnostics
    call output_diagnostics_index(mesh_i, mesh_o, mask_i, frac_o, &
         0, 3, glacier_region_i, glacier_region_o, 'Glacier Region ID', ndiag, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    ! Close the input file
    call pio_closefile(pioid_i)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made glacier region'
    end if
    call ESMF_VMLogMemInfo("At end of "//trim(subname))

  end subroutine mkglacierregion

end module mkglacierregionMod
