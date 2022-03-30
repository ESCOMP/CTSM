module mksoiltexnewMod

  !-----------------------------------------------------------------------
  ! Make soil data (texture)
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod     , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod      , only : shr_sys_abort
  use mkpioMod         , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod         , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkutilsMod       , only : chkerr
  use mkdiagnosticsMod , only : output_diagnostics_index
  use mkfileMod        , only : mkfile_output  
  use mkvarctl         , only : root_task, ndiag, spval
  use mkvarctl         , only : unsetsoil
  use mkvarpar         , only : nlevsoi

  implicit none
  private ! By default make data private

  public :: mksoiltexnew  ! Set soil texture

  integer                :: mapunit_value_max
  integer                :: num_soil_textures
  type(ESMF_DynamicMask) :: dynamicMask

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoiltexnew(file_mesh_i, file_data_i, mesh_o, pioid_o, pctlnd_pft_o, rc)
    !
    ! make %sand and %clay from IGBP soil data, which includes
    ! igbp soil 'mapunits' and their corresponding textures
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    real(r8)          , intent(in)    :: pctlnd_pft_o(:) ! PFT data: % of gridcell for PFTs
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Grid)        :: grid_i
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    type(ESMF_Field)       :: field_dstfrac
    type(file_desc_t)      :: pioid_i
    type(var_desc_t)       :: pio_varid
    integer                :: pio_vartype
    integer                :: dimid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: k,l,m,n
    integer                :: nlay         ! number of soil layers
    integer , allocatable  :: mask_i(:)
    real(r4), allocatable  :: sand_i(:,:)  ! input grid: percent sand
    real(r4), allocatable  :: clay_i(:,:)  ! input grid: percent clay
    real(r8), allocatable  :: sand_o(:,:)  ! % sand (output grid)
    real(r8), allocatable  :: clay_o(:,:)  ! % clay (output grid)
    integer , allocatable  :: mapunit_o(:)
    real(r4), pointer      :: dataptr(:)
    integer                :: mapunit      ! temporary igbp soil mapunit
    !
    integer                :: n_mapunits
    integer                :: lookup_index
    real(r4), allocatable  :: mapunit_i(:) ! input grid: igbp soil mapunits
    integer , allocatable  :: MapUnits(:)
    integer , allocatable  :: mapunit_lookup(:)
    type(var_desc_t)       :: pio_varid_sand
    type(var_desc_t)       :: pio_varid_clay
    integer                :: starts(3)    ! starting indices for reading lookup table
    integer                :: counts(3)    ! dimension counts for reading lookup table
    !
    integer                :: srcTermProcessing_Value = 0
    integer                :: rcode, ier   ! error status
    character(len=*), parameter :: subname = 'mksoiltexnew'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make %sand and %clay .....'
       write(ndiag,'(a)') ' Input file is '//trim(file_data_i)
    end if

    ! Determine ns_o and allocate output data
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(mapunit_o(ns_o))      ; mapunit_o(:) = 0
    allocate(sand_o(ns_o,nlevsoi)) ; sand_o(:,:) = spval
    allocate(clay_o(ns_o,nlevsoi)) ; clay_o(:,:) = spval

    !---------------------------------
    ! Determine mapunits on output grid
    !---------------------------------

    ! Open input data file
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_data_i))
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_data_i), pio_nowrite)

    ! Read in input mapunit mesh
    call ESMF_VMLogMemInfo("Before create grid_i in "//trim(subname))
    if (root_task) write(ndiag,*)"DEBUG: before create grid_i in "//trim(subname)
    grid_i = ESMF_GridCreate(filename=trim(file_data_i), &
         fileformat=ESMF_FILEFORMAT_GRIDSPEC, addCornerStagger=.true., addmask=.true., varname='MU', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (root_task) write(ndiag,*)"DEBUG: before create mesh_i in "//trim(subname)
    call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
    mesh_i = esmf_meshcreate(grid_i, rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))
    if (root_task) write(ndiag,*)"DEBUG: after create mesh_i in "//trim(subname)

    ! Determine ns_i (use the distgrid to the number of elements)
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read in mapunit data
    if (root_task) write(ndiag,*)"DEBUG: before mapunit read in "//trim(subname)
    allocate(mapunit_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'MU', mesh_i, mapunit_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))
    if (root_task) write(ndiag,*)"DEBUG: after mapunit read in "//trim(subname)

    ! Set mapunit values to zero where the input mask is 0
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    do ni = 1,ns_i
       if (mapunit_i(ni) == 0.) then
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

     ! Create ESMF fields that will be used below
    if (root_task) write(ndiag,*)"DEBUG: before field_i creation "//trim(subname)
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (root_task) write(ndiag,*)"DEBUG: before field_o creation "//trim(subname)
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_dstfrac = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle
    if (root_task) write(ndiag,*)"DEBUG: before route handle creation "//trim(subname)
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
         dstFracField= field_dstfrac, rc=rc)
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))
    if (root_task) write(ndiag,*)"DEBUG: after route handle creation "//trim(subname)

    ! Create a dynamic mask object
    ! The dynamic mask object further holds a pointer to the routine that will be called in order to
    ! handle dynamically masked elements - in this case its DynMaskProc (see below)
    if (root_task) write(ndiag,*)"DEBUG: before call to dynamic mask set creation "//trim(subname)
    call ESMF_DynamicMaskSetR4R8R4(dynamicMask, dynamicMaskRoutine=get_dominant_mapunit,  &
         handleAllElements=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (root_task) write(ndiag,*)"DEBUG: after call to dynamic mask set creation "//trim(subname)

    ! Determine dominant soil color in the field regrid call below
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = real(mapunit_i(:), kind=r4)

    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r4

    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, dynamicMask=dynamicMask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       mapunit_o(no) = int(dataptr(no))
    end do

    ! Write out mapunits
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out mapunits" 
    call mkfile_output(pioid_o,  mesh_o,  'mapunits', mapunit_o,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for mapunits')

    !---------------------------------
    ! Determine %sand and %clay on output grid - using above mapunits
    !---------------------------------

    ! rcode = pio_inq_dimid  (pioid_i, 'MapUnit', dimid)
    ! rcode = pio_inq_dimlen (pioid_i, dimid, n_mapunits)

    ! rcode = pio_inq_dimid  (pioid_i, 'soil_layer', dimid)
    ! rcode = pio_inq_dimlen (pioid_i, dimid, nlay)

    ! ! Read In MapUnits from the input file
    ! allocate(MapUnits(n_mapunits), stat=ier)
    ! rcode = pio_inq_varid(pioid_i, 'MapUnit', pio_varid)
    ! rcode = pio_get_var(pioid_i, pio_varid, MapUnits)

    ! ! Determine the mapunit lookup index from the value of the MapUnit variable
    ! mapunit_value_max = maxval(MapUnits)
    ! allocate(mapunit_lookup(mapunit_value_max))
    ! mapunit_lookup(:) = -999
    ! do n = 1,size(MapUnits)
    !    mapunit_lookup(MapUnits(n)) = n
    ! end do

    ! allocate(sand_i(mapunit_value_max,nlay), stat=ier)
    ! if (ier/=0) call shr_sys_abort()
    ! allocate(clay_i(mapunit_value_max,nlay), stat=ier)
    ! if (ier/=0) call shr_sys_abort()

    ! ! Get dimensions from input file and allocate memory for sand_i and clay_i

    ! rcode = pio_inq_varid(pioid_i, 'PCT_SAND', pio_varid_sand)
    ! rcode = pio_inq_varid(pioid_i, 'PCT_CLAY', pio_varid_clay)
    ! starts(1:3) = 1
    ! counts(1) = n_mapunits
    ! counts(2) = 1
    ! counts(3) = nlay

    ! allocate(sand_i(n_mapunits,nlay))
    ! allocate(clay_i(n_mapunits,nlay))
    ! rcode = pio_get_var(pioid_i, pio_varid_sand, starts, counts, sand_i)
    ! rcode = pio_get_var(pioid_i, pio_varid_clay, starts, counts, clay_i)

    ! do no = 1,ns_o
    !    if (mapunit_o(no) > 0) then
    !       ! valid value is obtained
    !       if (mapunit_o(no) > mapunit_value_max) then
    !          write(6,*)'mapunit_o is out of bounds ',mapunit_o(no)
    !          ! call shr_sys_abort("mapunit_o is out of bounds")
    !       end if
    !       lookup_index = mapunit_lookup(mapunit_o(no))
    !       do l = 1, nlay
    !          sand_o(no,l) = sand_i(lookup_index,l)
    !          clay_o(no,l) = clay_i(lookup_index,l)
    !       end do
    !    else
    !       ! use loam
    !       do l = 1, nlay
    !          sand_o(no,l) = 43.
    !          clay_o(no,l) = 18.
    !       end do
    !    end if
    ! end do

    ! ! Adjust pct sand and pct clay to be nearest integers and to be loam if pctlnd_pft is < 1.e-6
    ! ! Truncate all percentage fields on output grid. This is needed to insure that wt is zero
    ! ! (not a very small number such as 1e-16) where it really should be zero
    ! do no = 1,ns_o
    !    do k = 1,nlevsoi
    !       sand_o(no,k) = float(nint(sand_o(no,k)))
    !       clay_o(no,k) = float(nint(clay_o(no,k)))
    !    end do
    !    if (pctlnd_pft_o(no) < 1.e-6_r8) then
    !       sand_o(no,:) = 43._r8
    !       clay_o(no,:) = 18._r8
    !    end if
    ! end do

    ! if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil percent sand"
    ! call mkfile_output(pioid_o,  mesh_o,  'PCT_SAND', sand_o, lev1name='nlevsoi', rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

    ! if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil percent clay"
    ! call mkfile_output(pioid_o,  mesh_o,  'PCT_CLAY', clay_o, lev1name='nlevsoi', rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')
    ! call pio_syncfile(pioid_o)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %sand and %clay'
    end if

  end subroutine mksoiltexnew

  !================================================================================================
  subroutine get_dominant_mapunit(dynamicMaskList, dynamicSrcMaskValue, dynamicDstMaskValue, rc)

    use ESMF, only : ESMF_RC_ARG_BAD

    ! input/output arguments
    type(ESMF_DynamicMaskElementR4R8R4) , pointer              :: dynamicMaskList(:)
    real(ESMF_KIND_R4)                  , intent(in), optional :: dynamicSrcMaskValue
    real(ESMF_KIND_R4)                  , intent(in), optional :: dynamicDstMaskValue
    integer                             , intent(out)          :: rc

    ! local variables
    integer            :: ni, no, n
    real(ESMF_KIND_R4) :: wts_o(0:mapunit_value_max)
    integer            :: maxindex(1)
    real(ESMF_KIND_R4) :: maxvalue
    character(len=*), parameter :: subname = 'mksoiltex'
    !---------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (associated(dynamicMaskList)) then
       do no = 1, size(dynamicMaskList)
          dynamicMaskList(no)%dstElement = 0.d0
          wts_o(:) = 0.d0
          do ni = 1, size(dynamicMaskList(no)%factor)
             if (dynamicMaskList(no)%srcElement(ni) > 0.d0) then
                do n = 0,mapunit_value_max
                   if (dynamicMaskList(no)%srcElement(ni) == n) then
                      wts_o(n) = wts_o(n) + dynamicMaskList(no)%factor(ni)
                   end if
                enddo
             end if
          end do

          ! Determine the most dominant index of wts_o
          maxvalue = -999._r4
          maxindex = -999
          do n = 0,mapunit_value_max
             if (wts_o(n) > maxvalue) then
                maxindex(1) = n
                maxvalue = wts_o(n)
             end if
          end do
          if (maxindex(1) > mapunit_value_max) then
             write(6,*)'mapunit_o is out of bounds ',maxindex(1)
             call shr_sys_abort(subname//" mapunit_o is out of bounds")
          end if
          dynamicMaskList(no)%dstElement = real(maxindex(1), kind=r4)
       end do
    end if

  end subroutine get_dominant_mapunit

end module mksoiltexnewMod
