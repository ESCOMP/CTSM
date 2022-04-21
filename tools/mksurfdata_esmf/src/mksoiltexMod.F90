module mksoiltexMod

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

  public :: mksoiltex      ! Set soil texture

  integer                :: mapunit_value_max
  integer                :: num_soil_textures
  type(ESMF_DynamicMask) :: dynamicMask

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoiltex(file_mesh_i, file_mapunit_i, file_lookup_i, mesh_o, pioid_o, pctlnd_pft_o, rc)
    !
    ! make %sand and %clay from IGBP soil data, which includes
    ! igbp soil 'mapunits' and their corresponding textures
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i     ! input mesh/grid file name
    character(len=*)  , intent(in)    :: file_mapunit_i  ! input mapunit file name
    character(len=*)  , intent(in)    :: file_lookup_i   ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o          ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
    real(r8)          , intent(in)    :: pctlnd_pft_o(:) ! PFT data: % of gridcell for PFTs
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Grid)        :: grid_i
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    type(file_desc_t)      :: pioid_i
    type(var_desc_t)       :: pio_varid
    integer                :: pio_vartype
    integer                :: dimid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: k,l,m,n
    integer                :: nlay          ! number of soil layers
    integer                :: n_scid
    integer , allocatable  :: mask_i(:)
    real(r4), pointer      :: dataptr(:)
    integer                :: mapunit       ! temporary igbp soil mapunit
    integer,  allocatable  :: sand_i(:,:,:) ! input grid: percent sand
    integer,  allocatable  :: clay_i(:,:,:) ! input grid: percent clay
    real(r4), allocatable  :: sand_o(:,:)   ! % sand (output grid)
    real(r4), allocatable  :: clay_o(:,:)   ! % clay (output grid)
    integer                :: n_mapunits
    integer                :: lookup_index
    integer                :: SCID
    real(r4), allocatable  :: mapunit_i(:)  ! input grid: igbp soil mapunits
    integer , allocatable  :: mapunit_o(:)  ! output grid: igbp soil mapunits
    integer , allocatable  :: MapUnits(:)
    integer , allocatable  :: mapunit_lookup(:)
    type(var_desc_t)       :: pio_varid_sand
    type(var_desc_t)       :: pio_varid_clay
    integer                :: starts(3)     ! starting indices for reading lookup table
    integer                :: counts(3)     ! dimension counts for reading lookup table
    integer                :: srcTermProcessing_Value = 0
    integer                :: rcode, ier    ! error status
    character(len=*), parameter :: subname = 'mksoiltex'
    !-----------------------------------------------------------------------

    if (root_task) then
       write(ndiag,*)
       write(ndiag,'(1x,80a1)') ('=',k=1,80)
       write(ndiag,*)
       write(ndiag,'(a)') 'Attempting to make %sand and %clay .....'
       write(ndiag,'(a)') ' Input mapunit file is '//trim(file_mapunit_i)
       write(ndiag,'(a)') ' Input lookup table file is '//trim(file_lookup_i)
       write(ndiag,'(a)') ' Input mesh/grid file is '//trim(file_mesh_i)
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

    ! Determine input mesh
    if (trim(file_mesh_i) == trim(file_mapunit_i)) then
       ! input format is GRIDSPEC and read in grid and then create mesh
       if (root_task) write(ndiag,*)"reading grid_i and then creating mesh_i in "//trim(subname)
       call ESMF_VMLogMemInfo("Before create read in grid_i in "//trim(subname))
       grid_i = ESMF_GridCreate(filename=trim(file_mesh_i), &
            fileformat=ESMF_FILEFORMAT_GRIDSPEC, addCornerStagger=.true., addmask=.true., varname='MU', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("Before create mesh_i in "//trim(subname))
       mesh_i = esmf_meshcreate(grid_i, rc=rc) 
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMLogMemInfo("After create mesh_i in "//trim(subname))
    else
       ! Read in mesh directly 
       if (root_task) write(ndiag,*)"reading mesh_i directly in "//trim(subname)
       mesh_i = ESMF_MeshCreate(filename=trim(file_mesh_i), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Determine ns_i (use the distgrid to the number of elements)
    call ESMF_MeshGet(mesh_i, numOwnedElements=ns_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read in mapunit data
    if (root_task) write(ndiag,*)"Reading in mapunit data in "//trim(subname)
    call ESMF_VMLogMemInfo("Before pio_openfile for "//trim(file_mapunit_i))
    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_mapunit_i), pio_nowrite)
    allocate(mapunit_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid_i, 'MU', mesh_i, mapunit_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))
    call pio_closefile(pioid_i)

    ! Set mesh mask to zero where the mapunit values are 0
    if (root_task) write(ndiag,*)"Setting mask in mesh where mapunit data is 0 "//trim(subname)
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
    field_i = ESMF_FieldCreate(mesh_i, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_o = ESMF_FieldCreate(mesh_o, ESMF_TYPEKIND_R4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Create a route handle
    if (root_task) write(ndiag,*)" before route handle creation "//trim(subname)
    call ESMF_FieldRegridStore(field_i, field_o, routehandle=routehandle, &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After regridstore in "//trim(subname))
    if (root_task) write(ndiag,*)" after route handle creation "//trim(subname)

    ! Create a dynamic mask object
    ! The dynamic mask object further holds a pointer to the routine that will be called in order to
    ! handle dynamically masked elements - in this case its DynMaskProc (see below)
    if (root_task) write(ndiag,*)" before call to dynamic mask set creation "//trim(subname)
    call ESMF_DynamicMaskSetR4R8R4(dynamicMask, dynamicMaskRoutine=get_dominant_mapunit,  &
         handleAllElements=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (root_task) write(ndiag,*)" after call to dynamic mask set creation "//trim(subname)

    ! Determine values in field_i
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = real(mapunit_i(:), kind=r4)

    ! Determine mapunit_value_max (set it as a module variable so that it can be
    ! accessible to gen_dominant_mapunit) - this is needed in the dynamic mask routine
    mapunit_value_max = maxval(dataptr)

    ! Determine values in field_o
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dataptr(:) = 0._r4

    ! Determine mapunit_o
    call ESMF_FieldRegrid(field_i, field_o, routehandle=routehandle, dynamicMask=dynamicMask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_o, farrayptr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do no = 1,ns_o
       mapunit_o(no) = int(dataptr(no))
    end do

    do no = 1,ns_o
       if (mapunit_o(no) > mapunit_value_max) then
          write(6,*)'mapunit_o is out of bounds ',mapunit_o(no)
          ! call shr_sys_abort("mapunit_o is out of bounds")
       end if
    end do

    ! Write out mapunit_o
    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out mapunits "
    call mkfile_output(pioid_o,  mesh_o,  'mapunits', mapunit_o,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output for mapunits')
    call pio_syncfile(pioid_o)

    !---------------------------------
    ! Determine %sand and %clay on output grid - using above mapunits
    !---------------------------------

    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_lookup_i), pio_nowrite)

    rcode = pio_inq_dimid  (pioid_i, 'MapUnit', dimid)
    rcode = pio_inq_dimlen (pioid_i, dimid, n_mapunits)

    rcode = pio_inq_dimid  (pioid_i, 'soil_layer', dimid)
    rcode = pio_inq_dimlen (pioid_i, dimid, nlay)

    rcode = pio_inq_dimid  (pioid_i, 'SCID', dimid)
    rcode = pio_inq_dimlen (pioid_i, dimid, n_scid)

    ! Read In MapUnits from the input file
    allocate(MapUnits(n_mapunits), stat=ier)
    rcode = pio_inq_varid(pioid_i, 'MapUnit', pio_varid)
    rcode = pio_get_var(pioid_i, pio_varid, MapUnits)

    ! Determine the mapunit lookup index from the value of the MapUnit variable
    mapunit_value_max = maxval(MapUnits)
    allocate(mapunit_lookup(mapunit_value_max))
    mapunit_lookup(:) = -999
    do n = 1,size(MapUnits)
       mapunit_lookup(MapUnits(n)) = n
    end do

    allocate(sand_i(nlay,n_scid,n_mapunits))
    if (ier/=0) call shr_sys_abort()
    allocate(clay_i(nlay,n_scid,n_mapunits))
    if (ier/=0) call shr_sys_abort()

    ! Get dimensions from input file and allocate memory for sand_i and clay_i
    rcode = pio_inq_varid(pioid_i, 'PCT_SAND', pio_varid_sand)
    rcode = pio_inq_varid(pioid_i, 'PCT_CLAY', pio_varid_clay)

    rcode = pio_get_var(pioid_i, pio_varid_sand, sand_i)
    rcode = pio_get_var(pioid_i, pio_varid_clay, clay_i)

    do no = 1,ns_o

       if (pctlnd_pft_o(no) < 1.e-6_r8 .or. mapunit_o(no) == 0) then

          ! Adjust sand and clay be loam if pctlnd_pft is < 1.e-6 or mapunit is 0
          sand_o(no,:) = 43._r4
          clay_o(no,:) = 18._r4

       else

          ! Determine lookup_index
          lookup_index = mapunit_lookup(mapunit_o(no))

          ! Determine the top soil layer sand_o and clay_o 
          ! If its less than 0 search within the SCID array for the first index 
          ! that gives a pct sand that is greater than or  equal to 0
          ! Then determine the other soil layers sand_o
          sand_o(no,1) = float(sand_i(1,1,lookup_index))
          if (sand_o(no,1) < 0.) then
             do l = 2,n_scid
                if (float(sand_i(1,l,lookup_index)) >= 0.) then
                   sand_o(no,1) = float(sand_i(1,l,lookup_index))
                   exit
                end if
             end do
          end if
          if (sand_o(no,1) < 0.) then
             if (int(sand_o(no,1)) == -4) then
                write(6,'(a,i8)')'WARNING: changing sand_o from -4 to 99% at no = ',no
                sand_o(no,:) = 99._r4
             else
                write(6,'(a,i8,a,i8)')'WARNING: changing sand_o from ',int(sand_o(no,1)),' to 43 at no = ',no
                sand_o(no,:) = 43._r4
             end if
          end if
          do l = 2,nlay
             sand_o(no,l) = float(sand_i(l,1,lookup_index))
             if (sand_o(no,l) < 0. .and. l > 1) then
                sand_o(no,l) = sand_o(no,l-1)
             end if
          end do

          ! If its less than 0 search within the SCID array for the first index 
          ! that gives a pct clay that is greater than or  equal to 0
          ! Now determine the other soil layers clay_o
          clay_o(no,1) = float(clay_i(1,1,lookup_index))
          if (clay_o(no,1) < 0.) then
             do l = 2,n_scid
                if (float(clay_i(1,l,lookup_index)) >= 0.) then
                   clay_o(no,1) = float(clay_i(1,l,lookup_index))
                   exit
                end if
             end do
          end if
          if (clay_o(no,1) < 0.) then
             if (int(clay_o(no,1)) == -4) then
                write(6,'(a,i8)')'WARNING: changing clay_o from -4 to 1% at no = ',no
                clay_o(no,:) = 1._r4
             else
                write(6,'(a,i8,a,i8)')'WARNING: changing clay_o from ',int(clay_o(no,1)),' to 18 at no = ',no
                clay_o(no,:) = 18._r4
             end if
          end if
          if (clay_o(no,1) < 0.) then
             write(6,*)'ERROR: at no, lookup_index = ',no,lookup_index
             call shr_sys_abort('could not find a value >= 0 for clay_i') 
          end if
          do l = 2,nlay
             clay_o(no,l) = float(clay_i(l,1,lookup_index))
             if (clay_o(no,l) < 0. .and. l > 1) then
                clay_o(no,l) = clay_o(no,l-1)
             end if
          end do

       end if

    end do

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil percent sand"
    call mkfile_output(pioid_o,  mesh_o,  'PCT_SAND', sand_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil percent clay"
    call mkfile_output(pioid_o,  mesh_o,  'PCT_CLAY', clay_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error in calling mkfile_output')

    call pio_syncfile(pioid_o)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %sand and %clay'
    end if

  end subroutine mksoiltex

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

end module mksoiltexMod
