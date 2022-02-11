module mksoiltexMod

  !-----------------------------------------------------------------------
  ! Make soil data (texture)
  !-----------------------------------------------------------------------

  use ESMF
  use pio
  use shr_kind_mod   , only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod    , only : shr_sys_abort
  use mkpioMod       , only : mkpio_get_rawdata, mkpio_get_dimlengths
  use mkpioMod       , only : pio_iotype, pio_ioformat, pio_iosystem
  use mkesmfMod      , only : get_meshareas
  use mkutilsMod     , only : chkerr
  use mkvarctl
  use mkvarpar

  implicit none
  private ! By default make data private

  public :: mksoiltex      ! Set soil texture

  integer, parameter :: num=2  ! set soil mapunit number
  integer, parameter :: nlsm=4 ! number of soil textures

  integer :: num_soil_textures
  type(ESMF_DynamicMask) :: dynamicMask

  integer :: mapunit_value_max

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoiltex(file_mesh_i, file_data_i, mesh_o, sand_o, clay_o, mapunit_o, rc)
    !
    ! make %sand and %clay from IGBP soil data, which includes
    ! igbp soil 'mapunits' and their corresponding textures
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i ! input mesh file name
    character(len=*)  , intent(in)    :: file_data_i ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o      ! output mesh
    real(r8)          , intent(inout) :: sand_o(:,:) ! % sand (output grid)
    real(r8)          , intent(inout) :: clay_o(:,:) ! % clay (output grid)
    integer           , intent(inout) :: mapunit_o(:)
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_RouteHandle) :: routehandle
    type(ESMF_Mesh)        :: mesh_i
    type(ESMF_Field)       :: field_i
    type(ESMF_Field)       :: field_o
    type(ESMF_Field)       :: field_dstfrac
    type(file_desc_t)      :: pioid
    type(var_desc_t)       :: pio_varid
    integer                :: pio_vartype
    integer                :: dimid
    integer                :: ni,no
    integer                :: ns_i, ns_o
    integer                :: n,l,m
    character(len=38)      :: typ                     ! soil texture based on ...
    integer                :: nlay                    ! number of soil layers
    integer                :: mapunittemp             ! temporary igbp soil mapunit
    integer                :: maxovr
    integer , allocatable  :: mask_i(:)
    real(r8), allocatable  :: area_i(:)
    real(r8), allocatable  :: area_o(:)
    real(r4), allocatable  :: rmask_i(:)
    real(r4), allocatable  :: frac_o(:)
    real(r4), allocatable  :: sand_i(:,:)             ! input grid: percent sand
    real(r4), allocatable  :: clay_i(:,:)             ! input grid: percent clay
    real(r4), allocatable  :: mapunit_i(:)            ! input grid: igbp soil mapunits
    real(r4), pointer      :: dataptr(:)
    real(r8), pointer      :: dataptr_r8(:)
    character(len=38)      :: soil(0:nlsm)            ! name of each soil texture
    real(r4)               :: gast_i(0:nlsm)          ! global area, by texture type
    real(r4)               :: gast_o(0:nlsm)          ! global area, by texture type
    real(r4)               :: sum_fldi                ! global sum of dummy input fld
    real(r4)               :: sum_fldo                ! global sum of dummy output fld
    real(r4)               :: sumtex
    integer                :: rcode, ier              ! error status
    integer                :: srcTermProcessing_Value = 0
    character(len=*), parameter :: subname = 'mksoiltex'
    !-----------------------------------------------------------------------

    if (root_task) then
       write (ndiag,'(a)') 'Attempting to make %sand and %clay .....'
    end if

    if ( soil_clay_override /= unsetsoil )then
       write(6,*) 'Replace soil clay % for all points with: ', soil_clay_override
       if ( soil_sand_override == unsetsoil )then
          write (6,*) subname//':error: soil_clay set, but NOT soil_sand'
          call shr_sys_abort()
       end if
    end if
    if ( soil_sand_override /= unsetsoil )then
       write(6,*) 'Replace soil sand % for all points with: ', soil_sand_override
       if ( soil_clay_override == unsetsoil )then
          write (6,*) subname//':error: soil_sand set, but NOT soil_clay'
          call shr_sys_abort()
       end if
       sumtex = soil_sand_override + soil_clay_override
       if ( sumtex < 0.0_r4 .or. sumtex > 100.0_r4 )then
          write (6,*) subname//':error: soil_sand and soil_clay out of bounds: sand, clay = ', &
               soil_sand_override, soil_clay_override
          call shr_sys_abort()
       end if
    end if

    if (soil_sand_override /= unsetsoil .and. soil_clay_override /= unsetsoil) then
       if (root_task) then
          write(ndiag,'(a,i8)') ' Overriding soil color for all points with: ', soil_color_override
       end if
       sand_o(:,:) = soil_sand_override
       clay_o(:,:) = soil_clay_override
       RETURN
    end if

    ! Define the model surface types: 0:4
    soil(0) = 'no soil: ocean, glacier, lake, no data'
    soil(1) = 'clays                                 '
    soil(2) = 'sands                                 '
    soil(3) = 'loams                                 '
    soil(4) = 'silts                                 '

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

    ! Determine ns_o and allocate data_o
    call ESMF_MeshGet(mesh_o, numOwnedElements=ns_o, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the landmask from the file and reset the mesh mask based on that
    allocate(rmask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(mask_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'LANDMASK', mesh_i, rmask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do ni = 1,ns_i
       if (rmask_i(ni) > 0._r4) then
          mask_i(ni) = 1
       else
          mask_i(ni) = 0
       end if
    end do
    call ESMF_MeshSet(mesh_i, elementMask=mask_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Read in mapunit data
    allocate(mapunit_i(ns_i), stat=ier)
    if (ier/=0) call shr_sys_abort()
    call mkpio_get_rawdata(pioid, 'MAPUNITS', mesh_i, mapunit_i, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMLogMemInfo("After mkpio_getrawdata in "//trim(subname))

    ! Scale the input soil color by the input mask
    ! do ni = 1,ns_i
    !    if (mask_i(ni) == 0) then
    !       mapunit_i(ni) = 0._r4
    !    end if
    ! end do

    ! Determine mapunit_value_max (set it as a module variable so that it can be 
    ! accessible to gen_dominant_mapunit)
    rcode = pio_inq_dimid  (pioid, 'max_value_mapunit', dimid)
    rcode = pio_inq_dimlen (pioid, dimid, mapunit_value_max)

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

    ! Determin frac_o
    call ESMF_FieldGet(field_dstfrac, farrayptr=dataptr_r8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(frac_o(ns_o))
    frac_o(:) = real(dataptr_r8(:), kind=r4)

    ! Create a dynamic mask object
    ! The dynamic mask object further holds a pointer to the routine that will be called in order to
    ! handle dynamically masked elements - in this case its DynMaskProc (see below)
    call ESMF_DynamicMaskSetR4R8R4(dynamicMask, dynamicMaskRoutine=get_dominant_mapunit,  &
         handleAllElements=.true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine dominant soil color in the field regrid call below
    call ESMF_FieldGet(field_i, farrayptr=dataptr, rc=rc)
    dataptr(:) = mapunit_i(:)
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

    ! Get dimensions from input file and allocate memory for sand_i and clay_i
    rcode = pio_inq_dimid  (pioid, 'number_of_layers', dimid)
    rcode = pio_inq_dimlen (pioid, dimid, nlay)
    allocate(sand_i(mapunit_value_max,nlay), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(clay_i(mapunit_value_max,nlay), stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! read in sand_i and clay_i (they will read in total on all processors)
    rcode = pio_inq_varid(pioid, 'PCT_SAND', pio_varid)
    rcode = pio_get_var(pioid, pio_varid, sand_i)
    rcode = pio_inq_varid(pioid, 'PCT_CLAY', pio_varid)
    rcode = pio_get_var(pioid, pio_varid, clay_i)

    ! Set soil texture as follows:
    ! a. Use dominant igbp soil mapunit based on area of overlap unless 'no data' is dominant
    ! b. If this has no data, use loam for soil texture

    do no = 1,ns_o
       if (mapunit_o(no) > 0) then
          ! valid value is obtained
          if (mapunit_o(no) > mapunit_value_max) then
             write(6,*)'mapunit_o is out of bounds ',mapunit_o(no)
             ! call shr_sys_abort("mapunit_o is out of bounds")
          end if
          do l = 1, nlay
             sand_o(no,l) = sand_i(mapunit_o(no),l)
             clay_o(no,l) = clay_i(mapunit_o(no),l)
          end do
       else
          ! use loam
          do l = 1, nlay
             sand_o(no,l) = 43.
             clay_o(no,l) = 18.
          end do
       end if
    end do

    ! -----------------------------------------------------------------
    ! TODO:
    ! Compare global area of each soil type on input and output grids
    ! -----------------------------------------------------------------

    !     allocate(area_i(ns_i))
    !     call get_meshareas(mesh_i, area_i, rc)
    !     if (chkerr(rc,__LINE__,u_FILE_u)) return
    !     allocate(area_o(ns_o))
    !     call get_meshareas(mesh_o, area_o, rc)
    !     if (chkerr(rc,__LINE__,u_FILE_u)) return
    
    !     ! input grid: global areas by texture class
    !     gast_i(:) = 0.
    !     do l = 1, nlay
    !        do ni = 1,ns_i
    !           mapunittemp = nint(mapunit_i(ni))
    !           if (mapunittemp==0) then
    !              typ = 'no soil: ocean, glacier, lake, no data'
    !           else if (clay_i(mapunittemp,l) >= 40.) then
    !              typ = 'clays'
    !           else if (sand_i(mapunittemp,l) >= 50.) then
    !              typ = 'sands'
    !           else if (clay_i(mapunittemp,l)+sand_i(mapunittemp,l) < 50.) then
    !              if (rmask_i(ni) /= 0.) then
    !                 typ = 'silts'
    !              else            !if (mask(ni) == 0.) then no data
    !                 typ = 'no soil: ocean, glacier, lake, no data'
    !              end if
    !           else
    !              typ = 'loams'
    !           end if
    !           do m = 0, nlsm
    !              if (typ == soil(m)) go to 101
    !           end do
    !           write (6,*) 'MKSOILTEX error: sand = ',sand_i(mapunittemp,l), &
    !                ' clay = ',clay_i(mapunittemp,l), &
    !                ' not assigned to soil type for input grid lon,lat,layer = ',ni,l
    !           call shr_sys_abort()
    ! 101       continue
    !           gast_i(m) = gast_i(m) + area_i(ni)*mask_i(ni)*re**2
    !        end do
    !     end do
    
    !     ! output grid: global areas by texture class
    !     gast_o(:) = 0.
    !     do l = 1, nlay
    !        do no = 1,ns_o
    !           if (clay_o(no,l)==0. .and. sand_o(no,l)==0.) then
    !              typ = 'no soil: ocean, glacier, lake, no data'
    !           else if (clay_o(no,l) >= 40.) then
    !              typ = 'clays'
    !           else if (sand_o(no,l) >= 50.) then
    !              typ = 'sands'
    !           else if (clay_o(no,l)+sand_o(no,l) < 50.) then
    !              typ = 'silts'
    !           else
    !              typ = 'loams'
    !           end if
    !           do m = 0, nlsm
    !              if (typ == soil(m)) go to 102
    !           end do
    !           write (6,*) 'MKSOILTEX error: sand = ',sand_o(no,l), &
    !                ' clay = ',clay_o(no,l), &
    !                ' not assigned to soil type for output grid lon,lat,layer = ',no,l
    !           call shr_sys_abort()
    ! 102       continue
    !           gast_o(m) = gast_o(m) + area_o(no)*frac_o(no)*re**2
    !        end do
    !     end do
    
    !     ! Diagnostic output
    
    !     write (ndiag,*)
    !     write (ndiag,'(1x,70a1)') ('=',l=1,70)
    !     write (ndiag,*) 'Soil Texture Output'
    !     write (ndiag,'(1x,70a1)') ('=',l=1,70)
    !     write (ndiag,*)
    
    !     write (ndiag,*) 'The following table of soil texture classes is for comparison only.'
    !     write (ndiag,*) 'The actual data is continuous %sand, %silt and %clay not textural classes'
    !     write (ndiag,*)
    
    !     write (ndiag,*)
    !     write (ndiag,'(1x,70a1)') ('.',l=1,70)
    !     write (ndiag,1001)
    ! 1001 format (1x,'soil texture class',17x,' input grid area output grid area',/ &
    !              1x,33x,'     10**6 km**2','      10**6 km**2')
    !     write (ndiag,'(1x,70a1)') ('.',l=1,70)
    !     write (ndiag,*)
    
    !     do l = 0, nlsm
    !        write (ndiag,'(1x,a38,f16.3,f17.3)') soil(l),gast_i(l)*1.e-6,gast_o(l)*1.e-6
    !     end do
    
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write (ndiag,'(a)') 'Successfully made %sand and %clay'
       write (ndiag,*)
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
          maxindex = maxloc(wts_o(:)) 
          dynamicMaskList(no)%dstElement = real(maxindex(1)-1, kind=r4)
       end do
    end if

  end subroutine get_dominant_mapunit

end module mksoiltexMod
