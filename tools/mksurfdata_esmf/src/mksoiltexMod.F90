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
  use mkvarctl         , only : root_task, ndiag, spval, mpicom
  use mkvarctl         , only : unsetsoil
  use mkvarpar         , only : nlevsoi

  implicit none
  private ! By default make data private

#include <mpif.h>

  public :: mksoiltex      ! Set soil texture

  integer                :: mapunit_value_max
  integer                :: num_soil_textures
  type(ESMF_DynamicMask) :: dynamicMask

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  subroutine mksoiltex_check_layer_negative(data, no, lookup_index, name)
    !
    ! Arguments
    real(r4), intent(in) :: data
    integer, intent(in) :: no
    integer, intent(in) :: lookup_index
    character(*), intent(in) :: name
    if (data < 0._r4) then
       write(6,*)'ERROR: at no, lookup_index = ',no,lookup_index
       call shr_sys_abort('could not find a value >= 0 for '//name)
    end if
  end subroutine mksoiltex_check_layer_negative

  subroutine mksoiltex_i_to_o(no, lookup_index, n_scid, nlay, n_mapunits, name, val_neg_4, val_neg_other, data_i, data_o)
    !
    ! Fill output soil layers for a given texture component
    !
    ! Arguments
    integer, intent(in) :: no
    integer, intent(in) :: lookup_index
    integer, intent(in) :: n_scid, nlay, n_mapunits
    character(*), intent(in) :: name
    real(r4), intent(in) :: val_neg_4  ! Fallback value if read-in value is -4
    real(r4), intent(in) :: val_neg_other  ! Fallback value if read-in value is negative but not -4
    real(r4), intent(in) :: data_i(:,:,:)
    real(r4), intent(out) :: data_o(:,:)
    !
    ! Local variables
    integer :: l
    integer :: ier
    logical :: is_sand_dune

    ! Fill first layer of output array with first positive value on SCID dim of input array
    data_o(no,1) = data_i(1,1,lookup_index)
    if (data_o(no,1) < 0._r4) then
       do l = 2,n_scid
          if (data_i(1,l,lookup_index) >= 0._r4) then
             data_o(no,1) = data_i(1,l,lookup_index)
             exit
          end if
       end do
    end if

    ! Handle cases where no positive value was found
    if (data_o(no,1) < 0._r4) then
       is_sand_dune = int(data_o(no,1)) == -4
       if (is_sand_dune) then
          data_o(no,:) = val_neg_4
       else
          data_o(no,:) = val_neg_other
       end if
    end if

    ! Error on negative values in top layer
    call mksoiltex_check_layer_negative(data_o(no,1), no, lookup_index, name)

    ! Top soil layer is filled above. Here, we fill the other layers.
    do l = 2,nlay
       data_o(no,l) = data_i(l,1,lookup_index)

       ! If a layer is negative, fill it with the previous layer's value
       if (data_o(no,l) < 0._r4) then
          data_o(no,l) = data_o(no,l-1)
       end if
    end do

  end subroutine mksoiltex_i_to_o

  subroutine mksoiltex(file_mesh_i, file_mapunit_i, file_lookup_i, mesh_o, pioid_o, rc)
    !
    ! make %sand, %clay, organic carbon content, coarse fragments, bulk density,
    ! and pH measured in H2O
    !
    ! input/output variables
    character(len=*)  , intent(in)    :: file_mesh_i     ! input mesh/grid file name
    character(len=*)  , intent(in)    :: file_mapunit_i  ! input mapunit file name
    character(len=*)  , intent(in)    :: file_lookup_i   ! input data file name
    type(ESMF_Mesh)   , intent(in)    :: mesh_o          ! output mesh
    type(file_desc_t) , intent(inout) :: pioid_o
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
    integer                :: mapunit_value_max_local
    integer , allocatable  :: mask_i(:)
    real(r4), pointer      :: dataptr(:)
    integer                :: mapunit       ! temporary igbp soil mapunit
    integer,  allocatable  :: sand_i(:,:,:) ! input grid: percent sand
    integer,  allocatable  :: clay_i(:,:,:) ! input grid: percent clay
    integer,  allocatable  :: cfrag_i(:,:,:)  ! input grid: coarse fragments (vol% > 2 mm)
    real(r4), allocatable  :: bulk_i(:,:,:)  ! input grid: bulk density (g cm-3)
    real(r4), allocatable  :: orgc_i(:,:,:)  ! input grid: organic carbon content (gC kg-1)
    real(r4), allocatable  :: phaq_i(:,:,:)  ! input grid: soil pH measured in H2O (unitless)
    real(r4), allocatable  :: sand_o(:,:)   ! output grid: % sand
    real(r4), allocatable  :: clay_o(:,:)   ! output grid: % clay
    real(r4), allocatable  :: orgc_o(:,:)  ! output grid: organic carbon content (gC kg-1)
    real(r4), allocatable  :: cfrag_o(:,:)  ! output grid: coarse fragments (vol% > 2 mm)
    real(r4), allocatable  :: bulk_o(:,:)  ! output grid: bulk density (g cm-3)
    real(r4), allocatable  :: phaq_o(:,:)  ! output grid: soil pH measured in H2O (unitless)
    real(r4), allocatable  :: organic_o(:,:)  ! output grid: organic matter (kg m-3)
    integer                :: n_mapunits
    integer                :: lookup_index
    integer                :: SCID
    real(r4), allocatable  :: mapunit_i(:)  ! input grid: igbp soil mapunits
    integer , allocatable  :: mapunit_o(:)  ! output grid: igbp soil mapunits
    integer , allocatable  :: MapUnits(:)
    integer , allocatable  :: mapunit_lookup(:)
    type(var_desc_t)       :: pio_varid_sand
    type(var_desc_t)       :: pio_varid_clay
    type(var_desc_t)       :: pio_varid_orgc
    type(var_desc_t)       :: pio_varid_cfrag
    type(var_desc_t)       :: pio_varid_bulk
    type(var_desc_t)       :: pio_varid_phaq
    type(var_desc_t)       :: pio_varid_organic
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
       write(ndiag,'(a)') 'Attempting to make %sand, %clay, orgc, cfrag, bulk, phaq .....'
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
    allocate(orgc_o(ns_o,nlevsoi)) ; orgc_o(:,:) = spval
    allocate(cfrag_o(ns_o,nlevsoi)) ; cfrag_o(:,:) = spval
    allocate(bulk_o(ns_o,nlevsoi)) ; bulk_o(:,:) = spval
    allocate(phaq_o(ns_o,nlevsoi)) ; phaq_o(:,:) = spval
    allocate(organic_o(ns_o,nlevsoi)) ; organic_o(:,:) = spval

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
       else
          mask_i(ni) = 1
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
    !
    ! Note that dataptr (obtained from the input field) contains just a subset of the
    ! source data, based on the source data decomposition. So we need an mpi_allreduce to
    ! determine the global maximum value of mapunit.
    mapunit_value_max_local = maxval(dataptr)
    call mpi_allreduce(mapunit_value_max_local, mapunit_value_max, 1, MPI_INTEGER, MPI_MAX, mpicom, rcode)

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
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error calling mkfile_output for mapunits')
    call pio_syncfile(pioid_o)

    !---------------------------------
    ! Determine %sand, %clay, orgc, cfrag, bulk, phaq on output grid using
    ! mapunits
    !---------------------------------
    if (root_task) then
       write(ndiag,'(a)') 'WARNING: assigning       sand_o = -4 to 99%'
       write(ndiag,'(a)') 'WARNING: assigning other sand_o <  0 to 43%'
       write(ndiag,'(a)') 'WARNING: assigning       clay_o = -4 to 1%'
       write(ndiag,'(a)') 'WARNING: assigning other clay_o <  0 to 18%'
       write(ndiag,'(a)') 'WARNING: assigning       orgc_o = -4 to 1'
       write(ndiag,'(a)') 'WARNING: assigning other orgc_o <  0 to 0'
!      write(ndiag,'(a)') 'WARNING: same warnings for organic_o as for orgc_o'
       write(ndiag,'(a)') 'WARNING: same warnings for cfrag_o as for orgc_o'
       write(ndiag,'(a)') 'WARNING: assigning bulk_o < 0 to 1.5'
       write(ndiag,'(a)') 'WARNING: assigning phaq_o < 0 to 7'
    end if

    rcode = pio_openfile(pio_iosystem, pioid_i, pio_iotype, trim(file_lookup_i), pio_nowrite)

    rcode = pio_inq_dimid  (pioid_i, 'MapUnit', dimid)
    rcode = pio_inq_dimlen (pioid_i, dimid, n_mapunits)

    rcode = pio_inq_dimid  (pioid_i, 'soil_layer', dimid)
    rcode = pio_inq_dimlen (pioid_i, dimid, nlay)

    rcode = pio_inq_dimid  (pioid_i, 'SCID', dimid)
    rcode = pio_inq_dimlen (pioid_i, dimid, n_scid)

    ! Read In MapUnits from the input file
    allocate(MapUnits(n_mapunits), stat=ier)
    if (ier/=0) call shr_sys_abort()
    rcode = pio_inq_varid(pioid_i, 'MapUnit', pio_varid)
    rcode = pio_get_var(pioid_i, pio_varid, MapUnits)

    ! Determine the mapunit lookup index from the value of the MapUnit variable
    mapunit_value_max = maxval(MapUnits)
    allocate(mapunit_lookup(mapunit_value_max))
    mapunit_lookup(:) = -999
    do n = 1,size(MapUnits)
       mapunit_lookup(MapUnits(n)) = n
    end do

    allocate(sand_i(nlay,n_scid,n_mapunits), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(clay_i(nlay,n_scid,n_mapunits), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(orgc_i(nlay,n_scid,n_mapunits), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(cfrag_i(nlay,n_scid,n_mapunits), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(bulk_i(nlay,n_scid,n_mapunits), stat=ier)
    if (ier/=0) call shr_sys_abort()
    allocate(phaq_i(nlay,n_scid,n_mapunits), stat=ier)
    if (ier/=0) call shr_sys_abort()

    ! Get dimensions from input file and allocate memory for sand_i, clay_i,
    ! organic carbon content, coarse fragments, bulk density, pH measured in H2O
    rcode = pio_inq_varid(pioid_i, 'PCT_SAND', pio_varid_sand)
    rcode = pio_inq_varid(pioid_i, 'PCT_CLAY', pio_varid_clay)
    rcode = pio_inq_varid(pioid_i, 'ORGC', pio_varid_orgc)
    rcode = pio_inq_varid(pioid_i, 'CFRAG', pio_varid_cfrag)
    rcode = pio_inq_varid(pioid_i, 'BULK', pio_varid_bulk)
    rcode = pio_inq_varid(pioid_i, 'PHAQ', pio_varid_phaq)

    rcode = pio_get_var(pioid_i, pio_varid_sand, sand_i)
    rcode = pio_get_var(pioid_i, pio_varid_clay, clay_i)
    rcode = pio_get_var(pioid_i, pio_varid_orgc, orgc_i)
    rcode = pio_get_var(pioid_i, pio_varid_cfrag, cfrag_i)
    rcode = pio_get_var(pioid_i, pio_varid_bulk, bulk_i)
    rcode = pio_get_var(pioid_i, pio_varid_phaq, phaq_i)

    do no = 1,ns_o

       if (mapunit_o(no) == 0) then

          ! Set sand and clay to loam if mapunit is 0
          sand_o(no,:) = 43._r4
          clay_o(no,:) = 18._r4
          orgc_o(no,:) = 0._r4
          cfrag_o(no,:) = 0._r4
          bulk_o(no,:) = 1.5_r4  ! TODO Ok as a fill value?
          phaq_o(no,:) = 7._r4
          organic_o(no,:) = 0._r4

       else

          ! Determine lookup_index
          lookup_index = mapunit_lookup(mapunit_o(no))

          ! Fill output arrays
          call mksoiltex_i_to_o(no, lookup_index, n_scid, nlay, n_mapunits, &
               "sand", 99._r4, 43._r4, float(sand_i(:,:,:)), sand_o)
          call mksoiltex_i_to_o(no, lookup_index, n_scid, nlay, n_mapunits, &
               "clay", 1._r4, 18._r4, float(clay_i(:,:,:)), clay_o)
          call mksoiltex_i_to_o(no, lookup_index, n_scid, nlay, n_mapunits, &
               "cfrag", 1._r4, 0._r4, float(cfrag_i(:,:,:)), cfrag_o)
          call mksoiltex_i_to_o(no, lookup_index, n_scid, nlay, n_mapunits, &
               "bulk", 1.5_r4, 1.5_r4, bulk_i, bulk_o)  ! TODO: 1.5 ok for sand dunes and -7?
          call mksoiltex_i_to_o(no, lookup_index, n_scid, nlay, n_mapunits, &
               "phaq", 7._r4, 7._r4, phaq_i, phaq_o)
          call mksoiltex_i_to_o(no, lookup_index, n_scid, nlay, n_mapunits, &
               "orgc", 1._r4, 0._r4, orgc_i, orgc_o)

          ! ---------------------------------------------------------------
          ! Calculating organic_o
          ! ---------------------------------------------------------------
          ! Calculate organic from orgc_o, cfrag_o, and bulk_o, i.e. after
          ! these terms have been regridded. The plan is to move this
          ! calculation step from here to the CTSM. This approach keeps
          ! ORGC the same in fsurdat as in the raw data.
          !     Alternative approach previously available: Regrid organic_i
          ! (calculated from orgc_i, cfrag_i, and bulk_i) to organic_o. This
          ! approach first calculates organic_i and then regrids to organic_o
          ! rather than regridding all the terms first and then calculating
          ! organic_o. That would be organic_o_option 2, last available in
          ! commit d5f389a97.
          do l = 1, nlay
             organic_o(no,l) = orgc_o(no,l) * bulk_o(no,l) * &
                  (100._r4 - cfrag_o(no,l)) * 0.01_r4 / 0.58_r4
             ! Error on negative values in any layer
             call mksoiltex_check_layer_negative(organic_o(no,l), no, lookup_index, "organic")
          end do

       end if

    end do

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil percent sand"
    call mkfile_output(pioid_o,  mesh_o,  'PCT_SAND', sand_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error calling mkfile_output')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil percent clay"
    call mkfile_output(pioid_o,  mesh_o,  'PCT_CLAY', clay_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error calling mkfile_output')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil organic matter"
    call mkfile_output(pioid_o,  mesh_o,  'ORGANIC', organic_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error calling mkfile_output')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil organic carbon content"
    call mkfile_output(pioid_o,  mesh_o,  'ORGC', orgc_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error calling mkfile_output')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out coarse fragments in soil"
    call mkfile_output(pioid_o,  mesh_o,  'CFRAG', cfrag_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error calling mkfile_output')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil bulk density"
    call mkfile_output(pioid_o,  mesh_o,  'BULK', bulk_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error calling mkfile_output')

    if (root_task)  write(ndiag, '(a)') trim(subname)//" writing out soil pH measured in H2O"
    call mkfile_output(pioid_o,  mesh_o,  'PHAQ', phaq_o, lev1name='nlevsoi', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) call shr_sys_abort('error calling mkfile_output')

    call pio_syncfile(pioid_o)

    ! Release memory
    call ESMF_RouteHandleDestroy(routehandle, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_MeshDestroy(mesh_i, nogarbage = .true., rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) call shr_sys_abort()
    call ESMF_VMLogMemInfo("After destroy operations in "//trim(subname))

    if (root_task) then
       write(ndiag,'(a)') 'Successfully made %sand, %clay, orgc, cfrag, bulk, phaq .....'
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
