module lnd_set_decomp_and_domain

  use ESMF         , only : ESMF_VM, ESMF_MESH, ESMF_DistGrid, ESMF_Field
  use ESMF         , only : ESMF_RouteHandle, ESMF_SUCCESS, ESMF_MeshCreate
  use ESMF         , only : ESMF_FileFormat_ESMFMESH, ESMF_VMLogMemInfo
  use ESMF         , only : ESMF_FieldCreate, ESMF_FieldRedistStore, ESMF_FieldRedist
  use ESMF         , only : ESMF_FieldGet, ESMF_GridCreateNoPeriDimUfrm
  use ESMF         , only : ESMF_FieldRegrid, ESMF_FieldRegridStore
  use ESMF         , only : ESMF_Grid, ESMF_ARRAY, ESMF_DistGridCreate, ESMF_TYPEKIND_R8
  use ESMF         , only : ESMF_MESHLOC_ELEMENT, ESMF_FieldCreate, ESMF_STAGGERLOC_CORNER, ESMF_STAGGERLOC_CENTER
  use ESMF         , only : ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_DSTAREA
  use ESMF         , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_TERMORDER_SRCSEQ
  use ESMF         , only : ESMF_REGION_TOTAL, ESMF_REDUCE_SUM, ESMF_ARRAYCREATE
  use ESMF         , only : ESMF_MeshGet, ESMF_DistGridGet, ESMF_LOGFOUNDERROR
  use ESMF         , only : ESMF_LOGERR_PASSTHRU, ESMF_VMAllReduce, ESMF_FieldRegridGetArea
  use ESMF         , only : ESMF_FieldDestroy
  use shr_kind_mod , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod  , only : shr_sys_abort
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use spmdMod      , only : masterproc, mpicom
  use clm_varctl   , only : iulog, inst_suffix, FL => fname_len
  use abortutils   , only : endrun

  implicit none
  private ! except

  ! Module public routines
  public :: lnd_set_decomp_and_domain_from_readmesh
  public :: lnd_set_mesh_for_single_column
  public :: lnd_set_decomp_and_domain_for_single_column

  ! Module private routines
  private :: lnd_get_global_dims
  private :: lnd_set_lndmask_from_maskmesh
  private :: lnd_set_lndmask_from_lndmesh
  private :: lnd_set_lndmask_from_fatmlndfrc
  private :: lnd_set_ldomain_gridinfo_from_mesh
  private :: chkerr
  private :: pio_check_err

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lnd_set_decomp_and_domain_from_readmesh(driver, vm, meshfile_lnd, meshfile_mask, mesh_ctsm, &
       ni, nj, rc)

    use decompInitMod , only : decompInit_lnd
    use domainMod     , only : ldomain, domain_init
    use decompMod     , only : gindex_global, bounds_type, get_proc_bounds
    use clm_varpar    , only : nlevsoi
    use clm_varctl    , only : use_soil_moisture_streams

    ! input/output variables
    character(len=*)    , intent(in)    :: driver ! cmeps or lilac
    type(ESMF_VM)       , intent(in)    :: vm
    character(len=*)    , intent(in)    :: meshfile_lnd
    character(len=*)    , intent(in)    :: meshfile_mask
    type(ESMF_Mesh)     , intent(out)   :: mesh_ctsm
    integer             , intent(out)   :: ni,nj  ! global grid dimensions
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Mesh)        :: mesh_maskinput
    type(ESMF_Mesh)        :: mesh_lndinput
    type(ESMF_DistGrid)    :: distgrid_ctsm
    type(ESMF_Field)       :: field_lnd
    type(ESMF_Field)       :: field_ctsm
    type(ESMF_RouteHandle) :: rhandle_lnd2ctsm
    integer                :: g,n             ! indices
    integer                :: nlnd, nocn      ! local size of arrays
    integer                :: gsize           ! global size of grid
    logical                :: isgrid2d        ! true => grid is 2d
    type(bounds_type)      :: bounds          ! bounds
    integer                :: begg,endg       ! local bounds
    integer  , pointer     :: gindex_lnd(:)   ! global index space for just land points
    integer  , pointer     :: gindex_ocn(:)   ! global index space for just ocean points
    integer  , pointer     :: gindex_ctsm(:)  ! global index space for land and ocean points
    integer  , pointer     :: lndmask_glob(:)
    real(r8) , pointer     :: lndfrac_glob(:)
    real(r8) , pointer     :: lndfrac_loc_input(:)
    real(r8) , pointer     :: dataptr1d(:)
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Write diag info
    if (masterproc) then
       write(iulog,*)
       write(iulog,'(a)')' Input land mesh file '//trim(meshfile_lnd)
       write(iulog,'(a)')' Input mask mesh file '//trim(meshfile_mask)
       if (trim(meshfile_mask) /= trim(meshfile_lnd)) then
          write(iulog, '(a)') ' Obtaining land mask and fraction from mask file '//trim(meshfile_mask)
       else
          write(iulog, '(a)') ' Obtaining land mask and fraction from land mesh file '//trim(meshfile_lnd)
       end if
       write(iulog,*)
    end if

    ! Determine global 2d sizes from read of dimensions of surface dataset and allocate global memory
    call lnd_get_global_dims(ni, nj, gsize, isgrid2d)

    ! Read in the land mesh from the file
    mesh_lndinput = ESMF_MeshCreate(filename=trim(meshfile_lnd), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(driver) == 'cmeps') then
       ! Read in mask meshfile if needed
       if (trim(meshfile_mask) /= trim(meshfile_lnd)) then
#ifdef DEBUG
          ! This will get added to the ESMF PET files if DEBUG=TRUE and CREATE_ESMF_PET_FILES=TRUE
          call ESMF_VMLogMemInfo("clm: Before lnd mesh create in ")
#endif
          mesh_maskinput = ESMF_MeshCreate(filename=trim(meshfile_mask), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
#ifdef DEBUG
          ! This will get added to the ESMF PET files if DEBUG=TRUE and CREATE_ESMF_PET_FILES=TRUE
          call ESMF_VMLogMemInfo("clm: After lnd mesh create in ")
#endif
          ! Determine lndmask_glob and lndfrac_glob
          ! obain land mask and land fraction by mapping ocean mesh conservatively to land mesh
          ! Note that lndmask_glob and lndfrac_loc_input are allocated in lnd_set_lndmask_from_maskmesh
          call lnd_set_lndmask_from_maskmesh(mesh_lndinput, mesh_maskinput, vm, gsize, lndmask_glob, &
               lndfrac_loc_input, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
#ifdef DEBUG
          ! This will get added to the ESMF PET files if DEBUG=TRUE and CREATE_ESMF_PET_FILES=TRUE
          call ESMF_VMLogMemInfo("clm: After lnd_set_lndmask_from_maskmesh ")
#endif
       else
          ! obtain land mask from land mesh file - assume that land frac is identical to land mask
          call lnd_set_lndmask_from_lndmesh(mesh_lndinput, vm, gsize, lndmask_glob, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    else if (trim(driver) == 'lilac') then
       call lnd_set_lndmask_from_fatmlndfrc(lndmask_glob, lndfrac_glob, ni,nj)
    else
       call shr_sys_abort('driver '//trim(driver)//' is not supported, must be lilac or cmeps')
    end if

    ! Determine lnd decomposition that will be used by ctsm from lndmask_glob
    call decompInit_lnd(lni=ni, lnj=nj, amask=lndmask_glob)

    ! Determine ocn decomposition that will be used to create the full mesh
    ! note that the memory for gindex_ocn will be allocated in the following call
    ! but deallocated at the end of this routine
    call decompInit_ocn(ni=ni, nj=nj, amask=lndmask_glob, gindex_ocn=gindex_ocn)

    ! Get JUST gridcell processor bounds
    ! Remaining bounds (landunits, columns, patches) will be set after calling decompInit_glcp
    ! so get_proc_bounds is called twice and the gridcell information is just filled in twice
    call get_proc_bounds(bounds)
    begg = bounds%begg
    endg = bounds%endg

    ! Create ctsm gindex_lnd
    nlnd = endg - begg + 1
    allocate(gindex_lnd(nlnd))
    do g = begg, endg
       n = 1 + (g - begg)
       gindex_lnd(n) = gindex_global(g-begg+1)
    end do

    ! Initialize domain data structure
    call domain_init(domain=ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine ldomain%mask
    do g = begg, endg
       n = 1 + (g - begg)
       ldomain%mask(g) = lndmask_glob(gindex_lnd(n))
    end do

    ! Deallocate global pointer memory
    deallocate(lndmask_glob)

    ! Generate a ctsm global index that includes both land and ocean points
    nocn = size(gindex_ocn)
    allocate(gindex_ctsm(nlnd + nocn))
    do n = 1,nlnd+nocn
       if (n <= nlnd) then
          gindex_ctsm(n) = gindex_lnd(n)
       else
          gindex_ctsm(n) = gindex_ocn(n-nlnd)
       end if
    end do

    ! Generate a new mesh on the gindex decomposition
    distGrid_ctsm = ESMF_DistGridCreate(arbSeqIndexList=gindex_ctsm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    mesh_ctsm = ESMF_MeshCreate(mesh_lndinput, elementDistGrid=distgrid_ctsm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set ldomain%lonc, ldomain%latc and ldomain%area
    call lnd_set_ldomain_gridinfo_from_mesh(mesh_ctsm, vm, gindex_ctsm, begg, endg, isgrid2d, ni, nj, ldomain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set ldomain%lfrac
    ! Create fields on the input decomp and ctsm decomp
    ! Determine route handle to do a redist from input mesh read decomp to ctsm decomp
    ! Fill in the field on the input mesh read decomp with lndfrac_loc_input values
    ! Redistribute field_lnd to field_ctsm

    ! Determine ldomain%frac using ctsm decomposition
    if (trim(driver) == 'cmeps') then

       if (trim(meshfile_mask) /= trim(meshfile_lnd)) then
          field_lnd = ESMF_FieldCreate(mesh_lndinput, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          field_ctsm = ESMF_FieldCreate(mesh_ctsm, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldRedistStore(field_lnd, field_ctsm, routehandle=rhandle_lnd2ctsm, &
               ignoreUnmatchedIndices=.true., rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field_lnd, farrayptr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do n = 1,size(dataptr1d)
             dataptr1d(n) = lndfrac_loc_input(n)
          end do
          call ESMF_FieldRedist(field_lnd, field_ctsm, routehandle=rhandle_lnd2ctsm, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field_ctsm, farrayptr=dataptr1d, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          do g = begg, endg
             n = 1 + (g - begg)
             ldomain%frac(g) = dataptr1d(n)
          end do
       else
          ! ASSUME that land fraction is identical to land mask in this case
          do g = begg, endg
             ldomain%frac(g) = ldomain%mask(g)
          end do
       end if

    else

       do g = begg, endg
          n = 1 + (g - begg)
          ldomain%frac(g) = lndfrac_glob(gindex_lnd(n))
       end do
       deallocate(lndfrac_glob)

    end if

    ! Deallocate local pointer memory
    deallocate(gindex_lnd)
    deallocate(gindex_ocn)
    deallocate(gindex_ctsm)

  end subroutine lnd_set_decomp_and_domain_from_readmesh

  !===============================================================================
  subroutine lnd_set_mesh_for_single_column(scol_lon, scol_lat, mesh, rc)

    ! Generate a mesh for single column
    use clm_varcon, only : spval

    ! input/output variables
    real(r8)        , intent(in)  :: scol_lon
    real(r8)        , intent(in)  :: scol_lat
    type(ESMF_Mesh) , intent(out) :: mesh
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_Grid)        :: lgrid
    integer                :: maxIndex(2)
    real(r8)               :: mincornerCoord(2)
    real(r8)               :: maxcornerCoord(2)
    character(len=*), parameter :: subname= ' (lnd_set_mesh_for_single_column) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Use center and come up with arbitrary area delta lon and lat = .1 degree
    maxIndex(1)       = 1                ! number of lons
    maxIndex(2)       = 1                ! number of lats
    mincornerCoord(1) = scol_lon - .1_r8 ! min lon
    mincornerCoord(2) = scol_lat - .1_r8 ! min lat
    maxcornerCoord(1) = scol_lon + .1_r8 ! max lon
    maxcornerCoord(2) = scol_lat + .1_r8 ! max lat
    ! create the ESMF grid
    lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
         mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
         staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create the mesh from the lgrid
    mesh = ESMF_MeshCreate(lgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine lnd_set_mesh_for_single_column

  !===============================================================================
  subroutine lnd_set_decomp_and_domain_for_single_column(scol_lon, scol_lat, scol_mask, scol_frac)

    use decompInitMod , only : decompInit_lnd
    use decompMod     , only : bounds_type, get_proc_bounds
    use domainMod     , only : ldomain, domain_init
    use clm_varctl    , only : use_soil_moisture_streams
    use clm_varpar    , only : nlevsoi
    use clm_varcon    , only : spval

    ! input/output variables
    real(r8) , intent(in) :: scol_lon
    real(r8) , intent(in) :: scol_lat
    integer  , intent(in) :: scol_mask
    real(r8) , intent(in) :: scol_frac

    ! local variables
    type(bounds_type) :: bounds                ! bounds
    !-------------------------------------------------------------------------------

    ! Determine decomp and ldomain
    call decompInit_lnd(lni=1, lnj=1, amask=(/1/))

    ! Initialize processor bounds
    call get_proc_bounds(bounds)

    ! Initialize domain data structure
    call domain_init(domain=ldomain, isgrid2d=.false., ni=1, nj=1, nbeg=1, nend=1)

    ! Initialize ldomain attributes
    ldomain%lonc(1) = scol_lon
    ldomain%latc(1) = scol_lat
    ldomain%area(1) = spval
    ldomain%mask(1) = scol_mask
    ldomain%frac(1) = scol_frac

  end subroutine lnd_set_decomp_and_domain_for_single_column

  !===============================================================================
  subroutine lnd_get_global_dims(ni, nj, gsize, isgrid2d)

    ! Determine global 2d sizes from read of dimensions of surface dataset
    !
    ! Meshes do not indicate if the mesh can be represented as a logically rectangular
    ! grid. However, CTSM needs this information in the history file generation via the
    ! logical variable isgrid2d. Since for CMEPS and LILAC there is no longer the need for
    ! the fatmlndfrc file (where the isgrid2d variable was determined from before), the
    ! surface dataset is now used to determine if the underlying grid is 2d or not.

    use clm_varctl  , only : fsurdat, single_column
    use fileutils   , only : getfil
    use ncdio_pio   , only : ncd_io, file_desc_t, ncd_pio_openfile, ncd_pio_closefile, ncd_inqdlen, ncd_inqdid

    ! input/output variables
    integer, intent(out) :: ni
    integer, intent(out) :: nj
    integer, intent(out) :: gsize
    logical, intent(out) :: isgrid2d

    ! local variables
    character(len=CL) :: locfn
    type(file_desc_t) :: ncid    ! netcdf file id
    integer           :: dimid   ! netCDF dimension id
    logical           :: readvar ! read variable in or not
    logical           :: dim_exists
    logical           :: dim_found = .false.
    !-------------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read global dimensions from surface dataset'
       if (fsurdat == ' ') then
          write(iulog,*)'fsurdat must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    endif
    call getfil(fsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    dim_found = .false.
    call ncd_inqdid(ncid, 'lsmlon', dimid, dim_exists)
    if ( dim_exists ) then
       dim_found = .true.
       call ncd_inqdlen(ncid, dimid, ni, 'lsmlon')
       call ncd_inqdlen(ncid, dimid, nj, 'lsmlat')
    end if
    if (.not. dim_found) then
       call ncd_inqdid(ncid, 'gridcell', dimid, dim_exists)
       if ( dim_exists ) then
          dim_found = .true.
          call ncd_inqdlen(ncid, dimid, ni, 'gridcell')
          nj = 1
       end if
    end if
    if (.not. dim_found) then
       call shr_sys_abort('ERROR: surface dataset does not contain dims of lsmlon,lsmlat or gridcell')
    end if
    call ncd_pio_closefile(ncid)
    gsize = ni*nj
    if (single_column) then
       isgrid2d = .true.
    else if (nj == 1) then
       isgrid2d = .false.
    else
       isgrid2d = .true.
    end if
    if (masterproc) then
       write(iulog,'(a,2(i8,2x))') 'global ni,nj = ',ni,nj
       if (isgrid2d) then
          write(iulog,'(a)') 'model grid is 2-dimensional'
       else
          write(iulog,'(a)') 'model grid is not 2-dimensional'
       end if
    end if

  end subroutine lnd_get_global_dims

  !===============================================================================
  subroutine lnd_set_lndmask_from_maskmesh(mesh_lnd, mesh_mask, vm, gsize, lndmask_glob, lndfrac_loc, rc)

    ! If the landfrac/landmask file does not exists then determine the
    ! land fraction and land mask on the land grid by mapping the mask
    ! from the mesh_mask to the mesh_lnd mesh. Then write out the
    ! landfrac/landmesh file.  If the landfrac/landmask file does
    ! exist then simply read in the global land fraction and land mask
    ! from the file

    ! input/out variables
    type(ESMF_Mesh)     , intent(in)  :: mesh_lnd
    type(ESMF_Mesh)     , intent(in)  :: mesh_mask
    type(ESMF_VM)       , intent(in)  :: vm
    integer             , intent(in)  :: gsize
    integer             , pointer     :: lndmask_glob(:)
    real(r8)            , pointer     :: lndfrac_loc(:)
    integer             , intent(out) :: rc

    ! local variables:
    type(ESMF_DistGrid)    :: distgrid_lnd
    type(ESMF_RouteHandle) :: rhandle_mask2lnd
    type(ESMF_Field)       :: field_lnd
    type(ESMF_Field)       :: field_mask
    type(ESMF_DistGrid)    :: distgrid_mask
    integer  , pointer     :: gindex_input(:) ! global index space for land and ocean points
    integer  , pointer     :: lndmask_loc(:)
    integer  , pointer     :: itemp_glob(:)
    real(r8) , pointer     :: maskmask_loc(:) ! on ocean mesh
    real(r8) , pointer     :: maskfrac_loc(:) ! on land mesh
    real(r8) , pointer     :: dataptr1d(:)
    type(ESMF_Array)       :: elemMaskArray
    integer                :: lsize_lnd
    integer                :: lsize_mask
    integer                :: n, spatialDim
    integer                :: srcMaskValue = 0
    integer                :: dstMaskValue = -987987 ! spval for RH mask values
    integer                :: srcTermProcessing_Value = 0
    integer                :: klen
    real(r8)               :: fminval = 0.001_r8
    real(r8)               :: fmaxval = 1._r8
    logical                :: lexist
    logical                :: checkflag = .false.
    character(len=FL)      :: flandfrac
    character(len=CL)      :: flandfrac_status
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    flandfrac = './init_generated_files/ctsm_landfrac'//trim(inst_suffix)//'.nc'
    klen = len_trim(flandfrac) - 3 ! remove the .nc
    flandfrac_status = flandfrac(1:klen)//'.status'

    allocate(lndmask_glob(gsize)); lndmask_glob(:) = 0

    ! Determine if lndfrac/lndmask file exists
    inquire(file=trim(flandfrac), exist=lexist)

    if (lexist) then

       ! If file exists - read in lndmask and lndfrac
       if (masterproc) then
          write(iulog,*)
          write(iulog,'(a)')' Reading in land fraction and land mask from '//trim(flandfrac)
       end if
    else

       ! If file does not exist - compute lndmask and lndfrac and write to output file
       if (masterproc) then
          write(iulog,*)
          write(iulog,'(a)')' Computing land fraction and land mask by mapping mask from mesh_mask file'
       end if
       call ESMF_MeshGet(mesh_lnd, spatialDim=spatialDim, numOwnedElements=lsize_lnd, &
            elementDistGrid=distgrid_lnd, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(lndmask_loc(lsize_lnd))
       allocate(lndfrac_loc(lsize_lnd))

       ! create fields on land and ocean meshes
       field_lnd = ESMF_FieldCreate(mesh_lnd, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       field_mask = ESMF_FieldCreate(mesh_mask, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! create route handle to map ocean mask from mask mesh to land mesh
       call ESMF_FieldRegridStore(field_mask, field_lnd, routehandle=rhandle_mask2lnd, &
            srcMaskValues=(/srcMaskValue/), dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_DSTAREA, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! fill in values for field_mask with mask on mask mesh
       call ESMF_MeshGet(mesh_mask, elementdistGrid=distgrid_mask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_DistGridGet(distgrid_mask, localDe=0, elementCount=lsize_mask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(maskmask_loc(lsize_mask))
       elemMaskArray = ESMF_ArrayCreate(distgrid_mask, maskmask_loc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshGet(mesh_mask, elemMaskArray=elemMaskArray, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldGet(field_mask, farrayptr=dataptr1d, rc=rc)
       dataptr1d(:) = maskmask_loc(:)

       ! map mask mask to land mesh
       call ESMF_FieldRegrid(field_mask, field_lnd, routehandle=rhandle_mask2lnd, &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(mesh_lnd, spatialDim=spatialDim, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(maskfrac_loc(lsize_lnd))
       call ESMF_FieldGet(field_lnd, farrayptr=maskfrac_loc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,lsize_lnd
          lndfrac_loc(n) = 1._r8 - maskfrac_loc(n)
          if (lndfrac_loc(n) > fmaxval) lndfrac_loc(n) = 1._r8
          if (lndfrac_loc(n) < fminval) lndfrac_loc(n) = 0._r8
          if (lndfrac_loc(n) /= 0._r8) then
             lndmask_loc(n) = 1
          else
             lndmask_loc(n) = 0
          end if
       enddo
       call ESMF_FieldDestroy(field_lnd)
       call ESMF_FieldDestroy(field_mask)

       ! determine global landmask_glob - needed to determine the ctsm decomposition
       ! land frac, lats, lons and areas will be done below
       allocate(gindex_input(lsize_lnd))
       call ESMF_DistGridGet(distgrid_lnd, 0, seqIndexList=gindex_input, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,lsize_lnd
          lndmask_glob(gindex_input(n)) = lndmask_loc(n)
       end do
       allocate(itemp_glob(gsize))
       call ESMF_VMAllReduce(vm, sendData=lndmask_glob, recvData=itemp_glob, count=gsize, &
            reduceflag=ESMF_REDUCE_SUM, rc=rc)
       lndmask_glob(:) = int(itemp_glob(:))
       deallocate(itemp_glob)

       ! deallocate memory
       deallocate(maskmask_loc)
       deallocate(lndmask_loc)

    end if

  end subroutine lnd_set_lndmask_from_maskmesh

  !===============================================================================
  subroutine lnd_set_lndmask_from_lndmesh(mesh_lnd, vm, gsize, lndmask_glob, rc)

    ! input/out variables
    type(ESMF_Mesh)     , intent(in)  :: mesh_lnd
    type(ESMF_VM)       , intent(in)  :: vm
    integer             , intent(in)  :: gsize
    integer             , pointer     :: lndmask_glob(:)
    integer             , intent(out) :: rc

    ! local variables:
    integer             :: n
    integer             :: lsize
    integer , pointer   :: gindex(:)
    integer , pointer   :: lndmask_loc(:)
    integer , pointer   :: itemp_glob(:)
    type(ESMF_DistGrid) :: distgrid
    type(ESMF_Array)    :: elemMaskArray
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine lsize and distgrid_lnd
    call ESMF_MeshGet(mesh_lnd, elementdistGrid=distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distgrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine lndmask_loc
    ! The call to ESMF_MeshGet fills in the values of lndmask_loc
    allocate(lndmask_loc(lsize))
    elemMaskArray = ESMF_ArrayCreate(distgrid, lndmask_loc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_lnd, elemMaskArray=elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine global landmask_glob - needed to determine the ctsm decomposition
    ! land frac, lats, lons and areas will be done below
    allocate(gindex(lsize))
    allocate(itemp_glob(gsize))
    call ESMF_DistGridGet(distgrid, 0, seqIndexList=gindex, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(lndmask_glob(gsize)); lndmask_glob(:) = 0

    do n = 1,lsize
       lndmask_glob(gindex(n)) = lndmask_loc(n)
    end do
    call ESMF_VMAllReduce(vm, sendData=lndmask_glob, recvData=itemp_glob, count=gsize, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    lndmask_glob(:) = int(itemp_glob(:))
    deallocate(itemp_glob)
    deallocate(gindex)
    deallocate(lndmask_loc)

  end subroutine lnd_set_lndmask_from_lndmesh

  !===============================================================================
  subroutine lnd_set_lndmask_from_fatmlndfrc(mask, frac, ni, nj)

    ! Read the surface dataset grid related information
    ! This is used to set the domain decomposition - so global data is read here

    use clm_varctl , only : fatmlndfrc
    use fileutils  , only : getfil
    use ncdio_pio  , only : ncd_io, ncd_pio_openfile, ncd_pio_closefile, ncd_inqfdims, file_desc_t

    ! input/output variables
    integer         , pointer       :: mask(:)   ! grid mask
    real(r8)        , pointer       :: frac(:)   ! grid fraction
    integer         , intent(out)   :: ni, nj    ! global grid sizes

    ! local variables
    logical               :: isgrid2d
    integer               :: dimid,varid ! netCDF id's
    integer               :: ns          ! size of grid on file
    integer               :: n,i,j       ! index
    integer               :: ier         ! error status
    type(file_desc_t)     :: ncid        ! netcdf id
    character(len=256)    :: varname     ! variable name
    character(len=256)    :: locfn       ! local file name
    logical               :: readvar     ! read variable in or not
    integer , allocatable :: idata2d(:,:)
    real(r8), allocatable :: rdata2d(:,:)
    integer               :: unitn
    character(len=32)     :: subname = 'lnd_set_mask_from_fatmlndfrc' ! subroutine name
    !-----------------------------------------------------------------------

    ! Open file
    call getfil( fatmlndfrc, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Determine dimensions and if grid file is 2d or 1d
    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    if (masterproc) then
       write(iulog,*)'lat/lon grid flag (isgrid2d) is ',isgrid2d
    end if
    allocate(frac(ni*nj), mask(ni*nj))

    if (isgrid2d) then
       ! Grid is 2d
       allocate(idata2d(ni,nj))
       idata2d(:,:) = 1
       call ncd_io(ncid=ncid, varname='mask', data=idata2d, flag='read', readvar=readvar)
       if (readvar) then
          do j = 1,nj
             do i = 1,ni
                n = (j-1)*ni + i
                mask(n) = idata2d(i,j)
             enddo
          enddo
       else
          call endrun( msg=' ERROR: mask not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))
       end if
       deallocate(idata2d)
       allocate(rdata2d(ni,nj))
       rdata2d(:,:) = 1._r8
       call ncd_io(ncid=ncid, varname='frac', data=rdata2d, flag='read', readvar=readvar)
       if (readvar) then
          do j = 1,nj
             do i = 1,ni
                n = (j-1)*ni + i
                frac(n) = rdata2d(i,j)
             enddo
          enddo
       else
          call endrun( msg=' ERROR: mask not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))
       end if
       deallocate(rdata2d)
    else
       ! Grid is not 2d
       call ncd_io(ncid=ncid, varname='mask', data=mask, flag='read', readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: mask not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))
       end if
       call ncd_io(ncid=ncid, varname='frac', data=frac, flag='read', readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: frac not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))
       end if
    end if

    ! Close file
    call ncd_pio_closefile(ncid)

  end subroutine lnd_set_lndmask_from_fatmlndfrc

  !===============================================================================
  subroutine lnd_set_ldomain_gridinfo_from_mesh(mesh, vm, gindex, begg, endg, isgrid2d, ni, nj, ldomain, rc)

    use domainMod  , only : domain_type, lon1d, lat1d
    use clm_varcon , only : re

    use clm_varcon , only : grlnd
    use fileutils  , only : getfil
    use ncdio_pio  , only : ncd_io, file_desc_t, ncd_pio_openfile, ncd_pio_closefile

    ! input/output variables
    type(ESMF_Mesh)   , intent(in)    :: mesh
    type(ESMF_VM)     , intent(in)    :: vm
    integer           , intent(in)    :: gindex(:)
    integer           , intent(in)    :: begg,endg
    logical           , intent(in)    :: isgrid2d
    integer           , intent(in)    :: ni, nj
    type(domain_type) , intent(inout) :: ldomain
    integer           , intent(out)   :: rc

    ! local variables
    integer            :: g,n
    integer            :: gsize
    integer            :: numownedelements
    real(r8) , pointer :: ownedElemCoords(:)
    integer            :: spatialDim
    real(r8) , pointer :: dataptr1d(:)
    real(r8) , pointer :: lndlats_glob(:)
    real(r8) , pointer :: lndlons_glob(:)
    real(r8) , pointer :: rtemp_glob(:)
    type(ESMF_Field)   :: areaField
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine ldoman%latc and ldomain%lonc
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg,endg
       n = g - begg + 1
       ldomain%lonc(g) = ownedElemCoords(2*n-1)
       if (ldomain%lonc(g) == 360._r8) ldomain%lonc(g) = 0._r8 ! TODO: why the difference?
       ldomain%latc(g) = ownedElemCoords(2*n)
    end do

    ! Create ldomain%area by querying the mesh on the ctsm decomposition
    areaField = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(areaField, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(areaField, farrayPtr=dataptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = begg, endg
       ldomain%area(g) = dataptr1d(g-begg+1) * (re*re)
    end do
    call ESMF_FieldDestroy(areaField)

    ! If grid is 2d, determine lon1d and lat1d from mesh
    if (isgrid2d) then
       gsize = ni*nj
       allocate(rtemp_glob(gsize))

       ! Determine lon1d
       allocate(lndlons_glob(gsize))
       lndlons_glob(:) = 0._r8
       do n = 1,numownedelements
          if (ownedElemCoords(2*n-1) == 360._r8) then ! TODO: why is this needed?
             lndlons_glob(gindex(n)) = 0._r8
          else
             lndlons_glob(gindex(n)) = ownedElemCoords(2*n-1)
          end if
       end do
       call ESMF_VMAllReduce(vm, sendData=lndlons_glob, recvData=rtemp_glob, count=gsize, &
            reduceflag=ESMF_REDUCE_SUM, rc=rc)
       deallocate(lndlons_glob)
       allocate(lon1d(ni))
       do n = 1,ni
          lon1d(n) = rtemp_glob(n)
       end do

       ! Determine lat1d
       allocate(lndlats_glob(gsize))
       lndlats_glob(:) = 0._r8
       do n = 1,numownedelements
          lndlats_glob(gindex(n)) = ownedElemCoords(2*n)
       end do
       call ESMF_VMAllReduce(vm, sendData=lndlats_glob, recvData=rtemp_glob, count=gsize, &
            reduceflag=ESMF_REDUCE_SUM, rc=rc)
       deallocate(lndlats_glob)
       allocate(lat1d(nj))
       do n = 1,nj
          lat1d(n) = rtemp_glob((n-1)*ni + 1)
       end do
       deallocate(rtemp_glob)
    end if

  end subroutine lnd_set_ldomain_gridinfo_from_mesh

  !===============================================================================
  subroutine pio_check_err(ierror, description)
    use pio, only : PIO_NOERR
    integer     , intent(in) :: ierror
    character(*), intent(in) :: description
    if (ierror /= PIO_NOERR) then
       write (*,'(6a)') 'ERROR ', trim(description)
       call shr_sys_abort()
    endif
  end subroutine pio_check_err

  !===============================================================================
  logical function chkerr(rc, line, file)
    integer          , intent(in) :: rc
    integer          , intent(in) :: line
    character(len=*) , intent(in) :: file

    integer :: lrc
    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

  !===============================================================================
  subroutine lnd_set_read_write_landmask(flandfrac, flandfrac_status, write_file, read_file, &
       lndmask_glob, lndfrac_glob, gsize)

    ! Write or read landfrac and landmask to file so that mapping does not have to be done each time
    ! mapping the ocean mask to the land grid, it's possible that landfrac will be

    use ncdio_pio , only : ncd_io, file_desc_t, ncd_pio_openfile, ncd_pio_closefile
    use ncdio_pio , only : ncd_defdim, ncd_defvar, ncd_enddef, ncd_inqdlen
    use ncdio_pio , only : ncd_int, ncd_double, ncd_pio_createfile

    ! input/output variables
    character(len=*) , intent(in) :: flandfrac
    character(len=*) , intent(in) :: flandfrac_status
    logical          , intent(in) :: write_file
    logical          , intent(in) :: read_file
    integer          , pointer    :: lndmask_glob(:)
    real(r8)         , pointer    :: lndfrac_glob(:)
    integer          , intent(in) :: gsize

    ! local variables
    type(file_desc_t) :: pioid ! netcdf file id
    integer           :: dimid
    integer           :: iun
    integer           :: ioe
    integer           :: ier
    logical           :: lexists
    !-------------------------------------------------------------------------------

    if (write_file) then

       if (masterproc) then
          write(iulog,*)
          write(iulog,'(a)') 'lnd_set_decomp_and_domain: writing landmask and landfrac data to landfrac.nc'
          write(iulog,*)

          ! Remove file if it exists
          inquire(file=trim(flandfrac), exist=lexists)
          if (lexists) then
             open(unit=9876, file=flandfrac, status='old', iostat=ioe)
             if (ioe == 0) then
                close(9876, status='delete')
             end if
          end if
       end if
       call mpi_barrier(mpicom,ier)

       call ncd_pio_createfile(pioid, trim(flandfrac))
       call ncd_defdim (pioid, 'gridcell', gsize, dimid)
       call ncd_defvar(ncid=pioid, varname='landmask', xtype=ncd_int   , dim1name='gridcell')
       call ncd_defvar(ncid=pioid, varname='landfrac', xtype=ncd_double, dim1name='gridcell')
       call ncd_enddef(pioid)
       call ncd_io(ncid=pioid, varname='landmask', data=lndmask_glob, flag='write')
       call ncd_io(ncid=pioid, varname='landfrac', data=lndfrac_glob, flag='write')
       call ncd_pio_closefile(pioid)

       call mpi_barrier(mpicom,ier)
       if (masterproc) then
          open (newunit=iun, file=flandfrac_status, status='unknown',  iostat=ioe)
          if (ioe /= 0) then
             call endrun(msg='ERROR failed to open file '//trim(flandfrac_status)//errMsg(sourcefile, __LINE__))
          end if
          write(iun,'(a)')'Successfully wrote out '//trim(flandfrac_status)
          close(iun)
          write(iulog,'(a)')' Successfully wrote land fraction/mask status file '//trim(flandfrac_status)
       end if

    else if (read_file) then

       if (masterproc) then
          write(iulog,*)
          write(iulog,'(a)') 'lnd_set_decomp_and_domain: reading landmask and landfrac data from landfrac.nc'
          write(iulog,*)
       end if
       inquire(file=trim(flandfrac_status), exist=lexists)
       if (.not. lexists) then
          ! To read the file first check that the status file exists
          if (masterproc) then
             write(iulog,'(a)')' failed to find file '//trim(flandfrac_status)
             write(iulog,'(a)')' this indicates a problem in creating '//trim(flandfrac_status)
             write(iulog,'(a)')' remove '//trim(flandfrac)//' and try again'
          end if
          call endrun()
       else
          call ncd_pio_openfile (pioid, trim(flandfrac), 0)
          call ncd_io(ncid=pioid, varname='landmask', data=lndmask_glob, flag='read')
          call ncd_io(ncid=pioid, varname='landfrac', data=lndfrac_glob, flag='read')
          call ncd_pio_closefile(pioid)
       end if

    end if

  end subroutine lnd_set_read_write_landmask

  !===============================================================================
  subroutine decompInit_ocn(ni, nj, amask, gindex_ocn)

    ! !DESCRIPTION:
    ! calculate a decomposition of only ocn points (needed for the nuopc interface)

    ! !USES:
    use spmdMod  , only : npes, iam

    ! !ARGUMENTS:
    integer , intent(in) :: amask(:)
    integer , intent(in) :: ni,nj   ! domain global size
    integer , pointer, intent(out) :: gindex_ocn(:) ! this variable is allocated here, and is assumed to start unallocated

    ! !LOCAL VARIABLES:
    integer :: n,i,j,nocn
    integer :: nlnd_global
    integer :: nocn_global
    integer :: nocn_local
    integer :: my_ocn_start, my_ocn_end
    !------------------------------------------------------------------------------

    ! count total land and ocean gridcells
    nlnd_global = 0
    nocn_global = 0
    do n = 1,ni*nj
       if (amask(n) == 1) then
          nlnd_global = nlnd_global + 1
       else
          nocn_global = nocn_global + 1
       endif
    enddo

    ! create the a global index array for ocean points
    nocn_local = nocn_global / npes

    my_ocn_start = nocn_local*iam + min(iam, mod(nocn_global, npes)) + 1
    if (iam < mod(nocn_global, npes)) then
       nocn_local = nocn_local + 1
    end if
    my_ocn_end = my_ocn_start + nocn_local - 1

    allocate(gindex_ocn(nocn_local))
    nocn = 0
    do n = 1,ni*nj
       if (amask(n) == 0) then
          nocn = nocn + 1
          if (nocn >= my_ocn_start .and. nocn <= my_ocn_end) then
             gindex_ocn(nocn - my_ocn_start + 1) = n
          end if
       end if
    end do
  end subroutine decompInit_ocn

end module lnd_set_decomp_and_domain
