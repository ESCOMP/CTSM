module lnd_set_decomp_and_domain

  use ESMF
  use shr_kind_mod      , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use spmdMod           , only : masterproc
  use clm_varctl        , only : iulog
  use nuopc_shr_methods , only : chkerr

  implicit none
  private ! except

  ! Module public routines
  public :: lnd_set_decomp_and_domain_from_meshinfo
  public :: lnd_set_decomp_and_domain_from_newmesh

  ! Module private routines
  private :: clm_getlandmask_from_ocnmesh
  private :: clm_getlandmask_from_lndmesh
  private :: nc_check_err

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine lnd_set_decomp_and_domain_from_meshinfo(gcomp, mesh, ni, nj, rc)

    use NUOPC         , only : NUOPC_CompAttributeGet
    use decompInitMod , only : decompInit_ocn, decompInit_lnd, decompInit_lnd3D
    use domainMod     , only : ldomain, domain_init, lon1d, lat1d
    use decompMod     , only : ldecomp, bounds_type, get_proc_bounds
    use clm_varpar    , only : nlevsoi
    use clm_varctl    , only : use_soil_moisture_streams, single_column
    use clm_varcon    , only : re
    use lnd_comp_shr  , only : model_meshfile, model_clock

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(out)   :: mesh
    integer             , intent(out)   :: ni,nj ! global grid dimensions
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_Mesh)        :: mesh_lnd
    type(ESMF_Mesh)        :: mesh_ocn
    type(ESMF_DistGrid)    :: distgrid_mesh
    type(ESMF_DistGrid)    :: distgrid_lnd
    character(CL)          :: cvalue         ! config data
    integer                :: nlnd, nocn     ! local size of arrays
    integer                :: g,n            ! indices
    type(bounds_type)      :: bounds         ! bounds
    integer                :: begg,endg
    character(CL)          :: meshfile_ocn
    integer  , pointer     :: gindex_lnd(:)  ! global index space for just land points
    integer  , pointer     :: gindex_ocn(:)  ! global index space for just ocean points
    integer  , pointer     :: gindex(:)      ! global index space for land and ocean points
    integer  , pointer     :: mask(:)        ! local land/ocean mask
    integer  , pointer     :: lndmask_loc(:)
    real(r8) , pointer     :: lndfrac_loc(:)
    real(r8) , pointer     :: lndarea_loc(:)
    integer  , pointer     :: lndmask_glob(:)
    real(r8) , pointer     :: lndfrac_glob(:)
    real(r8) , pointer     :: lndarea_glob(:)
    real(r8) , pointer     :: lndlats_glob(:)
    real(r8) , pointer     :: lndlons_glob(:)
    real(r8) , pointer     :: rtemp_glob(:)
    integer  , pointer     :: itemp_glob(:)
    real(r8) , pointer     :: ownedElemCoords(:)
    real(r8) , pointer     :: dataptr1d(:)
    integer                :: lsize, gsize
    logical                :: isgrid2d
    integer                :: spatialDim
    type(ESMF_Field)       :: areaField
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! get vm
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global 2d sizes
    call NUOPC_CompAttributeGet(gcomp, name='lnd_ni', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ni
    call NUOPC_CompAttributeGet(gcomp, name='lnd_nj', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) nj
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

    ! read in the land mesh from the file
    mesh_lnd = ESMF_MeshCreate(filename=trim(model_meshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (masterproc) then
       write(iulog,'(a)')'land mesh file ',trim(model_meshfile)
    end if

    ! read in ocn mask meshfile
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ocnmask', value=meshfile_ocn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    mesh_ocn = ESMF_MeshCreate(filename=trim(meshfile_ocn), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (masterproc) then
       write(iulog,'(a)')'ocean mesh file ',trim(meshfile_ocn)
    end if

    ! set local land fraction and land mask for input read decomposition
    ! Note that lndmask_loc and lndfrac_loc are allocated in the following calls and lsize is returned
    if (trim(meshfile_ocn) == 'null') then
       ! obtain land mask from land mesh file - assume that land frac is identical to land mask
       call clm_getlandmask_from_lndmesh(mesh_lnd, lsize, lndmask_loc, lndfrac_loc, distgrid_lnd, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call clm_getlandmask_from_ocnmesh(mesh_lnd, mesh_ocn, lsize, lndmask_loc, lndfrac_loc, distgrid_lnd, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! determine global landmask_glob - needed to determine the ctsm decomposition
    ! land frac, lats, lons and areas will be done below
    allocate(gindex(lsize))
    call ESMF_DistGridGet(distgrid_lnd, 0, seqIndexList=gindex, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lndmask_glob(gsize)); lndmask_glob(:) = 0
    do n = 1,lsize
       lndmask_glob(gindex(n)) = lndmask_loc(n)
    end do
    allocate(itemp_glob(gsize))
    call ESMF_VMAllReduce(vm, sendData=lndmask_glob, recvData=itemp_glob, count=gsize, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    lndmask_glob(:) = int(itemp_glob(:))
    deallocate(itemp_glob)

    ! determine lnd decomposition that will be used by ctsm
    call decompInit_lnd(lni=ni, lnj=nj, amask=lndmask_glob)
    if (use_soil_moisture_streams) then
       call decompInit_lnd3D(lni=ni, lnj=nj, lnk=nlevsoi)
    end if

    ! Determine ocn decomposition that will be used to create the full mesh
    ! note that the memory for gindex_ocn will be allocated in the following call
    call decompInit_ocn(ni=ni, nj=nj, amask=lndmask_glob, gindex_ocn=gindex_ocn)

    ! *** Get JUST gridcell processor bounds ***
    ! Remaining bounds (landunits, columns, patches) will be set after calling decompInit_glcp
    ! so get_proc_bounds is called twice and the gridcell information is just filled in twice
    call get_proc_bounds(bounds)
    begg = bounds%begg
    endg = bounds%endg

    ! Create gindex_lnd 
    nlnd = endg - begg + 1
    allocate(gindex_lnd(nlnd))
    do g = begg, endg
       n = 1 + (g - begg)
       gindex_lnd(n) = ldecomp%gdc2glo(g)
    end do

    ! Initialize domain data structure
    call domain_init(domain=ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine ldomain%mask
    do g = begg, endg
       n = 1 + (g - begg)
       ldomain%mask(g) = lndmask_glob(gindex_lnd(n))
    end do
    deallocate(lndmask_glob)

    ! Determine ldomain%frac
    allocate(rtemp_glob(gsize))
    allocate(lndfrac_glob(gsize))
    lndfrac_glob(:) = 0._r8
    do n = 1,lsize
       lndfrac_glob(gindex(n)) = lndfrac_loc(n)
    end do
    call ESMF_VMAllReduce(vm, sendData=lndfrac_glob, recvData=rtemp_glob, count=gsize, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    lndfrac_glob(:) = rtemp_glob(:)
    do g = begg, endg
       ldomain%frac(g) = lndfrac_glob(gindex_lnd(g-begg+1))
    end do
    deallocate(lndfrac_glob)

    ! Get ownedElemCords from the mesh to be used to obtain ldoman%latc and ldomain%lonc
    call ESMF_MeshGet(mesh_lnd, spatialDim=spatialDim, rc=rc)
    allocate(ownedElemCoords(spatialDim*lsize))
    call ESMF_MeshGet(mesh_lnd, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ldomain%latc and global lat1d
    allocate(lndlats_glob(gsize))
    lndlats_glob(:) = 0._r8
    do n = 1,lsize
       lndlats_glob(gindex(n)) = ownedElemCoords(2*n)
    end do
    call ESMF_VMAllReduce(vm, sendData=lndlats_glob, recvData=rtemp_glob, count=gsize, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    lndlats_glob(:) = rtemp_glob(:)
    do g = begg, endg
       ldomain%latc(g) = lndlats_glob(gindex_lnd(g-begg+1))
    end do
    if (isgrid2d) then
       allocate(lat1d(nj))
       do n = 1,nj
          lat1d(n) = lndlats_glob((n-1)*ni + 1)
       end do
    end if
    deallocate(lndlats_glob)

    ! Determine ldomain%lonc and global lon1d
    allocate(lndlons_glob(gsize))
    lndlons_glob(:) = 0._r8
    do n = 1,lsize
       lndlons_glob(gindex(n)) = ownedElemCoords(2*n-1)
    end do
    call ESMF_VMAllReduce(vm, sendData=lndlons_glob, recvData=rtemp_glob, count=gsize, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    lndlons_glob(:) = rtemp_glob(:)
    do g = begg, endg
       ldomain%lonc(g) = lndlons_glob(gindex_lnd(g-begg+1))
    end do
    if (isgrid2d) then
       allocate(lon1d(ni))
       do n = 1,ni
          lon1d(n) = lndlons_glob(n)
       end do
    end if
    deallocate(lndlons_glob)
    deallocate(rtemp_glob)

    ! Create a global index that includes both land and ocean points
    nocn = size(gindex_ocn)
    allocate(gindex(nlnd + nocn))
    do n = 1,nlnd+nocn
       if (n <= nlnd) then
          gindex(n) = gindex_lnd(n)
       else
          gindex(n) = gindex_ocn(n-nlnd)
       end if
    end do

    ! Generate a new mesh on the gindex decomposition
    distGrid_mesh = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(gindex)
    mesh = ESMF_MeshCreate(mesh_lnd, elementDistGrid=distgrid_mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

  end subroutine lnd_set_decomp_and_domain_from_meshinfo

  !===============================================================================
  subroutine lnd_set_decomp_and_domain_from_newmesh(gcomp, mesh, ni, nj, rc)

    use NUOPC      , only : NUOPC_CompAttributeGet
    use clm_varctl , only : single_column
    use netcdf     , only : nf90_open, nf90_nowrite, nf90_noerr, nf90_close, nf90_strerror
    use netcdf     , only : nf90_inq_dimid, nf90_inq_varid, nf90_get_var
    use netcdf     , only : nf90_inquire_dimension, nf90_inquire_variable

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(out)   :: mesh
    integer             , intent(out)   :: ni,nj  ! global grid dimensions
    integer             , intent(out)   :: rc

    ! local variables
    integer                 :: ncid, ierr
    integer                 :: nv
    integer                 :: dimid_ni, dimid_nj, dimid_nv
    integer                 :: maxIndex(2)
    real(r8)                :: mincornerCoord(2)
    real(r8)                :: maxcornerCoord(2)
    type(ESMF_Grid)         :: lgrid
    real(r8), allocatable   :: xv(:,:,:), yv(:,:,:)
    integer                 :: varid_xv, varid_yv
    character(len=CL)       :: cvalue
    integer                 :: gsize
    logical                 :: isgrid2d
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! determine global 2d sizes
    call NUOPC_CompAttributeGet(gcomp, name='lnd_ni', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ni
    call NUOPC_CompAttributeGet(gcomp, name='lnd_nj', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) nj
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

    ! get the datm grid from the domain file
    call NUOPC_CompAttributeGet(gcomp, name='domain_lnd', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! open file
    ierr = nf90_open(cvalue, NF90_NOWRITE, ncid)
    call nc_check_err(ierr, 'nf90_open', trim(cvalue))
    ! get dimension ids
    ierr = nf90_inq_dimid(ncid, 'ni', dimid_ni)
    call nc_check_err(ierr, 'nf90_inq_dimid for ni', trim(cvalue))
    ierr = nf90_inq_dimid(ncid, 'nj', dimid_nj)
    call nc_check_err(ierr, 'nf90_inq_dimid for nj', trim(cvalue))
    ierr = nf90_inq_dimid(ncid, 'nv', dimid_nv)
    call nc_check_err(ierr, 'nf90_inq_dimid for nv', trim(cvalue))
    ! get dimension values
    ierr = nf90_inquire_dimension(ncid, dimid_ni, len=ni)
    call nc_check_err(ierr, 'nf90_inq_dimension for ni', trim(cvalue))
    ierr = nf90_inquire_dimension(ncid, dimid_nj, len=nj)
    call nc_check_err(ierr, 'nf90_inq_dimension for nj', trim(cvalue))
    ierr = nf90_inquire_dimension(ncid, dimid_nv, len=nv)
    call nc_check_err(ierr, 'nf90_inq_dimension for nv', trim(cvalue))
    ! get variable ids
    ierr = nf90_inq_varid(ncid, 'xv', varid_xv)
    call nc_check_err(ierr, 'nf90_inq_varid for xv', trim(cvalue))
    ierr = nf90_inq_varid(ncid, 'yv', varid_yv)
    call nc_check_err(ierr, 'nf90_inq_varid for yv', trim(cvalue))
    ! allocate memory for variables and get variable values
    allocate(xv(nv,ni,nj), yv(nv,ni,nj))
    ierr = nf90_get_var(ncid, varid_xv, xv)
    call nc_check_err(ierr, 'nf90_get_var for xv', trim(cvalue))
    ierr = nf90_get_var(ncid, varid_yv, yv)
    call nc_check_err(ierr, 'nf90_get_var for yv', trim(cvalue))
    ! close file
    ierr = nf90_close(ncid)
    call nc_check_err(ierr, 'nf90_close', trim(cvalue))
    ! create the grid
    maxIndex(1)       = ni          ! number of lons
    maxIndex(2)       = nj          ! number of lats
    mincornerCoord(1) = xv(1,1,1)   ! min lon
    mincornerCoord(2) = yv(1,1,1)   ! min lat
    maxcornerCoord(1) = xv(3,ni,nj) ! max lon
    maxcornerCoord(2) = yv(3,ni,nj) ! max lat
    deallocate(xv,yv)
    lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
         mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
         staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create the mesh from the grid
    mesh =  ESMF_MeshCreate(lgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO: initialize the decomposition
    ! initialize ldomain
    ! initialize the mask and mesh
    ! for created meshes assume the mask is 1
    ! create a pointer for mask and set it to 1
    ! call ESMF_MeshSet(mesh, elementMask=mask, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! deallocate(mask)

  end subroutine lnd_set_decomp_and_domain_from_newmesh

  !===============================================================================
  subroutine clm_getlandmask_from_ocnmesh(mesh_lnd, mesh_ocn, lsize_lnd, lndmask_loc, lndfrac_loc, distgrid_lnd, rc)

    ! input/out variables
    type(ESMF_Mesh)     , intent(in)  :: mesh_lnd
    type(ESMF_Mesh)     , intent(in)  :: mesh_ocn
    integer             , pointer     :: lndmask_loc(:)
    real(r8)            , pointer     :: lndfrac_loc(:)
    integer             , intent(out) :: lsize_lnd
    type(ESMF_DistGrid) , intent(out) :: distgrid_lnd
    integer             , intent(out) :: rc

    ! local variables:
    type(ESMF_RouteHandle) :: rhandle_ocn2lnd
    type(ESMF_Field)       :: field_lnd
    type(ESMF_Field)       :: field_ocn
    type(ESMF_DistGrid)    :: distgrid_ocn
    real(r8) , pointer     :: ocnmask_loc(:) ! on ocean mesh
    real(r8) , pointer     :: ocnfrac_loc(:) ! on land mesh
    real(r8) , pointer     :: dataptr1d(:)
    type(ESMF_Array)       :: elemMaskArray
    integer                :: lsize_ocn
    integer                :: n, spatialDim
    integer                :: srcMaskValue = 0
    integer                :: dstMaskValue = -987987 ! spval for RH mask values
    integer                :: srcTermProcessing_Value = 0
    real(r8)               :: fminval = 0.001_r8
    real(r8)               :: fmaxval = 1._r8
    logical                :: checkflag = .false.
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_MeshGet(mesh_lnd, spatialDim=spatialDim, numOwnedElements=lsize_lnd, &
         elementDistGrid=distgrid_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lndmask_loc(lsize_lnd))
    allocate(lndfrac_loc(lsize_lnd))

    ! create fields on land and ocean meshes
    field_lnd = ESMF_FieldCreate(mesh_lnd, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    field_ocn = ESMF_FieldCreate(mesh_ocn, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create route handle to map ocean mask from ocn mesh to land mesh
    call ESMF_FieldRegridStore(field_ocn, field_lnd, routehandle=rhandle_ocn2lnd, &
         srcMaskValues=(/srcMaskValue/), dstMaskValues=(/dstMaskValue/), &
         regridmethod=ESMF_REGRIDMETHOD_CONSERVE, normType=ESMF_NORMTYPE_DSTAREA, &
         srcTermProcessing=srcTermProcessing_Value, &
         ignoreDegenerate=.true., unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
          
    ! fill in values for field_ocn with mask on ocn mesh
    call ESMF_MeshGet(mesh_ocn, elementdistGrid=distgrid_ocn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distgrid_ocn, localDe=0, elementCount=lsize_ocn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ocnmask_loc(lsize_ocn))
    elemMaskArray = ESMF_ArrayCreate(distgrid_ocn, ocnmask_loc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_ocn, elemMaskArray=elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(field_ocn, farrayptr=dataptr1d, rc=rc)
    dataptr1d(:) = ocnmask_loc(:)

    ! map ocn mask to land mesh
    call ESMF_FieldRegrid(field_ocn, field_lnd, routehandle=rhandle_ocn2lnd, &
         termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=checkflag, zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(mesh_lnd, spatialDim=spatialDim, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ocnfrac_loc(lsize_lnd))
    call ESMF_FieldGet(field_lnd, farrayptr=ocnfrac_loc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,lsize_lnd
       lndfrac_loc(n) = 1._r8 - ocnfrac_loc(n)
       if (lndfrac_loc(n) > fmaxval) lndfrac_loc(n) = 1._r8
       if (lndfrac_loc(n) < fminval) lndfrac_loc(n) = 0._r8
       if (lndfrac_loc(n) /= 0._r8) then
          lndmask_loc(n) = 1
       else
          lndmask_loc(n) = 0
       end if
    enddo

    ! deallocate memory
    call ESMF_FieldDestroy(field_lnd)
    call ESMF_FieldDestroy(field_ocn)
    deallocate(ocnmask_loc)

  end subroutine clm_getlandmask_from_ocnmesh

  !===============================================================================
  subroutine clm_getlandmask_from_lndmesh(mesh_lnd, lsize, lndmask_loc, lndfrac_loc, distgrid_lnd, rc)

    ! input/out variables
    type(ESMF_Mesh)     , intent(in)  :: mesh_lnd
    integer             , intent(out) :: lsize
    integer             , pointer     :: lndmask_loc(:)
    real(r8)            , pointer     :: lndfrac_loc(:)
    type(ESMF_DistGrid) , intent(out) :: distgrid_lnd
    integer             , intent(out) :: rc

    ! local variables:
    type(ESMF_Array) :: elemMaskArray
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine lsize and distgrid_lnd
    call ESMF_MeshGet(mesh_lnd, elementdistGrid=distgrid_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distgrid_lnd, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine lndfrac_loc
    ! The call to ESMF_MeshGet fills in the values of lndmask_loc
    allocate(lndmask_loc(lsize))
    elemMaskArray = ESMF_ArrayCreate(distgrid_lnd, lndmask_loc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_lnd, elemMaskArray=elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine lndmask_loc
    ! ASSUME that land fraction is identical to land mask in this case
    allocate(lndfrac_loc(lsize))
    lndfrac_loc(:) = lndmask_loc(:)

  end subroutine clm_getlandmask_from_lndmesh

  !===============================================================================
  subroutine nc_check_err(ierror, description, filename)

    use shr_sys_mod , only : shr_sys_abort
    use netcdf      , only : nf90_noerr, nf90_strerror

    integer     , intent(in) :: ierror
    character(*), intent(in) :: description
    character(*), intent(in) :: filename

    if (ierror /= nf90_noerr) then
       write (*,'(6a)') 'ERROR ', trim(description),'. NetCDF file : "', trim(filename),&
            '". Error message:', trim(nf90_strerror(ierror))
       call shr_sys_abort()
    endif
  end subroutine nc_check_err

end module lnd_set_decomp_and_domain
