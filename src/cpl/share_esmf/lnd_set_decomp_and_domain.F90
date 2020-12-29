module lnd_set_decomp_and_domain

  use ESMF
  use shr_kind_mod      , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use spmdMod           , only : masterproc
  use clm_varctl        , only : iulog

  implicit none
  private ! except

  ! Module public routines
  public :: lnd_set_decomp_and_domain_from_readmesh
  public :: lnd_set_decomp_and_domain_from_newmesh

  ! Module private routines
  private :: lnd_get_global_dims
  private :: lnd_get_lndmask_from_ocnmesh
  private :: lnd_get_lndmask_from_lndmesh
  private :: lnd_set_ldomain_gridinfo
  private :: nc_check_err
  private :: chkerr 

  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
!===============================================================================
contains
!===============================================================================

  subroutine lnd_set_decomp_and_domain_from_readmesh(mode, vm, meshfile_lnd, meshfile_ocn, mesh_ctsm, ni, nj, rc)

    use decompInitMod , only : decompInit_ocn, decompInit_lnd, decompInit_lnd3D
    use domainMod     , only : ldomain, domain_init
    use decompMod     , only : ldecomp, bounds_type, get_proc_bounds
    use clm_varpar    , only : nlevsoi
    use clm_varctl    , only : fatmlndfrc, fsurdat
    use clm_varctl    , only : use_soil_moisture_streams, single_column
    !
    use ncdio_pio     , only : ncd_io, file_desc_t, ncd_pio_openfile, ncd_pio_closefile, ncd_inqdlen
    use abortutils    , only : endrun
    use shr_log_mod   , only : errMsg => shr_log_errMsg
    use fileutils     , only : getfil

    ! input/output variables
    character(len=*)    , intent(in)    :: mode  ! lilac or nuopc mode
    type(ESMF_VM)       , intent(in)    :: vm
    character(len=*)    , intent(in)    :: meshfile_lnd
    character(len=*)    , intent(in)    :: meshfile_ocn
    type(ESMF_Mesh)     , intent(out)   :: mesh_ctsm
    integer             , intent(out)   :: ni,nj ! global grid dimensions
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Mesh)        :: mesh_ocninput
    type(ESMF_Mesh)        :: mesh_lndinput
    type(ESMF_DistGrid)    :: distgrid_ctsm
    character(CL)          :: cvalue          ! config data
    integer                :: nlnd, nocn      ! local size of arrays
    integer                :: g,n             ! indices
    type(bounds_type)      :: bounds          ! bounds
    integer                :: begg,endg
    integer  , pointer     :: gindex_lnd(:)   ! global index space for just land points
    integer  , pointer     :: gindex_ocn(:)   ! global index space for just ocean points
    integer  , pointer     :: gindex_ctsm(:)  ! global index space for land and ocean points
    integer  , pointer     :: gindex_input(:) ! global index space for land and ocean points
    integer  , pointer     :: lndmask_glob(:)
    real(r8) , pointer     :: lndfrac_glob(:)
    integer                :: lsize_input
    integer                :: gsize
    logical                :: isgrid2d
    character(len=CL)      :: locfn 
    type(file_desc_t)      :: ncid    ! netcdf file id
    integer                :: dimid   ! netCDF dimension id
    logical                :: readvar ! read variable in or not
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine global 2d sizes from read of dimensions of surface dataset
    call lnd_get_global_dims(ni, nj, gsize, isgrid2d)

    ! Allocate global memory for land mask and land fraction
    allocate(lndmask_glob(gsize)); lndmask_glob(:) = 0
    allocate(lndfrac_glob(gsize)); lndfrac_glob(:) = 0._r8

    ! read in the land mesh from the file
    mesh_lndinput = ESMF_MeshCreate(filename=trim(meshfile_lnd), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (masterproc) then
       write(iulog,'(a)')'land mesh file ',trim(meshfile_lnd)
    end if

    ! Set global land fraction and global land mask across all processors
    if (mode == 'lilac' .and. trim(fatmlndfrc) /= 'null') then
       ! Note that is just for backwards compatibility 
       ! Read in global land mask and land fraction from fatmlndfrc
       call getfil( trim(fatmlndfrc), locfn, 0 )
       call ncd_pio_openfile (ncid, trim(locfn), 0)
       call ncd_io(ncid=ncid, varname='mask', data=lndmask_glob, flag='read', readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: variable mask not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))
       allocate(lndfrac_glob(ni*nj)); lndfrac_glob(:) = 0._r8
       call ncd_io(ncid=ncid, varname='frac', data=lndfrac_glob, flag='read', readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: variable frac not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))
       call ncd_pio_closefile(ncid)
    else
       if (trim(meshfile_ocn) /= 'null') then
          ! read in ocn mask meshfile
          mesh_ocninput = ESMF_MeshCreate(filename=trim(meshfile_ocn), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (masterproc) then
             write(iulog,'(a)')'ocean mesh file ',trim(meshfile_ocn)
          end if

          ! obain land mask and land fraction by mapping ocean mesh conservatively to land mesh
          call lnd_get_lndmask_from_ocnmesh(mesh_lndinput, mesh_ocninput, vm, gsize, lndmask_glob, lndfrac_glob, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          ! obtain land mask from land mesh file - assume that land frac is identical to land mask
          call lnd_get_lndmask_from_lndmesh(mesh_lndinput, vm, gsize, lndmask_glob, lndfrac_glob, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Determine lnd decomposition that will be used by ctsm
    call decompInit_lnd(lni=ni, lnj=nj, amask=lndmask_glob)
    if (use_soil_moisture_streams) then
       call decompInit_lnd3D(lni=ni, lnj=nj, lnk=nlevsoi)
    end if

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
       gindex_lnd(n) = ldecomp%gdc2glo(g)
    end do

    ! Initialize domain data structure
    call domain_init(domain=ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine ldomain%mask and ldomain%frac using ctsm decomposition
    do g = begg, endg
       n = 1 + (g - begg)
       ldomain%mask(g) = lndmask_glob(gindex_lnd(n))
       ldomain%frac(g) = lndfrac_glob(gindex_lnd(n))
    end do
    deallocate(lndmask_glob)
    deallocate(lndfrac_glob)

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

    ! Get ldomain%lonc, ldomain%latc and ldomain%area and optionally
    ! lon1d and lat1d if isgrid2d
    call lnd_set_ldomain_gridinfo(mesh_ctsm, vm, gindex_ctsm, bounds, isgrid2d, ni, nj, ldomain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Deallocate memory
    deallocate(gindex_lnd)
    deallocate(gindex_ocn)
    deallocate(gindex_ctsm)

  end subroutine lnd_set_decomp_and_domain_from_readmesh

  !===============================================================================
  subroutine lnd_set_decomp_and_domain_from_newmesh(domain_file, mesh, ni, nj, rc)

    ! Generate a new mesh from the input domain file and set the mask to 1

    use decompInitMod , only : decompInit_lnd, decompInit_lnd3D
    use decompMod     , only : ldecomp, bounds_type, get_proc_bounds
    use domainMod     , only : ldomain, domain_init
    use clm_varctl    , only : use_soil_moisture_streams, single_column
    use clm_varpar    , only : nlevsoi
    use netcdf        , only : nf90_open, nf90_nowrite, nf90_noerr, nf90_close, nf90_strerror
    use netcdf        , only : nf90_inq_dimid, nf90_inq_varid, nf90_get_var
    use netcdf        , only : nf90_inquire_dimension, nf90_inquire_variable

    ! input/output variables
    character(len=CL)   , intent(in)  :: domain_file
    type(ESMF_Mesh)     , intent(out) :: mesh
    integer             , intent(out) :: ni,nj  ! global grid dimensions
    integer             , intent(out) :: rc

    ! local variables
    logical               :: isgrid2d
    integer               :: g,n
    integer               :: nv
    integer               :: ncid, ierr
    integer               :: dimid_ni, dimid_nj, dimid_nv
    integer               :: maxIndex(2)
    real(r8)              :: mincornerCoord(2)
    real(r8)              :: maxcornerCoord(2)
    type(ESMF_Grid)       :: lgrid
    real(r8), allocatable :: xv(:,:,:), yv(:,:,:)
    integer               :: varid_xv, varid_yv
    integer               :: numownedelements
    integer, allocatable  :: lnd_mask(:)
    type(bounds_type)     :: bounds          ! bounds
    integer               :: begg,endg
    integer               :: nlnd
    integer, pointer      :: gindex_lnd(:)   ! global index space for just land points
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! open file
    ierr = nf90_open(domain_file, NF90_NOWRITE, ncid)
    call nc_check_err(ierr, 'nf90_open', trim(domain_file))
    ! get dimension ids
    ierr = nf90_inq_dimid(ncid, 'ni', dimid_ni)
    call nc_check_err(ierr, 'nf90_inq_dimid for ni', trim(domain_file))
    ierr = nf90_inq_dimid(ncid, 'nj', dimid_nj)
    call nc_check_err(ierr, 'nf90_inq_dimid for nj', trim(domain_file))
    ierr = nf90_inq_dimid(ncid, 'nv', dimid_nv)
    call nc_check_err(ierr, 'nf90_inq_dimid for nv', trim(domain_file))
    ! get dimension values
    ierr = nf90_inquire_dimension(ncid, dimid_ni, len=ni)
    call nc_check_err(ierr, 'nf90_inq_dimension for ni', trim(domain_file))
    ierr = nf90_inquire_dimension(ncid, dimid_nj, len=nj)
    call nc_check_err(ierr, 'nf90_inq_dimension for nj', trim(domain_file))
    ierr = nf90_inquire_dimension(ncid, dimid_nv, len=nv)
    call nc_check_err(ierr, 'nf90_inq_dimension for nv', trim(domain_file))
    ! get variable ids
    ierr = nf90_inq_varid(ncid, 'xv', varid_xv)
    call nc_check_err(ierr, 'nf90_inq_varid for xv', trim(domain_file))
    ierr = nf90_inq_varid(ncid, 'yv', varid_yv)
    call nc_check_err(ierr, 'nf90_inq_varid for yv', trim(domain_file))
    ! allocate memory for variables and get variable values
    allocate(xv(nv,ni,nj), yv(nv,ni,nj))
    ierr = nf90_get_var(ncid, varid_xv, xv)
    call nc_check_err(ierr, 'nf90_get_var for xv', trim(domain_file))
    ierr = nf90_get_var(ncid, varid_yv, yv)
    call nc_check_err(ierr, 'nf90_get_var for yv', trim(domain_file))
    ! close file
    ierr = nf90_close(ncid)
    call nc_check_err(ierr, 'nf90_close', trim(domain_file))
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

    ! Set the mesh mask to 1
    call ESMF_MeshGet(mesh, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lnd_mask(numownedelements))
    lnd_mask(:) = 1
    call ESMF_MeshSet(mesh, elementMask=lnd_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ldecomp and ldomain
    call decompInit_lnd(lni=ni, lnj=nj, amask=lnd_mask)
    if (use_soil_moisture_streams) then
       call decompInit_lnd3D(lni=ni, lnj=nj, lnk=nlevsoi)
    end if

    ! Initialize processor bounds
    call get_proc_bounds(bounds)
    begg = bounds%begg
    endg = bounds%endg

    ! Create ctsm gindex_lnd
    nlnd = endg - begg + 1
    allocate(gindex_lnd(nlnd))
    do g = begg, endg
       n = 1 + (g - begg)
       gindex_lnd(n) = ldecomp%gdc2glo(g)
    end do

    ! Initialize domain data structure
    isgrid2d = .true.
    call domain_init(domain=ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine ldomain%mask and ldomain%frac
    do g = begg, endg
       ldomain%mask(g) = 1
       ldomain%frac(g) = 1._r8
    end do
    deallocate(lnd_mask)

  end subroutine lnd_set_decomp_and_domain_from_newmesh

  !===============================================================================
  subroutine lnd_get_global_dims(ni, nj, gsize, isgrid2d)

    ! Determine global 2d sizes from read of dimensions of surface dataset

    use clm_varctl  , only : fsurdat, single_column
    use fileutils   , only : getfil
    use ncdio_pio   , only : ncd_io, file_desc_t, ncd_pio_openfile, ncd_pio_closefile, ncd_inqdlen
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg

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
    !-------------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to global dimensions from surface dataset'
       if (fsurdat == ' ') then
          write(iulog,*)'fsurdat must be specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    endif
    call getfil(fsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdlen(ncid, dimid, ni, 'lsmlon')
    call ncd_inqdlen(ncid, dimid, nj, 'lsmlat')
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
  subroutine lnd_get_lndmask_from_ocnmesh(mesh_lnd, mesh_ocn, vm, gsize, lndmask_glob, lndfrac_glob, rc)

    ! input/out variables
    type(ESMF_Mesh)     , intent(in)  :: mesh_lnd
    type(ESMF_Mesh)     , intent(in)  :: mesh_ocn
    type(ESMF_VM)       , intent(in)  :: vm
    integer             , intent(in)  :: gsize
    integer             , pointer     :: lndmask_glob(:)
    real(r8)            , pointer     :: lndfrac_glob(:)
    integer             , intent(out) :: rc

    ! local variables:
    type(ESMF_DistGrid)    :: distgrid_lnd
    type(ESMF_RouteHandle) :: rhandle_ocn2lnd
    type(ESMF_Field)       :: field_lnd
    type(ESMF_Field)       :: field_ocn
    type(ESMF_DistGrid)    :: distgrid_ocn
    integer  , pointer     :: gindex_input(:) ! global index space for land and ocean points
    integer  , pointer     :: lndmask_loc(:)
    integer  , pointer     :: itemp_glob(:)
    real(r8) , pointer     :: rtemp_glob(:)
    real(r8) , pointer     :: lndfrac_loc(:)
    real(r8) , pointer     :: ocnmask_loc(:) ! on ocean mesh
    real(r8) , pointer     :: ocnfrac_loc(:) ! on land mesh
    real(r8) , pointer     :: dataptr1d(:)
    type(ESMF_Array)       :: elemMaskArray
    integer                :: lsize_lnd
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
    call ESMF_FieldDestroy(field_lnd)
    call ESMF_FieldDestroy(field_ocn)

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

    ! Determine ldomain%frac using both input and ctsm decompositions
    ! lndfrac_glob is filled using the input decomposition and
    ! ldomin%frac is set using the ctsm decomposition
    allocate(rtemp_glob(gsize))
    do n = 1,lsize_lnd
       lndfrac_glob(gindex_input(n)) = lndfrac_loc(n)
    end do
    call ESMF_VMAllReduce(vm, sendData=lndfrac_glob, recvData=rtemp_glob, count=gsize, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    lndfrac_glob(:) = rtemp_glob(:)
    deallocate(rtemp_glob)

    ! deallocate memory
    deallocate(ocnmask_loc)
    deallocate(lndmask_loc)
    deallocate(lndfrac_loc)

  end subroutine lnd_get_lndmask_from_ocnmesh

  !===============================================================================
  subroutine lnd_get_lndmask_from_lndmesh(mesh_lnd, vm, gsize, lndmask_glob, lndfrac_glob, rc)

    ! input/out variables
    type(ESMF_Mesh)     , intent(in)  :: mesh_lnd
    type(ESMF_VM)       , intent(in)  :: vm
    integer             , intent(in)  :: gsize
    integer             , pointer     :: lndmask_glob(:)
    real(r8)            , pointer     :: lndfrac_glob(:)
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
    do n = 1,lsize
       lndmask_glob(gindex(n)) = lndmask_loc(n)
    end do
    call ESMF_VMAllReduce(vm, sendData=lndmask_glob, recvData=itemp_glob, count=gsize, &
         reduceflag=ESMF_REDUCE_SUM, rc=rc)
    lndmask_glob(:) = int(itemp_glob(:))
    deallocate(itemp_glob)
    deallocate(gindex)
    deallocate(lndmask_loc)

    ! ASSUME that land fraction is identical to land mask in this case
    lndfrac_glob(:) = lndmask_glob(:)

  end subroutine lnd_get_lndmask_from_lndmesh
 
  !===============================================================================
  subroutine lnd_set_ldomain_gridinfo(mesh, vm, gindex, bounds, isgrid2d, ni, nj, ldomain, rc)

    use domainMod  , only : domain_type, lon1d, lat1d
    use decompMod  , only : bounds_type, get_proc_bounds
    use clm_varcon , only : re

    ! input/output variables
    type(ESMF_Mesh)   , intent(in)    :: mesh
    type(ESMF_VM)     , intent(in)    :: vm
    integer           , intent(in)    :: gindex(:)
    type(bounds_type) , intent(in)    :: bounds
    logical           , intent(in)    :: isgrid2d
    integer           , intent(in)    :: ni,nj 
    type(domain_type) , intent(inout) :: ldomain
    integer           , intent(out)   :: rc 

    ! local variables
    integer                :: g,n             
    integer                :: gsize
    integer                :: begg,endg
    integer                :: numownedelements
    real(r8) , pointer     :: lndlats_glob(:)
    real(r8) , pointer     :: lndlons_glob(:)
    real(r8) , pointer     :: rtemp_glob(:)
    real(r8) , pointer     :: ownedElemCoords(:)
    integer                :: spatialDim
    real(r8) , pointer     :: dataptr1d(:)
    type(ESMF_Field)       :: areaField
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    begg = bounds%begg
    endg = bounds%endg

    ! Determine ldoman%latc and ldomain%lonc
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
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

    ! If grid is 2d, determine lon1d and lat1d
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

  end subroutine lnd_set_ldomain_gridinfo

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

end module lnd_set_decomp_and_domain
