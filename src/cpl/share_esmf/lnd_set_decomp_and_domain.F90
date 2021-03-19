module lnd_set_decomp_and_domain

  use ESMF
  use shr_kind_mod , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod  , only : shr_sys_abort
  use spmdMod      , only : masterproc
  use clm_varctl   , only : iulog

  implicit none
  private ! except

  ! Module public routines
  public :: lnd_set_decomp_and_domain_from_readmesh    ! nuopc/cmeps
  public :: lnd_set_decomp_and_domain_from_createmesh  ! nuopc/cmeps

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

    use decompInitMod , only : decompInit_ocn, decompInit_lnd, decompInit_lnd3D
    use domainMod     , only : ldomain, domain_init
    use decompMod     , only : ldecomp, bounds_type, get_proc_bounds
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
    type(ESMF_Mesh)     :: mesh_maskinput
    type(ESMF_Mesh)     :: mesh_lndinput
    type(ESMF_DistGrid) :: distgrid_ctsm
    integer             :: g,n             ! indices
    integer             :: nlnd, nocn      ! local size of arrays
    integer             :: gsize           ! global size of grid
    logical             :: isgrid2d        ! true => grid is 2d
    type(bounds_type)   :: bounds          ! bounds
    integer             :: begg,endg       ! local bounds
    integer  , pointer  :: gindex_lnd(:)   ! global index space for just land points
    integer  , pointer  :: gindex_ocn(:)   ! global index space for just ocean points
    integer  , pointer  :: gindex_ctsm(:)  ! global index space for land and ocean points
    integer  , pointer  :: lndmask_glob(:)
    real(r8) , pointer  :: lndfrac_glob(:)
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
    allocate(lndmask_glob(gsize)); lndmask_glob(:) = 0
    allocate(lndfrac_glob(gsize)); lndfrac_glob(:) = 0._r8

    ! Read in the land mesh from the file
    mesh_lndinput = ESMF_MeshCreate(filename=trim(meshfile_lnd), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (trim(driver) == 'cmeps') then
       ! Read in mask meshfile if needed
       if (trim(meshfile_mask) /= trim(meshfile_lnd)) then
          mesh_maskinput = ESMF_MeshCreate(filename=trim(meshfile_mask), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Determine lndmask_glob and lndfrac_glob
       if (trim(meshfile_mask) /= trim(meshfile_lnd)) then
          ! obain land mask and land fraction by mapping ocean mesh conservatively to land mesh
          call lnd_set_lndmask_from_maskmesh(mesh_lndinput, mesh_maskinput, vm, gsize, lndmask_glob, lndfrac_glob, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          ! obtain land mask from land mesh file - assume that land frac is identical to land mask
          call lnd_set_lndmask_from_lndmesh(mesh_lndinput, vm, gsize, lndmask_glob, lndfrac_glob, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    else if (trim(driver) == 'lilac') then
       call lnd_set_lndmask_from_fatmlndfrc(lndmask_glob, lndfrac_glob, ni,nj)
    else
       call shr_sys_abort('driver '//trim(driver)//' is not supported, must be lilac or cmeps')
    end if

    ! Determine lnd decomposition that will be used by ctsm from lndmask_glob
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

    ! Deallocate global pointer memory
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

    ! Set ldomain%lonc, ldomain%latc and ldomain%area
    call lnd_set_ldomain_gridinfo_from_mesh(mesh_ctsm, vm, gindex_ctsm, begg, endg, isgrid2d, ni, nj, ldomain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Deallocate local pointer memory
    deallocate(gindex_lnd)
    deallocate(gindex_ocn)
    deallocate(gindex_ctsm)

  end subroutine lnd_set_decomp_and_domain_from_readmesh

  !===============================================================================
  subroutine lnd_set_decomp_and_domain_from_createmesh(domain_file, vm, mesh_ctsm, ni, nj, rc)

    ! Generate a new mesh from the input domain file and set the mask to 1

    use decompInitMod , only : decompInit_lnd, decompInit_lnd3D
    use decompMod     , only : ldecomp, bounds_type, get_proc_bounds
    use domainMod     , only : ldomain, domain_init
    use clm_varctl    , only : use_soil_moisture_streams
    use clm_varctl    , only : scmlat, scmlon, single_column
    use clm_varpar    , only : nlevsoi
    use ncdio_pio     , only : pio_subsystem, io_type
    use pio

    ! input/output variables
    character(len=CL)   , intent(in)  :: domain_file
    type(ESMF_VM)       , intent(in)  :: vm
    type(ESMF_Mesh)     , intent(out) :: mesh_ctsm
    integer             , intent(out) :: ni,nj  ! global grid dimensions
    integer             , intent(out) :: rc

    ! local variables
    type(ESMF_Grid)       :: lgrid
    type(ESMF_Mesh)       :: mesh_lndcreate
    type(ESMF_DistGrid)   :: distgrid_ctsm
    integer, pointer      :: gindex_ctsm(:)       ! global index space for just land points
    logical               :: isgrid2d
    integer               :: i,j,g,n
    integer               :: nv
    integer               :: ierr
    integer               :: dimid
    integer               :: varid_xv, varid_yv
    integer               :: varid_xc, varid_yc
    integer               :: varid_area
    real(r8), allocatable :: xc(:,:), yc(:,:)     ! coordinates of centers
    real(r8), allocatable :: xv(:,:,:), yv(:,:,:) ! coordinates of corners
    integer               :: maxIndex(2)
    real(r8)              :: mincornerCoord(2)
    real(r8)              :: maxcornerCoord(2)
    integer               :: spatialDim
    integer               :: numownedelements
    real(r8) , pointer    :: ownedElemCoords(:)
    integer, allocatable  :: lnd_mask(:)
    type(bounds_type)     :: bounds               ! bounds
    integer               :: begg,endg
    integer               :: nlnd
    integer               :: start(2)             ! start index to read in for single column mode
    integer               :: count(2)             ! number of points to read in
    real(r8)              :: scol_data(1)         ! temporary
    integer , allocatable :: mask(:)              ! temporary
    real(r8), allocatable :: lats(:)              ! temporary
    real(r8), allocatable :: lons(:)              ! temporary
    real(r8), allocatable :: pos_lons(:)          ! temporary
    real(r8)              :: pos_scmlon           ! temporary
    real(r8)              :: scol_area            ! temporary
    type(file_desc_t)     :: pioid
    integer               :: rcode                ! error code
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    rcode = pio_openfile(pio_subsystem, pioid, io_type, trim(domain_file), pio_nowrite)
    call pio_check_err(rcode, 'error opening file '//trim(domain_file))
    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
    rcode = pio_inq_dimid(pioid, 'ni', dimid)
    call pio_check_err(rcode, 'pio_inq_dimid for ni in file '//trim(domain_file))
    rcode = pio_inquire_dimension(pioid, dimid, len=ni)
    call pio_check_err(rcode, 'pio_inq_dimension for ni in file '//trim(domain_file))
    rcode = pio_inq_dimid(pioid, 'nj', dimid)
    call pio_check_err(rcode, 'pio_inq_dimid for nj in file '//trim(domain_file))
    rcode = pio_inquire_dimension(pioid, dimid, len=nj)
    call pio_check_err(rcode, 'pio_inq_dimension for nj in file '//trim(domain_file))
    rcode = pio_inq_dimid(pioid, 'nv', dimid)
    call pio_check_err(rcode, 'pio_inq_dimid for nv in file '//trim(domain_file))
    rcode = pio_inquire_dimension(pioid, dimid, len=nv)
    call pio_check_err(rcode, 'pio_inq_dimension for nv in file '//trim(domain_file))
    rcode = pio_inq_varid(pioid, 'xc' , varid_xc)
    call pio_check_err(rcode, 'pio_inq_varid for yc in file '//trim(domain_file))
    rcode = pio_inq_varid(pioid, 'yc' , varid_yc)
    call pio_check_err(rcode, 'pio_inq_varid for yc in file '//trim(domain_file))
    rcode = pio_inq_varid(pioid, 'xv' , varid_xv)
    call pio_check_err(rcode, 'pio_inq_varid for xv in file '//trim(domain_file))
    rcode = pio_inq_varid(pioid, 'yv' , varid_yv)
    call pio_check_err(rcode, 'pio_inq_varid for yv in file '//trim(domain_file))
    rcode = pio_inq_varid(pioid, 'area', varid_area)
    call pio_check_err(rcode, 'pio_inq_varid for area in file '//trim(domain_file))

    if (single_column) then

       ! In this case the domain file is not a single point file - but normally a
       ! global domain file where a nearest neighbor search will be done to find
       ! the closest point in the domin file to scol_lon and scol_lat

       ! get center lats and lons from domain file
       allocate(xc(ni,nj))
       allocate(yc(ni,nj))
       rcode = pio_get_var(pioid, varid_xc, xc)
       call pio_check_err(rcode, 'pio_get_var for xc in file '//trim(domain_file))
       rcode = pio_get_var(pioid, varid_yc, yc)
       call pio_check_err(rcode, 'pio_get_var for yc in file '//trim(domain_file))

       ! find nearest neighbor indices of scmlon and scmlat in domain file
       allocate(lats(nj))
       allocate(lons(ni))
       allocate(pos_lons(ni))
       do i = 1,ni
          lons(i) = xc(i,1)
       end do
       do j = 1,nj
          lats(j) = yc(1,j)
       end do
       pos_lons(:)  = mod(lons(:)  + 360._r8, 360._r8)
       pos_scmlon = mod(scmlon + 360._r8, 360._r8)
       start(1) = (MINLOC(abs(pos_lons - pos_scmlon), dim=1))
       start(2) = (MINLOC(abs(lats      -scmlat    ), dim=1))
       count(:) = 1
       deallocate(lons)
       deallocate(lats)

       ! read in value of nearest neighbor lon and RESET scmlat
       rcode = pio_get_var(pioid, varid_xc, start, count, scol_data)
       call pio_check_err(rcode, 'pio_get_var for xc in file '//trim(domain_file))
       scmlon = scol_data(1)

       ! read in value of nearest neighbor lon and RESET scmlon
       rcode = pio_get_var(pioid, varid_yc, start, count, scol_data)
       call pio_check_err(rcode, 'pio_get_var for yc in file '//trim(domain_file))
       scmlat = scol_data(1)

       ! get area of gridcell
       rcode = pio_get_var(pioid, varid_area, start, count, scol_data)
       call pio_check_err(rcode, 'pio_get_var for area in file '//trim(domain_file))
       scol_area = scol_data(1)

       ! reset ni and nj to be single point values
       ni = 1
       nj = 1

       ! determine mincornerCoord and maxcornerCoord neede to create ESMF grid
       maxIndex(1)       = 1                        ! number of lons
       maxIndex(2)       = 1                        ! number of lats
       mincornerCoord(1) = scmlon - scol_area/2._r8 ! min lon
       mincornerCoord(2) = scmlat - scol_area/2._r8 ! min lat
       maxcornerCoord(1) = scmlon + scol_area/2._r8 ! max lon
       maxcornerCoord(2) = scmlat + scol_area/2._r8 ! max lat
       deallocate(xc,yc)

    else

       ! allocate xv and yv and read them in
       allocate(xv(nv,ni,nj))
       allocate(yv(nv,ni,nj))
       rcode = pio_get_var(pioid, varid_xv, xv)
       call pio_check_err(rcode, 'pio_get_var for xv in file '//trim(domain_file))
       rcode = pio_get_var(pioid, varid_yv, yv)
       call pio_check_err(rcode, 'pio_get_var for yv in file '//trim(domain_file))

       ! determine mincornerCoord and maxcornerCoord neede to create ESMF grid
       maxIndex(1)       = ni          ! number of lons
       maxIndex(2)       = nj          ! number of lats
       mincornerCoord(1) = xv(1,1,1)   ! min lon
       mincornerCoord(2) = yv(1,1,1)   ! min lat
       maxcornerCoord(1) = xv(3,ni,nj) ! max lon
       maxcornerCoord(2) = yv(3,ni,nj) ! max lat
       deallocate(xv,yv)

    end if

    ! close file
    call pio_seterrorhandling(pioid, PIO_INTERNAL_ERROR)
    call pio_closefile(pioid)

    ! create the ESMF grid
    lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
         mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
         staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create the mesh from the lgrid
    mesh_lndcreate =  ESMF_MeshCreate(lgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set the mesh mask to 1
    call ESMF_MeshGet(mesh_lndcreate, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    call ESMF_MeshGet(mesh_lndcreate, ownedElemCoords=ownedElemCoords, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lnd_mask(numownedelements))
    lnd_mask(:) = 1
    ! call ESMF_MeshSet(mesh_lndcreate, elementMask=lnd_mask, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine ldecomp and ldomain
    call decompInit_lnd(lni=ni, lnj=nj, amask=lnd_mask)
    if (use_soil_moisture_streams) then
       call decompInit_lnd3D(lni=ni, lnj=nj, lnk=nlevsoi)
    end if

    ! Initialize processor bounds
    call get_proc_bounds(bounds)
    begg = bounds%begg
    endg = bounds%endg

    ! Create gindex_ctsm
    nlnd = endg - begg + 1
    allocate(gindex_ctsm(nlnd))
    do g = begg, endg
       n = 1 + (g - begg)
       gindex_ctsm(n) = ldecomp%gdc2glo(g)
    end do

    ! Initialize domain data structure
    isgrid2d = .true.
    call domain_init(domain=ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine ldomain%mask and ldomain%frac
    do g = begg, endg
       ldomain%mask(g) = 1
       ldomain%frac(g) = 1._r8
    end do

    ! Generate a new mesh on the gindex decomposition
    distGrid_ctsm = ESMF_DistGridCreate(arbSeqIndexList=gindex_ctsm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    mesh_ctsm = ESMF_MeshCreate(mesh_lndcreate, elementDistGrid=distgrid_ctsm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set ldomain%lonc, ldomain%latc and ldomain%area
    call lnd_set_ldomain_gridinfo_from_mesh(mesh_ctsm, vm, gindex_ctsm, begg, endg, isgrid2d, ni, nj, ldomain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(lnd_mask)

  end subroutine lnd_set_decomp_and_domain_from_createmesh

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
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    use shr_sys_mod , only : shr_sys_abort

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
       write(iulog,*) 'Attempting to global dimensions from surface dataset'
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
  subroutine lnd_set_lndmask_from_maskmesh(mesh_lnd, mesh_mask, vm, gsize, lndmask_glob, lndfrac_glob, rc)

    ! input/out variables
    type(ESMF_Mesh)     , intent(in)  :: mesh_lnd
    type(ESMF_Mesh)     , intent(in)  :: mesh_mask
    type(ESMF_VM)       , intent(in)  :: vm
    integer             , intent(in)  :: gsize
    integer             , pointer     :: lndmask_glob(:)
    real(r8)            , pointer     :: lndfrac_glob(:)
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
    real(r8) , pointer     :: rtemp_glob(:)
    real(r8) , pointer     :: lndfrac_loc(:)
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
    deallocate(maskmask_loc)
    deallocate(lndmask_loc)
    deallocate(lndfrac_loc)

  end subroutine lnd_set_lndmask_from_maskmesh

  !===============================================================================
  subroutine lnd_set_lndmask_from_lndmesh(mesh_lnd, vm, gsize, lndmask_glob, lndfrac_glob, rc)

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

  end subroutine lnd_set_lndmask_from_lndmesh

  !===============================================================================
  subroutine lnd_set_lndmask_from_fatmlndfrc(mask, frac, ni, nj)

    ! Read the surface dataset grid related information
    ! This is used to set the domain decomposition - so global data is read here

    use clm_varctl , only : fatmlndfrc
    use fileutils  , only : getfil
    use ncdio_pio  , only : ncd_io, ncd_pio_openfile, ncd_pio_closefile, ncd_inqfdims, file_desc_t
    use abortutils , only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg

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
    character(len=32) :: subname = 'lnd_set_mask_from_fatmlndfrc' ! subroutine name
    !-----------------------------------------------------------------------

    ! Open file
    call getfil( fatmlndfrc, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Determine dimensions and if grid file is 2d or 1d
    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    if (masterproc) then
       write(iulog,*)'lat/lon grid flag (isgrid2d) is ',isgrid2d
    end if

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

    ! for reading in fatmlndfrc to override mesh data
    use clm_varctl , only : fatmlndfrc
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

    ! for sanity check - remove when this is done
    type(file_desc_t) :: ncid  ! netcdf id
    character(len=CL) :: locfn ! local file name
    real(r8), pointer :: lonc_atmlndfrc(:)
    real(r8), pointer :: latc_atmlndfrc(:)
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

    ! Sanity check- remove this when it is done
    call getfil( trim(fatmlndfrc), locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    allocate(lonc_atmlndfrc(numownedelements))
    allocate(latc_atmlndfrc(numownedelements))
    call ncd_io(ncid=ncid, varname= 'xc' , flag='read', data=lonc_atmlndfrc , dim1name=grlnd)
    call ncd_io(ncid=ncid, varname= 'yc' , flag='read', data=latc_atmlndfrc , dim1name=grlnd)
    do g = begg,endg
       n = g - begg + 1
       if ( abs(lonc_atmlndfrc(n) - ldomain%lonc(g)) > 1.e-11 .and. &
            abs(lonc_atmlndfrc(n) - ldomain%lonc(g)) /= 360._r8) then
          write(6,'(a,3(d20.13,2x))')'ERROR: lonc_atmlndfrac(n), ldomain%lonc(g), abs(diff) = ',&
               lonc_atmlndfrc(n), ldomain%lonc(g), abs(lonc_atmlndfrc(n) - ldomain%lonc(g))
          call shr_sys_abort()
       end if
       if (abs(latc_atmlndfrc(n) - ldomain%latc(g)) > 1.e-11) then
          write(6,'(a,3(d20.13,2x))')'ERROR: latc_atmlndfrac(n), ldomain%latc(g), abs(diff) = ',&
               latc_atmlndfrc(n), ldomain%latc(g), abs(latc_atmlndfrc(n) - ldomain%latc(g))
          call shr_sys_abort()
       end if
    end do
    deallocate(lonc_atmlndfrc)
    deallocate(latc_atmlndfrc)
    call ncd_pio_closefile(ncid)

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
  subroutine lnd_set_read_write_landmask(write_file, read_file, lndmask_glob, lndfrac_glob, gsize)

    ! This subroutine is currently unused (as of 2021-03-17), but it may be needed in the
    ! future. Its purpose is: Now that we get landmask and landfrac at runtime, from
    ! mapping the ocean mask to the land grid, it's possible that landfrac will be
    ! roundoff-level different with different processor counts. Mariana Vertenstein
    ! hasn't seen this happen yet, but if it does, then we can use this subroutine to
    ! solve this issue in tests that change processor count (ERP, PEM). I think Mariana's
    ! intent was: in the first run, we would write landmask and landfrac to a landfrac.nc
    ! file; then, in the second run (with different processor count), we would read that
    ! file rather than doing the mapping again. This way, both runs of the ERP or PEM
    ! test would use consistent landmask and landfrac values.

    use ncdio_pio , only : ncd_io, file_desc_t, ncd_pio_openfile, ncd_pio_closefile
    use ncdio_pio , only : ncd_defdim, ncd_defvar, ncd_enddef, ncd_inqdlen
    use ncdio_pio , only : ncd_int, ncd_double, ncd_pio_createfile

    ! input/output variables
    logical          , intent(in) :: write_file
    logical          , intent(in) :: read_file
    integer          , pointer    :: lndmask_glob(:)
    real(r8)         , pointer    :: lndfrac_glob(:)
    integer          , intent(in) :: gsize

    ! local variables
    type(file_desc_t) :: pioid ! netcdf file id
    integer           :: dimid
    character(len=CL) :: flandfrac = 'landfrac.nc'
    !-------------------------------------------------------------------------------

    if (write_file) then
       if (masterproc) then
          write(iulog,*)
          write(iulog,'(a)') 'lnd_set_decomp_and_domain: writing landmask and landfrac data to landfrac.nc'
          write(iulog,*)
       end if
       call ncd_pio_createfile(pioid, trim(flandfrac))
       call ncd_defdim (pioid, 'gridcell', gsize, dimid)
       call ncd_defvar(ncid=pioid, varname='landmask', xtype=ncd_int   , dim1name='gridcell')
       call ncd_defvar(ncid=pioid, varname='landfrac', xtype=ncd_double, dim1name='gridcell')
       call ncd_enddef(pioid)
       call ncd_io(ncid=pioid, varname='landmask', data=lndmask_glob, flag='write')
       call ncd_io(ncid=pioid, varname='landfrac', data=lndfrac_glob, flag='write')
       call ncd_pio_closefile(pioid)
    else if (read_file) then
       if (masterproc) then
          write(iulog,*)
          write(iulog,'(a)') 'lnd_set_decomp_and_domain: reading landmask and landfrac data from landfrac.nc'
          write(iulog,*)
       end if
       call ncd_pio_openfile (pioid, trim(flandfrac), 0)
       call ncd_io(ncid=pioid, varname='landmask', data=lndmask_glob, flag='read')
       call ncd_io(ncid=pioid, varname='landfrac', data=lndfrac_glob, flag='read')
       call ncd_pio_closefile(pioid)
    end if

  end subroutine lnd_set_read_write_landmask


end module lnd_set_decomp_and_domain
