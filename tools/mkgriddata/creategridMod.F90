module creategridMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: creategridMod
!
! !DESCRIPTION:
! Routines to create land model grid
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use fileutils   , only : getfil
  use mkvarctl   
  use domainMod   , only : domain_init, domain_type, domain_check
  use areaMod
  use ncdio
  use shr_const_mod, only : SHR_CONST_PI, SHR_CONST_REARTH, SHR_CONST_G
  use shr_sys_mod    , only : shr_sys_flush
!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: creategrid    ! Generate land model grid.
  public :: settopo       ! Generate topography
  public :: read_domain   ! read domain from netcdf file
  public :: write_domain  ! write domain to netcdf file
  public :: mkfile        ! create netcdf file

! !PRIVATE MEMBER FUNCTIONS:
  real(r8) :: flandmin = 0.001             ! minimum land frac for land cell
  real(r8) :: re = SHR_CONST_REARTH*0.001  ! radius of earth (km)
  type (gridmap_type), public :: gridmap_d2l
  type (gridmap_type), public :: gridmap_t2l
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: creategrid
!
! !INTERFACE:
  subroutine creategrid(fname,type)
!
! !DESCRIPTION:
! Generate land model grid.
! Surface grid edges -- Grids do not have to be global. To allow this, grids
! must define the north, east, south, and west edges:
! namelist variables
!    o fnavyoro : 20 min navy orography dataset
!    o edgen (edge(1)) : northern edge of grid (degrees): >  -90 and <= 90
!    o edgee (edge(2)) : eastern edge of grid (degrees) : see following notes
!    o edges (edge(3)) : southern edge of grid (degrees): >= -90 and <  90
!    o edgew (edge(4)) : western edge of grid (degrees) : see following notes
! For partial grids, northern and southern edges are any latitude
! between 90 (North Pole) and -90 (South Pole). Western and eastern
! edges are any longitude between -180 and 180, with longitudes
! west of Greenwich negative. That is, western edge >= -180 and < 180;
! eastern edge > western edge and <= 180.
! For global grids, northern and southern edges are 90 (North Pole)
! and -90 (South Pole). The western edge of the longitude grid starts
! at the dateline if the grid is generated (the longitudes for each grid
! cell correspond with the edges (i.e., range from -180 to 180)).
!
! !USES:
  use mkvarsur
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:

    integer  :: i,j                            !indices
    integer  :: lsmlon, lsmlat                 !local size
    real(r8) :: dx                             !land model cell width
    real(r8) :: dy                             !land model cell length
    character(len= 32) :: subname = 'create_grid'
    type(domain_type) :: ddomain
!-----------------------------------------------------------------------

    if (trim(type) == 'internal') then
       if (mksrf_lsmlon==-1 .and. mksrf_lsmlat==-1) then
          write(6,*) 'must specify mksrf_lsmlon/lat with internal type'
          stop
       endif
       if (mksrf_edgen ==-999. .or. mksrf_edges ==-999. .or. mksrf_edgee ==-999. .or. mksrf_edgew ==-999.) then
          write(6,*) 'must specify mksrf_edgen/edges/edgen/edgew with internal type'
          stop
       endif

       call read_domain(ddomain,fname)

       lsmlon     = mksrf_lsmlon
       lsmlat     = mksrf_lsmlat

       call domain_init(ldomain,lsmlon,lsmlat)

       ldomain%numlon(:) = lsmlon
       ldomain%edgen = mksrf_edgen
       ldomain%edgee = mksrf_edgee
       ldomain%edges = mksrf_edges
       ldomain%edgew = mksrf_edgew

       dx = (ldomain%edgee - ldomain%edgew) / lsmlon
       dy = (ldomain%edgen - ldomain%edges) / lsmlat

       if (dx <= 0._r8 .or. dy <= 0._r8) then
          write(6,*) 'creategrid ERROR edges wrong sign NESN:',ldomain%edgen,ldomain%edgee,ldomain%edges,ldomain%edgen
          stop
       endif

       do j = 1, lsmlat
       do i = 1, lsmlon
          ldomain%longxy(i,j) = ldomain%edgew + (2*i-1)*dx / 2.
          ldomain%latixy(i,j) = ldomain%edges + (2*j-1)*dy / 2
       end do
       end do

       ! Define edges and area of output land model grid cells
       call celledge (ldomain, ldomain%edgen , ldomain%edgee , &
                               ldomain%edges , ldomain%edgew)
       call cellarea (ldomain, ldomain%edgen , ldomain%edgee , &
                               ldomain%edges , ldomain%edgew)

       call ddomain_to_ldomain(ddomain,ldomain)

       write(6,*) ' '
       write(6,*) trim(subname),'- ddomain:'
       call domain_check(ddomain)

    elseif (trim(type) == 'external') then
       call read_domain(ldomain,fname,trim(type))
    else
       write(6,*) 'creategrid ERROR, type = ',trim(type)
       stop
    endif

    write (6,*) 'Successfully made land grid data'
    write (6,*)
    write(6,*) trim(subname),'- ldomain:'
    call domain_check(ldomain)

  end subroutine creategrid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: settopo
!
! !INTERFACE:
  subroutine settopo(fname)
!
! !DESCRIPTION:
! Generate topo data.
!
! !USES:
  use mkvarsur
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:

    integer  :: i,j,im,jm                      !indices
    integer  :: lsmlon, lsmlat                 !local size
    real(r8) :: dx                             !land model cell width
    real(r8) :: dy                             !land model cell length
    character(len= 32) :: subname = 'settopo'
    real(r8),allocatable :: fld_o(:,:)      !output grid: dummy field
    real(r8),allocatable :: fld_i(:,:)      !input grid: dummy field
    integer ,allocatable :: mask(:,:)       !temp for ldomain mask
    type(domain_type) :: tdomain

!-----------------------------------------------------------------------

! establish tdomain
    call read_domain(tdomain,fname)

    call domain_check(tdomain)

    allocate(fld_i(tdomain%ni,tdomain%nj))
    allocate(fld_o(ldomain%ni,ldomain%nj))
    allocate(mask (ldomain%ni,ldomain%nj))
    fld_i = 1.0
    fld_o = 1.0

    mask = ldomain%mask
    tdomain%mask = 1
    ldomain%mask = 1

    write(6,*) 'call areaini'

    call areaini(tdomain,ldomain,gridmap_t2l,fracin=fld_i,fracout=fld_o)

    write(6,*) 'areaini done'

    ldomain%mask = mask

    write(6,*) ' '
    write(6,*) trim(subname),':'
    call gridmap_checkmap(gridmap_t2l)

    write(6,*) 'call areaave'

    call areaave(tdomain%topo,ldomain%topo,gridmap_t2l)

    write(6,*) 'areaave done'

    write(6,*) trim(subname),': raw topo min/max = ',minval(tdomain%topo),maxval(tdomain%topo)
    write(6,*) trim(subname),': new topo min/max = ',minval(ldomain%topo),maxval(ldomain%topo)

    deallocate(fld_i,fld_o,mask)

  end subroutine settopo

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ddomain_to_ldomain
!
! !INTERFACE:
  subroutine ddomain_to_ldomain(ddomain,ldomain)
!
! !DESCRIPTION:
! Read a grid file

! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: ddomain
    type(domain_type),intent(out)   :: ldomain
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j                         !indices
    real(r8),allocatable :: fld_o(:,:)      !output grid: dummy field
    real(r8),allocatable :: fld_i(:,:)      !input grid: dummy field
    real(r8) :: sum_fldo                    !global sum of dummy output field
    real(r8) :: sum_fldi                    !global sum of dummy input field
    real(r8) :: relerr = 0.00001            !max error for sum weights ne 1
    character(len= 32) :: subname = 'ddomain_to_ldomain'
!-----------------------------------------------------------------

    ddomain%mask = 1
    ldomain%mask = 1

    ! Determine output grid longitudes and latitudes in increments of dx and dy
    ! Global latitude grid goes from south pole to north pole
    ! Global longitude grid starts at Dateline with western edge on Dateline

    allocate(fld_i(ddomain%ni,ddomain%nj))
    allocate(fld_o(ldomain%ni,ldomain%nj))
    fld_i = 1.0
    fld_o = 1.0

    call areaini(ddomain,ldomain,gridmap_d2l,fracin=fld_i,fracout=fld_o)

    write(6,*) ' '
    write(6,*) trim(subname),':'
    call gridmap_checkmap(gridmap_d2l)

    call areaave(ddomain%frac,ldomain%frac,gridmap_d2l)

    if (minval(ldomain%frac) < -0.000001_r8 .or. &
        maxval(ldomain%frac) >  1.000001_r8) then
       write (6,*) 'MKGRID error: fland out of bounds [0,1] ', &
          minval(ldomain%frac),maxval(ldomain%frac)
       stop
    end if

    where (ldomain%frac(:,:) > 1.0_r8)
       ldomain%frac(:,:) = 1.0_r8
    endwhere

    where (ldomain%frac(:,:) < 0.0_r8)
       ldomain%frac(:,:) = 0.0_r8
    endwhere

    do j=1,ddomain%nj
    do i=1,ddomain%ni
       fld_i(i,j) = ((j-1)*ddomain%ni + i)
    enddo
    enddo

    call areaave(fld_i,fld_o,gridmap_d2l)

    ! -----------------------------------------------------------------
    ! Error check1
    ! Compare global sum fld_o to global sum fld_i.
    ! -----------------------------------------------------------------

    ! This check is true only if both grids span the same domain.
    ! To obtain global sum of input field must multiply by
    ! fraction of input grid that is land as determined by input grid

    sum_fldi = 0.
    do j = 1, ddomain%nj
    do i = 1, ddomain%ni
       sum_fldi = sum_fldi + ddomain%area(i,j) * fld_i(i,j)
    end do
    end do

    sum_fldo = 0.
    do j = 1,ldomain%nj
    do i = 1,ldomain%ni
       sum_fldo = sum_fldo + ldomain%area(i,j) * fld_o(i,j)
    end do
    end do

    if ( abs(ldomain%edgen - ldomain%edges) == 180. .and. &
         abs(ldomain%edgee - ldomain%edgew) == 360. ) then
       if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
          write (6,*) 'MKGRID error: input field not conserved'
          write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
          write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
          stop
       end if
    end if

    ! Determine land mask

    where (ldomain%frac(:,:) < flandmin)
       ldomain%mask(:,:) = 0     !ocean
    elsewhere
       ldomain%mask(:,:) = 1     !land
    endwhere

    ! Reset landfrac to zero where landmask has been set to zero

    where (ldomain%mask(:,:) == 0)
       ldomain%frac(:,:) = 0
    endwhere

    ! deallocate dynamic memory

    deallocate(fld_i,fld_o)

  end subroutine ddomain_to_ldomain

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_domain
!
! !INTERFACE:
  subroutine read_domain(domain,fname,type)
!
! !DESCRIPTION:
! Read a grid file
!
! !USES:
    implicit none
    include 'netcdf.inc'
!
! !ARGUMENTS:
    type(domain_type),intent(inout) :: domain
    character(len=*) ,intent(in)    :: fname
    character(len=*) ,intent(in),optional :: type
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nlon,nlat                       !size
    logical :: etype                           !external type
    real(r8), allocatable :: lon1d(:)          !local array for 1d lon
    real(r8), allocatable :: lat1d(:)          !local array for 1d lat
    real(r8), allocatable :: xv(:,:,:)         !local array for corner lons
    real(r8), allocatable :: yv(:,:,:)         !local array for corner lats
    character(len=256) :: locfn                !local file name
    character(len=256) :: fnamel               !local file name
    integer :: i,j                             !indexes
    integer :: ncid                            !netCDF file id
    integer :: dimid                           !netCDF dimension id
    integer :: varid                           !netCDF variable id
    logical :: dimset                          !local ni,nj
    logical :: lonlatset                       !local lon(:,:), lat(:,:)
    logical :: edgeNESWset                     !local EDGE[NESW]
    logical :: llneswset                       !local lat[ns],lon[we]
    logical :: areaset                         !local area
    logical :: toposet                         !local topo
    logical :: landfracset                     !local landfrac
    logical :: maskset                         !local mask
    logical :: numlonset                       !local numlon
    integer :: ndims                           !number of dims for variable
    integer :: ier
    character(len= 32) :: subname = 'read_domain'
!-----------------------------------------------------------------

    fnamel = fname

    etype = .false.
    if (present(type)) then
       if (trim(type) == 'external') etype = .true.
    endif
    write(6,*) ' ' 

    dimset      = .false.
    lonlatset   = .false.
    edgeNESWset = .false.
    llneswset   = .false.
    areaset     = .false.
    toposet     = .false.
    landfracset = .false.
    maskset     = .false.
    numlonset   = .false.

    ! Read domain file and compute stuff as needed

    call getfil (fnamel, locfn, 0)
    call check_ret(nf_open(locfn, 0, ncid), subname)

    ier = nf_inq_dimid (ncid, 'lon', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read lon and lat dims'
       call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
    endif

    ier = nf_inq_dimid (ncid, 'ni', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read ni and nj dims'
       call check_ret(nf_inq_dimid  (ncid, 'ni', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'nj', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
    endif

    ier = nf_inq_dimid (ncid, 'lsmlon', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read lsmlon and lsmlat dims'
       call check_ret(nf_inq_dimid  (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'lsmlat', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
    endif

    if (dimset) then
       write(6,*) trim(subname),' initialized domain'
       call domain_init(domain,nlon,nlat)
    else
       write(6,*) trim(subname),' ERROR: nlon, nlat not set for domain_init'
       stop
    endif

    ier = nf_inq_varid (ncid, 'xc', varid)
    if (ier == NF_NOERR) then
       if (lonlatset) write(6,*) trim(subname),' WARNING, overwriting lat,lon'
       lonlatset = .true.
       write(6,*) trim(subname),' read xc and yc fields'
       call check_ret(nf_inq_varid (ncid, 'xc', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%longxy), subname)
       call check_ret(nf_inq_varid (ncid, 'yc', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%latixy), subname)
    endif

    ier = nf_inq_varid (ncid, 'lon', varid)
    if (ier == NF_NOERR) then
       if (lonlatset) write(6,*) trim(subname),' WARNING, overwriting lat,lon'
       lonlatset = .true.
       ier = nf_inq_varndims(ncid,varid,ndims)
       if (ndims == 1) then
          write(6,*) trim(subname),' read lon and lat 1d fields'
          allocate(lon1d(nlon),lat1d(nlat))
          call check_ret(nf_inq_varid (ncid, 'lon', varid), subname)
          call check_ret(nf_get_var_double (ncid, varid, lon1d), subname)
          call check_ret(nf_inq_varid (ncid, 'lat', varid), subname)
          call check_ret(nf_get_var_double (ncid, varid, lat1d), subname)
          do j = 1, nlat
          do i = 1, nlon
             domain%longxy(i,j) = lon1d(i)
             domain%latixy(i,j) = lat1d(j)
          enddo
          enddo
          deallocate(lon1d,lat1d)
       elseif (ndims == 2) then
          write(6,*) trim(subname),' read lon and lat 2d fields'
          call check_ret(nf_inq_varid (ncid, 'lat', varid), subname)
          call check_ret(nf_get_var_double (ncid, varid, domain%latixy), subname)
          call check_ret(nf_inq_varid (ncid, 'lon', varid), subname)
          call check_ret(nf_get_var_double (ncid, varid, domain%longxy), subname)
       else
          write(6,*) trim(subname),'ERROR: lon and lat illegal dim ',ndims
          stop
       endif
    endif

    ier = nf_inq_varid (ncid, 'LATIXY', varid)
    if (ier == NF_NOERR) then
       if (lonlatset) write(6,*) trim(subname),' WARNING, overwriting lat,lon'
       lonlatset = .true.
       write(6,*) trim(subname),' read LONGXY and LATIXY fields'
       call check_ret(nf_inq_varid (ncid, 'LONGXY', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%longxy), subname)
       call check_ret(nf_inq_varid (ncid, 'LATIXY', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%latixy), subname)
    endif

    ier = nf_inq_varid (ncid, 'NUMLON', varid)
    if (ier == NF_NOERR) then
       write(6,*) trim(subname),' check NUMLON for regular grid'
       numlonset = .true.
       call check_ret(nf_inq_varid (ncid, 'NUMLON', varid), subname)
       call check_ret(nf_get_var_int (ncid, varid, domain%numlon), subname)
       if (minval(domain%numlon) /= nlon .or. maxval(domain%numlon) /= nlon) then
          write(6,*) trim(subname),' ERROR not regular grid, stop', minval(domain%numlon),maxval(domain%numlon)
          stop
       endif
    endif

    ier = nf_inq_varid (ncid, 'EDGEN', varid)
    if (ier == NF_NOERR) then
       if (edgeNESWset) write(6,*) trim(subname),' WARNING, overwriting edges'
       edgeNESWset = .true.
       write(6,*) trim(subname),' read EDGE[NESW]'
       call check_ret(nf_inq_varid (ncid, 'EDGEN', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%edgen), subname)

       call check_ret(nf_inq_varid (ncid, 'EDGEE', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%edgee), subname)

       call check_ret(nf_inq_varid (ncid, 'EDGES', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%edges), subname)

       call check_ret(nf_inq_varid (ncid, 'EDGEW', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%edgew), subname)
       !--- check whether edges read are acceptable ---
       if ((abs(domain%edgen)+abs(domain%edges)+ &
            abs(domain%edgee)+abs(domain%edgew)) > 1.0e6 .or. &
            abs(domain%edgen-domain%edges) > 1.0e3 .or. &
            abs(domain%edgee-domain%edgew) > 1.0e3) then 
          edgeNESWset = .false.
       endif
    endif

    ier = nf_inq_varid (ncid, 'xv', varid)
    if (ier == NF_NOERR) then
       if (llneswset) write(6,*) trim(subname),' WARNING, overwriting lat[ns],lon[we]'
       llneswset = .true.
       write(6,*) trim(subname),' read xv and yv'
       allocate(xv(4,nlon,nlat),yv(4,nlon,nlat))
       !--- 4 pts are sw,se,ne,nw, careful of wraparound on w side
       call check_ret(nf_inq_varid (ncid, 'xv', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, xv), subname)
       call check_ret(nf_inq_varid (ncid, 'yv', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, yv), subname)
       do j = 1, nlat
       do i = 1, nlon
          if (xv(1,i,j) /= xv(4,i,j) .or. xv(2,i,j) /= xv(3,i,j) .or. &
              yv(1,i,j) /= yv(2,i,j) .or. yv(3,i,j) /= yv(4,i,j)) then
                 write(6,*) trim(subname),' ERROR not regular grid, stop', xv(:,i,j),yv(:,i,j)
                 stop
          endif
          domain%latn(i,j) = yv(3,i,j)
          domain%lone(i,j) = xv(3,i,j)
          domain%lats(i,j) = yv(1,i,j)
          domain%lonw(i,j) = xv(1,i,j)
       enddo
       enddo
       domain%edgee = maxval(domain%lone)
       domain%edgew = minval(domain%lonw)
       deallocate(xv,yv)
    endif

    ier = nf_inq_varid (ncid, 'LATN', varid)
    if (ier == NF_NOERR) then
       if (llneswset) write(6,*) trim(subname),' WARNING, overwriting lat[ns],lon[we]'
       llneswset = .true.
       write(6,*) trim(subname),' read LAT[NS],LON[WE]'
       call check_ret(nf_inq_varid (ncid, 'LATN', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%latn), subname)
       call check_ret(nf_inq_varid (ncid, 'LONE', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lone), subname)
       call check_ret(nf_inq_varid (ncid, 'LATS', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lats), subname)
       call check_ret(nf_inq_varid (ncid, 'LONW', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%lonw), subname)
    endif

    ier = nf_inq_varid (ncid, 'frac', varid)
    if (ier == NF_NOERR) then
       if (landfracset) write(6,*) trim(subname),' WARNING, overwriting frac'
       landfracset = .true.
       write(6,*) trim(subname),' read frac'
       call check_ret(nf_inq_varid (ncid, 'frac', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%frac), subname)
    endif

    ier = nf_inq_varid (ncid, 'LANDFRAC', varid)
    if (ier == NF_NOERR) then
       if (landfracset) write(6,*) trim(subname),' WARNING, overwriting frac'
       landfracset = .true.
       write(6,*) trim(subname),' read LANDFRAC'
       call check_ret(nf_inq_varid (ncid, 'LANDFRAC', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%frac), subname)
    endif

    ier = nf_inq_varid (ncid, 'landfract', varid)
    if (ier == NF_NOERR) then
       if (landfracset) write(6,*) trim(subname),' WARNING, overwriting frac'
       landfracset = .true.
       write(6,*) trim(subname),' read landfract'
       call check_ret(nf_inq_varid (ncid, 'landfract', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%frac), subname)
    endif

    ier = nf_inq_varid (ncid, 'mask', varid)
    if (ier == NF_NOERR) then
       if (maskset) write(6,*) trim(subname),' WARNING, overwriting mask'
       maskset = .true.
       write(6,*) trim(subname),' read mask'
       call check_ret(nf_inq_varid (ncid, 'mask', varid), subname)
       call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)
    endif

    ier = nf_inq_varid (ncid, 'LANDMASK', varid)
    if (ier == NF_NOERR) then
       if (maskset) write(6,*) trim(subname),' WARNING, overwriting mask'
       maskset = .true.
       write(6,*) trim(subname),' read LANDMASK'
       call check_ret(nf_inq_varid (ncid, 'LANDMASK', varid), subname)
       call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)
    endif

    ier = nf_inq_varid (ncid, 'area', varid)
    if (ier == NF_NOERR) then
       if (areaset) write(6,*) trim(subname),' WARNING, overwriting area'
       areaset = .true.
       write(6,*) trim(subname),' read area'
       call check_ret(nf_inq_varid (ncid, 'area', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%area), subname)
    endif

    ier = nf_inq_varid (ncid, 'AREA', varid)
    if (ier == NF_NOERR) then
       if (areaset) write(6,*) trim(subname),' WARNING, overwriting area'
       areaset = .true.
       write(6,*) trim(subname),' read AREA'
       call check_ret(nf_inq_varid (ncid, 'AREA', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%area), subname)
    endif

    ier = nf_inq_varid (ncid, 'PHIS', varid)
    if (ier == NF_NOERR) then
       if (toposet) write(6,*) trim(subname),' WARNING, overwriting topo'
       toposet = .true.
       write(6,*) trim(subname),' read TOPO'
       call check_ret(nf_inq_varid (ncid, 'PHIS', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%topo), subname)
       domain%topo = domain%topo/SHR_CONST_G
    endif

    ier = nf_inq_varid (ncid, 'TOPO', varid)
    if (ier == NF_NOERR) then
       if (toposet) write(6,*) trim(subname),' WARNING, overwriting topo'
       toposet = .true.
       write(6,*) trim(subname),' read TOPO'
       call check_ret(nf_inq_varid (ncid, 'TOPO', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%topo), subname)
    endif

    ier = nf_inq_varid (ncid, 'htopo', varid)
    if (ier == NF_NOERR) then
       if (toposet) write(6,*) trim(subname),' WARNING, overwriting topo'
       toposet = .true.
       write(6,*) trim(subname),' read htopo'
       call check_ret(nf_inq_varid (ncid, 'htopo', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%topo), subname)
    endif

    if (area_units == 1) then
       domain%area = domain%area * re * re
       area_units = 0
    endif

    call check_ret(nf_close(ncid), subname)

!----------------------
    if (etype) then
    if (mksrf_fcamfile /= '') then
    fnamel = mksrf_fcamfile
    call getfil (fnamel, locfn, 0)
    call check_ret(nf_open(locfn, 0, ncid), subname)

    ier = nf_inq_varid (ncid, 'lon', varid)
    if (ier == NF_NOERR) then
       if (lonlatset) write(6,*) trim(subname),' WARNING, overwriting lat,lon'
       lonlatset = .true.
       ier = nf_inq_varndims(ncid,varid,ndims)
       if (ndims == 1) then
          write(6,*) trim(subname),' read lon and lat 1d fields'
          allocate(lon1d(nlon),lat1d(nlat))
          call check_ret(nf_inq_varid (ncid, 'lon', varid), subname)
          call check_ret(nf_get_var_double (ncid, varid, lon1d), subname)
          call check_ret(nf_inq_varid (ncid, 'lat', varid), subname)
          call check_ret(nf_get_var_double (ncid, varid, lat1d), subname)
          do j = 1, nlat
          do i = 1, nlon
             domain%longxy(i,j) = lon1d(i)
             domain%latixy(i,j) = lat1d(j)
          enddo
          enddo
          deallocate(lon1d,lat1d)
       elseif (ndims == 2) then
          write(6,*) trim(subname),' read lon and lat 2d fields'
          call check_ret(nf_inq_varid (ncid, 'lat', varid), subname)
          call check_ret(nf_get_var_double (ncid, varid, domain%latixy), subname)
          call check_ret(nf_inq_varid (ncid, 'lon', varid), subname)
          call check_ret(nf_get_var_double (ncid, varid, domain%longxy), subname)
       else
          write(6,*) trim(subname),'ERROR: lon and lat illegal dim ',ndims
          stop
       endif
    endif

    call check_ret(nf_close(ncid), subname)

    endif
    endif
!----------------------
!----------------------
    if (etype) then
    if (mksrf_fcamtopo /= '') then
    fnamel = mksrf_fcamtopo
    call getfil (fnamel, locfn, 0)
    call check_ret(nf_open(locfn, 0, ncid), subname)

    ier = nf_inq_varid (ncid, 'PHIS', varid)
    if (ier == NF_NOERR) then
       if (toposet) write(6,*) trim(subname),' WARNING, overwriting topo'
       toposet = .true.
       write(6,*) trim(subname),' read TOPO'
       call check_ret(nf_inq_varid (ncid, 'PHIS', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%topo), subname)
       domain%topo = domain%topo/SHR_CONST_G
    endif

    ier = nf_inq_varid (ncid, 'TOPO', varid)
    if (ier == NF_NOERR) then
       if (toposet) write(6,*) trim(subname),' WARNING, overwriting topo'
       toposet = .true.
       write(6,*) trim(subname),' read TOPO'
       call check_ret(nf_inq_varid (ncid, 'TOPO', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%topo), subname)
    endif

    call check_ret(nf_close(ncid), subname)

    endif
    endif
!----------------------

    if (.not.lonlatset) then
       write(6,*) trim(subname),' ERROR: long and lati not set '
       stop
    endif

    if (.not.numlonset) then
       domain%numlon(:) = domain%ni
    endif

    if (.not.maskset.and.landfracset) then
       maskset = .true.
       where (domain%frac(:,:) < flandmin)
          domain%mask(:,:) = 0     !ocean
       elsewhere
          domain%mask(:,:) = 1     !land
       endwhere
    endif

    if (.not.landfracset.and.maskset) then
       landfracset = .true.
       where (domain%mask(:,:) == 0)
          domain%frac(:,:) = 0._r8     !ocean
       elsewhere
          domain%frac(:,:) = 1._r8     !land
       endwhere
    endif

    if (.not.llneswset) then
       llneswset = .true.
       if (edgeneswset) then
          write(6,*) trim(subname),' compute lat[ns],lon[we] from edge[nesw]'
          call celledge (domain, domain%edgen,domain%edgee,domain%edges,domain%edgew)
       else
          write(6,*) trim(subname),' compute lat[ns],lon[we]'
          call celledge(domain)
       endif
    endif

! check n/s/e/w/ consistent with center
    write(6,*) trim(subname),' check nesw consistent with center'
    call shr_sys_flush(6)
    do j = 1, nlat
    do i = 1, nlon
       if (domain%lone(i,j) < domain%longxy(i,j))  then
          domain%lone(i,j) = domain%lone(i,j) + 360.0_r8
          domain%edgee = max(domain%edgee,domain%lone(i,j))
       endif
       if (domain%lonw(i,j) > domain%longxy(i,j))  then
          domain%lonw(i,j) = domain%lonw(i,j) - 360.0_r8
          domain%edgew = min(domain%edgew,domain%lonw(i,j))
       endif
    enddo
    enddo

    if (.not.edgeneswset) then
       write(6,*) trim(subname),' set edges '
       call shr_sys_flush(6)
       edgeneswset = .true.
       domain%edgen = maxval(domain%latn)
       domain%edgee = maxval(domain%lone)
       domain%edges = minval(domain%lats)
       domain%edgew = minval(domain%lonw)
    endif

    if (.not.areaset .or. .not.area_valid) then
       areaset = .true.
       if (edgeneswset) then
          write(6,*) trim(subname),' compute cellarea with edge[nesw]'
          call shr_sys_flush(6)
          call cellarea (domain, domain%edgen,domain%edgee,domain%edges,domain%edgew)
       else
          write(6,*) trim(subname),' compute cellarea'
          call shr_sys_flush(6)
          call cellarea (domain)
       endif
    endif

    write(6,*) trim(subname),' done'
    call shr_sys_flush(6)

!    write(6,*) ' '
!    write(6,*) trim(subname),':'
!    call domain_check(domain)

  end subroutine read_domain
!----------------------------------------------------------------------------

  subroutine mkfile(lsmlon, lsmlat, fname, finfo, itype)

    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_sys_mod , only : shr_sys_getenv
    use fileutils   , only : get_filename
    use mkvarctl
    use ncdio

    implicit none
    integer, intent(in) :: lsmlon, lsmlat
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) :: finfo
    integer, intent(in), optional :: itype

    integer :: ncid
    integer :: j                    ! index
    integer :: pftsize              ! size of lsmpft dimension
    integer :: dimid                ! temporary
    integer :: values(8)            ! temporary
    character(len=256) :: str       ! global attribute string
    character(len=256) :: name      ! name of attribute
    character(len=256) :: unit      ! units of attribute
    character(len= 18) :: datetime  ! temporary
    character(len=  8) :: date      ! temporary
    character(len= 10) :: time      ! temporary
    character(len=  5) :: zone      ! temporary
    integer            :: ier       ! error status
    integer            :: omode     ! netCDF output mode
    character(len=32)  :: subname = 'mkfile'  ! subroutine name
    integer            :: type      ! 1=grid, 2=frac
!-----------------------------------------------------------------------

    type = 1
    if (present(itype)) then
       type = itype
    endif

    call check_ret(nf_create(trim(fname), nf_clobber, ncid), subname)
    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Define dimensions.

    call check_ret(nf_def_dim (ncid, 'lsmlon'    , lsmlon      , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'lsmlat'    , lsmlat      , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'nchar'  , 128         , dimid), subname)

    ! Create global attributes.

    str = 'NCAR-CSM'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Conventions', len_trim(str), trim(str)), subname)

    call date_and_time (date, time, zone, values)
    datetime(1:8) =        date(5:6) // '-' // date(7:8) // '-' // date(3:4)
    datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '
    str = 'created on: ' // datetime
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'History_Log', len_trim(str), trim(str)), subname)

    call shr_sys_getenv ('LOGNAME', str, ier)
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Logname', len_trim(str), trim(str)), subname)

    call shr_sys_getenv ('HOST', str, ier)
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Host', len_trim(str), trim(str)), subname)

    str = mksrf_fcamfile
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'fcamfile', len_trim(str), trim(str)), subname)

    str = mksrf_fccsmdom
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'fccsmdom', len_trim(str), trim(str)), subname)

    str = mksrf_fcamtopo
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'fcamtopo', len_trim(str), trim(str)), subname)

    str = finfo
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Other_Info', len_trim(str), trim(str)), subname)

    if (itype == 3) then
      str = mksrf_frawtopo
      call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
           'Topo_Input_Filename', len_trim(str), trim(str)), subname)
    endif

    str = 'Community Land Model: CLM3'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Source', len_trim(str), trim(str)), subname)

    str = &
'$HeadURL$'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Version', len_trim(str), trim(str)), subname)

    str = '$Id$'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Revision_Id', len_trim(str), trim(str)), subname)

    ! ----------------------------------------------------------------------
    ! Define variables
    ! ----------------------------------------------------------------------

    call ncd_defvar(ncid=ncid, varname='NUMLON', xtype=nf_int, &
         dim1name='lsmlat', &
         long_name='number of grid cells at each latitude', units='unitless')

    call ncd_defvar(ncid=ncid, varname='LONGXY', xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='longitude', units='degrees east')

    call ncd_defvar(ncid=ncid, varname='LATIXY', xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='latitude', units='degrees north')

    if (type == 1) then

       call ncd_defvar(ncid=ncid, varname='EDGEN', xtype=nf_double, &
         long_name='northern edge of surface grid', units='degrees north')
    
       call ncd_defvar(ncid=ncid, varname='EDGEE', xtype=nf_double, &
         long_name='eastern edge of surface grid', units='degrees east')
    
       call ncd_defvar(ncid=ncid, varname='EDGES', xtype=nf_double, &
         long_name='southern edge of surface grid', units='degrees north')
    
       call ncd_defvar(ncid=ncid, varname='EDGEW', xtype=nf_double, &
         long_name='western edge of surface grid', units='degrees east')

       call ncd_defvar(ncid=ncid, varname='LATN' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='latitude of north edge', units='degrees north')

       call ncd_defvar(ncid=ncid, varname='LONE' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='longitude of east edge', units='degrees east')

       call ncd_defvar(ncid=ncid, varname='LATS' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='latitude of south edge', units='degrees north')

       call ncd_defvar(ncid=ncid, varname='LONW' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='longitude of west edge', units='degrees east')

       call ncd_defvar(ncid=ncid, varname='AREA' , xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='area', units='km^2')

    elseif (type == 2) then

       call ncd_defvar(ncid=ncid, varname='LANDMASK', xtype=nf_int, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='land/ocean mask', units='0=ocean and 1=land')

       call ncd_defvar(ncid=ncid, varname='LANDFRAC', xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='land fraction', units='unitless')


    elseif (type == 3) then

       call ncd_defvar(ncid=ncid, varname='TOPO', xtype=nf_double, &
         dim1name='lsmlon', dim2name='lsmlat', &
         long_name='topography height', units='m')

    else

       write(6,*) 'ERROR: itype value invalid ',type
       stop

    endif

    ! End of define mode

    call check_ret(nf_enddef(ncid), subname)
    call check_ret(nf_close(ncid), subname)

  end subroutine mkfile

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_domain
!
! !INTERFACE:
  subroutine write_domain(domain,fname,itype)
!
! !USES:
!
! !DESCRIPTION:
! Write a domain to netcdf

! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout)       :: domain
    character(len=*) ,intent(in)          :: fname
    integer          ,intent(in),optional :: itype
!
! !REVISION HISTORY:
! Author: T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ncid                           !netCDF file id
    integer  :: omode                          !netCDF output mode
    character(len= 32) :: subname = 'write_domain'
    integer  :: type   ! 1=grid, 2=frac
!-----------------------------------------------------------------

     type = 1
     if (present(itype)) then
        type = itype
     endif

!    write(6,*) ' '
!    write(6,*) trim(subname),':'
!    call domain_check(domain)

    call check_ret(nf_open(trim(fname), nf_write, ncid), subname)
    ! File will be in define mode. Set fill mode to "no fill" to optimize performance

    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Write domain fields 

    call ncd_ioglobal(varname='NUMLON'  , data=domain%numlon, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LONGXY'  , data=domain%longxy, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LATIXY'  , data=domain%latixy, ncid=ncid, flag='write')

    if (type == 1) then

       call ncd_ioglobal(varname='EDGEN'   , data=domain%edgen, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='EDGEE'   , data=domain%edgee, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='EDGES'   , data=domain%edges, ncid=ncid, flag='write')
       call ncd_ioglobal(varname='EDGEW'   , data=domain%edgew, ncid=ncid, flag='write')

       call ncd_ioglobal(varname='LATN'    , data=domain%latn , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LONE'    , data=domain%lone , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LATS'    , data=domain%lats , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LONW'    , data=domain%lonw , ncid=ncid, flag='write')

       call ncd_ioglobal(varname='AREA'    , data=domain%area  , ncid=ncid, flag='write')

    elseif (type == 2) then

       call ncd_ioglobal(varname='LANDMASK', data=domain%mask  , ncid=ncid, flag='write')
       call ncd_ioglobal(varname='LANDFRAC', data=domain%frac  , ncid=ncid, flag='write')

    elseif (type == 3) then

       call ncd_ioglobal(varname='TOPO'    , data=domain%topo  , ncid=ncid, flag='write')

    else

       write(6,*) 'ERROR: itype value invalid ',type
       stop

    endif

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    ! Close grid data dataset

    call check_ret(nf_close(ncid), subname)

  end subroutine write_domain

!----------------------------------------------------------------------------
end module creategridMod
