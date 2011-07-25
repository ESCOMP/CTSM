module creategridMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkgridMod
!
! !DESCRIPTION:
! Routines to create land model grid
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use mkvarctl   
  use domainMod   , only : domain_init, domain_type, domain_check
  use areaMod
  use ncdio
!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: read_domain   ! read domain from netcdf file
  public :: write_domain  ! write domain to netcdf file

! !PRIVATE MEMBER FUNCTIONS:
  real(r8) :: flandmin = 0.001            !minimum land frac for land cell
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_domain
!
! !INTERFACE:
  subroutine read_domain(domain,fname,readmask)
!
! !DESCRIPTION:
! Read a grid file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain
    character(len=*) ,intent(in)    :: fname
    logical,optional, intent(in)    :: readmask ! true => read mask instead of landmask for urban parameters
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  include 'netcdf.inc'
    integer :: nlon,nlat                       !size
    real(r8), allocatable :: lon1d(:)          !local array for 1d lon
    real(r8), allocatable :: lat1d(:)          !local array for 1d lat
    real(r8), allocatable :: xv(:,:,:)         !local array for corner lons
    real(r8), allocatable :: yv(:,:,:)         !local array for corner lats
    integer :: i,j                             !indexes
    integer :: ncid                            !netCDF file id
    integer :: dimid                           !netCDF dimension id
    integer :: varid                           !netCDF variable id
    logical :: dimset                          !local ni,nj
    logical :: lonlatset                       !local lon(:,:), lat(:,:)
    logical :: edgeNESWset                     !local EDGE[NESW]
    logical :: llneswset                       !local lat[ns],lon[we]
    logical :: areaset                         !local area
    logical :: landfracset                     !local landfrac
    logical :: maskset                         !local mask
    logical :: numlonset                       !local numlon
    integer :: ndims                           !number of dims for variable
    integer :: ier
    logical :: lreadmask                       !local readmask
    character(len= 32) :: subname = 'read_domain'
!-----------------------------------------------------------------

    write(6,*) ' ' 

    dimset      = .false. 
    lonlatset   = .false. 
    edgeNESWset = .false. 
    llneswset   = .false. 
    areaset     = .false. 
    landfracset = .false. 
    maskset     = .false. 
    numlonset   = .false. 
    lreadmask   = .false.

    if (present(readmask)) then
       lreadmask = readmask
    end if

    ! Read domain file and compute stuff as needed

    call check_ret(nf_open(fname, 0, ncid), subname)

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
          domain%lats(i,j) = yv(2,i,j)
          domain%lonw(i,j) = xv(2,i,j)
       enddo
       enddo
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

    if (lreadmask) then
       ier = nf_inq_varid (ncid, 'mask', varid)
       if (ier == NF_NOERR) then
          if (maskset) write(6,*) trim(subname),' WARNING, overwriting mask'
          maskset = .true.
          write(6,*) trim(subname),' read mask'
          call check_ret(nf_inq_varid (ncid, 'mask', varid), subname)
          call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)
       endif
    else
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
    end if

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

    call check_ret(nf_close(ncid), subname)

    if (.not.lonlatset) then
       write(6,*) trim(subname),' ERROR: long and lati not set '
       stop
    endif

    if (.not.numlonset) then
       domain%numlon(:) = domain%ni
    endif

    if (.not.maskset.and.landfracset) then
       maskset = .true.
       where (domain%frac < flandmin)
          domain%mask = 0     !ocean
       elsewhere
          domain%mask = 1     !land
       endwhere
    endif

    if (.not.landfracset.and.maskset) then
       landfracset = .true.
       do j = 1, nlat
       do i = 1, nlon
          if ( domain%mask(i,j) == 0 )then
             domain%frac(i,j) = 0._r8     !ocean
          else
             domain%frac(i,j) = 1._r8     !land
          end if
       end do
       end do
    endif

    if (.not.llneswset) then
       if (edgeneswset) then
          write(6,*) trim(subname),' compute lat[ns],lon[we] from edge[nesw]'
          call celledge (domain, domain%edgen,domain%edgee,domain%edges,domain%edgew)
       else
          write(6,*) trim(subname),' compute lat[ns],lon[we]'
          call celledge(domain)
       endif
    endif

    if (.not.areaset) then
       if (edgeneswset) then
          write(6,*) trim(subname),' compute cellarea with edge[nesw]'
          call cellarea (domain, domain%edgen,domain%edgee,domain%edges,domain%edgew)
       else
          write(6,*) trim(subname),' compute cellarea'
          call cellarea (domain)
       endif
    endif

    if (.not.edgeneswset) then
        domain%edgen = maxval(domain%latn)
        domain%edgee = maxval(domain%lone)
        domain%edges = minval(domain%lats)
        domain%edgew = minval(domain%lonw)
    endif

!    write(6,*) ' '
!    write(6,*) trim(subname),':'
!    call domain_check(domain)

  end subroutine read_domain

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_domain
!
! !INTERFACE:
  subroutine write_domain(domain,fname)
!
! !USES:
!
! !DESCRIPTION:
! Write a domain to netcdf

! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain
    character(len=*) ,intent(in)    :: fname
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: ncid                           !netCDF file id
    integer  :: omode                          !netCDF output mode
    character(len= 32) :: subname = 'write_domain'
!-----------------------------------------------------------------

!    write(6,*) ' '
!    write(6,*) trim(subname),':'
!    call domain_check(domain)

    call check_ret(nf_open(trim(fname), nf_write, ncid), subname)
    ! File will be in define mode. Set fill mode to "no fill" to optimize performance

    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Write domain fields 

    call ncd_ioglobal(varname='EDGEN'   , data=domain%edgen, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EDGEE'   , data=domain%edgee, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EDGES'   , data=domain%edges, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='EDGEW'   , data=domain%edgew, ncid=ncid, flag='write')

    call ncd_ioglobal(varname='LATN'    , data=domain%latn , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LONE'    , data=domain%lone , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LATS'    , data=domain%lats , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LONW'    , data=domain%lonw , ncid=ncid, flag='write')

    call ncd_ioglobal(varname='NUMLON'  , data=domain%numlon, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='AREA'    , data=domain%area  , ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LONGXY'  , data=domain%longxy, ncid=ncid, flag='write')
    call ncd_ioglobal(varname='LATIXY'  , data=domain%latixy, ncid=ncid, flag='write')
!    call ncd_ioglobal(varname='LANDMASK', data=domain%mask  , ncid=ncid, flag='write')
!    call ncd_ioglobal(varname='LANDFRAC', data=domain%frac  , ncid=ncid, flag='write')

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    ! Close grid data dataset

    call check_ret(nf_close(ncid), subname)

  end subroutine write_domain

!----------------------------------------------------------------------------
end module creategridMod
