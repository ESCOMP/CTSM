!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: create_domain
!
! !INTERFACE:
   PROGRAM create_domain
!
! !DESCRIPTION:
! "off-line" tool to create a data model domain file from CLM land model input files.
!
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  implicit none
  include 'netcdf.inc'
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
!
! !LOCAL VARIABLES:
!EOP
  integer  :: ni,nj                         ! grid resolution
  integer  :: nv=4                          ! number of verticies
  integer  :: ni_grid,nj_grid               ! grid resolution
  integer  :: i,j,k                         ! indices
  integer  :: ier                           ! error status
  integer  :: omode                         ! netCDF output mode
  integer  :: ncid                          ! temporary
  integer  :: dimid(0:2)                    ! temporary
  integer  :: varid                         ! temporary
  character(len=SHR_KIND_CL) :: f_fracdata  ! clm fracdata file
  character(len=SHR_KIND_CL) :: f_griddata  ! clm griddata file
  character(len=SHR_KIND_CL) :: f_domain    ! output ocn domain file
  character(len=SHR_KIND_CL) :: str         ! global attribute string
  character(len=SHR_KIND_CL) :: name        ! name of attribute
  character(len=SHR_KIND_CL) :: unit        ! units of attribute
  character(len=  5) :: dtype               ! Type of data model

  real(r8), allocatable :: lndfrac(:,:)     ! Land fraction
  real(r8), allocatable :: ocnfrac(:,:)     ! Ocean fraction
  integer , allocatable :: ocnmask(:,:)     ! Ocean mask
  integer , allocatable :: lndmask(:,:)     ! Land mask
  real(r8), allocatable :: xc(:,:)          ! Longitudes
  real(r8), allocatable :: yc(:,:)          ! Latitudes
  real(r8), allocatable :: xv2(:,:)         ! Longitudes vertice in northern east
  real(r8), allocatable :: xv4(:,:)         ! Longitudes vertice in southern west
  real(r8), allocatable :: yv1(:,:)         ! Latitudes vertice in northern west
  real(r8), allocatable :: yv3(:,:)         ! Latitudes vertice in southern west
  real(r8), allocatable :: xv(:,:,:)        ! Longitude vertices
  real(r8), allocatable :: yv(:,:,:)        ! Latitude vertices
  real(r8), allocatable :: area(:,:)        ! Grid areas
  real(r8), parameter :: re = SHR_CONST_REARTH * 0.001_r8

  integer :: filename_position              ! Character position for last character of directory
  integer :: filename_length                ! Length of full pathname
  integer :: indx                           ! Index to mask and mask_values arrays
  integer, parameter :: ntype     = 2
  integer, parameter :: indx_docn = 1       ! MUST go with arrays set below
  integer, parameter :: indx_datm = 2       ! MUST go with arrays set below
  character(len=SHR_KIND_CL), parameter :: mask(ntype) = (/         &
            'ocean domain mask', 'land domain mask ' /)
  character(len=SHR_KIND_CL), parameter :: mask_values(ntype) = (/       &
         '1=ocean and 0=land, 0 indicates that cell is not active', &
         '0=ocean and 1=land, 0 indicates that cell is not active'  &
  /)

  character(len=32) :: subname = 'create_domain'
  
  namelist /domain_nl/ f_fracdata, f_griddata, f_domain, dtype

  !-----------------------------------------------------------------------
  ! Determine land fracdata, land griddata and output domain file
  !-----------------------------------------------------------------------

  dtype = "docn"
  write(6,*) 'Attempting to initialize control settings .....'
  read(5, domain_nl, iostat=ier)
  if (ier /= 0) then
     write(6,*)'error: namelist input resulted in error code ',ier
     stop
  endif
  if (      index(dtype, "docn") == 1 )then
     indx = indx_docn
     write(6,*)'Creating files for use with docn'
  else if ( index(dtype, "datm") == 1 )then
     indx = indx_datm
     write(6,*)'Creating files for use with datm for land model'
  else
     write(6,*)'error: dtype MUST be either docn or datm = ', trim(dtype)
     stop
  end if
  
  !-------------------------------------------------------------------
  ! Read clm fractional dataset ( for landfrac )
  !-------------------------------------------------------------------
  
  write(6,*) trim(subname),': reading in ',trim(f_fracdata)
  call check_ret(nf_open(f_fracdata, 0, ncid), subname)
  
  write(6,*) trim(subname),': reading ni and nj dims'
  call check_ret(nf_inq_dimid  (ncid, 'lsmlon', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, ni), subname)
  call check_ret(nf_inq_dimid  (ncid, 'lsmlat', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nj), subname)
  
  allocate(lndfrac(ni,nj))

  write(6,*) trim(subname),': reading LANDFRAC'
  call check_ret(nf_inq_varid (ncid, 'LANDFRAC', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, lndfrac), subname)
  
  write(6,*) trim(subname),': closing ',trim(f_fracdata)
  call check_ret(nf_close(ncid), subname)
    
  !-------------------------------------------------------------------
  ! Read clm grid dataset (for lats, longs, area)
  !-------------------------------------------------------------------

  write(6,*)
  write(6,*) trim(subname),': reading in ',trim(f_griddata)
  call check_ret(nf_open(f_griddata, 0, ncid), subname)
  
  write(6,*) trim(subname),': reading ni and nj dims'
  call check_ret(nf_inq_dimid  (ncid, 'lsmlon', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, ni_grid), subname)
  call check_ret(nf_inq_dimid  (ncid, 'lsmlat', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nj_grid), subname)
  
  if (ni_grid /= ni) then
     write(6,*)'ni_grid= ',ni_grid,' must be same as ni_frac= ',ni 
     stop
  end if

  if (nj_grid /= nj) then
     write(6,*)'nj_grid= ',ni_grid,' must be same as nj_frac= ',ni 
     stop
  end if

  allocate(xc(ni,nj), yc(ni,nj), xv(nv,ni,nj), yv(nv,ni,nj), area(ni,nj))
  allocate(yv1(ni,nj),xv2(ni,nj),yv3(ni,nj),xv4(ni,nj))

  write(6,*) trim(subname),': reading LONGXY'
  call check_ret(nf_inq_varid (ncid, 'LONGXY', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, xc), subname)
  
  write(6,*) trim(subname),': reading LATIXY'
  call check_ret(nf_inq_varid (ncid, 'LATIXY', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, yc), subname)

  write(6,*) trim(subname),': reading LATN'
  call check_ret(nf_inq_varid (ncid, 'LATN', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, yv3), subname)

  write(6,*) trim(subname),': reading LATS'
  call check_ret(nf_inq_varid (ncid, 'LATS', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, yv1), subname)

  write(6,*) trim(subname),': reading LONE'
  call check_ret(nf_inq_varid (ncid, 'LONE', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, xv2), subname)

  write(6,*) trim(subname),': reading LONW'
  call check_ret(nf_inq_varid (ncid, 'LONW', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, xv4), subname)

  xv(1,:,:) = xv4(:,:)
  xv(2,:,:) = xv2(:,:)
  xv(3,:,:) = xv2(:,:)
  xv(4,:,:) = xv4(:,:)

  ! Nix any negative longitudes
  do j = 1, nj
  do i = 1, ni
  do k = 1, nv
     if ( xv(k,i,j) < 0.0_r8 ) xv(k,i,j) = 360.0_r8 + xv(k,i,j)
  end do
  end do
  end do

  yv(1,:,:) = yv1(:,:)
  yv(2,:,:) = yv1(:,:)
  yv(3,:,:) = yv3(:,:)
  yv(4,:,:) = yv3(:,:)
  
  write(6,*) trim(subname),': reading AREA'
  call check_ret(nf_inq_varid (ncid, 'AREA', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, area), subname)
  
  write(6,*)' closing ',trim(f_griddata)
  call check_ret(nf_close(ncid), subname)

  !-------------------------------------------------------------------
  ! Create new domain file 
  !-------------------------------------------------------------------

  write(6,*)
  write(6,*)trim(subname),': creating f_domain = ',trim(f_domain)
  call check_ret(nf_create(trim(f_domain), nf_clobber, ncid), subname)
  call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

  call check_ret(nf_def_dim (ncid, 'nv', nv, dimid(0)), subname)
  call check_ret(nf_def_dim (ncid, 'ni', ni, dimid(1)), subname)
  call check_ret(nf_def_dim (ncid, 'nj', nj, dimid(2)), subname)

  call addglobal( ncid )
  
  filename_position = index( f_griddata, '/', back=.true. )
  filename_length  = len_trim(f_griddata)
  if ( filename_position == 0 )then
     str = f_griddata
  else
     str = f_griddata(filename_position+1:filename_length)
  end if
  write(6,*)'grid data is ',trim(str)
  call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
                 'Land_Grid_Dataset', len_trim(str), trim(str)), subname)
  
  filename_position = index( f_fracdata, '/', back=.true. )
  filename_length   = len_trim(f_fracdata)
  if ( filename_position == 0 )then
     str = f_fracdata
  else
     str = f_fracdata(filename_position+1:filename_length)
  end if
  write(6,*)'fractional data is ',trim(str)
  call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
                 'Land_Fraction_Dataset', len_trim(str), trim(str)), subname)
  
  name = 'xc'
  call check_ret(nf_def_var(ncid, trim(name), nf_double, 2, dimid(1:2), varid), subname)
  name = 'longitude of grid cell center'
  call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(name), trim(name)), subname)
  name = 'degrees_east'
  call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(name), trim(name)), subname)
  name = 'xv'
  call check_ret(nf_put_att_text(ncid, varid, 'bounds', len_trim(name), trim(name)), subname)

  name = 'yc'
  call check_ret(nf_def_var(ncid, trim(name), nf_double, 2, dimid(1:2), varid), subname)
  name = 'latitude of grid cell center'
  call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(name), trim(name)), subname)
  name = 'degrees_north'
  call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(name), trim(name)), subname)
  name = 'yv'
  call check_ret(nf_put_att_text(ncid, varid, 'bounds', len_trim(name), trim(name)), subname)

  name = 'xv'
  call check_ret(nf_def_var(ncid, trim(name), nf_double, 3, dimid(0:2), varid), subname)
  name = 'longitude of grid cell vertices'
  call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(name), trim(name)), subname)
  name = 'degrees_east'
  call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(name), trim(name)), subname)

  name = 'yv'
  call check_ret(nf_def_var(ncid, trim(name), nf_double, 3, dimid(0:2), varid), subname)
  name = 'latitude of grid cell vertices'
  call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(name), trim(name)), subname)
  name = 'degrees_north'
  call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(name), trim(name)), subname)

  name = 'mask'
  call check_ret(nf_def_var(ncid, trim(name), nf_int, 2, dimid(1:2), varid), subname)
  call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(mask(indx)), trim(mask(indx))), subname)
  name = 'xc yc'
  call check_ret(nf_put_att_text(ncid, varid, 'coordinate', len_trim(name), trim(name)), subname)
  name = 'unitless'
  call check_ret(nf_put_att_text(ncid, varid, 'note', len_trim(name), trim(name)), subname)
  call check_ret(nf_put_att_text(ncid, varid, 'comment', len_trim(mask_values(indx)), &
                                 trim(mask_values(indx))), subname)

  name = 'frac'
  call check_ret(nf_def_var(ncid, trim(name), nf_double, 2, dimid(1:2), varid), subname)
  name = 'fraction of grid cell that is active'
  call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(name), trim(name)), subname)
  name = 'xc yc'
  call check_ret(nf_put_att_text(ncid, varid, 'coordinate', len_trim(name), trim(name)), subname)
  name = 'unitless'
  call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(name), trim(name)), subname)
  name = "error if frac> 1.0+eps or frac < 0.0-eps; eps = 0.1000000E-11"
  call check_ret(nf_put_att_text(ncid, varid, 'filter1', len_trim(name), trim(name)), subname)
  name = "limit frac to [fminval,fmaxval]; fminval= 0.1000000E-02 fmaxval=  1.000000"
  call check_ret(nf_put_att_text(ncid, varid, 'filter2', len_trim(name), trim(name)), subname)

  name = 'area'
  call check_ret(nf_def_var(ncid, trim(name), nf_double, 2, dimid(1:2), varid), subname)
  name = "area of grid cell in radians squared" ;
  call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(name), trim(name)), subname)
  name = "xc yc"
  call check_ret(nf_put_att_text(ncid, varid, 'coordinate', len_trim(name), trim(name)), subname)
  name = 'radians2'
  call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(name), trim(name)), subname)

  ! End of define mode
  
  call check_ret(nf_enddef(ncid), subname)
  
  ! Write mask/frac fields
  
  allocate(ocnmask(ni,nj))
  allocate(lndmask(ni,nj))
  allocate(ocnfrac(ni,nj))

  do j = 1,nj
     do i = 1,ni
        ocnfrac(i,j) = 1._r8 - lndfrac(i,j)
        if ( lndfrac(i,j) > 0.0_r8 ) then
           lndmask(i,j) = 1
        else
           lndmask(i,j) = 0
        end if
        if ( ocnfrac(i,j) > 0.0_r8 ) then
           ocnmask(i,j) = 1
        else
           ocnmask(i,j) = 0
        end if
        ! Input area is assumed to be in km^2 - must convert to radians^2
        area(i,j) = area(i,j) / (re*re)
     end do
  end do
  
  write(6,*)'putting out xc'
  call check_ret(nf_inq_varid(ncid, 'xc', varid), subname)
  call check_ret(nf_put_var_double(ncid, varid, xc), subname) 
  
  write(6,*)'putting out yc'
  call check_ret(nf_inq_varid(ncid, 'yc', varid), subname)
  call check_ret(nf_put_var_double(ncid, varid, yc), subname) 

  write(6,*)'putting out xv'
  call check_ret(nf_inq_varid(ncid, 'xv', varid), subname)
  call check_ret(nf_put_var_double(ncid, varid, xv), subname) 
  
  write(6,*)'putting out yv'
  call check_ret(nf_inq_varid(ncid, 'yv', varid), subname)
  call check_ret(nf_put_var_double(ncid, varid, yv), subname) 
  
  if (      indx == indx_docn )then
     write(6,*)'putting out ocean frac'
     call check_ret(nf_inq_varid(ncid, 'frac', varid), subname)
     call check_ret(nf_put_var_double(ncid, varid, ocnfrac), subname) 
  
     write(6,*)'putting out ocean mask'
     call check_ret(nf_inq_varid(ncid, 'mask', varid), subname)
     call check_ret(nf_put_var_int(ncid, varid, ocnmask), subname) 
  else if ( indx == indx_datm )then
     write(6,*)'putting out land frac'
     call check_ret(nf_inq_varid(ncid, 'frac', varid), subname)
     call check_ret(nf_put_var_double(ncid, varid, lndfrac), subname) 
  
     write(6,*)'putting out land mask'
     call check_ret(nf_inq_varid(ncid, 'mask', varid), subname)
     call check_ret(nf_put_var_int(ncid, varid, lndmask), subname) 
  else
     write(6,*)'invalid data model type = ', trim(dtype)
     stop
  end if
  
  write(6,*)'putting out area'
  call check_ret(nf_inq_varid(ncid, 'area', varid), subname)
  call check_ret(nf_put_var_double(ncid, varid, area), subname) 
  
  ! Synchronize the disk copy of a netCDF dataset with in-memory buffers
  
  call check_ret(nf_sync(ncid), subname)
  
  ! Close grid data dataset
  
  call check_ret(nf_close(ncid), subname)

  write(6,'(/,A,//)') "Successfully created domain dataset for data model"
  
end program create_domain

!=============================================================================

subroutine check_ret(ret, calling)

  ! Check return status from netcdf call

  implicit none
  include 'netcdf.inc'
  integer, intent(in) :: ret
  character(len=*) :: calling

  if (ret /= NF_NOERR) then
     write(6,*)'netcdf error from ',trim(calling),' error= ',ret 
     write(6,*) NF_STRERROR(ret)
     stop
  end if

end subroutine check_ret


