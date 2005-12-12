program convert_pft

  implicit none
  include 'netcdf.inc'

!-----------------------------------------------------------------
! make surface type netcdf file
!-----------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer, parameter :: nlon = 720        !input grid : longitude points
  integer, parameter :: nlat = 360        !input grid : latitude  points
  integer, parameter :: numpft = 16       !number of plant types

  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)        
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid
  real(r8) :: dx,dy                       !grid increments

  real(r8) :: landmask(nlon,nlat)         !fraction of land
  real(r8) :: pct_pft(nlon,nlat,0:numpft) !percent pft

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id
  integer :: dimpft_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: landmask_id                  !landmask id
  integer :: pct_pft_id                   !pct_pft id

  integer :: i,j,m                        !indices
  integer :: ndata = 1                    !input unit
  integer :: ncid                         !netCDF file id
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: dim3_id(3)                   !netCDF dimension id for 3-d variables
  integer :: status                       !status

  character(len=80) :: filei, fileo       !input,output filenames
  character(len=80) :: name,unit          !netCDF attributes

!-----------------------------------------------------------------

! Determine input and output file names

  filei = '/ptmp/slevis/lsmv2_2/input/0.5x0.5/pft-igbp.5x.5'
  fileo = '/ptmp/slevis/lsmv2_2/input/0.5x0.5/mksrf_pft.nc'

! -----------------------------------------------------------------
! Determine grid for input data 
!
! Data are 0.5 x 0.5 degree, stored in latitude bands,
! from south to north. In a given latitude band, data begin
! at the dateline and proceed eastward. So first data
! point (x(1,1)) is a box centered at 89.75S, 179.75W
!
!   89.5S  ---------------------
!          |         |         |
!          |    x    |    x    |
!          |  (1,1)  |  (2,1)  |
!          |         |         |
!   90.0S  ---------------------
!        180.W     179.5W    179.0W
! -----------------------------------------------------------------

! Define North, East, South, West edges of grid

  edge(1) =   90.
  edge(2) =  180.
  edge(3) =  -90.
  edge(4) = -180.

! Make latitudes and longitudes at center of grid cell

  dx = (edge(2)-edge(4)) / nlon
  dy = (edge(1)-edge(3)) / nlat

  do j = 1, nlat
     do i = 1, nlon
        latixy(i,j) = (edge(3)+dy/2.) + (j-1)*dy
        longxy(i,j) = (edge(4)+dx/2.) + (i-1)*dx
       end do
  end do

  lat(:) = latixy(1,:)
  lon(:) = longxy(:,1)

! -----------------------------------------------------------------
! Create netcdf file
! -----------------------------------------------------------------

! Define global attributes

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'pft_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon', nlon    , dimlon_id)
  call wrap_def_dim (ncid, 'lat', nlat    , dimlat_id)
  call wrap_def_dim (ncid, 'pft', numpft+1, dimpft_id)

! Define grid variables

  name = 'lon'
  unit = 'degrees east'
  dim1_id(1) = dimlon_id
  call wrap_def_var (ncid,'LON', nf_float, 1, dim1_id, lon_id)
  call wrap_put_att_text (ncid, lon_id, 'long_name', name)
  call wrap_put_att_text (ncid, lon_id, 'units'    , unit)

  name = 'lat'
  unit = 'degrees north'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid,'LAT', nf_float, 1, dim1_id, lat_id)
  call wrap_put_att_text (ncid, lat_id, 'long_name', name)
  call wrap_put_att_text (ncid, lat_id, 'units'    , unit)

  name = 'longitude-2d'
  unit = 'degrees east'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'western edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

! Define pft variables

  name = 'land mask'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

  name = 'percent pft'
  unit = 'unitless'
  dim3_id(1) = dimlon_id
  dim3_id(2) = dimlat_id
  dim3_id(3) = dimpft_id
  call wrap_def_var (ncid ,'PCT_PFT' ,nf_float, 3, dim3_id, pct_pft_id)
  call wrap_put_att_text (ncid, pct_pft_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_pft_id, 'units'    , unit)

! End of definitions

  status = nf_enddef(ncid)

! Read in formatted surface data

  open (unit=ndata,file=trim(filei),status='unknown',&
       form='formatted',iostat=status)
  if (status .ne. 0) then
     write (6,*)'failed to open ',trim(filei),' on unit ',&
          ndata,' ierr=',status
     stop
  end if
  do j = 1, nlat
     do i = 1, nlon
        read (ndata,*) landmask(i,j),(pct_pft(i,j,m),m=0,numpft)
        if (landmask(i,j) == 100.) landmask(i,j) = 1.
     end do
  end do
  close(ndata)

! Write variables

  call wrap_put_var_realx (ncid, lon_id     , lon)
  call wrap_put_var_realx (ncid, lat_id     , lat)
  call wrap_put_var_realx (ncid, longxy_id  , longxy)
  call wrap_put_var_realx (ncid, latixy_id  , latixy)
  call wrap_put_var_realx (ncid, edgen_id   , edge(1))
  call wrap_put_var_realx (ncid, edgee_id   , edge(2))
  call wrap_put_var_realx (ncid, edges_id   , edge(3))
  call wrap_put_var_realx (ncid, edgew_id   , edge(4))
  call wrap_put_var_realx (ncid, landmask_id, landmask)
  call wrap_put_var_realx (ncid, pct_pft_id , pct_pft)

  call wrap_close(ncid)

end program convert_pft

!===============================================================================

subroutine wrap_create (path, cmode, ncid)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  character(len=*) path
  integer cmode, ncid, ret
  ret = nf_create (path, cmode, ncid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_create

!===============================================================================

subroutine wrap_def_dim (nfid, dimname, len, dimid)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, len, dimid
  character(len=*) :: dimname
  integer ret
  ret = nf_def_dim (nfid, dimname, len, dimid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_def_dim

!===============================================================================

subroutine wrap_def_var (nfid, name, xtype, nvdims, vdims, varid)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, xtype, nvdims, varid
  integer :: vdims(nvdims)
  character(len=*) :: name
  integer ret
  ret = nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_def_var

!===============================================================================

subroutine wrap_put_att_text (nfid, varid, attname, atttext)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, varid
  character(len=*) :: attname, atttext
  integer :: ret, siz
  siz = len_trim(atttext)
  ret = nf_put_att_text (nfid, varid, attname, siz, atttext)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_att_text

!===============================================================================

subroutine wrap_put_var_realx (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, varid
  real(r8) :: arr(*)
  integer :: ret
#ifdef CRAY
  ret = nf_put_var_real (nfid, varid, arr)
#else
  ret = nf_put_var_double (nfid, varid, arr)
#endif
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_var_realx

!===============================================================================

subroutine wrap_put_var_int (nfid, varid, arr)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: nfid, varid
  integer :: arr(*)
  integer :: ret
  ret = nf_put_var_int (nfid, varid, arr)
  if (ret.ne.NF_NOERR) call handle_error (ret)
end subroutine wrap_put_var_int
  
!===============================================================================

subroutine wrap_close (ncid)
  implicit none
  include 'netcdf.inc'
  integer, parameter :: r8 = selected_real_kind(12)
  integer :: ncid
  integer :: ret
  ret = nf_close (ncid)
  if (ret.ne.NF_NOERR) then
     write(6,*)'WRAP_CLOSE: nf_close failed for id ',ncid
     call handle_error (ret)
  end if
end subroutine wrap_close

!===============================================================================

subroutine handle_error(ret)
  implicit none
  include 'netcdf.inc'
  integer :: ret
  if (ret .ne. nf_noerr) then
     write(6,*) 'NCDERR: ERROR: ',nf_strerror(ret)
     call abort
  endif
end subroutine handle_error

!===============================================================================





