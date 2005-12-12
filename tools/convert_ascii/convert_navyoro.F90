program convert_navyoro

  implicit none
  include 'netcdf.inc'

!-----------------------------------------------------------------
! make surface type netcdf file
!-----------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

! File specific settings

  integer, parameter :: nlon = 1080      !number of navy oro longitudes
  integer, parameter :: nlat =  540      !number of navy oro latitudes
  real(r8) :: longxy(nlon,nlat)          !longitude dimension array (2d)        
  real(r8) :: latixy(nlon,nlat)          !longitude dimension array (2d)        
  real(r8) :: edge(4)                    !N,E,S,W edges of grid
  integer  :: numlon(nlat)               !fractional land                       
  real(r8) :: landfrac(nlon,nlat)        !number of longitudes for each latitude
  real(r8) :: dx,dy                      !grid increments
                                         
  integer :: dimlon_id                   !netCDF dimension id
  integer :: dimlat_id                   !netCDF dimension id
  integer :: longxy_id                   !2d longitude array id
  integer :: latixy_id                   !2d latitude array id
  integer :: edgen_id                    !northern edge of grid (edge(1)) id
  integer :: edgee_id                    !eastern  edge of grid (edge(2)) id
  integer :: edges_id                    !southern edge of grid (edge(3)) id
  integer :: edgew_id                    !western  edge of grid (edge(4)) id
  integer :: landfrac_id                 !fractional land id
  integer :: numlon_id                   !numlon id
                                         
  integer :: i,j,k                       !indicis
  integer :: ndata = 1                   !input unit
  integer :: ncid                        !netCDF file id
  integer :: dim1_id(1)                  !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                  !netCDF dimension id for 2-d variables
  integer :: status                      !status
                                         
  character(len=256) :: filei, fileo     !input,output filenames
  character(len=256) :: name,unit        !netCDF attributes

  real(r8) :: lat_s(nlat+1)              !grid cell latitude, southern edge (degrees)
  real(r8) :: lon_w(nlon+1)              !grid cell longitude, western edge (degrees)


!-----------------------------------------------------------------

! Determine input and output file names

  filei = 'navy_oro.fland.20min.dat'
  fileo = 'mksrf_navyoro_20min.nc'

! -----------------------------------------------------------------
! Create netcdf file
! -----------------------------------------------------------------

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'fractional_land_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

! Define grid variables 

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

  name = 'number of longitudes per latitude band'
  unit = 'unitless'
  dim1_id(1) = dimlat_id
  call wrap_def_var (ncid ,'NUMLON', nf_int, 1, dim1_id, numlon_id)
  call wrap_put_att_text (ncid, numlon_id, 'long_name', name)
  call wrap_put_att_text (ncid, numlon_id, 'units'    , unit)

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

! Define fractional land 

  name = 'fractional land'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'LANDFRAC', nf_float, 2, dim2_id, landfrac_id)
  call wrap_put_att_text (ncid, landfrac_id, 'long_name', name)
  call wrap_put_att_text (ncid, landfrac_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)

!----------------------------------------------------------------
! Read in formatted fractional land data and create output file
!----------------------------------------------------------------

! Read fractional land

  open (unit=ndata,file=trim(filei),status='unknown',form='formatted',iostat=status)
  if (status .ne. 0) then
     write (6,*)'failed to open ',trim(filei),' on unit ',ndata,' ierr=',status
     stop
  else
     write (6,*)'opened ',trim(filei),' on unit ',ndata
  end if

  do j = 1, nlat
     do i = 1, nlon
        read (ndata,*) landfrac(i,j)
     end do
  end do

  close(ndata)

! Input grid is regular grid starting at greenwich with the western edge 
! ON greenwich

  dy = 180./nlat
  do j = 1, nlat+1
     lat_s(j) = -90.0 + (j-1)*dy
  end do

  dx = 360./nlon
  do i = 1, nlon+1
     lon_w(i) = 0. + (i-1)*dx
  end do

  do j = 1, nlat
     do i = 1, nlon
        latixy(i,j) = (lat_s(j)+lat_s(j+1))/2.
        longxy(i,j) = (lon_w(i)+lon_w(i+1))/2.
     end do
  end do

  edge(1) =   90.
  edge(2) =  360.
  edge(3) =  -90.
  edge(4) =    0.

  numlon(:) = nlon

! Write output variables

  call wrap_put_var_realx (ncid, longxy_id   , longxy)
  call wrap_put_var_realx (ncid, latixy_id   , latixy)
  call wrap_put_var_realx (ncid, edgen_id    , edge(1))
  call wrap_put_var_realx (ncid, edgee_id    , edge(2))
  call wrap_put_var_realx (ncid, edges_id    , edge(3))
  call wrap_put_var_realx (ncid, edgew_id    , edge(4))
  call wrap_put_var_int   (ncid, numlon_id   , numlon)
  call wrap_put_var_realx (ncid, landfrac_id , landfrac)

! Close output file

  call wrap_close(ncid)

end program convert_navyoro

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





