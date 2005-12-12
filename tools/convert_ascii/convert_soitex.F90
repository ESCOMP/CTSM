program convert_soitex

  implicit none
  include 'netcdf.inc'

!-----------------------------------------------------------------
! make surface type netcdf file
!-----------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

! File specific settings

  integer, parameter :: nlon = 4320      !input grid : longitude points
  integer, parameter :: nlat = 2160      !input grid : latitude  points
  integer, parameter :: nlay = 10        !input grid : number of soil layers
  integer, parameter :: nmapunits = 4931 !input grid : # of igbp soil 'mapunits'
  integer, parameter :: mapunitmax = 6998!input grid : max value of 'mapunits'

  real(r8) :: dzsoi(10), zsoi(10)        !soil layer thickness and depth

  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid
  real(r8) :: dx,dy                       !grid increments

  integer  :: mu                          !current mapunit
  real(r8) :: pct_sand(mapunitmax,nlay)   !pct sand 
  real(r8) :: pct_clay(mapunitmax,nlay)   !pct clay
  real(r8) :: landmask(nlon,nlat)         !land mask
  real(r8) :: mapunit(nlon,nlat)          !global map of igbp soil mapunits
  real(r8) :: temp(nlon,nlat)             !same; used as temporary buffer

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id
  integer :: dimlay_id                    !netCDF dimension id
  integer :: dimmapunits_id               !netCDF dimension id
  integer :: dimmapunitmax_id             !netCDF dimension id

  integer :: dzsoi_id                     !soil thickness by layer
  integer :: zsoi_id                      !soil depth by layer
  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: lay_id                       !1d layer array id
  integer :: mapunit_id                   !2d mapunits array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: pct_sand_id                  !sand id
  integer :: pct_clay_id                  !clay id
  integer :: landmask_id                  !landmask id

  integer :: i,j,k                        !indices
  integer :: ndata = 1                    !input unit
  integer :: ncid                         !netCDF file id
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) ::  filei,   fileo   !input,output filenames
  character(len=256) ::  filei1,  filei2,  filei3,  filei4,  filei5,  filei6
  character(len=256) ::  filei7,  filei8,  filei9, filei10, filei11, filei12
  character(len=256) :: filei13, filei14, filei15, filei16, filei17, filei18
  character(len=256) :: filei19, filei20
  character(len=256) :: name,unit         !netCDF attributes

!-----------------------------------------------------------------

! Determine input and output file names

  filei = '/ptmp/slevis/lsminput/old/world.ascii'
  filei1= '/ptmp/slevis/lsminput/old/clm/Clay2.map'
  filei2= '/ptmp/slevis/lsminput/old/clm/Clay5.map'
  filei3= '/ptmp/slevis/lsminput/old/clm/Clay9.map'
  filei4= '/ptmp/slevis/lsminput/old/clm/Clay17.map'
  filei5= '/ptmp/slevis/lsminput/old/clm/Clay29.map'
  filei6= '/ptmp/slevis/lsminput/old/clm/Clay49.map'
  filei7= '/ptmp/slevis/lsminput/old/clm/Clay83.map'
  filei8= '/ptmp/slevis/lsminput/old/clm/Clay138.map'
  filei9= '/ptmp/slevis/lsminput/old/clm/Clay230.map'
  filei10= '/ptmp/slevis/lsminput/old/clm/Clay343.map'
  filei11= '/ptmp/slevis/lsminput/old/clm/Sand2.map'
  filei12= '/ptmp/slevis/lsminput/old/clm/Sand5.map'
  filei13= '/ptmp/slevis/lsminput/old/clm/Sand9.map'
  filei14= '/ptmp/slevis/lsminput/old/clm/Sand17.map'
  filei15= '/ptmp/slevis/lsminput/old/clm/Sand29.map'
  filei16= '/ptmp/slevis/lsminput/old/clm/Sand49.map'
  filei17= '/ptmp/slevis/lsminput/old/clm/Sand83.map'
  filei18= '/ptmp/slevis/lsminput/old/clm/Sand138.map'
  filei19= '/ptmp/slevis/lsminput/old/clm/Sand230.map'
  filei20= '/ptmp/slevis/lsminput/old/clm/Sand343.map'
  fileo = '/ptmp/slevis/lsminput/new/clm/5minx5min/mksrf_soitex.nc'

! -----------------------------------------------------------------
! Determine grid for input data 
!
! IGBP input data
!
! Data are igbp soil mapunits. Each mapunit can be thought of as a
! unique soil profile. Map resolution is 5min x 5min, stored in
! latitude bands, from north to south. In a given latitude band,
! data begin at the dateline (180W) and proceed eastward.

! Also from igbp we read the %sand and %clay that correspond to each
! mapunit by LSM soil layer. The raw data set created here, assembles all
! this data in one netcdf file. The final product goes from south to north.
! -----------------------------------------------------------------

! Define soil thicknesses and soil depths (model dependent)

  zsoi(1) = 0.0175
  zsoi(2) = 0.0451
  zsoi(3) = 0.0906
  zsoi(4) = 0.1656
  zsoi(5) = 0.2892
  zsoi(6) = 0.4930
  zsoi(7) = 0.8290
  zsoi(8) = 1.3829
  zsoi(9) = 2.2962
  zsoi(10) = 3.4332

  dzsoi(1) = zsoi(1) - 0.
  dzsoi(2) = zsoi(2) - zsoi(1)
  dzsoi(3) = zsoi(3) - zsoi(2)
  dzsoi(4) = zsoi(4) - zsoi(3)
  dzsoi(5) = zsoi(5) - zsoi(4)
  dzsoi(6) = zsoi(6) - zsoi(5)
  dzsoi(7) = zsoi(7) - zsoi(6)
  dzsoi(8) = zsoi(8) - zsoi(7)
  dzsoi(9) = zsoi(9) - zsoi(8)
  dzsoi(10) = zsoi(10) - zsoi(9)

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
        latixy(i,j) = (edge(1)-dy/2.) - (j-1)*dy
        longxy(i,j) = (edge(4)+dx/2.) + (i-1)*dx
       end do
  end do

  lat(:) = latixy(1,:)
  lon(:) = longxy(:,1)

! -----------------------------------------------------------------
! create netcdf file
! -----------------------------------------------------------------

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'igbp_soil_texture_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)
  call wrap_def_dim (ncid, 'number_of_layers'   , nlay      , dimlay_id)
  call wrap_def_dim (ncid, 'number_of_mapunits' , nmapunits , dimmapunits_id)
  call wrap_def_dim (ncid, 'max_value_mapunit'  , mapunitmax, dimmapunitmax_id)

! Define input file independent variables 

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
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! to possibly replace the next two variables
! find out about dimensioned variables
! (eg, see how pressure levels are treated)

  name = 'soil layer thickness'
  unit = 'm'
  dim1_id(1) = dimlay_id
  call wrap_def_var (ncid,'DZSOI', nf_float, 1, dim1_id, dzsoi_id)
  call wrap_put_att_text (ncid, dzsoi_id, 'long_name', name)
  call wrap_put_att_text (ncid, dzsoi_id, 'units'    , unit)

  name = 'soil layer depth'
  unit = 'm'
  dim1_id(1) = dimlay_id
  call wrap_def_var (ncid,'ZSOI', nf_float, 1, dim1_id, zsoi_id)
  call wrap_put_att_text (ncid, zsoi_id, 'long_name', name)
  call wrap_put_att_text (ncid, zsoi_id, 'units'    , unit)

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

! Define soil type and soil texture variables

  name = 'igbp soil mapunit'
  unit = 'unitless'
  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id
  call wrap_def_var (ncid ,'MAPUNITS' ,nf_float, 2, dim2_id, mapunit_id)
  call wrap_put_att_text (ncid, mapunit_id, 'long_name', name)
  call wrap_put_att_text (ncid, mapunit_id, 'units'    , unit)

  name = 'percent sand'
  unit = 'unitless'
  dim2_id(1) = dimmapunitmax_id
  dim2_id(2) = dimlay_id
  call wrap_def_var (ncid ,'PCT_SAND' ,nf_float, 2, dim2_id, pct_sand_id)
  call wrap_put_att_text (ncid, pct_sand_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_sand_id, 'units'    , unit)

  name = 'percent clay'
  unit = 'unitless'
  call wrap_def_var (ncid ,'PCT_CLAY' ,nf_float, 2, dim2_id, pct_clay_id)
  call wrap_put_att_text (ncid, pct_clay_id, 'long_name', name)
  call wrap_put_att_text (ncid, pct_clay_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)

! Read in formatted surface data

  open (unit=ndata,file=trim(filei),status='old',form='formatted',iostat=status)
  if (status .ne. 0) then
     write (6,*)'failed to open ',trim(filei),' on unit ',ndata,' ierr=',status
     stop
  end if
  do j = 1, nlat
     if (lat(j) <= 84. .and. lat(j) >= -56.5) then
        read (ndata,*) (mapunit(i,j),i=1,nlon)
        do i = 1, nlon
           if (mapunit(i,j) ==    0. .or. mapunit(i,j) ==  794. .or. &
               mapunit(i,j) == 1972. .or. mapunit(i,j) == 3214. .or. &
               mapunit(i,j) == 6997. .or. mapunit(i,j) == 6998.) then
              landmask(i,j) = 0. !ocean, no soil data, lakes, glaciers
           else
              landmask(i,j) = 1.
           end if
        end do
     else
        landmask(i,j) = 0.
        mapunit(i,j) = 0.
     end if
  end do
  close(ndata)

  open (unit=11,file=trim(filei1),status='old')
  open (unit=12,file=trim(filei2),status='old')
  open (unit=13,file=trim(filei3),status='old')
  open (unit=14,file=trim(filei4),status='old')
  open (unit=15,file=trim(filei5),status='old')
  open (unit=16,file=trim(filei6),status='old')
  open (unit=17,file=trim(filei7),status='old')
  open (unit=18,file=trim(filei8),status='old')
  open (unit=19,file=trim(filei9),status='old')
  open (unit=20,file=trim(filei10),status='old')
  open (unit=21,file=trim(filei11),status='old')
  open (unit=22,file=trim(filei12),status='old')
  open (unit=23,file=trim(filei13),status='old')
  open (unit=24,file=trim(filei14),status='old')
  open (unit=25,file=trim(filei15),status='old')
  open (unit=26,file=trim(filei16),status='old')
  open (unit=27,file=trim(filei17),status='old')
  open (unit=28,file=trim(filei18),status='old')
  open (unit=29,file=trim(filei19),status='old')
  open (unit=30,file=trim(filei20),status='old')

! initialize first

  do j = 1, nlay
     do i = 1, mapunitmax
        pct_clay(i,j) = 0.
        pct_sand(i,j) = 0.
     end do
  end do

! first clay

  do j = 1, nlay
     read(10+j,*) ! clear the first line
     do i = 1, nmapunits
        read(10+j,*) mu, pct_clay(mu,j)
     end do
     close(10+j)
  end do

! then sand

  do j = 1, nlay
     read(20+j,*) ! clear the first line
     do i = 1, nmapunits
        read(20+j,*) mu, pct_sand(mu,j)
     end do
     close(20+j)
  end do

! make north to south back to south to north

  do j = 1, nlat
     do i = 1, nlon
        temp(i,j) = mapunit(i,nlat-j+1)
     end do
  end do
  do j = 1, nlat
     do i = 1, nlon
        mapunit(i,j) = temp(i,j)
     end do
  end do
  do j = 1, nlat
     do i = 1, nlon
        temp(i,j) = landmask(i,nlat-j+1)
     end do
  end do
  do j = 1, nlat
     do i = 1, nlon
        landmask(i,j) = temp(i,j)
     end do
  end do
  do j = 1, nlat
     do i = 1, nlon
        temp(i,j) = latixy(i,nlat-j+1)
     end do
  end do
  do j = 1, nlat
     do i = 1, nlon
        latixy(i,j) = temp(i,j)
     end do
  end do
  do j = 1, nlat
     do i = 1, nlon
        temp(i,j) = longxy(i,nlat-j+1)
     end do
  end do
  do j = 1, nlat
     do i = 1, nlon
        longxy(i,j) = temp(i,j)
     end do
  end do

  lat(:) = latixy(1,:)
  lon(:) = longxy(:,1)

! Create output file

  call wrap_put_var_realx (ncid, lon_id     , lon)
  call wrap_put_var_realx (ncid, lat_id     , lat)
  call wrap_put_var_realx (ncid, longxy_id  , longxy)
  call wrap_put_var_realx (ncid, latixy_id  , latixy)
  call wrap_put_var_realx (ncid, landmask_id, landmask)
  call wrap_put_var_realx (ncid, edgen_id   , edge(1))
  call wrap_put_var_realx (ncid, edgee_id   , edge(2))
  call wrap_put_var_realx (ncid, edges_id   , edge(3))
  call wrap_put_var_realx (ncid, edgew_id   , edge(4))
  call wrap_put_var_realx (ncid, dzsoi_id   , dzsoi)
  call wrap_put_var_realx (ncid, zsoi_id    , zsoi)
  call wrap_put_var_realx (ncid, mapunit_id , mapunit)
  call wrap_put_var_realx (ncid, pct_sand_id, pct_sand)
  call wrap_put_var_realx (ncid, pct_clay_id, pct_clay)

  call wrap_close(ncid)

end program convert_soitex

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





