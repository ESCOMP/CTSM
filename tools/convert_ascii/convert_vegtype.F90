program make_surftype

  implicit none
#include <netcdf.inc>

!-----------------------------------------------------------------
! make surface type netcdf file
!-----------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

  integer, parameter :: nlon = 720  !input grid : longitude points
  integer, parameter :: nlat = 360  !input grid : latitude  points
  integer, parameter :: nlsm = 28   !number of LSM surface types

  real(r8) :: longxy(nlon,nlat)     !longitude dimension array (2d)        
  real(r8) :: latixy(nlon,nlat)     !longitude dimension array (2d)
  real(r8) :: lon(nlon)            
  real(r8) :: lat(nlat)
  real(r8) :: edge(4)
  real(r8) :: dx,dy

  integer :: surtyp(nlon,nlat)      !input surface type

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id
  integer :: nvegmax                !maximum value for input surface data

  integer :: lon_id                 !longitude array id
  integer :: lat_id                 !latitude array id
  integer :: longxy_id              !2d longitude array id
  integer :: latixy_id              !2d latitude array id
  integer :: edgen_id               !northern edge of grid (edge(1)) id
  integer :: edgee_id               !eastern  edge of grid (edge(2)) id
  integer :: edges_id               !southern edge of grid (edge(3)) id
  integer :: edgew_id               !western  edge of grid (edge(4)) id
  integer :: surftyp_id             !surface type id  

  integer :: i,j,k,m,n              !indices
  integer :: ndata = 1              !input unit
  integer :: ncid                   !netCDF file id
  integer :: dim1_id(1)             !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)             !netCDF dimension id for 2-d variables
  integer :: status                 !netCDF status

  integer :: miss=99999             !missing data indicator
  integer :: surftyp_i(nlon,nlat)   !input surface type (Olson or LSM)
  integer :: surftyp_o(nlon,nlat)   !output surface type (LSM)
  integer :: in2lsm(100)            !LSM surface type for each input type
  integer :: i_w                    !grid cell to west
  integer :: i_e                    !grid cell to east
  integer :: j_n                    !grid cell to north
  integer :: j_s                    !grid cell to south
  
  character(len=80) :: filei, fileo
  character(len=80) :: fvegtyp
  character(len=80) :: name,unit


!-----------------------------------------------------------------

! Determine input and output file names

  filei = 'olson.dat'
  fileo = 'olson.nc'

! -----------------------------------------------------------------
! Read in input data. 
! Data can be in one of two forms:
! OLSON vegetation types or LSM vegetation types. 
!  o The LSM vegetation types are from 0 to [nlsm]. 
!  o The OLSON vegetation types are from 0 to 73 and 
!    are mapped to LSM types.
!
! Data are 1/2 x 1/2 degree, stored in latitude bands, 
! from north to south. In a given latitude band, data begin 
! at Dateline (180W) and proceed eastward. So first data  
! point (x(1,1)) is a box centered at 89.75N, 179.75W. 
!
!   90.0N  ---------------------  
!          |         |         |
!          |    x    |    x    |
!          |  (1,1)  |  (2,1)  |
!          |         |         |
!   89.5N  --------------------- 
!        180.0W    179.5W    179.0W
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
! create netcdf file
! -----------------------------------------------------------------

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid,nf_global,'data_type','surface_type_data')

! Define dimensions

  call wrap_def_dim (ncid, 'LON' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'LAT' , nlat, dimlat_id)

! Define grid variables

  dim1_id(1) = dimlon_id

  dim2_id(1) = dimlon_id
  dim2_id(2) = dimlat_id

  name = 'lon'
  unit = 'degrees east'
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
  call wrap_def_var (ncid, 'LONGXY', nf_float, 2, dim2_id, longxy_id)
  call wrap_put_att_text (ncid, longxy_id, 'long_name', name)
  call wrap_put_att_text (ncid, longxy_id, 'units'    , unit)

  name = 'latitude-2d'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'LATIXY', nf_float, 2, dim2_id, latixy_id)
  call wrap_put_att_text (ncid, latixy_id, 'long_name', name)
  call wrap_put_att_text (ncid, latixy_id, 'units'    , unit)

  name = 'northern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGEN', nf_float, 0, 0, edgen_id)
  call wrap_put_att_text (ncid, edgen_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgen_id, 'units'    , unit)

  name = 'southern edge of surface grid'
  unit = 'degrees north'
  call wrap_def_var (ncid, 'EDGES', nf_float, 0, 0, edges_id)
  call wrap_put_att_text (ncid, edges_id, 'long_name', name)
  call wrap_put_att_text (ncid, edges_id, 'units'    , unit)

  name = 'eastern edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEE', nf_float, 0, 0, edgee_id)
  call wrap_put_att_text (ncid, edgee_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgee_id, 'units'    , unit)

  name = 'western edge of surface grid'
  unit = 'degrees east'
  call wrap_def_var (ncid, 'EDGEW', nf_float, 0, 0, edgew_id)
  call wrap_put_att_text (ncid, edgew_id, 'long_name', name)
  call wrap_put_att_text (ncid, edgew_id, 'units'    , unit)

  name = 'surface_type'
  unit = 'unitless'
  call wrap_def_var (ncid ,'SURFACE_TYPE' ,nf_int, 2, dim2_id, surftyp_id)
  call wrap_put_att_text (ncid, surftyp_id, 'long_name', name)
  call wrap_put_att_text (ncid, surftyp_id, 'units'    , unit)

  status = nf_enddef(ncid)

! -----------------------------------------------------------------
! read in formatted surface data
! -----------------------------------------------------------------

  open (unit=ndata,file=filei,status='unknown',form='formatted',iostat=status)
  if (status .ne. 0) then
     write (6,*)'failed to open ',trim(filei),' on unit ',ndata,' ierr=',status
     stop
  end if
  nvegmax = 0
  do j = 1, nlat
     do i = 1, nlon
        read (ndata,*) surftyp_i(i,j)
        nvegmax = max(nvegmax,surftyp_i(i,j))
     end do
  end do
  close(ndata)

! -----------------------------------------------------------------
! Olson data : Convert input surface types to LSM surface types
! -----------------------------------------------------------------

  if (nvegmax > nlsm) then

     do j = 1, nlat
        do i = 1, nlon

           k = surftyp_i(i,j)
           
! There are several values (2, 6, 8) in the data that are not defined. 
! Use neighboring cells in this order: west, east, north, south

           i_w = max(   1,i-1)
           i_e = min(nlon,i+1)
           j_n = min(   1,j-1)
           j_s = max(nlat,j+1)
           
           if (k==2 .or. k==6 .or. k==8) k = surftyp_i(i_w,j  )
           if (k==2 .or. k==6 .or. k==8) k = surftyp_i(i_e,j  )
           if (k==2 .or. k==6 .or. k==8) k = surftyp_i(i  ,j_n)
           if (k==2 .or. k==6 .or. k==8) k = surftyp_i(i  ,j_s)
           
! Split antarctica (17) into polar desert (69) and ice (70)
 
           if (k == 17) then
              if (j <= 313) then
                 k = 69
              else
                 k = 70
              end if
           end if

! 61 (eastern south taiga) will be classified as needleleaf deciduous tree. 
! Change 61 to 20 (main taiga = needleleaf evergreen tree) based on longitude

           if (k==61 .and. i<=576) k = 20 

! 61 (eastern south taiga) will be classified needleleaf deciduous tree. 
! Create additional needleleaf deciduous tree from 21 (main taiga) and
! 60 (southern taiga) based on longitude

           if (k==21 .and. i>=555) k = 61  
           if (k==60 .and. i>=582) k = 61   

! Change 26 (warm mixed) to broad-leaved humid forest based on latitude

           if (k==26 .and. j>=113) k = 29

! Split forest tundra (62, 63) into needleleaf evergreen forest tundra (62) 
! and needleleaf deciduous forest tundra (63) based on longitude

           if (k==63) k = 62
           if (k==62 .and. i>=490) k = 63   

! Error check

           if (k>100 .or. k<0) then
              write (6,*) 'ERROR: Olson surface type = ',k,' is undefined for lon,lat = ',i,j
              Stop
           end if

! Save modified OLSON type 

           surftyp_i(i,j) = k

        end do
     end do

! Assign each of the OLSON surface types to an LSM surface type.
! This mapping from OLSON to LSM is based on the BATS dataset code.
! Note: in2lsm(i) = OLSON type i

     in2lsm(1:19) = miss
     in2lsm(20) = 3 
     in2lsm(21) = 3                                                     
     in2lsm(22) = 3  
     in2lsm(23) = 6 
     in2lsm(24) = 8 
     in2lsm(25) = 9                                               
     in2lsm(26) = 9                                               
     in2lsm(27) = 7       
     in2lsm(28) = 10                                               
     in2lsm(29) = 10                                               
     in2lsm(30) = 24                                               
     in2lsm(31) = 26                                               
     in2lsm(32) = 12                                               
     in2lsm(33) = 10                                               
     in2lsm(34) = miss                                           
     in2lsm(35) = miss                                           
     in2lsm(36) = 28                                              
     in2lsm(37) = 25                                              
     in2lsm(38) = 23                                              
     in2lsm(39) = 23                                              
     in2lsm(40) = 17                                               
     in2lsm(41) = 18                                               
     in2lsm(42) = 17                                               
     in2lsm(43) = 12                                               
     in2lsm(44) = 28                                              
     in2lsm(45) = 28                                              
     in2lsm(46) = 20                                              
     in2lsm(47) = 20                                              
     in2lsm(48) = 20   
     in2lsm(49) = 22                                              
     in2lsm(50) = 2                                               
     in2lsm(51) = 22                                              
     in2lsm(52) = 22   
     in2lsm(53) = 19                                               
     in2lsm(54) = 19                                               
     in2lsm(55) = 15    
     in2lsm(56) = 16   
     in2lsm(57) = 15   
     in2lsm(58) = 16    
     in2lsm(59) = 21                                              
     in2lsm(60) = 6                                              
     in2lsm(61) = 4                                               
     in2lsm(62) = 13    
     in2lsm(63) = 14                                              
     in2lsm(64) = 20                                              
     in2lsm(65) = 0   
     in2lsm(66) = 0   
     in2lsm(67) = 0   
     in2lsm(68) = 0   
     in2lsm(69) = 2                                               
     in2lsm(70) = 1                                             
     in2lsm(71) = 22                                              
     in2lsm(72) = 27                                              
     in2lsm(73) = 0   
     in2lsm(74:100) = miss

  endif

! -----------------------------------------------------------------
! LSM input data : 1:1 correspondence between surface types
! -----------------------------------------------------------------

  if (nvegmax == nlsm) then
     do i = 1, nlsm
        in2lsm(i) = i
     end do
  end if

! -----------------------------------------------------------------
! Transform input surface types to LSM surface types
! -----------------------------------------------------------------

  surftyp_o(:,:) = miss 
  do j = 1 , nlat                                                 
     do i = 1, nlon                                                 
        if (surftyp_i(i,j) == 0) then
           surftyp_o(i,j) = 0
        else
           k = surftyp_i(i,j)
           surftyp_o(i,j) = in2lsm(k)
        end if
        if (surftyp_o(i,j)>nlsm .or. surftyp_o(i,j)<0) then
           write (6,*) 'ERROR: LSM surface type = ',surftyp_o(i,j),' is undefined for lon,lat = ',i,j
           stop
        end if
     end do
  end do

! Write output variables

  call wrap_put_var_realx (ncid, lon_id    , lon)
  call wrap_put_var_realx (ncid, lat_id    , lat)
  call wrap_put_var_realx (ncid, longxy_id , longxy)
  call wrap_put_var_realx (ncid, latixy_id , latixy)
  call wrap_put_var_realx (ncid, edgen_id  , edge(1))
  call wrap_put_var_realx (ncid, edgen_id  , edge(1))
  call wrap_put_var_realx (ncid, edgee_id  , edge(2))
  call wrap_put_var_realx (ncid, edges_id  , edge(3))
  call wrap_put_var_realx (ncid, edgew_id  , edge(4))
  call wrap_put_var_int   (ncid, surftyp_id, surftyp_o)

  call wrap_close(ncid)

end program make_surftype

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





