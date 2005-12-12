program convert_soicol

  implicit none
  include 'netcdf.inc'

!-----------------------------------------------------------------
! make soil colr netcdf file
!-----------------------------------------------------------------

  integer, parameter :: r8 = selected_real_kind(12)

! File specific settings

  integer, parameter :: nlon = 128        !input grid : longitude points
  integer, parameter :: nlat =  64        !input grid : latitude  points

  real(r8) :: gaulat(nlat)                !input grid: Gaussian latitudes
  real(r8) :: gauwt(nlat)                 !input grid: Gaussian weights

  real(r8) :: lon(nlon)                   !longitude dimension array (1d)
  real(r8) :: lat(nlat)                   !latitude dimension array (1d) 
  real(r8) :: longxy(nlon,nlat)           !longitude dimension array (2d)        
  real(r8) :: latixy(nlon,nlat)           !longitude dimension array (2d)
  real(r8) :: edge(4)                     !N,E,S,W edges of grid
  real(r8) :: dx,dy                       !grid increments

  real(r8) :: soil_color(nlon,nlat)       !lsm soil color
  real(r8) :: landmask(nlon,nlat)         !land mask derived from soil color

  integer :: dimlon_id                    !netCDF dimension id
  integer :: dimlat_id                    !netCDF dimension id

  integer :: lon_id                       !1d longitude array id
  integer :: lat_id                       !1d latitude array id
  integer :: longxy_id                    !2d longitude array id
  integer :: latixy_id                    !2d latitude array id
  integer :: edgen_id                     !northern edge of grid (edge(1)) id
  integer :: edgee_id                     !eastern  edge of grid (edge(2)) id
  integer :: edges_id                     !southern edge of grid (edge(3)) id
  integer :: edgew_id                     !western  edge of grid (edge(4)) id
  integer :: soil_color_id                !soil color id
  integer :: landmask_id                  !landmask id

  integer :: i,j                          !indicis
  integer :: ndata = 1                    !input unit
  integer :: ncid                         !netCDF file id
  integer :: dim1_id(1)                   !netCDF dimension id for 1-d variables
  integer :: dim2_id(2)                   !netCDF dimension id for 2-d variables
  integer :: status                       !status

  character(len=256) :: filei, fileo      !input,output filenames
  character(len=256) :: name,unit         !netCDF attributes

!-----------------------------------------------------------------

! Determine input and output file names

  filei = '/ptmp/slevis/lsmv2_2/input/0.5x0.5/bats.dat'
  fileo = '/ptmp/slevis/lsmv2_2/input/0.5x0.5/mksrf_soicol.nc'

! -----------------------------------------------------------------
! Determine grid for input data 
!
! BATS data are on T42 Gaussian grid, approximately 2.8 x 2.8 degrees,
! stored in latitude bands, from south to north. In a given latitude band, 
! data begin at Greenwich, centered on Greenwich, and proceed eastward. 
! -----------------------------------------------------------------

! Define North, East, South, West edges of grid

  edge(1) =  90.
  edge(2) = 358.59375
  edge(3) = -90.
  edge(4) =  -1.40625

! Make latitudes and longitudes at center of grid cell

  dx = (edge(2)-edge(4)) / nlon
  call mkgaulat (gaulat, gauwt, nlat)

  do j = 1, nlat
     do i = 1, nlon
        latixy(i,j) = -asin(gaulat(j)) * 180./(4.*atan(1.))
        longxy(i,j) = (edge(4)+dx/2.) + (i-1.)*dx
     end do
  end do

  lat(:) = latixy(1,:)
  lon(:) = longxy(:,1)

! -----------------------------------------------------------------
! create netcdf file
! -----------------------------------------------------------------

  call wrap_create (fileo, nf_clobber, ncid)
  call wrap_put_att_text (ncid, nf_global, 'data_type', 'soil_color_data')

! Define dimensions

  call wrap_def_dim (ncid, 'lon' , nlon, dimlon_id)
  call wrap_def_dim (ncid, 'lat' , nlat, dimlat_id)

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

! Define input file specific variables

  name = 'soil color'
  unit = 'unitless'
  dim2_id(1) = lon_id
  dim2_id(2) = lat_id
  call wrap_def_var (ncid ,'SOIL_COLOR' ,nf_float, 2, dim2_id, soil_color_id)
  call wrap_put_att_text (ncid, soil_color_id, 'long_name', name)
  call wrap_put_att_text (ncid, soil_color_id, 'units'    , unit)

  name = 'land mask'
  unit = 'unitless'
  call wrap_def_var (ncid ,'LANDMASK' ,nf_float, 2, dim2_id, landmask_id)
  call wrap_put_att_text (ncid, landmask_id, 'long_name', name)
  call wrap_put_att_text (ncid, landmask_id, 'units'    , unit)

! End of definition

  status = nf_enddef(ncid)

! Read in formatted surface data

  open (unit=ndata,file=trim(filei),status='unknown',form='formatted',iostat=status)
  if (status .ne. 0) then
     write (6,*)'failed to open ',trim(filei),' on unit ',ndata,' ierr=',status
     stop
  end if
  do j = 1, nlat
     do i = 1, nlon
        read (ndata,*) soil_color(i,j)
        if (soil_color(i,j)<0 .or. soil_color(i,j)>8) then
           write (6,*) 'ERROR: BATS soil color = ',soil_color(i,j), &
                ' is not valid for lon,lat = ',i,j
           stop
        end if
        if (soil_color(i,j)==0) then
           landmask(i,j) = 0.
        else
           landmask(i,j) = 1.
        end if
     end do
  end do
  close(ndata)

! Create output file

  call wrap_put_var_realx (ncid, lon_id       , lon)
  call wrap_put_var_realx (ncid, lat_id       , lat)
  call wrap_put_var_realx (ncid, longxy_id    , longxy)
  call wrap_put_var_realx (ncid, latixy_id    , latixy)
  call wrap_put_var_realx (ncid, edgen_id     , edge(1))
  call wrap_put_var_realx (ncid, edgee_id     , edge(2))
  call wrap_put_var_realx (ncid, edges_id     , edge(3))
  call wrap_put_var_realx (ncid, edgew_id     , edge(4))
  call wrap_put_var_realx (ncid, landmask_id  , landmask)
  call wrap_put_var_realx (ncid, soil_color_id, soil_color)

! Close output file

  call wrap_close(ncid)

end program convert_soicol

!===============================================================================

      subroutine mkgaulat (a, w, k)

! ------------------------ code history ------------------------------
! source file       : mkgaulat.F
! purpose           : Gaussian latitudes
! date first created: July 1995 - lsm version 1
! by whom           : Gordon Bonan
! date last revised : Jan 2000 - lsm version 2
! by whom           : Mariana Vertenstein
! --------------------------------------------------------------------

      implicit none

      integer, parameter :: r8 = selected_real_kind(12)

! ------------------------ arguments ---------------------------------
      integer , intent(in) :: k       !number of latitudes pole to pole
      real(r8), intent(out):: a(k)    !sine of latitudes
      real(r8), intent(out):: w(k)    !weights
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
      real(r8) pi       !value of pi
      real(r8) eps      !convergence criterion
      real(r8) c        !constant combination
      real(r8) fk       !real(r8) k
      real(r8) xz       !abscissa estimate
      real(r8) pkm1     !|
      real(r8) pkm2     !|-polynomials
      real(r8) pkmrk    !|
      real(r8) pk       !|
      real(r8) sp       !current iteration latitude increment
      real(r8) avsp     !|sp|
      real(r8) fn       !real n
      integer kk        !k/2 (number of latitudes in hemisphere)
      integer is        !latitude index
      integer iter      !iteration counter
      integer j         !index
      integer n         !index
      integer l         !index
      integer nn        !index
      real bz(50)       !table of first 50 zeros
      data bz          / 2.4048255577,   5.5200781103,                 &
         8.6537279129,  11.7915344391,  14.9309177086,  18.0710639679, &
        21.2116366299,  24.3524715308,  27.4934791320,  30.6346064684, &
        33.7758202136,  36.9170983537,  40.0584257646,  43.1997917132, &
        46.3411883717,  49.4826098974,  52.6240518411,  55.7655107550, &
        58.9069839261,  62.0484691902,  65.1899648002,  68.3314693299, &
        71.4729816036,  74.6145006437,  77.7560256304,  80.8975558711, & 
        84.0390907769,  87.1806298436,  90.3221726372,  93.4637187819, &
        96.6052679510,  99.7468198587, 102.8883742542, 106.0299309165, &
       109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309, &
       121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752, &
       134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737, &
       146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/
! --------------------------------------------------------------------

      eps = 1.e-6
      pi = 4.*atan(1.)
 
! The value eps, used for convergence tests in the iterations, 
! can be changed.  Newton iteration is used to find the abscissas.
 
      c = (1.-(2./pi)**2)*0.25
      fk = k
      kk = k/2

! Return n zeros (or if n>50, approximate zeros), of the Bessel function
! j0,in the array a. The first 50 zeros will be given exactly, and the
! remaining zeros are computed by extrapolation,and therefore not exact.

      n = kk
      nn = n
      if (n.gt.50) then
         a(50) = bz(50)
         do j=51,n
            a(j) = a(j-1) + pi
         end do
         nn = 49
      end if
      do j=1,nn
         a(j) = bz(j)
      end do

      do 30 is=1,kk
         xz = cos(a(is)/sqrt((fk+0.5)**2+c))
 
! This is the first approximation to xz
 
         iter = 0
   10    pkm2 = 1.
         pkm1 = xz
         iter = iter + 1

! Error exit

         if (iter.gt.10) then
            write(6,*) 'MKGAULAT error: no convergence in 10 iterations'
            Stop
         end if
 
! Computation of the legendre polynomial
 
         do 20 n=2,k
            fn = n
            pk = ((2.*fn-1.)*xz*pkm1-(fn-1.)*pkm2)/fn
            pkm2 = pkm1
            pkm1 = pk
   20    continue
         pkm1 = pkm2
         pkmrk = (fk*(pkm1-xz*pk))/(1.-xz**2)
         sp = pk/pkmrk
         xz = xz - sp
         avsp = abs(sp)
         if (avsp.gt.eps) go to 10
         a(is) = xz
         w(is) = (2.*(1.-xz**2))/(fk*pkm1)**2
   30 continue
      if (k.ne.kk*2) then
 
! For odd k computation of weight at the equator
 
         a(kk+1) = 0.
         pk = 2./fk**2
         do 40 n=2,k,2
            fn = n
            pk = pk*fn**2/(fn-1.)**2
   40    continue
         w(kk+1) = pk
      end if
 
! Complete the sets of abscissas and weights, using the symmetry.
 
      do 60 n=1,kk
         l = k + 1 - n
         a(l) = -a(n)
         w(l) = w(n)
   60 continue

      return
      end subroutine mkgaulat

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





