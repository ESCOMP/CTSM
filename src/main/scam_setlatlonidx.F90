subroutine scam_setlatlonidx(ncid ,targetlat ,targetlon, closelat, closelon, &
                             closelatidx, closelonidx)

!------------------------------------------------------------------------
! File: scam_setlatlon.F 
!     Returns the latitude and longitude index for use when CLM is running
!     with the Siongle Column Atmospheric Model (SCAM)
!
! Author: John Truesdale (jet@ucar.edu)
!
!------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use clm_varctl,   only: iulog
  use abortutils,   only: endrun
  use ncdio_pio
  
  implicit none

!------------------------------parameters-----------------------------------

  type(file_desc_t), intent(inout)  :: ncid     ! netcdf file id
  real(r8)         , intent(in)  :: targetlat,targetlon
  real(r8)         , intent(out) :: closelat,closelon
  integer          , intent(out) :: closelatidx,closelonidx

!------------------------------locals-----------------------------------
  
  real(r8),allocatable:: lats(:),lons(:),poslons(:)
  real(r8) prev, next,postargetlon
  integer  i                ! 
  integer  status
  integer  nlat_dimID, nlat, lat_varID 
  integer  nlon_dimID, nlon, lon_varID
  integer strt(2),cnt(2)
!
! Check mode: double or single precision.
! 
  call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
  status =  pio_inq_dimid ( ncid, 'lat', nlat_dimID )
  if ( status .ne.  PIO_noerr )  then
     status =  pio_inq_dimid ( ncid, 'lsmlat', nlat_dimID )
     if ( status .ne.  PIO_noerr )  then
        write(iulog,* )'ERROR - setlatlon.F:', &
             'Cant get variable dim for lat or lsmlat'
        call pio_closefile( ncid )
        call endrun
     endif
  endif

  status =  pio_inq_dimid ( ncid, 'lon', nlon_dimID )
  if ( status .ne.  pio_noerr )  then
     status =  pio_inq_dimid ( ncid, 'lsmlon', nlon_dimID )
     if ( status .ne.  pio_noerr )  then
        write(iulog,* )'ERROR - setlatlon.F:', &
             'Cant get variable dim for lon or lsmlon'
        call pio_closefile( ncid )
        call endrun
     endif
  endif
  call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)

  status = pio_inq_dimlen( ncid, nlat_dimid, nlat )
  status = pio_inq_dimlen( ncid, nlon_dimid, nlon )

  allocate(lats(nlat),lons(nlon),poslons(nlon))

  strt(1)=1
  strt(2)=1
  cnt(1)=nlat
  cnt(2)=1

  call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
  status = pio_inq_varid( ncid, 'lat', lat_varID )
  if ( status .ne. pio_noerr ) then
     status = pio_inq_varid( ncid, 'LATIXY', lat_varID )
     if ( status .ne. pio_noerr ) then
        write(iulog,* )'error - setlatlon.F:', &
             'Cant get variable ID for lat or LATIXY'
	call pio_closefile(ncid)
        return
     else
        strt(1) = nlon
        strt(2) = 1
        cnt(1)  = 1
        cnt(2)  = nlat
     endif
  endif
  call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
  status = pio_get_var (ncid, lat_varID, strt, cnt, lats)

  call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
  status = pio_inq_varid( ncid, 'lon', lon_varID )
  if ( status .ne. pio_noerr ) then
     status = pio_inq_varid( ncid, 'LONGXY', lon_varID )
     if ( status .ne. pio_noerr ) then
        write(iulog,* )'ERROR - setlatlon.F:', &
             'Cant get variable ID for lon'
	call pio_closefile(ncid)
        return
     endif
  endif
  call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
  strt(1) = 1
  strt(2) = 1
  cnt(1)  = nlon
  cnt(2)  = 1
  status = pio_get_var (ncid, lon_varID, strt, cnt, lons)

!!$ convert lons array and targetlon to 0,360

  poslons = mod(lons+360._r8,360._r8)
  postargetlon = mod(targetlon+360._r8,360._r8)

!!$ find index of value closest to 0.

  closelonidx = (MINLOC(abs(poslons-postargetlon),dim=1))
  closelatidx = (MINLOC(abs(lats-targetlat),dim=1))
  closelon = lons(closelonidx)
  closelat = lats(closelatidx)

  return
end subroutine scam_setlatlonidx


