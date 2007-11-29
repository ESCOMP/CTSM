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
  
  implicit none
#include <netcdf.inc>

!------------------------------parameters-----------------------------------

  integer, intent(in):: NCID
  real(r8), intent(in):: targetlat,targetlon
  real(r8), intent(out):: closelat,closelon
  integer, intent(out) :: closelatidx,closelonidx

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
  
  STATUS =  NF_INQ_DIMID ( NCID, 'lat', nlat_dimID )
  if ( STATUS .NE.  NF_NOERR )  then
     STATUS =  NF_INQ_DIMID ( NCID, 'lsmlat', nlat_dimID )
     if ( STATUS .NE.  NF_NOERR )  then
        write(iulog,* )'ERROR - setlatlon.F:', &
             'Cant get variable dim for lat or lsmlat'
        STATUS = NF_CLOSE( NCID )
        call endrun
     endif
  endif

  STATUS =  NF_INQ_DIMID ( NCID, 'lon', nlon_dimID )
  if ( STATUS .NE.  NF_NOERR )  then
     STATUS =  NF_INQ_DIMID ( NCID, 'lsmlon', nlon_dimID )
     if ( STATUS .NE.  NF_NOERR )  then
        write(iulog,* )'ERROR - setlatlon.F:', &
             'Cant get variable dim for lon or lsmlon'
        STATUS = NF_CLOSE( NCID )
        call endrun
     endif
  endif

  STATUS = NF_INQ_DIMLEN( NCID, nlat_dimID, nlat )
  STATUS = NF_INQ_DIMLEN( NCID, nlon_dimID, nlon )

  allocate(lats(nlat),lons(nlon),poslons(nlon))

  strt(1)=1
  strt(2)=1
  cnt(1)=nlat
  cnt(2)=1

  STATUS = NF_INQ_VARID( NCID, 'lat', lat_varID )
  if ( STATUS .NE. NF_NOERR ) then
     STATUS = NF_INQ_VARID( NCID, 'LATIXY', lat_varID )
     if ( STATUS .NE. NF_NOERR ) then
        write(iulog,* )'ERROR - setlatlon.F:', &
             'Cant get variable ID for lat or LATIXY'
        STATUS = NF_CLOSE( NCID )
        return
     else
        strt(1)=nlon
        strt(2)=1
        cnt(1)=1
        cnt(2)=nlat
     endif
  endif

  STATUS = nf_get_vara_double (ncid, lat_varID, strt, cnt, lats)
  if (STATUS/=NF_NOERR) then
     write(iulog,*)'Error: scam_setlatlonidx: error reading lat varid =', lat_varID
     call endrun
  end if

  STATUS = NF_INQ_VARID( NCID, 'lon', lon_varID )
  if ( STATUS .NE. NF_NOERR ) then
     STATUS = NF_INQ_VARID( NCID, 'LONGXY', lon_varID )
     if ( STATUS .NE. NF_NOERR ) then
        write(iulog,* )'ERROR - setlatlon.F:', &
             'Cant get variable ID for lon'
        STATUS = NF_CLOSE( NCID )
        return
     endif
  endif

  strt(1)=1
  strt(2)=1
  cnt(1)=nlon
  cnt(2)=1

  STATUS = nf_get_vara_double (ncid, lon_varID, strt, cnt, lons)
  if (STATUS/=NF_NOERR) then
     write(iulog,*)'Error: scam_setlatlonidx: error reading lon varid =', lon_varID
     call endrun
  end if

!!$ convert lons array and targetlon to 0,360

  poslons=mod(lons+360._r8,360._r8)
  postargetlon=mod(targetlon+360._r8,360._r8)

!!$ find index of value closest to 0.

  closelonidx=(MINLOC(abs(poslons-postargetlon),dim=1))
  closelatidx=(MINLOC(abs(lats-targetlat),dim=1))
  closelon=lons(closelonidx)
  closelat=lats(closelatidx)

  return
end subroutine scam_setlatlonidx


