module getnetcdfdata

! Description:
!   Routines for extracting a column from a netcdf file
!
! Author: 
!   
! Modules Used:
!
  use abortutils,   only: endrun

  implicit none
  include 'netcdf.inc'

  private
!
! Public Methods:
!
  public :: getncdata,setlatlonidx,handle_err
!
! !PRIVATE METHODS
!
  interface getncdata
     module procedure getncdata_real_1d
     module procedure getncdata_int_1d
     module procedure getncdata_real_scalar_2d
     module procedure getncdata_int_scalar_2d
     module procedure getncdata_real_scalar
     module procedure getncdata_int_scalar
  end interface

contains
subroutine setlatlonidx(ncid ,targetlat ,targetlon, closelat, closelon, &
                        closelatidx, closelonidx)
!------------------------------------------------------------------------
! File: setlatlon.F 
! Author: John Truesdale (jet@ucar.edu)
! $Id$
!
!------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use abortutils,   only: endrun
  implicit none

!------------------------------locals-----------------------------------
  
  integer, intent(in):: NCID
  real(r8), intent(in):: targetlat,targetlon
  real(r8), intent(out):: closelat,closelon
  integer, intent(out) :: closelatidx,closelonidx
  real(r8),allocatable:: lats(:),lons(:)
  real(r8) prev, next             !
  integer  i                ! 
  integer  status
  integer  nlat_dimID, nlat, lat_varID 
  integer  nlon_dimID, nlon, lon_varID
  integer strt(2),cnt(2)
  logical lnddata
!
! Check mode: double or single precision.
! 
  
  lnddata = .false.
  STATUS =  NF_INQ_DIMID ( NCID, 'lat', nlat_dimID )
  if ( STATUS .NE.  NF_NOERR )  then
     STATUS =  NF_INQ_DIMID ( NCID, 'lsmlat', nlat_dimID )
     if ( STATUS .NE.  NF_NOERR )  then
        write( 6,* )'ERROR - setlatlon.F:', &
             'Cant get variable dim for lat or lsmlat'
        STATUS = NF_CLOSE( NCID )
        call endrun
     else
        lnddata = .true.
     endif
  endif

  STATUS =  NF_INQ_DIMID ( NCID, 'lon', nlon_dimID )
  if ( STATUS .NE.  NF_NOERR )  then
     STATUS =  NF_INQ_DIMID ( NCID, 'lsmlon', nlon_dimID )
     if ( STATUS .NE.  NF_NOERR )  then
        write( 6,* )'ERROR - setlatlon.F:', &
             'Cant get variable dim for lon or lsmlon'
        STATUS = NF_CLOSE( NCID )
        call endrun
     endif
  endif

  STATUS = NF_INQ_DIMLEN( NCID, nlat_dimID, nlat )
  STATUS = NF_INQ_DIMLEN( NCID, nlon_dimID, nlon )

  allocate(lats(nlat),lons(nlon))

  STATUS = NF_INQ_VARID( NCID, 'lat', lat_varID )
  if ( STATUS .NE. NF_NOERR ) then
     STATUS = NF_INQ_VARID( NCID, 'LATIXY', lat_varID )
     if ( STATUS .NE. NF_NOERR ) then
        write( 6,* )'ERROR - setlatlon.F:', &
             'Cant get variable ID for lat or LATIXY'
        STATUS = NF_CLOSE( NCID )
        return
     endif
  endif

  if (lnddata) then
     strt(1)=1
     strt(2)=1
     cnt(1)=1
     cnt(2)=nlat
  else
     strt(1)=1
     strt(2)=1
     cnt(1)=nlat
     cnt(2)=1
  endif

  STATUS = nf_get_vara_double (ncid, lat_varID, strt, cnt, lats)
  if (STATUS/=NF_NOERR) then
     write(6,*)'Error: setlatlonidx: error reading lat varid =', lat_varID
     call endrun
  end if

  STATUS = NF_INQ_VARID( NCID, 'lon', lon_varID )
  if ( STATUS .NE. NF_NOERR ) then
     STATUS = NF_INQ_VARID( NCID, 'LONGXY', lon_varID )
     if ( STATUS .NE. NF_NOERR ) then
        write( 6,* )'ERROR - setlatlon.F:', &
             'Cant get variable ID for lon'
        STATUS = NF_CLOSE( NCID )
        return
     endif
  endif

  if (lnddata) then
     strt(1)=1
     strt(2)=1
     cnt(1)=nlon
     cnt(2)=1
  else
     strt(1)=1
     strt(2)=1
     cnt(1)=nlon
     cnt(2)=1
  endif

  STATUS = nf_get_vara_double (ncid, lon_varID, strt, cnt, lons)
  if (STATUS/=NF_NOERR) then
     write(6,*)'Error: setlatlonidx: error reading lon varid =', lon_varID
     call endrun
  end if

  closelatidx = 1
  do i = 1, nlat
     next = lats(i)
     prev = lats(closelatidx)
     if ( abs(targetlat - next) .LT.  &
          abs(targetlat - prev) ) then
        closelatidx = i
     endif
  enddo

  closelonidx = 1
  do i = 1, nlon
     if ( lons(i) .le. 180.0 ) then 
        next = lons(i)       
     else
        next = lons(i) - 360.0
     endif
     if ( lons(closelonidx) .le. 180.0 ) then 
        prev = lons(closelonidx)
     else
        prev = lons(closelonidx) - 360.0
     endif


     if ( abs(targetlon - next) .LT. &
          abs(targetlon - prev) ) then
        closelonidx = i
     endif
  enddo

  closelon=lons(closelonidx)
  closelat=lats(closelatidx)

  return
end subroutine setlatlonidx

subroutine getncdata_int_scalar_2d ( NCID, camlat, camlon, TimeIdx, &
                       varName, outData, STATUS )
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   implicit none
   integer     NCID, TimeIdx ,STATUS
   real(r8)   camlat,camlon
   character     varName*(*)
   integer outData(1,1)
   integer :: temp_array(1)

   call getncdata_int_1d  ( NCID, camlat, camlon, TimeIdx, &
                       varName, temp_array, STATUS )
   outData(1,1) = temp_array(1)
 end subroutine getncdata_int_scalar_2d

subroutine getncdata_real_scalar_2d ( NCID, camlat, camlon, TimeIdx, &
                       varName, outData, STATUS )
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   implicit none
   integer     NCID, TimeIdx ,STATUS
   character     varName*(*)
   real(r8)   camlat,camlon
   real(r8) outData(1,1)
   real(r8) :: temp_array(1)

   call getncdata_real_1d  ( NCID, camlat, camlon, TimeIdx, &
                       varName, temp_array, STATUS )
   outData(1,1) = temp_array(1)
 end subroutine getncdata_real_scalar_2d

subroutine getncdata_real_scalar ( NCID, camlat, camlon, TimeIdx, &
                       varName, outData, STATUS )
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   implicit none
#include <netcdf.inc>
   integer     NCID, TimeIdx ,STATUS
   real(r8)   camlat,camlon
   character     varName*(*)
   real(r8) outData
   real(r8) :: temp_array(1)

   call getncdata_real_1d  ( NCID, camlat, camlon, TimeIdx, &
                       varName, temp_array, STATUS )
   if (STATUS.eq.NF_NOERR) outData = temp_array(1)
 end subroutine getncdata_real_scalar

subroutine getncdata_int_scalar ( NCID, camlat, camlon, TimeIdx, &
                       varName, outData, STATUS )
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   implicit none
   integer     NCID, TimeIdx ,STATUS
   real(r8)   camlat,camlon
   character     varName*(*)
   integer outData
   integer :: temp_array(1)

   call getncdata_int_1d  ( NCID, camlat, camlon, TimeIdx, &
                       varName, temp_array, STATUS )
   outData = temp_array(1)
 end subroutine getncdata_int_scalar


subroutine getncdata_real_1d ( NCID, camlat, camlon, TimeIdx, &
                       varName, outData, STATUS )
!------------------------------------------------------------------------
! File: getncdata.F 
! Author: John Truesdale (jet@ucar.edu) 
! $Id$
!
!------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>

!-----------------------------------------------------------------------
! getncdata: extract the column of data for a particular
!     latitude, longitude, and timepoint from a netCDF file.
!
! INPUTS:
!
!	NCID		           NetCDF ID
!
!	camlat, camlon    target lat, lon of  column to extract.
!
!	TimeIdx		Time slice to extract from.  ( If a multi-timed
!				field.
!
!	varName	Desired field
!
!	outData		The pre-allocated destination where NCVGT will
!				store the data.
!
!       STATUS           Error code returned in argument list:
!                           zero if successful, non-zero otherwise
!
!-----------------------------------------------------------------------	
!
   real(r8),  intent (in)  :: camlat,camlon
   integer,   intent (in)  :: NCID, TimeIdx
   character, intent(in)   :: varName*(*)
   integer,   intent (out) :: STATUS
   real(r8),  intent(out)  :: outData(:)
!
!------------------------Local-----------------------------------------	
!
   real(r8) closelat,closelon
   integer     latidx, lonidx
   integer     totalloc
   integer     varID, var_ndims, dim_size, dims_set, i,var_type,j,varidtmp
   integer     var_dimIDs( NF_MAX_VAR_DIMS )
   integer     start(5 ), count(5 )
   integer     startloc(1), countloc(1)
   integer     mapnumland(300,300)
   integer , allocatable :: landvecixy(:)
   integer , allocatable :: landvecjxy(:)

   character     dim_name*( NF_MAX_NAME )

   logical usable_var
   logical use_nf_real,use_nf_int


   call setlatlonidx(ncid,camlat,camlon,closelat,closelon,latidx,lonidx)
!
! Check mode: double or single precision.
! 

#if USE_4BYTE_REAL
   use_nf_real = .true.
#else
   use_nf_real = .false.
#endif

!
! Get var ID.  Check to make sure it's there.
!
   STATUS = NF_INQ_VARID( NCID, varName, varID )

   if ( STATUS .NE. NF_NOERR )  then
      return
   endif
!
! Check the var variable's information with what we are expecting
! it to be.
!
   STATUS = NF_INQ_VARNDIMS( NCID, varID, var_ndims )
   if ( var_ndims .GT. 4 ) then
      write( 6,* ) 'ERROR - getncdata.F: The input var',varName, &
         'has', var_ndims, 'dimensions'
      STATUS = -1
      return
   endif

   STATUS =  NF_INQ_VARTYPE(NCID, varID, var_type)
   if ( var_type .NE. NF_FLOAT .and. var_type .NE. NF_DOUBLE .and. &
        var_type .NE. NF_INT ) then
      write( 6,* ) 'ERROR - getncdata.F: The input var',varName, &
         'has unknown type', var_type
      STATUS = -1
      return
   endif

   if (var_type==NF_INT) then
     use_nf_int = .true.
   else
     use_nf_int = .false.
   end if
   
!
!     surface variables
!
   if ( var_ndims .EQ. 0 ) then
      if (use_nf_int) then
         STATUS = NF_GET_VAR_INT( NCID, varID, outData )
      else
         if (use_nf_real) then
            STATUS = NF_GET_VAR_REAL( NCID, varID, outData )
         else
            STATUS = NF_GET_VAR_DOUBLE( NCID, varID, outData )
         endif
      endif
      call handle_err( STATUS )
      return
   endif

   STATUS = NF_INQ_VARDIMID( NCID, varID, var_dimIDs )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* ) 'ERROR - getncdata.F:Cant get dimension IDs for', varName
      call handle_err( STATUS )
      return
   endif
!     
!     Initialize the start and count arrays 
!     
   dims_set = 0
   do i =  var_ndims, 1, -1

      usable_var = .false.
      STATUS = NF_INQ_DIMNAME( NCID, var_dimIDs( i ), dim_name )
      if ( STATUS .NE. NF_NOERR ) then
         write(6,*) 'Error: getncdata.F - can''t get dim name', &
            'var_ndims = ', var_ndims, ' i = ',i
         call handle_err( STATUS )
         return
      endif
      if ( dim_name .EQ. 'lat' ) then
         start( i ) =  latIdx
         count( i ) = 1           ! Extract a single value 
         dims_set = dims_set + 1
         usable_var = .true.
      endif
      if ( dim_name .EQ. 'lsmlat' ) then
         start( i ) =  latIdx
         count( i ) = 1           ! Extract a single value 
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lon' ) then
         start( i ) = lonIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif
      if ( dim_name .EQ. 'lsmlon' ) then
         start( i ) = lonIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lev' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim len'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lsmpft' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim lsmpft'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'nlevsoi' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim nlevsoi'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'nlevsno' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim nlevso'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'nlevtot' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim nlevtot'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'numpatch' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim numpatch'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'maxpatch' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim maxpatch'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

! get correct landpoint out of numland array,  to do this first extract
! arrays landvecixy and landvecjxy and then build a tmp array for mapping
! lat lon to numland point(mapnumland array).

      if ( dim_name .EQ. 'numland' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim numland'
            call handle_err( STATUS )
            return
         endif
!
! Get var IDs for landvecixy and landvecjxy.
!
         STATUS = NF_INQ_VARID( NCID, 'LANDVECIXY', varidtmp )
         
         if ( STATUS .NE. NF_NOERR )  then
            write(6,*) 'Error: getncdata.F - can''t get var LANDVECIXY'
            call handle_err( STATUS )
         endif
         startloc( 1 ) = 1
         countloc( 1 ) = dim_size ! Extract all levels
         allocate (landvecixy(dim_size))
         STATUS = NF_GET_VARA_INT( NCID, varidtmp, startloc, countloc, landvecixy )

         STATUS = NF_INQ_VARID( NCID, 'LANDVECJXY', varidtmp )
         
         if ( STATUS .NE. NF_NOERR )  then
            write(6,*) 'Error: getncdata.F - can''t get var LANDVECJXY'
            call handle_err( STATUS )
         endif
         allocate (landvecjxy(dim_size))
         STATUS = NF_GET_VARA_INT( NCID, varidtmp, startloc, countloc, landvecjxy )
         mapnumland(:,:)=-1
         do j = 1,dim_size
            mapnumland(landvecixy(j),landvecjxy(j))=j
         end do
         if (mapnumland(lonIdx,latIdx)== -1) then
            write(6,*)'Error in getncdata mapping numland point to model lat and lon'
            call endrun
         end if
              
         start( i ) = mapnumland(lonIdx,latIdx)
         count( i ) = 1 ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
         deallocate(landvecixy)
         deallocate(landvecjxy)

      endif


      if ( dim_name .EQ. 'msub' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim msub'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'ilev' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim len'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'time' .OR. dim_name .EQ. 'tsec'  ) then
         start( i ) = TimeIdx
         count( i ) = 1      ! Extract a single value 
         dims_set = dims_set + 1   
         usable_var = .true.
      endif

      if ( usable_var .EQV. .false. ) then
         write( 6,* )'ERROR - getncdata.F: The input var ', &
            varName, ' has an unusable dimension ', dim_name
         STATUS = -1
      endif
   end do
   if ( dims_set .NE. var_ndims ) then
      write( 6,* )'ERROR - getncdata.F: Could not find all the', &
         'dimensions for input var ', varName
      write( 6,* )'Found ',dims_set, ' of ',var_ndims
      STATUS = -1
   endif

   if (use_nf_int) then
      STATUS = NF_GET_VARA_INT( NCID, varID, start, count, outData )
   else
      if (use_nf_real) then
         STATUS = NF_GET_VARA_REAL( NCID, varID, start, count, outData )
      else   
         STATUS = NF_GET_VARA_DOUBLE( NCID, varID, start, count, outData )
      endif
   end if

   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - getncdata.F: Could not get data for input var ', varName
      call handle_err( STATUS )
   endif

   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - getncdata.F: Could not get data for input var ', varName
      call handle_err( STATUS )
   endif

   return
end subroutine getncdata_real_1d

subroutine getncdata_int_1d ( NCID, camlat, camlon, TimeIdx, &
                       varName, outData, STATUS )
   use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
   implicit none
!-----------------------------------------------------------------------
#include <netcdf.inc>


!-----------------------------------------------------------------------
! getncdata: extract the column of data for a particular
!     latitude, longitude, and timepoint from a netCDF file.
!
! INPUTS:
!
!	NCID		           NetCDF ID
!
!	latIdx, lonIdx    lat, lon of  column to extract.
!
!	TimeIdx		Time slice to extract from.  ( If a multi-timed
!				field.
!
!	varName	Desired field
!
!	outData		The pre-allocated destination where NCVGT will
!				store the data.
!
!       STATUS           Error code returned in argument list:
!                           zero if successful, non-zero otherwise
!
!-----------------------------------------------------------------------	
!
!-----------------------------------------------------------------------	
!
   real(r8),  intent (in)  :: camlat,camlon
   integer,   intent (in)  :: NCID, TimeIdx
   character, intent(in)   :: varName*(*)
   integer,   intent (out) :: STATUS
   integer,   intent(out)  :: outData(:)
!
!------------------------Local-----------------------------------------	
!
   real(r8) closelat,closelon
   integer  latidx, lonidx
   integer  totalloc
   integer  varID, var_ndims, dim_size, dims_set, i,var_type,j,varidtmp
   integer  var_dimIDs( NF_MAX_VAR_DIMS )
   integer  start(5 ), count(5 )
   integer  startloc(1), countloc(1)
   integer  mapnumland(300,300)
   integer , allocatable :: landvecixy(:)
   integer , allocatable :: landvecjxy(:)

   character     dim_name*( NF_MAX_NAME )

   logical usable_var
   logical use_nf_real,use_nf_int


   call setlatlonidx(ncid,camlat,camlon,closelat,closelon,latidx,lonidx)

!
! Check mode: double or single precision.
! 

#if USE_4BYTE_REAL
   use_nf_real = .true.
#else
   use_nf_real = .false.
#endif

!
! Get var ID.  Check to make sure it's there.
!
   STATUS = NF_INQ_VARID( NCID, varName, varID )

   if ( STATUS .NE. NF_NOERR )  then
      return
   endif
!
! Check the var variable's information with what we are expecting
! it to be.
!
   STATUS = NF_INQ_VARNDIMS( NCID, varID, var_ndims )
   if ( var_ndims .GT. 4 ) then
      write( 6,* ) 'ERROR - getncdata.F: The input var',varName, &
         'has', var_ndims, 'dimensions'
      STATUS = -1
      return
   endif

   STATUS =  NF_INQ_VARTYPE(NCID, varID, var_type)
   if ( var_type .NE. NF_FLOAT .and. var_type .NE. NF_DOUBLE .and. &
        var_type .NE. NF_INT ) then
      write( 6,* ) 'ERROR - getncdata.F: The input var',varName, &
         'has unknown type', var_type
      STATUS = -1
      return
   endif

   if (var_type==NF_INT) then
      use_nf_int = .true.
   else
      use_nf_int = .false.
   end if

!
!     surface variables
!
   if ( var_ndims .EQ. 0 ) then
      if (use_nf_int) then
         STATUS = NF_GET_VAR_INT( NCID, varID, outData )
      else
         if (use_nf_real) then
           write(6,*)'Warning: getncdata is truncating real data into integer argument'
            STATUS = NF_GET_VAR_REAL( NCID, varID, outData )
         else
           write(6,*)'Warning: getncdata is truncating double data into integer argument'
            STATUS = NF_GET_VAR_DOUBLE( NCID, varID, outData )
         endif
      endif
      call handle_err( STATUS )
      return
   endif

   STATUS = NF_INQ_VARDIMID( NCID, varID, var_dimIDs )
   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* ) 'ERROR - getncdata.F:Cant get dimension IDs for', varName
      call handle_err( STATUS )
      return
   endif
!     
!     Initialize the start and count arrays 
!     
   dims_set = 0
   do i =  var_ndims, 1, -1

      usable_var = .false.
      STATUS = NF_INQ_DIMNAME( NCID, var_dimIDs( i ), dim_name )
      if ( STATUS .NE. NF_NOERR ) then
         write(6,*) 'Error: getncdata.F - can''t get dim name', &
            'var_ndims = ', var_ndims, ' i = ',i
         call handle_err( STATUS )
         return
      endif
      if ( dim_name .EQ. 'lat' ) then
         start( i ) =  latIdx
         count( i ) = 1           ! Extract a single value 
         dims_set = dims_set + 1
         usable_var = .true.
      endif
      if ( dim_name .EQ. 'lsmlat' ) then
         start( i ) =  latIdx
         count( i ) = 1           ! Extract a single value 
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lon' ) then
         start( i ) = lonIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif
      if ( dim_name .EQ. 'lsmlon' ) then
         start( i ) = lonIdx
         count( i ) = 1           ! Extract a single value
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lev' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim len'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'lsmpft' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim lsmpft'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'nlevsoi' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim nlevsoi'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'nlevsno' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim nlevso'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'nlevtot' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim nlevtot'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'numpatch' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim numpatch'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'maxpatch' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim maxpatch'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

! get correct landpoint out of numland array,  to do this first extract
! arrays landvecixy and landvecjxy and then build a tmp array for mapping
! lat lon to numland point(mapnumland array).

      if ( dim_name .EQ. 'numland' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim numland'
            call handle_err( STATUS )
            return
         endif
!
! Get var IDs for landvecixy and landvecjxy.
!
         STATUS = NF_INQ_VARID( NCID, 'LANDVECIXY', varidtmp )
         
         if ( STATUS .NE. NF_NOERR )  then
            write(6,*) 'Error: getncdata.F - can''t get var LANDVECIXY'
            call handle_err( STATUS )
         endif
         startloc( 1 ) = 1
         countloc( 1 ) = dim_size ! Extract all levels
         allocate (landvecixy(dim_size))
         STATUS = NF_GET_VARA_INT( NCID, varidtmp, startloc, countloc, landvecixy )

         STATUS = NF_INQ_VARID( NCID, 'LANDVECJXY', varidtmp )
         
         if ( STATUS .NE. NF_NOERR )  then
            write(6,*) 'Error: getncdata.F - can''t get var LANDVECJXY'
            call handle_err( STATUS )
         endif
         allocate (landvecjxy(dim_size))
         STATUS = NF_GET_VARA_INT( NCID, varidtmp, startloc, countloc, landvecjxy )
         mapnumland(:,:)=-1
         do j = 1,dim_size
            mapnumland(landvecixy(j),landvecjxy(j))=j
         end do
         if (mapnumland(lonIdx,latIdx)== -1) then
            write(6,*)'Error in getncdata mapping numland point to model lat and lon'
            call endrun
         end if
              
         start( i ) = mapnumland(lonIdx,latIdx)
         count( i ) = 1 ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
         deallocate(landvecixy)
         deallocate(landvecjxy)

      endif


      if ( dim_name .EQ. 'msub' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim msub'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'ilev' ) then
         STATUS = NF_INQ_DIMLEN( NCID, var_dimIDs( i ), dim_size )
         if ( STATUS .NE. NF_NOERR ) then
            write(6,*) 'Error: getncdata.F - can''t get dim len'
            call handle_err( STATUS )
            return
         endif
         start( i ) = 1
         count( i ) = dim_size ! Extract all levels
         dims_set = dims_set + 1
         usable_var = .true.
      endif

      if ( dim_name .EQ. 'time' .OR. dim_name .EQ. 'tsec'  ) then
         start( i ) = TimeIdx
         count( i ) = 1      ! Extract a single value 
         dims_set = dims_set + 1   
         usable_var = .true.
      endif

      if ( usable_var .EQV. .false. ) then
         write( 6,* )'ERROR - getncdata.F: The input var ', &
            varName, ' has an unusable dimension ', dim_name
         STATUS = -1
      endif
   end do
   if ( dims_set .NE. var_ndims ) then
      write( 6,* )'ERROR - getncdata.F: Could not find all the', &
         'dimensions for input var ', varName
      write( 6,* )'Found ',dims_set, ' of ',var_ndims
      STATUS = -1
   endif

   if (use_nf_int) then
      STATUS = NF_GET_VARA_INT( NCID, varID, start, count, outData )
   else
      if (use_nf_real) then
         STATUS = NF_GET_VARA_REAL( NCID, varID, start, count, outData )
      else   
         STATUS = NF_GET_VARA_DOUBLE( NCID, varID, start, count, outData )
      endif
   end if

   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - getncdata.F: Could not get data for input var ', varName
      call handle_err( STATUS )
   endif

   if ( STATUS .NE. NF_NOERR ) then
      write( 6,* )'ERROR - getncdata.F: Could not get data for input var ', varName
      call handle_err( STATUS )
   endif

   return
end subroutine getncdata_int_1d

subroutine handle_err( STATUS )

   implicit none

#include <netcdf.inc>

   integer STATUS

   if ( STATUS .NE. NF_NOERR ) then
      print *, 'NETCDF ERROR: ', NF_STRERROR( STATUS )
   endif
end subroutine handle_err


end module getnetcdfdata
