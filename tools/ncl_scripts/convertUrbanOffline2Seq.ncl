;
;  Convert the urban point datasets for Mexicocity and Vancouver from off-line to
;  format that datm8 can use.
;
;  Erik Kluzek
;  Jan/23/2008
;  $Id$
;  $HeadURL;
;
  ; Share constants copied from the src/csm_share/shr/shr_const_mod.F90 on 1/18/2007 EBK
  SHR_CONST_TKFRZ  = 273.15;
  SHR_CONST_STEBOL = 5.67e-8;
  SHR_CONST_PSTD   = 101325.0

undef( "eSat" );
function eSat( T [*][*]:float )
;
; Saturation vapor pressure from offline code on 1/18/2007
;
   local td, a0, a1, a2, a3, a4, a5, a6, b0, b1, b2, b3, b4, b5, b6;
begin
   ;--- coefficients for esat over water ---
   a0=6.107799961;
   a1=4.436518521e-01;
   a2=1.428945805e-02;
   a3=2.650648471e-04;
   a4=3.031240396e-06;
   a5=2.034080948e-08;
   a6=6.136820929e-11;

   ;--- coefficients for esat over ice ---
   b0=6.109177956;
   b1=5.034698970e-01;
   b2=1.886013408e-02;
   b3=4.176223716e-04;
   b4=5.824720280e-06;
   b5=4.838803174e-08;
   b6=1.838826904e-10;

   dsize = dimsizes(T);
   td  = new( dsize, "float" );
   sat = new( dsize, "float" );
   td  = T - SHR_CONST_TKFRZ;
   do i = 0, dsize(0)-1
      do j = 0, dsize(1)-1
         if (   td(i,j) .gt.  50. )then
              td(i,j) = 50.;
         else
           if ( td(i,j) .lt. -50. )then
              td(i,j) =-50.;
           end if
         end if
         if( td(i,j) .lt. 0.0 ) then
            sat = 100.0*(b0+td(i,j)*(b1+td(i,j)*(b2+td(i,j)*(b3+td(i,j)*(b4+td(i,j)*(b5+td(i,j)*b6))))));
         else
            sat = 100.0*(a0+td(i,j)*(a1+td(i,j)*(a2+td(i,j)*(a3+td(i,j)*(a4+td(i,j)*(a5+td(i,j)*a6))))));
         end if
      end do
   end do

   return( sat );
end

;
; Main script
;
begin
  ;
  ; Setup datasets to loop over
  ;
  print( "================================================================================================" );
  csmdata  = getenv("CSMDATA");
  if ( ismissing(csmdata) )then
     csmdata  = "/fs/cgd/csm/inputdata";
  end if
  cdate     = systemfunc( "date +%y%m%d" );
  dir1      = csmdata + "/atm/datm7/CLM1PT_data/mexicocityMEX.c"+cdate;
  dir2      = csmdata + "/atm/datm7/CLM1PT_data/vancouverCAN.c"+cdate;
  orig1     = csmdata + "/lnd/clm2/1pt_atmdata/mexicocityMEX.c051206/1993-12.nc";
  orig2     = csmdata + "/lnd/clm2/1pt_atmdata/vancouverCAN.c060920/1992-08.nc";
  system( "mkdir -p "+dir1 );
  system( "mkdir -p "+dir2 );
  filename1 = dir1+"/clm1pt-1993-12.nc";
  filename2 = dir2+"/clm1pt-1992-08.nc";
  date1     = "1993-12-01 02:00:00";
  date2     = "1992-08-12 20:00:00";
  ntime0    = (/   1, 283 /);
  ntimes    = (/ 160, 331 /);
  files     = (/ filename1, filename2 /);
  dates     = (/ date1,     date2     /);
  orig      = (/ orig1,     orig2     /);
  zbot      = (/ 28.4,      28.5      /);
  ldate     = systemfunc( "date" );
  user      = getenv("USER");

  ;
  ; Loop over the two datasets
  ;
  do ifile = 0, 1
     print( "read in file:   "+orig(ifile) );
     nc  = addfile( orig(ifile),  "r" );
     print( "write out file: "+files(ifile) );
     system( "/bin/rm "+files(ifile) );
     nco = addfile( files(ifile), "c" );
     ;
     ; Define dimensions
     ;
     dimnames  = (/ "scalar", "lon", "lat", "time" /);
     dims      = getfiledimsizes( nc );
     dims(3)   = ntimes(ifile);
     is_unlim  = (/ False, False, False, True /);
     filedimdef( nco, dimnames, dims, is_unlim );
     delete( dimnames );
     ;
     ; Define vars and add attributes from original file
     ;
     vars = getfilevarnames( nc );
     i = dimsizes(vars) - 1
     do while ( i .ge. 0 )
        dimlist = getfilevardims( nc, vars(i) )
        filevardef (    nco, vars(i), typeof(nc->$vars(i)$), dimlist );
        filevarattdef ( nco, vars(i), nc->$vars(i)$ );
        delete( dimlist );
        i = i - 1
     end do
     fileattdef ( nco, nc );
     ;
     ; Add attributes
     ;
     nco@conventions = "CF-1.0";
     nco@title       = "CLM single point datm input data";
     nco@history     = "Original data from "+orig(ifile)+"\n Time coords added by "+user+", "+ldate;
     ;
     ; Add time coordinate variable
     ;
     dims = getfiledimsizes(nc);
     filevardef( nco, "time", "float", (/"time"/) )
     nco->time@long_name = "observation time"
     nco->time@units     = " days since "+dates(ifile)
     nco->time@calendar  = "noleap";
     print( "Added attributes and time coordinate variable" );
     ;
     ; Add direct
     ;
     dimlist = getfilevardims( nc, "FSDS" )
     filevardef( nco, "FSDSdir", "float", dimlist );
     nco->FSDSdir@long_name = nc->FSDS@long_name + " direct";
     nco->FSDSdir@units     = nc->FSDS@units;
     nco->FSDSdir@mode      = nc->FSDS@mode;
     print( "Added FSDSdir" )
     ;
     ; Add diffuse
     ;
     filevardef( nco, "FSDSdif", "float", dimlist );
     nco->FSDSdif@long_name = nc->FSDS@long_name + " diffuse";
     nco->FSDSdif@units     = nc->FSDS@units;
     nco->FSDSdif@mode      = nc->FSDS@mode;
     print( "Added FSDSdif" )
     ;
     ; Add bottom layer height
     ;
     filevardef( nco, "ZBOT", "float", dimlist );
     nco->ZBOT@long_name = "observational height";
     nco->ZBOT@units     = "m";
     nco->ZBOT@mode      = nc->FSDS@mode;
     print( "Added ZBOT" )
     ;
     ; Add fill value to data
     ;
     vars           = getfilevarnames( nc );
     do i = 0, dimsizes(vars)-1
        nco->$vars(i)$@_FillValue = -1.;
     end do
     ;
     ; Add data
     ;
     print( "Adding data..." );
     ea   = new( dimsizes(nc->TBOT(0,:,:)), "float" );
     nend = ntimes(ifile)-1+ntime0(ifile);
     do v = 0, dimsizes(vars)-1
        dimlist = getfilevardims( nc, vars(v) )
        if ( dimlist(0) .eq. "time" ) then
           nco->$vars(v)$   = (/ nc->$vars(v)$(ntime0(ifile):nend,:,:) /);
        else
           nco->$vars(v)$   = (/ nc->$vars(v)$ /);
        end if
        delete( dimlist );
     end do
     nco->time(0)        = 0.0
     nco->ZBOT(:,:,:)    = (/ zbot(ifile) /);
     nco->FSDSdif(:,:,:) = 0.0;
     ;
     ; Add new data
     ;
     do i = 0, ntimes(ifile)-1
        n = ntime0(ifile) + i;
        if ( i .gt. 0 )then
           nco->time(i)     = (/nco->time(i-1)/) + 1.0/24.0;
        end if
        nco->FSDSdir(i,:,:) = (/nc->FSDS(n,:,:)/);
        ;
        ; If missing preasure...
        ;
        if ( any(ismissing(nco->PSRF(i,:,:))) )then
           nco->PSRF(i,:,:) = (/ SHR_CONST_PSTD /);
        end if
        ;
        ; If missing LW down...
        ;
        if ( any(ismissing(nco->FLDS(i,:,:))) )then
           T = nco->TBOT(i,:,:);
           if ( .not. any(ismissing(T) ) )then
              ; Prata (1996)
              P  = (/nco->PSRF(i,:,:)/);
              e  = P / ( P*100./( (/nco->RH(i,:,:)/)*eSat( (/T/) )) + 0.378);
              ea1 = 46.5*(0.01*e/(/T/) );
              ea = 1. - (1. + ea1 )*exp( -sqrt(1.2 + 3.*ea1 ) );
              nco->FLDS(i,:,:) = ea * SHR_CONST_STEBOL * (/T/)^4;
           end if
        end if
     end do
     print( nco->time );
     print( "Added data values" );
     ;
     ; Check that all data is valid
     ;
     do v = 0, dimsizes(vars)-1
        if ( any(ismissing(nco->$vars(v)$) ) )then
           print( "Found missing data on "+vars(v) );
           status_exit( 1 );
        end if
     end do
     delete( vars );
  end do

  print( " " );
  print( " " );
  print( "================================================================================================" );
  print( "Successfully converted datasets!!!" );
end
