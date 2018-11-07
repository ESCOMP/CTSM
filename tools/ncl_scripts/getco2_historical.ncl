;
; Take the greenhouse gas file used by CAM for historical (and future) representations of
; greenhouse gases, and convert it to a format that can be used by streams.
; So include domain data for a single point (or latitude bands) that covers the globe, as well
; as CO2 data over those latitude bands. In the process we also discard the other
; greenhouse gases, as the datm can only pass CO2.
;
;  Erik Kluzek
;  Mar/03/2010
;
begin
  ; ===========================================================================================================


  ; ===========================================================================================================
   ;
   ; Setup the namelist query script
   ;
   csmdata  = getenv("CSMDATA");
   clmroot  = getenv("CLM_ROOT");
   hgrid    = getenv("HGRID");          ; Get horizontal grid to use from env variable
   querynml = "bld/queryDefaultNamelist.pl -silent -justvalue ";
   if ( .not. ismissing(csmdata) )then
      querynml = querynml+" -csmdata "+csmdata;
   end if
   if ( ismissing(clmroot) )then
      querynml = "../../"+querynml;
   else
      querynml = clmroot+"/"+querynml;
   end if
   if ( ismissing(hgrid) )then
     hgrid = "lat-bands"
   end if
   ;
   ; Get input Greenhouse gas file and open it
   ;
   filetype  = "mkghg_bndtvghg";
   print( querynml+" -namelist clmexp -var "+filetype+" -options hgrid="+hgrid );
   ghgfile  = systemfunc( querynml+" -namelist clmexp -var "+filetype+" -options hgrid="+hgrid );
   print( "Use "+filetype+" file: "+ghgfile );
   if ( systemfunc("test -f "+ghgfile+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+ghgfile );
      exit
   end if
   ncg = addfile( ghgfile, "r" );

   ;
   ; Get date time-stamp to put on output CO2 file
   ;
   sdate     = systemfunc( "date +%y%m%d" );
   ldate     = systemfunc( "date" );

   sim_yr0 = ncg->date(0) / 10000;
   ntime   = dimsizes( ncg->date );
   sim_yr2 = ncg->date(ntime-1) / 10000;

   sim_yr_rng = "_simyr_"+sim_yr0 + "-" + sim_yr2;
   
   cmip_vers = "_CMIP6_";
   outco2filename = "fco2_datm_"+hgrid+sim_yr_rng+cmip_vers+"c"+sdate+".nc";
   system( "/bin/rm -f "+outco2filename );
   print( "output file: "+outco2filename );
   nco = addfile( outco2filename, "c" );
   ;
   ; Define dimensions
   ;
   if ( hgrid .eq. "lat-bands" )then
      nlat = dimsizes(ncg->lat);
   else
      if ( hgrid .eq. "global" )then
         nlat = 1
      else
         print( "hgrid type can only be global or lat-bands: "+hgrid )
         exit
      end if
   end if
   nlon = 1;
   nv   = 4;
   dimnames = (/ "time", "lat", "lon", "nv", "bounds" /);
   dsizes   = (/ ntime, nlat,  nlon, nv, 2 /);
   is_unlim = (/ True, False, False, False, False /);
   filedimdef( nco, dimnames, dsizes, is_unlim );
   ;
   ; Define variables
   ;
   vars = (/ "lonc", "latc", "lonv", "latv", "mask", "frac", "area", "CO2" /);
   units= (/ "degrees_east", "degrees_north", "degree_east", "degrees_north", "unitless", "unitless",          "radians^2",        "ppmv" /);
   lname= (/ "Longitude of grid cell center", "Latitude of grid cell center", "Longitudes of grid cell vertices", "Latitudes of grid cell vertices", "Mask of active cells: 1=active", "Fraction of grid cell that is active", "Area of grid cell", "CO2 concentration" /);
   print( "Define variables: "+vars );
   do i= 0, dimsizes(vars)-1
      if ( vars(i) .eq. "lonv" .or. vars(i) .eq. "latv" )then
         filevardef ( nco, vars(i), "double",  (/ "lat", "lon", "nv" /) );
      else
         if ( vars(i) .eq. "CO2" )then
            filevardef ( nco, vars(i),  "float",   (/ "time", "lat", "lon" /) );
            nco->$vars(i)$@coordinate  = "latc lonc time";
         else
            filevardef ( nco, vars(i), "double",  (/ "lat", "lon" /) );
         end if
      end if
      nco->$vars(i)$@units = units(i);
      nco->$vars(i)$@lname = lname(i);
   end do
   filevardef ( nco, "time",      "float",    (/ "time" /) );
   filevardef ( nco, "time_bnds", "float",    (/ "time", "bounds" /) );
   filevardef ( nco, "date",      "integer",  (/ "time" /) );
   varstatic = (/ "mask", "frac", "area" /);
   do i = 0, dimsizes(varstatic)-1
      nco->$varstatic(i)$@coordinate  = "latc lonc";
   end do
   nco->lonc@bounds      = "lonv";
   nco->latc@bounds      = "latv";
   ;
   ; Add attributes
   ;
   fileattdef ( nco, ncg );
   nco@history  = ldate+": Convert by getco2_historical.ncl";
   nco@source   = "Convert from:"+ghgfile;
   nco@Version  = systemfunc( "git describe" );
   filevarattdef( nco, "time", ncg->time );
   filevarattdef( nco, "date", ncg->date );
   nco->time_bnds@long_name = nco->time@long_name;
   nco->time_bnds@units     = nco->time@units;
   nco->time_bnds@calendar  = nco->time@calendar;
   ;
   ; Set static variables
   ;
   pi                      = 3.14159265358979323846d00;
   nco->mask               =   1;
   nco->frac               =   1.0;
   if ( nlat .gt. 1 )then
      nco->latc               = (/ ncg->lat/);
   else
      nco->latc               = (/ 0.0d00 /);
   end if
   nco->latv(nlat-1,0,0)   =  90.0d00;
   nco->latv(nlat-1,0,3)   =  90.0d00;
   if ( nlat .gt. 1 )then
     nco->latv(0:nlat-2,0,0) =  ( (/ ncg->lat(0:nlat-2) /) + (/ncg->lat(1:nlat-1) /) )*0.5d00
     nco->latv(0:nlat-2,0,3) =  (/ nco->latv(0:nlat-2,0,0) /);
     nco->latv(1:nlat-1,0,1) =  (/ nco->latv(0:nlat-2,0,0) /);
     nco->latv(1:nlat-1,0,2) =  (/ nco->latv(1:nlat-1,0,1) /);
   end if
   nco->latv(0,0,1)        = -90.0d00;
   nco->latv(0,0,2)        = -90.0d00;
   nco->lonv(:,0,0)        =   0.0d00;
   nco->lonv(:,0,3)        =   0.0d00;
   nco->lonc               = 180.0d00;
   nco->lonv(:,0,1)        = 360.0d00;
   nco->lonv(:,0,2)        = 360.0d00;
   clkws = gc_clkwise( nco->latv, nco->lonv );
   if ( any(clkws .eq. False) )then
      print( "Some varticies are NOT clockwise" );
      exit
   end if
   ; EBK -- NOTE The NCL function wasn't giving me the correct answer so I used the mathmatical expression
   ;nco->area = dble2flt( gc_qarea( nco->latv, nco->lonv ) );
   conv2rad = pi/180.0d00
   nco->area(:,0) = 2.0d00*pi*abs( sin((/nco->latv(:,0,0)/)*conv2rad) - sin((/nco->latv(:,0,1)/)*conv2rad) );
   if ( abs(sum(nco->area) - 4.0d00*pi) .gt. 1.d-14 )then
      print( "Area of globe does not sum to 4*pi as expected" );
      exit
   end if
   ;
   ; Time and date
   ;
   nco->date = (/ ncg->date /);
   nco->time = (/ ncg->time /);
   nco->date@comment = "This variable is NOT used when read by datm, the time coordinate is used";
   ;
   ; CO2
   ;
   print( "Copy CO2 for "+ntime+" time samples of data" );
   if ( nlat .gt. 1 )then
      do y = 0, nlat-1
        print( "latitude: "+ nco->latc(y,0) );
        nco->CO2(:,y,0) = (/ ncg->CO2_LBC(:,y) /) * 1.e6;
      end do
   else
      ; make sure all latitudes on file are the same for each time
      do itime = 0, ntime-1
         if ( max(ncg->CO2_LBC(itime,:)) .ne. min(ncg->CO2_LBC(itime,:)) )then
            print( "Global average, but latitudes are NOT constant" );
            exit
         end if
      end do
      nco->CO2(:,0,0) = (/ ncg->CO2_LBC(:,0) /) * 1.e6;
   end if
   print( "Average Global First CO2 ppmv value: Date="+nco->date(0)+" CO2="+avg(nco->CO2(0,:,0)     ) );
   print( "Average Global Last  CO2 ppmv value: Date="+nco->date(ntime-1)+" CO2="+avg(nco->CO2(ntime-1,:,0)) );

   print( "================================================================================================" );
   print( "Successfully created output historical CO2 file: "+outco2filename);

end
