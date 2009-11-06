;
; Extract a given year from the 20th Century transient Nitrogen file and put
; it in a constant ndepdat file.
;
; Erik Kluzek
; Nov/04/2009
; $Id$
; $HeadURL;
;
begin
   ; ===========================================================================================================

   res        = getenv("RES");        ; Get output resolution from env variable
   sim_yr     = getenv("SIM_YR");     ; Get the simulation year from env variable
   sim_yr_rng = getenv("SIM_YR_RNG"); ; Get the simulation year range from env variable

   if ( ismissing(res) )then
      res = "1.9x2.5";
   end if
   if ( ismissing(sim_yr) )then
      sim_yr = "1850";
   end if
   if ( ismissing(sim_yr_rng) )then
      sim_yr_rng = "1850-2000";
   end if

   ; ===========================================================================================================
   load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

   ;
   ; Setup the namelist query script
   ;
   csmdata  = getenv("CSMDATA");
   clmroot  = getenv("CLM_ROOT");
   querynml = "bld/queryDefaultNamelist.pl -silent -justvalue ";
   if ( .not. ismissing(csmdata) )then
      querynml = querynml+" -csmdata "+csmdata;
   end if
   if ( ismissing(clmroot) )then
      querynml = "../../"+querynml;
   else
      querynml = clmroot+"/models/lnd/clm/"+querynml;
   end if
   ;
   ; Use resolution to get input filenames and open them
   ;
   filetype = "fndepdyn";
   filename = systemfunc( querynml+" -res "+res+" -var "+filetype+" -options bgc=cn,sim_year_range="+sim_yr_rng );
   print( "Use "+filetype+" file: "+filename+" for sim year range ="+sim_yr_rng );
   if ( systemfunc("test -f "+filename+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+filename );
      exit
   end if
   nc = addfile( filename, "r" );

   ;
   ; Get date time-stamp to put on output file
   ;
   sdate     = systemfunc( "date +%y%m%d" );
   ldate     = systemfunc( "date" );

   outfilename = "ndep_clm_simyr"+sim_yr+"_"+res+"_c"+sdate+".nc";
   system( "/bin/rm -f "+outfilename );
   print( "output file: "+outfilename );
   nco = addfile( outfilename, "c" );
   ;
   ; Define dimensions
   ;
   dimnames = (/ "lat", "lon" /);
   nlat     = dimsizes( nc->lat );
   nlon     = dimsizes( nc->lon );
   dsizes   = (/ nlat,  nlon /);
   is_unlim = (/ False, False /);
   filedimdef( nco, dimnames, dsizes, is_unlim );
   ;
   ; Define variables
   ;
   vars  = (/ "NDEP_NHx_year", "NDEP_NOy_year", "NDEP_year", "lat", "lon" /);
   do i = 0, dimsizes(vars)-1
      print( "Add variable: "+vars(i) );
      if ( vars(i) .ne. "lat" .and. vars(i) .ne. "lon" )then
         dimlist = dimnames;
      else
         dimlist = getfilevardims( nc, vars(i) )
      end if
      filevardef ( nco, vars(i), typeof(nc->$vars(i)$), dimlist );
      filevarattdef ( nco, vars(i), nc->$vars(i)$ );
      delete( dimlist );
   end do
   ;
   ; Add global attributes
   ;
   print( "Add global attributes and time variables" );
   if ( isfilevaratt(nc,0,"Conventions") )then
      nco@Conventions = nc@Conventions;
   end if
   if ( isfilevaratt(nc,0,"source") )then
      nco@source      = nc@source;
   end if
   nco@history  = ldate+": getndepdatFrom20thCentury.ncl";
   nco@source   = "Year "+sim_yr+" from: "+filename;
   nco@Version  = "$HeadURL$";
   nco@Revision = "$Id$";
   ;
   ; Add coordinate vars
   ;
   nco->lon = (/nc->lon/);
   nco->lat = (/nc->lat/);
   ;
   ; Find time of desired year
   ;
   tim = closest_val( stringtoint( sim_yr ), (/nc->YEAR/) );
   print( "Time index closest is :"+tim+" year = "+nc->YEAR(tim)+" time="+nc->time(tim) );
   ;
   ; Add other vars
   ;
   do i = 0, dimsizes(vars)-1
      if ( vars(i) .ne. "lat" .and. vars(i) .ne. "lon" )then
          nco->$vars(i)$(:,:) = nc->$vars(i)$(tim,:,:);
      end if
   end do

   print( "================================================================================================" );
   print( "Successfully created output ndepdat file: "+outfilename );

end
