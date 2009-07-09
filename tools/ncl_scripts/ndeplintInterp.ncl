;
; Linearly interpolate nitrogen deposition between 1850 and 2000 from
; decadal averages to create a ndepdyn file.
;
; Erik Kluzek
; April/30/2009
; $Id$
; $HeadURL;
;
begin
   ; ===========================================================================================================

   res      = getenv("RES");   ; Get output resolution from env variable

   if ( ismissing(res) )then
      res = "1.9x2.5";
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
   filetype  = "fndepdat";
   sim_years = (/ 1850, 1855, 1865, 1875, 1885, 1895, 1905, 1915, 1925, 1935, 1945, 1955, 1965, 1975, 1985, 1995, 2000, 2005 /);
   nfiles    = dimsizes(sim_years)
   nc        = new( (/ nfiles /), "file"   );
   filenames = new( (/ nfiles /), "string" );
   do yr = 0, nfiles-1
      filenames(yr) = systemfunc( querynml+" -res "+res+" -var "+filetype+" -options bgc=cn,sim_year="+sim_years(yr) );
      print( "Use "+filetype+" file: "+filenames(yr)+" for sim_year="+sim_years(yr) );
      if ( systemfunc("test -f "+filenames(yr)+"; echo $?" ) .ne. 0 )then
         print( "Input "+filetype+" file does not exist or not found: "+filenames(yr) );
         exit
      end if
      nc(yr) = addfile( filenames(yr), "r" );
      ncy    = nc(yr);
      if ( yr .gt. 0 )then
         if ( dimsizes( ncy->lon ) .ne. nlon )then
             print( "Longitude is different size than previous file" );
             exit
         end if
         if ( dimsizes( ncy->lat ) .ne. nlat )then
             print( "Latitude is different size than previous file" );
             exit
         end if
      end if
      nlon = dimsizes( ncy->lon );
      nlat = dimsizes( ncy->lat );
   end do

   beg_sim_year   = sim_years(0)-1;
   end_sim_year   = sim_years(nfiles-1)+1;
   sim_year_range = beg_sim_year+"-"+end_sim_year;
   nyears         = end_sim_year - beg_sim_year + 1;
   ;
   ; Get date time-stamp to put on output file
   ;
   sdate     = systemfunc( "date +%y%m%d" );
   ldate     = systemfunc( "date" );

   outfilename = "fndep_clm_"+sim_year_range+"_"+res+"_c"+sdate+".nc";
   system( "/bin/rm -f "+outfilename );
   print( "output file: "+outfilename );
   nco = addfile( outfilename, "c" );
   ;
   ; Define dimensions
   ;
   dimnames = (/ "time", "lat", "lon" /);
   dsizes   = (/ nyears, nlat,  nlon /);
   is_unlim = (/ True, False, False /);
   filedimdef( nco, dimnames, dsizes, is_unlim );
   ;
   ; Define variables
   ;
   vars = (/ "lon", "lat", "time", "YEAR", "NDEP_year", "NDEP_AER_year" /);
   ncy0 = nc(0);
   do i = 0, dimsizes(vars)-1
      if ( vars(i) .ne. "YEAR" .and. vars(i) .ne. "time" )then
         print( "Add variable: "+vars(i) );
         if ( vars(i) .eq. "NDEP_year" .or. vars(i) .eq. "NDEP_AER_year")then
            dimlist = dimnames;
         else
            dimlist = getfilevardims( nc(0), vars(i) )
         end if
         filevardef (    nco, vars(i), typeof(ncy0->$vars(i)$), dimlist );
         filevarattdef ( nco, vars(i), ncy0->$vars(i)$ );
         delete( dimlist );
      end if
   end do
   filevardef ( nco, "YEAR", "integer", (/ "time" /) );
   nco->YEAR@long_name = "year";
   nco->YEAR@units     = "Year AD";
   filevardef ( nco, "time", "double", (/ "time" /) );
   nco->time@long_name = "time";
   nco->time@calendar  = "noleap";
   nco->time@units     = "days since 0000-01-01 00:00";
   daysperyear         = 365.0;
   ;
   ; Add global attributes
   ;
   print( "Add global attributes and time variables" );
   nco@Conventions = ncy0@Conventions;
   nco@source      = ncy0@source;
   do yr = 0, nfiles-1
      ncy           = nc(yr);
      history       = "history_"+sim_years(yr);
      nco@$history$ = ncy@history;
      case          = "case_"+sim_years(yr);
      nco@$case$    = ncy@case;
      source        = "source_"+sim_years(yr);
      nco@$source$  = "Input file:"+filenames(yr);
   end do
   nco@history  = ldate+": linearly interpolate in time between files by ndeplintInterp.ncl";
   nco@comment  = beg_sim_year+" is a repeat of "+sim_years(0)+" and "+end_sim_year+" is a repeat of "+sim_years(nfiles-1)
   nco@Version  = "$HeadURL$";
   nco@Revision = "$Id$";
   ;
   ; Add coordinate vars
   ;
   nco->lon = (/ncy0->lon/);
   nco->lat = (/ncy0->lat/);
   ;
   ; Year before first and first year is just first file
   ;
   print( "Copy first year to first year and also year before first year" );
   nco->NDEP_year(0,:,:)     = (/ ncy0->NDEP_year(:,:) /);
   nco->NDEP_AER_year(0,:,:) = (/ ncy0->NDEP_AER_year(:,:) /);
   nco->NDEP_year(1,:,:)     = (/ ncy0->NDEP_year(:,:) /);
   nco->NDEP_AER_year(1,:,:) = (/ ncy0->NDEP_AER_year(:,:) /);
   nco->YEAR(0:1)            = (/ beg_sim_year, beg_sim_year+1 /);
   nco->time(0:1)            = (/ int2dble(beg_sim_year)*daysperyear, int2dble(beg_sim_year+1)*daysperyear /);
   ;
   ; Year after last and last year is just from last file
   ;
   print( "Copy last year to last year and also year after last year as well as first year to first and year before first" );
   ncyn = nc(nfiles-1);
   nco->NDEP_year(nyears-2,:,:)     = (/ ncyn->NDEP_year(:,:) /);
   nco->NDEP_AER_year(nyears-2,:,:) = (/ ncyn->NDEP_AER_year(:,:) /);
   nco->NDEP_year(nyears-1,:,:)     = (/ ncyn->NDEP_year(:,:) /);
   nco->NDEP_AER_year(nyears-1,:,:) = (/ ncyn->NDEP_AER_year(:,:) /);
   nco->YEAR(nyears-2:)             = (/ end_sim_year-1, end_sim_year /);
   nco->time(nyears-2:)             = (/ int2dble(end_sim_year-1)*daysperyear, int2dble(end_sim_year)*daysperyear /);
   ;
   ; Loop through years in between now...
   ;
   yr   = 0;
   n    = 2;
   ncy0 = nc(0);
   ncy1 = nc(1);
   print( "low bound year="+sim_years(0)+"      avg(NDEP)="+avg(ncy0->NDEP_year(:,:)) );
   print( "hi  bound year="+sim_years(1)+"      avg(NDEP)="+avg(ncy1->NDEP_year(:,:)) );
   do year = sim_years(0)+1, sim_years(nfiles-1)-1
      if ( year .gt. sim_years(yr+1) ) then
         yr = yr + 1;
         ncy0 = nc(yr);
         ncy1 = nc(yr+1);
         print( "low bound year="+sim_years(yr)+"      avg(NDEP)="+avg(ncy0->NDEP_year(:,:)) );
         print( "hi  bound year="+sim_years(yr+1)+"      avg(NDEP)="+avg(ncy1->NDEP_year(:,:)) );
      end if

      t1                        = int2flt(year-sim_years(yr))/int2flt(sim_years(yr+1)-sim_years(yr));
      t0                        = 1.0 - t1;
      nco->NDEP_year(n,:,:)     = (/ t0*ncy0->NDEP_year(:,:)     + t1*ncy1->NDEP_year(:,:) /);
      nco->NDEP_AER_year(n,:,:) = (/ t0*ncy0->NDEP_AER_year(:,:) + t1*ncy1->NDEP_AER_year(:,:) /);
      nco->YEAR(n)              = year;
      nco->time(n)              = int2flt(year)*daysperyear;
      print( "year = "+year+" t0 ="+t0+" t1 = "+t1+" avg(NDEP)="+avg(nco->NDEP_year(n,:,:)) );
      n  = n  + 1;

   end do

   print( "================================================================================================" );
   print( "Successfully created output ndepdyn file: "+outfilename );

end
