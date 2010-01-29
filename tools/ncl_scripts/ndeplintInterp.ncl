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
   rcp      = getenv("RCP");   ; Get output RCP (representative concentration pathway) from env variable

   if ( ismissing(res) )then
      res = "1.9x2.5";
   end if
   if ( ismissing(rcp) )then
      rcp = "8.5";
   end if

   nyrsb4           = 6;   ; Number of years to copy first year before to...
   lastHistoricYear = 2005 ; Last year of historical data (could be 2005 or 1995)
                           ; if 2005 will use previous historical data until 2005 (which doesn't have data from 2006-2009)
                           ; if 1995 will use historical data until 1995 (and scenario data for 2005 2000-2005=historical, 2006-2009=scenario)
   copyLastHist     = False; If should copy the last historical year and then use it to interpolate to the first future year
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
   filetype  = "fatmgrid";
   gridfile  = systemfunc( querynml+" -res "+res+" -var "+filetype );
   print( "Use "+filetype+" file: "+gridfile );
   if ( systemfunc("test -f "+gridfile+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+gridfile );
      exit
   end if
   ncg       = addfile( gridfile, "r" );

   ;
   ; Set the simulation years to get data for...
   ;
   sim_years = (/ 1855, 1865, 1875, 1885, 1895, 1905, 1915, 1925, 1935, 1945, 1955, 1965, 1975, 1985, 1995, 2005, 2015, 2025, 2035, 2045, 2055, 2065, 2075, 2085, 2095, 2105 /);
   if ( copyLastHist )then
      ;
      ; Add in an extra year that's a copy of the last year, but labeled as the next year
      ;
      years = new( (/ dimsizes(sim_years)+1 /), typeof(sim_years) );
      histInd = ind( sim_years .eq. lastHistoricYear );
      if ( ismissing(histInd) )then
         print( "Last Historic Year of "+lastHistoicYear+" does NOT exist in sim_years array" );
         exit
      end if
      years(0:histInd)  = sim_years(0:histInd);
      years(histInd+1)  = sim_years(histInd)+1;
      repeatYear        = years(histInd+1);
      years(histInd+2:) = sim_years(histInd+1:);
      delete( sim_years );
      sim_years = years;
      delete( years );
   else
      repeatYear = 0;
   end if
   ;
   ; Setup arrays that will be used..
   ;
   filetype  = "fndepdat";
   nfiles    = dimsizes(sim_years)
   nc        = new( (/ nfiles /), "file"   );
   filenames = new( (/ nfiles /), "string" );
   do yr = 0, nfiles-1
      year = sim_years(yr);
      if ( year .eq. repeatYear )then
         year = repeatYear - 1;
      end if
      if ( year .gt. lastHistoricYear )then
         queryopts = ",rcp="+rcp;
         printopts = " rcp="+rcp;
      else
         printopts = " ";
         queryopts = " ";
      end if
      filenames(yr) = systemfunc( querynml+" -res "+res+" -var "+filetype+" -options bgc=cn,sim_year="+year+queryopts);
      print( "Use "+filetype+" file: "+filenames(yr)+" for sim_year="+sim_years(yr)+printopts );
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
   print( "Finished opening all files" );

   beg_sim_year   = sim_years(0)-nyrsb4;
   end_sim_year   = sim_years(nfiles-1)+1;
   sim_year_range = beg_sim_year+"-"+end_sim_year;
   nyears         = end_sim_year - beg_sim_year + 1;

   ;
   ; Get area from grid file and convert to m^2
   ;
   print( "Check that units of area are as expected" );
   vname = "AREA";
   units = ncg->$vname$@units;
   print( "Units of area:"+units );
   if ( units .ne. "km^2" )then
      print( "Units of area is NOT km^2 as expected = "+units );
      exit
   end if
   area = dble2flt( ncg->$vname$ ) * 1000000.0 ;  Area in m^2 from fatmgrid file which has it in kg^2
   pi   = 3.14159265358979323846d00
   re   = 6371.22d03
   print( "Sum of area:"+sum(area)+" Expected: "+4.0d00*pi*re*re );
   ;
   ; Check that units and variables are what is expected
   ;
   vnames = (/ "NDEP_year",      "NDEP_AER_year"  /);
   lnames = (/ "NOy deposition", "NHx deposition" /);
   print( "Check that units and variables are as expected" );
   ncy0   = nc(0);
   do i = 0, dimsizes(vnames)-1
      if ( .not. isfilevar( ncy0, vnames(i) ) )then
         print( "Variable: "+vnames(i)+" is NOT on the file :"+filenames(0)+" as expected" );
         exit
      end if
      if ( ncy0->$vnames(i)$@units .ne. "kg(N)/s" )then
         print( "Units of variables are NOT kg(N)/s as expected for "+vnames(i)+" = "+ncg->$vnames(i)$@units );
         exit
      end if
      if ( ncy0->$vnames(i)$@long_name .ne. lnames(i) )then
         print( "Long_name of variable "+vnames(i)+" does NOT equal the expected value of "+lnames(i) );
         exit
      end if
      print( "Variable: "+vnames(i)+" exists with long_name="+lnames(i)+" and units "+ncy0->$vnames(i)$@units );
   end do
   ;
   ; Get date time-stamp to put on output file
   ;
   sdate     = systemfunc( "date +%y%m%d" );
   ldate     = systemfunc( "date" );

   outfilename = "fndep_clm_rcp"+rcp+"_simyr"+sim_year_range+"_"+res+"_c"+sdate+".nc";
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
   vars  = (/ "lon",  "lat",  "time", "YEAR", "NDEP_year",                     "NDEP_NOy_year",  "NDEP_NHx_year" /);
   lname = (/ "file", "file", "file", "file", "Sum of NOy and NHx deposition", "NOy deposition", "NHx deposition" /);
   ncy0        = nc(0);
   daysperyear = 365.0;
   Factor      = 1000.0*3600.0*24.0*daysperyear;    ; Unit conversion to convert from Kg-N/sec to g-N/yr
   do i = 0, dimsizes(vars)-1
      if ( vars(i) .ne. "YEAR" .and. vars(i) .ne. "time" )then
         print( "Add variable: "+vars(i) );
         if ( vars(i) .eq. "NDEP_year" .or. vars(i) .eq. "NDEP_NOy_year" .or. vars(i) .eq. "NDEP_NHx_year")then
            dimlist = dimnames;
         else
            dimlist = getfilevardims( nc(0), vars(i) )
         end if
         if ( lname(i) .eq. "file" )then
            filevardef ( nco, vars(i), typeof(ncy0->$vars(i)$), dimlist );
            filevarattdef ( nco, vars(i), ncy0->$vars(i)$ );
         else
            filevardef ( nco, vars(i), "float", dimlist );
            nco->$vars(i)$@long_name = lname(i);
            nco->$vars(i)$@units     = "g(N)/m2/yr";
         end if
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
   ;
   ; Add global attributes
   ;
   print( "Add global attributes and time variables" );
   if ( isfilevaratt(ncy0,0,"Conventions") )then
      nco@Conventions = ncy0@Conventions;
   end if
   if ( isfilevaratt(ncy0,0,"source") )then
      nco@source      = ncy0@source;
   end if
   do yr = 0, nfiles-1
      ncy           = nc(yr);
      if ( isfilevaratt(ncy,0,"history") )then
         history       = "history_"+sim_years(yr);
         nco@$history$ = ncy@history;
      end if
      if ( isfilevaratt(ncy,0,"case") )then
         case          = "case_"+sim_years(yr);
         nco@$case$    = ncy@case;
      end if
      source        = "source_"+sim_years(yr);
      nco@$source$  = "Input file:"+filenames(yr);
   end do
   nco@history  = ldate+": linearly interpolate in time between files by ndeplintInterp.ncl";
   first_year   = beg_sim_year + nyrsb4 - 1;
   nco@comment  = beg_sim_year+" through "+first_year+" repeat "+sim_years(0)+" and "+end_sim_year+" is a repeat of "+sim_years(nfiles-1)
   nco@Version  = "$HeadURL$";
   nco@Revision = "$Id$";
   ;
   ; Add coordinate vars
   ;
   nco->lon = (/ncy0->lon/);
   nco->lat = (/ncy0->lat/);
   ;
   ; Years before first and first year is just first file
   ;
   print( "Copy first year to first year and also years before first year" );
   do i = 0, nyrsb4
      nco->NDEP_NOy_year(i,:,:) = (/ ncy0->NDEP_year(:,:) /) * Factor / area;
      nco->NDEP_NHx_year(i,:,:) = (/ ncy0->NDEP_AER_year(:,:) /) * Factor / area;
      nco->YEAR(i)              = beg_sim_year + i;
      nco->time(i)              = int2dble(beg_sim_year)*daysperyear + daysperyear*i;
      print( "year="+nco->YEAR(i)+" NOy = "+avg(nco->NDEP_NOy_year(i,:,:)) );
   end do
   ;
   ; Year after last and last year is just from last file
   ;
   print( "Copy last year to last year and also year after last year as well as first year to first and year before first" );
   ncyn = nc(nfiles-1);
   do i = -2, -1;
      nco->NDEP_NOy_year(nyears+i,:,:) = (/ ncyn->NDEP_year(:,:) /) * Factor / area;
      nco->NDEP_NHx_year(nyears+i,:,:) = (/ ncyn->NDEP_AER_year(:,:) /) * Factor / area;
      nco->YEAR(nyears+i)              =  end_sim_year+i+1;
      nco->time(nyears+i)              =  int2dble(end_sim_year+i)*daysperyear
      print( "year="+nco->YEAR(nyears+i)+" NOy = "+avg(nco->NDEP_NOy_year(nyears+i,:,:)) );
   end do
   ;
   ; Loop through years in between now...
   ;
   yr   = 0;
   n    = nyrsb4 + 1;
   ncy0 = nc(0);
   ncy1 = nc(1);
   print( "low bound year="+sim_years(0)+"      avg(NDEP)="+avg(ncy0->NDEP_year(:,:)*Factor/area) );
   print( "hi  bound year="+sim_years(1)+"      avg(NDEP)="+avg(ncy1->NDEP_year(:,:)*Factor/area) );
   do year = sim_years(0)+1, sim_years(nfiles-1)-1
      if ( year .gt. sim_years(yr+1) ) then
         yr = yr + 1;
         ncy0 = nc(yr);
         ncy1 = nc(yr+1);
         print( "low bound year="+sim_years(yr)+  "      avg(NDEP)="+avg(ncy0->NDEP_year(:,:)*Factor/area) );
         print( "hi  bound year="+sim_years(yr+1)+"      avg(NDEP)="+avg(ncy1->NDEP_year(:,:)*Factor/area) );
      end if

      t1                        = int2flt(year-sim_years(yr))/int2flt(sim_years(yr+1)-sim_years(yr));
      t0                        = 1.0 - t1;
      nco->NDEP_NOy_year(n,:,:) = (/ t0*ncy0->NDEP_year(:,:)     + t1*ncy1->NDEP_year(:,:) /);
      nco->NDEP_NOy_year(n,:,:) = (/ nco->NDEP_NOy_year(n,:,:) /) * Factor / area;
      nco->NDEP_NHx_year(n,:,:) = (/ t0*ncy0->NDEP_AER_year(:,:) + t1*ncy1->NDEP_AER_year(:,:) /);
      nco->NDEP_NHx_year(n,:,:) = (/ nco->NDEP_NHx_year(n,:,:) /) * Factor / area;
      nco->YEAR(n)              = year;
      nco->time(n)              = int2flt(year)*daysperyear;
      print( "year = "+year+" t0 ="+t0+" t1 = "+t1+" avg(NDEP)="+avg(nco->NDEP_NOy_year(n,:,:)) );
      n  = n  + 1;

   end do
   print( "Sum NOy and NHx together:" );

   nco->NDEP_year(:,:,:) = (/ nco->NDEP_NOy_year(:,:,:) /) + (/ nco->NDEP_NHx_year(:,:,:) /);
   print( "Sum NOy+NHx "+beg_sim_year+" = "+avg(nco->NDEP_year(0,:,:)) )
   print( "Sum NOy+NHx "+end_sim_year+" = "+avg(nco->NDEP_year(nyears-1,:,:)) )

   print( "================================================================================================" );
   print( "Successfully created output ndepdyn file: "+outfilename );

end
