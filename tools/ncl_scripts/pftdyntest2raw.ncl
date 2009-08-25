;
; Convert a pftdyn test dataset to raw datasets that can be run through mksurfdata.
;
; Erik Kluzek
; April/17/2009
; $Id$
; $HeadURL$
;
begin


   ; ===========================================================================================================

   res      = getenv("RES");   ; Get output resolution from env variable

   if ( ismissing(res) )then
      res = "1x1_tropicAtl";  ; resolution (10x15 or 1x1_tropicAtl)
   end if
   ; ===========================================================================================================
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
   ; Use resolution to get input filename and open it
   ;
   filetype = "fpftdyn";
   if ( res .eq. "1x1_tropicAtl" )then
     sim_yrs  = "1000-1004";
   else
     sim_yrs  = "1000-1002";
   end if
   filename = systemfunc( querynml+" -res "+res+" -options sim_year_range="+sim_yrs+" -var "+filetype );
   print( "Use file: "+filename );
   if ( systemfunc("test -f "+filename+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+filename );
      exit
   end if
   nc = addfile( filename, "r" );
   ;
   ; Use resolution to get input grid filename and open it
   ;
   filetype = "fatmlndfrc";
   gridfile = systemfunc( querynml+" -res "+res+" -var "+filetype );
   print( "Use file: "+gridfile );
   if ( systemfunc("test -f "+gridfile+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+gridfile );
      exit
   end if
   ncg = addfile( gridfile, "r" );
   ;
   ; Use resolution to get input dynamic PFT filename and open it for harvesting values
   ;
   filetype = "fpftdyn";
   dsim_yrs  = "1850-2000";
   dynfile = systemfunc( querynml+" -res "+res+" -options sim_year_range="+dsim_yrs+" -var "+filetype );
   print( "Use file: "+dynfile );
   if ( systemfunc("test -f "+dynfile+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+dynfile );
      exit
   end if
   ncd = addfile( dynfile, "r" );

   function getmyfile( filedesc:string )
   ;
   ; Return file handle for file described
   ;
   local myfile;
   begin
      if ( filedesc .eq. "nc" )then
         myfile = nc;      
      else
      if ( filedesc .eq. "ncg" )then
         myfile = ncg;      
      else
      if ( filedesc .eq. "ncd" )then
         myfile = ncd;      
      else
         print( "ERROR: bad file descriptor = "+filedesc );
         exit
      end if
      end if
      end if
      return( myfile );
   end

   ;
   ; Get date time-stamp to put on output filenames
   ;
   sdate     = systemfunc( "date +%y%m%d" );
   ldate     = systemfunc( "date" );
   ;
   ; Get dimension info.
   ;
   varname      = "PCT_PFT";
   dimnames     = (/ "pft", "lsmlat", "lsmlon" /);
   latgg        = ncg->LATIXY(lsmlat|:,lsmlon|0);
   longg        = ncg->LONGXY(lsmlat|0,lsmlon|:);
   nlat         = dimsizes( latgg );
   nlon         = dimsizes( longg );
   pft          = dimsizes( nc->$varname$(lsmpft|:,lsmlat|0,lsmlon|0,time|0) );
   numpft       = pft(0);
   dsizes       = (/ numpft, nlat, nlon /);
   is_unlim     = (/ False, False, False /);
   print( "dimensions:"+dimnames );

   ;
   ; Get variable info.
   ;
   vars = (/ "EDGEE", "EDGEN", "EDGES", "EDGEW", "LANDMASK", "LAT", "LATIXY", "LON", "LONGXY", "PCT_PFT", "GRAZING", "HARVEST_VH1", "HARVEST_VH2", "HARVEST_SH1", "HARVEST_SH2", "HARVEST_SH3" /);
   files= (/    "nc",    "nc",    "nc",    "nc",      "ncg", "ncg",     "nc", "ncg",     "nc",      "nc",     "ncd",         "ncd",         "ncd",        "ncd",         "ncd",         "ncd" /);
   ftype= (/  "file",  "file",  "file",  "file",     "file", "var",   "file", "var",  "file",   "notime",  "notime",      "notime",      "notime",       "notime",     "notime",      "notime" /);
   LAT      = nc->LATIXY( :, 0 );
   LON      = nc->LONGXY( 0, : );
   print( "vars on file:"+vars );
   ;
   ; Now loop over each year and output files for each year
   ;
   ntimes = dimsizes( nc->YEAR );
   do t = 0, ntimes - 1
      year = nc->YEAR(t);
      if ( year .lt. 10 ) then
         year        = 1000 + year;
      end if
      ; Open file for this year
      outfilename = "mksrf_pft_"+res+"_testyr"+year+"_c"+sdate+".nc";
      system( "/bin/rm -f "+outfilename );
      print( "output file: "+outfilename );
      nco = addfile( outfilename, "c" );
      ; Define dimensions
      filedimdef( nco, dimnames, dsizes, is_unlim );
      ;
      ; Define variables
      ;
      do i = 0, dimsizes(vars)-1
         ncf = getmyfile( files(i) );
         if ( ftype(i) .eq. "notime" ) then
            if ( vars(i) .eq. "PCT_PFT" ) then
               dimlist = dimnames;
            else
               dimlist = dimnames(1:);
            end if
            var     = ncf->$vars(i)$;
         else
            if ( ftype(i) .eq. "var" )then
                if ( vars(i) .eq. "LAT" ) then
                   dimlist = (/ "lsmlat" /);
                   var     = LAT;
                else
                if ( vars(i) .eq. "LON" ) then
                   dimlist = (/ "lsmlon" /);
                   var     = LON;
                end if
                end if
            else
               dimlist = getfilevardims( ncf, vars(i) )
               var     = ncf->$vars(i)$;
            end if
         end if
         filevardef (    nco, vars(i), typeof(var), dimlist );
         filevarattdef ( nco, vars(i), var );
         delete( dimlist );
         delete( var     );
      end do
      ;
      ; Add some attributes
      ;
      nco@Conventions   = nc@Conventions;
      nco@Logname       = nc@Logname;
      nco@creation_date = ldate;
      nco@history       = ldate + ": pftdyntest2raw.ncl res="+res;
      nco@version       = "$HeadURL$";
      nco@revision_id   = "$Id$";
      ;
      ; Now add the variables on the file
      ;
      do i = 0, dimsizes(vars)-1
         ncf = getmyfile( files(i) );
         if ( ftype(i) .eq. "notime" ) then
            if ( vars(i) .eq. "PCT_PFT" ) then
               nco->$vars(i)$ = (/ ncf->$vars(i)$(t,:,:,:) /);
            else
               nco->$vars(i)$ = (/ ncf->$vars(i)$(t,:,:) /);
            end if
         else
            if ( ftype(i) .eq. "var" )then
               if ( vars(i) .eq. "LAT" ) then
                  nco->$vars(i)$ = (/ LAT /);
               else
               if ( vars(i) .eq. "LON" ) then
                  nco->$vars(i)$ = (/ LON /);
               end if
               end if
            else
               nco->$vars(i)$ = (/ ncf->$vars(i)$ /);
            end if
         end if
      end do

   end do

   delete( vars  );
   delete( ftype );
   delete( files );

   print( "================================================================================================" );
   print( "Successfully created output "+filetype+" file: "+outfilename );

   ; ===========================================================================================================


end
