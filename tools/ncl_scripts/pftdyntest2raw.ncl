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
   filename = systemfunc( querynml+" -res "+res+" -var "+filetype );
   print( "Use file: "+filename );
   if ( systemfunc("test -f "+filename+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+filename );
      exit
   end if
   nc = addfile( filename, "r" );
   ;
   ; Get date time-stamp to put on output filenames
   ;
   sdate     = systemfunc( "date +%y%m%d" );
   ldate     = systemfunc( "date" );
   ;
   ; Get dimension info.
   ;
   varname   = "PCT_PFT";
   timdimnames  = getfilevardims(     nc,  varname );
   timdsizes    = getfilevardimsizes( nc,  varname );
   dimnames     = timdimnames( 1: );
   dimnames(0)  = "pft";
   dsizes       = timdsizes(   1: );
   is_unlim     = (/ False, False, False /);
   print( "dimensions:"+dimnames );

   ;
   ; Get variable info.
   ;
   vars = (/ "EDGEE", "EDGEN", "EDGES", "EDGEW", "LANDMASK", "LAT", "LATIXY", "LON", "LONGXY", "PCT_PFT" /);
   LANDFRAC = (/nc->LANDFRAC_PFT/);
   LANDMASK = new( dsizes(1:), "integer"  );
   LANDMASK = where( LANDFRAC > 0.0, 1, 0 );
   LANDMASK@long_name = "land mask";
   LANDMASK@units     = "unitless";
   LAT      = nc->LATIXY( :, 0 );
   LON      = nc->LONGXY( 0, : );
   print( "vars on file:"+vars );
   ;
   ; Now loop over each year and output files for each year
   ;
   ntimes = timdsizes( 0 );
   do t = 0, ntimes - 1
      year = nc->YEAR(t);
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
         if ( vars(i) .eq. "PCT_PFT" ) then
            dimlist = dimnames;
            vtype   = typeof(nc->$vars(i)$);
            var     = nc->$vars(i)$;
         else
            if ( vars(i) .eq. "LANDMASK" ) then
               dimlist = dimnames( 1: );
               vtype   = typeof(LANDMASK);
               var     = LANDMASK;
            else
               if ( vars(i) .eq. "LAT" ) then
                  dimlist = (/ "lsmlat" /);
                  vtype   = typeof(LAT);
                  var     = LAT;
               else
                  if ( vars(i) .eq. "LON" ) then
                     dimlist = (/ "lsmlon" /);
                     vtype   = typeof(LON);
                     var     = LON;
                  else
                     dimlist = getfilevardims( nc, vars(i) )
                     vtype   = typeof(nc->$vars(i)$);
                     var     = nc->$vars(i)$;
                  end if
               end if
            end if
         end if
         filevardef (    nco, vars(i), vtype, dimlist );
         filevarattdef ( nco, vars(i), var );
         delete( dimlist );
         delete( vtype   );
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
         if ( vars(i) .eq. "PCT_PFT" ) then
            nco->$vars(i)$ = (/ nc->$vars(i)$(t,:,:,:) /);
         else
            if ( vars(i) .eq. "LANDMASK" ) then
               nco->$vars(i)$ = (/ LANDMASK /);
            else
               if ( vars(i) .eq. "LAT" ) then
                  nco->$vars(i)$ = (/ LAT /);
               else
                  if ( vars(i) .eq. "LON" ) then
                     nco->$vars(i)$ = (/ LON /);
                  else
                     nco->$vars(i)$ = (/ nc->$vars(i)$ /);
                  end if
               end if
            end if
         end if
      end do

   end do

   delete( vars );

   print( "================================================================================================" );
   print( "Successfully created output "+filetype+" file: "+outfilename );

   ; ===========================================================================================================
end
