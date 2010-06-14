;
; Convert a pftdyn test dataset to raw datasets that can be run through mksurfdata.
;
; Erik Kluzek
; April/17/2009
; $Id$
; $HeadURL$
;
begin

   load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl";

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
   print( "Use input pftdyn file: "+filename );
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
   print( "Use frac file: "+gridfile );
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
   print( "Use harvest pftdyn file: "+dynfile );
   if ( systemfunc("test -f "+dynfile+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+dynfile );
      exit
   end if
   ncd = addfile( dynfile, "r" );

   function ReOrder360LONToNeg180( longitudes [*]:numeric )
   ;
   ; Reorder longitudes from 0 to 360 to -180 to 180
   ; and return indices to swap other arrays as well
   ;
   local swapIndex;
   begin

      nlons     = dimsizes(longitudes);
      swapIndex = ispan( 0, nlons-1, 1 );
      swapsum   = sum( swapIndex );
      midpt     = minind( abs(longitudes-180.0d00) );
      j         = nlons - midpt - 1;
      do i = 0, dimsizes(longitudes)-1
         j = j + 1;
         if ( i .eq. midpt )then
            j = 0;
         end if
         swapIndex(i) = j;
      end do
      if ( sum( swapIndex ) .ne. swapsum )then
         print( "swapIndex is NOT correct: "+swapIndex );
         exit
      end if

      newlongs = (/ longitudes /);
      do i = 0, dimsizes(longitudes)-1
         newlongs(i) = (/ longitudes( swapIndex(i) ) /);
         if ( newlongs(i) .ge. 180.0 )then
            newlongs(i) = newlongs(i) - 360.0;
         end if
      end do

      longitudes = (/ newlongs /);

      return( swapIndex );

   end

   function ReOrder2D( swapIndex, array[*][*]:numeric )
   ; 
   ; Reorder longitude dimension of input 2D array by input swapIndex
   ;
   local newarray;
   begin

      if ( dimsizes(swapIndex) .ne. dimsizes(array(0,:)) )then
         print( "ERROR: swapIndex and 2nd dim do NOT match" );
         exit
      end if
      if ( array!1 .ne. "lsmlon" )then
         print( "ERROR: 2nd dim is NOT lsmlon as expected" );
         exit
      end if
      newarray = (/ array /);
      do i = 0, dimsizes(swapIndex)-1
         newarray(:,i) = (/ array( :, swapIndex(i) ) /);
      end do

      return( newarray );
   end

   function ReOrder3D( swapIndex, array[*][*][*]:numeric )
   ; 
   ; Reorder longitude dimension of input 3D array by input swapIndex
   ;
   local newarray;
   begin
      if ( dimsizes(swapIndex) .ne. dimsizes(array(0,0,:)) )then
         print( "ERROR: swapIndex and 3rd dim do NOT match" );
         exit
      end if
      if ( array!2 .ne. "lsmlon" )then
         print( "ERROR: 3rd dim is NOT lsmlon as expected" );
         exit
      end if
      newarray = (/ array /);
      do i = 0, dimsizes(swapIndex)-1
         newarray(:,:,i) = (/ array( :, :, swapIndex(i) ) /);
      end do

      return( newarray );
   end

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
   files= (/    "nc", "nc",    "nc",    "nc",    "ncg",      "ncg", "nc",      "ncg", "nc",    "nc",      "ncd",     "ncd",         "ncd",        "ncd",         "ncd",          "ncd" /);
   ftype= (/  "file1D","file1D","file1D","file1D","file", "var",   "file",    "var", "file",   "notime",  "notime",   "notime",      "notime",       "notime",     "notime",      "notime" /);
   LAT       = nc->LATIXY( :, 0 );
   LON       = nc->LONGXY( 0, : );
   swapIndex = ReOrder360LONToNeg180( (/ LON /) );
   print( "vars on file:"+vars );
   LONE = nc->LONE( 0, : );
   LONW = nc->LONW( 0, : );
   s = ReOrder360LONToNeg180( (/LONE/)  );
   s = ReOrder360LONToNeg180( (/LONW/)  );
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
            var = ncf->$vars(i)$;
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
               nco->$vars(i)$ = ReOrder3D( swapIndex, ncf->$vars(i)$(t,:,:,:) );
            else
               nco->$vars(i)$ = ReOrder2D( swapIndex, ncf->$vars(i)$(t,:,:) );
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
               if ( ftype(i) .eq. "file1D" )then
                  nco->$vars(i)$ = (/ ncf->$vars(i)$ /);
               else
                  nco->$vars(i)$= ReOrder2D( swapIndex, ncf->$vars(i)$ );
                  if ( vars(i) .eq. "LONGXY" )then
                     var = where ( nco->$vars(i)$ .ge. 180.0, nco->$vars(i)$ - 360.0, nco->$vars(i)$ )
                     nco->$vars(i)$ = var;
                     delete(var);
                  end if
               end if
            end if
         end if
      end do
      nco->EDGEE = -LONE(nlon-1);
      nco->EDGEW = LONW(0)-360.0;

   end do

   delete( vars  );
   delete( ftype );
   delete( files );

   print( "================================================================================================" );
   print( "Successfully created output "+filetype+" file: "+outfilename );

   ; ===========================================================================================================


end
