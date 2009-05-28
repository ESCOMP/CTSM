;
; Regrid CLM Nitrogen-Deposition file (ndep) to a new resolution. 
; Works for both climatology monthly ndep files as well as yearly dynamic files.
; Edit the top few lines of the file for a different input simulation year, output
; resolution, or file type (climatology of yearly dynamic).
; Uses bld query to get input filenames. Also uses env variable CSMDATA for location
; of input files.
;
;  Erik Kluzek
;  Dec/14/2007
;  $Id$
;  $HeadURL;
;
begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE
  ; Edit the following as needed to interpolate to a new resolution.
  ;
  ; Input and output resolution and sim_year for input data
  ;
  resin    = "1.9x2.5";          ; Input resolution:  1.9x2.5 (normally NOT changed)
  res      = getenv("RES");      ; Get output resolution from env variable
  sim_year = getenv("SIM_YR");   ; Get input simulation year from env variable

  ; END OF SECTION TO EDIT
  ; ===========================================================================================================

  ; Defaults if env variable is NOT set

  if ( ismissing(res) )then
     res      = "48x96";      ; Output Resolution: 1.9x2.5, 48x96 etc.
  end if
  if ( ismissing(sim_year) )then
     sim_year = "2000";       ; Input/Output simulation year: 1850, 2000
  end if

  ; Figure out file type from sim-year

  if ( sim_year .eq. "1850-2000" )then
     filetype = "fndepdyn"
  else
     filetype = "fndepdat"
  end if

  ; END OF SECTION TO EDIT
  ; ===========================================================================================================
  ;
  ; Use above to get filenames, and various other information needed such as dates
  ;
  if ( filetype .eq. "fndepdyn" )then
     ftype    = "fndep_clm_";
  else
     ftype    = "ndep_clm_";
  end if
  csmdata  = getenv("CSMDATA");
  clmroot  = getenv("CLM_ROOT");
  querynml = "bld/queryDefaultNamelist.pl -silent -justvalue -options sim_year="+sim_year;
  querynml = querynml+",bgc=cn"
  if ( .not. ismissing(csmdata) )then
     querynml = querynml+" -csmdata "+csmdata;
  end if
  if ( ismissing(clmroot) )then
     querynml = "../../"+querynml;
  else
     querynml = clmroot+"/models/lnd/clm*/"+querynml;
  end if

  sdate     = systemfunc( "date +%y%m%d" );
  ldate     = systemfunc( "date" );
  print( "query string="+querynml )
  filename    = systemfunc( querynml+" -res "+resin+" -var "+filetype );
  gridinfile  = systemfunc( querynml+" -res "+resin+" -var fatmgrid" );
  gridfile    = systemfunc( querynml+" -res "+res  +" -var fatmgrid" );
  filenameout = ftype+sim_year+"_"+res+"_c"+sdate+".nc";
  ;
  ; Open files
  ;
  print( "Interpolate from file: "+filename );
  if ( systemfunc("test -f "+filename+"; echo $?" ) .ne. 0 )then
     print( "Input ndep file does not exist or not found: "+filename );
     exit
  end if
  nc       = addfile( filename,    "r" );
  print( "Use gridin file:       "+gridinfile );
  if ( systemfunc("test -f "+gridinfile+"; echo $?" ) .ne. 0 )then
     print( "Input gridinfile does not exist or not found: "+gridinfile );
     exit
  end if
  ncgi     = addfile( gridinfile,  "r" );
  print( "Use grid file:         "+gridfile );
  if ( systemfunc("test -f "+gridfile+"; echo $?" ) .ne. 0 )then
     print( "Input gridfile does not exist or not found: "+gridfile );
     exit
  end if
  ncg      = addfile( gridfile,    "r" );
  print( "Output file:           "+filenameout );
  if ( systemfunc("test -f "+filenameout+"; echo $?" ) .eq. 0 )then
     system( "/bin/rm -f "+filenameout );
  end if
  nco      = addfile( filenameout, "c" );
  ;
  ; Define dimensions
  ;
  if ( isfilevardim( nc, "NDEP_year", "time" ) )then
     dimnames  = (/ "lon", "lat", "time" /);
     dims      = getfiledimsizes( nc );
     ntime     = dims(2);
     ncgdimsiz = getfiledimsizes( ncg )
     dsizes    = (/ ncgdimsiz(1), ncgdimsiz(2), ntime /);
     is_unlim  = (/ False, False, True /);
     filedimdef( nco, dimnames, dsizes, is_unlim );
  else
     dimnames  = (/ "lon", "lat" /);
     ncgdimsiz = getfiledimsizes( ncg )
     dsizes    = (/ ncgdimsiz(1), ncgdimsiz(2) /);
     is_unlim  = (/ False, False /);
     filedimdef( nco, dimnames, dsizes, is_unlim );
  end if
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
  ; Add global attributes to output file
  ;
  nco@history            = ldate+": Regrid from "+resin+" resolution to "+res+" by ndepregrid.ncl";
  nco@source             = "Regrid from:"+filename;
  nco@gridfile           = gridfile;
  nco@ndepregridVersion  = "$HeadURL$";
  nco@ndepregridRevision = "$Id$";
  ;
  ; Copy coordinate variables over
  ;
  ; For climatology files...
  if ( filetype .eq. "fndepdat" )then
     if ( isfilevar( nc, "month" ) ) then
        nco->month = (/nc->month/);
     end if
  ; for yearly dynamic files...
  else
     nco->YEAR = (/nc->YEAR/);
     nco->time = (/nc->time/);
  end if
  nco->lon= doubletofloat( (/ncg->LONGXY(lsmlat|0,lsmlon|:) /) );
  nco->lat= doubletofloat( (/ncg->LATIXY(lsmlat|:,lsmlon|0) /) );
  lon     = (/ nc->lon /);
  lat     = (/ nc->lat /);
  lono    = (/ nco->lon /);
  lato    = (/ nco->lat /);
  areai   = doubletofloat( (/ ncgi->AREA /) );
  sumarea = sum(areai);
  areai   = areai / sumarea;
  areao   = doubletofloat( (/ ncg->AREA /) );
  sumarea = sum(areao);
  areao   = areao / sumarea;
  print( "areai sum = "+sum(areai) );
  print( "areao sum = "+sum(areao) );
  if ( abs(sum(areai)-sum(areao)) .gt. 1.e-4 ) then
     print( "Input and output areas are significantly different" );
     print( "Area-in:  "+sum(areai) );
     print( "Area-out: "+sum(areao) );
     exit
  end if
  ;
  ; loop over variables
  ;
  Cyclic = True;
  do i = 0, dimsizes( vars )-1

     vardimnames = getfilevardims( nc, vars(i) );
     ;
     ; If variable is not one of the dimensions -- regrid it and write to output file
     ;
     if ( (vars(i) .ne. "YEAR") .and. (vars(i) .ne. "month") .and. (vars(i) .ne. "lat") .and. (vars(i) .ne. "lon") .and. (vars(i) .ne. "LAT") .and. (vars(i) .ne. "LON") )then
       print( "Write variable: "+vars(i)+" to output file" );
       ;
       ; If time dimension
       ;
       if ( vardimnames(0) .eq. "time" )then
          do t = 0, ntime-1
             vart  = (/nc->$vars(i)$(time|t,lat|:,lon|:)/);
             varto = linint2 ( lon,  lat,  vart, Cyclic, lono, lato, 0 );
             nco->$vars(i)$(time|t,lat|:,lon|:) = (/varto/);
             vart  = vart*areai;
             varto = varto*areao;
             print( "Month: "+t+" Input w-sum: "+sum(vart)+" Output: "+sum(varto) );
          end do
       ;
       ; without time dimension
       ;
       else
          var  = (/nc->$vars(i)$/);
          varo = linint2 ( lon, lat,  var, Cyclic, lono, lato, 0 );
          nco->$vars(i)$ = (/varo/);
          var  = var*areai;
          varo = varo*areao;
          print( "Input w-sum: "+sum(var)+" Output: "+sum(varo) );
       end if
     end if
     delete( vardimnames );

  end do
  if ( isvar("varto") )then
     delete(varto);
  end if
  if ( isvar("vart") )then
     delete(vart);
  end if
  delete(areao);
  delete(areai);
  if ( isvar("varo") )then
     delete(varo);
  end if
  if ( isvar("var") )then
     delete(var)
  end if

  print( "================================================================================================" );
  print( "Successfully created output ndep file: "+filenameout );

end
