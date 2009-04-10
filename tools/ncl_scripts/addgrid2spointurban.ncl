;
; Add grid data to the single-point urban input datasets for specific
; cities (MexicoCity, Vancouver, etc.)
;
;  Erik Kluzek
;  Mar/27/2009
;  $Id$
;  $HeadURL$
;
begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE
  ; Edit the following as needed to work on the next dataset
  ;
  ; Input and output resolution
  ;
  res      = getenv("RES");   ; Get output resolution from env variable

  ; Defaults if env variable is NOT set

  if ( ismissing(res) )then
     res      = "1x1_mexicocityMEX";  ; Output resolution
  end if

  ; END OF SECTION TO EDIT
  ; ===========================================================================================================
  ;
  ; Get input and output filenames
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
     querynml = clmroot+"/models/lnd/clm*/"+querynml;
  end if

  furbfile    = "../../bld/urban_input/" + systemfunc( querynml+" -res "+res  +" -var furbinp -filenameonly" );
  gridfile    = systemfunc( querynml+" -res "+res  +" -var fatmgrid" );
  outfile     = res+".nc";
  ;
  ; Open files
  ;
  print( "Add grid data to file: "+furbfile );
  if ( systemfunc("test -f "+furbfile+"; echo $?" ) .ne. 0 )then
     print( "Input furbinp file does not exist or not found: "+furbfile );
     exit
  end if
  nc  = addfile( furbfile,    "r" );

  print( "Grid data from file: "+gridfile );
  if ( systemfunc("test -f "+gridfile+"; echo $?" ) .ne. 0 )then
     print( "Input gridfile file does not exist or not found: "+gridfile );
     exit
  end if
  ncg = addfile( gridfile,    "r" );

  system( "/bin/rm "+outfile );
  print( "Output to file: "+outfile );
  nco  = addfile( outfile,    "c" );

  ;
  ; Define dimensions
  ;
  print( "Define dimensions" );
  dimnames = (/ "lsmlat", "lsmlon", "numsolar", "numrad", "nlevurb" /);
  dims     = (/ 1,         1,        2,          2,        15        /);
  is_unlim = (/  False,    False,    False,      False,    False /);
  filedimdef( nco, dimnames, dims, is_unlim );
  ; 
  ; Define vars and add attributes from grid file
  ;
  print( "Define variables" );
  gvars = getfilevarnames( ncg )
  i = dimsizes(gvars) - 1
  do while ( i .ge. 0 )
     print( "Define gridfile variable: "+gvars(i) );
     dimlist = getfilevardims( ncg, gvars(i) )
     filevardef (    nco, gvars(i), typeof(ncg->$gvars(i)$), dimlist )
     filevarattdef ( nco, gvars(i), ncg->$gvars(i)$ );
     delete( dimlist );
     i = i - 1;
  end do
  nco@Other_Info  = ncg@Other_Info;
  nco@Conventions = ncg@Conventions;
  ; 
  ; Define vars and add attributes from original file
  ;
  vars = getfilevarnames( nc )
  i = dimsizes(vars) - 1
  do while ( i .ge. 0 )
     print( "Define variable: "+vars(i) );
     dimlist = getfilevardims( nc, vars(i) )
     filevardef (    nco, vars(i), typeof(nc->$vars(i)$), dimlist )
     filevarattdef ( nco, vars(i), nc->$vars(i)$ );
     delete( dimlist );
     i = i - 1;
  end do
  fileattdef(  nco, nc );

  ;
  ; Add PCT_URBAN and LANDMASK
  ;
  addvars   = (/ "PCT_URBAN",     "LANDMASK",  "mask" /);
  addvarsln = (/ "percent urban", "land mask", "mask" /);
  do i = 0, dimsizes( addvars ) - 1
     filevardef ( nco, addvars(i), "float", (/ "lsmlat", "lsmlon" /) )
     nco->$addvars(i)$@long_name = addvarsln(i)
     nco->$addvars(i)$@units     = "unitless";
  end do

  sdate = systemfunc( "date" );
  nco@History_Log = "Added grid data from "+gridfile+" on "+sdate + "grid file "+ncg@History_Log;

  print( "Now copy the actual variable data on the file" );
  ;
  ; Copy actual variable data from grid file
  ;
  do i = 0, dimsizes( gvars ) - 1
     nco->$gvars(i)$ = (/ ncg->$gvars(i)$ /);
  end do
  ;
  ; Copy actual variable data from original file
  ;
  do i = 0, dimsizes( vars ) - 1
     if ( typeof(nc->$vars(i)$) .ne. "int" )then
        nco->$vars(i)$ = (/ nc->$vars(i)$ /);
     end if
  end do

  nco->PCT_URBAN = 100.0;
  nco->LANDMASK  =   1.0;
  nco->mask      =   1.0;
  
  delete( vars  );
  delete( gvars );

  print( "================================================================================================" );
  print( "Successfully created output urban file: "+outfile );

end
