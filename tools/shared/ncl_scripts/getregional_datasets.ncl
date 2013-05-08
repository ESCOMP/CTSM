;
; Extract out regional datasets needed to run clm from the global datasets.
; NOTE: Requires at least NCL version 5.1.0 or later...
;
;  Erik Kluzek
;  Aug/28/2009
;  $Id$
;  $HeadURL;
;
begin
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl";
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ; Edit the following as needed to interpolate to a new resolution.
  ;
  ; Input resolution and position
  ;
  latS       = stringtodouble( getenv("S_LAT") );   ; Get south latitude from env variable
  latN       = stringtodouble( getenv("N_LAT") );   ; Get north latitude from env variable
  lonE       = stringtodouble( getenv("E_LON") );   ; Get east longitude from env variable
  lonW       = stringtodouble( getenv("W_LON") );   ; Get west longitude from env variable
  id         = getenv("CLM_USRDAT_NAME")            ; Get mydataid from env variable
  res        = getenv("RES");                       ; Get resolution from env variable
  rcp        = getenv("RCP");                       ; Get representative concentration pathway from env variable
  lmask      = getenv("MASK");                      ; Get input land/use mask to use from env variable
  sim_year   = stringtointeger( getenv("SIM_YR") ); ; Get input simulation year to use from env variable
  sim_yr_rng = getenv("SIM_YR_RNG");                ; Get input simulation year range to use
  debug_str  = getenv("DEBUG");                     ; Don't run just -- debug
  nomv_str   = getenv("NOMV");                      ; Don't move files -- leave in current directory
  print_str  = getenv("PRINT");                     ; Do Extra printing for debugging
  csmdata    = getenv("MYCSMDATA");                 ; Your personal CESM inputdata location

  if ( ismissing(res) )then
     res = "1.9x2.5";
  end if
  if ( ismissing(rcp) )then
     rcp = "-999.9";
  end if
  if ( ismissing(lmask) )then
     lmask = "gx1v6";
  end if
  if ( ismissing(latS) )then
     latS = 52.0d00;
  end if
  if ( ismissing(latN) )then
     latN = 73.0d00;
  end if
  if ( ismissing(lonW) )then
     lonW = 190.0d00;
  end if
  if ( ismissing(lonE) )then
     lonE = 220.0d00;
  end if
  if ( ismissing(sim_year) )then
     sim_year = 2000;
  end if
  if ( ismissing(sim_yr_rng) )then
     sim_yr_rng = "constant";
  end if
  if ( ismissing(print_str) )then
     printn = False;
  else
     if ( print_str .eq. "TRUE" )then
        printn = True;
     else
        printn = False;
     end if
  end if
  if ( ismissing(debug_str) )then
     debug = False;
  else
     if ( debug_str .eq. "TRUE" )then
        print( "DEBUG is TRUE do extra printing AND do NOT execute -- just print what WOULD happen" );
        debug  = True;
        printn = True;
     else
        debug = False;
     end if
  end if
  if ( ismissing(nomv_str) )then
     nomv = False;
  else
     if ( nomv_str .eq. "TRUE" )then
        print( "NOMV is TRUE do NOT move files leave them in current directory" );
        nomv = True;
     else
        nomv = False;
     end if
  end if
  if ( ismissing(id) )then
     id = "13x12pt_f19_alaskaUSA";
  end if
  if ( ismissing(csmdata) )then
     csmdata    = getenv("CSMDATA");   ; Standard CESM inputdata location
     if ( ismissing(csmdata) )then
        csmdata = "/fs/cgd/csm/inputdata";
     end if
  end if
  print( "Extract out regional datasets from global datasets" );
  if ( printn .eq. True )then
    print( "Global:   Resolution="+res+" mask="+lmask+" simyear="+sim_year );
    if ( sim_yr_rng .ne. "constant" )then
       print( "          sim_year_range="+sim_yr_rng+" rcp="+rcp );
    end if
    print( "Regional: id="+id+" Latitude="+latS+"-"+latN+" Longitude="+lonW+"-"+lonE );
  end if

  ;
  ; Setup the namelist query script
  ;
  ldate     = systemfunc( "date" );
  clmroot   = getenv("CLM_ROOT");
  querynml  = "bld/queryDefaultNamelist.pl -silent -justvalue ";
  querycfg  = "bld/config_query -silent -justvalue ";
  queryopts = " -options bgc=cn,mask="+lmask+",rcp="+rcp;
  querynyr  = ",sim_year="+sim_year+",sim_year_range="+sim_yr_rng;      ; Query for standard file years
  if ( sim_yr_rng .eq. "constant" )then
     queryyrs = ",sim_year=1850,sim_year_range=1850-2000";             ; Query for files that are always transient
  else
     queryyrs = querynyr;
  end if
  if ( .not. ismissing(csmdata) )then
     querynml = querynml+" -csmdata "+csmdata;
  end if
  if ( ismissing(clmroot) )then
     querynml = "../../"+querynml;
     querycfg = "../../"+querycfg;
  else
     querynml = clmroot+"/models/lnd/clm/"+querynml;
     querycfg = clmroot+"/models/lnd/clm/"+querycfg;
  end if
  GRID     = res+"_"+lmask;
  griddqry  = querycfg+GRID+" -id ATM_DOMAIN_PATH -justvalue";
  griddir  = systemfunc( griddqry );
  if ( ismissing(griddir) )then
     delete( griddqry )
     delete( griddir )
     gridqry  = querycfg+" ATM_DOMAIN_PATH -id value -justvalue -c config_definition.xml";
     griddir  = systemfunc( gridqry );
     if ( ismissing(griddir) )then
        print( "Input griddir does not exist or not found: query"+gridqry );
        exit
     end if
  end if
  gridqry  = querycfg+GRID+" -id ATM_DOMAIN_FILE -justvalue";
  gridfile = griddir + "/"+ systemfunc( gridqry );
  if ( printn .eq. True )then
     print( "          query = "+gridqry )
  end if
  ;
  ; Open file
  ;
  if ( systemfunc("test -f "+gridfile+"; echo $?" ) .ne. 0 )then
     print( "Input gridfile does not exist or not found: "+gridfile );
     exit
  end if
  if ( printn .eq. True )then
     print( "gridfile:"+gridfile );
  end if
  ncg     = addfile( gridfile,  "r" );

  if ( debug .eq. True )then
    print( "Env:      mycsmdata="+csmdata+" query="+querynml+queryopts+querynyr );
  end if
  indx = region_ind ( (/ncg->yc/), (/ncg->xc/), latS, latN, lonW, lonE );

  latdim = dimsizes(ncg->yc)
  londim = dimsizes(ncg->xc)
  if ( any( ismissing(indx)) )then
     print( "Indices:"+indx );
     print( "Missing indices found" );
     print( "LATIXY: "+ncg->yc );
     print( "LONGXY: "+ncg->xc );
     exit
  end if

  if ( debug .eq. True )then
     print( "Indices:"+indx );
  end if
  if ( printn .eq. True )then
     latdim = indx(3) - indx(2) + 1;
     londim = indx(1) - indx(0) + 1;
     print( "Grid size:"+latdim+"x"+londim );
     LOLAT = ncg->yc(indx(2),indx(0));
     HILAT = ncg->yc(indx(3),indx(1));
     print( "Actual grid span: Latitude="+LOLAT+"-"+HILAT );
     LOLON = ncg->xc(indx(2),indx(0));
     HILON = ncg->xc(indx(3),indx(1));
     print( "Actual grid span: Longitude="+LOLON+"-"+HILON );
  end if

  ;
  ; Setup filenames to process
  ;

  files     = (/ "fsurdat", "fpftdyn", "stream_fldfilename_ndep", "domainfile"      /);
  qryopt    = (/ querynyr,  querynyr,  queryyrs,                  "gridfile"        /);
  qnamelist = (/ " ",       " ",       "ndepdyn_nml",             " "               /);
  filelatnm = (/ "lsmlat",  "lsmlat",  "lat",                     "nj"              /);
  filelonnm = (/ "lsmlon",  "lsmlon",  "lon",                     "ni"              /);

  ;
  ; Loop over each of the files to process...
  ;
  do i = 0, dimsizes(files)-1
     ;
     ; Get the filename of the input global file and the output regional filename
     ;
     if ( qnamelist(i) .eq. " " ) then
        queryopt= " ";
     else
        queryopt= " -namelist "+qnamelist(i);
     end if
     if ( qryopt(i) .eq. "gridfile" )then
        globalfile = gridfile;
     else
        qry = querynml+queryopts+qryopt(i)+queryopt+" -res "+res+" -var "+files(i);
        globalfile = systemfunc( qry );
     end if
     if ( systemfunc("test -f "+globalfile+"; echo $?" ) .ne. 0 )then
        print( "Input global "+files(i)+" file does not exist or not found: "+globalfile );
        if ( printn .eq. True )then
           print( "qry was"+qry );
        end if
        continue;
     end if
     if ( debug .eq. True )then
        print( "Process file: "+globalfile );
     end if
     if ( qryopt(i) .eq. "gridfile" )then
        gridqry  = querycfg+"CLM_USRDAT -id ATM_DOMAIN_PATH -justvalue";
        griddir  = systemfunc( gridqry );
        qry  = querycfg+"CLM_USRDAT -id ATM_DOMAIN_FILE -justvalue";
        regfile = griddir + "/" + systemfunc( qry );
     else
        qry     = querynml+queryopts+qryopt(i)+queryopt+" -usrname "+id+" -var "+files(i);
        if ( nomv .eq. True )then
           qry  = qry+" -filenameonly";
        end if
        regfile = systemfunc( qry );
     end if
     if ( ismissing(regfile) )then
        print( "Output regional filename was NOT found: "+regfile );
        if ( printn .eq. True )then
           print( "qry was"+qry );
        end if
        continue;
     end if
     ;
     ; Run ncks on it over the region of interest
     ;
     cmd = "ncks -O -d "+filelatnm(i)+","+indx(0)+","+indx(1)+" -d "+filelonnm(i)+","+indx(2)+","+indx(3);
     cmd = cmd + " " + globalfile + " "+regfile;
     print( "Execute:"+cmd );
     if ( debug .eq. False )then
        if (  systemfunc( cmd+"; echo $?" ) .ne. 0 )then
           print( "Command did not complete successfully: " );
           exit
        end if
        ;
        ; Open up resultant file for writing
        ;
        nco = addfile( regfile, "w" );
        nco@history = nco@history + ":"+ldate + ": ";
        ;
        ; Add in coordinate variables for latitude/longitude
        ;
        var = filelatnm(i);
        if ( .not. isfilevar( nco, var ) .and. var .ne. "nj" )then
           if ( printn .eq. True )then
              print( "add "+var );
           end if
           nco@history = nco@history + " add "+var;
           nco->$var$  = (/ nco->LATIXY(:,0) /);
        end if
        var = filelonnm(i);
        if ( .not. isfilevar( nco, var ) .and. var .ne. "ni" )then
           if ( printn .eq. True )then
              print( "add "+var );
           end if
           nco@history = nco@history + " add "+var;
           nco->$var$  = (/ nco->LONGXY(0,:) /);
        end if
        ;
        ; Now check that internal variables are consistent
        ;
        var = "NUMLON";
        if( isfilevar( nco, var ) )then
           if ( printn .eq. True )then
              print( "modify "+var );
           end if
           nco@history = nco@history + " modify "+var;
           nco->$var$(:) = dimsizes( nco->LATIXY(0,:) );
        end if
        var = "EDGEN";
        if( isfilevar( nco, var ) )then
           if ( printn .eq. True )then
              print( "modify "+var );
           end if
           nco@history = nco@history + " modify "+var;
           nco->$var$ = max( nco->LATN );
        end if
        var = "EDGES";
        if( isfilevar( nco, var ) )then
           if ( printn .eq. True )then
              print( "modify "+var );
           end if
           nco@history = nco@history + " modify "+var;
           nco->$var$ = min( nco->LATS );
        end if
        var = "EDGEE";
        if( isfilevar( nco, var ) )then
           if ( printn .eq. True )then
              print( "modify "+var );
           end if
           nco@history = nco@history + " modify "+var;
           nco->$var$ = max( nco->LONE );
        end if
        var = "EDGEW";
        if( isfilevar( nco, var ) )then
           if ( printn .eq. True )then
              print( "modify "+var );
           end if
           nco@history = nco@history + " modify "+var;
           nco->$var$ = min( nco->LONW );
        end if
     end if
  end do

  print( "================================================================================================" );
  print( "Successfully created regional datasets from global datasets" );

end
