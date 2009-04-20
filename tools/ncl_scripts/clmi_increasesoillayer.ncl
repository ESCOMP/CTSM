;
; clmi_increasesoillayer.ncl: script to increase the number of soil layers on input
; clmi initial conditions file from 10 to 15. Copies the lowest layer to the extra
; layers created.
;
;  Erik Kluzek
;  Mar/17/2008
;  $Id$
;  $HeadURL;
;
begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE
  ; Edit the following as needed to interpolate to a new resolution.
  ;
  ; Input resolution for input data
  ;
  resin    = getenv("RES");         ; Get resolution from env variable
  maskin   = getenv("MASK");        ; Get land mask description from env variable
  bgc      = getenv("BGC");         ; Get BGC mode from env variable

  ; Defaults if env variable is NOT set

  if ( ismissing(resin) )then
     resin    = "1.9x2.5";       ; Input resolution
  end if
  if ( ismissing(maskin) )then
     maskin   = "gx1v5";         ; Input ocean mask
  end if
  if ( ismissing(bgc) )then
     bgc      = "none";          ; Biogeochemistry mode ("none, "cn" etc.)
  end if
  levgrnd  = 15;              ; Number of ground soil levels for output file


  ; END OF SECTION TO EDIT
  ; ===========================================================================================================

  print ( "resolution = "+resin   );
  print ( "mask       = "+maskin  );
  print ( "levgrnd    = "+levgrnd );
  ;
  ; Use above to get filenames, and various other information needed such as dates
  ;
  csmdata  = getenv("CSMDATA");
  clmroot  = getenv("CLM_ROOT");
  querynml = "bld/queryDefaultNamelist.pl -silent -justvalue -res "+resin+" -options mask="+maskin+",BGC="+bgc+" -var finidat";
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
  filename    = systemfunc( querynml );
  ;
  ; Open input file
  ;
  print( "Copy from file: "+filename );
  if ( systemfunc("test -f "+filename+"; echo $?" ) .ne. 0 )then
     print( "Input clmi file does not exist or not found: "+filename );
     exit
  end if
  nc       = addfile( filename,    "r" );
  ;
  ; Open output file
  ;
  curymd = nc->timemgr_rst_curr_ymd;
  year   = curymd / 10000;
  month  = (curymd - year*10000) / 100;
  day    = curymd - year*10000 - month*100;
  yearS  = sprinti( "%4.4d", year  );
  monthS = sprinti( "%2.2d", month );
  dayS   = sprinti( "%2.2d", day   );
  resinS = "fv"+resin;

  filenameout = "clmi."+"."+yearS+"-"+monthS+"-"+dayS+"_"+resinS+"_"+maskin+"_c"+sdate+".nc";
  system( "/bin/rm -f "+filenameout );
  print( "Write out to file: "+filenameout );
  if ( resin .eq. "0.47x0.63" ) then
     print( "Use large file format!" );
     setfileoption( "nc", "Format", "64BitOffset" );
  end if
  ncfo         = addfile( filenameout, "c" );
  ;
  ; Define dimensions
  ;
  dim_names  = getvardims( nc );
  ;; Rename levsoi to levgrnd
  dlvsoi = 5;
  if ( dim_names(dlvsoi) .eq. "levsoi" )then
     dim_names(dlvsoi) = "levgrnd";
  else
     print( "Input dimension number "+dlvsoi+" is NOT levsoi as expected" );
     exit
  end if
  dim_sizes    = getfiledimsizes( nc );
  dimUnlim  = new( dimsizes(dim_sizes), logical, "No_FillValue" );
  dimUnlim(:) = False;
  levsoi    = dim_sizes(dlvsoi);
  if ( levgrnd .le. levsoi )then
     print( "Output levgrnd less than or equal to input levsoi : "+levsoi + " \n" );
     exit
  end if
  dlvtot         = 8;
  if ( dim_names(dlvtot) .ne. "levtot" )then
     print( "Input dimension number "+dlvtot+" is NOT levtot as expected" );
     exit
  end if
  levtotin       = dim_sizes(dlvtot);
  levtot         = levtotin + levgrnd - levsoi;
  dim_sizes(dlvsoi) = levgrnd;
  dim_sizes(dlvtot) = levtot;
  ndims = dimsizes(dim_names) - 1;
  print( "Define dimensions" );
  filedimdef( ncfo, (/dim_names(1:ndims)/), (/dim_sizes(1:ndims)/), (/dimUnlim(1:ndims)/) );
  delete( dim_names );
  delete( dim_sizes   );
  delete( dimUnlim );
  ;
  ; Define vars and add attributes from original file
  ;
  vars   = getfilevarnames( nc );
  isGrnd = new( dimsizes(vars), logical, "No_FillValue" );
  isGrnd(:) = False;
  isSoi     = new( dimsizes(vars), logical, "No_FillValue" );
  isSoi(:)  = False;
  i      = dimsizes(vars) - 1
  print( "Define variables" );
  do while ( i .ge. 0 )
     dimlist = getfilevardims( nc, vars(i) )
     do j = 0, dimsizes(dimlist)-1
       if ( dimlist(j) .eq. "levtot" )then
          isGrnd(i)  = True;
       end if
       if ( dimlist(j) .eq. "levsoi" )then
          isSoi(i)   = True;
          dimlist(j) = "levgrnd";
       end if
     end do
     filevardef (    ncfo, vars(i), typeof(nc->$vars(i)$), dimlist );
     filevarattdef ( ncfo, vars(i), nc->$vars(i)$ );
     delete( dimlist );
     i = i - 1
  end do
  fileattdef ( ncfo, nc );
  ;
  ; Add global attributes to output file describing what happened
  ;
  ncfo@history = ldate+": increase soil layers from "+levsoi+" to "+levgrnd+" by clmi_increasesoillayer.ncl";
  ncfo@source  = "Copy from:"+filename+" increasing soil layers from "+levsoi+" to "+levgrnd;
  ;
  do i = 0, dimsizes( vars )-1

     if ( isGrnd(i) .eq. True )then
        print( "Copy levtot variable: "+vars(i) );
        ncfo->$vars(i)$(:,0:levtotin-1) = (/ nc->$vars(i)$(:,:) /);
        do j = levtotin, levtot-1
           ncfo->$vars(i)$(:,j) = (/ nc->$vars(i)$(:,levtotin-1) /);
        end do
     else
        if ( isSoi(i) .eq. True )then
           print( "Copy levsoi variable: "+vars(i) );
           ncfo->$vars(i)$(:,0:levsoi-1) = (/ nc->$vars(i)$(:,:) /);
           do j = levsoi, levgrnd-1
              ncfo->$vars(i)$(:,j) = (/ nc->$vars(i)$(:,levsoi-1) /);
           end do
        else
           print( "Copy ordinary variable: "+vars(i) );
           ncfo->$vars(i)$ = (/ nc->$vars(i)$ /);
        end if
     end if

  end do

  delete( vars   );
  delete( isGrnd );
  delete( isSoi  );

  print( "================================================================================================" );
  print( "Successfully created output clmi file: "+filenameout );

end
