;
; Soft link history files in a baseline directory to a filename with a date in the name
;
; NOTE: Requires at least NCL version 6.0.0 or later...
;
;  Erik Kluzek
;  Sep/14/2016
;

begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ;
  dir     = getenv("DIRECTORY");   ; Directory name to operate on
  newdir  = getenv("NEWDIR");      ; New directory name to have links for
  debug   = getenv("DEBUG");       ; DEBUG actually run softlink command or not
  verbose = getenv("VERBOSE");     ; VERBOSE printing or not...

  if ( ismissing(dir) )then
     print( "DIRECTORY is missing -- need to provide the directory to process" );
     dir = "/fs/cgd/csm/ccsm_baselines/clm4_5_12_r194/SMS_D_Ld1_Mmpi-serial.1x1_vancouverCAN.I1PTCLM50.hobart_nag.clm-default";
     ;status_exit( -1 );
  end if
  if ( ismissing(newdir) )then
     print( "NEWDIR is missing -- need to provide the new directory to create" );
     newdir = "/fs/cgd/csm/ccsm_baselines/clm4_5_12_r194.cime5/SMS_D_Ld1_Mmpi-serial.1x1_vancouverCAN.I1PTCLM50.hobart_nag.clm-default";
     ;status_exit( -1 );
  end if
  print( "Directory to operate over: "+dir );
  lvars  = (/ "DEBUG", "VERBOSE" /);
  varVal = (/ debug, verbose /);
  do i = 0, dimsizes(lvars)-1
     if ( ismissing(varVal(i)) )then
        ;print( lvars(i)+" not set" );
        varVal(i) = "FALSE";
     end if
     if ( varVal(i) .ne. "TRUE" .and. varVal(i) .ne. "FALSE" )then
        print( lvars(i)+" not set to TRUE or FALSE" );
        status_exit( -1 );
     end if
     ;print( lvars(i)+" is "+varVal(i) );
  end do
  debug   = varVal(0);
  verbose = varVal(1);
  hist_files = systemfunc("cd "+dir+"; ls *.h?.nc cism.h.nc");
  if ( verbose .eq. "TRUE" ) then
     print( hist_files );
  end if

  do i = 0, dimsizes(hist_files)-1
     nc = addfile( dir+"/"+hist_files(i), "r" );
     if ( verbose .eq. "TRUE" ) then
        print( "File:"+hist_files(i)+" Time:"+nc->time(0)+" "+nc->time@units );
     end if
     timearray = cd_calendar( nc->time, -5 );
     if ( verbose .eq. "TRUE" ) then
        print( "timearray = "+timearray );
     end if
     seconds = timearray(0,3)*3600 + timearray(0,4)*60 + timearray(0,5);
     if ( verbose .eq. "TRUE" ) then
        print( "seconds = "+seconds );
     end if
     ; Monthly files
     if ( dimsizes(nc->time) .eq. 1 .and. timearray(0,2) .eq. 1 )then
        date = sprinti( "%0.4d", timearray(0,0) )+sprinti("-%0.2d", timearray(0,1) )+sprinti("-%0.2d", timearray(0,2) );
     else
        date = sprinti( "%0.4d", timearray(0,0) )+sprinti("-%0.2d", timearray(0,1) )+sprinti("-%0.2d", timearray(0,2) )+sprinti("-%0.5d", seconds );
     end if
     if ( verbose .eq. "TRUE" ) then
        print( "date = "+date );
     end if
     newname = str_sub_str( hist_files(i), ".nc", "."+date+".nc" );
     if ( verbose .eq. "TRUE" ) then
        print( "newname= "+newname );
     end if
     cmd = "ln -s "+dir+"/"+hist_files(i)+" "+newdir+"/"+newname
     if ( verbose .eq. "TRUE" ) then
        print( "cmd= "+cmd );
     end if
     if ( debug .ne. "TRUE" ) then
        system( cmd );
     end if
     delete( timearray )
  end do
  print( "Successfully softlinked files to names with dates in the name" );
   
end
