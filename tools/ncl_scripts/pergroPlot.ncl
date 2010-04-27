;
; NCL program to plot up root mean square differences from two different
; sources. This is as a way to verify that a port to a new machine is valid.
; It can also be used to verify that code changes are merely roundoff level.
;
;  Erik Kluzek
;  Apr/15/2010
;  $Id: ndepregrid.ncl 22307 2010-04-08 06:10:02Z erik $
;  $HeadURL;
;
begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: Enter input using environment varibles
  ;
  ; RMS differences for input datasets
  ;
  rmsfile1 = getenv("RMSDAT");     ; Filename of first ASCII file with RMS differences
  rmsfile2 = getenv("RMSDAT2");    ; Filename of second ASCII file with RMS differences
  type     = getenv("TYPE");       ; Type of plot x11 to display to screen and ps for postscript
  ; ===========================================================================================================
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl";
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl";

  ; Defaults if env variable is NOT set
  if ( ismissing(rmsfile1) )then
     rmsfile1 = "RMSintrepid.dat";
  end if
  if ( ismissing(rmsfile2) )then
     rmsfile2 = "RMSbluefire.dat";
  end if
  if ( ismissing(type) )then
     type = "ps";
  end if
  ;
  ; Open the files and read the data
  ;
  ; First get the number of lines in each file
  if ( systemfunc("test -f "+rmsfile1+"; echo $?" ) .ne. 0 )then
     print( "Input RMS file does not exist or not found: "+rmsfile1);
     exit
  end if
  if ( systemfunc("test -f "+rmsfile2+"; echo $?" ) .ne. 0 )then
     print( "Input RMS file does not exist or not found: "+rmsfile2);
     exit
  end if
  ntimes1    = stringtointeger( systemfunc( "wc -l "+rmsfile1 ) );
  ntimes2    = stringtointeger( systemfunc( "wc -l "+rmsfile2 ) );
  ntimes     = min( (/ ntimes1, ntimes2 /) );
  if ( ntimes .le. 0 )then
     print( "One or both of the input RMS files is empty" );
     exit
  end if

  data       = new( (/ 2, ntimes /), "float" );

  print( "Read in "+ntimes+" time-steps from file: "+rmsfile1 );
  data(0,:)  = asciiread( rmsfile1, ntimes, "float" );
  if ( all(ismissing(data(0,:)) ) )then
     print( "Could NOT read any data in from this file" );
     exit
  end if

  print( "Read in "+ntimes+" time-steps from file: "+rmsfile2 );
  data(1,:)  = asciiread( rmsfile2, ntimes, "float" );
  if ( all(ismissing(data(1,:)) ) )then
     print( "Could NOT read any data in from this file" );
     exit
  end if

  data!0     = "files";
  data!1     = "time";
  data&files = (/ rmsfile1, rmsfile2 /);
  tics       = ispan( 0, ntimes-1, 1 );

  data@long_name = "RMS Difference";
  tics@long_name = "Time steps";
  tics@units     = "unitless";
  ;
  ; Do the plot
  ;
  wks   = gsn_open_wks ( type, "pergro");          ; open workstation

  res                   = True;                    ; plot mods desired
  res@tiMainString      = "Error Growth Plot";     ; add title
  res@xyLineThicknesses = (/1.0,2.0/);             ; make 2nd line thicker
  res@xyLineColors      = (/"red","green"/);       ; change line colors
  res@xyDashPatterns    = (/  1,   0 /);           ; dash patterns
  res@xyComputeYMin     = True;                    ; Compute the Y min
  res@xyYStyle          = "Log";                   ; Do a log plot in the Y axis

  print( "Do the plot without the legend" );
  plot  = gsn_csm_xy (wks,tics,data,res);          ; create plot

  ;
  ; Create legend
  ;
  lgres = True;
  lgres@lgLineColors      = res@xyLineColors;
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.10           ; height of legend (NDC)
  lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgnd = gsn_create_legend(wks,dimsizes(data&files),data&files,lgres);

  ;
  ; Use gsn_add_annotation to attach this legend to our existing plot.
  ; This way, if we resize the plot, the legend will stay with the
  ; plot and be resized automatically.
  ;
  amres                  = True;
  amres@amJust           = "TopLeft";   ; Use top left corner of box
                                        ; for determining its location.
  amres@amParallelPosF   = 0.25;        ; Move legend to right
  amres@amOrthogonalPosF = 0.25;        ; Move legend down.

  annoid = gsn_add_annotation(plot,lgnd,amres); ; add legend to plot

  print( "Now do a second plot with a legend added to it" );
  draw( plot );
  frame( wks );

end
