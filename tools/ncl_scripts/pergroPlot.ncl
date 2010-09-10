;
; NCL program to plot up root mean square differences from two different
; sources. This is as a way to verify that a port to a new machine is valid.
; It can also be used to verify that code changes are merely roundoff level.
;
;  Erik Kluzek
;  Apr/15/2010
;  $Id$
;  $HeadURL;
;
begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: Enter input using environment varibles
  ;
  ; RMS differences for input datasets
  ;
  nfiles     = 5;
  rmsfile    = new( (/ nfiles /), "string" );
  rmsfile(0) = getenv("RMSDAT");     ; Filename of first  ASCII file with RMS differences
  rmsfile(1) = getenv("RMSDAT2");    ; Filename of second ASCII file with RMS differences
  rmsfile(2) = getenv("RMSDAT3");    ; Filename of third  ASCII file with RMS differences
  rmsfile(3) = getenv("RMSDAT4");    ; Filename of fourth ASCII file with RMS differences
  rmsfile(4) = getenv("RMSDAT5");    ; Filename of fifth  ASCII file with RMS differences
  type       = getenv("TYPE");       ; Type of plot x11 to display to screen and ps for postscript
  ; ===========================================================================================================
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl";
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl";

  ; Defaults if env variables are NOT set
  if ( ismissing(rmsfile(0)) )then
     rmsfile(0) = "RMSintrepid.dat";
  end if
  if ( ismissing(rmsfile(1)) )then
     rmsfile(1) = "RMSbluefire.dat";
  end if
  if ( ismissing(type) )then
     type = "ps";
  end if
  ;
  ; Open the files and read the data
  ;
  print( "Check the input RMSDAT files and get the number of time samples in each" );
  ; First get the number of files and the # of lines in each file
  do i = 0, dimsizes(rmsfile)-1
     if ( ismissing(rmsfile(i)) )then
        nfiles = i
        break
     else
        if ( rmsfile(i) .eq. " " )then
           rmsfile(i) = rmsfile@_FillValue;
           nfiles = i
           break
        end if
        print( "Check: "+rmsfile(i) );
        if ( systemfunc("test -f "+rmsfile(i)+"; echo $?" ) .ne. 0 )then
           print( "Input RMS file does not exist or not found: "+rmsfile(i));
           exit
        end if
     end if
     ntimes2 = stringtointeger( systemfunc( "wc -l "+rmsfile(i) ) );
     if ( i .gt. 0 )then
        ntimes2 = min( (/ntimes,ntimes2/) )
     end if
     ntimes = ntimes2
  end do
  if ( ntimes .le. 0 .or. nfiles .le. 0 )then
     print( "One or more of the input RMS files is empty" );
     exit
  end if

  data = new( (/ nfiles, ntimes /), "float" );

  do i = 0, nfiles-1
     print( "Read in "+ntimes+" time-steps from file: "+rmsfile(i) );
     data(i,:)  = asciiread( rmsfile(i), ntimes, "float" );
     if ( all(ismissing(data(i,:)) ) )then
        print( "Could NOT read any data in from this file" );
        exit
     end if
  end do
  data = where( (data .eq. 0.0 ), data@_FillValue, data );

  data!0     = "files";
  data!1     = "time";
  data&files = rmsfile(0:nfiles-1);
  tics       = ispan( 0, ntimes-1, 1 );

  data@long_name = "RMS Difference";
  tics@long_name = "Time steps";
  tics@units     = "unitless";
  ;
  ; Do the plot
  ;
  plotfile = "pergro";
  if ( type .ne. "x11" )then
     print( "Plot out to file: "+plotfile+"."+type );
  end if
  wks   = gsn_open_wks ( type, plotfile );          ; open workstation

  lineThick = (/1.0,2.0,3.0,4.0,5.0/);
  lineColor = (/"red","green","navyblue","maroon","goldenrod"/);
  lineDash  = (/  1, 0, 2, 3, 4 /);

  res                   = True;                    ; plot mods desired
  res@tiMainString      = "Error Growth Plot";     ; add title
  res@xyComputeYMin     = True;                    ; Compute the Y min
  res@xyYStyle          = "Log";                   ; Do a log plot in the Y axis
  res@xyLineThicknesses = lineThick(0:nfiles-1)    ; line thickness
  res@xyLineColors      = lineColor(0:nfiles-1)    ; line colors
  res@xyDashPatterns    = lineDash(0:nfiles-1)     ; dash patterns

  print( "Do the plot without the legend" );
  plot  = gsn_csm_xy (wks,tics,data,res);          ; create plot

  ;
  ; Create legend
  ;
  lgres = True;
  lgres@lgLineColors        = res@xyLineColors;
  lgres@lgMonoLineThickness = False
  lgres@lgLineThicknesses   = res@xyLineThicknesses;
  lgres@lgMonoDashIndex     = False
  lgres@lgDashIndexes       = res@xyDashPatterns;
  lgres@lgItemType          = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF  = .08            ; legend label font thickness
  lgres@vpWidthF            = 0.13           ; width of legend (NDC)
  lgres@vpHeightF           = 0.10           ; height of legend (NDC)
  lgres@lgPerimThicknessF   = 2.0            ; thicken the box perimeter
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
  print( "Successfully created Perturbation Error Growth Plots" );

end
