;
;  mkscripgrid.ncl
;
;  Create SCRIP grid and mapping file for a land-only point or region.
;
;  Erik Kluzek
;  Dec/07/2011
;  $Id$
;  $HeadURL;
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
  ; ===========================================================================================================
  ; Set a few constants needed later
  cdate  = systemfunc( "date +%y%m%d" );
  ldate  = systemfunc( "date" );
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ; Edit the following as needed to interpolate to a new resolution.
  ;
  ; Input resolution and position
  ;
  name = getenv("PTNAME");                    ; Get name of this point

  latS = stringtodouble( getenv("S_LAT") );   ; Get south   latitude from env variable
  latN = stringtodouble( getenv("N_LAT") );   ; Get north   latitude from env variable
  lonE = stringtodouble( getenv("E_LON") );   ; Get east   longitude from env variable
  lonW = stringtodouble( getenv("W_LON") );   ; Get west   longitude from env variable

  nx   = stringtointeger( getenv("NX" )  );   ; Get number of grids along longitude lines
  ny   = stringtointeger( getenv("NY" )  );   ; Get number of grids along latitude lines

  imask = stringtointeger( getenv("IMASK") ); ; Get imask to use     from env variable

  print_str  = getenv("PRINT");               ; Do Extra printing for debugging

  outfilename = getenv("GRIDFILE");           ; Get filename from env variable

  if ( ismissing(nx) )then
     nx = 1;
  end if
  if ( ismissing(ny) )then
     ny = 1;
  end if
  if ( ismissing(imask) )then
     imask = 1;
  end if
  if ( ismissing(name) )then
     name = nx+"x"+ny+"pt_US-UMB";
  end if
  if ( ismissing(latS) )then
     latS = 45.5098;
  end if
  if ( ismissing(latN) )then
     latN = 45.6098;
  end if
  if ( ismissing(lonW) )then
     lonW = 275.2362;
  end if
  if ( ismissing(lonE) )then
     lonE = 275.3362;
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

  if ( ismissing(outfilename) )then
     if ( imask .eq. 1 )then
        outfilename = "SCRIPgrid_"+name+"_nomask_c"+cdate+".nc";
     else
        if ( imask .eq. 0 )then
           outfilename = "SCRIPgrid_"+name+"_noocean_c"+cdate+".nc";
        else
           outfilename = "SCRIPgrid_"+name+"_mask_c"+cdate+".nc";
        end if
     end if
  end if

  system( "/bin/rm -f "+outfilename );
  if ( printn )then
     print( "output file: "+outfilename );
  end if
  nc = addfile( outfilename, "c" );
  ;
  ; Define dimensions
  ;
  ns = nx * ny;
  dimnames = (/ "grid_size", "grid_corners", "grid_rank" /);
  dsizes   = (/ ns,           4,              2/);
  is_unlim = (/ False, False, False /);
  filedimdef( nc, dimnames, dsizes, is_unlim );

  ;
  ; Define variables
  ;
  dvars = (/ "grid_dims", "grid_center_lat", "grid_center_lon", "grid_corner_lat", "grid_corner_lon" /);

  filevardef ( nc, "grid_dims",       "integer", (/ "grid_rank" /) );
  filevardef ( nc, "grid_center_lat", "double",  (/ "grid_size" /) );
  filevardef ( nc, "grid_center_lon", "double",  (/ "grid_size" /) );
  filevardef ( nc, "grid_imask",      "integer", (/ "grid_size" /) );
  filevardef ( nc, "grid_corner_lat", "double",  (/ "grid_size", "grid_corners" /) );
  filevardef ( nc, "grid_corner_lon", "double",  (/ "grid_size", "grid_corners" /) );

  do i = 0, dimsizes(dvars)-1
    nc->$dvars(i)$@units = "degrees";
  end do
  nc->grid_imask@units = "unitless";
  ;
  ; Add global attributes to file
  ;
  nc@title   = "SCRIP grid file for "+name;
  nc@history = ldate+": create using mkscripgrid.ncl";
  nc@comment = "Ocean is assumed to non-existant at this point";
  nc@Version = "$HeadURL$";
  nc@Revision = "$Id$";
  ;
  ; Add data
  ;
  latxy      = new( (/ ns /), "double"  );
  lonxy      = new( (/ ns /), "double"  );
  latSxy     = new( (/ ns /), "double"  );
  lonExy     = new( (/ ns /), "double"  );
  latNxy     = new( (/ ns /), "double"  );
  lonWxy     = new( (/ ns /), "double"  );
  delX       = (lonE - lonW) / int2dble(nx);
  delY       = (latN - latS) / int2dble(ny);
  latSy      = fspan( latS,      latN-delY, ny );
  lonWx      = fspan( lonW,      lonE-delX, nx );
  latNy      = fspan( latS+delY, latN,      ny );
  lonEx      = fspan( lonW+delX, lonE,      nx );
  do i = 0, ns-1, nx
     j = i/nx;
     latSxy(i:i+nx-1) = latSy(j);
     latNxy(i:i+nx-1) = latNy(j);
     lonExy(i:i+nx-1) = lonEx;
     lonWxy(i:i+nx-1) = lonWx;
  end do
  ; Center latitude/longitudes
  latxy      = (latSxy + latNxy)/2.0d00
  lonxy      = (lonExy + lonWxy)/2.0d00
  ;
  nc->grid_dims       = (/ ns /);
  nc->grid_center_lat = (/ latxy /);
  nc->grid_center_lon = (/ lonxy /);
  nc->grid_imask      = (/ imask /);
  nc->grid_corner_lat(:,0) = (/ latSxy /);
  nc->grid_corner_lat(:,1) = (/ latSxy /);
  nc->grid_corner_lat(:,2) = (/ latNxy /);
  nc->grid_corner_lat(:,3) = (/ latNxy /);
  nc->grid_corner_lon(:,0) = (/ lonWxy /);
  nc->grid_corner_lon(:,1) = (/ lonExy /);
  nc->grid_corner_lon(:,2) = (/ lonExy /);
  nc->grid_corner_lon(:,3) = (/ lonWxy /);

  if ( printn )then
    print( "================================================================================================" );
    print( "Successfully created SCRIP grid file: "+outfilename);
  end if

  ; ===========================================================================================================

end
