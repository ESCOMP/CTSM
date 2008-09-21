;**************************************
; NCL script
;
; generate_ascii_avg_fv1_9x2_5_urbanparam_file_p7.ncl
;
; This program takes the 0.5 resolution netcdf file of averaged parameters
; (created from generate_netcdf_avg_halfdeg_urbanparam_file_p7.ncl),
; aggregates it to the resolution of the output grid and outputs to an ascii file.
;
; Keith Oleson, July 2008
;  $Id$
;  $HeadURL;
;
;**************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;========================================================================================
begin
 ;
 ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE
 ; Edit the following as needed
 ;

 ;
 ; General operation flags
 ;
 do_generate_ascii = 0   ; ------ if should create an output ascii file or not
 res     = "10x15"       ; ------ output resolution to run at

 ; END OF SECTION TO EDIT
 ;=====================================================================================
 ;
 ;
 ; Use above to get filenames, and various other information needed such as dates
 ;

 csmdata  = getenv("CSMDATA");
 clmroot  = getenv("CLM_ROOT");
 querynml = "bld/queryDefaultNamelist.pl -silent -justvalue ";
 if ( .not. ismissing(csmdata) )then
    querynml = querynml+" -csmdata "+csmdata;
 else
    csmdata = "/fs/cgd/csm/inputdata";
 end if
 if ( ismissing(clmroot) )then
    querynml = "../../"+querynml;
 else
    querynml = clmroot+"/models/lnd/clm*/"+querynml;
 end if

 sdate     = systemfunc( "date +%y%m%d" );
 ldate     = systemfunc( "date" );
 print( "query string="+querynml );
 ;
 ; Input data files
 ;
 furbinp = systemfunc( querynml+" -res "+res+" -filenameonly -var furbinp" );
 srf_file= systemfunc( querynml+" -res "+res+" -options furbinp="+furbinp+" -var fsurdat" );

 cdf_file = csmdata+"/lnd/clm2/rawdata/urban_data_Feddema:080410_0.5x0.5_avg_c080725.nc" ; netCDF file to read.
 ;
 ; Set output data filenames
 ;
 cdf_filename = "urbdata."+res+"_avg_c"+sdate+".nc";
 outfile      = "avg_densityclass_"+res+"_fluxes_p7_c"+sdate+".txt";

 ;=======================================================================================

 print ("=========================================")
 print ("Start Time: "+systemfunc("date") )
 print ("=========================================")


 ;
 ; Parameter sizes
 ;
 nlevurb  = 10    ; Number of thickness layers for urban
 numrad   = 2     ; Number or radiation bands
 numsolar = 2     ; Number of solar bands

; Get surface data
 
 print( "Open surface dataset at output resolution "+srf_file );
 surfdata_file  = addfile( srf_file, "r")
 latn           = dble2flt(surfdata_file->LATN)
 lats           = dble2flt(surfdata_file->LATS)
 lat_lr         = dble2flt(surfdata_file->LATIXY(:,0))
 lon_lr         = dble2flt(surfdata_file->LONGXY(0,:))
 nlat_lr        = dimsizes(lat_lr)
 nlon_lr        = dimsizes(lon_lr)
 pct_urban_lr   = surfdata_file->PCT_URBAN
 pct_urban_lr1d = ndtooned(pct_urban_lr)
 edgee          = surfdata_file->EDGEE
 edgen          = surfdata_file->EDGEN
 edges          = surfdata_file->EDGES
 edgew          = surfdata_file->EDGEW
 latixy         = surfdata_file->LATIXY
 longxy         = surfdata_file->LONGXY
 print((/"Done reading in "+res+" surface data"/))

 pct_urban_lr!0   = "lat"
 pct_urban_lr&lat = lat_lr
 pct_urban_lr!1   = "lon"
 pct_urban_lr&lon = lon_lr
 pct_urban_lr     = lonFlip(pct_urban_lr)
 lon_lr_flipped   = (/pct_urban_lr&lon/)
 print((/"Done assigning "+res+" lat/lons"/))

; Set up the edges so they correspond to -180 to 180 grid
; Special case where lonw is max and lone is lonmin is handled later
 lonw      = new((/nlat_lr,nlon_lr/),"float")
 lone      = new((/nlat_lr,nlon_lr/),"float")

 dellon    = new((/1/),"float")
 dellat    = new((/1/),"float")
 dellon    = (/ dble2flt( (/longxy(0,1)/) - (/longxy(0,0)/) ) /)
 dellat    = (/ dble2flt( (/latixy(1,0)/) - (/latixy(0,0)/) ) /)
 lonmax    =  180.0 - dellon + dellon/2.
 lonmin    = -180.0 + dellon - dellon/2.
 lonw(:,0) = lonmax
 lonw(:,1) = lonmin
 do i = 2,nlon_lr-1
   lonw(:,i) = lonw(:,i-1) + dellon
 end do
 lone(:,0) = lonmin
 do i = 1,nlon_lr-1
   lone(:,i) = lone(:,i-1) + dellon
 end do
 print((/"Done setting up "+res+" lat/lon edges"/))

; Get data
 print( "Open high resolution urban raw data file: "+cdf_file );
 data_file          = addfile(cdf_file,"r")
 nlat_hr            = dimsizes(data_file->LAT)
 nlon_hr            = dimsizes(data_file->LON)
 canyon_hwr_hr      = data_file->CANYON_HWR
 wtlunit_roof_hr    = data_file->WTLUNIT_ROOF
 wtroad_perv_hr     = data_file->WTROAD_PERV
 em_roof_hr         = data_file->EM_ROOF
 em_wall_hr         = data_file->EM_WALL
 em_improad_hr      = data_file->EM_IMPROAD
 em_perroad_hr      = data_file->EM_PERROAD
 alb_roof_hr        = data_file->ALB_ROOF
 alb_wall_hr        = data_file->ALB_WALL
 alb_improad_hr     = data_file->ALB_IMPROAD
 alb_perroad_hr     = data_file->ALB_PERROAD
 ht_roof_hr         = data_file->HT_ROOF
 wind_hgt_canyon_hr = data_file->WIND_HGT_CANYON
 tk_roof_hr         = data_file->TK_ROOF
 tk_wall_hr         = data_file->TK_WALL
 tk_improad_hr      = data_file->TK_IMPROAD
 cv_roof_hr         = data_file->CV_ROOF
 cv_wall_hr         = data_file->CV_WALL
 cv_improad_hr      = data_file->CV_IMPROAD
 thick_roof_hr      = data_file->THICK_ROOF
 thick_wall_hr      = data_file->THICK_WALL
 nlevimproad_hr     = data_file->NLEV_IMPROAD
 t_building_min_hr  = data_file->T_BUILDING_MIN
 t_building_max_hr  = data_file->T_BUILDING_MAX
 print((/"Done getting data"/))

; Arrays for low resolution data
 canyon_hwr_lr      = new((/nlat_lr,nlon_lr/),"float")
 wtlunit_roof_lr    = new((/nlat_lr,nlon_lr/),"float")
 wtroad_perv_lr     = new((/nlat_lr,nlon_lr/),"float")
 em_roof_lr         = new((/nlat_lr,nlon_lr/),"float")
 em_wall_lr         = new((/nlat_lr,nlon_lr/),"float")
 em_improad_lr      = new((/nlat_lr,nlon_lr/),"float")
 em_perroad_lr      = new((/nlat_lr,nlon_lr/),"float")
 ht_roof_lr         = new((/nlat_lr,nlon_lr/),"float")
 wind_hgt_canyon_lr = new((/nlat_lr,nlon_lr/),"float")
 thick_roof_lr      = new((/nlat_lr,nlon_lr/),"float")
 thick_wall_lr      = new((/nlat_lr,nlon_lr/),"float")
 nlevimproad_lr     = new((/nlat_lr,nlon_lr/),"float")
 t_building_min_lr  = new((/nlat_lr,nlon_lr/),"float")
 t_building_max_lr  = new((/nlat_lr,nlon_lr/),"float")
 alb_roof_lr        = new((/numsolar,numrad,nlat_lr,nlon_lr/),"float")
 alb_wall_lr        = new((/numsolar,numrad,nlat_lr,nlon_lr/),"float")
 alb_improad_lr     = new((/numsolar,numrad,nlat_lr,nlon_lr/),"float")
 alb_perroad_lr     = new((/numsolar,numrad,nlat_lr,nlon_lr/),"float")
 tk_roof_lr         = new((/nlevurb,nlat_lr,nlon_lr/),"float")
 tk_wall_lr         = new((/nlevurb,nlat_lr,nlon_lr/),"float")
 tk_improad_lr      = new((/nlevurb,nlat_lr,nlon_lr/),"float")
 cv_roof_lr         = new((/nlevurb,nlat_lr,nlon_lr/),"float")
 cv_wall_lr         = new((/nlevurb,nlat_lr,nlon_lr/),"float")
 cv_improad_lr      = new((/nlevurb,nlat_lr,nlon_lr/),"float")

; Assign coordinates for low resolution data
 canyon_hwr_lr!0   = "lat"
 canyon_hwr_lr&lat = lat_lr
 canyon_hwr_lr!1   = "lon"
 canyon_hwr_lr&lon = lon_lr_flipped
 copy_VarCoords(canyon_hwr_lr,wtlunit_roof_lr)
 copy_VarCoords(canyon_hwr_lr,wtroad_perv_lr)
 copy_VarCoords(canyon_hwr_lr,em_roof_lr)
 copy_VarCoords(canyon_hwr_lr,em_improad_lr)
 copy_VarCoords(canyon_hwr_lr,em_perroad_lr)
 copy_VarCoords(canyon_hwr_lr,em_wall_lr)
 copy_VarCoords(canyon_hwr_lr,ht_roof_lr)
 copy_VarCoords(canyon_hwr_lr,wind_hgt_canyon_lr)
 copy_VarCoords(canyon_hwr_lr,thick_wall_lr)
 copy_VarCoords(canyon_hwr_lr,thick_roof_lr)
 copy_VarCoords(canyon_hwr_lr,nlevimproad_lr)
 copy_VarCoords(canyon_hwr_lr,t_building_min_lr)
 copy_VarCoords(canyon_hwr_lr,t_building_max_lr)

 alb_roof_lr!0        = "numsolar"
 alb_roof_lr&numsolar = ispan(1,numsolar,1)
 alb_roof_lr!1        = "numrad"
 alb_roof_lr&numrad   = ispan(1,numrad,1)
 alb_roof_lr!2        = "lat"
 alb_roof_lr&lat      = lat_lr
 alb_roof_lr!3        = "lon"
 alb_roof_lr&lon      = lon_lr_flipped
 copy_VarCoords(alb_roof_lr,alb_wall_lr)
 copy_VarCoords(alb_roof_lr,alb_improad_lr)
 copy_VarCoords(alb_roof_lr,alb_perroad_lr)
 tk_roof_lr!0       = "nlevurb"
 tk_roof_lr&nlevurb = ispan(1,nlevurb,1)
 tk_roof_lr!1       = "lat"
 tk_roof_lr&lat     = lat_lr
 tk_roof_lr!2       = "lon"
 tk_roof_lr&lon     = lon_lr_flipped
 copy_VarCoords(tk_roof_lr,tk_wall_lr)
 copy_VarCoords(tk_roof_lr,tk_improad_lr)
 copy_VarCoords(tk_roof_lr,cv_roof_lr)
 copy_VarCoords(tk_roof_lr,cv_wall_lr)
 copy_VarCoords(tk_roof_lr,cv_improad_lr)

 print((/"Begin calculating average parameters for each low resolution gridcell"/))
 ;
 ; Expand search area by slightly more than half one high resolution gridcell in each direction
 ;
 delta = 0.26
 do j = 0,nlat_lr-1
    print((/"lat: "+lat_lr(j)/))
    do i = 0,nlon_lr-1
       if (pct_urban_lr(j,i) .gt. 0.) then
          ; Special case where lonw is max and lone is min
          if (i .eq. 0) then
             nl   = 3
             nb   = nlon_hr - nl -1
             ne   = nlon_hr - 1
             nb2  = nl+1
             ne2  = nl*2+1
             tmp1 = (/canyon_hwr_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
             tmp2 = (/canyon_hwr_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
             numer = new((/dimsizes(tmp1(:,0)),nl*2+2/),"float")
             numer(:,0:nl) = tmp1
             numer(:,nb2:ne2) = tmp2
             if ( .not. all(ismissing(numer)) ) then
                canyon_hwr_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/wtlunit_roof_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/wtlunit_roof_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                wtlunit_roof_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/wtroad_perv_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/wtroad_perv_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                wtroad_perv_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/em_roof_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/em_roof_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                em_roof_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/em_improad_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/em_improad_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                em_improad_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/em_perroad_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/em_perroad_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                em_perroad_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/em_wall_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/em_wall_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                em_wall_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/ht_roof_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/ht_roof_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                ht_roof_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/wind_hgt_canyon_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/wind_hgt_canyon_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                wind_hgt_canyon_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/thick_wall_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/thick_wall_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                thick_wall_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/thick_roof_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/thick_roof_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                thick_roof_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/t_building_min_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/t_building_min_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                t_building_min_lr(j,i) = avg(numer)
                numer(:,0:nl) = (/t_building_max_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/t_building_max_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                t_building_max_lr(j,i) = avg(numer)
                delete(tmp1)
                delete(tmp2)
                delete(numer)

                tmp1 = (/alb_roof_hr(:,:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                tmp2 = (/alb_roof_hr(:,:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                numer = new((/dimsizes(tmp1(:,0,0,0)),dimsizes(tmp1(0,:,0,0)),dimsizes(tmp1(0,0,:,0)),nl*2+2/),"float")
                numer(:,:,:,0:nl) = tmp1
                numer(:,:,:,nb2:ne2) = tmp2
                do ss = 0,numsolar-1
                do ll = 0,numrad-1
                   alb_roof_lr(ss,ll,j,i) = avg(numer(ss,ll,:,:))
                end do
                end do
                numer(:,:,:,0:nl) = (/alb_improad_hr(:,:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,:,:,nb2:ne2) = (/alb_improad_hr(:,:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                do ss = 0,numsolar-1
                do ll = 0,numrad-1
                  alb_improad_lr(ss,ll,j,i) = avg(numer(ss,ll,:,:))
                end do
                end do
                numer(:,:,:,0:nl) = (/alb_perroad_hr(:,:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,:,:,nb2:ne2) = (/alb_perroad_hr(:,:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                do ss = 0,numsolar-1
                do ll = 0,numrad-1
                  alb_perroad_lr(ss,ll,j,i) = avg(numer(ss,ll,:,:))
                end do
                end do
                numer(:,:,:,0:nl) = (/alb_wall_hr(:,:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,:,:,nb2:ne2) = (/alb_wall_hr(:,:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                do ss = 0,numsolar-1
                do ll = 0,numrad-1
                  alb_wall_lr(ss,ll,j,i) = avg(numer(ss,ll,:,:))
                end do
                end do
                delete(tmp1)
                delete(tmp2)
                delete(numer)

                tmp1 = (/tk_roof_hr(:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                tmp2 = (/tk_roof_hr(:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                numer = new((/dimsizes(tmp1(:,0,0)),dimsizes(tmp1(0,:,0)),nl*2+2/),"float")
                numer(:,:,0:nl) = tmp1
                numer(:,:,nb2:ne2) = tmp2
                do ll = 0,nlevurb-1
                  tk_roof_lr(ll,j,i) = avg(numer(ll,:,:))
                end do
                numer(:,:,0:nl) = (/tk_wall_hr(:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,:,nb2:ne2) = (/tk_wall_hr(:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                do ll = 0,nlevurb-1
                  tk_wall_lr(ll,j,i) = avg(numer(ll,:,:))
                end do
                numer(:,:,0:nl) = (/tk_improad_hr(:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,:,nb2:ne2) = (/tk_improad_hr(:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                do ll = 0,nlevurb-1
                   if ( .not. all(ismissing(numer(ll,:,:))) ) then
                      tk_improad_lr(ll,j,i) = avg(numer(ll,:,:))
                   end if
                end do
                numer(:,:,0:nl) = (/cv_roof_hr(:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,:,nb2:ne2) = (/cv_roof_hr(:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                do ll = 0,nlevurb-1
                  cv_roof_lr(ll,j,i) = avg(numer(ll,:,:))
                end do
                numer(:,:,0:nl) = (/cv_wall_hr(:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,:,nb2:ne2) = (/cv_wall_hr(:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                do ll = 0,nlevurb-1
                  cv_wall_lr(ll,j,i) = avg(numer(ll,:,:))
                end do
                numer(:,:,0:nl) = (/cv_improad_hr(:,{lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,:,nb2:ne2) = (/cv_improad_hr(:,{lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                do ll = 0,nlevurb-1
                   if ( .not. all(ismissing(numer(ll,:,:))) ) then
                      cv_improad_lr(ll,j,i) = avg(numer(ll,:,:))
                   end if
                end do
                delete(tmp1)
                delete(tmp2)
                delete(numer)

                tmp1 = (/nlevimproad_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                tmp2 = (/nlevimproad_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                numer = new((/dimsizes(tmp1(:,0)),nl*2+2/),"float")
                numer(:,0:nl) = (/nlevimproad_hr({lats(j,i)-delta:latn(j,i)+delta},0:nl)/)
                numer(:,nb2:ne2) = (/nlevimproad_hr({lats(j,i)-delta:latn(j,i)+delta},nb:ne)/)
                if (.not.(all(ismissing(numer)))) then
                  nlevimproad_lr(j,i) = 2
                else 
                  nlevimproad_lr(j,i) = 0
                end if
             end if
             delete(tmp1)
             delete(tmp2)
             delete(numer)
          else
             if (.not.(all(ismissing(nlevimproad_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta})))) ) then
                canyon_hwr_lr(j,i)      = avg(canyon_hwr_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                wtlunit_roof_lr(j,i)    = avg(wtlunit_roof_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                wtroad_perv_lr(j,i)     = avg(wtroad_perv_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                em_roof_lr(j,i)         = avg(em_roof_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                em_improad_lr(j,i)      = avg(em_improad_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                em_perroad_lr(j,i)      = avg(em_perroad_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                em_wall_lr(j,i)         = avg(em_wall_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                ht_roof_lr(j,i)         = avg(ht_roof_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                wind_hgt_canyon_lr(j,i) = avg(wind_hgt_canyon_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                thick_wall_lr(j,i)      = avg(thick_wall_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                thick_roof_lr(j,i)      = avg(thick_roof_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                t_building_min_lr(j,i)  = avg(t_building_min_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                t_building_max_lr(j,i)  = avg(t_building_max_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
    
                do ss = 0,numsolar-1
                do ll = 0,numrad-1
                  alb_roof_lr(ss,ll,j,i)    = avg(alb_roof_hr(ss,ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  alb_improad_lr(ss,ll,j,i) = avg(alb_improad_hr(ss,ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  alb_perroad_lr(ss,ll,j,i) = avg(alb_perroad_hr(ss,ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  alb_wall_lr(ss,ll,j,i)    = avg(alb_wall_hr(ss,ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                end do
                end do
    
                do ll = 0,nlevurb-1
                  tk_wall_lr(ll,j,i)    = avg(tk_wall_hr(ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  tk_roof_lr(ll,j,i)    = avg(tk_roof_hr(ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  if( .not. all(ismissing((tk_improad_hr(ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta})))) )then
                      tk_improad_lr(ll,j,i) = avg(tk_improad_hr(ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  end if
                  cv_wall_lr(ll,j,i)    = avg(cv_wall_hr(ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  cv_roof_lr(ll,j,i)    = avg(cv_roof_hr(ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  if( .not. all(ismissing((cv_improad_hr(ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta})))) )then
                     cv_improad_lr(ll,j,i) = avg(cv_improad_hr(ll,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}))
                  end if
                end do
                if ( .not.(all(nlevimproad_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta}) .eq. 0.))) then
                   nlevimproad_lr(j,i) = 2
    ;              print((/nlevimproad_hr({lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta})/))
    ;              print((/tk_improad_lr(0:floattointeger(nlevimproad_lr(j,i))-1,j,i)/))
    ;              print((/tk_improad_hr(0:floattointeger(nlevimproad_lr(j,i))-1,{lats(j,i)-delta:latn(j,i)+delta},{lonw(j,i)-delta:lone(j,i)+delta})/))
                else
                   nlevimproad_lr(j,i) = 0
                end if
             else
                nlevimproad_lr(j,i) = 0
             end if

          end if
       end if
    end do
 end do
 print((/"Finish calculating average parameters for each low resolution gridcell"/))

 ; Must lon flip back to original grid
 canyon_hwr_lr      = lonFlip(canyon_hwr_lr)
 wtlunit_roof_lr    = lonFlip(wtlunit_roof_lr)
 wtroad_perv_lr     = lonFlip(wtroad_perv_lr)
 em_roof_lr         = lonFlip(em_roof_lr)
 em_improad_lr      = lonFlip(em_improad_lr)
 em_perroad_lr      = lonFlip(em_perroad_lr)
 em_wall_lr         = lonFlip(em_wall_lr)
 alb_roof_lr        = lonFlip(alb_roof_lr)
 alb_improad_lr     = lonFlip(alb_improad_lr)
 alb_perroad_lr     = lonFlip(alb_perroad_lr)
 alb_wall_lr        = lonFlip(alb_wall_lr)
 ht_roof_lr         = lonFlip(ht_roof_lr)
 wind_hgt_canyon_lr = lonFlip(wind_hgt_canyon_lr)
 tk_wall_lr         = lonFlip(tk_wall_lr)
 tk_roof_lr         = lonFlip(tk_roof_lr)
 tk_improad_lr      = lonFlip(tk_improad_lr)
 cv_wall_lr         = lonFlip(cv_wall_lr)
 cv_roof_lr         = lonFlip(cv_roof_lr)
 cv_improad_lr      = lonFlip(cv_improad_lr)
 thick_wall_lr      = lonFlip(thick_wall_lr)
 thick_roof_lr      = lonFlip(thick_roof_lr)
 nlevimproad_lr     = lonFlip(nlevimproad_lr)
 t_building_min_lr  = lonFlip(t_building_min_lr)
 t_building_max_lr  = lonFlip(t_building_max_lr)

; Determine min/max
 print((/"canyon_hwr      Min: "+min(canyon_hwr_lr(:,:))     +" Max: "+max(canyon_hwr_lr(:,:))/))
 print((/"wtlunit_roof    Min: "+min(wtlunit_roof_lr(:,:))   +" Max: "+max(wtlunit_roof_lr(:,:))/))
 print((/"wtroad_perv     Min: "+min(wtroad_perv_lr(:,:))    +" Max: "+max(wtroad_perv_lr(:,:))/))
 print((/"em_roof         Min: "+min(em_roof_lr(:,:))        +" Max: "+max(em_roof_lr(:,:))/))
 print((/"em_improad      Min: "+min(em_improad_lr(:,:))     +" Max: "+max(em_improad_lr(:,:))/))
 print((/"em_perroad      Min: "+min(em_perroad_lr(:,:))     +" Max: "+max(em_perroad_lr(:,:))/))
 print((/"em_wall         Min: "+min(em_wall_lr(:,:))        +" Max: "+max(em_wall_lr(:,:))/))
 print((/"alb_roof        Min: "+min(alb_roof_lr(:,:,:,:))   +" Max: "+max(alb_roof_lr(:,:,:,:))/))
 print((/"alb_improad     Min: "+min(alb_improad_lr(:,:,:,:))+" Max: "+max(alb_improad_lr(:,:,:,:))/))
 print((/"alb_perroad     Min: "+min(alb_perroad_lr(:,:,:,:))+" Max: "+max(alb_perroad_lr(:,:,:,:))/))
 print((/"alb_wall        Min: "+min(alb_wall_lr(:,:,:,:))   +" Max: "+max(alb_wall_lr(:,:,:,:))/))
 print((/"ht_roof         Min: "+min(ht_roof_lr(:,:))        +" Max: "+max(ht_roof_lr(:,:))/))
 print((/"wind_hgt_canyon Min: "+min(wind_hgt_canyon_lr(:,:))+" Max: "+max(wind_hgt_canyon_lr(:,:))/))
 print((/"tk_wall         Min: "+min(tk_wall_lr(:,:,:))      +" Max: "+max(tk_wall_lr(:,:,:))/))
 print((/"tk_roof         Min: "+min(tk_roof_lr(:,:,:))      +" Max: "+max(tk_roof_lr(:,:,:))/))
 print((/"tk_improad      Min: "+min(tk_improad_lr(:,:,:))   +" Max: "+max(tk_improad_lr(:,:,:))/))
 print((/"cv_wall         Min: "+min(cv_wall_lr(:,:,:))      +" Max: "+max(cv_wall_lr(:,:,:))/))
 print((/"cv_roof         Min: "+min(cv_roof_lr(:,:,:))      +" Max: "+max(cv_roof_lr(:,:,:))/))
 print((/"cv_improad      Min: "+min(cv_improad_lr(:,:,:))   +" Max: "+max(cv_improad_lr(:,:,:))/))
 print((/"thick_wall      Min: "+min(thick_wall_lr(:,:))     +" Max: "+max(thick_wall_lr(:,:))/))
 print((/"thick_roof      Min: "+min(thick_roof_lr(:,:))     +" Max: "+max(thick_roof_lr(:,:))/))
 print((/"nlevimproad     Min: "+min(nlevimproad_lr(:,:))    +" Max: "+max(nlevimproad_lr(:,:))/))
 print((/"t_building_min  Min: "+min(t_building_min_lr(:,:)) +" Max: "+max(t_building_min_lr(:,:))/))
 print((/"t_building_max  Min: "+min(t_building_max_lr(:,:)) +" Max: "+max(t_building_max_lr(:,:))/))

 print( "number of levels of impervious road, and tk for each" );
 do j = 0,nlat_lr-1
 print((/"lat: "+lat_lr(j) /))
 do i = 0,nlon_lr-1
    if (.not.(ismissing(nlevimproad_lr(j,i))))
       if (nlevimproad_lr(j,i) .gt. 0.) then
          print((/nlevimproad_lr(j,i)/))
          print((/tk_improad_lr(0:floattointeger(nlevimproad_lr(j,i))-1,j,i)/))
       end if
    end if
 end do
 end do

 print( "Set attributes" );
 canyon_hwr_lr@long_name      = "canyon height to width ratio"
 canyon_hwr_lr@units          = "unitless"
 wtlunit_roof_lr@long_name    = "fraction of roof"
 wtlunit_roof_lr@units        = "unitless"
 wtroad_perv_lr@long_name     = "fraction of pervious road"
 wtroad_perv_lr@units         = "unitless"
 em_roof_lr@long_name         = "emissivity of roof"
 em_roof_lr@units             = "unitless"
 em_wall_lr@long_name         = "emissivity of wall"
 em_wall_lr@units             = "unitless"
 em_improad_lr@long_name      = "emissivity of impervious road"
 em_improad_lr@units          = "unitless"
 em_perroad_lr@long_name      = "emissivity of pervious road"
 em_perroad_lr@units          = "unitless"
 alb_roof_lr@long_name        = "albedo of roof"
 alb_roof_lr@units            = "unitless"
 alb_wall_lr@long_name        = "albedo of wall"
 alb_wall_lr@units            = "unitless"
 alb_improad_lr@long_name     = "albedo of impervious road"
 alb_improad_lr@units         = "unitless"
 alb_perroad_lr@long_name     = "albedo of pervious road"
 alb_perroad_lr@units         = "unitless"
 ht_roof_lr@long_name         = "height of roof"
 ht_roof_lr@units             = "meters"
 wind_hgt_canyon_lr@long_name = "height of wind in canyon"
 wind_hgt_canyon_lr@units     = "meters"
 nlevimproad_lr@long_name     = "number of impervious road layers"
 nlevimproad_lr@units         = "unitless"
 thick_wall_lr@long_name      = "thickness of wall"
 thick_wall_lr@units          = "m"
 thick_roof_lr@long_name      = "thickness of roof"
 thick_roof_lr@units          = "m"
 t_building_min_lr@long_name  = "minimum interior building temperature"
 t_building_min_lr@units      = "K"
 t_building_max_lr@long_name  = "maximum interior building temperature"
 t_building_max_lr@units      = "K"
 tk_wall_lr@long_name         = "thermal conductivity of wall"
 tk_wall_lr@units             = "W/m*K"
 tk_roof_lr@long_name         = "thermal conductivity of roof"
 tk_roof_lr@units             = "W/m*K"
 cv_wall_lr@long_name         = "volumetric heat capacity of wall"
 cv_wall_lr@units             = "J/m^3*K"
 cv_roof_lr@long_name         = "volumetric heat capacity of roof"
 cv_roof_lr@units             = "J/m^3*K"
 tk_improad_lr@long_name      = "thermal conductivity of impervious road"
 tk_improad_lr@units          = "W/m*K"
 cv_improad_lr@long_name      = "volumetric heat capacity of impervious road"
 cv_improad_lr@units          = "J/m^3*K"

 system("/bin/rm -f " + cdf_filename)
 out_file                  = addfile(cdf_filename,"c")
 out_file@title            = "Averaged urban parameters at "+res+" for TBD,HD, and MD classes"
 out_file@source           = "$HeadURL$"
 out_file@case_id          = "Feddema/Jackson BLDG_TABLE_RevINS_25072008.csv"
 out_file@date             = systemfunc("date")
 out_file->EDGEE           = edgee
 out_file->EDGEN           = edgen
 out_file->EDGES           = edges
 out_file->EDGEW           = edgew
 out_file->LAT             = lat_lr
 out_file->LATIXY          = latixy
 out_file->LON             = lon_lr
 out_file->LONGXY          = longxy
 out_file->CANYON_HWR      = canyon_hwr_lr
 out_file->WTLUNIT_ROOF    = wtlunit_roof_lr
 out_file->WTROAD_PERV     = wtroad_perv_lr
 out_file->EM_ROOF         = em_roof_lr
 out_file->EM_WALL         = em_wall_lr
 out_file->EM_IMPROAD      = em_improad_lr
 out_file->EM_PERROAD      = em_perroad_lr
 out_file->ALB_ROOF        = alb_roof_lr
 out_file->ALB_WALL        = alb_wall_lr
 out_file->ALB_IMPROAD     = alb_improad_lr
 out_file->ALB_PERROAD     = alb_perroad_lr
 out_file->HT_ROOF         = ht_roof_lr
 out_file->WIND_HGT_CANYON = wind_hgt_canyon_lr
 out_file->TK_ROOF         = tk_roof_lr
 out_file->TK_WALL         = tk_wall_lr
 out_file->TK_IMPROAD      = tk_improad_lr
 out_file->CV_ROOF         = cv_roof_lr
 out_file->CV_WALL         = cv_wall_lr
 out_file->CV_IMPROAD      = cv_improad_lr
 out_file->THICK_ROOF      = thick_roof_lr
 out_file->THICK_WALL      = thick_wall_lr
 out_file->NLEV_IMPROAD    = nlevimproad_lr
 out_file->T_BUILDING_MIN  = t_building_min_lr
 out_file->T_BUILDING_MAX  = t_building_max_lr

 if (do_generate_ascii .eq. 0) then
    ; Generate global ascii file.
    tmpoutfile1 = "avg_densityclass_"+res+"_fluxes_p7_c"+sdate+"_1.txt";
    tmpoutfile2 = "avg_densityclass_"+res+"_fluxes_p7_c"+sdate+"_2.txt";
    if(isfilepresent(tmpoutfile1))
       print("Warning: '" + tmpoutfile1 + "' exists. Will remove it.")
       system("/bin/rm " + tmpoutfile1)
    end if
    if(isfilepresent(tmpoutfile2))
       print("Warning: '" + tmpoutfile2 + "' exists. Will remove it.")
       system("/bin/rm " + tmpoutfile2)
    end if
    if(isfilepresent(outfile))
       print("Warning: '" + outfile + "' exists. Will remove it.")
       system("/bin/rm " + outfile)
    end if

    flt_fmt = "%11.3f";
    exp_fmt = "%11.3e";
    print((/"Begin generating ascii file: "+outfile/))
    do j = 0,nlat_lr-1
       print((/"lat: "+lat_lr(j)/))
       do i = 0,nlon_lr-1
          ;print((/"lon: "+i/))
          fg = (/ (/sprintf("canyon_hwr "+flt_fmt,canyon_hwr_lr(j,i))/), \
          (/sprintf("wtlunit_roof "+flt_fmt,wtlunit_roof_lr(j,i))/), \
          (/sprintf("wtroad_perv "+flt_fmt,wtroad_perv_lr(j,i))/), \
          (/sprintf("em_roof "+flt_fmt,em_roof_lr(j,i))/), \
          (/sprintf("em_improad "+flt_fmt,em_improad_lr(j,i))/), \
          (/sprintf("em_perroad "+flt_fmt,em_perroad_lr(j,i))/), \
          (/sprintf("em_wall "+flt_fmt,em_wall_lr(j,i))/), \
          (/sprintf("alb_roof_dir "+flt_fmt,alb_roof_lr(0,0,j,i))+sprintf(flt_fmt,alb_roof_lr(0,1,j,i))/), \
          (/sprintf("alb_roof_dif "+flt_fmt,alb_roof_lr(1,0,j,i))+sprintf(flt_fmt,alb_roof_lr(1,1,j,i))/), \
          (/sprintf("alb_improad_dir "+flt_fmt,alb_improad_lr(0,0,j,i))+sprintf(flt_fmt,alb_improad_lr(0,1,j,i))/), \
          (/sprintf("alb_improad_dif "+flt_fmt,alb_improad_lr(1,0,j,i))+sprintf(flt_fmt,alb_improad_lr(1,1,j,i))/), \
          (/sprintf("alb_perroad_dir "+flt_fmt,alb_perroad_lr(0,0,j,i))+sprintf(flt_fmt,alb_perroad_lr(0,1,j,i))/), \
          (/sprintf("alb_perroad_dif "+flt_fmt,alb_perroad_lr(1,0,j,i))+sprintf(flt_fmt,alb_perroad_lr(1,1,j,i))/), \
          (/sprintf("alb_wall_dir "+flt_fmt,alb_wall_lr(0,0,j,i))+sprintf(flt_fmt,alb_wall_lr(0,1,j,i))/), \
          (/sprintf("alb_wall_dif "+flt_fmt,alb_wall_lr(1,0,j,i))+sprintf(flt_fmt,alb_wall_lr(1,1,j,i))/), \
          (/sprintf("ht_roof "+flt_fmt,ht_roof_lr(j,i))/), \
          (/sprintf("wind_hgt_canyon "+flt_fmt,wind_hgt_canyon_lr(j,i))/), \
          (/sprintf("tk_wall "+flt_fmt,tk_wall_lr(0,j,i))+sprintf(flt_fmt,tk_wall_lr(1,j,i))+ \
          sprintf(flt_fmt,tk_wall_lr(2,j,i))+sprintf(flt_fmt,tk_wall_lr(3,j,i))+ \
          sprintf(flt_fmt,tk_wall_lr(4,j,i))+sprintf(flt_fmt,tk_wall_lr(5,j,i))+ \
          sprintf(flt_fmt,tk_wall_lr(6,j,i))+sprintf(flt_fmt,tk_wall_lr(7,j,i))+ \
          sprintf(flt_fmt,tk_wall_lr(8,j,i))+sprintf(flt_fmt,tk_wall_lr(9,j,i))/), \
          (/sprintf("tk_roof "+flt_fmt,tk_roof_lr(0,j,i))+sprintf(flt_fmt,tk_roof_lr(1,j,i))+ \
          sprintf(flt_fmt,tk_roof_lr(2,j,i))+sprintf(flt_fmt,tk_roof_lr(3,j,i))+ \
          sprintf(flt_fmt,tk_roof_lr(4,j,i))+sprintf(flt_fmt,tk_roof_lr(5,j,i))+ \
          sprintf(flt_fmt,tk_roof_lr(6,j,i))+sprintf(flt_fmt,tk_roof_lr(7,j,i))+ \
          sprintf(flt_fmt,tk_roof_lr(8,j,i))+sprintf(flt_fmt,tk_roof_lr(9,j,i))/), \
          (/sprintf("tk_improad "+flt_fmt,tk_improad_lr(0,j,i))+sprintf(flt_fmt,tk_improad_lr(1,j,i))+ \
          sprintf(flt_fmt,tk_improad_lr(2,j,i))+sprintf(flt_fmt,tk_improad_lr(3,j,i))+ \
          sprintf(flt_fmt,tk_improad_lr(4,j,i))+sprintf(flt_fmt,tk_improad_lr(5,j,i))+ \
          sprintf(flt_fmt,tk_improad_lr(6,j,i))+sprintf(flt_fmt,tk_improad_lr(7,j,i))+ \
          sprintf(flt_fmt,tk_improad_lr(8,j,i))+sprintf(flt_fmt,tk_improad_lr(9,j,i))/), \
          (/sprintf("cv_wall "+exp_fmt,cv_wall_lr(0,j,i))+sprintf(exp_fmt,cv_wall_lr(1,j,i))+ \
          sprintf(exp_fmt,cv_wall_lr(2,j,i))+sprintf(exp_fmt,cv_wall_lr(3,j,i))+ \
          sprintf(exp_fmt,cv_wall_lr(4,j,i))+sprintf(exp_fmt,cv_wall_lr(5,j,i))+ \
          sprintf(exp_fmt,cv_wall_lr(6,j,i))+sprintf(exp_fmt,cv_wall_lr(7,j,i))+ \
          sprintf(exp_fmt,cv_wall_lr(8,j,i))+sprintf(exp_fmt,cv_wall_lr(9,j,i))/), \
          (/sprintf("cv_roof "+exp_fmt,cv_roof_lr(0,j,i))+sprintf(exp_fmt,cv_roof_lr(1,j,i))+ \
          sprintf(exp_fmt,cv_roof_lr(2,j,i))+sprintf(exp_fmt,cv_roof_lr(3,j,i))+ \
          sprintf(exp_fmt,cv_roof_lr(4,j,i))+sprintf(exp_fmt,cv_roof_lr(5,j,i))+ \
          sprintf(exp_fmt,cv_roof_lr(6,j,i))+sprintf(exp_fmt,cv_roof_lr(7,j,i))+ \
          sprintf(exp_fmt,cv_roof_lr(8,j,i))+sprintf(exp_fmt,cv_roof_lr(9,j,i))/), \
          (/sprintf("cv_improad "+exp_fmt,cv_improad_lr(0,j,i))+sprintf(exp_fmt,cv_improad_lr(1,j,i))+ \
          sprintf(exp_fmt,cv_improad_lr(2,j,i))+sprintf(exp_fmt,cv_improad_lr(3,j,i))+ \
          sprintf(exp_fmt,cv_improad_lr(4,j,i))+sprintf(exp_fmt,cv_improad_lr(5,j,i))+ \
          sprintf(exp_fmt,cv_improad_lr(6,j,i))+sprintf(exp_fmt,cv_improad_lr(7,j,i))+ \
          sprintf(exp_fmt,cv_improad_lr(8,j,i))+sprintf(exp_fmt,cv_improad_lr(9,j,i))/), \
          (/sprintf("thick_wall "+flt_fmt,thick_wall_lr(j,i))/), \
          (/sprintf("thick_roof "+flt_fmt,thick_roof_lr(j,i))/), \
          (/sprinti("nlev_improad %2i",floattointeger(nlevimproad_lr(j,i)))/), \
          (/sprintf("t_building_min "+flt_fmt,t_building_min_lr(j,i))/), \
          (/sprintf("t_building_max "+flt_fmt,t_building_max_lr(j,i))/)/)
          if (j .eq. 0 .and. i .eq. 0) then
             asciiwrite(tmpoutfile1,fg)
          else
             asciiwrite(tmpoutfile2,fg)
             system("cat "+tmpoutfile1+" "+tmpoutfile2+">"+" "+outfile)
             system("mv "+outfile+" "+tmpoutfile1)
          end if
       end do  
    end do
    system("mv "+tmpoutfile1+" "+outfile)
    system("rm "+tmpoutfile2)

 end if

 print ("=========================================")
 print ("Finish Time: "+systemfunc("date") )
 print ("=========================================")


  print( "============================================================================");
  print( "Successfully created output urban input file: "+outfile);

end
