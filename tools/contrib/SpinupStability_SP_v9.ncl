; NCL script
; SpinupStability_SP_v9.ncl
; Script to examine stability of SP spinup simulation.
; This version operates on either monthly mean or multi-annual mean multi-variable history files 
; Keith Oleson, Sep 2018

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "~oleson/lnd_diag/run/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  print ("=========================================")
  print ("Start Time: "+systemfunc("date") )
  print ("=========================================")

;=======================================================================;
; This script produces a page of plots of various variables that are evaluated
; as to whether they are spunup or not.  A summary of variables in equilibrium
; and/or in disequilibrium is also produced to standard out. The postscript output
; is $caseid_spinup.ps in the pwd.
; The variables examined are FSH,EFLX_LH_TOT,FPSN,TWS,H2OSOI,TSOI.
; The percentage of land area in TSOI disequilibrium is also examined.
; Thresholds are defined below (i.e., global_thresh_*) and can be changed for individual
; variables.
; To run this script, just enter in your case name and username below. 
;  AND set the annual_hist flag to "True" if your case has annual mean history files or set 
;  annual_hist flag to "False" if your case has monthly mean history files.
; The script ; assumes that your history files are in /glade/scratch/$username/archive/$caseid/lnd/hist
;=======================================================================;

  caseid = "cesm20exp10j_1deg_CPLHST_nndep_1850pAD"
  username = "oleson"
  annual_hist = "False"
  cplot = "Global"
  subper = 20
  h2osoi_layer = 8               ; Desired soil layer (layer 8 is about 1m)
  tsoi_layer = 10                ; Desired soil layer (layer 10 is about 3m)

  do_plot = "True"
;=======================================================================;

  data_dir = "/glade/scratch/"+username+"/archive/"+caseid+"/lnd/hist/"

; Thresholds 
  glob_thresh_fsh = 0.02         ; global threshold for FSH equilibrium (delta W m-2 / yr)
  glob_thresh_lh = 0.02          ; global threshold for EFLX_LH_TOT equilibrium (delta W m-2 / yr)
  glob_thresh_gpp = 0.02         ; global threshold for FPSN equilibrium (delta PgC / yr)
  glob_thresh_tws = 0.001        ; global threshold for TWS equilibrium (delta m / yr)
  glob_thresh_h2osoi = 0.01      ; global threshold for H2OSOI equilibrium (delta mm mm-3 / yr)
  glob_thresh_tsoi = 0.02        ; global threshold for TSOI equilibrium (delta K / yr)
  glob_thresh_area = 3.0         ; global threshold percent area with TWS disequilibrium gt 0.01 m
  tws_thresh = 0.001             ; disequilibrium threshold for individual gridcells (m)

  if (annual_hist .eq. "True") then
    fls                       = systemfunc("ls " + data_dir + caseid+".clm2.h0.*-*-*-*"+".nc")
  else
    fls                       = systemfunc("ls " + data_dir + caseid+".clm2.h0.*-*"+".nc")
  end if
  flsdims                   = dimsizes(fls)

  if (annual_hist .eq. "True") then
    lstfile                   = addfile(fls(flsdims-1),"r")
  else
    lstfile                   = addfile(fls(flsdims-12),"r")  ;last month (DEC) of any year has mdate for next year
                                                              ;so grab JAN of last year
  end if

  if (annual_hist .eq. "True") then
    lstyrdim                = dimsizes(lstfile->mcdate)
    mcdate_lst              = lstfile->mcdate(lstyrdim-1)
  else
    mcdate_lst              = lstfile->mcdate
  end if
  lstyr                     = toint(mcdate_lst)/10000

  fstfile                   = addfile(fls(0),"r")
  if (annual_hist .eq. "True") then
    mcdate_fst              = fstfile->mcdate(0)
  else
    mcdate_fst              = fstfile->mcdate
  end if
  fstyr                     = toint(mcdate_fst)/10000

  yearplt                   = ispan(fstyr,lstyr,subper)
  yearpltrev                = yearplt(::-1)
  year                      = ispan(fstyr,lstyr,subper) - fstyr
  nyrs                      = dimsizes(year)
  if (subper .eq. 1) then
    yearpltmid                = ispan(fstyr+subper/2,yearplt(nyrs-2),subper)
  else
    yearpltmid                = ispan(fstyr+subper/2,yearplt(nyrs-1),subper)
  end if

; Build an array of monthly indices
  monstr                    = new(nyrs*12,"integer")
  do i = 0,nyrs-1
    monstr(i*12:i*12+11)    = ispan(year(i)*12,year(i)*12+11,1)
  end do

; Add the data files together
  if (annual_hist .eq. "True") then
    data                    = addfiles(fls, "r")
    ListSetType (data, "cat")
  else
    data                    = addfiles(fls(monstr), "r")
  end if

; Convert to annual means
  if (annual_hist .eq. "True") then
    fsh              = data[:]->FSH(year,:,:)
    lh               = data[:]->EFLX_LH_TOT(year,:,:)
    gpp              = data[:]->FPSN(year,:,:)
    tws              = data[:]->TWS(year,:,:)
    h2osoi           = data[:]->H2OSOI(year,h2osoi_layer-1,:,:)
    tsoi             = data[:]->TSOI(year,tsoi_layer-1,:,:)
  else
    fsh              = month_to_annual(data[:]->FSH,1)
    if (isfilevar(data[0],"EFLX_LH_TOT")) then
      lh               = month_to_annual(data[:]->EFLX_LH_TOT,1)
    else
      tmp1           = data[:]->FCTR
      tmp2           = data[:]->FCEV
      tmp3           = data[:]->FGEV
      tot            = tmp1+tmp2+tmp3
      copy_VarCoords(tmp1,tot)
      lh             = month_to_annual(tot,1)
      delete(tmp1)
      delete(tmp2)
      delete(tmp3)
      delete(tot)
    end if
    gpp              = month_to_annual(data[:]->FPSN,1)
    if (isfilevar(data[0],"TWS")) then
      tws              = month_to_annual(data[:]->TWS,1)
    else
      tmp1             = data[:]->H2OCAN
      tmp2             = data[:]->H2OSNO
      tmp3             = data[:]->WA
      tmp4             = data[:]->SOILLIQ
      tmp4s            = dim_sum_n(tmp4,1)
      tmp5             = data[:]->SOILICE
      tmp5s            = dim_sum_n(tmp5,1)
      tot              = tmp1+tmp2+tmp3+tmp4s+tmp5s
      copy_VarCoords(tmp1,tot)
      tws              = month_to_annual(tot,1)
      delete(tmp1)
      delete(tmp2)
      delete(tmp3)
      delete(tmp4)
      delete(tmp4s)
      delete(tmp5)
      delete(tmp5s)
    end if
    h2osoi           = month_to_annual(data[:]->H2OSOI(:,h2osoi_layer-1,:,:),1)
    tsoi             = month_to_annual(data[:]->TSOI(:,tsoi_layer-1,:,:),1)
  end if
  lat                       = data[0]->lat
  nlat                      = dimsizes(lat)
  lon                       = data[0]->lon
  nlon                      = dimsizes(lon)
  landfrac                  = data[0]->landfrac
  area                      = data[0]->area
; nbedrock                  = data[0]->nbedrock
  aream                     = area*1.e6
  landarea                  = landfrac*aream
  gtoPg                     = 1e-15
  secinyr                   = 60.*60.*24.*365.
  umoleCO2tomoleCO2         = 1./1.e6
  moleCO2togCO2             = 44./1.
  gCO2togC                  = 12./44.

; FSH
  if (isfilevar(data[0],"FSH")) then

  areasum            = sum(area*landfrac)
  areaL              = area*landfrac
  landareaL          = conform_dims(dimsizes(fsh),areaL,(/1,2/))  ; conforming dimensions of areaL to tlai
  fsh_glob           = dim_sum_n(fsh*landareaL/areasum,(/1,2/))   ; weighted global fsh (W m-2)
  fsh_glob!0         = "year"                                                
  fsh_glob&year      = yearplt                                               
  fsh_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    fsh_glob_del(i) = (fsh_glob(i+1) - fsh_glob(i))/subper
  end do
  fsh_glob_del!0     = "year"
  fsh_glob_del&year  = yearpltmid
  indx = where(abs(fsh_glob_del) .lt. glob_thresh_fsh,1,0)
  if (all(indx .eq. 1)) then
    fsh_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          fsh_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      fsh_glob_equil   = -999
    end if
  end if
  fsh_glob_equil@_FillValue = -999
  delete(indx)

  else
    fsh_glob_equil = new(1,typeof(yearplt),-999)
    fsh_glob_equil = -999
    fsh_glob       = new(dimsizes(dim_sum_n(fsh,(/1,2/))),typeof(fsh),-999)
    fsh_glob       = -999
    fsh_glob_del   = new(nyrs-1,"float",-999.)
    fsh_glob_del   = -999.
  end if

; EFLX_LH_TOT
; if (isfilevar(data[0],"EFLX_LH_TOT")) then

  lh_glob           = dim_sum_n(lh*landareaL/areasum,(/1,2/))    ; weighted global lh (W m-2)
  lh_glob!0         = "year"                                                
  lh_glob&year      = yearplt                                               
  lh_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    lh_glob_del(i) = (lh_glob(i+1) - lh_glob(i))/subper
  end do
  lh_glob_del!0     = "year"
  lh_glob_del&year  = yearpltmid
  indx = where(abs(lh_glob_del) .lt. glob_thresh_lh,1,0)
  if (all(indx .eq. 1)) then
    lh_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          lh_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      lh_glob_equil   = -999
    end if
  end if
  lh_glob_equil@_FillValue = -999
  delete(indx)

; else
;   lh_glob_equil = new(1,typeof(yearplt),-999)
;   lh_glob_equil = -999
;   lh_glob       = new(dimsizes(dim_sum_n(lh,(/1,2/))),typeof(lh),-999)
;   lh_glob       = -999
;   lh_glob_del   = new(nyrs-1,"float",-999.)
;   lh_glob_del   = -999.
; end if

; GPP
  landareaC             = conform_dims(dimsizes(gpp),landarea,(/1,2/))                    ; conforming dimensions of landarea to gpp
  gpp_area              = gpp*landareaC                                                   ; correcting gpp for total land area
  gpp_pg                = gpp_area*gtoPg*secinyr*umoleCO2tomoleCO2*moleCO2togCO2*gCO2togC ; umole CO2/s->Pg Carbon
  gpp_glob              = dim_sum_n(gpp_pg, (/1,2/))                                      ; sums over all latitudes
  gpp_glob!0            = "year"                                                
  gpp_glob&year         = yearplt                                               
  gpp_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    gpp_glob_del(i) = (gpp_glob(i+1) - gpp_glob(i))/subper
  end do
  gpp_glob_del!0     = "year"
  gpp_glob_del&year  = yearpltmid
  indx = where(abs(gpp_glob_del) .lt. glob_thresh_gpp,1,0)
  if (all(indx .eq. 1)) then
    gpp_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          gpp_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      gpp_glob_equil   = -999
    end if
  end if
  gpp_glob_equil@_FillValue = -999
  delete(indx)

; TWS
  tws_glob           = (dim_sum_n(tws*landareaL/areasum,(/1,2/)))/1.e3    ; weighted global tws (meters)
  tws_glob!0         = "year"                                                
  tws_glob&year      = yearplt                                               
  tws_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    tws_glob_del(i) = (tws_glob(i+1) - tws_glob(i))/subper
  end do
  tws_glob_del!0     = "year"
  tws_glob_del&year  = yearpltmid
  indx = where(abs(tws_glob_del) .lt. glob_thresh_tws,1,0)
  if (all(indx .eq. 1)) then
    tws_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          tws_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      tws_glob_equil   = -999
    end if
  end if
  tws_glob_equil@_FillValue = -999
  delete(indx)

; H2OSOI 
  if (isfilevar(data[0],"H2OSOI")) then

  h2osoi_glob           = dim_sum_n(h2osoi*landareaL/areasum,(/1,2/))    ; weighted global h2osoi (mm)
  h2osoi_glob!0         = "year"                                                
  h2osoi_glob&year      = yearplt                                               
  h2osoi_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    h2osoi_glob_del(i) = (h2osoi_glob(i+1) - h2osoi_glob(i))/subper
  end do
  h2osoi_glob_del!0     = "year"
  h2osoi_glob_del&year  = yearpltmid
  indx = where(abs(h2osoi_glob_del) .lt. glob_thresh_h2osoi,1,0)
  if (all(indx .eq. 1)) then
    h2osoi_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          h2osoi_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      h2osoi_glob_equil   = -999
    end if
  end if
  h2osoi_glob_equil@_FillValue = -999
  delete(indx)

  else
    h2osoi_glob_equil = -999
    h2osoi_glob_equil = new(1,typeof(yearplt),-999)
    h2osoi_glob_equil = -999
    h2osoi_glob       = new(dimsizes(dim_sum_n(h2osoi,(/1,2/))),typeof(h2osoi),-999)
    h2osoi_glob       = -999
    h2osoi_glob_del   = new(nyrs-1,"float",-999.)
    h2osoi_glob_del   = -999.
  end if

; TSOI 
  if (isfilevar(data[0],"TSOI")) then

  tsoi_glob           = dim_sum_n(tsoi*landareaL/areasum,(/1,2/))    ; weighted global tsoi (mm)
  tsoi_glob!0         = "year"                                                
  tsoi_glob&year      = yearplt                                               
  tsoi_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    tsoi_glob_del(i) = (tsoi_glob(i+1) - tsoi_glob(i))/subper
  end do
  tsoi_glob_del!0     = "year"
  tsoi_glob_del&year  = yearpltmid
  indx = where(abs(tsoi_glob_del) .lt. glob_thresh_tsoi,1,0)
  if (all(indx .eq. 1)) then
    tsoi_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          tsoi_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      tsoi_glob_equil   = -999
    end if
  end if
  tsoi_glob_equil@_FillValue = -999
  delete(indx)

  else
    h2osno_glob_equil = -999
    h2osno_glob_equil = new(1,typeof(yearplt),-999)
    h2osno_glob_equil = -999
    h2osno_glob       = new(dimsizes(dim_sum_n(h2osno,(/1,2/))),typeof(h2osno),-999)
    h2osno_glob       = -999
    h2osno_glob_del   = new(nyrs-1,"float",-999.)
    h2osno_glob_del   = -999.
  end if

; Land area not in TWS equilibrium
  landarea_noequil = new((/nyrs-1,nlat,nlon/),"float")
  do i = 0,nyrs-2
    landarea_noequil(i,:,:) = where(abs(tws(i+1,:,:)/1.e3 - tws(i,:,:)/1.e3)/subper .gt. tws_thresh, landarea, 0.)
  end do
  perc_landarea_noequil = 100.*(dim_sum_n(landarea_noequil,(/1,2/))/sum(landarea))
  indx = where(abs(perc_landarea_noequil) .lt. glob_thresh_area,1,0)
  if (all(indx .eq. 1)) then
    perc_landarea_glob_noequil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          perc_landarea_glob_noequil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      perc_landarea_glob_noequil   = -999
    end if
  end if
  perc_landarea_glob_noequil@_FillValue = -999
  delete(indx)
  tws_1_map = where(landarea_noequil(nyrs-2,:,:) .ne. 0.,(tws(nyrs-1,:,:)/1.e3-tws(nyrs-2,:,:)/1.e3)/subper,tws@_FillValue)
  tws_1_map!0 = "lat"
  tws_1_map&lat = lat
  tws_1_map!1 = "lon"
  tws_1_map&lon = lon
  tws_2_map = where(landarea_noequil(nyrs-3,:,:) .ne. 0.,(tws(nyrs-2,:,:)/1.e3-tws(nyrs-3,:,:)/1.e3)/subper,tws@_FillValue)
  copy_VarCoords(tws_1_map,tws_2_map)

;===============================Plotting====================================;
  if (do_plot .eq. "True") then

; wks_type = "png"
; wks_type@wkWidth  = 2500
; wks_type@wkHeight = 2500
; wks_type = "x11"
  wks_type = "ps"
; wks_type = "pdf"
  wks                        = gsn_open_wks (wks_type, caseid+"_SP_Spinup")
  gsn_define_colormap(wks, "ViBlGrWhYeOrRe")

  plot = new(13,"graphic")

  resP = True
; resP@gsnMaximize = True
  resP@gsnPaperOrientation = "portrait"
  resP@gsnFrame = False
  resP@gsnDraw = True
  resP@txString = caseid + " Annual Mean "
  resP@gsnPanelXWhiteSpacePercent = 2.
  resP@gsnPanelCenter = False
; resP@gsnPanelDebug = True

  res                        = True
  res@gsnFrame               = False
  res@gsnDraw                = False
  res@xyLineThicknessF       = 2

  polyres1 = True
  polyres1@gsLineDashPattern = 16
  polyres2 = True
  polyres2@gsLineDashPattern = 2

  res@tiXAxisString          = " "
  res@tiYAxisString          = "W m~S~-2~N~"
  res@tiMainString           = "FSH"
  plot(0)                    = gsn_csm_xy(wks,yearplt,fsh_glob,res)

  res@tiYAxisString          = "W m~S~-2~N~"
  res@tiMainString           = "Delta FSH " + "EqYr: "+fsh_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(1)                    = gsn_csm_xy(wks,yearpltmid,fsh_glob_del,res)
  prim1 = gsn_add_polyline(wks,plot(1),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim2 = gsn_add_polyline(wks,plot(1),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_fsh,-glob_thresh_fsh/),polyres2)
  prim3 = gsn_add_polyline(wks,plot(1),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_fsh,glob_thresh_fsh/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "W m~S~-2~N~"
  res@tiMainString           = "EFLX_LH_TOT"
  plot(2)                    = gsn_csm_xy(wks,yearplt,lh_glob,res)

  res@tiYAxisString          = "W m~S~-2~N~"
  res@tiMainString           = "Delta EFLX_LH_TOT " + "EqYr: "+lh_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(3)                    = gsn_csm_xy(wks,yearpltmid,lh_glob_del,res)
  prim4 = gsn_add_polyline(wks,plot(3),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim5 = gsn_add_polyline(wks,plot(3),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_lh,-glob_thresh_lh/),polyres2)
  prim6 = gsn_add_polyline(wks,plot(3),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_lh,glob_thresh_lh/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "Pg C yr~S~-1~N~"
  res@tiMainString           = "GPP"
  plot(4)                    = gsn_csm_xy(wks,yearplt,gpp_glob,res)

  res@tiYAxisString          = "Pg C yr~S~-1~N~"
  res@tiXAxisString          = "Spinup Year"
  res@tiMainString           = "Delta GPP " + "EqYr: "+gpp_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(5)                    = gsn_csm_xy(wks,yearpltmid,gpp_glob_del,res)
  prim7 = gsn_add_polyline(wks,plot(5),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim8 = gsn_add_polyline(wks,plot(5),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_gpp,-glob_thresh_gpp/),polyres2)
  prim9 = gsn_add_polyline(wks,plot(5),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_gpp,glob_thresh_gpp/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "m"
  res@tiMainString           = "TWS"
  plot(6)                   = gsn_csm_xy(wks,yearplt,tws_glob,res)

  res@tiYAxisString          = "m"
  res@tiMainString           = "Delta TWS " + "EqYr: "+tws_glob_equil
  res@trYMaxF                = 0.05
  res@trYMinF                = -0.05
  plot(7)                   = gsn_csm_xy(wks,yearpltmid,tws_glob_del,res)
  prim10 = gsn_add_polyline(wks,plot(7),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim11 = gsn_add_polyline(wks,plot(7),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_tws,-glob_thresh_tws/),polyres2)
  prim12 = gsn_add_polyline(wks,plot(7),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_tws,glob_thresh_tws/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "mm~S~3~N~ mm~S~3~N~"
  res@tiMainString           = "H2OSOI (layer "+h2osoi_layer+")"
  plot(8)                    = gsn_csm_xy(wks,yearplt,h2osoi_glob,res)

  res@tiYAxisString          = "mm~S~3~N~ mm~S~3~N~"
  res@tiXAxisString          = "Spinup Year"
  res@tiMainString           = "Delta H2OSOI " + "EqYr: "+h2osoi_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(9)                    = gsn_csm_xy(wks,yearpltmid,h2osoi_glob_del,res)
  prim13 = gsn_add_polyline(wks,plot(9),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim14 = gsn_add_polyline(wks,plot(9),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_h2osoi,-glob_thresh_h2osoi/),polyres2)
  prim15 = gsn_add_polyline(wks,plot(9),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_h2osoi,glob_thresh_h2osoi/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "K"
  res@tiMainString           = "TSOI (layer "+tsoi_layer+")"
  plot(10)                    = gsn_csm_xy(wks,yearplt,tsoi_glob,res)

  res@tiYAxisString          = "K"
  res@tiXAxisString          = "Spinup Year"
  res@tiMainString           = "Delta TSOI " + "EqYr: "+tsoi_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(11)                    = gsn_csm_xy(wks,yearpltmid,tsoi_glob_del,res)
  prim16 = gsn_add_polyline(wks,plot(11),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim17 = gsn_add_polyline(wks,plot(11),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_tsoi,-glob_thresh_tsoi/),polyres2)
  prim18 = gsn_add_polyline(wks,plot(11),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_tsoi,glob_thresh_tsoi/),polyres2)

  res@tiYAxisString          = "%"
  res@tiMainString           = "% Land Area in TWS Disequil " + "EqYr: "+perc_landarea_glob_noequil
  res@trYMaxF                = 80.0
  res@trYMinF                =  0.0
  plot(12)                   = gsn_csm_xy(wks,yearpltmid,perc_landarea_noequil,res)
  prim19 = gsn_add_polyline(wks,plot(12),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_area,glob_thresh_area/),polyres2)

  gsn_panel(wks,plot,(/4,4/),resP)

  delete(plot)
  resc                       = True             ; turn on plotting options
  resc@gsnSpreadColors       = True             ; spans all colors in colormap
  resc@cnFillMode            = "RasterFill"     ; raster mode
  resc@cnFillOn              = True             ; turn on color fill
  resc@cnLinesOn             = False            ; turn off contour lines
  resc@cnLineLabelsOn        = False            ; turn off contour line labels
  resc@cnLevelSelectionMode  = "ExplicitLevels" 
  resc@mpProjection          = "robinson"       ; Robinson grid projection
  if (cplot .eq. "Arctic") then
    resc@mpLimitMode   = "LatLon"
    resc@mpMinLatF     = 50.
    resc@mpMaxLatF     = 85.
    resc@mpMinLonF     = -180.
    resc@mpMaxLonF     = -95.
    resc@gsnAddCyclic  = False
    resc@mpProjection  = "CylindricalEquidistant"
  end if
  resc@gsnDraw               = True
  resc@gsnFrame              = False
  resc@lbAutoManage          = False
  resc@lbLabelFontHeightF    = 0.010
  resc@txFontHeightF         = 0.008
  resc@cnLevels              = (/-0.005,-0.004,-0.003,-0.002,-0.001,0.0,0.001,0.002,0.003,0.004,0.005/)
  resc@gsnLeftString         = "m"
  resc@gsnRightString        = ""

  resc@vpXF                  = 0.30                  ; position and sizes
  resc@vpYF                  = 0.28                  ; for XY plot
  resc@vpWidthF              = 0.30
  resc@vpHeightF             = 0.30
  resc@gsnCenterString       = "TWS Disequil Yr " + yearplt(nyrs-1) + " - " + yearplt(nyrs-2)
  plot                       = gsn_csm_contour_map(wks,tws_1_map,resc)

  resc@vpXF                  = 0.65                  ; position and sizes
  resc@vpYF                  = 0.28                  ; for XY plot
  resc@vpWidthF              = 0.30
  resc@vpHeightF             = 0.30
  resc@gsnCenterString       = "TWS Disequil Yr " + yearplt(nyrs-2) + " - " + yearplt(nyrs-3)
  plot                       = gsn_csm_contour_map(wks,tws_2_map,resc)

  frame(wks)

  end if   ; end do_plot

; Equilibrium summary
  print((/"======================================================================="/))
  print((/"======================================================================="/))
  print((/"EQUILIBRIUM SUMMARY"/))
  print((/"======================================================================="/))
  if (.not.(ismissing(fsh_glob_equil))) then
    print((/"FSH is in equilibrium. Eq. Yr. = "+fsh_glob_equil/))
  else
    print((/"FATAL: FSH is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(lh_glob_equil))) then
    print((/"EFLX_LH_TOT is in equilibrium. Eq. Yr. = "+lh_glob_equil/))
  else
    print((/"FATAL: EFLX_LH_TOT is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(gpp_glob_equil))) then
    print((/"GPP is in equilibrium. Eq. Yr. = "+gpp_glob_equil/))
  else
    print((/"FATAL: GPP is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(tws_glob_equil))) then
    print((/"TWS is in equilibrium. Eq. Yr. = "+tws_glob_equil/))
  else
    print((/"WARNING: TWS is NOT in equilibrium or is missing"/))
  end if
  if (.not.(ismissing(h2osoi_glob_equil))) then
    print((/"H2OSOI is in equilibrium. Eq. Yr. = "+h2osoi_glob_equil/))
  else
    print((/"FATAL: H2OSOI is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(tsoi_glob_equil))) then
    print((/"TSOI is in equilibrium. Eq. Yr. = "+tsoi_glob_equil/))
  else
    print((/"FATAL: TSOI is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(perc_landarea_glob_noequil))) then
    print((/"At least "+(100.-glob_thresh_area)+" percent of the land surface is in TWS equilibrium. Eq. Yr. = "+perc_landarea_glob_noequil/))
    print((/"Percent of the land surface not in equilibrium ("+sprintf("%6.2f",perc_landarea_noequil(nyrs-2))+"% )"/))
  else
    print((/"FATAL: Not enough of the land surface is in equilibrium ("+sprintf("%6.2f",perc_landarea_noequil(nyrs-2))+"% > "+sprintf("%6.2f",glob_thresh_area)+"%)"/))
  end if
  if (.not.(ismissing(fsh_glob_equil)) .and. \
      .not.(ismissing(lh_glob_equil))    .and. \
      .not.(ismissing(gpp_glob_equil))        .and. \
      .not.(ismissing(tws_glob_equil))        .and. \
      .not.(ismissing(h2osoi_glob_equil))        .and. \
      .not.(ismissing(tsoi_glob_equil))        .and. \
      .not.(ismissing(perc_landarea_glob_noequil))) then
    print((/"Congratulations! Your simulation is in equilibrium"/))
  else
    print((/"FATAL: Your simulation is not in equilibrium, 8 hours have been deducted from your PTO bank, try again"/))
  end if
  print((/"======================================================================="/))

  print ("=========================================")
  print ("Finish Time: "+systemfunc("date") )
  print ("=========================================")
  print ("Successfully ran the script")

end
