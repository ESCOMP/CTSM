; NCL script
; SpinupStability_v6.ncl
; Script to examine stability of spinup simulation.
; This version operates on either monthly mean or multi-annual mean multi-variable history files 
; Keith Oleson, May 2016

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
; The variables examined are TOTECOSYSC,TOTSOMC,TOTVEGC,TLAI,GPP,TWS.
; The percentage of land area in TOTECOSYSC disequilibrium is also examined.
; Thresholds are defined below (i.e., global_thresh_*) and can be changed for individual
; variables.
; To run this script, just enter in your case name and username below. 
;  AND set the annual_hist flag to "True" if your case has annual mean history files or set 
;  annual_hist flag to "False" if your case has monthly mean history files.
; The script ; assumes that your history files are in /glade/scratch/$username/archive/$caseid/lnd/hist
;=======================================================================;

  caseid = "clm5r221_2deg_calpv9_savetandg_1850ADspin"
  username = "oleson"
  annual_hist = "False"

; caseid = "clm5r221_2deg_calpv9_savetandg_c3c4_1850ADspin"
; username = "oleson"
; annual_hist = "False"

  do_plot = "True"
; do_plot = "False"
;=======================================================================;

  data_dir = "/glade/scratch/"+username+"/archive/"+caseid+"/lnd/hist/"
; data_dir = "/glade/scratch/"+username+"/"+caseid+"/run/"
; data_dir = "/glade/p/cesm/palwg_dev/pliocene/archive/F.e15.cam54_clm45_plio_fire_plio/lnd/hist/"
; data_dir = "/glade/p/cesm/palwg_dev/pliocene/archive/BGC_final_spinup_plio/lnd/hist/"

  subper = 20                    ; Subsampling period in years
; subper = 30                    ; Subsampling period in years
; subper = 5                     ; Subsampling period in years

; Thresholds 
  glob_thresh_totecosysc = 0.02  ; global threshold for TOTECOSYSC equilibrium (delta PgC / yr)
  glob_thresh_totsomc = 0.02     ; global threshold for TOTSOMC equilibrium (delta PgC / yr)
  glob_thresh_totvegc = 0.02     ; global threshold for TOTVEGC equilibrium (delta PgC / yr)
  glob_thresh_tlai = 0.02        ; global threshold for TLAI equilibrium (delta m2/m2 / yr)
  glob_thresh_gpp = 0.02         ; global threshold for GPP equilibrium (delta PgC / yr)
  glob_thresh_tws = 0.001        ; global threshold for TWS equilibrium (delta m / yr)
  glob_thresh_area = 3.0         ; global threshold percent area with TOTECOSYSC disequilibrium gt 1 g C/m2/yr
  totecosysc_thresh = 1.         ; disequilibrium threshold for individual gridcells (gC/m2/yr)

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
  yearpltmid                = ispan(fstyr+subper/2,yearplt(nyrs-1),subper)

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
    totecosysc              = data[:]->TOTECOSYSC(year,:,:)
    totsomc                 = data[:]->TOTSOMC(year,:,:)
    totvegc                 = data[:]->TOTVEGC(year,:,:)
    tlai                    = data[:]->TLAI(year,:,:)
    gpp                     = data[:]->GPP(year,:,:)
    tws                     = data[:]->TWS(year,:,:)
;   fert                    = data[:]->FERT(year,:,:)
  else
    totecosysc              = month_to_annual(data[:]->TOTECOSYSC,1)
    totsomc                 = month_to_annual(data[:]->TOTSOMC,1)
    totvegc                 = month_to_annual(data[:]->TOTVEGC,1)
    tlai                    = month_to_annual(data[:]->TLAI,1)
    gpp                     = month_to_annual(data[:]->GPP,1)
    tws                     = month_to_annual(data[:]->TWS,1)
;   fert                    = month_to_annual(data[:]->FERT,1)
  end if
  lat                       = data[0]->lat
  nlat                      = dimsizes(lat)
  lon                       = data[0]->lon
  nlon                      = dimsizes(lon)
  landfrac                  = data[0]->landfrac
  area                      = data[0]->area
  aream                     = area*1.e6
  landarea                  = landfrac*aream
  gtoPg                     = 1e-15
  secinyr                   = 60.*60.*24.*365.

; TOTECOSYSC
  landareaC                 = conform_dims(dimsizes(totecosysc),landarea,(/1,2/)) ; conforming dimensions of landarea to totecosysc
  totecosysc_area           = totecosysc*landareaC                                ; correcting totecosysc for total land area
  totecosysc_pg             = totecosysc_area*gtoPg                               ; g to Pg
  totecosysc_glob           = dim_sum_n(totecosysc_pg, (/1,2/))                   ; sums over all latitudes
  totecosysc_glob!0         = "year"                                              
  totecosysc_glob&year      = yearplt                                             
  totecosysc_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    totecosysc_glob_del(i)  = (totecosysc_glob(i+1) - totecosysc_glob(i))/subper
  end do
  totecosysc_glob_del!0     = "year"
  totecosysc_glob_del&year  = yearpltmid
  indx = where(abs(totecosysc_glob_del) .lt. glob_thresh_totecosysc,1,0)
  if (all(indx .eq. 1)) then
    totecosysc_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          totecosysc_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      totecosysc_glob_equil   = -999
    end if
  end if
  totecosysc_glob_equil@_FillValue = -999
  delete(indx)

; Land area not in TOTECOSYSC equilibrium
  landarea_noequil = new((/nyrs-1,nlat,nlon/),"float")
  do i = 0,nyrs-2
    landarea_noequil(i,:,:) = where((totecosysc(i+1,:,:) - totecosysc(i,:,:))/subper .gt. totecosysc_thresh, landarea, 0.)
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
  totecosysc_1_map = where(landarea_noequil(nyrs-2,:,:) .ne. 0.,(totecosysc(nyrs-1,:,:)-totecosysc(nyrs-2,:,:))/subper,totecosysc@_FillValue)
  totecosysc_1_map!0 = "lat"
  totecosysc_1_map&lat = lat
  totecosysc_1_map!1 = "lon"
  totecosysc_1_map&lon = lon
  totecosysc_2_map = where(landarea_noequil(nyrs-3,:,:) .ne. 0.,(totecosysc(nyrs-2,:,:)-totecosysc(nyrs-3,:,:))/subper,totecosysc@_FillValue)
  copy_VarCoords(totecosysc_1_map,totecosysc_2_map)

; TOTSOMC
  totsomc_area           = totsomc*landareaC                                   ; correcting totsomc for total land area
  totsomc_pg             = totsomc_area*gtoPg                                  ; g to Pg
  totsomc_glob           = dim_sum_n(totsomc_pg, (/1,2/))                      ; sums over all latitudes
  totsomc_glob!0         = "year"                                                
  totsomc_glob&year      = yearplt                                               
  totsomc_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    totsomc_glob_del(i)  = (totsomc_glob(i+1) - totsomc_glob(i))/subper
  end do
  totsomc_glob_del!0     = "year"
  totsomc_glob_del&year  = yearpltmid
  indx = where(abs(totsomc_glob_del) .lt. glob_thresh_totsomc,1,0)
  if (all(indx .eq. 1)) then
    totsomc_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          totsomc_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      totsomc_glob_equil   = -999
    end if
  end if
  totsomc_glob_equil@_FillValue = -999
  delete(indx)

; TOTVEGC
  totvegc_area              = totvegc*landareaC                                   ; correcting totvegc for total land area
  totvegc_pg                = totvegc_area*gtoPg                                  ; g to Pg
  totvegc_glob              = dim_sum_n(totvegc_pg, (/1,2/))                      ; sums over all latitudes
  totvegc_glob!0            = "year"                                                
  totvegc_glob&year         = yearplt                                               
  totvegc_glob_del          = new(nyrs-1,"float")
  do i = 0,nyrs-2
    totvegc_glob_del(i)     = (totvegc_glob(i+1) - totvegc_glob(i))/subper
  end do
  totvegc_glob_del!0        = "year"
  totvegc_glob_del&year     = yearpltmid
  indx = where(abs(totvegc_glob_del) .lt. glob_thresh_totvegc,1,0)
  if (all(indx .eq. 1)) then
    totvegc_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          totvegc_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      totvegc_glob_equil   = -999
    end if
  end if
  totvegc_glob_equil@_FillValue = -999
  delete(indx)

; TLAI
  areasum             = sum(area*landfrac)
  areaL               = area*landfrac
  landareaL           = conform_dims(dimsizes(tlai),areaL,(/1,2/))                ; conforming dimensions of areaL to tlai
  tlai_glob           = dim_sum_n(tlai*landareaL/areasum,(/1,2/))                 ; weighted global tlai
  tlai_glob!0         = "year"                                                
  tlai_glob&year      = yearplt                                               
  tlai_glob_del       = new(nyrs-1,"float")
  do i = 0,nyrs-2
    tlai_glob_del(i) = (tlai_glob(i+1) - tlai_glob(i))/subper
  end do
  tlai_glob_del!0     = "year"
  tlai_glob_del&year  = yearpltmid
  indx = where(abs(tlai_glob_del) .lt. glob_thresh_tlai,1,0)
  if (all(indx .eq. 1)) then
    tlai_glob_equil = yearplt(0)
  else
    if (.not.(all(indx .eq. 0)) .and. indx(dimsizes(indx)-1) .ne. 0) then
      indxrev = indx(::-1)
      do i = 0,dimsizes(indxrev)-1
        if (indxrev(i) .eq. 0) then
          tlai_glob_equil = yearpltrev(i-1)
          break 
        end if
      end do
      delete(indxrev)
    else
      tlai_glob_equil   = -999
    end if
  end if
  tlai_glob_equil@_FillValue = -999
  delete(indx)

; GPP
  gpp_area              = gpp*landareaC                                   ; correcting gpp for total land area
  gpp_pg                = gpp_area*gtoPg*secinyr                          ; g to Pg and sec to yrs
  gpp_glob              = dim_sum_n(gpp_pg, (/1,2/))                      ; sums over all latitudes
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

; FERT
; landareaC           = conform_dims(dimsizes(fert),landarea,(/1,2/)) ; conforming dimensions of landarea to totecosysc
; fert_area           = fert*landareaC                                ; correcting totecosysc for total land area
; fert_pg             = fert_area*gtoPg*secinyr                       ; g to Pg and sec to yrs
; fert_glob           = dim_sum_n(fert_pg, (/1,2/))                   ; sums over all latitudes
; print(fert_glob)

;===============================Plotting====================================;
  if (do_plot .eq. "True") then

; wks_type = "png"
; wks_type@wkWidth  = 2500
; wks_type@wkHeight = 2500
; wks_type = "x11"
  wks_type = "ps"
; wks_type = "pdf"
  wks                        = gsn_open_wks (wks_type, caseid+"_Spinup")
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
  res@tiYAxisString          = "Pg C"
  res@tiMainString           = "TOTECOSYSC"
  plot(0)                    = gsn_csm_xy(wks,yearplt,totecosysc_glob,res)

  res@tiYAxisString          = "Pg C"
  res@tiMainString           = "Delta TOTECOSYSC " + "EqYr: "+totecosysc_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(1)                    = gsn_csm_xy(wks,yearpltmid,totecosysc_glob_del,res)
  prim1 = gsn_add_polyline(wks,plot(1),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim2 = gsn_add_polyline(wks,plot(1),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_totecosysc,-glob_thresh_totecosysc/),polyres2)
  prim3 = gsn_add_polyline(wks,plot(1),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_totecosysc,glob_thresh_totecosysc/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "Pg C"
  res@tiMainString           = "TOTSOMC"
  plot(2)                    = gsn_csm_xy(wks,yearplt,totsomc_glob,res)

  res@tiYAxisString          = "Pg C"
  res@tiMainString           = "Delta TOTSOMC " + "EqYr: "+totsomc_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(3)                    = gsn_csm_xy(wks,yearpltmid,totsomc_glob_del,res)
  prim4 = gsn_add_polyline(wks,plot(3),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim5 = gsn_add_polyline(wks,plot(3),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_totsomc,-glob_thresh_totsomc/),polyres2)
  prim6 = gsn_add_polyline(wks,plot(3),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_totsomc,glob_thresh_totsomc/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "Pg C"
  res@tiMainString           = "TOTVEGC"
  plot(4)                    = gsn_csm_xy(wks,yearplt,totvegc_glob,res)

  res@tiYAxisString          = "Pg C"
  res@tiMainString           = "Delta TOTVEGC " + "EqYr: "+totvegc_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(5)                    = gsn_csm_xy(wks,yearpltmid,totvegc_glob_del,res)
  prim7 = gsn_add_polyline(wks,plot(5),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim8 = gsn_add_polyline(wks,plot(5),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_totvegc,-glob_thresh_totvegc/),polyres2)
  prim9 = gsn_add_polyline(wks,plot(5),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_totvegc,glob_thresh_totvegc/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "m ~S~2~N~ m~S~-2~N~"
  res@tiMainString           = "TLAI"
  plot(6)                    = gsn_csm_xy(wks,yearplt,tlai_glob,res)

  res@tiYAxisString          = "m ~S~2~N~ m~S~-2~N~"
  res@tiMainString           = "Delta TLAI " + "EqYr: "+tlai_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(7)                    = gsn_csm_xy(wks,yearpltmid,tlai_glob_del,res)
  prim10 = gsn_add_polyline(wks,plot(7),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim11 = gsn_add_polyline(wks,plot(7),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_tlai,-glob_thresh_tlai/),polyres2)
  prim12 = gsn_add_polyline(wks,plot(7),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_tlai,glob_thresh_tlai/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "Pg C yr~S~-1~N~"
  res@tiMainString           = "GPP"
  plot(8)                    = gsn_csm_xy(wks,yearplt,gpp_glob,res)

  res@tiYAxisString          = "Pg C yr~S~-1~N~"
  res@tiXAxisString          = "Spinup Year"
  res@tiMainString           = "Delta GPP " + "EqYr: "+gpp_glob_equil
  res@trYMaxF                = 0.2
  res@trYMinF                = -0.2
  plot(9)                    = gsn_csm_xy(wks,yearpltmid,gpp_glob_del,res)
  prim13 = gsn_add_polyline(wks,plot(9),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim14 = gsn_add_polyline(wks,plot(9),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_gpp,-glob_thresh_gpp/),polyres2)
  prim15 = gsn_add_polyline(wks,plot(9),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_gpp,glob_thresh_gpp/),polyres2)

  delete(res@trYMaxF)
  delete(res@trYMinF)
  res@tiYAxisString          = "m"
  res@tiMainString           = "TWS"
  plot(10)                   = gsn_csm_xy(wks,yearplt,tws_glob,res)

  res@tiYAxisString          = "m"
  res@tiMainString           = "Delta TWS " + "EqYr: "+tws_glob_equil
  res@trYMaxF                = 0.05
  res@trYMinF                = -0.05
  plot(11)                   = gsn_csm_xy(wks,yearpltmid,tws_glob_del,res)
  prim16 = gsn_add_polyline(wks,plot(11),(/0.,yearpltmid(nyrs-2)/),(/0.,0./),polyres1)
  prim17 = gsn_add_polyline(wks,plot(11),(/0.,yearpltmid(nyrs-2)/),(/-glob_thresh_tws,-glob_thresh_tws/),polyres2)
  prim18 = gsn_add_polyline(wks,plot(11),(/0.,yearpltmid(nyrs-2)/),(/glob_thresh_tws,glob_thresh_tws/),polyres2)

  res@tiYAxisString          = "%"
  res@tiMainString           = "% Land Area in TOTECOSYSC Disequil " + "EqYr: "+perc_landarea_glob_noequil
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
  resc@gsnDraw               = True
  resc@gsnFrame              = False
  resc@lbAutoManage          = False
  resc@lbLabelFontHeightF    = 0.010
  resc@txFontHeightF         = 0.008
  resc@cnLevels              = (/-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0/)
  resc@gsnLeftString         = "gC m~S~-2~N~"
  resc@gsnRightString        = ""

  resc@vpXF                  = 0.30                  ; position and sizes
  resc@vpYF                  = 0.28                  ; for XY plot
  resc@vpWidthF              = 0.30
  resc@vpHeightF             = 0.30
  resc@gsnCenterString       = "TOTECOSYSC Disequil Yr " + yearplt(nyrs-1) + " - " + yearplt(nyrs-2)
  plot                       = gsn_csm_contour_map(wks,totecosysc_1_map,resc)

  resc@vpXF                  = 0.65                  ; position and sizes
  resc@vpYF                  = 0.28                  ; for XY plot
  resc@vpWidthF              = 0.30
  resc@vpHeightF             = 0.30
  resc@gsnCenterString       = "TOTECOSYSC Disequil Yr " + yearplt(nyrs-2) + " - " + yearplt(nyrs-3)
  plot                       = gsn_csm_contour_map(wks,totecosysc_2_map,resc)

  frame(wks)

  end if   ; end do_plot

; Equilibrium summary
  print((/"======================================================================="/))
  print((/"======================================================================="/))
  print((/"EQUILIBRIUM SUMMARY"/))
  print((/"======================================================================="/))
  if (.not.(ismissing(totecosysc_glob_equil))) then
    print((/"TOTECOSYSC is in equilibrium. Eq. Yr. = "+totecosysc_glob_equil/))
  else
    print((/"FATAL: TOTECOSYSC is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(totsomc_glob_equil))) then
    print((/"TOTSOMC is in equilibrium. Eq. Yr. = "+totsomc_glob_equil/))
  else
    print((/"FATAL: TOTSOMC is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(totvegc_glob_equil))) then
    print((/"TOTVEGC is in equilibrium. Eq. Yr. = "+totvegc_glob_equil/))
  else
    print((/"FATAL: TOTVEGC is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(tlai_glob_equil))) then
    print((/"TLAI is in equilibrium. Eq. Yr. = "+tlai_glob_equil/))
  else
    print((/"FATAL: TLAI is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(gpp_glob_equil))) then
    print((/"GPP is in equilibrium. Eq. Yr. = "+gpp_glob_equil/))
  else
    print((/"FATAL: GPP is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(tws_glob_equil))) then
    print((/"TWS is in equilibrium. Eq. Yr. = "+tws_glob_equil/))
  else
    print((/"WARNING: TWS is NOT in equilibrium"/))
  end if
  if (.not.(ismissing(perc_landarea_glob_noequil))) then
    print((/"At least "+(100.-glob_thresh_area)+" percent of the land surface is in TOTECOSYSC equilibrium. Eq. Yr. = "+perc_landarea_glob_noequil/))
  else
    print((/"FATAL: Not enough of the land surface is in equilibrium ("+sprintf("%6.2f",perc_landarea_noequil(nyrs-2))+"% > "+sprintf("%6.2f",glob_thresh_area)+"%)"/))
  end if
  if (.not.(ismissing(totecosysc_glob_equil)) .and. \
      .not.(ismissing(totsomc_glob_equil))    .and. \
      .not.(ismissing(totvegc_glob_equil))    .and. \
      .not.(ismissing(tlai_glob_equil))       .and. \
      .not.(ismissing(gpp_glob_equil))        .and. \
      .not.(ismissing(perc_landarea_glob_noequil))) then
    print((/"Congratulations! Your simulation is in equilibrium"/))
  else
    print((/"FATAL: Your simulation is not in equilibrium, 8 hours have been deducted from your PTO bank, try again"/))
  end if
  print((/"======================================================================="/))

  print ("=========================================")
  print ("Finish Time: "+systemfunc("date") )
  print ("=========================================")

end
