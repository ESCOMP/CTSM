; NCL script
; create_PLUMBER2_sites_csv.ncl ; Keith Oleson, Sep 2023
; This script generates a csv file for the NEON site_and_regional infrastructure
;**************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  print ("=========================================")
  print ("Start Time: "+systemfunc("date") )
  print ("=========================================")


; Sites
  sim = systemfunc("ls -1 " + "/glade/work/oleson/PLUMBER2/datm_files/")

  nsim = dimsizes(sim)
  site_num       = new(nsim,"integer")
  atm_ncpl       = new(nsim,"string")
  startyear      = new(nsim,"string")
  endyear        = new(nsim,"string")
  localstartyear = new(nsim,"string")
  startmonth     = new(nsim,"string")
  startday       = new(nsim,"string")
  run_startdate  = new(nsim,"string")
  starttod       = new(nsim,"string")
  nsteps         = new(nsim,"string")
  latitude       = new(nsim,"double")
  longitude      = new(nsim,"double")

  pfttype = new((/nsim,2/),"integer")
  pftperc = new((/nsim,2/),"double")
  pctwet  = new(nsim,"double")

  do sel = 0,nsim - 1

  site_num(sel) = sel+1

; Specify input directory
  f_in_dir  = "/glade/work/oleson/PLUMBER2/datm_files/"+sim(sel)+"/CLM1PT_data/"
  fls = systemfunc("ls " + f_in_dir + "CTSM_DATM_" + sim(sel) + "_*-*"+".nc")
  command = "echo " + f_in_dir + "CTSM_DATM_" + sim(sel) + "_*.nc" + " | cut -c76-79"
  startyear(sel) = systemfunc(command)
  command = "echo " + f_in_dir + "CTSM_DATM_" + sim(sel) + "_*.nc" + " | cut -c81-84"
  endyear(sel) = systemfunc(command)
  print(fls)
  datmfile = addfile(fls,"r")
  time = datmfile->time
  printVarSummary(time)
  atm_ncpl(sel) = tostring(sprinti("%0.2i ",tointeger(86400./(time(1)-time(0)))))
  print(atm_ncpl(sel))
  if (isint64(time)) then
    tmp = toint(time)
    copy_VarCoords(time,tmp)
    copy_VarAtts(time,tmp)
  else
    tmp  = time
    copy_VarCoords(time,tmp)
    copy_VarAtts(time,tmp)
  end if
  printVarSummary(tmp)
  utc_date = cd_calendar(tmp(0),0)
  print(utc_date)
  localstartyear(sel) = tostring(tointeger(utc_date(0,0)))
  print(localstartyear(sel))
  startmonth(sel) = tostring(sprinti("%0.2i ",tointeger(utc_date(0,1))))
  print(startmonth(sel))
  startday(sel) = tostring(sprinti("%0.2i ",tointeger(utc_date(0,2))))
  print(startday(sel))
  starttod(sel) = tostring(tointeger(utc_date(0,3)/24 * 86400))
  print(starttod(sel))
  nsteps(sel) = dimsizes(tmp)

  latitude(sel) = datmfile->LATIXY(0,0)
  longitude(sel) = datmfile->LONGXY(0,0)

  if (longitude(sel) .gt. 180.d) then
    longitude(sel) = longitude(sel) - 360.d
  end if


  delete(time)
  delete(tmp)
  delete(utc_date)

; run_startdate(sel) = localstartyear(sel)+"-"+startmonth(sel)+"-"+startday(sel)
  run_startdate(sel) = tostring(sprinti("%0.4i",tointeger(localstartyear(sel))))+"-"+tostring(sprinti("%0.2i",tointeger(startmonth(sel))))+"-"+tostring(sprinti("%0.2i",tointeger(startday(sel))))

  f_in_dir  = "/glade/work/oleson/PLUMBER2/input_files/"+sim(sel)+"/"
  fls = systemfunc("ls " + f_in_dir + "surfdata_0.9x1.25_16pfts_" + sim(sel) + ".nc")
  surfdatafile = addfile(fls,"r")
  pct_natveg = surfdatafile->PCT_NATVEG(0,0)
  if (pct_natveg .gt. 0.) then
    pct_nat_pft = surfdatafile->PCT_NAT_PFT(:,0,0)
    indx = ind(pct_nat_pft(:) .gt. 0)
    npfts = dimsizes(indx)
    print(npfts)
    if (npfts .gt. 2) then
      print((/"Number of pfts is > 2"/))
      pfttype(sel,:) = -999
      pftperc(sel,:) = -999.
    else
      pfttype(sel,:) = -999
      pftperc(sel,:) = -999.
      do p = 0,npfts-1
         pfttype(sel,p) = indx(p)
         pftperc(sel,p) = pct_nat_pft(indx(p))
      end do
    end if
    pctwet(sel) = 0.
    delete(pct_nat_pft)
    delete(indx)
  else
    print((/"PCT_NATVEG is zero"/))
    pfttype(sel,:) = -999
    pftperc(sel,:) = -999.
    pct_wetland = surfdatafile->PCT_WETLAND(0,0)
    if (pct_wetland .gt. 0.) then
      print((/"Wetland Site"/))
      pctwet(sel) = pct_wetland
    else
      print((/"PCT_NATVEG and PCT_WETLAND are both zero"/))
      pctwet(sel) = 0.
      pct_crop = surfdatafile->PCT_CROP(0,0)
      if (pct_crop .gt. 0.) then
        print((/"Crop Site"/))
        pct_cft = surfdatafile->PCT_CFT(:,0,0)
        if (pct_cft(0) .gt. 0.) then
          pfttype(sel,0) = 15
          pftperc(sel,0) = (/pct_cft(0)/)
        else
          pfttype(sel,0) = -999
          pftperc(sel,0) = -999.
        end if
        if (pct_cft(1) .gt. 0.) then
          pfttype(sel,1) = 16
          pftperc(sel,1) = (/pct_cft(1)/)
        else
          pfttype(sel,1) = -999
          pftperc(sel,1) = -999.
        end if
      else
        print((/"PCT_NATVEG, PCT_WETLAND, PCT_CROP are all zero"/))
      end if
    end if
  end if



  end do

  print((/"Done with Loop"/))

  fname = "/glade/work/oleson/PLUMBER2/NEON_PLUMBER/PLUMBER2_sites.csv"
  header = (/"#start_year and end_year are used to define DATM_YR_START and DATM_YR_END.", \
             "#RUN_STARTDATE and START_TOD are specified because we are starting at GMT corresponding to local midnight.", \
             "#ATM_NCPL is specified so that the time step of the model matches the time interval specified by the atm forcing data.", \
             ",Site,Lat,Lon,pft1,pft1-%,pft2,pft2-%,start_year,end_year,RUN_STARTDATE,START_TOD,ATM_NCPL"/)
  hlist = [/header/]
  alist = [/site_num, sim, latitude, longitude, pfttype(:,0), pftperc(:,0), pfttype(:,1), pftperc(:,1), startyear, endyear, run_startdate, starttod, atm_ncpl/]
  print(alist)
  print((/"Done defining hlist and alist"/))
  write_table(fname, "w", hlist, "%s")
  print((/"Done with header"/))
  write_table(fname, "a", alist, "%d,%s,%9.6f,%9.6f,%d,%5.2f,%d,%5.2f,%s,%s,%s,%s,%s")

  print ("=========================================")
  print ("Finish Time: "+systemfunc("date") )
  print ("=========================================")

  end

