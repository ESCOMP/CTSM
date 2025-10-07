; Script written by Fang Li and emailed to slevis on 2025/07/24
; Script for generating the hdm (aka popdens or population density) stream file
; /glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/firedata/clmforc.Li_2025_CMIP7_hdm_0.5x0.5_simyr1850-2100_c250717.nc
; starting from preexisting file
; /glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/firedata/clmforc.Li_2018_SSP3_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc
; and raw data in /glade/work/fangli/hd/TRENDY/pop-dens_input4MIPs_population_CMIP_PIK-CMIP-1-0-0_gr_1850-2025.nc
; Is the latter the same as /glade/campaign/cesm/cesmdata/input4MIPs_raw/input4MIPs/CMIP7/CMIP/PIK/PIK-CMIP-1-0-0 ?

f=addfile("/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/firedata/clmforc.Li_2018_SSP3_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2100_c181205.nc","r")
hdm=f->hdm

f1=addfile("/glade/work/fangli/hd/TRENDY/pop-dens_input4MIPs_population_CMIP_PIK-CMIP-1-0-0_gr_1850-2025.nc","r")
hdmn=f1->pop_dens

f2=addfile("/glade/work/fangli/hd/TRENDY/ctl05.clm2.h0.1850-01.nc","r")
landfrac=f2->landfrac
landmask=f2->landmask
landf1=landfrac*landmask
landf=landf1
landf(:,0:359)=landf1(:,360:719)
landf(:,360:719)=landf1(:,0:359)

do ilat=0,359
 do ilon=0,719
  if(ismissing(landf(ilat,ilon)))then
      hdm(0:175,ilat,ilon)=(/hdmn(:,ilat,ilon)/)
  else
     hdm(0:175,ilat,ilon)=(/hdmn(:,ilat,ilon)/landf(ilat,ilon)/)
  end if 
  end do
end do

fout=addfile("clmforc.Li_2025_CMIP7_hdm_0.5x0.5_simyr1850-2100_c250717.nc","c")
 vNam=getfilevarnames(f)
    do i=0,dimsizes(vNam)-1
    if(vNam(dimsizes(vNam)-1-i).ne."hdm")then
       x=f->$vNam(dimsizes(vNam)-1-i)$
       fout->$vNam(dimsizes(vNam)-1-i)$=x   
       delete(x) 
      else
       fout->hdm=hdm
      end if
    end do

