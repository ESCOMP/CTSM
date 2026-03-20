; Script written by Fang Li and emailed to slevis on 2025/07/24
; Script for generating the abm (peak crop-fire month) raw dataset
; /glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/rawdata/mksrf_abm_0.5x0.5_simyr2000.c250715.nc
; starting from preexisting raw dataset
; /glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/rawdata/mksrf_abm_0.5x0.5_simyr2000.c240821.nc

f1=addfile("/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/rawdata/mksrf_abm_0.5x0.5_simyr2000.c240821.nc","r")
abm1=f1->abm
lon1=f1->lon
lat1=f1->lat

fils=systemfunc("ls /glade/work/fangli/obs/GFED5/crop05/BA*.nc")
fs=addfiles(fils,"r")
ListSetType(fs,"join")
baf=fs[:]->baf
bafm=dim_avg_n(baf,0)

f=addfile("fpc_crop05.nc","r") ; MCD12C1 12+14
fpc_crop=f->fpc_crop

fw=addfile("fpc_water05.nc","r")
fpc_water=fw->fpc_water

abm=abm1
do ilat=0, 359
 do ilon=0, 719
   if(fpc_crop(ilat,ilon).gt.0.005 .and. fpc_water(ilat,ilon).lt.0.5)then 
      if(sum(bafm(:,ilat,ilon)).gt.0.0)then
          abm(ilat,ilon)=maxind(bafm(:,ilat,ilon))+1
         else
            abm(ilat,ilon)=13
         end if
   else
      abm(ilat,ilon)=14
  end if
end do
end do

 fout1=addfile("abm05-raw.nc","c")
   fout1->abm=abm

