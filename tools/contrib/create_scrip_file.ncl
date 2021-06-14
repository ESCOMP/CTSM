load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin

nums = (/"00","01","02","03","04","05","06","07","08","09", \
         "10","11","12","13","14","15","16","17","18","19", \
         "20","21","22","23","24","25","26","27","28","29", \
         "30","31","32","33","34","35","36","37","38","39"/)

;************************************************
; get WRF fields
;************************************************

;change this to point to the relevant geo_em file. 
wrf_file  = addfile("geo_em.d01.nc", "r")


lat2d = wrf_file->XLAT_M(0,:,:)
lon2d = wrf_file->XLONG_M(0,:,:)
lon2d = where(lon2d.lt.0, lon2d+360.0, lon2d)
landmask = wrf_file->LANDMASK(0,:,:)

ns_corner_lat = wrf_file->XLAT_V(0,:,:)
ns_corner_lon = wrf_file->XLONG_U(0,:,:)
ew_corner_lat = wrf_file->XLAT_U(0,:,:)
ew_corner_lon = wrf_file->XLONG_U(0,:,:)

print(lat2d(0,0))
print(lon2d(0,0))

print(ns_corner_lat(0:1,0:1))
print(ns_corner_lon(0:1,0:1))
print(ew_corner_lat(0:1,0:1))
print(ew_corner_lon(0:1,0:1))

opt = True
opt@GridMask       = landmask
opt@ForceOverwrite = True

curvilinear_to_SCRIP("wrf2clm_land.nc",lat2d,lon2d,opt)

opt@GridMask       = where(landmask.eq.0,1,0)
curvilinear_to_SCRIP("wrf2clm_ocean.nc",lat2d,lon2d,opt)

infile = addfile("wrf2clm_land.nc","r")

lat1d = infile->grid_center_lat
lon1d = infile->grid_center_lon

lat1dc = infile->grid_corner_lat
lon1dc = infile->grid_corner_lon

print(lat1d(0))
print(lon1d(0))
print(lat1dc(0,:))
print(lon1dc(0,:))

end
