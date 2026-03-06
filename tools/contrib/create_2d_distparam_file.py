#! /usr/bin/env python
import subprocess
import numpy as np
import netCDF4 as netcdf4

# creates surface dataset with distributed parameters
dir_surfdata = '/fs/cgd/csm/inputdata/lnd/clm2/surfdata_esmf/'
dir_output='./'

# for coordinates, land mask, etc.
fsurf = f'{dir_surfdata}ctsm5.3.0/surfdata_0.9x1.25_hist_2000_78pfts_c240908.nc'

outfile = f'{dir_output}paramdata_0.9x1.25_ctsm60_params.c260204.nc'
paramfile = '/fs/cgd/csm/inputdata/lnd/clm2/paramdata/ctsm60_params.c260204.no_nan_fill.nc'

# calibrated parameters from Tang2025
param_list = ['e_ice','n_baseflow','baseflow_scalar','fff','d_max','FMAX','om_frac_sf','interception_fraction','watsat_sf','sucsat_sf','bsw_sf','zbedrock_sf','hksat_sf','slopebeta','liq_canopy_storage_scalar','maximum_leaf_wetted_fraction','n_melt_coef','accum_factor','upplim_destruct_metamorph','zsno','precip_repartition_nonglc_all_rain_t','medlynslope','jmaxb0','kmax','cv','leafcn']

# snow parameters
param_list += ['xdrdt','scvng_fct_mlt_sf','snw_rds_refrz','fresh_snw_rds_max']

# remove non-standard parameters and those already on surface data file
exclude_list = ['FMAX','zbedrock_sf','om_frac_sf']
for e in exclude_list:
    param_list.remove(e)

# namelist parameters (add all_snow_t)
pnamelist = ['baseflow_scalar','precip_repartition_nonglc_all_rain_t','precip_repartition_nonglc_all_snow_t']
param_list += ['precip_repartition_nonglc_all_snow_t']
# assign values from namelist_defaults_ctsm.xml
pnamelist_values = [1e-3,2,0]
pnamelist_attributes = [{'long_name':'baseflow scaling coefficient','units':'unitless'},{'long_name':'all rain temperature','units':'degrees C'},{'long_name':'all snow temperature','units':'degrees C'}]
    
nparams = len(param_list)

#--  Check whether file exists  ---------------------------------
command=['ls',outfile]
file_exists=subprocess.call(command,stderr=subprocess.PIPE)
if file_exists > 0:
    print('creating new file: ', outfile)
else:
    print('overwriting file: ', outfile)

#--  Open surface data file
f =  netcdf4.Dataset(fsurf, 'r')
global_attributes  = f.ncattrs()
variables = f.variables
dimensions = f.dimensions
jm, im = len(dimensions['lsmlat']),len(dimensions['lsmlon'])
lon = np.asarray(f.variables['LONGXY'][0,])
lat = np.asarray(f.variables['LATIXY'][:,0])

if 'PFTDATA_MASK' in f.variables.keys():            
    landmask = np.asarray(f.variables['PFTDATA_MASK'][:,])
if 'LANDFRAC_PFT' in f.variables.keys(): # ctsm >= 5.3.0
    landmask = np.where(f.variables['LANDFRAC_PFT'][:,]>0,1,0)

#--  Open parameter file
p1 =  netcdf4.Dataset(paramfile, 'r')

#--  Open output file
w =  netcdf4.Dataset(outfile, 'w', format='NETCDF3_64BIT')
w.title = 'CTSM parameters'
w.createDimension('lat',int(jm))
w.createDimension('lon',int(im))
w.createDimension('time',None)
w.createVariable('lon',float,('lon',))
w.variables['lon'].units = 'degrees east'
w.variables['lon'].long_name = 'longitude of grid cell center'
w.variables['lon'][:,] = lon
w.createVariable('lat',float,('lat',))
w.variables['lat'].units = 'degrees north'
w.variables['lat'].long_name = 'latitude of grid cell center'
w.variables['lat'][:,] = lat
w.createVariable('time',float,('time',))
w.variables['time'].units = 'days since 2000-01-01 00:00'
w.variables['time'].long_name = 'time'
w.variables['time'].calendar = 'noleap'
w.variables['time'][0,] = 0.5
for pname in param_list:
    print(pname)
    # set namelist variable values using specified values
    if pname in pnamelist:
        ip = pnamelist.index(pname)
        pvar = pnamelist_values[ip]
        #--  Create variable and populate array
        w.createVariable(pname,float,('time','lat','lon'))
        wvar = w.variables[pname][:,]
        #wvar[landmask>0,] = pvar # seems that valid data needed everywhere
        wvar[0,] = pvar

        # create spatial variation
        if 1==2 and pname == 'baseflow_scalar':
            wvar[0,:,lon < 60] = 1e-1
            wvar[0,:,lon > 250] = 1e-4
        if 1==2 and pname == 'precip_repartition_nonglc_all_rain_t':
            wvar[0,:,lon < 90] = 1
            wvar[0,:,lon > 250] = 3
        if 1==2 and pname == 'precip_repartition_nonglc_all_snow_t':
            wvar[0,:,lon < 90] = -3
            wvar[0,:,lon > 250] = -4
            
            
        w.variables[pname][0,] = wvar
        #--  Set attribute values
        for att in pnamelist_attributes[ip].keys():
            #print(pname,att,pnamelist_attributes[ip][att])
            w.variables[pname].setncattr(att,pnamelist_attributes[ip][att])
    else: # use parameter values from parameter file
        pvar = p1.variables[pname][:,]
        # only create 2d fields from scalars for now
        if len(pvar.shape) == 0:
            #--  Create variable and populate array

            #(shr_strdata_readstrm) ERROR: only double, real and short types are supported for stream read
            print('dtype ',pvar.dtype)
            if pvar.dtype in ['int','int32','int64']:
                #pvar = pvar.astype('short')
                pvar = pvar.astype(float)
            w.createVariable(pname,pvar.dtype,('time','lat','lon'))
            #w.variables[pname][landmask>0,] = pvar
            wvar = w.variables[pname][:,]
            #wvar[landmask>0,] = pvar # seems that valid data needed everywhere
            wvar[0,] = pvar

            if 1==2 and pname == 'd_max':
                wvar[0,:,lon < 60] = 5
                wvar[0,:,lon > 250] = 30
            if 1==2 and pname == 'xdrdt':
                wvar[0,lat < -60] = 1

            w.variables[pname][0,] = wvar
            #--  Set attribute values
            att = p1.variables[pname].ncattrs()
            for attname in att:
                w.variables[pname].setncattr(attname,p1.variables[pname].getncattr(attname))

w.close()
p1.close()
f.close()
print('created ',outfile)
