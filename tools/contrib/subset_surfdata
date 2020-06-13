#! /usr/bin/env python
import sys
import os
from getpass import getuser
import string
import subprocess
import numpy as np
import xarray as xr

def mprint(mstr):
    vnum=sys.version_info[0]
    if vnum == 3:
        print(mstr)
    if vnum == 2:
        print mstr
        
myname=getuser()
pwd=os.getcwd()
mprint(myname)
mprint(pwd)

#creates regional surface dataset and domain file

#--  Specify input and output directories
dir_output='/glade/scratch/'+myname+'/regional/'

#--  Create regional CLM domain file
create_domain   = True
#--  Create CLM surface data file
create_surfdata = True
#--  Create CLM surface data file
create_landuse  = False

tagnum=1
if tagnum == 1:
    tag='S.America'

    ln1=275.
    ln2=330.
    lt1=-40.
    lt2=15.

if tagnum == 2:
    tag='Western.US'

    ln1=284.
    ln2=296.
    lt1=44.
    lt2=53.

#--  Set time stamp
command='date "+%y%m%d"'
x2=subprocess.Popen(command,stdout=subprocess.PIPE,shell='True')
x=x2.communicate()
timetag = x[0].strip()

#--  Specify land domain file  ---------------------------------
fdomain  = '/glade/p/cesmdata/cseg/inputdata/share/domains/domain.lnd.fv1.9x2.5_gx1v7.170518.nc'
#fdomain2 = dir_output + 'domain.lnd.fv0.9x1.25_gx1v6.'+tag+'.090309.nc'
fdomain2 = dir_output + 'domain.lnd.fv1.9x2.5_gx1v7.'+tag+'_170518.nc'

#--  Specify surface data file  --------------------------------
fsurf    = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/surfdata_1.9x2.5_78pfts_CMIP6_simyr1850_c170824.nc'
#fsurf2   = dir_output + 'surfdata_0.9x1.25_16pfts_CMIP6_simyr2000_'+tag+'.c170706.nc'
fsurf2   = dir_output + 'surfdata_1.9x2.5_78pfts_CMIP6_simyr1850_'+tag+'_c170824.nc'

#--  Specify landuse file  -------------------------------------
fluse    = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/landuse.timeseries_1.9x2.5_hist_78pfts_CMIP6_simyr1850-2015_c170824.nc'
fluse2   = dir_output + 'landuse.timeseries_1.9x2.5_hist_78pfts_CMIP6_simyr1850-2015_'+tag+'.c170824.nc'

#--  Create CTSM domain file
if create_domain:
    f1  = xr.open_dataset(fdomain)
    # create 1d coordinate variables to enable sel() method
    lon0=np.asarray(f1['xc'][0,:])
    lat0=np.asarray(f1['yc'][:,0])
    lon=xr.DataArray(lon0,name='lon',dims='ni',coords={'ni':lon0})
    lat=xr.DataArray(lat0,name='lat',dims='nj',coords={'nj':lat0})
    # assign() not working on cheyenne
    #f2=f1.assign({'lon':lon,'lat':lat})
    f2=f1.assign()
    f2['lon'] = lon
    f2['lat'] = lat
    f2.reset_coords(['xc','yc'],inplace=True)

    # subset longitude and latitude arrays
    xind=np.where((lon >= ln1) & (lon <= ln2))[0]
    yind=np.where((lat >= lt1) & (lat <= lt2))[0]
    f3=f2.isel(nj=yind,ni=xind)

    wfile=fdomain2
    # mode 'w' overwrites file
    f3.to_netcdf(path=wfile, mode='w')
    mprint('created file '+fdomain2)
    f1.close(); f2.close(); f3.close()

#--  Create CTSM surface data file
if create_surfdata:
    f1  = xr.open_dataset(fsurf)
    # create 1d variables
    lon0=np.asarray(f1['LONGXY'][0,:])
    lon=xr.DataArray(lon0,name='lon',dims='lsmlon',coords={'lsmlon':lon0})
    lat0=np.asarray(f1['LATIXY'][:,0])
    lat=xr.DataArray(lat0,name='lat',dims='lsmlat',coords={'lsmlat':lat0})
    #f2=f1.assign({'lon':lon,'lat':lat})
    f2=f1.assign()
    f2['lon'] = lon
    f2['lat'] = lat
    # subset longitude and latitude arrays
    xind=np.where((lon >= ln1) & (lon <= ln2))[0]
    yind=np.where((lat >= lt1) & (lat <= lt2))[0]
    f3=f2.isel(lsmlat=yind,lsmlon=xind)

    # mode 'w' overwrites file
    f3.to_netcdf(path=fsurf2, mode='w')
    mprint('created file '+fsurf2)
    f1.close(); f2.close(); f3.close()

#--  Create CTSM transient landuse data file
if create_landuse:
    f1  = xr.open_dataset(fluse)
    # create 1d variables
    lon0=np.asarray(f1['LONGXY'][0,:])
    lon=xr.DataArray(lon0,name='lon',dims='lsmlon',coords={'lsmlon':lon0})
    lat0=np.asarray(f1['LATIXY'][:,0])
    lat=xr.DataArray(lat0,name='lat',dims='lsmlat',coords={'lsmlat':lat0})
    #f2=f1.assign({'lon':lon,'lat':lat})
    f2=f1.assign()
    f2['lon'] = lon
    f2['lat'] = lat
    # subset longitude and latitude arrays
    xind=np.where((lon >= ln1) & (lon <= ln2))[0]
    yind=np.where((lat >= lt1) & (lat <= lt2))[0]
    f3=f2.isel(lsmlat=yind,lsmlon=xind)
    # mode 'w' overwrites file
    f3.to_netcdf(path=fluse2, mode='w')
    mprint('created file '+fluse2)
    f1.close(); f2.close(); f3.close()



