#! /usr/bin/env python
import sys
import string
import subprocess
import numpy as np
from scipy.io import netcdf
import netCDF4 as netcdf4

'''
# load proper modules first, i.e.
 module load python/2.7.7
 module load all-python-libs
'''
#creates regional surface dataset and domain file

dir_input='/glade/p/cesm/cseg/inputdata/lnd/clm2/'
# must change this to a directory to which you have permissions
dir_output='/glade/scratch/'

#--  create regional CLM domain file
create_domain   = 0
#--  create ergional CLM surface data file
create_surfdata = 1
#--  create river direction file
create_rdirc = 1

tagnum=1
if tagnum == 1:
    tag='NE.US'
    rnum=2
    resnum=2
if tagnum == 2:
    tag='Western.US'
    rnum=2
    resnum=2
if tagnum == 3:
    tag='S.America'
    rnum=2
    resnum=1
if tagnum == 4:
    tag='California'
    rnum=2
    resnum=2

if resnum == 1:
    fdomain = '/glade/p/cesmdata/cseg/inputdata/share/domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc'

    fsurf = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/surfdata_0.9x1.25_simyr2000_c141219.nc'

    rdirc='/glade/p/cesmdata/cseg/inputdata/lnd/clm2/rtmdata/rdirc_0.5x0.5_simyr2000_slpmxvl_c120717.nc'

    fdomain2=dir_output+'domain.lnd.fv0.9x1.25_gx1v6.'+tag+'.nc'

    fsurf2=dir_output+'surfdata_0.9x1.25_16pfts_simyr2000_UHDMedian_'+tag+'.nc'
    rdirc2=dir_output+'rdirc_0.5x0.5_simyr2000_slpmxvl_'+tag+'.nc'

if resnum == 2:
    fdomain = '/glade/p/cesmdata/cseg/inputdata/share/domains/domain.clm/domain.lnd.0.125x0.125_tx0.1v2.140704.nc'
    fsurf = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/surfdata_0.125x0.125_simyr2000_c150114.nc'
    rdirc='/glade/p/cesmdata/cseg/inputdata/lnd/clm2/rtmdata/rdirc_0.1x0.1_simyr2000_c110712.nc'

    fdomain2=dir_output+'domain.lnd.0.125x0.125_tx0.1v2.'+tag+'.nc'
    fsurf2=dir_output+'surfdata_0.125x0.125_simyr2000.'+tag+'.nc'
    rdirc2=dir_output+'rdirc_0.1x0.1_'+tag+'.nc'



#--  read in coordinates from original domain file
f1 =  netcdf4.Dataset(fdomain, 'r', format='NETCDF4')
xc  = f1.variables['xc']
yc  = f1.variables['yc']
mask  = np.copy(f1.variables['mask'])

#--  convert coordinates to 1d
lon=np.asarray(xc[0,:])
lat=np.asarray(yc[:,0])
im=lon.size
jm=lat.size

if rnum == 2:
    if tagnum == 1: # NE.US
        ln1=284.
        ln2=296.
        lt1=44.
        lt2=53.
    if tagnum == 2: # W.US
        ln1=235.
        ln2=270.
        lt1=27.
        lt2=50.

    if tagnum == 3: # S.America
        ln1=275.
        ln2=330.
        lt1=-40.
        lt2=15.
    if tagnum == 4: # California
        ln1=235
        ln2=246.5
        lt1=31.5
        lt2=45.

    #--  create regional mask/area arrays
    #mask=np.zeros((jm,im))
    #where returns a tuple, extract list w/ '[0]'
    xind=np.where((lon >= ln1) & (lon <= ln2))[0]
    yind=np.where((lat >= lt1) & (lat <= lt2))[0]
    im_new=xind.size
    jm_new=yind.size

    mask=mask[yind,:]; mask=mask[:,xind]

    lonc=np.copy(xc) ; latc=np.copy(yc)
    lonc=lonc[yind,:]; lonc=lonc[:,xind]
    latc=latc[yind,:]; latc=latc[:,xind]

    #print lonc

dlon = np.abs(lon[0] - lon[1])
dlat = np.abs(lat[0] - lat[1])

#--  create domain file
##1
if create_domain == 1:
    frac  = np.copy(f1.variables['frac'])
    area  = np.copy(f1.variables['area'])

    frac=frac[yind,:]; frac=frac[:,xind]
    area=area[yind,:]; area=area[:,xind]


    nv = 4
#--  Create vertex matrices  ----------------------------
    lonv=np.empty((jm_new,im_new,nv), dtype='float64')
    latv=np.empty((jm_new,im_new,nv), dtype='float64')
#--  order is SW, SE, NE, NW
    lonv[:,:,0] = lonc - 0.5*dlon
    latv[:,:,0] = latc - 0.5*dlat
    lonv[:,:,1] = lonc + 0.5*dlon
    latv[:,:,1] = latc - 0.5*dlat
    lonv[:,:,2] = lonc + 0.5*dlon
    latv[:,:,2] = latc + 0.5*dlat
    lonv[:,:,3] = lonc - 0.5*dlon
    latv[:,:,3] = latc + 0.5*dlat

    #--  Check whether file exists  ---------------------------------
    command=['ls',fdomain2]
    file_exists=subprocess.call(command,stderr=subprocess.PIPE)
    
    w = netcdf4.Dataset(fdomain2, 'w', format='NETCDF4')
    
    if file_exists > 0:
        print 'creating new file: ', fdomain2
    else:
        print 'overwriting file: ', fdomain2
    command='date "+%y%m%d"'
    x2=subprocess.Popen(command,stdout=subprocess.PIPE,shell='True')
    x=x2.communicate()
    timetag = x[0].strip()
    w.creation_date = timetag 
    w.source_file = fdomain
    w.title = 'CESM domain data'
    w.createDimension('nv',int(nv))
    w.createDimension('ni',int(im_new))
    w.createDimension('nj',int(jm_new))
    
    wyc   = w.createVariable('yc','f8',('nj','ni'))
    wxc   = w.createVariable('xc','f8',('nj','ni'))
    wyv   = w.createVariable('yv','f8',('nj','ni','nv'))
    wxv   = w.createVariable('xv','f8',('nj','ni','nv'))
    wmask = w.createVariable('mask','i4',('nj','ni'))
    warea = w.createVariable('area','f8',('nj','ni'))
    wfrac = w.createVariable('frac','f8',('nj','ni'))
    
    wyc.units  = 'degrees north'
    wxc.units  = 'degrees east'
    wyv.units  = 'degrees north'
    wxv.units  = 'degrees east'
    wmask.units  = 'unitless'
    warea.units  = 'radians squared'
    wfrac.units  = 'unitless'
    
    wyc.long_name  = 'latitude of grid cell center'
    wxc.long_name  = 'longitude of grid cell center'
    wyv.long_name  = 'latitude of grid cell vertices'
    wxv.long_name  = 'longitude of grid cell vertices'
    wmask.long_name  = 'land domain mask'
    warea.long_name  = 'area of grid cell in radians squared'
    wfrac.long_name  = 'fraction of grid cell that is active'

    # write to file  --------------------------------------------
    wyc[:,:] = latc
    wxc[:,:] = lonc
    wyv[:,:,:] = latv
    wxv[:,:,:] = lonv
    wmask[:,:] = mask
    warea[:,:] = area
    wfrac[:,:] = frac
            
    w.close()

#--  create surface data file  --------------------------------------
##2
if create_surfdata == 1:
    f2 =  netcdf4.Dataset(fsurf, 'r', format='NETCDF4')
    global_attributes  = f2.ncattrs()
    variables = f2.variables
    dimensions = f2.dimensions

    #--  Check whether file exists  ---------------------------------
    command=['ls',fsurf2]
    file_exists=subprocess.call(command,stderr=subprocess.PIPE)
    if file_exists > 0:
        print 'creating new file: ', fsurf2
    else:
        print 'overwriting file: ', fsurf2

    #--  Open output file
    w = netcdf4.Dataset(fsurf2, 'w', format='NETCDF4')

    #--  Set global attributes
    for ga in global_attributes:
        setattr(w,ga,f2.getncattr(ga))
    #--  Set dimensions of output file
    for dim in dimensions.keys():
        print dim
        if dim == 'lsmlon':
            w.createDimension(dim,int(im_new))
        elif dim == 'lsmlat':
            w.createDimension(dim,int(jm_new))
        else:
            w.createDimension(dim,len(dimensions[dim]))

    for var in variables.keys():
        y=f2.variables[var].dimensions
        x2 = [x.encode('ascii') for x in y]
        vtype = f2.variables[var].datatype
        print var, vtype, x2
        wvar = w.createVariable(var, vtype, x2)

        if len(x2) > 0:
            fvar=np.copy(f2.variables[var])

        #--  Subset input variables
        for n in range(len(x2)):
            fdim = x2[n]
            if fdim == 'lsmlon':
                if n == 0:
                    fvar = fvar[xind,]
                if n == 1:
                    fvar = fvar[:,xind,]
                if n == 2:
                    fvar = fvar[:,:,xind,]
                if n == 3:
                    fvar = fvar[:,:,:,xind,]
            if fdim == 'lsmlat':
                if n == 0:
                    fvar = fvar[yind,]
                if n == 1:
                    fvar = fvar[:,yind,]
                if n == 2:
                    fvar = fvar[:,:,yind,]
                if n == 3:
                    fvar = fvar[:,:,:,yind,]

        #--  Set attribute values
        att=f2.variables[var].ncattrs()
        print att, '\n'
        km=len(att)
        for attname in att:
            print 'name: ',attname,' value: ',f2.variables[var].getncattr(attname)
            w.variables[var].setncattr(attname,f2.variables[var].getncattr(attname))

        #--  Write variable data to output file
        if len(x2) == 0:
            wvar[:] = np.copy(f2.variables[var][:])
        if len(wvar.shape) == 1:
            wvar[:] = fvar
        if len(wvar.shape) == 2:
            wvar[:,:] = fvar
        if len(wvar.shape) == 3:
            wvar[:,:,:] = fvar
        if len(wvar.shape) == 4:
            wvar[:,:,:,:] = fvar

   #--  Close output file
    w.close

#--  create river direction file
##3
if create_rdirc == 1:
    f2 =  netcdf4.Dataset(rdirc, 'r', format='NETCDF4')
    global_attributes  = f2.ncattrs()
    variables = f2.variables
    dimensions = f2.dimensions

    #--  create regional indices
    xc2  = np.copy(f2.variables['xc'])
    yc2  = np.copy(f2.variables['yc'])
    #--  convert coordinates to 1d
    rlon=np.asarray(xc2[0,:])
    rlon[rlon < 0] += 360.0
    rlat=np.asarray(yc2[:,0])
    rim=rlon.size
    rjm=rlat.size
    #where returns a tuple, extract list w/ '[0]'
    xind2=np.where((rlon >= ln1) & (rlon <= ln2))[0]
    yind2=np.where((rlat >= lt1) & (rlat <= lt2))[0]
    ni_new=xind2.size
    nj_new=yind2.size

    #--  Check whether file exists  ---------------------------------
    command=['ls',rdirc2]
    file_exists=subprocess.call(command,stderr=subprocess.PIPE)
    if file_exists > 0:
        print 'creating new file: ', rdirc2
    else:
        print 'overwriting file: ', rdirc2

    #--  Open output file
    w = netcdf4.Dataset(rdirc2, 'w', format='NETCDF4')

    #--  Set global attributes
    for ga in global_attributes:
        setattr(w,ga,f2.getncattr(ga))
    #--  Set dimensions of output file
    for dim in dimensions.keys():
        print dim
        if dim == 'ni':
            w.createDimension(dim,int(ni_new))
        elif dim == 'nj':
            w.createDimension(dim,int(nj_new))
        else:
            w.createDimension(dim,len(dimensions[dim]))

    for var in variables.keys():
        y=f2.variables[var].dimensions
        x2 = [x.encode('ascii') for x in y]
        vtype = f2.variables[var].datatype
        print var, vtype, x2
        wvar = w.createVariable(var, vtype, x2)

        if len(x2) > 1:
            fvar=np.copy(f2.variables[var])

        #--  Subset input variables
        for n in range(len(x2)):
            fdim = x2[n]
            if fdim == 'ni':
                if n == 0:
                    fvar = fvar[xind2,]
                if n == 1:
                    fvar = fvar[:,xind2,]
                if n == 2:
                    fvar = fvar[:,:,xind2,]
                if n == 3:
                    fvar = fvar[:,:,:,xind2,]
            if fdim == 'nj':
                if n == 0:
                    fvar = fvar[yind2,]
                if n == 1:
                    fvar = fvar[:,yind2,]
                if n == 2:
                    fvar = fvar[:,:,yind2,]
                if n == 3:
                    fvar = fvar[:,:,:,yind2,]

        #--  Set attribute values
        att=f2.variables[var].ncattrs()
        print att, '\n'
        km=len(att)
        for attname in att:
            print 'name: ',attname,' value: ',f2.variables[var].getncattr(attname)
            w.variables[var].setncattr(attname,f2.variables[var].getncattr(attname))

        #--  Set edges to zero
        if var == 'RTM_FLOW_DIRECTION':
            rjm2=fvar.shape[0]
            rim2=fvar.shape[1]
            fvar[:,0]    = 0
            fvar[:,rim2-1] = 0
            fvar[0,:]    = 0
            fvar[rjm2-1.:] = 0

        #--  Write variable data to output file
        if len(wvar.shape) == 1:
            wvar[:] = f2.variables[var][:]
        if len(wvar.shape) == 2:
            wvar[:,:] = fvar
        if len(wvar.shape) == 3:
            wvar[:,:,:] = fvar
        if len(wvar.shape) == 4:
            wvar[:,:,:,:] = fvar

   #--  Close output file
    w.close

#--  Close output file
f1.close
