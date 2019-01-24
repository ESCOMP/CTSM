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
# creates single point surface datasets, atmospheric forcing data, 
# and domain file by extracting data from global datasets.

'''
#------------------------------------------------------------------#
#---------------------  Instructions  -----------------------------#
#------------------------------------------------------------------#
After creating a case using a global compset, run preview_namelist.  
From the resulting lnd_in file in the run directory, find the name 
of the domain file, and the surface data file.  
From the datm streams files (e.g. datm.streams.txt.CLMGSWP3v1.Precip)
find the name of the datm forcing data domain file and forcing files.  
Use these file names as the sources for the single point files to 
be created (see below).

After running this script, point to the new CLM domain and surface 
dataset using the user_nl_clm file in the case directory.  In addition, 
copy the datm.streams files to the case directory, with the prefix 
'user_', e.g. user_datm.streams.txt.CLMGSWP3v1.Precip.  Change the 
information in the user_datm.streams* files to point to the single 
point datm data (domain and forcing files) created using this script.  

The domain file is not set via user_nl_clm, but requires changing 
LND_DOMAIN and ATM_DOMAIN (and their paths) in env_run.xml.  

Using single point forcing data requires specifying the nearest 
neighbor mapping algorithm for the datm streams (usually they are 
the first three in the list) in user_nl_datm: mapalgo = 'nn','nn','nn', 
..., where the '...' can still be 'bilinear', etc, depending on the 
other streams that are being used, e.g. aerosols, anomaly forcing, 
bias correction.

Single point simulations require a single processor.  The various 
"NTASKS_*" xml variables in env_mach_pes.xml should therefore be 
set to 1 and "ROOTPE_" should be set to zero.  The mpi-serial 
libraries should also be used, and can be set in env_build.xml 
by changing "MPILIB" to "mpi-serial" prior to setting up the case.  
PIO can be changed in env_run.xml to serial as well by changing 
PIO_TYPENAME to "netcdf".

River routing can be turned off either in env_run.xml or in user_nl_mosart 
by setting "MOSART_MODE" to "NULL".  

No initial conditions are available, so set use_init_interp = .true. 
in user_nl_clm, or CLM_FORCE_COLDSTART to "on" in env_run.xml.

Single point runs should be run on caldera, in that case change the 
JOB_QUEUE to caldera from regular in env_batch.xml

'''

# set 'create_' variables to 1 to create desired file type
#--  create single point CLM domain file
create_domain   = 0
#--  create single point CLM surface data file
create_surfdata = 0
#--  create single point DATM atmospheric forcing data
create_datm = 0
#--  create single point dynamic landuse file
create_landuse = 0

# change plon/plat to values of point to be extracted from global dataset
tagnum=1
if tagnum == 1: #Single point
   plon =270.0
   plat = 35.0
   rnum = 1
   tag=str(plon)+'_'+str(plat)

# First specify input data locations, then specify output data locations        
# Directory names should end in '/'

# specify CLM domain and surface dataset; this information can be 
# found in the lnd_in file

fdomain = '/glade/p/cesmdata/cseg/inputdata/share/domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc'
fsurf = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/surfdata_0.9x1.25_78pfts_CMIP6_simyr2000_c170824.nc'
flanduse = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/landuse.timeseries_0.9x1.25_rcp8.5_simyr1850-2100_c141219.nc'

# specify DATM domain and forcing dataset; this information can be
# found in the datm.streams.txt* files

dir_input_datm='/glade/p/cesm/lmwg/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/'
fdatmdomain = '/glade/p/cesmdata/cseg/lmwg/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/domain.lnd.360x720_gswp3.0v1.c170606.nc'

solrdir  = 'Solar/'
precdir  = 'Precip/'
tpqwldir = 'TPHWL/'

solrtag  = 'clmforc.GSWP3.c2011.0.5x0.5.Solr.'
prectag  = 'clmforc.GSWP3.c2011.0.5x0.5.Prec.'
tpqwtag  = 'clmforc.GSWP3.c2011.0.5x0.5.TPQWL.'

# specify locations for output; must change this to a directory to 
# which you have permissions

dir_output     ='/glade/scratch/'
dir_output_datm='/glade/scratch/'

# specify new name for CLM domain file and surface data file
fdomain2  = dir_output+'domain.lnd.fv0.9x1.25_'+tag+'.nc'
fsurf2    = dir_output+'surfdata_0.9x1.25_'+tag+'.nc'
flanduse2 = dir_output+'landuse.timeseries_0.9x1.25_rcp8.5_simyr1850-2100_'+tag+'.nc'

# specify new name for DATM domain file
fdatmdomain2  = dir_output_datm+'domain.lnd.360x720_'+tag+'.nc'

# specify years of atm data to extract
syr=1991
eyr=1991

''' ------------------------------------------------------

End of user specification section

 ------------------------------------------------------ '''


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

dlon = np.abs(lon[0] - lon[1])
dlat = np.abs(lat[0] - lat[1])

rnum = 1
if rnum == 1:
   
    #--  create single point mask/area arrays
    xind = np.argmin(abs(lon - plon))
    yind = np.argmin(abs(lat - plat))
    im_new=xind.size
    jm_new=yind.size

    mask=mask[yind,:]; mask=mask[xind]

    lonc=np.copy(xc) ; latc=np.copy(yc)
    lonc=lonc[yind,:]; lonc=lonc[xind]
    latc=latc[yind,:]; latc=latc[xind]

#--  create CLM domain file
##1
if create_domain == 1:
    frac  = np.copy(f1.variables['frac'])
    area  = np.copy(f1.variables['area'])

    frac=frac[yind,:]; frac=frac[xind]     
    area=area[yind,:]; area=area[xind]


    nv = 4
    #--  Create vertex matrices  ----------------------------
    lonv=np.empty((nv), dtype='float64')
    latv=np.empty((nv), dtype='float64')
    #--  order is SW, SE, NE, NW
    lonv[0] = lonc - 0.5*dlon
    latv[0] = latc - 0.5*dlat
    lonv[1] = lonc + 0.5*dlon
    latv[1] = latc - 0.5*dlat
    lonv[2] = lonc + 0.5*dlon
    latv[2] = latc + 0.5*dlat
    lonv[3] = lonc - 0.5*dlon
    latv[3] = latc + 0.5*dlat

    #--  Check whether file exists  ---------------------------------
    command=['ls',fdomain2]
    file_exists=subprocess.call(command,stderr=subprocess.PIPE)
    print file_exists

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

#--  create CLM surface data file  ----------------------------------
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
                    shp=np.asarray(fvar.shape)
                    shp[0] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[0,] = fvar[xind,]
                    fvar = tmp
                if n == 1:
                    shp=np.asarray(fvar.shape)
                    shp[1] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,0,] = fvar[:,xind,]
                    fvar = tmp
                if n == 2:
                    shp=np.asarray(fvar.shape)
                    shp[2] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,:,0,] = fvar[:,:,xind,]
                    fvar = tmp
                if n == 3:
                    shp=np.asarray(fvar.shape)
                    shp[3] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,:,:,0,] = fvar[:,:,:,xind,]
                    fvar = tmp
            if fdim == 'lsmlat':
                if n == 0:
                    shp=np.asarray(fvar.shape)
                    shp[0] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[0,] = fvar[yind,]
                    fvar = tmp
                if n == 1:
                    shp=np.asarray(fvar.shape)
                    shp[1] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,0,] = fvar[:,yind,]
                    fvar = tmp
                if n == 2:
                    shp=np.asarray(fvar.shape)
                    shp[2] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,:,0,] = fvar[:,:,yind,]
                    fvar = tmp
                if n == 3:
                    shp=np.asarray(fvar.shape)
                    shp[3] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,:,:,0,] = fvar[:,:,:,yind,]
                    fvar = tmp

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

#--  create datm data files  --------------------------------------
##3
if create_datm == 1:
   
#--  create datm domain file
   f2 =  netcdf4.Dataset(fdatmdomain, 'r', format='NETCDF4')
   xc  = f2.variables['xc']
   yc  = f2.variables['yc']
   mask  = np.copy(f2.variables['mask'])
   
#--  convert coordinates to 1d
   lon=np.asarray(xc[0,:])
   lat=np.asarray(yc[:,0])
   im=lon.size
   jm=lat.size
   
   dlon = np.abs(lon[0] - lon[1])
   dlat = np.abs(lat[0] - lat[1])
   xind = np.argmin(abs(lon - plon))
   yind = np.argmin(abs(lat - plat))
   im_new=xind.size
   jm_new=yind.size
   
   mask=mask[yind,:]; mask=mask[xind]
   
   lonc=np.copy(xc) ; latc=np.copy(yc)
   lonc=lonc[yind,:]; lonc=lonc[xind]
   latc=latc[yind,:]; latc=latc[xind]
   #frac  = np.copy(f2.variables['frac'])
   area  = np.copy(f2.variables['area'])

   # if no frac on dataset, use all 1s
   frac=np.ones((jm,im))

   frac=frac[yind,:]; frac=frac[xind]     
   area=area[yind,:]; area=area[xind]
   
   
   nv = 4
    #--  Create vertex matrices  ----------------------------
   lonv=np.empty((nv), dtype='float64')
   latv=np.empty((nv), dtype='float64')
    #--  order is SW, SE, NE, NW
   lonv[0] = lonc - 0.5*dlon
   latv[0] = latc - 0.5*dlat
   lonv[1] = lonc + 0.5*dlon
   latv[1] = latc - 0.5*dlat
   lonv[2] = lonc + 0.5*dlon
   latv[2] = latc + 0.5*dlat
   lonv[3] = lonc - 0.5*dlon
   latv[3] = latc + 0.5*dlat
   
    #--  Check whether file exists  ---------------------------------
   command=['ls',fdatmdomain2]
   file_exists=subprocess.call(command,stderr=subprocess.PIPE)
   print command
   print file_exists
   
   w = netcdf4.Dataset(fdatmdomain2, 'w', format='NETCDF4')
   
   if file_exists > 0:
      print 'creating new file: ', fdatmdomain2
   else:
      print 'overwriting file: ', fdatmdomain2
   command='date "+%y%m%d"'
   x2=subprocess.Popen(command,stdout=subprocess.PIPE,shell='True')
   x=x2.communicate()
   timetag = x[0].strip()
   w.creation_date = timetag
   w.source_file = fdatmdomain
   w.title = 'CESM datm domain data'
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
   print 'successfully created datm domain file: ', fdatmdomain2
   f2.close()

   infile=[]
   outfile=[]
   for y in range(syr,eyr+1):
      ystr=str(y)
      for m in range(1,13):
         mstr=str(m) 
         if m < 10:
            mstr='0'+mstr

         dtag=ystr+'-'+mstr

         fsolar=dir_input_datm+solrdir+solrtag+dtag+'.nc'
         fsolar2=dir_output_datm+solrtag+tag+'.'+dtag+'.nc'
         fprecip=dir_input_datm+precdir+prectag+dtag+'.nc'
         fprecip2=dir_output_datm+prectag+tag+'.'+dtag+'.nc'
         ftpqw=dir_input_datm+tpqwldir+tpqwtag+dtag+'.nc'
         ftpqw2=dir_output_datm+tpqwtag+tag+'.'+dtag+'.nc'
         
         infile+=[fsolar,fprecip,ftpqw]
         outfile+=[fsolar2,fprecip2,ftpqw2]

   nm=len(infile)
   for n in range(nm):
      print outfile[n], '\n'
      file_in = infile[n]
      file_out = outfile[n]

      f2 =  netcdf4.Dataset(file_in, 'r', format='NETCDF4')
      global_attributes  = f2.ncattrs()
      variables = f2.variables
      dimensions = f2.dimensions

      #--  Open output file
      w = netcdf4.Dataset(file_out, 'w', format='NETCDF4')
      
      #--  Set global attributes
      for ga in global_attributes:
         setattr(w,ga,f2.getncattr(ga))
      #--  Set dimensions of output file
      for dim in dimensions.keys():
         print dim
         if dim == 'lon':
            w.createDimension(dim,int(im_new))
         elif dim == 'lat':
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
            if fdim == 'lon':
               if n == 0:
                  shp=np.asarray(fvar.shape)
                  shp[0] = 1
                  tmp=np.empty(shp, dtype='float64')
                  tmp[0,] = fvar[xind,]
                  fvar = tmp
               if n == 1:
                  shp=np.asarray(fvar.shape)
                  shp[1] = 1
                  tmp=np.empty(shp, dtype='float64')
                  tmp[:,0,] = fvar[:,xind,]
                  fvar = tmp
               if n == 2:
                  shp=np.asarray(fvar.shape)
                  shp[2] = 1
                  tmp=np.empty(shp, dtype='float64')
                  tmp[:,:,0,] = fvar[:,:,xind,]
                  fvar = tmp
               if n == 3:
                  shp=np.asarray(fvar.shape)
                  shp[3] = 1
                  tmp=np.empty(shp, dtype='float64')
                  tmp[:,:,:,0,] = fvar[:,:,:,xind,]
                  fvar = tmp
            if fdim == 'lat':
               if n == 0:
                  shp=np.asarray(fvar.shape)
                  shp[0] = 1
                  tmp=np.empty(shp, dtype='float64')
                  tmp[0,] = fvar[yind,]
                  fvar = tmp
               if n == 1:
                  shp=np.asarray(fvar.shape)
                  shp[1] = 1
                  tmp=np.empty(shp, dtype='float64')
                  tmp[:,0,] = fvar[:,yind,]
                  fvar = tmp
               if n == 2:
                  shp=np.asarray(fvar.shape)
                  shp[2] = 1
                  tmp=np.empty(shp, dtype='float64')
                  tmp[:,:,0,] = fvar[:,:,yind,]
                  fvar = tmp
               if n == 3:
                  shp=np.asarray(fvar.shape)
                  shp[3] = 1
                  tmp=np.empty(shp, dtype='float64')
                  tmp[:,:,:,0,] = fvar[:,:,:,yind,]
                  fvar = tmp

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
      print 'datm files written to: ', dir_output_datm

#--  create landuse data file  --------------------------------------
#46
if create_landuse == 1:
    f2 =  netcdf4.Dataset(flanduse, 'r', format='NETCDF4')
    global_attributes  = f2.ncattrs()
    variables = f2.variables
    dimensions = f2.dimensions

    #--  Check whether file exists  ---------------------------------
    command=['ls',flanduse2]
    file_exists=subprocess.call(command,stderr=subprocess.PIPE)
    if file_exists > 0:
        print 'creating new file: ', flanduse2
    else:
        print 'overwriting file: ', flanduse2

    #--  Open output file
    w = netcdf4.Dataset(flanduse2, 'w', format='NETCDF4')

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
                    shp=np.asarray(fvar.shape)
                    shp[0] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[0,] = fvar[xind,]
                    fvar = tmp
                if n == 1:
                    shp=np.asarray(fvar.shape)
                    shp[1] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,0,] = fvar[:,xind,]
                    fvar = tmp
                if n == 2:
                    shp=np.asarray(fvar.shape)
                    shp[2] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,:,0,] = fvar[:,:,xind,]
                    fvar = tmp
                if n == 3:
                    shp=np.asarray(fvar.shape)
                    shp[3] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,:,:,0,] = fvar[:,:,:,xind,]
                    fvar = tmp
            if fdim == 'lsmlat':
                if n == 0:
                    shp=np.asarray(fvar.shape)
                    shp[0] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[0,] = fvar[yind,]
                    fvar = tmp
                if n == 1:
                    shp=np.asarray(fvar.shape)
                    shp[1] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,0,] = fvar[:,yind,]
                    fvar = tmp
                if n == 2:
                    shp=np.asarray(fvar.shape)
                    shp[2] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,:,0,] = fvar[:,:,yind,]
                    fvar = tmp
                if n == 3:
                    shp=np.asarray(fvar.shape)
                    shp[3] = 1
                    tmp=np.empty(shp, dtype='float64')
                    tmp[:,:,:,0,] = fvar[:,:,:,yind,]
                    fvar = tmp

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
    print 'successfully created landuse data file: ', flanduse2

