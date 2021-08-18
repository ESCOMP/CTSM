from base_case import BaseCase
import os
import numpy as np
import xarray as xr
from datetime import date

class SinglePointCase (BaseCase):
    """
    A case to encapsulate single point cases.
 
    ...
 
    Attributes
    ----------
    plat : float
        latitude
    plon : float
        longitude
    site_name: str -- default = None
        Site name
 
    Methods
    -------
    create_tag
        create a tag for single point which is the site name
        or the "lon-lat" format if the site name does not exist.
 
    create_domain_at_point
        Create domain file at a single point.
    create_landuse_at_point:
        Create landuse file at a single point.
    create_surfdata_at_point:
        Create surface dataset at a single point.
    create_datmdomain_at_point:
        Create DATM domain file at a single point.
    """
 
    def __init__(self, plat, plon,site_name,
                 create_domain, create_surfdata, create_landuse, create_datm,
                 overwrite_single_pft, dominant_pft, zero_nonveg_landunits,
                 uniform_snowpack, no_saturation_excess):
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm)
        self.plat = plat
        self.plon = plon
        self.site_name = site_name
        self.overwrite_single_pft = overwrite_single_pft
        self.dominant_pft = dominant_pft
        self.zero_nonveg_landunits = zero_nonveg_landunits
        self.uniform_snowpack = uniform_snowpack
        self.no_saturation_excess = no_saturation_excess

    def create_tag(self):
        if self.site_name:
            self.tag = self.site_name
        else:
             self.tag=str(self.plon)+'_'+str(self.plat)
 
    @staticmethod
    def create_fileout_name( filename,tag):
 
        basename = os.path.basename(filename)
        items = basename.split('_')
        today = date.today()
        today_string = today.strftime("%y%m%d")
        new_string = items[0]+"_"+items[2]+"_"+items[3]+"_"+ items[4] \
                    +"_"+items[5]+"_"+items[6]+"_"+tag+"_c"+today_string+".nc"
        return new_string
 
    def create_domain_at_point (self):
        print( "----------------------------------------------------------------------")
        print ("Creating domain file at ", self.plon, self.plat)
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fdomain_in, 'xc','yc','ni','nj')
        # extract gridcell closest to plon/plat
        f3 = f2.sel(ni=self.plon,nj=self.plat,method='nearest')
        # expand dimensions
        f3 = f3.expand_dims(['nj','ni'])
 
        #update attributes
        self.update_metadata(f3)
        f3.attrs['Created_from'] = self.fdomain_in
 
        wfile=self.fdomain_out
        f3.to_netcdf(path=wfile, mode='w')
        print('Successfully created file (fdomain_out)'+self.fdomain_out)
        f2.close(); f3.close()
 
 
 
    def create_landuse_at_point (self):
        print( "----------------------------------------------------------------------")
        print ("Creating landuse file at ", self.plon, self.plat, ".")
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fluse_in, 'LONGXY','LATIXY','lsmlon','lsmlat')
        # extract gridcell closest to plon/plat
        f3 = f2.sel(lsmlon=self.plon,lsmlat=self.plat,method='nearest')
 
        # expand dimensions
        f3 = f3.expand_dims(['lsmlat','lsmlon'])
        # specify dimension order 
        #f3 = f3.transpose('time','lat','lon')
        f3 = f3.transpose(u'time', u'cft', u'natpft', u'lsmlat', u'lsmlon')
        #f3['YEAR'] = f3['YEAR'].squeeze()
 
        # revert expand dimensions of YEAR
        year = np.squeeze(np.asarray(f3['YEAR']))
        x = xr.DataArray(year, coords={'time':f3['time']}, dims='time', name='YEAR')
        x.attrs['units']='unitless'
        x.attrs['long_name']='Year of PFT data'
        f3['YEAR'] = x 
 
        #update attributes
        self.update_metadata(f3)
        f3.attrs['Created_from'] = self.fluse_in
 
        wfile = self.fluse_out
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode='w')
        print('Successfully created file (luse_out)'+self.fluse_out,".")
        f2.close(); f3.close()


    def create_surfdata_at_point(self):
        print( "----------------------------------------------------------------------")
        print ("Creating surface dataset file at ", self.plon, self.plat, ".")
        # create 1d coordinate variables to enable sel() method
        filename = self.fsurf_in
        f2 = self.create_1d_coord(filename, 'LONGXY','LATIXY','lsmlon','lsmlat')
        # extract gridcell closest to plon/plat
        f3 = f2.sel(lsmlon=self.plon,lsmlat=self.plat,method='nearest')
        # expand dimensions
        f3 = f3.expand_dims(['lsmlat','lsmlon']).copy(deep=True)
 
        # modify surface data properties
        if self.overwrite_single_pft:
            f3['PCT_NAT_PFT'][:,:,:] = 0
            f3['PCT_NAT_PFT'][:,:,self.dominant_pft] = 100
        if self.zero_nonveg_landunits:
            f3['PCT_NATVEG'][:,:]  = 100
            f3['PCT_CROP'][:,:]    = 0
            f3['PCT_LAKE'][:,:]    = 0.
            f3['PCT_WETLAND'][:,:] = 0.
            f3['PCT_URBAN'][:,:,]   = 0.
            f3['PCT_GLACIER'][:,:] = 0.
        if self.uniform_snowpack:
            f3['STD_ELEV'][:,:] = 20.
        if self.no_saturation_excess:
            f3['FMAX'][:,:] = 0.
 
        # specify dimension order 
        #f3 = f3.transpose(u'time', u'cft', u'natpft', u'lsmlat', u'lsmlon')
        f3 = f3.transpose(u'time', u'cft', u'lsmpft', u'natpft', u'nglcec', u'nglcecp1', u'nlevsoi', u'nlevurb', u'numrad', u'numurbl', 'lsmlat', 'lsmlon')
        
        #update attributes
        self.update_metadata(f3)
        f3.attrs['Created_from'] = self.fsurf_in 
        del(f3.attrs['History_Log'])
        # mode 'w' overwrites file
        f3.to_netcdf(path=self.fsurf_out, mode='w')
        print('Successfully created file (fsurf_out) :'+self.fsurf_out)
        f2.close(); f3.close()

    def create_datmdomain_at_point(self): 
        print( "----------------------------------------------------------------------")
        print("Creating DATM domain file at ", self.plon, self.plat, ".") 
        # create 1d coordinate variables to enable sel() method
        filename = self.fdatmdomain_in
        f2 = self.create_1d_coord(filename, 'xc','yc','ni','nj')
        # extract gridcell closest to plon/plat
        f3 = f2.sel(ni=self.plon,nj=self.plat,method='nearest')
        # expand dimensions
        f3 = f3.expand_dims(['nj','ni'])
        wfile=self.fdatmdomain_out
        #update attributes
        self.update_metadata(f3)
        f3.attrs['Created_from'] = self.fdatmdomain_in 
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode='w')
        print('Successfully created file (fdatmdomain_out) :'+self.fdatmdomain_out)
        f2.close(); f3.close()

    def extract_datm_at(self, file_in, file_out):
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(file_in, 'LONGXY','LATIXY','lon','lat')
        # extract gridcell closest to plon/plat
        f3  = f2.sel(lon=self.plon,lat=self.plat,method='nearest')
        # expand dimensions
        f3 = f3.expand_dims(['lat','lon'])
        # specify dimension order 
        f3 = f3.transpose(u'scalar','time','lat','lon')
        
        #update attributes
        self.update_metadata(f3)
        f3.attrs['Created_from'] = file_in 
        # mode 'w' overwrites file
        f3.to_netcdf(path=file_out, mode='w')
        print('Successfully created file :'+ file_out)
        f2.close(); f3.close()
        
    def create_datm_at_point(self):
        print( "----------------------------------------------------------------------")
        print("Creating DATM files at ", self.plon, self.plat, ".")
        #--  specify subdirectory names and filename prefixes
        solrdir = 'Solar/'
        precdir = 'Precip/'
        tpqwldir = 'TPHWL/'
        prectag = 'clmforc.GSWP3.c2011.0.5x0.5.Prec.'
        solrtag = 'clmforc.GSWP3.c2011.0.5x0.5.Solr.'
        tpqwtag = 'clmforc.GSWP3.c2011.0.5x0.5.TPQWL.'
        
        #--  create data files  
        infile=[]
        outfile=[]
        for y in range(self.datm_syr,self.datm_eyr+1):
          ystr=str(y)
          for m in range(1,13):
             mstr=str(m) 
             if m < 10:
                mstr='0'+mstr
        
             dtag=ystr+'-'+mstr
        
             fsolar=self.dir_input_datm+solrdir+solrtag+dtag+'.nc'
             fsolar2=self.dir_output_datm+solrtag+self.tag+'.'+dtag+'.nc'
             fprecip=self.dir_input_datm+precdir+prectag+dtag+'.nc'
             fprecip2=self.dir_output_datm+prectag+self.tag+'.'+dtag+'.nc'
             ftpqw=self.dir_input_datm+tpqwldir+tpqwtag+dtag+'.nc'
             ftpqw2=self.dir_output_datm+tpqwtag+self.tag+'.'+dtag+'.nc'
        
             infile+=[fsolar,fprecip,ftpqw]
             outfile+=[fsolar2,fprecip2,ftpqw2]
        
        nm=len(infile)
        for n in range(nm):
            print(outfile[n])
            file_in = infile[n]
            file_out = outfile[n]
            self.extract_datm_at(file_in, file_out)
                                                                                                                               
        
        print('All DATM files are created in: '+self.dir_output_datm)



