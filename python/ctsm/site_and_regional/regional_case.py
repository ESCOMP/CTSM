from ctsm.site_and_regional.base_case import BaseCase

import numpy as np
import xarray as xr

class RegionalCase (BaseCase):
    """
    A case to encapsulate regional cases.
    """
 
    def __init__(self, lat1, lat2, lon1, lon2, reg_name,
                 create_domain, create_surfdata, create_landuse, create_datm):
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm)
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        self.reg_name = reg_name
 
    def create_tag(self):
        if self.reg_name:
            self.tag = self.reg_name
        else:
            self.tag=str(self.lon1)+'-'+str(self.lon2)+'_'+str(self.lat1)+'-'+str(self.lat2)
 
    def create_domain_at_reg (self):
        #print ("Creating domain file at region", self.lon1+"-"+self.lat2,self.lat1+"-"+self.lat2)
        print ("Creating domain file at region:", self.tag)
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fdomain_in, 'xc','yc','ni','nj')
        lat = f2['lat']
        lon = f2['lon']
        # subset longitude and latitude arrays
        xind=np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind=np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3=f2.isel(nj=yind,ni=xind)
 
        #update attributes
        self.update_metadata(f3)
        f3.attrs['Created_from'] = self.fdomain_in
 
        wfile=self.fdomain_out
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode='w')
        print('Successfully created file (fdomain_out)'+self.fdomain_out)
        f2.close(); f3.close()

 
 
    def create_surfdata_at_reg(self):
        #print ("Creating surface dataset file at region", self.lon1+"-"+self.lat2,self.lat1+"-"+self.lat2)
        print ("Creating surface dataset file at region:", self.tag)
        # create 1d coordinate variables to enable sel() method
        filename = self.fsurf_in
        f2 = self.create_1d_coord(filename, 'LONGXY','LATIXY','lsmlon','lsmlat')
        lat = f2['lat']
        lon = f2['lon']
        # subset longitude and latitude arrays
        xind=np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind=np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3=f2.isel(lsmlat=yind,lsmlon=xind)
 
        #update attributes
        self.update_metadata(f3)
        f3.attrs['Created_from'] = self.fsurf_in
 
        # mode 'w' overwrites file
        f3.to_netcdf(path=self.fsurf_out, mode='w')
        print('created file (fsurf_out)'+self.fsurf_out)
        #f1.close();
        f2.close(); f3.close()
 
 
    def create_landuse_at_reg (self):
        #print ("Creating surface dataset file at region", self.lon1+"-"+self.lat2,self.lat1+"-"+self.lat2)
        print ("Creating surface dataset file at region:",self.tag)
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fluse_in, 'LONGXY','LATIXY','lsmlon','lsmlat')
        lat = f2['lat']
        lon = f2['lon']
        # subset longitude and latitude arrays
        xind=np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind=np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3=f2.isel(lsmlat=yind,lsmlon=xind)
 
        #update attributes
        self.update_metadata(f3)
        f3.attrs['Created_from'] = self.fluse_in
 
        wfile=self.fluse_out
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode='w')
        print('Successfully created file (fdomain_out)'+self.fdomain_out)
        f2.close(); f3.close()

 
 



