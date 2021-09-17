#! /usr/bin/env python3
"""
Run this code by using the following wrapper script:
../../tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

#  Import libraries
import os
import subprocess

from datetime import date
from getpass import getuser

import xarray as xr


def get_git_sha():
    """
    Returns Git short SHA for the currect directory.
    """
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip().decode()


class ModifyFsurdat:
    """
    A case to encapsulate function that modifies fsurdat file.

    ...

    Attributes
    ----------

    Methods
    -------
    add_tag_to_filename(filename, tag)
       add a tag and timetag to a filename ending with 
       [._]cYYMMDD.nc or [._]YYMMDD.nc 

    update_metadata(nc, filename)

    modify:
        Modify input surface dataset
    """

    def __init__(self, overwrite_single_pft, dominant_pft,
                 zero_nonveg_landunits, uniform_snowpack, no_saturation_excess):
        self.overwrite_single_pft = overwrite_single_pft
        self.dominant_pft = dominant_pft
        self.zero_nonveg_landunits = zero_nonveg_landunits
        self.uniform_snowpack = uniform_snowpack
        self.no_saturation_excess = no_saturation_excess


    @staticmethod
    def add_tag_to_filename(filename, tag):
        """
        Add a tag and replace timetag of a filename
        # Expects file to end with [._]cYYMMDD.nc or [._]YYMMDD.nc
        # Add the tag to just before that ending part
        # and change the ending part to the current time tag
        """

        basename = os.path.basename(filename)
        cend = -10

        if ( basename[cend] == "c" ):
           cend = cend - 1
        if ( (basename[cend] != ".") and (basename[cend] != "_") ):
           print ( "Trouble figuring out where to add tag to filename:" + filename )
           os.abort()

        today = date.today()
        today_string = today.strftime("%y%m%d")

        return(basename[:cend] + "_" + tag + "_c" + today_string + '.nc')


    @staticmethod
    def update_metadata(nc, filename):

        #update attributes
        today = date.today()
        today_string = today.strftime("%Y-%m-%d")

        #get git hash
        sha = get_git_sha()

        nc.attrs['Created_on'] = today_string
        nc.attrs['Created_by'] = getuser()
        nc.attrs['Created_with'] = os.path.abspath(__file__) + " -- " + sha
        nc.attrs['Created_from'] = filename

        #delete unrelated attributes if they exist
        del_attrs = ['source_code', 'SVN_url', 'hostname', 'history'
                     'History_Log', 'Logname', 'Host', 'Version',
                     'Compiler_Optimized']
        attr_list = nc.attrs

        for attr in del_attrs:
            if attr in attr_list:
                del(nc.attrs[attr])


    def modify(self):

        print ("Creating surface dataset file")
        filename = self.fsurf_in
        print("Open file: " + filename)
        f1 = xr.open_dataset(filename)

        # modify surface data properties
        if self.overwrite_single_pft:
            f1['PCT_NAT_PFT'][:,:,:] = 0
            f1['PCT_NAT_PFT'][:,:,self.dominant_pft] = 100
        if self.zero_nonveg_landunits:
            f1['PCT_NATVEG'][:,:]  = 100
            f1['PCT_CROP'][:,:]    = 0
            f1['PCT_LAKE'][:,:]    = 0.
            f1['PCT_WETLAND'][:,:] = 0.
            f1['PCT_URBAN'][:,:,]   = 0.
            f1['PCT_GLACIER'][:,:] = 0.
        if self.uniform_snowpack:
            f1['STD_ELEV'][:,:] = 20.
        if self.no_saturation_excess:
            f1['FMAX'][:,:] = 0.

        # specify dimension order 
        f1= f1.transpose(u'time', u'cft', u'lsmpft', u'natpft', u'nglcec', u'nglcecp1', u'nlevsoi', u'nlevurb', u'numrad', u'numurbl', 'lsmlat', 'lsmlon')
        
        #update attributes
        self.update_metadata(f1, filename)

        # mode 'w' overwrites file if it exists
        f1.to_netcdf(path=self.fsurf_out, mode='w')
        print('Successfully created file (fsurf_out) :' + self.fsurf_out)
        f1.close()


