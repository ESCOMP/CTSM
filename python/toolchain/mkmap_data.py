#!/usr/bin/env python                                                                                                                                         
# -*- coding: utf-8 -*-
# 2020-12-13                Negin Sobhani
 
"""
Description:
This script includes mkmap_data.py ***
"""
 
from __future__ import print_function
 
# python standard libraries
import os
import sys
import re
import subprocess
import argparse
import logging
import shlex
import xarray as xr

 
from datetime import datetime
 
from gen_mksurf_namelist import get_parser, build_nl
 
__author__ = 'Negin Sobhani'
__email__ = 'negins@ucar.edu'

#def create_job (****):
    #"pbs script that gets submited" 

def submit_job (file_name):
    qsub file_name

def get_pair(line):
    line = re.sub('=', '', line)
    key, sep, value = line.strip( ).partition(" ")
    value = value.strip()
    return key, value

def name_weightfile(src_res, src_mask, dst_res,dst_mask):
    weight_fname = 'map_'+src_res+'_'+src_mask+'_to_'+dst_res+'_'+dst_mask+'.nc'
    return weight_fname 

def read_nl (nl_fname):
    nl_d={}
    with open("surfdata_4x5_hist_78pfts_CMIP6_2000-2005_c201215.namelist") as data:
        for line in data:
            key, value = get_pair(line)
            nl_d[key] = value
    print ('gitdescribe   : ', nl_d['gitdescribe'])
    print ('use_transient : ', nl_d['use_transient'].strip())

    print ('start_year    : ', nl_d['start_year'])
    print ('end_year      : ', nl_d['end_year'])
    return nl_d

def parse_mesh_fname(mesh_fname):
    """
    This function parse mesh filename and parse it.
    """
    mesh_fname = mesh_fname.rsplit('/')[-1]
    mesh_fname = mesh_fname.rsplit('_')

    res = mesh_fname[1]
    mask = mesh_fname[2]

    return res, mask


def main ():                                                                                                                                                  
    args = get_parser().parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    logging.debug('debugging in debug mode.')

    # TODO: should read the correct namelist
    nl_fname = "surfdata_4x5_hist_78pfts_CMIP6_2000-2005_c201215.namelist"

    # read the namelist
    nl_d = read_nl(nl_fname)
    input_path = "/glade/p/cesm/cseg/inputdata/"

    
    print ('-----------------------')
    dst_mesh_file = nl_d['dst_mesh_file']
    dst_res, dst_mask = parse_mesh_fname(dst_mesh_file)

    print ('dst_mesh_file :',  dst_mesh_file)


    src_flist = [
                'mksrf_fsoitex',
                'mksrf_forganic',
                #'mksrf_flakwat',
                'mksrf_fwetlnd',
                'mksrf_fmax',
                #'mksrf_fglacier',
                'mksrf_fvocef',
                'mksrf_furbtopo',
                'mksrf_fgdp',
                'mksrf_fpeat',
                #'mksrf_fsoildepth',
                'mksrf_fabm',
                'mksrf_furban',
                #'mksrf_ftopostats',
                'mksrf_fvegtyp',
                'mksrf_fsoicol',
                'mksrf_flai'
                ]

    for src_file in src_flist:
        src_fname = nl_d[src_file].strip()
        print (src_file, " : ", src_fname)
        ds = xr.open_dataset(src_fname)
        src_mesh_file = input_path+ ds.src_mesh_file
        src_res, src_mask = parse_mesh_fname(src_mesh_file)

        lrg_args = " "
        if (src_res =="3minx3min"):
            print ("src_res is: ", src_res)
            lrg_args=" --64bit_offset "

        print ('src_mesh_file is : ', src_mesh_file)


        w_fname = name_weightfile (src_res, src_mask, dst_res, dst_mask) 

        print ("weight_fname is : ", w_fname)

        esmf_bin_path = "/glade/work/dunlap/ESMF-INSTALL/8.0.0bs38/bin/bing/Linux.intel.64.mpt.default/"
        esmf_regrid = esmf_bin_path + "/ESMF_RegridWeightGen"

        cmd = esmf_regrid + " --ignore_unmapped"

        cmd = cmd + " -s " + src_mesh_file + " -d " +dst_mesh_file
        cmd = cmd + " -m conserve "+" -w "+ w_fname
        cmd = cmd + lrg_args
        cmd = cmd + "  --src_type SCRIP  --dst_type SCRIP"

        if os.path.isfile(w_fname):
            print ('weight file ', w_fname, 'already exists!')
        else:
            print (cmd)
            os.system(cmd)
        print ('-----------------------')

    '''
     /glade/work/dunlap/ESMF-INSTALL/8.0.0bs38/bin/bing/Linux.intel.64.mpt.default/ESMF_RegridWeightGen
     --ignore_unmapped -s /glade/p/cesm/cseg/inputdata/lnd/clm2/mappingdata/grids/SCRIPgrid_0.5x0.5_AVHRR_c110228.nc  -d
     /glade/p/cesm/cseg/inputdata/lnd/clm2/mappingdata/grids/SCRIPgrid_10x15_nomask_c110308.nc -m conserve -w
     map_0.5x0.5_AVHRR_to_10x15_nomask_aave_da_c201215.nc --src_type SCRIP  --dst_type SCRIP 

    /glade/work/dunlap/ESMF-INSTALL/8.0.0bs38/bin/bing/Linux.intel.64.mpt.default//ESMF_RegridWeightGen -s
    /glade/p/cesm/cseg/inputdata/lnd/clm2/mappingdata/grids/SCRIPgrid_0.25x0.25_MODIS_c170321.nc -d
    /glade/p/cesm/cseg/inputdata/lnd/clm2/mappingdata/grids/SCRIPgrid_10x15_nomask_c110308.nc -m conserve  -w
    map_0.25x0.25_MODIS_to_10x15nomask.nc  --src_type SCRIP  --dst_type SCRIP
    '''

if __name__ == "__main__":
    main()



