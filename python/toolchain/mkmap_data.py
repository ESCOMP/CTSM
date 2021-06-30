#!/usr/bin/env python                                                                                                                                         

# 2020-12-13                Negin Sobhani


"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This Python script is part of the simplified toolchain for creating
the surface dataset for ctsm cases.

After creating the namelist/control file using ./gen_mksurf_namelist.py
with options for your desired case, you should run :

./mkmap_data.py --namelist [namelist from ./gen_mksurf_namelist.py]

to create mapping files.


This python script is the alternative code to mkmapdata.sh with
 few differences:

1. It reads namelist file for DST mesh file.

2. For raw dataset in the namelist it finds their SRC mesh file and mask from
the necdf metadata.

3. Check if the weight (mapping file) already exists or not. 

4. If it does not exist, it creates the mapping (weight) file. 
 
-------------------------------------------------------------------
Instructions for running on Cheyenne/Casper:
 
load the following into your local environment:
 
    module load python
    ncar_pylib
-------------------------------------------------------------------
To see the available options:
    ./mkmap_data.py --help
 
To run the script:
    ./mkmap_data.py --namelist ${namelist from ./gen_mksurf_namelist.py}
 
To remove NPL from your environment on Cheyenne/Casper:
    deactivate
-------------------------------------------------------------------
"""

#TODO (NS):
# - [ ] Add docs and notes for functions.
# - [ ] Job submission classes for different machines. 
# - [x]  should read the correct namelist and check if it exist.
# - [ ] Check if the desired input_path exist
 
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
 
#from gen_mksurf_namelist import get_parser, build_nl
 
__author__ = 'Negin Sobhani'
__email__ = 'negins@ucar.edu'

def get_parser():
        """ 
            Get parser object for this script.
        """
        parser = argparse.ArgumentParser(description=__doc__,
                           formatter_class=argparse.RawDescriptionHelpFormatter)
 
        parser.print_usage = parser.print_help

        parser.add_argument('--namelist',
                    help='namelist for the ctsm case created by gen_mksurf_namelist.py ', 
                    action="store",
                    dest="namelist",
                    required=True)
        parser.add_argument('-d','--debug', 
                    help='Debug mode will print more information. ', 
                    action="store_true", 
                    dest="debug", 
                    default=False)
        parser.add_argument('--input_dir',
                    help='''
                    Path of your mesh files and input data.', 
                    [default: %(default)s]
                    ''',
                    action="store",
                    dest="input_path",
                    default="/glade/p/cesm/cseg/inputdata/")

        return parser                                                                                                        



def get_pair(line):
    line = re.sub('=', '', line)
    key, sep, value = line.strip( ).partition(" ")
    value = value.strip()
    return key, value

def name_weightfile(src_res, src_mask, dst_res,dst_mask):
    weight_fname = 'map_'+src_res+'_'+src_mask+'_to_'+dst_res+'_'+dst_mask+'.nc'
    return weight_fname 

def read_nl (namelist):
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

# I removed the sections of schedulers because they were not working.
#def create_job (****):
    #"pbs script that gets submited" 

#def submit_job (file_name):
#    qsub file_name


def main ():                                                                                                                                                  
    args = get_parser().parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    logging.debug('debugging in debug mode.')
    namelist = args.namelist
    input_path = args.input_path


    #-- check if the given namelist exist 
    if not os.path.isfile(namelist):
        print ("namelist filename :", namelist)
        print ("Error: namelist_fname :",namelist, "does not exist!")
        sys.exit('Please run ./gen_mksur_namelist.py to create the'+ 
                 'namelist for your case and make sure you are pointing '+ \
                 'to the created namelist using --namelist option.') 

    #-- check if the given input path exist
    if not os.path.exists(input_path):
        sys.exit('ERROR: \n'+
                 '\t input_path does not exist on this machine. \n'+
                 '\t Please point to the correct raw_dir using --input_path'+
                 'flag.')


    # read the namelist
    nl_d = read_nl(namelist)

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



