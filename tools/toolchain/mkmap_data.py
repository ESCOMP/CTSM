#!/usr/bin/env python                                                                                                                                         

# 2020-12-13                Negin Sobhani


"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This Python script is part of the simplified toolchain for creating
the surface dataset for ctsm cases.

After creating the namelist/control file using ./gen_mksurf_namelist.py
using the  options for your desired case, you should run :

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
# - [x] Should read the correct namelist and check if it exist.
# - [x] Check if the desired input_dir exist
# - [ ] Check if the given namelist exist
# - [ ] Download if the data does not exist: instead of download here download in gen_mksur_namelist


from __future__ import print_function
 
# python standard libraries
import os
import sys
import re
import tqdm
import subprocess
import argparse
import logging
import shlex
import xarray as xr

 
from datetime import datetime
from getpass import getuser
 
#from gen_mksurf_namelist import get_parser, build_nl
 
__author__ = 'Negin Sobhani'
__email__ = 'negins@ucar.edu'


myname = getuser()
host_name = subprocess.run(['hostname'], 
                    stdout=subprocess.PIPE).stdout.decode('utf-8')


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
        parser.add_argument('--input_dir',
                    help='''
                    Path of your mesh files and input data.', 
                    [default: %(default)s]
                    ''',
                    action="store",
                    dest="input_dir",
                    default="/glade/p/cesm/cseg/inputdata/")
        parser.add_argument('--map_dir',
                    help='''
                    Path of your weight (mapping) files.', 
                    [default: %(default)s]
                    ''',
                    action="store",
                    dest="map_dir",
                    default=os.getcwd())
        parser.add_argument('-d','--debug', 
                    help='''
                    Debug mode will print more information
                    [default: %(default)s]
                    ''',
                    action="store_true", 
                    dest="debug", 
                    default=True)


        # if host_name contains cheyenne
        if ("cheyenne" or "casper" in host_name):
            parser.add_argument('--account',
                        help='''
                        Cheyenne or Casper project account name.', 
                        [default: %(default)s]
                        ''',
                        action="store",
                        dest="account_name",
                        default=os.environ.get('PROJECT'))

        return parser                                                                                                        


def get_pair(line):
    """
    Function for reading namelist in pair format
    for creating a dictionary.

    Args:
        line (str): String that include = sign

    Returns:
        key (str) : field or variable name
        value (str) : field or variable value

    """

    line = re.sub('=', '', line)
    key, sep, value = line.strip( ).partition(" ")
    value = value.strip()
    return key, value

def name_weightfile(src_res, src_mask, dst_res,dst_mask):
    """
    Function for creating weight file name based on
    source and destination files resolutions and masks.

    Args:
        src_res (str) : Source file resolution
        src_mask (str) : Source file mask
        dst_res (str) : Destination file resolution
        dst_mask (str) : Destination file mask

    Returns:
        weight_fname (str): mapping/weight file name.
    """
    weight_fname = 'map_'+src_res+'_'+src_mask+'_to_'+dst_res+'_'+dst_mask+'.nc'
    return weight_fname 

def read_nl (namelist):
    """
    Function for reading the namelist /control file.

    Args:
        namelist (str): namelist file name

    Raises:
        Error and exit if the namelist file does not exist

    Returns:
        nl_d (dict) : dictionary including namelist
            control file information.
    """

    #================================================================================
    #TODO: 
    # When the files in raw input is updated we can use the newly created namelist
    # for now because the files are not updated I am pointing to temporary namelist
    # that points to Sam's directory where the updated files exist!.
    #================================================================================

    namelist = "/glade/scratch/negins/ctsm_toolchain/tools/toolchain/surfdata_4x5_hist_78pfts_CMIP6_2000-2005_c201215.namelist"

    #-- check if the given namelist exist 
    if not os.path.isfile(namelist):
        print ("namelist filename :", namelist)
        print ("Error: namelist_fname :",namelist, "does not exist!")
        sys.exit('Please run ./gen_mksur_namelist.py to create the'+ 
                 'namelist for your case and make sure you are pointing '+ \
                 'to the created namelist using --namelist option.')

    nl_d={}

    with open(namelist) as data:
        for line in data:
            key, value = get_pair(line)
            nl_d[key] = value

    #-- printing out namelist values in debug mode
    logging.debug(' namelist  = \n\t'+namelist +"\n")
    for k, v in nl_d.items():
        logging.debug(" "+k + " : " + v)

    return nl_d

def parse_mesh_fname(mesh_fname):
    """
    This function parse mesh filename
    to find the corresponding resolution and mask.
    It assumes this information in in the mesh_fname.

    Args:
        mesh_fname (str) : Mesh filename

    Raises:
        Error if mesh file does not exists.

    Returns:
        res (str) : mesh file resolution from mesh filename
        mask (str) : mesh mask from mesh filename
    """

    #-- check if the mesj file exists 
    if not os.path.isfile(mesh_fname):
        print ("mesh filename :", mesh_fname)
        print ("Error: mesh_fname :",mesh_fname, "does not exist!")

        #TODO: What should they do if mesh file does not exist

    print ("mesh_fname : ", mesh_fname)

    mesh_fname = mesh_fname.rsplit('/')[-1]
    mesh_fname = mesh_fname.rsplit('_')

    res = mesh_fname[1]
    mask = mesh_fname[2]

    #-- verbose print of res and mask in debug mode.
    logging.debug('  ==> res  : '+ res)
    logging.debug('  ==> mask : '+ mask)

    return res, mask

# TODO []: error checking grid and mask between global attributes and ****.
# for src files only and not dst files.
# If we include this erro check, we should require the users to add 
# grid and landmask attributes to their own raw dataset if they are 
# using their own data.

class PbsScheduler:
    """
    A class to encapsulate the scheduler object compatible
    with PBS on the NCAR Cheyenne and Casper systems.

    ...

    Attributes
    ----------
    account : str
        Cheyenne or Casper project code or account code

    nproc : str
        Number of processors to request

    nnodes: str
        Number of nodes to request

    ppn: str
        Number of processors per node to request

    mem: str
        Cheyenne node memory to request in GB

    queue: str
        Cheyenne or Casper queue to use
 
    walltime : str
        The wallclock time in HH:MM:SS format.

    email : str
        user email address

    mail_events : str
        Specifies the set of conditions under which
        a notification about the job is sent:

            'n': no email
            'b': email before each job begins
            'e': email after each job terminates
            'a': email after each job is aborted


    Methods
    -------

    """

    def __init__(self, account, nproc, nnodes, mem, ppn, queue, walltime, email=None, mail_events=None):

        self.account = account
        self.nproc = nproc
        self.nnodes = nnodes
        self.mem = mem
        self.ppn = ppn
        self.queue = queue
        self.walltime = walltime
        self.email = email
        self.mail_events = mail_events
        #self.mail_events = "abe" if mail_events is not None 

    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item])
                    for item in sorted(self.__dict__)))


    def write_job_pbs (self, job_fname, cmd):
        self.job_fname = job_fname

        job_file= open (self.job_fname, 'w')

        self.job_template = ( \
                "#!/bin/tcsh\n"
                "#PBS -N "+self.job_fname+"\n"
                "#PBS -A "+self.account+"\n"
                "#PBS -l walltime="+self.walltime+"\n"
                "#PBS -q "+self.queue+"\n"
                "#PBS -j oe"+"\n\n"
                #"#PBS -m "+ self.mail_events+"\n"
                "#PBS -M "+ self.email+"\n"
                #"#PBS -l select=4:ncpus=16:mpiprocs=16"+"\n"
                #"#PBS -l select=4:ncpus=2:mpiprocs=2:mem=109GB"+"\n"
                "#PBS -l select="+self.nnodes+":ncpus="+self.nproc+":mpiprocs="+self.ppn+":mem="+self.mem+"\n"
                "set MPI_SHEPHERD=true"+"\n"
                "mpiexec_mpt "+cmd+"\n"
                )
        job_file.write(self.job_template)
        job_file.close()

    def schedule (self):
        #cmd = "qsub "+self.job_fname

        #print (cmd)
        #khar = subprocess.check_output(["qsub", self.job_fname])
        khar = subprocess.run(["qsub", self.job_fname],stdout=subprocess.PIPE)
        #os.system(cmd)

        print ("------")
        test = khar.stdout.decode('utf-8')
        print (test)
        print (re.findall(r'^\D*(\d+)', test))

        print ("------")
        

def main ():                                                                                                                                                  
    args = get_parser().parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    logging.debug('debugging in debug mode.')

    namelist = args.namelist
    input_dir = args.input_dir
    map_dir = args.map_dir

    #-- check if the given input_dir exist
    if not os.path.exists(input_dir):
        sys.exit('ERROR: \n'+
                 '\t input_dir does not exist on this machine. \n'+
                 '\t Please point to the correct raw_dir using --input_dir'+
                 'flag.')

    #-- read the namelist
    nl_d = read_nl(namelist)

    #-- dst mesh file
    print ('-----------------------')
    dst_mesh_file = nl_d['dst_mesh_file']
    print ('dst_mesh_file :',  dst_mesh_file)
    print ('-----------------------')

    dst_res, dst_mask = parse_mesh_fname(dst_mesh_file)


    #src mesh files
    src_flist = [
                'mksrf_fsoitex',
                'mksrf_forganic',
                'mksrf_flakwat',
                'mksrf_fwetlnd',
                'mksrf_fmax',
                'mksrf_fglacier', 
                'mksrf_fvocef',
                'mksrf_furbtopo',
                'mksrf_fgdp',
                'mksrf_fpeat',
                'mksrf_fsoildepth', 
                'mksrf_fabm',
                'mksrf_furban',
                'mksrf_ftopostats', 
                'mksrf_fvegtyp',
                'mksrf_fsoicol',
                'mksrf_flai',
                'mksrf_dyn_lu'
                'mksrf_fvic'
                ]

    for src_file in tqdm.tqdm(src_flist):
        src_fname = nl_d[src_file].strip()
        print (src_file, " : ", src_fname)

        #-- check if src_file exist
        if not os.path.isfile(src_fname):
            print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print ("src_fname: ", src_fname, "does not exist!")
            print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        else:
            ds = xr.open_dataset(src_fname)
            src_mesh_file = input_dir+ ds.src_mesh_file
            src_res, src_mask = parse_mesh_fname(src_mesh_file)

            lrg_args = " "
            if (src_res =="3minx3min" or "1km-merge-10min"):
                lrg_args=" --64bit_offset "


            w_fname = name_weightfile (src_res, src_mask, dst_res, dst_mask) 

            print ("src_mesh_file is : ", src_mesh_file)
            print ("weight file is : ", w_fname)

            # TODO : How to get esmf binary path when not on cheyenne.

            esmf_bin_path = "/glade/work/dunlap/ESMF-INSTALL/8.0.0bs38/bin/bing/Linux.intel.64.mpt.default/"
            esmf_regrid = esmf_bin_path + "/ESMF_RegridWeightGen"

            cmd = esmf_regrid + " --ignore_unmapped"

            cmd = cmd + " -s " + src_mesh_file + " -d " +dst_mesh_file
            cmd = cmd + " -m conserve "+" -w "+ w_fname
            cmd = cmd + lrg_args
            cmd = cmd + "  --src_type SCRIP  --dst_type SCRIP"

            
            account = "P93300606"
            nproc = "2"
            nnodes = "2"
            mem = "109GB"
            ppn = "16"
            queue = "regular"
            walltime = "12:00:00"
            email = "negins@ucar.edu"



            #scheduler class:
            pbs_job = PbsScheduler(account, nproc, nnodes, mem, ppn, queue, walltime, email)

            job_fname = "regrid_submit"+w_fname +".sh"

            print ("cmd:", cmd)
            pbs_job.write_job_pbs(job_fname, cmd)
            #pbs_job.schedule()
            #print (pbs_job)
            #exit()
            
            #job_fname = "regrid_submit.sh"
            #job_file= open (job_fname, 'w')

            # job scheduler for myself:
            #job_template = ( \
            #        "#!/bin/tcsh\n"
            #        "#PBS -N ctsm_nwp"+"\n"
            #        "#PBS -A P93300606"+"\n"
            #        "#PBS -l walltime=12:00:00"+"\n"
            #        "#PBS -q regular"+"\n"
            #        "#PBS -j oe"+"\n\n"
            #        #"#PBS -l select=4:ncpus=16:mpiprocs=16"+"\n"
            #        "#PBS -l select=4:ncpus=2:mpiprocs=2:mem=109GB"+"\n"
            #        "set MPI_SHEPHERD=true"+"\n"
            #        "mpiexec_mpt "+cmd+"\n"
            #        )


            #job_file.write(job_template)
            #job_file.close()

            #cmd = "qsub "+job_fname



            #if os.path.isfile(w_fname):
            #    print ('weight file ', w_fname, 'already exists!')
            #else:
            #    print (cmd)
            #    os.system(cmd)
            #print ('-----------------------')

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



