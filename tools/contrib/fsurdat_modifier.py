#! /usr/bin/env python
"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|

Instructions for running on Cheyenne/Casper:

load the following into your local environment
    module load python
    ncar_pylib

-------------------------------------------------------------------
To see available options:
    ./fsurdat_modifier.py --help
-------------------------------------------------------------------

This script reads a surface dataset (fsurdat file) and outputs a
modified copy of the same file.

-------------------------------------------------------------------
To run the script:
    ./fsurdat_modifier.py
 
To remove NPL (ncar_pylib) from your environment on Cheyenne/Casper:
    deactivate
-------------------------------------------------------------------

"""

#  Import libraries
from __future__ import print_function

import sys
import os
import string
import logging
import subprocess
import argparse

import xarray as xr

from datetime import date
from getpass import getuser
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

myname = getuser()

def get_parser():                                                                                                                                                                   
        """
        Get parser object for this script.

        Command-line inputs
        - fsurdat_in: input file (str)
        - fsurdat_out: output file (str)
        - variable: variable to modify (str)
        - value: value assigned to the variable to be modified (float)
        """
        parser = ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)

        parser.print_usage = parser.print_help

        parser.add_argument('--dom_pft', 
                    help='Dominant PFT if overwrite_single_pft = .true. [default: %(default)s] ', 
                    action="store", 
                    dest="dom_pft",
                    type =int,
                    default=7)
        parser.add_argument('--uni_snow', 
                    help='Turn on the flag to create uniform snowpack. [default: %(default)s]', 
                    action="store", 
                    dest="uni_snow",
                    default=False)
        parser.add_argument('--overwrite_single_pft', 
                    help='Turn on the flag to make the whole grid 100%% single PFT. [default: %(default)s]', 
                    action="store", 
                    dest="overwrite_single_pft",
                    default=False)
        parser.add_argument('--zero_nonveg', 
                    help='Set all non-vegetation landunits to zero. [default: %(default)s]', 
                    action="store", 
                    dest="zero_nonveg",
                    type =bool,
                    default=False)
        parser.add_argument('--no_saturation_excess', 
                    help='Turn off saturation excess. [default: %(default)s]', 
                    action="store", 
                    dest="no_saturation_excess",
                    type =bool,
                    default=False)
        parser.add_argument('--fsurdat_in', 
                    help = 'Input surface dataset. [default: %(default)s]', 
                    action = "store", 
                    dest = "fsurdat_in",
                    type = str,
                    default = "/glade/p/cesmdata/cseg/inputdata/lnd/clm2/surfdata_map/release-clm5.0.18/surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr2000_c190214.nc")
        parser.add_argument('--fsurdat_out', 
                    help = 'Output surface dataset. [default: %(default)s]', 
                    action = "store", 
                    dest = "fsurdat_out",
                    type = str,
                    default = "/glade/scratch/" + myname + "/surfdata_0.9x1.25_hist_78pfts_CMIP6_simyr2000_c190214_modified.nc")

        return parser


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
        super().__init__()
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
        nc.attrs['Created_by'] = myname
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


def setup_logging(log_file, log_level):
    """
    Setup logging to log to console and log file.
    """

    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # setup log file
    one_mb = 1000000
    handler = logging.handlers.RotatingFileHandler(log_file, maxBytes=one_mb , backupCount=10)

    fmt = logging.Formatter(
            '%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')

    handler.setFormatter(fmt)
    root_logger.addHandler(handler)

    # setup logging to console
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(fmt)
    root_logger.addHandler(stream_handler)

    # redirect stdout/err to log file
    StreamToLogger.setup_stdout()
    StreamToLogger.setup_stderr()


class StreamToLogger(object):
    """
    Custom class to log all stdout and stderr streams.
    modified from:
    https://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/
    """
    def __init__(self, stream, logger, log_level=logging.INFO,
                 also_log_to_stream=False):
        self.logger = logger
        self.stream = stream
        self.log_level = log_level
        self.linebuf = ''
        self.also_log_to_stream = also_log_to_stream

    @classmethod
    def setup_stdout(cls, also_log_to_stream=True):
        """
        Setup logger for stdout
        """
        stdout_logger = logging.getLogger('STDOUT')
        sl = StreamToLogger(sys.stdout, stdout_logger, logging.INFO, also_log_to_stream)
        sys.stdout = sl

    @classmethod
    def setup_stderr(cls, also_log_to_stream=True):
        """
        Setup logger for stdout
        """
        stderr_logger = logging.getLogger('STDERR')
        sl = StreamToLogger(sys.stderr, stderr_logger, logging.ERROR, also_log_to_stream)
        sys.stderr = sl

    def write(self, buf):
        temp_linebuf = self.linebuf + buf 
        self.linebuf = ''
        for line in temp_linebuf.splitlines(True):
            if line[-1] == '\n':
                self.logger.log(self.log_level, line.rstrip())
            else:
                self.linebuf += line

    def flush(self):
        if self.linebuf != '': 
            self.logger.log(self.log_level, self.linebuf.rstrip()) 
        self.linebuf = ''


def main ():

    # Parse arguments from the command line
    args = get_parser().parse_args()

    # Set up logging
    today = date.today()
    today_string = today.strftime("%Y%m%d")
    pwd = os.getcwd()
    log_file = os.path.join(pwd, today_string + '.log')
    log_level = logging.DEBUG
    setup_logging(log_file, log_level)
    log = logging.getLogger(__name__)

    print("User = " + myname)
    print("Current directory = " + pwd)

    #--  Modify landunit structure
    overwrite_single_pft = args.overwrite_single_pft
    dominant_pft         = args.dom_pft
    zero_nonveg_landunits= args.zero_nonveg
    uniform_snowpack     = args.uni_snow
    no_saturation_excess = args.no_saturation_excess

    #--  Create ModifyFsurdat Object
    modify_fsurdat = ModifyFsurdat(overwrite_single_pft, dominant_pft,
                                   zero_nonveg_landunits, uniform_snowpack,
                                   no_saturation_excess)
    print(modify_fsurdat)

    #--  Set input and output filenames
    fsurf_in = args.fsurdat_in
    fsurf_out = args.fsurdat_out

    modify_fsurdat.fsurf_in = fsurf_in
    modify_fsurdat.fsurf_out = fsurf_out

    print ("fsurf_in   :",fsurf_in)  #
    print ("fsurf_out  :",fsurf_out) #

    #--  Create CTSM surface data file
    modify_fsurdat.modify()

    print( "Successful completion of script." )
    exit()


if __name__ == "__main__":
    main()

