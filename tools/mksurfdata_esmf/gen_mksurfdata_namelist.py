#!/usr/bin/env python3

import sys, os, shutil
import xml.etree.ElementTree as ET
import logging
import argparse, textwrap
import subprocess
from datetime import datetime

logger = logging.getLogger(__name__)

# valid options for SSP scenarios and pft years:
valid_opts = {"ssp_rcp": ["hist","SSP1-2.6","SSP3-7.0","SSP5-3.4","SSP2-4.5","SSP1-1.9","SSP4-3.4","SSP4-6.0","SSP5-8.5"]}

def get_parser():
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help

    parser.add_argument(
        "--start-year",
        help = textwrap.dedent('''\
               Simulation start year. 
               [Required]'''),
        action="store",
        dest="start_year",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--end-year",
        help = textwrap.dedent('''\
               Simulation end year. 
               [Required]'''),
        action="store",
        dest="end_year",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--res",
        help="""
            model resolution [default: %(default)s]
            To see available supported resolutions, simply invoke this command
            with a --res unknown opion
            """,
        action="store",
        dest="res",
        required=False,
        default="4x5",
    )
    parser.add_argument(
        "--glc-nec",
        help="""
            Number of glacier elevation classes to use. [default: %(default)s]
            """,
        action="store",
        dest="glc_nec",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--rundir",
        help="""
            Directory to run in.
            [default: %(default)s]
            """,
        action="store",
        dest="run_dir",
        required=False,
        default=os.getcwd(),
    )
    parser.add_argument(
        "--ssp-rcp",
        help="""
            Shared Socioeconomic Pathway and Representative
            Concentration Pathway Scenario name(s).
            [default: %(default)s]
            """,
        action="store",
        dest="ssp_rcp",
        required=False,
        choices=valid_opts["ssp_rcp"],
        default="hist",
    )
    parser.add_argument(
        "--rawdata-dir",
        help="""
            /path/of/root/of/input/data',
            [default: %(default)s]
            """,
        action="store",
        dest="input_path",
        default="/glade/p/cesm/cseg/inputdata/",
    )
    parser.add_argument(
        "--vic",
        help="""
            Flag for adding the fields required for the VIC model. 
            [default: %(default)s]
            """,
        action="store_true",
        dest="vic_flag",
        default=False,
    )
    parser.add_argument(
        "--glc",
        help="""
            Flag for adding the optional 3D glacier fields for verification of the glacier model.
            [default: %(default)s]
            """,
        action="store_true",
        dest="glc_flag",
        default=False,
    )
    parser.add_argument(
        "--hires_pft",
        help="""
            If you want to use the high-resolution pft dataset rather
            than the default lower resolution dataset.
            (Low resolution is at quarter-degree, high resolution at 3-minute)
            [Note: hires only available for 1850 and 2005.]
            [default: %(default)s]
            """,
        action="store_true",
        dest="hres_flag",
        default=False,
    )
    parser.add_argument(
        "--nocrop",
        help="""
            Create datasets with the extensive list of prognostic crop types.
            [default: %(default)s]
            """,
        action="store_false",
        dest="crop_flag",
        default=True,
    )
    parser.add_argument(
        "-f",
        "--fast",
        help="Toggle fast mode which does not user the large mapping file",
        action="store_true",
        dest="fast_flag",
        default=False,
    )
    parser.add_argument(
        "--potveg_flag",
        help="""
            Use Potential Vegetation for pft_years
            [default: %(default)s]
            """,
        action="store_false",
        dest="potveg_flag",
        default=False,
    )
    parser.add_argument(
        "--merge_gis",
        help="""
        If you want to use the glacier dataset that merges in
        the Greenland Ice Sheet data that CISM uses (typically
        used only if consistency with CISM is important)
        [default: %(default)s]
        """,
        action="store",
        dest="merge_gis",
        choices=["on","off"],
        default="off",
    )
    return parser

def main ():

    args = get_parser().parse_args()

    start_year = args.start_year
    end_year = args.end_year
    res = args.res
    ssp_rcp = args.ssp_rcp
    input_path = args.input_path
    crop_flag = args.crop_flag
    vic_flag = args.vic_flag
    glc_flag = args.glc_flag
    potveg = args.potveg_flag
    glc_nec = args.glc_nec
    merge_gis = args.merge_gis
    if args.hres_flag:
        hires_pft = 'on'
    else: 
        hires_pft = 'off'

    hostname = subprocess.check_output('hostname').strip().decode(encoding='UTF-8')
    print (f"hostname is {hostname}")
    logname = subprocess.check_output('logname').strip().decode(encoding='UTF-8')
    print (f"logname is {logname}")

    # determine pft_years - needed to parse xml file
    if int(start_year) == 1850 and int(end_year) == 1850:
        pft_years = "1850"
    elif int(start_year) == 2000 and int(end_year) == 2000:
        pft_years = "2000"
    elif int(start_year) == 2005 and int(end_year) == 2005:
        pft_years = "2005"
    elif int(start_year) >= 850 and int(end_year) <= 1849:
        pft_years = "0850-1849"
    elif int(start_year) >= 1850 and int(start_year) <= 2100 and int(end_year) <= 2005:
        pft_years = "1850-2015"
    elif int(start_year) >= 2016 and int(start_year) <= 2100 and int(end_year) <=2100:
        pft_years = "2016-2100"
    elif potveg:
        pft_years = "PtVg"
    else:
        print (f"start_year is {start_year} and end_year is {end_year}")
        print (f"ERROR: start and end years should be between 850 and 2105 or pot_veg flag needs to be set")
        sys.exit()

    # error check on glc_nec
    if (glc_nec <= 0) or (glc_nec >= 100):
        raise argparse.ArgumentTypeError("ERROR: glc_nec must be between 1 and 99.")

    # if args.debug:
    #     logging.basicConfig(level=logging.DEBUG)

    # create attribute list for parsing xml file
    attribute_list = {'hires_pft':hires_pft,
                      'pft_years':pft_years,
                      'ssp_rcp':ssp_rcp,
                      'mergeGIS':merge_gis}

    # create dictionary for raw data files names
    rawdata_files = {}

    tree1 = ET.parse('./gen_mksurfdata_namelist.xml')
    root = tree1.getroot()
    root.tag
    root.attrib
    for child1 in root:
        max_match_num = -1
        max_match_child = None
        for child2 in child1:
            if child2.tag == 'entry':
                num_match = 0
                for attrib in attribute_list:
                    # Get the value of the attrib for the entry
                    childval = child2.get(attrib, default=None)
                    if childval == attribute_list[attrib]:
                        num_match += 1
                    elif childval is not None:
                        num_match = -1
                        break
                if num_match > max_match_num: 
                    max_match_num = num_match
                    max_match_child = child2

        if max_match_child is None:
            print (f"{child1.tag} has no matches")
            raise "ERROR"

        for item in max_match_child:
            if item.tag == 'data_filename': 
                rawdata_files[child1.tag] = os.path.join(input_path, item.text)
                if '%y' in rawdata_files[child1.tag]:
                    # ERROR: keep %y here and do the replacement later
                    rawdata_files[child1.tag] = rawdata_files[child1.tag].replace("%y",str(start_year))
                if not os.path.isfile(rawdata_files[child1.tag]):
                    print(f"ERROR: intput rawdata file {rawdata_files[child1.tag]} does not exist")
                    sys.exit()

            if item.tag == 'mesh_filename':
                new_key = f"{child1.tag}_mesh"
                rawdata_files[new_key] = os.path.join(input_path, item.text)
                if not os.path.isfile(rawdata_files[new_key]):
                    print(f"ERROR: mesh file {rawdata_files[new_key]} does not exist")
                    sys.exit()

    tree2 = ET.parse('../../ccs_config/component_grids_nuopc.xml')
    root = tree2.getroot()
    model_mesh = ""
    for child1 in root:  # this is domain tag
        for name, value in child1.attrib.items():
            if value == res:
                for child2 in child1:
                    if child2.tag == 'mesh':
                        model_mesh = child2.text
                        rawdata_files["mksrf_fgrid_mesh"] = os.path.join(input_path,model_mesh.strip('$DIN_LOC_ROOT/'))
                    if child2.tag == 'nx':
                        rawdata_files["mksrf_fgrid_mesh_nx"] = child2.text
                    if child2.tag == 'ny':
                        rawdata_files["mksrf_fgrid_mesh_ny"] = child2.text

    if len(model_mesh) == 0:
        print (f"ERROR: input res {res} is invalid")
        valid_grids = []
        for child1 in root:  # this is domain tag
            for name, value in child1.attrib.items():
                valid_grids.append(value)
        print (f"valid grid values are {valid_grids}")
        sys.exit()

    # Determine num_pft
    if crop_flag:
        num_pft = "78"
    else:
        num_pft = "16"
    logger.debug(f" crop_flag = {str(crop_flag)} => num_pft = {num_pft}")

    # Create land-use txt file for a transient case.
    # Determine the run type and if a transient run create output landuse txt file
    if end_year > start_year:
        run_type = "transient"
    else:
        run_type = "timeslice"
    logger.info(f" run_type  = {run_type}")
    if run_type == 'transient':
        landuse_fname = f"transient_timeseries_hist_{num_pft}pfts_simyr{start_year}-{end_year}.txt"
        with open(landuse_fname, "w", encoding='utf-8') as landuse_file:
            for year in range(start_year, end_year + 1):
                file1 = rawdata_files["mksrf_fvegtyp"]
                landuse_input_fname = file1.replace("%y",str(start_year))
                if not os.path.isfile(landuse_input_fname):
                    logger.warning(f"landunit_input_fname: {landuse_input_fname}")
                landuse_line = f"{landuse_input_fname:<196}{str(year)}\n"
                # -- Each line is written twice in the original pl code:
                landuse_file.write(landuse_line)
                landuse_file.write(landuse_line)
                logger.debug(f"year : {year}")
                logger.debug(landuse_line)
        logger.info("Successfully created land use file : {landuse_fname}.")
        logger.info("-------------------------------------------------------")
    else:
        landuse_fname = ""

    time_stamp = datetime.today().strftime("%y%m%d")
    if end_year == start_year:
        nlfname = f"surfdata_{res}_{ssp_rcp}_{num_pft}pfts_CMIP6_{start_year}_c{time_stamp}.namelist"
        fsurdat = f"surfdata_{res}_{ssp_rcp}_{num_pft}pfts_CMIP6_{start_year}_c{time_stamp}.nc"
        fsurlog = f"surfdata_{res}_{ssp_rcp}_{num_pft}pfts_CMIP6_{start_year}_c{time_stamp}.log"
        fdyndat = ""
    else:
        nlfname = f"surfdata_{res}_{ssp_rcp}_{num_pft}pfts_CMIP6_{start_year}-{end_year}_c{time_stamp}.namelist"
        fsurdat = f"surfdata_{res}_{ssp_rcp}_{num_pft}pfts_CMIP6_{start_year}-{end_year}_c{time_stamp}.nc"
        fsurlog = f"surfdata_{res}_{ssp_rcp}_{num_pft}pfts_CMIP6_{start_year}-{end_year}_c{time_stamp}.log"
        fdyndat = f"landuse.timeseries_{res}_{ssp_rcp}_{num_pft}_CMIP6_{start_year}-{end_year}_c{time_stamp}.log"

    gitdescribe = subprocess.check_output('git describe', shell=True).strip()
    gitdescribe = gitdescribe.decode('utf-8')

    print (f"Creating input namelist file {nlfname}")
    with open(nlfname, "w",encoding='utf-8') as nlfile:
        nlfile.write("&mksurfdata_input \n")

        for key,value in rawdata_files.items():
            if key == 'mksrf_fgrid_mesh_nx' or key == 'mksrf_fgrid_mesh_ny':
                nlfile.write(f"  {key} = {value} \n")
            elif key != "mksrf_fvic" and key != "mksrf_fvic_mesh":
                nlfile.write(f"  {key} = \'{value}\' \n")
        
        mksrf_hrvtyp = rawdata_files["mksrf_fvegtyp"]
        nlfile.write( f"  mksrf_fhrvtyp = \'{mksrf_hrvtyp}\' \n")

        mksrf_hrvtyp_mesh = rawdata_files["mksrf_fvegtyp_mesh"]
        nlfile.write( f"  mksrf_fhrvtyp_mesh = \'{mksrf_hrvtyp_mesh}\' \n")

        nlfile.write(f"  numpft = {num_pft} \n")
        nlfile.write( "  no_inlandwet = .true. \n")
        nlfile.write( "  fdyndat = \' \' \n")
        nlfile.write(f"  fsurdat = \'{fsurdat}\' \n")
        nlfile.write(f"  fsurlog = \'{fsurlog}\' \n")
        nlfile.write(f"  mksrf_fdynuse = \'{landuse_fname} \' \n")
        nlfile.write(f"  gitdescribe = \'{gitdescribe}\' \n")
        nlfile.write( "  outnc_large_files = .false. \n")
        nlfile.write( "  outnc_double = .true. \n")

        if glc_flag:
            nlfile.write( "  outnc_3dglc = .true. \n")
        else:
            nlfile.write( "  outnc_3dglc = .false. \n")

        if vic_flag:
            nlfile.write( "  outnc_vic = .true. \n")
            mksrf_fvic = rawdata_files["mksrf_fvic"]
            nlfile.write(f"  mksrf_fvic = \'{mksrf_fvic}\' \n")
            mksrf_fvic_mesh = rawdata_files["mksrf_fvic_mesh"]
            nlfile.write(f"  mksrf_fvic_mesh = \'{mksrf_fvic_mesh}\' \n")
        else:
            nlfile.write("  outnc_vic = .false. \n")

        nlfile.write(f"  logname = \'{logname}\' \n")
        nlfile.write(f"  hostname = \'{hostname}\' \n")

        nlfile.write("/ \n")

if __name__ == "__main__":
    main()
