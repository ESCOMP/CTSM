#!/usr/bin/env python3
"""
gen_mksurfdata_namelist.py generates a namelist for use with the mksurfdata
executable. For detailed instructions, see README.
"""
import os
import sys
import xml.etree.ElementTree as ET
import logging
import argparse
import textwrap
import subprocess
from datetime import datetime
import netCDF4

_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            os.pardir,
                            os.pardir,
                            'python')
sys.path.insert(1, _CTSM_PYTHON)

from ctsm.ctsm_logging import setup_logging_pre_config, add_logging_args, process_logging_args

logger = logging.getLogger(__name__)

# valid options for SSP/RCP scenarios
valid_opts = {'ssp-rcp': ['SSP1-2.6', 'SSP3-7.0', 'SSP5-3.4', 'SSP2-4.5',
                          'SSP1-1.9', 'SSP4-3.4', 'SSP4-6.0', 'SSP5-8.5',
                          'none']}

def get_parser():
    """
    Get parser object for this script.
    """
    # set up logging allowing user control
    setup_logging_pre_config()

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help
    add_logging_args(parser)

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
            Model resolution (required) 
            To see available supported resolutions, simply invoke this command
            with a --res unknown option. For custom resolutions, provide a grid
            name of your choosing to be used in the name of the fsurdat file.
            """,
        action="store",
        dest="res",
        required=True,
    )
    parser.add_argument(
        "--model-mesh",
        help="""
            model mesh [default: %(default)s]
            Ignore --res and use --model-mesh to be this file
            """,
        action="store",
        dest="force_model_mesh_file",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--model-mesh-nx",
        help="""
            model mesh [default: %(default)s]
            Required when using --model-mesh: set nx to the grid's number of
            columns; expect nx x ny = elementCount for consistency with the
            model mesh
            """,
        action="store",
        dest="force_model_mesh_nx",
        required=False,
        default=None
    )
    parser.add_argument(
        "--model-mesh-ny",
        help="""
            model mesh [default: %(default)s]
            Required when using --model-mesh: set ny to the grid's number of
            rows; expect nx x ny = elementCount for consistency with the model
            mesh
            """,
        action="store",
        dest="force_model_mesh_ny",
        required=False,
        default=None
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
        choices=valid_opts["ssp-rcp"],
        default="none",
    )
    parser.add_argument(
        "--rawdata-dir",
        help="""
            /path/of/root/of/input/data
            on izumi use /fs/cgd/csm/inputdata
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
        "--inlandwet",
        help="""
            Flag for including inland wetlands.
            [default: %(default)s]
            """,
        action="store_true",
        dest="inlandwet",
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
        "--nosurfdata",
        help="""
            Do not output a surface datase
            This is useful if you only want a landuse_timeseries file
            [default: %(default)s]
            """,
        action="store_true",
        dest="surfdata_flag",
        default=False,
    )
    parser.add_argument(
        "--nocrop",
        help="""
            Do not create datasets with the extensive list of prognostic crop types.
            [default: %(default)s]
            """,
        action="store_true",
        dest="crop_flag",
        default=False,
    )
    parser.add_argument(
        "--potveg_flag",
        help="""
            Use Potential Vegetation for pft_years
            [default: %(default)s]
            """,
        action="store_true",
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
    """
    See docstring at the top.
    """
    args = get_parser().parse_args()
    process_logging_args(args)

    start_year = args.start_year
    end_year = args.end_year
    ssp_rcp = args.ssp_rcp
    res = args.res
    force_model_mesh_file = args.force_model_mesh_file
    force_model_mesh_nx = args.force_model_mesh_nx
    force_model_mesh_ny = args.force_model_mesh_ny
    input_path = args.input_path
    nocrop_flag = args.crop_flag
    nosurfdata_flag = args.surfdata_flag
    vic_flag = args.vic_flag
    inlandwet = args.inlandwet
    glc_flag = args.glc_flag
    potveg = args.potveg_flag
    glc_nec = args.glc_nec
    merge_gis = args.merge_gis
    if args.hres_flag:
        if (start_year == 1850 and end_year == 1850) or \
           (start_year == 2005 and end_year == 2005):
            hires_pft = 'on'
        else:
            error_msg = 'ERROR: for --hires_pft you must set both start-year ' \
                        'and end-year to 1850 or to 2005'
            sys.exit(error_msg)
    else:
        hires_pft = 'off'

    if force_model_mesh_file is not None:
        # open mesh_file to read element_count and, if available, orig_grid_dims
        mesh_file = netCDF4.Dataset(force_model_mesh_file, 'r')
        element_count = mesh_file.dimensions['elementCount'].size
        if 'origGridDims' in mesh_file.variables:
            orig_grid_dims = mesh_file.variables['origGridDims']
            force_model_mesh_nx = orig_grid_dims[0]
            force_model_mesh_ny = orig_grid_dims[1]
            mesh_file.close()
            important_msg = 'Found data for force_model_mesh_nx and ' \
                            'force_model_mesh_ny in the mesh file so ' \
                            'IGNORING ANY CORRESPONDING USER-ENTERED VALUES.'
            logger.info(important_msg)
        elif force_model_mesh_nx is None or force_model_mesh_ny is None:
            error_msg = 'ERROR: You set --model-mesh but the file does not ' \
                        'contain the variable origGridDims, so you MUST ALSO ' \
                        'SET --model-mesh-nx AND --model-mesh-ny'
            sys.exit(error_msg)

        # using force_model_mesh_nx and force_model_mesh_ny either from the
        # mesh file (see previous if statement) or the user-entered values
        if element_count != int(force_model_mesh_nx) * int(force_model_mesh_ny):
            error_msg = 'ERROR: The product of ' \
                        '--model-mesh-nx x --model-mesh-ny must equal ' \
                        'exactly elementCount in --model-mesh'
            sys.exit(error_msg)

    hostname = os.getenv("HOSTNAME")
    logname = os.getenv("LOGNAME")

    logger.info('hostname is %s', hostname)
    logger.info('logname is %s', logname)

    if ssp_rcp == 'none':
        if int(start_year) > 2015:
            error_msg = 'ERROR: if start-year > 2015 must add an --ssp_rcp ' \
                        'argument that is not none: valid opts for ssp-rcp ' \
                        f'are {valid_opts}'
            sys.exit(error_msg)
        elif int(end_year) > 2015:
            error_msg = 'ERROR: if end-year > 2015 must add an --ssp-rcp ' \
                        'argument that is not none: valid opts for ssp-rcp ' \
                        f'are {valid_opts}'
            sys.exit(error_msg)

    pft_years_ssp = "-999"

    # determine pft_years - needed to parse xml file
    if potveg:
        pft_years = "PtVg"
    elif int(start_year) == 1850 and int(end_year) == 1850:
        pft_years = "1850"
    elif int(start_year) == 2000 and int(end_year) == 2000:
        pft_years = "2000"
    elif int(start_year) == 2005 and int(end_year) == 2005:
        pft_years = "2005"
    elif int(start_year) >= 850 and int(end_year) <= 1849:
        pft_years = "0850-1849"
    elif int(start_year) >= 1850 and int(start_year) <= 2100 and int(end_year) <= 2015:
        pft_years = "1850-2015"
    elif int(start_year) >= 1850 and int(start_year) <= 2100 and int(end_year) <= 2100:
        pft_years = "1850-2015"
        pft_years_ssp = "2016-2100"
    elif int(start_year) >= 2016 and int(start_year) <= 2100 and int(end_year) <=2100:
        pft_years = "-999"
        pft_years_ssp = "2016-2100"
    else:
        error_msg = f'ERROR: start_year is {start_year} and end_year is ' \
                    f'{end_year}; start/end years should be between 850 and ' \
                     ' 2105 or pot_veg flag must be set'
        sys.exit(error_msg)

    logger.info('pft_years = %s', pft_years)

    # Create land-use txt file for a transient case.
    # Determine the run type and if a transient run create output landuse txt file
    if end_year > start_year:
        run_type = "transient"
    else:
        run_type = "timeslice"
    logger.info('run_type  = %s', run_type)

    # error check on glc_nec
    if (glc_nec <= 0) or (glc_nec >= 100):
        raise argparse.ArgumentTypeError("ERROR: glc_nec must be between 1 and 99.")

    # create attribute list for parsing xml file
    attribute_list = {'hires_pft':hires_pft,
                      'pft_years':pft_years,
                      'pft_years_ssp':pft_years_ssp,
                      'ssp_rcp':ssp_rcp,
                      'mergeGIS':merge_gis,
                      'res':res}

    # create dictionary for raw data files names
    rawdata_files = {}

    # determine input rawdata
    _must_run_download_input_data = False
    tool_path = os.path.dirname(os.path.abspath(__file__))
    xml_path = os.path.join(tool_path, 'gen_mksurfdata_namelist.xml')
    tree1 = ET.parse(xml_path)
    root = tree1.getroot()
    logger.info('root.tag: %s', root.tag)
    logger.info('root.attrib: %s', root.attrib)
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
            # For years greater than 2015 - mksrf_fvegtyp_ssp must have a match
            if start_year <= 2015:
                if 'mksrf_fvegtyp_ssp' not in child1.tag:
                    error_msg = f'ERROR: {child1.tag} has no matches'
                    sys.exit(error_msg)
                else:
                    continue
            else:
                # For years less than 2015 - mksrf_fvegtyp must have a match
                if 'mksrf_fvegtyp' not in child1.tag:
                    error_msg = f'ERROR: {child1.tag} has no matches'
                    sys.exit(error_msg)
                else:
                    continue

        for item in max_match_child:
            if item.tag == 'data_filename':
                rawdata_files[child1.tag] = os.path.join(input_path, item.text)
                if '%y' not in rawdata_files[child1.tag]:
                    if not os.path.isfile(rawdata_files[child1.tag]):
                        print('WARNING: input data file ' \
                              f'{rawdata_files[child1.tag]} for {child1.tag} ' \
                              'does not exist')
                        print('WARNING: run ./download_input_data to try TO ' \
                              'OBTAIN MISSING FILES')
                        _must_run_download_input_data = True

            if item.tag == 'mesh_filename':
                new_key = f"{child1.tag}_mesh"
                rawdata_files[new_key] = os.path.join(input_path, item.text)
                if not os.path.isfile(rawdata_files[new_key]):
                    print('WARNING: input mesh file ' \
                          f'{rawdata_files[new_key]} does not exist')
                    print('WARNING: run ./download_input_data to try TO ' \
                          'OBTAIN MISSING FILES')
                    _must_run_download_input_data = True

            if item.tag == 'lake_filename':
                new_key = f"{child1.tag}_lake"
                rawdata_files[new_key] = os.path.join(input_path, item.text)

            if item.tag == 'urban_filename':
                new_key = f"{child1.tag}_urban"
                rawdata_files[new_key] = os.path.join(input_path, item.text)

    # determine output mesh
    xml_path = os.path.join(tool_path, '../../ccs_config/component_grids_nuopc.xml')
    tree2 = ET.parse(xml_path)
    root = tree2.getroot()
    model_mesh = ""
    for child1 in root:  # this is domain tag
        for _, value in child1.attrib.items():
            if value == res:
                for child2 in child1:
                    if child2.tag == 'mesh':
                        model_mesh = child2.text
                        rawdata_files['mksrf_fgrid_mesh'] = \
                            os.path.join(input_path,
                                         model_mesh.strip('$DIN_LOC_ROOT/'))
                    if child2.tag == 'nx':
                        rawdata_files["mksrf_fgrid_mesh_nx"] = child2.text
                    if child2.tag == 'ny':
                        rawdata_files["mksrf_fgrid_mesh_ny"] = child2.text

    if not model_mesh and force_model_mesh_file is None:
        valid_grids = []
        for child1 in root:  # this is domain tag
            for _, value in child1.attrib.items():
                valid_grids.append(value)
        if res in valid_grids:
            error_msg = 'ERROR: You have requested a valid grid for which ' \
        '../../ccs_config/component_grids_nuopc.xml does not include a mesh ' \
        'file. For a regular regional or 1x1 grid, you may generate the ' \
        'fsurdat file using the subset_data tool instead. Alternatively ' \
        'and definitely for curvilinear grids, you may generate ' \
        'a mesh file using the workflow currently (2022/7/6) described in ' \
        'https://github.com/ESCOMP/CTSM/issues/1773#issuecomment-1163432584'
            sys.exit(error_msg)
        else:
            error_msg = f'ERROR: invalid input res {res}; ' \
                        f'valid grid values are {valid_grids}'
            sys.exit(error_msg)

    # Determine num_pft
    if nocrop_flag:
        num_pft = "16"
    else:
        num_pft = "78"
    logger.info('num_pft is %s', num_pft)

    # Write out if surface dataset will be created
    if nosurfdata_flag:
        logger.info('surface dataset will not be created')
    else:
        logger.info('surface dataset will be created')

    if run_type == 'transient':
        if ssp_rcp == 'none':
            landuse_fname = \
 f"landuse_timeseries_hist_{num_pft}pfts_simyr{start_year}-{end_year}.txt"
        else:
            landuse_fname = \
 f"landuse_timeseries_{ssp_rcp}_{num_pft}pfts_CMIP6_simyr{start_year}-{end_year}.txt"

        with open(landuse_fname, "w", encoding='utf-8') as landuse_file:
            for year in range(start_year, end_year + 1):
                year_str = str(year)
                if year <= 2015:
                    file1 = rawdata_files["mksrf_fvegtyp"]
                    file2 = rawdata_files["mksrf_fvegtyp_urban"]
                    file3 = rawdata_files["mksrf_fvegtyp_lake"]
                else:
                    file1 = rawdata_files["mksrf_fvegtyp_ssp"]
                    file2 = rawdata_files["mksrf_fvegtyp_ssp_urban"]
                    file3 = rawdata_files["mksrf_fvegtyp_ssp_lake"]

                landuse_input_fname = file1.replace("%y", year_str)
                landuse_input_fnam2 = file2.replace("%y", year_str)
                landuse_input_fnam3 = file3.replace("%y", year_str)
                if not os.path.isfile(landuse_input_fname):
                    print('WARNING: landunit_input_fname: ' \
                          f'{landuse_input_fname} does not exist')
                    print('WARNING: run ./download_input_data to try TO ' \
                          'OBTAIN MISSING FILES')
                    _must_run_download_input_data = True
                if not os.path.isfile(landuse_input_fnam2):
                    print('WARNING: landunit_input_fnam2: ' \
                          f'{landuse_input_fnam2} does not exist')
                    print('WARNING: run ./download_input_data to try TO ' \
                          'OBTAIN MISSING FILES')
                    _must_run_download_input_data = True
                if not os.path.isfile(landuse_input_fnam3):
                    print('WARNING: landunit_input_fnam3: ' \
                          f'{landuse_input_fnam3} does not exist')
                    print('WARNING: run ./download_input_data to try TO ' \
                          'OBTAIN MISSING FILES')
                    _must_run_download_input_data = True

                # -- Each line is written twice in the original perl code:
                landuse_line = f"{landuse_input_fname:<196}{year_str}\n"
                landuse_lin2 = f"{landuse_input_fnam2:<196}{year_str}\n"
                landuse_lin3 = f"{landuse_input_fnam3:<196}{year_str}\n"
                landuse_file.write(landuse_line)
                landuse_file.write(landuse_line)
                landuse_file.write(landuse_lin2)
                landuse_file.write(landuse_lin3)
                logger.debug('year : %s', year_str)
                logger.debug(landuse_line)
        print(f"Successfully created input landuse file {landuse_fname}")
    else:
        landuse_fname = ""

    time_stamp = datetime.today().strftime("%y%m%d")
    if ssp_rcp == 'none':
        ssp_rcp_name = 'hist'
    else:
        ssp_rcp_name = ssp_rcp
    if int(end_year) == int(start_year):
        fdyndat = ''
        prefix = \
 f'surfdata_{res}_{ssp_rcp_name}_{num_pft}pfts_CMIP6_{start_year}_c{time_stamp}.'
    else:
        fdyndat = \
 f'landuse.timeseries_{res}_{ssp_rcp_name}_{num_pft}_CMIP6_{start_year}-{end_year}_c{time_stamp}.nc'
        prefix = \
 f'surfdata_{res}_{ssp_rcp_name}_{num_pft}pfts_CMIP6_{start_year}-{end_year}_c{time_stamp}.'

    nlfname = f'{prefix}namelist'
    fsurdat = f'{prefix}nc'
    fsurlog = f'{prefix}log'

    git_desc_cmd = f'git -C {tool_path} describe'
    try:
        # The "git -C" option permits a system test to run this tool from
        # elsewhere while running the git command from the tool_path
        gitdescribe = subprocess.check_output(git_desc_cmd, shell=True).strip()
    except subprocess.CalledProcessError as e:
        # In case the "git -C" option is unavailable, as on casper (2022/5/24)
        # Still, this does NOT allow the system test to work on machines
        # without git -C
        logger.info('git -C option unavailable on casper as of 2022/7/2 %s', e)
        gitdescribe = subprocess.check_output('git describe', shell=True).strip()
    gitdescribe = gitdescribe.decode('utf-8')

    # The below two overrides are only used for testing an validation
    # it takes a long time to generate the mapping files
    # from 1km to the following two resolutions since the output mesh has so few points
    if res == "10x15":
        mksrf_ftopostats_override = os.path.join(input_path, 'lnd', 'clm2',
            'rawdata', 'surfdata_topo_10x15_c220303.nc')
        logger.info('will override mksrf_ftopostats with = %s',
                    mksrf_ftopostats_override)
    else:
        mksrf_ftopostats_override = ""

    # ----------------------------------------
    # Write output namelist file
    # ----------------------------------------

    with open(nlfname, "w",encoding='utf-8') as nlfile:
        nlfile.write("&mksurfdata_input \n")

        # -------------------
        # raw input data
        # -------------------
        if force_model_mesh_file is None:
            mksrf_fgrid_mesh_nx = rawdata_files["mksrf_fgrid_mesh_nx"]
            mksrf_fgrid_mesh_ny = rawdata_files["mksrf_fgrid_mesh_ny"]
            mksrf_fgrid_mesh    = rawdata_files["mksrf_fgrid_mesh"]
        else:
            mksrf_fgrid_mesh_nx = force_model_mesh_nx
            mksrf_fgrid_mesh_ny = force_model_mesh_ny
            mksrf_fgrid_mesh    = force_model_mesh_file
        nlfile.write( f"  mksrf_fgrid_mesh = \'{mksrf_fgrid_mesh}\' \n")
        nlfile.write( f"  mksrf_fgrid_mesh_nx = {mksrf_fgrid_mesh_nx} \n")
        nlfile.write( f"  mksrf_fgrid_mesh_ny = {mksrf_fgrid_mesh_ny} \n")

        for key,value in rawdata_files.items():
            if key == 'mksrf_ftopostats' and mksrf_ftopostats_override != '':
                nlfile.write(f"  mksrf_ftopostats_override = \'{mksrf_ftopostats_override}\' \n")
            elif '_fvic' not in key and 'mksrf_fvegtyp' not in key and 'mksrf_fgrid' not in key:
                # write everything else
                nlfile.write(f"  {key} = \'{value}\' \n")

        if start_year <= 2015:
            mksrf_fvegtyp       = rawdata_files["mksrf_fvegtyp"]
            mksrf_fvegtyp_mesh  = rawdata_files["mksrf_fvegtyp_mesh"]
            mksrf_fhrvtyp       = rawdata_files["mksrf_fvegtyp"]
            mksrf_fhrvtyp_mesh  = rawdata_files["mksrf_fvegtyp_mesh"]
        else:
            mksrf_fvegtyp       = rawdata_files["mksrf_fvegtyp_ssp"]
            mksrf_fvegtyp_mesh  = rawdata_files["mksrf_fvegtyp_ssp_mesh"]
            mksrf_fhrvtyp       = rawdata_files["mksrf_fvegtyp_ssp"]
            mksrf_fhrvtyp_mesh  = rawdata_files["mksrf_fvegtyp_ssp_mesh"]
        if '%y' in mksrf_fvegtyp:
            mksrf_fvegtyp = mksrf_fvegtyp.replace("%y",str(start_year))
        if '%y' in mksrf_fhrvtyp:
            mksrf_fhrvtyp = mksrf_fhrvtyp.replace("%y",str(start_year))
        if not os.path.isfile(mksrf_fvegtyp):
            print('WARNING: input mksrf_fvegtyp file ' \
                  f'{mksrf_fvegtyp} does not exist')
            print('WARNING: run ./download_input_data to try TO ' \
                  'OBTAIN MISSING FILES')
            _must_run_download_input_data = True
        if not os.path.isfile(mksrf_fhrvtyp):
            print('WARNING: input mksrf_fhrvtyp file ' \
                  f'{mksrf_fhrvtyp} does not exist')
            print('WARNING: run ./download_input_data to try TO ' \
                  'OBTAIN MISSING FILES')
            _must_run_download_input_data = True
        nlfile.write( f"  mksrf_fvegtyp = \'{mksrf_fvegtyp}\' \n")
        nlfile.write( f"  mksrf_fvegtyp_mesh = \'{mksrf_fvegtyp_mesh}\' \n")
        nlfile.write( f"  mksrf_fhrvtyp = \'{mksrf_fhrvtyp}\' \n")
        nlfile.write( f"  mksrf_fhrvtyp_mesh = \'{mksrf_fhrvtyp_mesh}\' \n")

        if vic_flag:
            mksrf_fvic = rawdata_files["mksrf_fvic"]
            nlfile.write(f"  mksrf_fvic = \'{mksrf_fvic}\' \n")
            mksrf_fvic_mesh = rawdata_files["mksrf_fvic_mesh"]
            nlfile.write(f"  mksrf_fvic_mesh = \'{mksrf_fvic_mesh}\' \n")

        nlfile.write( f"  mksrf_fdynuse = \'{landuse_fname} \' \n")

        # -------------------
        # output data files
        # -------------------
        if nosurfdata_flag:
            nlfile.write("  fsurdat = \' \' \n")
        else:
            nlfile.write(f"  fsurdat = \'{fsurdat}'\n")
        nlfile.write(f"  fsurlog = \'{fsurlog}\' \n")
        nlfile.write(f"  fdyndat = \'{fdyndat}\' \n")

        # -------------------
        # output data logicals
        # -------------------
        nlfile.write(f"  numpft = {num_pft} \n")
        nlfile.write(f"  no_inlandwet = .{str(not inlandwet).lower()}. \n")
        nlfile.write(f"  outnc_3dglc = .{str(glc_flag).lower()}. \n")
        nlfile.write(f"  outnc_vic = .{str(vic_flag).lower()}. \n")
        nlfile.write( "  outnc_large_files = .false. \n")
        nlfile.write( "  outnc_double = .true. \n")
        nlfile.write(f"  logname = \'{logname}\' \n")
        nlfile.write(f"  hostname = \'{hostname}\' \n")
        nlfile.write(f"  gitdescribe = \'{gitdescribe}\' \n")

        nlfile.write("/ \n")

    if _must_run_download_input_data:
        temp_nlfname = 'surfdata.namelist'
        os.rename(nlfname, temp_nlfname)
        nlfname = temp_nlfname

    print (f"Successfully created input namelist file {nlfname}")
    sys.exit(0)

if __name__ == "__main__":
    main()
