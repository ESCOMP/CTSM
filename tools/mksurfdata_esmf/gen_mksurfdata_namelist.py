#!/usr/bin/env python3

import sys, os, shutil
import xml.etree.ElementTree as ET
import logging
import argparse, textwrap
import subprocess
from datetime import datetime

logger = logging.getLogger(__name__)

# valid options for SSP scenarios and pft years:
valid_opts = {"ssp-rcp": ["none","SSP1-2.6","SSP3-7.0","SSP5-3.4","SSP2-4.5","SSP1-1.9","SSP4-3.4","SSP4-6.0","SSP5-8.5"]}

def get_parser():
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help

    parser.add_argument(
        '-v', '--verbose',
        help="ncrease output verbosity",
        action="store_true",
    )
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
        "--model-mesh",
        help="""
            model mesh [default: %(default)s]
            Ignore the --res option and force the model mesh file to be this input
            """,
        action="store",
        dest="force_model_mesh_file",
        required=False,
        default="none",
    )
    parser.add_argument(
        "--model-mesh-nx",
        help="""
            model mesh [default: %(default)s]
            Ignore the --res option and force the model mesh to have this nx 
            """,
        action="store",
        dest="force_model_mesh_nx",
        required=False,
        default="-999",
    )
    parser.add_argument(
        "--model-mesh-ny",
        help="""
            model mesh [default: %(default)s]
            Ignore the --res option and force the model mesh to have this ny 
            """,
        action="store",
        dest="force_model_mesh_ny",
        required=False,
        default="-999",
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

    args = get_parser().parse_args()

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
    glc_flag = args.glc_flag
    potveg = args.potveg_flag
    glc_nec = args.glc_nec
    merge_gis = args.merge_gis
    if args.hres_flag:
        if (start_year == 1850 and end_year == 1850) or \
           (start_year == 2005 and end_year == 2005):
            hires_pft = 'on'
        else:
            print(f"ERROR: for --hires_pft you must set both start-year & end-year to 1850 or to 2005")
            sys.exit(5)
    else:
        hires_pft = 'off'
    verbose = args.verbose

    if force_model_mesh_file != 'none':
        res = force_model_mesh_nx + 'x' + force_model_mesh_ny

    hostname = os.getenv("HOSTNAME")
    logname = os.getenv("LOGNAME")
    if args.verbose:
        print (f"hostname is {hostname}")
        print (f"logname is {logname}")

    if ssp_rcp == 'none':
        if int(start_year) > 2015:
            print(f"ERROR: if start-year is > 2015 must add an --ssp_rcp argument that is not 'none")
            print(f"  valid opts for ssp-rcp are {valid_opts}")
            sys.exit(10)
        elif int(end_year) > 2015:
            print(f"ERROR: if end-year is > 2015 must add an --ssp-rcp argument that is not 'none")
            print(f"  valid opts for ssp-rcp are {valid_opts}")
            sys.exit(10)

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
        print (f"start_year is {start_year} and end_year is {end_year}")
        print (f"ERROR: start and end years should be between 850 and 2105 or pot_veg flag needs to be set")
        sys.exit(10)

    if verbose:
        print (f"pft_years = {pft_years}")

    # Create land-use txt file for a transient case.
    # Determine the run type and if a transient run create output landuse txt file
    if end_year > start_year:
        run_type = "transient"
    else:
        run_type = "timeslice"
    if verbose:
        print(f"run_type  = {run_type}")

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
            # For years greater than 2015 - mksrf_fvegtyp_ssp must have a match
            if start_year <= 2015:
                if 'mksrf_fvegtyp_ssp' not in child1.tag:
                    print (f"ERROR: {child1.tag} has no matches")
                    sys.exit(15)
                else:
                    continue
            else:
                # For years less than 2015 - mksrf_fvegtyp must have a match
                if 'mksrf_fvegtyp' not in child1.tag:
                    print (f"ERROR: {child1.tag} has no matches")
                    sys.exit(15)
                else:
                    continue

        for item in max_match_child:
            if item.tag == 'data_filename':
                rawdata_files[child1.tag] = os.path.join(input_path, item.text)
                if '%y' not in rawdata_files[child1.tag]:
                    if not os.path.isfile(rawdata_files[child1.tag]):
                        print(f"ERROR: input data file {rawdata_files[child1.tag]} for {child1.tag} does not exist")
                        sys.exit(20)

            if item.tag == 'mesh_filename':
                new_key = f"{child1.tag}_mesh"
                rawdata_files[new_key] = os.path.join(input_path, item.text)
                if not os.path.isfile(rawdata_files[new_key]):
                    print(f"ERROR: input mesh file {rawdata_files[new_key]} does not exist")
                    sys.exit(30)

            if item.tag == 'lake_filename':
                new_key = f"{child1.tag}_lake"
                rawdata_files[new_key] = os.path.join(input_path, item.text)

            if item.tag == 'urban_filename':
                new_key = f"{child1.tag}_urban"
                rawdata_files[new_key] = os.path.join(input_path, item.text)

    # determine output mesh
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

    if force_model_mesh_file == 'none' and len(model_mesh) == 0:
        print (f"ERROR: input res {res} is invalid")
        valid_grids = []
        for child1 in root:  # this is domain tag
            for name, value in child1.attrib.items():
                valid_grids.append(value)
        print (f"valid grid values are {valid_grids}")
        sys.exit(40)

    # Determine num_pft
    if nocrop_flag:
        num_pft = "16"
    else:
        num_pft = "78"
    if verbose:
        print (f"num_pft is {num_pft}")

    # Write out if surface dataset will be created
    if verbose:
        if nosurfdata_flag:
            print(f"surface dataset will not be created")
        else:
            print(f"surface dataset will be created")

    if run_type == 'transient':
        if ssp_rcp == 'none':
            landuse_fname = f"landuse_timeseries_hist_{num_pft}pfts_simyr{start_year}-{end_year}.txt"
        else:
            landuse_fname = f"landuse_timeseries_{ssp_rcp}_{num_pft}pfts_CMIP6_simyr{start_year}-{end_year}.txt"

        with open(landuse_fname, "w", encoding='utf-8') as landuse_file:
            for year in range(start_year, end_year + 1):
                if year <= 2015:
                    file1 = rawdata_files["mksrf_fvegtyp"]
                    file2 = rawdata_files["mksrf_fvegtyp_urban"]
                    file3 = rawdata_files["mksrf_fvegtyp_lake"]
                else:
                    file1 = rawdata_files["mksrf_fvegtyp_ssp"]
                    file2 = rawdata_files["mksrf_fvegtyp_ssp_urban"]
                    file3 = rawdata_files["mksrf_fvegtyp_ssp_lake"]

                landuse_input_fname = file1.replace("%y",str(year))
                landuse_input_fnam2 = file2.replace("%y",str(year))
                landuse_input_fnam3 = file3.replace("%y",str(year))
                if not os.path.isfile(landuse_input_fname):
                     print(f"ERROR: landunit_input_fname: {landuse_input_fname} does not exist")
                     sys.exit(60)
                if not os.path.isfile(landuse_input_fnam2):
                     print(f"ERROR: landunit_input_fnam2: {landuse_input_fnam2} does not exist")
                     sys.exit(60)
                if not os.path.isfile(landuse_input_fnam3):
                     print(f"ERROR: landunit_input_fnam3: {landuse_input_fnam3} does not exist")
                     sys.exit(60)

                # -- Each line is written twice in the original perl code:
                landuse_line = f"{landuse_input_fname:<196}{str(year)}\n"
                landuse_lin2 = f"{landuse_input_fnam2:<196}{str(year)}\n"
                landuse_lin3 = f"{landuse_input_fnam3:<196}{str(year)}\n"
                landuse_file.write(landuse_line)
                landuse_file.write(landuse_line)
                landuse_file.write(landuse_lin2)
                landuse_file.write(landuse_lin3)
                logger.debug(f"year : {year}")
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
        nlfname = f"surfdata_{res}_{ssp_rcp_name}_{num_pft}pfts_CMIP6_{start_year}_c{time_stamp}.namelist"
        fsurdat = f"surfdata_{res}_{ssp_rcp_name}_{num_pft}pfts_CMIP6_{start_year}_c{time_stamp}.nc"
        fsurlog = f"surfdata_{res}_{ssp_rcp_name}_{num_pft}pfts_CMIP6_{start_year}_c{time_stamp}.log"
        fdyndat = ""
    else:
        nlfname = f"surfdata_{res}_{ssp_rcp_name}_{num_pft}pfts_CMIP6_{start_year}-{end_year}_c{time_stamp}.namelist"
        fsurdat = f"surfdata_{res}_{ssp_rcp_name}_{num_pft}pfts_CMIP6_{start_year}-{end_year}_c{time_stamp}.nc"
        fsurlog = f"surfdata_{res}_{ssp_rcp_name}_{num_pft}pfts_CMIP6_{start_year}-{end_year}_c{time_stamp}.log"
        fdyndat = f"landuse.timeseries_{res}_{ssp_rcp_name}_{num_pft}_CMIP6_{start_year}-{end_year}_c{time_stamp}.nc"

    gitdescribe = subprocess.check_output('git describe', shell=True).strip()
    gitdescribe = gitdescribe.decode('utf-8')

    # The below two overrides are only used for testing an validation
    # it takes a long time to generate the mapping files
    # from 1km to the following two resolutions since the output mesh has so few points
    if res == "10x15":
        mksrf_ftopostats_override = os.path.join(input_path,"lnd","clm2","rawdata","surfdata_topo_10x15_c220303.nc")
        if args.verbose:
            print (f"will override mksrf_ftopostats with = {mksrf_ftopostats_override}")
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
        if force_model_mesh_file != 'none':
            mksrf_fgrid_mesh_nx = force_model_mesh_nx
            mksrf_fgrid_mesh_ny = force_model_mesh_ny
            mksrf_fgrid_mesh    = force_model_mesh_file
        else:
            mksrf_fgrid_mesh_nx = rawdata_files["mksrf_fgrid_mesh_nx"] 
            mksrf_fgrid_mesh_ny = rawdata_files["mksrf_fgrid_mesh_ny"] 
            mksrf_fgrid_mesh    = rawdata_files["mksrf_fgrid_mesh"] 
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
            print(f"ERROR: input mksrf_fvegtyp file {mksrf_fvegtyp} does not exist")
            sys.exit(20)
        if not os.path.isfile(mksrf_fhrvtyp):
            print(f"ERROR: input mksrf_fhrvtyp file {mksrf_fhrvtyp} does not exist")
            sys.exit(20)
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
            nlfile.write(f"  fsurdat = \' \' \n")
        else:
            nlfile.write(f"  fsurdat = \'{fsurdat}'\n")
        nlfile.write(f"  fsurlog = \'{fsurlog}\' \n")
        nlfile.write(f"  fdyndat = \'{fdyndat}\' \n")

        # -------------------
        # output data logicals
        # -------------------
        nlfile.write(f"  numpft = {num_pft} \n")
        nlfile.write( "  no_inlandwet = .true. \n")
        if glc_flag:
            nlfile.write( "  outnc_3dglc = .true. \n")
        else:
            nlfile.write( "  outnc_3dglc = .false. \n")
        if vic_flag:
            nlfile.write( "  outnc_vic = .true. \n")
        else:
            nlfile.write("  outnc_vic = .false. \n")
        nlfile.write( "  outnc_large_files = .false. \n")
        nlfile.write( "  outnc_double = .true. \n")
        nlfile.write(f"  logname = \'{logname}\' \n")
        nlfile.write(f"  hostname = \'{hostname}\' \n")
        nlfile.write(f"  gitdescribe = \'{gitdescribe}\' \n")

        nlfile.write("/ \n")

    print (f"Successfully created input namelist file {nlfname}")
    sys.exit(0)

if __name__ == "__main__":
    main()
