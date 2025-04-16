"""
This module includes the definition for a RegionalCase classs.
"""

# -- Import libraries
# -- Import Python Standard Libraries
import logging
import os
import argparse
from datetime import datetime

# -- 3rd party libraries
import numpy as np
import xarray as xr

# -- import local classes for this script
from ctsm.site_and_regional.base_case import BaseCase, USRDAT_DIR
from ctsm.site_and_regional.mesh_type import MeshType
from ctsm.utils import add_tag_to_filename
from ctsm.utils import abort
from ctsm.config_utils import check_lon1_lt_lon2

logger = logging.getLogger(__name__)


class RegionalCase(BaseCase):
    """
    A class to encapsulate regional cases.

    ...
    Attributes
    ----------
    lat1 : float
        first (left) latitude of a region.
    lat1 : float
        second (right) latitude of a region.
    lon1 : float
        first (bottom) longitude of a region.
    lon2 : float
        second (top) longitude of a region.
    lon_type : int
        180 if longitudes are in [-180, 180], 360 if they're in [0, 360]
    reg_name: str -- default = None
        Region's name
    create_domain : bool
        flag for creating domain file
    create_mesh : bool
        flag for creating mesh file
    create_surfdata : bool
        flag for creating surface dataset
    create_landuse : bool
        flag for creating landuse file
    create_datm : bool
        flag for creating DATM files
    create_user_mods : bool
        flag for creating user mods files and folders
    overwrite : bool
        flag for over-writing files if they already exist


    Methods
    -------
    create_tag
        Create a tag for this region which is either
        region's name or a combination of bounds of this
        region lat1-lat2_lon1-lon2

    check_region_bounds
        Check for the regional bounds

    check_region_lats
        Check for the regional lats

    create_domain_at_reg
        Create domain file at this region

    create_surfdata_at_reg
        Create surface dataset at this region

    extract_mesh_at_reg
        Extract mesh from the domain dataset created by create_domain_at_reg

    create_landuse_at_reg
        Create landuse file at this region

    write_shell_commands(namelist)
        write out xml commands to a file for usermods (i.e. shell_commands) for regional settings.
    """

    # pylint: disable=too-many-instance-attributes
    # the ones we have are useful

    def __init__(
        self,
        *,
        lat1,
        lat2,
        lon1,
        lon2,
        reg_name,
        create_domain,
        create_surfdata,
        create_landuse,
        create_datm,
        create_user_mods,
        create_mesh,
        out_dir,
        overwrite,
    ):
        """
        Initializes RegionalCase with the given arguments.
        """
        super().__init__(
            create_domain=create_domain,
            create_surfdata=create_surfdata,
            create_landuse=create_landuse,
            create_datm=create_datm,
            create_user_mods=create_user_mods,
            overwrite=overwrite,
        )

        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        self.reg_name = reg_name
        self.create_mesh = create_mesh
        self.mesh = None
        self.out_dir = out_dir
        self.check_region_bounds()
        self.create_tag()
        self.ni = None
        self.nj = None

    def create_tag(self):
        """
        Create a tag for a region which is either the region name
        or
        the lat1-lat2_lon1-lon2 if the region name does not exist.
        """
        if self.reg_name:
            self.tag = self.reg_name
        else:
            self.tag = "{}-{}_{}-{}".format(
                str(self.lon1), str(self.lon2), str(self.lat1), str(self.lat2)
            )

    def check_region_bounds(self):
        """
        Check for the regional bounds
        """
        # If you're calling this, lat/lon bounds need to have been provided
        if any(x is None for x in [self.lon1, self.lon2, self.lat1, self.lat2]):
            raise argparse.ArgumentTypeError(
                "Latitude and longitude bounds must be provided and not None.\n"
                + f"   lon1: {self.lon1}\n"
                + f"   lon2: {self.lon2}\n"
                + f"   lat1: {self.lat1}\n"
                + f"   lat2: {self.lat2}"
            )
        # By now, you need to have already converted to longitude [0, 360]
        check_lon1_lt_lon2(self.lon1, self.lon2, 360)
        self.check_region_lats()

    def check_region_lats(self):
        """
        Check for the regional lat bound
        """
        if self.lat1 >= self.lat2:
            err_msg = """
            \n
            ERROR: lat1 is bigger than lat2.
            lat1 points to the westernmost longitude of the region. {}
            lat2 points to the easternmost longitude of the region. {}
            Please make sure lat1 is smaller than lat2.

            """.format(
                self.lat1, self.lat2
            )
            raise argparse.ArgumentTypeError(err_msg)

    def create_domain_at_reg(self, indir, file):
        """
        Create domain file for this RegionalCase class.
        """

        # specify files
        fdomain_in = os.path.join(indir, file)
        fdomain_out = add_tag_to_filename(fdomain_in, self.tag)
        logger.info("fdomain_in:  %s", fdomain_in)
        logger.info("fdomain_out: %s", os.path.join(self.out_dir, fdomain_out))
        logger.info("Creating domain file at region: %s", self.tag)

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fdomain_in, "xc", "yc", "ni", "nj")
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(nj=yind, ni=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fdomain_in

        # mode 'w' overwrites file
        wfile = os.path.join(self.out_dir, fdomain_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("Successfully created file (fdomain_out) %s", wfile)
        f_in.close()
        f_out.close()

        if self.create_mesh:
            mesh_out = os.path.join(
                self.out_dir,
                os.path.splitext(fdomain_out)[0] + "_ESMF_UNSTRUCTURED_MESH.nc",
            )
            self.mesh = mesh_out
            logger.info("creating mesh file from domain file: %s", wfile)
            ds = xr.open_dataset(wfile, mask_and_scale=False, decode_times=False).transpose()
            self.extract_mesh_at_reg(ds)

    def create_surfdata_at_reg(self, indir, file, user_mods_dir, specify_fsurf_out):
        """
        Create surface data file for this RegionalCase class.
        """

        logger.info("Creating surface dataset file at region: %s", self.tag)

        # specify files
        fsurf_in = os.path.join(indir, file)
        if specify_fsurf_out is None:
            fsurf_out = add_tag_to_filename(fsurf_in, self.tag, replace_res=True)
        else:
            fsurf_out = specify_fsurf_out

        logger.info("fsurf_in:  %s", fsurf_in)
        logger.info("fsurf_out: %s", os.path.join(self.out_dir, fsurf_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fsurf_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fsurf_in

        # mode 'w' overwrites file
        wfile = os.path.join(self.out_dir, fsurf_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("created file (fsurf_out) %s", wfile)
        f_in.close()
        f_out.close()

        # write to user_nl_clm if specified
        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                line = "fsurdat = '${}'".format(os.path.join(USRDAT_DIR, fsurf_out))
                self.write_to_file(line, nl_clm)

    def extract_mesh_at_reg(self, ds_in):
        """
        Create Mesh from Surface dataset netcdf file.
        """
        logger.info("Creating meshfile for  at region: %s", self.tag)

        lat_name = "yc"
        lon_name = "xc"

        lats = ds_in[lat_name].astype(np.float32)
        lons = ds_in[lon_name].astype(np.float32)

        self.ni = len(lats.ni)
        self.nj = len(lats.nj)

        mask = ds_in["mask"].astype(np.float32)
        if mask.max() > 1.0 or mask.min() < 0.0:
            abort("Mask variable is not within 0 to 1")

        this_mesh = MeshType(lats, lons, mask=mask)
        this_mesh.calculate_corners()
        this_mesh.calculate_nodes()
        this_mesh.create_esmf(self.mesh)

    def create_landuse_at_reg(self, indir, file, user_mods_dir):
        """
        Create land use data file for this RegionalCase class.
        """

        logger.info("Creating landuse file at region: %s", self.tag)

        # specify files
        fluse_in = os.path.join(indir, file)
        fluse_out = add_tag_to_filename(fluse_in, self.tag, replace_res=True)
        logger.info("fluse_in:  %s", fluse_in)
        logger.info("fluse_out: %s", os.path.join(self.out_dir, fluse_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fluse_in

        # mode 'w' overwrites file
        wfile = os.path.join(self.out_dir, fluse_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("Successfully created file (fluse_out) %s", wfile)
        f_in.close()
        f_out.close()

        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                # line = "landuse = '${}'".format(os.path.join(USRDAT_DIR, fluse_out))
                line = "flanduse_timeseries = '${}'".format(os.path.join(USRDAT_DIR, fluse_out))
                self.write_to_file(line, nl_clm)

    def create_mesh_at_reg(self, mesh_dir, mesh_surf):
        """
        Create a mesh subsetted for the RegionalCase class.
        """
        logger.info("----------------------------------------------------------------------")
        logger.info("Subsetting mesh file for region: %s", self.tag)

        today = datetime.today()
        today_string = today.strftime("%y%m%d")

        mesh_in = os.path.join(mesh_dir, mesh_surf)
        mesh_out = os.path.join(
            self.out_dir,
            os.path.splitext(mesh_surf)[0] + "_" + self.tag + "_c" + today_string + ".nc",
        )

        logger.info("mesh_in  :  %s", mesh_in)
        logger.info("mesh_out :  %s", mesh_out)

        self.mesh = mesh_out

        _, subset_element, subset_node, conn_dict = self.subset_mesh_at_reg(mesh_in)

        f_in = xr.open_dataset(mesh_in)
        self.write_mesh(f_in, subset_element, subset_node, conn_dict, mesh_out)

    def subset_mesh_at_reg(self, mesh_in):
        """
        This function subsets the mesh based on lat and lon bounds given by RegionalCase class.
        """
        f_in = xr.open_dataset(mesh_in)
        elem_count = len(f_in["elementCount"])
        elem_conn = f_in["elementConn"]
        num_elem_conn = f_in["numElementConn"]
        node_count = len(f_in["nodeCount"])
        node_coords = f_in["nodeCoords"]

        subset_element = []
        cnt = 0

        for n in range(elem_count):
            endx = elem_conn[n, : num_elem_conn[n].values].values
            endx[:,] -= 1  # convert to zero based index
            endx = [int(xi) for xi in endx]

            nlon = node_coords[endx, 0].values
            nlat = node_coords[endx, 1].values

            l1 = np.logical_or(nlon <= self.lon1, nlon >= self.lon2)
            l2 = np.logical_or(nlat <= self.lat1, nlat >= self.lat2)

            if np.any(np.logical_or(l1, l2)):
                pass
            else:
                subset_element.append(n)
                cnt += 1

        subset_node = []
        conn_dict = {}
        cnt = 1
        for n in range(node_count):
            nlon = node_coords[n, 0].values
            nlat = node_coords[n, 1].values

            l1 = np.logical_or(nlon <= self.lon1, nlon >= self.lon2)
            l2 = np.logical_or(nlat <= self.lat1, nlat >= self.lat2)

            if np.logical_or(l1, l2):
                conn_dict[n + 1] = -9999
            else:
                subset_node.append(n)
                conn_dict[n + 1] = cnt
                cnt += 1

            # -- reverse logic
            # l1 = np.logical_and(nlon >= self.lon1,nlon <= self.lon2)
            # l2 = np.logical_and(nlat >= self.lat1,nlat <= self.lat2)
            # if np.any(l1) and np.any(l2):
            #    subset_node.append(n)
            #    conn_dict[n+1] = cnt
            #    cnt+=1
            # else:
            #    conn_dict[n+1] = -9999

        return node_coords, subset_element, subset_node, conn_dict

    @staticmethod
    def write_mesh(f_in, subset_element, subset_node, conn_dict, mesh_out):
        """
        This function writes out the subsetted mesh file.
        """
        corner_pairs = f_in.variables["nodeCoords"][subset_node,]
        variables = f_in.variables
        global_attributes = f_in.attrs

        max_node_dim = len(f_in["maxNodePElement"])

        elem_count = len(subset_element)
        elem_conn_out = np.empty(shape=[elem_count, max_node_dim])
        elem_conn_index = f_in.variables["elementConn"][subset_element,]

        for n in range(elem_count):
            for m in range(max_node_dim):
                ndx = int(elem_conn_index[n, m])
                elem_conn_out[n, m] = conn_dict[ndx]

        num_elem_conn_out = np.empty(
            shape=[
                elem_count,
            ]
        )
        num_elem_conn_out[:] = f_in.variables["numElementConn"][subset_element,]

        center_coords_out = np.empty(shape=[elem_count, 2])
        center_coords_out[:, :] = f_in.variables["centerCoords"][subset_element, :]

        if "elementMask" in variables:
            elem_mask_out = np.empty(
                shape=[
                    elem_count,
                ]
            )
            elem_mask_out[:] = f_in.variables["elementMask"][subset_element,]

        if "elementArea" in variables:
            elem_area_out = np.empty(
                shape=[
                    elem_count,
                ]
            )
            elem_area_out[:] = f_in.variables["elementArea"][subset_element,]

        # -- create output dataset
        f_out = xr.Dataset()

        f_out["nodeCoords"] = xr.DataArray(
            corner_pairs, dims=("nodeCount", "coordDim"), attrs={"units": "degrees"}
        )

        f_out["elementConn"] = xr.DataArray(
            elem_conn_out,
            dims=("elementCount", "maxNodePElement"),
            attrs={"long_name": "Node indices that define the element connectivity"},
        )
        f_out.elementConn.encoding = {"dtype": np.int32}

        f_out["numElementConn"] = xr.DataArray(
            num_elem_conn_out,
            dims=("elementCount"),
            attrs={"long_name": "Number of nodes per element"},
        )
        f_out.numElementConn.encoding = {"dtype": np.int32}

        f_out["centerCoords"] = xr.DataArray(
            center_coords_out,
            dims=("elementCount", "coordDim"),
            attrs={"units": "degrees"},
        )

        # -- add mask
        if "elementMask" in variables:
            f_out["elementMask"] = xr.DataArray(
                elem_mask_out, dims=("elementCount"), attrs={"units": "unitless"}
            )
            f_out.elementMask.encoding = {"dtype": np.int32}

        if "elementArea" in variables:
            f_out["elementArea"] = xr.DataArray(
                elem_area_out, dims=("elementCount"), attrs={"units": "unitless"}
            )

        # -- setting fill values
        for var in variables:
            if "_FillValue" in f_in[var].encoding:
                f_out[var].encoding["_FillValue"] = f_in[var].encoding["_FillValue"]
            else:
                f_out[var].encoding["_FillValue"] = None

        # -- add global attributes
        for attr in global_attributes:
            if attr != "timeGenerated":
                f_out.attrs[attr] = global_attributes[attr]

        f_out.attrs = {
            "title": "ESMF unstructured grid file for a region",
            "created_by": "subset_data",
            "date_created": "{}".format(datetime.now()),
        }

        f_out.to_netcdf(mesh_out)
        logger.info("Successfully created file (mesh_out) %s", mesh_out)

    def write_shell_commands(self, namelist):
        """
        writes out xml commands commands to a file (i.e. shell_commands) for single-point runs
        """
        # write_to_file surrounds text with newlines
        with open(namelist, "w") as nl_file:
            self.write_to_file("# Change below line if you move the subset data directory", nl_file)
            self.write_to_file("./xmlchange {}={}".format(USRDAT_DIR, self.out_dir), nl_file)
            self.write_to_file("./xmlchange ATM_DOMAIN_MESH={}".format(str(self.mesh)), nl_file)
            self.write_to_file("./xmlchange LND_DOMAIN_MESH={}".format(str(self.mesh)), nl_file)
            self.write_to_file("./xmlchange MASK_MESH={}".format(str(str(self.mesh))), nl_file)
            self.write_to_file("./xmlchange ATM_NX={}".format(str(str(self.ni))), nl_file)
            self.write_to_file("./xmlchange LND_NX={}".format(str(str(self.ni))), nl_file)
            self.write_to_file("./xmlchange ATM_NY={}".format(str(str(self.nj))), nl_file)
            self.write_to_file("./xmlchange LND_NY={}".format(str(str(self.nj))), nl_file)
