"""
This module includes the definition for a parent class for SinglePointCase
and RegionalCase. The common functionalities of SinglePointCase and
RegionalCase are defined in this Class.
"""

# -- Import libraries

# -- standard libraries
import os.path
import logging
from collections import namedtuple

from datetime import date
from getpass import getuser

# -- 3rd party libraries
import numpy as np
import xarray as xr

# -- import local classes for this script
from ctsm.utils import abort
from ctsm.git_utils import get_ctsm_git_short_hash

USRDAT_DIR = "CLM_USRDAT_DIR"
logger = logging.getLogger(__name__)

# named tuple for datm input/output files and folder names
DatmFiles = namedtuple(
    "DatmFiles",
    "indir outdir fdomain_in dir_solar dir_prec dir_tpqw tag_solar tag_prec tag_tpqw name_solar "
    "name_prec name_tpqw ",
)


class BaseCase:
    """
    Parent class to SinglePointCase and RegionalCase
    ...
    Attributes
    ----------
    create_domain : bool
        flag for creating domain file
    create_surfdata : bool
        flag for creating surface dataset
    create_landuse : bool
        flag for creating landuse file
    create_datm : bool
        flag for creating DATM files
    create_user_mods
        flag for creating a user_mods directory
    overwrite : bool
        flag for overwriting if the file already exists

    Methods
    -------
    create_1d_coord(filename, lon_varname , lat_varname,x_dim , y_dim )
        create 1d coordinate variables to enable sel() method
    update_metadata(nc)
        Class method for adding some new attributes (such as date, username) and
        remove the old attributes from the netcdf file.
    write_to_file:
        Writes text to a file, surrounding text with \n characters
    write_to_netcdf:
        write xarray dataset to netcdf
    """

    def __init__(
        self,
        *,
        create_domain,
        create_surfdata,
        create_landuse,
        create_datm,
        create_user_mods,
        overwrite,
    ):
        """
        Initializes BaseCase with the given arguments.

        Parameters
        ----------
        create_domain : bool
            Flag for creating domain file a region/single point
        create_surfdata : bool
            Flag for creating domain file a region/single point
        create_landuse : bool
            Flag for creating landuse file a region/single point
        create_datmdata : bool
            Flag for creating datm files a region/single point
        create_user_mods : bool
            Flag for creating user mods directories and files for running CTSM
        overwrite : bool
            flag for overwriting if the file already exists
        """
        self.create_domain = create_domain
        self.create_surfdata = create_surfdata
        self.create_landuse = create_landuse
        self.create_datm = create_datm
        self.create_user_mods = create_user_mods
        self.overwrite = overwrite

    def __str__(self):
        """
        Converts ingredients of the BaseCase to string for printing.
        """
        return "{}\n{}".format(
            str(self.__class__),
            "\n".join(
                (
                    "{} = {}".format(str(key), str(self.__dict__[key]))
                    for key in sorted(self.__dict__)
                )
            ),
        )

    @staticmethod
    def create_1d_coord(filename, lon_varname, lat_varname, x_dim, y_dim):
        """
        Create 1d coordinate variables for a netcdf file to enable sel() method

        Parameters
        ----------
            filename (str) : name of the netcdf file
            lon_varname (str) : variable name that has 2d lon
            lat_varname (str) : variable name that has 2d lat
            x_dim (str) : dimension name in X -- lon
            y_dim (str): dimension name in Y -- lat

        Raises
        ------
            None

        Returns
        -------
            f_out (xarray Dataset): Xarray Dataset with 1-d coords

        """

        f_in = None
        if os.path.exists(filename):
            logger.debug("Open file: %s", filename)

            f_in = xr.open_dataset(filename)
        else:
            err_msg = "File not found : " + filename
            abort(err_msg)

        # create 1d coordinate variables to enable sel() method
        lon0 = np.asarray(f_in[lon_varname][0, :])
        lat0 = np.asarray(f_in[lat_varname][:, 0])
        lon = xr.DataArray(lon0, name="lon", dims=x_dim, coords={x_dim: lon0})
        lat = xr.DataArray(lat0, name="lat", dims=y_dim, coords={y_dim: lat0})

        f_out = f_in.assign({"lon": lon, "lat": lat})

        f_out.reset_coords([lon_varname, lat_varname])
        f_in.close()
        return f_out

    @staticmethod
    def update_metadata(nc_file):
        """
        Class method for adding some new attributes (such as date, username) and
        remove the old attributes from the netcdf file.

        Parameters
        ----------
            nc (xarray dataset) :
                Xarray dataset of netcdf file that we'd want to update it's metadata.

        Raises
        ------
            None

        Returns
        ------
            None

        """
        # update attributes
        today = date.today()
        today_string = today.strftime("%Y-%m-%d")

        # get git hash
        sha = get_ctsm_git_short_hash()

        nc_file.attrs["Created_on"] = today_string
        nc_file.attrs["Created_by"] = getuser()
        nc_file.attrs["Created_with"] = "./subset_data" + " -- " + sha

        # delete unrelated attributes if they exist
        del_attrs = [
            "source_code",
            "SVN_url",
            "hostname",
            "history",
            "History_Log",
            "Logname",
            "Host",
            "Version",
            "Compiler_Optimized",
        ]
        attr_list = nc_file.attrs

        for attr in del_attrs:
            if attr in attr_list:
                logger.debug("This attr should be deleted : %s", attr)
                del nc_file.attrs[attr]

    @staticmethod
    def write_to_file(text, file_out):
        """
        Writes text to a file, surrounding text with \n characters
        """
        file_out.write("\n{}\n".format(text))

    def write_to_netcdf(self, xr_ds, nc_fname):
        """
        Writes a netcdf file if
            - the file does not exist.
            or
            - overwrite flag is chosen.

        Args:
            xr_ds : Xarray Dataset
                The xarray dataset that we are write out to netcdf file.
            nc_fname : str
                Netcdf file name
        Raises:
            Error and aborts the code if the file exists and --overwrite is not used.
        """
        if not os.path.exists(nc_fname) or self.overwrite:
            # mode 'w' overwrites file
            xr_ds.to_netcdf(path=nc_fname, mode="w", format="NETCDF3_64BIT")
        else:
            err_msg = (
                "File "
                + nc_fname
                + " already exists."
                + "\n Either remove the file or use "
                + "--overwrite to overwrite the existing files."
            )
            abort(err_msg)
