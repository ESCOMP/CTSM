import os
import subprocess
import logging

import numpy as np
import xarray as xr

from datetime import date
from getpass import getuser

myname = getuser()
USRDAT_DIR = "CLM_USRDAT_DIR"

logger = logging.getLogger(__name__)

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
    Methods
    -------
    create_1d_coord(filename, lon_varname , lat_varname,x_dim , y_dim )
        create 1d coordinate variables to enable sel() method
    add_tag_to_filename(filename, tag)
       add a tag and timetag to a filename ending with
       [._]cYYMMDD.nc or [._]YYMMDD.nc
    update_metadata(self, nc)
       updates metadata for a netcdf file and removes attributes that should not be there
    """

    def __init__(self, create_domain, create_surfdata, create_landuse, create_datm, create_user_mods):
        self.create_domain = create_domain
        self.create_surfdata = create_surfdata
        self.create_landuse = create_landuse
        self.create_datm = create_datm
        self.create_user_mods = create_user_mods

    def __str__(self):
        return "{}\n{}".format(str(self.__class__), "\n".join(
            ("{} = {}".format(str(key), str(self.__dict__[key])) for key in sorted(self.__dict__))))

    @staticmethod
    def create_1d_coord(filename, lon_varname, lat_varname, x_dim, y_dim):
        """
        lon_varname : variable name that has 2d lon
        lat_varname : variable name that has 2d lat
        x_dim: dimension name in X -- lon
        y_dim: dimension name in Y -- lat
        """
        logging.debug("Open file: " + filename)
        f1 = xr.open_dataset(filename)

        # create 1d coordinate variables to enable sel() method
        lon0 = np.asarray(f1[lon_varname][0, :])
        lat0 = np.asarray(f1[lat_varname][:, 0])
        lon = xr.DataArray(lon0, name="lon", dims=x_dim, coords={x_dim: lon0})
        lat = xr.DataArray(lat0, name="lat", dims=y_dim, coords={y_dim: lat0})

        f2 = f1.assign({"lon": lon, "lat": lat})

        f2.reset_coords([lon_varname, lat_varname])
        f1.close()
        return f2

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
        if basename[cend] == "c":
            cend = cend - 1
        if (basename[cend] != ".") and (basename[cend] != "_"):
            logging.error("Trouble figuring out where to add tag to filename:" + filename)
            os.abort()
        today = date.today()
        today_string = today.strftime("%y%m%d")
        return basename[:cend] + "_" + tag + "_c" + today_string + ".nc"

    def update_metadata(self, nc):
        # update attributes
        today = date.today()
        today_string = today.strftime("%Y-%m-%d")

        # get git hash
        sha = self.get_git_sha()

        nc.attrs["Created_on"] = today_string
        nc.attrs["Created_by"] = myname
        nc.attrs["Created_with"] = os.path.abspath(__file__) + " -- " + sha

        # delete unrelated attributes if they exist
        del_attrs = [
            "source_code",
            "SVN_url",
            "hostname",
            "history" "History_Log",
            "Logname",
            "Host",
            "Version",
            "Compiler_Optimized",
        ]
        attr_list = nc.attrs

        for attr in del_attrs:
            if attr in attr_list:
                logging.debug ("This attr should be deleted : "+ attr)
                del nc.attrs[attr]

        # for attr, value in attr_list.items():
        #    print (attr + " = "+str(value))

    @staticmethod
    def get_git_sha():
        """
        Returns Git short SHA for the current directory.
        """
        try:
            sha = (subprocess.check_output(["git", "-C", os.path.dirname(__file__), "rev-parse", "--short", "HEAD"]).strip().decode())
        except subprocess.CalledProcessError:
            sha = "NOT-A-GIT-REPOSITORY"
        return sha

    @staticmethod
    def write_to_file(text, file):
        """
        Writes text to a file, surrounding text with \n characters
        """
        file.write("\n{}\n".format(text))
