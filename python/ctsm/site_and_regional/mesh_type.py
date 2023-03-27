"""
This module includes the definition and functions for defining a Grid or Mesh.
This enables creating ESMF mesh file (unstructured grid file)for valid 1D or 2D lats and lons.
"""
import os
import logging
import argparse
import datetime

import numpy as np
import xarray as xr
#import dask.array as da
#import dask.dataframe as dd
import pandas

# -- libraries for plotting mesh (make_mesh_plot)
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from ctsm.utils import abort

logger = logging.getLogger(__name__)


class MeshType:
    """
    An object for describing mesh or grid files.
    """

    def __init__(self, center_lats, center_lons, mask=None, mesh_name=None):
        """
        Construct a mesh object

        Attributes
        ----------
        center_lats : xarray DataArray
            array of latitudes (either 1D or 2D)
        center_lons : xarray DataArray
            array longitudes (either 1D or 2D)
        mesh_name : str or None, optional
            Name of the mesh
        mask : np array or None
            numpy array that include landmask

        Methods
        -------
        check_lat_lon_dims:
            check if dimensions of latitude and longitude is valid.

        create_artificial_mask:
            create an artificial mask of ones if mask does not exits.

        create_2d_coords:
            if 1d coords is provided, this method creates 2d coords from 1d coords.

        calculate_corners:
            calculate corner coordinates for each polygon in the grid

        calculate_node_coords:
            extract node coordiantes for nodeCoords variable on mesh file

        calculate_elem_conn
            calculate element connectivity for elementConn variable on mesh file

        create_esmf
            write mesh file to netcdf file.
        """

        self.mesh_name = mesh_name
        self.center_lats = center_lats
        self.center_lons = center_lons

        # -- dims of lat and lon (1d or 2d)
        self.lat_dims = len(self.center_lats.dims)
        self.lon_dims = len(self.center_lons.dims)
        self.check_lat_lon_dims()

        if mask is None:
            self.create_artificial_mask()
        else:
            #self.mask = da.from_array(np.array(mask.astype(np.int8)))
            self.mask = np.array(mask.astype(np.int8))

    def check_lat_lon_dims(self):
        """
        Check latitude and longitude dimensions to make sure they are valid.

        -------------
        Raises:
            Error (ArgumentTypeError):
                If the provided latitude has dimension >2.
            Error (ArgumentTypeError):
                If the provided longitude has dimension >2.
        """

        if self.lat_dims not in [1, 2]:
            err_msg = """
                Unrecognized grid! \n
                The dimension of latitude should be either 1 or 2 but it is {}.
                """.format(
                self.lat_dims
            )
            raise argparse.ArgumentTypeError(err_msg)

        if self.lon_dims not in [1, 2]:
            err_msg = """
                Unrecognized grid! \n
                The dimension of longitude should be either 1 or 2 but it is {}.
                """.format(
                self.lon_dims
            )
            raise argparse.ArgumentTypeError(err_msg)

    def create_artificial_mask(self):
        """
        create an artificial mask of 1 (i.e. all ocean) if no land mask is provided.
        """

        logger.info("Creating an artificial mask for this region...")

        if self.lat_dims == 1:
            # -- 1D mask (lat x lon)
            lats_size = self.center_lats.size
            lons_size = self.center_lons.size
            mask = np.ones([lons_size, lats_size], dtype=np.int8)
        elif self.lat_dims == 2:
            # -- 2D mask
            mask = np.ones(self.center_lats.shape, dtype=np.int8)
        #mask_da = da.from_array(mask)
        mask_da = mask
        self.mask = mask_da

    def create_2d_coords(self):
        """
        Create 2d center points for our mesh
        and convert them to Dask Array.
        """
        self.center_lats = self.center_lats.astype(np.float64, copy=False)
        self.center_lons = self.center_lons.astype(np.float64, copy=False)

        if self.lat_dims == 1:
            # -- 1D lats and lons
            lats_size = self.center_lats.size
            lons_size = self.center_lons.size

            # -- convert center points from 1d to 2d
            self.center_lat2d = np.broadcast_to(
                self.center_lats[:], (lons_size, lats_size)
            )
            self.center_lon2d = np.broadcast_to(
                self.center_lons[:], (lons_size, lats_size)
            )
        elif self.lat_dims == 2:
            # -- 2D lats and lons
            dims = self.center_lons.shape

            # -- convert 2D lats and lons to number x and y
            lons_size = dims[0]
            lats_size = dims[1]

            # -- convert to dask array
            #self.center_lat2d = da.from_array(np.array(self.center_lats))
            #self.center_lon2d = da.from_array(np.array(self.center_lons))
            self.center_lat2d = np.array(self.center_lats)
            self.center_lon2d = np.array(self.center_lons)

    def calculate_corners(self, unit="degrees"):
        """
        calculate corner coordinates by averaging adjacent cells

        Parameters
        ----------
        unit : {'degrees', 'radians'}, optional
            The unit of corner coordinates.
        """

        self.create_2d_coords()
        # -- pad center_lats for calculating edge gridpoints
        # -- otherwise we cannot calculate the corner coords
        # -- for the edge rows/columns.

        #padded_lat2d = np.pad(self.center_lat2d.compute(), (1, 1), mode="reflect", reflect_type="odd")
        padded_lat2d = np.pad(self.center_lat2d, (1, 1), mode="reflect", reflect_type="odd")

        # -- pad center_lons for calculating edge grids
        padded_lon2d = np.pad(self.center_lon2d, (1, 1), mode="reflect", reflect_type="odd")

        # -- calculate corner lats for each grid
        north_east = (
            padded_lat2d[1:-1, 1:-1]
            + padded_lat2d[0:-2, 1:-1]
            + padded_lat2d[1:-1, 2:]
            + padded_lat2d[0:-2, 2:]
        ) / 4.0
        north_west = (
            padded_lat2d[1:-1, 1:-1]
            + padded_lat2d[0:-2, 1:-1]
            + padded_lat2d[1:-1, 0:-2]
            + padded_lat2d[0:-2, 0:-2]
        ) / 4.0
        south_west = (
            padded_lat2d[1:-1, 1:-1]
            + padded_lat2d[1:-1, 0:-2]
            + padded_lat2d[2:, 1:-1]
            + padded_lat2d[2:, 0:-2]
        ) / 4.0
        south_east = (
            padded_lat2d[1:-1, 1:-1]
            + padded_lat2d[1:-1, 2:]
            + padded_lat2d[2:, 1:-1]
            + padded_lat2d[2:, 2:]
        ) / 4.0

        # -- order counter-clockwise
        self.corner_lats = np.stack(
            [
                north_west.T.reshape((-1,)).T,
                south_west.T.reshape((-1,)).T,
                south_east.T.reshape((-1,)).T,
                north_east.T.reshape((-1,)).T,
            ],
            axis=1,
        )

        # -- calculate corner lons for each grid
        north_east = (
            padded_lon2d[1:-1, 1:-1]
            + padded_lon2d[0:-2, 1:-1]
            + padded_lon2d[1:-1, 2:]
            + padded_lon2d[0:-2, 2:]
        ) / 4.0
        north_west = (
            padded_lon2d[1:-1, 1:-1]
            + padded_lon2d[0:-2, 1:-1]
            + padded_lon2d[1:-1, 0:-2]
            + padded_lon2d[0:-2, 0:-2]
        ) / 4.0
        south_west = (
            padded_lon2d[1:-1, 1:-1]
            + padded_lon2d[1:-1, 0:-2]
            + padded_lon2d[2:, 1:-1]
            + padded_lon2d[2:, 0:-2]
        ) / 4.0
        south_east = (
            padded_lon2d[1:-1, 1:-1]
            + padded_lon2d[1:-1, 2:]
            + padded_lon2d[2:, 1:-1]
            + padded_lon2d[2:, 2:]
        ) / 4.0

        # -- order counter-clockwise
        self.corner_lons = np.stack(
            [
                north_west.T.reshape((-1,)).T,
                south_west.T.reshape((-1,)).T,
                south_east.T.reshape((-1,)).T,
                north_east.T.reshape((-1,)).T,
            ],
            axis=1,
        )
        self.unit = unit

    def calculate_node_coords(self):
        """
        Calculates coordinates of each node (for 'nodeCoords' in ESMF mesh).
        In ESMF mesh, 'nodeCoords' is a two-dimensional array with dimension ('nodeCount','coordDim')
        """
        # -- create an array of corner pairs
        corner_pairs = np.stack(
            [self.corner_lons.T.reshape((-1,)).T, self.corner_lats.T.reshape((-1,)).T],
            axis=1,
        )

        # -- convert to float32 to find duplicates

        # -- remove coordinates that are shared between the elements
        #node_coords = dd.from_dask_array(corner_pairs).drop_duplicates().values
        #node_coords.compute_chunk_sizes()
        node_coords = pandas.DataFrame( data=corner_pairs ).drop_duplicates().values

        # -- check size of unique coordinate pairs
        dims = self.mask.shape
        nlon = dims[0]
        nlat = dims[1]
        elem_conn_size = nlon * nlat + nlon + nlat + 1

        self.node_coords = node_coords

        # -- error check to avoid issues later
        if self.node_coords.shape[0] != elem_conn_size:
            logger.warning(
                "The size of unique coordinate pairs is {} but expected size is {}!".format(
                    self.node_coords.shape[0], elem_conn_size
                )
            )
            logger.warning("This may result your simulation to crash later.")
            #abort( "Expected size for element connections is wrong" )

    def calculate_elem_conn(self):
        """
        Calculate element connectivity (for 'elementConn' in ESMF mesh).
        In ESMF mesh, 'elementConn' describes how the nodes are connected together.
        """
        #corners = dd.concat(
        #    [
        #        dd.from_dask_array(corner)
        #        for corner in [
        #            self.corner_lons.T.reshape((-1,)).T,
        #            self.corner_lats.T.reshape((-1,)).T,
        #        ]
        #    ],
        #    axis=1,
        #)
        #corners.columns = ["lon", "lat"]

        #elem_conn = corners.compute().groupby(["lon", "lat"], sort=False).ngroup() + 1
        #elem_conn = da.from_array(elem_conn.to_numpy())

        ## -- reshape to write to ESMF
        #self.elem_conn = elem_conn.T.reshape((4, -1)).T
        corners = pandas.concat(
            [
                pandas.DataFrame(corner)
                for corner in [
                    self.corner_lons.T.reshape((-1,)).T,
                    self.corner_lats.T.reshape((-1,)).T,
                ]
            ],
            axis=1,
        )
        corners.columns = ["lon", "lat"]

        elem_conn = corners.groupby(["lon", "lat"], sort=False).ngroup() + 1
        elem_conn = elem_conn.to_numpy()

        # -- reshape to write to ESMF
        self.elem_conn = elem_conn.T.reshape((4, -1)).T

    def create_esmf(self, mesh_fname, area=None):
        """
        Create an ESMF mesh file for the mesh

        Parameters
        ----------
        mesh_fname : str
            The path to write the ESMF meshfile

        area : numpy.ndarray or None
            Array containing element areas for the ESMF mesh file
            If None, ESMF calculates element areas internally.
        """
        # -- calculate node coordinates
        self.calculate_node_coords()

        # -- calculate element connections
        self.calculate_elem_conn()

        self.center_coords = np.stack(
            [
                self.center_lon2d.T.reshape((-1,)).T,
                self.center_lat2d.T.reshape((-1,)).T,
            ],
            axis=1,
        )
        # create output Xarray dataset
        ds_out = xr.Dataset()

        ds_out["origGridDims"] = xr.DataArray(
            np.array(self.center_lon2d.shape, dtype=np.int32), dims=("origGridRank")
        )
        ds_out["nodeCoords"] = xr.DataArray(
            self.node_coords, dims=("nodeCount", "coordDim"), attrs={"units": self.unit}
        )
        ds_out["elementConn"] = xr.DataArray(
            self.elem_conn,
            dims=("elementCount", "maxNodePElement"),
            attrs={
                "long_name": "Node indices that define the element connectivity",
                "_FillValue": -1,
            },
        )
        ds_out.elementConn.encoding = {"dtype": np.int32}

        ds_out["numElementConn"] = xr.DataArray(
            4 * np.ones(self.center_lon2d.size, dtype=np.int32),
            dims=("elementCount"),
            attrs={"long_name": "Number of nodes per element"},
        )
        ds_out["centerCoords"] = xr.DataArray(
            self.center_coords,
            dims=("elementCount", "coordDim"),
            attrs={"units": self.unit},
        )

        # -- add mask
        ds_out["elementMask"] = xr.DataArray(
            self.mask.T.reshape((-1,)).T,
            dims=("elementCount"),
            attrs={"units": "unitless", "_FillValue": -9999.0},
        )
        ds_out.elementMask.encoding = {"dtype": np.int32}

        # -- add area if provided
        if area is not None:
            da_area = np.array(area)
            ds_out["elementArea"] = xr.DataArray(
                da_area.T.reshape((-1,)).T,
                dims=("elementCount"),
                attrs={"units": "radians^2", "long_name": "area weights"},
            )

        # -- force no '_FillValue' if not specified (default Nan)
        for var in ds_out.variables:
            if "_FillValue" not in ds_out[var].encoding:
                ds_out[var].encoding["_FillValue"] = None

        # -- add global attributes
        if self.mesh_name:
            ds_out.attrs["title"] = "ESMF unstructured grid file  " + self.mesh_name
        else:
            ds_out.attrs["title"] = "ESMF unstructured grid file  "
        ds_out.attrs["gridType"] = "unstructured mesh"
        ds_out.attrs["version"] = "0.9"
        ds_out.attrs["conventions"] = "ESMFMESH"
        ds_out.attrs["date_created"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # -- write Xarray dataset to file
        if mesh_fname is not None:
            logger.info("Writing ESMF Mesh file to : %s", mesh_fname)
            ds_out.to_netcdf(mesh_fname)
            logger.info("Successfully created ESMF Mesh file : %s", mesh_fname)

        self.make_plot = True

        if self.make_plot:
            plot_regional = os.path.splitext(mesh_fname)[0] + "_regional" + ".png"

            plot_global = os.path.splitext(mesh_fname)[0] + "_global" + ".png"

            self.make_mesh_plot(plot_regional, plot_global)

    def make_mesh_plot(self, plot_regional, plot_global):
        """
        Create a plot for the ESMF mesh file

        Parameters
        ----------
        plot_regional : str
            The path to write the ESMF meshfile regional plot
        plot_global : str
            The path to write the ESMF meshfile global plot
        """

        # -- regional plot
        plt.figure(num=None, figsize=(15, 13), facecolor="w", edgecolor="k")
        ax = plt.axes(projection=ccrs.PlateCarree())

        ax.add_feature(cfeature.COASTLINE, edgecolor="black")
        ax.add_feature(cfeature.BORDERS, edgecolor="black")
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, edgecolor="black")
        ax.add_feature(cfeature.LAKES, edgecolor="black")
        ax.add_feature(cfeature.RIVERS)
        ax.gridlines(
            color="black",
            linestyle="dotted",
            draw_labels=True,
        )

        # -- plot corner coordinates
        #clats, clons = self.node_coords.T.compute()
        #elem_conn_vals = self.elem_conn.compute()
        clats, clons = self.node_coords.T
        elem_conn_vals = self.elem_conn
        element_counts = elem_conn_vals.shape[0]

        for index in range(element_counts):
            conns = [int(x) - 1 for x in elem_conn_vals[index]]

            lat_corners = clats[conns]
            lon_corners = clons[conns]
            poly_corners = np.zeros((len(lat_corners), 2))
            poly_corners[:, 1] = lon_corners
            poly_corners[:, 0] = lat_corners
            poly = mpatches.Polygon(
                poly_corners,
                closed=True,
                ec="black",
                lw=1,
                transform=ccrs.PlateCarree(),
                zorder=10,
                facecolor="none",
            )
            ax.add_patch(poly)

        # -- plot center coordinates
        #clon, clat = self.center_coords.T.compute()
        clon, clat = self.center_coords.T

        ax.scatter(
            clon,
            clat,
            color="tomato",
            marker="x",
            transform=ccrs.PlateCarree(),
            zorder=11,
        )
        lc_colors = {
            "Corner Coordinates": "black",  # value=0
            "Center Coordinates": "tomato",  # value=1
        }
        labels, handles = zip(
            *[(k, mpatches.Rectangle((0, 0), 1, 1, facecolor=v)) for k, v in lc_colors.items()]
        )

        ax.legend(handles, labels)

        plt.savefig(plot_regional, bbox_inches="tight")

        logger.info("Successfully created regional plots for ESMF Mesh file : %s", plot_regional)

        # -- global plot
        fig = plt.figure(num=None, figsize=(15, 10), facecolor="w", edgecolor="k")

        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

        ax.add_feature(cfeature.COASTLINE, edgecolor="black")
        ax.add_feature(cfeature.BORDERS, edgecolor="black")
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, edgecolor="black")
        ax.add_feature(cfeature.LAKES, edgecolor="black")
        ax.add_feature(cfeature.RIVERS)
        ax.gridlines(
            color="black",
            linestyle="dotted",
            draw_labels=True,
        )
        ax.set_global()

        # -- plot corner coordinates
        #clats, clons = self.node_coords.T.compute()
        #elem_conn_vals = self.elem_conn.compute()
        clats, clons = self.node_coords.T
        elem_conn_vals = self.elem_conn
        element_counts = elem_conn_vals.shape[0]

        for index in range(element_counts):
            conns = [int(x) - 1 for x in elem_conn_vals[index]]

            lat_corners = clats[conns]
            lon_corners = clons[conns]
            poly_corners = np.zeros((len(lat_corners), 2))
            poly_corners[:, 1] = lon_corners
            poly_corners[:, 0] = lat_corners
            poly = mpatches.Polygon(
                poly_corners,
                closed=True,
                ec="black",
                transform=ccrs.PlateCarree(),
                zorder=10,
                facecolor="none",
                linewidth=0.5,
            )
            ax.add_patch(poly)

        # -- plot center coordinates
        #clon, clat = self.center_coords.T.compute()
        clon, clat = self.center_coords.T

        ax.scatter(
            clon,
            clat,
            color="tomato",
            marker="o",
            s=1,
            transform=ccrs.PlateCarree(),
            zorder=11,
        )

        ax.legend(handles, labels)

        plt.savefig(plot_global, bbox_inches="tight")

        logger.info("Successfully created regional plots for ESMF Mesh file : %s", plot_global)
