"""
This module extends the mesh type for advanced plotting
"""

import logging

import numpy as np

# -- libraries for plotting mesh (make_mesh_plot)
# Turn import error off in case you are using a
# conda environment that doesn't include these.
# pylint: disable=import-error
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from ctsm.site_and_regional.mesh_type import MeshType

logger = logging.getLogger(__name__)


class MeshPlotType(MeshType):
    """
    Extend mesh type with some advanced plotting capability
    """

    def make_mesh_plot(self, plot_regional, plot_global, args):
        """
        Create plots for the ESMF mesh file

        Parameters
        ----------
        plot_regional : str
            The path to write the ESMF meshfile regional plot
        plot_global : str
            The path to write the ESMF meshfile global plot
        """

        self.mesh_plot(plot_regional, args, regional=True)
        self.mesh_plot(plot_global, args, regional=False)

    def mesh_plot(self, plot_file, args, regional):
        """Make a plot of a mesh file in either a regional or global grid"""
        # -- regional settings
        if regional:
            plt.figure(num=None, figsize=(15, 13), facecolor="w", edgecolor="k")
            # pylint: disable=abstract-class-instantiated
            ax = plt.axes(projection=ccrs.PlateCarree())
            plot_type = "regional"
            line_width = 1
            marker = "x"
            marker_size = 50
        # global settings
        else:
            fig = plt.figure(num=None, figsize=(15, 10), facecolor="w", edgecolor="k")
            # pylint: disable=abstract-class-instantiated
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
            plot_type = "global"
            line_width = 0.5
            marker = "o"
            marker_size = 0.1
        if args.no_center_coords:
            marker_size = 0

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
        if not regional:
            ax.set_global()

        # -- plot corner coordinates
        # clats, clons = self.node_coords.T.compute()
        # elem_conn_vals = self.elem_conn.compute()
        clats, clons = self.node_coords.T
        elem_conn_vals = self.elem_conn
        element_counts = elem_conn_vals.shape[0]
        # Scale marker and line size for a large number of points
        if element_counts > 60000:
            marker_size = 0.01
            line_width = line_width * 0.05

        for index in range(element_counts):
            conns = [int(x) - 1 for x in elem_conn_vals[index]]

            lat_corners = clats[conns]
            lon_corners = clons[conns]
            poly_corners = np.zeros((len(lat_corners), 2))
            poly_corners[:, 1] = lon_corners
            poly_corners[:, 0] = lat_corners
            # pylint: disable=abstract-class-instantiated
            poly = mpatches.Polygon(
                poly_corners,
                closed=True,
                ec="black",
                transform=ccrs.PlateCarree(),
                zorder=10,
                facecolor="none",
                linewidth=line_width,
            )
            ax.add_patch(poly)

        # -- plot center coordinates
        # clon, clat = self.center_coords.T.compute()
        clon, clat = self.center_coords.T

        # pylint: disable=abstract-class-instantiated
        ax.scatter(
            clon,
            clat,
            color="tomato",
            marker=marker,
            s=marker_size,
            transform=ccrs.PlateCarree(),
            zorder=11,
        )
        if regional:
            lc_colors = {
                "Corner Coordinates": "black",  # value=0
                "Center Coordinates": "tomato",  # value=1
            }
            labels, handles = zip(
                *[(k, mpatches.Rectangle((0, 0), 1, 1, facecolor=v)) for k, v in lc_colors.items()]
            )

            if not args.no_center_coords:
                ax.legend(handles, labels)

        plt.savefig(plot_file, bbox_inches="tight", dpi=args.dpi)

        logger.info("Successfully created %s plots for ESMF Mesh file : %s", plot_type, plot_file)
