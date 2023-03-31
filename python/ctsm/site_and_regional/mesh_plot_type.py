"""
This module extends the mesh type for advanced plotting
"""
import logging

import numpy as np

# -- libraries for plotting mesh (make_mesh_plot)
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
        # clats, clons = self.node_coords.T.compute()
        # elem_conn_vals = self.elem_conn.compute()
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
        # clon, clat = self.center_coords.T.compute()
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
        # clats, clons = self.node_coords.T.compute()
        # elem_conn_vals = self.elem_conn.compute()
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
        # clon, clat = self.center_coords.T.compute()
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
