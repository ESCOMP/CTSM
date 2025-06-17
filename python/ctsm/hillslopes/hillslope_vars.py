"""
A useful class for processing hillslope files
"""

import shutil
import numpy as np

# The below "pylint: disable" is because pylint complains that netCDF4 has no member Dataset, even
# though it does.
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm.hillslopes.hillslope_utils import add_variable_nc, add_longxy_latixy_nc, NETCDF_FORMAT


class HillslopeVars:
    """
    Fields to be added to hillslope_file
    """

    # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        ncolumns_per_gridcell,
        nhillslope,
        n_lat,
        n_lon,
        *,
        recurse=True,
        incl_latlon=False,
        incl_pftndx=False,
    ):

        # Variables that will actually be saved
        self.h_elev = np.zeros((ncolumns_per_gridcell, n_lat, n_lon))
        self.h_dist = np.zeros((ncolumns_per_gridcell, n_lat, n_lon))
        self.h_area = np.zeros((ncolumns_per_gridcell, n_lat, n_lon))
        self.h_slope = np.zeros((ncolumns_per_gridcell, n_lat, n_lon))
        self.h_aspect = np.zeros((ncolumns_per_gridcell, n_lat, n_lon))
        self.h_width = np.zeros((ncolumns_per_gridcell, n_lat, n_lon))
        self.h_bedrock = np.zeros((ncolumns_per_gridcell, n_lat, n_lon))
        self.h_stream_depth = np.zeros((n_lat, n_lon))
        self.h_stream_width = np.zeros((n_lat, n_lon))
        self.h_stream_slope = np.zeros((n_lat, n_lon))
        self.nhillcolumns = np.zeros((n_lat, n_lon), dtype=int)
        self.pct_hillslope = np.zeros((nhillslope, n_lat, n_lon))
        self.hillslope_index = np.zeros((ncolumns_per_gridcell, n_lat, n_lon), dtype=int)
        self.column_index = np.zeros((ncolumns_per_gridcell, n_lat, n_lon), dtype=int)
        self.downhill_column_index = np.zeros((ncolumns_per_gridcell, n_lat, n_lon), dtype=int)

        if incl_latlon:
            self.lon = np.zeros((n_lon))
            self.lat = np.zeros((n_lat))
            self.lon2d = np.zeros((n_lat, n_lon))
            self.lat2d = np.zeros((n_lat, n_lon))

        if incl_pftndx:
            self.hillslope_pftndx = np.zeros((ncolumns_per_gridcell, n_lat, n_lon), dtype=int)
        else:
            self.hillslope_pftndx = None

        # Placeholders for read-in data from each chunk
        self.chunk_mask = np.zeros((n_lat, n_lon))
        if recurse:
            self.this_chunk = HillslopeVars(
                ncolumns_per_gridcell,
                nhillslope,
                n_lat,
                n_lon,
                recurse=False,
                incl_latlon=incl_latlon,
            )

    def _read_samenames(self, chunk_ds):
        """
        Read hillslope variables from one chunk file with the same variable names
        in old and new versions

        Args:
            chunk_ds (xarray Dataset): Opened chunk file
        """
        self.this_chunk.nhillcolumns = chunk_ds.variables["nhillcolumns"][:,].astype(int)
        self.this_chunk.pct_hillslope = chunk_ds.variables["pct_hillslope"][:,]
        self.this_chunk.hillslope_index = chunk_ds.variables["hillslope_index"][:,].astype(int)
        self.this_chunk.column_index = chunk_ds.variables["column_index"][:,].astype(int)
        self.this_chunk.downhill_column_index = chunk_ds.variables["downhill_column_index"][
            :,
        ].astype(int)

    def _read_oldnames(self, chunk_file, read_bedrock, read_stream):
        """Read hillslope variables from one chunk file with old variable names

        Args:
            cfile (str): Path to chunk file
            read_bedrock (logical): Whether to read bedrock variable(s)
            read_stream (logical): Whether to read stream variable(s)
        """
        chunk_ds = Dataset(chunk_file, "r")
        self.chunk_mask = chunk_ds.variables["chunk_mask"][:]
        self.this_chunk.h_elev = chunk_ds.variables["h_height"][:]
        self.this_chunk.h_dist = chunk_ds.variables["h_length"][:]
        self.this_chunk.h_width = chunk_ds.variables["h_width"][:]
        self.this_chunk.h_area = chunk_ds.variables["h_area"][:]
        self.this_chunk.h_slope = chunk_ds.variables["h_slope"][:]
        self.this_chunk.h_aspect = chunk_ds.variables["h_aspect"][:]
        if read_bedrock:
            self.this_chunk.h_bedrock = chunk_ds.variables["h_bedrock"][:]
        if read_stream:
            self.this_chunk.h_stream_depth = chunk_ds.variables["h_stream_depth"][:]
            self.this_chunk.h_stream_width = chunk_ds.variables["h_stream_width"][:]
            self.h_stream_slope = chunk_ds.variables["hillslope_stream_slope"][:]

        self._read_samenames(chunk_ds)
        chunk_ds.close()

    def read(self, chunk_file, read_bedrock, read_stream, incl_latlon=False):
        """Read hillslope variables from one chunk file

        Args:
            cfile (str): Path to chunk file
            read_bedrock (logical): Whether to read bedrock variable(s)
            read_stream (logical): Whether to read stream variable(s)
        """

        chunk_ds = Dataset(chunk_file, "r")
        self.this_chunk.chunk_mask = chunk_ds.variables["chunk_mask"][:]

        if incl_latlon:
            self.this_chunk.lon = chunk_ds.variables["longitude"][:]
            self.this_chunk.lat = chunk_ds.variables["latitude"][:]
            self.this_chunk.lon2d = chunk_ds.variables["LONGXY"][:]
            self.this_chunk.lat2d = chunk_ds.variables["LATIXY"][:]

        try:
            self.this_chunk.h_elev = chunk_ds.variables["hillslope_elevation"][:]
        except KeyError:
            chunk_ds.close()
            self._read_oldnames(chunk_file, read_bedrock, read_stream)
            return
        self.this_chunk.h_dist = chunk_ds.variables["hillslope_distance"][:]
        self.this_chunk.h_width = chunk_ds.variables["hillslope_width"][:]
        self.this_chunk.h_area = chunk_ds.variables["hillslope_area"][:]
        self.this_chunk.h_slope = chunk_ds.variables["hillslope_slope"][:]
        self.this_chunk.h_aspect = chunk_ds.variables["hillslope_aspect"][:]
        if read_bedrock:
            self.this_chunk.h_bedrock = chunk_ds.variables["hillslope_bedrock_depth"][:]
        if read_stream:
            self.this_chunk.h_stream_depth = chunk_ds.variables["hillslope_stream_depth"][:]
            self.this_chunk.h_stream_width = chunk_ds.variables["hillslope_stream_width"][:]
            self.this_chunk.h_stream_slope = chunk_ds.variables["hillslope_stream_slope"][:]

        self._read_samenames(chunk_ds)
        chunk_ds.close()

    def update(
        self,
        i,
        j,
        add_bedrock,
        add_stream,
        *,
        landmask=None,
        incl_latlon=False,
        incl_chunkmask=False,
        this_chunk_1d=False,
        remove_if_too_few_aspects=True,
    ):
        """
        Update a gridcell in chunk
        """
        if landmask is not None and not (
            self.this_chunk.chunk_mask[j, i] > 0 and landmask[j, i] > 0
        ):
            return

        if incl_latlon:
            self.lon[i] = self.this_chunk.lon
            self.lat[j] = self.this_chunk.lat
            self.lon2d[j, i] = self.this_chunk.lon2d
            self.lat2d[j, i] = self.this_chunk.lat2d

        if incl_chunkmask:
            self.chunk_mask[j, i] = np.int32(self.this_chunk.chunk_mask)

        if this_chunk_1d:
            i_this = 0
            j_this = 0
        else:
            i_this = i
            j_this = j

        self.h_elev[:, j, i] = self.this_chunk.h_elev[:, j_this, i_this]
        self.h_dist[:, j, i] = self.this_chunk.h_dist[:, j_this, i_this]
        self.h_width[:, j, i] = self.this_chunk.h_width[:, j_this, i_this]
        self.h_area[:, j, i] = self.this_chunk.h_area[:, j_this, i_this]
        self.h_slope[:, j, i] = self.this_chunk.h_slope[:, j_this, i_this]
        self.h_aspect[:, j, i] = self.this_chunk.h_aspect[:, j_this, i_this]
        if add_bedrock:
            self.h_bedrock[:, j, i] = self.this_chunk.h_bedrock[:, j_this, i_this]
        if add_stream:
            self.h_stream_depth[j, i] = self.this_chunk.h_stream_depth[j_this, i_this]
            self.h_stream_width[j, i] = self.this_chunk.h_stream_width[j_this, i_this]
            self.h_stream_slope[j, i] = self.this_chunk.h_stream_slope[j_this, i_this]

        self.nhillcolumns[j, i] = self.this_chunk.nhillcolumns[j_this, i_this]
        self.pct_hillslope[:, j, i] = self.this_chunk.pct_hillslope[:, j_this, i_this]
        self.hillslope_index[:, j, i] = self.this_chunk.hillslope_index[:, j_this, i_this]
        self.column_index[:, j, i] = self.this_chunk.column_index[:, j_this, i_this]
        self.downhill_column_index[:, j, i] = self.this_chunk.downhill_column_index[
            :, j_this, i_this
        ]

        # if 2 or less valid aspects, remove all hillslope data
        if remove_if_too_few_aspects and self.this_chunk.nhillcolumns[j, i] > 0:
            # check number of hillslopes
            h_ndx = self.this_chunk.hillslope_index[:, j, i]
            nactual_hillslopes = np.unique(h_ndx[h_ndx > 0]).size
            if nactual_hillslopes < 3:
                self._remove_hillslope_data(i, j)

    def _remove_hillslope_data(self, i, j):
        self.h_elev[:, j, i] = 0
        self.h_dist[:, j, i] = 0
        self.h_width[:, j, i] = 0
        self.h_area[:, j, i] = 0
        self.h_slope[:, j, i] = 0
        self.h_aspect[:, j, i] = 0
        self.h_bedrock[:, j, i] = 0
        self.h_stream_depth[j, i] = 0
        self.h_stream_width[j, i] = 0
        self.h_stream_slope[j, i] = 0

        self.nhillcolumns[j, i] = 0
        self.pct_hillslope[:, j, i] = 0
        self.hillslope_index[:, j, i] = 0
        self.column_index[:, j, i] = 0
        self.downhill_column_index[:, j, i] = 0

    def save(
        self,
        *,
        input_file,
        output_file,
        ncolumns_per_gridcell,
        nhillslope,
        add_bedrock,
        add_stream,
        logger=None,
        n_lon=None,
        n_lat=None,
        incl_latlon=False,
        incl_chunkmask=False,
        save_fsurdat=False,
    ):
        """
        Save to netCDF
        """
        msg = f"Saving to {output_file}"
        if logger is None:
            print(msg)
        else:
            logger.info(msg)

        # Create and open file
        if save_fsurdat:
            shutil.copyfile(input_file, output_file)
            ds_out = Dataset(output_file, "a")
        else:
            ds_out = Dataset(output_file, "w", format=NETCDF_FORMAT)

        # Create dimensions
        if not save_fsurdat:
            if n_lon is None or n_lat is None:
                if n_lon != n_lat:
                    raise NotImplementedError("Unhandled: Only one of n_lon and n_lat being None")
                with Dataset(input_file, "r") as fsurdat:
                    n_lat = len(fsurdat.dimensions["lsmlat"])
                    n_lon = len(fsurdat.dimensions["lsmlon"])
            ds_out.createDimension("lsmlon", n_lon)
            ds_out.createDimension("lsmlat", n_lat)
        ds_out.createDimension("nhillslope", nhillslope)
        ds_out.createDimension("nmaxhillcol", ncolumns_per_gridcell)

        if incl_latlon:
            add_variable_nc(
                name="longitude",
                units="degrees",
                long_name="longitude - 1d",
                data=self.lon,
                dataset=ds_out,
                dims=["lsmlon"],
            )
            add_variable_nc(
                name="latitude",
                units="degrees",
                long_name="latitude - 1d",
                data=self.lat,
                dataset=ds_out,
                dims=["lsmlat"],
            )
            add_longxy_latixy_nc(self.lon2d, self.lat2d, ds_out)

        if incl_chunkmask:
            add_variable_nc(
                name="chunk_mask",
                units="unitless",
                long_name="chunk mask",
                data=self.chunk_mask,
                dataset=ds_out,
                dims=["lsmlat", "lsmlon"],
                data_type=np.int32,
            )

        add_variable_nc(
            name="hillslope_elevation",
            units="m",
            long_name="hillslope elevation above channel",
            data=self.h_elev,
            dataset=ds_out,
        )

        add_variable_nc(
            name="hillslope_distance",
            units="m",
            long_name="hillslope distance from channel",
            data=self.h_dist,
            dataset=ds_out,
        )

        add_variable_nc(
            name="hillslope_width",
            units="m",
            long_name="hillslope width",
            data=self.h_width,
            dataset=ds_out,
        )

        add_variable_nc(
            name="hillslope_area",
            units="m2",
            long_name="hillslope area",
            data=self.h_area,
            dataset=ds_out,
        )

        add_variable_nc(
            name="hillslope_slope",
            units="m/m",
            long_name="hillslope slope",
            data=self.h_slope,
            dataset=ds_out,
        )

        add_variable_nc(
            name="hillslope_aspect",
            units="radians",
            long_name="hillslope aspect (clockwise from North)",
            data=self.h_aspect,
            dataset=ds_out,
        )

        if add_bedrock:
            add_variable_nc(
                name="hillslope_bedrock_depth",
                units="meters",
                long_name="hillslope bedrock depth",
                data=self.h_bedrock,
                dataset=ds_out,
            )

        if add_stream:
            add_variable_nc(
                name="hillslope_stream_depth",
                units="meters",
                long_name="stream channel bankfull depth",
                data=self.h_stream_depth,
                dataset=ds_out,
                dims=["lsmlat", "lsmlon"],
            )
            add_variable_nc(
                name="hillslope_stream_width",
                units="meters",
                long_name="stream channel bankfull width",
                data=self.h_stream_width,
                dataset=ds_out,
                dims=["lsmlat", "lsmlon"],
            )
            add_variable_nc(
                name="hillslope_stream_slope",
                units="m/m",
                long_name="stream channel slope",
                data=self.h_stream_slope,
                dataset=ds_out,
                dims=["lsmlat", "lsmlon"],
            )

        add_variable_nc(
            name="nhillcolumns",
            units="unitless",
            long_name="number of columns per landunit",
            data_type=np.int32,
            data=self.nhillcolumns.astype(np.int32),
            dataset=ds_out,
            dims=["lsmlat", "lsmlon"],
        )

        add_variable_nc(
            name="pct_hillslope",
            units="per cent",
            long_name="percent hillslope of landunit",
            data=self.pct_hillslope,
            dataset=ds_out,
            dims=["nhillslope", "lsmlat", "lsmlon"],
        )

        add_variable_nc(
            name="hillslope_index",
            units="unitless",
            long_name="hillslope_index",
            data=self.hillslope_index.astype(np.int32),
            data_type=np.int32,
            dataset=ds_out,
        )

        add_variable_nc(
            name="column_index",
            units="unitless",
            long_name="column index",
            data=self.column_index.astype(np.int32),
            data_type=np.int32,
            dataset=ds_out,
        )

        add_variable_nc(
            name="downhill_column_index",
            units="unitless",
            long_name="downhill column index",
            data=self.downhill_column_index.astype(np.int32),
            dataset=ds_out,
            data_type=np.int32,
        )

        if self.hillslope_pftndx is not None:
            add_variable_nc(
                name="hillslope_pftndx",
                units="unitless",
                long_name="hillslope pft indices",
                data=self.hillslope_pftndx.astype(np.int32),
                dataset=ds_out,
                data_type=np.int32,
            )

        # Save
        ds_out.close()
