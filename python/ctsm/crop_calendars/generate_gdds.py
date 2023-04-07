# Import the CTSM Python utilities, functions for GDD generation
import utils
import generate_gdds_functions as gddfn
paramfile_dir = "/glade/p/cesmdata/cseg/inputdata/lnd/clm2/paramdata"
    
# Import other shared functions
import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 
import cropcal_module as cc
from cropcal_figs_module import *

# Import everything else
import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.transforms import Bbox
import warnings
import cartopy.crs as ccrs
import datetime as dt
import pickle
import datetime as dt
import argparse
import logging

# Figure settings
plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

# Info re: PFT parameter set
my_clm_ver = 51
my_clm_subver = "c211112"

# Suppress some warnings
import warnings
warnings.filterwarnings("ignore", message="__len__ for multi-part geometries is deprecated and will be removed in Shapely 2.0. Check the length of the `geoms` property instead to get the  number of parts of a multi-part geometry.")
warnings.filterwarnings("ignore", message="Iteration over multi-part geometries is deprecated and will be removed in Shapely 2.0. Use the `geoms` property to access the constituent parts of a multi-part geometry.")


def get_multicrop_maps(ds, theseVars, crop_fracs_yx, dummy_fill, gdd_units):
    
    # Get GDDs for these crops
    da_eachCFT = xr.concat((ds[x] for i, x in enumerate(theseVars)),
                            dim="cft")
    if "time" in ds.dims:
        da_eachCFT = da_eachCFT.isel(time=0, drop=True)
    da_eachCFT = da_eachCFT.where(da_eachCFT != dummy_fill)
    da_eachCFT.attrs['units'] = gdd_units
    
    # What are the maximum differences seen between different crop types?
    if len(theseVars) > 1:
        maxDiff = np.nanmax(da_eachCFT.max(dim="cft") - da_eachCFT.min(dim="cft"))
        if maxDiff > 0:
            print(f"   Max difference among crop types: {np.round(maxDiff)}")
    
    if crop_fracs_yx is None:
        return da_eachCFT.isel(cft=0, drop=True)
    
    # Warn if GDD is NaN anywhere that there is area
    da_eachCFT['cft'] = crop_fracs_yx['cft']
    gddNaN_areaPos = np.isnan(da_eachCFT) & (crop_fracs_yx > 0)
    if np.any(gddNaN_areaPos):
        total_bad_croparea = np.nansum(crop_fracs_yx.where(gddNaN_areaPos).values)
        total_croparea = np.nansum(crop_fracs_yx.values)
        print(f"   GDD reqt NaN but area positive ({np.round(total_bad_croparea/total_croparea*100, 1)}% of this crop's area)")
    
    # Get areas and weights, masking cell-crops with NaN GDDs
    crop_fracs_yx = crop_fracs_yx.where(~np.isnan(da_eachCFT))
    crop_area_yx = crop_fracs_yx.sum(dim="cft")
    weights_yx = crop_fracs_yx / crop_area_yx
    weights_sum_gt0 = weights_yx.sum(dim='cft').where(weights_yx > 0)
    assert(np.isclose(np.nanmin(weights_sum_gt0.values), 1.0))
    assert(np.isclose(np.nanmax(weights_sum_gt0.values), 1.0))
    
    # Mask GDDs and weights where there is no area
    da_eachCFT = da_eachCFT.where(crop_fracs_yx > 0)
    if len(theseVars)==1:
        return da_eachCFT.isel(cft=0, drop=True)
    weights_yx = weights_yx.where(crop_fracs_yx > 0)
    weights_sum = weights_yx.sum(dim='cft').where(crop_area_yx > 0)
    assert(np.isclose(np.nanmin(weights_sum.values), 1.0))
    assert(np.isclose(np.nanmax(weights_sum.values), 1.0))
    
    # Ensure grid match between GDDs and weights
    if not np.array_equal(da_eachCFT['lon'].values, weights_yx['lon'].values):
        raise RuntimeError("lon mismatch")
    if not np.array_equal(da_eachCFT['lat'].values, weights_yx['lat'].values):
        raise RuntimeError("lat mismatch")
    
    # Get area-weighted mean GDD requirements for all crops
    da = (da_eachCFT * weights_yx).sum(dim="cft")
    da.attrs['units'] = gdd_units
    da = da.where(crop_area_yx > 0)
    
    # Ensure that weighted mean is between each cell's min and max
    whereBad = (da < da_eachCFT.min(dim="cft")) | (da > da_eachCFT.max(dim="cft"))
    if np.any(whereBad):
        where_belowMin = da.where(da < da_eachCFT.min(dim="cft"))
        worst_belowMin = np.min((da_eachCFT.min(dim='cft') - where_belowMin).values)
        where_aboveMax = da.where(da > da_eachCFT.max(dim="cft"))
        worst_aboveMax = np.max((where_aboveMax - da_eachCFT.max(dim='cft')).values)
        worst = max(worst_belowMin, worst_aboveMax)
        tol = 1e-12
        if worst > 1e-12:
            raise RuntimeError(f"Some value is outside expected range by {worst} (exceeds tolerance {tol})")
    
    return da


def main(run_dir=None, first_season=None, last_season=None, sdates_file=None, hdates_file=None, output_dir=None, save_figs=True, only_make_figs=False, run1_name=None, run2_name=None, land_use_file=None, first_land_use_year=None, last_land_use_year=None, unlimited_season_length=False, logger=None):
    
    # Directories to save output files and figures
    if not output_dir:
        if only_make_figs:
            output_dir = run_dir
        else:
            output_dir = os.path.join(run_dir, "generate_gdds")
            if not unlimited_season_length:
                output_dir += ".mxmat"
            output_dir += "." + dt.datetime.now().strftime('%Y-%m-%d-%H%M%S')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    outdir_figs = os.path.join(output_dir, "figs")

    # Set up log file and function, if needed
    if logger is None:
        logging.basicConfig(level=logging.DEBUG,
                            format="",
                            filename=os.path.join(output_dir, 'generate_gdds.log'),
                            filemode='a')
        logger = logging.getLogger('')

    # Print some info
    gddfn.log(logger, f"Saving to {output_dir}")
    
    
    ##########################
    ### Import and process ###
    ##########################
    
    if not only_make_figs:
    
        # Keep 1 extra year to avoid incomplete final growing season for crops harvested after Dec. 31.
        y1_import_str = f"{first_season+1}-01-01"
        yN_import_str = f"{last_season+2}-01-01"

        gddfn.log(logger, f"Importing netCDF time steps {y1_import_str} through {yN_import_str} (years are +1 because of CTSM output naming)")

        pickle_file = os.path.join(output_dir, f'{first_season}-{last_season}.pickle')
        h1_ds_file = os.path.join(output_dir, f'{first_season}-{last_season}.h1_ds.nc')
        if os.path.exists(pickle_file):
            with open(pickle_file, 'rb') as f:
                first_season, last_season, pickle_year, gddaccum_yp_list, gddharv_yp_list, skip_patches_for_isel_nan_lastyear, lastYear_active_patch_indices_list, incorrectly_daily, gddharv_in_h3, save_figs, incl_vegtypes_str, incl_patches1d_itype_veg, mxsowings = pickle.load(f)
            print(f'Will resume import at {pickle_year+1}')
            h1_ds = None
        else:
            incorrectly_daily = False
            skip_patches_for_isel_nan_lastyear = np.ndarray([])
            gddharv_in_h3 = False
            pickle_year = -np.inf
            gddaccum_yp_list = []
            gddharv_yp_list = []
            incl_vegtypes_str = None
            lastYear_active_patch_indices_list = None
        sdates_rx = sdates_file
        hdates_rx = hdates_file
        
        if not unlimited_season_length:
            mxmats = cc.import_max_gs_length(paramfile_dir, my_clm_ver, my_clm_subver)
        else:
            mxmats = None
        
        for y, thisYear in enumerate(np.arange(first_season+1,last_season+3)):
            
            if thisYear <= pickle_year:
                continue
            
            h1_ds, sdates_rx, hdates_rx, gddaccum_yp_list, gddharv_yp_list, skip_patches_for_isel_nan_lastyear, lastYear_active_patch_indices_list, incorrectly_daily, gddharv_in_h3, incl_vegtypes_str, incl_patches1d_itype_veg, mxsowings = gddfn.import_and_process_1yr(first_season, last_season, y, thisYear, sdates_rx, hdates_rx, gddaccum_yp_list, gddharv_yp_list, skip_patches_for_isel_nan_lastyear, lastYear_active_patch_indices_list, incorrectly_daily, gddharv_in_h3, run_dir, incl_vegtypes_str, h1_ds_file, mxmats, cc.get_gs_len_da, logger)
            
            gddfn.log(logger, f'   Saving pickle file ({pickle_file})...')
            with open(pickle_file, 'wb') as f:
                pickle.dump([first_season, last_season, thisYear, gddaccum_yp_list, gddharv_yp_list, skip_patches_for_isel_nan_lastyear, lastYear_active_patch_indices_list, incorrectly_daily, gddharv_in_h3, save_figs, incl_vegtypes_str, incl_patches1d_itype_veg, mxsowings], f, protocol=-1)
                
        
        if isinstance(incl_vegtypes_str, list):
            incl_vegtypes_str = np.array(incl_vegtypes_str)
        plot_vegtypes_str = incl_vegtypes_str[[i for i,c in enumerate(gddaccum_yp_list) if not isinstance(c,type(None))]]
        
        gddfn.log(logger, "Done")
        
        if not h1_ds:
            h1_ds = xr.open_dataset(h1_ds_file)
    
    
    ######################################################
    ### Get and grid mean GDDs in GGCMI growing season ###
    ######################################################
    
    if not only_make_figs:
    
        longname_prefix = "GDD harvest target for "
        
        # Could skip this by saving sdates_rx['time_bounds']
        sdates_rx = gddfn.import_rx_dates("s", sdates_rx, incl_patches1d_itype_veg, mxsowings, logger)
        
        gddfn.log(logger, 'Getting and gridding mean GDDs...')
        gdd_maps_ds = gddfn.yp_list_to_ds(gddaccum_yp_list, h1_ds, incl_vegtypes_str, sdates_rx, longname_prefix, logger)
        gddharv_maps_ds = gddfn.yp_list_to_ds(gddharv_yp_list, h1_ds, incl_vegtypes_str, sdates_rx, longname_prefix, logger)
        
        # Fill NAs with dummy values
        dummy_fill = -1
        gdd_fill0_maps_ds = gdd_maps_ds.fillna(0)
        gdd_maps_ds = gdd_maps_ds.fillna(dummy_fill)
        gddfn.log(logger, 'Done getting and gridding means.')
        
        # Add dummy variables for crops not actually simulated
        gddfn.log(logger, "Adding dummy variables...")
        # Unnecessary?
        template_ds = xr.open_dataset(sdates_file, decode_times=True)
        all_vars = [v.replace("sdate","gdd") for v in template_ds if "sdate" in v]
        all_longnames = [template_ds[v].attrs["long_name"].replace("Planting day ", longname_prefix) + " (dummy)" for v in template_ds if "sdate" in v]
        dummy_vars = []
        dummy_longnames = []
        for v, thisVar in enumerate(all_vars):
            if thisVar not in gdd_maps_ds:
                dummy_vars.append(thisVar)
                dummy_longnames.append(all_longnames[v])
        
        def make_dummy(thisCrop_gridded, addend):
            dummy_gridded = thisCrop_gridded
            dummy_gridded.values = dummy_gridded.values*0 + addend
            return dummy_gridded
        for v in gdd_maps_ds:
            thisCrop_gridded = gdd_maps_ds[v].copy()
            thisCrop_fill0_gridded = gdd_fill0_maps_ds[v].copy()
            break
        dummy_gridded = make_dummy(thisCrop_gridded, -1)
        dummy_gridded0 = make_dummy(thisCrop_fill0_gridded, 0)
        
        for v, thisVar in enumerate(dummy_vars):
            if thisVar in gdd_maps_ds:
                gddfn.error(logger, f'{thisVar} is already in gdd_maps_ds. Why overwrite it with dummy?')
            dummy_gridded.name = thisVar
            dummy_gridded.attrs["long_name"] = dummy_longnames[v]
            gdd_maps_ds[thisVar] = dummy_gridded
            dummy_gridded0.name = thisVar
            dummy_gridded0.attrs["long_name"] = dummy_longnames[v]
            gdd_fill0_maps_ds[thisVar] = dummy_gridded0
        
        # Add lon/lat attributes
        def add_lonlat_attrs(ds):
            ds.lon.attrs = {\
                "long_name": "coordinate_longitude",
                "units": "degrees_east"}
            ds.lat.attrs = {\
                "long_name": "coordinate_latitude",
                "units": "degrees_north"}
            return ds
        gdd_maps_ds = add_lonlat_attrs(gdd_maps_ds)
        gdd_fill0_maps_ds = add_lonlat_attrs(gdd_fill0_maps_ds)
        gddharv_maps_ds = add_lonlat_attrs(gddharv_maps_ds)
        
        gddfn.log(logger, "Done.")
    
    
    ######################
    ### Save to netCDF ###
    ######################
    
    if not only_make_figs:
        gddfn.log(logger, "Saving...")
        
        # Get output file path
        datestr = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        outfile = os.path.join(output_dir, "gdds_" + datestr + ".nc")
        outfile_fill0 = os.path.join(output_dir, "gdds_fill0_" + datestr + ".nc")
        
        def save_gdds(sdates_file, hdates_file, outfile, gdd_maps_ds, sdates_rx):
            # Set up output file from template (i.e., prescribed sowing dates).
            template_ds = xr.open_dataset(sdates_file, decode_times=True)
            for v in template_ds:
                if "sdate" in v:
                    template_ds = template_ds.drop(v)
            template_ds.to_netcdf(path=outfile, format="NETCDF3_CLASSIC")
            template_ds.close()
        
            # Add global attributes
            comment = f"Derived from CLM run plus crop calendar input files {os.path.basename(sdates_file) and {os.path.basename(hdates_file)}}."
            gdd_maps_ds.attrs = {\
                "author": "Sam Rabin (sam.rabin@gmail.com)",
                "comment": comment,
                "created": dt.datetime.now().astimezone().isoformat()
                }
        
            # Add time_bounds
            gdd_maps_ds["time_bounds"] = sdates_rx.time_bounds
        
            # Save cultivar GDDs
            gdd_maps_ds.to_netcdf(outfile, mode="w", format="NETCDF3_CLASSIC")
        
        save_gdds(sdates_file, hdates_file, outfile, gdd_maps_ds, sdates_rx)
        save_gdds(sdates_file, hdates_file, outfile_fill0, gdd_fill0_maps_ds, sdates_rx)
        
        gddfn.log(logger, "Done saving.")
    
    
    ########################################
    ### Save things needed for mapmaking ###
    ########################################
    
    def add_attrs_to_map_ds(map_ds, incl_vegtypes_str, dummy_fill, outdir_figs, first_season, last_season):
        return map_ds.assign_attrs({'incl_vegtypes_str': incl_vegtypes_str,
                                    'dummy_fill': dummy_fill,
                                    'outdir_figs': outdir_figs,
                                    'y1': first_season,
                                    'yN': last_season})
    
    if not only_make_figs:
        if not os.path.exists(outdir_figs):
            os.makedirs(outdir_figs)

        gdd_maps_ds = add_attrs_to_map_ds(gdd_maps_ds, plot_vegtypes_str, dummy_fill, outdir_figs, first_season, last_season)
        gddharv_maps_ds = add_attrs_to_map_ds(gddharv_maps_ds, plot_vegtypes_str, dummy_fill, outdir_figs, first_season, last_season)
        
        gdd_maps_ds.to_netcdf(os.path.join(outdir_figs, "gdd_maps.nc"))
        gddharv_maps_ds.to_netcdf(os.path.join(outdir_figs, "gddharv_maps.nc"))
    
    
    #################################################
    ### Save before/after map and boxplot figures ###
    #################################################
    
    def get_bounds_ncolors(gdd_spacing, diff_map_yx):
        vmax = np.floor(np.nanmax(diff_map_yx.values)/gdd_spacing)*gdd_spacing
        vmin = -vmax
        epsilon = np.nextafter(0, 1)
        bounds = list(np.arange(vmin, vmax, gdd_spacing)) + [vmax-epsilon]
        if 0 in bounds:
            bounds.remove(0)
            bounds[bounds.index(-gdd_spacing)] /= 2
            bounds[bounds.index(gdd_spacing)] /= 2
        Ncolors = len(bounds) + 1
        return vmax, bounds, Ncolors    
    
    def make_map(ax, this_map, this_title, vmax, bin_width, fontsize_ticklabels, fontsize_titles, bounds=None, extend='both', cmap=None, cbar_ticks=None, vmin=None):
        
        if bounds:
            if not cmap:
                raise RuntimeError("Calling make_map() with bounds requires cmap to be specified")
            norm = mcolors.BoundaryNorm(bounds, cmap.N, extend=extend)
            im1 = ax.pcolormesh(this_map.lon.values, this_map.lat.values,
                                this_map, shading="auto",
                                norm=norm,
                                cmap=cmap)
        else:
            if np.any(this_map.values < 0):
                gdd_spacing = 500
                vmax = np.floor(np.nanmax(this_map.values)/gdd_spacing)*gdd_spacing
                if vmin is not None:
                    raise RuntimeError("Do not specify vmin in this call of make_map()")
                vmin = -vmax
                Ncolors = vmax/gdd_spacing
                if Ncolors % 2 == 0: Ncolors += 1
                if not cmap:
                    cmap = cm.get_cmap(cropcal_colors['div_other_nonnorm'], Ncolors)
                
                if np.any(this_map.values > vmax) and np.any(this_map.values < vmin):
                    extend = 'both'
                elif np.any(this_map.values > vmax):
                    extend = 'max'
                elif np.any(this_map.values < vmin):
                    extend = 'min'
                else:
                    extend = 'neither'
                
            else:
                if vmin is None:
                    vmin = 0
                else:
                    vmin = np.floor(vmin/500)*500
                vmax = np.floor(vmax/500)*500
                Ncolors = int(vmax/500)
                if not cmap:
                    cmap=cm.get_cmap(cropcal_colors['seq_other'], Ncolors+1)
                extend = 'max'
                extend_color = cmap.colors[-1]
                cmap = mcolors.ListedColormap(cmap.colors[:Ncolors])
                cmap.set_over(extend_color)
                
            im1 = ax.pcolormesh(this_map.lon.values, this_map.lat.values, 
                    this_map, shading="auto",
                    vmin=vmin, vmax=vmax,
                    cmap=cmap)
            
        ax.set_extent([-180,180,-63,90],crs=ccrs.PlateCarree())
        ax.coastlines(linewidth=0.3)
        ax.set_title(this_title, fontsize=fontsize_titles, fontweight="bold", y=0.96)
        cbar = plt.colorbar(im1, orientation="horizontal", fraction=0.1, pad=0.02,
                            aspect=40, extend=extend, spacing='proportional')
        cbar.ax.tick_params(labelsize=fontsize_ticklabels)
        cbar.ax.set_xlabel(this_map.attrs['units'],
                           fontsize=fontsize_ticklabels)
        cbar.ax.xaxis.set_label_coords(x=0.115, y=2.6)
        if cbar_ticks:
            cbar.ax.set_xticks(cbar_ticks)
        
        ticks = np.arange(-60, 91, bin_width)
        ticklabels = [str(x) for x in ticks]
        for i,x in enumerate(ticks):
            if x%2:
                ticklabels[i] = ''
        plt.yticks(np.arange(-60,91,15), labels=ticklabels,
                   fontsize=fontsize_ticklabels)
        plt.axis('off')
        
    def get_non_nans(in_da, fillValue):
        in_da = in_da.where(in_da != fillValue)
        return in_da.values[~np.isnan(in_da.values)]
    
    linewidth = 1.5
    def set_boxplot_props(bp, color, linewidth):
        linewidth = linewidth
        plt.setp(bp['boxes'], color=color, linewidth=linewidth)
        plt.setp(bp['whiskers'], color=color, linewidth=linewidth)
        plt.setp(bp['caps'], color=color, linewidth=linewidth)
        plt.setp(bp['medians'], color=color, linewidth=linewidth)
        plt.setp(bp['fliers'], markeredgecolor=color, markersize=6, linewidth=linewidth, markeredgewidth=linewidth/2)
    
    def make_plot(data, offset, linewidth):
        offset = 0.4*offset
        bpl = plt.boxplot(data, positions=np.array(range(len(data)))*2.0+offset, widths=0.6, 
                          boxprops=dict(linewidth=linewidth), whiskerprops=dict(linewidth=linewidth), 
                          capprops=dict(linewidth=linewidth), medianprops=dict(linewidth=linewidth),
                          flierprops=dict(markeredgewidth=0.5))
        return bpl
    
    def make_figures(first_land_use_year, last_land_use_year, land_use_file, run1_name, run2_name, thisDir=None, gdd_maps_ds=None, gddharv_maps_ds=None, outdir_figs=None, linewidth=1.5):
        if not gdd_maps_ds:
            if not thisDir:
                gddfn.error(logger, 'If not providing gdd_maps_ds, you must provide thisDir (location of gdd_maps.nc)')
            gdd_maps_ds = xr.open_dataset(thisDir + 'gdd_maps.nc')
        if not gddharv_maps_ds:
            if not thisDir:
                gddfn.error(logger, 'If not providing gddharv_maps_ds, you must provide thisDir (location of gddharv_maps.nc)')
            gddharv_maps_ds = xr.open_dataset(thisDir + 'gdd_maps.nc')
    
        # Get info
        incl_vegtypes_str = gdd_maps_ds.attrs['incl_vegtypes_str']
        dummy_fill = gdd_maps_ds.attrs['dummy_fill']
        if not outdir_figs:
            outdir_figs = gdd_maps_ds.attrs['outdir_figs']
        try:
            y1 = gdd_maps_ds.attrs['y1']
            yN = gdd_maps_ds.attrs['yN']
        # Backwards compatibility with a bug (fixed 2023-01-03)
        except:
            y1 = gdd_maps_ds.attrs['first_season']
            yN = gdd_maps_ds.attrs['last_season']
        # Import LU data, if doing so
        if land_use_file:
            y1_lu = y1 if first_land_use_year == None else first_land_use_year
            yN_lu = yN if last_land_use_year == None else last_land_use_year
            lu_ds = cc.open_lu_ds(land_use_file, y1_lu, yN_lu, gdd_maps_ds, ungrid=False)
            lu_years_text = f" (masked by {y1_lu}-{yN_lu} area)"
            lu_years_file = f"_mask{y1_lu}-{yN_lu}"
        else:
            lu_ds = None
    
        # layout = "3x1"
        # layout = "2x2"
        layout = "3x2"
        bin_width = 15
        lat_bin_edges = np.arange(0, 91, bin_width)
    
        fontsize_titles = 12
        fontsize_axislabels = 12
        fontsize_ticklabels = 12
    
        Nbins = len(lat_bin_edges)-1
        bin_names = ["All"]
        for b in np.arange(Nbins):
            lower = lat_bin_edges[b]
            upper = lat_bin_edges[b+1]
            bin_names.append(f"{lower}–{upper}")
        
        color_old = cropcal_colors_cases(run1_name)
        if color_old is None:
            color_old = '#beaed4'
        color_new = cropcal_colors_cases(run2_name)
        if color_new is None:
            color_new = '#7fc97f'
        gdd_units = 'GDD (°C • day)'
    
        # Maps
        ny = 3
        nx = 1
        gddfn.log(logger, "Making before/after maps...")
        for v, vegtype_str in enumerate(incl_vegtypes_str + ["Corn", "Cotton", "Rice", "Soybean", "Sugarcane", "Wheat"]):
            print(f"{vegtype_str}...")
            
            # Get component types
            if vegtype_str in incl_vegtypes_str:
                vegtypes_str = [vegtype_str]
            elif not lu_ds:
                raise RuntimeError(f"If mapping {vegtype_str}, you must provide land use dataset")
            else:
                vegtypes_str = [x for x in incl_vegtypes_str if vegtype_str.lower() in x]
            vegtypes_int = [utils.vegtype_str2int(x)[0] for x in vegtypes_str]
            
            # Crop fraction map (for masking and weighting)
            if lu_ds:
                crop_fracs_yx = (lu_ds.LANDFRAC_PFT * lu_ds.PCT_CROP * lu_ds.PCT_CFT.sel(cft=vegtypes_int)).sum(dim="time")
                if np.sum(crop_fracs_yx) == 0:
                    print(f"Skipping {vegtype_str} (no area)")
                    continue
            else:
                crop_fracs_yx = None

            theseVars = [f"gdd1_{x}" for x in vegtypes_int]
            gddharv_map_yx = get_multicrop_maps(gddharv_maps_ds, theseVars, crop_fracs_yx, dummy_fill, gdd_units)
            gdd_map_yx = get_multicrop_maps(gdd_maps_ds, theseVars, crop_fracs_yx, dummy_fill, gdd_units)
            
            # Get figure title
            if len(vegtypes_str) > 1:
                vegtype_str_title = vegtype_str
            else:                
                vegtype_str_title = vegtype_str.replace("_", " ")
                if "irrigated" not in vegtype_str:
                    vegtype_str_title = "rainfed " + vegtype_str_title
                vegtype_str_title = vegtype_str_title.capitalize()
                    
            vmin = min(np.min(gdd_map_yx), np.min(gddharv_map_yx)).values
            vmax = max(np.max(gdd_map_yx), np.max(gddharv_map_yx)).values
            
            # Set up figure and first subplot
            if layout == "3x1":
                fig = plt.figure(figsize=(7.5,14))
                ax = fig.add_subplot(ny,nx,1,projection=ccrs.PlateCarree())
            elif layout == "2x2":
                fig = plt.figure(figsize=(12,6))
                spec = fig.add_gridspec(nrows=2, ncols=2,
                                        width_ratios=[0.4,0.6])
                ax = fig.add_subplot(spec[0,0],projection=ccrs.PlateCarree())
            elif layout == "3x2":
                fig = plt.figure(figsize=(14,9))
                spec = fig.add_gridspec(nrows=3, ncols=2,
                                        width_ratios=[0.5,0.5],
                                        wspace=0.2)
                ax = fig.add_subplot(spec[0,0],projection=ccrs.PlateCarree())
            else:
                gddfn.error(logger, f"layout {layout} not recognized")
            
            thisMin = int(np.round(np.nanmin(gddharv_map_yx)))
            thisMax = int(np.round(np.nanmax(gddharv_map_yx)))
            thisTitle = f"{run1_name} (range {thisMin}–{thisMax})"
            make_map(ax, gddharv_map_yx, thisTitle, vmax, bin_width,
                     fontsize_ticklabels, fontsize_titles, vmin=vmin)
            
            if layout == "3x1":
                ax = fig.add_subplot(ny,nx,2,projection=ccrs.PlateCarree())
            elif layout in ["2x2", "3x2"]:
                ax = fig.add_subplot(spec[1,0],projection=ccrs.PlateCarree())
            else:
                gddfn.error(logger, f"layout {layout} not recognized")
            thisMin = int(np.round(np.nanmin(gdd_map_yx)))
            thisMax = int(np.round(np.nanmax(gdd_map_yx)))
            thisTitle = f"{run2_name} (range {thisMin}–{thisMax})"
            make_map(ax, gdd_map_yx, thisTitle, vmax, bin_width,
                     fontsize_ticklabels, fontsize_titles, vmin=vmin)
            
            # Difference
            if layout == "3x2":
                ax = fig.add_subplot(spec[2,0],projection=ccrs.PlateCarree())
                thisMin = int(np.round(np.nanmin(gdd_map_yx)))
                thisMax = int(np.round(np.nanmax(gdd_map_yx)))
                thisTitle = f"{run2_name} minus {run1_name}"
                diff_map_yx = gdd_map_yx - gddharv_map_yx
                diff_map_yx.attrs['units'] = gdd_units
                
                gdd_spacing = 500
                vmax, bounds, Ncolors = get_bounds_ncolors(gdd_spacing, diff_map_yx)
                if Ncolors < 9:
                    gdd_spacing = 250
                    vmax, bounds, Ncolors = get_bounds_ncolors(gdd_spacing, diff_map_yx)
                
                cmap = cm.get_cmap(cropcal_colors['div_other_nonnorm'], Ncolors)
                cbar_ticks = []
                include_0bin_ticks = Ncolors <= 13
                if vmax <= 3000:
                    tick_spacing = gdd_spacing*2
                elif vmax <= 5000:
                    tick_spacing = 1500
                else:
                    tick_spacing = 2000
                previous = -np.inf
                for x in bounds:
                    if (not include_0bin_ticks) and (x>0) and (previous<0):
                        cbar_ticks.append(0)
                    if x % tick_spacing == 0 or (include_0bin_ticks and abs(x)==gdd_spacing/2):
                        cbar_ticks.append(x)
                    previous = x
                
                make_map(ax, diff_map_yx, thisTitle, vmax, bin_width,
                        fontsize_ticklabels, fontsize_titles, bounds=bounds,
                        extend='both', cmap=cmap, cbar_ticks=cbar_ticks)
            
            # Boxplots #####################

            gdd_vector = get_non_nans(gdd_map_yx, dummy_fill)
            gddharv_vector = get_non_nans(gddharv_map_yx, dummy_fill)
            
            lat_abs = np.abs(gdd_map_yx.lat.values)
            gdd_bybin_old = [gddharv_vector]
            gdd_bybin_new = [gdd_vector]
            for b in np.arange(Nbins):
                lower = lat_bin_edges[b]
                upper = lat_bin_edges[b+1]
                lat_inds = np.where((lat_abs>=lower) & (lat_abs<upper))[0]
                gdd_vector_thisBin = get_non_nans(gdd_map_yx[lat_inds,:], dummy_fill)
                gddharv_vector_thisBin = get_non_nans(gddharv_map_yx[lat_inds,:], dummy_fill)
                gdd_bybin_old.append(gddharv_vector_thisBin)
                gdd_bybin_new.append(gdd_vector_thisBin)
                    
            if layout == "3x1":
                ax = fig.add_subplot(ny,nx,3)
            elif layout in ["2x2", "3x2"]:
                ax = fig.add_subplot(spec[:,1])
            else:
                gddfn.error(logger, f"layout {layout} not recognized")
    
            # Shift bottom of plot up to make room for legend
            ax_pos = ax.get_position()
            ax.set_position(Bbox.from_extents(ax_pos.x0, 0.19, ax_pos.x1, ax_pos.y1))
            # Define legend position
            legend_bbox_to_anchor = (0, -0.15, 1, 0.2)
            
            bpl = make_plot(gdd_bybin_old, -1, linewidth)
            bpr = make_plot(gdd_bybin_new, 1, linewidth)
            set_boxplot_props(bpl, color_old, linewidth)
            set_boxplot_props(bpr, color_new, linewidth)
            
            # draw temporary lines to create a legend
            plt.plot([], c=color_old, label=run1_name, linewidth=linewidth)
            plt.plot([], c=color_new, label=run2_name, linewidth=linewidth)
            plt.legend(fontsize=fontsize_titles,
                       bbox_to_anchor=legend_bbox_to_anchor,
                       ncol = 2,
                       loc='lower left',
                       mode = 'expand')
            
            plt.xticks(range(0, len(bin_names) * 2, 2), bin_names,
                       fontsize=fontsize_ticklabels)
            plt.yticks(fontsize=fontsize_ticklabels)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            
            plt.xlabel("Latitude zone (absolute value)", fontsize=fontsize_axislabels)
            plt.ylabel(gdd_units, fontsize=fontsize_axislabels)
            ax.yaxis.set_label_coords(-0.11, 0.5)
            plt.title(f"Zonal changes", fontsize=fontsize_titles, fontweight="bold")
            
            plt.suptitle(f"Maturity requirements: {vegtype_str_title}" + lu_years_text,
                         fontsize=fontsize_titles*1.2,
                         fontweight="bold",
                         y=0.95)
            
            if vegtype_str in incl_vegtypes_str:
                outfile = os.path.join(outdir_figs, f"{theseVars[0]}_{vegtype_str}_gs{y1}-{yN}{lu_years_file}.png")
            else:
                outfile = os.path.join(outdir_figs, f"{vegtype_str}_gs{y1}-{yN}{lu_years_file}.png")
            plt.savefig(outfile, dpi=300, transparent=False, facecolor='white',
                        bbox_inches='tight')
            plt.close()
    
        gddfn.log(logger, "Done.")
    
    if save_figs: 
        if only_make_figs:
            gdd_maps_ds = xr.open_dataset(os.path.join(run_dir, "figs", "gdd_maps.nc"))
            gddharv_maps_ds = xr.open_dataset(os.path.join(run_dir, "figs", "gddharv_maps.nc"))
        make_figures(first_land_use_year, last_land_use_year, land_use_file, run1_name, run2_name, gdd_maps_ds=gdd_maps_ds, gddharv_maps_ds=gddharv_maps_ds, outdir_figs=outdir_figs, linewidth=linewidth)


if __name__ == "__main__":
    
    ###############################
    ### Process input arguments ###
    ###############################
    
    # Set arguments
    parser = argparse.ArgumentParser(description="ADD DESCRIPTION HERE")
    parser.add_argument("-r", "--run-dir", 
                        help="Directory where run outputs can be found (and where outputs will go). If --only-make-figs, this is the directory with the preprocessed files (e.g., *.pickle file).",
                        required=True)
    parser.add_argument("-1", "--first-season", 
                        help="First growing season to include in calculation of mean",
                        required=True,
                        type=int)
    parser.add_argument("-n", "-N", "--last-season", 
                        help="Last growing season to include in calculation of mean",
                        required=True,
                        type=int)
    parser.add_argument("-o", "--output-dir",
                        help="Output directory. Default is auto-generated subdir of -r/--run-dir.")
    parser.add_argument("-sd", "--sdates-file", 
                        help="File of prescribed sowing dates",
                        required=True)
    parser.add_argument("-hd", "--hdates-file", 
                        help="File of prescribed harvest dates",
                        required=True)
    figsgroup = parser.add_mutually_exclusive_group()
    figsgroup.add_argument("--dont-save-figs", 
                           help="Do not save figures",
                           action="store_true", default=False)
    figsgroup.add_argument("--only-make-figs", 
                           help="Use preprocessed files to make figures only",
                           action="store_true", default=False)
    parser.add_argument("--run1-name", 
                        help="Name of original values to show in figures",
                        default="Old")
    parser.add_argument("--run2-name", 
                        help="Name of new values to show in figures",
                        default="New")
    parser.add_argument("-lu", "--land-use-file",
                        help="Path to CLM land use timeseries file, for masking figures",
                        default=None)
    parser.add_argument("--first-land-use-year",
                        help="First year in land use file to use for masking. Default --first-season.",
                        default=None,
                        type=int)
    parser.add_argument("--last-land-use-year",
                        help="Last year in land use file to use for masking. Default --last-season.",
                        default=None,
                        type=int)
    parser.add_argument("--unlimited-season-length", 
                        help="Limit mean growing season length based on CLM CFT parameter mxmat.",
                        action="store_true", default=False)
    
    # Get arguments
    args = parser.parse_args(sys.argv[1:])
    for k, v in sorted(vars(args).items()):
        print(f"{k}: {v}")
    save_figs = not args.dont_save_figs
    
    # Call main()
    main(run_dir=args.run_dir, first_season=args.first_season, last_season=args.last_season, sdates_file=args.sdates_file, hdates_file=args.hdates_file, output_dir=args.output_dir, save_figs=save_figs, only_make_figs=args.only_make_figs, run1_name=args.run1_name, run2_name=args.run2_name, land_use_file=args.land_use_file, first_land_use_year=args.first_land_use_year, last_land_use_year=args.last_land_use_year, unlimited_season_length=args.unlimited_season_length)

# main(run_dir="/Users/Shared/CESM_runs/tests_10x15_20230329_gddgen/202303301820",
#      sdates_file="/Users/Shared/CESM_work/crop_dates_mostrice/sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f10_f10_mg37.2000-2000.20230330_165301.nc",
#      hdates_file="/Users/Shared/CESM_work/crop_dates_mostrice/hdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-f10_f10_mg37.2000-2000.20230330_165301.nc",
#      first_season=1997, last_season=2003,
#      save_figs=False)