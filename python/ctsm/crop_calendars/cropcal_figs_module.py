import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
from matplotlib import cm
import matplotlib.collections as mplcol
import cartopy.feature as cfeature
import xarray as xr
import cftime
from scipy import stats

# Import the CTSM Python utilities
import utils

# Colormaps (maps)
cropcal_colors = {
    'seq_timeofyear': 'twilight_shifted',
    'seq_other': 'plasma', # magma_r? CMRmap_r?
    'div_yieldirr': 'BrBG',
    'div_timeofyear': 'twilight_shifted',
    'div_other_nonnorm': 'PuOr_r',
    'div_other_norm': 'RdBu_r',
    'underlay': [0.75, 0.75, 0.75, 1],
    'underlay_lighter': [0.85, 0.85, 0.85, 1],
    'underlay_lightest': [0.92, 0.92, 0.92, 1],
}

# Cases (line and scatter plots)
def cropcal_colors_cases(casename):
    case_color_dict = {
        'clm default': [x/255 for x in [92, 219, 219]],
        'prescribed calendars': [x/255 for x in [250, 102, 240]],
        'prescribed maturity': [x/255 for x in [128,0,0]],
        'prescribed sowing': [x/255 for x in [133, 92, 255]],
    }
    case_color_dict['5.0 lu'] = case_color_dict['clm default']
    case_color_dict['5.2 lu'] = case_color_dict['prescribed calendars']
    
    case_color = None
    casename_for_colors = casename.lower().replace(" (0)", "").replace(" (1)", "")
    if casename_for_colors in case_color_dict:
        case_color = case_color_dict[casename_for_colors]
    return case_color


def chunk_colorbar(this_map, cbar_spacing, cmap, crop, fontsize, pct_absdiffs_masked_before, sumdiff_beforemask, varInfo, vmin, vmax, posNeg=False, underlay=None, v=0):
    # Make a temporary plot with the same color axis and colorbar settings we would use in make_map().
    plt.pcolormesh(this_map, vmin=vmin, vmax=vmax)
    cb0 = plt.colorbar(location="bottom")
    cb0.ax.tick_params(labelsize=fontsize['ticklabels'])
    
    # Where did plt.colorbar() draw bin boundaries? These are referred to as "tick marks," but note that the extreme values might lie outside [vmin, vmax].
    ticklocations = cb0.get_ticks()
    bounds = ticklocations
    
    # In our plot, we will move vmin left and vmax right to ensure that the tick marks are the color bin boundaries.
    if cb0.vmin < bounds[0]:
        raise RuntimeError("Handle vmin < bounds[0]")
    elif cb0.vmax > bounds[-1]:
        raise RuntimeError("Handle vmax > bounds[-1]")
    elif 0 not in bounds:
        raise RuntimeError("Handle 0 not in bounds")
    vmin = bounds[0]
    vmax = bounds[-1]
    
    # Get number of color bins
    Nbins = len(bounds) - 1
    bottom_of_topbin = bounds[-2]
    bottom_of_2ndbin = bounds[-3]
    binwidth = bounds[-1] - bounds[-2]
    if Nbins < 8:
        Nbins *= 2
        bottom_of_2ndbin = bottom_of_topbin
        binwidth /= 2
        bottom_of_topbin += binwidth
    
    # Ensure that most extreme bin (on at least one side of 0) has at least one gridcell included. If not, remove the most extreme bins and check again.
    maxinmap = np.nanmax(np.abs(this_map.values))
    if maxinmap < bottom_of_topbin:
        if maxinmap < bottom_of_2ndbin:
            raise RuntimeError("How is maxinmap less than the bottom of the SECOND bin??")
        vmax -= binwidth
        vmin += binwidth
        Nbins-=2
        if ticklocations[0] < vmin:
            ticklocations = ticklocations[1:]
        if ticklocations[-1] > vmax:
            ticklocations = ticklocations[:-1]
        
    # Get new colormap with the right number of bins.
    this_cmap = cm.get_cmap(cmap, Nbins)
    if Nbins % 2:
        raise RuntimeError(f"Expected even number of color bins; got {Nbins}")
    
    # Special color for small-masked cells
    if 'maskcolorbar_near0' in varInfo and varInfo['maskcolorbar_near0'][v] is not None:
        
        # Get near-zero threshold
        nearzero_thresh = varInfo['maskcolorbar_near0'][v]
        if isinstance(nearzero_thresh, str):
            nearzero_parts = nearzero_thresh.split("|")
            if nearzero_parts[0] != "percentile":
                raise RuntimeError(f"Unable to parse maskcolorbar_near0 value {nearzero_thresh}")
            
            # The input value should be a percentile, 0-100
            nearzero_val_in = np.float(nearzero_parts[1])
            
            if len(nearzero_parts) > 2:
                if nearzero_parts[2] != "cumulative":
                    raise RuntimeError(f"Unable to parse maskcolorbar_near0 value {nearzero_thresh}")
                frac_to_include = 1 - nearzero_val_in/100
                nearzero_thresh = get_threshold_lowestpercentile_cumulative(this_map, frac_to_include)
            else:
                nearzero_thresh = np.nanpercentile(np.abs(this_map), nearzero_val_in)
        else:
            raise RuntimeError("Can only parse string values of maskcolorbar_near0")
        
        # Add near-zero bin
        bounds, cbar_spacing, pct_absdiffs_masked_before, ticklabels, ticklocations, vmax, vmin = maskcolorbar_near0(binwidth, bounds, cbar_spacing, crop, nearzero_thresh, pct_absdiffs_masked_before, posNeg, sumdiff_beforemask, this_map, ticklocations, vmax, vmin)
        
        # Add color for that bin
        if underlay is not None:
            raise RuntimeError("You need a different color to distinguish maskcolorbar_near0 cells from other-masked cells")
        this_cmap = get_ListedColormap(cmap, this_cmap, Nbins)
        if posNeg:
            if crop == "Crops decreasing":
                new_colors = np.concatenate((this_cmap.colors[:int(Nbins/2)],
                                                np.array([cropcal_colors['underlay']])),
                                            axis=0)
            elif crop == "Crops increasing":
                new_colors = np.concatenate((np.array([cropcal_colors['underlay']]),
                                                this_cmap.colors[int(Nbins/2)+1:]),
                                            axis=0)
            else:
                raise RuntimeError(f"posNeg: Crop {crop} not recognized for color bar (2)")
        else:
            new_colors = np.concatenate((this_cmap.colors[:int(Nbins/2)],
                                            np.array([cropcal_colors['underlay']]),
                                            this_cmap.colors[int(Nbins/2)+1:]),
                                        axis=0)
        
        this_cmap = mcolors.ListedColormap(new_colors)
    
    # Remove our temporary plot and its colorbar.
    plt.cla()
    cb0.remove()
    
    return bounds, cbar_spacing, pct_absdiffs_masked_before, this_cmap, ticklabels, ticklocations, vmin, vmax


def get_amount_masked(crop, this_map_timemean, sumdiff_beforemask, pct_absdiffs_masked_before, reason):
    sumdiff_aftermask = np.nansum(np.abs(this_map_timemean.values))
    pct_absdiffs_masked = 100 * (1 - sumdiff_aftermask / sumdiff_beforemask)
    pct_absdiffs_masked_here = pct_absdiffs_masked - pct_absdiffs_masked_before
    print(f"   Masked {crop} ({reason}): {round(pct_absdiffs_masked_here, 1)}%")
    pct_absdiffs_masked_before = pct_absdiffs_masked
    return pct_absdiffs_masked_before


def get_ListedColormap(cmap_name, this_cmap, Nbins):
    if isinstance(this_cmap, mcolors.LinearSegmentedColormap):
        this_cmap = cm.get_cmap(cmap_name)
        color_list = [this_cmap(x) for x in np.arange(0, 1, 1/Nbins)]
        color_list = []
        for i, x in enumerate(np.arange(0, 1+1e-9, 1/Nbins)):
            color_list.append(this_cmap(x))
            if i>0 and color_list[i] == color_list[i-1]:
                print(f"{prev_x} → color_list[{i-1}] = {color_list[i-1]}")
                print(f"{x} → color_list[{i}] = {color_list[i]}")
                raise RuntimeError("Repeated color!")
            prev_x = x
        this_cmap = mcolors.ListedColormap(color_list)
    elif not isinstance(this_cmap, mcolors.ListedColormap):
        raise RuntimeError(f"Not sure how to get list of colors from {type(this_cmap)}")
    
    return this_cmap


def get_threshold_lowestpercentile_cumulative(this_map_timemean, frac_to_include):
    flattened = np.abs(this_map_timemean.values).flatten()
    flattened_is_ok = np.where(~np.isnan(flattened))
    okflattened = flattened[flattened_is_ok]
    oksorted = np.flip(np.sort(okflattened))
    okcumsum = np.cumsum(oksorted)
    okcumprop = okcumsum / np.sum(oksorted)
    for i, x in enumerate(okcumprop):
        if x >= frac_to_include:
            break
    lowest_threshold = oksorted[i]
    return lowest_threshold


def get_non_rx_map(var_info, cases, casename, this_var, thisCrop_main, found_types, plot_y1, plot_yN, ref_casename):
    time_dim = var_info['time_dim']
    case = cases[casename]
    
    # Trim to included years
    try:
        this_ds = case['ds'].sel({time_dim: slice(plot_y1, plot_yN)})
    except:
        # Try converting years to cftime
        plot_y1 = cftime.DatetimeNoLeap(plot_y1, 1, 1)
        plot_yN = cftime.DatetimeNoLeap(plot_yN, 1, 1)
        this_ds = case['ds'].sel({time_dim: slice(plot_y1, plot_yN)})
    
    if this_var not in case['ds']:
        return xr.DataArray(), "continue"
    elif ref_casename and ref_casename!="rx" and cases[ref_casename]['res'] != case['res']:
        # Not bothering with regridding (for now?)
        return xr.DataArray(), "continue"
    this_map = this_ds[this_var]
    
    # Prepare to mask out patch-years with no area
    if "gs" in this_map.dims:
        croparea_ever_positive = this_ds['croparea_positive_wholeseason'].sum(dim="gs")
    elif "time" in this_map.dims or this_var in ["QIRRIG_DEMAND_PATCH_PKMTH"]:
        croparea_ever_positive = this_ds['croparea_positive_sowing'].sum(dim="time")
    else:
        raise RuntimeError(f"Unsure how to mask patch-years with no area for {this_var} with dims {this_map.dims}")
    this_ds['croparea_ever_positive'] = croparea_ever_positive
    
    # Grid the included vegetation types, if needed
    if "lon" not in this_map.dims:
        this_map = utils.grid_one_variable(this_ds, this_var, vegtype=found_types)
        croparea_ever_positive = utils.grid_one_variable(this_ds, 'croparea_ever_positive', vegtype=found_types) > 0
    # If not, select the included vegetation types
    else:
        this_map = this_map.sel(ivt_str=found_types)
        croparea_ever_positive = this_ds['croparea_ever_positive'].sel(ivt_str=found_types) > 0
        
    return this_map, croparea_ever_positive, time_dim


def make_map(ax, this_map, fontsize, bounds=None, cbar=None, cbar_labelpad=4.0, cbar_max=None, cbar_spacing='uniform', cmap=cropcal_colors['seq_other'], extend_bounds='both', extend_nonbounds='both', linewidth=1.0, lonlat_bin_width=None, show_cbar=False, subplot_label=None, this_title=None, ticklabels=None, ticklocations=None, underlay=None, underlay_color=None, units=None, vmax=None, vmin=None, vrange=None):
    
    
    if underlay is not None:
        if underlay_color is None:
            underlay_color = cropcal_colors['underlay']
        underlay_cmap = mcolors.ListedColormap(np.array([underlay_color, [1, 1, 1, 1]]))
        ax.pcolormesh(underlay.lon.values, underlay.lat.values,
                      underlay, cmap=underlay_cmap)
    
    if bounds is not None:
        norm = mcolors.BoundaryNorm(bounds, cmap.N, extend=extend_bounds)
        im = ax.pcolormesh(this_map.lon.values, this_map.lat.values,
                           this_map, shading="auto",
                           norm=norm,
                           cmap=cmap)
    else:
        im = ax.pcolormesh(this_map.lon.values, this_map.lat.values, 
                           this_map, shading="auto",
                           cmap=cmap,
                           vmin=vmin, vmax=vmax)
        if vrange:
            im.set_clim(vrange[0], vrange[1])
    ax.set_extent([-180,180,-63,90],crs=ccrs.PlateCarree())
    
    if subplot_label is not None:
        plt.text(0, 0.95, f"({subplot_label})", transform=ax.transAxes,
             fontsize=fontsize['axislabels'])
    
    # # Country borders
    # ax.add_feature(cfeature.BORDERS, linewidth=linewidth, edgecolor="white", alpha=0.5)
    # ax.add_feature(cfeature.BORDERS, linewidth=linewidth*0.6, alpha=0.3)
    
    # Coastlines
    ax.coastlines(linewidth=linewidth, color="white", alpha=0.5)
    ax.coastlines(linewidth=linewidth*0.6, alpha=0.3)
    
    if this_title:
        ax.set_title(this_title, fontsize=fontsize['titles'])
    if show_cbar:
        if cbar:
            cbar.remove()
        
        if bounds is not None:
            cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='horizontal', fraction=0.1, pad=0.02, spacing=cbar_spacing)
        else:
            cbar = plt.colorbar(im, ax=ax, orientation="horizontal", fraction=0.1, pad=0.02, extend=extend_nonbounds, spacing=cbar_spacing)
        
        deal_with_ticklabels(cbar, cbar_max, ticklabels, ticklocations, units, im)
        cbar.set_label(label=units, fontsize=fontsize['axislabels'], verticalalignment="center", labelpad=cbar_labelpad)
        cbar.ax.tick_params(labelsize=fontsize['ticklabels'])
        if units is not None and "month" in units.lower():
            cbar.ax.tick_params(length=0)
    
    
    if lonlat_bin_width:
        set_ticks(lonlat_bin_width, fontsize, "y")
        # set_ticks(lonlat_bin_width, fontsize, "x")
    else:
        # Need to do this for subplot row labels
        set_ticks(-1, fontsize, "y")
        plt.yticks([])
    for x in ax.spines:
        ax.spines[x].set_visible(False)
    
    if show_cbar:
        return im, cbar
    else:
        return im, None


# Note that we're not actually masking here; we're just setting a special color for a particular part of the colorbar
def maskcolorbar_near0(binwidth, bounds, cbar_spacing, crop, nearzero_thresh, pct_absdiffs_masked_before, posNeg, sumdiff_beforemask, this_map_timemean, ticklocations, vmax, vmin):
    ticklabels = None
    
    # Setting this to False will cause some weird results that I'm not going to try and fix for now.
    # For example:
        # Floating-point errors like 0.6000000000001
        # When there are 2 color bins between near-zero bin and first labeled tick,
        #    there should be an extra tick+label between those 2
    use_proportional_cbar = True
    if use_proportional_cbar:
        cbar_spacing = "proportional"
    else:
        cbar_spacing = "uniform"
            
    # Add near-zero bin
    if posNeg:
        if crop == "Crops decreasing":
            vmax = 0
            ticklocations = ticklocations[np.where(ticklocations <= 0)]
            bounds = np.concatenate((np.arange(vmin, -binwidth+1e-9, binwidth),
                                        np.array([-nearzero_thresh, 0])))
            if not use_proportional_cbar:
                ticklabels = np.concatenate(([str(x) for x in ticklocations[ticklocations<0]],
                                             ["-{:.2g}".format(nearzero_thresh)],
                                             ["0"]))
                ticklocations = np.concatenate((ticklocations[ticklocations < 0], [-nearzero_thresh, 0]))
        elif crop == "Crops increasing":
            vmin = 0
            ticklocations = ticklocations[np.where(ticklocations >= 0)]
            bounds = np.concatenate((np.array([0, nearzero_thresh]),
                                        np.arange(binwidth, vmax+1e-9, binwidth)))
            if not use_proportional_cbar:
                ticklabels = np.concatenate((["0"],
                                             ["{:.2g}".format(nearzero_thresh)],
                                             [str(x) for x in ticklocations[ticklocations>0]],))
                ticklocations = np.concatenate(([0, nearzero_thresh], ticklocations[ticklocations > 0]))
        else:
            raise RuntimeError(f"posNeg: Crop {crop} not recognized for color bar")
    else:
        bounds = np.concatenate((np.arange(vmin, -binwidth+1e-9, binwidth),
                                    np.array([-nearzero_thresh, nearzero_thresh]),
                                    np.arange(binwidth, vmax+1e-9, binwidth)))
        if not use_proportional_cbar:
            ticklabels = np.concatenate(([str(x) for x in ticklocations[ticklocations<0]],
                                         ["±{:.2g}".format(nearzero_thresh)],
                                         [str(x) for x in ticklocations[ticklocations>0]]))
    
    # Get stats
    this_map_timemean_fake = this_map_timemean.where(np.abs(this_map_timemean) >= nearzero_thresh)
    
    # Diagnostics
    pct_absdiffs_masked_before = get_amount_masked(crop, this_map_timemean_fake, sumdiff_beforemask, pct_absdiffs_masked_before, f"nearest to zero, threshold {nearzero_thresh}")
    
    # Because this is just a fake mask we must not do any additional masking, real or fake, after this. If we do, then our diagnostics for that will be messed up. To ensure we don't try to do any subsequent masking, set this to None, which should throw an error in get_amount_masked().
    pct_absdiffs_masked_before = None
    
    return bounds, cbar_spacing, pct_absdiffs_masked_before, ticklabels, ticklocations, vmax, vmin


def deal_with_ticklabels(cbar, cbar_max, ticklabels, ticklocations, units, im):
    if ticklocations is not None:
        cbar.set_ticks(ticklocations)
        if units is not None and units.lower() == "month":
            cbar.set_ticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
            units == "Month"
        elif ticklabels is not None:
            cbar.set_ticklabels(ticklabels)
    if isinstance(im, mplcol.QuadMesh):
        clim_max = im.get_clim()[1]
    else:
        clim_max = im
    if cbar_max is not None and clim_max > cbar_max:
        if ticklabels is not None:
            raise RuntimeError("How to handle this now that you are specifying ticklocations separate from ticklabels?")
        ticks = cbar.get_ticks()
        if ticks[-2] > cbar_max:
            raise RuntimeError(f"Specified cbar_max is {cbar_max} but highest bin BEGINS at {ticks[-2]}")
        ticklabels = ticks.copy()
        ticklabels[-1] = cbar_max
        for i, x in enumerate(ticklabels):
            if x == int(x):
                ticklabels[i] = str(int(x))
        cbar.set_ticks(ticks) # Calling this before set_xticklabels() avoids "UserWarning: FixedFormatter should only be used together with FixedLocator" (https://stackoverflow.com/questions/63723514/userwarning-fixedformatter-should-only-be-used-together-with-fixedlocator)
        cbar.set_ticklabels(ticklabels)


def set_ticks(lonlat_bin_width, fontsize, x_or_y):
    if x_or_y == "x":
        ticks = np.arange(-180, 181, lonlat_bin_width)
    else:
        ticks = np.arange(-60, 91, lonlat_bin_width)
        
    ticklabels = [str(x) for x in ticks]
    for i,x in enumerate(ticks):
        if x%2:
            ticklabels[i] = ''
    
    if x_or_y == "x":
        plt.xticks(ticks, labels=ticklabels,
                    fontsize=fontsize['ticklabels'])
    else:
        plt.yticks(ticks, labels=ticklabels,
                    fontsize=fontsize['ticklabels'])
