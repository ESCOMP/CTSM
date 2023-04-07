import re
import cf_units as cf
import cftime
import numpy as np
import xarray as xr
from cartopy.util import add_cyclic_point

#from xr_ds_ex import xr_ds_ex

# generate annual means, weighted by days / month
def weighted_annual_mean(array, time_in='time', time_out='time'):
    
    if isinstance(array[time_in].values[0], cftime.datetime):
        month_length = array[time_in].dt.days_in_month
        
        # After https://docs.xarray.dev/en/v0.5.1/examples/monthly-means.html
        group = f'{time_in}.year'
        weights = month_length.groupby(group) / month_length.groupby(group).sum()
        np.testing.assert_allclose(weights.groupby(group).sum().values, 1)
        array = (array * weights).groupby(group).sum(dim=time_in, skipna=True)
        if time_out != "year":
            array = array.rename({'year': time_out})
        
    else:
        mon_day  = xr.DataArray(np.array([31,28,31,30,31,30,31,31,30,31,30,31]), dims=['month'])
        mon_wgt  = mon_day/mon_day.sum()
        array = (array.rolling({time_in: 12}, center=False) # rolling
                .construct("month") # construct the array
                .isel({time_in: slice(11, None, 12)}) # slice so that the first element is [1..12], second is [13..24]
                .dot(mon_wgt, dims=["month"]))
        if time_in != time_out:
            array = array.rename({time_in: time_out})
        
    return array

def change_units(ds, variable_str, variable_bounds_str, target_unit_str):
    """ Applies unit conversion on an xarray DataArray """
    orig_units = cf.Unit(ds[variable_str].attrs["units"])
    target_units = cf.Unit(target_unit_str)
    variable_in_new_units = xr.apply_ufunc(
        orig_units.convert,
        ds[variable_bounds_str],
        target_units,
        output_dtypes=[ds[variable_bounds_str].dtype],
    )
    return variable_in_new_units

def clean_units(units):
    """replace some troublesome unit terms with acceptable replacements"""
    replacements = {'kgC':'kg', 'gC':'g', 'gC13':'g', 'gC14':'g', 'gN':'g',
                    'unitless':'1',
                    'years':'common_years', 'yr':'common_year',
                    'meq':'mmol', 'neq':'nmol'}
    units_split = re.split('( |\(|\)|\^|\*|/|-[0-9]+|[0-9]+)', units)
    units_split_repl = \
        [replacements[token] if token in replacements else token for token in units_split]
    return ''.join(units_split_repl)

def copy_fill_settings(da_in, da_out):
    """
    propagate _FillValue and missing_value settings from da_in to da_out
    return da_out
    """
    if '_FillValue' in da_in.encoding:
        da_out.encoding['_FillValue'] = da_in.encoding['_FillValue']
    else:
        da_out.encoding['_FillValue'] = None
    if 'missing_value' in da_in.encoding:
        da_out.attrs['missing_value'] = da_in.encoding['missing_value']
    return da_out

def dim_cnt_check(ds, varname, dim_cnt):
    """confirm that varname in ds has dim_cnt dimensions"""
    if len(ds[varname].dims) != dim_cnt:
        msg_full = 'unexpected dim_cnt=%d, varname=%s' % (len(ds[varname].dims), varname)
        raise ValueError(msg_full)

def time_set_mid(ds, time_name):
    """
    set ds[time_name] to midpoint of ds[time_name].attrs['bounds'], if bounds attribute exists
    type of ds[time_name] is not changed
    ds is returned
    """

    if 'bounds' not in ds[time_name].attrs:
        return ds

    # determine units and calendar of unencoded time values
    if ds[time_name].dtype == np.dtype('O'):
        units = 'days since 0000-01-01'
        calendar = 'noleap'
    else:
        units = ds[time_name].attrs['units']
        calendar = ds[time_name].attrs['calendar']

    # construct unencoded midpoint values, assumes bounds dim is 2nd
    tb_name = ds[time_name].attrs['bounds']
    if ds[tb_name].dtype == np.dtype('O'):
        tb_vals = cftime.date2num(ds[tb_name].values, units=units, calendar=calendar)
    else:
        tb_vals = ds[tb_name].values
    tb_mid = tb_vals.mean(axis=1)

    # set ds[time_name] to tb_mid
    if ds[time_name].dtype == np.dtype('O'):
        ds[time_name] = cftime.num2date(tb_mid, units=units, calendar=calendar)
    else:
        ds[time_name] = tb_mid

    return ds

def time_year_plus_frac(ds, time_name):
    """return time variable, as year plus fraction of year"""

    # this is straightforward if time has units='days since 0000-01-01' and calendar='noleap'
    # so convert specification of time to that representation

    # get time values as an np.ndarray of cftime objects
    if np.dtype(ds[time_name]) == np.dtype('O'):
        tvals_cftime = ds[time_name].values
    else:
        tvals_cftime = cftime.num2date(
            ds[time_name].values, ds[time_name].attrs['units'], ds[time_name].attrs['calendar'])

    # convert cftime objects to representation mentioned above
    tvals_days = cftime.date2num(tvals_cftime, 'days since 0000-01-01', calendar='noleap')

    return tvals_days / 365.0


# add cyclic point
def cyclic_dataarray(da, coord='lon'):
    """ Add a cyclic coordinate point to a DataArray along a specified
    named coordinate dimension.
    >>> from xray import DataArray
    >>> data = DataArray([[1, 2, 3], [4, 5, 6]],
    ...                      coords={'x': [1, 2], 'y': range(3)},
    ...                      dims=['x', 'y'])
    >>> cd = cyclic_dataarray(data, 'y')
    >>> print cd.data
    array([[1, 2, 3, 1],
           [4, 5, 6, 4]])
    """
    assert isinstance(da, xr.DataArray)

    lon_idx = da.dims.index(coord)
    cyclic_data, cyclic_coord = add_cyclic_point(da.values,
                                                 coord=da.coords[coord],
                                                 axis=lon_idx)

    # Copy and add the cyclic coordinate and data
    new_coords = dict(da.coords)
    new_coords[coord] = cyclic_coord
    new_values = cyclic_data

    new_da = xr.DataArray(new_values, dims=da.dims, coords=new_coords)

    # Copy the attributes for the re-constructed data and coords
    for att, val in da.attrs.items():
        new_da.attrs[att] = val
    for c in da.coords:
        for att in da.coords[c].attrs:
            new_da.coords[c].attrs[att] = da.coords[c].attrs[att]

    return new_da

# as above, but for a dataset
# doesn't work because dims are locked in a dataset
'''
def cyclic_dataset(ds, coord='lon'):
    assert isinstance(ds, xr.Dataset)

    lon_idx = ds.dims.index(coord)
    cyclic_data, cyclic_coord = add_cyclic_point(ds.values,
                                                 coord=ds.coords[coord],
                                                 axis=lon_idx)

    # Copy and add the cyclic coordinate and data
    new_coords = dict(ds.coords)
    new_coords[coord] = cyclic_coord
    new_values = cyclic_data

    new_ds = xr.DataSet(new_values, dims=ds.dims, coords=new_coords)

    # Copy the attributes for the re-constructed data and coords
    for att, val in ds.attrs.items():
        new_ds.attrs[att] = val
    for c in ds.coords:
        for att in ds.coords[c].attrs:
            new_ds.coords[c].attrs[att] = ds.coords[c].attrs[att]

    return new_ds
'''


# List of PFTs used in CLM
def define_pftlist():
    pftlist =  ["not_vegetated",
        "needleleaf_evergreen_temperate_tree",
        "needleleaf_evergreen_boreal_tree",
        "needleleaf_deciduous_boreal_tree",
        "broadleaf_evergreen_tropical_tree",
        "broadleaf_evergreen_temperate_tree",
        "broadleaf_deciduous_tropical_tree",
        "broadleaf_deciduous_temperate_tree",
        "broadleaf_deciduous_boreal_tree",
        "broadleaf_evergreen_shrub",
        "broadleaf_deciduous_temperate_shrub",
        "broadleaf_deciduous_boreal_shrub",
        "c3_arctic_grass",
        "c3_non-arctic_grass",
        "c4_grass",
        "unmanaged_c3_crop",
        "unmanaged_c3_irrigated",
        "temperate_corn",
        "irrigated_temperate_corn",
        "spring_wheat",
        "irrigated_spring_wheat",
        "winter_wheat",
        "irrigated_winter_wheat",
        "soybean",
        "irrigated_soybean",
        "barley",
        "irrigated_barley",
        "winter_barley",
        "irrigated_winter_barley",
        "rye",
        "irrigated_rye",
        "winter_rye",
        "irrigated_winter_rye",
        "cassava",
        "irrigated_cassava",
        "citrus",
        "irrigated_citrus",
        "cocoa",
        "irrigated_cocoa",
        "coffee",
        "irrigated_coffee",
        "cotton",
        "irrigated_cotton",
        "datepalm",
        "irrigated_datepalm",
        "foddergrass",
        "irrigated_foddergrass",
        "grapes",
        "irrigated_grapes",
        "groundnuts",
        "irrigated_groundnuts",
        "millet",
        "irrigated_millet",
        "oilpalm",
        "irrigated_oilpalm",
        "potatoes",
        "irrigated_potatoes",
        "pulses",
        "irrigated_pulses",
        "rapeseed",
        "irrigated_rapeseed",
        "rice",
        "irrigated_rice",
        "sorghum",
        "irrigated_sorghum",
        "sugarbeet",
        "irrigated_sugarbeet",
        "sugarcane",
        "irrigated_sugarcane",
        "sunflower",
        "irrigated_sunflower",
        "miscanthus",
        "irrigated_miscanthus",
        "switchgrass",
        "irrigated_switchgrass",
        "tropical_corn",
        "irrigated_tropical_corn",
        "tropical_soybean",
        "irrigated_tropical_soybean"]
    return pftlist


# Get CLM ivt number corresponding to a given name
def ivt_str2int(ivt_str):
    
    pftlist = define_pftlist()
    if isinstance(ivt_str, str):
        ivt_int = pftlist.index(ivt_str)
    elif isinstance(ivt_str, list) or isinstance(ivt_str, np.ndarray):
        ivt_int = [ivt_str2int(x) for x in ivt_str]
        if isinstance(ivt_str, np.ndarray):
            ivt_int = np.array(ivt_int)
    else:
        raise RuntimeError(f"Update ivt_str_to_int() to handle input of type {type(ivt_str)} (if possible)")
    
    return ivt_int

# Get CLM ivt name corresponding to a given number
def ivt_int2str(ivt_int):
    
    pftlist = define_pftlist()
    if np.issubdtype(type(ivt_int), np.integer) or int(ivt_int)==ivt_int:
        ivt_str = pftlist[int(ivt_int)]
    elif isinstance(ivt_int, list) or isinstance(ivt_int, np.ndarray):
        ivt_str = [ivt_int2str(x) for x in ivt_int]
        if isinstance(ivt_int, np.ndarray):
            ivt_str = np.array(ivt_str)
    elif isinstance(ivt_int, float):
        raise RuntimeError("List indices must be integers")
    else:
        raise RuntimeError(f"Update ivt_str_to_int() to handle input of type {type(ivt_int)} (if possible)")
    
    return ivt_str


# Does this vegetation type's name match (for a given comparison method) any member of a filtering list?
'''
Methods:
    ok_contains:    True if any member of this_filter is found in this_vegtype.
    notok_contains: True of no member of this_filter is found in this_vegtype.
    ok_exact:       True if this_vegtype matches any member of this_filter 
                    exactly.
    notok_exact:    True if this_vegtype does not match any member of 
                    this_filter exactly.
'''
def is_this_vegtype(this_vegtype, this_filter, this_method):

    # Make sure data type of this_vegtype is acceptable
    if isinstance(this_vegtype, float) and int(this_vegtype)==this_vegtype:
        this_vegtype = int(this_vegtype)
    data_type_ok = lambda x: isinstance(x, str) or isinstance(x, int) or isinstance(x, np.int64)
    ok_input = True
    if not data_type_ok(this_vegtype):
        if isinstance(this_vegtype, xr.core.dataarray.DataArray):
            this_vegtype = this_vegtype.values
        if isinstance(this_vegtype, (list, np.ndarray)):
            if len(this_vegtype) == 1 and data_type_ok(this_vegtype[0]):
                this_vegtype = this_vegtype[0]
            elif data_type_ok(this_vegtype[0]):
                raise TypeError("is_this_vegtype(): this_vegtype must be a single string or integer, not a list of them. Did you mean to call is_each_vegtype() instead?")
            else:
                ok_input = False
        else:
            ok_input = False
    if not ok_input:
        raise TypeError(f"is_this_vegtype(): First argument (this_vegtype) must be a string or integer, not {type(this_vegtype)}")
    
    # Make sure data type of this_filter is acceptable
    if not np.iterable(this_filter):
        raise TypeError(f"is_this_vegtype(): Second argument (this_filter) must be iterable (e.g., a list), not {type(this_filter)}")
    
    # Perform the comparison
    if this_method == "ok_contains":
        return any(n in this_vegtype for n in this_filter)
    elif this_method == "notok_contains":
        return not any(n in this_vegtype for n in this_filter)
    elif this_method == "ok_exact":
        return any(n == this_vegtype for n in this_filter)
    elif this_method == "notok_exact":
        return not any(n == this_vegtype for n in this_filter)
    else:
        raise ValueError(f"Unknown comparison method: '{this_method}'")


# Get boolean list of whether each vegetation type in list is a managed crop
'''
    this_vegtypelist: The list of vegetation types whose members you want to 
                      test.
    this_filter:      The list of strings against which you want to compare 
                      each member of this_vegtypelist.
    this_method:      How you want to do the comparison. See is_this_vegtype().
'''
def is_each_vegtype(this_vegtypelist, this_filter, this_method):
    
    if isinstance(this_vegtypelist, xr.DataArray):
        this_vegtypelist = this_vegtypelist.values
    
    return [is_this_vegtype(x, this_filter, this_method) for x in this_vegtypelist]

# Helper function to check that a list is strictly increasing
def is_strictly_increasing(L):
    # https://stackoverflow.com/a/4983359/2965321
    return all(x<y for x, y in zip(L, L[1:]))

# Ensure that longitude axis coordinates are monotonically increasing
def make_lon_increasing(xr_obj):
    if not "lon" in xr_obj.dims:
        return xr_obj
    
    lons = xr_obj.lon.values
    if is_strictly_increasing(lons):
        return xr_obj
    
    shift = 0
    while not is_strictly_increasing(lons) and shift < lons.size:
        shift = shift + 1
        lons = np.roll(lons, 1, axis=0)
    if not is_strictly_increasing(lons):
        raise RuntimeError("Unable to rearrange longitude axis so it's monotonically increasing")
        
    return xr_obj.roll(lon=shift, roll_coords=True)


# Convert a longitude axis that's -180 to 180 around the international date line to one that's 0 to 360 around the prime meridian. If you pass in a Dataset or DataArray, the "lon" coordinates will be changed. Otherwise, it assumes you're passing in numeric data.
def lon_idl2pm(lons_in, fail_silently = False):
    def check_ok(tmp, fail_silently):
        msg = ""
        
        if np.any(tmp > 180):
            msg = f"Maximum longitude is already > 180 ({np.max(tmp)})"
        elif np.any(tmp < -180):
            msg = f"Minimum longitude is < -180 ({np.min(tmp)})"
            
        if msg == "":
            return True
        elif fail_silently:
            return False
        else:
            raise ValueError(msg)
            
    def do_it(tmp):
        tmp = tmp + 360
        tmp = np.mod(tmp, 360)
        return tmp
    
    if isinstance(lons_in, (xr.DataArray, xr.Dataset)):
        if not check_ok(lons_in.lon.values, fail_silently):
            return lons_in
        lons_out = lons_in
        lons_out = lons_out.assign_coords(lon=do_it(lons_in.lon.values))
        lons_out = make_lon_increasing(lons_out)
    else:
        if not check_ok(lons_in, fail_silently):
            return lons_in
        lons_out = do_it(lons_in)
        if not is_strictly_increasing(lons_out):
            print("WARNING: You passed in numeric longitudes to lon_idl2pm() and these have been converted, but they're not strictly increasing.")
        print("To assign the new longitude coordinates to an Xarray object, use xarrayobject.assign_coordinates()! (Pass the object directly in to lon_idl2pm() in order to suppress this message.)")
        
    return lons_out


# Convert a longitude axis that's 0 to 360 around the prime meridian to one that's -180 to 180 around the international date line. If you pass in a Dataset or DataArray, the "lon" coordinates will be changed and the axis and data rolled---i.e., maps will be centered on the prime meridian, plus or minus any offset of your gridcell centers. Otherwise, this assumes you're passing in numeric data, and no rolling takes place.
def lon_pm2idl(lons_in, fail_silently = False):
    def check_ok(tmp, fail_silently):
        msg = ""
        if np.any(tmp < 0):
            msg = f"Minimum longitude is already < 0 ({np.min(tmp)})"
        elif np.any(tmp > 360):
            msg = f"Maximum longitude is > 360 ({np.max(tmp)})"
        
        if msg == "":
            return True
        elif fail_silently:
            return False
        else:
            raise ValueError(msg)
        
    def do_it(tmp):
        tmp = np.mod((tmp + 180),360)-180
        return tmp
    
    if isinstance(lons_in, (xr.DataArray, xr.Dataset)):
        if not check_ok(lons_in.lon.values, fail_silently):
            return lons_in
        lons_out = lons_in
        lons_out = lons_out.assign_coords(lon=do_it(lons_in.lon.values))
        lons_out = make_lon_increasing(lons_out)
    else:
        if not check_ok(lons_in, fail_silently):
            return lons_in
        lons_out = do_it(lons_in)
        if not is_strictly_increasing(lons_out):
            print("WARNING: You passed in numeric longitudes to lon_pm2idl() and these have been converted, but they're not strictly increasing.")
        print("To assign the new longitude coordinates to an Xarray object, use xarrayobject.assign_coordinates()! (Pass the object directly in to lon_pm2idl() in order to suppress this message.)")
        
    return lons_out


# List (strings) of managed crops in CLM.
def define_mgdcrop_list():
    notcrop_list = ["tree", "grass", "shrub", "unmanaged", "not_vegetated"]
    defined_pftlist = define_pftlist()
    is_crop = is_each_vegtype(defined_pftlist, notcrop_list, "notok_contains")
    return [defined_pftlist[i] for i, x in enumerate(is_crop) if x]


# Convert list of vegtype strings to integer index equivalents.
def vegtype_str2int(vegtype_str, vegtype_mainlist=None):
    
    convert_to_ndarray = not isinstance(vegtype_str, np.ndarray)
    if convert_to_ndarray:
        vegtype_str = np.array(vegtype_str)
    
    if isinstance(vegtype_mainlist, xr.Dataset):
        vegtype_mainlist = vegtype_mainlist.vegtype_str.values
    elif isinstance(vegtype_mainlist, xr.DataArray):
        vegtype_mainlist = vegtype_mainlist.values
    elif vegtype_mainlist == None:
        vegtype_mainlist = define_pftlist()
    if not isinstance(vegtype_mainlist, list) and isinstance(vegtype_mainlist[0], str):
        if isinstance(vegtype_mainlist, list):
            raise TypeError(f"Not sure how to handle vegtype_mainlist as list of {type(vegtype_mainlist[0])}")
        else:
            raise TypeError(f"Not sure how to handle vegtype_mainlist as type {type(vegtype_mainlist[0])}")

    if vegtype_str.shape == ():
        indices = np.array([-1])
    else:
        indices = np.full(len(vegtype_str), -1)
    for v in np.unique(vegtype_str):
        indices[np.where(vegtype_str == v)] = vegtype_mainlist.index(v)
    if convert_to_ndarray:
        indices = [int(x) for x in indices]
    return indices

# Flexibly subset time(s) and/or vegetation type(s) from an xarray Dataset or DataArray. Keyword arguments like dimension=selection. Selections can be individual values or slice()s. Optimize memory usage by beginning keyword argument list with the selections that will result in the largest reduction of object size. Use dimension "vegtype" to extract patches of designated vegetation type (can be string or integer).
# Can also do dimension=function---e.g., time=np.mean will take the mean over the time dimension.
def xr_flexsel(xr_object, patches1d_itype_veg=None, warn_about_seltype_interp=True, **kwargs):
    
    # Setup
    havewarned = False
    delimiter = "__"
    
    for key, selection in kwargs.items():
        
        if callable(selection):
            # It would have been really nice to do selection(xr_object, axis=key), but numpy methods and xarray methods disagree on "axis" vs. "dimension." So instead, just do this manually.
            if selection == np.mean:
                try:
                    xr_object = xr_object.mean(dim=key)
                except:
                    raise ValueError(f"Failed to take mean of dimension {key}. Try doing so outside of xr_flexsel().")
            else:
                raise ValueError(f"xr_flexsel() doesn't recognize function {selection}")

        elif key == "vegtype":
            
            # Convert to list, if needed
            if not isinstance(selection, list):
                selection = [selection]
            
            # Convert to indices, if needed
            if isinstance(selection[0], str):
                selection = vegtype_str2int(selection)
            
            # Get list of boolean(s)
            if isinstance(selection[0], int):
                if isinstance(patches1d_itype_veg, type(None)):
                    patches1d_itype_veg = xr_object.patches1d_itype_veg.values
                elif isinstance(patches1d_itype_veg, xr.core.dataarray.DataArray):
                    patches1d_itype_veg = patches1d_itype_veg.values
                is_vegtype = is_each_vegtype(patches1d_itype_veg, selection, "ok_exact")
            elif isinstance(selection[0], bool):
                if len(selection) != len(xr_object.patch):
                    raise ValueError(f"If providing boolean 'vegtype' argument to xr_flexsel(), it must be the same length as xr_object.patch ({len(selection)} vs. {len(xr_object.patch)})")
                is_vegtype = selection
            else:
                raise TypeError(f"Not sure how to handle 'vegtype' of type {type(selection[0])}")
            xr_object = xr_object.isel(patch=[i for i, x in enumerate(is_vegtype) if x])
            if "ivt" in xr_object:
                xr_object = xr_object.isel(ivt=is_each_vegtype(xr_object.ivt.values, selection, "ok_exact"))
        
        else:
            
            # Parse selection type, if provided
            if delimiter in key:
                key, selection_type = key.split(delimiter)
                
            # Check type of selection
            else:
                
                is_inefficient = False
                if isinstance(selection, slice):
                    slice_members = []
                    if selection == slice(0):
                        raise ValueError("slice(0) will be empty")
                    if selection.start != None:
                        slice_members = slice_members + [selection.start]
                    if selection.stop != None:
                        slice_members = slice_members + [selection.stop]
                    if selection.step != None:
                        slice_members = slice_members + [selection.step]
                    if slice_members==[]:
                        raise TypeError("slice is all None?")
                    this_type = int
                    for x in slice_members:
                        if x < 0 or not isinstance(x, int):
                            this_type = "values"
                            break
                elif isinstance(selection, np.ndarray):
                    if selection.dtype.kind in np.typecodes["AllInteger"]:
                        this_type = int
                    else:
                        is_inefficient = True
                        this_type = None
                        for x in selection:
                            if x < 0 or x%1 > 0:
                                if isinstance(x, int):
                                    this_type = "values"
                                else:
                                    this_type = type(x)
                                break
                        if this_type==None:
                            this_type = int
                            selection = selection.astype(int)
                else:
                    this_type = type(selection)
                
                warn_about_this_seltype_interp = warn_about_seltype_interp
                if this_type == list and isinstance(selection[0], str):
                    selection_type = "values"
                    warn_about_this_seltype_interp = False
                elif this_type == int:
                    selection_type = "indices"
                else:
                    selection_type = "values"
                
                if warn_about_this_seltype_interp:
                    # Suggest suppressing selection type interpretation warnings
                    if not havewarned:
                        print("xr_flexsel(): Suppress all 'selection type interpretation' messages by specifying warn_about_seltype_interp=False")
                        havewarned = True
                    if is_inefficient:
                        extra =  " This will also improve efficiency for large selections."
                    else:
                        extra = ""
                    print(f"xr_flexsel(): Selecting {key} as {selection_type} because selection was interpreted as {this_type}. If not correct, specify selection type ('indices' or 'values') in keyword like '{key}{delimiter}SELECTIONTYPE=...' instead of '{key}=...'.{extra}")
                    
            # Trim along relevant 1d axes
            if isinstance(xr_object, xr.Dataset) and key in ["lat","lon"]:
                if selection_type == "indices":
                    inclCoords = xr_object[key].values[selection]
                elif selection_type == "values":
                    if isinstance(selection, slice):
                        inclCoords = xr_object.sel({key: selection}, drop=False)[key].values
                    else:
                        inclCoords = selection
                else:
                    raise TypeError(f"selection_type {selection_type} not recognized")                    
                if key == "lat":
                    thisXY = "jxy"
                elif key=="lon":
                    thisXY = "ixy"
                else:
                    raise KeyError(f"Key '{key}' not recognized: What 1d_ suffix should I use for variable name?")
                pattern = re.compile(f"1d_{thisXY}")
                matches = [x for x in list(xr_object.keys()) if pattern.search(x) != None]
                for thisVar in matches:
                    if len(xr_object[thisVar].dims) != 1:
                        raise RuntimeError(f"Expected {thisVar} to have 1 dimension, but it has {len(xr_object[thisVar].dims)}: {xr_object[thisVar].dims}")
                    thisVar_dim = xr_object[thisVar].dims[0]
                    # print(f"Variable {thisVar} has dimension {thisVar_dim}")
                    thisVar_coords = xr_object[key].values[xr_object[thisVar].values.astype(int)-1]
                    # print(f"{thisVar_dim} size before: {xr_object.sizes[thisVar_dim]}")
                    ok_ind = []
                    new_1d_thisXY = []
                    for i, x in enumerate(thisVar_coords):
                        if x in inclCoords:
                           ok_ind = ok_ind + [i]
                           new_1d_thisXY = new_1d_thisXY + [(inclCoords==x).nonzero()[0] + 1]
                    xr_object = xr_object.isel({thisVar_dim: ok_ind})
                    new_1d_thisXY = np.array(new_1d_thisXY).squeeze()
                    xr_object[thisVar].values = new_1d_thisXY
                    # print(f"{thisVar_dim} size after: {xr_object.sizes[thisVar_dim]}")

            
            # Perform selection
            if selection_type == "indices":
                # Have to select like this instead of with index directly because otherwise assign_coords() will throw an error. Not sure why.
                if isinstance(selection, int):
                    # Single integer? Turn it into a slice.
                    selection = slice(selection,selection+1)
                elif isinstance(selection, np.ndarray) and not selection.dtype.kind in np.typecodes["AllInteger"]:
                    selection = selection.astype(int)
                xr_object = xr_object.isel({key: selection})
            elif selection_type == "values":
                xr_object = xr_object.sel({key: selection})
            else:
                raise TypeError(f"selection_type {selection_type} not recognized")
    
    return xr_object


# Get PFT of each patch, in both integer and string forms.
def get_patch_ivts(this_ds, this_pftlist):
    # First, get all the integer values; should be time*pft or pft*time. We will eventually just take the first timestep.
    vegtype_int = this_ds.patches1d_itype_veg
    vegtype_int.values = vegtype_int.values.astype(int)

    # Convert to strings.
    vegtype_str = list(np.array(this_pftlist)[vegtype_int.values])

    # Return a dictionary with both results
    return {"int": vegtype_int, "str": vegtype_str, "all_str": this_pftlist}


# Convert a list of strings with vegetation type names into a DataArray. Used to add vegetation type info in import_ds().
def get_vegtype_str_da(vegtype_str):
    nvt = len(vegtype_str)
    thisName = "vegtype_str"
    vegtype_str_da = xr.DataArray(\
        vegtype_str, 
        coords={"ivt": np.arange(0,nvt)}, 
        dims=["ivt"],
        name = thisName)
    return vegtype_str_da


# Function to drop unwanted variables in preprocessing of open_mfdataset(), making sure to NOT drop any unspecified variables that will be useful in gridding. Also adds vegetation type info in the form of a DataArray of strings.
# Also renames "pft" dimension (and all like-named variables, e.g., pft1d_itype_veg_str) to be named like "patch". This can later be reversed, for compatibility with other code, using patch2pft().
def mfdataset_preproc(ds, vars_to_import, vegtypes_to_import, timeSlice):
    
    # Rename "pft" dimension and variables to "patch", if needed
    if "pft" in ds.dims:
        pattern = re.compile("pft.*1d")
        matches = [x for x in list(ds.keys()) if pattern.search(x) != None]
        pft2patch_dict = {"pft": "patch"}
        for m in matches:
            pft2patch_dict[m] = m.replace("pft","patch").replace("patchs","patches")
        ds = ds.rename(pft2patch_dict)

    derived_vars = []
    if vars_to_import != None:
        # Split vars_to_import into variables that are vs. aren't already in ds
        derived_vars = [v for v in vars_to_import if v not in ds]
        present_vars = [v for v in vars_to_import if v in ds]
        vars_to_import = present_vars
        
        # Get list of dimensions present in variables in vars_to_import.
        dimList = []
        for thisVar in vars_to_import:
            # list(set(x)) returns a list of the unique items in x
            dimList = list(set(dimList + list(ds.variables[thisVar].dims)))
        
        # Get any _1d variables that are associated with those dimensions. These will be useful in gridding. Also, if any dimension is "pft", set up to rename it and all like-named variables to "patch"
        onedVars = []
        for thisDim in dimList:
            pattern = re.compile(f"{thisDim}.*1d")
            matches = [x for x in list(ds.keys()) if pattern.search(x) != None]
            onedVars = list(set(onedVars + matches))
        
        # Add dimensions and _1d variables to vars_to_import
        vars_to_import = list(set(vars_to_import \
            + list(ds.dims) + onedVars))
        
        # Add any _bounds variables
        bounds_vars = []
        for v in vars_to_import:
            bounds_var = v + "_bounds"
            if bounds_var in ds:
                bounds_vars = bounds_vars + [bounds_var]
        vars_to_import = vars_to_import + bounds_vars
        
        # Get list of variables to drop
        varlist = list(ds.variables)
        vars_to_drop = list(np.setdiff1d(varlist, vars_to_import))

        # Drop them
        ds = ds.drop_vars(vars_to_drop)
    
    # Add vegetation type info
    if "patches1d_itype_veg" in list(ds):
        this_pftlist = define_pftlist()
        get_patch_ivts(ds, this_pftlist) # Includes check of whether vegtype changes over time anywhere
        vegtype_da = get_vegtype_str_da(this_pftlist)
        patches1d_itype_veg_str = vegtype_da.values[ds.isel(time=0).patches1d_itype_veg.values.astype(int)]
        npatch = len(patches1d_itype_veg_str)
        patches1d_itype_veg_str = xr.DataArray( \
            patches1d_itype_veg_str,
            coords={"patch": np.arange(0,npatch)}, 
            dims=["patch"],
            name = "patches1d_itype_veg_str")
        ds = xr.merge([ds, vegtype_da, patches1d_itype_veg_str])

    # Restrict to veg. types of interest, if any
    if vegtypes_to_import != None:
        ds = xr_flexsel(ds, vegtype=vegtypes_to_import)
    
    # Restrict to time slice, if any
    if timeSlice:
        ds = ds.sel(time=timeSlice)

    # Finish import
    ds = xr.decode_cf(ds, decode_times = True)
    
    # Compute derived variables
    for v in derived_vars:
        if v == "HYEARS" and "HDATES" in ds and ds.HDATES.dims == ('time', 'mxharvests', 'patch'):
            yearList = np.array([np.float32(x.year - 1) for x in ds.time.values])
            hyears = ds["HDATES"].copy()
            hyears.values = np.tile(np.expand_dims(yearList, (1,2)), (1, ds.dims["mxharvests"], ds.dims["patch"]))
            hyears.values[ds.HDATES.values<=0] = ds.HDATES.values[ds.HDATES.values<=0]
            hyears.values[np.isnan(ds.HDATES.values)] = np.nan
            hyears.attrs["long_name"] = "DERIVED: actual crop harvest years"
            hyears.attrs["units"] = "year"
            ds["HYEARS"] = hyears
    
    return ds


# Rename "patch" dimension and any associated variables back to "pft". Uses a dictionary with the names of the dimensions and variables we want to rename. This allows us to do it all at once, which may be more efficient than one-by-one.
def patch2pft(xr_object):

    # Rename "patch" dimension
    patch2pft_dict = {}
    for thisDim in xr_object.dims:
        if thisDim == "patch":
            patch2pft_dict["patch"] = "pft"
            break
    
    # Rename variables containing "patch"
    if isinstance(xr_object, xr.Dataset):
        pattern = re.compile("patch.*1d")
        matches = [x for x in list(xr_object.keys()) if pattern.search(x) != None]
        if len(matches) > 0:
            for m in matches:
                patch2pft_dict[m] = m.replace("patches","patchs").replace("patch","pft")
    
    # Do the rename
    if len(patch2pft_dict) > 0:
        xr_object = xr_object.rename(patch2pft_dict)
    
    return xr_object


# Import a dataset that can be spread over multiple files, only including specified variables and/or vegetation types and/or timesteps, concatenating by time. DOES actually read the dataset into memory, but only AFTER dropping unwanted variables and/or vegetation types.
def import_ds(filelist, myVars=None, myVegtypes=None, timeSlice=None, myVars_missing_ok=[], only_active_patches=False, rename_lsmlatlon=False, chunks=None):
    
    # Convert myVegtypes here, if needed, to avoid repeating the process each time you read a file in xr.open_mfdataset().
    if myVegtypes != None:
        if not isinstance(myVegtypes, list):
            myVegtypes = [myVegtypes]
        if isinstance(myVegtypes[0], str):
            myVegtypes = vegtype_str2int(myVegtypes)
    
    # Same for these variables.
    if myVars != None:
        if not isinstance(myVars, list):
            myVars = [myVars]
    if myVars_missing_ok:
        if not isinstance(myVars_missing_ok, list):
            myVars_missing_ok = [myVars_missing_ok]
            
    # Make sure lists are actually lists
    if not isinstance(filelist, list):
        filelist = [filelist]
    if not isinstance(myVars_missing_ok, list):
        myVars_missing_ok = [myVars_missing_ok]
            
    # Remove files from list if they don't contain requested timesteps.
    # timeSlice should be in the format slice(start,end[,step]). start or end can be None to be unbounded on one side. Note that the standard slice() documentation suggests that only elements through end-1 will be selected, but that seems not to be the case in the xarray implementation.
    if timeSlice:
        new_filelist = []
        for file in sorted(filelist):
            if xr.open_dataset(file).time.sel(time=timeSlice).size:
                new_filelist.append(file)
            
             # If you found some matching files, but then you find one that doesn't, stop going through the list.
            elif new_filelist:
                break
        if not new_filelist:
            raise RuntimeError(f"No files found in timeSlice {timeSlice}")
        filelist = new_filelist

    # The xarray open_mfdataset() "preprocess" argument requires a function that takes exactly one variable (an xarray.Dataset object). Wrapping mfdataset_preproc() in this lambda function allows this. Could also just allow mfdataset_preproc() to access myVars and myVegtypes directly, but that's bad practice as it could lead to scoping issues.
    mfdataset_preproc_closure = \
        lambda ds: mfdataset_preproc(ds, myVars, myVegtypes, timeSlice)

    # Import
    if isinstance(filelist, list):
        this_ds = xr.open_mfdataset(sorted(filelist), \
            data_vars="minimal", 
            preprocess=mfdataset_preproc_closure,
            compat='override',
            coords='all',
            chunks=chunks)
    elif isinstance(filelist, str):
        this_ds = xr.open_dataset(filelist, chunks=chunks)
        this_ds = mfdataset_preproc(this_ds, myVars, myVegtypes, timeSlice)
        this_ds = this_ds.compute()
        
    # Include only active patches (or whatever)
    if only_active_patches:
        is_active = this_ds.patches1d_active.values
        p_active = np.where(is_active)[0]
        this_ds_active = this_ds.isel(patch=p_active)
    
    # Warn and/or error about variables that couldn't be imported or derived
    if myVars:
        missing_vars = [v for v in myVars if v not in this_ds]
        ok_missing_vars = [v for v in missing_vars if v in myVars_missing_ok]
        bad_missing_vars = [v for v in missing_vars if v not in myVars_missing_ok]
        if ok_missing_vars:
            print(f"Could not import some variables; either not present or not deriveable: {ok_missing_vars}")
        if bad_missing_vars:
            raise RuntimeError(f"Could not import some variables; either not present or not deriveable: {bad_missing_vars}")
    
    if rename_lsmlatlon:
        if "lsmlat" in this_ds.dims:
            this_ds = this_ds.rename({'lsmlat': 'lat'})
        if "lsmlon" in this_ds.dims:
            this_ds = this_ds.rename({'lsmlon': 'lon'})
    
    return this_ds


# Return a DataArray, with defined coordinates, for a given variable in a dataset.
def get_thisVar_da(thisVar, this_ds):

    # Make DataArray for this variable
    thisvar_da = np.array(this_ds.variables[thisVar])
    theseDims = this_ds.variables[thisVar].dims
    thisvar_da = xr.DataArray(thisvar_da, 
        dims = theseDims)

    # Define coordinates of this variable's DataArray
    dimsDict = dict()
    for thisDim in theseDims:
        dimsDict[thisDim] = this_ds[thisDim]
    thisvar_da = thisvar_da.assign_coords(dimsDict)
    thisvar_da.attrs = this_ds[thisVar].attrs

    return thisvar_da


# Given a DataArray, remove all patches except those planted with managed crops.
def trim_da_to_mgd_crop(thisvar_da, patches1d_itype_veg_str):

    # Handle input DataArray without patch dimension
    if not any(np.array(list(thisvar_da.dims)) == "patch"):
        print("Input DataArray has no patch dimension and therefore trim_to_mgd_crop() has no effect.")
        return thisvar_da
    
    # Throw error if patches1d_itype_veg_str isn't strings
    if isinstance(patches1d_itype_veg_str, xr.DataArray):
        patches1d_itype_veg_str = patches1d_itype_veg_str.values
    if not isinstance(patches1d_itype_veg_str[0], str):
        raise TypeError("Input patches1d_itype_veg_str is not in string form, and therefore trim_to_mgd_crop() cannot work.")
    
    # Get boolean list of whether each patch is planted with a managed crop
    notcrop_list = ["tree", "grass", "shrub", "unmanaged", "not_vegetated"]
    is_crop = is_each_vegtype(patches1d_itype_veg_str, notcrop_list, "notok_contains")

    # Warn if no managed crops were found, but still return the empty result
    if np.all(np.bitwise_not(is_crop)):
        print("No managed crops found! Returning empty DataArray.")
    return thisvar_da.isel(patch = [i for i, x in enumerate(is_crop) if x])


# Make a geographically gridded DataArray (with dimensions time, vegetation type [as string], lat, lon) of one variable within a Dataset. Optional keyword arguments will be passed to xr_flexsel() to select single steps or slices along the specified ax(ie)s.
#
# fillValue: Default None means grid will be filled with NaN, unless the variable in question already has a fillValue, in which case that will be used.
def grid_one_variable(this_ds, thisVar, fillValue=None, **kwargs):
    
    # Get this Dataset's values for selection(s), if provided
    this_ds = xr_flexsel(this_ds, \
        **kwargs)
    
    # Get DataArrays needed for gridding
    thisvar_da = get_thisVar_da(thisVar, this_ds)
    vt_da = None
    if "patch" in thisvar_da.dims:
        spatial_unit = "patch"
        xy_1d_prefix = "patches"
        if "patches1d_itype_veg" in this_ds:
            vt_da = get_thisVar_da("patches1d_itype_veg", this_ds)
    elif "gridcell" in thisvar_da.dims:
        spatial_unit = "gridcell"
        xy_1d_prefix = "grid"
    else:
        raise RuntimeError(f"What variables to use for _ixy and _jxy of variable with dims {thisvar_da.dims}?")
    ixy_da = get_thisVar_da(xy_1d_prefix + "1d_ixy", this_ds)
    jxy_da = get_thisVar_da(xy_1d_prefix + "1d_jxy", this_ds)
    
    
    if not fillValue and "_FillValue" in thisvar_da.attrs:
        fillValue = thisvar_da.attrs["_FillValue"]
    
    # Renumber vt_da to work as indices on new ivt dimension, if needed.
    ### Ensures that the unique set of vt_da values begins with 1 and
    ### contains no missing steps.
    if "ivt" in this_ds and vt_da is not None:
        vt_da.values = np.array([np.where(this_ds.ivt.values == x)[0][0] for x in vt_da.values])
    
    # Get new dimension list
    new_dims = list(thisvar_da.dims)
    ### Remove "[spatial_unit]".
    if spatial_unit in new_dims:
        new_dims.remove(spatial_unit)
    #  Add "ivt_str" (vegetation type, as string). This needs to go at the end, to avoid a possible situation where you wind up with multiple Ellipsis members of fill_indices.
    if "ivt" in this_ds and spatial_unit=="patch":
        new_dims.append("ivt_str")
    ### Add lat and lon to end of list
    new_dims = new_dims + ["lat", "lon"]

    # Set up empty array
    n_list = []
    for dim in new_dims:
        if dim == "ivt_str":
            n = this_ds.sizes["ivt"]
        elif dim in thisvar_da.coords:
            n = thisvar_da.sizes[dim]
        else:
            n = this_ds.sizes[dim]
        n_list = n_list + [n]
    thisvar_gridded = np.empty(n_list)
    if fillValue:
        thisvar_gridded[:] = fillValue
    else:
        thisvar_gridded[:] = np.NaN

    # Fill with this variable
    fill_indices = []
    for dim in new_dims:
        if dim == "lat":
            fill_indices.append(jxy_da.values.astype(int) - 1)
        elif dim == "lon":
            fill_indices.append(ixy_da.values.astype(int) - 1)
        elif dim == "ivt_str":
            fill_indices.append(vt_da)
        elif not fill_indices:
        # I.e., if fill_indices is empty. Could also do "elif len(fill_indices)==0".
            fill_indices.append(Ellipsis)
    try:
        thisvar_gridded[tuple(fill_indices[:len(fill_indices)])] = thisvar_da.values
    except:
        thisvar_gridded[tuple(fill_indices[:len(fill_indices)])] = thisvar_da.values.transpose()
    if not np.any(np.bitwise_not(np.isnan(thisvar_gridded))):
        if np.all(np.isnan(thisvar_da.values)):
            print('Warning: This DataArray (and thus map) is all NaN')
        else:
            raise RuntimeError("thisvar_gridded was not filled!")
    
    # Assign coordinates, attributes and name
    thisvar_gridded = xr.DataArray(thisvar_gridded, \
        dims=tuple(new_dims),
        attrs=thisvar_da.attrs)
    for dim in new_dims:
        if dim == "ivt_str":
            values = this_ds.vegtype_str.values
        elif dim in thisvar_da.coords:
            values = thisvar_da[dim]
        else:
            values = this_ds[dim].values
        thisvar_gridded = thisvar_gridded.assign_coords({dim: values})
    thisvar_gridded.name = thisVar
    
    # Add FillValue attribute
    if fillValue:
        thisvar_gridded.attrs["_FillValue"] = fillValue

    return thisvar_gridded


# Xarray's native resampler is nice, but it will result in ALL variables being resampled
# along the N specified dimension(s). This means that, e.g., variables that were
# supposed to be 1d will be 1+Nd afterwards. This function undoes that. The syntax is the
# same as for Xarray's resampler plus, at the beginning: 
#    (1) which object you want to resample, and
#    (2) the function you want to use for downsampling (e.g., "mean").
def resample(ds_in, thefunction, **kwargs):
    
    # This problem is not applicable to DataArrays, so just use Xarray's resampler.
    if isinstance(ds_in, xr.DataArray):
        da_resampler = ds_in.resample(**kwargs)
        da_out = getattr(da_resampler, thefunction)()
        return da_out
    
    # Get the original dimensions of each variable
    orig_dims = dict([(x, ds_in[x].dims) for x in ds_in])
    
    # Do the Xarray resampling
    ds_resampler = ds_in.resample(**kwargs)
    ds_out = getattr(ds_resampler, thefunction)()
    
    for v in ds_out:
        
        extra_dims = [d for d in ds_out[v].dims if d not in orig_dims[v]]
        if not extra_dims:
            print(f"Skipping {v}")
            continue
        
        # Xarray's resampler for now only supports resampling along one
        # dimension. However, I'm going to try and support a future
        # version that supports an arbitrary number of resampled dimensions;
        # this is, of course, untested.
        
        # For any newly-created dimensions, select the first value. Should
        # Be the same as all other values.
        indexers = dict([(d,0) for d in extra_dims])

        ds_out[v] = ds_out[v].isel(**indexers, drop=True)
    
    return ds_out


# Have a DataArray that's just for one time point? Repeat it for many.
def tile_over_time(da_in, years=None):
    
    if "time" not in da_in.dims:
        raise RuntimeError('Rework tile_over_time() to function with da_in lacking time dimension.')
    
    # Deal with Datasets
    if isinstance(da_in, xr.Dataset):
        new_attrs = {}
        for x in da_in.attrs:
            if x == "created":
                continue
            else:
                new_attrs[x] = da_in.attrs[x]
        ds_out = xr.Dataset(attrs=new_attrs)
        for v in da_in:
            if 'time' in da_in[v].dims and v != "time_bounds":
                ds_out[v] = tile_over_time(da_in[v], years=years)
            else:
                ds_out[v] = da_in[v].copy()
        return ds_out
    elif not isinstance(da_in, xr.DataArray):
        raise RuntimeError(f'tile_over_time() only works with xarray Datasets and DataArrays, not {type(da_in)}')
    
    # Get info about time in input DataArray
    dt0 = da_in.time.values[0]
    dt_type = type(dt0)
    has_year_zero = dt0.has_year_zero
    
    if type(years) != type(None):
        new_time = np.array([dt_type(x, dt0.month, dt0.day, has_year_zero=has_year_zero) for x in years])
        strftime_fmt = '%Y-%m-%d'
        Ntime = len(years)
    else:
        raise RuntimeError('Rework tile_over_time to work with something other than years')
    
    # Convert from cftime to int days since, for compatibility with NETCDF3-CLASSIC format.
    new_time_units = 'days since ' + new_time[0].strftime(strftime_fmt)
    new_time = cftime.date2num(new_time, new_time_units)
    
    # Set up time DataArray
    new_time_attrs = da_in.time.attrs
    new_time_attrs['units'] = new_time_units
    new_time_attrs['calendar'] = dt0.calendar
    new_time_da = xr.DataArray(new_time, dims=['time'], coords={'time': new_time}, attrs=new_time_attrs)
    
    # Get coordinates to be used in new DataArray
    new_coords = {}
    for x in da_in.coords:
        if x == "time":
            new_coords[x] = new_time_da
        else:
            new_coords[x] = da_in.coords[x]
    
    # Set up broadcasting DataArray
    new_shape = tuple([da_in.shape[i] if x != "time" else Ntime for i,x in enumerate(da_in.dims)])
    bc_da = xr.DataArray(np.ones(new_shape), dims=da_in.dims, coords=new_coords)
    
    # Get new DataArray
    da_out = (da_in.squeeze() * bc_da).transpose(*da_in.dims)
    da_out = da_out.assign_attrs(da_in.attrs)
    
    return da_out