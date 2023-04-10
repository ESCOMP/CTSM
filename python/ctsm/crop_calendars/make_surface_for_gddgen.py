import numpy as np
import xarray as xr
import argparse
import sys
import os

# Doing what I did before with flanduse_timeseries, except now with
# fsurdat.
def get_new_fsurdat_v0(surf, lu):
    
    # Where is each crop ever active?
    if "AREA" not in lu:
        lu["AREA"] = xr.DataArray(data=np.ones_like(lu["LANDFRAC_PFT"]),
                                  dims=lu["LANDFRAC_PFT"].dims,
                                  coords=lu["LANDFRAC_PFT"].coords)
    lu['AREA_CROP'] = (lu.AREA * lu.LANDFRAC_PFT * lu.PCT_CROP/100).transpose("time", "lsmlat", "lsmlon")
    lu['AREA_CFT'] = (lu['AREA_CROP'] * lu.PCT_CFT/100).transpose("time", "cft", "lsmlat", "lsmlon")
    ds_area_max = lu['AREA_CFT'].max(dim="time")
    lu['ever_active_bycft'] = (ds_area_max > 0).transpose("cft", "lsmlat", "lsmlon")
    lu['ever_active'] = lu['ever_active_bycft'].any(dim="cft")
    
    # For cells that EVER have cropland, when do they have 0 crop area? Change those 0% values to 1%.
    new_pct_crop_ar = lu['PCT_CROP'].values
    needs_pctcrop_0to1 = np.where((lu['AREA_CROP']==0) & lu['ever_active'])
    new_pct_crop_ar[needs_pctcrop_0to1] = 1.0
    lu['PCT_CROP'] = xr.DataArray(data = new_pct_crop_ar,
                                  coords = lu['PCT_CROP'].coords,
                                  attrs = lu['PCT_CROP'].attrs)
    lu['AREA_CROP'] = (lu.AREA * lu.LANDFRAC_PFT * lu.PCT_CROP/100).transpose("time", "lsmlat", "lsmlon")
    if np.any((lu['AREA_CROP']==0) & lu['ever_active']):
        raise RuntimeError("Failed to fill 0% CROP with 1% where needed.")

    # For cells that EVER have each CFT, when do they have 0 area of that CFT? Change those 0% values to something positive.
    new_pct_cft_tcyx = lu['PCT_CFT'].values
    new_pct_cft_tcyx[np.where((lu['AREA_CFT']==0) & lu['ever_active'])] = 1.0
    # Ensure sum to 100
    i = 0
    while np.any(~np.isclose(np.sum(new_pct_cft_tcyx, axis=1), 100.0)):
        i+=1
        if i > 10:
            raise RuntimeError('too many iterations')
        new_pct_cft_tcyx = 100 * (new_pct_cft_tcyx / np.expand_dims(np.sum(new_pct_cft_tcyx, axis=1), axis=1))
    lu['PCT_CFT'] = xr.DataArray(data = new_pct_cft_tcyx,
                                 coords = lu['PCT_CFT'].coords,
                                 attrs = lu['PCT_CFT'].attrs)
    
    # Just take the first timestep
    new_pct_crop_da = lu['PCT_CROP'].isel(time=0)
    new_pct_cft_da = lu['PCT_CFT'].isel(time=0)
    
    # flanduse_timeseries doesn't have PCT_NATVEG
    new_natveg_da = surf['PCT_NATVEG']
    
    # For cells where 
    
    return new_pct_crop_da, new_natveg_da, new_pct_cft_da


# Trying to minimize the number of crop PFTs that need to be simulated.
# CLM currently giving errors about wt_cft not summing to 1.
def get_new_fsurdat_v1(surf, lu, params):
    
    if params is None:
        raise RuntimeError("You must provide -p/--paramfile")
    
    ##########################################
    ### %% Get new PCT_CROP and PCT_NATVEG ###
    ##########################################

    new_pct_crop_da = lu['PCT_CROP'].max(dim="time")
    new_pct_crop_da.attrs = surf['PCT_CROP'].attrs

    surf_pct_crop_plus_natveg = surf['PCT_CROP'] + surf['PCT_NATVEG']
    if np.any(surf_pct_crop_plus_natveg < new_pct_crop_da):
        raise RuntimeError("Max CROP > CROP+NATVEG")

    new_natveg_da = surf_pct_crop_plus_natveg - new_pct_crop_da
    if np.any((new_natveg_da > 0) & (surf['PCT_NATVEG'] == 0)):
        print("You created some NATVEG area. Not necessarily a problem, but unexpected.")
    new_natveg_da.attrs = surf['PCT_NATVEG'].attrs
    
    
    ##################################################################
    ### Get new PCT_CFT (percentage of cropland that is each crop) ###
    ##################################################################
    
    # Sum all crops' max area, merging unrepresented types into their representative type
    cft_list_int = surf['cft'].values
    max_merged_pct_crop = np.full_like(surf['PCT_CFT'], 0.0)
    for i, c in enumerate(cft_list_int):
        mergetarget = params['mergetoclmpft'].sel(pft=c).values
        m = np.where(cft_list_int==mergetarget)[0]
        max_merged_pct_crop[m,:,:] += np.expand_dims(lu['PCT_CFT'].sel(cft=c).max(dim="time"), axis=0)
    max_merged_pct_crop_da = xr.DataArray(data=max_merged_pct_crop,
                                          dims=surf['PCT_CFT'].dims,
                                          attrs=surf['PCT_CFT'].attrs)
    
    # Ensure no area in merged-away crops
    for i, c in enumerate(cft_list_int):
        if (params['mergetoclmpft'].sel(pft=c) != c) and (max_merged_pct_crop_da.sel(cft=i).max() > 0):
            raise RuntimeError(f"Unexpected max_merged_pct_crop area for pft {c}")

    # Determine how many crops ever have any area
    ever_has_this_crop = np.full_like(max_merged_pct_crop, 0.0)
    ever_has_this_crop[max_merged_pct_crop > 0] = 1
    N_crops_ever_active = ever_has_this_crop.sum(axis=0)

    # Split crop area evenly among ever-included crops
    new_pct_cft = np.full_like(surf['PCT_CFT'].isel(cft=0), 0.0)
    new_pct_cft[N_crops_ever_active > 0] = 100 / N_crops_ever_active[N_crops_ever_active > 0]
    new_pct_cft = np.expand_dims(new_pct_cft, axis=0)
    new_pct_cft = np.tile(new_pct_cft, reps=[surf.dims['cft'], 1, 1])
    where_zero = np.where(max_merged_pct_crop_da)
    new_pct_cft[where_zero] = 0.0
    new_pct_cft_da = xr.DataArray(data=new_pct_cft,
                                  dims=surf['PCT_CFT'].dims,
                                  attrs=surf['PCT_CFT'].attrs)
    
    return new_pct_crop_da, new_natveg_da, new_pct_cft_da


def main(argv):
    
    ###############################
    ### Process input arguments ###
    ###############################
    
    # Set arguments
    parser = argparse.ArgumentParser(description="ADD DESCRIPTION HERE")
    parser.add_argument("-l", "--flanduse_timeseries", "--flanduse-timeseries", 
                        help="Land-use timeseries file (flanduse_timeseries) for CLM run",
                        required=True)
    parser.add_argument("-p", "--paramfile", 
                        help="Parameter file (paramfile) for CLM run")
    parser.add_argument("-s", "--fsurdat", 
                        help="Surface dataset (fsurdat) for CLM run",
                        required=True)
    parser.add_argument("-o", "--outfile", 
                        help="Output fsurdat file")
    args = parser.parse_args(argv)


    ##########################
    ### Import and process ###
    ##########################

    lu = xr.open_dataset(args.flanduse_timeseries)
    surf = xr.open_dataset(args.fsurdat)
    if args.paramfile is not None:
        params = xr.open_dataset(args.paramfile)

    # new_pct_crop_da, new_natveg_da, new_pct_cft_da = get_new_fsurdat_v1(surf, lu, params)
    new_pct_crop_da, new_natveg_da, new_pct_cft_da = get_new_fsurdat_v0(surf, lu)
    
    
    #############################
    ### Save to run directory ###
    #############################
    
    # Make new Dataset
    new_surf = surf
    new_surf['PCT_CROP'] = new_pct_crop_da
    new_surf['PCT_CFT'] = new_pct_cft_da
    new_surf['PCT_NATVEG'] = new_natveg_da
    
    # Save to new file
    if outfile is None:
        fsurdat_noext, ext = os.path.splitext(args.fsurdat)
        outfile = f"{fsurdat_noext}.GDDgen{ext}"
        outfile = os.path.basename(outfile)
    new_surf.to_netcdf(outfile, format="NETCDF3_64BIT")

    print(outfile)


if __name__ == "__main__":
    main(sys.argv[1:])
