"""
Replace values in default parameter set with values from Jackie's South America work
"""
import os
import numpy as np
import xarray as xr

PARAM_FILES_DIR = "/glade/u/home/samrabin/ctsm_fates-isforest/src/fates/parameter_files"
file_in = os.path.join(PARAM_FILES_DIR, "fates_params_shuman24.a8f9da07.nc")
file_out = os.path.join(PARAM_FILES_DIR, "fates_params_shuman24.a8f9da07.withchanges.nc")

ds_in = xr.open_dataset(file_in)
ds_out = ds_in.copy()

def replace_values(ds, var, vals_in, indices=None):
    """
    Replace values in DataArray
    """
    da_in = ds[var]
    print(var + ": " + da_in.attrs["long_name"])
    ndims = len(da_in.dims)
    if indices is None:
        vals_in = np.array(vals_in)
        if vals_in.size != da_in.size:
            raise RuntimeError(f"vals_in.size {vals_in.size} but da_in.size {da_in.size}")
        vals_new = vals_in
    else:
        vals_new = da_in.copy().values

        # Handle shortcut to replace all values
        if indices == "all":
            if isinstance(vals_in, (list, np.ndarray)):
                raise RuntimeError("For indices 'all', you must provide a single scalar vals_in")
            indices = np.arange(da_in.size)
            vals_in = np.full(shape=indices.shape, fill_value=vals_in)

        # Check input lists
        if not isinstance(vals_in, (list, np.ndarray)):
            vals_in = [vals_in]
        if not isinstance(indices, (list, np.ndarray)):
            indices = [indices]
        if max(indices) + 1 > da_in.size:
            raise RuntimeError(f"max(indices) {max(indices)} but da_in.size {da_in.size}")
        if len(vals_in) != len(indices):
            raise RuntimeError(f"len(vals_in) {len(vals_in)} but len(indices) {len(indices)}")

        # Convert to 1-d, if needed
        did_reshape = ndims > 1
        if did_reshape:
            orig_dims = vals_new.shape
            vals_new = vals_new.ravel()

        # Assign new values
        for i, val in enumerate(vals_in):
            vals_new[indices[i]] = val

        # Convert back from 1-d, if needed
        if did_reshape:
            vals_new = vals_new.reshape(orig_dims)

    # Print results
    vals_new_raveled = vals_new.ravel()
    this_str = ""
    any_diff = False
    for i, val_orig in enumerate(da_in.values.ravel()):
        val_new = vals_new_raveled[i]
        if val_orig != val_new:
            this_str += f"{val_orig} → {val_new}"
            any_diff = True
        elif indices is None or i in indices:
            this_str += "✔️"
        else:
            this_str += "_"
        this_str += ", "
    this_str = this_str[:-2]
    if not any_diff:
        this_str += " (Jackie's change[s] already incorporated into default param file)"
    else:
        print(da_in.values)
        print(vals_new)
    print(this_str)

    # Save to Dataset
    new_da = xr.DataArray(
        data=vals_new,
        attrs=da_in.attrs,
        dims=da_in.dims,
        coords=da_in.coords,
    )
    ds[var] = new_da
    print(" ")
    return ds

ds_out = replace_values(ds_out, "fates_fire_FBD", 0.95, 5)
ds_out = replace_values(ds_out, "fates_alloc_storage_cushion", 2.25, 2)
ds_out = replace_values(ds_out, "fates_allom_d2h1", [57.6, 57.6, 1])
ds_out = replace_values(ds_out, "fates_allom_d2h2", [0.74, 0.74, 1])
ds_out = replace_values(ds_out, "fates_allom_la_per_sa_int", 1000, 2)
ds_out = replace_values(ds_out, "fates_allom_sai_scaler", 0.0012, 2)
ds_out = replace_values(ds_out, "fates_fire_alpha_SH", [0.1487, 0.06, 1])
ds_out = replace_values(ds_out, "fates_fire_bark_scaler", [0.0301, 0.1085], [0, 1])
ds_out = replace_values(ds_out, "fates_fire_crown_kill", [1, 0.05, 1])
ds_out = replace_values(ds_out, "fates_grperc", [0.3, 0.3], [0, 1])
ds_out = replace_values(ds_out, "fates_leaf_vcmax25top", [41, 40], [0, 2])
ds_out = replace_values(ds_out, "fates_mort_hf_sm_threshold", [0.025, 0.025], [0, 1])
ds_out = replace_values(ds_out, "fates_wood_density", [0.6305, 0.695, 0.01])
ds_out = replace_values(ds_out, "fates_turnover_branch", [75, 75, 0.3208])
ds_out = replace_values(ds_out, "fates_allom_h2cd1", [0.33, 0.1], [0, 1])
ds_out = replace_values(ds_out, "fates_turnover_leaf", [[1.4025, 1.4025, 0.3208]])
ds_out = replace_values(ds_out, "fates_recruit_height_min", 0.5, 2)
ds_out = replace_values(ds_out, "fates_recruit_init_density", 20, 2)
ds_out = replace_values(ds_out, "fates_stoich_nitr", [0.02675, 0.02675, 0.16], [0, 1, 2])
ds_out = replace_values(ds_out, "fates_allom_d2bl1", [0.12668, 0.12668, 0.000964])
ds_out = replace_values(ds_out, "fates_allom_d2bl2", [1.2813, 1.2813, 1.9492])
ds_out = replace_values(ds_out, "fates_allom_d2ca_coefficient_max", [0.76865, 0.76865, 0.03])
ds_out = replace_values(ds_out, "fates_allom_d2ca_coefficient_min", [0.76865, 0.76865, 0.03])
ds_out = replace_values(ds_out, "fates_allom_l2fr", [0.4863, 0.4863], [0, 1])
ds_out = replace_values(ds_out, "fates_leaf_slamax", [0.03992, 0.03992, 0.0135])
ds_out = replace_values(ds_out, "fates_leaf_slatop", [0.01996, 0.01996, 0.0135])
ds_out = replace_values(ds_out, "fates_mort_scalar_cstarvation", [0.02956, 0.02956, 0.2])
ds_out = replace_values(ds_out, "fates_recruit_seed_alloc", [0.046801, 0.046801], [0, 1])
ds_out = replace_values(ds_out, "fates_fire_cg_strikes", 0.1)
ds_out = replace_values(ds_out, "fates_fire_threshold", 25.0)
ds_out = replace_values(ds_out, "fates_comp_excln", -1.0)
ds_out = replace_values(ds_out, "fates_mort_disturb_frac", 0.5)
ds_out = replace_values(ds_out, "fates_mort_understorey_death", 1.0)
ds_out = replace_values(ds_out, "fates_allom_agb1", [0.0673, 0.0673, 0.000964])
ds_out = replace_values(ds_out, "fates_allom_agb2", [0.976, 0.976, 1.9462])
ds_out = replace_values(ds_out, "fates_allom_amode", [3, 3], [0, 1])
ds_out = replace_values(ds_out, "fates_allom_agb3", [np.nan, np.nan, 0])
ds_out = replace_values(ds_out, "fates_allom_agb4", [np.nan, np.nan, 0])
ds_out = replace_values(ds_out, "fates_allom_fmode", 2.0, "all")
ds_out = replace_values(ds_out, "fates_allom_agb_frac", 1, 2)
ds_out = replace_values(ds_out, "fates_allom_dbh_maxheight", [200.0, 200.0, 1.0])
ds_out = replace_values(ds_out, "fates_allom_hmode", [5.0, 5.0, 3.0])
ds_out = replace_values(ds_out, "fates_allom_lmode", [3.0, 3.0, 1.0])
ds_out = replace_values(ds_out, "fates_allom_fnrt_prof_a", [7, 7], [1, 2])
ds_out = replace_values(ds_out, "fates_allom_fnrt_prof_b", [1, 1], [1, 2])
ds_out = replace_values(ds_out, "fates_rad_leaf_clumping_index", 0.85, 2)
ds_out = replace_values(ds_out, "fates_phen_evergreen", [1, 1], [1, 2])
ds_out = replace_values(ds_out, "fates_phen_stress_decid", [0, 0], [1, 2])
ds_out = replace_values(ds_out, "fates_leaf_stomatal_intercept", 10000, 2)
ds_out = replace_values(ds_out, "fates_leaf_stomatal_slope_ballberry", 4, 2)
ds_out = replace_values(ds_out, "fates_leaf_stomatal_slope_medlyn", [4.1, 4.1, 4.1])
ds_out = replace_values(ds_out, "fates_mort_bmort", [0.01303514, 0.01303514, 0.01303514])
ds_out = replace_values(ds_out, "fates_stoich_phos", [0.0033, 24], [1, 5])
ds_out = replace_values(ds_out, "fates_allom_d2h3", [21.6, 21.6], [0, 1])
ds_out = replace_values(ds_out, "fates_allom_d2bl3", [np.nan, np.nan, 0])
ds_out = replace_values(ds_out, "fates_cnp_turnover_nitr_retrans", np.nan, "all")
ds_out = replace_values(ds_out, "fates_cnp_turnover_phos_retrans", np.nan, "all")
ds_out = replace_values(ds_out, "fates_rad_leaf_xl", 0.1, 1)
ds_out = replace_values(ds_out, "fates_rad_leaf_rhonir", [0.46, 0.46, 0.28])
ds_out = replace_values(ds_out, "fates_rad_leaf_rhovis", [0.11, 0.11, 0.05])
ds_out = replace_values(ds_out, "fates_rad_stem_rhonir", [0.49, 0.49, 0.39])
ds_out = replace_values(ds_out, "fates_rad_stem_rhovis", [0.21, 0.21, 0.16])
ds_out = replace_values(ds_out, "fates_rad_leaf_taunir", [0.33, 0.33, 0.4])
ds_out = replace_values(ds_out, "fates_rad_leaf_tauvis", [0.06, 0.06], [0, 1])
ds_out = replace_values(ds_out, "fates_rad_stem_taunir", 0.001, 2)
ds_out = replace_values(ds_out, "fates_rad_stem_tauvis", 0.001, 2)
ds_out = replace_values(ds_out, "fates_recruit_seed_dbh_repro_threshold", [150, 0.5], [1, 2])
ds_out = replace_values(ds_out, "fates_phen_flush_fraction", np.nan, "all")

ds_out.to_netcdf(file_out)
