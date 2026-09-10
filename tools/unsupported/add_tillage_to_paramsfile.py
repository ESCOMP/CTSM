import xarray as xr
import numpy as np
import os
import subprocess
import argparse
import sys

def make_dataarray(np_array, decomp_model, ntill_intensities_max, ndecomp_pools_max, ntill_stages_max):
    intensities_dim = xr.IndexVariable(
        dims = "ntill_intensities_max",
        data = np.arange(ntill_intensities_max),
        )
    pools_dim = xr.IndexVariable(
        dims = "ndecomp_pools_max",
        data = np.arange(ndecomp_pools_max),
        )
    stages_dim = xr.IndexVariable(
        dims = "ntill_stages_max",
        data = np.arange(ntill_stages_max),
        )
    
    # Name DataArray
    if decomp_model.lower() == "mimics":
        da_name = "mimics"
    elif decomp_model.lower() == "century":
        da_name = "bgc"
    da_name += "_till_decompk_multipliers"
        
    da = xr.DataArray(
        data = np_array,
        dims = {
            "ntill_intensities_max": intensities_dim,
            "ndecomp_pools_max": pools_dim,
            "ntill_stages_max": stages_dim
            },
        name = da_name,
        attrs = {
            "long_name": f"Value by which decomp_k should be multiplied during tillage with {decomp_model} soil",
            "units": "unitless",
        }
    )
    
    # netCDF variable needs dimensions reversed from how they're specified in code
    da = da.transpose()
    
    return da

def main(file_in,
         file_out):
    # Get git info
    thisDir = os.path.dirname(__file__)
    git_status = subprocess.run(
        ["git", "status"],
        capture_output=True,
    )
    git_status = git_status.stdout.decode()
    repo_is_clean = "working tree clean" in git_status
    if not repo_is_clean:
        print("WARNING: Repo not clean; will not save params file.")
        print(git_status)
    git_log = subprocess.run(
        ["git", "log", "-1"],
        capture_output=True,
    )
    git_log = git_log.stdout.decode()


    # Set up dimensions
    ds0 = xr.open_dataset(file_in)
    ntill_intensities_max = 2
    ndecomp_pools_max = ds0.dims["ndecomp_pools_max"]
    ntill_stages_max = 3
    tillage_shape_ips = (ntill_intensities_max, ndecomp_pools_max, ntill_stages_max)
    tillage_shape_ps = (ndecomp_pools_max, ntill_stages_max)


    # Fill CENTURY array
    # Define pool indices
    i_litr_min = 0  # 1 in FORTRAN, but Python is 0-indexed
    i_met_lit = i_litr_min
    i_cel_lit = i_met_lit + 1
    i_lig_lit = i_cel_lit + 1
    i_act_som = i_lig_lit + 1
    i_slo_som = i_act_som + 1
    i_pas_som = i_slo_som + 1
    tillage_century = np.full(tillage_shape_ips, 1.0)
    tillage_century_lo = np.full(tillage_shape_ps, 1.0)
    tillage_century_lo[i_act_som,:] = np.array([1.0, 1.0, 1.0])
    tillage_century_lo[i_slo_som,:] = np.array([3.0, 1.6, 1.3])
    tillage_century_lo[i_pas_som,:] = np.array([3.0, 1.6, 1.3])
    tillage_century_lo[i_cel_lit,:] = np.array([1.5, 1.5, 1.1])
    tillage_century_lo[i_lig_lit,:] = np.array([1.5, 1.5, 1.1])
    tillage_century[0,:,:] = tillage_century_lo
    tillage_century_hi = np.full(tillage_shape_ps, 1.0)
    tillage_century_hi[i_act_som,:] = np.array([1.2, 1.0, 1.0])
    tillage_century_hi[i_slo_som,:] = np.array([4.8, 3.5, 2.5])
    tillage_century_hi[i_pas_som,:] = np.array([4.8, 3.5, 2.5])
    tillage_century_hi[i_cel_lit,:] = np.array([1.8, 1.5, 1.1])
    tillage_century_hi[i_lig_lit,:] = np.array([1.8, 1.5, 1.1])
    tillage_century[1,:,:] = tillage_century_hi


    # Fill MIMICS array
    i_litr_min = 1
    i_met_lit = i_litr_min
    i_str_lit = i_met_lit + 1
    i_avl_som = i_str_lit + 1
    i_chem_som = i_avl_som + 1
    i_phys_som = i_chem_som + 1
    tillage_mimics = np.full(tillage_shape_ips, 1.0)
    tillage_mimics[:,i_avl_som,:] = tillage_century[:,i_act_som,:]
    tillage_mimics[:,i_chem_som,:] = tillage_century[:,i_slo_som,:]
    tillage_mimics[:,i_phys_som,:] = tillage_century[:,i_pas_som,:]
    if not np.array_equal(tillage_century[:,i_cel_lit,:], tillage_century[:,i_lig_lit,:]):
        raise RuntimeError("How to combine 2 CENTURY litter pools into 1 MIMICS litter pool?")
    tillage_mimics[:,i_str_lit,:] = tillage_century[:,i_cel_lit,:]

    # Make DataArrays
    tillage_century_da = make_dataarray(
        tillage_century, "CENTURY",
        ntill_intensities_max, ndecomp_pools_max, ntill_stages_max)
    tillage_mimics_da = make_dataarray(
        tillage_mimics, "MIMICS",
        ntill_intensities_max, ndecomp_pools_max, ntill_stages_max)

    if not repo_is_clean:
        raise RuntimeError("Clean up git repo before trying to save!")

    ds1 = ds0.copy()
    ds0.close()

    ds1[tillage_century_da.name] = tillage_century_da
    ds1[tillage_mimics_da.name] = tillage_mimics_da
    ds1.attrs['latest_git_log'] = git_log

    ds1.to_netcdf(file_out, format="NETCDF3_CLASSIC")


if __name__ == "__main__":
    ###############################
    ### Process input arguments ###
    ###############################
    parser = argparse.ArgumentParser(
        description="Adds tillage parameters to a CLM parameter file (netCDF)."
    )

    # Define arguments
    parser.add_argument(
        "-i",
        "--input-file",
        help="Parameter file (netCDF) to which you wish to add tillage parameters.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-file",
        help="Output parameter file.",
        type=str,
        required=True,
    )

    # Get arguments
    args = parser.parse_args(sys.argv[1:])

    ###########
    ### Run ###
    ###########
    main(os.path.realpath(args.input_file),
         os.path.realpath(args.output_file),
    )