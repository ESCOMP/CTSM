These gridcell files were generated according to the procedure given in the `representative-hillslopes` submodule's `shell_scripts/README.md`, plus the following:
```bash
mv chunk_07_HAND_4_col_hillslope_geo_params_Trapezoidal_MERIT_j_001_i_001.nc chunk_06_HAND_4_col_hillslope_geo_params_Trapezoidal_MERIT_j_001_i_001.nc
mv chunk_25_HAND_4_col_hillslope_geo_params_Trapezoidal_MERIT_j_004_i_004.nc chunk_26_HAND_4_col_hillslope_geo_params_Trapezoidal_MERIT_j_004_i_004.nc
```

I.e., chunk 7's gridcell was lumped in with chunk 6, and chunk 25 was renamed to chunk 26. This change exercises the `combine_gridcell_files` script a bit more than having one gridcell per chunk and no skipped chunks.
