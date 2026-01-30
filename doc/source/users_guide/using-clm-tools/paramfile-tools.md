
# Tools for working with parameter files

This guide describes the features and usage of the `query_paramfile` and `set_paramfile` tools, located in `tools/param_utils/`. These utilities help users inspect and modify CLM parameter files.

Note that you need to have the `ctsm_pylib` conda environment activated to use these tools. See Sect. :numref:`using-ctsm-pylib` for more information.

## `query_paramfile`
**Purpose:** Print the values of one or more parameters from a CTSM parameter file (netCDF format).

**Features:**
- Print values for specified parameters or all.
- Optionally filter output by Plant Functional Types (PFTs) for PFT-specific parameters.

For more information, do `tools/param_utils/query_paramfile --help`.


### Example usage

Print all variables in a parameter file:
```bash
tools/param_utils/query_paramfile -i paramfile.nc
```

Print specific variables:
```bash
tools/param_utils/query_paramfile -i paramfile.nc jmaxha jmaxhd
```

Print values for specific PFTs:
```bash
tools/param_utils/query_paramfile -i paramfile.nc -p needleleaf_evergreen_temperate_tree,c4_grass medlynintercept medlynslope
```

## `set_paramfile`
**Purpose:** Change values of one or more parameters in a CTSM parameter file (netCDF format).

**Features:**
- Modify parameter values for all or selected PFTs.
- Optionally drop PFTs not specified.
- Set parameter values to fill (missing) values using `nan`.
- Ensures safe file handling and checks for argument validity.

Note that the output file must not already exist.

For more information, do `tools/param_utils/set_paramfile --help`.

### Example usage

Change a scalar parameter:
```bash
tools/param_utils/set_paramfile -i paramfile.nc -o output.nc jmaxha=51000
```

Change a one-dimensional parameter (`mimics_fmet` has the `segment` dimension, length 4):
```bash
tools/param_utils/set_paramfile -i paramfile.nc -o output.nc mimics_fmet=0.1,0.2,0.3,0.4
```

Change a one-dimensional parameter to be all one value (`mxmat` has the `pft` dimension, length 79):
```bash
tools/param_utils/set_paramfile -i paramfile.nc -o output.nc mxmat=360
```

Change a parameter for specific PFTs:
```bash
tools/param_utils/set_paramfile -i paramfile.nc -o output.nc -p needleleaf_evergreen_temperate_tree,c4_grass medlynintercept=99.9,100.1 medlynslope=2.99,1.99 mxmat=199
```

Set a parameter to the fill value:
```bash
tools/param_utils/set_paramfile -i paramfile.nc -o output.nc -p needleleaf_evergreen_temperate_tree,c4_grass fleafcn=nan,nan
```

## `compare_paramfiles`
**Purpose:** Compare values of one or more parameters in two CTSM parameter files (netCDF format). Compares both "raw" values and values after masking and scaling (applying `_FillValue` to convert missing values to NaN, multiplying by `scale_factor`, and adding `add_offset`).

### Example usage
Results truncated.

```shell
file0=/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/paramdata/clm_params78pftModDates.c130821.nc
file1=/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/paramdata/ctsm60_params.c251124.nc

tools/param_utils/compare_paramfiles $file0 $file1
```

```
File 0: /glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/paramdata/clm_params78pftModDates.c130821.nc
File 1: /glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/paramdata/ctsm60_params.c251124.nc

Variable(s) present in File 0 but not File 1:
   cn_s1_bgc
   cn_s2_bgc
   ...

Variable(s) present in File 1 but not File 0:
   C2_liq_Brun89
   FUN_fracfixers
   ...

...

allconss:
   Values differ:
      [pft 71 (miscanthus)] 0.0 → 2.0 (raw)
                            nan → 2.0 (masked/scaled)
      [pft 72 (irrigated_miscanthus)] 0.0 → 2.0 (raw)
                                      nan → 2.0 (masked/scaled)
      [pft 73 (switchgrass)] 0.0 → 2.0 (raw)
                             nan → 2.0 (masked/scaled)
      [pft 74 (irrigated_switchgrass)] 0.0 → 2.0 (raw)
                                       nan → 2.0 (masked/scaled)

...

atmch4:
   Attribute(s) present in File 1 but not File 0:
       _FillValue: nan
   Dimension names differ: File 0: ['allpfts']
                           File 1: []

baset:
   Attribute(s) present in File 1 but not File 0:
       _FillValue: nan
   Values differ:
      [pft 71 (miscanthus)] 0.0 → 8.0
      [pft 72 (irrigated_miscanthus)] 0.0 → 8.0
      [pft 73 (switchgrass)] 0.0 → 8.0
      [pft 74 (irrigated_switchgrass)] 0.0 → 8.0

...
```

For more information, do `tools/param_utils/compare_paramfiles --help`.