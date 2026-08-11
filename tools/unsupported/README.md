# $CTSMROOT/tools/unsupported/README.md

Updated 2026/08/11 slevis
Written 2019/01/24 ekluzek

## Purpose
This directory is for CTSM users to contribute scripts that others may find useful. Contributors should always add documentation within their scripts. We consider these scripts unsupported, and they may or may not work. Also settings within the scripts may be hardwired, e.g. paths may assume NCAR's directory structures.

## Basic instructions
Python scripts require the following settings before running on NCAR's supercomputers. The `py_env_create` step is required once, unless one needs to update their `ctsm_pylib` environment:

``` shell
module load conda
../../py_env_create
conda activate ctsm_pylib
```

ncl scripts work as follows on NCAR's supercomputers:
``` shell
module load ncl
ncl script_name.ncl
```

## Note
The separate repository `https://github.com/NCAR/CMIP7_inputdata_processing` includes contributions from the LMWG, each with README instructions.

