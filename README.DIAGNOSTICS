# QuickStart to Running Diagnostics of your Case with CUPiD (CESM Unified Postprocessing and Diagnostics)

[!NOTE]
Land Diagnostics Cheat Sheet is here:

This document has more details on more options and such

https://docs.google.com/document/d/1ziZWGgaj9FxgR6WRHyCAZnzIc6-4ezqfUDDpQv4MzUI/edit?tab=t.0

## Initial setup steps to do first time

This step is something you need to do the first time you want to run CUPiD for a case.
And you only need to do it again if the CUPiD environment changes.

``` shell
# Setup the conda environments
# TODO: Make this use py_env_create?
mamba env create --yes -f tools/cupid/environments/cupid-infrastructure.yml
mamba env create --yes -f tools/cupid/environments/cupid-analysis.yml
# Check that the environment is valid
conda activate cupid-infrastructure
which cupid-diagnostics
# Should return something like:
# $HOME/conda-envs/cupid-infrastructure/bin/cupid-diagnostics
# If it returns "Command not found." something is wrong with the environment
```

## Create your case using the CUPiD user-mod

This is similar to setting up any case, such as documented in the README and Quickstart guides.

``` shell
./create_newcase --case testIwCUPiD --res f09_t232 --compset I2000Clm60BgcCrop --user-mods-dirs clm-CUPiD
```

## Run CUPID in your case

After st_archive has run do the following. With RUN_POSTPROCESSING set to TRUE this will happen
with each case submission automatically. But, if you want to run it separately...


``` shell
./case.submit --only-job case.cupid
```