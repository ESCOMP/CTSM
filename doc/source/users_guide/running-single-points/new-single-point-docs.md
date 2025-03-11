# New single-point documentation

## 1. Subset the data

```shell
cd tools/site_and_regional
conda run -n ctsm_pylib ./subset_data point --lat 45.402252 --lon 267.201915 --site CC_C3_crujra --create-surface --create-datm --datm-syr 2000 --datm-eyr 2022 --dompft 13 --create-user-mods --outdir /glade/derecho/scratch/$USER/CC_crujra_subset_C3
```

- Ran into climate data issues. Had wanted to use CRU-JRA. Manually changed something to make that happen, but we can't figure out what. File an issue to make this easier (and to make default CRU-JRA)
- Note that `subset_data` doesn't work yet when sampling datm for a region: [Issue #2110: Don't allow users to try to subset datm for a region with subset_data: doesn't work yet](https://github.com/ESCOMP/CTSM/issues/2110)

## 2. Make case
```shell
./create_newcase --case /glade/derecho/scratch/krocci/CC_case_C4_test --res CLM_USRDAT --compset I2000Clm60Bgc --run-unsupported --user-mods-dirs /glade/derecho/scratch/krocci/CC_crujra_subset/user_mods/
```

`./case.setup` on Derecho will automatically set queue to `develop` and walltime to 1 hour. So you might need to modify that. The max walltime for develop is 1 hour.
