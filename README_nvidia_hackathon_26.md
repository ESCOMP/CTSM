# Instructions for running Profiling for the NVIDIA/OpenACC Hackathon 2026

### General things:

- We have changes here to cime and ccs_config
- The testlist to use is: nvidia-hackathon-26
- ccs_config changes make derecho_nvhpc add more compiler analyisis
- Using derecho_nvhpc also runs profiler
- Analyze the report file after you run your case

## How to profile a case:

There is extra output added when the build is done under derecho_nvhpc so that we can assess what optimization is being done. The profiler will be run automatically when you run a case with derecho_nvhpc. You then analyze the report generated when a case is run with the NVIDIA tool nsys.

``` shell
git clone --origin escomp -b nvidia-hackathon-26 https://github.com/ESCOMP/CTSM.git
./bin/git-fleximod update
cd cime/scripts
./create_test SMS_D.f10_f10_mg37.I2000Clm60BgcCrop.derecho_nvhpc.clm-crop --walltime 00:20:00 --no-build --test_id test_hackathon
cd $SCRATCH/SMS_D.f10_f10_mg37.I2000Clm60BgcCrop.derecho_nvhpc.clm-crop.test_hackathon
./case.build
# Get the module load and environment setup for later
source .env_mach_specific.csh # for cshell
. .env_mach_specific.sh # For bash
# Now examine the bld/lnd.bldlog.* file for the vectorization information
./case.submit
# After the case runs it will have run the profiler and created a report file
# Use nsys to analyze the report file
nsys analyze run/<report_file>
```

## Some things we'd like to do:

- Have a seperate compiler for the profiling option (nvhpc-prof)
- Fix some of the bugs that make changing the build for a case more difficult
- Bring in some of the longer term methods to make the build both documented, more flexible and easier to customize
- Bring in a way to do the profiling as a standard part of ccs_config
