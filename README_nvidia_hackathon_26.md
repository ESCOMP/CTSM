# Instructions for running Profiling for the NVIDIA/OpenACC Hackathon 2026

### General things:

- We have changes here to cime and ccs_config
- The testlist to use is: nvidia-hackathon-26
- ccs_config changes make derecho_nvhpc add more compiler analysis
- ccs_config changes also add derecho_nvhpc-prof to run the profiler in the run phase 
- You analyze the report file after you run your case

## How to profile a case:

There is extra output added when the build is done under derecho_nvhpc or derecho_nvhpc-prof so that we can assess what optimization is being done. The profiler will be run automatically when you run a case with derecho_nvhpc-prof. You then analyze the report generated when a case is run with the NVIDIA tool nsys.

``` shell
git clone --origin escomp -b nvidia-hackathon-26 https://github.com/ESCOMP/CTSM.git ctsm-nvidia-hackathon-26
cd ctsm-nvidia-hackathon-26
./bin/git-fleximod update
cd cime/scripts
./create_test SMS.f10_f10_mg37.I2000Clm60BgcCrop.derecho_nvhpc-prof.clm-crop --walltime 00:20:00 --no-build --test-id test_hackathon
cd $SCRATCH/SMS.f10_f10_mg37.I2000Clm60BgcCrop.derecho_nvhpc-prof.clm-crop.test_hackathon
# Build the model
./case.build
# Now examine the bld/lnd.bldlog.* file for the optimization information
# This will be at the end of it for all files

# Next submit the job into the batch queing system and wait for it to finish
./case.submit
#
# After the case runs it will have run the profiler and created a report file
#
# Get the module load and environment setup for later
source .env_mach_specific.csh # for cshell
. .env_mach_specific.sh # For bash
# Use the GUI nsys-ui to analyze the report file it generated
ls run/report1.nsys-rep
nsys-ui  # Open the above report file there
```
## Some things we'd like to do:

- Fix some of the bugs that make changing the build for a case more difficult
- Bring in some of the longer term methods to make the build both documented, more flexible and easier to customize
- Bring in a way to do the profiling as a standard part of ccs_config
