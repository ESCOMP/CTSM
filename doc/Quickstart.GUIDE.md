# Quick-Start to Using NUOPC Scripts for ctsm6_0
---

## Assumptions:

You want to use Derecho with ctsm6_0 BGC
to do a CTSM simulation with data atmosphere and the
latest CRUJRA atm forcing files and settings. You also want to cycle
the CRUJRA atm data between 1950 to 2010 and you want to run at
0.9x1.25 degree resolution.

## Process:

   ### Create the case

``` shell
   cd cime/scripts

   ./create_newcase --case <testcase> --mach derecho --res f09_t232 -compset I2000Clm60BgcCrop
   # (./create_newcase -help -- to get help on the script)

   # Setup the case

   cd <testcase>
   ./xmlchange id1=val1,id2=val2  # to make changes to any settings in the env_*.xml files
   ./case.setup
   #(./case.setup -help -- to get help on the script)

   # Add any namelist changes to the user_nl_* files

   $EDITOR user_nl_*

   # Compile the code

   ./case.build

   # Submit the run

   ./case.submit

```

## Information on Compsets:

     "I" compsets are the ones with CTSM and NUOPC driver and CDEPS data models without ice and ocean.
     Most of the "I" compsets are for clm6_0 physics and use the CRUJRA-2024 data with solar following
     the cosine of solar zenith angle, precipitation constant, and other
     variables linearly interpolated in time (and with appropriate time-stamps on
     the date). Previous CMIP6 simulations with clm5_0 physics used GSWP3 atmospheric forcing.

     To get a list of the "I" compsets use the `query_config` utility in `cime/scripts`

``` shell
    cd cime/scripts
    ./query_config --compsets clm
```

## Automatically resubmitting jobs:

   After doing a short simulation that you believe is correct

``` shell
   ./xmlchange CONTINUE_RUN=TRUE

   # Change RESUBMIT to number greater than 0, and CONTINUE_RUN to TRUE...

   ./case.submit
```
