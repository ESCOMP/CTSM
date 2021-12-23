# $CTSMROOT/README.NUOPC_driver

CTSM now by default uses the NUOPC based CMEPS driver!


## What's new?

MESH Files:
Mesh files to describe grids are new in both the driver namelist and for example in any
streams files.
Full ESMF Library is used:
The full ESMF Library is used and required to be built in order to run the model.
Single Point cases:
Single point cases can now set their location using PTS_LAT and PTS_LON.

## What's removed?

Domain files are no longer used. And mapping for regriding is created on the fly
rather than using fixed mapping files in almost all cases. Runoff mapping files
still need to be generated offline.

## What files change?

rpointer.drv becomes rpointer.cpl
cpl.log.* files get's split into med.log.* and drv.log.*
user_datm.streams.txt.* file changes goes into the user_nl_datm_streams files
datm.streams.txt.* files are all in one file called datm.streams.xml

## What XML variables change in your case?

DATM_CLMNCEP_YR_* variables change to DATM_YR_*

## New obscure options:

ESMF_AWARE_THREADING --- ESMF is aware of threading (can have differing number of threads in components)
CREATE_ESMF_PET_FILES -- Create output log files from ESMF for each Processor (PET)
ESMF_VERBOSITY_LEVEL --- Verbosity level for ESMF logging
ESMF_PROFILING_LEVEL --- Verbosity level for ESMF profiling

nuopc.runseq is a text file that determines how the driver operates. You can change the operation
by having an updated copy in your case directory.


## What if I want to use the previous MCT driver?

The MCT driver will be available for sometime going forward, but
new development won't go into it, and it will eventually be removed.
But, if you have to...
Use the "--driver mct" command line option to create_newcase
You can set COMP_INTERFACE in a case as well, but it won't create it with everything needed
so we recommend setting up a case from scratch.


For more notes see:

https://docs.google.com/presentation/d/1yjiKSEV53JDAJbYxhpY2T9GTxlWFzQAn
