How to run CTSM on Mistral?

1. Load extra modules
> module load python/2.7.12
> module load intel
> module load openmpi/2.0.2p2_hpcx-intel14
> module load esmf
> module load nco/4.7.5-gcc64
> module load ncl
> module load gcc/6.4.0
> module load cmake

2. Set up additionnal LD_LIBRARY_PATH paths
> export LD_LIBRARY_PATH="/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-openmpi2-intel14/lib:/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14/lib:$LD_LIBRARY_PATH"

3. Copy /pf/a/a271098/.cime/config_machines.xml and /pf/a/a271098/.cime/config_batch.xml into your home in a .cime directory

4. Use /work/aa0049/a271098/CTSM/cime/scripts/create_newcase in one of your home directory. 

From there, consult the documentation on this website: https://escomp.github.io/ctsm-docs/versions/master/html/index.html. You can use my script /pf/a/a271098/CTSM_runs/default/newcase_test.bash to help you.

All the results of your simulation will be saved in my work directory /work/aa0049/a271098/archive and /work/aa0049/a271098/scratch. Do not replace a case already existing! A request message will show up when your create the case:
Directory /work/aa0049/a271098/CTSM/scratch/I1850CLM50_001/bld already exists, (r)eplace, (a)bort, or (u)se existing?

If you want to store the results of your simulations in your work directory, follow the second part.

Enjoy!


How to install CTSM in your work directory /work/aa0049/$USER/

1. Clone CTSM github in your home directory
> git clone https://github.com/ESCOMP/CTSM.git
Watch out: only the version 5a0ba10 is working on mistral. More recent versions still need to be tested. To use this version do this following line in /work/aa0049/$USER/CTSM:
> git checkout 5a0ba100a094c7152fa0f7247d44f2e88a42823f

2. Unpack the model
> ./manage_externals/checkout_externals

3. Copy my CLMBuildNamelist.pm in your directory /work/aa0049/$USER/CTSM/bld/
> cp /work/aa0049/a271098/CTSM/bld/CLMBuildNamelist.pm .

5. Copy my config_inputdata.xml in your directory /work/aa0049/$USER/CTSM/cime/config/cesm
> cp /work/aa0049/a271098/CTSM/cime/config/cesm/config_inputdata.xml .

6. Start from "How to run CTSM on Mistral?"
- When you copy config_machines.xml, replace every a271098 by $USER

Enjoy!


If you need help or want to do a regional simulation, contact adamseau@awi.de.
