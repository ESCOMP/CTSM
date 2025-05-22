.. include:: ../substitutions.rst

.. _running-single-point-datasets:

******************************************
 Running Single Point Configurations
******************************************

In addition to running with the outputs of ``subset_data`` (Sect. :numref:`single_point_subset_data`), CLM supports running using single-point or regional datasets that are customized to a particular region. CLM supports a a small number of out-of-the-box single-point and regional datasets. However, users can create their own dataset.

To get the list of supported dataset resolutions do this:
::

   > cd $CTSMROOT/doc
   > ../bld/build-namelist -res list

Which results in the following:
::

   CLM build-namelist - valid values for res (Horizontal resolutions
   Note: 0.5x0.5, 5x5min, 10x10min, 3x3min and 0.33x0.33 are only used for CLM tools):
         Values: default 512x1024 360x720cru 128x256 64x128 48x96 32x64 8x16 94x192  \
                 0.23x0.31 0.47x0.63 0.9x1.25 1.9x2.5 2.5x3.33 4x5 10x15 5x5_amazon 1x1_tropicAtl \
                 1x1_vancouverCAN 1x1_mexicocityMEX 1x1_asphaltjungleNJ 1x1_brazil 1x1_urbanc_alpha 1x1_numaIA  \
                 1x1_smallvilleIA 0.5x0.5 3x3min 5x5min 10x10min 0.33x0.33 ne4np4 ne16np4 ne30np4 ne60np4  \
                 ne120np4 ne240np4 wus12 us20
         Default = 1.9x2.5
        (NOTE: resolution and mask and other settings may influence what the default is)

The resolution names that have an underscore in them ("_") are all single-point or regional resolutions.

.. note:: When running a single point, the number of processors is automatically set to one, which is the only value allowed.

.. warning::
   Just like running with the outputs from ``subset_data`` (Sect. :numref:`single_point_subset_data`), by default these setups sometimes run with ``MPILIB=mpi-serial`` (in the ``env_build.xml`` file) turned on, which allows you to run the model interactively. On some machines this mode is NOT supported and you may need to change it to FALSE before you are able to build.

.. _single-point-global-climate:

Single-point runs with global climate forcings
==============================================

Example: Use global forcings at a site without its own special forcings
-----------------------------------------------------------------------

This example uses the single-point site in Brazil.
::

   > cd cime/scripts
   > set SITE=1x1_brazil
   > ./create_newcase -case testSPDATASET -res $SITE -compset I2000Clm50SpGs
   > cd testSPDATASET

Then setup, build and run normally.

Example: Use global forcings at a site WITH its own special forcings
--------------------------------------------------------------------

The urban Mexico City test site has its own atmosphere forcing data (see Sect. :numref:`single-point-with-own-forcing`). To ignore that and run it with the default global forcing data, but over the period for which its own forcing data is provided, do the following:

::

   > cd cime/scripts
   # Set a variable to the site you want to use (as it's used several times below)
   > set SITE=1x1_mexicocityMEX
   > ./create_newcase -case testSPDATASET -res $SITE -compset I1PtClm50SpGs
   > cd testSPDATASET

(Note the use of ``I1Pt`` instead of ``I2000`` as in the example above.) Then setup, build and run normally.

.. _single-point-with-own-forcing:

Supported single-point runs for sites with their own atmospheric forcing
========================================================================

Of the supported single-point datasets we have three that also have atmospheric forcing data that go with them: Mexico City (Mexico), Vancouver, (Canada, British Columbia), and ``urbanc_alpha`` (test data for an Urban inter-comparison project). Mexico city and Vancouver also have namelist options in the source code for them to work with modified urban data parameters that are particular to these locations. To turn on the atmospheric forcing for these datasets, you set the ``env_run.xml DATM_MODE`` variable to ``CLM1PT``, and then the atmospheric forcing datasets will be used for the point picked. If you use one of the compsets that has "I1Pt" in the name that will be set automatically.

.. todo::
    Update the below, as ``queryDefaultNamelist.pl`` no longer exists.

When running with datasets that have their own atmospheric forcing you need to be careful to run over the period that data is available. If you have at least one year of forcing it will cycle over the available data over and over again no matter how long of a simulation you run. However, if you have less than a years worth of data (or if the start date doesn't start at the beginning of the year, or the end date doesn't end at the end of the year) then you won't be able to run over anything but the data extent. In this case you will need to carefully set the ``RUN_STARTDATE``, ``START_TOD`` and ``STOP_N/STOP_OPTION`` variables for your case to run over the entire time extent of your data. For the supported data points, these values are in the XML database and you can use the ``queryDefaultNamelist.pl`` script to query the values and set them for your case (they are set for the three urban test cases: Mexicocity, Vancouver, and urbanc_alpha).

Example: Use site-specific atmospheric forcings
-----------------------------------------------
In this example, we show how to use the atmospheric forcings specific to the Vancouver, Canada point.
::

   > cd cime/scripts

   # Set a variable to the site you want to use (as it's used several times below)
   > set SITE=1x1_vancouverCAN

   # Create a case at the single-point resolutions with their forcing
   > ./create_newcase -case testSPDATASETnAtmForcing -res $SITE -compset I1PtClm50SpGs
   > cd testSPDATASETnAtmForcing

   # Figure out the start and end date for this dataset
   # You can do this by examining the datafile.
   > set STOP_N=330
   > set START_YEAR=1992
   > set STARTDATE=${START_YEAR}-08-12
   > @ NDAYS = $STOP_N / 24
   > ./xmlchange RUN_STARTDATE=$STARTDATE,STOP_N=$STOP_N,STOP_OPTION=nsteps

   # Set the User namelist to set the output frequencies of the history files
   # Setting the stdurbpt use-case option create three history file streams
   # The frequencies and number of time-samples needs to be set
   > cat << EOF > user_nl_clm
   hist_mfilt = $NDAYS,$STOP_N,$STOP_N
   hist_nhtfrq = -1,1,1
   EOF

   > ./case.setup

.. warning:: If you don't set the start-year and run-length carefully as shown above the model will abort with a "dtlimit error" in the atmosphere model. Since, the forcing data for this site (and the MexicoCity site) is less than a year, the model won't be able to run for a full year. The ``1x1_urbanc_alpha`` site has data for more than a full year, but neither year is complete hence, it has the same problem (see the problem for this site above).

.. _creating-your-own-singlepoint-dataset:

Creating your own single-point dataset
===================================================

The following provides an example of setting up a case using ``CLM_USRDAT_NAME`` where you rename the files according to the ``CLM_USRDAT_NAME`` convention. We have an example of such datafiles in the repository for a specific region over Alaska (actually just a sub-set of the global f19 grid).

Example: Using CLM_USRDAT_NAME to run a simulation using user datasets for a specific region over Alaska
-----------------------------------------------------------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case my_userdataset_test -res CLM_USRDAT -compset I2000Clm50BgcCruGs
   > cd my_userdataset_test/
   > set GRIDNAME=13x12pt_f19_alaskaUSA
   > set LMASK=gx1v6
   > ./xmlchange CLM_USRDAT_NAME=$GRIDNAME,CLM_BLDNML_OPTS="-mask $LMASK"
   > ./xmlchange ATM_DOMAIN_FILE=domain.lnd.${GRIDNAME}_$LMASK.nc
   > ./xmlchange LND_DOMAIN_FILE=domain.lnd.${GRIDNAME}_$LMASK.nc

   # Make sure the file exists in your $CSMDATA or else use svn to download it there
   > ls $CSMDATA/lnd/clm2/surfdata_map/surfdata_${GRIDNAME}_simyr2000.nc

   # If it doesn't exist, comment out the following...
   #> setenv SVN_INP_URL https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/
   #> svn export $SVN_INP_URL/lnd/clm2/surfdata_map/surfdata_${GRIDNAME}_simyr2000.nc $CSMDATA/lnd/clm2/surfdata_map/surfdata_${GRIDNAME}_simyr2000.nc
   > ./case.setup

The first step is to create the domain and surface datasets using the process outlined in :ref:`using-clm-tools-section`. Below we show an example of the process.

Example: Creating a surface dataset for a single point
---------------------------------------------------------------------
.. todo::
    Update the below, as ``mksurfdata.pl`` no longer exists and domain files aren't needed with nuopc.

::

   # set the GRIDNAME and creation date that will be used later
   > setenv GRIDNAME 1x1_boulderCO
   > setenv CDATE    `date +%y%m%d`
   # Create the SCRIP grid file for the location and create a unity mapping file for it.
   > cd $CTSMROOT/tools/mkmapdata
   > ./mknoocnmap.pl -p 40,255 -n $GRIDNAME
   # Set pointer to MAPFILE just created that will be used later
   > setenv MAPFILE `pwd`/map_${GRIDNAME}_noocean_to_${GRIDNAME}_nomask_aave_da_${CDATE}.nc
   # create the mapping files needed by mksurfdata_esmf.
   > cd ../.././mkmapdata
   > setenv GRIDFILE ../mkmapgrids/SCRIPgrid_${GRIDNAME}_nomask_${CDATE}.nc
   > ./mkmapdata.sh -r $GRIDNAME -f $GRIDFILE -t regional
   # create the domain file
   > cd ../../../../tools/mapping/gen_domain_files/src
   > ../../../scripts/ccsm_utils/Machines/configure -mach cheyenne -compiler intel
   > gmake
   > cd ..
   > setenv OCNDOM domain.ocn_noocean.nc
   > setenv ATMDOM domain.lnd.{$GRIDNAME}_noocean.nc
   > ./gen_domain -m $MAPFILE -o $OCNDOM -l $ATMDOM
   # Save the location where the domain file was created
   > setenv GENDOM_PATH `pwd`
   # Finally create the surface dataset
   > cd ../../../../lnd/clm/tools/|version|/mksurfdata_esmf/src
   > gmake
   > cd ..
   > ./mksurfdata.pl -r usrspec -usr_gname $GRIDNAME -usr_gdate $CDATE

The next step is to create a case that points to the files you created above. We will still use the ``CLM_USRDAT_NAME`` option as a way to get a case setup without having to add the grid to scripts.

Example: Setting up a case from the single-point surface dataset just created
--------------------------------------------------------------------------------------------

.. todo::
    Change this to provide instructions for a CTSM checkout instead of a CESM one.

.. todo::
    Update the below, as domain files aren't needed with nuopc.

::

   # First setup an environment variable that points to the top of the CESM directory.
   > setenv CESMROOT <directory-of-path-to-main-cesm-directory>
   # Next make sure you have a inputdata location that you can write to
   # You only need to do this step once, so you won't need to do this in the future
   > setenv MYCSMDATA $HOME/inputdata     # Set env var for the directory for input data
   > ./link_dirtree $CSMDATA $MYCSMDATA
   # Copy the file you created above to your new $MYCSMDATA location following the CLMUSRDAT
   # naming convention (leave off the creation date)
   > cp $CESMROOT/$CTSMROOT/tools/mksurfdata_esmf/surfdata_${GRIDNAME}_simyr1850_$CDATE.nc \
   $MYCSMDATA/lnd/clm2/surfdata_map/surfdata_${GRIDNAME}_simyr1850.nc
   > cd $CESMROOT/cime/scripts
   > ./create_newcase -case my_usernldatasets_test -res CLM_USRDAT -compset I1850Clm50BgcCropCru \
   -mach cheyenne_intel
   > cd my_usernldatasets_test
   > ./xmlchange DIN_LOC_ROOT=$MYCSMDATA
   # Set the path to the location of gen_domain set in the creation step above
   > ./xmlchange ATM_DOMAIN_PATH=$GENDOM_PATH,LND_DOMAIN_PATH=$GENDOM_PATH
   > ./xmlchange ATM_DOMAIN_FILE=$ATMDOM,LND_DOMAIN_FILE=$ATMDOM
   > ./xmlchange CLM_USRDAT_NAME=$GRIDNAME
   > ./case.setup

.. note:: With this and previous versions of the model we recommended using ``CLM_USRDAT_NAME`` as a way to identify your own datasets without having to enter them into the XML database. This has the down-side that you can't include creation dates in your filenames, which means you can't keep track of different versions by date. It also means you HAVE to rename the files after you created them with ``mksurfdata.pl``. Now, since ``user_nl`` files are supported for ALL model components, and the same domain files are read by both CLM and DATM and set using the envxml variables: ``ATM_DOMAIN_PATH``, ``ATM_DOMAIN_FILE``, ``LND_DOMAIN_PATH``, and ``LND_DOMAIN_FILE`` -- you can use this mechanism (``user_nl_clm`` and ``user_nl_datm`` and those envxml variables) to point to your datasets in any location. In the future we will deprecate ``CLM_USRDAT_NAME`` and recommend ``user_nl_clm`` and ``user_nl_datm`` and the ``DOMAIN`` envxml variables.
