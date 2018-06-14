.. _creating-surface-datasets:

===========================
 Creating Surface Datasets
===========================

When just creating a replacement file for an existing one, the relevant tool should be used directly to create the file. When you are creating a set of files for a new resolution there are some dependencies between the tools that you need to keep in mind when creating them. The main dependency is that you MUST create a SCRIP grid file first as the SCRIP grid dataset is then input into the other tools. Also look at `Table 3-1 <CLM-URL>`_ which gives information on the files required and when. `Figure 2-1 <CLM-URL>`_ shows an overview of the general data-flow for creation of the fsurdat datasets.

Figure 2-1. Data Flow for Creation of Surface Datasets from Raw SCRIP Grid Files
--------------------------------------------------------------------------------
Insert figure 2-1

Starting from a SCRIP grid file that describes the grid you will run the model on, you first run **mkmapdata.sh** to create a list of mapping files. See `Figure 2-3 <CLM-URL>`_ for a more detailed view of how **mkmapdata.sh** works. The mapping files tell **mksurfdata_map** how to map between the output grid and the raw datasets that it uses as input. The output of **mksurfdata_map** is a surface dataset that you then use for running the model. See `Figure 2-6 <CLM-URL>`_ for a more detailed view of how **mksurfdata_map** works.

`Figure 2-2 <CLM-URL>`_ is the legend for this figure (`Figure 2-1 <CLM-URL>`_) and other figures in this chapter (`Figure 2-4 <CLM-URL>`_, `Figure 2-5 <CLM-URL>`_, and `Figure 2-6 <CLM-URL>`_).
Figure 2-2. Legend for Data Flow Figures
Insert figure 2-2

Green arrows define the input to a program, while red arrows define the output. Cylinders define files that are either created by a program or used as input for a program. Boxes are programs.

You start with a description of a SCRIP grid file for your output grid file and then create mapping files from the raw datasets to it. Once, the mapping files are created **mksurfdata_map** is run to create the surface dataset to run the model.

Creating a Complete Set of Files for Input to CLM
-------------------------------------------------

1. Create SCRIP grid datasets (if NOT already done)

   First you need to create a descriptor file for your grid, that includes the locations of cell centers and cell corners. There is also a "mask" field, but in this case the mask is set to one everywhere (i.e. all of the masks for the output model grid are "nomask"). An example SCRIP grid file is: $CSMDATA/lnd/clm2/mappingdata/grids/SCRIPgrid_10x15_nomask_c110308.nc. The mkmapgrids and mkscripgrid.ncl NCL script in the models/lnd/clm/tools/shared/mkmapgrids directory can help you with this. SCRIP grid files for all the standard CLM grids are already created for you. See the Section called Creating an output SCRIP grid file at a resolution to run the model on for more information on this.

2. Create domain dataset (if NOT already done)

   Next use gen_domain to create a domain file for use by DATM and CLM. This is required, unless a domain file was already created. See the Section called Creating a domain file for CLM and DATM for more information on this.

3. Create mapping files for mksurfdata_map (if NOT already done)

   Create mapping files for mksurfdata_map with mkmapdata.sh in models/lnd/clm/tools/shared/mkmapdata. See the Section called Creating mapping files that mksurfdata_map will use for more information on this.

4. Create surface datasets

   Next use mksurfdata_map to create a surface dataset, using the mapping datasets created on the previous step as input. There is a version for either clm4_0 or clm4_5 for this program. See the Section called Using mksurfdata_map to create surface datasets from grid datasets for more information on this.

5. Create some sort of initial condition dataset

   You then need to do one of the following three options to have an initial dataset to start from.

   a. Use spinup-procedures to create initial condition datasets

      The first option is to do the spinup procedures from arbitrary initial conditions to get good initial datasets. This is the most robust method to use. See the Section called Spinning up the Satellite Phenology Model (CLMSP spinup) in Chapter 4, the Section called Spinning up the CLM4.0 biogeochemistry Carbon-Nitrogen Model (CN spinup) in Chapter 4, or the Section called Spinning up the CLM4.0 Carbon-Nitrogen Dynamic Global Vegetation Model (CNDV spinup) in Chapter 4 for more information on this.

   b. Use interpinic to interpolate existing initial condition datasets

      The next option is to interpolate from spunup datasets at a different resolution, using interpinic. There is a version for either clm4_0 or clm4_5 for this program. See the Section called Using interpinic to interpolate initial conditions to different resolutions for more information on this.

   c. Start up from arbitrary initial conditions

      The last alternative is to run from arbitrary initial conditions without using any spun-up datasets. This is inappropriate when using CLM4.5-BGC or CLMCN (bgc=cn or cndv) as it takes a long time to spinup Carbon pools.

.. warning:: This is NOT recommended as many fields in CLM take a long time to equilibrate.

6. Enter the new datasets into the build-namelist XML database
   The last optional thing to do is to enter the new datasets into the build-namelist XML database. See Chapter 3 for more information on doing this. This is optional because the user may enter these files into their namelists manually. The advantage of entering them into the database is so that they automatically come up when you create new cases.

The ``models/lnd/clm/tools/README`` goes through the complete process for creating input files needed to run CLM. We repeat that file here:

.. include:: ../../clm5.0/tools/README
   :literal:

