.. include:: ../substitutions.rst

.. _creating-surface-datasets:

===========================
 Creating Surface Datasets
===========================

When just creating a replacement file for an existing one, the relevant tool should be used directly to create the file. When you are creating a set of files for a new resolution there are some dependencies between the tools that you need to keep in mind when creating them. The main dependency is that you MUST create a SCRIP grid file first as the SCRIP grid dataset is then input into the other tools. Also look at Table :numref:`reqd-files-table` which gives information on the files required and when. :numref:`Figure Data_Flow` shows an overview of the general data-flow for creation of the fsurdat datasets.

.. _Figure Data_Flow:

.. figure:: mkmapdata_mksurfdata.jpeg

  Data Flow for Creation of Surface Datasets from Raw SCRIP Grid Files

Starting from a SCRIP grid file that describes the grid you will run the model on, you first run ```mkmapdata.sh`` to create a list of mapping files. See :numref:`Figure mkmapdata.sh` for a more detailed view of how ``mkmapdata.sh`` works. The mapping files tell ``mksurfdata_esmf`` how to map between the output grid and the raw datasets that it uses as input. The output of ``mksurfdata_esmf`` is a surface dataset that you then use for running the model. See :numref:`Figure Workflow of CLM5 Land Use Data Tool and mksurfdata_esmf Tool` for a more detailed view of how ``mksurfdata_esmf`` works.

:numref:`Figure Data_Flow_Legend` is the legend for this figure (:numref:`Figure Data_Flow`) and other figures in this chapter (:numref:`Figure Global-Domain` and :numref:`Figure mknoocnmap.pl`).

.. _Figure Data_Flow_Legend:

.. figure:: LegendCLMToolDataFlow.jpeg

  Legend for Data Flow Figures

Green arrows define the input to a program, while red arrows define the output. Cylinders define files that are either created by a program or used as input for a program. Boxes are programs.

You start with a description of a SCRIP grid file for your output grid file and then create mapping files from the raw datasets to it. Once, the mapping files are created ``mksurfdata_esmf`` is run to create the surface dataset to run the model.

Creating a Complete Set of Files for Input to CLM
-------------------------------------------------

1. Create SCRIP grid datasets (if NOT already done)

   First you need to create a descriptor file for your grid, that includes the locations of cell centers and cell corners. There is also a "mask" field, but in this case the mask is set to one everywhere (i.e. all of the masks for the output model grid are "nomask"). An example SCRIP grid file is: ``$CSMDATA/lnd/clm2/mappingdata/grids/SCRIPgrid_10x15_nomask_c110308.nc``. The ``mkmapgrids`` and ``mkscripgrid.ncl`` NCL script in the ``$CTSMROOT/tools/mkmapgrids`` directory can help you with this. SCRIP grid files for all the standard CLM grids are already created for you. See the Section called Creating an output SCRIP grid file at a resolution to run the model on for more information on this.

.. todo::
    Update the below, as domain files aren't needed with nuopc.

2. Create domain dataset (if NOT already done)

   Next use ``gen_domain`` to create a domain file for use by DATM and CLM. This is required, unless a domain file was already created. See the Section called Creating a domain file for CLM and DATM for more information on this.

3. Create mapping files for ``mksurfdata_esmf`` (if NOT already done)

   Create mapping files for ``mksurfdata_esmf`` with ``mkmapdata.sh`` in ``$CTSMROOT/tools/mkmapdata``. See the Section called Creating mapping files that ``mksurfdata_esmf`` will use for more information on this.

4. Create surface datasets

   Next use ``mksurfdata_esmf`` to create a surface dataset, using the mapping datasets created on the previous step as input. There is a version for either clm4_0 or |version| for this program. See the Section called Using ``mksurfdata_esmf`` to create surface datasets from grid datasets for more information on this.

5. Enter the new datasets into the ``build-namelist`` XML database
   The last optional thing to do is to enter the new datasets into the ``build-namelist`` XML database. See Chapter 3 for more information on doing this. This is optional because the user may enter these files into their namelists manually. The advantage of entering them into the database is so that they automatically come up when you create new cases.

The ``$CTSMROOT/tools/README`` goes through the complete process for creating input files needed to run CLM. We repeat that file here:

.. include:: ../../../../tools/README
   :literal:

