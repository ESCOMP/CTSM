.. _managing-your-data-files:

.. include:: ../substitutions.rst

==============================
 Managing Your Data Own Files
==============================

If you are running on a supported machine (such as cheyenne or hopper) the standard input datasets will already be available and you won't have to check them out of the subversion inputdata server. However, you also will NOT be able to add your own datafiles to these standard inputdata directories -- because most likely you won't have permissions to do so. In order to add files to the XML database or to use ``CLM_USRDAT_NAME`` you need to put data in the standard locations so that they can be found. The recommended way to do this is to use the **link_dirtree** tool in the CESM scripts. Some information on **link_dirtree** is available in the `|cesmrelease| Scripts User's Guide <CLM-URL>`_. We also have some examples of it's use here and in other sections of this User's Guide.

Using **link_dirtree** is quite simple, you give the directory where data exists and then the directory that you want to create where datasets will point to the original source files. In the example below we use "$HOME/inputdata", but MYCSMDATA could be any directory you have access to where you want to put your data.

```
> cd scripts
# First make sure you have a inputdata location that you can write to 
# You only need to do this step once, so you won't need to do this in the future
# (except to bring in any updated files in the original $CSMDATA location).
> setenv MYCSMDATA $HOME/inputdata    # Set env var for the directory for input data
> ./link_dirtree $CSMDATA $MYCSMDATA
```

Then when you create a case you will change ``DIN_LOC_ROOT_CSMDATA`` to point to the location you linked to rather than the default location.

``> ./xmlchange DIN_LOC_ROOT_CSMDATA=$MYCSMDATA``

In order to list the files that you have created you merely need to use the UNIX command find to find the files that are NOT softlinks. So for example executing the following command:

``> find $MYCSMDATA -type f -print``

for me gives the following truncated list of CLM_USRDAT_NAME files that I have created.

```
/glade/p/work/erik/inputdata/lnd/clm2/pftdata/pft-physiology.c130503.nc
/glade/p/work/erik/inputdata/atm/datm7/CLM1PT_data/1x1pt_BE-Vie/1997-01.nc
/glade/p/work/erik/inputdata/atm/datm7/CLM1PT_data/1x1pt_BE-Vie/1997-02.nc
/glade/p/work/erik/inputdata/atm/datm7/CLM1PT_data/1x1pt_BE-Vie/1997-03.nc
/glade/p/work/erik/inputdata/atm/datm7/CLM1PT_data/1x1pt_BE-Vie/1997-04.nc
```

You can also use **find** to list files that have a particular pattern in the name as well (using the -name option with wildcards). Also you can always rerun the **link_dirtree** command if any new files are added that you need to be linked into your directory tree. Since, the files are soft-links -- it doesn't take up much space other than the files that you add there. This way all of the files are kept in one place, they are organized by usage according to CESM standards, and you can easily find your own files, and CLM can find them as well.



If you are running on a supported machine (such as cheyenne or hopper) the standard input datasets will already be available and you won't have to check them out of the subversion inputdata server. However, you also will NOT be able to add your own datafiles to these standard inputdata directories -- because most likely you won't have permissions to do so. In order to add files to the XML database or to use ``CLM_USRDAT_NAME`` you need to put data in the standard locations so that they can be found. The recommended way to do this is to use the **link_dirtree** tool in the CESM scripts. Some information on **link_dirtree** is available in the `|cesmrelease| Scripts User's Guide <CLM-URL>`_. We also have some examples of it's use here and in other sections of this User's Guide.

Using **link_dirtree** is quite simple, you give the directory where data exists and then the directory that you want to create where datasets will point to the original source files. In the example below we use "$HOME/inputdata", but MYCSMDATA could be any directory you have access to where you want to put your data.

```
> cd scripts
# First make sure you have a inputdata location that you can write to 
# You only need to do this step once, so you won't need to do this in the future
# (except to bring in any updated files in the original $CSMDATA location).
> setenv MYCSMDATA $HOME/inputdata    # Set env var for the directory for input data
> ./link_dirtree $CSMDATA $MYCSMDATA
```

Then when you create a case you will change ``DIN_LOC_ROOT_CSMDATA`` to point to the location you linked to rather than the default location.

``> ./xmlchange DIN_LOC_ROOT_CSMDATA=$MYCSMDATA``

In order to list the files that you have created you merely need to use the UNIX command find to find the files that are NOT softlinks. So for example executing the following command:

``> find $MYCSMDATA -type f -print``

for me gives the following truncated list of CLM_USRDAT_NAME files that I have created.

```
/glade/p/work/erik/inputdata/lnd/clm2/pftdata/pft-physiology.c130503.nc
/glade/p/work/erik/inputdata/atm/datm7/CLM1PT_data/1x1pt_BE-Vie/1997-01.nc
/glade/p/work/erik/inputdata/atm/datm7/CLM1PT_data/1x1pt_BE-Vie/1997-02.nc
/glade/p/work/erik/inputdata/atm/datm7/CLM1PT_data/1x1pt_BE-Vie/1997-03.nc
/glade/p/work/erik/inputdata/atm/datm7/CLM1PT_data/1x1pt_BE-Vie/1997-04.nc
```

You can also use **find** to list files that have a particular pattern in the name as well (using the -name option with wildcards). Also you can always rerun the **link_dirtree** command if any new files are added that you need to be linked into your directory tree. Since, the files are soft-links -- it doesn't take up much space other than the files that you add there. This way all of the files are kept in one place, they are organized by usage according to CESM standards, and you can easily find your own files, and CLM can find them as well.
