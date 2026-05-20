(comparing-history-files)=

# Comparing History Files

## `nccmp`

`nccmp` allows you to compare data and/or metadata between two netCDF files. It is available as a module (`module load nccmp`) on Derecho and Casper. See the project's [GitLab](https://gitlab.com/remikz/nccmp/-/blob/master/README.md) for more information and installation instructions.

## `cprnc`

`cprnc` ([GitHub](https://github.com/ESMCI/cprnc)) is a tool shared across CESM to compare two NetCDF history files. It differences every field that is shared on both files, and reports a summary of the difference. The summary includes the three largest differences, as well as the root mean square (RMS) difference. It also gives some summary information on the field as well. You have to enter at least one file, and up to two files. With one file it gives you summary information on the file, and with two it gives you information on the differences between the two. At the end it will give you a summary of the fields compared and how many fields were different and how many were identical.

`cprnc` is available as a module (`module load cprnc`) on Derecho and Casper. To install on a different machine, from the top level of a CTSM checkout, do:
```shell
cd cime/CIME/non_py/cprnc
mkdir bld
cd bld
cmake ../
make
```

The `cprnc` README is embedded below for your convenience. Note that this is for the version included with CTSM (via CIME), which might differ from the version (if any) available on your machine.

```{include} ../../../../cime/CIME/non_py/cprnc/README.md
```
