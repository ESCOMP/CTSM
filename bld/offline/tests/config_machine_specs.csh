#!/bin/csh
#-----------------------------------------------------------------------
#
# Set the variables needed for the configuration process that are specific
# to a given platform. The settings here are specific to NCAR platforms.
# The perl objects are able to provide default settings for several
# different sites as well as different platforms at those sites.
#
# This script sets the following ENV variables for the given platform.
#
#   LIB_NETCDF 	Location of directory with the netCDF library.
#   INC_NETCDF 	Location of directory with the netCDF includes.
#   LIB_MPI    	Location of directory with the MPI library.
#   INC_MPI    	Location of directory with the MPI includes.
#   SPMD       	If MPI distributed memory option should be turned on.
#   
# If these values are set outside this script those values will be retained.
#
# Run-time ENV variables are set in the CAM_run.pm Perl module.
#
#-----------------------------------------------------------------------
setenv OP_SYSTEM     `uname -s`

if ( ! $?ROOTDIR )then
  setenv ROOTDIR `pwd`/../src
endif

if ( $OP_SYSTEM == SunOS )then
#-----------------------------------------------------------------------
# SunOS (CGD Suns)
#-----------------------------------------------------------------------

  if ( ! $?LIB_NETCDF ) setenv LIB_NETCDF  /contrib/netcdf/lib
  if ( ! $?INC_NETCDF ) setenv INC_NETCDF  /contrib/include

  if ( ! $?LIB_MPI ) setenv LIB_MPI     /contrib/lib
  if ( ! $?INC_MPI ) setenv INC_MPI     /contrib/include

  if ( ! $?SPMD ) setenv SPMD TRUE

else if ( $OP_SYSTEM == AIX )then
#-----------------------------------------------------------------------
# AIX (bluesky)
#-----------------------------------------------------------------------

  if ( ! $?LIB_NETCDF ) setenv LIB_NETCDF  /usr/local/lib64/r4i4
  if ( ! $?INC_NETCDF ) setenv INC_NETCDF  /usr/local/include

  if ( ! $?SPMD ) setenv SPMD TRUE

else if ( $OP_SYSTEM == IRIX64 )then
#-----------------------------------------------------------------------
# SGI (tempest)
#-----------------------------------------------------------------------

  if ( ! $?LIB_NETCDF ) setenv LIB_NETCDF  /usr/local/lib64/r4i4
  if ( ! $?INC_NETCDF ) setenv INC_NETCDF  /usr/local/include

  if ( ! $?SPMD ) setenv SPMD TRUE

else if ( $OP_SYSTEM == Darwin )then
#-----------------------------------------------------------------------
# Darwin (MacOS)
#-----------------------------------------------------------------------

  if ( ! $?LIB_NETCDF ) setenv LIB_NETCDF  /usr/local/lib
  if ( ! $?INC_NETCDF ) setenv INC_NETCDF  /usr/local/include

  if ( ! $?SPMD ) setenv SPMD FALSE

else if ( $OP_SYSTEM == Linux )then
#-----------------------------------------------------------------------
# Linux (bangkok)
#-----------------------------------------------------------------------

  setenv MPI_ROOT "/usr/local/mpich-1.2.7p1-pgi-pgcc-pghf-6.1-3"
  if ( $?USER_FC ) setenv MPI_ROOT "/usr/local/mpich-1.2.7p1-gcc-g++-4.0.2-8-lf9562"
  if ( ! $?LIB_MPI ) setenv LIB_MPI $MPI_ROOT/lib
  if ( ! $?INC_MPI ) setenv INC_MPI $MPI_ROOT/include

  setenv NETCDF_ROOT "/usr/local/netcdf-3.6.1-beta3-pgi-hpf-cc-6.0-5"
  if ( $?USER_FC ) setenv NETCDF_ROOT "/usr/local/netcdf-3.6.1beta3-gcc-4.0.2-g77-lf9562"
  if ( ! $?LIB_NETCDF )  setenv LIB_NETCDF    $NETCDF_ROOT/lib
  if ( ! $?INC_NETCDF )  setenv INC_NETCDF    $NETCDF_ROOT/include

  if ( ! $?SPMD ) setenv SPMD    TRUE

else if ( $OP_SYSTEM == OSF1 )then
#-----------------------------------------------------------------------
# Compaq (prospect)
#-----------------------------------------------------------------------

  if ( ! $?LIB_NETCDF ) setenv LIB_NETCDF  /usr/local/lib32/r4i4
  if ( ! $?INC_NETCDF ) setenv INC_NETCDF  /usr/local/include

  if ( ! $?SPMD ) setenv SPMD TRUE
else if ( $OP_SYSTEM == UNICOS/mp )then
#-----------------------------------------------------------------------
# Cray X1 (phoenix)
#-----------------------------------------------------------------------

  if ( ! $?LIB_NETCDF ) setenv LIB_NETCDF  /apps/netcdf/prod/x1_r8/lib
  if ( ! $?INC_NETCDF ) setenv INC_NETCDF  /apps/netcdf/prod/x1_r8/include

  if ( ! $?SPMD ) setenv SPMD TRUE
else
  echo "Do not know how to handle this architecture: $OP_SYSTEM"
  exit
endif
