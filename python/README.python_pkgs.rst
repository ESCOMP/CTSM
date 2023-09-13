.. sectnum::

.. contents::

=====================================================================================
 Requirements to consider for python packages to be included for CTSM python tools
=====================================================================================

Requirements for CTSM python tools:

Any added dependencies should be discussed and approved by the CTSM-software team. Criteria for evaluation of third-party packages
include:
- How much value does this package provide to us beyond what we could get without it?
- How difficult it is to install this package (including: how much does this complicate / slow down the creation of the conda
environment)?
- How stable is this package?
- How well maintained is this package?
- Are there other packages that are more stable or better maintained that would provide nearly the same level of value?
- Tools that require extra packages should be done in a "contrib" type area out of the main part of tools (this would apply for
advanced plotting capability for example)
- We need to be able to reproduce working conda environments minimally on our test machines (currently cheyenne and izumi), but also
on any machine that we run CTSM on. If there is a machine that we can run CTSM on that we can't build the conda environments or run
the tools on -- that needs to be fixed.
- We need to tell the user how long to expect the conda environment to load, and give them options if the conda load is taking too
long
- Conda environments need to build robustly even for users who don't have ctsm_pylib loaded in their conda environment
- Currently dask will NOT be something we require for any of the main CTSM tools
- Currently we won't use conda-lock
- We specify the black version exactly so that black will function identically for all users
- We specify the pylint version exactly because pylint is finicky with version and we need it to work identically for all developers
- We might remove the need for python packages that aren't providing enough utility
Packages where specific versions seem to be required will have the version requirements in a least the >= form if not an exact
version
