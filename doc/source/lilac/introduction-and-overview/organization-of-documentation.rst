.. _organization-of-documentation:

.. highlight:: shell

===================================
 Organization of the documentation
===================================

This documentation is organized into the following high-level sections:

- :numref:`{number}. {name} <introduction-and-overview>`: This section gives a general
  introduction to LILAC and describes the organization of this documentation (you're
  reading this now!)

- :numref:`{number}. {name} <obtaining-building-and-running>`: This section provides
  instructions for building and running CTSM within an atmosphere model that has been set
  up to run with CTSM (e.g., WRF). If you are starting to use CTSM with an atmosphere
  model that does not yet have any calls to CTSM-LILAC, then you should start with section
  :numref:`{number}. {name} <calling-ctsm-from-atm>`.

- :numref:`{number}. {name} <calling-ctsm-from-atm>`: This section provides details on the
  Fortran code that needs to be added to an atmosphere model in order to call CTSM via
  LILAC as its land surface scheme. (In practice, this step comes before
  :numref:`{number}. {name} <obtaining-building-and-running>`, but it is included later in
  the documentation because it is of interest to fewer people.)

- :numref:`{number}. {name} <specific-atm-models>`: This section provides notes on running
  CTSM within specific atmosphere models.
