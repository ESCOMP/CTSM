module spmdMod
  ! Stub of spmdMod

  implicit none
  save
  private

  logical, parameter, public :: masterproc = .true.
  integer, parameter, public :: iam = 0
  integer, parameter, public :: mpicom = 0
  integer, parameter, public :: mpi_integer = 0
  integer, parameter, public :: mpi_logical = 0
  integer, parameter, public :: MPI_REAL8 = 1275070505

end module spmdMod
