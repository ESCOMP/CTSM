module mkvarctl

  !-----------------------------------------------------------------------
  ! Module containing control variables
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  integer, public  :: ndiag             ! output log unit
  logical, public  :: root_task         ! proc 0 logical for printing msgs
  integer, public  :: iam               ! processor number
  integer, public  :: npes              ! number of processors
  integer, public  :: mpicom            ! communicator group

  logical, public  :: outnc_large_files ! output files in 64-bit format for large files
  logical, public  :: outnc_double      ! output ALL data in files as 64-bit
  integer, public  :: outnc_dims = 2    ! only applicable to lat/lon grids
  logical, public  :: outnc_1d          ! true => output file is 1d
  logical, public  :: outnc_vic         ! true => output VIC fields
  logical, public  :: outnc_3dglc       ! true => output 3D glacier fields

  integer, public  :: numpft  = 16      ! number of plant types
  integer, public  :: nglcec = 10       ! number of elevation classes for glaciers

  integer , public :: soil_color_override
  real(r8), public :: soil_sand_override
  real(r8), public :: soil_clay_override
  real(r8), public :: soil_fmax_override

  ! Variables to override data read in with
  ! (all_urban is mostly for single-point mode, but could be used for sensitivity studies)
  logical  , public :: all_urban             ! output ALL data as 100% covered in urban
  logical  , public :: all_veg               ! if gridcell will be 100% vegetated land-cover
  logical  , public :: no_inlandwet          ! set wetland to 0% over land; wetland will only be used for ocean points
  real(r8) , public :: std_elev = -999.99_r8 ! Standard deviation of elevation (m) to use for ent

  real(r8), public, parameter :: spval     = 1.e36 ! special value
  integer,  public, parameter :: ispval    = -9999 ! special value
  integer , public, parameter :: unsetcol  = -999  ! flag to indicate soil color NOT set
  real(r8), public, parameter :: unsetsoil = -999.99_r8 ! Flag to signify soil texture override not set

end module mkvarctl
