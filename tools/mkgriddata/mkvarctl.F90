module mkvarctl

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkvarctl
!
! !DESCRIPTION:
! Module containing control variables
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
  private
!
  character(len=256), public :: mksrf_fnavyoro = ' ' ! filepath for 20 min navy orography dataset
  character(len=256), public :: mksrf_frawtopo = ' ' ! filepath for raw topo file
  character(len=256), public :: mksrf_fcamfile = ' ' ! filepath for cam grid file
  character(len=256), public :: mksrf_fcamtopo = ' ' ! filepath for cam topo file
  character(len=256), public :: mksrf_fccsmdom = ' ' ! filepath for cesm domain file
  character(len=256), public :: mksrf_fclmgrid = ' ' ! filepath for clm grid file
  integer           , public :: mksrf_lsmlon = 0     ! create create with this number of longitudes
  integer           , public :: mksrf_lsmlat = 0     ! create create with this number of latitudes
  real(r8)          , public :: mksrf_edgen          ! northern edge of grid (degrees): >  -90 and < 90
  real(r8)          , public :: mksrf_edgee          ! eastern edge of grid (degrees) : see following notes
  real(r8)          , public :: mksrf_edges          ! southern edge of grid (degrees): >= -90 and <  90
  real(r8)          , public :: mksrf_edgew          ! western edge of grid (degrees) : see following notes

  integer           , public :: area_units = 0       ! 0=km2, 1=rad2
  logical           , public :: area_valid = .true.  ! is area read valid
  integer           , public :: numpft     = 16      ! number of plant types

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 11/04
!
!EOP
!-----------------------------------------------------------------------

end module mkvarctl
