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
  character(len=256), public :: mksrf_fnavyoro = ' ' ! directory for 20 min navy orography dataset
  character(len=256), public :: mksrf_frawtopo = ' ' ! directory for raw topo file
  character(len=256), public :: mksrf_fcamfile = ' ' ! directory for cam grid file
  character(len=256), public :: mksrf_fcamtopo = ' ' ! directory for cam topo file
  character(len=256), public :: mksrf_fccsmdom = ' ' ! directory for ccsm domain file
  character(len=256), public :: mksrf_fclmgrid = ' ' ! directory for ccsm domain file
  integer           , public :: mksrf_lsmlon = 0     ! create create with this number of longitudes
  integer           , public :: mksrf_lsmlat = 0     ! create create with this number of latitudes
  real(r8)          , public :: mksrf_edgen          ! northern edge of grid (degrees): >  -90 and < 90
  real(r8)          , public :: mksrf_edgee          ! eastern edge of grid (degrees) : see following notes
  real(r8)          , public :: mksrf_edges          ! southern edge of grid (degrees): >= -90 and <  90
  real(r8)          , public :: mksrf_edgew          ! western edge of grid (degrees) : see following notes

  integer           , public :: area_units = 0       ! 0=km2, 1=rad2
  logical           , public :: area_valid = .true.  ! is area read valid

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 11/04
!
!EOP
!-----------------------------------------------------------------------

end module mkvarctl
