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
!
  character(len=256) :: mksrf_fnavyoro = ' ' ! directory for 20 min navy orography dataset
  character(len=256) :: mksrf_frawtopo = ' ' ! directory for raw topo file
  character(len=256) :: mksrf_fcamfile = ' ' ! directory for cam grid file
  character(len=256) :: mksrf_fcamtopo = ' ' ! directory for cam topo file
  character(len=256) :: mksrf_fccsmdom = ' ' ! directory for ccsm domain file
  character(len=256) :: mksrf_fclmgrid = ' ' ! directory for ccsm domain file
  integer            :: mksrf_lsmlon = 0     ! create create with this number of longitudes
  integer            :: mksrf_lsmlat = 0     ! create create with this number of latitudes
  real(r8)           :: mksrf_edgen          ! northern edge of grid (degrees): >  -90 and < 90
  real(r8)           :: mksrf_edgee          ! eastern edge of grid (degrees) : see following notes
  real(r8)           :: mksrf_edges          ! southern edge of grid (degrees): >= -90 and <  90
  real(r8)           :: mksrf_edgew          ! western edge of grid (degrees) : see following notes

  integer            :: area_units = 0       ! 0=km2, 1=rad2
  logical            :: area_valid = .true.  ! is area read valid

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 11/04
!
!EOP
!-----------------------------------------------------------------------

end module mkvarctl
