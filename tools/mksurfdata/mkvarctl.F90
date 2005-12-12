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
  character(len=256) :: mksrf_fgrid_global = ' '   ! land grid file name to use 
  character(len=256) :: mksrf_fgrid_regional = ' ' ! land grid file name to use 
  character(len=256) :: mksrf_fvegtyp              ! vegetation data file name
  character(len=256) :: mksrf_fsoitex              ! soil texture data file name
  character(len=256) :: mksrf_fsoicol              ! soil color data file name
  character(len=256) :: mksrf_flanwat              ! inland water data file name
  character(len=256) :: mksrf_furban               ! urban data file name
  character(len=256) :: mksrf_fglacier             ! glacier data file name
  character(len=256) :: mksrf_flai                 ! lai data filename
  character(len=256) :: mksrf_fdynuse              ! ascii file containing names of dynamic land use files
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 11/04
!
!EOP
!-----------------------------------------------------------------------

end module mkvarctl
