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
  logical            :: outnc_large_files     ! output files in 64-bit format for large files
  logical            :: outnc_double          ! output ALL data in files as 64-bit
  logical            :: all_urban             ! output ALL data as 100% covered in urban
  character(len=32)  :: mksrf_gridnm   = ' '  ! name of grid to use on output file
  character(len=256) :: mksrf_fgrid    = ' '  ! land grid file name to use 
  character(len=256) :: mksrf_gridtype = ' '  ! land gridtype, global or reg
  character(len=256) :: mksrf_fvegtyp  = ' '  ! vegetation data file name
  character(len=256) :: mksrf_fsoitex  = ' '  ! soil texture data file name
  character(len=256) :: mksrf_forganic = ' '  ! organic matter data file name
  character(len=256) :: mksrf_fsoicol  = ' '  ! soil color data file name
  character(len=256) :: mksrf_flanwat  = ' '  ! inland water data file name
  character(len=256) :: mksrf_furban   = ' '  ! urban data file name
  character(len=256) :: mksrf_firrig   = ' '  ! irrigated area data file name
  character(len=256) :: mksrf_fglacier = ' '  ! glacier data file name
  character(len=256) :: mksrf_ftopo    = ' '  ! topography data file name
  character(len=256) :: mksrf_ffrac    = ' '  ! grid fraction data file name
  character(len=256) :: mksrf_fmax     = ' '  ! fmax data file name
  character(len=256) :: mksrf_flai     = ' '  ! lai data filename
  character(len=256) :: mksrf_fdynuse  = ' '  ! ascii file containing names of dynamic land use files
  character(len=256) :: mksrf_fvocef   = ' '  ! VOC Emission Factor data file name
  integer, public    :: nglcec         = 10   ! number of elevation classes for glaciers
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 11/04
!
!EOP
!-----------------------------------------------------------------------

end module mkvarctl
