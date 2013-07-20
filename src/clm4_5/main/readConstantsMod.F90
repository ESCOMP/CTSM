module readConstantsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: readConstantsMod
!
! !DESCRIPTION:
! 
! 
!
! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none
  save
  private
! !PUBLIC MEMBER FUNCTIONS:
  public :: readConstants
!
! !REVISION HISTORY:
!   Dec 3 2012 : Created by S. Muszala
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readConstants
!
! !INTERFACE:
subroutine readConstants ()
!
! !DESCRIPTION:
! 
! 
!
! !USES:
   use clm_varctl, only : defCN
!
! !ARGUMENTS:
   implicit none
!
! !CALLED FROM:clm_initializeMod.F90::initialize1 
! 
!
! !REVISION HISTORY:
!   Dec 3 2012 : Created by S. Muszala
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
!
! !OTHER LOCAL VARIABLES:

!EOP
!-----------------------------------------------------------------------

#ifdef CN
      call CNConstReadFile()
#endif

end subroutine readConstants

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNConstReadFile
!
! !INTERFACE:
#ifdef CN
subroutine CNConstReadFile ()
!
! !DESCRIPTION:
! 
! 
!
! !USES:
!
   use CNAllocationMod  , only : readCNAllocConsts
   use CNDecompMod             , only : readCNDecompConsts
#ifdef CENTURY_DECOMP
   use CNDecompCascadeBGCMod   , only : readCNDecompBgcConsts
#else
   use CNDecompCascadeCNMod    , only : readCNDecompCnConsts
#endif
   use CNPhenologyMod   , only : readCNPhenolConsts
   use CNMRespMod       , only : readCNMRespConsts
   use CNNDynamicsMod   , only : readCNNDynamicsConsts
   use CNGapMortalityMod, only : readCNGapMortConsts   
#ifdef NITRIF_DENITRIF
   use CNNitrifDenitrifMod     , only : readCNNitrifDenitrifConsts
#endif
   use CNSoilLittVertTranspMod , only: readCNSoilLittVertTranspConsts
   use CNSharedConstsMod, only : CNConstReadShared,CNConstShareInst
   use clm_varctl       , only : fconsts,iulog
   use ncdio_pio        , only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t, &
                                 ncd_inqdid, ncd_inqdlen
   use spmdMod          , only : masterproc
   use fileutils        , only : getfil
   use abortutils       , only: endrun

! !ARGUMENTS:
   implicit none
!
! !CALLED FROM:  readConstantsMod.F90::readConstants
! 
!
! !REVISION HISTORY:
!   Dec 3 2012 : Created by S. Muszala
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
!
! !OTHER LOCAL VARIABLES:
   character(len=32)  :: subname = 'CNConstReadShared'
   character(len=100) :: errCode = 'Error reading in CN shared const file '
   character(len=256) :: locfn ! local file name
   type(file_desc_t)  :: ncid  ! pio netCDF file id
   integer            :: dimid ! netCDF dimension id
   integer            :: npft  ! number of pfts on pft-physiology file
   logical            :: readv ! has variable been read in or not
   real(r8)           :: tempr ! temporary to read in constant
   character(len=10) :: tString ! temp. var for reading

!EOP
!-----------------------------------------------------------------------

   !
   ! read CN shared consts
   !
   if (masterproc) then
      write(iulog,*) 'readConstantsMod.F90::CNConstReadFile:: reading const file'
   end if
   call getfil (fconsts, locfn, 0)
   call ncd_pio_openfile (ncid, trim(locfn), 0)
   call ncd_inqdid(ncid,'pft',dimid) 
   call ncd_inqdlen(ncid,dimid,npft) 

   !
   ! read CN shared consts
   !
   call CNConstReadShared(ncid)

   !
   ! populate each module with private constants
   !
   call readCNAllocConsts(ncid)
   call readCNDecompConsts(ncid)
#ifdef CENTURY_DECOMP
   call readCNDecompBgcConsts(ncid)
#else
   call readCNDecompCnConsts(ncid)
#endif
   call readCNPhenolConsts(ncid)
   call readCNMRespConsts (ncid)
   call readCNNDynamicsConsts (ncid)
   call readCNGapMortConsts (ncid)
#ifdef NITRIF_DENITRIF
   call readCNNitrifDenitrifConsts(ncid)
#endif
   call readCNSoilLittVertTranspConsts(ncid)
   !
   ! close CN const file
   !
   call ncd_pio_closefile(ncid)

end subroutine CNConstReadFile
#endif

end module readConstantsMod
