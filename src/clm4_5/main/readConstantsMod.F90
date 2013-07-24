module readConstantsMod

  !-----------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Read constants
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varctl   , only: use_cn, use_century_decomp, use_nitrif_denitrif
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readConstants
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readConstants ()
    !
    ! !ARGUMENTS:
    implicit none
    !-----------------------------------------------------------------------

    if (use_cn) then
       call CNConstReadFile()
    end if

  end subroutine readConstants

  !-----------------------------------------------------------------------
  subroutine CNConstReadFile ()
    !
    ! !DESCRIPTION:
    ! read CN shared constants
    !
    ! !USES:
    use CNAllocationMod         , only : readCNAllocConsts
    use CNDecompMod             , only : readCNDecompConsts
    use CNDecompCascadeBGCMod   , only : readCNDecompBgcConsts
    use CNDecompCascadeCNMod    , only : readCNDecompCnConsts
    use CNPhenologyMod          , only : readCNPhenolConsts
    use CNMRespMod              , only : readCNMRespConsts
    use CNNDynamicsMod          , only : readCNNDynamicsConsts
    use CNGapMortalityMod       , only : readCNGapMortConsts   
    use CNNitrifDenitrifMod     , only : readCNNitrifDenitrifConsts
    use CNSoilLittVertTranspMod , only : readCNSoilLittVertTranspConsts
    use CNSharedConstsMod       , only : CNConstReadShared,CNConstShareInst
    use clm_varctl              , only : fconsts, iulog
    use spmdMod                 , only : masterproc
    use fileutils               , only : getfil
    use abortutils              , only : endrun
    use ncdio_pio               , only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, &
                                         file_desc_t, ncd_inqdid, ncd_inqdlen
    !
    ! !ARGUMENTS:
    implicit none
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
    character(len=10)  :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    
    if (use_cn) then

       if (masterproc) then
          write(iulog,*) 'readConstantsMod.F90::CNConstReadFile:: reading const file'
       end if
       call getfil (fconsts, locfn, 0)
       call ncd_pio_openfile (ncid, trim(locfn), 0)
       call ncd_inqdid(ncid,'pft',dimid) 
       call ncd_inqdlen(ncid,dimid,npft) 

       call CNConstReadShared(ncid)
       !
       ! populate each module with private constants
       !
       call readCNAllocConsts(ncid)
       call readCNDecompConsts(ncid)
       if (use_century_decomp) then
          call readCNDecompBgcConsts(ncid)
       else
          call readCNDecompCnConsts(ncid)
       end if
       call readCNPhenolConsts(ncid)
       call readCNMRespConsts (ncid)
       call readCNNDynamicsConsts (ncid)
       call readCNGapMortConsts (ncid)
       if (use_nitrif_denitrif) then
          call readCNNitrifDenitrifConsts(ncid)
       end if
       call readCNSoilLittVertTranspConsts(ncid)
       !
       ! close CN const file
       !
       call ncd_pio_closefile(ncid)

    end if

  end subroutine CNConstReadFile

end module readConstantsMod
