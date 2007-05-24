module subgridMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: subgridMod
!
! !DESCRIPTION:
! sub-grid data and mapping types and modules
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use nanMod      , only : bigint,nan
  use abortutils  , only : endrun

  implicit none
  private	
  save

! !PUBLIC MEMBER FUNCTIONS:
  public subgrid_get_gcellinfo        ! Returns g,l,c,p properties from wtxy


! !REVISION HISTORY:
! 2006.07.04 T Craig, rename initSubgridMod
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !LOCAL MODULE VARIABLES:
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgrid_get_gcellinfo
!
! !INTERFACE:
  subroutine subgrid_get_gcellinfo (nw, &
                             nlunits, ncols, npfts, &
                             nveg, wtveg, &
                             ncrop, wtcrop, &
                             nlake, wtlake, &
                             nwetland, wtwetland, &
                             nglacier, wtglacier)
!
! !DESCRIPTION:
! Obtain gridcell properties
!
! !USES
  use clm_varpar  , only : numpft, maxpatch_pft, &
                           npatch_lake, npatch_glacier, npatch_wet, npatch_crop
  use clm_varctl  , only : allocate_all_vegpfts
  use clm_varsur  , only : wtxy

! !ARGUMENTS
    implicit none
    integer , intent(in)  :: nw                   ! wtxy cell index
    integer , optional, intent(out) :: nlunits    ! number of landunits
    integer , optional, intent(out) :: ncols      ! number of columns 
    integer , optional, intent(out) :: npfts      ! number of pfts 
    integer , optional, intent(out) :: nveg       ! number of vegetated pfts in naturally vegetated landunit
    real(r8), optional, intent(out) :: wtveg      ! weight (relative to gridcell) of naturally vegetated landunit
    integer , optional, intent(out) :: ncrop      ! number of crop pfts in crop landunit
    real(r8), optional, intent(out) :: wtcrop     ! weight (relative to gridcell) of crop landunit
    integer , optional, intent(out) :: nlake      ! number of lake pfts (columns) in lake landunit
    real(r8), optional, intent(out) :: wtlake     ! weight (relative to gridcell) of lake landunitof lake pfts (columns) in lake landunit
    integer , optional, intent(out) :: nwetland   ! number of wetland pfts (columns) in wetland landunit
    real(r8), optional, intent(out) :: wtwetland  ! weight (relative to gridcell) of wetland landunitof wetland pfts (columns) in wetland landunit
    integer , optional, intent(out) :: nglacier   ! number of glacier pfts (columns) in glacier landunit
    real(r8), optional, intent(out) :: wtglacier  ! weight (relative to gridcell) of glacier landunitof glacier pfts (columns) in glacier landunit
!
! !CALLED FROM:
! subroutines decomp_init, initGridCells
!
! !REVISION HISTORY:
! 2002.09.11  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                ! loop index
    integer  :: ipfts            ! number of pfts in gridcell
    integer  :: icols            ! number of columns in gridcell
    integer  :: ilunits          ! number of landunits in gridcell
    integer  :: npfts_per_lunit  ! number of pfts in landunit
    real(r8) :: wtlunit          ! weight (relative to gridcell) of landunit
!------------------------------------------------------------------------------

    ! Initialize pfts, columns and landunits counters for gridcell

    ipfts   = 0
    icols   = 0
    ilunits = 0

    ! Set naturally vegetated landunit assuming that 
    ! the vegetated landunit has one column

    npfts_per_lunit = 0
    wtlunit = 0._r8
    do m = 1, maxpatch_pft            
       if (wtxy(nw,m) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1
          wtlunit = wtlunit + wtxy(nw,m)
       end if
    end do
    if (npfts_per_lunit > 0) then
       if (allocate_all_vegpfts) npfts_per_lunit = numpft+1
       ilunits = ilunits + 1
       icols = icols + 1  
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nveg )) nveg  = npfts_per_lunit
    if (present(wtveg)) wtveg = wtlunit

    ! Set lake landunit

    npfts_per_lunit = 0
    wtlunit = 0._r8
    if (wtxy(nw,npatch_lake) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(nw,npatch_lake)
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nlake )) nlake  = npfts_per_lunit
    if (present(wtlake)) wtlake = wtlunit

    ! Set wetland landunit

    npfts_per_lunit = 0
    wtlunit = 0._r8
    if (wtxy(nw,npatch_wet) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(nw,npatch_wet)
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nwetland )) nwetland  = npfts_per_lunit
    if (present(wtwetland)) wtwetland = wtlunit

    ! Set glacier landunit

    npfts_per_lunit = 0
    wtlunit = 0._r8
    if (wtxy(nw,npatch_glacier) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(nw,npatch_glacier)
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nglacier )) nglacier  = npfts_per_lunit
    if (present(wtglacier)) wtglacier = wtlunit

    ! Set crop landunit if appropriate

    npfts_per_lunit = 0
    wtlunit = 0._r8
    do m = npatch_glacier+1, npatch_crop 
       if (wtxy(nw,m) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1
          wtlunit = wtlunit + wtxy(nw,m)
       end if
    end do
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(ncrop )) ncrop  = npfts_per_lunit
    if (present(wtcrop)) wtcrop = wtlunit

    ! Determine return arguments

    if (present(nlunits)) nlunits = ilunits
    if (present(ncols))   ncols   = icols
    if (present(npfts))   npfts   = ipfts

  end subroutine subgrid_get_gcellinfo

!-----------------------------------------------------------------------

end module subgridMod
