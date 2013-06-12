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
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog

  implicit none
  private	
  save

! !PUBLIC MEMBER FUNCTIONS:
  public subgrid_get_gcellinfo        ! Obtain gridcell properties


! !REVISION HISTORY:
! 2006.07.04 T Craig, rename initSubgridMod
!
!
! !PRIVATE MEMBER FUNCTIONS: None
!
! !PRIVATE DATA MEMBERS: None
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgrid_get_gcellinfo
!
! !INTERFACE:
  subroutine subgrid_get_gcellinfo (gi, &
                             nlunits, ncols, npfts, &
                             nveg, &
                             ncrop, &
                             nurban_tbd, &
                             nurban_hd, &
                             nurban_md, &
                             nlake, &
                             nwetland, &
                             nglacier, &
                             nglacier_mec,  &
                             glcmask)
!
! !DESCRIPTION:
! Obtain gridcell properties
!
! !USES
  use clm_varpar  , only : natpft_size, cft_size, maxpatch_urb, maxpatch_glcmec
  use clm_varctl  , only : allocate_all_vegpfts, create_crop_landunit
  use clm_varsur  , only : wt_lunit, wt_glc_mec
  use clm_varcon  , only : istsoil, istcrop, istice, istice_mec, istdlak, istwet, &
                           isturb_tbd, isturb_hd, isturb_md

! !ARGUMENTS
    implicit none
    integer , intent(in)  :: gi                   ! grid cell index
    integer , optional, intent(out) :: nlunits    ! number of landunits
    integer , optional, intent(out) :: ncols      ! number of columns 
    integer , optional, intent(out) :: npfts      ! number of pfts 
    integer , optional, intent(out) :: nveg       ! number of vegetated pfts in naturally vegetated landunit
    integer , optional, intent(out) :: ncrop      ! number of crop pfts in crop landunit
    integer , optional, intent(out) :: nurban_tbd ! number of urban pfts (columns) in urban TBD landunit
    integer , optional, intent(out) :: nurban_hd ! number of urban pfts (columns) in urban HD landunit
    integer , optional, intent(out) :: nurban_md ! number of urban pfts (columns) in urban MD landunit
    integer , optional, intent(out) :: nlake      ! number of lake pfts (columns) in lake landunit
    integer , optional, intent(out) :: nwetland   ! number of wetland pfts (columns) in wetland landunit
    integer , optional, intent(out) :: nglacier   ! number of glacier pfts (columns) in glacier landunit
    integer , optional, intent(out) :: nglacier_mec  ! number of glacier_mec pfts (columns) in glacier_mec landunit
    integer , optional, intent(in)  :: glcmask  ! = 1 if glc requires surface mass balance in this gridcell
!
! !CALLED FROM:
! subroutines decomp_init, initGridCells
!
! !REVISION HISTORY:
! 2002.09.11  Mariana Vertenstein  Creation.
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m                ! loop index
    integer  :: n                ! elevation class index
    integer  :: ipfts            ! number of pfts in gridcell
    integer  :: icols            ! number of columns in gridcell
    integer  :: ilunits          ! number of landunits in gridcell
    integer  :: npfts_per_lunit  ! number of pfts in landunit
!------------------------------------------------------------------------------

    ! Initialize pfts, columns and landunits counters for gridcell

    ipfts   = 0
    icols   = 0
    ilunits = 0

    ! Set naturally vegetated landunit

    npfts_per_lunit = 0
    if (allocate_all_vegpfts) then
       ! WJS (5-8-13): For some reason, we allocate space for the natural veg landunit
       ! whenever the crop landunit has area > 0, even if there is no natural veg here.
       ! I'm not sure why this is, but I'm maintaining this behavior for now.
       if (wt_lunit(gi, istsoil) > 0.0_r8 .or. wt_lunit(gi, istcrop) > 0.0_r8) then
          npfts_per_lunit = natpft_size
       end if
    else
       write(iulog,*) 'allocate_all_vegpfts=false is no longer supported'
       call endrun()
    end if
    if (npfts_per_lunit > 0) then ! true even when only crops are present
       ! Assume that the vegetated landunit has one column
       ilunits = ilunits + 1
       icols = icols + 1  
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nveg )) nveg  = npfts_per_lunit

    ! Set urban tall building district landunit

    npfts_per_lunit = 0
    if (wt_lunit(gi, isturb_tbd) > 0.0_r8) then
       npfts_per_lunit = maxpatch_urb
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nurban_tbd )) nurban_tbd  = npfts_per_lunit

    ! Set urban high density landunit

    npfts_per_lunit = 0
    if (wt_lunit(gi, isturb_hd) > 0.0_r8) then
       npfts_per_lunit = maxpatch_urb
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nurban_hd )) nurban_hd  = npfts_per_lunit

    ! Set urban medium density landunit

    npfts_per_lunit = 0
    if (wt_lunit(gi, isturb_md) > 0.0_r8) then
       npfts_per_lunit = maxpatch_urb
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nurban_md )) nurban_md  = npfts_per_lunit

    ! Set lake landunit

    npfts_per_lunit = 0
    if (wt_lunit(gi, istdlak) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nlake )) nlake  = npfts_per_lunit

    ! Set wetland landunit

    npfts_per_lunit = 0
    if (wt_lunit(gi, istwet) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nwetland )) nwetland  = npfts_per_lunit

    ! Set glacier landunit

    npfts_per_lunit = 0
    if (wt_lunit(gi, istice) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nglacier )) nglacier  = npfts_per_lunit

    ! Set glacier_mec landunit
    ! If glcmask = 1, we create a column for each elevation class even if the weight on
    ! the grid cell is 0.

    npfts_per_lunit = 0
    do m = 1, maxpatch_glcmec
       ! If the landunit has non-zero weight on the grid cell, and this column has
       ! non-zero weight on the landunit...
       if (wt_lunit(gi, istice_mec) > 0.0_r8 .and. wt_glc_mec(gi, m) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1

       elseif (present(glcmask)) then
          if (glcmask == 1) then      ! create a virtual column 
             npfts_per_lunit = npfts_per_lunit + 1
          endif  ! glcmask = 1 
       endif  ! wt > 0
    enddo   ! maxpatch_glcmec
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nglacier_mec )) nglacier_mec  = npfts_per_lunit

    ! Set crop landunit if appropriate

    npfts_per_lunit = 0
    if (create_crop_landunit) then
       if (allocate_all_vegpfts) then
          ! WJS (5-8-13): For some reason, we allocate space for the crop landunit whenever
          ! the natural veg landunit has area > 0, even if there is no crop here.
          ! I'm not sure why this is, but I'm maintaining this behavior for now.
          if (wt_lunit(gi, istsoil) > 0.0_r8 .or. wt_lunit(gi, istcrop) > 0.0_r8) then
             npfts_per_lunit = cft_size
          end if
       else
          write(iulog,*) 'allocate_all_vegpfts=false is not supported with create_crop_landunit'
          call endrun()
       end if
    end if
    if (npfts_per_lunit > 0) then ! true even if only natural veg is present
                                  ! (as long as create_crop_landunit is true)
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(ncrop )) ncrop  = npfts_per_lunit

    ! Determine return arguments

    if (present(nlunits)) nlunits = ilunits
    if (present(ncols))   ncols   = icols
    if (present(npfts))   npfts   = ipfts

  end subroutine subgrid_get_gcellinfo

!-----------------------------------------------------------------------

end module subgridMod
