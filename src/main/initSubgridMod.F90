module initSubgridMod

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar  , only : lsmlon, lsmlat, numpft, &
	                   maxpatch, maxpatch_pft, maxpatch_cft, &
                           npatch_lake, npatch_glacier, npatch_wet, npatch_crop
  use clm_varctl  , only : allocate_all_vegpfts
  use abortutils  , only : endrun

  implicit none
  private	
  save

  public get_gcell_info                ! returns gridcell landunits, columns and pfts properties
  public set_landunit_veg_compete      ! sets up gridcell naturally vegetated landunit
  public set_landunit_wet_ice_lake     ! sets up gridcell lake, wetland and glacier landunits
  public set_landunit_crop_noncompete  ! sets up gridcell for crop landunit

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gcell_info
!
! !INTERFACE:
  subroutine get_gcell_info (i, j, wtxy, &
	                     nlunits, ncols, npfts, &
                             nveg, wtveg, &
                             ncrop, wtcrop, &
                             nlake, wtlake, &
                             nwetland, wtwetland, &
                             nglacier, wtglacier )
!
! !DESCRIPTION:
! Obtain gridcell properties
!
! !USES

! !ARGUMENTS
    implicit none
    integer , intent(in)  :: i                    ! longitude index
    integer , intent(in)  :: j                    ! latitude index 
    real(r8), intent(in)  :: wtxy(lsmlon, lsmlat, maxpatch) ! subgrid pft weights
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
! subroutines initDecomp 
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
       if (wtxy(i,j,m) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1
          wtlunit = wtlunit + wtxy(i,j,m)
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
    if (wtxy(i,j,npatch_lake) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(i,j,npatch_lake)
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
    if (wtxy(i,j,npatch_wet) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(i,j,npatch_wet)
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
    if (wtxy(i,j,npatch_glacier) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(i,j,npatch_glacier)
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
       if (wtxy(i,j,m) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1
          wtlunit = wtlunit + wtxy(i,j,m)
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

  end subroutine get_gcell_info

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_veg_compete
!
! !INTERFACE:
  subroutine set_landunit_veg_compete (ltype, wtxy, vegxy, &
                                       i, j, gi, li, ci, pi, &
                                       i_sn, j_sn, gi_sn, li_sn, ci_sn, pi_sn, &
                                       nump, numc, numl, &
            	                       land1d_sn_gi, cols1d_sn_li, cols1d_sn_gi, &
                                       pfts1d_sn_ci, pfts1d_sn_li, pfts1d_sn_gi)
!
! !DESCRIPTION: 
! Initialize vegetated landunit with competition
!
! !USES
    use clmtype
    use clm_varcon, only : istsoil
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    real(r8), intent(in)    :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch weights
    integer , intent(in)    :: vegxy(lsmlon,lsmlat,maxpatch) ! PFT type 
    integer , intent(in)    :: i                 ! 2d longitude index
    integer , intent(in)    :: j                 ! 2d latitude index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    integer , intent(in)    :: i_sn              ! 2d longitude index
    integer , intent(in)    :: j_sn              ! 2d latitude index
    integer , intent(in)    :: gi_sn             ! gridcell s->n index
    integer , intent(inout) :: li_sn             ! landunit index
    integer , intent(inout) :: ci_sn             ! column index
    integer , intent(inout) :: pi_sn             ! pft index
    integer , intent(in)    :: nump
    integer , intent(in)    :: numc
    integer , intent(in)    :: numl
    integer , intent(inout) :: land1d_sn_gi(numl)
    integer , intent(inout) :: cols1d_sn_li(numc)
    integer , intent(inout) :: cols1d_sn_gi(numc) 
    integer , intent(inout) :: pfts1d_sn_ci(nump) 
    integer , intent(inout) :: pfts1d_sn_li(nump) 
    integer , intent(inout) :: pfts1d_sn_gi(nump) 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(i,j,m)
    integer  :: n                                ! loop index
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landunit
    real(r8) :: wtlunit2gcell                    ! weight relative to gridcell of landunit
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr         ! pointer to pft derived subtype
!------------------------------------------------------------------------

    ! Set decomposition properties

    call get_gcell_info(i, j, wtxy, nveg=npfts, wtveg=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       lptr => clm3%g%l
       cptr => clm3%g%l%c
       pptr => clm3%g%l%c%p

       ! Set landunit properties
       
       ncols = 1
       
       li = li + 1   
       lptr%ncolumns(li)  = 1
       lptr%coli(li)      = ci + 1
       lptr%colf(li)      = ci + ncols
       lptr%npfts(li)     = npfts
       lptr%pfti(li)      = pi + 1
       lptr%pftf(li)      = pi + npfts 
       lptr%gridcell(li)  = gi
       lptr%wtgcell(li)   = wtlunit2gcell
       lptr%ifspecial(li) = .false.
       lptr%lakpoi(li)    = .false.
       lptr%itype(li)     = ltype
       
       ! Set column properties for this landunit (only one column on landunit)
    
       ci = ci + 1
       cptr%npfts(ci)    = npfts
       cptr%pfti(ci)     = pi + 1
       cptr%pftf(ci)     = pi + npfts
       cptr%landunit(ci) = li
       cptr%gridcell(ci) = gi
       cptr%wtlunit(ci)  = 1.0_r8
       cptr%wtgcell(ci)  = wtlunit2gcell
       cptr%itype(ci)    = 1
       
       ! Set pft properties for this landunit

       if (allocate_all_vegpfts) then
          do n = 1,numpft+1
             pi = pi + 1
             pptr%column(pi)   = ci
             pptr%landunit(pi) = li
             pptr%gridcell(pi) = gi
             pptr%wtcol(pi)    = 0._r8
             pptr%wtlunit(pi)  = 0._r8
             pptr%wtgcell(pi)  = 0._r8
             pptr%mxy(pi)      = n
             pptr%itype(pi)    = n-1
             do m = 1,maxpatch_pft
                if (vegxy(i,j,m) == pptr%itype(pi) .and. wtxy(i,j,m) > 0._r8) then
                   pptr%wtcol(pi)    = pptr%wtcol(pi)   + wtxy(i,j,m) / wtlunit2gcell
                   pptr%wtlunit(pi)  = pptr%wtlunit(pi) + wtxy(i,j,m) / wtlunit2gcell
                   pptr%wtgcell(pi)  = pptr%wtgcell(pi) + wtxy(i,j,m)
                end if
             end do
          end do
       else
          do m = 1,maxpatch_pft
             if (wtxy(i,j,m) > 0._r8) then
                pi = pi + 1
                pptr%column(pi)   = ci
                pptr%landunit(pi) = li
                pptr%gridcell(pi) = gi
                pptr%wtcol(pi)    = wtxy(i,j,m) / wtlunit2gcell
                pptr%wtlunit(pi)  = wtxy(i,j,m) / wtlunit2gcell
                pptr%wtgcell(pi)  = wtxy(i,j,m)
                pptr%mxy(pi)      = m
                pptr%itype(pi)    = vegxy(i,j,m)
             end if
          end do
       end if

    end if

    ! Set south->north properties

    call get_gcell_info(i_sn, j_sn, wtxy, nveg=npfts, wtveg=wtlunit2gcell)

    if (npfts > 0) then

       li_sn = li_sn + 1   
       land1d_sn_gi(li_sn) = gi_sn

       ci_sn = ci_sn + 1
       cols1d_sn_li(ci_sn) = li_sn
       cols1d_sn_gi(ci_sn) = gi_sn

       if (allocate_all_vegpfts) then
          do n = 1,numpft+1
             pi_sn = pi_sn + 1
             pfts1d_sn_ci(pi_sn) = ci_sn
             pfts1d_sn_li(pi_sn) = li_sn
             pfts1d_sn_gi(pi_sn) = gi_sn
          end do
       else
          do m = 1,maxpatch_pft
             if (wtxy(i_sn,j_sn,m) > 0._r8) then
                pi_sn = pi_sn + 1
                pfts1d_sn_ci(pi_sn) = ci_sn
                pfts1d_sn_li(pi_sn) = li_sn
                pfts1d_sn_gi(pi_sn) = gi_sn
             end if
          end do
       end if

    end if

  end subroutine set_landunit_veg_compete
  
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_wet_ice_lake
!
! !INTERFACE:
  subroutine set_landunit_wet_ice_lake (ltype, wtxy, vegxy, &
                                   i, j, gi, li, ci, pi, &
                                   i_sn, j_sn, gi_sn, li_sn, ci_sn, pi_sn, &
                                   nump, numc, numl, &
                                   land1d_sn_gi, cols1d_sn_li, cols1d_sn_gi, &
                                   pfts1d_sn_ci, pfts1d_sn_li, pfts1d_sn_gi)
!
! !DESCRIPTION: 
! Initialize wet_ice_lake landunits that are non-urban (lake, wetland, glacier)
!
! !USES
    use clmtype
    use clm_varcon, only : istice, istwet, istdlak
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    real(r8), intent(in)    :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch weights
    integer , intent(in)    :: vegxy(lsmlon,lsmlat,maxpatch) ! PFT type 
    integer , intent(in)    :: i               ! 2-dim longitude index
    integer , intent(in)    :: j               ! 2-dim latitude index
    integer , intent(in)    :: gi              ! gridcell index
    integer , intent(inout) :: li              ! landunit index
    integer , intent(inout) :: ci              ! column index
    integer , intent(inout) :: pi              ! pft index
    integer , intent(in)    :: i_sn            ! 2-dim longitude index
    integer , intent(in)    :: j_sn            ! 2-dim latitude index
    integer , intent(in)    :: gi_sn           ! gridcell s->n index
    integer , intent(inout) :: li_sn           ! landunit index
    integer , intent(inout) :: ci_sn           ! column index
    integer , intent(inout) :: pi_sn           ! pft index
    integer , intent(in)    :: nump
    integer , intent(in)    :: numc
    integer , intent(in)    :: numl
    integer , intent(inout) :: land1d_sn_gi(numl)
    integer , intent(inout) :: cols1d_sn_li(numc)
    integer , intent(inout) :: cols1d_sn_gi(numc) 
    integer , intent(inout) :: pfts1d_sn_ci(nump) 
    integer , intent(inout) :: pfts1d_sn_li(nump) 
    integer , intent(inout) :: pfts1d_sn_gi(nump) 
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(i,j,m)
    integer  :: c                                ! column loop index
    integer  :: ctype                            ! column type
    integer  :: ier                              ! error status 
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landunit
    real(r8) :: wtlunit2gcell                    ! weight relative to gridcell of landunit
    real(r8) :: wtcol2lunit                      ! temporary weight with respect to landunit
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr         ! pointer to pft derived subtype
!------------------------------------------------------------------------

    ! Set decomposition properties

    if (ltype == istwet) then
       call get_gcell_info(i, j, wtxy, nwetland=npfts, wtwetland=wtlunit2gcell)
       m = npatch_wet
    else if (ltype == istdlak) then
       call get_gcell_info(i, j, wtxy, nlake=npfts, wtlake=wtlunit2gcell)
       m = npatch_lake
    else if (ltype == istice) then 
       call get_gcell_info(i, j, wtxy, nglacier=npfts, wtglacier=wtlunit2gcell)
       m = npatch_glacier
    else
       write(6,*)' set_landunit_wet_ice_lake: ltype of ',ltype,' not valid'
       write(6,*)' only istwet, istdlak or istice ltypes are valid'
       call endrun()
    end if

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       lptr => clm3%g%l
       cptr => clm3%g%l%c
       pptr => clm3%g%l%c%p
       
       if (npfts /=1) then
          write(6,*)' set_landunit_wet_ice_lake: compete landunit must have one column and one pft '
          write(6,*)' current values of ncols, pfts=',ncols,npfts
          call endrun()
       end if

       ncols = 1

       ! Currently assume that each landunit only has only one column 
       ! (of type 1) and that each column has its own pft
       
       wtcol2lunit = 1.0_r8/ncols
       ctype = 1

       ! Determine landunit properties 

       li = li + 1
       lptr%ncolumns(li)  = ncols
       lptr%coli(li)      = ci + 1
       lptr%colf(li)      = ci + ncols
       lptr%npfts(li)     = ncols
       lptr%pfti(li)      = pi + 1
       lptr%pftf(li)      = pi + npfts
       lptr%gridcell(li)  = gi
       lptr%wtgcell(li)   = wtlunit2gcell
       lptr%itype(li)     = ltype
       lptr%ifspecial(li) = .true.
       if (ltype == istdlak) then
          lptr%lakpoi(li) = .true.
       else
          lptr%lakpoi(li) = .false.
       end if
       
       ! Determine column and properties
       ! For the wet, ice or lake landunits it is assumed that each column has its own pft
       
       ci = ci + 1
       pi = pi + 1 
       
       cptr%npfts(ci)    = 1
       cptr%pfti(ci)     = pi 
       cptr%pftf(ci)     = pi 
       cptr%landunit(ci) = li
       cptr%gridcell(ci) = gi
       cptr%wtlunit(ci)  = wtcol2lunit
       cptr%wtgcell(ci)  = wtcol2lunit * wtlunit2gcell
       cptr%itype(ci)    = ctype
       
       pptr%column(pi)   = ci
       pptr%landunit(pi) = li
       pptr%gridcell(pi) = gi
       pptr%wtcol(pi)    = 1._r8
       pptr%wtlunit(pi)  = wtcol2lunit
       pptr%wtgcell(pi)  = wtcol2lunit * wtlunit2gcell
       pptr%mxy(pi)      = m
       pptr%itype(pi)    = vegxy(i,j,m)
       
    end if
       
    ! Set south->north properties

    if (ltype == istwet) then
       call get_gcell_info(i_sn, j_sn, wtxy, nwetland=npfts, wtwetland=wtlunit2gcell)
    else if (ltype == istdlak) then
       call get_gcell_info(i_sn, j_sn, wtxy, nlake=npfts, wtlake=wtlunit2gcell)
    else if (ltype == istice) then 
       call get_gcell_info(i_sn, j_sn, wtxy, nglacier=npfts, wtglacier=wtlunit2gcell)
    end if

    if (npfts > 0) then

       li_sn = li_sn + 1
       land1d_sn_gi(li_sn) = gi_sn
       
       ci_sn = ci_sn + 1
       cols1d_sn_gi(ci_sn)  = gi_sn
       cols1d_sn_li(ci_sn)  = li_sn
       
       pi_sn = pi_sn + 1
       pfts1d_sn_ci(pi_sn) = ci_sn
       pfts1d_sn_li(pi_sn) = li_sn
       pfts1d_sn_gi(pi_sn) = gi_sn

    end if

  end subroutine set_landunit_wet_ice_lake

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_crop_noncompete
!
! !INTERFACE:
  subroutine set_landunit_crop_noncompete (ltype, wtxy, vegxy, &
                                           i, j, gi, li, ci, pi, &
                                           i_sn, j_sn, gi_sn, li_sn, ci_sn, pi_sn, &
                                           nump, numc, numl, &
                                           land1d_sn_gi, cols1d_sn_li, cols1d_sn_gi, &
                                           pfts1d_sn_ci, pfts1d_sn_li, pfts1d_sn_gi)
!
! !DESCRIPTION: 
! Initialize crop landunit without competition
!
! !USES
    use clmtype
    use clm_varpar, only : npatch_crop, npatch_glacier
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    real(r8), intent(in)    :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch weights
    integer , intent(in)    :: vegxy(lsmlon,lsmlat,maxpatch) ! PFT types 
    integer , intent(in)    :: i                 ! 2d longitude index
    integer , intent(in)    :: j                 ! 2d latitude index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    integer , intent(in)    :: i_sn              ! 2d longitude index
    integer , intent(in)    :: j_sn              ! 2d latitude index
    integer , intent(in)    :: gi_sn             ! gridcell s->n index
    integer , intent(inout) :: li_sn             ! landunit index
    integer , intent(inout) :: ci_sn             ! column index
    integer , intent(inout) :: pi_sn             ! pft index
    integer , intent(in)    :: nump
    integer , intent(in)    :: numc
    integer , intent(in)    :: numl
    integer , intent(inout) :: land1d_sn_gi(numl)
    integer , intent(inout) :: cols1d_sn_li(numc)
    integer , intent(inout) :: cols1d_sn_gi(numc) 
    integer , intent(inout) :: pfts1d_sn_ci(nump) 
    integer , intent(inout) :: pfts1d_sn_li(nump) 
    integer , intent(inout) :: pfts1d_sn_gi(nump) 
!
! !REVISION HISTORY:
! Created by Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(i,j,m)
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landunit
    real(r8) :: wtlunit2gcell                    ! weight relative to gridcell of landunit
    type(landunit_type), pointer :: lptr         ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr         ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr         ! pointer to pft derived subtype
!------------------------------------------------------------------------

    ! Set decomposition properties

    call get_gcell_info(i, j, wtxy, ncrop=npfts, wtcrop=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module
       
       lptr => clm3%g%l
       cptr => clm3%g%l%c
       pptr => clm3%g%l%c%p
       
       ! Set landunit properties - each column has its own pft
       
       ncols = npfts
       
       li = li + 1   
       lptr%ncolumns(li)  = ncols
       lptr%coli(li)      = ci + 1
       lptr%colf(li)      = ci + ncols
       lptr%npfts(li)     = npfts
       lptr%pfti(li)      = pi + 1
       lptr%pftf(li)      = pi + npfts
       lptr%gridcell(li)  = gi
       lptr%wtgcell(li)   = wtlunit2gcell
       lptr%itype(li)     = ltype
       lptr%ifspecial(li) = .false.
       lptr%lakpoi(li)    = .false.
       
       ! Set column and pft properties for this landunit (each column has its own pft)

       do m = npatch_glacier+1, npatch_crop
          if (wtxy(i,j,m) > 0._r8) then
             ci = ci + 1
             pi = pi + 1
             
             cptr%gridcell(ci) = gi
             cptr%landunit(ci) = li
             cptr%npfts(ci)    = 1
             cptr%pfti(ci)     = pi 
             cptr%pftf(ci)     = pi 
             cptr%wtlunit(ci)  = wtxy(i,j,m) / wtlunit2gcell
             cptr%wtgcell(ci)  = wtxy(i,j,m)
             cptr%itype(ci)    = 1
             
             pptr%gridcell(pi) = gi
             pptr%landunit(pi) = li
             pptr%column(pi)   = ci
             pptr%wtcol(pi)    = 1.0_r8
             pptr%wtlunit(pi)  = wtxy(i,j,m) / wtlunit2gcell
             pptr%wtgcell(pi)  = wtxy(i,j,m)
             pptr%itype(pi)    = vegxy(i,j,m)
             pptr%mxy(pi)      = m
          end if
       end do

    end if
       
    ! Set south->north properties

    call get_gcell_info(i_sn, j_sn, wtxy, ncrop=npfts, wtcrop=wtlunit2gcell)

    if (npfts > 0) then

       li_sn = li_sn + 1
       land1d_sn_gi(li_sn) = gi_sn
       
       do m = npatch_glacier+1, npatch_crop
          if (wtxy(i_sn,j_sn,m) > 0._r8) then
             ci_sn = ci_sn + 1
             pi_sn = pi_sn + 1
             cols1d_sn_gi(ci_sn) = gi_sn
             cols1d_sn_li(ci_sn) = li_sn
             pfts1d_sn_ci(pi_sn) = ci_sn
             pfts1d_sn_li(pi_sn) = li_sn
             pfts1d_sn_gi(pi_sn) = gi_sn
          end if
       end do  

    end if

  end subroutine set_landunit_crop_noncompete

end module initSubgridMod
