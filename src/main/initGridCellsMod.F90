#include <misc.h>
#include <preproc.h>

module initGridcellsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initGridcellsMod
!
! !DESCRIPTION:
! Initializes sub-grid mapping for each land grid cell
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use abortutils  , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public initGridcells ! initialize sub-grid gridcell mapping 
  public get_sn_land1d ! returns s->n gridcell indices for 
                       ! each s->n landunit
  public get_sn_cols1d ! returns s->n gridcell or landunit indices for 
                       ! each s->n column
  public get_sn_pfts1d ! returns s->n gridcell, landunit or column indices for 
                       ! s->n pft column
!
! !PIVATE MEMBER FUNCTIONS:
  private set_gcell_otherprops
  private set_gcell_gptr
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL MODULE VARIABLES:
  integer, private, allocatable :: land1d_sn_gi(:)  ! corresponding s->n gridcell index for every s->n landunit
  integer, private, allocatable :: cols1d_sn_li(:)  ! corresponding s->n landunit index for every s->n column
  integer, private, allocatable :: cols1d_sn_gi(:)  ! corresponding s->n gridcell index for every s->n column
  integer, private, allocatable :: pfts1d_sn_ci(:)  ! corresponding s->n column   index for every s->n pft
  integer, private, allocatable :: pfts1d_sn_li(:)  ! corresponding s->n landunit index for every s->n pft
  integer, private, allocatable :: pfts1d_sn_gi(:)  ! corresponding s->n gridcell index for every s->n pft
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initGridcells
!
! !INTERFACE:
  subroutine initGridcells (vegxy, wtxy) 
!
! !DESCRIPTION: 
! Initialize sub-grid mapping and allocates space for derived type hierarchy.
! For each land gridcell determine landunit, column and pft properties.
!
! !USES
    use clmtype
    use clm_varpar      , only : lsmlon, lsmlat, maxpatch 
    use clm_varsur      , only : landmask, numlon 
    use clm_varcon      , only : istsoil, istice, istwet, istdlak, isturb
    use initGridIndexMod, only : set_index_xy, set_index_sn
    use initSubgridMod  , only : get_gcell_info, &
                                 set_landunit_veg_compete, &
                                 set_landunit_wet_ice_lake, &
    	                         set_landunit_crop_noncompete  
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) ! PFT type 
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch weights
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,i_sn,j_sn              ! indices 
    integer :: gi_sn,li_sn, ci_sn, pi_sn  ! indices
    integer :: gi,li,ci,pi,m              ! indices
    integer :: nveg                       ! number of pfts in naturally vegetated landunit
    real(r8):: wtveg                      ! weight (relative to gridcell) of naturally veg landunit
    integer :: ncrop                      ! number of crop pfts in crop landunit
    real(r8):: wtcrop                     ! weight (relative to gridcell) of crop landunit
    integer :: nlake                      ! number of pfts (columns) in lake landunit
    real(r8):: wtlake                     ! weight (relative to gridcell) of lake landunit
    integer :: nwetland                   ! number of pfts (columns) in wetland landunit
    real(r8):: wtwetland                  ! weight (relative to gridcell) of wetland landunit
    integer :: nglacier                   ! number of pfts (columns) in glacier landunit
    real(r8):: wtglacier                  ! weight (relative to gridcell) of glacier landunit
    integer :: ilunits                    ! temporary
    integer :: icols                      ! temporary
    integer :: ipfts                      ! temporary
    integer :: ier                        ! error status
    integer :: numg                       ! total number of gridcells across all processors
    integer :: numl                       ! total number of landunits across all processors 
    integer :: numc                       ! total number of columns across all processors 
    integer :: nump                       ! total number of pfts across all processors 
    integer, pointer :: ixy(:)            ! lon index for current decomposition
    integer, pointer :: jxy(:)            ! lat index for current decomposition
    integer, pointer :: ixy_sn(:)         ! lon index for s->n ordering
    integer, pointer :: jxy_sn(:)         ! lat index for s->n ordering
    integer, pointer :: gridcell_sn(:)    ! gridcell index for s->n ordering
    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
 !------------------------------------------------------------------------

    ! Set pointers into derived types for this module

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Find total global number of grid cells, landunits, columns and pfts 
    
    numg = 0
    numl = 0
    numc = 0
    nump = 0
    do j = 1, lsmlat
       do i = 1, numlon(j)
          if (landmask(i,j) == 1) then         
             call get_gcell_info (i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)
             numg = numg + 1
             numl = numl + ilunits
             numc = numc + icols
             nump = nump + ipfts
          end if
       end do
    end do

    ! Dynamic memory allocation

    allocate(land1d_sn_gi(numl), &
             cols1d_sn_gi(numc), cols1d_sn_li(numc), &
             pfts1d_sn_gi(nump), pfts1d_sn_li(nump), pfts1d_sn_ci(nump), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initGridCells(): allocation error for sn indices'
       call endrun()
    end if

    allocate(ixy(numg), jxy(numg), ixy_sn(numg), jxy_sn(numg), gridcell_sn(numg), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initGridCells(): allocation error for sn indices'
       call endrun()
    end if

    ! Determine i,j index for each land gridcell 
    ! (both in current decompsition and for s->n gridcell ordering)

    call set_index_xy(numg, ixy, jxy, ixy_sn, jxy_sn)

    ! Determine SN gridcell ordering

    call set_index_sn(numg, gridcell_sn, type1d=nameg)

    ! For each land gridcell on global grid determine landunit, column and pft properties

    li    = 0
    ci    = 0
    pi    = 0
    li_sn = 0
    ci_sn = 0
    pi_sn = 0

    do gi = 1,numg

       ! Determine indices

       i     = ixy(gi)
       j     = jxy(gi)
       i_sn  = ixy_sn(gi)
       j_sn  = jxy_sn(gi)
       gi_sn = gi

       ! Determine naturally vegetated landunit

       call set_landunit_veg_compete( &
            ltype=istsoil, wtxy=wtxy, vegxy=vegxy, &
            i=i      , j=j      , gi=gi      , li=li      , ci=ci      , pi=pi, &
            i_sn=i_sn, j_sn=j_sn, gi_sn=gi_sn, li_sn=li_sn, ci_sn=ci_sn, pi_sn=pi_sn, &
            nump=nump, numc=numc, numl=numl, &
            land1d_sn_gi=land1d_sn_gi, cols1d_sn_li=cols1d_sn_li, cols1d_sn_gi=cols1d_sn_gi, &
            pfts1d_sn_ci=pfts1d_sn_ci, pfts1d_sn_li=pfts1d_sn_li, pfts1d_sn_gi=pfts1d_sn_gi)

       ! Determine crop landunit

       call set_landunit_crop_noncompete( &
            ltype=istsoil, wtxy=wtxy, vegxy=vegxy,  &
            i=i      , j=j      , gi=gi      , li=li      , ci=ci      , pi=pi, &
            i_sn=i_sn, j_sn=j_sn, gi_sn=gi_sn, li_sn=li_sn, ci_sn=ci_sn, pi_sn=pi_sn, &
            nump=nump, numc=numc, numl=numl, &
            land1d_sn_gi=land1d_sn_gi, cols1d_sn_li=cols1d_sn_li, cols1d_sn_gi=cols1d_sn_gi, &
            pfts1d_sn_ci=pfts1d_sn_ci, pfts1d_sn_li=pfts1d_sn_li, pfts1d_sn_gi=pfts1d_sn_gi)

       ! Determine lake, wetland and glacier landunits 

       call set_landunit_wet_ice_lake( &
            ltype=istdlak, wtxy=wtxy, vegxy=vegxy,  &
            i=i      , j=j      , gi=gi      , li=li      , ci=ci      , pi=pi, &
            i_sn=i_sn, j_sn=j_sn, gi_sn=gi_sn, li_sn=li_sn, ci_sn=ci_sn, pi_sn=pi_sn, &
            nump=nump, numc=numc, numl=numl, &
            land1d_sn_gi=land1d_sn_gi, cols1d_sn_li=cols1d_sn_li, cols1d_sn_gi=cols1d_sn_gi, &
            pfts1d_sn_ci=pfts1d_sn_ci, pfts1d_sn_li=pfts1d_sn_li, pfts1d_sn_gi=pfts1d_sn_gi)
       
       call set_landunit_wet_ice_lake( &
            ltype=istwet, wtxy=wtxy, vegxy=vegxy,  &
            i=i      , j=j      , gi=gi      , li=li      , ci=ci      , pi=pi, &
            i_sn=i_sn, j_sn=j_sn, gi_sn=gi_sn, li_sn=li_sn, ci_sn=ci_sn, pi_sn=pi_sn, &
            nump=nump, numc=numc, numl=numl, &
            land1d_sn_gi=land1d_sn_gi, cols1d_sn_li=cols1d_sn_li, cols1d_sn_gi=cols1d_sn_gi, &
            pfts1d_sn_ci=pfts1d_sn_ci, pfts1d_sn_li=pfts1d_sn_li, pfts1d_sn_gi=pfts1d_sn_gi)
       
       call set_landunit_wet_ice_lake( &
            ltype=istice, wtxy=wtxy, vegxy=vegxy,  &
            i=i      , j=j      , gi=gi      , li=li      , ci=ci      , pi=pi, &
            i_sn=i_sn, j_sn=j_sn, gi_sn=gi_sn, li_sn=li_sn, ci_sn=ci_sn, pi_sn=pi_sn, &
            nump=nump, numc=numc, numl=numl, &
            land1d_sn_gi=land1d_sn_gi, cols1d_sn_li=cols1d_sn_li, cols1d_sn_gi=cols1d_sn_gi, &
            pfts1d_sn_ci=pfts1d_sn_ci, pfts1d_sn_li=pfts1d_sn_li, pfts1d_sn_gi=pfts1d_sn_gi)

    end do  

    ! Determine gridcell gptr properties

    call set_gcell_gptr(numg, ixy, jxy, wtxy)

    ! Set xy grid info (ixy, jxy, lat and lon) and areas for each subgrid type

    call set_gcell_otherprops(numg, ixy, jxy)

    ! Set south->north indices (the following are global quantities)

    call set_index_sn(numg, gptr%snindex, type1d=nameg)
    call set_index_sn(numg, lptr%snindex, type1d=namel)
    call set_index_sn(numg, cptr%snindex, type1d=namec)
    call set_index_sn(numg, pptr%snindex, type1d=namep)

    ! Deallocate dynamic memory

    deallocate(ixy, jxy, ixy_sn, jxy_sn, gridcell_sn)

  end subroutine initGridcells

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_gcell_gptr
!
! !INTERFACE:
  subroutine set_gcell_gptr(numg, ixy, jxy, wtxy)
!
! !DESCRIPTION: 
! Initialize clmtype gptr properties
!
! !USES
    use clmtype
    use clm_varpar    , only : lsmlon, lsmlat, maxpatch 
    use initSubGridMod, only : get_gcell_info	
    implicit none
!
! !ARGUMENTS:
    integer , intent(in) :: numg
    integer , intent(in) :: ixy(numg)
    integer , intent(in) :: jxy(numg)
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch weights
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 08/2004
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g       ! indices
    integer :: ilunits     ! temporary
    integer :: icols       ! temporary
    integer :: ipfts       ! temporary
    integer :: ngcells     ! temporary
    integer :: nlunits     ! temporary
    integer :: ncols       ! temporary
    integer :: npfts       ! temporary
    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
!------------------------------------------------------------------------

    ! Set pointers into derived types for this module

    gptr => clm3%g

    ! Determine gridcell properties (currently only one type of gridcell)

    ngcells = 0
    nlunits = 0
    ncols   = 0
    npfts   = 0

    do g = 1,numg

       i = ixy(g)
       j = jxy(g)

       gptr%luni(g) = nlunits + 1
       gptr%coli(g) = ncols   + 1
       gptr%pfti(g) = npfts   + 1

       call get_gcell_info(i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)

       ngcells = ngcells + 1
       nlunits = nlunits + ilunits
       ncols   = ncols   + icols
       npfts   = npfts   + ipfts

       gptr%lunf(g) = nlunits
       gptr%colf(g) = ncols
       gptr%pftf(g) = npfts

       gptr%nlandunits(g) = gptr%lunf(g) - gptr%luni(g) + 1
       gptr%ncolumns(g)   = gptr%colf(g) - gptr%coli(g) + 1
       gptr%npfts(g)      = gptr%pftf(g) - gptr%pfti(g) + 1

       gptr%itype(g) = 1

    end do

  end subroutine set_gcell_gptr

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_gcell_otherprops
!
! !INTERFACE:
  subroutine set_gcell_otherprops(numg, ixy, jxy)
!
! !DESCRIPTION: 
! Initialize clmtype gptr, lptr, cptr and pptr component arrays
! [ixy, jxy, latdeg, londeg, area]
!
! !USES
    use clmtype
    use shr_const_mod, only : SHR_CONST_PI
    use clm_varsur   , only : area, latixy, longxy, landfrac
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: numg
    integer, intent(in) :: ixy(numg)
    integer, intent(in) :: jxy(numg)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 08/2004
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,l,c,p
    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
!------------------------------------------------------------------------

    ! Set pointers into derived types for this module

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p
    
    ! Loop over all gridcells

    do g = 1,numg

       i = ixy(g)
       j = jxy(g)

       gptr%ixy(g) = i 
       gptr%jxy(g) = j 
       gptr%latdeg(g) = latixy(i,j) 
       gptr%londeg(g) = longxy(i,j) 
       gptr%lat(g)    = latixy(i,j) * SHR_CONST_PI/180._r8  
       gptr%lon(g)    = longxy(i,j) * SHR_CONST_PI/180._r8
       gptr%landfrac(g) = landfrac(i,j)
       gptr%area(g)  = area(i,j)

       do l = gptr%luni(g), gptr%lunf(g)

          lptr%ixy(l) = i
          lptr%jxy(l) = j
          lptr%latdeg(l) = latixy(i,j) 
          lptr%londeg(l) = longxy(i,j) 
          lptr%area(l)   = lptr%wtgcell(l) * area(i,j)

          do c = lptr%coli(l), lptr%colf(l)

             cptr%ixy(c) = i
             cptr%jxy(c) = j
             cptr%latdeg(c) = latixy(i,j) 
             cptr%londeg(c) = longxy(i,j) 
             cptr%area(c)   = cptr%wtgcell(c) * area(i,j)

             do p = cptr%pfti(c), cptr%pftf(c)

                pptr%ixy(p) = i
                pptr%jxy(p) = j
                pptr%latdeg(p) = latixy(i,j) 
                pptr%londeg(p) = longxy(i,j) 
                pptr%area(p)   = pptr%wtgcell(p) * area(i,j)

             end do
          end do
       end do
    end do

  end subroutine set_gcell_otherprops

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_sn_land1d
!
! !INTERFACE:
   subroutine get_sn_land1d(snindex, type1d, numl)
!
! !DESCRIPTION:
!  Creates south-> north indices at a given clm level
!
! !USES 
     use clmtype	
!
! !ARGUMENTS
     implicit none
     integer         , pointer    :: snindex(:)
     character(len=*), intent(in) :: type1d 
     integer         , intent(in) :: numl
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES
     integer :: l
!------------------------------------------------------------------------------

     do l = 1,numl
        if (type1d == nameg) then
           snindex(l) = land1d_sn_gi(l)
        end if
     end do

   end subroutine get_sn_land1d

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_sn_cols1d
!
! !INTERFACE:
   subroutine get_sn_cols1d(snindex, type1d, numc)
!
! !DESCRIPTION:
!  Creates south->north indices at a given clm level
!
! !USES 
     use clmtype	
!
! !ARGUMENTS
     implicit none
     integer         , pointer    :: snindex(:)
     character(len=*), intent(in) :: type1d 
     integer         , intent(in) :: numc
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES
     integer :: c
!------------------------------------------------------------------------------

     do c = 1,numc
        if (type1d == nameg) then
           snindex(c) = cols1d_sn_gi(c)
        else if (type1d == namel) then
           snindex(c) = cols1d_sn_li(c)
        end if
     end do

   end subroutine get_sn_cols1d

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_sn_pfts1d
!
! !INTERFACE:
   subroutine get_sn_pfts1d(snindex, type1d, nump)
!
! !DESCRIPTION:
!  Creates south-> north indices at a given clm level
!
! !USES 
     use clmtype	
!
! !ARGUMENTS
     implicit none
     integer         , pointer    :: snindex(:)
     character(len=*), intent(in) :: type1d 
     integer         , intent(in) :: nump
!
! !REVISION HISTORY:
! 2003.12.01  Mariana Vertenstein  Creation.
!
!EOP
!
! !LOCAL VARIABLES
     integer :: p
!------------------------------------------------------------------------------

     do p = 1,nump
        if (type1d == nameg) then
           snindex(p) = pfts1d_sn_gi(p)
        else if (type1d == namel) then
           snindex(p) = pfts1d_sn_li(p)
        else if (type1d == namec) then
           snindex(p) = pfts1d_sn_ci(p)
        end if
     end do

   end subroutine get_sn_pfts1d

end module initGridcellsMod
