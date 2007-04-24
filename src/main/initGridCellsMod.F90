#include <misc.h>
#include <preproc.h>

module initGridCellsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initGridCellsMod
!
! !DESCRIPTION:
! Initializes sub-grid mapping for each land grid cell
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc,iam,mpicom
  use abortutils  , only : endrun
  use clm_varsur  , only : wtxy, vegxy
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public initGridcells ! initialize sub-grid gridcell mapping 
!
! !PRIVATE MEMBER FUNCTIONS:
  private clm_ptrs_compdown
  private clm_ptrs_check
  private set_landunit_veg_compete
  private set_landunit_wet_ice_lake
  private set_landunit_crop_noncompete
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL MODULE VARIABLES:
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initGridcells
!
! !INTERFACE:
  subroutine initGridcells () 
!
! !DESCRIPTION: 
! Initialize sub-grid mapping and allocates space for derived type hierarchy.
! For each land gridcell determine landunit, column and pft properties.
!
! !USES
    use clmtype     , only : clm3, gridcell_type, landunit_type, &
                             column_type, pft_type
    use domainMod   , only : ldomain, adomain, gatm
    use decompMod   , only : ldecomp, adecomp, get_proc_global, get_proc_bounds
    use clm_varcon  , only : istsoil, istice, istwet, istdlak, isturb
    use subgridMod  , only : gcelldc, gcellsn, subgrid_alloc
    use subgridMod  , only : subgrid_get_gcellinfo
    use shr_const_mod,only : SHR_CONST_PI
    use spmdGathScatMod, only : gather_data_to_master
!
! !ARGUMENTS:
    implicit none
!    integer , intent(in) :: vegxy(:,:) ! PFT type 
!    real(r8), intent(in) :: wtxy(:,:)  ! subgrid patch weights
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: li,ci,pi,m,na,gdc,gsn,glo    ! indices
    integer :: nveg           ! number of pfts in naturally vegetated landunit
    real(r8):: wtveg          ! weight (gridcell) of naturally veg landunit
    integer :: ncrop          ! number of crop pfts in crop landunit
    real(r8):: wtcrop         ! weight (gridcell) of crop landunit
    integer :: nlake          ! number of pfts (columns) in lake landunit
    real(r8):: wtlake         ! weight (gridcell) of lake landunit
    integer :: nwetland       ! number of pfts (columns) in wetland landunit
    real(r8):: wtwetland      ! weight (gridcell) of wetland landunit
    integer :: nglacier       ! number of pfts (columns) in glacier landunit
    real(r8):: wtglacier      ! weight (gridcell) of glacier landunit
    integer :: ier            ! error status
    integer :: numg           ! total number of gridcells across all processors
    integer :: numl           ! total number of landunits across all processors 
    integer :: numc           ! total number of columns across all processors 
    integer :: nump           ! total number of pfts across all processors 
    integer :: begg,endg      ! local beg/end gridcells gdc
    integer :: begl,endl      ! local beg/end landunits
    integer :: begc,endc      ! local beg/end columns 
    integer :: begp,endp      ! local beg/end pfts
    logical :: my_gcell       ! is gdc gridcell on my pe
    integer :: nwtxy          ! wtxy cell index

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

    ! Get total global number of grid cells, landunits, columns and pfts 
    
    call get_proc_global(numg,numl,numc,nump)
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Dynamic memory allocation

    call subgrid_alloc(gcelldc,numg,numl,numc,nump,'gcelldc')
    call subgrid_alloc(gcellsn,numg,numl,numc,nump,'gcellsn')

    ! For each land gridcell on global grid determine landunit, column and pft properties

#if (1 == 1)
    li    = begl-1
    ci    = begc-1
    pi    = begp-1

    !----- Set clm3 variables -----
    do gdc = begg,endg
#else
    li    = 0
    ci    = 0
    pi    = 0

    !----- Set clm3 variables -----
    do gdc = 1,numg
#endif

       glo = ldecomp%gdc2glo(gdc)
!       nwtxy = glo
       nwtxy = gdc

       my_gcell = .false.
       if (gdc >= begg .and. gdc <= endg) then
          my_gcell = .true.
       endif

       ! Determine naturally vegetated landunit

       call set_landunit_veg_compete(               &
!            ltype=istsoil, wtxy=wtxy, vegxy=vegxy,  &
            ltype=istsoil, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine crop landunit

       call set_landunit_crop_noncompete(           &
!            ltype=istsoil, wtxy=wtxy, vegxy=vegxy,  &
            ltype=istsoil, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine lake, wetland and glacier landunits 

       call set_landunit_wet_ice_lake(              &
!            ltype=istdlak, wtxy=wtxy, vegxy=vegxy,  &
            ltype=istdlak, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       call set_landunit_wet_ice_lake(              &
!            ltype=istwet, wtxy=wtxy, vegxy=vegxy,   &
            ltype=istwet, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       call set_landunit_wet_ice_lake(              &
!            ltype=istice, wtxy=wtxy, vegxy=vegxy,   &
            ltype=istice, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Set clm3 lats/lons

       if (my_gcell) then
          gptr%latdeg(gdc) = ldomain%latc(gdc) 
          gptr%londeg(gdc) = ldomain%lonc(gdc) 
          gptr%lat(gdc)    = gptr%latdeg(gdc) * SHR_CONST_PI/180._r8  
          gptr%lon(gdc)    = gptr%londeg(gdc) * SHR_CONST_PI/180._r8
          gptr%area(gdc)   = ldomain%area(gdc)

!          na = adecomp%glo2gdc(ldomain%gatm(gdc))
          na = adecomp%glo2gdc(gatm(glo))
          gptr%londeg_a(gdc) = adomain%lonc(na)
          gptr%latdeg_a(gdc) = adomain%latc(na)
          gptr%lon_a   (gdc) = gptr%londeg_a(gdc) * SHR_CONST_PI/180._r8  
          gptr%lat_a   (gdc) = gptr%latdeg_a(gdc) * SHR_CONST_PI/180._r8  
       endif

    enddo

    ! Fill in subgrid datatypes

    call clm_ptrs_compdown()
    call clm_ptrs_check()

    ! Compute gcellsn indexes

    call gather_data_to_master (gptr%luni, gcelldc%g_li, clmlevel='gridcell')
    call gather_data_to_master (gptr%coli, gcelldc%g_ci, clmlevel='gridcell')
    call gather_data_to_master (gptr%pfti, gcelldc%g_pi, clmlevel='gridcell')
    call gather_data_to_master (gptr%nlandunits, gcelldc%g_ln, clmlevel='gridcell')
    call gather_data_to_master (gptr%ncolumns, gcelldc%g_cn, clmlevel='gridcell')
    call gather_data_to_master (gptr%npfts, gcelldc%g_pn, clmlevel='gridcell')

    li    = 1
    ci    = 1
    pi    = 1
    do gsn = 1,numg
       gdc = ldecomp%gsn2gdc(gsn)
       gcellsn%g_li(gsn) = li
       gcellsn%g_ci(gsn) = ci
       gcellsn%g_pi(gsn) = pi
       li = li + gcelldc%g_ln(gdc)
       ci = ci + gcelldc%g_cn(gdc)
       pi = pi + gcelldc%g_pn(gdc)
    enddo

  end subroutine initGridcells

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_ptrs_compdown
!
! !INTERFACE:
  subroutine clm_ptrs_compdown()
!
! !DESCRIPTION:
! Assumes the part of the subgrid pointing up has been set.  Fills 
! in the data pointing down.  Up is p_c, p_l, p_g, c_l, c_g, and l_g.
!
! This algorithm assumes all indices are monotonically increasing.
!
! Algorithm works as follows.  The p, c, and l loops march through
! the full arrays (nump, numc, and numl) checking the "up" indexes.
! As soon as the "up" index of the current (p,c,l) cell changes relative 
! to the previous (p,c,l) cell, the *i array will be set to point down 
! to that cell.  The *f array follows the same logic, so it's always the
! last "up" index from the previous cell when an "up" index changes.
!
! For example, a case where p_c(1:4) = 1 and p_c(5:12) = 2.  This 
! subroutine will set c_pi(1) = 1, c_pf(1) = 4, c_pi(2) = 5, c_pf(2) = 12.
!
! !USES
    use clmtype, only : clm3, gridcell_type, landunit_type, &
                        column_type, pft_type
    use decompMod , only : get_proc_bounds

! !ARGUMENTS
    implicit none
!
! !CALLED FROM:
! subroutines initGridCellsMod
!
! !REVISION HISTORY:
! 2005.11.15  T Craig Creation
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: begg,endg,begl,endl,begc,endc,begp,endp ! beg/end glcp
    integer :: g,l,c,p               ! loop counters
    integer :: curg,curl,curc,curp   ! tracks g,l,c,p indexes in arrays
    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
!------------------------------------------------------------------------------

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !--- Set the current c,l,g (curc, curl, curg) to zero for initialization,
    !---   these indices track the current "up" index.
    !--- Take advantage of locality of g/l/c/p cells
    !--- Loop p through full local begp:endp length
    !--- Separately check the p_c, p_l, and p_g indexes for a change in
    !---   the "up" index.
    !--- If there is a change, verify that the current c,l,g is within the 
    !---   valid range, and set c_pi, l_pi, or g_pi to that current c,l,g
    !--- Constantly update the c_pf, l_pf, and g_pf array.  When the
    !---   g, l, c index changes, the *_pf array will be set correctly
    !--- Do the same for cols setting c_li, c_gi, c_lf, c_gf and
    !---   lunits setting l_gi, l_gf.

    curc = 0
    curl = 0
    curg = 0
    do p = begp,endp
       if (pptr%column(p) /= curc) then
          curc = pptr%column(p)
          if (curc < begc .or. curc > endc) then
             write(6,*) 'clm_ptrs_compdown ERROR: pcolumn ',p,curc,begc,endc
             call endrun()
          endif
          cptr%pfti(curc) = p
       endif
       cptr%pftf(curc) = p
       cptr%npfts(curc) = cptr%pftf(curc) - cptr%pfti(curc) + 1
       if (pptr%landunit(p) /= curl) then
          curl = pptr%landunit(p)
          if (curl < begl .or. curl > endl) then
             write(6,*) 'clm_ptrs_compdown ERROR: plandunit ',p,curl,begl,endl
             call endrun()
          endif
          lptr%pfti(curl) = p
       endif
       lptr%pftf(curl) = p
       lptr%npfts(curl) = lptr%pftf(curl) - lptr%pfti(curl) + 1
       if (pptr%gridcell(p) /= curg) then
          curg = pptr%gridcell(p)
          if (curg < begg .or. curg > endg) then
             write(6,*) 'clm_ptrs_compdown ERROR: pgridcell ',p,curg,begg,endg
             call endrun()
          endif
          gptr%pfti(curg) = p
       endif
       gptr%pftf(curg) = p
       gptr%npfts(curg) = gptr%pftf(curg) - gptr%pfti(curg) + 1
    enddo

    curg = 0
    curl = 0
    do c = begc,endc
       if (cptr%landunit(c) /= curl) then
          curl = cptr%landunit(c)
          if (curl < begl .or. curl > endl) then
             write(6,*) 'clm_ptrs_compdown ERROR: clandunit ',c,curl,begl,endl
             call endrun()
          endif
          lptr%coli(curl) = c
       endif
       lptr%colf(curl) = c
       lptr%ncolumns(curl) = lptr%colf(curl) - lptr%coli(curl) + 1
       if (cptr%gridcell(c) /= curg) then
          curg = cptr%gridcell(c)
          if (curg < begg .or. curg > endg) then
             write(6,*) 'clm_ptrs_compdown ERROR: cgridcell ',c,curg,begg,endg
             call endrun()
          endif
          gptr%coli(curg) = c
       endif
       gptr%colf(curg) = c
       gptr%ncolumns(curg) = gptr%colf(curg) - gptr%coli(curg) + 1
    enddo

    curg = 0
    do l = begl,endl
       if (lptr%gridcell(l) /= curg) then
          curg = lptr%gridcell(l)
          if (curg < begg .or. curg > endg) then
             write(6,*) 'clm_ptrs_compdown ERROR: lgridcell ',l,curg,begg,endg
             call endrun()
          endif
          gptr%luni(curg) = l
       endif
       gptr%lunf(curg) = l
       gptr%nlandunits(curg) = gptr%lunf(curg) - gptr%luni(curg) + 1
    enddo

    end subroutine clm_ptrs_compdown
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_ptrs_check
!
! !INTERFACE:
  subroutine clm_ptrs_check()
!
! !DESCRIPTION:
! Checks and writes out a summary of subgrid data
!
! !USES
    use clmtype, only : clm3, gridcell_type, landunit_type, &
                        column_type, pft_type
    use decompMod , only : get_proc_bounds

! !ARGUMENTS
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 2005.11.15  T Craig Creation
!
!EOP
!
! !LOCAL VARIABLES:
    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
    integer :: begg,endg,begl,endl,begc,endc,begp,endp   ! beg/end indices
    integer :: g,l,c,p       ! loop counters
    logical :: error         ! error flag
!------------------------------------------------------------------------------

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p
    
    if (masterproc) write(6,*) ' '
    if (masterproc) write(6,*) '---clm_ptrs_check:'
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !--- check index ranges ---
    error = .false.
    if (minval(gptr%luni) < begl .or. maxval(gptr%luni) > endl) error=.true.
    if (minval(gptr%lunf) < begl .or. maxval(gptr%lunf) > endl) error=.true.
    if (minval(gptr%coli) < begc .or. maxval(gptr%coli) > endc) error=.true.
    if (minval(gptr%colf) < begc .or. maxval(gptr%colf) > endc) error=.true.
    if (minval(gptr%pfti) < begp .or. maxval(gptr%pfti) > endp) error=.true.
    if (minval(gptr%pftf) < begp .or. maxval(gptr%pftf) > endp) error=.true.
    if (error) then
       write(6,*) '   clm_ptrs_check: g index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(6,*) '   clm_ptrs_check: g index ranges - OK'

    error = .false.
    if (minval(lptr%gridcell) < begg .or. maxval(lptr%gridcell) > endg) error=.true.
    if (minval(lptr%coli) < begc .or. maxval(lptr%coli) > endc) error=.true.
    if (minval(lptr%colf) < begc .or. maxval(lptr%colf) > endc) error=.true.
    if (minval(lptr%pfti) < begp .or. maxval(lptr%pfti) > endp) error=.true.
    if (minval(lptr%pftf) < begp .or. maxval(lptr%pftf) > endp) error=.true.
    if (error) then
       write(6,*) '   clm_ptrs_check: l index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(6,*) '   clm_ptrs_check: l index ranges - OK'

    error = .false.
    if (minval(cptr%gridcell) < begg .or. maxval(cptr%gridcell) > endg) error=.true.
    if (minval(cptr%landunit) < begl .or. maxval(cptr%landunit) > endl) error=.true.
    if (minval(cptr%pfti) < begp .or. maxval(cptr%pfti) > endp) error=.true.
    if (minval(cptr%pftf) < begp .or. maxval(cptr%pftf) > endp) error=.true.
    if (error) then
       write(6,*) '   clm_ptrs_check: c index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(6,*) '   clm_ptrs_check: c index ranges - OK'

    error = .false.
    if (minval(pptr%gridcell) < begg .or. maxval(pptr%gridcell) > endg) error=.true.
    if (minval(pptr%landunit) < begl .or. maxval(pptr%landunit) > endl) error=.true.
    if (minval(pptr%column) < begc .or. maxval(pptr%column) > endc) error=.true.
    if (error) then
       write(6,*) '   clm_ptrs_check: p index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(6,*) '   clm_ptrs_check: p index ranges - OK'

    !--- check that indices in arrays are monotonically increasing ---
    error = .false.
    do g=begg+1,endg
      if (gptr%luni(g) < gptr%luni(g-1)) error = .true.
      if (gptr%lunf(g) < gptr%lunf(g-1)) error = .true.
      if (gptr%coli(g) < gptr%coli(g-1)) error = .true.
      if (gptr%colf(g) < gptr%colf(g-1)) error = .true.
      if (gptr%pfti(g) < gptr%pfti(g-1)) error = .true.
      if (gptr%pftf(g) < gptr%pftf(g-1)) error = .true.
      if (error) then
         write(6,*) '   clm_ptrs_check: g mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(6,*) '   clm_ptrs_check: g mono increasing - OK'

    error = .false.
    do l=begl+1,endl
      if (lptr%gridcell(l) < lptr%gridcell(l-1)) error = .true.
      if (lptr%coli(l) < lptr%coli(l-1)) error = .true.
      if (lptr%colf(l) < lptr%colf(l-1)) error = .true.
      if (lptr%pfti(l) < lptr%pfti(l-1)) error = .true.
      if (lptr%pftf(l) < lptr%pftf(l-1)) error = .true.
      if (error) then
         write(6,*) '   clm_ptrs_check: l mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(6,*) '   clm_ptrs_check: l mono increasing - OK'

    error = .false.
    do c=begc+1,endc
      if (cptr%gridcell(c) < cptr%gridcell(c-1)) error = .true.
      if (cptr%landunit(c) < cptr%landunit(c-1)) error = .true.
      if (cptr%pfti(c) < cptr%pfti(c-1)) error = .true.
      if (cptr%pftf(c) < cptr%pftf(c-1)) error = .true.
      if (error) then
         write(6,*) '   clm_ptrs_check: c mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(6,*) '   clm_ptrs_check: c mono increasing - OK'

    error = .false.
    do p=begp+1,endp
      if (pptr%gridcell(p) < pptr%gridcell(p-1)) error = .true.
      if (pptr%landunit(p) < pptr%landunit(p-1)) error = .true.
      if (pptr%column  (p) < pptr%column  (p-1)) error = .true.
      if (error) then
         write(6,*) '   clm_ptrs_check: p mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(6,*) '   clm_ptrs_check: p mono increasing - OK'

    !--- check that the tree is internally consistent ---
    error = .false.
    do g = begg, endg
       do l = gptr%luni(g),gptr%lunf(g)
          if (lptr%gridcell(l) /= g) error = .true.
          do c = lptr%coli(l),lptr%colf(l)
             if (cptr%gridcell(c) /= g) error = .true.
             if (cptr%landunit(c) /= l) error = .true.
             do p = cptr%pfti(c),cptr%pftf(c)
                if (pptr%gridcell(p) /= g) error = .true.
                if (pptr%landunit(p) /= l) error = .true.
                if (pptr%column(p)   /= c) error = .true.
                if (error) then
                   write(6,*) '   clm_ptrs_check: tree consistent - ERROR'
                   call endrun()
                endif
             enddo
          enddo
       enddo
    enddo
    if (masterproc) write(6,*) '   clm_ptrs_check: tree consistent - OK'
    if (masterproc) write(6,*) ' '

end subroutine clm_ptrs_check
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_veg_compete
!
! !INTERFACE:
!  subroutine set_landunit_veg_compete (ltype, wtxy, vegxy, &
  subroutine set_landunit_veg_compete (ltype, &
                           nw, gi, li, ci, pi, setdata)
!
! !DESCRIPTION: 
! Initialize vegetated landunit with competition
!
! !USES
    use clmtype   , only : clm3, model_type, gridcell_type, landunit_type, &
                           column_type,pft_type
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varcon, only : istsoil
    use clm_varpar, only : numpft, maxpatch_pft
    use clm_varctl, only : allocate_all_vegpfts
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
!    real(r8), intent(in)    :: wtxy(:,:)         ! subgrid patch weights
!    integer , intent(in)    :: vegxy(:,:)        ! PFT types 
    integer , intent(in)    :: nw                ! cell index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: n                                ! loop index
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    integer  :: pitype                           ! pft itype
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft

!------------------------------------------------------------------------

    ! Set decomposition properties

!    call subgrid_get_gcellinfo(nw, wtxy, nveg=npfts, wtveg=wtlunit2gcell)
    call subgrid_get_gcellinfo(nw, nveg=npfts, wtveg=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       lptr => clm3%g%l
       cptr => clm3%g%l%c
       pptr => clm3%g%l%c%p

       ncols = 1
       
       li = li + 1
       ci = ci + 1

       if (setdata) then
          ! Set landunit properties
          lptr%ifspecial(li) = .false.
          lptr%lakpoi(li)    = .false.
          lptr%itype(li)     = ltype
       
          lptr%gridcell (li) = gi
          lptr%wtgcell(li) = wtlunit2gcell

          ! Set column properties for this landunit (only one column on landunit)
          cptr%itype(ci)    = 1
      
          cptr%gridcell (ci) = gi
          cptr%wtgcell(ci) = wtlunit2gcell
          cptr%landunit (ci) = li
          cptr%wtlunit(ci) = 1.0_r8
       endif ! setdata

       ! Set pft properties for this landunit

       if (allocate_all_vegpfts) then
          do n = 1,numpft+1
             pi = pi + 1
             pitype = n-1
             if (setdata) then
                pptr%mxy(pi)      = n
                pptr%itype(pi)    = pitype
                pptr%gridcell (pi) = gi
                pptr%landunit (pi) = li
                pptr%column (pi) = ci
                pptr%wtgcell(pi) = 0.0_r8
                pptr%wtlunit(pi) = 0.0_r8
                pptr%wtcol(pi) = 0.0_r8
                do m = 1,maxpatch_pft
                   if (vegxy(nw,m) == pitype .and. wtxy(nw,m) > 0._r8) then
                      pptr%wtgcell(pi)  = pptr%wtgcell(pi) + wtxy(nw,m)
                      pptr%wtlunit(pi)  = pptr%wtlunit(pi) + wtxy(nw,m) / wtlunit2gcell
                      pptr%wtcol(pi)  = pptr%wtcol(pi) + wtxy(nw,m) / wtlunit2gcell
                   end if
                end do
             endif ! setdata
          end do
       else
          do m = 1,maxpatch_pft
             if (wtxy(nw,m) > 0._r8) then
                pi = pi + 1
                if (setdata) then
                   pptr%mxy(pi)      = m
                   pptr%itype(pi)    = vegxy(nw,m)
                   pptr%gridcell (pi) = gi
                   pptr%wtgcell(pi) = wtxy(nw,m)
                   pptr%landunit (pi) = li
                   pptr%wtlunit(pi) = wtxy(nw,m) / wtlunit2gcell
                   pptr%column (pi) = ci
                   pptr%wtcol(pi) = wtxy(nw,m) / wtlunit2gcell
                endif ! setdata
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
!  subroutine set_landunit_wet_ice_lake (ltype, wtxy, vegxy, &
  subroutine set_landunit_wet_ice_lake (ltype, &
                           nw, gi, li, ci, pi, setdata)
!
! !DESCRIPTION: 
! Initialize wet_ice_lake landunits that are non-urban (lake, wetland, glacier)
!
! !USES
    use clmtype   , only : clm3, model_type, gridcell_type, landunit_type, &
                           column_type,pft_type
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varcon, only : istice, istwet, istdlak
    use clm_varpar, only : npatch_lake, npatch_glacier, npatch_wet
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
!    real(r8), intent(in)    :: wtxy(:,:)         ! subgrid patch weights
!    integer , intent(in)    :: vegxy(:,:)        ! PFT types 
    integer , intent(in)    :: nw                ! cell index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: c                                ! column loop index
    integer  :: ctype                            ! column type
    integer  :: ier                              ! error status 
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    real(r8) :: wtcol2lunit                      ! col weight in landunit
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft
!------------------------------------------------------------------------

    ! Set decomposition properties

    if (ltype == istwet) then
!       call subgrid_get_gcellinfo(nw, wtxy, nwetland=npfts, wtwetland=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nwetland=npfts, wtwetland=wtlunit2gcell)
       m = npatch_wet
    else if (ltype == istdlak) then
!       call subgrid_get_gcellinfo(nw, wtxy, nlake=npfts, wtlake=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nlake=npfts, wtlake=wtlunit2gcell)
       m = npatch_lake
    else if (ltype == istice) then 
!       call subgrid_get_gcellinfo(nw, wtxy, nglacier=npfts, wtglacier=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nglacier=npfts, wtglacier=wtlunit2gcell)
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
          write(6,*)' set_landunit_wet_ice_lake: compete landunit must'// &
                    ' have one column and one pft '
          write(6,*)' current values of ncols, pfts=',ncols,npfts
          call endrun()
       end if

       ncols = 1

       ! Currently assume that each landunit only has only one column 
       ! (of type 1) and that each column has its own pft
       
       wtcol2lunit = 1.0_r8/ncols
       ctype = 1

       li = li + 1
       ci = ci + 1
       pi = pi + 1 

       if (setdata) then
       
          ! Determine landunit properties 

          lptr%itype(li)     = ltype
          lptr%ifspecial(li) = .true.
          if (ltype == istdlak) then
             lptr%lakpoi(li) = .true.
          else
             lptr%lakpoi(li) = .false.
          end if
       
          lptr%gridcell (li) = gi
          lptr%wtgcell(li) = wtlunit2gcell

          ! Determine column and properties
          ! For the wet, ice or lake landunits it is assumed that each 
          ! column has its own pft
       
          cptr%itype(ci)    = ctype
       
          cptr%gridcell (ci) = gi
          cptr%wtgcell(ci) = wtcol2lunit * wtlunit2gcell
          cptr%landunit (ci) = li
          cptr%wtlunit(ci) = wtcol2lunit

          ! Set pft properties

          pptr%mxy(pi)      = m
          pptr%itype(pi)    = vegxy(nw,m)
     
          pptr%gridcell (pi) = gi
          pptr%wtgcell(pi) = wtcol2lunit * wtlunit2gcell
          pptr%landunit (pi) = li
          pptr%wtlunit(pi) = wtcol2lunit
          pptr%column (pi) = ci
          pptr%wtcol(pi) = 1.0_r8
       endif ! setdata
    end if
       
  end subroutine set_landunit_wet_ice_lake

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_crop_noncompete
!
! !INTERFACE:
!  subroutine set_landunit_crop_noncompete (ltype, wtxy, vegxy, &
  subroutine set_landunit_crop_noncompete (ltype, &
                           nw, gi, li, ci, pi, setdata)
!
! !DESCRIPTION: 
! Initialize crop landunit without competition
!
! !USES
    use clmtype   , only : clm3, model_type, gridcell_type, landunit_type, &
                           column_type,pft_type
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varpar, only : npatch_crop, npatch_glacier
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
!    real(r8), intent(in)    :: wtxy(:,:)         ! subgrid patch weights
!    integer , intent(in)    :: vegxy(:,:)        ! PFT types 
    integer , intent(in)    :: nw                ! cell index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft
!------------------------------------------------------------------------

    ! Set decomposition properties

!    call subgrid_get_gcellinfo(nw, wtxy, ncrop=npfts, wtcrop=wtlunit2gcell)
    call subgrid_get_gcellinfo(nw, ncrop=npfts, wtcrop=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       lptr => clm3%g%l
       cptr => clm3%g%l%c
       pptr => clm3%g%l%c%p
       
       ! Set landunit properties - each column has its own pft
       
       ncols = npfts
       
       li = li + 1   

       if (setdata) then
          lptr%itype(li)     = ltype
          lptr%ifspecial(li) = .false.
          lptr%lakpoi(li)    = .false.
       
          lptr%gridcell (li) = gi
          lptr%wtgcell(li) = wtlunit2gcell
       endif ! setdata

       ! Set column and pft properties for this landunit 
       ! (each column has its own pft)

       do m = npatch_glacier+1, npatch_crop
          if (wtxy(nw,m) > 0._r8) then
             ci = ci + 1
             pi = pi + 1
             
             if (setdata) then
                cptr%itype(ci)    = 1
                pptr%itype(pi)    = vegxy(nw,m)
                pptr%mxy(pi)      = m
             
                cptr%gridcell (ci) = gi
                cptr%wtgcell(ci) = wtxy(nw,m)
                cptr%landunit (ci) = li
                cptr%wtlunit(ci) = wtxy(nw,m) / wtlunit2gcell

                pptr%gridcell (pi) = gi
                pptr%wtgcell(pi) = wtxy(nw,m)
                pptr%landunit (pi) = li
                pptr%wtlunit(pi) = wtxy(nw,m) / wtlunit2gcell
                pptr%column (pi) = ci
                pptr%wtcol(pi) = 1.0_r8
             endif ! setdata
          end if
       end do

    end if
       
  end subroutine set_landunit_crop_noncompete

!------------------------------------------------------------------------------

end module initGridCellsMod
