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
  use clm_varsur  , only : topoxy
  use clm_varctl  , only : iulog

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
  private set_landunit_urban
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !PRIVATE DATA MEMBERS: None
!EOP
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
    use clmtype 
    use domainMod   , only : ldomain
    use decompMod   , only : ldecomp, get_proc_global, get_proc_bounds
    use clm_varcon  , only : istsoil, istice, istwet, istdlak, isturb, istice_mec, &
                             udens_tbd, udens_hd, udens_md
    use clm_varctl  , only : create_glacier_mec_landunit
    use clm_varcon  , only : istcrop
    use subgridMod  , only : subgrid_get_gcellinfo
    use shr_const_mod,only : SHR_CONST_PI
    use surfrdMod   , only : crop_prog
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: li,ci,pi,m,na,gdc,gsn,glo    ! indices
    integer :: nveg           ! number of pfts in naturally vegetated landunit
    integer :: ltype          ! landunit type
    real(r8):: wtveg          ! weight (gridcell) of naturally veg landunit
    integer :: ncrop          ! number of crop pfts in crop landunit
    real(r8):: wtcrop         ! weight (gridcell) of crop landunit
    integer :: nlake          ! number of pfts (columns) in lake landunit
    real(r8):: wtlake         ! weight (gridcell) of lake landunit
    integer :: nwetland       ! number of pfts (columns) in wetland landunit
    real(r8):: wtwetland      ! weight (gridcell) of wetland landunit
    integer :: nglacier       ! number of pfts (columns) in glacier landunit
    real(r8):: wtglacier      ! weight (gridcell) of glacier landunit
    integer :: nglacier_mec   ! number of pfts (columns) in glacier landunit
    real(r8):: wtglacier_mec  ! weight (gridcell) of glacier_mec landunit
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

 !------------------------------------------------------------------------

    ! Set pointers into derived types for this module


    ! Get total global number of grid cells, landunits, columns and pfts 
    
    call get_proc_global(numg,numl,numc,nump)
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! For each land gridcell on global grid determine landunit, column and pft properties

    li    = begl-1
    ci    = begc-1
    pi    = begp-1

    if ( crop_prog )then
       ltype = istcrop
    else
       ltype = istsoil
    end if

    !----- Set clm3 variables -----
    do gdc = begg,endg

       glo = ldecomp%gdc2glo(gdc)
       nwtxy = gdc

       my_gcell = .false.
       if (gdc >= begg .and. gdc <= endg) then
          my_gcell = .true.
       endif

       ! Determine naturally vegetated landunit

       call set_landunit_veg_compete(               &
            ltype=istsoil, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine crop landunit

       call set_landunit_crop_noncompete(           &
            ltype=ltype, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine urban tall building district landunit

       call set_landunit_urban( &
!           ltype=isturb, wtxy=wtxy, vegxy=vegxy,   &
            ltype=isturb, udenstype=udens_tbd, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine urban high density landunit

       call set_landunit_urban( &
!           ltype=isturb, wtxy=wtxy, vegxy=vegxy,   &
            ltype=isturb, udenstype=udens_hd, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine urban medium density landunit

       call set_landunit_urban( &
!           ltype=isturb, wtxy=wtxy, vegxy=vegxy,   &
            ltype=isturb, udenstype=udens_md, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine lake, wetland and glacier landunits 

       call set_landunit_wet_ice_lake(              &
            ltype=istdlak, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       call set_landunit_wet_ice_lake(              &
            ltype=istwet, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       call set_landunit_wet_ice_lake(              &
            ltype=istice, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       if (create_glacier_mec_landunit) then
          call set_landunit_wet_ice_lake(              &
               ltype=istice_mec, &
               nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell, &
               glcmask = ldomain%glcmask(gdc))
       endif

       ! Make ice sheet masks

       grc%gris_mask(gdc) = 0._r8
       grc%gris_area(gdc) = 0._r8
       grc%aais_mask(gdc) = 0._r8
       grc%aais_area(gdc) = 0._r8
      
       ! Greenland mask
       if ( (ldomain%latc(gdc) >  58. .and. ldomain%latc(gdc) <= 67.  .and.   &
             ldomain%lonc(gdc) > 302. .and. ldomain%lonc(gdc) < 330.)         &
                                      .or.                                 &
            (ldomain%latc(gdc) >  67. .and. ldomain%latc(gdc) <= 70. .and.    &
             ldomain%lonc(gdc) > 300. .and. ldomain%lonc(gdc) < 345.)         &
                                      .or.                                 &
            (ldomain%latc(gdc) >  70. .and. ldomain%latc(gdc) <= 75. .and.    &
             ldomain%lonc(gdc) > 295. .and. ldomain%lonc(gdc) < 350.)         &
                                      .or.                                 &
            (ldomain%latc(gdc) >  75. .and. ldomain%latc(gdc) <= 79. .and.    &
             ldomain%lonc(gdc) > 285. .and. ldomain%lonc(gdc) < 350.)         &
                                      .or.                                 &
            (ldomain%latc(gdc) >  79. .and. ldomain%latc(gdc) <  85. .and.    &
             ldomain%lonc(gdc) > 290. .and. ldomain%lonc(gdc) < 355.) ) then
 
            grc%gris_mask(gdc) = 1.0_r8

      elseif (ldomain%latc(gdc) < -60.) then

            grc%aais_mask(gdc) = 1.0_r8

       endif  ! Greenland or Antarctic grid cell

       ! Set clm3 lats/lons

       if (my_gcell) then
          grc%gindex(gdc) = glo
          grc%latdeg(gdc) = ldomain%latc(gdc) 
          grc%londeg(gdc) = ldomain%lonc(gdc) 
          grc%lat(gdc)    = grc%latdeg(gdc) * SHR_CONST_PI/180._r8  
          grc%lon(gdc)    = grc%londeg(gdc) * SHR_CONST_PI/180._r8
          grc%area(gdc)   = ldomain%area(gdc)
       endif

    enddo

    ! Fill in subgrid datatypes

    call clm_ptrs_compdown()
    call clm_ptrs_check()

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
    use clmtype
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
!
! !LOCAL VARIABLES:
    integer :: begg,endg,begl,endl,begc,endc,begp,endp ! beg/end glcp
    integer :: g,l,c,p               ! loop counters
    integer :: curg,curl,curc,curp   ! tracks g,l,c,p indexes in arrays
!EOP
!------------------------------------------------------------------------------


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
       if (pft%column(p) /= curc) then
          curc = pft%column(p)
          if (curc < begc .or. curc > endc) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: pcolumn ',p,curc,begc,endc
             call endrun()
          endif
          col%pfti(curc) = p
       endif
       col%pftf(curc) = p
       col%npfts(curc) = col%pftf(curc) - col%pfti(curc) + 1
       if (pft%landunit(p) /= curl) then
          curl = pft%landunit(p)
          if (curl < begl .or. curl > endl) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: plandunit ',p,curl,begl,endl
             call endrun()
          endif
          lun%pfti(curl) = p
       endif
       lun%pftf(curl) = p
       lun%npfts(curl) = lun%pftf(curl) - lun%pfti(curl) + 1
       if (pft%gridcell(p) /= curg) then
          curg = pft%gridcell(p)
          if (curg < begg .or. curg > endg) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: pgridcell ',p,curg,begg,endg
             call endrun()
          endif
          grc%pfti(curg) = p
       endif
       grc%pftf(curg) = p
       grc%npfts(curg) = grc%pftf(curg) - grc%pfti(curg) + 1
    enddo

    curg = 0
    curl = 0
    do c = begc,endc
       if (col%landunit(c) /= curl) then
          curl = col%landunit(c)
          if (curl < begl .or. curl > endl) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: clandunit ',c,curl,begl,endl
             call endrun()
          endif
          lun%coli(curl) = c
       endif
       lun%colf(curl) = c
       lun%ncolumns(curl) = lun%colf(curl) - lun%coli(curl) + 1
       if (col%gridcell(c) /= curg) then
          curg = col%gridcell(c)
          if (curg < begg .or. curg > endg) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: cgridcell ',c,curg,begg,endg
             call endrun()
          endif
          grc%coli(curg) = c
       endif
       grc%colf(curg) = c
       grc%ncolumns(curg) = grc%colf(curg) - grc%coli(curg) + 1
    enddo

    curg = 0
    do l = begl,endl
       if (lun%gridcell(l) /= curg) then
          curg = lun%gridcell(l)
          if (curg < begg .or. curg > endg) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: lgridcell ',l,curg,begg,endg
             call endrun()
          endif
          grc%luni(curg) = l
       endif
       grc%lunf(curg) = l
       grc%nlandunits(curg) = grc%lunf(curg) - grc%luni(curg) + 1
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
    use clmtype
    use decompMod , only : get_proc_bounds

! !ARGUMENTS
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 2005.11.15  T Craig Creation
!
!
! !LOCAL VARIABLES:
    integer :: begg,endg,begl,endl,begc,endc,begp,endp   ! beg/end indices
    integer :: g,l,c,p       ! loop counters
    logical :: error         ! error flag
!EOP
!------------------------------------------------------------------------------

    
    if (masterproc) write(iulog,*) ' '
    if (masterproc) write(iulog,*) '---clm_ptrs_check:'
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !--- check index ranges ---
    error = .false.
    if (minval(grc%luni) < begl .or. maxval(grc%luni) > endl) error=.true.
    if (minval(grc%lunf) < begl .or. maxval(grc%lunf) > endl) error=.true.
    if (minval(grc%coli) < begc .or. maxval(grc%coli) > endc) error=.true.
    if (minval(grc%colf) < begc .or. maxval(grc%colf) > endc) error=.true.
    if (minval(grc%pfti) < begp .or. maxval(grc%pfti) > endp) error=.true.
    if (minval(grc%pftf) < begp .or. maxval(grc%pftf) > endp) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: g index ranges - ERROR'
       write(iulog,*)'minval,beg,maxval,end'
       write(iulog,*) minval(grc%luni),begl,maxval(grc%luni),endl
       write(iulog,*) minval(grc%lunf),begl,maxval(grc%lunf),endl
       write(iulog,*) minval(grc%coli),begc,maxval(grc%coli),endc
       write(iulog,*) minval(grc%colf),begc,maxval(grc%colf),endc
       write(iulog,*) minval(grc%pfti),begp,maxval(grc%pfti),endp
       write(iulog,*) minval(grc%pftf),begp,maxval(grc%pftf),endp
       call endrun()
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: g index ranges - OK'

    error = .false.
    if (minval(lun%gridcell) < begg .or. maxval(lun%gridcell) > endg) error=.true.
    if (minval(lun%coli) < begc .or. maxval(lun%coli) > endc) error=.true.
    if (minval(lun%colf) < begc .or. maxval(lun%colf) > endc) error=.true.
    if (minval(lun%pfti) < begp .or. maxval(lun%pfti) > endp) error=.true.
    if (minval(lun%pftf) < begp .or. maxval(lun%pftf) > endp) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: l index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: l index ranges - OK'

    error = .false.
    if (minval(col%gridcell) < begg .or. maxval(col%gridcell) > endg) error=.true.
    if (minval(col%landunit) < begl .or. maxval(col%landunit) > endl) error=.true.
    if (minval(col%pfti) < begp .or. maxval(col%pfti) > endp) error=.true.
    if (minval(col%pftf) < begp .or. maxval(col%pftf) > endp) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: c index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: c index ranges - OK'

    error = .false.
    if (minval(pft%gridcell) < begg .or. maxval(pft%gridcell) > endg) error=.true.
    if (minval(pft%landunit) < begl .or. maxval(pft%landunit) > endl) error=.true.
    if (minval(pft%column) < begc .or. maxval(pft%column) > endc) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: p index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: p index ranges - OK'

    !--- check that indices in arrays are monotonically increasing ---
    error = .false.
    do g=begg+1,endg
      if (grc%luni(g) < grc%luni(g-1)) error = .true.
      if (grc%lunf(g) < grc%lunf(g-1)) error = .true.
      if (grc%coli(g) < grc%coli(g-1)) error = .true.
      if (grc%colf(g) < grc%colf(g-1)) error = .true.
      if (grc%pfti(g) < grc%pfti(g-1)) error = .true.
      if (grc%pftf(g) < grc%pftf(g-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: g mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: g mono increasing - OK'

    error = .false.
    do l=begl+1,endl
      if (lun%gridcell(l) < lun%gridcell(l-1)) error = .true.
      if (lun%coli(l) < lun%coli(l-1)) error = .true.
      if (lun%colf(l) < lun%colf(l-1)) error = .true.
      if (lun%pfti(l) < lun%pfti(l-1)) error = .true.
      if (lun%pftf(l) < lun%pftf(l-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: l mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: l mono increasing - OK'

    error = .false.
    do c=begc+1,endc
      if (col%gridcell(c) < col%gridcell(c-1)) error = .true.
      if (col%landunit(c) < col%landunit(c-1)) error = .true.
      if (col%pfti(c) < col%pfti(c-1)) error = .true.
      if (col%pftf(c) < col%pftf(c-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: c mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: c mono increasing - OK'

    error = .false.
    do p=begp+1,endp
      if (pft%gridcell(p) < pft%gridcell(p-1)) error = .true.
      if (pft%landunit(p) < pft%landunit(p-1)) error = .true.
      if (pft%column  (p) < pft%column  (p-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: p mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: p mono increasing - OK'

    !--- check that the tree is internally consistent ---
    error = .false.
    do g = begg, endg
       do l = grc%luni(g),grc%lunf(g)
          if (lun%gridcell(l) /= g) error = .true.
          do c = lun%coli(l),lun%colf(l)
             if (col%gridcell(c) /= g) error = .true.
             if (col%landunit(c) /= l) error = .true.
             do p = col%pfti(c),col%pftf(c)
                if (pft%gridcell(p) /= g) error = .true.
                if (pft%landunit(p) /= l) error = .true.
                if (pft%column(p)   /= c) error = .true.
                if (error) then
                   write(iulog,*) '   clm_ptrs_check: tree consistent - ERROR'
                   call endrun()
                endif
             enddo
          enddo
       enddo
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: tree consistent - OK'
    if (masterproc) write(iulog,*) ' '

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
    use clmtype 
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varpar, only : numpft, maxpatch_pft, numcft
    use clm_varctl, only : allocate_all_vegpfts, create_crop_landunit
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
! Created by ?
! 2005.11.25 Updated by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: n                                ! loop index
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    integer  :: pitype                           ! pft itype
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell

!------------------------------------------------------------------------

    ! Set decomposition properties

!    call subgrid_get_gcellinfo(nw, wtxy, nveg=npfts, wtveg=wtlunit2gcell)
    call subgrid_get_gcellinfo(nw, nveg=npfts, wtveg=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module


       ncols = 1
       
       li = li + 1
       ci = ci + 1

       if (setdata) then
          ! Set landunit properties
          lun%ifspecial(li) = .false.
          lun%lakpoi(li)    = .false.
          lun%urbpoi(li)    = .false.
          lun%itype(li)     = ltype
       
          lun%gridcell (li) = gi
          lun%wtgcell(li) = wtlunit2gcell

          ! Set column properties for this landunit (only one column on landunit)
          col%itype(ci)    = 1
      
          col%gridcell (ci) = gi
          col%wtgcell(ci) = wtlunit2gcell
          col%landunit (ci) = li
          col%wtlunit(ci) = 1.0_r8
       endif ! setdata

       ! Set pft properties for this landunit

       if (create_crop_landunit) then
          do n = 1,numpft+1-numcft
             pi = pi + 1
             pitype = n-1
             if (setdata) then
                pft%mxy(pi)      = n
                pft%itype(pi)    = pitype
                pft%gridcell(pi) = gi
                pft%landunit(pi) = li
                pft%column (pi) = ci

                if (wtlunit2gcell > 0._r8) then
                   pft%wtgcell(pi) = 0.0_r8
                   pft%wtlunit(pi) = 0.0_r8
                   pft%wtcol(pi) = 0.0_r8
                   do m = 1,maxpatch_pft
                      if (vegxy(nw,m) == pitype) then
                         pft%wtgcell(pi)  = pft%wtgcell(pi) + wtxy(nw,m)
                         pft%wtlunit(pi)  = pft%wtlunit(pi) + wtxy(nw,m) / wtlunit2gcell
                         pft%wtcol(pi)  = pft%wtcol(pi) + wtxy(nw,m) / wtlunit2gcell
                      end if
                   end do
                else  ! wtlunit2gcell == 0._r8
                   ! TODO WJS: Temporarily setting this to equal weighting for all
                   ! pfts. In the future, we could potentially get some info about this
                   ! from the surface dataset, if it is changed to give pct_pft as % of
                   ! the pft on the landunit
                   pft%wtgcell(pi) = 0._r8
                   pft%wtlunit(pi) = 1._r8 / (numpft+1-numcft)
                   pft%wtcol(pi)   = 1._r8 / (numpft+1-numcft)
                end if
             endif ! setdata
          end do
       else if (allocate_all_vegpfts) then
          do n = 1,numpft+1
             pi = pi + 1
             pitype = n-1
             if (setdata) then
                pft%mxy(pi)      = n
                pft%itype(pi)    = pitype
                pft%gridcell(pi) = gi
                pft%landunit(pi) = li
                pft%column (pi) = ci

                if (wtlunit2gcell > 0._r8) then
                   pft%wtgcell(pi) = 0.0_r8
                   pft%wtlunit(pi) = 0.0_r8
                   pft%wtcol(pi) = 0.0_r8
                   do m = 1,maxpatch_pft
                      if (vegxy(nw,m) == pitype) then
                         pft%wtgcell(pi)  = pft%wtgcell(pi) + wtxy(nw,m)
                         pft%wtlunit(pi)  = pft%wtlunit(pi) + wtxy(nw,m) / wtlunit2gcell
                         pft%wtcol(pi)  = pft%wtcol(pi) + wtxy(nw,m) / wtlunit2gcell
                      end if
                   end do
                else  ! wtlunit2gcell == 0._r8
                   ! TODO WJS: Temporarily setting this to equal weighting for all
                   ! pfts. In the future, we could potentially get some info about this
                   ! from the surface dataset, if it is changed to give pct_pft as % of
                   ! the pft on the landunit
                   pft%wtgcell(pi) = 0._r8
                   pft%wtlunit(pi) = 1._r8 / (numpft+1)
                   pft%wtcol(pi)   = 1._r8 / (numpft+1)
                end if
             endif ! setdata
          end do
       else
          write(iulog,*) 'allocate_all_vegpfts=false is no longer supported'
          call endrun()
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
                           nw, gi, li, ci, pi, setdata, glcmask)
!
! !DESCRIPTION: 
! Initialize wet_ice_lake landunits that are non-urban (lake, wetland, glacier, glacier_mec)
!
! !USES
    use clmtype 
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varcon, only : istice, istwet, istdlak, istice_mec
    use clm_varpar, only : npatch_lake, npatch_glacier, npatch_wet
    use clm_varpar, only : npatch_glacier_mec, maxpatch_glcmec

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
    integer , intent(in), optional :: glcmask    ! = 1 where glc requires sfc mass balance
                                                 ! = 0 otherwise
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: c                                ! column loop index
    integer  :: ctype                            ! column type
    integer  :: ier                              ! error status 
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    real(r8) :: wtcol2lunit                      ! col weight in landunit

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
    else if (ltype == istice_mec) then
!       call subgrid_get_gcellinfo(nw, wtxy, nglacier_mec=npfts, wtglacier_mec=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nglacier_mec=npfts, wtglacier_mec=wtlunit2gcell, &
                                  glcmask = glcmask)
       ! NOTE: multiple columns per landunit, so m is not set here

    else
       write(iulog,*)' set_landunit_wet_ice_lake: ltype of ',ltype,' not valid'
       write(iulog,*)' only istwet, istdlak, istice and istice_mec ltypes are valid'
       call endrun()
    end if

    if (npfts > 0) then

       ! Set pointers into derived types for this module


       if (npfts /=1 .and. ltype /= istice_mec) then
          write(iulog,*)' set_landunit_wet_ice_lake: compete landunit must'// &
               ' have one column and one pft '
          write(iulog,*)' current values of ncols, pfts=',ncols,npfts
          call endrun()
       end if

       if (ltype==istice_mec) then   ! multiple columns per landunit

          ! Assume that columns are of type 1 and that each column has its own pft

          ctype = 1
          li = li + 1

          if (setdata) then

             ! Determine landunit properties

             lun%itype    (li) = ltype
             lun%ifspecial(li) = .true.
             lun%glcmecpoi(li) = .true.
             lun%lakpoi   (li) = .false.
             lun%urbpoi   (li) = .false.
             lun%gridcell (li) = gi
             lun%wtgcell  (li) = wtlunit2gcell

             ! Determine column and properties
             ! (Each column has its own pft)
             ! 
             ! For grid cells with glcmask = 1, make sure all the elevations classes
             !  are populated, even if some have zero fractional area.  This ensures that the 
             !  ice sheet component, glc, will receive a surface mass balance in each elevation 
             !  class wherever the SMB is needed.
             ! Columns with zero weight are referred to as "virtual" columns.
 
             do m = npatch_glacier+1, npatch_glacier_mec

                if (wtxy(nw,m) > 0._r8 .or. glcmask == 1) then

                   ci = ci + 1
                   pi = pi + 1
                   if (wtlunit2gcell > 0._r8) then
                      wtcol2lunit = wtxy(nw,m)/wtlunit2gcell
                   else   ! virtual landunit
                      ! TODO WJS: Temporarily setting this to equal weighting for all
                      ! columns on the landunit. In the future, we could potentially get
                      ! some info about this from the surface dataset, if it is changed
                      ! to give pct_glc_mec as % of the column on the landunit
                      wtcol2lunit = 1._r8 / maxpatch_glcmec
                   endif

                   col%itype    (ci) = ctype
                   col%gridcell (ci) = gi
                   col%wtgcell  (ci) = wtcol2lunit * wtlunit2gcell
                   col%landunit (ci) = li
                   col%wtlunit  (ci) = wtcol2lunit

                   ! Set sfc elevation too

                   cps%glc_topo(ci) = topoxy(nw,m)

                   ! Set pft properties

                   pft%mxy      (pi) = m
                   pft%itype    (pi) = vegxy(nw,m)
                   pft%gridcell (pi) = gi
                   pft%wtgcell  (pi) = wtcol2lunit * wtlunit2gcell
                   pft%landunit (pi) = li
                   pft%wtlunit  (pi) = wtcol2lunit
                   pft%column   (pi) = ci
                   pft%wtcol    (pi) = 1.0_r8

                endif   ! wtxy > 0 or glcmask = 1
             enddo      ! loop over columns
          endif         ! setdata

       else

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

             lun%itype    (li) = ltype
             lun%ifspecial(li) = .true.
             lun%urbpoi   (li) = .false.
             if (ltype == istdlak) then
                lun%lakpoi(li) = .true.
             else
                lun%lakpoi(li) = .false.
             end if
       
             lun%gridcell (li) = gi
             lun%wtgcell(li) = wtlunit2gcell

             ! Determine column and properties
             ! For the wet, ice or lake landunits it is assumed that each 
             ! column has its own pft
       
             col%itype(ci)    = ctype
       
             col%gridcell (ci) = gi
             col%wtgcell(ci) = wtcol2lunit * wtlunit2gcell
             col%landunit (ci) = li
             col%wtlunit(ci) = wtcol2lunit

             ! Set pft properties

             pft%mxy(pi)      = m
             pft%itype(pi)    = vegxy(nw,m)
     
             pft%gridcell (pi) = gi
             pft%wtgcell(pi) = wtcol2lunit * wtlunit2gcell
             pft%landunit (pi) = li
             pft%wtlunit(pi) = wtcol2lunit
             pft%column (pi) = ci
             pft%wtcol(pi) = 1.0_r8
          endif ! setdata
       end if   ! ltype = istice_mec
    endif       ! npfts > 0       

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
    use clmtype
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varctl, only : create_crop_landunit
    use clm_varpar, only : maxpatch_pft, numcft, npatch_glacier_mec
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
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
!------------------------------------------------------------------------

    ! Set decomposition properties

!    call subgrid_get_gcellinfo(nw, wtxy, ncrop=npfts, wtcrop=wtlunit2gcell)
    call subgrid_get_gcellinfo(nw, ncrop=npfts, wtcrop=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       
       ! Set landunit properties - each column has its own pft
       
       ncols = npfts
       
       li = li + 1   

       if (setdata) then
          lun%itype(li)     = ltype
          lun%ifspecial(li) = .false.
          lun%lakpoi(li)    = .false.
          lun%urbpoi(li)    = .false.
          lun%gridcell (li) = gi
          lun%wtgcell(li) = wtlunit2gcell
       endif ! setdata

       ! Set column and pft properties for this landunit 
       ! (each column has its own pft)

       if (create_crop_landunit) then
          do m = maxpatch_pft-numcft+1, maxpatch_pft
             ci = ci + 1
             pi = pi + 1
             
             if (setdata) then
                col%itype(ci)    = 1
                pft%itype(pi)    = m - 1
                pft%mxy(pi)      = m
          
                col%gridcell (ci) = gi
                col%wtgcell(ci) = wtxy(nw,m)
                col%landunit (ci) = li

                pft%gridcell (pi) = gi
                pft%wtgcell(pi) = wtxy(nw,m)
                pft%landunit (pi) = li
                pft%column (pi) = ci
                pft%wtcol(pi) = 1._r8
                if (wtlunit2gcell > 0) then
                   col%wtlunit(ci) = wtxy(nw,m) / wtlunit2gcell
                   pft%wtlunit(pi) = wtxy(nw,m) / wtlunit2gcell
                else
                   ! TODO WJS: Temporarily setting this to equal weighting for all crop
                   ! pfts. In the future, we could potentially get some info about this
                   ! from the surface dataset, if it is changed to give pct_cft as % of
                   ! the cft on the landunit
                   col%wtlunit(ci) = 1._r8 / numcft
                   pft%wtlunit(pi) = 1._r8 / numcft
                end if
             endif ! setdata
          end do
       end if

    end if
       
  end subroutine set_landunit_crop_noncompete

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_urban
!
! !INTERFACE:
!  subroutine set_landunit_urban (ltype, wtxy, vegxy, &
  subroutine set_landunit_urban (ltype, udenstype, &
                                 nw, gi, li, ci, pi, setdata)
!
! !DESCRIPTION: 
! Initialize urban landunits
!
! !USES
    use clm_varcon   , only : isturb, icol_roof, icol_sunwall, icol_shadewall, &
                              icol_road_perv, icol_road_imperv, &
                              udens_tbd, udens_hd, udens_md, udens_base
    use clm_varpar   , only : npatch_urban_tbd, npatch_urban_hd, npatch_urban_md, maxpatch_urb
    use clmtype
    use subgridMod   , only : subgrid_get_gcellinfo
    use UrbanInputMod, only : urbinp
    use decompMod    , only : ldecomp
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: udenstype         ! urban density type
!   real(r8), intent(in)    :: wtxy(:,:)         ! subgrid patch weights
!   integer , intent(in)    :: vegxy(:,:)        ! PFT types 
    integer , intent(in)    :: nw                ! cell index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: c             ! column loop index
    integer  :: m             ! m index in wtxy(nw,m)
    integer  :: n             ! urban density type index
    integer  :: ctype         ! column type
    integer  :: npfts         ! number of pfts in landunit
    integer  :: ncols         ! number of columns in landunit
    integer  :: npatch        ! npatch for the given urban density class
    real(r8) :: wtlunit2gcell ! weight relative to gridcell of landunit
    real(r8) :: wtcol2lunit   ! weight of column with respect to landunit
    real(r8) :: wtlunit_roof  ! weight of roof with respect to landunit
    real(r8) :: wtroad_perv   ! weight of pervious road column with respect to total road
    integer  :: ier           ! error status 
!------------------------------------------------------------------------

    ! Set decomposition properties, and set variables specific to urban density type

    select case (udenstype)
    case (udens_tbd)
       !   call subgrid_get_gcellinfo(nw, wtxy, nurban_tbd=npfts, wturban_tbd=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nurban_tbd=npfts, wturban_tbd=wtlunit2gcell)
       npatch = npatch_urban_tbd
    case (udens_hd)
       !   call subgrid_get_gcellinfo(nw, wtxy, nurban_hd=npfts, wturban_hd=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nurban_hd=npfts, wturban_hd=wtlunit2gcell)
       npatch = npatch_urban_hd
    case (udens_md)
       !   call subgrid_get_gcellinfo(nw, wtxy, nurban_md=npfts, wturban_md=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nurban_md=npfts, wturban_md=wtlunit2gcell)
       npatch = npatch_urban_md
    case default
       write(iulog,*)' set_landunit_urban: unknown udenstype: ', udenstype
       call endrun()
    end select

    n = udenstype - udens_base

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       
       ! Determine landunit properties - each columns has its own pft
       
       ncols = npfts

       li = li + 1
       if (setdata) then
          lun%itype    (li) = ltype
          lun%udenstype(li) = udenstype
          lun%ifspecial(li) = .true.
          lun%lakpoi   (li) = .false.
          lun%urbpoi   (li) = .true.

          lun%gridcell (li) = gi
          lun%wtgcell  (li) = wtlunit2gcell
       endif

       ! Loop through columns for this landunit and set the column and pft properties
       ! For the urban landunits it is assumed that each column has its own pft
       
       do m = npatch, npatch + maxpatch_urb - 1
          if (wtxy(nw,m) > 0._r8) then
                
             wtlunit_roof = urbinp%wtlunit_roof(nw,n)
             wtroad_perv  = urbinp%wtroad_perv(nw,n)
             
             if (m == npatch  ) then
                ctype = icol_roof
                wtcol2lunit = wtlunit_roof
             else if (m == npatch+1) then
                ctype = icol_sunwall
                wtcol2lunit = (1. - wtlunit_roof)/3
             else if (m == npatch+2) then
                ctype = icol_shadewall
                wtcol2lunit = (1. - wtlunit_roof)/3
             else if (m == npatch+3) then
                ctype = icol_road_imperv
                wtcol2lunit = ((1. - wtlunit_roof)/3) * (1.-wtroad_perv)
             else if (m == npatch+4) then
                ctype = icol_road_perv
                wtcol2lunit = ((1. - wtlunit_roof)/3) * (wtroad_perv)
             end if
             
             ci = ci + 1
             pi = pi + 1 
             
             if (setdata) then
                col%itype(ci)     = ctype

                col%gridcell (ci) = gi
                col%wtgcell  (ci) = wtcol2lunit * wtlunit2gcell
                col%landunit (ci) = li
                col%wtlunit  (ci) = wtcol2lunit

                pft%mxy     (pi)  = m
                pft%itype   (pi)  = vegxy(nw,m)
                
                pft%gridcell(pi)  = gi
                pft%wtgcell (pi)  = wtcol2lunit * wtlunit2gcell
                pft%landunit(pi)  = li
                pft%wtlunit (pi)  = wtcol2lunit
                pft%column  (pi)  = ci
                pft%wtcol   (pi)  = 1.0_r8
             end if
             
          end if
       end do   ! end of loop through urban columns-pfts
    end if

  end subroutine set_landunit_urban

!------------------------------------------------------------------------------

end module initGridCellsMod
