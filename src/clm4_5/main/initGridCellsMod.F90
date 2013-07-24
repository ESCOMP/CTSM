module initGridCellsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initializes sub-grid mapping for each land grid cell
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc,iam,mpicom
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  use decompMod   , only : bounds_type, ldecomp
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
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine initGridcells (bounds) 
    !
    ! !DESCRIPTION: 
    ! Initialize sub-grid mapping and allocates space for derived type hierarchy.
    ! For each land gridcell determine landunit, column and pft properties.
    !
    ! !USES
    use clmtype 
    use domainMod   , only : ldomain
    use clm_varcon  , only : istsoil, istice, istwet, istdlak, istice_mec, &
                             isturb_tbd, isturb_hd, isturb_md, istcrop
    use clm_varctl  , only : create_glacier_mec_landunit
    use shr_const_mod,only : SHR_CONST_PI
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: li,ci,pi,m,gdc    ! indices
    !------------------------------------------------------------------------

    ! For each land gridcell on global grid determine landunit, column and pft properties

    li = bounds%begl-1
    ci = bounds%begc-1
    pi = bounds%begp-1

    ! Determine naturally vegetated landunit
    do gdc = bounds%begg,bounds%endg
       call set_landunit_veg_compete(               &
            ltype=istsoil, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
    end do

    ! Determine crop landunit
    do gdc = bounds%begg,bounds%endg
       call set_landunit_crop_noncompete(           &
            ltype=istcrop, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
    end do

    ! Determine urban tall building district landunit
    do gdc = bounds%begg,bounds%endg
       call set_landunit_urban( &
            ltype=isturb_tbd, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)

    end do

    ! Determine urban high density landunit
    do gdc = bounds%begg,bounds%endg
       call set_landunit_urban( &
            ltype=isturb_hd, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
    end do

    ! Determine urban medium density landunit
    do gdc = bounds%begg,bounds%endg
       call set_landunit_urban( &
            ltype=isturb_md, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
    end do

    ! Determine lake, wetland and glacier landunits 
    do gdc = bounds%begg,bounds%endg
       call set_landunit_wet_ice_lake(              &
            ltype=istdlak, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
    end do

    do gdc = bounds%begg,bounds%endg
       call set_landunit_wet_ice_lake(              &
            ltype=istwet, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
    end do

    do gdc = bounds%begg,bounds%endg
       call set_landunit_wet_ice_lake(              &
            ltype=istice, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
    end do

    if (create_glacier_mec_landunit) then
       do gdc = bounds%begg,bounds%endg
          call set_landunit_wet_ice_lake(              &
               ltype=istice_mec, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true., &
               glcmask = ldomain%glcmask(gdc))
       end do
    endif

    ! Set some other gridcell-level variables

    do gdc = bounds%begg,bounds%endg

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

       grc%gindex(gdc) = ldecomp%gdc2glo(gdc)
       grc%latdeg(gdc) = ldomain%latc(gdc) 
       grc%londeg(gdc) = ldomain%lonc(gdc) 
       grc%lat(gdc)    = grc%latdeg(gdc) * SHR_CONST_PI/180._r8  
       grc%lon(gdc)    = grc%londeg(gdc) * SHR_CONST_PI/180._r8
       grc%area(gdc)   = ldomain%area(gdc)

    enddo

    ! Fill in subgrid datatypes

    call clm_ptrs_compdown(bounds)

    call clm_ptrs_check()

  end subroutine initGridcells

  !------------------------------------------------------------------------------
  subroutine clm_ptrs_compdown(bounds)
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
    !
    ! !ARGUMENTS
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c,p               ! loop counters
    integer :: curl,curc,curp      ! tracks l,c,p indexes in arrays
    !------------------------------------------------------------------------------

    !--- Set the current c,l (curc, curl) to zero for initialization,
    !---   these indices track the current "up" index.
    !--- Take advantage of locality of l/c/p cells
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
    do p = bounds%begp,bounds%endp
       if (pft%column(p) /= curc) then
          curc = pft%column(p)
          if (curc < bounds%begc .or. curc > bounds%endc) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: pcolumn ',p,curc,bounds%begc,bounds%endc
             call endrun()
          endif
          col%pfti(curc) = p
       endif
       col%pftf(curc) = p
       col%npfts(curc) = col%pftf(curc) - col%pfti(curc) + 1
       if (pft%landunit(p) /= curl) then
          curl = pft%landunit(p)
          if (curl < bounds%begl .or. curl > bounds%endl) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: plandunit ',p,curl,bounds%begl,bounds%endl
             call endrun()
          endif
          lun%pfti(curl) = p
       endif
       lun%pftf(curl) = p
       lun%npfts(curl) = lun%pftf(curl) - lun%pfti(curl) + 1
    enddo

    curl = 0
    do c = bounds%begc,bounds%endc
       if (col%landunit(c) /= curl) then
          curl = col%landunit(c)
          if (curl < bounds%begl .or. curl > bounds%endl) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: clandunit ',c,curl,bounds%begl,bounds%endl
             call endrun()
          endif
          lun%coli(curl) = c
       endif
       lun%colf(curl) = c
       lun%ncolumns(curl) = lun%colf(curl) - lun%coli(curl) + 1
    enddo

  end subroutine clm_ptrs_compdown

  !------------------------------------------------------------------------------
  subroutine clm_ptrs_check()
    !
    ! !DESCRIPTION:
    ! Checks and writes out a summary of subgrid data
    !
    ! !USES
    use clmtype
    use decompMod , only : get_proc_bounds
    !
    ! !ARGUMENTS
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg,begl,endl,begc,endc,begp,endp   ! beg/end indices
    integer :: g,l,c,p       ! loop counters
    integer :: l_prev        ! l value of previous point
    logical :: error         ! error flag
    !------------------------------------------------------------------------------

    
    if (masterproc) write(iulog,*) ' '
    if (masterproc) write(iulog,*) '---clm_ptrs_check:'
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !--- check index ranges ---
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
    do l=begl+1,endl
      if ((lun%itype(l) == lun%itype(l-1)) .and. &
           lun%gridcell(l) < lun%gridcell(l-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
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
      l = col%landunit(c)
      l_prev = col%landunit(c-1)
      if ((lun%itype(l) == lun%itype(l_prev)) .and. &
           col%gridcell(c) < col%gridcell(c-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
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
      l = pft%landunit(p)
      l_prev = pft%landunit(p-1)
      if ((lun%itype(l) == lun%itype(l_prev)) .and. &
           pft%gridcell(p) < pft%gridcell(p-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
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
    do l = begl, endl
       g = lun%gridcell(l)
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
    if (masterproc) write(iulog,*) '   clm_ptrs_check: tree consistent - OK'
    if (masterproc) write(iulog,*) ' '

  end subroutine clm_ptrs_check

  !------------------------------------------------------------------------
  subroutine set_landunit_veg_compete (ltype, gi, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize vegetated landunit with competition
    !
    ! !USES
    use clmtype
    use clm_varsur, only : wt_lunit, wt_nat_pft
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varpar, only : numpft, maxpatch_pft, numcft, natpft_lb, natpft_ub
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: m                                ! index
    integer  :: lb_offset                        ! offset between natpft_lb and 1
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: pitype                           ! pft itype
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    !------------------------------------------------------------------------

    ! Set decomposition properties

    call subgrid_get_gcellinfo(gi, nveg=npfts)
    wtlunit2gcell = wt_lunit(gi, ltype)

    if (npfts > 0) then
       li = li + 1
       ci = ci + 1

       ! Set landunit properties
       lun%ifspecial(li) = .false.
       lun%lakpoi(li)    = .false.
       lun%urbpoi(li)    = .false.
       lun%itype(li)     = ltype
       lun%gridcell(li)  = gi
       lun%wtgcell(li)   = wtlunit2gcell

       ! Set column properties for this landunit (only one column on landunit)
       col%itype(ci)    = 1
       col%gridcell(ci) = gi
       col%wtgcell(ci)  = wtlunit2gcell
       col%landunit(ci) = li
       col%wtlunit(ci)  = 1.0_r8

       ! Set pft properties for this landunit
       lb_offset = 1 - natpft_lb
       do m = natpft_lb,natpft_ub
          pi               = pi + 1
          pitype           = m
          pft%mxy(pi)      = m + lb_offset
          pft%itype(pi)    = pitype
          pft%gridcell(pi) = gi
          pft%landunit(pi) = li
          pft%column(pi)   = ci
          pft%wtcol(pi)    = wt_nat_pft(gi, m)
          pft%wtlunit(pi)  = pft%wtcol(pi)
          pft%wtgcell(pi)  = pft%wtlunit(pi) * wtlunit2gcell
       end do
    end if

  end subroutine set_landunit_veg_compete
  
  !------------------------------------------------------------------------
  subroutine set_landunit_wet_ice_lake (ltype, gi, li, ci, pi, setdata, glcmask)
    !
    ! !DESCRIPTION: 
    ! Initialize wet_ice_lake landunits that are non-urban (lake, wetland, glacier, glacier_mec)
    !
    ! !USES
    use clmtype
    use clm_varsur, only : wt_lunit, wt_glc_mec, topo_glc_mec
    use clm_varcon, only : istwet, istdlak, istice, istice_mec
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varpar, only : maxpatch_glcmec
    use pftvarcon , only : noveg

    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
    integer , intent(in), optional :: glcmask    ! = 1 where glc requires sfc mass balance
    !
    ! !LOCAL VARIABLES:
    integer  :: m                                ! index
    integer  :: c                                ! column loop index
    integer  :: ctype                            ! column type
    integer  :: ier                              ! error status 
    integer  :: npfts                            ! number of pfts in landunit
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    real(r8) :: wtcol2lunit                      ! col weight in landunit
    !------------------------------------------------------------------------

    ! Set decomposition properties

    if (ltype == istwet) then
       call subgrid_get_gcellinfo(gi, nwetland=npfts)
    else if (ltype == istdlak) then
       call subgrid_get_gcellinfo(gi, nlake=npfts)
    else if (ltype == istice) then 
       call subgrid_get_gcellinfo(gi, nglacier=npfts)
    else if (ltype == istice_mec) then
       call subgrid_get_gcellinfo(gi, nglacier_mec=npfts, glcmask = glcmask)
    else
       write(iulog,*)' set_landunit_wet_ice_lake: ltype of ',ltype,' not valid'
       write(iulog,*)' only istwet, istdlak, istice and istice_mec ltypes are valid'
       call endrun()
    end if

    wtlunit2gcell = wt_lunit(gi, ltype)

    if (npfts > 0) then

       if (npfts /=1 .and. ltype /= istice_mec) then
          write(iulog,*)' set_landunit_wet_ice_lake: compete landunit must'// &
               ' have one pft '
          write(iulog,*)' current value of npfts=',npfts
          call endrun()
       end if

       if (ltype==istice_mec) then   ! multiple columns per landunit

          ! Assume that columns are of type 1 and that each column has its own pft

          ctype = 1
          li = li + 1

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
 
          do m = 1, maxpatch_glcmec

             wtcol2lunit = wt_glc_mec(gi,m)

             if (wtcol2lunit > 0._r8 .or. glcmask == 1) then

                ci = ci + 1
                pi = pi + 1

                col%itype    (ci) = ctype
                col%gridcell (ci) = gi
                col%wtgcell  (ci) = wtcol2lunit * wtlunit2gcell
                col%landunit (ci) = li
                col%wtlunit  (ci) = wtcol2lunit

                ! Set sfc elevation too

                cps%glc_topo(ci) = topo_glc_mec(gi,m)

                ! Set pft properties

                pft%itype    (pi) = noveg
                pft%gridcell (pi) = gi
                pft%wtgcell  (pi) = wtcol2lunit * wtlunit2gcell
                pft%landunit (pi) = li
                pft%wtlunit  (pi) = wtcol2lunit
                pft%column   (pi) = ci
                pft%wtcol    (pi) = 1.0_r8

             endif   ! wtcol2lunit > 0 or glcmask = 1
          enddo      ! loop over columns

       else

          ! Currently assume that each landunit only has only one column 
          ! (of type 1) and that each column has its own pft
       
          wtcol2lunit = 1.0_r8
          ctype = 1

          li = li + 1
          ci = ci + 1
          pi = pi + 1

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
          lun%wtgcell(li)   = wtlunit2gcell

          ! Determine column and properties
          ! For the wet, ice or lake landunits it is assumed that each 
          ! column has its own pft

          col%itype(ci)     = ctype
          col%gridcell (ci) = gi
          col%wtgcell(ci)   = wtcol2lunit * wtlunit2gcell
          col%landunit (ci) = li
          col%wtlunit(ci)   = wtcol2lunit

          ! Set pft properties

          pft%itype(pi)     = noveg
          pft%gridcell (pi) = gi
          pft%wtgcell(pi)   =  wtcol2lunit * wtlunit2gcell
          pft%landunit (pi) = li
          pft%wtlunit(pi)   = wtcol2lunit
          pft%column (pi)   = ci
          pft%wtcol(pi)     = 1.0_r8

       end if   ! ltype = istice_mec
    endif       ! npfts > 0       

  end subroutine set_landunit_wet_ice_lake

  !------------------------------------------------------------------------
  subroutine set_landunit_crop_noncompete (ltype, gi, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize crop landunit without competition
    !
    ! Note about the ltype input argument: This provides the value for this landunit index
    ! (i.e., the crop landunit index). This may differ from the landunit's 'itype' value,
    ! since itype is istsoil if we are running with create_crop_landunit but crop_prog = false.
    !
    ! !USES
    use clmtype
    use clm_varsur, only : wt_lunit, wt_cft
    use clm_varcon, only : istcrop, istsoil
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varctl, only : create_crop_landunit
    use clm_varpar, only : maxpatch_pft, numcft, crop_prog, cft_lb, cft_ub
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: m                                ! index
    integer  :: npfts                            ! number of pfts in landunit
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    !------------------------------------------------------------------------

    ! Set decomposition properties

    call subgrid_get_gcellinfo(gi, ncrop=npfts)
    wtlunit2gcell = wt_lunit(gi, ltype)

    if (npfts > 0) then

       ! Set landunit properties - each column has its own pft
       
       li = li + 1   

       ! Note that we cannot simply use the 'ltype' argument to set itype here,
       ! because ltype will always indicate istcrop
       if ( crop_prog )then
          lun%itype(li) = istcrop
       else
          lun%itype(li) = istsoil
       end if
       lun%ifspecial(li) = .false.
       lun%lakpoi(li)    = .false.
       lun%urbpoi(li)    = .false.
       lun%gridcell (li) = gi
       lun%wtgcell(li) = wtlunit2gcell

       ! Set column and pft properties for this landunit 
       ! (each column has its own pft)

       if (create_crop_landunit) then
          do m = cft_lb, cft_ub
             ci = ci + 1
             pi = pi + 1
             
             col%itype(ci)    = 1
             pft%itype(pi)    = m
             pft%mxy(pi)      = m + 1

             col%gridcell (ci) = gi
             col%landunit (ci) = li
             col%wtlunit(ci) = wt_cft(gi,m)
             col%wtgcell(ci) = col%wtlunit(ci) * wtlunit2gcell

             pft%gridcell (pi) = gi
             pft%landunit (pi) = li
             pft%column (pi) = ci
             pft%wtcol(pi) = 1._r8
             pft%wtlunit(pi) = col%wtlunit(ci)
             pft%wtgcell(pi) = col%wtgcell(ci)
          end do
       end if

    end if
       
  end subroutine set_landunit_crop_noncompete

  !------------------------------------------------------------------------------
  subroutine set_landunit_urban (ltype, gi, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize urban landunits
    !
    ! !USES
    use clm_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, &
                              icol_road_perv, icol_road_imperv, &
                              isturb_tbd, isturb_hd, isturb_md, isturb_MIN
    use clm_varpar   , only : maxpatch_urb
    use clmtype
    use clm_varsur   , only : wt_lunit
    use subgridMod   , only : subgrid_get_gcellinfo
    use UrbanInputMod, only : urbinp
    use decompMod    , only : ldecomp
    use pftvarcon    , only : noveg
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: c             ! column loop index
    integer  :: m             ! index
    integer  :: n             ! urban density type index
    integer  :: ctype         ! column type
    integer  :: npfts         ! number of pfts in landunit
    real(r8) :: wtlunit2gcell ! weight relative to gridcell of landunit
    real(r8) :: wtcol2lunit   ! weight of column with respect to landunit
    real(r8) :: wtlunit_roof  ! weight of roof with respect to landunit
    real(r8) :: wtroad_perv   ! weight of pervious road column with respect to total road
    integer  :: ier           ! error status 
    !------------------------------------------------------------------------

    ! Set decomposition properties, and set variables specific to urban density type

    select case (ltype)
    case (isturb_tbd)
       call subgrid_get_gcellinfo(gi, nurban_tbd=npfts)
    case (isturb_hd)
       call subgrid_get_gcellinfo(gi, nurban_hd=npfts)
    case (isturb_md)
       call subgrid_get_gcellinfo(gi, nurban_md=npfts)
    case default
       write(iulog,*)' set_landunit_urban: unknown ltype: ', ltype
       call endrun()
    end select

    wtlunit2gcell = wt_lunit(gi, ltype)

    n = ltype - isturb_MIN + 1
    wtlunit_roof = urbinp%wtlunit_roof(gi,n)
    wtroad_perv  = urbinp%wtroad_perv(gi,n)

    if (npfts > 0) then

       ! Determine landunit properties - each columns has its own pft

       li = li + 1
       lun%itype    (li) = ltype
       lun%ifspecial(li) = .true.
       lun%lakpoi   (li) = .false.
       lun%urbpoi   (li) = .true.

       lun%gridcell (li) = gi
       lun%wtgcell  (li) = wtlunit2gcell

       ! Loop through columns for this landunit and set the column and pft properties
       ! For the urban landunits it is assumed that each column has its own pft
       
       do m = 1, maxpatch_urb
          
          if (m == 1) then
             ctype = icol_roof
             wtcol2lunit = wtlunit_roof
          else if (m == 2) then
             ctype = icol_sunwall
             wtcol2lunit = (1. - wtlunit_roof)/3
          else if (m == 3) then
             ctype = icol_shadewall
             wtcol2lunit = (1. - wtlunit_roof)/3
          else if (m == 4) then
             ctype = icol_road_imperv
             wtcol2lunit = ((1. - wtlunit_roof)/3) * (1.-wtroad_perv)
          else if (m == 5) then
             ctype = icol_road_perv
             wtcol2lunit = ((1. - wtlunit_roof)/3) * (wtroad_perv)
          end if

          ci = ci + 1
          pi = pi + 1 

          col%itype(ci)    = ctype
          col%gridcell(ci) = gi
          col%wtgcell(ci)  = wtcol2lunit * wtlunit2gcell
          col%landunit(ci) = li
          col%wtlunit(ci)  = wtcol2lunit

          pft%itype(pi)    = noveg
          pft%gridcell(pi) = gi
          pft%wtgcell(pi)  = wtcol2lunit * wtlunit2gcell
          pft%landunit(pi) = li
          pft%wtlunit(pi)  = wtcol2lunit
          pft%column(pi)   = ci
          pft%wtcol(pi)    = 1.0_r8

       end do   ! end of loop through urban columns-pfts
    end if

  end subroutine set_landunit_urban

end module initGridCellsMod
