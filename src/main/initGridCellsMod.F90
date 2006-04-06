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
!
! !PIVATE MEMBER FUNCTIONS:
  private set_clm_gptrs
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
  subroutine initGridcells (vegxy, wtxy) 
!
! !DESCRIPTION: 
! Initialize sub-grid mapping and allocates space for derived type hierarchy.
! For each land gridcell determine landunit, column and pft properties.
!
! !USES
    use clmtype         , only : clm3, gridcell_type, landunit_type, &
                                 column_type, pft_type
    use domainMod       , only : ldomain
    use decompMod       , only : ldecomp
    use clm_varpar      , only : lsmlon, lsmlat, maxpatch 
    use clm_varcon      , only : istsoil, istice, istwet, istdlak, isturb
    use initSubgridMod  , only : get_gcell_info, &
                                 set_landunit_veg_compete, &
                                 set_landunit_wet_ice_lake, &
    	                         set_landunit_crop_noncompete, &
                                 gcelldc, gcellsn, &
                                 subgrid_alloc,subgrid_compdown,subgrid_check
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
       do i = 1, lsmlon
!tcx fix this so numg,numl,numc,nump global are stored somewhere, not recomputed
          if (ldecomp%ij2gdc(i,j) > 0) then
             call get_gcell_info (i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)
             numg = numg + 1
             numl = numl + ilunits
             numc = numc + icols
             nump = nump + ipfts
          end if
       end do
    end do

    ! Dynamic memory allocation

    call subgrid_alloc(gcellsn,numg,numl,numc,nump,'gcellsn')
    call subgrid_alloc(gcelldc,numg,numl,numc,nump,'gcelldc')

    ! For each land gridcell on global grid determine landunit, column and pft properties

    li    = 0
    ci    = 0
    pi    = 0
    li_sn = 0
    ci_sn = 0
    pi_sn = 0

    !----- Set gcelldc and clm3 variables -----
    do gi = 1,numg

       i     = ldecomp%gdc2i(gi)
       j     = ldecomp%gdc2j(gi)

       ! derived values
       gptr%itype(gi) = 1

       ! Determine naturally vegetated landunit

       call set_landunit_veg_compete(               &
            ltype=istsoil, wtxy=wtxy, vegxy=vegxy,  &
            i=i, j=j, gi=gi, li=li, ci=ci, pi=pi,   &
            gcell=gcelldc, clm_input=clm3)

       ! Determine crop landunit

       call set_landunit_crop_noncompete(           &
            ltype=istsoil, wtxy=wtxy, vegxy=vegxy,  &
            i=i, j=j, gi=gi, li=li, ci=ci, pi=pi,   &
            gcell=gcelldc, clm_input=clm3)

       ! Determine lake, wetland and glacier landunits 

       call set_landunit_wet_ice_lake(              &
            ltype=istdlak, wtxy=wtxy, vegxy=vegxy,  &
            i=i, j=j, gi=gi, li=li, ci=ci, pi=pi,   &
            gcell=gcelldc, clm_input=clm3)

       call set_landunit_wet_ice_lake(              &
            ltype=istwet, wtxy=wtxy, vegxy=vegxy,   &
            i=i, j=j, gi=gi, li=li, ci=ci, pi=pi,   &
            gcell=gcelldc, clm_input=clm3)

       call set_landunit_wet_ice_lake(              &
            ltype=istice, wtxy=wtxy, vegxy=vegxy,   &
            i=i, j=j, gi=gi, li=li, ci=ci, pi=pi,   &
            gcell=gcelldc, clm_input=clm3)

    enddo

    !----- Set gcellsn -----
    do gi = 1,numg
       gi_sn = gi
       i_sn  = ldecomp%gsn2i(gi_sn)
       j_sn  = ldecomp%gsn2j(gi_sn)

       call set_landunit_veg_compete(                               &
            ltype=istsoil, wtxy=wtxy, vegxy=vegxy,                  &
            i=i_sn, j=j_sn, gi=gi_sn, li=li_sn, ci=ci_sn, pi=pi_sn, &
            gcell=gcellsn)

       call set_landunit_crop_noncompete(                           &
            ltype=istsoil, wtxy=wtxy, vegxy=vegxy,                  &
            i=i_sn, j=j_sn, gi=gi_sn, li=li_sn, ci=ci_sn, pi=pi_sn, &
            gcell=gcellsn)

       call set_landunit_wet_ice_lake(                              &
            ltype=istdlak, wtxy=wtxy, vegxy=vegxy,                  &
            i=i_sn, j=j_sn, gi=gi_sn, li=li_sn, ci=ci_sn, pi=pi_sn, &
            gcell=gcellsn)

       call set_landunit_wet_ice_lake(                              &
            ltype=istwet, wtxy=wtxy, vegxy=vegxy,                   &
            i=i_sn, j=j_sn, gi=gi_sn, li=li_sn, ci=ci_sn, pi=pi_sn, &
            gcell=gcellsn)

       call set_landunit_wet_ice_lake(                              &
            ltype=istice, wtxy=wtxy, vegxy=vegxy,                   &
            i=i_sn, j=j_sn, gi=gi_sn, li=li_sn, ci=ci_sn, pi=pi_sn, &
            gcell=gcellsn)

    end do  

    ! Fill in subgrid datatypes

    call subgrid_compdown(gcelldc)
    call subgrid_check(gcelldc)
    call subgrid_compdown(gcellsn)
    call subgrid_check(gcellsn)

    ! Set clm3 pointers for g,l,c,p indexes

    call set_clm_gptrs()

  end subroutine initGridcells

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_clm_gptrs
!
! !INTERFACE:
  subroutine set_clm_gptrs()
!
! !DESCRIPTION: 
! Initialize clmtype gptr properties
!
! !USES
  use clmtype       , only : clm3, gridcell_type, landunit_type, &
                             column_type, pft_type, model_type
  use domainMod     , only : ldomain,adomain
  use decompMod     , only : ldecomp
  use areaMod       , only : gridmap_setptrs
  use clm_atmlnd    , only : gridmap_a2l
  use initSubgridMod, only : gcelldc
  use decompMod     , only : get_proc_global
  use shr_const_mod , only : SHR_CONST_PI
  implicit none
!
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Created by T Craig 2005.11.15
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: i,j,g,l,c,p,ia,ja,n
  integer :: numg     ! global number of gridcells
  integer :: numl     ! global number of landunits
  integer :: numc     ! global number of columns
  integer :: nump     ! global number of pfts
  type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
  type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
  type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
  type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
  integer            ,pointer   :: n_a2l(:,:)   ! number of overlapping cells
  integer            ,pointer   :: i_a2l(:,:,:) ! i index of overlap input cell
  integer            ,pointer   :: j_a2l(:,:,:) ! j index of overlap input cell

!------------------------------------------------------------------------

    ! Set pointers into derived types for this module

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p
    
    call get_proc_global(numg,numl,numc,nump)

    ! pointers

    gptr%luni     => gcelldc%g_li
    gptr%lunf     => gcelldc%g_lf
    gptr%coli     => gcelldc%g_ci
    gptr%colf     => gcelldc%g_cf
    gptr%pfti     => gcelldc%g_pi
    gptr%pftf     => gcelldc%g_pf

    lptr%gridcell => gcelldc%l_g
    lptr%wtgcell  => gcelldc%l_gw
    lptr%coli     => gcelldc%l_ci
    lptr%colf     => gcelldc%l_cf
    lptr%pfti     => gcelldc%l_pi
    lptr%pftf     => gcelldc%l_pf

    cptr%gridcell => gcelldc%c_g
    cptr%wtgcell  => gcelldc%c_gw
    cptr%landunit => gcelldc%c_l
    cptr%wtlunit  => gcelldc%c_lw
    cptr%pfti     => gcelldc%c_pi
    cptr%pftf     => gcelldc%c_pf

    pptr%gridcell => gcelldc%p_g
    pptr%wtgcell  => gcelldc%p_gw
    pptr%landunit => gcelldc%p_l
    pptr%wtlunit  => gcelldc%p_lw
    pptr%column   => gcelldc%p_c
    pptr%wtcol    => gcelldc%p_cw

    gptr%nlandunits => gcelldc%g_ln
    gptr%ncolumns   => gcelldc%g_cn
    gptr%npfts      => gcelldc%g_pn
    lptr%ncolumns   => gcelldc%l_cn
    lptr%npfts      => gcelldc%l_pn
    cptr%npfts      => gcelldc%c_pn

    ! Set lats/lons

    do g = 1,numg

       i = ldecomp%gdc2i(g)
       j = ldecomp%gdc2j(g)

       gptr%ixy(g) = i 
       gptr%jxy(g) = j 
       gptr%latdeg(g) = ldomain%latc(i,j) 
       gptr%londeg(g) = ldomain%lonc(i,j) 
       gptr%lat(g)    = ldomain%latc(i,j) * SHR_CONST_PI/180._r8  
       gptr%lon(g)    = ldomain%lonc(i,j) * SHR_CONST_PI/180._r8
       gptr%area(g)   = ldomain%area(i,j)

    end do

    ! Set "atm" lats/lons in gridcell

    call gridmap_setptrs(gridmap_a2l,n_ovr=n_a2l,i_ovr=i_a2l,j_ovr=j_a2l)

    do g = 1,numg
       i = ldecomp%gdc2i(g)
       j = ldecomp%gdc2j(g)
       if (n_a2l(i,j) /= 1) then
          write(6,*) 'set_clm_gptrs ERROR: n_a2l must be one, ',n_a2l(i,j)
          call endrun()
       endif
       do n = 1,n_a2l(i,j)
          ia = i_a2l(i,j,n)
          ja = j_a2l(i,j,n)
          gptr%londeg_a(g) = adomain%lonc(ia,ja)
          gptr%latdeg_a(g) = adomain%latc(ia,ja)
       enddo
    enddo

    gptr%lon_a(:) = gptr%londeg_a(:) * SHR_CONST_PI/180._r8  
    gptr%lat_a(:) = gptr%latdeg_a(:) * SHR_CONST_PI/180._r8  

  end subroutine set_clm_gptrs

!------------------------------------------------------------------------

end module initGridcellsMod
