module initSubgridMod

  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar  , only : lsmlon, lsmlat, numpft, &
	                   maxpatch, maxpatch_pft, maxpatch_cft, &
                           npatch_lake, npatch_glacier, npatch_wet, npatch_crop
  use clm_varctl  , only : allocate_all_vegpfts
  use spmdMod     , only : masterproc
  use nanMod      , only : bigint,nan
  use abortutils  , only : endrun

  implicit none
  private	
  save

  public subgrid_alloc                ! allocates arrays in subgrid_type
  public subgrid_compdown             ! computes "down" data of subgrid_type
  public subgrid_check                ! check subgrid data
  public get_subgrid_size             ! get subgrid global numg,numl,numc,nump
  public get_gcell_info               ! returns g,l,c,p properties
  public set_landunit_veg_compete     ! sets up naturally vegetated landunits
  public set_landunit_wet_ice_lake    ! sets up lake, wetland and glacier lunts
  public set_landunit_crop_noncompete ! sets up crop landunit
  public get_sn_land1d                ! returns g     s->n indices for lunits
  public get_sn_cols1d                ! returns g,l   s->n indices for cols
  public get_sn_pfts1d                ! returns g,l,c s->n indices for pfts

  type subgrid_type
     character(len=32) :: name
     integer , pointer :: g_li(:)
     integer , pointer :: g_lf(:)
     integer , pointer :: g_ln(:)
     integer , pointer :: g_ci(:)
     integer , pointer :: g_cf(:)
     integer , pointer :: g_cn(:)
     integer , pointer :: g_pi(:)
     integer , pointer :: g_pf(:)
     integer , pointer :: g_pn(:)
     integer , pointer :: l_g(:)
     real(r8), pointer :: l_gw(:)
     integer , pointer :: l_ci(:)
     integer , pointer :: l_cf(:)
     integer , pointer :: l_cn(:)
     integer , pointer :: l_pi(:)
     integer , pointer :: l_pf(:)
     integer , pointer :: l_pn(:)
     integer , pointer :: c_g(:)
     real(r8), pointer :: c_gw(:)
     integer , pointer :: c_l(:)
     real(r8), pointer :: c_lw(:)
     integer , pointer :: c_pi(:)
     integer , pointer :: c_pf(:)
     integer , pointer :: c_pn(:)
     integer , pointer :: p_g(:)
     real(r8), pointer :: p_gw(:)
     integer , pointer :: p_l(:)
     real(r8), pointer :: p_lw(:)
     integer , pointer :: p_c(:)
     real(r8), pointer :: p_cw(:)
  end type subgrid_type
  public subgrid_type

  type (subgrid_type),public :: gcellsn
  type (subgrid_type),public :: gcelldc

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgrid_alloc
!
! !INTERFACE:
  subroutine subgrid_alloc(subgrid,ng,nl,nc,np,name)
!
! !DESCRIPTION:
! Allocate subgrid_type
!
! !USES

! !ARGUMENTS
    implicit none
    type(subgrid_type) :: subgrid               ! subgrid to init
    integer , intent(in)  :: ng,nl,nc,np           ! size to init
    character(len=*),intent(in),optional :: name   ! optional name
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
    integer :: ier
    logical :: error
!------------------------------------------------------------------------------

    error = .false.

    if (present(name)) then
       subgrid%name = trim(name)
    else
       subgrid%name = "unnamed"
    endif

    allocate(subgrid%g_li(ng),subgrid%g_lf(ng),subgrid%g_ln(ng),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%g_ci(ng),subgrid%g_cf(ng),subgrid%g_cn(ng),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%g_pi(ng),subgrid%g_pf(ng),subgrid%g_pn(ng),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%l_g (nl),subgrid%l_gw(nl),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%l_ci(nl),subgrid%l_cf(nl),subgrid%l_cn(nl),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%l_pi(nl),subgrid%l_pf(nl),subgrid%l_pn(nl),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%c_g (nc),subgrid%c_gw(nc),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%c_l (nc),subgrid%c_lw(nc),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%c_pi(nc),subgrid%c_pf(nc),subgrid%c_pn(nc),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%p_g (np),subgrid%p_gw(np),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%p_l (np),subgrid%p_lw(np),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%p_c (np),subgrid%p_cw(np),stat=ier)
    if (ier /= 0) error = .true.

    if (error) then
       write(6,*) 'subgrid_alloc ERROR: '
       call endrun()
    endif

    subgrid%g_li(:) = bigint; subgrid%g_lf(:) = bigint
    subgrid%g_ln(:) = bigint
    subgrid%g_ci(:) = bigint; subgrid%g_cf(:) = bigint
    subgrid%g_cn(:) = bigint
    subgrid%g_pi(:) = bigint; subgrid%g_pf(:) = bigint
    subgrid%g_pn(:) = bigint
    subgrid%l_g (:) = bigint; subgrid%l_gw(:) = nan
    subgrid%l_ci(:) = bigint; subgrid%l_cf(:) = bigint
    subgrid%l_cn(:) = bigint
    subgrid%l_pi(:) = bigint; subgrid%l_pf(:) = bigint
    subgrid%l_pn(:) = bigint
    subgrid%c_g (:) = bigint; subgrid%c_gw(:) = nan
    subgrid%c_l (:) = bigint; subgrid%c_lw(:) = nan
    subgrid%c_pi(:) = bigint; subgrid%c_pf(:) = bigint
    subgrid%c_pn(:) = bigint
    subgrid%p_g (:) = bigint; subgrid%p_gw(:) = nan
    subgrid%p_l (:) = bigint; subgrid%p_lw(:) = nan
    subgrid%p_c (:) = bigint; subgrid%p_cw(:) = nan

    end subroutine subgrid_alloc

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgrid_compdown
!
! !INTERFACE:
  subroutine subgrid_compdown(subgrid)
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

! !ARGUMENTS
    implicit none
    type(subgrid_type) :: subgrid
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
    integer :: numg,numl,numc,nump   ! sizes of arrays in subgrid
    integer :: g,l,c,p               ! loop counters
    integer :: curg,curl,curc,curp   ! tracks g,l,c,p indexes in arrays
!------------------------------------------------------------------------------

    call get_subgrid_size(subgrid,numg,numl,numc,nump)

    !--- Set the current c,l,g (curc, curl, curg) to zero for initialization,
    !---   these indices track the current "up" index.
    !--- Loop p through full nump length
    !--- Separately check the p_c, p_l, and p_g indexes for a change in
    !---   the "up" index.
    !--- If there is a change, verify that the current c,l,g is within the 
    !---   valid range, and set c_pi, l_pi, or g_pi to that current c,l,g
    !--- Constantly update the c_pf, l_pf, and g_pf array.  When the
    !---   g, l, c index changes, the *_pf array will be set correctly
    !--- Do the same for numc setting c_li, c_gi, c_lf, c_gf and
    !---   numl setting l_gi, l_gf.

    curc = 0
    curl = 0
    curg = 0
    do p = 1,nump
       if (subgrid%p_c(p) /= curc) then
          curc = subgrid%p_c(p)
          if (curc < 1 .or. curc > numc) then
             write(6,*) 'subgrid_compdown ERROR: p_c ',p,curc,numc
             call endrun()
          endif
          subgrid%c_pi(curc) = p
       endif
       subgrid%c_pf(curc) = p
       subgrid%c_pn(curc) = subgrid%c_pf(curc) - subgrid%c_pi(curc) + 1
       if (subgrid%p_l(p) /= curl) then
          curl = subgrid%p_l(p)
          if (curl < 1 .or. curl > numl) then
             write(6,*) 'subgrid_compdown ERROR: p_l ',p,curl,numl
             call endrun()
          endif
          subgrid%l_pi(curl) = p
       endif
       subgrid%l_pf(curl) = p
       subgrid%l_pn(curl) = subgrid%l_pf(curl) - subgrid%l_pi(curl) + 1
       if (subgrid%p_g(p) /= curg) then
          curg = subgrid%p_g(p)
          if (curg < 1 .or. curg > numg) then
             write(6,*) 'subgrid_compdown ERROR: p_g ',p,curg,numg
             call endrun()
          endif
          subgrid%g_pi(curg) = p
       endif
       subgrid%g_pf(curg) = p
       subgrid%g_pn(curg) = subgrid%g_pf(curg) - subgrid%g_pi(curg) + 1
    enddo

    curg = 0
    curl = 0
    do c = 1,numc
       if (subgrid%c_l(c) /= curl) then
          curl = subgrid%c_l(c)
          if (curl < 1 .or. curl > numl) then
             write(6,*) 'subgrid_compdown ERROR: c_l ',c,curl,numl
             call endrun()
          endif
          subgrid%l_ci(curl) = c
       endif
       subgrid%l_cf(curl) = c
       subgrid%l_cn(curl) = subgrid%l_cf(curl) - subgrid%l_ci(curl) + 1
       if (subgrid%c_g(c) /= curg) then
          curg = subgrid%c_g(c)
          if (curg < 1 .or. curg > numg) then
             write(6,*) 'subgrid_compdown ERROR: c_g ',c,curg,numg
             call endrun()
          endif
          subgrid%g_ci(curg) = c
       endif
       subgrid%g_cf(curg) = c
       subgrid%g_cn(curg) = subgrid%g_cf(curg) - subgrid%g_ci(curg) + 1
    enddo

    curg = 0
    do l = 1,numl
       if (subgrid%l_g(l) /= curg) then
          curg = subgrid%l_g(l)
          if (curg < 1 .or. curg > numg) then
             write(6,*) 'subgrid_compdown ERROR: l_g ',l,curg,numg
             call endrun()
          endif
          subgrid%g_li(curg) = l
       endif
       subgrid%g_lf(curg) = l
       subgrid%g_ln(curg) = subgrid%g_lf(curg) - subgrid%g_li(curg) + 1
    enddo

    end subroutine subgrid_compdown
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgrid_check
!
! !INTERFACE:
  subroutine subgrid_check(subgrid)
!
! !DESCRIPTION:
! Checks and writes out a summary of subgrid data
!
! !USES

! !ARGUMENTS
    implicit none
    type(subgrid_type) :: subgrid
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 2005.11.15  T Craig Creation
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: numg,numl,numc,nump   ! sizes of arrays in subgrid
    integer :: g,l,c,p       ! loop counters
    logical :: error         ! error flag
!------------------------------------------------------------------------------

    if (masterproc) write(6,*) ' '
    if (masterproc) write(6,*) '---subgrid_check: checking ',trim(subgrid%name)
    call get_subgrid_size(subgrid,numg,numl,numc,nump)

    !--- check index ranges ---
    error = .false.
    if (minval(subgrid%g_li) < 1 .or. maxval(subgrid%g_li) > numl) error=.true.
    if (minval(subgrid%g_lf) < 1 .or. maxval(subgrid%g_lf) > numl) error=.true.
    if (minval(subgrid%g_ci) < 1 .or. maxval(subgrid%g_ci) > numc) error=.true.
    if (minval(subgrid%g_cf) < 1 .or. maxval(subgrid%g_cf) > numc) error=.true.
    if (minval(subgrid%g_pi) < 1 .or. maxval(subgrid%g_pi) > nump) error=.true.
    if (minval(subgrid%g_pf) < 1 .or. maxval(subgrid%g_pf) > nump) error=.true.
    if (error) then
       write(6,*) '   subgrid_check: g index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(6,*) '   subgrid_check: g index ranges - OK'

    error = .false.
    if (minval(subgrid%l_g ) < 1 .or. maxval(subgrid%l_g ) > numg) error=.true.
    if (minval(subgrid%l_ci) < 1 .or. maxval(subgrid%l_ci) > numc) error=.true.
    if (minval(subgrid%l_cf) < 1 .or. maxval(subgrid%l_cf) > numc) error=.true.
    if (minval(subgrid%l_pi) < 1 .or. maxval(subgrid%l_pi) > nump) error=.true.
    if (minval(subgrid%l_pf) < 1 .or. maxval(subgrid%l_pf) > nump) error=.true.
    if (error) then
       write(6,*) '   subgrid_check: l index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(6,*) '   subgrid_check: l index ranges - OK'

    error = .false.
    if (minval(subgrid%c_g ) < 1 .or. maxval(subgrid%c_g ) > numg) error=.true.
    if (minval(subgrid%c_l ) < 1 .or. maxval(subgrid%c_l ) > numl) error=.true.
    if (minval(subgrid%c_pi) < 1 .or. maxval(subgrid%c_pi) > nump) error=.true.
    if (minval(subgrid%c_pf) < 1 .or. maxval(subgrid%c_pf) > nump) error=.true.
    if (error) then
       write(6,*) '   subgrid_check: c index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(6,*) '   subgrid_check: c index ranges - OK'

    error = .false.
    if (minval(subgrid%p_g ) < 1 .or. maxval(subgrid%p_g ) > numg) error=.true.
    if (minval(subgrid%p_l ) < 1 .or. maxval(subgrid%p_l ) > numl) error=.true.
    if (minval(subgrid%p_c ) < 1 .or. maxval(subgrid%p_c ) > numc) error=.true.
    if (error) then
       write(6,*) '   subgrid_check: p index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(6,*) '   subgrid_check: p index ranges - OK'

    !--- check that indices in arrays are monotonically increasing ---
    error = .false.
    do g=2,numg
      if (subgrid%g_li(g) < subgrid%g_li(g-1)) error = .true.
      if (subgrid%g_lf(g) < subgrid%g_lf(g-1)) error = .true.
      if (subgrid%g_ci(g) < subgrid%g_ci(g-1)) error = .true.
      if (subgrid%g_cf(g) < subgrid%g_cf(g-1)) error = .true.
      if (subgrid%g_pi(g) < subgrid%g_pi(g-1)) error = .true.
      if (subgrid%g_pf(g) < subgrid%g_pf(g-1)) error = .true.
      if (error) then
         write(6,*) '   subgrid_check: g mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(6,*) '   subgrid_check: g mono increasing - OK'

    error = .false.
    do l=2,numl
      if (subgrid%l_g (l) < subgrid%l_g (l-1)) error = .true.
      if (subgrid%l_ci(l) < subgrid%l_ci(l-1)) error = .true.
      if (subgrid%l_cf(l) < subgrid%l_cf(l-1)) error = .true.
      if (subgrid%l_pi(l) < subgrid%l_pi(l-1)) error = .true.
      if (subgrid%l_pf(l) < subgrid%l_pf(l-1)) error = .true.
      if (error) then
         write(6,*) '   subgrid_check: l mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(6,*) '   subgrid_check: l mono increasing - OK'

    error = .false.
    do c=2,numc
      if (subgrid%c_g (c) < subgrid%c_g (c-1)) error = .true.
      if (subgrid%c_l (c) < subgrid%c_l (c-1)) error = .true.
      if (subgrid%c_pi(c) < subgrid%c_pi(c-1)) error = .true.
      if (subgrid%c_pf(c) < subgrid%c_pf(c-1)) error = .true.
      if (error) then
         write(6,*) '   subgrid_check: c mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(6,*) '   subgrid_check: c mono increasing - OK'

    error = .false.
    do p=2,nump
      if (subgrid%p_g (p) < subgrid%p_g (p-1)) error = .true.
      if (subgrid%p_l (p) < subgrid%p_l (p-1)) error = .true.
      if (subgrid%p_c (p) < subgrid%p_c (p-1)) error = .true.
      if (error) then
         write(6,*) '   subgrid_check: p mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(6,*) '   subgrid_check: p mono increasing - OK'

    !--- check that the tree is internally consistent ---
    error = .false.
    do g = 1, numg
       do l = subgrid%g_li(g),subgrid%g_lf(g)
          if (subgrid%l_g(l) /= g) error = .true.
          do c = subgrid%l_ci(l),subgrid%l_cf(l)
             if (subgrid%c_g(c) /= g) error = .true.
             if (subgrid%c_l(c) /= l) error = .true.
             do p = subgrid%c_pi(c),subgrid%c_pf(c)
                if (subgrid%p_g(p) /= g) error = .true.
                if (subgrid%p_l(p) /= l) error = .true.
                if (subgrid%p_c(p) /= c) error = .true.
                if (error) then
                   write(6,*) '   subgrid_check: tree consistent - ERROR'
                   call endrun()
                endif
             enddo
          enddo
       enddo
    enddo
    if (masterproc) write(6,*) '   subgrid_check: tree consistent - OK'
    if (masterproc) write(6,*) ' '

end subroutine subgrid_check

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_subgrid_size
!
! !INTERFACE:
   subroutine get_subgrid_size(subgrid,numg,numl,numc,nump)
!
! !DESCRIPTION:
!  Returns global size of subgrid g,l,c,p
!
! !USES 
!
! !ARGUMENTS
     implicit none
     type(subgrid_type),intent(in)  :: subgrid
     integer           ,intent(out) :: numg,numl,numc,nump
!
! !REVISION HISTORY:
! 2005.11.15 T Craig Creation
!
!EOP
!
! !LOCAL VARIABLES
!------------------------------------------------------------------------------

   ! extract size directory from subgrid, pick appropriate arrays in subgrid
   numg = size(subgrid%g_pi)
   numl = size(subgrid%l_pi)
   numc = size(subgrid%c_pi)
   nump = size(subgrid%p_g)

   end subroutine get_subgrid_size

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
                           i, j, gi, li, ci, pi, gcell, clm_input)
!
! !DESCRIPTION: 
! Initialize vegetated landunit with competition
!
! !USES
    use clmtype   , only : model_type, gridcell_type, landunit_type, &
                           column_type,pft_type
    use clm_varcon, only : istsoil
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    real(r8), intent(in)    :: wtxy(:,:,:)       ! subgrid patch weights
    integer , intent(in)    :: vegxy(:,:,:)      ! PFT types 
    integer , intent(in)    :: i                 ! 2d longitude index
    integer , intent(in)    :: j                 ! 2d latitude index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    type(subgrid_type)       ,intent(inout) :: gcell
    type(model_type),optional,intent(inout),target :: clm_input  ! clm3
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(i,j,m)
    integer  :: n                                ! loop index
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    integer  :: pitype                           ! pft itype
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    logical  :: setclm                           ! clm present flag
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft

!------------------------------------------------------------------------

    setclm = .false.
    if (present(clm_input)) then
       setclm = .true.
    endif

    ! Set decomposition properties

    call get_gcell_info(i, j, wtxy, nveg=npfts, wtveg=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       if (setclm) then
          lptr => clm_input%g%l
          cptr => clm_input%g%l%c
          pptr => clm_input%g%l%c%p
       endif

       ! Set landunit properties
       
       ncols = 1
       
       li = li + 1
       if (setclm) then
          lptr%ifspecial(li) = .false.
          lptr%lakpoi(li)    = .false.
          lptr%itype(li)     = ltype
       endif
       
       gcell%l_g (li) = gi
       gcell%l_gw(li) = wtlunit2gcell

       ! Set column properties for this landunit (only one column on landunit)
    
       ci = ci + 1
       if (setclm) then
          cptr%itype(ci)    = 1
       endif
       
       gcell%c_g (ci) = gi
       gcell%c_gw(ci) = wtlunit2gcell
       gcell%c_l (ci) = li
       gcell%c_lw(ci) = 1.0_r8

       ! Set pft properties for this landunit

       if (allocate_all_vegpfts) then
          do n = 1,numpft+1
             pi = pi + 1
             pitype = n-1
             if (setclm) then
                pptr%mxy(pi)      = n
                pptr%itype(pi)    = pitype
             endif
             gcell%p_g (pi) = gi
             gcell%p_l (pi) = li
             gcell%p_c (pi) = ci
             gcell%p_gw(pi) = 0.0_r8
             gcell%p_lw(pi) = 0.0_r8
             gcell%p_cw(pi) = 0.0_r8
             do m = 1,maxpatch_pft
                if (vegxy(i,j,m) == pitype .and. wtxy(i,j,m) > 0._r8) then
                   gcell%p_gw(pi)  = gcell%p_gw(pi) + wtxy(i,j,m)
                   gcell%p_lw(pi)  = gcell%p_lw(pi) + wtxy(i,j,m) / wtlunit2gcell
                   gcell%p_cw(pi)  = gcell%p_cw(pi) + wtxy(i,j,m) / wtlunit2gcell
                end if
             end do

          end do
       else
          do m = 1,maxpatch_pft
             if (wtxy(i,j,m) > 0._r8) then
                pi = pi + 1
                if (setclm) then
                   pptr%mxy(pi)      = m
                   pptr%itype(pi)    = vegxy(i,j,m)
                endif
                gcell%p_g (pi) = gi
                gcell%p_gw(pi) = wtxy(i,j,m)
                gcell%p_l (pi) = li
                gcell%p_lw(pi) = wtxy(i,j,m) / wtlunit2gcell
                gcell%p_c (pi) = ci
                gcell%p_cw(pi) = wtxy(i,j,m) / wtlunit2gcell
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
                           i, j, gi, li, ci, pi, gcell, clm_input)
!
! !DESCRIPTION: 
! Initialize wet_ice_lake landunits that are non-urban (lake, wetland, glacier)
!
! !USES
    use clmtype   , only : model_type, gridcell_type, landunit_type, &
                           column_type,pft_type
    use clm_varcon, only : istice, istwet, istdlak
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    real(r8), intent(in)    :: wtxy(:,:,:)       ! subgrid patch weights
    integer , intent(in)    :: vegxy(:,:,:)      ! PFT types 
    integer , intent(in)    :: i                 ! 2d longitude index
    integer , intent(in)    :: j                 ! 2d latitude index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    type(subgrid_type)       ,intent(inout) :: gcell
    type(model_type),optional,intent(inout),target :: clm_input  ! clm3
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(i,j,m)
    integer  :: c                                ! column loop index
    integer  :: ctype                            ! column type
    integer  :: ier                              ! error status 
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    real(r8) :: wtcol2lunit                      ! col weight in landunit
    logical  :: setclm                           ! clm present flag
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft
!------------------------------------------------------------------------

    setclm = .false.
    if (present(clm_input)) then
       setclm = .true.
    endif

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

       if (setclm) then
          lptr => clm_input%g%l
          cptr => clm_input%g%l%c
          pptr => clm_input%g%l%c%p
       endif
       
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

       ! Determine landunit properties 

       li = li + 1
       if (setclm) then
          lptr%itype(li)     = ltype
          lptr%ifspecial(li) = .true.
          if (ltype == istdlak) then
             lptr%lakpoi(li) = .true.
          else
             lptr%lakpoi(li) = .false.
          end if
       endif
       
       gcell%l_g (li) = gi
       gcell%l_gw(li) = wtlunit2gcell

       ! Determine column and properties
       ! For the wet, ice or lake landunits it is assumed that each 
       ! column has its own pft
       
       ci = ci + 1
       pi = pi + 1 
       
       if (setclm) then
          cptr%itype(ci)    = ctype
       endif
       
       gcell%c_g (ci) = gi
       gcell%c_gw(ci) = wtcol2lunit * wtlunit2gcell
       gcell%c_l (ci) = li
       gcell%c_lw(ci) = wtcol2lunit

       if (setclm) then
          pptr%mxy(pi)      = m
          pptr%itype(pi)    = vegxy(i,j,m)
       endif
       
       gcell%p_g (pi) = gi
       gcell%p_gw(pi) = wtcol2lunit * wtlunit2gcell
       gcell%p_l (pi) = li
       gcell%p_lw(pi) = wtcol2lunit
       gcell%p_c (pi) = ci
       gcell%p_cw(pi) = 1.0_r8
    end if
       
  end subroutine set_landunit_wet_ice_lake

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_crop_noncompete
!
! !INTERFACE:
  subroutine set_landunit_crop_noncompete (ltype, wtxy, vegxy, &
                           i, j, gi, li, ci, pi, gcell, clm_input)
!
! !DESCRIPTION: 
! Initialize crop landunit without competition
!
! !USES
    use clmtype   , only : model_type, gridcell_type, landunit_type, &
                           column_type,pft_type
    use clm_varpar, only : npatch_crop, npatch_glacier
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    real(r8), intent(in)    :: wtxy(:,:,:)       ! subgrid patch weights
    integer , intent(in)    :: vegxy(:,:,:)      ! PFT types 
    integer , intent(in)    :: i                 ! 2d longitude index
    integer , intent(in)    :: j                 ! 2d latitude index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    type(subgrid_type)       ,intent(inout) :: gcell
    type(model_type),optional,intent(inout),target :: clm_input  ! clm3
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                                ! m index in wtxy(i,j,m)
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    logical  :: setclm                           ! clm present flag
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft
!------------------------------------------------------------------------

    setclm = .false.
    if (present(clm_input)) then
       setclm = .true.
    endif

    ! Set decomposition properties

    call get_gcell_info(i, j, wtxy, ncrop=npfts, wtcrop=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       if (setclm) then
          lptr => clm_input%g%l
          cptr => clm_input%g%l%c
          pptr => clm_input%g%l%c%p
       endif
       
       ! Set landunit properties - each column has its own pft
       
       ncols = npfts
       
       li = li + 1   
       if (setclm) then
          lptr%itype(li)     = ltype
          lptr%ifspecial(li) = .false.
          lptr%lakpoi(li)    = .false.
       endif
       
       gcell%l_g (li) = gi
       gcell%l_gw(li) = wtlunit2gcell

       ! Set column and pft properties for this landunit 
       ! (each column has its own pft)

       do m = npatch_glacier+1, npatch_crop
          if (wtxy(i,j,m) > 0._r8) then
             ci = ci + 1
             pi = pi + 1
             
             if (setclm) then
                cptr%itype(ci)    = 1
                pptr%itype(pi)    = vegxy(i,j,m)
                pptr%mxy(pi)      = m
             endif
             
             gcell%c_g (ci) = gi
             gcell%c_gw(ci) = wtxy(i,j,m)
             gcell%c_l (ci) = li
             gcell%c_lw(ci) = wtxy(i,j,m) / wtlunit2gcell

             gcell%p_g (pi) = gi
             gcell%p_gw(pi) = wtxy(i,j,m)
             gcell%p_l (pi) = li
             gcell%p_lw(pi) = wtxy(i,j,m) / wtlunit2gcell
             gcell%p_c (pi) = ci
             gcell%p_cw(pi) = 1.0_r8
          end if
       end do

    end if
       
  end subroutine set_landunit_crop_noncompete

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
     use clmtype, only : nameg,namel,namec,namep        
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
           snindex(l) = gcellsn%l_g(l)
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
     use clmtype, only : nameg,namel,namec,namep        
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
           snindex(c) = gcellsn%c_g(c)
        else if (type1d == namel) then
           snindex(c) = gcellsn%c_l(c)
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
     use clmtype, only : nameg,namel,namec,namep        
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
           snindex(p) = gcellsn%p_g(p)
        else if (type1d == namel) then
           snindex(p) = gcellsn%p_l(p)
        else if (type1d == namec) then
           snindex(p) = gcellsn%p_c(p)
        end if
     end do

   end subroutine get_sn_pfts1d

!------------------------------------------------------------------------------

end module initSubgridMod
