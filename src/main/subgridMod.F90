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
  public subgrid_alloc                ! allocates arrays in subgrid_type
  public subgrid_get_indexes          ! returns subgrid indices for gridcells
  public subgrid_get_gcellinfo        ! Returns g,l,c,p properties from wtxy

  type subgrid_type
     character(len=32) :: name
     integer , pointer :: g_li(:)
     integer , pointer :: g_ln(:)
     integer , pointer :: g_ci(:)
     integer , pointer :: g_cn(:)
     integer , pointer :: g_pi(:)
     integer , pointer :: g_pn(:)
!     integer , pointer :: g_lf(:)
!     integer , pointer :: g_cf(:)
!     integer , pointer :: g_pf(:)
!     integer , pointer :: l_g(:)
!     real(r8), pointer :: l_gw(:)
!     integer , pointer :: l_ci(:)
!     integer , pointer :: l_cf(:)
!     integer , pointer :: l_cn(:)
!     integer , pointer :: l_pi(:)
!     integer , pointer :: l_pf(:)
!     integer , pointer :: l_pn(:)
!     integer , pointer :: c_g(:)
!     real(r8), pointer :: c_gw(:)
!     integer , pointer :: c_l(:)
!     real(r8), pointer :: c_lw(:)
!     integer , pointer :: c_pi(:)
!     integer , pointer :: c_pf(:)
!     integer , pointer :: c_pn(:)
!     integer , pointer :: p_g(:)
!     real(r8), pointer :: p_gw(:)
!     integer , pointer :: p_l(:)
!     real(r8), pointer :: p_lw(:)
!     integer , pointer :: p_c(:)
!     real(r8), pointer :: p_cw(:)
  end type subgrid_type
  public subgrid_type

  type (subgrid_type),public :: gcelldc
  type (subgrid_type),public :: gcellsn

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

    allocate(subgrid%g_li(ng),subgrid%g_ln(ng),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%g_ci(ng),subgrid%g_cn(ng),stat=ier)
    if (ier /= 0) error = .true.
    allocate(subgrid%g_pi(ng),subgrid%g_pn(ng),stat=ier)
    if (ier /= 0) error = .true.

    if (error) then
       write(6,*) 'subgrid_alloc ERROR: '
       call endrun()
    endif

    subgrid%g_li(:) = bigint
    subgrid%g_ln(:) = bigint
    subgrid%g_ci(:) = bigint
    subgrid%g_cn(:) = bigint
    subgrid%g_pi(:) = bigint
    subgrid%g_pn(:) = bigint

    end subroutine subgrid_alloc

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgrid_get_indexes
!
! !INTERFACE:
   subroutine subgrid_get_indexes(gdc,gsn,type1d,dci,dcf,sni,snf)
!
! !DESCRIPTION:
!  Gets indices for dc2sn mapping routines
!
! !USES:
   use clmtype   , only : nameg, namel, namec, namep
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: gdc,gsn
    character(len=*), intent(in) :: type1d
    integer, intent(out) :: dci,dcf,sni,snf
!
! !REVISION HISTORY:
! 2005.11.15  T Craig Extracted from map_*_* subroutines
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: error
!------------------------------------------------------------------------------

     error = .false.
     if (type1d == nameg) then
        dci = gdc
        dcf = gdc
        sni = gsn
        snf = gsn
     else if (type1d == namel) then
        dci = gcelldc%g_li(gdc)
        dcf = dci + gcelldc%g_ln(gdc) - 1
        sni = gcellsn%g_li(gsn)
     else if (type1d == namec) then
        dci = gcelldc%g_ci(gdc)
        dcf = dci + gcelldc%g_cn(gdc) - 1
        sni = gcellsn%g_ci(gsn)
     else if (type1d == namep) then
        dci = gcelldc%g_pi(gdc)
        dcf = dci + gcelldc%g_pn(gdc) - 1
        sni = gcellsn%g_pi(gsn)
     end if
     snf = sni + dcf - dci

end subroutine subgrid_get_indexes

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
