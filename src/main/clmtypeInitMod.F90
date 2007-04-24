#include <misc.h>
#include <preproc.h>

module clmtypeInitMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clmtypeInitMod
!
! !DESCRIPTION:
! Allocate clmtype components and initialize them to signaling NaN.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use nanMod
  use clmtype
  use clm_varpar, only: maxpatch_pft, nlevsno, nlevsoi, numrad, nlevlak, numpft, ndst, nvoc
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initClmtype
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type
  private :: init_energy_balance_type
  private :: init_water_balance_type
  private :: init_pft_ecophys_constants
  private :: init_pft_DGVMecophys_constants
  private :: init_pft_pstate_type
  private :: init_pft_epv_type
  private :: init_pft_pdgvstate_type
  private :: init_pft_estate_type
  private :: init_pft_wstate_type
  private :: init_pft_cstate_type
  private :: init_pft_nstate_type
  private :: init_pft_eflux_type
  private :: init_pft_mflux_type
  private :: init_pft_wflux_type
  private :: init_pft_cflux_type
  private :: init_pft_nflux_type
  private :: init_pft_vflux_type
  private :: init_pft_dflux_type
  private :: init_column_pstate_type
  private :: init_column_estate_type
  private :: init_column_wstate_type
  private :: init_column_cstate_type
  private :: init_column_nstate_type
  private :: init_column_eflux_type
  private :: init_column_wflux_type
  private :: init_column_cflux_type
  private :: init_column_nflux_type
  private :: init_landunit_pstate_type
  private :: init_gridcell_pstate_type
  private :: init_gridcell_wflux_type
!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initClmtype
!
! !INTERFACE:
  subroutine initClmtype()
!
! !DESCRIPTION:
! Initialize clmtype components to signaling nan
! The following clmtype components should NOT be initialized here
! since they are set in routine clm_map which is called before this
! routine is invoked
!    *%area, *%wt, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
!    *%ifspecial, *%ityplun, *%itype
!    *%pfti, *%pftf, *%pftn
!    *%coli, *%colf, *%coln
!    *%luni, *%lunf, *%lunn
!
! !USES:
    use decompMod , only : get_proc_bounds, get_proc_global
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    call init_pft_type     (begp, endp, clm3%g%l%c%p)
    call init_column_type  (begc, endc, clm3%g%l%c)
    call init_landunit_type(begl, endl, clm3%g%l)
    call init_gridcell_type(begg, endg, clm3%g)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    ! pft DGVM-specific ecophysiological constants

    call init_pft_DGVMecophys_constants()

    ! energy balance structures (all levels)

    call init_energy_balance_type(begp, endp, clm3%g%l%c%p%pebal)
    call init_energy_balance_type(begc, endc, clm3%g%l%c%cebal)
    call init_energy_balance_type(begl, endl, clm3%g%l%lebal)
    call init_energy_balance_type(begg, endg, clm3%g%gebal)
    call init_energy_balance_type(1,       1, clm3%mebal)

    ! water balance structures (all levels)

    call init_water_balance_type(begp, endp, clm3%g%l%c%p%pwbal)
    call init_water_balance_type(begc, endc, clm3%g%l%c%cwbal)
    call init_water_balance_type(begl, endl, clm3%g%l%lwbal)
    call init_water_balance_type(begg, endg, clm3%g%gwbal)
    call init_water_balance_type(1,       1, clm3%mwbal)

    ! carbon balance structures (pft and column levels)

    call init_carbon_balance_type(begp, endp, clm3%g%l%c%p%pcbal)
    call init_carbon_balance_type(begc, endc, clm3%g%l%c%ccbal)

    ! nitrogen balance structures (pft and column levels)

    call init_nitrogen_balance_type(begp, endp, clm3%g%l%c%p%pnbal)
    call init_nitrogen_balance_type(begc, endc, clm3%g%l%c%cnbal)

    ! pft physical state variables at pft level and averaged to the column

    call init_pft_pstate_type(begp, endp, clm3%g%l%c%p%pps)
    call init_pft_pstate_type(begc, endc, clm3%g%l%c%cps%pps_a)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(begp, endp, clm3%g%l%c%p%pepv)

    ! pft DGVM state variables at pft level and averaged to column

    call init_pft_pdgvstate_type(begp, endp, clm3%g%l%c%p%pdgvs)
    call init_pft_pdgvstate_type(begc, endc, clm3%g%l%c%cdgvs%pdgvs_a)

    ! pft energy state variables at the pft level and averaged to the column

    call init_pft_estate_type(begp, endp, clm3%g%l%c%p%pes)
    call init_pft_estate_type(begc, endc, clm3%g%l%c%ces%pes_a)

    ! pft water state variables at the pft level and averaged to the column

    call init_pft_wstate_type(begp, endp, clm3%g%l%c%p%pws)
    call init_pft_wstate_type(begc, endc, clm3%g%l%c%cws%pws_a)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pcs)
    call init_pft_cstate_type(begc, endc, clm3%g%l%c%ccs%pcs_a)
    ! 4/14/05: PET
    ! Adding isotope code
    call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pc13s)
    call init_pft_cstate_type(begc, endc, clm3%g%l%c%cc13s%pcs_a)

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(begp, endp, clm3%g%l%c%p%pns)
    call init_pft_nstate_type(begc, endc, clm3%g%l%c%cns%pns_a)

    ! pft energy flux variables at pft level and averaged to column

    call init_pft_eflux_type(begp, endp, clm3%g%l%c%p%pef)
    call init_pft_eflux_type(begc, endc, clm3%g%l%c%cef%pef_a)

    ! pft momentum flux variables at pft level and averaged to the column

    call init_pft_mflux_type(begp, endp, clm3%g%l%c%p%pmf)
    call init_pft_mflux_type(begc, endc, clm3%g%l%c%cmf%pmf_a)

    ! pft water flux variables

    call init_pft_wflux_type(begp, endp, clm3%g%l%c%p%pwf)
    call init_pft_wflux_type(begc, endc, clm3%g%l%c%cwf%pwf_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pcf)
    call init_pft_cflux_type(begc, endc, clm3%g%l%c%ccf%pcf_a)
    ! 4/14/05: PET
    ! Adding isotope code
    call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pc13f)
    call init_pft_cflux_type(begc, endc, clm3%g%l%c%cc13f%pcf_a)
    

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(begp, endp, clm3%g%l%c%p%pnf)
    call init_pft_nflux_type(begc, endc, clm3%g%l%c%cnf%pnf_a)

    ! pft VOC flux variables at pft level and averaged to column

    call init_pft_vflux_type(begp, endp, clm3%g%l%c%p%pvf)
    call init_pft_vflux_type(begc, endc, clm3%g%l%c%cvf%pvf_a)

    ! pft dust flux variables at pft level and averaged to column

    call init_pft_dflux_type(begp, endp, clm3%g%l%c%p%pdf)
    call init_pft_dflux_type(begc, endc, clm3%g%l%c%cdf%pdf_a)

    ! column physical state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_pstate_type(begc, endc, clm3%g%l%c%cps)
    call init_column_pstate_type(begl, endl, clm3%g%l%lps%cps_a)
    call init_column_pstate_type(begg, endg, clm3%g%gps%cps_a)
    call init_column_pstate_type(1,       1, clm3%mps%cps_a)

    ! column energy state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_estate_type(begc, endc, clm3%g%l%c%ces)
    call init_column_estate_type(begl, endl, clm3%g%l%les%ces_a)
    call init_column_estate_type(begg, endg, clm3%g%ges%ces_a)
    call init_column_estate_type(1,       1, clm3%mes%ces_a)

    ! column water state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_wstate_type(begc, endc, clm3%g%l%c%cws)
    call init_column_wstate_type(begl, endl, clm3%g%l%lws%cws_a)
    call init_column_wstate_type(begg, endg, clm3%g%gws%cws_a)
    call init_column_wstate_type(1,       1, clm3%mws%cws_a)

    ! column carbon state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_cstate_type(begc, endc, clm3%g%l%c%ccs)
    call init_column_cstate_type(begl, endl, clm3%g%l%lcs%ccs_a)
    call init_column_cstate_type(begg, endg, clm3%g%gcs%ccs_a)
    call init_column_cstate_type(1,       1, clm3%mcs%ccs_a)
    ! 4/14/05: PET
    ! Adding isotope code
    call init_column_cstate_type(begc, endc, clm3%g%l%c%cc13s)

    ! column nitrogen state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_nstate_type(begc, endc, clm3%g%l%c%cns)
    call init_column_nstate_type(begl, endl, clm3%g%l%lns%cns_a)
    call init_column_nstate_type(begg, endg, clm3%g%gns%cns_a)
    call init_column_nstate_type(1,       1, clm3%mns%cns_a)

    ! column energy flux variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_eflux_type(begc, endc, clm3%g%l%c%cef)
    call init_column_eflux_type(begl, endl, clm3%g%l%lef%cef_a)
    call init_column_eflux_type(begg, endg, clm3%g%gef%cef_a)
    call init_column_eflux_type(1,       1, clm3%mef%cef_a)

    ! column water flux variables at column level and averaged to
    ! landunit, gridcell and model level

    call init_column_wflux_type(begc, endc, clm3%g%l%c%cwf)
    call init_column_wflux_type(begl, endl, clm3%g%l%lwf%cwf_a)
    call init_column_wflux_type(begg, endg, clm3%g%gwf%cwf_a)
    call init_column_wflux_type(1,       1, clm3%mwf%cwf_a)

    ! column carbon flux variables at column level

    call init_column_cflux_type(begc, endc, clm3%g%l%c%ccf)
    ! 4/14/05: PET
    ! Adding isotope code
    call init_column_cflux_type(begc, endc, clm3%g%l%c%cc13f)

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(begc, endc, clm3%g%l%c%cnf)

    ! land unit physical state variables

    call init_landunit_pstate_type(begl, endl, clm3%g%l%lps)

    ! gridcell DGVM variables
    call init_gridcell_dgvstate_type(begg, endg, clm3%g%gdgvs)

    ! gridcell physical state variables

    call init_gridcell_pstate_type(begg, endg, clm3%g%gps)

    ! gridcell: water flux variables

    call init_gridcell_wflux_type(begg, endg, clm3%g%gwf)

  end subroutine initClmtype

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_type
!
! !INTERFACE:
  subroutine init_pft_type (beg, end, p)
!
! !DESCRIPTION:
! Initialize components of pft_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(pft_type), intent(inout):: p
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(p%gridcell(beg:end),p%wtgcell(beg:end))
    allocate(p%landunit(beg:end),p%wtlunit(beg:end))
    allocate(p%column  (beg:end),p%wtcol  (beg:end))

    allocate(p%itype(beg:end))
    allocate(p%mxy(beg:end))

  end subroutine init_pft_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_type
!
! !INTERFACE:
  subroutine init_column_type (beg, end, c)
!
! !DESCRIPTION:
! Initialize components of column_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(column_type), intent(inout):: c
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(c%gridcell(beg:end),c%wtgcell(beg:end))
   allocate(c%landunit(beg:end),c%wtlunit(beg:end))

   allocate(c%pfti(beg:end),c%pftf(beg:end),c%npfts(beg:end))

   allocate(c%itype(beg:end))

  end subroutine init_column_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_type
!
! !INTERFACE:
  subroutine init_landunit_type (beg, end,l)
!
! !DESCRIPTION:
! Initialize components of landunit_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(landunit_type), intent(inout):: l
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(l%gridcell(beg:end),l%wtgcell(beg:end))

   allocate(l%coli(beg:end),l%colf(beg:end),l%ncolumns(beg:end))
   allocate(l%pfti(beg:end),l%pftf(beg:end),l%npfts   (beg:end))

   allocate(l%itype(beg:end))
   allocate(l%ifspecial(beg:end))
   allocate(l%lakpoi(beg:end))

  end subroutine init_landunit_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_type
!
! !INTERFACE:
  subroutine init_gridcell_type (beg, end,g)
!
! !DESCRIPTION:
! Initialize components of gridcell_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(gridcell_type), intent(inout):: g
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(g%luni(beg:end),g%lunf(beg:end),g%nlandunits(beg:end))
   allocate(g%coli(beg:end),g%colf(beg:end),g%ncolumns  (beg:end))
   allocate(g%pfti(beg:end),g%pftf(beg:end),g%npfts     (beg:end))

   allocate(g%area(beg:end))
   allocate(g%lat(beg:end))
   allocate(g%lon(beg:end))
   allocate(g%latdeg(beg:end))
   allocate(g%londeg(beg:end))
   allocate(g%lat_a(beg:end))
   allocate(g%lon_a(beg:end))
   allocate(g%latdeg_a(beg:end))
   allocate(g%londeg_a(beg:end))

  end subroutine init_gridcell_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_energy_balance_type
!
! !INTERFACE:
  subroutine init_energy_balance_type(beg, end, ebal)
!
! !DESCRIPTION:
! Initialize energy balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(energy_balance_type), intent(inout):: ebal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ebal%errsoi(beg:end))
    allocate(ebal%errseb(beg:end))
    allocate(ebal%errsol(beg:end))
    allocate(ebal%errlon(beg:end))

    ebal%errsoi(beg:end) = nan
    ebal%errseb(beg:end) = nan
    ebal%errsol(beg:end) = nan
    ebal%errlon(beg:end) = nan

  end subroutine init_energy_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_water_balance_type
!
! !INTERFACE:
  subroutine init_water_balance_type(beg, end, wbal)
!
! !DESCRIPTION:
! Initialize water balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(water_balance_type), intent(inout):: wbal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(wbal%begwb(beg:end))
    allocate(wbal%endwb(beg:end))
    allocate(wbal%errh2o(beg:end))

    wbal%begwb(beg:end) = nan
    wbal%endwb(beg:end) = nan
    wbal%errh2o(beg:end) = nan

  end subroutine init_water_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_carbon_balance_type
!
! !INTERFACE:
  subroutine init_carbon_balance_type(beg, end, cbal)
!
! !DESCRIPTION:
! Initialize carbon balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(carbon_balance_type), intent(inout):: cbal
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(cbal%begcb(beg:end))
    allocate(cbal%endcb(beg:end))
    allocate(cbal%errcb(beg:end))

    cbal%begcb(beg:end) = nan
    cbal%endcb(beg:end) = nan
    cbal%errcb(beg:end) = nan

  end subroutine init_carbon_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_nitrogen_balance_type
!
! !INTERFACE:
  subroutine init_nitrogen_balance_type(beg, end, nbal)
!
! !DESCRIPTION:
! Initialize nitrogen balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(nitrogen_balance_type), intent(inout):: nbal
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(nbal%begnb(beg:end))
    allocate(nbal%endnb(beg:end))
    allocate(nbal%errnb(beg:end))

    nbal%begnb(beg:end) = nan
    nbal%endnb(beg:end) = nan
    nbal%errnb(beg:end) = nan

  end subroutine init_nitrogen_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_ecophys_constants
!
! !INTERFACE:
  subroutine init_pft_ecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pftcon%ncorn(0:numpft))
    allocate(pftcon%nwheat(0:numpft))
    allocate(pftcon%noveg(0:numpft))
    allocate(pftcon%ntree(0:numpft))
    allocate(pftcon%smpso(0:numpft)) 
    allocate(pftcon%smpsc(0:numpft)) 
    allocate(pftcon%fnitr(0:numpft))
    allocate(pftcon%foln(0:numpft))
    allocate(pftcon%dleaf(0:numpft))
    allocate(pftcon%c3psn(0:numpft))
    allocate(pftcon%vcmx25(0:numpft))
    allocate(pftcon%mp(0:numpft))
    allocate(pftcon%qe25(0:numpft))
    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft,numrad))
    allocate(pftcon%rhos(0:numpft,numrad))
    allocate(pftcon%taul(0:numpft,numrad))
    allocate(pftcon%taus(0:numpft,numrad))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%roota_par(0:numpft))
    allocate(pftcon%rootb_par(0:numpft))
    allocate(pftcon%sla(0:numpft))
    allocate(pftcon%slatop(0:numpft))
    allocate(pftcon%dsladlai(0:numpft))
    allocate(pftcon%leafcn(0:numpft))
    allocate(pftcon%flnr(0:numpft))
    allocate(pftcon%woody(0:numpft))
    allocate(pftcon%lflitcn(0:numpft))
    allocate(pftcon%frootcn(0:numpft))
    allocate(pftcon%livewdcn(0:numpft))
    allocate(pftcon%deadwdcn(0:numpft))
    allocate(pftcon%froot_leaf(0:numpft))
    allocate(pftcon%stem_leaf(0:numpft))
    allocate(pftcon%croot_stem(0:numpft))
    allocate(pftcon%flivewd(0:numpft))
    allocate(pftcon%fcur(0:numpft))
    allocate(pftcon%lf_flab(0:numpft))
    allocate(pftcon%lf_fcel(0:numpft))
    allocate(pftcon%lf_flig(0:numpft))
    allocate(pftcon%fr_flab(0:numpft))
    allocate(pftcon%fr_fcel(0:numpft))
    allocate(pftcon%fr_flig(0:numpft))
    allocate(pftcon%dw_fcel(0:numpft))
    allocate(pftcon%dw_flig(0:numpft))
    allocate(pftcon%leaf_long(0:numpft))
    allocate(pftcon%evergreen(0:numpft))
    allocate(pftcon%stress_decid(0:numpft))
    allocate(pftcon%season_decid(0:numpft))
    allocate(pftcon%resist(0:numpft))

    pftcon%ncorn(:) = bigint
    pftcon%nwheat(:) = bigint
    pftcon%noveg(:) = bigint
    pftcon%ntree(:) = bigint
    pftcon%smpso(:) = nan
    pftcon%smpsc(:) = nan
    pftcon%fnitr(:) = nan
    pftcon%foln(:) = nan
    pftcon%dleaf(:) = nan
    pftcon%c3psn(:) = nan
    pftcon%vcmx25(:) = nan
    pftcon%mp(:) = nan
    pftcon%qe25(:) = nan
    pftcon%xl(:) = nan
    pftcon%rhol(:,:numrad) = nan
    pftcon%rhos(:,:numrad) = nan
    pftcon%taul(:,:numrad) = nan
    pftcon%taus(:,:numrad) = nan
    pftcon%z0mr(:) = nan
    pftcon%displar(:) = nan
    pftcon%roota_par(:) = nan
    pftcon%rootb_par(:) = nan
    pftcon%sla(:) = nan
    pftcon%slatop(:) = nan
    pftcon%dsladlai(:) = nan
    pftcon%leafcn(:) = nan
    pftcon%flnr(:) = nan
    pftcon%woody(:) = nan
    pftcon%lflitcn(:) = nan
    pftcon%frootcn(:) = nan
    pftcon%livewdcn(:) = nan
    pftcon%deadwdcn(:) = nan
    pftcon%froot_leaf(:) = nan
    pftcon%stem_leaf(:) = nan
    pftcon%croot_stem(:) = nan
    pftcon%flivewd(:) = nan
    pftcon%fcur(:) = nan
    pftcon%lf_flab(:) = nan
    pftcon%lf_fcel(:) = nan
    pftcon%lf_flig(:) = nan
    pftcon%fr_flab(:) = nan
    pftcon%fr_fcel(:) = nan
    pftcon%fr_flig(:) = nan
    pftcon%dw_fcel(:) = nan
    pftcon%dw_flig(:) = nan
    pftcon%leaf_long(:) = nan
    pftcon%evergreen(:) = nan
    pftcon%stress_decid(:) = nan
    pftcon%season_decid(:) = nan
    pftcon%resist(:) = nan

  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_DGVMecophys_constants
!
! !INTERFACE:
  subroutine init_pft_DGVMecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(dgv_pftcon%respcoeff(0:numpft))
    allocate(dgv_pftcon%flam(0:numpft))
    allocate(dgv_pftcon%resist(0:numpft))
    allocate(dgv_pftcon%l_turn(0:numpft))
    allocate(dgv_pftcon%l_long(0:numpft))
    allocate(dgv_pftcon%s_turn(0:numpft))
    allocate(dgv_pftcon%r_turn(0:numpft))
    allocate(dgv_pftcon%l_cton(0:numpft))
    allocate(dgv_pftcon%s_cton(0:numpft))
    allocate(dgv_pftcon%r_cton(0:numpft))
    allocate(dgv_pftcon%l_morph(0:numpft))
    allocate(dgv_pftcon%l_phen(0:numpft))
    allocate(dgv_pftcon%lmtorm(0:numpft))
    allocate(dgv_pftcon%crownarea_max(0:numpft))
    allocate(dgv_pftcon%init_lai(0:numpft))
    allocate(dgv_pftcon%x(0:numpft))
    allocate(dgv_pftcon%tcmin(0:numpft))
    allocate(dgv_pftcon%tcmax(0:numpft))
    allocate(dgv_pftcon%gddmin(0:numpft))
    allocate(dgv_pftcon%twmax(0:numpft))
    allocate(dgv_pftcon%lm_sapl(0:numpft))
    allocate(dgv_pftcon%sm_sapl(0:numpft))
    allocate(dgv_pftcon%hm_sapl(0:numpft))
    allocate(dgv_pftcon%rm_sapl(0:numpft))
    allocate(dgv_pftcon%tree(0:numpft))
    allocate(dgv_pftcon%summergreen(0:numpft))
    allocate(dgv_pftcon%raingreen(0:numpft))
    allocate(dgv_pftcon%reinickerp(0:numpft))
    allocate(dgv_pftcon%wooddens(0:numpft))
    allocate(dgv_pftcon%latosa(0:numpft))
    allocate(dgv_pftcon%allom1(0:numpft))
    allocate(dgv_pftcon%allom2(0:numpft))
    allocate(dgv_pftcon%allom3(0:numpft))

    dgv_pftcon%respcoeff(:) = nan
    dgv_pftcon%flam(:) = nan
    dgv_pftcon%resist(:) = nan
    dgv_pftcon%l_turn(:) = nan
    dgv_pftcon%l_long(:) = nan
    dgv_pftcon%s_turn(:) = nan
    dgv_pftcon%r_turn(:) = nan
    dgv_pftcon%l_cton(:) = nan
    dgv_pftcon%s_cton(:) = nan
    dgv_pftcon%r_cton(:) = nan
    dgv_pftcon%l_morph(:) = nan
    dgv_pftcon%l_phen(:) = nan
    dgv_pftcon%lmtorm(:) = nan
    dgv_pftcon%crownarea_max(:) = nan
    dgv_pftcon%init_lai(:) = nan
    dgv_pftcon%x(:) = nan
    dgv_pftcon%tcmin(:) = nan
    dgv_pftcon%tcmax(:) = nan
    dgv_pftcon%gddmin(:) = nan
    dgv_pftcon%twmax(:) = nan
    dgv_pftcon%lm_sapl(:) = nan
    dgv_pftcon%sm_sapl(:) = nan
    dgv_pftcon%hm_sapl(:) = nan
    dgv_pftcon%rm_sapl(:) = nan
    dgv_pftcon%tree(:) = .false.
    dgv_pftcon%summergreen(:) = .false.
    dgv_pftcon%raingreen(:) = .false.
    dgv_pftcon%reinickerp(:) = nan
    dgv_pftcon%wooddens(:) = nan
    dgv_pftcon%latosa(:) = nan
    dgv_pftcon%allom1(:) = nan
    dgv_pftcon%allom2(:) = nan
    dgv_pftcon%allom3(:) = nan

  end subroutine init_pft_DGVMecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pstate_type
!
! !INTERFACE:
  subroutine init_pft_pstate_type(beg, end, pps)
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !USES:
#if (defined CASA)
    use CASAMod   , only : npools, nresp_pools, nlive
    use clm_varcon, only : spval
#endif
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_pstate_type), intent(inout):: pps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pps%frac_veg_nosno(beg:end))
    allocate(pps%frac_veg_nosno_alb(beg:end))
    allocate(pps%emv(beg:end))
    allocate(pps%z0mv(beg:end))
    allocate(pps%z0hv(beg:end))
    allocate(pps%z0qv(beg:end))
    allocate(pps%rootfr(beg:end,1:nlevsoi))
    allocate(pps%rootr(beg:end,1:nlevsoi))
    allocate(pps%rresis(beg:end,1:nlevsoi))
    allocate(pps%dewmx(beg:end))
    allocate(pps%rssun(beg:end))
    allocate(pps%rssha(beg:end))
    allocate(pps%laisun(beg:end))
    allocate(pps%laisha(beg:end))
    allocate(pps%btran(beg:end))
    allocate(pps%fsun(beg:end))
    allocate(pps%tlai(beg:end))
    allocate(pps%tsai(beg:end))
    allocate(pps%elai(beg:end))
    allocate(pps%esai(beg:end))
    allocate(pps%fwet(beg:end))
    allocate(pps%fdry(beg:end))
    allocate(pps%dt_veg(beg:end))
    allocate(pps%htop(beg:end))
    allocate(pps%hbot(beg:end))
    allocate(pps%z0m(beg:end))
    allocate(pps%displa(beg:end))
    allocate(pps%albd(beg:end,1:numrad))
    allocate(pps%albi(beg:end,1:numrad))
    allocate(pps%fabd(beg:end,1:numrad))
    allocate(pps%fabi(beg:end,1:numrad))
    allocate(pps%ftdd(beg:end,1:numrad))
    allocate(pps%ftid(beg:end,1:numrad))
    allocate(pps%ftii(beg:end,1:numrad))
    allocate(pps%u10(beg:end))
    allocate(pps%fv(beg:end))
    allocate(pps%ram1(beg:end))
    allocate(pps%slasun(beg:end))
    allocate(pps%slasha(beg:end))
    allocate(pps%lncsun(beg:end))
    allocate(pps%lncsha(beg:end))
    allocate(pps%vcmxsun(beg:end))
    allocate(pps%vcmxsha(beg:end))
    allocate(pps%gdir(beg:end))
    allocate(pps%omega(beg:end,1:numrad))
    allocate(pps%eff_kid(beg:end,1:numrad))
    allocate(pps%eff_kii(beg:end,1:numrad))
    allocate(pps%sun_faid(beg:end,1:numrad))
    allocate(pps%sun_faii(beg:end,1:numrad))
    allocate(pps%sha_faid(beg:end,1:numrad))
    allocate(pps%sha_faii(beg:end,1:numrad))
    ! 4/14/05: PET
    ! Adding isotope code
    allocate(pps%cisun(beg:end))
    allocate(pps%cisha(beg:end))
    allocate(pps%alphapsnsun(beg:end))
    allocate(pps%alphapsnsha(beg:end))
    
#if (defined CASA)
    allocate(pps%Closs(beg:end,npools))  ! C lost to atm
    allocate(pps%Resp_C(beg:end,npools))
    allocate(pps%Tpool_C(beg:end,npools))! Total C pool size
    allocate(pps%eff(beg:end,nresp_pools))
    allocate(pps%frac_donor(beg:end,nresp_pools))
    allocate(pps%livefr(beg:end,nlive))  !live fraction
    allocate(pps%pet(beg:end))           !potential evaporation (mm h2o/s)
    allocate(pps%co2flux(beg:end))       ! net CO2 flux (g C/m2/sec) [+= atm]
    allocate(pps%fnpp(beg:end))          ! NPP  (g C/m2/sec)
    allocate(pps%soilt(beg:end))         !soil temp for top 30cm
    allocate(pps%smoist(beg:end))        !soil moisture for top 30cm
    allocate(pps%sz(beg:end))
    allocate(pps%watopt(beg:end))
    allocate(pps%watdry(beg:end))
    allocate(pps%soiltc(beg:end))         !soil temp for entire column
    allocate(pps%smoistc(beg:end))        !soil moisture for entire column
    allocate(pps%szc(beg:end))
    allocate(pps%watoptc(beg:end))
    allocate(pps%watdryc(beg:end))
    allocate(pps%Wlim(beg:end))
    allocate(pps%litterscalar(beg:end))
    allocate(pps%rootlitscalar(beg:end))
    allocate(pps%stressCD(beg:end))
    allocate(pps%excessC(beg:end))       ! excess Carbon (gC/m2/timestep)
    allocate(pps%bgtemp(beg:end))
    allocate(pps%bgmoist(beg:end))
    allocate(pps%plai(beg:end))          ! prognostic LAI (m2 leaf/m2 ground)
    allocate(pps%Cflux(beg:end))
    allocate(pps%XSCpool(beg:end))
    allocate(pps%tday(beg:end))     ! daily accumulated temperature (deg C)
    allocate(pps%tdayavg(beg:end))  ! daily averaged temperature (deg C)
    allocate(pps%tcount(beg:end))   ! counter for daily avg temp
    allocate(pps%degday(beg:end))   ! accumulated degree days (deg C)
    allocate(pps%ndegday(beg:end))  ! counter for number of degree days
    allocate(pps%stressT(beg:end))
    allocate(pps%stressW(beg:end))  ! water stress function for leaf loss
    allocate(pps%iseabeg(beg:end))  ! index for start of growing season
    allocate(pps%nstepbeg(beg:end)) ! nstep at start of growing season
    allocate(pps%lgrow(beg:end))    ! growing season index (0 or 1) to be
                                    ! passed daily to CASA to get NPP
    allocate(pps%sandfrac(beg:end))
    allocate(pps%clayfrac(beg:end))
#endif

    pps%frac_veg_nosno(beg:end) = bigint
    pps%frac_veg_nosno_alb(beg:end) = 0
    pps%emv(beg:end) = nan
    pps%z0mv(beg:end) = nan
    pps%z0hv(beg:end) = nan
    pps%z0qv(beg:end) = nan
    pps%rootfr(beg:end,:nlevsoi) = nan
    pps%rootr (beg:end,:nlevsoi) = nan
    pps%rresis(beg:end,:nlevsoi) = nan
    pps%dewmx(beg:end) = nan
    pps%rssun(beg:end) = nan
    pps%rssha(beg:end) = nan
    pps%laisun(beg:end) = nan
    pps%laisha(beg:end) = nan
    pps%btran(beg:end) = nan
    pps%fsun(beg:end) = nan
    pps%tlai(beg:end) = 0._r8
    pps%tsai(beg:end) = 0._r8
    pps%elai(beg:end) = 0._r8
    pps%esai(beg:end) = 0._r8
    pps%fwet(beg:end) = nan
    pps%fdry(beg:end) = nan
    pps%dt_veg(beg:end) = nan
    pps%htop(beg:end) = 0._r8
    pps%hbot(beg:end) = 0._r8
    pps%z0m(beg:end) = nan
    pps%displa(beg:end) = nan
    pps%albd(beg:end,:numrad) = nan
    pps%albi(beg:end,:numrad) = nan
    pps%fabd(beg:end,:numrad) = nan
    pps%fabi(beg:end,:numrad) = nan
    pps%ftdd(beg:end,:numrad) = nan
    pps%ftid(beg:end,:numrad) = nan
    pps%ftii(beg:end,:numrad) = nan
    pps%u10(beg:end) = nan
    pps%fv(beg:end) = nan
    pps%ram1(beg:end) = nan
    pps%slasun(beg:end) = nan
    pps%slasha(beg:end) = nan
    pps%lncsun(beg:end) = nan
    pps%lncsha(beg:end) = nan
    pps%vcmxsun(beg:end) = nan
    pps%vcmxsha(beg:end) = nan
    pps%gdir(beg:end) = nan
    pps%omega(beg:end,1:numrad) = nan
    pps%eff_kid(beg:end,1:numrad) = nan
    pps%eff_kii(beg:end,1:numrad) = nan
    pps%sun_faid(beg:end,1:numrad) = nan
    pps%sun_faii(beg:end,1:numrad) = nan
    pps%sha_faid(beg:end,1:numrad) = nan
    pps%sha_faii(beg:end,1:numrad) = nan
    ! 4/14/05: PET
    ! Adding isotope code
    pps%cisun(beg:end) = nan
    pps%cisha(beg:end) = nan
    pps%alphapsnsun(beg:end) = nan
    pps%alphapsnsha(beg:end) = nan

#if (defined CASA)
    pps%Closs(beg:end,:npools) = spval   !init w/ spval the variables that
    pps%Resp_C(beg:end,:npools) = nan    !go to history, because CASA
    pps%Tpool_C(beg:end,:npools) = spval !routines do not get called on
    pps%livefr(beg:end,:nlive) = spval   !first timestep of nsrest=0 and
    pps%pet(beg:end) = spval             !history would get nans
    pps%co2flux(beg:end) = nan           !in the first timestep
    pps%fnpp(beg:end) = nan
    pps%excessC(beg:end) = spval
    pps%bgtemp(beg:end) = spval
    pps%bgmoist(beg:end) = spval
    pps%plai(beg:end) = spval
    pps%Cflux(beg:end) = nan
    pps%XSCpool(beg:end) = spval
    pps%tdayavg(beg:end) = spval
    pps%degday(beg:end) = spval
    pps%stressT(beg:end) = spval
    pps%stressW(beg:end) = spval
    pps%stressCD(beg:end) = spval
    pps%iseabeg(beg:end) = spval
    pps%nstepbeg(beg:end) = spval
    pps%lgrow(beg:end) = spval
    pps%eff(beg:end,:nresp_pools) = nan
    pps%frac_donor(beg:end,:nresp_pools) = nan
    pps%soilt(beg:end) = spval                  ! on history file
    pps%smoist(beg:end) = spval                 ! on history file
    pps%sz(beg:end) = nan
    pps%watopt(beg:end) = nan
    pps%watdry(beg:end) = nan
    pps%soiltc(beg:end) = nan
    pps%smoistc(beg:end) = nan
    pps%szc(beg:end) = nan
    pps%watoptc(beg:end) = spval                ! on history file
    pps%watdryc(beg:end) = spval                ! on history file
    pps%Wlim(beg:end) = spval                   ! on history file
    pps%litterscalar(beg:end) = nan
    pps%rootlitscalar(beg:end) = nan
    pps%tday(beg:end) = nan
    pps%tcount(beg:end) = nan
    pps%ndegday(beg:end) = nan
    pps%sandfrac(beg:end) = nan
    pps%clayfrac(beg:end) = nan
#endif

  end subroutine init_pft_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_epv_type
!
! !INTERFACE:
  subroutine init_pft_epv_type(beg, end, pepv)
!
! !DESCRIPTION:
! Initialize pft ecophysiological variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_epv_type), intent(inout):: pepv
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pepv%dormant_flag(beg:end))
    allocate(pepv%days_active(beg:end))
    allocate(pepv%onset_flag(beg:end))
    allocate(pepv%onset_counter(beg:end))
    allocate(pepv%onset_gddflag(beg:end))
    allocate(pepv%onset_fdd(beg:end))
    allocate(pepv%onset_gdd(beg:end))
    allocate(pepv%onset_swi(beg:end))
    allocate(pepv%offset_flag(beg:end))
    allocate(pepv%offset_counter(beg:end))
    allocate(pepv%offset_fdd(beg:end))
    allocate(pepv%offset_swi(beg:end))
    allocate(pepv%lgsf(beg:end))
    allocate(pepv%bglfr(beg:end))
    allocate(pepv%bgtr(beg:end))
    allocate(pepv%dayl(beg:end))
    allocate(pepv%prev_dayl(beg:end))
    allocate(pepv%annavg_t2m(beg:end))
    allocate(pepv%tempavg_t2m(beg:end))
    allocate(pepv%gpp(beg:end))
    allocate(pepv%availc(beg:end))
    allocate(pepv%xsmrpool_recover(beg:end))
    allocate(pepv%xsmrpool_c13ratio(beg:end))
    allocate(pepv%alloc_pnow(beg:end))
    allocate(pepv%c_allometry(beg:end))
    allocate(pepv%n_allometry(beg:end))
    allocate(pepv%plant_ndemand(beg:end))
    allocate(pepv%tempsum_plant_ndemand(beg:end))
    allocate(pepv%annsum_plant_ndemand(beg:end))
    allocate(pepv%tempsum_retransn(beg:end))
    allocate(pepv%annsum_retransn(beg:end))
    allocate(pepv%avail_retransn(beg:end))
    allocate(pepv%plant_nalloc(beg:end))
    allocate(pepv%plant_calloc(beg:end))
    allocate(pepv%excess_cflux(beg:end))
    allocate(pepv%downreg(beg:end))
    allocate(pepv%prev_leafc_to_litter(beg:end))
    allocate(pepv%prev_frootc_to_litter(beg:end))
    allocate(pepv%tempsum_npp(beg:end))
    allocate(pepv%annsum_npp(beg:end))
    ! 4/21/05, PET
    ! Adding isotope code
	 allocate(pepv%rc13_canair(beg:end))
	 allocate(pepv%rc13_psnsun(beg:end))
	 allocate(pepv%rc13_psnsha(beg:end))

    pepv%dormant_flag(beg:end) = nan
    pepv%days_active(beg:end) = nan
    pepv%onset_flag(beg:end) = nan
    pepv%onset_counter(beg:end) = nan
    pepv%onset_gddflag(beg:end) = nan
    pepv%onset_fdd(beg:end) = nan
    pepv%onset_gdd(beg:end) = nan
    pepv%onset_swi(beg:end) = nan
    pepv%offset_flag(beg:end) = nan
    pepv%offset_counter(beg:end) = nan
    pepv%offset_fdd(beg:end) = nan
    pepv%offset_swi(beg:end) = nan
    pepv%lgsf(beg:end) = nan
    pepv%bglfr(beg:end) = nan
    pepv%bgtr(beg:end) = nan
    pepv%dayl(beg:end) = nan
    pepv%prev_dayl(beg:end) = nan
    pepv%annavg_t2m(beg:end) = nan
    pepv%tempavg_t2m(beg:end) = nan
    pepv%gpp(beg:end) = nan
    pepv%availc(beg:end) = nan
    pepv%xsmrpool_recover(beg:end) = nan
    pepv%xsmrpool_c13ratio(beg:end) = nan
    pepv%alloc_pnow(beg:end) = nan
    pepv%c_allometry(beg:end) = nan
    pepv%n_allometry(beg:end) = nan
    pepv%plant_ndemand(beg:end) = nan
    pepv%tempsum_plant_ndemand(beg:end) = nan
    pepv%annsum_plant_ndemand(beg:end) = nan
    pepv%tempsum_retransn(beg:end) = nan
    pepv%annsum_retransn(beg:end) = nan
    pepv%avail_retransn(beg:end) = nan
    pepv%plant_nalloc(beg:end) = nan
    pepv%plant_calloc(beg:end) = nan
    pepv%excess_cflux(beg:end) = nan
    pepv%downreg(beg:end) = nan
    pepv%prev_leafc_to_litter(beg:end) = nan
    pepv%prev_frootc_to_litter(beg:end) = nan
    pepv%tempsum_npp(beg:end) = nan
    pepv%annsum_npp(beg:end) = nan
    ! 4/21/05, PET
    ! Adding isotope code
    pepv%rc13_canair(beg:end) = nan
    pepv%rc13_psnsun(beg:end) = nan
    pepv%rc13_psnsha(beg:end) = nan
    
  end subroutine init_pft_epv_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pdgvstate_type
!
! !INTERFACE:
  subroutine init_pft_pdgvstate_type(beg, end, pdgvs)
!
! !DESCRIPTION:
! Initialize pft DGVM state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dgvstate_type), intent(inout):: pdgvs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdgvs%agdd0(beg:end))
    allocate(pdgvs%agdd5(beg:end))
    allocate(pdgvs%agddtw(beg:end))
    allocate(pdgvs%agdd(beg:end))
    allocate(pdgvs%t10(beg:end))
    allocate(pdgvs%t_mo(beg:end))
    allocate(pdgvs%t_mo_min(beg:end))
    allocate(pdgvs%fnpsn10(beg:end))
    allocate(pdgvs%prec365(beg:end))
    allocate(pdgvs%agdd20(beg:end))
    allocate(pdgvs%tmomin20(beg:end))
    allocate(pdgvs%t10min(beg:end))
    allocate(pdgvs%tsoi25(beg:end))
    allocate(pdgvs%annpsn(beg:end))
    allocate(pdgvs%annpsnpot(beg:end))
    allocate(pdgvs%present(beg:end))
    allocate(pdgvs%dphen(beg:end))
    allocate(pdgvs%leafon(beg:end))
    allocate(pdgvs%leafof(beg:end))
    allocate(pdgvs%nind(beg:end))
    allocate(pdgvs%lm_ind(beg:end))
    allocate(pdgvs%sm_ind(beg:end))
    allocate(pdgvs%hm_ind(beg:end))
    allocate(pdgvs%rm_ind(beg:end))
    allocate(pdgvs%lai_ind(beg:end))
    allocate(pdgvs%fpcinc(beg:end))
    allocate(pdgvs%fpcgrid(beg:end))
    allocate(pdgvs%crownarea(beg:end))
    allocate(pdgvs%bm_inc(beg:end))
    allocate(pdgvs%afmicr(beg:end))
    allocate(pdgvs%firelength (beg:end))
    allocate(pdgvs%litterag(beg:end))
    allocate(pdgvs%litterbg(beg:end))
    allocate(pdgvs%cpool_fast(beg:end))
    allocate(pdgvs%cpool_slow(beg:end))
    allocate(pdgvs%k_fast_ave(beg:end))
    allocate(pdgvs%k_slow_ave(beg:end))
    allocate(pdgvs%litter_decom_ave(beg:end))
    allocate(pdgvs%turnover_ind(beg:end))

    pdgvs%agdd0(beg:end) = nan
    pdgvs%agdd5(beg:end) = nan
    pdgvs%agddtw(beg:end) = nan
    pdgvs%agdd(beg:end) = nan
    pdgvs%t10(beg:end) = nan
    pdgvs%t_mo(beg:end) = nan
    pdgvs%t_mo_min(beg:end) = nan
    pdgvs%fnpsn10(beg:end) = nan
    pdgvs%prec365(beg:end) = nan
    pdgvs%agdd20(beg:end) = nan
    pdgvs%tmomin20(beg:end) = nan
    pdgvs%t10min(beg:end) = nan
    pdgvs%tsoi25(beg:end) = nan
    pdgvs%annpsn(beg:end) = nan
    pdgvs%annpsnpot(beg:end) = nan
    pdgvs%present(beg:end) = .false.
    pdgvs%dphen(beg:end) = nan
    pdgvs%leafon(beg:end) = nan
    pdgvs%leafof(beg:end) = nan
    pdgvs%nind(beg:end) = nan
    pdgvs%lm_ind(beg:end) = nan
    pdgvs%sm_ind(beg:end) = nan
    pdgvs%hm_ind(beg:end) = nan
    pdgvs%rm_ind(beg:end) = nan
    pdgvs%lai_ind(beg:end) = nan
    pdgvs%fpcinc(beg:end) = nan
    pdgvs%fpcgrid(beg:end) = nan
    pdgvs%crownarea(beg:end) = nan
    pdgvs%bm_inc(beg:end) = nan
    pdgvs%afmicr(beg:end) = nan
    pdgvs%firelength (beg:end) = nan
    pdgvs%litterag(beg:end) = nan
    pdgvs%litterbg(beg:end) = nan
    pdgvs%cpool_fast(beg:end) = nan
    pdgvs%cpool_slow(beg:end) = nan
    pdgvs%k_fast_ave(beg:end) = nan
    pdgvs%k_slow_ave(beg:end) = nan
    pdgvs%litter_decom_ave(beg:end) = nan
    pdgvs%turnover_ind(beg:end) = nan

  end subroutine init_pft_pdgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_estate_type
!
! !INTERFACE:
  subroutine init_pft_estate_type(beg, end, pes)
!
! !DESCRIPTION:
! Initialize pft energy state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_estate_type), intent(inout):: pes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pes%t_ref2m(beg:end))
    allocate(pes%t_ref2m_min(beg:end))
    allocate(pes%t_ref2m_max(beg:end))
    allocate(pes%t_ref2m_min_inst(beg:end))
    allocate(pes%t_ref2m_max_inst(beg:end))
    allocate(pes%q_ref2m(beg:end))
    allocate(pes%t_veg(beg:end))

    pes%t_ref2m(beg:end) = nan
    pes%t_ref2m_min(beg:end) = nan
    pes%t_ref2m_max(beg:end) = nan
    pes%t_ref2m_min_inst(beg:end) = nan
    pes%t_ref2m_max_inst(beg:end) = nan
    pes%q_ref2m(beg:end) = nan
    pes%t_veg(beg:end) = nan

  end subroutine init_pft_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wstate_type
!
! !INTERFACE:
  subroutine init_pft_wstate_type(beg, end, pws)
!
! !DESCRIPTION:
! Initialize pft water state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wstate_type), intent(inout):: pws !pft water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pws%h2ocan(beg:end))
    pws%h2ocan(beg:end) = nan

  end subroutine init_pft_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cstate_type
!
! !INTERFACE:
  subroutine init_pft_cstate_type(beg, end, pcs)
!
! !DESCRIPTION:
! Initialize pft carbon state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cstate_type), intent(inout):: pcs !pft carbon state
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pcs%leafc(beg:end))
    allocate(pcs%leafc_storage(beg:end))
    allocate(pcs%leafc_xfer(beg:end))
    allocate(pcs%frootc(beg:end))
    allocate(pcs%frootc_storage(beg:end))
    allocate(pcs%frootc_xfer(beg:end))
    allocate(pcs%livestemc(beg:end))
    allocate(pcs%livestemc_storage(beg:end))
    allocate(pcs%livestemc_xfer(beg:end))
    allocate(pcs%deadstemc(beg:end))
    allocate(pcs%deadstemc_storage(beg:end))
    allocate(pcs%deadstemc_xfer(beg:end))
    allocate(pcs%livecrootc(beg:end))
    allocate(pcs%livecrootc_storage(beg:end))
    allocate(pcs%livecrootc_xfer(beg:end))
    allocate(pcs%deadcrootc(beg:end))
    allocate(pcs%deadcrootc_storage(beg:end))
    allocate(pcs%deadcrootc_xfer(beg:end))
    allocate(pcs%gresp_storage(beg:end))
    allocate(pcs%gresp_xfer(beg:end))
    allocate(pcs%cpool(beg:end))
    allocate(pcs%xsmrpool(beg:end))
    allocate(pcs%pft_ctrunc(beg:end))
    allocate(pcs%dispvegc(beg:end))
    allocate(pcs%storvegc(beg:end))
    allocate(pcs%totvegc(beg:end))
    allocate(pcs%totpftc(beg:end))

    pcs%leafc(beg:end) = nan
    pcs%leafc_storage(beg:end) = nan
    pcs%leafc_xfer(beg:end) = nan
    pcs%frootc(beg:end) = nan
    pcs%frootc_storage(beg:end) = nan
    pcs%frootc_xfer(beg:end) = nan
    pcs%livestemc(beg:end) = nan
    pcs%livestemc_storage(beg:end) = nan
    pcs%livestemc_xfer(beg:end) = nan
    pcs%deadstemc(beg:end) = nan
    pcs%deadstemc_storage(beg:end) = nan
    pcs%deadstemc_xfer(beg:end) = nan
    pcs%livecrootc(beg:end) = nan
    pcs%livecrootc_storage(beg:end) = nan
    pcs%livecrootc_xfer(beg:end) = nan
    pcs%deadcrootc(beg:end) = nan
    pcs%deadcrootc_storage(beg:end) = nan
    pcs%deadcrootc_xfer(beg:end) = nan
    pcs%gresp_storage(beg:end) = nan
    pcs%gresp_xfer(beg:end) = nan
    pcs%cpool(beg:end) = nan
    pcs%xsmrpool(beg:end) = nan
    pcs%pft_ctrunc(beg:end) = nan
    pcs%dispvegc(beg:end) = nan
    pcs%storvegc(beg:end) = nan
    pcs%totvegc(beg:end) = nan
    pcs%totpftc(beg:end) = nan

  end subroutine init_pft_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nstate_type
!
! !INTERFACE:
  subroutine init_pft_nstate_type(beg, end, pns)
!
! !DESCRIPTION:
! Initialize pft nitrogen state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nstate_type), intent(inout):: pns !pft nitrogen state
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pns%leafn(beg:end))
    allocate(pns%leafn_storage(beg:end))
    allocate(pns%leafn_xfer(beg:end))
    allocate(pns%frootn(beg:end))
    allocate(pns%frootn_storage(beg:end))
    allocate(pns%frootn_xfer(beg:end))
    allocate(pns%livestemn(beg:end))
    allocate(pns%livestemn_storage(beg:end))
    allocate(pns%livestemn_xfer(beg:end))
    allocate(pns%deadstemn(beg:end))
    allocate(pns%deadstemn_storage(beg:end))
    allocate(pns%deadstemn_xfer(beg:end))
    allocate(pns%livecrootn(beg:end))
    allocate(pns%livecrootn_storage(beg:end))
    allocate(pns%livecrootn_xfer(beg:end))
    allocate(pns%deadcrootn(beg:end))
    allocate(pns%deadcrootn_storage(beg:end))
    allocate(pns%deadcrootn_xfer(beg:end))
    allocate(pns%retransn(beg:end))
    allocate(pns%npool(beg:end))
    allocate(pns%pft_ntrunc(beg:end))
    allocate(pns%dispvegn(beg:end))
    allocate(pns%storvegn(beg:end))
    allocate(pns%totvegn(beg:end))
    allocate(pns%totpftn(beg:end))

    pns%leafn(beg:end) = nan
    pns%leafn_storage(beg:end) = nan
    pns%leafn_xfer(beg:end) = nan
    pns%frootn(beg:end) = nan
    pns%frootn_storage(beg:end) = nan
    pns%frootn_xfer(beg:end) = nan
    pns%livestemn(beg:end) = nan
    pns%livestemn_storage(beg:end) = nan
    pns%livestemn_xfer(beg:end) = nan
    pns%deadstemn(beg:end) = nan
    pns%deadstemn_storage(beg:end) = nan
    pns%deadstemn_xfer(beg:end) = nan
    pns%livecrootn(beg:end) = nan
    pns%livecrootn_storage(beg:end) = nan
    pns%livecrootn_xfer(beg:end) = nan
    pns%deadcrootn(beg:end) = nan
    pns%deadcrootn_storage(beg:end) = nan
    pns%deadcrootn_xfer(beg:end) = nan
    pns%retransn(beg:end) = nan
    pns%npool(beg:end) = nan
    pns%pft_ntrunc(beg:end) = nan
    pns%dispvegn(beg:end) = nan
    pns%storvegn(beg:end) = nan
    pns%totvegn(beg:end) = nan
    pns%totpftn(beg:end) = nan

  end subroutine init_pft_nstate_type
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_eflux_type
!
! !INTERFACE:
  subroutine init_pft_eflux_type(beg, end, pef)
!
! !DESCRIPTION:
! Initialize pft energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_eflux_type), intent(inout):: pef
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pef%sabg(beg:end))
    allocate(pef%sabv(beg:end))
    allocate(pef%fsa(beg:end))
    allocate(pef%fsr(beg:end))
    allocate(pef%parsun(beg:end))
    allocate(pef%parsha(beg:end))
    allocate(pef%dlrad(beg:end))
    allocate(pef%ulrad(beg:end))
    allocate(pef%eflx_lh_tot(beg:end))
    allocate(pef%eflx_lh_grnd(beg:end))
    allocate(pef%eflx_soil_grnd(beg:end))
    allocate(pef%eflx_sh_tot(beg:end))
    allocate(pef%eflx_sh_grnd(beg:end))
    allocate(pef%eflx_sh_veg(beg:end))
    allocate(pef%eflx_lh_vege(beg:end))
    allocate(pef%eflx_lh_vegt(beg:end))
    allocate(pef%cgrnd(beg:end))
    allocate(pef%cgrndl(beg:end))
    allocate(pef%cgrnds(beg:end))
    allocate(pef%eflx_gnet(beg:end))
    allocate(pef%dgnetdT(beg:end))
    allocate(pef%eflx_lwrad_out(beg:end))
    allocate(pef%eflx_lwrad_net(beg:end))
    allocate(pef%fsds_vis_d(beg:end))
    allocate(pef%fsds_nir_d(beg:end))
    allocate(pef%fsds_vis_i(beg:end))
    allocate(pef%fsds_nir_i(beg:end))
    allocate(pef%fsr_vis_d(beg:end))
    allocate(pef%fsr_nir_d(beg:end))
    allocate(pef%fsr_vis_i(beg:end))
    allocate(pef%fsr_nir_i(beg:end))
    allocate(pef%fsds_vis_d_ln(beg:end))
    allocate(pef%fsds_nir_d_ln(beg:end))
    allocate(pef%fsr_vis_d_ln(beg:end))
    allocate(pef%fsr_nir_d_ln(beg:end))
    allocate(pef%sun_add(beg:end,1:numrad))
    allocate(pef%tot_aid(beg:end,1:numrad))
    allocate(pef%sun_aid(beg:end,1:numrad))
    allocate(pef%sun_aii(beg:end,1:numrad))
    allocate(pef%sha_aid(beg:end,1:numrad))
    allocate(pef%sha_aii(beg:end,1:numrad))
    allocate(pef%sun_atot(beg:end,1:numrad))
    allocate(pef%sha_atot(beg:end,1:numrad))
    allocate(pef%sun_alf(beg:end,1:numrad))
    allocate(pef%sha_alf(beg:end,1:numrad))
    allocate(pef%sun_aperlai(beg:end,1:numrad))
    allocate(pef%sha_aperlai(beg:end,1:numrad))

    pef%sabg(beg:end) = nan
    pef%sabv(beg:end) = nan
    pef%fsa(beg:end) = nan
    pef%fsr(beg:end) = nan
    pef%parsun(beg:end) = nan
    pef%parsha(beg:end) = nan
    pef%dlrad(beg:end) = nan
    pef%ulrad(beg:end) = nan
    pef%eflx_lh_tot(beg:end) = nan
    pef%eflx_lh_grnd(beg:end) = nan
    pef%eflx_soil_grnd(beg:end) = nan
    pef%eflx_sh_tot(beg:end) = nan
    pef%eflx_sh_grnd(beg:end) = nan
    pef%eflx_sh_veg(beg:end) = nan
    pef%eflx_lh_vege(beg:end) = nan
    pef%eflx_lh_vegt(beg:end) = nan
    pef%cgrnd(beg:end) = nan
    pef%cgrndl(beg:end) = nan
    pef%cgrnds(beg:end) = nan
    pef%eflx_gnet(beg:end) = nan
    pef%dgnetdT(beg:end) = nan
    pef%eflx_lwrad_out(beg:end) = nan
    pef%eflx_lwrad_net(beg:end) = nan
    pef%fsds_vis_d(beg:end) = nan
    pef%fsds_nir_d(beg:end) = nan
    pef%fsds_vis_i(beg:end) = nan
    pef%fsds_nir_i(beg:end) = nan
    pef%fsr_vis_d(beg:end) = nan
    pef%fsr_nir_d(beg:end) = nan
    pef%fsr_vis_i(beg:end) = nan
    pef%fsr_nir_i(beg:end) = nan
    pef%fsds_vis_d_ln(beg:end) = nan
    pef%fsds_nir_d_ln(beg:end) = nan
    pef%fsr_vis_d_ln(beg:end) = nan
    pef%fsr_nir_d_ln(beg:end) = nan
    pef%sun_add(beg:end,1:numrad) = nan
    pef%tot_aid(beg:end,1:numrad) = nan
    pef%sun_aid(beg:end,1:numrad) = nan
    pef%sun_aii(beg:end,1:numrad) = nan
    pef%sha_aid(beg:end,1:numrad) = nan
    pef%sha_aii(beg:end,1:numrad) = nan
    pef%sun_atot(beg:end,1:numrad) = nan
    pef%sha_atot(beg:end,1:numrad) = nan
    pef%sun_alf(beg:end,1:numrad) = nan
    pef%sha_alf(beg:end,1:numrad) = nan
    pef%sun_aperlai(beg:end,1:numrad) = nan
    pef%sha_aperlai(beg:end,1:numrad) = nan

  end subroutine init_pft_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_mflux_type
!
! !INTERFACE:
  subroutine init_pft_mflux_type(beg, end, pmf)
!
! !DESCRIPTION:
! Initialize pft momentum flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_mflux_type), intent(inout) :: pmf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pmf%taux(beg:end))
    allocate(pmf%tauy(beg:end))

    pmf%taux(beg:end) = nan
    pmf%tauy(beg:end) = nan

  end subroutine init_pft_mflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wflux_type
!
! !INTERFACE:
  subroutine init_pft_wflux_type(beg, end, pwf)
!
! !DESCRIPTION:
! Initialize pft water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wflux_type), intent(inout) :: pwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pwf%qflx_prec_intr(beg:end))
    allocate(pwf%qflx_prec_grnd(beg:end))
    allocate(pwf%qflx_rain_grnd(beg:end))
    allocate(pwf%qflx_snow_grnd(beg:end))
    allocate(pwf%qflx_snowcap(beg:end))
    allocate(pwf%qflx_evap_veg(beg:end))
    allocate(pwf%qflx_tran_veg(beg:end))
    allocate(pwf%qflx_evap_can(beg:end))
    allocate(pwf%qflx_evap_soi(beg:end))
    allocate(pwf%qflx_evap_tot(beg:end))
    allocate(pwf%qflx_evap_grnd(beg:end))
    allocate(pwf%qflx_dew_grnd(beg:end))
    allocate(pwf%qflx_sub_snow(beg:end))
    allocate(pwf%qflx_dew_snow(beg:end))

    pwf%qflx_prec_intr(beg:end) = nan
    pwf%qflx_prec_grnd(beg:end) = nan
    pwf%qflx_rain_grnd(beg:end) = nan
    pwf%qflx_snow_grnd(beg:end) = nan
    pwf%qflx_snowcap(beg:end) = nan
    pwf%qflx_evap_veg(beg:end) = nan
    pwf%qflx_tran_veg(beg:end) = nan
    pwf%qflx_evap_can(beg:end) = nan
    pwf%qflx_evap_soi(beg:end) = nan
    pwf%qflx_evap_tot(beg:end) = nan
    pwf%qflx_evap_grnd(beg:end) = nan
    pwf%qflx_dew_grnd(beg:end) = nan
    pwf%qflx_sub_snow(beg:end) = nan
    pwf%qflx_dew_snow(beg:end) = nan

  end subroutine init_pft_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cflux_type
!
! !INTERFACE:
  subroutine init_pft_cflux_type(beg, end, pcf)
!
! !DESCRIPTION:
! Initialize pft carbon flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cflux_type), intent(inout) :: pcf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pcf%psnsun(beg:end))
    allocate(pcf%psnsha(beg:end))
    allocate(pcf%fpsn(beg:end))
    allocate(pcf%frm(beg:end))
    allocate(pcf%frmf(beg:end))
    allocate(pcf%frms(beg:end))
    allocate(pcf%frmr(beg:end))
    allocate(pcf%frg(beg:end))
    allocate(pcf%dmi(beg:end))
    allocate(pcf%fco2(beg:end))
    allocate(pcf%fmicr(beg:end))
    allocate(pcf%m_leafc_to_litter(beg:end))
    allocate(pcf%m_frootc_to_litter(beg:end))
    allocate(pcf%m_leafc_storage_to_litter(beg:end))
    allocate(pcf%m_frootc_storage_to_litter(beg:end))
    allocate(pcf%m_livestemc_storage_to_litter(beg:end))
    allocate(pcf%m_deadstemc_storage_to_litter(beg:end))
    allocate(pcf%m_livecrootc_storage_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_litter(beg:end))
    allocate(pcf%m_leafc_xfer_to_litter(beg:end))
    allocate(pcf%m_frootc_xfer_to_litter(beg:end))
    allocate(pcf%m_livestemc_xfer_to_litter(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_litter(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_litter(beg:end))
    allocate(pcf%m_livestemc_to_litter(beg:end))
    allocate(pcf%m_deadstemc_to_litter(beg:end))
    allocate(pcf%m_livecrootc_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_to_litter(beg:end))
    allocate(pcf%m_gresp_storage_to_litter(beg:end))
    allocate(pcf%m_gresp_xfer_to_litter(beg:end))
    allocate(pcf%m_leafc_to_fire(beg:end))
    allocate(pcf%m_frootc_to_fire(beg:end))
    allocate(pcf%m_leafc_storage_to_fire(beg:end))
    allocate(pcf%m_frootc_storage_to_fire(beg:end))
    allocate(pcf%m_livestemc_storage_to_fire(beg:end))
    allocate(pcf%m_deadstemc_storage_to_fire(beg:end))
    allocate(pcf%m_livecrootc_storage_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_fire(beg:end))
    allocate(pcf%m_leafc_xfer_to_fire(beg:end))
    allocate(pcf%m_frootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_fire(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_litter_fire(beg:end))
    allocate(pcf%m_gresp_storage_to_fire(beg:end))
    allocate(pcf%m_gresp_xfer_to_fire(beg:end))
    allocate(pcf%leafc_xfer_to_leafc(beg:end))
    allocate(pcf%frootc_xfer_to_frootc(beg:end))
    allocate(pcf%livestemc_xfer_to_livestemc(beg:end))
    allocate(pcf%deadstemc_xfer_to_deadstemc(beg:end))
    allocate(pcf%livecrootc_xfer_to_livecrootc(beg:end))
    allocate(pcf%deadcrootc_xfer_to_deadcrootc(beg:end))
    allocate(pcf%leafc_to_litter(beg:end))
    allocate(pcf%frootc_to_litter(beg:end))
    allocate(pcf%leaf_mr(beg:end))
    allocate(pcf%froot_mr(beg:end))
    allocate(pcf%livestem_mr(beg:end))
    allocate(pcf%livecroot_mr(beg:end))
    allocate(pcf%leaf_curmr(beg:end))
    allocate(pcf%froot_curmr(beg:end))
    allocate(pcf%livestem_curmr(beg:end))
    allocate(pcf%livecroot_curmr(beg:end))
    allocate(pcf%leaf_xsmr(beg:end))
    allocate(pcf%froot_xsmr(beg:end))
    allocate(pcf%livestem_xsmr(beg:end))
    allocate(pcf%livecroot_xsmr(beg:end))
    allocate(pcf%psnsun_to_cpool(beg:end))
    allocate(pcf%psnshade_to_cpool(beg:end))
    allocate(pcf%cpool_to_xsmrpool(beg:end))
    allocate(pcf%cpool_to_leafc(beg:end))
    allocate(pcf%cpool_to_leafc_storage(beg:end))
    allocate(pcf%cpool_to_frootc(beg:end))
    allocate(pcf%cpool_to_frootc_storage(beg:end))
    allocate(pcf%cpool_to_livestemc(beg:end))
    allocate(pcf%cpool_to_livestemc_storage(beg:end))
    allocate(pcf%cpool_to_deadstemc(beg:end))
    allocate(pcf%cpool_to_deadstemc_storage(beg:end))
    allocate(pcf%cpool_to_livecrootc(beg:end))
    allocate(pcf%cpool_to_livecrootc_storage(beg:end))
    allocate(pcf%cpool_to_deadcrootc(beg:end))
    allocate(pcf%cpool_to_deadcrootc_storage(beg:end))
    allocate(pcf%cpool_to_gresp_storage(beg:end))
    allocate(pcf%cpool_leaf_gr(beg:end))
    allocate(pcf%cpool_leaf_storage_gr(beg:end))
    allocate(pcf%transfer_leaf_gr(beg:end))
    allocate(pcf%cpool_froot_gr(beg:end))
    allocate(pcf%cpool_froot_storage_gr(beg:end))
    allocate(pcf%transfer_froot_gr(beg:end))
    allocate(pcf%cpool_livestem_gr(beg:end))
    allocate(pcf%cpool_livestem_storage_gr(beg:end))
    allocate(pcf%transfer_livestem_gr(beg:end))
    allocate(pcf%cpool_deadstem_gr(beg:end))
    allocate(pcf%cpool_deadstem_storage_gr(beg:end))
    allocate(pcf%transfer_deadstem_gr(beg:end))
    allocate(pcf%cpool_livecroot_gr(beg:end))
    allocate(pcf%cpool_livecroot_storage_gr(beg:end))
    allocate(pcf%transfer_livecroot_gr(beg:end))
    allocate(pcf%cpool_deadcroot_gr(beg:end))
    allocate(pcf%cpool_deadcroot_storage_gr(beg:end))
    allocate(pcf%transfer_deadcroot_gr(beg:end))
    allocate(pcf%leafc_storage_to_xfer(beg:end))
    allocate(pcf%frootc_storage_to_xfer(beg:end))
    allocate(pcf%livestemc_storage_to_xfer(beg:end))
    allocate(pcf%deadstemc_storage_to_xfer(beg:end))
    allocate(pcf%livecrootc_storage_to_xfer(beg:end))
    allocate(pcf%deadcrootc_storage_to_xfer(beg:end))
    allocate(pcf%gresp_storage_to_xfer(beg:end))
    allocate(pcf%livestemc_to_deadstemc(beg:end))
    allocate(pcf%livecrootc_to_deadcrootc(beg:end))
    allocate(pcf%gpp(beg:end))
    allocate(pcf%mr(beg:end))
    allocate(pcf%current_gr(beg:end))
    allocate(pcf%transfer_gr(beg:end))
    allocate(pcf%storage_gr(beg:end))
    allocate(pcf%gr(beg:end))
    allocate(pcf%ar(beg:end))
    allocate(pcf%rr(beg:end))
    allocate(pcf%npp(beg:end))
    allocate(pcf%agnpp(beg:end))
    allocate(pcf%bgnpp(beg:end))
    allocate(pcf%litfall(beg:end))
    allocate(pcf%vegfire(beg:end))
    allocate(pcf%pft_cinputs(beg:end))
    allocate(pcf%pft_coutputs(beg:end))
    allocate(pcf%pft_fire_closs(beg:end))

    pcf%psnsun(beg:end) = nan
    pcf%psnsha(beg:end) = nan
    pcf%fpsn(beg:end) = nan
    pcf%frm(beg:end) = nan
    pcf%frmf(beg:end) = nan
    pcf%frms(beg:end) = nan
    pcf%frmr(beg:end) = nan
    pcf%frg(beg:end) = nan
    pcf%dmi(beg:end) = nan
    pcf%fco2(beg:end) = 0._r8
    pcf%fmicr(beg:end) = nan
    pcf%m_leafc_to_litter(beg:end) = nan
    pcf%m_frootc_to_litter(beg:end) = nan
    pcf%m_leafc_storage_to_litter(beg:end) = nan
    pcf%m_frootc_storage_to_litter(beg:end) = nan
    pcf%m_livestemc_storage_to_litter(beg:end) = nan
    pcf%m_deadstemc_storage_to_litter(beg:end) = nan
    pcf%m_livecrootc_storage_to_litter(beg:end) = nan
    pcf%m_deadcrootc_storage_to_litter(beg:end) = nan
    pcf%m_leafc_xfer_to_litter(beg:end) = nan
    pcf%m_frootc_xfer_to_litter(beg:end) = nan
    pcf%m_livestemc_xfer_to_litter(beg:end) = nan
    pcf%m_deadstemc_xfer_to_litter(beg:end) = nan
    pcf%m_livecrootc_xfer_to_litter(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_litter(beg:end) = nan
    pcf%m_livestemc_to_litter(beg:end) = nan
    pcf%m_deadstemc_to_litter(beg:end) = nan
    pcf%m_livecrootc_to_litter(beg:end) = nan
    pcf%m_deadcrootc_to_litter(beg:end) = nan
    pcf%m_gresp_storage_to_litter(beg:end) = nan
    pcf%m_gresp_xfer_to_litter(beg:end) = nan
    pcf%m_leafc_to_fire(beg:end) = nan
    pcf%m_frootc_to_fire(beg:end) = nan
    pcf%m_leafc_storage_to_fire(beg:end) = nan
    pcf%m_frootc_storage_to_fire(beg:end) = nan
    pcf%m_livestemc_storage_to_fire(beg:end) = nan
    pcf%m_deadstemc_storage_to_fire(beg:end) = nan
    pcf%m_livecrootc_storage_to_fire(beg:end) = nan
    pcf%m_deadcrootc_storage_to_fire(beg:end) = nan
    pcf%m_leafc_xfer_to_fire(beg:end) = nan
    pcf%m_frootc_xfer_to_fire(beg:end) = nan
    pcf%m_livestemc_xfer_to_fire(beg:end) = nan
    pcf%m_deadstemc_xfer_to_fire(beg:end) = nan
    pcf%m_livecrootc_xfer_to_fire(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_fire(beg:end) = nan
    pcf%m_livestemc_to_fire(beg:end) = nan
    pcf%m_deadstemc_to_fire(beg:end) = nan
    pcf%m_deadstemc_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_to_fire(beg:end) = nan
    pcf%m_deadcrootc_to_fire(beg:end) = nan
    pcf%m_deadcrootc_to_litter_fire(beg:end) = nan
    pcf%m_gresp_storage_to_fire(beg:end) = nan
    pcf%m_gresp_xfer_to_fire(beg:end) = nan
    pcf%leafc_xfer_to_leafc(beg:end) = nan
    pcf%frootc_xfer_to_frootc(beg:end) = nan
    pcf%livestemc_xfer_to_livestemc(beg:end) = nan
    pcf%deadstemc_xfer_to_deadstemc(beg:end) = nan
    pcf%livecrootc_xfer_to_livecrootc(beg:end) = nan
    pcf%deadcrootc_xfer_to_deadcrootc(beg:end) = nan
    pcf%leafc_to_litter(beg:end) = nan
    pcf%frootc_to_litter(beg:end) = nan
    pcf%leaf_mr(beg:end) = nan
    pcf%froot_mr(beg:end) = nan
    pcf%livestem_mr(beg:end) = nan
    pcf%livecroot_mr(beg:end) = nan
    pcf%leaf_curmr(beg:end) = nan
    pcf%froot_curmr(beg:end) = nan
    pcf%livestem_curmr(beg:end) = nan
    pcf%livecroot_curmr(beg:end) = nan
    pcf%leaf_xsmr(beg:end) = nan
    pcf%froot_xsmr(beg:end) = nan
    pcf%livestem_xsmr(beg:end) = nan
    pcf%livecroot_xsmr(beg:end) = nan
    pcf%psnsun_to_cpool(beg:end) = nan
    pcf%psnshade_to_cpool(beg:end) = nan
    pcf%cpool_to_xsmrpool(beg:end) = nan
    pcf%cpool_to_leafc(beg:end) = nan
    pcf%cpool_to_leafc_storage(beg:end) = nan
    pcf%cpool_to_frootc(beg:end) = nan
    pcf%cpool_to_frootc_storage(beg:end) = nan
    pcf%cpool_to_livestemc(beg:end) = nan
    pcf%cpool_to_livestemc_storage(beg:end) = nan
    pcf%cpool_to_deadstemc(beg:end) = nan
    pcf%cpool_to_deadstemc_storage(beg:end) = nan
    pcf%cpool_to_livecrootc(beg:end) = nan
    pcf%cpool_to_livecrootc_storage(beg:end) = nan
    pcf%cpool_to_deadcrootc(beg:end) = nan
    pcf%cpool_to_deadcrootc_storage(beg:end) = nan
    pcf%cpool_to_gresp_storage(beg:end) = nan
    pcf%cpool_leaf_gr(beg:end) = nan
    pcf%cpool_leaf_storage_gr(beg:end) = nan
    pcf%transfer_leaf_gr(beg:end) = nan
    pcf%cpool_froot_gr(beg:end) = nan
    pcf%cpool_froot_storage_gr(beg:end) = nan
    pcf%transfer_froot_gr(beg:end) = nan
    pcf%cpool_livestem_gr(beg:end) = nan
    pcf%cpool_livestem_storage_gr(beg:end) = nan
    pcf%transfer_livestem_gr(beg:end) = nan
    pcf%cpool_deadstem_gr(beg:end) = nan
    pcf%cpool_deadstem_storage_gr(beg:end) = nan
    pcf%transfer_deadstem_gr(beg:end) = nan
    pcf%cpool_livecroot_gr(beg:end) = nan
    pcf%cpool_livecroot_storage_gr(beg:end) = nan
    pcf%transfer_livecroot_gr(beg:end) = nan
    pcf%cpool_deadcroot_gr(beg:end) = nan
    pcf%cpool_deadcroot_storage_gr(beg:end) = nan
    pcf%transfer_deadcroot_gr(beg:end) = nan
    pcf%leafc_storage_to_xfer(beg:end) = nan
    pcf%frootc_storage_to_xfer(beg:end) = nan
    pcf%livestemc_storage_to_xfer(beg:end) = nan
    pcf%deadstemc_storage_to_xfer(beg:end) = nan
    pcf%livecrootc_storage_to_xfer(beg:end) = nan
    pcf%deadcrootc_storage_to_xfer(beg:end) = nan
    pcf%gresp_storage_to_xfer(beg:end) = nan
    pcf%livestemc_to_deadstemc(beg:end) = nan
    pcf%livecrootc_to_deadcrootc(beg:end) = nan
    pcf%gpp(beg:end) = nan
    pcf%mr(beg:end) = nan
    pcf%current_gr(beg:end) = nan
    pcf%transfer_gr(beg:end) = nan
    pcf%storage_gr(beg:end) = nan
    pcf%gr(beg:end) = nan
    pcf%ar(beg:end) = nan
    pcf%rr(beg:end) = nan
    pcf%npp(beg:end) = nan
    pcf%agnpp(beg:end) = nan
    pcf%bgnpp(beg:end) = nan
    pcf%litfall(beg:end) = nan
    pcf%vegfire(beg:end) = nan
    pcf%pft_cinputs(beg:end) = nan
    pcf%pft_coutputs(beg:end) = nan
    pcf%pft_fire_closs(beg:end) = nan

  end subroutine init_pft_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nflux_type
!
! !INTERFACE:
  subroutine init_pft_nflux_type(beg, end, pnf)
!
! !DESCRIPTION:
! Initialize pft nitrogen flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nflux_type), intent(inout) :: pnf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pnf%m_leafn_to_litter(beg:end))
    allocate(pnf%m_frootn_to_litter(beg:end))
    allocate(pnf%m_leafn_storage_to_litter(beg:end))
    allocate(pnf%m_frootn_storage_to_litter(beg:end))
    allocate(pnf%m_livestemn_storage_to_litter(beg:end))
    allocate(pnf%m_deadstemn_storage_to_litter(beg:end))
    allocate(pnf%m_livecrootn_storage_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_litter(beg:end))
    allocate(pnf%m_leafn_xfer_to_litter(beg:end))
    allocate(pnf%m_frootn_xfer_to_litter(beg:end))
    allocate(pnf%m_livestemn_xfer_to_litter(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_litter(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_litter(beg:end))
    allocate(pnf%m_livestemn_to_litter(beg:end))
    allocate(pnf%m_deadstemn_to_litter(beg:end))
    allocate(pnf%m_livecrootn_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_to_litter(beg:end))
    allocate(pnf%m_retransn_to_litter(beg:end))
    allocate(pnf%m_leafn_to_fire(beg:end))
    allocate(pnf%m_frootn_to_fire(beg:end))
    allocate(pnf%m_leafn_storage_to_fire(beg:end))
    allocate(pnf%m_frootn_storage_to_fire(beg:end))
    allocate(pnf%m_livestemn_storage_to_fire(beg:end))
    allocate(pnf%m_deadstemn_storage_to_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_fire(beg:end))
    allocate(pnf%m_leafn_xfer_to_fire(beg:end))
    allocate(pnf%m_frootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_fire(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_litter_fire(beg:end))
    allocate(pnf%m_retransn_to_fire(beg:end))
    allocate(pnf%leafn_xfer_to_leafn(beg:end))
    allocate(pnf%frootn_xfer_to_frootn(beg:end))
    allocate(pnf%livestemn_xfer_to_livestemn(beg:end))
    allocate(pnf%deadstemn_xfer_to_deadstemn(beg:end))
    allocate(pnf%livecrootn_xfer_to_livecrootn(beg:end))
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(beg:end))
    allocate(pnf%leafn_to_litter(beg:end))
    allocate(pnf%leafn_to_retransn(beg:end))
    allocate(pnf%frootn_to_litter(beg:end))
    allocate(pnf%retransn_to_npool(beg:end))
    allocate(pnf%sminn_to_npool(beg:end))
    allocate(pnf%npool_to_leafn(beg:end))
    allocate(pnf%npool_to_leafn_storage(beg:end))
    allocate(pnf%npool_to_frootn(beg:end))
    allocate(pnf%npool_to_frootn_storage(beg:end))
    allocate(pnf%npool_to_livestemn(beg:end))
    allocate(pnf%npool_to_livestemn_storage(beg:end))
    allocate(pnf%npool_to_deadstemn(beg:end))
    allocate(pnf%npool_to_deadstemn_storage(beg:end))
    allocate(pnf%npool_to_livecrootn(beg:end))
    allocate(pnf%npool_to_livecrootn_storage(beg:end))
    allocate(pnf%npool_to_deadcrootn(beg:end))
    allocate(pnf%npool_to_deadcrootn_storage(beg:end))
    allocate(pnf%leafn_storage_to_xfer(beg:end))
    allocate(pnf%frootn_storage_to_xfer(beg:end))
    allocate(pnf%livestemn_storage_to_xfer(beg:end))
    allocate(pnf%deadstemn_storage_to_xfer(beg:end))
    allocate(pnf%livecrootn_storage_to_xfer(beg:end))
    allocate(pnf%deadcrootn_storage_to_xfer(beg:end))
    allocate(pnf%livestemn_to_deadstemn(beg:end))
    allocate(pnf%livestemn_to_retransn(beg:end))
    allocate(pnf%livecrootn_to_deadcrootn(beg:end))
    allocate(pnf%livecrootn_to_retransn(beg:end))
    allocate(pnf%ndeploy(beg:end))
    allocate(pnf%pft_ninputs(beg:end))
    allocate(pnf%pft_noutputs(beg:end))
    allocate(pnf%pft_fire_nloss(beg:end))

    pnf%m_leafn_to_litter(beg:end) = nan
    pnf%m_frootn_to_litter(beg:end) = nan
    pnf%m_leafn_storage_to_litter(beg:end) = nan
    pnf%m_frootn_storage_to_litter(beg:end) = nan
    pnf%m_livestemn_storage_to_litter(beg:end) = nan
    pnf%m_deadstemn_storage_to_litter(beg:end) = nan
    pnf%m_livecrootn_storage_to_litter(beg:end) = nan
    pnf%m_deadcrootn_storage_to_litter(beg:end) = nan
    pnf%m_leafn_xfer_to_litter(beg:end) = nan
    pnf%m_frootn_xfer_to_litter(beg:end) = nan
    pnf%m_livestemn_xfer_to_litter(beg:end) = nan
    pnf%m_deadstemn_xfer_to_litter(beg:end) = nan
    pnf%m_livecrootn_xfer_to_litter(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_litter(beg:end) = nan
    pnf%m_livestemn_to_litter(beg:end) = nan
    pnf%m_deadstemn_to_litter(beg:end) = nan
    pnf%m_livecrootn_to_litter(beg:end) = nan
    pnf%m_deadcrootn_to_litter(beg:end) = nan
    pnf%m_retransn_to_litter(beg:end) = nan
    pnf%m_leafn_to_fire(beg:end) = nan
    pnf%m_frootn_to_fire(beg:end) = nan
    pnf%m_leafn_storage_to_fire(beg:end) = nan
    pnf%m_frootn_storage_to_fire(beg:end) = nan
    pnf%m_livestemn_storage_to_fire(beg:end) = nan
    pnf%m_deadstemn_storage_to_fire(beg:end) = nan
    pnf%m_livecrootn_storage_to_fire(beg:end) = nan
    pnf%m_deadcrootn_storage_to_fire(beg:end) = nan
    pnf%m_leafn_xfer_to_fire(beg:end) = nan
    pnf%m_frootn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_xfer_to_fire(beg:end) = nan
    pnf%m_deadstemn_xfer_to_fire(beg:end) = nan
    pnf%m_livecrootn_xfer_to_fire(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_to_fire(beg:end) = nan
    pnf%m_deadstemn_to_fire(beg:end) = nan
    pnf%m_deadstemn_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_to_fire(beg:end) = nan
    pnf%m_deadcrootn_to_fire(beg:end) = nan
    pnf%m_deadcrootn_to_litter_fire(beg:end) = nan
    pnf%m_retransn_to_fire(beg:end) = nan
    pnf%leafn_xfer_to_leafn(beg:end) = nan
    pnf%frootn_xfer_to_frootn(beg:end) = nan
    pnf%livestemn_xfer_to_livestemn(beg:end) = nan
    pnf%deadstemn_xfer_to_deadstemn(beg:end) = nan
    pnf%livecrootn_xfer_to_livecrootn(beg:end) = nan
    pnf%deadcrootn_xfer_to_deadcrootn(beg:end) = nan
    pnf%leafn_to_litter(beg:end) = nan
    pnf%leafn_to_retransn(beg:end) = nan
    pnf%frootn_to_litter(beg:end) = nan
    pnf%retransn_to_npool(beg:end) = nan
    pnf%sminn_to_npool(beg:end) = nan
    pnf%npool_to_leafn(beg:end) = nan
    pnf%npool_to_leafn_storage(beg:end) = nan
    pnf%npool_to_frootn(beg:end) = nan
    pnf%npool_to_frootn_storage(beg:end) = nan
    pnf%npool_to_livestemn(beg:end) = nan
    pnf%npool_to_livestemn_storage(beg:end) = nan
    pnf%npool_to_deadstemn(beg:end) = nan
    pnf%npool_to_deadstemn_storage(beg:end) = nan
    pnf%npool_to_livecrootn(beg:end) = nan
    pnf%npool_to_livecrootn_storage(beg:end) = nan
    pnf%npool_to_deadcrootn(beg:end) = nan
    pnf%npool_to_deadcrootn_storage(beg:end) = nan
    pnf%leafn_storage_to_xfer(beg:end) = nan
    pnf%frootn_storage_to_xfer(beg:end) = nan
    pnf%livestemn_storage_to_xfer(beg:end) = nan
    pnf%deadstemn_storage_to_xfer(beg:end) = nan
    pnf%livecrootn_storage_to_xfer(beg:end) = nan
    pnf%deadcrootn_storage_to_xfer(beg:end) = nan
    pnf%livestemn_to_deadstemn(beg:end) = nan
    pnf%livestemn_to_retransn(beg:end) = nan
    pnf%livecrootn_to_deadcrootn(beg:end) = nan
    pnf%livecrootn_to_retransn(beg:end) = nan
    pnf%ndeploy(beg:end) = nan
    pnf%pft_ninputs(beg:end) = nan
    pnf%pft_noutputs(beg:end) = nan
    pnf%pft_fire_nloss(beg:end) = nan

  end subroutine init_pft_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vflux_type
!
! !INTERFACE:
  subroutine init_pft_vflux_type(beg, end, pvf)
!
! !DESCRIPTION:
! Initialize pft VOC flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vflux_type), intent(inout) :: pvf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pvf%vocflx_tot(beg:end))
    allocate(pvf%vocflx(beg:end,1:nvoc))
    allocate(pvf%vocflx_1(beg:end))
    allocate(pvf%vocflx_2(beg:end))
    allocate(pvf%vocflx_3(beg:end))
    allocate(pvf%vocflx_4(beg:end))
    allocate(pvf%vocflx_5(beg:end))

    pvf%vocflx_tot(beg:end) = nan
    pvf%vocflx(beg:end,1:nvoc) = nan
    pvf%vocflx_1(beg:end) = nan
    pvf%vocflx_2(beg:end) = nan
    pvf%vocflx_3(beg:end) = nan
    pvf%vocflx_4(beg:end) = nan
    pvf%vocflx_5(beg:end) = nan

  end subroutine init_pft_vflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_dflux_type
!
! !INTERFACE:
  subroutine init_pft_dflux_type(beg, end, pdf)
!
! !DESCRIPTION:
! Initialize pft dust flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dflux_type), intent(inout):: pdf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdf%flx_mss_vrt_dst(beg:end,1:ndst))
    allocate(pdf%flx_mss_vrt_dst_tot(beg:end))
    allocate(pdf%vlc_trb(beg:end,1:ndst))
    allocate(pdf%vlc_trb_1(beg:end))
    allocate(pdf%vlc_trb_2(beg:end))
    allocate(pdf%vlc_trb_3(beg:end))
    allocate(pdf%vlc_trb_4(beg:end))

    pdf%flx_mss_vrt_dst(beg:end,1:ndst) = nan
    pdf%flx_mss_vrt_dst_tot(beg:end) = nan
    pdf%vlc_trb(beg:end,1:ndst) = nan
    pdf%vlc_trb_1(beg:end) = nan
    pdf%vlc_trb_2(beg:end) = nan
    pdf%vlc_trb_3(beg:end) = nan
    pdf%vlc_trb_4(beg:end) = nan

  end subroutine init_pft_dflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_pstate_type
!
! !INTERFACE:
  subroutine init_column_pstate_type(beg, end, cps)
!
! !DESCRIPTION:
! Initialize column physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_pstate_type), intent(inout):: cps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cps%snl(beg:end))      !* cannot be averaged up
    allocate(cps%isoicol(beg:end))  !* cannot be averaged up
    allocate(cps%bsw(beg:end,nlevsoi))
    allocate(cps%watsat(beg:end,nlevsoi))
    allocate(cps%watdry(beg:end,nlevsoi)) 
    allocate(cps%watopt(beg:end,nlevsoi)) 
    allocate(cps%hksat(beg:end,nlevsoi))
    allocate(cps%sucsat(beg:end,nlevsoi))
    allocate(cps%csol(beg:end,nlevsoi))
    allocate(cps%tkmg(beg:end,nlevsoi))
    allocate(cps%tkdry(beg:end,nlevsoi))
    allocate(cps%tksatu(beg:end,nlevsoi))
    allocate(cps%smpmin(beg:end))
    allocate(cps%hkdepth(beg:end))
    allocate(cps%wtfact(beg:end))
    allocate(cps%fracice(beg:end,nlevsoi))
    allocate(cps%gwc_thr(beg:end))
    allocate(cps%mss_frc_cly_vld(beg:end))
    allocate(cps%mbl_bsn_fct(beg:end))
    allocate(cps%do_capsnow(beg:end))
    allocate(cps%snowdp(beg:end))
    allocate(cps%snowage(beg:end))
    allocate(cps%frac_sno (beg:end))
    allocate(cps%zi(beg:end,-nlevsno+0:nlevsoi))
    allocate(cps%dz(beg:end,-nlevsno+1:nlevsoi))
    allocate(cps%z (beg:end,-nlevsno+1:nlevsoi))
    allocate(cps%frac_iceold(beg:end,-nlevsno+1:nlevsoi))
    allocate(cps%imelt(beg:end,-nlevsno+1:nlevsoi))
    allocate(cps%eff_porosity(beg:end,nlevsoi))
    allocate(cps%emg(beg:end))
    allocate(cps%z0mg(beg:end))
    allocate(cps%z0hg(beg:end))
    allocate(cps%z0qg(beg:end))
    allocate(cps%htvp(beg:end))
    allocate(cps%beta(beg:end))
    allocate(cps%zii(beg:end))
    allocate(cps%albgrd(beg:end,numrad))
    allocate(cps%albgri(beg:end,numrad))
    allocate(cps%rootr_column(beg:end,nlevsoi))
    allocate(cps%wf(beg:end))
    allocate(cps%bsw2(beg:end,nlevsoi))
    allocate(cps%psisat(beg:end,nlevsoi))
    allocate(cps%vwcsat(beg:end,nlevsoi))
    allocate(cps%soilpsi(beg:end,nlevsoi))
    allocate(cps%decl(beg:end))
    allocate(cps%coszen(beg:end))
    allocate(cps%fpi(beg:end))
    allocate(cps%fpg(beg:end))
    allocate(cps%annsum_counter(beg:end))
    allocate(cps%cannsum_npp(beg:end))
    allocate(cps%me(beg:end))
    allocate(cps%fire_prob(beg:end))
    allocate(cps%mean_fire_prob(beg:end))
    allocate(cps%fireseasonl(beg:end))
    allocate(cps%farea_burned(beg:end))
    allocate(cps%ann_farea_burned(beg:end))

    cps%isoicol(beg:end) = bigint
    cps%bsw(beg:end,1:nlevsoi) = nan
    cps%watsat(beg:end,1:nlevsoi) = nan
    cps%watdry(beg:end,1:nlevsoi) = nan  
    cps%watopt(beg:end,1:nlevsoi) = nan  
    cps%hksat(beg:end,1:nlevsoi) = nan
    cps%sucsat(beg:end,1:nlevsoi) = nan
    cps%csol(beg:end,1:nlevsoi) = nan
    cps%tkmg(beg:end,1:nlevsoi) = nan
    cps%tkdry(beg:end,1:nlevsoi) = nan
    cps%tksatu(beg:end,1:nlevsoi) = nan
    cps%smpmin(beg:end) = nan
    cps%hkdepth(beg:end) = nan
    cps%wtfact(beg:end) = nan
    cps%fracice(beg:end,1:nlevsoi) = nan
    cps%gwc_thr(beg:end) = nan
    cps%mss_frc_cly_vld(beg:end) = nan
    cps%mbl_bsn_fct(beg:end) = nan
    cps%do_capsnow (beg:end)= .false.
    cps%snowdp(beg:end) = nan
    cps%snowage(beg:end) = nan
    cps%frac_sno(beg:end) = nan
    cps%zi(beg:end,-nlevsno+0:nlevsoi) = nan
    cps%dz(beg:end,-nlevsno+1:nlevsoi) = nan
    cps%z (beg:end,-nlevsno+1:nlevsoi) = nan
    cps%frac_iceold(beg:end,-nlevsno+1:nlevsoi) = nan
    cps%imelt(beg:end,-nlevsno+1:nlevsoi) = bigint
    cps%eff_porosity(beg:end,1:nlevsoi) = nan
    cps%emg(beg:end) = nan
    cps%z0mg(beg:end) = nan
    cps%z0hg(beg:end) = nan
    cps%z0qg(beg:end) = nan
    cps%htvp(beg:end) = nan
    cps%beta(beg:end) = nan
    cps%zii(beg:end) = nan
    cps%albgrd(beg:end,:numrad) = nan
    cps%albgri(beg:end,:numrad) = nan
    cps%rootr_column(beg:end,1:nlevsoi) = nan
    cps%wf(beg:end) = nan
    cps%bsw2(beg:end,1:nlevsoi) = nan
    cps%psisat(beg:end,1:nlevsoi) = nan
    cps%vwcsat(beg:end,1:nlevsoi) = nan
    cps%soilpsi(beg:end,1:nlevsoi) = nan
    cps%decl(beg:end) = nan
    cps%coszen(beg:end) = nan
    cps%fpi(beg:end) = nan
    cps%fpg(beg:end) = nan
    cps%annsum_counter(beg:end) = nan
    cps%cannsum_npp(beg:end) = nan
    cps%me(beg:end) = nan
    cps%fire_prob(beg:end) = nan
    cps%mean_fire_prob(beg:end) = nan
    cps%fireseasonl(beg:end) = nan
    cps%farea_burned(beg:end) = nan
    cps%ann_farea_burned(beg:end) = nan

  end subroutine init_column_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_estate_type
!
! !INTERFACE:
  subroutine init_column_estate_type(beg, end, ces)
!
! !DESCRIPTION:
! Initialize column energy state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_estate_type), intent(inout):: ces
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(ces%t_grnd(beg:end))
    allocate(ces%dt_grnd(beg:end))
    allocate(ces%t_soisno(beg:end,-nlevsno+1:nlevsoi))
    allocate(ces%t_lake(beg:end,1:nlevlak))
    allocate(ces%tssbef(beg:end,-nlevsno+1:nlevsoi))
    allocate(ces%t_snow(beg:end))
    allocate(ces%thv(beg:end))
    allocate(ces%thm(beg:end))

    ces%t_grnd(beg:end) = nan
    ces%dt_grnd(beg:end) = nan
    ces%t_soisno(beg:end,-nlevsno+1:nlevsoi) = nan
    ces%t_lake(beg:end,1:nlevlak)= nan
    ces%tssbef(beg:end,-nlevsno+1:nlevsoi) = nan
    ces%t_snow(beg:end) = nan
    ces%thv(beg:end) = nan
    ces%thm(beg:end) = nan

  end subroutine init_column_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wstate_type
!
! !INTERFACE:
  subroutine init_column_wstate_type(beg, end, cws)
!
! !DESCRIPTION:
! Initialize column water state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wstate_type), intent(inout):: cws !column water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cws%h2osno(beg:end))
    allocate(cws%h2osoi_liq(beg:end,-nlevsno+1:nlevsoi))
    allocate(cws%h2osoi_ice(beg:end,-nlevsno+1:nlevsoi))
    allocate(cws%h2osoi_vol(beg:end,1:nlevsoi))
    allocate(cws%h2osno_old(beg:end))
    allocate(cws%qg(beg:end))
    allocate(cws%dqgdT(beg:end))
    allocate(cws%snowice(beg:end))
    allocate(cws%snowliq(beg:end))
    allocate(cws%soilalpha(beg:end))
    allocate(cws%zwt(beg:end))
    allocate(cws%fcov(beg:end))
    allocate(cws%wa(beg:end))
    allocate(cws%wt(beg:end))
    allocate(cws%qcharge(beg:end))

    cws%h2osno(beg:end) = nan
    cws%h2osoi_liq(beg:end,-nlevsno+1:nlevsoi)= nan
    cws%h2osoi_ice(beg:end,-nlevsno+1:nlevsoi) = nan
    cws%h2osoi_vol(beg:end,1:nlevsoi) = nan
    cws%h2osno_old(beg:end) = nan
    cws%qg(beg:end) = nan
    cws%dqgdT(beg:end) = nan
    cws%snowice(beg:end) = nan
    cws%snowliq(beg:end) = nan
    cws%soilalpha(beg:end) = nan
    cws%zwt(beg:end) = nan
    cws%fcov(beg:end) = nan
    cws%wa(beg:end) = nan
    cws%wt(beg:end) = nan
    cws%qcharge(beg:end) = nan

  end subroutine init_column_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cstate_type
!
! !INTERFACE:
  subroutine init_column_cstate_type(beg, end, ccs)
!
! !DESCRIPTION:
! Initialize column carbon state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cstate_type), intent(inout):: ccs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ccs%soilc(beg:end))
    allocate(ccs%cwdc(beg:end))
    allocate(ccs%litr1c(beg:end))
    allocate(ccs%litr2c(beg:end))
    allocate(ccs%litr3c(beg:end))
    allocate(ccs%soil1c(beg:end))
    allocate(ccs%soil2c(beg:end))
    allocate(ccs%soil3c(beg:end))
    allocate(ccs%soil4c(beg:end))
    allocate(ccs%col_ctrunc(beg:end))
    allocate(ccs%totlitc(beg:end))
    allocate(ccs%totsomc(beg:end))
    allocate(ccs%totecosysc(beg:end))
    allocate(ccs%totcolc(beg:end))

    ccs%soilc(beg:end) = nan
    ccs%cwdc(beg:end) = nan
    ccs%litr1c(beg:end) = nan
    ccs%litr2c(beg:end) = nan
    ccs%litr3c(beg:end) = nan
    ccs%soil1c(beg:end) = nan
    ccs%soil2c(beg:end) = nan
    ccs%soil3c(beg:end) = nan
    ccs%soil4c(beg:end) = nan
    ccs%col_ctrunc(beg:end) = nan
    ccs%totlitc(beg:end) = nan
    ccs%totsomc(beg:end) = nan
    ccs%totecosysc(beg:end) = nan
    ccs%totcolc(beg:end) = nan

  end subroutine init_column_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nstate_type
!
! !INTERFACE:
  subroutine init_column_nstate_type(beg, end, cns)
!
! !DESCRIPTION:
! Initialize column nitrogen state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nstate_type), intent(inout):: cns
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(cns%cwdn(beg:end))
    allocate(cns%litr1n(beg:end))
    allocate(cns%litr2n(beg:end))
    allocate(cns%litr3n(beg:end))
    allocate(cns%soil1n(beg:end))
    allocate(cns%soil2n(beg:end))
    allocate(cns%soil3n(beg:end))
    allocate(cns%soil4n(beg:end))
    allocate(cns%sminn(beg:end))
    allocate(cns%col_ntrunc(beg:end))
    allocate(cns%totlitn(beg:end))
    allocate(cns%totsomn(beg:end))
    allocate(cns%totecosysn(beg:end))
    allocate(cns%totcoln(beg:end))

    cns%cwdn(beg:end) = nan
    cns%litr1n(beg:end) = nan
    cns%litr2n(beg:end) = nan
    cns%litr3n(beg:end) = nan
    cns%soil1n(beg:end) = nan
    cns%soil2n(beg:end) = nan
    cns%soil3n(beg:end) = nan
    cns%soil4n(beg:end) = nan
    cns%sminn(beg:end) = nan
    cns%col_ntrunc(beg:end) = nan
    cns%totlitn(beg:end) = nan
    cns%totsomn(beg:end) = nan
    cns%totecosysn(beg:end) = nan
    cns%totcoln(beg:end) = nan

  end subroutine init_column_nstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_eflux_type
!
! !INTERFACE:
  subroutine init_column_eflux_type(beg, end, cef)
!
! !DESCRIPTION:
! Initialize column energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_eflux_type), intent(inout):: cef
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cef%eflx_snomelt(beg:end))
    allocate(cef%eflx_impsoil(beg:end))

    cef%eflx_snomelt(beg:end) = nan
    cef%eflx_impsoil(beg:end) = nan

  end subroutine init_column_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wflux_type
!
! !INTERFACE:
  subroutine init_column_wflux_type(beg, end, cwf)
!
! !DESCRIPTION:
! Initialize column water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wflux_type), intent(inout):: cwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cwf%qflx_infl(beg:end))
    allocate(cwf%qflx_surf(beg:end))
    allocate(cwf%qflx_drain(beg:end))
    allocate(cwf%qflx_top_soil(beg:end))
    allocate(cwf%qflx_snomelt(beg:end))
    allocate(cwf%qflx_qrgwl(beg:end))
    allocate(cwf%qmelt(beg:end))
    allocate(cwf%h2ocan_loss(beg:end))

    cwf%qflx_infl(beg:end) = nan
    cwf%qflx_surf(beg:end) = nan
    cwf%qflx_drain(beg:end) = nan
    cwf%qflx_top_soil(beg:end) = nan
    cwf%qflx_snomelt(beg:end) = nan
    cwf%qflx_qrgwl(beg:end) = nan
    cwf%qmelt(beg:end) = nan
    cwf%h2ocan_loss(beg:end) = nan

  end subroutine init_column_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cflux_type
!
! !INTERFACE:
  subroutine init_column_cflux_type(beg, end, ccf)
!
! !DESCRIPTION:
! Initialize column carbon flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cflux_type), intent(inout):: ccf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(ccf%m_leafc_to_litr1c(beg:end))
    allocate(ccf%m_leafc_to_litr2c(beg:end))
    allocate(ccf%m_leafc_to_litr3c(beg:end))
    allocate(ccf%m_frootc_to_litr1c(beg:end))
    allocate(ccf%m_frootc_to_litr2c(beg:end))
    allocate(ccf%m_frootc_to_litr3c(beg:end))
    allocate(ccf%m_leafc_storage_to_litr1c(beg:end))
    allocate(ccf%m_frootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_storage_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_storage_to_litr1c(beg:end))
    allocate(ccf%m_livecrootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_deadcrootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_leafc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_frootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livecrootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadcrootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_to_cwdc(beg:end))
    allocate(ccf%m_deadstemc_to_cwdc(beg:end))
    allocate(ccf%m_livecrootc_to_cwdc(beg:end))
    allocate(ccf%m_deadcrootc_to_cwdc(beg:end))
    allocate(ccf%m_gresp_storage_to_litr1c(beg:end))
    allocate(ccf%m_gresp_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_to_cwdc_fire(beg:end))
    allocate(ccf%m_deadcrootc_to_cwdc_fire(beg:end))
    allocate(ccf%m_litr1c_to_fire(beg:end))
    allocate(ccf%m_litr2c_to_fire(beg:end))
    allocate(ccf%m_litr3c_to_fire(beg:end))
    allocate(ccf%m_cwdc_to_fire(beg:end))
    allocate(ccf%leafc_to_litr1c(beg:end))
    allocate(ccf%leafc_to_litr2c(beg:end))
    allocate(ccf%leafc_to_litr3c(beg:end))
    allocate(ccf%frootc_to_litr1c(beg:end))
    allocate(ccf%frootc_to_litr2c(beg:end))
    allocate(ccf%frootc_to_litr3c(beg:end))
    allocate(ccf%cwdc_to_litr2c(beg:end))
    allocate(ccf%cwdc_to_litr3c(beg:end))
    allocate(ccf%litr1_hr(beg:end))
    allocate(ccf%litr1c_to_soil1c(beg:end))
    allocate(ccf%litr2_hr(beg:end))
    allocate(ccf%litr2c_to_soil2c(beg:end))
    allocate(ccf%litr3_hr(beg:end))
    allocate(ccf%litr3c_to_soil3c(beg:end))
    allocate(ccf%soil1_hr(beg:end))
    allocate(ccf%soil1c_to_soil2c(beg:end))
    allocate(ccf%soil2_hr(beg:end))
    allocate(ccf%soil2c_to_soil3c(beg:end))
    allocate(ccf%soil3_hr(beg:end))
    allocate(ccf%soil3c_to_soil4c(beg:end))
    allocate(ccf%soil4_hr(beg:end))
    allocate(ccf%lithr(beg:end))
    allocate(ccf%somhr(beg:end))
    allocate(ccf%hr(beg:end))
    allocate(ccf%sr(beg:end))
    allocate(ccf%er(beg:end))
    allocate(ccf%litfire(beg:end))
    allocate(ccf%somfire(beg:end))
    allocate(ccf%totfire(beg:end))
    allocate(ccf%nep(beg:end))
    allocate(ccf%nee(beg:end))
    allocate(ccf%col_cinputs(beg:end))
    allocate(ccf%col_coutputs(beg:end))
    allocate(ccf%col_fire_closs(beg:end))

    ccf%m_leafc_to_litr1c(beg:end) = nan
    ccf%m_leafc_to_litr2c(beg:end) = nan
    ccf%m_leafc_to_litr3c(beg:end) = nan
    ccf%m_frootc_to_litr1c(beg:end) = nan
    ccf%m_frootc_to_litr2c(beg:end) = nan
    ccf%m_frootc_to_litr3c(beg:end) = nan
    ccf%m_leafc_storage_to_litr1c(beg:end) = nan
    ccf%m_frootc_storage_to_litr1c(beg:end) = nan
    ccf%m_livestemc_storage_to_litr1c(beg:end) = nan
    ccf%m_deadstemc_storage_to_litr1c(beg:end) = nan
    ccf%m_livecrootc_storage_to_litr1c(beg:end) = nan
    ccf%m_deadcrootc_storage_to_litr1c(beg:end) = nan
    ccf%m_leafc_xfer_to_litr1c(beg:end) = nan
    ccf%m_frootc_xfer_to_litr1c(beg:end) = nan
    ccf%m_livestemc_xfer_to_litr1c(beg:end) = nan
    ccf%m_deadstemc_xfer_to_litr1c(beg:end) = nan
    ccf%m_livecrootc_xfer_to_litr1c(beg:end) = nan
    ccf%m_deadcrootc_xfer_to_litr1c(beg:end) = nan
    ccf%m_livestemc_to_cwdc(beg:end) = nan
    ccf%m_deadstemc_to_cwdc(beg:end) = nan
    ccf%m_livecrootc_to_cwdc(beg:end) = nan
    ccf%m_deadcrootc_to_cwdc(beg:end) = nan
    ccf%m_gresp_storage_to_litr1c(beg:end) = nan
    ccf%m_gresp_xfer_to_litr1c(beg:end) = nan
    ccf%m_deadstemc_to_cwdc_fire(beg:end) = nan
    ccf%m_deadcrootc_to_cwdc_fire(beg:end) = nan
    ccf%m_litr1c_to_fire(beg:end) = nan
    ccf%m_litr2c_to_fire(beg:end) = nan
    ccf%m_litr3c_to_fire(beg:end) = nan
    ccf%m_cwdc_to_fire(beg:end) = nan
    ccf%leafc_to_litr1c(beg:end) = nan
    ccf%leafc_to_litr2c(beg:end) = nan
    ccf%leafc_to_litr3c(beg:end) = nan
    ccf%frootc_to_litr1c(beg:end) = nan
    ccf%frootc_to_litr2c(beg:end) = nan
    ccf%frootc_to_litr3c(beg:end) = nan
    ccf%cwdc_to_litr2c(beg:end) = nan
    ccf%cwdc_to_litr3c(beg:end) = nan
    ccf%litr1_hr(beg:end) = nan
    ccf%litr1c_to_soil1c(beg:end) = nan
    ccf%litr2_hr(beg:end) = nan
    ccf%litr2c_to_soil2c(beg:end) = nan
    ccf%litr3_hr(beg:end) = nan
    ccf%litr3c_to_soil3c(beg:end) = nan
    ccf%soil1_hr(beg:end) = nan
    ccf%soil1c_to_soil2c(beg:end) = nan
    ccf%soil2_hr(beg:end) = nan
    ccf%soil2c_to_soil3c(beg:end) = nan
    ccf%soil3_hr(beg:end) = nan
    ccf%soil3c_to_soil4c(beg:end) = nan
    ccf%soil4_hr(beg:end) = nan
    ccf%lithr(beg:end) = nan
    ccf%somhr(beg:end) = nan
    ccf%hr(beg:end) = nan
    ccf%sr(beg:end) = nan
    ccf%er(beg:end) = nan
    ccf%litfire(beg:end) = nan
    ccf%somfire(beg:end) = nan
    ccf%totfire(beg:end) = nan
    ccf%nep(beg:end) = nan
    ccf%nee(beg:end) = nan
    ccf%col_cinputs(beg:end) = nan
    ccf%col_coutputs(beg:end) = nan
    ccf%col_fire_closs(beg:end) = nan

  end subroutine init_column_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nflux_type
!
! !INTERFACE:
  subroutine init_column_nflux_type(beg, end, cnf)
!
! !DESCRIPTION:
! Initialize column nitrogen flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nflux_type), intent(inout):: cnf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(cnf%ndep_to_sminn(beg:end))
    allocate(cnf%nfix_to_sminn(beg:end))
    allocate(cnf%m_leafn_to_litr1n(beg:end))
    allocate(cnf%m_leafn_to_litr2n(beg:end))
    allocate(cnf%m_leafn_to_litr3n(beg:end))
    allocate(cnf%m_frootn_to_litr1n(beg:end))
    allocate(cnf%m_frootn_to_litr2n(beg:end))
    allocate(cnf%m_frootn_to_litr3n(beg:end))
    allocate(cnf%m_leafn_storage_to_litr1n(beg:end))
    allocate(cnf%m_frootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_storage_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_storage_to_litr1n(beg:end))
    allocate(cnf%m_livecrootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_deadcrootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_leafn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_frootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livecrootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_deadcrootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_to_cwdn(beg:end))
    allocate(cnf%m_deadstemn_to_cwdn(beg:end))
    allocate(cnf%m_livecrootn_to_cwdn(beg:end))
    allocate(cnf%m_deadcrootn_to_cwdn(beg:end))
    allocate(cnf%m_retransn_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_to_cwdn_fire(beg:end))
    allocate(cnf%m_deadcrootn_to_cwdn_fire(beg:end))
    allocate(cnf%m_litr1n_to_fire(beg:end))
    allocate(cnf%m_litr2n_to_fire(beg:end))
    allocate(cnf%m_litr3n_to_fire(beg:end))
    allocate(cnf%m_cwdn_to_fire(beg:end))
    allocate(cnf%leafn_to_litr1n(beg:end))
    allocate(cnf%leafn_to_litr2n(beg:end))
    allocate(cnf%leafn_to_litr3n(beg:end))
    allocate(cnf%frootn_to_litr1n(beg:end))
    allocate(cnf%frootn_to_litr2n(beg:end))
    allocate(cnf%frootn_to_litr3n(beg:end))
    allocate(cnf%cwdn_to_litr2n(beg:end))
    allocate(cnf%cwdn_to_litr3n(beg:end))
    allocate(cnf%litr1n_to_soil1n(beg:end))
    allocate(cnf%sminn_to_soil1n_l1(beg:end))
    allocate(cnf%litr2n_to_soil2n(beg:end))
    allocate(cnf%sminn_to_soil2n_l2(beg:end))
    allocate(cnf%litr3n_to_soil3n(beg:end))
    allocate(cnf%sminn_to_soil3n_l3(beg:end))
    allocate(cnf%soil1n_to_soil2n(beg:end))
    allocate(cnf%sminn_to_soil2n_s1(beg:end))
    allocate(cnf%soil2n_to_soil3n(beg:end))
    allocate(cnf%sminn_to_soil3n_s2(beg:end))
    allocate(cnf%soil3n_to_soil4n(beg:end))
    allocate(cnf%sminn_to_soil4n_s3(beg:end))
    allocate(cnf%soil4n_to_sminn(beg:end))
    allocate(cnf%sminn_to_denit_l1s1(beg:end))
    allocate(cnf%sminn_to_denit_l2s2(beg:end))
    allocate(cnf%sminn_to_denit_l3s3(beg:end))
    allocate(cnf%sminn_to_denit_s1s2(beg:end))
    allocate(cnf%sminn_to_denit_s2s3(beg:end))
    allocate(cnf%sminn_to_denit_s3s4(beg:end))
    allocate(cnf%sminn_to_denit_s4(beg:end))
    allocate(cnf%sminn_to_denit_excess(beg:end))
    allocate(cnf%sminn_leached(beg:end))
    allocate(cnf%potential_immob(beg:end))
    allocate(cnf%actual_immob(beg:end))
    allocate(cnf%sminn_to_plant(beg:end))
    allocate(cnf%supplement_to_sminn(beg:end))
    allocate(cnf%gross_nmin(beg:end))
    allocate(cnf%net_nmin(beg:end))
    allocate(cnf%denit(beg:end))
    allocate(cnf%col_ninputs(beg:end))
    allocate(cnf%col_noutputs(beg:end))
    allocate(cnf%col_fire_nloss(beg:end))

    cnf%ndep_to_sminn(beg:end) = nan
    cnf%nfix_to_sminn(beg:end) = nan
    cnf%m_leafn_to_litr1n(beg:end) = nan
    cnf%m_leafn_to_litr2n(beg:end) = nan
    cnf%m_leafn_to_litr3n(beg:end) = nan
    cnf%m_frootn_to_litr1n(beg:end) = nan
    cnf%m_frootn_to_litr2n(beg:end) = nan
    cnf%m_frootn_to_litr3n(beg:end) = nan
    cnf%m_leafn_storage_to_litr1n(beg:end) = nan
    cnf%m_frootn_storage_to_litr1n(beg:end) = nan
    cnf%m_livestemn_storage_to_litr1n(beg:end) = nan
    cnf%m_deadstemn_storage_to_litr1n(beg:end) = nan
    cnf%m_livecrootn_storage_to_litr1n(beg:end) = nan
    cnf%m_deadcrootn_storage_to_litr1n(beg:end) = nan
    cnf%m_leafn_xfer_to_litr1n(beg:end) = nan
    cnf%m_frootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livestemn_xfer_to_litr1n(beg:end) = nan
    cnf%m_deadstemn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livecrootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_deadcrootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livestemn_to_cwdn(beg:end) = nan
    cnf%m_deadstemn_to_cwdn(beg:end) = nan
    cnf%m_livecrootn_to_cwdn(beg:end) = nan
    cnf%m_deadcrootn_to_cwdn(beg:end) = nan
    cnf%m_retransn_to_litr1n(beg:end) = nan
    cnf%m_deadstemn_to_cwdn_fire(beg:end) = nan
    cnf%m_deadcrootn_to_cwdn_fire(beg:end) = nan
    cnf%m_litr1n_to_fire(beg:end) = nan
    cnf%m_litr2n_to_fire(beg:end) = nan
    cnf%m_litr3n_to_fire(beg:end) = nan
    cnf%m_cwdn_to_fire(beg:end) = nan
    cnf%leafn_to_litr1n(beg:end) = nan
    cnf%leafn_to_litr2n(beg:end) = nan
    cnf%leafn_to_litr3n(beg:end) = nan
    cnf%frootn_to_litr1n(beg:end) = nan
    cnf%frootn_to_litr2n(beg:end) = nan
    cnf%frootn_to_litr3n(beg:end) = nan
    cnf%cwdn_to_litr2n(beg:end) = nan
    cnf%cwdn_to_litr3n(beg:end) = nan
    cnf%litr1n_to_soil1n(beg:end) = nan
    cnf%sminn_to_soil1n_l1(beg:end) = nan
    cnf%litr2n_to_soil2n(beg:end) = nan
    cnf%sminn_to_soil2n_l2(beg:end) = nan
    cnf%litr3n_to_soil3n(beg:end) = nan
    cnf%sminn_to_soil3n_l3(beg:end) = nan
    cnf%soil1n_to_soil2n(beg:end) = nan
    cnf%sminn_to_soil2n_s1(beg:end) = nan
    cnf%soil2n_to_soil3n(beg:end) = nan
    cnf%sminn_to_soil3n_s2(beg:end) = nan
    cnf%soil3n_to_soil4n(beg:end) = nan
    cnf%sminn_to_soil4n_s3(beg:end) = nan
    cnf%soil4n_to_sminn(beg:end) = nan
    cnf%sminn_to_denit_l1s1(beg:end) = nan
    cnf%sminn_to_denit_l2s2(beg:end) = nan
    cnf%sminn_to_denit_l3s3(beg:end) = nan
    cnf%sminn_to_denit_s1s2(beg:end) = nan
    cnf%sminn_to_denit_s2s3(beg:end) = nan
    cnf%sminn_to_denit_s3s4(beg:end) = nan
    cnf%sminn_to_denit_s4(beg:end) = nan
    cnf%sminn_to_denit_excess(beg:end) = nan
    cnf%sminn_leached(beg:end) = nan
    cnf%potential_immob(beg:end) = nan
    cnf%actual_immob(beg:end) = nan
    cnf%sminn_to_plant(beg:end) = nan
    cnf%supplement_to_sminn(beg:end) = nan
    cnf%gross_nmin(beg:end) = nan
    cnf%net_nmin(beg:end) = nan
    cnf%denit(beg:end) = nan
    cnf%col_ninputs(beg:end) = nan
    cnf%col_noutputs(beg:end) = nan
    cnf%col_fire_nloss(beg:end) = nan

  end subroutine init_column_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_pstate_type
!
! !INTERFACE:
  subroutine init_landunit_pstate_type(beg, end, lps)
!
! !DESCRIPTION:
! Initialize landunit physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (landunit_pstate_type), intent(inout):: lps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

  ! currently nothing is here - just a place holder

  end subroutine init_landunit_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_dgvstate_type
!
! !INTERFACE:
  subroutine init_gridcell_dgvstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell DGVM variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_dgvstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gps%afirefrac(beg:end))
    allocate(gps%acfluxfire(beg:end))
    allocate(gps%bmfm(beg:end,maxpatch_pft))
    allocate(gps%afmicr(beg:end,maxpatch_pft))
    allocate(gps%begwater(beg:end))
    allocate(gps%endwater(beg:end))
    allocate(gps%begenergy(beg:end))
    allocate(gps%endenergy(beg:end))
    gps%afirefrac(beg:end) = nan
    gps%acfluxfire(beg:end) = nan
    gps%bmfm(beg:end,1:maxpatch_pft) = nan
    gps%afmicr(beg:end,1:maxpatch_pft) = nan
    gps%begwater(beg:end) = nan
    gps%endwater(beg:end) = nan
    gps%begenergy(beg:end) = nan
    gps%endenergy(beg:end) = nan

  end subroutine init_gridcell_dgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_pstate_type
!
! !INTERFACE:
  subroutine init_gridcell_pstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_pstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

  end subroutine init_gridcell_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wflux_type
!
! !INTERFACE:
  subroutine init_gridcell_wflux_type(beg, end, gwf)
!
! !DESCRIPTION:
! Initialize gridcell water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wflux_type), intent(inout):: gwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(gwf%qchan2(beg:end))
    allocate(gwf%qchocn2(beg:end))

    gwf%qchan2(beg:end) = 0._r8
    gwf%qchocn2(beg:end) = 0._r8

  end subroutine init_gridcell_wflux_type

end module clmtypeInitMod
