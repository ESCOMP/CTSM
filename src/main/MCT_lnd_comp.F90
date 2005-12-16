#include <misc.h>
#include <preproc.h>

module MCT_lnd_comp

#if (defined COUP_CAM)

  use seq_fields_mod

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use abortutils    , only : endrun
  use m_GlobalSegMap, only : GlobalSegMap
  use m_AttrVect    , only : AttrVect

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: MCT_lnd_AvInit
  public :: MCT_lnd_ExportInit
  public :: MCT_lnd_Export
  public :: MCT_lnd_Import

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  type(GlobalSegMap), public :: GSMap_lnd

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  integer :: index_l2c_Sl_t            ! temperature
  integer :: index_l2c_Sl_tref         ! 2m reference temperature
  integer :: index_l2c_Sl_qref         ! 2m reference specific humidity
  integer :: index_l2c_Sl_avsdr        ! albedo: direct , visible
  integer :: index_l2c_Sl_anidr        ! albedo: direct , near-ir
  integer :: index_l2c_Sl_avsdf        ! albedo: diffuse, visible
  integer :: index_l2c_Sl_anidf        ! albedo: diffuse, near-ir
  integer :: index_l2c_Sl_snowh        ! snow height
  integer :: index_l2c_Sl_landfrac     ! land fraction
  integer :: index_l2c_Fall_taux       ! wind stress, zonal
  integer :: index_l2c_Fall_tauy       ! wind stress, meridional
  integer :: index_l2c_Fall_lat        ! latent          heat flux
  integer :: index_l2c_Fall_sen        ! sensible        heat flux
  integer :: index_l2c_Fall_lwup       ! upward longwave heat flux
  integer :: index_l2c_Fall_evap       ! evaporation    water flux
  integer :: index_l2c_Fall_swnet      ! 2m reference temperature
  integer :: index_l2c_Fall_nee = 0    ! co2 flux **For testing set to 0

  integer :: index_c2l_Sa_z            ! bottom atm level height
  integer :: index_c2l_Sa_u            ! bottom atm level zon wind
  integer :: index_c2l_Sa_v            ! bottom atm level mer wind
  integer :: index_c2l_Sa_tbot         ! bottom atm level temp
  integer :: index_c2l_Sa_ptem         ! bottom atm level pot temp
  integer :: index_c2l_Sa_shum         ! bottom atm level spec hum
  integer :: index_c2l_Sa_dens         ! bottom atm level air dens
  integer :: index_c2l_Sa_pbot         ! bottom atm level pressure
  integer :: index_c2l_Sa_pslv         ! sea level atm pressure     
  integer :: index_c2l_Faxa_lwdn       ! downward longwave heat flux
  integer :: index_c2l_Faxa_rainc      ! precip: liquid, convective
  integer :: index_c2l_Faxa_rainl      ! precip: liquid, large-scale
  integer :: index_c2l_Faxa_snowc      ! precip: frozen, convective
  integer :: index_c2l_Faxa_snowl      ! precip: frozen, large-scale
  integer :: index_c2l_Faxa_swndr      ! shortwave: nir direct  down
  integer :: index_c2l_Faxa_swvdr      ! shortwave: vis direct  down
  integer :: index_c2l_Faxa_swndf      ! shortwave: nir diffuse down
  integer :: index_c2l_Faxa_swvdf      ! shortwave: vis diffuse down
  integer :: index_c2l_Sa_co2prog = 0  ! bottom atm prognostic co2 **For testing set to 0
  integer :: index_c2l_Sa_co2diag = 0  ! bottom atm diagnostic co2 **For testing set to 0

  integer :: index_a2c_Sa_z            ! bottom atm level height
  integer :: index_a2c_Sa_u            ! bottom atm level zon wind
  integer :: index_a2c_Sa_v            ! bottom atm level mer wind
  integer :: index_a2c_Sa_tbot         ! bottom atm level temp
  integer :: index_a2c_Sa_ptem         ! bottom atm level pot temp
  integer :: index_a2c_Sa_shum         ! bottom atm level spec hum
  integer :: index_a2c_Sa_dens         ! bottom atm level air den
  integer :: index_a2c_Sa_pbot         ! bottom atm level pressure
  integer :: index_a2c_Sa_pslv         ! sea level atm pressure
  integer :: index_a2c_Faxa_lwdn       ! downward lw heat flux
  integer :: index_a2c_Faxa_rainc      ! prec: liquid "convective"
  integer :: index_a2c_Faxa_rainl      ! prec: liquid "large scale"
  integer :: index_a2c_Faxa_snowc      ! prec: frozen "convective"
  integer :: index_a2c_Faxa_snowl      ! prec: frozen "large scale"
  integer :: index_a2c_Faxa_swndr      ! sw: nir direct  downward
  integer :: index_a2c_Faxa_swvdr      ! sw: vis direct  downward
  integer :: index_a2c_Faxa_swndf      ! sw: nir diffuse downward
  integer :: index_a2c_Faxa_swvdf      ! sw: vis diffuse downward
  integer :: index_a2c_Faxa_swnet      ! sw: net
  integer :: index_a2c_Sa_co2prog =0      ! bottom atm level prognostic co2
  integer :: index_a2c_Sa_co2diag =0     ! bottom atm level diagnostic co2


  integer :: index_r2c_Forr_roff   ! runoff to ocean

  integer, dimension(:), allocatable :: perm  ! permutation array to reorder points

  logical :: noland = .false.  ! Flag if no land points here

!===============================================================
contains
!===============================================================

  subroutine MCT_lnd_setIndices( )

    use shr_string_mod, only : shr_string_listGetIndexF, shr_string_listGetNum
    implicit none

    ! Determine send indices

    index_l2c_Sl_landfrac= shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_landfrac')
    index_l2c_Sl_t       = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_t')
    index_l2c_Sl_snowh   = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_snowh')
    index_l2c_Sl_avsdr   = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_avsdr')
    index_l2c_Sl_anidr   = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_anidr')
    index_l2c_Sl_avsdf   = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_avsdf')
    index_l2c_Sl_anidf   = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_anidf')
    index_l2c_Sl_tref    = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_tref')
    index_l2c_Sl_qref    = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Sl_qref')
    index_l2c_Fall_taux  = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Fall_taux')
    index_l2c_Fall_tauy  = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Fall_tauy')
    index_l2c_Fall_lat   = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Fall_lat')
    index_l2c_Fall_sen   = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Fall_sen')
    index_l2c_Fall_lwup  = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Fall_lwup')
    index_l2c_Fall_evap  = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Fall_evap')
    index_l2c_Fall_swnet = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Fall_swnet')
    index_l2c_Fall_nee   = shr_string_listGetIndexF(cpl_fields_l2c_fields,'Fall_nee')

    ! Determine receive indices

    index_c2l_Sa_z       = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_z')
    index_c2l_Sa_u       = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_u')
    index_c2l_Sa_v       = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_v')
    index_c2l_Sa_ptem    = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_ptem')
    index_c2l_Sa_shum    = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_shum')
    index_c2l_Sa_pbot    = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_pbot')
    index_c2l_Sa_tbot    = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_tbot')
    index_c2l_Faxa_lwdn  = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_lwdn')
    index_c2l_Faxa_rainc = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_rainc')
    index_c2l_Faxa_rainl = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_rainl')
    index_c2l_Faxa_snowc = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_snowc')
    index_c2l_Faxa_snowl = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_snowl')
    index_c2l_Faxa_swndr = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_swndr')
    index_c2l_Faxa_swvdr = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_swvdr')
    index_c2l_Faxa_swndf = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_swndf')
    index_c2l_Faxa_swvdf = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Faxa_swvdf')
    index_c2l_Sa_co2prog = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_co2prog')
    index_c2l_Sa_co2diag = shr_string_listGetIndexF(cpl_fields_c2l_fields,'Sa_co2diag')

    index_a2c_Sa_z          = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_z')
    index_a2c_Sa_u          = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_u')
    index_a2c_Sa_v          = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_v')
    index_a2c_Sa_tbot       = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_tbot')
    index_a2c_Sa_ptem       = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_ptem')
    index_a2c_Sa_pbot       = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_pbot')
    index_a2c_Sa_pslv       = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_pslv')
    index_a2c_Sa_shum       = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_shum')
    index_a2c_Sa_dens       = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_dens')
    index_a2c_Faxa_swnet    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_swnet')
    index_a2c_Faxa_lwdn     = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_lwdn')
    index_a2c_Faxa_rainc    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_rainc')
    index_a2c_Faxa_rainl    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_rainl')
    index_a2c_Faxa_snowc    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_snowc')
    index_a2c_Faxa_snowl    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_snowl')
    index_a2c_Faxa_swndr    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_swndr')
    index_a2c_Faxa_swvdr    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_swvdr')
    index_a2c_Faxa_swndf    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_swndf')
    index_a2c_Faxa_swvdf    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Faxa_swvdf')
    index_a2c_Sa_co2prog    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_co2prog')
    index_a2c_Sa_co2diag    = shr_string_listGetIndexF(cpl_fields_a2c_fields,'Sa_co2diag')


    ! Determine runoff indices

    index_r2c_Forr_roff  = shr_string_listGetIndexF(cpl_fields_r2c_fields,'Forr_roff')

  end subroutine MCT_lnd_setIndices

!==================================================================================

  subroutine MCT_lnd_AvInit( a2c_l, c2l_l, l2c_l )

    !-------------------------------------------------------------------
    use m_AttrVect    , only : MCT_aVect_init => init
    use m_AttrVect    , only : MCT_aVect_zero => zero 
    use m_GlobalSegMap, only : MCT_GSMap_lsize => lsize
    use m_MergeSorts  , only : MCT_IndexSet => IndexSet
    use m_MergeSorts  , only : MCT_IndexSort => IndexSort
    use m_Permuter    , only : MCT_Permute => Permute

    use MCT_seq
    use decompMod     , only : get_proc_bounds
#if (defined SPMD)
    use mpishorthand  , only : mpicom
#endif
    use clmtype
    use clm_varpar    , only : lsmlon, lsmlat
    use decompMod     , only : get_proc_bounds


    implicit none
    type(AttrVect), intent(inout) :: a2c_l
    type(AttrVect), intent(inout) :: c2l_l
    type(AttrVect), intent(inout) :: l2c_l

#if !(defined SPMD)
    integer, parameter :: mpicom = 1
#endif

    integer,allocatable :: gindex(:)
    integer :: i, j, n, gi
    integer :: lsize,gsize
    integer :: ier
    integer :: begp, endp    ! beginning and ending pft indices
    integer :: begc, endc    ! beginning and ending column indices
    integer :: begl, endl    ! beginning and ending landunit indices
    integer :: begg, endg    ! beginning and ending gridcell indices
    !-------------------------------------------------------------------

    call MCT_lnd_SetIndices( )

    ! Build the land grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used
    ! in SCRIP
    
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    allocate(gindex(begg:endg),stat=ier)

    ! number the local grid

    do n = begg, endg
        i = clm3%g%ixy(n)
        j = clm3%g%jxy(n)
        gi = (j-1)*lsmlon + i
        gindex(n) = gi
    end do

    lsize = endg-begg+1
    gsize = lsmlon*lsmlat

    allocate(perm(lsize),stat=ier)

    ! reorder gindex to be in ascending order.

    ! initialize a permutation array

    call MCT_IndexSet(perm)

    ! derive a permutation that puts gindex in ascending order
    ! the default for IndexSort is Ascending.

    call MCT_IndexSort(lsize,perm,gindex)

    ! Sort gindex in-place

    call MCT_Permute(gindex,perm,lsize)

    call SetGSMap(lsize, gsize, gindex, mpicom, LNDID, GSMap_lnd)

    call MCT_aVect_init(a2c_l, rList=cpl_fields_a2c_fields,    &
         lsize=MCT_GSMap_lsize(GSMap_lnd, mpicom))

    call MCT_aVect_init(c2l_l, rList=cpl_fields_c2l_fields,    &
         lsize=MCT_GSMap_lsize(GSMap_lnd, mpicom))

    call MCT_aVect_init(l2c_l, rList=cpl_fields_l2c_fields,    &
         lsize=MCT_GSMap_lsize(GSMap_lnd, mpicom))

    call MCT_aVect_zero(a2c_l)
    call MCT_aVect_zero(c2l_l)
    call MCT_aVect_zero(l2c_l)

    deallocate(gindex)

  end subroutine MCT_lnd_AvInit

!=================================================================================

  subroutine MCT_lnd_ExportInit( l2as, l2af, l2c_l )

    use clmtype   , only : lnd2atm_state_type, lnd2atm_flux_type, clm3
    use clm_varsur, only : landfrac
    use clm_varcon, only : sb
    use decompMod , only : get_proc_bounds
    use m_AttrVect, only : MCT_Avt_Permute => Permute

    implicit none

    type(lnd2atm_state_type), intent(inout) :: l2as
    type(lnd2atm_flux_type) , intent(inout) :: l2af
    type(AttrVect)          , intent(inout) :: l2c_l

    integer :: g,i,ier,nstep
    integer :: begp, endp    ! beginning and ending pft indices
    integer :: begc, endc    ! beginning and ending column indices
    integer :: begl, endl    ! beginning and ending landunit indices
    integer :: begg, endg    ! beginning and ending gridcell indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    i=1
    do g = begg,endg
       l2c_l%rAttr(index_l2c_Sl_landfrac,i) =  landfrac(clm3%g%ixy(g), clm3%g%jxy(g))
       l2c_l%rAttr(index_l2c_Sl_t,i)        =  sqrt(sqrt(l2af%eflx_lwrad_out(g)/sb))
       l2c_l%rAttr(index_l2c_Sl_snowh,i)    =  l2as%h2osno(g)
       l2c_l%rAttr(index_l2c_Sl_avsdr,i)    =  l2as%albd(g,1)
       l2c_l%rAttr(index_l2c_Sl_anidr,i)    =  l2as%albd(g,2)
       l2c_l%rAttr(index_l2c_Sl_avsdf,i)    =  l2as%albi(g,1)
       l2c_l%rAttr(index_l2c_Sl_anidf,i)    =  l2as%albi(g,2)
       l2c_l%rAttr(index_l2c_Fall_lwup,i)   =  l2af%eflx_lwrad_out(g)
       i=i+1
    end do

    ! permute before using the Rearrange call

    call MCT_Avt_Permute(l2c_l,perm)
    
  end subroutine MCT_lnd_ExportInit

!====================================================================================

  subroutine MCT_lnd_Export( l2as, l2af, l2c_l )   

    use clmtype   , only : lnd2atm_state_type, lnd2atm_flux_type
    use decompMod , only : get_proc_bounds
    use m_AttrVect, only : MCT_Avt_Permute => Permute
    implicit none

    type(lnd2atm_state_type), intent(inout) :: l2as
    type(lnd2atm_flux_type) , intent(inout) :: l2af
    type(AttrVect)          , intent(inout) :: l2c_l

    integer :: g,i
    integer  :: begp, endp    ! beginning and ending pft indices
    integer  :: begc, endc    ! beginning and ending column indices
    integer  :: begl, endl    ! beginning and ending landunit indices
    integer  :: begg, endg    ! beginning and ending gridcell indices
    
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    i=1
    l2c_l%rAttr(:,:) = 0.0_r8
    do g = begg,endg
       l2c_l%rAttr(index_l2c_Sl_t,i)       =  l2as%t_rad(g)
       l2c_l%rAttr(index_l2c_Sl_snowh,i)   =  l2as%h2osno(g)
       l2c_l%rAttr(index_l2c_Sl_avsdr,i)   =  l2as%albd(g,1)
       l2c_l%rAttr(index_l2c_Sl_anidr,i)   =  l2as%albd(g,2)
       l2c_l%rAttr(index_l2c_Sl_avsdf,i)   =  l2as%albi(g,1)
       l2c_l%rAttr(index_l2c_Sl_anidf,i)   =  l2as%albi(g,2)
       l2c_l%rAttr(index_l2c_Sl_tref,i)    =  l2as%t_ref2m(g)
       l2c_l%rAttr(index_l2c_Sl_qref,i)    =  l2as%q_ref2m(g)
       l2c_l%rAttr(index_l2c_Fall_taux,i)  =  l2af%taux(g)
       l2c_l%rAttr(index_l2c_Fall_tauy,i)  =  l2af%tauy(g)
       l2c_l%rAttr(index_l2c_Fall_lat,i)   =  l2af%eflx_lh_tot(g)
       l2c_l%rAttr(index_l2c_Fall_sen,i)   =  l2af%eflx_sh_tot(g)
       l2c_l%rAttr(index_l2c_Fall_lwup,i)  =  l2af%eflx_lwrad_out(g)
       l2c_l%rAttr(index_l2c_Fall_evap,i)  =  l2af%qflx_evap_tot(g)
       l2c_l%rAttr(index_l2c_Fall_swnet,i) =  l2af%fsa(g)
       if (index_l2c_Fall_nee /= 0) then
          l2c_l%rAttr(index_l2c_Fall_nee,i) = l2af%nee(g)
       end if
       i=i+1
    end do

    ! permute before using the Rearrange call.

    call MCT_Avt_Permute(l2c_l,perm)

  end subroutine MCT_lnd_Export

!====================================================================================

  subroutine MCT_lnd_Import(a2ls, a2lf, a2c_l)

    use clmtype   , only : atm2lnd_state_type, atm2lnd_flux_type
    use clm_varctl, only : co2_type
    use clm_varcon, only : rair, o2_molar_const, co2_ppmv_const, c13ratio
    use m_AttrVectComms , only: AttrVect_gather => gather
    use m_AttrVect, only : MCT_Avt_lsize => lsize
    use m_AttrVect, only : MCT_Avt_Unpermute => Unpermute
    use decompMod , only : get_proc_bounds
    use time_manager,   only: get_nstep
    use spmd_utils    , only: masterproc, iam
#if (defined SPMD)
    use mpishorthand  , only: mpicom
#endif

    implicit none
    type(atm2lnd_state_type), intent(inout) :: a2ls
    type(atm2lnd_flux_type) , intent(inout) :: a2lf
    type(AttrVect)          , intent(inout) :: a2c_l

    type(AttrVect) :: Ga2c_l
    integer  :: g,i,nstep,ier
    real(r8) :: forc_rainc    ! rainxy Atm flux mm/s
    real(r8) :: forc_rainl    ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc    ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl    ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag ! temporary
    real(r8) :: co2_ppmv_prog ! temporary
    real(r8) :: co2_ppmv      ! temporary
    integer  :: begp, endp    ! beginning and ending pft indices
    integer  :: begc, endc    ! beginning and ending column indices
    integer  :: begl, endl    ! beginning and ending landunit indices
    integer  :: begg, endg    ! beginning and ending gridcell indices

    ! unpermute after rearrange call and before copying into local arrays.

    call MCT_Avt_Unpermute(a2c_l,perm)

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    if (co2_type == 'prognostic' .and. index_a2c_Sa_co2prog == 0) then
       write(6,*)' must have nonzero index_a2c_Sa_co2prog for co2_type equal to prognostic'
       call endrun()
    else if (co2_type == 'diagnostic' .and. index_a2c_Sa_co2diag == 0) then
       write(6,*)' must have nonzero index_a2c_Sa_co2diag for co2_type equal to diagnostic'
       call endrun()
    end if
    
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.
    
    i=1 
    do g = begg,endg
       
        ! Determine required receive fields

        a2ls%forc_hgt(g)     = a2c_l%rAttr(index_a2c_Sa_z,i)         ! zgcmxy  Atm state m
        a2ls%forc_u(g)       = a2c_l%rAttr(index_a2c_Sa_u,i)         ! forc_uxy  Atm state m/s
        a2ls%forc_v(g)       = a2c_l%rAttr(index_a2c_Sa_v,i)         ! forc_vxy  Atm state m/s
        a2ls%forc_th(g)      = a2c_l%rAttr(index_a2c_Sa_ptem,i)      ! forc_thxy Atm state K
        a2ls%forc_q(g)       = a2c_l%rAttr(index_a2c_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
        a2ls%forc_pbot(g)    = a2c_l%rAttr(index_a2c_Sa_pbot,i)      ! ptcmxy  Atm state Pa
        a2ls%forc_t(g)       = a2c_l%rAttr(index_a2c_Sa_tbot,i)      ! forc_txy  Atm state K
        a2lf%forc_lwrad(g)   = a2c_l%rAttr(index_a2c_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
        forc_rainc           = a2c_l%rAttr(index_a2c_Faxa_rainc,i)   ! mm/s
        forc_rainl           = a2c_l%rAttr(index_a2c_Faxa_rainl,i)   ! mm/s
        forc_snowc           = a2c_l%rAttr(index_a2c_Faxa_snowc,i)   ! mm/s
        forc_snowl           = a2c_l%rAttr(index_a2c_Faxa_snowl,i)   ! mm/s
        a2lf%forc_solad(g,2) = a2c_l%rAttr(index_a2c_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
        a2lf%forc_solad(g,1) = a2c_l%rAttr(index_a2c_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
        a2lf%forc_solai(g,2) = a2c_l%rAttr(index_a2c_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
        a2lf%forc_solai(g,1) = a2c_l%rAttr(index_a2c_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

        ! Determine optional receive fields

        if (index_a2c_Sa_co2prog /= 0) then
           co2_ppmv_prog = a2c_l%rAttr(index_a2c_Sa_co2prog,i)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv_const
        end if
 
        if (index_a2c_Sa_co2diag /= 0) then
           co2_ppmv_diag = a2c_l%rAttr(index_a2c_Sa_co2diag,i)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv_const
        end if

        ! Determine derived quantities for required fields

        a2ls%forc_hgt_u(g) = a2ls%forc_hgt(g)    !observational height of wind [m]
        a2ls%forc_hgt_t(g) = a2ls%forc_hgt(g)    !observational height of temperature [m]
        a2ls%forc_hgt_q(g) = a2ls%forc_hgt(g)    !observational height of humidity [m]
        a2ls%forc_vp(g)    = a2ls%forc_q(g) * a2ls%forc_pbot(g) &
                             / (0.622_r8 + 0.378_r8 * a2ls%forc_q(g))
        a2ls%forc_rho(g)   = (a2ls%forc_pbot(g) - 0.378_r8 * a2ls%forc_vp(g)) &
                             / (rair * a2ls%forc_t(g))
        a2ls%forc_po2(g)   = o2_molar_const * a2ls%forc_pbot(g)
        a2ls%forc_wind(g)  = sqrt(a2ls%forc_u(g)**2 + a2ls%forc_v(g)**2)
        a2lf%forc_solar(g) = a2lf%forc_solad(g,1) + a2lf%forc_solai(g,1) + &
                             a2lf%forc_solad(g,2) + a2lf%forc_solai(g,2)
        a2lf%forc_rain(g)  = forc_rainc + forc_rainl
        a2lf%forc_snow(g)  = forc_snowc + forc_snowl
        
        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        if (co2_type == 'prognostic') then
           co2_ppmv = co2_ppmv_prog
        else if (co2_type == 'diagnostic') then
           co2_ppmv = co2_ppmv_diag 
        else
           co2_ppmv = co2_ppmv_const      
        end if
        a2ls%forc_pco2(g) = co2_ppmv * 1.e-6_r8 * a2ls%forc_pbot(g) 
        a2ls%forc_pc13o2(g) = co2_ppmv * c13ratio * 1.e-6_r8 * a2ls%forc_pbot(g)
	 
        i=i+1

     end do

     ! debug write statements (remove)

     print *,'co2_type = ', co2_type, ' co2_ppmv = ', co2_ppmv

   end subroutine MCT_lnd_Import

#endif 

end module MCT_lnd_comp

