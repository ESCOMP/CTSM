module clm_cpl_indices
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_cpl_indices
!
! !DESCRIPTION:
!    Module containing the indices for the fields passed between CLM and
!    the driver. Includes the River Transport Model fields (RTM) and the
!    fields needed by the land-ice component (sno).
!
! !USES:

  implicit none

  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_cpl_indices_set        ! Set the coupler indices

!
! !PUBLIC DATA MEMBERS:
!
  ! lnd -> drv (required)

  integer, public ::index_l2x_Sl_t            ! temperature
  integer, public ::index_l2x_Sl_tref         ! 2m reference temperature
  integer, public ::index_l2x_Sl_qref         ! 2m reference specific humidity
  integer, public ::index_l2x_Sl_avsdr        ! albedo: direct , visible
  integer, public ::index_l2x_Sl_anidr        ! albedo: direct , near-ir
  integer, public ::index_l2x_Sl_avsdf        ! albedo: diffuse, visible
  integer, public ::index_l2x_Sl_anidf        ! albedo: diffuse, near-ir
  integer, public ::index_l2x_Sl_snowh        ! snow height
  integer, public ::index_l2x_Sl_u10          ! 10m wind
  integer, public ::index_l2x_Sl_ddvel        ! dry deposition velocities (optional)
  integer, public ::index_l2x_Sl_fv           ! friction velocity  
  integer, public ::index_l2x_Sl_ram1         ! aerodynamical resistance
  integer, public ::index_l2x_Fall_taux       ! wind stress, zonal
  integer, public ::index_l2x_Fall_tauy       ! wind stress, meridional
  integer, public ::index_l2x_Fall_lat        ! latent          heat flux
  integer, public ::index_l2x_Fall_sen        ! sensible        heat flux
  integer, public ::index_l2x_Fall_lwup       ! upward longwave heat flux
  integer, public ::index_l2x_Fall_evap       ! evaporation     water flux
  integer, public ::index_l2x_Fall_swnet      ! heat flux       shortwave net       
  integer, public ::index_l2x_Fall_fco2_lnd   ! co2 flux **For testing set to 0
  integer, public ::index_l2x_Fall_flxdst1    ! dust flux size bin 1    
  integer, public ::index_l2x_Fall_flxdst2    ! dust flux size bin 2    
  integer, public ::index_l2x_Fall_flxdst3    ! dust flux size bin 3    
  integer, public ::index_l2x_Fall_flxdst4    ! dust flux size bin 4
  integer, public ::index_l2x_Fall_flxvoc1    ! voc flux size bin 1    
  integer, public ::index_l2x_Fall_flxvoc2    ! voc flux size bin 2    
  integer, public ::index_l2x_Fall_flxvoc3    ! voc flux size bin 3    
  integer, public ::index_l2x_Fall_flxvoc4    ! voc flux size bin 4
  integer, public ::index_l2x_Fall_flxvoc5    ! voc flux size bin 5

  integer, public :: nflds_l2x = 0

  ! roff to driver (part of land for now) (optional if RTM is off)

  integer, public ::index_r2x_Forr_roff = 0   ! liquid runoff to ocean
  integer, public ::index_r2x_Forr_ioff = 0   ! ice runoff to ocean

  integer, public :: nflds_r2x = 0

  ! drv -> lnd (required)

  integer, public ::index_x2l_Sa_z            ! bottom atm level height
  integer, public ::index_x2l_Sa_u            ! bottom atm level zon wind
  integer, public ::index_x2l_Sa_v            ! bottom atm level mer wind
  integer, public ::index_x2l_Sa_ptem         ! bottom atm level pot temp
  integer, public ::index_x2l_Sa_shum         ! bottom atm level spec hum
  integer, public ::index_x2l_Sa_pbot         ! bottom atm level pressure
  integer, public ::index_x2l_Sa_tbot         ! bottom atm level temp
  integer, public ::index_x2l_Faxa_lwdn       ! downward lw heat flux
  integer, public ::index_x2l_Faxa_rainc      ! prec: liquid "convective"
  integer, public ::index_x2l_Faxa_rainl      ! prec: liquid "large scale"
  integer, public ::index_x2l_Faxa_snowc      ! prec: frozen "convective"
  integer, public ::index_x2l_Faxa_snowl      ! prec: frozen "large scale"
  integer, public ::index_x2l_Faxa_swndr      ! sw: nir direct  downward
  integer, public ::index_x2l_Faxa_swvdr      ! sw: vis direct  downward
  integer, public ::index_x2l_Faxa_swndf      ! sw: nir diffuse downward
  integer, public ::index_x2l_Faxa_swvdf      ! sw: vis diffuse downward
  integer, public ::index_x2l_Sa_co2prog      ! bottom atm level prognostic co2
  integer, public ::index_x2l_Sa_co2diag      ! bottom atm level diagnostic co2
  integer, public ::index_x2l_Faxa_bcphidry   ! flux: Black Carbon hydrophilic dry deposition
  integer, public ::index_x2l_Faxa_bcphodry   ! flux: Black Carbon hydrophobic dry deposition
  integer, public ::index_x2l_Faxa_bcphiwet   ! flux: Black Carbon hydrophilic wet deposition
  integer, public ::index_x2l_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer, public ::index_x2l_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer, public ::index_x2l_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer, public ::index_x2l_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer, public ::index_x2l_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer, public ::index_x2l_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer, public ::index_x2l_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer, public ::index_x2l_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer, public ::index_x2l_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer, public ::index_x2l_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer, public ::index_x2l_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition

  integer, public :: nflds_x2l = 0

  ! sno -> drv (only if land-ice model is NOT a stub model)
  integer, public ::index_s2x_Ss_tsrf01   = 0   ! glc MEC temperature
  integer, public ::index_s2x_Ss_tsrf02   = 0
  integer, public ::index_s2x_Ss_tsrf03   = 0
  integer, public ::index_s2x_Ss_tsrf04   = 0
  integer, public ::index_s2x_Ss_tsrf05   = 0
  integer, public ::index_s2x_Ss_tsrf06   = 0
  integer, public ::index_s2x_Ss_tsrf07   = 0
  integer, public ::index_s2x_Ss_tsrf08   = 0
  integer, public ::index_s2x_Ss_tsrf09   = 0
  integer, public ::index_s2x_Ss_tsrf10   = 0
  integer, public ::index_s2x_Ss_topo01   = 0   ! glc MEC topo height
  integer, public ::index_s2x_Ss_topo02   = 0
  integer, public ::index_s2x_Ss_topo03   = 0
  integer, public ::index_s2x_Ss_topo04   = 0
  integer, public ::index_s2x_Ss_topo05   = 0
  integer, public ::index_s2x_Ss_topo06   = 0
  integer, public ::index_s2x_Ss_topo07   = 0
  integer, public ::index_s2x_Ss_topo08   = 0
  integer, public ::index_s2x_Ss_topo09   = 0
  integer, public ::index_s2x_Ss_topo10   = 0
  integer, public ::index_s2x_Fgss_qice01 = 0   ! glc MEC ice flux
  integer, public ::index_s2x_Fgss_qice02 = 0
  integer, public ::index_s2x_Fgss_qice03 = 0
  integer, public ::index_s2x_Fgss_qice04 = 0
  integer, public ::index_s2x_Fgss_qice05 = 0
  integer, public ::index_s2x_Fgss_qice06 = 0
  integer, public ::index_s2x_Fgss_qice07 = 0
  integer, public ::index_s2x_Fgss_qice08 = 0
  integer, public ::index_s2x_Fgss_qice09 = 0
  integer, public ::index_s2x_Fgss_qice10 = 0

  integer, public :: nflds_s2x = 0

  ! drv -> sno (only if land-ice model is NOT a stub model)

  integer, public ::index_x2s_Sg_frac01   = 0   ! Fraction of glacier in glc MEC class 1
  integer, public ::index_x2s_Sg_topo01   = 0   ! Topo height in glc MEC class 1
  integer, public ::index_x2s_Fsgg_rofi01 = 0   ! Ice runoff from glc model
  integer, public ::index_x2s_Fsgg_rofl01 = 0   ! Liquid runoff from glc model
  integer, public ::index_x2s_Fsgg_hflx01 = 0
  integer, public ::index_x2s_Sg_frac02   = 0
  integer, public ::index_x2s_Sg_topo02   = 0
  integer, public ::index_x2s_Fsgg_rofi02 = 0
  integer, public ::index_x2s_Fsgg_rofl02 = 0
  integer, public ::index_x2s_Fsgg_hflx02 = 0
  integer, public ::index_x2s_Sg_frac03   = 0
  integer, public ::index_x2s_Sg_topo03   = 0
  integer, public ::index_x2s_Fsgg_rofi03 = 0
  integer, public ::index_x2s_Fsgg_rofl03 = 0
  integer, public ::index_x2s_Fsgg_hflx03 = 0
  integer, public ::index_x2s_Sg_frac04   = 0
  integer, public ::index_x2s_Sg_topo04   = 0
  integer, public ::index_x2s_Fsgg_rofi04 = 0
  integer, public ::index_x2s_Fsgg_rofl04 = 0
  integer, public ::index_x2s_Fsgg_hflx04 = 0
  integer, public ::index_x2s_Sg_frac05   = 0
  integer, public ::index_x2s_Sg_topo05   = 0
  integer, public ::index_x2s_Fsgg_rofi05 = 0
  integer, public ::index_x2s_Fsgg_rofl05 = 0
  integer, public ::index_x2s_Fsgg_hflx05 = 0
  integer, public ::index_x2s_Sg_frac06   = 0
  integer, public ::index_x2s_Sg_topo06   = 0
  integer, public ::index_x2s_Fsgg_rofi06 = 0
  integer, public ::index_x2s_Fsgg_rofl06 = 0
  integer, public ::index_x2s_Fsgg_hflx06 = 0
  integer, public ::index_x2s_Sg_frac07   = 0
  integer, public ::index_x2s_Sg_topo07   = 0
  integer, public ::index_x2s_Fsgg_rofi07 = 0
  integer, public ::index_x2s_Fsgg_rofl07 = 0
  integer, public ::index_x2s_Fsgg_hflx07 = 0
  integer, public ::index_x2s_Sg_frac08   = 0
  integer, public ::index_x2s_Sg_topo08   = 0
  integer, public ::index_x2s_Fsgg_rofi08 = 0
  integer, public ::index_x2s_Fsgg_rofl08 = 0
  integer, public ::index_x2s_Fsgg_hflx08 = 0
  integer, public ::index_x2s_Sg_frac09   = 0
  integer, public ::index_x2s_Sg_topo09   = 0
  integer, public ::index_x2s_Fsgg_rofi09 = 0
  integer, public ::index_x2s_Fsgg_rofl09 = 0
  integer, public ::index_x2s_Fsgg_hflx09 = 0
  integer, public ::index_x2s_Sg_frac10   = 0
  integer, public ::index_x2s_Sg_topo10   = 0
  integer, public ::index_x2s_Fsgg_rofi10 = 0
  integer, public ::index_x2s_Fsgg_rofl10 = 0
  integer, public ::index_x2s_Fsgg_hflx10 = 0

  integer, public :: nflds_x2s = 0
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 01/19/2011, Erik Kluzek:         Added protex headers
! 01/20/2011, Erik Kluzek:         Set nflds
!
!EOP
!-----------------------------------------------------------------------

!=======================================================================
contains
!=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_cpl_indices_set
!
! !CALLED FROM: lnd_comp_mct or lnd_comp_esmf
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
! !INTERFACE:
  subroutine clm_cpl_indices_set( )

!
! !DESCRIPTION: 
! Set the coupler indices needed by the land model coupler
! interface.
!
! !USES:
  use seq_flds_mod  , only: seq_flds_x2l_fields, seq_flds_l2x_fields,     &
                            seq_flds_x2s_fields, seq_flds_s2x_fields,     &
                                                 seq_flds_r2x_fields
  use mct_mod       , only: mct_aVect, mct_aVect_init, mct_avect_indexra, &
                            mct_aVect_clean, mct_avect_nRattr
  use seq_drydep_mod, only: drydep_fields_token, lnd_drydep
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 01/2011, Erik Kluzek:         Added protex headers
!
! !LOCAL VARIABLES:
    type(mct_aVect) :: l2x      ! temporary, land to coupler
    type(mct_aVect) :: x2l      ! temporary, coupler to land
    type(mct_aVect) :: r2x      ! temporary, runoff to coupler
    type(mct_aVect) :: s2x      ! temporary, glacier to coupler
    type(mct_aVect) :: x2s      ! temporary, coupler to glacier
    character(len=32) :: subname = 'clm_cpl_indices_set'  ! subroutine name
!EOP
!
!-----------------------------------------------------------------------

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2l, rList=seq_flds_x2l_fields, lsize=1)
    call mct_aVect_init(l2x, rList=seq_flds_l2x_fields, lsize=1)

    ! lnd -> drv 
    index_l2x_Sl_t          = mct_avect_indexra(l2x,'Sl_t')
    index_l2x_Sl_snowh      = mct_avect_indexra(l2x,'Sl_snowh')
    index_l2x_Sl_avsdr      = mct_avect_indexra(l2x,'Sl_avsdr')
    index_l2x_Sl_anidr      = mct_avect_indexra(l2x,'Sl_anidr')
    index_l2x_Sl_avsdf      = mct_avect_indexra(l2x,'Sl_avsdf')
    index_l2x_Sl_anidf      = mct_avect_indexra(l2x,'Sl_anidf')
    index_l2x_Sl_tref       = mct_avect_indexra(l2x,'Sl_tref')
    index_l2x_Sl_qref       = mct_avect_indexra(l2x,'Sl_qref')
    index_l2x_Sl_u10        = mct_avect_indexra(l2x,'Sl_u10')
    index_l2x_Fall_taux     = mct_avect_indexra(l2x,'Fall_taux')
    index_l2x_Fall_tauy     = mct_avect_indexra(l2x,'Fall_tauy')
    index_l2x_Fall_lat      = mct_avect_indexra(l2x,'Fall_lat')
    index_l2x_Fall_sen      = mct_avect_indexra(l2x,'Fall_sen')
    index_l2x_Fall_lwup     = mct_avect_indexra(l2x,'Fall_lwup')
    index_l2x_Fall_evap     = mct_avect_indexra(l2x,'Fall_evap')
    index_l2x_Fall_swnet    = mct_avect_indexra(l2x,'Fall_swnet')
    index_l2x_Sl_ram1       = mct_avect_indexra(l2x,'Sl_ram1')
    index_l2x_Sl_fv         = mct_avect_indexra(l2x,'Sl_fv')
    index_l2x_Fall_flxdst1  = mct_avect_indexra(l2x,'Fall_flxdst1')
    index_l2x_Fall_flxdst2  = mct_avect_indexra(l2x,'Fall_flxdst2')
    index_l2x_Fall_flxdst3  = mct_avect_indexra(l2x,'Fall_flxdst3')
    index_l2x_Fall_flxdst4  = mct_avect_indexra(l2x,'Fall_flxdst4')
    index_l2x_Fall_flxvoc1  = mct_avect_indexra(l2x,'Fall_flxvoc1' ,perrwith='quiet')
    index_l2x_Fall_flxvoc2  = mct_avect_indexra(l2x,'Fall_flxvoc2' ,perrwith='quiet')
    index_l2x_Fall_flxvoc3  = mct_avect_indexra(l2x,'Fall_flxvoc3' ,perrwith='quiet')
    index_l2x_Fall_flxvoc4  = mct_avect_indexra(l2x,'Fall_flxvoc4' ,perrwith='quiet')
    index_l2x_Fall_flxvoc5  = mct_avect_indexra(l2x,'Fall_flxvoc5' ,perrwith='quiet')
    index_l2x_Fall_fco2_lnd = mct_avect_indexra(l2x,'Fall_fco2_lnd',perrwith='quiet')
    if ( lnd_drydep )then
       index_l2x_Sl_ddvel = mct_avect_indexra(l2x, trim(drydep_fields_token))
    else
       index_l2x_Sl_ddvel = 0
    end if

    nflds_l2x = mct_avect_nRattr(l2x)

    ! drv -> lnd
    index_x2l_Sa_z          = mct_avect_indexra(x2l,'Sa_z')
    index_x2l_Sa_u          = mct_avect_indexra(x2l,'Sa_u')
    index_x2l_Sa_v          = mct_avect_indexra(x2l,'Sa_v')
    index_x2l_Sa_ptem       = mct_avect_indexra(x2l,'Sa_ptem')
    index_x2l_Sa_pbot       = mct_avect_indexra(x2l,'Sa_pbot')
    index_x2l_Sa_tbot       = mct_avect_indexra(x2l,'Sa_tbot')
    index_x2l_Sa_shum       = mct_avect_indexra(x2l,'Sa_shum')
    index_x2l_Faxa_lwdn     = mct_avect_indexra(x2l,'Faxa_lwdn')
    index_x2l_Faxa_rainc    = mct_avect_indexra(x2l,'Faxa_rainc')
    index_x2l_Faxa_rainl    = mct_avect_indexra(x2l,'Faxa_rainl')
    index_x2l_Faxa_snowc    = mct_avect_indexra(x2l,'Faxa_snowc')
    index_x2l_Faxa_snowl    = mct_avect_indexra(x2l,'Faxa_snowl')
    index_x2l_Faxa_swndr    = mct_avect_indexra(x2l,'Faxa_swndr')
    index_x2l_Faxa_swvdr    = mct_avect_indexra(x2l,'Faxa_swvdr')
    index_x2l_Faxa_swndf    = mct_avect_indexra(x2l,'Faxa_swndf')
    index_x2l_Faxa_swvdf    = mct_avect_indexra(x2l,'Faxa_swvdf')
    index_x2l_Faxa_bcphidry = mct_avect_indexra(x2l,'Faxa_bcphidry')
    index_x2l_Faxa_bcphodry = mct_avect_indexra(x2l,'Faxa_bcphodry')
    index_x2l_Faxa_bcphiwet = mct_avect_indexra(x2l,'Faxa_bcphiwet')
    index_x2l_Faxa_ocphidry = mct_avect_indexra(x2l,'Faxa_ocphidry')
    index_x2l_Faxa_ocphodry = mct_avect_indexra(x2l,'Faxa_ocphodry')
    index_x2l_Faxa_ocphiwet = mct_avect_indexra(x2l,'Faxa_ocphiwet')
    index_x2l_Faxa_dstdry1  = mct_avect_indexra(x2l,'Faxa_dstdry1')
    index_x2l_Faxa_dstdry2  = mct_avect_indexra(x2l,'Faxa_dstdry2')
    index_x2l_Faxa_dstdry3  = mct_avect_indexra(x2l,'Faxa_dstdry3')
    index_x2l_Faxa_dstdry4  = mct_avect_indexra(x2l,'Faxa_dstdry4')
    index_x2l_Faxa_dstwet1  = mct_avect_indexra(x2l,'Faxa_dstwet1')
    index_x2l_Faxa_dstwet2  = mct_avect_indexra(x2l,'Faxa_dstwet2')
    index_x2l_Faxa_dstwet3  = mct_avect_indexra(x2l,'Faxa_dstwet3')
    index_x2l_Faxa_dstwet4  = mct_avect_indexra(x2l,'Faxa_dstwet4')
    index_x2l_Sa_co2prog    = mct_avect_indexra(x2l,'Sa_co2prog',perrwith='quiet')
    index_x2l_Sa_co2diag    = mct_avect_indexra(x2l,'Sa_co2diag',perrwith='quiet')

    nflds_x2l = mct_avect_nRattr(x2l)

    call mct_aVect_clean(x2l)
    call mct_aVect_clean(l2x)

    ! runoff

    call mct_aVect_init(r2x, rList=seq_flds_r2x_fields, lsize=1)

    index_r2x_Forr_roff  = mct_avect_indexra(r2x,'Forr_roff')
    index_r2x_Forr_ioff  = mct_avect_indexra(r2x,'Forr_ioff')

    nflds_r2x = mct_avect_nRattr(r2x)

    call mct_aVect_clean(r2x)

    ! sno (for land-ice model)

    call mct_aVect_init(x2s, rList=seq_flds_x2s_fields, lsize=1)
    call mct_aVect_init(s2x, rList=seq_flds_s2x_fields, lsize=1)

#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
    index_x2s_Sg_frac01   = mct_avect_indexra(x2s,'Sg_frac01')
    index_x2s_Sg_topo01   = mct_avect_indexra(x2s,'Sg_topo01')
    index_x2s_Fsgg_rofi01 = mct_avect_indexra(x2s,'Fsgg_rofi01')
    index_x2s_Fsgg_rofl01 = mct_avect_indexra(x2s,'Fsgg_rofl01')
    index_x2s_Fsgg_hflx01 = mct_avect_indexra(x2s,'Fsgg_hflx01')
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
    index_x2s_Sg_frac02   = mct_avect_indexra(x2s,'Sg_frac02')
    index_x2s_Sg_topo02   = mct_avect_indexra(x2s,'Sg_topo02')
    index_x2s_Sg_frac03   = mct_avect_indexra(x2s,'Sg_frac03')
    index_x2s_Sg_topo03   = mct_avect_indexra(x2s,'Sg_topo03')
    index_x2s_Fsgg_rofi02 = mct_avect_indexra(x2s,'Fsgg_rofi02')
    index_x2s_Fsgg_rofl02 = mct_avect_indexra(x2s,'Fsgg_rofl02')
    index_x2s_Fsgg_hflx02 = mct_avect_indexra(x2s,'Fsgg_hflx02')
    index_x2s_Fsgg_rofi03 = mct_avect_indexra(x2s,'Fsgg_rofi03')
    index_x2s_Fsgg_rofl03 = mct_avect_indexra(x2s,'Fsgg_rofl03')
    index_x2s_Fsgg_hflx03 = mct_avect_indexra(x2s,'Fsgg_hflx03')
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
    index_x2s_Sg_frac04   = mct_avect_indexra(x2s,'Sg_frac04')
    index_x2s_Sg_topo04   = mct_avect_indexra(x2s,'Sg_topo04')
    index_x2s_Sg_frac05   = mct_avect_indexra(x2s,'Sg_frac05')
    index_x2s_Sg_topo05   = mct_avect_indexra(x2s,'Sg_topo05')
    index_x2s_Fsgg_rofi04 = mct_avect_indexra(x2s,'Fsgg_rofi04')
    index_x2s_Fsgg_rofl04 = mct_avect_indexra(x2s,'Fsgg_rofl04')
    index_x2s_Fsgg_hflx04 = mct_avect_indexra(x2s,'Fsgg_hflx04')
    index_x2s_Fsgg_rofi05 = mct_avect_indexra(x2s,'Fsgg_rofi05')
    index_x2s_Fsgg_rofl05 = mct_avect_indexra(x2s,'Fsgg_rofl05')
    index_x2s_Fsgg_hflx05 = mct_avect_indexra(x2s,'Fsgg_hflx05')
#endif
#if (defined GLC_NEC_10 )
    index_x2s_Sg_frac06   = mct_avect_indexra(x2s,'Sg_frac06')
    index_x2s_Sg_topo06   = mct_avect_indexra(x2s,'Sg_topo06')
    index_x2s_Sg_frac07   = mct_avect_indexra(x2s,'Sg_frac07')
    index_x2s_Sg_topo07   = mct_avect_indexra(x2s,'Sg_topo07')
    index_x2s_Sg_frac08   = mct_avect_indexra(x2s,'Sg_frac08')
    index_x2s_Sg_topo08   = mct_avect_indexra(x2s,'Sg_topo08')
    index_x2s_Sg_frac09   = mct_avect_indexra(x2s,'Sg_frac09')
    index_x2s_Sg_topo09   = mct_avect_indexra(x2s,'Sg_topo09')
    index_x2s_Sg_frac10   = mct_avect_indexra(x2s,'Sg_frac10')
    index_x2s_Sg_topo10   = mct_avect_indexra(x2s,'Sg_topo10')
    index_x2s_Fsgg_rofi06 = mct_avect_indexra(x2s,'Fsgg_rofi06')
    index_x2s_Fsgg_rofl06 = mct_avect_indexra(x2s,'Fsgg_rofl06')
    index_x2s_Fsgg_hflx06 = mct_avect_indexra(x2s,'Fsgg_hflx06')
    index_x2s_Fsgg_rofi07 = mct_avect_indexra(x2s,'Fsgg_rofi07')
    index_x2s_Fsgg_rofl07 = mct_avect_indexra(x2s,'Fsgg_rofl07')
    index_x2s_Fsgg_hflx07 = mct_avect_indexra(x2s,'Fsgg_hflx07')
    index_x2s_Fsgg_rofi08 = mct_avect_indexra(x2s,'Fsgg_rofi08')
    index_x2s_Fsgg_rofl08 = mct_avect_indexra(x2s,'Fsgg_rofl08')
    index_x2s_Fsgg_hflx08 = mct_avect_indexra(x2s,'Fsgg_hflx08')
    index_x2s_Fsgg_rofi09 = mct_avect_indexra(x2s,'Fsgg_rofi09')
    index_x2s_Fsgg_rofl09 = mct_avect_indexra(x2s,'Fsgg_rofl09')
    index_x2s_Fsgg_hflx09 = mct_avect_indexra(x2s,'Fsgg_hflx09')
    index_x2s_Fsgg_rofi10 = mct_avect_indexra(x2s,'Fsgg_rofi10')
    index_x2s_Fsgg_rofl10 = mct_avect_indexra(x2s,'Fsgg_rofl10')
    index_x2s_Fsgg_hflx10 = mct_avect_indexra(x2s,'Fsgg_hflx10')
#endif

    nflds_x2s = mct_avect_nRattr(x2s)

    ! sno -> drv

#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 || defined GLC_NEC_1)
    index_s2x_Ss_tsrf01   = mct_avect_indexra(s2x,'Ss_tsrf01')
    index_s2x_Ss_topo01   = mct_avect_indexra(s2x,'Ss_topo01')
    index_s2x_Fgss_qice01 = mct_avect_indexra(s2x,'Fgss_qice01')
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 || defined GLC_NEC_3 )
    index_s2x_Ss_tsrf02   = mct_avect_indexra(s2x,'Ss_tsrf02')
    index_s2x_Ss_topo02   = mct_avect_indexra(s2x,'Ss_topo02')
    index_s2x_Ss_tsrf03   = mct_avect_indexra(s2x,'Ss_tsrf03')
    index_s2x_Ss_topo03   = mct_avect_indexra(s2x,'Ss_topo03')
    index_s2x_Fgss_qice02 = mct_avect_indexra(s2x,'Fgss_qice02')
    index_s2x_Fgss_qice03 = mct_avect_indexra(s2x,'Fgss_qice03')
#endif
#if (defined GLC_NEC_10 || defined GLC_NEC_5 )
    index_s2x_Ss_tsrf04   = mct_avect_indexra(s2x,'Ss_tsrf04')
    index_s2x_Ss_topo04   = mct_avect_indexra(s2x,'Ss_topo04')
    index_s2x_Ss_tsrf05   = mct_avect_indexra(s2x,'Ss_tsrf05')
    index_s2x_Ss_topo05   = mct_avect_indexra(s2x,'Ss_topo05')
    index_s2x_Fgss_qice04 = mct_avect_indexra(s2x,'Fgss_qice04')
    index_s2x_Fgss_qice05 = mct_avect_indexra(s2x,'Fgss_qice05')
#endif
#if (defined GLC_NEC_10 )
    index_s2x_Ss_tsrf06   = mct_avect_indexra(s2x,'Ss_tsrf06')
    index_s2x_Ss_topo06   = mct_avect_indexra(s2x,'Ss_topo06')
    index_s2x_Ss_tsrf07   = mct_avect_indexra(s2x,'Ss_tsrf07')
    index_s2x_Ss_topo07   = mct_avect_indexra(s2x,'Ss_topo07')
    index_s2x_Ss_tsrf08   = mct_avect_indexra(s2x,'Ss_tsrf08')
    index_s2x_Ss_topo08   = mct_avect_indexra(s2x,'Ss_topo08')
    index_s2x_Ss_tsrf09   = mct_avect_indexra(s2x,'Ss_tsrf09')
    index_s2x_Ss_topo09   = mct_avect_indexra(s2x,'Ss_topo09')
    index_s2x_Ss_tsrf10   = mct_avect_indexra(s2x,'Ss_tsrf10')
    index_s2x_Ss_topo10   = mct_avect_indexra(s2x,'Ss_topo10')
    index_s2x_Fgss_qice06 = mct_avect_indexra(s2x,'Fgss_qice06')
    index_s2x_Fgss_qice07 = mct_avect_indexra(s2x,'Fgss_qice07')
    index_s2x_Fgss_qice08 = mct_avect_indexra(s2x,'Fgss_qice08')
    index_s2x_Fgss_qice09 = mct_avect_indexra(s2x,'Fgss_qice09')
    index_s2x_Fgss_qice10 = mct_avect_indexra(s2x,'Fgss_qice10')
#endif

    nflds_s2x = mct_avect_nRattr(s2x)

    call mct_aVect_clean(x2s)
    call mct_aVect_clean(s2x)

  end subroutine clm_cpl_indices_set

!=======================================================================

end module clm_cpl_indices
