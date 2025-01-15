module SnowSnicarMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate albedo of snow containing impurities 
  ! and the evolution of snow effective radius
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use clm_varctl      , only : iulog, snicar_numrad_snw, &
                               snicar_snw_shape, snicar_snobc_intmix, &
                               snicar_snodst_intmix, do_sno_oc
  use clm_varcon      , only : tfrz
  use shr_const_mod   , only : SHR_CONST_RHOICE
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type, subgrid_level_column
  use atm2lndType     , only : atm2lnd_type
  use WaterStateBulkType  , only : waterstatebulk_type
  use WaterDiagnosticBulkType  , only : waterdiagnosticbulk_type
  use WaterFluxBulkType   , only : waterfluxbulk_type
  use TemperatureType , only : temperature_type
  use GridcellType    , only : grc       
  use LandunitType    , only : lun       
  use ColumnType      , only : col       
  !
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SNICAR_RT        ! Snow albedo and vertically-resolved solar absorption
  public :: SnowAge_grain    ! Snow effective grain size evolution
  public :: SnowAge_init     ! Initial read in of snow-aging file
  public :: SnowOptics_init  ! Initial read in of snow-optics file

  type, private :: params_type
      real(r8) :: xdrdt         ! Arbitrary factor applied to snow aging rate (-)
      real(r8) :: snw_rds_refrz ! Effective radius of re-frozen snow (microns)
      real(r8) :: C2_liq_Brun89 ! Constant for liquid water grain growth [m3 s-1],
                                ! from Brun89: corrected for LWC in units of percent
      real(r8) :: fresh_snw_rds_max  ! maximum warm fresh snow effective radius [microns]
      real(r8) :: snw_rds_min  ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
  end type params_type
  type(params_type), private ::  params_inst
  !
  ! !PUBLIC DATA MEMBERS:
  integer,  public, parameter :: sno_nbr_aer =   8        ! number of aerosol species in snowpack
                                                          ! (indices described above) [nbr]
  logical,  public, parameter :: DO_SNO_AER =   .true.    ! parameter to include aerosols in snowpack radiative calculations

  ! !PRIVATE DATA MEMBERS:
  integer, parameter :: default_number_bands = 5  ! currently the only alternative is 480 bands
  integer, parameter :: highest_default_band = 5
  integer, parameter :: sec_highest_default_band = 4
  integer, parameter :: high_number_bands = 480

  integer,  parameter :: idx_Mie_snw_mx = 1471           ! number of effective radius indices used in Mie lookup table [idx]
  integer,  parameter :: idx_T_max      = 11             ! maximum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_T_min      = 1              ! minimum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_max   = 31             ! maximum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_min   = 1              ! minimum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_max   = 8              ! maximum snow density index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_min   = 1              ! minimum snow density index used in aging lookup table [idx]

  integer,  parameter :: snw_rds_max_tbl = 1500          ! maximum effective radius defined in Mie lookup table [microns]
  integer,  parameter :: snw_rds_min_tbl = 30            ! minimium effective radius defined in Mie lookup table [microns]
  real(r8), parameter :: snw_rds_max     = 1500._r8      ! maximum allowed snow effective radius [microns]

  real(r8), parameter :: min_snw = 1.0E-30_r8            ! minimum snow mass required for SNICAR RT calculation [kg m-2]

  real(r8), parameter :: C1_liq_Brun89 = 0._r8           ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: zeroed to accomodate dry snow aging, was 1.28E-17_r8

  real(r8), parameter :: tim_cns_bc_rmv  = 2.2E-8_r8     ! time constant for removal of BC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_oc_rmv  = 2.2E-8_r8     ! time constant for removal of OC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_dst_rmv = 2.2E-8_r8     ! time constant for removal of dust in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)

  ! snow and aerosol Mie parameters:
  ! (arrays declared here, but are set in iniTimeConst)
  ! (idx_Mie_snw_mx is number of snow radii with defined parameters (i.e. from 30um to 1500um))
  
  ! direct-beam weighted ice optical properties
  real(r8), allocatable :: ss_alb_snw_drc(:,:)  ! (idx_Mie_snw_mx, numrad_snw)
  real(r8), allocatable :: asm_prm_snw_drc(:,:)
  real(r8), allocatable :: ext_cff_mss_snw_drc(:,:)

  ! diffuse radiation weighted ice optical properties
  real(r8), allocatable :: ss_alb_snw_dfs(:,:)  ! (idx_Mie_snw_mx, numrad_snw)
  real(r8), allocatable :: asm_prm_snw_dfs(:,:)
  real(r8), allocatable :: ext_cff_mss_snw_dfs(:,:)

  ! hydrophilic BC
  real(r8), allocatable :: ss_alb_bc_hphil(:)  ! (numrad_snw)
  real(r8), allocatable :: asm_prm_bc_hphil(:)
  real(r8), allocatable :: ext_cff_mss_bc_hphil(:)

  ! hydrophobic BC
  real(r8), allocatable :: ss_alb_bc_hphob(:)  ! (numrad_snw)
  real(r8), allocatable :: asm_prm_bc_hphob(:)
  real(r8), allocatable :: ext_cff_mss_bc_hphob(:)

  ! hydrophilic OC
  real(r8), allocatable :: ss_alb_oc_hphil(:)  ! (numrad_snw)
  real(r8), allocatable :: asm_prm_oc_hphil(:)
  real(r8), allocatable :: ext_cff_mss_oc_hphil(:)

  ! hydrophobic OC
  real(r8), allocatable :: ss_alb_oc_hphob(:)  ! (numrad_snw)
  real(r8), allocatable :: asm_prm_oc_hphob(:)
  real(r8), allocatable :: ext_cff_mss_oc_hphob(:)

  ! dust species 1:
  real(r8), allocatable :: ss_alb_dst1(:)  ! (numrad_snw)
  real(r8), allocatable :: asm_prm_dst1(:)
  real(r8), allocatable :: ext_cff_mss_dst1(:)

  ! dust species 2:
  real(r8), allocatable :: ss_alb_dst2(:)  ! (numrad_snw)
  real(r8), allocatable :: asm_prm_dst2(:)
  real(r8), allocatable :: ext_cff_mss_dst2(:)

  ! dust species 3:
  real(r8), allocatable :: ss_alb_dst3(:)  ! (numrad_snw)
  real(r8), allocatable :: asm_prm_dst3(:)
  real(r8), allocatable :: ext_cff_mss_dst3(:)

  ! dust species 4:
  real(r8), allocatable :: ss_alb_dst4(:)  ! (numrad_snw)
  real(r8), allocatable :: asm_prm_dst4(:)
  real(r8), allocatable :: ext_cff_mss_dst4(:)

  ! downward solar radiation spectral weights for 5-band or 480-band
  real(r8), allocatable :: flx_wgt_dir(:)  ! (numrad_snw)  ! direct
  real(r8), allocatable :: flx_wgt_dif(:)  ! (numrad_snw)  ! diffuse

  ! best-fit parameters for snow aging defined over:
  !  11 temperatures from 225 to 273 K
  !  31 temperature gradients from 0 to 300 K/m
  !   8 snow densities from 0 to 350 kg/m3
  ! (arrays declared here, but are set in iniTimeConst)
  real(r8), allocatable :: snowage_tau(:,:,:)  ! (idx_rhos_max, idx_Tgrd_max, idx_T_max)
  real(r8), allocatable :: snowage_kappa(:,:,:)
  real(r8), allocatable :: snowage_drdt0(:,:,:)
  !
  ! !REVISION HISTORY:
  ! Created by Mark Flanner (Univ. of Michigan)
  ! Updated by Cenlin He (NCAR) based on Flanner et al. 2021 GMD

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !----------------------------------------------------------------------------
  subroutine readParams( ncid )
    !
    ! !USES:
    use ncdio_pio, only: file_desc_t
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'readParams_SnowSnicar'
    !--------------------------------------------------------------------

    ! Arbitrary factor applied to snow aging rate (-)
    call readNcdioScalar(ncid, 'xdrdt', subname, params_inst%xdrdt)
    ! Effective radius of re-frozen snow (microns)
    call readNcdioScalar(ncid, 'snw_rds_refrz', subname, params_inst%snw_rds_refrz)
    ! constant for liquid water grain growth [m3 s-1], from Brun89: corrected for LWC in units of percent
    call readNcdioScalar(ncid,  'C2_liq_Brun89', subname, params_inst%C2_liq_Brun89)
    ! maximum warm fresh snow effective radius [microns]
    call readNcdioScalar(ncid, 'fresh_snw_rds_max', subname, params_inst%fresh_snw_rds_max)
    ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
    call readNcdioScalar(ncid, 'snw_rds_min', subname, params_inst%snw_rds_min)

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SNICAR_RT (bounds, num_nourbanc, filter_nourbanc,  &
                        coszen, flg_slr_in, h2osno_liq, h2osno_ice, h2osno_total, snw_rds,   &
                        mss_cnc_aer_in, albsfc, albout, flx_abs, waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Determine reflectance of, and vertically-resolved solar absorption in, 
    ! snow with impurities.
    !
    ! Original references on physical models of snow reflectance include: 
    ! Wiscombe and Warren [1980] and Warren and Wiscombe [1980],
    ! Journal of Atmospheric Sciences, 37,
    !
    ! The multi-layer solution for multiple-scattering used here is from:
    ! Toon et al. [1989], Rapid calculation of radiative heating rates 
    ! and photodissociation rates in inhomogeneous multiple scattering atmospheres, 
    ! J. Geophys. Res., 94, D13, 16287-16301
    !
    ! The implementation of the SNICAR model in CLM/CSIM is described in:
    ! Flanner, M., C. Zender, J. Randerson, and P. Rasch [2007], 
    ! Present-day climate forcing and response from black carbon in snow,
    ! J. Geophys. Res., 112, D11202, doi: 10.1029/2006JD008003
    !
    ! Updated radiative transfer solver:
    !
    ! The multi-layer solution for multiple-scattering used here is from:
    ! Briegleb, P. and Light, B.: A Delta-Eddington mutiple scattering
    ! parameterization for solar radiation in the sea ice component of the
    ! community climate system model, 2007.
    !
    ! The implementation of the SNICAR-AD model in CLM is described in:
    ! Dang et al.2019, Inter-comparison and improvement of 2-stream shortwave
    ! radiative transfer models for unified treatment of cryospheric surfaces
    ! in ESMs; and Flanner et al. 2021, SNICAR-ADv3: a community tool for modeling 
    ! spectral snow albedo
    !
    ! !USES:
    use clm_varpar       , only : nlevsno, numrad, ivis, inir
    use clm_time_manager , only : get_nstep
    use shr_const_mod    , only : SHR_CONST_PI
    !
    ! !ARGUMENTS:
    type (bounds_type), intent(in)  :: bounds                                    
    integer           , intent(in)  :: num_nourbanc                                       ! number of columns in non-urban filter
    integer           , intent(in)  :: filter_nourbanc(:)                                 ! column filter for non-urban points
    real(r8)          , intent(in)  :: coszen         ( bounds%begc: )                    ! cosine of solar zenith angle for next time step (col) [unitless]
    integer           , intent(in)  :: flg_slr_in                                         ! flag: =1 for direct-beam incident flux,=2 for diffuse incident flux
    real(r8)          , intent(in)  :: h2osno_liq     ( bounds%begc: , -nlevsno+1: )      ! liquid water content (col,lyr) [kg/m2]
    real(r8)          , intent(in)  :: h2osno_ice     ( bounds%begc: , -nlevsno+1: )      ! ice content (col,lyr) [kg/m2]
    real(r8)          , intent(in)  :: h2osno_total   ( bounds%begc: )                    ! total snow content (col) [kg/m2] (may differ from sum of h2osno_liq and h2osno_ice if there is snow despite there being no explicit snow layers)
    integer           , intent(in)  :: snw_rds        ( bounds%begc: , -nlevsno+1: )      ! snow effective radius (col,lyr) [microns, m^-6]
    real(r8)          , intent(in)  :: mss_cnc_aer_in ( bounds%begc: , -nlevsno+1: , 1: ) ! mass concentration of all aerosol species (col,lyr,aer) [kg/kg]
    real(r8)          , intent(in)  :: albsfc         ( bounds%begc: , 1: )               ! albedo of surface underlying snow (col,bnd) [frc]
    real(r8)          , intent(out) :: albout         ( bounds%begc: , 1: )               ! snow albedo, averaged into 2 bands (=0 if no sun or no snow) (col,bnd) [frc]
    real(r8)          , intent(out) :: flx_abs        ( bounds%begc: , -nlevsno+1: , 1: ) ! absorbed flux in each layer per unit flux incident (col, lyr, bnd)
    type(waterdiagnosticbulk_type) , intent(in)  :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    !
    ! variables for snow radiative transfer calculations
    integer :: nir_bnd_bgn  ! first band index in near-IR spectrum [idx]
    integer :: nir_bnd_end  ! ending near-IR band index [idx]

    ! Local variables representing single-column values of arrays:
    integer :: snl_lcl                            ! negative number of snow layers [nbr]
    integer :: snw_rds_lcl(-nlevsno+1:0)          ! snow effective radius [m^-6]
    real(r8):: flx_slrd_lcl(1:snicar_numrad_snw)         ! direct beam incident irradiance [W/m2] (set to 1)
    real(r8):: flx_slri_lcl(1:snicar_numrad_snw)         ! diffuse incident irradiance [W/m2] (set to 1)
    real(r8):: mss_cnc_aer_lcl(-nlevsno+1:0,1:sno_nbr_aer) ! aerosol mass concentration (lyr,aer_nbr) [kg/kg]
    real(r8):: h2osno_lcl                         ! total column snow mass [kg/m2]
    real(r8):: h2osno_liq_lcl(-nlevsno+1:0)       ! liquid water mass [kg/m2]
    real(r8):: h2osno_ice_lcl(-nlevsno+1:0)       ! ice mass [kg/m2]
    real(r8):: albsfc_lcl(1:snicar_numrad_snw)           ! albedo of underlying surface [frc]
    real(r8):: ss_alb_snw_lcl(-nlevsno+1:0)       ! single-scatter albedo of ice grains (lyr) [frc]
    real(r8):: asm_prm_snw_lcl(-nlevsno+1:0)      ! asymmetry parameter of ice grains (lyr) [frc]
    real(r8):: ext_cff_mss_snw_lcl(-nlevsno+1:0)  ! mass extinction coefficient of ice grains (lyr) [m2/kg]
    real(r8):: ss_alb_aer_lcl(sno_nbr_aer)        ! single-scatter albedo of aerosol species (aer_nbr) [frc] 
    real(r8):: asm_prm_aer_lcl(sno_nbr_aer)       ! asymmetry parameter of aerosol species (aer_nbr) [frc]
    real(r8):: ext_cff_mss_aer_lcl(sno_nbr_aer)   ! mass extinction coefficient of aerosol species (aer_nbr) [m2/kg]

    ! Other local variables
    integer :: APRX_TYP                           ! two-stream approximation type
                                                  ! (1=Eddington, 2=Quadrature, 3=Hemispheric Mean) [nbr]
    integer :: DELTA                              ! flag to use Delta approximation (Joseph, 1976)
                                                  ! (1= use, 0= don't use)
    real(r8):: flx_wgt(1:snicar_numrad_snw)              ! weights applied to spectral bands,
                                                  ! specific to direct and diffuse cases (bnd) [frc] 
    integer :: flg_nosnl                          ! flag: =1 if there is snow, but zero snow layers,
                                                  ! =0 if at least 1 snow layer [flg]   
    integer :: trip                               ! flag: =1 to redo RT calculation if result is unrealistic
    integer :: flg_dover                          ! defines conditions for RT redo (explained below)
    real(r8):: albedo                             ! temporary snow albedo [frc]
    real(r8):: flx_sum                            ! temporary summation variable for NIR weighting
    real(r8):: albout_lcl(snicar_numrad_snw)             ! snow albedo by band [frc]
    real(r8):: flx_abs_lcl(-nlevsno+1:1,snicar_numrad_snw)! absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]
    real(r8):: L_snw(-nlevsno+1:0)                ! h2o mass (liquid+solid) in snow layer (lyr) [kg/m2]
    real(r8):: tau_snw(-nlevsno+1:0)              ! snow optical depth (lyr) [unitless]
    real(r8):: L_aer(-nlevsno+1:0,sno_nbr_aer)    ! aerosol mass in snow layer (lyr,nbr_aer) [kg/m2] 
    real(r8):: tau_aer(-nlevsno+1:0,sno_nbr_aer)  ! aerosol optical depth (lyr,nbr_aer) [unitless]
    real(r8):: tau_sum                            ! cumulative (snow+aerosol) optical depth [unitless]
    real(r8):: tau_clm(-nlevsno+1:0)              ! column optical depth from layer bottom to snowpack top (lyr) [unitless] 
    real(r8):: omega_sum                          ! temporary summation of single-scatter albedo of all aerosols [frc]
    real(r8):: g_sum                              ! temporary summation of asymmetry parameter of all aerosols [frc]
    real(r8):: tau(-nlevsno+1:0)                  ! weighted optical depth of snow+aerosol layer (lyr) [unitless]
    real(r8):: omega(-nlevsno+1:0)                ! weighted single-scatter albedo of snow+aerosol layer (lyr) [frc]
    real(r8):: g(-nlevsno+1:0)                    ! weighted asymmetry parameter of snow+aerosol layer (lyr) [frc]
    real(r8):: tau_star(-nlevsno+1:0)             ! transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer
                                                  ! (lyr) [unitless]
    real(r8):: omega_star(-nlevsno+1:0)           ! transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer (lyr) [frc]
    real(r8):: g_star(-nlevsno+1:0)               ! transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
                                                  ! (lyr) [frc]
    integer :: nstep                              ! current timestep [nbr] (debugging only)
    integer :: g_idx, c_idx, l_idx                ! gridcell, column, and landunit indices [idx]
    integer :: bnd_idx                            ! spectral band index (1 <= bnd_idx <= snicar_numrad_snw) [idx]
    integer :: rds_idx                            ! snow effective radius index for retrieving
                                                  ! Mie parameters from lookup table [idx]
    integer :: snl_btm                            ! index of bottom snow layer (0) [idx]
    integer :: snl_top                            ! index of top snow layer (-4 to 0) [idx]
    integer :: fc                                 ! column filter index
    integer :: i                                  ! layer index [idx]
    integer :: j                                  ! aerosol number index [idx]
    integer :: n                                  ! tridiagonal matrix index [idx]
    integer :: m                                  ! secondary layer index [idx]   
    real(r8):: F_direct(-nlevsno+1:0)             ! direct-beam radiation at bottom of layer interface (lyr) [W/m^2]
    real(r8):: F_net(-nlevsno+1:0)                ! net radiative flux at bottom of layer interface (lyr) [W/m^2]
    real(r8):: F_abs(-nlevsno+1:0)                ! net absorbed radiative energy (lyr) [W/m^2]
    real(r8):: F_abs_sum                          ! total absorbed energy in column [W/m^2]
    real(r8):: F_sfc_pls                          ! upward radiative flux at snowpack top [W/m^2]
    real(r8):: F_btm_net                          ! net flux at bottom of snowpack [W/m^2]                    
    real(r8):: F_sfc_net                          ! net flux at top of snowpack [W/m^2]
    real(r8):: energy_sum                         ! sum of all energy terms; should be 0.0 [W/m^2]
    real(r8):: F_direct_btm                       ! direct-beam radiation at bottom of snowpack [W/m^2]
    real(r8):: mu_not                             ! cosine of solar zenith angle (used locally) [frc]
    integer :: err_idx                            ! counter for number of times through error loop [nbr]
    real(r8):: lat_coord                          ! gridcell latitude (debugging only)
    real(r8):: lon_coord                          ! gridcell longitude (debugging only)
    integer :: sfctype                            ! underlying surface type (debugging only)
    real(r8):: pi                                 ! 3.1415...

    !-----------------------------------------------------------------------
    ! variables used for Toon et al. 1989 2-stream solver (Flanner et al. 2007):
    ! intermediate variables for radiative transfer approximation:
    real(r8):: gamma1(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(r8):: gamma2(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(r8):: gamma3(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(r8):: gamma4(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(r8):: lambda(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(r8):: GAMMA(-nlevsno+1:0)                ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(r8):: mu_one                             ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real(r8):: e1(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
    real(r8):: e2(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
    real(r8):: e3(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
    real(r8):: e4(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
    real(r8):: C_pls_btm(-nlevsno+1:0)            ! intermediate variable: upward flux at bottom interface (lyr) [W/m2]
    real(r8):: C_mns_btm(-nlevsno+1:0)            ! intermediate variable: downward flux at bottom interface (lyr) [W/m2]
    real(r8):: C_pls_top(-nlevsno+1:0)            ! intermediate variable: upward flux at top interface (lyr) [W/m2]
    real(r8):: C_mns_top(-nlevsno+1:0)            ! intermediate variable: downward flux at top interface (lyr) [W/m2]
    real(r8):: A(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(r8):: B(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(r8):: D(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(r8):: E(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(r8):: AS(-2*nlevsno+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(r8):: DS(-2*nlevsno+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(r8):: X(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real(r8):: Y(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)

    !-----------------------------------------------------------------------
    ! variables used for Adding-doubling 2-stream solver based on SNICAR-ADv3 version 
    ! (Dang et al. 2019; Flanner et al. 2021)
    real(r8):: trndir(-nlevsno+1:1)               ! solar beam down transmission from top
    real(r8):: trntdr(-nlevsno+1:1)               ! total transmission to direct beam for layers above
    real(r8):: trndif(-nlevsno+1:1)               ! diffuse transmission to diffuse beam for layers above
    real(r8):: rupdir(-nlevsno+1:1)               ! reflectivity to direct radiation for layers below
    real(r8):: rupdif(-nlevsno+1:1)               ! reflectivity to diffuse radiation for layers below
    real(r8):: rdndif(-nlevsno+1:1)               ! reflectivity to diffuse radiation for layers above
    real(r8):: dfdir(-nlevsno+1:1)                ! down-up flux at interface due to direct beam at top surface
    real(r8):: dfdif(-nlevsno+1:1)                ! down-up flux at interface due to diffuse beam at top surface
    real(r8):: dftmp(-nlevsno+1:1)                ! temporary variable for down-up flux at interface
    real(r8):: rdir(-nlevsno+1:0)                 ! layer reflectivity to direct radiation
    real(r8):: rdif_a(-nlevsno+1:0)               ! layer reflectivity to diffuse radiation from above
    real(r8):: rdif_b(-nlevsno+1:0)               ! layer reflectivity to diffuse radiation from below
    real(r8):: tdir(-nlevsno+1:0)                 ! layer transmission to direct radiation (solar beam + diffuse)
    real(r8):: tdif_a(-nlevsno+1:0)               ! layer transmission to diffuse radiation from above
    real(r8):: tdif_b(-nlevsno+1:0)               ! layer transmission to diffuse radiation from below
    real(r8):: trnlay(-nlevsno+1:0)               ! solar beam transm for layer (direct beam only)
    real(r8):: ts                                 ! layer delta-scaled extinction optical depth
    real(r8):: ws                                 ! layer delta-scaled single scattering albedo
    real(r8):: gs                                 ! layer delta-scaled asymmetry parameter
    real(r8):: extins                             ! extinction
    real(r8):: alp                                ! temporary for alpha
    real(r8):: gam                                ! temporary for agamm
    real(r8):: amg                                ! alp - gam
    real(r8):: apg                                ! alp + gam
    real(r8):: ue                                 ! temporary for u
    real(r8):: refk                               ! interface multiple scattering
    real(r8):: refkp1                             ! interface multiple scattering for k+1
    real(r8):: refkm1                             ! interface multiple scattering for k-1
    real(r8):: tdrrdir                            ! direct tran times layer direct ref
    real(r8):: tdndif                             ! total down diffuse = tot tran - direct tran
    real(r8):: taus                               ! scaled extinction optical depth
    real(r8):: omgs                               ! scaled single particle scattering albedo
    real(r8):: asys                               ! scaled asymmetry parameter
    real(r8):: lm                                 ! temporary for el
    real(r8):: mu                                 ! cosine solar zenith for either snow or water
    real(r8):: ne                                 ! temporary for n
    real(r8):: R1                                 ! perpendicular polarization reflection amplitude
    real(r8):: R2                                 ! parallel polarization reflection amplitude
    real(r8):: T1                                 ! perpendicular polarization transmission amplitude
    real(r8):: T2                                 ! parallel polarization transmission amplitude
    real(r8):: Rf_dir_a                           ! fresnel reflection to direct radiation
    real(r8):: Tf_dir_a                           ! fresnel transmission to direct radiation
    real(r8):: Rf_dif_a                           ! fresnel reflection to diff radiation from above
    real(r8):: Rf_dif_b                           ! fresnel reflection to diff radiation from below
    real(r8):: Tf_dif_a                           ! fresnel transmission to diff radiation from above
    real(r8):: Tf_dif_b                           ! fresnel transmission to diff radiation from below
    real(r8):: gwt                                ! gaussian weight
    real(r8):: swt                                ! sum of weights
    real(r8):: trn                                ! layer transmission
    real(r8):: rdr                                ! rdir for gaussian integration
    real(r8):: tdr                                ! tdir for gaussian integration
    real(r8):: smr                                ! accumulator for rdif gaussian integration
    real(r8):: smt                                ! accumulator for tdif gaussian integration
    real(r8):: exp_min                            ! minimum exponential value

    integer :: ng  ! gaussian integration index
    integer, parameter :: ngmax = 8  ! max gaussian integration index
    real(r8), parameter :: difgauspt(ngmax) = &  ! Gaussian integration angles (radians)
                      (/ 0.9894009_r8,  0.9445750_r8, &
                         0.8656312_r8,  0.7554044_r8, &
                         0.6178762_r8,  0.4580168_r8, &
                         0.2816036_r8,  0.0950125_r8/)
    real(r8), parameter :: difgauswt(ngmax) = &  ! Gaussian integration coefficients/weights
                      (/ 0.0271525_r8,  0.0622535_r8, &
                         0.0951585_r8,  0.1246290_r8, &
                         0.1495960_r8,  0.1691565_r8, &
                         0.1826034_r8,  0.1894506_r8/)

    integer :: snl_btm_itf                        ! index of bottom snow layer interfaces (1) [idx]
    ! constants used in algorithm
    real(r8), parameter :: c0 = 0.0_r8
    real(r8), parameter :: c1 = 1.0_r8
    real(r8), parameter :: c3 = 3.0_r8
    real(r8), parameter :: c4 = 4.0_r8
    real(r8), parameter :: c6 = 6.0_r8
    real(r8), parameter :: cp01 = 0.01_r8
    real(r8), parameter :: cp5 = 0.5_r8
    real(r8), parameter :: cp75 = 0.75_r8
    real(r8), parameter :: c1p5 = 1.5_r8
    real(r8), parameter :: trmin = 0.001_r8
    real(r8), parameter :: argmax = 10.0_r8  ! maximum argument of exponential
    ! constants and coefficients used for SZA parameterization
    real(r8), parameter :: sza_a0 =  0.085730_r8
    real(r8), parameter :: sza_a1 = -0.630883_r8
    real(r8), parameter :: sza_a2 =  1.303723_r8
    real(r8), parameter :: sza_b0 =  1.467291_r8
    real(r8), parameter :: sza_b1 = -3.338043_r8
    real(r8), parameter :: sza_b2 =  6.807489_r8
    real(r8), parameter :: puny =  1.0e-11_r8
    real(r8), parameter :: mu_75 =  0.2588_r8  ! cosine of 75 degree
    real(r8):: sza_c1                             ! coefficient, SZA parameteirzation
    real(r8):: sza_c0                             ! coefficient, SZA parameterization
    real(r8):: sza_factor                         ! factor used to adjust NIR direct albedo
    real(r8):: flx_sza_adjust                     ! direct NIR flux adjustment from sza_factor
    real(r8):: mu0                                ! incident solar zenith angle

    !-----------------------------------------------------------------------
    ! variables used for nonspherical snow grain treatment (He et al. 2017 J of Climate):
    character(len=15) :: sno_shp(-nlevsno+1:0)    ! Snow shape type: sphere, spheroid, hexagonal plate, koch snowflake
                                                  ! currently only assuming same shapes for all snow layers
    real(r8) :: sno_fs(-nlevsno+1:0)              ! Snow shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                                  ! only activated when snicar_snw_shape is nonspherical
                                                  ! 0=use recommended default value (He et al. 2017);
                                                  ! others(0<sno_fs<1)= user-specified value
    real(r8) :: sno_AR(-nlevsno+1:0)              ! Snow grain aspect ratio: ratio of grain width to length
                                                  ! only activated when snicar_snw_shape is nonspherical
                                                  ! 0=use recommended default value (He et al. 2017);
                                                  ! others(0.1<fs<20)= use user-specified value
    ! Constants and parameters for aspherical ice particles    
    ! asymmetry factor parameterization coefficients (6 bands) from Table 3 & Eqs. 6-7 in He et al. (2017)
    integer, parameter :: seven_bands = 7
    real(r8) :: g_wvl_ct(seven_bands)  ! center point for wavelength band (um)
    real(r8), parameter :: g_wvl(seven_bands+1) =  &  ! wavelength (um) division point
        (/ 0.25_r8, 0.70_r8, 1.41_r8, 1.90_r8, &
           2.50_r8, 3.50_r8, 4.00_r8, 5.00_r8 /)
    real(r8), parameter :: g_b0(seven_bands) =  &
        (/  9.76029E-1_r8,  9.67798E-1_r8,  1.00111_r8, 1.00224_r8, &
            9.64295E-1_r8,  9.97475E-1_r8,  9.97475E-1_r8 /)
    real(r8), parameter :: g_b1(seven_bands) =  &
        (/  5.21042E-1_r8,  4.96181E-1_r8,  1.83711E-1_r8,  1.37082E-1_r8, &
            5.50598E-2_r8,  8.48743E-2_r8,  8.48743E-2_r8 /)
    real(r8), parameter :: g_b2(seven_bands) =  &
        (/ -2.66792E-4_r8,  1.14088E-3_r8,  2.37011E-4_r8, -2.35905E-4_r8, &
            8.40449E-4_r8, -4.71484E-4_r8, -4.71484E-4_r8 /)
    ! Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007 JAS
    real(r8), parameter :: g_F07_c2(seven_bands) =  &
        (/  1.349959E-1_r8,  1.115697E-1_r8,  9.853958E-2_r8,  5.557793E-2_r8, &
            -1.233493E-1_r8,  0.0_r8        ,  0.0_r8         /)
    real(r8), parameter :: g_F07_c1(seven_bands) =  &
        (/ -3.987320E-1_r8, -3.723287E-1_r8, -3.924784E-1_r8, -3.259404E-1_r8, &
            4.429054E-2_r8, -1.726586E-1_r8, -1.726586E-1_r8 /)
    real(r8), parameter :: g_F07_c0(seven_bands) =  &
        (/  7.938904E-1_r8,  8.030084E-1_r8,  8.513932E-1_r8,  8.692241E-1_r8, &
            7.085850E-1_r8,  6.412701E-1_r8,  6.412701E-1_r8 /)
    real(r8), parameter :: g_F07_p2(seven_bands) =  &
        (/  3.165543E-3_r8,  2.014810E-3_r8,  1.780838E-3_r8,  6.987734E-4_r8, &
            -1.882932E-2_r8, -2.277872E-2_r8, -2.277872E-2_r8 /)
    real(r8), parameter :: g_F07_p1(seven_bands) =  &
        (/  1.140557E-1_r8,  1.143152E-1_r8,  1.143814E-1_r8,  1.071238E-1_r8, &
            1.353873E-1_r8,  1.914431E-1_r8,  1.914431E-1_r8 /)
    real(r8), parameter :: g_F07_p0(seven_bands) =  &
        (/  5.292852E-1_r8,  5.425909E-1_r8,  5.601598E-1_r8,  6.023407E-1_r8, &
            6.473899E-1_r8,  4.634944E-1_r8,  4.634944E-1_r8 /)

    ! other temporary variables
    real(r8), allocatable :: wvl_ct(:)            ! band center wavelength (um) for 5 or 480-band case
    real(r8) :: diam_ice                          ! effective snow grain diameter (SSA-equivalent) unit: microns
    real(r8) :: fs_sphd                           ! shape factor for spheroid snow
    real(r8), parameter :: fs_sphd_default = 0.929_r8  ! default; He et al. (2017), Table 1
    real(r8) :: fs_hex                            ! shape factor for reference hexagonal snow
    real(r8) :: fs_hex0                           ! shape factor for hexagonal plate
    real(r8), parameter :: fs_hex_ref = 0.788_r8  ! reference shape factor
    real(r8) :: fs_koch                           ! shape factor for Koch snowflake
    real(r8), parameter :: fs_koch_default = 0.712_r8  ! default; He et al. (2017), Table 1
    real(r8) :: AR_tmp                            ! aspect ratio temporary
    real(r8), parameter :: AR_tmp_default_1 = 0.5_r8  ! default; He et al. (2017), Table 1
    real(r8), parameter :: AR_tmp_default_2 = 2.5_r8  ! default; He et al. (2017), Table 1
    real(r8) :: g_ice_Cg_tmp(seven_bands)  ! temporary asymmetry factor correction coeff
    real(r8) :: gg_ice_F07_tmp(seven_bands)  ! temporary asymmetry factor related to geometric reflection & refraction
    real(r8) :: g_Cg_intp                         ! interpolated asymmetry factor correction coeff to target bands 
    real(r8) :: gg_F07_intp                       ! interpolated asymmetry factor related to geometric reflection & refraction
    real(r8) :: g_ice_F07                         ! asymmetry factor for Fu 2007 parameterization value
    integer  :: igb                               ! loop index

    !-----------------------------------------------------------------------
    ! variables used for BC-snow internal mixing (He et al. 2017 J of Climate):
    real(r8) :: enh_omg_bcint                     ! BC-induced enhancement in snow single-scattering co-albedo (1-omega)
    integer, parameter :: sixteen_bands = 16
    real(r8) :: enh_omg_bcint_tmp(1:sixteen_bands)  ! temporary BC-induced enhancement in snow 1-omega
    real(r8) :: enh_omg_bcint_tmp2(1:sixteen_bands)  ! temporary BC-induced enhancement in snow 1-omega
    real(r8) :: bcint_wvl_ct(1:sixteen_bands)  ! Parameterization band center wavelength (um)
    ! initialize for BC-snow internal mixing
    ! Eq. 8b & Table 4 in He et al., 2017 J. Climate (wavelength>1.2um, no BC-snow int mixing effect)
    real(r8), parameter :: bcint_wvl(sixteen_bands+1) =  &  ! Parameterization band (0.2-1.2um) for BC-induced enhancement in snow 1-omega
        (/ 0.20_r8, 0.25_r8, 0.30_r8, 0.33_r8, 0.36_r8, 0.40_r8, 0.44_r8, 0.48_r8, &
           0.52_r8, 0.57_r8, 0.64_r8, 0.69_r8, 0.75_r8, 0.78_r8, 0.87_r8, 1._r8, 1.2_r8 /)
    real(r8), parameter :: bcint_d0(sixteen_bands) =  &  ! Parameterization coefficients at each band center wavelength
        (/ 2.48045_r8   , 4.70305_r8   , 4.68619_r8   , 4.67369_r8   , 4.65040_r8   , &
           2.40364_r8   , 7.95408E-1_r8, 2.92745E-1_r8, 8.63396E-2_r8, 2.76299E-2_r8, &
           1.40864E-2_r8, 8.65705E-3_r8, 6.12971E-3_r8, 4.45697E-3_r8, 3.06648E-2_r8, &
           7.96544E-1_r8 /)
    real(r8), parameter :: bcint_d1(sixteen_bands) =  &  ! Parameterization coefficients at each band center wavelength
        (/ 9.77209E-1_r8, 9.73317E-1_r8, 9.79650E-1_r8, 9.84579E-1_r8, 9.93537E-1_r8, &
           9.95955E-1_r8, 9.95218E-1_r8, 9.74284E-1_r8, 9.81193E-1_r8, 9.81239E-1_r8, &
           9.55515E-1_r8, 9.10491E-1_r8, 8.74196E-1_r8, 8.27238E-1_r8, 4.82870E-1_r8, &
           4.36649E-2_r8 /)
    real(r8), parameter :: bcint_d2(sixteen_bands) =  &  ! Parameterization coefficients at each band center wavelength
        (/ 3.95960E-1_r8, 2.04820E-1_r8, 2.07410E-1_r8, 2.09390E-1_r8, 2.13030E-1_r8, &
           4.18570E-1_r8, 1.29682_r8   , 3.75514_r8   , 1.27372E+1_r8, 3.93293E+1_r8, &
           8.78918E+1_r8, 1.86969E+2_r8, 3.45600E+2_r8, 7.08637E+2_r8, 1.41067E+3_r8, &
           2.57288E+2_r8 /)
    real(r8), parameter :: den_bc = 1.7_r8  ! BC particle density (g/cm3)
    real(r8), parameter :: den_bc_target = 1.49_r8  ! target BC particle density (g/cm3) used in BC MAC adjustment
    real(r8), parameter :: Re_bc = 0.045_r8  ! target BC effective radius (um) used in BC MAC adjustment
    real(r8), parameter :: radius_1 = 0.1_r8  ! used with Re_bc (um)
    real(r8), parameter :: radius_2 = 0.05_r8  ! used with Re_bc (um)
    ! Eq. 1a,1b and Table S1 in He et al. 2018 GRL
    ! Parameterization coefficients for BC size adjustment in BC-snow int mix
    integer, parameter :: three_bands = 3
    real(r8), parameter :: bcint_m(three_bands) = (/ -0.8724_r8, -0.1866_r8, -0.0046_r8 /)
    real(r8), parameter :: bcint_n(three_bands) = (/ -0.0072_r8, -0.1918_r8, -0.5177_r8 /)

    real(r8) :: bcint_m_tmp                       ! temporary of bcint_m
    real(r8) :: bcint_n_tmp                       ! temporary of bcint_n
    real(r8) :: bcint_dd                          ! intermediate parameter
    real(r8) :: bcint_dd2                         ! intermediate parameter
    real(r8) :: bcint_f                           ! intermediate parameter
    real(r8) :: enh_omg_bcint_intp                ! BC-induced enhancement in snow 1-omega (logscale) interpolated to CLM wavelength
    real(r8) :: enh_omg_bcint_intp2               ! BC-induced enhancement in snow 1-omega interpolated to CLM wavelength
    real(r8) :: wvl_doint                         ! wavelength doing BC-snow int mixing (<=1.2um)
    integer  :: ibb                               ! loop index

    !-----------------------------------------------------------------------
    ! variables used for dust-snow internal mixing (He et al. 2019 JAMES):
    real(r8) :: enh_omg_dstint                    ! dust-induced enhancement in snow single-scattering co-albedo (1-omega)
    integer, parameter :: size_bins = 6
    real(r8) :: enh_omg_dstint_tmp(size_bins)  ! temporary dust-induced enhancement in snow 1-omega
    real(r8) :: enh_omg_dstint_tmp2(size_bins)  ! temporary dust-induced enhancement in snow 1-omega
    real(r8) :: dstint_wvl_ct(size_bins)  ! Parameterization band center wavelength (um)
    ! initialize for dust-snow internal mixing
    ! Eq. 1 and Table 1 in He et al. 2019 JAMES (wavelength>1.2um, no dust-snow int mixing effect)
    real(r8), parameter :: dstint_wvl(size_bins+1) =  &  ! Parameterization band (0.2-1.2um) for dust-induced enhancement in snow 1-omega
        (/ 0.2_r8, 0.2632_r8, 0.3448_r8, 0.4415_r8, 0.625_r8, 0.7782_r8, 1.2422_r8/)
    real(r8), parameter :: dstint_a1(size_bins) =  &  ! Parameterization coefficients at each band center wavelength
        (/ -2.1307E+1_r8, -1.5815E+1_r8, -9.2880_r8   , 1.1115_r8   , 1.0307_r8   , 1.0185_r8    /)
    real(r8), parameter :: dstint_a2(size_bins) =  &  ! Parameterization coefficients at each band center wavelength
        (/  1.1746E+2_r8,  9.3241E+1_r8,  4.0605E+1_r8, 3.7389E-1_r8, 1.4800E-2_r8, 2.8921E-4_r8 /)
    real(r8), parameter :: dstint_a3(size_bins) =  &  ! Parameterization coefficients at each band center wavelength
        (/  9.9701E-1_r8,  9.9781E-1_r8,  9.9848E-1_r8, 1.0035_r8   , 1.0024_r8   , 1.0356_r8    /)

    real(r8) :: enh_omg_dstint_intp               ! dust-induced enhancement in snow 1-omega (logscale) interpolated to CLM wavelength
    real(r8) :: enh_omg_dstint_intp2              ! dust-induced enhancement in snow 1-omega interpolated to CLM wavelength
    real(r8) :: tot_dst_snw_conc                  ! total dust content in snow across all size bins (ppm=ug/g)
    integer  :: idb                               ! loop index

    real(r8), parameter :: enh_omg_max = 1.e5_r8  ! reasonable maximum value for enh_omg_[bc,dst]int_intp2

    ! unit conversions
    real(r8), parameter :: kg_kg_to_ppm = 1.e6_r8  ! kg/kg to ppm
    real(r8), parameter :: kg_to_ug = 1.e9_r8  ! kg to micrograms

    character(len=*), parameter :: subname = 'SNICAR_RT'
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(coszen)         == (/bounds%endc/)),                 sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(h2osno_liq)     == (/bounds%endc, 0/)),              sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(h2osno_ice)     == (/bounds%endc, 0/)),              sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(h2osno_total)   == (/bounds%endc/)),                 sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(snw_rds)        == (/bounds%endc, 0/)),              sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(mss_cnc_aer_in) == (/bounds%endc, 0, sno_nbr_aer/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(albsfc)         == (/bounds%endc, numrad/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(albout)         == (/bounds%endc, numrad/)),         sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(flx_abs)        == (/bounds%endc, 1, numrad/)),      sourcefile, __LINE__)

    associate(& 
         snl         =>   col%snl                           , & ! Input:  [integer (:)]  negative number of snow layers (col) [nbr]
         frac_sno    =>   waterdiagnosticbulk_inst%frac_sno_eff_col    & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
         )

      ! initialize parameter and
      ! SNICAR/CLM snow band center wavelength (um)
      allocate(wvl_ct(snicar_numrad_snw))
      select case (snicar_numrad_snw)
      case (default_number_bands)
         nir_bnd_bgn = 2
         wvl_ct(:)  = (/ 0.5_r8, 0.85_r8, 1.1_r8, 1.35_r8, 3.25_r8 /)  ! 5-band
      case (high_number_bands)
         nir_bnd_bgn = 51
         do igb = 1, snicar_numrad_snw
            wvl_ct(igb) = 0.205_r8 + 0.01_r8 * (igb - 1._r8)  ! 480-band
         enddo
      case default
         write(iulog,*) subname//' ERROR: unknown snicar_numrad_snw value: ', snicar_numrad_snw
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end select
      nir_bnd_end    = snicar_numrad_snw 

      ! initialize for nonspherical snow grains
      sno_shp(:) = snicar_snw_shape ! currently only assuming same shapes for all snow layers
      sno_fs(:)  = 0._r8
      sno_AR(:)  = 0._r8

      g_wvl_ct(1:seven_bands) = g_wvl(2:seven_bands+1) * 0.5_r8 + g_wvl(1:seven_bands) * 0.5_r8
      dstint_wvl_ct(1:size_bins) = dstint_wvl(2:size_bins+1) * 0.5_r8 + dstint_wvl(1:size_bins) * 0.5_r8
      bcint_wvl_ct(1:sixteen_bands) = bcint_wvl(2:sixteen_bands+1) * 0.5_r8 + bcint_wvl(1:sixteen_bands) * 0.5_r8

      ! Define constants
      pi = SHR_CONST_PI

      ! always use Delta approximation for snow
      DELTA = 1

      ! Get current timestep
      nstep = get_nstep()

      ! Loop over all non-urban columns
      do fc = 1,num_nourbanc
         c_idx = filter_nourbanc(fc)

         ! Zero absorbed radiative fluxes:
         do i=-nlevsno+1,1,1
            flx_abs_lcl(i,:)   = 0._r8
            flx_abs(c_idx,i,:) = 0._r8
         enddo

         ! set snow/ice mass to be used for RT:
         h2osno_lcl = h2osno_total(c_idx)

         ! Qualifier for computing snow RT: 
         !  1) sunlight from atmosphere model 
         !  2) minimum amount of snow on ground. 
         !     Otherwise, set snow albedo to zero
         if ((coszen(c_idx) > 0._r8) .and. (h2osno_lcl > min_snw)) then     

            ! Set variables specific to CLM
            ! If there is snow, but zero snow layers, we must create a layer locally.
            ! This layer is presumed to have the fresh snow effective radius.
            if (snl(c_idx) > -1) then
               flg_nosnl         =  1
               snl_lcl           =  -1
               h2osno_ice_lcl(0) =  h2osno_lcl
               h2osno_liq_lcl(0) =  0._r8
               snw_rds_lcl(0)    =  nint(params_inst%snw_rds_min)
            else
               flg_nosnl         =  0
               snl_lcl           =  snl(c_idx)
               h2osno_liq_lcl(:) =  h2osno_liq(c_idx,:)
               h2osno_ice_lcl(:) =  h2osno_ice(c_idx,:)
               snw_rds_lcl(:)    =  snw_rds(c_idx,:)
            endif

            snl_btm   = 0
            snl_top   = snl_lcl+1

            ! for debugging only
            l_idx     = col%landunit(c_idx)
            g_idx     = col%gridcell(c_idx)
            sfctype   = lun%itype(l_idx)
            lat_coord = grc%latdeg(g_idx)
            lon_coord = grc%londeg(g_idx)


            ! Set local aerosol array
            do j=1,sno_nbr_aer
               mss_cnc_aer_lcl(:,j) = mss_cnc_aer_in(c_idx,:,j)
            enddo


            ! Set spectral underlying surface albedos to their corresponding VIS or NIR albedos
            albsfc_lcl(1:(nir_bnd_bgn-1))       = albsfc(c_idx,ivis)
            albsfc_lcl(nir_bnd_bgn:nir_bnd_end) = albsfc(c_idx,inir)
            

            ! Error check for snow grain size:
            do i=snl_top,snl_btm,1
               if ((snw_rds_lcl(i) < snw_rds_min_tbl) .or. (snw_rds_lcl(i) > snw_rds_max_tbl)) then
                  write (iulog,*) "SNICAR ERROR: snow grain radius of ", snw_rds_lcl(i), " out of bounds."
                  write (iulog,*) "NSTEP= ", nstep
                  write (iulog,*) "column: ", c_idx, " level: ", i, " snl(c)= ", snl_lcl
                  write (iulog,*) "lat= ", lat_coord, " lon= ", lon_coord
                  write (iulog,*) "h2osno_total(c)= ", h2osno_lcl
                  call endrun(subgrid_index=c_idx, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
               endif
            enddo


            ! Incident flux weighting parameters
            !  - sum of all VIS bands must equal 1
            !  - sum of all NIR bands must equal 1
            !
            ! Spectral bands (5-band case)
            !  Band 1: 0.3-0.7um (VIS)
            !  Band 2: 0.7-1.0um (NIR)
            !  Band 3: 1.0-1.2um (NIR)
            !  Band 4: 1.2-1.5um (NIR)
            !  Band 5: 1.5-5.0um (NIR)
            !
            ! Hyperspectral (10-nm) bands (480-band case)
            ! Bands 1~50  : 0.2-0.7um (VIS)
            ! Bands 51~480: 0.7~5.0um (NIR)
            !
            ! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
            !
!           ! works for both 5-band & 480-band, flux weights directly read from input data
            ! Direct:
            if (flg_slr_in == 1) then
               flx_wgt(1:snicar_numrad_snw) = flx_wgt_dir(1:snicar_numrad_snw)  ! VIS or NIR band sum is already normalized to 1.0 in input data
            ! Diffuse:
            elseif (flg_slr_in == 2) then
               flx_wgt(1:snicar_numrad_snw) = flx_wgt_dif(1:snicar_numrad_snw)  ! VIS or NIR band sum is already normalized to 1.0 in input data
            endif

            exp_min = exp(-argmax)

            ! Loop over snow spectral bands
            do bnd_idx = 1,snicar_numrad_snw
               ! flg_dover is not used since this algorithm is stable for mu_not > 0.01
               ! mu_not is cosine solar zenith angle above the fresnel level; make
               ! sure mu_not is large enough for stable and meaningful radiation
               ! solution: .01 is like sun just touching horizon with its lower edge
               ! equivalent to mu0 in sea-ice shortwave model ice_shortwave.F90
               mu_not = max(coszen(c_idx), cp01)

               flg_dover = 1    ! default is to redo
               err_idx   = 0    ! number of times through loop

               do while (flg_dover > 0)

                  ! Set direct or diffuse incident irradiance to 1
                  ! (This has to be within the bnd loop because mu_not is adjusted in rare cases)
                  if (flg_slr_in == 1) then
                     flx_slrd_lcl(bnd_idx) = 1._r8/(mu_not*pi) ! this corresponds to incident irradiance of 1.0
                     flx_slri_lcl(bnd_idx) = 0._r8
                  else
                     flx_slrd_lcl(bnd_idx) = 0._r8
                     flx_slri_lcl(bnd_idx) = 1._r8
                  endif

                  ! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
                  ! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
                  if (snicar_numrad_snw == default_number_bands .and. (bnd_idx == highest_default_band .or. bnd_idx == sec_highest_default_band)) then
                     mss_cnc_aer_lcl(:,:) = 0._r8
                  endif

                  if ( (snicar_numrad_snw == high_number_bands).and.(bnd_idx > 100) ) then ! >1.2um
                     mss_cnc_aer_lcl(:,:) = 0._r8
                  endif


                  !--------------------------- Start snow & aerosol optics --------------------------------
                  ! Define local Mie parameters based on snow grain size and aerosol species retrieved from a lookup table.

                  ! Spherical snow: single-scatter albedo, mass extinction coefficient, asymmetry factor
                  if (flg_slr_in == 1) then
                     do i=snl_top,snl_btm,1
                        rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                        ! snow optical properties (direct radiation)
                        ss_alb_snw_lcl(i)      = ss_alb_snw_drc(rds_idx,bnd_idx)
                        ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_drc(rds_idx,bnd_idx)
                        if (sno_shp(i) == 'sphere') asm_prm_snw_lcl(i) = asm_prm_snw_drc(rds_idx,bnd_idx)
                     enddo
                  elseif (flg_slr_in == 2) then
                     do i=snl_top,snl_btm,1
                        rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                        ! snow optical properties (diffuse radiation)
                        ss_alb_snw_lcl(i)      = ss_alb_snw_dfs(rds_idx,bnd_idx)
                        ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_dfs(rds_idx,bnd_idx)
                        if (sno_shp(i) == 'sphere') asm_prm_snw_lcl(i) = asm_prm_snw_dfs(rds_idx,bnd_idx)
                     enddo
                  endif

                  ! Nonspherical snow: shape-dependent asymmetry factors
                  do i=snl_top,snl_btm,1

                     select case (sno_shp(i))
                     case ('spheroid')
                        diam_ice = 2._r8 * snw_rds_lcl(i)   ! unit: microns
                        if (sno_fs(i) == 0._r8) then
                           fs_sphd = fs_sphd_default  ! default; He et al. (2017), Table 1
                        else
                           fs_sphd = sno_fs(i) ! user specified value
                        endif
                        fs_hex = fs_hex_ref  ! reference shape factor
                        if (sno_AR(i) == 0._r8) then
                           AR_tmp = AR_tmp_default_1  ! default; He et al. (2017), Table 1
                        else
                           AR_tmp = sno_AR(i)  ! user specified value
                        endif
                        do igb = 1, seven_bands
                           g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_sphd/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
                           gg_ice_F07_tmp(igb) = g_F07_c0(igb) + g_F07_c1(igb) * AR_tmp + g_F07_c2(igb) * (AR_tmp * AR_tmp)  ! Eqn. 3.1 in Fu (2007)
                        enddo

                     case ('hexagonal_plate')
                        diam_ice = 2._r8 * snw_rds_lcl(i)   ! unit: microns
                        if (sno_fs(i) == 0._r8) then
                           fs_hex0 = fs_hex_ref  ! default; He et al. (2017), Table 1
                        else
                           fs_hex0 = sno_fs(i) ! user specified value
                        endif
                        fs_hex = fs_hex_ref  ! reference shape factor
                        if (sno_AR(i) == 0._r8) then
                           AR_tmp = AR_tmp_default_2  ! default; He et al. (2017), Table 1
                        else
                           AR_tmp = sno_AR(i)  ! user specified value
                        endif
                        do igb = 1, seven_bands
                           g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_hex0/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
                           gg_ice_F07_tmp(igb) = g_F07_p0(igb) + g_F07_p1(igb) * log(AR_tmp) + g_F07_p2(igb) * (log(AR_tmp) * log(AR_tmp)) ! Eqn. 3.3 in Fu (2007)
                        enddo

                     case ('koch_snowflake')
                        diam_ice = 2._r8 * snw_rds_lcl(i) / 0.544_r8  ! unit: microns
                        if (sno_fs(i) == 0._r8) then
                           fs_koch = fs_koch_default  ! default; He et al. (2017), Table 1
                        else
                           fs_koch = sno_fs(i) ! user specified value
                        endif
                        fs_hex = fs_hex_ref  ! reference shape factor
                        if (sno_AR(i) == 0._r8) then
                           AR_tmp = AR_tmp_default_2  ! default; He et al. (2017), Table 1
                        else
                           AR_tmp = sno_AR(i)  ! user specified value
                        endif
                        do igb = 1, seven_bands
                           g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_koch/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
                           gg_ice_F07_tmp(igb) = g_F07_p0(igb) + g_F07_p1(igb) * log(AR_tmp) + g_F07_p2(igb) * (log(AR_tmp) * log(AR_tmp)) ! Eqn. 3.3 in Fu (2007)
                        enddo

                     case ('sphere')
                        ! DO NOTHING
                     case default
                        write(iulog,*) subname//' ERROR: unknown sno_shp for i: ', sno_shp(i), i
                        call endrun(msg=errMsg(sourcefile, __LINE__))
                     end select

                     ! compute nonspherical snow asymmetry factor
                     if (sno_shp(i) /= 'sphere') then
                        ! 7 wavelength bands for g_ice to be interpolated into targeted SNICAR bands here
                        ! use the piecewise linear interpolation subroutine created at the end of this module
                        ! tests showed the piecewise linear interpolation has similar results as pchip interpolation
                        call piecewise_linear_interp1d(seven_bands, g_wvl_ct, g_ice_Cg_tmp, wvl_ct(bnd_idx), g_Cg_intp)
                        call piecewise_linear_interp1d(seven_bands, g_wvl_ct, gg_ice_F07_tmp, wvl_ct(bnd_idx), gg_F07_intp)
                        g_ice_F07 = gg_F07_intp + 0.5_r8 * (1._r8 - gg_F07_intp) / ss_alb_snw_lcl(i)  ! Eq.2.2 in Fu (2007)
                        asm_prm_snw_lcl(i) = g_ice_F07 * g_Cg_intp     ! Eq.6, He et al. (2017)
                     endif

                     asm_prm_snw_lcl(i) = min(0.99_r8, asm_prm_snw_lcl(i))  !avoid unreasonable values (rarely occur in large-size spheroid cases)

                  enddo ! snow layer loop

                  ! aerosol species 2 optical properties, hydrophobic BC
                  ss_alb_aer_lcl(2)      = ss_alb_bc_hphob(bnd_idx)
                  asm_prm_aer_lcl(2)     = asm_prm_bc_hphob(bnd_idx)
                  ext_cff_mss_aer_lcl(2) = ext_cff_mss_bc_hphob(bnd_idx)
                  ! aerosol species 3 optical properties, hydrophilic OC
                  ss_alb_aer_lcl(3)      = ss_alb_oc_hphil(bnd_idx)
                  asm_prm_aer_lcl(3)     = asm_prm_oc_hphil(bnd_idx)
                  ext_cff_mss_aer_lcl(3) = ext_cff_mss_oc_hphil(bnd_idx)
                  ! aerosol species 4 optical properties, hydrophobic OC
                  ss_alb_aer_lcl(4)      = ss_alb_oc_hphob(bnd_idx)
                  asm_prm_aer_lcl(4)     = asm_prm_oc_hphob(bnd_idx)
                  ext_cff_mss_aer_lcl(4) = ext_cff_mss_oc_hphob(bnd_idx)

                  ! Optics for BC/dust-snow external mixing:
                  ! aerosol species 1 optical properties, hydrophilic BC
                  ss_alb_aer_lcl(1)      = ss_alb_bc_hphil(bnd_idx)
                  asm_prm_aer_lcl(1)     = asm_prm_bc_hphil(bnd_idx)
                  ext_cff_mss_aer_lcl(1) = ext_cff_mss_bc_hphil(bnd_idx)
                  ! aerosol species 5 optical properties, dust size1
                  ss_alb_aer_lcl(5)      = ss_alb_dst1(bnd_idx)
                  asm_prm_aer_lcl(5)     = asm_prm_dst1(bnd_idx)
                  ext_cff_mss_aer_lcl(5) = ext_cff_mss_dst1(bnd_idx)
                  ! aerosol species 6 optical properties, dust size2
                  ss_alb_aer_lcl(6)      = ss_alb_dst2(bnd_idx)
                  asm_prm_aer_lcl(6)     = asm_prm_dst2(bnd_idx)
                  ext_cff_mss_aer_lcl(6) = ext_cff_mss_dst2(bnd_idx)
                  ! aerosol species 7 optical properties, dust size3
                  ss_alb_aer_lcl(7)      = ss_alb_dst3(bnd_idx)
                  asm_prm_aer_lcl(7)     = asm_prm_dst3(bnd_idx)
                  ext_cff_mss_aer_lcl(7) = ext_cff_mss_dst3(bnd_idx)
                  ! aerosol species 8 optical properties, dust size4
                  ss_alb_aer_lcl(8)      = ss_alb_dst4(bnd_idx)
                  asm_prm_aer_lcl(8)     = asm_prm_dst4(bnd_idx)
                  ext_cff_mss_aer_lcl(8) = ext_cff_mss_dst4(bnd_idx)

                  ! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
                  ! 2. optical Depths (tau_snw, tau_aer)
                  ! 3. weighted Mie properties (tau, omega, g)

                  wvl_doint = wvl_ct(bnd_idx)

                  ! Weighted Mie parameters of each layer
                  do i=snl_top,snl_btm,1

                     ! Start BC/dust-snow internal mixing for wavelength<=1.2um
                     if (wvl_doint <= 1.2_r8) then
 
                        ! BC-snow internal mixing applied to hydrophilic BC if activated
                        ! BC-snow internal mixing primarily affect snow single-scattering albedo
                        if ( snicar_snobc_intmix .and. (mss_cnc_aer_lcl(i,1) > 0._r8) ) then
                           ! result from Eq.8b in He et al.(2017) is based on BC Re=0.1um &
                           ! MAC=6.81 m2/g (@550 nm) & BC density=1.7g/cm3 (den_bc).
                           ! To be consistent with Bond et al. 2006 recommeded value (BC MAC=7.5 m2/g @550nm)
                           ! we made adjustments on BC size & density as follows to get MAC=7.5m2/g:
                           ! (1) We use BC Re=0.045um [geometric mean diameter=0.06um (Dentener et al.2006, 
                           ! Yu and Luo,2009) & geometric std=1.5 (Flanner et al.2007;Aoki et al., 2011)].
                           ! (2) We tune BC density from 1.7 to 1.49 g/cm3 (den_bc_target) (Aoki et al., 2011).
                           ! These adjustments also lead to consistent results with Flanner et al. 2012 (ACP) lookup table
                           ! for BC-snow internal mixing enhancement in albedo reduction (He et al. 2018 ACP)
                           do ibb=1,sixteen_bands
                              enh_omg_bcint_tmp(ibb) = bcint_d0(ibb) * &
                                 ( (mss_cnc_aer_lcl(i,1) * kg_to_ug * den_bc / den_bc_target + bcint_d2(ibb))**bcint_d1(ibb) )
                              ! adjust enhancment factor for BC effective size from 0.1um to Re_bc (He et al. 2018 GRL Eqs.1a,1b)
                              if (ibb < 3) then ! near-UV
                                 bcint_m_tmp = bcint_m(1)
                                 bcint_n_tmp = bcint_n(1)
                              else if (ibb >= 3 .and. ibb <= 11) then ! visible
                                 bcint_m_tmp = bcint_m(2)
                                 bcint_n_tmp = bcint_n(2)
                              else  ! ibb > 11, NIR
                                 bcint_m_tmp = bcint_m(3)
                                 bcint_n_tmp = bcint_n(3)
                              endif
                              bcint_dd  = (Re_bc / radius_2)**bcint_m_tmp
                              bcint_dd2 = (radius_1 / radius_2)**bcint_m_tmp
                              bcint_f  = (Re_bc / radius_1)**bcint_n_tmp

                              enh_omg_bcint_tmp2(ibb)=LOG10(max(1._r8,bcint_dd*((enh_omg_bcint_tmp(ibb)/bcint_dd2)**bcint_f)))
                           enddo
                           ! piecewise linear interpolate into targeted SNICAR bands in a logscale space
                           call piecewise_linear_interp1d(sixteen_bands,bcint_wvl_ct,enh_omg_bcint_tmp2,wvl_doint,enh_omg_bcint_intp)
                           ! update snow single-scattering albedo
                           enh_omg_bcint_intp2 = 10._r8 ** enh_omg_bcint_intp                           
                           enh_omg_bcint_intp2 = min(enh_omg_max, max(enh_omg_bcint_intp2, 1._r8)) ! constrain enhancement to a reasonable range
                           ss_alb_snw_lcl(i)   = 1._r8 - (1._r8 - ss_alb_snw_lcl(i)) * enh_omg_bcint_intp2
                           ss_alb_snw_lcl(i)   = max(0.5_r8, min(ss_alb_snw_lcl(i),1._r8))
                           ! reset hydrophilic BC property to 0 since it is accounted by updated snow ss_alb above
                           ss_alb_aer_lcl(1)       = 0.0
                           asm_prm_aer_lcl(1)      = 0.0
                           ext_cff_mss_aer_lcl(1)  = 0.0
                        endif ! end if BC-snow mixing type

                        ! Dust-snow internal mixing applied to all size bins if activated
                        ! Dust-snow internal mixing primarily affect snow single-scattering albedo
                        ! default optics of externally mixed dust at 4 size bins based on effective
                        ! radius of 1.38um and sigma=2.0 with truncation to each size bin (Flanner et al. 2021 GMD)
                        ! parameterized dust-snow int mix results based on effective radius of 1.1um and sigma=2.0
                        ! from (He et al. 2019 JAMES). Thus, the parameterization can be approximately applied to
                        ! all dust size bins here.
                        tot_dst_snw_conc = (mss_cnc_aer_lcl(i,5) + mss_cnc_aer_lcl(i,6) + &
                                            mss_cnc_aer_lcl(i,7) + mss_cnc_aer_lcl(i,8)) * kg_kg_to_ppm
                        if ( snicar_snodst_intmix .and. (tot_dst_snw_conc > 0._r8) ) then
                           do idb=1, size_bins
                              enh_omg_dstint_tmp(idb) = dstint_a1(idb)+dstint_a2(idb)*(tot_dst_snw_conc**dstint_a3(idb))
                              enh_omg_dstint_tmp2(idb) = LOG10(max(enh_omg_dstint_tmp(idb),1._r8))
                           enddo
                           ! piecewise linear interpolate into targeted SNICAR bands in a logscale space
                           call piecewise_linear_interp1d(size_bins,dstint_wvl_ct,enh_omg_dstint_tmp2,wvl_doint,enh_omg_dstint_intp)
                           ! update snow single-scattering albedo
                           enh_omg_dstint_intp2 = 10._r8 ** enh_omg_dstint_intp
                           enh_omg_dstint_intp2 = min(enh_omg_max, max(enh_omg_dstint_intp2, 1._r8)) ! constrain enhancement to a reasonable range
                           ss_alb_snw_lcl(i) = 1._r8 - (1._r8 - ss_alb_snw_lcl(i)) * enh_omg_dstint_intp2
                           ss_alb_snw_lcl(i) = max(0.5_r8, min(ss_alb_snw_lcl(i),1._r8))
                           ! reset all dust optics to zero  since it is accounted by updated snow ss_alb above
                           ss_alb_aer_lcl(5:8)      = 0._r8
                           asm_prm_aer_lcl(5:8)     = 0._r8
                           ext_cff_mss_aer_lcl(5:8) = 0._r8
                        endif ! end if dust-snow internal mixing

                     endif ! end if BC/dust-snow internal mixing (bands<1.2um)

                     L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
                     tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)

                     do j=1,sno_nbr_aer
                        L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
                        tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
                     enddo

                     tau_sum   = 0._r8
                     omega_sum = 0._r8
                     g_sum     = 0._r8

                     do j=1,sno_nbr_aer
                        tau_sum    = tau_sum + tau_aer(i,j) 
                        omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
                        g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
                     enddo

                     tau(i)    = tau_sum + tau_snw(i)
                     omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
                     g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))

                  enddo  ! end do snow layers

                  ! DELTA transformations, if requested
                  if (DELTA == 1) then
                     do i=snl_top,snl_btm,1
                        g_star(i)     = g(i)/(1+g(i))
                        omega_star(i) = (1._r8 - g(i) * g(i)) * omega(i) / (1._r8 - omega(i) * (g(i) * g(i)))
                        tau_star(i)   = (1._r8 - omega(i) * (g(i) * g(i))) * tau(i)
                     enddo
                  else
                     do i=snl_top,snl_btm,1
                        g_star(i)     = g(i)
                        omega_star(i) = omega(i)
                        tau_star(i)   = tau(i)
                     enddo
                  endif
                 !--------------------------- End of snow & aerosol optics --------------------------------

                 !--------------------------- Start Adding-doubling RT solver  --------------------------------

                 ! Given input vertical profiles of optical properties, evaluate the
                 ! monochromatic Delta-Eddington adding-doubling solution

                 ! trndir, trntdr, trndif, rupdir, rupdif, rdndif are variables at the layer interface,
                 ! for snow with layers from snl_top to snl_btm there are snl_top to snl_btm+1 layer interface
                 snl_btm_itf = snl_btm + 1

                 ! initialization for layer interface
                 do i = snl_top,snl_btm_itf,1
                    trndir(i) = c0
                    trntdr(i) = c0
                    trndif(i) = c0
                    rupdir(i) = c0
                    rupdif(i) = c0
                    rdndif(i) = c0
                 enddo
                 ! initialize top interface of top layer
                 trndir(snl_top) = c1
                 trntdr(snl_top) = c1
                 trndif(snl_top) = c1
                 rdndif(snl_top) = c0

                 ! begin main level loop for snow layer interfaces except for the very bottom
                 do i = snl_top,snl_btm,1

                    ! initialize all layer apparent optical properties to 0
                    rdir  (i) = c0
                    rdif_a(i) = c0
                    rdif_b(i) = c0
                    tdir  (i) = c0
                    tdif_a(i) = c0
                    tdif_b(i) = c0
                    trnlay(i) = c0

                    ! compute next layer Delta-eddington solution only if total transmission
                    ! of radiation to the interface just above the layer exceeds trmin.
                    if (trntdr(i) > trmin ) then

                      ! delta-transformed single-scattering properties of this layer
                      ts = tau_star(i)
                      ws = omega_star(i)
                      gs = g_star(i)

                      ! Delta-Eddington solution expressions, Eq. 50: Briegleb and Light 2007
                      lm = sqrt(c3*(c1-ws)*(c1 - ws*gs))
                      ue = c1p5*(c1 - ws*gs)/lm
                      extins = max(exp_min, exp(-lm*ts))
                      ne = ((ue+c1)*(ue+c1)/extins) - ((ue-c1)*(ue-c1)*extins)

                      ! first calculation of rdif, tdif using Delta-Eddington formulas
                      ! Eq.: Briegleb 1992; alpha and gamma for direct radiation
                      rdif_a(i) = (ue * ue - c1) * (c1 / extins - extins) / ne
                      tdif_a(i) = c4*ue/ne

                      ! evaluate rdir,tdir for direct beam
                      trnlay(i) = max(exp_min, exp(-ts/mu_not))

                      ! Delta-Eddington solution expressions
                      ! Eq. 50: Briegleb and Light 2007; alpha and gamma for direct radiation
                      alp = cp75*ws*mu_not*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu_not*mu_not))
                      gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu_not*mu_not)/(c1-lm*lm*mu_not*mu_not))
                      apg = alp + gam
                      amg = alp - gam
                      rdir(i) = apg*rdif_a(i) +  amg*(tdif_a(i)*trnlay(i) - c1)
                      tdir(i) = apg*tdif_a(i) + (amg* rdif_a(i)-apg+c1)*trnlay(i)

                      ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
                      ! since Delta-Eddington rdif formula is not well-behaved (it is usually
                      ! biased low and can even be negative); use ngmax angles and gaussian
                      ! integration for most accuracy:
                      R1 = rdif_a(i) ! use R1 as temporary
                      T1 = tdif_a(i) ! use T1 as temporary
                      swt = c0
                      smr = c0
                      smt = c0
                      ! gaussian angles for the AD integral
                      do ng=1,ngmax
                         mu  = difgauspt(ng)
                         gwt = difgauswt(ng)
                         swt = swt + mu*gwt
                         trn = max(exp_min, exp(-ts/mu))
                         alp = cp75*ws*mu*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu*mu))
                         gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu*mu)/(c1-lm*lm*mu*mu))
                         apg = alp + gam
                         amg = alp - gam
                         rdr = apg*R1 + amg*T1*trn - amg
                         tdr = apg*T1 + amg*R1*trn - apg*trn + trn
                         smr = smr + mu*rdr*gwt
                         smt = smt + mu*tdr*gwt
                      enddo      ! ng
                      rdif_a(i) = smr/swt
                      tdif_a(i) = smt/swt

                      ! homogeneous layer
                      rdif_b(i) = rdif_a(i)
                      tdif_b(i) = tdif_a(i)

                    endif ! trntdr(k) > trmin

                    ! Calculate the solar beam transmission, total transmission, and
                    ! reflectivity for diffuse radiation from below at interface i,
                    ! the top of the current layer k:
                    !
                    !              layers       interface
                    !
                    !       ---------------------  i-1
                    !                i-1
                    !       ---------------------  i
                    !                 i
                    !       ---------------------

                    trndir(i+1) = trndir(i)*trnlay(i)            ! solar beam transmission from top
                    refkm1      = c1/(c1 - rdndif(i)*rdif_a(i))  ! interface multiple scattering for i-1
                    tdrrdir     = trndir(i)*rdir(i)              ! direct tran times layer direct ref
                    tdndif      = trntdr(i) - trndir(i)          ! total down diffuse = tot tran - direct tran
                    trntdr(i+1) = trndir(i)*tdir(i) + &          ! total transmission to direct beam for layers above
                                  (tdndif + tdrrdir*rdndif(i))*refkm1*tdif_a(i)
                    ! Eq. B4; Briegleb and Light 2007
                    rdndif(i+1) = rdif_b(i) + &                  ! reflectivity to diffuse radiation for layers above
                                  (tdif_b(i)*rdndif(i)*refkm1*tdif_a(i))
                    trndif(i+1) = trndif(i)*refkm1*tdif_a(i)     ! diffuse transmission to diffuse beam for layers above

                 enddo       ! end i main level loop

                 ! compute reflectivity to direct and diffuse radiation for layers
                 ! below by adding succesive layers starting from the underlying
                 ! ground and working upwards:
                 !
                 !              layers       interface
                 !
                 !       ---------------------  i
                 !                 i
                 !       ---------------------  i+1
                 !                i+1
                 !       ---------------------

                 ! set the underlying ground albedo == albedo of near-IR
                 ! unless bnd_idx < nir_bnd_bgn, for visible
                 rupdir(snl_btm_itf) = albsfc(c_idx,inir)
                 rupdif(snl_btm_itf) = albsfc(c_idx,inir)
                 if (bnd_idx < nir_bnd_bgn) then
                     rupdir(snl_btm_itf) = albsfc(c_idx,ivis)
                     rupdif(snl_btm_itf) = albsfc(c_idx,ivis)
                 endif

                 do i=snl_btm,snl_top,-1
                    ! interface scattering Eq. B5; Briegleb and Light 2007
                    refkp1 = c1/( c1 - rdif_b(i)*rupdif(i+1))
                    ! dir from top layer plus exp tran ref from lower layer, interface
                    ! scattered and tran thru top layer from below, plus diff tran ref
                    ! from lower layer with interface scattering tran thru top from below
                    rupdir(i) = rdir(i) &
                                + (        trnlay(i)  *rupdir(i+1) &
                                +  (tdir(i)-trnlay(i))*rupdif(i+1) ) * refkp1 * tdif_b(i)
                    ! dif from top layer from above, plus dif tran upwards reflected and
                    ! interface scattered which tran top from below
                    rupdif(i) = rdif_a(i) + tdif_a(i)*rupdif(i+1)*refkp1*tdif_b(i)
                 enddo       ! i

                 ! net flux (down-up) at each layer interface from the
                 ! snow top (i = snl_top) to bottom interface above land (i = snl_btm_itf)
                 ! the interface reflectivities and transmissivities required
                 ! to evaluate interface fluxes are returned from solution_dEdd;
                 ! now compute up and down fluxes for each interface, using the
                 ! combined layer properties at each interface:
                 !
                 !              layers       interface
                 !
                 !       ---------------------  i
                 !                 i
                 !       ---------------------

                 do i = snl_top, snl_btm_itf
                    ! interface scattering, Eq. 52; Briegleb and Light 2007
                    refk = c1/(c1 - rdndif(i)*rupdif(i))
                    ! dir tran ref from below times interface scattering, plus diff
                    ! tran and ref from below times interface scattering
                    ! fdirup(i) = (trndir(i)*rupdir(i) + &
                    !                 (trntdr(i)-trndir(i))  &
                    !                 *rupdif(i))*refk
                    ! dir tran plus total diff trans times interface scattering plus
                    ! dir tran with up dir ref and down dif ref times interface scattering
                    ! fdirdn(i) = trndir(i) + (trntdr(i) &
                    !               - trndir(i) + trndir(i)  &
                    !               *rupdir(i)*rdndif(i))*refk
                    ! diffuse tran ref from below times interface scattering
                    ! fdifup(i) = trndif(i)*rupdif(i)*refk
                    ! diffuse tran times interface scattering
                    ! fdifdn(i) = trndif(i)*refk

                    ! netflux, down - up
                    ! dfdir = fdirdn - fdirup
                    dfdir(i) = trndir(i) &
                              + (trntdr(i)-trndir(i)) * (c1 - rupdif(i)) * refk &
                              -  trndir(i)*rupdir(i)  * (c1 - rdndif(i)) * refk
                    if (dfdir(i) < puny) dfdir(i) = c0
                    ! dfdif = fdifdn - fdifup
                    dfdif(i) = trndif(i) * (c1 - rupdif(i)) * refk
                    if (dfdif(i) < puny) dfdif(i) = c0
                 enddo  ! k

                 ! SNICAR_AD_RT is called twice for direct and diffuse incident fluxes
                 ! direct incident
                 if (flg_slr_in == 1) then
                    albedo = rupdir(snl_top)
                    dftmp  = dfdir
                    refk   = c1/(c1 - rdndif(snl_top)*rupdif(snl_top))
                    F_sfc_pls = (trndir(snl_top)*rupdir(snl_top) + &
                                (trntdr(snl_top)-trndir(snl_top))  &
                                *rupdif(snl_top))*refk
                 !diffuse incident
                 else
                    albedo = rupdif(snl_top)
                    dftmp  = dfdif
                    refk   = c1/(c1 - rdndif(snl_top)*rupdif(snl_top))
                    F_sfc_pls = trndif(snl_top)*rupdif(snl_top)*refk
                 endif

                 ! Absorbed flux in each layer
                 do i=snl_top,snl_btm,1
                    F_abs(i) = dftmp(i)-dftmp(i+1)
                    flx_abs_lcl(i,bnd_idx) = F_abs(i)

                    ! ERROR check: negative absorption
                    if (flx_abs_lcl(i,bnd_idx) < -0.00001_r8) then
                      write (iulog,"(a,e13.6,a,i6,a,i6)") "SNICAR ERROR: negative absoption : ", &
                            flx_abs_lcl(i,bnd_idx), " at timestep: ", nstep, " at column: ", c_idx
                      write(iulog,*) "SNICAR_AD STATS: snw_rds(0)= ", snw_rds(c_idx,0)
                      write(iulog,*) "SNICAR_AD STATS: L_snw(0)= ", L_snw(0)
                      write(iulog,*) "SNICAR_AD STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
                      write(iulog,*) "SNICAR_AD STATS: soot1(0)= ", mss_cnc_aer_lcl(0,1)
                      write(iulog,*) "SNICAR_AD STATS: soot2(0)= ", mss_cnc_aer_lcl(0,2)
                      write(iulog,*) "SNICAR_AD STATS: dust1(0)= ", mss_cnc_aer_lcl(0,3)
                      write(iulog,*) "SNICAR_AD STATS: dust2(0)= ", mss_cnc_aer_lcl(0,4)
                      write(iulog,*) "SNICAR_AD STATS: dust3(0)= ", mss_cnc_aer_lcl(0,5)
                      write(iulog,*) "SNICAR_AD STATS: dust4(0)= ", mss_cnc_aer_lcl(0,6)
                      call endrun(subgrid_index=c_idx, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
                    endif
                 enddo

                 ! absobed flux by the underlying ground
                 F_btm_net = dftmp(snl_btm_itf)

                 ! note here, snl_btm_itf = 1 by snow column set up in CLM
                 flx_abs_lcl(1,bnd_idx) = F_btm_net

                 if (flg_nosnl == 1) then
                   ! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
                   !flx_abs_lcl(:,bnd_idx) = 0._r8
                   !flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net

                   ! changed on 20070408:
                   ! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
                   ! handles the case of no snow layers. Then, if a snow layer is addded between now and
                   ! SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
                    flx_abs_lcl(0,bnd_idx) = F_abs(0)
                    flx_abs_lcl(1,bnd_idx) = F_btm_net
                 endif

                 !Underflow check (we've already tripped the error condition above)
                 do i=snl_top,1,1
                    flx_abs_lcl(i,bnd_idx) = max(0._r8, flx_abs_lcl(i,bnd_idx))
                 enddo

                 F_abs_sum = 0._r8
                 do i=snl_top,snl_btm,1
                    F_abs_sum = F_abs_sum + F_abs(i)
                 enddo

                 ! no need to repeat calculations for adding-doubling solver
                 flg_dover = 0

                 !--------------------------- End of Adding-doubling RT solver  --------------------------------

               enddo !enddo while (flg_dover > 0)

               ! Energy conservation check:
               ! Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
               energy_sum = (mu_not*pi*flx_slrd_lcl(bnd_idx)) + flx_slri_lcl(bnd_idx) - (F_abs_sum + F_btm_net + F_sfc_pls)
               if (abs(energy_sum) > 0.00001_r8) then
                  write (iulog,"(a,e12.6,a,i6,a,i6)") "SNICAR ERROR: Energy conservation error of : ", energy_sum, &
                       " at timestep: ", nstep, " at column: ", c_idx
                  write(iulog,*) "F_abs_sum: ",F_abs_sum
                  write(iulog,*) "F_btm_net: ",F_btm_net
                  write(iulog,*) "F_sfc_pls: ",F_sfc_pls
                  write(iulog,*) "mu_not*pi*flx_slrd_lcl(bnd_idx): ", mu_not*pi*flx_slrd_lcl(bnd_idx)
                  write(iulog,*) "flx_slri_lcl(bnd_idx)", flx_slri_lcl(bnd_idx)
                  write(iulog,*) "bnd_idx", bnd_idx
                  write(iulog,*) "F_abs", F_abs
                  write(iulog,*) "albedo", albedo
                  call endrun(subgrid_index=c_idx, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
               endif

               albout_lcl(bnd_idx) = albedo

               ! Fail if albedo > 1
               if (albout_lcl(bnd_idx) > 1.0) then

                  write (iulog,*) "SNICAR ERROR: Albedo > 1.0 at c: ", c_idx, " NSTEP= ",nstep
                  write (iulog,*) "SNICAR STATS: bnd_idx= ",bnd_idx
                  write (iulog,*) "SNICAR STATS: albout_lcl(bnd)= ",albout_lcl(bnd_idx), &
                       " albsfc_lcl(bnd_idx)= ",albsfc_lcl(bnd_idx)
                  write (iulog,*) "SNICAR STATS: landtype= ", sfctype
                  write (iulog,*) "SNICAR STATS: h2osno_total= ", h2osno_lcl, " snl= ", snl_lcl
                  write (iulog,*) "SNICAR STATS: coszen= ", coszen(c_idx), " flg_slr= ", flg_slr_in

                  write (iulog,*) "SNICAR STATS: soot(-4)= ", mss_cnc_aer_lcl(-4,1)
                  write (iulog,*) "SNICAR STATS: soot(-3)= ", mss_cnc_aer_lcl(-3,1)
                  write (iulog,*) "SNICAR STATS: soot(-2)= ", mss_cnc_aer_lcl(-2,1)
                  write (iulog,*) "SNICAR STATS: soot(-1)= ", mss_cnc_aer_lcl(-1,1)
                  write (iulog,*) "SNICAR STATS: soot(0)= ", mss_cnc_aer_lcl(0,1)

                  write (iulog,*) "SNICAR STATS: L_snw(-4)= ", L_snw(-4)
                  write (iulog,*) "SNICAR STATS: L_snw(-3)= ", L_snw(-3)
                  write (iulog,*) "SNICAR STATS: L_snw(-2)= ", L_snw(-2)
                  write (iulog,*) "SNICAR STATS: L_snw(-1)= ", L_snw(-1)
                  write (iulog,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)

                  write (iulog,*) "SNICAR STATS: snw_rds(-4)= ", snw_rds(c_idx,-4)
                  write (iulog,*) "SNICAR STATS: snw_rds(-3)= ", snw_rds(c_idx,-3)
                  write (iulog,*) "SNICAR STATS: snw_rds(-2)= ", snw_rds(c_idx,-2)
                  write (iulog,*) "SNICAR STATS: snw_rds(-1)= ", snw_rds(c_idx,-1)
                  write (iulog,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)

                  call endrun(subgrid_index=c_idx, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
               endif

            enddo   ! loop over wvl bands


            ! Weight output NIR albedo appropriately
            select case (snicar_numrad_snw)
            case (default_number_bands)  ! 5-band case
              ! VIS band
              albout(c_idx,ivis) = albout_lcl(ivis)
            case (high_number_bands)  ! 480-band case
              ! average for VIS band
              flx_sum = 0._r8
              do bnd_idx= 1, (nir_bnd_bgn-1)
                 flx_sum = flx_sum + flx_wgt(bnd_idx)*albout_lcl(bnd_idx)
              end do
              albout(c_idx,ivis) = flx_sum / sum(flx_wgt(1:(nir_bnd_bgn-1)))
            end select

            ! average for NIR band (5 or 480-band case)
            flx_sum = 0._r8
            do bnd_idx = nir_bnd_bgn, nir_bnd_end
               flx_sum = flx_sum + flx_wgt(bnd_idx) * albout_lcl(bnd_idx)
            end do
            albout(c_idx,inir) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))

            ! Weight output NIR absorbed layer fluxes (flx_abs) appropriately
            select case (snicar_numrad_snw)
            case (default_number_bands)  ! 5-band case
              ! VIS band
              flx_abs(c_idx,:,1) = flx_abs_lcl(:,1)
            case (high_number_bands)  ! 480-band case
              ! average for VIS band
              do i=snl_top,1,1
                 flx_sum = 0._r8
                 do bnd_idx= 1,(nir_bnd_bgn-1)
                    flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_lcl(i,bnd_idx)
                 enddo
                 flx_abs(c_idx,i,ivis) = flx_sum / sum(flx_wgt(1:(nir_bnd_bgn-1)))
              end do
            end select

            ! average for NIR band (5 or 480-band case)
            do i = snl_top, 1, 1
               flx_sum = 0._r8
               do bnd_idx = nir_bnd_bgn, nir_bnd_end
                  flx_sum = flx_sum + flx_wgt(bnd_idx) * flx_abs_lcl(i,bnd_idx)
               end do
               flx_abs(c_idx,i,inir) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
            end do

            ! high solar zenith angle adjustment for Adding-doubling solver results
            ! near-IR direct albedo/absorption adjustment for high solar zenith angles
            ! solar zenith angle parameterization
            ! calculate the scaling factor for NIR direct albedo if SZA>75 degree
            if ((mu_not < mu_75) .and. (flg_slr_in == 1)) then
               sza_c1 = sza_a0 + sza_a1 * mu_not + sza_a2 * (mu_not * mu_not)
               sza_c0 = sza_b0 + sza_b1 * mu_not + sza_b2 * (mu_not * mu_not)
               sza_factor = sza_c1 * (log10(snw_rds_lcl(snl_top) * c1) - c6) + sza_c0
               flx_sza_adjust  = albout(c_idx,inir) * (sza_factor-c1) * sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
               albout(c_idx,inir) = albout(c_idx,inir) * sza_factor
               flx_abs(c_idx,snl_top,inir) = flx_abs(c_idx,snl_top,inir) - flx_sza_adjust
            endif


         ! If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
         elseif ( (coszen(c_idx) > 0._r8) .and. (h2osno_lcl < min_snw) .and. (h2osno_lcl > 0._r8) ) then
            albout(c_idx,ivis) = albsfc(c_idx,ivis)
            albout(c_idx,inir) = albsfc(c_idx,inir)

         ! There is either zero snow, or no sun
         else
            albout(c_idx,ivis) = 0._r8
            albout(c_idx,inir) = 0._r8
         endif    ! if column has snow and coszen > 0

      enddo    ! loop over all columns

    end associate

  end subroutine SNICAR_RT

  !-----------------------------------------------------------------------
  subroutine SnowAge_grain(bounds, &
       num_snowc, filter_snowc, num_nosnowc, filter_nosnowc, &
       waterfluxbulk_inst, waterstatebulk_inst, waterdiagnosticbulk_inst, temperature_inst, atm2lnd_inst)
    !
    ! !DESCRIPTION:
    ! Updates the snow effective grain size (radius). 
    ! Contributions to grain size evolution are from:
    !   1. vapor redistribution (dry snow) 
    !   2. liquid water redistribution (wet snow)
    !   3. re-freezing of liquid water
    ! 
    ! Vapor redistribution: Method is to retrieve 3 best-bit parameters that
    ! depend on snow temperature, temperature gradient, and density,
    ! that are derived from the microphysical model described in: 
    ! Flanner and Zender (2006), Linking snowpack microphysics and albedo
    ! evolution, J. Geophys. Res., 111, D12208, doi:10.1029/2005JD006834. 
    ! The parametric equation has the form: 
    ! dr/dt = drdt_0*(tau/(dr_fresh+tau))^(1/kappa), where:
    !   r is the effective radius,
    !   tau and kappa are best-fit parameters,
    !   drdt_0 is the initial rate of change of effective radius, and
    !   dr_fresh is the difference between the current and fresh snow states 
    !  (r_current - r_fresh).
    !
    ! Liquid water redistribution: Apply the grain growth function from:
    !   Brun, E. (1989), Investigation of wet-snow metamorphism in respect of 
    !   liquid-water content, Annals of Glaciology, 13, 22-26.
    !   There are two parameters that describe the grain growth rate as 
    !   a function of snow liquid water content (LWC). The "LWC=0" parameter
    !   is zeroed here because we are accounting for dry snowing with a 
    !   different representation
    !
    ! Re-freezing of liquid water: Assume that re-frozen liquid water clumps
    !   into an arbitrarily large effective grain size (snw_rds_refrz). 
    !   The phenomenon is observed (Grenfell), but so far unquantified, as far as 
    !   I am aware.
    !
    ! !USES:
    use clm_time_manager , only : get_step_size_real, get_nstep
    use clm_varpar       , only : nlevsno
    use clm_varcon, only: spval, secsphr
    use shr_const_mod    , only : SHR_CONST_RHOICE, SHR_CONST_PI
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_snowc         ! number of column snow points in column filter
    integer                , intent(in)    :: filter_snowc(:)   ! column filter for snow points
    integer                , intent(in)    :: num_nosnowc       ! number of column non-snow points in column filter
    integer                , intent(in)    :: filter_nosnowc(:) ! column filter for non-snow points
    type(waterfluxbulk_type)   , intent(in)    :: waterfluxbulk_inst
    type(waterstatebulk_type)  , intent(inout) :: waterstatebulk_inst
    type(waterdiagnosticbulk_type)  , intent(inout) :: waterdiagnosticbulk_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    !
    ! !LOCAL VARIABLES:
    integer :: snl_top                      ! top snow layer index [idx]
    integer :: snl_btm                      ! bottom snow layer index [idx]
    integer :: i                            ! layer index [idx]
    integer :: c_idx                        ! column index [idx]
    integer :: fc                           ! snow column filter index [idx]
    integer :: T_idx                        ! snow aging lookup table temperature index [idx]
    integer :: Tgrd_idx                     ! snow aging lookup table temperature gradient index [idx]
    integer :: rhos_idx                     ! snow aging lookup table snow density index [idx]
    real(r8) :: t_snotop                    ! temperature at upper layer boundary [K]
    real(r8) :: t_snobtm                    ! temperature at lower layer boundary [K]
    real(r8) :: dTdz(bounds%begc:bounds%endc,-nlevsno:0)    ! snow temperature gradient (col,lyr) [K m-1]
    real(r8) :: bst_tau                     ! snow aging parameter retrieved from lookup table [hour]
    real(r8) :: bst_kappa                   ! snow aging parameter retrieved from lookup table [unitless]
    real(r8) :: bst_drdt0                   ! snow aging parameter retrieved from lookup table [um hr-1]
    real(r8) :: dr                          ! incremental change in snow effective radius [um]
    real(r8) :: dr_wet                      ! incremental change in snow effective radius from wet growth [um]
    real(r8) :: dr_fresh                    ! difference between fresh snow r_e and current r_e [um]
    real(r8) :: newsnow                     ! fresh snowfall [kg m-2]
    real(r8) :: refrzsnow                   ! re-frozen snow [kg m-2]
    real(r8) :: frc_newsnow                 ! fraction of layer mass that is new snow [frc]
    real(r8) :: frc_oldsnow                 ! fraction of layer mass that is old snow [frc]
    real(r8) :: frc_refrz                   ! fraction of layer mass that is re-frozen snow [frc]
    real(r8) :: frc_liq                     ! fraction of layer mass that is liquid water[frc]    
    real(r8) :: dtime                       ! land model time step [sec]
    real(r8) :: rhos                        ! snow density [kg m-3]
    real(r8) :: h2osno_lyr                  ! liquid + solid H2O in snow layer [kg m-2]
    real(r8) :: cdz(-nlevsno+1:0)           ! column average layer thickness [m]
    real(r8) :: snw_rds_fresh               ! fresh snow radius [microns]
    !--------------------------------------------------------------------------!

    associate(                                                      & 
         snl                => col%snl                            , & ! Input:  [integer  (:)   ]  negative number of snow layers (col) [nbr]
         dz                 => col%dz                             , & ! Input:  [real(r8) (:,:) ]  layer thickness (col,lyr) [m]         

         qflx_snow_grnd_col => waterfluxbulk_inst%qflx_snow_grnd_col  , & ! Input:  [real(r8) (:)   ]  snow on ground after interception (col) [kg m-2 s-1]
         qflx_snofrz_lyr    => waterfluxbulk_inst%qflx_snofrz_lyr_col , & ! Input:  [real(r8) (:,:) ]  snow freezing rate (col,lyr) [kg m-2 s-1]

         frac_sno           => waterdiagnosticbulk_inst%frac_sno_eff_col   , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         h2osno_no_layers   => waterstatebulk_inst%h2osno_no_layers_col    , & ! Input:  [real(r8) (:)   ]  snow that is not resolved into layers (col) [mm H2O]
         h2osoi_liq         => waterstatebulk_inst%h2osoi_liq_col     , & ! Input:  [real(r8) (:,:) ]  liquid water content (col,lyr) [kg m-2]
         h2osoi_ice         => waterstatebulk_inst%h2osoi_ice_col     , & ! Input:  [real(r8) (:,:) ]  ice content (col,lyr) [kg m-2]        
         snw_rds            => waterdiagnosticbulk_inst%snw_rds_col        , & ! Output: [real(r8) (:,:) ]  effective grain radius (col,lyr) [microns, m-6]
         snw_rds_top        => waterdiagnosticbulk_inst%snw_rds_top_col    , & ! Output: [real(r8) (:)   ]  effective grain radius, top layer (col) [microns, m-6]
         sno_liq_top        => waterdiagnosticbulk_inst%sno_liq_top_col    , & ! Output: [real(r8) (:)   ]  liquid water fraction (mass) in top snow layer (col) [frc]

         t_soisno           => temperature_inst%t_soisno_col      , & ! Input:  [real(r8) (:,:) ]  soil and snow temperature (col,lyr) [K]
         t_grnd             => temperature_inst%t_grnd_col        , & ! Input:  [real(r8) (:)   ]  ground temperature (col) [K]            
         snot_top           => temperature_inst%snot_top_col      , & ! Output: [real(r8) (:)   ]  temperature in top snow layer (col) [K]            
         dTdz_top           => temperature_inst%dTdz_top_col        & ! Output: [real(r8) (:)   ]  temperature gradient in top layer (col) [K m-1]
         )
  

      ! set timestep and step interval
      dtime = get_step_size_real()

      ! loop over columns that have at least one snow layer
      do fc = 1, num_snowc
         c_idx = filter_snowc(fc)

         snl_btm = 0
         snl_top = snl(c_idx) + 1

         cdz(snl_top:snl_btm)=frac_sno(c_idx)*dz(c_idx,snl_top:snl_btm)

         ! loop over snow layers
         do i=snl_top,snl_btm,1
            !
            !**********  1. DRY SNOW AGING  ***********
            !
            h2osno_lyr = h2osoi_liq(c_idx,i) + h2osoi_ice(c_idx,i)

            ! temperature gradient
            if (i == snl_top) then 
               ! top layer
               t_snotop = t_soisno(c_idx,snl_top)
               t_snobtm = (t_soisno(c_idx,i+1)*dz(c_idx,i) &
                    + t_soisno(c_idx,i)*dz(c_idx,i+1)) &
                    / (dz(c_idx,i)+dz(c_idx,i+1))
            else
               t_snotop = (t_soisno(c_idx,i-1)*dz(c_idx,i) &
                    + t_soisno(c_idx,i)*dz(c_idx,i-1)) &
                    / (dz(c_idx,i)+dz(c_idx,i-1))
               t_snobtm = (t_soisno(c_idx,i+1)*dz(c_idx,i) &
                    + t_soisno(c_idx,i)*dz(c_idx,i+1)) &
                    / (dz(c_idx,i)+dz(c_idx,i+1))
            endif

            dTdz(c_idx,i) = abs((t_snotop - t_snobtm) / cdz(i))

            ! snow density
            rhos = (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i)) / cdz(i)

            ! make sure rhos doesn't drop below 50 (see rhos_idx below)
            rhos=max(50._r8,rhos)

            ! best-fit table indices
            T_idx    = nint((t_soisno(c_idx,i)-223) / 5) + 1
            Tgrd_idx = nint(dTdz(c_idx,i) / 10) + 1
            rhos_idx = nint((rhos-50) / 50) + 1

            ! boundary checks
            T_idx = max(T_idx, idx_T_min)
            T_idx = min(T_idx, idx_T_max)
            Tgrd_idx = max(Tgrd_idx, idx_Tgrd_min)
            Tgrd_idx = min(Tgrd_idx, idx_Tgrd_max)
            rhos_idx = max(rhos_idx, idx_rhos_min)
            rhos_idx = min(rhos_idx, idx_rhos_max)

            ! best-fit parameters
            bst_tau   = snowage_tau(rhos_idx,Tgrd_idx,T_idx)
            bst_kappa = snowage_kappa(rhos_idx,Tgrd_idx,T_idx)     
            bst_drdt0 = snowage_drdt0(rhos_idx,Tgrd_idx,T_idx)


            !LvK extra boundary check, to prevent when using old restart file with lower snw_rds_min than current run
            snw_rds(c_idx,i) = max(snw_rds(c_idx,i), params_inst%snw_rds_min)

            ! change in snow effective radius, using best-fit parameters
            dr_fresh = snw_rds(c_idx,i) - params_inst%snw_rds_min
            dr = (bst_drdt0 * (bst_tau / (dr_fresh + bst_tau))**(1._r8 / bst_kappa)) * (dtime / secsphr)

            !
            !**********  2. WET SNOW AGING  ***********
            !
            ! We are assuming wet and dry evolution occur simultaneously, and 
            ! the contributions from both can be summed. 
            ! This is justified by setting the linear offset constant C1_liq_Brun89 to zero [Brun, 1989]

            ! liquid water faction
            frc_liq = min(0.1_r8, (h2osoi_liq(c_idx,i) / (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i))))

            !dr_wet = 1E6_r8*(dtime*(C1_liq_Brun89 + C2_liq_Brun89*(frc_liq**(3))) / (4*SHR_CONST_PI*(snw_rds(c_idx,i)/1E6)**(2)))
            !simplified, units of microns:
            dr_wet = 1E18_r8*(dtime*(params_inst%C2_liq_Brun89*(frc_liq**(3))) / &
               (4._r8 * SHR_CONST_PI * (snw_rds(c_idx,i) * snw_rds(c_idx,i))))

            dr = dr + dr_wet

            !
            !**********  3. SNOWAGE SCALING  *************
            !
            ! Multiply rate of change of effective radius by some constant, xdrdt
            dr = dr*params_inst%xdrdt

            !
            !**********  4. INCREMENT EFFECTIVE RADIUS, ACCOUNTING FOR:  ***********
            !               DRY AGING
            !               WET AGING
            !               FRESH SNOW
            !               RE-FREEZING
            !
            ! new snowfall [kg/m2]
            newsnow = max(0._r8, (qflx_snow_grnd_col(c_idx)*dtime))

            ! snow that has re-frozen [kg/m2]
            refrzsnow = max(0._r8, (qflx_snofrz_lyr(c_idx,i)*dtime))

            ! fraction of layer mass that is re-frozen
            frc_refrz = refrzsnow / h2osno_lyr

            ! fraction of layer mass that is new snow
            if (i == snl_top) then
               frc_newsnow = newsnow / h2osno_lyr
            else
               frc_newsnow = 0._r8
            endif

            if ((frc_refrz + frc_newsnow) > 1._r8) then
               frc_refrz = frc_refrz / (frc_refrz + frc_newsnow)
               frc_newsnow = 1._r8 - frc_refrz
               frc_oldsnow = 0._r8
            else
               frc_oldsnow = 1._r8 - frc_refrz - frc_newsnow
            endif

            ! temperature dependent fresh grain size
            snw_rds_fresh = FreshSnowRadius(c_idx, atm2lnd_inst)

            ! mass-weighted mean of fresh snow, old snow, and re-frozen snow effective radius
            snw_rds(c_idx,i) = (snw_rds(c_idx,i)+dr)*frc_oldsnow + snw_rds_fresh*frc_newsnow + &
                               params_inst%snw_rds_refrz*frc_refrz
            !
            !**********  5. CHECK BOUNDARIES   ***********
            !
            ! boundary check
            snw_rds(c_idx,i) = max(snw_rds(c_idx,i), params_inst%snw_rds_min)
            snw_rds(c_idx,i) = min(snw_rds(c_idx,i), snw_rds_max)

            ! set top layer variables for history files
            if (i == snl_top) then
               snot_top(c_idx)    = t_soisno(c_idx,i)
               dTdz_top(c_idx)    = dTdz(c_idx,i)
               snw_rds_top(c_idx) = snw_rds(c_idx,i)
               sno_liq_top(c_idx) = h2osoi_liq(c_idx,i) / (h2osoi_liq(c_idx,i)+h2osoi_ice(c_idx,i))
            endif

         enddo
      enddo

      ! Special case: snow on ground, but not enough to have defined a snow layer:
      !   set snw_rds to fresh snow grain size:
      do fc = 1, num_nosnowc
         c_idx = filter_nosnowc(fc)
         if (h2osno_no_layers(c_idx) > 0._r8) then
            snw_rds(c_idx,0) = params_inst%snw_rds_min
         endif
      enddo

    end associate 

  end subroutine SnowAge_grain


  !-----------------------------------------------------------------------
  real(r8) function FreshSnowRadius(c_idx, atm2lnd_inst) 
    !
    ! !DESCRIPTION:
    ! Returns fresh snow grain radius, which is linearly dependent on temperature.
    ! This is implemented to remedy an outstanding bias that SNICAR has in initial
    ! grain size. See e.g. Sandells et al, 2017 for a discussion (10.5194/tc-11-229-2017).
    !
    ! Yang et al. (2017), 10.1016/j.jqsrt.2016.03.033
    !  discusses grain size observations, which suggest a temperature dependence. 
    ! 
    ! !REVISION HISTORY:
    ! Author: Leo VanKampenhout
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in)                :: c_idx         ! column index
    type(atm2lnd_type) , intent(in)    :: atm2lnd_inst  ! Forcing from atmosphere
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------
    real(r8), parameter :: tmin = tfrz - 30._r8       ! start of linear ramp
    real(r8), parameter :: tmax = tfrz - 0._r8        ! end of linear ramp
    real(r8) :: gs_min  ! minimum value
    real(r8) :: gs_max  ! maximum value

    associate( &
         forc_t      => atm2lnd_inst%forc_t_downscaled_col & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)
         )
       if ( params_inst%fresh_snw_rds_max <= params_inst%snw_rds_min )then
           FreshSnowRadius = params_inst%snw_rds_min
       else
           gs_max = params_inst%fresh_snw_rds_max
           gs_min = params_inst%snw_rds_min

           if (forc_t(c_idx) < tmin) then
               FreshSnowRadius = gs_min
           else if (forc_t(c_idx) > tmax) then
               FreshSnowRadius = gs_max
           else
               FreshSnowRadius = (tmax-forc_t(c_idx))/(tmax-tmin)*gs_min + &
                                 (forc_t(c_idx)-tmin)/(tmax-tmin)*gs_max
           end if
       end if

    end associate

  end function FreshSnowRadius



  !-----------------------------------------------------------------------
   subroutine SnowOptics_init( )
     
     use fileutils  , only : getfil
     use CLM_varctl , only : fsnowoptics, snicar_numrad_snw
     use CLM_varctl , only : snicar_solarspec, snicar_dust_optics
     use spmdMod    , only : masterproc
     use ncdio_pio  , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile

     type(file_desc_t)  :: ncid                        ! netCDF file id
     character(len=256) :: locfn                       ! local filename
     character(len= 32) :: subname = 'SnowOptics_init' ! subroutine name
     integer            :: ier                         ! error status
     logical :: readv  ! has variable been read in or not
     character(len=100) :: errCode = '-Error reading fsnowoptics file:'
     character(len=100) :: tString ! temp. var for reading
     character(len=3) :: short_case_dust_opt  ! subset of tString
     character(len=3) :: short_case_solarspec  ! subset of tString

     !
     ! Initialize optical variables
     allocate(ss_alb_snw_drc(idx_Mie_snw_mx,snicar_numrad_snw))
     allocate(asm_prm_snw_drc(idx_Mie_snw_mx,snicar_numrad_snw))
     allocate(ext_cff_mss_snw_drc(idx_Mie_snw_mx,snicar_numrad_snw))
     allocate(ss_alb_snw_dfs(idx_Mie_snw_mx,snicar_numrad_snw))
     allocate(asm_prm_snw_dfs(idx_Mie_snw_mx,snicar_numrad_snw))
     allocate(ext_cff_mss_snw_dfs(idx_Mie_snw_mx,snicar_numrad_snw))
     allocate(ss_alb_bc_hphil(snicar_numrad_snw))
     allocate(asm_prm_bc_hphil(snicar_numrad_snw))
     allocate(ext_cff_mss_bc_hphil(snicar_numrad_snw))
     allocate(ss_alb_bc_hphob(snicar_numrad_snw))
     allocate(asm_prm_bc_hphob(snicar_numrad_snw))
     allocate(ext_cff_mss_bc_hphob(snicar_numrad_snw))
     allocate(ss_alb_oc_hphil(snicar_numrad_snw))
     allocate(asm_prm_oc_hphil(snicar_numrad_snw))
     allocate(ext_cff_mss_oc_hphil(snicar_numrad_snw))
     allocate(ss_alb_oc_hphob(snicar_numrad_snw))
     allocate(asm_prm_oc_hphob(snicar_numrad_snw))
     allocate(ext_cff_mss_oc_hphob(snicar_numrad_snw))
     allocate(ss_alb_dst1(snicar_numrad_snw))
     allocate(asm_prm_dst1(snicar_numrad_snw))
     allocate(ext_cff_mss_dst1(snicar_numrad_snw))
     allocate(ss_alb_dst2(snicar_numrad_snw))
     allocate(asm_prm_dst2(snicar_numrad_snw))
     allocate(ext_cff_mss_dst2(snicar_numrad_snw))
     allocate(ss_alb_dst3(snicar_numrad_snw))
     allocate(asm_prm_dst3(snicar_numrad_snw))
     allocate(ext_cff_mss_dst3(snicar_numrad_snw))
     allocate(ss_alb_dst4(snicar_numrad_snw))
     allocate(asm_prm_dst4(snicar_numrad_snw))
     allocate(ext_cff_mss_dst4(snicar_numrad_snw))
     allocate(flx_wgt_dir(snicar_numrad_snw))
     allocate(flx_wgt_dif(snicar_numrad_snw))

     if (masterproc) write(iulog,*) 'Attempting to read snow optical properties...'
     call getfil (fsnowoptics, locfn, 0)
     call ncd_pio_openfile(ncid, locfn, 0)
     if(masterproc) write(iulog,*) subname,trim(fsnowoptics)

     select case (snicar_solarspec)
     case ('mid_latitude_winter')  ! mid-latitude winter spectrum
        short_case_solarspec = 'mlw'
     case ('mid_latitude_summer')  ! mid-latitude summer spectrum
        short_case_solarspec = 'mls'
     case ('sub_arctic_winter')  ! sub-Arctic winter spectrum
        short_case_solarspec = 'saw'
     case ('sub_arctic_summer')  ! sub-Arctic summer spectrum
        short_case_solarspec = 'sas'
     case ('summit_greenland_summer')  ! Summit,Greenland,summer spectrum
        short_case_solarspec = 'smm'
     case ('high_mountain_summer')  ! High Mountain summer spectrum
        short_case_solarspec = 'hmn'
     case default
        write(iulog,*) subname//' ERROR: unknown snicar_solarspec: ', snicar_solarspec
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end select

     select case (snicar_dust_optics)  ! dust optical properties
     case ('sahara')  ! Saharan dust (Balkanski et al., 2007, central hematite)
        short_case_dust_opt = 'sah'
     case ('san_juan_mtns_colorado')  ! San Juan Mountains, CO (Skiles et al, 2017)
        short_case_dust_opt = 'col'
     case ('greenland')  ! Greenland (Polashenski et al., 2015, central absorptivity)
        short_case_dust_opt = 'gre'
     case default
        write(iulog,*) subname//' ERROR: unknown snicar_dust_optics: ', snicar_dust_optics
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end select

     !--------------------- for 5-band data
     select case (snicar_numrad_snw)
     case (default_number_bands)  ! 5-band case

        ! The argument posNOTonfile=.true. is used here because this is a non-spatial file.
        ! This argument is relevant when running single_column.
        ! flux weights/spectrum
        tString = 'flx_wgt_dir5_'//short_case_solarspec
        call ncd_io(trim(tString), flx_wgt_dir, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'flx_wgt_dif5_'//short_case_solarspec
        call ncd_io(trim(tString), flx_wgt_dif, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        !
        ! THIS NOTE APPLIES TO ALL THE call ncd_io LINES BELOW WHERE
        ! bcphob AND ocphob GET ASSIGNED TO VARIABLES SUFFIXED bc_hphil/oc_hphil:
        !
        ! Assumption (1) applies here, in the input section.
        ! Assumption (2) applies later, in the snicar code.
        !
        ! 1) In this section, hydrophillic particles behave like hydrophobic
        ! particles. We assume bc_hphil/oc_hphil to have the same optics as bc_hphob/oc_hphob
        ! because sulfate coating on the bc_hphil/oc_hphil surface is assumed to be
        ! dissolved into the hydrometeo (i.e, snow grain here) during the
        ! deposition process. This is different from the assumption made in
        ! prior model versions, where bc_hphil/oc_hphil was coated by undissolved
        ! sulfate.
        ! 2) Later, in the snicar code, if the bc-snow internal mixing option
        ! is on, bc_hphil/oc_hphil (internally mixed within the snow grain) will be
        ! treated differently than bc_hphob/oc_hphob (mixed externally or outside the
        ! snow grain).
        !
        ! BC species 1 Mie parameters, uncoated BC, same as bc_hphob before BC-snow internal mixing
        tString = 'ss_alb_bcphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_bc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_bcphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), asm_prm_bc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_bcphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_bc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! BC species 2 Mie parameters, uncoated BC
        tString = 'ss_alb_bcphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_bc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_bcphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), asm_prm_bc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_bcphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_bc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! OC species 1 Mie parameters, uncoated OC, same as oc_hphob before OC-snow internal mixing
        tString = 'ss_alb_ocphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_oc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_ocphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), asm_prm_oc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_ocphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_oc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! OC species 2 Mie parameters, uncoated OC
        tString = 'ss_alb_ocphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_oc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_ocphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), asm_prm_oc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_ocphob_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_oc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! ice refractive index (Picard et al., 2016)
        tString = 'ss_alb_ice_pic16_dir_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_snw_drc, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_ice_pic16_dir_'//short_case_solarspec
        call ncd_io(trim(tString),asm_prm_snw_drc, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_ice_pic16_dir_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_snw_drc, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ss_alb_ice_pic16_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_snw_dfs, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_ice_pic16_dif_'//short_case_solarspec
        call ncd_io(trim(tString),asm_prm_snw_dfs, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_ice_pic16_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_snw_dfs, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

        ! dust species 1 Mie parameters
        tString = 'ss_alb_dust01_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_dst1, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_dust01_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), asm_prm_dst1, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_dust01_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_dst1, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! dust species 2 Mie parameters
        tString = 'ss_alb_dust02_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_dst2, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_dust02_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), asm_prm_dst2, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_dust02_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_dst2, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! dust species 3 Mie parameters
        tString = 'ss_alb_dust03_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_dst3, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_dust03_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), asm_prm_dst3, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_dust03_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_dst3, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! dust species 4 Mie parameters
        tString = 'ss_alb_dust04_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ss_alb_dst4, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_dust04_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), asm_prm_dst4, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_dust04_'//short_case_dust_opt//'_dif_'//short_case_solarspec
        call ncd_io(trim(tString), ext_cff_mss_dst4, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

     !-------------------- for 480-band data
     case (high_number_bands)

        ! BC species 1 Mie parameters, uncoated BC, same as bc_hphob before BC-snow internal mixing
        tString = 'ss_alb_bcphob'
        call ncd_io(trim(tString), ss_alb_bc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_bcphob'
        call ncd_io(trim(tString), asm_prm_bc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_bcphob'
        call ncd_io(trim(tString), ext_cff_mss_bc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! BC species 2 Mie parameters, uncoated BC
        tString = 'ss_alb_bcphob'
        call ncd_io(trim(tString), ss_alb_bc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_bcphob'
        call ncd_io(trim(tString), asm_prm_bc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_bcphob'
        call ncd_io(trim(tString), ext_cff_mss_bc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! OC species 1 Mie parameters, uncoated OC, same as oc_hphob before OC-snow internal mixing
        tString = 'ss_alb_ocphob'
        call ncd_io(trim(tString), ss_alb_oc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_ocphob'
        call ncd_io(trim(tString), asm_prm_oc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_ocphob'
        call ncd_io(trim(tString), ext_cff_mss_oc_hphil, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! OC species 2 Mie parameters, uncoated OC
        tString = 'ss_alb_ocphob'
        call ncd_io(trim(tString), ss_alb_oc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_ocphob'
        call ncd_io(trim(tString), asm_prm_oc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_ocphob'
        call ncd_io(trim(tString), ext_cff_mss_oc_hphob, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

        ! snow optical properties derived from different ice refractive index dataset
        ! same value for direct and diffuse due to high spectral res without spectra averaging in database (Picard et al., 2016)
        tString = 'ss_alb_ice_pic16'
        call ncd_io(trim(tString), ss_alb_snw_drc, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_ice_pic16'
        call ncd_io(trim(tString), asm_prm_snw_drc, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_ice_pic16'
        call ncd_io(trim(tString), ext_cff_mss_snw_drc, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ss_alb_ice_pic16'
        call ncd_io(trim(tString), ss_alb_snw_dfs, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_ice_pic16'
        call ncd_io(trim(tString), asm_prm_snw_dfs, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_ice_pic16'
        call ncd_io(trim(tString), ext_cff_mss_snw_dfs, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

        ! dust optical properties
        ! dust species 1 Mie parameters
        tString = 'ss_alb_dust01_'//short_case_dust_opt
        call ncd_io(trim(tString), ss_alb_dst1, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_dust01_'//short_case_dust_opt
        call ncd_io(trim(tString), asm_prm_dst1, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_dust01_'//short_case_dust_opt
        call ncd_io(trim(tString), ext_cff_mss_dst1, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! dust species 2 Mie parameters
        tString = 'ss_alb_dust02_'//short_case_dust_opt
        call ncd_io(trim(tString), ss_alb_dst2, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_dust02_'//short_case_dust_opt
        call ncd_io(trim(tString), asm_prm_dst2, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_dust02_'//short_case_dust_opt
        call ncd_io(trim(tString), ext_cff_mss_dst2, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! dust species 3 Mie parameters
        tString = 'ss_alb_dust03_'//short_case_dust_opt
        call ncd_io(trim(tString), ss_alb_dst3, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_dust03_'//short_case_dust_opt
        call ncd_io(trim(tString), asm_prm_dst3, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_dust03_'//short_case_dust_opt
        call ncd_io(trim(tString), ext_cff_mss_dst3, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        ! dust species 4 Mie parameters
        tString = 'ss_alb_dust04_'//short_case_dust_opt
        call ncd_io(trim(tString), ss_alb_dst4, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'asm_prm_dust04_'//short_case_dust_opt
        call ncd_io(trim(tString), asm_prm_dst4, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'ext_cff_mss_dust04_'//short_case_dust_opt
        call ncd_io(trim(tString), ext_cff_mss_dst4, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
 
        ! downward solar radiation spectral weights for 480-band
        tString = 'flx_wgt_dir480_'//short_case_solarspec
        call ncd_io(trim(tString), flx_wgt_dir, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
        tString = 'flx_wgt_dif480_'//short_case_solarspec
        call ncd_io(trim(tString), flx_wgt_dif, 'read', ncid, readv, posNOTonfile=.true.)
        if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

     case default
        write(iulog,*) subname//' ERROR: unknown snicar_numrad_snw: ', snicar_numrad_snw
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end select

     call ncd_pio_closefile(ncid)
     if (masterproc) then

        write(iulog,*) 'Successfully read snow optical properties'
        ! print some diagnostics:
        write (iulog,*) 'SNICAR: Mie single scatter albedos for direct-beam ice, rds=100um: ', &
             ss_alb_snw_drc(71,1), ss_alb_snw_drc(71,2), ss_alb_snw_drc(71,3),     &
             ss_alb_snw_drc(71,4), ss_alb_snw_drc(71,5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for diffuse ice, rds=100um: ',     &
             ss_alb_snw_dfs(71,1), ss_alb_snw_dfs(71,2), ss_alb_snw_dfs(71,3),     &
             ss_alb_snw_dfs(71,4), ss_alb_snw_dfs(71,5)
        if (do_sno_oc) then
           write (iulog,*) 'SNICAR: Including OC aerosols from snow radiative transfer calculations'
        else
           write (iulog,*) 'SNICAR: Excluding OC aerosols from snow radiative transfer calculations'
        endif
        write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic BC: ', &
             ss_alb_bc_hphil(1), ss_alb_bc_hphil(2), ss_alb_bc_hphil(3), ss_alb_bc_hphil(4), ss_alb_bc_hphil(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic BC: ', &
             ss_alb_bc_hphob(1), ss_alb_bc_hphob(2), ss_alb_bc_hphob(3), ss_alb_bc_hphob(4), ss_alb_bc_hphob(5)
        if (do_sno_oc) then
           write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic OC: ', &
                ss_alb_oc_hphil(1), ss_alb_oc_hphil(2), ss_alb_oc_hphil(3), ss_alb_oc_hphil(4), ss_alb_oc_hphil(5)
           write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic OC: ', &
                ss_alb_oc_hphob(1), ss_alb_oc_hphob(2), ss_alb_oc_hphob(3), ss_alb_oc_hphob(4), ss_alb_oc_hphob(5)
        endif
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 1: ', &
             ss_alb_dst1(1), ss_alb_dst1(2), ss_alb_dst1(3), ss_alb_dst1(4), ss_alb_dst1(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 2: ', &
             ss_alb_dst2(1), ss_alb_dst2(2), ss_alb_dst2(3), ss_alb_dst2(4), ss_alb_dst2(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 3: ', &
             ss_alb_dst3(1), ss_alb_dst3(2), ss_alb_dst3(3), ss_alb_dst3(4), ss_alb_dst3(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 4: ', &
             ss_alb_dst4(1), ss_alb_dst4(2), ss_alb_dst4(3), ss_alb_dst4(4), ss_alb_dst4(5)
        write(iulog,*)
     end if

   end subroutine SnowOptics_init

   !-----------------------------------------------------------------------
   subroutine SnowAge_init( )
     use CLM_varctl      , only : fsnowaging
     use fileutils       , only : getfil
     use spmdMod         , only : masterproc
     use ncdio_pio       , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile

     type(file_desc_t)  :: ncid                        ! netCDF file id
     character(len=256) :: locfn                       ! local filename
     character(len= 32) :: subname = 'SnowOptics_init' ! subroutine name
     integer            :: varid                       ! netCDF id's
     integer            :: ier                         ! error status
     logical :: readv  ! has variable been read in or not
     character(len=100) :: errCode = '-Error reading snow aging parameters:'
     character(len=100) :: tString ! temp. var for reading

     ! Open snow aging (effective radius evolution) file:
     allocate(snowage_tau(idx_rhos_max,idx_Tgrd_max,idx_T_max))
     allocate(snowage_kappa(idx_rhos_max,idx_Tgrd_max,idx_T_max))
     allocate(snowage_drdt0(idx_rhos_max,idx_Tgrd_max,idx_T_max))

     if(masterproc)  write(iulog,*) 'Attempting to read snow aging parameters .....'
     call getfil (fsnowaging, locfn, 0)
     call ncd_pio_openfile(ncid, locfn, 0)
     if(masterproc) write(iulog,*) subname,trim(fsnowaging)

     ! snow aging parameters

     tString = 'tau'
     call ncd_io(trim(tString), snowage_tau, 'read', ncid, readv, posNOTonfile=.true.)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
     tString = 'kappa'
     call ncd_io(trim(tString), snowage_kappa, 'read', ncid, readv, posNOTonfile=.true.)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
     tString = 'drdsdt0'
     call ncd_io(trim(tString), snowage_drdt0, 'read', ncid, readv, posNOTonfile=.true.)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

     call ncd_pio_closefile(ncid)
     if (masterproc) then

        write(iulog,*) 'Successfully read snow aging properties'

        ! print some diagnostics:
        write (iulog,*) 'SNICAR: snowage tau for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_tau(3,11,9)
        write (iulog,*) 'SNICAR: snowage kappa for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_kappa(3,11,9)
        write (iulog,*) 'SNICAR: snowage dr/dt_0 for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_drdt0(3,11,9)
     endif

   end subroutine SnowAge_init

   !-----------------------------------------------------------------------
   subroutine piecewise_linear_interp1d(nd, xd, yd, xi, yi)

     ! piecewise linear interpolation method for 1-dimensional data
     ! original author: John Burkardt, Florida State University, 09/22/2012
     ! Licencing: Original code distributed under the GNU LGPL license
     ! Original code: https://people.sc.fsu.edu/~jburkardt/f77_src/pwl_interp_1d/pwl_interp_1d.f
     ! Added and modified by Cenlin He (NCAR), 01/27/2022

     implicit none

     integer , intent(in)   :: nd         ! number of data points of (xd)
     real(r8), intent(in)   :: xd(1:nd)   ! x-value of data points
     real(r8), intent(in)   :: yd(1:nd)   ! y-value of data points
     real(r8), intent(in)   :: xi         ! x-value for to-be-interpolated point
     real(r8), intent(out)  :: yi         ! the interpolated value at xi

     ! local variables
     integer  :: i, k    ! loop index
     real(r8) :: t

     yi = 0._r8

     ! if only one data point
     if ( nd == 1 ) then
        yi = yd(1)
        return
     endif

     ! if multiple data points
     if ( xi < xd(1) ) then ! extrapolate
        t  = ( xi - xd(1) ) / ( xd(2) - xd(1) )
        yi = (1._r8 - t) * yd(1) + t * yd(2)
     elseif ( xi > xd(nd) ) then ! extrapolate
        t  = ( xi - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
        yi = (1._r8 - t) * yd(nd-1) + t * yd(nd)
     else  ! piecsewise interpolate
        do k = 2, nd
           if ( (xd(k-1) <= xi) .and. (xi <= xd(k)) ) then
              t  = ( xi - xd(k-1) ) / ( xd(k) - xd(k-1) )
              yi = (1._r8 - t) * yd(k-1) + t * yd(k)
              exit
           endif
        enddo
     endif

     return

   end subroutine piecewise_linear_interp1d
 
 end module SnowSnicarMod
