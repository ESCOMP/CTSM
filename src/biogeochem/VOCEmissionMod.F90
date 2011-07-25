module VOCEmissionMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: VOCEmissionMod
!
! !DESCRIPTION:
! Volatile organic compound emission
!
! !USES:
  use clm_varctl, only: iulog
  use abortutils, only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: VOCEmission
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: VOCEmission
!
! !INTERFACE:
  subroutine VOCEmission (lbp, ubp, num_soilp, filter_soilp )
!
! ! NEW DESCRIPTION
! Volatile organic compound emission
! This code simulates volatile organic compound emissions following:
! 1. Isoprene: Guenther et al., 2006 description of MEGAN emissions
!     following equations 2-9, 16-17, 20
! 2. Monoterpenes/OVOCs/ORVOCs/CO: algorithm presented in Guenther, A., 
!    1999: Modeling Biogenic Volatile Organic Compound Emissions to the 
!    Atmosphere. In Reactive Hydrocarbons in the Atmosphere, Ch. 3
!    With updates from MEGAN online user's guide 
!    ( http://acd.ucar.edu/~guenther/MEGAN/MEGANusersguide.pdf)
! This model relies on the assumption that 90% of isoprene and monoterpene
! emissions originate from canopy foliage:
!    E= epsilon * gamma * rho
! VOC flux (E) [ugC m-2 h-1] is calculated from baseline emission
! factors (epsilon) [ugC m-2 h-1] which are mapped for each PFT (isoprene)
! or constant for each PFT (others).  Note that for constant EFs the units
! of [ugC g-1 h-1] must be multiplied by the source density factor.
! The emission activity factor (gamma) [unitless] for isoprene includes 
! dependence on PPFT, temperature, LAI, leaf age and soil moisture.  
! The canopy environment constant was calculated offline for CLM+CAM at 
! standard conditions.
! The emission activity factor for the other emissions depends on temperature.
! We assume that the escape efficiency (rho) here is unity following
! Guenther et al., 2006.
! Subroutine written to operate at the patch level.
! IN FINAL IMPLEMENTATION, REMEMBER:
! 1. may wish to call this routine only as freq. as rad. calculations
! 2. may wish to place epsilon values directly in pft-physiology file
! Output: vocflx(nvoc) !VOC flux [ug C m-2 h-1]
!
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_atmlnd   , only : clm_a2l
    use clmtype
    use clm_varpar   , only : nvoc, numpft
    use clm_atmlnd   , only : clm_a2l
    use shr_const_mod, only : SHR_CONST_RGAS
    use clm_varcon   , only : denice
    use clm_varpar   , only : nlevsoi
    use pftvarcon    , only : ndllf_evr_tmp_tree,  ndllf_evr_brl_tree,    &
                              ndllf_dcd_brl_tree,  nbrdlf_evr_trp_tree,   &
                              nbrdlf_evr_tmp_tree, nbrdlf_dcd_brl_shrub,  &
                              nbrdlf_dcd_trp_tree, nbrdlf_dcd_tmp_tree,   &
                              nbrdlf_dcd_brl_tree, nbrdlf_evr_shrub,      &
                              nc3_arctic_grass,                           &
                              nc3crop,             nc4_grass,             &
                              noveg
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: num_soilp                   ! number of columns in soil pft filter
    integer, intent(in) :: filter_soilp(num_soilp)     ! pft filter for soil
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
! 2/1/02, Peter Thornton: migration to new data structure
! 4/15/06, Colette L. Heald: modify for updated MEGAN model (Guenther et al., 2006)
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pgridcell(:)     ! gridcell index of corresponding pft
    integer , pointer :: pcolumn(:)       ! column index of corresponding pft
    integer , pointer :: ivt(:)           ! pft vegetation type for current
    real(r8), pointer :: t_veg(:)         ! pft vegetation temperature (Kelvin)
    real(r8), pointer :: fsun(:)          ! sunlit fraction of canopy
    real(r8), pointer :: elai(:)          ! one-sided leaf area index with burying by snow
    real(r8), pointer :: clayfrac(:)      ! fraction of soil that is clay
    real(r8), pointer :: sandfrac(:)      ! fraction of soil that is sand
    real(r8), pointer :: forc_solad(:,:)  ! direct beam radiation (visible only)
    real(r8), pointer :: forc_solai(:,:)  ! diffuse radiation     (visible only)
    real(r8), pointer :: sla(:)           ! specific leaf area [m2 leaf g-1 C]
    real(r8), pointer :: h2osoi_vol(:,:)  ! volumetric soil water (m3/m3)
    real(r8), pointer :: h2osoi_ice(:,:)  ! ice soil content (kg/m3)
    real(r8), pointer :: dz(:,:)          ! depth of layer (m)
    real(r8), pointer :: coszen(:)        ! cosine of solar zenith angle
    real(r8), pointer :: efisop(:,:)      ! emission factors for isoprene for each pft [ug C m-2 h-1]
    real(r8), pointer :: elai_p(:)        ! one-sided leaf area index from previous timestep
    real(r8), pointer :: t_veg24(:)       ! avg pft vegetation temperature for last 24 hrs
    real(r8), pointer :: t_veg240(:)      ! avg pft vegetation temperature for last 240 hrs
    real(r8), pointer :: fsun24(:)        ! sunlit fraction of canopy last 24 hrs
    real(r8), pointer :: fsun240(:)       ! sunlit fraction of canopy last 240 hrs
    real(r8), pointer :: forc_solad24(:)  ! direct beam radiation last 24hrs  (visible only)
    real(r8), pointer :: forc_solai24(:)  ! diffuse radiation  last 24hrs     (visible only)
    real(r8), pointer :: forc_solad240(:) ! direct beam radiation last 240hrs (visible only)
    real(r8), pointer :: forc_solai240(:) ! diffuse radiation  last 240hrs    (visible only)
    real(r8), pointer :: bsw(:,:)         ! Clapp and Hornberger "b" (nlevgrnd)
    real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: sucsat(:,:)      ! minimum soil suction (mm) (nlevgrnd)

    real(r8), parameter :: smpmax = 2.57e5_r8 ! maximum soil matrix potential
!
! local pointers to original implicit out arrays
!
    real(r8), pointer :: vocflx(:,:)      ! VOC flux [ug C m-2 h-1]
    real(r8), pointer :: vocflx_tot(:)    ! VOC flux [ug C m-2 h-1]
    real(r8), pointer :: vocflx_1(:)      ! VOC flux(1) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_2(:)      ! VOC flux(2) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_3(:)      ! VOC flux(3) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_4(:)      ! VOC flux(4) [ug C m-2 h-1]
    real(r8), pointer :: vocflx_5(:)      ! VOC flux(5) [ug C m-2 h-1]
    real(r8), pointer :: Eopt_out(:)     
    real(r8), pointer :: topt_out(:)
    real(r8), pointer :: alpha_out(:)
    real(r8), pointer :: cp_out(:)
    real(r8), pointer :: paru_out(:)
    real(r8), pointer :: par24u_out(:)
    real(r8), pointer :: par240u_out(:)
    real(r8), pointer :: para_out(:)
    real(r8), pointer :: par24a_out(:)
    real(r8), pointer :: par240a_out(:)
    real(r8), pointer :: gamma_out(:)
    real(r8), pointer :: gammaT_out(:)
    real(r8), pointer :: gammaP_out(:)
    real(r8), pointer :: gammaL_out(:)
    real(r8), pointer :: gammaA_out(:)
    real(r8), pointer :: gammaS_out(:)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: fp,p,g,c,n,j            ! indices
    integer  :: ct_bad
    real(r8) :: epsilon(lbp:ubp)        ! emission factor [ugC m-2 h-1]
    real(r8) :: par                     ! temporary
    real(r8) :: par24                   ! temporary
    real(r8) :: par240                  ! temporary
    real(r8) :: density                 ! source density factor [g dry wgt foliar mass/m2 ground]
    real(r8) :: gamma(lbp:ubp)          ! activity factor (accounting for light, T, age, LAI conditions)
    real(r8) :: gamma_p                 ! activity factor for PPFD
    real(r8) :: gamma_l                 ! activity factor for PPFD & LAI
    real(r8) :: gamma_t                 ! activity factor for temperature
    real(r8) :: gamma_a                 ! activity factor for leaf age
    real(r8) :: gamma_sm                ! activity factor for soil moisture
    real(r8) :: x                       ! temporary 
    real(r8) :: Eopt                    ! temporary 
    real(r8) :: topt                    ! temporary 
    real(r8) :: cp                      ! temporary
    real(r8) :: alpha                   ! temporary
    real(r8) :: elai_prev               ! lai for previous timestep
    real(r8) :: fnew, fgro, fmat, fsen  ! fractions of leaves at different phenological stages
    real(r8) :: nl                      ! temporary number of soil levels
    real(r8) :: theta_ice               ! water content in ice in m3/m3
    real(r8) :: wilt                    ! wilting point in m3/m3
    real(r8) :: theta1                  ! temporary
!
! Constants
!
    real(r8), parameter :: R   = SHR_CONST_RGAS*0.001_r8   ! univ. gas constant [J K-1 mol-1]
    real(r8), parameter :: scale_mw =0.882_r8              ! conversion factor for isoprene -> carbon
    real(r8), parameter :: alpha_fix = 0.001_r8            ! empirical coefficient
    real(r8), parameter :: cp_fix = 1.21_r8                ! empirical coefficient
    real(r8), parameter :: ct1 = 95.0_r8                   ! empirical coefficient (70 in User's Guide)
    real(r8), parameter :: ct2 = 230.0_r8                  ! empirical coefficient  (200 in User's Guide)
    real(r8), parameter :: ct3 = 0.00831_r8                ! empirical coefficient (0.0083 in User's Guide)
    real(r8), parameter :: topt_fix = 317._r8              ! std temperature [K]
    real(r8), parameter :: Eopt_fix = 2.26_r8              ! empirical coefficient
    real(r8), parameter :: tstd = 303.15_r8                ! std temperature [K]
    real(r8), parameter :: bet = 0.09_r8                   ! beta empirical coefficient [K-1]
    real(r8), parameter :: clai1 = 0.49_r8                 ! empirical coefficient
    real(r8), parameter :: clai2 = 0.2_r8                  ! empirical coefficient
    real(r8), parameter :: clai3 = 5.0_r8                  ! empirical coefficient
    real(r8), parameter :: Anew = 0.01_r8                  ! relative emission factor for new plants
    real(r8), parameter :: Agro = 0.5_r8                   ! relative emission factor for new plants
    real(r8), parameter :: Amat = 1.0_r8                   ! relative emission factor for new plants
    real(r8), parameter :: Asen = 0.33_r8                  ! relative emission factor for new plants
    real(r8), parameter :: cce = 0.40_r8                   ! factor to set emissions to unity @ std
    real(r8), parameter :: cce1 = 0.47_r8                  ! same as Cce but for non-accumulated vars
    real(r8), parameter :: ca1 = 0.004_r8                  ! empirical coefficent for alpha
    real(r8), parameter :: ca2 = 0.0005_r8                 ! empirical coefficent for alpha
    real(r8), parameter :: ca3 = 0.0468_r8                 ! empirical coefficent for cp
    real(r8), parameter :: par0_sun = 200._r8              ! std conditions for past 24 hrs [umol/m2/s]
    real(r8), parameter :: par0_shade = 50._r8             ! std conditions for past 24 hrs [umol/m2/s]
    real(r8), parameter :: co1 = 313._r8                   ! empirical coefficient
    real(r8), parameter :: co2 = 0.6_r8                    ! empirical coefficient
    real(r8), parameter :: co3 = 2.034_r8                  ! empirical coefficient
    real(r8), parameter :: co4 = 0.05_r8                   ! empirical coefficient
    real(r8), parameter :: tstd0 = 297_r8                  ! std temperature [K]
    real(r8), parameter :: deltheta1=0.06_r8               ! empirical coefficient
    real(r8), parameter :: scaling_to_500_Tg = 5._r8/7._r8 ! J-F Larmaque's empirical scaling factor

!
! These are the values from version of genesis-ibis / 1000.
! CN calculates its own sla [m2 leaf g-1 C]
! Divide by 2 in the equation to get dry weight foliar mass from grams carbon
!
    real(r8) :: hardwire_sla(0:numpft)
    real(r8) :: slarea(lbp:ubp)           ! Specific leaf areas [m2 leaf g-1 C]
    real(r8) :: hardwire_droot(0:numpft)  ! Root depth [m]
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)
    forc_solad => clm_a2l%forc_solad
    forc_solai => clm_a2l%forc_solai
    efisop     => clm3%g%gve%efisop

    ! Assign local pointers to derived subtypes components (column-level)
    h2osoi_vol       => clm3%g%l%c%cws%h2osoi_vol
    h2osoi_ice       => clm3%g%l%c%cws%h2osoi_ice
    dz               => clm3%g%l%c%cps%dz
    bsw              => clm3%g%l%c%cps%bsw
    watsat           => clm3%g%l%c%cps%watsat
    sucsat           => clm3%g%l%c%cps%sucsat

    ! Assign local pointers to derived subtypes components (pft-level)
    pgridcell        => clm3%g%l%c%p%gridcell
    pcolumn          => clm3%g%l%c%p%column
    ivt              => clm3%g%l%c%p%itype
    t_veg            => clm3%g%l%c%p%pes%t_veg
    fsun             => clm3%g%l%c%p%pps%fsun
    elai             => clm3%g%l%c%p%pps%elai
    clayfrac         => clm3%g%l%c%p%pps%clayfrac
    sandfrac         => clm3%g%l%c%p%pps%sandfrac
    vocflx           => clm3%g%l%c%p%pvf%vocflx
    vocflx_tot       => clm3%g%l%c%p%pvf%vocflx_tot
    vocflx_1         => clm3%g%l%c%p%pvf%vocflx_1
    vocflx_2         => clm3%g%l%c%p%pvf%vocflx_2
    vocflx_3         => clm3%g%l%c%p%pvf%vocflx_3
    vocflx_4         => clm3%g%l%c%p%pvf%vocflx_4
    vocflx_5         => clm3%g%l%c%p%pvf%vocflx_5
    Eopt_out         => clm3%g%l%c%p%pvf%Eopt_out
    topt_out         => clm3%g%l%c%p%pvf%topt_out
    alpha_out        => clm3%g%l%c%p%pvf%alpha_out
    cp_out           => clm3%g%l%c%p%pvf%cp_out
    paru_out         => clm3%g%l%c%p%pvf%paru_out
    par24u_out       => clm3%g%l%c%p%pvf%par24u_out
    par240u_out      => clm3%g%l%c%p%pvf%par240u_out
    para_out         => clm3%g%l%c%p%pvf%para_out
    par24a_out       => clm3%g%l%c%p%pvf%par24a_out
    par240a_out      => clm3%g%l%c%p%pvf%par240a_out
    gammaL_out       => clm3%g%l%c%p%pvf%gammaL_out
    gammaT_out       => clm3%g%l%c%p%pvf%gammaT_out
    gammaP_out       => clm3%g%l%c%p%pvf%gammaP_out
    gammaA_out       => clm3%g%l%c%p%pvf%gammaA_out
    gammaS_out       => clm3%g%l%c%p%pvf%gammaS_out
    gamma_out        => clm3%g%l%c%p%pvf%gamma_out
    sla              => clm3%g%l%c%p%pps%slasha

    t_veg24          => clm3%g%l%c%p%pvs%t_veg24
    t_veg240         => clm3%g%l%c%p%pvs%t_veg240
    forc_solad24     => clm3%g%l%c%p%pvs%fsd24
    forc_solad240    => clm3%g%l%c%p%pvs%fsd240
    forc_solai24     => clm3%g%l%c%p%pvs%fsi24
    forc_solai240    => clm3%g%l%c%p%pvs%fsi240
    fsun24           => clm3%g%l%c%p%pvs%fsun24
    fsun240          => clm3%g%l%c%p%pvs%fsun240
    elai_p           => clm3%g%l%c%p%pvs%elai_p

    hardwire_sla(noveg)                                    = 0._r8     ! bare-soil

    hardwire_sla(ndllf_evr_tmp_tree)                       = 0.0125_r8 !needleleaf
    hardwire_sla(ndllf_evr_brl_tree)                       = 0.0125_r8 !Gordon Bonan suggests NET = 0.0076
    hardwire_sla(ndllf_dcd_brl_tree)                       = 0.0125_r8 !Gordon Bonan suggests NDT = 0.0200

    hardwire_sla(nbrdlf_evr_trp_tree)                      = 0.0250_r8 !broadleaf
    hardwire_sla(nbrdlf_evr_tmp_tree)                      = 0.0250_r8 !Gordon Bonan suggests BET = 0.0178
    hardwire_sla(nbrdlf_dcd_trp_tree)                      = 0.0250_r8 !Gordon Bonan suggests BDT = 0.0274
    hardwire_sla(nbrdlf_dcd_tmp_tree:nbrdlf_dcd_brl_shrub) = 0.0250_r8 

    hardwire_sla(nc3_arctic_grass:numpft)                  = 0.0200_r8 !grass/crop

! root depth (m) (defined based on Zeng et al., 2001, cf Guenther 2006)

    hardwire_droot(noveg)                                     = 0._r8   ! bare-soil
    hardwire_droot(ndllf_evr_tmp_tree:ndllf_evr_brl_tree)     = 1.8_r8  ! evergreen tree
    hardwire_droot(ndllf_dcd_brl_tree)                        = 2.0_r8  ! needleleaf deciduous boreal tree
    hardwire_droot(nbrdlf_evr_trp_tree:nbrdlf_evr_tmp_tree)   = 3.0_r8  ! broadleaf evergreen tree
    hardwire_droot(nbrdlf_dcd_trp_tree:nbrdlf_dcd_brl_tree)   = 2.0_r8  ! broadleaf deciduous tree
    hardwire_droot(nbrdlf_evr_shrub:nbrdlf_dcd_brl_shrub)     = 2.5_r8  ! shrub
    hardwire_droot(nc3_arctic_grass:numpft)                   = 1.5_r8  ! grass/crop

! initialize variables which get passed to the atmosphere
    vocflx(lbp:ubp, :)=0._r8

    ! Determine specific leaf array
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       slarea(p) = hardwire_sla(ivt(p))

    end do


    ! Begin loop through voc species
    !_______________________________________________________________________________

    do n = 1, nvoc
       select case (n)

       case(1)	

          do fp = 1,num_soilp
             p = filter_soilp(fp)
             g = pgridcell(p)


             ! epsilon: use gridded values for 6 PFTs specified by MEGAN following
             ! -------  Guenther et al. (2006).  Map the numpft CLM PFTs to these 6.
             !          Units: [ug C m-2 h-1] (convert input files from units of 
             !                 [ug isop m-2 h-1])
    	     epsilon(p) = 0._r8

             ! isoprenes:
             if (     ivt(p) == ndllf_evr_tmp_tree  &
             .or.     ivt(p) == ndllf_evr_brl_tree) then     !fineleaf evergreen
                	epsilon(p) = efisop(2,g)*scale_mw
             else if (ivt(p) == ndllf_dcd_brl_tree) then     !fineleaf deciduous
                	epsilon(p) = efisop(3,g)*scale_mw
             else if (ivt(p) >= nbrdlf_evr_trp_tree &
             .and.    ivt(p) <= nbrdlf_dcd_brl_tree) then    !broadleaf trees
                	epsilon(p) = efisop(1,g)*scale_mw
             else if (ivt(p) >= nbrdlf_evr_shrub &
             .and.    ivt(p) <= nbrdlf_dcd_brl_shrub) then   !shrubs
                	epsilon(p) = efisop(4,g)*scale_mw
             else if (ivt(p) >= nc3_arctic_grass &
             .and.    ivt(p) <= nc4_grass) then              !grass
                	epsilon(p) = efisop(5,g)*scale_mw
             else if (ivt(p) >= nc3crop) then                !crops
                	epsilon(p) =efisop(6,g)*scale_mw
             end if

          end do

       case(2)

          do fp = 1,num_soilp
             p = filter_soilp(fp)
             g = pgridcell(p)

             ! epsilon: use values from table 3 in Guenther (1997) which originate in
             ! -------  Guenther et al. (1995). In the comments below, I mention the pft
             !          category as described in table 3. Some values were taken directly
             !          from Guenther et al. (1995). Units: [ugC g-1 h-1]
             !          Values were updated on 1/2002 (Guenther, personal communication)

             ! monoterpenes:
             epsilon(p) = 0._r8
             ! monoterpenes:
             if (     ivt(p) >= ndllf_evr_tmp_tree &
             .and.    ivt(p) <= ndllf_evr_brl_tree) then     !needleleaf evergreen
                epsilon(p) = 2.0_r8
             else if (ivt(p) == ndllf_dcd_brl_tree) then     !needleleaf deciduous
                epsilon(p) = 1.6_r8
             else if (ivt(p) >= nbrdlf_evr_trp_tree  &
             .and.    ivt(p) <= nbrdlf_dcd_brl_tree) then    !broadleaf everg trop
                epsilon(p) = 0.4_r8
             else if (ivt(p) >= nbrdlf_evr_shrub &
             .and.    ivt(p) <= nbrdlf_dcd_brl_shrub) then   !other woody veg
                epsilon(p) = 0.8_r8
             else if (ivt(p) >= nc3_arctic_grass &
             .and.    ivt(p) <= numpft) then                 !grass & crop
                epsilon(p) = 0.1_r8
             end if
          end do

       case (3)
          do fp = 1,num_soilp
             p = filter_soilp(fp)
             g = pgridcell(p)

             ! other VOCs (OVOCs)
             epsilon(p) = 1.0_r8                 !Guenther (personal communication)
          end do

       case (4)
          do fp = 1,num_soilp
             p = filter_soilp(fp)
             g = pgridcell(p)

             ! other reactive VOCs (ORVOCs)
             epsilon(p) = 1.0_r8                 !Guenther (personal communication)
          end do

       case (5)
          do fp = 1,num_soilp
             p = filter_soilp(fp)
             g = pgridcell(p)

             ! CO
             epsilon(p) = 0.3_r8                 !Guenther (personal communication)
          end do


       case default

          write(iulog,*)'only nvocs up to index 5 are currently supported'
          call endrun()

       end select
       
       
       ct_bad=0

       select case (n)

       case (1)

          do fp = 1,num_soilp
             p = filter_soilp(fp)
             g = pgridcell(p)
             c = pcolumn(p)


             ! gamma: Activity factor. Units [dimensionless]
             ! =====  For isoprene include activity factors for LAI,PPFD, T, leaf age, and soil moisture

             ! Activity factor for LAI (Guenther et al., 2006)
             !------------------------
             ! Guenther et al., 2006 eq 3
             if ( (fsun240(p) > 0.0_r8) .and. (fsun240(p) < 1.e30_r8) ) then 
                 gamma_l = cce * elai(p)
             else
                 gamma_l = cce1 * elai(p)
             end if
	     gammaL_out(p)=gamma_l

             ! Activity factor for PPFD (Guenther et al., 2006)
             !-------------------------
	     ! With distinction between sunlit and shaded leafs, weight scalings by
             ! fsun and fshade 
             ! Scale total incident par by fraction of sunlit leaves (added on 1/2002)
             ! multiply w/m2 by 4.6 to get umol/m2/s for par (added 8/14/02)

             ! fvitt -- forc_solad240, forc_solai240 can be zero when CLM finidat is specified
             !          which will cause par240 to be zero and produce NaNs via log(par240)
             ! dml   -- fsun240 can be equal to or greater than one before 10 day averages are
             !           set on startup or if a new pft comes online during land cover change.
             !           Avoid this problem by only doing calculations with fsun240 when fsun240 is
             !           between 0 and 1
             if ( (fsun240(p) > 0._r8) .and. (fsun240(p) < 1._r8) .and.  (forc_solad240(p) > 0._r8) &
             .and. (forc_solai240(p) > 0._r8)) then
                ! With alpha and cp calculated based on eq 6 and 7:
                ! Note indexing for accumulated variables is all at pft level
                ! SUN:
                par = (forc_solad(g,1) + fsun(p) * forc_solai(g,1)) * 4.6_r8
                par24 = (forc_solad24(p) + fsun24(p) * forc_solai24(p)) * 4.6_r8
                par240 = (forc_solad240(p) + fsun240(p) * forc_solai240(p)) * 4.6_r8
                alpha = ca1 - ca2 * log(par240)
                cp = ca3 * exp(ca2 * (par24-par0_sun))*par240**(0.6_r8)
                gamma_p = fsun(p) * ( cp * alpha*par * (1._r8 + alpha*alpha*par*par)**(-0.5_r8) )
	        paru_out(p)=par
		par24u_out(p)=par24
                par240u_out(p)=par240
                ! SHADE:
                par = ((1._r8 - fsun(p)) * forc_solai(g,1)) * 4.6_r8
                par24 = ((1._r8 - fsun24(p)) * forc_solai24(p)) * 4.6_r8
                par240 = ((1._r8 - fsun240(p)) * forc_solai240(p)) * 4.6_r8
                alpha = ca1 - ca2 * log(par240)
                cp = ca3 * exp(ca2 * (par24-par0_shade))*par240**(0.6_r8)
                par = ((1._r8 - fsun(p)) * forc_solai(g,1)) * 4.6_r8
                gamma_p = gamma_p + (1-fsun(p)) * (cp*alpha*par*(1._r8 + alpha*alpha*par*par)**(-0.5_r8))
                para_out(p)=par
		par24a_out(p)=par24
 		par240a_out(p)=par240
             else
                ! With fixed alpha and cp (from MEGAN User's Guide):
                ! SUN: direct + diffuse  
                par = (forc_solad(g,1) + fsun(p) * forc_solai(g,1)) * 4.6_r8
                alpha = alpha_fix
                cp = cp_fix
                gamma_p = fsun(p) * ( cp * alpha*par * (1._r8 + alpha*alpha*par*par)**(-0.5_r8) )
		paru_out(p)=par
	        par24u_out(p)=-999
	        par240u_out(p)=-999
                ! SHADE: diffuse 
                par = ((1._r8 - fsun(p)) * forc_solai(g,1)) * 4.6_r8
                gamma_p = gamma_p + (1-fsun(p)) * (cp*alpha*par*(1._r8 + alpha*alpha*par*par)**(-0.5_r8))
		para_out(p)=par
                par24a_out(p)=-999
                par240a_out(p)=-999
             end if 
             alpha_out(p)=alpha
             cp_out(p)=cp
             gammaP_out(p)=gamma_p


             ! Activity factor for temperature (Guenther et al., 2006)
             !--------------------------------
             if ( (t_veg240(p) > 0.0_r8) .and. (t_veg240(p) < 1.e30_r8) ) then 
                ! topt and Eopt from eq 8 and 9:
                topt = co1 + (co2 * (t_veg240(p)-tstd0))
                Eopt = co3 * exp (co4 * (t_veg24(p)-tstd0)) * exp(co4 * (t_veg240(p) -tstd0))
	     else
                topt = topt_fix
                Eopt = Eopt_fix
             endif 
             x = ( (1._r8/topt) - (1._r8/(t_veg(p))) ) / ct3
             gamma_t = Eopt * ( ct2 * exp(ct1 * x)/(ct2 - ct1 * (1._r8 - exp(ct2 * x))) )
             topt_out(p)=topt
             Eopt_out(p)=Eopt
             gammaT_out(p)=gamma_t


             ! Activity factor for leaf age (Guenther et al., 2006)
             !-----------------------------
             ! If not CNDV elai is constant therefore gamma_a=1.0
             ! gamma_a set to unity for evergreens (PFTs 1, 2, 4, 5)
             ! Note that we assume here that the time step is shorter than the number of 
             !days after budbreak required to induce isoprene emissions (ti=12 days) and 
             ! the number of days after budbreak to reach peak emission (tm=28 days)
	     if ( (ivt(p) == ndllf_dcd_brl_tree) .or. (ivt(p) >= nbrdlf_dcd_trp_tree) ) then  ! non-evergreen

                if ( (elai_p(p) > 0.0_r8) .and. (elai_p(p) < 1.e30_r8) )then 
                   elai_prev = 2._r8*elai_p(p)-elai(p)  ! have accumulated average lai over last timestep
                   if (elai_prev == elai(p)) then
                      fnew = 0.0_r8
                      fgro = 0.0_r8
                      fmat = 1.0_r8
                      fsen = 0.0_r8
                   else if (elai_prev > elai(p)) then
                      fnew = 0.0_r8
                      fgro = 0.0_r8
                      fmat = 1.0_r8 - (elai_prev - elai(p))/elai_prev
                      fsen = (elai_prev - elai(p))/elai_prev
                   else if (elai_prev < elai(p)) then
                      fnew = 1 - (elai_prev / elai(p))
                      fgro = 0.0_r8
                      fmat = (elai_prev / elai(p))
                      fsen = 0.0_r8
                   end if             
                
                   gamma_a = fnew * Anew + fgro * Agro + fmat * Amat + fsen * Asen
	        else
                   gamma_a = 1.0_r8
                end if

             else
                gamma_a = 1.0_r8
             end if
             gammaA_out(p)=gamma_a


             ! Activity factor for soil moisture (Guenther et al., 2006) 
             !----------------------------------
             ! Calculate the mean scaling factor throughout the root depth.
             ! wilting point potential is in units of matric potential (mm) 
             ! (1 J/Kg = 0.001 MPa, approx = 0.1 m)
             ! convert to volumetric soil water using equation 7.118 of the CLM4 Technical Note
             if ((clayfrac(p) > 0) .and. (sandfrac(p) > 0)) then 
               gamma_sm = 0._r8
	       nl=0._r8

               do j = 1,nlevsoi
	         if  (sum(dz(c,1:j)) < hardwire_droot(ivt(p)))  then
                   theta_ice = h2osoi_ice(c,j)/(dz(c,j)*denice)
                   wilt = ((smpmax/sucsat(c,j))**(-1._r8/bsw(c,j))) * (watsat(c,j) - theta_ice)
                   theta1 = wilt + deltheta1
                   if (h2osoi_vol(c,j) >= theta1) then 
             	      gamma_sm = gamma_sm + 1._r8
                   else if ( (h2osoi_vol(c,j) > wilt) .and. (h2osoi_vol(c,j) < theta1) ) then
		      gamma_sm = gamma_sm + ( h2osoi_vol(c,j) - wilt ) / deltheta1
                   else
		      gamma_sm = gamma_sm + 0._r8
                   end if
		   nl=nl+1._r8
                 end if
 	       end do 

	       if (nl > 0) then
	         gamma_sm = gamma_sm/nl
	       endif
             else
	       gamma_sm = 1.0_r8
             end if
             gammaS_out(p)=gamma_sm


             ! Calculate total scaling factor
             !--------------------------------
	     gamma(p) = gamma_l * gamma_p * gamma_t * gamma_a * gamma_sm
             if ( (gamma(p) >=0.0_r8) .and. (gamma(p)< 100._r8) ) then
                gamma_out(p)=gamma(p)
	     else
                gamma_out(p)=gamma(p)
                write(iulog,*) 'clh GAMMA: ',gamma(p),gamma_l,gamma_p,gamma_t,gamma_a,gamma_sm
             end if

          end do

       case (2,3,4,5)

          do fp = 1,num_soilp
             p = filter_soilp(fp)
             g = pgridcell(p)

             ! gamma: Activity factor. Units [dimensionless]
             ! -----  For monoterpenes, OVOCs, ORVOCs, CO include simple activity factors 
             !        for LAI and T only (Guenther et al., 1995)
             gamma_t = exp(bet * (t_veg(p) - tstd))
	     gamma(p)=gamma_t

          end do

       end select

       do fp = 1,num_soilp
          p = filter_soilp(fp)
          g = pgridcell(p)

          ! density: Source density factor [g dry weight foliar mass m-2 ground]
          ! -------  Other than isoprene, need to convert EF units from 
          ! [ug g-1 h-1] to [ug m-2 h-1]
          if (ivt(p) > noveg) then
             density = elai(p) / (slarea(p) * 0.5_r8)
          else
             density = 0._r8
          end if

          ! calculate the voc flux
          ! ----------------------
	  select case (n)

          case(1)

              vocflx(p,n) = epsilon(p) * gamma(p) * scaling_to_500_Tg

          case(2,3,4,5)

              vocflx(p,n) = epsilon(p) * gamma(p) * density

          end select


       end do   ! end pft loop

    end do   ! end voc species loop
    !_______________________________________________________________________________

    ! Calculate total voc flux and individual components for history output

    do fp = 1,num_soilp
       p = filter_soilp(fp)
       vocflx_tot(p) = 0._r8
    end do
    do n = 1, nvoc
       do fp = 1,num_soilp
          p = filter_soilp(fp)
          vocflx_tot(p) = vocflx_tot(p) + vocflx(p,n)
       end do
    end do
    do fp = 1,num_soilp
       p = filter_soilp(fp)
       g = pgridcell(p)
       vocflx_1(p) = vocflx(p,1)
       vocflx_2(p) = vocflx(p,2)
       vocflx_3(p) = vocflx(p,3)
       vocflx_4(p) = vocflx(p,4)
       vocflx_5(p) = vocflx(p,5)
    end do

  end subroutine VOCEmission

end module VOCEmissionMod
