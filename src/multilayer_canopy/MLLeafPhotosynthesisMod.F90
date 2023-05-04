module MLLeafPhotosynthesisMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate leaf photosynthesis and stomatal conductance
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LeafPhotosynthesis       ! Leaf photosynthesis and stomatal conductance
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: ft                      ! Photosynthesis temperature response
  private :: fth                     ! Photosynthesis temperature inhibition
  private :: fth25                   ! Scaling factor for photosynthesis temperature inhibition
  private :: CiFunc                  ! Calculate An and gs for a specified Ci
  private :: GsFunc                  ! Calculate An for a specified gs
  private :: StomataOptimization     ! Photosynthesis and stomatal conductance with optimization
  private :: StomataEfficiency       ! Water-use efficiency and cavitation checks for optimal gs
  private :: StomataFluxes           ! Leaf calculations for a specified stomatal conductance
  private :: LeafWaterPotential      ! Leaf water potential for transpiration rate
  private :: LeafTranspiration       ! Leaf transpiration flux
  private :: C13Fractionation        ! 13C fractionation for photosynthesis
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  function ft (tl, ha) result(ans)
    !
    ! !DESCRIPTION:
    ! Photosynthesis temperature response
    !
    ! !USES:
    use clm_varcon, only : tfrz
    use MLclm_varcon, only : rgas
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: tl    ! Leaf temperature (K)
    real(r8), intent(in) :: ha    ! Activation energy (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans               ! Temperature function value
    !---------------------------------------------------------------------

    ans = exp( ha / (rgas * (tfrz+25._r8)) * (1._r8 - (tfrz+25._r8) / tl) )

  end function ft

  !-----------------------------------------------------------------------
  function fth (tl, hd, se, c) result(ans)
    !
    ! !DESCRIPTION:
    ! Photosynthesis temperature inhibition
    !
    ! !USES:
    use MLclm_varcon, only : rgas
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: tl    ! Leaf temperature (K)
    real(r8), intent(in) :: hd    ! Deactivation energy (J/mol)
    real(r8), intent(in) :: se    ! Entropy term (J/mol/K)
    real(r8), intent(in) :: c     ! Scaling factor for high temperature inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans               ! Temperature function value
    !---------------------------------------------------------------------

    ans = c / ( 1._r8 + exp( (-hd + se*tl) / (rgas*tl) ) )

  end function fth

  !-----------------------------------------------------------------------
  function fth25 (hd, se) result(ans)
    !
    ! !DESCRIPTION:
    ! Scaling factor for photosynthesis temperature inhibition
    !
    ! !USES:
    use clm_varcon, only : tfrz
    use MLclm_varcon, only : rgas
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: hd    ! Deactivation energy (J/mol)
    real(r8), intent(in) :: se    ! Entropy term (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans               ! Temperature function value
    !---------------------------------------------------------------------

    ans = 1._r8 + exp( (-hd + se * (tfrz+25._r8)) / (rgas * (tfrz+25._r8)) )

  end function fth25

  !-----------------------------------------------------------------------
  subroutine LeafPhotosynthesis (num_filter, filter, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance
    !
    ! !USES:
    use clm_varcon, only : tfrz
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only: kc25, ko25, cp25, kcha, koha, cpha
    use MLclm_varcon, only: vcmaxha_noacclim, vcmaxha_acclim, jmaxha_noacclim, jmaxha_acclim
    use MLclm_varcon, only: vcmaxhd_noacclim, vcmaxhd_acclim, jmaxhd_noacclim, jmaxhd_acclim
    use MLclm_varcon, only: vcmaxse_noacclim, vcmaxse_acclim, jmaxse_noacclim, jmaxse_acclim
    use MLclm_varcon, only: rdha, rdhd, rdse
    use MLclm_varcon, only: phi_psII, theta_j, vpd_min_MED, rh_min_BB
    use MLclm_varctl, only : gs_type, acclim_type
    use MLMathToolsMod, only : hybrid, quadratic
    use MLWaterVaporMod, only : SatVap
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    integer, intent(in) :: il           ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                      ! Aboveground layer index
    real(r8) :: vcmaxha                 ! Activation energy for vcmax (J/mol)
    real(r8) :: jmaxha                  ! Activation energy for jmax (J/mol)
    real(r8) :: vcmaxhd                 ! Deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd                  ! Deactivation energy for jmax (J/mol)
    real(r8) :: vcmaxse                 ! Entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse                  ! Entropy term for jmax (J/mol/K)
    real(r8) :: vcmaxc                  ! Scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc                   ! Scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: rdc                     ! Scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: qabs                    ! PAR utilized by PS II (umol photons/m2/s)
    real(r8) :: desat                   ! Derivative of saturation vapor pressure (Pa/K)
    real(r8) :: gs_err                  ! gs for error check
    real(r8) :: an_err                  ! An for error check
    real(r8) :: aquad,bquad,cquad       ! Terms for quadratic equations
    real(r8) :: r1,r2                   ! Roots of quadratic equation
    real(r8) :: ci0, ci1                ! Initial estimates for Ci
    real(r8) :: hs_term                 ! Leaf surface humidity term (-)
    real(r8) :: vpd_term                ! Leaf vapor pressure deficit term (kPa)
    real(r8) :: t1,t2,t3,t4             ! C4 temperature terms
    real(r8), parameter :: tol = 0.1_r8 ! Convergence tolerance for Ci (mol/mol)
    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    c3psn     => pftcon%c3psn                  , &  ! CLM: Photosynthetic pathway (1. = C3 plant, 0. = C4 plant)
    g0opt_BB  => pftcon%g0opt_BB               , &  ! CLM (new): Ball-Berry minimum leaf conductance, unstressed (mol H2O/m2/s)
    g1opt_BB  => pftcon%g1opt_BB               , &  ! CLM (new): Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    g0opt_MED => pftcon%g0opt_MED              , &  ! CLM (new): Medlyn minimum leaf conductance, unstressed (mol H2O/m2/s)
    g1opt_MED => pftcon%g1opt_MED              , &  ! CLM (new): Medlyn slope of conductance-photosynthesis relationship, unstressed
    tacclim   => mlcanopy_inst%tacclim_forcing , &  ! Average air temperature for acclimation (K)
    ncan      => mlcanopy_inst%ncan_canopy     , &  ! Number of layers
    btran     => mlcanopy_inst%btran_soil      , &  ! Soil wetness factor for stomatal conductance (-)
    dpai      => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    vcmax25   => mlcanopy_inst%vcmax25_profile , &  ! Canopy layer leaf maximum carboxylation rate at 25C (umol/m2/s)
    jmax25    => mlcanopy_inst%jmax25_profile  , &  ! Canopy layer C3 maximum electron transport rate at 25C (umol/m2/s)
    kp25      => mlcanopy_inst%kp25_profile    , &  ! Canopy layer C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    rd25      => mlcanopy_inst%rd25_profile    , &  ! Canopy layer leaf respiration rate at 25C (umol CO2/m2/s)
    eair      => mlcanopy_inst%eair_profile    , &  ! Canopy layer vapor pressure (Pa)
    cair      => mlcanopy_inst%cair_profile    , &  ! Canopy layer atmospheric CO2 (umol/mol)
    tleaf     => mlcanopy_inst%tleaf_leaf      , &  ! Leaf temperature (K)
    gbv       => mlcanopy_inst%gbv_leaf        , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gbc       => mlcanopy_inst%gbc_leaf        , &  ! Leaf boundary layer conductance: CO2 (mol CO2/m2 leaf/s)
    apar      => mlcanopy_inst%apar_leaf       , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
                                                    ! *** Output ***
    g0        => mlcanopy_inst%g0_canopy       , &  ! Ball-Berry or Medlyn minimum leaf conductance (mol H2O/m2/s)
    g1        => mlcanopy_inst%g1_canopy       , &  ! Ball-Berry or Medlyn slope parameter
    kc        => mlcanopy_inst%kc_leaf         , &  ! Leaf Michaelis-Menten constant for CO2 (umol/mol)
    ko        => mlcanopy_inst%ko_leaf         , &  ! Leaf Michaelis-Menten constant for O2 (mmol/mol)
    cp        => mlcanopy_inst%cp_leaf         , &  ! Leaf CO2 compensation point (umol/mol)
    vcmax     => mlcanopy_inst%vcmax_leaf      , &  ! Leaf maximum carboxylation rate (umol/m2/s)
    jmax      => mlcanopy_inst%jmax_leaf       , &  ! Leaf maximum electron transport rate (umol/m2/s)
    je        => mlcanopy_inst%je_leaf         , &  ! Leaf electron transport rate (umol/m2/s)
    kp        => mlcanopy_inst%kp_leaf         , &  ! Leaf C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    rd        => mlcanopy_inst%rd_leaf         , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    ci        => mlcanopy_inst%ci_leaf         , &  ! Leaf intercellular CO2 (umol/mol)
    hs        => mlcanopy_inst%hs_leaf         , &  ! Leaf fractional humidity at leaf surface (dimensionless)
    vpd       => mlcanopy_inst%vpd_leaf        , &  ! Leaf vapor pressure deficit (Pa)
    ceair     => mlcanopy_inst%ceair_leaf      , &  ! Leaf vapor pressure of air, constrained for stomatal conductance (Pa)
    leaf_esat => mlcanopy_inst%leaf_esat_leaf  , &  ! Leaf saturation vapor pressure (Pa)
    lwpleaf   => mlcanopy_inst%lwpleaf_leaf    , &  ! Leaf water potential (MPa)
                                                    ! *** Output from calls to CiFunc or StomataOptimization ***
    ac        => mlcanopy_inst%ac_leaf         , &  ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj        => mlcanopy_inst%aj_leaf         , &  ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ap        => mlcanopy_inst%ap_leaf         , &  ! Leaf product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    ag        => mlcanopy_inst%ag_leaf         , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an        => mlcanopy_inst%an_leaf         , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs        => mlcanopy_inst%cs_leaf         , &  ! Leaf surface CO2 (umol/mol)
    gs        => mlcanopy_inst%gs_leaf         , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
                                                    ! *** Output from call to C13Fractionation ***
    alphapsn  => mlcanopy_inst%alphapsn_leaf     &  ! Leaf 13C fractionation factor for photosynthesis (-)
    )

    do fp = 1, num_filter
       p = filter(fp)
       do ic = 1, ncan(p)

          ! Adjustments for temperature acclimation

          select case (acclim_type)
          case (0)
             ! No temperature acclimation
             vcmaxha = vcmaxha_noacclim
             jmaxha  = jmaxha_noacclim
             vcmaxhd = vcmaxhd_noacclim
             jmaxhd  = jmaxhd_noacclim
             vcmaxse = vcmaxse_noacclim
             jmaxse  = jmaxse_noacclim
          case (1)
             ! With temperature acclimation
             vcmaxha = vcmaxha_acclim
             jmaxha  = jmaxha_acclim
             vcmaxhd = vcmaxhd_acclim
             jmaxhd  = jmaxhd_acclim
             vcmaxse_acclim = 668.39_r8 - 1.07_r8 * min(max((tacclim(p)-tfrz),11._r8),35._r8)
             jmaxse_acclim  = 659.70_r8 - 0.75_r8 * min(max((tacclim(p)-tfrz),11._r8),35._r8)
             vcmaxse = vcmaxse_acclim
             jmaxse  = jmaxse_acclim
          case default
             call endrun (msg=' ERROR: LeafPhotosynthesis: acclim_type not valid')
          end select

          ! High temperature deactivation: 
          ! The factor "c" scales the deactivation to a value of 1.0 at 25C

          vcmaxc = fth25 (vcmaxhd, vcmaxse)
          jmaxc  = fth25 (jmaxhd, jmaxse)
          rdc    = fth25 (rdhd, rdse)

          ! Process each canopy layer

          if (dpai(p,ic) > 0._r8) then

             ! C3 photosynthetic temperature response

             kc(p,ic,il)    = kc25          * ft(tleaf(p,ic,il), kcha)
             ko(p,ic,il)    = ko25          * ft(tleaf(p,ic,il), koha)
             cp(p,ic,il)    = cp25          * ft(tleaf(p,ic,il), cpha)
             vcmax(p,ic,il) = vcmax25(p,ic) * ft(tleaf(p,ic,il), vcmaxha) * fth(tleaf(p,ic,il), vcmaxhd, vcmaxse, vcmaxc) 
             jmax(p,ic,il)  = jmax25(p,ic)  * ft(tleaf(p,ic,il), jmaxha)  * fth(tleaf(p,ic,il), jmaxhd, jmaxse, jmaxc)
             rd(p,ic,il)    = rd25(p,ic)    * ft(tleaf(p,ic,il), rdha)    * fth(tleaf(p,ic,il), rdhd, rdse, rdc)
             kp(p,ic,il)    = 0._r8

             ! C4 photosynthetic temperature response

             if (nint(c3psn(patch%itype(p))) == 0) then
                t1 = 2.0**( (tleaf(p,ic,il)-(tfrz+25._r8)) / 10._r8 ) 
                t2 = 1._r8 + exp(0.2_r8*((tfrz+15._r8)-tleaf(p,ic,il))) 
                t3 = 1._r8 + exp(0.3_r8*(tleaf(p,ic,il)-(tfrz+40._r8)))
                t4 = 1._r8 + exp(1.3_r8*(tleaf(p,ic,il)-(tfrz+55._r8)))
                vcmax(p,ic,il)  = vcmax25(p,ic) * t1 / (t2 * t3)
                rd(p,ic,il) = rd25(p,ic) * t1 / t4
                kp(p,ic,il) = kp25(p,ic) * t1
             end if

             ! Soil water effect

             select case (gs_type)
             case (0)
                ! Medlyn conductance
                vcmax(p,ic,il) = vcmax(p,ic,il) * btran(p)
                g0(p) = g0opt_MED(patch%itype(p))
                g1(p) = g1opt_MED(patch%itype(p))
             case (1)
                ! Ball-Berry conductance
                vcmax(p,ic,il) = vcmax(p,ic,il) * btran(p)
                g0(p) = max( g0opt_BB(patch%itype(p)) * btran(p), 1.e-06_r8 )
                g1(p) = g1opt_BB(patch%itype(p))
             end select

             ! Saturation vapor pressure at leaf temperature

             call SatVap (tleaf(p,ic,il), leaf_esat(p,ic,il), desat)

             ! Constrain canopy air vapor pressure to <= leaf_esat

             ceair(p,ic,il) = min(eair(p,ic), leaf_esat(p,ic,il))

             ! Constrain eair >= rh_min_BB * leaf_esat so that Ball-Berry solution
             ! does not blow up at low relative humidity. This ensures that hs does
             ! not go to zero. 

             select case (gs_type)
             case (1)
                ceair(p,ic,il) = max(ceair(p,ic,il), rh_min_BB*leaf_esat(p,ic,il))
             end select

             ! Electron transport rate for C3 plants

             qabs = 0.5_r8 * phi_psII * apar(p,ic,il)
             aquad = theta_j
             bquad = -(qabs + jmax(p,ic,il))
             cquad = qabs * jmax(p,ic,il)
             call quadratic (aquad, bquad, cquad, r1, r2)
             je(p,ic,il) = min(r1,r2)

             ! Calculate photosynthesis and stomatal conductance

             select case (gs_type)
             case (0, 1)

                ! Ball-Berry or Medlyn conductance

                if (nint(c3psn(patch%itype(p))) == 1) then
                   ci0 = 0.7_r8 * cair(p,ic)
                else if (nint(c3psn(patch%itype(p))) == 0) then
                   ci0 = 0.4_r8 * cair(p,ic)
                end if
                ci1 = ci0 * 0.99_r8

                ! Solve for Ci: Use CiFunc to iterate photosynthesis calculations
                ! until the change in Ci is < tol. Ci has units umol/mol

                ci(p,ic,il) = hybrid ('LeafPhotosynthesis', p, ic, il, mlcanopy_inst, CiFunc, ci0, ci1, tol)

             case (2)

                ! Use water-use efficiency optimization and cavitation check

                call StomataOptimization (p, ic, il, mlcanopy_inst)

             case default
                call endrun (msg=' ERROR: LeafPhotosynthesis: gs_type not valid')

             end select

             ! Relative humidity and vapor pressure at leaf surface

             hs(p,ic,il) = (gbv(p,ic,il)*eair(p,ic) + gs(p,ic,il)*leaf_esat(p,ic,il)) &
                         / ((gbv(p,ic,il) + gs(p,ic,il)) * leaf_esat(p,ic,il))
             vpd(p,ic,il) = max(leaf_esat(p,ic,il) - hs(p,ic,il)*leaf_esat(p,ic,il), 0.1_r8)

             ! Error checks

             if (gs(p,ic,il) < 0._r8) then
                call endrun (msg=' ERROR: LeafPhotosynthesis: negative stomatal conductance')
             end if

             ! Check solution for Ball-Berry or Medlyn conductance

             hs_term = (gbv(p,ic,il)*ceair(p,ic,il) + gs(p,ic,il)*leaf_esat(p,ic,il)) &
                     / ((gbv(p,ic,il) + gs(p,ic,il)) * leaf_esat(p,ic,il))
             vpd_term = (leaf_esat(p,ic,il) - hs_term * leaf_esat(p,ic,il)) * 0.001_r8

             select case (gs_type)
             case (1)

                ! Compare with Ball-Berry model: gs = g0 + g1 * An * hs/cs

                gs_err = g0(p) + g1(p) * max(an(p,ic,il), 0._r8) * hs_term / cs(p,ic,il)
                if (abs(gs(p,ic,il)-gs_err) > 1.e-06_r8) then
                   call endrun (msg=' ERROR: LeafPhotosynthesis: failed Ball-Berry error check')
                end if

             case (0)

                ! Compare with Medlyn model: gs = g0 + 1.6 * (1 + g1 / sqrt(Ds)) * An / cs
                ! Remember that the vpd term is constrained to >= vpd_min_MED to prevent
                ! infinite gs at low vpd

                if ((leaf_esat(p,ic,il) - ceair(p,ic,il)) > vpd_min_MED) then
                   gs_err = g0(p) + 1.6_r8 * (1._r8 + g1(p) / sqrt(vpd_term)) * max(an(p,ic,il),0._r8) / cs(p,ic,il)
                   if (abs(gs(p,ic,il)-gs_err) > 1.e-06_r8) then
                      call endrun (msg=' ERROR: LeafPhotosynthesis: failed Medlyn error check')
                   end if
                end if

             end select

             ! Compare with diffusion equation: An = (ca - ci) * gleaf

             an_err = (cair(p,ic) - ci(p,ic,il)) / (1._r8 / gbc(p,ic,il) + 1.6_r8 / gs(p,ic,il))
             if (an(p,ic,il) > 0._r8 .and. abs(an(p,ic,il)-an_err) > 0.01_r8) then
                call endrun (msg=' ERROR: LeafPhotosynthesis: failed diffusion error check')
             end if

             ! Leaf water potential

             select case (gs_type)
             case (0, 1)
                lwpleaf(p,ic,il) = 0._r8
             end select

          else

             rd(p,ic,il) = 0._r8
             select case (gs_type)
             case (0, 1)
                call CiFunc (p, ic, il, mlcanopy_inst, 0._r8, ci(p,ic,il))
             case (2)
                call GsFunc (p, ic, il, mlcanopy_inst, ci(p,ic,il))
             case default
                call endrun (msg=' ERROR: LeafPhotosynthesis: gs_type not valid')
             end select
             hs(p,ic,il) = 0._r8
             vpd(p,ic,il) = 0._r8
             lwpleaf(p,ic,il) = 0._r8

          end if

          ! 13C fractionation for photosynthesis

          call C13Fractionation (p, ic, il, mlcanopy_inst)

       end do
    end do

    end associate
  end subroutine LeafPhotosynthesis

  !-----------------------------------------------------------------------
  subroutine CiFunc (p, ic, il, mlcanopy_inst, ci_val, ci_dif)
    !
    ! !DESCRIPTION:
    ! Calculate leaf photosynthesis and stomatal conductance for a specified Ci
    ! (ci_val). Then calculate a new Ci from the diffusion equation. This
    ! function returns a value ci_dif = 0 when Ci has converged to the value that
    ! satisfies the metabolic, stomatal constraint, and diffusion equations.
    !
    ! !USES:
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : qe_c4, vpd_min_MED, colim_c3a, colim_c4a, colim_c4b
    use MLclm_varctl, only : colim_type, gs_type
    use MLMathToolsMod, only : quadratic
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: p       ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)  :: ic      ! Aboveground layer index
    integer, intent(in)  :: il      ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: ci_val  ! Input value for Ci (umol/mol)
    real(r8), intent(out) :: ci_dif ! Difference in Ci
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: aquad,bquad,cquad   ! Terms for quadratic equations
    real(r8) :: r1,r2               ! Roots of quadratic equation
    real(r8) :: ai                  ! Intermediate co-limited photosynthesis (umol CO2/m2/s)
    real(r8) :: gleaf               ! Leaf CO2 conductance (mol CO2/m2/s)
    real(r8) :: cinew               ! New value for Ci
    real(r8) :: term                ! Term for stomatal conductance
    real(r8) :: vpd_term            ! Vapor pressure deficit for Medlyn stomatal conductance (kPa)
    !---------------------------------------------------------------------

    associate ( &
                                                  ! *** Input ***
    c3psn     => pftcon%c3psn                , &  ! CLM: Photosynthetic pathway (1. = C3 plant, 0. = C4 plant)
    o2ref     => mlcanopy_inst%o2ref_forcing , &  ! Atmospheric O2 at reference height (mmol/mol)
    g0        => mlcanopy_inst%g0_canopy     , &  ! Ball-Berry or Medlyn minimum leaf conductance (mol H2O/m2/s)
    g1        => mlcanopy_inst%g1_canopy     , &  ! Ball-Berry or Medlyn slope parameter
    dpai      => mlcanopy_inst%dpai_profile  , &  ! Canopy layer plant area index (m2/m2)
    cair      => mlcanopy_inst%cair_profile  , &  ! Canopy layer atmospheric CO2 (umol/mol)
    gbv       => mlcanopy_inst%gbv_leaf      , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gbc       => mlcanopy_inst%gbc_leaf      , &  ! Leaf boundary layer conductance: CO2 (mol CO2/m2 leaf/s)
    apar      => mlcanopy_inst%apar_leaf     , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    kc        => mlcanopy_inst%kc_leaf       , &  ! Leaf Michaelis-Menten constant for CO2 (umol/mol)
    ko        => mlcanopy_inst%ko_leaf       , &  ! Leaf Michaelis-Menten constant for O2 (mmol/mol)
    cp        => mlcanopy_inst%cp_leaf       , &  ! Leaf CO2 compensation point (umol/mol)
    vcmax     => mlcanopy_inst%vcmax_leaf    , &  ! Leaf maximum carboxylation rate (umol/m2/s)
    je        => mlcanopy_inst%je_leaf       , &  ! Leaf electron transport rate (umol/m2/s)
    kp        => mlcanopy_inst%kp_leaf       , &  ! Leaf C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    rd        => mlcanopy_inst%rd_leaf       , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    ceair     => mlcanopy_inst%ceair_leaf    , &  ! Leaf vapor pressure of air, constrained for stomatal conductance (Pa)
    leaf_esat => mlcanopy_inst%leaf_esat_leaf, &  ! Leaf saturation vapor pressure (Pa)
                                                  ! *** Output ***
    ac        => mlcanopy_inst%ac_leaf       , &  ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj        => mlcanopy_inst%aj_leaf       , &  ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ap        => mlcanopy_inst%ap_leaf       , &  ! Leaf product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    ag        => mlcanopy_inst%ag_leaf       , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an        => mlcanopy_inst%an_leaf       , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs        => mlcanopy_inst%cs_leaf       , &  ! Leaf surface CO2 (umol/mol)
    gs        => mlcanopy_inst%gs_leaf         &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    )

    if (dpai(p,ic) > 0._r8) then

       ! First calculate the metabolic (demand-based) photosynthetic rate

       if (nint(c3psn(patch%itype(p))) == 1) then

          ! C3: Rubisco-limited photosynthesis
          ac(p,ic,il) = vcmax(p,ic,il) * max(ci_val-cp(p,ic,il),0._r8) / (ci_val + kc(p,ic,il)*(1._r8 + o2ref(p)/ko(p,ic,il)))
 
          ! C3: RuBP regeneration-limited photosynthesis
          aj(p,ic,il) = je(p,ic,il) * max(ci_val-cp(p,ic,il),0._r8) / (4._r8 * ci_val + 8._r8 * cp(p,ic,il))

          ! C3: Product-limited photosynthesis
          ap(p,ic,il) = 0._r8

       else if (nint(c3psn(patch%itype(p))) == 0) then

          ! C4: Rubisco-limited photosynthesis
          ac(p,ic,il) = vcmax(p,ic,il)
 
          ! C4: RuBP regeneration-limited photosynthesis
          aj(p,ic,il) = qe_c4 * apar(p,ic,il)

          ! C4: PEP carboxylase-limited (CO2-limited)
          ap(p,ic,il) = kp(p,ic,il) * max(ci_val, 0._r8)

       end if

       ! Photosynthesis as the minimum or co-limited rate

       select case (colim_type)
       case (0)

          ! No co-limitation - use minimum rate

          if (nint(c3psn(patch%itype(p))) == 1) then
             ag(p,ic,il) = min(ac(p,ic,il), aj(p,ic,il))
          else if (nint(c3psn(patch%itype(p))) == 0) then
             ag(p,ic,il) = min(ac(p,ic,il), aj(p,ic,il), ap(p,ic,il))
          end if

       case (1)

          ! First co-limit Ac and Aj

          if (nint(c3psn(patch%itype(p))) == 1) then
             aquad = colim_c3a
          else if (nint(c3psn(patch%itype(p))) == 0) then
             aquad = colim_c4a
          end if
          bquad = -(ac(p,ic,il) + aj(p,ic,il))
          cquad = ac(p,ic,il) * aj(p,ic,il)
          call quadratic (aquad, bquad, cquad, r1, r2)
          ai = min(r1,r2)

          ! Now co-limit again using Ap, but only for C4 plants

          if (nint(c3psn(patch%itype(p))) == 1) then
             ag(p,ic,il) = ai
          else if (nint(c3psn(patch%itype(p))) == 0) then
             aquad = colim_c4b
             bquad = -(ai + ap(p,ic,il))
             cquad = ai * ap(p,ic,il)
             call quadratic (aquad, bquad, cquad, r1, r2)
             ag(p,ic,il) = min(r1,r2)
          end if

       case default
          call endrun (msg=' ERROR: CiFunc: colim_type not valid')
       end select

       ! Prevent photosynthesis from being negative

       ac(p,ic,il) = max(ac(p,ic,il), 0._r8)
       aj(p,ic,il) = max(aj(p,ic,il), 0._r8)
       ap(p,ic,il) = max(ap(p,ic,il), 0._r8)
       ag(p,ic,il) = max(ag(p,ic,il), 0._r8)

       ! Net photosynthesis

       an(p,ic,il) = ag(p,ic,il) - rd(p,ic,il)

       ! CO2 at leaf surface

       cs(p,ic,il) = cair(p,ic) - an(p,ic,il) / gbc(p,ic,il)
       cs(p,ic,il) = max(cs(p,ic,il), 1._r8)

       ! Now use the stomatal constraint function to calculate gs given An

       select case (gs_type)
       case (1)

          ! Ball-Berry stomatal conductance: Solve the quadratic equation
          ! for gs given An: aquad*gs^2 + bquad*gs + cquad = 0. This is
          ! obtained by substituting hs = es / leaf_esat. The correct
          ! solution is the larger of the two roots. This solution is
          ! valid for An >= 0. With An <= 0, gs = g0.

          if (an(p,ic,il) > 0._r8) then
             term = an(p,ic,il) / cs(p,ic,il)
             aquad = 1._r8
             bquad = gbv(p,ic,il) - g0(p) - g1(p) * term
             cquad = -gbv(p,ic,il) * (g0(p) + g1(p) * term * ceair(p,ic,il) / leaf_esat(p,ic,il))
             call quadratic (aquad, bquad, cquad, r1, r2)
             gs(p,ic,il) = max(r1,r2)
          else
             gs(p,ic,il) = g0(p)
          end if

       case (0)

          ! Medlyn stomatal conductance is a similar quadratic equation.
          ! The vpd term is constrained to >= vpd_min_MED to prevent
          ! infinite gs at low vpd.

          if (an(p,ic,il) > 0._r8) then
             vpd_term = max((leaf_esat(p,ic,il) - ceair(p,ic,il)), vpd_min_MED) * 0.001_r8
             term = 1.6_r8 * an(p,ic,il) / cs(p,ic,il)
             aquad = 1._r8
             bquad = -(2._r8 * (g0(p) + term) + (g1(p) * term)**2 / (gbv(p,ic,il) * vpd_term))
             cquad = g0(p) * g0(p) + (2._r8 * g0(p) + term * (1._r8 - g1(p) * g1(p) / vpd_term)) * term
             call quadratic (aquad, bquad, cquad, r1, r2)
             gs(p,ic,il) = max(r1,r2)
          else
             gs(p,ic,il) = g0(p)
          end if

       end select

       ! Now use the diffusion (supply-based) photosynthetic rate to calculate Ci

       gleaf = 1._r8 / (1._r8/gbc(p,ic,il) + 1.6_r8/gs(p,ic,il))
       cinew = cair(p,ic) - an(p,ic,il) / gleaf

       ! ci_dif is the difference between the current Ci and the new Ci

       ci_dif = cinew - ci_val
       if (an(p,ic,il) < 0._r8) ci_dif = 0._r8

    else

       ac(p,ic,il) = 0._r8
       aj(p,ic,il) = 0._r8
       ap(p,ic,il) = 0._r8
       ag(p,ic,il) = 0._r8
       an(p,ic,il) = 0._r8
       cs(p,ic,il) = 0._r8
       gs(p,ic,il) = 0._r8
       ci_dif = 0._r8

    end if

    end associate
  end subroutine CiFunc

  !-----------------------------------------------------------------------
  subroutine GsFunc (p, ic, il, mlcanopy_inst, ci_val)
    !
    ! !DESCRIPTION:
    ! Calculate leaf photosynthesis for a specified stomatal conductance.
    ! Then calculate Ci from the diffusion equation. 
    !
    ! !USES:
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : qe_c4, colim_c3a, colim_c4a, colim_c4b
    use MLclm_varctl, only : colim_type
    use MLMathToolsMod, only : quadratic
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p        ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic       ! Aboveground layer index
    integer, intent(in) :: il       ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(out) :: ci_val ! Calculated value for Ci (umol/mol)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: gleaf               ! Leaf CO2 conductance (mol CO2/m2/s)
    real(r8) :: a0,b0               ! Terms for quadratic photosynthesis calculation
    real(r8) :: aquad,bquad,cquad   ! Terms for quadratic equations
    real(r8) :: r1,r2               ! Roots of quadratic equation
    real(r8) :: ai                  ! Intermediate co-limited photosynthesis (umol CO2/m2/s)
    !---------------------------------------------------------------------

    associate ( &
                                               ! *** Input ***
    c3psn  => pftcon%c3psn                , &  ! CLM: Photosynthetic pathway (1. = C3 plant, 0. = C4 plant)
    o2ref  => mlcanopy_inst%o2ref_forcing , &  ! Atmospheric O2 at reference height (mmol/mol)
    dpai   => mlcanopy_inst%dpai_profile  , &  ! Canopy layer plant area index (m2/m2)
    cair   => mlcanopy_inst%cair_profile  , &  ! Canopy layer atmospheric CO2 (umol/mol)
    gbc    => mlcanopy_inst%gbc_leaf      , &  ! Leaf boundary layer conductance: CO2 (mol CO2/m2 leaf/s)
    gs     => mlcanopy_inst%gs_leaf       , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    apar   => mlcanopy_inst%apar_leaf     , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    kc     => mlcanopy_inst%kc_leaf       , &  ! Leaf Michaelis-Menten constant for CO2 (umol/mol)
    ko     => mlcanopy_inst%ko_leaf       , &  ! Leaf Michaelis-Menten constant for O2 (mmol/mol)
    cp     => mlcanopy_inst%cp_leaf       , &  ! Leaf CO2 compensation point (umol/mol)
    vcmax  => mlcanopy_inst%vcmax_leaf    , &  ! Leaf maximum carboxylation rate (umol/m2/s)
    je     => mlcanopy_inst%je_leaf       , &  ! Leaf electron transport rate (umol/m2/s)
    kp     => mlcanopy_inst%kp_leaf       , &  ! Leaf C4 initial slope of CO2 response curve at 25C (mol/m2/s)
    rd     => mlcanopy_inst%rd_leaf       , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
                                               ! *** Output ***
    ac     => mlcanopy_inst%ac_leaf       , &  ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj     => mlcanopy_inst%aj_leaf       , &  ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ap     => mlcanopy_inst%ap_leaf       , &  ! Leaf product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    ag     => mlcanopy_inst%ag_leaf       , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an     => mlcanopy_inst%an_leaf       , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs     => mlcanopy_inst%cs_leaf         &  ! Leaf surface CO2 (umol/mol)
    )

    ! Calculate leaf photosynthesis for a specified stomatal conductance
    ! by substituting the diffusion equation for Ci into the metabolic
    ! (demand-based) photosynthetic equation. The resulting quadratic
    ! equation is: a An**2 + b An + c = 0

    if (dpai(p,ic) > 0._r8) then

       ! Leaf conductance: gbc has units mol CO2/m2/s, gs has units mol H2O/m2/s,
       ! gleaf has units mol CO2/m2/s

       gleaf = 1._r8 / (1._r8/gbc(p,ic,il) + 1.6_r8/gs(p,ic,il))

       ! Gross assimilation rates

       if (nint(c3psn(patch%itype(p))) == 1) then

          ! C3: Rubisco-limited photosynthesis

          a0 = vcmax(p,ic,il)
          b0 = kc(p,ic,il) * (1._r8 + o2ref(p) / ko(p,ic,il))

          aquad = 1._r8 / gleaf
          bquad = -(cair(p,ic) + b0) - (a0 - rd(p,ic,il)) / gleaf
          cquad = a0 * (cair(p,ic) - cp(p,ic,il)) - rd(p,ic,il) * (cair(p,ic) + b0)

          call quadratic (aquad, bquad, cquad, r1, r2)
          ac(p,ic,il) = min(r1,r2) + rd(p,ic,il)
 
          ! C3: RuBP regeneration-limited photosynthesis

          a0 = je(p,ic,il) / 4._r8
          b0 = 2._r8 * cp(p,ic,il)

          aquad = 1._r8 / gleaf
          bquad = -(cair(p,ic) + b0) - (a0 - rd(p,ic,il)) / gleaf
          cquad = a0 * (cair(p,ic) - cp(p,ic,il)) - rd(p,ic,il) * (cair(p,ic) + b0)

          call quadratic (aquad, bquad, cquad, r1, r2)
          aj(p,ic,il) = min(r1,r2) + rd(p,ic,il)

          ! C3: Product-limited photosynthesis

          ap(p,ic,il) = 0._r8

       else if (nint(c3psn(patch%itype(p))) == 0) then

          ! C4: Rubisco-limited photosynthesis
          ac(p,ic,il) = vcmax(p,ic,il)
 
          ! C4: RuBP-limited photosynthesis
          aj(p,ic,il) = qe_c4 * apar(p,ic,il)

          ! C4: PEP carboxylase-limited (CO2-limited)
          ap(p,ic,il) = kp(p,ic,il) * (cair(p,ic) * gleaf + rd(p,ic,il)) / (gleaf + kp(p,ic,il))

       end if

       ! Photosynthesis as the minimum or co-limited rate

       select case (colim_type)
       case (0)

          ! No co-limitation - use minimum rate

          if (nint(c3psn(patch%itype(p))) == 1) then
             ag(p,ic,il) = min(ac(p,ic,il), aj(p,ic,il))
          else if (nint(c3psn(patch%itype(p))) == 0) then
             ag(p,ic,il) = min(ac(p,ic,il), aj(p,ic,il), ap(p,ic,il))
          end if

       case (1)

          ! First co-limit Ac and Aj

          if (nint(c3psn(patch%itype(p))) == 1) then
             aquad = colim_c3a
          else if (nint(c3psn(patch%itype(p))) == 0) then
             aquad = colim_c4a
          end if
          bquad = -(ac(p,ic,il) + aj(p,ic,il))
          cquad = ac(p,ic,il) * aj(p,ic,il)
          call quadratic (aquad, bquad, cquad, r1, r2)
          ai = min(r1,r2)

          ! Now co-limit again using Ap, but only for C4 plants

          if (nint(c3psn(patch%itype(p))) == 1) then
             ag(p,ic,il) = ai
          else if (nint(c3psn(patch%itype(p))) == 0) then
             aquad = colim_c4b
             bquad = -(ai + ap(p,ic,il))
             cquad = ai * ap(p,ic,il)
             call quadratic (aquad, bquad, cquad, r1, r2)
             ag(p,ic,il) = min(r1,r2)
          end if

       case default
          call endrun (msg=' ERROR: GsFunc: colim_type not valid')
       end select

       ! Net photosynthesis

       an(p,ic,il) = ag(p,ic,il) - rd(p,ic,il)

       ! CO2 at leaf surface

       cs(p,ic,il) = cair(p,ic) - an(p,ic,il) / gbc(p,ic,il)
       cs(p,ic,il) = max(cs(p,ic,il), 1._r8)

       ! Intercelluar CO2

       ci_val = cair(p,ic) - an(p,ic,il) / gleaf

    else

       ac(p,ic,il) = 0._r8
       aj(p,ic,il) = 0._r8
       ap(p,ic,il) = 0._r8
       ag(p,ic,il) = 0._r8
       an(p,ic,il) = 0._r8
       cs(p,ic,il) = 0._r8
       ci_val = 0._r8

    end if

    end associate
  end subroutine GsFunc

  !-----------------------------------------------------------------------
  subroutine StomataOptimization (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Photosynthesis and stomatal conductance with optimization
    !
    ! !USES:
    use MLclm_varcon, only : gsmin_SPA
    use MLMathToolsMod, only : zbrent
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p              ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic             ! Aboveground layer index
    integer, intent(in) :: il             ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: gs1, gs2                  ! Initial guess for gs (mol H2O/m2/s)
    real(r8) :: check1, check2            ! Water-use efficiency and cavitation check for gs1 and gs2
    real(r8), parameter :: tol = 0.004_r8 ! gs is updated to accuracy tol (mol H2O/m2/s)
    !---------------------------------------------------------------------

    associate ( &
    dpai      => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    lwpveg    => mlcanopy_inst%lwpveg_profile  , &  ! Canopy layer leaf water potential (MPa)
    lwpleaf   => mlcanopy_inst%lwpleaf_leaf    , &  ! Leaf water potential (MPa)
    gs        => mlcanopy_inst%gs_leaf           &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    )

    ! Initialize leaf water potential (for sunlit or shaded leaf) to the layer
    ! value of the previous time step

    lwpleaf(p,ic,il) = lwpveg(p,ic)

    ! Low and high initial estimates for gs (mol H2O/m2/s)

    gs1 = gsmin_SPA
    gs2 = 2._r8

    ! Calculate gs

    if (dpai(p,ic) > 0._r8) then

       ! Check for minimum stomatal conductance linked to low light or drought stress
       ! based on the water-use efficiency and cavitation checks for gs1 and gs2

       call StomataEfficiency (p, ic, il, mlcanopy_inst, gs1, check1)
       call StomataEfficiency (p, ic, il, mlcanopy_inst, gs2, check2)

       if (check1 * check2 < 0._r8) then

          ! Calculate gs using the function StomataEfficiency to iterate gs
          ! to an accuracy of tol (mol H2O/m2/s)

          gs(p,ic,il) = zbrent ('StomataOptimization', p, ic, il, mlcanopy_inst, StomataEfficiency, gs1, gs2, tol)

       else

          ! Low light or drought stress. Set gs to minimum conductance

          gs(p,ic,il) = gsmin_SPA

       end if

    else

       gs(p,ic,il) = 0._r8

    end if

    ! Leaf fluxes and leaf water potential for this gs

    call StomataFluxes (p, ic, il, mlcanopy_inst, gs(p,ic,il), lwpleaf(p,ic,il))

    end associate
  end subroutine StomataOptimization

  !-----------------------------------------------------------------------
  subroutine StomataEfficiency (p, ic, il, mlcanopy_inst, gs_val, val)
    !
    ! !DESCRIPTION:
    ! Water-use efficiency check and cavitation check to determine optimal gs. 
    ! For the stomatal conductance gs_val, calculate photosynthesis and leaf
    ! water potential for an increase in stomatal conductance equal to "delta".
    ! The returned value is positive if this increase produces a change in
    ! photosynthesis > iota*vpd*delta or if the leaf water potential is > minlwp.
    ! The returned value is negative if the increase produces a change in
    ! photosynthesis < iota*vpd*delta or if the leaf water potential is < minlwp. 
    !
    ! !USES:
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p         ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic        ! Aboveground layer index
    integer, intent(in) :: il        ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: gs_val   ! Value for gs to use in calculations
    real(r8), intent(out) :: val     ! Returned minimum of the two stomatal checks
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: delta                ! Small difference for gs (mol H2O/m2/s)
    real(r8) :: leafwp               ! Current leaf water potential (MPa)
    real(r8) :: gs2                  ! Lower value for gs (mol H2O/m2/s)
    real(r8) :: an2                  ! Leaf photosynthesis at gs2 (umol CO2/m2/s)
    real(r8) :: gs1                  ! Higher value for gs (mol H2O/m2/s)
    real(r8) :: an1                  ! Leaf photosynthesis at gs1 (umol CO2/m2/s)
    real(r8) :: wue                  ! Water-use efficiency check
    real(r8) :: minpsi               ! Cavitation check
    !---------------------------------------------------------------------

    associate ( &
    minlwp_SPA  => pftcon%minlwp_SPA          , &  ! CLM (new): Minimum leaf water potential (MPa)
    iota_SPA    => pftcon%iota_SPA            , &  ! CLM (new): Stomatal water-use efficiency (umol CO2/ mol H2O)
    pref        => mlcanopy_inst%pref_forcing , &  ! Air pressure at reference height (Pa)
    lwpleaf     => mlcanopy_inst%lwpleaf_leaf , &  ! Leaf water potential (MPa)
    an          => mlcanopy_inst%an_leaf      , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    vpd         => mlcanopy_inst%vpd_leaf       &  ! Leaf vapor pressure deficit (Pa)
    )

    ! Specify "delta" as a small difference in gs (mol H2O/m2/s)

    delta = 0.001_r8

    ! Photosynthesis at lower gs (gs_val - delta)

    leafwp = lwpleaf(p,ic,il)
    gs2 = gs_val - delta
    call StomataFluxes (p, ic, il, mlcanopy_inst, gs2, leafwp)
    an2 = an(p,ic,il)

    ! Photosynthesis at higher gs (gs_val)

    leafwp = lwpleaf(p,ic,il)
    gs1 = gs_val
    call StomataFluxes (p, ic, il, mlcanopy_inst, gs1, leafwp)
    an1 = an(p,ic,il)

    ! Efficiency check: wue < 0 when d(An) / d(gs) < iota * vpd

    wue = (an1 - an2) - iota_SPA(patch%itype(p)) * delta * (vpd(p,ic,il) / pref(p))

    ! Cavitation check: minpsi < 0 when leafwp < minlwp_SPA

    minpsi = leafwp - minlwp_SPA(patch%itype(p))

    ! Return the minimum of the two checks

    val = min(wue, minpsi)

    end associate
  end subroutine StomataEfficiency

  !-----------------------------------------------------------------------
  subroutine StomataFluxes (p, ic, il, mlcanopy_inst, gs_val, leafwp)
    !
    ! !DESCRIPTION:
    ! Leaf calculations for a specified stomatal conductance (gs_val)
    !
    ! !USES:
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p             ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic            ! Aboveground layer index
    integer, intent(in) :: il            ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: gs_val       ! Value for gs to use in calculations
    real(r8), intent(inout) :: leafwp    ! Leaf water potential (MPa)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    associate ( &
    eair        => mlcanopy_inst%eair_profile    , &  ! Canopy layer vapor pressure (Pa)
    gs          => mlcanopy_inst%gs_leaf         , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    ci          => mlcanopy_inst%ci_leaf         , &  ! Leaf intercellular CO2
    hs          => mlcanopy_inst%hs_leaf         , &  ! Leaf fractional humidity at leaf surface (dimensionless)
    vpd         => mlcanopy_inst%vpd_leaf        , &  ! Leaf vapor pressure deficit (Pa)
    gbv         => mlcanopy_inst%gbv_leaf        , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    leaf_esat   => mlcanopy_inst%leaf_esat_leaf    &  ! Leaf saturation vapor pressure (Pa)
    )

    ! Use specified gs (gs_val)

    gs(p,ic,il) = gs_val

    ! Leaf photosynthesis for gs

    call GsFunc (p, ic, il, mlcanopy_inst, ci(p,ic,il))

    ! Relative humidity and vapor pressure at leaf surface

    hs(p,ic,il) = (gbv(p,ic,il)*eair(p,ic) + gs(p,ic,il)*leaf_esat(p,ic,il)) &
                / ((gbv(p,ic,il) + gs(p,ic,il)) * leaf_esat(p,ic,il))
    vpd(p,ic,il) = max(leaf_esat(p,ic,il) - hs(p,ic,il)*leaf_esat(p,ic,il), 0.1_r8)

    ! Leaf transpiration

    call LeafTranspiration (p, ic, il, mlcanopy_inst)

    ! Leaf water potential

    call LeafWaterPotential (p, ic, il, mlcanopy_inst, leafwp)

    end associate
  end subroutine StomataFluxes

  !-----------------------------------------------------------------------
  subroutine LeafWaterPotential (p, ic, il, mlcanopy_inst, leafwp)
    !
    ! !DESCRIPTION:
    ! Calculate leaf water potential for the current transpiration rate
    !
    ! !USES:
    use clm_varcon, only : denh2o, grav
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varctl, only : dtime_substep
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p                  ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic                 ! Aboveground layer index
    integer, intent(in) :: il                 ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(inout) :: leafwp         ! Leaf water potential (MPa)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: head                          ! Head of pressure  (MPa/m)
    real(r8) :: dtime                         ! Model time step (s)
    real(r8) :: y0                            ! Leaf water potential at beginning of timestep (MPa)
    real(r8) :: dy                            ! Change in leaf water potential (MPa)
    real(r8) :: a, b                          ! Intermediate calculation
    !---------------------------------------------------------------------

    associate ( &
    capac_SPA   => pftcon%capac_SPA           , &  ! CLM (new): Plant capacitance (mmol H2O/m2 leaf area/MPa)
    psis        => mlcanopy_inst%psis_soil    , &  ! Weighted soil water potential (MPa)
    dpai        => mlcanopy_inst%dpai_profile , &  ! Canopy layer plant area index (m2/m2)
    zs          => mlcanopy_inst%zs_profile   , &  ! Canopy layer height for scalar concentration and source (m)
    lsc         => mlcanopy_inst%lsc_profile  , &  ! Canopy layer leaf-specific conductance (mmol H2O/m2 leaf/s/MPa)
    trleaf      => mlcanopy_inst%trleaf_leaf    &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    )

    head = denh2o * grav * 1.e-06_r8

    ! Get step size

    dtime = dtime_substep

    ! Change in leaf water potential is: dy / dt = (a - y) / b. The integrated change 
    ! over a full model timestep is: dy = (a - y0) * (1 - exp(-dt/b))

    if (dpai(p,ic) > 0._r8) then
       y0 = leafwp
       a = psis(p) - head *  zs(p,ic) - 1000._r8 * trleaf(p,ic,il) / lsc(p,ic)
       b = capac_SPA(patch%itype(p)) / lsc(p,ic)
       dy = (a - y0) * (1._r8 - exp(-dtime/b))
       leafwp = y0 + dy
    else
       leafwp = 0._r8
    end if

    end associate
  end subroutine LeafWaterPotential

  !-----------------------------------------------------------------------
  subroutine LeafTranspiration (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Calculate leaf water transpiration flux
    !
    ! !USES:
    use MLWaterVaporMod, only : SatVap, LatVap
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p            ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic           ! Aboveground layer index
    integer, intent(in) :: il           ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: esat                    ! Saturation vapor pressure (Pa)
    real(r8) :: desat                   ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(r8) :: lambda                  ! Latent heat of vaporization (J/mol)
    real(r8) :: gleaf                   ! Leaf conductance for transpiration (mol H2O/m2 leaf/s)
    !---------------------------------------------------------------------

    associate ( &
                                                  ! *** Input ***
    tref      => mlcanopy_inst%tref_forcing  , &  ! Air temperature at reference height (K)
    pref      => mlcanopy_inst%pref_forcing  , &  ! Air pressure at reference height (Pa)
    dpai      => mlcanopy_inst%dpai_profile  , &  ! Canopy layer plant area index (m2/m2)
    eair      => mlcanopy_inst%eair_profile  , &  ! Canopy layer vapor pressure (Pa)
    fdry      => mlcanopy_inst%fdry_profile  , &  ! Canopy layer fraction of plant area index that is green and dry
    tleaf     => mlcanopy_inst%tleaf_leaf    , &  ! Leaf temperature (K)
    gs        => mlcanopy_inst%gs_leaf       , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv_leaf      , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
                                                  ! *** Output ***
    trleaf    => mlcanopy_inst%trleaf_leaf     &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    )

    if (dpai(p,ic) > 0._r8) then

       ! Saturation vapor pressure

       call SatVap (tleaf(p,ic,il), esat, desat)

       ! Latent heat of vaporization

       lambda = LatVap(tref(p))

       ! Leaf conductance for transpiration

       gleaf = gs(p,ic,il) * gbv(p,ic,il) / (gs(p,ic,il) + gbv(p,ic,il))

       ! Transpiration flux: mol H2O/m2/s

       trleaf(p,ic,il) = gleaf * fdry(p,ic) * (esat - eair(p,ic)) / pref(p)

    else

       trleaf(p,ic,il) = 0._r8

    end if

    end associate
  end subroutine LeafTranspiration

  !-----------------------------------------------------------------------
  subroutine C13Fractionation (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! 13C fractionation for photosynthesis
    !
    ! !USES:
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p        ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic       ! Aboveground layer index
    integer, intent(in) :: il       ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    associate ( &
                                                  ! *** Input ***
    c3psn     => pftcon%c3psn                , &  ! CLM: Photosynthetic pathway (1. = C3 plant, 0. = C4 plant)
    dpai      => mlcanopy_inst%dpai_profile  , &  ! Canopy layer plant area index (m2/m2)
    cair      => mlcanopy_inst%cair_profile  , &  ! Canopy layer atmospheric CO2 (umol/mol)
    apar      => mlcanopy_inst%apar_leaf     , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    ci        => mlcanopy_inst%ci_leaf       , &  ! Leaf intercellular CO2 (umol/mol)
                                                  ! *** Output ***
    alphapsn => mlcanopy_inst%alphapsn_leaf    &  ! Leaf 13C fractionation factor for photosynthesis (-)
    )

    if (dpai(p,ic) > 0._r8) then

       if (apar(p,ic,il) > 0._r8) then
          if (nint(c3psn(patch%itype(p))) == 1) then
             alphapsn(p,ic,il) = 1._r8 + (4.4_r8 + 22.6_r8 * ci(p,ic,il) / cair(p,ic)) / 1000._r8
          else if (nint(c3psn(patch%itype(p))) == 0) then
             alphapsn(p,ic,il) = 1._r8 + 4.4_r8 / 1000._r8
          end if
       else
          alphapsn(p,ic,il) = 1._r8
       end if

    else

       alphapsn(p,ic,il) = 0._r8

    end if

    end associate
  end subroutine C13Fractionation

end module MLLeafPhotosynthesisMod
