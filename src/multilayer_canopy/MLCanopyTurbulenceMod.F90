module MLCanopyTurbulenceMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Scalar source/sink fluxes and scalar profiles
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
  public :: CanopyTurbulence        ! Main routine for scalar source/sink fluxes and scalar profiles
  public :: LookupPsihatINI         ! Initialize the RSL psihat look-up tables
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: WellMixed              ! Canopy scalar profiles equal reference height values
  private :: HF2008                 ! Scalar source/sink fluxes and scalar profiles using H&F (2008) RSL theory
  private :: ObuFunc                ! Subroutine to solve for the Obukhov length
  private :: GetBeta                ! Calculate beta = u* / u at canopy top
  private :: GetPrSc                ! Prandlt number (Pr) and Schmidt number (Sc) at canopy top
  private :: GetPsiRSL              ! RSL-modified stability functions
  private :: phim_monin_obukhov     ! Monin-Obukhov phi stability function for momentum
  private :: phic_monin_obukhov     ! Monin-Obukhov phi stability function for scalars
  private :: psim_monin_obukhov     ! Monin-Obukhov psi stability function for momentum
  private :: psic_monin_obukhov     ! Monin-Obukhov psi stability function for scalars
  private :: LookupPsihat           ! Determines the RSL function psihat as provided through a look-up table
  private :: WindProfile            ! Wind speed profile above and within canopy
  private :: AerodynamicConductance ! Aerodynamic conductances above and within canopy
  private :: FluxProfileSolution    ! Scalar source/sink fluxes and concentration profiles using implicit solution

contains

  !-----------------------------------------------------------------------
  subroutine CanopyTurbulence (niter, num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy turbulence, scalar source/sink fluxes for leaves and soil, and
    ! scalar profiles above and within the canopy
    !
    ! !USES:
    use MLclm_varctl, only : turb_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: niter        ! Iteration index
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    select case (turb_type)
    case (0, -1)

       ! Use the well-mixed assumption, or read profile data from dataset

       call WellMixed (num_filter, filter, mlcanopy_inst)

    case (1)

       ! Use Harman & Finnigan (2008) roughness sublayer theory

       call HF2008 (niter, num_filter, filter, mlcanopy_inst)

    case default

       call endrun (msg=' ERROR: CanopyTurbulence: turb_type not valid')

    end select

  end subroutine CanopyTurbulence

  !-----------------------------------------------------------------------
  subroutine WellMixed (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy scalar profiles equal reference height values
    ! (well-mixed assumption) or are read in from dataset
    !
    ! !USES:
    use MLclm_varctl, only : turb_type
    use MLclm_varcon, only : cd
    use MLclm_varpar, only : isun, isha
    use MLLeafFluxesMod, only : LeafFluxes
    use MLMathToolsMod, only : hybrid
    use MLSoilFluxesMod, only : SoilFluxes
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fp                       ! Filter index
    integer :: p                        ! Patch index for CLM g/l/c/p hierarchy
    integer :: ic                       ! Aboveground layer index
    integer :: il                       ! Sunlit (1) or shaded (2) leaf index
    real(r8) :: obu0, obu1              ! Initial estimates for Obukhov length (m)
    real(r8) :: tol                     ! Accuracy tolerance for Obukhov length (m)
    real(r8) :: dummy                   ! Dummy argument
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *** Input ***
    uref      => mlcanopy_inst%uref_forcing     , &  ! Wind speed at reference height (m/s)
    tref      => mlcanopy_inst%tref_forcing     , &  ! Air temperature at reference height (K)
    eref      => mlcanopy_inst%eref_forcing     , &  ! Vapor pressure at reference height (Pa)
    co2ref    => mlcanopy_inst%co2ref_forcing   , &  ! Atmospheric CO2 at reference height (umol/mol)
    qref      => mlcanopy_inst%qref_forcing     , &  ! Specific humidity at reference height (kg/kg)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    ztop      => mlcanopy_inst%ztop_canopy      , &  ! Canopy foliage top height (m)
    lai       => mlcanopy_inst%lai_canopy       , &  ! Leaf area index of canopy (m2/m2)
    sai       => mlcanopy_inst%sai_canopy       , &  ! Stem area index of canopy (m2/m2)
    wind_data => mlcanopy_inst%wind_data_profile, &  ! Canopy layer wind speed FROM DATASET (m/s)
    tair_data => mlcanopy_inst%tair_data_profile, &  ! Canopy layer air temperature FROM DATASET (K)
    eair_data => mlcanopy_inst%eair_data_profile, &  ! Canopy layer vapor pressure FROM DATASET (Pa)
                                                     ! *** Output ***
    Lc        => mlcanopy_inst%Lc_canopy        , &  ! Canopy density length scale (m)
    ustar     => mlcanopy_inst%ustar_canopy     , &  ! Friction velocity (m/s)
    uaf       => mlcanopy_inst%uaf_canopy       , &  ! Wind speed at canopy top (m/s)
    taf       => mlcanopy_inst%taf_canopy       , &  ! Air temperature at canopy top (K)
    qaf       => mlcanopy_inst%qaf_canopy       , &  ! Specific humidity at canopy top (kg/kg)
    gac0      => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    wind      => mlcanopy_inst%wind_profile     , &  ! Canopy layer wind speed (m/s)
    tair      => mlcanopy_inst%tair_profile     , &  ! Canopy layer air temperature (K)
    eair      => mlcanopy_inst%eair_profile     , &  ! Canopy layer vapor pressure (Pa)
    cair      => mlcanopy_inst%cair_profile     , &  ! Canopy layer atmospheric CO2 (umol/mol)
    gac       => mlcanopy_inst%gac_profile      , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    shair     => mlcanopy_inst%shair_profile    , &  ! Canopy layer air sensible heat flux (W/m2)
    etair     => mlcanopy_inst%etair_profile    , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair     => mlcanopy_inst%stair_profile    , &  ! Canopy layer air storage heat flux (W/m2)
    mflx      => mlcanopy_inst%mflx_profile     , &  ! Canopy layer momentum flux (m2/s2)
    ! From LeafFluxes
    tleaf     => mlcanopy_inst%tleaf_leaf       , &  ! Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf_leaf      , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf_leaf      , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf    => mlcanopy_inst%lhleaf_leaf      , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf    => mlcanopy_inst%trleaf_leaf      , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf    => mlcanopy_inst%evleaf_leaf      , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    ! From SoilFluxes
    shsoi     => mlcanopy_inst%shsoi_soil       , &  ! Sensible heat flux: ground (W/m2)
    lhsoi     => mlcanopy_inst%lhsoi_soil       , &  ! Latent heat flux: ground (W/m2)
    etsoi     => mlcanopy_inst%etsoi_soil       , &  ! Water vapor flux: ground (mol H2O/m2/s)
    gsoi      => mlcanopy_inst%gsoi_soil        , &  ! Soil heat flux (W/m2)
    eg        => mlcanopy_inst%eg_soil          , &  ! Soil surface vapor pressure (Pa)
    tg        => mlcanopy_inst%tg_soil            &  ! Soil surface temperature (K)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Canopy density length scale

       Lc(p) = ztop(p) / (cd * (lai(p) + sai(p)))

       ! These are not used, but are needed to pass into the hybrid root solver

       ic = 0 ; il = 0

       ! Calculate Obukhov length (obu) using the subroutine ObuFunc to iterate until
       ! the change in obu is < tol. Do not use the returned value (dummy), and
       ! instead use the value used to calculate ustar

       obu0 = 100._r8          ! Initial estimate for Obukhov length (m)
       obu1 = -100._r8         ! Initial estimate for Obukhov length (m)
       tol = 0.1_r8            ! Accuracy tolerance for Obukhov length (m)

       dummy = hybrid ('WellMixed', p, ic, il, mlcanopy_inst, ObuFunc, obu0, obu1, tol)

       ! Scalar profiles and vertical fluxes

       do ic = 1, ncan(p)
          cair(p,ic) = co2ref(p)

          select case (turb_type)

          ! Use the well-mixed assumption

          case (0)
             wind(p,ic) = uref(p)
             tair(p,ic) = tref(p)
             eair(p,ic) = eref(p)

          ! Use profiles from dataset

          case (-1)
             wind(p,ic) = wind_data(p,ic)
             tair(p,ic) = tair_data(p,ic)
             eair(p,ic) = eair_data(p,ic)
             ! Set each profile individually to WMA if desired
             ! wind(p,ic) = uref(p)
             ! tair(p,ic) = tref(p)
             ! eair(p,ic) = eref(p)

          end select

          shair(p,ic) = 0._r8
          etair(p,ic) = 0._r8
          stair(p,ic) = 0._r8
          mflx(p,ic) = 0._r8
       end do

       ! Calculate leaf fluxes (per unit leaf area)

       do ic = 1, ncan(p)
          call LeafFluxes (p, ic, isun, mlcanopy_inst)
          call LeafFluxes (p, ic, isha, mlcanopy_inst)
       end do

       ! Calculate soil fluxes, but need gac0. Use a large resistance
       ! to approximate bare ground and so that soil fluxes are small.

       gac0(p) = (1._r8 / 100._r8) * 42.3_r8

       call SoilFluxes (p, mlcanopy_inst)

       ! Only needed for output files. These cannot be zero
       ! for analysis package to work.

       uaf(p) = uref(p)
       taf(p) = tref(p)
       qaf(p) = qref(p)
       do ic = 1, ncan(p)
          gac(p,ic) = (1._r8 / 10._r8) * 42.3_r8  ! small non-zero resistance
       end do

    end do

    end associate
  end subroutine WellMixed

  !-----------------------------------------------------------------------
  subroutine HF2008 (niter, num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy turbulence, aerodynamic conductances, wind/temperature/water vapor
    ! profiles, and scalar source/sink fluxes (leaves, soil)  using above- and 
    ! within-canopy coupling with the Harman and Finnigan (2008) roughness
    ! sublayer (RSL) parameterization 
    !
    ! !USES:
    use MLclm_varcon, only : mmh2o, mmdry, cd, eta_max
    use MLMathToolsMod, only : hybrid
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: niter        ! Iteration index
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                      ! Aboveground layer index
    integer  :: il                      ! Sunlit (1) or shaded (2) leaf index
    real(r8) :: obu0, obu1              ! Initial estimates for Obukhov length (m)
    real(r8) :: tol                     ! Accuracy tolerance for Obukhov length (m)
    real(r8) :: dummy                   ! Dummy argument
    real(r8) :: lm                      ! Length scale for momentum (m)
    real(r8) :: lm_over_beta            ! lm / beta
    real(r8) :: eta                     ! Used to limit value for lm / beta
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *************
                                                     ! *** Input ***
                                                     ! *************
    zref      => mlcanopy_inst%zref_forcing     , &  ! Atmospheric reference height (m)
    uref      => mlcanopy_inst%uref_forcing     , &  ! Wind speed at reference height (m/s)
    tref      => mlcanopy_inst%tref_forcing     , &  ! Air temperature at reference height (K)
    thref     => mlcanopy_inst%thref_forcing    , &  ! Atmospheric potential temperature at reference height (K)
    thvref    => mlcanopy_inst%thvref_forcing   , &  ! Atmospheric virtual potential temperature at reference height (K)
    eref      => mlcanopy_inst%eref_forcing     , &  ! Vapor pressure at reference height (Pa)
    qref      => mlcanopy_inst%qref_forcing     , &  ! Specific humidity at reference height (kg/kg)
    pref      => mlcanopy_inst%pref_forcing     , &  ! Air pressure at reference height (Pa)
    co2ref    => mlcanopy_inst%co2ref_forcing   , &  ! Atmospheric CO2 at reference height (umol/mol)
    rhomol    => mlcanopy_inst%rhomol_forcing   , &  ! Molar density at reference height (mol/m3)
    cpair     => mlcanopy_inst%cpair_forcing    , &  ! Specific heat of air (constant pressure) at reference height (J/mol/K)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    ntop      => mlcanopy_inst%ntop_canopy      , &  ! Index for top leaf layer
    ztop      => mlcanopy_inst%ztop_canopy      , &  ! Canopy foliage top height (m)
    lai       => mlcanopy_inst%lai_canopy       , &  ! Leaf area index of canopy (m2/m2)
    sai       => mlcanopy_inst%sai_canopy       , &  ! Stem area index of canopy (m2/m2)
    rhg       => mlcanopy_inst%rhg_soil         , &  ! Relative humidity of airspace at soil surface (fraction)
    rnsoi     => mlcanopy_inst%rnsoi_soil       , &  ! Net radiation: ground (W/m2)
    tg_bef    => mlcanopy_inst%tg_bef_soil      , &  ! Soil surface temperature for previous timestep (K)
    soilres   => mlcanopy_inst%soilres_soil     , &  ! Soil evaporative resistance (s/m)
    soil_t    => mlcanopy_inst%soil_t_soil      , &  ! Temperature of first snow/soil layer (K)
    soil_dz   => mlcanopy_inst%soil_dz_soil     , &  ! Depth to temperature of first snow/soil layer (m)
    soil_tk   => mlcanopy_inst%soil_tk_soil     , &  ! Thermal conductivity of first snow/soil layer (W/m/K)
    zw        => mlcanopy_inst%zw_profile       , &  ! Canopy height at interface between two adjacent layers (m)
    zs        => mlcanopy_inst%zs_profile       , &  ! Canopy layer height for scalar concentration and source (m)
    dz        => mlcanopy_inst%dz_profile       , &  ! Canopy layer thickness (m)
    dpai      => mlcanopy_inst%dpai_profile     , &  ! Canopy layer plant area index (m2/m2)
    fwet      => mlcanopy_inst%fwet_profile     , &  ! Canopy layer fraction of plant area index that is wet
    fdry      => mlcanopy_inst%fdry_profile     , &  ! Canopy layer fraction of plant area index that is green and dry
    fracsun   => mlcanopy_inst%fracsun_profile  , &  ! Canopy layer sunlit fraction (-)
    cpleaf    => mlcanopy_inst%cpleaf_profile   , &  ! Canopy layer leaf heat capacity (J/m2 leaf/K)
    tair_bef  => mlcanopy_inst%tair_bef_profile , &  ! Canopy layer air temperature for previous timestep (K)
    eair_bef  => mlcanopy_inst%eair_bef_profile , &  ! Canopy layer vapor pressure for previous timestep (Pa)
    gbh       => mlcanopy_inst%gbh_leaf         , &  ! Leaf boundary layer conductance: heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv_leaf         , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gs        => mlcanopy_inst%gs_leaf          , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    rnleaf    => mlcanopy_inst%rnleaf_leaf      , &  ! Leaf net radiation (W/m2 leaf)
    tleaf_bef => mlcanopy_inst%tleaf_bef_leaf   , &  ! Leaf temperature for previous timestep (K)

                                                     ! ********************
                                                     ! *** Input/Output ***
                                                     ! ********************
    taf       => mlcanopy_inst%taf_canopy       , &  ! Air temperature at canopy top (K)
    qaf       => mlcanopy_inst%qaf_canopy       , &  ! Specific humidity at canopy top (kg/kg)
    obuold    => mlcanopy_inst%obuold_canopy    , &  ! Obukhov length from previous iteration
    nmozsgn   => mlcanopy_inst%nmozsgn_canopy   , &  ! Number of times stability changes sign during iteration

                                                     ! **************
                                                     ! *** Output ***
                                                     ! **************
    ! From HF2008
    Lc        => mlcanopy_inst%Lc_canopy        , &  ! Canopy density length scale (m)
    cair      => mlcanopy_inst%cair_profile     , &  ! Canopy layer atmospheric CO2 (umol/mol)
    mflx      => mlcanopy_inst%mflx_profile     , &  ! Canopy layer momentum flux (m2/s2)
    ! From ObuFunc
    zdisp     => mlcanopy_inst%zdisp_canopy     , &  ! Displacement height (m)
    beta      => mlcanopy_inst%beta_canopy      , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy      , &  ! Prandtl (Schmidt) number at canopy top (-)
    ustar     => mlcanopy_inst%ustar_canopy     , &  ! Friction velocity (m/s)
    gac_to_hc => mlcanopy_inst%gac_to_hc_canopy , &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    obu       => mlcanopy_inst%obu_canopy       , &  ! Obukhov length (m)
    ! From WindProfile
    uaf       => mlcanopy_inst%uaf_canopy       , &  ! Wind speed at canopy top (m/s)
    wind      => mlcanopy_inst%wind_profile     , &  ! Canopy layer wind speed (m/s)
    ! From AerodynamicConductance
    gac0      => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    gac       => mlcanopy_inst%gac_profile      , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    ! From ScalarProfile
    tair      => mlcanopy_inst%tair_profile     , &  ! Canopy layer air temperature (K)
    eair      => mlcanopy_inst%eair_profile     , &  ! Canopy layer vapor pressure (Pa)
    shair     => mlcanopy_inst%shair_profile    , &  ! Canopy layer air sensible heat flux (W/m2)
    etair     => mlcanopy_inst%etair_profile    , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair     => mlcanopy_inst%stair_profile    , &  ! Canopy layer air storage heat flux (W/m2)
    ! From LeafFluxes (ScalarProfile)
    tleaf     => mlcanopy_inst%tleaf_leaf       , &  ! Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf_leaf      , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf_leaf      , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf    => mlcanopy_inst%lhleaf_leaf      , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf    => mlcanopy_inst%trleaf_leaf      , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf    => mlcanopy_inst%evleaf_leaf      , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    ! From SoilFluxes (ScalarProfile)
    shsoi     => mlcanopy_inst%shsoi_soil       , &  ! Sensible heat flux: ground (W/m2)
    lhsoi     => mlcanopy_inst%lhsoi_soil       , &  ! Latent heat flux: ground (W/m2)
    etsoi     => mlcanopy_inst%etsoi_soil       , &  ! Water vapor flux: ground (mol H2O/m2/s)
    gsoi      => mlcanopy_inst%gsoi_soil        , &  ! Soil heat flux (W/m2)
    eg        => mlcanopy_inst%eg_soil          , &  ! Soil surface vapor pressure (Pa)
    tg        => mlcanopy_inst%tg_soil            &  ! Soil surface temperature (K)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Initialization for first iteration

!      if (niter == 1) then
!         obuold(p) = 0._r8
!         nmozsgn(p) = 0
!      end if

       ! Canopy density length scale

       Lc(p) = ztop(p) / (cd * (lai(p) + sai(p)))

       ! These are not used, but are needed to pass into the hybrid root solver

       ic = 0 ; il = 0

       ! Calculate Obukhov length (obu) using the subroutine ObuFunc to iterate until
       ! the change in obu is < tol. Do not use the returned value (dummy), and
       ! instead use the value used to calculate ustar

       obu0 = 100._r8          ! Initial estimate for Obukhov length (m)
       obu1 = -100._r8         ! Initial estimate for Obukhov length (m)
       tol = 0.1_r8            ! Accuracy tolerance for Obukhov length (m)

       dummy = hybrid ('HF2008', p, ic, il, mlcanopy_inst, ObuFunc, obu0, obu1, tol)

       ! Check to see if Obukhov length is changing signs between iterations.
       ! If too many changes in sign, set it to a near-neutral value.

!      if (obuold(p)*obu(p) < 0._r8) nmozsgn(p) = nmozsgn(p) + 1
!      obuold(p) = obu(p)
!      if (nmozsgn(p) >= 4) then
!         obu(p) = -1000._r8
!         call ObuFunc (p, ic, il, mlcanopy_inst, obu(p), dummy)
!      end if

       ! The roughness sublayer parameterization uses the expression lm / beta
       ! to calculate wind speed and conductances. Use a constrained value for
       ! lm / beta based on maximum value for eta.

       lm = 2._r8 * beta(p)**3 * Lc(p)
!      lm_over_beta = lm / beta(p)

       eta = min (beta(p)/lm*ztop(p), eta_max)
       lm_over_beta = ztop(p) / eta

       ! Wind speed profile

       call WindProfile (p, lm_over_beta, mlcanopy_inst)

       ! Momentum flux profile: profile is defined at zw (similar to vertical fluxes)

       do ic = 1, ncan(p)
          if (zw(p,ic) > ztop(p)) then
             ! above canopy
             mflx(p,ic) = -ustar(p)**2
          else
             ! within canopy
             mflx(p,ic) = -(ustar(p)**2) * exp(2._r8*(zw(p,ic) - ztop(p)) / lm_over_beta)
          end if
       end do

       ! Aerodynamic conductances

       call AerodynamicConductance (p, lm_over_beta, mlcanopy_inst)

       ! Scalar source/sink fluxes for leaves/soil and concentration profiles.
       ! This uses an implicit solution for temperature and vapor pressure.
    
        call FluxProfileSolution (p, mlcanopy_inst)

       ! No profile for CO2

       do ic = 1, ncan(p)
          cair(p,ic) = co2ref(p)
       end do

       ! Temperature and vapor pressure at top of canopy

       taf(p) = tair(p,ntop(p))
       qaf(p) = mmh2o/mmdry * eair(p,ntop(p)) / (pref(p) - (1._r8-mmh2o/mmdry) * eair(p,ntop(p)))

    end do

    end associate
  end subroutine HF2008

  !-----------------------------------------------------------------------
  subroutine ObuFunc (p, ic, il, mlcanopy_inst, obu_val, obu_dif)
    !
    ! !DESCRIPTION:
    ! Solve for the Obukhov length. For the current estimate of the Obukhov length
    ! (obu_val), calculate u*, T*, and q* and then the new length (obu). The subroutine
    ! returns the change in Obukhov length (obu_dif), which equals zero when the
    ! Obukhov length does not change value between iterations.
    !
    ! !USES:
    use clm_varcon, only : grav, vkc
    use MLclm_varctl, only : sparse_canopy_type
    use MLclm_varcon, only : beta_neutral_max, cr, z0mg, zeta_min, zeta_max
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: p               ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)  :: ic              ! Aboveground layer index
    integer, intent(in)  :: il              ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: obu_val         ! Input value for Obukhov length (m)
    real(r8), intent(out) :: obu_dif        ! Difference in Obukhov length (m)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: obu_cur                     ! Current value for Obukhov length
    real(r8) :: obu_new                     ! New value for Obukhov length
    real(r8) :: c1                          ! Parameter to calculate beta_neutral
    real(r8) :: beta_neutral                ! Neutral value for beta = u*/u(h)
    real(r8) :: h_minus_d                   ! Displacement height below canopy top (ztop - zdisp)
    real(r8) :: zeta                        ! Monin-Obukhov stability parameter (z-d)/L
    real(r8) :: psim, psic                  ! psi functions for momentum and scalars
    real(r8) :: zlog                        ! log((zref-zdisp)/(ztop-zdisp))
    real(r8) :: tstar                       ! Temperature scale (K)
    real(r8) :: qstar                       ! Water vapor scale (kg/kg)
    real(r8) :: tvstar                      ! Virtual potential temperature scale (K)
    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    zref      => mlcanopy_inst%zref_forcing    , &  ! Atmospheric reference height (m)
    uref      => mlcanopy_inst%uref_forcing    , &  ! Wind speed at reference height (m/s)
    thref     => mlcanopy_inst%thref_forcing   , &  ! Atmospheric potential temperature at reference height (K)
    thvref    => mlcanopy_inst%thvref_forcing  , &  ! Atmospheric virtual potential temperature at reference height (K)
    qref      => mlcanopy_inst%qref_forcing    , &  ! Specific humidity at reference height (kg/kg)
    rhomol    => mlcanopy_inst%rhomol_forcing  , &  ! Molar density at reference height (mol/m3)
    ztop      => mlcanopy_inst%ztop_canopy     , &  ! Canopy foliage top height (m)
    lai       => mlcanopy_inst%lai_canopy      , &  ! Leaf area index of canopy (m2/m2)
    sai       => mlcanopy_inst%sai_canopy      , &  ! Stem area index of canopy (m2/m2)
    Lc        => mlcanopy_inst%Lc_canopy       , &  ! Canopy density length scale (m)
    taf       => mlcanopy_inst%taf_canopy      , &  ! Air temperature at canopy top (K)
    qaf       => mlcanopy_inst%qaf_canopy      , &  ! Specific humidity at canopy top (kg/kg)
                                                    ! *** Output ***
    zdisp     => mlcanopy_inst%zdisp_canopy    , &  ! Displacement height (m)
    beta      => mlcanopy_inst%beta_canopy     , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy     , &  ! Prandtl (Schmidt) number at canopy top (-)
    ustar     => mlcanopy_inst%ustar_canopy    , &  ! Friction velocity (m/s)
    gac_to_hc => mlcanopy_inst%gac_to_hc_canopy, &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    obu       => mlcanopy_inst%obu_canopy        &  ! Obukhov length (m)
    )

    ! Use this current value of Obukhov length

    obu_cur = obu_val

    ! Prevent near-zero value of Obukhov length

    if (abs(obu_cur) <= 0.1_r8) obu_cur = 0.1_r8

    ! Determine neutral value of beta

    c1 = (vkc / log((ztop(p) + z0mg)/z0mg))**2
    beta_neutral = min(sqrt(c1 + cr*(lai(p)+sai(p))), beta_neutral_max)

    ! Calculate beta = u* / u(h) for current Obukhov length

    call GetBeta (beta_neutral, Lc(p)/obu_cur, beta(p))

    ! Displacement height, and then adjust for canopy sparseness

    h_minus_d = beta(p)**2 * Lc(p)
    select case (sparse_canopy_type)
    case (1)
       h_minus_d = h_minus_d * (1._r8 - exp(-0.25_r8*(lai(p)+sai(p))/beta(p)**2))
    end select
    h_minus_d = min(ztop(p), h_minus_d)
    zdisp(p) = ztop(p) - h_minus_d

    if ((zref(p) - zdisp(p)) < 0._r8) then
       call endrun (msg=' ERROR: ObuFunc: zdisp height > zref')
    end if

    ! Calculate Prandlt number (Pr) and Schmidt number (Sc) at canopy
    ! top for current Obukhov length. Only need Pr because Sc = Pr.

    call GetPrSc (beta_neutral, beta_neutral_max, Lc(p)/obu_cur, PrSc(p))

    ! Calculate the stability functions psi for momentum and scalars. The
    ! returned function values (psim, psic) are the Monin-Obukhov psi functions
    ! and additionally include the roughness sublayer psi_hat functions.
    ! These are evaluated at the reference height and at the canopy height.
    ! Here, limit Obukhov length based on values of zeta so that extreme
    ! cases are excluded.

    zeta = (zref(p) - zdisp(p)) / obu_cur
    if (zeta >= 0._r8) then
       zeta = min(zeta_max, max(zeta,0.01_r8))
    else
       zeta = max(zeta_min, min(zeta,-0.01_r8))
    end if
    obu_cur = (zref(p) - zdisp(p)) / zeta

    call GetPsiRSL (zref(p), ztop(p), zdisp(p), obu_cur, beta(p), PrSc(p), psim, psic)

    ! Friction velocity

    zlog = log((zref(p)-zdisp(p)) / (ztop(p)-zdisp(p)))
    ustar(p) = uref(p) * vkc / (zlog + psim)

    ! Temperature scale

    tstar = (thref(p) - taf(p)) * vkc / (zlog + psic)

    ! Water vapor scale - use units of specific humidity (kg/kg)

    qstar = (qref(p) - qaf(p)) * vkc / (zlog + psic)

    ! Aerodynamic conductance to canopy height

    gac_to_hc(p) = rhomol(p) * vkc * ustar(p) / (zlog + psic)

    ! Save value for obu used to calculate ustar

    obu(p) = obu_cur

    ! Calculate new Obukhov length (m)

    tvstar = tstar + 0.61_r8 * thref(p) * qstar
    obu_new = ustar(p)**2 * thvref(p) / (vkc * grav * tvstar)

    ! Change in Obukhov length (m)

    obu_dif = obu_new - obu_val

    end associate
  end subroutine ObuFunc

  !-----------------------------------------------------------------------
  subroutine GetBeta (beta_neutral, LcL, beta)
    !
    ! !DESCRIPTION:
    ! Calculate beta = u* / u(h) for current Obukhov length
    !
    ! !USES:
    use MLclm_varcon, only : beta_min, beta_max
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: beta_neutral    ! Neutral value for beta = u*/u(h)
    real(r8), intent(in)  :: LcL             ! Canopy density scale (Lc) / Obukhov length (obu)
    real(r8), intent(out) :: beta            ! Value of u*/u(h) at canopy top
    !
    ! !LOCAL VARIABLES:
    real(r8) :: aa, bb, cc, dd, qq, rr       ! Terms for quadratic or cubic solution
    real(r8) :: LcL_val                      ! LcL with limits applied
    real(r8) :: y, fy, err                   ! Used for error check
    !---------------------------------------------------------------------

    ! Could place limits on LcL, but not used here. Instead, limits
    ! are placed on beta.

    LcL_val = LcL
    
    if (LcL_val <= 0._r8) then

       ! Unstable case: quadratic equation for beta^2 at LcL_val

       aa = 1._r8
       bb = 16._r8 * LcL_val * beta_neutral**4     
       cc = -beta_neutral**4
       beta = sqrt((-bb + sqrt(bb**2 - 4._r8 * aa * cc)) / (2._r8 * aa))

    else

       ! Stable case: cubic equation for beta at LcL_val

       aa = 5._r8 * LcL_val
       bb = 0._r8
       cc = 1._r8
       dd = -beta_neutral
       qq = (2._r8*bb**3 - 9._r8*aa*bb*cc + 27._r8*(aa**2)*dd)**2 - 4._r8*(bb**2 - 3._r8*aa*cc)**3
       qq = sqrt(qq)
       rr = 0.5_r8 * (qq + 2._r8*bb**3 - 9._r8*aa*bb*cc + 27._r8*(aa**2)*dd)
       rr = rr**(1._r8/3._r8)
       beta = -(bb+rr)/(3._r8*aa) - (bb**2 - 3._r8*aa*cc)/(3._r8*aa*rr)    

    end if

    ! Error check

    y = LcL_val * beta**2 
    fy = phim_monin_obukhov(y)
    err = beta * fy - beta_neutral
    if (abs(err) > 1.e-06_r8) then
       call endrun (msg=' ERROR: GetBeta: beta error')
    end if

    ! Place limits on beta

    beta = min(beta_max, max(beta,beta_min))

  end subroutine GetBeta

  !-----------------------------------------------------------------------
  subroutine GetPrSc (beta_neutral, beta_neutral_max, LcL, PrSc)
    !
    ! !DESCRIPTION:
    ! Calculate Prandlt number (Pr) and Schmidt number (Sc) at canopy
    ! top for current Obukhov length
    !
    ! !USES:
    use MLclm_varctl, only : sparse_canopy_type
    use MLclm_varcon, only : Pr0, Pr1, Pr2
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: beta_neutral     ! Neutral value for beta = u*/u(h)
    real(r8), intent(in)  :: beta_neutral_max ! Maximum value for beta in neutral conditions
    real(r8), intent(in)  :: LcL              ! Canopy density scale (Lc) / Obukhov length (obu)
    real(r8), intent(out) :: PrSc             ! Prandtl (Schmidt) number at canopy top
    !---------------------------------------------------------------------
    
    PrSc = Pr0 + Pr1 * tanh(Pr2*LcL)

    ! Adjust for canopy sparseness

    select case (sparse_canopy_type)
    case (1)
       PrSc = (1._r8 - beta_neutral/beta_neutral_max) * 1._r8 + (beta_neutral/beta_neutral_max) * PrSc
    end select

  end subroutine GetPrSc

  !-----------------------------------------------------------------------
  subroutine GetPsiRSL (za, hc, disp, obu, beta, PrSc, psim, psic)
    !
    ! !DESCRIPTION:
    ! Calculate the stability functions psi for momentum and scalars. The
    ! returned function values (psim, psic) are the Monin-Obukhov psi functions
    ! and additionally include the roughness sublayer psihat functions.
    ! These are evaluated between the height za and the canopy height hc.
    !
    ! !USES:
    use clm_varcon, only : vkc
    use MLclm_varcon, only : c2, dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: za        ! Atmospheric height (m)
    real(r8), intent(in)  :: hc        ! Canopy height (m)
    real(r8), intent(in)  :: disp      ! Displacement height (m)
    real(r8), intent(in)  :: obu       ! Obukhov length (m)
    real(r8), intent(in)  :: beta      ! Value of u* / u at canopy top
    real(r8), intent(in)  :: PrSc      ! Prandtl (Schmidt) number at canopy top
    real(r8), intent(out) :: psim      ! psi function for momentum including RSL influence
    real(r8), intent(out) :: psic      ! psi function for scalars  including RSL influence
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dt                     ! Displacement height below canopy top (dt = ztop - zdisp)
    real(r8) :: phim                   ! Monin-Obukhov phi function for momentum at canopy top
    real(r8) :: phic                   ! Monin-Obukhov phi function for scalars at canopy top
    real(r8) :: c1                     ! RSL magnitude multiplier
    real(r8) :: psihat1                ! RSL psihat function evaluated at za
    real(r8) :: psihat2                ! RSL psihat function evaluated at hc
    real(r8) :: psi1                   ! Monin-Obukhov psi function evaluated at za
    real(r8) :: psi2                   ! Monin-Obukhov psi function evaluated at hc
    !---------------------------------------------------------------------

    dt = hc - disp

    ! In the RSL theory, c1 and c2 represent the scaled magnitude and
    ! height over which the RSL theory modifies MOST via the psihat functions:
    !
    !
    !     z
    !     ^
    !     |                                    . ___
    !     |                                    .  ^
    !     |                                   .   |
    !     |                                  .    |
    !     |                                 .     |
    !     |                               .       |
    !     |                             .        c2
    !     |                          .            |
    !     |                       .               |
    !     |                  .                    |
    !     |           .                           |
    !     |   .                                 _\/_
    !     -------------------------------------------> u
    !         |<-------------- c1 --------------->|

    phim = phim_monin_obukhov((hc-disp)/obu)
    c1 = (1._r8 - vkc / (2._r8 * beta * phim)) * exp(0.5_r8*c2)

    ! Evaluate the roughness sublayer psihat function for momentum at
    ! the height za and at the canopy height hc. Values for psihat are obtained
    ! from a look-up table. Here, heights are above the canopy for compatibility
    ! with the supplied look-up table. These heights are also scaled by dt = hc-d
    ! so that the look-up table uses (za-hc)/dt and (hc-hc)/dt. Also the term
    ! dt in the integration of psihat is scaled by the Obukhov length L (dt/obu).
    ! This means that the returned psihat value needs to be scaled (multiplied) by
    ! c1 before it fully represents psihat as it appears in the RSL equations.

    call LookupPsihat ((za-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM, psihat1)
    call LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM, psihat2)
    psihat1 = psihat1 * c1
    psihat2 = psihat2 * c1

    ! Evaluate the Monin-Obukhov psi function for momentum at the height za
    ! and at the canopy height hc

    psi1 = psim_monin_obukhov((za-disp)/obu)
    psi2 = psim_monin_obukhov((hc-disp)/obu)

    ! psi function for momentum including RSL influence

    psim = -psi1 + psi2 + psihat1 - psihat2 + vkc / beta

    ! Now do the same for scalars

    phic = phic_monin_obukhov((hc-disp)/obu)
    c1 = (1._r8 - PrSc*vkc / (2._r8 * beta * phic)) * exp(0.5_r8*c2)

    call LookupPsihat ((za-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH, psihat1)
    call LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH, psihat2)
    psihat1 = psihat1 * c1
    psihat2 = psihat2 * c1

    psi1 = psic_monin_obukhov((za-disp)/obu)
    psi2 = psic_monin_obukhov((hc-disp)/obu)

    psic = -psi1 + psi2 + psihat1 - psihat2

  end subroutine GetPsiRSL

  !-----------------------------------------------------------------------
  function phim_monin_obukhov (zeta) result(phi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov phi stability function for momentum
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: phi               ! phi for momentum
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then           ! --- unstable
       phi = 1._r8 / sqrt(sqrt(1._r8 - 16._r8 * zeta))
    else                             ! --- stable
       phi = 1._r8 + 5._r8 * zeta
    end if

  end function phim_monin_obukhov

  !-----------------------------------------------------------------------
  function phic_monin_obukhov (zeta) result(phi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov phi stability function for scalars
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: phi               ! phi for scalars
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then           ! --- unstable
       phi = 1._r8 / sqrt(1._r8 - 16._r8 * zeta)
    else                             ! --- stable
       phi = 1._r8 + 5._r8 * zeta
    end if

  end function phic_monin_obukhov

  !-----------------------------------------------------------------------
  function psim_monin_obukhov (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov psi stability function for momentum
    !
    ! !USES:
    use clm_varcon, only : pi => rpi
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x                 ! (1 - 16*zeta)**1/4
    real(r8) :: psi               ! psi for momentum
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then           ! --- unstable
       x = sqrt(sqrt(1._r8 - 16._r8 * zeta))
       psi = 2._r8 * log((1._r8+x)/2._r8) + log((1._r8+x*x)/2._r8) - 2._r8*atan(x) + pi * 0.5_r8
    else                             ! --- stable
       psi = -5._r8 * zeta
    end if

  end function psim_monin_obukhov

  !-----------------------------------------------------------------------
  function psic_monin_obukhov (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov psi stability function for scalars
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x                 ! (1 - 16*zeta)**1/4
    real(r8) :: psi               ! psi for scalars
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then           ! --- unstable
       x = sqrt(sqrt(1._r8 - 16._r8 * zeta))
       psi = 2._r8 * log((1._r8+x*x)/2._r8)
    else                             ! --- stable
       psi = -5._r8 * zeta
    end if

  end function psic_monin_obukhov

  !-----------------------------------------------------------------------
  subroutine LookupPsihat (zdt, dtL, zdtgrid, dtLgrid, psigrid, psihat)
    !
    ! !DESCRIPTION:
    ! Determine the RSL function psihat as provided through a look-up table
    ! for input values of zdt and dtL. Linearly interpolate between values
    ! supplied on the look-up table grid defined by zdtgrid, dtLgrid, psigrid.
    !
    ! NOTE: The psihat presented in Harman and Finnigan (2007,2008) and Harman (2012)
    ! has been re-written in non-dimensional form such that it now appears as:
    !
    ! psihat(z) = c1 * A(z/(beta^2*Lc),(beta^2*Lc)/L)
    !
    ! This routine gets the value of A from a look-up table. Noting that dt=beta^2*Lc,
    ! this routine therefore requires values of z/dt and dt/L. In addition, this means
    ! that the returned psihat value needs to be scaled (multiplied) by c1 before it fully
    ! represents psihat as it appears in the RSL equations.
    !
    ! !USES:
    use MLclm_varcon, only : nZ, nL
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zdt            ! Height (above canopy) normalized by dt
    real(r8), intent(in) :: dtL            ! dt/L (displacement height/Obukhov length)
    real(r8), intent(in) :: zdtgrid(nZ,1)  ! Grid of zdt on which psihat is given
    real(r8), intent(in) :: dtLgrid(1,nL)  ! Grid of dtL on which psihat is given
    real(r8), intent(in) :: psigrid(nZ,nL) ! Grid of psihat values
    real(r8), intent(out):: psihat         ! Value of psihat
    !
    !LOCAL VARIABLES
    integer  :: ii, jj                     ! Looping indices
    integer  :: L1, L2, Z1, Z2             ! Grid indices for psihat sought
    real(r8) :: wL1, wL2, wZ1, wZ2         ! Weights for averaging
    !---------------------------------------------------------------------

    ! Find indices and weights for dtL values which bracket the specified dtL

    L1 = 0 ; L2 = 0
    if (dtL <= dtLgrid(1,1)) then
       L1 = 1
       L2 = 1
       wL1 = 0.5_r8
       wL2 = 0.5_r8
    else if (dtL > dtLgrid(1,nL)) then
       L1 = nL
       L2 = nL
       wL1 = 0.5_r8
       wL2 = 0.5_r8
    else
       do jj = 1, nL-1
          if ((dtL <= dtLgrid(1,jj+1)) .and. (dtL > dtLgrid(1,jj))) then
             L1 = jj
             L2 = jj + 1
             wL1 = (dtLgrid(1,L2) - dtL) / (dtLgrid(1,L2) - dtLgrid(1,L1))
             wL2 = 1._r8 - wL1
          end if
       end do
    end if

    if (L1 == 0 .or. L2 == 0) then
       call endrun (msg=' ERROR: LookupPsihat error: indices L1 and L2 not found')
    end if

    ! Find indices and weights for zdt values which bracket the specified zdt

    Z1 = 0 ; Z2 = 0
    if (zdt > zdtgrid(1,1)) then
       Z1 = 1
       Z2 = 1
       wZ1 = 0.5_r8
       wZ2 = 0.5_r8
    else if (zdt < zdtgrid(nZ,1)) then
       Z1 = nZ
       Z2 = nZ
       wZ1 = 0.5_r8
       wZ2 = 0.5_r8
    else
       do ii = 1, nZ-1
          if ((zdt >= zdtgrid(ii+1,1)) .and. (zdt < zdtgrid(ii,1))) then
             Z1 = ii
             Z2 = ii + 1
             wZ1 = (zdt - zdtgrid(ii+1,1)) / (zdtgrid(ii,1) - zdtgrid(ii+1,1))
             wZ2 = 1._r8 - wZ1
          end if
       end do
    end if

    if (Z1 == 0 .or. Z2 == 0) then
       call endrun (msg=' ERROR: LookupPsihat error: indices Z1 and Z2 not found')
    end if

    ! Calculate psihat as a weighted average of the values of psihat on the grid

    psihat = wZ1 * wL1 * psigrid(Z1,L1) + wZ2 * wL1 * psigrid(Z2,L1) &
           + wZ1 * wL2 * psigrid(Z1,L2) + wZ2 * wL2 * psigrid(Z2,L2)

  end subroutine LookupPsihat

  !-----------------------------------------------------------------------
  subroutine WindProfile (p, lm_over_beta, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Wind speed profile above and within canopy
    !
    ! !USES:
    use clm_varcon, only : vkc
    use MLclm_varcon, only : wind_min
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p              ! Patch index for CLM g/l/c/p hierarchy
    real(r8), intent(in) :: lm_over_beta  ! lm / beta
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                        ! Aboveground layer index
    real(r8) :: psim                      ! psi function for momentum
    real(r8) :: psic                      ! psi function for scalars
    real(r8) :: zlog_m                    ! log height
    real(r8) :: ave                       ! Average wind speed in canopy (m/s)
    !---------------------------------------------------------------------

    associate ( &
                                                 ! *** Input ***
    ncan      => mlcanopy_inst%ncan_canopy  , &  ! Number of aboveground layers
    ntop      => mlcanopy_inst%ntop_canopy  , &  ! Index for top leaf layer
    zs        => mlcanopy_inst%zs_profile   , &  ! Canopy layer height for scalar concentration and source (m)
    dz        => mlcanopy_inst%dz_profile   , &  ! Canopy layer thickness (m)
    ztop      => mlcanopy_inst%ztop_canopy  , &  ! Canopy foliage top height (m)
    zdisp     => mlcanopy_inst%zdisp_canopy , &  ! Displacement height (m)
    obu       => mlcanopy_inst%obu_canopy   , &  ! Obukhov length (m)
    beta      => mlcanopy_inst%beta_canopy  , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy  , &  ! Prandtl (Schmidt) number at canopy top (-)
    ustar     => mlcanopy_inst%ustar_canopy , &  ! Friction velocity (m/s)
                                                 ! *** Output ***
    uaf      => mlcanopy_inst%uaf_canopy    , &  ! Wind speed at canopy top (m/s)
    wind     => mlcanopy_inst%wind_profile    &  ! Canopy layer wind speed (m/s)
    )

    ! Above-canopy wind profile: wind speed is defined at zs

    do ic = ntop(p)+1, ncan(p)
       call GetPsiRSL (zs(p,ic), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim, psic)
       zlog_m = log((zs(p,ic)-zdisp(p)) / (ztop(p)-zdisp(p)))
       wind(p,ic) = ustar(p) / vkc * (zlog_m + psim)
    end do

    ! Wind speed at top of canopy

    uaf(p) = ustar(p) / beta(p)

    ! Within-canopy wind profile: wind speed is defined at zs

    do ic = 1, ntop(p)
       wind(p,ic) = uaf(p) * exp((zs(p,ic) - ztop(p)) / lm_over_beta)
       wind(p,ic) = max(wind(p,ic), wind_min)
    end do

    ! Average wind in the canopy, weighted by layer thickness

    ave = sum(wind(p,1:ntop(p))*dz(p,1:ntop(p)))/ sum(dz(p,1:ntop(p)))

    ! Replace wind speed at each layer with average wind in the canopy
    ! or with friction velocity

!   wind(p,1:ntop(p)) = ave
!   wind(p,1:ntop(p)) = ustar(p)

    end associate
  end subroutine WindProfile

  !-----------------------------------------------------------------------
  subroutine AerodynamicConductance (p, lm_over_beta, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Aerodynamic conductances above and within the canopy.
    ! Conductances are defined between zs(i) and zs(i+1).
    !
    ! !USES:
    use clm_varcon, only : vkc
    use MLclm_varctl, only : HF_extension_type
    use MLclm_varcon, only : z0mg, ra_max
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p              ! Patch index for CLM g/l/c/p hierarchy
    real(r8), intent(in) :: lm_over_beta  ! lm / beta
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                        ! Aboveground layer index
    real(r8) :: psim1, psim2              ! psi function for momentum
    real(r8) :: psic, psic1, psic2        ! psi function for scalars
    real(r8) :: zlog_m, zlog_c            ! log height
    real(r8) :: zu, zl                    ! Upper and lower heights for within canopy resistances (m)
    real(r8) :: res                       ! Resistance (s/m)
    real(r8) :: ustar_g                   ! Friction velocity at ground (m/s)
    real(r8) :: z0cg                      ! Roughness length of ground for scalars (m)
    real(r8) :: sumres                    ! Sum of aerodynamic resistances above canopy
    real(r8) :: gac_above_foliage         ! Aerodynamic conductance for top/bottom foliage layer
    real(r8) :: gac_below_foliage         ! Aerodynamic conductance for top/bottom foliage layer
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *** Input ***
    zref      => mlcanopy_inst%zref_forcing     , &  ! Atmospheric reference height (m)
    rhomol    => mlcanopy_inst%rhomol_forcing   , &  ! Molar density at reference height (mol/m3)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    ntop      => mlcanopy_inst%ntop_canopy      , &  ! Index for top leaf layer
    ztop      => mlcanopy_inst%ztop_canopy      , &  ! Canopy foliage top height (m)
    zdisp     => mlcanopy_inst%zdisp_canopy     , &  ! Displacement height (m)
    obu       => mlcanopy_inst%obu_canopy       , &  ! Obukhov length (m)
    beta      => mlcanopy_inst%beta_canopy      , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy      , &  ! Prandtl (Schmidt) number at canopy top (-)
    ustar     => mlcanopy_inst%ustar_canopy     , &  ! Friction velocity (m/s)
    gac_to_hc => mlcanopy_inst%gac_to_hc_canopy , &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    zs        => mlcanopy_inst%zs_profile       , &  ! Canopy layer height for scalar concentration and source (m)
    wind      => mlcanopy_inst%wind_profile     , &  ! Canopy layer wind speed (m/s)
                                                     ! *** Output ***
    gac0      => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    gac       => mlcanopy_inst%gac_profile        &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    )

    ! -------------------------------------
    ! Above-canopy aerodynamic conductances 
    ! conductance: mol/m3 * m/s = mol/m2/s
    ! -------------------------------------

    do ic = ntop(p)+1, ncan(p)-1
       call GetPsiRSL (zs(p,ic),   ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim1, psic1)
       call GetPsiRSL (zs(p,ic+1), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim2, psic2)
       ! equivalent to: -psi_c_z2 + psi_c_z1 + psi_c_rsl_z2 - psi_c_rsl_z1
       psic = psic2 - psic1
       zlog_c = log((zs(p,ic+1)-zdisp(p)) / (zs(p,ic)-zdisp(p)))
       gac(p,ic) = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)
    end do

    ! Special case for the top layer to the reference height

    ic = ncan(p)
    call GetPsiRSL (zs(p,ic), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim1, psic1)
    call GetPsiRSL (zref(p),  ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim2, psic2)
    psic = psic2 - psic1
    zlog_c = log((zref(p)-zdisp(p)) / (zs(p,ic)-zdisp(p)))
    gac(p,ic) = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)

    ! The top foliage layer includes terms for within and above the canopy. Here,
    ! calculate the conductance from top of foliage at height ztop to the layer
    ! immediately above it at height zs(ntop+1)

    ic = ntop(p)
    call GetPsiRSL (ztop(p),    ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim1, psic1)
    call GetPsiRSL (zs(p,ic+1), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim2, psic2)
    psic = psic2 - psic1
    zlog_c = log((zs(p,ic+1)-zdisp(p)) / (ztop(p)-zdisp(p)))
    gac_above_foliage = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)

    ! Make sure above-canopy aerodynamic resistances sum to 1/gac_to_hc

    sumres = 1._r8 / gac_above_foliage
    do ic = ntop(p)+1, ncan(p)
       sumres = sumres + 1._r8 / gac(p,ic)
    end do

    if (abs(1._r8/sumres - gac_to_hc(p)) > 1.e-06_r8) then
       call endrun (msg=' ERROR: AerodynamicConductance: above-canopy aerodynamic conductance error')
    end if

    ! --------------------------------------
    ! Within-canopy aerodynamic conductances
    ! res = resistance: s/m
    ! gac = conductance: (mol/m3) / (s/m) = mol/m2/s
    ! --------------------------------------

    do ic = 1, ntop(p)-1
       zl = zs(p,ic) - ztop(p)
       zu = zs(p,ic+1) - ztop(p)
       res = PrSc(p) / (beta(p) * ustar(p)) * (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
       gac(p,ic) = rhomol(p) / res
    end do

    ! Special case for top foliage layer: conductance from zs to ztop ...

    ic = ntop(p)
    zl = zs(p,ic) - ztop(p)
    zu = ztop(p) - ztop(p)
    res = PrSc(p) / (beta(p) * ustar(p)) * (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
    gac_below_foliage = rhomol(p) / res

    ! ... and now include additional conductance from top of foliage to next layer above

    gac(p,ic) = 1._r8 / (1._r8 / gac_below_foliage + 1._r8 / gac_above_foliage)

    ! Aerodynamic conductance at ground

    z0cg = 0.1_r8 * z0mg
    if (z0mg > zs(p,1) .or. z0cg > zs(p,1)) then
       call endrun (msg=' ERROR: AerodynamicConductance: soil roughness error')
    end if

    select case (HF_extension_type)

    case (1)

       ! Extend HF exponential profile to ground (taken as z0cg)

       zl = z0cg - ztop(p)
       zu = zs(p,1) - ztop(p)
       res = PrSc(p) / (beta(p) * ustar(p)) * (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
       gac0(p) = rhomol(p) / res

    case (2)

       ! Use log profile to ground (taken as z0mg, where wind speed is zero)

       zlog_m = log(zs(p,1)/z0mg)
       ustar_g = wind(p,1) * vkc / zlog_m
       gac0(p) = rhomol(p) * vkc * ustar_g / zlog_m

!      zlog_m = log(zs(p,1)/z0mg)                     !!! CLMml v0 CODE !!!
!      ustar_g = wind(p,1) * vkc / zlog_m             !!! CLMml v0 CODE !!!
!      ustar_g = max(ustar_g, 0.01_r8)                !!! CLMml v0 CODE !!!
!      z0cg = 0.1_r8 * z0mg                           !!! CLMml v0 CODE !!!
!      zlog_c = log(zs(p,1)/z0cg)                     !!! CLMml v0 CODE !!!
!      gac0(p) = rhomol(p) * vkc * ustar_g / zlog_c   !!! CLMml v0 CODE !!!

    end select

    ! Limit resistances to < 500 s/m

    res = min (rhomol(p)/gac0(p), ra_max)
    gac0(p) = rhomol(p) / res

    do ic = 1, ncan(p)
       res = min (rhomol(p)/gac(p,ic), ra_max)
       gac(p,ic) = rhomol(p) / res
    end do

    end associate
  end subroutine AerodynamicConductance

  !-----------------------------------------------------------------------
  subroutine FluxProfileSolution (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute scalar source/sink fluxes for leaves/soil and concentration
    ! profiles. This uses an implicit solution for temperature and vapor
    ! pressure. The boundary condition is the above-canopy scalar values at the
    ! reference height. Vegetation and ground temperatures and fluxes are
    ! calculated as part of the implicit solution.
    !
    ! !USES:
    use MLclm_varctl, only : dtime_substep
    use MLclm_varpar, only : isun, isha, nlevmlcan, nleaf
    use MLLeafFluxesMod, only : LeafFluxes
    use MLMathToolsMod, only: tridiag_2eq
    use MLSoilFluxesMod, only : SoilFluxes
    use MLWaterVaporMod, only : SatVap, LatVap
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p                  ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                            ! Aboveground layer index
    integer  :: il                            ! Sunlit (1) or shaded (2) leaf index
    real(r8) :: dtime                         ! Model time step (s)
    real(r8) :: lambda                        ! Latent heat of vaporization (J/mol)
    real(r8) :: esat                          ! Saturation vapor pressure (Pa)
    real(r8) :: desat                         ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(r8) :: qsat                          ! Saturation vapor pressure (mol/mol)
    real(r8) :: den                           ! Intermediate calculation
    real(r8) :: pai                           ! Plant area index of layer (m2/m2)
    real(r8) :: rho_dz_over_dt                ! Intermediate calculation for canopy air storage
    real(r8) :: gac_ic_minus_one              ! Special case for ic-1 = 0: use soil conductance not ic-1

    real(r8) :: qsat0                         ! Saturation vapor pressure at ground (mol/mol)
    real(r8) :: dqsat0                        ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: gsw                           ! Soil conductance for water vapor (mol H2O/m2/s)
    real(r8) :: gs0                           ! Total soil-to-air conductance for water vapor (mol H2O/m2/s)
    real(r8) :: alpha0                        ! Coefficient for ground temperature (dimensionless)
    real(r8) :: beta0                         ! Coefficient for ground temperature (K mol/mol)
    real(r8) :: delta0                        ! Coefficient for ground temperature (K)
    real(r8) :: c01                           ! Soil heat flux term (W/m2)
    real(r8) :: c02                           ! Soil heat flux term (W/m2/K)
    real(r8) :: t0                            ! Soil surface temperature (K)
    real(r8) :: sh0                           ! Ground sensible heat flux (W/m2)
    real(r8) :: et0                           ! Ground evaporation flux (mol H2O/m2/s)
    real(r8) :: g0                            ! Soil heat flux (W/m2)
    real(r8) :: e0                            ! Soil surface vapor pressure (Pa)

    real(r8) :: tleaf_implic(nlevmlcan,nleaf) ! Leaf temperature from implicit solution (K)
    real(r8) :: gleaf_sh(nlevmlcan,nleaf)     ! Leaf conductance for sensible heat (mol/m2/s)
    real(r8) :: gleaf_et(nlevmlcan,nleaf)     ! Leaf conductance for water vapor (mol/m2/s)
    real(r8) :: heatcap(nlevmlcan,nleaf)      ! Heat capacity of leaves (J/m2/K)
    real(r8) :: avail_energy(nlevmlcan,nleaf) ! Available energy for leaf (W/m2)
    real(r8) :: dqsat(nlevmlcan,nleaf)        ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: qsat_term(nlevmlcan,nleaf)    ! Intermediate calculation for saturation vapor pressure (mol/mol)
    real(r8) :: alpha(nlevmlcan,nleaf)        ! Coefficient for leaf temperature (dimensionless)
    real(r8) :: beta(nlevmlcan,nleaf)         ! Coefficient for leaf temperature (K mol/mol)
    real(r8) :: delta(nlevmlcan,nleaf)        ! Coefficient for leaf temperature (K)

    real(r8) :: a1(nlevmlcan)                 ! Coefficient for canopy air temperature
    real(r8) :: b11(nlevmlcan)                ! Coefficient for canopy air temperature
    real(r8) :: b12(nlevmlcan)                ! Coefficient for canopy air temperature
    real(r8) :: c1(nlevmlcan)                 ! Coefficient for canopy air temperature
    real(r8) :: d1(nlevmlcan)                 ! Coefficient for canopy air temperature

    real(r8) :: a2(nlevmlcan)                 ! Coefficient for canopy air water vapor mole fraction
    real(r8) :: b21(nlevmlcan)                ! Coefficient for canopy air water vapor mole fraction
    real(r8) :: b22(nlevmlcan)                ! Coefficient for canopy air water vapor mole fraction
    real(r8) :: c2(nlevmlcan)                 ! Coefficient for canopy air water vapor mole fraction
    real(r8) :: d2(nlevmlcan)                 ! Coefficient for canopy air water vapor mole fraction

    ! These are needed only for error checks
    real(r8) :: storage_sh(nlevmlcan)         ! Heat storage flux in air (W/m2)
    real(r8) :: storage_et(nlevmlcan)         ! Water vapor storage flux in air (mol H2O/m2/s)
    real(r8) :: stveg(nlevmlcan)              ! Canopy layer leaf storage heat flux (W/m2)
    real(r8) :: shsrc(nlevmlcan)              ! Canopy layer leaf sensible heat flux (W/m2)
    real(r8) :: etsrc(nlevmlcan)              ! Canopy layer leaf water vapor flux (mol H2O/m2/s)
    real(r8) :: stveg_leaf                    ! Canopy layer leaf storage heat flux from LeafFluxes (W/m2)
    real(r8) :: shsrc_leaf                    ! Canopy layer leaf sensible heat flux from LeafFluxes (W/m2)
    real(r8) :: etsrc_leaf                    ! Canopy layer leaf water vapor flux from LeafFluxes (mol H2O/m2/s)
    real(r8) :: err                           ! Energy imbalance (W/m2)
    real(r8) :: sum_src                       ! Sum of source flux over all layers
    real(r8) :: sum_storage                   ! Sum of storage flux over all layers
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *** Input ***
    tref      => mlcanopy_inst%tref_forcing     , &  ! Air temperature at reference height (K)
    thref     => mlcanopy_inst%thref_forcing    , &  ! Atmospheric potential temperature at reference height (K)
    eref      => mlcanopy_inst%eref_forcing     , &  ! Vapor pressure at reference height (Pa)
    pref      => mlcanopy_inst%pref_forcing     , &  ! Air pressure at reference height (Pa)
    rhomol    => mlcanopy_inst%rhomol_forcing   , &  ! Molar density at reference height (mol/m3)
    cpair     => mlcanopy_inst%cpair_forcing    , &  ! Specific heat of air (constant pressure) at reference height (J/mol/K)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    ntop      => mlcanopy_inst%ntop_canopy      , &  ! Index for top leaf layer
    gac0      => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    tg_bef    => mlcanopy_inst%tg_bef_soil      , &  ! Soil surface temperature for previous timestep (K)
    rhg       => mlcanopy_inst%rhg_soil         , &  ! Relative humidity of airspace at soil surface (fraction)
    rnsoi     => mlcanopy_inst%rnsoi_soil       , &  ! Net radiation: ground (W/m2)
    soilres   => mlcanopy_inst%soilres_soil     , &  ! Soil evaporative resistance (s/m)
    soil_t    => mlcanopy_inst%soil_t_soil      , &  ! Temperature of first snow/soil layer (K)
    soil_dz   => mlcanopy_inst%soil_dz_soil     , &  ! Depth to temperature of first snow/soil layer (m)
    soil_tk   => mlcanopy_inst%soil_tk_soil     , &  ! Thermal conductivity of first snow/soil layer (W/m/K)
    dz        => mlcanopy_inst%dz_profile       , &  ! Canopy layer thickness (m)
    dpai      => mlcanopy_inst%dpai_profile     , &  ! Canopy layer plant area index (m2/m2)
    fwet      => mlcanopy_inst%fwet_profile     , &  ! Canopy layer fraction of plant area index that is wet
    fdry      => mlcanopy_inst%fdry_profile     , &  ! Canopy layer fraction of plant area index that is green and dry
    fracsun   => mlcanopy_inst%fracsun_profile  , &  ! Canopy layer sunlit fraction (-)
    cpleaf    => mlcanopy_inst%cpleaf_profile   , &  ! Canopy layer leaf heat capacity (J/m2 leaf/K)
    gac       => mlcanopy_inst%gac_profile      , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    tair_bef  => mlcanopy_inst%tair_bef_profile , &  ! Canopy layer air temperature for previous timestep (K)
    eair_bef  => mlcanopy_inst%eair_bef_profile , &  ! Canopy layer vapor pressure for previous timestep (Pa)
    gbh       => mlcanopy_inst%gbh_leaf         , &  ! Leaf boundary layer conductance: heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv_leaf         , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gs        => mlcanopy_inst%gs_leaf          , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    rnleaf    => mlcanopy_inst%rnleaf_leaf      , &  ! Leaf net radiation (W/m2 leaf)
    tleaf_bef => mlcanopy_inst%tleaf_bef_leaf   , &  ! Leaf temperature for previous timestep (K)
                                                     ! *** Output ***
    tair      => mlcanopy_inst%tair_profile     , &  ! Canopy layer air temperature (K)
    eair      => mlcanopy_inst%eair_profile     , &  ! Canopy layer vapor pressure (Pa)
    shair     => mlcanopy_inst%shair_profile    , &  ! Canopy layer air sensible heat flux (W/m2)
    etair     => mlcanopy_inst%etair_profile    , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair     => mlcanopy_inst%stair_profile    , &  ! Canopy layer air storage heat flux (W/m2)
                                                     ! *** Output from LeafFluxes
    tleaf     => mlcanopy_inst%tleaf_leaf       , &  ! Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf_leaf      , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf_leaf      , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf    => mlcanopy_inst%lhleaf_leaf      , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf    => mlcanopy_inst%trleaf_leaf      , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf    => mlcanopy_inst%evleaf_leaf      , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
                                                     ! *** Output from SoilFluxes
    shsoi     => mlcanopy_inst%shsoi_soil       , &  ! Sensible heat flux: ground (W/m2)
    lhsoi     => mlcanopy_inst%lhsoi_soil       , &  ! Latent heat flux: ground (W/m2)
    etsoi     => mlcanopy_inst%etsoi_soil       , &  ! Water vapor flux: ground (mol H2O/m2/s)
    gsoi      => mlcanopy_inst%gsoi_soil        , &  ! Soil heat flux (W/m2)
    eg        => mlcanopy_inst%eg_soil          , &  ! Soil surface vapor pressure (Pa)
    tg        => mlcanopy_inst%tg_soil            &  ! Soil surface temperature (K)
    )

    ! Time step (s)

    dtime = dtime_substep

    ! Latent heat of vaporization

    lambda = LatVap(tref(p))

    ! Terms for ground temperature, which is calculated from the energy balance:
    !
    ! Rn0 - H0 - lambda*E0 - G = 0
    !
    ! and is rewritten in relation to T and q as:
    !
    ! T0 = alpha0*T(1) + beta0*q(1) + delta0
    !
    ! by substituting the flux equations for H0, E0, and G

    call SatVap (tg_bef(p), esat, desat)               ! Vapor pressure (Pa) at ground temperature
    qsat0 = esat / pref(p) ; dqsat0 = desat / pref(p)  ! Pa -> mol/mol

    gsw = (1._r8 / soilres(p)) * rhomol(p)             ! Soil conductance for water vapor: s/m -> mol H2O/m2/s
    gs0 = gac0(p) * gsw / (gac0(p) + gsw)              ! Total soil-to-air conductance, including aerodynamic term

    c02 = soil_tk(p) / soil_dz(p)                      ! Soil heat flux term (W/m2/K)
    c01 = -c02 * soil_t(p)                             ! Soil heat flux term (W/m2)

    den = cpair(p) * gac0(p) + lambda * rhg(p) * gs0 * dqsat0 + c02
    alpha0 = cpair(p) * gac0(p) / den
    beta0 = lambda * gs0 / den
    delta0 = (rnsoi(p) - lambda * rhg(p) * gs0 * (qsat0 - dqsat0 * tg_bef(p)) - c01) / den

    ! Similarly, re-arrange the leaf energy balance to calculate leaf
    ! temperature in relation to T and q:
    !
    ! Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
    ! Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)

    do ic = 1, ncan(p)

       ! Calculate terms for sunlit and shaded leaves

       if (dpai(p,ic) > 0._r8) then

          do il = 1, nleaf

             ! Leaf conductances

             gleaf_sh(ic,il) = 2._r8 * gbh(p,ic,il)
             gleaf_et(ic,il) = gs(p,ic,il)*gbv(p,ic,il)/(gs(p,ic,il)+gbv(p,ic,il)) * fdry(p,ic) + gbv(p,ic,il) * fwet(p,ic)

             ! Heat capacity of leaves

             heatcap(ic,il) = cpleaf(p,ic)

             ! Available energy: net radiation

             avail_energy(ic,il) = rnleaf(p,ic,il)

             ! Saturation vapor pressure and derivative for leaf temperature at time n: Pa -> mol/mol

             call SatVap (tleaf_bef(p,ic,il), esat, desat)
             qsat = esat / pref(p) ; dqsat(ic,il) = desat / pref(p)

             ! Term for linearized vapor pressure at leaf temperature:
             ! qsat(tleaf) = qsat(tleaf_bef) + dqsat * (tleaf - tleaf_bef)
             ! Here, qsat_term contains the terms with tleaf_bef

             qsat_term(ic,il) = qsat - dqsat(ic,il) * tleaf_bef(p,ic,il)

             ! alpha, beta, delta coefficients for leaf temperature

             den = heatcap(ic,il) / dtime + gleaf_sh(ic,il) * cpair(p) + gleaf_et(ic,il) * lambda * dqsat(ic,il)
             alpha(ic,il) = gleaf_sh(ic,il) * cpair(p) / den
             beta(ic,il) = gleaf_et(ic,il) * lambda / den
             delta(ic,il) = avail_energy(ic,il) / den &
                          - lambda * gleaf_et(ic,il) * qsat_term(ic,il) / den &
                          + heatcap(ic,il) / dtime * tleaf_bef(p,ic,il) / den

             ! Now scale flux terms for leaf area so that fluxes are for the canopy layer

             if (il == isun) then
                pai = dpai(p,ic) * fracsun(p,ic)
             else if (il == isha) then
                pai = dpai(p,ic) * (1._r8 - fracsun(p,ic))
             end if

             gleaf_sh(ic,il) = gleaf_sh(ic,il) * pai
             gleaf_et(ic,il) = gleaf_et(ic,il) * pai
             heatcap(ic,il) = heatcap(ic,il) * pai
             avail_energy(ic,il) = avail_energy(ic,il) * pai

          end do

       else

          ! Zero out terms

          do il = 1, nleaf
             gleaf_sh(ic,il) = 0._r8
             gleaf_et(ic,il) = 0._r8
             heatcap(ic,il) = 0._r8
             avail_energy(ic,il) = 0._r8
             dqsat(ic,il) = 0._r8
             qsat_term(ic,il) = 0._r8
             alpha(ic,il) = 0._r8
             beta(ic,il) = 0._r8
             delta(ic,il) = 0._r8
          end do

       end if

    end do

    ! The system of equations for air temperature (K) and water vapor (mol/mol)
    ! at each layer is:
    !
    ! a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
    ! a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
    !
    ! These equations are obtained by substituting the equations:
    !
    ! Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
    ! Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
    !
    ! and:
    !
    ! T0 = alpha0*T(1) + beta0*q(1) + delta0
    !
    ! into the one-dimensional scalar conservation equations for T and q.
    !
    ! Calculate the:
    ! a1, b11, b12, c1, d1 coefficients for air temperature
    ! a2, b21, b22, c2, d2 coefficients for water vapor mole fraction

    do ic = 1, ncan(p)

       ! Storage term

       rho_dz_over_dt = rhomol(p) * dz(p,ic) / dtime

       ! a1,b11,b12,c1,d1 coefficients for air temperature

       if (ic == 1) then
          gac_ic_minus_one = gac0(p)
       else
          gac_ic_minus_one = gac(p,ic-1)
       end if

       a1(ic) = -gac_ic_minus_one
       b11(ic) = rho_dz_over_dt + gac_ic_minus_one + gac(p,ic) &
               + gleaf_sh(ic,isun) * (1._r8 - alpha(ic,isun)) + gleaf_sh(ic,isha) * (1._r8 - alpha(ic,isha))
       b12(ic) = -gleaf_sh(ic,isun) * beta(ic,isun) - gleaf_sh(ic,isha) * beta(ic,isha)
       c1(ic) = -gac(p,ic)
       d1(ic) = rho_dz_over_dt * tair_bef(p,ic) + gleaf_sh(ic,isun) * delta(ic,isun) + gleaf_sh(ic,isha) * delta(ic,isha)

       ! Special case for top layer

       if (ic == ncan(p)) then
          c1(ic) = 0._r8
          d1(ic) = d1(ic) + gac(p,ic) * thref(p)
       end if

       ! Special case for first canopy layer (i.e., immediately above the ground)

       if (ic == 1) then
          a1(ic) = 0._r8
          b11(ic) = b11(ic) - gac0(p) * alpha0
          b12(ic) = b12(ic) - gac0(p) * beta0
          d1(ic) = d1(ic) + gac0(p) * delta0
       end if

       ! a2,b21,b22,c2,d2 coefficients for water vapor mole fraction

       if (ic == 1) then
          gac_ic_minus_one = gs0
       else
          gac_ic_minus_one = gac(p,ic-1)
       end if

       a2(ic) = -gac_ic_minus_one
       b21(ic) = -gleaf_et(ic,isun) * dqsat(ic,isun) * alpha(ic,isun) - gleaf_et(ic,isha) * dqsat(ic,isha) * alpha(ic,isha)
       b22(ic) = rho_dz_over_dt + gac_ic_minus_one + gac(p,ic) &
               + gleaf_et(ic,isun) * (1._r8 - dqsat(ic,isun) * beta(ic,isun)) &
               + gleaf_et(ic,isha) * (1._r8 - dqsat(ic,isha) * beta(ic,isha))
       c2(ic) = -gac(p,ic)
       d2(ic) = rho_dz_over_dt * (eair_bef(p,ic) / pref(p)) &
              + gleaf_et(ic,isun) * (dqsat(ic,isun) * delta(ic,isun) + qsat_term(ic,isun)) &
              + gleaf_et(ic,isha) * (dqsat(ic,isha) * delta(ic,isha) + qsat_term(ic,isha)) 

       ! Special case for top layer

       if (ic == ncan(p)) then
          c2(ic) = 0._r8
          d2(ic) = d2(ic) + gac(p,ic) * (eref(p) / pref(p))
       end if

       ! Special case for first canopy layer (i.e., immediately above the ground)

       if (ic == 1) then
          a2(ic) = 0._r8
          b21(ic) = b21(ic) - gs0 * rhg(p) * dqsat0 * alpha0
          b22(ic) = b22(ic) - gs0 * rhg(p) * dqsat0 * beta0
          d2(ic) = d2(ic) + gs0 * rhg(p) * (qsat0 + dqsat0 * (delta0 - tg_bef(p)))
       end if

    end do

    ! Solve for air temperature and water vapor (mol/mol):
    !
    ! a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
    ! a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
    !
    ! Note that as used here eair = mol/mol

    call tridiag_2eq (a1, b11, b12, c1, d1, a2, b21, b22, c2, d2, tair(p,:), eair(p,:), ncan(p))

    ! Soil surface temperature (K) and vapor pressure (mol/mol)

    t0 = alpha0 * tair(p,1) + beta0 * eair(p,1) + delta0
    e0 = rhg(p) * (qsat0 + dqsat0 * (t0 - tg_bef(p)))

    ! Leaf temperature

    do ic = 1, ncan(p)
       tleaf_implic(ic,isun) = alpha(ic,isun)*tair(p,ic) + beta(ic,isun)*eair(p,ic) + delta(ic,isun)
       tleaf_implic(ic,isha) = alpha(ic,isha)*tair(p,ic) + beta(ic,isha)*eair(p,ic) + delta(ic,isha)
    end do

    ! Convert water vapor from mol/mol to Pa

    do ic = 1, ncan(p)
       eair(p,ic) = eair(p,ic) * pref(p)
    end do
    e0 = e0 * pref(p)

    ! Use LeafFluxes to calculate leaf fluxes (per unit leaf area) for the
    ! current air temperature and vapor pressure profiles in the canopy. The
    ! flux and leaf temperature calculations there are the same as here, so
    ! the answers are the same in both routines. Compare answers to make sure
    ! the implicit flux-profile solution is correct.

    do ic = 1, ncan(p)
       call LeafFluxes (p, ic, isun, mlcanopy_inst)
       call LeafFluxes (p, ic, isha, mlcanopy_inst)
    end do

    ! Compare with leaf fluxes as calculated here. But remember that here the
    ! fluxes for sunlit/shaded leaves are multiplied by their leaf area

    do ic = 1, ncan(p)
       shsrc(ic) = 0._r8 ; etsrc(ic) = 0._r8 ; stveg(ic) = 0._r8
       if (dpai(p,ic) > 0._r8) then

          ! Layer fluxes as calculated here

          do il = 1, nleaf
             shsrc(ic) = shsrc(ic) + cpair(p) * (tleaf_implic(ic,il) - tair(p,ic)) * gleaf_sh(ic,il)
             call SatVap (tleaf_bef(p,ic,il), esat, desat)
             etsrc(ic) = etsrc(ic) + (esat + desat * (tleaf_implic(ic,il) - tleaf_bef(p,ic,il)) - eair(p,ic)) / pref(p) &
                       * gleaf_et(ic,il)
             stveg(ic) = stveg(ic) + heatcap(ic,il) * (tleaf_implic(ic,il) - tleaf_bef(p,ic,il)) / dtime
          end do

          ! Layer fluxes as calculated by LeafFluxes

          shsrc_leaf = (shleaf(p,ic,isun) * fracsun(p,ic) + shleaf(p,ic,isha) * (1._r8 - fracsun(p,ic))) * dpai(p,ic)
          etsrc_leaf = (trleaf(p,ic,isun) + evleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic) &
                     + (trleaf(p,ic,isha) + evleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
          stveg_leaf = (stleaf(p,ic,isun) * fracsun(p,ic) + stleaf(p,ic,isha) * (1._r8 - fracsun(p,ic))) * dpai(p,ic)

          ! Error checks

          if (abs(shsrc(ic)-shsrc_leaf) > 0.001_r8) then
             call endrun (msg=' ERROR: FluxProfileSolution: Leaf sensible heat flux error')
          end if

          if (abs(lambda*(etsrc(ic)-etsrc_leaf)) > 0.001_r8) then
             call endrun (msg=' ERROR: FluxProfileSolution: Leaf latent heat flux error')
          end if

          if (abs(stveg(ic)-stveg_leaf) > 0.001_r8) then
             call endrun (msg=' ERROR: FluxProfileSolution: Leaf heat storage error')
          end if

          if (abs(tleaf(p,ic,isun)-tleaf_implic(ic,isun)) > 1.e-06_r8) then
             call endrun (msg=' ERROR: FluxProfileSolution: Leaf temperature error (sunlit)')
          end if

          if (abs(tleaf(p,ic,isha)-tleaf_implic(ic,isha)) > 1.e-06_r8) then
             call endrun (msg=' ERROR: FluxProfileSolution: Leaf temperature error (shaded)')
          end if
       end if
    end do

    ! Use SoilFluxes to calculate soil fluxes. The flux and temperature
    ! calculations there are the same as here, so the answers are the
    ! same in both routines. Compare answers to make sure the
    ! implicit flux-profile solution is correct.

    call SoilFluxes (p, mlcanopy_inst)

    ! Compare with soil fluxes as calculated here

    sh0 = -cpair(p) * (tair(p,1) - t0) * gac0(p)
    et0 = -(eair(p,1) - e0) / pref(p) * gs0
    g0 = -soil_tk(p) / soil_dz(p) * soil_t(p) + soil_tk(p) / soil_dz(p) * t0

    if (abs(shsoi(p)-sh0) > 0.001_r8) then
       call endrun (msg=' ERROR: FluxProfileSolution: Soil sensible heat flux error')
    end if

    if (abs(lambda*(etsoi(p)-et0)) > 0.001_r8) then
       call endrun (msg=' ERROR: FluxProfileSolution: Soil latent heat flux error')
    end if

    if (abs(gsoi(p)-g0) > 0.001_r8) then
       call endrun (msg=' ERROR: FluxProfileSolution: Soil heat flux error')
    end if

    if (abs(tg(p)-t0) > 1.e-06_r8) then
       call endrun (msg=' ERROR: FluxProfileSolution: Soil surface temperature error')
    end if

    if (abs(eg(p)-e0) > 1.e-06_r8) then
       call endrun (msg=' ERROR: FluxProfileSolution: Soil surface vapor pressure error')
    end if

    ! Vertical sensible heat and water vapor fluxes between layers

    do ic = 1, ncan(p)-1
       shair(p,ic) = -cpair(p) * (tair(p,ic+1) - tair(p,ic)) * gac(p,ic)
       etair(p,ic) = -(eair(p,ic+1) - eair(p,ic)) / pref(p) * gac(p,ic)
    end do
    ic = ncan(p)
    shair(p,ic) = -cpair(p) * (thref(p) - tair(p,ic)) * gac(p,ic)
    etair(p,ic) = -(eref(p) - eair(p,ic)) / pref(p) * gac(p,ic)

    ! Canopy air storage flux (W/m2) and its component terms

    do ic = 1, ncan(p)
       storage_sh(ic) = rhomol(p) * cpair(p) * (tair(p,ic) - tair_bef(p,ic)) * dz(p,ic) / dtime
       storage_et(ic) = rhomol(p) * (eair(p,ic) - eair_bef(p,ic)) / pref(p) * dz(p,ic) / dtime
       stair(p,ic) = storage_sh(ic) + storage_et(ic) * lambda
    end do

    ! ---------------------------
    ! DO SOME CONSERVATION CHECKS
    ! ---------------------------

    ! Vegetation flux energy balance

    do ic = 1, ncan(p)
       err = avail_energy(ic,isun) + avail_energy(ic,isha) - shsrc(ic) - lambda * etsrc(ic) - stveg(ic)
       if (abs(err) > 0.001_r8) then
          call endrun (msg=' ERROR: FluxProfileSolution: Leaf energy balance error')
       end if
    end do

    ! Flux conservation at each layer

    do ic = 1, ncan(p)
       if (ic == 1) then
          err = storage_sh(ic) - (sh0 + shsrc(ic) - shair(p,ic))
       else
          err = storage_sh(ic) - (shair(p,ic-1) + shsrc(ic) - shair(p,ic))
       end if
       if (abs(err) > 0.001_r8) then
          call endrun (msg=' ERROR: FluxProfileSolution: Sensible heat layer conservation error')
       end if

       if (ic == 1) then
          err = storage_et(ic) - (et0 + etsrc(ic) - etair(p,ic))
       else
          err = storage_et(ic) - (etair(p,ic-1) + etsrc(ic) - etair(p,ic))
       end if
       err = err * lambda
       if (abs(err) > 0.001_r8) then
          call endrun (msg=' ERROR: FluxProfileSolution: Latent heat layer conservation error')
       end if
    end do

    ! Flux conservation for canopy sensible heat and latent heat. This is to
    ! check canopy conservation equation (so the sum is to ntop not ncan).

    sum_src = 0._r8 ; sum_storage = 0._r8
    do ic = 1, ntop(p)
       sum_src = sum_src + shsrc(ic)
       sum_storage = sum_storage + storage_sh(ic)
    end do

    err = (sh0 + sum_src - sum_storage) - shair(p,ntop(p))
    if (abs(err) > 0.001_r8) then
       call endrun (msg=' ERROR: FluxProfileSolution: Sensible heat canopy conservation error')
    end if

    sum_src = 0._r8 ; sum_storage = 0._r8
    do ic = 1, ntop(p)
       sum_src = sum_src + etsrc(ic)
       sum_storage = sum_storage + storage_et(ic)
    end do

    err = (et0 + sum_src - sum_storage) - etair(p,ntop(p))  ! mol H2O/m2/s
    err = err * lambda                                      ! W/m2
    if (abs(err) > 0.001_r8) then
       call endrun (msg=' ERROR: FluxProfileSolution: Latent heat canopy conservation error')
    end if

    ! Ground energy balance conservation

    err = rnsoi(p) - sh0 - lambda * et0 - g0
    if (abs(err) > 0.001_r8) then
       call endrun (msg=' ERROR: FluxProfileSolution: Ground temperature energy balance error')
    end if

    end associate
  end subroutine FluxProfileSolution

  !-----------------------------------------------------------------------
  subroutine LookupPsihatINI
    !
    ! !DESCRIPTION:
    ! Initialize the look-up tables needed to calculate the RSL psihat functions.
    ! Remember that in a netcdf file the dimensions appear in the opposite order:
    ! netcdf: psigridM_nc(nL,nZ) -> Fortran: psigridM(nZ,nL)
    ! netcdf: psigridH_nc(nL,nZ) -> Fortran: psigridH(nZ,nL)
    !
    ! !USES:
    use fileutils, only : getfil
    use ncdio_pio, only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t
    use ncdio_pio, only : ncd_inqdid, ncd_inqdlen
    use spmdMod, only : masterproc
    use clm_varctl, only : rslfile
    use MLclm_varcon, only : nZ, nL, dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
    !
    ! !ARGUMENTS:
    implicit none
    !
    !LOCAL VARIABLES
    character(len=256) :: locfn        ! Local file name
    type(file_desc_t) :: ncid          ! pio netCDF file id
    integer :: dimid                   ! netCDF dimension id
    logical :: readv                   ! read variable in or not

    real(r8) :: zdtgridM_nc(nZ)        ! netcdf data: Grid of zdt on which psihat is given for momentum
    real(r8) :: dtLgridM_nc(nL)        ! netcdf data: Grid of dtL on which psihat is given for momentum
    real(r8) :: psigridM_nc(nL,nZ)     ! netcdf data: Grid of psihat values for momentum
    real(r8) :: zdtgridH_nc(nZ)        ! netcdf data: Grid of zdt on which psihat is given for heat
    real(r8) :: dtLgridH_nc(nL)        ! netcdf data: Grid of dtL on which psihat is given for heat
    real(r8) :: psigridH_nc(nL,nZ)     ! netcdf data: Grid of psihat values for heat
    integer :: nZ_nc, nL_nc            ! netcdf data: dimensions
    integer :: ii, jj                  ! Looping indices
    !---------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read RSL look-up table .....'
    end if

    ! Get netcdf file

    call getfil (rslfile, locfn, 0)

    ! Open netcdf file

    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Check dimensions

    call ncd_inqdid (ncid, 'nZ', dimid)
    call ncd_inqdlen (ncid, dimid, nZ_nc)

    if (nZ_nc /= nZ) then
       call endrun (msg=' ERROR: LookupPsihatINI: nZ does not equal expected value')
    end if

    call ncd_inqdid (ncid, 'nL', dimid)
    call ncd_inqdlen (ncid, dimid, nL_nc)

    if (nL_nc /= nL) then
       call endrun (msg=' ERROR: LookupPsihatINI: nL does not equal expected value')
    end if

    ! Read variables

    call ncd_io('dtLgridM', dtLgridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading dtLgridM')

    call ncd_io('zdtgridM', zdtgridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading zdtgridM')

    call ncd_io('psigridM', psigridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading psigridM')

    call ncd_io('dtLgridH', dtLgridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading dtLgridH')

    call ncd_io('zdtgridH', zdtgridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading zdtgridH')

    call ncd_io('psigridH', psigridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading psigridH')

    ! Close netcdf file

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*) 'Successfully read RSL look-up table'
    end if

    ! Copy netcdf variables

    do jj = 1, nL
       dtLgridM(1,jj) = dtLgridM_nc(jj)
       dtLgridH(1,jj) = dtLgridH_nc(jj)
    end do

    do ii = 1, nZ
       zdtgridM(ii,1) = zdtgridM_nc(ii)
       zdtgridH(ii,1) = zdtgridH_nc(ii)
    end do

    do ii = 1, nZ
       do jj = 1, nL
          psigridM(ii,jj) = psigridM_nc(jj,ii)
          psigridH(ii,jj) = psigridH_nc(jj,ii)
       end do
    end do

    return

  end subroutine LookupPsihatINI

end module MLCanopyTurbulenceMod
