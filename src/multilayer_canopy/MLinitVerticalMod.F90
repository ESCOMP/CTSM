module MLinitVerticalMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialize multilayer canopy vertical structure and profiles
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use decompMod, only : bounds_type
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initVerticalStructure  ! Define canopy layer vertical structure
  public :: initVerticalProfiles   ! Initialize vertical profiles and canopy states
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initVerticalStructure (bounds, num_filter, filter, &
  canopystate_inst, frictionvel_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Define canopy layer vertical structure
    !
    ! !USES:
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : dpai_min, dht_tall, dht_short, dht_param
    use MLclm_varctl, only : pad_type
    use MLclm_varpar, only : nlevmlcan
    use MLMathToolsMod, only : beta_function
    use CanopyStateType, only : canopystate_type
    use FrictionVelocityMod, only : frictionvel_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filter            ! Number of patches in filter
    integer, intent(in) :: filter(:)             ! Patch filter
    type(canopystate_type), intent(in) :: canopystate_inst
    type(frictionvel_type), intent(in) :: frictionvel_inst
    type(mlcanopy_type), intent(out) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                               ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                               ! Aboveground layer index
    real(r8) :: dht                              ! Height increment (m)
    integer  :: n_zref                           ! Number of above-canopy layers
    real(r8) :: dz_zref                          ! Atmospheric reference height - canopy height (m)
    real(r8) :: beta                             ! Value of the beta function
    real(r8) :: zl, zu                           ! Bottom and top heights for a canopy layer (m)
    integer  :: num_int                          ! Number of layers for numerical integration of LAI between zl and zu
    real(r8) :: dz_int                           ! Height increment for numerical integration of LAI between zl and zu (m)
    integer  :: ic_int                           ! Do loop index for numerical integration of LAI between zl and zu
    real(r8) :: z_int                            ! Height for numerical integration of LAI between zl and zu (m)
    real(r8) :: zrel                             ! Height relative to canopy top (z_int/ztop)
    real(r8) :: beta_pdf                         ! Value for beta probability density function
    real(r8) :: lad                              ! Layer leaf area density (m2/m3)
    real(r8) :: lai_err, sai_err, pai_err        ! Sum of leaf, stem, or plant area index over all layers (m2/m2)
    real(r8) :: lai_miss, sai_miss               ! Missing leaf (stem) area after imposing dpai_min constraint on layers (m2/m2)

    real(r8) :: dlai(bounds%begp:bounds%endp,1:nlevmlcan)  ! Layer leaf area index (m2/m2)
    real(r8) :: dsai(bounds%begp:bounds%endp,1:nlevmlcan)  ! Layer stem area index (m2/m2)
    real(r8) ::   zw(bounds%begp:bounds%endp,0:nlevmlcan)  ! Height at interface between two adjacent layers (m)
    !---------------------------------------------------------------------

    associate ( &
                                                          ! *** Input ***
    pbeta       => pftcon%pbeta                      , &  ! CLM (new): Parameter for the leaf area density beta distribution
    qbeta       => pftcon%qbeta                      , &  ! CLM (new): Parameter for the leaf area density beta distribution
    forc_hgt_u  => frictionvel_inst%forc_hgt_u_patch , &  ! CLM: Atmospheric reference height (m)
    htop        => canopystate_inst%htop_patch       , &  ! CLM: Canopy height (m)
    elai        => canopystate_inst%elai_patch       , &  ! CLM: Leaf area index of canopy (m2/m2)
    esai        => canopystate_inst%esai_patch       , &  ! CLM: Stem area index of canopy (m2/m2)
                                                          ! *** Output ***
    zref        => mlcanopy_inst%zref_forcing        , &  ! Atmospheric reference height (m)
    ztop        => mlcanopy_inst%ztop_canopy         , &  ! Canopy height (m)
    ncan        => mlcanopy_inst%ncan_canopy         , &  ! Number of layers
    ntop        => mlcanopy_inst%ntop_canopy         , &  ! Index for top leaf layer
    nbot        => mlcanopy_inst%nbot_canopy         , &  ! Index for bottom leaf layer
    dlai_frac   => mlcanopy_inst%dlai_frac_profile   , &  ! Canopy layer leaf area index (fraction of canopy total)
    dsai_frac   => mlcanopy_inst%dsai_frac_profile   , &  ! Canopy layer stem area index (fraction of canopy total)
    zs          => mlcanopy_inst%zs_profile          , &  ! Canopy layer height for scalar concentration and source (m)
    dz          => mlcanopy_inst%dz_profile            &  ! Canopy layer thickness (m)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Top height of canopy

       ztop(p) = htop(p)

       ! Determine number of within-canopy layers by specifying height increment
       ! based on tall or short canopies as specified by dht_param

       if (ztop(p) > dht_param) then
          dht = dht_tall
       else
          dht = dht_short
       end if

       ! Define ntop as the number of within-canopy layers and adjust dht if needed

       ntop(p) = nint(ztop(p) / dht)
       dht = ztop(p) / float(ntop(p))

       ! Calculate heights at layer interfaces (zw). These are the heights
       ! for the conductance between two scalar concentrations. They are
       ! defined for ic = 0 (ground) to ic = ncan (top above-canopy layer).
       ! First calculate within-canopy heights (for ic = 0 to ntop)

       ic = ntop(p)
       zw(p,ic) = ztop(p)
       do ic = ntop(p)-1, 0, -1
          zw(p,ic) = zw(p,ic+1) - dht
       end do

       if (zw(p,0) < 0._r8 .and. zw(p,0) >= -1.e-9_r8) then
          zw(p,0) = 0._r8
       else if (abs(zw(p,0)) > 1.e-9_r8) then
          call endrun (msg=' ERROR: initVerticalStructure: zw(p,0) improperly defined')
       end if

       ! Now calculate the above-canopy layers and their heights
       ! (for ic = ntop+1 to ncan)

       zref(p) = forc_hgt_u(p)
       dz_zref = zref(p) - ztop(p)
       n_zref = nint(dz_zref / dht)
       dht = dz_zref / float(n_zref)
       ncan(p) = ntop(p) + n_zref

       if (ncan(p) > nlevmlcan) then
          call endrun (msg=' ERROR: initVerticalStructure: ncan > nlevmlcan')
       end if

       ic = ncan(p)
       zw(p,ic) = zref(p)
       do ic = ncan(p)-1, ntop(p)+1, -1
          zw(p,ic) = zw(p,ic+1) - dht
       end do

       ! Thickness of each layer

       do ic = 1, ncan(p)
          dz(p,ic) = zw(p,ic) - zw(p,ic-1)
       end do

       ! Determine heights of the scalar concentration and scalar source
       ! (zs). These are physically centered between the conductance points
       ! (i.e., in the middle of the layer).

       do ic = 1, ncan(p)
          zs(p,ic) = 0.5_r8 * (zw(p,ic) + zw(p,ic-1))
       end do

       ! Calculate leaf area index at each height by numerically integrating
       ! the leaf area density distribution between the bottom and top
       ! heights for that layer

       beta = beta_function (pbeta(patch%itype(p)), qbeta(patch%itype(p)))

       do ic = 1, ntop(p)
          zl = zw(p,ic-1) ; zu = zw(p,ic)
          dlai(p,ic) = 0._r8 ; dsai(p,ic) = 0._r8

          ! Use beta distribution with numerical integration between zl and zu

          if (pad_type == 1) then
             num_int = 100                         ! 100 sublayers for numerical integration
             dz_int = (zu - zl) / float(num_int)   ! dz for numerical integration
             do ic_int = 1, num_int

                if (ic_int == 1) then
                   z_int = zl + 0.5_r8 * dz_int
                else
                   z_int = z_int + dz_int
                end if
                zrel = min(z_int/ztop(p), 1._r8)
                beta_pdf = (1._r8 / beta) * zrel**(pbeta(patch%itype(p))-1._r8) * &
                           (1._r8 - zrel)**(qbeta(patch%itype(p))-1._r8)

                ! Leaf area

                lad = (elai(p) / ztop(p)) * beta_pdf
                dlai(p,ic) = dlai(p,ic) + lad * dz_int

                ! Stem area

                lad = (esai(p) / ztop(p)) * beta_pdf
                dsai(p,ic) = dsai(p,ic) + lad * dz_int

             end do
          end if

          ! Use a uniform profile

          if (pad_type == 2) then
             dlai(p,ic) = (elai(p) / ztop(p)) * (zu - zl)
             dsai(p,ic) = (esai(p) / ztop(p)) * (zu - zl)
          end if

       end do

       ! Check to make sure sum of numerical lai+sai matches canopy lai+sai

       pai_err = sum(dlai(p,1:ntop(p))) + sum(dsai(p,1:ntop(p)))
       if (abs(pai_err - (elai(p)+esai(p))) > 1.e-5_r8) then
          call endrun (msg=' ERROR: initVerticalStructure: plant area does not match CLM input')
       end if

       ! Now determine the number of layers with plant area. Set the layers with
       ! plant area < dpai_min to zero and sum the "missing" leaf and stem area
       ! so that this can be distributed back across the vegetation layers.

       lai_miss = 0._r8 ; sai_miss = 0._r8
       do ic = 1, ntop(p)
          if ((dlai(p,ic)+dsai(p,ic)) < dpai_min) then
             lai_miss = lai_miss + dlai(p,ic)
             sai_miss = sai_miss + dsai(p,ic)
             dlai(p,ic) = 0._r8
             dsai(p,ic) = 0._r8
          end if
       end do

       ! Distribute the missing leaf area across vegetation layers
       ! in proportion to the leaf area profile

       if (lai_miss > 0._r8) then
          lai_err = sum(dlai(p,1:ntop(p)))
          do ic = 1, ntop(p)
             dlai(p,ic) = dlai(p,ic) + lai_miss * (dlai(p,ic) / lai_err)
          end do
       end if

       ! Now do the same for stem area

       if (sai_miss > 0._r8) then
          sai_err = sum(dsai(p,1:ntop(p)))
          do ic = 1, ntop(p)
             dsai(p,ic) = dsai(p,ic) + sai_miss * (dsai(p,ic) / sai_err)
          end do
       end if

       ! Find the lowest leaf layer

       nbot(p) = 0
       do ic = ntop(p), 1, -1
          if ((dlai(p,ic)+dsai(p,ic)) > 0._r8) nbot(p) = ic
       end do

       if (nbot(p) == 0) then
          write (iulog,*) 'initVerticalStructure error: nbot not defined'
          call endrun()
       end if

       ! Error check

       lai_err = sum(dlai(p,1:ntop(p)))
       if (abs(lai_err - elai(p)) > 1.e-06_r8) then
          call endrun (msg=' ERROR: initVerticalStructure: leaf area does not match CLM input after redistribution')
       end if

       sai_err = sum(dsai(p,1:ntop(p)))
       if (abs(sai_err - esai(p)) > 1.e-5_r8) then
          call endrun (msg=' ERROR: initVerticalStructure: stem area does not match CLM input after redistribution')
       end if

       ! Zero out above-canopy layers

       do ic = ntop(p)+1, ncan(p)
          dlai(p,ic) = 0._r8
          dsai(p,ic) = 0._r8
       end do

       ! Normalize profiles

       do ic = 1, ncan(p)
          dlai_frac(p,ic) = dlai(p,ic) / elai(p)
          dsai_frac(p,ic) = dsai(p,ic) / esai(p)
       end do

    end do

    end associate
  end subroutine initVerticalStructure

  !-----------------------------------------------------------------------
  subroutine initVerticalProfiles (num_filter, filter, &
  atm2lnd_inst, mlcanopy_inst, wateratm2lndbulk_inst)
    !
    ! !DESCRIPTION:
    ! Initialize vertical profiles and canopy states
    !
    ! !USES:
    use PatchType, only : patch
    use MLclm_varcon, only : mmh2o, mmdry
    use MLclm_varpar, only : isun, isha
    use atm2lndType, only : atm2lnd_type
    use Wateratm2lndBulkType, only : wateratm2lndbulk_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter            ! Number of patches in filter
    integer, intent(in) :: filter(:)             ! Patch filter
    type(atm2lnd_type) , intent(in)    :: atm2lnd_inst
    type(wateratm2lndbulk_type), intent(in) :: wateratm2lndbulk_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                               ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                                ! Grid cell index for CLM g/l/c/p hierarchy
    integer  :: ic                               ! Aboveground layer index
    !---------------------------------------------------------------------

    associate ( &
                                                              ! *** Input ***
    forc_u      => atm2lnd_inst%forc_u_grc               , &  ! CLM: Atmospheric wind speed in east direction (m/s)
    forc_v      => atm2lnd_inst%forc_v_grc               , &  ! CLM: Atmospheric wind speed in north direction (m/s)
    forc_pco2   => atm2lnd_inst%forc_pco2_grc            , &  ! CLM: Atmospheric CO2 partial pressure (Pa)
    forc_t      => atm2lnd_inst%forc_t_downscaled_col    , &  ! CLM: Atmospheric temperature (K)
    forc_q      => wateratm2lndbulk_inst%forc_q_downscaled_col, &  ! CLM: Atmospheric specific humidity (kg/kg)
    forc_pbot   => atm2lnd_inst%forc_pbot_downscaled_col , &  ! CLM: Atmospheric pressure (Pa)
    ncan        => mlcanopy_inst%ncan_canopy             , &  ! Number of layers
                                                              ! *** Output ***
    taf         => mlcanopy_inst%taf_canopy              , &  ! Air temperature at canopy top (K)
    qaf         => mlcanopy_inst%qaf_canopy              , &  ! Specific humidity at canopy top (kg/kg)
    tg          => mlcanopy_inst%tg_soil                 , &  ! Soil surface temperature (K)
    gac0        => mlcanopy_inst%gac0_soil               , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    wind        => mlcanopy_inst%wind_profile            , &  ! Canopy layer wind speed (m/s)
    tair        => mlcanopy_inst%tair_profile            , &  ! Canopy layer air temperature (K)
    eair        => mlcanopy_inst%eair_profile            , &  ! Canopy layer vapor pressure (Pa)
    cair        => mlcanopy_inst%cair_profile            , &  ! Canopy layer atmospheric CO2 (umol/mol)
    h2ocan      => mlcanopy_inst%h2ocan_profile          , &  ! Canopy layer intercepted water (kg H2O/m2)
    lwpveg      => mlcanopy_inst%lwpveg_profile          , &  ! Canopy layer leaf water potential (MPa)
    tleaf       => mlcanopy_inst%tleaf_leaf                &  ! Leaf temperature (K)
    )

    ! Initialization of vertical profiles and canopy states

    do fp = 1, num_filter
       p = filter(fp)
       c = patch%column(p)
       g = patch%gridcell(p)

       do ic = 1, ncan(p)
          wind(p,ic) = max (1.0_r8, sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
          tair(p,ic) = forc_t(c)
          eair(p,ic) = forc_q(c) * forc_pbot(c) / (mmh2o/mmdry + (1._r8-mmh2o/mmdry) * forc_q(c))
          cair(p,ic) = forc_pco2(g) / forc_pbot(c) * 1.e06_r8

          tleaf(p,ic,isun) = forc_t(c)
          tleaf(p,ic,isha) = forc_t(c)
          h2ocan(p,ic) = 0._r8
          lwpveg(p,ic) = -0.1_r8
       end do

       taf(p) = forc_t(c)
       qaf(p) = forc_q(c)
       tg(p) = forc_t(c)
       gac0(p) = 1._r8 / 10._r8 * 42.3_r8
    end do

    end associate
  end subroutine initVerticalProfiles

end module MLinitVerticalMod
