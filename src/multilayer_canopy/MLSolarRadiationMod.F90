module MLSolarRadiationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar radiation transfer through canopy
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use decompMod, only : bounds_type
  use PatchType, only : patch
  use pftconMod, only : pftcon
  use shr_kind_mod, only : r8 => shr_kind_r8
  use MLCanopyFluxesType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SolarRadiation                 ! Main driver for radiative transfer
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: NormanRadiation               ! Norman radiative transfer
  private :: TwoStreamRadiation            ! Two-stream approximation radiative transfer
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SolarRadiation (bounds, num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Solar radiation transfer through canopy
    !
    ! !USES:
    use clm_varcon, only : pi => rpi
    use clm_varpar, only : numrad, ivis
    use MLclm_varctl, only : light_type
    use MLclm_varpar, only : nlevmlcan, isun, isha
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filter                    ! Number of patches in filter
    integer, intent(in) :: filter(:)                     ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                       ! Filter index
    integer  :: p                                        ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                        ! Column index for CLM g/l/c/p hierarchy
    integer  :: ic                                       ! Aboveground layer index
    integer  :: ib                                       ! Waveband index
    integer  :: j                                        ! Sky angle index
    real(r8) :: angle                                    ! Sky angle (5, 15, 25, 35, 45, 55, 65, 75, 85 degrees)
    real(r8) :: gdirj                                    ! Relative projected area of leaf elements in the direction of sky angle
    real(r8) :: chil                                     ! Departure of leaf angle from spherical orientation (-0.4 <= xl <= 0.6)
    real(r8) :: phi1, phi2                               ! Term in Ross-Goudriaan function for gdir
    real(r8) :: gdir                                     ! Relative projected area of leaf elements in the direction of solar beam
    real(r8) :: wl, ws                                   ! Leaf and stem fraction of canopy layer
    real(r8) :: rho(bounds%begp:bounds%endp,1:nlevmlcan,1:numrad)    ! Leaf/stem reflectance
    real(r8) :: tau(bounds%begp:bounds%endp,1:nlevmlcan,1:numrad)    ! Leaf/stem transmittance
    real(r8) :: omega(bounds%begp:bounds%endp,1:nlevmlcan,1:numrad)  ! Leaf/stem scattering coefficient
    real(r8) :: kb(bounds%begp:bounds%endp,1:nlevmlcan)              ! Direct beam extinction coefficient

    ! For Norman radiation
    real(r8) :: tbi(bounds%begp:bounds%endp,0:nlevmlcan) ! Exponential transmittance of direct beam onto canopy layer

    ! For two-stream radiation
    real(r8) :: asu                                       ! Single scattering albedo
    real(r8) :: tmp0,tmp1,tmp2                            ! Intermediate variables
    real(r8) :: avmu(bounds%begp:bounds%endp,1:nlevmlcan) ! Average inverse diffuse optical depth per unit leaf area
    real(r8) :: betad(bounds%begp:bounds%endp,1:nlevmlcan,1:numrad)   ! Upscatter parameter for diffuse radiation
    real(r8) :: betab(bounds%begp:bounds%endp,1:nlevmlcan,1:numrad)   ! Upscatter parameter for direct beam radiation
    !---------------------------------------------------------------------

    associate ( &
                                                          ! *** Input ***
    xl          => pftcon%xl                       , &    ! CLM: Departure of leaf angle from spherical orientation (-)
    rhol        => pftcon%rhol                     , &    ! CLM: Leaf reflectance (-)
    taul        => pftcon%taul                     , &    ! CLM: Leaf transmittance (-)
    rhos        => pftcon%rhos                     , &    ! CLM: Stem reflectance (-)
    taus        => pftcon%taus                     , &    ! CLM: Stem transmittance (-)
    clump_fac   => pftcon%clump_fac                , &    ! CLM (new): Foliage clumping index (-)
    solar_zen   => mlcanopy_inst%solar_zen_forcing , &    ! Solar zenith angle (radians)
!   lai         => mlcanopy_inst%lai_canopy        , &    ! Leaf area index of canopy (m2/m2)
!   sai         => mlcanopy_inst%sai_canopy        , &    ! Stem area index of canopy (m2/m2)
    ncan        => mlcanopy_inst%ncan_canopy       , &    ! Number of layers
    ntop        => mlcanopy_inst%ntop_canopy       , &    ! Index for top leaf layer
    nbot        => mlcanopy_inst%nbot_canopy       , &    ! Index for bottom leaf layer
    dlai        => mlcanopy_inst%dlai_profile      , &    ! Canopy layer leaf area index (m2/m2)
    dsai        => mlcanopy_inst%dsai_profile      , &    ! Canopy layer stem area index (m2/m2)
    dpai        => mlcanopy_inst%dpai_profile      , &    ! Canopy layer plant area index (m2/m2)
    sumpai      => mlcanopy_inst%sumpai_profile    , &    ! Canopy layer cumulative plant area index (m2/m2)
                                                          ! *** Output ***
    fracsun     => mlcanopy_inst%fracsun_profile   , &    ! Canopy layer sunlit fraction (-)
    tb          => mlcanopy_inst%tb_profile        , &    ! Canopy layer transmittance of direct beam radiation (-)
    td          => mlcanopy_inst%td_profile        , &    ! Canopy layer transmittance of diffuse radiation (-)
    swleaf      => mlcanopy_inst%swleaf_leaf       , &    ! Leaf absorbed solar radiation (W/m2 leaf)
    apar        => mlcanopy_inst%apar_leaf           &    ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    )

    ! Calculate canopy layer optical properties

    do fp = 1, num_filter
       p = filter(fp)

       ! Zero out variables for all layers

       rho(p,1:ncan(p),1:numrad) = 0._r8
       tau(p,1:ncan(p),1:numrad) = 0._r8
       omega(p,1:ncan(p),1:numrad) = 0._r8
       kb(p,1:ncan(p)) = 0._r8
       fracsun(p,1:ncan(p)) = 0._r8
       tb(p,1:ncan(p)) = 0._r8
       td(p,1:ncan(p)) = 0._r8
       tbi(p,0:ncan(p)) = 0._r8
       avmu(p,1:ncan(p)) = 0._r8
       betab(p,1:ncan(p),1:numrad) = 0._r8
       betad(p,1:ncan(p),1:numrad) = 0._r8

       do ic = ntop(p), 1 , -1

          if (dpai(p,ic) > 0._r8) then

             ! Weight reflectance and transmittance by lai and sai and calculate
             ! leaf scattering coefficient

             wl = dlai(p,ic) / dpai(p,ic) ; ws = dsai(p,ic) / dpai(p,ic)
             do ib = 1, numrad
                rho(p,ic,ib) = max(rhol(patch%itype(p),ib)*wl + rhos(patch%itype(p),ib)*ws, 1.e-06_r8)
                tau(p,ic,ib) = max(taul(patch%itype(p),ib)*wl + taus(patch%itype(p),ib)*ws, 1.e-06_r8)
                omega(p,ic,ib) = rho(p,ic,ib) + tau(p,ic,ib)
             end do

             ! Direct beam extinction coefficient

             chil = min(max(xl(patch%itype(p)), -0.4_r8), 0.6_r8)
             if (abs(chil) <= 0.01_r8) chil = 0.01_r8

             phi1 = 0.5_r8 - 0.633_r8 * chil - 0.330_r8 * chil * chil
             phi2 = 0.877_r8 * (1._r8 - 2._r8 * phi1)

             gdir = phi1 + phi2 * cos(solar_zen(p))
             kb(p,ic) = gdir / cos(solar_zen(p))
             kb(p,ic) = min(kb(p,ic), 40._r8)

             ! Sunlit fraction of layer

             fracsun(p,ic) = clump_fac(patch%itype(p)) * exp(-kb(p,ic) * sumpai(p,ic) * clump_fac(patch%itype(p)))

             ! Direct beam transmittance (tb) through a single layer

             tb(p,ic) = exp(-kb(p,ic) * dpai(p,ic) * clump_fac(patch%itype(p)))

             ! Diffuse transmittance through a single layer (also needed for longwave
             ! radiation). Estimated for nine sky angles in increments of 10 degrees.

             td(p,ic) = 0._r8
             do j = 1, 9
                angle = (5._r8 + (j - 1) * 10._r8) * pi / 180._r8
                gdirj = phi1 + phi2 * cos(angle)
                td(p,ic) = td(p,ic) &
                         + exp(-gdirj / cos(angle) * dpai(p,ic) * clump_fac(patch%itype(p))) * sin(angle) * cos(angle)
             end do
             td(p,ic) = td(p,ic) * 2._r8 * (10._r8 * pi / 180._r8)

             !-------------------------------------------------
             ! Special parameters for Norman radiative transfer
             !-------------------------------------------------

             ! Transmittance (tbi) of unscattered direct beam onto layer i

             if (ic == ntop(p)) then
                tbi(p,ntop(p)) = 1._r8
             else
                tbi(p,ic) = tbi(p,ic+1) * exp(-kb(p,ic+1) * dpai(p,ic+1) * clump_fac(patch%itype(p)))
             end if

             !-----------------------------------------------------
             ! Special parameters for two-stream radiative transfer
             !-----------------------------------------------------

             ! avmu - average inverse diffuse optical depth per unit leaf area

             avmu(p,ic) = (1._r8 - phi1/phi2 * log((phi1+phi2)/phi1)) / phi2

             ! betad - upscatter parameter for diffuse radiation
             ! betab - upscatter parameter for direct beam radiation

             do ib = 1, numrad

                ! upscatter parameter for diffuse radiation

                betad(p,ic,ib) = 0.5_r8 / omega(p,ic,ib) * ( rho(p,ic,ib) + tau(p,ic,ib) &
                               + (rho(p,ic,ib)-tau(p,ic,ib)) * ((1._r8+chil)/2._r8)**2 )

                ! upscatter parameter for direct beam radiation

                tmp0 = gdir + phi2 * cos(solar_zen(p))
                tmp1 = phi1 * cos(solar_zen(p))
                tmp2 = 1._r8 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1)
                asu = 0.5_r8 * omega(p,ic,ib) * gdir / tmp0 * tmp2
                betab(p,ic,ib) = (1._r8 + avmu(p,ic)*kb(p,ic)) / (omega(p,ic,ib)*avmu(p,ic)*kb(p,ic)) * asu
             end do

          end if

       end do

       ! Direct beam transmittance onto ground

       tbi(p,0) = tbi(p,nbot(p)) * exp(-kb(p,nbot(p)) * dpai(p,nbot(p)) * clump_fac(patch%itype(p)))

    end do

    ! Calculate radiative transfer through canopy

    select case (light_type)
    case (1)
       call NormanRadiation (bounds, num_filter, filter, &
       rho, tau, omega, tbi, mlcanopy_inst)
    case (2)
       call TwoStreamRadiation (bounds, num_filter, filter, &
       omega, kb, avmu, betad, betab, mlcanopy_inst)
    case default
       call endrun (msg=' ERROR: SolarRadiation: light_type not valid')
    end select

    ! APAR per unit sunlit and shaded leaf area

    do fp = 1, num_filter
       p = filter(fp)
       do ic = 1, ncan(p)
          apar(p,ic,isun) = swleaf(p,ic,isun,ivis) * 4.6_r8
          apar(p,ic,isha) = swleaf(p,ic,isha,ivis) * 4.6_r8
       end do
    end do

    end associate
  end subroutine SolarRadiation

  !-----------------------------------------------------------------------
  subroutine NormanRadiation (bounds, num_filter, filter, &
  rho, tau, omega, tbi, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute solar radiation transfer through canopy using Norman (1979)
    !
    ! !USES:
    use clm_varpar, only : numrad
    use MLclm_varpar, only : nlevmlcan, isun, isha
    use MLMathToolsMod, only : tridiag
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer , intent(in) :: num_filter                                           ! Number of patches in filter
    integer , intent(in) :: filter(:)                                            ! Patch filter
    real(r8), intent(in) :: rho(bounds%begp:bounds%endp,1:nlevmlcan,1:numrad)    ! Leaf/stem reflectance
    real(r8), intent(in) :: tau(bounds%begp:bounds%endp,1:nlevmlcan,1:numrad)    ! Leaf/stem transmittance
    real(r8), intent(in) :: omega(bounds%begp:bounds%endp,1:nlevmlcan,1:numrad)  ! Leaf/stem scattering coefficient
    real(r8), intent(in) :: tbi(bounds%begp:bounds%endp,0:nlevmlcan) ! Exponential transmittance of direct beam onto canopy layer
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                               ! Filter index
    integer  :: p                                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                                               ! Aboveground layer index
    integer  :: icm1                                             ! Layer below ic (ic-1)
    integer  :: ib                                               ! Waveband index
    real(r8) :: suminc                                           ! Incident radiation for energy conservation check
    real(r8) :: sumref                                           ! Reflected radiation for energy conservation check
    real(r8) :: sumabs                                           ! Absorbed radiation for energy conservation check
    real(r8) :: err                                              ! Error check (W/m2)
    real(r8) :: refld                                            ! Term for diffuse radiation reflected by layer
    real(r8) :: trand                                            ! Term for diffuse radiation transmitted by layer
    integer  :: m                                                ! Index to the tridiagonal matrix
    real(r8) :: aic, bic                                         ! Intermediate terms for tridiagonal matrix
    real(r8) :: eic, fic                                         ! Intermediate terms for tridiagonal matrix
    integer, parameter :: neq = (nlevmlcan+1)*2                  ! Number of tridiagonal equations to solve
    real(r8) :: atri(neq), btri(neq)                             ! Entries in tridiagonal matrix
    real(r8) :: ctri(neq), dtri(neq)                             ! Entries in tridiagonal matrix
    real(r8) :: utri(neq)                                        ! Tridiagonal solution
    real(r8) :: swbeam                                           ! Direct beam solar flux onto canopy layer (W/m2 ground)
    real(r8) :: swabsb                                           ! Absorbed direct beam solar flux for canopy layer (W/m2 ground)
    real(r8) :: swabsd                                           ! Absorbed diffuse solar flux for canopy layer (W/m2 ground)
    real(r8) :: swsun                                            ! Absorbed solar radiation, sunlit fraction of layer (W/m2 ground)
    real(r8) :: swsha                                            ! Absorbed solar radiation, shaded fraction of layer (W/m2 ground)
    real(r8) :: swup(bounds%begp:bounds%endp,0:nlevmlcan,numrad) ! Upward diffuse solar flux above canopy layer (W/m2 ground)
    real(r8) :: swdn(bounds%begp:bounds%endp,0:nlevmlcan,numrad) ! Downward diffuse solar flux onto canopy layer (W/m2 ground)
    !---------------------------------------------------------------------

    associate ( &
                                                      ! *** Input ***
    swskyb      => mlcanopy_inst%swskyb_forcing  , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd_forcing  , &  ! Atmospheric diffuse solar radiation (W/m2)
    ncan        => mlcanopy_inst%ncan_canopy     , &  ! Number of layers
    ntop        => mlcanopy_inst%ntop_canopy     , &  ! Index for top leaf layer
    nbot        => mlcanopy_inst%nbot_canopy     , &  ! Index for bottom leaf layer
    albsoib     => mlcanopy_inst%albsoib_soil    , &  ! Direct beam albedo of ground (-)
    albsoid     => mlcanopy_inst%albsoid_soil    , &  ! Diffuse albedo of ground (-)
    dpai        => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    fracsun     => mlcanopy_inst%fracsun_profile , &  ! Canopy layer sunlit fraction (-)
    tb          => mlcanopy_inst%tb_profile      , &  ! Canopy layer transmittance of direct beam radiation (-)
    td          => mlcanopy_inst%td_profile      , &  ! Canopy layer transmittance of diffuse radiation (-)
                                                      ! *** Output ***
    swveg       => mlcanopy_inst%swveg_canopy    , &  ! Absorbed solar radiation: vegetation (W/m2)
    swvegsun    => mlcanopy_inst%swvegsun_canopy , &  ! Absorbed solar radiation: sunlit canopy (W/m2)
    swvegsha    => mlcanopy_inst%swvegsha_canopy , &  ! Absorbed solar radiation: shaded canopy (W/m2)
    albcan      => mlcanopy_inst%albcan_canopy   , &  ! Albedo above canopy (-)
    swsoi       => mlcanopy_inst%swsoi_soil      , &  ! Absorbed solar radiation: ground (W/m2)
    swleaf      => mlcanopy_inst%swleaf_leaf       &  ! Leaf absorbed solar radiation (W/m2 leaf)
    )

    do ib = 1, numrad
       do fp = 1, num_filter
          p = filter(fp)

          ! Zero out radiative fluxes for all layers

          swup(p,0,ib) = 0._r8
          swdn(p,0,ib) = 0._r8

          do ic = 1, ncan(p)
             swup(p,ic,ib) = 0._r8
             swdn(p,ic,ib) = 0._r8
             swleaf(p,ic,isun,ib) = 0._r8
             swleaf(p,ic,isha,ib) = 0._r8
          end do

          ! Set up and solve tridiagonal system of equations for upward and
          ! downward fluxes. There are two equations for each layer and the
          ! soil. These equations are referenced by "m". The first equation
          ! is the upward flux and the second equation is the downward flux.

          m = 0

          ! Soil: upward flux

          m = m + 1 
          atri(m) = 0._r8 
          btri(m) = 1._r8 
          ctri(m) = -albsoid(p,ib) 
          dtri(m) = swskyb(p,ib) * tbi(p,0) * albsoib(p,ib) 

          ! Soil: downward flux

          refld = (1._r8 - td(p,nbot(p))) * rho(p,nbot(p),ib) 
          trand = (1._r8 - td(p,nbot(p))) * tau(p,nbot(p),ib) + td(p,nbot(p)) 
          aic = refld - trand * trand / refld 
          bic = trand / refld 

          m = m + 1 
          atri(m) = -aic 
          btri(m) = 1._r8 
          ctri(m) = -bic 
          dtri(m) = swskyb(p,ib) * tbi(p,nbot(p)) * (1._r8 - tb(p,nbot(p))) * (tau(p,nbot(p),ib) - rho(p,nbot(p),ib) * bic)

          ! Leaf layers, excluding top layer

          do ic = nbot(p), ntop(p)-1

             ! Upward flux

             refld = (1._r8 - td(p,ic)) * rho(p,ic,ib)
             trand = (1._r8 - td(p,ic)) * tau(p,ic,ib) + td(p,ic)
             fic = refld - trand * trand / refld
             eic = trand / refld

             m = m + 1
             atri(m) = -eic
             btri(m) = 1._r8
             ctri(m) = -fic
             dtri(m) = swskyb(p,ib) * tbi(p,ic) * (1._r8 - tb(p,ic)) * (rho(p,ic,ib) - tau(p,ic,ib) * eic)

             ! Downward flux

             refld = (1._r8 - td(p,ic+1)) * rho(p,ic+1,ib)
             trand = (1._r8 - td(p,ic+1)) * tau(p,ic+1,ib) + td(p,ic+1)
             aic = refld - trand * trand / refld
             bic = trand / refld

             m = m + 1
             atri(m) = -aic
             btri(m) = 1._r8
             ctri(m) = -bic
             dtri(m) = swskyb(p,ib) * tbi(p,ic+1) * (1._r8 - tb(p,ic+1)) * (tau(p,ic+1,ib) - rho(p,ic+1,ib) * bic)

          end do

          ! Top canopy layer: upward flux

          ic = ntop(p)
          refld = (1._r8 - td(p,ic)) * rho(p,ic,ib)
          trand = (1._r8 - td(p,ic)) * tau(p,ic,ib) + td(p,ic)
          fic = refld - trand * trand / refld
          eic = trand / refld

          m = m + 1
          atri(m) = -eic
          btri(m) = 1._r8
          ctri(m) = -fic
          dtri(m) = swskyb(p,ib) * tbi(p,ic) * (1._r8 - tb(p,ic)) * (rho(p,ic,ib) - tau(p,ic,ib) * eic)

          ! Top canopy layer: downward flux

          m = m + 1
          atri(m) = 0._r8
          btri(m) = 1._r8
          ctri(m) = 0._r8
          dtri(m) = swskyd(p,ib)

          ! Solve tridiagonal system of equations for upward and downward fluxes

          call tridiag (atri, btri, ctri, dtri, utri, m)

          ! Now copy the solution (utri) to the upward (swup) and downward (swdn)
          ! fluxes for each layer
          ! swup =  Upward diffuse flux above layer
          ! swdn =  Downward diffuse flux onto layer

          m = 0

          ! Soil fluxes

          m = m + 1
          swup(p,0,ib) = utri(m)
          m = m + 1
          swdn(p,0,ib) = utri(m)

          ! Leaf layer fluxes

          do ic = nbot(p), ntop(p)
             m = m + 1
             swup(p,ic,ib) = utri(m)
             m = m + 1
             swdn(p,ic,ib) = utri(m)
          end do

       end do             ! end patch loop
    end do                ! end waveband loop

    ! Compute fluxes

    do ib = 1, numrad
       do fp = 1, num_filter
          p = filter(fp)

          ! Solar radiation absorbed by ground (soil)

          swbeam = tbi(p,0) * swskyb(p,ib)
          swabsb = swbeam * (1._r8 - albsoib(p,ib))
          swabsd = swdn(p,0,ib) * (1._r8 - albsoid(p,ib))
          swsoi(p,ib) = swabsb + swabsd

          ! Leaf layer fluxes

          swveg(p,ib) = 0._r8
          swvegsun(p,ib) = 0._r8
          swvegsha(p,ib) = 0._r8

          do ic = nbot(p), ntop(p)

             ! Downward direct beam incident on layer and absorbed direct
             ! beam and diffuse for layer. Note the special case for first
             ! leaf layer, where the upward flux from below is from the ground.
             ! The ground is ic=0, but nbot-1 will not equal 0 if there are lower
             ! canopy layers without leaves. This coding accommodates open
             ! layers beneath the bottom of the canopy.

             swbeam = tbi(p,ic) * swskyb(p,ib)
             swabsb = swbeam * (1._r8 - tb(p,ic)) * (1._r8 - omega(p,ic,ib))
             if (ic == nbot(p)) then
                icm1 = 0
             else
                icm1 = ic - 1
             end if
             swabsd = (swdn(p,ic,ib) + swup(p,icm1,ib)) * (1._r8 - td(p,ic)) * (1._r8 - omega(p,ic,ib))

             ! Absorbed radiation for shaded and sunlit portions of layer

             swsha = swabsd * (1._r8 - fracsun(p,ic))
             swsun = swabsd * fracsun(p,ic) + swabsb 

             ! Per unit sunlit and shaded leaf area

             swleaf(p,ic,isun,ib) = swsun / (fracsun(p,ic) * dpai(p,ic))
             swleaf(p,ic,isha,ib) = swsha / ((1._r8 - fracsun(p,ic)) * dpai(p,ic))

             ! Sum solar radiation absorbed by vegetation and sunlit/shaded leaves

             swveg(p,ib) = swveg(p,ib) + (swabsb + swabsd)
             swvegsun(p,ib) = swvegsun(p,ib) + swsun
             swvegsha(p,ib) = swvegsha(p,ib) + swsha

          end do  

          ! Albedo

          suminc = swskyb(p,ib) + swskyd(p,ib)
          if (suminc > 0._r8) then
             albcan(p,ib) = swup(p,ntop(p),ib) / suminc
          else
             albcan(p,ib) = 0._r8
          end if

          ! Conservation check for total radiation balance: absorbed = incoming - outgoing

          sumref = albcan(p,ib) * (swskyb(p,ib) + swskyd(p,ib))
          sumabs = suminc - sumref
          err = sumabs - (swveg(p,ib) + swsoi(p,ib))
          if (abs(err) > 1.e-03_r8) then
             call endrun (msg='ERROR: NormanRadiation: total solar conservation error')
          end if

          ! Sunlit and shaded absorption

          err = (swvegsun(p,ib) + swvegsha(p,ib)) - swveg(p,ib) 
          if (abs(err) > 1.e-03_r8) then
             call endrun (msg='ERROR: NormanRadiation: sunlit/shade solar conservation error')
          end if

       end do             ! end patch loop
    end do                ! end waveband loop

    end associate
  end subroutine NormanRadiation

  !-----------------------------------------------------------------------
  subroutine TwoStreamRadiation (bounds, num_filter, filter, &
  omega, kb, avmu, betad, betab, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute solar radiation transfer through canopy using the two-stream
    ! approximation. This solution uses the leaf-scale fluxes (per unit
    ! leaf area) and so requires thin canopy layers.
    !
    ! !USES:
    use clm_varpar, only : numrad
    use MLclm_varpar, only : nlevmlcan, isun, isha
    !
    ! Arguments
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filter                             ! Number of patches in filter
    integer, intent(in) :: filter(:)                              ! Patch filter
    real(r8), intent(in) :: omega(bounds%begp:bounds%endp,1:nlevmlcan,numrad) ! Leaf/stem scattering coefficient
    real(r8), intent(in) :: kb(bounds%begp:bounds%endp,1:nlevmlcan)           ! Direct beam extinction coefficient
    real(r8), intent(in) :: avmu(bounds%begp:bounds%endp,1:nlevmlcan)         ! Average inverse diffuse optical depth per unit leaf area
    real(r8), intent(in) :: betad(bounds%begp:bounds%endp,1:nlevmlcan,numrad) ! Upscatter parameter for diffuse radiation
    real(r8), intent(in) :: betab(bounds%begp:bounds%endp,1:nlevmlcan,numrad) ! Upscatter parameter for direct beam radiation
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! Local variable declarations
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ib                      ! Waveband index
    integer  :: ic                      ! Aboveground layer index
    real(r8) :: b,c,d,h,u,v             ! Intermediate two-stream parameters
    real(r8) :: g1,g2                   ! Intermediate two-stream parameters
    real(r8) :: s1,s2                   ! Intermediate two-stream parameters
    real(r8) :: num1,num2               ! Intermediate two-stream parameters
    real(r8) :: den1,den2               ! Intermediate two-stream parameters
    real(r8) :: n1b,n2b                 ! Two-stream parameters
    real(r8) :: n1d,n2d                 ! Two-stream parameters
    real(r8) :: a1b,a2b                 ! Parameter for sunlit/shaded leaf radiation absorption
    real(r8) :: a1d,a2d                 ! Parameter for sunlit/shaded leaf radiation absorption

    real(r8) :: iupwb0                  ! Direct beam flux scattered upward (reflected) above canopy (W/m2)
    real(r8) :: iupwb                   ! Direct beam flux scattered upward at the canopy depth (W/m2)
    real(r8) :: idwnb                   ! Direct beam flux scattered downward below canopy (W/m2)
    real(r8) :: iabsb                   ! Direct beam flux absorbed by canopy (W/m2)
    real(r8) :: iabsb_sun               ! Direct beam flux absorbed by sunlit canopy (W/m2)
    real(r8) :: iabsb_sha               ! Direct beam flux absorbed by shaded canopy (W/m2)

    real(r8) :: iupwd0                  ! Diffuse flux scattered upward (reflected) above canopy (W/m2)
    real(r8) :: iupwd                   ! Diffuse flux scattered upward at the canopy depth (W/m2)
    real(r8) :: idwnd                   ! Diffuse flux scattered downward below canopy (W/m2)
    real(r8) :: iabsd                   ! Diffuse flux absorbed by canopy (W/m2)
    real(r8) :: iabsd_sun               ! Diffuse flux absorbed by sunlit canopy (W/m2)
    real(r8) :: iabsd_sha               ! Diffuse flux absorbed by shaded canopy (W/m2)

    real(r8) :: diupwb, didwnb          ! Flux derivatives
    real(r8) :: diupwd, didwnd          ! Flux derivatives
    real(r8) :: ilbs, ild               ! Leaf fluxes (per unit leaf area)

    real(r8) :: dir,dif                 ! Direct beam and diffuse fluxes
    real(r8) :: sun,sha                 ! Sunlit and shaded fluxes
    real(r8) :: suminc                  ! Incident radiation for energy conservation check
    real(r8) :: sumref                  ! Reflected radiation for energy conservation check
    real(r8) :: sumabs                  ! Absorbed radiation for energy conservation check

    integer :: ik
    !---------------------------------------------------------------------

    associate ( &
                                                      ! *** Input ***
    clump_fac   => pftcon%clump_fac              , &  ! CLM (new): Foliage clumping index (-)
    swskyb      => mlcanopy_inst%swskyb_forcing  , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd_forcing  , &  ! Atmospheric diffuse solar radiation (W/m2)
    ncan        => mlcanopy_inst%ncan_canopy     , &  ! Number of layers
    ntop        => mlcanopy_inst%ntop_canopy     , &  ! Index for top leaf layer
    nbot        => mlcanopy_inst%nbot_canopy     , &  ! Index for bottom leaf layer
    lai         => mlcanopy_inst%lai_canopy      , &  ! Leaf area index of canopy (m2/m2)
    sai         => mlcanopy_inst%sai_canopy      , &  ! Stem area index of canopy (m2/m2)
    albsoib     => mlcanopy_inst%albsoib_soil    , &  ! Direct beam albedo of ground (-)
    albsoid     => mlcanopy_inst%albsoid_soil    , &  ! Diffuse albedo of ground (-)
    dpai        => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    sumpai      => mlcanopy_inst%sumpai_profile  , &  ! Canopy layer cumulative plant area index (m2/m2)
    fracsun     => mlcanopy_inst%fracsun_profile , &  ! Canopy layer sunlit fraction (-)
                                                      ! *** Output ***
    swveg       => mlcanopy_inst%swveg_canopy    , &  ! Absorbed solar radiation: vegetation (W/m2)
    swvegsun    => mlcanopy_inst%swvegsun_canopy , &  ! Absorbed solar radiation: sunlit canopy (W/m2)
    swvegsha    => mlcanopy_inst%swvegsha_canopy , &  ! Absorbed solar radiation: shaded canopy (W/m2)
    albcan      => mlcanopy_inst%albcan_canopy   , &  ! Albedo above canopy (-)
    swsoi       => mlcanopy_inst%swsoi_soil      , &  ! Absorbed solar radiation: ground (W/m2)
    swleaf      => mlcanopy_inst%swleaf_leaf       &  ! Leaf absorbed solar radiation (W/m2 leaf)
    )

    do ib = 1, numrad
       do fp = 1, num_filter
          p = filter(fp)

    ik = ntop(p)

          !----------------------------------------------------------------
          ! Canopy fluxes using cumulative lai+sai
          !----------------------------------------------------------------

          ! Common terms

          b = (1._r8 - (1._r8 - betad(p,ik,ib)) * omega(p,ik,ib)) / avmu(p,ik)
          c = betad(p,ik,ib) * omega(p,ik,ib) / avmu(p,ik)
          h = sqrt(b*b - c*c)
          u = (h - b - c) / (2._r8 * h)
          v = (h + b + c) / (2._r8 * h)
          d = omega(p,ik,ib) * kb(p,ik) * swskyb(p,ib) / (h*h - kb(p,ik)*kb(p,ik))
          g1 = (betab(p,ik,ib) * kb(p,ik) - b * betab(p,ik,ib) - c * (1._r8 - betab(p,ik,ib))) * d
          g2 = ((1._r8 - betab(p,ik,ib)) * kb(p,ik) + c * betab(p,ik,ib) + b * (1._r8 - betab(p,ik,ib))) * d
          s1 = exp(-h * (lai(p)+sai(p)) * clump_fac(patch%itype(p)))
          s2 = exp(-kb(p,ik) * (lai(p)+sai(p)) * clump_fac(patch%itype(p)))

          ! Terms for direct beam radiation

          num1 = v * (g1 + g2 * albsoid(p,ib) + albsoib(p,ib) * swskyb(p,ib)) * s2
          num2 = g2 * (u + v * albsoid(p,ib)) * s1
          den1 = v * (v + u * albsoid(p,ib)) / s1
          den2 = u * (u + v * albsoid(p,ib)) * s1
          n2b = (num1 - num2) / (den1 - den2)
          n1b = (g2 - n2b * u) / v

          a1b = -g1 *      (1._r8 - s2*s2) / (2._r8 * kb(p,ik)) &
              +  n1b * u * (1._r8 - s2*s1) / (kb(p,ik) + h) + n2b * v * (1._r8 - s2/s1) / (kb(p,ik) - h)
          a2b =  g2 *      (1._r8 - s2*s2) / (2._r8 * kb(p,ik)) &
              -  n1b * v * (1._r8 - s2*s1) / (kb(p,ik) + h) - n2b * u * (1._r8 - s2/s1) / (kb(p,ik) - h)

          ! Direct beam radiative fluxes

          iupwb0 = -g1 + n1b * u + n2b * v
          iupwb = -g1 * s2 + n1b * u * s1 + n2b * v / s1
          idwnb =  g2 * s2 - n1b * v * s1 - n2b * u / s1
          iabsb = swskyb(p,ib) - iupwb0 - (1._r8 - albsoid(p,ib)) * idwnb - (1._r8 - albsoib(p,ib)) * swskyb(p,ib) * s2
          iabsb_sun = (1._r8 - omega(p,ik,ib)) * ((1._r8 - s2) * swskyb(p,ib) + 1._r8 / avmu(p,ik) * (a1b + a2b) * clump_fac(patch%itype(p)))
          iabsb_sha = iabsb - iabsb_sun

          ! Terms for diffuse radiation
 
          num1 = swskyd(p,ib) * (u + v * albsoid(p,ib)) * s1
          den1 = v * (v + u * albsoid(p,ib)) / s1
          den2 = u * (u + v * albsoid(p,ib)) * s1
          n2d = num1 / (den1 - den2)
          n1d = -(swskyd(p,ib) + n2d * u) / v

          a1d =  n1d * u * (1._r8 - s2*s1) / (kb(p,ik) + h) + n2d * v * (1._r8 - s2/s1) / (kb(p,ik) - h)
          a2d = -n1d * v * (1._r8 - s2*s1) / (kb(p,ik) + h) - n2d * u * (1._r8 - s2/s1) / (kb(p,ik) - h)

          ! Diffuse radiative fluxes

          iupwd0 = n1d * u + n2d * v
          iupwd =  n1d * u * s1 + n2d * v / s1
          idwnd = -n1d * v * s1 - n2d * u / s1
          iabsd = swskyd(p,ib) - iupwd0 - (1._r8 - albsoid(p,ib)) * idwnd
          iabsd_sun = (1._r8 - omega(p,ik,ib)) / avmu(p,ik) * (a1d + a2d) * clump_fac(patch%itype(p))
          iabsd_sha = iabsd - iabsd_sun

          !----------------------------------------------------------------
          ! Save necessary radiative fluxes
          !----------------------------------------------------------------

          ! Albedo

          suminc = swskyb(p,ib) + swskyd(p,ib)
          sumref = iupwb0 + iupwd0
          if (suminc > 0._r8) then
             albcan(p,ib) = sumref / suminc
          else
             albcan(p,ib) = 0._r8
          end if

          ! Solar radiation absorbed by canopy

          swveg(p,ib) = iabsb +  iabsd
          swvegsun(p,ib) = iabsb_sun + iabsd_sun
          swvegsha(p,ib) = iabsb_sha + iabsd_sha

          ! Solar radiation absorbed by ground (soil)

          dir = swskyb(p,ib) * s2 * (1._r8 - albsoib(p,ib))
          dif = (idwnb + idwnd) * (1._r8 - albsoid(p,ib))
          swsoi(p,ib) = dir + dif

          ! Conservation check: total incident = total reflected + total absorbed

          suminc = swskyb(p,ib) + swskyd(p,ib)
          sumref = iupwb0 + iupwd0
          sumabs = swveg(p,ib) + swsoi(p,ib)

          if (abs(suminc - (sumabs+sumref)) >= 1.e-06_r8) then
             call endrun (msg='ERROR: TwoStreamRadiation: total solar radiation conservation error')
          end if

          !----------------------------------------------------------------
          ! Repeat two-stream calculations for each leaf layer to
          ! calculate leaf fluxes
          !----------------------------------------------------------------

          ! Zero out fluxes for all layers

          do ic = 1, ncan(p)
             swleaf(p,ic,isun,ib) = 0._r8
             swleaf(p,ic,isha,ib) = 0._r8
          end do

          ! Calculate fluxes for leaf layers

          do ic = nbot(p), ntop(p)

             ! s1 and s2 depend on cumulative plant area index

             s1 = exp(-h * sumpai(p,ic) * clump_fac(patch%itype(p)))
             s2 = exp(-kb(p,ik) * sumpai(p,ic) * clump_fac(patch%itype(p)))

             ! ilbs - absorbed direct beam flux (scattered direct component) per unit leaf area
             ! at cumulative LAI, average for all leaves (J / m2 leaf / s)

             diupwb =  kb(p,ik) * g1 * s2 - h * n1b * u * s1 + h * n2b * v / s1
             didwnb = -kb(p,ik) * g2 * s2 + h * n1b * v * s1 - h * n2b * u / s1
             ilbs = (omega(p,ik,ib) * kb(p,ik) * swskyb(p,ib) * s2 + (diupwb - didwnb)) * clump_fac(patch%itype(p))

             ! ild - absorbed diffuse flux per unit leaf area at cumulative LAI,
             ! average for all leaves (J / m2 leaf / s)

             diupwd = -h * n1d * u * s1 + h * n2d * v / s1
             didwnd =  h * n1d * v * s1 - h * n2d * u / s1
             ild = (diupwd - didwnd) * clump_fac(patch%itype(p))

             ! Save leaf fluxes per unit sunlit and shaded leaf area (W/m2 leaf)

             swleaf(p,ic,isun,ib) = (1._r8 - omega(p,ik,ib)) * kb(p,ik) * swskyb(p,ib) + (ilbs + ild)
             swleaf(p,ic,isha,ib) =  ilbs + ild

          end do

       end do      ! end grid point loop
    end do         ! end waveband loop

    !---------------------------------------------------------------------
    ! Adjust leaf fluxes as needed. The sum of the fluxes for sunlit and
    ! shaded leaves should equal the total absorbed by the canopy, but may
    ! not because of inaccuracies in the flux derivatives (this is a small
    ! error if the dpai increment is small). Normalize these fluxes to sum
    ! to the canopy absorption.
    !---------------------------------------------------------------------

    do ib = 1, numrad
       do fp = 1, num_filter
          p = filter(fp)

          ! Sum canopy absorption (W/m2 ground) using leaf fluxes per unit sunlit
          ! and shaded leaf area (W/m2 leaf)

          sumabs = 0._r8
          do ic = nbot(p), ntop(p)
             sun = swleaf(p,ic,isun,ib) * fracsun(p,ic) * dpai(p,ic)
             sha = swleaf(p,ic,isha,ib) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             sumabs = sumabs + sun + sha
          end do

          ! Normalize profile

          if (sumabs > 0.0_r8) then
             do ic = nbot(p), ntop(p)
                swleaf(p,ic,isun,ib) = swleaf(p,ic,isun,ib) * swveg(p,ib) / sumabs
                swleaf(p,ic,isha,ib) = swleaf(p,ic,isha,ib) * swveg(p,ib) / sumabs
             end do
          end if

       end do
    end do

    end associate
  end subroutine TwoStreamRadiation

end module MLSolarRadiationMod
