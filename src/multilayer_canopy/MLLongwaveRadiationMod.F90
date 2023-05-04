module MLLongwaveRadiationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate longwave radiation transfer through canopy
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
  public :: LongwaveRadiation
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LongwaveRadiation (bounds, num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Longwave radiation transfer through canopy using Norman (1979)
    !
    ! !USES:
    use clm_varcon, only : sb
    use decompMod, only : bounds_type
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : emg
    use MLclm_varpar, only : isun, isha, nlevmlcan
    use MLMathToolsMod, only : tridiag
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filter           ! Number of patches in filter
    integer, intent(in) :: filter(:)            ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                              ! Filter index
    integer  :: p                               ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                              ! Aboveground layer index
    integer  :: icm1                            ! Layer below ic (ic-1)
    real(r8) :: sumabs                          ! Absorbed radiation for energy conservation check
    real(r8) :: err                             ! Error check
    real(r8) :: omega                           ! Leaf scattering coefficient
    real(r8) :: rho                             ! Leaf reflectance
    real(r8) :: tau                             ! Leaf transmittance
    real(r8) :: trand                           ! Term for longwave radiation transmitted by layer
    real(r8) :: refld                           ! Term for longwave radiation reflected by layer
    real(r8) :: lw_source_sun                   ! Longwave radiation emitted by sunlit leaf (W/m2)
    real(r8) :: lw_source_sha                   ! Longwave radiation emitted by shaded leaf (W/m2)
    real(r8) :: lw_source(nlevmlcan)            ! Longwave radiation emitted by leaf layer (W/m2)
    integer  :: m                               ! Index to the tridiagonal matrix
    real(r8) :: aic, bic                        ! Intermediate terms for tridiagonal matrix
    real(r8) :: eic, fic                        ! Intermediate terms for tridiagonal matrix
    integer, parameter :: neq = (nlevmlcan+1)*2 ! Number of tridiagonal equations to solve
    real(r8) :: atri(neq), btri(neq)            ! Entries in tridiagonal matrix
    real(r8) :: ctri(neq), dtri(neq)            ! Entries in tridiagonal matrix
    real(r8) :: utri(neq)                       ! Tridiagonal solution
    real(r8) :: lwabs                           ! Absorbed longwave flux (W/m2 ground)
    real(r8) :: lwup_layer(bounds%begp:bounds%endp,0:nlevmlcan) ! Upward longwave flux above canopy layer (W/m2 ground)
    real(r8) :: lwdn_layer(bounds%begp:bounds%endp,0:nlevmlcan) ! Downward longwave flux onto canopy layer (W/m2 ground)
    !---------------------------------------------------------------------

    associate ( &
                                                   ! *** Input ***
    emleaf   => pftcon%emleaf                 , &  ! CLM (new): Leaf emissivity (-)
    lwsky    => mlcanopy_inst%lwsky_forcing   , &  ! Atmospheric longwave radiation (W/m2)
    ncan     => mlcanopy_inst%ncan_canopy     , &  ! Number of layers
    ntop     => mlcanopy_inst%ntop_canopy     , &  ! Index for top leaf layer
    nbot     => mlcanopy_inst%nbot_canopy     , &  ! Index for bottom leaf layer
    tg       => mlcanopy_inst%tg_soil         , &  ! Soil surface temperature (K)
    dpai     => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    fracsun  => mlcanopy_inst%fracsun_profile , &  ! Canopy layer sunlit fraction (-)
    td       => mlcanopy_inst%td_profile      , &  ! Canopy layer transmittance of diffuse radiation (-)
    tleaf    => mlcanopy_inst%tleaf_leaf      , &  ! Leaf temperature (K)
                                                   ! *** Output ***
    lwup     => mlcanopy_inst%lwup_canopy     , &  ! Upward longwave radiation above canopy (W/m2)
    lwveg    => mlcanopy_inst%lwveg_canopy    , &  ! Absorbed longwave radiation: vegetation (W/m2)
    lwsoi    => mlcanopy_inst%lwsoi_soil      , &  ! Absorbed longwave radiation: ground (W/m2)
    lwleaf   => mlcanopy_inst%lwleaf_leaf       &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Zero out radiative fluxes for all layers

       lwup_layer(p,0) = 0._r8
       lwdn_layer(p,0) = 0._r8

       do ic = 1, ncan(p)
          lwup_layer(p,ic) = 0._r8
          lwdn_layer(p,ic) = 0._r8
          lwleaf(p,ic,isun) = 0._r8
          lwleaf(p,ic,isha) = 0._r8
       end do

       ! Leaf scattering coefficient

       omega = 1._r8 - emleaf(patch%itype(p))

       ! Terms for longwave radiation reflected and transmitted by a layer:
       ! intercepted radiation is reflected

       rho = omega 
       tau = 0._r8

       ! Terms for longwave radiation reflected and transmitted by a layer:
       ! intercepted radiation is both reflected and transmitted

!      rho = omega * 0.5_r8
!      tau = omega * 0.5_r8

       ! Emitted longwave radiation is weighted average of sunlit and shaded leaves

       do ic = nbot(p), ntop(p)
          lw_source_sun = emleaf(patch%itype(p)) * sb * tleaf(p,ic,isun)**4
          lw_source_sha = emleaf(patch%itype(p)) * sb * tleaf(p,ic,isha)**4
          lw_source(ic) = (lw_source_sun * fracsun(p,ic) + lw_source_sha * (1._r8 - fracsun(p,ic))) &
                        * (1._r8 - td(p,ic))
       end do

       ! Set up and solve tridiagonal system of equations for upward and
       ! downward fluxes. There are two equations for each leaf layer and
       ! the soil. These equations are referenced by "m". The first
       ! equation is the upward flux and the second equation is the downward flux.

       m = 0

       ! Soil: upward flux

       m = m + 1
       atri(m) = 0._r8
       btri(m) = 1._r8
       ctri(m) = -(1._r8 - emg)
       dtri(m) = emg * sb * tg(p)**4

       ! Soil: downward flux

       refld = (1._r8 - td(p,nbot(p))) * rho
       trand = (1._r8 - td(p,nbot(p))) * tau + td(p,nbot(p))
       aic = refld - trand * trand / refld
       bic = trand / refld

       m = m + 1
       atri(m) = -aic
       btri(m) = 1._r8
       ctri(m) = -bic
       dtri(m) = (1._r8 - bic) * lw_source(nbot(p))

       ! Leaf layers, excluding top layer

       do ic = nbot(p), ntop(p)-1

          ! Upward flux

          refld = (1._r8 - td(p,ic)) * rho
          trand = (1._r8 - td(p,ic)) * tau + td(p,ic)
          fic = refld - trand * trand / refld
          eic = trand / refld

          m = m + 1
          atri(m) = -eic
          btri(m) = 1._r8
          ctri(m) = -fic
          dtri(m) = (1._r8 - eic) * lw_source(ic)

          ! Downward flux

          refld = (1._r8 - td(p,ic+1)) * rho
          trand = (1._r8 - td(p,ic+1)) * tau + td(p,ic+1)
          aic = refld - trand * trand / refld
          bic = trand / refld

          m = m + 1
          atri(m) = -aic
          btri(m) = 1._r8
          ctri(m) = -bic
          dtri(m) = (1._r8 - bic) * lw_source(ic+1)

       end do

       ! Top canopy layer: upward flux

       ic = ntop(p)
       refld = (1._r8 - td(p,ic)) * rho
       trand = (1._r8 - td(p,ic)) * tau + td(p,ic)
       fic = refld - trand * trand / refld
       eic = trand / refld

       m = m + 1
       atri(m) = -eic
       btri(m) = 1._r8
       ctri(m) = -fic
       dtri(m) = (1._r8 - eic) * lw_source(ic)

       ! Top canopy layer: downward flux

       m = m + 1
       atri(m) = 0._r8
       btri(m) = 1._r8
       ctri(m) = 0._r8
       dtri(m) = lwsky(p)

       ! Solve tridiagonal system of equations for upward and downward fluxes

       call tridiag (atri, btri, ctri, dtri, utri, m)

       ! Now copy the solution (utri) to the upward (lwup_layer) and downward (lwdn_layer)
       ! fluxes for each layer:
       ! lwup_layer = Upward longwave flux above layer
       ! lwdn_layer = Downward longwave flux onto layer

       m = 0

       ! Soil fluxes

       m = m + 1
       lwup_layer(p,0) = utri(m)
       m = m + 1
       lwdn_layer(p,0) = utri(m)

       ! Leaf layer fluxes

       do ic = nbot(p), ntop(p)
          m = m + 1
          lwup_layer(p,ic) = utri(m)
          m = m + 1
          lwdn_layer(p,ic) = utri(m)
       end do

       ! Absorbed longwave radiation for ground (soil)

       lwsoi(p) = lwdn_layer(p,0) - lwup_layer(p,0)

       ! Leaf layer fluxes

       lwveg(p) = 0._r8

       do ic = nbot(p), ntop(p)

          ! Absorbed longwave radiation for layer. Note special case for first
          ! leaf layer, where the upward flux from below is from the ground.
          ! The ground is ic=0, but nbot-1 will not equal 0 if there are lower
          ! canopy layers without leaves.

          if (ic == nbot(p)) then
             icm1 = 0
          else
             icm1 = ic - 1
          end if
          lwabs = emleaf(patch%itype(p)) * (lwdn_layer(p,ic)+lwup_layer(p,icm1)) * (1._r8 - td(p,ic)) &
                - 2._r8 * lw_source(ic)
          lwleaf(p,ic, isun) = lwabs / dpai(p,ic)
          lwleaf(p,ic, isha) = lwabs / dpai(p,ic)

          ! Sum longwave radiation absorbed by vegetation

          lwveg(p) = lwveg(p) + lwabs

       end do

       ! Canopy emitted longwave radiation

       lwup(p) = lwup_layer(p,ntop(p))

       ! Conservation check for total radiation balance: absorbed = incoming - outgoing

       sumabs = lwsky(p) - lwup(p)
       err = sumabs - (lwveg(p) + lwsoi(p))
       if (abs(err) > 1.e-03_r8) then
          call endrun (msg='ERROR: LongwaveRadiation: total longwave conservation error')
       end if

    end do

    end associate
  end subroutine LongwaveRadiation

end module MLLongwaveRadiationMod
