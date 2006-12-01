#include <misc.h>
#include <preproc.h>

module Biogeophysics2Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Biogeophysics2Mod
!
! !DESCRIPTION:
! Performs the calculation of soil/snow and ground temperatures
! and updates surface fluxes based on the new ground temperature.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Biogeophysics2   ! Calculate soil/snow and ground temperatures
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
! !IROUTINE: Biogeophysics2
!
! !INTERFACE:
  subroutine Biogeophysics2 (lbc, ubc, lbp, ubp, num_nolakec, &
             filter_nolakec, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! This is the main subroutine to execute the calculation of soil/snow and
! ground temperatures and update surface fluxes based on the new ground
! temperature
!
! Calling sequence is:
! Biogeophysics2:             surface biogeophysics driver
!    -> SoilTemperature:      soil/snow and ground temperatures
!          -> SoilTermProp    thermal conductivities and heat capacities
!          -> Tridiagonal     tridiagonal matrix solution
!          -> PhaseChange     phase change of liquid/ice contents
!
! (1) Snow and soil temperatures
!     o The volumetric heat capacity is calculated as a linear combination
!       in terms of the volumetric fraction of the constituent phases.
!     o The thermal conductivity of soil is computed from
!       the algorithm of Johansen (as reported by Farouki 1981), and the
!       conductivity of snow is from the formulation used in
!       SNTHERM (Jordan 1991).
!     o Boundary conditions:
!       F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
!     o Soil / snow temperature is predicted from heat conduction
!       in 10 soil layers and up to 5 snow layers.
!       The thermal conductivities at the interfaces between two
!       neighboring layers (j, j+1) are derived from an assumption that
!       the flux across the interface is equal to that from the node j
!       to the interface and the flux from the interface to the node j+1.
!       The equation is solved using the Crank-Nicholson method and
!       results in a tridiagonal system equation.
!
! (2) Phase change (see PhaseChange.F90)
!
! !USES:
    use clmtype
    use clm_atmlnd        , only : clm_a2l
    use clm_time_manager  , only : get_step_size
    use clm_varcon        , only : hvap, cpair, grav, vkc, tfrz, sb
    use clm_varpar        , only : nlevsno, nlevsoi, max_pft_per_col
    use SoilTemperatureMod, only : SoilTemperature
    use subgridAveMod     , only : p2c
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! Migrated to clm2.0 by Keith Oleson and Mariana Vertenstein
! Migrated to clm2.1 new data structures by Peter Thornton and M. Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)         ! pft's column index
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    real(r8), pointer :: pwtgcell(:)        ! pft's weight relative to corresponding column
    integer , pointer :: npfts(:)           ! column's number of pfts 
    integer , pointer :: pfti(:)            ! column's beginning pft index 
    integer , pointer :: snl(:)             ! number of snow layers
    logical , pointer :: do_capsnow(:)      ! true => do snow capping
    real(r8), pointer :: forc_lwrad(:)      ! downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: emg(:)             ! ground emissivity
    real(r8), pointer :: htvp(:)            ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: t_grnd(:)          ! ground temperature (Kelvin)
    integer , pointer :: frac_veg_nosno(:)  ! fraction of vegetation not covered by snow (0 OR 1 now) [-]
    real(r8), pointer :: cgrnds(:)          ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl(:)          ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: dlrad(:)           ! downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad(:)           ! upward longwave radiation above the canopy [W/m2]
    real(r8), pointer :: eflx_sh_veg(:)     ! sensible heat flux from leaves (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_veg(:)   ! vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:)   ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can(:)   ! evaporation from leaves and stems (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtcol(:)           ! pft weight relative to column
    real(r8), pointer :: tssbef(:,:)        ! soil/snow temperature before update
    real(r8), pointer :: t_soisno(:,:)      ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)    ! ice lens (kg/m2) (new)
    real(r8), pointer :: h2osoi_liq(:,:)    ! liquid water (kg/m2) (new)
! 
! local pointers to implicit inout arguments
!
    real(r8), pointer :: eflx_sh_grnd(:)    ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:)   ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_snowcap(:)    ! excess precipitation due to snow capping (mm H2O /s)
! 
! local pointers to implicit out arguments
! 
    real(r8), pointer :: dt_grnd(:)         ! change in t_grnd, last iteration (Kelvin)
    real(r8), pointer :: eflx_soil_grnd(:)  ! soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_sh_tot(:)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_tot(:)   ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: eflx_lh_tot(:)     ! total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: qflx_evap_grnd(:)  ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow(:)   ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow(:)   ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)   ! ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: eflx_lwrad_out(:)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lh_vege(:)    ! veg evaporation heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vegt(:)    ! veg transpiration heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_grnd(:)    ! ground evaporation heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: errsoi_pft(:)      ! pft-level soil/lake energy conservation error (W/m**2)
    real(r8), pointer :: errsoi_col(:)      ! column-level soil/lake energy conservation error (W/m**2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: p,c,g,j,pi,l         ! indices
    integer  :: fc,fp                ! lake filtered column and pft indices
    real(r8) :: dtime                ! land model time step (sec)
    real(r8) :: egsmax(lbc:ubc)      ! max. evaporation which soil can provide at one time step
    real(r8) :: egirat(lbc:ubc)      ! ratio of topsoil_evap_tot : egsmax
    real(r8) :: tinc(lbc:ubc)        ! temperature difference of two time step
    real(r8) :: xmf(lbc:ubc)         ! total latent heat of phase change of ground water
    real(r8) :: sumwt(lbc:ubc)       ! temporary
    real(r8) :: evaprat(lbp:ubp)     ! ratio of qflx_evap_soi/topsoil_evap_tot
    real(r8) :: save_qflx_evap_soi   ! temporary storage for qflx_evap_soi
    real(r8) :: topsoil_evap_tot(lbc:ubc)          ! column-level total evaporation from top soil layer
    real(r8) :: fact(lbc:ubc, -nlevsno+1:nlevsoi)  ! used in computing tridiagonal matrix
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_lwrad => clm_a2l%forc_lwrad

    ! Assign local pointers to derived subtypes components (column-level)

    npfts      => clm3%g%l%c%npfts
    pfti       => clm3%g%l%c%pfti
    snl        => clm3%g%l%c%cps%snl
    do_capsnow => clm3%g%l%c%cps%do_capsnow
    htvp       => clm3%g%l%c%cps%htvp
    emg        => clm3%g%l%c%cps%emg
    t_grnd     => clm3%g%l%c%ces%t_grnd
    dt_grnd    => clm3%g%l%c%ces%dt_grnd
    t_soisno   => clm3%g%l%c%ces%t_soisno
    tssbef     => clm3%g%l%c%ces%tssbef
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    errsoi_col => clm3%g%l%c%cebal%errsoi

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn        => clm3%g%l%c%p%column
    pgridcell      => clm3%g%l%c%p%gridcell
    pwtgcell       => clm3%g%l%c%p%wtgcell
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    sabg           => clm3%g%l%c%p%pef%sabg
    dlrad          => clm3%g%l%c%p%pef%dlrad
    ulrad          => clm3%g%l%c%p%pef%ulrad
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    eflx_sh_veg    => clm3%g%l%c%p%pef%eflx_sh_veg
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_evap_veg  => clm3%g%l%c%p%pwf%qflx_evap_veg
    qflx_tran_veg  => clm3%g%l%c%p%pwf%qflx_tran_veg
    qflx_evap_can  => clm3%g%l%c%p%pwf%qflx_evap_can
    qflx_snowcap   => clm3%g%l%c%p%pwf%qflx_snowcap
    qflx_evap_tot  => clm3%g%l%c%p%pwf%qflx_evap_tot
    qflx_evap_grnd => clm3%g%l%c%p%pwf%qflx_evap_grnd
    qflx_sub_snow  => clm3%g%l%c%p%pwf%qflx_sub_snow
    qflx_dew_snow  => clm3%g%l%c%p%pwf%qflx_dew_snow
    qflx_dew_grnd  => clm3%g%l%c%p%pwf%qflx_dew_grnd
    eflx_soil_grnd => clm3%g%l%c%p%pef%eflx_soil_grnd
    eflx_sh_tot    => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_lh_tot    => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out
    eflx_lwrad_net => clm3%g%l%c%p%pef%eflx_lwrad_net
    eflx_lh_vege   => clm3%g%l%c%p%pef%eflx_lh_vege
    eflx_lh_vegt   => clm3%g%l%c%p%pef%eflx_lh_vegt
    eflx_lh_grnd   => clm3%g%l%c%p%pef%eflx_lh_grnd
    cgrnds         => clm3%g%l%c%p%pef%cgrnds
    cgrndl         => clm3%g%l%c%p%pef%cgrndl
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    errsoi_pft     => clm3%g%l%c%p%pebal%errsoi
    wtcol          => clm3%g%l%c%p%wtcol

    ! Get step size

    dtime = get_step_size()

    ! Determine soil temperatures including surface soil temperature

    call SoilTemperature(lbc, ubc, num_nolakec, filter_nolakec, xmf , fact)

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       j = snl(c)+1

       ! Calculate difference in soil temperature from last time step, for
       ! flux corrections

       tinc(c) = t_soisno(c,j) - tssbef(c,j)

       ! Determine ratio of topsoil_evap_tot

       egsmax(c) = (h2osoi_ice(c,j)+h2osoi_liq(c,j)) / dtime

       ! added to trap very small negative soil water,ice

       if (egsmax(c) < 0._r8) then
          egsmax(c) = 0._r8
       end if
    end do

    ! A preliminary pft loop to determine if corrections are required for
    ! excess evaporation from the top soil layer... Includes new logic
    ! to distribute the corrections between pfts on the basis of their
    ! evaporative demands.
    ! egirat holds the ratio of demand to availability if demand is
    ! greater than availability, or 1.0 otherwise.
    ! Correct fluxes to present soil temperature

!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       eflx_sh_grnd(p) = eflx_sh_grnd(p) + tinc(c)*cgrnds(p)
       qflx_evap_soi(p) = qflx_evap_soi(p) + tinc(c)*cgrndl(p)
    end do

    ! Set the column-average qflx_evap_soi as the weighted average over all pfts
    ! but only count the pfts that are evaporating

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       topsoil_evap_tot(c) = 0._r8
       sumwt(c) = 0._r8
    end do

    do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             if (pwtgcell(p)>0._r8) then
                topsoil_evap_tot(c) = topsoil_evap_tot(c) + qflx_evap_soi(p) * wtcol(p)
             end if
          end if
       end do
    end do

    ! Calculate ratio for rescaling pft-level fluxes to meet availability

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       if (topsoil_evap_tot(c) > egsmax(c)) then
          egirat(c) = (egsmax(c)/topsoil_evap_tot(c))
       else
          egirat(c) = 1.0_r8
       end if
    end do

!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)
       j = snl(c)+1

       ! Correct soil fluxes for possible evaporation in excess of top layer water
       ! excess energy is added to the sensible heat flux from soil

       if (egirat(c) < 1.0_r8) then
          save_qflx_evap_soi = qflx_evap_soi(p)
          qflx_evap_soi(p) = qflx_evap_soi(p) * egirat(c)
          eflx_sh_grnd(p) = eflx_sh_grnd(p) + (save_qflx_evap_soi - qflx_evap_soi(p))*htvp(c)
       end if

       ! Ground heat flux

       eflx_soil_grnd(p) = sabg(p) + dlrad(p) + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) &
            - emg(c)*sb*tssbef(c,j)**3*(tssbef(c,j) + 4._r8*tinc(c)) &
            - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))

       ! Total fluxes (vegetation + ground)

       eflx_sh_tot(p) = eflx_sh_veg(p) + eflx_sh_grnd(p)
       qflx_evap_tot(p) = qflx_evap_veg(p) + qflx_evap_soi(p)
       eflx_lh_tot(p)= hvap*qflx_evap_veg(p) + htvp(c)*qflx_evap_soi(p)

       ! Assign ground evaporation to sublimation from soil ice or to dew
       ! on snow or ground

       qflx_evap_grnd(p) = 0._r8
       qflx_sub_snow(p) = 0._r8
       qflx_dew_snow(p) = 0._r8
       qflx_dew_grnd(p) = 0._r8

       if (qflx_evap_soi(p) >= 0._r8) then
          ! for evaporation partitioning between liquid evap and ice sublimation, 
	  ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
	  if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0.) then
             qflx_evap_grnd(p) = max(qflx_evap_soi(p)*(h2osoi_liq(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
	  else
	     qflx_evap_grnd(p) = 0.
	  end if
          qflx_sub_snow(p) = qflx_evap_soi(p) - qflx_evap_grnd(p)
       else
          if (t_grnd(c) < tfrz) then
             qflx_dew_snow(p) = abs(qflx_evap_soi(p))
          else
             qflx_dew_grnd(p) = abs(qflx_evap_soi(p))
          end if
       end if

       ! Update the pft-level qflx_snowcap
       ! This was moved in from Hydrology2 to keep all pft-level
       ! calculations out of Hydrology2

       if (snl(c) < 0 .and. do_capsnow(c)) then
          qflx_snowcap(p) = qflx_snowcap(p) + qflx_dew_snow(p) + qflx_dew_grnd(p)
       end if

       ! Outgoing long-wave radiation from vegetation + ground
       ! For conservation we put the increase of ground longwave to outgoing

       eflx_lwrad_out(p) = ulrad(p) &
            + (1-frac_veg_nosno(p))*(1._r8-emg(c))*forc_lwrad(g) &
            + (1-frac_veg_nosno(p))*emg(c)*sb * tssbef(c,j)**4 &
            + 4._r8*emg(c)*sb*tssbef(c,j)**3*tinc(c)

       ! Variables needed by history tape

       qflx_evap_can(p)  = qflx_evap_veg(p) - qflx_tran_veg(p)
       eflx_lh_vege(p)   = (qflx_evap_veg(p) - qflx_tran_veg(p)) * hvap
       eflx_lh_vegt(p)   = qflx_tran_veg(p) * hvap
       eflx_lh_grnd(p)   = qflx_evap_soi(p) * htvp(c)
       eflx_lwrad_net(p) = eflx_lwrad_out(p) - forc_lwrad(g)

    end do

    ! Soil Energy balance check

!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       errsoi_pft(p) = eflx_soil_grnd(p) - xmf(c)
    end do
    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fp = 1,num_nolakep
          p = filter_nolakep(fp)
          c = pcolumn(p)
          if (j >= snl(c)+1) then
             errsoi_pft(p) = errsoi_pft(p) - (t_soisno(c,j)-tssbef(c,j))/fact(c,j)
          end if
       end do
    end do

    ! lake balance for errsoi is not over pft
    ! therefore obtain column-level radiative temperature

    call p2c(num_nolakec, filter_nolakec, errsoi_pft, errsoi_col)

  end subroutine Biogeophysics2

end module Biogeophysics2Mod
