#include <misc.h>
#include <preproc.h>

module Biogeophysics1Mod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Biogeophysics1Mod
!
! !DESCRIPTION:
! Performs calculation of leaf temperature and surface fluxes.
! Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: Biogeophysics1   ! Calculate leaf temperature and surface fluxes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Biogeophysics1
!
! !INTERFACE:
  subroutine Biogeophysics1(lbg, ubg, lbc, ubc, lbp, ubp, &
       num_nolakec, filter_nolakec, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! This is the main subroutine to execute the calculation of leaf temperature
! and surface fluxes. Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! Calling sequence is:
! Biogeophysics1:           surface biogeophysics driver
!  -> QSat:                 saturated vapor pressure, specific humidity, and
!                           derivatives at ground surface and derivatives at
!                           leaf surface using updated leaf temperature
! Leaf temperature
! Foliage energy conservation is given by the foliage energy budget
! equation:
!                Rnet - Hf - LEf = 0
! The equation is solved by Newton-Raphson iteration, in which this
! iteration includes the calculation of the photosynthesis and
! stomatal resistance, and the integration of turbulent flux profiles.
! The sensible and latent heat transfer between foliage and atmosphere
! and ground is linked by the equations:
!                Ha = Hf + Hg and Ea = Ef + Eg
!
! !USES:
    use clmtype
    use clm_atmlnd         , only : clm_a2l
    use clm_varcon         , only : denh2o, denice, roverg, hvap, hsub, &
                                    istice, istwet, zlnd, zsno, spval
    use clm_varpar         , only : nlevsoi, nlevsno
    use QSatMod            , only : QSat
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg                    ! gridcell-index bounds
    integer, intent(in) :: lbc, ubc                    ! column-index bounds
    integer, intent(in) :: lbp, ubp                    ! pft-index bounds
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
    integer , pointer :: ivt(:)           !pft vegetation type
    integer , pointer :: ityplun(:)       !landunit type
    integer , pointer :: clandunit(:)     !column's landunit index
    integer , pointer :: cgridcell(:)     !column's gridcell index
    real(r8), pointer :: forc_pbot(:)     !atmospheric pressure (Pa)
    real(r8), pointer :: forc_q(:)        !atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_t(:)        !atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_hgt_t(:)    !observational height of temperature [m]
    real(r8), pointer :: forc_th(:)       !atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_u(:)        !atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v(:)        !atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: smpmin(:)        !restriction for min of soil potential (mm)
    integer , pointer :: snl(:)           !number of snow layers
    real(r8), pointer :: frac_sno(:)      !fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: h2osno(:)        !snow water (mm H2O)
    real(r8), pointer :: elai(:)          !one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)          !one-sided stem area index with burying by snow
    real(r8), pointer :: z0mr(:)          !ratio of momentum roughness length to canopy top height (-)
    real(r8), pointer :: displar(:)       !ratio of displacement height to canopy top height (-)
    real(r8), pointer :: htop(:)          !canopy top (m)
    real(r8), pointer :: dz(:,:)          !layer depth (m)
    real(r8), pointer :: t_soisno(:,:)    !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_liq(:,:)  !liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)  !ice lens (kg/m2)
    real(r8), pointer :: watsat(:,:)      !volumetric soil water at saturation (porosity)
    real(r8), pointer :: sucsat(:,:)      !minimum soil suction (mm)
    real(r8), pointer :: bsw(:,:)         !Clapp and Hornberger "b"
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: t_grnd(:)        !ground temperature (Kelvin)
    real(r8), pointer :: qg(:)            !ground specific humidity [kg/kg]
    real(r8), pointer :: dqgdT(:)         !d(qg)/dT
    real(r8), pointer :: emg(:)           !ground emissivity
    real(r8), pointer :: htvp(:)          !latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: beta(:)          !coefficient of convective velocity [-]
    real(r8), pointer :: zii(:)           !convective boundary height [m]
    real(r8), pointer :: thm(:)           !intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8), pointer :: thv(:)           !virtual potential temperature (kelvin)
    real(r8), pointer :: z0mg(:)          !roughness length over ground, momentum [m]
    real(r8), pointer :: z0hg(:)          !roughness length over ground, sensible heat [m]
    real(r8), pointer :: z0qg(:)          !roughness length over ground, latent heat [m]
    real(r8), pointer :: emv(:)           !vegetation emissivity
    real(r8), pointer :: z0m(:)           !momentum roughness length (m)
    real(r8), pointer :: displa(:)        !displacement height (m)
    real(r8), pointer :: z0mv(:)          !roughness length over vegetation, momentum [m]
    real(r8), pointer :: z0hv(:)          !roughness length over vegetation, sensible heat [m]
    real(r8), pointer :: z0qv(:)          !roughness length over vegetation, latent heat [m]
    real(r8), pointer :: eflx_sh_tot(:)   !total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot(:)   !total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: eflx_sh_veg(:)   !sensible heat flux from leaves (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_tot(:) !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_veg(:) !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:) !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: cgrnd(:)         !deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: cgrnds(:)        !deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl(:)        !deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8) ,pointer :: tssbef(:,:)      !soil/snow temperature before update
    real(r8) ,pointer :: soilalpha(:)     !factor that reduces ground saturated specific humidity (-)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: g,l,c,p !indices
    integer  :: j       !soil/snow level index
    integer  :: fp      !lake filter pft index
    integer  :: fc      !lake filter column index
    real(r8) :: qred    !soil surface relative humidity
    real(r8) :: avmuir  !ir inverse optical depth per unit leaf area
    real(r8) :: eg      !water vapor pressure at temperature T [pa]
    real(r8) :: qsatg   !saturated humidity [kg/kg]
    real(r8) :: degdT   !d(eg)/dT
    real(r8) :: qsatgdT !d(qsatg)/dT
    real(r8) :: fac     !soil wetness of surface layer
    real(r8) :: psit    !negative potential of soil
    real(r8) :: hr      !relative humidity
    real(r8) :: wx      !partial volume of ice and water of surface layer
!------------------------------------------------------------------------------

   ! Assign local pointers to derived type members (gridcell-level)

    forc_hgt_t    => clm_a2l%forc_hgt_t
    forc_pbot     => clm_a2l%forc_pbot
    forc_q        => clm_a2l%forc_q
    forc_t        => clm_a2l%forc_t
    forc_th       => clm_a2l%forc_th
    forc_u        => clm_a2l%forc_u
    forc_v        => clm_a2l%forc_v

    ! Assign local pointers to derived type members (landunit-level)

    ityplun       => clm3%g%l%itype

    ! Assign local pointers to derived type members (column-level)

    cgridcell     => clm3%g%l%c%gridcell
    clandunit     => clm3%g%l%c%landunit
    beta          => clm3%g%l%c%cps%beta
    dqgdT         => clm3%g%l%c%cws%dqgdT
    emg           => clm3%g%l%c%cps%emg
    frac_sno      => clm3%g%l%c%cps%frac_sno
    h2osno        => clm3%g%l%c%cws%h2osno
    htvp          => clm3%g%l%c%cps%htvp
    qg            => clm3%g%l%c%cws%qg
    smpmin        => clm3%g%l%c%cps%smpmin
    snl           => clm3%g%l%c%cps%snl
    t_grnd        => clm3%g%l%c%ces%t_grnd
    thm           => clm3%g%l%c%ces%thm
    thv           => clm3%g%l%c%ces%thv
    z0hg          => clm3%g%l%c%cps%z0hg
    z0mg          => clm3%g%l%c%cps%z0mg
    z0qg          => clm3%g%l%c%cps%z0qg
    zii           => clm3%g%l%c%cps%zii
    bsw           => clm3%g%l%c%cps%bsw
    dz            => clm3%g%l%c%cps%dz
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    soilalpha     => clm3%g%l%c%cws%soilalpha
    sucsat        => clm3%g%l%c%cps%sucsat
    t_soisno      => clm3%g%l%c%ces%t_soisno
    tssbef        => clm3%g%l%c%ces%tssbef
    watsat        => clm3%g%l%c%cps%watsat

    ! Assign local pointers to derived type members (pft-level)

    ivt           => clm3%g%l%c%p%itype
    elai          => clm3%g%l%c%p%pps%elai
    esai          => clm3%g%l%c%p%pps%esai
    htop          => clm3%g%l%c%p%pps%htop
    emv           => clm3%g%l%c%p%pps%emv
    z0m           => clm3%g%l%c%p%pps%z0m
    displa        => clm3%g%l%c%p%pps%displa
    z0mv          => clm3%g%l%c%p%pps%z0mv
    z0hv          => clm3%g%l%c%p%pps%z0hv
    z0qv          => clm3%g%l%c%p%pps%z0qv
    eflx_sh_tot   => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_lh_tot   => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_sh_veg   => clm3%g%l%c%p%pef%eflx_sh_veg
    qflx_evap_tot => clm3%g%l%c%p%pwf%qflx_evap_tot
    qflx_evap_veg => clm3%g%l%c%p%pwf%qflx_evap_veg
    qflx_tran_veg => clm3%g%l%c%p%pwf%qflx_tran_veg
    cgrnd         => clm3%g%l%c%p%pef%cgrnd
    cgrnds        => clm3%g%l%c%p%pef%cgrnds
    cgrndl        => clm3%g%l%c%p%pef%cgrndl

    ! Assign local pointers to derived type members (ecophysiological)

    z0mr          => pftcon%z0mr
    displar       => pftcon%displar

    do j = -nlevsno+1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          tssbef(c,j) = t_soisno(c,j)
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       g = cgridcell(c)

       ! begin calculations that relate only to the column level
       ! Ground and soil temperatures from previous time step

       t_grnd(c) = t_soisno(c,snl(c)+1)

       ! Saturated vapor pressure, specific humidity and their derivatives
       ! at ground surface

       qred = 1._r8
       if (ityplun(l)/=istwet .AND. ityplun(l)/=istice) then
          wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)
          fac  = min(1._r8, wx/watsat(c,1))
          fac  = max( fac, 0.01_r8 )
          psit = -sucsat(c,1) * fac ** (-bsw(c,1))
          psit = max(smpmin(c), psit)
          hr   = exp(psit/roverg/t_grnd(c))
          qred = (1._r8-frac_sno(c))*hr + frac_sno(c)
          soilalpha(c) = qred
       else
          soilalpha(c) = spval
       end if

       call QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg, qsatgdT)

       qg(c) = qred*qsatg
       dqgdT(c) = qred*qsatgdT

       if (qsatg > forc_q(g) .and. forc_q(g) > qred*qsatg) then
          qg(c) = forc_q(g)
          dqgdT(c) = 0._r8
       end if

       ! Ground Emissivity

       if (h2osno(c)>0._r8 .or. ityplun(l)==istice) then
          emg(c) = 0.97_r8
       else
          emg(c) = 0.96_r8
       end if

       ! Latent heat. We arbitrarily assume that the sublimation occurs
       ! only as h2osoi_liq = 0

       htvp(c) = hvap
       if (h2osoi_liq(c,snl(c)+1) <= 0._r8 .and. h2osoi_ice(c,snl(c)+1) > 0._r8) htvp(c) = hsub

       ! Switch between vaporization and sublimation causes rapid solution
       ! separation in perturbation growth test

#if (defined PERGRO)
       htvp(c) = hvap
#endif

       ! Ground roughness lengths over non-lake columns (includes bare ground, ground
       ! underneath canopy, wetlands, etc.)

       if (frac_sno(c) > 0._r8) then
          z0mg(c) = zsno
       else
          z0mg(c) = zlnd
       end if
       z0hg(c) = z0mg(c)            ! initial set only
       z0qg(c) = z0mg(c)            ! initial set only

       ! Potential, virtual potential temperature, and wind speed at the
       ! reference height

       beta(c) = 1._r8
       zii(c)  = 1000._r8
       thm(c)  = forc_t(g) + 0.0098_r8*forc_hgt_t(g)
       thv(c)  = forc_th(g)*(1._r8+0.61_r8*forc_q(g))

    end do ! (end of columns loop)

    ! Initialization

!dir$ concurrent
!cdir nodep
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)

       ! Initial set (needed for history tape fields)

       eflx_sh_tot(p) = 0._r8
       eflx_lh_tot(p) = 0._r8
       eflx_sh_veg(p) = 0._r8
       qflx_evap_tot(p) = 0._r8
       qflx_evap_veg(p) = 0._r8
       qflx_tran_veg(p) = 0._r8

       ! Initial set for calculation

       cgrnd(p)  = 0._r8
       cgrnds(p) = 0._r8
       cgrndl(p) = 0._r8

       ! Vegetation Emissivity

       avmuir = 1._r8
       emv(p) = 1._r8-exp(-(elai(p)+esai(p))/avmuir)

       ! Roughness lengths over vegetation

       z0m(p)    = z0mr(ivt(p)) * htop(p)
       displa(p) = displar(ivt(p)) * htop(p)

       z0mv(p)   = z0m(p)
       z0hv(p)   = z0mv(p)
       z0qv(p)   = z0mv(p)

    end do

  end subroutine Biogeophysics1

end module Biogeophysics1Mod
