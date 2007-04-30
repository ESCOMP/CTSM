#include <misc.h>
#include <preproc.h>

module SoilHydrologyMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SoilHydrologyMod
!
! !DESCRIPTION:
! Calculate soil hydrology
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRunoff  ! Calculate surface runoff
  public :: Infiltration   ! Calculate infiltration into surface soil layer
  public :: SoilWater      ! Calculate soil hydrology
  public :: Drainage       ! Calculate subsurface drainage
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 04/25/07 Keith Oleson: CLM3.5 hydrology
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SurfaceRunoff
!
! !INTERFACE:
  subroutine SurfaceRunoff (lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
       vol_liq, icefrac)
!
! !DESCRIPTION:
! Calculate surface runoff
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon, only : denice, denh2o, wimp
    use clm_varpar, only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                   ! column bounds
    integer , intent(in)  :: lbp, ubp                   ! pft bounds   
    integer , intent(in)  :: num_soilc                  ! number of column soil points in column filter
    integer , intent(in)  :: filter_soilc(ubc-lbc+1)    ! column filter for soil points
    real(r8), intent(out) :: vol_liq(lbc:ubc,1:nlevsoi) ! partial volume of liquid water in layer
    real(r8), intent(out) :: icefrac(lbc:ubc,1:nlevsoi) ! fraction of ice in layer (-)
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/26/02, Peter Thornton: Migrated to new data structures.
! 4/26/05, David Lawrence: Made surface runoff for dry soils a function
!   of rooting fraction in top three soil layers.
! 04/25/07  Keith Oleson: Completely new routine for CLM3.5 hydrology
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    real(r8), pointer :: qflx_top_soil(:)  !net water input into soil from top (mm/s)
    real(r8), pointer :: watsat(:,:)       !volumetric soil water at saturation (porosity)
    real(r8), pointer :: hkdepth(:)        !decay factor (m)
    real(r8), pointer :: zwt(:)            !water table depth (m)
    real(r8), pointer :: fcov(:)           !fractional area with water table at surface
    real(r8), pointer :: dz(:,:)           !layer depth (m)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
    real(r8), pointer :: wtfact(:)         !maximum saturated fraction for a gridcell
    real(r8), pointer :: hksat(:,:)        !hydraulic conductivity at saturation (mm H2O /s)
    real(r8), pointer :: bsw(:,:)          !Clapp and Hornberger "b"
    real(r8), pointer :: sucsat(:,:)       !minimum soil suction (mm)
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: qflx_surf(:)      !surface runoff (mm H2O /s)
    real(r8), pointer :: eff_porosity(:,:) !effective porosity = porosity - vol_ice
    real(r8), pointer :: fracice(:,:)      !fractional impermeability (-)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c,j,fc,g                   !indices
    real(r8) :: vol_ice(lbc:ubc,1:nlevsoi) !partial volume of ice lens in layer
    real(r8) :: fff(lbc:ubc)               !decay factor (m-1)
    real(r8) :: s1                         !variable to calculate qinmax
    real(r8) :: su                         !variable to calculate qinmax
    real(r8) :: v                          !variable to calculate qinmax
    real(r8) :: qinmax                     !maximum infiltration capacity (mm/s)

!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    qflx_top_soil => clm3%g%l%c%cwf%qflx_top_soil
    qflx_surf     => clm3%g%l%c%cwf%qflx_surf
    watsat        => clm3%g%l%c%cps%watsat
    hkdepth       => clm3%g%l%c%cps%hkdepth
    dz            => clm3%g%l%c%cps%dz
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    fcov          => clm3%g%l%c%cws%fcov
    eff_porosity  => clm3%g%l%c%cps%eff_porosity
    wtfact        => clm3%g%l%c%cps%wtfact
    zwt           => clm3%g%l%c%cws%zwt
    fracice       => clm3%g%l%c%cps%fracice
    hksat         => clm3%g%l%c%cps%hksat
    bsw           => clm3%g%l%c%cps%bsw
    sucsat        => clm3%g%l%c%cps%sucsat

    do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)

          ! Porosity of soil, partial volume of ice and liquid, fraction of ice in each layer,
          ! fractional impermeability
   
          vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
          eff_porosity(c,j) = watsat(c,j)-vol_ice(c,j)
          vol_liq(c,j) = min(eff_porosity(c,j), h2osoi_liq(c,j)/(dz(c,j)*denh2o))

          icefrac(c,j) = min(1._r8,h2osoi_ice(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j)))

          fracice(c,j) = max(0._r8,exp(-3._r8*(1._r8-icefrac(c,j)))- exp(-3._r8))

       end do
    end do

    ! Saturated fraction

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       fff(c)  = 1._r8 / hkdepth(c)
       fcov(c) = (1._r8 - fracice(c,1)) * wtfact(c) * exp(-0.5_r8*fff(c)*zwt(c)) + fracice(c,1)
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)

       ! Maximum infiltration capacity
       s1        = max(0.01_r8,vol_liq(c,1)/max(wimp,eff_porosity(c,1)))
       su        = max(0._r8,(s1-fcov(c)*1._r8) / (1._r8-fcov(c)))
       v         = -bsw(c,1)*sucsat(c,1)/(0.5_r8*dz(c,1)*1000._r8)
       qinmax    = (1._r8+v*(su-1._r8))*hksat(c,1)

       ! Surface runoff
       qflx_surf(c) =  fcov(c) * qflx_top_soil(c) + &
                       (1._r8-fcov(c)) * max(0._r8, qflx_top_soil(c)-qinmax)

    end do

  end subroutine SurfaceRunoff

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Infiltration
!
! !INTERFACE:
  subroutine Infiltration(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! Calculate infiltration into surface soil layer (minus the evaporation)
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                   ! column bounds
    integer, intent(in) :: num_soilc                  ! number of column soil points in column filter
    integer, intent(in) :: filter_soilc(ubc-lbc+1)    ! column filter for soil points
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/27/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: snl(:)            !minus number of snow layers
    real(r8), pointer :: qflx_top_soil(:)  !net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_surf(:)      !surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_evap_grnd(:) !ground surface evaporation rate (mm H2O/s) [+]
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: qflx_infl(:)      !infiltration (mm H2O /s)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: c, fc    !indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    snl            => clm3%g%l%c%cps%snl
    qflx_top_soil  => clm3%g%l%c%cwf%qflx_top_soil
    qflx_surf      => clm3%g%l%c%cwf%qflx_surf
    qflx_infl      => clm3%g%l%c%cwf%qflx_infl
    qflx_evap_grnd => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd

    ! Infiltration into surface soil layer (minus the evaporation)

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       if (snl(c) >= 0) then
          qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) - qflx_evap_grnd(c)
       else
          qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c)
       end if
    end do

  end subroutine Infiltration

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilWater
!
! !INTERFACE:
  subroutine SoilWater(lbc, ubc, num_soilc, &
       filter_soilc, vol_liq, dwat, hk, dhkdw)
!
! !DESCRIPTION:
! Soil hydrology
! Soil moisture is predicted from a 10-layer model (as with soil
! temperature), in which the vertical soil moisture transport is governed
! by infiltration, runoff, gradient diffusion, gravity, and root
! extraction through canopy transpiration.  The net water applied to the
! surface layer is the snowmelt plus precipitation plus the throughfall
! of canopy dew minus surface runoff and evaporation.
! CLM3.5 uses a zero-flow bottom boundary condition.
!
! The vertical water flow in an unsaturated porous media is described by
! Darcy's law, and the hydraulic conductivity and the soil negative
! potential vary with soil water content and soil texture based on the work
! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
! integrated over the layer thickness, in which the time rate of change in
! water mass must equal the net flow across the bounding interface, plus the
! rate of internal source or sink. The terms of water flow across the layer
! interfaces are linearly expanded by using first-order Taylor expansion.
! The equations result in a tridiagonal system equation.
!
! Note: length units here are all millimeter
! (in temperature subroutine uses same soil layer
! structure required but lengths are m)
!
! Richards equation:
!
! d wat      d     d wat d psi
! ----- = - -- [ k(----- ----- - 1) ] + S
!   dt      dz       dz  d wat
!
! where: wat = volume of water per volume of soil (mm**3/mm**3)
! psi = soil matrix potential (mm)
! dt  = time step (s)
! z   = depth (mm)
! dz  = thickness (mm)
! qin = inflow at top (mm h2o /s)
! qout= outflow at bottom (mm h2o /s)
! s   = source/sink flux (mm h2o /s)
! k   = hydraulic conductivity (mm h2o /s)
!
!                       d qin                  d qin
! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
!                       d wat(j-1)             d wat(j)
!                ==================|=================
!                                  < qin
!
!                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
!
!                                  > qout
!                ==================|=================
!                        d qout               d qout
! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
!                        d wat(j)             d wat(j+1)
!
!
! Solution: linearize k and psi about d wat and use tridiagonal
! system of equations to solve for d wat,
! where for layer j
!
!
! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon    , only : wimp
    use clm_varpar    , only : nlevsoi, max_pft_per_col
    use shr_const_mod , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use TridiagonalMod, only : Tridiagonal
    use clm_time_manager  , only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                   ! column bounds
    integer , intent(in)  :: num_soilc                  ! number of column soil points in column filter
    integer , intent(in)  :: filter_soilc(ubc-lbc+1)    ! column filter for soil points
    real(r8), intent(in)  :: vol_liq(lbc:ubc,1:nlevsoi) ! soil water per unit volume [mm/mm]
    real(r8), intent(out) :: dwat(lbc:ubc,1:nlevsoi)    ! change of soil water [m3/m3]
    real(r8), intent(out) :: hk(lbc:ubc,1:nlevsoi)      ! hydraulic conductivity [mm h2o/s]
    real(r8), intent(out) :: dhkdw(lbc:ubc,1:nlevsoi)   ! d(hk)/d(vol_liq)
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/27/02, Peter Thornton: Migrated to new data structures. Includes
! treatment of multiple PFTs on a single soil column.
! 04/25/07 Keith Oleson: CLM3.5 hydrology
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: npfts(:)             ! column's number of pfts - ADD
    real(r8), pointer :: pwtcol(:)            ! weight relative to column for each pft
    real(r8), pointer :: pwtgcell(:)          ! weight relative to gridcell for each pft
    real(r8), pointer :: z(:,:)               ! layer depth (m)
    real(r8), pointer :: dz(:,:)              ! layer thickness (m)
    real(r8), pointer :: smpmin(:)            ! restriction for min of soil potential (mm)
    real(r8), pointer :: qflx_infl(:)         ! infiltration (mm H2O /s)
    real(r8), pointer :: qflx_tran_veg_pft(:) ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg_col(:) ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: eff_porosity(:,:)    ! effective porosity = porosity - vol_ice
    real(r8), pointer :: watsat(:,:)          ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: hksat(:,:)           ! hydraulic conductivity at saturation (mm H2O /s)
    real(r8), pointer :: bsw(:,:)             ! Clapp and Hornberger "b"
    real(r8), pointer :: sucsat(:,:)          ! minimum soil suction (mm)
    real(r8), pointer :: t_soisno(:,:)        ! soil temperature (Kelvin)
    real(r8), pointer :: rootr_pft(:,:)       ! effective fraction of roots in each soil layer
    integer , pointer :: pfti(:)              ! beginning pft index for each column
    real(r8), pointer :: fracice(:,:)         ! fractional impermeability (-)
    real(r8), pointer :: h2osoi_vol(:,:)      ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
!
! local pointers to original implicit inout arguments
!
    real(r8), pointer :: h2osoi_liq(:,:)      ! liquid water (kg/m2)
!
! local pointer s to original implicit out arguments
!
    real(r8), pointer :: rootr_col(:,:)       ! effective fraction of roots in each soil layer
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: p,c,fc,j                  ! do loop indices
    integer  :: jtop(lbc:ubc)             ! top level at each column
    real(r8) :: dtime                     ! land model time step (sec)
    real(r8) :: amx(lbc:ubc,1:nlevsoi)    ! "a" left off diagonal of tridiagonal matrix
    real(r8) :: bmx(lbc:ubc,1:nlevsoi)    ! "b" diagonal column for tridiagonal matrix
    real(r8) :: cmx(lbc:ubc,1:nlevsoi)    ! "c" right off diagonal tridiagonal matrix
    real(r8) :: rmx(lbc:ubc,1:nlevsoi)    ! "r" forcing term of tridiagonal matrix
    real(r8) :: zmm(lbc:ubc,1:nlevsoi)    ! layer depth [mm]
    real(r8) :: dzmm(lbc:ubc,1:nlevsoi)   ! layer thickness [mm]
    real(r8) :: den                       ! used in calculating qin, qout
    real(r8) :: dqidw0                    ! d(qin)/d(vol_liq(i-1))
    real(r8) :: dqidw1                    ! d(qin)/d(vol_liq(i))
    real(r8) :: dqodw1                    ! d(qout)/d(vol_liq(i))
    real(r8) :: dqodw2                    ! d(qout)/d(vol_liq(i+1))
    real(r8) :: dsmpdw(lbc:ubc,1:nlevsoi) ! d(smp)/d(vol_liq)
    real(r8) :: num                       ! used in calculating qin, qout
    real(r8) :: qin                       ! flux of water into soil layer [mm h2o/s]
    real(r8) :: qout                      ! flux of water out of soil layer [mm h2o/s]
    real(r8) :: s_node                    ! soil wetness
    real(r8) :: s1                        ! "s" at interface of layer
    real(r8) :: s2                        ! k*s**(2b+2)
    real(r8) :: smp(lbc:ubc,1:nlevsoi)    ! soil matrix potential [mm]
    real(r8) :: sdamp                     ! extrapolates soiwat dependence of evaporation
    integer  :: pi                        ! pft index
    real(r8) :: temp(lbc:ubc)             ! accumulator for rootr weighting
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    npfts             => clm3%g%l%c%npfts
    z                 => clm3%g%l%c%cps%z
    dz                => clm3%g%l%c%cps%dz
    smpmin            => clm3%g%l%c%cps%smpmin
    watsat            => clm3%g%l%c%cps%watsat
    hksat             => clm3%g%l%c%cps%hksat
    bsw               => clm3%g%l%c%cps%bsw
    sucsat            => clm3%g%l%c%cps%sucsat
    eff_porosity      => clm3%g%l%c%cps%eff_porosity
    rootr_col         => clm3%g%l%c%cps%rootr_column
    t_soisno          => clm3%g%l%c%ces%t_soisno
    h2osoi_liq        => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol        => clm3%g%l%c%cws%h2osoi_vol
    qflx_infl         => clm3%g%l%c%cwf%qflx_infl
    fracice           => clm3%g%l%c%cps%fracice
    qflx_tran_veg_col => clm3%g%l%c%cwf%pwf_a%qflx_tran_veg
    pfti              => clm3%g%l%c%pfti

    ! Assign local pointers to derived type members (pft-level)

    qflx_tran_veg_pft => clm3%g%l%c%p%pwf%qflx_tran_veg
    rootr_pft         => clm3%g%l%c%p%pps%rootr
    pwtcol            => clm3%g%l%c%p%wtcol
    pwtgcell          => clm3%g%l%c%p%wtgcell

    ! Get time step

    dtime = get_step_size()

    ! Because the depths in this routine are in mm, use local
    ! variable arrays instead of pointers

    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          zmm(c,j) = z(c,j)*1.e3_r8
          dzmm(c,j) = dz(c,j)*1.e3_r8
       end do
    end do

    ! First step is to calculate the column-level effective rooting
    ! fraction in each soil layer. This is done outside the usual
    ! PFT-to-column averaging routines because it is not a simple
    ! weighted average of the PFT level rootr arrays. Instead, the
    ! weighting depends on both the per-unit-area transpiration
    ! of the PFT and the PFTs area relative to all PFTs.

    temp(:) = 0._r8

    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          rootr_col(c,j) = 0._r8
       end do
    end do

    do pi = 1,max_pft_per_col
       do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
          do fc = 1, num_soilc
             c = filter_soilc(fc)
             if (pi <= npfts(c)) then
                p = pfti(c) + pi - 1
                if (pwtgcell(p)>0._r8) then
                   rootr_col(c,j) = rootr_col(c,j) + rootr_pft(p,j) * qflx_tran_veg_pft(p) * pwtcol(p)
                end if
             end if
          end do
       end do
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          if (pi <= npfts(c)) then
             p = pfti(c) + pi - 1
             if (pwtgcell(p)>0._r8) then
                temp(c) = temp(c) + qflx_tran_veg_pft(p) * pwtcol(p)
             end if
          end if
       end do
    end do

    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          if (temp(c) /= 0._r8) then
             rootr_col(c,j) = rootr_col(c,j)/temp(c)
          end if
       end do
    end do

    ! Hydraulic conductivity and soil matric potential and their derivatives

    sdamp = 0._r8
    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)

          s1 = 0.5_r8*(h2osoi_vol(c,j) + h2osoi_vol(c,min(nlevsoi, j+1))) / &
               (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
          s1 = min(1._r8, s1)
          s2 = hksat(c,j)*s1**(2._r8*bsw(c,j)+2._r8)

          hk(c,j) = (1._r8-0.5_r8*(fracice(c,j)+fracice(c,min(nlevsoi, j+1))))*s1*s2

          dhkdw(c,j) = (1._r8-0.5_r8*(fracice(c,j)+fracice(c,min(nlevsoi, j+1))))* &
                       (2._r8*bsw(c,j)+3._r8)*s2*0.5_r8/watsat(c,j)
          if(j == nlevsoi) dhkdw(c,j) = dhkdw(c,j) * 2._r8

          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          s_node = min(1.0_r8, s_node)

          smp(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
          smp(c,j) = max(smpmin(c), smp(c,j))

          dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/(s_node*watsat(c,j))

       end do
    end do

    ! Set up r, a, b, and c vectors for tridiagonal solution

    ! Node j=1 (top)

    j = 1
!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       qin    = qflx_infl(c)
       den    = (zmm(c,j+1)-zmm(c,j))
       num    = (smp(c,j+1)-smp(c,j)) - den
       qout   = -hk(c,j)*num/den
       dqodw1 = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
       dqodw2 = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
       rmx(c,j) =  qin - qout - qflx_tran_veg_col(c) * rootr_col(c,j)
       amx(c,j) =  0._r8
       bmx(c,j) =  dzmm(c,j)*(sdamp+1._r8/dtime) + dqodw1
       cmx(c,j) =  dqodw2
    end do

    ! Nodes j=2 to j=nlevsoi-1

    do j = 2, nlevsoi - 1
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          den    = (zmm(c,j) - zmm(c,j-1))
          num    = (smp(c,j)-smp(c,j-1)) - den
          qin    = -hk(c,j-1)*num/den
          dqidw0 = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
          dqidw1 = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
          den    = (zmm(c,j+1)-zmm(c,j))
          num    = (smp(c,j+1)-smp(c,j)) - den
          qout   = -hk(c,j)*num/den
          dqodw1 = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
          dqodw2 = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
          rmx(c,j) =  qin - qout - qflx_tran_veg_col(c)*rootr_col(c,j)
          amx(c,j) = -dqidw0
          bmx(c,j) =  dzmm(c,j)/dtime - dqidw1 + dqodw1
          cmx(c,j) =  dqodw2
       end do
    end do

    ! Node j=nlevsoi (bottom)

    j = nlevsoi
!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       den    = (zmm(c,j) - zmm(c,j-1))
       num    = (smp(c,j)-smp(c,j-1)) - den
       qin    = -hk(c,j-1)*num/den
       dqidw0 = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
       dqidw1 = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
       qout   =  0._r8  ! zero-flow bottom boundary condition
       dqodw1 =  0._r8  ! zero-flow bottom boundary condition
       rmx(c,j) =  qin - qout - qflx_tran_veg_col(c)*rootr_col(c,j)
       amx(c,j) = -dqidw0
       bmx(c,j) =  dzmm(c,j)/dtime - dqidw1 + dqodw1
       cmx(c,j) =  0._r8
    end do

    ! Solve for dwat

    jtop(:) = 1
    call Tridiagonal(lbc, ubc, 1, nlevsoi, jtop, num_soilc, filter_soilc, &
                     amx, bmx, cmx, rmx, dwat(lbc:ubc,1:nlevsoi))

    ! Renew the mass of liquid water

    do j= 1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          h2osoi_liq(c,j) = h2osoi_liq(c,j) + dwat(c,j)*dzmm(c,j)
       end do
    end do

  end subroutine SoilWater

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Drainage
!
! !INTERFACE:
  subroutine Drainage(lbc, ubc, num_soilc, filter_soilc, vol_liq, hk, &
          icefrac)
!
! !DESCRIPTION:
! Calculate subsurface drainage
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_time_manager, only : get_step_size
    use clm_varcon  , only: pondmx, tfrz
    use clm_varpar  , only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                   ! column bounds
    integer , intent(in) :: num_soilc                  ! number of column soil points in column filter
    integer , intent(in) :: filter_soilc(ubc-lbc+1)    ! column filter for soil points
    real(r8), intent(in) :: vol_liq(lbc:ubc,1:nlevsoi) ! partial volume of liquid water in layer
    real(r8), intent(in) :: hk(lbc:ubc,1:nlevsoi)      ! hydraulic conductivity (mm h2o/s)
    real(r8), intent(in) :: icefrac(lbc:ubc,1:nlevsoi) ! fraction of ice in layer
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 4/26/05, Peter Thornton and David Lawrence: Turned off drainage from
! middle soil layers for both wet and dry fractions.
! 04/25/07  Keith Oleson: Completely new routine for CLM3.5 hydrology
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: snl(:)            !number of snow layers
    real(r8), pointer :: qflx_snowcap(:)   !excess precipitation due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)  !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow(:)  !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow(:)  !sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: dz(:,:)           !layer depth (m)
    real(r8), pointer :: bsw(:,:)          !Clapp and Hornberger "b"
    real(r8), pointer :: eff_porosity(:,:) !effective porosity = porosity - vol_ice
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: hksat(:,:)        !hydraulic conductivity at saturation (mm H2O /s)
    real(r8), pointer :: sucsat(:,:)       !minimum soil suction (mm)
    real(r8), pointer :: z(:,:)            !layer depth (m)
    real(r8), pointer :: zi(:,:)           !interface level below a "z" level (m)
    real(r8), pointer :: watsat(:,:)       !volumetric soil water at saturation (porosity)
    real(r8), pointer :: hkdepth(:)        !decay factor (m)
    real(r8), pointer :: zwt(:)            !water table depth (m)
    real(r8), pointer :: wa(:)             !water in the unconfined aquifer (mm)
    real(r8), pointer :: wt(:)             !total water storage (unsaturated soil water + groundwater) (mm)
    real(r8), pointer :: qcharge(:)        !aquifer recharge rate (mm/s)
!
! local pointers to original implicit inout arguments
!
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: qflx_drain(:)     !sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)     !qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
    real(r8), pointer :: eflx_impsoil(:)   !implicit evaporation for soil temperature equation
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c,j,fc                   !indices
    real(r8) :: dtime                    !land model time step (sec)
    real(r8) :: xs(lbc:ubc)              !water needed to bring soil moisture to watmin (mm)
    real(r8) :: dzmm(lbc:ubc,1:nlevsoi)  !layer thickness (mm)
    real(r8) :: watmin                   !minimum soil moisture (mm)
    integer  :: jwt(lbc:ubc)             !index of the soil layer right above the water table (-)
    real(r8) :: rsub_bot(lbc:ubc)        !subsurface runoff - bottom drainage (mm/s)
    real(r8) :: rsub_sat(lbc:ubc)        !subsurface runoff - saturation excess (mm/s)
    real(r8) :: rsub_top(lbc:ubc)        !subsurface runoff - topographic control (mm/s)
    real(r8) :: fff(lbc:ubc)             !decay factor (m-1)
    real(r8) :: xsi(lbc:ubc)             !excess soil water above saturation at layer i (mm)
    real(r8) :: xs1(lbc:ubc)             !excess soil water above saturation at layer 1 (mm)
    real(r8) :: smpfz(1:nlevsoi)         !matric potential of layer right above water table (mm)
    real(r8) :: wtsub                    !summation of hk*dzmm for layers below water table (mm**2/s)
    real(r8) :: rous                     !aquifer yield (-)
    real(r8) :: wh                       !smpfz(jwt)-z(jwt) (mm)
    real(r8) :: wh_zwt                   !water head at the water table depth (mm)
    real(r8) :: ws                       !summation of pore space of layers below water table (mm)
    real(r8) :: s_node                   !soil wetness (-)
    real(r8) :: dzsum                    !summation of dzmm of layers below water table (mm)
    real(r8) :: icefracsum               !summation of icefrac*dzmm of layers below water table (-)
    real(r8) :: fracice_rsub(lbc:ubc)    !fractional impermeability of soil layers (-)
    real(r8) :: ka                       !hydraulic conductivity of the aquifer (mm/s)
    real(r8) :: dza                      !fff*(zwt-z(jwt)) (-)
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    snl           => clm3%g%l%c%cps%snl
    dz            => clm3%g%l%c%cps%dz
    bsw           => clm3%g%l%c%cps%bsw
    t_soisno      => clm3%g%l%c%ces%t_soisno
    hksat         => clm3%g%l%c%cps%hksat
    sucsat        => clm3%g%l%c%cps%sucsat
    z             => clm3%g%l%c%cps%z
    zi            => clm3%g%l%c%cps%zi
    watsat        => clm3%g%l%c%cps%watsat
    hkdepth       => clm3%g%l%c%cps%hkdepth
    zwt           => clm3%g%l%c%cws%zwt
    wa            => clm3%g%l%c%cws%wa
    wt            => clm3%g%l%c%cws%wt
    qcharge       => clm3%g%l%c%cws%qcharge
    eff_porosity  => clm3%g%l%c%cps%eff_porosity
    qflx_snowcap  => clm3%g%l%c%cwf%pwf_a%qflx_snowcap
    qflx_dew_grnd => clm3%g%l%c%cwf%pwf_a%qflx_dew_grnd
    qflx_dew_snow => clm3%g%l%c%cwf%pwf_a%qflx_dew_snow
    qflx_sub_snow => clm3%g%l%c%cwf%pwf_a%qflx_sub_snow
    qflx_drain    => clm3%g%l%c%cwf%qflx_drain
    qflx_qrgwl    => clm3%g%l%c%cwf%qflx_qrgwl
    eflx_impsoil  => clm3%g%l%c%cef%eflx_impsoil
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice

    ! Get time step

    dtime = get_step_size()

    ! Convert layer thicknesses from m to mm

    do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          dzmm(c,j) = dz(c,j)*1.e3_r8
       end do
    end do

    ! Initial set

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       qflx_drain(c) = 0._r8 
       rsub_bot(c)   = 0._r8
       rsub_sat(c)   = 0._r8
       rsub_top(c)   = 0._r8
       fracice_rsub(c) = 0._r8
    end do

    ! The layer index of the first unsaturated layer, i.e., the layer right above
    ! the water table

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       jwt(c) = nlevsoi
       do j = 2,nlevsoi
          if(zwt(c) <= zi(c,j)) then
             jwt(c) = j-1
             exit
          end if
       enddo
    end do

    ! Topographic runoff
!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       fff(c)         = 1._r8/ hkdepth(c)
       dzsum = 0._r8
       icefracsum = 0._r8
       do j = jwt(c), nlevsoi
          dzsum  = dzsum + dzmm(c,j)
          icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
       end do
       fracice_rsub(c) = max(0._r8,exp(-3._r8*(1._r8-(icefracsum/dzsum)))- exp(-3._r8))
       rsub_top(c)    = (1._r8 - fracice_rsub(c)) * 4.5e-4_r8 * exp(-fff(c)*zwt(c))
    end do

    rous = 0.2_r8

    ! Water table calculation

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)

       ! Matric potential at the layer above the water table
       s_node = vol_liq(c,jwt(c))/eff_porosity(c,jwt(c))
       s_node = max(s_node, 0.01_r8)
       s_node = min(1.0_r8, s_node)
       smpfz(jwt(c)) = -sucsat(c,jwt(c))*s_node**(-bsw(c,jwt(c)))
       smpfz(jwt(c)) = max(-80000.0_r8, smpfz(jwt(c)))

       ! Water head at the water table and hydraulic conductivity of the aquifer
       if(jwt(c) == nlevsoi) then
          dza = fff(c)*(zwt(c)-z(c,jwt(c)))
          ka = hk(c,jwt(c)) * (1.0_r8-exp(-dza))/dza
          wh_zwt  = -zwt(c) * 1000._r8
       else
          ka = hk(c,jwt(c))
          wh_zwt  = -sucsat(c,jwt(c)+1) - zwt(c) * 1000._r8
       endif

       ! Recharge rate qcharge to groundwater (positive to aquifer)
       wh      = smpfz(jwt(c))  - z(c,jwt(c))*1.e3_r8
       qcharge(c) = -ka * (wh_zwt-wh)  /((zwt(c)-z(c,jwt(c)))*1000._r8)

       ! To limit qcharge  (for the first several timesteps)
       qcharge(c) = max(-10.0_r8/dtime,qcharge(c))
       qcharge(c) = min( 10.0_r8/dtime,qcharge(c))

       ! Water storage in aquifer + soil
       wt(c)  = wt(c) + (qcharge(c) - rsub_top(c)) * dtime

       if(jwt(c) == nlevsoi) then             ! water table is below the soil column
          wa(c)  = wa(c) + (qcharge(c) -rsub_top(c)) * dtime
          wt(c)  = wa(c)
          zwt(c)     = (zi(c,nlevsoi) + 25._r8) - wa(c)/1000._r8/rous
          h2osoi_liq(c,nlevsoi) = h2osoi_liq(c,nlevsoi) - qcharge(c) * dtime
          h2osoi_liq(c,nlevsoi) = h2osoi_liq(c,nlevsoi) + max(0._r8,(wa(c)-5000._r8))
          wa(c)  = min(wa(c), 5000._r8)
       else                                ! water table within soil layers
          if (jwt(c) == nlevsoi-1) then       ! water table within bottom soil layer

             zwt(c) = zi(c,nlevsoi)- (wt(c)-rous*1000._r8*25._r8) /eff_porosity(c,nlevsoi)/1000._r8

          else                                   ! water table within soil layers 1-9
             ws = 0._r8   ! water used to fill soil air pores regardless of water content
             do j = jwt(c)+2,nlevsoi
               ws = ws + eff_porosity(c,j) * 1000._r8 * dz(c,j)
             enddo
             zwt(c) = zi(c,jwt(c)+1)-(wt(c)-rous*1000_r8*25._r8-ws) /eff_porosity(c,jwt(c)+1)/1000._r8
          endif

          wtsub = 0._r8
          do j = jwt(c)+1, nlevsoi
             wtsub = wtsub + hk(c,j)*dzmm(c,j)
          end do

          ! Remove subsurface runoff
          do j = jwt(c)+1, nlevsoi 
             h2osoi_liq(c,j) = h2osoi_liq(c,j) - rsub_top(c)*dtime*hk(c,j)*dzmm(c,j)/wtsub
          end do
       end if

       zwt(c) = max(0.05_r8,zwt(c))
       zwt(c) = min(80._r8,zwt(c))

    end do

    !  excessive water above saturation added to the above unsaturated layer like a bucket
    !  if column fully saturated, excess water goes to runoff

    do j = nlevsoi,2,-1
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          xsi(c)            = max(h2osoi_liq(c,j)-eff_porosity(c,j)*dzmm(c,j),0._r8)
          h2osoi_liq(c,j)   = min(eff_porosity(c,j)*dzmm(c,j), h2osoi_liq(c,j))
          h2osoi_liq(c,j-1) = h2osoi_liq(c,j-1) + xsi(c)
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       xs1(c)          = max(h2osoi_liq(c,1)-(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice(c,1)),0._r8)
       h2osoi_liq(c,1) = min(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice(c,1), h2osoi_liq(c,1))
       rsub_sat(c)     = xs1(c) / dtime
    end do

    ! Limit h2osoi_liq to be greater than or equal to watmin.
    ! Get water needed to bring h2osoi_liq equal watmin from lower layer.
    ! If insufficient water in soil layers, get from aquifer water

    watmin = 0.01_r8

    do j = 1, nlevsoi-1
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          if (h2osoi_liq(c,j) < 0._r8) then
             xs(c) = watmin - h2osoi_liq(c,j)
          else
             xs(c) = 0._r8
          end if
          h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
          h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
       end do
    end do

    j = nlevsoi
!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       if (h2osoi_liq(c,j) < watmin) then
          xs(c) = watmin-h2osoi_liq(c,j)
       else
          xs(c) = 0._r8
       end if
       h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
       wa(c) = wa(c) - xs(c)
       wt(c) = wt(c) - xs(c)
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)

       ! Sub-surface runoff and drainage

       qflx_drain(c) = rsub_sat(c) + rsub_top(c)

       ! Set imbalance for snow capping

       qflx_qrgwl(c) = qflx_snowcap(c)

       ! Implicit evaporation term is now zero

       eflx_impsoil(c) = 0._r8

       ! Renew the ice and liquid mass due to condensation

       if (snl(c)+1 >= 1) then
          h2osoi_liq(c,1) = h2osoi_liq(c,1) + qflx_dew_grnd(c) * dtime
          h2osoi_ice(c,1) = h2osoi_ice(c,1) + (qflx_dew_snow(c) * dtime)
          if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
             qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
             h2osoi_ice(c,1) = 0._r8
          else
             h2osoi_ice(c,1) = h2osoi_ice(c,1) - (qflx_sub_snow(c) * dtime)
          end if
       end if
    end do

  end subroutine Drainage

end module SoilHydrologyMod
