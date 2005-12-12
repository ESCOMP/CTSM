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
       zwice, vol_liq, s, zwt, fcov)
!
! !DESCRIPTION:
! Calculate surface runoff
! The original code was provide by Robert E. Dickinson based on
! following clues:  exponential decrease of Ksat, a water table
! level determination level including highland and lowland levels
! and fractional area of wetland (water table above the surface).
! Runoff is parameterized from the lowlands in terms of precip
! incident on wet areas and a base flow, where these are estimated
! using ideas from TOPMODEL.
! The original scheme was modified by Z.-L. Yang and G.-Y. Niu,
! o  using a new method to determine water table depth and
!    the fractional wet area (fcov)
! o  computing runoff (surface and subsurface) from this
!    fraction and the remaining fraction (i.e. 1-fcov)
! o  for the 1-fcov part, using BATS1e method to compute
!    surface and subsurface runoff.
! The original code on soil moisture and runoff were provided by
! R. E. Dickinson in July 1996.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon, only : denice, denh2o, wimp
    use clm_varpar, only : nlevsoi, maxpatch_pft
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                   ! column bounds
    integer , intent(in)  :: lbp, ubp                   ! pft bounds   
    integer , intent(in)  :: num_soilc                  ! number of column soil points in column filter
    integer , intent(in)  :: filter_soilc(ubc-lbc+1)    ! column filter for soil points
    real(r8), intent(out) :: zwice(lbc:ubc)             ! the sum of ice mass of soil (kg/m2)
    real(r8), intent(out) :: vol_liq(lbc:ubc,1:nlevsoi) ! partial volume of liquid water in layer
    real(r8), intent(out) :: s(lbc:ubc,1:nlevsoi)       ! wetness of soil (including ice)
    real(r8), intent(out) :: zwt(lbc:ubc)               ! water table depth
    real(r8), intent(out) :: fcov(lbc:ubc)              ! fractional area with water table at surface
!
! !CALLED FROM:
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/26/02, Peter Thornton: Migrated to new data structures.
! 4/26/05, David Lawrence: Made surface runoff for dry soils a function
!   of rooting fraction in top three soil layers.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: cgridcell(:)      !gridcell index for each column
    real(r8), pointer :: wtfact(:)         !Fraction of model area with high water table
    real(r8), pointer :: qflx_top_soil(:)  !net water input into soil from top (mm/s)
    real(r8), pointer :: watsat(:,:)       !volumetric soil water at saturation (porosity)
    real(r8), pointer :: dz(:,:)           !layer depth (m)
    real(r8), pointer :: zi(:,:)           !interface level below a "z" level (m)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
    integer , pointer :: npfts(:)          !column's number of pfts - ADD   
    real(r8), pointer :: wtcol(:)          !weight relative to column for each pft 
    integer , pointer :: pfti(:)           !beginning pft index for each column 
    real(r8), pointer :: rootfr(:,:)       !fraction of roots in each soil layer 
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: qflx_surf(:)      !surface runoff (mm H2O /s)
    real(r8), pointer :: eff_porosity(:,:) !effective porosity = porosity - vol_ice
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: p,c,g,j,fc                 ! indices
    integer  :: pi                        ! pft index  
    real(r8):: vol_ice(lbc:ubc,1:nlevsoi) ! partial volume of ice lens in layer
    real(r8):: zmean(lbc:ubc)             ! The surface soil layers contributing to runoff
    real(r8):: wmean(lbc:ubc)             ! The averaged soil wetness in surface soil layers
    real(r8):: infil_fact_col(lbc:ubc)    ! Infiltration enhancement factor for pft's 
    real(r8):: infil_fact(lbp:ubp)        ! Infiltration enhancement factor 
    real(r8), parameter:: fz = 1._r8      ! coefficient for water table depth
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (gridcell-level)

    wtfact        => clm3%g%gps%wtfact

    ! Assign local pointers to derived subtype components (column-level)

    cgridcell     => clm3%g%l%c%gridcell
    qflx_top_soil => clm3%g%l%c%cwf%qflx_top_soil
    qflx_surf     => clm3%g%l%c%cwf%qflx_surf
    watsat        => clm3%g%l%c%cps%watsat
    dz            => clm3%g%l%c%cps%dz
    zi            => clm3%g%l%c%cps%zi
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    eff_porosity  => clm3%g%l%c%cps%eff_porosity
    npfts         => clm3%g%l%c%npfts   
    wtcol         => clm3%g%l%c%p%wtcol   
    pfti          => clm3%g%l%c%pfti   
    rootfr        => clm3%g%l%c%p%pps%rootfr

    ! Porosity of soil, partial volume of ice and liquid and
    ! water table depth

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       zwice(c) = 0._r8
       wmean(c) = 0._r8
    end do

    do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)

          ! Porosity of soil, partial volume of ice and liquid

          zwice(c) = zwice(c) + h2osoi_ice(c,j)
          vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
          eff_porosity(c,j) = watsat(c,j)-vol_ice(c,j)
          vol_liq(c,j) = min(eff_porosity(c,j), h2osoi_liq(c,j)/(dz(c,j)*denh2o))

          ! Wetness of soil

          s(c,j) = min(1._r8,(vol_ice(c,j)+vol_liq(c,j))/watsat(c,j))

          ! Averaged soil wetness in surface soil layers

          wmean(c) = wmean(c) + s(c,j)*dz(c,j)
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       g = cgridcell(c)

       ! Determine water table depth

       zwt(c) = fz * (zi(c,nlevsoi) - wmean(c))

       ! Determine saturation fraction

       fcov(c) = wtfact(g) * min(1._r8,exp(-zwt(c)))

       ! Re-initialize wmean and initialize zmean to zero for use below

       wmean(c) = 0._r8
       zmean(c) = 0._r8
       infil_fact_col(c) = 0._r8  
    end do

    do pi = 1,maxpatch_pft    
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc       
          c = filter_soilc(fc)    
          if (pi <= npfts(c)) then    
             p = pfti(c) + pi - 1    
             infil_fact(p) = 0._r8
             do j = 1, 3
                infil_fact(p)=rootfr(p,j) + infil_fact(p)
             end do
             infil_fact_col(c) = infil_fact_col(c) + infil_fact(p) * wtcol(p)    
          end if    
       end do    
    end do    

    do j = 1, 3
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          zmean(c) = zmean(c) + dz(c,j)
          wmean(c) = wmean(c) + s(c,j) * dz(c,j)
       end do
    end do

    ! If top soil layer is impermeable then all qflx_top_soil goes to surface runoff

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       wmean(c) = wmean(c) / zmean(c)
       if (eff_porosity(c,1) < wimp) then
          qflx_surf(c) =  max(0._r8, fcov(c) * qflx_top_soil(c)) + &
                          max(0._r8, (1._r8-fcov(c)) * qflx_top_soil(c))
       else
#if (defined CN)
          ! keep the saturated fraction qflx_over, but turn off the
          ! unsaturated fraction
			 qflx_surf(c) =  max(0._r8, fcov(c) * qflx_top_soil(c))
#else
          qflx_surf(c) =  max(0._r8, fcov(c) * qflx_top_soil(c)) + &
            max(0._r8,(1._r8-fcov(c))*  &
            min(1._r8,wmean(c)**(4+16*infil_fact_col(c)))*qflx_top_soil(c))  
#endif
       end if
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
! The original code on soil moisture and runoff were provided by
! R. E. Dickinson in July 1996.
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
    integer , pointer :: snl(:)           ! minus number of snow layers
    real(r8), pointer :: qflx_top_soil(:) ! net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_surf(:)     ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_evap_grnd(:)! ground surface evaporation rate (mm H2O/s) [+]
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: qflx_infl(:)     !infiltration (mm H2O /s)
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
    use time_manager  , only : get_step_size
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
    real(r8) :: dzmm(lbc:ubc,1:nlevsoi)   !layer thickness [mm]
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
    qflx_infl         => clm3%g%l%c%cwf%qflx_infl
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

    ! set initial values

    sdamp = 0._r8
    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)

          ! Set hydraulic conductivity to zero if effective porosity 5% in any of
          ! two neighbor layers or liquid content (theta) less than 0.001

          if (      (eff_porosity(c,j) < wimp) &
               .or. (eff_porosity(c,min(nlevsoi,j+1)) < wimp) &
               .or. (vol_liq(c,j) <= 1.e-3_r8)) then
             hk(c,j) = 0._r8
             dhkdw(c,j) = 0._r8
          else
             s1 = 0.5_r8*(vol_liq(c,j)+vol_liq(c,min(nlevsoi,j+1))) / &
                  (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi,j+1))))
             s2 = hksat(c,j)*s1**(2._r8*bsw(c,j)+2._r8)
             hk(c,j) = s1*s2
             dhkdw(c,j) = (2._r8*bsw(c,j)+3._r8)*s2*0.5_r8/watsat(c,j)
             if (j == nlevsoi) dhkdw(c,j) = dhkdw(c,j) * 2._r8
          end if

          ! Evaluate hydraulic conductivity, soil matric potential,
          ! d(smp)/d(vol_liq), and d(hk)/d(vol_liq).

          if (t_soisno(c,j) > SHR_CONST_TKFRZ) then

             s_node = max(vol_liq(c,j)/watsat(c,j), 0.01_r8)
             s_node = min(1.0_r8, s_node)
             smp(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
             smp(c,j) = max(smpmin(c), smp(c,j))        ! Limit soil suction
             dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/(s_node*watsat(c,j))

          else

             ! When ice is present, the matric potential is only related to temperature
             ! by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
             ! Unit 1 Joule = 1 (kg m2/s2), J/kg /(m/s2) ==> m ==> 1e3 mm

             smp(c,j) = 1.e3_r8 * SHR_CONST_LATICE/SHR_CONST_G * &
                  (t_soisno(c,j)-SHR_CONST_TKFRZ)/t_soisno(c,j)
             smp(c,j) = max(smpmin(c), smp(c,j))        ! Limit soil suction
             dsmpdw(c,j) = 0._r8

          end if
       end do
    end do

    ! Set up r, a, b, and c vectors for tridiagonal solution

    ! Node j=1

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

    ! Node j=nlevsoi

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
       qout   =  hk(c,j)
       dqodw1 =  dhkdw(c,j)
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
  subroutine Drainage(lbc, ubc, num_soilc, filter_soilc, &
       zwice, vol_liq, s, zwt, fcov, hk, dhkdw, dwat)
!
! !DESCRIPTION:
! Calculate subsurface drainage
! The original code was provide by Robert E. Dickinson based on
! following clues:  exponential decrease of Ksat, a water table
! level determination level including highland and lowland levels
! and fractional area of wetland (water table above the surface).
! Runoff is parameterized from the lowlands in terms of precip
! incident on wet areas and a base flow, where these are estimated
! using ideas from TOPMODEL.
! The original scheme was modified by Z.-L. Yang and G.-Y. Niu,
! *  using a new method to determine water table depth and
!    the fractional wet area (fcov)
! *  computing runoff (surface and subsurface) from this
!    fraction and the remaining fraction (i.e. 1-fcov)
! *  for the 1-fcov part, using BATS1e method to compute
!    surface and subsurface runoff.
! The original code on soil moisture and runoff were provided by
! R. E. Dickinson in July 1996.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use time_manager, only : get_step_size
    use clm_varcon  , only: pondmx
    use clm_varpar  , only : nlevsoi
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                   ! column bounds
    integer , intent(in) :: num_soilc                  ! number of column soil points in column filter
    integer , intent(in) :: filter_soilc(ubc-lbc+1)    ! column filter for soil points
    real(r8), intent(in) :: zwice(lbc:ubc)             ! the sum of ice mass of soil (kg/m2)
    real(r8), intent(in) :: vol_liq(lbc:ubc,1:nlevsoi) ! partial volume of liquid water in layer
    real(r8), intent(in) :: s(lbc:ubc,1:nlevsoi)       ! wetness of soil (including ice)
    real(r8), intent(in) :: zwt(lbc:ubc)               ! water table depth
    real(r8), intent(in) :: fcov(lbc:ubc)              ! fractional area with water table at surface
    real(r8), intent(in) :: hk(lbc:ubc,1:nlevsoi)      ! hydraulic conductivity (mm h2o/s)
    real(r8), intent(in) :: dhkdw(lbc:ubc,1:nlevsoi)   ! d(hk)/d(vol_liq)
    real(r8), intent(in) :: dwat(lbc:ubc,1:nlevsoi)    ! change in soil water
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 4/26/05, Peter Thornton and David Lawrence: Turned off drainage from
! middle soil layers for both wet and dry fractions.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer:: snl(:)            !number of snow layers
    real(r8), pointer:: qflx_snowcap(:)   !excess precipitation due to snow capping (mm H2O /s) [+]
    real(r8), pointer:: qflx_dew_grnd(:)  !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer:: qflx_dew_snow(:)  !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer:: qflx_sub_snow(:)  !sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer:: dz(:,:)           !layer depth (m)
    real(r8), pointer:: bsw(:,:)          !Clapp and Hornberger "b"
    real(r8), pointer:: eff_porosity(:,:) !effective porosity = porosity - vol_ice
!
! local pointers to original implicit inout arguments
!
    real(r8), pointer:: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer:: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to original implicit out arguments
!
    real(r8), pointer:: qflx_drain(:)     !sub-surface runoff (mm H2O /s)
    real(r8), pointer:: qflx_qrgwl(:)     !qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer:: eflx_impsoil(:)   !implicit evaporation for soil temperature equation
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c,j,fc                   !indices
    real(r8) :: dtime                    !land model time step (sec)
    real(r8) :: xs(lbc:ubc)              !excess soil water above saturation
    real(r8) :: dzmm(lbc:ubc,1:nlevsoi)  !layer thickness (mm)
    real(r8) :: watmin                   !minimum soil moisture
    real(r8) :: hksum(lbc:ubc)           !summation of hydraulic cond for layers 6->9
    real(r8) :: zsat(lbc:ubc)            !hydraulic conductivity weighted soil thickness
    real(r8) :: wsat(lbc:ubc)            !hydraulic conductivity weighted soil wetness
    real(r8) :: qflx_drain_wet(lbc:ubc)  !subsurface runoff from "wet" part (mm h2o/s)
    real(r8) :: qflx_drain_dry(lbc:ubc)  !subsurface runoff from "dry" part (mm h2o/s)
    real(r8) :: dzksum(lbc:ubc)          !hydraulic conductivity weighted soil thickness
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    snl           => clm3%g%l%c%cps%snl
    dz            => clm3%g%l%c%cps%dz
    bsw           => clm3%g%l%c%cps%bsw
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

    ! Streamflow and total runoff

    ! Convert layer thicknesses from m to mm

    do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          dzmm(c,j) = dz(c,j)*1.e3_r8
       end do
    end do

    ! The amount of streamflow is assumed maintained by flow from the
    ! lowland water table with different levels contributing according to
    ! their thickness and saturated hydraulic conductivity, i.e. a given
    ! layer below the water table interface loses water at a rate per unit
    ! depth given by qflx_drain*hk/(sum over all layers below this water table
    ! of hk*dz). Because this is a slow smooth process, and not strongly
    ! coupled to water in any one layer, it should remain stable for
    ! explicit time differencing. Hence, for simplicity it is removed
    ! explicitly prior to the main soil water calculation.
    ! Another assumption: no subsurface runoff for ice mixed soil
    ! Zong-Liang Yang & G.-Y. Niu

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       qflx_drain(c) = 0._r8          ! subsurface runoff
       qflx_drain_wet(c) = 0._r8      ! subsurface runoff
       qflx_drain_dry(c) = 0._r8      ! subsurface runoff
       hksum(c) = 0._r8
    end do

    do j = 6,nlevsoi-1
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          hksum(c) = hksum(c) + hk(c,j)
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       if (zwice(c) <= 0._r8 .AND. hksum(c) > 0._r8) then
          zsat(c) = 0._r8
          wsat(c) = 0._r8
          dzksum(c) = 0._r8
       end if
    end do

    do j = 6,nlevsoi-1
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          if (zwice(c) <= 0._r8 .AND. hksum(c) > 0._r8) then
             zsat(c) = zsat(c) + dz(c,j)*hk(c,j)
             wsat(c) = wsat(c) + s(c,j)*dz(c,j)*hk(c,j)
             dzksum(c) = dzksum(c) + hk(c,j)*dz(c,j)
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       if (zwice(c) <= 0._r8 .AND. hksum(c) > 0._r8) then
          wsat(c) = wsat(c) / zsat(c)
          qflx_drain_dry(c) = 0.0_r8
          qflx_drain_wet(c) = 0.0_r8
          qflx_drain(c) = qflx_drain_dry(c) + qflx_drain_wet(c)
       end if
    end do

    do j = 6, nlevsoi-1
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          if (zwice(c) <= 0._r8 .AND. hksum(c) > 0._r8) then
             h2osoi_liq(c,j) = h2osoi_liq(c,j) - dtime*qflx_drain(c)*dz(c,j)*hk(c,j)/dzksum(c)
          end if
       end do
    end do

    ! Limit h2osoi_liq to be greater than or equal to watmin.
    ! Get water needed to bring h2osoi_liq equal watmin from lower layer.

    watmin = 0.0_r8

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
       qflx_drain(c) = qflx_drain(c) - xs(c)/dtime
    end do

    ! Determine water in excess of saturation

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       xs(c) = max(0._r8, h2osoi_liq(c,1)-(pondmx + eff_porosity(c,1)*dzmm(c,1)))
       if (xs(c) > 0._r8) h2osoi_liq(c,1) = pondmx + eff_porosity(c,1)*dzmm(c,1)
    end do

    do j = 2,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          xs(c) = xs(c) + max(h2osoi_liq(c,j) - eff_porosity(c,j)*dzmm(c,j), 0._r8)  ! [mm]
          h2osoi_liq(c,j) = min(eff_porosity(c,j)*dzmm(c,j), h2osoi_liq(c,j))
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_soilc
       c = filter_soilc(fc)

       ! Sub-surface runoff and drainage

       qflx_drain(c) = qflx_drain(c) + xs(c)/dtime &
            + hk(c,nlevsoi) + dhkdw(c,nlevsoi)*dwat(c,nlevsoi) ! [mm/s]

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
