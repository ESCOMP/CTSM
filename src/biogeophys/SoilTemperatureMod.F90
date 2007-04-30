#include <misc.h>
#include <preproc.h>

module SoilTemperatureMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SoilTemperatureMod
!
! !DESCRIPTION:
! Calculates snow and soil temperatures including phase change
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilTemperature  ! Snow and soil temperatures including phase change
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilThermProp   ! Set therm conductivities and heat cap of snow/soil layers
  private :: PhaseChange     ! Calculation of the phase change within snow and soil layers
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
! !IROUTINE: SoilTemperature
!
! !INTERFACE:
  subroutine SoilTemperature(lbc, ubc, num_nolakec, filter_nolakec, &
                             xmf, fact)
!
! !DESCRIPTION:
! Snow and soil temperatures including phase change
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), and the
!   conductivity of snow is from the formulation used in
!   SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction
!   in 10 soil layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two
!   neighboring layers (j, j+1) are derived from an assumption that
!   the flux across the interface is equal to that from the node j
!   to the interface and the flux from the interface to the node j+1.
!   The equation is solved using the Crank-Nicholson method and
!   results in a tridiagonal system equation.
!
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd    , only : clm_a2l
    use clm_time_manager  , only : get_step_size
    use clm_varcon    , only : sb, capr, cnfac
    use clm_varpar    , only : nlevsno, nlevsoi, max_pft_per_col
    use TridiagonalMod, only : Tridiagonal
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                    ! column bounds
    integer , intent(in)  :: num_nolakec                 ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    real(r8), intent(out) :: xmf(lbc:ubc)                ! total latent heat of phase change of ground water
    real(r8), intent(out) :: fact(lbc:ubc, -nlevsno+1:nlevsoi) ! used in computing tridiagonal matrix
!
! !CALLED FROM:
! subroutine Biogeophysics2 in module Biogeophysics2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! Changed references for tg to t_grnd, for consistency with the
! rest of the code (tg eliminated as redundant)
! 2/14/02, Peter Thornton: Migrated to new data structures. Added pft loop
! in calculation of net ground heat flux.
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    integer , pointer :: npfts(:)           ! column's number of pfts 
    integer , pointer :: pfti(:)            ! column's beginning pft index 
    real(r8), pointer :: pwtcol(:)          ! weight of pft relative to column
    real(r8), pointer :: pwtgcell(:)        ! weight of pft relative to corresponding gridcell
    real(r8), pointer :: forc_lwrad(:)      ! downward infrared (longwave) radiation (W/m**2)
    integer , pointer :: snl(:)             ! number of snow layers
    real(r8), pointer :: htvp(:)            ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: emg(:)             ! ground emissivity
    real(r8), pointer :: cgrnd(:)           ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: dlrad(:)           ! downward longwave radiation blow the canopy [W/m2]
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    integer , pointer :: frac_veg_nosno(:)  ! fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
    real(r8), pointer :: eflx_sh_grnd(:)    ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:)   ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: zi(:,:)            ! interface level below a "z" level (m)
    real(r8), pointer :: dz(:,:)            ! layer depth (m)
    real(r8), pointer :: z(:,:)             ! layer thickness (m)
    real(r8), pointer :: t_soisno(:,:)      ! soil temperature (Kelvin)
! 
! local pointers to  original implicit inout arguments
!
    real(r8), pointer :: t_grnd(:)          ! ground surface temperature [K]
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: eflx_gnet(:)       ! net ground heat flux into the surface (W/m**2)
    real(r8), pointer :: dgnetdT(:)         ! temperature derivative of ground net heat flux
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: j,c,p,g,pi                       !indices
    integer  :: fc                               !lake filtered column indices
    integer  :: jtop(lbc:ubc)                    !top level at each column
    real(r8) :: dtime                            !land model time step (sec)
    real(r8) :: at (lbc:ubc,-nlevsno+1:nlevsoi)  !"a" vector for tridiagonal matrix
    real(r8) :: bt (lbc:ubc,-nlevsno+1:nlevsoi)  !"b" vector for tridiagonal matrix
    real(r8) :: ct (lbc:ubc,-nlevsno+1:nlevsoi)  !"c" vector for tridiagonal matrix
    real(r8) :: rt (lbc:ubc,-nlevsno+1:nlevsoi)  !"r" vector for tridiagonal solution
    real(r8) :: cv (lbc:ubc,-nlevsno+1:nlevsoi)  !heat capacity [J/(m2 K)]
    real(r8) :: tk (lbc:ubc,-nlevsno+1:nlevsoi)  !thermal conductivity [W/(m K)]
    real(r8) :: fn (lbc:ubc,-nlevsno+1:nlevsoi)  !heat diffusion through the layer interface [W/m2]
    real(r8) :: fn1(lbc:ubc,-nlevsno+1:nlevsoi)  !heat diffusion through the layer interface [W/m2]
    real(r8) :: brr(lbc:ubc,-nlevsno+1:nlevsoi)  !temporary
    real(r8) :: dzm                              !used in computing tridiagonal matrix
    real(r8) :: dzp                              !used in computing tridiagonal matrix
    real(r8) :: hs(lbc:ubc)                      !net energy flux into the surface (w/m2)
    real(r8) :: dhsdT(lbc:ubc)                   !d(hs)/dT
    real(r8) :: temp1(lbc:ubc)                   !temporary variable
    real(r8) :: temp2(lbc:ubc)                   !temporary variable
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_lwrad     => clm_a2l%forc_lwrad

    ! Assign local pointers to derived subtypes components (column-level)

    npfts          => clm3%g%l%c%npfts
    pfti           => clm3%g%l%c%pfti
    snl            => clm3%g%l%c%cps%snl
    htvp           => clm3%g%l%c%cps%htvp
    emg            => clm3%g%l%c%cps%emg
    t_grnd         => clm3%g%l%c%ces%t_grnd
    zi             => clm3%g%l%c%cps%zi
    dz             => clm3%g%l%c%cps%dz
    z              => clm3%g%l%c%cps%z
    t_soisno       => clm3%g%l%c%ces%t_soisno

    ! Assign local pointers to derived subtypes components (pft-level)

    pgridcell      => clm3%g%l%c%p%gridcell
    pwtcol         => clm3%g%l%c%p%wtcol
    pwtgcell       => clm3%g%l%c%p%wtgcell  
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    cgrnd          => clm3%g%l%c%p%pef%cgrnd
    dlrad          => clm3%g%l%c%p%pef%dlrad
    sabg           => clm3%g%l%c%p%pef%sabg
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    eflx_gnet      => clm3%g%l%c%p%pef%eflx_gnet
    dgnetdT        => clm3%g%l%c%p%pef%dgnetdT

    ! Get step size

    dtime = get_step_size()

    ! Compute ground surface and soil temperatures

    ! Thermal conductivity and Heat capacity

    call SoilThermProp(lbc, ubc, num_nolakec, filter_nolakec, tk, cv)

    ! Net ground heat flux into the surface and its temperature derivative
    ! Added a pfts loop here to get the average of hs and dhsdT over
    ! all PFTs on the column. Precalculate the terms that do not depend on PFT.

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       temp1(c) =    emg(c) * sb * t_grnd(c)**4
       temp2(c) = 4._r8*emg(c) * sb * t_grnd(c)**3
    end do

    hs(lbc:ubc) = 0._r8
    dhsdT(lbc:ubc) = 0._r8
    do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             if (pwtgcell(p)>0._r8) then
                g = pgridcell(p)
                eflx_gnet(p) = sabg(p) + dlrad(p) + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) &
                               - temp1(c) - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
                dgnetdT(p) = - cgrnd(p) - temp2(c)
                hs(c) = hs(c) + eflx_gnet(p) * pwtcol(p)
                dhsdT(c) = dhsdT(c) + dgnetdT(p) * pwtcol(p)
             end if
          end if
       end do
    end do

    ! Determine heat diffusion through the layer interface and factor used in computing
    ! tridiagonal matrix and set up vector r and vectors a, b, c that define tridiagonal
    ! matrix and solve system

    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then
             if (j == snl(c)+1) then
                fact(c,j) = dtime/cv(c,j) * dz(c,j) / (0.5_r8*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))
                fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
             else if (j <= nlevsoi-1) then
                fact(c,j) = dtime/cv(c,j)
                fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                dzm     = (z(c,j)-z(c,j-1))
             else if (j == nlevsoi) then
                fact(c,j) = dtime/cv(c,j)
                fn(c,j) = 0._r8
             end if
          end if
       enddo
    end do

    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then
             if (j == snl(c)+1) then
                dzp     = z(c,j+1)-z(c,j)
                at(c,j) = 0._r8
                bt(c,j) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                ct(c,j) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                rt(c,j) = t_soisno(c,j) +  fact(c,j)*( hs(c) - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
             else if (j <= nlevsoi-1) then
                dzm     = (z(c,j)-z(c,j-1))
                dzp     = (z(c,j+1)-z(c,j))
                at(c,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                bt(c,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                ct(c,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp
                rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
             else if (j == nlevsoi) then
                dzm     = (z(c,j)-z(c,j-1))
                at(c,j) =   - (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                bt(c,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                ct(c,j) = 0._r8
                rt(c,j) = t_soisno(c,j) - cnfac*fact(c,j)*fn(c,j-1)
             end if
          end if
       enddo
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       jtop(c) = snl(c) + 1
    end do
    call Tridiagonal(lbc, ubc, -nlevsno+1, nlevsoi, jtop, num_nolakec, filter_nolakec, &
                     at, bt, ct, rt, t_soisno(lbc:ubc,-nlevsno+1:nlevsoi))

    ! Melting or Freezing

    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then
             if (j <= nlevsoi-1) then
                fn1(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
             else if (j == nlevsoi) then
                fn1(c,j) = 0._r8
             end if
          end if
       end do
    end do

    do j = -nlevsno+1,nlevsoi
!dir$ prefervector
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then
             if (j == snl(c)+1) then
                brr(c,j) = cnfac*fn(c,j) + (1._r8-cnfac)*fn1(c,j)
             else
                brr(c,j) = cnfac*(fn(c,j)-fn(c,j-1)) + (1._r8-cnfac)*(fn1(c,j)-fn1(c,j-1))
             end if
          end if
       end do
    end do

    call PhaseChange (lbc, ubc, num_nolakec, filter_nolakec, fact, brr, hs, dhsdT, xmf)

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       t_grnd(c) = t_soisno(c,snl(c)+1)
    end do

  end subroutine SoilTemperature

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilThermProp
!
! !INTERFACE:
  subroutine SoilThermProp (lbc, ubc,  num_nolakec, filter_nolakec, tk, cv)
!
! !DESCRIPTION:
! Calculation of thermal conductivities and heat capacities of
! snow/soil layers
! (1) The volumetric heat capacity is calculated as a linear combination
!     in terms of the volumetric fraction of the constituent phases.
!
! (2) The thermal conductivity of soil is computed from the algorithm of
!     Johansen (as reported by Farouki 1981), and of snow is from the
!     formulation used in SNTHERM (Jordan 1991).
! The thermal conductivities at the interfaces between two neighboring
! layers (j, j+1) are derived from an assumption that the flux across
! the interface is equal to that from the node j to the interface and the
! flux from the interface to the node j+1.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clmtype
    use clm_varcon  , only : denh2o, denice, tfrz, tkwat, tkice, tkair, &
                             cpice,  cpliq,  istice, istwet
    use clm_varpar  , only : nlevsno, nlevsoi
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc                       ! column bounds
    integer , intent(in)  :: num_nolakec                    ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(ubc-lbc+1)      ! column filter for non-lake points
    real(r8), intent(out) :: cv(lbc:ubc,-nlevsno+1:nlevsoi) ! heat capacity [J/(m2 K)]
    real(r8), intent(out) :: tk(lbc:ubc,-nlevsno+1:nlevsoi) ! thermal conductivity [W/(m K)]
!
! !CALLED FROM:
! subroutine SoilTemperature in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/13/02, Peter Thornton: migrated to new data structures
! 7/01/03, Mariana Vertenstein: migrated to vector code
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer , pointer :: clandunit(:)     ! column's landunit
    integer , pointer :: ityplun(:)       ! landunit type
    integer , pointer :: snl(:)           ! number of snow layers
    real(r8), pointer :: h2osno(:)        ! snow water (mm H2O)
!
! local pointers to original implicit in arrays
!
    real(r8), pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
    real(r8), pointer :: tksatu(:,:)      ! thermal conductivity, saturated soil [W/m-K]
    real(r8), pointer :: tkmg(:,:)        ! thermal conductivity, soil minerals  [W/m-K]
    real(r8), pointer :: tkdry(:,:)       ! thermal conductivity, dry soil (W/m/Kelvin)
    real(r8), pointer :: csol(:,:)        ! heat capacity, soil solids (J/m**3/Kelvin)
    real(r8), pointer :: dz(:,:)          ! layer depth (m)
    real(r8), pointer :: zi(:,:)          ! interface level below a "z" level (m)
    real(r8), pointer :: z(:,:)           ! layer thickness (m)
    real(r8), pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: l,c,j                     ! indices
    integer  :: fc                        ! lake filtered column indices
    real(r8) :: bw                        ! partial density of water (ice + liquid)
    real(r8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
    real(r8) :: dke                       ! kersten number
    real(r8) :: fl                        ! fraction of liquid or unfrozen water to total water
    real(r8) :: satw                      ! relative total water content of soil.
    real(r8) :: thk(lbc:ubc,-nlevsno+1:nlevsoi) ! thermal conductivity of layer
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (landunit-level)

    ityplun    => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit  => clm3%g%l%c%landunit
    snl        => clm3%g%l%c%cps%snl
    h2osno     => clm3%g%l%c%cws%h2osno
    watsat     => clm3%g%l%c%cps%watsat
    tksatu     => clm3%g%l%c%cps%tksatu
    tkmg       => clm3%g%l%c%cps%tkmg
    tkdry      => clm3%g%l%c%cps%tkdry
    csol       => clm3%g%l%c%cps%csol
    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    t_soisno   => clm3%g%l%c%ces%t_soisno
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice

    ! Thermal conductivity of soil from Farouki (1981)

    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)

          ! Only examine levels from 1->nlevsoi
          if (j >= 1) then
             l = clandunit(c)
             if (ityplun(l) /= istwet .AND. ityplun(l) /= istice) then
                satw = (h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
                satw = min(1._r8, satw)
                if (satw > .1e-6_r8) then
                   fl = h2osoi_liq(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j))
                   if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                      dke = max(0._r8, log10(satw) + 1.0_r8)
                      dksat = tksatu(c,j)
                   else                               ! Frozen soil
                      dke = satw
                      dksat = tkmg(c,j)*0.249_r8**(fl*watsat(c,j))*2.29_r8**watsat(c,j)
                   endif
                   thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
                else
                   thk(c,j) = tkdry(c,j)
                endif
             else
                thk(c,j) = tkwat
                if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
             endif
          endif

          ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
          ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
          if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then
             bw = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/dz(c,j)
             thk(c,j) = tkair + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkair)
          end if

       end do
    end do

    ! Thermal conductivity at the layer interface

    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1 .AND. j <= nlevsoi-1) then
             tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                  /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
          else if (j == nlevsoi) then
             tk(c,j) = 0._r8
          end if
       end do
    end do

    ! Soil heat capacity, from de Vires (1963)

    do j = 1, nlevsoi
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          if (ityplun(l) /= istwet .AND. ityplun(l) /= istice) then
             cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) +   &
               (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
          else
             cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
          endif
          if (j == 1) then
             if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8) then
                cv(c,j) = cv(c,j) + cpice*h2osno(c)
             end if
          end if
       enddo
    end do

    ! Snow heat capacity

    do j = -nlevsno+1,0
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
             cv(c,j) = cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j)
          end if
       end do
    end do

  end subroutine SoilThermProp

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: PhaseChange
!
! !INTERFACE:
  subroutine PhaseChange (lbc, ubc, num_nolakec, filter_nolakec, fact, &
                          brr, hs, dhsdT, xmf)
!
! !DESCRIPTION:
! Calculation of the phase change within snow and soil layers:
! (1) Check the conditions for which the phase change may take place,
!     i.e., the layer temperature is great than the freezing point
!     and the ice mass is not equal to zero (i.e. melting),
!     or the layer temperature is less than the freezing point
!     and the liquid water mass is greater than the allowable supercooled 
!     liquid water calculated from freezing point depression (i.e. freezing).
! (2) Assess the rate of phase change from the energy excess (or deficit)
!     after setting the layer temperature to freezing point.
! (3) Re-adjust the ice and liquid mass, and the layer temperature
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_time_manager, only : get_step_size
    use clm_varcon  , only : tfrz, hfus, grav
    use clm_varpar  , only : nlevsno, nlevsoi
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lbc, ubc                             ! column bounds
    integer , intent(in) :: num_nolakec                          ! number of column non-lake points in column filter
    integer , intent(in) :: filter_nolakec(ubc-lbc+1)            ! column filter for non-lake points
    real(r8), intent(in) :: brr   (lbc:ubc, -nlevsno+1:nlevsoi)  ! temporary
    real(r8), intent(in) :: fact  (lbc:ubc, -nlevsno+1:nlevsoi)  ! temporary
    real(r8), intent(in) :: hs    (lbc:ubc)                      ! net ground heat flux into the surface
    real(r8), intent(in) :: dhsdT (lbc:ubc)                      ! temperature derivative of "hs"
    real(r8), intent(out):: xmf   (lbc:ubc)                      ! total latent heat of phase change
!
! !CALLED FROM:
! subroutine SoilTemperature in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/14/02, Peter Thornton: Migrated to new data structures.
! 7/01/03, Mariana Vertenstein: Migrated to vector code
! 04/25/07 Keith Oleson: CLM3.5 Hydrology
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
    integer , pointer :: snl(:)           !number of snow layers
    real(r8), pointer :: h2osno(:)        !snow water (mm H2O)
!
! local pointers to original implicit inout scalars
!
    real(r8), pointer :: snowdp(:)        !snow height (m)
!
! local pointers to original implicit out scalars
!
    real(r8), pointer :: qflx_snomelt(:)  !snow melt (mm H2O /s)
    real(r8), pointer :: eflx_snomelt(:)  !snow melt heat flux (W/m**2)
!
! local pointers to original implicit in arrays
!
    real(r8), pointer :: h2osoi_liq(:,:)  !liquid water (kg/m2) (new)
    real(r8), pointer :: h2osoi_ice(:,:)  !ice lens (kg/m2) (new)
    real(r8), pointer :: tssbef(:,:)      !temperature at previous time step [K]
    real(r8), pointer :: sucsat(:,:)      !minimum soil suction (mm)
    real(r8), pointer :: watsat(:,:)      !volumetric soil water at saturation (porosity)
    real(r8), pointer :: bsw(:,:)         !Clapp and Hornberger "b"
    real(r8), pointer :: dz(:,:)          !layer thickness (m)
!
! local pointers to original implicit inout arrays
!
    real(r8), pointer :: t_soisno(:,:)    !soil temperature (Kelvin)
!
! local pointers to original implicit out arrays
!
    integer, pointer :: imelt(:,:)        !flag for melting (=1), freezing (=2), Not=0 (new)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: j,c,g                              !do loop index
    integer  :: fc                                 !lake filtered column indices
    real(r8) :: dtime                              !land model time step (sec)
    real(r8) :: heatr                              !energy residual or loss after melting or freezing
    real(r8) :: temp1                              !temporary variables [kg/m2]
    real(r8) :: hm(lbc:ubc,-nlevsno+1:nlevsoi)     !energy residual [W/m2]
    real(r8) :: xm(lbc:ubc,-nlevsno+1:nlevsoi)     !melting or freezing within a time step [kg/m2]
    real(r8) :: wmass0(lbc:ubc,-nlevsno+1:nlevsoi) !initial mass of ice and liquid (kg/m2)
    real(r8) :: wice0 (lbc:ubc,-nlevsno+1:nlevsoi) !initial mass of ice (kg/m2)
    real(r8) :: wliq0 (lbc:ubc,-nlevsno+1:nlevsoi) !initial mass of liquid (kg/m2)
    real(r8) :: supercool(lbc:ubc,nlevsoi)         !supercooled water in soil (kg/m2) 
    real(r8) :: propor                             !proportionality constant (-)
    real(r8) :: tinc                               !t(n+1)-t(n) (K)
    real(r8) :: smp                                !frozen water potential (mm)
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    snl          => clm3%g%l%c%cps%snl
    h2osno       => clm3%g%l%c%cws%h2osno
    snowdp       => clm3%g%l%c%cps%snowdp
    qflx_snomelt => clm3%g%l%c%cwf%qflx_snomelt
    eflx_snomelt => clm3%g%l%c%cef%eflx_snomelt
    h2osoi_liq   => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice   => clm3%g%l%c%cws%h2osoi_ice
    imelt        => clm3%g%l%c%cps%imelt
    t_soisno     => clm3%g%l%c%ces%t_soisno
    tssbef       => clm3%g%l%c%ces%tssbef
    bsw          => clm3%g%l%c%cps%bsw
    sucsat       => clm3%g%l%c%cps%sucsat
    watsat       => clm3%g%l%c%cps%watsat
    dz           => clm3%g%l%c%cps%dz

    ! Get step size

    dtime = get_step_size()

    ! Initialization

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       qflx_snomelt(c) = 0._r8
       xmf(c) = 0._r8
    end do

    do j = -nlevsno+1,nlevsoi       ! all layers
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then

             ! Initialization
             imelt(c,j) = 0
             hm(c,j) = 0._r8
             xm(c,j) = 0._r8
             wice0(c,j) = h2osoi_ice(c,j)
             wliq0(c,j) = h2osoi_liq(c,j)
             wmass0(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

    do j = -nlevsno+1,0             ! snow layers 
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then

             ! Melting identification
             ! If ice exists above melt point, melt some to liquid.
             if (h2osoi_ice(c,j) > 0._r8 .AND. t_soisno(c,j) > tfrz) then
                imelt(c,j) = 1
                t_soisno(c,j) = tfrz
             endif

             ! Freezing identification
             ! If liquid exists below melt point, freeze some to ice.
             if (h2osoi_liq(c,j) > 0._r8 .AND. t_soisno(c,j) < tfrz) then
                imelt(c,j) = 2
                t_soisno(c,j) = tfrz
             endif
          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

    do j = 1,nlevsoi             ! soil layers 
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (h2osoi_ice(c,j) > 0. .AND. t_soisno(c,j) > tfrz) then
             imelt(c,j) = 1
             t_soisno(c,j) = tfrz
          endif

          ! from Zhao (1997) and Koren (1999)
          supercool(c,j) = 0.0_r8
          if(t_soisno(c,j) < tfrz) then
             smp = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
             supercool(c,j) = watsat(c,j)*(smp/sucsat(c,j))**(-1._r8/bsw(c,j))
             supercool(c,j) = supercool(c,j)*dz(c,j)*1000._r8       ! (mm)
          endif

          if (h2osoi_liq(c,j) > supercool(c,j) .AND. t_soisno(c,j) < tfrz) then
             imelt(c,j) = 2
             t_soisno(c,j) = tfrz
          endif

          ! If snow exists, but its thickness is less than the critical value (0.01 m)
          if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. j == 1) then
             if (t_soisno(c,j) > tfrz) then
                imelt(c,j) = 1
                t_soisno(c,j) = tfrz
             endif
          endif
       end do
    enddo

    do j = -nlevsno+1,nlevsoi       ! all layers
!dir$ concurrent
!cdir nodep
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)

          if (j >= snl(c)+1) then

             ! Calculate the energy surplus and loss for melting and freezing
             if (imelt(c,j) > 0) then
                tinc = t_soisno(c,j)-tssbef(c,j)
                if (j > snl(c)+1) then
                   hm(c,j) = brr(c,j) - tinc/fact(c,j)
                else
                   hm(c,j) = hs(c) + dhsdT(c)*tinc + brr(c,j) - tinc/fact(c,j)
                endif
             endif

             ! These two errors were checked carefully (Y. Dai).  They result from the
             ! computed error of "Tridiagonal-Matrix" in subroutine "thermal".
             if (imelt(c,j) == 1 .AND. hm(c,j) < 0._r8) then
                hm(c,j) = 0._r8
                imelt(c,j) = 0
             endif
             if (imelt(c,j) == 2 .AND. hm(c,j) > 0._r8) then
                hm(c,j) = 0._r8
                imelt(c,j) = 0
             endif

             ! The rate of melting and freezing

             if (imelt(c,j) > 0 .and. abs(hm(c,j)) > 0._r8) then
                xm(c,j) = hm(c,j)*dtime/hfus                           ! kg/m2

                ! If snow exists, but its thickness is less than the critical value
                ! (1 cm). Note: more work is needed to determine how to tune the
                ! snow depth for this case
                if (j == 1) then
                   if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. xm(c,j) > 0._r8) then
                      temp1 = h2osno(c)                           ! kg/m2
                      h2osno(c) = max(0._r8,temp1-xm(c,j))
                      propor = h2osno(c)/temp1
                      snowdp(c) = propor * snowdp(c)
                      heatr = hm(c,j) - hfus*(temp1-h2osno(c))/dtime   ! W/m2
                      if (heatr > 0._r8) then
                         xm(c,j) = heatr*dtime/hfus                    ! kg/m2
                         hm(c,j) = heatr                               ! W/m2
                      else
                         xm(c,j) = 0._r8
                         hm(c,j) = 0._r8
                      endif
                      qflx_snomelt(c) = max(0._r8,(temp1-h2osno(c)))/dtime   ! kg/(m2 s)
                      xmf(c) = hfus*qflx_snomelt(c)
                   endif
                endif

                heatr = 0._r8
                if (xm(c,j) > 0._r8) then
                   h2osoi_ice(c,j) = max(0._r8, wice0(c,j)-xm(c,j))
                   heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                else if (xm(c,j) < 0._r8) then
                   if (j <= 0) then
                      h2osoi_ice(c,j) = min(wmass0(c,j), wice0(c,j)-xm(c,j))  ! snow
                   else
                      if (wmass0(c,j) < supercool(c,j)) then
                         h2osoi_ice(c,j) = 0._r8
                      else
                         h2osoi_ice(c,j) = min(wmass0(c,j) - supercool(c,j),wice0(c,j)-xm(c,j))
                      endif
                   endif
                   heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                endif

                h2osoi_liq(c,j) = max(0._r8,wmass0(c,j)-h2osoi_ice(c,j))

                if (abs(heatr) > 0._r8) then
                   if (j > snl(c)+1) then
                      t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr
                   else
                      t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr/(1._r8-fact(c,j)*dhsdT(c))
                   endif
                   if (j <= 0) then    ! snow
                      if (h2osoi_liq(c,j)*h2osoi_ice(c,j)>0._r8) t_soisno(c,j) = tfrz
                   end if
                endif

                xmf(c) = xmf(c) + hfus * (wice0(c,j)-h2osoi_ice(c,j))/dtime

                if (imelt(c,j) == 1 .AND. j < 1) then
                   qflx_snomelt(c) = qflx_snomelt(c) + max(0._r8,(wice0(c,j)-h2osoi_ice(c,j)))/dtime
                endif
             endif

          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

    ! Needed for history file output

!dir$ concurrent
!cdir nodep
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       eflx_snomelt(c) = qflx_snomelt(c) * hfus
    end do

  end subroutine PhaseChange

end module SoilTemperatureMod
