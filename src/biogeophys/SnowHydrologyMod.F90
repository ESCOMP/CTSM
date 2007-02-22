#include <misc.h>
#include <preproc.h>

module SnowHydrologyMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SnowHydrologyMod
!
! !DESCRIPTION:
! Calculate snow hydrology.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar  , only : nlevsno
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SnowWater         ! Change of snow mass and the snow water onto soil
  public :: SnowCompaction    ! Change in snow layer thickness due to compaction
  public :: CombineSnowLayers ! Combine snow layers less than a min thickness
  public :: DivideSnowLayers  ! Subdivide snow layers if they exceed maximum thickness
  public :: BuildSnowFilter   ! Construct snow/no-snow filters
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: Combo            ! Returns the combined variables: dz, t, wliq, wice.
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
! !IROUTINE: SnowWater
!
! !INTERFACE:
  subroutine SnowWater(lbc, ubc, num_snowc, filter_snowc, &
                       num_nosnowc, filter_nosnowc)
!
! !DESCRIPTION:
! Evaluate the change of snow mass and the snow water onto soil.
! Water flow within snow is computed by an explicit and non-physical
! based scheme, which permits a part of liquid water over the holding
! capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
! percolate into the underlying layer.  Except for cases where the
! porosity of one of the two neighboring layers is less than 0.05, zero
! flow is assumed. The water flow out of the bottom of the snow pack will
! participate as the input of the soil water and runoff.  This subroutine
! uses a filter for columns containing snow which must be constructed prior
! to being called.
!
! !USES:
    use clmtype
    use clm_varcon  , only : denh2o, denice, wimp, ssi
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: num_snowc                   ! number of snow points in column filter
    integer, intent(in) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
    integer, intent(in) :: num_nosnowc                 ! number of non-snow points in column filter
    integer, intent(in) :: filter_nosnowc(ubc-lbc+1)   ! column filter for non-snow points
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 15 November 2000: Mariana Vertenstein
! 2/26/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: snl(:)              !number of snow layers
    logical , pointer :: do_capsnow(:)       !true => do snow capping
    real(r8), pointer :: qflx_snomelt(:)     !snow melt (mm H2O /s)
    real(r8), pointer :: qflx_rain_grnd(:)   !rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow(:)    !sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_grnd(:)   !ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_dew_snow(:)    !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)    !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: dz(:,:)             !layer depth (m)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: qflx_top_soil(:)     !net water input into soil from top (mm/s)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: h2osoi_ice(:,:)     !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)     !liquid water (kg/m2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c, j, fc                           !do loop/array indices
    real(r8) :: dtime                              !land model time step (sec)
    real(r8) :: qin(lbc:ubc)                       !water flow into the elmement (mm/s)
    real(r8) :: qout(lbc:ubc)                      !water flow out of the elmement (mm/s)
    real(r8) :: wgdif                              !ice mass after minus sublimation
    real(r8) :: vol_liq(lbc:ubc,-nlevsno+1:0)      !partial volume of liquid water in layer
    real(r8) :: vol_ice(lbc:ubc,-nlevsno+1:0)      !partial volume of ice lens in layer
    real(r8) :: eff_porosity(lbc:ubc,-nlevsno+1:0) !effective porosity = porosity - vol_ice
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    snl            => clm3%g%l%c%cps%snl
    do_capsnow     => clm3%g%l%c%cps%do_capsnow
    qflx_snomelt   => clm3%g%l%c%cwf%qflx_snomelt
    qflx_rain_grnd => clm3%g%l%c%cwf%pwf_a%qflx_rain_grnd
    qflx_sub_snow  => clm3%g%l%c%cwf%pwf_a%qflx_sub_snow
    qflx_evap_grnd => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd
    qflx_dew_snow  => clm3%g%l%c%cwf%pwf_a%qflx_dew_snow
    qflx_dew_grnd  => clm3%g%l%c%cwf%pwf_a%qflx_dew_grnd
    qflx_top_soil  => clm3%g%l%c%cwf%qflx_top_soil
    dz             => clm3%g%l%c%cps%dz
    h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq

    ! Determine model time step

    dtime = get_step_size()

    ! Renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
    ! surface snow layer resulting from sublimation (frost) / evaporation (condense)

!dir$ concurrent
!cdir nodep
    do fc = 1,num_snowc
       c = filter_snowc(fc)
       if (do_capsnow(c)) then
          wgdif = h2osoi_ice(c,snl(c)+1) - qflx_sub_snow(c)*dtime
          h2osoi_ice(c,snl(c)+1) = wgdif
          if (wgdif < 0._r8) then
             h2osoi_ice(c,snl(c)+1) = 0._r8
             h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
          end if
          h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) - qflx_evap_grnd(c) * dtime
       else
          wgdif = h2osoi_ice(c,snl(c)+1) + (qflx_dew_snow(c) - qflx_sub_snow(c)) * dtime
          h2osoi_ice(c,snl(c)+1) = wgdif
          if (wgdif < 0._r8) then
             h2osoi_ice(c,snl(c)+1) = 0._r8
             h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
          end if
          h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) +  &
               (qflx_rain_grnd(c) + qflx_dew_grnd(c) - qflx_evap_grnd(c)) * dtime
       end if
       h2osoi_liq(c,snl(c)+1) = max(0._r8, h2osoi_liq(c,snl(c)+1))
    end do

    ! Porosity and partial volume

    do j = -nlevsno+1, 0
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             vol_ice(c,j) = min(1._r8, h2osoi_ice(c,j)/(dz(c,j)*denice))
             eff_porosity(c,j) = 1._r8 - vol_ice(c,j)
             vol_liq(c,j) = min(eff_porosity(c,j),h2osoi_liq(c,j)/(dz(c,j)*denh2o))
          end if
       end do
    end do

    ! Capillary forces within snow are usually two or more orders of magnitude
    ! less than those of gravity. Only gravity terms are considered.
    ! the genernal expression for water flow is "K * ss**3", however,
    ! no effective parameterization for "K".  Thus, a very simple consideration
    ! (not physically based) is introduced:
    ! when the liquid water of layer exceeds the layer's holding
    ! capacity, the excess meltwater adds to the underlying neighbor layer.

    qin(:) = 0._r8

    do j = -nlevsno+1, 0
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             h2osoi_liq(c,j) = h2osoi_liq(c,j) + qin(c)
             if (j <= -1) then
                ! No runoff over snow surface, just ponding on surface
                if (eff_porosity(c,j) < wimp .OR. eff_porosity(c,j+1) < wimp) then
                   qout(c) = 0._r8
                else
                   qout(c) = max(0._r8,(vol_liq(c,j)-ssi*eff_porosity(c,j))*dz(c,j))
                   qout(c) = min(qout(c),(1._r8-vol_ice(c,j+1)-vol_liq(c,j+1))*dz(c,j+1))
                end if
             else
                qout(c) = max(0._r8,(vol_liq(c,j) - ssi*eff_porosity(c,j))*dz(c,j))
             end if
             qout(c) = qout(c)*1000._r8
             h2osoi_liq(c,j) = h2osoi_liq(c,j) - qout(c)
             qin(c) = qout(c)
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       ! Qout from snow bottom
       qflx_top_soil(c) = qout(c) / dtime
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)
       qflx_top_soil(c) = qflx_rain_grnd(c) + qflx_snomelt(c)
    end do

  end subroutine SnowWater

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowCompaction
!
! !INTERFACE:
  subroutine SnowCompaction(lbc, ubc, num_snowc, filter_snowc)
!
! !DESCRIPTION:
! Determine the change in snow layer thickness due to compaction and
! settling.
! Three metamorphisms of changing snow characteristics are implemented,
! i.e., destructive, overburden, and melt. The treatments of the former
! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution
! due to melt metamorphism is simply taken as a ratio of snow ice
! fraction after the melting versus before the melting.
!
! !USES:
    use clmtype
    use clm_time_manager, only : get_step_size
    use clm_varcon  , only : denice, denh2o, tfrz
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                ! column bounds
    integer, intent(in) :: num_snowc               ! number of column snow points in column filter
    integer, intent(in) :: filter_snowc(ubc-lbc+1) ! column filter for snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    integer,  pointer :: snl(:)             !number of snow layers
!
! local pointers to implicit in arguments
!
    integer,  pointer :: imelt(:,:)        !flag for melting (=1), freezing (=2), Not=0
    real(r8), pointer :: frac_iceold(:,:)  !fraction of ice relative to the tot water
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: dz(:,:)           !layer depth (m)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: j, c, fc                   ! indices
    real(r8):: dtime                      ! land model time step (sec)
    real(r8), parameter :: c2 = 23.e-3_r8    ! [m3/kg]
    real(r8), parameter :: c3 = 2.777e-6_r8  ! [1/s]
    real(r8), parameter :: c4 = 0.04_r8      ! [1/K]
    real(r8), parameter :: c5 = 2.0_r8       !
    real(r8), parameter :: dm = 100.0_r8     ! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
    real(r8), parameter :: eta0 = 9.e+5_r8   ! The Viscosity Coefficient Eta0 [kg-s/m2]
    real(r8) :: burden(lbc:ubc) ! pressure of overlying snow [kg/m2]
    real(r8) :: ddz1   ! Rate of settling of snowpack due to destructive metamorphism.
    real(r8) :: ddz2   ! Rate of compaction of snowpack due to overburden.
    real(r8) :: ddz3   ! Rate of compaction of snowpack due to melt [1/s]
    real(r8) :: dexpf  ! expf=exp(-c4*(273.15-t_soisno)).
    real(r8) :: fi     ! Fraction of ice relative to the total water content at current time step
    real(r8) :: td     ! t_soisno - tfrz [K]
    real(r8) :: pdzdtc ! Nodal rate of change in fractional-thickness due to compaction [fraction/s]
    real(r8) :: void   ! void (1 - vol_ice - vol_liq)
    real(r8) :: wx     ! water mass (ice+liquid) [kg/m2]
    real(r8) :: bi     ! partial density of ice [kg/m3]

!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (column-level)

    snl         => clm3%g%l%c%cps%snl
    dz          => clm3%g%l%c%cps%dz
    imelt       => clm3%g%l%c%cps%imelt
    frac_iceold => clm3%g%l%c%cps%frac_iceold
    t_soisno    => clm3%g%l%c%ces%t_soisno
    h2osoi_ice  => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq  => clm3%g%l%c%cws%h2osoi_liq

    ! Get time step

    dtime = get_step_size()

    ! Begin calculation - note that the following column loops are only invoked if snl(c) < 0

    burden(:) = 0._r8

    do j = -nlevsno+1, 0
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then

             wx = h2osoi_ice(c,j) + h2osoi_liq(c,j)
             void = 1._r8 - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o) / dz(c,j)

             ! Allow compaction only for non-saturated node and higher ice lens node.
             if (void > 0.001_r8 .and. h2osoi_ice(c,j) > .1_r8) then
                bi = h2osoi_ice(c,j) / dz(c,j)
                fi = h2osoi_ice(c,j) / wx
                td = tfrz-t_soisno(c,j)
                dexpf = exp(-c4*td)

                ! Settling as a result of destructive metamorphism

                ddz1 = -c3*dexpf
                if (bi > dm) ddz1 = ddz1*exp(-46.0e-3_r8*(bi-dm))

                ! Liquid water term

                if (h2osoi_liq(c,j) > 0.01_r8*dz(c,j)) ddz1=ddz1*c5

                ! Compaction due to overburden

                ddz2 = -burden(c)*exp(-0.08_r8*td - c2*bi)/eta0

                ! Compaction occurring during melt

                if (imelt(c,j) == 1) then
                   ddz3 = - 1._r8/dtime * max(0._r8,(frac_iceold(c,j) - fi)/frac_iceold(c,j))
                else
                   ddz3 = 0._r8
                end if

                ! Time rate of fractional change in dz (units of s-1)

                pdzdtc = ddz1 + ddz2 + ddz3

                ! The change in dz due to compaction

                dz(c,j) = dz(c,j) * (1._r8+pdzdtc*dtime)
             end if

             ! Pressure of overlying snow

             burden(c) = burden(c) + wx

          end if
       end do
    end do

  end subroutine SnowCompaction

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CombineSnowLayers
!
! !INTERFACE:
  subroutine CombineSnowLayers(lbc, ubc, num_snowc, filter_snowc)
!
! !DESCRIPTION:
! Combine snow layers that are less than a minimum thickness or mass
! If the snow element thickness or mass is less than a prescribed minimum,
! then it is combined with a neighboring element.  The subroutine
! clm\_combo.f90 then executes the combination of mass and energy.
!
! !USES:
    use clmtype
    use clm_varcon, only : istsoil
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: lbc, ubc                    ! column bounds
    integer, intent(inout) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(inout) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer, pointer :: clandunit(:)       !landunit index for each column
    integer, pointer :: ityplun(:)         !landunit type
!
! local pointers to implicit inout arguments
!
    integer , pointer :: snl(:)            !number of snow layers
    real(r8), pointer :: h2osno(:)         !snow water (mm H2O)
    real(r8), pointer :: snowdp(:)         !snow height (m)
    real(r8), pointer :: dz(:,:)           !layer depth (m)
    real(r8), pointer :: zi(:,:)           !interface level below a "z" level (m)
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: z(:,:)            !layer thickness (m)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: c, fc                 ! column indices
    integer :: i,k                   ! loop indices
    integer :: j,l                   ! node indices
    integer :: msn_old(lbc:ubc)      ! number of top snow layer
    integer :: mssi(lbc:ubc)         ! node index
    integer :: neibor                ! adjacent node selected for combination
    real(r8):: zwice(lbc:ubc)        ! total ice mass in snow
    real(r8):: zwliq (lbc:ubc)       ! total liquid water in snow
    real(r8):: dzmin(5)              ! minimum of top snow layer

    data dzmin /0.010_r8, 0.015_r8, 0.025_r8, 0.055_r8, 0.115_r8/
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (landunit-level)

    ityplun    => clm3%g%l%itype

    ! Assign local pointers to derived subtypes (column-level)

    clandunit  => clm3%g%l%c%landunit
    snl        => clm3%g%l%c%cps%snl
    snowdp     => clm3%g%l%c%cps%snowdp
    h2osno     => clm3%g%l%c%cws%h2osno
    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    t_soisno   => clm3%g%l%c%ces%t_soisno
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq

    ! Check the mass of ice lens of snow, when the total is less than a small value,
    ! combine it with the underlying neighbor.

!dir$ concurrent
!cdir nodep
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       msn_old(c) = snl(c)
    end do

    ! The following loop is NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       l = clandunit(c)
       do j = msn_old(c)+1,0
          if (h2osoi_ice(c,j) <= .1_r8) then
             if (ityplun(l) == istsoil) then
                h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
                h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)
             else if (ityplun(l) /= istsoil .and. j /= 0) then
                h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
                h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)
             end if

             ! shift all elements above this down one.
             if (j > snl(c)+1 .and. snl(c) < -1) then
                do i = j, snl(c)+2, -1
                   t_soisno(c,i)   = t_soisno(c,i-1)
                   h2osoi_liq(c,i) = h2osoi_liq(c,i-1)
                   h2osoi_ice(c,i) = h2osoi_ice(c,i-1)
                   dz(c,i)         = dz(c,i-1)
                end do
             end if
             snl(c) = snl(c) + 1
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       h2osno(c) = 0._r8
       snowdp(c) = 0._r8
       zwice(c)  = 0._r8
       zwliq(c)  = 0._r8
    end do

    do j = -nlevsno+1,0
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
             snowdp(c) = snowdp(c) + dz(c,j)
             zwice(c)  = zwice(c) + h2osoi_ice(c,j)
             zwliq(c)  = zwliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Check the snow depth - all snow gone
    ! The liquid water assumes ponding on soil surface.

!dir$ concurrent
!cdir nodep
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       l = clandunit(c)
       if (snowdp(c) < 0.01_r8 .and. snowdp(c) > 0._r8) then
          snl(c) = 0
          h2osno(c) = zwice(c)
          if (h2osno(c) <= 0._r8) snowdp(c) = 0._r8
          if (ityplun(l) == istsoil) h2osoi_liq(c,1) = h2osoi_liq(c,1) + zwliq(c)
       end if
    end do

    ! Check the snow depth - snow layers combined
    ! The following loop IS NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       ! Two or more layers

       if (snl(c) < -1) then

          msn_old(c) = snl(c)
          mssi(c) = 1

          do i = msn_old(c)+1,0
             if (dz(c,i) < dzmin(mssi(c))) then

                if (i == snl(c)+1) then
                   ! If top node is removed, combine with bottom neighbor.
                   neibor = i + 1
                else if (i == 0) then
                   ! If the bottom neighbor is not snow, combine with the top neighbor.
                   neibor = i - 1
                else
                   ! If none of the above special cases apply, combine with the thinnest neighbor
                   neibor = i + 1
                   if ((dz(c,i-1)+dz(c,i)) < (dz(c,i+1)+dz(c,i))) neibor = i-1
                end if

                ! Node l and j are combined and stored as node j.
                if (neibor > i) then
                   j = neibor
                   l = i
                else
                   j = i
                   l = neibor
                end if

                call Combo (dz(c,j), h2osoi_liq(c,j), h2osoi_ice(c,j), &
                   t_soisno(c,j), dz(c,l), h2osoi_liq(c,l), h2osoi_ice(c,l), t_soisno(c,l) )

                ! Now shift all elements above this down one.
                if (j-1 > snl(c)+1) then
                   do k = j-1, snl(c)+2, -1
                      t_soisno(c,k) = t_soisno(c,k-1)
                      h2osoi_ice(c,k) = h2osoi_ice(c,k-1)
                      h2osoi_liq(c,k) = h2osoi_liq(c,k-1)
                      dz(c,k) = dz(c,k-1)
                   end do
                end if

                ! Decrease the number of snow layers
                snl(c) = snl(c) + 1
                if (snl(c) >= -1) EXIT

             else

                ! The layer thickness is greater than the prescribed minimum value
                mssi(c) = mssi(c) + 1

             end if
          end do

       end if

    end do

    ! Reset the node depth and the depth of layer interface

    do j = 0, -nlevsno+1, -1
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c) + 1) then
             z(c,j) = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

  end subroutine CombineSnowLayers

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DivideSnowLayers
!
! !INTERFACE:
  subroutine DivideSnowLayers(lbc, ubc, num_snowc, filter_snowc)
!
! !DESCRIPTION:
! Subdivides snow layers if they exceed their prescribed maximum thickness.
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: lbc, ubc                    ! column bounds
    integer, intent(inout) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(inout) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit inout arguments
!
    integer , pointer :: snl(:)            !number of snow layers
    real(r8), pointer :: dz(:,:)           !layer depth (m)
    real(r8), pointer :: zi(:,:)           !interface level below a "z" level (m)
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: z(:,:)            !layer thickness (m)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: j, c, fc               ! indices
    real(r8) :: drr                    ! thickness of the combined [m]
    integer  :: msno                   ! number of snow layer 1 (top) to msno (bottom)
    real(r8) :: dzsno(lbc:ubc,nlevsno) ! Snow layer thickness [m]
    real(r8) :: swice(lbc:ubc,nlevsno) ! Partial volume of ice [m3/m3]
    real(r8) :: swliq(lbc:ubc,nlevsno) ! Partial volume of liquid water [m3/m3]
    real(r8) :: tsno(lbc:ubc ,nlevsno) ! Nodel temperature [K]
    real(r8) :: zwice                  ! temporary
    real(r8) :: zwliq                  ! temporary
    real(r8) :: propor                 ! temporary
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    snl        => clm3%g%l%c%cps%snl
    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    t_soisno   => clm3%g%l%c%ces%t_soisno
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq

    ! Begin calculation - note that the following column loops are only invoked
    ! for snow-covered columns

    do j = 1,nlevsno
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= abs(snl(c))) then
             dzsno(c,j) = dz(c,j+snl(c))
             swice(c,j) = h2osoi_ice(c,j+snl(c))
             swliq(c,j) = h2osoi_liq(c,j+snl(c))
             tsno(c,j)  = t_soisno(c,j+snl(c))
          end if
       end do
    end do

!dir$ concurrent
!cdir nodep
    do fc = 1, num_snowc
       c = filter_snowc(fc)

       msno = abs(snl(c))

       if (msno == 1) then
          ! Specify a new snow layer
          if (dzsno(c,1) > 0.03_r8) then
             msno = 2
             dzsno(c,1) = dzsno(c,1)/2._r8
             swice(c,1) = swice(c,1)/2._r8
             swliq(c,1) = swliq(c,1)/2._r8
             dzsno(c,2) = dzsno(c,1)
             swice(c,2) = swice(c,1)
             swliq(c,2) = swliq(c,1)
             tsno(c,2)  = tsno(c,1)
          end if
       end if

       if (msno > 1) then
          if (dzsno(c,1) > 0.02_r8) then
             drr = dzsno(c,1) - 0.02_r8
             propor = drr/dzsno(c,1)
             zwice = propor*swice(c,1)
             zwliq = propor*swliq(c,1)
             propor = 0.02_r8/dzsno(c,1)
             swice(c,1) = propor*swice(c,1)
             swliq(c,1) = propor*swliq(c,1)
             dzsno(c,1) = 0.02_r8

             call Combo (dzsno(c,2), swliq(c,2), swice(c,2), tsno(c,2), drr, &
                  zwliq, zwice, tsno(c,1))

             ! Subdivide a new layer
             if (msno <= 2 .and. dzsno(c,2) > 0.07_r8) then
                msno = 3
                dzsno(c,2) = dzsno(c,2)/2._r8
                swice(c,2) = swice(c,2)/2._r8
                swliq(c,2) = swliq(c,2)/2._r8
                dzsno(c,3) = dzsno(c,2)
                swice(c,3) = swice(c,2)
                swliq(c,3) = swliq(c,2)
                tsno(c,3)  = tsno(c,2)
             end if
          end if
       end if

       if (msno > 2) then
          if (dzsno(c,2) > 0.05_r8) then
             drr = dzsno(c,2) - 0.05_r8
             propor = drr/dzsno(c,2)
             zwice = propor*swice(c,2)
             zwliq = propor*swliq(c,2)
             propor = 0.05_r8/dzsno(c,2)
             swice(c,2) = propor*swice(c,2)
             swliq(c,2) = propor*swliq(c,2)
             dzsno(c,2) = 0.05_r8

             call Combo (dzsno(c,3), swliq(c,3), swice(c,3), tsno(c,3), drr, &
                  zwliq, zwice, tsno(c,2))

             ! Subdivided a new layer
             if (msno <= 3 .and. dzsno(c,3) > 0.18_r8) then
                msno =  4
                dzsno(c,3) = dzsno(c,3)/2._r8
                swice(c,3) = swice(c,3)/2._r8
                swliq(c,3) = swliq(c,3)/2._r8
                dzsno(c,4) = dzsno(c,3)
                swice(c,4) = swice(c,3)
                swliq(c,4) = swliq(c,3)
                tsno(c,4)  = tsno(c,3)
             end if
          end if
       end if

       if (msno > 3) then
          if (dzsno(c,3) > 0.11_r8) then
             drr = dzsno(c,3) - 0.11_r8
             propor = drr/dzsno(c,3)
             zwice = propor*swice(c,3)
             zwliq = propor*swliq(c,3)
             propor = 0.11_r8/dzsno(c,3)
             swice(c,3) = propor*swice(c,3)
             swliq(c,3) = propor*swliq(c,3)
             dzsno(c,3) = 0.11_r8

             call Combo (dzsno(c,4), swliq(c,4), swice(c,4), tsno(c,4), drr, &
                  zwliq, zwice, tsno(c,3))

             ! Subdivided a new layer
             if (msno <= 4 .and. dzsno(c,4) > 0.41_r8) then
                msno = 5
                dzsno(c,4) = dzsno(c,4)/2._r8
                swice(c,4) = swice(c,4)/2._r8
                swliq(c,4) = swliq(c,4)/2._r8
                dzsno(c,5) = dzsno(c,4)
                swice(c,5) = swice(c,4)
                swliq(c,5) = swliq(c,4)
                tsno(c,5)  = tsno(c,4)
             end if
          end if
       end if

       if (msno > 4) then
          if (dzsno(c,4) > 0.23_r8) then
             drr = dzsno(c,4) - 0.23_r8
             propor = drr/dzsno(c,4)
             zwice = propor*swice(c,4)
             zwliq = propor*swliq(c,4)
             propor = 0.23_r8/dzsno(c,4)
             swice(c,4) = propor*swice(c,4)
             swliq(c,4) = propor*swliq(c,4)
             dzsno(c,4) = 0.23_r8

             call Combo (dzsno(c,5), swliq(c,5), swice(c,5), tsno(c,5), drr, &
                  zwliq, zwice, tsno(c,4))
          end if
       end if

       snl(c) = -msno

    end do

    do j = -nlevsno+1,0
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             dz(c,j)         = dzsno(c,j-snl(c))
             h2osoi_ice(c,j) = swice(c,j-snl(c))
             h2osoi_liq(c,j) = swliq(c,j-snl(c))
             t_soisno(c,j)   = tsno(c,j-snl(c))
          end if
       end do
    end do

    do j = 0, -nlevsno+1, -1
!dir$ concurrent
!cdir nodep
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

  end subroutine DivideSnowLayers

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Combo
!
! !INTERFACE:
  subroutine Combo(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
!
! !DESCRIPTION:
! Combines two elements and returns the following combined
! variables: dz, t, wliq, wice.
! The combined temperature is based on the equation:
! the sum of the enthalpies of the two elements =
! that of the combined element.
!
! !USES:
    use clm_varcon,  only : cpice, cpliq, tfrz, hfus
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)    :: dz2   ! nodal thickness of 2 elements being combined [m]
    real(r8), intent(in)    :: wliq2 ! liquid water of element 2 [kg/m2]
    real(r8), intent(in)    :: wice2 ! ice of element 2 [kg/m2]
    real(r8), intent(in)    :: t2    ! nodal temperature of element 2 [K]
    real(r8), intent(inout) :: dz    ! nodal thickness of 1 elements being combined [m]
    real(r8), intent(inout) :: wliq  ! liquid water of element 1
    real(r8), intent(inout) :: wice  ! ice of element 1 [kg/m2]
    real(r8), intent(inout) :: t     ! nodel temperature of elment 1 [K]
!
! !CALLED FROM:
! subroutine CombineSnowLayers in this module
! subroutine DivideSnowLayers in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8) :: dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
    real(r8) :: wliqc ! Combined liquid water [kg/m2]
    real(r8) :: wicec ! Combined ice [kg/m2]
    real(r8) :: tc    ! Combined node temperature [K]
    real(r8) :: h     ! enthalpy of element 1 [J/m2]
    real(r8) :: h2    ! enthalpy of element 2 [J/m2]
    real(r8) :: hc    ! temporary
!-----------------------------------------------------------------------

    dzc = dz+dz2
    wicec = (wice+wice2)
    wliqc = (wliq+wliq2)
    h = (cpice*wice+cpliq*wliq) * (t-tfrz)+hfus*wliq
    h2= (cpice*wice2+cpliq*wliq2) * (t2-tfrz)+hfus*wliq2

    hc = h + h2
    if(hc < 0._r8)then
       tc = tfrz + hc/(cpice*wicec + cpliq*wliqc)
    else if (hc.le.hfus*wliqc) then
       tc = tfrz
    else
       tc = tfrz + (hc - hfus*wliqc) / (cpice*wicec + cpliq*wliqc)
    end if

    dz = dzc
    wice = wicec
    wliq = wliqc
    t = tc

  end subroutine Combo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BuildSnowFilter
!
! !INTERFACE:
  subroutine BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec, &
                             num_snowc, filter_snowc, &
                             num_nosnowc, filter_nosnowc)
!
! !DESCRIPTION:
! Constructs snow filter for use in vectorized loops for snow hydrology.
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: lbc, ubc                    ! column bounds
    integer, intent(in)  :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in)  :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(out) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(out) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
    integer, intent(out) :: num_nosnowc                 ! number of column non-snow points in column filter
    integer, intent(out) :: filter_nosnowc(ubc-lbc+1)   ! column filter for non-snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in Hydrology2Mod
! subroutine CombineSnowLayers in this module
!
! !REVISION HISTORY:
! 2003 July 31: Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in arguments
    integer , pointer :: snl(:)                        ! number of snow layers
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer  :: fc, c
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    snl => clm3%g%l%c%cps%snl

    ! Build snow/no-snow filters for other subroutines

    num_snowc = 0
    num_nosnowc = 0
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       if (snl(c) < 0) then
          num_snowc = num_snowc + 1
          filter_snowc(num_snowc) = c
       else
          num_nosnowc = num_nosnowc + 1
          filter_nosnowc(num_nosnowc) = c
       end if
    end do

  end subroutine BuildSnowFilter

end module SnowHydrologyMod
