#include <misc.h>
#include <preproc.h>

module mkarbinitMod

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkarbinit
!
! !INTERFACE:
  subroutine mkarbinit()
!
! !DESCRIPTION:
! Initializes the following time varying variables:
! water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
! snow       : snowdp, snowage, snl, dz, z, zi
! temperature: t_soisno, t_veg, t_grnd
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use shr_const_mod, only : SHR_CONST_TKFRZ
    use clmtype
    use clm_varpar   , only : nlevsoi, nlevsno, nlevlak
    use clm_varcon   , only : bdsno, istice, istwet, istsoil, isturb, &
                              denice, denh2o, spval, sb, icol_road_perv, &
                              icol_road_imperv, icol_roof, icol_sunwall, &
                              icol_shadewall
    use clm_varctl   , only : iulog, pertlim
    use spmdMod      , only : masterproc
    use decompMod    , only : get_proc_bounds
    use shr_sys_mod  , only : shr_sys_flush
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)        ! column index associated with each pft
    integer , pointer :: ctype(:)          ! column type
    integer , pointer :: clandunit(:)      ! landunit index associated with each column
    integer , pointer :: ltype(:)          ! landunit type
    logical , pointer :: lakpoi(:)         ! true => landunit is a lake point
    real(r8), pointer :: dz(:,:)           ! layer thickness depth (m)
    real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity) (nlevsoi)
    real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(r8), pointer :: bsw2(:,:)         ! Clapp and Hornberger "b" for CN code
    real(r8), pointer :: psisat(:,:)       ! soil water potential at saturation for CN code (MPa)
    real(r8), pointer :: vwcsat(:,:)       ! volumetric water content at saturation for CN code (m3/m3)
    real(r8), pointer :: zi(:,:)           ! interface level below a "z" level (m)
    real(r8), pointer :: wa(:)             ! water in the unconfined aquifer (mm)
    real(r8), pointer :: wt(:)             ! total water storage (unsaturated soil water + groundwater) (mm)
    real(r8), pointer :: zwt(:)            ! water table depth (m)
!
! local pointers to implicit out arguments
!
    integer , pointer :: snl(:)            ! number of snow layers
    real(r8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
    real(r8), pointer :: t_lake(:,:)       ! lake temperature (Kelvin)  (1:nlevlak)
    real(r8), pointer :: t_grnd(:)         ! ground temperature (Kelvin)
    real(r8), pointer :: t_veg(:)          ! vegetation temperature (Kelvin)
    real(r8), pointer :: t_ref2m(:)        ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: h2osoi_vol(:,:)   ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: h2ocan_col(:)     ! canopy water (mm H2O) (column-level)
    real(r8), pointer :: h2ocan_pft(:)     ! canopy water (mm H2O) (pft-level)
    real(r8), pointer :: h2osno(:)         ! snow water (mm H2O)
    real(r8), pointer :: snowdp(:)         ! snow height (m)
    real(r8), pointer :: snowage(:)        ! non dimensional snow age [-] (new)
    real(r8), pointer :: eflx_lwrad_out(:) ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: soilpsi(:,:)      ! soil water potential in each soil layer (MPa)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer :: j,l,c,p      ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8):: vwc,psi      ! for calculating soilpsi
    real(r8):: pertval      ! for calculating temperature perturbation
!-----------------------------------------------------------------------

    if ( masterproc )then
        write(iulog,*) 'Setting initial data to non-spun up values'
        if ( pertlim /= 0.0_r8 ) &
        write(iulog,*) 'Applying perturbation to initial surface temperature'
    end if

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype
    lakpoi     => clm3%g%l%lakpoi

    ! Assign local pointers to derived subtypes components (column-level)

    ctype      => clm3%g%l%c%itype
    clandunit  => clm3%g%l%c%landunit
    snl        => clm3%g%l%c%cps%snl
    dz         => clm3%g%l%c%cps%dz
    watsat     => clm3%g%l%c%cps%watsat
    bsw2       => clm3%g%l%c%cps%bsw2
    vwcsat     => clm3%g%l%c%cps%vwcsat
    psisat     => clm3%g%l%c%cps%psisat
    soilpsi    => clm3%g%l%c%cps%soilpsi
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol
    h2ocan_col => clm3%g%l%c%cws%pws_a%h2ocan
    snowage    => clm3%g%l%c%cps%snowage
    snowdp     => clm3%g%l%c%cps%snowdp
    h2osno     => clm3%g%l%c%cws%h2osno
    t_soisno   => clm3%g%l%c%ces%t_soisno
    t_lake     => clm3%g%l%c%ces%t_lake
    t_grnd     => clm3%g%l%c%ces%t_grnd
    zi         => clm3%g%l%c%cps%zi
    wa         => clm3%g%l%c%cws%wa
    wt         => clm3%g%l%c%cws%wt
    zwt        => clm3%g%l%c%cws%zwt

    ! Assign local pointers to derived subtypes components (pft-level)

    pcolumn        => clm3%g%l%c%p%column
    h2ocan_pft     => clm3%g%l%c%p%pws%h2ocan
    t_veg          => clm3%g%l%c%p%pes%t_veg
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
    eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out  

    ! Determine subgrid bounds on this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! NOTE: h2ocan, h2osno, snowdp and snowage has valid values everywhere
    ! canopy water (pft level)

    do p = begp, endp
       h2ocan_pft(p) = 0._r8
       
       ! added for canopy water mass balance under dynamic pft weights
       !clm3%g%l%c%p%pps%tlai(p) = 0._r8
       !clm3%g%l%c%p%pps%tsai(p) = 0._r8
       !clm3%g%l%c%p%pps%elai(p) = 0._r8
       !clm3%g%l%c%p%pps%esai(p) = 0._r8
       !clm3%g%l%c%p%pps%htop(p) = 0._r8
       !clm3%g%l%c%p%pps%hbot(p) = 0._r8
       !clm3%g%l%c%p%pps%frac_veg_nosno_alb(p) = 0._r8
    end do

!dir$ concurrent
!cdir nodep
    do c = begc,endc

       ! canopy water (column level)

       h2ocan_col(c) = 0._r8

       ! snow water

       l = clandunit(c)
       if (ltype(l) == istice) then
          h2osno(c) = 1000._r8
       else
          h2osno(c) = 0._r8
       endif

       ! snow depth

       snowdp(c)  = h2osno(c) / bdsno

       ! snow age

       snowage(c) = 0._r8

    end do

    ! Set snow layer number, depth and thickiness

    call snowdp2lev(begc, endc)

    ! Set snow/soil temperature, note:
    ! t_soisno only has valid values over non-lake
    ! t_lake   only has valid values over lake
    ! t_grnd has valid values over all land
    ! t_veg  has valid values over all land

!dir$ concurrent
!cdir nodep
    do c = begc,endc

       t_soisno(c,-nlevsno+1:nlevsoi) = spval
       t_lake(c,1:nlevlak) = spval

       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake
          t_soisno(c,-nlevsno+1:0) = spval
          if (snl(c) < 0) then    !snow layer temperatures
             do j = snl(c)+1, 0
                t_soisno(c,j) = 250._r8
             enddo
          endif
          if (ltype(l) == istice) then
             do j = 1, nlevsoi
                t_soisno(c,j) = 250._r8
             end do
          else if (ltype(l) == istwet) then
             do j = 1, nlevsoi
                t_soisno(c,j) = 277._r8
             end do
          else
#if (defined VANCOUVER)
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
             ! Set road top layer to initial air temperature and interpolate other
             ! layers down to 20C in bottom layer
               do j = 1, nlevsoi
                  t_soisno(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevsoi-1)) 
               end do
             ! Set wall and roof layers to initial air temperature
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
               do j = 1, nlevsoi
                  t_soisno(c,j) = 297.56
               end do
             else
               do j = 1, nlevsoi
                  t_soisno(c,j) = 283._r8
               end do
             end if
#elif (defined MEXICOCITY)
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
             ! Set road top layer to initial air temperature and interpolate other
             ! layers down to 22C in bottom layer
               do j = 1, nlevsoi
                  t_soisno(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevsoi-1)) 
               end do
             ! Set wall and roof layers to initial air temperature
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_roof) then
               do j = 1, nlevsoi
                  t_soisno(c,j) = 289.46
               end do
             else
               do j = 1, nlevsoi
                  t_soisno(c,j) = 283._r8
               end do
             end if
#elif (defined GRANDVIEW)
             if (ctype(c) == icol_road_perv .or. ctype(c) == icol_road_imperv) then 
             ! Set road top layer to 18.5C and interpolate other
             ! layers down to 16.5C in bottom layer
               do j = 1, nlevsoi
                  t_soisno(c,j) = 291.66 - (j-1) * ((291.66-289.66)/(nlevsoi-1)) 
               end do
             ! Set wall top layer to 18.35C and interpolate other
             ! layers down to 17.0C in bottom layer
             else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall) then
               do j = 1, nlevsoi
                  t_soisno(c,j) = 291.51 - (j-1) * ((291.51-290.16)/(nlevsoi-1)) 
               end do
             else if (ctype(c) == icol_roof) then
               do j = 1, nlevsoi
                  t_soisno(c,j) = 291.56_r8
               end do
             end if
#else
             do j = 1, nlevsoi
                t_soisno(c,j) = 283._r8
             end do
#endif
          endif
          t_grnd(c) = t_soisno(c,snl(c)+1)
       else                     !lake
          t_lake(c,1:nlevlak) = 277._r8
          t_grnd(c) = t_lake(c,1)
       endif
       if ( pertlim /= 0.0_r8 )then
          if (.not. lakpoi(l)) then  !not lake
             do j = 1, nlevsoi
                call random_number (pertval)
                pertval       = 2._r8*pertlim*(0.5_r8 - pertval)
                t_soisno(c,j) = t_soisno(c,j)*(1._r8 + pertval)
             end do
             t_grnd(c) = t_soisno(c,snl(c)+1)
          else                       !lake
             do j = 1, nlevlak
                call random_number (pertval)
                pertval     = 2._r8*pertlim*(0.5_r8 - pertval)
                t_lake(c,j) = t_lake(c,j)*(1._r8 + pertval)
             end do
             t_grnd(c) = t_lake(c,1)
          endif
       end if

    end do

!dir$ concurrent
!cdir nodep
    do p = begp, endp
       c = pcolumn(p)
#if (defined VANCOUVER)
       t_veg(p) = 297.56
       t_ref2m(p) = 297.56
#elif (defined MEXICOCITY)
       t_veg(p) = 289.46
       t_ref2m(p) = 289.46
#elif (defined GRANDVIEW)
       t_veg(p) = 291.56
       t_ref2m(p) = 291.56
#else
       t_veg(p) = 283._r8
       t_ref2m(p) = 283._r8
#endif
       eflx_lwrad_out(p) = sb * (t_grnd(c))**4
    end do

    ! Set snow/soil ice and liquid mass

    ! volumetric water is set first and liquid content and ice lens are obtained
    ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
    ! and urban pervious road (other urban columns have zero soil water)

    h2osoi_vol(begc:endc,         1:nlevsoi) = spval
    h2osoi_liq(begc:endc,-nlevsno+1:nlevsoi) = spval
    h2osoi_ice(begc:endc,-nlevsno+1:nlevsoi) = spval

    wa(begc:endc)  = 5000._r8
    wt(begc:endc)  = 5000._r8
    zwt(begc:endc) = 0._r8

!dir$ concurrent
!cdir nodep
    do c = begc,endc
       l = clandunit(c)
       if (.not. lakpoi(l)) then  !not lake
          wa(c)  = 4800._r8
          wt(c)  = wa(c)
          zwt(c) = (25._r8 + zi(c,nlevsoi)) - wa(c)/0.2_r8 /1000._r8  ! One meter below soil column
       end if
    end do

    do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          l = clandunit(c)
          if (.not. lakpoi(l)) then  !not lake

             ! volumetric water
             if (ltype(l) == istsoil) then
                h2osoi_vol(c,j) = 0.4_r8
!KO                h2osoi_vol(c,j) = 0.05_r8
             else if (ltype(l) == isturb) then 
                if (ctype(c) == icol_road_perv) then
                   h2osoi_vol(c,j) = 0.3_r8
                else
                   h2osoi_vol(c,j) = 0.0_r8
                end if
             else
                h2osoi_vol(c,j) = 1.0_r8
             endif
             h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
           
             ! soil layers
             if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                h2osoi_ice(c,j)  = dz(c,j)*denice*h2osoi_vol(c,j)
                h2osoi_liq(c,j) = 0._r8
             else
                h2osoi_ice(c,j) = 0._r8
                h2osoi_liq(c,j) = dz(c,j)*denh2o*h2osoi_vol(c,j)
             endif

#if (defined CN)
             ! soil water potential (added 10/21/03, PET)
             ! required for CN code
             if (ltype(l) == istsoil) then
                if (h2osoi_liq(c,j) > 0._r8) then
                   vwc = h2osoi_liq(c,j)/(dz(c,j)*denh2o)
                   psi = psisat(c,j) * (vwc/vwcsat(c,j))**bsw2(c,j)
                   soilpsi(c,j) = max(psi, -15.0_r8)
                   soilpsi(c,j) = min(soilpsi(c,j),0.0_r8)
                end if
             end if
#endif
          end if

       end do

    end do

    ! Set snow

    do j = -nlevsno+1, 0
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          l = clandunit(c)
          if (.not. lakpoi(l)) then  !not lake
             if (j > snl(c)) then
                h2osoi_ice(c,j) = dz(c,j)*250._r8
                h2osoi_liq(c,j) = 0._r8
             end if
          end if
       end do
    end do

  end subroutine mkarbinit

end module mkarbinitMod
