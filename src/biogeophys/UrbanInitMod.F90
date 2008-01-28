#include <misc.h>
#include <preproc.h>

module UrbanInitMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: UrbanInitMod
! 
! !DESCRIPTION: 
! Initialize urban data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun  
  use shr_sys_mod , only : shr_sys_flush 
  use clm_varctl  , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  save

  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanInitTimeVar   ! Initialize urban time varying variables
  public :: UrbanInitTimeConst ! Initialize urban time constant variables
  public :: UrbanInitAero      ! Calculate urban landunit aerodynamic constants
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInitAero
!
! !INTERFACE:
  subroutine UrbanInitAero( )
!
! !DESCRIPTION: 
! Calculate urban land unit aerodynamic constants using Macdonald (1998) as used in
! Grimmond and Oke (1999)
!
! !USES:
    use clmtype   , only : clm3
    use clm_varcon, only : isturb, vkc
    use decompMod , only : get_proc_bounds
!
! !ARGUMENTS:
    implicit none
!
! local pointers to original implicit in arguments (urban clump)
!
    real(r8), pointer :: ht_roof(:)    ! height of urban roof (m)
    real(r8), pointer :: canyon_hwr(:) ! ratio of building height to street width (-)
    integer , pointer :: ltype(:)      ! landunit type
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: z_0_town(:)   ! urban landunit momentum roughness length (m)
    real(r8), pointer :: z_d_town(:)   ! urban landunit displacement height (m)
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Keith Oleson January 2005
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8), parameter :: alpha = 4.43_r8 ! coefficient used to calculate z_d_town
    real(r8), parameter :: beta = 1.0_r8   ! coefficient used to calculate z_d_town
    real(r8), parameter :: C_d = 1.2_r8    ! drag coefficient as used in Grimmond and Oke (1999)
    real(r8) :: plan_ai                    ! plan area index - ratio building area to plan area (-)
    real(r8) :: frontal_ai                 ! frontal area index of buildings (-)
    real(r8) :: build_lw_ratio             ! building short/long side ratio (-)
    integer  :: l,g                        ! indices
    integer  :: begp, endp                 ! clump beginning and ending pft indices
    integer  :: begc, endc                 ! clump beginning and ending column indices
    integer  :: begl, endl                 ! clump beginning and ending landunit indices
    integer  :: begg, endg                 ! clump beginning and ending gridcell indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (landunit level)

    ltype      => clm3%g%l%itype
    z_0_town   => clm3%g%l%z_0_town
    z_d_town   => clm3%g%l%z_d_town
    ht_roof    => clm3%g%l%ht_roof
    canyon_hwr => clm3%g%l%canyon_hwr

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do l = begl, endl 
      if (ltype(l) == isturb) then 

         ! Calculate plan area index 
         plan_ai = canyon_hwr(l)/(canyon_hwr(l) + 1._r8)

         ! Building shape shortside/longside ratio (e.g. 1 = square )
         ! This assumes the building occupies the entire canyon length
         build_lw_ratio = plan_ai

         ! Calculate frontal area index
         frontal_ai = (1._r8 - plan_ai) * canyon_hwr(l)

         ! Adjust frontal area index for different building configuration
         frontal_ai = frontal_ai * sqrt(1/build_lw_ratio) * sqrt(plan_ai)
         
         ! Calculate displacement height
         
#if (defined VANCOUVER)
         z_d_town(l) = 3.5_r8
#elif (defined MEXICOCITY)
         z_d_town(l) = 10.9_r8
#else
         z_d_town(l) = (1._r8 + alpha**(-plan_ai) * (plan_ai - 1._r8)) * ht_roof(l)
#endif
         
         ! Calculate the roughness length
         
#if (defined VANCOUVER)
         z_0_town(l) = 0.35_r8
#elif (defined MEXICOCITY)
         z_0_town(l) = 2.2_r8
#else
         z_0_town(l) = ht_roof(l) * (1._r8 - z_d_town(l) / ht_roof(l)) * &
                       exp(-1.0_r8 * (0.5_r8 * beta * C_d / vkc**2 * &
                       (1 - z_d_town(l) / ht_roof(l)) * frontal_ai)**(-0.5_r8))
#endif
      end if
   end do

 end subroutine UrbanInitAero

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInitTimeConst
!
! !INTERFACE:
  subroutine UrbanInitTimeConst()
!
! !DESCRIPTION: 
! Initialize urban time-constant variables
!
! !USES:
    use clmtype      , only : clm3
    use clm_varcon   , only : isturb, icol_roof, icol_sunwall, icol_shadewall, &
                              icol_road_perv, icol_road_imperv, spval
    use decompMod    , only : get_proc_bounds, ldecomp
    use UrbanInputMod, only : urbinp
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    integer , pointer :: gdc(:)                 ! grid index for landunit 
    integer , pointer :: coli(:)                ! beginning column index for landunit 
    integer , pointer :: colf(:)                ! ending column index for landunit
    integer , pointer :: ctype(:)               ! column type
    integer , pointer :: ltype(:)               ! landunit type index
    integer , pointer :: lgridcell(:)           ! gridcell of corresponding landunit
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: canyon_hwr(:)          ! urban canyon height to width ratio
    real(r8), pointer :: emg(:)                 ! ground emissivity
    real(r8), pointer :: wtroad_perv(:)         ! weight of pervious column to total road
    real(r8), pointer :: ht_roof(:)             ! height of urban roof (m)
    real(r8), pointer :: wtlunit_roof(:)        ! weight of roof with respect to landunit
    real(r8), pointer :: wind_hgt_canyon(:)     ! height above road at which wind in canyon is to be computed (m)
    real(r8), pointer :: eflx_traffic(:)        ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_traffic_factor(:) ! multiplicative factor for sensible heat flux from urban traffic
    real(r8), pointer :: eflx_wasteheat(:)      ! sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: t_building(:)          ! internal building temperature (K)
    real(r8), pointer :: t_building_max(:)      ! maximum internal building temperature (K)
    real(r8), pointer :: t_building_min(:)      ! minimum internal building temperature (K)
    real(r8), pointer :: tk_wall(:,:)           ! thermal conductivity of urban wall (W/m/K)
    real(r8), pointer :: tk_roof(:,:)           ! thermal conductivity of urban roof (W/m/K)
    real(r8), pointer :: tk_improad(:,:)        ! thermal conductivity of urban impervious road (W/m/K)
    real(r8), pointer :: cv_wall(:,:)           ! thermal conductivity of urban wall (J/m^3/K)
    real(r8), pointer :: cv_roof(:,:)           ! thermal conductivity of urban roof (J/m^3/K)
    real(r8), pointer :: cv_improad(:,:)        ! thermal conductivity of urban impervious road (J/m^3/K)
    real(r8), pointer :: sandfrac_road(:,:)     ! sand fraction of urban road
    real(r8), pointer :: clayfrac_road(:,:)     ! clay fraction of urban road
    real(r8), pointer :: scalez_wall(:)         ! layer thickness discretization of urban wall
    real(r8), pointer :: scalez_roof(:)         ! layer thickness discretization of urban roof
    real(r8), pointer :: thick_wall(:)          ! thickness of urban wall (m)
    real(r8), pointer :: thick_roof(:)          ! thickness of urban roof (m)
!
!EOP
!
! !OTHER LOCAL VARIABLES
    integer  :: nc,fl,ib,l,c,p,g          ! indices
    integer  :: ier                       ! error status
    integer  :: begp, endp                ! clump beginning and ending pft indices
    integer  :: begc, endc                ! clump beginning and ending column indices
    integer  :: begl, endl                ! clump beginning and ending landunit indices
    integer  :: begg, endg                ! clump beginning and ending gridcell indices

    ! Assign local pointers to derived type members (landunit-level)

    ltype               => clm3%g%l%itype
    lgridcell           => clm3%g%l%gridcell
    coli                => clm3%g%l%coli
    colf                => clm3%g%l%colf
    canyon_hwr          => clm3%g%l%canyon_hwr
    wtroad_perv         => clm3%g%l%wtroad_perv 
    ht_roof             => clm3%g%l%ht_roof
    wtlunit_roof        => clm3%g%l%wtlunit_roof
    wind_hgt_canyon     => clm3%g%l%wind_hgt_canyon
    eflx_traffic        => clm3%g%l%lef%eflx_traffic
    eflx_traffic_factor => clm3%g%l%lef%eflx_traffic_factor
    eflx_wasteheat      => clm3%g%l%lef%eflx_wasteheat
    t_building          => clm3%g%l%lps%t_building
    t_building_max      => clm3%g%l%lps%t_building_max
    t_building_min      => clm3%g%l%lps%t_building_min
    canyon_hwr          => clm3%g%l%canyon_hwr
    tk_wall             => clm3%g%l%lps%tk_wall
    tk_roof             => clm3%g%l%lps%tk_roof
    tk_improad          => clm3%g%l%lps%tk_improad
    cv_wall             => clm3%g%l%lps%cv_wall
    cv_roof             => clm3%g%l%lps%cv_roof
    cv_improad          => clm3%g%l%lps%cv_improad
    sandfrac_road       => clm3%g%l%lps%sandfrac_road
    clayfrac_road       => clm3%g%l%lps%clayfrac_road
    scalez_wall         => clm3%g%l%lps%scalez_wall
    scalez_roof         => clm3%g%l%lps%scalez_roof
    thick_wall          => clm3%g%l%lps%thick_wall
    thick_roof          => clm3%g%l%lps%thick_roof

    ! Assign local pointers to derived type members (column-level)

    ctype               => clm3%g%l%c%itype
    emg                 => clm3%g%l%c%cps%emg
    
   ! Initialize time constant urban variables

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do l = begl, endl
       if (ltype(l) == isturb) then
          g =  lgridcell(l)

          canyon_hwr(l)         = urbinp%canyon_hwr(g)
          wtroad_perv(l)        = urbinp%wtroad_perv(g)
          ht_roof(l)            = urbinp%ht_roof(g)
          wtlunit_roof(l)       = urbinp%wtlunit_roof(g)
          wind_hgt_canyon(l)    = urbinp%wind_hgt_canyon(g)
          tk_wall(l,:)          = urbinp%tk_wall(g,:)
          tk_roof(l,:)          = urbinp%tk_roof(g,:)
          tk_improad(l,:)       = urbinp%tk_improad(g,:)
          cv_wall(l,:)          = urbinp%cv_wall(g,:)
          cv_roof(l,:)          = urbinp%cv_roof(g,:)
          cv_improad(l,:)       = urbinp%cv_improad(g,:)
          sandfrac_road(l,:)    = urbinp%sandfrac_road(g,:)
          clayfrac_road(l,:)    = urbinp%clayfrac_road(g,:)
          scalez_wall(l)        = urbinp%scalez_wall(g)
          scalez_roof(l)        = urbinp%scalez_roof(g)
          thick_wall(l)         = urbinp%thick_wall(g)
          thick_roof(l)         = urbinp%thick_roof(g)

          do c = coli(l),colf(l)
             if (ctype(c) == icol_roof       ) emg(c) = urbinp%em_roof(g)
             if (ctype(c) == icol_sunwall    ) emg(c) = urbinp%em_wall(g)
             if (ctype(c) == icol_shadewall  ) emg(c) = urbinp%em_wall(g)
             if (ctype(c) == icol_road_imperv) emg(c) = urbinp%em_improad(g)
             if (ctype(c) == icol_road_perv  ) emg(c) = urbinp%em_perroad(g)
          end do

          ! Inferred from Sailor and Lu 2004
          eflx_traffic_factor(l) = 0.0_r8
!KO          eflx_traffic_factor(l) = 3.6_r8 * (canyon_hwr(l)-0.5_r8) + 1.0_r8

#if (defined VANCOUVER || defined MEXICOCITY || defined GRANDVIEW)
          ! Minimal heating or air conditioning
          t_building_max(l) = 340.00_r8
          t_building_min(l) = 230.00_r8
#else
          ! Arbitrary comfort values
          t_building_max(l) = 297.60_r8  !~75F
          t_building_min(l) = 291.16_r8  !~65F
#endif
       else
          eflx_traffic(l) = spval
          eflx_traffic_factor(l) = spval
          eflx_wasteheat(l) = spval
          t_building(l)     = spval
          t_building_max(l) = spval
          t_building_min(l) = spval
       end if
    end do

  end subroutine UrbanInitTimeConst

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInitTimeVar
!
! !INTERFACE:
  subroutine UrbanInitTimeVar( )
!
! !DESCRIPTION: 
! Initialize urban time-varying variables
!
! !USES:
    use clmtype   , only : clm3
    use clm_varcon, only : isturb, spval, icol_road_perv
    use decompMod , only : get_proc_bounds
!
! !ARGUMENTS:
    implicit none
!
! local pointers to original implicit in arguments (urban clump)
!
    integer , pointer :: ltype(:)      ! landunit type
    integer , pointer :: lgridcell(:)  ! gridcell of corresponding landunit
    integer , pointer :: clandunit(:)  ! landunit index of corresponding column
    integer , pointer :: ctype(:)      ! column type
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: taf(:)                ! urban canopy air temperature (K)
    real(r8), pointer :: qaf(:)                ! urban canopy air specific humidity (kg/kg)
    real(r8), pointer :: eflx_building_heat(:) ! heat flux from urban building interior to walls, roof (W/m**2)
    real(r8), pointer :: fcov(:)               ! fractional area with water table at surface
    real(r8), pointer :: qcharge(:)            !aquifer recharge rate (mm/s)
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Keith Oleson February 2005
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: l,g,c         ! indices
    integer :: begp, endp    ! clump beginning and ending pft indices
    integer :: begc, endc    ! clump beginning and ending column indices
    integer :: begl, endl    ! clump beginning and ending landunit indices
    integer :: begg, endg    ! clump beginning and ending gridcell indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (landunit level)

    taf                => clm3%g%l%taf
    qaf                => clm3%g%l%qaf
    ltype              => clm3%g%l%itype
    lgridcell          => clm3%g%l%gridcell
    clandunit          => clm3%g%l%c%landunit
    eflx_building_heat => clm3%g%l%c%cef%eflx_building_heat
    fcov               => clm3%g%l%c%cws%fcov
    qcharge            => clm3%g%l%c%cws%qcharge
    ctype              => clm3%g%l%c%itype

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do l = begl, endl 
       g = lgridcell(l)
       if (ltype(l) == isturb) then 
#if (defined VANCOUVER)
          taf(l) = 297.56_r8
          qaf(l) = 0.0111_r8
#elif (defined MEXICOCITY)
          taf(l) = 289.46_r8
          qaf(l) = 0.00248_r8
#elif (defined GRANDVIEW)
          taf(l) = 291.56_r8
          qaf(l) = 0.002_r8
#else
          taf(l) = 283._r8
          ! Arbitrary set since forc_q is not yet available
          qaf(l) = 1.e-4_r8
#endif
       end if
    end do
    do c = begc, endc 
       l = clandunit(c)
       if (ltype(l) == isturb) then 
          eflx_building_heat(c) = 0._r8
          !
          ! Set hydrology variables for urban to spvalue -- as only valid for pervious road
          !
          if (ctype(c) /= icol_road_perv  )then
             fcov(c)    = spval
             qcharge(c) = spval
          end if
       else
          eflx_building_heat(c) = spval
       end if
    end do
    
  end subroutine UrbanInitTimeVar
  
end module UrbanInitMod
