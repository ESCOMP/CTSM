module SoilMoistStressMod
#include "shr_assert.h"
  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates soil moisture stress for plant gpp and transpiration
  !
  ! After discussion with other developers, I have now removed all functions that
  ! return array, and decalared all variables that will be modified as intent(inout).
  ! The initialization will be done whenever the variable is initialized. This avoids
  ! code crash when initialization is not done appropriately, and make the code safer
  ! during the long-term maintenance
  !
  ! Created by Jinyun Tang, Feb., 2014 
  implicit none
  
  save
  private
  
  integer ::   root_moist_stress_method
  integer, parameter :: moist_stress_clm_default  = 0  !default method for calculating root moisture stress
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: calc_root_moist_stress
  public :: calc_effective_soilporosity
  public :: calc_effective_snowporosity
  public :: calc_volumetric_h2oliq
  public :: set_perchroot_opt
  public :: init_root_moist_stress
  
  ! !PUBLIC DATA MEMBERS:
  logical,  private :: perchroot     = .false.  ! true => btran is based only on unfrozen soil levels
  logical,  private :: perchroot_alt = .false.  ! true => btran is based on active layer (defined over two years); 
  
  contains

  subroutine init_root_moist_stress()
  !
  !DESCRIPTION
  !specify the method to compute root soil moisture stress
  !
  implicit none
  
  root_moist_stress_method = moist_stress_clm_default   
  end subroutine init_root_moist_stress

!--------------------------------------------------------------------------------
  subroutine set_perchroot_opt(perchroot_global, perchroot_alt_global)
  !
  !DESCRIPTIONS
  !set up local perchroot logical switches, in the future, this wil be
  !read in as namelist
  implicit none
  logical, intent(in) :: perchroot_global
  logical, intent(in) :: perchroot_alt_global
  
  perchroot = perchroot_global
  perchroot_alt = perchroot_alt_global
  
  end subroutine set_perchroot_opt
  
!--------------------------------------------------------------------------------
  
  subroutine calc_effective_soilporosity(bounds, ubj, numf, filter, watsat, &
     h2osoi_ice, dz, denice, eff_por)
  !
  ! !DESCRIPTIONS
  ! compute the effective soil porosity
  !
  ! !USES
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  
  implicit none
  type(bounds_type), intent(in) :: bounds                         !bounds
  integer , intent(in) :: ubj                                     !lbinning level indices
  integer , intent(in) :: numf                                    !filter dimension
  integer , intent(in) :: filter(:)                               !filter  
  real(r8), intent(in) :: watsat( bounds%begc: , 1: )             !soil porosity
  real(r8), intent(in) :: h2osoi_ice( bounds%begc: , 1: )         !ice water content, kg H2o/m2
  real(r8), intent(in) :: dz(bounds%begc: , 1: )                  !layer thickness, m
  real(r8), intent(in) :: denice                                  !ice density, kg/m3

  real(r8), intent(inout) :: eff_por( bounds%begc: ,1: )  ! effective porosity

  !return variable

  !local variables
  integer :: c, j, fc                                             !indices
  real(r8):: vol_ice    !volumetric ice
  !obtain the array size
  
  ! Enforce expected array sizes
  SHR_ASSERT_ALL((ubound(watsat)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(h2osoi_ice)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(dz)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(eff_por)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))

  
  !main calculation loop
  !it assumes the soil layers start from 1
  do j = 1, ubj
    do fc = 1, numf
      c = filter(fc)
      !compute the volumetric ice content
      vol_ice=min(watsat(c,j), h2osoi_ice(c,j)/(denice*dz(c,j)))
      
      !compute the maximum soil space to fill liquid water and air
      eff_por(c,j) = watsat(c,j) - vol_ice
    enddo
  enddo
  end subroutine calc_effective_soilporosity
!--------------------------------------------------------------------------------

  subroutine calc_effective_snowporosity(bounds, lbj, jtop, numf, filter, &
     h2osoi_ice, dz, denice, eff_por)
  !
  ! !DESCRIPTIONS
  ! compute the effective porosity snow
  !
  ! !USES
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg    
  implicit none
  type(bounds_type), intent(in) :: bounds                         !bounds
  integer , intent(in) :: lbj                                     !ubing level indices
  integer , intent(in) :: jtop( bounds%begc: )                    !top level for each column [col]    
  integer , intent(in) :: numf                                    !filter dimension
  integer , intent(in) :: filter(:)                               !filter  
  real(r8), intent(in) :: h2osoi_ice( bounds%begc: , lbj: )       !ice water content, kg H2o/m2
  real(r8), intent(in) :: dz(bounds%begc: , lbj: )                !layer thickness, m
  real(r8), intent(in) :: denice                                  !ice density, kg/m3

  real(r8), intent(inout) :: eff_por( bounds%begc: ,lbj: )        !returning effective porosity
  
  
  !local variables
  integer :: c, j, fc                                             !indices
  integer :: ubj
  real(r8) :: vol_ice     !volumetric ice
  
  ubj = 0
 ! Enforce expected array sizes
  SHR_ASSERT_ALL((ubound(jtop)    == (/bounds%endc/)), errMsg(__FILE__, __LINE__)) 
  SHR_ASSERT_ALL((ubound(h2osoi_ice)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(dz)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(eff_por) == (/bounds%endc,0/)), errMsg(__FILE__, __LINE__))
  
  !main calculation loop

  !it assumes snow layer ends at 0
  do j = lbj,0
    do fc = 1, numf
      c = filter(fc)
      if(j>=jtop(c))then
        !compute the volumetric ice content       
        vol_ice=min(1._r8, h2osoi_ice(c,j)/(denice*dz(c,j)))
         
        !compute the maximum snow void space to fill liquid water and air         
        eff_por(c,j) = 1._r8 - vol_ice
      endif
    enddo
  enddo

  end subroutine calc_effective_snowporosity
  
!--------------------------------------------------------------------------------
  subroutine calc_volumetric_h2oliq(bounds, jtop, lbj, ubj, numf, filter,&
     eff_porosity, h2osoi_liq, dz, denh2o, vol_liq)
  !
  ! !DESCRIPTIONS
  ! compute the volumetric liquid water content
  !
  !
  ! !USES
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg  
  implicit none
  type(bounds_type), intent(in) :: bounds                       !bounds
  integer , intent(in)    :: jtop( bounds%begc: )               ! top level for each column [col]  
  integer , intent(in) :: lbj, ubj                              !lbinning and ubing level indices
  integer , intent(in) :: numf                                  !filter dimension
  integer , intent(in) :: filter(:)                             !filter    
  real(r8), intent(in) :: eff_porosity(bounds%begc: , lbj: )    !effective soil porosity
  real(r8), intent(in) :: h2osoi_liq(bounds%begc: , lbj: )      !liquid water content [kg H2o/m2]
  real(r8), intent(in) :: dz(bounds%begc: , lbj: )              !layer thickness [m]
  real(r8), intent(in) :: denh2o                                !water density [kg/m3]

  real(r8), intent(inout) :: vol_liq(bounds%begc: , lbj: )         !volumetric liquid water content  
  !local variables    
  integer :: c, j, fc                                             !indices  
  
  ! Enforce expected array sizes  
  SHR_ASSERT_ALL((ubound(jtop)    == (/bounds%endc/)), errMsg(__FILE__, __LINE__)) 
  SHR_ASSERT_ALL((ubound(h2osoi_liq)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(eff_porosity)  == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(dz)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))  
  SHR_ASSERT_ALL((ubound(vol_liq)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))  


  !main calculation loop
  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(j>=jtop(c))then
         !volume of liquid is no greater than effective void space
         vol_liq(c,j) = min(eff_porosity(c,j), h2osoi_liq(c,j)/(dz(c,j)*denh2o))
      endif
    enddo
  enddo

  end subroutine calc_volumetric_h2oliq
!--------------------------------------------------------------------------------

  subroutine normalize_unfrozen_rootfr(bounds, ubj, fn, filterp, pcolumn, cps_data,&
     pps_data, ces_data, rootfr_unf)
  !
  ! !DESCRIPTIONS
  ! normalize root fraction for total unfrozen depth 
  !
  !
  ! !USES
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment 
  use SimpleMathMod  , only : array_normalization
  use clmtype        , only : column_pstate_type
  use clmtype        , only : column_estate_type
  use clmtype        , only : pft_pstate_type
  implicit none
  type(bounds_type), intent(in) :: bounds                    !bounds
  type(column_pstate_type), intent(in) :: cps_data        !derived type to pass input column physical state variables
  type(pft_pstate_type),    intent(in) :: pps_data        !pft physical state
  type(column_estate_type), intent(in) :: ces_data        !column energy state
  integer , intent(in) :: ubj                             !ubinning level indices
  integer , intent(in) :: fn                              !filter dimension
  integer , intent(in) :: filterp(:)                      !filter
  integer , intent(in) :: pcolumn(bounds%begp: )          !column info related to pft

  real(r8), intent(inout):: rootfr_unf(bounds%begp:bounds%endp, 1:ubj)  !normalized root fraction in unfrozen layers


  !local variables
  !real(r8) :: rootsum(bounds%begp:bounds%endp)  
  integer :: p, c, j, f                                   !indices  
 
  ! Enforce expected array sizes   
  SHR_ASSERT_ALL((ubound(pcolumn)     == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
  
  associate(                                               &
    rootfr               => pps_data%rootfr              , & !fraction of roots in each soil layer
    t_soisno             => ces_data%t_soisno            , & !soil temperature    [K]
    altmax_lastyear_indx => cps_data%altmax_lastyear_indx, & !prior year maximum annual depth of thaw
    altmax_indx          => cps_data%altmax_indx           & !maximum annual depth of thaw
  )
  !main calculation loop  
  !Initialize rootfr_unf to zero.
  !I found it necessary to ensure the pgi compiler not
  !to complain with float point exception. However, it raises a question how
  !to make sure those values that are initialized with nan or spval are not reset
  !to zero within similar coding style. Jinyun Tang, May 23, 2014.
  
  ! Define rootfraction for unfrozen soil only
  if (perchroot .or. perchroot_alt) then
    if (perchroot_alt) then
      ! use total active layer (defined ass max thaw depth for current and prior year)
      do j = 1, ubj
        do f = 1, fn
          p = filterp(f)
          c = pcolumn(p)
          
          if ( j <= max(altmax_lastyear_indx(c), altmax_indx(c), 1) )then
            rootfr_unf(p,j) = rootfr(p,j)
          else
            rootfr_unf(p,j) = 0._r8
          end if
        end do
      end do
    else
      ! use instantaneous temperature
      do j = 1, ubj
        do f = 1, fn
          p = filterp(f)
          c = pcolumn(p)
               
          if (t_soisno(c,j) >= tfrz) then
            rootfr_unf(p,j) = rootfr(p,j)
          else
            rootfr_unf(p,j) = 0._r8
          end if
        end do
      end do
         
    end if ! perchroot_alt          
  end if ! perchroot
  
  !normalize the root fraction for each pft
  call array_normalization(bounds%begp, bounds%endp, 1, ubj, &
            fn, filterp, rootfr_unf(bounds%begp:bounds%endp, 1:ubj))


  end associate        
  
  end subroutine normalize_unfrozen_rootfr
  
!--------------------------------------------------------------------------------
  subroutine calc_root_moist_stress_clm45default(bounds, &
     nlevgrnd, fn, filterp, rootfr_unf, pft_data,  &
     pftcon_data,  cws_data, cps_data, ces_data, &
     pps_data)
  !
  ! DESCRIPTIONS
  ! compute the root water stress using the default clm45 approach
  
  !
  ! USES
  use shr_kind_mod   , only : r8 => shr_kind_r8  
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment
  use clmtype        , only : pft_epc_type
  use clmtype        , only : pft_type
  use clmtype        , only : pft_pstate_type
  use clmtype        , only : column_pstate_type
  use clmtype        , only : column_wstate_type
  use clmtype        , only : column_estate_type
  use SoiWatRetCurveParMod,only : soil_suction
  implicit none
  type(bounds_type),        intent(in) :: bounds                        !bounds
  real(r8)                , intent(in) :: rootfr_unf(bounds%begp: , 1: )!                  !
  type(pft_type)          , intent(in) :: pft_data                      !data structure to pass input pft relevant variables
  type(pft_epc_type)      , intent(in) :: pftcon_data                   !data structure to pass input pft ecophysiological constants structur
  type(column_pstate_type), intent(in) :: cps_data                      !data structure to pass input column physical state variables
  type(column_wstate_type), intent(in) :: cws_data
  type(column_estate_type), intent(in) :: ces_data  
  integer,                  intent(in) :: nlevgrnd                      !number of vertical layers
  integer,                  intent(in) :: fn                            !number of filters
  integer,                  intent(in) :: filterp(:)                    !filter array          
  type(pft_pstate_type)   , intent(inout) :: pps_data                      !data structure to pass input from pft pstate variables

  !local parameter
  real(r8), parameter :: btran0 = 0.0_r8  ! initial value
  !local variables
  real(r8) :: smp_node, s_node  !temporary variables
  real(r8) :: smp_node_lf       !temporary variable
  integer :: p, f, j, c, l      !indices

  ! Enforce expected array sizes   
  SHR_ASSERT_ALL((ubound(rootfr_unf)     == (/bounds%endp, nlevgrnd/)), errMsg(__FILE__, __LINE__))  
  associate(                                                &
   ivt                   => pft_data%itype                , & !vegetation type
   pcolumn               => pft_data%column               , & !column indices of the given pfts
   plandunit             => pft_data%landunit             , & !land unit of the given pfts
   rootfr                => pps_data%rootfr               , & !fraction of roots in each soil layer
   btran                 => pps_data%btran                , & !integrated soil water stress
   btran2                => pps_data%btran2               , & !integrated soil water stress square
   rootr                 => pps_data%rootr                , & !active root fraction for uptaking water 
   rresis                => pps_data%rresis               , & !root soil water stress ([0, 1])   
   smpso                 => pftcon_data%smpso             , & !soil water pontential at full stomatal opening (mm)
   smpsc                 => pftcon_data%smpsc             , & !soil water pontential at full stomatal closure (mm)
   h2osoi_liqvol         => cws_data%h2osoi_liqvol        , & !liquid volumetric moisture, will be used for BeTR
   h2osoi_vol            => cws_data%h2osoi_vol           , & !ltotal volumetric moisture
   eff_porosity          => cps_data%eff_porosity         , & !effective soil porosity, will be used for BeTR
   watsat                => cps_data%watsat               , & !soil porosity
   sucsat                => cps_data%sucsat               , & !saturated soil water potential (mm)
   bsw                   => cps_data%bsw                  , & !clapp-hornberg shape parameter
   t_soisno              => ces_data%t_soisno               & !soil temperature
  ) 
  do j = 1,nlevgrnd
    do f = 1, fn
      p = filterp(f)
      c = pcolumn(p)
      l = plandunit(p)

      ! Root resistance factors
      ! rootr effectively defines the active root fraction in each layer      
      if (h2osoi_liqvol(c,j) .le. 0._r8 .or. t_soisno(c,j) .le. tfrz-2._r8) then
        rootr(p,j) = 0._r8
      else
        s_node = max(h2osoi_liqvol(c,j)/eff_porosity(c,j),0.01_r8)

        !smp_node = max(smpsc(ivt(p)), -sucsat(c,j)*s_node**(-bsw(c,j)))
        call soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_node)
        smp_node = max(smpsc(ivt(p)), smp_node)

        rresis(p,j) = min( (eff_porosity(c,j)/watsat(c,j))* &
                          (smp_node - smpsc(ivt(p))) / (smpso(ivt(p)) - smpsc(ivt(p))), 1._r8)
                          
                          
        if (.not. (perchroot .or. perchroot_alt) ) then
           rootr(p,j) = rootfr(p,j)*rresis(p,j)
        else
           rootr(p,j) = rootfr_unf(p,j)*rresis(p,j)
        end if
        !it is possible to further separate out a btran function, but I will leave it for the moment, jyt
        btran(p)    = btran(p) + max(rootr(p,j),0._r8)
        !smp_node_lf = max(smpsc(ivt(p)), -sucsat(c,j)*(h2osoi_vol(c,j)/watsat(c,j))**(-bsw(c,j)))
        s_node = h2osoi_vol(c,j)/watsat(c,j)
        call soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_node_lf)
        !smp_node_lf =  -sucsat(c,j)*(h2osoi_vol(c,j)/watsat(c,j))**(-bsw(c,j))
        smp_node_lf = max(smpsc(ivt(p)), smp_node_lf) 
        btran2(p)   = btran2(p) +rootfr(p,j)*min((smp_node_lf - smpsc(ivt(p))) / (smpso(ivt(p)) - smpsc(ivt(p))), 1._r8)
      endif 
    end do
  end do
  
  ! Normalize root resistances to get layer contribution to ET
   do j = 1,nlevgrnd
      do f = 1, fn
         p = filterp(f)
         if (btran(p) > btran0) then
           rootr(p,j) = rootr(p,j)/btran(p)
         else
           rootr(p,j) = 0._r8
         end if
      end do
   end do  
  end associate
  end subroutine calc_root_moist_stress_clm45default

!--------------------------------------------------------------------------------
  subroutine calc_root_moist_stress(bounds, nlevgrnd, fn, filterp, &
     pft_data, pftcon_data, cws_data, cps_data, ces_data, pps_data)
  !
  ! DESCRIPTIONS
  ! compute the root water stress using different approaches
  
  !
  ! USES
  use shr_kind_mod   , only : r8 => shr_kind_r8  
  use decompMod      , only : bounds_type
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clmtype        , only : pft_type
  use clmtype        , only : pft_epc_type
  use clmtype        , only : pft_pstate_type
  use clmtype        , only : column_wstate_type
  use clmtype        , only : column_estate_type  
  use clmtype        , only : column_pstate_type
  use clm_varcon     , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment 
  use abortutils    , only : endrun       
  implicit none
  type(bounds_type),        intent(in)    :: bounds                       !bounds
  integer,                  intent(in)    :: nlevgrnd
  integer,                  intent(in)    :: fn
  integer,                  intent(in)    :: filterp(:)
  type(pft_type),           intent(in)    :: pft_data                   !data structure to pass input pft relevant variables 
  type(pft_epc_type),       intent(in)    :: pftcon_data
  type(column_wstate_type), intent(in)    :: cws_data
  type(column_pstate_type), intent(in)    :: cps_data
  type(column_estate_type), intent(in)    :: ces_data
  type(pft_pstate_type)   , intent(inout) :: pps_data    
  
  real(r8) :: smp_node, s_node  !temporary variables
  real(r8) :: rootfr_unf(bounds%begp:bounds%endp,1:nlevgrnd) ! Rootfraction defined for unfrozen layers only.
  
  integer :: p, f, j, c, l      !indices

   character(len=32) :: subname = 'calc_root_moist_stress'  ! subroutine name

  !define normalized rootfraction for unfrozen soil
  !define normalized rootfraction for unfrozen soil
  rootfr_unf(bounds%begp:bounds%endp,1:nlevgrnd) = 0._r8
  
  call normalize_unfrozen_rootfr(bounds, &
      ubj = nlevgrnd,                                          &
      fn = fn,                                                 &
      filterp = filterp,                                       &
      pcolumn = pft_data%column(bounds%begp:bounds%endp),      &
      cps_data =  cps_data,                                    &
      pps_data =  pps_data,                                    &
      ces_data = ces_data,                                     &
      rootfr_unf=rootfr_unf(bounds%begp:bounds%endp,1:nlevgrnd))

  !suppose h2osoi_liq, eff_porosity are already computed somewhere else
  
  select case (root_moist_stress_method)
  !add other methods later
  case (moist_stress_clm_default)
  
    call calc_root_moist_stress_clm45default(bounds,                    &
       nlevgrnd = nlevgrnd,                                             &
       fn = fn,                                                         &
       filterp = filterp,                                               &
       rootfr_unf=rootfr_unf(bounds%begp:bounds%endp,1:nlevgrnd),       &
       pft_data = pft_data,                                             &
       pftcon_data = pftcon_data,                                       &
       cws_data = cws_data,                                             &
       cps_data = cps_data,                                             &
       ces_data = ces_data,                                             &       
       pps_data = pps_data                                              )
  case default
      call endrun(subname // ':: a root moisture stress function must be specified!')     
  end select

  end subroutine calc_root_moist_stress  
end module SoilMoistStressMod
