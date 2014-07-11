module SurfaceResistanceMod
#include "shr_assert.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SurfaceResistanceMod
!
! !DESCRIPTION:
! Module holding routines for calculation of surface resistances of the different tracers
! transported with BeTR. The surface here refers to water and soil, not including canopy
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_TKFRZ
  use clm_varctl  , only: iulog
  
  implicit none
  save
  private
  integer :: soil_stress_method   !choose the method for soil resistance calculation
  
  integer, parameter :: leepielke_1992 = 0 !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: calc_soilevap_stress
  public :: do_soilevap_beta
  public :: init_soil_stress
  !
  ! !REVISION HISTORY:
  ! 6/25/2013 Created by Jinyun Tang
  !
  !EOP
  !-----------------------------------------------------------------------
  
contains

   subroutine init_soil_stress()
   !
   !DESCRIPTIONS
   ! initialize method for soil stress calculation
   implicit none
   
   
   soil_stress_method = leepielke_1992
   
   end subroutine init_soil_stress
   
  !------------------------------------------------------------------------------   
   
   subroutine calc_soilevap_stress(bounds, num_nolakec, filter_nolakec, col, lun, cps, cws)
   !
   ! DESCRIPTIONS
   ! compute the stress factor for soil evaporation calculation
   !
   use shr_kind_mod   , only : r8 => shr_kind_r8     
   use shr_const_mod      , only : SHR_CONST_PI  
   use decompMod,      only : bounds_type
   use clmtype,        only : column_type, column_wstate_type, column_pstate_type
   use clmtype,        only : landunit_type
   use abortutils    , only : endrun      
   implicit none
   type(bounds_type), intent(in) :: bounds    ! bounds   
   integer, intent(in) :: num_nolakec
   integer, intent(in) :: filter_nolakec(:)
   type(column_type),        intent(in) :: col
   type(landunit_type),      intent(in) :: lun
   type(column_wstate_type), intent(inout) :: cws
   type(column_pstate_type), intent(in) :: cps

   character(len=32) :: subname = 'calc_soilevap_stress'  ! subroutine name
   associate(                &
   soilbeta                  =>    cws%soilbeta                  & ! Output: [real(r8) (:)] factor that reduces ground evaporation
   )
   
   
   !select the right method and do the calculation
   select case (soil_stress_method)
   
   case (leepielke_1992)
     call calc_beta_leepielke1992(bounds, num_nolakec, filter_nolakec,lun, col, cws, cps, soilbeta(bounds%begc:bounds%endc))
   case default
      call endrun(subname // ':: a soilevap stress function must be specified!')     
   end select
   
   end associate
   end subroutine calc_soilevap_stress
   
   
  !------------------------------------------------------------------------------   
   subroutine calc_beta_leepielke1992(bounds, num_nolakec, filter_nolakec, lun, col, cws, cps, soilbeta)
   !
   !
   ! DESCRIPTION
   ! compute the lee-pielke beta factor to scal actual soil evaporation from potential evaporation
   !
   ! USES
   use shr_kind_mod   , only : r8 => shr_kind_r8     
   use shr_const_mod      , only : SHR_CONST_PI
   use shr_log_mod    , only : errMsg => shr_log_errMsg   
   use decompMod,      only : bounds_type
   use clmtype,        only : column_type, column_wstate_type, column_pstate_type
   use clmtype,        only : landunit_type
   use shr_infnan_mod, only : &
       nan => shr_infnan_nan, assignment(=)
   use clm_varcon         , only : denh2o, denice
   use landunit_varcon    , only : istice, istice_mec, istwet, istsoil, istcrop
   use column_varcon,       only :  icol_roof, icol_sunwall, icol_shadewall, &
         icol_road_imperv, icol_road_perv
   implicit none
   type(bounds_type), intent(in) :: bounds    ! bounds   
   integer, intent(in) :: num_nolakec
   integer, intent(in) :: filter_nolakec(:)
   type(column_type),        intent(in) :: col
   type(landunit_type),      intent(in) :: lun
   type(column_wstate_type), intent(in) :: cws
   type(column_pstate_type), intent(in) :: cps
   real(r8), intent(inout) :: soilbeta(bounds%begc:bounds%endc)
   
   !local variables
   real(r8) :: fac, fac_fc, wx      !temporary variables
   integer  :: c, l, fc     !indices
   

   SHR_ASSERT_ALL((ubound(soilbeta)    == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
  
   associate(                                                     &
    ityplun                   =>    lun%itype                    , & ! Input:  [integer (:)] landunit type                             
    clandunit                 =>   col%landunit                  , & ! Input:  [integer (:)] column's landunit index                   
    ctype                     =>    col%itype                    , & ! Input:  [integer (:)] column type                               
    h2osoi_ice                =>    cws%h2osoi_ice               , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2)                       
    h2osoi_liq                =>    cws%h2osoi_liq               , & ! Input:  [real(r8) (:,:)] liquid water (kg/m2)                   
    watsat                    =>    cps%watsat                   , & ! Input:  [real(r8) (:,:)] volumetric soil water at saturation (porosity)
    watfc                     =>    cps%watfc                    , & ! Input:  [real(r8) (:,:)] volumetric soil water at field capacity
    dz                        =>    cps%dz                       , & ! Input:  [real(r8) (:,:)] layer depth (m)                        
    frac_sno                  =>    cps%frac_sno                 , & ! Input:  [real(r8) (:)] fraction of ground covered by snow (0 to 1)
    frac_h2osfc               =>    cps%frac_h2osfc                & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   )
   
   do fc = 1,num_nolakec
      c = filter_nolakec(fc)
      l = clandunit(c)   
      if (ityplun(l)/=istwet .AND. ityplun(l)/=istice  &
                              .AND. ityplun(l)/=istice_mec) then
         if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
            wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)
            fac  = min(1._r8, wx/watsat(c,1))
            fac  = max( fac, 0.01_r8 )
            !! Lee and Pielke 1992 beta, added by K.Sakaguchi
            if (wx < watfc(c,1) ) then  !when water content of ths top layer is less than that at F.C.
               fac_fc  = min(1._r8, wx/watfc(c,1))  !eqn5.66 but divided by theta at field capacity
               fac_fc  = max( fac_fc, 0.01_r8 )
               ! modify soil beta by snow cover. soilbeta for snow surface is one
               soilbeta(c) = (1._r8-frac_sno(c)-frac_h2osfc(c)) &
                     *0.25_r8*(1._r8 - cos(SHR_CONST_PI*fac_fc))**2._r8 &
                              + frac_sno(c)+ frac_h2osfc(c)
            else   !when water content of ths top layer is more than that at F.C.
               soilbeta(c) = 1._r8
            end if             
         else if (ctype(c) == icol_road_perv) then
            soilbeta(c) = 0._r8
         else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall) then
            soilbeta(c) = 0._r8          
         else if (ctype(c) == icol_roof .or. ctype(c) == icol_road_imperv) then
            soilbeta(c) = 0._r8
         endif   
      else
         soilbeta(c) =   1._r8
      endif
   enddo   
   end associate
   end subroutine calc_beta_leepielke1992
   
  !------------------------------------------------------------------------------   
   function do_soilevap_beta()result(lres)
   !
   !DESCRIPTION
   ! return true if the moisture stress for soil evaporation is computed as beta factor
   ! otherwise false
   implicit none
   logical :: lres
   
   
   if(soil_stress_method==leepielke_1992)then
      lres=.true.
   else
      lres=.false.
   endif   
   return
   end function do_soilevap_beta
end module SurfaceResistanceMod
