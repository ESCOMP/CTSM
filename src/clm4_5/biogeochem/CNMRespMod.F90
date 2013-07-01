module CNMRespMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNMRespMod
!
! !DESCRIPTION:
! Module holding maintenance respiration routines for coupled carbon
! nitrogen code.
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use clm_varpar   , only: nlevgrnd
   use shr_const_mod, only: SHR_CONST_TKFRZ
   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public :: CNMResp
!
! !REVISION HISTORY:
! 8/14/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: Migrated all subroutines to vector data structures.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNMResp
!
! !INTERFACE:
subroutine CNMResp(lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
!
! !USES:
   use clmtype
   use pftvarcon    , only : npcropmin
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc                    ! column-index bounds
   integer, intent(in) :: num_soilc                 ! number of soil points in column filter
   integer, intent(in) :: filter_soilc(:)   ! column filter for soil points
   integer, intent(in) :: num_soilp                 ! number of soil points in pft filter
   integer, intent(in) :: filter_soilp(:)   ! pft filter for soil points
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 8/14/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
   real(r8), pointer :: grainn(:)     ! (kgN/m2) grain N
   integer :: c,p,j          ! indices
   integer :: fp             ! soil filter pft index
   integer :: fc             ! soil filter column index
   real(r8):: mr             ! maintenance respiration (gC/m2/s)
   real(r8):: br             ! base rate (gC/gN/s)
   real(r8):: q10            ! temperature dependence
   real(r8):: tc             ! temperature correction, 2m air temp (unitless)
   real(r8):: tcsoi(lbc:ubc,nlevgrnd) ! temperature correction by soil layer (unitless)
!EOP
!-----------------------------------------------------------------------

   associate(&    
   t_soisno       =>    ces%t_soisno         , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   t_ref2m        =>    pes%t_ref2m          , & ! Input:  [real(r8) (:)]  2 m height surface air temperature (Kelvin)       
   leafn          =>    pns%leafn            , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
   frootn         =>    pns%frootn           , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
   livestemn      =>    pns%livestemn        , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
   livecrootn     =>    pns%livecrootn       , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   rootfr         =>    pps%rootfr           , & ! Input:  [real(r8) (:,:)]  fraction of roots in each soil layer  (nlevgrnd)
   leaf_mr        =>    pcf%leaf_mr          , & ! InOut:  [real(r8) (:)]                                                    
   froot_mr       =>    pcf%froot_mr         , & ! InOut:  [real(r8) (:)]                                                    
   livestem_mr    =>    pcf%livestem_mr      , & ! InOut:  [real(r8) (:)]                                                    
   livecroot_mr   =>    pcf%livecroot_mr     , & ! InOut:  [real(r8) (:)]                                                    
   grain_mr       =>    pcf%grain_mr         , & ! InOut:  [real(r8) (:)]                                                    
   lmrsun         =>    pcf%lmrsun           , & ! InOut:  [real(r8) (:)]  sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
   lmrsha         =>    pcf%lmrsha           , & ! InOut:  [real(r8) (:)]  shaded leaf maintenance respiration rate (umol CO2/m**2/s)
   laisun         =>    pps%laisun           , & ! InOut:  [real(r8) (:)]  sunlit projected leaf area index                  
   laisha         =>    pps%laisha           , & ! InOut:  [real(r8) (:)]  shaded projected leaf area index                  
   frac_veg_nosno =>    pps%frac_veg_nosno   , & ! InOut:  [integer (:)]  fraction of vegetation not covered by snow (0 OR 1) [-]
   ivt            =>    pft%itype            , & ! Input:  [integer (:)]  pft vegetation type                                
   pcolumn        =>    pft%column           , & ! Input:  [integer (:)]  index into column level quantities                 
   plandunit      =>    pft%landunit         , & ! Input:  [integer (:)]  index into landunit level quantities               
   clandunit      =>    col%landunit         , & ! Input:  [integer (:)]  index into landunit level quantities               
   itypelun       =>    lun%itype            , & ! Input:  [integer (:)]  landunit type                                      
   woody          =>    pftcon%woody           & ! Input:  [real(r8) (:)]  binary flag for woody lifeform (1=woody, 0=not woody)
   )
   grainn         => pns%grainn

   ! base rate for maintenance respiration is from:
   ! M. Ryan, 1991. Effects of climate change on plant respiration.
   ! Ecological Applications, 1(2), 157-167.
   ! Original expression is br = 0.0106 molC/(molN h)
   ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
   br = 2.525e-6_r8
   ! Peter Thornton: 3/13/09 
   ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
   ! to improve seasonal cycle of atmospheric CO2 concentration in global
   ! simulatoins
   q10 = 1.5_r8

   ! column loop to calculate temperature factors in each soil layer
   do j=1,nlevgrnd
      do fc = 1, num_soilc
         c = filter_soilc(fc)

         ! calculate temperature corrections for each soil layer, for use in
         ! estimating fine root maintenance respiration with depth

         tcsoi(c,j) = q10**((t_soisno(c,j)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
      end do
   end do

   ! pft loop for leaves and live wood
   do fp = 1, num_soilp
      p = filter_soilp(fp)

      ! calculate maintenance respiration fluxes in
      ! gC/m2/s for each of the live plant tissues.
      ! Leaf and live wood MR

      tc = q10**((t_ref2m(p)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
      if (frac_veg_nosno(p) == 1) then
         leaf_mr(p) = lmrsun(p) * laisun(p) * 12.011e-6_r8 + &
                      lmrsha(p) * laisha(p) * 12.011e-6_r8
      else
         leaf_mr(p) = 0._r8
      end if

      if (woody(ivt(p)) == 1) then
         livestem_mr(p) = livestemn(p)*br*tc
         livecroot_mr(p) = livecrootn(p)*br*tc
      else if (ivt(p) >= npcropmin) then
         livestem_mr(p) = livestemn(p)*br*tc
         grain_mr(p) = grainn(p)*br*tc
      end if
   end do

   ! soil and pft loop for fine root
   do j = 1,nlevgrnd
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = pcolumn(p)

         ! Fine root MR
         ! rootfr(j) sums to 1.0 over all soil layers, and
         ! describes the fraction of root mass that is in each
         ! layer.  This is used with the layer temperature correction
         ! to estimate the total fine root maintenance respiration as a
         ! function of temperature and N content.

         froot_mr(p) = froot_mr(p) + frootn(p)*br*tcsoi(c,j)*rootfr(p,j)
      end do
   end do

    end associate 
 end subroutine CNMResp

#endif

end module CNMRespMod
