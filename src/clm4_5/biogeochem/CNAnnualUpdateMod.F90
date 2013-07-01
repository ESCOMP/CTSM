module CNAnnualUpdateMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNAnnualUpdateMod
!
! !DESCRIPTION:
! Module for updating annual summation variables
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CNAnnualUpdate
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNAnnualUpdate
!
! !INTERFACE:
subroutine CNAnnualUpdate(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                          num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update annual summation variables
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size, get_days_per_year
   use clm_varcon      , only: secspday
   use pft2colMod      , only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: lbp, ubp        ! pft bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 10/1/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
   integer :: c,p          ! indices
   integer :: fp,fc        ! lake filter indices
   real(r8):: dt           ! radiation time step (seconds)
!EOP
!-----------------------------------------------------------------------
   associate(& 
   annsum_counter            => cps%annsum_counter         , & ! InOut:  [real(r8) (:)]  seconds since last annual accumulator turnover
   tempsum_potential_gpp     => pepv%tempsum_potential_gpp , & ! InOut:  [real(r8) (:)]  temporary annual sum of potential GPP   
   annsum_potential_gpp      => pepv%annsum_potential_gpp  , & ! InOut:  [real(r8) (:)]  annual sum of potential GPP             
   tempmax_retransn          => pepv%tempmax_retransn      , & ! InOut:  [real(r8) (:)]  temporary annual max of retranslocated N pool (gN/m2)
   annmax_retransn           => pepv%annmax_retransn       , & ! InOut:  [real(r8) (:)]  annual max of retranslocated N pool (gN/m2)
   tempavg_t2m               => pepv%tempavg_t2m           , & ! InOut:  [real(r8) (:)]  temporary average 2m air temperature (K)
   annavg_t2m                => pepv%annavg_t2m            , & ! InOut:  [real(r8) (:)]  annual average 2m air temperature (K)   
   tempsum_npp               => pepv%tempsum_npp           , & ! InOut:  [real(r8) (:)]  temporary sum NPP (gC/m2/yr)            
   annsum_npp                => pepv%annsum_npp            , & ! InOut:  [real(r8) (:)]  annual sum NPP (gC/m2/yr)               
   cannsum_npp               => cps%cannsum_npp            , & ! InOut:  [real(r8) (:)]  column annual sum NPP (gC/m2/yr)        
   cannavg_t2m               => cps%cannavg_t2m            , & ! InOut:  [real(r8) (:)] annual average of 2m air temperature, averaged from pft-level (K)
#if (defined CNDV)
   tempsum_litfall           => pepv%tempsum_litfall       , & ! InOut:  [real(r8) (:)]  temporary sum litfall (gC/m2/yr)        
   annsum_litfall            => pepv%annsum_litfall        , & ! InOut:  [real(r8) (:)]  annual sum litfall (gC/m2/yr)           
#endif
   pcolumn                   => pft%column                 & ! Input:  [integer (:)]  index into column level                  
   )

   ! set time steps
   dt = real( get_step_size(), r8 )

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      annsum_counter(c) = annsum_counter(c) + dt
   end do

   if (num_soilc .gt. 0) then

   if (annsum_counter(filter_soilc(1)) >= get_days_per_year() * secspday) then
      ! pft loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
            ! update annual plant ndemand accumulator
            annsum_potential_gpp(p)  = tempsum_potential_gpp(p)
            tempsum_potential_gpp(p) = 0._r8

            ! update annual total N retranslocation accumulator
            annmax_retransn(p)  = tempmax_retransn(p)
            tempmax_retransn(p) = 0._r8

            ! update annual average 2m air temperature accumulator
            annavg_t2m(p)  = tempavg_t2m(p)
            tempavg_t2m(p) = 0._r8

            ! update annual NPP accumulator, convert to annual total
            annsum_npp(p) = tempsum_npp(p) * dt
            tempsum_npp(p) = 0._r8

#if (defined CNDV)
            ! update annual litfall accumulator, convert to annual total
            annsum_litfall(p) = tempsum_litfall(p) * dt
            tempsum_litfall(p) = 0._r8
#endif
      end do

      ! use p2c routine to get selected column-average pft-level fluxes and states
      call p2c(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, annsum_npp, cannsum_npp)
      call p2c(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, annavg_t2m, cannavg_t2m)
   end if

   end if

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      if (annsum_counter(c) >= get_days_per_year() * secspday) annsum_counter(c) = 0._r8
   end do

    end associate 
 end subroutine CNAnnualUpdate
!-----------------------------------------------------------------------
#endif

end module CNAnnualUpdateMod
