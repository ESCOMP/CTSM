 module EDAccumulateFluxesMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This routine accumulates NPP, GPP and respiration of each cohort over the course of each 24 hour period. 
  ! The fluxes are stored per cohort, and the npp_clm (etc) fluxes are calcualted in EDPhotosynthesis
  ! This routine cannot be in EDPhotosynthesis because EDPhotosynthesis is a loop and therefore would
  ! erroneously add these things up multiple times. 
  ! Rosie Fisher. March 2014. 
  
  
  ! !USES:
  use EDtypesMod, only : patch, cohort, gridCellEdState
  use clmtype,    only : pcf, pft

  implicit none
  private
  save
  
  public :: AccumulateFluxes_ED

  type(cohort), pointer  :: currentCohort ! current cohort
  type(patch) , pointer  :: currentPatch ! current patch

  contains

  !------------------------------------------------------------------------------
  subroutine AccumulateFluxes_ED(p)
      
   use shr_kind_mod,   only : r8 => shr_kind_r8
   use EDClmtype,      only : EDpft
   
   implicit none

   integer, intent(in) :: p  !patch/'p'

   integer :: iv !leaf layer
   integer :: g  !gridcell
  
   associate(& 
   fpsn                => pcf%fpsn        , & ! Output: [real(r8) (:)]   photosynthesis (umol CO2 /m**2 /s)
   psncanopy           => pcf%psncanopy   , & ! Output: [real(r8) (:,:)] canopy scale photosynthesis umol CO2 /m**2/ s
   ED_patch            => EDpft%ED_patch  , & ! Input:  [real(r8) (:,:)] does this 'p' have any vegetation associated with it?
   pgridcell           => pft%gridcell      & ! Input:  [real(r8) (:,:)] gridcell of 'p'
   )
    
    fpsn(p) = psncanopy(p)
         if(ED_patch(p) == 1)then
            g = pgridcell(p)
            currentPatch => gridCellEdState(g)%spnt%oldest_patch   
            do while(p /= currentPatch%clm_pno)
               currentPatch => currentPatch%younger
            enddo

            currentCohort => currentPatch%shortest

            do while(associated(currentCohort))

               ! Accumulate fluxes from hourly to daily values. 
               ! _clm fluxes are KgC/indiv/timestep _acc are KgC/indiv/day
               currentCohort%npp_acc  = currentCohort%npp_acc  + currentCohort%npp_clm 
               currentCohort%gpp_acc  = currentCohort%gpp_acc  + currentCohort%gpp_clm 
               currentCohort%resp_acc = currentCohort%resp_acc + currentCohort%resp_clm

               do iv=1,currentCohort%nv
                  if(currentCohort%year_net_uptake(iv) == 999._r8)then ! note that there were leaves in this layer this year. 
                     currentCohort%year_net_uptake(iv) = 0._r8
                  end if
                  currentCohort%year_net_uptake(iv) = currentCohort%year_net_uptake(iv) + currentCohort%ts_net_uptake(iv)
               enddo

               currentCohort => currentCohort%taller
            enddo ! while(associated(currentCohort)

         end if !ED_patch
     end associate
end subroutine AccumulateFluxes_ED

end module EDAccumulateFluxesMod
