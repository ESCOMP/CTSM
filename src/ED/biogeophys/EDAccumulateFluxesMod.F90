module EDAccumulateFluxesMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This routine accumulates NPP, GPP and respiration of each cohort over the course of each 24 hour period. 
  ! The fluxes are stored per cohort, and the npp_tstep (etc) fluxes are calcualted in EDPhotosynthesis
  ! This routine cannot be in EDPhotosynthesis because EDPhotosynthesis is a loop and therefore would
  ! erroneously add these things up multiple times. 
  ! Rosie Fisher. March 2014. 
  !
  ! !USES:
  implicit none
  private
  !
  public :: AccumulateFluxes_ED

  logical :: DEBUG = .false.  ! for debugging this module
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine AccumulateFluxes_ED(nsites, sites, bc_in, bc_out)
    !
    ! !DESCRIPTION:
    ! see above
    !
    ! !USES:
    use shr_kind_mod      , only : r8 => shr_kind_r8
    use clm_varctl        , only : iulog
    use EDTypesMod        , only : ed_patch_type, ed_cohort_type, ed_site_type
    use FatesInterfaceMod , only : bc_in_type,bc_out_type
    !
    ! !ARGUMENTS    
    integer,            intent(in)            :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_in_type),   intent(in)            :: bc_in(nsites)
    type(bc_out_type),  intent(inout)         :: bc_out(nsites)
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer  :: ccohort ! current cohort
    type(ed_patch_type) , pointer  :: cpatch ! current patch
    integer :: iv !leaf layer
    integer :: c  ! clm/alm column
    integer :: s  ! ed site
    integer :: ifp ! index fates patch
    !----------------------------------------------------------------------
    
    do s = 1, nsites
       
       ifp = 0
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))                 
          ifp = ifp+1

          if( bc_in(s)%filter_photo_pa(ifp) == 3 ) then
             ccohort => cpatch%shortest
             do while(associated(ccohort))
                
                ! Accumulate fluxes from hourly to daily values. 
                ! _tstep fluxes are KgC/indiv/timestep _acc are KgC/indiv/day
                
                if ( DEBUG ) then
                   write(iulog,*) 'EDAccumFlux 64 ',ccohort%npp_acc, &
                         ccohort%npp_tstep
                   write(iulog,*) 'EDAccumFlux 66 ',ccohort%gpp_tstep
                   write(iulog,*) 'EDAccumFlux 67 ',ccohort%resp_tstep
                endif
                
                ccohort%npp_acc  = ccohort%npp_acc  + ccohort%npp_tstep 
                ccohort%gpp_acc  = ccohort%gpp_acc  + ccohort%gpp_tstep 
                ccohort%resp_acc = ccohort%resp_acc + ccohort%resp_tstep
                
                do iv=1,ccohort%nv
                   if(ccohort%year_net_uptake(iv) == 999._r8)then ! note that there were leaves in this layer this year. 
                      ccohort%year_net_uptake(iv) = 0._r8
                   end if
                   ccohort%year_net_uptake(iv) = ccohort%year_net_uptake(iv) + ccohort%ts_net_uptake(iv)
                enddo
                
                ccohort => ccohort%taller
             enddo ! while(associated(ccohort))
          end if
          cpatch => cpatch%younger
       end do  ! while(associated(cpatch))
    end do
    return
    
 end subroutine AccumulateFluxes_ED

end module EDAccumulateFluxesMod

