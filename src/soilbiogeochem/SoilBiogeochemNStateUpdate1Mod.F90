module SoilBiogeochemNStateUpdate1Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                       , only: r8 => shr_kind_r8
  use clm_time_manager                   , only : get_step_size_real
  use clm_varpar                         , only : nlevdecomp, ndecomp_cascade_transitions
  use clm_varctl                         , only : iulog, use_nitrif_denitrif, use_crop
  use SoilBiogeochemDecompCascadeConType , only : use_soil_matrixcn
  use clm_varcon                         , only : nitrif_n2o_loss_frac
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemNitrogenStateType    , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenfluxType     , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use CNSharedParamsMod                  , only : use_fun
  use ColumnType                         , only : col 
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: SoilBiogeochemNStateUpdate1
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemNStateUpdate1(num_bgc_soilc, filter_bgc_soilc,  &
       soilbiogeochem_state_inst, soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                                 , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
    type(soilbiogeochem_state_type)         , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k ! indices
    integer :: fc        ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)

    !-----------------------------------------------------------------------

    associate(                                                                   & 
         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool        , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool     , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step

         ndep_prof             => soilbiogeochem_state_inst%ndep_prof_col      , & ! Input:  [real(r8) (:,:)   ]  profile over which N deposition is distributed through column (1/m)
         nfixation_prof        => soilbiogeochem_state_inst%nfixation_prof_col , & ! Input:  [real(r8) (:,:)   ]  profile over which N fixation is distributed through column (1/m)

         nf                    => soilbiogeochem_nitrogenflux_inst             , & ! Output:
         ns                    => soilbiogeochem_nitrogenstate_inst              & ! Output:
         )

      ! set time steps
      dt = get_step_size_real()

      do j = 1, nlevdecomp
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            if(use_fun)then !RF in FUN logic, the fixed N goes straight into the plant, and not into the SMINN pool. 
 	               ! N deposition and fixation (put all into NH4 pool)
	               ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%ndep_to_sminn_col(c)*dt * ndep_prof(c,j)
	               ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%ffix_to_sminn_col(c)*dt * nfixation_prof(c,j)
	          else
	            if (.not. use_nitrif_denitrif) then

	               ! N deposition and fixation
	               ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + nf%ndep_to_sminn_col(c)*dt * ndep_prof(c,j)
	               ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + nf%nfix_to_sminn_col(c)*dt * nfixation_prof(c,j)

	            else

	               ! N deposition and fixation (put all into NH4 pool)
	               ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%ndep_to_sminn_col(c)*dt * ndep_prof(c,j)
	               ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%nfix_to_sminn_col(c)*dt * nfixation_prof(c,j)
                       
                end if
           end if
          
         end do
         
      end do

      ! repeating N dep and fixation for crops
      if ( use_crop )then
         do j = 1, nlevdecomp

            ! column loop
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (.not. use_nitrif_denitrif) then

                  ! N deposition and fixation
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) &
                       + nf%fert_to_sminn_col(c)*dt * ndep_prof(c,j)
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) &
                       + nf%soyfixn_to_sminn_col(c)*dt * nfixation_prof(c,j)

               else

                  ! N deposition and fixation (put all into NH4 pool)
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) &
                       + nf%fert_to_sminn_col(c)*dt * ndep_prof(c,j)
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) &
                       + nf%soyfixn_to_sminn_col(c)*dt * nfixation_prof(c,j)

               end if
            end do
         end do
      end if

      ! decomposition fluxes
      ! TODO slevis: improve indentation
   if (.not. use_soil_matrixcn) then
      do k = 1, ndecomp_cascade_transitions
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)

               nf%decomp_npools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                    nf%decomp_npools_sourcesink_col(c,j,cascade_donor_pool(k)) - &
                    nf%decomp_cascade_ntransfer_vr_col(c,j,k) * dt
            end do
         end do
      end do


      do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)

                  nf%decomp_npools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                       nf%decomp_npools_sourcesink_col(c,j,cascade_receiver_pool(k)) + &
                       (nf%decomp_cascade_ntransfer_vr_col(c,j,k) + &
                        nf%decomp_cascade_sminn_flux_vr_col(c,j,k)) * dt
               end do
            end do
         else  ! terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  nf%decomp_npools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                       nf%decomp_npools_sourcesink_col(c,j,cascade_donor_pool(k)) - &
                       nf%decomp_cascade_sminn_flux_vr_col(c,j,k) * dt
               end do
            end do
         end if
      end do
   end if  ! 

      if (.not. use_nitrif_denitrif) then

         !--------------------------------------------------------
         !-------------    NITRIF_DENITRIF OFF -------------------
         !--------------------------------------------------------

         ! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes and denitrification fluxes
         do k = 1, ndecomp_cascade_transitions
            if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               do j = 1, nlevdecomp
                  ! column loop
                  do fc = 1,num_bgc_soilc
                     c = filter_bgc_soilc(fc)
                     ns%sminn_vr_col(c,j)  = ns%sminn_vr_col(c,j) - &
                          (nf%sminn_to_denit_decomp_cascade_vr_col(c,j,k) + &
                          nf%decomp_cascade_sminn_flux_vr_col(c,j,k))* dt
                  end do
               end do
            else
               do j = 1, nlevdecomp
                  ! column loop
                  do fc = 1,num_bgc_soilc
                     c = filter_bgc_soilc(fc)
                     ns%sminn_vr_col(c,j)  = ns%sminn_vr_col(c,j) - &
                          nf%sminn_to_denit_decomp_cascade_vr_col(c,j,k)* dt

                     ns%sminn_vr_col(c,j)  = ns%sminn_vr_col(c,j) + &
                          nf%decomp_cascade_sminn_flux_vr_col(c,j,k)* dt

                  end do
               end do
            endif
         end do

         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               ! "bulk denitrification"
               ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) - nf%sminn_to_denit_excess_vr_col(c,j) * dt

               ! total plant uptake from mineral N
               if ( .not. use_fun ) then
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) - nf%sminn_to_plant_vr_col(c,j)*dt
               else
                  ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) - nf%sminn_to_plant_fun_vr_col(c,j)*dt
               end if
               ! flux that prevents N limitation (when Carbon_only is set)
               ns%sminn_vr_col(c,j) = ns%sminn_vr_col(c,j) + nf%supplement_to_sminn_vr_col(c,j)*dt
            end do
         end do

      else   

         !--------------------------------------------------------
         !-------------    NITRIF_DENITRIF ON --------------------
         !--------------------------------------------------------

         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)

               ! mineralization fluxes (divert a fraction of this stream to nitrification flux, add the rest to NH4 pool)
               ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%gross_nmin_vr_col(c,j)*dt

               ! immobilization fluxes
               ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) - nf%actual_immob_nh4_vr_col(c,j)*dt

               ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) - nf%actual_immob_no3_vr_col(c,j)*dt

               ! plant uptake fluxes
               if ( .not. use_fun )then
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) - nf%smin_nh4_to_plant_vr_col(c,j)*dt

                  ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) - nf%smin_no3_to_plant_vr_col(c,j)*dt
               else 
                  ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) -  nf%sminn_to_plant_fun_nh4_vr_col(c,j)*dt

                  ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) -  nf%sminn_to_plant_fun_no3_vr_col(c,j)*dt
               end if
             

               ! Account for nitrification fluxes
               ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) - nf%f_nit_vr_col(c,j) * dt

               ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) + nf%f_nit_vr_col(c,j) * dt &
                    * (1._r8 - nitrif_n2o_loss_frac)

               ! Account for denitrification fluxes
               ns%smin_no3_vr_col(c,j) = ns%smin_no3_vr_col(c,j) - nf%f_denit_vr_col(c,j) * dt

               ! flux that prevents N limitation (when Carbon_only is set; put all into NH4)
               ns%smin_nh4_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + nf%supplement_to_sminn_vr_col(c,j)*dt

               ! update diagnostic total
               ns%sminn_vr_col(c,j) = ns%smin_nh4_vr_col(c,j) + ns%smin_no3_vr_col(c,j)
               
            end do ! end of column loop
         end do
              
      end if

    end associate

  end subroutine SoilBiogeochemNStateUpdate1

end module SoilBiogeochemNStateUpdate1Mod
