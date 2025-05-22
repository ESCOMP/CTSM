module CNC14DecayMod

  !-----------------------------------------------------------------------
  ! Module for 14-carbon flux variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use clm_time_manager                   , only : get_step_size_real, get_average_days_per_year
  use clm_varpar                         , only : nlevdecomp, ndecomp_pools
  use clm_varcon                         , only : secspday
  use clm_varctl                         , only : spinup_state
  use CNSharedParamsMod                  , only : use_matrixcn
  use decompMod                          , only : bounds_type
  use pftconMod                          , only : npcropmin
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, use_soil_matrixcn
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use PatchType                          , only : patch
  use ColumnType                         , only : col
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: C14Decay
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine C14Decay( bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst, &
       c14_cnveg_carbonflux_inst,  c14_soilbiogeochem_carbonflux_inst )
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the radioactive decay of C14
    !
    implicit none
    ! !ARGUMENTS:
    type(bounds_type)                     , intent(in)    :: bounds        
    integer                               , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                               , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                               , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                               , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(CNVeg_carbonstate_type)          , intent(inout) :: c14_cnveg_carbonstate_inst
    type(soilbiogeochem_carbonstate_type) , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(CNVeg_carbonflux_type)           , intent(inout) :: c14_cnveg_carbonflux_inst
    type(soilbiogeochem_carbonflux_type)  , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,j,g,l,p,fc,c,i
    real(r8) :: dt            ! radiation time step (seconds)
    real(r8) :: half_life
    real(r8) :: decay_const
    real(r8) :: days_per_year ! days per year
    real(r8) :: spinup_term   ! spinup accelerated decomposition factor, used to accelerate transport as well
    !-----------------------------------------------------------------------

    associate(                                                                               & 
         spinup_factor      =>    decomp_cascade_con%spinup_factor                         , & ! Input:   [real(r8) (:)     ]  factor for AD spinup associated with each pool

         decomp_cpools_vr   =>    c14_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col , & ! Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools

         cropseedc_deficit  =>    c14_cnveg_carbonstate_inst%cropseedc_deficit_patch       , & ! Output:  [real(r8) (:)     ]                                          
         seedc              =>    c14_cnveg_carbonstate_inst%seedc_grc                     , & ! Output:  [real(r8) (:)     ]                                          
         cpool              =>    c14_cnveg_carbonstate_inst%cpool_patch                   , & ! Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool  
         xsmrpool           =>    c14_cnveg_carbonstate_inst%xsmrpool_patch                , & ! Output:  [real(r8) (:)     ]  (gC/m2) execss maint resp C pool        
         deadcrootc         =>    c14_cnveg_carbonstate_inst%deadcrootc_patch              , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C              
         deadcrootc_storage =>    c14_cnveg_carbonstate_inst%deadcrootc_storage_patch      , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage      
         deadcrootc_xfer    =>    c14_cnveg_carbonstate_inst%deadcrootc_xfer_patch         , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer     
         deadstemc          =>    c14_cnveg_carbonstate_inst%deadstemc_patch               , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                     
         deadstemc_storage  =>    c14_cnveg_carbonstate_inst%deadstemc_storage_patch       , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage             
         deadstemc_xfer     =>    c14_cnveg_carbonstate_inst%deadstemc_xfer_patch          , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer            
         frootc             =>    c14_cnveg_carbonstate_inst%frootc_patch                  , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C                     
         frootc_storage     =>    c14_cnveg_carbonstate_inst%frootc_storage_patch          , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage             
         frootc_xfer        =>    c14_cnveg_carbonstate_inst%frootc_xfer_patch             , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer            
         gresp_storage      =>    c14_cnveg_carbonstate_inst%gresp_storage_patch           , & ! Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage      
         gresp_xfer         =>    c14_cnveg_carbonstate_inst%gresp_xfer_patch              , & ! Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer     
         leafc              =>    c14_cnveg_carbonstate_inst%leafc_patch                   , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C                          
         leafc_storage      =>    c14_cnveg_carbonstate_inst%leafc_storage_patch           , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                  
         leafc_xfer         =>    c14_cnveg_carbonstate_inst%leafc_xfer_patch              , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                 
         livecrootc         =>    c14_cnveg_carbonstate_inst%livecrootc_patch              , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C              
         livecrootc_storage =>    c14_cnveg_carbonstate_inst%livecrootc_storage_patch      , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage      
         livecrootc_xfer    =>    c14_cnveg_carbonstate_inst%livecrootc_xfer_patch         , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer     
         livestemc          =>    c14_cnveg_carbonstate_inst%livestemc_patch               , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C                     
         livestemc_storage  =>    c14_cnveg_carbonstate_inst%livestemc_storage_patch       , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage             
         livestemc_xfer     =>    c14_cnveg_carbonstate_inst%livestemc_xfer_patch          , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer            
         pft_ctrunc         =>    c14_cnveg_carbonstate_inst%ctrunc_patch                  , & ! Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation 

         ileaf_to_iout_fic        => c14_cnveg_carbonflux_inst%ileaf_to_iout_fi                  , & ! Input: [integer (:)] Index of fire related C transfer from leaf pool to outside of vegetation pools
         ileafst_to_iout_fic      => c14_cnveg_carbonflux_inst%ileafst_to_iout_fi                , & ! Input: [integer (:)] Index of fire related C transfer from leaf storage pool to outside of vegetation pools
         ileafxf_to_iout_fic      => c14_cnveg_carbonflux_inst%ileafxf_to_iout_fi                , & ! Input: [integer (:)] Index of fire related C transfer from leaf transfer pool to outside of vegetation pools
         ifroot_to_iout_fic       => c14_cnveg_carbonflux_inst%ifroot_to_iout_fi                 , & ! Input: [integer (:)] Index of fire related C transfer from fine root pool to outside of vegetation pools
         ifrootst_to_iout_fic     => c14_cnveg_carbonflux_inst%ifrootst_to_iout_fi               , & ! Input: [integer (:)] Index of fire related C transfer from fine root storage pool to outside of vegetation pools
         ifrootxf_to_iout_fic     => c14_cnveg_carbonflux_inst%ifrootxf_to_iout_fi               , & ! Input: [integer (:)] Index of fire related C transfer from fine root transfer pool to outside of vegetation pools
         ilivestem_to_iout_fic    => c14_cnveg_carbonflux_inst%ilivestem_to_iout_fi              , & ! Input: [integer (:)] Index of fire related C transfer from live stem pool to outside of vegetation pools
         ilivestemst_to_iout_fic  => c14_cnveg_carbonflux_inst%ilivestemst_to_iout_fi            , & ! Input: [integer (:)] Index of fire related C transfer from live stem storage pool to outside of vegetation pools
         ilivestemxf_to_iout_fic  => c14_cnveg_carbonflux_inst%ilivestemxf_to_iout_fi            , & ! Input: [integer (:)] Index of fire related C transfer from live stem transfer pool to outside of vegetation pools
         ideadstem_to_iout_fic    => c14_cnveg_carbonflux_inst%ideadstem_to_iout_fi              , & ! Input: [integer (:)] Index of fire related C transfer from dead stem pool to outside of vegetation pools
         ideadstemst_to_iout_fic  => c14_cnveg_carbonflux_inst%ideadstemst_to_iout_fi            , & ! Input: [integer (:)] Index of fire related C transfer from dead stem storage pool to outside of vegetation pools
         ideadstemxf_to_iout_fic  => c14_cnveg_carbonflux_inst%ideadstemxf_to_iout_fi            , & ! Input: [integer (:)] Index of fire related C transfer from dead stem transfer pool to outside of vegetation pools
         ilivecroot_to_iout_fic   => c14_cnveg_carbonflux_inst%ilivecroot_to_iout_fi             , & ! Input: [integer (:)] Index of fire related C transfer from live coarse root pool to outside of vegetation pools
         ilivecrootst_to_iout_fic => c14_cnveg_carbonflux_inst%ilivecrootst_to_iout_fi           , & ! Input: [integer (:)] Index of fire related C transfer from live coarse root storage pool to outside of vegetation pools
         ilivecrootxf_to_iout_fic => c14_cnveg_carbonflux_inst%ilivecrootxf_to_iout_fi           , & ! Input: [integer (:)] Index of fire related C transfer from live coarse root transfer pool to outside of vegetation pools
         ideadcroot_to_iout_fic   => c14_cnveg_carbonflux_inst%ideadcroot_to_iout_fi             , & ! Input: [integer (:)] Index of fire related C transfer from dead coarse root pool to outside of vegetation pools
         ideadcrootst_to_iout_fic => c14_cnveg_carbonflux_inst%ideadcrootst_to_iout_fi           , & ! Input: [integer (:)] Index of fire related C transfer from dead coarse root storage pool to outside of vegetation pools
         ideadcrootxf_to_iout_fic => c14_cnveg_carbonflux_inst%ideadcrootxf_to_iout_fi             & ! Input: [integer (:)] Index of fire related C transfer from dead coarse root transfer pool to outside of vegetation pools
         )

      ! set time steps
      dt = get_step_size_real()
      days_per_year = get_average_days_per_year()

      half_life = 5730._r8 * secspday * days_per_year
      decay_const = - log(0.5_r8) / half_life

      do g = bounds%begg, bounds%endg
         seedc(g) = seedc(g) * (1._r8 - decay_const * dt)
      end do

      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( spinup_state >= 1) then
                  ! speed up radioactive decay by the same factor as decomposition so tat SOM ages prematurely in all respects
                  spinup_term = spinup_factor(l)
                  if ( abs(spinup_factor(l) - 1._r8) .gt. .000001_r8 ) then
                     spinup_term = spinup_term  * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
                  endif
               else
                  spinup_term = 1._r8
               endif
               if(.not. use_soil_matrixcn)then
                  decomp_cpools_vr(c,j,l) = decomp_cpools_vr(c,j,l) * (1._r8 - decay_const * spinup_term * dt)
               else
                  associate( &
                    matrix_decomp_fire_k => c14_soilbiogeochem_carbonflux_inst%matrix_decomp_fire_k_col   & ! Output: [real(r8) (:,:)   ]  (gC/m3/step) VR deomp. C fire loss in matrix representation
                  )
                  matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) = matrix_decomp_fire_k(c,j+nlevdecomp*(l-1)) &
                                                             - spinup_term * decay_const * dt
                  end associate
               end if
            end do
         end do
      end do ! end of columns loop

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         cpool(p)              = cpool(p)               * (1._r8 - decay_const * dt)
         xsmrpool(p)           = xsmrpool(p)            * (1._r8 - decay_const * dt)
         if(.not. use_matrixcn)then
            ! NOTE: Any changes here need to be applied below
            deadcrootc(p)         = deadcrootc(p)          * (1._r8 - decay_const * dt)
            deadcrootc_storage(p) = deadcrootc_storage(p)  * (1._r8 - decay_const * dt)
            deadcrootc_xfer(p)    = deadcrootc_xfer(p)     * (1._r8 - decay_const * dt)
            deadstemc(p)          = deadstemc(p)           * (1._r8 - decay_const * dt)
            deadstemc_storage(p)  = deadstemc_storage(p)   * (1._r8 - decay_const * dt)
            deadstemc_xfer(p)     = deadstemc_xfer(p)      * (1._r8 - decay_const * dt)
            frootc(p)             = frootc(p)              * (1._r8 - decay_const * dt)
            frootc_storage(p)     = frootc_storage(p)      * (1._r8 - decay_const * dt)
            frootc_xfer(p)        = frootc_xfer(p)         * (1._r8 - decay_const * dt)
            leafc(p)              = leafc(p)               * (1._r8 - decay_const * dt)
            leafc_storage(p)      = leafc_storage(p)       * (1._r8 - decay_const * dt)
            leafc_xfer(p)         = leafc_xfer(p)          * (1._r8 - decay_const * dt)
            livecrootc(p)         = livecrootc(p)          * (1._r8 - decay_const * dt)
            livecrootc_storage(p) = livecrootc_storage(p)  * (1._r8 - decay_const * dt)
            livecrootc_xfer(p)    = livecrootc_xfer(p)     * (1._r8 - decay_const * dt)
            livestemc(p)          = livestemc(p)           * (1._r8 - decay_const * dt)
            livestemc_storage(p)  = livestemc_storage(p)   * (1._r8 - decay_const * dt)
            livestemc_xfer(p)     = livestemc_xfer(p)      * (1._r8 - decay_const * dt)
         else
            associate( &
              matrix_fitransfer => c14_cnveg_carbonflux_inst%matrix_fitransfer_patch & ! Output: [real(r8) (:,:)   ] (gC/m2/s) C transfer rate from fire processes
            )
            ! Each of these MUST correspond to the code above. Any changes in
            ! code above need to apply here as well
            matrix_fitransfer(p,ideadcroot_to_iout_fic)   = decay_const
            matrix_fitransfer(p,ideadcrootst_to_iout_fic) = decay_const
            matrix_fitransfer(p,ideadcrootxf_to_iout_fic) = decay_const
            matrix_fitransfer(p,ideadstem_to_iout_fic)    = decay_const
            matrix_fitransfer(p,ideadstemst_to_iout_fic)  = decay_const
            matrix_fitransfer(p,ideadstemxf_to_iout_fic)  = decay_const
            matrix_fitransfer(p,ifroot_to_iout_fic)       = decay_const
            matrix_fitransfer(p,ifrootst_to_iout_fic)     = decay_const
            matrix_fitransfer(p,ifrootxf_to_iout_fic)     = decay_const
            matrix_fitransfer(p,ileaf_to_iout_fic)        = decay_const
            matrix_fitransfer(p,ileafst_to_iout_fic)      = decay_const
            matrix_fitransfer(p,ileafxf_to_iout_fic)      = decay_const
            matrix_fitransfer(p,ilivecroot_to_iout_fic)   = decay_const
            matrix_fitransfer(p,ilivecrootst_to_iout_fic) = decay_const
            matrix_fitransfer(p,ilivecrootxf_to_iout_fic) = decay_const
            matrix_fitransfer(p,ilivestem_to_iout_fic)    = decay_const
            matrix_fitransfer(p,ilivestemst_to_iout_fic)  = decay_const
            matrix_fitransfer(p,ilivestemxf_to_iout_fic)  = decay_const
            end associate
         end if
         gresp_storage(p)      = gresp_storage(p)       * (1._r8 - decay_const * dt)
         gresp_xfer(p)         = gresp_xfer(p)          * (1._r8 - decay_const * dt)
         pft_ctrunc(p)         = pft_ctrunc(p)          * (1._r8 - decay_const * dt)

         ! NOTE(wjs, 2017-02-02) This isn't a completely robust way to check if this is a
         ! prognostic crop patch (at the very least it should also check if <= npcropmax;
         ! ideally it should use a prognostic_crop flag that doesn't seem to exist
         ! currently). But I'm just being consistent with what's done elsewhere (e.g., in
         ! CStateUpdate1).
         if (patch%itype(p) >= npcropmin) then ! skip 2 generic crops
            cropseedc_deficit(p)  = cropseedc_deficit(p)   * (1._r8 - decay_const * dt)
         end if
      end do

    end associate

  end subroutine C14Decay

end module CNC14DecayMod
