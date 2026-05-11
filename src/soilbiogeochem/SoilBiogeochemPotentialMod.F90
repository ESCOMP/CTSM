module SoilBiogeochemPotentialMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate potential decomp rates and total immobilization demand.
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use decompMod                          , only : bounds_type
  use clm_varpar                         , only : nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar                         , only : i_cop_mic, i_oli_mic
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, mimics_decomp, decomp_method
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType    , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType     , only : soilbiogeochem_nitrogenflux_type
  use clm_varctl                         , only : use_fates, iulog
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams
  public :: SoilBiogeochemPotential
  !
  type, private :: params_type
     real(r8) :: dnp         !denitrification proportion
  end type Params_type
  !
  type(params_type), private :: params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read parameters
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io
    use abortutils   , only: endrun
    use shr_log_mod  , only: errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNDecompParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    tString='dnp'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%dnp=tempr 

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemPotential (bounds, num_bgc_soilc, filter_bgc_soilc, &
       soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,  &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst, &
       cn_decomp_pools, p_decomp_cpool_loss, p_decomp_cn_gain, &
       pmnf_decomp_cascade, p_decomp_npool_to_din)
    !
    ! !USES:
    use shr_log_mod                       , only : errMsg => shr_log_errMsg
    use SoilBiogeochemDecompCascadeConType, only : i_atm
    !
    ! !ARGUMENT:
    type(bounds_type)                       , intent(in)    :: bounds   
    integer                                 , intent(in)    :: num_bgc_soilc          ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:)    ! filter for soil columns
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonstate_type)   , intent(in)    :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    real(r8)                                , intent(out)   :: cn_decomp_pools(bounds%begc:,1:,1:)     ! c:n ratios of applicable pools
    real(r8)                                , intent(out)   :: p_decomp_cpool_loss(bounds%begc:,1:,1:) ! potential C loss from one pool to another
    real(r8)                                , intent(out)   :: pmnf_decomp_cascade(bounds%begc:,1:,1:) ! potential mineral N flux, from one pool to another
    real(r8)                                , intent(out)   :: p_decomp_npool_to_din(bounds%begc:,1:,1:)  ! potential flux to dissolved inorganic N
    real(r8)                                , intent(out)   :: p_decomp_cn_gain(bounds%begc:,1:,1:)  ! C:N ratio of the flux gained by the receiver pool
    !
    ! !LOCAL VARIABLES:
    integer :: c,j,k,l,m                                    !indices
    integer :: fc                                           !filter column index
    integer :: begc,endc                                    !bounds 
    real(r8):: immob(bounds%begc:bounds%endc,1:nlevdecomp)  !potential N immobilization
    real(r8):: p_decomp_cpool_gain(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions)  ! potential C gain by receiver pool
    real(r8):: p_decomp_npool_gain(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions)  ! potential N gain by receiver pool
    real(r8):: p_decomp_cpool_gain_sum(1:ndecomp_pools)  ! total potential C gain by receiver pool (only microbial pools)
    real(r8):: p_decomp_npool_gain_sum(1:ndecomp_pools)  ! total potential N gain by receiver pool (only microbial pools)
    real(r8):: decomp_nc_loss_donor  ! N:C ratio of donor pool
    real(r8):: p_decomp_cn_diff_ratio  ! relative change in receiver pool C:N
    real(r8):: p_decomp_npool_loss  ! potential N flux out of donor pool
    real(r8):: ratio                                        !temporary variable
    !-----------------------------------------------------------------------
   
    begc = bounds%begc; endc = bounds%endc

    SHR_ASSERT_ALL_FL((ubound(cn_decomp_pools)     == (/endc,nlevdecomp,ndecomp_pools/))               , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(p_decomp_cn_gain)    == (/endc,nlevdecomp,ndecomp_pools/))               , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(p_decomp_cpool_loss) == (/endc,nlevdecomp,ndecomp_cascade_transitions/)) , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(pmnf_decomp_cascade) == (/endc,nlevdecomp,ndecomp_cascade_transitions/)) , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(p_decomp_npool_to_din) == (/endc,nlevdecomp,ndecomp_cascade_transitions/)) , sourcefile, __LINE__)

    associate(                                                                                                          & 
         cascade_donor_pool               => decomp_cascade_con%cascade_donor_pool                                 , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool            => decomp_cascade_con%cascade_receiver_pool                              , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step
         floating_cn_ratio_decomp_pools   => decomp_cascade_con%floating_cn_ratio_decomp_pools                     , & ! Input:  [logical  (:)     ]  TRUE => pool has fixed C:N ratio                   
         initial_cn_ratio                 => decomp_cascade_con%initial_cn_ratio                                   , & ! Input:  [real(r8) (:)     ]  c:n ratio for initialization of pools             

         nue_decomp_cascade               => soilbiogeochem_state_inst%nue_decomp_cascade_col                      , & ! Input:  [real(r8) (:)     ]  N use efficiency for a given transition (gN going into microbe / gN decomposed)
         rf_decomp_cascade                => soilbiogeochem_carbonflux_inst%rf_decomp_cascade_col                  , & ! Input:  [real(r8) (:,:,:) ]  respired fraction in decomposition step (frac)
         pathfrac_decomp_cascade          => soilbiogeochem_carbonflux_inst%pathfrac_decomp_cascade_col            , & ! Input:  [real(r8) (:,:,:) ]  what fraction of C passes from donor to receiver pool through a given transition (frac)

         decomp_npools_vr                 => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col                , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools

         decomp_cpools_vr                 => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col                  , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools

         potential_immob_vr               => soilbiogeochem_nitrogenflux_inst%potential_immob_vr_col               , & ! Output: [real(r8) (:,:)   ]                                                  
         gross_nmin_vr                    => soilbiogeochem_nitrogenflux_inst%gross_nmin_vr_col                    , & ! Output: [real(r8) (:,:)   ]                                                  
         
         cn_col                           => soilbiogeochem_carbonflux_inst%cn_col                                 , & ! Input: [real(r8) (:,:)   ]  C:N ratio
         decomp_k                         => soilbiogeochem_carbonflux_inst%decomp_k_col                           , & ! Input: [real(r8) (:,:,:) ]  decomposition rate coefficient (1./sec)
         phr_vr                           => soilbiogeochem_carbonflux_inst%phr_vr_col                               & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)                           
         )

      
      ! set initial values for potential C and N fluxes
      p_decomp_cpool_loss(begc:endc, :, :) = 0._r8
      pmnf_decomp_cascade(begc:endc, :, :) = 0._r8

      ! column loop to calculate potential decomp rates and total immobilization demand

      !! calculate c:n ratios of applicable pools
      do l = 1, ndecomp_pools
         if ( floating_cn_ratio_decomp_pools(l) ) then
            do j = 1,nlevdecomp
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if ( decomp_npools_vr(c,j,l) > 0._r8 ) then
                     cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
                  end if
               end do
            end do
         else
            do j = 1,nlevdecomp
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
               end do
            end do
         end if
      end do

      ! calculate the non-nitrogen-limited fluxes
      ! these fluxes include the  "/ dt" term to put them on a
      ! per second basis, since the rate constants have been
      ! calculated on a per timestep basis.

      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)

               if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. &
                    decomp_k(c,j,cascade_donor_pool(k)) > 0._r8 ) then
                  p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) &
                       * decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)

                  if ( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)) ) then  !! not transition of cwd to litter

                     if (cascade_receiver_pool(k) /= i_atm ) then  ! not 100% respiration
                        ratio = 0._r8

                        if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8) then
                           ratio = cn_decomp_pools(c,j,cascade_receiver_pool(k))/cn_decomp_pools(c,j,cascade_donor_pool(k))
                        endif

                        pmnf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(c,j,k) * (1.0_r8 - rf_decomp_cascade(c,j,k) - ratio) &
                             / cn_decomp_pools(c,j,cascade_receiver_pool(k)) )

                     else   ! 100% respiration
                        pmnf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
                     endif
                  else   ! CWD -> litter OR mimics_decomp is true
                     pmnf_decomp_cascade(c,j,k) = 0._r8

                     if (decomp_method == mimics_decomp) then
                        ! N:C ratio of donor pools (N:C instead of C:N because
                        ! already checked that we're not dividing by zero)
                        decomp_nc_loss_donor = &
                           decomp_npools_vr(c,j,cascade_donor_pool(k)) / &
                           decomp_cpools_vr(c,j,cascade_donor_pool(k))
                        ! calculate N fluxes from donor pools
                        p_decomp_npool_loss = &
                           p_decomp_cpool_loss(c,j,k) * decomp_nc_loss_donor
                        ! Track N lost to DIN from messy eating, ie the
                        ! diff between p_decomp_npool_loss (pertains to
                        ! donor) and p_decomp_npool_gain (receiver)
                        p_decomp_cpool_gain(c,j,k) = &
                           p_decomp_cpool_loss(c,j,k) * (1.0_r8 - rf_decomp_cascade(c,j,k))
                        p_decomp_npool_gain(c,j,k) = &
                           p_decomp_npool_loss * nue_decomp_cascade(k)
                        p_decomp_npool_to_din(c,j,k) = &
                           p_decomp_npool_loss - p_decomp_npool_gain(c,j,k)
                     end if  ! using mimics
                  end if  ! floating_cn_ratio_decomp_pools = .true.
               else  ! donors not donating
                  p_decomp_cpool_gain(c,j,k) = 0.0_r8
                  p_decomp_npool_gain(c,j,k) = 0.0_r8
                  p_decomp_npool_to_din(c,j,k) = 0.0_r8
               end if  ! donors donating? (ie decomp_cpools_vr & decomp_k > 0)
            end do  ! column loop

         end do  ! soil level loop
      end do  ! transitions loop

      ! Calculate cn_gain into microbial biomass (in first do k
      ! transitions loop).
      ! Compare cn_gain to target C:N ratio of microbial biomass pools
      ! to determine immobilization vs. mineralization (in second do k
      ! transitions loop).
      if (decomp_method == mimics_decomp) then
         do j = 1,nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               ! Sum C & N fluxes from all transitions into m1 & m2 pools.
               ! Had to form a new loop for the summation due to the order
               ! necessary, ie do k as the innermost loop.
               p_decomp_cpool_gain_sum(:) = 0.0_r8
               p_decomp_npool_gain_sum(:) = 0.0_r8
               do k = 1, ndecomp_cascade_transitions
                  if (cascade_receiver_pool(k) == i_cop_mic .or. &
                      cascade_receiver_pool(k) == i_oli_mic) then
                     if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. &
                         decomp_k(c,j,cascade_donor_pool(k)) > 0._r8 ) then
                        p_decomp_cpool_gain_sum(cascade_receiver_pool(k)) = &
                           p_decomp_cpool_gain_sum(cascade_receiver_pool(k)) + &
                           p_decomp_cpool_gain(c,j,k)
                        p_decomp_npool_gain_sum(cascade_receiver_pool(k)) = &
                           p_decomp_npool_gain_sum(cascade_receiver_pool(k)) + &
                           p_decomp_npool_gain(c,j,k)
                        ! Once k-loop completes, the left hand side should end
                        ! up with the correct cn ratio
                        if (p_decomp_npool_gain_sum(cascade_receiver_pool(k)) > 0.0_r8) then
                           p_decomp_cn_gain(c,j,cascade_receiver_pool(k)) = &
                              p_decomp_cpool_gain_sum(cascade_receiver_pool(k)) / &
                              p_decomp_npool_gain_sum(cascade_receiver_pool(k))
                        else
                           p_decomp_cn_gain(c,j,cascade_receiver_pool(k)) = 0.0_r8
                        end if  ! denominator check
                     end if  ! donors donating (decomp_cpools_vr & decomp_k > 0)
                  end if  ! microbes receiving
               end do  ! transitions loop
               do k = 1, ndecomp_cascade_transitions
                  if (cascade_receiver_pool(k) == i_cop_mic .or. &
                      cascade_receiver_pool(k) == i_oli_mic) then
                     if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. &
                         decomp_k(c,j,cascade_donor_pool(k)) > 0._r8 ) then
                        ! if p_decomp_cn_diff < 0  N mineralization
                        !                     > 0  immobilization
                        ! "min" in next line turns off immobilization flux
                        p_decomp_cn_diff_ratio = min(0.0_r8,  &
                           (p_decomp_cn_gain(c,j,cascade_receiver_pool(k)) - &
                            cn_col(c,cascade_receiver_pool(k))) / cn_col(c,cascade_receiver_pool(k)))
                        ! Actual amount of N that's mineralized or that would
                        ! need to be immobilized
                        ! negative=mineralization: add to the DIN pool
                        ! positive=immobilizaiton: compete for N with plants to
                        !                          see how much we get
                        pmnf_decomp_cascade(c,j,k) = p_decomp_cn_diff_ratio * p_decomp_npool_gain(c,j,k)
                     end if  ! donors donating (decomp_cpools_vr & decomp_k > 0)
                  end if  ! microbes receiving
               end do  ! transitions loop
            end do  ! column loop
         end do  ! soil levels loop
      end if  ! using mimics

      ! Sum up all the potential immobilization fluxes (positive pmnf flux)
      ! and all the mineralization fluxes (negative pmnf flux)
      do j = 1,nlevdecomp
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            immob(c,j) = 0._r8
         end do
      end do
      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (pmnf_decomp_cascade(c,j,k) > 0._r8) then
                  immob(c,j) = immob(c,j) + pmnf_decomp_cascade(c,j,k)
               else
                  gross_nmin_vr(c,j) = gross_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
               end if
               if (decomp_method == mimics_decomp) then
                  gross_nmin_vr(c,j) = gross_nmin_vr(c,j) + p_decomp_npool_to_din(c,j,k)
               end if
            end do
         end do
      end do

      do j = 1,nlevdecomp
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            potential_immob_vr(c,j) = immob(c,j)
         end do
      end do

      ! Add up potential hr for methane calculations
      do j = 1,nlevdecomp
         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            phr_vr(c,j) = 0._r8
         end do
      end do
      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               phr_vr(c,j) = phr_vr(c,j) + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
            end do
         end do
      end do


    end associate

  end subroutine SoilBiogeochemPotential
 
end module SoilBiogeochemPotentialMod
