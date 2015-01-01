module CNBalanceCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon/nitrogen mass balance checking.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_infnan_mod                  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use decompMod                       , only : bounds_type
  use abortutils                      , only : endrun
  use clm_varctl                      , only : iulog, use_nitrif_denitrif
  use clm_time_manager                , only : get_step_size
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use SoilBiogeochemNitrogenfluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemCarbonfluxType    , only : soilbiogeochem_carbonflux_type
  use ColumnType                      , only : col                
  use GridcellType                    , only : grc
  !
  implicit none
  private
  ! 
  ! !PUBLIC TYPES:
  type, public :: cn_balance_type
     private
     real(r8), pointer :: begcb_col(:)        ! (gC/m2) carbon mass, beginning of time step 
     real(r8), pointer :: endcb_col(:)        ! (gC/m2) carbon mass, end of time step
     real(r8), pointer :: begnb_col(:)        ! (gN/m2) nitrogen mass, beginning of time step 
     real(r8), pointer :: endnb_col(:)        ! (gN/m2) nitrogen mass, end of time step 
     logical , pointer :: beg_vals_set_col(:) ! Whether begcb/begnb have been set for this column in this time step
   contains
     procedure , public  :: Init
     procedure , public  :: BeginCNBalance
     procedure , public  :: CBalanceCheck
     procedure , public  :: NBalanceCheck
     procedure , private :: InitAllocate 
  end type cn_balance_type
  !
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    class(cn_balance_type)         :: this
    type(bounds_type) , intent(in) :: bounds  

    call this%InitAllocate(bounds)
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    class(cn_balance_type)         :: this
    type(bounds_type) , intent(in) :: bounds  

    integer :: begc, endc

    begc = bounds%begc; endc= bounds%endc

    allocate(this%begcb_col(begc:endc)) ; this%begcb_col(:) = nan
    allocate(this%endcb_col(begc:endc)) ; this%endcb_col(:) = nan
    allocate(this%begnb_col(begc:endc)) ; this%begnb_col(:) = nan
    allocate(this%endnb_col(begc:endc)) ; this%endnb_col(:) = nan
    allocate(this%beg_vals_set_col(begc:endc)) ; this%beg_vals_set_col(:) = .false.
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine BeginCNBalance(this, bounds, num_soilc, filter_soilc, &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Calculate beginning column-level carbon/nitrogen balance, for mass conservation check
    !
    ! !ARGUMENTS:
    class(cn_balance_type)         , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds          
    integer                        , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_carbonstate_type)   , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c
    !-----------------------------------------------------------------------

    associate(                                            & 
         col_begcb    => this%begcb_col                  , & ! Output: [real(r8) (:)]  (gC/m2) carbon mass, beginning of time step
         col_begnb    => this%begnb_col                  , & ! Output: [real(r8) (:)]  (gN/m2) nitrogen mass, beginning of time step 
         beg_vals_set => this%beg_vals_set_col           , & ! Output: [logical  (:)]  Whether begcb/begnb have been set
         totcolc      => cnveg_carbonstate_inst%totc_col , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
         totcoln      => cnveg_nitrogenstate_inst%totn_col & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
         )
    
    beg_vals_set(bounds%begc:bounds%endc) = .false.

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       col_begcb(c) = totcolc(c)
       col_begnb(c) = totcoln(c)
       beg_vals_set(c) = .true.
    end do

    end associate

  end subroutine BeginCNBalance
 
  !-----------------------------------------------------------------------
  subroutine CBalanceCheck(this, bounds, num_soilc, filter_soilc, &
       soilbiogeochem_carbonflux_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Perform carbon mass conservation check for column and patch
    !
    ! !ARGUMENTS:
    class(cn_balance_type)               , intent(inout) :: this
    type(bounds_type)                    , intent(in)    :: bounds          
    integer                              , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                              , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_carbonflux_type) , intent(in)    :: soilbiogeochem_carbonflux_inst
    type(cnveg_carbonflux_type)          , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: cnveg_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,err_index    ! indices
    integer  :: fc             ! lake filter indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    real(r8) :: col_cinputs
    real(r8) :: col_coutputs
    real(r8) :: col_errcb(bounds%begc:bounds%endc) 
    !-----------------------------------------------------------------------

    associate(                                                                            & 
         col_begcb               =>    this%begcb_col                                   , & ! Input:  [real(r8) (:) ]  (gC/m2) carbon mass, beginning of time step 
         col_endcb               =>    this%endcb_col                                   , & ! Output: [real(r8) (:) ]  (gC/m2) carbon mass, end of time step 
         beg_vals_set            =>    this%beg_vals_set_col                            , & ! Input:  [logical  (:) ]  Whether begcb/begnb have been set in this time step
         dwt_closs               =>    cnveg_carbonflux_inst%dwt_closs_col              , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total carbon loss from product pools and conversion
         product_closs           =>    cnveg_carbonflux_inst%product_closs_col          , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total wood product carbon loss
         gpp                     =>    cnveg_carbonflux_inst%gpp_col                    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) gross primary production      
         er                      =>    cnveg_carbonflux_inst%er_col                     , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs          =>    cnveg_carbonflux_inst%fire_closs_col             , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total column-level fire C loss
         col_hrv_xsmrpool_to_atm =>    cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_col    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool harvest mortality 

         som_c_leached           =>    soilbiogeochem_carbonflux_inst%som_c_leached_col , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total SOM C loss from vertical transport 

         totcolc                 =>    cnveg_carbonstate_inst%totc_col                    & ! Input:  [real(r8) (:) ]  (gC/m2) total column carbon, incl veg and cpool
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (.not. beg_vals_set(c)) then
            ! Skip the check if the beginning values weren't set for this column. This
            ! can happen, for example, if this is a newly-active column.
            cycle
         end if

         ! calculate the total column-level carbon storage, for mass conservation check
         col_endcb(c) = totcolc(c)

         ! calculate total column-level inputs
         col_cinputs = gpp(c)

         ! calculate total column-level outputs
         ! er = ar + hr, col_fire_closs includes patch-level fire losses
         col_coutputs = er(c) + col_fire_closs(c) + dwt_closs(c) + product_closs(c) + col_hrv_xsmrpool_to_atm(c)

         ! subtract leaching flux
         col_coutputs = col_coutputs - som_c_leached(c)

         ! calculate the total column-level carbon balance error for this time step
         col_errcb(c) = (col_cinputs - col_coutputs)*dt - (col_endcb(c) - col_begcb(c))

         ! check for significant errors
         if (abs(col_errcb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
         end if

      end do ! end of columns loop

      if (err_found) then
         c = err_index
         write(iulog,*)'column cbalance error = ', col_errcb(c), c
         write(iulog,*)'Latdeg,Londeg=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
         write(iulog,*)'begcb       = ',col_begcb(c)
         write(iulog,*)'endcb       = ',col_endcb(c)
         write(iulog,*)'delta store = ',col_endcb(c)-col_begcb(c)
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine CBalanceCheck

  !-----------------------------------------------------------------------
  subroutine NBalanceCheck(this, bounds, num_soilc, filter_soilc, &
       soilbiogeochem_nitrogenflux_inst, cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Perform nitrogen mass conservation check
    !
    ! !USES:
    use clm_varpar, only : crop_prog
    !
    ! !ARGUMENTS:
    class(cn_balance_type)                  , intent(inout) :: this
    type(bounds_type)                       , intent(in)    :: bounds          
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc                                      (:) ! filter for soil columns
    type(soilbiogeochem_nitrogenflux_type)  , intent(in)    :: soilbiogeochem_nitrogenflux_inst
    type(cnveg_nitrogenflux_type)           , intent(in)    :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,err_index,j  ! indices
    integer :: fc             ! lake filter indices
    logical :: err_found      ! error flag
    real(r8):: dt             ! radiation time step (seconds)
    real(r8):: col_ninputs(bounds%begc:bounds%endc) 
    real(r8):: col_noutputs(bounds%begc:bounds%endc) 
    real(r8):: col_errnb(bounds%begc:bounds%endc) 
    !-----------------------------------------------------------------------

    associate(                                                                             & 
         col_begnb           => this%begnb_col                                           , & ! Input:  [real(r8) (:) ]  (gN/m2) nitrogen mass, beginning of time step   
         col_endnb           => this%endnb_col                                           , & ! Output: [real(r8) (:) ]  (gN/m2) nitrogen mass, end of time step         
         beg_vals_set        => this%beg_vals_set_col                                    , & ! Input:  [logical  (:) ]  Whether begcb/begnb have been set in this time step
         ndep_to_sminn       => soilbiogeochem_nitrogenflux_inst%ndep_to_sminn_col       , & ! Input:  [real(r8) (:) ]  (gN/m2/s) atmospheric N deposition to soil mineral N        
         nfix_to_sminn       => soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col       , & ! Input:  [real(r8) (:) ]  (gN/m2/s) symbiotic/asymbiotic N fixation to soil mineral N 
         fert_to_sminn       => soilbiogeochem_nitrogenflux_inst%fert_to_sminn_col       , & ! Input:  [real(r8) (:) ]  (gN/m2/s)                                         
         soyfixn_to_sminn    => soilbiogeochem_nitrogenflux_inst%soyfixn_to_sminn_col    , & ! Input:  [real(r8) (:) ]  (gN/m2/s)                                         
         supplement_to_sminn => soilbiogeochem_nitrogenflux_inst%supplement_to_sminn_col , & ! Input:  [real(r8) (:) ]  (gN/m2/s) supplemental N supply                           
         denit               => soilbiogeochem_nitrogenflux_inst%denit_col               , & ! Input:  [real(r8) (:) ]  (gN/m2/s) total rate of denitrification           
         sminn_leached       => soilbiogeochem_nitrogenflux_inst%sminn_leached_col       , & ! Input:  [real(r8) (:) ]  (gN/m2/s) soil mineral N pool loss to leaching   
         smin_no3_leached    => soilbiogeochem_nitrogenflux_inst%smin_no3_leached_col    , & ! Input:  [real(r8) (:) ]  (gN/m2/s) soil mineral NO3 pool loss to leaching 
         smin_no3_runoff     => soilbiogeochem_nitrogenflux_inst%smin_no3_runoff_col     , & ! Input:  [real(r8) (:) ]  (gN/m2/s) soil mineral NO3 pool loss to runoff   
         f_n2o_nit           => soilbiogeochem_nitrogenflux_inst%f_n2o_nit_col           , & ! Input:  [real(r8) (:) ]  (gN/m2/s) flux of N2o from nitrification 
         som_n_leached       => soilbiogeochem_nitrogenflux_inst%som_n_leached_col       , & ! Input:  [real(r8) (:) ]  (gN/m2/s) total SOM N loss from vertical transport

         col_fire_nloss      => cnveg_nitrogenflux_inst%fire_nloss_col                   , & ! Input:  [real(r8) (:) ]  (gN/m2/s) total column-level fire N loss 
         dwt_nloss           => cnveg_nitrogenflux_inst%dwt_nloss_col                    , & ! Input:  [real(r8) (:) ]  (gN/m2/s) total nitrogen loss from product pools and conversion
         product_nloss       => cnveg_nitrogenflux_inst%product_nloss_col                , & ! Input:  [real(r8) (:) ]  (gN/m2/s) total wood product nitrogen loss

         totcoln             => cnveg_nitrogenstate_inst%totn_col                          & ! Input:  [real(r8) (:) ]  (gN/m2) total column nitrogen, incl veg 
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         if (.not. beg_vals_set(c)) then
            ! Skip the check if the beginning values weren't set for this column. This
            ! can happen, for example, if this is a newly-active column.
            cycle
         end if

         ! calculate the total column-level nitrogen storage, for mass conservation check
         col_endnb(c) = totcoln(c)

         ! calculate total column-level inputs
         col_ninputs(c) = ndep_to_sminn(c) + nfix_to_sminn(c) + supplement_to_sminn(c)
         if (crop_prog) then
            col_ninputs(c) = col_ninputs(c) + fert_to_sminn(c) + soyfixn_to_sminn(c)
         end if

         ! calculate total column-level outputs
         col_noutputs(c) = denit(c) + col_fire_nloss(c) + dwt_nloss(c) + product_nloss(c)

         if (.not. use_nitrif_denitrif) then
            col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
         else
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
         end if

         col_noutputs(c) = col_noutputs(c) - som_n_leached(c)

         ! calculate the total column-level nitrogen balance error for this time step
         col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - (col_endnb(c) - col_begnb(c))

         if (abs(col_errnb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
         end if

      end do ! end of columns loop

      if (err_found) then
         c = err_index
         write(iulog,*)'column nbalance error = ',col_errnb(c), c
         write(iulog,*)'Latdeg,Londeg         = ',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
         write(iulog,*)'begnb                 = ',col_begnb(c)
         write(iulog,*)'endnb                 = ',col_endnb(c)
         write(iulog,*)'delta store           = ',col_endnb(c)-col_begnb(c)
         write(iulog,*)'input mass            = ',col_ninputs(c)*dt
         write(iulog,*)'output mass           = ',col_noutputs(c)*dt
         write(iulog,*)'net flux              = ',(col_ninputs(c)-col_noutputs(c))*dt
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine NBalanceCheck

end module CNBalanceCheckMod
