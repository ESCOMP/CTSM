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
  use CNSharedParamsMod               , only : use_fun

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
   contains
     procedure , public  :: Init
     procedure , public  :: BeginCNBalance
     procedure , public  :: CBalanceCheck
     procedure , public  :: NBalanceCheck
     procedure , private :: InitAllocate 
  end type cn_balance_type
  !

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
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
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine BeginCNBalance(this, bounds, num_soilc, filter_soilc, &
       cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Calculate beginning column-level carbon/nitrogen balance, for mass conservation check
    !
    ! Should be called after the CN state summaries have been computed for this time step
    ! (which should be after the dynamic landunit area updates and the associated filter
    ! updates - i.e., using the new version of the filters)
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
         totcolc      => cnveg_carbonstate_inst%totc_col , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
         totcoln      => cnveg_nitrogenstate_inst%totn_col & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
         )
    
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       col_begcb(c) = totcolc(c)
       col_begnb(c) = totcoln(c)
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
         wood_harvestc           =>    cnveg_carbonflux_inst%wood_harvestc_col          , & ! Input:  [real(r8) (:) ]  (gC/m2/s) wood harvest (to product pools)
         grainc_to_cropprodc     =>    cnveg_carbonflux_inst%grainc_to_cropprodc_col    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) grain C to 1-year crop product pool
         gpp                     =>    cnveg_carbonflux_inst%gpp_col                    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) gross primary production
         er                      =>    cnveg_carbonflux_inst%er_col                     , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs          =>    cnveg_carbonflux_inst%fire_closs_col             , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total column-level fire C loss
         col_hrv_xsmrpool_to_atm =>    cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_col    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool harvest mortality 
         col_xsmrpool_to_atm     =>   cnveg_carbonflux_inst%xsmrpool_to_atm_col         , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool crop harvest loss to atm
         som_c_leached           =>    soilbiogeochem_carbonflux_inst%som_c_leached_col , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total SOM C loss from vertical transport 

         totcolc                 =>    cnveg_carbonstate_inst%totc_col                    & ! Input:  [real(r8) (:) ]  (gC/m2) total column carbon, incl veg and cpool
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! calculate the total column-level carbon storage, for mass conservation check
         col_endcb(c) = totcolc(c)

         ! calculate total column-level inputs
         col_cinputs = gpp(c)

         ! calculate total column-level outputs
         ! er = ar + hr, col_fire_closs includes patch-level fire losses
         col_coutputs = er(c) + col_fire_closs(c) + col_hrv_xsmrpool_to_atm(c) + &
              col_xsmrpool_to_atm(c)

         ! Fluxes to product pools are included in column-level outputs: the product
         ! pools are not included in totcolc, so are outside the system with respect to
         ! these balance checks. (However, the dwt flux to product pools is NOT included,
         ! since col_begcb is initialized after the dynamic area adjustments - i.e.,
         ! after the dwt term has already been taken out.)
         col_coutputs = col_coutputs + &
              wood_harvestc(c) + &
              grainc_to_cropprodc(c)

         ! subtract leaching flux
         col_coutputs = col_coutputs - som_c_leached(c)

         ! calculate the total column-level carbon balance error for this time step
         col_errcb(c) = (col_cinputs - col_coutputs)*dt - &
              (col_endcb(c) - col_begcb(c))

         ! check for significant errors
         if (abs(col_errcb(c)) > 1e-7_r8) then
            err_found = .true.
            err_index = c
         end if
          if (abs(col_errcb(c)) > 1e-8_r8) then
            write(iulog,*) 'cbalance warning',c,col_errcb(c),col_endcb(c)
         end if



      end do ! end of columns loop

      if (err_found) then
         c = err_index
         write(iulog,*)'column cbalance error    = ', col_errcb(c), c
         write(iulog,*)'Latdeg,Londeg=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
         write(iulog,*)'begcb                    = ',col_begcb(c)
         write(iulog,*)'endcb                    = ',col_endcb(c)
         write(iulog,*)'delta store              = ',col_endcb(c)-col_begcb(c)
         write(iulog,*)'--- Inputs ---'
         write(iulog,*)'gpp                      = ',gpp(c)*dt
         write(iulog,*)'--- Outputs ---'
         write(iulog,*)'er                       = ',er(c)*dt
         write(iulog,*)'col_fire_closs           = ',col_fire_closs(c)*dt
         write(iulog,*)'col_hrv_xsmrpool_to_atm  = ',col_hrv_xsmrpool_to_atm(c)*dt
         write(iulog,*)'col_xsmrpool_to_atm      = ',col_xsmrpool_to_atm(c)*dt
         write(iulog,*)'wood_harvestc            = ',wood_harvestc(c)*dt
         write(iulog,*)'grainc_to_cropprodc      = ',grainc_to_cropprodc(c)*dt
         write(iulog,*)'-1*som_c_leached         = ',som_c_leached(c)*dt
         call endrun(msg=errMsg(sourcefile, __LINE__))
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
    use clm_varctl, only : use_crop
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
         ndep_to_sminn       => soilbiogeochem_nitrogenflux_inst%ndep_to_sminn_col       , & ! Input:  [real(r8) (:) ]  (gN/m2/s) atmospheric N deposition to soil mineral N        
         nfix_to_sminn       => soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col       , & ! Input:  [real(r8) (:) ]  (gN/m2/s) symbiotic/asymbiotic N fixation to soil mineral N 
         ffix_to_sminn       => soilbiogeochem_nitrogenflux_inst%ffix_to_sminn_col       , & ! Input:  [real(r8) (:) ]  (gN/m2/s) free living N fixation to soil mineral N         
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
         wood_harvestn       => cnveg_nitrogenflux_inst%wood_harvestn_col                , & ! Input:  [real(r8) (:) ]  (gN/m2/s) wood harvest (to product pools)
         grainn_to_cropprodn => cnveg_nitrogenflux_inst%grainn_to_cropprodn_col          , & ! Input:  [real(r8) (:) ]  (gN/m2/s) grain N to 1-year crop product pool

         totcoln             => cnveg_nitrogenstate_inst%totn_col                          & ! Input:  [real(r8) (:) ]  (gN/m2) total column nitrogen, incl veg 
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         ! calculate the total column-level nitrogen storage, for mass conservation check
         col_endnb(c) = totcoln(c)

         ! calculate total column-level inputs
         col_ninputs(c) = ndep_to_sminn(c) + nfix_to_sminn(c) + supplement_to_sminn(c)
         
         if(use_fun)then
            col_ninputs(c) = col_ninputs(c) + ffix_to_sminn(c) ! for FUN, free living fixation is a seprate flux. RF. 
         endif
     
         if (use_crop) then
            col_ninputs(c) = col_ninputs(c) + fert_to_sminn(c) + soyfixn_to_sminn(c)
         end if

         ! calculate total column-level outputs
         col_noutputs(c) = denit(c) + col_fire_nloss(c)

         ! Fluxes to product pools are included in column-level outputs: the product
         ! pools are not included in totcoln, so are outside the system with respect to
         ! these balance checks. (However, the dwt flux to product pools is NOT included,
         ! since col_begnb is initialized after the dynamic area adjustments - i.e.,
         ! after the dwt term has already been taken out.)
         col_noutputs(c) = col_noutputs(c) + &
              wood_harvestn(c) + &
              grainn_to_cropprodn(c)

         if (.not. use_nitrif_denitrif) then
            col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
         else
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
         end if

         col_noutputs(c) = col_noutputs(c) - som_n_leached(c)

         ! calculate the total column-level nitrogen balance error for this time step
         col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - &
              (col_endnb(c) - col_begnb(c))

         if (abs(col_errnb(c)) > 1e-3_r8) then
            err_found = .true.
            err_index = c
         end if
         
         if (abs(col_errnb(c)) > 1e-7_r8) then
            write(iulog,*) 'nbalance warning',c,col_errnb(c),col_endnb(c)
            write(iulog,*)'inputs,ffix,nfix,ndep = ',ffix_to_sminn(c)*dt,nfix_to_sminn(c)*dt,ndep_to_sminn(c)*dt
            write(iulog,*)'outputs,lch,roff,dnit = ',smin_no3_leached(c)*dt, smin_no3_runoff(c)*dt,f_n2o_nit(c)*dt
         end if

      end do ! end of columns loop

      if (err_found) then
         c = err_index
         write(iulog,*)'column nbalance error    = ',col_errnb(c), c
         write(iulog,*)'Latdeg,Londeg            = ',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
         write(iulog,*)'begnb                    = ',col_begnb(c)
         write(iulog,*)'endnb                    = ',col_endnb(c)
         write(iulog,*)'delta store              = ',col_endnb(c)-col_begnb(c)
         write(iulog,*)'input mass               = ',col_ninputs(c)*dt
         write(iulog,*)'output mass              = ',col_noutputs(c)*dt
         write(iulog,*)'net flux                 = ',(col_ninputs(c)-col_noutputs(c))*dt
         write(iulog,*)'inputs,ffix,nfix,ndep    = ',ffix_to_sminn(c)*dt,nfix_to_sminn(c)*dt,ndep_to_sminn(c)*dt
         write(iulog,*)'outputs,ffix,nfix,ndep   = ',smin_no3_leached(c)*dt, smin_no3_runoff(c)*dt,f_n2o_nit(c)*dt
        
         
         
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

    end associate

  end subroutine NBalanceCheck

end module CNBalanceCheckMod
