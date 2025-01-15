module CNBalanceCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon/nitrogen mass balance checking.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_infnan_mod                  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use decompMod                       , only : bounds_type, subgrid_level_gridcell, subgrid_level_column
  use abortutils                      , only : endrun
  use clm_varctl                      , only : iulog, use_nitrif_denitrif, use_fates_bgc
  use clm_time_manager                , only : get_step_size_real
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenfluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemCarbonfluxType    , only : soilbiogeochem_carbonflux_type
  use CNProductsMod                   , only : cn_products_type
  use ColumnType                      , only : col                
  use GridcellType                    , only : grc
  use CNSharedParamsMod               , only : use_fun
  use CLMFatesInterfaceMod            , only : hlm_fates_interface_type
  use clm_varpar                      , only : nlevdecomp
  
  !
  implicit none
  private
  ! 
  ! !PUBLIC TYPES:
  type, public :: cn_balance_type
     private
     real(r8), pointer :: begcb_col(:)  ! (gC/m2) column carbon mass, beginning of time step
     real(r8), pointer :: endcb_col(:)  ! (gC/m2) column carbon mass, end of time step
     real(r8), pointer :: begnb_col(:)  ! (gN/m2) column nitrogen mass, beginning of time step
     real(r8), pointer :: endnb_col(:)  ! (gN/m2) column nitrogen mass, end of time step
     real(r8), pointer :: begcb_grc(:)  ! (gC/m2) gridcell carbon mass, beginning of time step
     real(r8), pointer :: endcb_grc(:)  ! (gC/m2) gridcell carbon mass, end of time step
     real(r8), pointer :: begnb_grc(:)  ! (gN/m2) gridcell nitrogen mass, beginning of time step
     real(r8), pointer :: endnb_grc(:)  ! (gN/m2) gridcell nitrogen mass, end of time step
     real(r8)          :: cwarning      ! (gC/m2) For a Carbon balance warning
     real(r8)          :: nwarning      ! (gN/m2) For a Nitrogen balance warning
     real(r8)          :: cerror        ! (gC/m2) For a Carbon balance error
     real(r8)          :: nerror        ! (gN/m2) For a Nitrogen balance error
   contains
     procedure , public  :: Init
     procedure , public  :: BeginCNGridcellBalance
     procedure , public  :: BeginCNColumnBalance
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

    ! Set warning and error levels for Carbon and Nitrogen balance
    ! These could become namelist items if we want them to change for different
    ! types of cases
    this%cwarning = 1.e-8_r8
    this%nwarning = 1.e-7_r8
    this%nerror   = 1.e-3_r8
    this%cerror   = 1.e-7_r8
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    class(cn_balance_type)         :: this
    type(bounds_type) , intent(in) :: bounds  

    integer :: begc, endc
    integer :: begg, endg

    begg = bounds%begg; endg = bounds%endg

    allocate(this%begcb_grc(begg:endg)) ; this%begcb_grc(:) = nan
    allocate(this%endcb_grc(begg:endg)) ; this%endcb_grc(:) = nan
    allocate(this%begnb_grc(begg:endg)) ; this%begnb_grc(:) = nan
    allocate(this%endnb_grc(begg:endg)) ; this%endnb_grc(:) = nan

    begc = bounds%begc; endc= bounds%endc

    allocate(this%begcb_col(begc:endc)) ; this%begcb_col(:) = nan
    allocate(this%endcb_col(begc:endc)) ; this%endcb_col(:) = nan
    allocate(this%begnb_col(begc:endc)) ; this%begnb_col(:) = nan
    allocate(this%endnb_col(begc:endc)) ; this%endnb_col(:) = nan
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine BeginCNGridcellBalance(this, bounds, cnveg_carbonflux_inst, &
       soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst, &
       c_products_inst, n_products_inst)
    !
    ! !DESCRIPTION:
    ! Calculate beginning gridcell-level carbon/nitrogen balance
    ! for mass conservation check
    !
    ! Should be called after CN state summaries have been computed
    ! and before the dynamic landunit area updates
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_balance_type)         , intent(inout)    :: this
    type(bounds_type)              , intent(in)       :: bounds
    type(soilbiogeochem_carbonstate_type), intent(in) :: soilbiogeochem_carbonstate_inst
    type(cnveg_carbonflux_type)    , intent(in)       :: cnveg_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(in) :: soilbiogeochem_nitrogenstate_inst
    type(cn_products_type)         , intent(in)       :: c_products_inst
    type(cn_products_type)         , intent(in)       :: n_products_inst
    !
    ! !LOCAL VARIABLES:
    integer :: g
    integer :: begg, endg
    real(r8) :: hrv_xsmrpool_amount_left_to_dribble(bounds%begg:bounds%endg)
    real(r8) :: gru_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg)
    real(r8) :: dwt_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg)
    !-----------------------------------------------------------------------

    associate(                                                &
         begcb          => this%begcb_grc                   , &  ! Output: [real(r8) (:)]  (gC/m2) gridcell carbon mass, beginning of time step
         begnb          => this%begnb_grc                   , &  ! Output: [real(r8) (:)]  (gN/m2) gridcell nitrogen mass, beginning of time step
         totc           => soilbiogeochem_carbonstate_inst%totc_grc  , &  ! Input:  [real(r8) (:)]  (gC/m2) total gridcell carbon, incl veg and cpool
         totn           => soilbiogeochem_nitrogenstate_inst%totn_grc, &  ! Input:  [real(r8) (:)]  (gN/m2) total gridcell nitrogen, incl veg
         c_cropprod1    => c_products_inst%cropprod1_grc    , &  ! Input:  [real(r8) (:)]  (gC/m2) carbon in crop products
         n_cropprod1    => n_products_inst%cropprod1_grc    , &  ! Input:  [real(r8) (:)]  (gC/m2) nitrogen in crop products
         c_tot_woodprod => c_products_inst%tot_woodprod_grc , &  ! Input:  [real(r8) (:)]  (gC/m2) total carbon in wood products
         n_tot_woodprod => n_products_inst%tot_woodprod_grc   &  ! Input:  [real(r8) (:)]  (gC/m2) total nitrogen in wood products
         )

    begg = bounds%begg; endg = bounds%endg
    
    if(.not.use_fates_bgc)then
       call cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_dribbler%get_amount_left_to_dribble_beg( &
         bounds, hrv_xsmrpool_amount_left_to_dribble(bounds%begg:bounds%endg))
       call cnveg_carbonflux_inst%dwt_conv_cflux_dribbler%get_amount_left_to_dribble_beg( &
            bounds, dwt_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg))
       call cnveg_carbonflux_inst%gru_conv_cflux_dribbler%get_amount_left_to_dribble_beg( &
            bounds, gru_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg))
    else
       hrv_xsmrpool_amount_left_to_dribble(bounds%begg:bounds%endg) = 0._r8
       dwt_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg) = 0._r8
       gru_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg) = 0._r8
    end if

    do g = begg, endg
       begcb(g) = totc(g) + c_tot_woodprod(g) + c_cropprod1(g) + &
                  hrv_xsmrpool_amount_left_to_dribble(g) + &
                  gru_conv_cflux_amount_left_to_dribble(g) + &
                  dwt_conv_cflux_amount_left_to_dribble(g)
       begnb(g) = totn(g) + n_tot_woodprod(g) + n_cropprod1(g)
    end do

    end associate

  end subroutine BeginCNGridcellBalance

  !-----------------------------------------------------------------------
  subroutine BeginCNColumnBalance(this, bounds, num_soilc, filter_soilc, &
       soilbiogeochem_carbonstate_inst,soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Calculate beginning column-level carbon/nitrogen balance, for mass conservation check
    !
    ! Should be called after CN state summaries have been recomputed for this time step
    ! (which should be after the dynamic landunit area updates and the associated filter
    ! updates - i.e., using the new version of the filters)
    !
    ! !ARGUMENTS:
    class(cn_balance_type)         , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds          
    integer                        , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_carbonstate_type), intent(in) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type), intent(in) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c
    !-----------------------------------------------------------------------

    associate(                                            & 
         col_begcb    => this%begcb_col                  , & ! Output: [real(r8) (:)]  (gC/m2) column carbon mass, beginning of time step
         col_begnb    => this%begnb_col                  , & ! Output: [real(r8) (:)]  (gN/m2) column nitrogen mass, beginning of time step
         totcolc      => soilbiogeochem_carbonstate_inst%totc_col , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
         totcoln      => soilbiogeochem_nitrogenstate_inst%totn_col & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
         )
    
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       col_begcb(c) = totcolc(c)
       col_begnb(c) = totcoln(c)

    end do

    end associate

  end subroutine BeginCNColumnBalance
 
  !-----------------------------------------------------------------------
  subroutine CBalanceCheck(this, bounds, num_soilc, filter_soilc, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst, c_products_inst, &
       clm_fates)
    !
    ! !USES:
    use subgridAveMod, only: c2g
    
    !
    ! !DESCRIPTION:
    ! Perform carbon mass conservation check for column and patch
    !
    ! Note on FATES: On fates colums, there is no vegetation biomass
    !                and no gpp flux. There is a litter input flux.
    
    !
    ! !ARGUMENTS:
    class(cn_balance_type)               , intent(inout) :: this
    type(bounds_type)                    , intent(in)    :: bounds          
    integer                              , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                              , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_carbonflux_type) , intent(in)    :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type), intent(inout) :: soilbiogeochem_carbonstate_inst
    type(cnveg_carbonflux_type)          , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: cnveg_carbonstate_inst
    type(cn_products_type)               , intent(in)    :: c_products_inst
    type(hlm_fates_interface_type)       , intent(inout) :: clm_fates
    
    !
    ! !LOCAL VARIABLES:
    integer :: c, g, err_index ! indices
    integer  :: s              ! fates site index (follows c)
    integer  :: fc             ! lake filter indices
    integer  :: ic             ! index of the current clump
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    real(r8) :: col_cinputs, grc_cinputs
    real(r8) :: col_coutputs, grc_coutputs
    real(r8) :: col_errcb(bounds%begc:bounds%endc) 
    real(r8) :: grc_errcb(bounds%begg:bounds%endg)
    real(r8) :: som_c_leached_grc(bounds%begg:bounds%endg)
    real(r8) :: hrv_xsmrpool_amount_left_to_dribble(bounds%begg:bounds%endg)
    real(r8) :: gru_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg)
    real(r8) :: dwt_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg)

    !-----------------------------------------------------------------------

    associate(                                                                            & 
         grc_begcb               =>    this%begcb_grc                                   , & ! Input:  [real(r8) (:) ]  (gC/m2) gridcell-level carbon mass, beginning of time step
         grc_endcb               =>    this%endcb_grc                                   , & ! Output: [real(r8) (:) ]  (gC/m2) gridcell-level carbon mass, end of time step
         totgrcc                 =>    soilbiogeochem_carbonstate_inst%totc_grc         , & ! Output:  [real(r8) (:)]  (gC/m2) total gridcell carbon, incl veg and cpool
         nbp_grc                 =>    cnveg_carbonflux_inst%nbp_grc                    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) net biome production (positive for sink)
         cropprod1_grc           =>    c_products_inst%cropprod1_grc                    , & ! Input:  [real(r8) (:)]  (gC/m2) carbon in crop products
         tot_woodprod_grc        =>    c_products_inst%tot_woodprod_grc                 , & ! Input:  [real(r8) (:)]  (gC/m2) total carbon in wood products
         dwt_seedc_to_leaf_grc   =>    cnveg_carbonflux_inst%dwt_seedc_to_leaf_grc      , & ! Input:  [real(r8) (:)]  (gC/m2/s) seed source sent to leaf
         dwt_seedc_to_deadstem_grc =>  cnveg_carbonflux_inst%dwt_seedc_to_deadstem_grc  , & ! Input:  [real(r8) (:)]  (gC/m2/s) seed source sent to deadstem
         col_begcb               =>    this%begcb_col                                   , & ! Input:  [real(r8) (:) ]  (gC/m2) carbon mass, beginning of time step 
         col_endcb               =>    this%endcb_col                                   , & ! Output: [real(r8) (:) ]  (gC/m2) carbon mass, end of time step 
         wood_harvestc           =>    cnveg_carbonflux_inst%wood_harvestc_col          , & ! Input:  [real(r8) (:) ]  (gC/m2/s) wood harvest (to product pools)
         gru_conv_cflux          =>    cnveg_carbonflux_inst%gru_conv_cflux_col         , & ! Input:  [real(r8) (:) ]  (gC/m2/s) wood harvest (to product pools)
         gru_wood_productc_gain  =>    cnveg_carbonflux_inst%gru_wood_productc_gain_col , & ! Input:  [real(r8) (:) ]  (gC/m2/s) wood harvest (to product pools)
         crop_harvestc_to_cropprodc     =>    cnveg_carbonflux_inst%crop_harvestc_to_cropprodc_col    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) crop harvest C to 1-year crop product pool
         gpp                     =>    cnveg_carbonflux_inst%gpp_col                    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) gross primary production
         er                      =>    cnveg_carbonflux_inst%er_col                     , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs          =>    cnveg_carbonflux_inst%fire_closs_col             , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total column-level fire C loss
         col_hrv_xsmrpool_to_atm =>    cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_col    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool harvest mortality 
         col_xsmrpool_to_atm     =>    cnveg_carbonflux_inst%xsmrpool_to_atm_col         , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool crop harvest loss to atm
         som_c_leached           =>    soilbiogeochem_carbonflux_inst%som_c_leached_col , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total SOM C loss from vertical transport 

         totcolc                 =>    soilbiogeochem_carbonstate_inst%totc_col          , & ! Input:  [real(r8) (:) ]  (gC/m2) total column carbon, incl veg and cpool
         fates_litter_flux       =>    soilbiogeochem_carbonflux_inst%fates_litter_flux  &   ! Total carbon litter flux from FATES to CLM [gC/m2/s]
         )

      ! set time steps
      dt = get_step_size_real()

      ! clump index
      ic = bounds%clump_index
      
      err_found = .false.
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! calculate the total column-level carbon storage, for mass conservation check
         ! for bigleaf, totcolc includes soil and all of the veg c pools including cpool, xfer, etc
         ! for fates, totcolc only includes soil and non-fates litter carbon,
         ! see soibiogeochem_carbonstate_inst%summary for calculations
         col_endcb(c) = totcolc(c)
         
         
         if( col%is_fates(c) ) then

            s = clm_fates%f2hmap(ic)%hsites(c)
            
            col_cinputs = fates_litter_flux(c)
            
            ! calculate total column-level outputs
            ! fates has already exported burn losses and fluxes to the atm
            ! So they are irrelevant here
            ! (gC/m2/s) total heterotrophic respiration
            col_coutputs = soilbiogeochem_carbonflux_inst%hr_col(c)

         else
            
            ! calculate total column-level inputs
            col_cinputs = gpp(c)
            
            ! calculate total column-level outputs
            ! er = ar + hr, col_fire_closs includes patch-level fire losses
            col_coutputs = er(c) + col_fire_closs(c) + col_hrv_xsmrpool_to_atm(c) + &
                 col_xsmrpool_to_atm(c) + gru_conv_cflux(c)
         
            ! Fluxes to product pools are included in column-level outputs: the product
            ! pools are not included in totcolc, so are outside the system with respect to
            ! these balance checks. (However, the dwt flux to product pools is NOT included,
            ! since col_begcb is initialized after the dynamic area adjustments - i.e.,
            ! after the dwt term has already been taken out.)
            col_coutputs = col_coutputs + &
                 wood_harvestc(c) + &
                 gru_wood_productc_gain(c) + &
                 crop_harvestc_to_cropprodc(c)
         
         end if

         ! subtract leaching flux
         col_coutputs = col_coutputs - som_c_leached(c)

         ! calculate the total column-level carbon balance error for this time step
         col_errcb(c) = (col_cinputs - col_coutputs)*dt - &
              (col_endcb(c) - col_begcb(c))

         ! check for significant errors
         if (abs(col_errcb(c)) > this%cerror) then
            err_found = .true.
            err_index = c
         end if
          if (abs(col_errcb(c)) > this%cwarning) then
            write(iulog,*) 'cbalance warning at c =', c, col_errcb(c), col_endcb(c)
         end if

      end do ! end of columns loop

      if (err_found) then
         c = err_index
         write(iulog,*)'column cbalance error    = ', col_errcb(c), c
         write(iulog,*)'is fates column?         = ', col%is_fates(c)
         write(iulog,*)'Latdeg,Londeg=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
         write(iulog,*)'begcb                    = ',col_begcb(c)
         write(iulog,*)'endcb                    = ',col_endcb(c)
         write(iulog,*)'delta store              = ',col_endcb(c)-col_begcb(c)
         write(iulog,*)'--- Inputs ---'
         if( col%is_fates(c) ) then
            write(iulog,*)'fates litter_flux        = ',fates_litter_flux(c)*dt
         else
            write(iulog,*)'gpp                      = ',gpp(c)*dt
         end if
         write(iulog,*)'--- Outputs ---'
         if( .not.col%is_fates(c) ) then
            write(iulog,*)'er                       = ',er(c)*dt
            write(iulog,*)'col_fire_closs           = ',col_fire_closs(c)*dt
            write(iulog,*)'col_hrv_xsmrpool_to_atm  = ',col_hrv_xsmrpool_to_atm(c)*dt
            write(iulog,*)'col_xsmrpool_to_atm      = ',col_xsmrpool_to_atm(c)*dt
            write(iulog,*)'wood_harvestc            = ',wood_harvestc(c)*dt
            write(iulog,*)'crop_harvestc_to_cropprodc = ', crop_harvestc_to_cropprodc(c)*dt
         else
            write(iulog,*)'hr                       = ',soilbiogeochem_carbonflux_inst%hr_col(c)*dt
         end if
         write(iulog,*)'-1*som_c_leached         = ',som_c_leached(c)*dt
         call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg=errMsg(sourcefile, __LINE__))
      end if

      ! Repeat error check at the gridcell level
      call c2g( bounds = bounds, &
         carr = totcolc(bounds%begc:bounds%endc), &
         garr = totgrcc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')
      call c2g( bounds = bounds, &
         carr = som_c_leached(bounds%begc:bounds%endc), &
         garr = som_c_leached_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

      err_found = .false.
      do g = bounds%begg, bounds%endg
         ! calculate gridcell-level carbon storage for mass conservation check
         ! Notes:
         ! totgrcc = totcolc = totc_p2c_col(c) + soilbiogeochem_cwdc_col(c) + soilbiogeochem_totlitc_col(c) + soilbiogeochem_totmicc_col(c) + soilbiogeochem_totsomc_col(c) + soilbiogeochem_ctrunc_col(c)
         ! totc_p2c_col = totc_patch = totvegc_patch(p) + xsmrpool_patch(p) + ctrunc_patch(p) + cropseedc_deficit_patch(p)
         ! Not including seedc_grc in grc_begcb and grc_endcb because
         ! seedc_grc forms out of thin air, for now, and equals
         ! -1 * (dwt_seedc_to_leaf_grc(g) + dwt_seedc_to_deadstem_grc(g))
         ! We account for the latter fluxes as inputs below; the same
         ! fluxes have entered the pools earlier in the timestep. For true
         ! conservation we would need to add a flux out of npp into seed.

         if(.not.use_fates_bgc)then
            call cnveg_carbonflux_inst%hrv_xsmrpool_to_atm_dribbler%get_amount_left_to_dribble_end( &
                 bounds, hrv_xsmrpool_amount_left_to_dribble(bounds%begg:bounds%endg))
            call cnveg_carbonflux_inst%dwt_conv_cflux_dribbler%get_amount_left_to_dribble_end( &
                 bounds, dwt_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg))
            call cnveg_carbonflux_inst%gru_conv_cflux_dribbler%get_amount_left_to_dribble_end( &
                 bounds, gru_conv_cflux_amount_left_to_dribble(bounds%begg:bounds%endg))

            grc_endcb(g) = totgrcc(g) + tot_woodprod_grc(g) + cropprod1_grc(g) + &
                 hrv_xsmrpool_amount_left_to_dribble(g) + &
                 gru_conv_cflux_amount_left_to_dribble(g) + &
                 dwt_conv_cflux_amount_left_to_dribble(g)
            
            ! calculate total gridcell-level inputs
            ! slevis notes:
            ! nbp_grc = nep_grc - fire_closs_grc - hrv_xsmrpool_to_atm_dribbled_grc - &
            !           dwt_conv_cflux_dribbled_grc - gru_conv_cflux_dribbled_grc - product_closs_grc
            grc_cinputs = nbp_grc(g) + & 
                 dwt_seedc_to_leaf_grc(g) + dwt_seedc_to_deadstem_grc(g)
            
            ! calculate total gridcell-level outputs
            grc_coutputs = - som_c_leached_grc(g)
            
            ! calculate the total gridcell-level carbon balance error
            ! for this time step
            grc_errcb(g) = (grc_cinputs - grc_coutputs) * dt - &
                 (grc_endcb(g) - grc_begcb(g))
            
         else
            
            ! Totally punt on this for now. We just don't track these gridscale variables yet (RGK)
            grc_cinputs  = 0._r8
            grc_endcb(g) = grc_begcb(g)
            grc_coutputs = 0._r8
            grc_errcb(g) = 0._r8
            
         end if
         
         ! check for significant errors
         if (abs(grc_errcb(g)) > this%cerror) then
            err_found = .true.
            err_index = g
         end if
         if (abs(grc_errcb(g)) > this%cwarning) then
            write(iulog,*) 'cbalance warning at g =', g, grc_errcb(g), grc_endcb(g)
         end if
      end do ! end of gridcell loop

      if (err_found) then
         g = err_index
         write(iulog,*)'gridcell cbalance error =', grc_errcb(g), g
         write(iulog,*)'latdeg, londeg          =', grc%latdeg(g), grc%londeg(g)
         write(iulog,*)'begcb                   =', grc_begcb(g)
         write(iulog,*)'endcb                   =', grc_endcb(g)
         write(iulog,*)'delta store             =', grc_endcb(g) - grc_begcb(g)
         write(iulog,*)'--- Inputs ---'
         write(iulog,*)'nbp_grc                 =', nbp_grc(g) * dt
         write(iulog,*)'dwt_seedc_to_leaf_grc   =', dwt_seedc_to_leaf_grc(g) * dt
         write(iulog,*)'dwt_seedc_to_deadstem_grc =', dwt_seedc_to_deadstem_grc(g) * dt
         write(iulog,*)'--- Outputs ---'
         write(iulog,*)'-1*som_c_leached_grc    = ', som_c_leached_grc(g) * dt
         call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, msg=errMsg(sourcefile, __LINE__))
      end if

    end associate

  end subroutine CBalanceCheck

  !-----------------------------------------------------------------------
  subroutine NBalanceCheck(this, bounds, num_soilc, filter_soilc, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst, &
       cnveg_nitrogenflux_inst, &
       cnveg_nitrogenstate_inst, n_products_inst, atm2lnd_inst, clm_fates)
    !
    ! !DESCRIPTION:
    ! Perform nitrogen mass conservation check
    !
    ! !USES:
    use clm_varctl, only : use_crop
    use subgridAveMod, only: c2g
    use atm2lndType, only: atm2lnd_type
    !
    ! !ARGUMENTS:
    class(cn_balance_type)                  , intent(inout) :: this
    type(bounds_type)                       , intent(in)    :: bounds          
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc                                      (:) ! filter for soil columns
    type(soilbiogeochem_nitrogenflux_type)  , intent(in)    :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)           , intent(in)    :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(cn_products_type)                  , intent(in)    :: n_products_inst
    type(atm2lnd_type)                      , intent(in)    :: atm2lnd_inst
    type(hlm_fates_interface_type)          , intent(inout) :: clm_fates
    
    !
    ! !LOCAL VARIABLES:
    integer :: c,err_index,j,s ! indices
    integer :: ic              ! index of clump
    integer :: g               ! gridcell index
    integer :: fc              ! lake filter indices
    logical :: err_found       ! error flag
    real(r8):: dt              ! radiation time step (seconds)
    real(r8):: col_ninputs(bounds%begc:bounds%endc) 
    real(r8):: col_noutputs(bounds%begc:bounds%endc) 
    real(r8):: col_errnb(bounds%begc:bounds%endc) 
    real(r8):: col_ninputs_partial(bounds%begc:bounds%endc)
    real(r8):: col_noutputs_partial(bounds%begc:bounds%endc)
    real(r8):: grc_ninputs_partial(bounds%begg:bounds%endg)
    real(r8):: grc_noutputs_partial(bounds%begg:bounds%endg)
    real(r8):: grc_ninputs(bounds%begg:bounds%endg)
    real(r8):: grc_noutputs(bounds%begg:bounds%endg)
    real(r8):: grc_errnb(bounds%begg:bounds%endg)
    !-----------------------------------------------------------------------

    associate(                                                                             & 
         grc_begnb           => this%begnb_grc                                           , & ! Input:  [real(r8) (:) ]  (gN/m2) gridcell nitrogen mass, beginning of time step
         grc_endnb           => this%endnb_grc                                           , & ! Output: [real(r8) (:) ]  (gN/m2) gridcell nitrogen mass, end of time step
         totgrcn             => soilbiogeochem_nitrogenstate_inst%totn_grc                        , & ! Input:  [real(r8) (:) ]  (gN/m2) total gridcell nitrogen, incl veg
         cropprod1_grc       => n_products_inst%cropprod1_grc                            , & ! Input:  [real(r8) (:)]  (gN/m2) nitrogen in crop products
         product_loss_grc    => n_products_inst%product_loss_grc                         , & ! Input:  [real(r8) (:)]  (gN/m2) losses from wood & crop products
         tot_woodprod_grc    => n_products_inst%tot_woodprod_grc                         , & ! Input:  [real(r8) (:)]  (gN/m2) total nitrogen in wood products
         dwt_seedn_to_leaf_grc   =>    cnveg_nitrogenflux_inst%dwt_seedn_to_leaf_grc     , & ! Input:  [real(r8) (:)]  (gN/m2/s) seed source sent to leaf
         dwt_seedn_to_deadstem_grc =>  cnveg_nitrogenflux_inst%dwt_seedn_to_deadstem_grc , & ! Input:  [real(r8) (:)]  (gN/m2/s) seed source sent to deadstem
         dwt_conv_nflux_grc  =>  cnveg_nitrogenflux_inst%dwt_conv_nflux_grc              , & ! Input:  [real(r8) (:)]  (gN/m2/s) dwt_conv_nflux_patch summed to the gridcell-level
         col_begnb           => this%begnb_col                                           , & ! Input:  [real(r8) (:) ]  (gN/m2) column nitrogen mass, beginning of time step
         col_endnb           => this%endnb_col                                           , & ! Output: [real(r8) (:) ]  (gN/m2) column nitrogen mass, end of time step
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
         gru_conv_nflux_grc  => cnveg_nitrogenflux_inst%gru_conv_nflux_grc               , & ! Input:  [real(r8) (:) ]  (gC/m2/s) wood harvest (to product pools) summed to the gridcell level
         gru_conv_nflux      => cnveg_nitrogenflux_inst%gru_conv_nflux_col               , & ! Input:  [real(r8) (:) ]  (gC/m2/s) wood harvest (to product pools)
         gru_wood_productn_gain => cnveg_nitrogenflux_inst%gru_wood_productn_gain_col    , & ! Input:  [real(r8) (:) ]  (gC/m2/s) wood harvest (to product pools)
         gru_wood_productn_gain_grc => cnveg_nitrogenflux_inst%gru_wood_productn_gain_grc, & ! Input:  [real(r8) (:) ]  (gC/m2/s) wood harvest (to product pools) summed to the gridcell level
         crop_harvestn_to_cropprodn => cnveg_nitrogenflux_inst%crop_harvestn_to_cropprodn_col          , & ! Input:  [real(r8) (:) ]  (gN/m2/s) crop harvest N to 1-year crop product pool

         totcoln             => soilbiogeochem_nitrogenstate_inst%totn_col               , & ! Input:  [real(r8) (:) ]  (gN/m2) total column nitrogen, incl veg
         sminn_to_plant      => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_col,       &
         fates_litter_flux   => soilbiogeochem_nitrogenflux_inst%fates_litter_flux  &   ! Total nitrogen litter flux from FATES to CLM [gN/m2/s]
         )


      ! set time steps
      dt = get_step_size_real()

      ! initialize local arrays
      col_ninputs_partial(:) = 0._r8
      col_noutputs_partial(:) = 0._r8

      ! clump index
      ic = bounds%clump_index
      
      err_found = .false.
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         ! calculate the total column-level nitrogen storage, for mass conservation check
         col_endnb(c) = totcoln(c)

         ! calculate total column-level inputs
         col_ninputs(c) = ndep_to_sminn(c) + nfix_to_sminn(c) + supplement_to_sminn(c)

         ! If using fates, pass in the decomposition flux
         if( col%is_fates(c) ) then
            col_ninputs(c) = col_ninputs(c)  + fates_litter_flux(c)
         end if
         
         if(use_fun)then
            col_ninputs(c) = col_ninputs(c) + ffix_to_sminn(c) ! for FUN, free living fixation is a seprate flux. RF. 
         endif
     
         if (use_crop) then
            col_ninputs(c) = col_ninputs(c) + fert_to_sminn(c) + soyfixn_to_sminn(c)
         end if

         col_ninputs_partial(c) = col_ninputs(c)
         
         ! calculate total column-level outputs

         col_noutputs(c) = denit(c)

         if( .not.col%is_fates(c) ) then
            
            col_noutputs(c) = col_noutputs(c) + col_fire_nloss(c) + gru_conv_nflux(c)

            ! Fluxes to product pools are included in column-level outputs: the product
            ! pools are not included in totcoln, so are outside the system with respect to
            ! these balance checks. (However, the dwt flux to product pools is NOT included,
            ! since col_begnb is initialized after the dynamic area adjustments - i.e.,
            ! after the dwt term has already been taken out.)
            col_noutputs(c) = col_noutputs(c) + &
                 wood_harvestn(c) + &
                 gru_wood_productn_gain(c) + &
                 crop_harvestn_to_cropprodn(c)
         else
            
            ! If we are using fates, remove plant uptake
            col_noutputs(c) = col_noutputs(c) +  sminn_to_plant(c)
            
         end if

         if (.not. use_nitrif_denitrif) then
            col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
         else
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
         end if

         col_noutputs(c) = col_noutputs(c) - som_n_leached(c)
         
         col_noutputs_partial(c) = col_noutputs(c)

         if( .not.col%is_fates(c) ) then
            col_noutputs_partial(c) = col_noutputs_partial(c) - &
                 wood_harvestn(c) - &
                 crop_harvestn_to_cropprodn(c)
         end if
         
         ! calculate the total column-level nitrogen balance error for this time step
         col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - &
              (col_endnb(c) - col_begnb(c))

         if (abs(col_errnb(c)) > this%nerror) then
            err_found = .true.
            err_index = c
         end if
         
         if (abs(col_errnb(c)) > this%nwarning) then
            write(iulog,*) 'nbalance warning at c =', c, col_errnb(c), col_endnb(c)
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
         if(col%is_fates(c))then
            write(iulog,*)'inputs,ndep,nfix,suppn= ',ndep_to_sminn(c)*dt,nfix_to_sminn(c)*dt,supplement_to_sminn(c)*dt
         else
            write(iulog,*)'inputs,ffix,nfix,ndep = ',ffix_to_sminn(c)*dt,nfix_to_sminn(c)*dt,ndep_to_sminn(c)*dt
         end if
         if(col%is_fates(c))then
            write(iulog,*)'outputs,lch,roff,dnit,plnt = ',smin_no3_leached(c)*dt, smin_no3_runoff(c)*dt,f_n2o_nit(c)*dt,sminn_to_plant(c)*dt
         else
            write(iulog,*)'outputs,lch,roff,dnit    = ',smin_no3_leached(c)*dt, smin_no3_runoff(c)*dt,f_n2o_nit(c)*dt
         end if
         call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg=errMsg(sourcefile, __LINE__))
      end if

      if_notfates: if(.not.use_fates_bgc)then

         ! Repeat error check at the gridcell level
         call c2g( bounds = bounds, &
              carr = totcoln(bounds%begc:bounds%endc), &
              garr = totgrcn(bounds%begg:bounds%endg), &
              c2l_scale_type = 'unity', &
              l2g_scale_type = 'unity')
         call c2g( bounds = bounds, &
              carr = col_ninputs_partial(bounds%begc:bounds%endc), &
              garr = grc_ninputs_partial(bounds%begg:bounds%endg), &
              c2l_scale_type = 'unity', &
              l2g_scale_type = 'unity')
         call c2g( bounds = bounds, &
              carr = col_noutputs_partial(bounds%begc:bounds%endc), &
              garr = grc_noutputs_partial(bounds%begg:bounds%endg), &
              c2l_scale_type = 'unity', &
              l2g_scale_type = 'unity')

         err_found = .false.
         do g = bounds%begg, bounds%endg
            ! calculate the total gridcell-level nitrogen storage, for mass conservation check
            ! Notes:
            ! Not including seedn_grc in grc_begnb and grc_endnb because
            ! seedn_grc forms out of thin air, for now, and equals
            ! -1 * (dwt_seedn_to_leaf_grc(g) + dwt_seedn_to_deadstem_grc(g))
            ! We account for the latter fluxes as inputs below; the same
            ! fluxes have entered the pools earlier in the timestep. For true
            ! conservation we would need to add a flux out of nfix into seed.
            grc_endnb(g) = totgrcn(g) + tot_woodprod_grc(g) + cropprod1_grc(g)

            ! calculate total gridcell-level inputs
            grc_ninputs(g) = grc_ninputs_partial(g) + &
                 dwt_seedn_to_leaf_grc(g) + &
                 dwt_seedn_to_deadstem_grc(g)
            
            ! calculate total gridcell-level outputs
            grc_noutputs(g) = grc_noutputs_partial(g) + &
                 dwt_conv_nflux_grc(g) + &
                 product_loss_grc(g) - &
                 ! Subtract the next one because it is present in
                 ! grc_noutputs_partial but not needed at the
                 ! gridcell level
                 gru_wood_productn_gain_grc(g)

            ! calculate the total gridcell-level nitrogen balance error for this time step
            grc_errnb(g) = (grc_ninputs(g) - grc_noutputs(g)) * dt - &
                 (grc_endnb(g) - grc_begnb(g))
      
            if (abs(grc_errnb(g)) > this%nerror) then
               err_found = .true.
               err_index = g
            end if
            
            if (abs(grc_errnb(g)) > this%nwarning) then
               write(iulog,*) 'nbalance warning at g =', g, grc_errnb(g), grc_endnb(g)
            end if
         end do
         if (err_found) then
            g = err_index
            write(iulog,*) 'gridcell nbalance error  =', grc_errnb(g), g
            write(iulog,*) 'latdeg, londeg           =', grc%latdeg(g), grc%londeg(g)
            write(iulog,*) 'begnb                    =', grc_begnb(g)
            write(iulog,*) 'endnb                    =', grc_endnb(g)
            write(iulog,*) 'delta store              =', grc_endnb(g) - grc_begnb(g)
            write(iulog,*) 'input mass               =', grc_ninputs(g) * dt
            write(iulog,*) 'output mass              =', grc_noutputs(g) * dt
            write(iulog,*) 'net flux                 =', (grc_ninputs(g) - grc_noutputs(g)) * dt
            write(iulog,*) '--- Inputs ---'
            write(iulog,*) 'grc_ninputs_partial      =', grc_ninputs_partial(g) * dt
            write(iulog,*) 'dwt_seedn_to_leaf_grc    =', dwt_seedn_to_leaf_grc(g) * dt
            write(iulog,*) 'dwt_seedn_to_deadstem_grc =', dwt_seedn_to_deadstem_grc(g) * dt
            write(iulog,*) '--- Outputs ---'
            write(iulog,*) 'grc_noutputs_partial     =', grc_noutputs_partial(g) * dt
            write(iulog,*) 'dwt_conv_nflux_grc       =', dwt_conv_nflux_grc(g) * dt
            write(iulog,*) '-gru_wood_productn_gain_grc =', -gru_wood_productn_gain_grc(g) * dt
            write(iulog,*) 'product_loss_grc         =', product_loss_grc(g) * dt
            call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, msg=errMsg(sourcefile, __LINE__))
         end if
         
      end if if_notfates

    end associate
    
  end subroutine NBalanceCheck

end module CNBalanceCheckMod
