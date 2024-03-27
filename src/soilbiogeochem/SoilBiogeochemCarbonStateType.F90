module SoilBiogeochemCarbonStateType

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use decompMod                          , only : bounds_type
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp, nlevsoi
  use clm_varcon                         , only : spval, ispval, dzsoi_decomp, zisoi, zsoi, c3_r2
  use clm_varctl                         , only : iulog, spinup_state, use_fates_bgc
  use landunit_varcon                    , only : istcrop, istsoil
  use abortutils                         , only : endrun
  use spmdMod                            , only : masterproc 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, mimics_decomp, decomp_method, use_soil_matrixcn
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  use SparseMatrixMultiplyMod            , only : sparse_matrix_type, vector_type
  use CNVegCarbonStateType               , only : cnveg_carbonstate_type
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: soilbiogeochem_carbonstate_type
     
     ! all c pools involved in decomposition
     real(r8), pointer :: decomp_cpools_vr_col (:,:,:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: decomp0_cpools_vr_col(:,:,:) ! (gC/m3) vertically-resolved C baseline (initial value of this year) in decomposing (litter, cwd, soil) pools in dimension (col,nlev,npools)
     real(r8), pointer :: decomp_cpools_vr_SASUsave_col(:,:,:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: decomp_soilc_vr_col  (:,:)   ! (gC/m3) vertically-resolved decomposing total soil c pool
     real(r8), pointer :: ctrunc_vr_col        (:,:)   ! (gC/m3) vertically-resolved column-level sink for C truncation

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: ctrunc_col              (:)   ! (gC/m2) column-level sink for C truncation
     real(r8), pointer :: totmicc_col             (:)   ! (gC/m2) total microbial carbon
     real(r8), pointer :: totlitc_col             (:)   ! (gC/m2) total litter carbon
     real(r8), pointer :: totlitc_1m_col          (:)   ! (gC/m2) total litter carbon to 1 meter
     real(r8), pointer :: totsomc_col             (:)   ! (gC/m2) total soil organic matter carbon
     real(r8), pointer :: totsomc_1m_col          (:)   ! (gC/m2) total soil organic matter carbon to 1 meter
     real(r8), pointer :: cwdc_col                (:)   ! (gC/m2) coarse woody debris C (diagnostic)
     real(r8), pointer :: decomp_cpools_1m_col    (:,:) ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
     real(r8), pointer :: decomp_cpools_col       (:,:) ! (gC/m2)  decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: dyn_cbal_adjustments_col(:)   ! (gC/m2) adjustments to each column made in this timestep via dynamic column
                                                        ! area adjustments (note: this variable only makes sense at the column-level:
                                                        ! it is meaningless if averaged to the gridcell-level)
     integer           :: restart_file_spinup_state     ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
     real(r8)          :: totvegcthresh                 ! threshold for total vegetation carbon to zero out decomposition pools


     ! Carbon totals, includes soil, cpool and vegetation
     real(r8), pointer :: totc_col                            (:) ! (gC/m2) total column carbon, incl veg and cpool
     real(r8), pointer :: totecosysc_col                      (:) ! (gC/m2) total ecosystem carbon, incl veg but excl cpool 
     real(r8), pointer :: totc_grc                            (:) ! (gC/m2) total gridcell carbon

     
     ! Matrix-cn
     real(r8), pointer :: matrix_cap_decomp_cpools_col    (:,:)   ! (gC/m2) C capacity in decomposing (litter, cwd, soil) N pools in dimension (col,npools)
     real(r8), pointer :: matrix_cap_decomp_cpools_vr_col (:,:,:) ! (gC/m3) vertically-resolved C capacity in decomposing (litter, cwd, soil) pools in dimension(col,nlev,npools)
     real(r8), pointer :: in_acc                          (:,:)   ! (gC/m3/yr) accumulated litter fall C input per year in dimension(col,nlev*npools)
     real(r8), pointer :: in_acc_2d                       (:,:,:) ! (gC/m3/yr) accumulated litter fall C input per year in dimension(col,nlev,npools)
     real(r8), pointer :: tran_acc                        (:,:,:) ! (gC/m3/yr) accumulated C transfers from j to i (col,i,j) per year in dimension(col,nlev*npools,nlev*npools)
     real(r8), pointer :: vert_up_tran_acc                (:,:,:) ! (gC/m3/yr) accumulated upward vertical C transport in dimension(col,nlev,npools)
     real(r8), pointer :: vert_down_tran_acc              (:,:,:) ! (gC/m3/yr) accumulated downward vertical C transport in dimension(col,nlev,npools)
     real(r8), pointer :: exit_acc                        (:,:,:) ! (gC/m3/yr) accumulated exit C in dimension(col,nlev,npools)
     real(r8), pointer :: hori_tran_acc                   (:,:,:) ! (gC/m3/yr) accumulated C transport between pools at the same level in dimension(col,nlev,ntransfers)
     type(sparse_matrix_type) :: AKXcacc                          ! (gC/m3/yr) accumulated N transfers from j to i (col,i,j) per year in dimension(col,nlev*npools,nlev*npools) in sparse matrix type
     type(vector_type) :: matrix_Cinter                           ! (gC/m3)    vertically-resolved decomposing (litter, cwd, soil) N pools in dimension(col,nlev*npools) in vector type

   contains

     procedure , public  :: Init   
     procedure , public  :: SetValues 
     procedure , public  :: Restart
     procedure , public  :: Summary
     procedure , public  :: SetTotVgCThresh
     procedure , public  :: DynamicColumnAdjustments  ! adjust state variables when column areas change
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     


  end type soilbiogeochem_carbonstate_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type, ratio, c12_soilbiogeochem_carbonstate_inst)

    class(soilbiogeochem_carbonstate_type)                       :: this
    type(bounds_type)                     , intent(in)           :: bounds  
    character(len=3)                      , intent(in)           :: carbon_type
    real(r8)                              , intent(in)           :: ratio
    type(soilbiogeochem_carbonstate_type) , intent(in), optional :: c12_soilbiogeochem_carbonstate_inst

    this%totvegcthresh = nan
    call this%InitAllocate ( bounds)
    call this%InitHistory ( bounds, carbon_type )
    if (present(c12_soilbiogeochem_carbonstate_inst)) then
       call this%InitCold  ( bounds, ratio, c12_soilbiogeochem_carbonstate_inst )
    else
       call this%InitCold  ( bounds, ratio) 
    end if

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate( this%decomp_cpools_col    (begc :endc,1:ndecomp_pools))   ; this%decomp_cpools_col    (:,:) = nan
    allocate( this%decomp_cpools_1m_col (begc :endc,1:ndecomp_pools))   ; this%decomp_cpools_1m_col (:,:) = nan
    if(use_soil_matrixcn)then
       allocate( this%matrix_cap_decomp_cpools_col    (begc :endc,1:ndecomp_pools))   ; this%matrix_cap_decomp_cpools_col    (:,:) = nan
    end if

    allocate( this%ctrunc_vr_col(begc :endc,1:nlevdecomp_full)) ; 
    this%ctrunc_vr_col        (:,:) = nan

    allocate(this%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
    this%decomp_cpools_vr_col(:,:,:)= nan
    !matrix-spinup
    if(use_soil_matrixcn)then
       allocate(this%matrix_cap_decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
       this%matrix_cap_decomp_cpools_vr_col(:,:,:)= nan
       allocate(this%decomp0_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
       this%decomp0_cpools_vr_col(:,:,:)= nan
       allocate(this%decomp_cpools_vr_SASUsave_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
       this%decomp_cpools_vr_SASUsave_col(:,:,:)= nan
       allocate(this%in_acc(begc:endc,1:nlevdecomp*ndecomp_pools))
       this%in_acc(:,:)= nan
       allocate(this%tran_acc(begc:endc,1:nlevdecomp*ndecomp_pools,1:nlevdecomp*ndecomp_pools))
       this%tran_acc(:,:,:)= nan

       allocate(this%in_acc_2d(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%in_acc_2d(:,:,:)= nan
       allocate(this%vert_up_tran_acc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%vert_up_tran_acc(:,:,:)= nan
       allocate(this%vert_down_tran_acc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%vert_down_tran_acc(:,:,:)= nan
       allocate(this%exit_acc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%exit_acc(:,:,:)= nan
       allocate(this%hori_tran_acc(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
       this%hori_tran_acc(:,:,:)= nan
       call this%AKXcacc%InitSM(ndecomp_pools*nlevdecomp,begc,endc,decomp_cascade_con%n_all_entries)
       call this%matrix_Cinter%InitV        (ndecomp_pools*nlevdecomp,begc,endc)
    end if
    allocate(this%decomp_soilc_vr_col(begc:endc,1:nlevdecomp_full))  
    this%decomp_soilc_vr_col(:,:)= nan

    allocate(this%ctrunc_col     (begc :endc)) ; this%ctrunc_col     (:) = nan
    allocate(this%cwdc_col       (begc :endc)) ; this%cwdc_col       (:) = nan

    allocate(this%totmicc_col    (begc :endc)) ; this%totmicc_col    (:) = nan
    allocate(this%totlitc_col    (begc :endc)) ; this%totlitc_col    (:) = nan
    allocate(this%totsomc_col    (begc :endc)) ; this%totsomc_col    (:) = nan
    allocate(this%totlitc_1m_col (begc :endc)) ; this%totlitc_1m_col (:) = nan
    allocate(this%totsomc_1m_col (begc :endc)) ; this%totsomc_1m_col (:) = nan
    allocate(this%dyn_cbal_adjustments_col (begc:endc)) ; this%dyn_cbal_adjustments_col (:) = nan

    allocate(this%totc_col                 (begc:endc)) ; this%totc_col                 (:) = nan
    allocate(this%totecosysc_col           (begc:endc)) ; this%totecosysc_col           (:) = nan
    allocate(this%totc_grc                 (begg:endg)) ; this%totc_grc                 (:) = nan
    
    this%restart_file_spinup_state = huge(1)

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, carbon_type) 
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_carbonstate_type) :: this
    type(bounds_type) , intent(in)          :: bounds 
    character(len=3)  , intent(in)          :: carbon_type
    !
    ! !LOCAL VARIABLES:
    integer           :: l
    integer           :: begc ,endc
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    !-------------------------------
    ! C12 state variables - column
    !-------------------------------

    if (carbon_type == 'c12') then

       if ( nlevdecomp_full > 1 ) then
          this%decomp_soilc_vr_col(begc:endc,:) = spval
          call hist_addfld2d (fname='SOILC_vr', units='gC/m^3',  type2d='levsoi', &
               avgflag='A', long_name='SOIL C (vertically resolved)', &
               ptr_col=this%decomp_soilc_vr_col)
       end if

       this%decomp_cpools_col(begc:endc,:) = spval
       do l  = 1, ndecomp_pools
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_cpools_vr_col(:,1:nlevsoi,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'_C_vr'
             longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC/m^3',  type2d='levsoi', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr)
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'_C'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)

          if ( nlevdecomp_full > 1 ) then
             data1dptr => this%decomp_cpools_1m_col(:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'_C_1m'
             longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr, default='inactive')
          endif
       end do
 
       if (decomp_method == mimics_decomp) then
          this%totmicc_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTMICC', units='gC/m^2', &
            avgflag='A', long_name='total microbial carbon', &
            ptr_col=this%totmicc_col)
       end if

       ! Matrix solution history fields
       if(use_soil_matrixcn)then
          this%matrix_cap_decomp_cpools_col(begc:endc,:) = spval
          do l  = 1, ndecomp_pools
             if ( nlevdecomp_full > 1 ) then
                data2dptr => this%matrix_cap_decomp_cpools_vr_col(:,1:nlevsoi,l)
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_Cap_vr'
                longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C capacity (vertically resolved)'
                call hist_addfld2d (fname=fieldname, units='gC/m^3',  type2d='levsoi', &
                     avgflag='I', long_name=longname, &
                     ptr_col=data2dptr)
             endif

             if ( nlevdecomp_full .eq. 1)then
                data1dptr => this%matrix_cap_decomp_cpools_col(:,l)
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_Cap'
                longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C capacity'
                call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                     avgflag='I', long_name=longname, &
                     ptr_col=data1dptr)
             end if

          end do
 
       end if
 
       this%totlitc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTLITC', units='gC/m^2', &
            avgflag='A', long_name='total litter carbon', &
            ptr_col=this%totlitc_col)

       this%totsomc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTSOMC', units='gC/m^2', &
            avgflag='A', long_name='total soil organic matter carbon', &
            ptr_col=this%totsomc_col)

       if ( nlevdecomp_full > 1 ) then
          this%totlitc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTLITC_1m', units='gC/m^2', &
               avgflag='A', long_name='total litter carbon to 1 meter depth', &
               ptr_col=this%totlitc_1m_col)
       end if

       if ( nlevdecomp_full > 1 ) then
          this%totsomc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTSOMC_1m', units='gC/m^2', &
               avgflag='A', long_name='total soil organic matter carbon to 1 meter depth', &
               ptr_col=this%totsomc_1m_col)
       end if

       this%ctrunc_col(begc:endc) = spval
       call hist_addfld1d (fname='COL_CTRUNC', units='gC/m^2',  &
            avgflag='A', long_name='column-level sink for C truncation', &
            ptr_col=this%ctrunc_col, default='inactive')

       this%dyn_cbal_adjustments_col(begc:endc) = spval
       call hist_addfld1d (fname='DYN_COL_SOIL_ADJUSTMENTS_C', units='gC/m^2', &
            avgflag='SUM', &
            long_name='Adjustments in soil carbon due to dynamic column areas; &
            &only makes sense at the column level: should not be averaged to gridcell', &
            ptr_col=this%dyn_cbal_adjustments_col, default='inactive')

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
            avgflag='A', long_name='total column carbon, incl veg and cpool but excl product pools', &
            ptr_col=this%totc_col)

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
            avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool and product pools', &
            ptr_col=this%totecosysc_col)
       
   end if

    !-------------------------------
    ! C13 state variables - column
    !-------------------------------

    if ( carbon_type == 'c13' ) then

       if ( nlevdecomp_full > 1 ) then
          this%decomp_soilc_vr_col(begc:endc,:) = spval
          call hist_addfld2d (fname='C13_SOILC_vr', units='gC13/m^3',  type2d='levsoi', &
               avgflag='A', long_name='C13 SOIL C (vertically resolved)', &
               ptr_col=this%decomp_soilc_vr_col, default='inactive')
       end if

       this%decomp_cpools_vr_col(begc:endc,:,:) = spval
       do l = 1, ndecomp_pools
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_cpools_vr_col(:,1:nlevsoi,l)
             fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'_C_vr'
             longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC13/m^3',  type2d='levsoi', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'_C'
          longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC13/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')
       end do

       if (decomp_method == mimics_decomp) then
          this%totmicc_col(begc:endc) = spval
          call hist_addfld1d (fname='C13_TOTMICC', units='gC/m^2', &
            avgflag='A', long_name='C13 total microbial carbon', &
            ptr_col=this%totmicc_col)
       end if

       ! Matrix solution history fields
       if(use_soil_matrixcn)then
          this%matrix_cap_decomp_cpools_vr_col(begc:endc,:,:) = spval
          do l = 1, ndecomp_pools
             if ( nlevdecomp_full > 1 ) then
                data2dptr => this%matrix_cap_decomp_cpools_vr_col(:,1:nlevsoi,l)
                fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_Cap_vr'
                longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C capacity (vertically resolved)'
                call hist_addfld2d (fname=fieldname, units='gC13/m^3',  type2d='levsoi', &
                     avgflag='I', long_name=longname, &
                     ptr_col=data2dptr, default='inactive')
             endif

             if ( nlevdecomp_full .eq. 1)then
                data1dptr => this%matrix_cap_decomp_cpools_col(:,l)
                fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_Cap'
                longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C capacity'
                call hist_addfld1d (fname=fieldname, units='gC13/m^2', &
                     avgflag='I', long_name=longname, &
                     ptr_col=data1dptr)
             end if
          end do
       end if

       this%totlitc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTLITC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total litter carbon', &
            ptr_col=this%totlitc_col)

       this%totsomc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTSOMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total soil organic matter carbon', &
            ptr_col=this%totsomc_col)

       if ( nlevdecomp_full > 1 ) then
          this%totlitc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C13_TOTLITC_1m', units='gC13/m^2', &
               avgflag='A', long_name='C13 total litter carbon to 1 meter', &
               ptr_col=this%totlitc_1m_col, default='inactive')
       end if

       if ( nlevdecomp_full > 1 ) then
          this%totsomc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C13_TOTSOMC_1m', units='gC13/m^2', &
               avgflag='A', long_name='C13 total soil organic matter carbon to 1 meter', &
               ptr_col=this%totsomc_1m_col, default='inactive')
       endif

       this%ctrunc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_COL_CTRUNC', units='gC13/m^2',  &
            avgflag='A', long_name='C13 column-level sink for C truncation', &
            ptr_col=this%ctrunc_col, default='inactive')

       this%dyn_cbal_adjustments_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_DYN_COL_SOIL_ADJUSTMENTS_C', units='gC13/m^2', &
            avgflag='SUM', &
            long_name='C13 adjustments in soil carbon due to dynamic column areas; &
            &only makes sense at the column level: should not be averaged to gridcell', &
            ptr_col=this%dyn_cbal_adjustments_col, default='inactive')

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total column carbon, incl veg and cpool but excl product pools', &
            ptr_col=this%totc_col, default='inactive')

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTECOSYSC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool and product pools', &
            ptr_col=this%totecosysc_col)
       
    endif

    !-------------------------------
    ! C14 state variables - column
    !-------------------------------

    if ( carbon_type == 'c14' ) then

       if ( nlevdecomp_full > 1 ) then
          this%decomp_soilc_vr_col(begc:endc,:) = spval
          call hist_addfld2d (fname='C14_SOILC_vr', units='gC14/m^3',  type2d='levsoi', &
               avgflag='A', long_name='C14 SOIL C (vertically resolved)', &
               ptr_col=this%decomp_soilc_vr_col)
       end if

       this%decomp_cpools_vr_col(begc:endc,:,:) = spval
       do l = 1, ndecomp_pools
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_cpools_vr_col(:,1:nlevsoi,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'_C_vr'
             longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC14/m^3',  type2d='levsoi', &
                  avgflag='A', long_name=longname, ptr_col=data2dptr, default='inactive')
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'_C'
          longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC14/m^2', &
               avgflag='A', long_name=longname, ptr_col=data1dptr, default='inactive')

          if ( nlevdecomp_full > 1 ) then
             data1dptr => this%decomp_cpools_1m_col(:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'_C_1m'
             longname =  'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                  avgflag='A', long_name=longname, ptr_col=data1dptr, default='inactive')
          endif
       end do

       if (decomp_method == mimics_decomp) then
          this%totmicc_col(begc:endc) = spval
          call hist_addfld1d (fname='C14_TOTMICC', units='gC/m^2', &
            avgflag='A', long_name='C14 total microbial carbon', &
            ptr_col=this%totmicc_col)
       end if

       if(use_soil_matrixcn)then
          this%matrix_cap_decomp_cpools_vr_col(begc:endc,:,:) = spval
          do l = 1, ndecomp_pools
             if ( nlevdecomp_full > 1 ) then
                data2dptr => this%matrix_cap_decomp_cpools_vr_col(:,1:nlevsoi,l)
                fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_Cap_vr'
                longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C capacity (vertically resolved)'
                call hist_addfld2d (fname=fieldname, units='gC14/m^3',  type2d='levsoi', &
                     avgflag='I', long_name=longname, ptr_col=data2dptr, default='inactive')
             endif

             if ( nlevdecomp_full .eq. 1)then
                data1dptr => this%matrix_cap_decomp_cpools_col(:,l)
                fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_Cap'
                longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C capacity'
                call hist_addfld1d (fname=fieldname, units='gC14/m^2', &
                     avgflag='I', long_name=longname, ptr_col=data1dptr)
             end if
          end do
       end if

       this%totlitc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTLITC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total litter carbon', &
            ptr_col=this%totlitc_col)

       this%totsomc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTSOMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total soil organic matter carbon', &
            ptr_col=this%totsomc_col)

       if ( nlevdecomp_full > 1 ) then       
          this%totlitc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C14_TOTLITC_1m', units='gC14/m^2', &
               avgflag='A', long_name='C14 total litter carbon to 1 meter', &
               ptr_col=this%totlitc_1m_col, default='inactive')

          this%totsomc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C14_TOTSOMC_1m', units='gC14/m^2', &
               avgflag='A', long_name='C14 total soil organic matter carbon to 1 meter', &
               ptr_col=this%totsomc_1m_col, default='inactive')
       endif

       this%ctrunc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_COL_CTRUNC', units='gC14/m^2', &
            avgflag='A', long_name='C14 column-level sink for C truncation', &
            ptr_col=this%ctrunc_col, default='inactive')

       this%dyn_cbal_adjustments_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_DYN_COL_SOIL_ADJUSTMENTS_C', units='gC14/m^2', &
            avgflag='SUM', &
            long_name='C14 adjustments in soil carbon due to dynamic column areas; &
            &only makes sense at the column level: should not be averaged to gridcell', &
            ptr_col=this%dyn_cbal_adjustments_col, default='inactive')

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTCOLC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total column carbon, incl veg and cpool but excl product pools', &
            ptr_col=this%totc_col, default='inactive')

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTECOSYSC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total ecosystem carbon, incl veg but excl cpool and product pools', &
            ptr_col=this%totecosysc_col)

    endif

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, ratio, c12_soilbiogeochem_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonstate_type) :: this 
    type(bounds_type) , intent(in)         :: bounds  
    real(r8)          , intent(in)         :: ratio
    type(soilbiogeochem_carbonstate_type), intent(in), optional :: c12_soilbiogeochem_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,j,k,g
    integer :: fc                                        ! filter index
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    !-----------------------------------------------------------------------

    ! initialize column-level variables

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       ! matrix spinup
       if(use_soil_matrixcn)then
          this%in_acc(c,:) = 0._r8
          this%AKXcacc%M(c,:) = 0._r8
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          if (.not. present(c12_soilbiogeochem_carbonstate_inst)) then !c12

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   if (zsoi(j) < decomp_cascade_con%initial_stock_soildepth ) then  !! only initialize upper soil column
                      this%decomp_cpools_vr_col(c,j,k) = decomp_cascade_con%initial_stock(k)
                      if(use_soil_matrixcn)then
                         this%matrix_cap_decomp_cpools_vr_col(c,j,k) = decomp_cascade_con%initial_stock(k)
                      end if
                   else
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                      if(use_soil_matrixcn)then
                         this%matrix_cap_decomp_cpools_vr_col(c,j,k) = 0._r8
                      end if
                   endif
                end do
                this%ctrunc_vr_col(c,j) = 0._r8
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                      if(use_soil_matrixcn)then
                         this%matrix_cap_decomp_cpools_vr_col(c,j,k) = 0._r8
                      end if
                   end do
                   this%ctrunc_vr_col(c,j) = 0._r8
                end do
             end if
             this%decomp_cpools_col(c,1:ndecomp_pools)    = decomp_cascade_con%initial_stock(1:ndecomp_pools)
             this%decomp_cpools_1m_col(c,1:ndecomp_pools) = decomp_cascade_con%initial_stock(1:ndecomp_pools)
             if(use_soil_matrixcn)then
                this%matrix_cap_decomp_cpools_col(c,1:ndecomp_pools)    = decomp_cascade_con%initial_stock(1:ndecomp_pools)
             end if

          else

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   this%decomp_cpools_vr_col(c,j,k) = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(c,j,k) * ratio
                   if(use_soil_matrixcn)then
                      this%matrix_cap_decomp_cpools_vr_col(c,j,k) = c12_soilbiogeochem_carbonstate_inst%matrix_cap_decomp_cpools_vr_col(c,j,k) * ratio
                   end if
                end do
                this%ctrunc_vr_col(c,j) = c12_soilbiogeochem_carbonstate_inst%ctrunc_vr_col(c,j) * ratio
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                      if(use_soil_matrixcn)then
                         this%matrix_cap_decomp_cpools_vr_col(c,j,k) = 0._r8
                      end if
                   end do
                   this%ctrunc_vr_col(c,j) = 0._r8
                end do
             end if
             do k = 1, ndecomp_pools
                this%decomp_cpools_col(c,k)    = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_col(c,k) * ratio
                this%decomp_cpools_1m_col(c,k) = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_1m_col(c,k) * ratio
                if(use_soil_matrixcn)then
                   this%matrix_cap_decomp_cpools_col(c,k)    = c12_soilbiogeochem_carbonstate_inst%matrix_cap_decomp_cpools_col(c,k) * ratio
                end if
             end do

          endif
          if(use_soil_matrixcn)then
             do j = 1, nlevdecomp_full
                do k = 1, ndecomp_pools
                   this%in_acc_2d(c,j,k) = 0._r8
                   this%vert_up_tran_acc(c,j,k) = 0._r8
                   this%vert_down_tran_acc(c,j,k) = 0._r8
                   this%exit_acc(c,j,k) = 0._r8
                   this%decomp0_cpools_vr_col(c,j,k) = max(this%decomp_cpools_vr_col(c,j,k),1.e-30_r8)
                   this%decomp_cpools_vr_SASUsave_col(c,j,k) = 0._r8
                end do
                do k = 1, ndecomp_cascade_transitions
                   this%hori_tran_acc(c,j,k) = 0._r8
                end do
             end do
             do j = 1,decomp_cascade_con%n_all_entries
                this%AKXcacc%M(c,j) = 0._r8
             end do
          end if
       end if


       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          if (present(c12_soilbiogeochem_carbonstate_inst) .and. (.not.col%is_fates(c)) ) then
             this%cwdc_col(c)    = c12_soilbiogeochem_carbonstate_inst%cwdc_col(c) * ratio
          else
             this%cwdc_col(c)    = 0._r8
          end if
          this%ctrunc_col(c)     = 0._r8
          this%totmicc_col(c)    = 0._r8
          this%totlitc_col(c)    = 0._r8
          this%totsomc_col(c)    = 0._r8
          this%totlitc_1m_col(c) = 0._r8
          this%totsomc_1m_col(c) = 0._r8

          this%totc_col(c)       = 0._r8
          this%totecosysc_col(c) = 0._r8
       end if
       
    end do

    do g = bounds%begg, bounds%endg
       this%totc_grc(g)  = 0._r8
    end do
    
    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics
    
    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! initialize fields for special filters
    call this%SetValues (num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag, carbon_type, totvegc_col, c12_soilbiogeochem_carbonstate_inst )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod       , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_time_manager     , only : is_restart, get_nstep
    use shr_const_mod        , only : SHR_CONST_PDB
    use clm_varcon           , only : c14ratio
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_carbonstate_type)                      :: this
    type(bounds_type)                     , intent(in)           :: bounds 
    type(file_desc_t)                     , intent(inout)        :: ncid   ! netcdf id
    character(len=*)                      , intent(in)           :: flag   !'read' or 'write'
    character(len=3)                      , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    real(r8)                              , intent(in)           :: totvegc_col(bounds%begc:bounds%endc) ! (gC/m2) total 
                                                                                                         ! vegetation carbon
    type(soilbiogeochem_carbonstate_type) , intent(in), optional :: c12_soilbiogeochem_carbonstate_inst

    !
    ! !LOCAL VARIABLES:
    integer  :: i,j,k,l,c,fc
    real(r8) :: m                   ! multiplier for the exit_spinup code
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname   ! temporary
    logical  :: readvar
    integer  :: idata
    logical  :: exit_spinup  = .false.
    logical  :: enter_spinup = .false.
    logical  :: found = .false.
    integer  :: i_decomp,j_decomp,i_lev,j_lev
    !------------------------------------------------------------------------

    if (carbon_type == 'c12') then

       do k = 1, ndecomp_pools
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
          ptr2d => this%decomp_cpools_vr_col(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='g/m3', fill_value=spval, &
               scale_by_thickness=.false., &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
          if (flag=='read' .and. .not. readvar) then
             call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
                  errMsg(sourcefile, __LINE__))
          end if
       end do

       if (use_soil_matrixcn)then
          do k = 1, ndecomp_pools
             varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
             ptr2d => this%matrix_cap_decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_Cap_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
             ptr2d => this%decomp0_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"0_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end do
          if(flag=='write')then
             do i = 1,ndecomp_pools
                do j = 1,nlevdecomp
                   this%in_acc_2d(:,j,i) = this%in_acc(:,j+(i-1)*nlevdecomp)
                end do
             end do
             do i = 1,decomp_cascade_con%n_all_entries
                found = .false.
                j_lev    = mod(decomp_cascade_con%all_j(i) - 1,nlevdecomp)  + 1
                j_decomp = (decomp_cascade_con%all_j(i) - j_lev)/nlevdecomp + 1
                i_lev    = mod(decomp_cascade_con%all_i(i) - 1,nlevdecomp)  + 1
                i_decomp = (decomp_cascade_con%all_i(i) - i_lev)/nlevdecomp + 1
                if(i_decomp .eq. j_decomp .and. j_lev - i_lev .eq. 1)then
                   this%vert_up_tran_acc(:,i_lev,i_decomp) = this%AKXcacc%M(:,i)
                   found = .true.
                else
                   if(i_decomp .eq. j_decomp .and. i_lev - j_lev .eq. 1)then
                      this%vert_down_tran_acc(:,i_lev,i_decomp) =  this%AKXcacc%M(:,i)
                      found = .true.
                   else
                      if(i_decomp .eq. j_decomp .and. i_lev .eq. j_lev)then
                         this%exit_acc(:,i_lev,i_decomp) = this%AKXcacc%M(:,i)
                         found = .true.
                      else
                         do k=1,ndecomp_cascade_transitions
                            if(i_decomp .ne. j_decomp .and. i_lev .eq. j_lev .and. &
                               i_decomp .eq. decomp_cascade_con%cascade_receiver_pool(k) .and. &
                               j_decomp .eq. decomp_cascade_con%cascade_donor_pool(k) .and. .not. found)then
                               this%hori_tran_acc(:,i_lev,k) = this%AKXcacc%M(:,i)
                               found = .true.
                            end if
                         end do
                      end if
                   end if
                end if
                if(.not. found) write(iulog,*) 'Error in storing matrix restart variables',i
             end do
          end if
          do k = 1, ndecomp_pools
             varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
             ptr2d => this%in_acc_2d(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_input_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%vert_up_tran_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vert_up_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%vert_down_tran_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vert_down_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%exit_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_exit_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end do
          do i = 1, ndecomp_cascade_transitions
             varname=trim(decomp_cascade_con%cascade_step_name(i))//'c'
             ptr2d => this%hori_tran_acc(:,:,i)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_hori_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end do
          if(flag=='read')then
             do i = 1,ndecomp_pools
                do j = 1,nlevdecomp
                   this%in_acc(:,j+(i-1)*nlevdecomp) = this%in_acc_2d(:,j,i)
                end do
             end do
             do i = 1,decomp_cascade_con%n_all_entries
                found = .false.
                j_lev    = mod(decomp_cascade_con%all_j(i) - 1,nlevdecomp)  + 1
                j_decomp = (decomp_cascade_con%all_j(i) - j_lev)/nlevdecomp + 1
                i_lev    = mod(decomp_cascade_con%all_i(i) - 1,nlevdecomp)  + 1
                i_decomp = (decomp_cascade_con%all_i(i) - i_lev)/nlevdecomp + 1
                if(i_decomp .eq. j_decomp .and. j_lev - i_lev .eq. 1)then
                   this%AKXcacc%M(:,i) = this%vert_up_tran_acc(:,i_lev,i_decomp) 
                   found = .true.
                else
                   if(i_decomp .eq. j_decomp .and. i_lev - j_lev .eq. 1)then
                      this%AKXcacc%M(:,i) = this%vert_down_tran_acc(:,i_lev,i_decomp) 
                      found = .true.
                   else
                      if(i_decomp .eq. j_decomp .and. i_lev .eq. j_lev)then
                         this%AKXcacc%M(:,i) = this%exit_acc(:,i_lev,i_decomp) 
                         found = .true.
                      else
                         do k=1,ndecomp_cascade_transitions
                            if(i_decomp .ne. j_decomp .and. i_lev .eq. j_lev .and. &
                               i_decomp .eq. decomp_cascade_con%cascade_receiver_pool(k) .and. &
                               j_decomp .eq. decomp_cascade_con%cascade_donor_pool(k) .and. .not. found)then
                               this%AKXcacc%M(:,i) = this%hori_tran_acc(:,i_lev,k) 
                               found = .true.
                            end if
                         end do
                      end if
                   end if
                end if
                if(.not. found) write(iulog,*) 'Error in storing matrix restart variables',i
             end do
          end if
       end if

       ptr2d => this%ctrunc_vr_col
       call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc_vr', xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='gC/m3', fill_value=spval, &
            scale_by_thickness=.false., &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
               errMsg(sourcefile, __LINE__))
       end if

    end if

    !--------------------------------
    ! C13 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c13' ) then

       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
          ptr2d => this%decomp_cpools_vr_col(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='g/m3', fill_value=spval, &
               scale_by_thickness=.false., &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
          if(use_soil_matrixcn)then
             ptr2d => this%matrix_cap_decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_Cap_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
             ptr2d => this%decomp0_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"0_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col' &
                  // ' with atmospheric c13 value for: '//trim(varname)
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (this%decomp_cpools_vr_col(i,j,k) /= spval .and. .not. isnan(this%decomp_cpools_vr_col(i,j,k)) ) then
                      this%decomp_cpools_vr_col(i,j,k) = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(i,j,k) * c3_r2
                   endif
                   if(use_soil_matrixcn)then
                      if (this%matrix_cap_decomp_cpools_vr_col(i,j,k) /= spval .and. .not. isnan(this%matrix_cap_decomp_cpools_vr_col(i,j,k)) ) then
                         this%matrix_cap_decomp_cpools_vr_col(i,j,k) = c12_soilbiogeochem_carbonstate_inst%matrix_cap_decomp_cpools_vr_col(i,j,k) * c3_r2
                      endif
                   end if
                end do
             end do
          end if
       end do

       if (use_soil_matrixcn)then
          if(flag=='write')then
             do i = 1,ndecomp_pools
                do j = 1,nlevdecomp
                   this%in_acc_2d(:,j,i) = this%in_acc(:,j+(i-1)*nlevdecomp)
                end do
             end do
             do i = 1,decomp_cascade_con%n_all_entries
                found = .false.
                j_lev    = mod(decomp_cascade_con%all_j(i) - 1,nlevdecomp)  + 1
                j_decomp = (decomp_cascade_con%all_j(i) - j_lev)/nlevdecomp + 1
                i_lev    = mod(decomp_cascade_con%all_i(i) - 1,nlevdecomp)  + 1
                i_decomp = (decomp_cascade_con%all_i(i) - i_lev)/nlevdecomp + 1
                if(i_decomp .eq. j_decomp .and. j_lev - i_lev .eq. 1)then
                   this%vert_up_tran_acc(:,i_lev,i_decomp) = this%AKXcacc%M(:,i)
                   found = .true.
                else
                   if(i_decomp .eq. j_decomp .and. i_lev - j_lev .eq. 1)then
                      this%vert_down_tran_acc(:,i_lev,i_decomp) =  this%AKXcacc%M(:,i)
                      found = .true.
                   else
                      if(i_decomp .eq. j_decomp .and. i_lev .eq. j_lev)then
                         this%exit_acc(:,i_lev,i_decomp) = this%AKXcacc%M(:,i)
                         found = .true.
                      else
                         do k=1,ndecomp_cascade_transitions
                            if(i_decomp .ne. j_decomp .and. i_lev .eq. j_lev .and. &
                               i_decomp .eq. decomp_cascade_con%cascade_receiver_pool(k) .and. &
                               j_decomp .eq. decomp_cascade_con%cascade_donor_pool(k) .and. .not. found)then
                               this%hori_tran_acc(:,i_lev,k) = this%AKXcacc%M(:,i)
                               found = .true.
                            end if
                         end do
                      end if
                   end if
                end if
                if(.not. found) write(iulog,*) 'Error in storing matrix restart variables',i
             end do
          end if
          do k = 1, ndecomp_pools
             varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
             ptr2d => this%in_acc_2d(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_input_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%vert_up_tran_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vert_up_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%vert_down_tran_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vert_down_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%exit_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_exit_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end do
          do i = 1, ndecomp_cascade_transitions
             varname=trim(decomp_cascade_con%cascade_step_name(i))//'c_13'
             ptr2d => this%hori_tran_acc(:,:,i)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_hori_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end do
          if(flag=='read')then
             do i = 1,ndecomp_pools
                do j = 1,nlevdecomp
                   this%in_acc(:,j+(i-1)*nlevdecomp) = this%in_acc_2d(:,j,i)
                end do
             end do
             do i = 1,decomp_cascade_con%n_all_entries
                found = .false.
                j_lev    = mod(decomp_cascade_con%all_j(i) - 1,nlevdecomp)  + 1
                j_decomp = (decomp_cascade_con%all_j(i) - j_lev)/nlevdecomp + 1
                i_lev    = mod(decomp_cascade_con%all_i(i) - 1,nlevdecomp)  + 1
                i_decomp = (decomp_cascade_con%all_i(i) - i_lev)/nlevdecomp + 1
                if(i_decomp .eq. j_decomp .and. j_lev - i_lev .eq. 1)then
                   this%AKXcacc%M(:,i) = this%vert_up_tran_acc(:,i_lev,i_decomp) 
                   found = .true.
                else
                   if(i_decomp .eq. j_decomp .and. i_lev - j_lev .eq. 1)then
                      this%AKXcacc%M(:,i) = this%vert_down_tran_acc(:,i_lev,i_decomp) 
                      found = .true.
                   else
                      if(i_decomp .eq. j_decomp .and. i_lev .eq. j_lev)then
                         this%AKXcacc%M(:,i) = this%exit_acc(:,i_lev,i_decomp) 
                         found = .true.
                      else
                         do k=1,ndecomp_cascade_transitions
                            if(i_decomp .ne. j_decomp .and. i_lev .eq. j_lev .and. &
                               i_decomp .eq. decomp_cascade_con%cascade_receiver_pool(k) .and. &
                               j_decomp .eq. decomp_cascade_con%cascade_donor_pool(k) .and. .not. found)then
                               this%AKXcacc%M(:,i) = this%hori_tran_acc(:,i_lev,k) 
                               found = .true.
                            end if
                         end do
                      end if
                   end if
                end if
                if(.not. found) write(iulog,*) 'Error in storing matrix restart variables',i
             end do
          end if
       end if

       ptr2d => this%ctrunc_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='gC/m3', fill_value=spval, &
            scale_by_thickness=.false., &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    end if

    !--------------------------------
    ! C14 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c14' ) then

       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
          ptr2d => this%decomp_cpools_vr_col(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='g/m3', fill_value=spval, &
               scale_by_thickness=.false., &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
          if(use_soil_matrixcn)then
             ptr2d => this%matrix_cap_decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_Cap_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
             ptr2d => this%decomp0_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"0_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col with atmospheric c14 value for: '//&
                  trim(varname)
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (this%decomp_cpools_vr_col(i,j,k) /= spval .and. .not. isnan(this%decomp_cpools_vr_col(i,j,k)) ) then
                      this%decomp_cpools_vr_col(i,j,k) = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(i,j,k) * c3_r2
                   endif
                   if(use_soil_matrixcn)then
                      if (this%matrix_cap_decomp_cpools_vr_col(i,j,k) /= spval .and. .not. isnan(this%matrix_cap_decomp_cpools_vr_col(i,j,k)) ) then
                         this%matrix_cap_decomp_cpools_vr_col(i,j,k) = c12_soilbiogeochem_carbonstate_inst%matrix_cap_decomp_cpools_vr_col(i,j,k) * c3_r2
                      endif
                   end if
                end do
             end do
          end if
       end do

       if (use_soil_matrixcn)then
          if(flag=='write')then
             do i = 1,ndecomp_pools
                do j = 1,nlevdecomp
                   this%in_acc_2d(:,j,i) = this%in_acc(:,j+(i-1)*nlevdecomp)
                end do
             end do
             do i = 1,decomp_cascade_con%n_all_entries
                found = .false.
                j_lev    = mod(decomp_cascade_con%all_j(i) - 1,nlevdecomp)  + 1
                j_decomp = (decomp_cascade_con%all_j(i) - j_lev)/nlevdecomp + 1
                i_lev    = mod(decomp_cascade_con%all_i(i) - 1,nlevdecomp)  + 1
                i_decomp = (decomp_cascade_con%all_i(i) - i_lev)/nlevdecomp + 1
                if(i_decomp .eq. j_decomp .and. j_lev - i_lev .eq. 1)then
                   this%vert_up_tran_acc(:,i_lev,i_decomp) = this%AKXcacc%M(:,i)
                   found = .true.
                else
                   if(i_decomp .eq. j_decomp .and. i_lev - j_lev .eq. 1)then
                      this%vert_down_tran_acc(:,i_lev,i_decomp) =  this%AKXcacc%M(:,i)
                      found = .true.
                   else
                      if(i_decomp .eq. j_decomp .and. i_lev .eq. j_lev)then
                         this%exit_acc(:,i_lev,i_decomp) = this%AKXcacc%M(:,i)
                         found = .true.
                      else
                         do k=1,ndecomp_cascade_transitions
                            if(i_decomp .ne. j_decomp .and. i_lev .eq. j_lev .and. &
                               i_decomp .eq. decomp_cascade_con%cascade_receiver_pool(k) .and. &
                               j_decomp .eq. decomp_cascade_con%cascade_donor_pool(k) .and. .not. found)then
                               this%hori_tran_acc(:,i_lev,k) = this%AKXcacc%M(:,i)
                               found = .true.
                            end if
                         end do
                      end if
                   end if
                end if
                if(.not. found) write(iulog,*) 'Error in storing matrix restart variables',i
             end do
          end if
          do k = 1, ndecomp_pools
             varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
             ptr2d => this%in_acc_2d(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_input_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%vert_up_tran_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vert_up_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%vert_down_tran_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vert_down_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
                ptr2d => this%exit_acc(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_exit_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end do
          do i = 1, ndecomp_cascade_transitions
             varname=trim(decomp_cascade_con%cascade_step_name(i))//'c_14'
             ptr2d => this%hori_tran_acc(:,:,i)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_hori_tran_acc_vr", xtype=ncd_double,  &
                dim1name='column', dim2name='levgrnd', switchdim=.true., &
                long_name='',  units='', fill_value=spval, scale_by_thickness=.false., &
                interpinic_flag='interp', readvar=readvar, data=ptr2d)
          end do
          if(flag=='read')then
             do i = 1,ndecomp_pools
                do j = 1,nlevdecomp
                   this%in_acc(:,j+(i-1)*nlevdecomp) = this%in_acc_2d(:,j,i)
                end do
             end do
             do i = 1,decomp_cascade_con%n_all_entries
                found = .false.
                j_lev    = mod(decomp_cascade_con%all_j(i) - 1,nlevdecomp)  + 1
                j_decomp = (decomp_cascade_con%all_j(i) - j_lev)/nlevdecomp + 1
                i_lev    = mod(decomp_cascade_con%all_i(i) - 1,nlevdecomp)  + 1
                i_decomp = (decomp_cascade_con%all_i(i) - i_lev)/nlevdecomp + 1
                if(i_decomp .eq. j_decomp .and. j_lev - i_lev .eq. 1)then
                   this%AKXcacc%M(:,i) = this%vert_up_tran_acc(:,i_lev,i_decomp) 
                   found = .true.
                else
                   if(i_decomp .eq. j_decomp .and. i_lev - j_lev .eq. 1)then
                      this%AKXcacc%M(:,i) = this%vert_down_tran_acc(:,i_lev,i_decomp) 
                      found = .true.
                   else
                      if(i_decomp .eq. j_decomp .and. i_lev .eq. j_lev)then
                         this%AKXcacc%M(:,i) = this%exit_acc(:,i_lev,i_decomp) 
                         found = .true.
                      else
                         do k=1,ndecomp_cascade_transitions
                            if(i_decomp .ne. j_decomp .and. i_lev .eq. j_lev .and. &
                               i_decomp .eq. decomp_cascade_con%cascade_receiver_pool(k) .and. &
                               j_decomp .eq. decomp_cascade_con%cascade_donor_pool(k) .and. .not. found)then
                               this%AKXcacc%M(:,i) = this%hori_tran_acc(:,i_lev,k) 
                               found = .true.
                            end if
                         end do
                      end if
                   end if
                end if
                if(.not. found) write(iulog,*) 'Error in storing matrix restart variables',i
             end do
          end if
       end if

       ptr2d => this%ctrunc_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='gC/m3', fill_value=spval, &
            scale_by_thickness=.false., &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)

    end if

    !--------------------------------
    ! Spinup state
    !--------------------------------

        if (carbon_type == 'c12') then
           if (flag == 'write') idata = spinup_state
           call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
             long_name='Spinup state of the model that wrote this restart file: ' &
             // ' 0 = normal model mode, 1 = AD spinup', units='', &
             interpinic_flag='copy', readvar=readvar,  data=idata)
           if (flag == 'read') then
              if (readvar) then
                 this%restart_file_spinup_state = idata
              else
                 call endrun(msg=' CNRest: spinup_state was not on the restart file and is required' // &
                   errMsg(sourcefile, __LINE__))
              end if
           end if
        else
           this%restart_file_spinup_state = c12_soilbiogeochem_carbonstate_inst%restart_file_spinup_state
        endif

        ! now compare the model and restart file spinup states, and either take the 
        ! model into spinup mode or out of it if they are not identical
        ! taking model out of spinup mode requires multiplying each decomposing pool 
        ! by the associated AD factor.
        ! putting model into spinup mode requires dividing each decomposing pool 
        ! by the associated AD factor.
        ! only allow this to occur on first timestep of model run.
        
        if (flag == 'read' .and. spinup_state /= this%restart_file_spinup_state ) then
           if (spinup_state == 0 .and. this%restart_file_spinup_state >= 1 ) then
              if ( masterproc ) write(iulog,*) ' CNRest: taking ',carbon_type,' SOM pools out of AD spinup mode'
              exit_spinup = .true.
           else if (spinup_state >= 1 .and. this%restart_file_spinup_state == 0 ) then
              if ( masterproc ) write(iulog,*) ' CNRest: taking ',carbon_type,' SOM pools into AD spinup mode'
              enter_spinup = .true.
           else
              call endrun(msg=' CNRest: error in entering/exiting spinup.  spinup_state ' &
                   // ' != restart_file_spinup_state, but do not know what to do'//&
                   errMsg(sourcefile, __LINE__))
           end if
           if (get_nstep() >= 2) then
              call endrun(msg=' CNRest: error in entering/exiting spinup - should occur only when nstep = 1'//&
                   errMsg(sourcefile, __LINE__))
           endif
           if ( exit_spinup .and. isnan(this%totvegcthresh) )then
              call endrun(msg=' CNRest: error in exit spinup - totvegcthresh was not set with SetTotVgCThresh'//&
                   errMsg(sourcefile, __LINE__))
           end if
           do k = 1, ndecomp_pools
              if ( exit_spinup ) then
                 m = decomp_cascade_con%spinup_factor(k)
              else if ( enter_spinup ) then
                 m = 1. / decomp_cascade_con%spinup_factor(k)
              end if
              do c = bounds%begc, bounds%endc
                 l = col%landunit(c)
                 do j = 1, nlevdecomp_full
                    if ( abs(m - 1._r8) .gt. 0.000001_r8 .and. exit_spinup) then
                       this%decomp_cpools_vr_col(c,j,k) = this%decomp_cpools_vr_col(c,j,k) * m * &
                            get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
                       ! If there is no vegetation carbon, implying that all vegetation has died, then 
                       ! reset decomp pools to near zero during exit_spinup to avoid very 
                       ! large and inert soil carbon stocks; note that only pools with spinup factor > 1 
                       ! will be affected, which means that total SOMC and LITC pools will not be set to 0.
                       if (totvegc_col(c) <= this%totvegcthresh .and. lun%itype(l) /= istcrop) then 
                         this%decomp_cpools_vr_col(c,j,k) = 0.0_r8
                       endif
                    elseif ( abs(m - 1._r8) .gt. 0.000001_r8 .and. enter_spinup) then
                       this%decomp_cpools_vr_col(c,j,k) = this%decomp_cpools_vr_col(c,j,k) * m / &
                            get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
                    else
                       this%decomp_cpools_vr_col(c,j,k) = this%decomp_cpools_vr_col(c,j,k) * m 
                    endif
                 end do
              end do
           end do
        end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state variables
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_carbonstate_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_column
       i = filter_column(fi)
       if ( .not. col%is_fates(i) ) then
          this%cwdc_col(i)       = value_column
       end if
       this%ctrunc_col(i)     = value_column
       this%totmicc_col(i)    = value_column
       this%totlitc_col(i)    = value_column
       this%totlitc_1m_col(i) = value_column
       this%totsomc_col(i)    = value_column
       this%totsomc_1m_col(i) = value_column
       this%totc_col(i)       = value_column
       this%totecosysc_col(i) = value_column
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%ctrunc_vr_col(i,j) = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cpools_col(i,k) = value_column
          this%decomp_cpools_1m_col(i,k) = value_column
          if(use_soil_matrixcn)then
             this%matrix_cap_decomp_cpools_col(i,k) = value_column
          end if
       end do
    end do

    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cpools_vr_col(i,j,k) = value_column
             if(use_soil_matrixcn)then
                this%matrix_cap_decomp_cpools_vr_col(i,j,k) = value_column
                this%decomp0_cpools_vr_col(i,j,k) = value_column
             end if
          end do
       end do
    end do

    if(use_soil_matrixcn)then
       do j = 1,nlevdecomp
          do k = 1, ndecomp_pools
             do fi = 1, num_column
                i = filter_column(fi)
                this%in_acc_2d(i,j,k)          = value_column
                this%vert_up_tran_acc(i,j,k)   = value_column
                this%vert_down_tran_acc(i,j,k) = value_column
                this%exit_acc(i,j,k) = value_column
             end do
          end do
          do k = 1, ndecomp_cascade_transitions
             do fi = 1, num_column
                i = filter_column(fi)
                this%hori_tran_acc(i,j,k)   = value_column
             end do
          end do
       end do
    end if
    
    if(use_soil_matrixcn)then
       do j = 1,decomp_cascade_con%n_all_entries
          do fi = 1, num_column
             i = filter_column(fi)
             this%AKXcacc%M(i,j) = value_column
          end do
       end do
    end if

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_allc, filter_allc, num_bgc_soilc, filter_bgc_soilc,cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Perform column-level carbon summary calculations
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonstate_type)          :: this
    type(bounds_type)               , intent(in)    :: bounds          
    integer                         , intent(in)    :: num_allc       ! number of columns in soil filter
    integer                         , intent(in)    :: filter_allc(:) ! filter for all active columns
    integer                         , intent(in)    :: num_bgc_soilc       ! number of columns in soil filter
    integer                         , intent(in)    :: filter_bgc_soilc(:) ! filter for all active columns
    type(cnveg_carbonstate_type)    , intent(inout) :: cnveg_carbonstate_inst
    
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l       ! indices
    integer  :: fc            ! filter indices
    real(r8) :: maxdepth      ! depth to integrate soil variables
    integer  :: num_local     ! Either num_bgc_soilc or num_allc, depending
                              ! on if its a fates run, its different because
                              ! the cnveg variables are not allocated w/ fates
    real(r8) :: ecovegc_col
    real(r8) :: totvegc_col
    !-----------------------------------------------------------------------

    ! vertically integrate each of the decomposing C pools
    do l = 1, ndecomp_pools
       do fc = 1,num_allc
          c = filter_allc(fc)
          this%decomp_cpools_col(c,l) = 0._r8
          if(use_soil_matrixcn)then
             this%matrix_cap_decomp_cpools_col(c,l) = 0._r8
          end if
       end do
    end do
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%decomp_cpools_col(c,l) = &
                  this%decomp_cpools_col(c,l) + &
                  this%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
             if(use_soil_matrixcn)then
                this%matrix_cap_decomp_cpools_col(c,l) = &
                     this%matrix_cap_decomp_cpools_col(c,l) + &
                     this%matrix_cap_decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
             end if
          end do
       end do
    end do

    if ( nlevdecomp > 1) then

       ! vertically integrate each of the decomposing C pools to 1 meter
       maxdepth = 1._r8
       do l = 1, ndecomp_pools
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%decomp_cpools_1m_col(c,l) = 0._r8
          end do
       end do
       do l = 1, ndecomp_pools
          do j = 1, nlevdecomp
             if ( zisoi(j) <= maxdepth ) then
                do fc = 1,num_allc
                   c = filter_allc(fc)
                   this%decomp_cpools_1m_col(c,l) = &
                        this%decomp_cpools_1m_col(c,l) + &
                        this%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
                end do
             elseif ( zisoi(j-1) < maxdepth ) then
                do fc = 1,num_allc
                   c = filter_allc(fc)
                   this%decomp_cpools_1m_col(c,l) = &
                        this%decomp_cpools_1m_col(c,l) + &
                        this%decomp_cpools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
                end do
             endif
          end do
       end do

    endif

    ! Add soil carbon pools together to produce vertically-resolved decomposing total soil c pool
    if ( nlevdecomp_full > 1 ) then
       do j = 1, nlevdecomp
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%decomp_soilc_vr_col(c,j) = 0._r8
          end do
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_soil(l) ) then
             do j = 1, nlevdecomp
                do fc = 1,num_allc
                   c = filter_allc(fc)
                   this%decomp_soilc_vr_col(c,j) = this%decomp_soilc_vr_col(c,j) + &
                        this%decomp_cpools_vr_col(c,j,l)
                end do
             end do
          end if
       end do
    end if

    ! truncation carbon
    do fc = 1,num_allc
       c = filter_allc(fc)
       this%ctrunc_col(c) = 0._r8
    end do
    do j = 1, nlevdecomp
       do fc = 1,num_allc
          c = filter_allc(fc)
          this%ctrunc_col(c) = &
               this%ctrunc_col(c) + &
               this%ctrunc_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! total litter carbon in the top meter (TOTLITC_1m)
    if ( nlevdecomp > 1) then
       do fc = 1,num_allc
          c = filter_allc(fc)
          this%totlitc_1m_col(c) = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_litter(l) ) then
             do fc = 1,num_allc
                c = filter_allc(fc)
                this%totlitc_1m_col(c) = this%totlitc_1m_col(c) + &
                     this%decomp_cpools_1m_col(c,l)
             end do
          endif
       end do
    end if

    ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
    if ( nlevdecomp > 1) then
       do fc = 1,num_allc
          c = filter_allc(fc)
          this%totsomc_1m_col(c) = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_soil(l) ) then
             do fc = 1,num_allc
                c = filter_allc(fc)
                this%totsomc_1m_col(c) = this%totsomc_1m_col(c) + this%decomp_cpools_1m_col(c,l)
             end do
          end if
       end do
    end if

    ! total microbial carbon (TOTMICC)
    do fc = 1,num_allc
       c = filter_allc(fc)
       this%totmicc_col(c) = 0._r8
    end do
    do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_microbe(l) ) then
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%totmicc_col(c) = this%totmicc_col(c) + this%decomp_cpools_col(c,l)
          end do
       endif
    end do

    ! total litter carbon (TOTLITC)
    do fc = 1,num_allc
       c = filter_allc(fc)
       this%totlitc_col(c) = 0._r8
    end do
    do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(l) ) then
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%totlitc_col(c) = this%totlitc_col(c) + this%decomp_cpools_col(c,l)
          end do
       endif
    end do

    ! total soil organic matter carbon (TOTSOMC)
    do fc = 1,num_allc
       c = filter_allc(fc)
       this%totsomc_col(c) = 0._r8
    end do
    do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_soil(l) ) then
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%totsomc_col(c) = this%totsomc_col(c) + this%decomp_cpools_col(c,l)
          end do
       end if
    end do

    do fc = 1,num_allc
       c = filter_allc(fc)
       ! coarse woody debris carbon
       this%cwdc_col(c) = 0._r8
    end do
    
    if (use_fates_bgc) then
       num_local = num_bgc_soilc
    else
       num_local = num_allc
    end if
    do fc = 1,num_local
       if(use_fates_bgc) then
          c = filter_bgc_soilc(fc)
       else
          c = filter_allc(fc)
       end if
       if(col%is_fates(c)) then
          totvegc_col = 0._r8
          ecovegc_col = 0._r8
       else
          do l = 1, ndecomp_pools
             if ( decomp_cascade_con%is_cwd(l) ) then
                this%cwdc_col(c) = this%cwdc_col(c) + this%decomp_cpools_col(c,l)
             end if
          end do
          totvegc_col = cnveg_carbonstate_inst%totc_p2c_col(c)
          ecovegc_col = cnveg_carbonstate_inst%totvegc_col(c)
       end if
       
       ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
       this%totecosysc_col(c) =   &
            this%cwdc_col(c)    + &
            this%totmicc_col(c) + &
            this%totlitc_col(c) + &
            this%totsomc_col(c) + &
            ecovegc_col
       ! total column carbon, including veg and cpool (TOTCOLC)
       this%totc_col(c) =         &
            this%cwdc_col(c)    + &
            this%totmicc_col(c) + &
            this%totlitc_col(c) + &
            this%totsomc_col(c) + &
            this%ctrunc_col(c)  + &
            totvegc_col
    end do
       
  end subroutine Summary

  !------------------------------------------------------------------------
  subroutine SetTotVgCThresh(this, totvegcthresh)

    class(soilbiogeochem_carbonstate_type)                       :: this
    real(r8)                              , intent(in)           :: totvegcthresh

    if ( totvegcthresh <= 0.0_r8 )then
        call endrun(msg=' ERROR totvegcthresh is zero or negative and should be > 0'//&
               errMsg(sourcefile, __LINE__))
    end if
    this%totvegcthresh = totvegcthresh

  end subroutine SetTotVgCThresh


  !-----------------------------------------------------------------------
  subroutine DynamicColumnAdjustments(this, bounds, clump_index, column_state_updater)
    !
    ! !DESCRIPTION:
    ! Adjust state variables when column areas change due to dynamic landuse
    !
    ! !USES:
    use dynColumnStateUpdaterMod, only : column_state_updater_type
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonstate_type) , intent(inout) :: this
    type(bounds_type)                      , intent(in)    :: bounds

    ! Index of clump on which we're currently operating. Note that this implies that this
    ! routine must be called from within a clump loop.
    integer                                , intent(in)    :: clump_index

    type(column_state_updater_type)        , intent(in)    :: column_state_updater
    !
    ! !LOCAL VARIABLES:
    integer :: j  ! level
    integer :: l  ! decomp pool
    real(r8) :: adjustment_one_level(bounds%begc:bounds%endc)
    integer :: begc, endc

    character(len=*), parameter :: subname = 'DynamicColumnAdjustments'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    this%dyn_cbal_adjustments_col(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call column_state_updater%update_column_state_no_special_handling( &
               bounds = bounds, &
               clump_index = clump_index, &
               var    = this%decomp_cpools_vr_col(begc:endc, j, l), &
               adjustment = adjustment_one_level(begc:endc))

               this%dyn_cbal_adjustments_col(begc:endc) = &
               this%dyn_cbal_adjustments_col(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       end do
    end do

    do j = 1, nlevdecomp
       call column_state_updater%update_column_state_no_special_handling( &
            bounds = bounds, &
            clump_index = clump_index, &
            var    = this%ctrunc_vr_col(begc:endc, j), &
            adjustment = adjustment_one_level(begc:endc))
       this%dyn_cbal_adjustments_col(begc:endc) = &
            this%dyn_cbal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
    end do

  end subroutine DynamicColumnAdjustments


end module SoilBiogeochemCarbonStateType
