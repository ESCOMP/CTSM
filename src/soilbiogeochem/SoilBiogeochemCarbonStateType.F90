module SoilBiogeochemCarbonStateType

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use decompMod                          , only : bounds_type
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp, nlevsoi
  use clm_varcon                         , only : spval, ispval, dzsoi_decomp, zisoi, zsoi, c3_r2
  use clm_varctl                         , only : iulog, use_vertsoilc, spinup_state, use_fates 
  use landunit_varcon                    , only : istcrop, istsoil
  use abortutils                         , only : endrun
  use spmdMod                            , only : masterproc 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: soilbiogeochem_carbonstate_type
     
     ! all c pools involved in decomposition
     real(r8), pointer :: decomp_cpools_vr_col (:,:,:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: decomp_soilc_vr_col  (:,:)   ! (gC/m3) vertically-resolved decomposing total soil c pool
     real(r8), pointer :: ctrunc_vr_col        (:,:)   ! (gC/m3) vertically-resolved column-level sink for C truncation

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: ctrunc_col               (:) ! (gC/m2) column-level sink for C truncation
     real(r8), pointer :: totlitc_col          (:)     ! (gC/m2) total litter carbon
     real(r8), pointer :: totlitc_1m_col       (:)     ! (gC/m2) total litter carbon to 1 meter
     real(r8), pointer :: totsomc_col          (:)     ! (gC/m2) total soil organic matter carbon
     real(r8), pointer :: totsomc_1m_col       (:)     ! (gC/m2) total soil organic matter carbon to 1 meter
     real(r8), pointer :: cwdc_col             (:)     ! (gC/m2) coarse woody debris C (diagnostic)
     real(r8), pointer :: decomp_cpools_1m_col (:,:)   ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
     real(r8), pointer :: decomp_cpools_col    (:,:)   ! (gC/m2)  decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: dyn_cbal_adjustments_col (:) ! (gC/m2) adjustments to each column made in this timestep via dynamic column area adjustments (note: this variable only makes sense at the column-level: it is meaningless if averaged to the gridcell-level)
     integer  :: restart_file_spinup_state             ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
     real(r8)          :: totvegcthresh                ! threshold for total vegetation carbon to zero out decomposition pools

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
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    allocate( this%decomp_cpools_col    (begc :endc,1:ndecomp_pools))   ; this%decomp_cpools_col    (:,:) = nan
    allocate( this%decomp_cpools_1m_col (begc :endc,1:ndecomp_pools))   ; this%decomp_cpools_1m_col (:,:) = nan

    allocate( this%ctrunc_vr_col(begc :endc,1:nlevdecomp_full)) ; 
    this%ctrunc_vr_col        (:,:) = nan

    allocate(this%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
    this%decomp_cpools_vr_col(:,:,:)= nan
    allocate(this%decomp_soilc_vr_col(begc:endc,1:nlevdecomp_full))  
    this%decomp_soilc_vr_col(:,:)= nan

    allocate(this%ctrunc_col     (begc :endc)) ; this%ctrunc_col     (:) = nan
    if ( .not. use_fates ) then
       allocate(this%cwdc_col       (begc :endc)) ; this%cwdc_col       (:) = nan
    endif
    allocate(this%totlitc_col    (begc :endc)) ; this%totlitc_col    (:) = nan
    allocate(this%totsomc_col    (begc :endc)) ; this%totsomc_col    (:) = nan
    allocate(this%totlitc_1m_col (begc :endc)) ; this%totlitc_1m_col (:) = nan
    allocate(this%totsomc_1m_col (begc :endc)) ; this%totsomc_1m_col (:) = nan
    allocate(this%dyn_cbal_adjustments_col (begc:endc)) ; this%dyn_cbal_adjustments_col (:) = nan

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
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC/m^3',  type2d='levsoi', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr)
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)

          if ( nlevdecomp_full > 1 ) then
             data1dptr => this%decomp_cpools_1m_col(:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_1m'
             longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr, default='inactive')
          endif
       end do
 
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
             fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC13/m^3',  type2d='levsoi', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC13/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')
       end do

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
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC14/m^3',  type2d='levsoi', &
                  avgflag='A', long_name=longname, ptr_col=data2dptr, default='inactive')
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC14/m^2', &
               avgflag='A', long_name=longname, ptr_col=data1dptr, default='inactive')

          if ( nlevdecomp_full > 1 ) then
             data1dptr => this%decomp_cpools_1m_col(:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_1m'
             longname =  'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                  avgflag='A', long_name=longname, ptr_col=data1dptr, default='inactive')
          endif
       end do

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
    integer :: p,c,l,j,k
    integer :: fc                                        ! filter index
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    !-----------------------------------------------------------------------

    ! initialize column-level variables

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          if (.not. present(c12_soilbiogeochem_carbonstate_inst)) then !c12

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   if (zsoi(j) < decomp_cascade_con%initial_stock_soildepth ) then  !! only initialize upper soil column
                      this%decomp_cpools_vr_col(c,j,k) = decomp_cascade_con%initial_stock(k)
                   else
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                   endif
                end do
                this%ctrunc_vr_col(c,j) = 0._r8
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                   end do
                   this%ctrunc_vr_col(c,j) = 0._r8
                end do
             end if
             this%decomp_cpools_col(c,1:ndecomp_pools)    = decomp_cascade_con%initial_stock(1:ndecomp_pools)
             this%decomp_cpools_1m_col(c,1:ndecomp_pools) = decomp_cascade_con%initial_stock(1:ndecomp_pools)

          else

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   this%decomp_cpools_vr_col(c,j,k) = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(c,j,k) * ratio
                end do
                this%ctrunc_vr_col(c,j) = c12_soilbiogeochem_carbonstate_inst%ctrunc_vr_col(c,j) * ratio
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                   end do
                   this%ctrunc_vr_col(c,j) = 0._r8
                end do
             end if
             do k = 1, ndecomp_pools
                this%decomp_cpools_col(c,k)    = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_col(c,k) * ratio
                this%decomp_cpools_1m_col(c,k) = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_1m_col(c,k) * ratio
             end do

          endif
       end if

       if ( .not. use_fates ) then
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
             if (present(c12_soilbiogeochem_carbonstate_inst)) then
                this%cwdc_col(c)    = c12_soilbiogeochem_carbonstate_inst%cwdc_col(c) * ratio
             else
                this%cwdc_col(c)    = 0._r8
             end if
             this%ctrunc_col(c)     = 0._r8
             this%totlitc_col(c)    = 0._r8
             this%totsomc_col(c)    = 0._r8
             this%totlitc_1m_col(c) = 0._r8
             this%totsomc_1m_col(c) = 0._r8
          end if
       end if
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
    integer  :: i,j,k,l,c
    real(r8) :: m                   ! multiplier for the exit_spinup code
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname   ! temporary
    logical  :: readvar
    integer  :: idata
    logical  :: exit_spinup  = .false.
    logical  :: enter_spinup = .false.
    ! flags for comparing the model and restart decomposition cascades
    integer  :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !------------------------------------------------------------------------

    if (carbon_type == 'c12') then

       do k = 1, ndecomp_pools
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
          if (use_vertsoilc) then
             ptr2d => this%decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_cpools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
                  errMsg(sourcefile, __LINE__))
          end if
       end do

       if (use_vertsoilc) then
          ptr2d => this%ctrunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc_vr', xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%ctrunc_vr_col(:,1) ! nlevdecomp = 1; so treat as 1D variable
          call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc', xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
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
          if (use_vertsoilc) then
             ptr2d => this%decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_cpools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col' &
                  // ' with atmospheric c13 value for: '//trim(varname)
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (this%decomp_cpools_vr_col(i,j,k) /= spval .and. .not. isnan(this%decomp_cpools_vr_col(i,j,k)) ) then
                      this%decomp_cpools_vr_col(i,j,k) = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do

       if (use_vertsoilc) then
          ptr2d => this%ctrunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%ctrunc_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13", xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
    end if

    !--------------------------------
    ! C14 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c14' ) then

       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
          if (use_vertsoilc) then
             ptr2d => this%decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_cpools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col with atmospheric c14 value for: '//&
                  trim(varname)
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (this%decomp_cpools_vr_col(i,j,k) /= spval .and. .not. isnan(this%decomp_cpools_vr_col(i,j,k)) ) then
                      this%decomp_cpools_vr_col(i,j,k) = c12_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do

       if (use_vertsoilc) then 
          ptr2d => this%ctrunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%ctrunc_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14", xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if

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
       if ( .not. use_fates ) then
          this%cwdc_col(i)       = value_column
       end if
       this%ctrunc_col(i)     = value_column
       this%totlitc_col(i)    = value_column
       this%totlitc_1m_col(i) = value_column
       this%totsomc_col(i)    = value_column
       this%totsomc_1m_col(i) = value_column
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
       end do
    end do

    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cpools_vr_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_allc, filter_allc)
    !
    ! !DESCRIPTION:
    ! Perform column-level carbon summary calculations
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonstate_type)          :: this
    type(bounds_type)               , intent(in)    :: bounds          
    integer                         , intent(in)    :: num_allc       ! number of columns in allc filter
    integer                         , intent(in)    :: filter_allc(:) ! filter for all active columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l       ! indices
    integer  :: fc            ! filter indices
    real(r8) :: maxdepth      ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    ! vertically integrate each of the decomposing C pools
    do l = 1, ndecomp_pools
       do fc = 1,num_allc
          c = filter_allc(fc)
          this%decomp_cpools_col(c,l) = 0._r8
       end do
    end do
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%decomp_cpools_col(c,l) = &
                  this%decomp_cpools_col(c,l) + &
                  this%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
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

    ! coarse woody debris carbon
    if (.not. use_fates ) then
       do fc = 1,num_allc
          c = filter_allc(fc)
          this%cwdc_col(c) = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_cwd(l) ) then
             do fc = 1,num_allc
                c = filter_allc(fc)
                this%cwdc_col(c) = this%cwdc_col(c) + this%decomp_cpools_col(c,l)
             end do
          end if
       end do
       
    end if

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
