module WaterBalanceType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water balance-related variables that apply to both
  ! bulk water and water tracers.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use clm_varctl     , only : iulog
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterbalance_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: h2osno_old_col         (:)   ! col snow mass for previous time step (kg/m2) (new)
     real(r8), pointer :: liq1_grc               (:)   ! grc initial gridcell total h2o liq content
     real(r8), pointer :: liq2_grc               (:)   ! grc post land cover change total liq content
     real(r8), pointer :: ice1_grc               (:)   ! grc initial gridcell total h2o ice content
     real(r8), pointer :: ice2_grc               (:)   ! grc post land cover change total ice content


     ! Balance Checks

     real(r8), pointer :: begwb_col              (:)   ! water mass begining of the time step
     real(r8), pointer :: endwb_col              (:)   ! water mass end of the time step
     real(r8), pointer :: errh2o_patch           (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2o_col             (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2osno_col          (:)   ! snow water conservation error(mm H2O)

   contains

     procedure          :: Init         
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  

  end type waterbalance_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info)

    class(waterbalance_type)            :: this
    type(bounds_type) , intent(in)    :: bounds  
    class(water_info_base_type), intent(in), target :: info

    this%info => info

    call this%InitAllocate(bounds) 

    call this%InitHistory(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(waterbalance_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate(this%h2osno_old_col         (begc:endc))                     ; this%h2osno_old_col         (:)   = nan   
    allocate(this%liq1_grc               (begg:endg))                     ; this%liq1_grc               (:)   = nan
    allocate(this%liq2_grc               (begg:endg))                     ; this%liq2_grc               (:)   = nan
    allocate(this%ice1_grc               (begg:endg))                     ; this%ice1_grc               (:)   = nan
    allocate(this%ice2_grc               (begg:endg))                     ; this%ice2_grc               (:)   = nan




    allocate(this%begwb_col              (begc:endc))                     ; this%begwb_col              (:)   = nan
    allocate(this%endwb_col              (begc:endc))                     ; this%endwb_col              (:)   = nan
    allocate(this%errh2o_patch           (begp:endp))                     ; this%errh2o_patch           (:)   = nan
    allocate(this%errh2o_col             (begc:endc))                     ; this%errh2o_col             (:)   = nan
    allocate(this%errh2osno_col          (begc:endc))                     ; this%errh2osno_col          (:)   = nan
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_lch4
    use clm_varctl     , only : hist_wrtch4diag
    use clm_varpar     , only : nlevsno, nlevsoi
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal, no_snow_zero
    !
    ! !ARGUMENTS:
    class(waterbalance_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    this%liq1_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('LIQUID_CONTENT1'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('initial gridcell total liq content'), &
         ptr_lnd=this%liq1_grc)

    this%liq2_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('LIQUID_CONTENT2'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('post landuse change gridcell total liq content'), &
         ptr_lnd=this%liq2_grc, default='inactive')     

    this%ice1_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('ICE_CONTENT1'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('initial gridcell total ice content'), &
         ptr_lnd=this%ice1_grc)     

    this%ice2_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('ICE_CONTENT2'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('post land cover change total ice content'), &
         ptr_lnd=this%ice2_grc, default='inactive')


    this%errh2o_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('ERRH2O'), &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('total water conservation error'), &
         ptr_col=this%errh2o_col)

    this%errh2osno_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('ERRH2OSNO'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('imbalance in snow depth (liquid water)'), &
         ptr_col=this%errh2osno_col, c2l_scale_type='urbanf')
  end subroutine InitHistory

  !-----------------------------------------------------------------------


end module WaterBalanceType
