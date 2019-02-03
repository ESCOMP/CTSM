module WaterBalanceType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water balance-related variables that apply to both
  ! bulk water and water tracers.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use decompMod      , only : BOUNDS_SUBGRID_PATCH, BOUNDS_SUBGRID_COLUMN, BOUNDS_SUBGRID_GRIDCELL
  use clm_varcon     , only : spval
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  use WaterTracerUtils, only : AllocateVar1d
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

     real(r8), pointer :: snow_sources_col         (:)   ! col snow sources (mm H2O/s)
     real(r8), pointer :: snow_sinks_col           (:)   ! col snow sinks (mm H2O/s)

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
  subroutine Init(this, bounds, info, tracer_vars)

    class(waterbalance_type), intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds  
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: tracer_vars

    this%info => info

    call this%InitAllocate(bounds, tracer_vars)

    call this%InitHistory(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, tracer_vars)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterbalance_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  
    type(water_tracer_container_type), intent(inout) :: tracer_vars
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------

    call AllocateVar1d(var = this%h2osno_old_col, name = 'h2osno_old_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%liq1_grc, name = 'liq1_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL)
    call AllocateVar1d(var = this%liq2_grc, name = 'liq2_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL)
    call AllocateVar1d(var = this%ice1_grc, name = 'ice1_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL)
    call AllocateVar1d(var = this%ice2_grc, name = 'ice2_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL)

    call AllocateVar1d(var = this%snow_sources_col, name = 'snow_sources_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%snow_sinks_col, name = 'snow_sinks_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)

    call AllocateVar1d(var = this%begwb_col, name = 'begwb_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%endwb_col, name = 'endwb_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%errh2o_patch, name = 'errh2o_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)
    call AllocateVar1d(var = this%errh2o_col, name = 'errh2o_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%errh2osno_col, name = 'errh2osno_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(waterbalance_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at any
    ! given time step but only if there is at least one snow layer (for all landunits 
    ! except lakes).  Also note that monthly average files of snow_sources and snow_sinks
    ! sinks must be weighted by number of days in the month to diagnose, for example, an 
    ! annual value of the change in h2osno. 

    this%snow_sources_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOW_SOURCES'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('snow sources (liquid water)'), &
         ptr_col=this%snow_sources_col, c2l_scale_type='urbanf')

    this%snow_sinks_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOW_SINKS'),  &
         units='mm/s',  &
         avgflag='A', &
         long_name=this%info%lname('snow sinks (liquid water)'), &
         ptr_col=this%snow_sinks_col, c2l_scale_type='urbanf')

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

end module WaterBalanceType
