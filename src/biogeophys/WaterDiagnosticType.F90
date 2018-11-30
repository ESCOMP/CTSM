module WaterDiagnosticType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water diagnostic variables that apply to both bulk
  ! water and water tracers. Diagnostic variables are neither fundamental state variables
  ! nor fluxes between those fundamental states, but are typically derived from those
  ! states and/or fluxes.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use decompMod      , only : BOUNDS_SUBGRID_PATCH, BOUNDS_SUBGRID_COLUMN, BOUNDS_SUBGRID_LANDUNIT, BOUNDS_SUBGRID_GRIDCELL
  use clm_varctl     , only : use_vancouver, use_mexicocity
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  use WaterTracerUtils, only : AllocateVar1d
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterdiagnostic_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: snowice_col            (:)   ! col average snow ice lens
     real(r8), pointer :: snowliq_col            (:)   ! col average snow liquid water

     real(r8), pointer :: total_plant_stored_h2o_col(:) ! col water that is bound in plants, including roots, sapwood, leaves, etc
                                                        ! in most cases, the vegetation scheme does not have a dynamic
                                                        ! water storage in plants, and thus 0.0 is a suitable for the trivial case.
                                                        ! When FATES is coupled in with plant hydraulics turned on, this storage
                                                        ! term is set to non-zero. (kg/m2 H2O)

     real(r8), pointer :: h2osoi_liqice_10cm_col (:)   ! col liquid water + ice lens in top 10cm of soil (kg/m2)
     real(r8), pointer :: tws_grc                (:)   ! grc total water storage (mm H2O)
     real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)
     real(r8), pointer :: qg_snow_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_soil_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_h2osfc_col          (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_col                 (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qaf_lun                (:)   ! lun urban canopy air specific humidity (kg/kg)

   contains

     procedure          :: Init         
     procedure          :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type waterdiagnostic_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info, tracer_vars)

    class(waterdiagnostic_type), intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds  
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: tracer_vars

    this%info => info

    call this%InitAllocate(bounds, tracer_vars)

    call this%InitHistory(bounds)

    call this%InitCold(bounds)

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
    class(waterdiagnostic_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  
    type(water_tracer_container_type), intent(inout) :: tracer_vars
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------

    call AllocateVar1d(var = this%snowice_col, name = 'snowice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%snowliq_col, name = 'snowliq_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%total_plant_stored_h2o_col, name = 'total_plant_stored_h2o_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%h2osoi_liqice_10cm_col, name = 'h2osoi_liqice_10cm_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%tws_grc, name = 'tws_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_GRIDCELL)
    call AllocateVar1d(var = this%qg_snow_col, name = 'qg_snow_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qg_soil_col, name = 'qg_soil_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qg_h2osfc_col, name = 'qg_h2osfc_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qg_col, name = 'qg_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_COLUMN)
    call AllocateVar1d(var = this%qaf_lun, name = 'qaf_lun', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_LANDUNIT)
    call AllocateVar1d(var = this%q_ref2m_patch, name = 'q_ref2m_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = BOUNDS_SUBGRID_PATCH)

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
    class(waterdiagnostic_type), intent(in) :: this
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

    this%h2osoi_liqice_10cm_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SOILWATER_10CM'),  &
         units='kg/m2', &
         avgflag='A', &
         long_name=this%info%lname('soil liquid water + ice in top 10cm of soil (veg landunits only)'), &
         ptr_col=this%h2osoi_liqice_10cm_col, set_urb=spval, set_lake=spval, l2g_scale_type='veg')

    this%tws_grc(begg:endg) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('TWS'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('total water storage'), &
         ptr_lnd=this%tws_grc)

    ! (rgk 02-02-2017) There is intentionally no entry  here for stored plant water
    !                  I think that since the value is zero in all cases except
    !                  for FATES plant hydraulics, it will be confusing for users
    !                  when they see their plants have no water in output files.
    !                  So it is not useful diagnostic information. The information
    !                  can be provided through FATES specific history diagnostics
    !                  if need be.


    this%q_ref2m_patch(begp:endp) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('Q2M'), &
         units='kg/kg',  &
         avgflag='A', &
         long_name=this%info%lname('2m specific humidity'), &
         ptr_patch=this%q_ref2m_patch)



    ! Snow properties - these will be vertically averaged over the snow profile

    this%snowliq_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOWLIQ'),  &
         units='kg/m2',  &
         avgflag='A', &
         long_name=this%info%lname('snow liquid water'), &
         ptr_col=this%snowliq_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('SNOWLIQ_ICE'),  &
         units='kg/m2',  &
         avgflag='A', &
         long_name=this%info%lname('snow liquid water (ice landunits only)'), &
         ptr_col=this%snowliq_col, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%snowice_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOWICE'),  &
         units='kg/m2', &
         avgflag='A', &
         long_name=this%info%lname('snow ice'), &
         ptr_col=this%snowice_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('SNOWICE_ICE'),  &
         units='kg/m2', &
         avgflag='A', &
         long_name=this%info%lname('snow ice (ice landunits only)'), &
         ptr_col=this%snowice_col, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')


  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    use ncdio_pio       , only : file_desc_t
    !
    ! !ARGUMENTS:
    class(waterdiagnostic_type), intent(in) :: this
    type(bounds_type)     , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: l
    real(r8)           :: ratio
    !-----------------------------------------------------------------------

    ratio = this%info%get_ratio()

    ! Water Stored in plants is almost always a static entity, with the exception
    ! of when FATES-hydraulics is used. As such, this is trivially set to 0.0 (rgk 03-2017)
    this%total_plant_stored_h2o_col(bounds%begc:bounds%endc) = 0.0_r8


    do l = bounds%begl, bounds%endl 
       if (lun%urbpoi(l)) then
          if (use_vancouver) then
             this%qaf_lun(l) = 0.0111_r8 * ratio
          else if (use_mexicocity) then
             this%qaf_lun(l) = 0.00248_r8 * ratio
          else
             this%qaf_lun(l) = 1.e-4_r8 * ratio ! Arbitrary set since forc_q is not yet available
          end if
       end if
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use clm_varcon       , only : nameg
    use ncdio_pio        , only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterdiagnostic_type), intent(in) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,j
    logical  :: readvar
    !------------------------------------------------------------------------



    ! TWS is needed when methane is on and the TWS_inversion is used to get exact
    ! restart.
    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('TWS'), &
         xtype=ncd_double,  &
         dim1name=nameg, &
         long_name=this%info%lname('Total Water Storage'), &
         units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%tws_grc)


    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('qaf'), &
         xtype=ncd_double, dim1name='landunit', &
         long_name=this%info%lname('urban canopy specific humidity'), &
         units='kg/kg', &
         interpinic_flag='interp', readvar=readvar, data=this%qaf_lun)

  end subroutine Restart

end module WaterDiagnosticType
