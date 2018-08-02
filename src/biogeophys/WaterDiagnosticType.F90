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
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use clm_varctl     , only : use_vancouver, use_mexicocity, iulog
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
  type, public :: waterdiagnostic_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: snowice_col            (:)   ! col average snow ice lens
     real(r8), pointer :: snowliq_col            (:)   ! col average snow liquid water

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
  subroutine Init(this, bounds, info)

    class(waterdiagnostic_type)            :: this
    type(bounds_type) , intent(in)    :: bounds  
    class(water_info_base_type), intent(in), target :: info

    this%info => info

    call this%InitAllocate(bounds) 

    call this%InitHistory(bounds)

    call this%InitCold(bounds)

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
    class(waterdiagnostic_type) :: this
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

    allocate(this%snowice_col            (begc:endc))                     ; this%snowice_col            (:)   = nan   
    allocate(this%snowliq_col            (begc:endc))                     ; this%snowliq_col            (:)   = nan   
    allocate(this%h2osoi_liqice_10cm_col (begc:endc))                     ; this%h2osoi_liqice_10cm_col (:)   = nan   
    allocate(this%tws_grc                (begg:endg))                     ; this%tws_grc                (:)   = nan



    allocate(this%qg_snow_col            (begc:endc))                     ; this%qg_snow_col            (:)   = nan   
    allocate(this%qg_soil_col            (begc:endc))                     ; this%qg_soil_col            (:)   = nan   
    allocate(this%qg_h2osfc_col          (begc:endc))                     ; this%qg_h2osfc_col          (:)   = nan   
    allocate(this%qg_col                 (begc:endc))                     ; this%qg_col                 (:)   = nan   
    allocate(this%qaf_lun                (begl:endl))                     ; this%qaf_lun                (:)   = nan
    allocate(this%q_ref2m_patch          (begp:endp))                     ; this%q_ref2m_patch          (:)   = nan


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
    class(waterdiagnostic_type) :: this
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
    use shr_const_mod   , only : shr_const_pi
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use shr_spfn_mod    , only : shr_spfn_erf
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use landunit_varcon , only : istwet, istsoil, istdlak, istcrop, istice_mec  
    use column_varcon   , only : icol_shadewall, icol_road_perv
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use clm_varcon      , only : spval, sb, bdsno 
    use clm_varcon      , only : zlnd, tfrz, spval, pc
    use clm_varctl      , only : fsurdat, iulog
    use clm_varctl        , only : use_bedrock
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    class(waterdiagnostic_type)                :: this
    type(bounds_type)     , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev
    real(r8)           :: maxslope, slopemax, minslope
    real(r8)           :: d, fd, dfdd, slope0,slopebeta
    real(r8) ,pointer  :: std (:)     
    logical            :: readvar 
    type(file_desc_t)  :: ncid        
    character(len=256) :: locfn       
    !-----------------------------------------------------------------------



    do l = bounds%begl, bounds%endl 
       if (lun%urbpoi(l)) then
          if (use_vancouver) then
             this%qaf_lun(l) = 0.0111_r8
          else if (use_mexicocity) then
             this%qaf_lun(l) = 0.00248_r8
          else
             this%qaf_lun(l) = 1.e-4_r8 ! Arbitrary set since forc_q is not yet available
          end if
       end if
    end do


    associate(snl => col%snl) 



    end associate

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod          , only : masterproc
    use clm_varcon       , only : pondmx, watmin, spval, nameg
    use landunit_varcon  , only : istcrop, istdlak, istsoil  
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager , only : is_first_step
    use clm_varctl       , only : bound_h2osoi
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterdiagnostic_type) :: this
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
