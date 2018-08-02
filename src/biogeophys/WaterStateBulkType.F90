module WaterStateBulkType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water state variables that just apply to bulk
  ! water. Note that this type extends the base waterstate_type, so the full
  ! waterstatebulk_type contains the union of the fields defined here and the fields
  ! defined in waterstate_type.
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
  use WaterStateType , only : waterstate_type
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, extends(waterstate_type), public :: waterstatebulk_type

     real(r8), pointer :: snow_persistence_col   (:)   ! col length of time that ground has had non-zero snow thickness (sec)
     real(r8), pointer :: int_snow_col           (:)   ! col integrated snowfall (mm H2O)


   contains

     procedure          :: InitBulk         
     procedure          :: RestartBulk      
     procedure, private :: InitBulkAllocate 
     procedure, private :: InitBulkHistory  
     procedure, private :: InitBulkCold     

  end type waterstatebulk_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine InitBulk(this, bounds, info, &
       h2osno_input_col, watsat_col, t_soisno_col)

    class(waterstatebulk_type)     :: this
    type(bounds_type) , intent(in) :: bounds
    class(water_info_base_type), intent(in), target :: info
    real(r8)          , intent(in) :: h2osno_input_col(bounds%begc:)
    real(r8)          , intent(in) :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)          , intent(in) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)


    call this%Init(bounds, info, &
       h2osno_input_col, watsat_col, t_soisno_col)

    call this%InitBulkAllocate(bounds) 

    call this%InitBulkHistory(bounds)

    call this%InitBulkCold(bounds, &
       h2osno_input_col)

  end subroutine InitBulk

  !------------------------------------------------------------------------
  subroutine InitBulkAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(waterstatebulk_type) :: this
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

    allocate(this%snow_persistence_col   (begc:endc))                     ; this%snow_persistence_col   (:)   = nan
    allocate(this%int_snow_col           (begc:endc))                     ; this%int_snow_col           (:)   = nan   



  end subroutine InitBulkAllocate

  !------------------------------------------------------------------------
  subroutine InitBulkHistory(this, bounds)
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
    class(waterstatebulk_type) :: this
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


    ! Snow properties - these will be vertically averaged over the snow profile

    this%int_snow_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('INT_SNOW'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('accumulated swe (vegetated landunits only)'), &
         ptr_col=this%int_snow_col, l2g_scale_type='veg', &
         default='inactive')

    call hist_addfld1d ( &
         fname=this%info%fname('INT_SNOW_ICE'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('accumulated swe (ice landunits only)'), &
         ptr_col=this%int_snow_col, l2g_scale_type='ice', &
         default='inactive')

    this%snow_persistence_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('SNOW_PERSISTENCE'),  &
         units='seconds',  &
         avgflag='I', &
         long_name=this%info%lname('Length of time of continuous snow cover (nat. veg. landunits only)'), &
         ptr_col=this%snow_persistence_col, l2g_scale_type='natveg') 


  end subroutine InitBulkHistory

  !-----------------------------------------------------------------------
  subroutine InitBulkCold(this, bounds, &
       h2osno_input_col)
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
    class(waterstatebulk_type)                :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: h2osno_input_col(bounds%begc:)
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

    SHR_ASSERT_ALL((ubound(h2osno_input_col)     == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))

    do c = bounds%begc,bounds%endc
       this%int_snow_col(c)           = h2osno_input_col(c) 
       this%snow_persistence_col(c)   = 0._r8
    end do


    associate(snl => col%snl) 


    end associate

  end subroutine InitBulkCold

  !------------------------------------------------------------------------
  subroutine RestartBulk(this, bounds, ncid, flag, &
       watsat_col)
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
    class(waterstatebulk_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    real(r8)         , intent(in)    :: watsat_col (bounds%begc:, 1:)  ! volumetric soil water at saturation (porosity)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,j
    logical  :: readvar
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(watsat_col) == (/bounds%endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))

    call this%restart (bounds, ncid, flag=flag, &
         watsat_col=watsat_col(bounds%begc:bounds%endc,:)) 


    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('INT_SNOW'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('accuumulated snow'), &
         units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%int_snow_col)
    if (flag=='read' .and. .not. readvar) then
       this%int_snow_col(:) = 0.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('SNOW_PERS'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('continuous snow cover time'), &
         units='sec', &
         interpinic_flag='interp', readvar=readvar, data=this%snow_persistence_col)    
    if (flag=='read' .and. .not. readvar) then
         this%snow_persistence_col(:) = 0.0_r8
    end if

  end subroutine RestartBulk


end module WaterStateBulkType
