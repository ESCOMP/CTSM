module dynlakeFileMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the dataset that specifies transient areas of the lake landunit 
  !
  ! !USES:
  
#include "shr_assert.h"
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type, bounds_level_proc
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use clm_varctl            , only : iulog
  use clm_varcon            , only : grlnd
  use abortutils            , only : endrun
  use spmdMod               , only : masterproc
  
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynlake_init     ! initialize information read from landuse.timeseries dataset
  public :: dynlake_interp   ! get landuse data for the current time step, if needed
  !
  ! ! PRIVATE TYPES
  type(dyn_file_type), target      :: dynlake_file ! information for the file containing transient lake data
  type(dyn_var_time_uninterp_type) :: wtlake       ! weight of the lake landunit

  ! Names of variables on file
  character(len=*), parameter :: lake_varname = 'PCT_LAKE'

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynlake_init(bounds, dynlake_filename)
    !
    ! !DESCRIPTION:
    ! Initialize dataset containing transient lake info (position it to the right time
    ! samples that bound the initial model date)
    !
    ! !USES:
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds           ! proc-level bounds
    character(len=*)  , intent(in) :: dynlake_filename ! name of file containing transient lake information
    !
    ! !LOCAL VARIABLES:
    integer :: num_points     ! number of spatial points
	 
    character(len=*), parameter :: subname = 'dynlake_init'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT(bounds%level == bounds_level_proc, subname // ': argument must be PROC-level bounds')

    if (masterproc) then
       write(iulog,*) 'Attempting to read lake dynamic landuse data .....'
    end if

    ! Get the year from the START of the timestep; this way, we'll update lake areas
    ! starting after the year boundary. This is consistent with the timing of glacier
    ! updates, and will likely be consistent with the timing of lake updates determined
    ! prognostically, if lake areas are ever determined prognostically rather than
    ! prescribed ahead of time.
    dynlake_file = dyn_file_type(dynlake_filename, YEAR_POSITION_START_OF_TIMESTEP)

    ! read data PCT_LAKE 
    !
    ! Note: if you want to change transient lakes so that they are interpolated, rather
    ! than jumping to each year's value on Jan 1 of that year, simply change wtlake and
    ! to be of type dyn_var_time_interp_type (rather than
    ! dyn_var_time_uninterp_type), and change the following constructors to construct
    ! variables of dyn_var_time_interp_type. That's all you need to do.
    num_points = (bounds%endg - bounds%begg + 1)
    wtlake = dyn_var_time_uninterp_type( &
         dyn_file = dynlake_file, varname=lake_varname, &
         dim1name=grlnd, conversion_factor=100._r8, &
         do_check_sums_equal_1=.false., data_shape=[num_points])
 
  end subroutine dynlake_init
  
  
  !-----------------------------------------------------------------------
  subroutine dynlake_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get lake cover for model time, when needed.
    !
    ! Sets col%wtlunit and lun%wtgcell for lake landunits.
    !
    ! Note that lake cover currently jumps to its new value at the start of the year.
    ! However, as mentioned above, this behavior can be changed to time interpolation
    ! simply by making wtlake and wtcft dyn_var_time_interp_type variables rather than
    ! dyn_var_time_uninterp_type. 
    !
    ! !USES:
    use landunit_varcon   , only : istdlak
    use subgridWeightsMod , only : set_landunit_weight
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: g              ! indices
    real(r8), allocatable :: wtlake_cur(:)  ! current weight of the lake landunit
    
    character(len=*), parameter :: subname = 'dynlake_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == bounds_level_proc, subname // ': argument must be PROC-level bounds')

    call dynlake_file%time_info%set_current_year()

    ! Set new landunit area
    allocate(wtlake_cur(bounds%begg:bounds%endg))
    call wtlake%get_current_data(wtlake_cur)
    do g = bounds%begg, bounds%endg
       call set_landunit_weight(g, istdlak, wtlake_cur(g))
    end do
    deallocate(wtlake_cur)

  end subroutine dynlake_interp

end module dynlakeFileMod
