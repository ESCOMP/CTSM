module dynurbanFileMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the dataset that specifies transient areas of the urban landunits
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
  use landunit_varcon       , only : numurbl
  use UrbanParamsType       , only : CheckUrban
  
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynurban_init    ! initialize information read from landuse.timeseries dataset
  public :: dynurban_interp  ! get landuse data for the current time step, if needed
  !
  ! ! PRIVATE TYPES
  type(dyn_file_type), target      :: dynurban_file ! information for the file containing transient urban data
  type(dyn_var_time_uninterp_type) :: wturban       ! weight of the urban landunits

  ! Names of variables on file
  character(len=*), parameter :: urban_varname = 'PCT_URBAN'

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynurban_init(bounds, dynurban_filename)
    !
    ! !DESCRIPTION:
    ! Initialize dataset containing transient urban info (position it to the right time
    ! samples that bound the initial model date)
    !
    ! !USES:
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    use ncdio_pio      , only : check_dim_size
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds            ! proc-level bounds
    character(len=*)  , intent(in) :: dynurban_filename ! name of file containing transient urban information
    !
    ! !LOCAL VARIABLES:
    integer :: num_points     ! number of spatial points
    integer :: wturb_shape(2) ! shape of the wturb data
	 
    character(len=*), parameter :: subname = 'dynurban_init'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT(bounds%level == bounds_level_proc, subname // ': argument must be PROC-level bounds')

    if (masterproc) then
       write(iulog,*) 'Attempting to read urban dynamic landuse data .....'
    end if

    ! Get the year from the START of the timestep; this way, we'll update urban areas
    ! starting after the year boundary. This is consistent with the timing of glacier
    ! updates, and will likely be consistent with the timing of urban updates determined
    ! prognostically, if urban areas are ever determined prognostically rather than
    ! prescribed ahead of time.
    dynurban_file = dyn_file_type(dynurban_filename, YEAR_POSITION_START_OF_TIMESTEP)
    call check_dim_size(dynurban_file, 'numurbl', numurbl)

    ! read data PCT_URBAN
    !
    ! Note: if you want to change transient urban landunits so that they are interpolated, rather
    ! than jumping to each year's value on Jan 1 of that year, simply change wturban
    ! to be of type dyn_var_time_interp_type (rather than
    ! dyn_var_time_uninterp_type), and change the following constructors to construct
    ! variables of dyn_var_time_interp_type. That's all you need to do.
    num_points = (bounds%endg - bounds%begg + 1)
    wturb_shape = [num_points, numurbl]
    wturban = dyn_var_time_uninterp_type( &
          dyn_file = dynurban_file, varname=urban_varname, &
          dim1name=grlnd, conversion_factor=100._r8, &
          do_check_sums_equal_1=.false., data_shape=wturb_shape)

 
  end subroutine dynurban_init
  
  
  !-----------------------------------------------------------------------
  subroutine dynurban_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get urban cover for model time, when needed.
    !
    ! Sets lun%wtgcell for urban landunits.
    !
    ! Note that urban cover currently jumps to its new value at the start of the year.
    ! However, as mentioned above, this behavior can be changed to time interpolation
    ! simply by making wturban and wtcft dyn_var_time_interp_type variables rather than
    ! dyn_var_time_uninterp_type. 
    !
    ! !USES:
    use landunit_varcon   , only : isturb_tbd, isturb_hd, isturb_md
    use subgridWeightsMod , only : set_landunit_weight
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: g              ! indices
    real(r8), allocatable :: wturban_cur(:,:) ! current weight of the urban landunits
    
    character(len=*), parameter :: subname = 'dynurban_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == bounds_level_proc, subname // ': argument must be PROC-level bounds')

    call dynurban_file%time_info%set_current_year()

    ! Set new landunit areas
    allocate(wturban_cur(bounds%begg:bounds%endg, numurbl))
    call wturban%get_current_data(wturban_cur)
    do g = bounds%begg, bounds%endg
       call set_landunit_weight(g, isturb_tbd, wturban_cur(g,1))
       call set_landunit_weight(g, isturb_hd, wturban_cur(g,2))
       call set_landunit_weight(g, isturb_md, wturban_cur(g,3))
    end do

    ! Check that urban data is valid
    call CheckUrban(bounds%begg, bounds%endg, wturban_cur(bounds%begg:bounds%endg,:), subname)

    deallocate(wturban_cur)

  end subroutine dynurban_interp

end module dynurbanFileMod