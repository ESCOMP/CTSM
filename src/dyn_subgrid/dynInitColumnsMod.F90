module dynInitColumnsMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle initialization of columns that just switched from inactive to active
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use decompMod            , only : bounds_type
  use abortutils           , only : endrun
  use clm_varctl           , only : iulog  
  use clm_varcon           , only : namec
  use TemperatureType      , only : temperature_type
  use WaterStateBulkType       , only : waterstatebulk_type
  use SoilHydrologyType    , only : soilhydrology_type
  use GridcellType         , only : grc
  use LandunitType         , only : lun
  use ColumnType           , only : col
  use dynColumnTemplateMod , only : template_col_from_landunit, TEMPLATE_NONE_FOUND
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  ! The following is the public interface to the routines in this module:
  public :: initialize_new_columns  ! Do initialization for all columns that are newly-active in this time step
  
  ! The following are public only for unit testing purposes, and should not be called
  ! directly by application code:
  public :: initial_template_col_crop ! Find column to use as a template for a crop column that has newly become active
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: initial_template_col_dispatcher ! Find column to use as a template; dispatcher to the appropriate routine based on landunit type
  private :: initial_template_col_soil       ! Find column to use as a template for a vegetated column that has newly become active
  private :: copy_state                      ! Copy a subset of state variables from template column to newly-active column

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initialize_new_columns(bounds, cactive_prior, &
       temperature_inst, waterstatebulk_inst, soilhydrology_inst)
    !
    ! !DESCRIPTION:
    ! Do initialization for all columns that are newly-active in this time step
    !
    ! !USES:
    use GetGlobalValuesMod , only : GetGlobalWrite
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                        ! bounds
    logical                , intent(in)    :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    type(temperature_type)   , intent(inout) :: temperature_inst
    type(waterstatebulk_type)    , intent(inout) :: waterstatebulk_inst
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c          ! column index
    integer :: c_template ! index of template column

    character(len=*), parameter :: subname = 'initialize_new_columns'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_ALL((ubound(cactive_prior) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    do c = bounds%begc, bounds%endc
       ! If this column is newly-active, then we need to initialize it using the routines in this module
       if (col%active(c) .and. .not. cactive_prior(c)) then
          c_template = initial_template_col_dispatcher(bounds, c, cactive_prior(bounds%begc:bounds%endc))
          if (c_template /= TEMPLATE_NONE_FOUND) then
             call copy_state(c, c_template, &
                  temperature_inst, waterstatebulk_inst, soilhydrology_inst)
          else
             write(iulog,*) subname// ' WARNING: No template column found to initialize newly-active column'
             write(iulog,*) '-- keeping the state that was already in memory, possibly from arbitrary initialization'
             call GetGlobalWrite(decomp_index=c, clmlevel=namec)
          end if
       end if
    end do

  end subroutine initialize_new_columns


  !-----------------------------------------------------------------------
  function initial_template_col_dispatcher(bounds, c_new, cactive_prior) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for the given column that has newly become active;
    ! this is a dispatcher that calls the appropriate routine based on the landunit type of c_new.
    !
    ! Returns TEMPLATE_NONE_FOUND if there is no column to use for initialization
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop, istice_mec, istdlak, istwet, isturb_MIN, isturb_MAX
    !
    ! !ARGUMENTS:
    integer :: c_template  ! function result
    type(bounds_type) , intent(in) :: bounds                        ! bounds
    integer           , intent(in) :: c_new                         ! column index that needs initialization
    logical           , intent(in) :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    !
    ! !LOCAL VARIABLES:
    integer :: l     ! landunit index
    integer :: ltype ! landunit type

    character(len=*), parameter :: subname = 'initial_template_col_dispatcher'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_ALL((ubound(cactive_prior) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    l = col%landunit(c_new)
    ltype = lun%itype(l)
    select case(ltype)
    case(istsoil)
       c_template = initial_template_col_soil(c_new)
    case(istcrop)
       c_template = initial_template_col_crop(bounds, c_new, cactive_prior(bounds%begc:bounds%endc))
    case(istice_mec)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active glacier mec column not yet implemented'
       write(iulog,*) 'Expectation is that glacier mec columns should be active from the start of the run wherever they can grow'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(sourcefile, __LINE__))
    case(istdlak)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active lake column not yet implemented'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(sourcefile, __LINE__))
    case(istwet)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active wetland column not yet implemented'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(sourcefile, __LINE__))
    case(isturb_MIN:isturb_MAX)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active urban column not yet implemented'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(sourcefile, __LINE__))
    case default
       write(iulog,*) subname// ' ERROR: Unknown landunit type: ', ltype
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(sourcefile, __LINE__))
    end select

  end function initial_template_col_dispatcher


  !-----------------------------------------------------------------------
  function initial_template_col_soil(c_new) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for a vegetated column that has newly become active.
    !
    ! For now, we assume that the only vegetated columns that can newly become active are
    ! ones with 0 weight on the grid cell (i.e., virtual columns). For these, we simply
    ! keep the state at the current value (likely arbitrary initial conditions), and so
    ! return TEMPLATE_NONE_FOUND from this function. Within this function, we check this assumption.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer              :: c_template ! function result
    integer , intent(in) :: c_new        ! column index that needs initialization
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initial_template_col_soil'
    !-----------------------------------------------------------------------

    if (col%wtgcell(c_new) > 0._r8) then
       write(iulog,*) subname// ' ERROR: Expectation is that the only vegetated columns that&
            & can newly become active are ones with 0 weight on the grid cell'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(sourcefile, __LINE__))
    end if

    c_template = TEMPLATE_NONE_FOUND
    
  end function initial_template_col_soil

  !-----------------------------------------------------------------------
  function initial_template_col_crop(bounds, c_new, cactive_prior) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for a crop column that has newly become active
    !
    ! Returns TEMPLATE_NONE_FOUND if there is no column to use for initialization
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    integer :: c_template  ! function result
    type(bounds_type) , intent(in) :: bounds                        ! bounds
    integer           , intent(in) :: c_new                         ! column index that needs initialization
    logical           , intent(in) :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initial_template_col_crop'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(cactive_prior) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    
    ! First try to find an active column on the vegetated landunit; if there is none, then
    ! find the first active column on the crop landunit; if there is none, then
    ! template_col will be TEMPLATE_NONE_FOUND
    c_template = template_col_from_landunit(bounds, c_new, istsoil, cactive_prior(bounds%begc:bounds%endc))
    if (c_template == TEMPLATE_NONE_FOUND) then
       c_template = template_col_from_landunit(bounds, c_new, istcrop, cactive_prior(bounds%begc:bounds%endc))
    end if

  end function initial_template_col_crop


  !-----------------------------------------------------------------------
  subroutine copy_state(c_new, c_template, &
       temperature_inst, waterstatebulk_inst, soilhydrology_inst)
    !
    ! !DESCRIPTION:
    ! Copy a subset of state variables from a template column (c_template) to a newly-
    ! active column (c_new)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer, intent(in) :: c_new      ! index of newly-active column
    integer, intent(in) :: c_template ! index of column to use as a template
    type(temperature_type)  , intent(inout) :: temperature_inst
    type(waterstatebulk_type)   , intent(inout) :: waterstatebulk_inst
    type(soilhydrology_type), intent(inout) :: soilhydrology_inst
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'copy_state'
    !-----------------------------------------------------------------------

    ! For now, just copy a few key variables
    ! TODO(wjs, 2016-08-31) Figure out what else should be copied here

    ! We only copy the below-ground portion of these multi-level variables, not the
    ! above-ground (snow) portion. This is because it is challenging to initialize the
    ! snow pack in a consistent state, requiring copying many more state variables - and
    ! if you initialize it in a partly-inconsistent state, you get balance errors. So, for
    ! now at least, we (Dave Lawrence, Keith Oleson, Bill Sacks) have decided that it's
    ! safest to just let the snow pack in the new column start at cold start conditions.

    temperature_inst%t_soisno_col(c_new,1:) = temperature_inst%t_soisno_col(c_template,1:)

    ! TODO(wjs, 2016-08-31) If we had more general uses of this initial template col
    ! infrastructure (copying state between very different landunits), then we might need
    ! to handle bedrock layers - e.g., zeroing out any water that would be added to a
    ! bedrock layer(?). But for now we just use this initial template col infrastructure
    ! for nat veg -> crop, for which the bedrock will be the same, so we're not dealing
    ! with that complexity for now.
    waterstatebulk_inst%h2osoi_liq_col(c_new,1:) = waterstatebulk_inst%h2osoi_liq_col(c_template,1:)
    waterstatebulk_inst%h2osoi_ice_col(c_new,1:) = waterstatebulk_inst%h2osoi_ice_col(c_template,1:)
    waterstatebulk_inst%h2osoi_vol_col(c_new,1:) = waterstatebulk_inst%h2osoi_vol_col(c_template,1:)

    soilhydrology_inst%wa_col(c_new) = soilhydrology_inst%wa_col(c_template)

  end subroutine copy_state



end module dynInitColumnsMod
