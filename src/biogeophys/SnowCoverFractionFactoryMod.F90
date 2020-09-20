module SnowCoverFractionFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of snow_cover_fraction_base_type. This module figures
  ! out the particular type to return.
  !
  ! !USES:
  use shr_log_mod                             , only : errMsg => shr_log_errMsg
  use decompMod                               , only : bounds_type
  use ColumnType                              , only : column_type
  use glcBehaviorMod                          , only : glc_behavior_type
  use ncdio_pio                               , only : file_desc_t
  use clm_varctl                              , only : iulog, use_subgrid_fluxes
  use abortutils                              , only : endrun
  use SnowCoverFractionBaseMod                , only : snow_cover_fraction_base_type
  use SnowCoverFractionNiuYang2007Mod         , only : snow_cover_fraction_niu_yang_2007_type
  use SnowCoverFractionSwensonLawrence2012Mod , only : snow_cover_fraction_swenson_lawrence_2012_type
  implicit none
  save
  private
  !
  ! !PUBLIC ROUTINES:
  public :: CreateAndInitSnowCoverFraction

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
contains

  !-----------------------------------------------------------------------
  function CreateAndInitSnowCoverFraction(snow_cover_fraction_method, &
       bounds, col, glc_behavior, NLFilename, params_ncid) &
       result(scf_method)
    !
    ! !DESCRIPTION:
    ! Create an instance of the appropriate snow_cover_fraction_base_type (based on
    ! snow_cover_fraction_method), initialize it and return it
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_base_type), allocatable :: scf_method  ! function result

    character(len=*)        , intent(in)    :: snow_cover_fraction_method
    type(bounds_type)       , intent(in)    :: bounds
    type(column_type)       , intent(in)    :: col
    type(glc_behavior_type) , intent(in)    :: glc_behavior
    character(len=*)        , intent(in)    :: NLFilename ! Namelist filename
    type(file_desc_t)       , intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CreateAndInitSnowCoverFraction'
    !-----------------------------------------------------------------------

    select case (snow_cover_fraction_method)
    case ('SwensonLawrence2012')
       allocate(scf_method, &
            source = snow_cover_fraction_swenson_lawrence_2012_type())

       ! NOTE(wjs, 2019-08-05) This would be more straightforward if we could just do all
       ! initialization in the constructor: then we could avoid an Init method, and the
       ! consequent need for a 'select type' construct. But we have had some compilers
       ! (e.g., intel17) that don't seem to work right with function constructors; Init
       ! methods seem to be more reliable.
       select type (scf_method)
       class is (snow_cover_fraction_swenson_lawrence_2012_type)
          call scf_method%Init( &
               bounds       = bounds, &
               col          = col, &
               glc_behavior = glc_behavior, &
               NLFilename   = NLFilename, &
               params_ncid  = params_ncid)
       class default
          call endrun(msg = "Unexpected class", &
               additional_msg = errMsg(sourcefile, __LINE__))
       end select

    case ('NiuYang2007')
       allocate(scf_method, &
            source = snow_cover_fraction_niu_yang_2007_type())

       ! NOTE(wjs, 2019-08-05) This would be more straightforward if we could just do all
       ! initialization in the constructor: then we could avoid an Init method, and the
       ! consequent need for a 'select type' construct. But we have had some compilers
       ! (e.g., intel17) that don't seem to work right with function constructors; Init
       ! methods seem to be more reliable.
       select type (scf_method)
       class is (snow_cover_fraction_niu_yang_2007_type)
          call scf_method%Init( &
               params_ncid = params_ncid)
       class default
          call endrun(msg = "Unexpected class", &
               additional_msg = errMsg(sourcefile, __LINE__))
       end select

    case default
       write(iulog,*) subname//' ERROR: unknown snow_cover_fraction_method: ', &
            snow_cover_fraction_method
       call endrun(msg = 'unknown snow_cover_fraction_method', &
            additional_msg = errMsg(sourcefile, __LINE__))
    end select

  end function CreateAndInitSnowCoverFraction

end module SnowCoverFractionFactoryMod
