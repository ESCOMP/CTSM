module NutrientCompetitionFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of nutrient_competition_method_type. This module figures
  ! out the particular type to return.
  !
  ! !USES:
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use clm_varctl          , only : iulog

  implicit none
  save
  private
  !
  ! !PUBLIC ROUTINES:
  public :: create_nutrient_competition_method  ! create an object of class nutrient_competition_method_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  function create_nutrient_competition_method(bounds) result(nutrient_competition_method)
    !
    ! !DESCRIPTION:
    ! Create and return an object of nutrient_competition_method_type. The particular type
    ! is determined based on a namelist parameter.
    !
    ! !USES:
    use shr_kind_mod                      , only : SHR_KIND_CL
    use NutrientCompetitionMethodMod      , only : nutrient_competition_method_type
    use NutrientCompetitionCLM45defaultMod, only : nutrient_competition_clm45default_type
    use NutrientCompetitionFlexibleCNMod  , only : nutrient_competition_FlexibleCN_type
    use decompMod                         , only : bounds_type

    ! FIXME(bja, 2015-06) need to pass method control in as a parameter
    ! instead of relying on a global!
    use clm_varctl, only : use_flexibleCN   

    !
    ! !ARGUMENTS:
    class(nutrient_competition_method_type), allocatable :: nutrient_competition_method  ! function result
    type(bounds_type),                       intent(in)  :: bounds
    !
    ! !LOCAL VARIABLES:

    ! For now, hard-code the method. Eventually this will be set from namelist, either by
    ! this routine (appropriate if the 'method' is in its own namelist group), or do the
    ! namelist read outside this module and pass the method in as a parameter (appropriate
    ! if the 'method' is part of a larger namelist group).
    character(len=SHR_KIND_CL) :: method
    
    character(len=*), parameter :: subname = 'create_nutrient_competition_method'
    !-----------------------------------------------------------------------
    
    ! FIXME(bja, 2015-06) flexible_cn may need to be
    ! merged with other nitrogen code, so a more robust method of
    ! selecting the competition method will depend on how the science
    ! is merged.
    method = "clm45default"
    if (use_flexibleCN) then
       method = "flexible_cn"
    end if
    
    select case (trim(method))
       
    case ("clm45default")
       allocate(nutrient_competition_method, &
            source=nutrient_competition_clm45default_type())

    case ("flexible_cn")
       allocate(nutrient_competition_method, &
            source=nutrient_competition_FlexibleCN_type())

    case default
       write(iulog,*) subname//' ERROR: unknown method: ', method
       call endrun(msg=errMsg(sourcefile, __LINE__))

    end select
    call nutrient_competition_method%Init(bounds)

  end function create_nutrient_competition_method

end module NutrientCompetitionFactoryMod
