module NutrientCompetitionMethodMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Abstract base class for functions to calculate nutrient competition
  !
  ! Created by Jinyun Tang, following Bill Sack's implementation of polymorphism
  ! !USES:
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: nutrient_competition_method_type

  type, abstract :: nutrient_competition_method_type
     private
   contains
     ! initialize the nutrient competition module
     procedure(init_nutrient_competition_interface), deferred :: init_nutrient_competition

     ! compute plant nutrient demand
     procedure(calc_plant_nutrient_demand_interface), deferred :: calc_plant_nutrient_demand

     ! read in nutrient competition kinetic parameters
     procedure(read_nutrient_competition_params_interface), deferred :: read_nutrient_competition_params
     
     ! compute the nutrient yield for different components
     procedure(calc_nutrient_competition_interface), deferred :: calc_nutrient_competition
   
  end type nutrient_competition_method_type

  abstract interface

     ! Note: The following code is adapted based on what Bill Scaks has done for soil water retention curve
     ! polymorphism. Therefore, I also keep some suggestions he gave there.
     !
     ! - Make the interfaces contain all possible inputs that are needed by any
     !   implementation; each implementation will then ignore the inputs it doesn't need.
     !
     ! - For inputs that are needed only by particular implementations - and particularly
     !   for inputs that are constant in time 
     !   pass these into the constructor, and save pointers to these inputs as components
     !   of the child type that needs them. Then they aren't needed as inputs to the
     !   individual routines, allowing the interfaces for these routines to remain more
     !   consistent between different implementations.

     subroutine init_nutrient_competition_interface(this, bounds, crop_prog)
       ! !DESCRIPTION:
       ! initialize nutrient competition calculation
       !
       ! !USES:
       use decompMod           , only : bounds_type       
       import :: nutrient_competition_method_type
       !
       ! !ARGUMENTS:
       class(nutrient_competition_method_type), intent(in) :: this

       type(bounds_type), intent(in) :: bounds    
       logical, intent(in) :: crop_prog


     end subroutine init_nutrient_competition_interface
  !---------------------------------------------------------------------------

     subroutine calc_plant_nutrient_demand_interface(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         photosyns_vars, canopystate_vars, crop_vars, carbonstate_vars, cnstate_vars, &
         carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, c13_carbonflux_vars, c14_carbonflux_vars)
       ! !DESCRIPTION:
       ! Compute plant nutrient demand
       ! right now phosporus is not included by that can be done easily by add another two varaibles such as
       ! phosphorusstate_vars and phosphorusflux_vars
       ! !USES:
       use decompMod           , only : bounds_type  
       use CNCarbonFluxType    , only : carbonflux_type
       use CNCarbonStateType   , only : carbonstate_type
       use CNNitrogenFluxType  , only : nitrogenflux_type
       use CNNitrogenStateType , only : nitrogenstate_type
       use CNStateType         , only : cnstate_type
       use PhotosynthesisType  , only : photosyns_type
       use CropType            , only : crop_type 
       use CanopyStateType     , only : canopystate_type
       import :: nutrient_competition_method_type
       !
       ! !ARGUMENTS:
       class(nutrient_competition_method_type), intent(in) :: this
       type(bounds_type)        , intent(in)    :: bounds
       integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
       integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
       integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
       integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
       type(photosyns_type)     , intent(in)    :: photosyns_vars
       type(crop_type)          , intent(in)    :: crop_vars
       type(canopystate_type)   , intent(in)    :: canopystate_vars
       type(cnstate_type)       , intent(inout) :: cnstate_vars
       type(carbonstate_type)   , intent(in)    :: carbonstate_vars
       type(carbonflux_type)    , intent(inout) :: carbonflux_vars
       type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
       type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
       type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
       type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
       
     end subroutine calc_plant_nutrient_demand_interface
  !---------------------------------------------------------------------------
     subroutine read_nutrient_competition_params_interface(this, ncid)
       ! !DESCRIPTION:
       ! read in kinetic parameters that are needed for doing nutrient competition
       !
       ! !USES:
       use ncdio_pio , only : file_desc_t
       import :: nutrient_competition_method_type
       !
       ! !ARGUMENTS:
       class(nutrient_competition_method_type), intent(in) :: this
       type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
       
     end subroutine read_nutrient_competition_params_interface

  !---------------------------------------------------------------------------     
     subroutine calc_nutrient_competition_interface(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
     photosyns_vars, cnstate_vars, carbonstate_vars, carbonflux_vars, c13_carbonflux_vars, &
     c14_carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars)
     !
     !DESCRIPTION
     !calculate nutrient yield after considering competition between different components
     !
     ! USES
     use decompMod           , only : bounds_type       
     use CNCarbonFluxType    , only : carbonflux_type
     use CNCarbonStateType   , only : carbonstate_type
     use CNNitrogenFluxType  , only : nitrogenflux_type
     use CNNitrogenStateType , only : nitrogenstate_type
     use CNStateType         , only : cnstate_type
     use PhotosynthesisType  , only : photosyns_type
            
     
     import :: nutrient_competition_method_type
     
  ! !ARGUMENTS:
     class(nutrient_competition_method_type), intent(in) :: this
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
     integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
     integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
     integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
     type(photosyns_type)     , intent(in)    :: photosyns_vars   
     type(cnstate_type)       , intent(inout) :: cnstate_vars
     type(carbonstate_type)   , intent(in)    :: carbonstate_vars
     type(carbonflux_type)    , intent(inout) :: carbonflux_vars
     type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
     type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
     type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
     type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
  
     end subroutine calc_nutrient_competition_interface
  end interface

end module NutrientCompetitionMethodMod
