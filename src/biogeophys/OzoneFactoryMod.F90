module OzoneFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of ozone_base_type. This module figures out the
  ! particular type to return.
  !
  ! !USES:
  use decompMod   , only : bounds_type

  implicit none
  save
  private

  !
  ! !PUBLIC ROUTINES:
  public :: create_and_init_ozone_type  ! create an object of class ozone_base_type

contains

  !-----------------------------------------------------------------------
  function create_and_init_ozone_type(bounds) result(ozone)
    !
    ! !DESCRIPTION:
    ! Create and initialize an object of ozone_base_type, and return this object. The
    ! particular type is determined based on the o3_veg_stress_method namelist parameter.
    !
    ! !USES:
    use clm_varctl   , only : o3_veg_stress_method
    use OzoneBaseMod , only : ozone_base_type
    use OzoneOffMod  , only : ozone_off_type
    use OzoneMod     , only : ozone_type
    !
    ! !ARGUMENTS:
    class(ozone_base_type), allocatable :: ozone  ! function result
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'create_and_init_ozone_type'
    !-----------------------------------------------------------------------
    
    if (o3_veg_stress_method=='unset') then 
       allocate(ozone_off_type :: ozone)
    else 
       allocate(ozone_type :: ozone)
    endif

    call ozone%Init(bounds, o3_veg_stress_method)
    
  end function create_and_init_ozone_type

end module OzoneFactoryMod
