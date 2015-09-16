module SoilWaterRetentionCurveVanGenuchten1980Mod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Implementation of soil_water_retention_curve_type using the Clapp-Hornberg 1978
  ! parameterizations.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: soil_water_retention_curve_vangenuchten_1980_type
  
  type, extends(soil_water_retention_curve_type) :: &
       soil_water_retention_curve_vangenuchten_1980_type
     private
   contains
     procedure :: soil_hk              ! compute hydraulic conductivity
     procedure :: soil_suction         ! compute soil suction potential
     procedure :: soil_suction_inverse ! compute relative saturation at which soil suction is equal to a target value
  end type soil_water_retention_curve_vangenuchten_1980_type

  interface soil_water_retention_curve_vangenuchten_1980_type
     ! initialize a new soil_water_retention_curve_vangenuchten_1980_type object
     module procedure constructor  
  end interface soil_water_retention_curve_vangenuchten_1980_type

contains

  !-----------------------------------------------------------------------
  type(soil_water_retention_curve_vangenuchten_1980_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type soil_water_retention_curve_vangenuchten_1980_type.
    ! For now, this is simply a place-holder.
    !-----------------------------------------------------------------------

  end function constructor

  !-----------------------------------------------------------------------
  subroutine soil_hk(this, c, j, s, imped, soilstate_inst, hk, dhkds)
    !
    ! !DESCRIPTION:
    ! Compute hydraulic conductivity
    !
    ! !USES:
    use SoilStateType  , only : soilstate_type
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_vangenuchten_1980_type), intent(in) :: this
    integer,  intent(in)             :: c        !column index
    integer,  intent(in)             :: j        !level index
    real(r8), intent(in)             :: s        !relative saturation, [0, 1]
    real(r8), intent(in)             :: imped    !ice impedance
    type(soilstate_type), intent(in) :: soilstate_inst
    real(r8), intent(out)            :: hk       !hydraulic conductivity [mm/s]
    real(r8), optional, intent(out)  :: dhkds    !d[hk]/ds   [mm/s]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_hk'
    !-----------------------------------------------------------------------
    
    associate(& 
         hksat             =>    soilstate_inst%hksat_col(c,j)          , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         bsw               =>    soilstate_inst%bsw_col(c,j)              & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
         )


    !compute hydraulic conductivity
    hk=imped*hksat*s**(2._r8*bsw+3._r8)

    !compute the derivative
    if(present(dhkds))then
       dhkds=(2._r8*bsw+3._r8)*hk/s
    endif

    end associate 

  end subroutine soil_hk

  !-----------------------------------------------------------------------
  subroutine soil_suction(this, c, j, s, soilstate_inst, smp, dsmpds)
    !j, 
    ! !DESCRIPTION:
    ! Compute soil suction potential
    !
    ! !USES:
    use SoilStateType  , only : soilstate_type
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_vangenuchten_1980_type), intent(in) :: this
    integer,  intent(in)             :: c       !column index
    integer,  intent(in)             :: j        !level index
    real(r8), intent(in)             :: s        !relative saturation, [0, 1]
    type(soilstate_type), intent(in) :: soilstate_inst
    real(r8), intent(out)            :: smp      !soil suction, negative, [mm]
    real(r8), optional, intent(out)  :: dsmpds   !d[smp]/ds, [mm]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_suction'
    !-----------------------------------------------------------------------
    
    associate(& 
         bsw               =>    soilstate_inst%bsw_col(c,j)            , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                       
         sucsat            =>    soilstate_inst%sucsat_col(c,j)           & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         )

    !compute soil suction potential, negative
    smp = -sucsat*s**(-bsw)

    !compute derivative
    if(present(dsmpds))then
       dsmpds=-bsw*smp/s
    endif

    end associate 

  end subroutine soil_suction

  !-----------------------------------------------------------------------
  subroutine soil_suction_inverse(this, c, j, smp_target, soilstate_inst, s_target)
    !
    ! !DESCRIPTION:
    ! Compute relative saturation at which soil suction is equal to a target value.
    ! This is done by inverting the soil_suction equation to solve for s.
    !
    ! !USES:
    use SoilStateType  , only : soilstate_type
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_vangenuchten_1980_type), intent(in) :: this
    integer,  intent(in)             :: c       !column index
    integer,  intent(in)             :: j        !level index
    type(soilstate_type), intent(in) :: soilstate_inst
    real(r8) , intent(in)  :: smp_target ! target soil suction, negative [mm]
    real(r8) , intent(out) :: s_target   ! relative saturation at which smp = smp_target [0,1]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_suction_inverse'
    !-----------------------------------------------------------------------
    
    associate(& 
         bsw               =>    soilstate_inst%bsw_col(c,j)            , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
         sucsat            =>    soilstate_inst%sucsat_col(c,j)           & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         )

    s_target = (-smp_target/sucsat)**(-1._r8/bsw)

    end associate 

  end subroutine soil_suction_inverse

end module SoilWaterRetentionCurveVanGenuchten1980Mod


