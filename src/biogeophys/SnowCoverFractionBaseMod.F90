module SnowCoverFractionBaseMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Abstract base class for methods of computing snow cover fraction
  !
  ! !USES:
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: snow_cover_fraction_base_type

  type, abstract :: snow_cover_fraction_base_type
     private
   contains
     ! Update snow depth and snow fraction
     procedure(UpdateSnowDepthAndFrac_Interface), deferred :: UpdateSnowDepthAndFrac

     ! Add new snow to integrated snow fall
     procedure(AddNewsnowToIntsnow_Interface), deferred :: AddNewsnowToIntsnow

     ! Single-point function: return fractional snow cover during melt
     procedure(FracSnowDuringMelt_Interface), deferred :: FracSnowDuringMelt
  end type snow_cover_fraction_base_type

  abstract interface

     subroutine UpdateSnowDepthAndFrac_Interface(this, bounds, num_c, filter_c, &
          urbpoi, h2osno_total, snowmelt, int_snow, newsnow, bifall, &
          snow_depth, frac_sno)
       ! Update snow depth and snow fraction
       use decompMod, only : bounds_type
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: snow_cover_fraction_base_type

       class(snow_cover_fraction_base_type), intent(in) :: this
       type(bounds_type), intent(in) :: bounds
       integer, intent(in) :: num_c       ! number of columns in filter_c
       integer, intent(in) :: filter_c(:) ! column filter to operate over

       logical  , intent(in)    :: urbpoi( bounds%begc: )       ! true if the given column is urban
       real(r8) , intent(in)    :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
       real(r8) , intent(in)    :: snowmelt( bounds%begc: )     ! total snow melt in the time step (mm H2O)
       real(r8) , intent(in)    :: int_snow( bounds%begc: )     ! integrated snowfall (mm H2O)
       real(r8) , intent(in)    :: newsnow( bounds%begc: )      ! total new snow in the time step (mm H2O)
       real(r8) , intent(in)    :: bifall( bounds%begc: )       ! bulk density of newly fallen dry snow (kg/m3)

       real(r8) , intent(inout) :: snow_depth( bounds%begc: )   ! snow height (m)
       real(r8) , intent(inout) :: frac_sno( bounds%begc: )     ! fraction of ground covered by snow (0 to 1)
     end subroutine UpdateSnowDepthAndFrac_Interface

     subroutine AddNewsnowToIntsnow_Interface(this, bounds, num_c, filter_c, &
          newsnow, h2osno_total, frac_sno, &
          int_snow)
       ! Add new snow to integrated snow fall
       use decompMod, only : bounds_type
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: snow_cover_fraction_base_type

       class(snow_cover_fraction_base_type), intent(in) :: this
       type(bounds_type), intent(in) :: bounds
       integer, intent(in) :: num_c       ! number of columns in filter_c
       integer, intent(in) :: filter_c(:) ! column filter to operate over

       real(r8) , intent(in)    :: newsnow( bounds%begc: )      ! total new snow in the time step (mm H2O)
       real(r8) , intent(in)    :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
       real(r8) , intent(in)    :: frac_sno( bounds%begc: )     ! fraction of ground covered by snow (0 to 1)
       real(r8) , intent(inout) :: int_snow( bounds%begc: )     ! integrated snowfall (mm H2O)
     end subroutine AddNewsnowToIntsnow_Interface

     pure function FracSnowDuringMelt_Interface(this, c, h2osno_total, int_snow) result(frac_sno)
       ! Single-point function: return fractional snow cover during melt
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: snow_cover_fraction_base_type

       real(r8) :: frac_sno  ! function result
       class(snow_cover_fraction_base_type), intent(in) :: this
       integer , intent(in) :: c            ! column we're operating on
       real(r8), intent(in) :: h2osno_total ! total snow water (mm H2O)
       real(r8), intent(in) :: int_snow     ! integrated snowfall (mm H2O)
     end function FracSnowDuringMelt_Interface
  end interface

end module SnowCoverFractionBaseMod
