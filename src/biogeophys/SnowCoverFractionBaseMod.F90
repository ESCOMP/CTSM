module SnowCoverFractionBaseMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Abstract base class for methods of computing snow cover fraction
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use clm_varctl      , only : use_subgrid_fluxes
  use landunit_varcon , only : istdlak

  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: snow_cover_fraction_base_type

  type, abstract :: snow_cover_fraction_base_type
     private
   contains
     ! ------------------------------------------------------------------------
     ! Common subroutines, implemented here
     ! ------------------------------------------------------------------------

     ! Calculate frac_sno_eff given frac_sno
     procedure :: CalcFracSnoEff

     ! ------------------------------------------------------------------------
     ! Subroutines that must be implemented by classes that extend this base class
     ! ------------------------------------------------------------------------

     ! Update snow depth and snow fraction
     procedure(UpdateSnowDepthAndFrac_Interface), deferred :: UpdateSnowDepthAndFrac

     ! Add new snow to integrated snow fall
     procedure(AddNewsnowToIntsnow_Interface), deferred :: AddNewsnowToIntsnow

     ! Single-point function: return fractional snow cover during melt
     procedure(FracSnowDuringMelt_Interface), deferred :: FracSnowDuringMelt
  end type snow_cover_fraction_base_type

  abstract interface

     subroutine UpdateSnowDepthAndFrac_Interface(this, bounds, num_c, filter_c, &
          lun_itype_col, urbpoi, h2osno_total, snowmelt, int_snow, newsnow, bifall, &
          snow_depth, frac_sno, frac_sno_eff)
       ! Update snow depth and snow fraction
       use decompMod, only : bounds_type
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: snow_cover_fraction_base_type

       class(snow_cover_fraction_base_type), intent(in) :: this
       type(bounds_type), intent(in) :: bounds
       integer, intent(in) :: num_c       ! number of columns in filter_c
       integer, intent(in) :: filter_c(:) ! column filter to operate over

       integer  , intent(in)    :: lun_itype_col( bounds%begc: ) ! landunit type for each column
       logical  , intent(in)    :: urbpoi( bounds%begc: )        ! true if the given column is urban
       real(r8) , intent(in)    :: h2osno_total( bounds%begc: )  ! total snow water (mm H2O)
       real(r8) , intent(in)    :: snowmelt( bounds%begc: )      ! total snow melt in the time step (mm H2O)
       real(r8) , intent(in)    :: int_snow( bounds%begc: )      ! integrated snowfall (mm H2O)
       real(r8) , intent(in)    :: newsnow( bounds%begc: )       ! total new snow in the time step (mm H2O)
       real(r8) , intent(in)    :: bifall( bounds%begc: )        ! bulk density of newly fallen dry snow (kg/m3)

       real(r8) , intent(inout) :: snow_depth( bounds%begc: )   ! snow height (m)
       real(r8) , intent(inout) :: frac_sno( bounds%begc: )     ! fraction of ground covered by snow (0 to 1)
       real(r8) , intent(inout) :: frac_sno_eff( bounds%begc: ) ! eff. fraction of ground covered by snow (0 to 1)
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

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
contains

  !-----------------------------------------------------------------------
  subroutine CalcFracSnoEff(this, bounds, num_c, filter_c, &
       lun_itype_col, urbpoi, frac_sno, &
       frac_sno_eff)
    !
    ! !DESCRIPTION:
    ! Calculate frac_sno_eff given frac_sno
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_base_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_c       ! number of columns in filter_c
    integer, intent(in) :: filter_c(:) ! column filter to operate over

    integer  , intent(in)    :: lun_itype_col( bounds%begc: ) ! landunit type for each column
    logical  , intent(in)    :: urbpoi( bounds%begc: )        ! true if the given column is urban
    real(r8) , intent(in)    :: frac_sno( bounds%begc: )      ! fraction of ground covered by snow (0 to 1)
    real(r8) , intent(inout) :: frac_sno_eff( bounds%begc: )  ! eff. fraction of ground covered by snow (0 to 1)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'CalcFracSnoEff'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(lun_itype_col, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(urbpoi, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno_eff, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_c
       c = filter_c(fc)

       if (urbpoi(c) .or. lun_itype_col(c) == istdlak) then
          frac_sno_eff(c) = 1._r8
       else
          if (use_subgrid_fluxes) then
             frac_sno_eff(c) = frac_sno(c)
          else
             frac_sno_eff(c) = 1._r8
          end if
       end if
    end do

  end subroutine CalcFracSnoEff

end module SnowCoverFractionBaseMod
