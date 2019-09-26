module SnowCoverFractionNiuYang2007Mod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Implementation of snow_cover_fraction_base_type using the Niu & Yang 2007
  ! parameterization.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use ncdio_pio      , only : file_desc_t
  use clm_varctl     , only : iulog, use_subgrid_fluxes
  use paramUtilMod   , only : readNcdioScalar
  use SnowCoverFractionBaseMod, only : snow_cover_fraction_base_type

  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: snow_cover_fraction_niu_yang_2007_type

  type, extends(snow_cover_fraction_base_type) :: snow_cover_fraction_niu_yang_2007_type
     private
     real(r8) :: zlnd = -999._r8  ! Roughness length for soil (m)
   contains
     procedure, public :: UpdateSnowDepthAndFrac
     procedure, public :: AddNewsnowToIntsnow
     procedure, public :: FracSnowDuringMelt
     procedure, public :: Init

     procedure, private :: ReadParams
  end type snow_cover_fraction_niu_yang_2007_type

  interface snow_cover_fraction_niu_yang_2007_type
     module procedure Constructor
  end interface snow_cover_fraction_niu_yang_2007_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
contains

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine UpdateSnowDepthAndFrac(this, bounds, num_c, filter_c, &
       lun_itype_col, urbpoi, h2osno_total, snowmelt, int_snow, newsnow, bifall, &
       snow_depth, frac_sno, frac_sno_eff)
    !
    ! !DESCRIPTION:
    ! Update snow depth and snow fraction using the NiuYang2007 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_niu_yang_2007_type), intent(in) :: this
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

    real(r8) , intent(inout) :: snow_depth( bounds%begc: )    ! snow height (m)
    real(r8) , intent(inout) :: frac_sno( bounds%begc: )      ! fraction of ground covered by snow (0 to 1)
    real(r8) , intent(inout) :: frac_sno_eff( bounds%begc: )  ! eff. fraction of ground covered by snow (0 to 1)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'UpdateSnowDepthAndFrac'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(lun_itype_col, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(urbpoi, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snowmelt, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(newsnow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bifall, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snow_depth, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno_eff, 1) == bounds%endc), sourcefile, __LINE__)

    associate( &
         begc => bounds%begc, &
         endc => bounds%endc  &
         )

    do fc = 1, num_c
       c = filter_c(fc)

       if (h2osno_total(c) == 0.0_r8) then
          ! NOTE(wjs, 2019-08-03) This resetting of snow_depth to 0 when h2osno_total is
          ! 0 may already be done elsewhere; if it isn't, it probably *should* be done
          ! elsewhere rather than here.
          snow_depth(c) = 0._r8
       end if

       snow_depth(c) = snow_depth(c) + newsnow(c) / bifall(c)

       if (snow_depth(c) > 0.0_r8) then
          frac_sno(c) = tanh(snow_depth(c) / (2.5_r8 * this%zlnd * &
               (min(800._r8,(h2osno_total(c)+ newsnow(c))/snow_depth(c))/100._r8)**1._r8) )
       else
          frac_sno(c) = 0._r8
       end if

       ! NOTE(wjs, 2019-08-03) I'm not sure what the intent is of the following block, and
       ! it feels to me like this should be looking at (h2osno_total+newsnow), both in the
       ! condition and in the min, but for now I'm keeping the pre-existing logic.
       if (h2osno_total(c) > 0.0_r8 .and. h2osno_total(c) < 1.0_r8) then
          frac_sno(c) = min(frac_sno(c), h2osno_total(c))
       end if
    end do

    call this%CalcFracSnoEff(bounds, num_c, filter_c, &
         lun_itype_col = lun_itype_col(begc:endc), &
         urbpoi        = urbpoi(begc:endc), &
         frac_sno      = frac_sno(begc:endc), &
         frac_sno_eff  = frac_sno_eff(begc:endc))

    end associate
  end subroutine UpdateSnowDepthAndFrac

  !-----------------------------------------------------------------------
  subroutine AddNewsnowToIntsnow(this, bounds, num_c, filter_c, &
       newsnow, h2osno_total, frac_sno, &
       int_snow)
    !
    ! !DESCRIPTION:
    ! Add new snow to integrated snow fall
    !
    ! This is straightforward for this NiuYang2007 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_niu_yang_2007_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_c       ! number of columns in filter_c
    integer, intent(in) :: filter_c(:) ! column filter to operate over

    real(r8) , intent(in)    :: newsnow( bounds%begc: )      ! total new snow in the time step (mm H2O)
    real(r8) , intent(in)    :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
    real(r8) , intent(in)    :: frac_sno( bounds%begc: )     ! fraction of ground covered by snow (0 to 1)
    real(r8) , intent(inout) :: int_snow( bounds%begc: )     ! integrated snowfall (mm H2O)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'AddNewsnowToIntsnow'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(newsnow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_c
       c = filter_c(fc)
       int_snow(c) = int_snow(c) + newsnow(c)
    end do

  end subroutine AddNewsnowToIntsnow

  !-----------------------------------------------------------------------
  pure function FracSnowDuringMelt(this, c, h2osno_total, int_snow) result(frac_sno)
    !
    ! !DESCRIPTION:
    ! Single-point function giving frac_snow during times when the snow pack is melting
    !
    ! This is currently not implemented for the NiuYang07 method, so simply returns NaN. (We
    ! can't abort, since this is a pure function.) That's okay since this is only called
    ! if use_subgrid_fluxes is true, and use_subgrid_fluxes cannot currently be used with
    ! the NiuYang07 method.
    !
    ! !ARGUMENTS:
    real(r8) :: frac_sno  ! function result
    class(snow_cover_fraction_niu_yang_2007_type), intent(in) :: this
    integer , intent(in) :: c            ! column we're operating on
    real(r8), intent(in) :: h2osno_total ! total snow water (mm H2O)
    real(r8), intent(in) :: int_snow     ! integrated snowfall (mm H2O)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'FracSnowDuringMelt'
    !-----------------------------------------------------------------------

    frac_sno = nan

  end function FracSnowDuringMelt

  ! ========================================================================
  ! Infrastructure routines (for initialization)
  ! ========================================================================

  function Constructor()
    type(snow_cover_fraction_niu_yang_2007_type) :: Constructor
    ! DO NOTHING (simply return a variable of the appropriate type)
  end function Constructor

  !-----------------------------------------------------------------------
  subroutine Init(this, params_ncid)
    !
    ! !DESCRIPTION:
    ! Initialize this instance of the NiuYang2007 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_niu_yang_2007_type) , intent(inout) :: this
    type(file_desc_t)       , intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter             :: subname = 'Init'
    !-----------------------------------------------------------------------

    if (use_subgrid_fluxes) then
       write(iulog,*) 'ERROR: Attempt to use NiuYang07 snow cover fraction parameterization with use_subgrid_fluxes.'
       write(iulog,*) 'These two options are incompatible.'
       call endrun('NiuYang07 snow cover fraction parameterization incompatible with use_subgrid_fluxes')
    end if

    call this%ReadParams( &
         params_ncid = params_ncid)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine ReadParams(this, params_ncid)
    !
    ! !DESCRIPTION:
    ! Read netCDF parameters needed for the NiuYang2007 method
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_niu_yang_2007_type) , intent(inout) :: this
    type(file_desc_t) , intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'ReadParams'
    !-----------------------------------------------------------------------

    ! Roughness length for soil (m)
    call readNcdioScalar(params_ncid, 'zlnd', subname, this%zlnd)

  end subroutine ReadParams

end module SnowCoverFractionNiuYang2007Mod
