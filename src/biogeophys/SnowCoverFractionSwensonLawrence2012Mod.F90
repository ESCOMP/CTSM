module SnowCoverFractionSwensonLawrence2012Mod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Implementation of snow_cover_fraction_base_type using the Swenson & Lawrence 2007
  ! parameterization.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_nl_mod     , only : shr_nl_find_group_name
  use shr_mpi_mod    , only : shr_mpi_bcast
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use ncdio_pio      , only : file_desc_t
  use clm_varctl     , only : iulog
  use spmdMod        , only : masterproc, mpicom
  use fileutils      , only : getavu, relavu, opnfil
  use clm_varcon     , only : rpi
  use ColumnType     , only : column_type
  use glcBehaviorMod , only : glc_behavior_type
  use landunit_varcon, only : istice_mec
  use paramUtilMod   , only : readNcdioScalar
  use SnowCoverFractionBaseMod, only : snow_cover_fraction_base_type

  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: snow_cover_fraction_swenson_lawrence_2012_type

  type, extends(snow_cover_fraction_base_type) :: snow_cover_fraction_swenson_lawrence_2012_type
     private
     real(r8) :: int_snow_max = -999._r8  ! limit applied to integrated snowfall when determining changes in snow-covered fraction during melt (mm H2O)
     real(r8) :: accum_factor = -999._r8  ! Accumulation constant for fractional snow covered area (unitless)
     real(r8), allocatable :: n_melt(:)   ! SCA shape parameter
   contains
     procedure, public :: UpdateSnowDepthAndFrac
     procedure, public :: AddNewsnowToIntsnow
     procedure, public :: FracSnowDuringMelt
     procedure, public :: Init

     procedure, private :: ReadNamelist
     procedure, private :: ReadParams
     procedure, private :: CheckValidInputs
     procedure, private :: SetDerivedParameters
  end type snow_cover_fraction_swenson_lawrence_2012_type

  interface snow_cover_fraction_swenson_lawrence_2012_type
     module procedure Constructor
  end interface snow_cover_fraction_swenson_lawrence_2012_type

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
    ! Update snow depth and snow fraction using the SwensonLawrence2012 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(in) :: this
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
    integer  :: fc, c
    real(r8) :: z_avg                                 ! grid cell average snow depth

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

    ! ------------------------------------------------------------------------
    ! Update frac_sno
    ! ------------------------------------------------------------------------

    do fc = 1, num_c
       c = filter_c(fc)

       !======================  FSCA PARAMETERIZATIONS  ======================
       ! fsca parameterization based on *changes* in swe
       if (h2osno_total(c) == 0._r8) then
          if (newsnow(c) > 0._r8) then
             frac_sno(c) = tanh(this%accum_factor * newsnow(c))
          else
             ! NOTE(wjs, 2019-08-07) This resetting of frac_sno to 0 when h2osno_total is 0
             ! may already be done elsewhere; if it isn't, it possibly *should* be done
             ! elsewhere rather than here.
             frac_sno(c) = 0._r8
          end if

       else ! h2osno_total(c) > 0
          if (snowmelt(c) > 0._r8) then
             ! first compute change from melt during previous time step
             frac_sno(c) = this%FracSnowDuringMelt( &
                  c            = c, &
                  h2osno_total = h2osno_total(c), &
                  int_snow     = int_snow(c))
          end if

          if (newsnow(c) > 0._r8) then
             ! Update fsca by new snow event, add to previous fsca

             ! The form in Swenson & Lawrence 2012 (eqn. 3) is:
             ! 1._r8 - (1._r8 - tanh(this%accum_factor * newsnow(c))) * (1._r8 - frac_sno(c))
             !
             ! This form is algebraically equivalent, but simpler and less prone to
             ! roundoff errors (see https://github.com/ESCOMP/ctsm/issues/784)
             frac_sno(c) = frac_sno(c) + tanh(this%accum_factor * newsnow(c)) * (1._r8 - frac_sno(c))

          end if
       end if
    end do

    call this%CalcFracSnoEff(bounds, num_c, filter_c, &
         lun_itype_col = lun_itype_col(begc:endc), &
         urbpoi        = urbpoi(begc:endc), &
         frac_sno      = frac_sno(begc:endc), &
         frac_sno_eff  = frac_sno_eff(begc:endc))

    ! ------------------------------------------------------------------------
    ! Update snow_depth
    ! ------------------------------------------------------------------------

    do fc = 1, num_c
       c = filter_c(fc)

       if (h2osno_total(c) > 0.0_r8) then
          if (frac_sno_eff(c) > 0._r8)then
             snow_depth(c)=snow_depth(c) + newsnow(c)/(bifall(c) * frac_sno_eff(c))
          else
             snow_depth(c)=0._r8
          end if

       else ! h2osno_total == 0
          if (newsnow(c) > 0._r8) then
             z_avg = newsnow(c)/bifall(c)
             snow_depth(c) = z_avg/frac_sno_eff(c)
          else
             snow_depth(c) = 0._r8
          end if
       end if
    end do

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
    ! For this SwensonLawrence2012 parameterization, this involves some care to ensure consistency
    ! between int_snow and other terms.
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_c       ! number of columns in filter_c
    integer, intent(in) :: filter_c(:) ! column filter to operate over

    real(r8) , intent(in)    :: newsnow( bounds%begc: )      ! total new snow in the time step (mm H2O)
    real(r8) , intent(in)    :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
    real(r8) , intent(in)    :: frac_sno( bounds%begc: )     ! fraction of ground covered by snow (0 to 1)
    real(r8) , intent(inout) :: int_snow( bounds%begc: )     ! integrated snowfall (mm H2O)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: temp_intsnow    ! temporary version of int_snow

    character(len=*), parameter :: subname = 'AddNewsnowToIntsnow'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(newsnow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_c
       c = filter_c(fc)

       if (newsnow(c) > 0._r8) then
          ! reset int_snow after accumulation events: make int_snow consistent with new
          ! fsno, h2osno_total
          temp_intsnow= (h2osno_total(c) + newsnow(c)) &
               / (0.5*(cos(rpi*(1._r8-max(frac_sno(c),1e-6_r8))**(1./this%n_melt(c)))+1._r8))
          int_snow(c) = min(1.e8_r8,temp_intsnow)
       end if

       ! NOTE(wjs, 2019-07-25) Sean Swenson and Bill Sacks aren't sure whether this extra
       ! addition of new_snow is correct: it seems to be double-adding newsnow, but we're
       ! not positive that it's wrong. This seems to have been in place ever since the
       ! clm45 branch came to the trunk.
       int_snow(c) = int_snow(c) + newsnow(c)
    end do

  end subroutine AddNewsnowToIntsnow

  !-----------------------------------------------------------------------
  pure function FracSnowDuringMelt(this, c, h2osno_total, int_snow) result(frac_sno)
    !
    ! !DESCRIPTION:
    ! Single-point function giving frac_snow during times when the snow pack is melting
    !
    ! !ARGUMENTS:
    real(r8) :: frac_sno  ! function result
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(in) :: this
    integer , intent(in) :: c            ! column we're operating on
    real(r8), intent(in) :: h2osno_total ! total snow water (mm H2O)
    real(r8), intent(in) :: int_snow     ! integrated snowfall (mm H2O)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: int_snow_limited  ! integrated snowfall, limited to be no greater than int_snow_max [mm]
    real(r8) :: smr

    character(len=*), parameter :: subname = 'FracSnowDuringMelt'
    !-----------------------------------------------------------------------

    int_snow_limited = min(int_snow, this%int_snow_max)
    smr = min(1._r8, h2osno_total/int_snow_limited)

    frac_sno = 1. - (acos(min(1._r8,(2.*smr - 1._r8)))/rpi)**(this%n_melt(c))

  end function FracSnowDuringMelt

  ! ========================================================================
  ! Infrastructure routines (for initialization)
  ! ========================================================================

  function Constructor()
    type(snow_cover_fraction_swenson_lawrence_2012_type) :: Constructor
    ! DO NOTHING (simply return a variable of the appropriate type)
  end function Constructor

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, col, glc_behavior, NLFilename, params_ncid)
    !
    ! !DESCRIPTION:
    ! Initialize this instance of the SwensonLawrence2012 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    type(column_type)       , intent(in)    :: col
    type(glc_behavior_type) , intent(in)    :: glc_behavior
    character(len=*)        , intent(in)    :: NLFilename  ! Namelist filename
    type(file_desc_t)       , intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:
    real(r8)                                :: n_melt_coef    ! n_melt parameter (unitless)
    real(r8) :: n_melt_glcmec  ! SCA shape parameter for glc_mec columns

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    call this%ReadNamelist( &
         NLFilename = NLFilename, &
         n_melt_glcmec = n_melt_glcmec)

    call this%ReadParams( &
         params_ncid = params_ncid, &
         n_melt_coef = n_melt_coef)

    if (masterproc) then
       call this%CheckValidInputs( &
            n_melt_glcmec = n_melt_glcmec)
    end if

    call this%SetDerivedParameters( &
         bounds        = bounds, &
         col           = col, &
         glc_behavior  = glc_behavior, &
         n_melt_coef   = n_melt_coef, &
         n_melt_glcmec = n_melt_glcmec)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename, n_melt_glcmec)
    ! !DESCRIPTION:
    ! Read namelist values needed for SwensonLawrence2012 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(inout) :: this
    character(len=*) , intent(in)  :: NLFilename    ! Namelist filename
    real(r8)         , intent(out) :: n_melt_glcmec ! SCA shape parameter for glc_mec columns
    !
    ! !LOCAL VARIABLES:
    ! local variables corresponding to the namelist variables we read here
    real(r8) :: int_snow_max

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: nmlname = 'scf_swenson_lawrence_2012_inparm'

    character(len=*), parameter :: subname = 'ReadNamelist'
    !-----------------------------------------------------------------------

    namelist /scf_swenson_lawrence_2012_inparm/ int_snow_max, n_melt_glcmec

    int_snow_max = nan
    n_melt_glcmec = nan

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=scf_swenson_lawrence_2012_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast(int_snow_max, mpicom)
    call shr_mpi_bcast(n_melt_glcmec, mpicom)

    this%int_snow_max = int_snow_max

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname // ' settings:'
       write(iulog, nml=scf_swenson_lawrence_2012_inparm)
       write(iulog,*) ' '
    end if

  end subroutine ReadNamelist

  !-----------------------------------------------------------------------
  subroutine ReadParams(this, params_ncid, n_melt_coef)
    !
    ! !DESCRIPTION:
    ! Read netCDF parameters needed for the SwensonLawrence2012 method
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(inout) :: this
    type(file_desc_t) , intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    real(r8)          , intent(out)   :: n_melt_coef ! n_melt parameter (unitless)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'ReadParams'
    !-----------------------------------------------------------------------

    ! Accumulation constant for fractional snow covered area (unitless)
    call readNcdioScalar(params_ncid, 'accum_factor', subname, this%accum_factor)

    ! n_melt parameter (unitless)
    call readNcdioScalar(params_ncid, 'n_melt_coef', subname, n_melt_coef)

  end subroutine ReadParams

  !-----------------------------------------------------------------------
  subroutine CheckValidInputs(this, n_melt_glcmec)
    !
    ! !DESCRIPTION:
    ! Check for validity of parameters read from namelist and netCDF
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(in) :: this
    real(r8), intent(in) :: n_melt_glcmec  ! SCA shape parameter for glc_mec columns
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CheckValidInputs'
    !-----------------------------------------------------------------------

    if (this%int_snow_max <= 0.0_r8) then
       write(iulog,*)'ERROR: int_snow_max = ', this%int_snow_max,' is not supported, must be greater than 0.0.'
       call endrun(msg=' ERROR: invalid value for int_snow_max in CLM namelist. '//&
            errMsg(sourcefile, __LINE__))
    end if

    if (n_melt_glcmec <= 0.0_r8) then
       write(iulog,*)'ERROR: n_melt_glcmec = ', n_melt_glcmec,' is not supported, must be greater than 0.0.'
       call endrun(msg=' ERROR: invalid value for n_melt_glcmec in CLM namelist. '//&
            errMsg(sourcefile, __LINE__))
    end if

  end subroutine CheckValidInputs

  !-----------------------------------------------------------------------
  subroutine SetDerivedParameters(this, bounds, col, glc_behavior, n_melt_coef, n_melt_glcmec)
    !
    ! !DESCRIPTION:
    ! Set parameters that are derived from other inputs
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(inout) :: this
    type(bounds_type)       , intent(in) :: bounds
    type(column_type)       , intent(in) :: col
    type(glc_behavior_type) , intent(in) :: glc_behavior
    real(r8)                , intent(in) :: n_melt_coef
    real(r8)                , intent(in) :: n_melt_glcmec ! SCA shape parameter for glc_mec columns
    !
    ! !LOCAL VARIABLES:
    integer :: c, g

    character(len=*), parameter :: subname = 'SetDerivedParameters'
    !-----------------------------------------------------------------------

    allocate(this%n_melt(bounds%begc:bounds%endc))

    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)

       if (col%lun_itype(c) == istice_mec .and. glc_behavior%allow_multiple_columns_grc(g)) then
          ! ice_mec columns already account for subgrid topographic variability through
          ! their use of multiple elevation classes; thus, to avoid double-accounting for
          ! topographic variability in these columns, we ignore topo_std and use a fixed
          ! value of n_melt.
          this%n_melt(c) = n_melt_glcmec
       else
          this%n_melt(c) = n_melt_coef / max(10._r8, col%topo_std(c))
       end if
    end do

  end subroutine SetDerivedParameters

end module SnowCoverFractionSwensonLawrence2012Mod
