module SnowCoverFractionMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Classes for various methods of computing snow cover fraction
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_nl_mod     , only : shr_nl_find_group_name
  use shr_mpi_mod    , only : shr_mpi_bcast
  use fileutils      , only : getavu, relavu, opnfil
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use spmdMod        , only : masterproc, mpicom
  use clm_varcon     , only : rpi
  use clm_varctl     , only : iulog, use_subgrid_fluxes
  use ncdio_pio      , only : file_desc_t
  use paramUtilMod   , only : readNcdioScalar
  use ColumnType     , only : column_type
  use glcBehaviorMod , only : glc_behavior_type
  use landunit_varcon, only : istice_mec
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: CreateAndInitScfMethod
  public :: snow_cover_fraction_type

  type, abstract :: snow_cover_fraction_type
     private
   contains
     ! Initialize this instance
     procedure(Init_Interface), deferred :: Init

     ! Update snow depth and snow fraction
     procedure(UpdateSnowDepthAndFrac_Interface), deferred :: UpdateSnowDepthAndFrac

     ! Add new snow to integrated snow fall
     procedure(AddNewsnowToIntsnow_Interface), deferred :: AddNewsnowToIntsnow

     ! Single-point function: return fractional snow cover during melt
     procedure(FracSnowDuringMelt_Interface), deferred :: FracSnowDuringMelt
  end type snow_cover_fraction_type

  abstract interface

     subroutine Init_Interface(this, bounds, col, glc_behavior, NLFilename, params_ncid)
       ! Initialize this instance
       use decompMod, only : bounds_type
       use ColumnType, only : column_type
       use glcBehaviorMod, only : glc_behavior_type
       use ncdio_pio, only : file_desc_t
       import :: snow_cover_fraction_type

       class(snow_cover_fraction_type), intent(inout) :: this
       type(bounds_type), intent(in) :: bounds
       type(column_type), intent(in) :: col
       type(glc_behavior_type), intent(in) :: glc_behavior
       character(len=*), intent(in) :: NLFilename ! Namelist filename
       type(file_desc_t), intent(inout) :: params_ncid ! pio netCDF file id for parameter file
     end subroutine Init_Interface

     subroutine UpdateSnowDepthAndFrac_Interface(this, bounds, num_c, filter_c, &
          urbpoi, h2osno_total, snowmelt, int_snow, newsnow, bifall, &
          snow_depth, frac_sno)
       ! Update snow depth and snow fraction
       use decompMod, only : bounds_type
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: snow_cover_fraction_type

       class(snow_cover_fraction_type), intent(in) :: this
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
       import :: snow_cover_fraction_type

       class(snow_cover_fraction_type), intent(in) :: this
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
       import :: snow_cover_fraction_type

       real(r8) :: frac_sno  ! function result
       class(snow_cover_fraction_type), intent(in) :: this
       integer , intent(in) :: c            ! column we're operating on
       real(r8), intent(in) :: h2osno_total ! total snow water (mm H2O)
       real(r8), intent(in) :: int_snow     ! integrated snowfall (mm H2O)
     end function FracSnowDuringMelt_Interface
  end interface

  type, extends(snow_cover_fraction_type) :: snow_cover_fraction_niu_yang_2007_type
     ! This implements the N&Y07 (Niu & Yang 2007) parameterization of snow cover fraction
     private
     real(r8) :: zlnd = -999._r8  ! Roughness length for soil (m)
   contains
     procedure, public :: Init => InitNiuYang2007
     procedure, public :: UpdateSnowDepthAndFrac => UpdateSnowDepthAndFracNiuYang2007
     procedure, public :: AddNewsnowToIntsnow => AddNewsnowToIntsnowNiuYang2007
     procedure, public :: FracSnowDuringMelt => FracSnowDuringMeltNiuYang2007

     procedure, private :: ReadParamsNiuYang2007
  end type snow_cover_fraction_niu_yang_2007_type

  interface snow_cover_fraction_niu_yang_2007_type
     module procedure ConstructorNiuYang2007
  end interface snow_cover_fraction_niu_yang_2007_type

  type, extends(snow_cover_fraction_type) :: snow_cover_fraction_swenson_lawrence_2012_type
     ! This implements the S&L12 (Swenson & Lawrence 2012) parameterization of snow cover fraction
     private
     real(r8) :: int_snow_max = -999._r8  ! limit applied to integrated snowfall when determining changes in snow-covered fraction during melt (mm H2O)
     real(r8) :: accum_factor = -999._r8  ! Accumulation constant for fractional snow covered area (unitless)
     real(r8), allocatable :: n_melt(:)   ! SCA shape parameter
   contains
     procedure, public :: Init => InitSwensonLawrence2012
     procedure, public :: UpdateSnowDepthAndFrac => UpdateSnowDepthAndFracSwensonLawrence2012
     procedure, public :: AddNewsnowToIntsnow => AddNewsnowToIntsnowSwensonLawrence2012
     procedure, public :: FracSnowDuringMelt => FracSnowDuringMeltSwensonLawrence2012

     procedure, private :: ReadNamelistSwensonLawrence2012
     procedure, private :: ReadParamsSwensonLawrence2012
     procedure, private :: CheckValidInputsSwensonLawrence2012
     procedure, private :: SetDerivedParametersSwensonLawrence2012
  end type snow_cover_fraction_swenson_lawrence_2012_type

  interface snow_cover_fraction_swenson_lawrence_2012_type
     module procedure ConstructorSwensonLawrence2012
  end interface snow_cover_fraction_swenson_lawrence_2012_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
contains

  ! ========================================================================
  ! Factory method
  ! ========================================================================

  !-----------------------------------------------------------------------
  function CreateAndInitScfMethod(snow_cover_fraction_method, &
       bounds, col, glc_behavior, NLFilename, params_ncid) &
       result(scf_method)
    !
    ! !DESCRIPTION:
    ! Create an instance of the appropriate snow_cover_fraction_type (based on
    ! snow_cover_fraction_method), initialize it and return it
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_type), allocatable :: scf_method  ! function result

    character(len=*)        , intent(in)    :: snow_cover_fraction_method
    type(bounds_type)       , intent(in)    :: bounds
    type(column_type)       , intent(in)    :: col
    type(glc_behavior_type) , intent(in)    :: glc_behavior
    character(len=*)        , intent(in)    :: NLFilename ! Namelist filename
    type(file_desc_t)       , intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CreateAndInitScfMethod'
    !-----------------------------------------------------------------------

    select case (snow_cover_fraction_method)
    case ('SwensonLawrence2012')
       allocate(scf_method, &
            source = snow_cover_fraction_swenson_lawrence_2012_type())
    case ('NiuYang2007')
       allocate(scf_method, &
            source = snow_cover_fraction_niu_yang_2007_type())
    case default
       write(iulog,*) subname//' ERROR: unknown snow_cover_fraction_method: ', &
            snow_cover_fraction_method
       call endrun(msg = 'unknown snow_cover_fraction_method', &
            additional_msg = errMsg(sourcefile, __LINE__))
    end select

    call scf_method%Init( &
         bounds       = bounds, &
         col          = col, &
         glc_behavior = glc_behavior, &
         NLFilename   = NLFilename, &
         params_ncid  = params_ncid)

  end function CreateAndInitScfMethod

  ! ========================================================================
  ! Methods for N&Y07 (Niu & Yang 2007) parameterization
  ! ========================================================================

  function ConstructorNiuYang2007()
    type(snow_cover_fraction_niu_yang_2007_type) :: ConstructorNiuYang2007
    ! DO NOTHING (simply return a variable of the appropriate type)
  end function ConstructorNiuYang2007

  !-----------------------------------------------------------------------
  subroutine InitNiuYang2007(this, bounds, col, glc_behavior, NLFilename, params_ncid)
    !
    ! !DESCRIPTION:
    ! Initialize this instance of the NiuYang2007 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_niu_yang_2007_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    type(column_type)       , intent(in)    :: col
    type(glc_behavior_type) , intent(in)    :: glc_behavior
    character(len=*)        , intent(in)    :: NLFilename  ! Namelist filename
    type(file_desc_t)       , intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter             :: subname = 'InitNiuYang2007'
    !-----------------------------------------------------------------------

    if (use_subgrid_fluxes) then
       write(iulog,*) 'ERROR: Attempt to use N&Y07 snow cover fraction parameterization with use_subgrid_fluxes.'
       write(iulog,*) 'These two options are incompatible.'
       call endrun('N&Y07 snow cover fraction parameterization incompatible with use_subgrid_fluxes')
    end if

    call this%ReadParamsNiuYang2007( &
         params_ncid = params_ncid)

  end subroutine InitNiuYang2007

  !-----------------------------------------------------------------------
  subroutine ReadParamsNiuYang2007(this, params_ncid)
    !
    ! !DESCRIPTION:
    ! Read netCDF parameters needed for the NiuYang2007 method
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_niu_yang_2007_type) , intent(inout) :: this
    type(file_desc_t) , intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'ReadParamsNiuYang2007'
    !-----------------------------------------------------------------------

    ! Roughness length for soil (m)
    call readNcdioScalar(params_ncid, 'zlnd', subname, this%zlnd)

  end subroutine ReadParamsNiuYang2007

  !-----------------------------------------------------------------------
  subroutine UpdateSnowDepthAndFracNiuYang2007(this, bounds, num_c, filter_c, &
       urbpoi, h2osno_total, snowmelt, int_snow, newsnow, bifall, &
       snow_depth, frac_sno)
    !
    ! !DESCRIPTION:
    ! Update snow depth and snow fraction using the NiuYang2007 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_niu_yang_2007_type), intent(in) :: this
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
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c

    character(len=*), parameter :: subname = 'UpdateSnowDepthAndFracNiuYang2007'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(urbpoi, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snowmelt, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(newsnow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bifall, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snow_depth, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)

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

  end subroutine UpdateSnowDepthAndFracNiuYang2007

  !-----------------------------------------------------------------------
  subroutine AddNewsnowToIntsnowNiuYang2007(this, bounds, num_c, filter_c, &
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

    character(len=*), parameter :: subname = 'AddNewsnowToIntsnowNiuYang2007'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(newsnow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_c
       c = filter_c(fc)
       int_snow(c) = int_snow(c) + newsnow(c)
    end do

  end subroutine AddNewsnowToIntsnowNiuYang2007

  !-----------------------------------------------------------------------
  pure function FracSnowDuringMeltNiuYang2007(this, c, h2osno_total, int_snow) result(frac_sno)
    !
    ! !DESCRIPTION:
    ! Single-point function giving frac_snow during times when the snow pack is melting
    !
    ! This is currently not implemented for the N&Y07 method, so simply returns NaN. (We
    ! can't abort, since this is a pure function.) That's okay since this is only called
    ! if use_subgrid_fluxes is true, and use_subgrid_fluxes cannot currently be used with
    ! the N&Y07 method.
    !
    ! !ARGUMENTS:
    real(r8) :: frac_sno  ! function result
    class(snow_cover_fraction_niu_yang_2007_type), intent(in) :: this
    integer , intent(in) :: c            ! column we're operating on
    real(r8), intent(in) :: h2osno_total ! total snow water (mm H2O)
    real(r8), intent(in) :: int_snow     ! integrated snowfall (mm H2O)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'FracSnowDuringMeltNiuYang2007'
    !-----------------------------------------------------------------------

    frac_sno = nan

  end function FracSnowDuringMeltNiuYang2007

  ! ========================================================================
  ! Methods for SwensonLawrence2012 parameterization
  ! ========================================================================

  function ConstructorSwensonLawrence2012()
    type(snow_cover_fraction_swenson_lawrence_2012_type) :: ConstructorSwensonLawrence2012
    ! DO NOTHING (simply return a variable of the appropriate type)
  end function ConstructorSwensonLawrence2012

  !-----------------------------------------------------------------------
  subroutine InitSwensonLawrence2012(this, bounds, col, glc_behavior, NLFilename, params_ncid)
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

    character(len=*), parameter :: subname = 'InitSwensonLawrence2012'
    !-----------------------------------------------------------------------

    call this%ReadNamelistSwensonLawrence2012( &
         NLFilename = NLFilename, &
         n_melt_glcmec = n_melt_glcmec)

    call this%ReadParamsSwensonLawrence2012( &
         params_ncid = params_ncid, &
         n_melt_coef = n_melt_coef)

    if (masterproc) then
       call this%CheckValidInputsSwensonLawrence2012( &
            n_melt_glcmec = n_melt_glcmec)
    end if

    call this%SetDerivedParametersSwensonLawrence2012( &
         bounds        = bounds, &
         col           = col, &
         glc_behavior  = glc_behavior, &
         n_melt_coef   = n_melt_coef, &
         n_melt_glcmec = n_melt_glcmec)

  end subroutine InitSwensonLawrence2012

  !-----------------------------------------------------------------------
  subroutine ReadNamelistSwensonLawrence2012(this, NLFilename, n_melt_glcmec)
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

    character(len=*), parameter :: subname = 'ReadNamelistSwensonLawrence2012'
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

  end subroutine ReadNamelistSwensonLawrence2012

  !-----------------------------------------------------------------------
  subroutine ReadParamsSwensonLawrence2012(this, params_ncid, n_melt_coef)
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

    character(len=*), parameter :: subname = 'ReadParamsSwensonLawrence2012'
    !-----------------------------------------------------------------------

    ! Accumulation constant for fractional snow covered area (unitless)
    call readNcdioScalar(params_ncid, 'accum_factor', subname, this%accum_factor)

    ! n_melt parameter (unitless)
    call readNcdioScalar(params_ncid, 'n_melt_coef', subname, n_melt_coef)

  end subroutine ReadParamsSwensonLawrence2012

  !-----------------------------------------------------------------------
  subroutine CheckValidInputsSwensonLawrence2012(this, n_melt_glcmec)
    !
    ! !DESCRIPTION:
    ! Check for validity of parameters read from namelist and netCDF
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(in) :: this
    real(r8), intent(in) :: n_melt_glcmec  ! SCA shape parameter for glc_mec columns
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CheckValidInputsSwensonLawrence2012'
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

  end subroutine CheckValidInputsSwensonLawrence2012

  !-----------------------------------------------------------------------
  subroutine SetDerivedParametersSwensonLawrence2012(this, bounds, col, glc_behavior, n_melt_coef, n_melt_glcmec)
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

    character(len=*), parameter :: subname = 'SetDerivedParametersSwensonLawrence2012'
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

  end subroutine SetDerivedParametersSwensonLawrence2012

  !-----------------------------------------------------------------------
  subroutine UpdateSnowDepthAndFracSwensonLawrence2012(this, bounds, num_c, filter_c, &
       urbpoi, h2osno_total, snowmelt, int_snow, newsnow, bifall, &
       snow_depth, frac_sno)
    !
    ! !DESCRIPTION:
    ! Update snow depth and snow fraction using the SwensonLawrence2012 parameterization
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_swenson_lawrence_2012_type), intent(in) :: this
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
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: z_avg                                 ! grid cell average snow depth

    character(len=*), parameter :: subname = 'UpdateSnowDepthAndFracSwensonLawrence2012'
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL((ubound(urbpoi, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snowmelt, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(newsnow, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(bifall, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(snow_depth, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(frac_sno, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_c
       c = filter_c(fc)
       if (h2osno_total(c) > 0.0_r8) then
          !======================  FSCA PARAMETERIZATIONS  ======================
          ! fsca parameterization based on *changes* in swe
          ! first compute change from melt during previous time step
          if(snowmelt(c) > 0._r8) then
             frac_sno(c) = this%FracSnowDuringMelt( &
                  c            = c, &
                  h2osno_total = h2osno_total(c), &
                  int_snow     = int_snow(c))
          end if

          ! update fsca by new snow event, add to previous fsca
          if (newsnow(c) > 0._r8) then
             frac_sno(c) = 1._r8 - (1._r8 - tanh(this%accum_factor * newsnow(c))) * (1._r8 - frac_sno(c))
          end if

          !====================================================================

          ! for subgrid fluxes
          if (use_subgrid_fluxes .and. .not. urbpoi(c)) then
             if (frac_sno(c) > 0._r8)then
                snow_depth(c)=snow_depth(c) + newsnow(c)/(bifall(c) * frac_sno(c))
             else
                snow_depth(c)=0._r8
             end if
          else
             ! for uniform snow cover
             snow_depth(c)=snow_depth(c)+newsnow(c)/bifall(c)
          end if

       else ! h2osno_total == 0
          ! initialize frac_sno and snow_depth when no snow present initially
          if (newsnow(c) > 0._r8) then
             z_avg = newsnow(c)/bifall(c)
             frac_sno(c) = tanh(this%accum_factor * newsnow(c))

             ! update snow_depth to be consistent with frac_sno, z_avg
             if (use_subgrid_fluxes .and. .not. urbpoi(c)) then
                snow_depth(c)=z_avg/frac_sno(c)
             else
                snow_depth(c)=newsnow(c)/bifall(c)
             end if
          else
             snow_depth(c) = 0._r8
             frac_sno(c) = 0._r8
          end if
       end if
    end do

  end subroutine UpdateSnowDepthAndFracSwensonLawrence2012

  !-----------------------------------------------------------------------
  subroutine AddNewsnowToIntsnowSwensonLawrence2012(this, bounds, num_c, filter_c, &
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

    character(len=*), parameter :: subname = 'AddNewsnowToIntsnowSwensonLawrence2012'
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

  end subroutine AddNewsnowToIntsnowSwensonLawrence2012

  !-----------------------------------------------------------------------
  pure function FracSnowDuringMeltSwensonLawrence2012(this, c, h2osno_total, int_snow) result(frac_sno)
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

    character(len=*), parameter :: subname = 'FracSnowDuringMeltSwensonLawrence2012'
    !-----------------------------------------------------------------------

    int_snow_limited = min(int_snow, this%int_snow_max)
    smr = min(1._r8, h2osno_total/int_snow_limited)

    frac_sno = 1. - (acos(min(1._r8,(2.*smr - 1._r8)))/rpi)**(this%n_melt(c))

  end function FracSnowDuringMeltSwensonLawrence2012

end module SnowCoverFractionMod
