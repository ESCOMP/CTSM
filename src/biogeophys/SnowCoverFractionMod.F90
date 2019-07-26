module SnowCoverFractionMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Classes for various methods of computing snow cover fraction
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use spmdMod        , only : masterproc, mpicom
  use clm_varcon     , only : rpi
  use clm_varctl     , only : iulog, use_subgrid_fluxes
  use ncdio_pio      , only : file_desc_t
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
       use decompMod, only : bounds_type
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: snow_cover_fraction_type

       class(snow_cover_fraction_type), intent(in) :: this
       type(bounds_type), intent(in) :: bounds
       integer, intent(in) :: num_c       ! number of columns in filter_c
       integer, intent(in) :: filter_c(:) ! column filter to operate over

       logical, intent(in) :: urbpoi( bounds%begc: ) ! true if the given column is urban
       real(r8), intent(in) :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
       ! FIXME(wjs, 2019-07-22) document the following arguments
       real(r8), intent(in) :: snowmelt( bounds%begc: )
       real(r8), intent(in) :: int_snow( bounds%begc: )
       real(r8), intent(in) :: newsnow( bounds%begc: )
       real(r8), intent(in) :: bifall( bounds%begc: )

       real(r8), intent(inout) :: snow_depth( bounds%begc: )
       real(r8), intent(inout) :: frac_sno( bounds%begc: )
     end subroutine UpdateSnowDepthAndFrac_Interface

     subroutine AddNewsnowToIntsnow_Interface(this, bounds, num_c, filter_c, &
          newsnow, h2osno_total, frac_sno, &
          int_snow)
       use decompMod, only : bounds_type
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: snow_cover_fraction_type

       class(snow_cover_fraction_type), intent(in) :: this
       type(bounds_type), intent(in) :: bounds
       integer, intent(in) :: num_c       ! number of columns in filter_c
       integer, intent(in) :: filter_c(:) ! column filter to operate over

       ! FIXME(wjs, 2019-07-22) document the following arguments
       real(r8), intent(in) :: newsnow( bounds%begc: )
       real(r8), intent(in) :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
       real(r8), intent(in) :: frac_sno( bounds%begc: )
       real(r8), intent(inout) :: int_snow( bounds%begc: )
     end subroutine AddNewsnowToIntsnow_Interface

     pure function FracSnowDuringMelt_Interface(this, c, h2osno_total, int_snow) result(frac_sno)
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: snow_cover_fraction_type

       real(r8) :: frac_sno  ! function result
       class(snow_cover_fraction_type), intent(in) :: this
       integer , intent(in) :: c            ! column we're operating on
       real(r8), intent(in) :: h2osno_total ! total snow water (mm H2O)
       real(r8), intent(in) :: int_snow     ! integrated snowfall (mm H2O)
     end function FracSnowDuringMelt_Interface
  end interface

  type, extends(snow_cover_fraction_type) :: snow_cover_fraction_clm5_type
     private
     real(r8) :: int_snow_max = -999._r8  ! limit applied to integrated snowfall when determining changes in snow-covered fraction during melt (mm H2O)
     real(r8) :: accum_factor = -999._r8  ! Accumulation constant for fractional snow covered area (unitless)
     real(r8), allocatable :: n_melt(:)   ! SCA shape parameter
   contains
     procedure, public :: Init => InitClm5
     procedure, public :: UpdateSnowDepthAndFrac => UpdateSnowDepthAndFracClm5
     procedure, public :: AddNewsnowToIntsnow => AddNewsnowToIntsnowClm5
     procedure, public :: FracSnowDuringMelt => FracSnowDuringMeltClm5

     procedure, private :: ReadNamelistClm5
     procedure, private :: ReadParamsClm5
     procedure, private :: CheckValidInputsClm5
     procedure, private :: SetDerivedParametersClm5
  end type snow_cover_fraction_clm5_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
contains

  ! ========================================================================
  ! Factory method
  ! ========================================================================

  !-----------------------------------------------------------------------
  function CreateAndInitScfMethod(bounds, col, glc_behavior, NLFilename, params_ncid) result(scf_method)
    !
    ! !DESCRIPTION:
    ! Create an instance of the appropriate snow_cover_fraction_type, initialize it and
    ! return it
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_type), allocatable :: scf_method  ! function result
    type(bounds_type), intent(in) :: bounds
    type(column_type), intent(in) :: col
    type(glc_behavior_type), intent(in) :: glc_behavior
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    type(file_desc_t), intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CreateAndInitScfMethod'
    !-----------------------------------------------------------------------

    ! FIXME(wjs, 2019-07-26) Have some logic to create the appropriate snow cover
    ! fraction method (based on something read from controlMod/clm_varctl???)
    allocate(scf_method, &
         source = snow_cover_fraction_clm5_type())

    call scf_method%Init( &
         bounds       = bounds, &
         col          = col, &
         glc_behavior = glc_behavior, &
         NLFilename   = NLFilename, &
         params_ncid  = params_ncid)

  end function CreateAndInitScfMethod

  ! ========================================================================
  ! Methods for CLM5 parameterization
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine InitClm5(this, bounds, col, glc_behavior, NLFilename, params_ncid)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(column_type), intent(in) :: col
    type(glc_behavior_type), intent(in) :: glc_behavior
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    type(file_desc_t), intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:
    real(r8) :: n_melt_coef    ! n_melt parameter (unitless)
    real(r8) :: n_melt_glcmec  ! SCA shape parameter for glc_mec columns

    character(len=*), parameter :: subname = 'InitClm5'
    !-----------------------------------------------------------------------

    call this%ReadNamelistClm5( &
         NLFilename = NLFilename, &
         n_melt_glcmec = n_melt_glcmec)

    call this%ReadParamsClm5( &
         params_ncid = params_ncid, &
         n_melt_coef = n_melt_coef)

    if (masterproc) then
       call this%CheckValidInputsClm5( &
            n_melt_glcmec = n_melt_glcmec)
    end if

    call this%SetDerivedParametersClm5( &
         bounds        = bounds, &
         col           = col, &
         glc_behavior  = glc_behavior, &
         n_melt_coef   = n_melt_coef, &
         n_melt_glcmec = n_melt_glcmec)

  end subroutine InitClm5

  !-----------------------------------------------------------------------
  subroutine ReadNamelistClm5(this, NLFilename, n_melt_glcmec)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use shr_mpi_mod    , only : shr_mpi_bcast
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    real(r8), intent(out) :: n_melt_glcmec  ! SCA shape parameter for glc_mec columns
    !
    ! !LOCAL VARIABLES:
    ! local variables corresponding to the namelist variables we read here
    real(r8) :: int_snow_max

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: nmlname = 'scf_clm5_inparm'

    character(len=*), parameter :: subname = 'ReadNamelistClm5'
    !-----------------------------------------------------------------------

    namelist /scf_clm5_inparm/ int_snow_max, n_melt_glcmec

    int_snow_max = nan
    n_melt_glcmec = nan

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=scf_clm5_inparm, iostat=ierr)
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
       write(iulog, nml=scf_clm5_inparm)
       write(iulog,*) ' '
    end if

  end subroutine ReadNamelistClm5

  !-----------------------------------------------------------------------
  subroutine ReadParamsClm5(this, params_ncid, n_melt_coef)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(inout) :: this
    type(file_desc_t), intent(inout) :: params_ncid  ! pio netCDF file id for parameter file
    real(r8), intent(out) :: n_melt_coef ! n_melt parameter (unitless)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'ReadParamsClm5'
    !-----------------------------------------------------------------------

    ! Accumulation constant for fractional snow covered area (unitless)
    call readNcdioScalar(params_ncid, 'accum_factor', subname, this%accum_factor)

    ! n_melt parameter (unitless)
    call readNcdioScalar(params_ncid, 'n_melt_coef', subname, n_melt_coef)

  end subroutine ReadParamsClm5

  !-----------------------------------------------------------------------
  subroutine CheckValidInputsClm5(this, n_melt_glcmec)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(in) :: this
    real(r8), intent(out) :: n_melt_glcmec  ! SCA shape parameter for glc_mec columns
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CheckValidInputsClm5'
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

  end subroutine CheckValidInputsClm5

  !-----------------------------------------------------------------------
  subroutine SetDerivedParametersClm5(this, bounds, col, glc_behavior, n_melt_coef, n_melt_glcmec)
    !
    ! !DESCRIPTION:
    ! Set parameters that are derived from other inputs
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    type(column_type), intent(in) :: col
    type(glc_behavior_type), intent(in) :: glc_behavior
    real(r8), intent(in) :: n_melt_coef
    real(r8), intent(in) :: n_melt_glcmec ! SCA shape parameter for glc_mec columns
    !
    ! !LOCAL VARIABLES:
    integer :: c, g

    character(len=*), parameter :: subname = 'SetDerivedParametersClm5'
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

  end subroutine SetDerivedParametersClm5

  !-----------------------------------------------------------------------
  subroutine UpdateSnowDepthAndFracClm5(this, bounds, num_c, filter_c, &
       urbpoi, h2osno_total, snowmelt, int_snow, newsnow, bifall, &
       snow_depth, frac_sno)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_c       ! number of columns in filter_c
    integer, intent(in) :: filter_c(:) ! column filter to operate over

    logical, intent(in) :: urbpoi( bounds%begc: ) ! true if the given column is urban
    real(r8), intent(in) :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
    ! FIXME(wjs, 2019-07-22) document the following arguments
    real(r8), intent(in) :: snowmelt( bounds%begc: )
    real(r8), intent(in) :: int_snow( bounds%begc: )
    real(r8), intent(in) :: newsnow( bounds%begc: )
    real(r8), intent(in) :: bifall( bounds%begc: )

    real(r8), intent(inout) :: snow_depth( bounds%begc: )
    real(r8), intent(inout) :: frac_sno( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: fsno_new
    real(r8) :: fmelt
    real(r8) :: z_avg                                 ! grid cell average snow depth

    character(len=*), parameter :: subname = 'UpdateSnowDepthAndFracClm5'
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
             fsno_new = 1._r8 - (1._r8 - tanh(this%accum_factor * newsnow(c))) * (1._r8 - frac_sno(c))
             frac_sno(c) = fsno_new
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
             fmelt=newsnow(c)
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

  end subroutine UpdateSnowDepthAndFracClm5

  ! FIXME(wjs, 2019-07-25) For the old formulation (n&y) we'll just have the single line,
  !   int_snow(c) = int_snow(c) + newsnow(c)
  !-----------------------------------------------------------------------
  subroutine AddNewsnowToIntsnowClm5(this, bounds, num_c, filter_c, &
       newsnow, h2osno_total, frac_sno, &
       int_snow)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_c       ! number of columns in filter_c
    integer, intent(in) :: filter_c(:) ! column filter to operate over

    ! FIXME(wjs, 2019-07-22) document the following arguments
    real(r8), intent(in) :: newsnow( bounds%begc: )
    real(r8), intent(in) :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
    real(r8), intent(in) :: frac_sno( bounds%begc: )
    real(r8), intent(inout) :: int_snow( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    real(r8) :: temp_intsnow    ! temporary version of int_snow

    character(len=*), parameter :: subname = 'AddNewsnowToIntsnowClm5'
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

  end subroutine AddNewsnowToIntsnowClm5

  !-----------------------------------------------------------------------
  pure function FracSnowDuringMeltClm5(this, c, h2osno_total, int_snow) result(frac_sno)
    !
    ! !DESCRIPTION:
    ! Single-point function giving frac_snow during times when the snow pack is melting
    !
    ! !ARGUMENTS:
    real(r8) :: frac_sno  ! function result
    class(snow_cover_fraction_clm5_type), intent(in) :: this
    integer , intent(in) :: c            ! column we're operating on
    real(r8), intent(in) :: h2osno_total ! total snow water (mm H2O)
    real(r8), intent(in) :: int_snow     ! integrated snowfall (mm H2O)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: int_snow_limited  ! integrated snowfall, limited to be no greater than int_snow_max [mm]
    real(r8) :: smr

    character(len=*), parameter :: subname = 'FracSnowDuringMeltClm5'
    !-----------------------------------------------------------------------

    int_snow_limited = min(int_snow, this%int_snow_max)
    smr = min(1._r8, h2osno_total/int_snow_limited)

    frac_sno = 1. - (acos(min(1._r8,(2.*smr - 1._r8)))/rpi)**(this%n_melt(c))

  end function FracSnowDuringMeltClm5


end module SnowCoverFractionMod
