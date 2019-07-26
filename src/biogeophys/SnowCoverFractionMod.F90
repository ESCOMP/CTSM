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
  use clm_varcon     , only : rpi
  use clm_varctl     , only : iulog, use_subgrid_fluxes
  use ncdio_pio      , only : file_desc_t
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: CreateAndInitScfMethod
  ! FIXME(wjs, 2019-07-26) this specific type should not be public: only the generic type
  ! should be public
  public :: snow_cover_fraction_clm5_type

  type :: snow_cover_fraction_clm5_type
     private
     real(r8) :: int_snow_max = -999._r8  ! limit applied to integrated snowfall when determining changes in snow-covered fraction during melt (mm H2O)
     real(r8) :: accum_factor = -999._r8  ! Accumulation constant for fractional snow covered area (unitless)
   contains
     procedure, public :: Init => InitClm5
     procedure, public :: UpdateSnowDepthAndFracClm5
     procedure, public :: AddNewsnowToIntsnowClm5
     procedure, public :: FracSnowDuringMeltClm5

     procedure, private :: ReadNamelistClm5
     procedure, private :: ReadParamsClm5
     procedure, private :: CheckValidInputsClm5
  end type snow_cover_fraction_clm5_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
contains

  ! ========================================================================
  ! Factory method
  ! ========================================================================

  !-----------------------------------------------------------------------
  function CreateAndInitScfMethod(NLFilename, params_ncid) result(scf_method)
    !
    ! !DESCRIPTION:
    ! Create an instance of the appropriate snow_cover_fraction_type, initialize it and
    ! return it
    !
    ! !ARGUMENTS:
    ! FIXME(wjs, 2019-07-26) Change this to generic snow_cover_fraction_type
    class(snow_cover_fraction_clm5_type), allocatable :: scf_method  ! function result
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
         NLFilename  = NLFilename, &
         params_ncid = params_ncid)

  end function CreateAndInitScfMethod

  ! ========================================================================
  ! Methods for CLM5 parameterization
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine InitClm5(this, NLFilename, params_ncid)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    type(file_desc_t), intent(inout) :: params_ncid ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitClm5'
    !-----------------------------------------------------------------------

    call this%ReadNamelistClm5(NLFilename)
    call this%ReadParamsClm5(params_ncid)
    call this%CheckValidInputsClm5()

  end subroutine InitClm5

  !-----------------------------------------------------------------------
  subroutine ReadNamelistClm5(this, NLFilename)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    ! local variables corresponding to the namelist variables we read here
    real(r8) :: int_snow_max

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: nmlname = 'scf_clm5_inparm'

    character(len=*), parameter :: subname = 'ReadNamelistClm5'
    !-----------------------------------------------------------------------

    namelist /scf_clm5_inparm/ int_snow_max

    int_snow_max = nan

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

    this%int_snow_max = int_snow_max

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname // ' settings:'
       write(iulog, nml=scf_clm5_inparm)
       write(iulog,*) ' '
    end if

  end subroutine ReadNamelistClm5

  !-----------------------------------------------------------------------
  subroutine ReadParamsClm5(this, params_ncid)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use paramUtilMod, only: readNcdioScalar
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(inout) :: this
    type(file_desc_t), intent(inout) :: params_ncid  ! pio netCDF file id for parameter file
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'ReadParamsClm5'
    !-----------------------------------------------------------------------

    ! Accumulation constant for fractional snow covered area (unitless)
    call readNcdioScalar(params_ncid, 'accum_factor', subname, this%accum_factor)

  end subroutine ReadParamsClm5

  !-----------------------------------------------------------------------
  subroutine CheckValidInputsClm5(this)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS:
    class(snow_cover_fraction_clm5_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CheckValidInputsClm5'
    !-----------------------------------------------------------------------

    if (this%int_snow_max <= 0.0_r8) then
       write(iulog,*)'ERROR: int_snow_max = ', this%int_snow_max,' is not supported, must be greater than 0.0.'
       call endrun(msg=' ERROR: invalid value for int_snow_max in CLM namelist. '//&
            errMsg(sourcefile, __LINE__))
    end if

  end subroutine CheckValidInputsClm5


  !-----------------------------------------------------------------------
  subroutine UpdateSnowDepthAndFracClm5(this, bounds, num_c, filter_c, &
       urbpoi, n_melt, h2osno_total, snowmelt, int_snow, newsnow, bifall, &
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
    ! FIXME(wjs, 2019-07-25) Move n_melt to being class-local
    real(r8), intent(in) :: n_melt( bounds%begc: )
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
    SHR_ASSERT_FL((ubound(n_melt, 1) == bounds%endc), sourcefile, __LINE__)
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
             frac_sno(c) = this%FracSnowDuringMeltClm5( &
                  h2osno_total = h2osno_total(c), &
                  int_snow = int_snow(c), &
                  n_melt = n_melt(c))
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
       newsnow, h2osno_total, frac_sno, n_melt, &
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
    ! FIXME(wjs, 2019-07-25) Move n_melt to being class-local
    real(r8), intent(in) :: n_melt( bounds%begc: )
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
    SHR_ASSERT_FL((ubound(n_melt, 1) == bounds%endc), sourcefile, __LINE__)
    SHR_ASSERT_FL((ubound(int_snow, 1) == bounds%endc), sourcefile, __LINE__)

    do fc = 1, num_c
       c = filter_c(fc)

       if (newsnow(c) > 0._r8) then
          ! reset int_snow after accumulation events: make int_snow consistent with new
          ! fsno, h2osno_total
          temp_intsnow= (h2osno_total(c) + newsnow(c)) &
               / (0.5*(cos(rpi*(1._r8-max(frac_sno(c),1e-6_r8))**(1./n_melt(c)))+1._r8))
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
  pure function FracSnowDuringMeltClm5(this, h2osno_total, int_snow, n_melt) result(frac_sno)
    !
    ! !DESCRIPTION:
    ! Single-point function giving frac_snow during times when the snow pack is melting
    !
    ! !ARGUMENTS:
    real(r8) :: frac_sno  ! function result
    class(snow_cover_fraction_clm5_type), intent(in) :: this
    real(r8), intent(in) :: h2osno_total ! total snow water (mm H2O)
    real(r8), intent(in) :: int_snow     ! integrated snowfall (mm H2O)
    ! FIXME(wjs, 2019-07-25) move this into the class
    real(r8), intent(in) :: n_melt
    !
    ! !LOCAL VARIABLES:
    real(r8) :: int_snow_limited  ! integrated snowfall, limited to be no greater than int_snow_max [mm]
    real(r8) :: smr

    character(len=*), parameter :: subname = 'FracSnowDuringMeltClm5'
    !-----------------------------------------------------------------------

    int_snow_limited = min(int_snow, this%int_snow_max)
    smr = min(1._r8, h2osno_total/int_snow_limited)

    frac_sno = 1. - (acos(min(1._r8,(2.*smr - 1._r8)))/rpi)**(n_melt)

  end function FracSnowDuringMeltClm5


end module SnowCoverFractionMod
