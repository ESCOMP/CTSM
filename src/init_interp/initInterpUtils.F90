module initInterpUtils

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains various utilities used by initInterp
  !
  ! !USES:

  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use ncdio_pio    , only: file_desc_t, ncd_inqdlen, ncd_inqdid, ncd_io
  use abortutils   , only: endrun

  implicit none
  private
  save

  ! Public methods

  public :: glc_elevclasses_are_same  ! Function that determines whether the glacier elevation classes are the same in the input and output files

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  logical function glc_elevclasses_are_same(ncidi, ncido)
    !
    ! !DESCRIPTION:
    ! Determines if the glacier elevation classes are the same in ncidi and ncido
    !
    ! Returns .true. if they are the same (i.e., same number and bounds of elevation
    ! classes, within roundoff), false if not.
    !
    ! !ARGUMENTS:
    type(file_desc_t) , intent(inout) :: ncidi
    type(file_desc_t) , intent(inout) :: ncido
    !
    ! !LOCAL VARIABLES:
    integer :: dimid_dummy
    logical :: dimexist
    integer :: glc_nec_input
    integer :: glc_nec_output
    logical :: readvar
    real(r8), pointer :: elevclass_bounds_input(:)
    real(r8), pointer :: elevclass_bounds_output(:)
    integer :: elevclass

    real(r8), parameter :: bounds_tol = 1.e-4_r8  ! tolerance for checking equality of elevclass bounds

    character(len=*), parameter :: subname = 'glc_elevclasses_are_same'
    !-----------------------------------------------------------------------

    ! BACKWARDS_COMPATIBILITY(wjs, 2018-03-19) Old restart files generated from
    ! configurations with istice rather than istice don't have a 'glc_nec' dimension.
    ! Users may still be using files generated like that. The value of this function
    ! should be irrelevant in that case. We can remove this code once we can rely on all
    ! users' finidat files having been generated from configurations with istice.
    call ncd_inqdid(ncidi, 'glc_nec', dimid_dummy, dimexist=dimexist)
    if (.not. dimexist) then
       glc_elevclasses_are_same = .false.
       return
    end if

    call ncd_inqdlen(ncido, dimid_dummy, glc_nec_output, name='glc_nec')
    call ncd_inqdlen(ncidi, dimid_dummy, glc_nec_input, name='glc_nec')

    if (glc_nec_input == glc_nec_output) then
       allocate(elevclass_bounds_input(0:glc_nec_input))
       allocate(elevclass_bounds_output(0:glc_nec_output))
       call ncd_io(ncid=ncido, varname='glc_elevclass_bounds', &
            data=elevclass_bounds_output, flag='read', readvar=readvar)
       if (.not. readvar) then
          call endrun('glc_elevclass_bounds not found on output file ' // &
               errMsg(sourcefile, __LINE__))
       end if
       call ncd_io(ncid=ncidi, varname='glc_elevclass_bounds', &
            data=elevclass_bounds_input, flag='read', readvar=readvar)
       if (.not. readvar) then
          ! BACKWARDS_COMPATIBILITY(wjs, 2018-03-19) Older restart files don't have this
          ! variable, but it's safe to assume that any old restart file was generated
          ! with the current elevation class bounds, as given below. Once we can rely on
          ! old restart files having the glc_elevclass_bounds variable, we should replace
          ! this hard-coded setting with a call to endrun, as we have for ncido.
          if (glc_nec_input == 10) then
             elevclass_bounds_input = [0._r8, 200._r8, 400._r8, 700._r8, 1000._r8, &
                  1300._r8, 1600._r8, 2000._r8, 2500._r8, 3000._r8, 10000._r8]
          else
             call endrun('glc_elevclass_bounds not found on input file ' // &
                  errMsg(sourcefile, __LINE__))
          end if
       end if

       glc_elevclasses_are_same = .true.
       do elevclass = 0, glc_nec_input
          if (abs(elevclass_bounds_input(elevclass) - elevclass_bounds_output(elevclass)) &
               > bounds_tol) then
             glc_elevclasses_are_same = .false.
          end if
       end do

       deallocate(elevclass_bounds_input)
       deallocate(elevclass_bounds_output)

    else  ! glc_nec_input /= glc_nec_output
       glc_elevclasses_are_same = .false.
    end if

  end function glc_elevclasses_are_same

end module initInterpUtils
