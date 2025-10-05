module SelfTestDriver

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains a top-level driver to the self-test code.
  !
  ! See the README file in this directory for a high-level overview of these self-tests.

  use decompMod, only : bounds_type
  use TestNcdioPio, only : test_ncdio_pio
  use abortutils, only : endrun
  use clm_varctl, only : iulog

  implicit none
  private
  save

  ! Public routines

  public :: self_test_driver ! Run the self-tests asked for
  public :: self_test_readnml  ! Read in the general self testing options for overall code flow
  public :: for_testing_bypass_init_after_self_tests ! For testing bypass the rest of the initialization after the self test driver was run
  public :: for_testing_bypass_run_except_clock_advance ! For testing bypass most of the run phase other than the clock advance

  ! Private module data
  logical :: for_testing_bypass_init ! For testing bypass the initialization phase after the self-test driver
  logical :: for_testing_bypass_run ! For testing bypass most of the run phase except the time advance
  logical :: for_testing_run_ncdiopio_tests  ! true => run tests of ncdio_pio
  logical :: for_testing_run_decomp_init_tests ! true => run tests of decompInit
  logical, public :: for_testing_exit_after_self_tests ! true => exit after running self tests

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine self_test_driver(bounds)
    !
    ! !DESCRIPTION:
    ! Top-level driver to the self-test code
    !
    ! This subroutine should be called all the time, but each set of self tests is only
    ! run if the appropriate flag is set.
    !
    ! !USES:
    use decompMod, only : bounds_type
    use TestNcdioPio, only : test_ncdio_pio
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_Finalize
    use shr_sys_mod, only : shr_sys_flush
    use spmdMod, only : masterproc
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'self_test_driver'
    integer :: ntests = 0
    !-----------------------------------------------------------------------

    if (for_testing_run_ncdiopio_tests) then
       ntests = ntests + 1
       call test_ncdio_pio(bounds)
    end if
    if (for_testing_run_decomp_init_tests) then
       ntests = ntests + 1
    end if
    if (for_testing_exit_after_self_tests) then
       ! Print out some messaging if we are exiting after self tests.
       if ( masterproc ) then 
          if ( ntests == 0 )then
             write(iulog,*) 'WARNING: You are exiting after self tests were run -- but no self tests were run.'
          else
             write(iulog,*) 'Exiting after running ', ntests, ' self tests.'
          end if
          call shr_sys_flush(iulog)
          call ESMF_LogWrite(' exiting after running self tests', ESMF_LOGMSG_INFO)
       end if
    end if

  end subroutine self_test_driver

  !-----------------------------------------------------------------------
  subroutine self_test_readnml(NLFileName)
    !
    ! !DESCRIPTION:
    ! Namelist read for the self-test driver. This includes bypass options
    ! that will be used in other parts of the code to bypass bits of the code
    ! for testing purposes.
    !
    ! !USES:
    use shr_nl_mod , only : shr_nl_find_group_name
    use spmdMod, only : masterproc, mpicom
    use shr_mpi_mod, only : shr_mpi_bcast
    use clm_varctl, only : iulog
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename  ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr ! error code
    integer :: unitn ! unit for namelist file

    ! Namelist name: this has to be matched with the name in the read stqatement
    character(len=*), parameter :: nmlname = 'for_testing_options'
    !-----------------------------------------------------------------------

    namelist /for_testing_options/ for_testing_bypass_init, for_testing_bypass_run, &
         for_testing_run_ncdiopio_tests, for_testing_run_decomp_init_tests, &
         for_testing_exit_after_self_tests

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       write(iulog,*) 'Read in '//nmlname//' namelist'
       open(newunit=unitn, status='old', file=NLFilename)
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unit=unitn, nml=for_testing_options, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist", file=sourcefile, line=__LINE__)
          end if
       else
          call endrun(msg="ERROR finding "//nmlname//"namelist", file=sourcefile, line=__LINE__)
       end if
       close(unitn)
    end if

    call shr_mpi_bcast (for_testing_bypass_init, mpicom)
    call shr_mpi_bcast (for_testing_bypass_run, mpicom)
    call shr_mpi_bcast(for_testing_run_ncdiopio_tests, mpicom)
    call shr_mpi_bcast(for_testing_run_decomp_init_tests, mpicom)
    call shr_mpi_bcast(for_testing_exit_after_self_tests, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=for_testing_options)
       write(iulog,*) ' '
    end if

  end subroutine self_test_readnml

  !-----------------------------------------------------------------------

  logical function for_testing_bypass_init_after_self_tests()
    ! Determine if should exit initialization early after having run the self tests
    if ( for_testing_bypass_init ) then
       for_testing_bypass_init_after_self_tests = .true.
    else
       for_testing_bypass_init_after_self_tests = .false.
    end if
  end function for_testing_bypass_init_after_self_tests

  !-----------------------------------------------------------------------

  logical function for_testing_bypass_run_except_clock_advance()
    ! Determine if should skip most of the run phase other than the clock advance
    if ( for_testing_bypass_init ) then
       for_testing_bypass_run_except_clock_advance = .true.
    else
       for_testing_bypass_run_except_clock_advance = .false.
    end if
  end function for_testing_bypass_run_except_clock_advance

  !-----------------------------------------------------------------------

end module SelfTestDriver
