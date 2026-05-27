module debugMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Provides a place to define gridcells, landunits, columns, and/or patches that we want debug
  ! output for, as well as functions for checking whether a given g/l/c/p is being debugged.
  !
  ! !USES:
#include "shr_assert.h"
  use abortutils       , only : endrun
  use decompMod        , only : subgrid_level_gridcell, subgrid_level_landunit
  use decompMod        , only : subgrid_level_column, subgrid_level_patch
  use decompMod        , only : get_global_index
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER SUBROUTINES AND FUNCTIONS:
  public :: debugMod_init
  public :: do_debug_patch
  public :: do_debug_column
  public :: do_debug_landunit
  public :: do_debug_gridcell
  !
  ! !PRIVATE MEMBER DATA:

  ! Global indices for objects to debug.
  integer :: debug_p
  integer :: debug_c
  integer :: debug_l
  integer :: debug_g

  character(len=*), parameter, private :: sourcefile = __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine debugMod_init()
    !
    ! Handle namelist variables
    !
    ! !USES:
    use clm_varctl       , only : iulog
    use spmdMod          , only : masterproc, mpicom
    use controlMod       , only : NLFilename
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use shr_mpi_mod      , only : shr_mpi_bcast
    !
    ! !LOCAL VARIABLES:
    integer                 :: nu_nml                     ! unit for namelist file
    integer                 :: nml_error                  ! namelist i/o error flag
    integer :: nml_debug
    !-----------------------------------------------------------------------

    namelist /debug/                  &
         debug_p, &
         debug_c, &
         debug_l, &
         debug_g

    ! Default namelist variable values
    debug_p = -999
    debug_c = -999
    debug_l = -999
    debug_g = -999

    ! Read namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'debug', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=debug,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun('ERROR reading debug namelist')
          end if
       else
          call endrun('ERROR finding debug namelist')
       end if
       close(nu_nml)

       ! Print values
       write(iulog,*)
       write(iulog,*) 'Debug settings:'
       write(iulog, '(a,i8)') '  debug_p = ', debug_p
       write(iulog, '(a,i8)') '  debug_c = ', debug_c
       write(iulog, '(a,i8)') '  debug_l = ', debug_l
       write(iulog, '(a,i8)') '  debug_g = ', debug_g
    endif
    call shr_mpi_bcast(debug_p, mpicom)
    call shr_mpi_bcast(debug_c, mpicom)
    call shr_mpi_bcast(debug_l, mpicom)
    call shr_mpi_bcast(debug_g, mpicom)


  end subroutine debugMod_init

  !-----------------------------------------------------------------------
  function do_debug_patch(i_local, i_global) result(do_debug)
    !
    ! !DESCRIPTION:
    ! Return .true. if the given global patch index is the one we're debugging; .false. otherwise.
    ! Provide only i_local or i_global.
    !
    ! !ARGUMENTS:
    integer, optional :: i_local   ! Local patch index
    integer, optional :: i_global  ! Global patch index
    !
    ! !RESULT:
    logical :: do_debug
    integer :: i  ! Global patch index
    !-----------------------------------------------------------------------

    ! This makes it so that compilers should build this function as a simple boolean if we haven't
    ! requested any patch to be debugged, which reduces performance cost in production runs.
    if (debug_p < 1) then
       do_debug = .false.
       return
    end if

    if (present(i_global)) then
       call shr_assert(.not. present(i_local), file=sourcefile, line=__LINE__)
       call shr_assert(i_global >= 1, file=sourcefile, line=__LINE__)
       i = i_global
    else
       call shr_assert(present(i_local), file=sourcefile, line=__LINE__)
       i = get_global_index(subgrid_index=i_local, subgrid_level=subgrid_level_patch, &
            donot_abort_on_badindex=.false.)
    end if

    do_debug = i == debug_p
    return

  end function do_debug_patch

  !-----------------------------------------------------------------------
  function do_debug_column(i_local, i_global) result(do_debug)
    !
    ! !DESCRIPTION:
    ! Return .true. if the given global column index is the one we're debugging; .false. otherwise.
    ! Provide only i_local or i_global.
    !
    ! !ARGUMENTS:
    integer, optional :: i_local   ! Local column index
    integer, optional :: i_global  ! Global column index
    !
    ! !RESULT:
    logical :: do_debug
    integer :: i  ! Global column index
    !-----------------------------------------------------------------------

    ! This makes it so that compilers should build this function as a simple boolean if we haven't
    ! requested any column to be debugged, which reduces performance cost in production runs.
    if (debug_c < 1) then
       do_debug = .false.
       return
    end if

    if (present(i_global)) then
       call shr_assert(.not. present(i_local), file=sourcefile, line=__LINE__)
       call shr_assert(i_global >= 1, file=sourcefile, line=__LINE__)
       i = i_global
    else
       call shr_assert(present(i_local), file=sourcefile, line=__LINE__)
       i = get_global_index(subgrid_index=i_local, subgrid_level=subgrid_level_column, &
            donot_abort_on_badindex=.false.)
    end if

    do_debug = i == debug_c
    return

  end function do_debug_column

  !-----------------------------------------------------------------------
  function do_debug_landunit(i_local, i_global) result(do_debug)
    !
    ! !DESCRIPTION:
    ! Return .true. if the given global landunit index is the one we're debugging; .false. otherwise.
    ! Provide only i_local or i_global.
    !
    ! !ARGUMENTS:
    integer, optional :: i_local   ! Local landunit index
    integer, optional :: i_global  ! Global landunit index
    !
    ! !RESULT:
    logical :: do_debug
    integer :: i  ! Global landunit index
    !-----------------------------------------------------------------------

    ! This makes it so that compilers should build this function as a simple boolean if we haven't
    ! requested any landunit to be debugged, which reduces performance cost in production runs.
    if (debug_l < 1) then
       do_debug = .false.
       return
    end if

    if (present(i_global)) then
       call shr_assert(.not. present(i_local), file=sourcefile, line=__LINE__)
       call shr_assert(i_global >= 1, file=sourcefile, line=__LINE__)
       i = i_global
    else
       call shr_assert(present(i_local), file=sourcefile, line=__LINE__)
       i = get_global_index(subgrid_index=i_local, subgrid_level=subgrid_level_landunit, &
            donot_abort_on_badindex=.false.)
    end if

    do_debug = i == debug_l
    return

  end function do_debug_landunit

  !-----------------------------------------------------------------------
  function do_debug_gridcell(i_local, i_global) result(do_debug)
    !
    ! !DESCRIPTION:
    ! Return .true. if the given global gridcell index is the one we're debugging; .false. otherwise.
    ! Provide only i_local or i_global.
    !
    ! !ARGUMENTS:
    integer, optional :: i_local   ! Local gridcell index
    integer, optional :: i_global  ! Global gridcell index
    !
    ! !RESULT:
    logical :: do_debug
    integer :: i  ! Global gridcell index
    !-----------------------------------------------------------------------

    ! This makes it so that compilers should build this function as a simple boolean if we haven't
    ! requested any gridcell to be debugged, which reduces performance cost in production runs.
    if (debug_g < 1) then
       do_debug = .false.
       return
    end if

    if (present(i_global)) then
       call shr_assert(.not. present(i_local), file=sourcefile, line=__LINE__)
       call shr_assert(i_global >= 1, file=sourcefile, line=__LINE__)
       i = i_global
    else
       call shr_assert(present(i_local), file=sourcefile, line=__LINE__)
       i = get_global_index(subgrid_index=i_local, subgrid_level=subgrid_level_gridcell, &
            donot_abort_on_badindex=.false.)
    end if

    do_debug = i == debug_g
    return

  end function do_debug_gridcell

end module debugMod



