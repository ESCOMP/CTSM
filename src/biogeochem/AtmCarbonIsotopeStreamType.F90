module AtmCarbonIsotopeStreamType


#include "shr_assert.h"
  !
  ! Description:
  !
  ! This extends the stream base type to implement streams for atmospheric
  ! Carbon isotope ratios that are read in from streams datasets (delta C13 and delta C14).
  !
  use shr_kind_mod , only : r8 => shr_kind_r8
  use clm_varctl , only : iulog
  use abortutils , only : endrun
  use decompMod , only : bounds_type
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use CTSMForce2DStreamBaseType, only : ctsm_force_2DStream_base_type

  implicit none
  private

  !-----------------------------------------------------------------------
  ! Atmospheric Delta C13 Stream Type
  !-----------------------------------------------------------------------
  character(len=*), parameter :: varname_c13 = 'delta13co2_in_air'
  type, public, extends(ctsm_force_2DStream_base_type) :: atm_delta_c13_stream_type
     private
     real(r8), public, allocatable :: atm_delta_c13(:) ! delta C13 data array
  contains

      ! Public Methods
      procedure, public :: C13Init ! C13 initialization
      procedure, public :: Init => C13Init ! Generic name for the initialization
      procedure, public :: C13Interp  ! C13 Interp method to fill the local data array
      procedure, public :: Interp => C13Interp  ! Generic name for the Interp method
      procedure, public :: C13ClassClean ! C13 clean method as a class method
      procedure, public :: Clean => C13ClassClean ! Generic name for the clean method
      final :: C13TypeClean  ! This clean method may be called by the compiler when the type goes out of scope
      ! Private methods
      procedure, private :: C13InitAllocate ! Allocate the local C13 data

  end type atm_delta_c13_stream_type

  !-----------------------------------------------------------------------
  ! Atmospheric Delta C14 Stream Type
  !-----------------------------------------------------------------------
  character(len=*), parameter :: varname_c14 = 'Delta14co2_in_air'
  type, public, extends(ctsm_force_2DStream_base_type) :: atm_delta_c14_stream_type
     private
     real(r8), public, allocatable :: atm_delta_c14(:) ! delta c14 data array
  contains

      ! Public Methods
      procedure, public :: C14Init ! C14 initialization
      procedure, public :: Init => C14Init ! Generic name for the initialization
      procedure, public :: C14Interp  ! C14 Interp method to fill the local data array
      procedure, public :: Interp => C14Interp  ! Generic name for the Interp method
      procedure, public :: C14ClassClean ! C14 clean method as a class method
      procedure, public :: Clean => C14ClassClean ! Generic name for the clean method
      final :: C14TypeClean ! This clean method may be called by the compiler when the type goes out of scope
      ! Private methods
      procedure, private :: C14InitAllocate ! Allocate the local C14 data

  end type atm_delta_c14_stream_type

  character(len=*), parameter, private :: sourcefile = &
  __FILE__

  !-----------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------

    !-------------------------------------------------------------------------------------

    subroutine C13Init( this, bounds, fldfilename, meshfile, mapalgo, tintalgo, taxmode, &
                        year_first, year_last, model_year_align )
         !
         ! Initialize the atmospheric delta C13 stream type
         !
         ! Arguments:
         class(atm_delta_c13_stream_type), intent(inout) :: this
         type(bounds_type), intent(in) :: bounds
         character(*), intent(in) :: fldfilename ! stream data filename (full pathname) (single file)
         character(*), intent(in) :: meshfile ! full pathname to stream mesh file (none for global data)
         character(*), intent(in) :: mapalgo ! stream mesh -> model mesh mapping type
         character(*), intent(in) :: tintalgo ! time interpolation algorithm
         character(*), intent(in) :: taxMode ! time axis mode
         integer, intent(in) :: year_first ! first year to use
         integer, intent(in) :: year_last ! last  year to use
         integer, intent(in) :: model_year_align ! align yearFirst with this model year

         ! Since C13 data is a single global value mapalgo and meshfile must both be none
         call shr_assert( trim(mapalgo) == "none", "mapalgo MUST be none for C13 streams"//errMsg( file=sourcefile, line=__LINE__) )
         call shr_assert( trim(meshfile) == "none", "meshfile MUST be none for C13 streams"//errMsg( file=sourcefile, line=__LINE__) )
         call this%InitBase( bounds, varnames = (/ varname_c13 /), fldfilename=fldfilename, meshfile=meshfile, &
                             mapalgo=mapalgo, tintalgo=tintalgo, taxmode=taxmode, name=varname_c13, &
                             year_first=year_first, year_last=year_last, model_year_align=model_year_align )
         call this%C13InitAllocate( bounds )
         call this%Advance( )
         call this%Check1DPtrSize( bounds )

    end subroutine C13Init

    !-------------------------------------------------------------------------------------

    subroutine C13InitAllocate( this, bounds )
         ! Allocate memory for the delta C13 data array
         use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
         class(atm_delta_c13_stream_type), intent(inout) :: this
         type(bounds_type), intent(in) :: bounds

         integer :: begg, endg

         begg = bounds%begg; endg = bounds%endg
         allocate( this%atm_delta_c13( bounds%begg : bounds%endg ) ); this%atm_delta_c13 = nan
    end subroutine C13InitAllocate

    !-------------------------------------------------------------------------------------

    subroutine C13ClassClean( this )
         ! Clean up memory for the C13 stream type as a class method
         class(atm_delta_c13_stream_type), intent(inout) :: this

         call C13TypeClean( this )

    end subroutine C13ClassClean

    !-------------------------------------------------------------------------------------

    subroutine C13TypeClean( this )
         ! Clean up memory for the C13 stream type for this specific type
         type(atm_delta_c13_stream_type), intent(inout) :: this

         deallocate( this%atm_delta_c13 )
         call this%CleanBase()

    end subroutine C13TypeClean

    !-------------------------------------------------------------------------------------

    subroutine C13Interp( this, bounds )
         !
         ! Fill the local CTSM grid delta C13 array with data from the stream
         !
         ! Arguments
         class(atm_delta_c13_stream_type), intent(inout) :: this
         type(bounds_type), intent(in) :: bounds

         ! Local Variables
         integer :: g
         real(r8), pointer :: dataptr1d(:)
         integer :: rc ! error return code

         ! Get pointer for stream data that is time and spatially interpolated to model time and grid
         call this%GetPtr1D( varname_c13, dataptr1d )

         do g = bounds%begg, bounds%endg
            this%atm_delta_c13(g) = dataptr1d(g)
         end do
    end subroutine C13Interp

    !-------------------------------------------------------------------------------------

    subroutine C14Init( this, bounds, fldfilename, meshfile, mapalgo, tintalgo, taxmode, &
                         year_first, year_last, model_year_align )
         !
         ! Initialize the atmospheric delta C14 stream type
         !
         ! Arguments:
         class(atm_delta_c14_stream_type), intent(inout) :: this
         type(bounds_type), intent(in) :: bounds
         character(*), intent(in) :: fldfilename ! stream data filename (full pathname) (single file)
         character(*), intent(in) :: meshfile ! full pathname to stream mesh file (none for global data)
         character(*), intent(in) :: mapalgo ! stream mesh -> model mesh mapping type
         character(*), intent(in) :: tintalgo ! time interpolation algorithm
         character(*), intent(in) :: taxMode ! time axis mode
         integer, intent(in) :: year_first ! first year to use
         integer, intent(in) :: year_last ! last  year to use
         integer, intent(in) :: model_year_align ! align yearFirst with this model year

         character(len=len(mapalgo)) :: str_mapalgo ! Temporary copy of mapalgo so can change input if meshfile is none

         ! Since C14 data has latitude bands mapalgo and meshfile can neither be set to none
         call shr_assert( trim(mapalgo) /= "none", "mapalgo MUST NOT be none for C14 streams"//errMsg( file=sourcefile, line=__LINE__) )
         ! TODO: Uncomment this error check when we are ready for the test for this to change answers
         !call shr_assert( trim(meshfile) /= "none", "meshfile MUST NOT be none for C14 streams"//errMsg( file=sourcefile, line=__LINE__) )
         ! TOD: Remove this tempoary bit at the same time
         if ( trim(meshfile) == "none" )then
            str_mapalgo = "none"
         else
            str_mapalgo = mapalgo
         end if
         ! TODO: End of temporary bit

         call this%InitBase( bounds, varnames = (/ varname_c14 /), fldfilename=fldfilename, meshfile=meshfile, &
                             mapalgo=str_mapalgo, tintalgo=tintalgo, taxmode=taxmode, name=varname_c14, &
                             year_first=year_first, year_last=year_last, model_year_align=model_year_align )
         call this%C14InitAllocate( bounds )
         call this%Advance( )
         call this%Check1DPtrSize( bounds )

    end subroutine C14Init

    !-------------------------------------------------------------------------------------

    subroutine C14InitAllocate( this, bounds )
         ! Allocate memory for the delta C14 data array
         use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
         ! Arguments
         class(atm_delta_c14_stream_type), intent(inout) :: this
         type(bounds_type), intent(in) :: bounds

         integer :: begg, endg

         begg = bounds%begg; endg = bounds%endg
         allocate( this%atm_delta_c14( bounds%begg : bounds%endg ) ); this%atm_delta_c14 = nan
    end subroutine C14InitAllocate

    !-------------------------------------------------------------------------------------

    subroutine C14Interp( this, bounds )
         ! Fill the local CTSM grid delta C13 array with data from the stream
         !
         ! Arguments:
         class(atm_delta_c14_stream_type), intent(inout) :: this
         type(bounds_type), intent(in) :: bounds

         ! Local Variables
         integer :: g
         real(r8), pointer :: dataptr1d(:)
         integer :: rc ! error return code

         ! Get pointer for stream data that is time and spatially interpolated to model time and grid
         call this%GetPtr1D( varname_c14, dataptr1d )

         do g = bounds%begg, bounds%endg
            this%atm_delta_c14(g) = dataptr1d(g)
         end do
    end subroutine C14Interp

    !-------------------------------------------------------------------------------------

    subroutine C14ClassClean( this )
         ! Clean up memory for the C14 stream type as a class method
         class(atm_delta_c14_stream_type), intent(inout) :: this

         call C14TypeClean( this )

    end subroutine C14ClassClean

    !-------------------------------------------------------------------------------------

    subroutine C14TypeClean( this )
         ! Clean up memory for the C14 stream type for this specific type
         type(atm_delta_c14_stream_type), intent(inout) :: this

         deallocate( this%atm_delta_c14 )
         call this%CleanBase()

    end subroutine C14TypeClean

    !-------------------------------------------------------------------------------------

end module AtmCarbonIsotopeStreamType