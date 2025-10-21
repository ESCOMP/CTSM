module AtmCarbonIsotopeStreamType
  use shr_kind_mod , only : r8 => shr_kind_r8
  use clm_varctl , only : iulog
  use abortutils , only : endrun
  use decompMod , only : bounds_type
  use CTSMForce2DStreamBaseType, only : ctsm_force_2DStream_base_type

  implicit none
  private

  type, public, extends(ctsm_force_2DStream_base_type) :: atm_delta_c13_stream_type
     private
     real(r8), allocatable :: atm_delta_c13(:) ! delta C13 data array
  contains

      ! Public Methods
      procedure, public :: C13Init
      procedure, public :: Init => C13Init
      procedure, public :: C13Interp 
      procedure, public :: Interp => C13Interp 
      procedure, public :: C13ClassClean
      procedure, public :: Clean => C13ClassClean
      final :: C13TypeClean
      ! Private methods
      procedure, private :: C13InitAllocate

  end type atm_delta_c13_stream_type

  type, public, extends(ctsm_force_2DStream_base_type) :: atm_delta_c14_stream_type
     private
     real(r8), allocatable :: atm_delta_c14(:) ! delta c14 data array
  contains

      ! Public Methods
      procedure, public :: C14Init
      procedure, public :: Init => C14Init
      procedure, public :: C14Interp 
      procedure, public :: Interp => C14Interp 
      procedure, public :: C14ClassClean
      procedure, public :: Clean => C14ClassClean
      final :: C14TypeClean
      ! Private methods
      procedure, private :: C14InitAllocate

  end type atm_delta_c14_stream_type

  contains

    subroutine C13Init( this, bounds, varnames, fldfilename, meshfile, mapalgo, tintalgo, taxmode, name, &
                        year_first, year_last, model_year_align )
         ! Uses:
         ! Arguments:
         class(atm_delta_c13_stream_type), intent(inout) :: this 
         type(bounds_type), intent(in) :: bounds
         character(*), intent(in) :: varnames(:) ! variable names to read from stream file
         character(*), intent(in) :: fldfilename ! stream data filename (full pathname) (single file)
         character(*), intent(in) :: meshfile ! full pathname to stream mesh file (none for global data)
         character(*), intent(in) :: mapalgo ! stream mesh -> model mesh mapping type
         character(*), intent(in) :: tintalgo ! time interpolation algorithm
         character(*), intent(in) :: taxMode ! time axis mode
         character(*), intent(in) :: name ! name of stream
         integer, intent(in) :: year_first ! first year to use
         integer, intent(in) :: year_last ! last  year to use
         integer, intent(in) :: model_year_align ! align yearFirst with this model year

         call this%InitBase( bounds, varnames, fldfilename, meshfile, mapalgo, tintalgo, taxmode, name, &
                           year_first, year_last, model_year_align )
         call this%C13InitAllocate( bounds )

     end subroutine C13Init

    subroutine C13InitAllocate( this, bounds )
         use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
         class(atm_delta_c13_stream_type), intent(inout) :: this 
         type(bounds_type), intent(in) :: bounds

         integer :: begg, endg

         begg = bounds%begg; endg = bounds%endg
         allocate( this%atm_delta_c13( bounds%begg : bounds%endg ) ); this%atm_delta_c13 = nan
    end subroutine C13InitAllocate

    subroutine C13ClassClean( this )
         class(atm_delta_c13_stream_type), intent(inout) :: this 

         call C13TypeClean( this )

    end subroutine C13ClassClean

    subroutine C13TypeClean( this )
         type(atm_delta_c13_stream_type), intent(inout) :: this 

         deallocate( this%atm_delta_c13 )
         call this%CleanBase()

    end subroutine C13TypeClean

    subroutine C13Interp( this, bounds )
         class(atm_delta_c13_stream_type), intent(inout) :: this 
         type(bounds_type), intent(in) :: bounds
    end subroutine C13Interp

     subroutine C14Init( this, bounds, varnames, fldfilename, meshfile, mapalgo, tintalgo, taxmode, name, &
                         year_first, year_last, model_year_align )
         ! Uses:
         ! Arguments:
         class(atm_delta_c14_stream_type), intent(inout) :: this 
         type(bounds_type), intent(in) :: bounds
         character(*), intent(in) :: varnames(:) ! variable names to read from stream file
         character(*), intent(in) :: fldfilename ! stream data filename (full pathname) (single file)
         character(*), intent(in) :: meshfile ! full pathname to stream mesh file (none for global data)
         character(*), intent(in) :: mapalgo ! stream mesh -> model mesh mapping type
         character(*), intent(in) :: tintalgo ! time interpolation algorithm
         character(*), intent(in) :: taxMode ! time axis mode
         character(*), intent(in) :: name ! name of stream
         integer, intent(in) :: year_first ! first year to use
         integer, intent(in) :: year_last ! last  year to use
         integer, intent(in) :: model_year_align ! align yearFirst with this model year

         call this%InitBase( bounds, varnames, fldfilename, meshfile, mapalgo, tintalgo, taxmode, name, &
                           year_first, year_last, model_year_align )
         call this%C14InitAllocate( bounds )

     end subroutine C14Init

    subroutine C14InitAllocate( this, bounds )
         use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
         class(atm_delta_c14_stream_type), intent(inout) :: this 
         type(bounds_type), intent(in) :: bounds

         integer :: begg, endg

         begg = bounds%begg; endg = bounds%endg
         allocate( this%atm_delta_c14( bounds%begg : bounds%endg ) ); this%atm_delta_c14 = nan
    end subroutine C14InitAllocate

    subroutine C14Interp( this, bounds )
         class(atm_delta_c14_stream_type), intent(inout) :: this 
         type(bounds_type), intent(in) :: bounds
    end subroutine C14Interp

    subroutine C14ClassClean( this )
         class(atm_delta_c14_stream_type), intent(inout) :: this 

         call C14TypeClean( this )

    end subroutine C14ClassClean


    subroutine C14TypeClean( this )
         type(atm_delta_c14_stream_type), intent(inout) :: this 

         deallocate( this%atm_delta_c14 )
         call this%CleanBase()

    end subroutine C14TypeClean

end module AtmCarbonIsotopeStreamType