module CTSMForce2DStreamBaseType

  use shr_kind_mod , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use abortutils , only : endrun
  use decompMod , only : bounds_type
  use clm_varctl, only : FL => fname_len

  implicit none
  private

  type, abstract, public :: ctsm_force_2DStream_base_type
     private
     character(len=FL) :: stream_filename ! The stream data filename (also in sdat)
     character(len=CL) :: stream_name ! The stream name (also in sdat)
  contains

      ! PUBLIC METHODS
      procedure(Init_interface) , public, deferred :: Init
      procedure, public, non_overridable :: InitBase ! Initialize and read data in the streams, , store the g_to_ig index array
      procedure(Clean_interface), public, deferred :: Clean  ! Clean and deallocate the object class method
      procedure, public, non_overridable :: CleanBase ! Clean method for the base type
      procedure, public, non_overridable :: Advance   ! Advance the streams data to the current model date
      procedure, public :: GetPtr1D  ! Get pointer to the 1D data array
      procedure(Interp_interface), public, deferred :: Interp  ! method in extensions to turn stream data into output data

  end type ctsm_force_2DStream_base_type

  abstract interface

     subroutine Init_interface( this, bounds, fldfilename, meshfile, mapalgo, tintalgo, taxmode, &
                           year_first, year_last, model_year_align )
         ! Uses:
         use decompMod , only : bounds_type
         import :: ctsm_force_2DStream_base_type

         ! Arguments:
         class(ctsm_force_2DStream_base_type), intent(inout) :: this 
         type(bounds_type), intent(in) :: bounds
         character(*), intent(in) :: fldfilename ! stream data filename (full pathname) (single file)
         ! NOTE: fldfilename could be expanded to an array if needed, but currently we only have one file
         character(*), intent(in) :: meshfile ! full pathname to stream mesh file (none for global data)
         character(*), intent(in) :: mapalgo ! stream mesh -> model mesh mapping type
         character(*), intent(in) :: tintalgo ! time interpolation algorithm
         character(*), intent(in) :: taxMode ! time axis mode
         integer, intent(in) :: year_first ! first year to use
         integer, intent(in) :: year_last ! last  year to use
         integer, intent(in) :: model_year_align ! align yearFirst with this model year
     end subroutine Init_interface

     subroutine Clean_interface(this)
       ! Uses:
       import :: ctsm_force_2DStream_base_type
       !
       ! Arguments:
       class(ctsm_force_2DStream_base_type), intent(inout) :: this
     end subroutine Clean_interface 

     subroutine Interp_interface(this, bounds)
       ! Uses:
       use decompMod , only : bounds_type
       import :: ctsm_force_2DStream_base_type
       !
       ! Arguments:
       class(ctsm_force_2DStream_base_type), intent(inout) :: this
       type(bounds_type), intent(in) :: bounds
     end subroutine Interp_interface

  end interface

    character(len=*), parameter, private :: sourcefile = &
       __FILE__

   contains

      subroutine InitBase( this, bounds, varnames, fldfilename, meshfile, mapalgo, tintalgo, taxmode, name, &
                           year_first, year_last, model_year_align )
         ! Uses:
         ! Arguments:
         class(ctsm_force_2DStream_base_type), intent(inout) :: this 
         type(bounds_type), intent(in) :: bounds
         character(*), intent(in) :: varnames(:) ! variable names to read from stream file
         character(*), intent(in) :: fldfilename ! stream data filename (full pathname) (single file)
         ! NOTE: fldfilename could be expanded to an array if needed, but currently we only have one file
         character(*), intent(in) :: meshfile ! full pathname to stream mesh file (none for global data)
         character(*), intent(in) :: mapalgo ! stream mesh -> model mesh mapping type
         character(*), intent(in) :: tintalgo ! time interpolation algorithm
         character(*), intent(in) :: taxMode ! time axis mode
         character(*), intent(in) :: name ! name of stream
         integer, intent(in) :: year_first ! first year to use
         integer, intent(in) :: year_last ! last  year to use
         integer, intent(in) :: model_year_align ! align yearFirst with this model year

      end subroutine InitBase

      subroutine CleanBase( this )
         class(ctsm_force_2DStream_base_type) , intent(inout) :: this 

         call endrun('CTSMForce2DStreamBaseType: CleanBase method not implemented in stub')

      end subroutine CleanBase

      subroutine Advance(this)
         ! Arguments:
         class(ctsm_force_2DStream_base_type), intent(inout) :: this

         call endrun('CTSMForce2DStreamBaseType: Advance method not implemented in stub')

      end subroutine Advance

      subroutine GetPtr1D(this, fldname, dataptr1d)
         ! Get the pointer to the 1D data array for the given field name
         ! Uses:
         ! Arguments:
         class(ctsm_force_2DStream_base_type), intent(inout) :: this
         character(*), intent(in) :: fldname  ! field name to get pointer for
         real(r8), pointer :: dataptr1d(:)  ! Pointer to the 1D data

      end subroutine GetPtr1D

end module CTSMForce2DStreamBaseType
