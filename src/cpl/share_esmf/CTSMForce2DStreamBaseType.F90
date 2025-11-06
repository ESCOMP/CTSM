module CTSMForce2DStreamBaseType

#include "shr_assert.h"

  use ESMF, only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
  use dshr_strdata_mod , only : shr_strdata_type 
  use shr_kind_mod , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use clm_varctl , only : iulog
  use spmdMod , only : masterproc, mpicom, iam
  use abortutils , only : endrun
  use decompMod , only : bounds_type
  use clm_varctl, only : FL => fname_len

  implicit none
  private

  type, abstract, public :: ctsm_force_2DStream_base_type
     private
     type(shr_strdata_type), public :: sdat  ! Stream data type
     character(len=FL) :: stream_filename ! The stream data filename (also in sdat)
     character(len=CL) :: stream_name ! The stream name (also in sdat)
  contains

      ! PUBLIC METHODS
      procedure(Init_interface) , public, deferred :: Init
      procedure, public, non_overridable :: InitBase ! Initialize and read data in the streams, , store the g_to_ig index array
      procedure(Clean_interface), public, deferred :: Clean  ! Clean and deallocate the object class method
      procedure, public, non_overridable :: CleanBase ! Clean method for the base type
      procedure, public, non_overridable :: Advance   ! Advance the streams data to the current model date
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
         use lnd_comp_shr , only : mesh, model_clock
         use dshr_strdata_mod , only : shr_strdata_init_from_inline
         use decompMod , only : bounds_level_proc
         use shr_log_mod , only : errMsg => shr_log_errMsg
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

         ! Local variables
         integer, parameter :: offset = 0 ! time offset in seconds of stream data
         integer :: rc ! error return code

         ! Some error checking...
         SHR_ASSERT( bounds%level == bounds_level_proc, "InitBase should have a processor bounds, so we can do some checking"//errMsg( sourcefile, __LINE__) )
         SHR_ASSERT( bounds%begg == 1, "Make sure the starting bounds index is 1 so we know the mapping to gridcells is correct"//errMsg( sourcefile, __LINE__) )
         if ( len(fldfilename) >= FL )then
            call endrun( 'stream field filename is too long:'//trim(fldfilename), file=sourcefile, line=__LINE__ )
         end if
         this%stream_filename = fldfilename
         this%stream_name = name
         call shr_strdata_init_from_inline(this%sdat, &
               my_task             = iam, &
               logunit             = iulog, &
               compname            = 'LND', &
               model_clock         = model_clock,&
               model_mesh          = mesh, &
               stream_meshfile     = trim(meshfile), &
               stream_lev_dimname  = 'null', &
               stream_mapalgo      = mapalgo, &
               stream_filenames    = (/trim(fldfilename)/), &
               stream_fldlistFile  = varnames, &
               stream_fldListModel = varnames, &
               stream_yearFirst    = year_first, &
               stream_yearLast     = year_last, &
               stream_yearAlign    = model_year_align, &
               stream_offset       = offset, &
               stream_taxmode      = taxmode, &
               stream_dtlimit      = 1.0e30_r8, &
               stream_tintalgo     = tintalgo, &
               stream_name         = name, &
               rc                  = rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=sourcefile)) then
            write(iulog,*) ' Streams initialization failing for ', trim(name), ' stream file = ', trim(fldfilename)
            call endrun( 'CTSM forcing Streams initialization failing', file=sourcefile, line=__LINE__ )
         end if
      end subroutine InitBase

      subroutine CleanBase( this )
         class(ctsm_force_2DStream_base_type) , intent(inout) :: this 

         integer :: ierr ! error code

         ! Currently no data to deallocate other than the stream data type
         ! The stream data type doesn't have a clean method right now
         ! So doing a few things manually here
      end subroutine CleanBase

      subroutine Advance(this)
         ! Uses:
         use clm_time_manager , only : get_curr_date
         use dshr_strdata_mod , only : shr_strdata_advance
         !
         ! Arguments:
         class(ctsm_force_2DStream_base_type), intent(inout) :: this
         ! !LOCAL VARIABLES:
         integer :: year    ! year (0, ...) for nstep+1
         integer :: mon     ! month (1, ..., 12) for nstep+1
         integer :: day     ! day of month (1, ..., 31) for nstep+1
         integer :: sec     ! seconds into current date for nstep+1
         integer :: mcdate  ! Current model date (yyyymmdd)
         integer :: rc      ! Error return code

         ! Advance sdat stream
         call get_curr_date(year, mon, day, sec)
         mcdate = year*10000 + mon*100 + day
         call shr_strdata_advance(this%sdat, ymd=mcdate, tod=sec, logunit=iulog, istr='CTSMForce2DStreamBase', rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            write(iulog,*) ' Streams advance failing for ', trim(this%stream_name), ' stream file = ', trim(this%stream_filename)
            call endrun( 'CTSM forcing Streams advance failing', file=sourcefile, line=__LINE__ )
         end if
      end subroutine Advance

end module CTSMForce2DStreamBaseType
