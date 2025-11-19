module CTSMForce2DStreamBaseType

!
! Description:
!
! Base module to handle 2D streams in CTSM. Specific streams extend this object
! for the details needed to handle a specific stream file.
!
! Having this base type allows the ESMF specific streams implementation to be isolated
! from the CTSM code. This allows the streams code this is based on to change in one place.
! It also makes it easier to unit-test extensions of this type as they become pretty standard
! CTSM code and there is a unit-tester stub for this code.
!

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

  !-----------------------------------------------------------------------
  ! Base 2D streams type
  !-----------------------------------------------------------------------
  type, abstract, public :: ctsm_force_2DStream_base_type
     private
     type(shr_strdata_type) :: sdat  ! Stream data type
     character(len=FL) :: stream_filename ! The stream data filename (also in sdat)
     character(len=CL) :: stream_name ! The stream name (also in sdat)
  contains

     ! PUBLIC METHODS
     procedure(Init_interface) , public, deferred :: Init ! Initiale the extended type
     procedure, public, non_overridable :: InitBase ! Initialize and read data in the streams
     procedure(Clean_interface), public, deferred :: Clean  ! Clean and deallocate the object class method
     procedure, public, non_overridable :: CleanBase ! Clean method for the base type
     procedure, public, non_overridable :: Advance   ! Advance the streams data to the current model date
     procedure, public :: GetPtr1D  ! Get pointer to the 1D data array
     procedure, public :: Check1DPtrSize  ! Check that the data pointers are the expected size
     procedure(Interp_interface), public, deferred :: Interp  ! method in extensions to turn stream data into CTSM  data

  end type ctsm_force_2DStream_base_type
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Interfaces that will be deferred to the extended type
  !-----------------------------------------------------------------------
  abstract interface

     !-----------------------------------------------------------------------

     subroutine Init_interface( this, bounds, fldfilename, meshfile, mapalgo, tintalgo, taxmode, &
                           year_first, year_last, model_year_align )
         ! Description:
         !
         ! Initialize the specific stream type that extends the base type
         ! Normally the extended type will call the InitBase as well as doing other initialization needed
         !
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

     !-----------------------------------------------------------------------

     subroutine Clean_interface(this)
       ! Description:
       ! Clean up any memory allocated in the specific stream type that extends the base type
       ! Normally the extended type will call the CleanBase method as well as other things needed.
       ! Uses:
       import :: ctsm_force_2DStream_base_type
       !
       ! Arguments:
       class(ctsm_force_2DStream_base_type), intent(inout) :: this
     end subroutine Clean_interface

     !-----------------------------------------------------------------------

     subroutine Interp_interface(this, bounds)
       ! Description:
       ! Get the current time data from the streams and put it into the data of the extension.
       ! What this looks like may vary with the streams extension, but in general it will use
       ! The GetPtr1D method to get the streams data.
       ! Uses:
       use decompMod , only : bounds_type
       import :: ctsm_force_2DStream_base_type
       !
       ! Arguments:
       class(ctsm_force_2DStream_base_type), intent(inout) :: this
       type(bounds_type), intent(in) :: bounds
     end subroutine Interp_interface

  end interface
  !-----------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
     __FILE__

  !-----------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------

     !-----------------------------------------------------------------------

     subroutine InitBase( this, bounds, varnames, fldfilename, meshfile, mapalgo, tintalgo, taxmode, name, &
                          year_first, year_last, model_year_align )
         !
         ! Description:
         !
         ! Initialization of the base type. Extended types will normally call this as part of their initialization.
         !
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

     !-----------------------------------------------------------------------

     subroutine Check1DPtrSize( this, bounds )
         !
         ! Description:
         !
         ! Check that the stream data pointer size is as expected
         !
         use shr_kind_mod , only : CS => shr_kind_CS
         use dshr_stream_mod, only : shr_stream_streamType, shr_stream_getStreamFieldList
         use dshr_stream_mod, only : shr_stream_getMeshFileName
         use shr_log_mod , only : errMsg => shr_log_errMsg
         ! Arguments:
         class(ctsm_force_2DStream_base_type), intent(inout) :: this
         type(bounds_type), intent(in) :: bounds
         ! Local variables
         integer :: n ! Indices
         integer :: nvars ! Number of variables
         real(r8), pointer :: dataptr1d(:)  ! Pointer to the 1D data
         character(len=CS), allocatable :: varnames(:)
         type(shr_stream_streamType), pointer :: stream => NULL()
         character(len=CL) :: meshname

         ! Loop through the list of varnames

         stream => this%sdat%stream(1)
         nvars = stream%nvars
         allocate( varnames(nvars) )
         call shr_stream_getStreamFieldList( stream, varnames )
         call shr_stream_getMeshFileName( stream, meshname )
         do n = 1, nvars
             call this%GetPtr1D( varnames(n), dataptr1d )
             if ( trim(meshname) == 'none' ) then
                call shr_assert( size(dataptr1d) == 1, "Expect stream data to be 1 when no mesh given"//errMsg( file=sourcefile, line=__LINE__) )
             else
                call shr_assert( size(dataptr1d) == bounds%endg-bounds%begg + 1, "Expect stream data to be the size of grid bounds"//errMsg( file=sourcefile, line=__LINE__) )
             end if
         end do
         deallocate( varnames )

     end subroutine Check1DPtrSize

     !-----------------------------------------------------------------------

     subroutine CleanBase( this )
         ! Description:
         ! Clean up any memory in the base type as needed.
         ! Normally types that extend this base type will call this as part of their clean operation
         !
         ! Arguments:
         class(ctsm_force_2DStream_base_type) , intent(inout) :: this

         integer :: ierr ! error code

         ! Currently no data to deallocate other than the stream data type
         ! The stream data type doesn't have a clean method right now
         ! So doing a few things manually here
     end subroutine CleanBase

     !-----------------------------------------------------------------------

     subroutine Advance(this)
         !
         ! Description:
         !
         ! Advance the stream to the current time-step
         !
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

     !-----------------------------------------------------------------------

     subroutine GetPtr1D(this, fldname, dataptr1d)
         !
         ! Description:
         !
         ! Get the pointer to the 1D data array for the given field name
         ! Normally stream extensions will use this in the Interp method to
         ! save the stream data locally.
         !
         ! Uses:
         use dshr_methods_mod , only : dshr_fldbun_getfldptr
         ! Arguments:
         class(ctsm_force_2DStream_base_type), intent(inout) :: this
         character(*), intent(in) :: fldname  ! field name to get pointer for
         real(r8), pointer :: dataptr1d(:)  ! Pointer to the 1D data

         ! Local variables
         integer :: rc ! error return code

         ! Get pointer for stream data that is time and spatially interpolated to model time and grid
         call dshr_fldbun_getFldPtr(this%sdat%pstrm(1)%fldbun_model, fldname=fldname, fldptr1=dataptr1d, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=sourcefile)) then
            call endrun( 'Error getting field pointer for '//trim(fldname)//' from stream data', file=sourcefile, line=__LINE__ )
         end if

     end subroutine GetPtr1D

     !-----------------------------------------------------------------------

end module CTSMForce2DStreamBaseType
