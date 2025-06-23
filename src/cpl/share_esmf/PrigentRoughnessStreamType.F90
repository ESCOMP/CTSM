module PrigentRoughnessStreamType


  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in the Prigent et al. (1997) roughness length streams file
  ! Created by Danny M. Leung 22 Nov 2022
  ! !USES
  use ESMF
  use dshr_strdata_mod , only : shr_strdata_type
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use spmdMod          , only : mpicom, masterproc
  use clm_varctl       , only : iulog, FL => fname_len
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: prigent_roughness_stream_type
     real(r8), pointer, public :: prigent_rghn  (:)         ! Prigent et al. (1997) roughness length (m)
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: Init            ! Initialize and read data in
      procedure, public :: UseStreams      ! If Prigent rougness streams will be used
      procedure, public :: IsStreamInit    ! If the streams have been initialized and read in, so data can be used
      procedure, public :: Clean           ! Clean and deallocate the object

      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: InitAllocate   ! Allocate data

  end type prigent_roughness_stream_type

  ! ! PRIVATE DATA:
  type, private :: streamcontrol_type
     character(len=FL)  :: stream_fldFileName_prigentroughness   ! data Filename
     character(len=FL)  :: stream_meshfile_prigentroughness      ! mesh Filename
     character(len=CL)  :: prigentroughnessmapalgo               ! map algo
  contains
     procedure, private :: ReadNML      ! Read in control namelist
  end type streamcontrol_type

  type(streamcontrol_type), private :: control             ! Stream control data
  logical                 , private :: NMLRead = .false.   ! If namelist has been read
  logical                 , private :: InitDone = .false.  ! If initialization of streams has been done

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine Init(this, bounds, NLFilename)
   !
   ! Initialize the prigent roughness stream object
   !
   ! Uses:
   use spmdMod          , only : iam
   use lnd_comp_shr     , only : mesh, model_clock
   use dshr_strdata_mod , only : shr_strdata_init_from_inline, shr_strdata_print
   use dshr_strdata_mod , only : shr_strdata_advance
   use dshr_methods_mod , only : dshr_fldbun_getfldptr
   !
   ! arguments
   implicit none
   class(prigent_roughness_stream_type) :: this
   type(bounds_type), intent(in)   :: bounds
   character(len=*),  intent(in)   :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer                        :: ig, g, n           ! Indices
   integer                        :: year               ! year (0, ...) for nstep+1
   integer                        :: mon                ! month (1, ..., 12) for nstep+1
   integer                        :: day                ! day of month (1, ..., 31) for nstep+1
   integer                        :: sec                ! seconds into current date for nstep+1
   integer                        :: mcdate             ! Current model date (yyyymmdd)
   type(shr_strdata_type)         :: sdat_rghn      ! input data stream
   character(len=16), allocatable :: stream_varnames(:) ! array of stream field names
   integer                        :: rc                 ! error code
   real(r8), pointer              :: dataptr1d(:)       ! temporary pointer
   character(len=*), parameter    :: stream_name = 'prigent_roughness'
   character(len=*), parameter    :: subname     = 'PrigentRoughnessStream::Init'
   !-----------------------------------------------------------------------

      call control%ReadNML( bounds, NLFileName )

      call this%InitAllocate( bounds )

      if ( this%useStreams() )then


         allocate(stream_varnames(1))
         stream_varnames = (/"Z0a"/)  ! varname in cdf5_Z0a_Prigent-Globe-025x025-09262022.nc, in centimeter

         if (masterproc) then
            write(iulog,*) '  stream_varnames                  = ',stream_varnames
         end if

         ! Initialize the cdeps data type sdat_rghn
         call shr_strdata_init_from_inline(sdat_rghn,                                   &
              my_task             = iam,                                                &
              logunit             = iulog,                                              &
              compname            = 'LND',                                              &
              model_clock         = model_clock,                                        &
              model_mesh          = mesh,                                               &
              stream_meshfile     = control%stream_meshfile_prigentroughness,           &
              stream_lev_dimname  = 'null',                                             &
              stream_mapalgo      = control%prigentroughnessmapalgo,                    &
              stream_filenames    = (/trim(control%stream_fldFileName_prigentroughness)/), &
              stream_fldlistFile  = stream_varnames,                                    &
              stream_fldListModel = stream_varnames,                                    &
              stream_yearFirst    = 1997,                                               &
              stream_yearLast     = 1997,                                               &
              stream_yearAlign    = 1,                                                  &
              stream_offset       = 0,                                                  &
              stream_taxmode      = 'extend',                                           &
              stream_dtlimit      = 1.0e30_r8,                                          &
              stream_tintalgo     = 'linear',                                           &
              stream_name         = 'Prigent roughness',                                &
              rc                  = rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! Explicitly set current date to a hardcoded constant value. Otherwise
         ! using the real date can cause roundoff differences that are
         ! detrected as issues with exact restart.  EBK M05/20/2017
         ! call get_curr_date(year, mon, day, sec)
         year = 1997
         mon  = 12
         day  = 31
         sec  = 0
         mcdate = year*10000 + mon*100 + day

         call shr_strdata_advance(sdat_rghn, ymd=mcdate, tod=sec, logunit=iulog, istr='prigentrghn', rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! Get pointer for stream data that is time and spatially interpolate to model time and grid
         do n = 1,size(stream_varnames)
            call dshr_fldbun_getFldPtr(sdat_rghn%pstrm(1)%fldbun_model, stream_varnames(n), fldptr1=dataptr1d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
            end if
            if (trim(stream_varnames(n)) == 'Z0a') then
               ig = 0
               do g = bounds%begg,bounds%endg
                  ig = ig+1
                  this%prigent_rghn(g) = dataptr1d(ig)
               end do

            end if

         end do
         deallocate(stream_varnames)
         InitDone = .true.
      end if

  end subroutine Init

  !==============================================================================
  logical function UseStreams(this)
    !
    ! !DESCRIPTION:
    ! Return true if the Prigent Roughness stream is being used
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(prigent_roughness_stream_type) :: this
    !
    character(len=*), parameter :: subname = 'PrigentRoughnessStream::UseStreams'
    !
    if ( this%IsStreamInit() )then
       if ( trim(control%stream_fldFileName_prigentroughness) == '' )then
          UseStreams = .false.  ! Prigent streams are off without a filename given
       else
          UseStreams = .true.
       end if
    else
       UseStreams = .true.
    end if
  end function UseStreams

  !==============================================================================
  logical function IsStreamInit(this)
    !
    ! !DESCRIPTION:
    ! Return true if the streams have been initialized
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(prigent_roughness_stream_type) :: this
    !
    character(len=*), parameter :: subname = 'PrigentRoughnessStream::IsStreamInit'
    !
    if ( .not. NMLRead )then
       call endrun(msg=subname//' ERROR Namelist has NOT been read first, call Init before this')
    end if
    if ( InitDone )then
       IsStreamInit = .true.
    else
       IsStreamInit = .false.
    end if
  end function IsStreamInit

  !==============================================================================
  subroutine Clean(this)
   !
   ! Deallocate and clean the object
   !
   ! Uses:
   !
   ! arguments
   implicit none
   class(prigent_roughness_stream_type) :: this
   !
   ! local variables
   !-----------------------------------------------------------------------
   deallocate(this%prigent_rghn)
   this%prigent_rghn => NULL()
   InitDone = .false.

  end subroutine Clean

  !==============================================================================
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    implicit none
    class(prigent_roughness_stream_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    if ( this%useStreams() )then
       allocate(this%prigent_rghn(begg:endg))
    else
       allocate(this%prigent_rghn(0))
    end if
    this%prigent_rghn(:) = nan

  end subroutine InitAllocate

  !==============================================================================
  subroutine ReadNML(this, bounds, NLFilename)
   !
   ! Read the namelist data stream information.
   !
   ! Uses:
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   !
   ! arguments
   implicit none
   class(streamcontrol_type)     :: this
   type(bounds_type), intent(in) :: bounds
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   logical            :: use_prigent_roughness = .true.
   character(len=FL)  :: stream_fldFileName_prigentroughness = ' '
   character(len=FL)  :: stream_meshfile_prigentroughness = ' '
   character(len=CL)  :: prigentroughnessmapalgo = 'bilinear'
   character(len=*), parameter :: namelist_name = 'prigentroughness'    ! MUST agree with group name in namelist definition to read.
   character(len=*), parameter :: subName = "('prigentroughness::ReadNML')"
   !-----------------------------------------------------------------------

   namelist /prigentroughness/ &               ! MUST agree with namelist_name above
        prigentroughnessmapalgo,  stream_fldFileName_prigentroughness, stream_meshfile_prigentroughness, &
        use_prigent_roughness

   ! Default values for namelist

   ! Read prigentroughness namelist
   if (masterproc) then
      open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, namelist_name, status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=prigentroughness,iostat=nml_error)   ! MUST agree with namelist_name above
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
   endif

   call shr_mpi_bcast(use_prigent_roughness               , mpicom)
   call shr_mpi_bcast(prigentroughnessmapalgo             , mpicom)
   call shr_mpi_bcast(stream_fldFileName_prigentroughness , mpicom)
   call shr_mpi_bcast(stream_meshfile_prigentroughness    , mpicom)

   ! Error checking
   if ( .not. use_prigent_roughness )then
      if ( len_trim(stream_fldFileName_prigentroughness) /= 0 )then
         call endrun(msg=' ERROR stream_fldFileName_prigentroughness is set, but use_prigent_roughness is FALSE' &
                     //errMsg(sourcefile, __LINE__))
      end if
      if ( len_trim(stream_meshfile_prigentroughness) /= 0 )then
         call endrun(msg=' ERROR stream_meshfile_prigentroughness is set, but use_prigent_roughness is FALSE' &
                     //errMsg(sourcefile, __LINE__))
      end if
   else
      if ( len_trim(stream_fldFileName_prigentroughness) == 0 )then
         call endrun(msg=' ERROR stream_fldFileName_prigentroughness is NOT set, but use_prigent_roughness is TRUE' &
                     //errMsg(sourcefile, __LINE__))
      end if
      if ( len_trim(stream_meshfile_prigentroughness) == 0 )then
         call endrun(msg=' ERROR stream_meshfile_prigentroughness is NOT set, but use_prigent_roughness is TRUE' &
                    //errMsg(sourcefile, __LINE__))
      end if
   end if

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) namelist_name, ' stream settings:'
      write(iulog,*) '  use_prigent_roughness               = ',use_prigent_roughness
      if ( use_prigent_roughness )then
         write(iulog,*) '  stream_fldFileName_prigentroughness = ',trim(stream_fldFileName_prigentroughness)
         write(iulog,*) '  stream_meshfile_prigentroughness    = ',trim(stream_meshfile_prigentroughness)
         write(iulog,*) '  prigentroughnessmapalgo             = ',trim(prigentroughnessmapalgo)
      end if
   endif
   this%stream_fldFileName_prigentroughness = stream_fldFileName_prigentroughness
   this%stream_meshfile_prigentroughness    = stream_meshfile_prigentroughness
   this%prigentroughnessmapalgo             = prigentroughnessmapalgo

   ! Mark namelist read as having been done
   NMLRead = .true.

 end subroutine ReadNML

end module PrigentRoughnessStreamType
