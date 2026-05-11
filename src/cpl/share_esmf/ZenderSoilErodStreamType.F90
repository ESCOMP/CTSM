module ZenderSoilErodStreamType
#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in the Zender et al. (2003b) Dust source function streams file that has been read in from CAM instead of CLM. dmleung 11 Mar 2023
  ! pathname in CAM: /glade/p/cesmdata/cseg/inputdata/atm/cam/dst/
  ! relevant filenames: (CAM6) dst_source2x2tunedcam6-2x2-04062017.nc (default)
  !                     (CAM5) dst_source2x2_cam5.4_c150327.nc
  !                     (CAM4) dst_source2x2tuned-cam4-06132012.nc
  ! These files are largely similar and the differences are mainly only a little tuning.
  ! This .F90 file for now only deals with the CAM6 source function, which can be used in CAM5 and CAM4 too. Not sure if we will expand the code to include a namelist control to deal
  ! with the other files.
  !
  ! !USES
  use ESMF             , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_Finalize, ESMF_END_ABORT
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

  type, public :: soil_erod_stream_type
     real(r8), pointer, private :: soil_erodibility  (:)         ! Zender et al. (2003b) dust source function (or soil erodibility)
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: Init            ! Initialize and read data in
      procedure, public :: CalcDustSource ! Calculate dust source spatial filter (basically truncating stream data value smaller than 0.1 following CAM's practice) based on input streams
      procedure, public :: UseStreams      ! If streams will be used

      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: InitAllocate   ! Allocate data

  end type soil_erod_stream_type

  ! ! PRIVATE DATA:
  type, private :: streamcontrol_type
     character(len=FL)  :: stream_fldFileName_zendersoilerod   ! data Filename
     character(len=FL)  :: stream_meshfile_zendersoilerod      ! mesh Filename
     character(len=CL)  :: zendersoilerod_mapalgo              ! map algo
     logical            :: namelist_set = .false.              ! if namelist was set yet
  contains
     procedure, private :: ReadNML     ! Read in namelist
  end type streamcontrol_type

  type(streamcontrol_type), private :: control        ! Stream control data

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine Init(this, bounds, NLFilename)
   !
   ! Initialize the Zender soil eroditability stream object
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
   class(soil_erod_stream_type) :: this
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
   type(shr_strdata_type)         :: sdat_erod      ! input data stream
   character(len=16), allocatable :: stream_varnames(:) ! array of stream field names
   integer                        :: rc                 ! error code
   real(r8), pointer              :: dataptr1d(:)       ! temporary pointer
   character(len=*), parameter    :: stream_name = 'zendersoilerod'
   !-----------------------------------------------------------------------

      call control%ReadNML( bounds, NLFileName )
      call this%InitAllocate( bounds )

      if ( this%useStreams() )then     ! is this a namelist input and is it set in namelist default

         allocate(stream_varnames(1))
         stream_varnames = (/"mbl_bsn_fct_geo"/)  ! varname in the dust source file; the variable is dimensionless

         if (masterproc) then
            write(iulog,*) '  stream_varnames                  = ',stream_varnames
            flush(iulog)
         end if

         ! Initialize the cdeps data type sdat_erod
         call shr_strdata_init_from_inline(sdat_erod,                                  &     ! what is this function and where does it come from?
              my_task             = iam,                                                &
              logunit             = iulog,                                              &
              compname            = 'LND',                                              &
              model_clock         = model_clock,                                        &
              model_mesh          = mesh,                                               &
              stream_meshfile     = control%stream_meshfile_zendersoilerod,              &
              stream_lev_dimname  = 'null',                                             &
              stream_mapalgo      = control%zendersoilerod_mapalgo,                       &
              stream_filenames    = (/trim(control%stream_fldFileName_zendersoilerod)/), &
              stream_fldlistFile  = stream_varnames,                                    &
              stream_fldListModel = stream_varnames,                                    &
              stream_yearFirst    = 2003,                                               &
              stream_yearLast     = 2003,                                               &
              stream_yearAlign    = 1,                                                  &
              stream_offset       = 0,                                                  &
              stream_taxmode      = 'extend',                                           &
              stream_dtlimit      = 1.0e30_r8,                                          &
              stream_tintalgo     = 'linear',                                           &
              stream_name         = 'Zender soil erodibility',                          &
              rc                  = rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            write(iulog,*) 'Error on stream initialize -- see PET*.ESMF_LogFile(s)'
            call endrun("ESMF log error")
         end if

         ! Explicitly set current date to a hardcoded constant value. Otherwise
         ! using the real date can cause roundoff differences that are
         ! detrected as issues with exact restart.  EBK M05/20/2017
         ! call get_curr_date(year, mon, day, sec)
         year = 2003
         mon  = 12
         day  = 31
         sec  = 0
         mcdate = year*10000 + mon*100 + day

         call shr_strdata_advance(sdat_erod, ymd=mcdate, tod=sec, logunit=iulog, istr='zendersoilerod', rc=rc)   ! what is istr and do I need to change elsewhere because the change of istr here
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            write(iulog,*) 'Error on stream advance -- see PET*.ESMF_LogFile(s)'
            call endrun("ESMF log error")
         end if

         ! Get pointer for stream data that is time and spatially interpolate to model time and grid
         do n = 1,size(stream_varnames)
            call dshr_fldbun_getFldPtr(sdat_erod%pstrm(1)%fldbun_model, stream_varnames(n), fldptr1=dataptr1d, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
               write(iulog,*) 'Error on get field pointer -- see PET*.ESMF_LogFile(s)'
               call endrun("ESMF log error")
            end if
            if (trim(stream_varnames(n)) == 'mbl_bsn_fct_geo') then
               ig = 0
               do g = bounds%begg,bounds%endg
                  ig = ig+1
                  this%soil_erodibility(g) = dataptr1d(ig)
               end do

            end if

         end do
         ! TODO: EBK 03/25/2024: When shr_strdata adds a clean method we should invoke it here to save memory
         ! This is talked about in https://github.com/ESCOMP/CDEPS/issues/261

      end if

  end subroutine Init

  !==============================================================================
  logical function UseStreams(this)
    !
    ! !DESCRIPTION:
    ! Return true if the Zender method is being used and the soil erodability
    ! file is being used with it
    !
    ! !USES:
    use shr_dust_emis_mod, only : is_dust_emis_zender, is_zender_soil_erod_from_land
    !
    ! !ARGUMENTS:
    implicit none
    class(soil_erod_stream_type) :: this
    !
    ! !LOCAL VARIABLES:
    if ( .not. control%namelist_set )then
       call endrun(msg=' ERROR namelist NOT set before being used'//errMsg(sourcefile, __LINE__))
    end if
    if ( is_dust_emis_zender() .and. is_zender_soil_erod_from_land() )then
       UseStreams = .true.
    else
       UseStreams = .false.
    end if
  end function UseStreams

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
    class(soil_erod_stream_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    if ( this%useStreams() ) then
       allocate(this%soil_erodibility  (begg:endg))
    else
       allocate(this%soil_erodibility  (0))
    end if
    this%soil_erodibility  (:)   = nan

  end subroutine InitAllocate

  !==============================================================================
  subroutine CalcDustSource(this, bounds, soil_erod)
    !
    ! !DESCRIPTION:
    !  Calculate the soil eroditability for the Zender dust method.
    !
    ! !USES:
    use ColumnType              , only : col
    !use PatchType               , only : patch
    !USES
    use landunit_varcon         , only : istdlak
    use LandunitType            , only : lun
    !
    ! !ARGUMENTS:
    implicit none
    class(soil_erod_stream_type)             :: this
    type(bounds_type)              , intent(in)    :: bounds
    real(r8)                       , intent(inout) :: soil_erod(bounds%begc:)      ! [fraction] rock drag partition factor (roughness effect)
    !
    ! !LOCAL VARIABLES:
    !integer  :: g, c, fc    ! Indices
    integer  :: g, p, fp, l, c    ! Indices
    !real(r8) :: z0s         ! smooth roughness length (m)

    ! constants
    real(r8),parameter :: soil_erod_threshold = 0.1_r8   ! CAM soil erodibility threshold; below threshold -> soil_erod = 0_r8     11 Mar 2023
    !---------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(soil_erod)        == (/bounds%endc/)), sourcefile, __LINE__)

    !associate(                                                     &
         !z           =>   col%z                                   & ! Input:  [real(r8) (:,:) ]  layer depth (m) (-nlevsno+1:nlevsoi)
    !)


    ! dmleung: this loop truncates soil erodibility values smaller than a threshold value (set as 0.1). We save the drag partition factor as a grid level quantity.
    do c = bounds%begc,bounds%endc
        g = col%gridcell(c)
        l = col%landunit(c)
       if (lun%itype(l) /= istdlak) then ! not lake (can only be used during initialization)

          if (this%soil_erodibility(g) .lt. soil_erod_threshold ) then
             soil_erod(c) = 0._r8
          else
             soil_erod(c) = this%soil_erodibility(g)
          end if

       end if
    end do

    !end associate

  end subroutine CalcDustSource

  !==============================================================================
  subroutine ReadNML(this, bounds, NLFilename)
   !
   ! Read the namelist data stream information for the Zender method soil
   ! eroditability file
   !
   ! Uses:
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   use shr_dust_emis_mod, only : is_zender_soil_erod_from_land
   !
   ! arguments
   implicit none
   class(streamcontrol_type)     :: this
   type(bounds_type), intent(in) :: bounds
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: i         ! Indices
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   character(len=FL)  :: stream_fldFileName_zendersoilerod = ' '
   character(len=FL)  :: stream_meshfile_zendersoilerod = ' '
   character(len=CL)  :: zendersoilerod_mapalgo = ' '
   character(len=FL)  :: tmp_file_array(3)
   character(len=*), parameter :: namelist_name = 'zendersoilerod'    ! MUST agree with group name in namelist definition to read.
   character(len=*), parameter :: subName = "('zendersoilerod::ReadNML')"
   !-----------------------------------------------------------------------

   namelist /zendersoilerod/ &               ! MUST agree with namelist_name above
        zendersoilerod_mapalgo,  stream_fldFileName_zendersoilerod, &
        stream_meshfile_zendersoilerod

   ! Default values for namelist

   ! Read zenderdustsource namelist
   if (masterproc) then
      open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, namelist_name, status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=zendersoilerod, iostat=nml_error)   ! MUST agree with namelist_name above
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
   endif

   call shr_mpi_bcast(zendersoilerod_mapalgo            , mpicom)
   call shr_mpi_bcast(stream_fldFileName_zendersoilerod , mpicom)
   call shr_mpi_bcast(stream_meshfile_zendersoilerod    , mpicom)

   if (masterproc .and. is_zender_soil_erod_from_land() ) then
      write(iulog,*) ' '
      write(iulog,*) namelist_name, ' stream settings:'
      write(iulog,*) '  stream_fldFileName_zendersoilerod = ',stream_fldFileName_zendersoilerod
      write(iulog,*) '  stream_meshfile_zendersoilerod    = ',stream_meshfile_zendersoilerod
      write(iulog,*) '  zendersoilerod_mapalgo            = ',zendersoilerod_mapalgo
   endif

   tmp_file_array(1) = stream_fldFileName_zendersoilerod
   tmp_file_array(2) = stream_meshfile_zendersoilerod
   tmp_file_array(3) = zendersoilerod_mapalgo
   if ( is_zender_soil_erod_from_land() ) then
      do i = 1, size(tmp_file_array)
         if ( len_trim(tmp_file_array(i)) == 0 )then
            call endrun(msg=' ERROR '//trim(tmp_file_array(i))//' must be set when Zender_2003 is being used and zender_soil_erod_source is lnd'//errMsg(sourcefile, __LINE__))
         end if
      end do
   else
      do i = 1, size(tmp_file_array)
         if ( len_trim(tmp_file_array(i)) > 0 )then
            call endrun(msg=' ERROR '//trim(tmp_file_array(i))//' is set and MUST iNOT be when Zender_2003 is NOT being used or zender_soil_erod_source is atm'//errMsg(sourcefile, __LINE__))
         end if
      end do
   end if
   this%stream_fldFileName_zendersoilerod = stream_fldFileName_zendersoilerod
   this%stream_meshfile_zendersoilerod    = stream_meshfile_zendersoilerod
   this%zendersoilerod_mapalgo            = zendersoilerod_mapalgo

   this%namelist_set = .true.

 end subroutine ReadNML

end module ZenderSoilErodStreamType
