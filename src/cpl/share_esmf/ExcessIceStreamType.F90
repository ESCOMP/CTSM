module ExcessIceStreamType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in excess ice initial bulk values from data stream.
  ! Needed in parameterization for excess ice in soil (Lee et al., 2014).
  ! Used when use_excess_ice is true for initialization:
  ! startup type runs starting from coldstart and initial datasets
  ! that do not have required variables
  ! or hybrid runs from cases with use_excess_ice was false.
  ! Dataset is interpolated to 0.125x0.125 degrees grid from Brown et al., 1997
  ! with values derived from permafrost types.
  ! Values represent fraction of excess ice within soil column
  ! and are distributed within it later in initialization
  !
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

  public  :: UseExcessIceStreams      ! If streams will be used

  type, public :: excessicestream_type
     real(r8), pointer, private :: exice_bulk  (:)         ! excess ice bulk value (-)
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public  :: Init            ! Initialize and read data in
      procedure, public  :: CalcExcessIce   ! Calculate excess ice ammount

      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: InitAllocate   ! Allocate data

  end type excessicestream_type
    ! ! PRIVATE DATA:
  type, private :: streamcontrol_type
     character(len=FL)  :: stream_fldFileName_exice   ! data Filename
     character(len=FL)  :: stream_meshfile_exice      ! mesh Filename
     character(len=CL)  :: stream_mapalgo_exice       ! map algo
  contains
     procedure, private :: ReadNML     ! Read in namelist
  end type streamcontrol_type

  logical :: namelist_read = .false.
  type(streamcontrol_type), private :: control        ! Stream control data

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine Init(this, bounds, NLFilename)
    !
    use spmdMod          , only : iam
    use lnd_comp_shr     , only : mesh, model_clock
    use dshr_strdata_mod , only : shr_strdata_init_from_inline, shr_strdata_print
    use dshr_strdata_mod , only : shr_strdata_advance
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    !
    ! arguments
    implicit none
    class(excessicestream_type) :: this
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
    type(shr_strdata_type)         :: sdat_exice         ! input data stream
    character(len=16), allocatable :: stream_varnames(:) ! array of stream field names
    integer                        :: rc                 ! error code
    real(r8), pointer              :: dataptr1d(:)       ! temporary pointer
    character(len=*), parameter    :: stream_name = 'excess ice'

    call this%InitAllocate( bounds )
    call control%ReadNML( bounds, NLFileName )
    if ( UseExcessIceStreams() )then
      allocate(stream_varnames(1))
            stream_varnames = (/"EXICE"/)
      
      if (masterproc) then
        write(iulog,*) '  stream_varnames                  = ',stream_varnames
        write(iulog,*) '  Values will be used if the variable is not on the initial conditions dataset'
      end if

      call shr_strdata_init_from_inline(sdat_exice,                                    &
      my_task             = iam,                                                &
      logunit             = iulog,                                              &
      compname            = 'LND',                                              &
      model_clock         = model_clock,                                        &
      model_mesh          = mesh,                                               &
      stream_meshfile     = control%stream_meshfile_exice,                      &
      stream_lev_dimname  = 'null',                                             &
      stream_mapalgo      = control%stream_mapalgo_exice,                       &
      stream_filenames    = (/trim(control%stream_fldFileName_exice)/),         &
      stream_fldlistFile  = stream_varnames,                                    &
      stream_fldListModel = stream_varnames,                                    &
      stream_yearFirst    = 1996,                                               &
      stream_yearLast     = 1996,                                               &
      stream_yearAlign    = 1,                                                  &
      stream_offset       = 0,                                                  &
      stream_taxmode      = 'extend',                                           &
      stream_dtlimit      = 1.0e30_r8,                                          &
      ! in ch4FinundatedStreamType it is set to linear but we have a single date dataset
      stream_tintalgo     = 'nearest',                                          &
      stream_name         = 'excess ice ',                                      &
      rc                  = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
      
      !TODO
      ! Explicitly set current date to a hardcoded constant value. Otherwise
      ! using the real date can cause roundoff differences that are
      ! detrected as issues with exact restart.  EBK M05/20/2017
      ! call get_curr_date(year, mon, day, sec)
      year = 1996
      mon  = 12
      day  = 31
      sec  = 0
      mcdate = year*10000 + mon*100 + day
      
      call shr_strdata_advance(sdat_exice, ymd=mcdate, tod=sec, logunit=iulog, istr='exice', rc=rc) 
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
        
        ! Get pointer for stream data that is time and spatially interpolate to model time and grid
        do n = 1,size(stream_varnames)
          call dshr_fldbun_getFldPtr(sdat_exice%pstrm(1)%fldbun_model, stream_varnames(n), fldptr1=dataptr1d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
          if (trim(stream_varnames(n)) == 'EXICE') then
            ig = 0
            do g = bounds%begg,bounds%endg
              ig = ig+1
              this%exice_bulk(g) = dataptr1d(ig)
            end do
          end if
        end do
    end if
  end subroutine Init

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
    class(excessicestream_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begc, endc
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg


    allocate(this%exice_bulk(begg:endg))            ;  this%exice_bulk(:)   = nan

  end subroutine InitAllocate

  subroutine CalcExcessIce(this,bounds,exice_bulk_init)
  
  ! only transfers grid values to columns
   use shr_const_mod   , only : SHR_CONST_TKFRZ
   use landunit_varcon , only : istwet, istsoil, istcrop, istice
   use column_varcon   , only : icol_road_perv, icol_road_imperv
   use clm_varcon      , only : denice 
   use clm_varcon      , only : tfrz
   use ColumnType      , only : col
   use LandunitType    , only : lun
   implicit none
   class(excessicestream_type)        :: this
   type(bounds_type),  intent(in)     :: bounds
   real(r8)         ,  intent(inout)  :: exice_bulk_init(bounds%begc:bounds%endc) 
   !
   ! !LOCAL VARIABLES:
   integer  :: begc, endc
   integer  :: begg, endg
   integer  :: c, l, g  !counters
   
   exice_bulk_init(bounds%begc:bounds%endc)=0.0_r8

   do c = bounds%begc,bounds%endc
      g = col%gridcell(c)
      l = col%landunit(c)
      if ((.not. lun%lakpoi(l)) .and. (.not. lun%urbpoi(l)) .and. (.not. lun%itype(l) == istwet) .and. (.not. lun%itype(l) == istice)) then  !not lake
         if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
            exice_bulk_init(c)=this%exice_bulk(g)
         else 
            exice_bulk_init(c) = 0.0_r8
         endif
      else
         exice_bulk_init(c)=0.0_r8
      endif
   enddo

  end subroutine CalcExcessIce

  logical function UseExcessIceStreams()
  !
  ! !DESCRIPTION:
  ! Return true if
  !
  ! !USES:
  !
  ! !ARGUMENTS:
  implicit none
  !
  ! !LOCAL VARIABLES:
  if ( .not. namelist_read ) then
      call endrun(msg=' ERROR UseExcessIceStreams being called, but namelist has not been read yet'//errMsg(sourcefile, __LINE__))
  end if
  if ( trim(control%stream_fldFileName_exice) == '' )then
     UseExcessIceStreams = .false.
  else
     UseExcessIceStreams = .true.
  end if
end function UseExcessIceStreams

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
  logical            :: use_excess_ice_streams = .false.         ! logical to turn on use of excess ice streams
  character(len=FL)  :: stream_fldFileName_exice = ' '
  character(len=FL)  :: stream_meshfile_exice = ' '
  character(len=CL)  :: stream_mapalgo_exice = 'bilinear'
  character(len=*), parameter :: namelist_name = 'exice_streams'    ! MUST agree with name in namelist and read
  character(len=*), parameter :: subName = "('exice_streams::ReadNML')"
  !-----------------------------------------------------------------------

  namelist /exice_streams/ &               ! MUST agree with namelist_name above
       stream_mapalgo_exice,  stream_fldFileName_exice, stream_meshfile_exice, use_excess_ice_streams

  ! Default values for namelist

  ! Read excess ice namelist
  if (masterproc) then
     open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
     call shr_nl_find_group_name(nu_nml, namelist_name, status=nml_error)
     if (nml_error == 0) then
        read(nu_nml, nml=exice_streams,iostat=nml_error)   ! MUST agree with namelist_name above
        if (nml_error /= 0) then
           call endrun(msg=' ERROR reading '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
        end if
     else
        call endrun(msg=' ERROR finding '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
     end if
     close(nu_nml)
  endif

  call shr_mpi_bcast(use_excess_ice_streams   , mpicom)
  call shr_mpi_bcast(stream_mapalgo_exice     , mpicom)
  call shr_mpi_bcast(stream_fldFileName_exice , mpicom)
  call shr_mpi_bcast(stream_meshfile_exice    , mpicom)

  if (masterproc) then
     write(iulog,*) ' '
     if ( use_excess_ice_streams ) then
        write(iulog,*) 'excess ice streams are enabled: '
        write(iulog,*) namelist_name, ' stream settings:'
        write(iulog,*) '  stream_fldFileName_exice = ',stream_fldFileName_exice
        write(iulog,*) '  stream_meshfile_exice    = ',stream_meshfile_exice
        write(iulog,*) '  stream_mapalgo_exice     = ',stream_mapalgo_exice
        if ( trim(stream_fldFileName_exice) == '' )then
            call endrun(msg=' ERROR excess ice streams are on, but stream_fldFileName_exice is NOT set'//errMsg(sourcefile, __LINE__))
        end if
     else
        write(iulog,*) 'excess ice streams are off'
        if ( trim(stream_fldFileName_exice) /= '' )then
            call endrun(msg=' ERROR excess ice streams are off, but stream_fldFileName_exice is set'//errMsg(sourcefile, __LINE__))
        end if
     end if
  endif
  this%stream_fldFileName_exice = stream_fldFileName_exice
  this%stream_meshfile_exice    = stream_meshfile_exice
  this%stream_mapalgo_exice     = stream_mapalgo_exice
  namelist_read                 = .true.

end subroutine ReadNML


end module ExcessIceStreamType
