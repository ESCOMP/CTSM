module FanStreamMod

#include "shr_assert.h"

  !----------------------------------------------------------------------- 
  ! Contains methods for reading in FAN nitrogen deposition (in the form of
  ! manure) data file
  ! Also includes functions for fan stream file handling and 
  ! interpolation.
  !
  ! uses:
  use ESMF
  use dshr_strdata_mod , only : shr_strdata_type
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_cl
  use clm_varcon       , only : ispval
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use spmdMod          , only : mpicom, masterproc, comp_id, iam
  use clm_varctl       , only : iulog
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type

  !
  implicit none
  private

  ! Public interfaces
  public :: fanstream_init            ! Initialize FAN streams data
  public :: fanstream_interp          ! interpolates between two years of FAN file data
  public :: set_bcast_fanstream_pars  ! Set teh namelist parameters for the FAN streams

  ! Private module data
  type(shr_strdata_type)  :: sdat_fan            ! input data streams
  integer :: stream_year_first_fan = ispval      ! first year in stream to use
  integer :: stream_year_last_fan = ispval       ! last year in stream to use
  integer :: model_year_align_fan = ispval       ! align year to align model years with FAN streams
  character(len=CL)  :: stream_fldFileName_fan   ! FAN stream filename
  character(len=CL)  :: stream_meshFile_fan      ! FAN mesh filename
  character(len=CL)  :: fan_mapalgo = 'bilinear' ! FAN stream mapping algorithm
  integer, parameter :: nFields = 6
  character(len=16)  :: stream_varnames(nFields) ! array of stream field names
  logical :: crop_manure_per_crop                ! If manure is per crop or per land area, changes the variables read in
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================

  subroutine set_bcast_fanstream_pars(str_yr_first, str_yr_last, mdl_yr_align, mapalgo, str_filename, str_meshname, crop_man_is_percrop)
    !-----------------------------------------------------------------------
    !    
    ! Set the FAN stream namelist parameters
    !
    !-----------------------------------------------------------------------
    ! uses:
    use shr_mpi_mod, only : shr_mpi_bcast
    ! Arguments:
    implicit none
    integer, intent(in) :: str_yr_first, str_yr_last, mdl_yr_align
    ! whether manure_sgrz and manure_ngrz are per crop or land area:
    logical, intent(in) :: crop_man_is_percrop 
    character(len=*), intent(in) :: str_filename, str_meshname, mapalgo
    !-----------------------------------------------------------------------

    stream_year_first_fan = str_yr_first
    stream_year_last_fan = str_yr_last
    model_year_align_fan = mdl_yr_align
    stream_fldFileName_fan = str_filename
    stream_meshFile_fan = str_meshname
    crop_manure_per_crop = crop_man_is_percrop
    fan_mapalgo = mapalgo

    call shr_mpi_bcast(stream_year_first_fan, mpicom)
    call shr_mpi_bcast(stream_year_last_fan, mpicom)
    call shr_mpi_bcast(model_year_align_fan, mpicom)
    call shr_mpi_bcast(stream_fldFileName_fan, mpicom)
    call shr_mpi_bcast(stream_meshFile_fan, mpicom)
    call shr_mpi_bcast(crop_manure_per_crop, mpicom)
    call shr_mpi_bcast(fan_mapalgo, mpicom)
    
  end subroutine set_bcast_fanstream_pars

  !************************************************************************************
  
  subroutine fanstream_init(bounds, NLFilename)
   !-----------------------------------------------------------------------
   !    
   ! Initialize data stream information for FAN
   !
   !-----------------------------------------------------------------------
   ! uses:
   use clm_varctl                , only : inst_name
   use ncdio_pio                 , only : pio_subsystem
   use shr_pio_mod               , only : shr_pio_getiotype
   use shr_nl_mod                , only : shr_nl_find_group_name
   use lnd_comp_shr              , only : mesh, model_clock
   use dshr_strdata_mod          , only : shr_strdata_init_from_inline, shr_strdata_print
   use dshr_methods_mod          , only : dshr_fldbun_getfldptr
   !
   ! Arguments:
   implicit none
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! Local variables:
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   integer            :: rc        ! error code
   character(len=16)  :: streamvar, streamvar2    ! Specific stream variable names
   character(*), parameter :: subName = "('fanstream_init')"
   !-----------------------------------------------------------------------

   if (stream_year_first_fan == ispval) then
      call endrun(msg=subName//'ERROR stream_year_first_fan not set at '//errMsg(sourcefile, __LINE__))
   end if

   if (crop_manure_per_crop) then
      streamvar  = 'manure_sgrz_crop'
      streamvar2 = 'manure_ngrz_crop'
   else
      streamvar  = 'manure_sgrz'
      streamvar2 = 'manure_ngrz'
   end if
   stream_varnames = (/ &
        'manure_grz      ',   &
        streamvar,      &
        streamvar2,     &
        'fract_urea      ',   &
        'fract_nitr      ',   &
        'soilph          '    &
   /)
   
   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'FAN stream settings:'
      write(iulog,*) '  stream_year_first_fan  = ',stream_year_first_fan
      write(iulog,*) '  stream_year_last_fan   = ',stream_year_last_fan   
      write(iulog,*) '  model_year_align_fan   = ',model_year_align_fan   
      write(iulog,*) '  stream_fldFileName_fan = ',stream_fldFileName_fan
      write(iulog,*) '  stream_meshFile_fan    = ',stream_meshFile_fan
      write(iulog,*) '  stream_varnames        = ',stream_varnames
      write(iulog,*) ' '
   endif
   !
   ! Initialize the cdeps data type sdat_fan
   call shr_strdata_init_from_inline(sdat_fan, &
              my_task             = iam, &
              logunit             = iulog, &
              compname            = 'LND', &
              model_clock         = model_clock, &
              model_mesh          = mesh, &
              stream_meshfile     = stream_meshFile_fan, &
              stream_lev_dimname  = 'null', &
              stream_mapalgo      = fan_mapalgo, &
              stream_filenames    = (/trim(stream_fldFileName_fan)/), &
              stream_fldlistFile  = stream_varnames, &
              stream_fldListModel = stream_varnames, &
              stream_yearFirst    = stream_year_first_fan, &
              stream_yearLast     = stream_year_last_fan, &
              stream_yearAlign    = model_year_align_fan, &
              stream_offset       = 0, &
              stream_taxmode      = 'extend', &
              stream_dtlimit      = 1.0e30_r8, &
              stream_tintalgo     = 'linear', &
              stream_name         = 'FAN manure file', &
              rc                  = rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if

   if (masterproc) then
      call shr_strdata_print(sdat_fan,'CLMFAN data')
   endif

 end subroutine fanstream_init
  
 !================================================================

 subroutine fanstream_interp(bounds, atm2lnd_inst)
   !-----------------------------------------------------------------------
   !
   ! Interpoalte the FAN data to the current simulation time
   !
   !-----------------------------------------------------------------------
   use clm_time_manager, only : get_curr_date, get_curr_days_per_year
   use clm_varcon      , only : secspday
   use atm2lndType     , only : atm2lnd_type
   use shr_infnan_mod  , only : isinf => shr_infnan_isinf
   use dshr_methods_mod, only : dshr_fldbun_getfldptr
   use dshr_strdata_mod, only : shr_strdata_advance
   !
   ! Arguments
   implicit none
   type(bounds_type) , intent(in)    :: bounds  
   type(atm2lnd_type), intent(inout) :: atm2lnd_inst
   !
   ! Local variables
   integer :: g, ig , n ! Indices
   integer :: year    ! year (0, ...) for nstep+1
   integer :: mon     ! month (1, ..., 12) for nstep+1
   integer :: day     ! day of month (1, ..., 31) for nstep+1
   integer :: sec     ! seconds into current date for nstep+1
   integer :: mcdate  ! Current model date (yyyymmdd)
   integer :: dayspyr ! days per year
   integer :: rc      ! error code
   real(r8), pointer :: dataptr1d(:)  ! Temporary data array to put stream data into
   character(*), parameter :: subName = "('fanstream_interp')"
   !-----------------------------------------------------------------------
   SHR_ASSERT_FL( (lbound(atm2lnd_inst%forc_ndep_grz_grc,1) == bounds%begg ), sourcefile, __LINE__)
   SHR_ASSERT_FL( (ubound(atm2lnd_inst%forc_ndep_grz_grc,1) == bounds%endg ), sourcefile, __LINE__)

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day
   dayspyr = get_curr_days_per_year( )

   call shr_strdata_advance(sdat_fan, ymd=mcdate, tod=sec, logunit=iulog, istr='clmfan', rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if

   do n = 1, nFields
       dataptr1d => NULL()
       call dshr_fldbun_getFldPtr(sdat_fan%pstrm(1)%fldbun_model, trim(stream_varnames(n)), &
            fldptr1=dataptr1d, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       select case( trim(stream_varnames(n)) )
       case( 'manure_grz')
          do g = bounds%begg,bounds%endg
             atm2lnd_inst%forc_ndep_grz_grc(g) = dataptr1d(g-bounds%begg+1) / (secspday * dayspyr)
          end do
       case( 'manure_sgrz_crop')
          do g = bounds%begg,bounds%endg
             if ( isinf(dataptr1d(g)) ) then
                atm2lnd_inst%forc_ndep_sgrz_grc(g) = 9.0e99_r8
             else
                atm2lnd_inst%forc_ndep_sgrz_grc(g) = dataptr1d(g-bounds%begg+1) / (secspday * dayspyr)
             end if
          end do
       case( 'manure_ngrz_crop')
          do g = bounds%begg,bounds%endg
             if ( isinf(dataptr1d(g)) ) then
                atm2lnd_inst%forc_ndep_ngrz_grc(g) = 9.0e99_r8
             else
                atm2lnd_inst%forc_ndep_ngrz_grc(g) = dataptr1d(g-bounds%begg+1) / (secspday * dayspyr)
             end if
          end do
       case( 'manure_sgrz'     )
          do g = bounds%begg,bounds%endg
             if ( isinf(dataptr1d(g)) ) then
                atm2lnd_inst%forc_ndep_sgrz_grc(g) = 9.0e99_r8
             else
                atm2lnd_inst%forc_ndep_sgrz_grc(g) = dataptr1d(g-bounds%begg+1) / (secspday * dayspyr)
             end if
          end do
       case( 'manure_ngrz'     )
          do g = bounds%begg,bounds%endg
             if ( isinf(dataptr1d(g)) ) then
                atm2lnd_inst%forc_ndep_ngrz_grc(g) = 9.0e99_r8
             else
                atm2lnd_inst%forc_ndep_ngrz_grc(g) = dataptr1d(g-bounds%begg+1) / (secspday * dayspyr)
             end if
          end do
       case( 'fract_urea'      )
          do g = bounds%begg,bounds%endg
             atm2lnd_inst%forc_ndep_urea_grc(g) = dataptr1d(g-bounds%begg+1)
          end do
       case( 'fract_nitr'      )
          do g = bounds%begg,bounds%endg
             atm2lnd_inst%forc_ndep_nitr_grc(g) = dataptr1d(g-bounds%begg+1)
          end do
       case( 'soilph'          )
          do g = bounds%begg,bounds%endg
             atm2lnd_inst%forc_soilph_grc(g)= dataptr1d(g-bounds%begg+1)
          end do
       case default
          call endrun(msg=subName//'ERROR FAN stream variable is not handled'//trim(stream_varnames(n))//errMsg(sourcefile, __LINE__))
       end select
   end do
   dataptr1d => NULL()

 end subroutine fanstream_interp
    
end module FanStreamMod
