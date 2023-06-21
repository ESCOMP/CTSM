module IrrigationStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read in irrigation-related variables from data stream
  !
  ! !USES:
  use ESMF               , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_Finalize, ESMF_END_ABORT, ESMF_SUCCESS
  use dshr_strdata_mod   , only : shr_strdata_type, shr_strdata_print
  use dshr_strdata_mod   , only : shr_strdata_init_from_inline, shr_strdata_advance
  use dshr_methods_mod   , only : dshr_fldbun_getfldptr
  use dshr_stream_mod    , only : shr_stream_file_null
  use shr_kind_mod       , only : r8 => shr_kind_r8, cl => shr_kind_CL, cxx => shr_kind_CXX
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use shr_mpi_mod        , only : shr_mpi_bcast
  use decompMod          , only : bounds_type, subgrid_level_column
  use abortutils         , only : endrun
  use clm_varctl         , only : iulog, use_soil_moisture_streams
  use controlMod         , only : NLFilename
  use LandunitType       , only : lun
  use ColumnType         , only : col                
  use SoilStateType      , only : soilstate_type
  use WaterStateBulkType , only : waterstatebulk_type
  use perf_mod           , only : t_startf, t_stopf
  use spmdMod            , only : masterproc, mpicom, iam
  use clm_time_manager   , only : get_calendar, get_curr_date
  use clm_nlUtilsMod     , only : find_nlgroup_name
  use clm_varpar         , only : nlevsoi
  use clm_varcon         , only : denh2o, denice, watmin, spval
  use landunit_varcon    , only : istsoil, istcrop
  use lnd_comp_shr       , only : mesh, model_meshfile, model_clock
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: PrescribedIrrigationInit    ! position datasets for irrigation
  public :: PrescribedIrrigationAdvance ! Advance the irrigation stream (outside of Open-MP loops)
  public :: PrescribedIrrigationInterp  ! interpolates between two periods of irrigation data

  ! !PRIVATE MEMBER DATA:
  type(shr_strdata_type)  :: sdat_irrig                    ! soil moisture input data stream
  
  ! ADDED
  integer :: ism                          ! irrigation stream index
  integer :: nfields                      ! number of fields in irrigation stream
  integer, public, parameter :: nirrig = 2        ! number of fields in irrigation stream
  character(SHR_KIND_CXX), allocatable    :: fldNames(:)
  
  ! UNSURE
  character(*), parameter :: stream_var_name = "H2OSOI"    ! base string for field string
  character(len=CL)       :: stream_lev_dimname = 'levsoi' ! name of vertical layer dimension
  
  ! SURE
  integer, allocatable    :: g_to_ig(:)                    ! Array matching gridcell index to data index
  logical                 :: irrig_ignore_data_if_missing  ! If should ignore overridding a point with soil moisture data
                                                           ! from the streams file, if the streams file shows that point 
                                                           ! as missing (namelist item)
  ! !PRIVATE TYPES:
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  character(*), parameter :: u_FILE_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================
  
  subroutine PrescribedIrrigationInit(bounds)
    !
    ! Initialize data stream information for soil moisture.
    !
    ! !USES:
    ! !USES: UNSURE
    use clm_varctl       , only : inst_name
    use clm_time_manager , only : get_calendar
    use ncdio_pio        , only : pio_subsystem
    use shr_pio_mod      , only : shr_pio_getiotype
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use ndepStreamMod    , only : clm_domain_mct
    use shr_stream_mod   , only : shr_stream_file_null
    use shr_string_mod   , only : shr_string_listCreateField
    use clm_varpar       , only : nlevsoi
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)  :: bounds          ! bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: rc                         ! error coe
    integer            :: i                          ! index
    integer            :: stream_year_first_irrig    ! first year in Ustar stream to use
    integer            :: stream_year_last_irrig     ! last year in Ustar stream to use
    integer            :: model_year_align_irrig     ! align stream_year_first_soilm with 
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    integer            :: irrig_offset               ! Offset in time for dataset (sec)
    character(len=CL)  :: stream_fldfilename_irrig   ! ustar stream filename to read
    character(len=CL)  :: soilm_tintalgo = 'nearest'  ! Time interpolation alogrithm
    !character(len=CL)  :: stream_mapalgo = 'bilinear'
    !real(r8)           :: stream_dtlimit = 15._r8
    !character(len=CL)  :: stream_taxmode = 'cycle'
    character(*), parameter :: subName = "('PrescribedIrrigationInit')"
    character(*), parameter :: F00 = "('(PrescribedIrrigationInit) ',4a)"
	
	! ADDED
	character(SHR_KIND_CXX)    :: fldList            ! field string
    !-----------------------------------------------------------------------
    !
    ! deal with namelist variables here in init
    !
    namelist /irrigation_streams/   &
         stream_year_first_irrig,      &
         stream_year_last_irrig,       &
         model_year_align_irrig,       &
         irrig_tintalgo,               &
         irrig_offset,                 &
         irrig_ignore_data_if_missing, &
         stream_fldfilename_irrig

    ! Default values for namelist
    stream_year_first_irrig      = 1      ! first year in stream to use
    stream_year_last_irrig       = 1      ! last  year in stream to use
    model_year_align_irrig       = 1      ! align stream_year_first_soilm with this model year
    stream_fldfilename_irrig     = shr_stream_file_null
    irrig_offset                 = 0
    irrig_ignore_data_if_missing = .false.

    ! Read irrig_streams namelist
    if (masterproc) then
	   ! nu_nml = getavu()
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'irrigation_streams', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=irrigation_streams,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading irrigation_streams namelist')
          end if
       else
          call endrun(subname // ':: ERROR finding irrig_streams namelist')
       end if
       close(nu_nml)
	   ! call relavu( nu_nml )
    endif

    call shr_mpi_bcast(stream_year_first_irrig, mpicom)
    call shr_mpi_bcast(stream_year_last_irrig, mpicom)
    call shr_mpi_bcast(model_year_align_irrig, mpicom)
    call shr_mpi_bcast(stream_fldfilename_irrig, mpicom)
    call shr_mpi_bcast(irrig_tintalgo, mpicom)
    call shr_mpi_bcast(irrig_offset, mpicom)
    call shr_mpi_bcast(irrig_ignore_data_if_missing, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'irrigation_stream settings:'
       write(iulog,*) '  stream_year_first_irrig  = ',stream_year_first_irrig  
       write(iulog,*) '  stream_year_last_irrig   = ',stream_year_last_irrig   
       write(iulog,*) '  model_year_align_irrig   = ',model_year_align_irrig   
       write(iulog,*) '  stream_fldfilename_irrig = ',trim(stream_fldfilename_irrig)
       write(iulog,*) '  irrig_tintalgo           = ',trim(irrig_tintalgo)
       write(iulog,*) '  irrig_offset             = ',irrig_offset
       if ( irrig_ignore_data_if_missing ) then
          write(iulog,*) '  Do NOT override a point with streams data if the streams data is missing'
       else
          write(iulog,*) '  Abort, if you find a model point where the input streams data is set to missing value'
       end if
    endif

    ! Initialize the cdeps data type sdat_irrig
    ! TODO: for now stream_meshfile is the same as the model meshfile - must generalize this if want to have
    ! stream be at a different resolution

    call shr_strdata_init_from_inline(sdat_irrig,                  &
         my_task             = iam,                                &
         logunit             = iulog,                              &
         compname            = 'LND',                              &
         model_clock         = model_clock,                        &
         model_mesh          = mesh,                               &
         stream_meshfile     = model_meshfile,                     &
         stream_lev_dimname  = trim(stream_lev_dimname),           & 
         stream_mapalgo      = trim(stream_mapalgo),               &
         stream_filenames    = (/trim(stream_fldfilename_irrig)/), &
         stream_fldlistFile  = (/trim(stream_var_name)/),          &
         stream_fldListModel = (/trim(stream_var_name)/),          &
         stream_yearFirst    = stream_year_first_irrig,            &
         stream_yearLast     = stream_year_last_irrig,             &
         stream_yearAlign    = model_year_align_irrig,             &
         stream_offset       = irrig_offset,                       &
         stream_taxmode      = stream_taxmode,                     &
         stream_dtlimit      = stream_dtlimit,                     &
         stream_tintalgo     = irrig_tintalgo,                     &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    if (masterproc) then
       call shr_strdata_print(sdat_soilm,'irrigation data')
    endif

  end subroutine PrescribedIrrigationInit

  !==============================================================================
  subroutine PrescribedIrrigationAdvance( bounds )
    !
    ! Advanace the prescribed irrigation stream
    !
    ! !USES:
    ! use clm_time_manager, only : get_curr_date
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: rc
    integer :: g, ig
    integer :: ier     ! error code
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day

    ! Advance sdat stream
    call shr_strdata_advance(sdat_irrig, ymd=mcdate, tod=sec, logunit=iulog, istr='irrigation', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Map gridcell to 1->local_size (g_to_ig is a module variable)
    ier = 0
    if ( .not. allocated(g_to_ig) )then
       allocate (g_to_ig(bounds%begg:bounds%endg), stat=ier)
       if (ier /= 0) then
          write(iulog,*) 'Prescribed irrigation allocation error'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          g_to_ig(g) = ig
       end do
    end if

  end subroutine PrescribedIrrigationAdvance

  !==============================================================================
  subroutine PrescribedIrrigationInterp(bounds, irrigation_inst)
    !
    ! Assign data stream information for prescribed irrigation.
    !
    ! !USES:
	use clm_time_manager, only : get_curr_date, get_curr_time
    use clm_varpar      , only : numcft
	use clm_varcon      , only : denh2o, denice, watmin, spval
	use landunit_varcon , only : istsoil, istcrop
	use dshr_methods_mod , only : dshr_fldbun_getfldptr
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)         , intent(in)    :: bounds
    type(irrigation_type)      , intent(in)    :: irrigation_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: rc
    integer  :: c, g, j, ig, n
	
	integer :: p, g, j, ig, ip, n, ivt
	integer :: mcsec                      ! seconds of current date
    integer :: mdcur                      ! current day
    integer :: mscur                      ! seconds of current day
    integer :: mcdate                     ! current date
    integer :: yr,mon,day,nbsec           ! year,month,day,seconds components of a date
    integer :: hours,minutes,secs         ! hours,minutes,seconds of hh:mm:ss
	
	real(r8), allocatable :: sfc_irrig_rate_patch_prescribed (:,:)
	real(r8), allocatable :: irrig_rate_demand_patch_prescribed (:,:)
	
	real(r8), allocatable :: irrig_method (:,:)    ! prescribed type of irrigation
    real(r8), allocatable :: irrig_rate_duration   (:,:)    ! prescribed duration of irrigation
    real(r8), allocatable :: irrig_start_time  (:,:)
    integer,  allocatable :: irrig_crop_type   (:,:)
    real(r8), allocatable :: irrig_lon   (:,:)
    real(r8), allocatable :: irrig_lat   (:,:)
    real(r8), parameter :: eps = 1e-3_r8 
	
	
    !real(r8) :: soilm_liq_frac            ! liquid fraction of soil moisture
    !real(r8) :: soilm_ice_frac            ! ice fraction of soil moisture
    !real(r8) :: moisture_increment        ! soil moisture adjustment increment
    !real(r8) :: h2osoi_vol_initial        ! initial vwc value
    !real(r8), pointer :: dataptr2d(:,:)   ! first dimension is level, second is data on that level
    character(*), parameter    :: subName = "('PrescribedIrrigationInterp')"
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)
	
    associate( &
         sfc_irrig_rate_patch  =>    irrigation_inst%sfc_irrig_rate_patch       & 		! Input:  [real(r8) (:,:) ] patch irrigation flux (mm H2O/s) from surface water
		 irrig_rate_demand_patch  =>    irrigation_inst%irrig_rate_demand_patch       & ! Input:  [real(r8) (:,:) ] patch irrigation flux (mm H2O/s) demand
         )
      SHR_ASSERT_FL( (lbound(sfc_irrig_rate_patch,1) <= bounds%begp ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(sfc_irrig_rate_patch,1) >= bounds%endp ), sourcefile, __LINE__)
	  
	  SHR_ASSERT_FL( (lbound(irrig_rate_demand_patch,1) <= bounds%begp ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(irrig_rate_demand_patch,1) >= bounds%endp ), sourcefile, __LINE__)

      ! Get pointer for stream data that is time and spatially interpolate to model time and grid
      call dshr_fldbun_getFldPtr(sdat_soilm%pstrm(1)%fldbun_model, trim(stream_var_name), fldptr2=dataptr2d,  rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
	  

      ! Check that inner most dimensino of dataptr2d is equal to nlevsoi
      !if (size(dataptr2d,dim=1) /= nlevsoi) then
      !   if (masterproc) then
      !      write(iulog,*) 'ERROR: dataptr2d(dim=1) = ',size(dataptr2d,dim=1),&
      !           ' and  nlevsoi = ',nlevsoi,' must match '
      !   end if
      !   call endrun(trim(subname) // &
      !        ' ERROR:: The input soil moisture stream does not have levels equal to nlevsoi')
      !end if

      ! Set the prescribed soil moisture read from the file everywhere
	  call get_curr_time (mdcur, mscur)
      call get_curr_date (yr, mon, day, mcsec)
	  
	  allocate(sfc_irrig_rate_patch_prescribed(bounds%begg:bounds%endg,nirrig))
	  allocate(irrig_rate_demand_patch_prescribed(bounds%begg:bounds%endg,nirrig))
	  
	  allocate(irrig_method(bounds%begg:bounds%endg,nirrig))
      allocate(irrig_rate_duration(bounds%begg:bounds%endg,nirrig))
      allocate(irrig_start_time(bounds%begg:bounds%endg,nirrig))
      allocate(irrig_crop_type(bounds%begg:bounds%endg,nirrig))
      allocate(irrig_lon(bounds%begg:bounds%endg,nirrig))
      allocate(irrig_lat(bounds%begg:bounds%endg,nirrig))
	  
	  ! Read data from data streams
      do g = bounds%begg, bounds%endg
         ig = g_to_ig(g)

         do j = 1, nirrig
		 
		    !n = ig + (j-1)*size(g_to_ig)
            n = ig + (j-1)*size(g_to_ig)
			
			ip = mct_aVect_indexRA(sdat_irrig%avs(1),'sfc_irrig_rate_patch')
            sfc_irrig_rate_patch_prescribed(g,j) = sdat_irrig%avs(1)%rAttr(ip,n)
			
			ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_rate_demand_patch')
            irrig_rate_demand_patch_prescribed(g,j) = sdat_irrig%avs(1)%rAttr(ip,n)
			
            ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_method')
            irrig_method(g,j) = sdat_irrig%avs(1)%rAttr(ip,n)

            ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_duration')
            irrig_rate_duration(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

            ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_start_time')
            irrig_start_time(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

            ip = mct_aVect_indexRA(sdat_irrig%avs(1),'crop_type')
            irrig_crop_type(g,j)   = int(sdat_irrig%avs(1)%rAttr(ip,n))

            ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_longitude')
            irrig_lon(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

            ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_latitude')
            irrig_lat(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

         end do
      end do
	  
	  do p = bounds%begp, bounds%endp
         if (lun%itype(patch%landunit(p)) == istcrop) then
            g = patch%gridcell(p)
            ivt = patch%itype(p)


            do j = 1, nirrig

              ! if((abs(irrig_lon(g) - grc%londeg(g)) < eps) .and. (abs(irrig_lat(g) - grc%latdeg(g)) < eps) .and. (ivt == irrig_crop_type(g))) then
              ! if(ivt == int(irrig_crop_type(g,j))) then
			  if((abs(irrig_lon(g,j) - grc%londeg(g)) < eps) .and. (abs(irrig_lat(g,j) - grc%latdeg(g)) < eps) .and. (ivt == irrig_crop_type(g,j))) then

                  if ((mscur >=  irrig_start_time(g,j)) .and. (mscur <= (irrig_start_time(g,j)+irrig_rate_duration(g,j)))) then

                     sfc_irrig_rate_patch(p) = sfc_irrig_rate_patch_prescribed(g,j)
					 irrig_rate_demand_patch(p) = irrig_rate_demand_patch_prescribed(g,j)

                  else
                     sfc_irrig_rate_patch(p) = 0._r8
					 irrig_rate_demand_patch(p) = 0._r8
                  endif
               endif
            enddo
         endif
      enddo

      deallocate(sfc_irrig_rate_patch_prescribed, irrig_rate_demand_patch_prescribed, irrig_method,irrig_rate_duration,irrig_start_time,irrig_crop_type,irrig_lon,irrig_lat)
	  
	end associate
	
  end subroutine PrescribedIrrigationInterp

end module IrrigationStreamMod
