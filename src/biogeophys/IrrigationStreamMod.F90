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
  use PatchType          , only : patch
  use IrrigationMod      , only : irrigation_type,irrig_method_drip
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
  type(shr_strdata_type)  :: sdat_irrig                 ! irrigation input data stream
  
  integer :: ism                                        ! irrigation stream index
  integer     , parameter :: numFields = 8              ! number of fields to build field string
  character(cl)           :: stream_varnames(numFields) ! field string
  character(len=CL)       :: stream_lev_dimname = 'cft' ! name of vertical layer dimension
  integer :: num_irrig_cft                              ! number of irrigated cfts on stream file

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
    character(len=CL)  :: irrig_tintalgo = 'nearest'  ! Time interpolation alogrithm
    character(len=CL)  :: stream_mapalgo = 'nn'
    real(r8)           :: stream_dtlimit = 15._r8
    character(len=CL)  :: stream_taxmode = 'cycle'
    character(*), parameter :: subName = "('PrescribedIrrigationInit')"
    character(*), parameter :: F00 = "('(PrescribedIrrigationInit) ',4a)"
	
    ! ADDED
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

    stream_varnames = (/'irrig_depth_column','irrig_target_smp_column','irrig_fraction_column','irrig_method_patch','irrig_duration_column','irrig_start_time_column','surface_water_ponding_column','crop_type'/)

    if (masterproc) then
       do i = 1,numFields
          write(iulog,'(a,a)' ) '  stream_varname         = ',trim(stream_varnames(i))
       end do
    endif

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
         stream_fldlistFile  = stream_varnames,          &
         stream_fldListModel = stream_varnames,          &
         stream_yearFirst    = stream_year_first_irrig,            &
         stream_yearLast     = stream_year_last_irrig,             &
         stream_yearAlign    = model_year_align_irrig,             &
         stream_offset       = irrig_offset,                       &
         stream_taxmode      = stream_taxmode,                     &
         stream_dtlimit      = stream_dtlimit,                     &
         stream_tintalgo     = irrig_tintalgo,                     &
         stream_name         = 'Irrig',                            &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if	

    ! get actual number of cfts on stream file
    num_irrig_cft = sdat_irrig%pstrm(1)%stream_nlev

    if (masterproc) then
       call shr_strdata_print(sdat_irrig,'irrigation data')
       write(iulog,*) 'number of irrigated cfts on stream file ',num_irrig_cft
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
	use clm_varcon      , only : denh2o, denice, watmin, spval
	use landunit_varcon , only : istsoil, istcrop
	use dshr_methods_mod , only : dshr_fldbun_getfldptr
	
	use GridcellType     , only : grc
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)         , intent(in)    :: bounds
    type(irrigation_type)      , intent(in)    :: irrigation_inst
    ! !LOCAL VARIABLES:
    integer :: rc
    integer :: lsize
    integer :: c	
    integer :: p, g, j, ig, ip, n, ivt
    integer :: mcsec                      ! seconds of current date
    integer :: mdcur                      ! current day
    integer :: mscur                      ! seconds of current day
    integer :: mcdate                     ! current date
    integer :: yr,mon,day,nbsec           ! year,month,day,seconds components of a date
    integer :: hours,minutes,secs         ! hours,minutes,seconds of hh:mm:ss
	
    real(r8), allocatable :: irrig_rate_prescribed  (:,:)
    integer, allocatable :: irrig_method_prescribed (:,:)    
	integer, allocatable :: irrig_method	 (:,:)												 
    integer, allocatable :: irrig_duration    (:,:)
    integer, allocatable :: irrig_start_time  (:,:)
    integer,  allocatable :: irrig_crop_type   (:,:)
    real(r8), allocatable :: irrig_lon     (:,:)
    real(r8), allocatable :: irrig_lat     (:,:)
    real(r8), allocatable :: irrig_target  (:,:)
    real(r8), allocatable :: irrig_depth   (:,:)
    real(r8), allocatable :: irrig_fraction(:,:)
	real(r8), allocatable :: surface_water_ponding(:,:)
    real(r8), parameter :: eps = 1e-3_r8 
	
	
    real(r8), pointer :: dataptr1d(:)
    real(r8), pointer :: dataptr2d(:,:)
	
	integer, pointer :: dataptr2d_int(:,:)
	
    character(*), parameter    :: subName = "('PrescribedIrrigationInterp')"
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)
	
    associate( &
         sfc_irrig_rate_patch     			 =>    irrigation_inst%sfc_irrig_rate_patch      		, & 		! Input:  [real(r8) (:,:) ] patch irrigation flux (mm H2O/s) from surface water
         irrig_rate_demand_patch  			 =>    irrigation_inst%irrig_rate_demand_patch   		, & ! Input:  [real(r8) (:,:) ] patch irrigation flux (mm H2O/s) demand
         irrig_method_patch 				 =>    irrigation_inst%irrig_method_patch             	, & ! Input:  [real(r8) (:,:) ] patch irrigation method (drip,sprinkler,flood)
         irrig_target_smp_column  			 =>    irrigation_inst%irrig_target_smp_column  		, & ! Input:  [real(r8) (:,:) ] column irrigation target soil matric potential (mm)
         irrig_depth_column       			 =>    irrigation_inst%irrig_depth_column       		, & ! Input:  [real(r8) (:,:) ] column irrigation depth (m)
         irrig_threshold_fraction_column     =>    irrigation_inst%irrig_threshold_fraction_column  ,& ! Input:  [real(r8) (:,:) ] column irrigation threshold fraction
		 irrig_start_time_column  			 =>    irrigation_inst%irrig_start_time_column  		, & ! Input:  [real(r8) (:,:) ] column irrigation target soil matric potential (mm)
         irrig_duration_column       		 =>    irrigation_inst%irrig_duration_column       		, & ! Input:  [real(r8) (:,:) ] column irrigation depth (m)
		 surface_water_ponding_column    	 =>    irrigation_inst%surface_water_ponding_column		  & ! Input:  [real(r8) (:,:) ] column irrigation threshold fraction
		 )
      SHR_ASSERT_FL( (lbound(sfc_irrig_rate_patch,1) <= bounds%begp ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(sfc_irrig_rate_patch,1) >= bounds%endp ), sourcefile, __LINE__)	
      SHR_ASSERT_FL( (lbound(sfc_irrig_rate_patch,2) == 1 ), sourcefile, __LINE__)
!      SHR_ASSERT_FL( (ubound(sfc_irrig_rate_patch,2) >= 64 ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(irrig_rate_demand_patch,1) >= bounds%endp ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(irrig_rate_demand_patch,1) <= bounds%begp ), sourcefile, __LINE__)	
      SHR_ASSERT_FL( (lbound(irrig_rate_demand_patch,2) == 1 ),                        sourcefile, __LINE__)
!      SHR_ASSERT_FL( (ubound(irrig_rate_demand_patch,2) >= 64 ),          sourcefile, __LINE__)
	  
      ! Get pointer for stream data that is time and spatially interpolate to model time and grid

      ! Place all data from each type into a temporary 2d array
      lsize = bounds%endg - bounds%begg + 1
	
	  ! Read variables based on stream_varnames
      do n = 1,numFields
         call dshr_fldbun_getFldPtr(sdat_irrig%pstrm(1)%fldbun_model, trim(stream_varnames(n)), &
              fldptr2=dataptr2d,  rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! Read data from data streams
         select case (stream_varnames(n))
         case ('sfc_irrig_rate_patch')
            allocate(irrig_rate_prescribed(bounds%begg:bounds%endg,num_irrig_cft))    
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_rate_prescribed(g,:) = dataptr2d(:,ig)
            enddo
         case ('irrig_method_patch')
            allocate(irrig_method_prescribed(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_method_prescribed(g,:) = dataptr2d(:,ig)
            enddo
         case ('irrig_duration_column')
            allocate(irrig_duration(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_duration(g,:) = INT(dataptr2d(:,ig))
            enddo
         case ('irrig_start_time_column')
            allocate(irrig_start_time(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_start_time(g,:) = INT(dataptr2d(:,ig))
            enddo
         case ('crop_type')
            allocate(irrig_crop_type(bounds%begg:bounds%endg,num_irrig_cft))
	        do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_crop_type(g,:) = dataptr2d(:,ig)
            enddo
         case ('irrig_longitude')
            allocate(irrig_lon(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_lon(g,:) = dataptr2d(:,ig)
            enddo
         case ('irrig_latitude')
            allocate(irrig_lat(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_lat(g,:) = dataptr2d(:,ig)
            enddo
         case ('irrig_target_smp_column')
            allocate(irrig_target(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_target(g,:) = dataptr2d(:,ig)
            enddo
         case ('irrig_depth_column')
            allocate(irrig_depth(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_depth(g,:) = dataptr2d(:,ig)
            enddo
         case ('irrig_fraction_column')
            allocate(irrig_fraction(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               irrig_fraction(g,:) = dataptr2d(:,ig) 
            enddo
		 case ('surface_water_ponding_column')
            allocate(surface_water_ponding(bounds%begg:bounds%endg,num_irrig_cft))
            do g = bounds%begg, bounds%endg
               ig = g_to_ig(g)
               surface_water_ponding(g,:) = dataptr2d(:,ig) 
            enddo
         end select

      end do
      call get_curr_time (mdcur, mscur)
      call get_curr_date (yr, mon, day, mcsec)
	  
	  
      do p = bounds%begp, bounds%endp
         if (lun%itype(patch%landunit(p)) == istcrop) then
            g = patch%gridcell(p)
            c = patch%column(p)
            ivt = patch%itype(p)

            !set all to zero as test                     
            !sfc_irrig_rate_patch(p)    = 0._r8
            !irrig_rate_demand_patch(p) = 0._r8

            do j = 1, num_irrig_cft
              if(ivt == int(irrig_crop_type(g,j))) then
                     ! check if these params on stream file before assign them to arrays in irrigationMod.F90
                     if (allocated(irrig_method_prescribed)) then
                        irrig_method_patch(p)      = irrig_method_prescribed(g,j)
                     endif
                     if (allocated(irrig_target)) then
                        irrig_target_smp_column(c) = irrig_target(g,j)
					 end if
					 if (allocated(irrig_depth)) then
                        irrig_depth_column(c)      = irrig_depth(g,j)
					 end if
                     if (allocated(irrig_fraction)) then
                        irrig_threshold_fraction_column(c)   = irrig_fraction(g,j)
					 end if
                     if (allocated(irrig_start_time)) then
                        irrig_start_time_column(c)      = INT(irrig_start_time(g,j))
					 end if
					 if (allocated(irrig_duration)) then
						irrig_duration_column(c)      = INT(irrig_duration(g,j))
                     endif
					 if (allocated(surface_water_ponding)) then 
						surface_water_ponding_column(c)      = surface_water_ponding(g,j)
                     endif
                     ! prescribed irrigation rate on stream file
                     if (allocated(irrig_rate_prescribed)) then 
                        if ((mscur >=  irrig_start_time(g,j)) .and. (mscur <= (irrig_start_time(g,j)+irrig_duration(g,j)))) then
                           
                           sfc_irrig_rate_patch(p) = irrig_rate_prescribed(g,j)
                           irrig_rate_demand_patch(p) = irrig_rate_prescribed(g,j)
                           
                        else
                           sfc_irrig_rate_patch(p)    = 0._r8
                           irrig_rate_demand_patch(p) = 0._r8
                        endif
                     endif

               !endif

               endif
            enddo
         endif
      enddo

      if (allocated(irrig_rate_prescribed)) then
         deallocate(irrig_rate_prescribed)
      endif
      if (allocated(irrig_method_prescribed)) then
         deallocate(irrig_method_prescribed)
      endif
      if (allocated(irrig_duration)) then
         deallocate(irrig_duration)
      endif
      if (allocated(irrig_start_time)) then
         deallocate(irrig_start_time)
      endif
      if (allocated(irrig_crop_type)) then
         deallocate(irrig_crop_type)
      endif
      if (allocated(irrig_lon)) then
         deallocate(irrig_lon)
      endif
      if (allocated(irrig_lat)) then
         deallocate(irrig_lat)
      endif
      if (allocated(irrig_target)) then
         deallocate(irrig_target)
      endif
      if (allocated(irrig_depth)) then
         deallocate(irrig_depth)
      endif
      if (allocated(irrig_fraction)) then
         deallocate(irrig_fraction)
      endif
	  if (allocated(irrig_method)) then
         deallocate(irrig_method)
      endif
	  if (allocated(surface_water_ponding)) then
         deallocate(surface_water_ponding)
      endif

    end associate
  end subroutine PrescribedIrrigationInterp

end module IrrigationStreamMod
