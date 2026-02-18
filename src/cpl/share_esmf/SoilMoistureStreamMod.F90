module SoilMoistureStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read in soil moisture from data stream
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
  use clm_varctl         , only : iulog, use_soil_moisture_streams, FL => fname_len
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
  public :: PrescribedSoilMoistureInit    ! position datasets for soil moisture
  public :: PrescribedSoilMoistureAdvance ! Advance the soil moisture stream (outside of Open-MP loops)
  public :: PrescribedSoilMoistureInterp  ! interpolates between two periods of soil moisture data

  ! !PRIVATE MEMBER DATA:
  type(shr_strdata_type)  :: sdat_soilm                    ! soil moisture input data stream
  character(*), parameter :: stream_var_name = "H2OSOI"    ! base string for field string
  character(len=CL)       :: stream_lev_dimname = 'levsoi' ! name of vertical layer dimension
  integer, allocatable    :: g_to_ig(:)                    ! Array matching gridcell index to data index
  logical                 :: soilm_ignore_data_if_missing  ! If should ignore overridding a point with soil moisture data
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
  
  subroutine PrescribedSoilMoistureInit(bounds)
    !
    ! Initialize data stream information for soil moisture.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)  :: bounds          ! bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: rc                         ! error coe
    integer            :: i                          ! index
    integer            :: stream_year_first_soilm    ! first year in Ustar stream to use
    integer            :: stream_year_last_soilm     ! last year in Ustar stream to use
    integer            :: model_year_align_soilm     ! align stream_year_first_soilm with 
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    integer            :: soilm_offset               ! Offset in time for dataset (sec)
    character(len=FL)  :: stream_fldfilename_soilm   ! ustar stream filename to read
    character(len=CL)  :: soilm_tintalgo = 'linear'  ! Time interpolation alogrithm
    character(len=CL)  :: stream_mapalgo = 'bilinear'
    real(r8)           :: stream_dtlimit = 15._r8
    character(len=CL)  :: stream_taxmode = 'cycle'
    character(*), parameter :: subName = "('PrescribedSoilMoistureInit')"
    character(*), parameter :: F00 = "('(PrescribedSoilMoistureInit) ',4a)"
    !-----------------------------------------------------------------------
    !
    ! deal with namelist variables here in init
    !
    namelist /soil_moisture_streams/   &
         stream_year_first_soilm,      &
         stream_year_last_soilm,       &
         model_year_align_soilm,       &
         soilm_tintalgo,               &
         soilm_offset,                 &
         soilm_ignore_data_if_missing, &
         stream_fldfilename_soilm

    ! Default values for namelist
    stream_year_first_soilm      = 1      ! first year in stream to use
    stream_year_last_soilm       = 1      ! last  year in stream to use
    model_year_align_soilm       = 1      ! align stream_year_first_soilm with this model year
    stream_fldfilename_soilm     = shr_stream_file_null
    soilm_offset                 = 0
    soilm_ignore_data_if_missing = .false.

    ! Read soilm_streams namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'soil_moisture_streams', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=soil_moisture_streams,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading soil_moisture_streams namelist')
          end if
       else
          call endrun(subname // ':: ERROR finding soilm_streams namelist')
       end if
       close(nu_nml)
    endif

    call shr_mpi_bcast(stream_year_first_soilm, mpicom)
    call shr_mpi_bcast(stream_year_last_soilm, mpicom)
    call shr_mpi_bcast(model_year_align_soilm, mpicom)
    call shr_mpi_bcast(stream_fldfilename_soilm, mpicom)
    call shr_mpi_bcast(soilm_tintalgo, mpicom)
    call shr_mpi_bcast(soilm_offset, mpicom)
    call shr_mpi_bcast(soilm_ignore_data_if_missing, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'soil_moisture_stream settings:'
       write(iulog,*) '  stream_year_first_soilm  = ',stream_year_first_soilm  
       write(iulog,*) '  stream_year_last_soilm   = ',stream_year_last_soilm   
       write(iulog,*) '  model_year_align_soilm   = ',model_year_align_soilm   
       write(iulog,*) '  stream_fldfilename_soilm = ',trim(stream_fldfilename_soilm)
       write(iulog,*) '  soilm_tintalgo           = ',trim(soilm_tintalgo)
       write(iulog,*) '  soilm_offset             = ',soilm_offset
       if ( soilm_ignore_data_if_missing ) then
          write(iulog,*) '  Do NOT override a point with streams data if the streams data is missing'
       else
          write(iulog,*) '  Abort, if you find a model point where the input streams data is set to missing value'
       end if
    endif

    ! Initialize the cdeps data type sdat_soilm
    ! TODO: for now stream_meshfile is the same as the model meshfile - must generalize this if want to have
    ! stream be at a different resolution

    call shr_strdata_init_from_inline(sdat_soilm,                  &
         my_task             = iam,                                &
         logunit             = iulog,                              &
         compname            = 'LND',                              &
         model_clock         = model_clock,                        &
         model_mesh          = mesh,                               &
         stream_meshfile     = model_meshfile,                     &
         stream_lev_dimname  = trim(stream_lev_dimname),           & 
         stream_mapalgo      = trim(stream_mapalgo),               &
         stream_filenames    = (/trim(stream_fldfilename_soilm)/), &
         stream_fldlistFile  = (/trim(stream_var_name)/),          &
         stream_fldListModel = (/trim(stream_var_name)/),          &
         stream_yearFirst    = stream_year_first_soilm,            &
         stream_yearLast     = stream_year_last_soilm,             &
         stream_yearAlign    = model_year_align_soilm,             &
         stream_offset       = soilm_offset,                       &
         stream_taxmode      = stream_taxmode,                     &
         stream_dtlimit      = stream_dtlimit,                     &
         stream_tintalgo     = soilm_tintalgo,                     &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    if (masterproc) then
       call shr_strdata_print(sdat_soilm,'soil moisture data')
    endif

  end subroutine PrescribedSoilMoistureInit

  !==============================================================================
  subroutine PrescribedSoilMoistureAdvance( bounds )
    !
    ! Advanace the prescribed soil moisture stream
    !
    ! !USES:
    !
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
    call shr_strdata_advance(sdat_soilm, ymd=mcdate, tod=sec, logunit=iulog, istr='soil_moisture', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Map gridcell to 1->local_size (g_to_ig is a module variable)
    ier = 0
    if ( .not. allocated(g_to_ig) )then
       allocate (g_to_ig(bounds%begg:bounds%endg), stat=ier)
       if (ier /= 0) then
          write(iulog,*) 'Prescribed soil moisture allocation error'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ig = 0
       do g = bounds%begg,bounds%endg
          ig = ig+1
          g_to_ig(g) = ig
       end do
    end if

  end subroutine PrescribedSoilMoistureAdvance

  !==============================================================================
  subroutine PrescribedSoilMoistureInterp(bounds, soilstate_inst, waterstatebulk_inst)
    !
    ! Assign data stream information for prescribed soil moisture.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)         , intent(in)    :: bounds
    type(soilstate_type)      , intent(in)    :: soilstate_inst
    type(waterstatebulk_type) , intent(inout) :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: rc
    integer  :: c, g, j, ig, n
    real(r8) :: soilm_liq_frac            ! liquid fraction of soil moisture
    real(r8) :: soilm_ice_frac            ! ice fraction of soil moisture
    real(r8) :: moisture_increment        ! soil moisture adjustment increment
    real(r8) :: h2osoi_vol_initial        ! initial vwc value
    real(r8), pointer :: dataptr2d(:,:)   ! first dimension is level, second is data on that level
    character(*), parameter    :: subName = "('PrescribedSoilMoistureInterp')"
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)
    associate( &
         dz               =>    col%dz                             ,   & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         watsat           =>    soilstate_inst%watsat_col          ,   & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col ,   & ! Input/Output:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col ,   & ! Input/Output:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_vol       =>    waterstatebulk_inst%h2osoi_vol_col ,   & ! Output: volumetric soil water (m3/m3)
         h2osoi_vol_prs   =>    waterstatebulk_inst%h2osoi_vol_prs_grc & ! Output: prescribed volumetric soil water (m3/m3)
         )
      SHR_ASSERT_FL( (lbound(h2osoi_vol,1) <= bounds%begc ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_vol,1) >= bounds%endc ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_vol,2) == 1 ),               sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_vol,2) >= nlevsoi ),         sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(dz,1) <= bounds%begc ),             sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(dz,1) >= bounds%endc ),             sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(dz,2) <= 1 ),                       sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(dz,2) >= nlevsoi ),                 sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(watsat,1) <= bounds%begc ),         sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(watsat,1) >= bounds%endc ),         sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(watsat,2) <= 1 ),                   sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(watsat,2) >= nlevsoi ),             sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_liq,1) <= bounds%begc ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_liq,1) >= bounds%endc ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_liq,2) <= 1 ),               sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_liq,2) >= nlevsoi ),         sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_ice,1) <= bounds%begc ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_ice,1) >= bounds%endc ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_ice,2) <= 1 ),               sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_ice,2) >= nlevsoi ),         sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_vol_prs,1) <= bounds%begg ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_vol_prs,1) >= bounds%endg ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_vol_prs,2) == 1 ),           sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_vol_prs,2) >= nlevsoi ),     sourcefile, __LINE__)

      ! Get pointer for stream data that is time and spatially interpolate to model time and grid
      call dshr_fldbun_getFldPtr(sdat_soilm%pstrm(1)%fldbun_model, trim(stream_var_name), fldptr2=dataptr2d,  rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      ! Check that inner most dimensino of dataptr2d is equal to nlevsoi
      if (size(dataptr2d,dim=1) /= nlevsoi) then
         if (masterproc) then
            write(iulog,*) 'ERROR: dataptr2d(dim=1) = ',size(dataptr2d,dim=1),&
                 ' and  nlevsoi = ',nlevsoi,' must match '
         end if
         call endrun(trim(subname) // &
              ' ERROR:: The input soil moisture stream does not have levels equal to nlevsoi')
      end if

      ! Set the prescribed soil moisture read from the file everywhere
      do g = bounds%begg, bounds%endg
         ig = g_to_ig(g)
         do j = 1, nlevsoi
             h2osoi_vol_prs(g,j) = dataptr2d(j,ig)

             ! If soil moiture is being interpolated in time and the result is
             ! large that probably means one of the two data points is missing (set to spval)
             if ( h2osoi_vol_prs(g,j) > 10.0_r8 .and. (h2osoi_vol_prs(g,j) /= spval) )then
                h2osoi_vol_prs(g,j) = spval
             end if
         end do
      end do
      
      do c = bounds%begc, bounds%endc
         ! Set variable for each gridcell/column combination
         g = col%gridcell(c)
         ig = g_to_ig(g)
            
         ! EBK Jan/2020, also check weights on gridcell (See https://github.com/ESCOMP/CTSM/issues/847)
         if ( (lun%itype(col%landunit(c)) == istsoil) .or. &
              (lun%itype(col%landunit(c)) == istcrop) .and. (col%wtgcell(c) /= 0._r8) ) then

            ! this is a 2d field (gridcell/nlevsoi) !
            do j = 1, nlevsoi
               ! if soil water is zero, liq/ice fractions cannot be calculated
               if((h2osoi_liq(c, j) + h2osoi_ice(c, j)) > 0._r8) then
                  
                  ! save original soil moisture value
                  h2osoi_vol_initial = h2osoi_vol(c,j)
            
                  ! Check if the vegetated land mask from the dataset on the
                  ! file is different
                  if ( (h2osoi_vol_prs(g,j) == spval) .and. (h2osoi_vol_initial /= spval) )then
                     if ( soilm_ignore_data_if_missing )then
                        cycle
                     else
                        write(iulog,*) 'Input soil moisture dataset is not vegetated as expected: gridcell=', &
                                        g, ' active = ', col%active(c)
                        call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                             msg = subname //&
                             ' ERROR:: The input soil moisture stream is NOT vegetated for one of the land points' )
                     end if
                  end if

                  ! update volumetric soil moisture from data prescribed from the file
                  h2osoi_vol(c,j) = h2osoi_vol_prs(g,j)

                  ! calculate liq/ice mass fractions
                  soilm_liq_frac  = h2osoi_liq(c, j) /(h2osoi_liq(c, j) + h2osoi_ice(c, j))
                  soilm_ice_frac  = h2osoi_ice(c, j) /(h2osoi_liq(c, j) + h2osoi_ice(c, j))

                  ! calculate moisture increment
                  moisture_increment = h2osoi_vol(c,j) - h2osoi_vol_initial
                  ! add limitation check
                  moisture_increment = min((watsat(c,j) - h2osoi_vol_initial),max(-(h2osoi_vol_initial-watmin),moisture_increment))

                  ! update liq/ice water mass due to (volumetric) moisture increment
                  h2osoi_liq(c,j) = h2osoi_liq(c,j) + (soilm_liq_frac * moisture_increment * dz(c, j) * denh2o)
                  h2osoi_ice(c,j) = h2osoi_ice(c,j) + (soilm_ice_frac * moisture_increment * dz(c, j) * denice)

               else
                  call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                       msg = subname // ':: ERROR h2osoil liquid plus ice is zero')
               endif
            enddo
         endif      
      end do

    end associate

  end subroutine PrescribedSoilMoistureInterp

end module SoilMoistureStreamMod
