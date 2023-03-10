module IrrigationStreamMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read in irrigation from data stream
  !
  ! !USES:
  use shr_strdata_mod , only : shr_strdata_type, shr_strdata_create
  use shr_strdata_mod , only : shr_strdata_print, shr_strdata_advance
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_kind_mod    , only : CL => shr_kind_CL
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : scmlat,scmlon,single_column
  use clm_varctl      , only : iulog, use_irrigation_streams
  use clm_varcon      , only : grlnd
  use controlMod      , only : NLFilename
  use decompMod       , only : gsMap_lnd2Dsoi_gdc2glo
  use domainMod       , only : ldomain
  use fileutils       , only : getavu, relavu
  use GridcellType    , only : grc
  use LandunitType    , only : lun
  use PatchType       , only : patch
  use IrrigationMod   , only : irrigation_type
  use perf_mod        , only : t_startf, t_stopf
  use spmdMod         , only : masterproc
  use spmdMod         , only : mpicom, comp_id
  use mct_mod
  use ncdio_pio   
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
  type(shr_strdata_type) :: sdat_irrig    ! irrigation input data stream
  integer :: ism                          ! irrigation stream index
  integer :: nfields                      ! number of fields in irrigation stream
  integer, parameter :: nirrig = 2        ! number of crops on streams file
  character(SHR_KIND_CXX), allocatable    :: fldNames(:)
  integer, allocatable :: g_to_ig(:)      ! Array matching gridcell index to data index
  logical :: irrig_ignore_data_if_missing ! If should ignore overridding a point with irrigation data
                                          ! from the streams file, if the streams file shows that point 
                                          ! as missing (namelist item)
  !
  ! !PRIVATE TYPES:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains
  
  !-----------------------------------------------------------------------
  !
  ! irrigation_init
  !
  !-----------------------------------------------------------------------
  subroutine PrescribedIrrigationInit(bounds)
    !
    ! Initialize data stream information for irrigation.
    !
    !
    ! !USES:
    use clm_varctl       , only : inst_name
    use clm_time_manager , only : get_calendar
    use ncdio_pio        , only : pio_subsystem
    use shr_pio_mod      , only : shr_pio_getiotype
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use ndepStreamMod    , only : clm_domain_mct
    use shr_stream_mod   , only : shr_stream_file_null
    use shr_string_mod   , only : shr_string_listCreateField
    use clm_varpar       , only : nlevsoi
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds          ! bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: i                          ! index
    integer            :: stream_year_first_irrig    ! first year in Ustar stream to use
    integer            :: stream_year_last_irrig     ! last year in Ustar stream to use
    integer            :: model_year_align_irrig     ! align stream_year_first_irrig with 
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    integer            :: irrig_offset               ! Offset in time for dataset (sec)
    type(mct_ggrid)    :: dom_clm                    ! domain information 
    character(len=CL)  :: stream_fldfilename_irrig   ! ustar stream filename to read
    character(len=CL)  :: irrig_tintalgo = 'nearest' ! Time interpolation alogrithm

    character(*), parameter    :: subName = "('PrescribedIrrigationInit')"
    character(*), parameter    :: F00 = "('(PrescribedIrrigationInit) ',4a)"
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
    model_year_align_irrig       = 1      ! align stream_year_first_irrig with this model year
    stream_fldfilename_irrig     = shr_stream_file_null
    irrig_offset                 = 0
    irrig_ignore_data_if_missing = .false.

    ! Read irrig_streams namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
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
       call relavu( nu_nml )
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
       write(iulog,*) '  irrig_tintalgo = ',trim(irrig_tintalgo)
       write(iulog,*) '  irrig_offset = ',irrig_offset
       if ( irrig_ignore_data_if_missing )then
          write(iulog,*) '  Do NOT override a point with streams data if the streams data is missing'
       else
          write(iulog,*) '  Abort, if you find a model point where the input streams data is set to missing value'
       end if

    endif

    !scs
    call clm_domain_mct (bounds, dom_clm, nlevels=nirrig)
    !call clm_domain_mct (bounds, dom_clm)

    !
    ! create the field list for these fields...use in shr_strdata_create
    !
    fldList = 'irrig_rate:irrig_duration:irrig_start_time:crop_type:irrig_longitude:irrig_latitude'
    nfields = 6

    allocate(fldNames(nfields))
    fldNames = (/'irrig_rate','irrig_duration','irrig_start_time','crop_type','irrig_longitude','irrig_latitude'/)

    if (masterproc) write(iulog,*) 'fieldlist: ', trim(fldList)
    
    call shr_strdata_create(sdat_irrig,name="irrigation",    &
         pio_subsystem=pio_subsystem,                  & 
         pio_iotype=shr_pio_getiotype(inst_name),      &
         mpicom=mpicom, compid=comp_id,                &
         gsmap=gsMap_lnd2Dsoi_gdc2glo, ggrid=dom_clm,  &
         nxg=ldomain%ni, nyg=ldomain%nj,               &
         nzg=nirrig,                                  &
         yearFirst=stream_year_first_irrig,            &
         yearLast=stream_year_last_irrig,              &
         yearAlign=model_year_align_irrig,             &
         offset=irrig_offset,                          &
         domFilePath='',                               &
         domFileName=trim(stream_fldFileName_irrig),   &
         domTvarName='time',                           &
         domXvarName='lon' ,                           &
         domYvarName='lat' ,                           &  
         domZvarName='numcft' ,                        &  
         domAreaName='area',                           &
         domMaskName='mask',                           &
         filePath='',                                  &
         filename=(/stream_fldFileName_irrig/),        &
         fldListFile=fldList,                          &
         fldListModel=fldList,                         &
         fillalgo='none',                              &
         mapalgo='none',                               &
         tintalgo=irrig_tintalgo,                      &
         calendar=get_calendar(),                      &
!         dtlimit = 15._r8,                             &
         taxmode='cycle'                               )

    if (masterproc) then
       call shr_strdata_print(sdat_irrig,'irrigation data')
    endif

  end subroutine PrescribedIrrigationInit


  !-----------------------------------------------------------------------
  !
  ! PrescribedIrrigationAdvance
  !
  !-----------------------------------------------------------------------
  subroutine PrescribedIrrigationAdvance( bounds )
    !
    ! Advanace the prescribed irrigation stream
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    character(len=CL)  :: stream_var_name
    integer :: g, ig
    integer :: ier    ! error code
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)

    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day

!    stream_var_name = 'H2OSOI'

    ! Determine variable index
!    ism = mct_aVect_indexRA(sdat_irrig%avs(1),trim(stream_var_name))

    call shr_strdata_advance(sdat_irrig, mcdate, sec, mpicom, 'irrigation')

!!$    do g = 1, nfields
!!$       call shr_strdata_advance(sdat_irrig, mcdate, sec, mpicom, fldNames(g))
!!$    enddo

    ! Map gridcell to AV index
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

  !-----------------------------------------------------------------------
  !
  ! PrescribedIrrigationInterp
  !
  !-----------------------------------------------------------------------
  subroutine PrescribedIrrigationInterp(bounds, irrigation_inst)
    !
    ! Assign data stream information for prescribed irrigation.
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date, get_curr_time
    use clm_varpar      , only : numcft
    use clm_varcon      , only : denh2o, denice, watmin, spval
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)         , intent(in)    :: bounds
    type(irrigation_type)     , intent(inout) :: irrigation_inst
    !
    ! !LOCAL VARIABLES:
    integer :: p, g, j, ig, ip, n, ivt
    integer :: mcsec                      ! seconds of current date
    integer :: mdcur                      ! current day
    integer :: mscur                      ! seconds of current day
    integer :: mcdate                     ! current date
    integer :: yr,mon,day,nbsec           ! year,month,day,seconds components of a date
    integer :: hours,minutes,secs         ! hours,minutes,seconds of hh:mm:ss
    real(r8), allocatable :: irrig_rate_prescribed (:,:)    ! prescribed rate of irrigation
    real(r8), allocatable :: irrig_rate_duration   (:,:)    ! prescribed duration of irrigation
    real(r8), allocatable :: irrig_start_time   (:,:)
    real(r8), allocatable :: irrig_crop_type   (:,:)
    real(r8), allocatable :: irrig_lon   (:,:)
    real(r8), allocatable :: irrig_lat   (:,:)
    real(r8), parameter :: eps = 1e-3_r8 

    character(*), parameter    :: subName = "('PrescribedIrrigationInterp')"

    !-----------------------------------------------------------------------

!    SHR_ASSERT_FL( (lbound(sdat_irrig%avs(1)%rAttr,1) == ism ), sourcefile, __LINE__)
!    SHR_ASSERT_FL( (ubound(sdat_irrig%avs(1)%rAttr,1) == ism ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)

    associate( &
         irrig_rate_patch  =>    irrigation_inst%irrig_rate_patch       & ! Input:  [real(r8) (:,:) ]  current irrigation rate [mm/s]
         )
      SHR_ASSERT_FL( (lbound(irrig_rate_patch,1) <= bounds%begp ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(irrig_rate_patch,1) >= bounds%endp ), sourcefile, __LINE__)

      !
      ! Set the prescribed irrigation read from the file everywhere
      !

      call get_curr_time (mdcur, mscur)
      call get_curr_date (yr, mon, day, mcsec)

      allocate(irrig_rate_prescribed(bounds%begg:bounds%endg,nirrig))
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

             ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_rate')
             irrig_rate_prescribed(g,j) = sdat_irrig%avs(1)%rAttr(ip,n)

             ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_duration')
             irrig_rate_duration(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

             ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_start_time')
             irrig_start_time(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

             ip = mct_aVect_indexRA(sdat_irrig%avs(1),'crop_type')
             irrig_crop_type(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

             ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_longitude')
             irrig_lon(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

             ip = mct_aVect_indexRA(sdat_irrig%avs(1),'irrig_latitude')
             irrig_lat(g,j)   = sdat_irrig%avs(1)%rAttr(ip,n)

         end do
      end do
      

      do p = bounds%begp, bounds%endp

         g = patch%gridcell(p)
         ivt = patch%itype(p)

         do j = 1, nirrig

            if((abs(irrig_lon(g,j) - grc%londeg(g)) < eps) .and. (abs(irrig_lat(g,j) - grc%latdeg(g)) < eps) .and. (ivt == irrig_crop_type(g,j))) then


               if ((mscur >=  irrig_start_time(g,j)) .and. (mscur <=  (irrig_start_time(g,j)+irrig_rate_duration(g,j)))) then

                  irrig_rate_patch(p) = irrig_rate_prescribed(g,j)
               else
                  irrig_rate_patch(p) = 0._r8
               endif
            endif
         enddo
      enddo

      deallocate(irrig_rate_prescribed,irrig_rate_duration,irrig_start_time,irrig_crop_type,irrig_lon,irrig_lat)

    end associate

  end subroutine PrescribedIrrigationInterp

end module IrrigationStreamMod
