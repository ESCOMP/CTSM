module SoilMoistureStreamMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read in soil moisture from data stream
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
  use clm_varctl      , only : iulog, use_soil_moisture_streams
  use clm_varcon      , only : grlnd
  use controlMod      , only : NLFilename
  use decompMod       , only : gsMap_lnd2Dsoi_gdc2glo
  use domainMod       , only : ldomain
  use fileutils       , only : getavu, relavu
  use ColumnType      , only : col                
  use SoilStateType   , only : soilstate_type
  use WaterStateBulkType, only : waterstatebulk_type
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
  public :: PrescribedSoilMoistureInit    ! position datasets for soil moisture
  public :: PrescribedSoilMoistureInterp  ! interpolates between two periods of soil moisture data

  ! !PRIVATE MEMBER DATA:
  type(shr_strdata_type) :: sdat_soilm  ! soil moisture input data stream
  !
  ! !PRIVATE TYPES:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains
  
  !-----------------------------------------------------------------------
  !
  ! soil_moisture_init
  !
  !-----------------------------------------------------------------------
  subroutine PrescribedSoilMoistureInit(bounds)
    !
    ! Initialize data stream information for soil moisture.
    !
    !
    ! !USES:
    use clm_varctl       , only : inst_name
    use clm_time_manager , only : get_calendar
    use ncdio_pio        , only : pio_subsystem
    use shr_pio_mod      , only : shr_pio_getiotype
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use ndepStreamMod    , only : clm_domain_mct
    use histFileMod      , only : hist_addfld1d
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
    integer            :: stream_year_first_soilm    ! first year in Ustar stream to use
    integer            :: stream_year_last_soilm     ! last year in Ustar stream to use
    integer            :: model_year_align_soilm     ! align stream_year_first_soilm with 
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    type(mct_ggrid)    :: dom_clm                    ! domain information 
    character(len=CL)  :: stream_fldfilename_soilm   ! ustar stream filename to read
    character(len=CL)  :: soilm_mapalgo = 'bilinear' ! Mapping alogrithm

    character(*), parameter    :: subName = "('soil_moisture_init')"
    character(*), parameter    :: F00 = "('(soil_moisture_init) ',4a)"
    character(*), parameter    :: soilmString = "H2OSOI"  ! base string for field string
    character(SHR_KIND_CXX)    :: fldList            ! field string
    !-----------------------------------------------------------------------
    !
    ! deal with namelist variables here in init
    !
    namelist /soil_moisture_streams/         &
         stream_year_first_soilm,    &
         stream_year_last_soilm,     &
         model_year_align_soilm,     &
         soilm_mapalgo,              &
         stream_fldfilename_soilm

    ! Default values for namelist
    stream_year_first_soilm     = 1      ! first year in stream to use
    stream_year_last_soilm      = 1      ! last  year in stream to use
    model_year_align_soilm      = 1      ! align stream_year_first_soilm with this model year
    stream_fldfilename_soilm    = shr_stream_file_null

    ! Read soilm_streams namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
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
       call relavu( nu_nml )
    endif

    call shr_mpi_bcast(stream_year_first_soilm, mpicom)
    call shr_mpi_bcast(stream_year_last_soilm, mpicom)
    call shr_mpi_bcast(model_year_align_soilm, mpicom)
    call shr_mpi_bcast(stream_fldfilename_soilm, mpicom)

    if (masterproc) then

       write(iulog,*) ' '
       write(iulog,*) 'soil_moisture_stream settings:'
       write(iulog,*) '  stream_year_first_soilm  = ',stream_year_first_soilm  
       write(iulog,*) '  stream_year_last_soilm   = ',stream_year_last_soilm   
       write(iulog,*) '  model_year_align_soilm   = ',model_year_align_soilm   
       write(iulog,*) '  stream_fldfilename_soilm = ',trim(stream_fldfilename_soilm)

    endif

    call clm_domain_mct (bounds, dom_clm)

    !
    ! create the field list for these fields...use in shr_strdata_create
    !
    fldList = trim(soilmString)

    write(iulog,*) 'fieldlist: ', trim(fldList), nlevsoi
    
    call shr_strdata_create(sdat_soilm,name="soil_moisture",    &
         pio_subsystem=pio_subsystem,                  & 
         pio_iotype=shr_pio_getiotype(inst_name),      &
         mpicom=mpicom, compid=comp_id,                &
         gsmap=gsMap_lnd2Dsoi_gdc2glo, ggrid=dom_clm,  &
         nxg=ldomain%ni, nyg=ldomain%nj,               &
         nzg=nlevsoi,                                  &
         yearFirst=stream_year_first_soilm,            &
         yearLast=stream_year_last_soilm,              &
         yearAlign=model_year_align_soilm,             &
         offset=0,                                     &
         domFilePath='',                               &
         domFileName=trim(stream_fldFileName_soilm),   &
         domTvarName='time',                           &
         domXvarName='lon' ,                           &
         domYvarName='lat' ,                           &  
         domZvarName='levsoi' ,                        &  
         domAreaName='area',                           &
         domMaskName='mask',                           &
         filePath='',                                  &
         filename=(/stream_fldFileName_soilm/),        &
         fldListFile=fldList,                          &
         fldListModel=fldList,                         &
         fillalgo='none',                              &
         mapalgo=soilm_mapalgo,                        &
         calendar=get_calendar(),                      &
         dtlimit = 15._r8,                             &
         taxmode='cycle'                               )

    write(*,*) 'initializedstream'
    if (masterproc) then
       call shr_strdata_print(sdat_soilm,'soil moisture data')
    endif

  end subroutine PrescribedSoilMoistureInit

  !-----------------------------------------------------------------------
  !
  ! PrescribedSoilMoistureInterp
  !
  !-----------------------------------------------------------------------
  subroutine PrescribedSoilMoistureInterp(bounds, soilstate_inst, &
       waterstatebulk_inst)
    !
    ! Assign data stream information for prescribed soil moisture.
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date
    use clm_varpar      , only : nlevsoi
    use clm_varcon      , only : denh2o, denice 
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)         , intent(in)    :: bounds
    type(soilstate_type)      , intent(in)    :: soilstate_inst
    type(waterstatebulk_type) , intent(inout) :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c, g, j, ism, ig, n
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    integer :: ier    ! error code
    integer, allocatable :: g_to_ig(:)
    real(r8) :: soilm_liq_frac            ! liquid fraction of soil moisture
    real(r8) :: soilm_ice_frac            ! ice fraction of soil moisture
    real(r8) :: moisture_increment        ! soil moisture adjustment increment
    character(len=CL)  :: stream_var_name
    !-----------------------------------------------------------------------

    associate( &
         dz               =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         watsat           =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col    , & ! Input/Output:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col    , & ! Input/Output:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_vol       =>    waterstatebulk_inst%h2osoi_vol_col          & ! Output: volumetric soil water (m3/m3)
         )

      call get_curr_date(year, mon, day, sec)
      mcdate = year*10000 + mon*100 + day
      
      stream_var_name = 'H2OSOI'
      
      ! Determine variable index
      ism = mct_aVect_indexRA(sdat_soilm%avs(1),trim(stream_var_name))
      
      call shr_strdata_advance(sdat_soilm, mcdate, sec, mpicom, trim(stream_var_name))
      
      ! Map gridcell to AV index
      ier = 0
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
      
!      write(iulog,*) 'sdatavs: ',sdat_soilm%avs
      write(iulog,*) ' '
!      write(iulog,*) 'sdatsize: ',sdat_soilm%avs(1)
      call mct_aVect_info(2, sdat_soilm%avs(1))

      ! Read data from stream into column level variable
      
      do c = bounds%begc, bounds%endc
         !
         ! Set variable for each gridcell/column combination
         !
         ig = g_to_ig(col%gridcell(c))
         
         !       this is a 2d field (gridcell/nlevsoi) !
         do j = 1, nlevsoi
            ! initialize increment to original soil moisture value
            moisture_increment = h2osoi_vol(c,j)
            
            ! update volumetric soil moisture
            n = ig + (j-1)*ldomain%ni*ldomain%nj
            h2osoi_vol(c,j) = sdat_soilm%avs(1)%rAttr(ism,n)
            h2osoi_vol(c,j) = min(h2osoi_vol(c,j), watsat(c, j))
            
            ! calculate liq/ice mass fractions
            soilm_liq_frac  = h2osoi_liq(c, j) /(h2osoi_liq(c, j) + h2osoi_ice(c, j))
            soilm_ice_frac  = h2osoi_ice(c, j) /(h2osoi_liq(c, j) + h2osoi_ice(c, j))
            ! update moisture increment (relative to original value)
            moisture_increment = h2osoi_vol(c,j) - moisture_increment 
            
            ! update liq/ice water mass due to (volumetric) moisture increment
            h2osoi_liq(c,j) = h2osoi_liq(c,j) + (soilm_liq_frac * moisture_increment * dz(c, j) * denh2o)
            h2osoi_ice(c,j) = h2osoi_ice(c,j) + (soilm_ice_frac * moisture_increment * dz(c, j) * denice)
         enddo
         
      end do

    end associate

  end subroutine PrescribedSoilMoistureInterp

end module SoilMoistureStreamMod
