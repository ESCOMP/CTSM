module SoilMoistureStreamMod

  ! **********************************************************************
  ! --------------------------- IMPORTANT NOTE ---------------------------
  !
  ! In cases using the NUOPC driver/mediator, we use a different version of this module,
  ! based on CDEPS, which resides in src/cpl/nuopc/. Changes to the science here should
  ! also be made in the similar file in src/cpl/nuopc. Once we start using CDEPS by
  ! default, we can remove this version and move the CDEPS-based version into its place.
  ! **********************************************************************

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read in soil moisture from data stream
  !
  ! !USES:
  use shr_strdata_mod , only : shr_strdata_type, shr_strdata_create
  use shr_strdata_mod , only : shr_strdata_print, shr_strdata_advance
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_kind_mod    , only : CL => shr_kind_CL, CXX => shr_kind_CXX
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : iulog, use_soil_moisture_streams, inst_name
  use clm_varcon      , only : grlnd
  use controlMod      , only : NLFilename
  use decompMod       , only : gsMap_lnd2Dsoi_gdc2glo
  use domainMod       , only : ldomain
  use fileutils       , only : getavu, relavu
  use LandunitType    , only : lun
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
  public :: PrescribedSoilMoistureAdvance ! Advance the soil moisture stream (outside of Open-MP loops)
  public :: PrescribedSoilMoistureInterp  ! interpolates between two periods of soil moisture data

  ! !PRIVATE MEMBER DATA:
  type(shr_strdata_type) :: sdat_soilm    ! soil moisture input data stream
  integer :: ism                          ! Soil moisture steram index
  integer, allocatable :: g_to_ig(:)      ! Array matching gridcell index to data index
  logical :: soilm_ignore_data_if_missing ! If should ignore overridding a point with soil moisture data
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
  ! soil_moisture_init
  !
  !-----------------------------------------------------------------------
  subroutine PrescribedSoilMoistureInit(bounds)
    !
    ! Initialize data stream information for soil moisture.
    !
    !
    ! !USES:
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
    integer            :: stream_year_first_soilm    ! first year in Ustar stream to use
    integer            :: stream_year_last_soilm     ! last year in Ustar stream to use
    integer            :: model_year_align_soilm     ! align stream_year_first_soilm with
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    integer            :: soilm_offset               ! Offset in time for dataset (sec)
    type(mct_ggrid)    :: dom_clm                    ! domain information
    character(len=CL)  :: stream_fldfilename_soilm   ! ustar stream filename to read
    character(len=CL)  :: soilm_tintalgo = 'linear'  ! Time interpolation alogrithm

    character(*), parameter    :: subName = "('PrescribedSoilMoistureInit')"
    character(*), parameter    :: F00 = "('(PrescribedSoilMoistureInit) ',4a)"
    character(*), parameter    :: soilmString = "H2OSOI"  ! base string for field string
    character(CXX)             :: fldList            ! field string
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
       write(iulog,*) '  soilm_tintalgo = ',trim(soilm_tintalgo)
       write(iulog,*) '  soilm_offset = ',soilm_offset
       if ( soilm_ignore_data_if_missing )then
          write(iulog,*) '  Do NOT override a point with streams data if the streams data is missing'
       else
          write(iulog,*) '  Abort, if you find a model point where the input streams data is set to missing value'
       end if

    endif

    call clm_domain_mct (bounds, dom_clm, nlevels=nlevsoi)

    !
    ! create the field list for these fields...use in shr_strdata_create
    !
    fldList = trim(soilmString)

    if (masterproc) write(iulog,*) 'fieldlist: ', trim(fldList)

    call shr_strdata_create(sdat_soilm,name="soil_moisture",    &
         pio_subsystem=pio_subsystem,                  &
         pio_iotype=shr_pio_getiotype(inst_name),          &
         mpicom=mpicom, compid=comp_id,                &
         gsmap=gsMap_lnd2Dsoi_gdc2glo, ggrid=dom_clm,  &
         nxg=ldomain%ni, nyg=ldomain%nj,               &
         nzg=nlevsoi,                                  &
         yearFirst=stream_year_first_soilm,            &
         yearLast=stream_year_last_soilm,              &
         yearAlign=model_year_align_soilm,             &
         offset=soilm_offset,                          &
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
         mapalgo='none',                               &
         tintalgo=soilm_tintalgo,                      &
         calendar=get_calendar(),                      &
         dtlimit = 15._r8,                             &
         taxmode='cycle'                               )

    if (masterproc) then
       call shr_strdata_print(sdat_soilm,'soil moisture data')
    endif

  end subroutine PrescribedSoilMoistureInit


  !-----------------------------------------------------------------------
  !
  ! PrescribedSoilMoistureAdvance
  !
  !-----------------------------------------------------------------------
  subroutine PrescribedSoilMoistureAdvance( bounds )
    !
    ! Advanace the prescribed soil moisture stream
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

    stream_var_name = 'H2OSOI'

    ! Determine variable index
    ism = mct_aVect_indexRA(sdat_soilm%avs(1),trim(stream_var_name))

    call shr_strdata_advance(sdat_soilm, mcdate, sec, mpicom, trim(stream_var_name))

    ! Map gridcell to AV index
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
    use clm_varcon      , only : denh2o, denice, watmin, spval
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)         , intent(in)    :: bounds
    type(soilstate_type)      , intent(in)    :: soilstate_inst
    type(waterstatebulk_type) , intent(inout) :: waterstatebulk_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c, g, j, ig, n
    real(r8) :: soilm_liq_frac            ! liquid fraction of soil moisture
    real(r8) :: soilm_ice_frac            ! ice fraction of soil moisture
    real(r8) :: moisture_increment        ! soil moisture adjustment increment
    real(r8) :: h2osoi_vol_initial        ! initial vwc value
    character(*), parameter    :: subName = "('PrescribedSoilMoistureInterp')"

    !-----------------------------------------------------------------------

    SHR_ASSERT_FL( (lbound(sdat_soilm%avs(1)%rAttr,1) == ism ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(sdat_soilm%avs(1)%rAttr,1) == ism ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (lbound(g_to_ig,1) <= bounds%begg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(g_to_ig,1) >= bounds%endg ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (lbound(sdat_soilm%avs(1)%rAttr,2) <= g_to_ig(bounds%begg) ), sourcefile, __LINE__)
    SHR_ASSERT_FL( (ubound(sdat_soilm%avs(1)%rAttr,2) >= g_to_ig(bounds%endg)+(nlevsoi-1)*size(g_to_ig) ), sourcefile, __LINE__)
    associate( &
         dz               =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         watsat           =>    soilstate_inst%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col        , & ! Input/Output:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col        , & ! Input/Output:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_vol       =>    waterstatebulk_inst%h2osoi_vol_col        , & ! Output: volumetric soil water (m3/m3)
         h2osoi_vol_prs   =>    waterstatebulk_inst%h2osoi_vol_prs_grc      & ! Output: prescribed volumetric soil water (m3/m3)
         )
      SHR_ASSERT_FL( (lbound(h2osoi_vol,1) <= bounds%begc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_vol,1) >= bounds%endc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_vol,2) == 1 ),           sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_vol,2) >= nlevsoi ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(dz,1) <= bounds%begc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(dz,1) >= bounds%endc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(dz,2) <= 1 ),           sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(dz,2) >= nlevsoi ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(watsat,1) <= bounds%begc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(watsat,1) >= bounds%endc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(watsat,2) <= 1 ),           sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(watsat,2) >= nlevsoi ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_liq,1) <= bounds%begc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_liq,1) >= bounds%endc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_liq,2) <= 1 ),           sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_liq,2) >= nlevsoi ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_ice,1) <= bounds%begc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_ice,1) >= bounds%endc ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_ice,2) <= 1 ),           sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_ice,2) >= nlevsoi ),     sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_vol_prs,1) <= bounds%begg ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_vol_prs,1) >= bounds%endg ), sourcefile, __LINE__)
      SHR_ASSERT_FL( (lbound(h2osoi_vol_prs,2) == 1 ),           sourcefile, __LINE__)
      SHR_ASSERT_FL( (ubound(h2osoi_vol_prs,2) >= nlevsoi ),     sourcefile, __LINE__)
      !
      ! Set the prescribed soil moisture read from the file everywhere
      !
      do g = bounds%begg, bounds%endg
         ig = g_to_ig(g)
         do j = 1, nlevsoi

             !n = ig + (j-1)*size(g_to_ig)
             n = ig + (j-1)*size(g_to_ig)

             h2osoi_vol_prs(g,j) = sdat_soilm%avs(1)%rAttr(ism,n)

             ! If soil moiture is being interpolated in time and the result is
             ! large that probably means one of the two data points is missing (set to spval)
             if ( h2osoi_vol_prs(g,j) > 10.0_r8 .and. (h2osoi_vol_prs(g,j) /= spval) )then
                h2osoi_vol_prs(g,j) = spval
             end if

         end do
      end do

      do c = bounds%begc, bounds%endc
         !
         ! Set variable for each gridcell/column combination
         !
         g = col%gridcell(c)
         ig = g_to_ig(g)

         ! EBK Jan/2020, also check weights on gridcell (See https://github.com/ESCOMP/CTSM/issues/847)
         if ( (lun%itype(col%landunit(c)) == istsoil) .or. (lun%itype(col%landunit(c)) == istcrop) .and. &
              (col%wtgcell(c) /= 0._r8) ) then
            !       this is a 2d field (gridcell/nlevsoi) !
            do j = 1, nlevsoi

               n = ig + (j-1)*size(g_to_ig)

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
                        call endrun(subname // ' ERROR:: The input soil moisture stream is NOT vegetated for one of the land points' )
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
                  call endrun(subname // ':: ERROR h2osoil liquid plus ice is zero')
               endif
            enddo
         endif
      end do

    end associate

  end subroutine PrescribedSoilMoistureInterp

end module SoilMoistureStreamMod
