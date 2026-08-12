module WaterStateType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water state variables that apply to both bulk water
  ! and water tracers.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use decompMod      , only : subgrid_level_patch, subgrid_level_column, subgrid_level_landunit, subgrid_level_gridcell
  use clm_varctl     , only : use_bedrock, use_excess_ice, iulog
  use spmdMod        , only : masterproc
  use clm_varctl     , only : use_fates, use_hillslope
  use clm_varpar     , only : nlevgrnd, nlevsoi, nlevurb, nlevmaxurbgrnd, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  use WaterTracerUtils, only : AllocateVar1d, AllocateVar2d
  use ExcessIceStreamType, only : excessicestream_type, UseExcessIceStreams
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterstate_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: h2osno_no_layers_col   (:)   ! col snow that is not resolved into layers; this is non-zero only if there is too little snow for there to be explicit snow layers (mm H2O)
     real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_vol_col         (:,:) ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: h2osoi_vol_prs_grc     (:,:) ! grc volumetric soil water prescribed (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: h2osfc_col             (:)   ! col surface water (mm H2O)
     real(r8), pointer :: snocan_patch           (:)   ! patch canopy snow water (mm H2O)
     real(r8), pointer :: liqcan_patch           (:)   ! patch canopy liquid water (mm H2O)

     real(r8), pointer :: wa_col                 (:)   ! col water in the unconfined aquifer (mm)

     ! For the following dynbal baseline variables: positive values are subtracted to
     ! avoid counting liquid water content of "virtual" states; negative values are added
     ! to account for missing states in the model.
     real(r8), pointer :: dynbal_baseline_liq_col(:)    ! baseline liquid water content subtracted from each column's total liquid water calculation (mm H2O)
     real(r8), pointer :: dynbal_baseline_ice_col(:)    ! baseline ice content subtracted from each column's total ice calculation (mm H2O)

     real(r8) :: aquifer_water_baseline                 ! baseline value for water in the unconfined aquifer (wa_col) for this bulk / tracer (mm)

     real(r8), pointer :: excess_ice_col         (:,:)  ! col excess ice (kg/m2) (new) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: exice_bulk_init        (:)    ! inital value for excess ice (new) (unitless)

     ! Hillslope stream variables
     real(r8), pointer :: stream_water_volume_lun(:)   ! landunit volume of water in the streams (m3)

   contains

     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: CalculateTotalH2osno
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private :: CheckSnowConsistency

  end type waterstate_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info, tracer_vars, &
       h2osno_input_col, watsat_col, t_soisno_col, use_aquifer_layer, exice_coldstart_depth, exice_init_conc_col)

    class(waterstate_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds  
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: tracer_vars
    real(r8)          , intent(in) :: h2osno_input_col(bounds%begc:)
    real(r8)          , intent(in) :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)          , intent(in) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    logical           , intent(in) :: use_aquifer_layer ! whether an aquifer layer is used in this run
    real(r8)          , intent(in) :: exice_coldstart_depth ! depth below which excess ice will be present
    real(r8)          , intent(in) :: exice_init_conc_col(bounds%begc:bounds%endc) ! initial coldstart excess ice concentration (from the stream file)

    this%info => info

    call this%InitAllocate(bounds, tracer_vars)

    call this%InitHistory(bounds, use_aquifer_layer)
    call this%InitCold(bounds = bounds, &
      h2osno_input_col = h2osno_input_col, &
      watsat_col = watsat_col, &
      t_soisno_col = t_soisno_col, &
      use_aquifer_layer = use_aquifer_layer, &
      exice_coldstart_depth = exice_coldstart_depth , exice_init_conc_col = exice_init_conc_col)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, tracer_vars)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterstate_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  
    type(water_tracer_container_type), intent(inout) :: tracer_vars
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------

    call AllocateVar1d(var = this%h2osno_no_layers_col, name = 'h2osno_no_layers_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column)
    call AllocateVar2d(var = this%h2osoi_vol_col, name = 'h2osoi_vol_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column, &
         dim2beg = 1, dim2end = nlevmaxurbgrnd)
    call AllocateVar2d(var = this%h2osoi_vol_prs_grc, name = 'h2osoi_vol_prs_grc', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_gridcell, &
         dim2beg = 1, dim2end = nlevgrnd)
    call AllocateVar2d(var = this%h2osoi_ice_col, name = 'h2osoi_ice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column, &
         dim2beg = -nlevsno+1, dim2end = nlevmaxurbgrnd)
    call AllocateVar2d(var = this%h2osoi_liq_col, name = 'h2osoi_liq_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column, &
         dim2beg = -nlevsno+1, dim2end = nlevmaxurbgrnd)
    call AllocateVar1d(var = this%snocan_patch, name = 'snocan_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_patch)
    call AllocateVar1d(var = this%liqcan_patch, name = 'liqcan_patch', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_patch)
    call AllocateVar1d(var = this%h2osfc_col, name = 'h2osfc_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column)
    call AllocateVar1d(var = this%wa_col, name = 'wa_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column)
    call AllocateVar1d(var = this%dynbal_baseline_liq_col, name = 'dynbal_baseline_liq_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column)
    call AllocateVar1d(var = this%dynbal_baseline_ice_col, name = 'dynbal_baseline_ice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column)
    call AllocateVar1d(var = this%stream_water_volume_lun, name = 'stream_water_volume_lun', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_landunit)
    !excess ice vars
    call AllocateVar2d(var = this%excess_ice_col, name = 'excess_ice_col', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column, &
         dim2beg = -nlevsno+1, dim2end = nlevmaxurbgrnd)
    call AllocateVar1d(var = this%exice_bulk_init, name = 'exice_bulk_init', &
         container = tracer_vars, &
         bounds = bounds, subgrid_level = subgrid_level_column)

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, use_aquifer_layer)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    use clm_varctl     , only : use_soil_moisture_streams
    use GridcellType   , only : grc
    !
    ! !ARGUMENTS:
    class(waterstate_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds  
    logical          , intent(in) :: use_aquifer_layer ! whether an aquifer layer is used in this run
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begl, endl
    integer           :: begg, endg
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    data2dptr => this%h2osoi_liq_col(:,-nlevsno+1:0)
    call hist_addfld2d ( &
         fname=this%info%fname('SNO_LIQH2O'), &
         units='kg/m2', type2d='levsno',  &
         avgflag='A', &
         long_name=this%info%lname('Snow liquid water content'), &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    data2dptr => this%h2osoi_ice_col(:,-nlevsno+1:0)
    call hist_addfld2d ( &
         fname=this%info%fname('SNO_ICE'), &
         units='kg/m2', type2d='levsno',  &
         avgflag='A', &
         long_name=this%info%lname('Snow ice content'), &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    data2dptr => this%h2osoi_vol_col(begc:endc,1:nlevsoi)
    call hist_addfld2d ( &
         fname=this%info%fname('H2OSOI'),  &
         units='mm3/mm3', type2d='levsoi', &
         avgflag='A', &
         long_name=this%info%lname('volumetric soil water (natural vegetated and crop landunits only)'), &
         ptr_col=this%h2osoi_vol_col, l2g_scale_type='veg')

    if ( use_soil_moisture_streams )then
       call hist_addfld2d ( &
            fname=this%info%fname('H2OSOI_PRESCRIBED_GRC'), &
            units='mm3/mm3', type2d='levsoi',  &
            avgflag='A', &
            long_name=this%info%lname('volumetric soil water prescribed (vegetated landunits only)'), &
            ptr_gcell=this%h2osoi_vol_prs_grc, l2g_scale_type='veg',  default='inactive')
    end if

    ! this%h2osoi_liq_col(begc:endc,:) = spval
    ! call hist_addfld2d ( &
    !      fname=this%info%fname('SOILLIQ'),  &
    !      units='kg/m2', type2d='levgrnd', &
    !      avgflag='A', &
    !      long_name=this%info%lname('soil liquid water (natural vegetated and crop landunits only)'), &
    !      ptr_col=this%h2osoi_liq_col, l2g_scale_type='veg')

    data2dptr => this%h2osoi_liq_col(begc:endc,1:nlevsoi) 
    call hist_addfld2d ( &
         fname=this%info%fname('SOILLIQ'),  &
         units='kg/m2', type2d='levsoi', &
         avgflag='A', &
         long_name=this%info%lname('soil liquid water (natural vegetated and crop landunits only)'), &
         ptr_col=data2dptr, l2g_scale_type='veg')

    data2dptr => this%h2osoi_ice_col(begc:endc,1:nlevsoi)
    call hist_addfld2d ( &
         fname=this%info%fname('SOILICE'),  &
         units='kg/m2', type2d='levsoi', &
         avgflag='A', &
         long_name=this%info%lname('soil ice (natural vegetated and crop landunits only)'), &
         ptr_col=data2dptr, l2g_scale_type='veg')

    this%snocan_patch(begp:endp) = spval 
    call hist_addfld1d ( &
         fname=this%info%fname('SNOCAN'), &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('intercepted snow'), &
         ptr_patch=this%snocan_patch, set_lake=0._r8)

    this%liqcan_patch(begp:endp) = spval 
    call hist_addfld1d ( &
         fname=this%info%fname('LIQCAN'), &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('intercepted liquid water'), &
         ptr_patch=this%liqcan_patch, set_lake=0._r8)

    this%h2osfc_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('H2OSFC'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('surface water depth'), &
         ptr_col=this%h2osfc_col)

    if (use_aquifer_layer) then
       this%wa_col(begc:endc) = spval
       call hist_addfld1d (fname=this%info%fname('WA'),  units='mm',  &
            avgflag='A', &
            long_name=this%info%lname('water in the unconfined aquifer (natural vegetated and crop landunits only)'), &
            ptr_col=this%wa_col, l2g_scale_type='veg')
    end if

    if (use_hillslope) then
       this%stream_water_volume_lun(begl:endl) = spval
       call hist_addfld1d (fname=this%info%fname('STREAM_WATER_VOLUME'),  units='m3',  &
            avgflag='A', &
            long_name=this%info%lname('volume of water in stream channel (hillslope hydrology only)'), &
            ptr_lunit=this%stream_water_volume_lun, l2g_scale_type='natveg',  default='inactive')
    end if

    ! Add excess ice fields to history

    if (use_excess_ice) then
       data2dptr => this%excess_ice_col(begc:endc,1:nlevsoi)
       call hist_addfld2d (fname='EXCESS_ICE',  units='kg/m2', type2d='levsoi', &
           avgflag='A', long_name='excess soil ice (vegetated landunits only)', &
           ptr_col=this%excess_ice_col, l2g_scale_type='veg', default = 'inactive')
    end if

    ! (rgk 02-02-2017) There is intentionally no entry  here for stored plant water
    !                  I think that since the value is zero in all cases except
    !                  for FATES plant hydraulics, it will be confusing for users
    !                  when they see their plants have no water in output files.
    !                  So it is not useful diagnostic information. The information
    !                  can be provided through FATES specific history diagnostics
    !                  if need be.


  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       h2osno_input_col, watsat_col, t_soisno_col, use_aquifer_layer, exice_coldstart_depth, exice_init_conc_col)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use landunit_varcon , only : istwet, istsoil, istcrop, istice
    use column_varcon   , only : icol_road_perv, icol_road_imperv
    use clm_varcon      , only : denice, denh2o, bdsno , zisoi
    use clm_varcon      , only : tfrz, aquifer_water_baseline
    use initVerticalMod , only : find_soil_layer_containing_depth
    !
    ! !ARGUMENTS:
    class(waterstate_type), intent(inout) :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: h2osno_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: watsat_col(bounds%begc:, 1:)            ! volumetric soil water at saturation (porosity)
    real(r8)              , intent(in)    :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    logical               , intent(in)    :: use_aquifer_layer                       ! whether an aquifer layer is used in this run
    real(r8)              , intent(in)    :: exice_coldstart_depth ! depth below which excess ice will be present
    real(r8)              , intent(in)    :: exice_init_conc_col(bounds%begc:bounds%endc) ! initial coldstart excess ice concentration (from the stream file)
    !
    ! !LOCAL VARIABLES:
    integer            :: c,j,l,nlevs,g 
    integer            :: nbedrock, nexice ! layer containing 0.5 m
    real(r8)           :: ratio
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(h2osno_input_col)     == (/bounds%endc/))          , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(watsat_col)           == (/bounds%endc,nlevmaxurbgrnd/)) , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_soisno_col)         == (/bounds%endc,nlevmaxurbgrnd/)) , sourcefile, __LINE__)

    ratio = this%info%get_ratio()

    associate(snl => col%snl) 

      this%h2osfc_col(bounds%begc:bounds%endc) = 0._r8
      this%snocan_patch(bounds%begp:bounds%endp) = 0._r8
      this%liqcan_patch(bounds%begp:bounds%endp) = 0._r8
      this%stream_water_volume_lun(bounds%begl:bounds%endl) = 0._r8

      !--------------------------------------------
      ! Set soil water
      !--------------------------------------------

      ! volumetric water is set first and liquid content and ice lens are obtained
      ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
      ! and urban pervious road (other urban columns have zero soil water)

      this%h2osoi_vol_col(bounds%begc:bounds%endc,         1:) = spval
      this%h2osoi_vol_prs_grc(bounds%begg:bounds%endg,     1:) = spval
      this%h2osoi_liq_col(bounds%begc:bounds%endc,-nlevsno+1:) = spval
      this%h2osoi_ice_col(bounds%begc:bounds%endc,-nlevsno+1:) = spval
      do c = bounds%begc,bounds%endc
         l = col%landunit(c)
         if (.not. lun%lakpoi(l)) then  !not lake

            ! volumetric water
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  if (use_bedrock .and. col%nbedrock(c) <=nlevsoi) then
                     nbedrock = col%nbedrock(c)
                  else
                     nbedrock = nlevsoi
                  endif
                  if (j > nbedrock) then
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  else
                     if(use_fates) then
                         this%h2osoi_vol_col(c,j) = 0.75_r8*watsat_col(c,j)*ratio
                     else
                         this%h2osoi_vol_col(c,j) = 0.15_r8*ratio
                     end if
                  endif
               end do
            else if (lun%urbpoi(l)) then
               if (col%itype(c) == icol_road_perv) then
                  nlevs = nlevgrnd
                  do j = 1, nlevs
                     if (j <= nlevsoi) then
                        this%h2osoi_vol_col(c,j) = 0.3_r8 * ratio
                     else
                        this%h2osoi_vol_col(c,j) = 0.0_r8
                     end if
                  end do
               else if (col%itype(c) == icol_road_imperv) then
                  nlevs = nlevgrnd
                  do j = 1, nlevs
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  end do
               else
                  nlevs = nlevurb
                  do j = 1, nlevs
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  end do
               end if
            else if (lun%itype(l) == istwet) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  if (j > nlevsoi) then
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  else
                     this%h2osoi_vol_col(c,j) = 1.0_r8 * ratio
                  endif
               end do
            else if (lun%itype(l) == istice) then
               nlevs = nlevgrnd 
               do j = 1, nlevs
                  this%h2osoi_vol_col(c,j) = 1.0_r8 * ratio
               end do
            else
               write(iulog,*) 'water_state_type InitCold: unhandled landunit type ', lun%itype(l)
               call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg = 'unhandled landunit type', &
                    additional_msg = errMsg(sourcefile, __LINE__))
            endif
            do j = 1, nlevs
               this%h2osoi_vol_col(c,j) = min(this%h2osoi_vol_col(c,j), watsat_col(c,j)*ratio)
               if (t_soisno_col(c,j) <= SHR_CONST_TKFRZ) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*denice*this%h2osoi_vol_col(c,j) ! ratio already applied
                  this%h2osoi_liq_col(c,j) = 0._r8
               else
                  this%h2osoi_ice_col(c,j) = 0._r8
                  this%h2osoi_liq_col(c,j) = col%dz(c,j)*denh2o*this%h2osoi_vol_col(c,j) ! ratio already applied
               endif
            end do
            if (snl(c) == 0) then
               this%h2osno_no_layers_col(c) = h2osno_input_col(c) * ratio
            else
               this%h2osno_no_layers_col(c) = 0._r8
            end if
            do j = -nlevsno+1, 0
               if (j > snl(c)) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*250._r8 * ratio
                  this%h2osoi_liq_col(c,j) = 0._r8
               end if
            end do
         end if
      end do


      !--------------------------------------------
      ! Set Lake water
      !--------------------------------------------

      do c = bounds%begc, bounds%endc
         l = col%landunit(c)

         if (lun%lakpoi(l)) then
            if (snl(c) == 0) then
               this%h2osno_no_layers_col(c) = h2osno_input_col(c) * ratio
            else
               this%h2osno_no_layers_col(c) = 0._r8
            end if
            do j = -nlevsno+1, 0
               if (j > snl(c)) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*bdsno * ratio
                  this%h2osoi_liq_col(c,j) = 0._r8
               end if
            end do
            do j = 1,nlevgrnd
               if (j <= nlevsoi) then ! soil
                  this%h2osoi_vol_col(c,j) = watsat_col(c,j) * ratio
                  this%h2osoi_liq_col(c,j) = spval
                  this%h2osoi_ice_col(c,j) = spval
               else                  ! bedrock
                  this%h2osoi_vol_col(c,j) = 0._r8
               end if
            end do
         end if
      end do

      !--------------------------------------------
      ! For frozen layers !TODO - does the following make sense ???? it seems to overwrite everything
      !--------------------------------------------

      do c = bounds%begc, bounds%endc
         do j = 1,nlevmaxurbgrnd
            if (this%h2osoi_vol_col(c,j) /= spval) then
               if (t_soisno_col(c,j) <= tfrz) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*denice*this%h2osoi_vol_col(c,j) ! ratio already applied
                  this%h2osoi_liq_col(c,j) = 0._r8
               else
                  this%h2osoi_ice_col(c,j) = 0._r8
                  this%h2osoi_liq_col(c,j) = col%dz(c,j)*denh2o*this%h2osoi_vol_col(c,j) ! ratio already applied
               endif
            end if
         end do
      end do


      this%aquifer_water_baseline = aquifer_water_baseline * ratio
      this%wa_col(bounds%begc:bounds%endc)  = this%aquifer_water_baseline
      if (use_aquifer_layer) then
         ! NOTE(wjs, 2018-11-27) There is no fundamental reason why wa_col should be
         ! initialized differently based on use_aquifer_layer, but we (Bill Sacks and Sean
         ! Swenson) want to change the cold start initialization of wa_col to be
         ! aquifer_water_baseline everywhere for use_aquifer_layer .false., and we aren't
         ! sure of the implications of this change for use_aquifer_layer .true., so are
         ! maintaining the old cold start initialization in the latter case.
         do c = bounds%begc,bounds%endc
            l = col%landunit(c)
            if (.not. lun%lakpoi(l)) then  !not lake
               if (lun%urbpoi(l)) then
                  if (col%itype(c) == icol_road_perv) then
                     ! Note that the following hard-coded constant (on the next line)
                     ! seems implicitly related to aquifer_water_baseline 
                     this%wa_col(c)  = 4800._r8 * ratio
                  else
                     this%wa_col(c)  = spval
                  end if
               else
                  ! Note that the following hard-coded constant (on the next line) seems
                  ! implicitly related to aquifer_water_baseline
                  this%wa_col(c)  = 4000._r8 * ratio
               end if
            end if
         end do
      end if

      ! Initialize dynbal_baseline_liq_col and dynbal_baseline_ice_col: for some columns,
      ! these are set elsewhere in initialization, but we need them to be 0 for columns
      ! for which they are not explicitly set.
      this%dynbal_baseline_liq_col(bounds%begc:bounds%endc) = 0._r8
      this%dynbal_baseline_ice_col(bounds%begc:bounds%endc) = 0._r8

      !Initialize excess ice
      this%exice_bulk_init(bounds%begc:bounds%endc) = exice_init_conc_col(bounds%begc:bounds%endc)
      this%excess_ice_col(bounds%begc:bounds%endc,:) = 0.0_r8
      if (use_excess_ice) then
         do c = bounds%begc,bounds%endc
            g = col%gridcell(c)
            l = col%landunit(c)
            if (.not. lun%lakpoi(l)) then  !not lake
               if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
                  if (zisoi(nlevsoi) >= exice_coldstart_depth) then
                     call find_soil_layer_containing_depth(exice_coldstart_depth,nexice)
                  else
                     nexice=nlevsoi-1
                  endif
                  if (use_bedrock .and. col%nbedrock(c) <=nlevsoi) then
                     nbedrock = col%nbedrock(c)
                  else
                     nbedrock = nlevsoi
                  endif
                  do j = 2, nlevmaxurbgrnd ! ignore first layer
                     if (nexice<nbedrock) then ! bedrock below 1 m
                        if (j >= nexice .and. j<nbedrock .and. t_soisno_col(c,j) <= tfrz ) then
                           this%excess_ice_col(c,j) = col%dz(c,j)*denice*(this%exice_bulk_init(c))
                        else
                           this%excess_ice_col(c,j) = 0.0_r8
                        endif
                     else 
                        this%excess_ice_col(c,j) = 0.0_r8
                     end if
                  end do
               endif
            else ! just in case zeros for lakes and other columns
               this%excess_ice_col(c,-nlevsno+1:nlevmaxurbgrnd) = 0.0_r8
            endif
         enddo
      end if
    end associate

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       watsat_col, t_soisno_col, altmax_lastyear_indx)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use clm_varcon       , only : denice, denh2o, pondmx, watmin, tfrz
    use landunit_varcon  , only : istcrop, istdlak, istsoil  
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager , only : is_first_step, is_restart
    use clm_varctl       , only : bound_h2osoi, nsrest, nsrContinue
    use ncdio_pio        , only : file_desc_t, ncd_double
    use ExcessIceStreamType, only : UseExcessIceStreams
    use restUtilMod        , only : restartvar, RestartExcessIceIssue
    !
    ! !ARGUMENTS:
    class(waterstate_type), intent(in) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid                                    ! netcdf id
    character(len=*) , intent(in)    :: flag                                    ! 'read' or 'write'
    real(r8)         , intent(in)    :: watsat_col (bounds%begc:, 1:)           ! volumetric soil water at saturation (porosity)
    real(r8)         , intent(in)    :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    integer          , intent(in)    :: altmax_lastyear_indx(bounds%begc:)      !col active layer index last year
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,l,j,nlevs,nbedrock
    logical  :: readvar
    real(r8) :: maxwatsat    ! maximum porosity    
    real(r8) :: excess       ! excess volumetric soil water
    real(r8) :: totwat       ! total soil water (mm)
    logical :: excess_ice_on_restart ! Excess ice fields are on the restart file
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(watsat_col) == (/bounds%endc,nlevmaxurbgrnd/)) , sourcefile, __LINE__)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('H2OSFC'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('surface water'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osfc_col)
    if (flag=='read' .and. .not. readvar) then
       this%h2osfc_col(bounds%begc:bounds%endc) = 0.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('H2OSNO_NO_LAYERS')//':'//this%info%fname('H2OSNO'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('snow that is not resolved into layers'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osno_no_layers_col)
    ! BACKWARDS_COMPATIBILITY(wjs, 2019-06-06) If h2osno_no_layers is read from the old
    ! h2osno, then it will be non-zero for an explicit-layered snow pack. We fix that
    ! here. We can (and should) remove this backwards compatibility code at the same time
    ! as we remove ":H2OSNO" from the restart variable name above.
    if (flag == 'read' .and. .not. is_restart()) then
       do c = bounds%begc, bounds%endc
          if (col%snl(c) < 0) then
             this%h2osno_no_layers_col(c) = 0._r8
          end if
       end do
    end if

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('H2OSOI_LIQ'), &
         xtype=ncd_double,  &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name=this%info%lname('liquid water'), &
         units='kg/m2', &
         scale_by_thickness=.true., &
         interpinic_flag='interp', readvar=readvar, data=this%h2osoi_liq_col)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('H2OSOI_ICE'), &
         xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name=this%info%lname('ice lens'), &
         units='kg/m2', &
         scale_by_thickness=.true., &
         interpinic_flag='interp', readvar=readvar, data=this%h2osoi_ice_col)
         
    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('SNOCAN'), &
         xtype=ncd_double,  &
         dim1name='pft', &
         long_name=this%info%lname('canopy snow water'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%snocan_patch)

    ! NOTE(wjs, 2015-07-01) In old restart files, there was no LIQCAN variable. However,
    ! H2OCAN had similar meaning. So if we can't find LIQCAN, use H2OCAN to initialize
    ! liqcan_patch.
    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('LIQCAN')//':'//this%info%fname('H2OCAN'), &
         xtype=ncd_double,  &
         dim1name='pft', &
         long_name=this%info%lname('canopy liquid water'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%liqcan_patch)

    call restartvar(ncid=ncid, flag=flag, varname=this%info%fname('WA'), xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('water in the unconfined aquifer'), units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%wa_col)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('DYNBAL_BASELINE_LIQ'), &
         xtype=ncd_double, &
         dim1name='column', &
         long_name=this%info%lname("baseline liquid water mass subtracted from each column's total water calculation"), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%dynbal_baseline_liq_col)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('DYNBAL_BASELINE_ICE'), &
         xtype=ncd_double, &
         dim1name='column', &
         long_name=this%info%lname("baseline ice mass subtracted from each column's total ice calculation"), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%dynbal_baseline_ice_col)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('STREAM_WATER_VOLUME'), &
         xtype=ncd_double,  &
         dim1name='landunit', &
         long_name=this%info%lname('water in stream channel'), &
         units='m3', &
         interpinic_flag='interp', readvar=readvar, data=this%stream_water_volume_lun)
    ! Restart excess ice vars
    if (.not. use_excess_ice) then
       ! no need to even define the restart vars
       call RestartExcessIceIssue( &
            ncid = ncid, &
            flag = flag, &
            excess_ice_on_restart = excess_ice_on_restart)
       if( excess_ice_on_restart ) then
          if (masterproc) then
             write(iulog,*) '--WARNING-- Starting from initial conditions with excess ice present.'
             write(iulog,*) 'But use_excess_ice=.false.'
             write(iulog,*) 'This will cause soil moisture and temperature not being in equilibrium'
          endif
       endif
       this%excess_ice_col(bounds%begc:bounds%endc,-nlevsno+1:nlevmaxurbgrnd)=0.0_r8
    else
       call RestartExcessIceIssue( &
            ncid = ncid, &
            flag = flag, &
            excess_ice_on_restart = excess_ice_on_restart)
       ! have to at least define them 
       call restartvar(ncid=ncid, flag=flag, varname=this%info%fname('EXCESS_ICE'), xtype=ncd_double,  &
            dim1name='column', dim2name='levtot', switchdim=.true., &
            long_name=this%info%lname('excess soil ice (vegetated landunits only)'), &
            units='kg/m2', scale_by_thickness=.true., &
            interpinic_flag='interp', readvar=readvar, data=this%excess_ice_col)
       if (flag == 'read' .and. ((.not. readvar) .or. (.not.  excess_ice_on_restart)) ) then ! when reading restart that does not have excess ice in it
          if (nsrest == nsrContinue) then
             call endrun(msg = "On a continue run, excess ice fields MUST be on the restart file "// & 
             errMsg(sourcefile, __LINE__))
          else if ( .not. UseExcessIceStreams() )then
             call endrun(msg = "This input initial conditions file does NOT include excess ice fields" // &
                         ", and use_excess_ice_streams is off, one or the other needs to be changed  "// & 
                         errMsg(sourcefile, __LINE__))
          end if
          if (masterproc) then
             write(iulog,*) 'Excess ice data is read from the stream and not from restart file!'
          endif
          do c = bounds%begc,bounds%endc
          l = col%landunit(c)
          if (.not. lun%lakpoi(l)) then  !not lake
             if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
                if (use_bedrock .and. col%nbedrock(c)>nlevsoi) then
                   nbedrock = col%nbedrock(c)
                else
                   nbedrock = nlevsoi
                end if
                do j = 2, nlevmaxurbgrnd ! ignore first layer
                   if(altmax_lastyear_indx(c) < nbedrock) then
                      if (j>altmax_lastyear_indx(c) .and. j<nbedrock &
                           .and. t_soisno_col(c,j) <= tfrz) then
                         this%excess_ice_col(c,j) = col%dz(c,j)*denice*(this%exice_bulk_init(c)) ! exice_bulk_init should be already read from the stream during InitCold
                      else
                         this%excess_ice_col(c,j) = 0.0_r8
                      end if
                   else
                      this%excess_ice_col(c,j) = 0.0_r8
                   end if
                end do
             else
                this%excess_ice_col(c,-nlevsno+1:nlevmaxurbgrnd) = 0.0_r8
             end if
          end if
          end do ! end of column loop
       endif! end of old file restart
    endif ! end of exice restart

    ! Determine volumetric soil water (for read only)
    if (flag == 'read' ) then
       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          if ( col%itype(c) == icol_sunwall   .or. &
               col%itype(c) == icol_shadewall .or. &
               col%itype(c) == icol_roof )then
             nlevs = nlevurb
          else
             nlevs = nlevgrnd
          end if
          if ( lun%itype(l) /= istdlak ) then ! This calculation is now done for lakes in initLake.
             do j = 1,nlevs
                this%h2osoi_vol_col(c,j) = this%h2osoi_liq_col(c,j)/(col%dz(c,j)*denh2o) &
                                         + this%h2osoi_ice_col(c,j)/(col%dz(c,j)*denice)
             end do
          end if
       end do
    end if

    ! If initial run -- ensure that water is properly bounded (read only)
    if (flag == 'read' ) then
       if ( is_first_step() .and. bound_h2osoi) then
          do c = bounds%begc, bounds%endc
             l = col%landunit(c)
             if ( col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall .or. &
                  col%itype(c) == icol_roof )then
                nlevs = nlevurb
             else
                nlevs = nlevgrnd
             end if
             do j = 1,nlevs
                l = col%landunit(c)
                if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
                   this%h2osoi_liq_col(c,j) = max(0._r8,this%h2osoi_liq_col(c,j))
                   this%h2osoi_ice_col(c,j) = max(0._r8,this%h2osoi_ice_col(c,j))
                   this%h2osoi_vol_col(c,j) = this%h2osoi_liq_col(c,j)/(col%dz(c,j)*denh2o) &
                                       + this%h2osoi_ice_col(c,j)/(col%dz(c,j)*denice)
                   if (j == 1) then
                      maxwatsat = (watsat_col(c,j)*col%dz(c,j)*1000.0_r8 + pondmx) / (col%dz(c,j)*1000.0_r8)
                   else
                      maxwatsat =  watsat_col(c,j)
                   end if
                   if (this%h2osoi_vol_col(c,j) > maxwatsat) then 
                      excess = (this%h2osoi_vol_col(c,j) - maxwatsat)*col%dz(c,j)*1000.0_r8
                      totwat = this%h2osoi_liq_col(c,j) + this%h2osoi_ice_col(c,j)
                      this%h2osoi_liq_col(c,j) = this%h2osoi_liq_col(c,j) - &
                                           (this%h2osoi_liq_col(c,j)/totwat) * excess
                      this%h2osoi_ice_col(c,j) = this%h2osoi_ice_col(c,j) - &
                                           (this%h2osoi_ice_col(c,j)/totwat) * excess
                   end if
                   this%h2osoi_liq_col(c,j) = max(watmin,this%h2osoi_liq_col(c,j))
                   this%h2osoi_ice_col(c,j) = max(watmin,this%h2osoi_ice_col(c,j))
                   this%h2osoi_vol_col(c,j) = this%h2osoi_liq_col(c,j)/(col%dz(c,j)*denh2o) &
                                             + this%h2osoi_ice_col(c,j)/(col%dz(c,j)*denice)
                end if
             end do
          end do
       end if

    endif   ! end if if-read flag

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine CalculateTotalH2osno(this, &
       bounds, num_c, filter_c, caller, &
       h2osno_total)
    !
    ! !DESCRIPTION:
    ! Calculate h2osno_total over the given column filter
    !
    ! If running in debug mode, also assert that we don't have any unresolved snow if snl
    ! < 0, and that we don't have any resolved snow if snl == 0.
    !
    ! !ARGUMENTS:
    class(waterstate_type) , intent(in)    :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_c                        ! number of columns in filter
    integer                , intent(in)    :: filter_c(:)                  ! filter for columns to operate over
    character(len=*)       , intent(in)    :: caller                       ! name of caller (used in error messages)
    real(r8)               , intent(inout) :: h2osno_total( bounds%begc: ) ! total snow water (mm H2O)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    integer :: j

    character(len=*), parameter :: subname = 'CalculateTotalH2osno'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(h2osno_total, 1) == bounds%endc), sourcefile, __LINE__)

#ifndef NDEBUG
    call this%CheckSnowConsistency(num_c, filter_c, caller)
#endif

    do fc = 1, num_c
       c = filter_c(fc)
       h2osno_total(c) = this%h2osno_no_layers_col(c)

       do j = col%snl(c)+1, 0
          h2osno_total(c) = &
               h2osno_total(c) + &
               this%h2osoi_ice_col(c,j) + &
               this%h2osoi_liq_col(c,j)
       end do
    end do

  end subroutine CalculateTotalH2osno

  !-----------------------------------------------------------------------
  subroutine CheckSnowConsistency(this, num_c, filter_c, caller)
    !
    ! !DESCRIPTION:
    ! Make sure we only have unresolved snow where we should, and that we only have
    ! resolved snow where we should.
    !
    ! !ARGUMENTS:
    class(waterstate_type) , intent(in) :: this
    integer                , intent(in) :: num_c       ! number of columns in filter
    integer                , intent(in) :: filter_c(:) ! filter for columns to operate over
    character(len=*)       , intent(in) :: caller      ! name of caller (used in error messages)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c
    integer :: j
    logical :: ice_bad
    logical :: liq_bad

    character(len=*), parameter :: subname = 'CheckSnowConsistency'
    !-----------------------------------------------------------------------

    do fc = 1, num_c
       c = filter_c(fc)
       if (col%snl(c) < 0) then
          if (this%h2osno_no_layers_col(c) /= 0._r8) then
             write(iulog,*) subname//' ERROR: col has snow layers but non-zero h2osno_no_layers'
             write(iulog,*) '(Called from: ', trim(caller), ')'
             write(iulog,*) 'c, snl, h2osno_no_layers = ', c, col%snl(c), &
                  this%h2osno_no_layers_col(c)
             call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                  msg = subname//' ERROR: col has snow layers but non-zero h2osno_no_layers')
          end if
       end if

       do j = -nlevsno+1, col%snl(c)
          ice_bad = (this%h2osoi_ice_col(c,j) /= 0._r8 .and. this%h2osoi_ice_col(c,j) /= spval)
          liq_bad = (this%h2osoi_liq_col(c,j) /= 0._r8 .and. this%h2osoi_liq_col(c,j) /= spval)
          if (ice_bad .or. liq_bad) then
             write(iulog,*) subname//' ERROR: col has non-zero h2osoi_ice or h2osoi_liq outside resolved snow layers'
             write(iulog,*) '(Called from: ', trim(caller), ')'
             write(iulog,*) 'c, j, snl, h2osoi_ice, h2osoi_liq = ', c, j, col%snl(c), &
                  this%h2osoi_ice_col(c,j), this%h2osoi_liq_col(c,j)
             call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, &
                  msg = subname//' ERROR: col has non-zero h2osoi_ice or h2osoi_liq outside resolved snow layers')
          end if
       end do
    end do

  end subroutine CheckSnowConsistency

end module WaterStateType
