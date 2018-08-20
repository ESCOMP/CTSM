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
  use decompMod      , only : bounds_type
  use clm_varctl     , only : iulog
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use WaterInfoBaseType, only : water_info_base_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterstate_type

     class(water_info_base_type), pointer :: info

     real(r8), pointer :: h2osno_col             (:)   ! col snow water (mm H2O)
     real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_vol_col         (:,:) ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: h2ocan_patch           (:)   ! patch canopy water (mm H2O)
     real(r8), pointer :: h2osfc_col             (:)   ! col surface water (mm H2O)
     real(r8), pointer :: snocan_patch           (:)   ! patch canopy snow water (mm H2O)
     real(r8), pointer :: liqcan_patch           (:)   ! patch canopy liquid water (mm H2O)


   contains

     procedure          :: Init         
     procedure          :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type waterstate_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, info, &
       h2osno_input_col, watsat_col, t_soisno_col)

    class(waterstate_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    class(water_info_base_type), intent(in), target :: info
    real(r8)          , intent(in) :: h2osno_input_col(bounds%begc:)
    real(r8)          , intent(in) :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)          , intent(in) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)

    this%info => info

    call this%InitAllocate(bounds)

    call this%InitHistory(bounds)

    call this%InitCold(bounds, &
       h2osno_input_col, watsat_col, t_soisno_col)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate(this%h2osno_col             (begc:endc))                     ; this%h2osno_col             (:)   = nan   
    allocate(this%h2osoi_vol_col         (begc:endc, 1:nlevgrnd))         ; this%h2osoi_vol_col         (:,:) = nan
    allocate(this%h2osoi_ice_col         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_col         (:,:) = nan
    allocate(this%h2osoi_liq_col         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_col         (:,:) = nan
    allocate(this%h2ocan_patch           (begp:endp))                     ; this%h2ocan_patch           (:)   = nan  
    allocate(this%snocan_patch           (begp:endp))                     ; this%snocan_patch           (:)   = nan  
    allocate(this%liqcan_patch           (begp:endp))                     ; this%liqcan_patch           (:)   = nan  
    allocate(this%h2osfc_col             (begc:endc))                     ; this%h2osfc_col             (:)   = nan   




  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_lch4
    use clm_varctl     , only : hist_wrtch4diag
    use clm_varpar     , only : nlevsno, nlevsoi
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal, no_snow_zero
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    ! h2osno also includes snow that is part of the soil column (an 
    ! initial snow layer is only created if h2osno > 10mm). 

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
         long_name=this%info%lname('volumetric soil water (vegetated landunits only)'), &
         ptr_col=this%h2osoi_vol_col, l2g_scale_type='veg')

    ! this%h2osoi_liq_col(begc:endc,:) = spval
    ! call hist_addfld2d ( &
    !      fname=this%info%fname('SOILLIQ'),  &
    !      units='kg/m2', type2d='levgrnd', &
    !      avgflag='A', &
    !      long_name=this%info%lname('soil liquid water (vegetated landunits only)'), &
    !      ptr_col=this%h2osoi_liq_col, l2g_scale_type='veg')

    data2dptr => this%h2osoi_liq_col(begc:endc,1:nlevsoi) 
    call hist_addfld2d ( &
         fname=this%info%fname('SOILLIQ'),  &
         units='kg/m2', type2d='levsoi', &
         avgflag='A', &
         long_name=this%info%lname('soil liquid water (vegetated landunits only)'), &
         ptr_col=data2dptr, l2g_scale_type='veg')

    data2dptr => this%h2osoi_ice_col(begc:endc,1:nlevsoi)
    call hist_addfld2d ( &
         fname=this%info%fname('SOILICE'),  &
         units='kg/m2', type2d='levsoi', &
         avgflag='A', &
         long_name=this%info%lname('soil ice (vegetated landunits only)'), &
         ptr_col=data2dptr, l2g_scale_type='veg')

    this%h2ocan_patch(begp:endp) = spval 
    call hist_addfld1d ( &
         fname=this%info%fname('H2OCAN'), &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('intercepted water'), &
         ptr_patch=this%h2ocan_patch, set_lake=0._r8)

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

    call hist_addfld1d ( &
         fname=this%info%fname('H2OSNO'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('snow depth (liquid water)'), &
         ptr_col=this%h2osno_col, c2l_scale_type='urbanf')

    call hist_addfld1d ( &
         fname=this%info%fname('H2OSNO_ICE'), &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('snow depth (liquid water, ice landunits only)'), &
         ptr_col=this%h2osno_col, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')


    this%h2osfc_col(begc:endc) = spval
    call hist_addfld1d ( &
         fname=this%info%fname('H2OSFC'),  &
         units='mm',  &
         avgflag='A', &
         long_name=this%info%lname('surface water depth'), &
         ptr_col=this%h2osfc_col)


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
       h2osno_input_col, watsat_col, t_soisno_col)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    use shr_const_mod   , only : shr_const_pi
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use shr_spfn_mod    , only : shr_spfn_erf
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use landunit_varcon , only : istwet, istsoil, istdlak, istcrop, istice_mec  
    use column_varcon   , only : icol_shadewall, icol_road_perv
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use clm_varcon      , only : denice, denh2o, spval, sb, bdsno 
    use clm_varcon      , only : zlnd, tfrz, spval, pc
    use clm_varctl      , only : fsurdat, iulog
    use clm_varctl        , only : use_bedrock
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    class(waterstate_type)                :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: h2osno_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)              , intent(in)    :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev,nlevs 
    real(r8)           :: maxslope, slopemax, minslope
    real(r8)           :: d, fd, dfdd, slope0,slopebeta
    real(r8) ,pointer  :: std (:)     
    logical            :: readvar 
    type(file_desc_t)  :: ncid        
    character(len=256) :: locfn       
    integer            :: nbedrock
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osno_input_col)     == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(watsat_col)           == (/bounds%endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno_col)         == (/bounds%endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))

    do c = bounds%begc,bounds%endc
       this%h2osno_col(c)             = h2osno_input_col(c) 
    end do



    associate(snl => col%snl) 

      this%h2osfc_col(bounds%begc:bounds%endc) = 0._r8
      this%h2ocan_patch(bounds%begp:bounds%endp) = 0._r8
      this%snocan_patch(bounds%begp:bounds%endp) = 0._r8
      this%liqcan_patch(bounds%begp:bounds%endp) = 0._r8


      !--------------------------------------------
      ! Set soil water
      !--------------------------------------------

      ! volumetric water is set first and liquid content and ice lens are obtained
      ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
      ! and urban pervious road (other urban columns have zero soil water)

      this%h2osoi_vol_col(bounds%begc:bounds%endc,         1:) = spval
      this%h2osoi_liq_col(bounds%begc:bounds%endc,-nlevsno+1:) = spval
      this%h2osoi_ice_col(bounds%begc:bounds%endc,-nlevsno+1:) = spval
      do c = bounds%begc,bounds%endc
         l = col%landunit(c)
         if (.not. lun%lakpoi(l)) then  !not lake

            ! volumetric water
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  if (use_bedrock) then
                     nbedrock = col%nbedrock(c)
                  else
                     nbedrock = nlevsoi
                  endif
                  if (j > nbedrock) then
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  else
                     this%h2osoi_vol_col(c,j) = 0.15_r8
                  endif
               end do
            else if (lun%urbpoi(l)) then
               if (col%itype(c) == icol_road_perv) then
                  nlevs = nlevgrnd
                  do j = 1, nlevs
                     if (j <= nlevsoi) then
                        this%h2osoi_vol_col(c,j) = 0.3_r8
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
                     this%h2osoi_vol_col(c,j) = 1.0_r8
                  endif
               end do
            else if (lun%itype(l) == istice_mec) then
               nlevs = nlevgrnd 
               do j = 1, nlevs
                  this%h2osoi_vol_col(c,j) = 1.0_r8
               end do
            endif
            do j = 1, nlevs
               this%h2osoi_vol_col(c,j) = min(this%h2osoi_vol_col(c,j), watsat_col(c,j))
               if (t_soisno_col(c,j) <= SHR_CONST_TKFRZ) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*denice*this%h2osoi_vol_col(c,j)
                  this%h2osoi_liq_col(c,j) = 0._r8
               else
                  this%h2osoi_ice_col(c,j) = 0._r8
                  this%h2osoi_liq_col(c,j) = col%dz(c,j)*denh2o*this%h2osoi_vol_col(c,j)
               endif
            end do
            do j = -nlevsno+1, 0
               if (j > snl(c)) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*250._r8
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
            do j = -nlevsno+1, 0
               if (j > snl(c)) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*bdsno
                  this%h2osoi_liq_col(c,j) = 0._r8
               end if
            end do
            do j = 1,nlevgrnd
               if (j <= nlevsoi) then ! soil
                  this%h2osoi_vol_col(c,j) = watsat_col(c,j)
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
         do j = 1,nlevgrnd
            if (this%h2osoi_vol_col(c,j) /= spval) then
               if (t_soisno_col(c,j) <= tfrz) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*denice*this%h2osoi_vol_col(c,j)
                  this%h2osoi_liq_col(c,j) = 0._r8
               else
                  this%h2osoi_ice_col(c,j) = 0._r8
                  this%h2osoi_liq_col(c,j) = col%dz(c,j)*denh2o*this%h2osoi_vol_col(c,j)
               endif
            end if
         end do
      end do

    end associate

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       watsat_col)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod          , only : masterproc
    use clm_varcon       , only : denice, denh2o, pondmx, watmin, spval, nameg
    use landunit_varcon  , only : istcrop, istdlak, istsoil  
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager , only : is_first_step
    use clm_varctl       , only : bound_h2osoi
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    real(r8)         , intent(in)    :: watsat_col (bounds%begc:, 1:)  ! volumetric soil water at saturation (porosity)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,j,nlevs
    logical  :: readvar
    real(r8) :: maxwatsat    ! maximum porosity    
    real(r8) :: excess       ! excess volumetric soil water
    real(r8) :: totwat       ! total soil water (mm)
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(watsat_col) == (/bounds%endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))

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
         varname=this%info%fname('H2OSNO'), &
         xtype=ncd_double,  &
         dim1name='column', &
         long_name=this%info%lname('snow water'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osno_col)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('H2OSOI_LIQ'), &
         xtype=ncd_double,  &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name=this%info%lname('liquid water'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osoi_liq_col)

    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('H2OSOI_ICE'), &
         xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name=this%info%lname('ice lens'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osoi_ice_col)
         
    call restartvar(ncid=ncid, flag=flag, &
         varname=this%info%fname('H2OCAN'), &
         xtype=ncd_double,  &
         dim1name='pft', &
         long_name=this%info%lname('canopy water'), &
         units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2ocan_patch)

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


end module WaterStateType
