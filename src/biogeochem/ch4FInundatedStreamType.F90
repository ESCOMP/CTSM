
module ch4FInundatedStreamType

#include "shr_assert.h"

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading in finundated streams file for methane code.
  !
  ! !USES
  use shr_kind_mod   , only: r8 => shr_kind_r8, CL => shr_kind_cl
  use spmdMod        , only: mpicom, masterproc
  use clm_varctl     , only: iulog
  use abortutils     , only: endrun
  use decompMod      , only: bounds_type
  use ch4varcon      , only: finundation_mtd

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: ch4finundatedstream_type
     real(r8), pointer, private :: zwt0_gdc  (:)         ! col coefficient for determining finundated (m)
     real(r8), pointer, private :: f0_gdc    (:)         ! col maximum inundated fraction for a gridcell (for methane code)
     real(r8), pointer, private :: p3_gdc    (:)         ! col coefficient for determining finundated (m)
     real(r8), pointer, private :: fws_slope_gdc     (:) ! col slope in fws = slope * tws + intercept (A coefficient)
     real(r8), pointer, private :: fws_intercept_gdc (:) ! col slope in fws = slope * tws + intercept (B coefficient)
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: Init            ! Initialize and read data in
      procedure, public :: CalcFinundated  ! Calculate finundated based on input streams
      procedure, public :: UseStreams      ! If streams will be used

      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: InitAllocate   ! Allocate data

  end type ch4finundatedstream_type


  ! ! PRIVATE DATA:

  type, private :: streamcontrol_type
     character(len=CL)  :: stream_fldFileName_ch4finundated   ! Filename
     character(len=CL)  :: ch4finundatedmapalgo               ! map algo
     character(len=CL)  :: fldList                            ! List of fields to read
  contains
     procedure, private :: ReadNML     ! Read in namelist
  end type streamcontrol_type

  type(streamcontrol_type), private :: control        ! Stream control data

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================

  subroutine Init(this, bounds, NLFilename)
   !    
   ! Initialize the ch4 finundated stream object
   !
   ! Uses:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar, get_curr_date
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   use ndepStreamMod    , only : clm_domain_mct
   use domainMod        , only : ldomain
   use decompMod        , only : bounds_type, gsmap_lnd_gdc2glo
   use mct_mod          , only : mct_ggrid, mct_avect_indexra
   use shr_strdata_mod  , only : shr_strdata_type, shr_strdata_create
   use shr_strdata_mod  , only : shr_strdata_print, shr_strdata_advance
   use spmdMod          , only : comp_id, iam
   use ch4varcon        , only : finundation_mtd_h2osfc
   use ch4varcon        , only : finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
   !
   ! arguments
   implicit none
   class(ch4finundatedstream_type) :: this
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: ig, g            ! Indices
   type(mct_ggrid)    :: dom_clm          ! domain information 
   type(shr_strdata_type) :: sdat         ! input data stream
   integer            :: index_ZWT0       = 0 ! Index of ZWT0 field
   integer            :: index_F0         = 0 ! Index of F0 field
   integer            :: index_P3         = 0 ! Index of P3 field
   integer            :: index_FWS_TWS_A  = 0 ! Index of FWS_TWS_A field
   integer            :: index_FWS_TWS_B  = 0 ! Index of FWS_TWS_B field
   integer            :: year                 ! year (0, ...) for nstep+1
   integer            :: mon                  ! month (1, ..., 12) for nstep+1
   integer            :: day                  ! day of month (1, ..., 31) for nstep+1
   integer            :: sec                  ! seconds into current date for nstep+1
   integer            :: mcdate               ! Current model date (yyyymmdd)
   character(len=*), parameter :: stream_name = 'ch4finundated'
   character(*), parameter :: subName = "('ch4finundatedstream::Init')"
   !-----------------------------------------------------------------------
   if ( finundation_mtd /= finundation_mtd_h2osfc )then
      call this%InitAllocate( bounds )
      call control%ReadNML( bounds, NLFileName )

      if ( this%useStreams() )then
         call clm_domain_mct (bounds, dom_clm)

         call shr_strdata_create(sdat,name=stream_name,&
           pio_subsystem=pio_subsystem,               & 
           pio_iotype=shr_pio_getiotype(inst_name),   &
           mpicom=mpicom, compid=comp_id,             &
           gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,    &
           nxg=ldomain%ni, nyg=ldomain%nj,            &
           yearFirst=1996,                            &
           yearLast=1996,                             &
           yearAlign=1,                               &
           offset=0,                                  &
           domFilePath='',                            &
           domFileName=trim(control%stream_fldFileName_ch4finundated), &
           domTvarName='time',                        &
           domXvarName='LONGXY' ,                     &
           domYvarName='LATIXY' ,                     &  
           domAreaName='AREA',                        &
           domMaskName='LANDMASK',                    &
           filePath='',                               &
           filename=(/trim(control%stream_fldFileName_ch4finundated)/),&
           fldListFile=control%fldList,               &
           fldListModel=control%fldList,              &
           fillalgo='none',                           &
           mapalgo=control%ch4finundatedmapalgo,      &
           calendar=get_calendar(),                   &
           taxmode='extend'                           )

         if (masterproc) then
            call shr_strdata_print(sdat,'CLM '//stream_name//' data')
         endif

         if( finundation_mtd == finundation_mtd_ZWT_inversion )then
            index_ZWT0      = mct_avect_indexra(sdat%avs(1),'ZWT0')
            index_F0        = mct_avect_indexra(sdat%avs(1),'F0' )
            index_P3        = mct_avect_indexra(sdat%avs(1),'P3' )
         else if( finundation_mtd == finundation_mtd_TWS_inversion )then
            index_FWS_TWS_A = mct_avect_indexra(sdat%avs(1),'FWS_TWS_A')
            index_FWS_TWS_B = mct_avect_indexra(sdat%avs(1),'FWS_TWS_B')
         end if


         ! Explicitly set current date to a hardcoded constant value. Otherwise
         ! using the real date can cause roundoff differences that are
         ! detrected as issues with exact restart.  EBK M05/20/2017
         !call get_curr_date(year, mon, day, sec)
         year = 1996
         mon  = 12
         day  = 31
         sec  = 0
         mcdate = year*10000 + mon*100 + day

         call shr_strdata_advance(sdat, mcdate, sec, mpicom, 'ch4finundated')

         ! Get the data
         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
            if ( index_ZWT0 > 0 )then
               this%zwt0_gdc(g) = sdat%avs(1)%rAttr(index_ZWT0,ig)
            end if
            if ( index_F0 > 0 )then
               this%f0_gdc(g) = sdat%avs(1)%rAttr(index_F0,ig)
            end if
            if ( index_P3 > 0 )then
               this%p3_gdc(g) = sdat%avs(1)%rAttr(index_P3,ig)
            end if
            if ( index_FWS_TWS_A > 0 )then
               this%fws_slope_gdc(g) = sdat%avs(1)%rAttr(index_FWS_TWS_A,ig)
            end if
            if ( index_FWS_TWS_B > 0 )then
               this%fws_intercept_gdc(g) = sdat%avs(1)%rAttr(index_FWS_TWS_B,ig)
            end if
         end do
      end if
   end if

  end subroutine Init

  !-----------------------------------------------------------------------
  logical function UseStreams(this)
    !
    ! !DESCRIPTION:
    ! Return true if 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(ch4finundatedstream_type) :: this
    !
    ! !LOCAL VARIABLES:
    if ( trim(control%stream_fldFileName_ch4finundated) == '' )then
       UseStreams = .false. 
    else
       UseStreams = .true. 
    end if
  end function UseStreams

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use ch4varcon     , only: finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
    !
    ! !ARGUMENTS:
    implicit none
    class(ch4finundatedstream_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begc, endc
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    if( finundation_mtd == finundation_mtd_ZWT_inversion )then
       allocate(this%zwt0_gdc         (begg:endg))            ;  this%zwt0_gdc         (:)   = nan
       allocate(this%f0_gdc           (begg:endg))            ;  this%f0_gdc           (:)   = nan
       allocate(this%p3_gdc           (begg:endg))            ;  this%p3_gdc           (:)   = nan
    else if( finundation_mtd == finundation_mtd_TWS_inversion )then
       allocate(this%fws_slope_gdc    (begg:endg))            ;  this%fws_slope_gdc    (:)   = nan
       allocate(this%fws_intercept_gdc(begg:endg))            ;  this%fws_intercept_gdc(:)   = nan
    end if

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine CalcFinundated(this, bounds, num_soilc, filter_soilc, soilhydrology_inst, &
                            waterdiagnosticbulk_inst, qflx_surf_lag_col, finundated )
    !
    ! !DESCRIPTION:
    ! 
    ! Calculate finundated according to the appropriate methodology
    !
    ! !USES:
    use ColumnType       , only : col
    use ch4varcon        , only : finundation_mtd_h2osfc, finundation_mtd_ZWT_inversion
    use ch4varcon        , only : finundation_mtd_TWS_inversion
    use clm_varpar       , only : nlevsoi
    use SoilHydrologyType, only : soilhydrology_type
    use WaterDiagnosticBulkType   , only : waterdiagnosticbulk_type
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    implicit none
    class(ch4finundatedstream_type)                        :: this
    type(bounds_type)                      , intent(in)    :: bounds
    integer                                , intent(in)    :: num_soilc          ! number of column soil points in column filter
    integer                                , intent(in)    :: filter_soilc(:)    ! column filter for soil points
    type(soilhydrology_type)               , intent(in)    :: soilhydrology_inst
    type(waterdiagnosticbulk_type)                  , intent(in)    :: waterdiagnosticbulk_inst
    real(r8)                               , intent(in)    :: qflx_surf_lag_col(bounds%begc:) !time-lagged surface runoff (mm H2O /s)
    real(r8)                               , intent(inout) :: finundated(bounds%begc:)        ! fractional inundated area in soil column (excluding dedicated wetland columns)
    !
    ! !LOCAL VARIABLES:
    integer  :: g, c, fc    ! Indices
    real(r8) :: zwt_actual  ! Total water storage (ZWT) to use either perched or total depending on conditions

    SHR_ASSERT_ALL((ubound(qflx_surf_lag_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(finundated)        == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate(                                                                 &
         z                    =>   col%z                                     , & ! Input:  [real(r8) (:,:) ]  layer depth (m) (-nlevsno+1:nlevsoi)
         zwt                  =>   soilhydrology_inst%zwt_col                , & ! Input:  [real(r8) (:)   ]  water table depth (m)
         zwt_perched          =>   soilhydrology_inst%zwt_perched_col        , & ! Input:  [real(r8) (:)   ]  perched water table depth (m)
         tws                  =>   waterdiagnosticbulk_inst%tws_grc                   , & ! Input:  [real(r8) (:)   ]  total water storage (kg m-2)
         frac_h2osfc          =>   waterdiagnosticbulk_inst%frac_h2osfc_col             & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
    )

    ! Calculate finundated
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       g = col%gridcell(c)
       select case( finundation_mtd )
          case ( finundation_mtd_h2osfc )
             finundated(c) = frac_h2osfc(c)
          case ( finundation_mtd_ZWT_inversion )
             if (this%zwt0_gdc(g) > 0._r8) then
                if (zwt_perched(c) < z(c,nlevsoi)-1.e-5_r8 .and. zwt_perched(c) < zwt(c)) then
                   zwt_actual = zwt_perched(c)
                else
                   zwt_actual = zwt(c)
                end if
                finundated(c) = this%f0_gdc(g) * exp(-zwt_actual/this%zwt0_gdc(g)) + this%p3_gdc(g)*qflx_surf_lag_col(c)
             else
                finundated(c) = this%p3_gdc(g)*qflx_surf_lag_col(c)
             end if
          case ( finundation_mtd_TWS_inversion )
             finundated(c) = this%fws_slope_gdc(g) * tws(g) + this%fws_intercept_gdc(g)
       end select
       finundated(c) = min( 1.0_r8, max( 0.0_r8, finundated(c) ) )
    end do
    end associate

  end subroutine CalcFinundated
  !==============================================================================

  subroutine ReadNML(this, bounds, NLFilename)
   !    
   ! Read the namelist data stream information.  
   !
   ! Uses:
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   use fileutils        , only : getavu, relavu
   use ch4varcon        , only : finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
   !
   ! arguments
   implicit none
   class(streamcontrol_type) :: this
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   character(len=CL)  :: stream_fldFileName_ch4finundated = ' '
   character(len=CL)  :: ch4finundatedmapalgo = 'bilinear'
   character(len=*), parameter :: namelist_name = 'ch4finundated'    ! MUST agree with name in namelist and read
   character(len=*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(len=*), parameter :: subName = "('ch4finundated::ReadNML')"
   character(len=*), parameter :: F00 = "('(ch4finundated_readnml) ',4a)"
   !-----------------------------------------------------------------------

   namelist /ch4finundated/        &               ! MUST agree with namelist_name above
        ch4finundatedmapalgo,  stream_fldFileName_ch4finundated

   ! Default values for namelist

   ! Read ch4finundateddyn_nml namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, namelist_name, status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=ch4finundated,iostat=nml_error)   ! MUST agree with namelist_name above
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_fldFileName_ch4finundated, mpicom)
   call shr_mpi_bcast(ch4finundatedmapalgo            , mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) namelist_name, ' stream settings:'
      write(iulog,*) '  stream_fldFileName_ch4finundated = ',stream_fldFileName_ch4finundated
      write(iulog,*) '  ch4finundatedmapalgo             = ',ch4finundatedmapalgo
      write(iulog,*) ' '
   endif
   this%stream_fldFileName_ch4finundated = stream_fldFileName_ch4finundated
   this%ch4finundatedmapalgo             = ch4finundatedmapalgo
   if (      finundation_mtd == finundation_mtd_ZWT_inversion  )then
       this%fldList = "ZWT0:F0:P3"
   else if ( finundation_mtd == finundation_mtd_TWS_inversion  )then
       this%fldList = "FWS_TWS_A:FWS_TWS_B"
   else
       call endrun(msg=' ERROR do NOT know what list of variables to read for this finundation_mtd type'// &
                   errMsg(sourcefile, __LINE__))
   end if

 end subroutine ReadNML

end module ch4FInundatedStreamType
