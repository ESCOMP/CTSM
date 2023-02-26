module PrigentRoughnessStreamType
!dmleung modified based on ch4FInundatedStreamType on 17 Nov 2022
#include "shr_assert.h"             ! What is this? In many modules but not dust module

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in the Prigent et al. (1997) roughness length streams file. dmleung 22 Nov 2022
  !
  ! !USES
  use ESMF                         
  use dshr_strdata_mod , only : shr_strdata_type
  use shr_kind_mod     , only : r8 => shr_kind_r8, CL => shr_kind_cl
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use spmdMod          , only : mpicom, masterproc
  use clm_varctl       , only : iulog
  use abortutils       , only : endrun
  use decompMod        , only : bounds_type
  !use ch4varcon        , only : finundation_mtd   ! dmleung commented out. need to replace anything here?

  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: prigentroughnessstream_type
     ! dmleung commented the next five lines. Does it need any replacement here?
     real(r8), pointer, private :: prigent_rghn  (:)         ! Prigent et al. (1997) roughness length (m)
     !real(r8), pointer, private :: zwt0_gdc  (:)         ! col coefficient for determining finundated (m)
     !real(r8), pointer, private :: f0_gdc    (:)         ! col maximum inundated fraction for a gridcell (for methane code)
     !real(r8), pointer, private :: p3_gdc    (:)         ! col coefficient for determining finundated (m)
     !real(r8), pointer, private :: fws_slope_gdc     (:) ! col slope in fws = slope * tws + intercept (A coefficient)
     !real(r8), pointer, private :: fws_intercept_gdc (:) ! col slope in fws = slope * tws + intercept (B coefficient)
  contains

      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: Init            ! Initialize and read data in
      !procedure, public :: CalcFinundated  ! Calculate finundated based on input streams
      procedure, public :: CalcDragPartition ! Calculate drag partitioning based on input streams
      procedure, public :: UseStreams      ! If streams will be used -- dmleung set default as yes

      ! !PRIVATE MEMBER FUNCTIONS:
      procedure, private :: InitAllocate   ! Allocate data

  end type prigentroughnessstream_type

  ! ! PRIVATE DATA:
  type, private :: streamcontrol_type
     character(len=CL)  :: stream_fldFileName_prigentroughness   ! data Filename
     character(len=CL)  :: stream_meshfile_prigentroughness      ! mesh Filename
     character(len=CL)  :: prigentroughnessmapalgo               ! map algo
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
   ! Initialize the prigent roughness stream object
   !
   ! Uses:
   use spmdMod          , only : iam                !what is this and do I need this
   !use ch4varcon        , only : finundation_mtd_h2osfc
   !use ch4varcon        , only : finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
   use lnd_comp_shr     , only : mesh, model_clock
   use dshr_strdata_mod , only : shr_strdata_init_from_inline, shr_strdata_print
   use dshr_strdata_mod , only : shr_strdata_advance
   use dshr_methods_mod , only : dshr_fldbun_getfldptr                !what is this and do I need this (read ncdata?)
   !
   ! arguments
   implicit none
   class(prigentroughnessstream_type) :: this
   type(bounds_type), intent(in)   :: bounds
   character(len=*),  intent(in)   :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer                        :: ig, g, n           ! Indices
   integer                        :: year               ! year (0, ...) for nstep+1
   integer                        :: mon                ! month (1, ..., 12) for nstep+1
   integer                        :: day                ! day of month (1, ..., 31) for nstep+1
   integer                        :: sec                ! seconds into current date for nstep+1
   integer                        :: mcdate             ! Current model date (yyyymmdd)
   type(shr_strdata_type)         :: sdat_rghn      ! input data stream
   character(len=16), allocatable :: stream_varnames(:) ! array of stream field names
   integer                        :: rc                 ! error code
   real(r8), pointer              :: dataptr1d(:)       ! temporary pointer
   character(len=*), parameter    :: stream_name = 'prigentroughness'
   !-----------------------------------------------------------------------

   !if ( finundation_mtd /= finundation_mtd_h2osfc )then     ! how should I change this? comment out for now
      call this%InitAllocate( bounds )
      call control%ReadNML( bounds, NLFileName )

      if ( this%useStreams() )then     ! is this a namelist input and is it set in namelist default

         !if (finundation_mtd == finundation_mtd_ZWT_inversion  )then
         !   allocate(stream_varnames(3))
         !   stream_varnames = (/"ZWT0","F0  ","P3  "/)
         !else if ( finundation_mtd == finundation_mtd_TWS_inversion  )then
            allocate(stream_varnames(1))
            stream_varnames = (/"Z0a"/)  ! varname in cdf5_Z0a_Prigent-Globe-025x025-09262022.nc, in centimeter
         !else
         !   call endrun(msg=' ERROR do NOT know what list of variables to read for this finundation_mtd type'// &
         !        errMsg(sourcefile, __LINE__))
         !end if

         if (masterproc) then
            write(iulog,*) '  stream_varnames                  = ',stream_varnames
      end if

      ! Initialize the cdeps data type sdat_rghn
      call shr_strdata_init_from_inline(sdat_rghn,                                  &     ! what is this function and where does it come from?
              my_task             = iam,                                                &
              logunit             = iulog,                                              &
              compname            = 'LND',                                              &
              model_clock         = model_clock,                                        &
              model_mesh          = mesh,                                               &
              stream_meshfile     = control%stream_meshfile_prigentroughness,              &
              stream_lev_dimname  = 'null',                                             &
              stream_mapalgo      = control%prigentroughnessmapalgo,                       &
              stream_filenames    = (/trim(control%stream_fldFileName_prigentroughness)/), &
              stream_fldlistFile  = stream_varnames,                                    &
              stream_fldListModel = stream_varnames,                                    &
              stream_yearFirst    = 1997,                                               &
              stream_yearLast     = 1997,                                               &
              stream_yearAlign    = 1,                                                  &
              stream_offset       = 0,                                                  &
              stream_taxmode      = 'extend',                                           &
              stream_dtlimit      = 1.0e30_r8,                                          &
              stream_tintalgo     = 'linear',                                           &
              stream_name         = 'Prigent roughness',                                 &
              rc                  = rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! Explicitly set current date to a hardcoded constant value. Otherwise
         ! using the real date can cause roundoff differences that are
         ! detrected as issues with exact restart.  EBK M05/20/2017
         ! call get_curr_date(year, mon, day, sec)
         year = 1997
         mon  = 12
         day  = 31
         sec  = 0
         mcdate = year*10000 + mon*100 + day

         call shr_strdata_advance(sdat_rghn, ymd=mcdate, tod=sec, logunit=iulog, istr='prigentrghn', rc=rc)   ! what is istr and do I need to change elsewhere because the change of istr here
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
         end if

         ! Get pointer for stream data that is time and spatially interpolate to model time and grid
         do n = 1,size(stream_varnames) 
            call dshr_fldbun_getFldPtr(sdat_rghn%pstrm(1)%fldbun_model, stream_varnames(n), fldptr1=dataptr1d, rc=rc)   ! is the ncdf datasets read here, looping across all available variables in the files? Also, not sure about pstrm and fldbun_model
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
            end if
            if (trim(stream_varnames(n)) == 'Z0a') then
               ig = 0
               do g = bounds%begg,bounds%endg
                  ig = ig+1                   ! not sure why +1 is needed but it's okay
                  this%prigent_rghn(g) = dataptr1d(ig)
               end do

            end if

         end do
      end if
   !end if         !comment out for now

  end subroutine Init

  !==============================================================================
  logical function UseStreams(this)
    !
    ! !DESCRIPTION:
    ! Return true if
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    class(prigentroughnessstream_type) :: this
    !
    ! !LOCAL VARIABLES:
    if ( trim(control%stream_fldFileName_prigentroughness) == '' )then
       UseStreams = .false.  ! dmleung: this won't happen and UseStreams will always be true
    else
       UseStreams = .true.
    end if
  end function UseStreams

  !==============================================================================
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    !use ch4varcon     , only: finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
    !
    ! !ARGUMENTS:
    implicit none
    class(prigentroughnessstream_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !integer  :: begc, endc
    integer  :: begg, endg
    !---------------------------------------------------------------------

    !begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    !if( finundation_mtd == finundation_mtd_ZWT_inversion )then
       allocate(this%prigent_rghn     (begg:endg))            ;  this%prigent_rghn     (:)   = nan
    !   allocate(this%zwt0_gdc         (begg:endg))            ;  this%zwt0_gdc         (:)   = nan
    !   allocate(this%f0_gdc           (begg:endg))            ;  this%f0_gdc           (:)   = nan
    !   allocate(this%p3_gdc           (begg:endg))            ;  this%p3_gdc           (:)   = nan
    !else if( finundation_mtd == finundation_mtd_TWS_inversion )then
    !   allocate(this%fws_slope_gdc    (begg:endg))            ;  this%fws_slope_gdc    (:)   = nan
    !   allocate(this%fws_intercept_gdc(begg:endg))            ;  this%fws_intercept_gdc(:)   = nan
    !end if

  end subroutine InitAllocate

  !==============================================================================
  !subroutine CalcFinundated(this, bounds, num_soilc, filter_soilc, soilhydrology_inst, &
  !                          waterdiagnosticbulk_inst, qflx_surf_lag_col, finundated )
  !subroutine CalcDragPartition(this, bounds, num_nolakep, filter_nolakep, dpfct_rock)
  subroutine CalcDragPartition(this, bounds, dpfct_rock)
    !
    ! !DESCRIPTION:
    ! Commented below by dmleung 31 Dec 2022
    ! Calculate the drag partition effect of friction velocity due to surface roughness following 
    ! Leung et al. (2022).  This module is used in the dust emission module DUSTMod.F90 for 
    ! calculating drag partitioning. The drag partition equation comes from Marticorena and 
    ! Bergametti (1995) with constants modified by Darmenova et al. (2009). Here it is assumed 
    ! that this equation is used only over arid/desertic regions, such that Catherine Prigent's
    ! roughness measurements represents mostly rocks. For more vegetated areas, the vegetation
    ! roughness and drag partitioning are calculated in the DustEmission subroutine. This 
    ! subroutine is used in the InitCold subroutine of DUSTMod.F90.
    !
    ! !USES:
    !use ColumnType              , only : col
    use PatchType               , only : patch
    !USES dmleung added 31 Dec 2022
    use landunit_varcon         , only : istdlak
    use LandunitType            , only : lun
    !use ch4varcon               , only : finundation_mtd_h2osfc, finundation_mtd_ZWT_inversion
    !use ch4varcon               , only : finundation_mtd_TWS_inversion
    !use clm_varpar              , only : nlevsoi
    !use SoilHydrologyType       , only : soilhydrology_type
    !use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
    !
    ! !ARGUMENTS:
    implicit none
    class(prigentroughnessstream_type)             :: this
    type(bounds_type)              , intent(in)    :: bounds
    !integer                        , intent(in)    :: num_soilc                       ! number of column soil points in column filter
    !integer                        , intent(in)    :: filter_soilc(:)                 ! column filter for soil points
    !integer                        , intent(in)    :: num_nolakep                      !
    !integer                        , intent(in)    :: filter_nolakep(num_nolakep)      !
    !type(soilhydrology_type)       , intent(in)    :: soilhydrology_inst
    !type(waterdiagnosticbulk_type) , intent(in)    :: waterdiagnosticbulk_inst
    !real(r8)                       , intent(in)    :: qflx_surf_lag_col(bounds%begc:) !time-lagged surface runoff (mm H2O /s)
    !real(r8)                       , intent(inout) :: finundated(bounds%begc:)        ! fractional inundated area in soil column (excluding dedicated wetland columns)
    real(r8)                       , intent(inout) :: dpfct_rock(bounds%begp:)      ! [fraction] rock drag partition factor (roughness effect)
    !
    ! !LOCAL VARIABLES:
    !integer  :: g, c, fc    ! Indices
    integer  :: g, p, fp, l    ! Indices
    real(r8) :: z0s         ! smooth roughness length (m)
    !real(r8) :: zwt_actual  ! Total water storage (ZWT) to use either perched or total depending on conditions
    
    ! constants
    real(r8), parameter :: D_p = 130e-6_r8           ! [m] Medium soil particle diameter, assuming a global constant of ~130 um following Leung et al. (2022)
    real(r8), parameter :: X = 10_r8                 ! [m] distance downwind of the roughness element (rock). Assume estimating roughness effect at a distance of 10 m following Leung et al. (2022)
    !---------------------------------------------------------------------

    !SHR_ASSERT_ALL_FL((ubound(qflx_surf_lag_col) == (/bounds%endc/)), sourcefile, __LINE__)    ! not sure what this is
    !SHR_ASSERT_ALL_FL((ubound(finundated)        == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dpfct_rock)        == (/bounds%endp/)), sourcefile, __LINE__) !what's the use of this

    !associate(                                                     &
         !z           =>   col%z                                   & ! Input:  [real(r8) (:,:) ]  layer depth (m) (-nlevsno+1:nlevsoi)
         !zwt         =>   soilhydrology_inst%zwt_col         ,     & ! Input:  [real(r8) (:)   ]  water table depth (m)
         !zwt_perched =>   soilhydrology_inst%zwt_perched_col ,     & ! Input:  [real(r8) (:)   ]  perched water table depth (m)
         !tws         =>   waterdiagnosticbulk_inst%tws_grc   ,     & ! Input:  [real(r8) (:)   ]  total water storage (kg m-2)
         !frac_h2osfc =>   waterdiagnosticbulk_inst%frac_h2osfc_col & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
    !)

    ! Calculate finundated
    !do fc = 1, num_soilc
    !   c = filter_soilc(fc)
    !   g = col%gridcell(c)
    !   select case( finundation_mtd )
    !      case ( finundation_mtd_h2osfc )
    !         finundated(c) = frac_h2osfc(c)
    !      case ( finundation_mtd_ZWT_inversion )
    !         if (this%zwt0_gdc(g) > 0._r8) then
    !            if (zwt_perched(c) < z(c,nlevsoi)-1.e-5_r8 .and. zwt_perched(c) < zwt(c)) then
    !               zwt_actual = zwt_perched(c)
    !            else
    !               zwt_actual = zwt(c)
    !            end if
    !            finundated(c) = this%f0_gdc(g) * exp(-zwt_actual/this%zwt0_gdc(g)) + this%p3_gdc(g)*qflx_surf_lag_col(c)
    !         else
    !            finundated(c) = this%p3_gdc(g)*qflx_surf_lag_col(c)
    !         end if
    !      case ( finundation_mtd_TWS_inversion )
    !         finundated(c) = this%fws_slope_gdc(g) * tws(g) + this%fws_intercept_gdc(g)
    !   end select
    !   finundated(c) = min( 1.0_r8, max( 0.0_r8, finundated(c) ) )
    !end do

    ! dmleung: this loop calculates the drag partition effect (or roughness effect) of rocks. We save the drag partition factor as a patch level quantity.
    z0s = 2_r8 * D_p / 30_r8 ! equation from Frank M. White (2006). Here we assume soil medium size is a global constant, and so is smooth roughness length.
    !do fp = 1,num_nolakep
       !p = filter_nolakep(fp)
    do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       l = patch%landunit(p)
       if (lun%itype(l) /= istdlak) then
          dpfct_rock(p) = 1._r8 - ( log(this%prigent_rghn(g)*0.01_r8/z0s) / log(0.7_r8*(X/z0s)**0.8_r8) ) ! Calculating rock drag partition factor using Marticorena and Bergametti (1995). 0.01 is used to convert Z0a from centimeter to meter.
       end if
    end do

    !end associate

  end subroutine CalcDragPartition

  !==============================================================================
  subroutine ReadNML(this, bounds, NLFilename)
   !
   ! Read the namelist data stream information.
   !
   ! Uses:
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   !use ch4varcon        , only : finundation_mtd_ZWT_inversion, finundation_mtd_TWS_inversion
   !
   ! arguments
   implicit none
   class(streamcontrol_type)     :: this
   type(bounds_type), intent(in) :: bounds
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   character(len=CL)  :: stream_fldFileName_prigentroughness = ' '
   character(len=CL)  :: stream_meshfile_prigentroughness = ' '
   character(len=CL)  :: prigentroughnessmapalgo = 'bilinear'
   character(len=*), parameter :: namelist_name = 'prigentroughness'    ! MUST agree with group name in namelist definition to read. dmleung commented
   character(len=*), parameter :: subName = "('prigentroughness::ReadNML')"
   !-----------------------------------------------------------------------

   namelist /prigentroughness/ &               ! MUST agree with namelist_name above
        prigentroughnessmapalgo,  stream_fldFileName_prigentroughness, stream_meshfile_prigentroughness

   ! Default values for namelist

   ! Read prigentroughness namelist
   if (masterproc) then
      open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, namelist_name, status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=prigentroughness,iostat=nml_error)   ! MUST agree with namelist_name above
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding '//namelist_name//' namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
   endif

   call shr_mpi_bcast(prigentroughnessmapalgo             , mpicom)
   call shr_mpi_bcast(stream_fldFileName_prigentroughness , mpicom)
   call shr_mpi_bcast(stream_meshfile_prigentroughness    , mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) namelist_name, ' stream settings:'
      write(iulog,*) '  stream_fldFileName_prigentroughness = ',stream_fldFileName_prigentroughness
      write(iulog,*) '  stream_meshfile_prigentroughness    = ',stream_meshfile_prigentroughness
      write(iulog,*) '  prigentroughnessmapalgo             = ',prigentroughnessmapalgo
   endif
   this%stream_fldFileName_prigentroughness = stream_fldFileName_prigentroughness
   this%stream_meshfile_prigentroughness    = stream_meshfile_prigentroughness
   this%prigentroughnessmapalgo             = prigentroughnessmapalgo

 end subroutine ReadNML

end module PrigentRoughnessStreamType
