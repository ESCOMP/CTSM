module Wateratm2lndBulkType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water atm2lnd variables that just apply to bulk
  ! water. Note that this type extends the base wateratm2lnd_type, so the full
  ! wateratm2lndbulk_type contains the union of the fields defined here and the fields
  ! defined in wateratm2lnd_type.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use PatchType      , only : patch
  use clm_varctl     , only : iulog, use_fates, use_cn, use_cndv
  use clm_varcon     , only : spval
  use WaterAtm2lndType , only : wateratm2lnd_type
  use WaterInfoBaseType, only : water_info_base_type
  use WaterTracerContainerType, only : water_tracer_container_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, extends(wateratm2lnd_type), public :: wateratm2lndbulk_type

     real(r8), pointer :: volrmch_grc                   (:)   ! rof volr main channel (m3)
     real(r8), pointer :: volr_grc                      (:)   ! rof volr total volume (m3)
     real(r8), pointer :: forc_rh_grc                   (:)   ! atmospheric relative humidity (%)
     real(r8) , pointer :: prec365_col                  (:)   ! col 365-day running mean of tot. precipitation (see comment in UpdateAccVars regarding why this is col-level despite other prec accumulators being patch-level)
     real(r8) , pointer :: prec60_patch                 (:)   ! patch 60-day running mean of tot. precipitation (mm/s)
     real(r8) , pointer :: prec10_patch                 (:)   ! patch 10-day running mean of tot. precipitation (mm/s)
     real(r8) , pointer :: rh30_patch                   (:)   ! patch 30-day running mean of relative humidity
     real(r8) , pointer :: prec24_patch                 (:)   ! patch 24-hour running mean of tot. precipitation (mm/s)
     real(r8) , pointer :: rh24_patch                   (:)   ! patch 24-hour running mean of relative humidity

   contains

     procedure, public  :: InitBulk
     procedure, public  :: InitForTesting  ! Should only be used in unit tests
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars
     procedure, public  :: Clean
     procedure, private :: InitBulkAllocate
     procedure, private :: InitBulkHistory
     procedure, private :: InitBulkCold

  end type wateratm2lndbulk_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine InitBulk(this, bounds, info, vars)

    class(wateratm2lndbulk_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds
    class(water_info_base_type), intent(in), target :: info
    type(water_tracer_container_type), intent(inout) :: vars


    call this%Init(bounds, info, vars)

    call this%InitBulkAllocate(bounds)

    call this%InitBulkHistory(bounds)

    call this%InitBulkCold(bounds)

  end subroutine InitBulk

  !------------------------------------------------------------------------
  subroutine InitForTesting(this, bounds, info)
    ! Init routine only for unit testing
    class(wateratm2lndbulk_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds
    class(water_info_base_type), intent(in), target :: info

    type(water_tracer_container_type) :: vars

    ! In unit tests, we don't care about the vars structure, so don't force tests to
    ! create it
    call vars%init()
    call this%InitBulk(bounds, info, vars)
  end subroutine InitForTesting

  !------------------------------------------------------------------------
  subroutine InitBulkAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(wateratm2lndbulk_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%volr_grc                      (begg:endg))        ; this%volr_grc    (:)   = ival
    allocate(this%volrmch_grc                   (begg:endg))        ; this%volrmch_grc (:)   = ival
    allocate(this%forc_rh_grc                   (begg:endg))        ; this%forc_rh_grc (:)   = ival
    allocate(this%prec365_col                   (begc:endc))        ; this%prec365_col (:)   = nan
    allocate(this%prec60_patch                  (begp:endp))        ; this%prec60_patch(:)   = nan
    allocate(this%prec10_patch                  (begp:endp))        ; this%prec10_patch(:)   = nan
    allocate(this%rh30_patch                    (begp:endp))        ; this%rh30_patch  (:)   = nan
    if (use_fates) then
       allocate(this%prec24_patch               (begp:endp))        ; this%prec24_patch(:)   = nan
       allocate(this%rh24_patch                 (begp:endp))        ; this%rh24_patch  (:)   = nan
    end if


  end subroutine InitBulkAllocate

  !------------------------------------------------------------------------
  subroutine InitBulkHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod    , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(wateratm2lndbulk_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begg = bounds%begg; endg= bounds%endg


    this%volr_grc(begg:endg) = spval
    call hist_addfld1d (fname=this%info%fname('VOLR'),  units='m3',  &
         avgflag='A', long_name=this%info%lname('river channel total water storage'), &
         ptr_lnd=this%volr_grc)

    this%volrmch_grc(begg:endg) = spval
    call hist_addfld1d (fname=this%info%fname('VOLRMCH'),  units='m3',  &
         avgflag='A', long_name=this%info%lname('river channel main channel water storage'), &
         ptr_lnd=this%volrmch_grc)

    this%forc_rh_grc(begg:endg) = spval
    call hist_addfld1d (fname=this%info%fname('RH'), units='%',  &
         avgflag='A', long_name=this%info%lname('atmospheric relative humidity'), &
         ptr_gcell=this%forc_rh_grc, default='inactive')

    if (use_cn) then
       this%rh30_patch(begp:endp) = spval
       call hist_addfld1d (fname=this%info%fname('RH30'), units='%',  &
            avgflag='A', long_name=this%info%lname('30-day running mean of relative humidity'), &
            ptr_patch=this%rh30_patch, default='inactive')

       this%prec10_patch(begp:endp) = spval
       call hist_addfld1d (fname=this%info%fname('PREC10'), units='MM H2O/S',  &
            avgflag='A', long_name=this%info%lname('10-day running mean of PREC'), &
            ptr_patch=this%prec10_patch, default='inactive')

       this%prec60_patch(begp:endp) = spval
       call hist_addfld1d (fname=this%info%fname('PREC60'), units='MM H2O/S',  &
            avgflag='A', long_name=this%info%lname('60-day running mean of PREC'), &
            ptr_patch=this%prec60_patch, default='inactive')
    end if

  end subroutine InitBulkHistory

  !-----------------------------------------------------------------------
  subroutine InitBulkCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions
    !
    ! !ARGUMENTS:
    class(wateratm2lndbulk_type), intent(inout) :: this
    type(bounds_type)     , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    ! Nothing to do for now

  end subroutine InitBulkCold

  !------------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES
    use clm_varcon  , only : spval
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(wateratm2lndbulk_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    !---------------------------------------------------------------------

    if (use_cn) then
       call init_accum_field (name='PREC10', units='MM H2O/S', &
            desc='10-day running mean of total precipitation', accum_type='runmean', accum_period=-10, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='PREC60', units='MM H2O/S', &
            desc='60-day running mean of total precipitation', accum_type='runmean', accum_period=-60, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='RH30', units='%', &
            desc='30-day running mean of relative humidity', accum_type='runmean', accum_period=-30, &
            subgrid_type='pft', numlev=1, init_value=100._r8)
    end if

    if (use_cndv) then
       ! The following is a running mean with the accumulation period is set to -365 for a 365-day running mean.
       call init_accum_field (name='PREC365', units='MM H2O/S', &
            desc='365-day running mean of total precipitation', accum_type='runmean', accum_period=-365, &
            subgrid_type='column', numlev=1, init_value=0._r8)
    end if

    if ( use_fates ) then
       call init_accum_field (name='PREC24', units='m', &
            desc='24hr sum of precipitation', accum_type='runmean', accum_period=-1, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       ! Fudge - this neds to be initialized from the restat file eventually.
       call init_accum_field (name='RH24', units='m', &
            desc='24hr average of RH', accum_type='runmean', accum_period=-1, &
            subgrid_type='pft', numlev=1, init_value=100._r8)
    end if

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(wateratm2lndbulk_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: begc, endc
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    real(r8), pointer :: rbufslc(:)  ! temporary
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="InitAccVars allocation error for rbufslp"//&
            errMsg(sourcefile, __LINE__))
    endif
    ! Allocate needed dynamic memory for single level col field
    allocate(rbufslc(begc:endc), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="InitAccVars allocation error for rbufslc"//&
            errMsg(sourcefile, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    if (use_cn) then
       call extract_accum_field ('PREC10', rbufslp, nstep)
       this%prec10_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('PREC60', rbufslp, nstep)
       this%prec60_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('RH30', rbufslp, nstep)
       this%rh30_patch(begp:endp) = rbufslp(begp:endp)
    end if

    if (use_cndv) then
       call extract_accum_field ('PREC365' , rbufslc, nstep)
       this%prec365_col(begc:endc) = rbufslc(begc:endc)
    end if

    if (use_fates) then
       call extract_accum_field ('PREC24', rbufslp, nstep)
       this%prec24_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('RH24', rbufslp, nstep)
       this%rh24_patch(begp:endp) = rbufslp(begp:endp)
    end if

    deallocate(rbufslp)
    deallocate(rbufslc)

  end subroutine InitAccVars

  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(wateratm2lndbulk_type), intent(in) :: this
    type(bounds_type)      , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p                     ! indices
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: rbufslp(:)      ! temporary single level - patch level
    real(r8), pointer :: rbufslc(:)      ! temporary single level - column level
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'UpdateAccVars allocation error for rbufslp'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif
    ! Allocate needed dynamic memory for single level col field
    allocate(rbufslc(begc:endc), stat=ier)
    if (ier/=0) then
       write(iulog,*)'UpdateAccVars allocation error for rbufslc'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif


    ! Precipitation accumulators
    !
    ! For CNDV, we use a column-level accumulator. We cannot use a patch-level
    ! accumulator for CNDV because this is used for establishment, so must be available
    ! for inactive patches. In principle, we could/should switch to column-level for the
    ! other precip accumulators, too; we'd just need to be careful about backwards
    ! compatibility with old restart files.

    do p = begp,endp
       c = patch%column(p)
       rbufslp(p) = this%forc_rain_downscaled_col(c) + this%forc_snow_downscaled_col(c)
       rbufslc(c) = this%forc_rain_downscaled_col(c) + this%forc_snow_downscaled_col(c)
    end do

    if (use_cn) then
       ! Accumulate and extract PREC60 (accumulates total precipitation as 60-day running mean)
       call update_accum_field  ('PREC60', rbufslp, nstep)
       call extract_accum_field ('PREC60', this%prec60_patch, nstep)

       ! Accumulate and extract PREC10 (accumulates total precipitation as 10-day running mean)
       call update_accum_field  ('PREC10', rbufslp, nstep)
       call extract_accum_field ('PREC10', this%prec10_patch, nstep)
    end if

    if (use_cndv) then
       ! Accumulate and extract PREC365 (accumulates total precipitation as 365-day running mean)
       ! See above comment regarding why this is at the column-level despite other prec
       ! accumulators being at the patch level.
       call update_accum_field  ('PREC365', rbufslc, nstep)
       call extract_accum_field ('PREC365', this%prec365_col, nstep)

    end if

    if (use_fates) then
       call update_accum_field  ('PREC24', rbufslp, nstep)
       call extract_accum_field ('PREC24', this%prec24_patch, nstep)

       do p = bounds%begp,bounds%endp
          g = patch%gridcell(p)
          rbufslp(p) = this%forc_rh_grc(g)
       end do
       call update_accum_field  ('RH24', rbufslp, nstep)
       call extract_accum_field ('RH24', this%rh24_patch, nstep)
    end if

    if (use_cn) then
       do p = begp,endp
          g = patch%gridcell(p)
          rbufslp(p) = this%forc_rh_grc(g)
       end do
       ! Accumulate and extract RH30 (accumulates RH as 30-day running mean)
       call update_accum_field  ('RH30', rbufslp, nstep)
       call extract_accum_field ('RH30', this%rh30_patch, nstep)
    endif

    deallocate(rbufslp)
    deallocate(rbufslc)

  end subroutine UpdateAccVars

  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Finalize this instance
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(wateratm2lndbulk_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    ! atm->lnd
    deallocate(this%forc_rh_grc)

    ! atm->lnd not downscaled
    deallocate(this%forc_q_not_downscaled_grc)
    deallocate(this%forc_rain_not_downscaled_grc)
    deallocate(this%forc_snow_not_downscaled_grc)

    ! atm->lnd downscaled
    deallocate(this%forc_q_downscaled_col)
    deallocate(this%forc_rain_downscaled_col)
    deallocate(this%forc_snow_downscaled_col)

    ! rof->lnd
    deallocate(this%forc_flood_grc)
    deallocate(this%volr_grc)
    deallocate(this%volrmch_grc)

    ! anomaly forcing
    deallocate(this%prec365_col)
    deallocate(this%prec60_patch)
    deallocate(this%prec10_patch)
    if (use_fates) then
       deallocate(this%prec24_patch)
       deallocate(this%rh24_patch)
    end if

  end subroutine Clean

end module Wateratm2lndBulkType
