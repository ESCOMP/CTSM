#include <misc.h>
#include <preproc.h>

module accFldsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: accFldsMod
!
! !DESCRIPTION:
! This module contains subroutines that initialize, update and extract
! the user-specified fields over user-defined intervals. Each interval
! and accumulation type is unique to each field processed.
! Subroutine [initAccumFlds] defines the fields to be processed
! and the type of accumulation. Subroutine [updateAccumFlds] does
! the actual accumulation for a given field. Fields are accumulated
! by calls to subroutine [update_accum_field]. To accumulate a field,
! it must first be defined in subroutine [initAccumFlds] and then
! accumulated by calls to [updateAccumFlds].
! Four types of accumulations are possible:
!   o average over time interval
!   o running mean over time interval
!   o running accumulation over time interval
! Time average fields are only valid at the end of the averaging interval.
! Running means are valid once the length of the simulation exceeds the
! averaging interval. Accumulated fields are continuously accumulated.
! The trigger value "-99999." resets the accumulation to zero.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initAccFlds     ! Initialization accumulator fields
  public :: initAccClmtype  ! Initialize clmtype variables obtained from accum fields
  public :: updateAccFlds   ! Update accumulator fields
!
! !REVISION HISTORY:
! Created by M. Vertenstein 03/2003
!
!EOP

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initAccFlds()
!
! !INTERFACE:
  subroutine initAccFlds()
!
! !DESCRIPTION:
! Initializes accumulator and sets up array of accumulated fields
!
! !USES:
    use accumulMod   , only : init_accum_field, print_accum_fields
    use clm_time_manager , only : get_step_size
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use nanMod       , only : bigint
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY::
! Created by M. Vertenstein 03/2003
!
!EOP
!
! LOCAL VARIABLES:
!
    integer :: dtime                     !time step size
    integer, parameter :: not_used = bigint
!------------------------------------------------------------------------

    ! Hourly average of 2m temperature.

    dtime = get_step_size()
    call init_accum_field(name='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

#if (defined DGVM)
    ! 30-day average of 2m temperature.

    call init_accum_field (name='TDA', units='K', &
         desc='30-day average of 2-m temperature', &
         accum_type='timeavg', accum_period=-30, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following are running means.
    ! The accumulation period is set to 10 days for a 10-day running mean.

    call init_accum_field (name='T10', units='K', &
         desc='10-day running mean of 2-m temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=SHR_CONST_TKFRZ+20._r8)

    call init_accum_field (name='FNPSN10', units='UMOL/M2S', &
         desc='10-day running mean net cpy photosynth', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='PREC365', units='MM H2O/S', &
         desc='365-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-365, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following are accumulated fields.
    ! These types of fields are accumulated until a trigger value resets
    ! the accumulation to zero (see subroutine update_accum_field).
    ! Hence, [accper] is not valid.

    call init_accum_field (name='AGDD0', units='K', &
         desc='growing degree-days base 0C', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='AGDD5', units='K', &
         desc='growing degree-days base -5C', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='AGDDTW', units='K', &
         desc='growing degree-days base twmax', &
         accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='AGDD', units='K', &
         desc='growing degree-days base 5C', &
         accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)
#endif

    ! Print output of accumulated fields

    call print_accum_fields()

  end subroutine initAccFlds

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateAccFlds
!
! !INTERFACE:
  subroutine updateAccFlds()
!
! !DESCRIPTION:
! Update and/or extract accumulated fields
!
! !USES:
    use clmtype
    use clm_atmlnd   , only : clm_a2l
    use decompMod    , only : get_proc_bounds
    use clm_varcon   , only : spval
    use pftvarcon    , only : pftpar
    use shr_const_mod, only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod   , only : update_accum_field, extract_accum_field
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by M. Vertenstein 03/2003
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: itype(:)            ! pft vegetation
    integer , pointer :: pgridcell(:)        ! index into gridcell level quantities
    real(r8), pointer :: forc_t(:)           ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_rain(:)        ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)        ! snow rate [mm/s]
    real(r8), pointer :: frmf(:)             ! leaf maintenance respiration  (umol CO2 /m**2 /s)
    real(r8), pointer :: fpsn(:)             ! photosynthesis (umol CO2 /m**2 /s)
    real(r8), pointer :: t_ref2m(:)          ! 2 m height surface air temperature (Kelvin)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: t_ref2m_min(:)      ! daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max(:)      ! daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_inst(:) ! instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst(:) ! instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t10(:)              ! 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: t_mo(:)             ! 30-day average temperature (Kelvin)
    real(r8), pointer :: t_mo_min(:)         ! annual min of t_mo (Kelvin)
    real(r8), pointer :: fnpsn10(:)          ! 10-day running mean net photosynthesis
    real(r8), pointer :: prec365(:)          ! 365-day running mean of tot. precipitation
    real(r8), pointer :: agdd0(:)            ! accumulated growing degree days above 0 deg C
    real(r8), pointer :: agdd5(:)            ! accumulated growing degree days above -5
    real(r8), pointer :: agddtw(:)           ! accumulated growing degree days above twmax
    real(r8), pointer :: agdd(:)             ! accumulated growing degree days above 5
!
!EOP
!
! OTHER LOCAL VARIABLES:
    integer :: g,l,c,p                   ! indices
    integer :: itypveg                   ! vegetation type
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: year                      ! year (0, ...) for nstep
    integer :: month                     ! month (1, ..., 12) for nstep
    integer :: day                       ! day of month (1, ..., 31) for nstep
    integer :: secs                      ! seconds into current date for nstep
    logical :: end_cd                    ! temporary for is_end_curr_day() value
    integer :: ier                       ! error status
    integer :: begp, endp                !  per-proc beginning and ending pft indices
    integer :: begc, endc                !  per-proc beginning and ending column indices
    integer :: begl, endl                !  per-proc beginning and ending landunit indices
    integer :: begg, endg                !  per-proc gridcell ending gridcell indices
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_t => clm_a2l%forc_t
    forc_rain => clm_a2l%forc_rain
    forc_snow => clm_a2l%forc_snow

    ! Assign local pointers to derived subtypes components (pft-level)

    itype            => clm3%g%l%c%p%itype
    pgridcell        => clm3%g%l%c%p%gridcell
    t_ref2m          => clm3%g%l%c%p%pes%t_ref2m
    fpsn             => clm3%g%l%c%p%pcf%fpsn
    frmf             => clm3%g%l%c%p%pcf%frmf
    t_ref2m_max_inst => clm3%g%l%c%p%pes%t_ref2m_max_inst
    t_ref2m_min_inst => clm3%g%l%c%p%pes%t_ref2m_min_inst
    t_ref2m_max      => clm3%g%l%c%p%pes%t_ref2m_max
    t_ref2m_min      => clm3%g%l%c%p%pes%t_ref2m_min
    t_mo             => clm3%g%l%c%p%pdgvs%t_mo
    t_mo_min         => clm3%g%l%c%p%pdgvs%t_mo_min
    t10              => clm3%g%l%c%p%pdgvs%t10
    fnpsn10          => clm3%g%l%c%p%pdgvs%fnpsn10
    prec365          => clm3%g%l%c%p%pdgvs%prec365
    agdd0            => clm3%g%l%c%p%pdgvs%agdd0
    agdd5            => clm3%g%l%c%p%pdgvs%agdd5
    agddtw           => clm3%g%l%c%p%pdgvs%agddtw
    agdd             => clm3%g%l%c%p%pdgvs%agdd

    ! Determine calendar information

    dtime = get_step_size()
    nstep = get_nstep()
    call get_curr_date (year, month, day, secs)

    ! Don't do any accumulation if nstep is zero
    ! (only applies to coupled or cam mode)

    if (nstep == 0) return

    ! NOTE: currently only single level pft fields are used below
    ! Variables are declared above that should make it easy to incorporate
    ! multi-level or single-level fields of any subgrid type

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(6,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun
    endif

    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV', t_ref2m, nstep)
    call extract_accum_field ('TREFAV', rbufslp, nstep)
    end_cd = is_end_curr_day()
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          t_ref2m_max_inst(p) = max(rbufslp(p), t_ref2m_max_inst(p))
          t_ref2m_min_inst(p) = min(rbufslp(p), t_ref2m_min_inst(p))
       endif
       if (end_cd) then
          t_ref2m_max(p) = t_ref2m_max_inst(p)
          t_ref2m_min(p) = t_ref2m_min_inst(p)
          t_ref2m_max_inst(p) = -spval
          t_ref2m_min_inst(p) =  spval
       else if (secs == int(dtime)) then
          t_ref2m_max(p) = spval
          t_ref2m_min(p) = spval
       endif
    end do

#if (defined DGVM)
    ! Accumulate and extract TDA
    ! (accumulates TBOT as 30-day average)
    ! Also determine t_mo_min

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_t(g)
    end do
    call update_accum_field  ('TDA', rbufslp, nstep)
    call extract_accum_field ('TDA', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       t_mo(p) = rbufslp(p)
       t_mo_min(p) = min(t_mo_min(p), rbufslp(p))
    end do

    ! Accumulate and extract T10
    !(acumulates TSA as 10-day running mean)

    call update_accum_field  ('T10', t_ref2m, nstep)
    call extract_accum_field ('T10', t10, nstep)

    ! Accumulate and extract FNPSN10
    !(accumulates fpsn-frmf as 10-day running mean)

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       rbufslp(p) = fpsn(p) - frmf(p)
    end do
    call update_accum_field  ('FNPSN10', rbufslp, nstep)
    call extract_accum_field ('FNPSN10', fnpsn10, nstep)

    ! Accumulate and extract PREC365
    ! (accumulates total precipitation as 365-day running mean)

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       g = pgridcell(p)
       rbufslp(p) = forc_rain(g) + forc_snow(g)
    end do
    call update_accum_field  ('PREC365', rbufslp, nstep)
    call extract_accum_field ('PREC365', prec365, nstep)

    ! Accumulate growing degree days based on 10-day running mean temperature.
    ! Accumulate GDD above 0C and -5C using extracted t10 from accumulated variable.
    ! The trigger to reset the accumulated values to zero is -99999.
    ! agddtw is currently reset at the end of each year in subr. lpj

    ! Accumulate and extract AGDDO

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       rbufslp(p) = (t10(p) - SHR_CONST_TKFRZ) * dtime / SHR_CONST_CDAY
       if (rbufslp(p) < 0._r8) rbufslp(p) = -99999._r8
    end do
    call update_accum_field  ('AGDD0', rbufslp, nstep)
    call extract_accum_field ('AGDD0', agdd0, nstep)

    ! Accumulate and extract AGDD5

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       rbufslp(p) = (t10(p) - (SHR_CONST_TKFRZ - 5.0_r8))*dtime / SHR_CONST_CDAY
       if (rbufslp(p) < 0._r8) rbufslp(p) = -99999._r8
    end do
    call update_accum_field  ('AGDD5', rbufslp, nstep)
    call extract_accum_field ('AGDD5', agdd5, nstep)

    ! Accumulate and extract AGDDTW

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       itypveg = itype(p)
       rbufslp(p) = max(0.0_r8, (t10(p) - (SHR_CONST_TKFRZ+pftpar(itypveg,31))) &
            * dtime/SHR_CONST_CDAY)
    end do
    call update_accum_field  ('AGDDTW', rbufslp, nstep)
    call extract_accum_field ('AGDDTW', agddtw, nstep)

    ! Accumulate and extract AGDD

!dir$ concurrent
!cdir nodep
    do p = begp,endp
       rbufslp(p) = max(0.0_r8, (t_ref2m(p) - (SHR_CONST_TKFRZ + 5.0_r8)) &
            * dtime/SHR_CONST_CDAY)
    end do
    call update_accum_field  ('AGDD', rbufslp, nstep)
    call extract_accum_field ('AGDD', agdd, nstep)
#endif

    ! Deallocate dynamic memory

    deallocate(rbufslp)

  end subroutine updateAccFlds

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initAccClmtype
!
! !INTERFACE:
  subroutine initAccClmtype
!
! !DESCRIPTION:
! Initialize clmtype variables that are associated with
! time accumulated fields. This routine is called in an initial run
! at nstep=0 for cam and csm mode and at nstep=1 for offline mode.
! This routine is also always called for a restart run and
! therefore must be called after the restart file is read in
! and the accumulated fields are obtained.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use decompMod   , only : get_proc_bounds, get_proc_global
    use accumulMod  , only : extract_accum_field
    use clm_time_manager, only : get_nstep
    use clm_varctl  , only : nsrest
    use clm_varcon  , only : spval
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: t_ref2m_min(:)      ! daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max(:)      ! daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_inst(:) ! instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst(:) ! instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t10(:)              ! 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: t_mo(:)             ! 30-day average temperature (Kelvin)
    real(r8), pointer :: fnpsn10(:)          ! 10-day running mean net photosynthesis
    real(r8), pointer :: prec365(:)          ! 365-day running mean of tot. precipitation
    real(r8), pointer :: agdd0(:)            ! accumulated growing degree days above 0 deg C
    real(r8), pointer :: agdd5(:)            ! accumulated growing degree days above -5
    real(r8), pointer :: agddtw(:)           ! accumulated growing degree days above twmax
    real(r8), pointer :: agdd(:)             ! accumulated growing degree days above 5
!
! !LOCAL VARIABLES:
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: nstep        ! time step
    integer :: ier          ! error status
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    real(r8), pointer :: rbufslp(:)  ! temporary
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (pft-level)

    t_ref2m_max_inst => clm3%g%l%c%p%pes%t_ref2m_max_inst
    t_ref2m_min_inst => clm3%g%l%c%p%pes%t_ref2m_min_inst
    t_ref2m_max      => clm3%g%l%c%p%pes%t_ref2m_max
    t_ref2m_min      => clm3%g%l%c%p%pes%t_ref2m_min
    t10              => clm3%g%l%c%p%pdgvs%t10
    t_mo             => clm3%g%l%c%p%pdgvs%t_mo
    fnpsn10          => clm3%g%l%c%p%pdgvs%fnpsn10
    prec365          => clm3%g%l%c%p%pdgvs%prec365
    agdd0            => clm3%g%l%c%p%pdgvs%agdd0
    agdd5            => clm3%g%l%c%p%pdgvs%agdd5
    agddtw           => clm3%g%l%c%p%pdgvs%agddtw
    agdd             => clm3%g%l%c%p%pdgvs%agdd

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Determine time step

    nstep = get_nstep()

    ! Initialize 2m ref temperature max and min values

    if (nsrest == 0) then
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          t_ref2m_max(p) = spval
          t_ref2m_min(p) = spval
          t_ref2m_max_inst(p) = -spval
          t_ref2m_min_inst(p) =  spval
       end do
    end if

#if (defined DGVM)

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(6,*)'update_accum_hist allocation error for rbufslp'
       call endrun
    endif

    ! Initialize clmtype variables that are to be time accumulated

    call extract_accum_field ('T10', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       t10(p) = rbufslp(p)
    end do

    call extract_accum_field ('TDA', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       t_mo(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDD0', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       agdd0(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDD5', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       agdd5(p) = rbufslp(p)
    end do

    call extract_accum_field ('FNPSN10', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       fnpsn10(p) = rbufslp(p)
    end do

    call extract_accum_field ('PREC365', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       prec365(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDDTW', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       agddtw(p) = rbufslp(p)
    end do

    call extract_accum_field ('AGDD', rbufslp, nstep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       agdd(p) = rbufslp(p)
    end do

    deallocate(rbufslp)
#endif

  end subroutine initAccClmtype

end module accFldsMod
