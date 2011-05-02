
module CASAPhenologyMod

#if (defined CASA)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CASAPhenologyMod
!
! !DESCRIPTION:
! Initialize and run the CASA vegetation phenology.
!
! replaces pheno_params.h, pheno_paramsi.F, phenotci.F, phenoinit.F,
!          phenology.F, phenotcdyn.h, phenotvdyn.h from LSM-CASA code
!
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8
   use clm_varpar  , only : numpft
!
! !PUBLIC TYPES:
   implicit none
   save

  ! define parameters and constants used in Phenology

  integer , parameter :: nsecbeg  = 0    ! seconds at start of a day
  real(r8), parameter :: ngrowmin = 30.0_r8 ! minimum growing season of 30 days

  integer  :: evergreen(0:numpft) ! evergreen flag (0 or 1), 1 = evergreen 
  real(r8) :: tbase(0:numpft)     ! base temperature (deg C)
  real(r8) :: ddcrit(0:numpft)    ! degree days for start of growing season (deg C)
  real(r8) :: tdaycrit(0:numpft)  ! daily temp threshold for dropping leaves (deg C)

!  tbase     min temperature for starting growth (degrees C)
!        = 5C if not crops, use 10C if crops
!        = 40F for wheat, barley, rye, oats...              (not used)
!        = 45F for sunflower, potato...                     (not used)
!        = 50F for corn, sorghym, rice, soybeans, tomato ...(not used)

!  ddcrit    cumulative degree days threshold for start of
!            growing season (deg C)
!        = tbase until better info

!  tdaycrit  daily temperature threshold for dropping leaves (deg C)
!        = tbase until better info

    data tbase/ 999.00_r8,                                                  &
                  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8, &
                  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8, 10.00_r8, 10.00_r8/

    data ddcrit/ 999.00_r8,                                                  &
                   5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8, &
                   5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8, 10.00_r8, 10.00_r8/

    data tdaycrit/ 999.00_r8,                                                  &
                     5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8, &
                     5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8,  5.00_r8, 10.00_r8, 10.00_r8/
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initCASAPhenology              ! Initialize CASA phenology
  public :: CASAPhenology                  ! Compute CASA phenology
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initCASAPhenology
!
! !INTERFACE:
  subroutine initCASAPhenology()
!
! !DESCRIPTION:
! Initialize CASA vegetation phenology.
!
! !USES:
    use clmtype
    use decompMod , only : get_proc_bounds
    use clm_varctl, only : nsrest, nsrStartup
!
! !ARGUMENTS:
    implicit none
!
! !LOCAL VARIABLES:
    integer , pointer :: ivt(:)       ! pft vegetation type
    real(r8), pointer :: tday(:)      ! daily accumulated temperature (deg C)
    real(r8), pointer :: tdayavg(:)   ! daily averaged temperature (deg C)
    real(r8), pointer :: tcount(:)    ! counter for daily avg temp
    real(r8), pointer :: degday(:)    ! accumulated degree days (deg C) 
    real(r8), pointer :: ndegday(:)   ! counter for number of degree days
    real(r8), pointer :: stressT(:)   ! temperature stress function for leaf
                                      ! loss apply to Litterfall of decid veg
    real(r8), pointer :: stressW(:)   ! water stress function for leaf loss
    real(r8), pointer :: stressCD(:)  ! cold and drought stress function (sec-1)
    real(r8), pointer :: iseabeg(:)   ! index for start of growing season
    real(r8), pointer :: nstepbeg(:)  ! nstep at start of growing season
    real(r8), pointer :: lgrow(:)     ! growing season index (0 or 1) to be 
                                      ! passed daily to CASA to get NPP
    integer p                ! pft index
    integer begp, endp       ! per-proc beginning and ending pft indices
    integer begc, endc       ! per-proc beginning and ending column indices
    integer begl, endl       ! per-proc beginning and ending landunit indices
    integer begg, endg       ! per-proc gridcell ending gridcell indices
!
! !CALLED FROM:
! initialize in initializeMod
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    ! implicit intent in
    !============================================================
    ivt       => clm3%g%l%c%p%itype

    ! implicit intent out
    !============================================================
    tday      => clm3%g%l%c%p%pps%tday
    tdayavg   => clm3%g%l%c%p%pps%tdayavg
    tcount    => clm3%g%l%c%p%pps%tcount
    degday    => clm3%g%l%c%p%pps%degday
    ndegday   => clm3%g%l%c%p%pps%ndegday
    stressT   => clm3%g%l%c%p%pps%stressT
    stressW   => clm3%g%l%c%p%pps%stressW
    stressCD  => clm3%g%l%c%p%pps%stressCD
    iseabeg   => clm3%g%l%c%p%pps%iseabeg
    nstepbeg  => clm3%g%l%c%p%pps%nstepbeg
    lgrow     => clm3%g%l%c%p%pps%lgrow

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! set up array of crop types used in phenology
    ! set up array of evergreen types used in phenology

    do p = 0,numpft
       if (p == 1 .or. p == 2 .or. p == 4 .or. p == 5 .or. p == 9) then
          evergreen(p) = 1
       else
          evergreen(p) = 0
       end if
    end do

    do p = begp, endp   
!      stressT(p)  = 0.0
!      stressW(p)  = 0.0
       if (nsrest == nsrStartup) then
          if (evergreen(ivt(p)) == 1 ) then
             lgrow(p) = 1._r8
          else
             lgrow(p) = 0._r8
          end if
          iseabeg(p)  = 0.0_r8
          nstepbeg(p) = 0.0_r8
          ndegday(p)  = 0.0_r8
          tday(p)     = 0.0_r8
          tcount(p)   = 0.0_r8
          tdayavg(p)  = 0.0_r8
          degday(p)   = 0.0_r8
          stressCD(p) = 0.0_r8
          stressT(p)  = 0.0_r8    !! added to restart
          stressW(p)  = 0.0_r8    !! added to restart
       end if
    end do
    
  end subroutine initCASAPhenology

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CASAPhenology
!
! !INTERFACE:
  subroutine CASAPhenology(lbp, ubp, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Compute CASA vegetation phenology.
!
! !USES:
    use clmtype
    use clm_varcon  , only : tfrz, secspday
    use clm_time_manager, only : get_step_size, get_nstep, get_curr_date
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp   ! pft-index bounds
    integer, intent(in) :: num_soilp  ! number of soil points in pft filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! pft filter for soil points
!
! !LOCAL VARIABLES:
    integer , pointer :: ivt(:)       ! pft vegetation type
    integer , pointer :: pgridcell(:) ! pft's gridcell index
    real(r8), pointer :: latdeg(:)    ! latitude (degrees)
    real(r8), pointer :: tday(:)      ! daily accumulated temperature (deg C)
    real(r8), pointer :: tdayavg(:)   ! daily averaged temperature (deg C)
    real(r8), pointer :: tcount(:)    ! counter for daily avg temp
    real(r8), pointer :: degday(:)    ! accumulated degree days (deg C) 
    real(r8), pointer :: ndegday(:)   ! counter for number of degree days
    real(r8), pointer :: stressT(:)   ! temperature stress function for leaf
                                      ! loss apply to Litterfall of decid veg
    real(r8), pointer :: stressW(:)   ! water stress function for leaf loss
    real(r8), pointer :: stressCD(:)  ! cold and drought stress function (sec-1)
    real(r8), pointer :: iseabeg(:)   ! index for start of growing season
    real(r8), pointer :: nstepbeg(:)  ! nstep at start of growing season
    real(r8), pointer :: lgrow(:)     ! growing season index (0 or 1) to be 
                                      ! passed daily to CASA to get NPP
    real(r8), pointer :: t_ref2m(:)   ! 2m height surface air temperature (K)
    real(r8), pointer :: btran(:)     ! transpiration factor (0 to 1)
    !==============================================================
    integer  :: g          ! gridcell index
    integer  :: p          ! pft index
    integer  :: f          ! filter index
    integer  :: nstepmin
    integer  :: mcsec_n    !current seconds of next timestep (0, ..., seconds/day)
    integer  :: mcsec      !current seconds of current date (0, ..., seconds/day)
    integer  :: kyr        !year
    integer  :: kmo        !month (1, ..., 12)
    integer  :: kda        !day of month (1, ..., 31)
    integer  :: nstep      !time step number
    real(r8) :: dtime      !land model time step (sec) 
    real(r8) :: tnorm
!
! !CALLED FROM:
! casa_ecosystemDyn
!
! !REVISION HISTORY:
! 2004.06.08 Vectorized and reformatted by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

    ! implicit intent in
    !============================================================
    ivt       => clm3%g%l%c%p%itype
    pgridcell => clm3%g%l%c%p%gridcell
    latdeg    => clm3%g%latdeg

    ! implicit intent inout
    !============================================================
    tday     => clm3%g%l%c%p%pps%tday    ! daily accumulated temperature (deg C)
    tdayavg  => clm3%g%l%c%p%pps%tdayavg ! daily averaged temperature (deg C)
    tcount   => clm3%g%l%c%p%pps%tcount  ! counter for daily avg temp
    degday   => clm3%g%l%c%p%pps%degday  ! accumulated degree days (deg C) 
    ndegday  => clm3%g%l%c%p%pps%ndegday ! counter for number of degree days

    ! implicit intent out
    !============================================================
    stressT  => clm3%g%l%c%p%pps%stressT  ! temperature stress function for leaf loss apply to Litterfall of deciduous veg
    stressW  => clm3%g%l%c%p%pps%stressW  ! water stress function for leaf loss
    stressCD => clm3%g%l%c%p%pps%stressCD ! cold and drought stress function (sec-1)
    iseabeg  => clm3%g%l%c%p%pps%iseabeg  ! index for start of growing season
    nstepbeg => clm3%g%l%c%p%pps%nstepbeg ! nstep at start of growing season
    lgrow    => clm3%g%l%c%p%pps%lgrow    ! growing season index (0 or 1) to be passed daily to CASA to get NPP
    t_ref2m  => clm3%g%l%c%p%pes%t_ref2m  ! 2m surface air temperature (K)
    btran    => clm3%g%l%c%p%pps%btran    ! transpiration factor (0 to 1)

! -----------------------------------------------------------------------

    ! Get step size and current timestep

    dtime = get_step_size()
    nstep = get_nstep()

    ! convert ngrowmin (days) to timesteps

    nstepmin = ngrowmin*secspday/dtime

    ! ----------------------------------------------------------------------
    ! initialize arrays at start of the winter season
    ! NOTE: kda is not a julian day or this will not work.

    ! removed kmobeg (start of winter season) and hardwired to 1 or 7
    ! easy to reintroduce if expect it to vary (slevis)
    ! assign 0 deg latitude to the Northern hemisphere
    ! as the thermal equator is south of the equator.

    call get_curr_date(kyr, kmo, kda, mcsec)
    if (mcsec == nsecbeg .and. kda == 1) then
       do f = 1, num_soilp
          p = filter_soilp(f)
          g = pgridcell(p)
          if ((latdeg(g) >= 0._r8 .and. kmo == 1) .or. & ! Jan start in NH
              (latdeg(g) <  0._r8 .and. kmo == 7)) then  ! Jul start in SH
             degday(p)  = 0.0_r8
             ndegday(p) = 0.0_r8
             iseabeg(p) = 0.0_r8
          end if
       end do
    end if

    ! ----------------------------------------------------------------------
    ! DAILY AVERAGING 
    
    ! initialize arrays on first timestep of each new day 

    if (mcsec == nsecbeg) then
       do f = 1, num_soilp
          p = filter_soilp(f)
          ! tdayavg is needed for degday calculation and possibly for 
          ! monthly history accumulations - don't reset to zero.
          ! tdayavg(k) = 0.0

          tday(p) = 0.0_r8
          tcount(p) = 0.0_r8
       end do
    end if

    do f = 1, num_soilp
       p = filter_soilp(f)

       ! get min/max daily temp for dtr
       ! accumulate temperature during the day

       tday(p) = tday(p) + (t_ref2m(p)-tfrz)      !deg C
       tcount(p) = tcount(p) + 1.0_r8
    end do

    ! get daily averaged temperature (at end of day)
    ! mcsec_n = 0 is start of next day

    call get_curr_date(kyr, kmo, kda, mcsec_n, offset=int(dtime))
    if (mcsec_n /= nsecbeg) return   ! exit if not end of day

    ! reset counters for start of next day
    ! tdayavg is needed for degday calculation and possibly for 
    ! monthly history accumulations - don't reset to zero.

    do f = 1, num_soilp
       p = filter_soilp(f)
       if (tcount(p) /= 0.0_r8) then
          tdayavg(p) = tday(p)/tcount(p)
          tday(p) = 0.0_r8
          tcount(p) = 0.0_r8
       end if
    end do

    ! END OF DAILY AVERAGING
    ! IF NOT END OF DAY, EXIT FROM SUBROUTINE
    ! THIS EXITING WILL MISS THE FIRST DAY OF THE GROWING SEASON
    ! ----------------------------------------------------------------------
    
    ! accumulate degree days since daily Temp > base temperature for vegtype
    ! for corn, tbase = 50F
    ! Foley (1996) essentially has Tbase = 5C for winter deciduous trees
    
    do f = 1, num_soilp
       p = filter_soilp(f)
       if (tdayavg(p) > tbase(ivt(p))) then
          if (iseabeg(p) == 0.0_r8) nstepbeg(p) = real(nstep)
          iseabeg(p) = 1.0_r8
          degday(p) = degday(p) + tdayavg(p)
          ndegday(p) = ndegday(p)+1.0_r8
       end if
    end do

  ! ----------------------------------------------------------------------
  ! start growing season if degday > threshold
  ! here we assume that ddcrit(k) = tbase(k)
  !
  ! lgrow is either 0 or 1.  
  ! It is to be passed daily to the GPP subroutine.
  ! April02: we are passing this to the first call to CASA
  ! where we go from GPP to NPP.
  ! If lgrow=0 even if GPP > 0, we set NPP to zero, and
  ! pretend that all GPP went into autotrophic respiration.

    do f = 1, num_soilp
       p = filter_soilp(f)
       if (evergreen(ivt(p)) == 0) then
          if (degday(p) > ddcrit(ivt(p))) lgrow(p) = 1.0_r8
       end if
    end do

    ! ----------------------------------------------------------------------
    ! end of growing season
    ! if daily air temperature < critical value for veg
    ! Foley et al. (1996) used a critical temp of 5C for all veg
    ! drop leaves gradually by stress index
    ! ala Equation 7 and 8 of Dickinson et al J Climate 1998.
    ! Note Dickinson used canopy temperature rather than air temperature
    !  6/23/02
    !  Dickinson et al. (JClim 1998) Eqn(7) has T as argument in exponential
    !  and the exp goes to infinity very fast
    !              stressT(k)=exp(-(tdayavg(k)-tdaycrit(i)))
    !  modify the equation to use normalized T as argument
    !
    ! - Dickinson et al. Equation 6:
    !   stressCD is cold and drought stress, unit is sec-1
    !   cross-check:  1/stressCD ~1-2 months
    !   stressCD(k) to be added to "annK(m,LEAF)" and "annK(m,FROOT)"
    !       in casa_bgfluxes.F
    !----------------------------------------------------------------------
    ! Notes from iyf 02/07/11
    ! 
    ! Dickinson 1998 Interactive canopies for a climate model, J Climate
    ! 
    ! Cold and drought stress for phenology - Eqn 7-8
    ! 
    ! Here we modify the arguments of Dickinson's equation
    !      StressT = exp(-(Tavg-Tcrit)/Tcrit))
    !      StressW = exp((1-btran)
    ! No change to
    !      StressCD = (StressT+StressD)x2e-7    (sec-1)
    ! 
    ! Typical values:
    ! ---------------
    ! let stressT = 0.  btran = 0.1 --> tau=1/stressCD = 26 days.
    !                   btran = 0.9 -->                  56 days
    ! 
    ! Default values:
    ! ---------------
    ! evergreens:  stressT = 0, stressW = 0
    ! 
    ! Notes:
    ! ------
    ! In Dickinson's formulation of StressT, arg=-(T-Tcrit).  The normalization
    !  to Tcrit is necessary.
    ! 
    ! In Dickinson's formulation of StressW, arg=100(W-1), where W is a water
    ! stress
    ! term, est from the reduction by soil drying of the max evapotranspiration
    ! allowed by roots.  W>1 is permanent wilting.
    ! 
    ! In the LSM (documentation page 107),
    !    BTRAN = \sum (w_i*r_i)
    !          where w_i = (theta_i - theta_dry)/(theta_opt - theta-dry)
    !          theta_i is the volumetric soil moisture for layer i,
    !          and r_i is the relative root abundance.
    !    BTRAN varies between 0 (dry) and 1 (wet).
    !
    !----------------------------------------------------------------------
    
    ! To distinguish the autumn from the spring,
    ! put in a minimum growing season of 30 days
    
    do f = 1, num_soilp
       p = filter_soilp(f)

       ! Initialize every time step
       stressT(p)  = 0._r8
       stressW(p)  = 0._r8
       stressCD(p) = 0._r8

       if (evergreen(ivt(p)) == 0 .and. nstep-nstepbeg(p) > nstepmin &
          .and. tdayavg(p) < tdaycrit(ivt(p))) then
          lgrow(p) = 0.0_r8
          tnorm = tdaycrit(ivt(p))
          if (abs(tnorm) <= 1.e-5_r8) tnorm = 1.0_r8
          stressT(p) = exp(-(tdayavg(p)-tdaycrit(ivt(p)))/tnorm)
          stressW(p) = exp(1.0_r8-btran(p))
          stressCD(p) = (stressT(p)+stressW(p))*2.e-7_r8
       end if

    end do

  ! We apply StressT to litterfall of deciduous vegetation.
  ! Need to adjust seasonally invariant turnover times from before.
  ! April 02:  Dickinson has a water stress as well.  Not yet implemnted.
  ! July 02 - implemented water stress and cold/drought stress.


  end subroutine CASAPhenology

#endif

end module CASAPhenologyMod
