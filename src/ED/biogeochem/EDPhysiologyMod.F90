module EDPhysiologyMod

#include "shr_assert.h"

  ! ============================================================================
  ! Miscellaneous physiology routines from ED. 
  ! ============================================================================

  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_varctl          , only : iulog 
  use spmdMod             , only : masterproc
  use TemperatureType     , only : temperature_type
  use SoilStateType       , only : soilstate_type
  use WaterstateType      , only : waterstate_type
  use pftconMod           , only : pftcon
  use EDEcophysContype    , only : EDecophyscon
  use EDCohortDynamicsMod , only : allocate_live_biomass, zero_cohort
  use EDCohortDynamicsMod , only : create_cohort, fuse_cohorts, sort_cohorts
  use EDTypesMod          , only : dg_sf, dinc_ed, external_recruitment
  use EDTypesMod          , only : ncwd, cp_nlevcan, numpft_ed, senes
  use EDTypesMod          , only : ed_site_type, ed_patch_type, ed_cohort_type

  implicit none
  private

  public :: canopy_derivs
  public :: non_canopy_derivs
  public :: trim_canopy
  public :: phenology
  public :: phenology_leafonoff
  public :: Growth_Derivatives
  public :: recruitment
  public :: cwd_input
  public :: cwd_out
  public :: fragmentation_scaler
  public :: seeds_in
  public :: seed_decay
  public :: seed_germination
  public :: flux_into_litter_pools

  logical, parameter :: DEBUG  = .false. ! local debug flag

  ! ============================================================================

contains

  ! ============================================================================
  subroutine canopy_derivs( currentPatch )
    !
    ! !DESCRIPTION:
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type) , intent(inout), target :: currentPatch
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer ::currentCohort
    !----------------------------------------------------------------------

    ! call plant growth functions

    currentCohort => currentPatch%shortest

    do while(associated(currentCohort))
       call Growth_Derivatives(currentCohort)
       currentCohort => currentCohort%taller
    enddo

  end subroutine canopy_derivs

  ! ============================================================================
  subroutine non_canopy_derivs( currentPatch, temperature_inst )
    !
    ! !DESCRIPTION:
    ! Returns time differentials of the state vector
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type)    , intent(inout) :: currentPatch
    type(temperature_type) , intent(in)    :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer c,p
    !----------------------------------------------------------------------

    currentPatch%leaf_litter_in(:)   = 0.0_r8
    currentPatch%root_litter_in(:)   = 0.0_r8
    currentPatch%dleaf_litter_dt(:)  = 0.0_r8
    currentPatch%droot_litter_dt(:)  = 0.0_r8
    currentPatch%leaf_litter_out(:)  = 0.0_r8
    currentPatch%root_litter_out(:)  = 0.0_r8
    currentPatch%cwd_AG_in(:)        = 0.0_r8
    currentPatch%cwd_BG_in(:)        = 0.0_r8
    currentPatch%cwd_AG_out(:)       = 0.0_r8
    currentPatch%cwd_BG_out(:)       = 0.0_r8
    currentPatch%seeds_in(:)         = 0.0_r8  
    currentPatch%seed_decay(:)       = 0.0_r8
    currentPatch%seed_germination(:) = 0.0_r8

    ! update seed fluxes 
    call seeds_in(currentPatch)
    call seed_decay(currentPatch)
    call seed_germination(currentPatch)

    ! update fragmenting pool fluxes
    call cwd_input(currentPatch)
    call cwd_out( currentPatch, temperature_inst)

    do p = 1,numpft_ed
       currentPatch%dseed_dt(p) = currentPatch%seeds_in(p) - currentPatch%seed_decay(p) - currentPatch%seed_germination(p)
    enddo   
    
    do c = 1,ncwd
       currentPatch%dcwd_AG_dt(c) = currentPatch%cwd_AG_in(c) - currentPatch%cwd_AG_out(c) 
       currentPatch%dcwd_BG_dt(c) = currentPatch%cwd_BG_in(c) - currentPatch%cwd_BG_out(c) 
    enddo

    do p = 1,numpft_ed
       currentPatch%dleaf_litter_dt(p) = currentPatch%leaf_litter_in(p) - currentPatch%leaf_litter_out(p) 
       currentPatch%droot_litter_dt(p) = currentPatch%root_litter_in(p) - currentPatch%root_litter_out(p) 
    enddo

    ! currentPatch%leaf_litter_in(:)  = 0.0_r8
    ! currentPatch%root_litter_in(:)  = 0.0_r8
    ! currentPatch%leaf_litter_out(:) = 0.0_r8
    ! currentPatch%root_litter_out(:) = 0.0_r8
    ! currentPatch%CWD_AG_in(:)       = 0.0_r8
    ! currentPatch%cwd_bg_in(:)       = 0.0_r8
    ! currentPatch%CWD_AG_out(:)      = 0.0_r8
    ! currentPatch%cwd_bg_out(:)      = 0.0_r8

  end subroutine non_canopy_derivs

  ! ============================================================================
  subroutine trim_canopy( currentSite )
    !
    ! !DESCRIPTION:
    ! Canopy trimming / leaf optimisation. Removes leaves in negative annual carbon balance. 
    !
    ! !USES:
    !
    use EDParamsMod,          only : ED_val_grperc
    use EDGrowthFunctionsMod, only : tree_lai
    !
    ! !ARGUMENTS    
    type (ed_site_type),intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type) , pointer :: currentCohort
    type (ed_patch_type)  , pointer :: currentPatch

    real(r8) :: inc        ! rate at which canopy acclimates to uptake 
    real(r8) :: trim_limit ! this is the limit of the canopy trimming routine, so that trees 
                           ! can't just lose all their leaves and have no reproductive costs.
    integer  :: z          ! leaf layer
    integer  :: trimmed    ! was this layer trimmed in this year? If not expand the canopy. 

    trim_limit = 0.3_r8    ! Arbitrary limit to reductions in leaf area with stress. Without this nothing ever dies.  
    inc = 0.03_r8          ! Arbitrary incremental change in trimming function. Controls 
                           ! rate at which leaves are optimised to their environment. 
    !----------------------------------------------------------------------

    currentPatch => currentSite%youngest_patch

    do while(associated(currentPatch))
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort)) 
          trimmed = 0    
          currentCohort%treelai = tree_lai(currentCohort)    
          currentCohort%nv = ceiling((currentCohort%treelai+currentCohort%treesai)/dinc_ed)
          if (currentCohort%nv > cp_nlevcan)then
             write(iulog,*) 'nv > cp_nlevcan',currentCohort%nv,currentCohort%treelai,currentCohort%treesai, &
                  currentCohort%c_area,currentCohort%n,currentCohort%bl
          endif

          !Leaf cost vs netuptake for each leaf layer. 
          do z = 1,cp_nlevcan
             if (currentCohort%year_net_uptake(z) /= 999._r8)then !there was activity this year in this leaf layer. 
                !Leaf Cost kgC/m2/year-1
                !decidous costs. 
                if (pftcon%season_decid(currentCohort%pft) == 1.or.pftcon%stress_decid(currentCohort%pft) == 1)then 
                   currentCohort%leaf_cost =  1._r8/(pftcon%slatop(currentCohort%pft)*1000.0_r8)
                   currentCohort%leaf_cost = currentCohort%leaf_cost + 1.0_r8/(pftcon%slatop(currentCohort%pft)*1000.0_r8) * &
                        pftcon%froot_leaf(currentCohort%pft) / EDecophyscon%root_long(currentCohort%pft)
                   currentCohort%leaf_cost = currentCohort%leaf_cost * (ED_val_grperc(currentCohort%pft) + 1._r8)
                else !evergreen costs
                   currentCohort%leaf_cost = 1.0_r8/(pftcon%slatop(currentCohort%pft)* &
                        pftcon%leaf_long(currentCohort%pft)*1000.0_r8) !convert from sla in m2g-1 to m2kg-1
                   currentCohort%leaf_cost = currentCohort%leaf_cost + 1.0_r8/(pftcon%slatop(currentCohort%pft)*1000.0_r8) * &
                        pftcon%froot_leaf(currentCohort%pft) / EDecophyscon%root_long(currentCohort%pft)
                   currentCohort%leaf_cost = currentCohort%leaf_cost * (ED_val_grperc(currentCohort%pft) + 1._r8)
                endif
                if (currentCohort%year_net_uptake(z) < currentCohort%leaf_cost)then
                   if (currentCohort%canopy_trim > trim_limit)then

                      if ( DEBUG ) then
                         write(iulog,*) 'trimming leaves',currentCohort%canopy_trim,currentCohort%leaf_cost
                      endif

                      ! keep trimming until none of the canopy is in negative carbon balance.              
                      if (currentCohort%hite > EDecophyscon%hgt_min(currentCohort%pft))then
                         currentCohort%canopy_trim = currentCohort%canopy_trim - inc    
                         if (pftcon%evergreen(currentCohort%pft) /= 1)then
                            currentCohort%laimemory = currentCohort%laimemory*(1.0_r8 - inc) 
                         endif
                         trimmed = 1
                      endif
                   endif
                endif
             endif !leaf activity? 
          enddo !z
          if (currentCohort%NV.gt.2)then
             ! leaf_cost may be uninitialized, removing its diagnostic from the log
             ! to allow checking with fpe_traps (RGK)
             write(iulog,*) 'nv>4',currentCohort%year_net_uptake(1:6),currentCohort%canopy_trim
          endif

          currentCohort%year_net_uptake(:) = 999.0_r8
          if (trimmed == 0.and.currentCohort%canopy_trim < 1.0_r8)then
             currentCohort%canopy_trim = currentCohort%canopy_trim + inc
          endif 

          if ( DEBUG ) then
             write(iulog,*) 'trimming',currentCohort%canopy_trim
          endif
         
          ! currentCohort%canopy_trim = 1.0_r8 !FIX(RF,032414) this turns off ctrim for now. 
          currentCohort => currentCohort%shorter
       enddo
       currentPatch => currentPatch%older
    enddo

  end subroutine trim_canopy

  ! ============================================================================
  subroutine phenology( currentSite, temperature_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    ! Phenology. 
    !
    ! !USES:
    use clm_varcon, only : tfrz
    use clm_time_manager, only : get_curr_date
    use clm_time_manager, only : get_ref_date, timemgr_datediff 
    use EDTypesMod, only : udata
    use PatchType , only : patch   
    !
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), target :: currentSite
    type(temperature_type)  , intent(in)            :: temperature_inst
    type(waterstate_type)   , intent(in)            :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: t_veg24(:) 
    integer  :: t            ! day of year
    integer  :: ncolddays    ! no days underneath the threshold for leaf drop
    integer  :: ncolddayslim ! critical no days underneath the threshold for leaf drop
    integer  :: i
    integer  :: timesincedleafon,timesincedleafoff,timesinceleafon,timesinceleafoff
    integer  :: refdate
    integer  :: curdate
    
    integer  :: yr                       ! year (0, ...)
    integer  :: mon                      ! month (1, ..., 12)
    integer  :: day                      ! day of month (1, ..., 31)
    integer  :: sec                      ! seconds of the day
    integer  :: patchi                   ! the first CLM/ALM patch index of the associated column
    integer  :: coli                     ! the CLM/ALM column index of the associated site

    real(r8) :: gdd_threshold
    real(r8) :: a,b,c        ! params of leaf-pn model from botta et al. 2000. 
    real(r8) :: cold_t       ! threshold below which cold days are counted 
    real(r8) :: coldday      ! definition of a 'chilling day' for botta model 
    integer  :: ncdstart     ! beginning of counting period for chilling degree days.
    integer  :: gddstart     ! beginning of counting period for growing degree days.
    real(r8) :: drought_threshold
    real(r8) :: off_time     ! minimum number of days between leaf off and leaf on for drought phenology 
    real(r8) :: temp_in_C    ! daily averaged temperature in celcius
    real(r8) :: mindayson 
    real(r8) :: modelday

    !------------------------------------------------------------------------

    ! INTERF-TODO: THIS IS A BAND-AID, AS I WAS HOPING TO REMOVE CLM_PNO
    ! ALREADY REMOVED currentSite%clmcolumn, hence the need for these

    patchi = currentSite%oldest_patch%clm_pno-1
    coli   = patch%column(patchi)

    t_veg24       => temperature_inst%t_veg24_patch ! Input:  [real(r8) (:)]  avg pft vegetation temperature for last 24 hrs    

    call get_curr_date(yr, mon, day, sec)
    curdate = yr*10000 + mon*100 + day
    
    call get_ref_date(yr, mon, day, sec)
    refdate = yr*10000 + mon*100 + day
  
    call timemgr_datediff(refdate, 0, curdate, sec, modelday)
    if ( masterproc ) write(iulog,*) 'modelday',modelday

    ! Parameter of drought decid leaf loss in mm in top layer...FIX(RF,032414) 
    ! - this is arbitrary and poorly understood. Needs work. ED_
    drought_threshold = 0.15 
    off_time = 100.0_r8

    !Parameters of Botta et al. 2000 GCB,6 709-725 
    a = -68.0_r8
    b = 638.0_r8
    c = -0.001_r8
    coldday = 5.0_r8    !ed_ph_chiltemp

    mindayson = 30
     
    !Parameters from SDGVM model of senesence
    ncolddayslim = 5
    cold_t   = 7.5_r8  ! ed_ph_coldtemp

    t  = udata%time_period
    temp_in_C = t_veg24(patchi) - tfrz

    !-----------------Cold Phenology--------------------!              

    !Zero growing degree and chilling day counters
    if (currentSite%lat > 0)then
       ncdstart = 270  !Northern Hemisphere begining November
       gddstart = 1    !Northern Hemisphere begining January
    else
       ncdstart = 120  !Southern Hemisphere beginning May
       gddstart = 181  !Northern Hemisphere begining July
    endif
    
    ! FIX(SPM,032414) - this will only work for the first year, no?
    if (t == ncdstart)then
       currentSite%ncd = 0._r8
    endif

    !Accumulate growing/chilling days after start of counting period
    if (temp_in_C  <  coldday)then
       currentSite%ncd = currentSite%ncd + 1.0_r8
    endif

    gdd_threshold = a + b*exp(c*currentSite%ncd) !GDD accumulation function, which also depends on chilling days.

    !Accumulate temperature of last 10 days.
    currentSite%last_n_days(2:senes) =  currentSite%last_n_days(1:senes-1)
    currentSite%last_n_days(1) = temp_in_C                                      
    !count number of days for leaves off
    ncolddays = 0
    do i = 1,senes
       if (currentSite%last_n_days(i) < cold_t)then
          ncolddays = ncolddays + 1
       endif
    enddo

    ! Here is where we do the GDD accumulation calculation
    !
    ! reset GDD on set dates
    if (t == gddstart)then
       currentSite%ED_GDD_site = 0._r8
    endif
    !
    ! accumulate the GDD using daily mean temperatures
    if (t_veg24(patchi) .gt. tfrz) then
       currentSite%ED_GDD_site = currentSite%ED_GDD_site + t_veg24(currentSite%oldest_patch%clm_pno-1) - tfrz
    endif
    

    timesinceleafoff = modelday - currentSite%leafoffdate
    !LEAF ON: COLD DECIDUOUS. Needs to
    !1) have exceeded the growing degree day threshold 
    !2) The leaves should not be on already
    !3) There should have been at least on chilling day in the counting period.  
    if (currentSite%ED_GDD_site > gdd_threshold)then
       if (currentSite%status == 1) then
          if (currentSite%ncd >= 1) then
             currentSite%status = 2     !alter status of site to 'leaves on'
             ! NOTE(bja, 2015-01) should leafondate = modelday to be consistent with leaf off?
             currentSite%leafondate = t !record leaf on date   
             if ( DEBUG ) write(iulog,*) 'leaves on'
          endif !ncd
       endif !status
    endif !GDD

    timesinceleafon = modelday - currentSite%leafondate


    !LEAF OFF: COLD THRESHOLD
    !Needs to:
    !1) have exceeded the number of cold days threshold
    !2) have exceeded the minimum leafon time.
    !3) The leaves should not be off already
    !4) The day of the year should be larger than the counting period. (not sure if we need this/if it will break the restarting)
    
    if (ncolddays > ncolddayslim)then
     if (timesinceleafon > mindayson)then
       if (currentSite%status == 2)then
          currentSite%status = 1        !alter status of site to 'leaves on'
          currentSite%leafoffdate = modelday   !record leaf off date   
          if ( DEBUG ) write(iulog,*) 'leaves off'
       endif
    endif
    endif

    !LEAF OFF: COLD LIFESPAN THRESHOLD
    if(timesinceleafoff > 400)then !remove leaves after a whole year when there is no 'off' period.  
       if(currentSite%status == 2)then
          currentSite%status = 1        !alter status of site to 'leaves on'
          currentSite%leafoffdate = modelday   !record leaf off date   
          if ( DEBUG ) write(iulog,*) 'leaves off'
       endif
    endif

    !-----------------Drought Phenology--------------------!
    ! Principles of drought-deciduos phenology model...
    ! The 'dstatus' flag is 2 when leaves are on, and 1 when leaves area off. 
    ! The following sets those site-level flags, which are acted on in phenology_deciduos. 
    ! A* The leaves live for either the length of time the soil moisture is over the threshold 
    ! or the lifetime of the leaves, whichever is shorter. 
    ! B*: If the soil is only wet for a very short time, then the leaves stay on for 100 days
    ! C*: The leaves are only permitted to come ON for a 60 day window around when they last came on, 
    ! to prevent 'flickering' on in response to wet season storms
    ! D*: We don't allow anything to happen in the first ten days to allow the water memory window to come into equlibirum. 
    ! E*: If the soil is always wet, the leaves come on at the beginning of the window, and then last for their lifespan. 
    ! ISSUES
    ! 1. It's not clear what water content we should track. Here we are tracking the top layer, 
    ! but we probably should track something like BTRAN,
    ! but BTRAN is defined for each PFT, and there could potentially be more than one stress-dec PFT.... ?
    ! 2. In the beginning, the window is set at an arbitrary time of the year, so the leaves might come on 
    ! in the dry season, using up stored reserves
    ! for the stress-dec plants, and potentially killing them. To get around this, we need to read in the 
    ! 'leaf on' date from some kind of start-up file
    ! but we would need that to happen for every resolution, etc. 
    ! 3. Will this methodology properly kill off the stress-dec trees where there is no water stress? 
    ! What about where the wet period coincides with the
    ! warm period? We would just get them overlapping with the cold-dec trees, even though that isn't appropriate.... 
    ! Why don't the drought deciduous trees grow
    ! in the North? Is cold decidousness maybe even the same as drought deciduosness there (and so does this 
    ! distinction actually matter??).... 

    !Accumulate surface water memory of last 10 days.
    currentSite%water_memory(1) = waterstate_inst%h2osoi_vol_col(coli,1) 
    do i = 1,9 !shift memory along one
       currentSite%water_memory(11-i) = currentSite%water_memory(10-i)
    enddo

    !In drought phenology, we often need to force the leaves to stay on or off as moisture fluctuates...     
    timesincedleafoff = 0
    if (currentSite%dstatus == 1)then !the leaves are off. How long have they been off? 
       !leaves have come on, but last year, so at a later date than now.
       if (currentSite%dleafoffdate > 0.and.currentSite%dleafoffdate > t)then 
          timesincedleafoff = t + (360 - currentSite%dleafoffdate)
       else
          timesincedleafoff = t - currentSite%dleafoffdate    
       endif
    endif

    timesincedleafon = 0
    !the leaves are on. How long have they been on? 
    if (currentSite%dstatus == 2)then  
       !leaves have come on, but last year, so at a later date than now.
       if (currentSite%dleafondate > 0.and.currentSite%dleafondate > t)then 
          timesincedleafon = t + (360 - currentSite%dleafondate)
       else
          timesincedleafon = t - currentSite%dleafondate      
       endif
    endif

    !LEAF ON: DROUGHT DECIDUOUS WETNESS
    !Here, we used a window of oppurtunity to determine if we are close to the time when then leaves came on last year
    if ((t >= currentSite%dleafondate - 30.and.t <= currentSite%dleafondate + 30).or.(t > 360 - 15.and. &
         currentSite%dleafondate < 15))then ! are we in the window?
       if (sum(currentSite%water_memory(1:10)/10._r8) >= drought_threshold.and.currentSite%dstatus == 1.and.t >= 10)then 
          ! leave some minimum time between leaf off and leaf on to prevent 'flickering'.  
          if (timesincedleafoff > off_time)then  
             currentSite%dstatus = 2     !alter status of site to 'leaves on'
             currentSite%dleafondate = t   !record leaf on date
          endif
       endif
    endif

   !we still haven't done budburst by end of window
    if (t == currentSite%dleafondate+30.and.currentSite%dstatus == 1)then 
       currentSite%dstatus = 2    ! force budburst!
       currentSite%dleafondate = t   ! record leaf on date
    endif

    !LEAF OFF: DROUGHT DECIDUOUS LIFESPAN - if the leaf gets to the end of its useful life. A*, E*
    if (currentSite%dstatus == 2.and.t >= 10)then  !D*
       !Are the leaves at the end of their lives? !FIX(RF,0401014)- this is hardwiring....
       if (timesincedleafon > 365.0*pftcon%leaf_long(7))then 
          currentSite%dstatus = 1         !alter status of site to 'leaves on'
          currentSite%dleafoffdate = t    !record leaf on date          
       endif
    endif

    !LEAF OFF: DROUGHT DECIDUOUS DRYNESS - if the soil gets too dry, and the leaves have already been on a while... 
    if (currentSite%dstatus == 2.and.t >= 10)then  !D*
       if (sum(currentSite%water_memory(1:10)/10._r8) <= drought_threshold)then 
          if (timesincedleafon > 100)then !B* Have the leaves been on for some reasonable length of time? To prevent flickering. 
             currentSite%dstatus = 1      !alter status of site to 'leaves on'
             currentSite%dleafoffdate = t !record leaf on date           
          endif
       endif
    endif

    call phenology_leafonoff(currentSite)

  end subroutine phenology

  ! ============================================================================
  subroutine phenology_leafonoff(currentSite)
    !
    ! !DESCRIPTION:
    ! Controls the leaf on and off economics
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type) , pointer :: currentPatch     
    type(ed_cohort_type), pointer :: currentCohort  

    real(r8)           :: store_output ! the amount of the store to put into leaves - is a barrier against negative storage and C starvation. 

    !------------------------------------------------------------------------

    currentPatch => CurrentSite%oldest_patch   

    store_output  = 0.5_r8

    do while(associated(currentPatch))    
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))        
                
          !COLD LEAF ON
          if (pftcon%season_decid(currentCohort%pft) == 1)then
             if (currentSite%status == 2)then !we have just moved to leaves being on . 
                if (currentCohort%status_coh == 1)then !Are the leaves currently off?        
                   currentCohort%status_coh = 2    !Leaves are on, so change status to stop flow of carbon out of bstore. 
                   if (currentCohort%laimemory <= currentCohort%bstore)then
                      currentCohort%bl = currentCohort%laimemory !extract stored carbon to make new leaves.
                   else
                      ! we can only put on as much carbon as there is in the store...
                      ! nb. Putting all of bstore into leaves is C-starvation suicidal. 
                      ! The tendency for this could be parameterized
                      currentCohort%bl = currentCohort%bstore * store_output
                   endif

                   ! Add deployed carbon to alive biomass pool
                   currentCohort%balive = currentCohort%balive + currentCohort%bl

                   if ( DEBUG ) write(iulog,*) 'EDPhysMod 1 ',currentCohort%bstore

                   currentCohort%bstore = currentCohort%bstore - currentCohort%bl  ! Drain store

                   if ( DEBUG ) write(iulog,*) 'EDPhysMod 2 ',currentCohort%bstore

                   currentCohort%laimemory = 0.0_r8

                endif !pft phenology
             endif ! growing season 

             !COLD LEAF OFF
             currentCohort%leaf_litter = 0.0_r8 !zero leaf litter for today. 
             if (currentSite%status == 1)then !past leaf drop day? Leaves still on tree?  
                if (currentCohort%status_coh == 2)then ! leaves have not dropped
                   currentCohort%status_coh      = 1                  
                   !remember what the lai was this year to put the same amount back on in the spring... 
                   currentCohort%laimemory   = currentCohort%bl  
                   ! decrement balive for leaf litterfall        
                   currentCohort%balive      = currentCohort%balive - currentCohort%bl 
                   ! add lost carbon to litter
                   currentCohort%leaf_litter = currentCohort%bl 
                   currentCohort%bl          = 0.0_r8                            
                endif !leaf status
             endif !currentSite status
          endif  !season_decid

          !DROUGHT LEAF ON
          if (pftcon%stress_decid(currentCohort%pft) == 1)then
             if (currentSite%dstatus == 2)then !we have just moved to leaves being on . 
                if (currentCohort%status_coh == 1)then !is it the leaf-on day? Are the leaves currently off?       
                   currentCohort%status_coh = 2    !Leaves are on, so change status to stop flow of carbon out of bstore. 
                   if (currentCohort%laimemory <= currentCohort%bstore)then
                      currentCohort%bl = currentCohort%laimemory !extract stored carbon to make new leaves.
                   else
                    currentCohort%bl = currentCohort%bstore * store_output    !we can only put on as much carbon as there is in the store...
                    endif
                   currentCohort%balive = currentCohort%balive + currentCohort%bl

                   if ( DEBUG ) write(iulog,*) 'EDPhysMod 3 ',currentCohort%bstore

                   currentCohort%bstore = currentCohort%bstore - currentCohort%bl ! empty store

                   if ( DEBUG ) write(iulog,*) 'EDPhysMod 4 ',currentCohort%bstore

                   currentCohort%laimemory = 0.0_r8

                endif !currentCohort status again?
             endif   !currentSite status

             !DROUGHT LEAF OFF
             if (currentSite%dstatus == 1)then        
                if (currentCohort%status_coh == 2)then ! leaves have not dropped
                   currentCohort%status_coh      = 1   
                   currentCohort%laimemory   = currentCohort%bl
                   ! decrement balive for leaf litterfall  
                   currentCohort%balive      = currentCohort%balive - currentCohort%bl   
                   ! add retranslocated carbon (very small) to store.      
                   currentCohort%bstore      = currentCohort%bstore        
                   ! add falling leaves to litter pools . convert to KgC/m2                    
                   currentCohort%leaf_litter = currentCohort%bl  
                   currentCohort%bl          = 0.0_r8                                        
                endif
             endif !status
          endif !drought dec.
          currentCohort => currentCohort%shorter
       enddo !currentCohort

       currentPatch => currentPatch%younger

    enddo !currentPatch

  end subroutine phenology_leafonoff


  ! ============================================================================
  subroutine seeds_in( cp_pnt )
    !
    ! !DESCRIPTION:
    !  Flux from plants into seed pool. 
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), target :: cp_pnt ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type),  pointer :: currentPatch
    type(ed_site_type),   pointer :: currentSite
    type(ed_cohort_type), pointer :: currentCohort
    integer :: p
    !----------------------------------------------------------------------

    currentPatch => cp_pnt
    currentSite  => currentPatch%siteptr
   
    currentPatch%seeds_in(:) = 0.0_r8
    currentPatch%seed_rain_flux(:) = 0.0_r8
    
    currentCohort => currentPatch%tallest
    do while (associated(currentCohort))
       p = currentCohort%pft
       currentPatch%seeds_in(p) = currentPatch%seeds_in(p) +  currentCohort%seed_prod * currentCohort%n/currentPatch%area
       currentCohort => currentCohort%shorter
    enddo !cohort loop

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))
       if (EXTERNAL_RECRUITMENT == 1) then !external seed rain - needed to prevent extinction  
          do p = 1,numpft_ed
           currentPatch%seeds_in(p) = currentPatch%seeds_in(p) + EDecophyscon%seed_rain(p) !KgC/m2/year
           currentPatch%seed_rain_flux(p) = currentPatch%seed_rain_flux(p) + EDecophyscon%seed_rain(p) !KgC/m2/year
          enddo
       endif
       currentPatch => currentPatch%younger
    enddo

  end subroutine seeds_in
  
  ! ============================================================================
  subroutine seed_decay( currentPatch )
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into leaf litter pool    
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type),intent(inout) :: currentPatch ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    integer  ::  p
    real(r8) :: seed_turnover !complete seed turnover rate in yr-1. 
    !----------------------------------------------------------------------

    seed_turnover = 0.51_r8  ! from Liscke and Loffler 2006  
    ! decays the seed pool according to exponential model
    ! sd_mort is in yr-1
    do p = 1,numpft_ed 
       currentPatch%seed_decay(p) =  currentPatch%seed_bank(p) * seed_turnover
    enddo
 
  end subroutine seed_decay

  ! ============================================================================
  subroutine seed_germination( currentPatch ) 
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into sapling pool    
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type),intent(inout) :: currentPatch ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    integer :: p
    real(r8) max_germination !cap on germination rates. KgC/m2/yr Lishcke et al. 2009
    real(r8) germination_timescale !yr-1
    !----------------------------------------------------------------------

    germination_timescale = 0.5_r8 !this is arbitrary
    max_germination = 1.0_r8 !this is arbitrary

    do p = 1,numpft_ed
       currentPatch%seed_germination(p) =  min(currentPatch%seed_bank(p) * germination_timescale,max_germination)
    enddo

  end subroutine seed_germination

  ! ============================================================================
  subroutine Growth_Derivatives( currentCohort)
    !
    ! !DESCRIPTION:
    !  Main subroutine controlling growth and allocation derivatives    
    !
    ! !USES:
    use EDGrowthFunctionsMod , only : Bleaf, dDbhdBd, dhdbd, hite, mortality_rates,dDbhdBl
    use EDTypesMod           , only : udata
    !
    ! !ARGUMENTS    
    type(ed_cohort_type),intent(inout), target :: currentCohort
    !
    ! !LOCAL VARIABLES:
    type(ed_site_type),  pointer :: currentSite
    real(r8) :: dbldbd   !rate of change of dead biomass per unit dbh 
    real(r8) :: dbrdbd   !rate of change of root biomass per unit dbh
    real(r8) :: dbswdbd  !rate of change of sapwood biomass per unit dbh
    real(r8) :: dhdbd_fn !rate of change of height per unit dbh
    real(r8) :: va       !fraction of growth going to alive biomass
    real(r8) :: vs       !fraction of growth going to structural biomass
    real(r8) :: u,h      !intermediates 
    real(r8) :: frac     !fraction the stored carbon is of target store amount
    real(r8) :: f_store  !fraction of NPP allocated to storage in this timestep (functionf of stored pool)
    real(r8) :: gr_fract !fraction of carbon balance that is allocated to growth (not reproduction)
    real(r8) :: target_balive  !target leaf biomass under allometric optimum.  
    real(r8) :: cmort    ! starvation mortality rate (fraction per year)
    real(r8) :: bmort    ! background mortality rate (fraction per year)
    real(r8) :: hmort    ! hydraulic failure mortality rate (fraction per year)
    real(r8) :: balive_loss
    !----------------------------------------------------------------------

    currentSite => currentCohort%siteptr

    ! Mortality for trees in the understorey. 
    !if trees are in the canopy, then their death is 'disturbance'. This probably needs a different terminology
    if (currentCohort%canopy_layer > 1)then 
       call mortality_rates(currentCohort,cmort,hmort,bmort)
       currentCohort%dndt = -1.0_r8 * (cmort+hmort+bmort) * currentCohort%n
    else
       currentCohort%dndt = 0._r8
    endif

    ! Height
    currentCohort%hite = Hite(currentCohort) 
    h = currentCohort%hite
                       
    call allocate_live_biomass(currentCohort,0)

   ! calculate target size of living biomass compartment for a given dbh.   
    target_balive = Bleaf(currentCohort) * (1.0_r8 + pftcon%froot_leaf(currentCohort%pft) + &
         EDecophyscon%sapwood_ratio(currentCohort%pft)*h)
    !target balive without leaves. 
    if (currentCohort%status_coh == 1)then 
       target_balive = Bleaf(currentCohort) * (pftcon%froot_leaf(currentCohort%pft) + &
            EDecophyscon%sapwood_ratio(currentCohort%pft) * h)
    endif

    ! NPP 
    if ( DEBUG ) write(iulog,*) 'EDphys 716 ',currentCohort%npp_acc

    currentCohort%npp  = currentCohort%npp_acc  * udata%n_sub   !Link to CLM. convert from kgC/indiv/day into kgC/indiv/year
    currentCohort%gpp  = currentCohort%gpp_acc  * udata%n_sub   !Link to CLM. convert from kgC/indiv/day into kgC/indiv/year
    currentCohort%resp = currentCohort%resp_acc * udata%n_sub   !Link to CLM. convert from kgC/indiv/day into kgC/indiv/year

    currentSite%flux_in = currentSite%flux_in + currentCohort%npp_acc * currentCohort%n

    ! Maintenance demands     
    if (pftcon%evergreen(currentCohort%pft) == 1)then !grass and EBT
       currentCohort%leaf_md = currentCohort%bl / pftcon%leaf_long(currentCohort%pft)
       currentCohort%root_md = currentCohort%br / EDecophyscon%root_long(currentCohort%pft)
       currentCohort%md      = currentCohort%root_md + currentCohort%leaf_md
    endif

    !FIX(RF,032414) - I took out the stem turnover demand as it seemed excesively high and caused odd size-reated 
    ! decline affect
    !with which I am not especially comfortable, particularly as the concept of sapwood turnover is unclear for trees that 
    !are still in an expansion phase. 

    if (pftcon%season_decid(currentCohort%pft) == 1)then 
       currentCohort%root_md = currentCohort%br /EDecophyscon%root_long(currentCohort%pft)
       currentCohort%leaf_md = 0._r8
       currentCohort%md = currentCohort%root_md + currentCohort%leaf_md
    endif

    if (pftcon%stress_decid(currentCohort%pft) == 1)then 
       currentCohort%root_md = currentCohort%br /EDecophyscon%root_long(currentCohort%pft)
       currentCohort%leaf_md = 0._r8
       currentCohort%md = currentCohort%root_md + currentCohort%leaf_md
    endif

    if (pftcon%stress_decid(currentCohort%pft) /= 1.and.pftcon%season_decid(currentCohort%pft) /= 1.and. &
         pftcon%evergreen(currentCohort%pft) /= 1)then
       write(iulog,*) 'problem with phenology definitions',currentCohort%pft,pftcon%stress_decid(currentCohort%pft), &
            pftcon%season_decid(currentCohort%pft),pftcon%evergreen(currentCohort%pft)
    endif

    ! FIX(RF,032414) -turned off for now as it makes balive go negative....
    ! FIX(RF,032414) jan2012 0.01_r8 * currentCohort%bdead
    currentCohort%woody_turnover = 0.0_r8
    currentCohort%md = currentCohort%md + currentCohort%woody_turnover

    ! Calculate carbon balance 
    ! this is the fraction of maintenance demand we -have- to do...

    if ( DEBUG ) write(iulog,*) 'EDphys 760 ',currentCohort%npp, currentCohort%md, &
                   EDecophyscon%leaf_stor_priority(currentCohort%pft)

    currentCohort%carbon_balance = currentCohort%npp - currentCohort%md *  EDecophyscon%leaf_stor_priority(currentCohort%pft)

    ! Allowing only carbon from NPP pool to account for npp flux into the maintenance pools
    ! ie this does not include any use of storage carbon or balive to make up for missing carbon balance in the transfer
    currentCohort%npp_leaf  = min(currentCohort%npp*currentCohort%leaf_md/currentCohort%md, &
                                  currentCohort%leaf_md*EDecophyscon%leaf_stor_priority(currentCohort%pft))
    currentCohort%npp_froot = min(currentCohort%npp*currentCohort%root_md/currentCohort%md, &
                                  currentCohort%root_md*EDecophyscon%leaf_stor_priority(currentCohort%pft))


    if (Bleaf(currentCohort) > 0._r8)then

       if ( DEBUG ) write(iulog,*) 'EDphys A ',currentCohort%carbon_balance

       if (currentCohort%carbon_balance > 0._r8)then !spend C on growing and storing

          !what fraction of the target storage do we have? 
          frac = max(0.0_r8,currentCohort%bstore/(Bleaf(currentCohort) * EDecophyscon%cushion(currentCohort%pft)))
          ! FIX(SPM,080514,fstore never used ) 
          f_store = max(exp(-1.*frac**4._r8) - exp( -1.0_r8 ),0.0_r8)  
          !what fraction of allocation do we divert to storage?
          !what is the flux into the store?
          currentCohort%storage_flux = currentCohort%carbon_balance * f_store                     

          if ( DEBUG ) write(iulog,*) 'EDphys B ',f_store

          !what is the tax on the carbon available for growth? 
          currentCohort%carbon_balance = currentCohort%carbon_balance * (1.0_r8 - f_store)  
       else  !cbalance is negative. Take C out of store to pay for maintenance respn.
          currentCohort%storage_flux = currentCohort%carbon_balance 
          currentCohort%carbon_balance = 0._r8 
       endif

    else

       currentCohort%storage_flux = 0._r8
       currentCohort%carbon_balance = 0._r8
       write(iulog,*) 'ED: no leaf area in gd', currentCohort%indexnumber,currentCohort%n,currentCohort%bdead, &
             currentCohort%dbh,currentCohort%balive

    endif

    !Do we have enough carbon left over to make up the rest of the turnover demand? 
    balive_loss = 0._r8
    if (currentCohort%carbon_balance > currentCohort%md*(1.0_r8- EDecophyscon%leaf_stor_priority(currentCohort%pft)))then ! Yes...
       currentCohort%carbon_balance = currentCohort%carbon_balance - currentCohort%md * (1.0_r8 - &
             EDecophyscon%leaf_stor_priority(currentCohort%pft))

       currentCohort%npp_leaf  = currentCohort%npp_leaf  + &
            currentCohort%leaf_md *  (1.0_r8-EDecophyscon%leaf_stor_priority(currentCohort%pft))
       currentCohort%npp_froot = currentCohort%npp_froot + &
            currentCohort%root_md *  (1.0_r8-EDecophyscon%leaf_stor_priority(currentCohort%pft))

    else ! we can't maintain constant leaf area and root area. Balive is reduced

       currentCohort%npp_leaf  = currentCohort%npp_leaf  + &
             max(0.0_r8,currentCohort%carbon_balance*(currentCohort%leaf_md/currentCohort%md))
       currentCohort%npp_froot = currentCohort%npp_froot + &
             max(0.0_r8,currentCohort%carbon_balance*(currentCohort%root_md/currentCohort%md))

       balive_loss = currentCohort%md *(1.0_r8- EDecophyscon%leaf_stor_priority(currentCohort%pft))- currentCohort%carbon_balance
       currentCohort%carbon_balance = 0._r8
    endif

    !********************************************/
    ! Allometry & allocation of remaining carbon*/
    !********************************************/
    !Use remaining carbon to refill balive or to get larger. 

    !only if carbon balance is +ve
    if ((currentCohort%balive >= target_balive).AND.(currentCohort%carbon_balance >  0._r8))then 
       ! fraction of carbon going into active vs structural carbon        
       if (currentCohort%dbh <= EDecophyscon%max_dbh(currentCohort%pft))then ! cap on leaf biomass
          dbldbd = dDbhdBd(currentCohort)/dDbhdBl(currentCohort) 
          dbrdbd = pftcon%froot_leaf(currentCohort%pft) * dbldbd
          dhdbd_fn = dhdbd(currentCohort)
          dbswdbd = EDecophyscon%sapwood_ratio(currentCohort%pft) * (h*dbldbd + currentCohort%bl*dhdbd_fn)
          u  = 1.0_r8 / (dbldbd + dbrdbd + dbswdbd)     
          va = 1.0_r8 / (1.0_r8 + u)
          vs = u / (1.0_r8 + u)
          gr_fract = 1.0_r8 - EDecophyscon%seed_alloc(currentCohort%pft)
       else
          dbldbd = 0._r8; dbrdbd = 0._r8 ;dbswdbd = 0._r8      
          va = 0.0_r8
          vs = 1.0_r8
          gr_fract = 1.0_r8 - (EDecophyscon%seed_alloc(currentCohort%pft) + EDecophyscon%clone_alloc(currentCohort%pft))
       endif

       !FIX(RF,032414) - to fix high bl's. needed to prevent numerical errors without the ODEINT.  
       if (currentCohort%balive > target_balive*1.1_r8)then  
          va = 0.0_r8; vs = 1._r8
          write(iulog,*) 'using high bl cap',target_balive,currentCohort%balive                        
       endif

    else         
       dbldbd = 0._r8; dbrdbd = 0._r8; dbswdbd = 0._r8
       va = 1.0_r8; vs = 0._r8                        
       gr_fract = 1.0_r8
    endif

    ! calculate derivatives of living and dead carbon pools  
    currentCohort%dbalivedt = gr_fract * va * currentCohort%carbon_balance - balive_loss
    currentCohort%dbdeaddt  = gr_fract * vs * currentCohort%carbon_balance
    currentCohort%dbstoredt = currentCohort%storage_flux

    if ( DEBUG ) write(iulog,*) 'EDPhys dbstoredt I ',currentCohort%dbstoredt

    currentCohort%seed_prod = (1.0_r8 - gr_fract) * currentCohort%carbon_balance
    if (abs(currentCohort%npp-(currentCohort%dbalivedt+currentCohort%dbdeaddt+currentCohort%dbstoredt+ &
         currentCohort%seed_prod+currentCohort%md)) > 0.0000000001_r8)then
       write(iulog,*) 'error in carbon check growth derivs',currentCohort%npp- &
            (currentCohort%dbalivedt+currentCohort%dbdeaddt+currentCohort%dbstoredt+currentCohort%seed_prod+currentCohort%md)
       write(iulog,*) 'cohort fluxes',currentCohort%pft,currentCohort%canopy_layer,currentCohort%n, &
            currentCohort%npp,currentCohort%dbalivedt,balive_loss, &
            currentCohort%dbdeaddt,currentCohort%dbstoredt,currentCohort%seed_prod,currentCohort%md * &
            EDecophyscon%leaf_stor_priority(currentCohort%pft)
       write(iulog,*) 'proxies' ,target_balive,currentCohort%balive,currentCohort%dbh,va,vs,gr_fract
    endif

    ! prevent negative leaf pool (but not negative store pool). This is also a numerical error prevention, 
    ! but it shouldn't happen actually... 
    if (-1.0_r8*currentCohort%dbalivedt * udata%deltat > currentCohort%balive*0.99)then 
       write(iulog,*) 'using non-neg leaf mass cap',currentCohort%balive , currentCohort%dbalivedt,currentCohort%dbstoredt, &
            currentCohort%carbon_balance
       currentCohort%dbstoredt = currentCohort%dbstoredt + currentCohort%dbalivedt

       if ( DEBUG ) write(iulog,*) 'EDPhys dbstoredt II ',currentCohort%dbstoredt

       currentCohort%dbalivedt = 0._r8 
    endif

    currentCohort%npp_bseed = currentCohort%seed_prod
    currentCohort%npp_store = max(0.0_r8,currentCohort%storage_flux)

    ! calculate change in diameter and height 
    currentCohort%ddbhdt = currentCohort%dbdeaddt * dDbhdBd(currentCohort)
    currentCohort%dhdt   = currentCohort%dbdeaddt * dHdBd(currentCohort)

    ! If the cohort has grown, it is not new
    currentCohort%isnew=.false.

  end subroutine Growth_Derivatives

  ! ============================================================================
  subroutine recruitment( t, currentPatch )
    !
    ! !DESCRIPTION:
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    use EDGrowthFunctionsMod, only : bdead,dbh, Bleaf
    use EDTypesMod, only : udata
    !
    ! !ARGUMENTS    
    integer, intent(in) :: t
    type(ed_patch_type), intent(inout), pointer :: currentPatch
    !
    ! !LOCAL VARIABLES:
    integer :: ft
    type (ed_cohort_type) , pointer :: temp_cohort
    integer :: cohortstatus
    !----------------------------------------------------------------------

    allocate(temp_cohort) ! create temporary cohort
    call zero_cohort(temp_cohort)

    do ft = 1,numpft_ed

       temp_cohort%canopy_trim = 0.8_r8  !starting with the canopy not fully expanded 
       temp_cohort%pft         = ft
       temp_cohort%hite        = EDecophyscon%hgt_min(ft)
       temp_cohort%dbh         = Dbh(temp_cohort)
       temp_cohort%bdead       = Bdead(temp_cohort)
       temp_cohort%balive      = Bleaf(temp_cohort)*(1.0_r8 + pftcon%froot_leaf(ft) &
            + EDecophyscon%sapwood_ratio(ft)*temp_cohort%hite)
       temp_cohort%bstore      = EDecophyscon%cushion(ft)*(temp_cohort%balive/ (1.0_r8 + pftcon%froot_leaf(ft) &
            + EDecophyscon%sapwood_ratio(ft)*temp_cohort%hite))
       temp_cohort%n           = currentPatch%area * currentPatch%seed_germination(ft)*udata%deltat &
            / (temp_cohort%bdead+temp_cohort%balive+temp_cohort%bstore)
 
       if (t == 1)then
          write(iulog,*) 'filling in cohorts where there are none left; this will break carbon balance', &
               currentPatch%patchno,currentPatch%area
          temp_cohort%n = 0.1_r8*currentPatch%area
          write(iulog,*) 'cohort n',ft,temp_cohort%n
       endif

       temp_cohort%laimemory = 0.0_r8     
       if (pftcon%season_decid(temp_cohort%pft) == 1.and.currentPatch%siteptr%status == 1)then
         temp_cohort%laimemory = (1.0_r8/(1.0_r8 + pftcon%froot_leaf(ft) + &
              EDecophyscon%sapwood_ratio(ft)*temp_cohort%hite))*temp_cohort%balive
       endif
       if (pftcon%stress_decid(temp_cohort%pft) == 1.and.currentPatch%siteptr%dstatus == 1)then
         temp_cohort%laimemory = (1.0_r8/(1.0_r8 + pftcon%froot_leaf(ft) + &
            EDecophyscon%sapwood_ratio(ft)*temp_cohort%hite))*temp_cohort%balive
       endif

       cohortstatus = currentPatch%siteptr%status
       if (pftcon%stress_decid(ft) == 1)then !drought decidous, override status. 
          cohortstatus = currentPatch%siteptr%dstatus
       endif

       if (temp_cohort%n > 0.0_r8 )then
           if ( DEBUG ) write(iulog,*) 'EDPhysiologyMod.F90 call create_cohort '
           call create_cohort(currentPatch, temp_cohort%pft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
                temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore,  &
                temp_cohort%laimemory, cohortstatus, temp_cohort%canopy_trim, currentPatch%NCL_p)
       endif

    enddo  !pft loop

    deallocate(temp_cohort) ! delete temporary cohort

  end subroutine recruitment

  ! ============================================================================
  subroutine CWD_Input( currentPatch)
    !
    ! !DESCRIPTION:
    ! Generate litter fields from turnover.  
    !
    ! !USES:
    use SFParamsMod , only : SF_val_CWD_frac
    use EDParamsMod , only : ED_val_ag_biomass
    use EDTypesMod  , only : udata
    !
    ! !ARGUMENTS    
    type(ed_patch_type),intent(inout), target :: currentPatch
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: currentCohort
    integer  :: c,p
    real(r8) :: not_dead_n !projected remaining number of trees in understorey cohort after turnover
    real(r8) :: dead_n !understorey dead tree density
    integer  :: pft
    !----------------------------------------------------------------------

    ! ================================================        
    ! Other direct litter fluxes happen in phenology and in spawn_patches. 
    ! ================================================   

    currentCohort => currentPatch%shortest

    do while(associated(currentCohort))
      pft = currentCohort%pft        
      ! ================================================        
      ! Litter from tissue turnover. KgC/m2/year
      ! ================================================   
      currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
               currentCohort%leaf_md * currentCohort%n/currentPatch%area !turnover

      currentPatch%root_litter_in(pft) = currentPatch%root_litter_in(pft) + &
               currentCohort%root_md * currentCohort%n/currentPatch%area !turnover
      currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
         currentCohort%leaf_litter * currentCohort%n/currentPatch%area/udata%deltat

      !daily leaf loss needs to be scaled up to the annual scale here. 
      
      do c = 1,ncwd
         currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + currentCohort%woody_turnover * &
              SF_val_CWD_frac(c) * currentCohort%n/currentPatch%area *ED_val_ag_biomass
         currentPatch%cwd_BG_in(c) = currentPatch%cwd_BG_in(c) + currentCohort%woody_turnover * &
              SF_val_CWD_frac(c) * currentCohort%n/currentPatch%area *(1.0_r8-ED_val_ag_biomass)
      enddo

      if (currentCohort%canopy_layer > 1)then   

          ! ================================================        
          ! Litter fluxes for understorey  mortality. KgC/m2/year
          ! ================================================
          dead_n = -1.0_r8 * currentCohort%dndt / currentPatch%area

          currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
               (currentCohort%bl+currentCohort%leaf_litter/udata%deltat)* dead_n          
          currentPatch%root_litter_in(pft) = currentPatch%root_litter_in(pft) + &
               (currentCohort%br+currentCohort%bstore)     * dead_n

          do c = 1,ncwd
             currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + (currentCohort%bdead+currentCohort%bsw) * &
                   SF_val_CWD_frac(c) * dead_n * ED_val_ag_biomass
             currentPatch%cwd_BG_in(c) = currentPatch%cwd_BG_in(c) + (currentCohort%bdead+currentCohort%bsw) * &
                  SF_val_CWD_frac(c) * dead_n * (1.0_r8-ED_val_ag_biomass)

             if (currentPatch%cwd_AG_in(c) < 0.0_r8)then
                write(iulog,*) 'negative CWD in flux',currentPatch%cwd_AG_in(c), &
                     (currentCohort%bdead+currentCohort%bsw), dead_n
             endif
          enddo

       endif !canopy layer

       currentCohort => currentCohort%taller

    enddo  ! end loop over cohorts 

    do p = 1,numpft_ed
       currentPatch%leaf_litter_in(p) = currentPatch%leaf_litter_in(p) + currentPatch%seed_decay(p) !KgC/m2/yr
    enddo

  end subroutine CWD_Input

  ! ============================================================================
  subroutine fragmentation_scaler( currentPatch, temperature_inst )
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! FIX(SPM, 091914) this should be a function as it returns a value in currentPatch%fragmentation_scaler
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_PI, SHR_CONST_TKFRZ
    use EDSharedParamsMod  , only : EDParamsShareInst

    !
    ! !ARGUMENTS    
    type(ed_patch_type)    , intent(inout) :: currentPatch
    type(temperature_type) , intent(in)    :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    logical  :: use_century_tfunc = .false.
    type(ed_site_type), pointer :: currentSite
    integer  :: p,j
    real(r8) :: t_scalar
    real(r8) :: w_scalar
    real(r8) :: catanf                ! hyperbolic temperature function from CENTURY
    real(r8) :: catanf_30             ! hyperbolic temperature function from CENTURY
    real(r8) :: t1                    ! temperature argument
    real(r8) :: Q10                   ! temperature dependence
    real(r8) :: froz_q10              ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
    real(r8), pointer :: t_veg24(:)
    !----------------------------------------------------------------------

    catanf(t1) = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))

    t_veg24 => temperature_inst%t_veg24_patch      ! Input:  [real(r8) (:)]  avg pft vegetation temperature for last 24 hrs

    catanf_30 = catanf(30._r8)
    
!    c = currentPatch%siteptr%clmcolumn
    p = currentPatch%clm_pno
    
    ! set "froz_q10" parameter
    froz_q10  = EDParamsShareInst%froz_q10  
    Q10       = EDParamsShareInst%Q10

    if ( .not. use_century_tfunc ) then
    !calculate rate constant scalar for soil temperature,assuming that the base rate constants 
    !are assigned for non-moisture limiting conditions at 25C. 
      if (t_veg24(p)  >=  SHR_CONST_TKFRZ) then
        t_scalar = Q10**((t_veg24(p)-(SHR_CONST_TKFRZ+25._r8))/10._r8)
                 !  Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8)
      else
        t_scalar = (Q10**(-25._r8/10._r8))*(froz_q10**((t_veg24(p)-SHR_CONST_TKFRZ)/10._r8))
                  !Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8)
      endif
    else
      ! original century uses an arctangent function to calculate the temperature dependence of decomposition      
      t_scalar = max(catanf(t_veg24(p)-SHR_CONST_TKFRZ)/catanf_30,0.01_r8)
    endif    
   
    !Moisture Limitations   
    !BTRAN APPROACH - is quite simple, but max's out decomp at all unstressed soil moisture values, which is not realistic.  
    !litter decomp is proportional to water limitation on average... 
    w_scalar = sum(currentPatch%btran_ft(1:numpft_ed))/numpft_ed 

    currentPatch%fragmentation_scaler =  min(1.0_r8,max(0.0_r8,t_scalar * w_scalar))
    
  end subroutine fragmentation_scaler
  
  ! ============================================================================
  subroutine cwd_out( currentPatch, temperature_inst )
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    use SFParamsMod, only : SF_val_max_decomp
    use EDTypesMod , only : udata
    !
    ! !ARGUMENTS    
    type(ed_patch_type)    , intent(inout), target :: currentPatch
    type(temperature_type) , intent(in)            :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_site_type), pointer :: currentSite
    integer :: c,ft
    !----------------------------------------------------------------------

    currentSite => currentPatch%siteptr
    currentPatch%root_litter_out(:) = 0.0_r8
    currentPatch%leaf_litter_out(:) = 0.0_r8

    call fragmentation_scaler(currentPatch, temperature_inst)

    !Flux of coarse woody debris into decomposing litter pool. 

    currentPatch%cwd_ag_out(1:ncwd) = 0.0_r8
    currentPatch%cwd_bg_out(1:ncwd) = 0.0_r8
    currentPatch%leaf_litter_out(1:numpft_ed) = 0.0_r8
    currentPatch%root_litter_out(1:numpft_ed) = 0.0_r8
    
    do c = 1,ncwd  
       currentPatch%cwd_ag_out(c)      = max(0.0_r8,   currentPatch%cwd_ag(c) * &
            SF_val_max_decomp(c+1) * currentPatch%fragmentation_scaler )  
       currentPatch%cwd_bg_out(c)      = max(0.0_r8,   currentPatch%cwd_bg(c) * &
            SF_val_max_decomp(c+1) * currentPatch%fragmentation_scaler )
    enddo

    ! this is the rate at which dropped leaves stop being part of the burnable pool and begin to be part of the 
    ! decomposing pool. This should probably be highly sensitive to moisture, but also to the type of leaf 
    ! thick leaves can dry out before they are decomposed, for example. 
    ! this section needs further scientific input. 

    do ft = 1,numpft_ed
       currentPatch%leaf_litter_out(ft) = max(0.0_r8,currentPatch%leaf_litter(ft)* SF_val_max_decomp(dg_sf) * &
            currentPatch%fragmentation_scaler )
       currentPatch%root_litter_out(ft) = max(0.0_r8,currentPatch%root_litter(ft)* SF_val_max_decomp(dg_sf) * &
            currentPatch%fragmentation_scaler )
       if ( currentPatch%leaf_litter_out(ft)<0.0_r8.or.currentPatch%root_litter_out(ft)<0.0_r8)then
         write(iulog,*) 'root or leaf out is negative?',SF_val_max_decomp(dg_sf),currentPatch%fragmentation_scaler
       endif
    enddo

    !add up carbon going into fragmenting pools
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%leaf_litter_out) * &
         currentPatch%area *udata%deltat!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%root_litter_out) * &
         currentPatch%area *udata%deltat!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%cwd_ag_out) * &
         currentPatch%area *udata%deltat!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%cwd_bg_out) * &
         currentPatch%area *udata%deltat!kgC/site/day

  end subroutine cwd_out



  subroutine flux_into_litter_pools(nsites, sites, bc_in, bc_out)
    ! Created by Charlie Koven and Rosie Fisher, 2014-2015
    ! take the flux out of the fragmenting litter pools and port into the decomposing litter pools. 
    ! in this implementation, decomposing pools are assumed to be humus and non-flammable, whereas fragmenting pools
    ! are assumed to be physically fragmenting but not respiring. This is a simplification, but allows us to 
    ! a) reconcile the need to track both chemical fractions (lignin, cellulose, labile) and size fractions (trunk, branch, etc.)
    ! b) to impose a realistic delay on the surge of nutrients into the litter pools when large CWD is added to the system via mortality
    
    ! because of the different subgrid structure, this subroutine includes the functionality that in the big-leaf BGC model, is calculated in SoilBiogeochemVerticalProfileMod
    
    ! The ED code is resolved at a daily timestep, but all of the CN-BGC fluxes are passed in as derivatives per second, 
    ! and then accumulated in the CNStateUpdate routines. One way of doing this is to pass back the CN fluxes per second, 
    ! and keep them constant for the whole day (making sure they are not overwritten.
    ! This means that the carbon gets passed back and forth between the photosynthesis code (fast timestepping) to the ED code (slow timestepping), back to the BGC code (fast timestepping).
    ! This means that the state update for the litter pools and for the CWD pools occurs at different timescales. 
    

    use EDTypesMod, only : AREA, numpft_ed, cp_numlevdecomp_full, cp_numlevdecomp
    use SoilBiogeochemVerticalProfileMod, only: surfprof_exp
    use EDCLMLinkMod, only: cwd_fcel_ed, cwd_flig_ed
    use pftconMod, only : pftcon
    use shr_const_mod, only: SHR_CONST_CDAY
    use clm_varcon, only : zisoi, dzsoi_decomp, zsoi
    use EDParamsMod, only : ED_val_ag_biomass
    use FatesInterfaceMod, only : bc_in_type, bc_out_type
    use clm_varctl, only : use_vertsoilc
    use abortutils  , only : endrun

    ! INTERF-TODO: remove the control parameters: exponential_rooting_profile, pftspecific_rootingprofile, rootprof_exp, surfprof_exp, zisoi, dzsoi_decomp, zsoi
    !
    implicit none   
    !
    ! !ARGUMENTS    
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(:)
    type(bc_out_type)       , intent(inout)           :: bc_out(:)
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type)  , pointer :: currentPatch
    type (ed_cohort_type) , pointer :: currentCohort
    type(ed_site_type), pointer :: cs
    integer p,ci,j,s
    real(r8) time_convert    ! from year to seconds
    real(r8) mass_convert    ! ED uses kg, CLM uses g
    integer           :: begp,endp
    integer           :: begc,endc                                    !bounds 
    !------------------------------------------------------------------------
    real(r8) :: cinput_rootfr(1:numpft_ed, 1:cp_numlevdecomp_full)      ! column by pft root fraction used for calculating inputs
    real(r8) :: croot_prof_perpatch(1:cp_numlevdecomp_full)
    real(r8) :: surface_prof(1:cp_numlevdecomp_full)
    integer  :: ft
    real(r8) :: rootfr_tot(1:numpft_ed), biomass_bg_ft(1:numpft_ed)
    real(r8) :: surface_prof_tot, leaf_prof_sum, stem_prof_sum, froot_prof_sum, biomass_bg_tot
    real(r8) :: delta

    ! NOTE(bja, 201608) these were removed from clm in clm4_5_10_r187
    logical, parameter :: exponential_rooting_profile = .true.
    logical, parameter :: pftspecific_rootingprofile = .true.

    ! NOTE(bja, 201608) as of clm4_5_10_r187 rootprof_exp is now a
    ! private function level parameter in RootBiophysMod.F90::exponential_rootfr()
    real(r8), parameter :: rootprof_exp  = 3.  ! how steep profile is
    ! for root C inputs (1/ e-folding depth) (1/m)

    ! NOTE(bja, 201608) as of clm4_5_10_r187 rootprof_beta is now a
    ! two dimensional array with the second dimension being water,1,
    ! or carbon,2,. These are currently hard coded, but may be
    ! overwritten by the namelist.

    ! Note cdk 2016/08 we actually want to use the carbon index here rather than the water index.  
    ! Doing so will be answer changing though so perhaps easiest to do this in steps.
    integer, parameter :: rooting_profile_varindex_water = 1

    real(r8) :: leaf_prof(1:nsites, 1:cp_numlevdecomp)
    real(r8) :: froot_prof(1:nsites,  1:numpft_ed, 1:cp_numlevdecomp)
    real(r8) :: croot_prof(1:nsites, 1:cp_numlevdecomp)
    real(r8) :: stem_prof(1:nsites, 1:cp_numlevdecomp)
    
    
    delta = 0.001_r8    
    !no of seconds in a year. 
    time_convert =  365.0_r8*SHR_CONST_CDAY

    ! number of grams in a kilogram
    mass_convert = 1000._r8
    
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! first calculate vertical profiles
      ! define two types of profiles: 
      ! (1) a surface profile, for leaves and stem inputs, which is the same for each pft but differs from one site to the next to avoid inputting any C into permafrost or bedrock
      ! (2) a fine root profile, which is indexed by both site and pft, differs for each pft and also from one site to the next to avoid inputting any C into permafrost or bedrock
      ! (3) a coarse root profile, which is the root-biomass=weighted average of the fine root profiles
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (use_vertsoilc) then
         
         ! define a single shallow surface profile for surface additions (leaves, stems, and N deposition)
         surface_prof(:) = 0._r8
         do j = 1, cp_numlevdecomp
            surface_prof(j) = exp(-surfprof_exp * zsoi(j)) / dzsoi_decomp(j)
         end do
         
         ! initialize profiles to zero
         leaf_prof(1:nsites, :)      = 0._r8
         froot_prof(1:nsites, 1:numpft_ed, :)     = 0._r8
         croot_prof(1:nsites, :)     = 0._r8
         stem_prof(1:nsites, :)      = 0._r8
         
         cinput_rootfr(1:numpft_ed, :)     = 0._r8
            
         ! calculate pft-specific rooting profiles in the absence of permafrost or bedrock limitations
         if ( exponential_rooting_profile ) then
            if ( .not. pftspecific_rootingprofile ) then
               ! define rooting profile from exponential parameters
               do ft = 1, numpft_ed
                  do j = 1, cp_numlevdecomp
                     cinput_rootfr(ft,j) = exp(-rootprof_exp * zsoi(j)) / dzsoi_decomp(j)
                  end do
               end do
            else
               ! use beta distribution parameter from Jackson et al., 1996
               do ft = 1, numpft_ed
                  do j = 1, cp_numlevdecomp
                     cinput_rootfr(ft,j) = ( pftcon%rootprof_beta(ft, rooting_profile_varindex_water) ** (zisoi(j-1)*100._r8) - &
                          pftcon%rootprof_beta(ft, rooting_profile_varindex_water) ** (zisoi(j)*100._r8) ) &
                          / dzsoi_decomp(j)
                  end do
               end do
            endif
         else
            do ft = 1,numpft_ed 
               do j = 1, cp_numlevdecomp
                  ! use standard CLM root fraction profiles;
                  cinput_rootfr(ft,j) =  ( .5_r8*( &
                       exp(-pftcon%roota_par(ft) * zisoi(j-1))  &
                       + exp(-pftcon%rootb_par(ft) * zisoi(j-1))  &
                       - exp(-pftcon%roota_par(ft) * zisoi(j))    &
                       - exp(-pftcon%rootb_par(ft) * zisoi(j))))  / dzsoi_decomp(j)
               end do
            end do
         endif
         !

         do s = 1,nsites
            !
            ! now add permafrost constraint: integrate rootfr over active layer of soil site,
            ! truncate below permafrost or bedrock table where present, and rescale so that integral = 1
            do ft = 1,numpft_ed 
               rootfr_tot(ft) = 0._r8
            end do
            surface_prof_tot = 0._r8
            !
            do j = 1, min(max(bc_in(s)%max_rooting_depth_index_col, 1), cp_numlevdecomp)
               surface_prof_tot = surface_prof_tot + surface_prof(j)  * dzsoi_decomp(j)
            end do
            do ft = 1,numpft_ed 
               do j = 1, min(max(bc_in(s)%max_rooting_depth_index_col, 1), cp_numlevdecomp)
                  rootfr_tot(ft) = rootfr_tot(ft) + cinput_rootfr(ft,j) * dzsoi_decomp(j)
               end do
            end do
            !
            ! rescale the fine root profile
            do ft = 1,numpft_ed 
               if ( (bc_in(s)%max_rooting_depth_index_col > 0) .and. (rootfr_tot(ft) > 0._r8) ) then
                  ! where there is not permafrost extending to the surface, integrate the profiles over the active layer
                  ! this is equivalent to integrating over all soil layers outside of permafrost regions
                  do j = 1, min(max(bc_in(s)%max_rooting_depth_index_col, 1), cp_numlevdecomp)
                     froot_prof(s,ft,j) = cinput_rootfr(ft,j) / rootfr_tot(ft)
                  end do
               else
                  ! if fully frozen, or no roots, put everything in the top layer
                  froot_prof(s,ft,1) = 1._r8/dzsoi_decomp(1)
               endif
            end do
            !
            ! rescale the shallow profiles
            if ( (bc_in(s)%max_rooting_depth_index_col > 0) .and. (surface_prof_tot > 0._r8) ) then
               ! where there is not permafrost extending to the surface, integrate the profiles over the active layer
               ! this is equivalent to integrating over all soil layers outside of permafrost regions
               do j = 1, min(max(bc_in(s)%max_rooting_depth_index_col, 1), cp_numlevdecomp)
                  ! set all surface processes to shallower profile
                  leaf_prof(s,j) = surface_prof(j)/ surface_prof_tot
                  stem_prof(s,j) = surface_prof(j)/ surface_prof_tot
               end do
            else
               ! if fully frozen, or no roots, put everything in the top layer
               leaf_prof(s,1) = 1._r8/dzsoi_decomp(1)
               stem_prof(s,1) = 1._r8/dzsoi_decomp(1)
               do j = 2, cp_numlevdecomp
                  leaf_prof(s,j) = 0._r8
                  stem_prof(s,j) = 0._r8
               end do
            endif
         end do
         
      else
         
         ! for one layer decomposition model, set profiles to unity
         leaf_prof(1:nsites, :) = 1._r8
         froot_prof(1:nsites, 1:numpft_ed, :) = 1._r8
         stem_prof(1:nsites, :) = 1._r8
         
      end if
      
      ! sanity check to ensure they integrate to 1
      do s = 1, nsites
         ! check the leaf and stem profiles
         leaf_prof_sum = 0._r8
         stem_prof_sum = 0._r8
         do j = 1, cp_numlevdecomp
            leaf_prof_sum = leaf_prof_sum + leaf_prof(s,j) *  dzsoi_decomp(j)
            stem_prof_sum = stem_prof_sum + stem_prof(s,j) *  dzsoi_decomp(j)
         end do
         if ( ( abs(stem_prof_sum - 1._r8) > delta ) .or.  ( abs(leaf_prof_sum - 1._r8) > delta ) ) then
            write(iulog, *) 'profile sums: ',  leaf_prof_sum, stem_prof_sum
            write(iulog, *) 'surface_prof: ', surface_prof
            write(iulog, *) 'surface_prof_tot: ', surface_prof_tot
            write(iulog, *) 'leaf_prof: ',  leaf_prof(s,:)
            write(iulog, *) 'stem_prof: ',  stem_prof(s,:)
            write(iulog, *) 'max_rooting_depth_index_col: ', bc_in(s)%max_rooting_depth_index_col
            write(iulog, *) 'dzsoi_decomp: ',  dzsoi_decomp            
            call endrun()
         endif
         ! now check each fine root profile
         do ft = 1,numpft_ed 
            froot_prof_sum = 0._r8
            do j = 1, cp_numlevdecomp
               froot_prof_sum = froot_prof_sum + froot_prof(s,ft,j) *  dzsoi_decomp(j)
            end do
            if ( ( abs(froot_prof_sum - 1._r8) > delta ) ) then
               write(iulog, *) 'profile sums: ', froot_prof_sum
               call endrun()
            endif
         end do
      end do
      
      ! zero the site-level C input variables
      do s = 1, nsites
         do j = 1, cp_numlevdecomp
            bc_out(s)%FATES_c_to_litr_lab_c_col(j) = 0._r8
            bc_out(s)%FATES_c_to_litr_cel_c_col(j) = 0._r8
            bc_out(s)%FATES_c_to_litr_lig_c_col(j) = 0._r8
            croot_prof(s,j)         = 0._r8
         end do
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! now disaggregate the inputs vertically, using the vertical profiles
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do s = 1,nsites
         
         !      do g = bounds%begg,bounds%endg
         !         if (firstsoilpatch(g) >= 0 .and. ed_allsites_inst(g)%istheresoil) then 
         currentPatch => sites(s)%oldest_patch
         
         do while(associated(currentPatch))
            
            ! the CWD pools lose information about which PFT they came from; for the stems this doesn't matter as they all have the same profile, 
            ! however for the coarse roots they may have different profiles.  to approximately recover this information, loop over all cohorts in patch 
            ! to calculate the total root biomass in that patch of each pft, and then rescale the croot_prof as the weighted average of the froot_prof
            biomass_bg_ft(1:numpft_ed) = 0._r8
            currentCohort => currentPatch%tallest
            do while(associated(currentCohort))      
               biomass_bg_ft(currentCohort%pft) = biomass_bg_ft(currentCohort%pft) + &
                    currentCohort%b * (currentCohort%n / currentPatch%area) * (1.0_r8-ED_val_ag_biomass)
               currentCohort => currentCohort%shorter
            enddo !currentCohort
            ! 
            biomass_bg_tot = 0._r8
            do ft = 1,numpft_ed 
               biomass_bg_tot = biomass_bg_tot + biomass_bg_ft(ft)
            end do
            !         
            do j = 1, cp_numlevdecomp
               ! zero this for each patch
               croot_prof_perpatch(j) = 0._r8
            end do
            !
            if ( biomass_bg_tot .gt. 0._r8) then
               do ft = 1,numpft_ed 
                  do j = 1, cp_numlevdecomp
                     croot_prof_perpatch(j) = croot_prof_perpatch(j) + froot_prof(s,ft,j) * biomass_bg_ft(ft) / biomass_bg_tot
                  end do
               end do
            else ! no biomass
               croot_prof_perpatch(1) = 1./dzsoi_decomp(1)
            end if

            !
            ! add croot_prof as weighted average (weighted by patch area) of croot_prof_perpatch
            do j = 1, cp_numlevdecomp
               croot_prof(s, j) = croot_prof(s, j) + croot_prof_perpatch(j) * currentPatch%area / AREA
            end do
            !
            ! now disaggregate, vertically and by decomposition substrate type, the actual fluxes from CWD and litter pools
            !
            ! do c = 1, ncwd
            !    write(iulog,*)'cdk CWD_AG_out', c, currentpatch%CWD_AG_out(c), cwd_fcel_ed, currentpatch%area/AREA
            !    write(iulog,*)'cdk CWD_BG_out', c, currentpatch%CWD_BG_out(c), cwd_fcel_ed, currentpatch%area/AREA
            ! end do
            ! do ft = 1,numpft_ed
            !    write(iulog,*)'cdk leaf_litter_out', ft, currentpatch%leaf_litter_out(ft), cwd_fcel_ed, currentpatch%area/AREA
            !    write(iulog,*)'cdk root_litter_out', ft, currentpatch%root_litter_out(ft), cwd_fcel_ed, currentpatch%area/AREA
            ! end do
            ! !
            ! CWD pools fragmenting into decomposing litter pools. 
            do ci = 1, ncwd
               do j = 1, cp_numlevdecomp
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%CWD_AG_out(ci) * cwd_fcel_ed * currentpatch%area/AREA * stem_prof(s,j)  
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%CWD_AG_out(ci) * cwd_flig_ed * currentpatch%area/AREA * stem_prof(s,j)
                  !
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%CWD_BG_out(ci) * cwd_fcel_ed * currentpatch%area/AREA * croot_prof_perpatch(j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%CWD_BG_out(ci) * cwd_flig_ed * currentpatch%area/AREA * croot_prof_perpatch(j)
               end do
            end do
            
            ! leaf and fine root pools. 
            do ft = 1,numpft_ed
               do j = 1, cp_numlevdecomp
                  bc_out(s)%FATES_c_to_litr_lab_c_col(j) = bc_out(s)%FATES_c_to_litr_lab_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * pftcon%lf_flab(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * pftcon%lf_fcel(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * pftcon%lf_flig(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  !
                  bc_out(s)%FATES_c_to_litr_lab_c_col(j) = bc_out(s)%FATES_c_to_litr_lab_c_col(j) + &
                       currentpatch%root_litter_out(ft) * pftcon%fr_flab(ft) * currentpatch%area/AREA * froot_prof(s,ft,j)
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%root_litter_out(ft) * pftcon%fr_fcel(ft) * currentpatch%area/AREA * froot_prof(s,ft,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%root_litter_out(ft) * pftcon%fr_flig(ft) * currentpatch%area/AREA * froot_prof(s,ft,j)
                  !
                  !! and seed_decay too.  for now, use the same lability fractions as for leaf litter
                  bc_out(s)%FATES_c_to_litr_lab_c_col(j) = bc_out(s)%FATES_c_to_litr_lab_c_col(j) + &
                       currentpatch%seed_decay(ft) * pftcon%lf_flab(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%seed_decay(ft) * pftcon%lf_fcel(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%seed_decay(ft) * pftcon%lf_flig(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  !
               enddo
            end do
              
              currentPatch => currentPatch%younger
           end do !currentPatch

        end do  ! do sites(s)
     
        do s = 1, nsites
           do j = 1, cp_numlevdecomp                    
              ! time unit conversion
              bc_out(s)%FATES_c_to_litr_lab_c_col(j)=bc_out(s)%FATES_c_to_litr_lab_c_col(j) * mass_convert / time_convert
              bc_out(s)%FATES_c_to_litr_cel_c_col(j)=bc_out(s)%FATES_c_to_litr_cel_c_col(j) * mass_convert / time_convert
              bc_out(s)%FATES_c_to_litr_lig_c_col(j)=bc_out(s)%FATES_c_to_litr_lig_c_col(j) * mass_convert / time_convert
              
           end do
        end do
        
        ! write(iulog,*)'cdk FATES_c_to_litr_lab_c: ', FATES_c_to_litr_lab_c
        ! write_col(iulog,*)'cdk FATES_c_to_litr_cel_c: ', FATES_c_to_litr_cel_c    
        ! write_col(iulog,*)'cdk FATES_c_to_litr_lig_c: ', FATES_c_to_litr_lig_c
        ! write_col(iulog,*)'cdk cp_numlevdecomp_full,  bounds%begc, bounds%endc: ', cp_numlevdecomp_full, bounds%begc, bounds%endc
        ! write(iulog,*)'cdk leaf_prof: ', leaf_prof
        ! write(iulog,*)'cdk stem_prof: ', stem_prof    
        ! write(iulog,*)'cdk froot_prof: ', froot_prof
        ! write(iulog,*)'cdk croot_prof_perpatch: ', croot_prof_perpatch
        ! write(iulog,*)'cdk croot_prof: ', croot_prof

    end subroutine flux_into_litter_pools

end module EDPhysiologyMod
