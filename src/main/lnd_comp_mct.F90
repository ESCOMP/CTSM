#include <misc.h>
#include <preproc.h>

module lnd_comp_mct
  
#if (defined SEQ_MCT)

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: lnd_comp_mct
!
! !DESCRIPTION:
!
! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use mct_mod          , only : mct_aVect
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: lnd_init_mct
  public :: lnd_run_mct
  public :: lnd_final_mct
  SAVE
  private                              ! By default make data private
!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
! !PRIVATE MEMBER FUNCTIONS:
  private :: lnd_SetgsMap_mct
  private :: lnd_chkAerDep_mct
  private :: lnd_domain_mct
  private :: lnd_export_mct
  private :: lnd_import_mct
#ifdef RTM
  private :: rof_SetgsMap_mct
  private :: rof_domain_mct
  private :: rof_export_mct
#endif
!
! !PRIVATE VARIABLES
!
! Time averaged flux fields
!  
  type(mct_aVect)   :: l2x_l_SNAP
  type(mct_aVect)   :: l2x_l_SUM
!
! Time averaged counter for flux fields
!
  integer :: avg_count
!
! Atmospheric mode  
!
  logical :: atm_prognostic

!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_init_mct
!
! !INTERFACE:
  subroutine lnd_init_mct( EClock, cdata_l, x2l_l, l2x_l, &
                                   cdata_r, r2x_r,        &
                                   cdata_s, x2s_s, s2x_s, &
                                   NLFilename )
!
! !DESCRIPTION:
! Initialize land surface model and obtain relevant atmospheric model arrays
! back from (i.e. albedos, surface temperature and snow cover over land).
!
! !USES:
    use shr_kind_mod     , only : r8 => shr_kind_r8
    use clm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, &
                                  set_nextsw_cday
    use clm_atmlnd       , only : clm_mapl2a, clm_l2a, atm_l2a
    use clm_comp         , only : clm_init0, clm_init1, clm_init2
    use clm_varctl       , only : finidat,single_column, set_clmvarctl
    use controlMod       , only : control_setNL
    use domainMod        , only : amask, adomain
    use clm_varpar       , only : rtmlon, rtmlat
    use clm_varorb       , only : eccen, obliqr, lambm0, mvelpp
    use abortutils       , only : endrun
    use esmf_mod         , only : ESMF_Clock
    use clm_varctl       , only : iulog
    use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                  shr_file_getLogUnit, shr_file_getLogLevel, &
                                  shr_file_getUnit, shr_file_setIO
    use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
    use spmdMod          , only : masterproc, spmd_init
    use seq_timemgr_mod  , only : seq_timemgr_EClockGetData
    use seq_infodata_mod , only : seq_infodata_type, seq_infodata_GetData, seq_infodata_PutData, &
                                  seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                  seq_infodata_start_type_brnch
    use mct_mod          , only : mct_aVect, mct_gsMap, mct_gGrid, mct_aVect_init, mct_aVect_zero
    use seq_flds_mod
    use seq_flds_indices
!
! !ARGUMENTS:
    type(ESMF_Clock),             intent(in)    :: EClock
    type(seq_cdata),              intent(inout) :: cdata_l
    type(mct_aVect),              intent(inout) :: x2l_l, l2x_l
    type(seq_cdata),              intent(inout) :: cdata_r
    type(mct_aVect),              intent(inout) ::        r2x_r
    type(seq_cdata),  optional,   intent(inout) :: cdata_s
    type(mct_aVect),  optional,   intent(inout) :: x2s_s, s2x_s
    character(len=*), optional,   intent(in)    :: NLFilename 
!
! !LOCAL VARIABLES:
    integer                                     :: LNDID	
    integer                                     :: mpicom_lnd       	
    type(mct_gsMap),              pointer       :: GSMap_lnd
    type(mct_gGrid),              pointer       :: dom_l
    type(mct_gsMap),              pointer       :: GSMap_rof
    type(mct_gGrid),              pointer       :: dom_r
    type(seq_infodata_type),      pointer       :: infodata
    integer  :: lsize           ! size of attribute vector
    integer  :: i,j             ! indices
    integer  :: dtime_sync
    integer  :: dtime_clm
    logical  :: exists               ! true if file exists
    real(r8) :: scmlat
    real(r8) :: scmlon
    real(r8) :: nextsw_cday     ! calday from clock of next radiation computation
    character(len=SHR_KIND_CL) :: caseid
    character(len=SHR_KIND_CL) :: ctitle
    character(len=SHR_KIND_CL) :: starttype
    character(len=SHR_KIND_CL) :: calendar
    character(len=SHR_KIND_CL) :: hostname     ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version      ! Model version
    character(len=SHR_KIND_CL) :: username     ! user running the model
    integer  :: nsrest
    integer :: startype
    integer :: perpetual_ymd
    integer :: ref_ymd
    integer :: ref_tod
    integer :: start_ymd
    integer :: start_tod
    integer :: stop_ymd
    integer :: stop_tod
    logical :: brnch_retain_casename
    logical :: perpetual_run
    integer :: lbnum
    integer  :: shrlogunit,shrloglev ! old values
    character(len=32), parameter :: sub = 'lnd_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------
    ! Set cdata data

    call seq_cdata_setptrs(cdata_l, ID=LNDID, mpicom=mpicom_lnd, &
         gsMap=GSMap_lnd, dom=dom_l, infodata=infodata)

    call seq_cdata_setptrs(cdata_r, &
         gsMap=gsMap_rof, dom=dom_r) 

    ! Initialize clm MPI communicator 

    call spmd_init( mpicom_lnd )

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_init_mct:start::',lbnum)
    endif
#endif                      

    ! Initialize io log unit

    if (masterproc) then
       inquire(file='lnd_modelio.nml',exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('lnd_modelio.nml',iulog)
       end if
       write(iulog,format) "CLM land model initialization"
    end if

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Use infodata to set orbital values

    call seq_infodata_GetData( infodata, orb_eccen=eccen, orb_mvelpp=mvelpp, &
         orb_lambm0=lambm0, orb_obliqr=obliqr )

    ! Consistency check on namelist filename	

    call control_setNL( 'lnd_in' )

    ! Initialize clm
    ! clm_init0 reads namelist, grid and surface data
    ! clm_init1 and clm_init2 performs rest of initialization	
    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )
    call seq_infodata_GetData(infodata, perpetual=perpetual_run,                &
                              perpetual_ymd=perpetual_ymd, case_name=caseid,    &
                              case_desc=ctitle, single_column=single_column,    &
                              scmlat=scmlat, scmlon=scmlon,                     &
                              brnch_retain_casename=brnch_retain_casename,      &
                              start_type=starttype, model_version=version,      &
                              hostname=hostname, username=username              &
                                )
    call set_timemgr_init( calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
                           ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd,         &
                           stop_tod_in=stop_tod,  perpetual_run_in=perpetual_run,                &
                           perpetual_ymd_in=perpetual_ymd )
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = 0
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = 1
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = 3
    else
       call endrun( sub//' ERROR: unknown starttype' )
    end if

    call set_clmvarctl(    caseid_in=caseid, ctitle_in=ctitle,                     &
                           brnch_retain_casename_in=brnch_retain_casename,         &
                           single_column_in=single_column, scmlat_in=scmlat,       &
                           scmlon_in=scmlon, nsrest_in=nsrest, version_in=version, &
                           hostname_in=hostname, username_in=username )

    call clm_init0( )

    ! If in SCM mode and no land then exit out of initialization

    if ( single_column .and. amask(1)==0) then
       call seq_infodata_PutData(infodata, lnd_present=.false.)
       return
    end if

    call clm_init1( )
    call clm_init2()

    ! Check that clm internal dtime aligns with clm coupling interval

    call seq_timemgr_EClockGetData(EClock, dtime=dtime_sync )
    dtime_clm = get_step_size()
    if(masterproc) write(iulog,*)'dtime_sync= ',dtime_sync,' dtime_clm= ',dtime_clm,' mod = ',mod(dtime_sync,dtime_clm)
    if (mod(dtime_sync,dtime_clm) /= 0) then
       write(iulog,*)'clm dtime ',dtime_clm,' and Eclock dtime ',dtime_sync,' never align'
       call endrun( sub//' ERROR: time out of sync' )
    end if

    ! Initialize lnd gsMap

    call lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_lnd ) 	
    lsize = mct_gsMap_lsize(gsMap_lnd, mpicom_lnd)

    ! Initialize lnd domain

    call lnd_domain_mct( lsize, gsMap_lnd, dom_l )

    ! Initialize lnd attribute vectors

    call mct_aVect_init(x2l_l, rList=seq_flds_x2l_fields, lsize=lsize)
    call mct_aVect_zero(x2l_l)

    call mct_aVect_init(l2x_l, rList=seq_flds_l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x_l)

    call mct_aVect_init(l2x_l_SNAP, rList=seq_flds_l2x_fluxes, lsize=lsize)
    call mct_aVect_zero(l2x_l_SNAP)

    call mct_aVect_init(l2x_l_SUM , rList=seq_flds_l2x_fluxes, lsize=lsize)
    call mct_aVect_zero(l2x_l_SUM )

    if (masterproc) then
       write(iulog,format)'time averaging the following flux fields over the coupling interval'
       write(iulog,format) trim(seq_flds_l2x_fluxes)
    end if

    ! Create mct land export state

    call clm_mapl2a(clm_l2a, atm_l2a)
    call lnd_export_mct( atm_l2a, l2x_l )

#ifdef RTM
    ! Initialize rof gsMap

    call rof_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_rof ) 
    lsize = mct_gsMap_lsize(gsMap_rof, mpicom_lnd)

    ! Initialize rof domain

    call rof_domain_mct( lsize, gsMap_rof, dom_r )

    ! Initialize rtm attribute vectors		

    call mct_aVect_init(r2x_r, rList=seq_flds_r2x_fields, lsize=lsize)
    call mct_aVect_zero(r2x_r)

    ! Create mct river runoff export state

    call rof_export_mct( r2x_r )
#endif

    ! Initialize averaging counter

    avg_count = 0

    ! Set land modes

    call seq_infodata_PutData( infodata, lnd_prognostic=.true.)
    call seq_infodata_PutData( infodata, lnd_nx = adomain%ni, lnd_ny = adomain%nj)
#ifdef RTM
    call seq_infodata_PutData( infodata, rof_present=.true.)
    call seq_infodata_PutData( infodata, rof_nx = rtmlon, rof_ny = rtmlat)
#else
    call seq_infodata_PutData( infodata, rof_present=.false.)
#endif
    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    call set_nextsw_cday( nextsw_cday )

    ! Determine atmosphere modes

    call seq_infodata_GetData(infodata, atm_prognostic=atm_prognostic)
    if (masterproc) then
       if ( atm_prognostic )then
          write(iulog,format) 'Atmospheric input is from a prognostic model'
       else
          write(iulog,format) 'Atmospheric input is from a data model'
       end if
    end if

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_int_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine lnd_init_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_run_mct
!
! !INTERFACE:
  subroutine lnd_run_mct( EClock, cdata_l, x2l_l, l2x_l, &
                                  cdata_r,        r2x_r, &
                                  cdata_s, x2s_s, s2x_s )
!
! !DESCRIPTION:
! Run clm model
!
! !USES:
    use shr_kind_mod    ,only : r8 => shr_kind_r8
    use clm_atmlnd      ,only : clm_mapl2a, clm_mapa2l
    use clm_atmlnd      ,only : clm_l2a, atm_l2a, atm_a2l, clm_a2l
    use clm_comp        ,only : clm_run1, clm_run2
    use clm_time_manager,only : get_curr_date, get_nstep, get_curr_calday, get_step_size, &
                                advance_timestep, set_nextsw_cday
    use domainMod       ,only : adomain
    use decompMod       ,only : get_proc_bounds_atm
    use abortutils      ,only : endrun
    use esmf_mod        ,only : ESMF_Clock
    use clm_varctl      ,only : iulog
    use shr_file_mod    ,only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel
    use seq_cdata_mod   ,only : seq_cdata, seq_cdata_setptrs
    use seq_timemgr_mod ,only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use seq_infodata_mod,only : seq_infodata_type, seq_infodata_GetData
    use spmdMod         ,only : masterproc, mpicom
    use perf_mod        ,only : t_startf, t_stopf, t_barrierf
    use mct_mod         ,only : mct_aVect, mct_aVect_accum, mct_aVect_copy, mct_aVect_avg, &
                                mct_aVect_zero
    use mct_mod        , only : mct_gGrid, mct_gGrid_exportRAttr, mct_gGrid_lsize
    use aerdepMod      , only : aerdepini
!
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(in)    :: cdata_l
    type(mct_aVect)             , intent(inout) :: x2l_l
    type(mct_aVect)             , intent(inout) :: l2x_l
    type(seq_cdata)             , intent(in)    :: cdata_r
    type(mct_aVect)             , intent(inout) :: r2x_r
    type(seq_cdata),  optional,   intent(in)    :: cdata_s
    type(mct_aVect),  optional,   intent(inout) :: x2s_s
    type(mct_aVect),  optional,   intent(inout) :: s2x_s
!
! !LOCAL VARIABLES:
    integer :: ymd_sync        ! Sync date (YYYYMMDD)
    integer :: yr_sync         ! Sync current year
    integer :: mon_sync        ! Sync current month
    integer :: day_sync        ! Sync current day
    integer :: tod_sync        ! Sync current time of day (sec)
    integer :: ymd             ! CLM current date (YYYYMMDD)
    integer :: yr              ! CLM current year
    integer :: mon             ! CLM current month
    integer :: day             ! CLM current day
    integer :: tod             ! CLM current time of day (sec)
    integer :: dtime           ! time step increment (sec)
    integer :: nstep           ! time step index
    logical :: rstwr_sync      ! .true. ==> write restart file before returning
    logical :: rstwr           ! .true. ==> write restart file before returning
    logical :: nlend_sync      ! Flag signaling last time-step
    logical :: nlend           ! .true. ==> last time-step
    logical :: dosend          ! true => send data back to driver
    logical :: doalb           ! .true. ==> do albedo calculation on this time step
    real(r8):: nextsw_cday     ! calday from clock of next radiation computation
    real(r8):: caldayp1        ! clm calday plus dtime offset
    integer :: shrlogunit,shrloglev       ! old values
    integer :: begg, endg    
    integer :: lbnum
    type(seq_infodata_type),pointer :: infodata
    type(mct_gGrid),        pointer :: dom_l
    real(r8),               pointer :: data(:)  ! temporary
    integer :: g,i,lsize       ! counters
    logical,save :: first_call = .true.   ! first call work
    logical :: never_doAlb                ! if doalb never set
    character(len=32)            :: rdate ! date char string for restart file names
    character(len=32), parameter :: sub = "lnd_run_mct"
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_mct:start::',lbnum)
    endif
#endif

    ! Reset shr logging to my log file

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Determine time of next atmospheric shortwave calculation

    call seq_cdata_setptrs(cdata_l, infodata=infodata, dom=dom_l)
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)
    call seq_infodata_GetData(infodata, nextsw_cday=nextsw_cday )

    call set_nextsw_cday( nextsw_cday )

    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend_sync = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr_sync = seq_timemgr_RestartAlarmIsOn( EClock )

    lsize = mct_gGrid_lsize(dom_l)
    call get_proc_bounds_atm(begg, endg)
    if (first_call) then
       allocate(data(lsize))
       call mct_gGrid_exportRattr(dom_l,"ascale",data,lsize) 
       do g = begg,endg
          i = 1 + (g - begg)
           adomain%asca(g) = data(i)
       end do
       deallocate(data)

       call lnd_chkAerDep_mct( x2l_l )

       call aerdepini( )   ! Will be removed...

    endif

    
    ! Map MCT to land data type

    call t_startf ('lc_lnd_import')
    call lnd_import_mct( x2l_l, atm_a2l )
    call t_stopf ('lc_lnd_import')

    call t_startf ('lc_clm_mapa2l')
    call clm_mapa2l(atm_a2l, clm_a2l)
    call t_stopf ('lc_clm_mapa2l')
    
    ! Loop over time steps in coupling interval

    dosend      = .false.
    never_doAlb = .true.
    do while(.not. dosend)

       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated

       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))

       ! Determine doalb based on nextsw_cday sent from atm model

       dtime = get_step_size()
       caldayp1 = get_curr_calday(offset=dtime)
       doalb = abs(nextsw_cday- caldayp1) < 1.e-10_r8
       if ( doalb ) never_doAlb = .false.

       ! Determine if time to write cam restart and stop

       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.

       ! Run clm 

       call t_barrierf('sync_clm_run1', mpicom)
       call t_startf ('clm_run1')
       call clm_run1( doalb )
       call t_stopf ('clm_run1')

       call t_barrierf('sync_clm_run2', mpicom)
       call t_startf ('clm_run2')
       call clm_run2( rstwr, nlend, rdate )
       call t_stopf ('clm_run2')

       ! Map land data type to MCT
       
       call t_startf ('lc_clm_mapl2a')
       call clm_mapl2a(clm_l2a, atm_l2a)
       call t_stopf ('lc_clm_mapl2a')
       
       call t_startf ('lc_lnd_export')
       call lnd_export_mct( atm_l2a, l2x_l )
       call t_stopf ('lc_lnd_export')
       
       ! Compute snapshot attribute vector for accumulation
       
! don't accumulate on first coupling freq ts0 and ts1
! for consistency with ccsm3 when flxave is off
       nstep = get_nstep()
       if (nstep <= 1) then
          call mct_aVect_copy( l2x_l, l2x_l_SUM )
          avg_count = 1
       else
          call mct_aVect_copy( l2x_l, l2x_l_SNAP )
          call mct_aVect_accum( aVin=l2x_l_SNAP, aVout=l2x_l_SUM )
          avg_count = avg_count + 1
       endif
       
       ! Advance clm time step
       
       call t_startf ('lc_clm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_clm2_adv_timestep')

    end do

    ! Finish accumulation of attribute vector and average and zero out partial sum and counter
    
    call mct_aVect_avg ( l2x_l_SUM, avg_count)
    call mct_aVect_copy( l2x_l_SUM, l2x_l )
    call mct_aVect_zero( l2x_l_SUM) 
    avg_count = 0                   

#ifdef RTM
    ! Create river runoff output state

    call t_startf ('lc_rof_export')
    call rof_export_mct( r2x_r )
    call t_stopf ('lc_rof_export')
#endif
       
    ! Check that internal clock is in sync with master clock
    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' clm ymd=',ymd     ,'  clm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call endrun( sub//":: CLM clock not in sync with Master Sync clock" )
    end if
    if ( never_doAlb .and. (nextsw_cday > 0) .and. (nstep > 2) )then
       write(iulog,*)'nstep=',nstep, 'nextsw_cday=', nextsw_cday, 'caldayp1=', caldayp1
       call endrun( sub//":: doalb never set to true over coupling interval -- something is wrong" )
    end if

    
    ! Reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    first_call  = .false.

  end subroutine lnd_run_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_final_mct
!
! !INTERFACE:
  subroutine lnd_final_mct( )
!
! !DESCRIPTION:
! Finalize land surface model
!
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

   ! fill this in
  end subroutine lnd_final_mct

!=================================================================================

  subroutine lnd_SetgsMap_mct( mpicom_lnd, LNDID, gsMap_lnd )

    !-------------------------------------------------------------------
    !
    ! Uses
    !
    use shr_kind_mod , only : r8 => shr_kind_r8
    use decompMod    , only : get_proc_bounds_atm, adecomp
    use domainMod    , only : adomain
    use mct_mod      , only : mct_gsMap, mct_gsMap_init
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_lnd
    integer        , intent(in)  :: LNDID
    type(mct_gsMap), intent(out) :: gsMap_lnd
    !
    ! Local Variables
    !
    integer,allocatable :: gindex(:)
    integer :: i, j, n, gi
    integer :: lsize,gsize
    integer :: ier
    integer :: begg, endg    
    !-------------------------------------------------------------------

    ! Build the land grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    call get_proc_bounds_atm(begg, endg)

    allocate(gindex(begg:endg),stat=ier)

    ! number the local grid

    do n = begg, endg
        gindex(n) = adecomp%gdc2glo(n)
    end do
    lsize = endg-begg+1
    gsize = adomain%ni*adomain%nj

    call mct_gsMap_init( gsMap_lnd, gindex, mpicom_lnd, LNDID, lsize, gsize )

    deallocate(gindex)

  end subroutine lnd_SetgsMap_mct

!=================================================================================

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd_chkAerDep_mct
!
! !INTERFACE:
  subroutine lnd_chkAerDep_mct( x2l_l )
!
! !DESCRIPTION:
! Check aerosol deposition values sent from atmosphere to make sure data is filled.
! If any data is set to special-value or if all data is equal to zero than mark
! data as NOT filled. Model will abort if the data sent from the atmosphere is NOT
! filled.
!
!------------------------------------------------------------------------------
!BOP
! !USES:
    use shr_const_mod    , only : spval => SHR_CONST_SPVAL
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varctl       , only : iulog
    use clm_varctl       , only : set_caerdep_from_file, set_dustdep_from_file  ! This will be removed
    use abortutils       , only : endrun
    use seq_flds_indices , only : index_x2l_Faxa_bcphidry, index_x2l_Faxa_bcphodry, &
                                  index_x2l_Faxa_bcphiwet,                          &
                                  index_x2l_Faxa_ocphidry, index_x2l_Faxa_ocphodry, &
                                  index_x2l_Faxa_ocphiwet,                          &
                                  index_x2l_Faxa_dstdry1, index_x2l_Faxa_dstdry2,   &
                                  index_x2l_Faxa_dstdry3, index_x2l_Faxa_dstdry4,   &
                                  index_x2l_Faxa_dstwet1, index_x2l_Faxa_dstwet2,   &
                                  index_x2l_Faxa_dstwet3, index_x2l_Faxa_dstwet4
    use spmdMod          , only : masterproc
!
! !ARGUMENTS:
    implicit none

    type(mct_aVect)             , intent(inout) :: x2l_l
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!EOP
!---------------------------------------------------------------------------
    !
    ! Local Variables
    !
    logical :: caerdep_filled = .true.     ! Flag if carbon aerosol deposition is filled
    logical :: dustdep_filled = .true.     ! Flag if dust deposition is filled

    ! If ANY values set to special value -- then mark data as NOT filled
    if ( any(abs(x2l_l%rAttr(index_x2l_Faxa_bcphidry,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_bcphodry,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_bcphiwet,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_ocphidry,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_ocphodry,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_ocphiwet,:) - spval)/spval < 0.0001_r8 ) &
    )then
        caerdep_filled = .false.
    end if

    ! If ANY values set to special value -- then mark dust data as NOT filled
    if ( any(abs(x2l_l%rAttr(index_x2l_Faxa_dstdry1,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_dstdry2,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_dstdry3,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_dstdry4,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_dstwet1,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_dstwet2,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_dstwet3,:) - spval)/spval < 0.0001_r8 ) &
    .or. any(abs(x2l_l%rAttr(index_x2l_Faxa_dstwet4,:) - spval)/spval < 0.0001_r8 ) &
    )then
        dustdep_filled = .false.
    end if

    ! If ALL values are set to zero -- then mark carbon aerosol dep. as NOT filled
    if (  all(x2l_l%rAttr(index_x2l_Faxa_bcphidry,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_bcphodry,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_bcphiwet,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_ocphidry,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_ocphodry,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_ocphiwet,:) == 0.0_r8) &
    )then
        caerdep_filled = .false.
    end if

    ! If ALL values are set to zero -- then mark dust dep. as NOT filled
    if (  all(x2l_l%rAttr(index_x2l_Faxa_dstdry1,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_dstdry2,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_dstdry3,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_dstdry4,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_dstwet1,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_dstwet2,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_dstwet3,:) == 0.0_r8) &
    .and. all(x2l_l%rAttr(index_x2l_Faxa_dstwet4,:) == 0.0_r8) &
    )then
        dustdep_filled = .false.
    end if

    if ( caerdep_filled )then
       if ( masterproc ) &
       write(iulog,*) "Using aerosol deposition sent from atmosphere model"
    else
       if ( masterproc )then
          write(iulog,*) "WARNING: Reading carbon aerosol deposition from CLM input file"
          write(iulog,*) "WARNING: aerosol deposition from atm is either spval or all zero"
       end if
       !call endrun( "Aerosol deposition data is sent but NOT filled from the atmosphere model" )
    end if
    if ( dustdep_filled )then
       if ( masterproc ) &
       write(iulog,*) "Using dust deposition sent from atmosphere model"
    else
       if ( masterproc )then
          write(iulog,*) "WARNING: Reading dust deposition from CLM input file"
          write(iulog,*) "WARNING: Dust deposition from atm is either spval or all zero"
       end if
       !call endrun( "Dust deposition data is sent but NOT filled from the atmosphere model" )
    end if
    call shr_sys_flush( iulog )

    ! This will be removed...........
    set_caerdep_from_file = .not. caerdep_filled
    set_dustdep_from_file = .not. dustdep_filled
    ! To here........................

  end subroutine lnd_chkAerDep_mct

!====================================================================================

  subroutine lnd_export_mct( l2a, l2x_l )   

    !-----------------------------------------------------
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_time_manager, only : get_nstep  
    use clm_atmlnd      , only : lnd2atm_type
    use domainMod       , only : adomain
    use decompMod       , only : get_proc_bounds_atm, adecomp
    use seq_flds_indices

    type(lnd2atm_type), intent(inout) :: l2a
    type(mct_aVect)   , intent(inout) :: l2x_l

    integer :: g,i
    integer :: begg, endg    ! beginning and ending gridcell indices
    !-----------------------------------------------------
    
    call get_proc_bounds_atm(begg, endg)

    l2x_l%rAttr(:,:) = 0.0_r8

    ! ccsm sign convention is that fluxes are positive downward

!dir$ concurrent
    do g = begg,endg
       i = 1 + (g-begg)
       l2x_l%rAttr(index_l2x_Sl_landfrac,i) =  adomain%frac(g)
       l2x_l%rAttr(index_l2x_Sl_t,i)        =  l2a%t_rad(g)
       l2x_l%rAttr(index_l2x_Sl_snowh,i)    =  l2a%h2osno(g)
       l2x_l%rAttr(index_l2x_Sl_avsdr,i)    =  l2a%albd(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidr,i)    =  l2a%albd(g,2)
       l2x_l%rAttr(index_l2x_Sl_avsdf,i)    =  l2a%albi(g,1)
       l2x_l%rAttr(index_l2x_Sl_anidf,i)    =  l2a%albi(g,2)
       l2x_l%rAttr(index_l2x_Sl_tref,i)     =  l2a%t_ref2m(g)
       l2x_l%rAttr(index_l2x_Sl_qref,i)     =  l2a%q_ref2m(g)
       l2x_l%rAttr(index_l2x_Fall_taux,i)   = -l2a%taux(g)
       l2x_l%rAttr(index_l2x_Fall_tauy,i)   = -l2a%tauy(g)
       l2x_l%rAttr(index_l2x_Fall_lat,i)    = -l2a%eflx_lh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_sen,i)    = -l2a%eflx_sh_tot(g)
       l2x_l%rAttr(index_l2x_Fall_lwup,i)   = -l2a%eflx_lwrad_out(g)
       l2x_l%rAttr(index_l2x_Fall_evap,i)   = -l2a%qflx_evap_tot(g)
       l2x_l%rAttr(index_l2x_Fall_swnet,i)  =  l2a%fsa(g)
       if (index_l2x_Fall_nee /= 0) then
          l2x_l%rAttr(index_l2x_Fall_nee,i) = -l2a%nee(g)  
       end if

       ! optional fields for dust.  The index = 0 is a good way to flag it,
       ! but I have set it up so that l2a doesn't have ram1,fv,flxdst[1-4] if
       ! progsslt or dust aren't running. 
#if ( defined DUST || defined PROGSSLT )
       if (index_l2x_Sl_ram1 /= 0 )  l2x_l%rAttr(index_l2x_Sl_ram1,i) = l2a%ram1(g)
       if (index_l2x_Sl_fv   /= 0 )  l2x_l%rAttr(index_l2x_Sl_fv,i)   = l2a%fv(g)
#endif
#if ( defined DUST )
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst1,i)= -l2a%flxdst(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst2,i)= -l2a%flxdst(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst3,i)= -l2a%flxdst(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x_l%rAttr(index_l2x_Fall_flxdst4,i)= -l2a%flxdst(g,4)
#endif
    end do

  end subroutine lnd_export_mct

!====================================================================================

  subroutine lnd_import_mct( x2l_l, a2l )

! 27 February 2008: Keith Oleson; Forcing height change
    !-----------------------------------------------------
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use clm_atmlnd      , only: atm2lnd_type
    use clm_varctl      , only: co2_type, co2_ppmv
    use clm_varcon      , only: rair, o2_molar_const, c13ratio
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    use decompMod       , only: get_proc_bounds_atm
    use abortutils      , only: endrun
    use clm_varctl      , only : set_caerdep_from_file, set_dustdep_from_file ! will be removed
    use clm_varctl      , only: iulog
    use mct_mod         , only: mct_aVect
    use seq_flds_indices
    !
    ! Arguments
    !
    type(mct_aVect)   , intent(inout) :: x2l_l
    type(atm2lnd_type), intent(inout) :: a2l
    !
    ! Local Variables
    !
    integer  :: g,i,nstep,ier
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e                    !vapor pressure (Pa)
    real(r8) :: qsat                 !saturation specific humidity (kg/kg)
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag        ! temporary
    real(r8) :: co2_ppmv_prog        ! temporary
    real(r8) :: co2_ppmv_val         ! temporary
    integer  :: begg, endg           ! beginning and ending gridcell indices
    integer  :: co2_type_idx         ! integer flag for co2_type options
    real(r8) :: esatw                !saturation vapor pressure over water (Pa)
    real(r8) :: esati                !saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 !coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 !coefficients for esat over ice
    real(r8) :: tdc, t               !Kelvins to Celcius function and its input
    character(len=32), parameter :: sub = 'lnd_import_mct'

    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
               a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
               a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
               a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
               b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
               b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
               b6=1.838826904e-10_r8)
!
! function declarations
!
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))

    !-----------------------------------------------------

    call get_proc_bounds_atm(begg, endg)

    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if
    
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.
    

!dir$ concurrent
    do g = begg,endg
        i = 1 + (g - begg)
       
        ! Determine required receive fields

        a2l%forc_hgt(g)     = x2l_l%rAttr(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
        a2l%forc_u(g)       = x2l_l%rAttr(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
        a2l%forc_v(g)       = x2l_l%rAttr(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
        a2l%forc_th(g)      = x2l_l%rAttr(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
        a2l%forc_q(g)       = x2l_l%rAttr(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
        a2l%forc_pbot(g)    = x2l_l%rAttr(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
        a2l%forc_t(g)       = x2l_l%rAttr(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
        a2l%forc_lwrad(g)   = x2l_l%rAttr(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
        forc_rainc          = x2l_l%rAttr(index_x2l_Faxa_rainc,i)   ! mm/s
        forc_rainl          = x2l_l%rAttr(index_x2l_Faxa_rainl,i)   ! mm/s
        forc_snowc          = x2l_l%rAttr(index_x2l_Faxa_snowc,i)   ! mm/s
        forc_snowl          = x2l_l%rAttr(index_x2l_Faxa_snowl,i)   ! mm/s
        a2l%forc_solad(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
        a2l%forc_solad(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
        a2l%forc_solai(g,2) = x2l_l%rAttr(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
        a2l%forc_solai(g,1) = x2l_l%rAttr(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

        ! atmosphere coupling, if using prognostic aerosols
        ! This if will be removed so always on....
        if ( .not. set_caerdep_from_file ) then
           a2l%forc_aer(g,1) =  x2l_l%rAttr(index_x2l_Faxa_bcphidry,i)
           a2l%forc_aer(g,2) =  x2l_l%rAttr(index_x2l_Faxa_bcphodry,i)
           a2l%forc_aer(g,3) =  x2l_l%rAttr(index_x2l_Faxa_bcphiwet,i)
           a2l%forc_aer(g,4) =  x2l_l%rAttr(index_x2l_Faxa_ocphidry,i)
           a2l%forc_aer(g,5) =  x2l_l%rAttr(index_x2l_Faxa_ocphodry,i)
           a2l%forc_aer(g,6) =  x2l_l%rAttr(index_x2l_Faxa_ocphiwet,i)
        endif
        ! This if will be removed so always on....
        if ( .not. set_dustdep_from_file ) then
           a2l%forc_aer(g,7)  =  x2l_l%rAttr(index_x2l_Faxa_dstwet1,i)
           a2l%forc_aer(g,8)  =  x2l_l%rAttr(index_x2l_Faxa_dstdry1,i)
           a2l%forc_aer(g,9)  =  x2l_l%rAttr(index_x2l_Faxa_dstwet2,i)
           a2l%forc_aer(g,10) =  x2l_l%rAttr(index_x2l_Faxa_dstdry2,i)
           a2l%forc_aer(g,11) =  x2l_l%rAttr(index_x2l_Faxa_dstwet3,i)
           a2l%forc_aer(g,12) =  x2l_l%rAttr(index_x2l_Faxa_dstdry3,i)
           a2l%forc_aer(g,13) =  x2l_l%rAttr(index_x2l_Faxa_dstwet4,i)
           a2l%forc_aer(g,14) =  x2l_l%rAttr(index_x2l_Faxa_dstdry4,i)
        endif

        ! Determine optional receive fields

        if (index_x2l_Sa_co2prog /= 0) then
           co2_ppmv_prog = x2l_l%rAttr(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
        else
           co2_ppmv_prog = co2_ppmv
        end if
 
        if (index_x2l_Sa_co2diag /= 0) then
           co2_ppmv_diag = x2l_l%rAttr(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
        else
           co2_ppmv_diag = co2_ppmv
        end if

        ! Determine derived quantities for required fields
        a2l%forc_hgt_u(g) = a2l%forc_hgt(g)    !observational height of wind [m]
        a2l%forc_hgt_t(g) = a2l%forc_hgt(g)    !observational height of temperature [m]
        a2l%forc_hgt_q(g) = a2l%forc_hgt(g)    !observational height of humidity [m]
        a2l%forc_vp(g)    = a2l%forc_q(g) * a2l%forc_pbot(g) &
                            / (0.622_r8 + 0.378_r8 * a2l%forc_q(g))
        a2l%forc_rho(g)   = (a2l%forc_pbot(g) - 0.378_r8 * a2l%forc_vp(g)) &
                            / (rair * a2l%forc_t(g))
        a2l%forc_po2(g)   = o2_molar_const * a2l%forc_pbot(g)
        a2l%forc_wind(g)  = sqrt(a2l%forc_u(g)**2 + a2l%forc_v(g)**2)
        a2l%forc_solar(g) = a2l%forc_solad(g,1) + a2l%forc_solai(g,1) + &
                            a2l%forc_solad(g,2) + a2l%forc_solai(g,2)
        a2l%forc_rain(g)  = forc_rainc + forc_rainl
        a2l%forc_snow(g)  = forc_snowc + forc_snowl
        a2l%rainf    (g)  = a2l%forc_rain(g) + a2l%forc_snow(g)

        if (a2l%forc_t(g) > SHR_CONST_TKFRZ) then
           e = esatw(tdc(a2l%forc_t(g)))
        else
           e = esati(tdc(a2l%forc_t(g)))
        end if
        qsat           = 0.622_r8*e / (a2l%forc_pbot(g) - 0.378_r8*e)
        a2l%forc_rh(g) = 100.0_r8*(a2l%forc_q(g) / qsat)
        ! Make sure relative humidity is properly bounded
        ! a2l%forc_rh(g) = min( 100.0_r8, a2l%forc_rh(g) )
        ! a2l%forc_rh(g) = max(   0.0_r8, a2l%forc_rh(g) )
        
        ! Determine derived quantities for optional fields
        ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
        ! Note that forc_pbot is in Pa

        if (co2_type_idx == 1) then
           co2_ppmv_val = co2_ppmv_prog
        else if (co2_type_idx == 2) then
           co2_ppmv_val = co2_ppmv_diag 
        else
           co2_ppmv_val = co2_ppmv
        end if
        a2l%forc_pco2(g)   = co2_ppmv_val * 1.e-6_r8 * a2l%forc_pbot(g) 
        a2l%forc_pc13o2(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * a2l%forc_pbot(g)
	 
     end do

   end subroutine lnd_import_mct

!===============================================================================

  subroutine lnd_domain_mct( lsize, gsMap_l, dom_l )

    !-------------------------------------------------------------------
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varcon  , only : re
    use domainMod   , only : adomain
    use decompMod   , only : get_proc_bounds_atm, adecomp
    use spmdMod     , only : iam
    use mct_mod     , only : mct_gsMap, mct_gGrid, mct_gGrid_importIAttr, &
                             mct_gGrid_importRAttr, mct_gGrid_init,       &
                             mct_gsMap_orderedPoints
    use seq_flds_mod
    !
    ! Arguments
    !
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(inout) :: gsMap_l
    type(mct_ggrid), intent(out)   :: dom_l      
    !
    ! Local Variables
    !
    integer :: g,i,j              ! index
    integer :: begg, endg         ! beginning and ending gridcell indices
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    call mct_gGrid_init( GGrid=dom_l, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_l, iam, idata)
    call mct_gGrid_importIAttr(dom_l,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_l,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_l,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_l,"mask" ,data,lsize) 
    !
    ! Determine bounds
    !
    call get_proc_bounds_atm(begg, endg)
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = adomain%lonc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lon",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = adomain%latc(g)
    end do
    call mct_gGrid_importRattr(dom_l,"lat",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = adomain%area(g)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_l,"area",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = real(adomain%mask(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"mask",data,lsize) 

    do g = begg,endg
       i = 1 + (g - begg)
       data(i) = real(adomain%frac(g), r8)
    end do
    call mct_gGrid_importRattr(dom_l,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine lnd_domain_mct
    
!===============================================================================
    
#ifdef RTM
  subroutine rof_SetgsMap_mct( mpicom_l, LNDID, gsMap_r )

    !-------------------------------------------------------------------
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varpar  , only : rtmlon, rtmlat
    use RunoffMod   , only : runoff
    use abortutils  , only : endrun
    use clm_varctl  , only : iulog
    use mct_mod     , only : mct_gsMap, mct_gsMap_init
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_l
    integer        , intent(in)  :: LNDID
    type(mct_gsMap), intent(out) :: gsMap_r
    !
    ! Local Variables
    !
    integer,allocatable :: gindex(:)
    integer :: n, ni
    integer :: lsize,gsize
    integer :: ier
    character(len=32), parameter :: sub = 'rof_SetgsMap_mct'
    !-------------------------------------------------------------------

    ! Build the rof grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    gsize = rtmlon*rtmlat
    lsize = runoff%lnumro
    allocate(gindex(lsize),stat=ier)

    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          if (ni > runoff%lnumro) then
             write(iulog,*) sub, ' : ERROR runoff count',n,ni,runoff%lnumro
             call endrun( sub//' ERROR: runoff > expected' )
          endif
          gindex(ni) = runoff%gindex(n)
       endif
    end do
    if (ni /= runoff%lnumro) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumro
       call endrun( sub//' ERROR: runoff not equal to expected' )
    endif

    call mct_gsMap_init( gsMap_r, gindex, mpicom_l, LNDID, lsize, gsize )

    deallocate(gindex)

  end subroutine rof_SetgsMap_mct

!===============================================================================

  subroutine rof_domain_mct( lsize, gsMap_r, dom_r )

    !-------------------------------------------------------------------
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varcon  , only : re
    use RunoffMod   , only : runoff
    use abortutils  , only : endrun
    use clm_varctl  , only : iulog
    use spmdMod     , only : iam
    use mct_mod     , only : mct_gsMap, mct_gGrid, mct_gGrid_importIAttr, &
                             mct_gGrid_importRAttr, mct_gGrid_init, mct_gsMap_orderedPoints
    use seq_flds_mod
    use seq_flds_indices
    !
    ! Arguments
    !
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(inout) :: gsMap_r
    type(mct_ggrid), intent(out)   :: dom_r      
    !
    ! Local Variables
    !
    integer :: n, ni              ! index
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    character(len=32), parameter :: sub = 'rof_domain_mct'
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    call mct_gGrid_init( GGrid=dom_r, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_r, iam, idata)
    call mct_gGrid_importIAttr(dom_r,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_r,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_r,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_r,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_r,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_r,"mask" ,data,lsize) 
    !
    ! Determine bounds numbering consistency
    !
    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          if (ni > runoff%lnumro) then
             write(iulog,*) sub, ' : ERROR runoff count',n,ni,runoff%lnumro
             call endrun( sub//' ERROR: runoff > expected' )
          endif
       end if
    end do
    if (ni /= runoff%lnumro) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumro
       call endrun( sub//' ERROR: runoff not equal to expected' )
    endif
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          data(ni) = runoff%lonc(n)
       end if
    end do
    call mct_gGrid_importRattr(dom_r,"lon",data,lsize) 

    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          data(ni) = runoff%latc(n)
       end if
    end do
    call mct_gGrid_importRattr(dom_r,"lat",data,lsize) 

    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          data(ni) = runoff%area(n)*1.0e-6_r8/(re*re)
       end if
    end do
    call mct_gGrid_importRattr(dom_r,"area",data,lsize) 

    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          data(ni) = 1.0_r8
       end if
    end do
    call mct_gGrid_importRattr(dom_r,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_r,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine rof_domain_mct

!====================================================================================

  subroutine rof_export_mct( r2x_r)

    use shr_kind_mod, only : r8 => shr_kind_r8
    use RunoffMod   , only : runoff, nt_rtm, rtm_tracers
    use abortutils  , only : endrun
    use clm_varctl  , only : iulog
    use mct_mod     , only : mct_aVect
    use seq_flds_indices
    !-----------------------------------------------------
    !
    ! Arguments
    ! 
    type(mct_aVect), intent(inout) :: r2x_r
    !
    ! Local variables
    !
    integer :: ni, n, nt, nliq, nfrz
    character(len=32), parameter :: sub = 'rof_export_mct'
    !-----------------------------------------------------
    
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*)'RtmUpdateInput: ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call endrun()
    endif

    ni = 0
    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 2) then
          ni = ni + 1
          r2x_r%rAttr(index_r2x_Forr_roff,ni) = runoff%runoff(n,nliq)/(runoff%area(n)*1.0e-6_r8*1000._r8)
          r2x_r%rAttr(index_r2x_Forr_ioff,ni) = runoff%runoff(n,nfrz)/(runoff%area(n)*1.0e-6_r8*1000._r8)
          if (ni > runoff%lnumro) then
             write(iulog,*) sub, ' : ERROR runoff count',n,ni
             call endrun( sub//' : ERROR runoff > expected' )
          endif
       endif
    end do
    if (ni /= runoff%lnumro) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumro
       call endrun( sub//' : ERROR runoff not equal to expected' )
    endif

  end subroutine rof_export_mct
#endif

!====================================================================================

#endif

end module lnd_comp_mct
