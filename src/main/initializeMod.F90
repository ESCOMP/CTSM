#include <misc.h>
#include <preproc.h>

module initializeMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initializeMod
!
! !DESCRIPTION:
! Performs land model initialization
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use spmdMod         , only : masterproc
  use shr_sys_mod     , only : shr_sys_flush
  use abortutils      , only : endrun
  use clm_varctl      , only : nsrest
  use clm_varsur      , only : wtxy,vegxy
! !PUBLIC TYPES:
  implicit none
  save
!
  private    ! By default everything is private

! !PUBLIC MEMBER FUNCTIONS:
  public :: initialize1
  public :: initialize2
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private header    ! echo version numbers
  private do_restread
!-----------------------------------------------------------------------
! !PRIVATE DATA:

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize1
!
! !INTERFACE:
  subroutine initialize1( CCSMInit )
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clm_inparm] namelist.
! o Reads surface data on model grid.
! o Defines the multiple plant types and fraction areas for each surface type.
! o Builds the appropriate subgrid <-> grid mapping indices and weights.
! o Set up parallel processing.
! o Initializes time constant variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an initial run.
! o Initializes history file output.
! o Initializes river routing model.
! o Initializes accumulation variables.
!
! !USES:
    use clm_varpar, only : lsmlon, lsmlat, maxpatch
    use clm_varpar, only : clm_varpar_init
    use decompMod , only : decomp_atm_init, decomp_domg2l
    use decompMod , only : adecomp,ldecomp
    use decompMod , only : get_proc_clumps, get_clump_bounds, &
                           get_proc_bounds, get_proc_bounds_atm, &
                           decomp_lnd_init, decomp_glcp_init
    use domainMod , only : domain_check,domain_setsame, domain_clean, domain_type
    use domainMod , only : adomain,ldomain
    use domainMod , only : alatlon,llatlon,gatm,amask,pftm
    use domainMod , only : latlon_check, latlon_setsame
    use areaMod   , only : map_setgatm
    use surfrdMod , only : surfrd,surfrd_get_grid,surfrd_get_frac,&
                           surfrd_get_topo, surfrd_get_latlon
    use clm_varctl, only : fsurdat, fatmgrid, fatmlndfrc, &
                           fatmtopo, flndtopo
    use controlMod, only : control_init, control_print
    use shr_InputInfo_mod, only : shr_InputInfo_initType
!
! !ARGUMENTS:
    type(shr_InputInfo_initType), intent(in), optional :: CCSMInit
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ier                       ! error status
    integer  :: i,j,n1,n2,n,k             ! loop indices
    integer  :: nl,nlg                ! gdc and glo lnd indices
    real(r8) :: rmaxlon,rmaxlat           ! local min/max vars
    logical  :: samegrids             ! are atm and lnd grids same?
    integer  :: begg, endg            ! clump beg and ending gridcell indices
    integer  :: begg_atm, endg_atm    ! proc beg and ending gridcell indices
    type(domain_type) :: lgdomain   ! generic "global" domain for i/o
    type(domain_type) :: agdomain   ! generic "global" domain for i/o
!-----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Initialize run control variables, timestep
    ! ------------------------------------------------------------------------

    call header()

    if (masterproc) then
       write (6,*) 'Attempting to initialize the land model .....'
       write (6,*)
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif

    call clm_varpar_init ()
    if ( present(CCSMInit) )then
       call control_init ( CCSMInit )
    else
       call control_init ( )
    end if

    if (masterproc) call control_print()

    ! ------------------------------------------------------------------------
    ! Initialize the subgrid hierarchy
    ! ------------------------------------------------------------------------

    if (masterproc) then
       write (6,*) 'Attempting to read alatlon from fatmgrid'
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif

    call surfrd_get_latlon(alatlon, fatmgrid, amask, fatmlndfrc)
    call latlon_check(alatlon)
    if (masterproc) then
       write (6,*) 'amask size/min/max ',size(amask),minval(amask),maxval(amask)
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif

    if (masterproc) then
       write (6,*) 'Attempting to read llatlon from fsurdat'
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif

    call surfrd_get_latlon(llatlon, fsurdat, pftm, pftmflag=.true.)
    call latlon_check(llatlon)

    lsmlon = llatlon%ni
    lsmlat = llatlon%nj

    if (llatlon%ni < alatlon%ni .or. llatlon%nj < alatlon%nj) then
       if (masterproc) write(6,*) 'ERROR llatlon size > alatlon size: ', &
          n,llatlon%ni, llatlon%nj, alatlon%ni, alatlon%nj
       call endrun()
    endif

    !--- check if grids are "close", adjust, continue, or end  ---
    !--- set llatlon== alatlon if lats/lons < 0.001 different ---
    !--- exit if lats/lon > 0.001 and < 1.0 different          ---
    !--- continue if lat/lons > 1.0 different                  ---
    samegrids = .false.
    if (alatlon%ni == llatlon%ni .and. alatlon%nj == llatlon%nj) then
       rmaxlon = 0.0_r8
       rmaxlat = 0.0_r8
       do n = 1,alatlon%ni
          rmaxlon = max(rmaxlon,abs(alatlon%lonc(n)-llatlon%lonc(n)))
       enddo
       do n = 1,alatlon%nj
          rmaxlat = max(rmaxlat,abs(alatlon%latc(n)-llatlon%latc(n)))
       enddo
       if (rmaxlon < 0.001_r8 .and. rmaxlat < 0.001_r8) then
          if (masterproc) write(6,*) 'initialize1: set llatlon =~ alatlon', &
             ':continue',rmaxlon,rmaxlat
          call latlon_setsame(alatlon,llatlon)
          samegrids = .true.
       elseif (rmaxlon < 1.0_r8 .and. rmaxlat < 1.0_r8) then
          if (masterproc) write(6,*) 'initialize1: alatlon/llatlon mismatch', &
             ':error',rmaxlon,rmaxlat
          call endrun()
       else
          if (masterproc) write(6,*) 'initialize1: alatlon/llatlon different', &
              ':continue',rmaxlon,rmaxlat
       endif
    else
       if (masterproc) write(6,*) 'initialize1: alatlon/llatlon different ', &
          'sizes:continue'
    endif

    call decomp_atm_init(alatlon,amask)

#if (1 == 0)
    !--- set gatm "unity".  doesn't support finemesh grids yet
    if (llatlon%ni /= alatlon%ni .or. llatlon%nj /= alatlon%nj) then
       if (masterproc) write(6,*) 'ERROR llatlon size /= alatlon size: ', &
          n,llatlon%ni, llatlon%nj, alatlon%ni, alatlon%nj
       call endrun()
    else
       call map_setgatm_UNITY(gatm,alatlon,amask)
    endif
#endif
    call map_setgatm(gatm,alatlon,llatlon,amask,pftm)

    ! Initialize clump and processor decomposition 
    call decomp_lnd_init(alatlon%ns,alatlon%ni,alatlon%nj,llatlon%ns,llatlon%ni,llatlon%nj)

    !--- Read atm grid -------------------------------------------------------------------

    if (masterproc) then
       write (6,*) 'Attempting to read agdomain from fatmgrid'
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif
    call surfrd_get_grid(agdomain, fatmgrid)

    if (masterproc) then
       write (6,*) 'Attempting to read atm landfrac from fatmlndfrc'
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif
    call surfrd_get_frac(agdomain, fatmlndfrc)

    if (fatmtopo /= " ") then
    if (masterproc) then
       write (6,*) 'Attempting to read atm topo from fatmtopo'
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif
    call surfrd_get_topo(agdomain, fatmtopo)
    endif

    if (masterproc) then
       write(6,*) 'agdomain status:'
       call domain_check(agdomain)
    endif

    ! Set "local" domains
    call get_proc_bounds_atm(begg_atm, endg_atm)
    call decomp_domg2l(agdomain,adomain,adecomp,begg_atm,endg_atm)
    if (masterproc) then
       write(6,*) 'adomain status:'
       call domain_check(adomain)
    endif

    call domain_clean(agdomain)

    !--- Read lnd grid -------------------------------------------------------------------

    if (masterproc) then
       write (6,*) 'Attempting to read lgdomain from fsurdat ',trim(fsurdat)
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif
    call surfrd_get_grid(lgdomain, fsurdat)

    if (flndtopo /= " ") then
    if (masterproc) then
       write (6,*) 'Attempting to read lnd topo from flndtopo ',trim(flndtopo)
#ifndef UNICOSMP
       call shr_sys_flush(6)
#endif
    endif
    call surfrd_get_topo(lgdomain, flndtopo)
    endif

    if (masterproc) then
       call domain_check(lgdomain)
    endif

    call get_proc_bounds    (begg    , endg)
    call decomp_domg2l(lgdomain,ldomain,ldecomp,begg,endg)
    if (masterproc) then
       write(6,*) 'ldomain status:'
       call domain_check(ldomain)
    endif

    call domain_clean(lgdomain)

    !--- overwrite ldomain if same grids --------------------------------------------

    if (samegrids) then
       if (masterproc) write(6,*) 'initialize1: samegrids true, set ldomain =~ adomain'
       call domain_setsame(adomain,ldomain)
    endif

    ! Allocate surface grid dynamic memory (for wtxy and vegxy arrays)

    call get_proc_bounds    (begg    , endg)
    allocate (vegxy(begg:endg,maxpatch), &
              wtxy(begg:endg,maxpatch),  &
              stat=ier)   
    if (ier /= 0) then
       write(6,*)'initialize allocation error'; call endrun()
    endif

    ! Read surface dataset and set up vegetation type [vegxy] and 
    ! weight [wtxy] arrays for [maxpatch] subgrid patches.
    
    call surfrd (fsurdat, ldomain)

    call decomp_glcp_init(alatlon%ns,alatlon%ni,alatlon%nj,llatlon%ns,llatlon%ni,llatlon%nj)

  end subroutine initialize1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize2
!
! !INTERFACE:
  subroutine initialize2( SyncClock )
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clm_inparm] namelist.
! o Reads surface data on model grid.
! o Defines the multiple plant types and fraction areas for each surface type.
! o Builds the appropriate subgrid <-> grid mapping indices and weights.
! o Set up parallel processing.
! o Initializes time constant variables.
! o Reads restart data for a restart or branch run.
! o Reads initial data and initializes the time variant variables for an initial run.
! o Initializes history file output.
! o Initializes river routing model.
! o Initializes accumulation variables.
!
! !USES:
    use eshr_timemgr_mod, only : eshr_timemgr_clockType, eshr_timemgr_clockInfoType, &
                                 eshr_timemgr_clockGet
    use clm_atmlnd      , only : init_atm2lnd_type, init_lnd2atm_type, &
                                 clm_a2l, clm_l2a, atm_a2l, atm_l2a
    use initGridCellsMod, only : initGridCells
    use clm_varctl      , only : finidat, fpftdyn, fndepdyn
    use clmtypeInitMod  , only : initClmtype
    use domainMod       , only : gatm
    use domainMod       , only : ldomain, adomain
    use decompMod       , only : adecomp,ldecomp
    use areaMod         , only : map1dl_a2l, map1dl_l2a
    use areaMod         , only : map_setmapsFM
    use decompMod       , only : get_proc_clumps, get_clump_bounds, &
                                 get_proc_bounds, get_proc_bounds_atm, &
                                 decomp_lnd_init
    use filterMod       , only : allocFilters, setFilters
    use pftdynMod       , only : pftdyn_init, pftdyn_interp
    use histFldsMod     , only : initHistFlds
    use histFileMod     , only : htapes_build
    use restFileMod     , only : restFile_getfile, &
                                 restFile_open, restFile_close, &
                                 restFile_read, restFile_read_binary
    use accFldsMod      , only : initAccFlds, initAccClmtype
    use mkarbinitMod    , only : mkarbinit
    use ndepFileMod     , only : ndepdyn_init, ndepdyn_interp
#if (defined DGVM)
    use DGVMMod            , only : resetTimeConstDGVM, resetWeightsDGVM
    use DGVMEcosystemDynMod, only : DGVMEcosystemDynini
#else
    use STATICEcosysDynMod , only : EcosystemDynini
#endif
#if (defined DUST) 
    use DustMod         , only : Dustini
#endif
#if (defined CASA)
    use CASAMod         , only : initCASA
    use CASAPhenologyMod, only : initCASAPhenology
#endif
#if (defined RTM) 
    use RtmMod          , only : Rtmini
#endif
    use clm_time_manager, only : get_curr_date, get_nstep, advance_timestep, &
                                 timemgr_init, timemgr_restart_io, timemgr_restart
    use fileutils       , only : getfil
!
! !ARGUMENTS:
   type(eshr_timeMgr_clockType), optional, intent(IN) :: SyncClock ! Synchronization clock
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,k,n1,n2           ! indices
    integer  :: yr                    ! current year (0, ...)
    integer  :: mon                   ! current month (1 -> 12)
    integer  :: day                   ! current day (1 -> 31)
    integer  :: ncsec                 ! current time of day [seconds]
    integer  :: nc                    ! clump index
    integer  :: nclumps               ! number of clumps on this processor
    integer  :: begp, endp            ! clump beg and ending pft indices
    integer  :: begc, endc            ! clump beg and ending column indices
    integer  :: begl, endl            ! clump beg and ending landunit indices
    integer  :: begg, endg            ! clump beg and ending gridcell indices
    integer  :: begg_atm, endg_atm    ! proc beg and ending gridcell indices
    character(len=256) :: fnamer             ! name of netcdf restart file 
    character(len=256) :: pnamer             ! full pathname of netcdf restart file
    character(len=256) :: fnamer_bin         ! name of binary restart file
    character(len=256) :: pnamer_bin         ! full pathname of binary restart file
    integer  :: ncid                         ! netcdf id
    logical ,parameter :: a2ltrue = .true.   ! local
    logical ,parameter :: a2lfalse = .false. ! local
    !
    character(len=80) :: calendar           ! Calendar type
    integer           :: start_ymd          ! Start date (YYYYMMDD)
    integer           :: start_tod          ! Start time of day (sec)
    integer           :: ref_ymd            ! Reference date (YYYYMMDD)
    integer           :: ref_tod            ! Reference time of day (sec)
    integer           :: stop_ymd           ! Stop date (YYYYMMDD)
    integer           :: stop_tod           ! Stop time of day (sec)
    logical           :: log_print          ! Flag to print out log information or not
    logical           :: perpetual_run      ! If in perpetual mode or not
    integer           :: perpetual_ymd      ! Perpetual date (YYYYMMDD)
!----------------------------------------------------------------------

    ! Set the a2l and l2a maps

    call map_setmapsFM(adomain,ldomain,gatm,map1dl_a2l,map1dl_l2a)

    ! Allocate memory and initialize values of clmtype data structures

    call initClmtype()

    call get_proc_bounds    (begg    , endg)
    call init_atm2lnd_type  (begg    , endg    , clm_a2l)
    call init_lnd2atm_type  (begg    , endg    , clm_l2a)

    call get_proc_bounds_atm(begg_atm, endg_atm)
    call init_atm2lnd_type  (begg_atm, endg_atm, atm_a2l)
    call init_lnd2atm_type  (begg_atm, endg_atm, atm_l2a)

    ! Build hierarchy and topological info for derived typees

    call initGridCells()

    ! Deallocate surface grid dynamic memory (for wtxy and vegxy arrays)

    deallocate (vegxy, wtxy)

    ! ------------------------------------------------------------------------
    ! Initialize time constant variables 
    ! ------------------------------------------------------------------------

    ! Initialize Ecosystem Dynamics 

#if (defined DGVM)
    call DGVMEcosystemDynini()
#elif defined (CN)
    ! currently no call required
#else
    call EcosystemDynini()
#endif

    ! Initialize dust emissions model 

#if (defined DUST) 
    call Dustini()
#endif

    ! Initialize time constant variables (this must be called after
    ! DGVMEcosystemDynini() and before initCASA())

    call iniTimeConst()
    
    ! ------------------------------------------------------------------------
    ! Obtain restart file if appropriate
    ! ------------------------------------------------------------------------

    if (do_restread()) then
       call restFile_getfile( file=fnamer, path=pnamer )
    end if

    ! ------------------------------------------------------------------------
    ! Initialize time manager
    ! ------------------------------------------------------------------------

    if (nsrest == 0) then  
       if (present(SyncClock)) then	
          call eshr_timemgr_clockGet(                                          &
               SyncClock, start_ymd=start_ymd,                                 &
               start_tod=start_tod, ref_ymd=ref_ymd,                           &
               ref_tod=ref_tod, stop_ymd=stop_ymd, stop_tod=stop_tod,          &
               perpetual_run=perpetual_run, perpetual_ymd=perpetual_ymd,       &
               calendar=calendar )
          call timemgr_init(                                                   &
               calendar_in=calendar, start_ymd_in=start_ymd,                   &
               start_tod_in=start_tod, ref_ymd_in=ref_ymd,                     &
               ref_tod_in=ref_tod, stop_ymd_in=stop_ymd, stop_tod_in=stop_tod, &
               perpetual_run_in=perpetual_run, perpetual_ymd_in=perpetual_ymd )
       else
          call timemgr_init()
       end if
    else
       call restFile_open( flag='read', file=fnamer, ncid=ncid )
       call timemgr_restart_io( ncid=ncid, flag='read' )
       call restFile_close( ncid=ncid )
       if (present(SyncClock)) then	
          call eshr_timemgr_clockGet( SyncClock, stop_ymd=stop_ymd, stop_tod=stop_tod )
          call timemgr_restart( stop_ymd, stop_tod )
       else
          call timemgr_restart()
       end if
    end if

    ! Initialize river routing model, after time manager because ts needed

#if (defined RTM)
    if (masterproc) write(6,*)'Attempting to initialize RTM'
    call Rtmini()
    if (masterproc) write(6,*)'Successfully initialized RTM'
#ifndef UNICOSMP
    call shr_sys_flush(6)
#endif
#endif

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! Initialize accumulator fields to be time accumulated for various purposes.
    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed.

    call initAccFlds()

    ! ------------------------------------------------------------------------
    ! Set arbitrary initial conditions for time varying fields 
    ! used in coupled carbon-nitrogen code
    ! ------------------------------------------------------------------------
    
#if (defined CN)
    if (nsrest == 0) then
       call CNiniTimeVar()
    end if
#endif

    ! ------------------------------------------------------------------------
    ! Read restart/initial info
    ! ------------------------------------------------------------------------

    ! No weight related information can be contained in the routines,  
    ! "mkarbinit, inicfile and restFile". 

    if (do_restread()) then
       if (masterproc) write(6,*)'reading restart file ',fnamer
       call restFile_read( fnamer )
    else if (nsrest == 0 .and. finidat == ' ') then
       call mkarbinit()
    end if
       
    ! For restart run, read binary history restart

    if (nsrest == 1) then
       if (masterproc)then 
          pnamer_bin = pnamer(1:(len_trim(pnamer)-3))
          call getfil( pnamer_bin, fnamer_bin )
       end if
       call restFile_read_binary( fnamer_bin )
    end if

    ! ------------------------------------------------------------------------
    ! Initialization of dynamic pft weights
    ! ------------------------------------------------------------------------

    ! Determine correct pft weights (interpolate pftdyn dataset if initial run)
    ! Otherwise these are read in for a restart run
    
    if (fpftdyn /= ' ') then
       call pftdyn_init()
       if (nsrest == 0) call pftdyn_interp()
    end if

    ! ------------------------------------------------------------------------
    ! Initialize dynamic nitrogen deposition
    ! ------------------------------------------------------------------------

    if (fndepdyn /= ' ') then
       call ndepdyn_init()
       call ndepdyn_interp()
    end if

    ! ------------------------------------------------------------------------
    ! Start offline run at nstep = 1
    ! ------------------------------------------------------------------------

#if (defined OFFLINE)
    if (nsrest == 0) call advance_timestep()
#endif

    ! ------------------------------------------------------------------------
    ! Initialization of model parameterizations that are needed after
    ! restart file is read in
    ! ------------------------------------------------------------------------

#if (defined CASA)
    ! Initialize CASA 

    call initCASA()
    call initCASAPhenology()
#endif

    ! ------------------------------------------------------------------------
    ! Initialize history and accumator buffers
    ! ------------------------------------------------------------------------

    ! Initialize master history list. 

    call initHistFlds()

    ! Initialize active history fields. This is only done if not a restart run. 
    ! If a restart run, then this information has already been obtained from the 
    ! restart data read above. Note that routine htapes_build needs time manager 
    ! information, so this call must be made after the restart information has been read.

    if (nsrest == 0 .or. nsrest == 3) call htapes_build ()

    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! This routine is called in an initial run at nstep=0 for cam and csm mode
    ! and at nstep=1 for offline mode. This routine is also always called for a
    ! restart run and must therefore be called after the restart file is read in

    call initAccClmtype()

    ! --------------------------------------------------------------
    ! Note - everything below this point needs updated weights
    ! --------------------------------------------------------------

    ! Initialize filters
    
    call allocFilters()
    call setFilters()

#if (defined DGVM)
    ! Determine new subgrid weights and areas (obtained from fpcgrid) and
    ! reset DGVM time constant variables for input pft types

    nclumps = get_proc_clumps()
!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
#if !defined (USE_OMP)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
#endif
    do nc = 1,nclumps
       call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
       call resetWeightsDGVM(begg, endg, begc, endc, begp, endp)
       call resetTimeConstDGVM(begp, endp)
    end do
#if !defined (USE_OMP)
!CSD$ END PARALLEL DO
#endif
!$OMP END PARALLEL DO
#endif

    ! End initialization

    if (masterproc) then
       write (6,*) 'Successfully initialized the land model'
       if (nsrest == 0) then
          write (6,*) 'begin initial run at: '
       else
          write (6,*) 'begin continuation run at:'
       end if
       call get_curr_date(yr, mon, day, ncsec)
       write (6,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
       write (6,*)
       write (6,'(72a1)') ("*",i=1,60)
       write (6,*)
    endif

  end subroutine initialize2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: header
!
! !INTERFACE:
  subroutine header()
!
! !DESCRIPTION:
! Echo and save model version number
!
! !USES:
    use clm_varctl  , only : version
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!
!EOP
!-----------------------------------------------------------------------

    version = 'CLM MODEL version 3.5'
    if ( masterproc )then
      write (6,*) trim(version)
      write (6,*)
    end if

  end subroutine header

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_restread
!
! !INTERFACE:
  logical function do_restread( )
!
! !DESCRIPTION:
! Determine if restart file will be read
!
! !USES:
    use clm_varctl, only : finidat
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    do_restread = .false.
    if (nsrest == 0 .and. finidat /= ' ') then
       do_restread = .true.
    end if
    if (nsrest == 1 .or. nsrest == 3) then
       do_restread = .true.
    end if
  end function do_restread
  
end module initializeMod
