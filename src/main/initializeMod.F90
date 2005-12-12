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
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initialize
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

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize
!
! !INTERFACE:
  subroutine initialize( )
!
! !DESCRIPTION:
! Land model initialization.
! o Initializes run control variables via the [clmexp] namelist.
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
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_sys_mod     , only : shr_sys_flush
    use spmdMod         , only : masterproc
    use clm_varpar      , only : lsmlon, lsmlat, maxpatch
    use clm_varctl      , only : finidat, fpftdyn, fndepdyn, nsrest 
    use clmtypeInitMod  , only : initClmtype
    use initGridCellsMod, only : initGridCells
    use controlMod      , only : control_init, control_print
    use filterMod       , only : allocFilters, setFilters
    use pftdynMod       , only : pftdyn_init, pftdyn_interp
    use decompMod       , only : initDecomp, get_proc_clumps, get_clump_bounds, get_proc_bounds
    use histFldsMod     , only : initHistFlds
    use histFileMod     , only : htapes_build
    use restFileMod     , only : restFile_getfile, restFile_open, restFile_close, &
                                 restFile_read, restFile_read_binary
    use accFldsMod      , only : initAccFlds, initAccClmtype
    use surfFileMod     , only : surfrd
    use mkarbinitMod    , only : mkarbinit
    use ndepFileMod     , only : ndepdyn_init, ndepdyn_interp
#if (defined DGVM)
    use DGVMMod            , only : resetTimeConstDGVM, resetWeightsDGVM, gatherWeightsDGVM
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
#if (defined COUP_CAM)
    use time_manager    , only : get_curr_date, get_nstep 
#else
    use time_manager    , only : get_curr_date, get_nstep, advance_timestep, &
                                 timemgr_init, timemgr_restart
#endif
    use abortutils      , only : endrun
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine program_off if cpp token OFFLINE is defined
! routine program_csm if cpp token COUP_CSM is defined
! routine atmlnd_ini in module atm_lndMod if cpp token COUP_CAM is defined
!
! !REVISION HISTORY:
! Created by Gordon Bonan, Sam Levis and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,k                 ! indices
    integer  :: yr                    ! current year (0, ...)
    integer  :: mon                   ! current month (1 -> 12)
    integer  :: day                   ! current day (1 -> 31)
    integer  :: ncsec                 ! current time of day [seconds]
    integer  :: nc                    ! clump index
    integer  :: nclumps               ! number of clumps on this processor
    integer  :: begp, endp            ! clump beginning and ending pft indices
    integer  :: begc, endc            ! clump beginning and ending column indices
    integer  :: begl, endl            ! clump beginning and ending landunit indices
    integer  :: begg, endg            ! clump beginning and ending gridcell indices
    integer  :: ier                   ! error status
    integer , allocatable :: vegxy(:,:,:) ! vegetation type
    real(r8), allocatable :: wtxy(:,:,:)  ! subgrid weights
    character(len=256) :: fnamer          ! name of netcdf restart file 
    character(len=256) :: pnamer          ! full pathname of netcdf restart file
    character(len=256) :: fnamer_bin      ! name of binary restart file
    integer  :: ncid                      ! netcdf id
!-----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Initialize run control variables, timestep
    ! ------------------------------------------------------------------------

    call header()

    if (masterproc) then
       write (6,*) 'Attempting to initialize the land model .....'
       write (6,*)
    endif
    call control_init ()
    if (masterproc) call control_print()

    ! ------------------------------------------------------------------------
    ! Initialize the subgrid hierarchy
    ! ------------------------------------------------------------------------

    ! Allocate surface grid dynamic memory (for wtxy and vegxy arrays)

    allocate (vegxy(lsmlon,lsmlat,maxpatch), wtxy(lsmlon,lsmlat,maxpatch), stat=ier)   
    if (ier /= 0) then
       write(6,*)'initialize allocation error'; call endrun()
    endif

    ! Read surface dataset and set up vegetation type [vegxy] and weight [wtxy] arrays
    ! for [maxpatch] subgrid patches.
    
    call surfrd (vegxy, wtxy)
    
    ! Initialize clump and processor decomposition 
    
    call initDecomp(wtxy)   
    
    ! Allocate memory and initialize values of clmtype data structures

    call initClmtype()

    ! Build hierarchy and topological info for derived typees

    call initGridCells(vegxy, wtxy)

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
    
    ! Initialize river routing model 

#if (defined RTM)
    if (masterproc) write(6,*)'Attempting to initialize RTM'
    call Rtmini()
    if (masterproc) write(6,*)'Successfully initialized RTM'
#endif

    ! ------------------------------------------------------------------------
    ! Obtain restart file if appropriate
    ! ------------------------------------------------------------------------

    if (do_restread()) then
       call restFile_getfile( file=fnamer, path=pnamer )
    end if

    ! ------------------------------------------------------------------------
    ! Initialize time manager
    ! ------------------------------------------------------------------------

#if (! defined COUP_CAM)
    if (nsrest == 0) then  
       call timemgr_init()
    else
       call restFile_open( flag='read', file=fnamer, ncid=ncid )
       call timemgr_restart( ncid=ncid, flag='read' )
       call restFile_close( ncid=ncid )
    end if
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
       
#if (!defined SCAM)
    ! For restart run, read binary history restart

    if (nsrest == 1) then
       if (masterproc) fnamer_bin = fnamer(1:(len_trim(fnamer)-3))
       call restFile_read_binary( fnamer_bin )
    end if
#endif

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

#if (defined CASA) && (!defined SCAM) 
    ! Initialize CASA 

    call initCASA()
    call initCASAPhenology()
#endif

    ! ------------------------------------------------------------------------
    ! Initialize history and accumator buffers
    ! ------------------------------------------------------------------------

#if (!defined SCAM) 
    ! Initialize master history list. 

    call initHistFlds()

    ! Initialize active history fields. This is only done if not a restart run. 
    ! If a restart run, then this information has already been obtained from the 
    ! restart data read above. Note that routine htapes_build needs time manager 
    ! information, so this call must be made after the restart information has been read.

    if (nsrest == 0 .or. nsrest == 3) call htapes_build ()
#endif

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
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
    do nc = 1,nclumps
       call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
       call resetWeightsDGVM(begg, endg, begc, endc, begp, endp)
       call resetTimeConstDGVM(begp, endp)
    end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO
    call gatherWeightsDGVM()
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

  end subroutine initialize

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
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_varctl  , only : version
    use spmdMod     , only : masterproc
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

    version = 'CLM MODEL version 3.1'
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
    use clm_varctl, only : finidat, nsrest
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
#if (!defined SCAM)
    if (nsrest == 1 .or. nsrest == 3) then
       do_restread = .true.
    end if
#endif
  end function do_restread
  
end module initializeMod
