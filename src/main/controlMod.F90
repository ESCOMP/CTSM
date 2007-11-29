#include <misc.h>
#include <preproc.h>

module controlMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: controlMod
!
! !DESCRIPTION:
! Module which initializes run control variables. The following possible
! namelist variables are set default values and possibly read in on startup
!
! === define run ======================= >>>> Only if mode does NOT = ccsm_seq
!
!    o caseid                = 256 character case name
!    o ctitle                = 256 character case title
!    o nsrest                = integer flag. 0: initial run. 1: restart: 3: branch
!    o version               = 256 character model version description
!    o username              = 256 character user username
!    o hostname              = 256 character hostname of machine running on
!    o brnch_retain_casename = logical if a branch simulation can use the same casename 
!                              as the case it branches from
!
! === model time =======================
!
!    o dtime      = integer model time step (s)
!
! ====================================== >>>> Only if mode does NOT = ccsm_seq
!
!    o calendar       = Calendar to use in date calculations.
!                      'no_leap' (default) or 'gregorian'
!    o start_ymd      = Starting date for run encoded in yearmmdd format.
!                       Default value is read from initial conditions file.
!    o start_tod      = Starting time of day for run in seconds since 0Z.
!                       Default value is read from initial conditions file.
!    o stop_ymd       = Stopping date for run encoded in yearmmdd format.
!                       No default.
!    o stop_tod       = Stopping time of day for run in seconds since 0Z.
!                       Default: 0.
!    o stop_final_ymd = Final stopping time. (>>>> only in offline mode)
!                       Default: 99991231
!    o nelapse        = nnn, Specify the ending time for the run as an interval
!                       starting at the current time in either timesteps
!                       (if positive) or days (if negative).
!                       Either nestep or (stop_ymd,stop_tod) take precedence.
!    o nestep         = nnnn, Specify the ending time for the run as an interval
!                       starting at (start_ymd,start_tod) in either timesteps
!                       (if positive) or days (if negative).
!                       (stop_ymd,stop_tod) takes precedence if set.
!    o ref_ymd        = Reference date for time coordinate encoded in yearmmdd format.
!                       Default value is start_ymd.
!    o ref_tod        = Reference time of day for time coordinate in seconds since 0Z.
!                       Default value is start_tod.
!
! === input data ===
!
!    o finidat         = 256 character initial conditions file name
!    o fsurdat         = 256 character surface data file name
!    o flndtopo        = 256 character land topography file name
!    o fatmgrid        = 256 character atmosphere grid data file name
!    o fatmlndfrc      = 256 character landfrac (on atm grid) file name
!    o fatmtopo        = 256 character atmosphere topography file name (on atm grid)
!    o flndtopo        = 256 character land topography file name
!    o fndepdat        = 254 character nitrogen deposition data file name (netCDF)
!    o fndepdyn        = 254 character nitrogen deposition data file name (netCDF) changing dynamically in time
!    o fpftcon         = 256 character data file with PFT physiological constants
!    o fpftdyn         = 256 character data file with PFT physiological constants changing dynamically in time
!    o frivinp_rtm     = 256 character input data file for rtm
!    o nrevsn          = 256 character restart file name for use with branch run
!
! === offline forcing data === >>> Only for mode=offline
!
!    o offline_atmdir  = 256 character directory for input atm data files (can be Mass Store)
!    o cycle_begyr     = integer, first year of offline atm data (e.g. 1948)
!    o cycle_nyr       = integer, number of years of offline atm data to cycle
!
! === history and restart files ===
!
!    o hist_empty_htapes = true if do not want ANY output fields on files by default
!    o hist_avgflag_pertape = averaging flag to use by default for all fields on a file series
!    o hist_ndens    = integer, can have value of 1 (nc_double) or 2 (nf_float)
!    o hist_dov2xy   = true if want grid-average history field (false = vector)
!    o hist_nhtfrq   = integer history interval (+ = iterations,  - = hours, 0=monthly ave)
!    o hist_mfilt    = integer number of time samples per history file
!    o hist_fincl1   = 10 character name of fields for first  auxillary history file
!    o hist_fincl2   = 10 character name of fields for second auxillary history file
!    o hist_fincl3   = 10 character name of fields for first  auxillary history file
!    o hist_fincl4   = 10 character name of fields for second auxillary history file
!    o hist_fincl5   = 10 character name of fields for first  auxillary history file
!    o hist_fincl6   = 10 character name of fields for second auxillary history file
!    o hist_fexcl1   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl2   = 8  character name of fields for second auxillary history file
!    o hist_fexcl3   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl4   = 8  character name of fields for second auxillary history file
!    o hist_fexcl5   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl6   = 8  character name of fields for second auxillary history file
!    o hist_crtinic  = 8  character frequency to generate initial dataset
!                         ['6-HOURLY','DAILY','MONTHLY','YEARLY','NONE']
!    o rest_flag     = logical, turns off restart file writing [.TRUE., .FALSE.]
! === Biogeochem=CASA =======
!    o lnpp        = 1=gpp*gppfact,2=fn(lgrow)*gppfact
!    o lalloc      = 0=fixed allocation, 1=dynamic allocation
!    o q10         = temperature dependence
!    o spunup      = 0=no, 1=yes (used with nsrest/=1 only)
!    o fcpool      = Carbon Pool initial state filename
!
! === Decomposition   =======
!
!    o clump_pproc = clumps per processor
!    o nsegspc     = number of segments per clump for decomposition
!
! === model physics ===
!
!    o irad         = integer solar radiation frequency (+ = iteration. - = hour)
!    o wrtdia       = true if want output written
!    o csm_doflxave = true => flux averaging is to be performed (only used for csm mode)
!    o co2_ppmv     = CO2 volume mixing ratio
!    o pertlim      = perturbation limit
!    o create_croplandunit = logical on if to create crop as separate landunits
!
!    >>>>>>>>>>>>> Only for mode=ext_ccsm_con or seq_ccsm
!    o co2_type     = type of CO2 feedback, choices are constant, prognostic or diagnostic
!
! === rtm control variables ===
!
!    o rtm_nsteps  = if > 1, average rtm over rtm_nsteps time steps
!
! When coupled to CAM: base calendar info, nstep, nestep, nsrest, and time
! step are input to the land model from CAM. The values in the clm_inparm namelist
! are not used. 
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, SHR_KIND_CL
  use clm_varpar   , only : maxpatch_pft, numpft
  use clm_varctl
  use spmdMod
  use decompMod    , only : clump_pproc
  use histFileMod  , only : max_tapes, max_namlen, &
                            hist_empty_htapes, hist_dov2xy, &
                            hist_avgflag_pertape, hist_type1d_pertape, &
                            hist_nhtfrq, hist_ndens, hist_mfilt, &
                            hist_fincl1, hist_fincl2, hist_fincl3, &
                            hist_fincl4, hist_fincl5, hist_fincl6, &
                            hist_fexcl1, hist_fexcl2, hist_fexcl3, &
                            hist_fexcl4, hist_fexcl5, hist_fexcl6
  use restFileMod  , only : rest_flag
  use shr_const_mod, only : SHR_CONST_CDAY
  use abortutils   , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: control_setNL ! Set namelist filename
  public :: control_init  ! initial run control information
  public :: control_print ! print run control information
!
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! PRIVATE TYPES:
! Namelist variables only used locally
  character(len=256) :: rpntpath                         ! full UNIX pathname of restart pointer file
  character(len=  7) :: runtyp(4)                        ! run type
  character(len=SHR_KIND_CL) :: NLFilename = 'lnd.stdin' ! Namelist filename
#if (defined _OPENMP)
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
                                             ! concurrently in a single parallel region
#endif
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_setNL
!
! !INTERFACE:
  subroutine control_setNL( NLfile )

    implicit none
!
! !DESCRIPTION:
! Set the namelist filename to use
!
!
! !ARGUMENTS:
  character(len=*), intent(IN) :: NLFile ! Namelist filename
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=32) :: subname = 'control_setNL'  ! subroutine name
    logical :: lexist                               ! File exists

    ! Error checking...
    if ( len_trim(NLFile) == 0 )then
       call endrun( subname//' error: nlfilename entered is not set' )
    end if
    inquire (file = trim(NLFile), exist = lexist)
    if ( .not. lexist )then
       call endrun( subname//' error: NLfilename entered does NOT exist:'//trim(NLFile) )
    end if
    if ( len_trim(NLFile) > len(NLFilename) )then
       call endrun( subname//' error: entered NLFile is too long' )
    end if
    ! Set the filename
    NLFilename = NLFile
  end subroutine control_setNL

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_init
!
! !INTERFACE:
  subroutine control_init( CCSMInit )
!
! !DESCRIPTION:
! Initialize CLM run control information
!
! !USES:
    use clm_time_manager
#if (defined CASA)
    use CASAMod          , only : lnpp, lalloc, q10, spunup, fcpool
#endif
    use shr_inputinfo_mod, only : shr_inputInfo_initType,       &
                                  shr_inputInfo_initGetData,    &
                                  shr_inputInfo_initIsBranch,   &
                                  shr_inputInfo_initIsContinue, &
                                  shr_inputInfo_initIsStartup
    use fileutils        , only : getavu, relavu
    use shr_string_mod   , only : shr_string_getParentDir

    implicit none
!
! !ARGUMENTS:
    type(shr_InputInfo_initType), intent(in), optional :: CCSMInit   ! Input CCSM info

    include 'netcdf.inc'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,n                ! loop indices
    integer :: iundef               ! integer undefined value
    real(r8):: rundef               ! real undefined value
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'control_init'  ! subroutine name
!------------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! clm time manager info

#if (defined OFFLINE) || (defined COUP_CSM)
    namelist /clm_inparm/  &
         ctitle, caseid, nsrest, version, hostname, username, &
         calendar, nelapse, nestep, start_ymd, start_tod,     &
         stop_ymd, stop_tod, ref_ymd, ref_tod,                &
         brnch_retain_casename, stop_final_ymd
#endif

    ! clm input datasets

    namelist / clm_inparm/ &
	 dtime	

    namelist /clm_inparm/  &
         finidat, fsurdat, fatmgrid, fatmlndfrc, fatmtopo, flndtopo, &
         fpftcon, frivinp_rtm,  &
         fpftdyn, fndepdat, fndepdyn, nrevsn

#if (defined OFFLINE)
    namelist /clm_inparm/  &
          offline_atmdir, cycle_begyr, cycle_nyr
#endif

    ! clm history, restart options

    namelist /clm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, hist_mfilt, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6, &
         hist_crtinic, rest_flag

    ! clm bgc info

#if (defined CASA)
    namelist /clm_inparm/  &
         lnpp, lalloc, q10, spunup, fcpool
#endif

#if (defined SEQ_MCT) || (defined SEQ_ESMF) || (defined COUP_CSM)
    namelist /clm_inparm / &
         co2_type
#endif

    ! clm other options

    namelist /clm_inparm/  &
         clump_pproc, irad, wrtdia, csm_doflxave, rtm_nsteps, pertlim, &
         create_crop_landunit, nsegspc, co2_ppmv
         
    ! ----------------------------------------------------------------------
    ! Default values
    ! ----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to initialize run control settings .....'
    endif

    runtyp(0 + 1) = 'initial'
    runtyp(1 + 1) = 'restart'
    runtyp(3 + 1) = 'branch '

    iundef = -9999999
    rundef = -9999999._r8

    ! control variables

    caseid      = ' '
    ctitle      = ' '
    nsrest      = iundef
    username    = ' '
    hostname    = ' '
    source      = "Community Land Model CLM3"
    version     = source
    conventions = "CF-1.0"

    ! initial data

    fsurdat     = ' '
    fatmgrid    = ' '
    fatmlndfrc  = ' '
    fatmtopo    = ' '
    flndtopo    = ' '
    fndepdat    = ' '
    fndepdyn    = ' '
    finidat     = ' '
    fpftcon     = ' '
    frivinp_rtm = ' '
    fpftdyn     = ' '
    nrevsn      = ' '

    ! offline mode

    offline_atmdir   = ' '
    cycle_begyr = iundef
    cycle_nyr   = iundef

    ! landunit generation

    create_crop_landunit = .false.
    if (maxpatch_pft == numpft+1) then
       allocate_all_vegpfts = .true.
    else
       allocate_all_vegpfts = .false.
    end if

    ! bgc

    co2_type = 'constant'

    ! history file variables

    hist_crtinic = 'NONE'

    ! other namelist variables

    irad          = -1
    wrtdia        = .false.
    csm_doflxave  = .true.
    pertlim       = 0._r8
    single_column = .false.
    scmlat        = -999.
    scmlon        = -999.
    nsegspc       = 20
    co2_ppmv      = 355._r8

#if (defined RTM)
    ! If rtm_nsteps is not set in the namelist then
    ! will be given default value below

    rtm_nsteps = -999
#endif

#if (defined CASA)
    lnpp = 2
    lalloc = 1
    q10 = 2.0_r8          ! set Q10 to 2.0  03/11/19
    spunup = 0
    fcpool = ' '
#endif

    ! Set clumps per procoessor

    clump_pproc = 1
#if (defined _OPENMP)
    clump_pproc = omp_get_max_threads()
#else
#if (defined UNICOSMP)
#if (defined SSP)
    clump_pproc = 1
#else
    clump_pproc = 1 ! 4 when using CSDs in driver.F90; 1 otherwise
#endif
#endif
#endif

    if (masterproc) then

       ! ----------------------------------------------------------------------
       ! Read namelist from standard input. 
       ! ----------------------------------------------------------------------

       if ( len_trim(NLFilename) == 0  )then
          call endrun( subname//' error: nlfilename not set' )
       end if
       unitn = getavu()
       write(iulog,*) 'Read in clm_inparm namelist from: ', trim(NLFilename)
       open( unitn, file=trim(NLFilename), status='old' )
       ierr = 1
       do while ( ierr /= 0 )
          read(unitn, clm_inparm, iostat=ierr)
          if (ierr < 0) then
             call endrun( subname//' encountered end-of-file on namelist read' )
          endif
       end do
       call relavu( unitn )

       ! ----------------------------------------------------------------------
       ! Consistency checks on input namelist.
       ! ----------------------------------------------------------------------

       ! Consistency settings for co2 type

       if (co2_type /= 'constant' .and. co2_type /= 'prognostic' .and. co2_type /= 'diagnostic') then
          write(iulog,*)'co2_type = ',co2_type,' is not supported'
          write(iulog,*)'choices are constant, prognostic or diagnostic'
          call endrun()
       end if
#if (defined OFFLINE)
       if (co2_type /= 'constant') then
          write(iulog,*)'co2_type = ',co2_type,' is not supported in offline mode'
          write(iulog,*)'choice is only constant'
          call endrun()
       end if
#endif

       ! Consistency settings for dynamic land use, etc.

       if (fpftdyn /= ' ' .and. create_crop_landunit) then
          write(iulog,*)'dynamic landuse is currently not supported with create_crop_landunit option'
          call endrun()
       end if
       if (fpftdyn /= ' ') then
#if (defined DGVM)
          write(iulog,*)'dynamic landuse is currently not supported with DGVM option'
          call endrun()
#elif (defined CASA)          
          write(iulog,*)'dynamic landuse is currently not supported with CASA option'
          call endrun()
#endif       
       end if

#if (defined SEQ_MCT) || (defined SEQ_ESMF)
       ! Override select set of namelist values with sequential driver input

       if ( .not. present(CCSMInit) )then
          call endrun( subname//' error CCSMInit not present but is '// &
                       'required when linking with CAM' )
       end if
       call shr_inputInfo_initGetData( CCSMInit, case_name=caseid,         &
                                       case_desc=ctitle,                   &
                                       version=version, username=username, &
                                       hostname=hostname,                  &
                                       brnch_retain_casename=brnch_retain_casename,  &
                                       single_column=single_column ,       &
				       scmlat=scmlat,scmlon=scmlon)
       if (      shr_inputInfo_initIsStartup(  CCSMInit ) )then
          nsrest = 0
       else if ( shr_inputInfo_initIsContinue( CCSMInit ) )then
          nsrest = 1
       else if ( shr_inputInfo_initIsBranch(   CCSMInit ) )then
          nsrest = 3
       end if
#if (defined RTM) || (defined DGVM)
       if (is_perpetual()) then
          write(iulog,*)'RTM or DGVM cannot be defined in perpetual mode'
          call endrun()
       end if
#endif
       if (is_perpetual()) then
          if (finidat == ' ') then
             write(iulog,*)'must specify initial dataset for perpetual mode'
             call endrun()
          end if
       end if
#endif

#if (defined RTM)
       ! If rtm_nsteps was not entered in the namelist, give it the following default value

       if (rtm_nsteps == -999) then
          rtm_nsteps = (3600*3)/dtime ! 3 hours
       endif
#endif

       ! Check that hist_type_1d is not set for primary tape

       if (hist_type1d_pertape(1) /= ' ') then
          write(iulog,*)'CONTROL_INIT error: hist_type1d_pertape can only be set for tapes 2-6'
          call endrun()
       end if

       ! Check on run type

       if (nsrest == iundef) then
          write(iulog,*) 'error: must set nsrest'
          call endrun()
       end if
       if (nsrest == 3 .and. nrevsn == ' ') then
          write(iulog,*) 'error: need to set restart data file name'
          call endrun()
       end if

       ! Check on offline mode 

#if (defined OFFLINE)
       if (offline_atmdir == ' ') then
          write(iulog,*)'error: atmos  input data file must be specified'; call endrun()
       end if
       if (cycle_begyr /= iundef) then
          if (cycle_nyr == iundef) then
	     write(iulog,*)'error: if cycle_begyr is set, cycle_nyr must also be set'; call endrun()
	  end if
       end if
       if (cycle_nyr /= iundef) then
          if (cycle_begyr == iundef) then
	     write(iulog,*)'error: if cycle_nyr is set, cycle_begyr must also be set'; call endrun()
	  end if
       end if
#endif

       ! Check on ccsm mode 

#if (defined COUP_CSM)
       if (csm_doflxave .and. irad ==1 ) then
          write(iulog,*)'error: irad must be greater that one if', &
            ' flux averaging option is enabled'
          call endrun
       end if
#endif

       ! Check on nitrogen deposition dataset
       
       if (fndepdat /= ' ' .and. fndepdyn /= ' ') then
          write(iulog,*)'namelist error: only one of fndepdat or fndepdyn can be defined'
          call endrun()
       end if
       
       ! Model physics

       if (irad < 0) irad = nint(-irad*3600._r8/dtime)

       if ( (co2_ppmv <= 0.0_r8) .or. (co2_ppmv > 3000.0_r8) )then
          write(iulog,*)'namelist error: co2_ppmv is out of a reasonable range'
          call endrun()
       end if
       ! History and restart files

       do i = 1, max_tapes
          if (hist_nhtfrq(i) == 0) then
             hist_mfilt(i) = 1
          else if (hist_nhtfrq(i) < 0) then
             hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*dtime))
          endif
       end do
       
       rpntpath = 'rpointer.lnd'
       
       if (nsrest == 0) nrevsn = ' '
       if (nsrest == 1) nrevsn = 'set by restart pointer file file'
       
#if (defined DGVM)
       hist_crtinic = 'YEARLY'
#else
       if (trim(hist_crtinic) /= 'MONTHLY'  .and. trim(hist_crtinic) /= 'YEARLY' .and. &
            trim(hist_crtinic) /= '6-HOURLY' .and. trim(hist_crtinic) /= 'DAILY'  ) then
          hist_crtinic = 'NONE'
       endif
#endif

       ! Split the full pathname of the restart pointer file into a 
       ! directory name and a file name
       
       rpntdir = '.'        ! set path = "."
       rpntfil = rpntpath   ! use whole input string.
100    continue
       
    endif   ! end of if-masterproc if-block

    ! ----------------------------------------------------------------------
    ! Broadcast all control information if appropriate
    ! ----------------------------------------------------------------------

    call control_spmd()
    
    if (masterproc) then
       write(iulog,*) 'Successfully initialized run control settings'
       write(iulog,*)
    endif

  end subroutine control_init


!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_spmd
!
! !INTERFACE:
  subroutine control_spmd()
!
! !DESCRIPTION:
! Distribute namelist data all processors. The cpp SPMD definition
! provides for the funnelling of all program i/o through the master
! processor. Processor 0 either reads restart/history data from the
! disk and distributes it to all processors, or collects data from
! all processors and writes it to disk.
!
! !USES:
!
    use clm_time_manager
#if (defined CASA)
    use CASAMod, only : lnpp, lalloc, q10, spunup
#endif
    use spmdMod, only : mpicom
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer ier       !error code
!-----------------------------------------------------------------------

    ! run control variables

    call mpi_bcast (caseid,   len(caseid),   MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (ctitle,   len(ctitle),   MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (version,  len(version),  MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hostname, len(hostname), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (username, len(username), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nsrest,              1,  MPI_INTEGER  , 0, mpicom, ier)

    call mpi_bcast (dtime    , 1, MPI_INTEGER  , 0, mpicom, ier)

#if (defined OFFLINE) || (defined COUP_CSM)
    call mpi_bcast (nestep   , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (nelapse  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (start_ymd, 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (start_tod, 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (stop_ymd , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (stop_tod , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (ref_ymd  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (ref_tod  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (calendar ,len(calendar), MPI_CHARACTER, 0, mpicom, ier)
#endif

    ! initial file variables

    call mpi_bcast (nrevsn  , len(nrevsn)  , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat , len(finidat) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsurdat , len(fsurdat) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmgrid, len(fatmgrid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmlndfrc,len(fatmlndfrc),MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmtopo, len(fatmtopo) ,MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (flndtopo, len(flndtopo) ,MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fndepdat, len(fndepdat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fndepdyn, len(fndepdyn), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fpftcon , len(fpftcon) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fpftdyn , len(fpftdyn) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (frivinp_rtm, len(frivinp_rtm), MPI_CHARACTER, 0, mpicom, ier)

    ! Landunit generation

    call mpi_bcast(create_crop_landunit, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(allocate_all_vegpfts, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! BGC

    call mpi_bcast (co2_type, len(co2_type), MPI_CHARACTER, 0, mpicom, ier)

    ! physics variables

    call mpi_bcast (irad        , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (csm_doflxave, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (rtm_nsteps  , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nsegspc     , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (wrtdia      , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (single_column,1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (scmlat,       1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (scmlon,       1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (co2_ppmv    , 1, MPI_REAL8,   0, mpicom, ier)

    ! history file variables

    call mpi_bcast (hist_empty_htapes, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, size(hist_ndens), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_crtinic, len(hist_crtinic), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_avgflag_pertape, size(hist_avgflag_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_type1d_pertape, max_namlen*size(hist_type1d_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl1, max_namlen*size(hist_fexcl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl2, max_namlen*size(hist_fexcl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl3, max_namlen*size(hist_fexcl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl4, max_namlen*size(hist_fexcl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl5, max_namlen*size(hist_fexcl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl6, max_namlen*size(hist_fexcl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl1, (max_namlen+2)*size(hist_fincl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl2, (max_namlen+2)*size(hist_fincl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl3, (max_namlen+2)*size(hist_fincl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl4, (max_namlen+2)*size(hist_fincl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl5, (max_namlen+2)*size(hist_fincl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl6, (max_namlen+2)*size(hist_fincl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (rest_flag, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! restart file variables

    call mpi_bcast (rpntpath, len(rpntpath), MPI_CHARACTER, 0, mpicom, ier)

    ! clump decomposition variables

    call mpi_bcast (clump_pproc, 1, MPI_INTEGER, 0, mpicom, ier)

    ! error growth perturbation limit
    call mpi_bcast (pertlim, 1, MPI_REAL8, 0, mpicom, ier)

#if (defined CASA)
    call mpi_bcast (lnpp  , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (lalloc, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (spunup, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (q10   , 1, MPI_REAL8  , 0, mpicom, ier)
#endif

  end subroutine control_spmd

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_print
!
! !INTERFACE:
  subroutine control_print ()
!
! !DESCRIPTION:
! Write out run control variables
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i  !loop index
!------------------------------------------------------------------------

    write(iulog,*) 'define run:'
    write(iulog,*) '   source                = ',trim(source)
    write(iulog,*) '   version               = ',trim(version)
    write(iulog,*) '   run type              = ',runtyp(nsrest+1)
    write(iulog,*) '   case title            = ',trim(ctitle)
    write(iulog,*) '   username              = ',trim(username)
    write(iulog,*) '   hostname              = ',trim(hostname)
    write(iulog,*) 'input data files:'
    write(iulog,*) '   PFT physiology = ',trim(fpftcon)
    if (fsurdat == ' ') then
       write(iulog,*) '   fsurdat, surface dataset not set'
    else
       write(iulog,*) '   surface data   = ',trim(fsurdat)
    end if
    if (flndtopo == ' ') then
       write(iulog,*) '   flndtopo not set'
    else
       write(iulog,*) '   land topographic data = ',trim(flndtopo)
    end if
    if (fatmgrid == ' ') then
       write(iulog,*) '   fatmgrid not set, using fsurdat'
       fatmgrid = fsurdat
       write(iulog,*) '   atm grid data  = ',trim(fatmgrid)
    else
       write(iulog,*) '   atm grid data  = ',trim(fatmgrid)
    end if
    if (fatmlndfrc == ' ') then
       write(iulog,*) '   fatmlndfrc not set, using fatmgrid'
       fatmlndfrc = fatmgrid
       write(iulog,*) '   land frac data = ',trim(fatmlndfrc)
    else
       write(iulog,*) '   land frac data = ',trim(fatmlndfrc)
    end if
    if (fatmtopo == ' ') then
       write(iulog,*) '   fatmtopo not set'
    else
       write(iulog,*) '   atm topographic data = ',trim(fatmtopo)
    end if
    if (fndepdat == ' ') then
        write(iulog,*) '   NOT using input data for nitrogen deposition'
    else
        write(iulog,*) '   nitrogen deposition data = ',trim(fndepdat)
    endif
    if (fndepdyn == ' ') then
        write(iulog,*) '   NOT using dynamic input data for nitrogen deposition'
    else
        write(iulog,*) '   dynamic nitrogen deposition data = ',trim(fndepdyn)
    endif
    if (nsrest == 0 .and. finidat == ' ') write(iulog,*) '   initial data created by model'
    if (nsrest == 0 .and. finidat /= ' ') write(iulog,*) '   initial data   = ',trim(finidat)
    if (nsrest /= 0) write(iulog,*) '   restart data   = ',trim(nrevsn)
#if (defined OFFLINE)
    if (offline_atmdir /= ' ') then
       write(iulog,*) '   atmospheric forcing data    = ',trim(offline_atmdir)
    end if
    if (cycle_begyr /= -9999999) then
       write(iulog,*) '   first year of atmospheric forcing data for cycling = ',cycle_begyr
       write(iulog,*) '   number of years in atmospheric forcing data cycle  = ',cycle_nyr
    end if
#elif (defined SEQ_MCT) || (defined SEQ_ESMF)
    write(iulog,*) '   atmospheric forcing data is from sequential ccsm model'
#elif (defined COUP_CSM)
    write(iulog,*) '   atmospheric forcint data is from ccsm flux coupler'
#endif
#if (defined RTM)
    if (frivinp_rtm /= ' ') write(iulog,*) '   RTM river data       = ',trim(frivinp_rtm)
#endif
    write(iulog,*) 'Restart parameters:'
    write(iulog,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(iulog,*)'   restart pointer file name          = ',trim(rpntfil)
    if (hist_crtinic == 'MONTHLY') then
       write(iulog,*)'initial datasets will be written monthly'
    else if (hist_crtinic == 'YEARLY') then
       write(iulog,*)'initial datasets will be written yearly'
    else if (hist_crtinic == 'DAILY') then
       write(iulog,*)'initial datasets will be written daily'
    else if (hist_crtinic == '6-HOURLY') then
       write(iulog,*)'initial datasets will be written 6-hourly'
    else
       write(iulog,*)'initial datasets will not be produced'
    endif
    write(iulog,*) 'model physics parameters:'
#if (defined PERGRO)
    write(iulog,*) '   flag for random perturbation test is set'
#else
    write(iulog,*) '   flag for random perturbation test is not set'
#endif
    write(iulog,*) '   solar radiation frequency (iterations) = ',irad
    write(iulog,*) '   CO2 volume mixing ratio   (umol/mol)   = ', co2_ppmv
#if (defined COUP_CSM)
    write(iulog,*) 'communication with the flux coupler'
    if (csm_doflxave) then
       write(iulog,*)'    data will be sent to the flux coupler ', &
            'only when an albedo calculation is performed '
       write(iulog,*)'     fluxes will be averaged on steps where ', &
            'communication with the flux coupler does not occur'
    else
       write(iulog,*)'    data will be sent and received to/from ', &
            'the flux coupler at every time step except for nstep=1'
    endif
#endif
#if (defined RTM)
    if (rtm_nsteps > 1) then
       write(iulog,*)'river runoff calculation performed only every ',rtm_nsteps,' nsteps'
    else
       write(iulog,*)'river runoff calculation performed every time step'
    endif
#endif
    if (nsrest == 1) then
       write(iulog,*) 'restart warning:'
       write(iulog,*) '   Namelist not checked for agreement with initial run.'
       write(iulog,*) '   Namelist should not differ except for ending time step and run type'
    end if
    if (nsrest == 3) then
       write(iulog,*) 'branch warning:'
       write(iulog,*) '   Namelist not checked for agreement with initial run.'
       write(iulog,*) '   Surface data set and reference date should not differ from initial run'
    end if
#if (defined COUP_CSM)
    write(iulog,*) '   last time step determined by flux coupler'
#endif
    if ( pertlim /= 0.0_r8 ) &
    write(iulog,*) '   perturbation limit = ',pertlim
    write(iulog,*) '   maxpatch_pft         = ',maxpatch_pft
    write(iulog,*) '   allocate_all_vegpfts = ',allocate_all_vegpfts
    write(iulog,*) '   nsegspc              = ',nsegspc

  end subroutine control_print

end module controlMod
