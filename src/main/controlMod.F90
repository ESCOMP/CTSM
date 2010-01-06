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
! Note: For definitions of namelist variablses see
!       ../../bld/namelist_files/namelist_definition.xml
!       Display the file in a browser to see it neatly formatted in html.
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, SHR_KIND_CL
  use clm_varpar   , only : maxpatch_pft
  use clm_varctl   , only : caseid, ctitle, nsrest, brnch_retain_casename, hostname, model_version=>version,    &
                            iulog, hist_crtinic, outnc_large_files, finidat, fsurdat, fatmgrid, fatmlndfrc,     &
                            fatmtopo, flndtopo, fndepdat, fndepdyn, fpftdyn, fpftcon, nrevsn, frivinp_rtm,      &
                            create_crop_landunit, allocate_all_vegpfts, fget_archdev, &
                            co2_type, wrtdia, co2_ppmv, rtm_nsteps, nsegspc, pertlim,       &
                            hist_pioflag, ncd_lowmem2d, ncd_pio_def, ncd_pio_UseRearranger, username,           &
                            ncd_pio_UseBoxRearr, ncd_pio_SerialCDF, ncd_pio_IODOF_rootonly, ncd_pio_DebugLevel, &
                            ncd_pio_num_iotasks, fsnowaging, fsnowoptics, &
                            faerdep
  use spmdMod      , only : masterproc
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
  use UrbanMod     , only : urban_hac, urban_traffic
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
!
! !PRIVATE TYPES:
! Namelist variables only used locally
  character(len=  7) :: runtyp(4)                        ! run type
  character(len=SHR_KIND_CL) :: NLFilename = 'lnd.stdin' ! Namelist filename
#if (defined _OPENMP)
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
                                             ! concurrently in a single parallel region
#endif
!EOP
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
!
! !LOCAL VARIABLES:
!EOP
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
  subroutine control_init( )
!
! !DESCRIPTION:
! Initialize CLM run control information
!
! !USES:
    use clm_time_manager , only : set_timemgr_init, is_perpetual, get_timemgr_defaults
#if (defined CASA)
    use CASAMod          , only : lnpp, lalloc, q10, spunup, fcpool
#endif
    use fileutils        , only : getavu, relavu
    use shr_string_mod   , only : shr_string_getParentDir
    use clm_varctl       , only : clmvarctl_init

    implicit none
!
    include 'netcdf.inc'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=32)  :: starttype ! infodata start type
    integer :: i,j,n                ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'control_init'  ! subroutine name
!------------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! clm input datasets

    integer :: dtime    ! Integer time-step
    namelist / clm_inparm/ &
	 dtime	

    namelist /clm_inparm/  &
         finidat, fsurdat, fatmgrid, fatmlndfrc, fatmtopo, flndtopo, &
         fpftcon, frivinp_rtm,  &
         fpftdyn, fndepdat, fndepdyn, nrevsn, &
         fsnowoptics, fsnowaging

    namelist /clm_inparm/ faerdep

    ! clm history, restart options

    namelist /clm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, hist_mfilt, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6, &
         hist_crtinic, rest_flag, outnc_large_files, &
         hist_pioflag, ncd_lowmem2d, ncd_pio_def, &
         ncd_pio_UseRearranger, ncd_pio_UseBoxRearr, ncd_pio_SerialCDF, &
         ncd_pio_IODOF_rootonly, ncd_pio_DebugLevel, ncd_pio_num_iotasks

    ! clm bgc info

#if (defined CASA)
    namelist /clm_inparm/  &
         lnpp, lalloc, q10, spunup, fcpool
#endif

    namelist /clm_inparm / &
         co2_type

    ! clm other options

    namelist /clm_inparm/  &
         clump_pproc, wrtdia, rtm_nsteps, pertlim, &
         create_crop_landunit, nsegspc, co2_ppmv

    ! clm urban options

    namelist /clm_inparm/  &
         urban_hac, urban_traffic
         
    ! ----------------------------------------------------------------------
    ! Default values
    ! ----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to initialize run control settings .....'
    endif

    runtyp(0 + 1) = 'initial'
    runtyp(1 + 1) = 'restart'
    runtyp(3 + 1) = 'branch '

#if (defined CASA)
    lnpp = 2
    lalloc = 1
    q10 = 2.0_r8          ! set Q10 to 2.0  03/11/19
    spunup = 0
    fcpool = ' '
#endif

    ! Set clumps per procoessor

#if (defined _OPENMP)
    clump_pproc = omp_get_max_threads()
#else
    clump_pproc = 1
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

       call set_timemgr_init( dtime_in=dtime )

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

       if (urban_traffic) then
          write(iulog,*)'Urban traffic fluxes are not implemented currently'
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
       
    endif   ! end of if-masterproc if-block

    call clmvarctl_init( masterproc, dtime )

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
! Distribute namelist data all processors. All program i/o is 
! funnelled through the master processor. Processor 0 either 
! reads restart/history data from the disk and distributes 
! it to all processors, or collects data from
! all processors and writes it to disk.
!
! !USES:
!
#if (defined CASA)
    use CASAMod,    only : lnpp, lalloc, q10, spunup
#endif
    use spmdMod,    only : mpicom, MPI_CHARACTER, MPI_INTEGER, MPI_LOGICAL, MPI_REAL8
    use clm_varctl, only : single_column, scmlat, scmlon, rpntfil
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier       !error code
!-----------------------------------------------------------------------

    ! run control variables

    call mpi_bcast (caseid,         len(caseid),        MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (ctitle,         len(ctitle),        MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (model_version,  len(model_version), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hostname,       len(hostname),      MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (username,       len(username),      MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nsrest,                     1,      MPI_INTEGER  , 0, mpicom, ier)

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
#if (defined RTM)
    call mpi_bcast (frivinp_rtm, len(frivinp_rtm), MPI_CHARACTER, 0, mpicom, ier)
#endif
    call mpi_bcast (fsnowoptics,  len(fsnowoptics),  MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsnowaging,   len(fsnowaging),   MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fget_archdev, len(fget_archdev), MPI_CHARACTER, 0, mpicom, ier)
    
    call mpi_bcast (faerdep,      len(faerdep),      MPI_CHARACTER, 0, mpicom, ier)

    ! Landunit generation

    call mpi_bcast(create_crop_landunit, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(allocate_all_vegpfts, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! BGC

    call mpi_bcast (co2_type, len(co2_type), MPI_CHARACTER, 0, mpicom, ier)

    ! physics variables

    call mpi_bcast (urban_hac     , len(urban_hac), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (urban_traffic , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (rtm_nsteps  , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nsegspc     , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (wrtdia      , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (single_column,1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (scmlat,       1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (scmlon,       1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (co2_ppmv    , 1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (hist_pioflag, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ncd_lowmem2d, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ncd_pio_def , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ncd_pio_UseRearranger , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ncd_pio_UseBoxRearr   , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ncd_pio_SerialCDF     , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ncd_pio_IODOF_rootonly, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (ncd_pio_DebugLevel    , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (ncd_pio_num_iotasks   , 1, MPI_INTEGER, 0, mpicom, ier)

    ! history file variables

    call mpi_bcast (outnc_large_files, 1, MPI_LOGICAL, 0, mpicom, ier)
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

    call mpi_bcast (rpntfil, len(rpntfil), MPI_CHARACTER, 0, mpicom, ier)

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
! Write out the clm namelist run control variables
!
! !USES:
!
    use clm_varctl, only : source, rpntdir, rpntfil
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer i  !loop index
!------------------------------------------------------------------------

    write(iulog,*) 'define run:'
    write(iulog,*) '   source                = ',trim(source)
    write(iulog,*) '   model_version         = ',trim(model_version)
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
    if (fsnowoptics == ' ') then
       write(iulog,*) '   snow optical properties file NOT set'
    else
       write(iulog,*) '   snow optical properties file = ',trim(fsnowoptics)
    endif
    if (fsnowaging == ' ') then
       write(iulog,*) '   snow aging parameters file NOT set'
    else
       write(iulog,*) '   snow aging parameters file = ',trim(fsnowaging)
    endif
    if (faerdep == ' ') then
       write(iulog,*) '   aerosol deposition file NOT set'
    else
       write(iulog,*) '   aerosol deposition file = ',trim(faerdep)
    endif

    if (nsrest == 0 .and. finidat == ' ') write(iulog,*) '   initial data created by model'
    if (nsrest == 0 .and. finidat /= ' ') write(iulog,*) '   initial data   = ',trim(finidat)
    if (nsrest /= 0) write(iulog,*) '   restart data   = ',trim(nrevsn)
    write(iulog,*) '   atmospheric forcing data is from sequential ccsm model'
#if (defined RTM)
    if (frivinp_rtm /= ' ') write(iulog,*) '   RTM river data       = ',trim(frivinp_rtm)
#endif
    write(iulog,*) 'Restart parameters:'
    write(iulog,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(iulog,*)'   restart pointer file name          = ',trim(rpntfil)
    if ( outnc_large_files ) then
       write(iulog,*)'Large file support for output files is ON'
    end if
    if ( trim(fget_archdev) /= "null:" ) then
       write(iulog,*)'try to retreive input files that do NOT exist from archival device: ', trim(fget_archdev)
    end if
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
    write(iulog,*) '   CO2 volume mixing ratio   (umol/mol)   = ', co2_ppmv
    write(iulog,*) '   urban air conditioning/heating and wasteheat   = ', urban_hac
    write(iulog,*) '   urban traffic flux   = ', urban_traffic
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
    if ( pertlim /= 0.0_r8 ) &
    write(iulog,*) '   perturbation limit = ',pertlim
    write(iulog,*) '   maxpatch_pft         = ',maxpatch_pft
    write(iulog,*) '   allocate_all_vegpfts = ',allocate_all_vegpfts
    write(iulog,*) '   nsegspc              = ',nsegspc

!tcx for debugging
    write(iulog,*) 'history/PIO parameters:'
    write(iulog,*) '   hist_pioflag           = ', hist_pioflag
    write(iulog,*) '   ncd_lowmem2d           = ', ncd_lowmem2d
    write(iulog,*) '   ncd_pio_def            = ', ncd_pio_def
    write(iulog,*) '   ncd_pio_UseRearranger  = ', ncd_pio_UseRearranger
    write(iulog,*) '   ncd_pio_UseBoxRearr    = ', ncd_pio_UseBoxRearr
    write(iulog,*) '   ncd_pio_SerialCDF      = ', ncd_pio_SerialCDF
    write(iulog,*) '   ncd_pio_IODOF_rootonly = ', ncd_pio_IODOF_rootonly
    write(iulog,*) '   ncd_pio_DebugLevel     = ', ncd_pio_DebugLevel
    write(iulog,*) '   ncd_pio_num_iotasks    = ', ncd_pio_num_iotasks

  end subroutine control_print

end module controlMod
