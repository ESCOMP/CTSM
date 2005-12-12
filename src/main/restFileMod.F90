#include <misc.h>
#include <preproc.h>

module restFileMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: restFileMod
!
! !DESCRIPTION:
! Reads from or writes to/ the CLM restart file.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils,   only : endrun
  use ncdio       
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: restFile_read
  public :: restFile_write
  public :: restFile_read_binary
  public :: restFile_write_binary
  public :: restFile_open
  public :: restFile_close
  public :: restFile_getfile
  public :: restFile_filename       ! Sets restart filename
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: restFile_read_pfile     
  private :: restFile_write_pfile    ! Writes restart pointer file
  private :: restFile_archive        ! Close restart file and write restart pointer file
  private :: restFile_dimset
  private :: restFile_dimcheck
  private :: restFile_enddef
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !PRIVATE TYPES:
  private
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write
!
! !INTERFACE:
  subroutine restFile_write( file )
!
! !DESCRIPTION:
! Read/write CLM restart file.
!
! !USES:
#if (defined OFFLINE) || (defined COUP_CSM)
    use time_manager     , only : timemgr_restart, get_nstep
#else
    use time_manager     , only : get_nstep
#endif
    use subgridRestMod   , only : SubgridRest
    use BiogeophysRestMod, only : BiogeophysRest
#if (defined CN)
    use CNRestMod        , only : CNRest
#endif
#if (defined DGVM)
    use DGVMRestMod      , only : DGVMRest
#endif
#if (defined RTM)
    use RtmMod           , only : RTMRest
#endif
#if (defined CASA)
    use CASAMod          , only : CASARest
#endif
#if (defined COUP_CSM)
    use clm_csmMod       , only : csm_restart 
#endif
    use accumulMod       , only : accumulRest
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: file  ! output netcdf restart file
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ncid    ! netcdf id
    integer :: i       ! index
!-----------------------------------------------------------------------

    ! Open restart file

    call restFile_open( flag='write', file=file, ncid=ncid )

    ! Define dimensions

    call restFile_dimset ( ncid )

    ! Define restart file variables

#if (defined OFFLINE) || (defined COUP_CSM)
    call timemgr_restart(ncid, flag='define')
#endif
    call SubgridRest( ncid, flag='define' )
    call BiogeophysRest( ncid, flag='define' )
#if (defined CN)
    call CNRest( ncid, flag='define' )
#endif
#if (defined CASA)
    call CASARest( ncid, flag='define' )
#endif
#if (defined DGVM)
    call DGVMRest( ncid, flag='define' )
#endif
#if (defined RTM)
    call RtmRest( ncid, flag='define' )
#endif
#if (defined COUP_CSM)
    call csm_restart(ncid, flag='define')
#endif
    call accumulRest( ncid, flag='define' )
    call restFile_enddef( ncid )

    ! Write restart file variables
    
#if (defined OFFLINE) || (defined COUP_CSM)
    call timemgr_restart(ncid, flag='write')
#endif
    call SubgridRest( ncid, flag='write' )
    call BiogeophysRest( ncid, flag='write' )
#if (defined CN)
    call CNRest( ncid, flag='write' )
#endif
#if (defined CASA)
    call CASARest( ncid, flag='write' )
#endif
#if (defined DGVM)
    call DGVMRest( ncid, flag='write' )
#endif
#if (defined RTM)
    call RtmRest( ncid, flag='write' )
#endif
#if (defined COUP_CSM)
    call csm_restart(ncid, flag='write')
#endif
    call accumulRest( ncid, flag='write' )
    
    ! Close and archive restart file
    
    call restFile_close( ncid )
    call restFile_archive( file )
    
    ! Write restart pointer file
    
    call restFile_write_pfile( file )
    
    ! Write out diagnostic info

    if (masterproc) then
       write (6,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
       write (6,'(72a1)') ("-",i=1,60)
    end if
    
  end subroutine restFile_write

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read
!
! !INTERFACE:
  subroutine restFile_read( file )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
    use BiogeophysRestMod, only : BiogeophysRest
#if (defined CN)
    use CNRestMod        , only : CNRest
#endif
#if (defined DGVM)
    use DGVMRestMod      , only : DGVMRest
#endif
#if (defined RTM)
    use RtmMod           , only : RTMRest
#endif
#if (defined CASA)
    use CASAMod          , only : CASARest
#endif
#if (defined COUP_CSM)
    use clm_csmMod       , only : csm_restart 
#endif
    use accumulMod       , only : accumulRest
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: file  ! output netcdf restart file
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ncid                        ! netcdf id
    integer :: i                           ! index
!-----------------------------------------------------------------------

    ! Open file

    call restFile_open( flag='read', file=file, ncid=ncid )

    ! Read file

    call restFile_dimcheck( ncid )
    call BiogeophysRest( ncid, flag='read' )
#if (defined CN)
    call CNRest( ncid, flag='read' )
#endif
#if (defined CASA)
    call CASARest( ncid, flag='read' )
#endif
#if (defined DGVM)
    call DGVMRest( ncid, flag='read' )
#endif
#if (defined RTM)
    call RtmRest( ncid, flag='read' )
#endif
#if (defined COUP_CSM)
    call csm_restart( ncid, flag='read' )
#endif
    call accumulRest( ncid, flag='read' )
    
    ! Close file 

    call restFile_close( ncid )

    ! Write out diagnostic info

    if (masterproc) then
       write(6,'(72a1)') ("-",i=1,60)
       write(6,*) 'Successfully read restart data for restart run'
       write(6,*)
    end if

  end subroutine restFile_read

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read_binary
!
! !INTERFACE:
  subroutine restFile_read_binary( file )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
    use fileutils  , only : relavu, opnfil, getfil, getavu
    use histFileMod, only : restart_history 
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file   ! binary restart file
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nio                         ! Fortran unit number
!-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       call opnfil (file, nio, 'u')
    end if
    call restart_history(nio, flag='read')
    if (masterproc) then
       call relavu (nio)
    end if

  end subroutine restFile_read_binary

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write_binary
!
! !INTERFACE:
  subroutine restFile_write_binary( file )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
    use fileutils  , only : relavu, opnfil, getfil, getavu
    use histFileMod, only : restart_history 
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file   ! binary restart file
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nio                         ! Fortran unit number
!-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       call opnfil (file, nio, 'u')
    end if
    call restart_history(nio, flag='write')
    if (masterproc) then
       call relavu (nio)
       call restFile_archive( file )
    end if

  end subroutine restFile_write_binary

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_getfile
!
! !INTERFACE:
  subroutine restFile_getfile( file, path )
!
! !DESCRIPTION:
! Determine and obtain netcdf restart file
!
! !USES:
    use clm_varctl, only : caseid, finidat, nrevsn, nsrest, brnch_retain_casename
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: status                      ! return status
    integer :: length                      ! temporary          
!-----------------------------------------------------------------------

    if (masterproc) then

       ! Restart run:
       ! Restart file pathname is read restart pointer file 

       if (nsrest==1) then
          call restFile_read_pfile( path )
          call getfil( path, file, 0 )
       end if
       
       ! Branch run: 
       ! Restart file pathname is obtained from namelist "nrevsn"
       ! Check case name consistency (case name must be different for branch run, 
       ! unless namelist specification states otherwise)

       if (nsrest==3) then
          length = len_trim(nrevsn)
          if (nrevsn(length-2:length) == '.nc') then
             path = trim(nrevsn) 
          else
             path = trim(nrevsn) // '.nc'
          end if
          call getfil( path, file, 0 )

          status = verify(caseid, file)
          if (status == 0 .and. .not.(brnch_retain_casename)) then
             write(6,*) 'Must change case name on branch run if ',&
                  'brnch_retain_casename namelist is not set'
             write(6,*) 'previous case filename= ',trim(file),&
                  ' current case = ',trim(caseid)
             call endrun()
          end if
       end if

       ! Initial run: 
       ! Restart file pathname is obtained from namelist "finidat"

       if (nsrest==0) then
          call getfil( finidat, file, 0 )
       end if

    end if

  end subroutine restFile_getfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read_pfile
!
! !INTERFACE:
  subroutine restFile_read_pfile( pnamer )
!
! !DESCRIPTION:
! Setup restart file and perform necessary consistency checks
!
! !USES:
    use fileutils , only : opnfil, getavu, relavu
    use clm_varctl, only : rpntfil, rpntdir
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: pnamer ! full path of binary restart file
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
!-----------------------------------------------------------------------

    ! Obtain the restart file from the restart pointer file. 
    ! For restart runs, the restart pointer file contains the full pathname 
    ! of the restart file. For branch runs, the namelist variable 
    ! [nrevsn] contains the full pathname of the restart file. 
    ! New history files are always created for branch runs.
       
    if (masterproc) then
       write (6,*) 'Reading restart pointer file....'
       nio = getavu()
       locfn = trim(rpntdir) //'/'// trim(rpntfil)
       call opnfil (locfn, nio, 'f')
       read (nio,'(a256)') pnamer
       call relavu (nio)
       write (6,*) 'Reading restart data.....'
       write (6,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_read_pfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_archive
!
! !INTERFACE:
  subroutine restFile_archive( file )
!
! !DESCRIPTION:
! Close and archive restart file and write restart pointer file if
! in write mode, otherwise just close restart file if in read mode
!
! !USES:
#if (defined COUP_CSM)
    use clm_csmMod  , only : csmstop_next
#endif
    use time_manager, only : is_last_step
    use clm_varctl  , only : mss_irt, mss_wpass, archive_dir
    use fileutils   , only : putfil, set_filename
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file  ! local output filename
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i                   !index
    logical :: lremove             !true => remove file after archive
    character(len=256) :: rem_fn   !remote (archive) filename
    character(len=256) :: rem_dir  !remote (archive) directory
!-----------------------------------------------------------------------

   if (masterproc) then

      lremove = .true.
#if (defined OFFLINE) || (defined COUP_CAM)
      if (is_last_step()) lremove = .false.
#elif (defined COUP_CSM)
      if (csmstop_next) lremove = .false.
#endif
      if (mss_irt > 0) then
         rem_dir = trim(archive_dir) // '/rest/'
         rem_fn = set_filename(rem_dir, file)
         call putfil( file, rem_fn, mss_wpass, mss_irt, lremove )
      endif
      write(6,*) 'Successfully wrote local restart file ',trim(file)
      write(6,'(72a1)') ("-",i=1,60)
      write(6,*)

   end if

 end subroutine restFile_archive

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write_pfile
!
! !INTERFACE:
  subroutine restFile_write_pfile( fnamer )
!
! !DESCRIPTION:
! Open restart pointer file. Write names of current binary and netcdf
! restart files. If using mass store, these are the mass store names
! except if mss_irt=0 (no mass store files written). 
!
! !USES:
    use clm_varctl, only : rpntdir, mss_irt, archive_dir, rpntfil
    use fileutils , only : set_filename, relavu
    use fileutils , only : getavu, opnfil
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fnamer
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
    character(len=256) :: rem_dir   ! remote directory
!-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)
       call opnfil( filename, nio, 'f' )
       
       if (mss_irt == 0) then
          write(nio,'(a)') fnamer
       else
          rem_dir = trim(archive_dir) // '/rest/'
          write(nio,'(a)') set_filename( rem_dir, fnamer )
       endif
       
       call relavu( nio )
       write(6,*)'Successfully wrote local restart pointer file'
    end if

  end subroutine restFile_write_pfile

!-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

    use time_manager, only : get_nstep
    
    implicit none
    character(len=*), intent(in) :: flag ! flag to specify read or write
    character(len=*), intent(in) :: file ! filename
    integer, intent(out)         :: ncid ! netcdf id

    integer :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (masterproc) then
       if (flag == 'write') then

          ! Create new netCDF file (in define mode) and set fill mode
          ! to "no fill" to optimize performance

          write(6,*)
          write(6,*)'restFile_open: writing restart dataset at ',&
               trim(file), 'at nstep = ',get_nstep()
          write(6,*)
          call check_ret( nf_create(trim(file), nf_clobber, ncid), subname )
          call check_ret( nf_set_fill(ncid, nf_nofill, omode), subname )

       else if (flag == 'read') then
       
          ! Open netcdf restart file

          write (6,*) 'Reading restart dataset'
          call check_ret( nf_open(file, nf_nowrite, ncid), subname )

       end if
    end if
  
  end subroutine restFile_open

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_filename
!
! !INTERFACE:
  character(len=256) function restFile_filename( type, offset )
!
! !DESCRIPTION:
!
! !USES:
    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date, get_step_size
!
! !ARGUMENTS:
    implicit none
    character(*), intent(in) :: type         ! output type "binary" or "netcdf"
    integer, optional, intent(in) :: offset  ! offset from current time in seconds
                                             ! positive for future times and 
                                             ! negative for previous times
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       ! date char string
    integer :: day                    ! day (1 -> 31)
    integer :: mon                    ! month (1 -> 12)
    integer :: yr                     ! year (0 -> ...)
    integer :: sec                    ! seconds into current day
!-----------------------------------------------------------------------

    ! Note - the only difference between a restart and an initial file
    ! is that an initial file is written one time step before the date
    ! stamp associated with the file name. Consequently it can be used
    ! for an initial run and have "restart" type of capabilities for that run

    if (masterproc) then
       if (present(offset)) then
          call get_curr_date (yr, mon, day, sec, offset=offset)
       else
          call get_curr_date (yr, mon, day, sec)
       end if
       write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
       
       if (trim(type) == 'binary') then
          restFile_filename = "./"//trim(caseid)//".clm2.r."//trim(cdate)
       else if (trim(type) == 'netcdf') then
          if (present(offset)) then
             restFile_filename = "./"//trim(caseid)//".clm2.i."//trim(cdate)//".nc"
             write(6,*)'writing initial file ',trim(restFile_filename),' for model date = ',cdate
          else
             restFile_filename = "./"//trim(caseid)//".clm2.r."//trim(cdate)//".nc"
             write(6,*)'writing restart file ',trim(restFile_filename),' for model date = ',cdate
          end if
       else
          write(6,*)'restart file type ',trim(type),' is not supported'; call endrun()
       end if
    else
      restfile_filename = 'not_defined'
    end if
    
  end function restFile_filename

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_dimset
!
! !INTERFACE:
  subroutine restFile_dimset( ncid )
!
! !DESCRIPTION:
! Read/Write initial data from/to netCDF instantaneous initial data file
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use time_manager, only : get_nstep, get_curr_date
#ifdef SPMD
    use spmdMod     , only : mpicom, MPI_LOGICAL
#endif
    use clm_varctl  , only : caseid, ctitle, version, allocate_all_vegpfts, fsurdat
    use clm_varpar  , only : numrad, rtmlon, rtmlat, nlevlak, nlevsno, nlevsoi
    use shr_sys_mod , only : shr_sys_getenv
    use decompMod   , only : get_proc_bounds, get_proc_global
#ifdef RTM
    use RunoffMod   , only : get_proc_rof_global
#endif
#if (defined CASA)
  use CASAMod       , only : nlive, npools
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid           ! netCDF dataset id
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: yr                  ! current year (0 -> ...)
    integer :: mon                 ! current month (1 -> 12)
    integer :: day                 ! current day (1 -> 31)
    integer :: mcsec               ! seconds of current date
    integer :: mcdate              ! current date
    integer :: dimid               ! netCDF dimension id
    integer :: numg                ! total number of gridcells across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: nump                ! total number of pfts across all processors
    integer :: nrof_lnd            ! total number of land runoff points across all procs
    integer :: nrof_ocn            ! total number of ocean runoff points across all procs
    integer :: ier                 ! error status
    integer :: strlen_dimid        ! string dimension id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: str
    character(len= 32) :: subname='restFile_dimset' ! subroutine name
!------------------------------------------------------------------------

    if (masterproc) then

       call get_proc_global(numg, numl, numc, nump)
#if (defined RTM)
       call get_proc_rof_global(nrof_lnd, nrof_ocn)
#endif

       ! Define dimensions

       call check_ret( nf_def_dim(ncid, 'gridcell', numg           , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'landunit', numl           , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'column'  , numc           , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'pft'     , nump           , dimid), subname )
       
       call check_ret( nf_def_dim(ncid, 'levsoi'  , nlevsoi        , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'levlak'  , nlevlak        , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'levsno'  , nlevsno        , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'levtot'  , nlevsno+nlevsoi, dimid), subname )
       call check_ret( nf_def_dim(ncid, 'numrad'  , numrad         , dimid), subname )
#if (defined CASA)
       call check_ret(nf_def_dim (ncid, 'nlive'   , nlive          , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'npools'  , npools         , dimid), subname)
#endif
#if (defined RTM)
       call check_ret( nf_def_dim(ncid, 'ocnrof'  , nrof_ocn       , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'lndrof'  , nrof_lnd       , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'rtmlon'  , rtmlon         , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'rtmlat'  , rtmlat         , dimid), subname )
#endif
       call check_ret( nf_def_dim(ncid, 'string_length', 64        , dimid), subname)
       
       ! Define global attributes
       
       str = 'CF-1.0'
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'conventions', len_trim(str), trim(str)), subname)

       call getdatetime(curdate, curtime)
       str = 'created on ' // curdate // ' ' // curtime
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'history', len_trim(str), trim(str)), subname)

       call shr_sys_getenv ('LOGNAME', str, ier)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'logname', len_trim(str), trim(str)), subname)

       call shr_sys_getenv ('HOST', str, ier)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'host', len_trim(str), trim(str)), subname)

       str = 'Community Land Model: CLM3'
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'source', len_trim(str), trim(str)), subname)

       str = '$Name$'
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'version', len_trim(str), trim(str)), subname)

       str = '$Id$'
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'revision_id', len_trim(str), trim(str)), subname)

       str = ctitle
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'case_title', len_trim(str), trim(str)), subname)

       str = caseid
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'case_id', len_trim(str), trim(str)), subname)

       str = fsurdat
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'surface_dataset', len_trim(str), trim(str)), subname)
    end if

  end subroutine restFile_dimset
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_dimcheck
!
! !INTERFACE:
  subroutine restFile_dimcheck( ncid )
!
! !DESCRIPTION:
! Check dimensions of restart file
!
! !USES:
    use decompMod, only : get_proc_bounds, get_proc_global
#if (defined RTM)
    use RunoffMod, only : get_proc_rof_global
#endif
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: ncid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: numg     ! total number of gridcells across all processors
    integer :: numl     ! total number of landunits across all processors
    integer :: numc     ! total number of columns across all processors
    integer :: nump     ! total number of pfts across all processors
    integer :: nrof_lnd ! total number of land runoff points across all procs
    integer :: nrof_ocn ! total number of ocean runoff points across all procs
    character(len=32) :: subname='restFile_dimcheck' ! subroutine name
!-----------------------------------------------------------------------

    ! Get relevant sizes

    if (masterproc) then
       call get_proc_global(numg, numl, numc, nump)
#if (defined RTM)
       call get_proc_rof_global(nrof_lnd, nrof_ocn)
#endif

#if ( !defined SCAM )

       call check_dim(ncid, 'gridcell', numg)
       call check_dim(ncid, 'landunit', numl)
       call check_dim(ncid, 'column'  , numc)
       call check_dim(ncid, 'pft'     , nump)
       call check_dim(ncid, 'levsno'  , nlevsno)
       call check_dim(ncid, 'levsoi'  , nlevsoi)
       call check_dim(ncid, 'levlak'  , nlevlak) 
#if (defined CASA)
       ! Dimensions should be checked, but this will only work for initial
       ! datasets created with CASA enabled so do not normally do this.
       ! call check_dim(ncid, 'nlive'   , nlive)
       ! call check_dim(ncid, 'npools'  , npools)
#endif

#endif
    end if

  end subroutine restFile_dimcheck

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_enddef
!
! !INTERFACE:
  subroutine restFile_enddef( ncid )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nio    ! Fortran unit number
    character(len=32) :: subname='restFile_enddef' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       call check_ret(nf_enddef(ncid), subname)
    end if

  end subroutine restFile_enddef

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_close
!
! !INTERFACE:
  subroutine restFile_close( ncid )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=32) :: subname='restFile_close' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       call check_ret(nf_close(ncid), subname)
    end if

  end subroutine restFile_close

end module restFileMod



