program cprnc

  ! $Id$

  use header
  use nldat
  implicit none

  include 'netcdf.inc'

  ! Local workspace

  character*128 cvsid
  character*180 arg
  character*180 file(2)

  integer, dimension(2) :: nlon  = -1
  integer, dimension(2) :: nlat  = -1 
  integer, dimension(2) :: nlev  = -1

  integer n, nn, nargs
  integer numcases
  integer ntimemax
  integer ret(2)

  logical twocases
  logical itexists
  logical lnorm             ! print level by level L2 and L-INF norm stats
  logical verbose

#ifdef UNICOSMP
  integer ilen, ierr
#endif

  namelist/nlcpr/nsteps, nstepe, iprs, ipre, kprs, kpre, latprs, &
       latpre, matchts

  ! Externals

#ifdef UNICOSMP
  integer lenchr,ipxfargc
  external lenchr,ipxfargc
#else
  integer lenchr,iargc
  external lenchr,iargc
#endif

  !----------------------------------------------------------------------------

  cvsid = '$Id$'
  write(6,*)'You are using cprlndnc with cvsid:',cvsid

  ! Default settings before parsing argument list

  file(:) = ' '
  nsteps = 0
  nstepe = 99999999
  lnorm = .false.
  matchts = .true.
  twocases = .false.
  verbose = .false.

#ifdef UNICOSMP
  nargs = ipxfargc()
#else
  nargs = iargc()
#endif
  n = 1
  do while (n .le. nargs)
     arg = ' '
#ifdef UNICOSMP
     call pxfgetarg(n, arg, ilen, ierr)
#else
     call getarg (n, arg)
#endif
     n = n + 1
     if (arg .eq. '-e') then
#ifdef UNICOSMP
        call pxfgetarg(n, arg, ilen, ierr)
#else
        call getarg (n, arg)
#endif
        n = n + 1
        if (n.gt.nargs .or. arg(1:1).lt.'0' .or. arg(1:1).gt.'9') then
           call usage_exit ('The -e option requires an input ending timestep')
        end if
        read (arg, '(f15.7)') nstepe
        write (6,*) 'Stopping analysis after NSTEP=', nstepe
     else if (arg .eq. '-m') then
        matchts = .true.
     else if (arg .eq. '-n') then
        lnorm = .true.
     else if (arg .eq. '-s') then
#ifdef UNICOSMP
        call pxfgetarg(n, arg, ilen, ierr)
#else
        call getarg (n, arg)
#endif
        n = n + 1
        if (n.gt.nargs .or. arg(1:1).lt.'0' .or. arg(1:1).gt.'9') then
           call usage_exit ('The -s option requires an input start timestep')
        end if
        read (arg, '(f15.7)') nsteps
        write (6,*) 'Starting analysis at NSTEP=', nsteps
     else if (arg .eq. '-v') then
        verbose = .true.
     else
        if (file(1).eq.' ') then
           file(1) = arg(1:lenchr(arg))//char(0)
           write (6,*) 'file 1=',file(1)(1:lenchr(file(1)))
        else if (file(2).eq.' ') then
           file(2) = arg(1:lenchr(arg))//char(0)
           write (6,*) 'file 2=',file(2)(1:lenchr(file(2)))
           twocases = .true.
        else
           call usage_exit (' ')
        end if
     end if
  end do

  if (file(1).eq.' ') then
     call usage_exit ('You must enter at least 1 input file')
  else if (n.lt.nargs) then
     call usage_exit ('Either -d or -g must be specified')
  end if
  !
  ! Default namelist settings
  !
  iprs = 0
  ipre = 0
  kprs = 0
  kpre = 0
  latprs = 0
  latpre = 0

  inquire (file='nl.cpr', exist=itexists)
  if (itexists) then
     open (unit=9, file='nl.cpr', status='OLD')
     read (9,nlcpr)
     close (9)
     if (iprs.gt.0) then

        if (ipre.le.0 .or. ipre.lt.iprs) then
           write(6,*)'cprtps: IPRS set but IPRE invalid or not set'
           stop 99
        end if

      if (latprs.le.0 .or. latpre.lt.latprs) then
        write(6,*)'cprlndnc: LATPRS and/or LATPRE invalid or not set'
        stop 99
      end if

      if (kprs.le.0 .or. kpre.lt.kprs) then
        write(6,*)'cprlndnc: KPRS and/or KPRE invalid or not set'
        stop 99
      end if
    end if
  end if

  numcases = 1
  if (twocases) numcases = 2
!
! Open files, then get dimension and variable id info
!
  do n=1,numcases
    inquire (file=file(n), exist=itexists)
    if (.not.itexists) then
      write(6,*)'cprnc: Unable to find input file ',file(n)
      stop 99
    end if

    call wrap_open (file(n)(1:lenchr(file(n))), NF_NOWRITE, ncid(n))

    call wrap_get_att_text (ncid(n), NF_GLOBAL, 'case_id', case(n))
    call wrap_get_att_text (ncid(n), NF_GLOBAL, 'case_title', title(n))

    call wrap_inq_dimid (ncid(n), 'lon'   , londimid(n))
    call wrap_inq_dimid (ncid(n), 'lat'   , latdimid(n))
    call wrap_inq_dimid (ncid(n), 'levsoi', levdimid(n))
    call wrap_inq_dimid (ncid(n), 'time'  , unlimdimid(n))

    call wrap_inq_dimlen (ncid(n), londimid(n)  , nlon(n))
    call wrap_inq_dimlen (ncid(n), latdimid(n)  , nlat(n))
    call wrap_inq_dimlen (ncid(n), levdimid(n)  , nlev(n))
    call wrap_inq_dimlen (ncid(n), unlimdimid(n), ntime(n))

    call wrap_inq_varid (ncid(n), 'lon'   , lonid(n))
    call wrap_inq_varid (ncid(n), 'lat'   , latid(n))
    call wrap_inq_varid (ncid(n), 'levsoi', levid(n))
    call wrap_inq_varid (ncid(n), 'mcdate', ncdateid(n))
    call wrap_inq_varid (ncid(n), 'mcsec' , ncsecid(n))
    call wrap_inq_varid (ncid(n), 'area'  , areaid(n)) 
!
! Variables that should not result in a fatal error if not found
!
    ret(1) = nf_inq_varid (ncid(n), 'nstep', nstephid(n))
    if (ret(1).ne.NF_NOERR) nstephid(n) = -1

    ret(1) = nf_inq_varid (ncid(n), 'date_written', date_writtenid(n))
    if (ret(1).ne.NF_NOERR) date_writtenid(n) = -1

    ret(2) = nf_inq_varid (ncid(n), 'time_written', time_writtenid(n))
    if (ret(2).ne.NF_NOERR) time_writtenid(n) = -1
  end do
!
! Allocate space for header variables
!
  ntimemax = ntime(1)
  if (twocases .and. ntime(2).gt.ntime(1)) then
    ntimemax = ntime(2)
  end if

  allocate(date_written(ntimemax,numcases))
  allocate(time_written(ntimemax,numcases))
  allocate(nsteph(ntimemax,numcases))
  allocate(ncdate(ntimemax,numcases))
  allocate(ncsec(ntimemax,numcases))
  allocate(time(ntimemax))

  allocate(lon(nlon(1),numcases))
  allocate(lat(nlat(1),numcases))
  allocate(lev(nlev(1),numcases))
  allocate(area(nlon(1),nlat(1),numcases))

  if (twocases) allocate (time2(ntime(2)))
!
! Read in header variables
!
  do n=1,numcases
    ret(1) = nf_get_var_text (ncid(n), date_writtenid(n), date_written(1,n))
    if (ret(1).ne.NF_NOERR) then
      do nn=1,ntimemax
        date_written(nn,n) = 'UNKNOWN'
      end do
    end if

    ret(2) = nf_get_var_text (ncid(n), time_writtenid(n), time_written(1,n))
    if (ret(2).ne.NF_NOERR) then
      do nn=1,ntimemax
        time_written(nn,n) = 'UNKNOWN'
      end do
    end if

    ret(1) = nf_get_var_int (ncid(n), nstephid(n), nsteph(1,n))
    if (ret(1).ne.NF_NOERR) then
      do nn=1,ntimemax
        nsteph(nn,n) = -9999
      end do
    end if

    call wrap_get_var_int (ncid(n), ncdateid(n), ncdate(1,n))
    call wrap_get_var_int (ncid(n), ncsecid(n),  ncsec(1,n))

    call wrap_get_var_realx (ncid(n), lonid(n), lon(1,n))
    call wrap_get_var_realx (ncid(n), latid(n), lat(1,n))
    call wrap_get_var_realx (ncid(n), levid(n), lev(1,n))
    call wrap_get_var_realx (ncid(n), areaid(n), area(1,1,n))

  end do

! Currently no time variable so comment out for now
! call wrap_get_var_realx (ncid(1), timeid(1), time)
! if (twocases) call wrap_get_var_realx (ncid(2), timeid(2), time2)
  
  if (iprs.gt.0) then
    if (ipre.gt.nlon(1)) then
      write(6,*)'Bad ipre'
      write(6,*)'ipre,nlon=  ', ipre, nlon(1)
      stop 99
    end if

    if (iprs.lt.1) then
      write(6,*)'Bad iprs'
      write(6,*)'iprs=  ', iprs
      stop 99
    end if
  end if

! Print header diffs

  if (twocases) then
    call prhddiff (nlev(1), nlon(1), nlat(1))
  end if

! Print 1-time header

  write(6,'(/,a)')' SUMMARY OF FIELD DIFFERENCES:'
  write(6,802)'FIELD','# DIFFS','# POSS','MAX','MIN','DIFFMAX','VALUES','RDIFMAX','VALUES'

! Call the routine which loops through the time variables

  call cpr (nlon(1), nlat(1), nlev(1), numcases, lnorm, verbose)

802 format(1x ,a5,2x,a7,1x,a6,3x,a3,19x,a3,20x,a7,1x,a7,16x,a7,1x,a6)
end program cprnc

!===========================================================================

subroutine usage_exit (arg)
  implicit none
  character*(*) arg
  
  if (arg.ne.' ') write (6,*) arg
  write (6,*) 'Usage: cprlndnc [-c] [-s nsteps] [-e nstepe] [-m] ', &
              'file1 file2'
  stop 999
end subroutine usage_exit
