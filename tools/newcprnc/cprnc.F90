program cprnc
!-------------------------------------------------------------------------------------------
!
! Purpose: Driver program for executable named "cprnc". Analyze variable info on 1 or 2
!          CAM or CLM history files.  When 2 files are being analyzed, analyze difference
!          info as well.
! 
! Method:  Parse argument list. Check for rquired variables and attributes. Then call routine
!          which loops through time and gathers statistics.
!
! Primary call tree:
!
!   cprnc
!     parsepr <-------------------- Parse cmd-line arguments about printing
!     chkdims <-------------------- Check validity of dimensions when analyzing 2 files
!     timeloop <------------------- Driver routine for looping through time and variables
!       prhddiff <----------------- Print header differences when analyzing 2 files
!       prheader <----------------- Print header information
!       gather_stats <------------- Gather single-variable stats
!       print_stats <-------------- Print single-variable stats
!       gather_comparison_stats <-- Gather per-field per-timestep difference info
!       print_comparison_stats <--- Print per-field per-timestep difference info
!
! $Id$
!-------------------------------------------------------------------------------------------
   use utils,         only: endrun, lenchr
   use specialvalues, only: uninit_int
   use netcdfids,     only: unlimdimid, dateid, datesecid, nstephid, date_writtenid, time_writtenid, timeid
   use cntlvars,      only: iprs, ipre, jprs, jpre, kprs, kpre, verbose, lnorm, constps, matchts
   use nfwrappers,    only: wrap_open, wrap_inq_nvars, wrap_get_att_text, wrap_inq_dimid, &
                            wrap_inq_dimlen, wrap_inq_varid, wrap_inq_unlimdim, wrap_inq_var
   use driver,        only: timeloop

   implicit none
      
   include 'netcdf.inc'
!
! Local workspace
!
   character(len=128) :: cvsid                   ! for printing CVS version info
   character(len=180) :: arg                     ! cmd-line argument
   character(len=180) :: fname(2) = ' '          ! input filenames
   character(len=NF_MAX_NAME) :: caseid(2) = ' ' ! case ids (from namelist)
   character(len=NF_MAX_NAME) :: title(2) = ' '  ! case title (from namelist)
   character(len=NF_MAX_NAME) :: name            ! variable name

   integer :: nfid(2)  = -1                      ! file ids
   integer :: nvars(2) = uninit_int              ! number of variables on each of possibly 2 files
   integer :: ntime(2) = uninit_int              ! length of unlimited dimension
   integer :: n                                  ! index
   integer :: nargs                              ! number of cmd line args
   integer :: numcases = 1                       ! number of files (convenient to initialize to 1)
   integer :: nchars                             ! used in stripping nulls
   integer :: xtype                              ! variable type from netcdf
   integer :: ndims(2)                           ! number of dimensions in variable
   integer :: dimids(NF_MAX_DIMS,2)              ! dimension ids
   integer :: natts                              ! number of attributes
   
   logical :: itexists                           ! flag indicates 2 files are present
#ifdef UNICOSMP
   integer :: ilen, ierr                         ! length and error flag for pxfgetarg()
#endif
!
! Externals
!
#ifdef UNICOSMP
   integer, external :: ipxfargc                 ! number of cmd-line args
#else
   integer, external :: iargc                    ! number of cmd-line args
#endif

   cvsid = '$Id$'
   write(6,*)'You are using cprnc with cvsid:', cvsid
!
! Parse arg list
!
#ifdef UNICOSMP
   nargs = ipxfargc ()
#else
   nargs = iargc ()
#endif
   n = 1
   do while (n <= nargs)
      arg = ' '
#ifdef UNICOSMP
      call pxfgetarg (n, arg, ilen, ierr)
#else
      call getarg (n, arg)
#endif
      n = n + 1
      select case (arg)
      case ('-c')
         constps = .true.
         write (6,*)'Assuming constant PS everywhere (currently this changes nothing)'
      case ('-m')
         matchts = .false.
      case ('-n')
         lnorm = .true.
      case ('-v')
         verbose = .true.
      case ('-ipr')
#ifdef UNICOSMP
         call pxfgetarg (n, arg, ilen, ierr)
#else
         call getarg (n, arg)
#endif
         n = n + 1
         call parsepr (arg, iprs, ipre)
      case ('-jpr')
#ifdef UNICOSMP
         call pxfgetarg (n, arg, ilen, ierr)
#else
         call getarg (n, arg)
#endif
         n = n + 1
         call parsepr (arg, jprs, jpre)
      case ('-kpr')
#ifdef UNICOSMP
         call pxfgetarg (n, arg, ilen, ierr)
#else
         call getarg (n, arg)
#endif
         n = n + 1
         call parsepr (arg, kprs, kpre)
      case default
         if (fname(1) == ' ') then
            fname(1) = arg(1:len_trim(arg))//char(0)
            nchars = lenchr (fname(1))
            write (6,*) 'file 1=',fname(1)(1:nchars)
         else if (fname(2)==' ') then
            fname(2) = arg(1:len_trim(arg))//char(0)
            nchars = lenchr (fname(2))
            write (6,*) 'file 2=',fname(2)(1:nchars)
            numcases = 2
         else
            call usage_exit (' ')
         end if
      end select
   end do
!
! Must have at least 1 file input
!      
   if (fname(1) == ' ') then
      call usage_exit ('You must enter at least 1 input file')
   end if
!
! Open files, then get dimension and variable id info
!
   do n=1,numcases
      inquire (file=fname(n), exist=itexists)
      if (.not.itexists) then
         call endrun ('cprnc: Unable to find input file '//fname(n))
      end if

      call wrap_open (fname(n), NF_NOWRITE, nfid(n))
      call wrap_inq_nvars (nfid(n), nvars(n))
!
! A 5th argument to get_att_text, or a 4th argument to inq_varid results in
! the variable receiving that value if the item is not available on the
! netcdf file.  If the extra argument is not present and the item is not
! available, the code will abort.
!
      call wrap_get_att_text (nfid(n), NF_GLOBAL, 'case', caseid(n), 'Unspecified')
      call wrap_get_att_text (nfid(n), NF_GLOBAL, 'title', title(n), 'Unspecified')
      call wrap_inq_unlimdim (nfid(n), unlimdimid(n))
      call wrap_inq_dimlen (nfid(n), unlimdimid(n), ntime(n))

      call wrap_inq_varid (nfid(n), 'date', dateid(n), -1)
      call wrap_inq_varid (nfid(n), 'datesec', datesecid(n), -1)
      call wrap_inq_varid (nfid(n), 'nsteph', nstephid(n), -1)
      call wrap_inq_varid (nfid(n), 'date_written', date_writtenid(n), -1)
      call wrap_inq_varid (nfid(n), 'time_written', time_writtenid(n), -1)
      call wrap_inq_varid (nfid(n), 'time', timeid(n))
!
! Ensure variable named "time" is what we think it is
!
      call wrap_inq_var (nfid(n), timeid(n), name, xtype, ndims(n), dimids(:,n), natts)
      if (ndims(n) /= 1 .or. dimids(1,n) /= unlimdimid(n)) then
         call endrun ('cprnc: '//name//' is not a 1-d field of dimension NF_UNLIMITED')
      end if
   end do
!
! Check validity of print args.
!  
   if (iprs > 0 .or. jprs > 0 .or. kprs > 0) then
      if (ipre <= 0 .or. ipre < iprs) then
         write(6,*)'iprs,ipre=', iprs, ipre
         call endrun ('cprnc: unspecified or invalid "-ipr" input. syntax is -ipr start[:end]')
      end if
      
      if (jprs <= 0 .or. jpre < jprs) then
         write(6,*)'jprs,jpre=', jprs, jpre
         call endrun ('cprnc: unspecified or invalid "-jpr" input. syntax is -jpr start[:end]')
      end if
      
      if (kprs <= 0 .or. kpre < kprs) then
         write(6,*)'kprs,kpre=', kprs, kpre
         call endrun ('cprnc: unspecified or invalid "-kpr" input. syntax is -kpr start[:end]')
      end if
   end if
!
! Ensure that dimensionality matches between cases
!
   if (numcases > 1) then
      call chkdims (nfid)
   end if
!
! Print 1-time header info
!
   do n=1,numcases
      nchars = lenchr (caseid(n))
      write(6,*) 'CASE  ', n, ':', caseid(n)(1:nchars)
      nchars = lenchr (title(n))
      write(6,*) 'TITLE ', n, ':', title(n)(1:nchars)
   end do
   write(6,*)
!
! Call the routine which loops through variables in time and prints stats
!
   call timeloop (nfid, nvars, numcases, ntime)
   stop 0
end program cprnc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine usage_exit (arg)
   implicit none

   character(len=*), intent(in) :: arg
  
   if (arg /= ' ') write (6,*) arg
   write(6,*)'Usage: cprnc [-c ] [-m] [-n] [-v] [-ipr iprs[:ipre]] [-jpr jprs[:jpre]] [-kpr kprs[:kpre]] file1 [file2]'
   write(6,*)'-c: Assume constant surface pressure for mass-weighting (mass-weighting not yet implemented)'
   write(6,*)'-m: Compare each time sample. Default is false, i.e. match "time" coordinate values before comparing'
   write(6,*)'-n: Compute various norms, e.g. L-2, L-INF (not yet implemented)'
   write(6,*)'-v: Verbose output'
   write(6,*)'-ipr iprs[:ipre]: Print variable values for the specified 1st index subrange. If present,'
   write(6,*)'                  -jpr jprs[:jpre] and -kpr kprs[:kpre] must also be specified for 2nd and'
   write(6,*)'                  3rd index subranges, respectively. For variables with only 1 or 2 non-time'
   write(6,*)'                  dimensions, -jpr and -kpr settings which do not apply will be ignored.'

   stop 999
end subroutine usage_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine parsepr (arg, v1, v2)
!-------------------------------------------------------------------------------------------
! Purpose: Parse cmd line args about printing.
!
! Method:  Input is expected in the form: number1[:number2].  number1 is converted to v1.
!          If present, number2 is converted to v2.  Otherwise set v2 = v1
!-------------------------------------------------------------------------------------------
   implicit none

   character(len=*), intent(in) :: arg    ! cmd line arg expected of the form 'num1:num2' or 'num1'

   integer, intent(out) :: v1             ! e.g. num1 from above example
   integer, intent(out) :: v2             ! e.g. num2 or num1 (default) from above example
   
   integer :: i, j                        ! indices through arg
   integer :: ierr                        ! io error status
!
! First, try to get an integer number for everything up to ":"
!
   i = 1
   do while (i < len(arg) .and. arg(i:i) >= '0' .and. arg(i:i) <= '9' .and. arg(i:i) /= ':')
      i = i + 1
   end do
   read (arg(1:i-1), '(i5)') v1
!
! Next, if ":" comes after the number, look for the next number
!
   if (arg(i:i) == ':') then
      j = i + 1
      do while (j < len(arg) .and. arg(j:j) >= '0' .and. arg(j:j) <= '9')
         j = j + 1
      end do
      read (arg(i+1:j-1), '(i5)', iostat=ierr) v2
!
! On unexpected input set v2 = v1, e.g. "-ipr 2:blah" will mean the same
! as "-ipr 2:2"
!
      if (ierr /= 0) then
         v2 = v1
      end if
   else
!
! ":" not present. Interpret for example '-ipr 2' to mean '-ipr 2:2'
!
      v2 = v1
   end if

   return
end subroutine parsepr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chkdims (nfid)
!-------------------------------------------------------------------------------------------
! Purpose: Ensure conformability of dimensions when 2 files are being analyzed.
!
! Method:  Loop through dimensions of 1st file, ensure that the same dimension name
!          with the same size exists on the 2nd file.  The unlimited dimension
!          can have a different size.
!-------------------------------------------------------------------------------------------
   use utils,      only: endrun, lenchr
   use nfwrappers, only: wrap_inq_ndims, wrap_inq_dimname, wrap_inq_dimlen, wrap_inq_dimid
   use netcdfids,  only: unlimdimid
      
   implicit none
!
! Arguments
!
   integer, intent(in) :: nfid(2)          ! file ids
!
! Local workspace
!
   include 'netcdf.inc'

   character(len=NF_MAX_NAME) :: dimname   ! dimension name

   integer :: dimlen1, dimlen2             ! dimension lengths
   integer :: dimid                        ! dimension id
   integer :: ndims                        ! number of dimensions
   integer :: n                            ! index
!
! Sanity check
!
   if (any (nfid(:) < 0)) then
      call endrun ('chkdims: must have 2 active file ids')
   end if
!
! Loop over dimensions in file 1, ensure that dimensions of the same name and length
! exist on file 2.  Allow different lengths (but not names) for the unlimited dimension.
! Note that if file2 contains only a subset of the dimensions of file1 the code will
! fail, but if file1 contains only a subset of the dimensions of file2 it will work.
!
   call wrap_inq_ndims (nfid(1), ndims)
   do n=1,ndims
      dimlen1 = -1
      dimlen2 = -1
      call wrap_inq_dimname (nfid(1), n, dimname)
      call wrap_inq_dimlen (nfid(1), n, dimlen1)
      call wrap_inq_dimid (nfid(2), dimname, dimid)
      call wrap_inq_dimlen (nfid(2), dimid, dimlen2)
      if (dimlen1 == dimlen2) then
         write(6,*)'Dimension ',dimname(1:lenchr(dimname)),' has length ', dimlen1, ' on both files'
      else
         if (dimname == 'time') then
            write(6,*)'chkdims note: time dimension lengths differ:', dimlen1, dimlen2
         elseif (dimlen1 == -1 .or. dimlen2 == -1) then
            write(6,*)'Dimension ',dimname(1:lenchr(dimname)),' only exists on one file ', dimlen1,dimlen2
         else
            write(6,*)'Dimension ',dimname(1:lenchr(dimname)), ' has lengths ', dimlen1, dimlen2
!            call endrun ('chkdims')
         end if
      end if
   end do

   return
end subroutine chkdims
