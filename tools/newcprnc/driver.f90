module driver

   use prec,          only: r4, r8
   use types,         only: varinfo, diffinfo
   use cntlvars,      only: matchts, lnorm, constps
   use netcdfids,     only: unlimdimid, timeid
   use specialvalues, only: uninit_r8
   use gatherstats,   only: gather_stats, gather_comparison_stats
   use nfwrappers,    only: wrap_inq_var, wrap_get_vara_double, wrap_inq_unlimdim, wrap_inq_dimlen, &
                            wrap_inq_dimname, wrap_get_var_double
   use printers,      only: prhddiff, prheader, print_stats, print_comparison_stats
   use utils,         only: endrun

   implicit none

PRIVATE

   include 'netcdf.inc'

   public :: timeloop       ! driver routine which loops through time and calls stats-gathering routines
   private :: kosher        ! logical function determines whether field should be analyzed
   private :: get_fillvalue ! retrieve _FillValue attribute if present
   private :: setsizes      ! sets start, kount arrays for netcdf calls, and total array size

CONTAINS

   subroutine timeloop (nfid, nvars, numcases, ntime)
!-------------------------------------------------------------------------------------------
! Purpose: Loop through time samples, computing and printing stats for each variable
!
! Method:  For each time index, loop through the variables and print stats (max, min, etc.).
!          Print comparison stats for fields which are on both files
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid(2)             ! file id(s)
      integer, intent(in) :: nvars(2)            ! number of variables on each file
      integer, intent(in) :: numcases            ! number of files (1 or 2)
      integer, intent(in) :: ntime(2)            ! length of time dimension
!
! Local workspace
!
      real(r8), parameter :: timeepsilon = 1.e-9 ! time diff less than this considered insignificant
 
      character(len=NF_MAX_NAME) :: name         ! variable name
      character(len=NF_MAX_NAME) :: dimnames(3)  ! dimension names

      integer :: t(2)                   ! time indices
      integer :: t1, t2
      equivalence (t(1),t1)             ! equivalence for convenience
      equivalence (t(2),t2)

      integer :: n,nn                   ! file indices
      integer :: v                      ! variable index
      integer :: xtype                  ! variable type from netcdf
      integer :: numfcpr                ! number of fields compared (total)
      integer :: numfdif                ! number of fields which differed (total)
      integer :: ndims(2)               ! number of dimensions in variable
      integer :: dimsizes(3)            ! sizes of non-time (probably spatial) dimensions
      integer :: length(2)              ! size needed for buf allocation

      integer :: varid(2)               ! variable ids
      integer :: idum                   ! required arg to subroutine that is unused
      integer :: start(NF_MAX_DIMS,2)   ! shape input to get_vara routines
      integer :: kount(NF_MAX_DIMS,2)   ! shape input to get_vara routines
      integer :: dimids(NF_MAX_DIMS,2)  ! dimension ids
      integer :: natts                  ! number of attributes (required for netcdf calls)
      
      logical :: twocases               ! comparison or single file analysis?
      logical :: dodiffer               ! flag indicates differences were found

      real(r8) :: fillvalue(2)          ! _FillValue (if applicable).  Points are ignored.
      real(r8) :: time1(ntime(1))       ! time variable file 1
      real(r8) :: time2(ntime(2))       ! time variable file 2
      real(r8) :: timediff              ! time difference between samples when there are 2 files
!
! Declare data array as double precision instead of r8 for netcdf consistency
!
      double precision, allocatable :: buf(:,:)    ! field values

      type (varinfo) :: varstats(numcases)         ! variable statistics
      type (diffinfo) :: diffstats                 ! difference statistics
!
! Retrieve time variable(s) and print header diffs if 2 files are being analyzed
!
      call wrap_get_var_double (nfid(1), timeid(1), time1)
      twocases = numcases > 1
      if (twocases) then
         call wrap_get_var_double (nfid(2), timeid(2), time2)
         call prhddiff (nfid, numcases)
      end if
!
! Set initial values before looping
!
      t(:)    = 1
      numfcpr = 0
      numfdif = 0
!
! Top of loop over time samples: first match up times
!
1     continue
      if (twocases) then
         if (matchts) then
            do
               timediff = abs (time1(t1) - time2(t2))
!
! If times match then exit the loop
!
               if (timediff < timeepsilon) exit
               write(6,'(1x,a,f12.4,/1x,a,f12.4)') 'Current times differ:time1=', time1(t1), &
                                                   '                     time2=', time2(t2)
               if (time1(t1) < time2(t2)) then
                  write(6,*) 'Skipping a time sample on file 1'
                  t1 = t1 + 1
                  if (t1 > ntime(1)) then
                     write(6,*)'End of file 1 reached: stopping'
                     stop 0
                  end if
               else
                  write(6,*) 'Skipping a time sample on file 2'
                  t2 = t2 + 1
                  if (t2 > ntime(2)) then
                     write(6,*)'End of file 2 reached: stopping'
                     stop 0
                  end if
               end if
            end do
         end if

         call prheader (nfid, numcases, t, ntime)
         write(6,800) t1, t2, time1(t1), time2(t2)
         write(6,801)
         write(6,802)
         write(6,803)
         write(6,804)
         write(6,805)
         write(6,*)

800      format('SUMMARY of matching fields time samples ', 2i6, ' times ', 1p, 2e15.8, ':')
801      format(/'                        ( dim1, dim2, dim3)')
802      format( '                        (indx1,indx2,indx3) file 1')
803      format( ' FIELD    DIFFS NFILL1   MAX1          ',8x,'MIN1',19x,'DIFFMAX  VALUES',16x,'RDIFMAX VALUES')
804      format( '        ARRSIZE NFILL2   MAX2          ',8x,'MIN2')
805      format( '                        (indx1,indx2,indx3) file 2')
!
! Loop over file 1 variables, print comparison stats for those fields which are also on file 2
!
         do v=1,nvars(1)
            varid(1) = v
            dimids(:,1) = -1
            call wrap_inq_var (nfid(1), varid(1), name, xtype, ndims(1), dimids(:,1), natts)
!
! Skip to the next variable if type and shape criteria are not met (floating point, something by time)
!
            if (.not. kosher (nfid(1), name, ndims(1), dimids(:,1), unlimdimid(1), xtype)) then
               cycle
            end if
!
! Compare only if field also exists on file 2, and has the same characteristics as the file 1 field
!
            if (nf_inq_varid (nfid(2), name, varid(2)) == NF_NOERR) then
               dimids(:,2) = -1
               call wrap_inq_var (nfid(2), varid(2), name, xtype, ndims(2), dimids(:,2), natts)
               if (kosher (nfid(2), name, ndims(2), dimids(:,2), unlimdimid(2), xtype) .and. &
                   shapesmatch (nfid, ndims, dimids)) then
!
! Determine sizes, perform sanity check, then allocate buffer
!
                  do n=1,numcases
                     call setsizes (nfid(n), t(n), ndims(n), dimids(:,n), start(:,n), &
                                    kount(:,n), length(n), dimsizes, dimnames)
                  end do

                  if (length(1) /= length(2)) then
                     write(6,*)'lengths for buffer allocation do not match:', length(:)
                     call endrun ('driver')
                  end if

                  allocate (buf(length(1),numcases))
!
! Compute non-comparison stats for this variable on both files
!
                  do n=1,numcases
                     call wrap_get_vara_double (nfid(n), varid(n), start(:,n), kount(:,n), buf(:,n))
                     fillvalue(n) = get_fillvalue (nfid(n), varid(n))

                     call gather_stats (name, ndims(n), dimsizes(1), dimsizes(2), dimsizes(3), &
                                        buf(:,n), fillvalue(n), .false., dimnames, varstats(n))
                  end do
!
! Gather comparison stats, then release buffer space
!
                  call gather_comparison_stats (name, ndims(1), dimsizes(1), dimsizes(2), dimsizes(3), &
                                                buf, fillvalue, dimnames, diffstats, dodiffer)
                  deallocate (buf)
!
! Accumulate stats on number of fields compared and number that differ
! Then print the comparison stats
!
                  numfcpr = numfcpr + 1
                  if (dodiffer) then
                     numfdif = numfdif + 1
                  end if
                  call print_comparison_stats (varstats, diffstats)
               else
                  write(6,*)'driver:', name, ' is on both files but dimensions do not match'
               end if    ! variable shape and characteristics make it analyzable
            end if       ! variable exists on both files
         end do          ! v=1,nvars(1)
      end if             ! twocases
!
! Gather and print non-comparison stats (e.g. when analyzing only 1 file, or variable does
! not exist on both files)
!
      do n=1,numcases
         nn = mod(n,2) + 1
         if (n == 1) then
            write(6,900) n, t(n), time1(t(n))
         else if (n == 2) then
            write(6,900) n, t(n), time2(t(n))
         end if
         write(6,801)
         write(6,802)
         write(6,904)
         write(6,*)

900      format('SUMMARY of fields only on file ', i1, ' time sample', i6, ' time ', 1p, e15.8, ':')
904      format(' FIELD  ARRSIZE NFILL    MAX                   MIN')

         do v=1,nvars(n)
            varid(n) = v
            dimids(:,n) = -1
            call wrap_inq_var (nfid(n), varid(n), name, xtype, ndims(n), dimids(:,n), natts)
!
! If data are floating point and shape is something by time, then analyze if field is not 
! also on other the file. If it is on the other file, stats have already been printed in 
! the "do v=1,nvars(1)" loop above
!
            if (kosher (nfid(n), name, ndims(n), dimids(:,n), unlimdimid(n), xtype)) then
               if (.not. twocases .or. nf_inq_varid (nfid(nn), name, idum) /= NF_NOERR) then
!
! Determine sizes, allocate buffer, read in the variable, then gather and print stats
!                               
                  call setsizes (nfid(n), t(n), ndims(n), dimids(:,n), start(:,n), &
                                 kount(:,n), length(n), dimsizes, dimnames)

                  allocate (buf(length(n),1))
                  call wrap_get_vara_double (nfid(n), varid(n), start(:,n), kount(:,n), buf)
                  fillvalue(n) = get_fillvalue (nfid(n), varid(n))

                  call gather_stats (name, ndims(n), dimsizes(1), dimsizes(2), dimsizes(3), &
                                     buf, fillvalue(n), .true., dimnames, varstats(n))
                  deallocate (buf)

                  call print_stats (varstats(n))
               end if
            end if
         end do     ! v=1,nvars(n)
      end do        ! n=1,numcases

      write(6,806)
806   format(132('*'))
!
! Clean up if at end of file(s), or increment time counter and continue
!
      if (twocases .and. (t1 >= ntime(1) .or. t2 >= ntime(2))) then
         do n=1,numcases
            if (t(n) < ntime(n)) then
               write(6,*) 'NOTE: ', ntime(n)-t(n), ' time samples at the end of file ', n, ' were skipped'
            end if
         end do

         write(6,'(a,i8,a,i8,a)') 'A total of ', numfcpr, ' fields were compared, of which', &
                                  numfdif,' had non-zero diffs'
         return
      end if

      if (.not. twocases .and. t1 >= ntime(1)) then
         return
      end if
!
! Move to next time index in file(s)
!
      t(:) = t(:) + 1
      goto 1

   end subroutine timeloop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine setsizes (nfid, t, ndims, dimids, start, &
                        kount, length, dimsizes, dimnames)
!-------------------------------------------------------------------------------------------
! Purpose: Determine dimension sizes and total length, and set start and kount arrays
!          for input to nf_get_vara_double
!
! Method:  Call appropriate netcdf inquiry functions
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid                  ! file id
      integer, intent(in) :: t                     ! time index
      integer, intent(in) :: ndims                 ! number of dimensions in variable
      integer, intent(in) :: dimids(NF_MAX_DIMS)   ! dimension ids of variable

      integer, intent(out) :: start(NF_MAX_DIMS)   ! for input to nf_get_vara routines
      integer, intent(out) :: kount(NF_MAX_DIMS)   ! for input to nf_get_vara routines
      integer, intent(out) :: length               ! total array size (excluding time index)
      integer, intent(out) :: dimsizes(3)          ! individual non-time dimension sizes

      character(len=*), intent(out) :: dimnames(3) ! dimension names
!
! Local workspace
!
      integer :: n                                 ! dimension index
      integer :: unlimdimid_loc                    ! to distinguish from global variable unlimdimid
!
! First ensure that last dimension id on variable is unlimited (sanity check)
!
      call wrap_inq_unlimdim (nfid, unlimdimid_loc)
      if (dimids(ndims) /= unlimdimid_loc) then
         write(6,*) 'setsizes: dimids(', ndims, ') does not equal unlimited dimension id=', unlimdimid_loc
         call endrun ('setsizes')
      end if
!
! ndims > 4 should have been skipped and a msg printed by function kosher
!
      if (ndims > 4) then
         call endrun ('setsizes: ndims too big')
      end if

      start(:) = -1
      kount(:) = -1
      length = 1
!
! Loop through dimensions, setting start, kount, dimnames, dimsizes, and length for output
!
      do n=1,ndims-1
         call wrap_inq_dimlen (nfid, dimids(n), kount(n))
         call wrap_inq_dimname (nfid, dimids(n), dimnames(n))
         start(n)    = 1
         length      = length * kount(n)
         dimsizes(n) = kount(n)
      end do
!
! Unlimited dimension: want 1 slice at time index t
!
      start(ndims) = t
      kount(ndims) = 1
!
! Set inapplicable dimension sizes to 1 (convenience for gather_stats routines)
! Set inapplicable dimension names to dashes for later printing (print_stats routines)
!
      dimsizes(ndims:3) = 1
      dimnames(ndims:3) = '-----'

      return
   end subroutine setsizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function kosher (nfid, name, ndims, dimids, unlimdimid_loc, xtype)
!-------------------------------------------------------------------------------------------
! Purpose: Determine whether input field has appropriate characterisics for analysis
!
! Method:  Examine input vars and make appropriate netcdf calls to see if the field
!          is floating point, dimensioned something by time, where "something" may
!          be 1, 2, or 3 dimensions
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid                ! file id
      character(len=*), intent(in) :: name       ! variable name
      integer, intent(in) :: ndims               ! number of dimensions associated with variable
      integer, intent(in) :: dimids(NF_MAX_DIMS) ! dimension ids
      integer, intent(in) :: unlimdimid_loc      ! _loc distinguishes from global variable
      integer, intent(in) :: xtype               ! variable type (netcdf)

      integer :: n                               ! dimension index
      integer :: dimlen                          ! dimension length

      if (unlimdimid_loc < 1 .or. xtype < 0) then
         call endrun ('kosher: bad input')
      end if
!
! No more than 3 dimensions plus time allowed
!
      if (ndims > 4) then
         kosher = .false.
         write(6,*)'kosher: variable ', name, ' has too many dimensions'
         return
      end if
!
! Ensure that all dimensions have length > 0
!
      do n=1,ndims
         call wrap_inq_dimlen (nfid, dimids(n), dimlen)
         if (dimlen == 0) then
            kosher = .false.
            return
         end if
      end do
!
! Variable must be of type float, and dimensioned something by time
!
      if (ndims > 0 .and. (xtype == NF_FLOAT .or. xtype == NF_DOUBLE) .and. &
          dimids(1) /= unlimdimid_loc .and. dimids(ndims) == unlimdimid_loc) then
         kosher = .true.
      else
         kosher = .false.
      end if

      return
   end function kosher

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function shapesmatch (nfid, ndims, dimids)
!-------------------------------------------------------------------------------------------
! Purpose: Determine whether shapes of 2 netcdf variables match
!
! Method:  Make netcdf calls to get dimension names and lengths. Both names AND sizes
!          must all match in order for shapes to match
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid(2)               ! file ids
      integer, intent(in) :: ndims(2)              ! number of dimensions
      integer, intent(in) :: dimids(NF_MAX_DIMS,2) ! dimension ids
!
! Local workspace
!
      integer :: n                                 ! dimension index
      integer :: dimlen1, dimlen2                  ! dimension lengths
      
      character(len=NF_MAX_NAME) :: dimname1, dimname2 ! dimension names
!
! Check on sanity of input
!
      if (ndims(1) < 1 .or. ndims(2) < 1) then
         call endrun ('shapesmatch: bad ndims input')
      end if

      if (ndims(1) /= ndims(2)) then
         shapesmatch = .false.
         return
      end if

      do n=1,ndims(1)
         call wrap_inq_dimname (nfid(1), dimids(n,1), dimname1)
         call wrap_inq_dimname (nfid(2), dimids(n,2), dimname2)
!
! If dimension names do not match return false
! One might allow the unlimited dimension name to differ, but we do not
!
         if (dimname1 /= dimname2) then
            shapesmatch = .false.
            return
         end if
!
! If dimension lengths other than the unlimited dimension do not match return false
!
         call wrap_inq_dimlen (nfid(1), dimids(n,1), dimlen1)
         call wrap_inq_dimlen (nfid(2), dimids(n,2), dimlen2)
         if (dimlen1 /= dimlen2) then
            if (dimids(n,1) /= unlimdimid(1) .or. dimids(n,2) /= unlimdimid(2)) then
               shapesmatch = .false.
               return
            end if
         end if
      end do

      shapesmatch = .true.
      return
   end function shapesmatch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   real(r8) function get_fillvalue (nfid, varid)
!-------------------------------------------------------------------------------------------
! Purpose: Determine _FillValue attribute for a given variable (if present)
!
! Method:  Make appropriate netcdf calls.  If _FillValue is not present, return uninit_r8
!-------------------------------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nfid   ! file id
      integer, intent(in) :: varid  ! variable id
!
! Local workspace
!
      integer :: xtype               ! variable type (netcdf)
      integer :: numelem             ! number of elements in attribute
      integer :: ret                 ! return code

      double precision :: fillvalue  ! value of _FillValue attribute
!
! If fillvalue is there, determine its value and return it. Otherwise return uninit_r8
!
      get_fillvalue = uninit_r8
      if (nf_inq_att (nfid, varid, '_FillValue', xtype, numelem) == NF_NOERR) then
         if ((xtype == NF_REAL .or. xtype == NF_DOUBLE) .and. numelem == 1) then
            ret = nf_get_att_double (nfid, varid, '_FillValue', fillvalue)
            get_fillvalue = fillvalue
         end if
      end if

      return
   end function get_fillvalue
end module driver
